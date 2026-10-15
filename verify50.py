"""50-year, tightened-solver verification of every one-shot candidate.
Requires FLAT ZERO virus (not merely <50) from year 2 to year 50,
reservoir extinct, CD4 healthy. NaN/artifacts = FAIL."""
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from scipy.integrate import solve_ivp
from model import P, _events, make_deriv

Y = 365.25
ART_NEVER = 1e9

def sim_strict(par, t_end, n_out, rtol=1e-7, atol=1e-9, timeout=None):
    """Stricter piecewise integrator (tighter tolerances, smaller max step)."""
    import signal as _sig
    t0 = 0.0
    y0 = np.array((par["lam"]/par["d"], 0.,0.,0.,50.,0.,0.))
    teval = np.linspace(t0, t_end, n_out)
    if timeout:
        class Kit(BaseException): pass
        old = _sig.getsignal(_sig.SIGALRM)
        _sig.signal(_sig.SIGALRM, lambda s,f:(_ for _ in ()).throw(Kit()))
        _sig.setitimer(_sig.ITIMER_REAL, timeout)
    times, states = [t0], [y0.copy()]
    Yside = y0.copy()
    try:
        ev = _events(par, t0, t_end)
        for a, b in zip(ev[:-1], ev[1:]):
            if b <= a: continue
            sub = teval[(teval >= a) & (teval <= b)]
            if len(sub) < 2: sub = np.array([a, b])
            res = solve_ivp(make_deriv(par), (a, b), Yside, method="LSODA",
                            t_eval=sub, rtol=rtol, atol=atol,
                            min_step=1e-9, max_step=1.0)
            if not res.success:
                return None, None
            Yside = res.y[:, -1]
            np.maximum(Yside, 0.0, out=Yside)
            if par["strike_at"] < 1e9 and abs(b - par["strike_at"]) < 1e-6:
                s = par
                Yside[3] *= (1 - s["strike_erase"])
                Yside[1] *= (1 - s["strike_kill"])
                Yside[2] *= (1 - s["strike_kill"])
                Yside[4] *= (1 - s["strike_neutral"])
                Yside[5] *= (1 - s["strike_neutral"])
                Yside[6] *= (1 - s["strike_healD"])
                if s["strike_boostx"] != 1.0:
                    Yside[0] = min(Yside[0]*s["strike_boostx"], 1.5*par["lam"]/par["d"])
                np.maximum(Yside, 0.0, out=Yside)
            times.extend(res.t[1:])
            states.extend(res.y.T[1:])
    except BaseException:
        return None, None
    finally:
        if timeout:
            _sig.setitimer(_sig.ITIMER_REAL, 0)
            _sig.signal(_sig.SIGALRM, old)
    return np.array(times), np.array(states)

def flatzero(t, YY, from_year=2.0, to_year=50.0):
    """Virus is numerically ZERO (<=1e-3 copies) for the whole flat window,
       reservoir extinct, CD4 healthy at end. No NaN allowed."""
    if t is None or not np.isfinite(YY).all():
        return False
    V = YY[:,4] + YY[:,5]
    L = YY[:,3]
    m = (t >= from_year*Y) & (t <= to_year*Y)
    if m.sum() < 10:
        return False
    return bool((np.abs(V[m]) <= 1e-3).all()
                and (L[m[-50:]] <= 1e-6).all()
                and YY[m[-1], 0] >= 500)

def banner(name, ov, stage="chronic"):
    day = {"early": 14, "chronic": 365}[stage]
    p = dict(P); p.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                          strike_at=day, strike_kill=1.0, strike_neutral=0.999999,
                          strike_healD=0.9, **ov)
    t, YY = sim_strict(p, 50*Y, 1000, timeout=8.0)
    ok = flatzero(t, YY)
    tag = "FLAT-ZERO 50y CURE" if ok else ("numerics/NaN", "REBOUND or partial")[0 if t is None else 1]
    print(f"[{tag:22s}] {name}")
    return ok

print("======== 50-year flat-zero verification (tight solver) ========")
print("--- reference designs that ALREADY proved robust ---")
banner("REF G2: chronic ART 3y + exhaust purge (finite, no permanents)",
       dict(), stage="chronic")
# G2 uses ART+shocks, not strike: handle specially
def g2_ref():
    p = dict(P); p.update(t_art_on=365, t_art_off=1460,
                          shocks=list(range(430,1460,30)), lra=0.20)
    t, YY = sim_strict(p, 50*Y, 1000, timeout=12.0)
    ok = flatzero(t, YY, from_year=4.0)
    print(f"[{'FLAT-ZERO 50y CURE' if ok else 'FAIL':22s}] REF G2: chronic ART 3y + exhaust purge")
    return ok
g2_ref()
banner("REF M1b: one-shot CCR5+permanent Ab factory", dict(ccr5=0.98, ccr5_start=365,
       bnab=(365,0.9,0.0), bnab_w=45*Y, strike_at=1e9), stage="chronic")

print("\n--- PURGE-INFUSION one-shots (prior 30y 'cures') ---")
banner("P1 chronic: erase.999 win30 pur1.5 cl10",
       dict(strike_erase=0.999, strike_refr_len=30, strike_purge=1.5, strike_cleanup=10))
banner("P2 chronic: erase.999 win90 pur0.5 cl10",
       dict(strike_erase=0.999, strike_refr_len=90, strike_purge=0.5, strike_cleanup=10))
banner("P3 chronic: erase.999 win60 pur3 cl10",
       dict(strike_erase=0.999, strike_refr_len=60, strike_purge=3.0, strike_cleanup=10))
banner("P4 chronic: erase.99 win90 pur1.5 cl100",
       dict(strike_erase=0.99, strike_refr_len=90, strike_purge=1.5, strike_cleanup=100))
banner("P5 chronic: erase.999 win30 pur1.5 cl100",
       dict(strike_erase=0.999, strike_refr_len=30, strike_purge=1.5, strike_cleanup=100))
banner("P6 early: erase.9 win30 pur1.5 cl100", dict(strike_erase=0.9,
       strike_refr_len=30, strike_purge=1.5, strike_cleanup=100), stage="early")
banner("P7 early: erase.9 win15 pur1.5 cl100", dict(strike_erase=0.9,
       strike_refr_len=15, strike_purge=1.5, strike_cleanup=100), stage="early")

print("\n--- NUCLEAR ERASE one-shot ---")
banner("N1 chronic: erase0.99999999 + 365d lockdown",
       dict(strike_erase=0.99999999, strike_refr_len=365, strike_refr=0.999,
            strike_purge=0.0, strike_cleanup=1.0))