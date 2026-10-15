"""Debug: 50y behavior of reference designs (G2, M1b) and one-shot designs,
tight solver, trajectory diagnostics on decade marks."""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import solve_ivp
from model import P, _events, make_deriv

Y = 365.25
ART_NEVER = 1e9

def sim(par, t_end, n_out, rtol=1e-7, atol=1e-9):
    t0 = 0.0
    y0 = np.array((par["lam"]/par["d"],0.,0.,0.,50.,0.,0.))
    teval = np.linspace(t0, t_end, n_out)
    times, states = [t0], [y0.copy()]
    Yside = y0.copy()
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
            Yside[0] = min(Yside[0]*s["strike_boostx"], 1.5*par["lam"]/par["d"])
            np.maximum(Yside, 0.0, out=Yside)
        times.extend(res.t[1:])
        states.extend(res.y.T[1:])
    return np.array(times), np.array(states)

def report(name, p, from_yr=5.0):
    t, YY = sim(p, 50*Y, 1500)
    if t is None:
        print(f"{name}: integration FAILED"); return
    V = YY[:,4]+YY[:,5]; L = YY[:,3]; X = YY[:,0]; yr = t/Y
    print(f"{name}")
    for yy in [5,10,15,20,25,30,35,40,45,50]:
        i = np.argmin(np.abs(yr-yy))
        print(f"   yr {yy:2d}: X={X[i]:7.1f}  V={V[i]:11.3e}  L/ml={L[i]*1000:9.3e}")
    m = yr > from_yr
    over = np.where(V[m] > 1e-3)[0]
    if len(over):
        print(f"   -> first VL>1e-3 off-window at yr {yr[m][over[0]]:.1f} "
              f"(V={V[m][over[0]]:.2e}), and max={V[m].max():.2e}")
    else:
        print(f"   -> V flat-zero entire post window {from_yr}-50y")
    print()

ART = dict(P); ART.update(t_art_on=365, t_art_off=1460,
                          shocks=list(range(430,1460,30)), lra=0.20)
report("G2  chronic ART 3y + exhaust purge", ART)

M1b = dict(P); M1b.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                          ccr5=0.98, ccr5_start=365, bnab=(365,0.9,0.0), bnab_w=45*Y)
report("M1b one-shot CCR5 98% + Ab factory 45y", M1b)

N1 = dict(P); N1.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                        strike_at=365, strike_erase=0.99999999, strike_kill=1.0,
                        strike_neutral=0.999999, strike_refr=0.999, strike_refr_len=365,
                        strike_healD=0.9, strike_purge=0.0, strike_cleanup=1.0)
report("N1 chronic one-time erase 1e-8 residual + 365d lock", N1)

P1 = dict(P); P1.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                        strike_at=365, strike_erase=0.999, strike_kill=1.0,
                        strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=30,
                        strike_healD=0.9, strike_purge=1.5, strike_cleanup=10)
report("P1 chronic one-time purge-infusion erase.999 win30 pur1.5", P1, from_yr=2.0)

P6 = dict(P); P6.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                        strike_at=14, strike_erase=0.9, strike_kill=1.0,
                        strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=30,
                        strike_healD=0.9, strike_purge=1.5, strike_cleanup=100)
report("P6 early (d14) one-time purge-infusion", P6, from_yr=2.0)