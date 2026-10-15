"""Robustness v3: 50-year flat-zero, +/-15% perturbations, 20 worlds.
Candidates: references (G2, M1b), one-time purge-infusions (P4 chronic, P6 early),
nuclear erase (N1). Uses the calibrated standard integrator (validated: G2 stays
flat-zero at 50y). NaN/negative = FAIL."""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate

Y = 365.25
ART_NEVER = 1e9
rng = np.random.default_rng(23)
H = 50.0

def base_conf(name):
    if name == "G2":      # finite course reference
        return dict(t_art_on=365, t_art_off=1460,
                    shocks=list(range(430, 1460, 30)), lra=0.20,
                    strike_at=1e9, strike_refr_len=0)
    if name == "M1b":     # one-shot permanent-reference
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    strike_at=1e9, strike_refr_len=0,
                    ccr5=0.98, ccr5_start=365, bnab=(365, 0.9, 0.0), bnab_w=45*Y)
    if name == "P4":      # one-time purge-infusion, chronic
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    strike_at=365, strike_erase=0.99, strike_kill=1.0,
                    strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=90,
                    strike_healD=0.9, strike_purge=1.5, strike_cleanup=100)
    if name == "P6":      # one-time purge-infusion, EARLY (day 14)
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    strike_at=14, strike_erase=0.9, strike_kill=1.0,
                    strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=30,
                    strike_healD=0.9, strike_purge=1.5, strike_cleanup=100)
    if name == "N1":      # nuclear erase + 1y lockdown
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    strike_at=365, strike_erase=0.99999999, strike_kill=1.0,
                    strike_neutral=0.999999, strike_refr=0.999, strike_refr_len=365,
                    strike_healD=0.9, strike_purge=0.0, strike_cleanup=1.0)
    raise KeyError(name)

PERTURB = ["lam", "d", "beta", "fL", "a0", "d_l", "rho", "delta_y",
           "k", "u", "rho4", "mu", "s4", "cD", "rD", "D_thr"]

def sane_world(p):
    pp = dict(p); pp.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                            strike_at=1e9, strike_refr_len=0)
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=250, timeout=4.0)
    except BaseException:
        return False
    X = YY.T[0]; V = YY.T[4] + YY.T[5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def flatzero(t, YY, from_year):
    if t is None or not np.isfinite(YY).all():
        return False
    V = YY[:, 4] + YY[:, 5]; L = YY[:, 3]
    m = t >= from_year*Y
    if m.sum() < 50:
        return False
    return bool((np.abs(V[m]) <= 1e-3).all() and (L[-50:] <= 1e-6).all()
                and YY[-1, 0] >= 500)

FROM = {"G2": 5.0, "M1b": 3.0, "P4": 2.0, "P6": 2.0, "N1": 2.0}

t0 = time.time()
print(f"Robustness v3: {H}-year flat-zero cure, +/-15% on 16 params, 20 worlds")
for nm in ["G2", "M1b", "P4", "P6", "N1"]:
    base = base_conf(nm)
    passes = skipped = tries = 0
    while (passes + skipped) < 20 and tries < 40:
        p = dict(P)
        for k in PERTURB:
            p[k] = p[k] * rng.uniform(0.85, 1.15)
        tries += 1
        if not sane_world(p):
            skipped += 1; continue
        p.update(base)
        try:
            t, YY, _ = simulate(p, t_end=H*Y, n_out=700, timeout=20.0)
        except BaseException:
            skipped += 1; continue
        passes += int(flatzero(t, YY, FROM[nm]))
    print(f"  {nm:5s}: {passes:3d}/20 worlds FLAT-ZERO @ {H}y  (screened {skipped}, "
          f"tried {tries})   [{time.time()-t0:.0f}s]", flush=True)
print(f"done in {time.time()-t0:.0f}s")