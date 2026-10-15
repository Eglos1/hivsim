"""Miracle hunt: chronic-stage one-shot protocols that hit 20/20 flat-zero
at 50 years, with the FIXED interpolation. 20 perturbed worlds per design."""
import warnings, time
warnings.filterwarnings("ignore")
import numpy as np
from model import P, simulate

Y = 365.25
ART_NEVER = 1e9
rng = np.random.default_rng(41)
PERTURB = ["lam", "d", "beta", "fL", "a0", "d_l", "rho", "delta_y",
           "k", "u", "rho4", "mu", "s4", "cD", "rD", "D_thr"]

def base_LRA(p):  # shared strike scaffolding
    p.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
             strike_kill=1.0, strike_neutral=0.999999, strike_refr=0.9999,
             strike_healD=0.9)

def conf(name):
    p = dict(P)
    if name == "N1x":      # nuclear erase, more headroom
        base_LRA(p); p.update(strike_at=365, strike_erase=0.999999999,
                              strike_refr=0.999, strike_refr_len=365,
                              strike_purge=0.0, strike_cleanup=1.0)
    elif name == "N1xp":   # + purge/cleanup inside the 365d lock
        base_LRA(p); p.update(strike_at=365, strike_erase=0.99999999,
                              strike_refr=0.999, strike_refr_len=365,
                              strike_purge=1.5, strike_cleanup=100)
    elif name == "P4d":    # 99-day purge + ART folded in, cleanup 1000
        base_LRA(p); p.update(t_art_on=365, t_art_off=455,
                              strike_at=365, strike_erase=0.999,
                              strike_refr_len=90, strike_purge=2.0,
                              strike_cleanup=1000)
    elif name == "P4e":    # 120-day burn + ART 4 months
        base_LRA(p); p.update(t_art_on=365, t_art_off=485,
                              strike_at=365, strike_erase=0.9999,
                              strike_refr_len=120, strike_purge=3.0,
                              strike_cleanup=1000)
    elif name == "DT1":    # double-tap: two 45d burns, ART 4.5 months
        base_LRA(p); p.update(t_art_on=365, t_art_off=500,
                              strike_at=[365, 410],
                              strike_erase=[0.99, 0.99],
                              strike_refr_len=[45, 45],
                              strike_purge=[1.5, 1.5],
                              strike_cleanup=[1000, 1000])
    elif name == "DT2":    # double-tap, heavier purge
        base_LRA(p); p.update(t_art_on=365, t_art_off=500,
                              strike_at=[365, 410],
                              strike_erase=[0.99, 0.995],
                              strike_refr_len=[60, 60],
                              strike_purge=[3.0, 3.0],
                              strike_cleanup=[1000, 1000])
    else:
        raise KeyError(name)
    return p

def sane(p):
    pp = dict(p); pp.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER,
                            shocks=[], strike_at=1e9, strike_refr_len=0)
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=250, timeout=4)
    except BaseException:
        return False
    X = YY[:, 0]; V = YY[:, 4] + YY[:, 5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def flatzero(t, YY, fy=2.0):
    if t is None or not np.isfinite(YY).all():
        return False
    V = YY[:, 4] + YY[:, 5]; L = YY[:, 3]
    m = t >= fy*Y
    return bool((np.abs(V[m]) <= 1e-3).all() and (L[-50:] <= 1e-6).all()
                and YY[-1, 0] >= 500)

t0 = time.time()
print("MIRACLE HUNT: chronic one-shots, 50y flat-zero, 20 worlds each (fixed solver)")
for nm in ["N1x", "N1xp", "P4d", "P4e", "DT1", "DT2"]:
    passc = skip = tries = 0
    while passc + skip < 20 and tries < 45:
        p = dict(P)
        for k in PERTURB:
            p[k] *= rng.uniform(0.85, 1.15)
        tries += 1
        if not sane(p):
            skip += 1; continue
        p.update(conf(nm))
        try:
            t, YY, _ = simulate(p, t_end=50*Y, n_out=800, timeout=25)
        except BaseException:
            skip += 1; continue
        passc += int(flatzero(t, YY))
    print(f"  {nm:5s}: {passc:2d}/20  (screened {skip}, tried {tries}) "
          f"[{time.time()-t0:.0f}s]", flush=True)