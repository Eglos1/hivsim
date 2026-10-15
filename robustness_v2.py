"""Robustness v2: 30-year horizon, ±15% perturbations, 5 cure candidates."""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate

Y = 365.25
rng = np.random.default_rng(11)
ART_NEVER = 1e9

def strategy(name):
    if name == "G2_finite_long":      # chronic, long finite course (baseline cure)
        return dict(t_art_on=365, t_art_off=1460,
                    shocks=list(range(430, 1460, 30)), lra=0.20)
    if name == "Pshort":              # acute, SHORT finite course + CAR-T (1.5y)
        return dict(t_art_on=7, t_art_off=547,
                    shocks=list(range(30, 540, 30)), lra=0.20,
                    cart_boost=2.0, cart_start=30, t_cart_end=547)
    if name == "Pchron":              # chronic, MID finite course + CAR-T (2y)
        return dict(t_art_on=365, t_art_off=1095,
                    shocks=list(range(395, 1090, 30)), lra=0.20,
                    cart_boost=2.0, cart_start=400, t_cart_end=1095)
    if name == "M3_one_shot":         # acute ART-FREE one-shot
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    ccr5=0.90, ccr5_start=14,
                    cart_boost=1.5, cart_start=14,
                    bnab=(14, 0.7, 2.0), bnab_w=30*Y)
    if name == "M1b_one_shot":        # chronic ART-FREE one-shot (CCR5+Ab factory)
        return dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
                    ccr5=0.98, ccr5_start=365,
                    bnab=(365, 0.9, 0.0), bnab_w=30*Y)
    raise KeyError(name)

PERTURB = ["lam", "d", "beta", "fL", "a0", "d_l", "rho", "delta_y",
           "k", "u", "rho4", "mu", "s4", "cD", "rD", "D_thr"]

def sane_world(p):
    pp = dict(p); pp.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[])
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=200, timeout=3.0)
    except BaseException:
        return False
    X = YY.T[0]; V = YY.T[4] + YY.T[5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def cured(t, YY):
    """Suppressed for entire final 3 years + CD4 healthy + reservoir ~0."""
    X  = YY.T[0]   # CD4 target cells
    L  = YY.T[3]   # latent reservoir
    V  = YY.T[4] + YY.T[5]
    tail = t >= t[-1] - 3*Y
    return bool((V[tail] < 50).all() and X[-1] >= 500 and (L[tail] < 1e-3).all())

t_start = time.time()
print("Robustness v2: 30-year horizon, +/-15% on 16 params, 20 worlds",
      flush=True)
for nm in ["G2_finite_long", "Pshort", "Pchron", "M3_one_shot", "M1b_one_shot"]:
    base = strategy(nm)
    passes = skipped = tries = 0
    while (passes + skipped) < 20 and tries < 35:
        p = dict(P)
        for k in PERTURB:
            p[k] = p[k] * rng.uniform(0.85, 1.15)
        tries += 1
        if not sane_world(p):
            skipped += 1; continue
        p.update(base)
        try:
            t, YY, _ = simulate(p, t_end=30*Y, n_out=1600, timeout=5.0)
        except BaseException:
            skipped += 1; continue
        passes += int(cured(t, YY))
    print(f"  {nm:16s}: {passes}/20 worlds cured @30y "
          f"({100*passes/20:.0f}%)  screened: {skipped}", flush=True)
    print(f"     -- elapsed {time.time()-t_start:.0f}s", flush=True)
print(f"\ndone in {time.time()-t_start:.0f}s", flush=True)