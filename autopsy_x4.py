"""Deep autopsy: X4 ember life-cycle in boundary-gang relapses."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from model import P, simulate
Y = 365.25; ART_NEVER = 1e9
PERTURB = ["lam","d","beta","fL","a0","d_l","rho","delta_y","k","u",
           "rho4","mu","s4","cD","rD","D_thr"]

def sane(p):
    pp = dict(p); pp.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER,
                            shocks=[], strike_at=1e9, strike_refr_len=0)
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=250, timeout=4)
    except BaseException:
        return False
    X = YY[:,0]; V = YY[:,4] + YY[:,5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def PAN(pur=3.0):
    return dict(t_art_on=365, t_art_off=1460, shocks=[],
                strike_at=[365], strike_erase=[0.999999999],
                strike_kill=[1.0], strike_neutral=[0.999999],
                strike_refr=[0.9999], strike_refr_len=[1095],
                strike_healD=[0.9], strike_purge=[pur],
                strike_cleanup=[1000], ccr5=0.999, ccr5_start=365)

worlds = []
rng = np.random.default_rng(777)
while len(worlds) < 60:
    p = dict(P)
    for k in PERTURB:
        p[k] *= rng.uniform(0.85, 1.15)
    if sane(p):
        worlds.append(p)

for wi in (2, 49):
    w = dict(worlds[wi])
    conf = PAN(pur=3.0 if wi == 2 else 1.5)
    pp = dict(w); pp.update(conf)
    t, YY, _ = simulate(pp, t_end=60*Y, n_out=4000, timeout=240,
                        atol=1e-8, rtol=1e-6)
    X, Y5, Y4, L, V5, V4, D = YY.T
    print(f"\n=== world {wi} params ===")
    for k in ("a0","d_l","rho","s4","rho4","mu","beta","delta_y","k","u"):
        print(f"   {k}={w[k]:.6e}")
    print("   trajectory (X, y4, l/ml, v4):")
    for yy in [1, 2, 3, 4, 5, 6, 8, 10, 15, 20]:
        i = np.argmin(np.abs(t - yy*Y))
        print(f"    yr{yy:2d}: X={X[i]:6.1f}  y4={Y4[i]:.4e}  y5={Y5[i]:.4e}  "
              f"L/ml={L[i]*1000:.4e}  v4={V4[i]:.4e}  v5={V5[i]:.4e}")
    m = t >= 2*Y
    if np.abs(V4[m]).max() > 1e-3:
        t1 = t[m][np.where(np.abs(V4[m]) > 1e-3)[0][0]]/Y
        print(f"   X4 viremia >1e-3 at yr {t1:.2f}")