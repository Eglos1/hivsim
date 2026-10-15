"""Even-better hunting: aligned 60-world comparison for permanent-CCR5 variants."""
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

def eval_world(p, conf, atol=1e-6, rtol=1e-4):
    pp = dict(p); pp.update(conf)
    try:
        t, YY, _ = simulate(pp, t_end=60*Y, n_out=600, timeout=40,
                            atol=atol, rtol=rtol)
    except BaseException:
        return "ERR"
    V5, V4 = YY[:,4], YY[:,5]; V = V5 + V4; m = t >= 2*Y
    ok = bool((np.abs(V[m]).max() <= 1e-3) and (YY[-1,0] >= 500)
              and (np.abs(YY[-50:,3]).max()*1000 <= 1e-6))
    if ok:
        return "C"
    st = "R5" if np.abs(V5[m]).max() >= np.abs(V4[m]).max() else "X4"
    # failing fast (relapse before yr 10) or slow?
    if t[m][np.where(np.abs(V[m]) > 1e-3)[0][0]]/Y < 10:
        return st + "-fast"
    return st + "-slow"

def NT(ccr5=0.0, pur=1.5, win=365, art_end=1095, erase=0.999999999,
       extra=None):
    d = dict(t_art_on=365, t_art_off=art_end, shocks=[],
             strike_at=[365], strike_erase=[erase], strike_kill=[1.0],
             strike_neutral=[0.999999], strike_refr=[0.9999],
             strike_refr_len=[win], strike_healD=[0.9],
             strike_purge=[pur], strike_cleanup=[1000],
             ccr5=ccr5, ccr5_start=365)
    if extra:
        d.update(extra)
    return d

# second-strike variant (both strikes effective, inside the ART tail)
dbl = dict(strike_at=[365, 730], strike_erase=[0.999999999, 0.99],
           strike_purge=[1.5, 2.0], strike_refr_len=[365, 180],
           strike_refr=[0.9999, 0.9999], strike_kill=[1.0, 1.0],
           strike_neutral=[0.999999, 0.999999],
           strike_cleanup=[1000, 1000], strike_healD=[0.9, 0.9])

cfgs = [
    ("NT(ccr5=0)            ", NT()),
    ("NT + ccr5 99%         ", NT(0.99)),
    ("NT + ccr5 99.9%       ", NT(0.999)),
    ("NT + ccr5 99.99%      ", NT(0.99999)),
    ("NT + ccr5.999 + pur2.0", NT(0.999, pur=2.0)),
    ("NT + ccr5.999 + dbl@730", NT(0.999, extra=dbl)),
]

rng = np.random.default_rng(777)
worlds = []
while len(worlds) < 60:
    p = dict(P)
    for k in PERTURB:
        p[k] *= rng.uniform(0.85, 1.15)
    if sane(p):
        worlds.append(p)
print(f"cached {len(worlds)} worlds\n")

for name, cf in cfgs:
    res = [eval_world(w, cf) for w in worlds]
    nc = sum(1 for r in res if r == "C")
    nerr = sum(1 for r in res if r == "ERR")
    bad = [(i, r) for i, r in enumerate(res) if r != "C"]
    print(f"{name}: {nc}/60 cured (ERR={nerr})  failures: {bad}")