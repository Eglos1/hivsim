"""Frontier (fast): one pass per world, judge STRICT + HUMAN endpoints,
stream results as they compute."""
import warnings, sys; warnings.filterwarnings("ignore")
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

PAN = dict(t_art_on=365, t_art_off=1460, shocks=[],
           strike_at=[365], strike_erase=[0.999999999], strike_kill=[1.0],
           strike_neutral=[0.999999], strike_refr=[0.9999],
           strike_refr_len=[1095], strike_healD=[0.9],
           strike_purge=[1.5], strike_cleanup=[1000],
           ccr5=0.999, ccr5_start=365)

def run(w, atol):
    pp = dict(w); pp.update(PAN)
    try:
        t, YY, _ = simulate(pp, t_end=60*Y, n_out=500, timeout=60,
                            atol=atol, rtol=1e-6)
        return t, YY
    except BaseException:
        return None, None

def judge(t, YY):
    if t is None or not np.isfinite(YY).all():
        return "ERR", "ERR"
    V = YY[:,4] + YY[:,5]
    V5, V4 = YY[:,4], YY[:,5]
    xs = YY[:,0]
    mS = t >= 2*Y
    strict = (np.abs(V[mS]).max() <= 1e-3 and YY[-1,0] >= 500
              and np.abs(YY[-50:,3]).max()*1000 <= 1e-6)
    mH = t >= 3*Y
    human = (V[mH].max() < 50 and ys_min() >= 0 and YY[-1,0] >= 500)
    # dominant strain of any relapse blip
    st = "-"
    if np.abs(V[mS]).max() > 1e-3:
        st = "R5" if np.abs(V5[mS]).max() >= np.abs(V4[mS]).max() else "X4"
    return ("C" if strict else "R") + st, ("C" if human else "R") + st

def ys_min():
    return 0

def collect(seed):
    rng = np.random.default_rng(seed)
    worlds = []
    while len(worlds) < 60:
        p = dict(P)
        for k in PERTURB:
            p[k] *= rng.uniform(0.85, 1.15)
        if sane(p):
            worlds.append(p)
    return worlds

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 777
worlds = collect(seed)
print(f"seed {seed}: {len(worlds)} worlds cached", flush=True)
per_world = []
for i, w in enumerate(worlds):
    tL, YYL = run(w, 1e-8)
    tH, YYH = run(w, 1e-6)
    sL_s, sL_h = judge(tL, YYL)
    sH_s, sH_h = judge(tH, YYH)
    per_world.append((sL_s, sH_s, sL_h, sH_h))
    if i % 10 == 9:
        print(f"  ...{i+1}/60 done", flush=True)

def census(idx_lo, idx_hi):
    bc = sum(1 for x in per_world if x[idx_lo].startswith("C") and x[idx_hi].startswith("C"))
    br = [i for i,x in enumerate(per_world) if not x[idx_lo].startswith("C") and not x[idx_hi].startswith("C")]
    fl = [i for i,x in enumerate(per_world) if x[idx_lo].startswith("C") != x[idx_hi].startswith("C")]
    return bc, br, fl

for name, lo_i, hi_i in [("STRICT",0,1), ("HUMAN",2,3)]:
    bc, br, fl = census(lo_i, hi_i)
    print(f"[{name}] cured_both={bc}/60  relapsed_both={br}  flips={fl}", flush=True)

print("knife-edge biology (STRICT flips):", flush=True)
for (sL_s, sH_s, *_ ), w in zip(per_world, worlds):
    if sL_s.startswith("C") != sH_s.startswith("C"):
        slack = w["a0"] + w["d_l"] - w["rho"]
        gain5 = (1-w["fL"])*w["beta"]*w["lam"]/w["d"]*w["k"]/w["u"]/w["delta_y"]
        gain4 = (1-w["fL"])*w["beta"]*w["rho4"]*w["lam"]/w["d"]*w["k"]/w["u"]/w["delta_y"]
        print(f"   world slack_res={slack:.3e} R5gain={gain5:.2f} X4gain={gain4:.2f} "
              f"s4={w['s4']:.2e} d_l={w['d_l']:.3e} rho={w['rho']:.3e} "
              f"(lo={sL_s} hi={sH_s})", flush=True)