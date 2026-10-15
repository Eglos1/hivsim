"""Does 'block-and-lock' (permanent X4 silencing) close the last wall?
PAN + silence4 from various starts, judged STRICT + HUMAN, two tolerances."""
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

def PAN(sil4=0.0, sil4_start=None, sil5=0.0):
    return dict(t_art_on=365, t_art_off=1460, shocks=[],
                strike_at=[365], strike_erase=[0.999999999],
                strike_kill=[1.0], strike_neutral=[0.999999],
                strike_refr=[0.9999], strike_refr_len=[1095],
                strike_healD=[0.9], strike_purge=[1.5],
                strike_cleanup=[1000], ccr5=0.999, ccr5_start=365,
                silence4=sil4, silence5=sil5,
                silence_start=sil4_start if sil4_start is not None else 1e9)

def run(w, conf, atol):
    pp = dict(w); pp.update(conf)
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
    mS = t >= 2*Y
    strict = (np.abs(V[mS]).max() <= 1e-3 and YY[-1,0] >= 500
              and np.abs(YY[-50:,3]).max()*1000 <= 1e-6)
    mH = t >= 3*Y
    human = (V[mH].max() < 50 and YY[-1,0] >= 500)
    return ("C" if strict else "R"), ("C" if human else "R")

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

cfgs = [
    ("PAN alone (recap)                ", PAN()),
    ("PAN + X4-lock@1460 (0.999999999)", PAN(sil4=0.999999999, sil4_start=1460)),
    ("PAN + X4-lock@1095              ", PAN(sil4=0.999999999, sil4_start=1095)),
    ("PAN + X4-lock@365               ", PAN(sil4=0.999999999, sil4_start=365)),
    ("PAN + both-lock@1460            ", PAN(sil4=0.999999999, sil4_start=1460,
                                             sil5=0.999999999)),
]
for seed in (777, 4242):
    worlds = collect(seed)
    print(f"\n===== seed {seed} =====", flush=True)
    for name, cf in cfgs:
        out = []
        for w in worlds:
            tL, YYL = run(w, cf, 1e-8)
            tH, YYH = run(w, cf, 1e-6)
            out.append((judge(tL, YYL), judge(tH, YYH)))
        # census (STRICT index 0 of each judge tuple)
        bc = sum(1 for o in out if o[0][0] == "C" and o[1][0] == "C")
        br = sorted(i for i, o in enumerate(out)
                    if o[0][0] == "R" and o[1][0] == "R")
        fl = sorted(i for i, o in enumerate(out)
                    if o[0][0] != o[1][0])
        hu = sum(1 for o in out if o[0][1] == "C" and o[1][1] == "C")
        print(f"  {name}: STRICT both-cured={bc}/60 relapsers={br} flips={fl} "
              f"| HUMAN both-cured={hu}/60", flush=True)