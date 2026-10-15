"""Final champion sweep: purge-through-tail protocols, aligned 60 worlds."""
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
        t, YY, _ = simulate(pp, t_end=60*Y, n_out=700, timeout=90,
                            atol=atol, rtol=rtol)
    except BaseException:
        return "ERR"
    V = YY[:,4] + YY[:,5]; m = t >= 2*Y
    ok = bool((np.abs(V[m]).max() <= 1e-3) and (YY[-1,0] >= 500)
              and (np.abs(YY[-50:,3]).max()*1000 <= 1e-6))
    if ok:
        return "C"
    st = "X4" if np.abs(YY[m,5]).max() > np.abs(YY[m,4]).max() else "R5"
    t1 = t[m][np.where(np.abs(V[m]) > 1e-3)[0][0]]/Y
    return "%s@%.1fy" % (st, t1)

def prot(erase, pur, win, art_end, ccr5, refr_win):
    return dict(t_art_on=365, t_art_off=art_end, shocks=[],
                strike_at=[365], strike_erase=[erase], strike_kill=[1.0],
                strike_neutral=[0.999999], strike_refr=[0.9999],
                strike_refr_len=[refr_win], strike_healD=[0.9],
                strike_purge=[pur], strike_cleanup=[1000],
                ccr5=ccr5, ccr5_start=365)

rng = np.random.default_rng(777)
worlds = []
while len(worlds) < 60:
    p = dict(P)
    for k in PERTURB:
        p[k] *= rng.uniform(0.85, 1.15)
    if sane(p):
        worlds.append(p)
print(f"cached {len(worlds)} worlds\n")

cfgs = [
    ("NT            (recap 59/60) ", prot(0.999999999, 1.5, 365, 1095, 0.0, 365)),
    ("B  purge-2y-tail + ccr5.999 ", prot(0.999999999, 1.5, 1095, 1460, 0.999, 1095)),
    ("D  purge3/545d + ART3y +ccr5", prot(0.999999999, 3.0, 545, 1460, 0.999, 545)),
    ("E  purge+lock whole ART     ", prot(0.999999999, 1.5, 1095, 1460, 0.999, 1095)),
]
for name, cf in cfgs:
    res = [eval_world(w, cf) for w in worlds]
    nc = sum(1 for r in res if r == "C")
    bad = [(i, r) for i, r in enumerate(res) if r != "C"]
    print(f"{name}: {nc}/60 cured   failures: {bad}")

print("\n== champion B at two tolerances (flipper census) ==")
champ = cfgs[1][1]
lo = [eval_world(w, champ, atol=1e-8, rtol=1e-6) for w in worlds]
hi = [eval_world(w, champ, atol=1e-6, rtol=1e-4) for w in worlds]
both_fail = [i for i in range(60) if lo[i] != "C" and hi[i] != "C"]
flip = [i for i in range(60) if (lo[i] == "C") != (hi[i] == "C")]
print(f"cured(lo)={sum(1 for r in lo if r=='C')}  cured(hi)={sum(1 for r in hi if r=='C')}")
print(f"DEFINITIVE failures (fail both): {both_fail}")
print(f"tolerance-dependent (flip)     : {flip}")