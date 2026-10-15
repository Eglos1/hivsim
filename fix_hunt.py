"""Autopsy + fix hunt (fast): seed 777 worlds 2, 15, 49, 16 under fixes."""
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

def PAN(pur=1.5, win=1095, art_end=1460, erase=0.999999999, crush=None):
    s_at, s_er, s_ki, s_ne, s_re, s_rl, s_he, s_pu, s_cl = \
        [365], [erase], [1.0], [0.999999], [0.9999], [win], [0.9], [pur], [1000]
    if crush is not None:
        s_at.append(1095); s_er.append(0.0); s_ki.append(0.0); s_ne.append(0.0)
        s_re.append(0.0);  s_rl.append(365); s_he.append(0.0)
        s_pu.append(max(crush - pur, 0.0)); s_cl.append(1.0)
    return dict(t_art_on=365, t_art_off=art_end, shocks=[],
                strike_at=s_at, strike_erase=s_er, strike_kill=s_ki,
                strike_neutral=s_ne, strike_refr=s_re, strike_refr_len=s_rl,
                strike_healD=s_he, strike_purge=s_pu, strike_cleanup=s_cl,
                ccr5=0.999, ccr5_start=365)

def run(p, conf):
    pp = dict(p); pp.update(conf)
    try:
        t, YY, _ = simulate(pp, t_end=60*Y, n_out=600, timeout=60,
                            atol=1e-8, rtol=1e-6)
        return t, YY
    except BaseException:
        return None, None

def diagnose(t, YY):
    if t is None or not np.isfinite(YY).all():
        return "ERR"
    V = YY[:,4] + YY[:,5]; V5 = YY[:,4]; V4 = YY[:,5]; m = t >= 2*Y
    if np.abs(V[m]).max() <= 1e-3:
        return "CURED"
    st = "R5" if np.abs(V5[m]).max() >= np.abs(V4[m]).max() else "X4"
    t1 = t[m][np.where(np.abs(V[m]) > 1e-3)[0][0]]/Y
    return f"R@{t1:.1f}y/{st}"

worlds = []
rng = np.random.default_rng(777)
while len(worlds) < 60:
    p = dict(P)
    for k in PERTURB:
        p[k] *= rng.uniform(0.85, 1.15)
    if sane(p):
        worlds.append(p)
T = {i: worlds[i] for i in (2, 15, 49, 16)}

fixes = [
    ("baseline PAN       ", PAN()),
    ("A purge3.0 whole    ", PAN(pur=3.0)),
    ("C final crush to 5  ", PAN(crush=5.0)),
    ("E erase 1e-11       ", PAN(erase=0.99999999999)),
    ("F erase1e-11+purge3 ", PAN(erase=0.99999999999, pur=3.0)),
    ("F2 erase1e-11+crush5", PAN(erase=0.99999999999, crush=5.0))]
for name, cf in fixes:
    cells = []
    for i in (2, 15, 49, 16):
        cells.append("%s:%s" % (i, diagnose(*run(T[i], cf))))
    print(name, cells, flush=True)