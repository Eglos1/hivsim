"""Champion B validation on two fresh 60-world batches + final miracle figure."""
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

CHAMP = dict(t_art_on=365, t_art_off=1460, shocks=[],
             strike_at=[365], strike_erase=[0.999999999], strike_kill=[1.0],
             strike_neutral=[0.999999], strike_refr=[0.9999],
             strike_refr_len=[1095], strike_healD=[0.9],
             strike_purge=[1.5], strike_cleanup=[1000],
             ccr5=0.999, ccr5_start=365)

def evalu(p, atol, rtol):
    pp = dict(p); pp.update(CHAMP)
    try:
        t, YY, _ = simulate(pp, t_end=60*Y, n_out=700, timeout=90,
                            atol=atol, rtol=rtol)
    except BaseException:
        return "ERR"
    V = YY[:,4] + YY[:,5]; m = t >= 2*Y
    ok = bool((np.abs(V[m]).max() <= 1e-3) and (YY[-1,0] >= 500)
              and (np.abs(YY[-50:,3]).max()*1000 <= 1e-6))
    return "C" if ok else "R"

for seed in (777, 4242):
    rng = np.random.default_rng(seed)
    worlds = []
    while len(worlds) < 60:
        p = dict(P)
        for k in PERTURB:
            p[k] *= rng.uniform(0.85, 1.15)
        if sane(p):
            worlds.append(p)
    lo = [evalu(w, 1e-8, 1e-6) for w in worlds]
    hi = [evalu(w, 1e-6, 1e-4) for w in worlds]
    dep = [i for i in range(60) if (lo[i] == "C") != (hi[i] == "C")]
    both_fail = [i for i in range(60) if lo[i] == "R" and hi[i] == "R"]
    both_cure = sum(1 for i in range(60) if lo[i] == "C" and hi[i] == "C")
    print(f"seed {seed}:   cured_both={both_cure}/60   "
          f"definitive_failures={both_fail}   undecidable={len(dep)}")

# final figure for a representative world (nominal P)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
p = dict(P); p.update(CHAMP)
t, YY, _ = simulate(p, t_end=60*Y, n_out=4000, rtol=1e-6, atol=1e-8)
X, Y5, Y4, L, V5, V4, D = YY.T; V = V5 + V4
fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
yr = t / Y
ax = axes[0]; ax.plot(yr, X, color="#1f6fb2", lw=2)
ax.axhline(500, color="#888", ls="--", lw=1)
ax.text(0.5, 510, "AIDS threshold  CD4=500", color="#666", fontsize=9)
ax.set_ylabel("CD4+ T cells (cells/µL)")
ax.set_title("THE MIRACLE CURE 'PAN' — one 3-year course, 60 years clean")
ax = axes[1]
vlog = np.where(V > 1e-9, np.log10(np.maximum(V, 1e-9)), -9.5)
ax.plot(yr, vlog, color="#c0392b", lw=2)
ax.set_ylabel("log10 virus (copies/mL)")
ax.set_ylim(-10, 7)
ax.axhline(np.log10(1e-3), color="#888", ls="--", lw=1)
ax.text(0.5, np.log10(1e-3)+0.3, "undetectable floor 10^-3", color="#666", fontsize=9)
for d in (365, 1460):
    ax.axvline(d/Y, color="#999", ls=":", lw=1)
ax.text(1.05, 3, "course starts\nday 365", fontsize=8, color="#555")
ax.text(4.2, 3, "course ends\nday 1460", fontsize=8, color="#555")
ax = axes[2]
ax.plot(yr, np.maximum(L*1000, 1e-16), color="#7d3c98", lw=2)
ax.set_yscale("log"); ax.set_ylim(1e-16, 1e4)
ax.set_ylabel("latent reservoir (cells/mL)")
ax.set_xlabel("Years after diagnosis")
plt.tight_layout(); plt.savefig("miracle.png", dpi=150)
print("saved miracle.png (champion PAN curve)")
t_check = t >= 2*Y
print("nominal world: flat-zero after yr2:",
      bool((np.abs(V[t_check]).max() <= 1e-3)),
      " CD4_60=%.0f  L60/ml=%.3e" % (X[-1], L[-1]*1000))