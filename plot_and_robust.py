"""Robustness check + summary figure (optimized)."""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate

Y = 365.25
rng = np.random.default_rng(7)

def run_cure_fn(name):
    if name == "G1":
        return dict(t_art_on=7, t_art_off=1460,
                    shocks=list(range(60, 1470, 30)), lra=0.20)
    if name == "G2":
        return dict(t_art_on=365, t_art_off=1460,
                    shocks=list(range(430, 1460, 30)), lra=0.20)
    if name == "H1":
        return dict(t_art_on=7, t_art_off=1460,
                    shocks=list(range(90, 1470, 90)), lra=0.12,
                    ccr5=0.90, ccr5_start=150,
                    cart_boost=1.5, cart_start=150,
                    bnab=(200, 0.7, 2.0), bnab_w=15*Y)
    raise KeyError(name)

PERTURB = ["lam", "d", "beta", "fL", "a0", "d_l", "rho", "delta_y",
           "k", "u", "rho4", "mu", "s4", "cD", "rD", "D_thr"]

def sane_world(p):
    """Only keep parameter worlds that produce a plausible untreated course."""
    pp = dict(p); pp.update(t_art_on=1e9, t_art_off=1e9, shocks=[])
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=250, timeout=3.0)
    except BaseException:
        return False
    X = YY.T[0]; V = YY.T[4] + YY.T[5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def cured_3y_off_art(p, t, YY):
    stop = p["t_art_off"]
    tail = (t >= stop + 3*Y) & (t <= t[-1])
    v = YY.T[4] + YY.T[5]
    return bool(tail.any() and (v[tail] < 50).all())

print("Robustness: perturb 16 params by +/-15%; test cure strategies", flush=True)
for name in ["G1", "G2", "H1"]:
    base = run_cure_fn(name)
    passes = skipped = tries = 0
    while (passes + skipped) < 25 and tries < 40:
        p = dict(P)
        for key in PERTURB:
            p[key] = p[key] * rng.uniform(0.85, 1.15)
        tries += 1
        if not sane_world(p):
            skipped += 1
            continue
        p.update(base)
        try:
            t, YY, _ = simulate(p, t_end=8.5*Y, n_out=1100, timeout=3.0)
        except BaseException:
            skipped += 1
            continue
        passes += int(cured_3y_off_art(p, t, YY))
    print(f"  {name}: off-ART cure in {passes}/25 tested worlds "
          f"({100*passes/25:.0f}%)  [worlds screened out: {skipped}]", flush=True)

print("\nContrast - simple ART only, same worlds:", flush=True)
passes = skipped = tries = 0
while (passes + skipped) < 25 and tries < 40:
    p = dict(P)
    for key in PERTURB:
        p[key] = p[key] * rng.uniform(0.85, 1.15)
    tries += 1
    if not sane_world(p):
        skipped += 1; continue
    p.update(t_art_on=365, t_art_off=1460, shocks=[])
    try:
        t, YY, _ = simulate(p, t_end=8.5*Y, n_out=1100, timeout=3.0)
    except BaseException:
        skipped += 1; continue
    passes += int(cured_3y_off_art(p, t, YY))
print(f"  ART only: off-ART cure in {passes}/25 tested worlds "
      f"({100*passes/25:.0f}%)", flush=True)

# ------------------------------------------------------------------ figure
from matplotlib import pyplot as plt
plt.rcParams.update({"figure.dpi": 110})

def path(label, ov, color, lw=2, ls="-"):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    return t/Y, YY.T[0], YY.T[4] + YY.T[5]

scen = [
    ("untreated -> AIDS by ~yr 9",
     dict(t_art_on=1e9, t_art_off=1e9, shocks=[]), "#c62828", 2.2, "-"),
    ("ART only -> rebound day 4",
     dict(t_art_on=365, t_art_off=1460, shocks=[]), "#1565c0", 2.0, "-"),
    ("ART + LRA purge -> rebound day 10",
     dict(t_art_on=365, t_art_off=1460, shocks=list(range(430, 1470, 60))),
     "#f9a825", 2.0, "-."),
    ("CCR5 mono-therapy -> X4 escape",
     dict(t_art_on=365, t_art_off=1460, ccr5=0.9, ccr5_start=500),
     "#8e24aa", 2.0, ":"),
    ("CURE: early ART + max purge",
     dict(t_art_on=7, t_art_off=1460, shocks=list(range(60, 1470, 30)), lra=0.20),
     "#2e7d32", 2.2, "-"),
    ("CURE: chronic ART + max purge",
     dict(t_art_on=365, t_art_off=1460, shocks=list(range(430, 1460, 30)), lra=0.20),
     "#00695c", 2.0, "--"),
]
t0 = time.time()
fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(10, 7))
for lab, ov, col, lw, ls in scen:
    tY, X, V = path(lab, ov, col, lw, ls)
    a1.plot(tY, X, color=col, lw=lw, ls=ls, label=lab)
    a2.plot(tY, np.maximum(V, 1), color=col, lw=lw, ls=ls, label=lab)
a1.axhline(200, color="gray", lw=0.8, ls=":")
a1.axhspan(0, 200, color="red", alpha=0.05)
a2.axhline(50, color="gray", lw=0.8, ls=":")
a1.set_ylim(0, 1000); a1.set_ylabel("CD4+ T cells / uL")
a2.set_yscale("log"); a2.set_ylim(1, 1e7); a2.set_ylabel("viral load copies/mL")
for ax in (a1, a2):
    ax.legend(fontsize=8, loc="center right"); ax.grid(alpha=0.3)
a2.set_xlabel("years")
fig.suptitle("Poking at HIV until it stops: strategy outcomes", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("outcomes.png")
print(f"\nsaved outcomes.png  ({time.time()-t0:.0f}s)", flush=True)