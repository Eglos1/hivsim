"""Figure: finding a cure that requires nothing lifelong."""
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from matplotlib import pyplot as plt

Y = 365.25
plt.rcParams.update({"figure.dpi": 110})
ART_NEVER = 1e9

def path(label, ov, color, lw=2, ls="-", span=10):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=span*Y, n_out=span*365)
    return t/Y, YY.T[0], YY.T[4] + YY.T[5]

scen = [
    ("untreated -> AIDS", dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[]),
     "#c62828", 2.0, "-", 10),
    ("ART only -> rebound d4", dict(t_art_on=365, t_art_off=1460, shocks=[]),
     "#8d6e63", 2.0, "-", 10),
    ("CURE  finite: chronic ART 3y + purge", 
     dict(t_art_on=365, t_art_off=1460, shocks=list(range(430, 1460, 30)), lra=0.20),
     "#00695c", 2.2, "-", 30),
    ("CURE  finite short: early ART 1.5y + purge + CAR-T",
     dict(t_art_on=7, t_art_off=547, shocks=list(range(30, 540, 30)), lra=0.20,
          cart_boost=2.0, cart_start=30, t_cart_end=547),
     "#2e7d32", 2.2, "-", 30),
    ("CURE  one-shot ART-free (chronic)",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          ccr5=0.98, ccr5_start=365, bnab=(365, 0.9, 0.0), bnab_w=30*Y),
     "#4527a0", 2.2, "-", 30),
    ("NEAR MISS  CCR5+60% bNAb (X4 relapse ~yr 16)",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          ccr5=0.98, ccr5_start=365, bnab=(365, 0.6, 0.0), bnab_w=30*Y),
     "#ef6c00", 2.0, ":", 30),
]

fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(11, 7.5))
for lab, ov, col, lw, ls, span in scen:
    tY, X, V = path(lab, ov, col, lw, ls, span)
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
fig.suptitle("Cures that end — and one near-miss (X4 escape at low antibody)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("notlifelong.png")
print("saved notlifelong.png")