"""Final figure: one-time cures designed by the simulation (50-y horizon)."""
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from matplotlib import pyplot as plt

Y = 365.25
ART_NEVER = 1e9
plt.rcParams.update({"figure.dpi": 110, "font.size": 9})

def run(ov, span=50, n=900):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=span*Y, n_out=n, timeout=30)
    return t/Y, YY

scen = [
    ("untreated -> AIDS", dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[]),
     "#c62828", "-", 2.0),
    ("NEAR-MISS  chronic one-shot, 30-d purge window",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          strike_at=365, strike_erase=0.999, strike_kill=1.0,
          strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=30,
          strike_healD=0.9, strike_purge=1.5, strike_cleanup=10),
     "#ef6c00", ":", 2.2),
    ("CURE  one-time: chronic purge-infusion 99% erase, 90-d window",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          strike_at=365, strike_erase=0.99, strike_kill=1.0,
          strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=90,
          strike_healD=0.9, strike_purge=1.5, strike_cleanup=100),
     "#00897b", "-", 2.2),
    ("CURE  one-time, EARLY (day 14) purge-infusion 90% erase",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          strike_at=14, strike_erase=0.9, strike_kill=1.0,
          strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=30,
          strike_healD=0.9, strike_purge=1.5, strike_cleanup=100),
     "#2e7d32", "-", 2.4),
    ("CURE  one-time: near-total erase + 1y lockdown",
     dict(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
          strike_at=365, strike_erase=0.99999999, strike_kill=1.0,
          strike_neutral=0.999999, strike_refr=0.999, strike_refr_len=365,
          strike_healD=0.9, strike_purge=0.0, strike_cleanup=1.0),
     "#6a1b9a", "-", 2.0),
    ("CURE  finite course (G2, blank + purge 3y, then stop)",
     dict(t_art_on=365, t_art_off=1460,
          shocks=list(range(430, 1460, 30)), lra=0.20),
     "#1565c0", "-", 2.0),
]

fig, (a1, a2) = plt.subplots(2, 1, sharex=True, figsize=(11, 7.5))
for lab, ov, col, ls, lw in scen:
    tY, YY = run(ov)
    a1.plot(tY, YY[:, 0], color=col, lw=lw, ls=ls, label=lab)
    a2.plot(tY, np.maximum(YY[:, 4] + YY[:, 5], 1e-15), color=col, lw=lw, ls=ls)

a1.axhline(200, color="gray", lw=0.8, ls=":")
a1.axhspan(0, 200, color="red", alpha=0.05)
a2.axhline(50, color="gray", lw=0.8, ls=":")
a2.axhline(1e-3, color="gray", lw=0.8, ls=":")
a2.text(48.1, 1.5, "flat-zero floor", fontsize=7, color="gray")
a1.set_ylim(0, 1000); a1.set_ylabel("CD4+ cells / uL")
a2.set_yscale("log"); a2.set_ylim(1e-15, 1e7); a2.set_ylabel("viral load copies/mL")
for ax in (a1, a2):
    ax.legend(fontsize=7.5, loc="center right", ncol=2)
    ax.grid(alpha=0.3)
a2.set_xlabel("years")
fig.suptitle("Designing a ONE-TIME cure (no lifelong therapy) - 50-year horizon",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("onetime.png")
print("saved onetime.png")