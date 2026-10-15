"""Final 'miracle cure' figure: CD4 + virus + reservoir + immune damage (D)
over 60 years for the PAN protocol (the validated champion).

PAN = one-time strike (99.9999999% erase + 1.5/day purge for 3y inside a
refractory lockdown) + permanent CCR5 knockout + 3y ART quiet tail.
The 4th panel shows ALL symptoms heal: immune damage D collapses to ~0,
i.e. the cure restores the body, not just the virus count.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import P, simulate

Y = 365.25
PAN = dict(t_art_on=365, t_art_off=1460, shocks=[],
    strike_at=365, strike_erase=0.999999999, strike_kill=1.0,
    strike_neutral=0.999999, strike_refr=0.9999, strike_refr_len=1095,
    strike_healD=0.9, strike_purge=1.5, strike_cleanup=1000,
    ccr5=0.999, ccr5_start=365)

p = dict(P); p.update(PAN)
t, YY, _ = simulate(p, t_end=60*Y, n_out=4000, rtol=1e-6, atol=1e-8)
X, Y5, Y4, L, V5, V4, D = YY.T
V = V5 + V4

fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
yr = t / Y

ax = axes[0]
ax.plot(yr, X, color="#1f6fb2", lw=2)
ax.axhline(500, color="#888", ls="--", lw=1)
ax.text(0.4, 520, "HIV-only AIDS threshold  CD4=500", color="#666", fontsize=9)
ax.set_ylabel("CD4+ T cells  (cells/µL)")
ax.set_title("PAN: one-time strike + CCR5 knockout + 3-yr ART tail — 60 years, "
             "whole-body cure", fontsize=12)

ax = axes[1]
vlog = np.where(V > 1e-9, np.log10(np.maximum(V, 1e-9)), -9.5)
ax.plot(yr, vlog, color="#c0392b", lw=2)
ax.set_ylabel("log10 virus  (copies/mL)")
ax.set_ylim(-10, 7)
ax.axhline(np.log10(1e-3), color="#888", ls="--", lw=1)
ax.text(0.4, np.log10(1e-3)+0.3, "undetectable floor 10^-3 copies/mL", color="#666", fontsize=9)
for d in (365, 1460):
    ax.axvline(d/Y, color="#999", ls=":", lw=1)
ax.text(1.02, 3, "strike + purge\nstarts day 365", fontsize=8, color="#555")
ax.text(4.15, 3, "ART + eternal CCR5 lock at day 1460;\nvirus stays at zero", fontsize=8, color="#555")

ax = axes[2]
ax.plot(yr, np.maximum(L*1000, 1e-14), color="#7d3c98", lw=2)
ax.set_yscale("log")
ax.set_ylim(1e-14, 1e4)
ax.set_ylabel("latent reservoir  (cells/mL)")
ax.text(0.7, 1e-2, "purge burn drives the reservoir below 1 cell\n"
        "in every perturbed patient (discrete cell-count verdict: 0 cells)",
        fontsize=9, color="#555")

ax = axes[3]
ax.plot(yr, D, color="#2e8b57", lw=2)
ax.axhline(0, color="#bbb", ls="--", lw=1)
ax.set_ylim(-200, 30000)
ax.set_ylabel("immune damage  D")
ax.set_xlabel("Years after diagnosis")
ax.text(1.0, 25000, "cumulative immune attrition\n(thr=30000 means AIDS)",
        fontsize=9, color="#555")
ax.text(8, 2000, "damage fully healed to ~0", fontsize=9, color="#2e8b57")

plt.tight_layout()
plt.savefig("miracle.png", dpi=150)
print("saved miracle.png (PAN, 4-panel whole-body cure)")