"""Experiment suite: poke and prod the HIV model until it stops."""
import numpy as np
from model import P, simulate, summarize, rebound_info, cured_check

YEARS = 365.25

def run(label, **overrides):
    p = dict(P)
    p.update(overrides)
    t, Y, _ = simulate(p, t_end=10*YEARS, n_out=3650)
    summarize(t, Y, p, label)
    return t, Y, p

# ---------------------------------------------------------------- EXPERIMENT A
print("="*72)
print("EXPERIMENT A - do nothing (control)")
print("="*72)
tA, YA, pA = run("A: untreated control",
                 t_art_on=1e9, t_art_off=1e9, shocks=[])

# ---------------------------------------------------------------- EXPERIMENT B
print()
print("="*72)
print("EXPERIMENT B - standard ART, then stop (analytical treatment interruption)")
print("="*72)
tB, YB, pB = run("B: ART days 30-760, then stop", shocks=[])

# ---------------------------------------------------------------- EXPERIMENT C
print()
print("="*72)
print("EXPERIMENT C - ART + shock-and-kill cycles (latency reversal)",)
print("="*72)
pC = dict(P)
pC.update(shocks=[130, 190, 250, 310, 370, 430, 490, 550, 610, 670, 730])
tC, YC, pC = run("C: ART + 11 latency-reversal shocks", shocks=pC["shocks"])

# ---------------------------------------------------------------- EXPERIMENT D
print()
print("="*72)
print("EXPERIMENT D - ART + CCR5 gene-editing protection (90% block)")
print("="*72)
tD, YD, pD = run("D: ART + CCR5 90% protection", ccr5=0.90)

# ---------------------------------------------------------------- EXPERIMENT E
print()
print("="*72)
print("EXPERIMENT E - ART + CAR-T boost + shocks (kill + purge, no CCR5)")
print("="*72)
pE = dict(P)
pE.update(cart_boost=3.0, shocks=[130, 220, 310, 400, 490, 580, 670, 730])
tE, YE, pE = run("E: ART + CAR-T(3x) + 8 shocks", **pE)

# ---------------------------------------------------------------- EXPERIMENT F
print()
print("="*72)
print("EXPERIMENT F - combination test-and-treat EARLY: ART week 1 + shocks + CCR5-bets + CAR-T")
print("="*72)
pF = dict(P)
pF.update(t_art_on=7.0, ccr5=0.90, cart_boost=3.0,
          shocks=[120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 740],
          bnab=(400, 0.70, 2.0))
tF, YF, pF = run("F: early ART + CCR5(90%) + CAR-T(3x) + 12 shocks + bNAb", **pF)