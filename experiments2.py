"""Experiment suite v2: establish chronic infection, THEN treat."""
from model import P, simulate, summarize, rebound_info, cured_check

Y = 365.25

def run(label, **ov):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    summarize(t, YY, p, label)
    return t, YY, p

# Every scenario: infection day 0, untreated until ART starts day 365 (chronic),
# ART for 3 years, then interrupted at day 1460.

print("="*72)
print("EXPERIMENT A2 - untreated control (re-baseline)")
print("="*72)
run("A2: untreated", t_art_on=1e9, t_art_off=1e9, shocks=[])

print()
print("="*72)
print("EXPERIMENT B2 - ART alone at chronic stage, then interrupt")
print("="*72)
run("B2: ART start d365 stop d1460", t_art_on=365, t_art_off=1460, shocks=[])

print()
print("="*72)
print("EXPERIMENT C2 - ART + shock-and-kill (latency reversal). lra=0.12/day,")
print("               12 shocks every ~60d across the 3 ART years")
print("="*72)
shk = [430, 495, 560, 625, 690, 755, 820, 885, 950, 1015, 1080, 1145, 1210, 1275, 1340, 1405]
run("C2: ART + 16 LRA shocks", t_art_on=365, t_art_off=1460, shocks=shk, lra=0.12)

print()
print("="*72)
print("EXPERIMENT D2 - ART + CCR5 gene therapy (90%) at chronic stage")
print("="*72)
run("D2: ART + CCR5 90% (d500)", t_art_on=365, t_art_off=1460,
    ccr5=0.90, ccr5_start=500, shocks=[])

print()
print("="*72)
print("EXPERIMENT E2 - ART + CAR-T (5x killing boost) at chronic stage")
print("="*72)
run("E2: ART + CAR-T 5x (d500)", t_art_on=365, t_art_off=1460,
    cart_boost=5.0, cart_start=500, shocks=[])

print()
print("="*72)
print("EXPERIMENT F2 - early test-and-treat (ART at week 1) + shocks + CCR5")
print("="*72)
shkF = [90, 150, 210, 270, 330, 390, 450, 510, 570, 630, 690, 750]
run("F2: early ART (d7) + 12 LRA + CCR5 90% (d120)",
    t_art_on=7, t_art_off=1460, shocks=shkF,
    ccr5=0.90, ccr5_start=120, lra=0.12)