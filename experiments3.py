"""Experiment suite v3: chronic infection THEN treat, two-strain HIV."""
from model import P, simulate, summarize

Y = 365.25

def run(label, **ov):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    summarize(t, YY, p, label)
    return t, YY, p

print("="*74)
print("A3  UNTREATED control")
print("="*74)
run("A3 untreated", t_art_on=1e9, t_art_off=1e9, shocks=[])

print()
print("="*74)
print("B3  ART only (start d365, stop d1460)")
print("="*74)
run("B3 ART only", t_art_on=365, t_art_off=1460, shocks=[])

print()
print("="*74)
print("C3  ART + 16 shock-and-kill / latency-reversal pulses")
print("="*74)
shk = list(range(430, 1470, 60))
run("C3 ART + LRA x%d" % len(shk), t_art_on=365, t_art_off=1460, shocks=shk)

print()
print("="*74)
print("D3  ART + CCR5 gene therapy 90%  <- now with X4 escape channel")
print("="*74)
run("D3 ART + CCR5 90%", t_art_on=365, t_art_off=1460,
    ccr5=0.90, ccr5_start=500, shocks=[])

print()
print("="*74)
print("E3  ART + CAR-T modest 1.5x  <- realistic power, not magic 5x")
print("="*74)
run("E3 ART + CAR-T 1.5x", t_art_on=365, t_art_off=1460,
    cart_boost=1.5, cart_start=500, shocks=[])

print()
print("="*74)
print("F3  KITCHEN SINK: early ART (d7) + LRA x12 + CCR5 90% + CAR-T 1.5x + bNAb")
print("="*74)
shkF = list(range(90, 1470, 120))
run("F3 early ART + LRA + CCR5 + CAR-T + bNAb",
    t_art_on=7, t_art_off=1460, shocks=shkF,
    ccr5=0.90, ccr5_start=150, cart_boost=1.5, cart_start=150,
    bnab=(400, 0.7, 2.0))