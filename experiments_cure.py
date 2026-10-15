"""The search: what actually stops HIV in the simulation?"""
from model import P, simulate, summarize, rebound_info, cured_check

Y = 365.25

def run(label, **ov):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    summarize(t, YY, p, label)
    return t, YY, p

print("="*74)
print("G1  Sterilizing route: EARLY ART + maximal purge")
print("    ART during acute (d7), 47 latency-reversal pulses lra=0.20")
print("="*74)
shockG1 = list(range(60, 1470, 30))
run("G1 early ART + MAX purge", t_art_on=7, t_art_off=1460,
    shocks=shockG1, lra=0.20)

print()
print("="*74)
print("G2  Sterilizing route attempted LATE (chronic infection)")
print("    ART d365, 35 pulses lra=0.20 over 3 ART years")
print("="*74)
shockG2 = list(range(430, 1460, 30))
run("G2 chronic ART + MAX purge", t_art_on=365, t_art_off=1460,
    shocks=shockG2, lra=0.20)

print()
print("="*74)
print("H1  Functional route: 'cell & gene factory' combo")
print("    early ART + purge + CCR5 + CAR-T + PERMANENT bNAb factory")
print("="*74)
run("H1 early ART + factory combo",
    t_art_on=7, t_art_off=1460,
    shocks=list(range(90, 1470, 90)), lra=0.12,
    ccr5=0.90, ccr5_start=150,
    cart_boost=1.5, cart_start=150,
    bnab=(200, 0.7, 2.0), bnab_w=15*Y)

print()
print("="*74)
print("H2  Functional route WITHOUT permanent bNAb (all finite)")
print("    same as H1 but bNAb window only 120 days")
print("="*74)
run("H2 finite antibody course",
    t_art_on=7, t_art_off=1460,
    shocks=list(range(90, 1470, 90)), lra=0.12,
    ccr5=0.90, ccr5_start=150,
    cart_boost=1.5, cart_start=150,
    bnab=(200, 0.7, 2.0), bnab_w=120)

print()
print("="*74)
print("SWEEP  purge intensity vs rebound latency (early ART)")
print("       early ART d7, lra=0.12, N latency-reversal pulses")
print("="*74)
for n in [0, 4, 8, 16, 24, 32]:
    sh = list(range(60, 1470, max(1, int(1410/max(n,1))))) if n else []
    p = dict(P); p.update(t_art_on=7, t_art_off=1460, shocks=sh)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    rb = rebound_info(t, YY, p)
    pct = YY.T[4][-1]  # v4 at end
    tot = YY.T[4][-1] + YY.T[5][-1]
    print(f"  pulses={n:>2}: rebound = "
          f"{'NONE within 10y (cure)' if rb is None else f'{rb:5.0f} days after ART stop'}"
          f"   | end X4% = {100*pct/max(tot,1e-9):5.1f}")

print()
print("="*74)
print("SWEEP2 CCR5 protection alone at chronic stage: how far can it go?")
print("       CCR5 start d500, ART stop d1460, shocks none")
print("="*74)
for eff in [0.5, 0.7, 0.8, 0.9, 0.98]:
    p = dict(P); p.update(t_art_on=365, t_art_off=1460, ccr5=eff, ccr5_start=500)
    t, YY, _ = simulate(p, t_end=10*Y, n_out=3650)
    rb = rebound_info(t, YY, p)
    v4 = YY.T[5][-1]; vt = YY.T[4][-1]+YY.T[5][-1]
    print(f"  ccr5={eff:5.2f}: rebound = "
          f"{'NONE' if rb is None else f'{rb:5.0f} days'}"
          f"   | end VL={vt:10.0f}  X4%={100*v4/max(vt,1e-9):5.1f}")