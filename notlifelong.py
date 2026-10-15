"""Hunting for a cure that requires nothing lifelong:
  - completely ART-FREE single-shot cell & gene therapy
  - shortest possible ART course (acute and chronic starts)
"""
from model import P, simulate, rebound_info, sustained_remission

Y = 365.25
ART_NEVER = 1e9

def run(label, **ov):
    p = dict(P); p.update(ov)
    t, YY, _ = simulate(p, t_end=12*Y, n_out=4380)
    X = YY.T[0]; V = YY.T[4] + YY.T[5]
    L = YY.T[3]
    cured = sustained_remission(t, YY)
    rb = rebound_info(t, YY, p)
    rbstr = ("REMAINS suppressed (no ART to stop)" if rb is None
             else f"rebound {rb:.0f}d off ART")
    status = " CURE (no lifelong anything)" if cured else " NOT cured"
    print(f"[{'*' if cured else ' '}] {label}")
    print(f"      CD4_final={X[-1]:6.0f}  VL_final={V[-1]:10.0f}  "
          f"reservoir_final={L[-1]*1000:9.4f}/mL  -> {rbstr}"
          f"{status}")
    return t, YY, p

print("================ SANITY: helpers on known scenarios ================")
run("sanity untreated (AIDS by yr ~9)", t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[])
# quick changed signature
p0 = dict(P)
t0_, Y0_, _ = simulate(p0, t_end=12*Y, n_out=4380)
print("sanity G1 early+max purge:", sustained_remission(t0_, Y0_), flush=True)
print()
print("=================== ONE-SHOT ART-FREE EXPLORATION ===================")
print("(therapy applied once at chronic stage day 365; NO ART ever)")
run("M1  CCR5 98% + permanent bNAb 90% + CAR-T 2x   (day 365)",
    t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
    ccr5=0.98, ccr5_start=365, cart_boost=2.0, cart_start=365,
    bnab=(365, 0.9, 0.0), bnab_w=40*Y)
run("M1b CCR5 98% + permanent bNAb 90%  (day 365, no CAR-T)",
    t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
    ccr5=0.98, ccr5_start=365,
    bnab=(365, 0.9, 0.0), bnab_w=40*Y)
run("M2  CCR5 98% + WEAK bNAb 60% (day 365) -- borderline",
    t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
    ccr5=0.98, ccr5_start=365,
    bnab=(365, 0.60, 0.0), bnab_w=40*Y)
run("M3  one-shot at ACUTE (day 14): CCR5 90% + bNAb 70% + CAR-T 1.5x",
    t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
    ccr5=0.90, ccr5_start=14, cart_boost=1.5, cart_start=14,
    bnab=(14, 0.7, 2.0), bnab_w=40*Y)

print()
print("=============== SHORTEST ART COURSE (early treatment) ===============")
print("ART starts day 7; purge pulses every 60d during ART; stop drug after D days")
for D in [60, 90, 120, 180, 270, 365]:
    sh = [d for d in range(45, D-5, 60)]
    run(f"N early ART {D}d + {len(sh)} purge pulses",
        t_art_on=7, t_art_off=7+D, shocks=sh, lra=0.20)

print()
print("=========== SHORTEST COURSE STARTED AT CHRONIC (day 365) ===========")
for D in [90, 180, 365, 720]:
    sh = [d for d in range(400, 365+D-5, 60)]
    run(f"O chronic ART {D}d + {len(sh)} purge pulses",
        t_art_on=365, t_art_off=365+D, shocks=sh, lra=0.20)