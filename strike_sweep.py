"""Design-space search: the minimal ONE-TIME STRIKE that durably cures
(chronic & early infection, no ART ever, 30-yr horizon)."""
import warnings
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate, sustained_remission

Y = 365.25
ART_NEVER = 1e9

def strike(stage, erase, refr, refr_len, kill=1.0, neutral=0.999999,
           healD=0.9, boostx=1.2):
    p = dict(P)
    day = {"early": 14, "chronic": 365}[stage]
    p.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER, shocks=[],
             strike_at=day, strike_erase=erase, strike_kill=kill,
             strike_neutral=neutral, strike_refr=refr, strike_refr_len=refr_len,
             strike_healD=healD, strike_boostx=boostx)
    return p

def cured30(stage, erase, refr, refr_len, timeout=4.0):
    p = strike(stage, erase, refr, refr_len)
    t, YY, _ = simulate(p, t_end=30*Y, n_out=1500, timeout=timeout)
    return sustained_remission(t, YY)

# ---------------- 1) any cure at refr 0.999 with moderate/long lockdown? ----
print("== chronic, day-365 strike | erase sweep (lockdown: refr .999, len 365d) ==")
for erase in [0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999]:
    ok = cured30("chronic", erase, 0.999, 365)
    print(f"   erase {erase:9.8f}: {'CURED' if ok else 'rebound'}")

print("\n== chronic, day-365 strike | longer lockdown 730d ==")
for erase in [0.99, 0.999, 0.9999, 0.99999, 0.999999]:
    ok = cured30("chronic", erase, 0.999, 730)
    print(f"   erase {erase:9.8f}: {'CURED' if ok else 'rebound'}")

print("\n== earliest stage (d14) | lockdown 365d ==")
for erase in [0.9, 0.99, 0.999, 0.9999, 0.99999]:
    ok = cured30("early", erase, 0.999, 365)
    print(f"   erase {erase:9.8f}: {'CURED' if ok else 'rebound'}")

print("\n== chronic | weaker 90% lockdown, 730d ==")
for erase in [0.999, 0.9999, 0.99999, 0.999999]:
    ok = cured30("chronic", erase, 0.90, 730)
    print(f"   erase {erase:9.8f}: {'CURED' if ok else 'rebound'}")

# ---------------- 2) minimal lockdown length at threshold erase --------------
print("\n== chronic | erase 0.9999999 | how short can lockdown be? ==")
for ln in [30, 60, 90, 180, 365]:
    ok = cured30("chronic", 0.9999999, 0.999, ln)
    print(f"   lockdown {ln:4d}d: {'CURED' if ok else 'rebound'}")

# ---------------- 3) what does a near-wipe trajectory look like? -------------
p = strike("chronic", 0.99999999, 0.999, 365)
t, YY, _ = simulate(p, t_end=30*Y, n_out=1200)
X, y5, y4, L, v5, v4, V, D = YY.T[0], YY.T[1], YY.T[2], YY.T[3], YY.T[4], YY.T[5], YY.T[4]+YY.T[5], YY.T[6]
for yr in [1, 2, 5, 10, 20, 30]:
    i = np.argmin(np.abs(t - yr*Y))
    print(f"   yr {yr:2d}: CD4={X[i]:6.1f}  L={L[i]*1000:10.4f}/mL  "
          f"VL={V[i]:10.3e}  D={D[i]:8.1f}")