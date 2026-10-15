"""stoch_demo.py — discrete-cell verdict for PAN vs the old G2 protocol.

Takes every perturbed world from seeds 777 & 4242, runs each protocol through
the end of the finite course, converts the reservoir at ART-stop into integer
latent cell counts (R5 + X4), then runs stoch_tail.Pcure (Monte-Carlo of the
last cells) to get P(cure) instead of the ODE's "cured / flip / undecidable".

Headline hypothesis: worlds where the ODE saw a 1e-90/mL proviral "ember" are
ARTIFACTS — less than one cell is zero cells. Discrete physics makes those
worlds cures.
"""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from stoch_tail import stoch_Pcure, R0_strain, latent_cells

Y = 365.25
ART_NEVER = 1e9
PERTURB = ["lam","d","beta","fL","a0","d_l","rho","delta_y","k","u","rho4","mu","s4","cD","rD","D_thr"]

def sane(p):
    pp = dict(p); pp.update(t_art_on=ART_NEVER, t_art_off=ART_NEVER,
                            shocks=[], strike_at=1e9, strike_refr_len=0)
    try:
        t, YY, _ = simulate(pp, t_end=2.5*Y, n_out=250, timeout=4)
    except BaseException:
        return False
    X = YY[:,0]; V = YY[:,4] + YY[:,5]
    if not (np.isfinite(X).all() and np.isfinite(V).all()):
        return False
    i = int(len(t)*0.98)
    return 120 < X[i] < 900 and 1e3 < V[i] < 1e7

def make_worlds(seed, count=60):
    rng = np.random.default_rng(seed); out = []
    while len(out) < count:
        p = dict(P)
        for k in PERTURB:
            p[k] *= rng.uniform(0.85, 1.15)
        if sane(p):
            out.append(p)
    return out

def proto(name):
    if name == "PAN":
        return dict(t_art_on=365, t_art_off=1460, shocks=[],
            strike_at=[365], strike_erase=[0.999999999], strike_kill=[1.0],
            strike_neutral=[0.999999], strike_refr=[0.9999], strike_refr_len=[1095],
            strike_healD=[0.9], strike_purge=[1.5], strike_cleanup=[1000],
            ccr5=0.999, ccr5_start=365)
    if name == "G2":
        return dict(t_art_on=365, t_art_off=1460,
                    shocks=list(range(430, 1460, 30)), lra=0.20)
    if name == "PAN_lite":  # weak variant: shallow erase + short purge, no CCR5
        return dict(t_art_on=365, t_art_off=1460, shocks=[],
            strike_at=[365], strike_erase=[0.999], strike_kill=[1.0],
            strike_neutral=[0.999], strike_refr=[0.99], strike_refr_len=[90],
            strike_healD=[0.5], strike_purge=[0.3], strike_cleanup=[10],
            ccr5=0.0, ccr5_start=1e9)
    raise KeyError(name)

def end_of_course(w, cf, time=0.0):
    """ODE state right at ART-stop (+time). Returns L/ml, x, and free virus at stop."""
    p = dict(w); p.update(cf)
    STOP = p["t_art_off"] + time
    t, YY, _ = simulate(p, t_end=STOP, n_out=240, timeout=240)
    i = np.argmin(np.abs(t - p["t_art_off"]))
    # l is cells/uL -> x1000 = cells/mL; free virus already in copies/mL
    return YY[i, 3] * 1000.0, YY[i, 0], YY[i, 4], YY[i, 5]

def main():
    V = 1e6          # body volume scale: cells per (reservoir cell/mL)  [headline]
    R0_MC = 600      # Monte-Carlo trajectories per world
    t0 = time.time()
    worlds = make_worlds(777) + make_worlds(4242)
    print(f"worlds: {len(worlds)}  (seeds 777 + 4242)  [continuous-time tail MC, R={R0_MC}]", flush=True)

    for name in ["PAN", "PAN_lite", "G2"]:
        cf = proto(name)
        cured = []; Lml = []; vstops = []
        for wi, w in enumerate(worlds):
            L_ml, x, v5s, v4s = end_of_course(w, cf)
            Lml.append(max(L_ml, 0.0)); vstops.append(max(v5s + v4s, 0.0))
            n = latent_cells(L_ml, volume=V)
            s4 = w["s4"]
            n4 = int(round(n * s4)); n5 = n - n4
            p = stoch_Pcure(n5, n4, w, a_react=1.0, purge_extra=0.0, purge_until=0,
                            art_until=w["t_art_off"], t0=w["t_art_off"],   # clock starts AT ART-stop
                            ccr5_eff=cf.get("ccr5", 0.0), R=R0_MC,
                            rng=np.random.default_rng(1000 + wi))
            cured.append(p)
            if wi % 20 == 0:
                print(f"  {name} world {wi}: L_stop={L_ml:.2e}/mL -> {n} cells "
                      f"({n4} X4)  P(cure)={p:.3f}", flush=True)
        arr = np.array(cured); Lm = np.array(Lml); vs = np.array(vstops)
        nz = arr == 1.0
        print(f"\n=== {name}: end-of-course reservoir across {len(worlds)} worlds ===")
        print(f"  reservoir at ART-stop : median {np.median(Lm):.1e}/mL, max {np.max(Lm):.1e}/mL")
        print(f"  free virus at stop    : median {np.median(vs):.1e} copies/mL "
              f"(expected seeding clone = beta*x/u * v ~ {np.median(vs)*1.5e-6:.1e} cells)")
        print(f"  worlds where leftover reservoir < 1 cell : {int((Lm*V < 1).sum())}/{len(Lm)}")
        print(f"  discrete P(cure): min={arr.min():.3f}  mean={arr.mean():.3f}  median={np.median(arr):.3f}")
        print(f"  guaranteed cures (0 cells at stop)        : {int(nz.sum())}/{len(arr)}")
        print(f"  expected late relapses (discrete)         : {sum(1-p for p in arr):.2f} of {len(arr)}",
              flush=True)

    # ------- immune recovery: does PAN heal the whole body, not just zero virus?
    cf = proto("PAN")
    cd4_stop, cd4_60, D60, V60 = [], [], [], []
    for wi, w in enumerate(worlds):
        p = dict(w); p.update(cf)
        try:
            t, YY, _ = simulate(p, t_end=60*Y, n_out=500, timeout=120, atol=1e-8, rtol=1e-6)
        except BaseException:
            try:
                t, YY, _ = simulate(p, t_end=60*Y, n_out=400, timeout=300, atol=1e-6, rtol=1e-4)
            except BaseException:
                t, YY = None, None
        if YY is None:
            cd4_stop.append(np.nan); cd4_60.append(np.nan); D60.append(np.nan); V60.append(np.nan)
            print(f"  PAN immune world {wi}: SKIPPED (timeout)")
            continue
        i = np.argmin(np.abs(t - 1460))
        cd4_stop.append(YY[i, 0])
        cd4_60.append(YY[-1, 0])
        D60.append(YY[-1, 6])
        V60.append(abs(YY[-1, 4]) + abs(YY[-1, 5]))
        if wi % 25 == 0:
            print(f"  PAN immune world {wi}: CD4 stop={YY[i,0]:.0f} -> 60y={YY[-1,0]:.0f}  D60={YY[-1,6]:.0f}", flush=True)
    cd4s = np.array(cd4_stop); cd4e = np.array(cd4_60); De = np.array(D60); Ve = np.array(V60)
    print("\n=== PAN: immune reconstitution (cure as a whole, not just virus) ===")
    print(f"  CD4 at ART-stop : median {np.median(cd4s):.0f} cells/uL (AIDS-range ~{np.median([w['lam']/w['d'] for w in worlds[:1]]):.0f} healthy)")
    print(f"  CD4 at 60y      : median {np.median(cd4e):.0f}  min {cd4e.min():.0f}"
          f"  | worlds recovered to >=500 : {(cd4e>=500).sum()}/120")
    print(f"  immune damage D : median {np.median(De):.0f} (healthy ~0; AIDS threshold {30000})"
          f"  max {De.max():.0f}")
    print(f"  virus at 60y    : median {np.median(Ve):.1e} copies/mL, max {Ve.max():.1e}")

    # ------- risk table: P(cure) vs surviving cells for a median world -------
    w0 = worlds[0]
    R0s = R0_strain(w0, 1.0, 0.994)    # effective R5 with CCR5 at 0.994
    R0x = R0_strain(w0, w0["rho4"], 0.0)
    print(f"\n(median world: R0_R5(eff with ccr5 .994)={R0s:.2f}, R0_X4={R0x:.2f})")
    print("residual at ART-stop | P(cure) no-CCR5 | P(cure) CCR5-pen (0.994)")
    for n4, n5 in [(1,0),(1,1),(2,0),(5,5),(20,2),(50,50),(200,0),(1000,0)]:
        pc1 = stoch_Pcure(n5, n4, w0, art_until=1460, t0=1460, ccr5_eff=0.0, R=2000,
                          rng=np.random.default_rng(5))
        pc2 = stoch_Pcure(n5, n4, w0, art_until=1460, t0=1460, ccr5_eff=0.994, R=2000,
                          rng=np.random.default_rng(6))
        print(f"   n4={n4:5d} n5={n5:5d}      |    {pc1:.3f}        |      {pc2:.3f}")
    print(f"\n[done in {time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()