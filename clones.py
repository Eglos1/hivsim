"""clones.py — clone-level stochastic reservoir dynamics.

The survival question isn't just "how many cells" but "how many distinct
clonal proviruses" — and critically whether a reactivation burst can RE-SEED
new latent clones (the "reservoir feeds itself" persistence loop).

Each surviving latent cell belongs to a clone (one integration event). Daily:
  reactivate (removed; if resettable viremia exists, with prob f_reseed the
              burst seeds a NEW clone),
  die, or divide (clone grows, stays same clone).
Strains R5/X4 tracked separately; CCR5 knockout blocks R5 reactivation spread.

Deliverables:
  S1 3 X4 clones x1 cell, no reseed          -> P(cure) ~ baseline die-roll
  S2 same, WITH reseed 4e-4 / reactivation   -> persistence loop appears
  S3 40 R5 clones x5 cells, no CCR5          -> R5 bombs, P(cure) low
  S4 40 R5 clones x5 cells, CCR5 locked      -> R5 neutered, HIGHLY cured
  + diversity (Shannon entropy over clones) trajectory for S1 vs S2.
"""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P
from stoch_tail import R0_strain
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

Y = 365.25

def clone_sim(sizes5, sizes4, w, art_until=1460.0, t0=1460.0,
              ccr5_eff=0.0, f_reseed=0.0, R=1500, horizon=60*Y, seed=9):
    """MC over R trajectories of clone-structured reservoirs.

    sizes5/sizes4: (list of clone sizes) per strain, per initial trajectory.
    Returns (P_cure, meanH_t) where H_t = Shannon entropy of clone sizes
    averaged over live trajectories at sampled years.
    """
    rng = np.random.default_rng(seed)
    a0, d_l, rho = w["a0"], w["d_l"], w["rho"]
    p_r, p_d, p_div = 1-np.exp(-a0), 1-np.exp(-d_l), 1-np.exp(-rho)
    R05 = R0_strain(w, 1.0, ccr5_eff)
    R04 = R0_strain(w, w["rho4"], 0.0)
    p5t = 1 - 1/max(R05, 1.0) if R05 > 1 else 0.0
    p4t = 1 - 1/max(R04, 1.0) if R04 > 1 else 0.0

    # per trajectory: row = (clone_id, size, strain) simulation via compact
    # representation: for each trajectory keep two lists of clone sizes
    traj = []
    for _ in range(R):
        c5 = list(sizes5); c4 = list(sizes4)
        traj.append((c5, c4))
    live = np.ones(R, dtype=bool)
    cured = np.zeros(R, dtype=bool)
    t = float(t0)
    stem = []; H = []
    n_clone_max = 0   # global clone id counter (for reseeded clones)
    while t < horizon and live.any():
        spread = t >= art_until
        for i in np.where(live)[0]:
            c5, c4 = traj[i]
            for lst, strain, pt in ((c5, "5", p5t), (c4, "4", p4t)):
                n = len(lst)
                if n == 0 or sum(lst) == 0:
                    lst[:] = []
                    continue
                total = sum(lst)
                # reactivations and remove from specific clones (multinomial)
                r = int(rng.binomial(total, p_r))
                if r:
                    alloc = rng.multinomial(r, np.array(lst)/total)
                    nlst = [max(lst[k]-alloc[k], 0) for k in range(n)]
                    lst[:] = [s for s in nlst if s > 0]
                    # takeoff if this strain can spread (relapse -> kill traj)
                    if spread and pt > 0:
                        if rng.random() < (1-(1-pt)**r):
                            live[i] = False
                            lst[:] = []; c4[:] = []
                            continue
                    # reseeding: fresh integration sites from the burst
                    if f_reseed > 0 and rng.random() < (1-(1-f_reseed)**r):
                        lst.append(1)                 # brand-new clone
            # deaths (multinomial across clones, no strain change)
            for lst in (c5, c4):
                if not lst: continue
                n = len(lst); total = sum(lst)
                if total == 0:
                    lst[:] = []
                    continue
                dd = int(rng.binomial(total, p_d))
                if dd:
                    alloc = rng.multinomial(dd, np.array(lst)/total)
                    nlst = [max(lst[k]-alloc[k], 0) for k in range(n)]
                    lst[:] = [s for s in nlst if s > 0]
            # divisions (within clone: binomial per clone)
            for lst in (c5, c4):
                nb = len(lst)
                if nb == 0: continue
                add = rng.binomial(np.array(lst, dtype=np.int64), p_div)
                for k in range(nb):
                    lst[k] += int(add[k])
        # extinction within live trajectories
        for i in np.where(live)[0]:
            c5, c4 = traj[i]
            if not c5 and not c4:
                cured[i] = True; live[i] = False
        # sample entropy every ~2 years
        if int(t) % int(2*Y) == 0 and live.any():
            Hs = []
            for i in np.where(live)[0]:
                sizes = traj[i][0] + traj[i][1]
                tot = sum(sizes)
                if tot > 0:
                    pks = np.array(sizes)/tot
                    Hs.append(-(pks*np.log(pks)).sum())
            stem.append(t/Y); H.append(np.mean(Hs) if Hs else 0.0)
        t += 1.0
    H_t = (np.array(stem), np.array(H))
    return cured.sum()/R, H_t

def main():
    w = dict(P)
    print("clone-level stochastic dynamics (median world)")
    print(f"  R0_R5(no ccr5)={R0_strain(w,1.0,0.0):.2f}  "
          f"R0_R5(ccr5 .999)={R0_strain(w,1.0,0.999):.2f}  "
          f"R0_X4={R0_strain(w,w['rho4'],0.0):.2f}")
    S1, h1 = clone_sim([], [1,1,1], w, ccr5_eff=0.999, f_reseed=0.0)
    S2, h2 = clone_sim([], [1,1,1], w, ccr5_eff=0.999, f_reseed=4e-4)
    S3, _ = clone_sim([5]*40, [], w, ccr5_eff=0.0,  f_reseed=0.0)
    S4, _ = clone_sim([5]*40, [], w, ccr5_eff=0.999, f_reseed=0.0)
    print("\nS1  3 X4 clones x1, NO reseed          : P(cure) = %.3f" % S1)
    print("S2  3 X4 clones x1, reseed f=4e-4       : P(cure) = %.3f  <- persistence loop" % S2)
    print("S3  40 R5 clones x5, no CCR5 lock       : P(cure) = %.3f" % S3)
    print("S4  40 R5 clones x5, CCR5 lock 0.999    : P(cure) = %.3f" % S4)

    # diversity figure: Shannon entropy of surviving clone pool over time
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.plot(h1[0], h1[1], color="#1f6fb2", marker="o", ms=3, lw=2, label="no reseed")
    ax1.plot(h2[0], h2[1], color="#c0392b", marker="s", ms=3, lw=2, label="reseed 4e-4")
    ax1.set_xlabel("years"); ax1.set_ylabel("Shannon entropy of clones (H)")
    ax1.set_title("Clone diversity among survivors\n(3 X4 clones start)")
    ax1.legend(); ax1.grid(alpha=.3)
    ax2.bar(["S1 no-reseed", "S2 reseed", "S3 R5 no-lock", "S4 R5 +lock"],
            [S1, S2, S3, S4], color=["#1f6fb2", "#c0392b", "#7d3c98", "#2e8b57"])
    ax2.set_ylim(0, 1.05); ax2.set_ylabel("P(cure) 60y")
    ax2.set_title("Clone-level cure probability")
    plt.tight_layout(); plt.savefig("clones.png", dpi=150)
    print("\nsaved clones.png")

if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"[done in {time.time()-t0:.0f}s]")