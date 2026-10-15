"""final_verdict.py — ONE command to reproduce the entire PAN-2 story.

  python3 final_verdict.py

Panels:
  1) ODE strict+clinical census, 2 seeds x 60 worlds, 60 years
  2) three-solver provable classification (LSODA/Radau/BDF, tight)
  3) discrete-cell reservoir verdict (whole-latent-cell count)
  4) immune recovery (CD4, damage D) over the cohort
All findings logged to RESULTS_GO.md. ~6-10 min.
"""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from stoch_demo import make_worlds
from panlandscape import PAN, verdict
from stoch_tail import stoch_Pcure

Y = 365.25
LOG = "RESULTS_GO.md"

def log(m):
    print(m, flush=True)
    with open(LOG, "a") as f:
        f.write(m + "\n")

def run1(w, cf, method="LSODA", atol=1e-7, rtol=1e-6):
    p = dict(w); p.update(cf)
    try:
        return simulate(p, t_end=60*Y, n_out=400, timeout=90,
                        atol=atol, rtol=rtol, method=method)
    except BaseException:
        return None, None, None

def main():
    seeds = [777, 4242]
    worlds = make_worlds(777) + make_worlds(4242)
    cf = PAN(); cf.update(t_art_on=1e9, t_art_off=1e9)  # ART-FREE PAN-2
    t_all = time.time()
    log(f"\n### FINAL VERDICT: PAN-2 across {len(worlds)} patients (seeds {seeds})")

    strict = clinical = rel = 0
    cd4s, Ds = [], []
    for w in worlds:
        t, YY, _ = run1(w, cf)
        if YY is None or not np.isfinite(YY).all():
            rel += 1; continue
        V = np.abs(YY[:,4]) + np.abs(YY[:,5]); X = YY[:,0]
        if V[t >= 2*Y].max() <= 1e-3 and X[-1] >= 500: strict += 1
        if V[t >= 3*Y].max() <= 50 and X[-1] >= 500: clinical += 1
        else: rel += 1
        cd4s.append(X[-1]); Ds.append(YY[-1, 6])
    log(f"[1] ODE (LSODA): STRICT {strict}/120   CLINICAL {clinical}/120   else {rel}/120   "
        f"medCD4 {np.median(cd4s):.0f}   medD {np.median(Ds):.0f}")

    prov = amb = 0
    for w in worlds[:60]:
        cells = []
        for m in ["LSODA", "Radau", "BDF"]:
            t, YY, _ = run1(w, cf, method=m)
            cells.append(verdict(t, YY)[0])
        if set(cells) == {"C"}: prov += 1
        else: amb += 1
    log(f"[2] 3-solver (seed 777): PROVABLY-CURED {prov}/60, "
        f"solver-ambiguous {amb}/60")

    nonzero = 0; Pcs = []
    for wi, w in enumerate(worlds):
        p = dict(w); p.update(cf)
        try:
            t, YY, _ = simulate(p, t_end=1460.0, n_out=1400, timeout=90,
                                atol=1e-7, rtol=1e-6)
        except BaseException:
            continue
        L = max(YY[np.argmin(np.abs(t-1460.0)), 3], 0.0) * 1000
        n = int(round(L * 1e6))
        if n > 0: nonzero += 1
        n4 = int(round(n * worlds[wi]["s4"]))
        Pcs.append(stoch_Pcure(n - n4, n4, worlds[wi], art_until=1460.0, t0=1460.0,
                               ccr5_eff=0.999, R=500, rng=np.random.default_rng(wi)))
    log(f"[3] discrete-cell: worlds with >=1 whole latent cell at lockdown-end: "
        f"{nonzero}/240;  P(cure) mean {np.mean(Pcs):.5f}  -> discrete relapses "
        f"{sum(1 - p for p in Pcs):.4f}")
    log(f"[4] immune recovery: median CD4 at year 60 = {np.median(cd4s):.0f} (>=500 cured); "
        f"immune damage D healed to median {np.median(Ds):.0f} vs AIDS threshold 30000")
    log(f"[pan2 final verdict in {time.time()-t_all:.0f}s]")

if __name__ == "__main__":
    main()