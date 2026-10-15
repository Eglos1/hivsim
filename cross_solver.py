"""cross_solver.py — definitive verdict: a world is 'provably cured' only if
it cures under THREE independent stiff ODE solvers (LSODA, Radau, BDF) at two
tolerance families. This resolves the knife-edge 'flips' once and for all by
refusing to call a coin-flip a cure. Logs to RESULTS_GO.md.
"""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from stoch_demo import make_worlds
from panlandscape import PAN, verdict

Y = 365.25
METHODS = ["LSODA", "Radau", "BDF"]
TOLS = [("tight", 1e-7, 1e-6), ("loose", 1e-5, 1e-4)]  # (name, atol, rtol)

def run(w, cf, method, atol, rtol):
    p = dict(w); p.update(cf)
    try:
        t, YY, _ = simulate(p, t_end=60*Y, n_out=400, timeout=100,
                            atol=atol, rtol=rtol, method=method)
        return t, YY
    except BaseException:
        return None, None

def classify(w, cf):
    """Returns 'C' if cured in >=5 of the 6 (method,tol) cells and NO 'R'
    observed anywhere; 'R' if any cell says relapse and majority says R;
    'A' otherwise (ambiguous knife-edge)."""
    cells = []
    for method in METHODS:
        for name, atol, rtol in TOLS:
            t, YY = run(w, cf, method, atol, rtol)
            cells.append(verdict(t, YY)[0])
    nC, nR, nE = cells.count("C"), cells.count("R"), cells.count("E")
    if nE:  # a solver refused to integrate: conservative -> ambiguous
        return "A"
    if nR == 0 and nC > 0:
        return "C"
    if nC == 0 and nR > 0:
        return "R"
    return "A"

def main():
    worlds = make_worlds(777)
    log = open("RESULTS_GO.md", "a")
    msg = (f"\n### E6 cross-integrator definitive verdict (3 solvers x 2 tol families, "
           f"{len(worlds)} worlds, seed 777)")
    print(msg); log.write(msg + "\n")
    for tag, mode in [("PAN (ART 3y + CCR5 .999)", "art"),
                      ("PAN ART-FREE + CCR5 .999", "artfree")]:
        cf = PAN()
        if mode == "artfree":
            cf.update(t_art_on=1e9, t_art_off=1e9)
        t0 = time.time()
        out = [classify(w, cf) for w in worlds]
        q = dict(C=out.count("C"), R=out.count("R"), A=out.count("A"))
        msg = (f"{tag}: PROVABLY-CURED={q['C']}/60  definitive-relapse={q['R']}  "
               f"ambiguous={q['A']}  ({time.time()-t0:.0f}s)")
        print(msg); log.write(msg + "\n")
    log.close()

if __name__ == "__main__":
    main()