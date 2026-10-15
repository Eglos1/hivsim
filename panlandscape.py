"""panlandscape.py — autonomous protocol-hardening sweeps around PAN.

Subcommands (each logged quietly to RESULTS_GO.md):
  erase    strike-efficiency ladder  erase in {0.99 ... 0.999999999}
  start    therapy-start-day ladder  day 7 (acute) ... 1460 (late chronic)
  pw       purge x lockdown-window landscape (isolate 'window beats intensity')
  artfree  PAN without the ART tail
  purge    purge-strength ladder
Shared: 60 perturbed worlds (seed 777), tight-ish solver (atol 1e-7, rtol 1e-6).
"""
import sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from stoch_demo import make_worlds

Y = 365.25
LOG = "RESULTS_GO.md"

def log(txt):
    print(txt, flush=True)
    with open(LOG, "a") as f:
        f.write(txt + "\n")

def PAN(strike_at=365, erase=0.999999999, purge=1.5, win=1095,
        art_off=1460, ccr5=0.999):
    return dict(t_art_on=strike_at, t_art_off=art_off, shocks=[],
        strike_at=[strike_at], strike_erase=[erase], strike_kill=[1.0],
        strike_neutral=[0.999999], strike_refr=[0.9999], strike_refr_len=[win],
        strike_healD=[0.9], strike_purge=[purge], strike_cleanup=[1000],
        ccr5=ccr5, ccr5_start=strike_at)

def run(w, cf):
    p = dict(w); p.update(cf)
    for tol in [(1e-6, 1e-4, 350), (1e-4, 1e-3, 250)]:
        try:
            t, YY, _ = simulate(p, t_end=60*Y, n_out=tol[2], timeout=90,
                                atol=tol[0], rtol=tol[1])
            return t, YY
        except BaseException:
            continue
    return None, None

def verdict(t, YY, rem_start=2*Y):
    if YY is None or not np.isfinite(YY).all(): return "E", None
    V = np.abs(YY[:,4]) + np.abs(YY[:,5])
    X = YY[:,0]
    m = t >= rem_start
    ok = (V[m].max() <= 1e-3) and (X[-1] >= 500)
    return ("C" if ok else "R"), dict(cd4_final=X[-1], cd4_nadir=X[t < rem_start].min() if (t < rem_start).any() else X[-1], D_final=YY[-1,6])

def census(worlds, cf, rem_start=2*Y):
    stats = []
    for w in worlds:
        t, YY = run(w, cf)
        v, s = verdict(t, YY, rem_start)
        stats.append((v, s))
    C = sum(1 for v, _ in stats if v == "C")
    R = len(stats) - C - sum(1 for v, _ in stats if v == "E")
    E = sum(1 for v, _ in stats if v == "E")
    return dict(C=C, R=R, E=E, stats=stats)

def do_erase(worlds):
    steps = [0.99, 0.999, 0.9999, 0.99999, 0.999999,
             0.9999999, 0.99999999, 0.999999999]
    log("\n### E1 strike-efficiency ladder  (erase@365, purge1.5~win1095, ccr5 .999)")
    for e in steps:
        t0 = time.time()
        c = census(worlds, PAN(erase=e))
        log(f"erase={e:.9f}  C={c['C']}/60  R={c['R']}  E={c['E']}  ({time.time()-t0:.0f}s)")

def do_start(worlds):
    log("\n### E2 therapy-start-day ladder  (PAN core, running ART starts at strike)")
    for day in [7, 30, 150, 365, 730, 1095, 1460]:
        t0 = time.time()
        c = census(worlds, PAN(strike_at=day, art_off=day+1095),
                   rem_start=day + 365)
        nadirs = [s["cd4_nadir"] for v, s in c["stats"] if s]
        log(f"start day {day:5d}  C={c['C']}/60  R={c['R']}  E={c['E']}  "
            f"| pre-tx CD4 nadir median={np.median(nadirs):5.0f} "
            f"(sick while waiting)  ({time.time()-t0:.0f}s)")

def do_pw(worlds, subset=20):
    log("\n### E3 purge x lockdown-window landscape  (purge/win; die = window beats intensity)")
    wins = [180, 365, 730, 1095]
    purges = [0.5, 1.0, 1.5, 2.5]
    ws = worlds[:subset]
    log("purge\\win | " + " | ".join(f"w{int(w)}" for w in wins))
    for p in purges:
        row = [f"p{p:<5.2f}"]
        for w_ in wins:
            c = census(ws, PAN(purge=p, win=w_))
            row.append(f"C{c['C']}/{subset}")
        log("  " + " | ".join(row))

def do_artfree(worlds):
    log("\n### E4 PAN minus the ART tail  (strike+purge+CCR5, no ART at all)")
    t0 = time.time()
    na = dict(P); na.update(PAN())
    na.update(t_art_on=1e9, t_art_off=1e9)
    c = census(worlds, na)
    log(f"ART-free PAN   C={c['C']}/60  R={c['R']}  E={c['E']}  ({time.time()-t0:.0f}s)")

def do_purge(worlds):
    log("\n### E5 purge-strength ladder  (win=1095 fixed)")
    for p in [0.0, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0]:
        t0 = time.time()
        c = census(worlds, PAN(purge=p))
        log(f"purge={p:<4.1f} (win1095)  C={c['C']}/60  R={c['R']}  E={c['E']}  ({time.time()-t0:.0f}s)")

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    worlds = make_worlds(777)
    log(f"\n=====================================================================\n"
        f"panlandscape: {cmd}  ({len(worlds)} worlds seed 777, 60y, atol 1e-7)")
    t_all = time.time()
    if cmd in ("all", "erase"): do_erase(worlds)
    if cmd in ("all", "start"): do_start(worlds)
    if cmd in ("all", "pw"):    do_pw(worlds)
    if cmd in ("all", "artfree"): do_artfree(worlds)
    if cmd in ("all", "purge"): do_purge(worlds)
    log(f"[panlandscape:{cmd} done in {time.time()-t_all:.0f}s]")