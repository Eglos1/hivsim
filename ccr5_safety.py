"""ccr5_safety.py — what happens when you block CCR5 without killing the
reservoir? The clinical fear: permanent CCR5 knockout is a patient-assisted
selection experiment for X4 escape. We measure that fear across 60 perturbed
worlds (seed 777) under five regimens:

  bas    untreated chronic HIV
  c5     permanent CCR5 knockout only (0.99 @ day 365, no ART)
  art    3-yr ART only (365 -> 1460)
  artc5  ART + permanent CCR5 knockout (PAN WITHOUT the strike+purge)
  PAN    full champion (strike+purge runs through ART tail + CCR5 permanent)

Metrics per regimen: X4 escape fraction, median days to X4 escape, dominance
of X4 at escape, AIDS fraction by 60y, cure fraction.
"""
import warnings, time
import numpy as np
warnings.filterwarnings("ignore")
from model import P, simulate
from stoch_demo import make_worlds   # reuse exact same perturbed roster

Y = 365.25

def cond(name):
    return {
        "bas":   dict(),
        "c5":    dict(t_art_on=1e9, t_art_off=1e9, shocks=[],
                      strike_at=1e9, strike_refr_len=0,
                      ccr5=0.99, ccr5_start=365),
        "art":   dict(t_art_on=365, t_art_off=1460, shocks=[],
                      strike_at=1e9, strike_refr_len=0),
        "artc5": dict(t_art_on=365, t_art_off=1460, shocks=[],
                      strike_at=1e9, strike_refr_len=0,
                      ccr5=0.99, ccr5_start=365),
        "PAN":   dict(t_art_on=365, t_art_off=1460, shocks=[],
                      strike_at=[365], strike_erase=[0.999999999], strike_kill=[1.0],
                      strike_neutral=[0.999999], strike_refr=[0.9999], strike_refr_len=[1095],
                      strike_healD=[0.9], strike_purge=[1.5], strike_cleanup=[1000],
                      ccr5=0.999, ccr5_start=365),
    }[name]

def run60(w, cf):
    p = dict(w); p.update(cf)
    try:
        t, YY, _ = simulate(p, t_end=60*Y, n_out=1500, timeout=120, atol=1e-8, rtol=1e-6)
        return t, YY
    except BaseException:
        try:
            t, YY, _ = simulate(p, t_end=60*Y, n_out=1200, timeout=300, atol=1e-6, rtol=1e-4)
            return t, YY
        except BaseException:
            return None, None

def analyze(t, YY):
    if YY is None or not np.isfinite(YY).all():
        return dict(x4_esc=np.nan, esc_t=np.nan, x4_dom=np.nan,
                    aids=np.nan, cd4_60=np.nan, cure=np.nan)
    v4, v5, X = YY[:,5], YY[:,4], YY[:,0]
    V = v4 + v5
    # X4 escape: v4 above 1e-3 copies/mL at any point after day 400
    m = t >= 400
    esc = (v4[m] > 1e-3).any()
    esc_t = t[m][v4[m] > 1e-3].min() if esc else None
    dom = None
    if esc:
        i = np.where(v4[m] > 1e-3)[0][0]
        dom = v4[m][i] > v5[m][i]
    aids = (X[t >= 30*Y] < 200).any() or not np.isfinite(X[-1])
    cured = (np.abs(V[t >= 2*Y]).max() <= 1e-3) and (X[-1] >= 500)
    return dict(x4_esc=esc, esc_t=esc_t, x4_dom=dom,
                aids=bool(aids), cd4_60=X[-1], cure=bool(cured))

def main():
    worlds = make_worlds(777)
    print(f"safety study: {len(worlds)} perturbed worlds (seed 777) x 5 regimens")
    hdr = ("regimen  | X4-escape | median day | X4-dominates | AIDS 60y | cure  | median CD4 60y")
    print(hdr); print("-" * len(hdr))
    for nm in ["bas", "c5", "art", "artc5", "PAN"]:
        cf = cond(nm)
        esc_t, n_esc, n_dom, n_aids, n_cure, cd4 = [], 0, 0, 0, 0, []
        for wi, w in enumerate(worlds):
            t, YY = run60(w, cf)
            r = analyze(t, YY)
            if r["x4_esc"] is True:
                n_esc += 1; esc_t.append(r["esc_t"])
                if r["x4_dom"]: n_dom += 1
            if r["aids"]: n_aids += 1
            if r["cure"]: n_cure += 1
            cd4.append(r["cd4_60"])
        print(f"{nm:8s} | {n_esc:2d}/{len(worlds)}        | "
              f"{np.median(esc_t)/Y:6.1f}y | {n_dom:3d}/{n_esc}         | "
              f"{n_aids:2d}/{len(worlds)}     | {n_cure:2d}/{len(worlds)}  | {np.nanmedian(cd4):7.0f}", flush=True)

if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"[done in {time.time()-t0:.0f}s]")