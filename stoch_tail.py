"""stoch_tail.py — discrete-cell latent-reservoir engine.

When the reservoir falls below ~100 cells the ODE becomes invalid (you cannot
have a fraction of a cell, and single provirus events decide everything).
This engine takes the end-of-course reservoir from the ODE, converts it to an
integer number of latent cells (R5 and X4), and runs a continuous-time/master
Monte-Carlo of the last cells:

  per latent cell, per day:
    reactivate   @ a(t)      -> cell gone; if the strain can spread at this
                                moment, a takeoff roll happens with prob
                                (1 - 1/R0) — on a hit the patient relapses
    die          @ d_l       -> cell gone
    divide       @ rho       -> cell doubles (homeostatic proliferation)

  extinction = all latent cells gone  =>  CURE
  first successful takeoff            =>  RELAPSE

R5 spread is permanently blocked when a CCR5 knockout is active (this is why
PAN's permanent CCR5 knockout is the key move). X4 is never blockable.
"""
import numpy as np
from math import exp

def R0_strain(w, rho4=1.0, ccr5_eff=0.0):
    """Replication number of the reactivation loop at full CD4, per world."""
    x = w["lam"] / w["d"]
    b = w["beta"] * rho4 * (1 - ccr5_eff)
    return (1 - w["fL"]) * b * x * w["k"] / (w["u"] * w["delta_y"])

def latent_cells(L_per_ml, volume=1e6):
    """Map ODE-ish reservoir (cells/mL) to integer latent cell count."""
    n = int(round(max(L_per_ml, 0.0) * volume))
    return n

def stoch_Pcure(n5, n4, w, a_react=1.0, purge_extra=0.0, purge_until=0.0,
                art_until=1e12, t0=0.0, ccr5_eff=0.999, R=2000,
                t_horizon=60*365.25, rng=None):
    """Monte-Carlo P(cure) from (n5, n4) latent cells.

    a_react      = multiplier on baseline reactivation a0
    purge_extra  = activation added/day while t < purge_until (the tail purge)
    art_until    = day after which virus can spread again (ART off)
    ccr5_eff     = permanent CCR5 knockout efficiency (blocks R5 spread)
    """
    if rng is None:
        rng = np.random.default_rng(7)
    a0 = w["a0"] * a_react
    d_l = w["d_l"]
    rho = w["rho"]
    R05 = R0_strain(w, 1.0, ccr5_eff)
    R04 = R0_strain(w, w["rho4"], 0.0)

    n5v = np.full(R, n5, dtype=np.float64)
    n4v = np.full(R, n4, dtype=np.float64)
    dead = np.zeros(R, dtype=bool)     # True once outcome decided (cure/relapse)
    cured = np.zeros(R, dtype=bool)
    t = float(t0)          # clock starts where the measured reservoir lives
    dt = 1.0  # day; fine when per-day probs are small
    while t < t_horizon and not dead.all():
        a = a0 + (purge_extra if t < purge_until else 0.0)
        # per-day event probabilities (Poisson thinning)
        p_r   = 1 - exp(-a)
        p_d   = 1 - exp(-d_l)
        p_div = 1 - exp(-rho)
        # cells of settled trajectories are zero -> binomial gives zero events
        n5i = np.rint(np.maximum(n5v, 0)).astype(np.int64)
        n4i = np.rint(np.maximum(n4v, 0)).astype(np.int64)
        r5  = rng.binomial(n5i, p_r).astype(np.float64)
        r4  = rng.binomial(n4i, p_r).astype(np.float64)
        n5v = n5v - r5
        n4v = n4v - r4
        d5  = rng.binomial(n5i - r5.astype(np.int64), p_d).astype(np.float64)
        d4  = rng.binomial(n4i - r4.astype(np.int64), p_d).astype(np.float64)
        n5v = n5v - d5
        n4v = n4v - d4
        n5r = np.rint(np.maximum(n5v, 0)).astype(np.int64)
        n4r = np.rint(np.maximum(n4v, 0)).astype(np.int64)
        v5  = rng.binomial(n5r, p_div).astype(np.float64)
        v4  = rng.binomial(n4r, p_div).astype(np.float64)
        n5v = n5v + v5
        n4v = n4v + v4
        n5v = np.clip(n5v, 0, None); n4v = np.clip(n4v, 0, None)
        # takeoff rolls for reactivated cells (per reactivation event)
        if t >= art_until:
            live = ~dead
            if r4.any():
                p4 = 1 - 1/max(R04, 1.0)
                hit = rng.random(R) < (1 - (1 - p4) ** r4)
                dead[hit & live] = True
            if ccr5_eff < 1 and R05 > 1 and r5.any():
                p5 = 1 - 1/max(R05, 1.0)
                hit = rng.random(R) < (1 - (1 - p5) ** r5)
                dead[hit & live] = True
        # extinction verdict
        live = ~dead
        settled = live & (n5v + n4v < 1.0)
        cured[settled] = True
        dead[settled] = True
        t += dt
    # unresolved trajectories at horizon: only zero-cell ones count as cured
    live_zero = (~dead) & (n5v + n4v < 1.0)
    cured[live_zero] = True
    return float(cured.sum()) / R

def risk_table(w, cells_sweep, ccr5_eff=0.999, R=1200, **kw):
    """P(cure) over (x4, x5) cell-count grid → table for the report."""
    out = {}
    for n4, n5 in cells_sweep:
        p = stoch_Pcure(n5, n4, w, ccr5_eff=ccr5_eff, R=R, rng=np.random.default_rng(1234), **kw)
        out[(n4, n5)] = p
    return out