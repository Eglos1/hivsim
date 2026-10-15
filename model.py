"""
HIV-1 in-host dynamics model  (stable Nowak-May form + latent reservoir
                              + R5/X4 two-strain coreceptor escape)
========================================================================
States
  x   : uninfected CD4+ T cells        [cells/uL]
  y5  : productively infected, R5 (CCR5-tropic)  [cells/uL]
  y4  : productively infected, X4 (CXCR4-tropic) [cells/uL]
  l   : latently infected cells (the reservoir)  [cells/uL]
  v5  : R5 free virus                  [copies/mL]
  v4  : X4 free virus                  [copies/mL]
  D   : cumulative immune damage (drives AIDS progression, heals on
        deep suppression)

Why two strains? It is THE reason single-target CCR5 knock-out therapies
keep failing in the clinic: block CCR5 and the virus switches to CXCR4.
ART/PI, bNAbs and CAR-T act on both strains; CCR5 protection does not.

Interventions
  ART   : blocks new infection (RTI) + virion production (PI)
  SHOCK : latency-reversing-agent pulses raise activation a (purge)
  CCR5  : gene-therapy protection, blocks the R5 strain only
  CART  : armed T-cells, kill infected cells of either strain
  BNAB  : neutralizing-antibody window (both strains)

Units: x,y,l in cells/uL; v in copies/mL.
This is a mathematical model for study, not medical advice.
"""

import numpy as np
from scipy.integrate import solve_ivp

# ------------------------------------------------------------------ params
P = dict(
    lam = 10.0,          # CD4 cell source             (cells/uL/day)
    d   = 0.010,         # CD4 natural death           (1/day)
    beta= 1.0e-7,        # R5 infectivity              (mL/(copy*day))
    fL  = 1e-4,          # fraction of infections seeding latent pool
    a0  = 5e-4,          # baseline latency activation (1/day)
    d_l = 0.0030,        # latent cell death           (1/day)
    rho = 0.0011,        # latent cell growth          (1/day)
    delta_y = 1.0,       # infected-cell death (immune) (1/day)
    k   = 6.5e5,         # virions/mL per infected cell/day
    u   = 23.0,          # virus clearance             (1/day)

    # X4 strain / coreceptor escape
    rho4 = 0.95,         # X4 transmissibility vs R5
    mu   = 1e-6,         # R5->X4 mutation per R5 infection
    s4   = 1e-3,         # fraction of reservoir that is X4 provinces

    # damage -> AIDS progression
    D_thr = 30_000.0,
    D_cap = 300_000.0,
    rD    = 1.0e-5,
    kD    = 4.0e-6,
    cD    = 1.5e-4,      # immune attrition (CTL exhaustion)
    repair= 1.0e-4,

    # ART
    e_art_infect = 0.99999,
    e_art_prod   = 0.99999,

    # shock-and-kill
    lra       = 0.12,    # activation added per shock   (1/day)
    shock_d   = 5.0,     # shock duration               (days)
    shocks    = [],

    # immunotherapies (with start days so they can't cheat)
    ccr5       = 0.0,
    ccr5_start = 1e9,
    cart_boost = 0.0,
    cart_start = 1e9,
    t_cart_end = 1e9,   # day CAR-T cells are removed/expire (so finite,
                          # non-lifelong courses can be modeled)
    bnab       = (None, 0.0, 0.0),   # (start, neutralization frac, kill boost)
    bnab_w     = 120.0,              # window duration (days) -- long for
                                      # gene-therapy 'antibody factories'

    # "block-and-lock": permanent silencing of proviral REACTIVATION, per
    # strain, from silence_start on (0 = no effect). Counterfactual analog of
    # Tat-inhibitor 'block and lock' (e.g. didehydro-Cortistatin A): locked
    # proviruses can never wake up, even if they remain in the reservoir.
    silence5     = 0.0,   # fraction of R5-provirus reactivation silenced
    silence4     = 0.0,   # fraction of X4-provirus reactivation silenced
    silence_start= 1e9,

    # treatment schedule
    t_art_on  = 30.0,
    t_art_off = 760.0,

    # ---- "ONE-TIME STRIKE": a single, instantaneous intervention ----
    # A hypothetical mechanism applied at ONE time point (no lifelong drugs,
    # no permanent genetic change). Delivered in one session it:
    #   * erases a fraction of ALL latent proviruses in place  (strike_erase)
    #   * kills a fraction of every currently-infected cell    (strike_kill)
    #   * clears a fraction of free virus                      (strike_neutral)
    #   * briefly makes target cells refractory to infection   (strike_refr
    #     for strike_refr_len days afterwards) so leftover virus cannot reseed
    #   * resets accumulated immune damage & boosts CD4        (strike_healD,
    #     strike_boostx)
    strike_at      = 1e9,     # day of the strike (1e9 = never)
    strike_erase   = 0.0,     # fraction of latent reservoir eliminated in place
    strike_kill    = 0.0,     # fraction of productively-infected cells killed
    strike_neutral = 0.0,     # fraction of free virus cleared
    strike_refr    = 0.0,     # fractional reduction of infectivity after strike
    strike_refr_len= 0.0,     # how many days the target lockdown lasts
    strike_healD   = 0.0,     # fraction of accumulated damage reset
    strike_boostx  = 1.0,     # CD4 boost factor (capped)
    strike_purge   = 0.0,     # extra latency-activation rate DURING the lockdown
                              #   window (1/day) -> burn reservoir out fast
    strike_cleanup = 1.0,     # kill-rate multiplier for infected cells DURING
                              #   the lockdown window (clean up reactivators)
)

UNDET = 50.0   # copies/mL detection limit

# ---------------------------------------------------------------- dynamics
def _strike_list(par, wall=1e9):
    """Normalize the strike parameters into a list of strike descriptors.
    Each descriptor holds its own timing and strengths, enabling coordinated
    multi-strike ('double-tap') protocols. Scalars -> single strike; lists ->
    one strike per entry."""
    at = par["strike_at"]
    if np.isscalar(at):
        if at >= wall:
            return []
        n = 1
    else:
        n = len(at)
    def get(k):
        v = par[k]
        return [v]*n if np.isscalar(v) else list(v)
    return [dict(at=float(at if np.isscalar(at) else at[i]),
                 erase=float(get("strike_erase")[i]),
                 kill=float(get("strike_kill")[i]),
                 neutral=float(get("strike_neutral")[i]),
                 refr=float(get("strike_refr")[i]),
                 refr_len=float(get("strike_refr_len")[i]),
                 healD=float(get("strike_healD")[i]),
                 boostx=float(get("strike_boostx")[i]),
                 purge=float(get("strike_purge")[i]),
                 cleanup=float(get("strike_cleanup")[i]))
            for i in range(n)]

def make_deriv(par):
    beta0 = par["beta"]
    a0    = par["a0"]
    shocks = list(par["shocks"])
    shD    = par["shock_d"]
    lra    = par["lra"]
    bnab_start, bnab_n, bnab_k = par["bnab"]
    bnab_w = par["bnab_w"]
    strikes = _strike_list(par)

    def deriv(t, y):
        x, y5, y4, l, v5, v4, D = y
        art_on = par["t_art_on"] <= t < par["t_art_off"]
        v = v5 + v4

        # strain infectivities
        b5 = beta0
        if t >= par["ccr5_start"]:
            b5 *= (1 - par["ccr5"])          # CCR5 therapy only blocks R5
        b4 = beta0 * par["rho4"]             # X4 ignores CCR5 protection
        if art_on:
            b5 *= (1 - par["e_art_infect"])
            b4 *= (1 - par["e_art_infect"])
        if bnab_start is not None and bnab_start <= t <= bnab_start + bnab_w:
            b5 *= (1 - bnab_n)
            b4 *= (1 - bnab_n)
        # strike lockdown window (whichever strike is active right now)
        in_strike = None
        for s in strikes:
            if s["at"] <= t < s["at"] + s["refr_len"]:
                in_strike = s
                break
        if in_strike is not None:
            b5 *= (1 - in_strike["refr"])
            b4 *= (1 - in_strike["refr"])

        # infected-cell killing (immune attrition + CAR-T on both strains)
        dy = par["delta_y"]/(1 + par["cD"]*max(D,0))
        if par["cart_start"] <= t < par["t_cart_end"]:
            dy *= (1 + par["cart_boost"])
        if in_strike is not None and in_strike["cleanup"] != 1.0:
            dy *= in_strike["cleanup"]
        if bnab_start is not None and bnab_start <= t <= bnab_start + bnab_w:
            dy *= (1 + bnab_k)

        # virion production (PI block on ART)
        k_eff = par["k"] * (1 - (par["e_art_prod"] if art_on else 0.0))
        a = a0 + sum(lra for s in shocks if s <= t <= s + shD)
        if in_strike is not None:
            a += in_strike["purge"]

        # infection events
        e5 = b5*x*v5
        e4 = b4*x*v4

        # block-and-lock: silenced proviruses cannot reactivate
        sil5 = (1 - par["silence5"]) if t >= par["silence_start"] else 1.0
        sil4 = (1 - par["silence4"]) if t >= par["silence_start"] else 1.0

        dx  = par["lam"] - par["d"]*x - (e5 + e4) - par["kD"]*D*x
        dy5 = (1 - par["fL"])*e5*(1 - par["mu"]) + a*l*(1 - par["s4"])*sil5 - dy*y5
        dy4 = (1 - par["fL"])*(e4 + e5*par["mu"]) + a*l*par["s4"]*sil4    - dy*y4
        dl  = par["fL"]*(e5 + e4) + par["rho"]*l - (a + par["d_l"])*l
        dv5 = k_eff*y5 - par["u"]*v5
        dv4 = k_eff*y4 - par["u"]*v4

        if v > par["D_thr"]:
            dD = par["rD"]*(min(v, par["D_cap"]) - par["D_thr"])
        elif v < UNDET:
            dD = -par["repair"]*max(D,0.0)
        else:
            dD = 0.0

        return (dx, dy5, dy4, dl, dv5, dv4, dD)

    return deriv

# ---------------------------------------------------------------- simulate
def _events(par, t0, t_end):
    """Every time at which the vector field discontinuously changes."""
    ev = {t0, t_end, par["t_art_on"], par["t_art_off"],
          par["ccr5_start"], par["cart_start"], par["t_cart_end"],
          par["silence_start"]}
    for s in _strike_list(par):
        ev.add(s["at"])
        if s["refr_len"] > 0:
            ev.add(s["at"] + s["refr_len"])
    bs, *_ = par["bnab"]
    if bs is not None:
        ev.add(bs); ev.add(bs + par["bnab_w"])
    for s in par["shocks"]:
        ev.add(s); ev.add(s + par["shock_d"])
    ev = sorted(float(e) for e in ev if t0 <= e <= t_end)
    return ev

_ALARM_HANDLER = None

def _arm(timeout):
    """SIGALRM-based watchdog on Linux main thread. Returns disarm func."""
    global _ALARM_HANDLER
    import signal
    class Timeout(BaseException): pass
    def kill(sig, frm): raise Timeout("simulation watchdog")
    try:
        _ALARM_HANDLER = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, kill)
        signal.setitimer(signal.ITIMER_REAL, timeout)
    except (ValueError, TypeError):   # not main thread / bad input
        return lambda: None
    return lambda: signal.setitimer(signal.ITIMER_REAL, 0)

def _disarm():
    import signal
    signal.setitimer(signal.ITIMER_REAL, 0)
    if _ALARM_HANDLER is not None:
        signal.signal(signal.SIGALRM, _ALARM_HANDLER)

def simulate(par, t0=0.0, t_end=365.25*12, y0=None, n_out=None, timeout=None,
             rtol=1e-4, atol=1e-6, max_step=2.0, min_step=1e-6, method="LSODA"):
    if y0 is None:
        y0 = np.array((par["lam"]/par["d"], 0.0, 0.0, 0.0, 50.0, 0.0, 0.0))
    if n_out is None:
        n_out = int(t_end - t0) + 1
    teval = np.linspace(t0, t_end, n_out)
    disarm = _arm(timeout) if timeout else None

    times, states = [t0], [np.asarray(y0, float)]
    Yside = np.asarray(y0, float)
    try:
        for a, b in zip(_events(par, t0, t_end)[:-1], _events(par, t0, t_end)[1:]):
            if b <= a:
                continue
            # always sample the segment endpoints explicitly so sharp
            # discontinuities (strikes) are captured exactly once
            sub = np.unique(np.concatenate([
                teval[(teval > a) & (teval < b)], [a, b]]))
            res = solve_ivp(make_deriv(par), (a, b), Yside, method=method,
                            t_eval=sub, rtol=rtol, atol=atol,
                            min_step=min_step, max_step=max_step)
            if not res.success:
                raise RuntimeError(f"segment [{a:.1f},{b:.1f}]: {res.message}")
            Yside = res.y[:, -1]
            np.maximum(Yside, 0.0, out=Yside)   # populations cannot go negative
            # apply any one-time strikes whose time this segment ends at
            for s in _strike_list(par):
                if abs(b - s["at"]) < 1e-6:
                    Yside[3] *= (1.0 - s["erase"])
                    Yside[1] *= (1.0 - s["kill"])
                    Yside[2] *= (1.0 - s["kill"])
                    Yside[4] *= (1.0 - s["neutral"])
                    Yside[5] *= (1.0 - s["neutral"])
                    Yside[6] *= (1.0 - s["healD"])
                    if s["boostx"] != 1.0:
                        Yside[0] = min(Yside[0]*s["boostx"],
                                       1.5*par["lam"]/par["d"])
                    np.maximum(Yside, 0.0, out=Yside)
            times.extend(res.t[1:])
            states.extend(res.y.T[1:])
    except BaseException:
        if disarm: disarm()
        raise
    if disarm:
        disarm()

    t = np.asarray(times)
    Y = np.asarray(states)
    # dedupe timestamps (segments overlap at boundary event times; keep the
    # LAST state at each time, i.e. the post-transform one)
    if t.shape[0] != np.unique(t).shape[0]:
        _, idx = np.unique(t, return_index=True)
        t, Y = t[idx], Y[idx]
    # map onto the requested teval grid with clamped, monotone interpolation
    if t.shape[0] != teval.shape[0] or not np.allclose(t, teval):
        Te = np.asarray(teval, float)
        Ye = np.empty((teval.shape[0], Y.shape[1]))
        for c in range(Y.shape[1]):
            Ye[:, c] = np.interp(Te, t, Y[:, c], left=Y[0, c], right=Y[-1, c])
        t, Y = Te, Ye
    return t, Y, y0

# ----------------------------------------------------------------- metrics
def unpack(t, Y):
    X, Y5, Y4, L, V5, V4, D = Y.T
    return X, Y5, Y4, L, V5, V4, V5+V4, D

def rebound_info(t, Y, par):
    X, Y5, Y4, L, V5, V4, V, D = unpack(t, Y)
    if par["t_art_off"] is None or par["t_art_off"] > t[-1]:
        return None
    m = t >= par["t_art_off"]
    if not m.any():
        return None
    over = np.where(V[m] >= UNDET)[0]
    return t[m][over[0]] - par["t_art_off"] if len(over) else None

def cured_check(t, Y, par, window_years=3.0):
    """Off-ART viremia stays <50 copies/mL for the window -> remission."""
    X, Y5, Y4, L, V5, V4, V, D = unpack(t, Y)
    stop = par["t_art_off"]
    tail = (t >= stop + window_years*365.25) & (t < t[-1])
    return bool(tail.any() and (V[tail] < UNDET).all() and X[tail][-1] >= 500)

def sustained_remission(t, Y, window_years=3.0):
    """Cure criterion that works whether or not ART was ever stopped:
       during the FINAL window_years of the sim the virus stays
       <50 copies/mL and CD4 is healthy (>=500) and the reservoir is gone."
    """
    X, Y5, Y4, L, V5, V4, V, D = unpack(t, Y)
    tail = t >= t[-1] - window_years*365.25
    if not tail.any() or V[tail].size == 0:
        return False
    return bool((V[tail] < UNDET).all()
                and X[tail][-1] >= 500
                and (L[tail] < 1e-3).all())

def summarize(t, Y, par, label):
    X, Y5, Y4, L, V5, V4, V, D = unpack(t, Y)
    print(f"--- {label} ---")
    for yr in [1, 2, 5, 10]:
        i = np.argmin(np.abs(t - yr*365.25))
        print(f"  yr {yr:>2}: CD4={X[i]:6.1f}  VL={V[i]:10.0f} copies/mL  "
              f"X4={100*V4[i]/max(V[i],1e-9):5.1f}%  "
              f"reservoir={L[i]*1000:8.2f}/mL  D={D[i]:6.1f}")
    if V[-1] >= UNDET:
        print(f"  end : CD4={X[-1]:6.1f}  VL={V[-1]:10.0f}  "
              f"X4={100*V4[-1]/max(V[-1],1e-9):5.1f}%  D={D[-1]:6.1f}")
    rb = rebound_info(t, Y, par)
    if rb is not None:
        print(f"  ART interrupted day {par['t_art_off']:.0f} -> rebound in "
              f"{rb:.0f} days")
    else:
        print(f"  ART interrupted day {par['t_art_off']:.0f} -> "
              f"NO rebound within simulation window")
        print("  ** off-ART remission" if cured_check(t, Y, par)
              else "  (partial: not yet confirmed for 3-year window)")