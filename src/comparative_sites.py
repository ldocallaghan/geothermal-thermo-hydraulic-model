"""
comparative_sites.py
===================
Run the site-evaluation tool across four European candidate provinces and
tabulate. The point: the SAME models score very different verdicts, driven by
real geotherm + stress data -- which is the pitch-ready output.

Sites & sources (gradients/stress from cited literature; deep-stress ratios are
literature-typical estimates where not directly measured at superhot depth --
flagged as the key uncertainty):
  - Upper Rhine Graben / Soultz (FR): 200C@5km, Shmin~0.54Sv, SHmax~Sv. Genter;
    Cornet/Valley "Stress State at Soultz".
  - Larderello (IT): vapour-dominated, 350C@2.2km, supercritical K-horizon 3-7km
    (= brittle-ductile transition); post-collisional EXTENSION, normal faulting.
    Bertini/Gianelli; DESCRAMBLE Venelle-2.
  - United Downs / Carnmenellis (UK): 190C@5km, ~33-35C/km radiogenic granite;
    strike-slip stress (Cornubian), high horizontal anisotropy (Pine&Batchelor,
    Rosemanowes). Reinecker/Ledingham UDDGP.
  - Pannonian Basin (HU): ~45-50C/km, heat flow 90-100 mW/m2, Miocene back-arc
    EXTENSION, thin crust. Lenkey/Horvath; Toth geothermal atlas.
"""
import numpy as np
import geo_constants as C
import site_evaluation as se
from site_evaluation import (SiteProfile, SOULTZ, StressProfile, evaluate,
                             transitional_bounds)


def make_geotherm(segments):
    """segments: list of (z_top, T_top, gradient) cumulative; returns vectorized fn."""
    def gt(z):
        z = np.asarray(z, dtype=float)
        conds, vals = [], []
        for (z0, T0, g), (z1, *_ ) in zip(segments, segments[1:] + [(np.inf, 0, 0)]):
            conds.append((z >= z0) & (z < z1))
            vals.append(T0 + g * (z - z0))
        return np.select(conds, vals, default=vals[-1])
    return gt


LARDERELLO = SiteProfile(
    name="Larderello (Italy)",
    geotherm=make_geotherm([(0, 15, 0.152), (2200, 350, 0.077)]),
    target_depth=2850.0, target_T=400.0,            # just below the ~450C K-horizon
    Sv_grad=2600 * 9.81, K0_min=0.55, SHmax_over_Sv=0.75,   # normal-fault, low aniso
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.6, E_rock=45e9, UCS=140e6,
    # No measured magnitudes found. Regime: normal faulting alongside left-
    # lateral strike-slip, the intermediate stress switching between vertical
    # and horizontal (Brogi et al., field data and focal mechanisms in the Lago
    # Basin; full citation to confirm), i.e. SHmax ~ Sv. Hydrostatic Pp is the
    # default; the vapour-dominated reservoir is likely underpressured.
    stress_cases=transitional_bounds(
        2600 * 9.81, source="regime: normal/strike-slip transition (Brogi et al.)"),
    stress_basis="regime bounds")

CORNWALL = SiteProfile(
    name="United Downs / Carnmenellis (UK)",
    geotherm=make_geotherm([(0, 15, 0.035), (5000, 190, 0.028)]),
    target_depth=12500.0, target_T=400.0,
    Sv_grad=2630 * 9.81, K0_min=0.55, SHmax_over_Sv=1.40,   # strike-slip, HIGH aniso
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=3.3, E_rock=60e9, UCS=180e6,
    # Reinecker et al. (2021) section 6.3; see data/ud1/reinecker2021.md. Shmin
    # is a fit to the Rosemanowes hydrofracture tests of Pine et al. (1983b),
    # "to depths of 2000 m", so that is where the data end; SHmax is derived
    # from it at mu 0.8.
    # Pore pressure 9.494 MPa/km below a static fluid level 61 m down.
    stress_cases=(StressProfile(
        label="Reinecker et al. (2021)",
        Sv_grad=25.275e3, Shmin_grad=13.21e3, Shmin_int=3.0e6,
        SHmax_grad=25.99e3, SHmax_int=5.9e6, Pp_grad=9.494e3, Pp_datum=61.0,
        source="Reinecker et al. (2021), Geothermics 97, 102226",
        z_data=(0.0, 2000.0)),),
    stress_basis="measured",
    # C8, src/calibrate_ud1.py against the BGS UD-1 image log, wall cooling
    # bracketed 0-40 K (decided 3 October): intact rock (no breakout to 4 km)
    # and the breakout zones' median, at both ends of the bracket. UCS=180e6
    # above stays as v1.0's uncalibrated value.
    strength_cases=(("intact, 0 K", 203e6), ("intact, 40 K", 172e6),
                    ("weak zones, 0 K", 147e6), ("weak zones, 40 K", 118e6)))

PANNONIAN = SiteProfile(
    name="Pannonian Basin (Hungary)",
    geotherm=make_geotherm([(0, 15, 0.04625), (4000, 200, 0.035)]),
    target_depth=9714.0, target_T=400.0,
    Sv_grad=2550 * 9.81, K0_min=0.60, SHmax_over_Sv=0.90,   # back-arc, moderate
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.7, E_rock=50e9, UCS=160e6,
    # No measured magnitudes found. The basin interior is strike-slip, locally
    # transtensional, today (Bada et al. 2007, Foldtani Kozlony 137(3)); the
    # Miocene extension behind v1.0's normal-fault ratios is no longer active.
    # Bounded as transtensional, SHmax ~ Sv (decided 3 October).
    stress_cases=transitional_bounds(
        2550 * 9.81, source="regime: strike-slip, locally transtensional (Bada et al. 2007)"),
    stress_basis="regime bounds")

SITES = [SOULTZ, LARDERELLO, CORNWALL, PANNONIAN]


def classify(o):
    """(survives, stability, verdict).

    v1.0 adapter: the v1.0 rule, breakout fully suppressed below Shmin.
    v1.1 (C6): the mud-window verdict on breakout width, worst case in the
    site's stress range, with the range appended when the cases disagree. A
    stress state over the frictional cap (C2) can't be GO.
    """
    surv = o["m1"]["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP
    if o.get("v10"):
        return (surv,) + _classify_v10(o)
    w = o["stress"]["window"]
    if w["verdict"] == "GO":
        stab = f"OK ({w['width_hydro']:.0f} deg)"
    elif w["verdict"] == "CONDITIONAL":
        stab = f"+{w['overbalance_MPa']:.0f}MPa ({w['SG_lo']:.2f} SG)"
    else:
        stab = "NO (window shut)"
    verdict = w["verdict"]
    if verdict != "NO-GO" and o["T_rock"] >= 420:
        verdict = "GO* (ductile-drill)"
    verdicts = {c["window"]["verdict"] for c in o["stress_cases"]}
    if len(verdicts) > 1:
        lo, hi = min(verdicts, key=se.VERDICT_ORDER.index), se.worst_verdict(verdicts)
        verdict += f" [{lo}..{hi}]"
    if not o["stress"]["admissible"]["admissible"]:
        verdict = verdict.replace("GO", "CONDITIONAL", 1) if verdict.startswith("GO") else verdict
        verdict += " (stress inadmissible)"
    return surv, stab, verdict


def classify_sensitivity(o, key):
    """A sensitivity verdict: "C7 off" (no wall thermal stress) or the
    site-calibrated "W_max 63" (D2). None under the v1.0 adapter."""
    s = (o.get("sensitivity") or {}).get(key)
    if s is None:
        return None
    return classify(dict(o, stress=s["stress"], stress_cases=s["stress_cases"]))[2]


def classify_no_thermal(o):
    """The C7 sensitivity: the verdict with the wall thermal stress off."""
    return classify_sensitivity(o, "C7 off")


def _classify_v10(o):
    b = o["breakout"]
    if b["frac_limited"]:
        stab = "NO (frac-limited)"
    elif not b["breaks"]:
        stab = "OK"
    else:
        stab = f"+{b['overbalance_MPa']:.0f}MPa mud"
    bdt = o["T_rock"] >= 400  # at/over brittle-ductile edge
    if b["frac_limited"]:
        verdict = "NO-GO (stability)"
    elif bdt and o["T_rock"] >= 420:
        verdict = "GO* (ductile-drill)"
    elif b["breaks"]:
        verdict = "CONDITIONAL"
    else:
        verdict = "GO"
    return stab, verdict


def print_stress_cases(o):
    """Stress basis, data coverage, admissibility (C2), mud window (C6) and the
    range of stress cases."""
    s, ref = o["site"], o["stress"]
    beyond = ref["beyond_data"]
    cover = ("no magnitude data" if beyond is None else
             f"{beyond/1000:.1f} km below the deepest data" if beyond > 0 else
             "within the data")
    a = ref["admissible"]
    print(f"    stress: {s.stress_basis}, {ref['profile'].source}; target {cover}")
    print(f"    admissibility at mu {a['mu']}: effective S1/S3 {a['ratio']:.2f} vs cap "
          f"{a['cap']:.2f} -> {'admissible' if a['admissible'] else 'INADMISSIBLE'}"
          + ("; SHmax capped at the frictional limit" if ref["cap_binds"] else ""))
    print(f"    wall {o['T_wall']:.0f} C (Model 1, {o['m_min']} kg/s) against rock "
          f"{o['T_rock']:.0f} C: thermal hoop stress {o['dsigma_T']/1e6:+.0f} MPa"
          + ("" if o["thermal"] else " (C7 OFF)"))
    w = ref["window"]
    print(f"    mud window (W_max {w['W_max']:.0f} deg, wall {w['T_wall']:.0f} C): "
          f"{w['SG_lo']:.2f}-{w['SG_hi']:.2f} SG "
          f"({w['Pw_lo']/1e6:.0f}-{w['Pw_hi']/1e6:.0f} MPa){'' if w['open'] else ', SHUT'}; "
          f"breakout {w['width_hydro']:.0f} deg at hydrostatic, "
          f"{w['width_hi']:.0f} deg at the heaviest mud")
    tf = ("at every mud weight above pore pressure" if w["Pw_tensile"] <= ref["Pp"] else
          f"above {w['Pw_tensile']/1e6:.0f} MPa ({w['Pw_tensile']/(C.MUD_SG_GRAD*o['z']):.2f} SG)")
    print(f"    tensile fractures (T0 {C.WALL_T0/1e6:.0f} MPa) initiate {tf}; reported, not a bound")
    if len(o["stress_cases"]) > 1:
        for c in o["stress_cases"]:
            cw = c["window"]
            mark = "*" if c is ref else " "
            print(f"     {mark} {c['label']:<44} Shmin {c['Shmin']/1e6:5.0f} "
                  f"SHmax {c['SHmax']/1e6:5.0f}  window {cw['SG_lo']:.2f}-{cw['SG_hi']:.2f} SG  "
                  f"{cw['width_hydro']:3.0f} deg at hydro  {cw['verdict']}")
        print("       (* = worst case in range, used for the table row)")


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    results = [evaluate(s) for s in SITES]

    hdr = (f"{'Site':<26}{'depth':>7}{'T_rock':>7}{'K0':>5}{'aniso':>6}"
           f"{'bitC':>6}{'MW':>5}{'ROP':>5}{'gain':>5}{'stability':>20}  {'verdict':<34}"
           f"{'C7 off':<34}{'W_max ' + format(C.BREAKOUT_W_MAX_SITE_DEG, '.0f')}")
    print("=" * len(hdr))
    print("COMPARATIVE SITING TABLE  (target 400 C, vacuum tubing, quench-assist)")
    print("=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for s, o in zip(SITES, results):
        surv, stab, verdict = classify(o)
        d = o["drill"]
        short = s.name.split("(")[0].strip()[:25]
        print(f"{short:<26}{o['z']/1000:6.1f}k{o['T_rock']:7.0f}{o['K0']:5.2f}"
              f"{o['anisotropy']:6.2f}{o['m1']['T_bottom_delivered']:6.0f}"
              f"{o['MW_prod']:5.1f}{d['ROP']*3600:5.1f}{o['rop_gain']:5.1f}"
              f"{stab:>20}  {verdict:<34}{classify_no_thermal(o) or '':<34}"
              f"{classify_sensitivity(o, 'W_max ' + format(C.BREAKOUT_W_MAX_SITE_DEG, '.0f')) or ''}")
    print("-" * len(hdr))
    print("notes: MW at 10 kg/s early-life; ROP m/hr; gain = quench ROP multiplier;")
    print(f"       stability: breakout lobe width at hydrostatic mud, or the overbalance")
    print(f"       (and mud SG) that brings it within {C.BREAKOUT_W_MAX_DEG:.0f} deg below Shmin less")
    print(f"       {C.MUD_MARGIN_SG} SG; [a..b] = verdict range across the site's stress cases.")
    print("       C7 off = sensitivity: the verdict without the wall's thermal hoop stress.")
    print(f"       W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f} = sensitivity: the widest breakout UD-1 logged in a")
    print("       trouble-free section (site-calibrated, D2), in place of 90 deg.")
    print()
    # one-line readout per site
    for s, o in zip(SITES, results):
        b = o["breakout"]
        print(f"* {s.name}:")
        print(f"    depth {o['z']/1000:.1f} km to {o['T_rock']:.0f} C; "
              f"Sv {o['Sv']/1e6:.0f} / SHmax {o['SHmax']/1e6:.0f} / Shmin {o['Shmin']/1e6:.0f} MPa; "
              f"breakout sig_th {b['sigma_theta']/1e6:.0f} vs MC {b['mc_cold']/1e6:.0f} MPa, "
              f"P_need {b['P_need']/1e6:.0f} vs Shmin {o['Shmin']/1e6:.0f} MPa")
        print_stress_cases(o)

