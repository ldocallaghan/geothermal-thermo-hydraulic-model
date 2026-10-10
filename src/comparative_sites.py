"""
comparative_sites.py
===================
Run the site-evaluation tool across four European candidate provinces for a
400 C target, and report them in two tiers by what their verdicts rest on.

Evidence-based: stress magnitudes and rock strength measured or calibrated at
the site, and the stability model checked against a well there.
  - Upper Rhine Graben / Soultz (FR): stresses and lab UCS from Valley & Evans
    (2007); checked against GPK3/GPK4 (crosscheck_soultz.py).
  - United Downs / Carnmenellis (UK): stresses and pore pressure from
    Reinecker et al. (2021); strength calibrated on UD-1 (calibrate_ud1.py).
Speculative: stresses bounded by the faulting regime, strength unsourced.
  - Larderello (IT): temperature measured in Venelle-2 (Bertani et al. 2018);
    regime from Liotta & Brogi.
  - Pannonian Basin (HU): regime from Bada et al. (2007); regional gradient.

Every input's basis and source is in its SiteProfile.data_basis; the sources
are written up in data/sites/stress_sources.md.
"""
import sys

import numpy as np
import fluids
import geo_constants as C
import model1_coupled as m1
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
    geotherm=make_geotherm([(0, 15, 0.152), (2200, 350, 0.077)]),   # v1.0; see geotherm_v11
    target_depth=2850.0, target_T=400.0,
    Sv_grad=2600 * 9.81, K0_min=0.55, SHmax_over_Sv=0.75,   # normal-fault, low aniso
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.6, E_rock=45e9, UCS=140e6,
    # No measured magnitudes found. Regime: normal faulting alongside left-
    # lateral strike-slip, the intermediate stress switching between vertical
    # and horizontal (Liotta & Brogi, fault-slip data and focal mechanisms in
    # the Lago Basin), i.e. SHmax ~ Sv. Pore pressure is assumed hydrostatic;
    # the evidence points both ways (a vapour-dominated reservoir, but
    # over-pressured fluids along active faults).
    stress_cases=transitional_bounds(
        2600 * 9.81, source="regime: normal/strike-slip transition (Liotta & Brogi)"),
    stress_basis="regime bounds",
    # Venelle-2 (Bertani et al. 2018, SGP-TR-213): 350 C at 2.2 km, >= 504 C
    # at 2,815 m, 507-517 C at ~2.9 km. v1.0's profile read 400 C at 2.85 km,
    # more than 100 K too cold there; 400 C is at ~2.4 km.
    geotherm_v11=make_geotherm([(0, 15, 0.152), (2200, 350, 0.2504), (2815, 504, 0.1)]),
    target_depth_v11=2400.0,
    temperature_data_to=2900.0,
    data_basis={
        "stress": ("regime bounds", "Liotta & Brogi (manuscript, Geothermics); no "
                   "magnitudes published"),
        "strength": ("unsourced", "v1.0 140 MPa; only an Elba micaschist analogue found"),
        "pore pressure": ("assumed", "hydrostatic; vapour-dominated field, with "
                          "over-pressured fluids reported along active faults"),
        "temperature": ("measured", "Venelle-2 to 2.9 km (Bertani et al. 2018)"),
        "well check": ("none", "Venelle-2 drilled to 2.9 km (losses, stuck pipe), but no "
                       "breakout data published"),
    })

CORNWALL = SiteProfile(
    name="United Downs / Carnmenellis (UK)",
    geotherm=make_geotherm([(0, 15, 0.035), (5000, 190, 0.028)]),   # v1.0; see geotherm_v11
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
    # Calibrated by src/calibrate_ud1.py against the BGS UD-1 image log:
    # intact rock (no breakout to 4 km) and the breakout zones' median, each at
    # both ends of the 0-40 K bracket on wall cooling while UD-1 was drilled
    # (no mud temperatures are published). UCS=180e6 above stays as v1.0's
    # uncalibrated value.
    strength_cases=(("intact, 0 K", 203e6), ("intact, 40 K", 172e6),
                    ("weak zones, 0 K", 147e6), ("weak zones, 40 K", 118e6)),
    # Reinecker et al. (2021) give ~180 C at 5 km; v1.0's 190 C was 10 K
    # hot. v1.0's 28 C/km below 5 km is kept (unsourced), which puts 400 C at
    # ~12.86 km instead of 12.5 km.
    geotherm_v11=make_geotherm([(0, 15, 0.033), (5000, 180, 0.028)]),
    target_depth_v11=5000.0 + (400.0 - 180.0) / 0.028,
    temperature_data_to=5058.0,
    data_basis={
        "stress": ("measured", "Reinecker et al. (2021); Shmin data to 2.0 km, "
                   "SHmax derived at mu 0.8"),
        "strength": ("calibrated", "UD-1 breakouts, BGS image log (calibrate_ud1.py)"),
        "pore pressure": ("measured", "Reinecker et al. (2021)"),
        "temperature": ("extrapolated", "~180 C at 5 km (Reinecker et al. 2021); "
                        "28 C/km assumed below"),
        "well check": ("calibrated", "UD-1 to 5,058 m TVD (calibrate_ud1.py)"),
    })

PANNONIAN = SiteProfile(
    name="Pannonian Basin (Hungary)",
    geotherm=make_geotherm([(0, 15, 0.04625), (4000, 200, 0.035)]),
    target_depth=9714.0, target_T=400.0,
    Sv_grad=2550 * 9.81, K0_min=0.60, SHmax_over_Sv=0.90,   # back-arc, moderate
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.7, E_rock=50e9, UCS=160e6,
    # No measured magnitudes found. The basin interior is strike-slip, locally
    # transtensional, today (Bada et al. 2007, Foldtani Kozlony 137(3)); the
    # Miocene extension behind v1.0's normal-fault ratios is no longer active.
    # Bounded as transtensional, SHmax ~ Sv, as the interior is described. The
    # full strike-slip range (SHmax up to the friction cap) is wider: at 9.7 km
    # it spans GO to CONDITIONAL with the wall cooled, and CONDITIONAL to NO-GO
    # without cooling.
    stress_cases=transitional_bounds(
        2550 * 9.81, source="regime: strike-slip, locally transtensional (Bada et al. 2007)"),
    stress_basis="regime bounds",
    data_basis={
        "stress": ("regime bounds", "Bada et al. (2007); orientation maps only "
                   "(Bekesi et al. 2023), no magnitudes"),
        "strength": ("unsourced", "v1.0 160 MPa; no basement data found"),
        "pore pressure": ("assumed", "hydrostatic, but the deep Great Hungarian Plain "
                          "is reported overpressured by 1-35 MPa (not yet verified)"),
        "temperature": ("regional", "45-50 C/km regional gradient, extrapolated to "
                        "9.7 km; no deep well cited"),
        "well check": ("none", "no well with wellbore-failure data"),
    })

import crosscheck_newberry as nb                           # noqa: E402

NEWBERRY = SiteProfile(
    name="Newberry / NWG 55-29 (USA)",
    geotherm=nb.geotherm, target_depth=nb.target_depth(), target_T=400.0,
    Sv_grad=nb.profile().Sv_grad, K0_min=0.53, SHmax_over_Sv=0.85,   # ratios unused
    rho_fluid_grad=C.HYDROSTATIC_GRAD,
    # mean of the corrected conductivities over the open hole, 6,500-10,000 ft
    # (AltaRock, "Well 55-29 Heat Flow Values", GDR 271)
    k_rock=2.14, E_rock=50e9, UCS=79e6,
    # Davatzes and Hickman (2011) at both their friction coefficients; the data
    # cover the imaged open hole, 6,435-8,860 ft MD
    stress_cases=(nb.profile(0.55), nb.profile(0.70)),
    stress_basis="measured range",
    strength_cases=nb.granodiorite_strengths(),
    temperature_data_to=nb.Z_DATA_T,
    mu=0.55,
    data_basis={
        "stress": ("measured range", "Davatzes & Hickman (2011): Sv and pore pressure "
                   "measured; Shmin at frictional equilibrium (mu 0.55, 0.70); SHmax from "
                   "breakout widths, to 8,860 ft MD"),
        "strength": ("log-derived", "UCS-porosity relation of Davatzes & Hickman on the "
                     "neutron log, granodiorite 8,807 ft to TD (GDR 271)"),
        "pore pressure": ("measured", "static survey, October 2008 (GDR 271); "
                          "underpressured, Pf/Sv ~0.34"),
        "temperature": ("extrapolated", f"{nb.T_DATA_END:.0f} C at {nb.Z_DATA_T:.0f} m TVD "
                        f"(static survey 2008); {nb.GRAD_DEEP*1000:.0f} C/km below"),
        "well check": ("not reproduced", "NWG 55-29 (crosscheck_newberry.py): the volcanic "
                       "breakout width is reproduced at 25-35 K of wall cooling, but "
                       "breakouts are predicted in the granodiorite, where none were logged"),
    })

SITES = [SOULTZ, LARDERELLO, CORNWALL, PANNONIAN, NEWBERRY]


def classify(o):
    """(survives, stability, verdict).

    v1.0 adapter: the v1.0 rule, breakout fully suppressed below Shmin.
    v1.1: the mud-window verdict on breakout width, for the worst of the
    site's stress and strength cases, with the range appended when the cases
    disagree. A stress state over the frictional cap can't be GO.
    """
    surv = o["m1"]["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP
    if o.get("v10"):
        return (surv,) + _classify_v10(o)
    w = o["stress"]["window"]
    if w["verdict"] == "GO":
        stab = f"OK ({w['width_hydro']:.0f} deg)"
    elif w["verdict"] == "CONDITIONAL":
        ob = w["overbalance_MPa"]
        stab = f"+{ob:.1f}MPa ({w['SG_lo']:.2f} SG)" if ob < 10 else f"+{ob:.0f}MPa ({w['SG_lo']:.2f} SG)"
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
    """A sensitivity verdict: "no wall cooling" (the wall at rock
    temperature: no thermal hoop stress and no strength gain from cooling),
    "W_max 63" (UD-1's site-calibrated breakout limit in place of 90 deg), or
    "no wall cooling, W_max 63" (both). None under the v1.0 adapter."""
    s = (o.get("sensitivity") or {}).get(key)
    if s is None:
        return None
    return classify(dict(o, stress=s["stress"], stress_cases=s["stress_cases"]))[2]


def classify_no_thermal(o):
    """The verdict for an uncooled wall."""
    return classify_sensitivity(o, "no wall cooling")


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
    """Stress basis, data coverage, frictional admissibility, the mud window,
    and the range of stress and strength cases."""
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
          + ("" if o["thermal"] else " (switched off)"))
    w = ref["window"]
    print(f"    mud window (W_max {w['W_max']:.0f} deg, wall {w['T_wall']:.0f} C): "
          f"{w['SG_lo']:.2f}-{w['SG_hi']:.2f} SG "
          f"({w['Pw_lo']/1e6:.0f}-{w['Pw_hi']/1e6:.0f} MPa){'' if w['open'] else ', SHUT'}; "
          f"at hydrostatic mud {se.describe_width(w['width_hydro'])}, "
          f"at the heaviest {se.describe_width(w['width_hi'])}")
    tf = ("at every mud weight above pore pressure" if w["Pw_tensile"] <= ref["Pp"] else
          f"above {w['Pw_tensile']/1e6:.0f} MPa ({w['Pw_tensile']/(C.MUD_SG_GRAD*o['z']):.2f} SG)")
    print(f"    tensile fractures (T0 {C.WALL_T0/1e6:.0f} MPa) initiate {tf}; reported, not a bound")
    if len(o["stress_cases"]) > 1:
        for c in o["stress_cases"]:
            cw = c["window"]
            mark = "*" if c is ref else " "
            print(f"     {mark} {c['label']:<44} Shmin {c['Shmin']/1e6:5.0f} "
                  f"SHmax {c['SHmax']/1e6:5.0f}  window {cw['SG_lo']:.2f}-{cw['SG_hi']:.2f} SG  "
                  f"{se.describe_width(cw['width_hydro']):>21} at hydro  {cw['verdict']}")
        print("       (* = worst case in range, used for the table row)")


def beyond_data_km(o):
    """(stress, temperature): km by which the target lies below the deepest
    data; None where there is no data at all."""
    s = o["site"]
    stress = o["stress"]["beyond_data"]
    temp = (None if s.temperature_data_to is None
            else max(o["z"] - s.temperature_data_to, 0.0))
    return (None if stress is None else stress / 1000.0,
            None if temp is None else temp / 1000.0)


def print_table(title, sites, results):
    w63 = "W_max " + format(C.BREAKOUT_W_MAX_SITE_DEG, ".0f")
    hdr = (f"{'Site':<26}{'depth':>7}{'T_rock':>7}{'K0':>5}{'aniso':>6}"
           f"{'bitC':>6}{'MW':>5}{'ROP':>5}{'gain':>5}{'stability':>20}  {'verdict':<34}"
           f"{'no wall cooling':<34}{w63:<34}{'no cooling, ' + w63:<34}"
           f"{'beyond data (stress / T)'}")
    print("=" * len(hdr))
    print(title)
    print("=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for s, o in zip(sites, results):
        surv, stab, verdict = classify(o)
        d = o["drill"]
        short = s.name.split("(")[0].strip()[:25]
        bs, bt = beyond_data_km(o)
        beyond = (f"{'no data' if bs is None else f'{bs:.1f} km'} / "
                  f"{'no data' if bt is None else f'{bt:.1f} km'}")
        print(f"{short:<26}{o['z']/1000:6.1f}k{o['T_rock']:7.0f}{o['K0']:5.2f}"
              f"{o['anisotropy']:6.2f}{o['m1']['T_bottom_delivered']:6.0f}"
              f"{o['MW_prod']:5.1f}{d['ROP']*3600:5.1f}{o['rop_gain']:5.1f}"
              f"{stab:>20}  {verdict:<34}{classify_no_thermal(o) or '':<34}"
              f"{classify_sensitivity(o, w63) or '':<34}"
              f"{classify_sensitivity(o, 'no wall cooling, ' + w63) or '':<34}{beyond}")
        if s.tier == "speculative":
            gaps = ", ".join(f"{k} ({b})" for k, b in s.missing())
            print(f"{'':<26}   missing for evidence-based: {gaps}")
    print("-" * len(hdr))


def classify_circ_sensitivity(o, key):
    """A verdict with the wall temperature from a circulation sensitivity:
    "dual-wall pipe (not manufactured)", "static fracture bound" (no
    circulating friction on the upper bound), "rock exposed 1 yr", "film x0.75"
    (FORGE's fitted annulus film) or "vacuum-insulated pipe"."""
    s = (o.get("circ_sensitivity") or {}).get(key)
    if s is None:
        return None
    return classify(dict(o, stress=s["stress"], stress_cases=s["stress_cases"]))[2]


def change_reason(o, old):
    """Why the verdict moved from the earlier circulation model, in words."""
    new_v, old_v = classify(o)[2], classify(old)[2]
    c = o["circ"]
    t = next(t for t in c["tried"] if t["pipe"] == c["pipe"].name)
    set_by = ("hole cleaning sets the flow" if (t["cleaning_flow"] or 0) > (t["survival_flow"] or 0)
              else "tool survival sets the flow")
    wall = f"wall {old['T_wall']:.0f} -> {o['T_wall']:.0f} C"
    if new_v == old_v:
        return f"unchanged; {wall}"
    return f"{wall}: {c['pipe'].name} pipe at {c['m_dot']} kg/s ({set_by})"


def print_circulation(title, sites, results, olds):
    lim, vmin = m1.PUMP_LIMIT / 1e6, m1.HOLE_CLEANING_V
    hdr = (f"{'Site':<26}{'pipe':>17}{'flow':>6}{'gpm':>6}{'SPP':>6}{'ECD':>6}{'v_ann':>6}"
           f"{'BHCT':>6}{'1 yr':>6}{'heat':>6}  {'verdict':<30}{'GO at':<26}{'earlier model':<30}")
    print("=" * len(hdr))
    print(title)
    print("=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for s, o, old in zip(sites, results, olds):
        c = o["circ"]
        h = c["hyd"]
        short = s.name.split("(")[0].strip()[:25]
        spp = f"{h['spp']/1e6:.1f}" + ("!" if h["over_pump_limit"] else "")
        va = f"{h['v_ann_min']:.2f}" + ("" if h["cleans_hole"] else "!")
        g = c["go"]
        go = ("already GO" if classify(o)[2].startswith("GO") else
              "none within the limit" if g is None else
              f"{g['m_dot']:.0f} kg/s, {g['gpm']:.0f} gpm, {g['spp']/1e6:.0f} MPa"
              + (", within rig" if g["within_rig"] else ", beyond rig"))
        print(f"{short:<26}{c['pipe'].name:>17}{c['m_dot']:6.0f}{c['gpm']:6.0f}{spp:>6}"
              f"{c['ecd_SG']:6.3f}{va:>6}{c['bhct']:6.0f}{c['bhct_one_year']:6.0f}"
              f"{o['MW_prod']:6.1f}  {classify(o)[2]:<30}{go:<26}{classify(old)[2]:<30}")
        print(f"{'':<26}   against the earlier model: {change_reason(o, old)}")
        if c["conflict"]:
            print(f"{'':<26}   conflict: {c['conflict']}")
    print("-" * len(hdr))
    print(f"  flow kg/s and gal/min; SPP standpipe pressure, MPa (limit {lim:.1f}, ! over); ECD the")
    print("  annular friction as an equivalent density, SG, added to the mud on the fracture bound;")
    print(f"  v_ann lowest annular velocity, m/s (hole-cleaning reference {vmin:.3f}, ! under); BHCT")
    print("  bottom-hole circulating temperature, C, the wall temperature for the verdict; 1 yr: the")
    print("  same with the rock exposed for a year; heat: heat returned while drilling, early life,")
    print("  single loop, MW. GO at: the lowest flow at which the same pipe gives GO within the pump")
    print(f"  limit, against what {se.RIG_PUMPS} pumps rated as NOV's 2,200 hp 14-P-220 deliver at that")
    print("  standpipe pressure. Earlier model: the single 0.02 W/m K pipe at the lowest survivable flow.")
    print()
    keys = ("dual-wall pipe (not manufactured)", "static fracture bound", "rock exposed 1 yr",
            f"film x{se.FILM_MULT_FORGE}", "vacuum-insulated pipe")
    print("Circulation sensitivities (verdict, wall temperature)")
    for s, o in zip(sites, results):
        short = s.name.split("(")[0].strip()[:25]
        print(f"  {short}")
        for k in keys:
            v = o["circ_sensitivity"][k]
            print(f"    {k:<36}{classify_circ_sensitivity(o, k)} ({v['T_wall']:.0f} C)")
    print()


def print_cycle(title, sites, results):
    """The verdict through the drilling cycle: while drilling, at a connection
    and through a trip, with the safe pause, the trip and its fluid weight."""
    hdr = (f"{'Site':<26}{'drilling':>13}{'connection':>13}{'fresh rock':>13}{'trip':>13}{'drill SG':>9}"
           f"{'limit':>7}{'safe pause':>11}{'trip h':>7}{'trip SG':>8}{'limit':>7}  {'over the cycle':<15}"
           f"{'while circulating (v1.2)'}")
    print("=" * len(hdr))
    print(title)
    print("=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    fmt_pause = lambda h: "> 7 d" if h == float("inf") else f"{h * 60:.0f} min" if h < 1 else f"{h:.1f} h"
    for s, o in zip(sites, results):
        c = o["cycle"]
        short = s.name.split("(")[0].strip()[:25]
        st = c["verdicts"]
        print(f"{short:<26}{st['drilling']:>13}{st['connection']:>13}{st['connection, fresh rock']:>13}"
              f"{st['trip']:>13}{c['drill_SG']:9.2f}{c['drill_SG_hi']:7.2f}{fmt_pause(c['safe_pause_h']):>11}"
              f"{c['trip_h']:7.0f}"
              f"{c['trip_SG']:8.2f}{c['trip_SG_hi']:7.2f}  {c['site_verdict']:<15}{classify(o)[2]}")
    print("-" * len(hdr))
    print("  fresh rock: a connection on the rock just drilled, after FORGE's 10th-percentile circulation")
    print("  between the last new hole and the pumps stopping. drill SG: the static fluid weight that")
    print("  holds the wall as the bit exposes it (with the circulating friction) and through both")
    print("  connections, against its limit (the fracture limit less the friction); safe pause: how")
    print("  long from the end of drilling the hole holds at that weight with no circulation; trip h:")
    print("  the bottom of the hole without circulation for a trip, from FORGE 16B's tripping speeds")
    print("  and routine surface time scaled to depth; trip SG: the static fluid the hole needs for")
    print("  the trip, against the fracture limit (Shmin less 0.05 SG). With the pressure budget the")
    print("  SGs are equivalent static densities at the target, the trip SG includes the swab")
    print("  allowance and its limit the surge allowance, and the safe pause lets the column heat.")
    if any(o.get("budget") for o in results):
        print_budget(sites, results)
        print_windows(sites, results)
        print_tools(sites, results)
        print_status(sites, results)
    print()
    print("Cycle sensitivities (over the cycle; drill SG; safe pause; trip SG)")
    for s, o in zip(sites, results):
        short = s.name.split("(")[0].strip()[:25]
        print(f"  {short}")
        for k, c in o["cycle_sensitivity"].items():
            print(f"    {k:<40}{c['site_verdict']:<28}{c['drill_SG']:6.2f}{fmt_pause(c['safe_pause_h']):>10}"
                  f"{c['trip_SG']:8.2f}")
        print(f"    {'staged circulation':<40}no change at the bottom: staging cools only the hole above the bit")
        print(f"    {'double bit run':<40}no change at the bottom: the element one stand up is the same")
    print()
    pct = lambda x: "already NO-GO" if x is None else ("none within 20%" if x == float("inf") else f"{100 * x:.1f}%")
    print("Margins over the cycle: the flow range (from the flow used to the highest flow at which the")
    print("verdict holds), and the fall in Shmin that closes the window, with SHmax held or scaled with")
    print("it, with the frictional cap on SHmax at the site's friction coefficient and at 1.0; mu: the")
    print("friction the deciding case's stresses need at the first closure with SHmax held")
    for s_, o in zip(sites, results):
        short = s_.name.split("(")[0].strip()[:25]
        fc, sm = o["flow_ceiling"], o["shmin_margin"]
        mu = lambda d: "" if d["friction_needed"] is None else f" (mu {d['friction_needed']:.2f})"
        unit = "kg/s" if o.get("fluid") is None else "L/s"
        if fc["m_dot"] is None:
            print(f"  {short:<26}flow {o['circ']['m_dot']:.0f} {unit}; no flow ceiling: already NO-GO")
        else:
            print(f"  {short:<26}flow {o['circ']['m_dot']:.0f} to {fc['m_dot']:.0f} {unit}, ceiling set by"
                  f" {fc['by']} ({fc['gpm']:.0f} gpm, {'within' if fc['within_rig'] else 'beyond'} rig capacity)")
        for tag, name in (("site", f"cap at mu {s_.mu:.2f}"), ("mu_hi", "cap at mu 1.00")):
            d = sm[tag]
            print(f"  {'':<26}Shmin fall, {name}: SHmax held {pct(d['SHmax_held'])}{mu(d)}, "
                  f"scaled {pct(d['SHmax_scaled'])}")
    print()
    print("Failed depth at the Shmin azimuth from the full stress field, minutes after the bit exposes")
    print("the rock, at the drilling fluid with its circulating friction (deciding case), against the")
    print("uncooled depth whose temperature the width check uses")
    for s_, o in zip(sites, results):
        c = o["cycle"]
        short = s_.name.split("(")[0].strip()[:25]
        cells = "  ".join(f"{k['minutes']:.0f} min {100 * k['depth']:.1f} cm" for k in c["skin"])
        print(f"  {short:<26}uncooled {100 * c['skin_depth_uncooled']:.1f} cm;  {cells}")
    print()
    print("Trip verdict up the open hole (height above the bit; trip SG needed)")
    for s, o in zip(sites, results):
        short = s.name.split("(")[0].strip()[:25]
        cells = "  ".join(f"{r['height']:.0f} m {r['SG_needed']:.2f}" for r in o["trip_profile"])
        print(f"  {short:<26}{cells}")
    print()


def print_budget(sites, results):
    """The pressure budget: the drilling fluid each site converged on, the
    water column in each state, and what the trip fluid loses as the column
    heats through the trip."""
    print()
    print("Pressure budget: SGs are equivalent static densities at the target; the fluid is a")
    print("barite-weighted water-based mud, given by its density at surface (20 C)")
    hdr = (f"  {'Site':<26}{'fluid':>8}{'passes':>7}{'ECD':>7}{'water':>7}{'at trip':>8}"
           f"{'trip SG':>8}{'at start':>9}{'trip fluid':>11}{'lost':>9}{'return':>8}")
    print(hdr)
    for s, o in zip(sites, results):
        c, b = o["cycle"], o["cycle"]["budget"]
        short = s.name.split("(")[0].strip()[:25]
        fluid = "water" if o["fluid"] is None else f"{o['fluid'].sg:.3f}"
        trip_f = "water" if b["trip_rho_surface"] is None else f"{b['trip_rho_surface'] / 1000:.3f}"
        conv = "" if o["fluid_converged"] else "*"
        print(f"  {short:<26}{fluid:>8}{len(o['fluid_passes']):>6}{conv:1}{c['ecd_SG']:7.3f}"
              f"{b['water_SG']:7.3f}{b['water_SG_trip']:8.3f}{c['trip_SG']:8.3f}{b['trip_SG_start']:9.3f}"
              f"{trip_f:>11}{b['trip_dP_heating'] / 1e6:6.1f} MPa{o['circ']['T_return']:6.0f} C")
    print("  fluid: the drilling fluid's surface weight, iterated with the circulation it sets to")
    print(f"  within 0.002 (passes; * not converged). water / at trip: water's column at the")
    print("  circulating temperatures and at the end of the trip. trip SG / at start: the trip fluid")
    print("  at the end of the trip (with the 0.02 swab allowance) and at its start, against the")
    print("  fracture limit less the 0.02 surge allowance; lost: the bottomhole pressure the trip")
    print("  fluid loses as the column heats. return: the surface return; at 100 C the annulus")
    print("  would need backpressure.")
    hot = [s.name.split("(")[0].strip() for s, o in zip(sites, results) if o["circ"].get("mud_too_hot")]
    if hot:
        print(f"  flag: the circulating bottomhole temperature exceeds the {fluids.MUD_T_LIMIT:.0f} C at which")
        print("  conventional water-based muds start to break down: " + ", ".join(hot))


def print_windows(sites, results):
    """The fluid that suits every scenario in every state down the hole, as a
    density at surface, with the bound and depth that bind it."""
    print()
    print("Common windows: the fluid's density at surface (20 C) that suits every stress and strength")
    print("case in every state, down the open hole to the casing shoe; the bound that binds, the state,")
    print("case and measured depth that set it")
    f = lambda x: "none" if x == float("inf") else f"{x / 1000:.3f}"
    for s, o in zip(sites, results):
        short = s.name.split("(")[0].strip()[:25]
        for k, name in (("drilling_fluid", "drilling fluid"), ("trip_fluid", "trip fluid")):
            w = o["windows"][k]
            lo, hi = w["lo_by"], w["hi_by"]
            head = f"{f(w['fluid'])}" if w["open"] else w["note"]
            print(f"  {short:<26}{name:<15}{head}")
            print(f"  {'':<26}{'':<15}lower {f(w['lo'])}: {lo['state']}, {lo['case']}, {lo['md']:.0f} m")
            print(f"  {'':<26}{'':<15}upper {f(w['hi'])}: {hi['state']}, {hi['case']}, {hi['md']:.0f} m")
            short = ""
    print("  Each state's own window per case is in the result (windows/states). The site verdict")
    print("  above is unchanged: it is the worst case's; this is the prescription, a separate output.")


def print_status(sites, results):
    print()
    print("Status: INDETERMINATE where a prerequisite failed or is unknown, NO-GO with the failing parts,")
    print("CONDITIONAL with its conditions")
    for s, o in zip(sites, results):
        short = s.name.split("(")[0].strip()[:25]
        st = o["status"]
        print(f"  {short:<26}{st['label']:<15}" + "; ".join(st["reasons"]))


def print_tools(sites, results):
    """The tool (the BHA's internal fluid) and the coating through the cycle,
    the schedule each needs, and the site's status."""
    print()
    print(f"Tools and coating: the tool against {C.BHA_SURVIVAL_TEMP:.0f} C while drilling, through the connection")
    print("on fresh rock, and running back in with FORGE's staging (and without); the coating's")
    print("steel-side temperature against its rating")
    hdr = (f"  {'Site':<26}{'bit':>6}{'conn.':>7}{'in':>6}{'unstaged':>9}  {'tool':<34}"
           f"{'coating max':>12}  {'coating':<15}{'status'}")
    print(hdr)
    for s, o in zip(sites, results):
        tl = o["tools"]
        short = s.name.split("(")[0].strip()[:25]
        co = tl["coating"]
        cmax = max(co["circulating_max"], co["connection"]["max"], co["trip_in"]["max"])
        print(f"  {short:<26}{tl['drilling']['T']:6.0f}{tl['connection']['peak']:7.0f}"
              f"{tl['trip_in']['staged']['peak']:6.0f}{tl['trip_in']['unstaged']['peak']:9.0f}  "
              f"{tl['tool']:<34}{cmax:12.0f}  {co['status']:<15}{o['status']['label']}")
    print("  bit: the fluid delivered to the bit while drilling; conn.: the tool's peak through the")
    print("  connection; in / unstaged: its peak running back in with FORGE's staging and without.")
    print("Schedules the tool needs, where FORGE's practice fails")
    for s, o in zip(sites, results):
        tl = o["tools"]
        short = s.name.split("(")[0].strip()[:25]
        c, tr = tl["connection"], tl["trip_in"]
        if not tl["drilling"]["ok"]:
            print(f"  {short:<26}none: the tool is above {C.BHA_SURVIVAL_TEMP:.0f} C while circulating at the bit")
            continue
        if not c["ok"]:
            sc = c["schedule"]
            print(f"  {short:<26}connection: " + ("none short of continuous circulation" if sc is None else
                  f"{sc['precirc_min']:.1f} min of circulation before it, pumps off "
                  f"{sc['off_min']:.1f} min (tool {sc['peak']:.0f} C)"))
        if not tr["ok"]:
            sc = tr["schedule"]
            print(f"  {short:<26}running in: " + ("none short of continuous circulation" if sc is None else
                  f"staging every {sc['spacing_m']:.0f} m for {sc['minutes']:.0f} min, {sc['stages']} stages, "
                  f"{sc['staging_h']:.0f} h circulating ({100 * sc['circulating_fraction']:.0f}% of the run in); "
                  f"tool {sc['peak']:.0f} C"))
    print("  A schedule counts as found if it falls short of continuous circulation; the share of the")
    print("  run in spent circulating is given for the reader to judge whether it is practical.")
    print("Coating: hottest while circulating; through the connection; running in (minutes above rating)")
    for s, o in zip(sites, results):
        co = o["tools"]["coating"]
        short = s.name.split("(")[0].strip()[:25]
        above = lambda d: "" if d["minutes_above"] is None else f" ({d['minutes_above']:.0f} min above)"
        rating = "no rating" if co["rating"] is None else f"rating {co['rating']:.0f} C"
        print(f"  {short:<26}{co['pipe']}, {rating}: {co['circulating_max']:.0f} C; "
              f"{co['connection']['max']:.0f} C{above(co['connection'])}; "
              f"{co['trip_in']['max']:.0f} C{above(co['trip_in'])}")
    print("  Running in, staging is taken to cool the fluid in the string but not the column; its")
    print("  circulation also flushes the annulus around the coated pipe above the bit, so the time")
    print("  above the rating on the way in is overstated.")


def golden_rows(results):
    rows = {}
    for o in results:
        c = o["circ"]
        rows[o["site"].name] = dict(
            verdict=classify(o)[2], pipe=c["pipe"].name, m_dot=float(c["m_dot"]),
            bhct=float(c["bhct"]), bhct_one_year=float(c["bhct_one_year"]),
            spp_MPa=float(c["hyd"]["spp"]) / 1e6, v_ann=float(c["hyd"]["v_ann_min"]),
            ecd_SG=float(c["ecd_SG"]), MW=float(o["MW_prod"]),
            go_m_dot=None if c["go"] is None else float(c["go"]["m_dot"]),
            sensitivities={k: classify_sensitivity(o, k) for k in o["sensitivity"]}
            | {k: classify_circ_sensitivity(o, k) for k in o["circ_sensitivity"]})
    return rows


def golden_rows_v13(results):
    rows = {}
    for o in results:
        c = o["cycle"]
        rows[o["site"].name] = dict(
            verdicts=c["verdicts"], site_verdict=c["site_verdict"], drill_SG=c["drill_SG"],
            drill_SG_hi=c["drill_SG_hi"],
            safe_pause_h=None if c["safe_pause_h"] == float("inf") else c["safe_pause_h"],
            trip_h=c["trip_h"], trip_SG=c["trip_SG"], trip_SG_hi=c["trip_SG_hi"],
            T_ref=c["T_ref"],
            sensitivities={k: dict(site_verdict=v["site_verdict"], trip_SG=v["trip_SG"],
                                         drill_SG=v["drill_SG"])
                           for k, v in o["cycle_sensitivity"].items()},
            trip_profile=[(r["height"], r["verdict"], r["SG_needed"]) for r in o["trip_profile"]],
            flow_ceiling=dict(m_dot=o["flow_ceiling"]["m_dot"], by=o["flow_ceiling"]["by"]),
            shmin_margin={k: {kk: (None if vv is None else (-1.0 if vv == float("inf") else vv))
                              for kk, vv in d.items()} for k, d in o["shmin_margin"].items()})
    return rows


def print_basis(site):
    """Each input's data basis and source."""
    for k in ("stress", "strength", "pore pressure", "temperature", "well check"):
        basis, src = site.data_basis.get(k, ("none", ""))
        print(f"    {k + ':':<15}{basis:<16}{src}")


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    results = [evaluate(s) for s in SITES]
    olds = [evaluate(s, circulation="v1.1") for s in SITES]
    if "--write-golden" in sys.argv:
        import json, os
        gold = os.path.join(os.path.dirname(__file__), "..", "tests", "golden")
        with open(os.path.join(gold, "v12_site_table.json"), "w") as fh:
            json.dump(golden_rows(results), fh, indent=1)
        with open(os.path.join(gold, "v13_site_table.json"), "w") as fh:
            json.dump(golden_rows_v13(results), fh, indent=1)
    tiers = {t: [(s, o) for s, o in zip(SITES, results) if s.tier == t]
             for t in ("evidence-based", "speculative")}
    tiers_old = {t: [o for s, o in zip(SITES, olds) if s.tier == t]
                 for t in ("evidence-based", "speculative")}

    print_table("EVIDENCE-BASED SITES  (target 400 C; stress and strength measured or "
                "calibrated at the site, stability checked against a well there)",
                *zip(*tiers["evidence-based"]))
    print()
    print_table("SPECULATIVE SITES  (inputs missing, so the verdict rests on assumptions; "
                "NOT comparable with the table above)", *zip(*tiers["speculative"]))
    print("notes: MW = heat returned while drilling, early life, single loop; ROP m/hr;")
    print("       gain = quench ROP multiplier; bitC = the wall temperature for the verdict;")
    print(f"       stability: breakout lobe width at hydrostatic mud, or the overbalance")
    print(f"       (and mud SG) that brings it within {C.BREAKOUT_W_MAX_DEG:.0f} deg below Shmin less")
    print(f"       {C.MUD_MARGIN_SG} SG; [a..b] = verdict range across the site's stress and strength cases.")
    print("       no wall cooling = sensitivity: the wall at rock temperature, so neither the")
    print("       thermal hoop stress nor the strength gain from cooling.")
    print(f"       W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f} = sensitivity: the widest breakout UD-1 logged in a")
    print("       trouble-free section (site-calibrated), in place of 90 deg.")
    print("       beyond data: how far the target lies below the deepest stress / temperature data.")
    print()
    for tier in ("evidence-based", "speculative"):
        if tiers[tier]:
            print_cycle(f"THE DRILLING CYCLE, {tier.upper()} SITES: drilling, connection and trip",
                        *zip(*tiers[tier]))
    for tier in ("evidence-based", "speculative"):
        if tiers[tier]:
            print_circulation(f"CIRCULATION, {tier.upper()} SITES: the flow, pipe and wall "
                              "temperature behind each verdict, beside the earlier model's verdict",
                              *zip(*tiers[tier]), tiers_old[tier])
    for tier in ("evidence-based", "speculative"):
        for s, o in tiers[tier]:
            b = o["breakout"]
            print(f"* {s.name}  [{tier}]:")
            print_basis(s)
            print(f"    depth {o['z']/1000:.1f} km to {o['T_rock']:.0f} C; "
                  f"Sv {o['Sv']/1e6:.0f} / SHmax {o['SHmax']/1e6:.0f} / Shmin {o['Shmin']/1e6:.0f} MPa; "
                  f"breakout sig_th {b['sigma_theta']/1e6:.0f} vs MC {b['mc_cold']/1e6:.0f} MPa, "
                  f"P_need {b['P_need']/1e6:.0f} vs Shmin {o['Shmin']/1e6:.0f} MPa")
            print_stress_cases(o)
