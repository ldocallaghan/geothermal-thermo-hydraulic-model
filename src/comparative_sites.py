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
    # Bounded as transtensional, SHmax ~ Sv, as the interior is described; the
    # full strike-slip range (SHmax up to the friction cap) spans GO to NO-GO
    # at 9.7 km and says nothing about the site.
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

SITES = [SOULTZ, LARDERELLO, CORNWALL, PANNONIAN]


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
    """A sensitivity verdict: "no wall cooling" (the wall's thermal hoop
    stress switched off) or "W_max 63" (UD-1's site-calibrated breakout limit
    in place of 90 deg). None under the v1.0 adapter."""
    s = (o.get("sensitivity") or {}).get(key)
    if s is None:
        return None
    return classify(dict(o, stress=s["stress"], stress_cases=s["stress_cases"]))[2]


def classify_no_thermal(o):
    """The verdict with the wall's thermal hoop stress switched off."""
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
           f"{'no wall cooling':<34}{w63:<34}{'beyond data (stress / T)'}")
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
              f"{classify_sensitivity(o, w63) or '':<34}{beyond}")
        if s.tier == "speculative":
            gaps = ", ".join(f"{k} ({b})" for k, b in s.missing())
            print(f"{'':<26}   missing for evidence-based: {gaps}")
    print("-" * len(hdr))


def print_basis(site):
    """Each input's data basis and source."""
    for k in ("stress", "strength", "pore pressure", "temperature", "well check"):
        basis, src = site.data_basis.get(k, ("none", ""))
        print(f"    {k + ':':<15}{basis:<16}{src}")


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    results = [evaluate(s) for s in SITES]
    tiers = {t: [(s, o) for s, o in zip(SITES, results) if s.tier == t]
             for t in ("evidence-based", "speculative")}

    print_table("EVIDENCE-BASED SITES  (target 400 C; stress and strength measured or "
                "calibrated at the site, stability checked against a well there)",
                *zip(*tiers["evidence-based"]))
    print()
    print_table("SPECULATIVE SITES  (inputs missing, so the verdict rests on assumptions; "
                "NOT comparable with the table above)", *zip(*tiers["speculative"]))
    print("notes: MW at 10 kg/s early-life; ROP m/hr; gain = quench ROP multiplier;")
    print(f"       stability: breakout lobe width at hydrostatic mud, or the overbalance")
    print(f"       (and mud SG) that brings it within {C.BREAKOUT_W_MAX_DEG:.0f} deg below Shmin less")
    print(f"       {C.MUD_MARGIN_SG} SG; [a..b] = verdict range across the site's stress and strength cases.")
    print("       no wall cooling = sensitivity: the verdict without the wall's thermal hoop stress.")
    print(f"       W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f} = sensitivity: the widest breakout UD-1 logged in a")
    print("       trouble-free section (site-calibrated), in place of 90 deg.")
    print("       beyond data: how far the target lies below the deepest stress / temperature data.")
    print()
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
