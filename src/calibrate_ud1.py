"""
calibrate_ud1.py
===============
Calibrate rock-mass strength against the breakouts logged in UD-1, the
5,275 m MD well at United Downs.

Run from the repo root:   python src/calibrate_ud1.py

Inputs (all in data/ud1/):
  - stresses and pore pressure: Reinecker et al. (2021), via the United Downs
    profile in comparative_sites.CORNWALL
  - breakout depths, widths and heights: the BGS image-log interpretation,
    UD1_Borehole_Imaging_Interpretation.csv (provenance in bgs_image_log.md)
  - mud: below 1.05 SG in the 12.25" section (Reinecker et al. 2021, s. 5);
    1.05 SG is used, the upper bound, with 1.00 SG as a sensitivity
  - formation temperature: linear from the model's 15 C surface to the paper's
    "around 180 C at 5 km" (v1.0's site profile read ~10 K hot there)

The breakout check is the one the site runs use: Kirsch hoop stress, failure
in effective stress, with the thermal hoop stress of the cooled wall. No mud
or circulating temperatures are published for UD-1, so wall cooling is
bracketed at 0-40 K below the formation temperature and the spread is carried
into the strength range. 40 K is roughly the gap between static and
circulating bottom-hole temperature at 4 km in conventional drilling.

The calibration is closed-form. Each logged breakout of width W at depth z
fixes the strength of the rock it formed in: the failure criterion inverted
for sigma_cm, then un-derated to UCS at 25 C. Rock that did not break out puts a
lower bound on intact strength: enough to keep the width at zero at the
deepest such depth. The fit therefore separates weak or fractured zones (the
breakouts) from intact rock (everything else), and reports both.

Only the near-vertical 12.25" section (900-4,000 m MD, inclination <= 16 deg)
is fitted. The deviated 8.5" section is a qualitative check: a calibration
that makes it undrillable is wrong.
"""
import csv
import os
import sys

import numpy as np

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(ROOT, "data", "ud1", "UD1_Borehole_Imaging_Interpretation.csv")
FIG = os.path.join(ROOT, "figures", "ud1_calibration.png")

PROFILE = CORNWALL.stress_cases[0]          # Reinecker et al. (2021)
MUD_SG = 1.05                               # 12.25" section, upper bound
SECTION_MD = (900.0, 4000.0)                # 12.25" section
TD_TVD = 5058.0
T_SURF, T_5KM = C.SURFACE_TEMP, 180.0       # Reinecker et al. (2021)
COOLING_K = (0.0, 10.0, 20.0, 30.0, 40.0)   # wall below formation temperature
PHI_DEG = 35.0                              # model5.KMC; 30 and 40 as sensitivity


# ------------------------------------------------------------------ inputs
def load_log(path=LOG):
    """Image-log features as a list of dicts (two header rows: names, units)."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    names = rows[0]
    out = []
    for r in rows[2:]:
        d = dict(zip(names, r))
        out.append(dict(md=float(d["Measured Depth"]),
                        dev=float(d["Hole Deviation"]),
                        cls=d["Dip Classification"],
                        width=float(d["Breakout_Width_N"]),
                        height=float(d["Breakout_Height_N"])))
    return sorted(out, key=lambda f: f["md"])


def md_to_tvd(features):
    """MD -> TVD by the average-angle method over every logged station.

    The hole is taken as vertical above the first station (907 m MD, 1.5 deg).
    Over the 12.25" section the correction is at most 10 m (at 4,000 m MD).
    """
    md = np.array([f["md"] for f in features])
    dev = np.radians([f["dev"] for f in features])
    tvd = md[0] + np.concatenate(
        [[0.0], np.cumsum(np.diff(md) * np.cos(0.5 * (dev[1:] + dev[:-1])))])
    return lambda m: float(np.interp(m, md, tvd))


def T_formation(z):
    return T_SURF + (T_5KM - T_SURF) * z / 5000.0


def kmc(phi_deg):
    s = np.sin(np.radians(phi_deg))
    return (1 + s) / (1 - s)


# ------------------------------------------------------- forward / inverse
def state(z, cooling, SG=MUD_SG):
    """Stresses, pressures and temperatures at TVD z with the wall cooled.

    The wall is never colder than the surface (the mud's inlet temperature),
    so shallow depths get less than the nominal cooling."""
    T_rock = T_formation(z)
    T_wall = max(T_rock - cooling, min(T_rock, T_SURF))
    return dict(Sh=PROFILE.Shmin(z), SH=PROFILE.SHmax(z), Pp=float(PROFILE.Pp(z)),
                Pw=SG * C.MUD_SG_GRAD * z, T_rock=T_rock, T_wall=T_wall,
                dT=se.thermal_hoop_stress(T_rock, T_wall))


def implied_ucs(z, width_deg, cooling, SG=MUD_SG, phi_deg=PHI_DEG):
    """UCS [Pa, at 25 C] of rock that breaks out to exactly width_deg at z.

    The failure criterion around the hole: failure where cos 2phi > -a/d with
    a = SH + Sh - Pw + dT - Pp - sigma_cm - K (Pw - Pp), d = 2 (SH - Sh).
    A lobe of width W has -a/d = cos W, so
    sigma_cm = SH + Sh - Pw + dT - Pp - K (Pw - Pp) + d cos W.
    width 0 gives the strength that just holds the wall intact.
    """
    s = state(z, cooling, SG)
    K = kmc(phi_deg)
    scm = (s["SH"] + s["Sh"] - s["Pw"] + s["dT"] - s["Pp"] - K * (s["Pw"] - s["Pp"])
           + 2 * (s["SH"] - s["Sh"]) * np.cos(np.radians(width_deg)))
    return scm / m5.sigma_cm(s["T_wall"], 1.0)


def predicted_width(z, ucs, cooling, SG=MUD_SG):
    s = state(z, cooling, SG)
    return m5.breakout_width(s["SH"], s["Sh"], s["Pw"], s["Pp"], s["T_wall"], ucs, s["dT"])


def onset_depth(ucs, cooling, z_lo=500.0, z_hi=6000.0, SG=MUD_SG):
    """Shallowest TVD at which rock of this UCS breaks out (bisection)."""
    f = lambda z: implied_ucs(z, 0.0, cooling, SG) - ucs
    if f(z_hi) <= 0:
        return None
    if f(z_lo) > 0:
        return z_lo
    for _ in range(60):
        mid = 0.5 * (z_lo + z_hi)
        z_lo, z_hi = (mid, z_hi) if f(mid) <= 0 else (z_lo, mid)
    return z_hi


# ------------------------------------------------------------- calibration
def section_breakouts(features, tvd):
    lo, hi = SECTION_MD
    return [dict(f, tvd=tvd(f["md"])) for f in features
            if f["cls"] == "Breakout" and lo <= f["md"] < hi]


def intact_depths(breakouts, tvd, step=5.0):
    """TVDs in the 12.25" section outside every logged breakout."""
    lo, hi = SECTION_MD
    spans = [(tvd(b["md"] - b["height"] / 2), tvd(b["md"] + b["height"] / 2))
             for b in breakouts]
    zs = np.arange(tvd(lo), tvd(hi) + 1e-9, step)
    return [z for z in zs if not any(a <= z <= b for a, b in spans)]


def calibrate(features=None, SG=MUD_SG, phi_deg=PHI_DEG, cooling=COOLING_K):
    """Per cooling level: weak-zone UCS per breakout and the intact lower bound."""
    features = load_log() if features is None else features
    tvd = md_to_tvd(features)
    bo = section_breakouts(features, tvd)
    intact_z = intact_depths(bo, tvd)
    out = dict(breakouts=bo, tvd=tvd, SG=SG, phi_deg=phi_deg, levels={})
    for dK in cooling:
        weak = np.array([implied_ucs(b["tvd"], b["width"], dK, SG, phi_deg) for b in bo])
        need = np.array([implied_ucs(z, 0.0, dK, SG, phi_deg) for z in intact_z])
        out["levels"][dK] = dict(
            weak=weak, weak_p10=np.percentile(weak, 10), weak_p50=np.median(weak),
            weak_p90=np.percentile(weak, 90),
            intact_min=need.max(), intact_binding_tvd=intact_z[int(need.argmax())])
    return out


def strength_cases(cal):
    """The United Downs strength cases for the site runs: weak-zone median
    and intact lower bound, each at both ends of the cooling bracket, so the
    site verdict carries the calibration's full spread. Values in Pa."""
    lv = cal["levels"]
    lo, hi = min(lv), max(lv)
    return (("intact, 0 K cooling", lv[lo]["intact_min"]),
            ("intact, 40 K cooling", lv[hi]["intact_min"]),
            ("weak zones, 0 K cooling", lv[lo]["weak_p50"]),
            ("weak zones, 40 K cooling", lv[hi]["weak_p50"]))


# ------------------------------------------------------------------ checks
def td_window(ucs, cooling, SG=MUD_SG):
    """Mud window at UD-1's TD, as if vertical (the hole is ~35 deg there)."""
    z = TD_TVD
    s = state(z, cooling, SG)
    return se.mud_window(s["Sh"], s["SH"], s["Pp"], z, s["Pw"], s["T_rock"],
                         T_wall=s["T_wall"], UCS=ucs, dsigma_T=s["dT"])


def tensile_threshold_SG(z, cooling, T0=C.WALL_T0):
    """Mud SG above which tensile fractures initiate at the SHmax azimuth."""
    s = state(z, cooling)
    return (3 * s["Sh"] - s["SH"] + s["dT"] - s["Pp"] + T0) / (C.MUD_SG_GRAD * z)


def depth_table(cal, features, cooling, bins=np.append(np.arange(900.0, 4000.0, 250.0), 4000.0)):
    """Observed against predicted breakouts by 250 m MD bin."""
    lv = cal["levels"][cooling]
    tvd = cal["tvd"]
    rows = []
    for a, b in zip(bins[:-1], bins[1:]):
        obs = [x for x in cal["breakouts"] if a <= x["md"] < b]
        z = tvd(0.5 * (a + b))
        rows.append(dict(
            md=(a, b), n=len(obs), length=sum(x["height"] for x in obs),
            width_obs=float(np.median([x["width"] for x in obs])) if obs else None,
            w_weak=predicted_width(z, lv["weak_p50"], cooling),
            w_intact=predicted_width(z, lv["intact_min"], cooling)))
    return rows


# ------------------------------------------------------------------ report
def report(cal, features):
    MPa = 1e6
    bo = cal["breakouts"]
    print("=" * 88)
    print("UD-1 STRENGTH CALIBRATION  (Reinecker et al. 2021 stresses + BGS image log)")
    print("=" * 88)
    print(f"12.25\" section {SECTION_MD[0]:.0f}-{SECTION_MD[1]:.0f} m MD: {len(bo)} breakouts, "
          f"{sum(b['height'] for b in bo):.0f} m, widths {min(b['width'] for b in bo):.0f}-"
          f"{max(b['width'] for b in bo):.0f} deg (median {np.median([b['width'] for b in bo]):.0f})")
    tvd = cal["tvd"]
    last = features[-1]
    td = tvd(last["md"]) + (5275.0 - last["md"]) * np.cos(np.radians(last["dev"]))
    print(f"MD->TVD: average-angle over {len(features)} stations; "
          f"4,000 m MD = {tvd(4000.0):.0f} m TVD; extended to TD, {td:.0f} m TVD "
          f"(paper: {TD_TVD:.0f})")
    print(f"mud {cal['SG']:.2f} SG, friction angle {cal['phi_deg']:.0f} deg, "
          f"formation {T_SURF:.0f} C + {(T_5KM - T_SURF) / 5:.1f} C/km, thermal hoop stress ON")
    print("-" * 88)
    print(f"{'wall cooling':>12} | {'weak-zone UCS P10 / P50 / P90 (MPa)':^38} | "
          f"{'intact UCS >= (MPa)':^20} | binding at")
    for dK, lv in cal["levels"].items():
        print(f"{dK:>10.0f} K | {lv['weak_p10']/MPa:>11.0f} {lv['weak_p50']/MPa:>8.0f} "
              f"{lv['weak_p90']/MPa:>8.0f}{'':>9} | {lv['intact_min']/MPa:>14.0f}{'':>6} | "
              f"{lv['intact_binding_tvd']:.0f} m TVD")
    lv = cal["levels"]
    lo, hi = min(lv), max(lv)
    print("-" * 88)
    print(f"STRENGTH RANGE (cooling 0-40 K): weak zones {lv[hi]['weak_p50']/MPa:.0f}-"
          f"{lv[lo]['weak_p50']/MPa:.0f} MPa (median), "
          f"{lv[hi]['weak_p10']/MPa:.0f}-{lv[lo]['weak_p90']/MPa:.0f} MPa (P10-P90);")
    print(f"                                 intact rock at least "
          f"{lv[hi]['intact_min']/MPa:.0f}-{lv[lo]['intact_min']/MPa:.0f} MPa")
    for dK in (lo, hi):
        for name, u in (("weak P50", lv[dK]["weak_p50"]), ("intact", lv[dK]["intact_min"])):
            z0 = onset_depth(u, dK)
            print(f"  {dK:.0f} K, {name:<8}: breakouts start at "
                  f"{'none above 6 km' if z0 is None else f'{z0:.0f} m TVD'}")
    print(f"  logged: first breakout {min(b['md'] for b in bo):.0f} m MD")

    print("-" * 88)
    print("PREDICTED vs OBSERVED by depth (20 K cooling, mid-bracket)")
    print(f"{'MD bin':>12} {'n':>4} {'length':>8} {'obs width':>10} | "
          f"{'pred width, weak P50':>21} {'pred width, intact':>19}")
    for r in depth_table(cal, features, 20.0):
        wo = "-" if r["width_obs"] is None else f"{r['width_obs']:.0f}"
        print(f"{r['md'][0]:>5.0f}-{r['md'][1]:<6.0f} {r['n']:>4} {r['length']:>7.1f}m {wo:>10} | "
              f"{r['w_weak']:>19.0f} d {r['w_intact']:>17.0f} d")

    print("-" * 88)
    print("CHECKS")
    w12 = max(b["width"] for b in bo)
    deep = [f for f in features if f["cls"] == "Breakout" and f["md"] >= SECTION_MD[1]]
    print(f"  W_max: widest breakout in the 12.25\" section (drilled without trouble) "
          f"{w12:.0f} deg; in the deviated 8.5\" section {max(f['width'] for f in deep):.0f} deg")
    print(f"      -> site-calibrated W_max = {w12:.0f} deg, against the "
          f"{C.BREAKOUT_W_MAX_DEG:.0f} deg default (Zoback 2007)")
    for dK in (lo, hi):
        for name, u in (("weak P50", lv[dK]["weak_p50"]), ("intact", lv[dK]["intact_min"])):
            w = td_window(u, dK)
            print(f"  TD {TD_TVD:.0f} m TVD, {dK:.0f} K, {name:<8}: {w['verdict']:<11} "
                  f"breakout {w['width_hydro']:.0f} deg at {MUD_SG} SG mud; "
                  f"window {w['SG_lo']:.2f}-{w['SG_hi']:.2f} SG")
    for dK in (lo, hi):
        sg = [tensile_threshold_SG(z, dK) for z in (1000.0, 2500.0, 3700.0)]
        print(f"  tensile initiation threshold, {dK:.0f} K: {sg[0]:.2f} / {sg[1]:.2f} / "
              f"{sg[2]:.2f} SG at 1.0 / 2.5 / 3.7 km -> predicted at {MUD_SG} SG "
              f"{'everywhere' if max(sg) < MUD_SG else 'in part'}; logged at 2,665 and "
              f"3,667 m MD (BGS), c. 2,500-3,700 m (paper)")
    print("-" * 88)
    print("SENSITIVITY (20 K cooling): weak-zone P50 / intact bound, MPa")
    for label, kw in (("mud 1.00 SG", dict(SG=1.00)), ("friction angle 30 deg", dict(phi_deg=30.0)),
                      ("friction angle 40 deg", dict(phi_deg=40.0))):
        c = calibrate(features, cooling=(20.0,), **kw)["levels"][20.0]
        print(f"  {label:<22} {c['weak_p50']/MPa:6.0f} / {c['intact_min']/MPa:6.0f}")
    print("=" * 88)


# ------------------------------------------------------------------ figure
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"     # reference slots 1-3
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def figure(cal, features, path=FIG):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    MPa = 1e6
    lv = cal["levels"]
    lo, hi = min(lv), max(lv)
    tvd = cal["tvd"]
    bo = cal["breakouts"]
    deep = [dict(f, tvd=tvd(f["md"])) for f in features
            if f["cls"] == "Breakout" and f["md"] >= SECTION_MD[1]]
    zs = np.linspace(900.0, TD_TVD, 300)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 6.5), sharey=True,
                                 facecolor=SURF)
    for a in (a1, a2):
        a.set_facecolor(SURF)
        a.grid(True, color=GRID, lw=0.8)
        a.tick_params(colors=INK2)
        for s in a.spines.values():
            s.set_color(GRID)
        a.axhline(tvd(SECTION_MD[1]), color=INK2, lw=0.8, ls=":")

    # left: breakout width with depth
    for dK, alpha in ((lo, 1.0), (hi, 0.55)):
        a1.plot([predicted_width(z, lv[dK]["weak_p50"], dK) for z in zs], zs,
                color=ORANGE, lw=2, alpha=alpha)
        a1.plot([predicted_width(z, lv[dK]["intact_min"], dK) for z in zs], zs,
                color=AQUA, lw=2, alpha=alpha)
    a1.scatter([b["width"] for b in bo], [b["tvd"] for b in bo], s=36, color=BLUE,
               edgecolor=SURF, linewidth=1.5, zorder=3)
    a1.scatter([b["width"] for b in deep], [b["tvd"] for b in deep], s=36,
               facecolor="none", edgecolor=BLUE, linewidth=1.5, zorder=3)
    for W, lab, y in ((C.BREAKOUT_W_MAX_DEG, "W_max 90", 950),
                      (63.0, "widest in 12.25\"\n(63, site-calibrated)", 1250)):
        a1.axvline(W, color=INK2, lw=0.8, ls="--")
        a1.text(W - 1, y, lab, color=INK2, fontsize=8, va="top", ha="right")
    a1.text(98, tvd(SECTION_MD[1]) - 30, "12.25\" above | 8.5\" deviated below",
            color=INK2, fontsize=8, ha="right")
    a1.text(76, 4980, "weak-zone\nmedian", color=INK, fontsize=8, va="center")
    a1.text(4, 4650, "intact\nbound", color=INK, fontsize=8, va="center")
    a1.text(4, 1900, "logged breakouts:\n12.25\" filled, 8.5\" hollow\n"
            "curves: 0 K and 40 K cooling\n(nearly identical)", color=INK, fontsize=8)
    a1.set_xlim(0, 100)
    a1.set_xlabel("breakout width (deg)", color=INK)
    a1.set_ylabel("TVD (m)", color=INK)
    a1.set_title("Breakout width: logged and predicted at 1.05 SG", color=INK,
                 fontsize=10, loc="left")

    # right: implied strength per breakout across the cooling bracket
    for b in bo:
        u = [implied_ucs(b["tvd"], b["width"], dK) / MPa for dK in (lo, hi)]
        a2.plot(u, [b["tvd"]] * 2, color=BLUE, lw=2, solid_capstyle="round")
    a2.axvspan(lv[hi]["intact_min"] / MPa, lv[lo]["intact_min"] / MPa, color=AQUA, alpha=0.18)
    a2.axvspan(lv[hi]["weak_p50"] / MPa, lv[lo]["weak_p50"] / MPa, color=ORANGE, alpha=0.18)
    a2.axvline(CORNWALL.UCS / MPa, color=INK2, lw=0.8, ls="--")
    a2.text(CORNWALL.UCS / MPa + 2, 950, "v1.0 UCS 180", color=INK2, fontsize=8, va="top")
    a2.text(lv[hi]["intact_min"] / MPa + 1, 1250, "intact\nlower bound", color=INK, fontsize=8)
    a2.text(lv[hi]["weak_p50"] / MPa + 1, 1250, "weak-zone\nmedian", color=INK, fontsize=8)
    a2.text(35, 2050, "each line: one breakout,\n40 K (left) to 0 K (right)",
            color=INK, fontsize=8)
    a2.set_xlim(30, 230)
    a2.set_xlabel("implied UCS at 25 C (MPa)", color=INK)
    a2.set_title("Strength implied by each breakout", color=INK, fontsize=10, loc="left")
    a1.set_ylim(TD_TVD + 100, 800)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=SURF)
    plt.close(fig)
    return path


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    feats = load_log()
    cal = calibrate(feats)
    report(cal, feats)
    print("wrote", figure(cal, feats))
