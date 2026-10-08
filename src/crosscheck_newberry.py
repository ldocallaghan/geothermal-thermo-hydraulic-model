"""
crosscheck_newberry.py
======================
A check of the stability model against the Newberry EGS well NWG 55-29,
where the wall was imaged at up to 277 C. Nothing is fitted to Newberry.

Run from the repo root:   python src/crosscheck_newberry.py

Source: Davatzes, N.C., Hickman, S.H. (2011), "Preliminary Analysis of Stress
in the Newberry EGS Well NWG 55-29", GRC Transactions 35 (data/newberry/).
Logs: GDR submission 271 (AltaRock Energy), extracted by
data/newberry/build_extracts.py.

  * Stresses, as Davatzes and Hickman built them, at their friction
    coefficient (0.55; 0.70 as their alternative), not the 0.8 used for United
    Downs: their Sv at 8,420 ft (8,852.9 psi) carried linearly in TVD; pore
    pressure from the static survey of October 2008 (the well is
    underpressured, Pf/Sv ~ 0.34); Shmin at frictional failure for normal
    faulting; SHmax from their fit to the breakout widths, 0.772 z + 1,112.6 psi
    (0.76336 z + 1,094.1 at mu 0.70), z the TVD in feet.
  * Strength from their UCS-porosity relation, UCS = 13,800 exp(-0.04744 phi)
    psi, on the neutron porosity log, depth by depth.
  * Borehole fluid at formation pressure (Pm = Pf), as in their model.
  * The wall was logged during inject-to-cool, at up to 277 C, but its
    temperature by depth isn't published. So the wall cooling below the static
    temperature is scanned, and the range consistent with the log is reported.

Logged (their Fig. 3): breakouts throughout the volcanics above about 8,610 ft
MD, modal width 35.86 deg; none in the granodiorite below. The image log ends
at 8,860 ft, so the granodiorite it saw is 8,807-8,860 ft (massive), with
intruded transition rock from 8,610 ft. "Throughout" is taken here as more
than half the logged volcanics broken out.

Their model differs from this one in two ways, both tested here: it sets the
thermal stress to zero, and it compares the hoop stress with UCS directly,
with no loss of strength with temperature.
"""
import csv
import os
import sys

import numpy as np

import model5_convergence_confinement as m5
import site_evaluation as se

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "newberry")
FT, PSI = 0.3048, 6894.757

SV_8420_PSI, MD_8420 = 8852.9, 8420.0
SHMAX_FIT = {0.55: (0.772, 1112.6), 0.70: (0.76336, 1094.1)}   # psi/ft, psi
MU = 0.55
SHOE_MD, VOLCANIC_BASE_MD, GRANODIORITE_TOP_MD, LOG_BASE_MD = 6435.0, 8610.0, 8807.0, 8860.0
OBS_MODAL_WIDTH = 35.86
LOG_MAX_T = 277.0
COOLING_K = np.arange(0.0, 101.0, 5.0)


def _read(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        return list(csv.DictReader(fh))


_dev = _read("nwg5529_deviation.csv")
_DEV_MD = np.array([float(r["md_ft"]) for r in _dev])
_DEV_TVD = np.array([float(r["tvd_ft"]) for r in _dev])
_st = _read("nwg5529_static_2008.csv")
ST_MD = np.array([float(r["md_ft"]) for r in _st])
ST_TVD = np.array([float(r["tvd_ft"]) for r in _st])
ST_T = np.array([float(r["temp_C"]) for r in _st])
ST_P = np.array([float(r["press_psia"]) for r in _st])


def tvd_ft(md_ft):
    return np.interp(md_ft, _DEV_MD, _DEV_TVD)


def static_T(md_ft):
    return np.interp(md_ft, ST_MD, ST_T)


def pore_pressure_fit():
    """Linear fit of the static pressure against TVD over the open hole:
    gradient [Pa/m] and the depth [m] where it reaches zero."""
    m = ST_MD >= SHOE_MD
    a, b = np.polyfit(ST_TVD[m] * FT, ST_P[m] * PSI, 1)
    return float(a), float(-b / a)


def profile(mu=MU):
    """Davatzes and Hickman's stress model as a StressProfile in TVD [m]."""
    Sv_grad = SV_8420_PSI * PSI / (float(tvd_ft(MD_8420)) * FT)
    Pp_grad, Pp_datum = pore_pressure_fit()
    q = (np.sqrt(mu ** 2 + 1) + mu) ** 2
    Sh_grad = Sv_grad / q + Pp_grad * (1 - 1 / q)
    Sh_int = -Pp_grad * Pp_datum * (1 - 1 / q)
    g, c = SHMAX_FIT[mu]
    return se.StressProfile(
        label=f"Davatzes & Hickman (2011), mu {mu}", Sv_grad=Sv_grad,
        Shmin_grad=Sh_grad, Shmin_int=Sh_int, SHmax_grad=g * PSI / FT, SHmax_int=c * PSI,
        Pp_grad=Pp_grad, Pp_datum=Pp_datum,
        source="Davatzes and Hickman (2011); static survey Oct 2008 (GDR 271)", mu=mu,
        z_data=(SHOE_MD * FT, float(tvd_ft(LOG_BASE_MD)) * FT))


def ucs_from_porosity(nphi_pct):
    return 13800.0 * np.exp(-0.04744 * np.asarray(nphi_pct, float)) * PSI


_lg = _read("nwg5529_open_hole_logs.csv")
LOG_MD = np.array([float(r["md_ft_gl"]) for r in _lg])
LOG_NPHI = np.array([float(r["nphi_pct"]) for r in _lg])


def zones():
    return {"volcanics": (SHOE_MD, VOLCANIC_BASE_MD),
            "transition": (VOLCANIC_BASE_MD, GRANODIORITE_TOP_MD),
            "granodiorite (logged)": (GRANODIORITE_TOP_MD, LOG_BASE_MD)}


# The cooled zone for the check that credits cooling only as deep as the
# breakout reaches (model5.breakout_with_skin). The duration of the
# inject-to-cool period before the 2008 log isn't in GDR 271 or the paper, so
# it is bracketed from one hour to a day. 8.5" hole (the LAS bit size); the
# film coefficient as for UD-1.
A_HOLE = 8.5 * 0.0254 / 2
CIRC_HOURS = (1.0, 24.0)
_shapes = {}
_depths = {}


def _shape(hours):
    import calibrate_ud1 as cal
    import wall_thermal as wt
    if hours not in _shapes:
        _shapes[hours] = wt.cooling_shape(hours * 3600.0, A_HOLE, 2.14, 2.7e6, cal.H_ANN_UD1)
    return _shapes[hours]


def widths(p, cooling, ours=True, hours=None):
    """Breakout width at each 2-ft log depth in the imaged open hole.
    ours=False is Davatzes and Hickman's criterion: no thermal stress and no
    loss of strength with temperature."""
    m = (LOG_MD >= SHOE_MD) & (LOG_MD <= LOG_BASE_MD)
    md = LOG_MD[m]
    z = tvd_ft(md) * FT
    ucs = ucs_from_porosity(LOG_NPHI[m])
    out = np.empty(len(md))
    for i, (zi, mdi, u) in enumerate(zip(z, md, ucs)):
        Pp = float(p.Pp(zi))
        T_rock = float(static_T(mdi))
        if ours and hours is not None:
            r_f, g = _shape(hours)
            thermo = -se.thermal_hoop_stress(T_rock, T_rock - 1.0)
            key = (p.label, i)
            if key not in _depths:
                _depths[key] = m5.uncooled_depth(p.SHmax(zi), p.Shmin(zi), Pp, Pp, A_HOLE,
                                                 T_rock, u, thermo)
            out[i] = m5.breakout_with_skin(p.SHmax(zi), p.Shmin(zi), Pp, Pp, A_HOLE, r_f,
                                           T_rock - cooling * g, T_rock, u, thermo,
                                           depth_ref=_depths[key], detail=False)["width"]
            continue
        if ours:
            T_wall = T_rock - cooling
            dT = se.thermal_hoop_stress(T_rock, T_wall)
        else:
            T_wall, dT = 25.0, 0.0          # sigma_cm(25 C) = UCS, no thermal term
        out[i] = m5.breakout_width(p.SHmax(zi), p.Shmin(zi), Pp, Pp, T_wall, u, dT)
    return md, out


def summary(md, w):
    s = {}
    for name, (a, b) in zones().items():
        m = (md >= a) & (md < b)
        wz = w[m]
        s[name] = dict(n=int(m.sum()), frac=float(np.mean(wz > 0)),
                       median=float(np.median(wz[wz > 0])) if np.any(wz > 0) else 0.0)
    return s


def volcanics_match(s, tol=10.0):
    v = s["volcanics"]
    return v["frac"] > 0.5 and abs(v["median"] - OBS_MODAL_WIDTH) <= tol


def consistent(s, tol=10.0):
    return volcanics_match(s, tol) and s["granodiorite (logged)"]["frac"] < 0.1


def crosscheck(mu=MU):
    p = profile(mu)
    rows = [dict(cooling=None, label="their criterion", s=summary(*widths(p, 0.0, ours=False)))]
    for dK in COOLING_K:
        rows.append(dict(cooling=float(dK), label=f"{dK:.0f} K", s=summary(*widths(p, dK))))
    ok = [r["cooling"] for r in rows[1:] if consistent(r["s"])]
    vol = [r["cooling"] for r in rows[1:] if volcanics_match(r["s"])]
    skin = {}
    for h in CIRC_HOURS:
        srows = [dict(cooling=float(dK), s=summary(*widths(p, dK, True, h))) for dK in COOLING_K]
        skin[h] = dict(rows=srows, volcanics_cooling=[r["cooling"] for r in srows if volcanics_match(r["s"])],
                       consistent_cooling=[r["cooling"] for r in srows if consistent(r["s"])])
    return dict(mu=mu, profile=p, rows=rows, consistent_cooling=ok, volcanics_cooling=vol, skin=skin)


def report(res):
    p = res["profile"]
    z = float(tvd_ft(MD_8420)) * FT
    print("=" * 96)
    print(f"NEWBERRY NWG 55-29 CROSS-CHECK  (Davatzes & Hickman 2011 stresses at mu {res['mu']}; "
          "nothing fitted)")
    print("=" * 96)
    print(f"at 8,420 ft MD ({z/FT:.0f} ft TVD): Sv {p.Sv(z)/PSI:.0f}, Pp {p.Pp(z)/PSI:.0f}, "
          f"Shmin {p.Shmin(z)/PSI:.0f}, SHmax {p.SHmax(z)/PSI:.0f} psi "
          f"(paper: 8,853 / 3,122 / 4,675 / 7,305)")
    cap = p.cap_depth(res["mu"])
    print(f"frictional admissibility at mu {res['mu']}: "
          + ("admissible over the open hole" if cap is None or cap > p.z_data[1]
             else f"SHmax over the cap below {cap:.0f} m"))
    for name, (a, b) in zones().items():
        m = (LOG_MD >= a) & (LOG_MD < b)
        u = ucs_from_porosity(LOG_NPHI[m]) / 1e6
        print(f"  {name:<22} {a:.0f}-{b:.0f} ft MD: static {static_T(a):.0f}-{static_T(b):.0f} C, "
              f"UCS from porosity {np.percentile(u, 10):.0f}-{np.percentile(u, 90):.0f} MPa "
              f"(median {np.median(u):.0f})")
    print("-" * 96)
    print(f"{'wall cooling':<16}" + "".join(f"{n:>26}" for n in zones()) + "   consistent")
    print(f"{'':<16}" + "".join(f"{'broken / median width':>26}" for _ in zones()))
    for r in res["rows"]:
        cells = "".join(f"{100*v['frac']:>15.0f}% / {v['median']:>4.0f} deg" for v in r["s"].values())
        flag = "" if r["cooling"] is None else ("  yes" if consistent(r["s"]) else "")
        print(f"{r['label']:<16}{cells}{flag}")
    print("-" * 96)
    print(f"logged: breakouts throughout the volcanics, modal width {OBS_MODAL_WIDTH} deg; "
          "none in the granodiorite")
    print(f"consistent here: volcanics > 50% broken with median within 10 deg of "
          f"{OBS_MODAL_WIDTH:.0f}, granodiorite < 10%")
    ok = res["consistent_cooling"]
    print("wall cooling consistent with the log: "
          + ("none in 0-100 K" if not ok else f"{min(ok):.0f}-{max(ok):.0f} K"))
    vol = res["volcanics_cooling"]
    print("  with the volcanics alone: "
          + ("none" if not vol else f"{min(vol):.0f}-{max(vol):.0f} K"))
    g = [r for r in res["rows"] if r["cooling"] is not None and r["cooling"] in vol]
    if g:
        fr = [100 * r["s"]["granodiorite (logged)"]["frac"] for r in g]
        print(f"  but at that cooling the logged granodiorite breaks out over "
              f"{min(fr):.0f}-{max(fr):.0f}% of its length, where none were logged:")
        print("  the porosity relation gives it 60-82 MPa, little more than the volcanics,")
        print("  so the model can't tell the two apart (their criterion can't either).")
    for h, sk in res["skin"].items():
        v_ = sk["volcanics_cooling"]
        print(f"  with the cooling credited as deep as the breakout reaches (circulated {h:.0f} h): "
              f"volcanics alone {'none' if not v_ else f'{min(v_):.0f}-{max(v_):.0f} K'}; "
              f"consistent with the log {'none' if not sk['consistent_cooling'] else sk['consistent_cooling']}")
    print(f"for scale: the log's maximum temperature, {LOG_MAX_T:.0f} C, against "
          f"{static_T(LOG_BASE_MD):.0f} C static at its base, is "
          f"{static_T(LOG_BASE_MD) - LOG_MAX_T:.0f} K of cooling in the fluid there")
    print("=" * 96)


# ------------------------------------------------------------------ the site
# Static temperature against TVD [m], extrapolated below the deepest reading
# at the gradient of the deepest 1,500 ft (about 460 m) of the survey.
_deep = ST_MD >= ST_MD[-1] - 1500.0
GRAD_DEEP = float(np.polyfit(ST_TVD[_deep] * FT, ST_T[_deep], 1)[0])     # C/m
Z_DATA_T = float(ST_TVD[-1] * FT)
T_DATA_END = float(ST_T[-1])


def geotherm(z):
    z = np.asarray(z, float)
    inside = np.interp(z, ST_TVD * FT, ST_T)
    return np.where(z > Z_DATA_T, T_DATA_END + GRAD_DEEP * (z - Z_DATA_T), inside)


def target_depth(T=400.0):
    return Z_DATA_T + (T - T_DATA_END) / GRAD_DEEP


def granodiorite_strengths():
    """UCS from the porosity relation in the granodiorite, 8,807 ft to TD:
    10th percentile, median and 90th percentile."""
    m = LOG_MD >= GRANODIORITE_TOP_MD
    u = ucs_from_porosity(LOG_NPHI[m])
    return tuple((f"log-derived granodiorite {k}", float(np.percentile(u, q)))
                 for k, q in (("p10", 10), ("median", 50), ("p90", 90)))


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    for mu in (0.55, 0.70):
        report(crosscheck(mu))
