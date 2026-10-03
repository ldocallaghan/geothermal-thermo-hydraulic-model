"""
crosscheck_soultz.py
===================
C9 of spec v1.1: a sanity check of the v1.1 stability model against the 5 km
Soultz wells. Nothing is fitted to Soultz.

Run from the repo root:   python src/crosscheck_soultz.py

Stresses: Valley & Evans (2007), the SOULTZ stress cases (C1), valid 1.5-5 km.
Mud at hydrostatic ("the annulus pressure was maintained at near-hydrostatic
conditions", Valley & Evans). Formation temperature: the Upper Rhine geotherm,
200 C at 5 km. Wall cooling bracketed 0-40 K, as in the UD-1 calibration.

Strengths tried, from three sources:
  - the model's Soultz UCS (170 MPa, v1.0, uncalibrated)
  - UCS of unaltered Soultz granite, 100-130 MPa (ten lab tests, Valley &
    Evans 2007)
  - the UD-1 calibration (C8): weak-zone median and intact bound, 118-203 MPa

Observed in GPK3 and GPK4 (Valley & Evans 2007), the comparison targets:
  - breakouts at 5 km in both wells, drilled to TD
  - GPK4: no high-confidence breakouts above 3,000 m TVD, sparse to 3,670 m,
    dense below; about 10% of the logged length broken out in each well
  - tensile fractures almost continuous to 2,180 m TVD in GPK4, sporadic below
  - wall thermal stress at the time of logging -17 to -31 MPa
"""
import sys

import numpy as np

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL

SOULTZ = se.SOULTZ
Z_CHECK = 5000.0
COOLING_K = (0.0, 40.0)
OBS_ONSET_FIRST, OBS_ONSET_DENSE = 3000.0, 3670.0     # GPK4, m TVD
OBS_DITF_CONTINUOUS_TO = 2180.0                        # GPK4, m TVD

STRENGTHS = (
    ("site UCS (v1.0)", (SOULTZ.UCS,)),
    ("Soultz lab UCS 100-130", (100e6, 115e6, 130e6)),
    ("UD-1 weak zones (C8)", tuple(u for l, u in CORNWALL.strength_cases if "weak" in l)),
    ("UD-1 intact bound (C8)", tuple(u for l, u in CORNWALL.strength_cases if "intact" in l)),
)


def state(p, z, cooling):
    T_rock = float(se.urg_geotherm(z))
    T_wall = T_rock - cooling
    return dict(Sh=p.Shmin(z), SH=p.SHmax(z), Pp=float(p.Pp(z)),
                Pw=C.HYDROSTATIC_GRAD * z, T_rock=T_rock, T_wall=T_wall,
                dT=se.thermal_hoop_stress(T_rock, T_wall))


def width(p, z, ucs, cooling):
    s = state(p, z, cooling)
    return m5.breakout_width(s["SH"], s["Sh"], s["Pw"], s["Pp"], s["T_wall"], ucs, s["dT"])


def window(p, z, ucs, cooling):
    s = state(p, z, cooling)
    return se.mud_window(s["Sh"], s["SH"], s["Pp"], z, s["Pw"], s["T_rock"],
                         T_wall=s["T_wall"], UCS=ucs, dsigma_T=s["dT"])


def onset(p, ucs, cooling, z_lo=1500.0, z_hi=Z_CHECK, step=10.0):
    """Shallowest TVD in the data's range where the wall breaks out."""
    for z in np.arange(z_lo, z_hi + step / 2, step):
        if width(p, z, ucs, cooling) > 0.0:
            return float(z)
    return None


def tensile_depths(p, cooling, z_lo=1500.0, z_hi=Z_CHECK, step=50.0):
    """TVDs at which hydrostatic mud initiates tensile fractures (T0 = 0)."""
    out = []
    for z in np.arange(z_lo, z_hi + step / 2, step):
        s = state(p, z, cooling)
        if s["Pw"] > 3 * s["Sh"] - s["SH"] + s["dT"] - s["Pp"] + C.WALL_T0:
            out.append(float(z))
    return out


def crosscheck():
    rows = []
    for name, ucss in STRENGTHS:
        for ucs in ucss:
            for dK in COOLING_K:
                rows.append(dict(
                    source=name, ucs=ucs, cooling=dK,
                    width=[width(p, Z_CHECK, ucs, dK) for p in SOULTZ.stress_cases],
                    window=[window(p, Z_CHECK, ucs, dK) for p in SOULTZ.stress_cases],
                    onset=[onset(p, ucs, dK) for p in SOULTZ.stress_cases]))
    return rows


def report(rows):
    MPa = 1e6
    labels = [p.label.split(" (")[0].replace("SHmax ", "") for p in SOULTZ.stress_cases]
    print("=" * 100)
    print("C9  SOULTZ CROSS-CHECK at 5 km  (Valley & Evans 2007 stresses; hydrostatic mud; nothing fitted)")
    print("=" * 100)
    print(f"SHmax cases: {', '.join(labels)}; wall cooling {COOLING_K[0]:.0f} and "
          f"{COOLING_K[1]:.0f} K; widths in deg")
    print(f"{'strength source':<24}{'UCS':>5}{'cool':>6} | {'width at 5 km':^17} | "
          f"{'window at 5 km (SG)':^26} | breakout onset, m TVD")
    for r in rows:
        w = "/".join(f"{x:.0f}" for x in r["width"])
        win = " ".join(f"{x['SG_lo']:.2f}-{x['SG_hi']:.2f}" for x in r["window"][1:2])
        verdicts = "/".join(x["verdict"][:4] for x in r["window"])
        on = "/".join("-" if x is None else f"{x:.0f}" for x in r["onset"])
        print(f"{r['source']:<24}{r['ucs']/MPa:>5.0f}{r['cooling']:>5.0f}K | {w:^17} | "
              f"{win:>11} {verdicts:<14} | {on}")
    print("-" * 100)
    print(f"observed (GPK4): first high-confidence breakouts ~{OBS_ONSET_FIRST:.0f} m, "
          f"dense below {OBS_ONSET_DENSE:.0f} m TVD; breakouts at 5 km in GPK3 and GPK4")
    tens = {dK: tensile_depths(SOULTZ.stress_cases[1], dK) for dK in COOLING_K}
    for dK, zs in tens.items():
        span = "none" if not zs else f"{zs[0]:.0f}-{zs[-1]:.0f} m TVD"
        print(f"tensile initiation, mid SHmax, {dK:.0f} K: {span} "
              f"(observed: continuous to {OBS_DITF_CONTINUOUS_TO:.0f} m, sporadic below)")
    ok = all(x["open"] for r in rows for x in r["window"])
    bo = sum(1 for r in rows for w in r["width"] if w > 0)
    n = sum(len(r["width"]) for r in rows)
    print("-" * 100)
    print(f"breakouts predicted at 5 km in {bo} of {n} strength x cooling x SHmax cases; "
          f"mud window open in {'all' if ok else 'not all'} of them")
    print("=" * 100)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    report(crosscheck())
