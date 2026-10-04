"""Print the v1.0 regression anchors as a table, beside the published profile.

    python tests/report_v10_anchors.py

Stability arithmetic only -- no Models 1-4, so no water table needed. Every row
is covered by an assertion in tests/test_v10_stability_anchors.py or
tests/test_published_profile_v10_physics.py.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import geo_constants as C                                   # noqa: E402
import model5_convergence_confinement as m5                 # noqa: E402
import site_evaluation as se                                # noqa: E402
from comparative_sites import CORNWALL                      # noqa: E402

MPa = 1.0e6
SV, SHMIN, SHMAX = 25.275, (13.21, 3.0), (25.99, 5.9)       # published, MPa/km + MPa


def pub(z):
    return (SV * z / 1000.0 * MPa,
            (SHMIN[0] * z / 1000.0 + SHMIN[1]) * MPa,
            (SHMAX[0] * z / 1000.0 + SHMAX[1]) * MPa)


def v10(z):
    return se.stresses_v10(CORNWALL, z)[:3]


def cap(mu):
    return (np.sqrt(1 + mu ** 2) + mu) ** 2


def eff_ratio(Shmin, SHmax, z):
    Pp = C.HYDROSTATIC_GRAD * z
    return (SHmax - Pp) / (Shmin - Pp)


def breakout_depths(stresses, mud_SG, T_of_z):
    hits = []
    for z in np.arange(900.0, 4001.0, 100.0):
        _, Shmin, SHmax = stresses(z)
        b = se.breakout_v10(Shmin, SHmax, mud_SG * 9.81e3 * z, T_of_z(z))
        if bool(b["breaks"]):
            hits.append(z)
    return hits


def main():
    w = 100
    print("=" * w)
    print("v1.0 REGRESSION ANCHORS  (commit d927bc4 stability arithmetic)")
    print("=" * w)
    print(f"{'check':<46}{'v1.0 inputs':>26}{'published UDDGP':>28}")
    print("-" * w)

    # anisotropy
    _, sh_v, sH_v = v10(5058.0)
    _, sh_p, sH_p = pub(5058.0)
    print(f"{'SHmax/Shmin at United Downs':<46}"
          f"{f'{sH_v / sh_v:.2f} (any z)':>26}{f'{sH_p / sh_p:.2f} at 5 km':>28}")
    _, sh_p12, sH_p12 = pub(12500.0)
    print(f"{'   same, at 12.5 km':<46}{sH_v / sh_v:>26.2f}{sH_p12 / sh_p12:>28.2f}")

    # effective stress ratio
    print(f"{'effective S1/S3, hydrostatic Pp (5 km)':<46}"
          f"{eff_ratio(*v10(5000.0)[1:], 5000.0):>26.2f}"
          f"{eff_ratio(*pub(5000.0)[1:], 5000.0):>28.2f}")
    print(f"{'   same, at 12.5 km':<46}"
          f"{eff_ratio(*v10(12500.0)[1:], 12500.0):>26.2f}"
          f"{eff_ratio(*pub(12500.0)[1:], 12500.0):>28.2f}")
    caps = f"{cap(0.6):.2f} / {cap(0.8):.2f} / {cap(0.85):.2f}"
    print(f"{'frictional cap R(mu) at 0.6 / 0.8 / 0.85':<46}{caps:>26}{caps:>28}")
    print("-" * w)

    # mud to suppress breakout
    for z, label, T_wall in ((5058.0, "UD-1 at 5,058 m TVD", 180.0),
                             (12500.0, "target depth 12.5 km", None)):
        cells = []
        for stresses in (v10, z and pub):
            _, Shmin, SHmax = stresses(z)
            T_rock = float(CORNWALL.geotherm(z))
            b = se.breakout_v10(Shmin, SHmax, C.HYDROSTATIC_GRAD * z, T_rock,
                                T_wall=T_wall)
            tag = "NO-GO" if bool(b["frac_limited"]) else "CONDITIONAL"
            cells.append(f"{b['P_need'] / MPa:.0f} vs {Shmin / MPa:.0f} MPa, {tag}")
        print(f"{label + ': mud vs Shmin':<46}{cells[0]:>26}{cells[1]:>28}")
    print("-" * w)

    # logged interval
    for mud in (1.05, 1.10):
        hv = breakout_depths(v10, mud, lambda z: float(CORNWALL.geotherm(z)))
        hp = breakout_depths(pub, mud, lambda z: float(CORNWALL.geotherm(z)))
        fmt = lambda h: ("none" if not h else
                         f"{h[0] / 1000:.1f}-{h[-1] / 1000:.1f} km ({len(h)} of 32)")
        print(f"{f'breakout in 900-4,000 m at {mud:.2f} SG':<46}{fmt(hv):>26}{fmt(hp):>28}")
    print(f"{'breakouts logged in UD-1 over that interval':<46}"
          f"{'27, totalling 139 m':>26}{'27, totalling 139 m':>28}")
    print("-" * w)
    print(f"{'rock-mass strength used (global C.UCS)':<46}"
          f"{f'{m5.sigma_cm(200.0) / MPa:.0f} MPa at 200 C':>26}"
          f"{f'site UCS {CORNWALL.UCS / MPa:.0f} MPa ignored':>28}")
    print("=" * w)
    print("Published gradients verified against Reinecker et al. (2021), section 6.3;")
    print("see data/ud1/reinecker2021.md. UD-1 wall at the paper's ~180 C at 5 km.")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    main()
