"""The published United Downs stress profile under v1.0 physics (spec test 2).

This is the comparison that motivates C1: swap only the stress inputs, keep the
v1.0 breakout arithmetic, and United Downs stops being frac-limited.

SOURCE WARNING: the three gradients below are the ones written into the v1.1
spec, attributed to Reinecker et al. (2021), Geothermics. They reached the spec
through a summarising tool and are NOT yet checked against the paper (no access
from this environment). C1 must verify them before they become site data; until
then this file locks the spec's target numbers, not the paper's.
"""
import numpy as np
import pytest

import geo_constants as C
import site_evaluation as se

MPa = 1.0e6

# Reinecker et al. (2021) profile as quoted in the v1.1 spec -- UNVERIFIED.
SV_GRAD = 25.275 * MPa / 1000.0                  # Pa/m
SHMIN_GRAD, SHMIN_INTERCEPT = 13.21 * MPa / 1000.0, 3.0 * MPa
SHMAX_GRAD, SHMAX_INTERCEPT = 25.99 * MPa / 1000.0, 5.9 * MPa

UD1_TD_TVD = 5058.0      # m, UD-1 total depth (5,275 m MD)
TARGET_DEPTH = 12500.0   # m, the v1.0 400 C target depth for the site


def published(z):
    """(Sv, Shmin, SHmax) [Pa] from the published profile, linearly extrapolated."""
    return (SV_GRAD * z,
            SHMIN_GRAD * z + SHMIN_INTERCEPT,
            SHMAX_GRAD * z + SHMAX_INTERCEPT)


@pytest.mark.parametrize("z", [5000.0, 12500.0])
def test_published_profile_anisotropy_is_about_1p97(z):
    """C1's done-when: ~1.97 at both 5 and 12.5 km, against v1.0's 2.55."""
    _, Shmin, SHmax = published(z)
    assert SHmax / Shmin == pytest.approx(1.97, abs=0.02)


def test_published_profile_at_ud1_total_depth():
    """36 MPa of mud needed against Shmin of 70 MPa -- comfortably drillable,
    where v1.0 inputs give 65 against 72."""
    Sv, Shmin, SHmax = published(UD1_TD_TVD)
    P_fluid = C.HYDROSTATIC_GRAD * UD1_TD_TVD
    # 185 C wall, the measured UD-1 temperature at TD and the value the spec
    # uses here. At the 200 C cap instead, P_need reads 37.0 MPa.
    b = se.breakout_v10(Shmin, SHmax, P_fluid, 185.0, T_wall=185.0)

    assert Shmin / MPa == pytest.approx(69.8, abs=0.5)
    assert b["P_need"] / MPa == pytest.approx(36.4, abs=0.5)
    assert bool(b["frac_limited"]) is False


def test_published_profile_at_12p5km_is_conditional_not_frac_limited():
    """Spec test 2: P_need ~140 MPa against Shmin ~168 MPa, so the v1.0 NO-GO
    becomes CONDITIONAL on stress inputs alone -- before any physics change."""
    Sv, Shmin, SHmax = published(TARGET_DEPTH)
    P_fluid = C.HYDROSTATIC_GRAD * TARGET_DEPTH
    b = se.breakout_v10(Shmin, SHmax, P_fluid, 400.0)

    assert Shmin / MPa == pytest.approx(168.1, abs=0.5)
    assert b["P_need"] / MPa == pytest.approx(139.8, abs=0.5)
    assert bool(b["frac_limited"]) is False
    assert bool(b["breaks"]) is True             # breaks at hydrostatic -> CONDITIONAL


@pytest.mark.parametrize("z, expected", [(5000.0, 4.51), (12500.0, 4.78)])
def test_published_profile_effective_stress_ratio(z, expected):
    """Effective S1/S3 at hydrostatic Pp: 4.5 to 4.8, i.e. just under the 4.68
    frictional cap at mu 0.85 at 5 km and just over it at 12.5 km (C2)."""
    _, Shmin, SHmax = published(z)
    Pp = C.HYDROSTATIC_GRAD * z
    assert (SHmax - Pp) / (Shmin - Pp) == pytest.approx(expected, abs=0.02)


def test_published_profile_predicts_no_breakout_over_the_logged_interval():
    """The other half of the v1.0 failure: with the published stresses, v1.0
    physics predicts no breakout anywhere in 900-4,000 m, where 27 were logged.
    So the stress inputs alone do not explain the log -- the physics must change
    too (C3 to C8)."""
    for mud_SG in (1.05, 1.10):
        for z in np.arange(900.0, 4001.0, 100.0):
            _, Shmin, SHmax = published(z)
            Pw = mud_SG * 9.81e3 * z
            T_rock = 15.0 + 0.035 * z
            b = se.breakout_v10(Shmin, SHmax, Pw, T_rock)
            assert bool(b["breaks"]) is False, f"{z} m at {mud_SG} SG"


def test_frictional_caps_bracket_the_two_profiles():
    """The caps themselves (pure arithmetic; C2 adds the function that uses
    them): 3.12 at mu 0.6 and 4.68 at mu 0.85."""
    def cap(mu):
        return (np.sqrt(1 + mu ** 2) + mu) ** 2
    assert cap(0.6) == pytest.approx(3.12, abs=0.01)
    assert cap(0.85) == pytest.approx(4.68, abs=0.01)
