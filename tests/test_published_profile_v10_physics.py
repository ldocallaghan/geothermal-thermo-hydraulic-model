"""The published United Downs stress profile under v1.0 physics (spec test 2).

This is the comparison that motivates C1: swap only the stress inputs, keep the
v1.0 breakout arithmetic, and United Downs stops being frac-limited.

SOURCE: Reinecker et al. (2021), Geothermics 97, 102226, section 6.3. The three
gradients below were checked against the paper and are correct as the spec gives
them. Two things around them were NOT: the paper reports a pore pressure gradient
(9.494 MPa/km below a fluid level at ~61 m), and it derives SHmax by ASSUMING a
friction coefficient of 0.8. See data/ud1/reinecker2021.md.
"""
import numpy as np
import pytest

import geo_constants as C
import site_evaluation as se

MPa = 1.0e6

# Reinecker et al. (2021), section 6.3 -- verified against the paper.
SV_GRAD = 25.275 * MPa / 1000.0                  # Pa/m
SHMIN_GRAD, SHMIN_INTERCEPT = 13.21 * MPa / 1000.0, 3.0 * MPa
SHMAX_GRAD, SHMAX_INTERCEPT = 25.99 * MPa / 1000.0, 5.9 * MPa

# Pore pressure: 9.494 MPa/km below a static fluid level at ~61 m below ground.
PP_GRAD, PP_FLUID_LEVEL = 9.494 * MPa / 1000.0, 61.0
MU_PAPER = 0.8           # the friction coefficient the paper's SHmax is built on

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
    # 180 C wall: the paper gives "around 180 C at 5 km". At the 200 C cap
    # instead, P_need reads 37.0 MPa.
    b = se.breakout_v10(Shmin, SHmax, P_fluid, 180.0, T_wall=180.0)

    assert Shmin / MPa == pytest.approx(69.8, abs=0.5)
    assert b["P_need"] / MPa == pytest.approx(36.3, abs=0.5)
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
    """Effective S1/S3 at 10 MPa/km hydrostatic Pp: 4.5 to 4.8, i.e. just under
    the 4.68 cap at mu 0.85 at 5 km and just over it at 12.5 km. Under the
    paper's own pore pressure it is 4.0 to 4.25, under the 4.33 cap at the
    default mu 0.8 (D1) -- see the tests further down."""
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
    them): 3.12 at mu 0.6, 4.33 at the default 0.8 (D1), 4.68 at 0.85."""
    def cap(mu):
        return (np.sqrt(1 + mu ** 2) + mu) ** 2
    assert cap(0.6) == pytest.approx(3.12, abs=0.01)
    assert cap(0.8) == pytest.approx(4.33, abs=0.01)
    assert cap(0.85) == pytest.approx(4.68, abs=0.01)


# --------------------------------------- the paper's own pore pressure (6.3)
def Pp_paper(z):
    return PP_GRAD * max(z - PP_FLUID_LEVEL, 0.0)


def cap(mu):
    return (np.sqrt(1 + mu ** 2) + mu) ** 2


def implied_mu(ratio):
    """Friction coefficient a given effective S1/S3 sits at equilibrium with."""
    return (ratio - 1) / (2 * np.sqrt(ratio))


def test_paper_pore_pressure_is_not_the_models_hydrostatic():
    """9.494 MPa/km below 61 m, against C.HYDROSTATIC_GRAD's 10 MPa/km from
    surface: 3 MPa lighter at 5 km, 7 MPa at 12.5 km. C3 needs the per-site
    override for exactly this."""
    assert Pp_paper(5000.0) / MPa == pytest.approx(46.9, abs=0.1)
    assert C.HYDROSTATIC_GRAD * 5000.0 / MPa == pytest.approx(50.0, abs=0.1)
    assert Pp_paper(12500.0) / MPa == pytest.approx(118.1, abs=0.1)


@pytest.mark.parametrize("z, ratio, mu", [
    (5000.0, 4.01, 0.75),
    (12500.0, 4.25, 0.79),
])
def test_published_profile_sits_at_mu_0p8_under_its_own_pore_pressure(z, ratio, mu):
    """The profile implies the 0.8 it was CONSTRUCTED with, not 0.83-0.86.
    The spec's higher figures come from substituting 10 MPa/km hydrostatic."""
    _, Shmin, SHmax = published(z)
    R = (SHmax - Pp_paper(z)) / (Shmin - Pp_paper(z))
    assert R == pytest.approx(ratio, abs=0.01)
    assert implied_mu(R) == pytest.approx(mu, abs=0.01)


@pytest.mark.parametrize("mu", [0.8, 0.85])
def test_published_profile_is_admissible_at_12p5km_under_the_papers_pore_pressure(mu):
    """CONTRADICTS the C2 done-when ("the cap binds on it at 12.5 km"). That
    holds only with 10 MPa/km hydrostatic Pp and mu = 0.85 (4.77 against 4.68).
    With the paper's pore pressure the extrapolated profile is 4.25, inside the
    cap at both 0.8 (4.33) and 0.85 (4.68), so nothing binds."""
    _, Shmin, SHmax = published(12500.0)
    Pp = Pp_paper(12500.0)
    assert SHmax <= Pp + cap(mu) * (Shmin - Pp)


def test_frictional_cap_at_the_papers_friction_coefficient():
    assert cap(MU_PAPER) == pytest.approx(4.33, abs=0.01)
