"""Regression anchors for the v1.0 stability arithmetic.

These hold v1.0's behaviour fixed so that v1.1 can be compared with it. They
run against the real code (`site_evaluation.stresses_v10`,
`site_evaluation.breakout_v10`, `model5`) and deliberately avoid Models 1-4,
so they need neither the IAPWS water table nor a BVP solve.

The targets are v1.0's numbers at commit d927bc4, checked against UD-1 where
the well gives a comparison. Tolerances are tight because the arithmetic is
deterministic; they are stated in MPa so a failure reads in physical units.
"""
import numpy as np
import pytest

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, SITES

MPa = 1.0e6


# --------------------------------------------------------------- primitives
def test_mohr_coulomb_coefficient_is_35_degree_friction_angle():
    phi = np.radians(35.0)
    assert m5.KMC == pytest.approx((1 + np.sin(phi)) / (1 - np.sin(phi)))
    assert m5.KMC == pytest.approx(3.6902, abs=1e-4)


@pytest.mark.parametrize("T_C, expected_MPa", [
    (25.0, 200.0),      # no derating at the reference temperature
    (180.0, 172.1),     # UD-1 wall at 5,058 m
    (200.0, 168.5),     # the v1.0 cooled-wall assumption
    (400.0, 132.5),     # rock temperature at the 400 C targets
])
def test_sigma_cm_temperature_derating(T_C, expected_MPa):
    """With no UCS passed, sigma_cm reads the GLOBAL C.UCS = 200 MPa: the v1.0
    behaviour, kept as the default so the v1.0 adapter reproduces it."""
    assert C.UCS == 200.0 * MPa
    assert m5.sigma_cm(T_C) / MPa == pytest.approx(expected_MPa, abs=0.1)


# ------------------------------------------------- Cornwall v1.0 stress state
def test_cornwall_v10_anisotropy_is_constant_with_depth():
    """v1.0 builds stresses from constant ratios, so SHmax/Shmin is 2.55 at
    every depth. v1.1 replaces it with the measured profile (~1.97)."""
    assert CORNWALL.anisotropy == pytest.approx(2.545, abs=0.005)
    for z in (900.0, 4000.0, 5058.0, 12500.0):
        _, Shmin, SHmax, _ = se.stresses_v10(CORNWALL, z)
        assert SHmax / Shmin == pytest.approx(2.545, abs=0.005)


def test_cornwall_v10_effective_stress_ratio_exceeds_every_frictional_cap():
    """Effective S1/S3 = 6.23 at hydrostatic Pp, against frictional caps of
    3.12 (mu 0.6) and 4.68 (mu 0.85): v1.0's stress input was frictionally
    impossible. v1.0 had no such check; v1.1's admissibility check flags it."""
    z = 5058.0
    _, Shmin, SHmax, _ = se.stresses_v10(CORNWALL, z)
    Pp = C.HYDROSTATIC_GRAD * z
    assert (SHmax - Pp) / (Shmin - Pp) == pytest.approx(6.23, abs=0.01)


# -------------------------------------------------------- UD-1 at 5,058 m TVD
def test_ud1_measured_depth_mud_to_suppress_breakout_vs_shmin():
    """v1.0 inputs at UD-1's TD: 65 MPa needed against Shmin of 72 MPa.
    The published profile gives 36 vs 70 (test_published_profile_v10_physics.py)."""
    z = 5058.0
    _, Shmin, SHmax, P_fluid = se.stresses_v10(CORNWALL, z)
    T_rock = float(CORNWALL.geotherm(z))
    # Wall at 180 C, Reinecker et al. (2021)'s "around 180 C at 5 km". The repo
    # geotherm reads 191.6 C here; at that temperature P_need is 65.3 MPa, so
    # the choice moves nothing at the precision these anchors use.
    b = se.breakout_v10(Shmin, SHmax, P_fluid, T_rock, T_wall=180.0)

    assert T_rock == pytest.approx(191.6, abs=0.1)
    assert Shmin / MPa == pytest.approx(71.8, abs=0.5)
    assert b["P_need"] / MPa == pytest.approx(64.9, abs=0.5)
    assert bool(b["breaks"]) is True             # breaks out at hydrostatic mud
    assert bool(b["frac_limited"]) is False      # but 65 < 72: still mud-controllable


# ------------------------------------------------------- United Downs target
def test_united_downs_target_depth_is_frac_limited():
    """THE v1.0 verdict driver: at 12.5 km the mud needed to suppress breakout
    (215 MPa) exceeds Shmin (177 MPa), so v1.0 calls the site NO-GO."""
    z = CORNWALL.target_depth
    assert z == 12500.0
    _, Shmin, SHmax, P_fluid = se.stresses_v10(CORNWALL, z)
    T_rock = float(CORNWALL.geotherm(z))
    b = se.breakout_v10(Shmin, SHmax, P_fluid, T_rock)

    assert Shmin / MPa == pytest.approx(177.4, abs=0.5)
    assert b["P_need"] / MPa == pytest.approx(215.0, abs=0.5)
    assert bool(b["frac_limited"]) is True
    assert bool(b["breaks"]) is True


# ------------------------------- breakouts over the logged UD-1 interval
def _breakout_depths(site, mud_SG, depths):
    """Depths [m] where the v1.0 check predicts breakout at this mud weight."""
    out = []
    for z in depths:
        _, Shmin, SHmax, _ = se.stresses_v10(site, z)
        Pw = mud_SG * 9.81e3 * z
        T_rock = float(site.geotherm(z))
        if bool(se.breakout_v10(Shmin, SHmax, Pw, T_rock)["breaks"]):
            out.append(float(z))
    return out


LOGGED_INTERVAL = np.arange(900.0, 4001.0, 100.0)


def test_v10_predicts_breakout_only_at_the_bottom_of_the_logged_interval():
    """UD-1 broke out repeatedly between ~900 and ~4,000 m MD (27 breakouts
    totalling 139 m per Reinecker et al. 2021; 24 totalling 46 m in the BGS log).

    v1.0 inputs reproduce almost none of that: at the 12.25" section's mud
    weight (1.05 SG) breakout appears only at 3,900-4,000 m, and at 1.10 SG,
    the shallower sections' mud weight, nowhere at all.
    """
    assert _breakout_depths(CORNWALL, 1.05, LOGGED_INTERVAL) == [3900.0, 4000.0]
    assert _breakout_depths(CORNWALL, 1.10, LOGGED_INTERVAL) == []


# ------------------------------------------------------------- the four sites
@pytest.mark.parametrize("site, frac_limited", [
    ("Upper Rhine Graben / Soultz-sous-Forets (France)", False),
    ("Larderello (Italy)", False),
    ("United Downs / Carnmenellis (UK)", True),
    ("Pannonian Basin (Hungary)", False),
])
def test_v10_frac_limited_flag_per_site(site, frac_limited):
    """The stability half of the v1.0 comparative table, without Models 1-4.
    United Downs is the only site v1.0 rules out, and frac_limited is why."""
    s = next(x for x in SITES if x.name == site)
    _, Shmin, SHmax, P_fluid = se.stresses_v10(s, s.target_depth)
    T_rock = float(s.geotherm(s.target_depth))
    b = se.breakout_v10(Shmin, SHmax, P_fluid, T_rock)
    assert bool(b["frac_limited"]) is frac_limited


def test_per_site_ucs_is_ignored_by_the_v10_default():
    """Each SiteProfile carries a UCS; v1.0 never used it. v1.1 site runs pass
    their own strengths explicitly (test_site_strength.py); the default stays
    v1.0's global value."""
    assert CORNWALL.UCS == 180.0 * MPa
    assert m5.sigma_cm(200.0) == pytest.approx(C.UCS * (1 - 0.0009 * 175))
