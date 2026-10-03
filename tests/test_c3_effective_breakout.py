"""C3: effective-stress breakout check (spec v1.1).

sigma_theta(theta) = SHmax + Shmin - 2 (SHmax - Shmin) cos 2theta - Pw + dsigma_T,
failing when sigma_theta - Pp > sigma_cm(T_wall, UCS) + KMC (Pw - Pp).
Done when: with C7 off and Pp = 0 it reproduces the v1.0 inequality.
"""
import numpy as np
import pytest

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, SITES

MPa = 1.0e6
UD = CORNWALL.stress_cases[0]


# ------------------------------------------------------------- Kirsch hoop
def test_kirsch_extremes():
    """Spec test 3: 3 SHmax - Shmin - Pw at the Shmin azimuth (90 deg from
    SHmax), 3 Shmin - SHmax - Pw along SHmax (0 deg)."""
    SH, Sh, Pw = 150 * MPa, 80 * MPa, 50 * MPa
    assert m5.hoop_stress(90.0, SH, Sh, Pw) == pytest.approx(3 * SH - Sh - Pw)
    assert m5.hoop_stress(0.0, SH, Sh, Pw) == pytest.approx(3 * Sh - SH - Pw)
    th = np.arange(0.0, 360.0, 0.5)
    s = m5.hoop_stress(th, SH, Sh, Pw)
    assert s.max() == pytest.approx(3 * SH - Sh - Pw)
    assert s.min() == pytest.approx(3 * Sh - SH - Pw)
    assert th[np.argmax(s)] in (90.0, 270.0)


def test_kirsch_isotropic_is_2s_minus_pw():
    s = m5.hoop_stress(np.arange(0.0, 180.0, 15.0), 100 * MPa, 100 * MPa, 40 * MPa)
    assert np.allclose(s, 160 * MPa)


def test_thermal_term_adds_directly():
    base = m5.hoop_stress(90.0, 150 * MPa, 80 * MPa, 50 * MPa)
    assert m5.hoop_stress(90.0, 150 * MPa, 80 * MPa, 50 * MPa, -30 * MPa) == \
        pytest.approx(base - 30 * MPa)


# ------------------------------------------------ reduction to v1.0 (done-when)
@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
@pytest.mark.parametrize("z", [1000.0, 4000.0, 12500.0])
def test_pp_zero_reproduces_v10(site, z):
    _, Shmin, SHmax, Pw = se.stresses_v10(site, z)
    T = float(site.geotherm(z))
    for UCS in (None, site.UCS):
        new = se.breakout_eff(Shmin, SHmax, Pw, 0.0, T, UCS=UCS)
        old = se.breakout_v10(Shmin, SHmax, Pw, T, UCS=UCS)
        assert new.keys() == old.keys()
        for k in old:
            assert new[k] == pytest.approx(old[k], rel=1e-12), k


# --------------------------------------------------------- effective stress
def test_failure_is_mohr_coulomb_in_effective_stress():
    Sh, SH, Pw, Pp, T = 80 * MPa, 150 * MPa, 50 * MPa, 40 * MPa, 150.0
    b = se.breakout_eff(Sh, SH, Pw, Pp, T, UCS=180 * MPa)
    lhs = (3 * SH - Sh - Pw) - Pp
    rhs = m5.sigma_cm(T, 180 * MPa) + m5.KMC * (Pw - Pp)
    assert b["breaks"] is bool(lhs > rhs)
    assert b["mc_cold"] - Pp == pytest.approx(rhs)


def test_p_need_is_where_the_check_flips():
    Sh, SH, Pp, T = 120 * MPa, 260 * MPa, 90 * MPa, 200.0
    P = se.breakout_eff(Sh, SH, 100 * MPa, Pp, T)["P_need"]
    assert se.breakout_eff(Sh, SH, P - 1e3, Pp, T)["breaks"] is True
    assert se.breakout_eff(Sh, SH, P + 1e3, Pp, T)["breaks"] is False


def test_pore_pressure_raises_the_mud_needed():
    """d(P_need)/d(Pp) = (KMC - 1)/(1 + KMC) ~ 0.57: strength is lost faster
    than the hoop stress is relieved."""
    Sh, SH, Pw, T = 120 * MPa, 260 * MPa, 100 * MPa, 200.0
    a = se.breakout_eff(Sh, SH, Pw, 0.0, T)["P_need"]
    b = se.breakout_eff(Sh, SH, Pw, 50 * MPa, T)["P_need"]
    assert (b - a) / (50 * MPa) == pytest.approx((m5.KMC - 1) / (m5.KMC + 1))


def test_site_runs_use_the_profile_pore_pressure():
    z = CORNWALL.target_depth
    ref, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, thermal=False)
    assert ref["Pp"] == pytest.approx(float(UD.Pp(z)))
    expect = se.breakout_eff(ref["Shmin"], ref["SHmax"], C.HYDROSTATIC_GRAD * z,
                             ref["Pp"], 400.0, UCS=CORNWALL.UCS)
    assert ref["breakout"] == expect


def test_united_downs_at_12p5_km():
    """Published stresses, paper Pp, site UCS, 200 C wall, no thermal stress:
    the mud needed rises by (KMC-1)/(KMC+1) * Pp over the total-stress check,
    from 143 to 210 MPa, past Shmin (168). C5-C7 revisit this."""
    z = CORNWALL.target_depth
    ref, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, thermal=False)
    b = ref["breakout"]
    assert b["P_need"] / MPa == pytest.approx(210.9, abs=0.5)
    assert b["frac_limited"] is True


def test_ud1_breakouts_now_appear_in_the_lower_logged_interval():
    """Observation, not calibration (that is C8). At the 12.25" section's mud
    weight (1.05 SG) and the site UCS, with the wall at rock temperature, the
    published profile predicts breakouts at 3.6-4.0 km. Reinecker et al. log
    them over 900-4,000 m MD, 'notably in the lower part'. The total-stress
    check predicted none (test_published_profile_v10_physics.py)."""
    hits = []
    for z in np.arange(900.0, 4001.0, 100.0):
        T = float(CORNWALL.geotherm(z))
        b = se.breakout_eff(UD.Shmin(z), UD.SHmax(z), 1.05 * 9.81e3 * z,
                            float(UD.Pp(z)), T, T_wall=T, UCS=CORNWALL.UCS)
        if b["breaks"]:
            hits.append(float(z))
    assert hits == [3600.0, 3700.0, 3800.0, 3900.0, 4000.0]


def test_ud1_at_td_breaks_out_but_is_not_frac_limited():
    """UD-1 reached 5,058 m TVD with breakouts. At a 180 C wall the check
    needs 67 MPa of mud against Shmin of 70: breakouts at the 52 MPa used,
    but still suppressible below Shmin."""
    z = 5058.0
    b = se.breakout_eff(UD.Shmin(z), UD.SHmax(z), 1.05 * 9.81e3 * z,
                        float(UD.Pp(z)), 180.0, T_wall=180.0, UCS=CORNWALL.UCS)
    assert b["P_need"] / MPa == pytest.approx(67.2, abs=0.3)
    assert b["breaks"] is True and b["frac_limited"] is False
