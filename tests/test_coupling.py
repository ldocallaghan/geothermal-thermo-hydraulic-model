"""The evaluation's own coupling: the site's stiffness in the thermal stress,
Model 3 on the site's circulation, and the rate of penetration and exposure
iterated together."""
from dataclasses import replace

import pytest

import geo_constants as C
import model3_optimiser as m3
import site_evaluation as se


def test_a_coupled_site_takes_its_own_stiffness():
    s = replace(se.SOULTZ, coupled=True)
    base = se.thermal_hoop_stress(400.0, 200.0)
    assert se.thermal_hoop_stress(400.0, 200.0, se.SOULTZ) == base
    assert se.thermal_hoop_stress(400.0, 200.0, s) == pytest.approx(base * s.E_rock / C.E_ROCK)


def test_the_face_model_takes_the_coolant_it_is_given():
    """A hotter coolant at the face quenches less and drills slower; without
    the quench the coolant does not matter."""
    cold = m3.evaluate_at(130.0, 10700.0, 400.0, 0.6)
    hot = m3.evaluate_at(185.0, 10700.0, 400.0, 0.6)
    assert hot["ROP"] < cold["ROP"]
    assert hot["T_cold"] == 185.0
    assert (m3.evaluate_at(130.0, 10700.0, 400.0, 0.6, quench=False)["ROP"]
            == m3.evaluate_at(185.0, 10700.0, 400.0, 0.6, quench=False)["ROP"])


@pytest.mark.slow
def test_the_rate_of_penetration_is_consistent_with_the_circulation(evaluated):
    o = evaluated(se.SOULTZ)
    assert o["site"].coupled
    last = o["fluid_passes"][-1]
    assert o["fluid_converged"]
    assert o["drill"]["ROP"] * 3600 == pytest.approx(last["rop_m_h"], rel=se.ROP_TOL)
    assert o["drill"]["T_cold"] == pytest.approx(o["circ"]["bhct"])


@pytest.mark.slow
def test_thermal_off_reaches_every_state_of_the_cycle(evaluated):
    o = evaluated(se.SOULTZ)
    c = se.site_cycle(o["site"], o, thermal=False)
    for state, cases in c["cases"].items():
        assert all(k["dsigma_T"] == 0.0 for k in cases), state
