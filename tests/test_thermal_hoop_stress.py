"""Thermal hoop stress at the cooled wall.

dsigma_T = -E alpha / (1 - nu) * (T_rock - T_wall), with E and alpha from
Model 2, the wall at Model 1's circulating bottom-hole temperature. On by
default in the verdicts; it can be switched off, the v1.0 adapter ignores it,
and every site run also carries the verdict without it as a sensitivity.
"""
import numpy as np
import pytest

import geo_constants as C
import model2_spallation as m2
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, SITES

MPa = 1.0e6


def test_formula_uses_model2_properties_at_rock_temperature():
    T_rock, T_wall = 400.0, 200.0
    expect = -m2.E_of_T(T_rock) * m2.alpha_of_T(T_rock) / (1 - C.NU_ROCK) * 200.0
    assert se.thermal_hoop_stress(T_rock, T_wall) == pytest.approx(expect)


def test_cooling_is_tensile_heating_compressive_none_is_zero():
    assert se.thermal_hoop_stress(400.0, 200.0) < 0.0
    assert se.thermal_hoop_stress(200.0, 250.0) > 0.0
    assert se.thermal_hoop_stress(300.0, 300.0) == 0.0


def test_order_100_mpa_for_200_k_at_12p5_km():
    """A 200 K drop at 400 C rock gives about 100 MPa: comparable to the
    strength, which is why a calibration without it would be biased."""
    assert se.thermal_hoop_stress(400.0, 200.0) / MPa == pytest.approx(-102.3, abs=0.5)


def test_cooling_narrows_breakouts_and_promotes_tensile_fractures():
    SH, Sh, Pw, Pp, T_wall = 260 * MPa, 130 * MPa, 100 * MPa, 95 * MPa, 150.0
    dT = se.thermal_hoop_stress(400.0, T_wall)
    assert (m5.breakout_width(SH, Sh, Pw, Pp, T_wall, dsigma_T=dT)
            < m5.breakout_width(SH, Sh, Pw, Pp, T_wall))
    on = se.mud_window(Sh, SH, Pp, 10000.0, Pw, 400.0, T_wall=T_wall, dsigma_T=dT)
    off = se.mud_window(Sh, SH, Pp, 10000.0, Pw, 400.0, T_wall=T_wall)
    assert on["Pw_tensile"] == pytest.approx(off["Pw_tensile"] + dT)


def test_switch_off_restores_step_5():
    z = CORNWALL.target_depth
    on, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, T_wall=164.0)
    off, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, T_wall=164.0,
                                 thermal=False)
    assert off["dsigma_T"] == 0.0
    assert on["dsigma_T"] == pytest.approx(se.thermal_hoop_stress(400.0, 164.0))
    assert off["window"] == se.mud_window(
        off["Shmin"], off["SHmax"], off["Pp"], z, C.HYDROSTATIC_GRAD * z, 400.0,
        T_wall=164.0, UCS=CORNWALL.UCS)


def test_united_downs_go_with_cooling_conditional_without():
    """At the Model 1 wall temperature (~164 C at 2 kg/s), cooling turns the
    12.5 km target from CONDITIONAL to GO."""
    z = CORNWALL.target_depth
    on, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, T_wall=164.0)
    off, _ = se.stability_inputs(CORNWALL, z, 400.0, UCS=CORNWALL.UCS, T_wall=164.0,
                                 thermal=False)
    assert on["window"]["verdict"] == "GO"
    assert on["window"]["width_hydro"] == pytest.approx(82.0, abs=1.0)
    assert off["window"]["verdict"] == "CONDITIONAL"


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_v10_adapter_ignores_the_thermal_term(site):
    z = site.target_depth
    T = float(site.geotherm(z))
    ref, _ = se.stability_inputs(site, z, T, v10=True, T_wall=100.0, thermal=True)
    Sv, Shmin, SHmax, Pw = se.stresses_v10(site, z)
    assert ref["dsigma_T"] == 0.0
    assert ref["breakout"] == se.breakout_v10(Shmin, SHmax, Pw, T)


# --------------------------------------------- full runs (Model 1, slow)
@pytest.mark.slow
@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_evaluate_uses_model1_wall_and_reports_the_sensitivity(site):
    from comparative_sites import classify, classify_no_thermal
    o = se.evaluate(site)
    assert o["thermal"] is True
    assert o["T_wall"] == pytest.approx(o["m1"]["T_bottom_delivered"])
    assert o["dsigma_T"] == pytest.approx(se.thermal_hoop_stress(o["T_rock"], o["T_wall"]))
    nt = o["no_thermal"]["stress"]
    assert nt["dsigma_T"] == 0.0
    # an uncooled wall: at rock temperature, so no cold-strength gain either
    assert nt["T_wall"] == o["T_rock"]
    assert all(c["T_wall"] == o["T_rock"] for c in o["no_thermal"]["stress_cases"])
    assert classify_no_thermal(o) is not None
    # thermal=False removes only the thermal stress; the strength still sees
    # the cooled wall, so it is not the uncooled sensitivity
    off = se.evaluate(site, thermal=False)
    assert off["dsigma_T"] == 0.0 and off["T_wall"] < off["T_rock"]


@pytest.mark.slow
def test_v10_evaluate_has_no_sensitivity_line():
    from comparative_sites import classify_no_thermal
    o = se.evaluate(CORNWALL, v10=True)
    assert o["thermal"] is False and o["no_thermal"] is None
    assert classify_no_thermal(o) is None


@pytest.mark.slow
@pytest.mark.parametrize("v10", [False, True])
def test_creep_closure_uses_the_same_wall_temperature(v10):
    """Model 4's cooled creep rate takes the wall temperature the stability
    checks use: Model 1's in v1.1, v1.0's 200 C under the adapter."""
    import model4_hole_stability as m4
    o = se.evaluate(CORNWALL, v10=v10)
    T_wall = 200.0 if v10 else o["m1"]["T_bottom_delivered"]
    assert o["T_wall"] == pytest.approx(T_wall)
    expect = m4.closure_rate(o["z"], o["T_rock"], T_wall, 7 * 86400, cooled=True) * se.YEAR * 100
    assert o["creep_cold"] == pytest.approx(expect)


@pytest.mark.slow
def test_sensitivities_cover_uncooled_and_the_calibrated_breakout_limit():
    """The uncooled wall alone, UD-1's 63 deg limit alone, and the two
    together; at United Downs the combination closes the window (NO-GO)."""
    from comparative_sites import classify_sensitivity
    o = se.evaluate(CORNWALL)
    w63 = f"W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f}"
    assert set(o["sensitivity"]) == {"no wall cooling", w63, f"no wall cooling, {w63}"}
    both = o["sensitivity"][f"no wall cooling, {w63}"]["stress"]
    assert both["T_wall"] == o["T_rock"] and both["window"]["W_max"] == 63.0
    assert classify_sensitivity(o, f"no wall cooling, {w63}") == "NO-GO"
