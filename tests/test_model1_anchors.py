"""Model 1 results as they stand before the circulation model is extended.

These pin the coupled counterflow solve at its current settings: a single
vertical well, one insulation layer, and rock exposed for a year. Any later
change to the geometry, insulation or exposure must reproduce them when run
with these settings.

The boundary-value solve moves at about the 1e-3 level with numpy/scipy
versions, so temperatures are held to 0.5 C and heat to 1%.
"""
import pytest

import model1_coupled as m1
import site_evaluation as se
import well_geometry as wg

pytestmark = pytest.mark.slow


@pytest.mark.parametrize("m_dot, k_ins, bit, ret, mw", [
    (2.0, 0.10, 317.747, 57.819, 0.1486),
    (2.0, 0.02, 170.557, 81.713, 0.3483),
    (10.0, 0.02, 59.154, 135.068, 3.9855),
])
def test_450c_well_at_12p4_km(m_dot, k_ins, bit, ret, mw):
    """The README's Fig. 1: 35 C/km to 450 C at 12.4 km, 40 C inlet."""
    r = m1.solve(m_dot=m_dot, pipe=wg.uniform_wall(k_ins), verbose=False)
    assert r["L"] == pytest.approx(12428.6, abs=0.5)
    assert r["T_bottom_delivered"] == pytest.approx(bit, abs=0.5)
    assert r["T_return_surface"] == pytest.approx(ret, abs=0.5)
    assert r["Q_product"] / 1e6 == pytest.approx(mw, rel=0.01)


def test_soultz_at_10_kg_s():
    """The Soultz evaluation's production run: 10.7 km, 10 kg/s."""
    s = se.SOULTZ
    r = m1.solve(m_dot=10.0, T_inj=s.T_inj, pipe=wg.LEGACY_VACUUM, k_rock=s.k_rock,
                 geotherm=s.temperature(), target_depth=s.target(),
                 Q_face=30000.0, verbose=False)
    assert r["T_bottom_delivered"] == pytest.approx(54.937, abs=0.5)
    assert r["T_return_surface"] == pytest.approx(136.750, abs=0.5)
    assert r["Q_product"] / 1e6 == pytest.approx(4.0567, rel=0.01)


def test_soultz_wall_temperature_for_the_verdict():
    """The wall temperature the stability verdict uses: the bit temperature
    at the minimum flow that keeps the tools below 200 C."""
    o = se.evaluate(se.SOULTZ)
    assert o["m_min"] == 2
    assert o["T_wall"] == pytest.approx(148.088, abs=0.5)
