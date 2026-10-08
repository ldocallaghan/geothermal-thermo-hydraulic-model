"""Exposure time of the wall to circulation, varying with depth."""
import numpy as np
import pytest

import model1_coupled as m1
import well_geometry as wg


def test_constant_rate_of_penetration():
    t = m1.exposure_from_rop(md_bit=3000.0, rop_m_per_h=10.0, t_circulating=7200.0)
    assert t(np.array([3000.0, 2900.0, 0.0])) == pytest.approx(
        [7200.0, 7200.0 + 36000.0, 7200.0 + 1.08e6])


def test_drilling_history_uses_first_arrival():
    """The bit drills to 100 m, trips out to 40 m, then drills on to 200 m:
    the wall from 40 to 100 m keeps its first exposure."""
    time = [0.0, 100.0, 150.0, 200.0, 300.0]
    md = [0.0, 100.0, 40.0, 100.0, 200.0]
    t = m1.exposure_from_history(time, md, t_now=400.0)
    assert t(np.array([50.0, 100.0, 150.0, 200.0])) == pytest.approx(
        [350.0, 300.0, 150.0, 100.0])


def test_floor_of_one_hour():
    z = np.array([0.0, 10.0])
    t = m1.exposure_seconds(z, None, lambda s: np.array([10.0, 1e6]))
    assert t == pytest.approx([m1.T_EXPOSURE_FLOOR, 1e6])


def test_t_years_path_is_unchanged():
    t = m1.exposure_seconds(np.zeros(3), 1.0, None)
    assert t == pytest.approx(np.full(3, 3.1536e7))


@pytest.mark.slow
def test_constant_one_year_reproduces_t_years():
    a = m1.solve(m_dot=10.0, pipe=wg.LEGACY_VACUUM, verbose=False)
    b = m1.solve(m_dot=10.0, pipe=wg.LEGACY_VACUUM, verbose=False,
                 exposure=lambda s: np.full_like(np.asarray(s, float), 3.1536e7))
    assert b["T_bottom_delivered"] == pytest.approx(a["T_bottom_delivered"], abs=1e-6)
    assert b["T_return_surface"] == pytest.approx(a["T_return_surface"], abs=1e-6)


@pytest.mark.slow
def test_freshly_drilled_hole_is_hotter_at_the_bit():
    """Rock that has had less time to cool supplies more heat, so the bit runs hotter."""
    year = m1.solve(m_dot=10.0, pipe=wg.LEGACY_VACUUM, verbose=False)
    fresh = m1.solve(m_dot=10.0, pipe=wg.LEGACY_VACUUM, verbose=False,
                     exposure=m1.exposure_from_rop(year["L"], 5.0, 3600.0))
    assert fresh["T_bottom_delivered"] > year["T_bottom_delivered"] + 1.0
