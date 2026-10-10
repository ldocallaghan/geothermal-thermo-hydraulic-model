"""The new-site template and the validation of a site's inputs."""
import os
import sys
from dataclasses import replace

import numpy as np
import pytest

import site_evaluation as se

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples"))
from site_template import make_site   # noqa: E402


def test_the_template_rebuilds_soultz():
    s = make_site()
    assert s == se.SOULTZ
    assert s.tier == "evidence-based"
    se.validate_site(s)


def test_a_unit_slip_is_caught():
    with pytest.raises(ValueError, match="target_depth = 10.7 m is outside"):
        se.validate_site(replace(se.SOULTZ, target_depth=10.7))
    with pytest.raises(ValueError, match="E_rock"):
        se.validate_site(replace(se.SOULTZ, E_rock=55.0))


def test_temperature_must_rise_down_the_open_hole():
    def geo(z):
        z = np.asarray(z, dtype=float)
        return np.where((z > 6000) & (z < 6500), 250.0, se.urg_geotherm(z))
    with pytest.raises(ValueError, match="does not rise with depth"):
        se.validate_site(replace(se.SOULTZ, geotherm=geo))


def test_stresses_must_keep_the_stated_regime():
    p = replace(se.SOULTZ.stress_cases[0], regime="reverse")
    with pytest.raises(ValueError, match="reverse order"):
        se.validate_site(replace(se.SOULTZ, stress_cases=(p,)))
    p = replace(se.SOULTZ.stress_cases[0], Pp_grad=16e3)
    with pytest.raises(ValueError, match="pore pressure .* is not below Shmin"):
        se.validate_site(replace(se.SOULTZ, stress_cases=(p,)))


def test_a_basis_needs_a_source():
    basis = dict(se.SOULTZ.data_basis, stress=("measured", ""))
    with pytest.raises(ValueError, match="no source"):
        se.validate_site(replace(se.SOULTZ, data_basis=basis))


def test_a_layered_geotherm_follows_its_points():
    g = se.layered_geotherm([(0, 10.0), (2000, 90.0), (5000, 200.0)], 0.03)
    assert float(g(1000)) == pytest.approx(50.0)
    assert float(g(6000)) == pytest.approx(230.0)
