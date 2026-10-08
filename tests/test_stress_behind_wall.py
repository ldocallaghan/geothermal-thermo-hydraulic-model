"""Stresses and failure behind the wall face, against closed forms."""
import numpy as np
import pytest

import model5_convergence_confinement as m5

MPa = 1e6


def test_kirsch_at_the_wall_is_the_hoop_stress_and_the_fluid_pressure():
    th = np.linspace(0, 180, 37)
    sr, st, tau = m5.kirsch(0.1, th, 150 * MPa, 90 * MPa, 50 * MPa, 0.1)
    assert sr == pytest.approx(np.full_like(th, 50 * MPa))
    assert st == pytest.approx(m5.hoop_stress(th, 150 * MPa, 90 * MPa, 50 * MPa))
    assert tau == pytest.approx(np.zeros_like(th), abs=1e-6)


def test_kirsch_hand_checked():
    """SHmax 60, Shmin 40, Pw 20 MPa, a = 1, r = 2: by hand, along SHmax
    sigma_r 44.375, sigma_theta 45.625; at 45 deg 42.5, 57.5 and tau -13.125."""
    sr, st, tau = m5.kirsch(2.0, 0.0, 60.0, 40.0, 20.0, 1.0)
    assert (sr, st, tau) == pytest.approx((44.375, 45.625, 0.0))
    sr, st, tau = m5.kirsch(2.0, 45.0, 60.0, 40.0, 20.0, 1.0)
    assert (sr, st, tau) == pytest.approx((42.5, 57.5, -13.125))


def test_kirsch_tends_to_the_far_field():
    sr, st, tau = m5.kirsch(1e4, 0.0, 60.0, 40.0, 20.0, 1.0)
    assert (sr, st) == pytest.approx((60.0, 40.0), rel=1e-6)


def test_symmetric_about_the_shmin_azimuth():
    for phi in (5.0, 20.0, 40.0):
        a = m5.kirsch(0.15, 90 + phi, 150.0, 90.0, 50.0, 0.1)
        b = m5.kirsch(0.15, 90 - phi, 150.0, 90.0, 50.0, 0.1)
        assert a[0] == pytest.approx(b[0]) and a[1] == pytest.approx(b[1])
        assert a[2] == pytest.approx(-b[2])


def test_thermal_stress_at_the_wall_is_the_existing_hoop_term():
    r = np.geomspace(0.1, 10.0, 400)
    dsr, dst = m5.thermal_stresses(r, np.full_like(r, -100.0), 0.1, 0.5 * MPa)
    assert dsr[0] == pytest.approx(0.0)
    assert dst[0] == pytest.approx(-50 * MPa)


def test_thermal_stress_for_a_cooled_annulus():
    """dT constant from a to b, zero beyond; K = E alpha dT / (1 - nu)."""
    a, b, dT, thermo = 0.1, 0.3, -80.0, 0.5 * MPa
    K = thermo * dT
    r = np.geomspace(a, 3.0, 20001)
    dsr, dst = m5.thermal_stresses(r, np.where(r <= b, dT, 0.0), a, thermo)
    for x in (0.15, 0.25):
        i = np.argmin(abs(r - x))
        assert dsr[i] == pytest.approx(K * (r[i] ** 2 - a ** 2) / (2 * r[i] ** 2), rel=1e-3)
        assert dst[i] == pytest.approx(K * (r[i] ** 2 + a ** 2) / (2 * r[i] ** 2), rel=1e-3)
    for x in (0.5, 1.5):
        i = np.argmin(abs(r - x))
        assert dsr[i] == pytest.approx(K * (b ** 2 - a ** 2) / (2 * r[i] ** 2), rel=1e-3)
        assert dst[i] == pytest.approx(-K * (b ** 2 - a ** 2) / (2 * r[i] ** 2), rel=1e-3)


@pytest.mark.parametrize("cooling", [0.0, 60.0, 150.0])
def test_uniform_temperature_reproduces_the_wall_width(cooling):
    """With the temperature uniform behind the wall, the width judged over the
    rock equals the closed-form width at the wall, within 0.5 deg."""
    import site_evaluation as se
    SH, Sh, Pw, Pp, a, Tr, ucs = 150 * MPa, 90 * MPa, 50 * MPa, 48 * MPa, 0.11, 250.0, 120 * MPa
    Tw = Tr - cooling
    thermo = -se.thermal_hoop_stress(Tr, Tr - 1.0)
    # cooled uniformly over the rock out to well beyond 2a
    b = m5.breakout_behind_wall(SH, Sh, Pw, Pp, a, np.array([a, 1e3]), np.array([Tw, Tw]),
                                Tr, ucs, thermo, r_out=2.0)
    w = m5.breakout_width(SH, Sh, Pw, Pp, Tw, ucs, se.thermal_hoop_stress(Tr, Tw))
    if cooling == 0.0:
        assert b["width"] == pytest.approx(w, abs=0.5)
    assert b["width_wall"] == pytest.approx(w, abs=0.5)


def test_with_no_cooling_failure_starts_at_the_wall_at_every_site():
    """Failure starts at the wall. (The failed region reaches wider just behind
    the wall than at it, so its angular extent is not the breakout width.)"""
    import site_evaluation as se
    from comparative_sites import SITES
    for site in SITES:
        z = site.target()
        T = float(site.temperature()(z))
        for p in site.stress_cases:
            for _, ucs in site.strength_cases:
                Pw = site.rho_fluid_grad * z
                b = m5.breakout_behind_wall(p.SHmax(z), p.Shmin(z), Pw, float(p.Pp(z)), 0.11,
                                            np.array([0.11, 1e3]), np.array([T, T]), T, ucs,
                                            -se.thermal_hoop_stress(T, T - 1.0))
                assert b["r_first"] == pytest.approx(0.11), (site.name, p.label)
                assert b["width"] >= b["width_wall"] - 0.2


def _skin_case(cooling, b):
    """A skin cooled by `cooling` K out to radius b, rock temperature beyond."""
    import site_evaluation as se
    SH, Sh, Pw, Pp, a, Tr, ucs = 285 * MPa, 149 * MPa, 105 * MPa, 105 * MPa, 0.11, 400.0, 100 * MPa
    thermo = -se.thermal_hoop_stress(Tr, Tr - 1.0)
    r = np.array([a, b, b * 1.0001, 1e3])
    T = np.array([Tr - cooling, Tr - cooling, Tr, Tr])
    out = m5.breakout_with_skin(SH, Sh, Pw, Pp, a, r, T, Tr, ucs, thermo)
    wall = lambda Tw: m5.breakout_width(SH, Sh, Pw, Pp, Tw, ucs, se.thermal_hoop_stress(Tr, Tw))
    return out, wall


def test_skin_width_with_no_cooling_is_the_wall_width():
    out, wall = _skin_case(0.0, 0.2)
    assert out["width"] == pytest.approx(wall(400.0))
    assert out["width"] == pytest.approx(out["width_wall"])


def test_a_skin_thicker_than_the_breakout_keeps_its_benefit():
    out, wall = _skin_case(200.0, 0.3)
    assert out["depth_ref"] < 0.3 - 0.11
    assert out["width"] == pytest.approx(wall(200.0))


def test_a_skin_thinner_than_the_breakout_is_judged_by_the_rock_behind_it():
    out, wall = _skin_case(200.0, 0.1101)
    assert out["depth_ref"] > 0.0001
    assert out["width"] == pytest.approx(wall(400.0))
    assert out["width_wall"] == pytest.approx(wall(200.0))
