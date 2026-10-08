"""Pressure losses, bit nozzles, standpipe pressure and annular velocity."""
import numpy as np
import pytest

import model1_coupled as m1
import well_geometry as wg

PSI = 6894.757
GPM = 6.30902e-5   # m^3/s


def test_laminar_friction_is_hagen_poiseuille():
    assert m1.friction_factor(1000.0) == pytest.approx(0.064)


def test_turbulent_friction_is_blasius():
    assert m1.friction_factor(5e4) == pytest.approx(0.316 * 5e4 ** -0.25)


def test_blasius_against_smooth_colebrook():
    """Blasius holds to a few percent up to Re ~1e5 and under-reads above."""
    for Re in (1e4, 1e5):
        assert m1.friction_factor(Re) == pytest.approx(m1.colebrook_smooth(Re), rel=0.03)
    assert m1.friction_factor(1e6) < 0.9 * m1.colebrook_smooth(1e6)


def test_colebrook_smooth_known_value():
    # Moody chart, smooth pipe, Re = 1e5: f = 0.0180
    assert m1.colebrook_smooth(1e5) == pytest.approx(0.0180, abs=2e-4)


def test_bit_nozzles_reproduce_the_forge_rig():
    """Eight 14/32-inch nozzles, 600 gal/min of 8.40 lb/gal mud: the rig's
    daily report gives 193 psi."""
    rho = 8.40 * 119.826
    dp = m1.nozzle_dp(600 * GPM * rho, rho)
    assert dp / PSI == pytest.approx(193, abs=2)


def _pipe_in_hole(L=3000.0):
    return wg.WellGeometry.single(L, wg.CONVENTIONAL, r_hole=0.108)


@pytest.mark.slow
def test_string_loss_laminar_and_turbulent():
    import water_table as wt
    g = _pipe_in_hole()
    T, P = 60.0, m1.P_of_z(1500.0)
    rho, _, mu, _ = (float(x[0]) for x in wt.props(np.array([T]), np.array([P])))
    D, L = 2 * wg.CONVENTIONAL.r_bore, 3000.0
    A = np.pi * D ** 2 / 4

    m_dot = 0.05   # laminar
    Q = m_dot / rho
    hp = 128 * mu * L * Q / (np.pi * D ** 4)
    assert m1.hydraulics(m_dot, g, T_avg=T)["dp_string"] == pytest.approx(hp, rel=1e-9)

    m_dot = 30.0   # turbulent
    v = m_dot / (rho * A)
    blasius = 0.316 * (rho * v * D / mu) ** -0.25 * L / D * 0.5 * rho * v ** 2
    assert m1.hydraulics(m_dot, g, T_avg=T)["dp_string"] == pytest.approx(blasius, rel=1e-9)


@pytest.mark.slow
def test_standpipe_pressure_rises_with_flow():
    g = _pipe_in_hole()
    spp = [m1.hydraulics(md, g)["spp"] for md in (5, 10, 20, 40, 60)]
    assert np.all(np.diff(spp) > 0)


@pytest.mark.slow
def test_annular_velocity_and_hole_cleaning():
    g = _pipe_in_hole()
    h = m1.hydraulics(40.0, g, T_avg=60.0)
    A = np.pi * (0.108 ** 2 - wg.CONVENTIONAL.r_out ** 2)
    rho = 40.0 / (h["v_ann_min"] * A)
    assert 970 < rho < 1000
    assert h["cleans_hole"] == (h["v_ann_min"] >= m1.HOLE_CLEANING_V)


@pytest.mark.slow
def test_hot_annulus_lowers_standpipe_pressure():
    """With a solved temperature profile the lighter, hotter annulus
    returns some of the pump's work: the buoyancy term is negative."""
    r = m1.solve(m_dot=10.0, pipe=wg.LEGACY_VACUUM, verbose=False)
    h = m1.hydraulics(10.0, r["geometry"], result=r)
    assert h["dp_buoyancy"] < 0
    assert h["spp"] == pytest.approx(h["dp_string"] + h["dp_annulus"] + h["dp_bit"]
                                     + h["dp_buoyancy"])


@pytest.mark.slow
def test_friction_heat_in_an_insulated_string():
    """With a near-adiabatic string, the bore warms by its friction loss over
    rho cp, give or take what crosses the wall."""
    import benchmark_utaustin as b
    on, off = b.run(wg.DUAL_WALL), b.run(wg.DUAL_WALL, friction_heat=False)
    expect = on["hyd"]["dp_string"] / (b.MUD["rho"] * b.MUD["cp"])
    assert b.bhct(on) - b.bhct(off) == pytest.approx(expect, rel=0.15)
