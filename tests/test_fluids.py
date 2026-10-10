"""The drilling fluid (fluids.py): a barite-weighted water-based mud, its
water limit, and its Bingham-plastic pressure loss."""
import numpy as np
import pytest

import fluids

PSI_PER_FT = 6894.757 / 0.3048     # Pa/m
CP = 1e-3                           # Pa s
LBF100 = fluids.LBF_100FT2


def test_bingham_annular_loss_by_hand():
    """Laminar flow in one annulus section, in field units (Bourgoyne et al.
    1986, ch. 4): an 8.5 in hole around 5 in pipe, mean velocity 1 ft/s,
    PV 27 cP, YP 28 lbf/100 ft2:
        dp/dL = 27 x 1 / (1000 x 3.5^2) + 28 / (200 x 3.5)
              = 0.0022041 + 0.0400000 = 0.0422041 psi/ft."""
    hand = (27 * 1.0 / (1000 * 3.5 ** 2) + 28 / (200 * 3.5)) * PSI_PER_FT
    grad, Re_a = fluids.bingham_gradient(1300.0, 27 * CP, 28 * LBF100, 0.3048, 3.5 * 0.0254,
                                         annulus=True)
    assert Re_a < fluids.RE_CRITICAL
    assert float(grad) == pytest.approx(hand, rel=1e-3)


def test_bingham_pipe_loss_by_hand():
    """Laminar flow in a 3 in bore at 2 ft/s, PV 27 cP, YP 28 lbf/100 ft2:
        dp/dL = 27 x 2 / (1500 x 3^2) + 28 / (225 x 3) psi/ft."""
    hand = (27 * 2.0 / (1500 * 3.0 ** 2) + 28 / (225 * 3.0)) * PSI_PER_FT
    grad, Re_a = fluids.bingham_gradient(1300.0, 27 * CP, 28 * LBF100, 2 * 0.3048, 3.0 * 0.0254,
                                         annulus=False)
    assert Re_a < fluids.RE_CRITICAL
    assert float(grad) == pytest.approx(hand, rel=1e-3)


def test_without_a_yield_point_the_turbulent_loss_is_blasius():
    rho, mu, v, D = 1000.0, 1e-3, 2.0, 0.1
    grad, Re_a = fluids.bingham_gradient(rho, mu, 0.0, v, D, annulus=True)
    Re = rho * v * D / mu
    assert Re_a == pytest.approx(Re)
    assert float(grad) == pytest.approx(0.316 * Re ** -0.25 * rho * v ** 2 / (2 * D))


def test_a_yield_point_keeps_slow_flow_laminar():
    """At the same velocity, the yield point raises the apparent viscosity and
    can hold the flow laminar where water would be turbulent."""
    _, Re_w = fluids.bingham_gradient(1300.0, 5 * CP, 0.0, 1.5, 0.09, annulus=True)
    _, Re_m = fluids.bingham_gradient(1300.0, 5 * CP, 28 * LBF100, 1.5, 0.09, annulus=True)
    assert Re_w > fluids.RE_CRITICAL > Re_m


# ---------------------------------------------------------- with water properties
@pytest.mark.slow
def test_the_mud_without_barite_or_yield_point_is_water():
    import water_table as wt
    T = np.array([20.0, 150.0, 300.0])
    P = np.array([1e5, 50e6, 100e6])
    w = fluids.water()
    assert w.phi == pytest.approx(0.0, abs=1e-12)
    for a, b in zip(w.props(T, P), wt.props(T, P)):
        assert np.allclose(a, b, rtol=1e-12)


@pytest.mark.slow
def test_a_mud_has_its_surface_weight_and_expands_less_than_water():
    m = fluids.Mud(1400.0)
    rho0 = m.props(np.array([fluids.T_REF]), np.array([fluids.P_ATM]))[0][0]
    assert rho0 == pytest.approx(1400.0, rel=1e-3)
    w = fluids.water()
    T, P = np.array([fluids.T_REF, 300.0]), np.array([fluids.P_ATM, 100e6])
    exp_m = m.props(T, P)[0]
    exp_w = w.props(T, P)[0]
    assert exp_m[1] / exp_m[0] > exp_w[1] / exp_w[0]       # loses a smaller fraction
    with pytest.raises(ValueError):
        fluids.Mud(900.0)


@pytest.mark.slow
def test_the_plastic_viscosity_follows_water():
    """PV at the reference temperature is the measured 27 cP, and it falls
    with temperature as water's viscosity does."""
    m = fluids.Mud(1300.0)
    mu = m.props(np.array([fluids.PV_REF_T, 93.3]), np.array([fluids.P_ATM, fluids.P_ATM]))[2]
    assert mu[0] == pytest.approx(fluids.PV_REF, rel=1e-3)
    # Anawe and Folayan (2018): 9 cP at 200 F
    assert mu[1] == pytest.approx(9 * CP, rel=0.1)


@pytest.mark.slow
def test_water_limit_reproduces_the_water_hydraulics():
    """With water's density and no yield point, the fluid gives v1.3.2's
    hydraulics to within 1%."""
    import model1_coupled as m1
    import site_evaluation as se
    g = se.site_well(10700.0, se.SITE_PIPES[0])
    a = m1.hydraulics(40.0, g, T_avg=80.0)
    b = m1.hydraulics(40.0, g, T_avg=80.0, fluid=fluids.water())
    for k in ("dp_string", "dp_annulus", "dp_bit", "spp"):
        assert b[k] == pytest.approx(a[k], rel=0.01), k
