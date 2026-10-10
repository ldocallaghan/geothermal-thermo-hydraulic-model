"""The pressure budget (pressure.py): one wall pressure for every state of the
drilling cycle."""
import numpy as np
import pytest
from iapws.iapws97 import _PSat_T

import geo_constants as C
import pressure as pr
import site_evaluation as se
from comparative_sites import CORNWALL

Z = 10700.0
TVD = np.linspace(0.0, Z, 501)


def constant(rho):
    return dict(rho=rho, cp=4180.0, mu=1e-3, k=0.6)


def test_constant_gradient_reproduces_the_earlier_column():
    b = pr.bottomhole_pressure("connection", Z, grad=se.SOULTZ.rho_fluid_grad)
    assert b["P_static"] == se.SOULTZ.rho_fluid_grad * Z
    assert b["SG"] == pytest.approx(C.HYDROSTATIC_GRAD / C.MUD_SG_GRAD)


def test_a_constant_density_column_integrates_to_rho_g_z():
    """With a constant 1,000 kg/m3 fluid, no friction and an open annulus, the
    budget is 1000 g z, an equivalent density of 1.000; at the density of the
    constant gradient it reproduces rho_fluid_grad * z."""
    b = pr.bottomhole_pressure("connection", Z, (TVD, np.full_like(TVD, 20.0)), fluid=constant(1000.0))
    assert b["SG"] == pytest.approx(1.0, abs=1e-9)
    rho = se.SOULTZ.rho_fluid_grad / pr.G
    b = pr.bottomhole_pressure("connection", Z, (TVD, np.full_like(TVD, 20.0)), fluid=constant(rho))
    assert b["P_static"] == pytest.approx(se.SOULTZ.rho_fluid_grad * Z, rel=1e-12)


def fake_run(dp=(1.0e6, 2.0e6, 3.0e6), ends=(0.0, 4000.0, 8000.0, Z)):
    return dict(hyd=dict(pieces=[dict(md_top=a, md_bottom=b, dp_ann=d)
                                 for a, b, d in zip(ends[:-1], ends[1:], dp)]))


def test_annular_friction_accumulates_from_the_surface():
    run = fake_run()
    assert pr.annulus_friction(run, 0.0) == 0.0
    assert pr.annulus_friction(run, 2000.0) == pytest.approx(0.5e6)
    assert pr.annulus_friction(run, 6000.0) == pytest.approx(2.0e6)
    assert pr.annulus_friction(run, Z) == pytest.approx(6.0e6)


def test_drilling_adds_friction_and_a_trip_its_allowance():
    temps = (TVD, np.full_like(TVD, 20.0))
    c = pr.bottomhole_pressure("connection", Z, temps, fluid=constant(1000.0), run=fake_run())
    d = pr.bottomhole_pressure("drilling", Z, temps, fluid=constant(1000.0), run=fake_run())
    t = pr.bottomhole_pressure("trip", Z, temps, fluid=constant(1000.0), allowance_SG=0.02)
    assert c["friction"] == 0.0 and c["P_lo"] == c["P_hi"] == c["P_static"]
    assert d["P_lo"] == d["P_hi"] == pytest.approx(c["P_static"] + 6.0e6)
    assert d["SG"] == c["SG"]          # the equivalent density is the static column
    assert t["P_lo"] == pytest.approx(c["P_static"] - 0.02 * C.MUD_SG_GRAD * Z)
    assert t["P_hi"] == pytest.approx(c["P_static"] + 0.02 * C.MUD_SG_GRAD * Z)


def test_the_equivalent_density_leaves_out_the_wellhead_pressure():
    temps = (TVD, np.full_like(TVD, 20.0))
    b = pr.bottomhole_pressure("connection", Z, temps, fluid=constant(1000.0), P_surface=2.0e6)
    assert b["P_static"] == pytest.approx(1000.0 * pr.G * Z + 2.0e6)
    assert b["SG"] == pytest.approx(1.0)


def test_backpressure_only_when_the_return_would_boil():
    assert pr.backpressure(73.0) == 0.0
    assert pr.backpressure(99.0) == 0.0
    assert pr.backpressure(120.0) == pytest.approx(_PSat_T(393.15) * 1e6 - pr.P_ATM)
    assert pr.backpressure(120.0) > 0.0


def test_the_window_in_equivalent_density_leaves_out_the_wellhead_pressure():
    """A wellhead pressure moves the wall pressure, not the bounds: the window
    in equivalent density shifts down by it."""
    z, T = 12500.0, 400.0
    ref0, _ = se.stability_inputs(CORNWALL, z, T, UCS=CORNWALL.UCS)
    Pw = ref0["window"]["Pw_hydro"]
    ref1, _ = se.stability_inputs(CORNWALL, z, T, UCS=CORNWALL.UCS, Pw=Pw + 1.0e6, P_surface=1.0e6)
    w0, w1 = ref0["window"], ref1["window"]
    shift = 1.0e6 / (C.MUD_SG_GRAD * z)
    assert w1["SG_hydro"] == pytest.approx(w0["SG_hydro"])
    assert w1["SG_lo"] == pytest.approx(w0["SG_lo"] - shift)
    assert w1["SG_hi"] == pytest.approx(w0["SG_hi"] - shift)


# ---------------------------------------------------------------- with Model 1
@pytest.mark.slow
@pytest.mark.parametrize("site,circ,formation", [(se.SOULTZ, 0.960, 0.859),
                                                 (CORNWALL, 0.961, 0.876)])
def test_water_columns_match_the_review(evaluated, site, circ, formation):
    """Water at the annulus temperatures of the base flow, with the rock
    exposed for a year, and at formation temperature (a bound on a long trip),
    as the review integrated them."""
    o = evaluated(site)
    z = o["z"]
    water = se.circulate(site, z, site.temperature(), o["circ"]["pipe"], float(o["circ"]["m_dot"]))
    b = pr.bottomhole_pressure("connection", z, pr.circulating_temperatures(water))
    assert b["SG"] == pytest.approx(circ, abs=0.005)
    tvd = np.linspace(0.0, z, 1001)
    b = pr.bottomhole_pressure("trip", z, (tvd, [site.temperature()(x) for x in tvd]))
    assert b["SG"] == pytest.approx(formation, abs=0.005)


@pytest.mark.slow
def test_the_trip_column_is_hotter_than_the_circulating_one(evaluated):
    o = evaluated(se.SOULTZ)
    c = o["cycle"]["budget"]
    tvd_t, T_t = c["temperatures"]["trip"]
    tvd_c, T_c = c["temperatures"]["circulating"]
    shoe = 4000.0
    assert np.all(np.diff(tvd_t) > 0)
    below = tvd_t > shoe
    assert np.all(T_t[below] >= np.interp(tvd_t[below], tvd_c, T_c) - 1e-6)
    assert np.max(np.diff(tvd_t[tvd_t >= shoe])) <= se.TRIP_COLUMN_DZ + 1e-6
    assert c["water_SG_trip"] < c["water_SG"]
    assert c["trip_SG_start"] > o["cycle"]["trip_SG"]


@pytest.mark.slow
def test_a_boiling_return_is_held_with_backpressure(monkeypatch):
    """If the return would boil at an open wellhead, the run is repeated with
    the backpressure that keeps it liquid, and the cycle carries it."""
    monkeypatch.setattr(pr, "backpressure", lambda T: 2.0e6)
    z = 5000.0
    geotherm = lambda d: 15.0 + 0.06 * d
    r = se.circulate(se.SOULTZ, z, geotherm, se.SITE_PIPES[0], 20.0)
    assert r["backpressure"] == 2.0e6
    assert r["P_surface"] == pytest.approx(pr.P_ATM + 2.0e6)
    Pw, bp = se.water_column(se.SOULTZ, r, z)
    assert bp == 2.0e6
    assert (Pw - bp) / (C.MUD_SG_GRAD * z) < 1.0
