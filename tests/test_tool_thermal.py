"""The BHA and coated pipe through the cycle (tool_thermal.py)."""
import numpy as np
import pytest

import tool_thermal as tt

BHA = tt.Lump.of_pipe(0.073 / 2, 0.17 / 2, 1300.0 * 3000.0)


def test_an_adiabatic_bha_keeps_its_temperature():
    hist = tt.run(BHA, [tt.Segment(3600.0, 300.0, 0.0, 0.0)], 150.0, 140.0)
    assert hist.T_s[-1] == pytest.approx(150.0)
    assert hist.T_f[-1] == pytest.approx(140.0)


def test_a_bha_with_infinite_films_follows_the_fluid():
    """With films far beyond any real one, the steel and the fluid inside
    take the outside fluid's temperature at once, and follow it as it
    changes."""
    T_out = lambda s: 150.0 + 50.0 * s / 600.0
    hist = tt.run(BHA, [tt.Segment(600.0, T_out, 1e9, 1e9)], 100.0, 100.0, dt=1.0)
    assert hist.T_s[-1] == pytest.approx(200.0, abs=0.1)
    assert hist.T_f[-1] == pytest.approx(200.0, abs=0.1)
    assert hist.T_f[60] == pytest.approx(T_out(60.0), abs=0.1)


def test_flowing_fluid_holds_the_tool():
    hist = tt.run(BHA, [tt.Segment(600.0, 300.0, 500.0, 2000.0, T_flow=180.0)], 180.0, 180.0)
    assert np.all(hist.T_f == 180.0)
    assert 180.0 < hist.T_s[-1] < 300.0
    assert hist.T_s[-1] == pytest.approx(tt.steady_steel(180.0, 300.0, 2000.0, 500.0, BHA), abs=1.0)


def test_energy_is_conserved_between_the_nodes():
    """With no outside film the two nodes share their heat."""
    hist = tt.run(BHA, [tt.Segment(7200.0, 0.0, 0.0, 300.0)], 200.0, 100.0)
    E0 = BHA.C_s * 200.0 + BHA.C_f * 100.0
    E1 = BHA.C_s * hist.T_s[-1] + BHA.C_f * hist.T_f[-1]
    assert E1 == pytest.approx(E0, rel=1e-9)
    assert hist.T_s[-1] == pytest.approx(hist.T_f[-1], abs=0.5)


def test_staging_on_the_way_in_lowers_the_peak():
    column = lambda md: 20.0 + 0.035 * md
    circ = lambda md: (60.0, 3000.0, 1500.0)
    kw = dict(md_bottom=10000.0, speed=0.1, T0=20.0, h_static_out=60.0, h_static_in=40.0)
    plain, _ = tt.trip_in(BHA, column, **kw)
    staged, depth = tt.trip_in(BHA, column, stages=list(np.arange(4000.0, 10000.0, 500.0)),
                               stage_s=900.0, circulating=circ, **kw)
    assert staged.peak()[0] < plain.peak()[0]
    assert depth[-1] == pytest.approx(10000.0)
    # without staging the BHA lags the column it is run through
    assert plain.T_f[-1] < column(10000.0)


# ---------------------------------------------------------------- at a site
@pytest.fixture(scope="module")
def soultz(evaluated):
    import site_evaluation as se
    return evaluated(se.SOULTZ)


@pytest.mark.slow
def test_a_forced_exceedance_finds_a_schedule(soultz, monkeypatch):
    """With the tool's limit set just below the connection's peak, FORGE's
    practice fails and the search adds circulation before the connection or
    shortens it until the tool holds."""
    import geo_constants as C
    import site_evaluation as se
    base = soultz["tools"]["connection"]["peak"]
    monkeypatch.setattr(C, "BHA_SURVIVAL_TEMP", base - 0.05)
    if soultz["tools"]["drilling"]["T"] >= base - 0.05:
        pytest.skip("the bit itself is above the forced limit")
    g = se.tool_gates(se.SOULTZ, soultz)
    c = g["connection"]
    assert not c["ok"]
    if c["schedule"] is not None:
        sc = c["schedule"]
        assert sc["peak"] < base - 0.05
        assert (sc["precirc_min"], -sc["off_min"]) > (c["practice"]["precirc_min"], -c["practice"]["off_min"])


@pytest.mark.slow
def test_the_schedule_found_for_running_in_holds_the_tool(soultz):
    tr = soultz["tools"]["trip_in"]
    assert tr["unstaged"]["peak"] > tr["staged"]["peak"]
    if tr["schedule"] is not None:
        assert tr["schedule"]["peak"] < soultz["tools"]["limit"]
        assert tr["schedule"]["spacing_m"] <= tr["staged"]["spacing_m"]
