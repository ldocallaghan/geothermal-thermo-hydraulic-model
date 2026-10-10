"""Common windows across scenarios and down the hole."""
import pytest

import comparative_sites as cs
import site_evaluation as se

pytestmark = pytest.mark.slow


def test_pannonian_has_no_common_drilling_fluid(evaluated):
    """Each of Pannonian's Shmin cases has a drilling window, but no fluid
    suits them all: the deciding pair is two different cases."""
    o = evaluated(cs.PANNONIAN)
    w = o["windows"]["drilling_fluid"]
    assert not w["open"] and w["fluid"] is None
    a, b = w["deciding_pair"]
    assert a != b
    assert w["note"] == "none: the scenarios need different fluids"
    for lab in (a, b):
        lo = max(o["windows"]["states"][st][lab]["lo"] for st in ("drilling", "connection", "fresh rock"))
        hi = min(o["windows"]["states"][st][lab]["hi"] for st in ("drilling", "connection", "fresh rock"))
        assert lo <= hi, lab


def test_soultz_reports_where_the_trip_upper_bound_binds(evaluated):
    o = evaluated(se.SOULTZ)
    w = o["windows"]["trip_fluid"]
    shoe = max(h.md_top for h in o["m1"]["geometry"].hole)
    assert w["hi_by"]["state"] == "trip"
    assert shoe <= w["hi_by"]["md"] <= o["z"]
    assert len(o["windows"]["depths"]) >= (o["z"] - shoe) / se.WINDOW_DZ


def test_the_bounds_hold_the_cycle_at_the_bottom(evaluated):
    """The drilling fluid the cycle requires at the bottom of the deciding
    case is no heavier at surface than the common lower bound."""
    o = evaluated(se.SOULTZ)
    rho = o["cycle"]["budget"]["drill_rho_surface"] or 998.0
    assert o["windows"]["drilling_fluid"]["lo"] >= rho - 5.0


def test_the_cycle_table_prints(evaluated, capsys):
    rs = [evaluated(se.SOULTZ), evaluated(cs.PANNONIAN)]
    cs.print_cycle("test", [se.SOULTZ, cs.PANNONIAN], rs)
    out = capsys.readouterr().out
    assert "Common windows" in out and "Status:" in out and "Tools and coating" in out
