"""Verdicts through the drilling cycle."""
import json
import os

import numpy as np
import pytest

import site_evaluation as se
from comparative_sites import SITES, golden_rows_v13
from conftest import GOLDEN

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def soultz():
    return se.evaluate(se.SOULTZ)


def test_safe_pause_rises_with_the_fluid_weight(soultz):
    c = soultz["cycle"]
    pauses = [c["safe_pause_at"](sg) for sg in (c["drill_SG"], c["drill_SG"] + 0.05, c["drill_SG"] + 0.1)]
    assert pauses == sorted(pauses)


def test_zero_pause_is_the_drilling_field_at_static_pressure(soultz):
    c = soultz["cycle"]
    T0 = float(np.interp(c["r_ref"], c["result"].r, c["result"].fields[3]))
    assert c["after_pause"](0.0) == pytest.approx(T0, abs=0.5)


def test_an_uncooled_state_gives_the_no_wall_cooling_verdict(soultz):
    """The trip state with the wall back at rock temperature is the "no wall
    cooling" sensitivity."""
    ref, _ = se.stability_inputs(se.SOULTZ, soultz["z"], soultz["T_rock"], T_wall=soultz["T_rock"])
    nwc = soultz["sensitivity"]["no wall cooling"]["stress"]["window"]
    assert ref["window"]["verdict"] == nwc["verdict"]
    assert ref["window"]["SG_lo"] == pytest.approx(nwc["SG_lo"], abs=0.01)


def test_the_trip_needs_more_than_drilling(soultz):
    c = soultz["cycle"]
    assert c["T_ref"]["trip"] > c["T_ref"]["drilling"]
    assert c["trip_SG"] >= c["drill_SG"] - 1e-6


def test_the_safe_pause_is_bracketed(soultz):
    c = soultz["cycle"]
    h = c["safe_pause_h"]
    if 0 < h < se.SAFE_PAUSE_MAX_H:
        assert c["width_ok"](c["after_pause"](0.98 * h), c["drill_SG"] * se.C.MUD_SG_GRAD * c["z"])
        assert not c["width_ok"](c["after_pause"](1.02 * h), c["drill_SG"] * se.C.MUD_SG_GRAD * c["z"])


def test_site_table_v13_is_unchanged():
    with open(os.path.join(GOLDEN, "v13_site_table.json")) as fh:
        pinned = json.load(fh)
    now = golden_rows_v13([se.evaluate(s) for s in SITES])
    assert now.keys() == pinned.keys()
    for site, row in pinned.items():
        assert now[site]["verdicts"] == row["verdicts"], site
        assert now[site]["site_verdict"] == row["site_verdict"], site
        for k in ("drill_SG", "trip_SG", "trip_SG_hi"):
            assert now[site][k] == pytest.approx(row[k], abs=0.01), (site, k)
        assert now[site]["trip_h"] == pytest.approx(row["trip_h"], abs=0.1)
