"""The site table with the validated circulation model, pinned, and the rules
that choose each site's flow and pipe."""
import json
import os

import pytest

import model1_coupled as m1
import site_evaluation as se
from comparative_sites import SITES, golden_rows
from conftest import GOLDEN

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def results():
    return [se.evaluate(s) for s in SITES]


def test_site_table_is_unchanged(results):
    with open(os.path.join(GOLDEN, "v12_site_table.json")) as fh:
        pinned = json.load(fh)
    now = golden_rows(results)
    assert now.keys() == pinned.keys()
    for site, row in pinned.items():
        for k in ("verdict", "pipe", "m_dot", "sensitivities"):
            assert now[site][k] == row[k], (site, k)
        for k in ("bhct", "bhct_one_year", "spp_MPa", "v_ann", "MW"):
            assert now[site][k] == pytest.approx(row[k], abs=0.5), (site, k)


def test_flow_rule(results):
    """The flow meets survival, hole cleaning and the pump limit, unless the
    run reports a conflict."""
    for o in results:
        c = o["circ"]
        if c["conflict"]:
            continue
        assert c["bhct"] < se.C.BHA_SURVIVAL_TEMP
        assert c["hyd"]["v_ann_min"] >= m1.HOLE_CLEANING_V
        assert c["hyd"]["spp"] <= m1.PUMP_LIMIT


def test_wall_temperature_is_the_chosen_runs_bhct(results):
    for o in results:
        assert o["T_wall"] == pytest.approx(o["circ"]["bhct"])
        assert o["m1"]["T_bottom_delivered"] == pytest.approx(o["circ"]["bhct"])


def test_fresh_wall_runs_hotter_than_a_years_exposure(results):
    for o in results:
        assert o["circ"]["bhct"] >= o["circ"]["bhct_one_year"] - 0.5
