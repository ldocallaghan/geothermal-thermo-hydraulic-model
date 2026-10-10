"""The v1.3.3 site table, pinned (made in the CI environment, on the lock
file's versions): the status, the cycle and its sensitivities, the common
windows and their binding depths, the margins, and the tools and coating."""
import json
import os

import pytest

from comparative_sites import SITES, golden_rows_v133
from conftest import GOLDEN

pytestmark = pytest.mark.slow
PATH = os.path.join(GOLDEN, "v133_site_table.json")


@pytest.fixture(scope="module")
def rows(evaluated):
    if not os.path.exists(PATH):
        pytest.skip("no v1.3.3 golden yet")
    with open(PATH) as fh:
        pinned = json.load(fh)
    return pinned, golden_rows_v133([evaluated(s) for s in SITES])


def close(a, b, tol):
    if a is None or b is None or isinstance(a, str):
        return a == b
    return abs(a - b) <= tol


def test_statuses_and_verdicts_are_unchanged(rows):
    pinned, now = rows
    assert now.keys() == pinned.keys()
    for site, row in pinned.items():
        n = now[site]
        assert n["status"] == row["status"], site
        assert n["verdicts"] == row["verdicts"], site
        assert n["site_verdict"] == row["site_verdict"], site
        assert n["deciding_case"] == row["deciding_case"], site
        assert n["tools"]["tool"] == row["tools"]["tool"], site
        assert n["tools"]["coating"] == row["tools"]["coating"], site
        assert n["flow_ceiling"] == row["flow_ceiling"], site
        assert n["flow"] == row["flow"] and n["fluid_passes"] == row["fluid_passes"], site


def test_fluids_bounds_and_pauses_are_unchanged(rows):
    pinned, now = rows
    for site, row in pinned.items():
        n = now[site]
        for k in ("drill_SG", "drill_SG_hi", "trip_SG", "trip_SG_start", "trip_SG_hi", "ecd_SG"):
            assert close(n[k], row[k], 0.005), (site, k)
        assert close(n["fluid_SG"], row["fluid_SG"], 0.005), site
        assert close(n["safe_pause_h"], row["safe_pause_h"], max(1.0, 0.05 * (row["safe_pause_h"] or 0))), site
        assert close(n["trip_h"], row["trip_h"], 0.1), site
        assert close(n["rop_m_h"], row["rop_m_h"], 0.05), site
        for k, v in row["sensitivities"].items():
            w = n["sensitivities"][k]
            assert w["site_verdict"] == v["site_verdict"], (site, k)
            assert close(w["drill_SG"], v["drill_SG"], 0.01), (site, k)
        for (h, v, sg), (h2, v2, sg2) in zip(row["trip_profile"], n["trip_profile"]):
            assert h == h2 and v == v2 and close(sg, sg2, 0.01), (site, h)


def test_common_windows_and_binding_depths_are_unchanged(rows):
    pinned, now = rows
    for site, row in pinned.items():
        for k, w in row["windows"].items():
            x = now[site]["windows"][k]
            assert x["open"] == w["open"] and x["note"] == w["note"], (site, k)
            assert close(x["lo"], w["lo"], 5.0) and close(x["hi"], w["hi"], 5.0), (site, k)
            for b in ("lo_by", "hi_by"):
                assert x[b]["state"] == w[b]["state"] and x[b]["case"] == w[b]["case"], (site, k, b)
                assert close(x[b]["md"], w[b]["md"], 101.0), (site, k, b)


def test_margins_and_tools_are_unchanged(rows):
    pinned, now = rows
    for site, row in pinned.items():
        n = now[site]
        for tag, d in row["shmin_margin"].items():
            for k, v in d.items():
                assert close(n["shmin_margin"][tag][k], v, 0.01), (site, tag, k)
        for k in ("bit", "connection", "running_in", "running_in_unstaged", "coating_max"):
            assert close(n["tools"][k], row["tools"][k], 2.0), (site, k)
