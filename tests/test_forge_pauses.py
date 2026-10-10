"""The FORGE 16B pause check, pinned. As for the other validations, the tests
hold the results fixed; none asserts the fit is good."""
import json

import pytest

import validate_forge_pauses as vp

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def results():
    return vp.check(True), vp.check(False)


def test_pause_check_is_unchanged(results):
    with open(vp.GOLDEN) as fh:
        pinned = json.load(fh)
    c, c0 = results
    assert c["f"] == pytest.approx(pinned["f"], abs=0.02)
    # The no-node factor depends on the Python version: 0.357 on Python 3.14.4
    # (with the lock file's packages, and with numpy 2.4.4 and scipy 1.17.1),
    # 0.391 on Python 3.11, 3.12 and 3.13. The tolerance covers that spread.
    assert c0["f"] == pytest.approx(pinned["f_no_node"], abs=0.04)
    for now, was in zip(c["rows"] + c0["rows"], pinned["rows"] + pinned["rows_no_node"]):
        assert now["obs"] == pytest.approx(was["obs"], abs=0.5)
        assert now["mod"] == pytest.approx(was["mod"], abs=1.0)


def test_only_trips_are_fitted(results):
    c, _ = results
    for r in c["rows"]:
        assert (r["pred"] is None) == (r["kind"] != "trip")
