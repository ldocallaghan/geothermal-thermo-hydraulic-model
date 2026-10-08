"""The FORGE 16B validation, pinned. These hold the fitted values, the BHA 12
blind prediction and the decomposition fixed so that any change shows; no
test asserts that the fit is good."""
import json

import pytest

import validate_forge16b as v

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def results():
    return v.golden(v.run_all())


@pytest.fixture(scope="module")
def pinned():
    with open(v.GOLDEN) as fh:
        return json.load(fh)


def test_fitted_values(results, pinned):
    assert results["film_mult"] == pytest.approx(pinned["film_mult"], abs=0.02)
    assert results["k_idp"] == pytest.approx(pinned["k_idp"], rel=0.03)


def test_run_means_and_blind_prediction(results, pinned):
    for run, t in pinned["run_means"].items():
        assert results["run_means"][run] == pytest.approx(t, abs=1.0), run
    assert results["bha12_late"] == pytest.approx(pinned["bha12_late"], abs=1.0)


def test_decomposition(results, pinned):
    for name, d in pinned["decomposition"].items():
        assert results["decomposition"][name] == pytest.approx(d, abs=1.0), name


def test_mwd_sensor_sits_in_the_bha():
    for run in v.RUNS:
        bha = float(v.RUNS[run]["bha_nmdc_ft"])
        assert 0 < v.MWD_ABOVE_BIT_FT < bha
