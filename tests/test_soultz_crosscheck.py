"""Cross-check of the stability model against the Soultz wells
(crosscheck_soultz.py).

The model must predict breakouts at 5 km with an open mud window, as in GPK3
and GPK4. Nothing is fitted to Soultz; the comparison targets are the
observations in Valley & Evans (2007).
"""
import pytest

import crosscheck_soultz as cs
from comparative_sites import CORNWALL

MPa = 1.0e6


@pytest.fixture(scope="module")
def rows():
    return cs.crosscheck()


def test_breakouts_at_5_km_with_an_open_window(rows):
    """For v1.0's Soultz UCS and the measured lab range, every case breaks
    out at 5 km; for every strength tried, the window stays open."""
    for r in rows:
        if r["source"] in ("v1.0 UCS (unsourced)", "Soultz lab UCS 100-130"):
            assert all(w > 0 for w in r["width"]), (r["source"], r["ucs"], r["cooling"])
        assert all(x["open"] for x in r["window"])


def test_only_the_strongest_intact_rock_escapes_breakout(rows):
    misses = [(r["source"], r["cooling"]) for r in rows for w in r["width"] if w == 0.0]
    assert misses and all(src == "UD-1 intact bound" and dK == 40.0
                          for src, dK in misses)


def test_lab_strength_reproduces_the_gpk4_onset_at_the_shmax_lower_bound(rows):
    """Valley & Evans needed SHmax >= 0.9 Sv to put the breakout onset at
    3,670 m TVD with the lab UCS; the model, with the same strength and the
    lower-bound SHmax, puts it within 100 m of that."""
    r = next(r for r in rows if r["ucs"] == 115e6 and r["cooling"] == 40.0)
    assert r["onset"][0] == pytest.approx(cs.OBS_ONSET_DENSE, abs=100.0)


def test_ud1_weak_zone_strength_brackets_the_observed_onset(rows):
    """Transferred from UD-1 unchanged, the weak-zone strengths put the
    mid-range onset between 2.6 and 3.9 km, around GPK4's 3.0-3.67 km."""
    onsets = [r["onset"][1] for r in rows if r["source"] == "UD-1 weak zones"]
    assert min(onsets) < cs.OBS_ONSET_FIRST < cs.OBS_ONSET_DENSE < max(onsets)


def test_strength_sources_are_the_documented_ones():
    src = dict(cs.STRENGTHS)
    assert src["Soultz lab UCS 100-130"] == (100e6, 115e6, 130e6)
    assert set(src["UD-1 weak zones"]) == {u for l, u in CORNWALL.strength_cases
                                               if "weak" in l}


def test_tensile_initiation_does_not_discriminate_with_t0_zero():
    """Observed tensile fractures are continuous only to 2,180 m TVD; with
    zero tensile strength the model initiates them at every depth checked."""
    zs = cs.tensile_depths(cs.SOULTZ.stress_cases[1], 0.0)
    assert zs[0] == 1500.0 and zs[-1] == 5000.0


def test_results_are_unchanged_from_v122():
    """Every row of the cross-check (widths at 5 km, onsets, verdicts), as
    v1.2.2 gave it; the wall-only stability check is held fixed."""
    import json, os
    from conftest import GOLDEN
    with open(os.path.join(GOLDEN, "soultz_crosscheck_v122.json")) as fh:
        pinned = json.load(fh)
    rows = cs.crosscheck()
    assert len(rows) == len(pinned)
    for r, p in zip(rows, pinned):
        assert (r["source"], r["ucs"] / 1e6, r["cooling"]) == (p["source"], p["ucs"], p["cooling"])
        assert r["width"] == pytest.approx(p["width"], abs=1e-3)
        assert r["onset"] == p["onset"]
        assert [x["verdict"] for x in r["window"]] == p["verdicts"]
