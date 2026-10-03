"""C9: Soultz cross-check (spec v1.1), and spec test 6.

Test 6: breakouts at 5 km with a non-empty window. Nothing is fitted to Soultz;
the comparison targets are the GPK3/GPK4 observations in Valley & Evans (2007).
"""
import pytest

import crosscheck_soultz as cs
from comparative_sites import CORNWALL

MPa = 1.0e6


@pytest.fixture(scope="module")
def rows():
    return cs.crosscheck()


def test_breakouts_at_5_km_with_an_open_window(rows):
    """Spec test 6, for the model's own Soultz UCS and the measured range."""
    for r in rows:
        if r["source"] in ("site UCS (v1.0)", "Soultz lab UCS 100-130"):
            assert all(w > 0 for w in r["width"]), (r["source"], r["ucs"], r["cooling"])
        assert all(x["open"] for x in r["window"])


def test_only_the_strongest_intact_rock_escapes_breakout(rows):
    misses = [(r["source"], r["cooling"]) for r in rows for w in r["width"] if w == 0.0]
    assert misses and all(src == "UD-1 intact bound (C8)" and dK == 40.0
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
    onsets = [r["onset"][1] for r in rows if r["source"] == "UD-1 weak zones (C8)"]
    assert min(onsets) < cs.OBS_ONSET_FIRST < cs.OBS_ONSET_DENSE < max(onsets)


def test_strength_sources_are_the_documented_ones():
    src = dict(cs.STRENGTHS)
    assert src["Soultz lab UCS 100-130"] == (100e6, 115e6, 130e6)
    assert set(src["UD-1 weak zones (C8)"]) == {u for l, u in CORNWALL.strength_cases
                                               if "weak" in l}


def test_tensile_initiation_does_not_discriminate_with_t0_zero():
    """Observed tensile fractures are continuous only to 2,180 m TVD; with
    T0 = 0 (D4) the model initiates them at every depth checked."""
    zs = cs.tensile_depths(cs.SOULTZ.stress_cases[1], 0.0)
    assert zs[0] == 1500.0 and zs[-1] == 5000.0
