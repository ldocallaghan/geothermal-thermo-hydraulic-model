"""The UT Austin (Wu et al. 2025) benchmark, pinned. These hold Model 1's
results fixed so that any change shows; they don't assert agreement with Wu."""
import json

import pytest

import benchmark_utaustin as b

pytestmark = pytest.mark.slow


def test_benchmark_table_is_unchanged():
    with open(b.GOLDEN) as fh:
        golden = json.load(fh)
    now = b.golden_rows()
    assert now.keys() == golden.keys()
    for key, row in golden.items():
        for col, val in row.items():
            assert now[key][col] == pytest.approx(val, abs=0.5), (key, col)


def test_pipe_ranking_matches_wu():
    """Coldest to hottest at the bit: dual-wall, internal, external, conventional."""
    for orient in ("vertical", "horizontal"):
        t = [b.bhct(b.run(b.PIPES[n], horizontal=orient == "horizontal"))
             for n in ("dual-wall", "internal", "external", "conventional")]
        assert t == sorted(t)
