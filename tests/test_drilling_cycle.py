"""Drilling-cycle timings measured from the FORGE 16B record, pinned."""
import pytest

import drilling_cycle as dc


@pytest.fixture(scope="module")
def t():
    return dc.forge_timings()


def test_trips_without_circulation_match_the_record(t):
    """The three trips between the trial runs left the bottom of the hole
    without circulation for about 24, 34 and 36 h."""
    assert t.trips_h[:3] == pytest.approx((24.0, 34.0, 36.0), abs=1.0)


def test_measured_timings_are_unchanged(t):
    assert t.trip_out_ft_h == pytest.approx(1694, abs=1)
    assert t.trip_in_ft_h == pytest.approx(1833, abs=1)
    assert t.surface_h == pytest.approx(3.93, abs=0.01)
    assert t.connection_median_min == pytest.approx(2.0, abs=0.05)
    assert t.connection_p90_min == pytest.approx(5.0, abs=0.05)
    assert t.stand_ft == pytest.approx(94, abs=1)
    assert t.bit_run_h == pytest.approx(11.4, abs=0.05)


def test_a_routine_trip_at_forge_depth_takes_as_long_as_forge_s_routine_trips(t):
    """The two trips with only a BHA change at surface took 13.5 and 15.0 h."""
    assert dc.trip_hours(9800 * dc.FT, t) == pytest.approx(14.2, abs=1.0)


def test_element_history(t):
    h = dc.element_history(10000.0, 5.0, 2.0, t)
    assert [s.kind for s in h] == ["circulating", "static", "circulating", "circulating", "static"]
    assert sum(s.hours for s in h if s.label == "drilling the last stand") == pytest.approx(
        t.stand_ft * dc.FT / 5.0)
    assert h[-1].hours == pytest.approx(dc.trip_hours(10000.0, t))
