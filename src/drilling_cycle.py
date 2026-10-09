"""
drilling_cycle.py
=================
The drilling cycle a point on the borehole wall goes through, with timings
measured from the FORGE 16B record (data/forge16b/cycle_timings.csv and
trip_surface_operations.csv, built by data/forge16b/build_extracts.py).

The connections come from connections_10s.csv, the 10-second record over the
whole well: how long the pumps were off, and how long the fluid circulated
between the last new hole and the pumps stopping.

A wall element at depth z, one stand above the end of a bit run (the element
that has circulated least before the trip, and so reheats fastest), goes
through:

  1. drilled past, by the bit at the rate of penetration;
  2. circulating while the bit drills the last stand, with one connection pause;
  3. circulating for one bottoms-up before the trip;
  4. static for the trip: out at the measured speed, the routine time at
     surface, back in at the measured speed;
  5. circulating again.

Timings (D3): the medians of the trips' speeds out and in, with stops over
ten minutes excluded, scaled to the site's depth; the median and 90th
percentile connection, and the median and 10th percentile circulation after
the last new hole; the routine surface time (a bit and BHA change and
nothing else); bit runs of the FORGE length in hours. FORGE 16B is the only
drilling record in the repository, drilled in hot granite at 2.5 km TVD.
"""
import csv
import os
import statistics
import sys
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "forge16b")
FT = 0.3048


def _read(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        return list(csv.DictReader(fh))


def _pct(xs, q):
    xs = sorted(xs)
    return float(np.percentile(xs, q))


@dataclass(frozen=True)
class Timings:
    trip_out_ft_h: float
    trip_in_ft_h: float
    surface_h: float
    connection_median_min: float
    connection_p90_min: float
    fresh_circ_median_min: float     # circulation between the last new hole and the pumps stopping
    fresh_circ_p10_min: float
    stand_ft: float
    bit_run_h: float
    stage_spacing_ft: float
    stage_minutes: float
    trips_h: tuple            # bottom of the hole without circulation, each trip


def forge_timings():
    rows = _read("cycle_timings.csv")
    by = lambda k: [r for r in rows if r["kind"] == k]
    f = lambda r, c: float(r[c])
    conn10 = _read("connections_10s.csv")
    conn = [float(r["pumps_off_min"]) for r in conn10]
    fresh = [float(r["circulated_after_drilling_min"]) for r in conn10]
    # stand length: the hole drilled between connections where it was a whole
    # stand (connections are also taken for surveys and checks in between)
    drilled = [f(r, "drilled_since_ft") for r in by("connection")
               if r["drilled_since_ft"] and 80.0 <= f(r, "drilled_since_ft") <= 110.0]
    surf = [float(r["surface_hours"]) for r in _read("trip_surface_operations.csv")
            if r["routine"] == "yes"]
    # staged circulation on the way in, below 5,500 ft (where staging began)
    stages = [r for r in by("staged circulation") if f(r, "depth_ft") > 5500.0]
    spacing = []
    for a, b in zip(stages, stages[1:]):
        d = f(b, "depth_ft") - f(a, "depth_ft")
        if 100.0 < d < 1500.0:
            spacing.append(d)
    return Timings(
        trip_out_ft_h=statistics.median(f(r, "speed_ft_h") for r in by("trip out")),
        trip_in_ft_h=statistics.median(f(r, "speed_ft_h") for r in by("trip in")),
        surface_h=statistics.median(surf),
        connection_median_min=statistics.median(conn),
        connection_p90_min=_pct(conn, 90),
        fresh_circ_median_min=statistics.median(fresh),
        fresh_circ_p10_min=_pct(fresh, 10),
        stand_ft=statistics.median(drilled),
        bit_run_h=statistics.median(f(r, "hours") for r in by("bit run")),
        stage_spacing_ft=statistics.median(spacing),
        stage_minutes=statistics.median(60 * f(r, "hours") for r in stages),
        trips_h=tuple(f(r, "hours") for r in by("trip")),
    )


def trip_hours(md_m, t=None, scale=1.0):
    """Time the bottom of a hole at measured depth md_m spends without
    circulation for a trip: out, routine surface time, back in."""
    t = t or forge_timings()
    d = md_m / FT
    return scale * (d / t.trip_out_ft_h + t.surface_h + d / t.trip_in_ft_h)


@dataclass(frozen=True)
class Step:
    kind: str           # "circulating" or "static"
    hours: float
    label: str


def element_history(md_m, rop_m_h, bottoms_up_h, t=None, connection="median",
                    trip_scale=1.0):
    """The history of the wall element one stand above the end of a bit run at
    measured depth md_m, as steps of circulating and static time."""
    t = t or forge_timings()
    stand_h = t.stand_ft * FT / rop_m_h
    conn = (t.connection_median_min if connection == "median" else t.connection_p90_min) / 60.0
    return [
        Step("circulating", stand_h / 2, "drilling the last stand"),
        Step("static", conn, "connection"),
        Step("circulating", stand_h / 2, "drilling the last stand"),
        Step("circulating", bottoms_up_h, "bottoms-up"),
        Step("static", trip_hours(md_m, t, trip_scale), "trip"),
    ]


def report(t=None):
    t = t or forge_timings()
    print("=" * 72)
    print("Drilling-cycle timings, FORGE 16B: trips from the Pason 1-minute record, 21-29 May 2023;")
    print("connections from the 10-second record over the whole well")
    print("=" * 72)
    print("bottom of the hole without circulation, each trip [h]: "
          + ", ".join(f"{h:.1f}" for h in t.trips_h))
    print(f"tripping speed, stops over 10 min excluded: out {t.trip_out_ft_h:.0f} ft/h, "
          f"in {t.trip_in_ft_h:.0f} ft/h (medians)")
    print(f"routine surface time (bit and BHA change only): {t.surface_h:.1f} h")
    print(f"connection, pumps off (10-second record, whole well): median "
          f"{t.connection_median_min:.1f} min, 90th percentile {t.connection_p90_min:.1f} min")
    print(f"circulation between the last new hole and the pumps stopping: median "
          f"{60 * t.fresh_circ_median_min:.0f} s, 10th percentile {60 * t.fresh_circ_p10_min:.0f} s")
    print(f"stand length: {t.stand_ft:.0f} ft")
    print(f"bit run: median {t.bit_run_h:.1f} h on bottom")
    print(f"staged circulation on the way in: every {t.stage_spacing_ft:.0f} ft, "
          f"{t.stage_minutes:.0f} min each (medians)")
    print("-" * 72)
    for km in (2.5, 10.0, 12.9):
        print(f"a trip from {km:4.1f} km MD leaves the bottom without circulation for "
              f"{trip_hours(km * 1000, t):5.1f} h")
    print("=" * 72)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    report()
