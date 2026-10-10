"""Rebuild the FORGE 16B extracts in this folder from the raw GDR downloads.

    python data/forge16b/build_extracts.py

Expects the raw files in data/forge16b/raw/ (not committed; see README.md):

  - "16B(78)-32 Well Survey.zip"  (GDR 1516)        -> survey_16b.csv
  - "Utah FORGE Deep Well Temperature Profiles_Sept 2022.xlsx"
                                  (GDR 1421)        -> formation_temperature_16a.csv
  - "16B_Pason.zip", unzipped to raw/16B_Pason/
                                  (GDR 1516)        -> pason_trial_1min.csv,
                                                       drilling_history_10min.csv,
                                                       connections_10s.csv

run_observations.csv is transcribed by hand from the Eavor trial report, so it
is not rebuilt here. cycle_timings.csv is built from pason_trial_1min.csv.
"""
import csv
import io
import os
import re
import zipfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw")


def survey():
    """MD, inclination, azimuth and TVD from the gyro survey's text report."""
    z = zipfile.ZipFile(os.path.join(RAW, "16B(78)-32 Well Survey.zip"))
    name = next(n for n in z.namelist() if n.endswith("Final Survey Report.txt"))
    rows = []
    for line in io.TextIOWrapper(z.open(name), encoding="latin-1"):
        f = line.split()
        if len(f) >= 4 and re.fullmatch(r"-?\d+\.\d+", f[0] or ""):
            rows.append((float(f[0]), float(f[1]), float(f[2]), float(f[3])))
    with open(os.path.join(HERE, "survey_16b.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["md_ft", "inc_deg", "azi_deg", "tvd_ft"])
        w.writerows(rows)
    return len(rows)


_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _sheet(path, sheet_name):
    """Rows of one worksheet as lists of strings (a minimal xlsx reader)."""
    z = zipfile.ZipFile(path)
    shared = [("".join(t.text or "" for t in si.iter("{%s}t" % _NS["m"])))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS)]
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target")
            for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    sheet = next(s for s in wb.find("m:sheets", _NS) if s.get("name") == sheet_name)
    part = "xl/" + rels[sheet.get("{%s}id" % _NS["r"])].replace("xl/", "").lstrip("/")
    for row in ET.fromstring(z.read(part)).iter("{%s}row" % _NS["m"]):
        cells = []
        for c in row.findall("m:c", _NS):
            v = c.find("m:v", _NS)
            val = None if v is None else (shared[int(v.text)] if c.get("t") == "s" else v.text)
            col = 0
            for ch in re.match(r"[A-Z]+", c.get("r")).group():
                col = col * 26 + ord(ch) - 64
            while len(cells) < col - 1:
                cells.append(None)
            cells.append(val)
        yield cells


def formation_temperature(step_ft=25.0):
    """16A(78)-32 cement-bond-log temperature of 16 Aug 2021 (~216 days of
    recovery), against TVD, thinned to one point per step_ft of TVD."""
    path = os.path.join(RAW, "Utah FORGE Deep Well Temperature Profiles_Sept 2022.xlsx")
    out, last = [], None
    for r in _sheet(path, "16A(78)-32_2nd"):
        try:
            md, tvd, t = float(r[0]), float(r[1]), float(r[2])
        except (TypeError, ValueError, IndexError):
            continue
        if last is None or tvd - last >= step_ft:
            out.append((md, round(tvd, 1), t))
            last = tvd
    with open(os.path.join(HERE, "formation_temperature_16a.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["md_ft", "tvd_ft", "temp_F"])
        w.writerows(out)
    return len(out)


PASON_COLS = {
    "Hole Depth (feet)": "hole_depth_ft",
    "Bit Depth (feet)": "bit_depth_ft",
    "Standpipe Pressure (psi)": "spp_psi",
    "Total Pump Output (gal_per_min)": "flow_gpm",
    "TEMP IN MANIFOLD (DEGREES)": "temp_in_F",
    "TEMP OUT FLOW (DEGREES)": "temp_out_F",
    "Rate Of Penetration (ft_per_hr)": "rop_ft_h",
    "Differential Pressure (psi)": "diff_pressure_psi",
}


def pason_trial(start="2023/05/21", end="2023/05/29", minutes=1):
    """The Pason 10-second record over the trial (BHA 10 to 13), averaged over
    each minute. -999.25 is Pason's null and is dropped before averaging."""
    path = os.path.join(RAW, "16B_Pason", "10 Second Data.csv")
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        idx = {PASON_COLS[h]: i for i, h in enumerate(head) if h in PASON_COLS}
        names = list(PASON_COLS.values())
        buckets, order = {}, []
        for r in rd:
            if r[0] < start:
                continue
            if r[0] > end:
                break
            hh, mm, _ = r[1].split(":")
            key = f"{r[0].replace('/', '-')} {hh}:{int(mm) // minutes * minutes:02d}"
            if key not in buckets:
                buckets[key] = {n: [] for n in names}
                order.append(key)
            for n in names:
                try:
                    v = float(r[idx[n]])
                except ValueError:
                    continue
                if v != -999.25:
                    buckets[key][n].append(v)
    with open(os.path.join(HERE, "pason_trial_1min.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time"] + names)
        for k in order:
            b = buckets[k]
            w.writerow([k] + [f"{sum(b[n]) / len(b[n]):.2f}" if b[n] else "" for n in names])
    return len(order)


def drilling_history(minutes=10, circulating_gpm=100.0):
    """The whole Pason record in 10-minute bins: deepest hole depth reached and
    the fraction of the bin spent circulating (pump output above 100 gal/min),
    for the wall's exposure to circulation."""
    path = os.path.join(RAW, "16B_Pason", "10 Second Data.csv")
    out, key0, hmax, n, nc = [], None, 0.0, 0, 0
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        ih, iq = head.index("Hole Depth (feet)"), head.index("Total Pump Output (gal_per_min)")
        for r in rd:
            hh, mm, _ = r[1].split(":")
            key = f"{r[0].replace('/', '-')} {hh}:{int(mm) // minutes * minutes:02d}"
            if key != key0 and key0 is not None:
                out.append((key0, f"{hmax:.2f}", f"{nc / n:.3f}"))
                n = nc = 0
            key0 = key
            try:
                h, q = float(r[ih]), float(r[iq])
            except ValueError:
                continue
            if h != -999.25:
                hmax = max(hmax, h)
            n += 1
            nc += q != -999.25 and q > circulating_gpm
        out.append((key0, f"{hmax:.2f}", f"{nc / max(n, 1):.3f}"))
    with open(os.path.join(HERE, "drilling_history_10min.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time", "max_hole_depth_ft", "circulating_fraction"])
        w.writerows(out)
    return len(out)


def _pason_minutes():
    from datetime import datetime
    with open(os.path.join(HERE, "pason_trial_1min.csv"), newline="") as fh:
        rows = list(csv.DictReader(fh))
    num = lambda r, k: float(r[k]) if r[k] else float("nan")
    return [dict(t=datetime.strptime(r["time"], "%Y-%m-%d %H:%M"),
                 hole=num(r, "hole_depth_ft"), bit=num(r, "bit_depth_ft"),
                 flow=num(r, "flow_gpm")) for r in rows]


def _runs(flags):
    """(start, end) index pairs of consecutive True values."""
    out, i, n = [], 0, len(flags)
    while i < n:
        if flags[i]:
            j = i
            while j + 1 < n and flags[j + 1]:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


ON_BOTTOM_FT = 100.0      # bit within this of the hole bottom
PUMPS_OFF_GPM = 100.0     # pump output below this is no circulation
FULL_FLOW_GPM = 500.0
# The Pason bit depth stops updating while the BHA and heavy-weight pipe
# (about 1,300 ft together) are handled at surface; a bit shallower than this
# is taken as out of the hole.
SURFACE_FT = 1000.0
STOP_MIN = 10             # a stationary spell longer than this is a stop on a trip


def cycle_timings():
    """Drilling-cycle timings of FORGE 16B, 21-29 May 2023, from the 1-minute
    Pason extract: every gap in circulation with the bit on bottom (trips and
    pauses), the trips' speeds out and in with stops excluded, the time at
    surface, staged-circulation stops on the way in, connections while
    drilling, and the stand length (hole drilled between connections)."""
    m = _pason_minutes()
    n = len(m)
    hrs = lambda i, j: (m[j]["t"] - m[i]["t"]).total_seconds() / 3600.0
    on_bottom_circ = [m[i]["flow"] >= PUMPS_OFF_GPM and m[i]["hole"] - m[i]["bit"] <= ON_BOTTOM_FT
                      for i in range(n)]
    rows = []

    # gaps in circulation at the bottom of the hole
    idx = [i for i in range(n) if on_bottom_circ[i]]
    for i, j in zip(idx, idx[1:]):
        if j - i < 20:
            continue
        bits = [m[k]["bit"] for k in range(i, j + 1)]
        tripped = min(bits) < SURFACE_FT
        rows.append(dict(kind="trip" if tripped else "pause on bottom",
                         start=m[i]["t"], end=m[j]["t"], hours=hrs(i, j),
                         depth_ft=m[i]["hole"], min_bit_ft=min(bits)))

    # bit runs: on bottom between the end of one trip and the start of the next
    # (only spans in which hole was drilled; the first starts when the bit
    # first reached bottom in the record)
    trips = sorted([r for r in rows if r["kind"] == "trip"], key=lambda r: r["start"])
    hole_at = {m[i]["t"]: m[i]["hole"] for i in range(n)}
    starts = [m[idx[0]]["t"]] + [r["end"] for r in trips[:-1]]
    for t0, r1 in zip(starts, trips):
        drilled = hole_at[r1["start"]] - hole_at[t0]
        if drilled > 50.0:
            rows.append(dict(kind="bit run", start=t0, end=r1["start"],
                             hours=(r1["start"] - t0).total_seconds() / 3600.0,
                             depth_ft=r1["depth_ft"], drilled_since_ft=drilled))

    # trips out and in: from leaving bottom to surface, and from surface to bottom
    for r in [r for r in rows if r["kind"] == "trip"]:
        i0 = next(k for k in range(n) if m[k]["t"] == r["start"])
        i1 = next(k for k in range(n) if m[k]["t"] == r["end"])
        k_surf = next(k for k in range(i0, i1) if m[k]["bit"] < SURFACE_FT)
        k_leave = max(k for k in range(i0, i1) if m[k]["bit"] < SURFACE_FT)
        for kind, a, b in (("trip out", i0, k_surf), ("trip in", k_leave, i1)):
            moving = 0
            still = 0
            for k in range(a + 1, b + 1):
                if abs(m[k]["bit"] - m[k - 1]["bit"]) < 5.0:
                    still += 1
                else:
                    moving += 1 + (still if still <= STOP_MIN else 0)
                    still = 0
            dist = abs(m[b]["bit"] - m[a]["bit"])
            rows.append(dict(kind=kind, start=m[a]["t"], end=m[b]["t"], hours=hrs(a, b),
                             depth_ft=r["depth_ft"], speed_ft_h=dist / (moving / 60.0),
                             moving_hours=moving / 60.0))
        rows.append(dict(kind="at surface", start=m[k_surf]["t"], end=m[k_leave]["t"],
                         hours=hrs(k_surf, k_leave), depth_ft=r["depth_ft"]))
        # staged circulation on the way in: circulating with the bit off bottom
        circ = [m[k]["flow"] >= PUMPS_OFF_GPM and m[k]["hole"] - m[k]["bit"] > ON_BOTTOM_FT
                for k in range(k_leave, i1)]
        for a, b in _runs(circ):
            if b - a + 1 >= 5:
                rows.append(dict(kind="staged circulation", start=m[k_leave + a]["t"],
                                 end=m[k_leave + b]["t"], hours=(b - a + 1) / 60.0,
                                 depth_ft=m[k_leave + a]["bit"]))

    # connections while drilling: pumps off for under an hour with the bit on
    # bottom on both sides, during a run (hole depth rising)
    off = [m[i]["flow"] < PUMPS_OFF_GPM for i in range(n)]
    last_conn_hole = None
    for a, b in _runs(off):
        if a == 0 or b == n - 1 or b - a + 1 > 60:
            continue
        before, after = m[a - 1], m[b + 1]
        if (before["hole"] - before["bit"] <= ON_BOTTOM_FT and after["hole"] - after["bit"] <= ON_BOTTOM_FT
                and before["flow"] >= FULL_FLOW_GPM):
            drilled = None if last_conn_hole is None else before["hole"] - last_conn_hole
            rows.append(dict(kind="connection", start=m[a]["t"], end=m[b]["t"],
                             hours=(b - a + 1) / 60.0, depth_ft=before["hole"],
                             drilled_since_ft=drilled))
            last_conn_hole = before["hole"]

    cols = ["kind", "start", "end", "hours", "depth_ft", "min_bit_ft", "speed_ft_h",
            "moving_hours", "drilled_since_ft"]
    with open(os.path.join(HERE, "cycle_timings.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in sorted(rows, key=lambda r: r["start"]):
            w.writerow([r.get(c) if not isinstance(r.get(c), float) else f"{r[c]:.3f}"
                        for c in cols])
    return len(rows)


def connections_10s(full_gpm=500.0):
    """Connections while drilling over the whole well, from the 10-second
    record: the pumps off (below 100 gal/min) for under an hour, at full flow
    in the minute before, with the bit within 100 ft of bottom on both sides
    and new hole made in the half hour before. For each, how long the pumps
    were off, and how long the fluid circulated between the last new hole and
    the pumps stopping (the cooling the freshly drilled rock gets before the
    connection)."""
    from datetime import datetime, timezone
    path = os.path.join(RAW, "16B_Pason", "10 Second Data.csv")
    t, hole, bit, q = [], [], [], []
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        ih, ib = head.index("Hole Depth (feet)"), head.index("Bit Depth (feet)")
        iq = head.index("Total Pump Output (gal_per_min)")
        for r in rd:
            try:
                h, b, f = float(r[ih]), float(r[ib]), float(r[iq])
            except ValueError:
                continue
            if -999.25 in (h, b, f):
                continue
            t.append(datetime.strptime(r[0] + " " + r[1], "%Y/%m/%d %H:%M:%S").replace(tzinfo=timezone.utc))
            hole.append(h); bit.append(b); q.append(f)
    ts = [x.timestamp() for x in t]
    n = len(t)
    off = [f < PUMPS_OFF_GPM for f in q]
    rows = []
    for a, b in _runs(off):
        if a == 0 or b == n - 1:
            continue
        dur = ts[b] - ts[a] + 10.0
        before = a - 1
        if (dur > 3600.0 or max(q[max(0, before - 6):before + 1]) < full_gpm
                or hole[before] - bit[before] > ON_BOTTOM_FT
                or hole[b + 1] - bit[b + 1] > ON_BOTTOM_FT):
            continue
        k = before
        while k > 0 and ts[before] - ts[k] < 1800.0 and not hole[k] > hole[k - 1] + 0.01:
            k -= 1
        if not (hole[k] > hole[k - 1] + 0.01 and ts[before] - ts[k] < 1800.0):
            continue
        rows.append((t[a].strftime("%Y-%m-%d %H:%M:%S"), f"{hole[before]:.1f}",
                     f"{dur / 60.0:.2f}", f"{(ts[a] - ts[k]) / 60.0:.2f}"))
    with open(os.path.join(HERE, "connections_10s.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["start", "hole_depth_ft", "pumps_off_min", "circulated_after_drilling_min"])
        w.writerows(rows)
    return len(rows)


if __name__ == "__main__":
    print("survey_16b.csv:", survey(), "stations")
    print("formation_temperature_16a.csv:", formation_temperature(), "points")
    print("pason_trial_1min.csv:", pason_trial(), "minutes")
    print("drilling_history_10min.csv:", drilling_history(), "bins")
    print("cycle_timings.csv:", cycle_timings(), "periods")
    print("connections_10s.csv:", connections_10s(), "connections")
