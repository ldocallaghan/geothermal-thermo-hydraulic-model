"""Rebuild the FORGE 16B extracts in this folder from the raw GDR downloads.

    python data/forge16b/build_extracts.py

Expects the raw files in data/forge16b/raw/ (not committed; see README.md):

  - "16B(78)-32 Well Survey.zip"  (GDR 1516)        -> survey_16b.csv
  - "Utah FORGE Deep Well Temperature Profiles_Sept 2022.xlsx"
                                  (GDR 1421)        -> formation_temperature_16a.csv
  - "16B_Pason.zip", unzipped to raw/16B_Pason/
                                  (GDR 1516)        -> pason_trial_1min.csv,
                                                       drilling_history_10min.csv

run_observations.csv is transcribed by hand from the Eavor trial report, so it
is not rebuilt here.
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


if __name__ == "__main__":
    print("survey_16b.csv:", survey(), "stations")
    print("formation_temperature_16a.csv:", formation_temperature(), "points")
    print("pason_trial_1min.csv:", pason_trial(), "minutes")
    print("drilling_history_10min.csv:", drilling_history(), "bins")
