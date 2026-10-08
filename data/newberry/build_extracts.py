"""Rebuild the NWG 55-29 extracts in this folder from the raw GDR downloads.

    python data/newberry/build_extracts.py

Expects the raw files of GDR submission 271 in data/newberry/raw/gdr271/ (not
committed; see README.md):

  - Static_Profile_Oct2008_PT.csv   -> nwg5529_static_2008.csv
  - Deviation_Data.csv              -> nwg5529_deviation.csv
  - 55_29_RUN_3.las (open hole)     -> nwg5529_open_hole_logs.csv
"""
import csv
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw", "gdr271")
KB_ABOVE_GL_FT = 31.0      # RUN_3.las: "APD 31 ft, depth above permanent datum (GL)"


def _deviation():
    with open(os.path.join(RAW, "Deviation_Data.csv"), newline="") as fh:
        rows = [(float(r["MD_ft"]), float(r["TVD_ft"])) for r in csv.DictReader(fh)]
    return rows


def _tvd(md, dev):
    for (m0, t0), (m1, t1) in zip(dev, dev[1:]):
        if m0 <= md <= m1:
            return t0 + (t1 - t0) * (md - m0) / (m1 - m0)
    return dev[-1][1] + (md - dev[-1][0])


def deviation():
    dev = _deviation()
    with open(os.path.join(HERE, "nwg5529_deviation.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["md_ft", "tvd_ft"])
        w.writerows(dev)
    return len(dev)


def static_profile(step_ft=25.0):
    """The static pressure and temperature survey of October 2008, thinned to
    one point per step_ft of MD, with TVD from the deviation survey."""
    dev = _deviation()
    out, last = [], None
    with open(os.path.join(RAW, "Static_Profile_Oct2008_PT.csv"), newline="") as fh:
        rd = csv.reader(fh)
        next(rd)
        for r in rd:
            md, T, P = float(r[2]), float(r[1]), float(r[4])
            if last is None or md - last >= step_ft:
                out.append((md, round(_tvd(md, dev), 1), round(T, 2), P))
                last = md
    with open(os.path.join(HERE, "nwg5529_static_2008.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["md_ft", "tvd_ft", "temp_C", "press_psia"])
        w.writerows(out)
    return len(out)


def open_hole_logs(bin_ft=2.0):
    """Neutron porosity, bulk density, caliper and gamma ray from the open-hole
    run (Halliburton, 13 July 2008), as medians over bin_ft of depth, with
    depth moved from KB to ground level to match Davatzes and Hickman (2011)."""
    path = os.path.join(RAW, "55_29_RUN_3.las")
    with open(path, encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("~A"))
    names = lines[i].split()[1:]
    col = {n: names.index(n) for n in ("DEPT", "NPHI", "RHOB", "CALI", "GR")}
    bins = {}
    for l in lines[i + 1:]:
        f = l.split()
        if len(f) < len(names):
            continue
        md = float(f[col["DEPT"]]) - KB_ABOVE_GL_FT
        key = round(md / bin_ft) * bin_ft
        b = bins.setdefault(key, {k: [] for k in ("NPHI", "RHOB", "CALI", "GR")})
        for k in b:
            v = float(f[col[k]])
            if v != -999.25:
                b[k].append(v)
    with open(os.path.join(HERE, "nwg5529_open_hole_logs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["md_ft_gl", "nphi_pct", "rhob_gcc", "cali_in", "gr_api"])
        for key in sorted(bins):
            b = bins[key]
            if not b["NPHI"]:
                continue
            w.writerow([key] + [f"{statistics.median(b[k]):.3f}" if b[k] else ""
                                for k in ("NPHI", "RHOB", "CALI", "GR")])
    return len(bins)


if __name__ == "__main__":
    print("nwg5529_deviation.csv:", deviation(), "stations")
    print("nwg5529_static_2008.csv:", static_profile(), "points")
    print("nwg5529_open_hole_logs.csv:", open_hole_logs(), "bins")
