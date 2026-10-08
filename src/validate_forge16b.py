"""
validate_forge16b.py
====================
Validation of Model 1 against the insulated drill pipe (IDP) trial in Utah
FORGE well 16B(78)-32, May 2023. Data and their sources are in data/forge16b/
(README.md there).

The trial ran four bits through 9-1/2-inch hole at 60-70 deg inclination:

  BHA 10  7,584-8,085 ft MD  regular drill pipe            MWD average 180 F
  BHA 11  8,085-8,585 ft     full IDP string                MWD average 149 F
  BHA 12  8,585-9,255 ft     ~70% IDP above regular pipe    MWD average 164 F,
                                                            levelling at 150-160 F
  BHA 13  9,255-9,800 ft     regular drill pipe            MWD average 220 F;
                                                            steady 211 F at 700 gal/min,
                                                            224 F at 600 gal/min

The model is the steady counterflow model, run quasi-steadily: at regular times
while each bit was on bottom, it is solved with that minute's bit depth, flow
and inlet temperature (Pason), the wall's exposure to circulation from the
drilling record, and the string in the hole at the time. The MWD temperature
is compared with the fluid temperature in the bore at the sensor, 84 ft above
the bit, because the tool sits in the downgoing flow; the annulus temperature
there is reported beside it.

Protocol:
  1. Fit the annulus film multiplier on the runs without IDP (BHA 10, and BHA 13
     at 600 and 700 gal/min), with the rock's properties at their measured
     values. If the unfitted model is within the data's scatter, fit nothing.
  2. With that fixed, fit one effective conductivity for Eavor's IDP on BHA 11.
     It lumps coating, steel and joints together; their specification is not
     published, so the IDP is modelled as one wall layer across the 5-1/2-inch
     pipe body.
  3. Predict BHA 12 with nothing refitted.
  4. Rebuild Eavor's breakdown of the BHA 10 to BHA 11 change, one factor at a time.
"""
import csv
import os
import sys
from datetime import datetime

import numpy as np
from scipy.optimize import brentq, minimize_scalar

import model1_coupled as m1
import well_geometry as wg

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "forge16b")
FT, INCH, GPM, PSI = 0.3048, 0.0254, 6.30902e-5, 6894.757


def F_to_C(f):
    return (np.asarray(f, float) - 32.0) / 1.8


def C_to_F(c):
    return np.asarray(c, float) * 1.8 + 32.0


def _read(name):
    with open(os.path.join(DATA, name), newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- rock
# Means of three FORGE granitoid cores at 100 C (MetaRock 2021, tabulated by
# Xing et al. 2026; data/forge16b/README.md). The range across samples and
# temperatures is K_ROCK_RANGE.
K_ROCK = (2.149 + 2.341 + 3.002) / 3
RHOCP_ROCK = (2680.0 + 2600.0 + 2650.0) / 3 * (864.8 + 862.3 + 1091.1) / 3
K_ROCK_RANGE = (2.041, 3.087)


# ---------------------------------------------------------------- survey and temperature
def survey():
    rows = _read("survey_16b.csv")
    md = np.array([float(r["md_ft"]) for r in rows]) * FT
    tvd = np.array([float(r["tvd_ft"]) for r in rows]) * FT
    if md[0] > 0:
        md, tvd = np.concatenate([[0.0], md]), np.concatenate([[0.0], tvd])
    return wg.Survey(md, tvd)


def formation_temperature():
    """Static temperature of 16A(78)-32 against TVD [m] -> [C]."""
    rows = _read("formation_temperature_16a.csv")
    tvd = np.array([float(r["tvd_ft"]) for r in rows]) * FT
    T = F_to_C([float(r["temp_F"]) for r in rows])
    return lambda d: np.interp(d, tvd, T)


# ---------------------------------------------------------------- hole
def hole(md_bit_m):
    """16-inch surface casing to 1,136 ft in 22-inch hole; 11-3/4-inch 65 lb/ft
    intermediate casing to 4,837 ft in 14-3/4-inch hole; 9-1/2-inch open hole
    below (End of Well Report). The 11-3/4-inch casing's ID is the API nominal
    for 65 lb/ft. The 16-inch casing's weight isn't given; a 1/2-inch wall is
    assumed (its resistance is negligible next to the cement and rock)."""
    st, ce = (lambda a, b: wg.Layer(a, b, wg.K_STEEL)), (lambda a, b: wg.Layer(a, b, wg.K_CEMENT))
    r = lambda d_in: d_in * INCH / 2
    i_id, i_od, h2 = r(10.682), r(11.75), r(14.75)
    s_id, s_od, h1 = r(15.0), r(16.0), r(22.0)
    return [
        wg.HoleInterval(0.0, 1136 * FT, i_id,
                        (st(i_id, i_od), ce(i_od, s_id), st(s_id, s_od), ce(s_od, h1))),
        wg.HoleInterval(1136 * FT, 4837 * FT, i_id, (st(i_id, i_od), ce(i_od, h2))),
        wg.open_hole(4837 * FT, max(md_bit_m, 4837 * FT + 1.0), 9.5 * INCH),
    ]


# ---------------------------------------------------------------- string
def _steel(name, od_in, id_in):
    return wg.PipeType(name, (wg.Layer(id_in * INCH / 2, od_in * INCH / 2, wg.K_STEEL),))


def _bha_dimensions():
    """Length-weighted OD and ID of the BHA from motor to filter sub (EOW Fig. 90)."""
    rows = [r for r in _read("string_components.csv") if r["bha"] == "17" and 2 <= int(r["item"]) <= 10]
    L = np.array([float(r["length_ft"]) for r in rows])
    od = np.array([float(r["od_in"]) for r in rows])
    idd = np.array([float(r["id_in"]) for r in rows])
    return float(np.sum(L * od) / L.sum()), float(np.sum(L * idd) / L.sum())


BHA_OD, BHA_ID = _bha_dimensions()
BHA = _steel("BHA", BHA_OD, BHA_ID)
COLLAR = _steel("drill collars", 6.813, 2.875)
HWDP = _steel("heavy-weight drill pipe", 5.5, 3.625)
HWDP_JOINT_FT = 913.42 / 30            # EOW Fig. 90
DRILL_PIPE = _steel("5-1/2 in 24.7 lb/ft drill pipe", 5.5, 4.670)
MWD_ABOVE_BIT_FT = 84.0                # EOW Fig. 102, survey sensor offset


def idp(k_eff):
    """Eavor's IDP as one effective wall on the 5-1/2-inch body."""
    return wg.PipeType(f"IDP, k_eff {k_eff:.3g}",
                       (wg.Layer(4.670 * INCH / 2, 5.5 * INCH / 2, k_eff),))


RUNS = {r["run"]: r for r in _read("run_observations.csv") if r["run"] != "BHA13"}
RUNS["BHA13"] = next(r for r in _read("run_observations.csv") if r["run"] == "BHA13")
DEPTHS = {k: (float(r["depth_in_ft_md"] or 7584), float(r["depth_out_ft_md"])) for k, r in RUNS.items()}
DEPTHS["BHA10"] = (7584.0, 8085.0)     # EOW bit record; the trial report gives only depth out
BHA12_REGULAR_FT = 2378.5
NOZZLES = {"BHA10": (14,) * 8, "BHA11": (14,) * 8, "BHA12": (14, 14, 14, 16, 16, 16),
           "BHA13": (14,) * 8}         # rig_hydraulics.csv


def tfa(run):
    return float(sum(np.pi / 4 * (n / 32 * INCH) ** 2 for n in NOZZLES[run]))


def string(run, bit_ft, k_idp=None):
    r = RUNS[run]
    bha_ft = float(r["bha_nmdc_ft"])
    dc_ft = float(r["drill_collars_ft"] or 0)
    hw_ft = int(r["hwdp_joints"]) * HWDP_JOINT_FT
    dp_ft = bit_ft - bha_ft - dc_ft - hw_ft
    segs = []
    if run == "BHA11":
        segs.append((dp_ft, idp(k_idp)))
    elif run == "BHA12":
        segs += [(dp_ft - BHA12_REGULAR_FT, idp(k_idp)), (BHA12_REGULAR_FT, DRILL_PIPE)]
    else:
        segs.append((dp_ft, DRILL_PIPE))
    segs.append((hw_ft, HWDP))
    if dc_ft:
        segs.append((dc_ft, COLLAR))
    segs.append((bha_ft, BHA))
    return [wg.StringSegment(L * FT, p) for L, p in segs]


def geometry(run, bit_ft, k_idp=None):
    return wg.WellGeometry(string(run, bit_ft, k_idp), hole(bit_ft * FT), SURVEY)


# ---------------------------------------------------------------- drilling record
def _time(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M")


def circulation_clock():
    """Cumulative circulating time [s] against wall-clock time, and the deepest
    hole [m] reached by then, from the whole Pason record."""
    rows = _read("drilling_history_10min.csv")
    t = np.array([_time(r["time"]).timestamp() for r in rows])
    frac = np.array([float(r["circulating_fraction"]) for r in rows])
    depth = np.array([float(r["max_hole_depth_ft"]) for r in rows]) * FT
    return t, np.cumsum(frac * 600.0), depth


def exposure(t_now):
    """The wall at each depth has been exposed for the circulating time since
    the bit first reached it. Time without circulation is not counted; the
    rock's recovery during it is not modelled either."""
    t, circ, depth = CLOCK
    c_now = float(np.interp(t_now, t, circ))
    return m1.exposure_from_history(circ, depth, c_now)


def on_bottom_samples(run, every_min=30):
    """Times when the bit was drilling in this run, every `every_min` minutes."""
    lo, hi = DEPTHS[run]
    out = []
    for r in PASON:
        try:
            h, b, q = float(r["hole_depth_ft"]), float(r["bit_depth_ft"]), float(r["flow_gpm"])
            tin = float(r["temp_in_F"])
        except ValueError:
            continue
        rop = float(r["rop_ft_h"] or 0)
        if lo + 1 < h <= hi + 0.5 and b > h - 5 and q > 300 and rop > 1:
            dp = float(r["diff_pressure_psi"]) if r["diff_pressure_psi"] else 0.0
            out.append(dict(time=_time(r["time"]).timestamp(), bit_ft=b, flow_gpm=q,
                            T_in_F=tin, motor_dp_psi=max(dp, 0.0)))
    keep, last = [], -1e18
    for s in out:
        if s["time"] - last >= every_min * 60:
            keep.append(s)
            last = s["time"]
    return keep


# ---------------------------------------------------------------- model
def model_at(run, s, film_mult=1.0, k_idp=None, k_rock=K_ROCK, **override):
    """Model 1 at one on-bottom sample; override any of bit_ft, flow_gpm,
    T_in_F, time, motor_dp_psi to move one factor at a time."""
    s = {**s, **override}
    g = geometry(run, s["bit_ft"], k_idp)
    T_in = float(F_to_C(s["T_in_F"]))
    rho_in = float(m1.fluid_props(np.array([T_in]), np.array([m1.P_SURF]))[0][0])
    Q = s["flow_gpm"] * GPM
    exp = m1.smooth_exposure(exposure(s["time"]), g.md_bit)
    kw = dict(m_dot=Q * rho_in, T_inj=T_in, geotherm=T_FORM, geometry=g,
              target_rock_T=float(T_FORM(SURVEY.tvd(g.md_bit))),
              Q_face=s["motor_dp_psi"] * PSI * Q, k_rock=k_rock, rhocp_rock=RHOCP_ROCK,
              friction_heat=True, tfa=tfa(run), film_mult=film_mult, verbose=False)
    # seed with a uniform one-day exposure, which always converges
    r = m1.solve(**kw, t_years=1 / 365.0)
    r = m1.solve(**kw, exposure=exp, breaks=m1.exposure_breaks(exp, g.md_bit), init=r)
    z_mwd = g.md_bit - MWD_ABOVE_BIT_FT * FT
    return dict(mwd_F=float(C_to_F(np.interp(z_mwd, r["z"], r["Td"]))),
                ann_F=float(C_to_F(np.interp(z_mwd, r["z"], r["Tu"]))),
                ok=bool(r["success"]), **s)


def run_series(run, film_mult=1.0, k_idp=None, every_min=30, k_rock=K_ROCK):
    return [model_at(run, s, film_mult, k_idp, k_rock) for s in on_bottom_samples(run, every_min)]


def mean_mwd(series, flow=None):
    """Mean modelled MWD temperature over converged samples (optionally those
    within 25 gal/min of `flow`)."""
    xs = [p["mwd_F"] for p in series
          if p["ok"] and (flow is None or abs(p["flow_gpm"] - flow) < 25)]
    return float(np.mean(xs)) if xs else float("nan")


MEASURED = {"BHA10": 180.0, "BHA11": 149.0, "BHA12": 164.0, "BHA13": 220.0,
            "BHA12 levelled": (150.0, 160.0), "BHA13 700": 211.0, "BHA13 600": 224.0}


def step1_targets(film_mult, every_min=60):
    s10, s13 = run_series("BHA10", film_mult, every_min=every_min), \
        run_series("BHA13", film_mult, every_min=every_min)
    return {"BHA10": mean_mwd(s10), "BHA13 700": mean_mwd(s13, 700),
            "BHA13 600": mean_mwd(s13, 600)}


SURVEY = survey()
T_FORM = formation_temperature()
CLOCK = circulation_clock()
PASON = _read("pason_trial_1min.csv")


# ---------------------------------------------------------------- protocol
def fit_film(every_min=60):
    """Step 1: the annulus film multiplier that best matches BHA 10 and BHA 13
    (600 and 700 gal/min), by least squares. Its uncertainty is the range of
    multipliers that match each of the three on its own."""
    keys = ("BHA10", "BHA13 700", "BHA13 600")
    cache = {}

    def model(f):
        if f not in cache:
            cache[f] = step1_targets(f, every_min)
        return cache[f]

    def sse(f):
        m = model(f)
        return sum((m[k] - MEASURED[k]) ** 2 for k in keys)

    best = minimize_scalar(sse, bounds=(0.2, 3.0), method="bounded",
                           options=dict(xatol=0.01)).x
    singles = []
    for k in keys:
        g = lambda f: model(f)[k] - MEASURED[k]
        try:
            singles.append(brentq(g, 0.2, 3.0, xtol=0.01))
        except ValueError:
            singles.append(float("nan"))
    return dict(unfitted=model(1.0), film_mult=float(best), fitted=model(best),
                single_fits=dict(zip(keys, singles)))


def fit_idp(film_mult, target=MEASURED["BHA11"], every_min=30):
    """Step 2: the IDP's effective wall conductivity that reproduces BHA 11's
    average MWD temperature."""
    g = lambda lk: mean_mwd(run_series("BHA11", film_mult, 10 ** lk, every_min)) - target
    return 10 ** brentq(g, np.log10(0.01), np.log10(wg.K_STEEL), xtol=0.005)


def decomposition(film_mult, k_idp):
    """Step 4: Eavor's breakdown of the BHA 10 to BHA 11 change, rebuilt one
    factor at a time from the middle of BHA 10 (flows and inlets are the run
    averages in the trial report)."""
    s10 = on_bottom_samples("BHA10", 30)
    s11 = on_bottom_samples("BHA11", 30)
    b10, b11 = s10[len(s10) // 2], s11[len(s11) // 2]
    base = dict(flow_gpm=600.0, T_in_F=98.0)
    m = lambda run, s, k=None, **o: model_at(run, s, film_mult, k, **{**base, **o})["mwd_F"]
    t0 = m("BHA10", b10)
    rows = [
        ("flow 600 to 700 gal/min", m("BHA10", b10, flow_gpm=700.0) - t0, "-12 to -14"),
        ("inlet 98 to 121 F (run averages)", m("BHA10", b10, T_in_F=121.0) - t0, "+23"),
        ("inlet 130 to 110 F (chiller on)",
         m("BHA10", b10, T_in_F=110.0) - m("BHA10", b10, T_in_F=130.0), "-20 (1:1 rule)"),
        ("depth, mid BHA 10 to mid BHA 11", m("BHA11", b11, wg.K_STEEL) - t0, "+7"),
        ("IDP on, at mid BHA 11",
         m("BHA11", b11, k_idp) - m("BHA11", b11, wg.K_STEEL), "-47 to -49"),
    ]
    return t0, rows


def run_all(every_min=30):
    out = {}
    step1 = fit_film()
    f = step1["film_mult"]
    k_idp = fit_idp(f)
    lo, hi = np.nanmin(list(step1["single_fits"].values())), \
        np.nanmax(list(step1["single_fits"].values()))
    k_range = (fit_idp(lo), fit_idp(hi)) if np.isfinite(lo) and np.isfinite(hi) else (np.nan, np.nan)
    series = {run: run_series(run, f, k_idp, every_min) for run in ("BHA10", "BHA11", "BHA12", "BHA13")}
    s12 = series["BHA12"]
    late = [p for p in s12 if p["bit_ft"] > 8700]
    out.update(step1=step1, k_idp=k_idp, k_idp_range=k_range, series=series,
               bha12_mean=mean_mwd(s12), bha12_late=mean_mwd(late),
               decomposition=decomposition(f, k_idp))
    # rock conductivity sensitivity on the fitted model
    out["rock_sensitivity"] = {
        k: {run: mean_mwd(run_series(run, f, k_idp, 60, k_rock=k)) for run in ("BHA10", "BHA11", "BHA12")}
        for k in K_ROCK_RANGE}
    return out


def report(o):
    s1 = o["step1"]
    print("=" * 78)
    print("Model 1 against the FORGE 16B IDP trial: MWD temperature [F]")
    print("=" * 78)
    print("\n1. Rock and annulus side, runs without IDP (rock k %.2f W/m K, measured)" % K_ROCK)
    print(f"{'case':>12}{'measured':>10}{'unfitted':>10}{'fitted':>9}")
    for k in ("BHA10", "BHA13 700", "BHA13 600"):
        print(f"{k:>12}{MEASURED[k]:10.0f}{s1['unfitted'][k]:10.1f}{s1['fitted'][k]:9.1f}")
    sf = s1["single_fits"]
    print(f"  annulus film multiplier: {s1['film_mult']:.2f}  (fitting each case alone: "
          + ", ".join(f"{k} {v:.2f}" for k, v in sf.items()) + ")")
    print(f"\n2. IDP effective wall conductivity, fitted on BHA 11 (149 F): "
          f"{o['k_idp']:.3f} W/m K  (with the film multiplier at the ends of its range: "
          f"{o['k_idp_range'][0]:.3f} to {o['k_idp_range'][1]:.3f})")
    print("\n3. All runs with the fitted values; BHA 12 is a blind prediction")
    print(f"{'run':>8}{'measured':>16}{'model':>8}{'error':>7}{'annulus':>9}{'samples':>9}")
    rows = [("BHA10", 180.0), ("BHA11", 149.0), ("BHA12", 164.0), ("BHA13", 220.0)]
    for run, meas in rows:
        ser = o["series"][run]
        mm = mean_mwd(ser)
        ann = float(np.mean([p["ann_F"] for p in ser if p["ok"]]))
        tag = " (blind)" if run == "BHA12" else ""
        print(f"{run:>8}{meas:16.0f}{mm:8.1f}{mm-meas:+7.1f}{ann:9.1f}{len(ser):9d}{tag}")
    print(f"{'BHA12':>8}{'150-160 (late)':>16}{o['bha12_late']:8.1f}"
          f"{o['bha12_late']-155:+7.1f}{'':9}{'':9} (blind, below 8,700 ft)")
    t0, rows = o["decomposition"]
    print(f"\n4. BHA 10 to BHA 11, one factor at a time (from {t0:.0f} F at mid BHA 10)")
    print(f"{'factor':>36}{'Model 1':>9}{'Eavor':>13}")
    for name, d, eav in rows:
        print(f"{name:>36}{d:+9.1f}{eav:>13}")
    print("\nRock conductivity sensitivity (fitted model, run means)")
    for k, v in o["rock_sensitivity"].items():
        print(f"  k = {k:.2f}: " + ", ".join(f"{r} {x:.0f}" for r, x in v.items()))


def figure(o, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    colours = {"BHA10": "k", "BHA11": "tab:blue", "BHA12": "tab:green", "BHA13": "tab:red"}
    meas = {"BHA10": 180, "BHA11": 149, "BHA12": 164, "BHA13": 220}
    for run, ser in o["series"].items():
        d = [p["bit_ft"] for p in ser if p["ok"]]
        t = [p["mwd_F"] for p in ser if p["ok"]]
        ax[0].plot(d, t, "o-", color=colours[run], ms=4,
                   label=f"{run} model" + (" (blind)" if run == "BHA12" else ""))
        lo, hi = DEPTHS[run]
        ax[0].hlines(meas[run], lo, hi, color=colours[run], lw=3, alpha=0.4)
    ax[0].fill_between([8700, 9255], 150, 160, color="tab:green", alpha=0.15, lw=0)
    ax[0].set(xlabel="Bit depth [ft MD]", ylabel="MWD temperature [F]",
              title="Modelled (points) and measured run averages (bars)")
    ax[0].legend(fontsize=8)
    t0, rows = o["decomposition"]
    eav = [-13, 23, -20, 7, -48]
    y = np.arange(len(rows))
    ax[1].barh(y + 0.2, [r[1] for r in rows], 0.4, label="Model 1")
    ax[1].barh(y - 0.2, eav, 0.4, label="Eavor (empirical)", color="0.6")
    ax[1].set_yticks(y, [r[0] for r in rows], fontsize=8)
    ax[1].axvline(0, color="k", lw=0.8)
    ax[1].set(xlabel="Change in MWD temperature [F]", title="BHA 10 to BHA 11, one factor at a time")
    ax[1].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


GOLDEN = os.path.join(HERE, "..", "tests", "golden", "validate_forge16b.json")


def golden(o):
    s1 = o["step1"]
    t0, rows = o["decomposition"]
    return dict(film_mult=s1["film_mult"], k_idp=o["k_idp"],
                run_means={r: mean_mwd(o["series"][r]) for r in o["series"]},
                bha12_late=o["bha12_late"],
                decomposition={name: d for name, d, _ in rows})


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    o = run_all()
    report(o)
    figure(o, os.path.join(HERE, "..", "figures", "validate_forge16b.png"))
    print("\nsaved validate_forge16b.png")
    if "--write-golden" in sys.argv:
        import json
        with open(GOLDEN, "w") as fh:
            json.dump(golden(o), fh, indent=1)
        print("wrote", os.path.basename(GOLDEN))
