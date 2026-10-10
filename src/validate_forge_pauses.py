"""
validate_forge_pauses.py
========================
A check of the wall solver (wall_thermal.py) against the return temperature
FORGE 16B recorded when circulation resumed after a pause.

When circulation stops, the fluid in the hole drifts towards the temperature of
the rock around it: it warms in the deep, hot part of the hole and cools near
surface. When circulation resumes, the annulus column is pushed to surface,
shallowest fluid first, and the return temperature over the first bottoms-up
traces it. The bottom of the hole alone can't be compared with the returns: on
the way back in from a trip, the bit circulates from stops above bottom, and
the fluid exactly at the bottom is a few hundred gallons in a bottoms-up of
25,000. So the whole column is modelled.

Every gap in circulation at the bottom of the hole of more than 20 minutes in
the trial record (21-29 May 2023) is used: pauses with the bit on bottom, and
trips.

The model:
  * At points every 200 ft from surface to the bottom, the wall solver is run
    through that depth's history, minute by minute from the Pason record:
    circulating when the pumps are on and the bit is below it (at the annulus
    temperature and film coefficient Model 1 gives there for that flow, inlet
    and bit depth, with the run's string and the annulus film unfitted, as the
    sites use it), static otherwise, with the fluid at that depth a well-mixed
    node coupled to the wall by the laminar film (Nu = 4.36), around the string
    or filling the hole when the string is out. Depths drilled before the
    record begins have their earlier circulating time, from the 10-minute
    drilling history, as one circulating spell.
  * At the end of a pause the column of fluid temperatures is pushed to surface
    as plug flow by the measured pump output.
  * The anomaly is the excess over the return temperature steady circulation
    would give at the resumed flow and inlet. For the model, the excess of the
    column over Model 1's annulus profile at the resumed state. For the record,
    the measured return less Model 1's return at the resumed state, shifted by
    Model 1's error on the return just before the pause.
  * The returning fluid exchanges heat with the string and the rock on its way
    up, which plug flow ignores, so one attenuation factor f scales the model's
    anomaly. The pauses shorter than 6 h can't be used to fit it: most were
    gyro pump-downs or circulation at 110 to 340 gal/min, so their first
    bottoms-up runs for hours into the next trip. f is fitted on the trips
    instead, each trip predicted with f from the others (leave one out). The
    short pauses are reported, excluded from the fit.
  * D2's no-node sensitivity (the fluid at the wall's temperature) is run
    the same way.

Simplifications, each stated in the README: plug flow; casing and cement are
not resolved (the cased section is treated as rock behind the casing's inner
radius); the fluid inside the string and its displacement by the string on
the way in are not modelled, so only the first bottoms-up (the annulus column)
is compared.
"""
import csv
import os
import sys
from datetime import datetime, timezone

import numpy as np

import model1_coupled as m1
import validate_forge16b as v
import wall_thermal as wt

HERE = os.path.dirname(os.path.abspath(__file__))
PUMPS_ON_GPM = 100.0
ON_BOTTOM_FT = 100.0
K_IDP_UNFITTED = 4.703      # validate_forge16b, the chain with the film unfitted
DZ_FT = 200.0
SHOE_FT = 4837.0
R_OPEN = 9.5 * v.INCH / 2
R_CASED = 10.682 * v.INCH / 2


def _minutes():
    out = []
    for r in v.PASON:
        try:
            out.append(dict(t=datetime.strptime(r["time"], "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc).timestamp(),
                            hole=float(r["hole_depth_ft"]), bit=float(r["bit_depth_ft"]),
                            flow=float(r["flow_gpm"]), tin=float(r["temp_in_F"]),
                            tout=float(r["temp_out_F"]),
                            dp=float(r["diff_pressure_psi"] or 0.0)))
        except ValueError:
            continue
    return out


M = _minutes()
T_INDEX = {m["t"]: i for i, m in enumerate(M)}


# The Pason clock carries no zone. The times are read as UTC so that the epochs
# do not depend on the machine's time zone and print as the record reads.
def _ts(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc).timestamp()


CYCLE = list(csv.DictReader(open(os.path.join(v.DATA, "cycle_timings.csv"))))
_runs = [r for r in CYCLE if r["kind"] == "bit run"]
RUN_ENDS = {"BHA10": _ts(_runs[1]["start"]), "BHA11": _ts(_runs[2]["start"]),
            "BHA12": _ts(_runs[3]["start"])}


def run_for(ts):
    """The string in the hole at time ts (the gyro and coring runs after BHA 13
    used regular pipe with no collars, as BHA 13 did)."""
    for run in ("BHA10", "BHA11", "BHA12"):
        if ts < RUN_ENDS[run]:
            return run
    return "BHA13"


def pauses():
    """Every gap in circulation at the bottom of the hole over 20 minutes."""
    return [dict(kind=r["kind"], start=_ts(r["start"]), end=_ts(r["end"]),
                 hours=float(r["hours"]), depth_ft=float(r["depth_ft"]))
            for r in CYCLE if r["kind"] in ("pause on bottom", "trip")
            and float(r["hours"]) * 60 > 20]


# ------------------------------------------------------------ circulating states
_cache = {}


def circulating_profile(ts, bit_ft, flow, tin, dp):
    """Model 1 for a circulating spell: measured depth [ft] along the hole,
    annulus temperature [C] and annulus film coefficient, and the return
    temperature [C]. Cached on the run, bit depth to 200 ft, flow to 25 gal/min
    and inlet to 2 F."""
    run = run_for(ts)
    bit = max(round(bit_ft / 200.0) * 200.0, 1200.0)
    key = (run, bit, round(flow / 25) * 25, round(tin / 2) * 2)
    if key not in _cache:
        s = dict(time=ts, bit_ft=bit, flow_gpm=key[2], T_in_F=key[3], motor_dp_psi=max(dp, 0.0))
        r, g = v.solve_at(run, s, 1.0, K_IDP_UNFITTED)
        z = r["z"]
        ga = g.at(np.clip(z, 0.5, g.md_bit - 0.5))
        h = m1.h_dittus(r["m_dot"], ga["A_ann"], ga["Dh_ann"], r["Tu"], m1.P_of_z(ga["tvd"]))
        _cache[key] = (z / v.FT, r["Tu"], h, float(r["Tu"][0]), r["success"])
    return _cache[key]


def _mean(rows, k):
    return float(np.mean([m[k] for m in rows]))


def static_node(z_ft, string_in):
    """Heat capacity [J/m/K] of the water at depth z_ft and its laminar film
    coefficient to the wall, around the string or filling the hole."""
    tvd = v.SURVEY.tvd(z_ft * v.FT)
    P, T = m1.P_of_z(np.array([tvd])), np.array([float(v.T_FORM(tvd))])
    rho, cp, mu, k = (float(x[0]) for x in m1.fluid_props(T, P))
    a = R_CASED if z_ft < SHOE_FT else R_OPEN
    if string_in:
        r_out = 5.5 * v.INCH / 2
        area, Dh = np.pi * (a ** 2 - r_out ** 2), 2 * (a - r_out)
    else:
        area, Dh = np.pi * a ** 2, 2 * a
    return rho * cp * area, wt.NU_LAMINAR * k / Dh


NODE = True     # False: the no-node sensitivity of D2 (fluid at the wall's temperature)


def element(z_ft, t_end):
    """The wall solver at depth z_ft through its history to t_end. Returns the
    times and fluid temperatures [C] at the end of each state."""
    first = next((i for i, m in enumerate(M) if m["hole"] >= z_ft), None)
    if first is None:
        return None
    i_end = T_INDEX[t_end]
    flags = [M[i]["flow"] >= PUMPS_ON_GPM and M[i]["bit"] >= z_ft for i in range(len(M))]
    states, kinds, i = [], [], first
    # depths drilled before the record: their earlier circulation as one spell
    t_clock, circ, depth = v.CLOCK
    if first == 0 or M[first]["hole"] - z_ft > 50.0:
        t_reach = float(np.interp(z_ft * v.FT, depth, t_clock))
        pre = float(np.interp(M[0]["t"], t_clock, circ) - np.interp(t_reach, t_clock, circ))
        first_circ = next(k for k in range(len(M)) if flags[k])
        prof = circulating_profile(M[first_circ]["t"], M[first_circ]["bit"], M[first_circ]["flow"],
                                   M[first_circ]["tin"], M[first_circ]["dp"])
        if pre > 0:
            states.append(wt.State("circulating", pre, T_fluid=float(np.interp(z_ft, prof[0], prof[1])),
                                   h=float(np.interp(z_ft, prof[0], prof[2]))))
            kinds.append(("circulating", M[0]["t"]))
        i = 0
    while i < i_end:
        c = flags[i]
        j = i
        while j + 1 < i_end and flags[j + 1] == c:
            j += 1
        rows = M[i:j + 1]
        dur = (j - i + 1) * 60.0
        if c:
            prof = circulating_profile(rows[0]["t"], _mean(rows, "bit"), _mean(rows, "flow"),
                                       _mean(rows, "tin"), _mean(rows, "dp"))
            states.append(wt.State("circulating", dur, T_fluid=float(np.interp(z_ft, prof[0], prof[1])),
                                   h=float(np.interp(z_ft, prof[0], prof[2]))))
        else:
            string_in = float(np.median([m["bit"] for m in rows])) > z_ft
            C, hs = static_node(z_ft, string_in)
            states.append(wt.State("static", dur, h=hs, C=C if NODE else 0.0))
        kinds.append(("circulating" if c else "static", M[j]["t"] + 60.0))
        i = j + 1
    tvd = v.SURVEY.tvd(z_ft * v.FT)
    a = R_CASED if z_ft < SHOE_FT else R_OPEN
    res = wt.run(states, float(v.T_FORM(tvd)), a, v.K_ROCK, v.RHOCP_ROCK, n=100, growth=1.25)
    out = []
    for (kind, t_state_end), s, (_, T_w, T_n) in zip(kinds, states, res.ends):
        out.append((t_state_end, s.T_fluid if kind == "circulating" else (T_n if NODE else T_w)))
    return out


def column(depths, t_end):
    """Fluid temperature [C] at each depth when circulation resumes at t_end."""
    col = []
    for z in depths:
        e = ELEMENTS[z]
        val = None
        for t, T in e:
            if t <= t_end + 61.0:
                val = T
        col.append(val)
    return np.array(col, float)


# --------------------------------------------------------------- the comparison
def anomalies(p, f=1.0):
    """Observed and modelled return-temperature anomalies [F] over the first
    bottoms-up after circulation resumes, minute by minute."""
    depth = p["depth_ft"]
    depths = np.arange(DZ_FT / 2, depth, DZ_FT)
    col = column(depths, p["end"])
    i0, i1 = T_INDEX[p["start"]], T_INDEX[p["end"]]
    before = [m for m in M[max(i0 - 10, 0):i0 + 1] if m["flow"] >= PUMPS_ON_GPM]
    after = [m for m in M[i1:i1 + 30] if m["flow"] >= PUMPS_ON_GPM]
    pb = circulating_profile(before[0]["t"], _mean(before, "bit"), _mean(before, "flow"),
                             _mean(before, "tin"), _mean(before, "dp"))
    pa = circulating_profile(after[0]["t"], depth, _mean(after, "flow"),
                             _mean(after, "tin"), _mean(after, "dp"))
    offset = _mean(before, "tout") - (1.8 * pb[3] + 32)        # Model 1's error, F
    ref_F = 1.8 * pa[3] + 32 + offset
    anom_col = 1.8 * (col - np.interp(depths, pa[0], pa[1]))   # F, by depth
    # cumulative annulus volume from surface to each depth, at the bit on bottom
    g = v.geometry(run_for(p["end"]), depth, K_IDP_UNFITTED)
    zz = np.linspace(0, depth, 400)
    a = g.at(np.clip(zz * v.FT, 0.5, g.md_bit - 0.5))["A_ann"]
    vol = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * np.diff(zz * v.FT))]) / 3.785e-3
    v_total = float(vol[-1])
    obs, mod, pumped, k = [], [], 0.0, i1
    while pumped < v_total and k < len(M) - 1:
        pumped += max(M[k]["flow"], 0.0)
        z_here = float(np.interp(pumped, vol, zz))
        obs.append(M[k]["tout"] - ref_F)
        mod.append(f * float(np.interp(z_here, depths, anom_col)))
        k += 1
    return dict(obs=np.array(obs), mod=np.array(mod), minutes=len(obs), v_total=v_total,
                ref_F=ref_F, offset=offset, anom_col=anom_col, depths=depths)


ELEMENTS = {}


def build(t_end):
    """Run the wall solver at every depth through the record to t_end."""
    zmax = max(m["hole"] for m in M)
    for z in np.arange(DZ_FT / 2, zmax, DZ_FT):
        if z not in ELEMENTS:
            e = element(float(z), t_end)
            if e is not None:
                ELEMENTS[float(z)] = e


def _f(rows, key):
    return sum(r["obs"] * r[key] for r in rows) / sum(r[key] ** 2 for r in rows)


def check(node=True):
    global NODE
    NODE = node
    ELEMENTS.clear()
    ps = pauses()
    build(max(p["end"] for p in ps))
    rows = []
    for p in ps:
        a = anomalies(p)
        rows.append(dict(p, obs=float(np.mean(a["obs"])), mod=float(np.mean(a["mod"])),
                         minutes=a["minutes"], offset=a["offset"], curve=a))
    trips = [r for r in rows if r["kind"] == "trip"]
    for r in rows:
        if r["kind"] == "trip":
            r["f_loo"] = _f([o for o in trips if o is not r], "mod")
            r["pred"] = r["f_loo"] * r["mod"]
        else:
            r["f_loo"], r["pred"] = None, None
    NODE = True
    return dict(rows=rows, f=_f(trips, "mod"), node=node)


def report(c, c0=None):
    print("=" * 112)
    print("FORGE 16B: return-temperature anomaly over the first bottoms-up after a pause, "
          "observed and modelled [F]")
    print("=" * 112)
    print(f"{'start':<14}{'kind':<17}{'hours':>6}{'depth':>7}{'BU min':>8}{'model err':>10}"
          f"{'observed':>10}{'plug flow':>10}{'x f':>8}{'error':>8}{'no node':>9}")
    for i, r in enumerate(c["rows"]):
        nn = "" if c0 is None else f"{c0['rows'][i]['mod']:9.1f}"
        if r["kind"] == "trip":
            print(f"{datetime.fromtimestamp(r['start'], timezone.utc).strftime('%d May %H:%M'):<14}{r['kind']:<17}"
                  f"{r['hours']:6.1f}{r['depth_ft']:7.0f}{r['minutes']:8d}{r['offset']:+10.1f}"
                  f"{r['obs']:10.1f}{r['mod']:10.1f}{r['pred']:8.1f}{r['pred'] - r['obs']:+8.1f}{nn}")
        else:
            print(f"{datetime.fromtimestamp(r['start'], timezone.utc).strftime('%d May %H:%M'):<14}{r['kind']:<17}"
                  f"{r['hours']:6.1f}{r['depth_ft']:7.0f}{r['minutes']:8d}{r['offset']:+10.1f}"
                  f"{r['obs']:10.1f}{r['mod']:10.1f}{'':8}{'excluded':>8}{nn}")
    print("-" * 112)
    trips = [r for r in c["rows"] if r["kind"] == "trip"]
    rms = np.sqrt(np.mean([(r["pred"] - r["obs"]) ** 2 for r in trips]))
    print(f"attenuation f on the trips: {c['f']:.3f} (each trip left out of its own fit: "
          + ", ".join(f"{r['f_loo']:.2f}" for r in trips) + f"); RMS error {rms:.1f} F")
    if c0 is not None:
        t0 = [r for r in c0["rows"] if r["kind"] == "trip"]
        rms0 = np.sqrt(np.mean([(r["pred"] - r["obs"]) ** 2 for r in t0]))
        print(f"no fluid node: f {c0['f']:.3f}, RMS error {rms0:.1f} F; predicted minus observed "
              + ", ".join(f"{r['pred'] - r['obs']:+.1f}" for r in t0))
    print("excluded: pauses with the bit on bottom, whose first bottoms-up runs for hours at low "
          "flow (gyro pump-downs and the like)")
    print("model err: Model 1's return temperature less the measured one before the pause, "
          "removed from the record")


def figure(c, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    trips = [r for r in c["rows"] if r["kind"] == "trip"]
    fig, ax = plt.subplots(1, len(trips) + 1, figsize=(4 * (len(trips) + 1), 4.2))
    for r in trips:
        ax[0].plot(r["obs"], r["pred"], "s", color="tab:red")
        ax[0].annotate(f"{r['hours']:.0f} h", (r["obs"], r["pred"]), fontsize=7,
                       xytext=(4, 2), textcoords="offset points")
    vals = [x for r in trips for x in (r["obs"], r["pred"])]
    lim = (min(vals) - 3, max(vals) + 3)
    ax[0].plot(lim, lim, "k:", lw=1)
    ax[0].set(xlim=lim, ylim=lim, xlabel="observed mean anomaly [F]", ylabel="modelled [F]",
              title="trips, each predicted with f\nfitted on the others")
    for a, r in zip(ax[1:], trips):
        cu = r["curve"]
        t = np.arange(cu["minutes"])
        a.plot(t, cu["obs"], "k-", lw=1.5, label="observed")
        a.plot(t, r["f_loo"] * cu["mod"], "tab:red", lw=1.5, label="model x f")
        a.plot(t, cu["mod"], "tab:red", lw=0.8, ls=":", label="plug flow")
        a.set(xlabel="minutes after circulation resumed",
              title=f"trip of {r['hours']:.0f} h, {datetime.fromtimestamp(r['start'], timezone.utc).strftime('%d May')}")
        a.legend(fontsize=7)
    ax[1].set_ylabel("return temperature anomaly [F]")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


GOLDEN = os.path.join(HERE, "..", "tests", "golden", "forge_pauses.json")

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    c = check(True)
    c0 = check(False)
    report(c, c0)
    figure(c, os.path.join(HERE, "..", "figures", "forge_pauses.png"))
    print("saved forge_pauses.png")
    if "--write-golden" in sys.argv:
        import json
        keep = lambda cc: [{k: r[k] for k in ("start", "kind", "hours", "obs", "mod", "pred")}
                           for r in cc["rows"]]
        with open(GOLDEN, "w") as fh:
            json.dump(dict(f=c["f"], f_no_node=c0["f"], rows=keep(c), rows_no_node=keep(c0)),
                      fh, indent=1)
        print("wrote", os.path.basename(GOLDEN))
