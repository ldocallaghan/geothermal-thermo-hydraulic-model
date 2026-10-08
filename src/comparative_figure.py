"""Comparative siting figure (v1.1): what decides each verdict, and how much of
it rests on data.

Left: breakout width at hydrostatic mud across each site's stress and strength
cases, against depth to 400 C, with the 90 deg verdict limit (Zoback 2007)
and UD-1's site-calibrated 63 deg. Right: how deep each site's stress and
temperature data reach against its target. Filled = evidence-based, hollow =
speculative."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geo_constants as C
from comparative_sites import SITES, classify
from site_evaluation import evaluate

BLUE = "#2a78d6"                                    # reference palette slot 1
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

results = [evaluate(s) for s in SITES]
short = {s.name: s.name.split("(")[0].split("/")[0].strip() for s in SITES}
# the x-axis of the right panel has room for one short word per site
tick = {s.name: {"Upper Rhine Graben": "Soultz", "United Downs": "United\nDowns",
                 "Pannonian Basin": "Pannonian"}.get(short[s.name], short[s.name])
        for s in SITES}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.5), sharey=True,
                             gridspec_kw=dict(width_ratios=[1.35, 1]), facecolor=SURF)
for a in (a1, a2):
    a.set_facecolor(SURF)
    a.grid(True, color=GRID, lw=0.8)
    a.tick_params(colors=INK2)
    for sp in a.spines.values():
        sp.set_color(GRID)

# left: breakout width range at hydrostatic mud
for s, o in zip(SITES, results):
    ws = [c["window"]["width_hydro"] for c in o["stress_cases"]]
    y = o["z"] / 1000
    filled = s.tier == "evidence-based"
    a1.plot([min(ws), max(ws)], [y, y], color=BLUE, lw=2, solid_capstyle="round", zorder=2)
    a1.scatter([o["stress"]["window"]["width_hydro"]], [y], s=90, zorder=3,
               color=BLUE if filled else SURF, edgecolor=BLUE, linewidth=2)
    verdict = classify(o)[2]
    a1.annotate(f"{short[s.name]}  ({s.tier})\n{verdict}", (max(ws), y),
                xytext=(8, 0), textcoords="offset points", va="center",
                fontsize=8, color=INK)
for W, lab, yy in ((C.BREAKOUT_W_MAX_DEG, "W_max 90", 0.8),
                   (C.BREAKOUT_W_MAX_SITE_DEG, "UD-1 calibrated 63", 1.6)):
    a1.axvline(W, color=INK2, lw=0.8, ls="--")
    a1.text(W - 2, yy, lab, color=INK2, fontsize=8, ha="right", va="top")
a1.set_xlim(-5, 235)
a1.set_xlabel("breakout width at hydrostatic mud (deg), range across cases", color=INK)
a1.set_ylabel("depth to 400 C (km)", color=INK)
a1.set_title("Breakout width decides the verdict (marker = the row's case)",
             color=INK, fontsize=10, loc="left")
a1.text(0, 14.1, "filled: evidence-based   hollow: speculative", color=INK2, fontsize=8)

# right: how deep the data reach
for i, (s, o) in enumerate(zip(SITES, results)):
    x = i + 1
    filled = s.tier == "evidence-based"
    zd = o["stress"]["profile"].z_data
    if zd is not None:
        a2.plot([x - 0.12] * 2, [zd[0] / 1000, zd[1] / 1000], color=BLUE, lw=6,
                solid_capstyle="butt")
    if s.temperature_data_to is not None:
        a2.plot([x + 0.12] * 2, [0, s.temperature_data_to / 1000], color=BLUE, lw=6,
                alpha=0.35, solid_capstyle="butt")
    deepest = max([zd[1] if zd else 0.0, s.temperature_data_to or 0.0]) / 1000
    a2.plot([x, x], [deepest, o["z"] / 1000], color=INK2, lw=1, ls=":")
    a2.scatter([x], [o["z"] / 1000], s=90, zorder=3,
               color=BLUE if filled else SURF, edgecolor=BLUE, linewidth=2)
    a2.text(x, o["z"] / 1000 + 0.5, "target", color=INK2, fontsize=7, ha="center", va="top")
    if zd is None:
        a2.text(x - 0.12, 0.6, "no stress\ndata", color=INK2, fontsize=7, ha="center")
a2.set_xticks(range(1, len(SITES) + 1))
a2.set_xticklabels([tick[s.name] for s in SITES], fontsize=8, color=INK)
a2.set_xlim(0.4, len(SITES) + 0.6)
a2.text(1.55, 7.0, "dark: stress data\nlight: temperature data\n"
        "dotted: extrapolation\nto the target", color=INK2, fontsize=8)
a2.set_title("How deep the data reach", color=INK, fontsize=10, loc="left")

a1.set_ylim(14.5, 0)
fig.tight_layout()
_figdir = os.path.join(os.path.dirname(__file__), "..", "figures")
fig.savefig(os.path.join(_figdir, "comparative_sites.png"), dpi=130, facecolor=SURF)
print("saved comparative_sites.png")
