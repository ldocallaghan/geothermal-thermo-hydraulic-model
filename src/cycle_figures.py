"""Figures for the drilling cycle at Soultz: the wall through a bit run and the
trip that ends it (Fig. A), and the failure zone around the hole at the end of
drilling and the end of the trip (Fig. B)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
import wall_thermal

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")


def main():
    site = se.SOULTZ
    o = se.evaluate(site)
    c = o["cycle"]
    z, T_rock = o["z"], o["T_rock"]
    gov = next(pc for pc in c["per_case"] if pc["label"] == c["deciding_case"])
    p, ucs = gov["profile"], gov["UCS"]
    sg = C.MUD_SG_GRAD * z
    thermo = -se.thermal_hoop_stress(T_rock, T_rock - 1.0)

    # the history again, with the temperature at the breakout's depth probed
    states, kinds = [], []
    for st in c["history"]:
        if st.kind == "circulating":
            states.append(wall_thermal.State("circulating", st.hours * 3600, T_fluid=c["T_fluid"], h=c["h"]))
        else:
            r_str = 0.07 if st.label == "connection" else 0.0
            cap, hs = se._static_node(T_rock, z, c["a"], r_str)
            states.append(wall_thermal.State("static", st.hours * 3600, h=hs, C=cap))
        kinds.append(st)
    res = wall_thermal.run(states, T_rock, c["a"], site.k_rock, C.RHO_ROCK * C.CP_ROCK, n=150,
                           growth=1.2, probe_r=c["r_ref"])
    t_h = res.t / 3600.0
    # which state each step is in (for circulating friction on the bounds)
    bounds = np.cumsum([0.0] + [st.hours for st in kinds])
    circ = np.array([kinds[min(np.searchsorted(bounds, t, side="right") - 1, len(kinds) - 1)].kind
                     == "circulating" for t in t_h])
    Pw_d = c["drill_SG"] * sg
    width = np.array([m5.breakout_width(p.SHmax(z), p.Shmin(z), Pw_d + (c["ecd_SG"] * sg if ci else 0.0),
                                        float(p.Pp(z)), T, ucs, thermo * (T - T_rock))
                      for T, ci in zip(res.T_probe, circ)])
    lo = np.array([max(m5.mud_for_width(p.SHmax(z), p.Shmin(z), float(p.Pp(z)), T, C.BREAKOUT_W_MAX_DEG,
                                        ucs, thermo * (T - T_rock)), float(p.Pp(z))) / sg
                   for T in res.T_probe]) - np.where(circ, c["ecd_SG"], 0.0)
    hi = (p.Shmin(z) - C.MUD_MARGIN_SG * sg) / sg - np.where(circ, c["ecd_SG"], 0.0)

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(3, 2, width_ratios=(3, 1))
    ax = [fig.add_subplot(gs[0, 0])]
    ax += [fig.add_subplot(gs[i, 0], sharex=ax[0]) for i in (1, 2)]
    axp = fig.add_subplot(gs[:, 1])
    ax[0].plot(t_h, res.T_wall, "tab:blue", label="wall")
    ax[0].plot(t_h, res.T_probe, "tab:red", label=f"rock {100 * (c['r_ref'] - c['a']):.1f} cm behind the wall")
    ax[0].axhline(T_rock, color="k", ls=":", lw=1, label="formation")
    ax[0].set_ylabel("temperature [C]")
    ax[0].legend(fontsize=8)
    ax[1].plot(t_h, width, "k")
    ax[1].axhline(C.BREAKOUT_W_MAX_DEG, color="tab:red", ls="--", lw=1, label="90 deg limit")
    ax[1].set_ylabel(f"breakout width at {c['drill_SG']:.2f} SG [deg]")
    ax[1].legend(fontsize=8)
    ax[2].fill_between(t_h, lo, hi, where=hi >= lo, color="tab:green", alpha=0.3,
                       label="workable static fluid weight")
    ax[2].plot(t_h, lo, "tab:green")
    ax[2].plot(t_h, hi, "tab:red", label="fracture limit")
    ax[2].axhline(c["drill_SG"], color="k", lw=1, label="drilling fluid")
    ax[2].set_ylabel("fluid weight [SG]")
    ax[2].set_xlabel("hours since the bit passed this depth")
    ax[2].legend(fontsize=8)
    t_trip = bounds[-2]
    holds = c["safe_pause_h"] >= bounds[-1] - t_trip
    for a in ax:
        a.axvspan(t_trip, bounds[-1], color="0.9", zorder=0)
        if not holds:
            a.axvline(t_trip + c["safe_pause_h"], color="tab:orange", lw=1.2)
        a.grid(alpha=0.3)
    note = (f"the drilling fluid holds the wall through the trip (safe pause {c['safe_pause_h']:.0f} h)"
            if holds else f"safe pause {c['safe_pause_h']:.1f} h")
    ax[0].text(t_trip + 1, ax[0].get_ylim()[0] + 5, note, color="tab:orange", fontsize=8)
    ax[0].text(t_trip + 1, T_rock - 15, f"trip, {c['trip_h']:.0f} h without circulation", fontsize=8)
    # the trip fluid needed up the open hole
    tp = o["trip_profile"]
    hgt = [r["height"] / 1000 for r in tp]
    axp.plot([r["SG_needed"] for r in tp], hgt, "o-", color="tab:green", label="needed at the end of the trip")
    axp.plot([r["SG_hi"] for r in tp], hgt, "--", color="tab:red", label="fracture limit, pumps off")
    axp.set_xlabel("fluid weight [SG]")
    axp.set_ylabel("height above the bit [km]")
    axp.set_title("The trip, up the open hole", fontsize=9)
    axp.legend(fontsize=8, loc="upper right")
    axp.grid(alpha=0.3)
    for a in ax[:2]:
        plt.setp(a.get_xticklabels(), visible=False)
    fig.suptitle(f"Soultz at {z / 1000:.1f} km: the wall one stand above the end of a bit run "
                 f"({c['deciding_case']})", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "cycle_soultz.png"), dpi=130)

    # Fig. B: the failure margin around the hole at the end of drilling and of
    # the trip. Elastic Mohr-Coulomb fails a wider region just behind the wall
    # than at it, so the width is judged at the wall with the temperature at
    # the uncooled failure depth (the dashed circle).
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw=dict(projection="polar"),
                            gridspec_kw=dict(wspace=0.35))
    lev = np.linspace(-150, 150, 13)
    for axp, (name, field_i, Pw) in zip(axs, (("end of drilling", 3, Pw_d + c["ecd_SG"] * sg),
                                             ("end of the trip", 4, Pw_d))):
        b = m5.breakout_behind_wall(p.SHmax(z), p.Shmin(z), Pw, float(p.Pp(z)), c["a"], res.r,
                                    res.fields[field_i], T_rock, ucs, thermo, r_out=2.0, d_theta=1.0)
        T_r = float(np.interp(c["r_ref"], res.r, res.fields[field_i]))
        w = m5.breakout_width(p.SHmax(z), p.Shmin(z), Pw, float(p.Pp(z)), T_r, ucs,
                              thermo * (T_r - T_rock))
        ph = np.radians(b["phi"])
        # the quarter computed (0 to 90 deg from the Shmin azimuth), mirrored
        th = np.concatenate([ph, np.pi - ph[::-1], np.pi + ph, 2 * np.pi - ph[::-1]])
        M = b["margin"] / 1e6
        M = np.concatenate([M, M[:, ::-1], M, M[:, ::-1]], axis=1)
        R, TH = np.meshgrid(b["r"] / c["a"], th, indexing="ij")
        cf = axp.contourf(TH, R, np.clip(M, lev[0], lev[-1]), levels=lev, cmap="RdBu_r")
        axp.contour(TH, R, M, levels=[0.0], colors="k", linewidths=0.8)
        tt = np.linspace(0, 2 * np.pi, 200)
        axp.plot(tt, np.full_like(tt, c["r_ref"] / c["a"]), "k--", lw=1)
        for base in (0.0, np.pi):
            arc = np.radians(np.linspace(-w / 2, w / 2, 50)) + base
            axp.plot(arc, np.full_like(arc, 1.0), color="tab:orange", lw=4)
        axp.set_ylim(0, 2)
        axp.set_yticks([1.0, 1.5, 2.0])
        axp.set_xticks(np.radians([0, 90, 180, 270]))
        axp.set_xticklabels(["Shmin", "SHmax", "Shmin", "SHmax"])
        axp.tick_params(axis="x", pad=10)
        axp.set_rlabel_position(112)
        axp.set_title(f"{name}\n{T_r:.0f} C at the dashed circle, width {w:.0f} deg", fontsize=9, pad=18)
    cb = fig.colorbar(cf, ax=axs, shrink=0.7, pad=0.06)
    cb.set_label("Mohr-Coulomb margin [MPa], failure where positive")
    fig.suptitle(f"Soultz at {z / 1000:.1f} km, {c['drill_SG']:.2f} SG ({c['deciding_case']}): "
                 "\nthe rock out to twice the hole radius; the breakout width (orange) is judged at the wall",
                 fontsize=9)
    fig.savefig(os.path.join(FIG, "failure_zone_soultz.png"), dpi=130, bbox_inches="tight")
    print("saved cycle_soultz.png, failure_zone_soultz.png")


if __name__ == "__main__":
    main()
