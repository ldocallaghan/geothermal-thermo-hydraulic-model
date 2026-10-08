"""
benchmark_utaustin.py
=====================
Code-to-code benchmark of Model 1 against Wu et al. (2025), "A Comprehensive
Evaluation of Drill Pipe Insulation for Downhole Temperature Management Using
Physics-Based Models" (Stanford Geothermal Workshop, SGP-TR-202). Their inputs
and results are in data/benchmarks/utaustin2025.md.

Their wells are rebuilt from the published inputs and Model 1 is run with
nothing fitted:

  * the vertical and horizontal base cases (their Fig. 2), each pipe type;
  * the vertical well against reservoir temperature at 600 gal/min (Fig. 5);
  * the vertical 400 C well against flow rate (Fig. 6);
  * standpipe pressure for the vertical well at 600 gal/min (Figs. 7 and 9).

Model 1 is run with their constant mud properties, for a like-for-like
comparison, and with IAPWS-95 water at the same volumetric flow. Heat from
friction in the string, annulus and bit nozzles is included: their dual-wall
pipe loses about 70 MPa to friction, which warms the mud by about 18 C.

Two readings of their coated pipes are run. Their text describes a 1 mm coating
(the model's named pipe types). Their Table 5 lists a single conductivity for
each pipe, which their model may apply across the whole wall; the "whole wall"
column runs it that way, as a diagnostic only.

Assumptions where the paper is silent, each stated here:

  * The mud is treated as Newtonian, with its plastic viscosity (14 mPa s).
  * Their transient run lasts 1,500 minutes with the hole already drilled, so
    the whole wall has been exposed for 1,500 minutes.
  * No heat is added at the bit.
  * Casings run to surface, the gap between two casings is cement, and rock
    starts at the outermost casing's OD (the drilled hole sizes aren't given).
  * The bit's nozzles aren't given; the standpipe pressure uses the model's
    default (the FORGE 16B trial bits' 1.203 in^2).
  * The horizontal well's 400 m curve builds from 0 to 90 deg at constant
    curvature and lands at 8,600 m TVD; the drill pipe takes up the extra length.
"""
import os
import sys

import numpy as np

import model1_coupled as m1
import well_geometry as wg

MUD = dict(rho=1078.0, cp=3750.0, mu=0.014, k=0.75)
RHO_WATER_INLET = 992.2   # kg/m^3, water at 40 C, 1 atm (IAPWS-95)
GPM = 6.30902e-5          # m^3/s
T_SURF, K_ROCK, RHOCP_ROCK = 17.0, 2.31, 2800.0 * 930.0
TVD = 8600.0
T_INLET = 40.0
EXPOSURE = 1500 * 60.0    # s

HWDP, COLLAR = wg.WU_HWDP, wg.WU_COLLAR

PIPES = {"conventional": wg.CONVENTIONAL, "internal": wg.INTERNALLY_COATED,
         "external": wg.EXTERNALLY_COATED, "dual-wall": wg.DUAL_WALL}
WHOLE_WALL = {
    "conventional": wg.CONVENTIONAL,
    "internal": wg.PipeType("internal, whole wall", (wg.Layer(0.1194 / 2, 0.1397 / 2, 0.47),)),
    "external": wg.PipeType("external, whole wall", (wg.Layer(0.1214 / 2, 0.1417 / 2, 1.31),)),
    "dual-wall": wg.PipeType("dual-wall, whole wall", (wg.Layer(0.040, 0.1397 / 2, 0.092),)),
}

# Their results (data/benchmarks/utaustin2025.md); None where not given
FIG2 = {"vertical": {"conventional": 179, "internal": 76, "external": 121, "dual-wall": 59},
        "horizontal": {"conventional": 206, "internal": 100, "external": 153, "dual-wall": 65}}
FIG5_T = (200, 250, 300, 350, 400)
FIG5_G = (21.3, 27.1, 32.9, 38.7, 44.5)   # C/km, their Table 6
FIG5 = {"internal": (71, 80, 89, 98, 107), "external": (110, 114, 132, 148, 165),
        "dual-wall": (58, 59, 60, 61, 62)}
FIG6_GPM = (400, 500, 600, 700, 800)
FIG6 = {"internal": (145, 123, 107, 95, 87), "external": (212, 186, 165, 148, 135),
        "dual-wall": (59, 60, 61, 63, 66)}
SPP_600 = {"internal": (16, 18), "external": (16, 18), "dual-wall": (67, 67)}  # MPa, read off Fig. 7


def well(pipe, horizontal=False):
    if horizontal:
        curve, lateral = 400.0, 2000.0
        kop = TVD - 2 * curve / np.pi
        md = kop + curve + lateral
        survey = wg.Survey.from_stations(
            np.concatenate([[0.0, kop], kop + np.linspace(0, curve, 41)[1:], [md]]),
            np.concatenate([[0.0, 0.0], np.linspace(0, 90, 41)[1:], [90.0]]))
    else:
        md, survey = TVD, None
    return wg.wu2025_well(md, pipe, survey)


def run(pipe, G_per_km=24.0, gpm=600.0, fluid=MUD, horizontal=False, friction_heat=True):
    Q = gpm * GPM
    m_dot = Q * (fluid["rho"] if fluid is not None else RHO_WATER_INLET)
    g = well(pipe, horizontal)
    G = G_per_km / 1000.0
    r = m1.solve(m_dot=m_dot, T_inj=T_INLET, G=G, T_surf=T_SURF,
                 target_rock_T=T_SURF + G * TVD, Q_face=0.0, k_rock=K_ROCK,
                 geometry=g, fluid=fluid, rhocp_rock=RHOCP_ROCK, verbose=False,
                 friction_heat=friction_heat,
                 exposure=lambda s: np.full_like(np.asarray(s, float), EXPOSURE))
    r["hyd"] = m1.hydraulics(m_dot, g, result=r, fluid=fluid)
    return r


def bhct(r):
    return float(r["Td"][-1])


def tables():
    out = {}
    print("=" * 78)
    print("Model 1 against Wu et al. (2025): BHCT [C], nothing fitted")
    print("=" * 78)
    print("\nBase case, 600 gal/min, 24 C/km, 1,500 min (their Fig. 2)")
    print("  mud: their properties; water: IAPWS-95; no fric.: mud without friction")
    print("  heating; whole wall: Table 5 conductivity across the whole pipe wall, mud")
    print(f"{'well':>11}{'pipe':>14}{'Wu':>6}{'mud':>7}{'diff':>6}{'water':>7}{'diff':>6}"
          f"{'no fric.':>9}{'whole wall':>11}{'diff':>6}")
    for orient in ("vertical", "horizontal"):
        h = orient == "horizontal"
        for name, pipe in PIPES.items():
            a = bhct(run(pipe, horizontal=h))
            b = bhct(run(pipe, fluid=None, horizontal=h))
            c = bhct(run(pipe, horizontal=h, friction_heat=False))
            d = bhct(run(WHOLE_WALL[name], horizontal=h))
            ref = FIG2[orient][name]
            out[("fig2", orient, name)] = (ref, a, b, c, d)
            print(f"{orient:>11}{name:>14}{ref:6.0f}{a:7.0f}{a-ref:+6.0f}{b:7.0f}{b-ref:+6.0f}"
                  f"{c:9.0f}{d:11.0f}{d-ref:+6.0f}")

    print("\nVertical, 600 gal/min, against reservoir temperature (their Fig. 5), mud")
    print("  Model 1 / Model 1 whole wall (Wu)")
    print(f"{'pipe':>14}" + "".join(f"{t:>15}" for t in FIG5_T))
    for name, refs in FIG5.items():
        row = []
        for G, ref in zip(FIG5_G, refs):
            a = bhct(run(PIPES[name], G_per_km=G))
            d = bhct(run(WHOLE_WALL[name], G_per_km=G))
            out[("fig5", name, G)] = (ref, a, d)
            row.append(f"{a:3.0f}/{d:3.0f} ({ref:3d})")
        print(f"{name:>14}" + "".join(f"{c:>15}" for c in row))

    print("\nVertical 400 C well against flow (their Fig. 6), mud")
    print("  Model 1 / Model 1 whole wall (Wu)")
    print(f"{'pipe':>14}" + "".join(f"{q:>15}" for q in FIG6_GPM))
    for name, refs in FIG6.items():
        row = []
        for q, ref in zip(FIG6_GPM, refs):
            a = bhct(run(PIPES[name], G_per_km=44.5, gpm=q))
            d = bhct(run(WHOLE_WALL[name], G_per_km=44.5, gpm=q))
            out[("fig6", name, q)] = (ref, a, d)
            row.append(f"{a:3.0f}/{d:3.0f} ({ref:3d})")
        print(f"{name:>14}" + "".join(f"{c:>15}" for c in row))

    print("\nStandpipe pressure, vertical, 600 gal/min, mud [MPa] (their Figs. 7, 9)")
    print("  bit: the model's default nozzles; theirs aren't given")
    print(f"{'pipe':>14}{'Wu':>9}{'M1':>8}{'string':>8}{'annulus':>8}{'bit':>6}{'buoy':>6}"
          f"{'v_ann m/s':>10}")
    for name, pipe in PIPES.items():
        h = run(pipe)["hyd"]
        ref = SPP_600.get(name)
        out[("spp", name)] = (ref, h["spp"] / 1e6)
        rs = "  n/a" if ref is None else f"{ref[0]}-{ref[1]}" if ref[0] != ref[1] else f"{ref[0]}"
        print(f"{name:>14}{rs:>9}{h['spp']/1e6:8.1f}{h['dp_string']/1e6:8.1f}"
              f"{h['dp_annulus']/1e6:8.1f}{h['dp_bit']/1e6:6.1f}{h['dp_buoyancy']/1e6:6.1f}"
              f"{h['v_ann_min']:10.2f}")
    return out


def figure(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colours = {"conventional": "k", "internal": "tab:blue", "external": "tab:orange",
               "dual-wall": "tab:green"}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    for (kind, orient, name), (ref, a, b, c, d) in ((k, v) for k, v in out.items()
                                                     if k[0] == "fig2"):
        mk = "o" if orient == "vertical" else "s"
        ax[0].plot(ref, a, mk, color=colours[name], ms=8,
                   label=f"{name}, {orient}")
        ax[0].plot(ref, d, mk, color=colours[name], ms=8, mfc="none")
    lim = (40, 240)
    ax[0].plot(lim, lim, "k:", lw=1)
    ax[0].fill_between(lim, np.subtract(lim, 15), np.add(lim, 15), color="0.9")
    ax[0].set(xlim=lim, ylim=lim, xlabel="Wu et al. BHCT [C]", ylabel="Model 1 BHCT [C]",
              title="Base case (filled: 1 mm coating; open: whole wall)")
    ax[0].legend(fontsize=7)
    for name, refs in FIG5.items():
        ax[1].plot(FIG5_T, refs, "--", color=colours[name])
        ax[1].plot(FIG5_T, [out[("fig5", name, G)][1] for G in FIG5_G], "-o",
                   color=colours[name], label=name)
        if name != "dual-wall":
            ax[1].plot(FIG5_T, [out[("fig5", name, G)][2] for G in FIG5_G], ":o",
                       color=colours[name], mfc="none")
    ax[1].set(xlabel="Reservoir temperature at 8.6 km [C]", ylabel="BHCT [C]",
              title="Vertical, 600 gal/min\n(solid: 1 mm coating; dotted: whole wall; dashed: Wu)")
    ax[1].legend(fontsize=8)
    for name, refs in FIG6.items():
        ax[2].plot(FIG6_GPM, refs, "--", color=colours[name])
        ax[2].plot(FIG6_GPM, [out[("fig6", name, q)][1] for q in FIG6_GPM], "-o",
                   color=colours[name], label=name)
        if name != "dual-wall":
            ax[2].plot(FIG6_GPM, [out[("fig6", name, q)][2] for q in FIG6_GPM], ":o",
                       color=colours[name], mfc="none")
    ax[2].set(xlabel="Flow [gal/min]", ylabel="BHCT [C]",
              title="Vertical, 400 C reservoir\n(solid: 1 mm coating; dotted: whole wall; dashed: Wu)")
    ax[2].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


GOLDEN = os.path.join(os.path.dirname(__file__), "..", "tests", "golden",
                      "benchmark_utaustin.json")


def golden_rows():
    """The base-case comparison and standpipe pressures, as pinned by the tests."""
    rows = {}
    for orient in ("vertical", "horizontal"):
        for name, pipe in PIPES.items():
            h = orient == "horizontal"
            rows[f"{orient}/{name}"] = dict(
                wu=FIG2[orient][name], mud=bhct(run(pipe, horizontal=h)),
                whole_wall=bhct(run(WHOLE_WALL[name], horizontal=h)))
    for name, pipe in PIPES.items():
        rows[f"spp/{name}"] = dict(spp_MPa=run(pipe)["hyd"]["spp"] / 1e6)
    return rows


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if "--write-golden" in sys.argv:
        import json
        with open(GOLDEN, "w") as fh:
            json.dump(golden_rows(), fh, indent=1)
        print("wrote", os.path.basename(GOLDEN))
        sys.exit()
    out = tables()
    path = os.path.join(os.path.dirname(__file__), "..", "figures", "benchmark_utaustin.png")
    figure(out, path)
    print(f"\nsaved {os.path.basename(path)}")
