"""Depth-profile figures for Model 1 (coupled counterflow borehole).

The 450 C well at 12.4 km (35 C/km) on the site well design (Wu et al. 2025's
casing and string to the target), with friction heating and the loop held at
its operating pressure, the rock exposed for a year: conventional pipe,
NOV's TK-Drakon coated pipe and the dual-wall pipe of Xiao et al. (2022), at
the 30 kg/s that cleans an 8.67-inch hole.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geo_constants as C
import model1_coupled as m1
import well_geometry as wg

L = (C.TARGET_ROCK_TEMP - C.SURFACE_TEMP) / C.GEOTHERM_GRADIENT
M_DOT = 30.0
shoes = tuple(min(s, f * L) for s, f in zip(wg.WU_SHOES, (0.05, 0.25, 0.5)))

cases = [("conventional pipe", wg.CONVENTIONAL, "tab:red"),
         ("TK-Drakon coated pipe", wg.TK_DRAKON, "tab:orange"),
         ("dual-wall pipe", wg.DUAL_WALL, "tab:green")]


def run(pipe):
    return m1.solve(m_dot=M_DOT, geometry=wg.wu2025_well(L, pipe, shoes=shoes),
                    friction_heat=True, P_surface=m1.P_OP, Q_face=28000.0,
                    t_years=1.0, verbose=False)


if __name__ == "__main__":
    fig, axes = plt.subplots(1, len(cases), figsize=(15, 5.2), sharey=True)
    for ax, (label, pipe, colour) in zip(axes, cases):
        r = run(pipe)
        zkm = r["z"] / 1000
        ax.plot(r["Trock"], zkm, "k--", lw=1.3, label="rock")
        ax.plot(r["Td"], zkm, color="tab:blue", lw=2, label="down the pipe (to the bit)")
        ax.plot(r["Tu"], zkm, color="tab:red", lw=2, label="up the annulus")
        ax.axvline(C.BHA_SURVIVAL_TEMP, color="gray", ls=":", lw=1)
        ax.text(C.BHA_SURVIVAL_TEMP + 4, 1.0, "200 C\ntool limit", fontsize=8, color="gray")
        ax.invert_yaxis()
        ax.set_xlabel("Temperature [C]")
        ax.set_title(f"{label}, {M_DOT:.0f} kg/s\nbit {r['T_bottom_delivered']:.0f} C | "
                     f"return {r['T_return_surface']:.0f} C | {r['Q_product']/1e6:.1f} MW", fontsize=9)
        ax.grid(alpha=0.3)
        print(f"{label:<24} bit {r['T_bottom_delivered']:6.1f} C  return {r['T_return_surface']:6.1f} C"
              f"  heat {r['Q_product']/1e6:5.2f} MW  converged {r['success']}")
    axes[0].set_ylabel("Depth [km]")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Model 1: temperatures along a 12.4 km well to 450 C rock (35 C/km)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(__file__), "..", "figures", "model1_profiles.png"), dpi=130)
    print("saved model1_profiles.png")
