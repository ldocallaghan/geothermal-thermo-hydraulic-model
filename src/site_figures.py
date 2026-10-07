"""Dashboard figure for the site evaluation (Soultz / Upper Rhine Graben), v1.1:
measured stress profiles (Valley & Evans 2007), the breakout-width verdict
and the data basis."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geo_constants as C
import model1_coupled as m1
from comparative_sites import classify
from site_evaluation import SOULTZ, evaluate
import well_geometry as wg

o = evaluate(SOULTZ)
s = SOULTZ
z = o["z"]
surv, stab, verdict = classify(o)

# Model 1 profiles with the layered geotherm at a production flow
r = m1.solve(m_dot=10.0, T_inj=s.T_inj, pipe=wg.LEGACY_VACUUM, k_rock=s.k_rock,
             geotherm=s.temperature(), target_depth=z, Q_face=30000.0, verbose=False)
zkm = r["z"] / 1000

fig, ax = plt.subplots(1, 3, figsize=(15, 6))

# panel 1: layered geotherm + fluid profiles
ax[0].plot(r["Trock"], zkm, "k--", lw=1.5, label="rock (layered geotherm)")
ax[0].plot(r["Td"], zkm, color="tab:blue", lw=2, label="downcomer (to bit)")
ax[0].plot(r["Tu"], zkm, color="tab:red", lw=2, label="annulus (return)")
ax[0].axvline(200, color="gray", ls=":", lw=1); ax[0].text(205, 1, "200C\nceiling", fontsize=7, color="gray")
ax[0].axvline(374, color="purple", ls=":", lw=1); ax[0].text(330, 9.5, "supercrit\n374C", fontsize=7, color="purple")
ax[0].invert_yaxis(); ax[0].set_xlabel("Temperature [C]"); ax[0].set_ylabel("Depth [km]")
ax[0].set_title(f"Thermal: bit {r['T_bottom_delivered']:.0f}C, "
                f"return {r['T_return_surface']:.0f}C, {r['Q_product']/1e6:.1f} MW", fontsize=9)
ax[0].legend(fontsize=7, loc="lower left"); ax[0].grid(alpha=0.3)

# panel 2: measured stress profiles (Valley & Evans 2007), extrapolated below 5 km
zz = np.linspace(0, z, 50)
lo, mid, hi = s.stress_cases
ax[1].axhspan(lo.z_data[0] / 1000, lo.z_data[1] / 1000, color="0.9", zorder=0)
ax[1].text(5, lo.z_data[1] / 1000 - 0.2, "stress data 1.5-5 km", fontsize=7, color="dimgray")
ax[1].plot(lo.Sv(zz)/1e6, zz/1000, "k-", lw=2, label="Sv")
ax[1].fill_betweenx(zz/1000, lo.SHmax(zz)/1e6, hi.SHmax(zz)/1e6, color="tab:red", alpha=0.3,
                    label="SHmax, 0.90-1.05 Sv")
ax[1].plot(lo.Shmin(zz)/1e6, zz/1000, color="tab:green", lw=1.6, label="Shmin")
ax[1].plot(lo.Pp(zz)/1e6, zz/1000, color="tab:blue", lw=1.6, ls="--", label="Pp and mud (hydrostatic)")
ax[1].invert_yaxis(); ax[1].set_xlabel("Stress [MPa]"); ax[1].set_ylabel("Depth [km]")
w = o["stress"]["window"]
ax[1].set_title(f"Stress: worst-case breakout {w['width_hydro']:.0f} deg at hydrostatic\n"
                f"mud window {w['SG_lo']:.2f}-{w['SG_hi']:.2f} SG", fontsize=9)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].grid(alpha=0.3)

# panel 3: scorecard
ax[2].axis("off")
b = o["breakout"]; d = o["drill"]
rows = [
    ("SITE", SOULTZ.name.split("/")[1].split("(")[0].strip(), "k"),
    ("Target", f"{o['T_rock']:.0f} C @ {o['z']/1000:.1f} km", "k"),
    ("", "", "k"),
    ("Tool survival", f"OK ({o['m1']['T_bottom_delivered']:.0f} C bit)", "g"),
    ("Energy (10 kg/s)", f"{o['MW_prod']:.1f} MW_th", "k"),
    ("Drillability", f"{d['regime']}, {d['ROP']*3600:.1f} m/hr, {o['rop_gain']:.1f}x", "g"),
    ("Creep closure", f"controlled ({o['creep_hot']/max(o['creep_cold'],1e-30):.0e}x)", "g"),
    ("Isotropic stab.", f"{o['grc']['reg_c']}, {o['grc']['u_cold']*1000:.1f} mm", "g"),
    ("Breakout (width)", f"{min(c['window']['width_hydro'] for c in o['stress_cases']):.0f}-"
     f"{max(c['window']['width_hydro'] for c in o['stress_cases']):.0f} deg; {stab}", "orange"),
    ("Wall cooling", f"{o['T_wall']:.0f} C wall, {o['dsigma_T']/1e6:+.0f} MPa hoop", "k"),
    ("", "", "k"),
    ("VERDICT", verdict, "orange"),
    ("Data tier", f"{s.tier}; target {o['stress']['beyond_data']/1000:.1f} km below data", "k"),
    ("Strength", "lab UCS 100-130 MPa, not in-situ", "k"),
]
y = 0.95
for k, v, c in rows:
    ax[2].text(0.02, y, k, fontsize=9, fontweight="bold" if k in ("SITE","VERDICT") else "normal")
    ax[2].text(0.45, y, v, fontsize=9, color={"g":"green","orange":"darkorange","k":"black"}[c])
    y -= 0.075
ax[2].set_title("Scorecard", fontsize=10)

fig.suptitle("Site evaluation dashboard: Upper Rhine Graben / Soultz-sous-Forets", fontsize=12)
fig.tight_layout()
import os
_figdir = os.path.join(os.path.dirname(__file__), "..", "figures")
fig.savefig(os.path.join(_figdir, "site_dashboard.png"), dpi=130)
print("saved site_dashboard.png")
