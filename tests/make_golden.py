"""Regenerate the v1.0 four-site golden table.

Run from the repo root with the project venv:

    python tests/make_golden.py

This runs the FULL v1.0 pipeline (Models 1-5 per site), so it needs the IAPWS
water table (src/water_table.npz) and takes minutes. The output is the frozen
record of v1.0's comparative table; nothing in v1.1 should change it while the
v1.0 adapter is in use.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from comparative_sites import SITES, classify          # noqa: E402
from site_evaluation import evaluate                    # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden",
                   "v10_site_table.json")


def row(site):
    o = evaluate(site, v10=True)
    surv, stab, verdict = classify(o)
    b = o["breakout"]
    return {
        "name": site.name,
        "z_m": float(o["z"]),
        "T_rock_C": float(o["T_rock"]),
        "anisotropy": float(site.anisotropy),
        "Sv_MPa": float(o["Sv"]) / 1e6,
        "Shmin_MPa": float(o["Shmin"]) / 1e6,
        "SHmax_MPa": float(o["SHmax"]) / 1e6,
        "P_fluid_MPa": float(o["P_fluid"]) / 1e6,
        "m_min_kg_s": float(o["m_min"]),
        "T_bit_C": float(o["m1"]["T_bottom_delivered"]),
        "MW_prod": float(o["MW_prod"]),
        "ROP_m_hr": float(o["drill"]["ROP"]) * 3600.0,
        "rop_gain": float(o["rop_gain"]),
        "regime": o["drill"]["regime"],
        "grc_regime_cold": o["grc"]["reg_c"],
        "sigma_theta_MPa": float(b["sigma_theta"]) / 1e6,
        "mc_cold_MPa": float(b["mc_cold"]) / 1e6,
        "P_need_MPa": float(b["P_need"]) / 1e6,
        "overbalance_MPa": float(b["overbalance_MPa"]),
        "breaks": bool(b["breaks"]),
        "frac_limited": bool(b["frac_limited"]),
        "survives": bool(surv),
        "stability": stab,
        "verdict": verdict,
    }


if __name__ == "__main__":
    table = [row(s) for s in SITES]
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump({"note": "v1.0 pipeline output at commit d927bc4; "
                           "numpy/scipy versions affect the BVP outputs at the "
                           "1e-3 level, which is why the test tolerances are "
                           "relative and loose for Model 1 quantities.",
                   "sites": table}, f, indent=2)
        f.write("\n")
    for r in table:
        print(f"{r['name'][:34]:<36}{r['verdict']:>20}  "
              f"P_need {r['P_need_MPa']:6.1f} vs Shmin {r['Shmin_MPa']:6.1f} MPa")
    print("wrote", OUT)
