"""
A site's inputs, as a template.

Copy this file, replace every value with your site's, and give every value
its source. The values below are Soultz-sous-Forets, the reference case,
so that examples/new_site.py can check the template reproduces its row of
the site table. docs/new_site.md walks through the steps.

Units are SI throughout: metres, pascals, degrees Celsius. Each field's
allowed range is the one site_evaluation.validate_site checks
(site_evaluation.INPUT_RANGES); a value outside it is most often a unit
slip.

Sources: every input needs one (the data-basis tiers). The stress, strength,
pore pressure, temperature and well-check entries in data_basis decide
whether the site is evidence-based or speculative: stress "measured" or
"measured range", strength "calibrated", "measured (lab)" or "log-derived",
and a well check "calibrated" or "checked", or the site is reported as
speculative. The labels allowed are site_evaluation.BASES.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import geo_constants as C                            # noqa: E402
from site_evaluation import SiteProfile, StressProfile, urg_geotherm   # noqa: E402

# Valley & Evans (2007), Stanford Geothermal Workshop SGP-TR-183, valid 1.5 to
# 5.0 km: Sv = -1.30 + 25.50 z, Shmin = -1.78 + 14.06 z and
# -1.17 + 22.95 z <= SHmax < -1.37 + 26.78 z [MPa, z in km].
STRESS_SOURCE = "Valley & Evans (2007)"
SHMAX_RANGE = (("SHmax lower bound (0.90 Sv)", 22.95, -1.17),
               ("SHmax mid-range", 24.865, -1.27),
               ("SHmax upper bound (1.05 Sv)", 26.78, -1.37))


def make_site():
    return SiteProfile(
        # a name; the command line finds a site by a word of it
        name="Upper Rhine Graben / Soultz-sous-Forets (France)",
        # geotherm: any function from depth [m] to rock temperature [C], rising
        # with depth through the open hole. For measured points use
        # site_evaluation.layered_geotherm([(z, T), ...], gradient_below).
        geotherm=urg_geotherm,
        # the depth [m] of the target, 500 to 20,000, and its temperature
        # [C], 50 to 600, which the geotherm must give there within 5 C
        target_depth=10700.0,
        target_T=400.0,
        # the overburden gradient [Pa/m], 15,000 to 35,000 (rho_rock g); used
        # by the v1.0 adapter, the stress cases below set the stresses
        Sv_grad=2650.0 * 9.81,
        # the v1.0 adapter's ratios (Shmin/Sv, SHmax/Sv), dimensionless
        K0_min=0.54,
        SHmax_over_Sv=1.0,
        # the earlier model's fluid column [Pa/m], 8,000 to 20,000; v1.3.3 and
        # later integrate the column at its temperature
        rho_fluid_grad=C.HYDROSTATIC_GRAD,
        # rock thermal conductivity [W/m/K], 0.5 to 7; stiffness [Pa], 5e9 to
        # 1.5e11, which the thermal stress uses; the v1.0 strength [Pa]
        k_rock=2.9,
        E_rock=55e9,
        UCS=170e6,
        # the fluid's temperature at the pumps [C], 0 to 100
        T_inj=40.0,
        # stress cases: each a linear profile S = grad z + intercept [Pa/m, Pa],
        # with pore pressure Pp_grad (z - Pp_datum), the depths its data cover
        # (z_data, m), its source, and optionally its faulting regime
        # ("normal", "strike-slip" or "reverse"), whose stress order is checked
        stress_cases=tuple(
            StressProfile(label=lab, Sv_grad=25.50e3, Sv_int=-1.30e6,
                          Shmin_grad=14.06e3, Shmin_int=-1.78e6,
                          SHmax_grad=g * 1e3, SHmax_int=i * 1e6,
                          source=STRESS_SOURCE, z_data=(1500.0, 5000.0))
            for lab, g, i in SHMAX_RANGE),
        stress_basis="measured range",
        # strength cases: (label, UCS [Pa]), 10e6 to 500e6 each
        strength_cases=(("lab UCS 100", 100e6), ("lab UCS 115", 115e6), ("lab UCS 130", 130e6)),
        # the deepest temperature data [m]
        temperature_data_to=5000.0,
        # the friction coefficient of the frictional cap on SHmax below the
        # data and of the admissibility check, 0.2 to 1.5
        mu=C.FRICTION_MU,
        # the basis and source of each input
        data_basis={
            "stress": ("measured range", "Valley & Evans (2007), 1.5-5.0 km"),
            "strength": ("measured (lab)", "Valley & Evans (2007): 10 samples, 100-130 MPa; "
                         "not calibrated in situ"),
            "pore pressure": ("measured", "near-hydrostatic (Valley & Evans 2007)"),
            "temperature": ("extrapolated", "200 C at 5 km measured (GPK wells); 35 C/km "
                            "assumed below"),
            "well check": ("checked", "GPK3/GPK4 breakout onset and occurrence "
                           "(crosscheck_soultz.py)"),
        },
    )
