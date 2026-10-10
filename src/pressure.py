"""
pressure.py
===========
The pressure on the wall through the drilling cycle, from one budget:

    P_w(z) = P_surface + integral_0^z rho_f(T, P) g dTVD + dP_friction(z) + dP_allowance

  * P_surface: the pressure held at the wellhead. The annulus is open while
    drilling, so it is zero (gauge) unless the return would boil at
    atmospheric pressure (backpressure()).
  * the column: the fluid's density integrated down the hole at the column's
    temperature in that state of the cycle. While drilling and through a
    connection the column has Model 1's annulus temperatures; at the end of a
    trip the open hole has the static fluid's temperature (site_evaluation).
  * dP_friction: the annular friction loss between the depth and the surface,
    while circulating.
  * dP_allowance: surge and swab during a trip, taken off the lower bound and
    added to the upper.

Pressures here are gauge, relative to the atmosphere, as the stresses and
pore pressures they are compared with are; the fluid's properties are taken at
the absolute pressure. An SG is the equivalent static density at the depth,
(P_static - P_surface) / (1000 g z).
"""
import numpy as np
from scipy.optimize import brentq
from iapws.iapws97 import _PSat_T, _TSat_P

import fluids
import geo_constants as C
import model1_coupled as m1

P_ATM = m1.P_SURF        # Pa absolute, at the wellhead
G = m1.G_GRAV

# Surge and swab allowance on each bound through a trip [SG]: an assumption,
# no cited value for deep wells has been found.
TRIP_ALLOWANCE_SG = 0.02


def sg_per_pa(z):
    """Pa per unit SG at true vertical depth z."""
    return C.MUD_SG_GRAD * z


def constant_column(z, grad, P_surface=0.0):
    """Static pressure [Pa gauge] at depth z of a column of constant gradient."""
    return P_surface + grad * np.asarray(z, float)


def column(tvd, T, fluid=None, P_surface=0.0):
    """Static pressure [Pa gauge] down a fluid column whose temperature is T
    [C] at true vertical depths tvd (increasing from the surface, tvd[0] = 0).
    fluid as for model1_coupled.fluid_props (None for water). The density
    depends on the pressure, so the integral is repeated with the pressure
    from the pass before; two passes settle it to well under 1 kPa."""
    tvd = np.asarray(tvd, float)
    T = np.clip(np.asarray(T, float), 1.0, 480.0)
    P = P_surface + C.HYDROSTATIC_GRAD * tvd
    for _ in range(3):
        rho = m1.fluid_props(T, P + P_ATM, fluid)[0]
        dP = 0.5 * (rho[1:] + rho[:-1]) * G * np.diff(tvd)
        P = P_surface + np.concatenate([[0.0], np.cumsum(dP)])
    return P


def circulating_temperatures(run):
    """True vertical depths and annulus temperatures of a Model 1 run, from
    the surface to the bit."""
    tvd = run["geometry"].survey.tvd(run["z"])
    return tvd, np.asarray(run["Tu"], float)


def annulus_friction(run, md):
    """Friction loss in the annulus [Pa] between measured depth md and the
    surface, from a run's hydraulics (model1_coupled.hydraulics)."""
    out = 0.0
    for p in run["hyd"]["pieces"]:
        if md <= p["md_top"]:
            break
        frac = min(1.0, (md - p["md_top"]) / (p["md_bottom"] - p["md_top"]))
        out += frac * float(p["dp_ann"])
    return out


def backpressure(T_return):
    """The wellhead pressure [Pa gauge] that keeps a return at T_return [C]
    liquid: zero while it is below the boiling point at atmospheric pressure,
    else its saturation pressure less the atmosphere."""
    if T_return + 273.15 < _TSat_P(P_ATM / 1e6):
        return 0.0
    return _PSat_T(T_return + 273.15) * 1e6 - P_ATM


def bottomhole_pressure(state, z, temperatures=None, fluid=None, run=None,
                        P_surface=0.0, allowance_SG=0.0, grad=None):
    """The pressure on the wall at true vertical depth z in one state of the
    cycle.

    state: "drilling" (circulating: the static column plus the annular
    friction from run), "connection" or "trip" (static). temperatures:
    (tvd, T) of the column in that state; grad instead gives a column of
    constant gradient [Pa/m] (the earlier model). allowance_SG: surge and swab.

    Returns P_static, friction and allowance [Pa], and the wall pressure on
    the lower bound P_lo (swab) and on the upper bound P_hi (surge); with
    SG, the static column as an equivalent density."""
    if grad is not None:
        P_static = float(constant_column(z, grad, P_surface))
    else:
        tvd, T = temperatures
        P_static = float(np.interp(z, tvd, column(tvd, T, fluid, P_surface)))
    friction = annulus_friction(run, z) if state == "drilling" else 0.0
    allowance = allowance_SG * sg_per_pa(z)
    return dict(state=state, z=z, P_surface=P_surface, P_static=P_static, friction=friction,
                allowance=allowance, P_lo=P_static + friction - allowance,
                P_hi=P_static + friction + allowance,
                SG=(P_static - P_surface) / sg_per_pa(z))


def surface_density_for(SG, z, temperatures, make_fluid, P_surface=0.0, hi=3000.0):
    """The surface density [kg/m3] of the fluid make_fluid(rho) whose static
    column in the given temperatures reaches the equivalent density SG at z;
    None if water's column already does."""
    tvd, T = temperatures
    f = lambda rho: (np.interp(z, tvd, column(tvd, T, make_fluid(rho), P_surface)) - P_surface) \
        / sg_per_pa(z) - SG
    lo = fluids.water().rho_surface
    if f(lo) >= 0.0:
        return None
    return brentq(f, lo, hi, xtol=0.01)
