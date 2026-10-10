"""
fluids.py
=========
The drilling fluid: water, or a water-based mud weighted with barite.

Density. The mud is water and barite. The barite's volume fraction phi is set
by the mud's density at surface conditions (T_REF, atmospheric):

    phi = (rho_surface - rho_w(T_REF)) / (rho_barite(T_REF) - rho_w(T_REF))

and down the hole the water follows IAPWS-95 at the local temperature and
pressure while the barite, incompressible, expands only thermally:

    rho(T, P) = phi rho_barite(T) + (1 - phi) rho_w(T, P)

so a weighted mud expands less with temperature than water does. Other
solids (clay, polymers) are left out: at these weights barite is most of the
solids by volume.

Rheology. Bingham plastic, with the plastic viscosity and yield point of a
barite-weighted water-based mud measured by Anawe and Folayan (2018) at
80 F (PV 27 cP, YP 28 lbf/100 ft2). They give no density, so the rheology
does not vary with it. At the site flows the annulus is turbulent, where the
loss follows the plastic viscosity and not the yield point, so the
sensitivity halves the plastic viscosity. Between 80 and 200 F their plastic
viscosity falls from 27 to 9 cP, as water's viscosity does to within about
7%, so the plastic viscosity is taken as the water's at the local temperature
and pressure times its ratio at T_REF. Their yield point halves over the same
range; with no law beyond 200 F it is held at its 80 F value, which gives
the larger annular friction, the conservative side for the fracture bound.

Thermal properties. Heat capacity by mass fraction; conductivity by the
Maxwell relation for dispersed spheres.

Hydraulics (bingham_gradient): the laminar pressure loss of a Bingham plastic
in pipe and in the narrow-slot annulus, with the apparent-viscosity Reynolds
number for the change to turbulence, after Bourgoyne et al. (1986), ch. 4;
in turbulence the plastic viscosity in the Blasius friction factor on the
hydraulic diameter, the correlation model1_coupled uses for water.

Sources:
  Anawe, P.A.L. and Folayan, J.A. (2018). Data analyses on temperature-
    dependent behaviour of water based drilling fluid rheological models.
    Data in Brief 21, 289-298. doi:10.1016/j.dib.2018.09.100
  API Specification 13A: drilling-grade barite, specific gravity 4.20.
  Robie, R.A. and Hemingway, B.S. (1995). Thermodynamic properties of
    minerals. USGS Bulletin 2131: barite heat capacity at 298 K.
  Bourgoyne, A.T., Millheim, K.K., Chenevert, M.E. and Young, F.S. (1986).
    Applied Drilling Engineering. SPE Textbook Series 2, ch. 4.
"""
import numpy as np

import water_table as wt

T_REF = 20.0                  # C, the surface conditions a mud weight is quoted at
P_ATM = 1.0e5                 # Pa

BARITE_RHO = 4200.0           # kg/m3 at T_REF, API 13A drilling-grade barite
BARITE_CP = 436.0             # J/kg/K, 101.8 J/mol/K / 0.2334 kg/mol (Robie and Hemingway 1995)
BARITE_K = 1.5                # W/m/K; to be confirmed against a cited value
BARITE_BETA = 5.0e-5          # 1/K, volumetric thermal expansion; to be confirmed

PV_REF = 27.0e-3              # Pa s at 80 F (26.7 C), Anawe and Folayan (2018)
PV_REF_T = 26.7               # C
LBF_100FT2 = 0.4788026        # Pa
YP_REF = 28.0 * LBF_100FT2    # Pa (28 lbf/100 ft2 at 80 F), held in temperature

# conventional water-based muds start to break down above about 300 F
# (secondary source; to be confirmed against a drilling-fluids text)
MUD_T_LIMIT = 150.0           # C

RE_CRITICAL = 2100.0          # apparent-viscosity Reynolds number (Bourgoyne et al. 1986)


def _water(T, P):
    return wt.props(np.asarray(T, float), np.asarray(P, float))


def _rho_water_ref():
    return float(np.ravel(_water([T_REF], [P_ATM])[0])[0])


class Mud:
    """A barite-weighted water-based mud of density rho_surface [kg/m3] at
    T_REF and atmospheric pressure. With rho_surface at water's density and
    yield_point 0 it is water (pv_ref None takes water's own viscosity)."""

    def __init__(self, rho_surface, pv_ref=PV_REF, yield_point=YP_REF, pv_ref_T=PV_REF_T):
        rho_w = _rho_water_ref()
        if rho_surface < rho_w - 1e-6:
            raise ValueError("a barite mud is no lighter than water")
        self.rho_surface = float(rho_surface)
        self.phi = (rho_surface - rho_w) / (BARITE_RHO - rho_w)
        self.yield_point = float(yield_point)
        mu_w_ref = float(np.ravel(_water([pv_ref_T], [P_ATM])[2])[0])
        self.mu_ratio = 1.0 if pv_ref is None else pv_ref / mu_w_ref

    @property
    def sg(self):
        return self.rho_surface / 1000.0

    def props(self, T, P):
        """rho, cp, mu (the plastic viscosity), k at temperature T [C] and
        pressure P [Pa absolute], as model1_coupled.fluid_props."""
        rho_w, cp_w, mu_w, k_w = _water(T, P)
        T = np.asarray(T, float)
        rho_b = BARITE_RHO / (1.0 + BARITE_BETA * (T - T_REF))
        phi = self.phi
        rho = phi * rho_b + (1.0 - phi) * rho_w
        x_b = phi * rho_b / rho                   # mass fraction of barite
        cp = x_b * BARITE_CP + (1.0 - x_b) * cp_w
        # Maxwell: spheres of conductivity k_b at volume fraction phi in water
        k_b = BARITE_K
        k = k_w * (k_b + 2 * k_w + 2 * phi * (k_b - k_w)) / (k_b + 2 * k_w - phi * (k_b - k_w))
        return rho, cp, mu_w * self.mu_ratio, k


def water():
    """Water as a Mud: no barite, no yield point, water's viscosity."""
    return Mud(_rho_water_ref(), pv_ref=None, yield_point=0.0)


def bingham_gradient(rho, mu_p, tau_y, v, D, annulus):
    """Frictional pressure gradient [Pa/m] of a Bingham plastic at mean
    velocity v in a pipe of diameter D, or an annulus of hydraulic diameter D
    (outer less inner diameter). Laminar, from the plug-flow limits of
    Bourgoyne et al. (1986), ch. 4:

        pipe:    32 mu_p v / D^2 + 16 tau_y / (3 D)
        annulus: 48 mu_p v / D^2 +  6 tau_y / D        (narrow slot)

    (in field units mu_p v / (1500 d^2) + tau_y / (225 d) and
    mu_p v / (1000 (d2 - d1)^2) + tau_y / (200 (d2 - d1)), psi/ft). The flow
    is turbulent when the Reynolds number on the apparent viscosity, the
    Newtonian viscosity giving the same laminar loss, exceeds RE_CRITICAL;
    then the Blasius friction factor on the plastic viscosity."""
    rho, mu_p, v, D = (np.asarray(x, float) for x in (rho, mu_p, v, D))
    v = np.maximum(np.abs(v), 1e-9)
    if annulus:
        lam = 48.0 * mu_p * v / D ** 2 + 6.0 * tau_y / D
        mu_a = mu_p + tau_y * D / (8.0 * v)
    else:
        lam = 32.0 * mu_p * v / D ** 2 + 16.0 * tau_y / (3.0 * D)
        mu_a = mu_p + tau_y * D / (6.0 * v)
    Re_a = rho * v * D / mu_a
    turb = 0.316 * (rho * v * D / mu_p) ** -0.25 * rho * v ** 2 / (2.0 * D)
    return np.where(Re_a > RE_CRITICAL, turb, lam), Re_a
