"""
site_evaluation.py
=================
Site-evaluation tool: wires Models 1-5 into a single verdict for a candidate
site. Feed in a site's temperature profile, in-situ stress and rock properties;
get back drillability, tool survival, hole stability and energy, each with its
operating envelope.

Hole stability is judged on breakout width: an effective-stress Kirsch check
with the thermal stress of the cooled wall, run over every stress and strength
case the site's data support, with a mud window bounded by Shmin. Each site
also records what every input rests on (data_basis), which decides whether its
verdict is evidence-based or speculative.

The site defined here is the Upper Rhine Graben at Soultz-sous-Forets (France),
the best-characterised deep-geothermal site in Europe (GPK-1..4 wells to
~5 km). Temperature: ~90 C/km to 1.4 km, then near-isothermal convective
granite to 200 C at 5 km (Genter et al.; MDPI Geosciences 2020). Stress and
rock strength: Valley & Evans (2007); see data/sites/stress_sources.md.
"""
from dataclasses import dataclass, field
import numpy as np
import geo_constants as C
import model1_coupled as m1
import model2_spallation as m2
import model3_optimiser as m3
import model4_hole_stability as m4
import model5_convergence_confinement as m5
import wall_thermal
import well_geometry as wg

YEAR = 3.156e7


# --------------------------------------------------------------- site profile
def urg_geotherm(z):
    """Layered Upper Rhine Graben geotherm [degC], vectorized. z in metres.
    0-1400 m: ~90 C/km (sediments); 1400-5000 m: convective granite -> 200 C at
    5 km; >5000 m: ~35 C/km conductive (regional high heat flow)."""
    z = np.asarray(z, dtype=float)
    return np.select(
        [z <= 1400, z <= 5000],
        [12.0 + 0.090 * z,
         138.0 + 0.0172 * (z - 1400.0)],
        default=200.0 + 0.035 * np.clip(z - 5000.0, 0, None))


# ------------------------------------------------------------ stress profiles
@dataclass(frozen=True)
class StressProfile:
    """Linear in-situ stress with depth: S(z) = grad * z + intercept [Pa, z in m].

    Pore pressure is Pp_grad * (z - Pp_datum) below a static fluid level at
    Pp_datum [m], zero above it. z_data is the depth range [m] the magnitude data
    actually cover, or None when the profile rests on regime bounds alone.
    """
    label: str
    Sv_grad: float
    Shmin_grad: float
    SHmax_grad: float
    Sv_int: float = 0.0
    Shmin_int: float = 0.0
    SHmax_int: float = 0.0
    Pp_grad: float = C.HYDROSTATIC_GRAD
    Pp_datum: float = 0.0
    source: str = ""
    z_data: tuple = None
    # the friction coefficient the case was derived with, if it has its own;
    # its admissibility and cap are then checked at that value
    mu: float = None

    @classmethod
    def from_ratios(cls, Sv_grad, K0_min, SHmax_over_Sv, label="v1.0 ratios",
                    Pp_grad=C.HYDROSTATIC_GRAD, **kw):
        """v1.0's constant ratios against a linear overburden (no intercepts)."""
        return cls(label=label, Sv_grad=Sv_grad, Shmin_grad=K0_min * Sv_grad,
                   SHmax_grad=SHmax_over_Sv * Sv_grad, Pp_grad=Pp_grad, **kw)

    def Sv(self, z):
        return self.Sv_grad * z + self.Sv_int

    def Shmin(self, z):
        return self.Shmin_grad * z + self.Shmin_int

    def SHmax(self, z):
        return self.SHmax_grad * z + self.SHmax_int

    def Pp(self, z):
        return self.Pp_grad * np.maximum(np.asarray(z, dtype=float) - self.Pp_datum, 0.0)

    def beyond_data(self, z):
        """Metres by which z lies below the deepest data; None with no data."""
        return None if self.z_data is None else max(z - self.z_data[1], 0.0)

    def cap_depth(self, mu):
        """Depth [m] below which SHmax exceeds the frictional limit at mu, or None.

        Below the static fluid level every term is linear in z, so the crossing
        of SHmax - Pp = R (Shmin - Pp) is solved in closed form.
        """
        R = m5.frictional_cap(mu)
        a = self.SHmax_grad - self.Pp_grad - R * (self.Shmin_grad - self.Pp_grad)
        b = (self.SHmax_int - R * self.Shmin_int
             - (R - 1) * self.Pp_grad * self.Pp_datum)
        a_tol = m5.CAP_RTOL * (abs(self.SHmax_grad) + R * abs(self.Shmin_grad))
        if abs(a) <= a_tol:              # parallel to the cap: on it, or over it everywhere
            return 0.0 if b > 1.0 else None
        if a < 0:
            return None                  # moves away from the cap with depth
        return float(max(-b / a, self.Pp_datum))


def transitional_bounds(Sv_grad, source, Sv_int=0.0, Pp_grad=C.HYDROSTATIC_GRAD,
                        mu=C.FRICTION_MU, fractions=(0.0, 0.25, 0.5, 0.75, 1.0)):
    """Profiles for a normal/strike-slip transition regime (SHmax ~ Sv), for a
    site where the faulting regime is known but no magnitudes are measured.

    SHmax = Sv; Shmin runs from the frictional floor Pp + (Sv - Pp)/R(mu), at
    fraction 0, to Sv at fraction 1. Hydrostatic Pp from surface keeps every
    profile linear. The floor case is the most anisotropic one the regime allows.
    """
    R = m5.frictional_cap(mu)
    floor_g = Pp_grad + (Sv_grad - Pp_grad) / R
    floor_i = Sv_int / R
    return tuple(
        StressProfile(label=f"Shmin {f:.0%} of floor-to-Sv", Sv_grad=Sv_grad,
                      Sv_int=Sv_int, SHmax_grad=Sv_grad, SHmax_int=Sv_int,
                      Shmin_grad=floor_g + f * (Sv_grad - floor_g),
                      Shmin_int=floor_i + f * (Sv_int - floor_i),
                      Pp_grad=Pp_grad, source=source, z_data=None)
        for f in fractions)


@dataclass
class SiteProfile:
    name: str
    geotherm: callable          # z[m] -> T_rock[degC]
    target_depth: float         # m
    target_T: float             # degC (consistency check vs geotherm)
    Sv_grad: float              # Pa/m, vertical (overburden) stress gradient
    K0_min: float               # Shmin / Sv  (low => extensional, favourable)
    SHmax_over_Sv: float        # SHmax / Sv
    rho_fluid_grad: float       # Pa/m, in-hole fluid pressure gradient
    k_rock: float; E_rock: float; UCS: float
    T_inj: float = 40.0
    # The stress cases the site's data support: one for a fully measured
    # profile, several where a magnitude is a range or the stresses rest on
    # regime bounds. Empty means v1.0's ratios, via from_ratios.
    stress_cases: tuple = ()
    stress_basis: str = "v1.0 ratios"     # "measured", "measured range" or "regime bounds"
    # Rock-strength cases, (label, UCS at 25 C [Pa]), each run against every
    # stress case: a calibrated range, or a measured lab range. Empty means
    # the single value UCS.
    strength_cases: tuple = ()
    # Data basis: data_basis maps each input -- "stress",
    # "strength", "pore pressure", "temperature", "well check" -- to
    # (basis, source), basis one of BASES. geotherm_v11 / target_depth_v11
    # replace the v1.0 temperature profile where measurements contradict it
    # (the v1.0 adapter keeps the old ones); temperature_data_to is the
    # deepest measured temperature [m], None if there is none.
    data_basis: dict = field(default_factory=dict)
    geotherm_v11: callable = None
    target_depth_v11: float = None
    temperature_data_to: float = None
    # friction coefficient for frictional admissibility and the cap: the
    # United Downs value unless the site's stresses were derived with another
    mu: float = C.FRICTION_MU

    def __post_init__(self):
        if not self.strength_cases:
            self.strength_cases = (("site UCS", self.UCS),)
        if not self.stress_cases:
            self.stress_cases = (StressProfile.from_ratios(
                self.Sv_grad, self.K0_min, self.SHmax_over_Sv,
                Pp_grad=self.rho_fluid_grad, source="v1.0 site ratios"),)

    @property
    def anisotropy(self):
        return self.SHmax_over_Sv / self.K0_min   # SHmax / Shmin

    def temperature(self, v10=False):
        """The temperature profile a run uses: v1.0's under the adapter."""
        return self.geotherm if v10 or self.geotherm_v11 is None else self.geotherm_v11

    def target(self, v10=False):
        return (self.target_depth if v10 or self.target_depth_v11 is None
                else self.target_depth_v11)

    @property
    def tier(self):
        """"evidence-based" only if stress magnitudes and rock strength are
        both measured or calibrated at the site and the stability model has
        been checked against a real well there; otherwise "speculative"."""
        b = {k: v[0] for k, v in self.data_basis.items()}
        ok = (b.get("stress") in ("measured", "measured range")
              and b.get("strength") in ("calibrated", "measured (lab)", "log-derived")
              and b.get("well check") in ("calibrated", "checked"))
        return "evidence-based" if ok else "speculative"

    def missing(self):
        """Inputs that keep a site speculative, as (input, basis)."""
        need = {"stress": ("measured", "measured range"),
                "strength": ("calibrated", "measured (lab)", "log-derived"),
                "well check": ("calibrated", "checked")}
        return [(k, self.data_basis.get(k, ("none", ""))[0]) for k, good in need.items()
                if self.data_basis.get(k, ("none", ""))[0] not in good]


# Data-basis labels, strongest first
BASES = ("calibrated", "measured", "measured range", "measured (lab)", "log-derived",
         "checked", "regime bounds", "extrapolated", "regional", "assumed",
         "not reproduced", "unsourced", "none")


SOULTZ = SiteProfile(
    name="Upper Rhine Graben / Soultz-sous-Forets (France)",
    geotherm=urg_geotherm,
    target_depth=10700.0,       # ~400 C per layered geotherm
    target_T=400.0,
    Sv_grad=2650.0 * 9.81,      # URG granite ~2650 kg/m3 -> 26.0 kPa/m
    K0_min=0.54,                # Shmin/Sv ~0.54 (graben extension) -- measured
    SHmax_over_Sv=1.0,          # SHmax ~ Sv (strike-slip/normal transition)
    rho_fluid_grad=C.HYDROSTATIC_GRAD,
    k_rock=2.9,                 # URG biotite granite ~2.5-3.2 W/m/K
    E_rock=55e9,                # granite ~50-60 GPa
    UCS=170e6,                  # v1.0: "crystalline basement ~150-200 MPa", unsourced
    # Valley & Evans (2007), Stanford Geothermal Workshop SGP-TR-183, valid
    # 1.5-5.0 km: Sv = -1.30 + 25.50z, Shmin = -1.78 + 14.06z,
    # -1.17 + 22.95z <= SHmax < -1.37 + 26.78z  [MPa, z in km], i.e.
    # 0.90 Sv <= SHmax <= 1.05 Sv. Near-hydrostatic pore pressure.
    stress_cases=tuple(
        StressProfile(label=lab, Sv_grad=25.50e3, Sv_int=-1.30e6,
                      Shmin_grad=14.06e3, Shmin_int=-1.78e6,
                      SHmax_grad=g * 1e3, SHmax_int=i * 1e6,
                      source="Valley & Evans (2007)", z_data=(1500.0, 5000.0))
        for lab, g, i in (("SHmax lower bound (0.90 Sv)", 22.95, -1.17),
                          ("SHmax mid-range", 24.865, -1.27),
                          ("SHmax upper bound (1.05 Sv)", 26.78, -1.37))),
    stress_basis="measured range",
    # Lab UCS of ten samples of unaltered Soultz granite, 100-130 MPa
    # (Valley & Evans 2007), in place of the unsourced 170 MPa above
    strength_cases=(("lab UCS 100", 100e6), ("lab UCS 115", 115e6), ("lab UCS 130", 130e6)),
    temperature_data_to=5000.0,
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


def thermal_hoop_stress(T_rock, T_wall):
    """Thermal hoop stress [Pa] at a wall cooled from T_rock to T_wall.

    dsigma_T = -E alpha / (1 - nu) * (T_rock - T_wall): the fully constrained
    thermoelastic stress Model 2 applies at the cutting face, with E and alpha
    from Model 2 at the rock temperature, as Model 2 evaluates them. Negative
    (less compressive) when the wall is cooled.
    """
    E, al = m2.E_of_T(T_rock), m2.alpha_of_T(T_rock)
    return -E * al / (1.0 - C.NU_ROCK) * (T_rock - T_wall)


def stability_inputs(site, z, T_rock, v10=False, mu=None, UCS=None,
                     T_wall=None, thermal=True, W_max=C.BREAKOUT_W_MAX_DEG, ecd_SG=0.0,
                     Pw_add=0.0):
    """Stress state, frictional admissibility, breakout and mud window for
    every stress and strength case at depth z.

    v10=True reproduces v1.0: the site ratios and the total-stress breakout
    check, nothing else. Otherwise each of the site's stress cases is
    evaluated with its own pore pressure, the breakout check is in effective
    stress, and SHmax is capped at the frictional limit below the depth the
    case's data cover. The reference case -- the one a site's table row
    reports -- is the one with the worst verdict, then the narrowest mud
    window (closest to NO-GO), then the widest breakout at hydrostatic mud.
    Under v10 it is the one needing the most mud to suppress breakout outright.

    Strength: UCS if given, else each of the site's strength cases, run
    against every stress case; the v1.0 adapter uses the global C.UCS.

    T_wall is the wall temperature for strength and, with thermal=True, for the
    thermal hoop stress; evaluate() passes Model 1's circulating bottom-hole
    temperature. None means v1.0's min(T_rock, 200 C). The v1.0 adapter
    ignores both. Needs no Model 1-4 calls.
    """
    if mu is None:
        mu = site.mu
    # Pw_add: pressure added to the fluid on both bounds (circulating friction
    # while drilling, when the breakout bound is checked circulating too)
    P_mud = site.rho_fluid_grad * z + Pw_add
    if T_wall is None:
        T_wall = min(T_rock, T_WALL_CAP)
    dsigma_T = thermal_hoop_stress(T_rock, T_wall) if thermal and not v10 else 0.0
    if v10:
        Sv, Shmin, SHmax, _ = stresses_v10(site, z)
        profiles = [(site.stress_cases[0], Sv, Shmin, SHmax, C.HYDROSTATIC_GRAD * z)]
    else:
        profiles = [(p, p.Sv(z), p.Shmin(z), p.SHmax(z), float(p.Pp(z)))
                    for p in site.stress_cases]
    strengths = ([("", UCS)] if v10 or UCS is not None else
                 list(site.strength_cases))
    cases = []
    for (p, Sv, Shmin, SHmax, Pp), (s_label, UCS) in (
            (pr, st) for pr in profiles for st in strengths):
        window = None
        beyond = None if v10 else p.beyond_data(z)
        mu_c = p.mu if p.mu is not None and not v10 else mu
        SHmax_lim = m5.SHmax_frictional_limit(Shmin, Pp, mu_c)
        cap_binds = bool(beyond is not None and beyond > 0
                         and SHmax > SHmax_lim * (1 + m5.CAP_RTOL))
        if cap_binds:
            SHmax = SHmax_lim
        if not v10:
            window = mud_window(Shmin, SHmax, Pp, z, P_mud, T_rock, T_wall=T_wall,
                                UCS=UCS, dsigma_T=dsigma_T, W_max=W_max, ecd_SG=ecd_SG)
        cases.append(dict(window=window, T_wall=T_wall, dsigma_T=dsigma_T,
            UCS=C.UCS if UCS is None else UCS, strength_label=s_label,
            label=p.label + (f", {s_label}" if len(strengths) > 1 else ""),
            profile=p, Sv=Sv, Shmin=Shmin, SHmax=SHmax, Pp=Pp,
            beyond_data=beyond, cap_binds=cap_binds, cap_depth=p.cap_depth(mu_c),
            admissible=m5.stress_admissible(Sv, Shmin, SHmax, Pp, mu_c),
            breakout=(breakout_v10(Shmin, SHmax, P_mud, T_rock, UCS=UCS) if v10 else
                      breakout_eff(Shmin, SHmax, P_mud, Pp, T_rock, T_wall=T_wall,
                                   UCS=UCS, dsigma_T=dsigma_T))))
    if v10:
        ref = max(cases, key=lambda c: c["breakout"]["P_need"])
    else:   # worst verdict, then narrowest window, then widest breakout at hydro
        ref = max(cases, key=lambda c: (VERDICT_ORDER.index(c["window"]["verdict"]),
                                        c["window"]["Pw_lo"] - c["window"]["Pw_hi"],
                                        c["window"]["width_hydro"]))
    return ref, cases


def stresses_v10(site, z):
    """v1.0 stress state at depth z [m]: (Sv, Shmin, SHmax, P_fluid) in Pa.

    Constant ratios against a linear overburden gradient, so SHmax/Shmin is the
    same at every depth. v1.1 replaces these with measured, depth-dependent
    profiles (StressProfile).
    """
    Sv = site.Sv_grad * z
    return Sv, site.K0_min * Sv, site.SHmax_over_Sv * Sv, site.rho_fluid_grad * z


# ------------------------------------------------- v1.0 breakout arithmetic
# Wall temperature cap [degC]: the v1.0 checks assume the circulating fluid
# holds the wall at or below this, so strength is evaluated at the lower of the
# rock temperature and this cap.
T_WALL_CAP = 200.0


def breakout_v10(Shmin, SHmax, P_fluid, T_rock, T_wall=None, UCS=None):
    """v1.0 anisotropic breakout check at the Shmin azimuth of a vertical hole.

    Kirsch hoop stress sigma_theta = 3*SHmax - Shmin - P_i, failed against a
    Mohr-Coulomb limit KMC*P_i + sigma_cm(T). NOTE (and the reason v1.1 exists):
    pore pressure never enters -- the in-hole fluid pressure P_i is used both as
    the hoop-relieving term and as the confining term on the strength side.

    Returns the dict v1.0 reports as out["breakout"]. Extracted verbatim from
    evaluate() so the v1.0 numbers can be regression-tested without building the
    water table or running Models 1-4. UCS=None is v1.0's global C.UCS.
    """
    if T_wall is None:
        T_wall = min(T_rock, T_WALL_CAP)
    sth = 3 * SHmax - Shmin - P_fluid
    mc_cold = m5.KMC * P_fluid + m5.sigma_cm(T_wall, UCS)
    mc_hot = m5.KMC * P_fluid + m5.sigma_cm(T_rock, UCS)
    # mud weight needed to suppress breakout even when cold
    P_need = (3 * SHmax - Shmin - m5.sigma_cm(T_wall, UCS)) / (1 + m5.KMC)
    # If the mud weight needed to stop breakout exceeds Shmin, you would
    # hydraulically fracture the formation (lose returns) -> NOT mud-controllable.
    return dict(sigma_theta=sth, mc_hot=mc_hot, mc_cold=mc_cold,
                P_need=P_need, over_hydro=P_need / P_fluid,
                overbalance_MPa=(P_need - P_fluid) / 1e6,
                breaks=sth > mc_cold, frac_limited=P_need > Shmin)


def breakout_eff(Shmin, SHmax, Pw, Pp, T_rock, T_wall=None, UCS=None, dsigma_T=0.0):
    """Effective-stress breakout check at the Shmin azimuth.

    Fails when sigma_theta - Pp > sigma_cm(T_wall, UCS) + KMC (Pw - Pp): the
    Kirsch hoop stress from model5.hoop_stress, Mohr-Coulomb in effective
    stress, with no filter-cake credit beyond Pw - Pp. With Pp = 0 and
    dsigma_T = 0 it is breakout_v10 exactly. Same keys as breakout_v10, with
    sigma_theta the total hoop stress and mc_cold / mc_hot the total-stress
    limits (Pp plus the effective strength), so the two stay comparable.
    """
    if T_wall is None:
        T_wall = min(T_rock, T_WALL_CAP)
    K = m5.KMC
    sth = float(m5.hoop_stress(90.0, SHmax, Shmin, Pw, dsigma_T))
    mc_cold = Pp + m5.sigma_cm(T_wall, UCS) + K * (Pw - Pp)
    mc_hot = Pp + m5.sigma_cm(T_rock, UCS) + K * (Pw - Pp)
    # the Pw at which sigma_theta - Pp equals the cold limit
    P_need = (3 * SHmax - Shmin + dsigma_T - m5.sigma_cm(T_wall, UCS)
              + (K - 1) * Pp) / (1 + K)
    return dict(sigma_theta=sth, mc_hot=mc_hot, mc_cold=mc_cold,
                P_need=P_need, over_hydro=P_need / Pw,
                overbalance_MPa=(P_need - Pw) / 1e6,
                breaks=bool(sth > mc_cold), frac_limited=bool(P_need > Shmin))


VERDICT_ORDER = ("GO", "CONDITIONAL", "NO-GO")


def describe_width(width_deg):
    """A breakout width for printing. 180 deg means the wall fails all the way
    round (with equal horizontal stresses, all or nothing), which is general
    yielding of the wall rather than a wide breakout."""
    return "wall yields all round" if width_deg >= 180.0 else f"{width_deg:.0f} deg breakout"


def mud_window(Shmin, SHmax, Pp, z, Pw_hydro, T_rock, T_wall=None, UCS=None,
               dsigma_T=0.0, W_max=C.BREAKOUT_W_MAX_DEG, margin_SG=C.MUD_MARGIN_SG,
               T0=C.WALL_T0, ecd_SG=0.0):
    """Mud window and drillability verdict, judged on breakout width.

    The two bounds are checked under the conditions that set them. Breakout,
    the lower bound, at static mud weight: with the pumps off at a connection,
    for minutes while the wall is still cool, it sees only the static column.
    A trip takes the hole out of circulation for a day or more at these
    depths, and the wall reheats towards rock temperature; that is closer to
    the uncooled sensitivity, and reheating is not modelled. Fracture, the
    upper bound, at circulating pressure: the static mud plus ecd_SG, the
    annular friction loss as an equivalent density, must stay below Shmin less
    the margin. ecd_SG = 0 is the static check on both.

    Lower bound: the lightest mud keeping the breakout lobe at or below W_max,
    and never below pore pressure (no underbalanced drilling). Upper bound:
    Shmin less margin_SG of mud (lost returns). Verdict: GO if the lobe is
    within W_max at hydrostatic mud, CONDITIONAL if a mud weight inside the
    window gets it there, NO-GO if none does. Tensile-fracture initiation at
    the wall (effective 3 Shmin - SHmax - Pw - Pp + dsigma_T < -T0) is
    reported and does not bound the window: UD-1 logged tensile fractures and
    still reached TD, so initiation alone doesn't stop a well.
    """
    if T_wall is None:
        T_wall = min(T_rock, T_WALL_CAP)
    sg = C.MUD_SG_GRAD * z                       # Pa per unit SG at this depth
    width = lambda Pw: m5.breakout_width(SHmax, Shmin, Pw, Pp, T_wall, UCS, dsigma_T)
    Pw_lo = max(m5.mud_for_width(SHmax, Shmin, Pp, T_wall, W_max, UCS, dsigma_T), Pp)
    Pw_hi = Shmin - margin_SG * sg - ecd_SG * sg
    w_hydro = width(Pw_hydro)
    if w_hydro <= W_max and Pw_hydro <= Pw_hi:
        verdict, Pw_mud = "GO", Pw_hydro
    elif Pw_lo <= Pw_hi:
        verdict, Pw_mud = "CONDITIONAL", Pw_lo
    else:
        verdict, Pw_mud = "NO-GO", None
    # tensile initiation once Pw exceeds this
    Pw_tensile = 3 * Shmin - SHmax + dsigma_T - Pp + T0
    return dict(
        verdict=verdict, W_max=W_max, T_wall=T_wall,
        Pw_hydro=Pw_hydro, Pw_lo=Pw_lo, Pw_hi=Pw_hi, open=bool(Pw_lo <= Pw_hi),
        SG_hydro=Pw_hydro / sg, SG_lo=Pw_lo / sg, SG_hi=Pw_hi / sg, ecd_SG=ecd_SG,
        SG_hi_static=(Shmin - margin_SG * sg) / sg,
        width_hydro=w_hydro, width_hi=width(Pw_hi),
        Pw_mud=Pw_mud,
        overbalance_MPa=None if Pw_mud is None else (Pw_mud - Pw_hydro) / 1e6,
        overbalance_SG=None if Pw_mud is None else (Pw_mud - Pw_hydro) / sg,
        Pw_tensile=Pw_tensile,
        tensile_at_hydro=bool(Pw_hydro > Pw_tensile),
        tensile_at_mud=None if Pw_mud is None else bool(Pw_mud > Pw_tensile))


def worst_verdict(verdicts):
    return max(verdicts, key=VERDICT_ORDER.index)


# ------------------------------------------------------------- the evaluation
# ------------------------------------------------------------ circulation
# The site well: Wu et al. (2025)'s deep-well design (well_geometry.wu2025_well)
# to the target, with the casing shoes moved up for shallow targets. The pipe
# is the best-insulated one that is commercially available, trying the
# commercial types from best insulated to worst. The dual-wall pipe's figures
# come from a modelling study and no such drill pipe is on sale
# (data/materials.md), so it is a sensitivity, as is vacuum-insulated tubing,
# which is not a drill pipe.
SITE_PIPES = (wg.TK_DRAKON, wg.INTERNALLY_COATED, wg.EXTERNALLY_COATED, wg.CONVENTIONAL)
COMMERCIAL_PIPES = SITE_PIPES
SITE_FLOWS = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100)   # kg/s
# The highest flow FORGE 16B was drilled at (rig_hydraulics.csv).
FORGE_MAX_FLOW_GPM = 700.0
# What the rig's pumps can deliver at a standpipe pressure: three triplex mud
# pumps (FORGE 16B's rig had at least three, from its daily reports and Pason
# channels), each rated as NOV's 2,200 hp 14-P-220: 1,980 hydraulic hp at 105
# strokes/min, up to 1,215 gal/min on 9-inch liners and 7,475 psi on 5.5-inch
# liners (NOV technical marketing sheet; data/materials.md).
RIG_PUMPS = 3
PUMP_HYD_HP = 1980.0
PUMP_MAX_GPM = 1215.0
PUMP_MAX_PSI = 7475.0


def rig_capacity_gpm(spp_Pa, pumps=RIG_PUMPS):
    """Flow [gal/min] the rig's pumps deliver at a standpipe pressure, choosing
    liners to suit; 0 above the pumps' pressure rating."""
    psi = spp_Pa / 6894.757
    if psi > PUMP_MAX_PSI:
        return 0.0
    return pumps * min(PUMP_HYD_HP * 1714.0 / max(psi, 1.0), PUMP_MAX_GPM)
GPM = 6.30902e-5   # m^3/s
Q_FACE_SITE = 30000.0
FILM_MULT_FORGE = 0.75   # fitted on FORGE 16B (validate_forge16b.py); a sensitivity


def site_well(z, pipe):
    shoes = tuple(min(s, f * z) for s, f in zip(wg.WU_SHOES, (0.05, 0.25, 0.5)))
    return wg.wu2025_well(z, pipe, shoes=shoes)


def circulate(site, z, geotherm, pipe, m_dot, exposure=None, film_mult=1.0):
    """Model 1 on the site well, with friction heating, the loop held at the
    model's operating pressure (5 MPa at the wellhead) so the return stays liquid. Solved first with the
    rock exposed for a year (the earlier setting, kept as `old`), then, if an
    exposure function is given, with that, seeded from the first."""
    g = site_well(z, pipe)
    kw = dict(m_dot=m_dot, T_inj=site.T_inj, k_rock=site.k_rock, geotherm=geotherm,
              geometry=g, target_rock_T=float(geotherm(z)), Q_face=Q_FACE_SITE,
              friction_heat=True, film_mult=film_mult, P_surface=m1.P_OP, verbose=False)
    old = m1.solve(**kw, t_years=1.0)
    if exposure is None:
        r = old
    else:
        r = m1.solve(**kw, exposure=exposure, init=old)
        if not r["success"]:
            r = m1.solve_seeded(exposure, **kw)
    r["old"] = old
    r["hyd"] = m1.hydraulics(m_dot, g, result=r, P_surface=m1.P_OP)
    return r


def ecd_sg(run, z):
    """Annular friction loss of a circulating run as an equivalent density [SG]."""
    return run["hyd"]["dp_annulus"] / (C.MUD_SG_GRAD * z)


def flow_gpm(run):
    """Volumetric flow [gal/min] at the inlet."""
    rho_in = float(m1.fluid_props(np.array([run["Td"][0]]), np.array([m1.P_OP]))[0][0])
    return run["m_dot"] / rho_in / GPM


def flow_for_go(site, z, geotherm, pipe, exposure, start, thermal=True):
    """The lowest flow, from `start` up, at which the verdict is GO within the
    pump limit and with the bit below its survival limit; None if there is
    none up to SITE_FLOWS[-1]."""
    T_rock = float(geotherm(z))
    for md in (f for f in SITE_FLOWS if f >= start):
        r = circulate(site, z, geotherm, pipe, md, exposure)
        if r["hyd"]["over_pump_limit"]:
            return None
        if r["T_bottom_delivered"] >= C.BHA_SURVIVAL_TEMP:
            continue
        ref, cases = stability_inputs(site, z, T_rock, T_wall=float(r["T_bottom_delivered"]),
                                      thermal=thermal, ecd_SG=ecd_sg(r, z))
        if (ref["window"]["verdict"] == "GO"
                and all(c["admissible"]["admissible"] for c in cases)):
            gpm = flow_gpm(r)
            cap = rig_capacity_gpm(float(r["hyd"]["spp"]))
            return dict(m_dot=md, bhct=float(r["T_bottom_delivered"]),
                        spp=float(r["hyd"]["spp"]), ecd_SG=ecd_sg(r, z), gpm=gpm,
                        above_forge=gpm > FORGE_MAX_FLOW_GPM, rig_capacity_gpm=cap,
                        within_rig=gpm <= cap)
    return None


def choose_circulation(site, z, geotherm, exposure=None, pipes=SITE_PIPES, film_mult=1.0):
    """The flow for the wall temperature: the larger of the lowest flow that
    keeps the bit below its survival limit and the lowest that meets the
    hole-cleaning annular velocity, with the standpipe pressure within the pump
    limit; and the best-insulated pipe for which such a flow exists.

    If no pipe meets all three, the conflict is reported and the verdict uses
    the lowest survivable flow within the pump limit, as before; if nothing
    survives, the best pipe's highest flow within the limit."""
    tried = []
    for pipe in pipes:
        t = dict(pipe=pipe, survival_flow=None, cleaning_flow=None, runs={})
        for md in SITE_FLOWS:
            r = circulate(site, z, geotherm, pipe, md, exposure, film_mult)
            if r["hyd"]["over_pump_limit"]:
                t["pump_limited_at"] = md
                break
            t["runs"][md] = r
            if t["survival_flow"] is None and r["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP:
                t["survival_flow"] = md
            if t["cleaning_flow"] is None and r["hyd"]["cleans_hole"]:
                t["cleaning_flow"] = md
            if t["survival_flow"] is not None and t["cleaning_flow"] is not None:
                return dict(pipe=pipe, m_dot=md, run=r, conflict=None, tried=tried + [t])
        tried.append(t)
    for t in tried:
        if t["survival_flow"] is not None:
            md = t["survival_flow"]
            why = ("no flow within the pump limit both keeps the bit below "
                   f"{C.BHA_SURVIVAL_TEMP:.0f} C and cleans the hole")
            return dict(pipe=t["pipe"], m_dot=md, run=t["runs"][md], conflict=why, tried=tried)
    t = next((t for t in tried if t["runs"]), tried[0])
    md = max(t["runs"]) if t["runs"] else SITE_FLOWS[0]
    r = t["runs"].get(md) or circulate(site, z, geotherm, t["pipe"], md, exposure, film_mult)
    return dict(pipe=t["pipe"], m_dot=md, run=r, tried=tried,
                conflict=f"no pipe keeps the bit below {C.BHA_SURVIVAL_TEMP:.0f} C "
                         "within the pump limit")


def evaluate(site: SiteProfile, v10=False, thermal=True, circulation="v1.3"):
    """Run Models 1-5 for one site.

    circulation="v1.3" adds the verdict through the drilling cycle (site_cycle)
    to circulation="v1.2", the circulation model as validated: the site well and
    string, a named pipe type, friction heating, the wall's exposure from Model
    3's rate of penetration, and the flow rule above. circulation="v1.1" keeps
    the earlier single-pipe model (a 0.02 W/m K wall, a year's exposure, the
    lowest survivable flow) for the comparison with v1.1.1; the v1.0 adapter
    always uses it.

    v10=True is the v1.0 adapter: it reproduces v1.0 exactly (the global C.UCS
    instead of the site's own), for the regression tests. Each v1.1 change that
    alters the stability numbers is switched off by it. thermal=False switches
    off the wall's thermal hoop stress alone. Either way the result carries
    three sensitivities: an uncooled wall (at rock temperature, so neither
    thermal stress nor cold-strength gain), the site-calibrated breakout
    limit, and the two together.
    """
    z = site.target(v10)
    geotherm = site.temperature(v10)
    T_rock = float(geotherm(z))

    if circulation == "v1.3" and not v10:
        return _evaluate_v13(site, thermal)
    if circulation == "v1.2" and not v10:
        return _evaluate_v12(site, thermal)

    # --- (1) SURVIVAL & ENERGY: Model 1 with the real LAYERED geotherm ---
    # scan flow for the minimum that keeps the bit < survival ceiling
    surv = None
    for md in (2, 4, 7, 10, 14, 20):
        r = m1.solve(m_dot=md, T_inj=site.T_inj, pipe=wg.LEGACY_VACUUM, k_rock=site.k_rock,
                     geotherm=geotherm, target_depth=z, Q_face=30000.0,
                     verbose=False)
        if r["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP:
            surv = (md, r); break
    if surv is None:
        surv = (20, r)
    m_min, m1_run = surv

    # wall temperature: the circulating bottom-hole temperature at the
    # minimum survivable flow, for strength and thermal stress alike; v1.0's
    # 200 C cap under the adapter
    T_wall = None if v10 else float(m1_run["T_bottom_delivered"])
    out = _stability(site, z, T_rock, T_wall, v10, thermal)
    out["m_min"], out["m1"] = m_min, m1_run
    K0, G_deep = out["K0"], out["G_deep"]

    # --- (2/3) DRILLABILITY: regime + quench ROP gain (spallation uses K0=Shmin)
    drill = m3.evaluate(G_deep, K0, out["m_min"], pipe=wg.LEGACY_VACUUM,
                        target_rock_T=T_rock, quench=True)
    drill_nq = m3.evaluate(G_deep, K0, out["m_min"], pipe=wg.LEGACY_VACUUM,
                           target_rock_T=T_rock, quench=False)
    out["drill"] = drill
    out["rop_gain"] = drill["ROP"] / drill_nq["ROP"] if drill_nq["ROP"] > 0 else np.nan

    # production-flow energy (the survival run uses min flow, which minimises MW)
    rp = m1.solve(m_dot=10.0, T_inj=site.T_inj, pipe=wg.LEGACY_VACUUM, k_rock=site.k_rock,
                  geotherm=geotherm, target_depth=z, Q_face=30000.0, verbose=False)
    out["MW_prod"] = rp["Q_product"] / 1e6
    out["Tret_prod"] = rp["T_return_surface"]
    out["circulation"] = "v1.1"
    return out


def _stability(site, z, T_rock, T_wall, v10=False, thermal=True, ecd_SG=0.0):
    """Stability (Models 4 and 5) for a wall temperature, with the three
    sensitivities; ecd_SG is the circulating friction on the fracture bound."""
    ref, cases = stability_inputs(site, z, T_rock, v10=v10, T_wall=T_wall,
                                  thermal=thermal, ecd_SG=ecd_SG)
    # sensitivities, each a full re-run of the stability check (no Model 1).
    # "No wall cooling" puts the wall at rock temperature: no thermal hoop
    # stress and no strength gain from cooling, i.e. an uncooled hole.
    w63 = f"W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f}"
    sens = {} if v10 else {
        "no wall cooling": stability_inputs(site, z, T_rock, T_wall=T_rock,
                                            thermal=False, ecd_SG=ecd_SG),
        w63: stability_inputs(site, z, T_rock, T_wall=T_wall, thermal=thermal,
                              W_max=C.BREAKOUT_W_MAX_SITE_DEG, ecd_SG=ecd_SG),
        f"no wall cooling, {w63}": stability_inputs(
            site, z, T_rock, T_wall=T_rock, thermal=False,
            W_max=C.BREAKOUT_W_MAX_SITE_DEG, ecd_SG=ecd_SG)}
    UCS = ref["UCS"]
    Sv, Shmin, SHmax = ref["Sv"], ref["Shmin"], ref["SHmax"]
    P_fluid = site.rho_fluid_grad * z
    # Shmin/Sv at the target from the reference case; v1.0 used the site ratio
    K0 = site.K0_min if v10 else Shmin / Sv
    # effective gradient that reproduces the real target depth in scalar-G models
    G_deep = (T_rock - C.SURFACE_TEMP) / z
    out = dict(site=site, z=z, T_rock=T_rock, Sv=Sv, Shmin=Shmin, SHmax=SHmax,
               P_fluid=P_fluid, UCS=UCS, v10=v10,
               K0=K0, anisotropy=SHmax / Shmin, Pp=ref["Pp"],
               stress=ref, stress_cases=cases, thermal=thermal and not v10,
               T_wall=ref["T_wall"], dsigma_T=ref["dsigma_T"],
               sensitivity={k: dict(stress=r, stress_cases=c) for k, (r, c) in sens.items()},
               no_thermal=None if v10 else dict(stress=sens["no wall cooling"][0],
                                                stress_cases=sens["no wall cooling"][1]),
               G_deep=G_deep)

    # --- (4) CREEP closure: hot vs cooled wall ---
    # cooled wall at the same temperature as the stability checks; 200 C
    # under the v1.0 adapter
    out["creep_hot"] = m4.closure_rate(z, T_rock, ref["T_wall"], 7*86400, cooled=False) * YEAR * 100
    out["creep_cold"] = m4.closure_rate(z, T_rock, ref["T_wall"], 7*86400, cooled=True) * YEAR * 100

    # --- (5) STABILITY: isotropic GRC (uses Shmin as p0) + anisotropic breakout ---
    u_hot, rp_hot, pcr_h, reg_h = m5.grc(P_fluid, Shmin, T_rock, UCS=UCS)
    u_cold, rp_cold, pcr_c, reg_c = m5.grc(P_fluid, Shmin, ref["T_wall"], UCS=UCS)
    out["grc"] = dict(u_hot=u_hot, u_cold=u_cold, rp_hot=rp_hot, rp_cold=rp_cold,
                      reg_h=reg_h, reg_c=reg_c)
    out["breakout"] = ref["breakout"]
    return out


def _evaluate_v12(site, thermal=True):
    z = site.target()
    geotherm = site.temperature()
    T_rock = float(geotherm(z))
    # Shmin/Sv for Model 3 doesn't depend on the wall temperature
    pre = stability_inputs(site, z, T_rock, T_wall=T_rock, thermal=False)[0]
    K0 = pre["Shmin"] / pre["Sv"]
    G_deep = (T_rock - C.SURFACE_TEMP) / z

    # pass 1, rock exposed for a year: pipe and flow, then Model 3's rate of
    # penetration for them; pass 2, the wall's exposure from that rate
    first = choose_circulation(site, z, geotherm)
    rop = m3.evaluate(G_deep, K0, first["m_dot"], pipe=first["pipe"],
                      target_rock_T=T_rock, quench=True)["ROP"]
    exposure = m1.exposure_from_rop(z, rop * 3600.0, 0.0)
    circ = choose_circulation(site, z, geotherm, exposure)
    run = circ["run"]
    T_wall = float(run["T_bottom_delivered"])
    ecd = ecd_sg(run, z)

    out = _stability(site, z, T_rock, T_wall, thermal=thermal, ecd_SG=ecd)
    out["m_min"], out["m1"] = circ["m_dot"], run
    out["circulation"] = "v1.2"
    out["circ"] = dict(pipe=circ["pipe"], m_dot=circ["m_dot"], conflict=circ["conflict"],
                       rop=rop, hyd=run["hyd"], bhct=T_wall, ecd_SG=ecd,
                       gpm=flow_gpm(run), above_forge=flow_gpm(run) > FORGE_MAX_FLOW_GPM,
                       bhct_one_year=float(run["old"]["T_bottom_delivered"]),
                       tried=[dict(pipe=t["pipe"].name, survival_flow=t["survival_flow"],
                                   cleaning_flow=t["cleaning_flow"],
                                   pump_limited_at=t.get("pump_limited_at"))
                              for t in circ["tried"]])

    drill = m3.evaluate(G_deep, K0, circ["m_dot"], pipe=circ["pipe"],
                        target_rock_T=T_rock, quench=True)
    drill_nq = m3.evaluate(G_deep, K0, circ["m_dot"], pipe=circ["pipe"],
                           target_rock_T=T_rock, quench=False)
    out["drill"] = drill
    out["rop_gain"] = drill["ROP"] / drill_nq["ROP"] if drill_nq["ROP"] > 0 else np.nan
    # heat returned while drilling, early life, single loop
    out["MW_prod"] = run["Q_product"] / 1e6
    out["Tret_prod"] = run["T_return_surface"]

    # the lowest flow at which the same pipe gives GO, within the pump limit
    out["circ"]["go"] = flow_for_go(site, z, geotherm, circ["pipe"], exposure,
                                    circ["m_dot"], thermal)

    # circulation sensitivities, each with its own wall temperature and
    # circulating friction
    film = circulate(site, z, geotherm, circ["pipe"], circ["m_dot"], exposure, FILM_MULT_FORGE)
    alt = {
        "rock exposed 1 yr": (run["old"]["T_bottom_delivered"], ecd),
        f"film x{FILM_MULT_FORGE}": (film["T_bottom_delivered"], ecd_sg(film, z)),
        "static fracture bound": (T_wall, 0.0),
    }
    dual = choose_circulation(site, z, geotherm, exposure, pipes=(wg.DUAL_WALL,))
    out["circ"]["dual_wall"] = dict(pipe=dual["pipe"].name, m_dot=dual["m_dot"],
                                    conflict=dual["conflict"],
                                    bhct=float(dual["run"]["T_bottom_delivered"]),
                                    spp=float(dual["run"]["hyd"]["spp"]))
    alt["dual-wall pipe (not manufactured)"] = (dual["run"]["T_bottom_delivered"],
                                                ecd_sg(dual["run"], z))
    vit = choose_circulation(site, z, geotherm, exposure, pipes=(wg.VACUUM_INSULATED,))
    out["circ"]["vacuum"] = dict(m_dot=vit["m_dot"], conflict=vit["conflict"],
                                 bhct=float(vit["run"]["T_bottom_delivered"]))
    alt["vacuum-insulated pipe"] = (vit["run"]["T_bottom_delivered"], ecd_sg(vit["run"], z))
    out["circ_sensitivity"] = {}
    for k, (Tw, e) in alt.items():
        r, c = stability_inputs(site, z, T_rock, T_wall=float(Tw), thermal=thermal, ecd_SG=e)
        out["circ_sensitivity"][k] = dict(stress=r, stress_cases=c, T_wall=float(Tw), ecd_SG=e)
    return out


# ------------------------------------------------------- the drilling cycle
SAFE_PAUSE_MAX_H = 168.0     # a week; a longer safe pause is reported as "> 7 d"


def _bottom_conditions(run, z):
    """Annulus temperature [C], film coefficient, hole radius and string OD at
    measured depth z, and the bottoms-up time [h], from a circulating run."""
    g = run["geometry"]
    zz = min(z, g.md_bit - 0.5)
    ga = g.at(np.array([zz]))
    Tu = float(np.interp(zz, run["z"], run["Tu"]))
    h = float(m1.h_dittus(run["m_dot"], ga["A_ann"], ga["Dh_ann"], np.array([Tu]),
                          m1.P_of_z(ga["tvd"]))[0])
    top, bot = m1._pieces(g)
    A = g.at(0.5 * (top + bot))["A_ann"]
    rho = float(m1.fluid_props(np.array([run["Td"][0]]), np.array([m1.P_OP]))[0][0])
    bu_h = float(np.sum(A * (bot - top))) * rho / run["m_dot"] / 3600.0
    return dict(Tu=Tu, h=h, r_wall=float(ga["r_wall"][0]), r_string=float(ga["r_out"][0]),
                bu_h=bu_h)


def _static_node(T_rock, z_tvd, r_wall, r_string):
    """Heat capacity [J/m/K] and laminar film of the water in a static hole,
    around the string (r_string > 0) or filling it."""
    P, T = m1.P_of_z(np.array([z_tvd])), np.array([T_rock])
    rho, cp, mu, k = (float(x[0]) for x in m1.fluid_props(T, P))
    area = np.pi * (r_wall ** 2 - r_string ** 2)
    Dh = 2 * (r_wall - r_string)
    return rho * cp * area, wall_thermal.NU_LAMINAR * k / Dh


def _ref_depth(site, z, T_rock, Pw):
    """The deepest reach of the uncooled breakout across the site's stress and
    strength cases [m beyond the wall]: the depth whose temperature the wall
    check uses (model5.breakout_with_skin)."""
    thermo = -thermal_hoop_stress(T_rock, T_rock - 1.0)
    d = 0.0
    for p in site.stress_cases:
        for _, ucs in site.strength_cases:
            d = max(d, m5.uncooled_depth(p.SHmax(z), p.Shmin(z), Pw, float(p.Pp(z)), 0.11,
                                         T_rock, ucs, thermo))
    return d


def site_cycle(site, out, connection="p90", trip_scale=1.0, run=None, label="base", node=True,
               form_h=0.0, fresh="p10"):
    """Verdicts through the drilling cycle for the wall element one stand above
    the end of a bit run at the target: drilling (circulating pressure, the
    wall as the bit passes it), the connection (static, the wall at the end of
    the 90th percentile connection) and the trip (static, the wall at the end
    of the trip), with the safe pause and the trip fluid weight.

    Drilling is judged when the bit exposes the rock. Elastic failure is
    immediate, and the rock at the depth the breakout reaches has not yet been
    cooled then, so the fluid must hold the wall at formation temperature.
    form_h: a time [h] the breakout is assumed to take to form, during which
    circulation cools that rock (the sensitivity; no data support one).

    The connection is also checked on the rock the bit has just drilled: it
    is circulated against only from the last new hole until the pumps stop
    (fresh: FORGE's 10th percentile, "p10", or its median), then sits through
    the connection at static pressure. The drilling fluid must hold the wall
    in every state and stay below the fracture limit while circulating."""
    import drilling_cycle as dc
    t = dc.forge_timings()
    z = out["z"]
    T_rock = out["T_rock"]
    run = run or out["m1"]
    bc = _bottom_conditions(run, z)
    a = bc["r_wall"]
    rop = out["drill"]["ROP"] * 3600.0
    hist = dc.element_history(z, rop, bc["bu_h"], t, connection, trip_scale)
    Pw = site.rho_fluid_grad * z
    d_ref = _ref_depth(site, z, T_rock, Pw)
    states = []
    for st in hist:
        if st.kind == "circulating":
            states.append(wall_thermal.State("circulating", st.hours * 3600, T_fluid=bc["Tu"], h=bc["h"]))
        else:
            r_str = bc["r_string"] if st.label == "connection" else 0.0
            cap, hs = _static_node(T_rock, z, a, r_str)
            states.append(wall_thermal.State("static", st.hours * 3600, h=hs, C=cap if node else 0.0))
    total = sum(st.hours for st in hist) * 3600 + SAFE_PAUSE_MAX_H * 3600
    r_max = wall_thermal.r_max_for(a, site.k_rock / (C.RHO_ROCK * C.CP_ROCK), total)
    res = wall_thermal.run(states, T_rock, a, site.k_rock, C.RHO_ROCK * C.CP_ROCK, n=150,
                           r_max=r_max, growth=1.2)
    T_ref = lambda field: float(np.interp(a + d_ref, res.r, field))
    T_conn, T_trip, T_end = T_ref(res.fields[1]), T_ref(res.fields[4]), T_ref(res.fields[3])
    # the connection on freshly drilled rock at the bottom of the hole
    fresh_min = t.fresh_circ_p10_min if fresh == "p10" else t.fresh_circ_median_min
    conn_h = (t.connection_median_min if connection == "median" else t.connection_p90_min) / 60.0
    cap_f, hs_f = _static_node(T_rock, z, a, bc["r_string"])
    st_f = [wall_thermal.State("circulating", fresh_min * 60, T_fluid=bc["Tu"], h=bc["h"]),
            wall_thermal.State("static", conn_h * 3600, h=hs_f, C=cap_f if node else 0.0,
                               T_node0=bc["Tu"])]
    # with time allowed for breakouts to form, the fresh rock is judged when
    # they would form, circulated against again after the connection
    rest = form_h * 3600 - (fresh_min * 60 + conn_h * 3600)
    if rest > 0:
        st_f.append(wall_thermal.State("circulating", rest, T_fluid=bc["Tu"], h=bc["h"]))
    r_f = wall_thermal.run(st_f, T_rock, a, site.k_rock, C.RHO_ROCK * C.CP_ROCK, n=150, r_max=r_max,
                           growth=1.2)
    T_fresh = T_ref(r_f.fields[-1])
    if form_h > 0:
        r0 = wall_thermal.run([wall_thermal.State("circulating", form_h * 3600, T_fluid=bc["Tu"], h=bc["h"])],
                              T_rock, a, site.k_rock, C.RHO_ROCK * C.CP_ROCK, n=150, r_max=r_max,
                              growth=1.2)
        T_drill = T_ref(r0.fields[-1])
    else:
        T_drill = T_rock
    ecd = ecd_sg(run, z)
    sg = C.MUD_SG_GRAD * z
    ref_d, cases_d = stability_inputs(site, z, T_rock, T_wall=T_drill, Pw_add=ecd * sg)
    ref_c, cases_c = stability_inputs(site, z, T_rock, T_wall=T_conn)
    ref_t, cases_t = stability_inputs(site, z, T_rock, T_wall=T_trip)
    ref_f, cases_f = stability_inputs(site, z, T_rock, T_wall=T_fresh)
    states_out = {"drilling": (ref_d, T_drill), "connection": (ref_c, T_conn),
                  "connection, fresh rock": (ref_f, T_fresh), "trip": (ref_t, T_trip)}

    # Each stress and strength case is followed through the cycle on its own:
    # its drilling fluid weight (the lightest static weight that holds the wall
    # while circulating, less the friction it was checked with, and through a
    # connection, never below water), its trip fluid weight and its verdict.
    # The site takes the worst case, with the range across cases; the weights
    # and safe pause reported are the deciding case's (the worst verdict, then
    # the heaviest trip fluid).
    water = site.rho_fluid_grad / C.MUD_SG_GRAD
    need = lambda w: w["SG_lo"] if w["verdict"] != "GO" else 0.0
    per_case = []
    for cd, cc, cf, ct in zip(cases_d, cases_c, cases_f, cases_t):
        drill = max(water, need(cd["window"]) - ecd, need(cc["window"]), need(cf["window"]))
        # the drilling fluid, set by whichever state needs most, must not
        # fracture the rock while circulating
        drill_hi = cd["window"]["SG_hi"] - ecd
        v_d = cd["window"]["verdict"] if drill <= drill_hi + 1e-9 else "NO-GO"
        vs = [v_d, cc["window"]["verdict"], cf["window"]["verdict"], ct["window"]["verdict"]]
        v = "GO" if all(x == "GO" for x in vs) else ("CONDITIONAL" if all(x != "NO-GO" for x in vs)
                                                       else "NO-GO")
        per_case.append(dict(label=ct["label"], profile=ct["profile"], UCS=ct["UCS"], verdict=v,
                             states=dict(zip(("drilling", "connection", "connection, fresh rock",
                                              "trip"), vs)),
                             drill_SG=drill, drill_SG_hi=drill_hi,
                             trip_SG=max(water, need(ct["window"])), trip_SG_hi=ct["window"]["SG_hi"]))
    gov = max(per_case, key=lambda c: (VERDICT_ORDER.index(c["verdict"]), c["trip_SG"]))
    w_drill_SG = gov["drill_SG"]
    Pw_drill = w_drill_SG * sg
    thermo = -thermal_hoop_stress(T_rock, T_rock - 1.0)
    p_g, ucs_g = gov["profile"], gov["UCS"]

    def width_ok(T_r, Pw_drill=Pw_drill):
        w = m5.breakout_width(p_g.SHmax(z), p_g.Shmin(z), Pw_drill, float(p_g.Pp(z)), T_r, ucs_g,
                              thermo * (T_r - T_rock))
        return w <= C.BREAKOUT_W_MAX_DEG + 1e-6

    field_drill = res.fields[3]
    C_s, h_s = _static_node(T_rock, z, a, bc["r_string"])
    C_s = C_s if node else 0.0

    def after_pause(hours):
        rr = wall_thermal.run([wall_thermal.State("static", hours * 3600, h=h_s, C=C_s, T_node0=bc["Tu"])],
                              T_rock, a, site.k_rock, C.RHO_ROCK * C.CP_ROCK, n=150, r_max=r_max,
                              T0=field_drill, growth=1.2)
        return T_ref(rr.fields[-1])

    def safe_pause_at(SG):
        """The longest pause [h] from the end of drilling for which the
        deciding case's width stays within W_max at a static fluid weight SG
        (0 if it fails at once, inf beyond a week)."""
        P = SG * sg
        if not width_ok(T_ref(field_drill), P):
            return 0.0
        if width_ok(after_pause(SAFE_PAUSE_MAX_H), P):
            return float("inf")
        lo, hi = 0.0, SAFE_PAUSE_MAX_H
        for _ in range(18):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if width_ok(after_pause(mid), P) else (lo, mid)
        return 0.5 * (lo + hi)

    # The heuristic (the wall check at the temperature of the rock at the
    # uncooled breakout's depth) against the full stress field in the first
    # minutes after exposure, when the contracting skin pushes hoop stress into
    # the rock behind it: the failed depth at the Shmin azimuth, at the
    # drilling fluid with the circulating friction, for the deciding case.
    Pw_circ = Pw_drill + ecd * sg
    skin = []
    for minutes in (0.0, 1.0, 2.0, 5.0, 10.0):
        if minutes == 0.0:
            field = np.full_like(res.r, T_rock)
        else:
            rr = wall_thermal.run([wall_thermal.State("circulating", minutes * 60, T_fluid=bc["Tu"],
                                                      h=bc["h"])], T_rock, a, site.k_rock,
                                  C.RHO_ROCK * C.CP_ROCK, n=150, r_max=r_max, growth=1.2)
            field = rr.fields[-1]
        b = m5.breakout_behind_wall(p_g.SHmax(z), p_g.Shmin(z), Pw_circ, float(p_g.Pp(z)), a, res.r,
                                    field, T_rock, ucs_g, thermo, d_theta=1.0)
        skin.append(dict(minutes=minutes, depth=b["depth"], width_wall=b["width_wall"],
                         T_ref=T_ref(field)))
    depth_unc = m5.uncooled_depth(p_g.SHmax(z), p_g.Shmin(z), Pw_circ, float(p_g.Pp(z)), a, T_rock,
                                  ucs_g, thermo)

    safe = safe_pause_at(w_drill_SG)
    trip_h = hist[-1].hours
    order = [c["verdict"] for c in per_case]
    site_v = worst_verdict(order)
    if len(set(order)) > 1:
        site_v += f" [{min(order, key=VERDICT_ORDER.index)}..{worst_verdict(order)}]"
    verdicts = {k: worst_verdict([c["states"][k] for c in per_case])
                for k in ("drilling", "connection", "connection, fresh rock", "trip")}
    return dict(label=label, z=z, T_rock=T_rock, rop_m_h=rop, bottoms_up_h=bc["bu_h"],
                T_fluid=bc["Tu"], h=bc["h"], depth_ref=d_ref, ecd_SG=ecd,
                T_ref={**{k: v[1] for k, v in states_out.items()}, "end of drilling": T_end},
                states={k: v[0] for k, v in states_out.items()},
                cases={"drilling": cases_d, "connection": cases_c, "trip": cases_t},
                per_case=per_case, deciding_case=gov["label"],
                verdicts=verdicts,
                drill_SG=w_drill_SG, safe_pause_h=safe, trip_h=trip_h,
                trip_SG=gov["trip_SG"], trip_SG_hi=gov["trip_SG_hi"], drill_SG_hi=gov["drill_SG_hi"],
                fresh_circ_min=fresh_min,
                trip_window_open=bool(gov["trip_SG"] <= gov["trip_SG_hi"]),
                site_verdict=site_v, history=hist, skin=skin, skin_depth_uncooled=depth_unc, result=res, r_ref=a + d_ref, a=a,
                safe_pause_at=safe_pause_at, after_pause=after_pause, width_ok=width_ok)


def trip_profile(site, out, heights=(0.0, 100.0, 250.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0,
                                      5000.0, 6000.0)):
    """The trip verdict up the open hole: wall elements at heights above the
    end of the bit run, up to the casing shoe, each circulated since the bit
    passed at the rate of penetration (connections ignored) and through the
    bottoms-up, then the trip."""
    import drilling_cycle as dc
    t = dc.forge_timings()
    run = out["m1"]
    z_bit = out["z"]
    geotherm = site.temperature()
    shoe = max(h.md_top for h in run["geometry"].hole)
    rop = out["drill"]["ROP"] * 3600.0
    rows = []
    for H in heights:
        z = z_bit - t.stand_ft * dc.FT - H
        if z <= shoe:
            break
        T_rock = float(geotherm(z))
        bc = _bottom_conditions(run, z)
        trip = dc.trip_hours(z_bit, t)
        cap, hs = _static_node(T_rock, z, bc["r_wall"], 0.0)
        states = [wall_thermal.State("circulating", ((H + t.stand_ft * dc.FT) / rop + bc["bu_h"]) * 3600,
                                     T_fluid=bc["Tu"], h=bc["h"]),
                  wall_thermal.State("static", trip * 3600, h=hs, C=cap)]
        res = wall_thermal.run(states, T_rock, bc["r_wall"], site.k_rock, C.RHO_ROCK * C.CP_ROCK,
                               n=120, growth=1.2)
        d_ref = _ref_depth(site, z, T_rock, site.rho_fluid_grad * z)
        T_ref = float(np.interp(bc["r_wall"] + d_ref, res.r, res.fields[-1]))
        ref, _ = stability_inputs(site, z, T_rock, T_wall=T_ref)
        w = ref["window"]
        rows.append(dict(height=H, z=z, T_rock=T_rock, T_ref=T_ref, verdict=w["verdict"],
                         SG_lo=w["SG_lo"], SG_hi=w["SG_hi"],
                         SG_needed=w["SG_hydro"] if w["verdict"] == "GO" else w["SG_lo"]))
    return rows


def _evaluate_v13(site, thermal=True):
    """The v1.2 evaluation, with the verdict through the drilling cycle."""
    out = _evaluate_v12(site, thermal)
    out["circulation"] = "v1.3"
    out["cycle"] = site_cycle(site, out)
    z, geotherm = out["z"], site.temperature()
    rop = out["circ"]["rop"] * 3600.0
    exposure = m1.exposure_from_rop(z, rop, 0.0)
    sens = {
        "half trip time": site_cycle(site, out, trip_scale=0.5, label="half trip time"),
        "double trip time": site_cycle(site, out, trip_scale=2.0, label="double trip time"),
        "no fluid node": site_cycle(site, out, label="no fluid node", node=False),
        "breakouts form over 1 h": site_cycle(site, out, label="breakouts form over 1 h", form_h=1.0),
        "median circulation before a connection": site_cycle(
            site, out, label="median circulation before a connection", fresh="median"),
    }
    dual = out["circ"]["dual_wall"]
    run_dual = circulate(site, z, geotherm, wg.DUAL_WALL, dual["m_dot"], exposure)
    sens["dual-wall pipe"] = site_cycle(site, out, run=run_dual, label="dual-wall pipe")
    # the highest flow within the pump limit on the base pipe
    best = None
    for md in (f for f in SITE_FLOWS if f > out["circ"]["m_dot"]):
        r = circulate(site, z, geotherm, out["circ"]["pipe"], md, exposure)
        if r["hyd"]["over_pump_limit"]:
            break
        best = r
    if best is not None:
        sens[f"highest flow, {best['m_dot']:.0f} kg/s"] = site_cycle(
            site, out, run=best, label=f"highest flow, {best['m_dot']:.0f} kg/s")
    out["cycle_sensitivity"] = sens
    out["trip_profile"] = trip_profile(site, out)
    return out


def report(o):
    s = o["site"]
    print("=" * 80)
    print(f"SITE EVALUATION:  {s.name}")
    print("=" * 80)
    print(f"  Target: {o['T_rock']:.0f} C rock at {o['z']/1000:.1f} km "
          f"(layered geotherm; supercritical-class)")
    print(f"  In-situ stress: Sv={o['Sv']/1e6:.0f}  SHmax={o['SHmax']/1e6:.0f}  "
          f"Shmin={o['Shmin']/1e6:.0f} MPa  | K0={o['K0']:.2f}  anisotropy={o['anisotropy']:.2f}")
    print(f"  Stress basis: {s.stress_basis} ({o['stress']['profile'].source}); "
          f"pore pressure {o['Pp']/1e6:.0f} MPa")
    a = o["stress"]["admissible"]
    print(f"  Admissibility at mu {a['mu']}: effective S1/S3 {a['ratio']:.2f} vs cap "
          f"{a['cap']:.2f} -> {'admissible' if a['admissible'] else 'INADMISSIBLE'}")
    print(f"  Fluid pressure (hydrostatic): {o['P_fluid']/1e6:.0f} MPa")
    print("-" * 80)
    md, r = o["m_min"], o["m1"]
    print(f"  [1] TOOL SURVIVAL & ENERGY (layered geotherm, vacuum tubing)")
    print(f"      min flow to survive: {md} kg/s -> bit {r['T_bottom_delivered']:.0f} C "
          f"(ceiling {C.BHA_SURVIVAL_TEMP:.0f}) {'OK' if r['T_bottom_delivered']<200 else 'FAIL'}")
    print(f"      surface return {r['T_return_surface']:.0f} C | early-life {r['Q_product']/1e6:.1f} MW_th")
    print(f"      at production flow 10 kg/s: {o['MW_prod']:.1f} MW_th, return {o['Tret_prod']:.0f} C")
    print("-" * 80)
    d = o["drill"]
    print(f"  [2/3] DRILLABILITY (quench-assist; spallation uses Shmin/Sv={o['K0']:.2f})")
    print(f"      regime: {d['regime']} | effective MSE {d['MSE_eff']/1e6:.0f} MPa | "
          f"ROP {d['ROP']*3600:.1f} m/hr | quench gain {o['rop_gain']:.1f}x")
    print("-" * 80)
    print(f"  [4] CREEP CLOSURE (time-dependent)")
    print(f"      hot wall: {o['creep_hot']:.2e} %/yr  ->  cooled {o['T_wall']:.0f}C: {o['creep_cold']:.2e} %/yr")
    print(f"      cooling factor ~ {o['creep_hot']/max(o['creep_cold'],1e-30):.1e}x  (squeezing controlled)")
    print("-" * 80)
    g = o["grc"]; b = o["breakout"]
    print(f"  [5] HOLE STABILITY (rock-mass UCS {o['UCS']/1e6:.0f} MPa at 25 C)")
    print(f"      isotropic GRC (p0=Shmin): {g['reg_c']}, convergence {g['u_cold']*1000:.1f} mm, "
          f"plastic r/a {g['rp_cold']:.2f}")
    print(f"      ANISOTROPIC breakout: sigma_theta={b['sigma_theta']/1e6:.0f} MPa vs "
          f"MC-limit cold {b['mc_cold']/1e6:.0f} MPa -> "
          f"{'BREAKS OUT' if b['sigma_theta']>b['mc_cold'] else 'stable'}")
    print(f"      mud weight to suppress breakout: {b['P_need']/1e6:.0f} MPa vs "
          f"hydrostatic {o['P_fluid']/1e6:.0f} MPa = +{b['overbalance_MPa']:.0f} MPa "
          f"({b['over_hydro']:.2f}x hydrostatic)")
    w = o["stress"]["window"]
    if w is not None:
        print(f"      wall {o['T_wall']:.0f} C from Model 1: thermal hoop stress "
              f"{o['dsigma_T']/1e6:+.0f} MPa{'' if o['thermal'] else ' (switched off)'}")
        print(f"      at hydrostatic mud: {describe_width(w['width_hydro'])}, "
              f"limit {w['W_max']:.0f} deg; mud window {w['SG_lo']:.2f}-{w['SG_hi']:.2f} SG "
              f"-> {w['verdict']}")
    print("=" * 80)
    print("  VERDICT")
    verdict_lines(o)
    print("=" * 80)


def verdict_lines(o):
    s = o["site"]; b = o["breakout"]
    print(f"   + Extensional graben (K0={o['K0']:.2f}): low confinement aids quench & limits")
    print(f"     the minimum-stress; tool survival closes with active cooling.")
    print(f"   + Cooling controls time-dependent creep by ~{o['creep_hot']/max(o['creep_cold'],1e-30):.0e}x.")
    print(f"   - HIGH anisotropy (SHmax/Shmin={o['anisotropy']:.2f}) drives breakout, but only marginally:")
    print(f"     +{b['overbalance_MPa']:.0f} MPa overbalance ({b['over_hydro']:.2f}x hydrostatic) suppresses it,")
    print(f"     plus cooling margin. The binding constraint -- as in the real GPK wells, which")
    print(f"     broke out yet were drilled to 5 km -- not a showstopper.")
    print(f"   ~ At {o['T_rock']:.0f} C the rock is at the brittle-ductile edge; a 374 C/~10 km target")
    print(f"     sits more safely in the brittle field.")
    print(f"   => CONDITIONAL GO: drillable & survivable with quench+cooling; stability is mud-")
    print(f"      weight-limited by anisotropy, not temperature. A lower-anisotropy URG segment")
    print(f"      would score higher -- which is exactly what the tool is for.")


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    o = evaluate(SOULTZ)
    report(o)
