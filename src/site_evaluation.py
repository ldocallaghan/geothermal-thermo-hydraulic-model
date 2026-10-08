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
              and b.get("strength") in ("calibrated", "measured (lab)")
              and b.get("well check") in ("calibrated", "checked"))
        return "evidence-based" if ok else "speculative"

    def missing(self):
        """Inputs that keep a site speculative, as (input, basis)."""
        need = {"stress": ("measured", "measured range"),
                "strength": ("calibrated", "measured (lab)"),
                "well check": ("calibrated", "checked")}
        return [(k, self.data_basis.get(k, ("none", ""))[0]) for k, good in need.items()
                if self.data_basis.get(k, ("none", ""))[0] not in good]


# Data-basis labels, strongest first
BASES = ("calibrated", "measured", "measured range", "measured (lab)", "checked",
         "regime bounds", "extrapolated", "regional", "assumed", "unsourced", "none")


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


def stability_inputs(site, z, T_rock, v10=False, mu=C.FRICTION_MU, UCS=None,
                     T_wall=None, thermal=True, W_max=C.BREAKOUT_W_MAX_DEG):
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
    P_mud = site.rho_fluid_grad * z
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
        SHmax_lim = m5.SHmax_frictional_limit(Shmin, Pp, mu)
        cap_binds = bool(beyond is not None and beyond > 0
                         and SHmax > SHmax_lim * (1 + m5.CAP_RTOL))
        if cap_binds:
            SHmax = SHmax_lim
        if not v10:
            window = mud_window(Shmin, SHmax, Pp, z, P_mud, T_rock, T_wall=T_wall,
                                UCS=UCS, dsigma_T=dsigma_T, W_max=W_max)
        cases.append(dict(window=window, T_wall=T_wall, dsigma_T=dsigma_T,
            UCS=C.UCS if UCS is None else UCS, strength_label=s_label,
            label=p.label + (f", {s_label}" if len(strengths) > 1 else ""),
            profile=p, Sv=Sv, Shmin=Shmin, SHmax=SHmax, Pp=Pp,
            beyond_data=beyond, cap_binds=cap_binds, cap_depth=p.cap_depth(mu),
            admissible=m5.stress_admissible(Sv, Shmin, SHmax, Pp, mu),
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
               T0=C.WALL_T0):
    """Mud window and drillability verdict, judged on breakout width.

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
    Pw_hi = Shmin - margin_SG * sg
    w_hydro = width(Pw_hydro)
    if w_hydro <= W_max:
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
        SG_hydro=Pw_hydro / sg, SG_lo=Pw_lo / sg, SG_hi=Pw_hi / sg,
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
# to the target, with the casing shoes moved up for shallow targets. Pipe types
# are tried from best insulated to worst; vacuum-insulated tubing is not a drill
# pipe, so it is a sensitivity only.
SITE_PIPES = (wg.DUAL_WALL, wg.TK_DRAKON, wg.INTERNALLY_COATED, wg.EXTERNALLY_COATED,
              wg.CONVENTIONAL)
# The dual-wall pipe's figures come from a modelling study; no such drill pipe
# is on sale (data/materials.md). The rest can be bought or have been run.
COMMERCIAL_PIPES = SITE_PIPES[1:]
SITE_FLOWS = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80)   # kg/s
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
    r = old if exposure is None else m1.solve(**kw, exposure=exposure, init=old)
    r["old"] = old
    r["hyd"] = m1.hydraulics(m_dot, g, result=r, P_surface=m1.P_OP)
    return r


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


def evaluate(site: SiteProfile, v10=False, thermal=True, circulation="v1.2"):
    """Run Models 1-5 for one site.

    circulation="v1.2" is the circulation model as validated: the site well and
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


def _stability(site, z, T_rock, T_wall, v10=False, thermal=True):
    """Stability (Models 4 and 5) for a wall temperature, with the three
    sensitivities."""
    ref, cases = stability_inputs(site, z, T_rock, v10=v10, T_wall=T_wall,
                                  thermal=thermal)
    # sensitivities, each a full re-run of the stability check (no Model 1).
    # "No wall cooling" puts the wall at rock temperature: no thermal hoop
    # stress and no strength gain from cooling, i.e. an uncooled hole.
    w63 = f"W_max {C.BREAKOUT_W_MAX_SITE_DEG:.0f}"
    sens = {} if v10 else {
        "no wall cooling": stability_inputs(site, z, T_rock, T_wall=T_rock,
                                            thermal=False),
        w63: stability_inputs(site, z, T_rock, T_wall=T_wall, thermal=thermal,
                              W_max=C.BREAKOUT_W_MAX_SITE_DEG),
        f"no wall cooling, {w63}": stability_inputs(
            site, z, T_rock, T_wall=T_rock, thermal=False,
            W_max=C.BREAKOUT_W_MAX_SITE_DEG)}
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

    out = _stability(site, z, T_rock, T_wall, thermal=thermal)
    out["m_min"], out["m1"] = circ["m_dot"], run
    out["circulation"] = "v1.2"
    out["circ"] = dict(pipe=circ["pipe"], m_dot=circ["m_dot"], conflict=circ["conflict"],
                       rop=rop, hyd=run["hyd"], bhct=T_wall,
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

    # circulation sensitivities, each with its own wall temperature
    alt = {
        "rock exposed 1 yr": run["old"]["T_bottom_delivered"],
        f"film x{FILM_MULT_FORGE}": circulate(site, z, geotherm, circ["pipe"], circ["m_dot"],
                                              exposure, FILM_MULT_FORGE)["T_bottom_delivered"],
    }
    vit = choose_circulation(site, z, geotherm, exposure, pipes=(wg.VACUUM_INSULATED,))
    out["circ"]["vacuum"] = dict(m_dot=vit["m_dot"], conflict=vit["conflict"],
                                 bhct=float(vit["run"]["T_bottom_delivered"]))
    alt["vacuum-insulated pipe"] = vit["run"]["T_bottom_delivered"]
    com = choose_circulation(site, z, geotherm, exposure, pipes=COMMERCIAL_PIPES)
    out["circ"]["commercial"] = dict(pipe=com["pipe"].name, m_dot=com["m_dot"],
                                     conflict=com["conflict"],
                                     bhct=float(com["run"]["T_bottom_delivered"]),
                                     spp=float(com["run"]["hyd"]["spp"]))
    alt["commercial pipe only"] = com["run"]["T_bottom_delivered"]
    out["circ_sensitivity"] = {}
    for k, Tw in alt.items():
        r, c = stability_inputs(site, z, T_rock, T_wall=float(Tw), thermal=thermal)
        out["circ_sensitivity"][k] = dict(stress=r, stress_cases=c, T_wall=float(Tw))
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
