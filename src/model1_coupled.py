"""
model1_coupled.py
=================
MODEL 1, STAGE 2 -- coupled 1-D counterflow borehole heat exchanger.

Coaxial well (Eavor-like), z measured DOWNWARD along the hole from surface
(z=0) to bit (z=L). The string and hole come from a WellGeometry
(well_geometry.py); rock temperature and pressure are taken at true vertical
depth through its survey.
  * COLD fluid flows DOWN an insulated central pipe:   T_d(z)
  * HOT  fluid returns UP the annulus, against the rock: T_u(z)
  * rock wall temperature: T_rock(z) = T_surf + G*TVD(z), or a site profile

This routing is the crux of the "one loop does both" thesis: insulate the
downcomer so fluid reaches the bit COLD (tool survival + quench), while the
annulus harvests the rock's heat on the way up (the energy product).

Energy balances (per unit length), m_dot = mass flow in each leg:

  m_dot*cp_d * dT_d/dz =  UAi(z) * (T_u - T_d)
  m_dot*cp_u * dT_u/dz = -UAo(z) * (T_rock(z) - T_u) + UAi(z) * (T_u - T_d)

  UAi(z) [W/m/K] = centre<->annulus conductance per length: bore film, the
                   pipe wall's layers and the annulus film in series, weighted
                   by length between the pipe body and its bare tool joints
  UAo(z) [W/m/K] = rock<->annulus conductance per length: annulus film, any
                   casing and cement, and rock conduction in series

BCs:
  T_d(0) = T_inj
  T_u(L) = T_d(L) + Q_face/(m_dot*cp_L)          (face heat load at turnaround)

Water properties: IAPWS-95 via a precomputed table (water_table.py),
fully vectorized. P(z) hydrostatic (+1 atm surface).
"""
import numpy as np
from scipy.integrate import solve_bvp
import geo_constants as C
import water_table as wt
import well_geometry as wg
from well_geometry import WellGeometry

P_SURF = 1.0e5  # Pa, atmospheric at wellhead
P_OP = 5.0e6    # Pa, loop operating pressure for the surface heat-product calc
                # (keeps the hot return liquid; a pumped loop is pressurized)

def P_of_z(z, P_surface=P_SURF):
    return P_surface + C.HYDROSTATIC_GRAD * z


def fluid_props(T_C, P_Pa, fluid=None):
    """rho, cp, mu, k of the circulating fluid: IAPWS-95 water, or constant
    values from a dict with those keys (a drilling mud, say)."""
    if fluid is None:
        return wt.props(T_C, P_Pa)
    shape = np.broadcast(np.asarray(T_C), np.asarray(P_Pa)).shape
    return tuple(np.full(shape, float(fluid[q])) for q in ("rho", "cp", "mu", "k"))


def _cp(T_C, P_Pa, fluid=None):
    if fluid is None:
        return wt.CP(np.column_stack([np.clip(T_C, 1, 480), np.clip(P_Pa, 1e5, 220e6)]))
    return fluid_props(T_C, P_Pa, fluid)[1]


def h_dittus(m_dot, area, Dh, T_C, P_Pa, n=0.4, fluid=None):
    """Vectorized convective coefficient [W/m^2/K]. Laminar floor Nu=4.36."""
    rho, cp, mu, k = fluid_props(T_C, P_Pa, fluid)
    Gflux = m_dot / area
    Re = Gflux * Dh / mu
    Pr = cp * mu / k
    Nu = np.where(Re < 2300.0, 4.36, 0.023 * np.power(Re, 0.8) * np.power(Pr, n))
    return Nu * k / Dh


G_GRAV = 9.81
PUMP_LIMIT = 51.7e6          # Pa (7,500 psi), the limit used by Wu et al. (2025)
HOLE_CLEANING_V = 1.245      # m/s, around the drill pipe at FORGE 16B, 600 gal/min
                             # in 9-1/2-inch hole at 60-70 deg (data/materials.md)
BIT_TFA_DEFAULT = 7.759e-4   # m^2, eight 14/32-inch nozzles (FORGE 16B trial bits)
NOZZLE_CD = 0.95
SECONDS_PER_YEAR = 3.1536e7
BLEND = 1.0                 # m, width over which a change of geometry is ramped
T_EXPOSURE_FLOOR = 3600.0   # s


def exposure_from_rop(md_bit, rop_m_per_h, t_circulating=0.0):
    """Exposure time t(s) [s] when the hole was drilled at a constant rate of
    penetration and has since circulated for t_circulating seconds: the wall at
    s has been exposed since the bit passed it."""
    v = rop_m_per_h / 3600.0
    return lambda s: (md_bit - np.asarray(s, float)) / v + t_circulating


def exposure_from_history(time_s, bit_md, t_now):
    """Exposure time t(s) [s] from a drilling record: times and bit depths. The
    wall at s was first exposed when the bit first reached s, which is taken
    from the running maximum of bit depth (reaming and trips don't reset it)."""
    time_s = np.asarray(time_s, float)
    deepest = np.maximum.accumulate(np.asarray(bit_md, float))

    def t(s):
        s = np.asarray(s, float)
        # first record at or below s, interpolated from the record before it
        i = np.clip(np.searchsorted(deepest, s, side="left"), 1, len(deepest) - 1)
        d0, d1 = deepest[i - 1], deepest[i]
        frac = np.where(d1 > d0, (s - d0) / np.where(d1 > d0, d1 - d0, 1.0), 1.0)
        return t_now - (time_s[i - 1] + np.clip(frac, 0, 1) * (time_s[i] - time_s[i - 1]))
    return t


def smooth_exposure(exposure, md_bit, dz=25.0):
    """A drilling record gives an exposure time with a kink at every record,
    which the collocation solve refines without end. Resample ln t every dz
    metres and join the points with a monotone cubic."""
    from scipy.interpolate import PchipInterpolator
    s = np.linspace(0.0, md_bit, max(int(np.ceil(md_bit / dz)), 2) + 1)
    f = PchipInterpolator(s, np.log(np.maximum(exposure(s), T_EXPOSURE_FLOOR)))
    return lambda z: np.exp(f(np.clip(np.asarray(z, float), 0.0, md_bit)))


def exposure_breaks(exposure, md_bit, dz=1.0, jump=0.3):
    """Depths where an exposure-time function jumps (where the bit paused on a
    trip or a long stop): ln t changes by more than `jump` within dz metres.
    Passed to solve(breaks=...) so the solve is split there."""
    s = np.arange(0.0, md_bit + dz, dz)
    lt = np.log(np.maximum(exposure(s), T_EXPOSURE_FLOOR))
    i = np.nonzero(np.abs(np.diff(lt)) > jump)[0]
    return s[i] + 0.5 * dz


def exposure_seconds(z, t_years, exposure):
    """The time the rock term uses at each node. With exposure=None it is the
    constant t_years (floored at 1e-3 yr, as the model always has); otherwise
    exposure(z), floored at one hour.

    The rock term is quasi-steady line-source conduction, r_inf = r + 2 sqrt(alpha t).
    Its logarithm only makes sense once 2 sqrt(alpha t) is well beyond the hole
    radius; at alpha ~1e-6 m^2/s it is about 12 cm after an hour, comparable to
    a 10 cm hole. Shorter times would need the transient cylindrical-source
    solution, so they are floored instead."""
    if exposure is None:
        return np.full_like(np.asarray(z, float), max(t_years, 1e-3) * SECONDS_PER_YEAR)
    return np.maximum(exposure(z), T_EXPOSURE_FLOOR)


def conductances(z, T_d, T_u, m_dot, t_years, k_rock, geometry, exposure=None,
                 fluid=None, rhocp_rock=C.RHO_ROCK * C.CP_ROCK, geom_z=None,
                 film_mult=1.0, P_surface=P_SURF):
    """Vectorized per-length UAi (centre<->annulus), UAo (rock<->annulus).
    geom_z, if given, is where the string and hole are looked up (a point inside
    the same piece), so that nodes on a piece boundary take that piece's geometry."""
    g = geometry.at(z if geom_z is None else geom_z)
    g["tvd"] = geometry.survey.tvd(z)
    P = P_of_z(g["tvd"], P_surface)
    h_d = h_dittus(m_dot, g["A_bore"], g["Dh_bore"], T_d, P, fluid=fluid)
    h_a = film_mult * h_dittus(m_dot, g["A_ann"], g["Dh_ann"], T_u, P, fluid=fluid)

    R_in = 1.0 / (h_d * 2 * np.pi * g["r_bore"])
    R_ao = 1.0 / (h_a * 2 * np.pi * g["r_out"])
    f = g["f_joint"]
    UAi = (1 - f) / (R_in + g["R_body"] + R_ao) + f / (R_in + g["R_joint"] + R_ao)

    alpha = k_rock / rhocp_rock
    t = exposure_seconds(z, t_years, exposure)
    r_inf = g["r_rock"] + 2.0 * np.sqrt(alpha * t)
    R_aw = 1.0 / (h_a * 2 * np.pi * g["r_wall"])
    R_rock = np.log(r_inf / g["r_rock"]) / (2 * np.pi * k_rock)
    UAo = 1.0 / (R_aw + g["R_hole"] + R_rock)
    return UAi, UAo, P


def solve(m_dot=2.0, T_inj=40.0, G=C.GEOTHERM_GRADIENT, T_surf=C.SURFACE_TEMP,
          target_rock_T=C.TARGET_ROCK_TEMP, Q_face=28000.0, t_years=1.0,
          pipe=None, k_rock=C.K_ROCK, n_nodes=160, verbose=True,
          geotherm=None, target_depth=None, geometry=None, exposure=None,
          fluid=None, rhocp_rock=C.RHO_ROCK * C.CP_ROCK, friction_heat=False,
          tfa=BIT_TFA_DEFAULT, film_mult=1.0, breaks=(), init=None, P_surface=P_SURF):
    # P_surface: pressure at the wellhead. Atmospheric by default; a closed,
    # pressurized loop (P_OP) keeps a return above 100 C liquid.
    # init: a previous solve() result whose profiles seed this one (continuation).
    # breaks: extra measured depths at which to split the solve, where a
    # coefficient jumps for a reason other than geometry (exposure_breaks).
    # film_mult: multiplier on the annulus film coefficient (a calibration knob).
    # friction_heat: add the heat dissipated by friction in the bore and annulus,
    # and across the bit nozzles (total flow area tfa), to the fluid. Off by
    # default, which is how the model was first calibrated.
    # fluid: constant fluid properties (see fluid_props); default IAPWS-95 water.
    # rhocp_rock: rock volumetric heat capacity [J/m^3/K], for its diffusivity.
    # exposure: optional callable MD[m]->seconds the wall has been exposed to
    # circulation (exposure_from_rop, exposure_from_history); else t_years everywhere.
    # geotherm: optional callable TVD[m]->T_rock[degC] (e.g. a layered site profile).
    # The well is either a full WellGeometry, or one pipe type in a vertical open
    # hole of the original radius, whose length is target_depth if a geotherm is
    # given, else the depth where T_surf + G z reaches target_rock_T.
    if geometry is None:
        if pipe is None:
            raise ValueError("name a pipe type, or give a WellGeometry")
        L = target_depth if geotherm is not None else (target_rock_T - T_surf) / G
        geometry = WellGeometry.single(L, pipe)
    elif pipe is not None:
        raise ValueError("give a pipe type or a WellGeometry, not both")
    L = geometry.md_bit
    T_tvd = geotherm if geotherm is not None else (lambda d: T_surf + G * d)
    Trock = lambda zz: T_tvd(geometry.survey.tvd(zz))
    z = np.linspace(0, L, n_nodes)
    P_L = float(P_of_z(geometry.survey.tvd(L), P_surface))

    def derivs(zz, Td, Tu, geom_z=None):
        UAi, UAo, P = conductances(zz, Td, Tu, m_dot, t_years, k_rock, geometry, exposure,
                                   fluid, rhocp_rock, geom_z, film_mult, P_surface)
        cpd = _cp(Td, P, fluid)
        cpu = _cp(Tu, P, fluid)
        dTd = UAi * (Tu - Td) / (m_dot * cpd)
        dTu = (-UAo * (Trock(zz) - Tu) + UAi * (Tu - Td)) / (m_dot * cpu)
        if friction_heat:
            g = geometry.at(zz if geom_z is None else geom_z)
            dTd = dTd + _friction_gradient(m_dot, Td, P, g["A_bore"], g["Dh_bore"], fluid) / cpd
            dTu = dTu - _friction_gradient(m_dot, Tu, P, g["A_ann"], g["Dh_ann"], fluid) / cpu
        return dTd, dTu

    def dT_face(Td_L):
        cp_L = float((wt.cp if fluid is None else
                      lambda T, P: fluid_props(T, P, fluid)[1])(np.array([Td_L]), np.array([P_L]))[0])
        dT = Q_face / (m_dot * cp_L)
        if friction_heat:
            rho_L = float(np.ravel(fluid_props(np.array([Td_L]), np.array([P_L]), fluid)[0])[0])
            dT += nozzle_dp(m_dot, rho_L, tfa) / (rho_L * cp_L)
        return dT

    # Where the string or hole changes, the coefficients would jump, which the
    # collocation solve can't resolve. Across each change they are ramped
    # linearly over BLEND metres instead: within half of it either side, the
    # derivatives are blended between the geometry above and below.
    top, _ = _pieces(geometry, breaks)
    changes = top[1:]

    def odes(zz, y):
        dTd, dTu = derivs(zz, y[0], y[1])
        for b in changes:
            m = np.abs(zz - b) < 0.5 * BLEND
            if m.any():
                zm = zz[m]
                aL = derivs(zm, y[0][m], y[1][m], np.full_like(zm, b - BLEND))
                aR = derivs(zm, y[0][m], y[1][m], np.full_like(zm, b + BLEND))
                w = (zm - b) / BLEND + 0.5
                dTd[m] = (1 - w) * aL[0] + w * aR[0]
                dTu[m] = (1 - w) * aL[1] + w * aR[1]
        return np.vstack([dTd, dTu])

    def bc(ya, yb):
        return np.array([ya[0] - T_inj, yb[1] - yb[0] - dT_face(yb[0])])

    if len(changes):
        z = np.unique(np.concatenate([z, changes - 0.5 * BLEND, changes + 0.5 * BLEND]))
    y0 = np.vstack([np.interp(z, [0, L], [T_inj, T_inj + 30]),
                    np.interp(z, [0, L], [T_inj + 30, target_rock_T * 0.7])]) \
        if init is None else np.vstack([np.interp(z, init["z"], init["Td"]),
                                        np.interp(z, init["z"], init["Tu"])])
    sol = solve_bvp(odes, bc, z, y0, max_nodes=20000, tol=1e-4)
    sol_at = sol.sol

    zz = np.linspace(0, L, 400)
    Td, Tu = sol_at(zz)
    # Heat delivered to the surface plant. The loop is PRESSURIZED, so the hot
    # return stays liquid (evaluating cp at 1 atm would wrongly treat >100 C
    # return as steam, cp~2000, halving the result). Use the loop operating
    # pressure P_OP and the mean temperature of the product stream.
    cp_prod = float(wt.cp(np.array([0.5 * (Tu[0] + Td[0])]), np.array([P_OP]))[0]) \
        if fluid is None else float(fluid["cp"])
    Q_product = m_dot * cp_prod * (Tu[0] - Td[0])

    Trock_arr = np.array([Trock(zi) for zi in zz]) if geotherm is not None else Trock(zz)
    res = dict(sol=sol, z=zz, Td=Td, Tu=Tu, Trock=Trock_arr, L=L, geometry=geometry,
               T_bottom_delivered=Td[-1], T_return_surface=Tu[0],
               Q_product=Q_product, m_dot=m_dot,
               # a solution outside the water table's range (above 480 C) is spurious
               success=bool(sol.success and max(Td.max(), Tu.max()) < 480.0),
               P_bottom=P_L)
    if verbose:
        _print(res, T_inj, Q_face, t_years if exposure is None else "by depth")
    return res


def colebrook_smooth(Re, n_iter=30):
    """Darcy friction factor for a smooth pipe (Colebrook at zero roughness)."""
    Re = np.asarray(Re, float)
    f = 0.316 * Re ** -0.25
    for _ in range(n_iter):
        f = (-2.0 * np.log10(2.51 / (Re * np.sqrt(f)))) ** -2
    return f


def friction_factor(Re, model="blasius"):
    """Darcy friction factor: 64/Re below Re 2300, else Blasius (or smooth
    Colebrook as a check)."""
    Re = np.asarray(Re, float)
    turb = colebrook_smooth(np.maximum(Re, 2300.0)) if model == "colebrook" \
        else 0.316 * np.maximum(Re, 1.0) ** -0.25
    return np.where(Re < 2300.0, 64.0 / np.maximum(Re, 1e-12), turb)


def nozzle_dp(m_dot, rho, tfa=BIT_TFA_DEFAULT, cd=NOZZLE_CD):
    """Pressure drop across the bit nozzles [Pa]."""
    Q = m_dot / rho
    return rho * Q ** 2 / (2 * cd ** 2 * tfa ** 2)


def _friction_gradient(m_dot, T, P, area, Dh, fluid=None):
    """Friction pressure gradient divided by density, |dp/dz| / rho [J/kg/m]:
    the heat dissipated per unit mass of fluid per metre. The circulating
    fluid is liquid, so its properties are taken at no less than the loop
    pressure P_OP: otherwise a solver iterate above 100 C near the surface
    picks up steam's density and viscosity, and the friction heat runs away."""
    rho, _, mu, _ = fluid_props(np.clip(T, 1, 480), np.maximum(P, P_OP), fluid)
    v = m_dot / (rho * area)
    f = friction_factor(rho * v * Dh / mu)
    return f * v ** 2 / (2 * Dh)


def _pieces(geometry, breaks=()):
    """Measured-depth pieces of constant geometry: cut at segment ends and
    hole intervals, and at any extra breaks."""
    L = geometry.md_bit
    cuts = np.unique(np.concatenate([[0.0, L], geometry._seg_ends,
                                     [h.md_top for h in geometry.hole],
                                     np.asarray(breaks, float)]))
    cuts = cuts[(cuts >= 0) & (cuts <= L)]
    cuts = cuts[np.concatenate([[True], np.diff(cuts) > 1e-6])]
    return cuts[:-1], cuts[1:]


def hydraulics(m_dot, geometry, result=None, T_avg=80.0, tfa=BIT_TFA_DEFAULT,
               friction="blasius", bit=True, fluid=None, P_surface=P_SURF):
    """Pressure losses along the string and annulus, bit nozzle loss,
    standpipe pressure and annular velocity.

    With a solve() result, water properties are taken at each piece's local
    bore or annulus temperature and pressure, and the standpipe pressure
    includes the buoyancy of the hot annulus against the cold bore. Without
    one, properties are at T_avg and mid-depth pressure throughout.

      SPP = dP_string + dP_bit + dP_annulus + integral (rho_ann - rho_bore) g dTVD
    """
    top, bot = _pieces(geometry)
    mid = 0.5 * (top + bot)
    g = geometry.at(mid)
    dtvd = geometry.survey.tvd(bot) - geometry.survey.tvd(top)
    if result is None:
        P = P_of_z(geometry.survey.tvd(geometry.md_bit / 2), P_surface)
        rho, _, mu, _ = (float(np.ravel(x)[0]) for x in
                         fluid_props(np.array([T_avg]), np.array([P]), fluid))
        rho_b = rho_a = np.full_like(mid, rho)
        mu_b = mu_a = np.full_like(mid, mu)
        rho_bit = rho
    else:
        P = P_of_z(g["tvd"], P_surface)
        T_b = np.interp(mid, result["z"], result["Td"])
        T_a = np.interp(mid, result["z"], result["Tu"])
        rho_b, _, mu_b, _ = fluid_props(T_b, P, fluid)
        rho_a, _, mu_a, _ = fluid_props(T_a, P, fluid)
        rho_bit = float(rho_b[-1])

    def loss(rho, mu, area, Dh):
        v = m_dot / (rho * area)
        f = friction_factor(rho * v * Dh / mu, friction)
        return f * ((bot - top) / Dh) * 0.5 * rho * v ** 2, v

    dp_bore, _ = loss(rho_b, mu_b, g["A_bore"], g["Dh_bore"])
    dp_ann, v_ann = loss(rho_a, mu_a, g["A_ann"], g["Dh_ann"])
    dp_bit = nozzle_dp(m_dot, rho_bit, tfa) if bit else 0.0
    dp_buoy = float(np.sum((rho_a - rho_b) * G_GRAV * dtvd))
    spp = float(dp_bore.sum() + dp_ann.sum()) + dp_bit + dp_buoy
    return dict(
        pieces=[dict(md_top=a, md_bottom=b, dp_bore=x, dp_ann=y, v_ann=v)
                for a, b, x, y, v in zip(top, bot, dp_bore, dp_ann, v_ann)],
        dp_string=float(dp_bore.sum()), dp_annulus=float(dp_ann.sum()),
        dp_bit=dp_bit, dp_buoyancy=dp_buoy, spp=spp,
        over_pump_limit=spp > PUMP_LIMIT,
        v_ann_min=float(v_ann.min()), cleans_hole=float(v_ann.min()) >= HOLE_CLEANING_V,
        pump_power=spp * m_dot / float(np.mean(rho_b)))


def pump_power(m_dot, L=None, T_avg=80.0, geometry=None):
    """Friction power and pressure drop down the bore and up the annulus, at
    T_avg and mid-depth pressure, without the bit. geometry=None is the
    original coaxial layout over length L."""
    if geometry is None:
        geometry = WellGeometry.single(L, wg.LEGACY_VACUUM)  # only the radii matter
    h = hydraulics(m_dot, geometry, T_avg=T_avg, bit=False)
    dP = h["dp_string"] + h["dp_annulus"]
    P = P_of_z(geometry.survey.tvd(geometry.md_bit / 2))
    rho = float(wt.props(np.array([T_avg]), np.array([P]))[0][0])
    return dP * (m_dot / rho), dP


def _print(r, T_inj, Q_face, t_years):
    print("=" * 72)
    print("COUPLED 1-D COUNTERFLOW BOREHOLE  (Model 1, stage 2)")
    print("=" * 72)
    print(f"  Depth L                 : {r['L']/1000:.1f} km   (bottom rock {r['Trock'][-1]:.0f} C)")
    print(f"  Mass flow (per leg)     : {r['m_dot']:.2f} kg/s")
    pipes = ", ".join(dict.fromkeys(seg.pipe.name for seg in r["geometry"].string))
    print(f"  Injection temp          : {T_inj:.0f} C ; pipe: {pipes} ; t={t_years}{' yr' if isinstance(t_years, (int, float)) else ''}")
    print(f"  Face heat load Q_face   : {Q_face/1000:.1f} kW")
    print(f"  BVP converged           : {r['success']}  (P_bottom {r['P_bottom']/1e6:.0f} MPa)")
    print("-" * 72)
    pp, dP = pump_power(r['m_dot'], geometry=r['geometry'])
    ok = r['T_bottom_delivered'] < C.BHA_SURVIVAL_TEMP
    print(f"  >> Bit-delivered temp T_d(L) : {r['T_bottom_delivered']:6.1f} C   "
          f"(ceiling {C.BHA_SURVIVAL_TEMP:.0f} C)  {'OK' if ok else 'FAIL'}")
    print(f"  >> Return temp at surface    : {r['T_return_surface']:6.1f} C")
    print(f"  >> Thermal power exported    : {r['Q_product']/1e6:6.2f} MW_th")
    print(f"  >> Pump (friction) power     : {pp/1e3:6.1f} kW   (dP~{dP/1e6:.1f} MPa)")
    if pp > 0:
        print(f"  >> Thermal/pump ratio        : {r['Q_product']/max(pp,1):6.0f}x")
    print("=" * 72)


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    solve(m_dot=2.0, T_inj=40.0, Q_face=28000.0, t_years=1.0, pipe=wg.LEGACY_VACUUM)
    print("\n### mass-flow sweep (legacy vacuum tubing, t=1yr) ###")
    print(f"{'m_dot':>7}{'T_bit[C]':>10}{'T_surf[C]':>11}{'Q_MWth':>9}{'pump_kW':>9}{'conv':>6}")
    for md in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        r = solve(m_dot=md, t_years=1.0, pipe=wg.LEGACY_VACUUM, verbose=False)
        pp, _ = pump_power(md, geometry=r["geometry"])
        print(f"{md:7.1f}{r['T_bottom_delivered']:10.1f}{r['T_return_surface']:11.1f}"
              f"{r['Q_product']/1e6:9.2f}{pp/1e3:9.1f}{str(r['success']):>6}")
    print("\n### pipe type (m_dot=2 kg/s) ###")
    print(f"{'pipe':>32}{'T_bit[C]':>10}{'Q_MWth':>9}")
    for p in wg.PIPE_TYPES.values():
        r = solve(m_dot=2.0, pipe=p, verbose=False)
        print(f"{p.name:>32}{r['T_bottom_delivered']:10.1f}{r['Q_product']/1e6:9.2f}")
    print("\n### operating-time sensitivity (rock cooldown, m_dot=2 kg/s) ###")
    print(f"{'years':>8}{'T_bit[C]':>10}{'Q_MWth':>9}")
    for ty in (0.1, 1.0, 5.0, 30.0):
        r = solve(m_dot=2.0, t_years=ty, pipe=wg.LEGACY_VACUUM, verbose=False)
        print(f"{ty:8.1f}{r['T_bottom_delivered']:10.1f}{r['Q_product']/1e6:9.2f}")
