"""
model5_convergence_confinement.py
================================
Does a COLD LOAD-BEARING SHELL + fluid support keep a deep hole open even where
the elastic wall stress exceeds rock strength? Breakout != collapse: the rock
yields into a plastic annulus that converges to a finite displacement, and a
support pressure stabilises it. This is the tunnelling 'convergence-confinement'
method (Ground Reaction Curve), applied here with TEMPERATURE-DEPENDENT strength
and stiffness so the cooling enters physically.

Model (Duncan Fama 1993 / Hoek, circular opening, hydrostatic far field p0,
Mohr-Coulomb elastic-perfectly-plastic):

  k     = (1+sin phi)/(1-sin phi)
  p_cr  = (2 p0 - sigma_cm)/(1 + k)            critical support pressure
  if p_i >= p_cr   (ELASTIC):
        u_i = a (p0 - p_i)(1+nu)/E
        r_p = a
  if p_i <  p_cr   (PLASTIC annulus radius r_p, wall convergence u_i):
        r_p = a [ 2(p0(k-1)+sigma_cm) / ((1+k)((k-1)p_i+sigma_cm)) ]^(1/(k-1))
        u_i = a(1+nu)/E [ 2(1-nu)(p0-p_cr)(r_p/a)^2 - (1-2nu)(p0-p_i) ]

Cooling enters via sigma_cm(T) (strength rises as rock cools) and E(T). The
cold-rock GRC is valid out to the cooled-shell radius; we check self-consistency
(is the plastic zone r_p contained within the shell the dwell has cooled?).

KEY CORRECTION over Model 4(E): proper Mohr-Coulomb confines the wall via the
RADIAL stress (the fluid pressure), through the (1+k) factor -- so the simple
'sigma_theta > UCS' test in 4(E) was too pessimistic. Fluid support does a lot
of the work; cooling EXTENDS the stable envelope and cuts convergence.
"""
import numpy as np
import geo_constants as C

ALPHA = C.K_ROCK / (C.RHO_ROCK * C.CP_ROCK)
PHI = np.radians(35.0)
KMC = (1 + np.sin(PHI)) / (1 - np.sin(PHI))


def sigma_cm(T_C, UCS=None):
    """Rock-mass compressive strength [Pa], declining with T (cooling -> stronger).

    UCS [Pa] is the strength at 25 C; None means the global C.UCS, read at call
    time. Site runs pass their own SiteProfile.UCS.
    """
    if UCS is None:
        UCS = C.UCS
    return UCS * max(1.0 - 0.0009 * max(T_C - 25.0, 0.0), 0.4)


def hoop_stress(theta_deg, SHmax, Shmin, Pw, dsigma_T=0.0):
    """Kirsch hoop stress [Pa] at the wall of a vertical hole.

    theta is measured from the SHmax azimuth, so the compressive peak,
    3 SHmax - Shmin - Pw, is at 90 deg (the Shmin azimuth, where breakouts
    form) and the minimum, 3 Shmin - SHmax - Pw, at 0 deg. dsigma_T is the
    thermal hoop stress from cooling the wall (site_evaluation.
    thermal_hoop_stress); zero when that is off. Vectorised over theta.
    """
    th = np.radians(np.asarray(theta_deg, dtype=float))
    return SHmax + Shmin - 2.0 * (SHmax - Shmin) * np.cos(2.0 * th) - Pw + dsigma_T


def _width_threshold(SHmax, Shmin, Pw, Pp, T_wall, UCS, dsigma_T):
    """Effective failure margin at the Shmin azimuth, and the anisotropy term.

    Breakout where sigma_theta(theta) - Pp > sigma_cm + KMC (Pw - Pp). Writing
    phi = theta - 90 deg, sigma_theta = SHmax + Shmin + 2 (SHmax - Shmin) cos 2phi
    - Pw + dsigma_T, so the wall fails where cos 2phi > -a / d, with
    a = SHmax + Shmin - Pw + dsigma_T - Pp - limit and d = 2 (SHmax - Shmin).
    """
    limit = sigma_cm(T_wall, UCS) + KMC * (Pw - Pp)
    a = SHmax + Shmin - Pw + dsigma_T - Pp - limit
    return a, 2.0 * (SHmax - Shmin)


def breakout_width(SHmax, Shmin, Pw, Pp, T_wall, UCS=None, dsigma_T=0.0):
    """Width [deg] of one breakout lobe, centred on the Shmin azimuth.

    Closed form of the effective-stress Mohr-Coulomb criterion (see
    _width_threshold) around the circumference: the lobe spans
    |phi| < arccos(-a/d) / 2, so its full width is arccos(-a/d). 0 when the
    wall holds everywhere, 180 when it fails all round (two lobes meeting).
    """
    a, d = _width_threshold(SHmax, Shmin, Pw, Pp, T_wall, UCS, dsigma_T)
    if d < 0.0:
        raise ValueError("SHmax < Shmin: the horizontal stresses are mislabelled")
    if d == 0.0:                         # equal horizontal stresses: all or nothing
        return 180.0 if a > 0 else 0.0
    c = -a / d
    if c >= 1.0:
        return 0.0
    if c <= -1.0:
        return 180.0
    return float(np.degrees(np.arccos(c)))


def mud_for_width(SHmax, Shmin, Pp, T_wall, W_deg, UCS=None, dsigma_T=0.0):
    """Lightest Pw [Pa] that keeps the lobe at or below W_deg: the lower bound
    of the mud window.

    Width <= W  <=>  -a/d >= cos W  <=>  Pw >= (SHmax + Shmin + dsigma_T
    + (KMC - 1) Pp - sigma_cm + d cos W) / (1 + KMC). W = 0 gives the v1.0
    full-suppression pressure.
    """
    d = 2.0 * (SHmax - Shmin)
    return (SHmax + Shmin + dsigma_T + (KMC - 1.0) * Pp - sigma_cm(T_wall, UCS)
            + d * np.cos(np.radians(W_deg))) / (1.0 + KMC)


def frictional_cap(mu):
    """Largest effective S1/S3 a crust of optimally oriented faults with
    friction mu can sustain: ((1 + mu^2)^0.5 + mu)^2 (Jaeger & Cook)."""
    return (np.sqrt(1.0 + mu ** 2) + mu) ** 2


# Relative tolerance for "at the cap": regime-bound profiles are built exactly on
# it, and rounding must not tip them over.
CAP_RTOL = 1e-9


def stress_admissible(Sv, Shmin, SHmax, Pp, mu):
    """Frictional admissibility of a stress state.

    Compares effective S1/S3 with frictional_cap(mu). Returns ratio, cap,
    margin (cap - ratio, negative when over) and the admissible flag. Works on
    any regime: S1 and S3 are the largest and smallest of the three.
    """
    S1, S3 = max(Sv, Shmin, SHmax), min(Sv, Shmin, SHmax)
    cap = frictional_cap(mu)
    ratio = (S1 - Pp) / (S3 - Pp) if S3 > Pp else np.inf
    return dict(ratio=ratio, cap=cap, margin=cap - ratio,
                admissible=bool(ratio <= cap * (1 + CAP_RTOL)), mu=mu)


def SHmax_frictional_limit(Shmin, Pp, mu):
    """Largest SHmax the frictional cap allows given Shmin: Pp + R(mu)(Shmin - Pp)."""
    return Pp + frictional_cap(mu) * (Shmin - Pp)


def E_of_T(T_C):
    return max(C.E_ROCK * (1.0 - 0.0007 * max(T_C - 25.0, 0.0)), 0.3 * C.E_ROCK)


def grc(p_i, p0, T_wall, a=C.CAVITY_RADIUS, nu=C.NU_ROCK, UCS=None):
    """Ground reaction: returns (u_i [m] wall convergence, r_p/a, p_cr [Pa], regime)."""
    scm = sigma_cm(T_wall, UCS)
    E = E_of_T(T_wall)
    p_cr = (2 * p0 - scm) / (1 + KMC)
    if p_i >= p_cr:
        u_i = a * (p0 - p_i) * (1 + nu) / E
        return u_i, 1.0, p_cr, "elastic"
    rp_a = (2 * (p0 * (KMC - 1) + scm) /
            ((1 + KMC) * ((KMC - 1) * p_i + scm))) ** (1.0 / (KMC - 1))
    u_i = a * (1 + nu) / E * (2 * (1 - nu) * (p0 - p_cr) * rp_a ** 2
                              - (1 - 2 * nu) * (p0 - p_i))
    return u_i, rp_a, p_cr, "PLASTIC"


def dwell_to_cool(radius_m, a=C.CAVITY_RADIUS):
    """Time [s] for the conduction cold front (2 sqrt(alpha t)) to reach a radius."""
    delta = max(radius_m - a, 0.0)
    return (delta / 2.0) ** 2 / ALPHA


def fmt_t(s):
    if s < 86400: return f"{s/3600:.1f} h"
    if s < 3.156e7: return f"{s/86400:.1f} d"
    return f"{s/3.156e7:.1f} yr"


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    G = 0.035
    a = C.CAVITY_RADIUS

    print("=" * 80)
    print("(A) Is the hole even plastic? Fluid support vs critical pressure, by K0")
    print("    p_i = hydrostatic fluid column; p0 = K0 * lithostatic")
    print("=" * 80)
    print(f"{'z[km]':>6}{'T_rock':>7}{'p_i[MPa]':>9}{'K0':>5}{'p0[MPa]':>9}"
          f"{'p_cr_hot':>9}{'p_cr_cold':>10}{'regime(hot/cold)':>18}")
    for zkm in (10, 12, 14):
        z = zkm * 1000.0; T_rock = C.SURFACE_TEMP + G * z
        p_i = C.HYDROSTATIC_GRAD * z
        for K0 in (0.7, 1.0, 1.3):
            p0 = K0 * C.LITHOSTATIC_GRAD * z
            _, _, pcr_h, rh = grc(p_i, p0, T_rock, a)
            _, _, pcr_c, rc = grc(p_i, p0, 200.0, a)
            print(f"{zkm:6.0f}{T_rock:7.0f}{p_i/1e6:9.0f}{K0:5.1f}{p0/1e6:9.0f}"
                  f"{pcr_h/1e6:9.0f}{pcr_c/1e6:10.0f}{rh+'/'+rc:>18}")
    print("  -> hydrostatic fluid alone keeps K0<~0.9 ELASTIC even hot; compressional")
    print("     crust (K0>=1) goes plastic -> that is where cooling/shell/mud-weight earn")
    print("     their keep. (Corrects Model 4E, which ignored fluid confinement.)")

    print("\n" + "=" * 80)
    print("(B) Compressional case (K0=1.3, z=12 km): convergence HOT vs COOLED")
    print("    vs support pressure p_i  (the Ground Reaction Curve)")
    print("=" * 80)
    z = 12000.0; T_rock = C.SURFACE_TEMP + G * z
    p0 = 1.3 * C.LITHOSTATIC_GRAD * z
    print(f"  p0={p0/1e6:.0f} MPa, T_rock={T_rock:.0f} C")
    print(f"{'p_i[MPa]':>9}{'  | HOT: u[mm] r_p/a regime':<30}{'| COOLED200: u[mm] r_p/a regime':<32}")
    for pi_frac in (0.30, 0.40, 0.50, 0.60, 0.70):
        p_i = pi_frac * C.LITHOSTATIC_GRAD * z
        uh, rph, _, regh = grc(p_i, p0, T_rock, a)
        uc, rpc, _, regc = grc(p_i, p0, 200.0, a)
        print(f"{p_i/1e6:9.0f}  | {uh*1000:7.1f} {rph:5.2f} {regh:<8}"
              f"| {uc*1000:7.1f} {rpc:5.2f} {regc:<8}")
    print("  -> cooling raises strength+stiffness: smaller plastic zone, less convergence,")
    print("     and pushes the elastic threshold to lower support pressure.")

    print("\n" + "=" * 80)
    print("(C) The cold-shell ratchet: dwell to cool out to the plastic radius r_p")
    print("    (so the strong cold rock contains the yielded zone) -- K0=1.3, z=12km")
    print("=" * 80)
    print(f"{'p_i[MPa]':>9}{'r_p (hot)[m]':>13}{'shell needed[m]':>16}{'dwell':>10}"
          f"{'u_cooled[mm]':>13}")
    for pi_frac in (0.40, 0.50, 0.60):
        p_i = pi_frac * C.LITHOSTATIC_GRAD * z
        uh, rph, _, _ = grc(p_i, p0, T_rock, a)
        uc, rpc, _, regc = grc(p_i, p0, 200.0, a)
        r_p_m = rph * a
        shell = (rpc * a) if regc == "PLASTIC" else (rph * a)  # cool out to contain yield
        t = dwell_to_cool(shell)
        print(f"{p_i/1e6:9.0f}{r_p_m:13.2f}{shell:16.2f}{fmt_t(t):>10}{uc*1000:13.1f}")
    print("  -> plastic zones are sub-metre; cooling out to contain them takes hours-days,")
    print("     well within a drilling dwell. The cold shell as a contained, self-")
    print("     supporting yielded annulus is GEOMECHANICALLY CREDIBLE here.")

    print("\n" + "=" * 80)
    print("(D) Verdict band: deepest stable depth vs stress regime (hydrostatic fluid)")
    print("=" * 80)
    print(f"{'K0':>5}{'cooled? ':>9}{'max stable depth (u<1% radius)':>34}")
    for K0 in (0.7, 1.0, 1.3):
        for Tlab, Tw in (("hot", None), ("cooled200", 200.0)):
            zmax = 0
            for zkm in np.arange(6, 25, 0.5):
                z = zkm * 1000.0
                Tr = C.SURFACE_TEMP + G * z
                Tw_use = Tr if Tw is None else Tw
                p_i = C.HYDROSTATIC_GRAD * z
                p0 = K0 * C.LITHOSTATIC_GRAD * z
                u, rp, _, _ = grc(p_i, p0, Tw_use, a)
                if u < 0.01 * a:
                    zmax = zkm
                else:
                    break
            print(f"{K0:5.1f}{Tlab:>9}{zmax:>30.1f} km")

    print("\n" + "=" * 80)
    print("(E) The honest hard part: STRESS ANISOTROPY drives localized breakout")
    print("    Vertical hole: max hoop stress at the min-horizontal azimuth:")
    print("       sigma_theta = 3*S_H - S_h - P_i   (S_H=max, S_h=min horizontal)")
    print("    Mohr-Coulomb breakout when sigma_theta > k*P_i + sigma_cm.")
    print("    Aniso = S_H/S_h. z=12 km, S_h = K0*sigma_v, hydrostatic P_i.")
    print("=" * 80)
    z = 12000.0; T_rock = C.SURFACE_TEMP + G * z
    P_i = C.HYDROSTATIC_GRAD * z
    sig_v = C.LITHOSTATIC_GRAD * z
    print(f"{'K0':>5}{'Aniso':>7}{'sig_theta[MPa]':>15}{'MC limit hot':>14}"
          f"{'MC limit cold':>15}{'breakout?':>20}")
    for K0 in (0.7, 1.0):
        S_h = K0 * sig_v
        for aniso in (1.0, 1.3, 1.6):
            S_H = aniso * S_h
            sth = 3 * S_H - S_h - P_i
            lim_hot = KMC * P_i + sigma_cm(T_rock)
            lim_cold = KMC * P_i + sigma_cm(200.0)
            v = ("breaks hot&cold" if sth > lim_cold else
                 "breaks hot only" if sth > lim_hot else "stable")
            print(f"{K0:5.1f}{aniso:7.1f}{sth/1e6:15.0f}{lim_hot/1e6:14.0f}"
                  f"{lim_cold/1e6:15.0f}{v:>20}")
    print("  -> anisotropy is the real breakout driver: at Aniso~1.6 the wall breaks out")
    print("     even cold -> needs higher mud weight or accepts a breakout-widened hole.")
    print("     Cooling buys one strength tier (a 'breaks hot only' band becomes stable);")
    print("     it helps but does NOT license arbitrary anisotropy. The defensible claim:")
    print("     in moderate-anisotropy, low-K0 crust, cooling + fluid support keep a")
    print("     stable hole several km below the conventional ductile limit.")


# ------------------------------------------------- stresses behind the wall face
def kirsch(r, theta_deg, SHmax, Shmin, Pw, a):
    """Stresses around a vertical hole of radius a at radius r, angle theta from
    the SHmax azimuth, with wellbore pressure Pw; compression positive (Jaeger,
    Cook and Zimmerman 2007; Zoback 2007). Returns sigma_r, sigma_theta, tau_rt.
    Broadcasts over r and theta."""
    r = np.asarray(r, float)
    th = np.radians(np.asarray(theta_deg, float))
    q2, q4 = (a / r) ** 2, (a / r) ** 4
    s, d = 0.5 * (SHmax + Shmin), 0.5 * (SHmax - Shmin)
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    sr = s * (1 - q2) + d * (1 - 4 * q2 + 3 * q4) * c2 + Pw * q2
    st = s * (1 + q2) - d * (1 + 3 * q4) * c2 - Pw * q2
    tau = -d * (1 + 2 * q2 - 3 * q4) * s2
    return sr, st, tau


def thermal_stresses(r, dT, a, thermo):
    """Plane-strain thermoelastic stresses [Pa] at radii r (from a outwards)
    for an axisymmetric temperature change dT(r) = T(r) - T_rock, with the inner
    boundary free and the rock unbounded (Timoshenko and Goodier 1970, section
    150), compression positive. thermo = E alpha / (1 - nu), held at its value
    at the rock temperature. Returns dsigma_r, dsigma_theta.

      dsigma_r     =  thermo (1/r^2) int_a^r dT s ds
      dsigma_theta =  thermo [dT(r) - (1/r^2) int_a^r dT s ds]
    """
    r = np.asarray(r, float)
    dT = np.asarray(dT, float)
    f = dT * r
    integral = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(r))])
    I = integral / r ** 2
    return thermo * I, thermo * (dT - I)


def failure_margin(sr, st, tau, Pp, T, UCS):
    """Mohr-Coulomb margin [Pa] in effective stress, in the plane of the hole:
    sigma_1' - (sigma_cm(T) + KMC sigma_3'); failure where it is positive."""
    m = 0.5 * (sr + st) - Pp
    rad = np.sqrt((0.5 * (st - sr)) ** 2 + tau ** 2)
    s1, s3 = m + rad, m - rad
    cm = np.vectorize(lambda t: sigma_cm(t, UCS))(T) if np.ndim(T) else sigma_cm(T, UCS)
    return s1 - (cm + KMC * s3)


def breakout_behind_wall(SHmax, Shmin, Pw, Pp, a, r_field, T_field, T_rock, UCS, thermo,
                         r_out=2.0, n_r=81, d_theta=0.1):
    """Failure over the rock from the wall out to r_out * a, with the stresses
    from Kirsch and the thermal field, and the strength at the local
    temperature. In elastic Mohr-Coulomb the failed region reaches wider just
    behind the wall than at it, so "width" here is reported, not used for the
    verdict (see breakout_with_skin). T_field(r_field) is the temperature around the hole (T_rock
    beyond it). Returns the width [deg] (the angular extent, about the Shmin
    azimuth, of failure at any radius), the width at the wall alone, the depth
    of the failed zone at the Shmin azimuth [m beyond the wall], and the
    radius at which the failure margin is greatest there."""
    r = a * np.geomspace(1.0, r_out, n_r)
    T = np.interp(r, r_field, T_field, right=T_rock)
    dsr, dst = thermal_stresses(r, T - T_rock, a, thermo)
    # angles from the Shmin azimuth (90 deg from SHmax), 0 to 90
    phi = np.arange(0.0, 90.0 + d_theta / 2, d_theta)
    R, PH = np.meshgrid(r, phi, indexing="ij")
    sr, st, tau = kirsch(R, 90.0 - PH, SHmax, Shmin, Pw, a)
    sr, st = sr + dsr[:, None], st + dst[:, None]
    marg = failure_margin(sr, st, tau, Pp, np.broadcast_to(T[:, None], R.shape), UCS)
    fails = marg > 0
    def extent(col):
        # contiguous failed angles from the Shmin azimuth
        if not col[0]:
            return 0.0
        k = np.argmin(col) if not col.all() else len(col)
        return 2 * phi[k - 1] if k > 0 else 0.0
    any_r = fails.any(axis=0)
    width = extent(any_r)
    width_wall = extent(fails[0])
    at_shmin = fails[:, 0]
    depth = float(r[np.max(np.nonzero(at_shmin))] - a) if at_shmin.any() else 0.0
    r_first = float(r[np.argmax(marg[:, 0])])
    return dict(width=width, width_wall=width_wall, depth=depth, r_first=r_first,
                r=r, phi=phi, margin=marg)


def breakout_with_skin(SHmax, Shmin, Pw, Pp, a, r_field, T_field, T_rock, UCS, thermo):
    """Breakout width with the cooling credited only as deep as the breakout
    would reach.

    The width is judged at the wall, as an image log measures it and as the
    90 deg limit and the UD-1 calibration assume. Elastic Mohr-Coulomb gives
    failure behind the wall that reaches wider than the wall's breakout (wings
    just behind its edges), so the extent at any radius is not a measure of
    breakout width. What the cooled skin changes is the rock the breakout grows
    into. So: find the depth the uncooled breakout's failure reaches at the
    Shmin azimuth, and apply the wall check with the temperature (strength and
    thermal stress) of the rock at that depth. A skin thicker than that keeps
    its benefit; a thin or reheated one is judged by the warmer rock behind it.
    With no cooling this is the wall width exactly.

    Returns the width, the wall-only width (the v1.2.2 measure), the reference
    depth and the temperature used, and the behind-wall failure analysis of
    the cooled field (failed depth at the Shmin azimuth, where failure starts)."""
    uncooled = breakout_behind_wall(SHmax, Shmin, Pw, Pp, a, np.array([a, 1e3]),
                                    np.array([T_rock, T_rock]), T_rock, UCS, thermo)
    r_ref = a + uncooled["depth"]
    T_ref = float(np.interp(r_ref, r_field, T_field, right=T_rock))
    T_wall = float(np.interp(a, r_field, T_field, right=T_rock))
    width = breakout_width(SHmax, Shmin, Pw, Pp, T_ref, UCS, thermo * (T_ref - T_rock))
    width_wall = breakout_width(SHmax, Shmin, Pw, Pp, T_wall, UCS, thermo * (T_wall - T_rock))
    cooled = breakout_behind_wall(SHmax, Shmin, Pw, Pp, a, r_field, T_field, T_rock, UCS, thermo)
    return dict(width=width, width_wall=width_wall, depth_ref=uncooled["depth"], T_ref=T_ref,
                T_wall=T_wall, failed_depth=cooled["depth"], r_first=cooled["r_first"],
                behind=cooled)
