"""C1 (stress profiles as functions of depth) and C2 (frictional admissibility).

Spec v1.1 done-when:
- C1: v1.0 inputs still reproduce the v1.0 table exactly; the United Downs
  profile gives SHmax/Shmin of about 1.97 at 5 and 12.5 km.
- C2: the cap is 3.12 at 0.6, 4.33 at 0.8 and 4.68 at 0.85; the v1.0 Cornwall
  input is flagged; the published profile passes at 5 and 12.5 km under its own
  pore pressure at mu 0.8; the cap binds on a synthetic profile, and on the
  published profile with hydrostatic Pp at mu 0.85, with the depth reported.
All fast: stability_inputs() needs no Model 1-4 calls.
"""
import numpy as np
import pytest

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, LARDERELLO, PANNONIAN, SITES
from site_evaluation import SOULTZ, StressProfile

MPa = 1.0e6
UD = CORNWALL.stress_cases[0]


def at_target(site, **kw):
    z = site.target_depth
    return se.stability_inputs(site, z, float(site.geotherm(z)), UCS=site.UCS, **kw)


# ----------------------------------------------------------------------- C1
@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
@pytest.mark.parametrize("z", [1000.0, 5058.0, 12500.0])
def test_ratio_constructor_reproduces_v10_stresses(site, z):
    """To rounding: K0*Sv_grad*z against K0*(Sv_grad*z). The v1.0 table itself is
    bit-exact because the v10 adapter calls stresses_v10 directly."""
    p = StressProfile.from_ratios(site.Sv_grad, site.K0_min, site.SHmax_over_Sv)
    Sv, Shmin, SHmax, _ = se.stresses_v10(site, z)
    assert p.Sv(z) == Sv
    assert p.Shmin(z) == pytest.approx(Shmin, rel=1e-14)
    assert p.SHmax(z) == pytest.approx(SHmax, rel=1e-14)
    assert p.Pp(z) == C.HYDROSTATIC_GRAD * z


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_v10_adapter_ignores_the_new_profiles(site):
    z = site.target_depth
    T = float(site.geotherm(z))
    ref, cases = se.stability_inputs(site, z, T, v10=True)
    Sv, Shmin, SHmax, P_fluid = se.stresses_v10(site, z)
    assert len(cases) == 1
    assert (ref["Sv"], ref["Shmin"], ref["SHmax"]) == (Sv, Shmin, SHmax)
    assert ref["breakout"] == se.breakout_v10(Shmin, SHmax, P_fluid, T)
    assert ref["cap_binds"] is False and ref["beyond_data"] is None


def test_united_downs_profile_is_reinecker_2021():
    z = 5000.0
    assert UD.Sv(z) / MPa == pytest.approx(126.375)
    assert UD.Shmin(z) / MPa == pytest.approx(69.05)
    assert UD.SHmax(z) / MPa == pytest.approx(135.85)
    assert UD.Pp(z) / MPa == pytest.approx(9.494 * (5.0 - 0.061), abs=1e-9)
    assert UD.Pp(30.0) == 0.0                      # above the static fluid level
    assert CORNWALL.stress_basis == "measured"


@pytest.mark.parametrize("z", [5000.0, 12500.0])
def test_united_downs_anisotropy_is_about_1p97(z):
    assert UD.SHmax(z) / UD.Shmin(z) == pytest.approx(1.97, abs=0.02)


def test_target_depth_beyond_the_data_is_reported():
    ref, _ = at_target(CORNWALL)
    assert ref["beyond_data"] == pytest.approx(12500.0 - 2800.0)
    ref, _ = at_target(SOULTZ)
    assert ref["beyond_data"] == pytest.approx(10700.0 - 5000.0)
    for site in (LARDERELLO, PANNONIAN):
        ref, _ = at_target(site)
        assert ref["beyond_data"] is None           # regime bounds, no data


def test_soultz_profile_is_valley_and_evans_2007():
    """Shmin ~0.54 Sv and 0.90 Sv <= SHmax <= 1.05 Sv at 5 km (the spec's
    'Shmin about 0.54 Sv, SHmax about Sv')."""
    lo, mid, hi = SOULTZ.stress_cases
    z = 5000.0
    assert lo.Shmin(z) / lo.Sv(z) == pytest.approx(0.54, abs=0.01)
    assert lo.SHmax(z) / lo.Sv(z) == pytest.approx(0.90, abs=0.01)
    assert hi.SHmax(z) / hi.Sv(z) == pytest.approx(1.05, abs=0.01)
    assert lo.SHmax(z) < mid.SHmax(z) < hi.SHmax(z)
    assert lo.z_data == (1500.0, 5000.0)


@pytest.mark.parametrize("site", [LARDERELLO, PANNONIAN], ids=["Larderello", "Pannonian"])
def test_regime_bounds_run_from_the_frictional_floor_to_sv(site):
    R = m5.frictional_cap(C.FRICTION_MU)
    floor, *_, iso = site.stress_cases
    for z in (2000.0, site.target_depth):
        Sv, Pp = floor.Sv(z), float(floor.Pp(z))
        assert (Sv - Pp) / (floor.Shmin(z) - Pp) == pytest.approx(R)
        assert floor.SHmax(z) == Sv
        assert iso.Shmin(z) == pytest.approx(Sv) and iso.SHmax(z) == Sv
    assert site.stress_basis == "regime bounds"


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_reference_case_is_the_worst_in_range(site):
    ref, cases = at_target(site)
    assert ref["breakout"]["P_need"] == max(c["breakout"]["P_need"] for c in cases)


def test_regime_bound_reference_is_the_floor_case():
    for site in (LARDERELLO, PANNONIAN):
        ref, _ = at_target(site)
        assert ref["profile"] is site.stress_cases[0]


# ----------------------------------------------------------------------- C2
@pytest.mark.parametrize("mu, cap", [(0.6, 3.12), (0.8, 4.33), (0.85, 4.68)])
def test_frictional_cap(mu, cap):
    assert m5.frictional_cap(mu) == pytest.approx(cap, abs=0.01)


def test_default_friction_is_d1():
    assert C.FRICTION_MU == 0.8
    assert C.FRICTION_MU_SENSITIVITY == (0.6, 0.85)


@pytest.mark.parametrize("z", [5058.0, 12500.0])
def test_v10_cornwall_input_is_flagged(z):
    Sv, Shmin, SHmax, _ = se.stresses_v10(CORNWALL, z)
    a = m5.stress_admissible(Sv, Shmin, SHmax, C.HYDROSTATIC_GRAD * z, C.FRICTION_MU)
    assert a["ratio"] == pytest.approx(6.23, abs=0.01)
    assert a["admissible"] is False and a["margin"] < 0


@pytest.mark.parametrize("z, ratio", [(5000.0, 4.01), (12500.0, 4.25)])
def test_published_profile_passes_under_its_own_pore_pressure(z, ratio):
    a = m5.stress_admissible(UD.Sv(z), UD.Shmin(z), UD.SHmax(z), float(UD.Pp(z)),
                             C.FRICTION_MU)
    assert a["ratio"] == pytest.approx(ratio, abs=0.01)
    assert a["admissible"] is True and a["margin"] > 0


def test_admissibility_works_in_a_normal_faulting_regime():
    """S1 is Sv there, not SHmax."""
    a = m5.stress_admissible(Sv=100 * MPa, Shmin=30 * MPa, SHmax=50 * MPa,
                             Pp=20 * MPa, mu=0.8)
    assert a["ratio"] == pytest.approx(80.0 / 10.0)
    assert a["admissible"] is False


def test_cap_does_not_bind_on_united_downs_at_d1():
    ref, _ = at_target(CORNWALL)
    assert ref["cap_binds"] is False
    assert ref["admissible"]["admissible"] is True
    assert ref["cap_depth"] > 20000.0           # ~22 km, well below the target


def test_cap_binds_on_the_published_profile_with_hydrostatic_pp_at_0p85():
    hydro = StressProfile(label="UD, hydrostatic Pp", Sv_grad=UD.Sv_grad,
                          Shmin_grad=UD.Shmin_grad, Shmin_int=UD.Shmin_int,
                          SHmax_grad=UD.SHmax_grad, SHmax_int=UD.SHmax_int,
                          z_data=UD.z_data)
    # (25.99 - 10 - R (13.21 - 10)) z = 3R - 5.9  [MPa, km], R(0.85) = 4.676
    R = m5.frictional_cap(0.85)
    z_bind = hydro.cap_depth(0.85)
    assert z_bind / 1000.0 == pytest.approx((3 * R - 5.9) / (15.99 - 3.21 * R))
    assert z_bind / 1000.0 == pytest.approx(8.30, abs=0.01)
    z = z_bind
    assert (hydro.SHmax(z) - hydro.Pp(z)) / (hydro.Shmin(z) - hydro.Pp(z)) == pytest.approx(R)

    site = _site_with(hydro, target_depth=12500.0)
    ref, _ = se.stability_inputs(site, 12500.0, 400.0, mu=0.85)
    assert ref["cap_binds"] is True
    assert ref["SHmax"] == pytest.approx(m5.SHmax_frictional_limit(ref["Shmin"], ref["Pp"], 0.85))
    assert ref["admissible"]["ratio"] == pytest.approx(R)


def test_cap_binds_on_a_synthetic_profile_only_below_its_data():
    p = StressProfile(label="synthetic", Sv_grad=26e3, Shmin_grad=14e3,
                      SHmax_grad=40e3, z_data=(0.0, 3000.0))
    z_bind = p.cap_depth(0.8)
    assert z_bind == 0.0                    # over the cap at every depth
    within = _site_with(p, target_depth=3000.0)
    below = _site_with(p, target_depth=6000.0)
    ref_in, _ = se.stability_inputs(within, 3000.0, 150.0)
    ref_out, _ = se.stability_inputs(below, 6000.0, 250.0)
    assert ref_in["cap_binds"] is False and ref_in["admissible"]["admissible"] is False
    assert ref_out["cap_binds"] is True and ref_out["admissible"]["admissible"] is True


def test_cap_depth_where_it_starts_partway_down():
    # SHmax - Pp = R (Shmin - Pp) at 4 km, over the cap below it
    R = m5.frictional_cap(0.8)
    Shmin_g, Pp_g = 15e3, 10e3
    SHmax_g = Pp_g + R * (Shmin_g - Pp_g) + 1e3
    p = StressProfile(label="t", Sv_grad=26e3, Shmin_grad=Shmin_g,
                      SHmax_grad=SHmax_g, SHmax_int=-4.0e6)
    assert p.cap_depth(0.8) == pytest.approx(4000.0)


def _site_with(profile, target_depth):
    return se.SiteProfile(
        name="synthetic", geotherm=lambda z: 15.0 + 0.03 * np.asarray(z),
        target_depth=target_depth, target_T=400.0, Sv_grad=profile.Sv_grad,
        K0_min=0.6, SHmax_over_Sv=1.0, rho_fluid_grad=C.HYDROSTATIC_GRAD,
        k_rock=3.0, E_rock=50e9, UCS=180e6, stress_cases=(profile,),
        stress_basis="measured")
