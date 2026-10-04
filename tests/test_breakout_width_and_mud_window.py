"""Breakout width, the mud window, and the GO / CONDITIONAL / NO-GO verdict.

Width: the closed form matches brute-force sampling and a hand-worked Kirsch
case, and equal horizontal stresses give all or nothing (0 or 180 deg).
Window and verdict: each verdict is reachable from synthetic inputs, the
bounds are what they claim, tensile initiation is reported without bounding
the window, and every site and case gets a window.
"""
import numpy as np
import pytest

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, SITES

MPa = 1.0e6


def width_by_sampling(SH, Sh, Pw, Pp, T, UCS=None, dT=0.0):
    """Breakout width by brute force: 0.5 deg steps over one half-turn."""
    th = np.arange(0.0, 180.0, 0.5)
    fails = (m5.hoop_stress(th, SH, Sh, Pw, dT) - Pp
             > m5.sigma_cm(T, UCS) + m5.KMC * (Pw - Pp))
    return 0.5 * fails.sum()


# ------------------------------------------------------------ breakout width
def test_equal_horizontal_stresses_give_0_or_180():
    S, Pp, T = 100 * MPa, 40 * MPa, 200.0
    assert m5.breakout_width(S, S, 90 * MPa, Pp, T) == 0.0
    assert m5.breakout_width(S, S, 41 * MPa, Pp, T, UCS=20 * MPa) == 180.0


def test_hand_checked_kirsch_case():
    """SHmax 150, Shmin 80, Pw 50, Pp 40 MPa, sigma_cm 100 MPa (UCS 100 MPa at
    25 C). Limit = 100 + KMC*10 = 136.9 MPa effective; hoop - Pp at phi from the
    Shmin azimuth = 230 + 140 cos 2phi - 90 = 140 + 140 cos 2phi. Fails where
    cos 2phi > -3.1/140, i.e. |phi| < 45.63 deg: width 91.27 deg."""
    K = m5.KMC
    expected = np.degrees(np.arccos(-(140.0 - 100.0 - 10.0 * K) / 140.0))
    assert expected == pytest.approx(91.27, abs=0.01)
    w = m5.breakout_width(150 * MPa, 80 * MPa, 50 * MPa, 40 * MPa, 25.0, UCS=100 * MPa)
    assert w == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("Pw_MPa", [30.0, 50.0, 70.0, 90.0, 110.0])
def test_closed_form_matches_half_degree_sampling(Pw_MPa):
    args = (260 * MPa, 130 * MPa, Pw_MPa * MPa, 95 * MPa, 200.0)
    assert m5.breakout_width(*args, UCS=180 * MPa) == pytest.approx(
        width_by_sampling(*args, UCS=180 * MPa), abs=1.0)


def test_width_shrinks_as_mud_weight_rises():
    w = [m5.breakout_width(260 * MPa, 130 * MPa, P * MPa, 95 * MPa, 200.0)
         for P in range(90, 200, 10)]
    assert all(a >= b for a, b in zip(w, w[1:]))
    assert w[0] > 0.0 and w[-1] == 0.0


def test_width_zero_threshold_is_the_full_suppression_pressure():
    """mud_for_width at 0 deg is breakout_eff's P_need."""
    SH, Sh, Pp, T = 260 * MPa, 130 * MPa, 95 * MPa, 200.0
    P0 = m5.mud_for_width(SH, Sh, Pp, T, 0.0)
    assert P0 == pytest.approx(se.breakout_eff(Sh, SH, 100 * MPa, Pp, T)["P_need"])


@pytest.mark.parametrize("W", [30.0, 60.0, 90.0, 120.0])
def test_mud_for_width_gives_that_width(W):
    SH, Sh, Pp, T = 260 * MPa, 130 * MPa, 95 * MPa, 200.0
    P = m5.mud_for_width(SH, Sh, Pp, T, W)
    assert m5.breakout_width(SH, Sh, P, Pp, T) == pytest.approx(W, abs=1e-6)


def test_mislabelled_stresses_raise():
    with pytest.raises(ValueError):
        m5.breakout_width(80 * MPa, 150 * MPa, 50 * MPa, 40 * MPa, 25.0)


# ---------------------------------------------------- mud window and verdict
def window(SH, Sh, Pp, z=5000.0, Pw=None, **kw):
    Pw = C.HYDROSTATIC_GRAD * z if Pw is None else Pw
    return se.mud_window(Sh, SH, Pp, z, Pw, 200.0, **kw)


def test_go_when_hydrostatic_mud_keeps_the_lobe_within_w_max():
    w = window(110 * MPa, 80 * MPa, 50 * MPa)
    assert w["verdict"] == "GO"
    assert w["width_hydro"] <= C.BREAKOUT_W_MAX_DEG
    assert w["overbalance_MPa"] == 0.0


def test_conditional_when_heavier_mud_inside_the_window_does_it():
    w = window(150 * MPa, 85 * MPa, 50 * MPa, UCS=120 * MPa)
    assert w["width_hydro"] > C.BREAKOUT_W_MAX_DEG
    assert w["verdict"] == "CONDITIONAL" and w["open"]
    assert w["Pw_lo"] <= w["Pw_hi"]
    assert w["overbalance_MPa"] == pytest.approx((w["Pw_lo"] - w["Pw_hydro"]) / MPa)
    assert w["overbalance_SG"] == pytest.approx(
        (w["Pw_lo"] - w["Pw_hydro"]) / (C.MUD_SG_GRAD * 5000.0))
    width_at_lo = m5.breakout_width(150 * MPa, 85 * MPa, w["Pw_lo"], 50 * MPa, 200.0,
                                    UCS=120 * MPa)
    assert width_at_lo == pytest.approx(C.BREAKOUT_W_MAX_DEG, abs=1e-6)


def test_no_go_when_the_window_is_shut():
    w = window(200 * MPa, 70 * MPa, 50 * MPa, UCS=60 * MPa)
    assert w["verdict"] == "NO-GO" and not w["open"]
    assert w["Pw_lo"] > w["Pw_hi"]
    assert w["width_hi"] > C.BREAKOUT_W_MAX_DEG
    assert w["overbalance_MPa"] is None


def test_upper_bound_is_shmin_less_the_margin():
    z = 5000.0
    w = window(150 * MPa, 85 * MPa, 50 * MPa, z=z)
    assert w["Pw_hi"] == pytest.approx(85 * MPa - C.MUD_MARGIN_SG * C.MUD_SG_GRAD * z)


def test_lower_bound_never_below_pore_pressure():
    w = window(60 * MPa, 59 * MPa, 50 * MPa, UCS=400 * MPa)
    assert w["Pw_lo"] == 50 * MPa


def test_tensile_initiation_is_reported_not_a_bound():
    """With zero tensile strength, tensile initiation at the SHmax azimuth once
    Pw > 3 Shmin - SHmax - Pp; the window and verdict ignore it."""
    SH, Sh, Pp = 150 * MPa, 85 * MPa, 50 * MPa
    w = window(SH, Sh, Pp, UCS=120 * MPa)
    assert w["Pw_tensile"] == pytest.approx(3 * Sh - SH - Pp)
    assert w["tensile_at_mud"] is True
    assert w["verdict"] == "CONDITIONAL"


def test_published_profile_initiates_tensile_fractures_at_every_mud_weight_at_12p5_km():
    """With zero tensile strength the published profile initiates tensile fractures at any
    mud above pore pressure at 12.5 km, which is why initiation can't bound the
    window. Holds under the paper's pore pressure as well as hydrostatic."""
    ref, _ = se.stability_inputs(CORNWALL, 12500.0, 400.0, UCS=CORNWALL.UCS)
    assert ref["window"]["Pw_tensile"] <= ref["Pp"]


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_window_reported_for_every_site_and_case(site):
    z = site.target_depth
    ref, cases = se.stability_inputs(site, z, float(site.geotherm(z)), UCS=site.UCS)
    for c in cases:
        w = c["window"]
        assert w["verdict"] in se.VERDICT_ORDER
        for k in ("Pw_lo", "Pw_hi", "SG_lo", "SG_hi", "width_hydro", "width_hi",
                  "Pw_tensile"):
            assert np.isfinite(w[k])


def test_v10_adapter_has_no_window():
    z = CORNWALL.target_depth
    ref, _ = se.stability_inputs(CORNWALL, z, 400.0, v10=True)
    assert ref["window"] is None


@pytest.mark.parametrize("site, verdict", [
    ("Upper Rhine Graben / Soultz-sous-Forets (France)", "CONDITIONAL"),
    ("Larderello (Italy)", "GO"),
    ("United Downs / Carnmenellis (UK)", "CONDITIONAL"),
    ("Pannonian Basin (Hungary)", "CONDITIONAL"),
])
def test_width_verdicts_without_wall_cooling(site, verdict):
    """The width verdict on its own: v1.0's depths, strengths and 200 C wall,
    and no thermal hoop stress. Judging on width rather than full
    suppression already brings United Downs back from v1.0's NO-GO."""
    s = next(x for x in SITES if x.name == site)
    z = s.target_depth
    ref, _ = se.stability_inputs(s, z, float(s.geotherm(z)), UCS=s.UCS, thermal=False)
    assert ref["window"]["verdict"] == verdict


# ------------------------------------------------------------- classify()
def _o(verdicts, ref_index=0, admissible=True, T_rock=400.0):
    cases = [dict(window=dict(verdict=v, width_hydro=50.0, overbalance_MPa=5.0,
                              SG_lo=1.1)) for v in verdicts]
    ref = dict(cases[ref_index], admissible=dict(admissible=admissible))
    return dict(m1=dict(T_bottom_delivered=150.0), stress=ref, stress_cases=cases,
                T_rock=T_rock, v10=False)


def test_classify_single_case():
    from comparative_sites import classify
    assert classify(_o(["GO"]))[2] == "GO"
    assert classify(_o(["CONDITIONAL"]))[1] == "+5.0MPa (1.10 SG)"
    assert classify(_o(["NO-GO"]))[2] == "NO-GO"


def test_classify_reports_the_range_across_cases():
    from comparative_sites import classify
    assert classify(_o(["GO", "CONDITIONAL", "NO-GO"], ref_index=2))[2] == "NO-GO [GO..NO-GO]"


def test_classify_inadmissible_stress_is_never_go():
    from comparative_sites import classify
    v = classify(_o(["GO"], admissible=False))[2]
    assert v == "CONDITIONAL (stress inadmissible)"
    assert classify(_o(["NO-GO"], admissible=False))[2] == "NO-GO (stress inadmissible)"
