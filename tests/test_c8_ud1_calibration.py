"""C8: UD-1 strength calibration (spec v1.1), and spec test 5.

Test 5: predicted breakouts fall inside 900-4,000 m at section mud weights,
tensile fractures are predicted where they were logged, and the hole stays
drillable at 5,058 m.
"""
import numpy as np
import pytest

import calibrate_ud1 as cu
import geo_constants as C
import model5_convergence_confinement as m5
from comparative_sites import CORNWALL

MPa = 1.0e6


@pytest.fixture(scope="module")
def feats():
    return cu.load_log()


@pytest.fixture(scope="module")
def cal(feats):
    return cu.calibrate(feats)


# ------------------------------------------------------------------ inputs
def test_log_contents(feats):
    classes = [f["cls"] for f in feats]
    assert len(feats) == 12031
    assert classes.count("Breakout") == 75
    assert classes.count("Induced fracture") == 2


def test_section_breakouts(cal):
    bo = cal["breakouts"]
    assert len(bo) == 24
    assert sum(b["height"] for b in bo) == pytest.approx(46.0, abs=0.1)
    assert min(b["width"] for b in bo) == pytest.approx(30.0, abs=0.5)
    assert max(b["width"] for b in bo) == pytest.approx(63.3, abs=0.1)


def test_md_to_tvd(cal, feats):
    tvd = cal["tvd"]
    assert tvd(907.15) == pytest.approx(907.15)
    assert tvd(4000.0) == pytest.approx(3990.0, abs=1.0)    # <= 10 m correction
    # independent check: extend the last station to TD (5,275 m MD) at its
    # inclination and land on the paper's 5,058 m TVD
    last = feats[-1]
    td = tvd(last["md"]) + (5275.0 - last["md"]) * np.cos(np.radians(last["dev"]))
    assert td == pytest.approx(5058.0, abs=15.0)


def test_formation_temperature_is_the_papers():
    assert cu.T_formation(5000.0) == pytest.approx(180.0)
    assert cu.T_formation(0.0) == C.SURFACE_TEMP


def test_wall_never_colder_than_the_surface():
    s = cu.state(1000.0, 40.0)
    assert s["T_wall"] == C.SURFACE_TEMP
    assert cu.state(4000.0, 40.0)["T_wall"] == pytest.approx(cu.T_formation(4000.0) - 40.0)


# ---------------------------------------------------------- inversion
@pytest.mark.parametrize("W", [0.0, 30.0, 52.0, 63.0, 89.0])
@pytest.mark.parametrize("cooling", [0.0, 40.0])
def test_inversion_round_trips_through_model5(W, cooling):
    """implied_ucs inverts model5.breakout_width exactly at 35 deg friction."""
    z = 3800.0
    u = cu.implied_ucs(z, W, cooling)
    assert cu.predicted_width(z, u, cooling) == pytest.approx(W, abs=1e-6)


def test_friction_angle_35_is_model5():
    assert cu.kmc(35.0) == pytest.approx(m5.KMC)


# --------------------------------------------------------- calibration
def test_weak_zones_are_weaker_than_intact_rock(cal):
    for lv in cal["levels"].values():
        assert lv["weak_p90"] < lv["intact_min"]
        assert lv["weak_p10"] < lv["weak_p50"] < lv["weak_p90"]


def test_more_cooling_means_weaker_rock(cal):
    """Cooling relieves the hoop stress, so the same breakout needs weaker
    rock: about 0.7 MPa of UCS per kelvin."""
    lv = cal["levels"]
    p50 = [lv[k]["weak_p50"] for k in sorted(lv)]
    assert all(a > b for a, b in zip(p50, p50[1:]))
    slope = (lv[0.0]["weak_p50"] - lv[40.0]["weak_p50"]) / 40.0 / MPa
    assert slope == pytest.approx(0.72, abs=0.05)


def test_intact_bound_binds_at_the_foot_of_the_section(cal):
    for lv in cal["levels"].values():
        assert lv["intact_binding_tvd"] > 3950.0


def test_site_strength_cases_are_the_calibration(cal):
    got = dict(CORNWALL.strength_cases)
    for (label, u), want in zip(cu.strength_cases(cal),
                                ("intact, 0 K", "intact, 40 K",
                                 "weak zones, 0 K", "weak zones, 40 K")):
        assert got[want] == pytest.approx(u, abs=0.5 * MPa), label


def test_site_calibrated_w_max_is_the_widest_12p25_breakout(cal):
    assert C.BREAKOUT_W_MAX_SITE_DEG == round(max(b["width"] for b in cal["breakouts"]))


# ------------------------------------------------------------- test 5
@pytest.mark.parametrize("cooling", [0.0, 40.0])
def test_predicted_breakouts_fall_inside_the_logged_interval(cal, cooling):
    """At 1.05 SG, weak-zone rock starts breaking out at ~2.9 km TVD, inside
    900-4,000 m (logged: one at 2,210 m, the rest from 3,095 m MD), and
    intact rock not until the foot of the section."""
    lv = cal["levels"][cooling]
    z_weak = cu.onset_depth(lv["weak_p50"], cooling)
    z_intact = cu.onset_depth(lv["intact_min"], cooling)
    assert 900.0 < z_weak < 4000.0
    assert z_weak == pytest.approx(2900.0, abs=50.0)
    assert z_intact == pytest.approx(3987.0, abs=5.0)


def test_predicted_widths_match_the_lower_section(cal, feats):
    rows = cu.depth_table(cal, feats, 20.0)
    for r in rows:
        if r["md"][0] >= 3650.0:
            assert r["w_weak"] == pytest.approx(r["width_obs"], abs=8.0)


@pytest.mark.parametrize("cooling", [0.0, 40.0])
@pytest.mark.parametrize("z", [2500.0, 2665.0, 3667.0, 3700.0])
def test_tensile_fractures_predicted_where_logged(cooling, z):
    """At 1.05 SG, tensile initiation (T0 = 0) is predicted at the logged
    depths. It is also predicted above and below them: with T0 = 0 the model
    does not discriminate where tensile fractures form."""
    assert cu.tensile_threshold_SG(z, cooling) < cu.MUD_SG


@pytest.mark.parametrize("label, ucs", CORNWALL.strength_cases)
@pytest.mark.parametrize("cooling", [0.0, 40.0])
def test_hole_drillable_at_td(label, ucs, cooling):
    w = cu.td_window(ucs, cooling)
    assert w["open"] and w["verdict"] in ("GO", "CONDITIONAL")


@pytest.mark.parametrize("z", np.arange(4000.0, 5058.0, 100.0))
def test_deviated_section_not_made_undrillable(cal, z):
    """Qualitative check on the 8.5" section: at the weakest calibrated case
    the window stays open (vertical-hole approximation)."""
    ucs = cal["levels"][40.0]["weak_p50"]
    s = cu.state(z, 40.0)
    from site_evaluation import mud_window
    w = mud_window(s["Sh"], s["SH"], s["Pp"], z, s["Pw"], s["T_rock"],
                   T_wall=s["T_wall"], UCS=ucs, dsigma_T=s["dT"])
    assert w["open"]


def test_figure_is_written(cal, feats, tmp_path):
    p = cu.figure(cal, feats, path=str(tmp_path / "ud1.png"))
    assert (tmp_path / "ud1.png").stat().st_size > 10_000
