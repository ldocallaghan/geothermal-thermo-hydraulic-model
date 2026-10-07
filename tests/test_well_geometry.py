"""Well and string geometry, and the pipe wall as layers in series."""
import numpy as np
import pytest

import well_geometry as wg


def test_vertical_survey_gives_tvd_equal_to_md():
    s = np.array([0.0, 500.0, 3000.0, 12000.0])
    assert np.allclose(wg.Survey.vertical().tvd(s), s)
    st = wg.Survey.from_stations([0, 1000, 5000], [0, 0, 0], [0, 45, 90])
    assert np.allclose(st.tvd(s), s)


def test_straight_inclined_hole():
    """Constant inclination: TVD = MD cos(inc), and minimum curvature adds nothing."""
    st = wg.Survey.from_stations([0, 1000, 2000], [30, 30, 30], [10, 10, 10])
    assert st.tvd(1500.0) == pytest.approx(1500 * np.cos(np.radians(30)))


def test_minimum_curvature_build_section():
    """A 0-to-90 degree build over an arc of length L has radius 2L/pi, so the
    TVD gained is that radius."""
    L = 1000.0
    st = wg.Survey.from_stations([0, L], [0, 90])
    assert st.tvd(L) == pytest.approx(2 * L / np.pi)


def test_series_resistance_against_hand_calculation():
    p = wg.INTERNALLY_COATED
    hand = (np.log(60.7 / 59.7) / (2 * np.pi * 0.47)
            + np.log(69.85 / 60.7) / (2 * np.pi * 45.0))
    assert p.body_resistance() == pytest.approx(hand)
    assert p.r_bore == pytest.approx(0.0597)
    assert p.r_out == pytest.approx(0.06985)


def test_insulation_ranks_as_published():
    R = {name: p.body_resistance() for name, p in wg.PIPE_TYPES.items()}
    assert (R["conventional"] < R["externally coated"] < R["internally coated"]
            < R["dual-wall"] < R["vacuum-insulated (sensitivity)"])


def test_joint_is_bare_steel_across_the_pipe():
    p = wg.DUAL_WALL
    assert p.joint_resistance() == pytest.approx(
        np.log(p.r_out / p.r_bore) / (2 * np.pi * wg.K_STEEL))


def test_layers_must_be_concentric():
    with pytest.raises(ValueError):
        wg.PipeType("gap", (wg.Layer(0.05, 0.06, 1.0), wg.Layer(0.061, 0.07, 1.0)))


def test_annulus_changes_at_the_casing_shoe():
    g = wg.WellGeometry(
        [wg.StringSegment(3000.0, wg.CONVENTIONAL)],
        [wg.cased(0, 1000, 0.22, 0.24, 0.31), wg.open_hole(1000, 3000, 0.216)])
    a = g.at(np.array([500.0, 1500.0]))
    r_out = wg.CONVENTIONAL.r_out
    assert a["A_ann"] == pytest.approx(np.pi * (np.array([0.11, 0.108]) ** 2 - r_out ** 2))
    assert a["r_rock"] == pytest.approx([0.155, 0.108])
    cement = np.log(0.155 / 0.12) / (2 * np.pi * wg.K_CEMENT)
    steel = np.log(0.12 / 0.11) / (2 * np.pi * wg.K_STEEL)
    assert a["R_hole"] == pytest.approx([steel + cement, 0.0])


def test_string_segments_map_to_measured_depth():
    g = wg.WellGeometry(
        [wg.StringSegment(1000.0, wg.DUAL_WALL), wg.StringSegment(500.0, wg.CONVENTIONAL)],
        [wg.open_hole(0, 1500, 0.216)])
    assert g.md_bit == 1500.0
    r = g.at(np.array([0.0, 999.0, 1001.0, 1500.0]))["r_bore"]
    assert r == pytest.approx([0.040, 0.040, 0.0607, 0.0607])


def test_string_must_fit_and_hole_must_cover_it():
    with pytest.raises(ValueError):
        wg.WellGeometry([wg.StringSegment(100.0, wg.CONVENTIONAL)],
                        [wg.open_hole(0, 100, 0.12)])
    with pytest.raises(ValueError):
        wg.WellGeometry([wg.StringSegment(100.0, wg.CONVENTIONAL)],
                        [wg.open_hole(0, 50, 0.216)])


@pytest.mark.slow
def test_n_identical_segments_equal_one():
    import model1_coupled as m1
    one = m1.solve(m_dot=4.0, pipe=wg.LEGACY_VACUUM, verbose=False)
    L = one["L"]
    g = wg.WellGeometry([wg.StringSegment(L / 4, wg.LEGACY_VACUUM)] * 4,
                        [wg.HoleInterval(0.0, L, wg.C.R_WELL)])
    four = m1.solve(m_dot=4.0, geometry=g, verbose=False)
    assert four["T_bottom_delivered"] == pytest.approx(one["T_bottom_delivered"], abs=1e-3)
    assert four["T_return_surface"] == pytest.approx(one["T_return_surface"], abs=1e-3)


@pytest.mark.slow
@pytest.mark.parametrize("pipe", [wg.INTERNALLY_COATED, wg.DUAL_WALL])
def test_joint_fraction_limits(pipe):
    """Fraction 0 is the uniform pipe body; fraction 1 is bare steel throughout."""
    import model1_coupled as m1
    z = np.linspace(0.0, 10000.0, 7)
    T_d, T_u = np.full_like(z, 60.0), np.full_like(z, 150.0)
    bare = wg.PipeType("bare", (wg.Layer(pipe.r_bore, pipe.r_out, wg.K_STEEL),))

    def UAi(p):
        g = wg.WellGeometry.single(10000.0, p, r_hole=0.108)
        return m1.conductances(z, T_d, T_u, 10.0, 1.0, 2.5, g)[0]

    assert UAi(pipe.with_joints(0.0)) == pytest.approx(UAi(pipe))
    assert UAi(pipe.with_joints(1.0)) == pytest.approx(UAi(bare))
    half = UAi(pipe.with_joints(0.5))
    assert half == pytest.approx(0.5 * (UAi(pipe) + UAi(bare)))
