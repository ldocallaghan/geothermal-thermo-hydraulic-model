"""The radial conduction solver for the wall, against analytic solutions."""
import numpy as np
import pytest

import wall_thermal as wt

A, K, RC, T_ROCK = 0.11, 2.5, 2.5e6, 200.0
ALPHA = K / RC


def test_talbot_inversion_matches_jaeger_short_time_series():
    """For small tau the constant-temperature flux is 1/sqrt(pi tau) + 1/2 -
    sqrt(tau/pi)/4 + tau/8 (Carslaw and Jaeger 1959)."""
    tau = 0.01
    series = 1 / np.sqrt(np.pi * tau) + 0.5 - 0.25 * np.sqrt(tau / np.pi) + tau / 8
    assert wt.jaeger_flux(tau)[0] == pytest.approx(series, rel=1e-4)


def test_wall_flux_matches_jaeger_from_10_minutes_to_30_days():
    res = wt.run([wt.State("fixed", 30 * 86400, T_wall=100.0)], T_ROCK, A, K, RC)
    for t in (600, 3600, 86400, 10 * 86400, 30 * 86400):
        i = np.argmin(abs(res.t - t))
        q = -res.q_wall[i] / (2 * np.pi * A)
        exact = K * (T_ROCK - 100.0) / A * wt.jaeger_flux(ALPHA * res.t[i] / A ** 2)[0]
        assert q == pytest.approx(exact, rel=0.02), t


def test_recovery_after_a_constant_rate_period():
    """With no fluid node the wall recovers as the cylinder solution, by
    superposition, to within 2% once the pause exceeds 5 a^2/alpha; the line
    source is 6 to 13% away over the same times and closes on it slowly."""
    q, t_on = -1500.0, 2 * 86400
    res = wt.run([wt.State("flux", t_on, q=q), wt.State("static", 3 * 86400, C=0.0, h=1.0)],
                 T_ROCK, A, K, RC)
    t5 = 5 * A ** 2 / ALPHA
    gaps = []
    for dt in (t5, 2 * t5, 3 * 86400):
        i = np.argmin(abs(res.t - (t_on + dt)))
        model = res.T_wall[i] - T_ROCK
        exact = wt.cylinder_recovery(q, K, ALPHA, A, t_on, res.t[i] - t_on)[0]
        assert model == pytest.approx(exact, rel=0.02), dt
        gaps.append(abs(model / wt.line_source_recovery(q, K, ALPHA, t_on, res.t[i] - t_on) - 1))
    assert gaps == sorted(gaps, reverse=True)


def _cycle(C):
    return [wt.State("circulating", 3 * 86400, T_fluid=80.0, h=1500.0),
            wt.State("static", 36 * 3600, h=28.0, C=C)]


def test_grid_refinement_changes_the_wall_by_less_than_half_a_kelvin():
    C = 4e6 * np.pi * (0.11 ** 2 - 0.07 ** 2)
    walls = [wt.run(_cycle(C), T_ROCK, A, K, RC, n=n).fields[-1][0] for n in (100, 200, 400)]
    assert max(walls) - min(walls) < 0.5


def test_energy_is_conserved_over_a_cycle():
    C = 4e6 * np.pi * (0.11 ** 2 - 0.07 ** 2)
    h = _cycle(C) + [wt.State("circulating", 86400, T_fluid=80.0, h=1500.0)]
    res = wt.run(h, T_ROCK, A, K, RC)
    assert res.energy_in + res.energy_far == pytest.approx(res.stored, rel=0.01)


def test_fluid_node_with_no_capacity_is_the_no_node_case():
    a = wt.run(_cycle(0.0), T_ROCK, A, K, RC)
    b = wt.run(_cycle(1e-6), T_ROCK, A, K, RC)
    assert b.fields[-1][0] == pytest.approx(a.fields[-1][0], abs=0.05)


def test_fluid_node_with_huge_capacity_holds_the_wall():
    res = wt.run([wt.State("static", 36 * 3600, h=1e5, C=1e15, T_node0=80.0)],
                 T_ROCK, A, K, RC)
    assert res.T_wall[-1] == pytest.approx(80.0, abs=0.5)


def test_a_fluid_node_slows_the_reheating():
    C = 4e6 * np.pi * (0.11 ** 2 - 0.07 ** 2)
    assert (wt.run(_cycle(C), T_ROCK, A, K, RC).fields[-1][0]
            < wt.run(_cycle(0.0), T_ROCK, A, K, RC).fields[-1][0])
