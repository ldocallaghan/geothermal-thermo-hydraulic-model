"""Grid and time-step convergence of the wall solver at the breakout's
reference depth over the first minutes: 20 s of cooling as the bit exposes
the rock, then a 90th-percentile connection (6.3 min) with the fluid node.

The fresh-rock check reads the end of the connection, which halving both the
grid spacing and the time steps moves by 0.7 K at the site resolution (150
nodes) and refining eight times by 0.8 K. After the 20 s of cooling alone the
same halving moves it by 1.3 K: the first cell, about 4.5 mm, is as wide as
the heat has diffused, so that instant is not resolved to 1 K. The time step
is converged (0.04 K)."""
import numpy as np
import pytest

import wall_thermal as wt

A, K, RHOCP = 0.1207, 2.9, 2.7e6          # Soultz: hole radius, conductivity, rho c
T_ROCK, T_FLUID, H = 400.0, 190.0, 2000.0  # rock, annulus fluid [C], film [W/m2/K]
C_NODE, H_STATIC = 3.0e5, 60.0             # node [J/m/K] and its laminar film
D_REF = 0.016                              # the uncooled breakout's depth [m]


def at_reference(n, dt0, growth):
    hist = [wt.State("circulating", 20.0, T_fluid=T_FLUID, h=H),
            wt.State("static", 6.3 * 60, h=H_STATIC, C=C_NODE, T_node0=T_FLUID)]
    res = wt.run(hist, T_ROCK, A, K, RHOCP, n=n, r_max=30.0, dt0=dt0, growth=growth)
    return [float(np.interp(A + D_REF, res.r, f)) for f in res.fields]


def test_the_end_of_the_connection_converges_within_a_kelvin():
    coarse = at_reference(150, 1.0, 1.2)
    fine = at_reference(300, 0.5, np.sqrt(1.2))
    assert coarse[-1] == pytest.approx(fine[-1], abs=1.0)
    assert coarse[-1] == pytest.approx(at_reference(1200, 0.125, 1.2 ** 0.125)[-1], abs=1.0)


def test_the_time_step_is_converged():
    assert at_reference(150, 1.0, 1.2)[-1] == pytest.approx(at_reference(150, 0.1, 1.02)[-1], abs=0.1)


def test_the_first_twenty_seconds_are_resolved_to_a_kelvin_and_a_half():
    coarse, fine = at_reference(150, 1.0, 1.2), at_reference(300, 0.5, np.sqrt(1.2))
    assert coarse[0] == pytest.approx(fine[0], abs=1.5)
    assert coarse[0] < T_ROCK
