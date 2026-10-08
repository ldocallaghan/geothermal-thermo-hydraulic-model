"""
wall_thermal.py
===============
Temperature of the rock around a borehole through the drilling cycle.

One-dimensional axisymmetric conduction in the rock,

    dT/dt = alpha (1/r) d/dr (r dT/dr),        a <= r <= r_max,

solved by implicit (backward Euler) finite volumes on a log-spaced grid, with
the far field held at the rock temperature. The boundary at the wall depends
on the state of the well:

  * circulating: convection to the annulus fluid, T_fluid and h given (from
    Model 1 for the site's pipe and flow);
  * static: the fluid in the hole is one well-mixed node of heat capacity C
    [J/m/K], coupled to the wall by a film coefficient h (the laminar floor,
    Nu = 4.36), starting from its temperature at the end of circulation; with
    C = 0 the wall is insulated (the no-node sensitivity);
  * fixed wall temperature, and fixed heat flux: for the analytic checks.

A history is a list of State; run() returns the wall temperature through time
and the field at the end of each state.

Checks (tests/test_wall_thermal.py): the wall flux against the constant-
temperature cylinder (Jaeger 1956; Carslaw and Jaeger 1959), the late
recovery after a constant-rate period against the line source (Bullard 1947),
grid convergence and energy conservation.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded
from scipy.special import kv

NU_LAMINAR = 4.36


@dataclass
class State:
    """One state of the well. kind: "circulating" (T_fluid, h), "static"
    (h, C; the fluid node starts at T_node0, or at the last circulating fluid
    temperature if None), "fixed" (T_wall) or "flux" (q, W/m into the wall,
    negative when heat is drawn out)."""
    kind: str
    duration: float                 # s
    T_fluid: float = None
    h: float = None
    C: float = 0.0
    T_node0: float = None
    T_wall: float = None
    q: float = None
    label: str = ""


@dataclass
class Grid:
    a: float
    r_max: float
    n: int = 200

    def build(self):
        """Node radii (log-spaced, the first at the wall), face radii, and the
        volume per unit length of each control volume."""
        r = self.a * (self.r_max / self.a) ** np.linspace(0.0, 1.0, self.n)
        faces = np.concatenate([[r[0]], np.sqrt(r[:-1] * r[1:]), [r[-1]]])
        vol = np.pi * (faces[1:] ** 2 - faces[:-1] ** 2)
        return r, faces, vol


def r_max_for(a, alpha, t_total):
    """a + 10 sqrt(alpha t), and at least 30 m."""
    return max(a + 10.0 * np.sqrt(alpha * t_total), 30.0)


def time_steps(duration, dt0=1.0, growth=1.08, dt_max=None):
    """Steps that start small after a change of state and grow geometrically."""
    dt_max = dt_max or duration / 50.0
    out, t, dt = [], 0.0, dt0
    while t < duration - 1e-9:
        dt = min(dt, dt_max, duration - t)
        out.append(dt)
        t += dt
        dt *= growth
    return out


@dataclass
class Result:
    t: np.ndarray
    T_wall: np.ndarray
    T_node: np.ndarray
    q_wall: np.ndarray              # W/m into the rock at the wall
    r: np.ndarray
    fields: list = field(default_factory=list)   # T(r) at the end of each state
    energy_in: float = 0.0          # J/m through the wall, total
    energy_far: float = 0.0         # J/m through the far boundary, total
    stored: float = 0.0             # J/m change in the rock's stored heat
    vol: np.ndarray = None
    ends: list = field(default_factory=list)     # (t, T_wall, T_node) at the end of each state


def run(history, T_rock, a, k, rhocp, n=200, r_max=None, T0=None, dt0=1.0, growth=1.08):
    """Conduct heat through the history of states. T_rock: far-field and
    initial temperature (or T0, an initial field on the grid)."""
    alpha = k / rhocp
    total = sum(s.duration for s in history)
    r, faces, vol = Grid(a, r_max or r_max_for(a, alpha, total), n).build()
    T = np.full(n, float(T_rock)) if T0 is None else np.array(T0, float)
    T_start = T.copy()
    # conductance between neighbouring nodes per unit length
    G = 2 * np.pi * k / np.log(r[1:] / r[:-1])
    cap = rhocp * vol
    ts, Tw, Tn, qw = [0.0], [T[0]], [np.nan], [0.0]
    fields, ends, t = [], [], 0.0
    e_in = e_far = 0.0
    T_fluid_last = None
    for s in history:
        T_node = None
        if s.kind == "static" and s.C > 0:
            T_node = s.T_node0 if s.T_node0 is not None else (
                T_fluid_last if T_fluid_last is not None else T[0])
        for dt in time_steps(s.duration, dt0, growth):
            m = n + (1 if T_node is not None else 0)
            ab = np.zeros((3, m))
            rhs = np.zeros(m)
            # rock nodes: cap/dt (T - T_old) = sum G (T_j - T_i)
            diag = cap / dt
            ab[1, :n] = diag
            rhs[:n] = diag * T
            ab[1, :n - 1] += G
            ab[1, 1:n] += G
            ab[0, 1:n] = -G                  # upper diagonal: A[i, i+1]
            ab[2, :n - 1] = -G               # lower diagonal: A[i+1, i]
            # far field fixed
            ab[1, n - 1], ab[2, n - 2], rhs[n - 1] = 1.0, 0.0, T_rock
            # the wall
            A_w = 2 * np.pi * a
            if s.kind == "circulating":
                hA = s.h * A_w
                ab[1, 0] += hA
                rhs[0] += hA * s.T_fluid
            elif s.kind == "flux":
                rhs[0] += s.q
            elif s.kind == "fixed":
                ab[1, 0], ab[0, 1], rhs[0] = 1.0, 0.0, s.T_wall
            elif s.kind == "static" and T_node is not None:
                hA = s.h * A_w
                ab[1, 0] += hA
                ab[0, n] = 0.0
                # couple the wall (0) and the node (n): only adjacent in a full
                # matrix, so solve this case densely below
            T_old_wall = T[0]
            if s.kind == "static" and T_node is not None:
                A = np.zeros((m, m))
                for i in range(n):
                    A[i, i] = ab[1, i]
                    if i + 1 < n:
                        A[i, i + 1] = ab[0, i + 1]
                        A[i + 1, i] = ab[2, i]
                A[0, n] = -hA
                A[n, 0] = -hA
                A[n, n] = s.C / dt + hA
                rhs[n] = s.C / dt * T_node
                sol = _solve_bordered(A, rhs, n)
                T_new, T_node = sol[:n], sol[n]
            else:
                T_new = solve_banded((1, 1), ab[:, :n], rhs[:n])
            # heat into the rock at the wall over the step, from the wall
            # node's balance: cap0 dT0/dt = q_in + G0 (T1 - T0)
            q_in = cap[0] * (T_new[0] - T[0]) / dt - G[0] * (T_new[1] - T_new[0])
            q_far = G[-1] * (T_new[-1] - T_new[-2])     # into the far node
            e_in += q_in * dt
            e_far += q_far * dt
            T = T_new
            t += dt
            ts.append(t); Tw.append(T[0]); qw.append(q_in)
            Tn.append(T_node if T_node is not None else np.nan)
        if s.kind == "circulating":
            T_fluid_last = s.T_fluid
        fields.append(T.copy())
        ends.append((t, float(T[0]), float(T_node) if T_node is not None else np.nan))
    stored = float(np.sum(cap[:-1] * (T[:-1] - T_start[:-1])))
    return Result(np.array(ts), np.array(Tw), np.array(Tn), np.array(qw), r, fields,
                  e_in, e_far, stored, vol, ends)


def _solve_bordered(A, rhs, n):
    """Tridiagonal rock block bordered by the fluid node; a dense solve is
    small enough (n ~ 200) to need nothing cleverer."""
    return np.linalg.solve(A, rhs)


# ------------------------------------------------------------- analytic checks
def _talbot(F, t, M=24):
    """Inverse Laplace transform of F at t, fixed Talbot method (Abate and
    Valko 2004)."""
    rr = 2.0 * M / (5.0 * t)
    acc = 0.5 * np.real(F(rr + 0j)) * np.exp(rr * t)
    for k_ in range(1, M):
        th = k_ * np.pi / M
        cot = np.cos(th) / np.sin(th)
        s = rr * th * (cot + 1j)
        sig = th + (th * cot - 1.0) * cot
        acc += np.real(np.exp(t * s) * F(s) * (1.0 + 1j * sig))
    return rr / M * acc


def cylinder_flux_temperature(tau):
    """Dimensionless wall temperature change 2 pi k dT / q of a cylinder
    carrying a constant heat flow q from tau = 0, at tau = alpha t / a^2: the
    inverse transform of K0(sqrt p) / (p^1.5 K1(sqrt p)) (Carslaw and Jaeger
    1959, section 13.5)."""
    F = lambda p: kv(0, np.sqrt(p)) / (p ** 1.5 * kv(1, np.sqrt(p)))
    return np.array([_talbot(F, t) for t in np.atleast_1d(np.asarray(tau, float))])


def cylinder_recovery(q, k, alpha, a, t_on, dt_off):
    """Wall temperature change dt_off after a period t_on of constant heat
    flow q [W/m into the rock] ends, for the cylinder of radius a, by
    superposition. The line source (below) is its limit when both times are
    long against a^2/alpha."""
    f = cylinder_flux_temperature
    tau_on, tau = alpha * t_on / a ** 2, alpha * np.asarray(dt_off, float) / a ** 2
    return q / (2 * np.pi * k) * (f(tau_on + tau) - f(tau))


def jaeger_flux(tau, n_terms=24):
    """Dimensionless wall flux q a / (k (T_rock - T_wall)) for a cylinder held
    at a constant temperature in an infinite medium, at tau = alpha t / a^2:
    the inverse Laplace transform of K1(sqrt p) / (sqrt p K0(sqrt p)) (Carslaw
    and Jaeger 1959, section 13.5), by the fixed Talbot method (Abate and
    Valko 2004)."""
    tau = np.atleast_1d(np.asarray(tau, float))
    M = n_terms
    out = np.empty_like(tau)
    for j, t in enumerate(tau):
        rr = 2.0 * M / (5.0 * t)
        F = lambda p: kv(1, np.sqrt(p)) / (np.sqrt(p) * kv(0, np.sqrt(p)))
        acc = 0.5 * np.real(F(rr + 0j)) * np.exp(rr * t)
        for k_ in range(1, M):
            th = k_ * np.pi / M
            cot = np.cos(th) / np.sin(th)
            s = rr * th * (cot + 1j)
            sig = th + (th * cot - 1.0) * cot
            acc += np.real(np.exp(t * s) * F(s) * (1.0 + 1j * sig))
        out[j] = rr / M * acc
    return out


def line_source_recovery(q, k, alpha, t_on, dt_off):
    """Wall temperature change (from rock temperature) a time dt_off after a
    period t_on of constant heat flow q [W/m into the rock] ends, by
    superposition of line sources (Bullard 1947; the Horner form at late
    times): q/(4 pi k) [E1(r^2/4 alpha dt_off) - E1(r^2/4 alpha (t_on + dt_off))]
    -> q/(4 pi k) ln((t_on + dt_off) / dt_off)."""
    return q / (4.0 * np.pi * k) * np.log((t_on + dt_off) / dt_off)


def cooling_shape(t, a, k, rhocp, h, n=120):
    """The cooled zone after circulating for t seconds, per kelvin of fluid
    cooling: g(r) = (T_rock - T(r)) / (T_rock - T_fluid). Conduction is linear,
    so T(r) = T_rock - cooling g(r) for any rock temperature and cooling."""
    res = run([State("circulating", t, T_fluid=-1.0, h=h)], 0.0, a, k, rhocp, n=n, growth=1.2)
    return res.r, -res.fields[-1]
