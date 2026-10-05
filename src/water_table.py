"""
water_table.py
==============
Fast, VECTORIZED water properties for the borehole models.

IAPWS-95 is accurate but ~1 ms/call -- far too slow to call per node inside a
BVP solve. Here we evaluate IAPWS-95 ONCE on a (T, P) grid, cache it to disk,
and expose bilinear RegularGridInterpolators that take and return numpy arrays.
This is standard practice (build a steam-table once, interpolate thereafter) and
is smooth enough for these models. The table is loaded, or built, on first use,
so importing the models costs nothing until water properties are needed.

Properties: rho [kg/m^3], cp [J/kg/K], mu [Pa.s], k [W/m/K].
Inputs to the lookups: T [degC], P [Pa] (arrays ok).
"""
import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from iapws import IAPWS95

KELVIN = 273.15
_CACHE = os.path.join(os.path.dirname(__file__), "water_table.npz")

# Grid spans the whole borehole envelope with margin.
T_GRID = np.arange(1.0, 481.0, 2.0)            # degC
P_GRID = np.arange(0.1, 220.1, 2.0) * 1e6      # Pa


def _build():
    nT, nP = len(T_GRID), len(P_GRID)
    rho = np.empty((nT, nP)); cp = np.empty((nT, nP))
    mu = np.empty((nT, nP)); k = np.empty((nT, nP))
    for i, Tc in enumerate(T_GRID):
        for j, Pp in enumerate(P_GRID):
            w = IAPWS95(T=Tc + KELVIN, P=Pp / 1e6)
            if w.cp is None or w.mu is None or w.k is None:
                w = IAPWS95(T=Tc + KELVIN, P=Pp / 1e6 * 1.001 + 1e-3)
            rho[i, j] = w.rho
            cp[i, j] = w.cp * 1000.0
            mu[i, j] = w.mu
            k[i, j] = w.k
    np.savez(_CACHE, T=T_GRID, P=P_GRID, rho=rho, cp=cp, mu=mu, k=k)
    return rho, cp, mu, k


def _load():
    if os.path.exists(_CACHE):
        d = np.load(_CACHE)
        if (len(d["T"]) == len(T_GRID) and len(d["P"]) == len(P_GRID)):
            return d["rho"], d["cp"], d["mu"], d["k"]
    return _build()


_INTERP = {}


def _interpolators():
    """RHO, CP, MU, K interpolators, loading or building the table once."""
    if not _INTERP:
        kw = dict(bounds_error=False, fill_value=None)  # None -> extrapolate at edges
        for name, arr in zip(("RHO", "CP", "MU", "K"), _load()):
            _INTERP[name] = RegularGridInterpolator((T_GRID, P_GRID), arr, **kw)
    return _INTERP


def __getattr__(name):
    """Lazy module attributes: wt.RHO, wt.CP, wt.MU, wt.K."""
    if name in ("RHO", "CP", "MU", "K"):
        return _interpolators()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def _clip(T_C, P_Pa):
    T_C = np.clip(T_C, T_GRID[0], T_GRID[-1])
    P_Pa = np.clip(P_Pa, P_GRID[0], P_GRID[-1])
    return np.column_stack([np.atleast_1d(T_C), np.atleast_1d(P_Pa)])


def props(T_C, P_Pa):
    """Vectorized: returns (rho, cp, mu, k) arrays for T[degC], P[Pa]."""
    pts = _clip(T_C, P_Pa)
    f = _interpolators()
    return f["RHO"](pts), f["CP"](pts), f["MU"](pts), f["K"](pts)


def cp(T_C, P_Pa):
    return _interpolators()["CP"](_clip(T_C, P_Pa))


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    print("table shape:", _interpolators()["CP"].values.shape, "cached at", _CACHE)
    for T in (15, 200, 450):
        for P in (5e6, 120e6):
            r, c, m, kk = props(np.array([T]), np.array([P]))
            print(f"T={T:4} P={P/1e6:5.0f}MPa rho={r[0]:7.1f} cp={c[0]:7.0f} mu={m[0]*1e6:6.1f} k={kk[0]:.3f}")
