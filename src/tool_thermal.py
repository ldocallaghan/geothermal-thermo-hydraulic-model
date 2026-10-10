"""
tool_thermal.py
===============
The temperature of the bottom-hole assembly and of the coated drill pipe
through the drilling cycle.

The BHA is lumped per metre into two nodes: the steel of the collars, and
the fluid inside them, where the MWD and LWD electronics sit; the tool
temperature is that fluid's. The steel exchanges heat with the fluid
outside it and with the fluid inside through film coefficients:

    C_s dT_s/dt = H_out (T_out - T_s) + H_in (T_f - T_s)
    C_f dT_f/dt = H_in (T_s - T_f)            (static inside)
    T_f = T_flow                              (fluid flowing through)

with H = h x perimeter. While circulating, the films are Model 1's (forced
convection in the bore and the annulus); with the pumps off, the laminar
floor of the static fluid node (wall_thermal.NU_LAMINAR). Outside, the
fluid is the annulus while circulating, the static fluid node of the
connection, and the column at the end of the trip on the way back in.
Integrated by backward Euler.

Assumptions: the fluid inside the string does not exchange with the column
while the string is run in (a float valve in the BHA); the BHA's own heat
capacity is not fed back into the static node, which therefore runs hotter
than it would, the conservative side for the tool; staging circulation
delivers the bore temperature Model 1 gives at that depth, and the column
below the bit is not cooled by it.

Steel: plain carbon steel at 300 K (Incropera, F.P., DeWitt, D.P., Bergman,
T.L. and Lavine, A.S. (2007), Fundamentals of Heat and Mass Transfer, 6th
ed., Table A.1). Its heat capacity rises with temperature, so the value at
300 K heats the steel fastest.
"""
from dataclasses import dataclass, field

import numpy as np

STEEL_RHO = 7854.0      # kg/m3
STEEL_CP = 434.0        # J/kg/K


@dataclass
class Lump:
    """A length of string per metre: steel heat capacity C_s and internal
    fluid heat capacity C_f [J/m/K], and the inner and outer perimeters [m]."""
    C_s: float
    C_f: float
    P_in: float
    P_out: float
    D_in: float

    @classmethod
    def of_pipe(cls, r_bore, r_out, rhocp_fluid):
        steel = STEEL_RHO * STEEL_CP * np.pi * (r_out ** 2 - r_bore ** 2)
        return cls(C_s=steel, C_f=rhocp_fluid * np.pi * r_bore ** 2, P_in=2 * np.pi * r_bore,
                   P_out=2 * np.pi * r_out, D_in=2 * r_bore)


@dataclass
class Segment:
    """A stretch of time: T_out the fluid outside (a value, or a function of
    the time [s] from the segment's start), h_out and h_in the films, and
    T_flow the fluid flowing through the string (None: static inside)."""
    duration: float
    T_out: object
    h_out: float
    h_in: float
    T_flow: float = None
    label: str = ""


@dataclass
class History:
    t: np.ndarray
    T_s: np.ndarray
    T_f: np.ndarray
    labels: list = field(default_factory=list)

    def peak(self):
        i = int(np.argmax(self.T_f))
        return float(self.T_f[i]), float(self.t[i])

    def time_above(self, T, which="T_f"):
        """Time [s] the node spends above T."""
        x = getattr(self, which)
        dt = np.diff(self.t)
        return float(np.sum(dt[(x[1:] > T) | (x[:-1] > T)]))


def run(lump, segments, T_s0, T_f0, dt=10.0):
    """Integrate the two nodes through the segments."""
    t, Ts, Tf = 0.0, float(T_s0), float(T_f0)
    ts, Tss, Tfs, labels = [t], [Ts], [Tf], [""]
    for seg in segments:
        n = max(1, int(np.ceil(seg.duration / dt)))
        h = seg.duration / n
        Ho, Hi = seg.h_out * lump.P_out, seg.h_in * lump.P_in
        for k in range(1, n + 1):
            To = seg.T_out(k * h) if callable(seg.T_out) else seg.T_out
            if seg.T_flow is not None:
                Tf = seg.T_flow
                Ts = (lump.C_s / h * Ts + Ho * To + Hi * Tf) / (lump.C_s / h + Ho + Hi)
            else:
                # [C_s/h + Ho + Hi, -Hi; -Hi, C_f/h + Hi] [Ts; Tf] = [C_s/h Ts + Ho To; C_f/h Tf]
                a, b = lump.C_s / h + Ho + Hi, -Hi
                c, d = -Hi, lump.C_f / h + Hi
                r1, r2 = lump.C_s / h * Ts + Ho * To, lump.C_f / h * Tf
                det = a * d - b * c
                Ts, Tf = (r1 * d - b * r2) / det, (a * r2 - c * r1) / det
            t += h
            ts.append(t); Tss.append(Ts); Tfs.append(Tf); labels.append(seg.label)
    return History(np.array(ts), np.array(Tss), np.array(Tfs), labels)


def steady_steel(T_flow, T_out, h_in, h_out, lump):
    """The steel's temperature with fluid flowing through at T_flow and
    outside at T_out, at steady state."""
    Hi, Ho = h_in * lump.P_in, h_out * lump.P_out
    return (Hi * T_flow + Ho * T_out) / (Hi + Ho)


def trip_in(lump, column, md_bottom, speed, T0, h_static_out, h_static_in, stages=(),
            stage_s=0.0, circulating=None, dt=30.0):
    """Run the string into the hole at speed [m/s] to md_bottom, through a
    column whose temperature is column(md), from T0 at surface, holding for
    stage_s of circulation at each of the measured depths in stages;
    circulating(md) gives (T_flow, h_in, h_out) there. Returns the History
    and the depth [m] at each of its times."""
    T_s, T_f = T0, T0
    depth, t_all, Ts_all, Tf_all, labels = [0.0], [0.0], [T0], [T0], [""]
    stops = sorted(d for d in stages if 0.0 < d < md_bottom)
    md, t = 0.0, 0.0
    for target in stops + [md_bottom]:
        while md < target - 1e-9:
            step = min(speed * dt, target - md)
            h = step / speed
            md += step
            seg = Segment(h, float(column(md)), h_static_out, h_static_in)
            hist = run(lump, [seg], T_s, T_f, dt=h)
            T_s, T_f = float(hist.T_s[-1]), float(hist.T_f[-1])
            t += h
            depth.append(md); t_all.append(t); Ts_all.append(T_s); Tf_all.append(T_f)
            labels.append("running in")
        if target < md_bottom and stage_s > 0 and circulating is not None:
            T_flow, h_in, h_out = circulating(md)
            hist = run(lump, [Segment(stage_s, T_flow, h_out, h_in, T_flow=T_flow)], T_s, T_f)
            T_s, T_f = float(hist.T_s[-1]), float(hist.T_f[-1])
            t += stage_s
            depth.append(md); t_all.append(t); Ts_all.append(T_s); Tf_all.append(T_f)
            labels.append("staging")
    return History(np.array(t_all), np.array(Ts_all), np.array(Tf_all), labels), np.array(depth)
