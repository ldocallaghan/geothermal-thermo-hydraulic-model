"""
well_geometry.py
================
The well and drill string that Model 1 circulates through.

Positions are measured depth s along the hole, from surface (s = 0) to the bit.
A survey maps s to true vertical depth, where rock temperature and hydrostatic
pressure are taken.

  * The STRING is a list of segments from surface to bit (drill pipe, insulated
    pipe, heavy-weight pipe, collars). Each has a pipe type: the wall from bore
    to outside as concentric layers, and the fraction of its length taken up by
    bare tool joints.
  * The HOLE is a list of intervals. Each gives the radius of the annulus's
    outer boundary (casing ID, or the open hole) and the layers between it and
    the rock (casing steel, cement). In open hole there are none.

Resistances are per unit length [K m/W]: a layer from r_in to r_out with
conductivity k is ln(r_out/r_in) / (2 pi k), and layers add in series.
"""
from dataclasses import dataclass

import numpy as np

import geo_constants as C

K_STEEL = 45.0   # W/m/K, drill pipe and casing (Wu et al. 2025, Tables 3 to 5)
K_CEMENT = 1.0   # W/m/K, cement sheath (Dufrene et al. 2025)
INCH = 0.0254


@dataclass(frozen=True)
class Layer:
    r_in: float
    r_out: float
    k: float

    def resistance(self):
        return np.log(self.r_out / self.r_in) / (2 * np.pi * self.k)


def _check_concentric(layers):
    for a, b in zip(layers, layers[1:]):
        if not np.isclose(a.r_out, b.r_in):
            raise ValueError(f"layers not concentric: {a} then {b}")


@dataclass(frozen=True)
class PipeType:
    """A drill pipe's wall, inside out. A tool joint is modelled as bare steel
    across the same bore and outside radius."""
    name: str
    layers: tuple
    joint_fraction: float = 0.0
    source: str = ""

    def __post_init__(self):
        _check_concentric(self.layers)
        if not 0.0 <= self.joint_fraction <= 1.0:
            raise ValueError("joint_fraction must be in [0, 1]")

    @property
    def r_bore(self):
        return self.layers[0].r_in

    @property
    def r_out(self):
        return self.layers[-1].r_out

    def body_resistance(self):
        return sum(l.resistance() for l in self.layers)

    def joint_resistance(self):
        return Layer(self.r_bore, self.r_out, K_STEEL).resistance()

    def with_joints(self, fraction):
        return PipeType(self.name, self.layers, fraction, self.source)


def uniform_wall(k, r_in=C.R_INNER_PIPE_IN, r_out=C.R_INNER_PIPE_OUT, name=None):
    """One wall layer of conductivity k, on the model's original 100 mm bore
    and 120 mm outside diameter unless given."""
    return PipeType(name or f"uniform wall, {k} W/m/K", (Layer(r_in, r_out, k),))


# Named pipe types. Dimensions are Wu et al. (2025), Table 5: 5-1/2 inch pipe,
# 139.7 mm OD and 121.4 mm ID. Their conductivities are as cited there.
# Joint geometry is not published for any of them, so joint_fraction is 0 and
# joints are a sensitivity (with_joints).
_R_ID, _R_OD = 0.1214 / 2, 0.1397 / 2
_STEEL_BODY = Layer(_R_ID, _R_OD, K_STEEL)

CONVENTIONAL = PipeType(
    "conventional", (_STEEL_BODY,),
    source="Wu et al. (2025), Table 5")
INTERNALLY_COATED = PipeType(
    "internally coated", (Layer(0.1194 / 2, _R_ID, 0.47), _STEEL_BODY),
    source="1 mm composite at 0.47 W/m/K, Vetsak et al. (2024) via Wu et al. (2025)")
EXTERNALLY_COATED = PipeType(
    "externally coated", (_STEEL_BODY, Layer(_R_OD, 0.1417 / 2, 1.31)),
    source="1 mm composite at 1.31 W/m/K, Vetsak et al. (2024) via Wu et al. (2025)")
# Wu et al. give the dual-wall pipe an 80 mm bore and a phenolic fill at
# 0.092 W/m/K. The fill is taken to run from the bore to the outer tube, which
# has the conventional pipe's wall; the inner tube's steel adds nothing measurable.
DUAL_WALL = PipeType(
    "dual-wall", (Layer(0.080 / 2, _R_ID, 0.092), _STEEL_BODY),
    source="phenolic fill at 0.092 W/m/K, Xiao et al. (2022) via Wu et al. (2025)")
# Vacuum-insulated tubing is a production tubular; no source describes a drill
# pipe built that way. Its apparent conductivity, averaged over body and
# connections, is placed across the dual-wall pipe's envelope, as a sensitivity.
VACUUM_INSULATED = PipeType(
    "vacuum-insulated (sensitivity)", (Layer(0.080 / 2, _R_ID, 0.04), _STEEL_BODY),
    source="apparent 0.04 W/m/K with connections, range 0.02-0.08, Dufrene et al. (2025)")

PIPE_TYPES = {p.name: p for p in
              (CONVENTIONAL, INTERNALLY_COATED, EXTERNALLY_COATED, DUAL_WALL,
               VACUUM_INSULATED)}

# The wall the model used before it had named pipe types: 10 mm at 0.02 W/m/K
# on a 100 mm bore, the best-case end of the vacuum-tubing range. Kept so the
# v1.0 site table and the regression anchors reproduce.
LEGACY_VACUUM = uniform_wall(0.02, name="legacy vacuum tubing, 0.02 W/m/K")


@dataclass(frozen=True)
class StringSegment:
    length: float
    pipe: PipeType


@dataclass(frozen=True)
class HoleInterval:
    """From md_top to md_bottom: the annulus's outer radius r_wall and the
    layers outside it, up to the rock. layers=() is open hole."""
    md_top: float
    md_bottom: float
    r_wall: float
    layers: tuple = ()

    def __post_init__(self):
        if self.layers:
            if not np.isclose(self.layers[0].r_in, self.r_wall):
                raise ValueError("first layer must start at r_wall")
            _check_concentric(self.layers)

    @property
    def r_rock(self):
        return self.layers[-1].r_out if self.layers else self.r_wall

    def resistance(self):
        return sum(l.resistance() for l in self.layers)


def open_hole(md_top, md_bottom, d_hole):
    return HoleInterval(md_top, md_bottom, d_hole / 2)


def cased(md_top, md_bottom, casing_id, casing_od, d_hole, k_cement=K_CEMENT):
    """A single casing string cemented into a hole of diameter d_hole."""
    return HoleInterval(md_top, md_bottom, casing_id / 2,
                        (Layer(casing_id / 2, casing_od / 2, K_STEEL),
                         Layer(casing_od / 2, d_hole / 2, k_cement)))


class Survey:
    """Measured depth to true vertical depth. Built from survey stations
    (minimum curvature) or from a table of TVD against MD; interpolated
    linearly between stations."""

    def __init__(self, md, tvd):
        self.md = np.asarray(md, float)
        self.tvd_at = np.asarray(tvd, float)
        if np.any(np.diff(self.md) <= 0):
            raise ValueError("survey MD must increase")

    @classmethod
    def vertical(cls):
        return cls([0.0, 1.0], [0.0, 1.0])

    @classmethod
    def from_stations(cls, md, inc_deg, azi_deg=None, tvd0=0.0):
        md = np.asarray(md, float)
        I = np.radians(np.asarray(inc_deg, float))
        A = np.radians(np.zeros_like(md) if azi_deg is None
                       else np.asarray(azi_deg, float))
        I1, I2, A1, A2 = I[:-1], I[1:], A[:-1], A[1:]
        cos_dl = np.cos(I2 - I1) - np.sin(I1) * np.sin(I2) * (1 - np.cos(A2 - A1))
        dl = np.arccos(np.clip(cos_dl, -1.0, 1.0))
        rf = np.ones_like(dl)
        big = dl > 1e-9
        rf[big] = 2 / dl[big] * np.tan(dl[big] / 2)
        dtvd = np.diff(md) / 2 * (np.cos(I1) + np.cos(I2)) * rf
        return cls(md, tvd0 + np.concatenate([[0.0], np.cumsum(dtvd)]))

    def tvd(self, s):
        s = np.asarray(s, float)
        # beyond the last station, carry on along the last station's slope
        out = np.interp(s, self.md, self.tvd_at)
        slope = (self.tvd_at[-1] - self.tvd_at[-2]) / (self.md[-1] - self.md[-2])
        return np.where(s > self.md[-1], self.tvd_at[-1] + slope * (s - self.md[-1]), out)


class WellGeometry:
    """The string and hole, with per-node lookups vectorised over s."""

    def __init__(self, string, hole, survey=None):
        self.string = tuple(string)
        self.hole = tuple(hole)
        self.survey = survey or Survey.vertical()
        self.md_bit = float(sum(seg.length for seg in self.string))
        self._seg_ends = np.cumsum([seg.length for seg in self.string])
        self._hole_tops = np.array([h.md_top for h in self.hole])
        if self._hole_tops[0] > 0 or self.hole[-1].md_bottom < self.md_bit - 1e-6:
            raise ValueError("hole intervals must cover the string")
        for a, b in zip(self.hole, self.hole[1:]):
            if not np.isclose(a.md_bottom, b.md_top):
                raise ValueError("hole intervals must be contiguous")
        for h in self.hole:
            for seg in self.string:
                if seg.pipe.r_out >= h.r_wall:
                    raise ValueError(f"{seg.pipe.name} does not fit inside r_wall {h.r_wall}")

        p = [seg.pipe for seg in self.string]
        self._r_bore = np.array([q.r_bore for q in p])
        self._r_out = np.array([q.r_out for q in p])
        self._R_body = np.array([q.body_resistance() for q in p])
        self._R_joint = np.array([q.joint_resistance() for q in p])
        self._f_joint = np.array([q.joint_fraction for q in p])
        self._r_wall = np.array([h.r_wall for h in self.hole])
        self._r_rock = np.array([h.r_rock for h in self.hole])
        self._R_hole = np.array([h.resistance() for h in self.hole])

    @classmethod
    def single(cls, length, pipe, r_hole=C.R_WELL, survey=None):
        """One string segment in open hole: the model's original layout."""
        return cls([StringSegment(length, pipe)], [HoleInterval(0.0, length, r_hole)], survey)

    def segment_index(self, s):
        i = np.searchsorted(self._seg_ends, np.asarray(s, float), side="left")
        return np.clip(i, 0, len(self.string) - 1)

    def hole_index(self, s):
        i = np.searchsorted(self._hole_tops, np.asarray(s, float), side="right") - 1
        return np.clip(i, 0, len(self.hole) - 1)

    def at(self, s):
        """Per-node geometry at measured depths s."""
        i, j = self.segment_index(s), self.hole_index(s)
        r_bore, r_out, r_wall = self._r_bore[i], self._r_out[i], self._r_wall[j]
        return dict(
            tvd=self.survey.tvd(s),
            r_bore=r_bore, r_out=r_out, r_wall=r_wall, r_rock=self._r_rock[j],
            R_body=self._R_body[i], R_joint=self._R_joint[i], f_joint=self._f_joint[i],
            R_hole=self._R_hole[j],
            A_bore=np.pi * r_bore ** 2, Dh_bore=2 * r_bore,
            A_ann=np.pi * (r_wall ** 2 - r_out ** 2), Dh_ann=2 * (r_wall - r_out),
        )


# The deep vertical well of Wu et al. (2025), Tables 3 and 4: three casings
# cemented to surface, 0.22 m (8.67-inch) open hole, and a string of drill pipe,
# heavy-weight pipe and collars. Casings run to surface, the gap between two
# casings is cement, and rock starts at the outermost casing (the drilled hole
# sizes aren't given).
WU_HWDP = PipeType("heavy-weight drill pipe", (Layer(0.083 / 2, 0.14 / 2, K_STEEL),))
WU_COLLAR = PipeType("drill collar", (Layer(0.073 / 2, 0.17 / 2, K_STEEL),))
WU_L_HWDP, WU_L_COLLAR, WU_L_BIT = 223.0, 143.0, 0.8
WU_SHOES = (100.0, 1000.0, 4000.0)   # conductor, surface, intermediate [m]


def wu2025_hole(md_bottom, shoes=WU_SHOES):
    st = lambda a, b: Layer(a, b, K_STEEL)
    ce = lambda a, b: Layer(a, b, K_CEMENT)
    inter, surf, cond = (0.11, 0.12), (0.16, 0.17), (0.24, 0.255)
    s_cond, s_surf, s_int = shoes
    return [
        HoleInterval(0.0, s_cond, inter[0], (st(*inter), ce(inter[1], surf[0]), st(*surf),
                                            ce(surf[1], cond[0]), st(*cond))),
        HoleInterval(s_cond, s_surf, inter[0], (st(*inter), ce(inter[1], surf[0]), st(*surf))),
        HoleInterval(s_surf, s_int, inter[0], (st(*inter),)),
        open_hole(s_int, md_bottom, 0.22),
    ]


def wu2025_string(md_bit, pipe):
    l_dp = md_bit - WU_L_HWDP - WU_L_COLLAR - WU_L_BIT
    return [StringSegment(l_dp, pipe), StringSegment(WU_L_HWDP, WU_HWDP),
            StringSegment(WU_L_COLLAR + WU_L_BIT, WU_COLLAR)]


def wu2025_well(md_bit, pipe, survey=None, shoes=WU_SHOES):
    return WellGeometry(wu2025_string(md_bit, pipe), wu2025_hole(md_bit, shoes), survey)
