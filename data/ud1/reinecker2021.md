# Reinecker et al. (2021) — verified source values

Reinecker, J., Gutmanis, J., Foxford, A., Cotton, L., Dalby, C., Law, R. (2021).
"Geothermal exploration and reservoir modelling of the United Downs deep
geothermal project, Cornwall (UK)." *Geothermics* 97, 102226.
https://doi.org/10.1016/j.geothermics.2021.102226 — open access, CC BY 4.0.

Every figure below was read off the paper directly (section numbers given), not
from a summary. Values the v1.1 spec carried that turned out to be wrong are
flagged **CORRECTION**.

## Stress (section 6.3)

| Quantity | Paper | How it was obtained |
|---|---|---|
| `Sv` | **25.275 MPa/km** | best fit to wireline density logs |
| `Shmin` | **13.21 MPa/km + 3 MPa** | linear fit to refracture pressures from the Rosemanowes hydrofracture tests (Pine et al. 1983b, "to depths of 2000 m"), *not* a UD-1 measurement |
| `SHmax` | **25.99 MPa/km + 5.9 MPa** | derived from Shmin, formation pressure and the fracture fluid pressure needed to initiate slip, **assuming a coefficient of friction μ = 0.8** |
| Pore pressure | **9.494 MPa/km below a static fluid level at ~61 m below ground level** | reported directly |
| Regime | strike-slip, "critically stressed for shearing on appropriately oriented fractures" | |
| SH azimuth | 134° ± 25° (breakouts), 134° ± 12° (tensile fractures) | |

The paper's own slip-tendency analysis (Fig. 12) likewise uses "a reasonable
static friction coefficient of 0.8" (citing Byerlee 1978).

**CORRECTION — pore pressure.** The spec assumed hydrostatic unless reported. It
is reported, and 9.494 MPa/km below 61 m is not the same as the model's
10 MPa/km from surface: at 5 km it is 46.9 MPa against 50.0 MPa.

**CORRECTION — μ.** The published SHmax is not independent evidence for μ ≈ 0.85.
It was *constructed* at frictional equilibrium with μ = 0.8, so reading μ back
out of it is circular. Under the paper's own pore pressure the profile implies
μ = 0.75 at 5 km and 0.79 at 12.5 km, i.e. the 0.8 it was built with. The
spec's "0.83 to 0.86" comes from using 10 MPa/km hydrostatic instead.

### Effective S1/S3 against the frictional cap

Caps R(μ) = (√(1+μ²)+μ)²: **3.12** at 0.6, **4.33** at 0.8, **4.68** at 0.85.

| Depth | Pp = paper | S1/S3 | implied μ | Pp = 10 MPa/km | S1/S3 | implied μ |
|---|---|---|---|---|---|---|
| 5,000 m | 46.9 MPa | 4.01 | 0.75 | 50.0 MPa | 4.51 | 0.83 |
| 5,058 m | 47.4 MPa | 4.02 | 0.75 | 50.6 MPa | 4.51 | 0.83 |
| 12,500 m | 118.1 MPa | 4.25 | 0.79 | 125.0 MPa | 4.77 | 0.86 |

SHmax/Shmin is 1.97 at every depth (both intercepts are small), as the spec says.

## Well (sections 5, 6.3; Figs. 5, 8)

- UD-1: 5,275 m MD, **5,058 m TVD**. UD-2: 2,393 m MD, 2,214 m TVD.
- Hole sections: 24" to 247 m, 17.5" to 900 m, 12.25" to 4,000 m, 8.5" to 5,275 m.
- Kick-off at 3,390 m MD, building to ~35° inclination below 4 km. The 12.25"
  section is **nearly vertical** — max 15° deviation, only in its lowest 200 m.
  So MD ≈ TVD over 900–4,000 m, but the 8.5" section is strongly deviated and a
  vertical-hole Kirsch solution does not apply there.
- Temperature "around 180 °C at 5 km depth" (abstract, conclusions).

**CORRECTION — mud weight by section.** The paper (section 5): mud weight was kept
"as low as possible (**1.10 SG in the 24" and 17.5" sections and below 1.05 SG in
the 12.25" and 8.5" sections**) to reduce losses but still enable effective hole
cleaning." The spec had this the other way round. The logged breakout interval
(900–4,000 m) is the 12.25" section, so its mud weight is **below 1.05 SG**, not
1.10 SG.

## Borehole wall failure (section 6.3, Fig. 8)

- **27 breakouts, 139 m total, between c. 900 and c. 4,000 m MD** (12.25" section).
- **13 drilling-induced tensile fractures, 168 m total, between c. 2,500 and
  3,700 m MD.** (The spec gave the counts but not this depth range.)
- Breakouts are "locally developed in the 12.25" section, notably in the lower
  part, and **extensively in the 8.5" section**" — i.e. the 4,000–5,275 m
  interval broke out extensively and the well still reached TD. Orientations
  there were not used because the inclination exceeds 15°.
- Breakouts "are narrow and weak indicating high horizontal stress anisotropy".
- Scatter in breakout orientation rises above c. 3,600 m MD, attributed to
  increasing inclination and active fault structures.

## Mud losses and hydraulics (sections 5, 6.4)

- Dynamic losses up to **11 m³/h** in the 8.5" section below 4,000 m, especially
  4,850–5,080 m MD, around an open fracture zone at 4,890 m MD.
- Losses occurred in both wells; no measures beyond reducing circulation rate.

## Decision taken on the back of this (D1, revised 3 October)

The friction coefficient defaults to **μ = 0.8**, the value Reinecker et al. use
to derive SHmax and to assess slip tendency, with **0.6 and 0.85 reported as
sensitivities**. Caps: 3.12, 4.33, 4.68.

Consequence for C2's "done when": the published profile, extrapolated to
12.5 km, is **admissible at every one of those μ values** under the paper's pore
pressure (effective S1/S3 = 4.25). The spec's expectation that the cap binds
there holds only with 10 MPa/km hydrostatic Pp *and* μ = 0.85 (4.77 against
4.68). C2 should still implement the cap and report where it binds; it just does
not bind on United Downs.
