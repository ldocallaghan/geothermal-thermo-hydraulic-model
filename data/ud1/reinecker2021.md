# Reinecker et al. (2021): United Downs source values

Reinecker, J., Gutmanis, J., Foxford, A., Cotton, L., Dalby, C., Law, R. (2021).
"Geothermal exploration and reservoir modelling of the United Downs deep
geothermal project, Cornwall (UK)." *Geothermics* 97, 102226.
https://doi.org/10.1016/j.geothermics.2021.102226. Open access, CC BY 4.0.

The values the model takes from the paper, with section numbers.

## Stress (section 6.3)

| Quantity | Paper | How it was obtained |
|---|---|---|
| `Sv` | **25.275 MPa/km** | best fit to wireline density logs |
| `Shmin` | **13.21 MPa/km + 3 MPa** | linear fit to refracture pressures from the Rosemanowes hydrofracture tests (Pine et al. 1983b, "to depths of 2000 m"), *not* a UD-1 measurement |
| `SHmax` | **25.99 MPa/km + 5.9 MPa** | derived from Shmin, formation pressure and the fracture fluid pressure needed to initiate slip, **assuming a coefficient of friction μ = 0.8** |
| Pore pressure | **9.494 MPa/km below a static fluid level at ~61 m below ground level** | reported directly |
| Regime | strike-slip, "critically stressed for shearing on appropriately oriented fractures" | |
| SH azimuth | 134° ± 25° (breakouts), 134° ± 12° (tensile fractures) | |

The paper's slip-tendency analysis (Fig. 12) also uses "a reasonable static
friction coefficient of 0.8" (citing Byerlee 1978).

**Pore pressure.** 9.494 MPa/km below 61 m is lighter than the model's default
10 MPa/km from surface: 46.9 MPa against 50.0 MPa at 5 km. The United Downs
stress profile carries the paper's value.

**Friction.** Because SHmax was derived at frictional equilibrium with μ = 0.8,
the profile can't be used as independent evidence for any friction value:
under the paper's own pore pressure it implies μ = 0.75 at 5 km and 0.79 at
12.5 km, the 0.8 it was built with. Substituting 10 MPa/km hydrostatic pore
pressure would give 0.83 to 0.86 instead, close to Byerlee's 0.85, but that is
an artefact of the substitution. The model's default μ is the paper's 0.8,
with 0.6 and 0.85 reported as sensitivities.

### Effective S1/S3 against the frictional cap

Caps R(μ) = (√(1+μ²)+μ)²: **3.12** at 0.6, **4.33** at 0.8, **4.68** at 0.85.

| Depth | Pp = paper | S1/S3 | implied μ | Pp = 10 MPa/km | S1/S3 | implied μ |
|---|---|---|---|---|---|---|
| 5,000 m | 46.9 MPa | 4.01 | 0.75 | 50.0 MPa | 4.51 | 0.83 |
| 5,058 m | 47.4 MPa | 4.02 | 0.75 | 50.6 MPa | 4.51 | 0.83 |
| 12,500 m | 118.1 MPa | 4.25 | 0.79 | 125.0 MPa | 4.77 | 0.86 |

SHmax/Shmin is 1.97 at every depth (both intercepts are small). Extrapolated
to 12.5 km under the paper's pore pressure, the profile stays inside the cap at
μ = 0.8 and 0.85, so the frictional cap does not bind on United Downs; it would
only with 10 MPa/km hydrostatic pore pressure and μ = 0.85.

## Well (sections 5, 6.3; Figs. 5, 8)

- UD-1: 5,275 m MD, **5,058 m TVD**. UD-2: 2,393 m MD, 2,214 m TVD.
- Hole sections: 24" to 247 m, 17.5" to 900 m, 12.25" to 4,000 m, 8.5" to 5,275 m.
- Kick-off at 3,390 m MD, building to ~35° inclination below 4 km. The 12.25"
  section is **nearly vertical**, deviating at most 15° and only in its lowest
  200 m. So MD ≈ TVD over 900–4,000 m, but the 8.5" section is strongly
  deviated and a vertical-hole Kirsch solution does not apply there.
- Temperature "around 180 °C at 5 km depth" (abstract, conclusions).
- **Mud weight by section** (section 5): kept "as low as possible (1.10 SG in
  the 24" and 17.5" sections and below 1.05 SG in the 12.25" and 8.5"
  sections) to reduce losses but still enable effective hole cleaning". The
  logged breakout interval (900–4,000 m) is the 12.25" section, so its mud
  weight is **below 1.05 SG**.

## Borehole wall failure (section 6.3, Fig. 8)

- **27 breakouts, 139 m total, between c. 900 and c. 4,000 m MD** (12.25" section).
- **13 drilling-induced tensile fractures, 168 m total, between c. 2,500 and
  3,700 m MD.**
- Breakouts are "locally developed in the 12.25" section, notably in the lower
  part, and extensively in the 8.5" section": the 4,000–5,275 m interval broke
  out extensively and the well still reached TD. Orientations there were not
  used because the inclination exceeds 15°.
- The breakouts "are narrow and deep indicating high horizontal stress
  anisotropy".
- Scatter in breakout orientation rises above c. 3,600 m MD, attributed to
  increasing inclination and active fault structures.

The BGS interpretation of the same image log (`bgs_image_log.md`) gives
different counts.

## Mud losses and hydraulics (sections 5, 6.4)

- Dynamic losses up to **11 m³/h** in the 8.5" section below 4,000 m, especially
  4,850–5,080 m MD, around an open fracture zone at 4,890 m MD.
- Losses occurred in both wells; no measures beyond reducing circulation rate.
