# Wu et al. (2025): benchmark inputs and results

Wu, Y., Zhang, Y., Bettir, N., Ashok, P., van Oort, E. (2025). "A Comprehensive
Evaluation of Drill Pipe Insulation for Downhole Temperature Management Using
Physics-Based Models." *Proceedings, 50th Workshop on Geothermal Reservoir
Engineering*, Stanford University, SGP-TR-202.
https://pangea.stanford.edu/ERE/db/GeoConf/papers/SGW/2025/Wu1.pdf

Their model is the transient thermo-hydraulic model of Fallah et al. (2021),
Gu et al. (2022) and Ma et al. (2016), run for 1,500 minutes to reach steady
state. They don't validate it against field data in this paper.

Values in the tables below come from the paper's text and tables where it gives
them. Values marked (fig.) were read off its figures (Figs. 2, 5, 6, 7 and 9),
to about ±3 °C and ±2 MPa.

## Inputs (Tables 2 to 5)

| | Vertical well | Horizontal well |
|---|---|---|
| Hole size | 0.22 m (8.67 in) | 0.22 m (8.67 in) |
| True vertical depth | 8,600 m | 8,600 m |
| Inclination | 0° | 90° in the lateral |
| Lateral | none | 2,000 m, after a 400 m curve |

Drilling fluid (water-based mud): density 1,078 kg/m³, plastic viscosity
14 mPa·s, cp 3,750 J/kg·K, k 0.75 W/m·K. Flow 136 m³/h (600 gal/min); inlet
40 °C.

Formation: surface 17 °C, gradient 24 °C/km (about 223 °C at 8.6 km), density
2,800 kg/m³, cp 930 J/kg·K, k 2.31 W/m·K.

Casing (all steel, 45 W/m·K):

| | ID | OD | Length |
|---|---|---|---|
| Conductor | 0.48 m | 0.51 m | 100 m |
| Surface | 0.32 m | 0.34 m | 1,000 m |
| Intermediate | 0.22 m | 0.24 m | 4,000 m |
| Open hole | | 0.22 m | 4,600 m |

String (steel, 45 W/m·K):

| | ID | OD | Length |
|---|---|---|---|
| Drill pipe | 0.121 m | 0.14 m | 8,224 m |
| Heavy-weight drill pipe | 0.083 m | 0.14 m | 223 m |
| Drill collars | 0.073 m | 0.17 m | 143 m |
| Bit | | 0.22 m | 0.8 m |

Pipe types (Table 5):

| Type | ID | OD | Insulation | Conductivity | Their source |
|---|---|---|---|---|---|
| Conventional | 121.4 mm | 139.7 mm | none | 45 W/m·K | Xiao et al. (2024) |
| Internally coated | 119.4 mm | 139.7 mm | 1 mm composite fibre-resin | 0.470 W/m·K | Vetsak et al. (2024) |
| Externally coated | 121.4 mm | 141.7 mm | 1 mm composite fibre-resin | 1.310 W/m·K | Vetsak et al. (2024) |
| Dual-wall | 80.0 mm | 139.7 mm | phenolic resin | 0.092 W/m·K | Xiao et al. (2022) |

Pump limit: 7,500 psi (51.7 MPa), after Drillmec (2025).

## Results

Steady BHCT at 1,500 minutes, base case (Fig. 2):

| Pipe | Vertical | Horizontal |
|---|---|---|
| Conventional | 179 °C (fig.) | 206 °C (fig.) |
| Internally coated | 76 °C (fig.) | 100 °C (text) |
| Externally coated | 121 °C (fig.) | 149–153 °C (fig.; text "153") |
| Dual-wall | 59 °C (fig.) | 65 °C (text) |

BHCT against reservoir temperature, 600 gal/min (Fig. 5; Table 6 gives the
gradients: 21.3, 27.1, 32.9, 38.7 and 44.5 °C/km for 200 to 400 °C):

| Reservoir | 200 °C | 250 °C | 300 °C | 350 °C | 400 °C |
|---|---|---|---|---|---|
| Vertical, internally coated | 71 (text) | 80 | 89 | 98 | 107 (text) |
| Vertical, externally coated | 110 | 114 | 132 | 148 | 165 |
| Vertical, dual-wall | 58 | 59 | 60 | 61 | 62 (text: "steady at 59") |
| Horizontal, internally coated | 92 | 101 | 109 | 122 | 134 |
| Horizontal, externally coated | 135 | 145 | 170 | 194 | 232 |
| Horizontal, dual-wall | 62 | 63 | 63 | 64 | 65 |

The 400 °C case against flow rate (Fig. 6), vertical well:

| Flow, gal/min | 400 | 500 | 600 | 700 | 800 |
|---|---|---|---|---|---|
| Internally coated | 145 | 123 | 107 | 95 | 87 |
| Externally coated | 212 | 186 | 165 | 148 | 135 |
| Dual-wall | 59 | 60 | 61 | 63 | 66 |

The paper gives no conventional-pipe result at 400 °C.

Standpipe pressure, vertical well (Figs. 7 and 9): about 16 to 18 MPa for both
coated pipes and about 67 MPa for full dual-wall pipe at 600 gal/min (fig.).
Dual-wall pipe exceeds the 51.7 MPa limit above 500 gal/min (text); the coated
pipes stay below it from 400 to 800 gal/min (text). No conventional-pipe value
is given.

## Corrections to what was expected

The vertical-well BHCTs had been summarised as about 115, 100, 135 and 85 °C
(conventional, internal, external, dual-wall). The figure gives about 179, 76,
121 and 59 °C. The horizontal conventional-pipe value is about 206 °C, not
170 °C.
