# Newberry well NWG 55-29

Data for the stability cross-check against NWG 55-29 (`src/crosscheck_newberry.py`)
and the Newberry site in the site table.

## Sources and licence

- **Davatzes, N.C., Hickman, S.H. (2011).** "Preliminary Analysis of Stress in
  the Newberry EGS Well NWG 55-29." *GRC Transactions* 35. GDR submission 737.
  Notes in `davatzes_hickman_2011.md`.
- **Newberry EGS Demonstration: Well 55-29 Stimulation Data.** Cladouhos, T.T.
  (AltaRock Energy). GDR submission 271, https://doi.org/10.15121/1357916.
  CC-BY 4.0. The extracts here are derived from it.

The raw files are not committed. To rebuild the extracts, download these into
`data/newberry/raw/gdr271/` and run `python data/newberry/build_extracts.py`:

| File (GDR 271) | SHA-256 (first 16) |
|---|---|
| Static_Profile_Oct2008_PT.csv | ac30b6ab5050c8fc |
| Deviation_Data.csv | 21f150d8e1d11eb8 |
| 55_29_RUN_3.las | 86c34f536a948864 |
| 55-29_Heat_Flow.pdf (read, not extracted) | be4d8fb98da2725a |

The paper itself is `raw/davatzes-and-hickman-2011.pdf` (8d357bfc805d07f4).

## Extracts

- `nwg5529_static_2008.csv`: the static pressure and temperature survey of
  October 2008, every 25 ft of MD, with TVD from the deviation survey. It
  reaches 331.4 °C at 9,993 ft MD (2,989 m TVD). Its pressures are the
  formation's: 3,120 psi at 8,420 ft against the paper's Pf of 3,121.6 psi.
- `nwg5529_deviation.csv`: MD and TVD. The well deviates 10–15° over the
  imaged interval (the paper); 8,420 ft MD is 8,285 ft TVD.
- `nwg5529_open_hole_logs.csv`: neutron porosity, bulk density, caliper and
  gamma ray from the open-hole run (Halliburton, 13 July 2008), as medians over
  2 ft, with depth moved from the kelly bushing to ground level (31 ft, the
  LAS header's depth above permanent datum) to match the paper's MD GL.

## Thermal conductivity

AltaRock's "Well 55-29 Heat Flow Values" gives conductivities by depth
interval, corrected for temperature after Birch and Clark (1940): 1.99 to
2.47 W/m·K, and 1.99 to 2.30 W/m·K over the open hole below 6,500 ft. The site
uses the mean over the open hole, 2.14 W/m·K.

## What the cross-check found in the data

- The paper's Shmin at 8,420 ft, 4,675.3 psi, is what its frictional-equilibrium
  formula gives at μ = 0.70; at μ = 0.55, which the text calls the best
  estimate, the formula gives 5,118 psi. Both cases are run.
- The paper's SHmax fit, 0.772 z + 1,112.6 psi, gives 7,508 psi at 8,420 ft
  (z in TVD ft) against the 7,305 psi it quotes there, within the scatter
  of SHmax from individual breakouts.
- From the porosity relation, the logged granodiorite (8,807–8,860 ft) has
  60–82 MPa (median 70), the volcanics 51–74 MPa (median 61). That is too
  little contrast to explain breakouts in one and none in the other, under
  either the paper's criterion or this model's.
