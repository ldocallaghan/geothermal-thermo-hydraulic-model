# FORGE well 16B(78)-32: insulated drill pipe trial

Data for validating the circulation model (Model 1) against the insulated
drill pipe (IDP) trial Eavor ran in Utah FORGE well 16B(78)-32 in May 2023.

## Sources and licence

- **Utah FORGE: Well 16B(78)-32 Drilling Data.** McLennan, J., Mock, B.,
  Swearingen, L., Baldwin, R., Hodder, M., Vetsak, A., Kuhns, A.T., Breland,
  J., England, K. (2023). Geothermal Data Repository, submission 1516.
  https://doi.org/10.15121/1998591. CC-BY 4.0.
- **Utah FORGE: Deep Wells Temperature Surveys as of September 2022.** Jones, C.
  (2022). Geothermal Data Repository, submission 1421.
  https://doi.org/10.15121/1893533. CC-BY 4.0.

Both are used under CC-BY 4.0. The extracts here are derived from them.

The raw files are not committed. To rebuild the extracts, download these into
`data/forge16b/raw/` (unzipping `16B_Pason.zip` to `raw/16B_Pason/`) and run
`python data/forge16b/build_extracts.py`:

| File (GDR) | SHA-256 (first 16) |
|---|---|
| Eavor FORGE Insulated Drill Pipe Trial Results - 23 June 2023 (003).pdf (1516) | b392f98198f3b661 |
| 16B(78)-32 Well Survey.zip (1516) | 38c59c9452da709f |
| 16B Daily Reports.zip (1516) | 2c2475d95be35a1f |
| 16B mud temp logs.zip (1516) | ecefb1fb69ecda06 |
| 16B_Pason.zip (1516) | 858a6eaf56c2cb3c |
| End of Well Report-16B78-32-May_2024.pdf (1516) | b35f25eda60fcec9 |
| Utah FORGE Deep Well Temperature Profiles_Sept 2022.xlsx (1421) | c1a35f7914b5f591 |
| papers/xing_sgw2026.pdf (Stanford SGW 2026, Xing et al., for the rock properties) | 6d345bd2869d25fd |

## Extracts

- `survey_16b.csv`: the gyro survey (minimum curvature), 428 stations. TVD is
  from the kelly bushing, 32 ft above ground level (KB 5,447.65 ft, GL
  5,415.65 ft). TD is 10,947 ft MD, 8,391 ft TVD below KB, 8,359 ft below GL.
- `formation_temperature_16a.csv`: the cement-bond-log temperature in the
  neighbouring well 16A(78)-32 on 16 August 2021, about 216 days after the rig
  left, every 25 ft of TVD. The spreadsheet notes that its TVD comes from a
  polynomial fit to MD, usually within 30 to 40 ft of the directional survey.
  At 16B's TD (about 8,357 ft TVD) it reads 425.7 °F (219 °C). It agrees with
  the 58-32 survey after 1,371 days of recovery to within 3 °F at 7,500 ft, so
  it is close to equilibrium. 16B runs parallel to 16A, about 300 ft above it.
- `run_observations.csv`: the four trial runs, transcribed from the Eavor trial
  report (text and its Table 1): depths, string make-up, flow, inlet temperature
  and MWD temperatures. Blank cells are not given in the report; flow and inlet
  temperature for BHA 12 and 13 will come from the Pason data.

- `rig_hydraulics.csv`: the rig's hydraulics for each trial run, transcribed
  from the daily drilling reports reproduced in the End of Well Report (reports
  36 to 42, 22 to 28 May 2023): flow (average and maximum over the interval),
  standpipe pressure, mud weight, bit nozzles in 32nds of an inch, and the
  rig's own bit pressure drop and annular velocities around the drill pipe and
  the BHA. The End of Well Report numbers these BHAs 17 to 20; the trial report
  and the drilling contractors call them 10 to 13. Blank cells are not given.
  The mud was a low-solids water-based mud at 8.40 lb/gal (99.5% water).
  The standpipe pressure includes the mud motor and MWD, which a friction and
  nozzle model does not, so it is an upper bound for that model; the bit
  pressure drop and annular velocities are direct checks. The rig's annular
  velocity around the drill pipe is 24.51 Q / (9.5² − 5.5²) ft/min, Q in gal/min.
  The drill pipe is 5-1/2-inch, 24.7 lb/ft, S-135 (daily reports).

- `string_components.csv`: the BHA, collars and heavy-weight pipe, with OD, ID
  and length, from the BHA tables in the End of Well Report (Fig. 90 for BHA
  #17, trial BHA 10, and Fig. 102 for BHA #19, trial BHA 12; the trial report
  says BHA 11 used the same BHA as BHA 10). BHA #19's table puts the
  directional survey sensor 84 ft above the bit; that is taken as the MWD
  temperature sensor's position in all four runs. It also lists "25 STD Drill
  pipe" between the heavy-weight pipe and a crossover to the IDP, which
  settles how BHA 12's string was made up: the trial report's 8,264.49 ft of
  "IDP" is the whole drill pipe (it closes the depth tally with the BHA and
  heavy-weight pipe), of which the bottom 2,378.5 ft was regular pipe, leaving
  5,886 ft (71%) of IDP. The trial report's "~70% IDP" agrees.
  The End of Well Report lists 30 joints of heavy-weight pipe for BHA #17
  where the trial report gives 28; the trial report's tally closes on the run
  depths, so its lengths are used, with the End of Well Report's dimensions.
- `pason_trial_1min.csv`: the Pason record from 21 to 29 May 2023 averaged
  per minute: hole and bit depth, standpipe pressure, pump output, inlet and
  return mud temperature, rate of penetration and motor differential pressure.
- `drilling_history_10min.csv`: the whole Pason record in 10-minute bins: the
  deepest hole reached so far, and the fraction of the bin spent circulating
  (pump output above 100 gal/min).

- `cycle_timings.csv`: the drilling cycle measured from the 1-minute Pason
  extract (built by `cycle_timings()` in `build_extracts.py`): every gap in
  circulation with the bit on bottom (trips, and pauses with the bit on
  bottom), each trip's speed out and in with stationary spells over ten
  minutes excluded, the time at surface, staged-circulation stops on the way
  in, connections while drilling (pumps off for under an hour with the bit on
  bottom on both sides), and the bit runs. The Pason bit depth stops updating
  at about 900 ft while the BHA is handled, so a bit shallower than 1,000 ft
  counts as out of the hole. In the Pason record the first trip runs from
  22:38 on 21 May to 23:17 on 22 May.
- `connections_10s.csv`: every connection while drilling over the whole well,
  from the 10-second record (built by `connections_10s()`): the pumps off for
  under an hour, at full flow in the minute before, with the bit within 100 ft
  of bottom on both sides and new hole made in the half hour before. For each,
  how long the pumps were off, and how long the fluid circulated between the
  last new hole and the pumps stopping. The 10-second record resolves both
  better than the 1-minute extract, which also covers only the trial week.
- `trip_surface_operations.csv`: what was done at surface on each trip, from
  the daily reports in the End of Well Report, and whether it was routine (a
  bit and BHA change and nothing else).

## Rock thermal properties

Measured by MetaRock Laboratories (2021) on three FORGE granitoid cores, as
tabulated by Xing, P., Jones, C., Damjanac, B., Simmons, S., Moore, J., Deo, M.,
McLennan, J. (2026), "Influence of Thermal Properties Heterogeneity on
Geothermal Production and Fracture Spacing for Enhanced Geothermal System",
*Proc. 51st Workshop on Geothermal Reservoir Engineering*, Stanford,
SGP-TR-230, Table 1 (the underlying dataset is GDR submission 1430):

| Sample | Density | k at 35 / 100 / 200 °C | cp at 35 / 100 / 200 °C |
|---|---|---|---|
| 58-32A (7,440 ft) | 2.68 g/cc | 2.165 / 2.149 / 2.041 W/m·K | 798.7 / 864.8 / 981.3 J/kg·K |
| 16A(78)-32 | 2.60 g/cc | 2.389 / 2.341 / 2.206 W/m·K | 785.8 / 862.3 / 944.5 J/kg·K |
| 58-32B (6,804 ft) | 2.65 g/cc | 3.087 / 3.002 / 2.672 W/m·K | 988.2 / 1091.1 / 1175.2 J/kg·K |

The validation uses the means at 100 °C, the temperature of the cooled rock
near the wall: k = 2.50 W/m·K and ρc = 2,643 × 939 = 2.48 MJ/m³·K. The
spread across samples and temperatures, 2.04 to 3.09 W/m·K, is run as a
sensitivity. Gwynn et al. (2019, UGS Misc. Pub. 169-L) report 2.0 to
3.9 W/m·K on 58-32 cuttings at room temperature, depending on quartz content
(as cited by Xing et al.).

## What the dataset does and doesn't contain

- The Pason 10-second record (1.8 GB) has flow rate, hole and bit depth,
  standpipe pressure, inclination, and inlet ("TEMP IN MANIFOLD") and return
  ("TEMP OUT FLOW") mud temperature. It has no downhole temperature.
- The mud-temperature LAS files give inlet and return temperature against depth.
- The daily reports give tool temperatures only while tripping, during
  circulate-and-cool stops.
- **The MWD temperature while drilling is therefore available only from the
  Eavor trial report**: run averages, BHA 12's levelled range and BHA 13's
  steady values in the text, and plotted points in its figures. The validation
  uses the text values.
- The trial report gives no dimensions or conductivity for the IDP ("different
  versions, grades and coating types"). Its coating conductivities are reported
  by Vetsak et al. (2024) as 0.47 W/m·K (internal) and 1.31 W/m·K (external),
  1 mm thick, as cited by Wu et al. (2025).

## Casing and hole (End of Well Report)

16-inch surface casing to 1,136 ft MD; 11-3/4-inch 65 lb/ft intermediate
casing to 4,837 ft MD (14-3/4-inch hole to 4,845 ft); 9-1/2-inch open hole
below, to TD at 10,947 ft MD. All four trial runs used 9-1/2-inch PDC bits.

## Notes on the trial report

- Eavor's history match "closely aligned with the measured values" for both
  insulated and non-insulated runs. Its sentence about "drastically
  over-estimating the observed MWD temperatures when using non-insulated
  regular drill pipe" describes a counterfactual: BHA 11 and 12 re-simulated
  as if run on regular pipe, to show the IDP's effect. It is not a failure of
  their model on the runs without IDP.
- Eavor's model is described in "Enablement of High-Temperature Well Drilling
  for Multilateral Closed-Loop Geothermal Systems" (Stanford Geothermal
  Workshop, 2023).
- The "maximum formation temperature 450 °F" sometimes quoted for this well is
  not stated in the trial report or the End of Well Report; the 16A survey
  gives about 426 °F at 16B's TD.
