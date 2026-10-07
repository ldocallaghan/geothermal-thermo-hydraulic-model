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
`data/forge16b/raw/` and run `python data/forge16b/build_extracts.py`:

| File (GDR) | SHA-256 (first 16) |
|---|---|
| Eavor FORGE Insulated Drill Pipe Trial Results - 23 June 2023 (003).pdf (1516) | b392f98198f3b661 |
| 16B(78)-32 Well Survey.zip (1516) | 38c59c9452da709f |
| 16B Daily Reports.zip (1516) | 2c2475d95be35a1f |
| 16B mud temp logs.zip (1516) | ecefb1fb69ecda06 |
| 16B_Pason.zip (1516) | 858a6eaf56c2cb3c |
| End of Well Report-16B78-32-May_2024.pdf (1516) | b35f25eda60fcec9 |
| Utah FORGE Deep Well Temperature Profiles_Sept 2022.xlsx (1421) | c1a35f7914b5f591 |

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
