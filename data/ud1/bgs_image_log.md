# UD-1 borehole image interpretation (BGS)

`UD1_Borehole_Imaging_Interpretation.csv` is the British Geological Survey's
interpretation of the UD-1 acoustic image log, unmodified.

- **Dataset:** "United Downs 1 borehole image interpretation", BGS, 10 July 2024.
  NGDC accession item 185507 (BGSID 21246966). DOI
  [10.5285/a040a246-8691-428b-b146-13874e743b0a](https://doi.org/10.5285/a040a246-8691-428b-b146-13874e743b0a).
  Listed on data.gov.uk as `united-downs-1-borehole-image-interpretation`.
- **Downloaded:** 3 October 2026 from
  `https://webservices.bgs.ac.uk/accessions/download/185507?fileName=UD1_Borehole_Imaging_Interpretation.csv`.
  SHA-256 `4508b9a5fadb7473714eab87d72611bedba6f0b46af837567899af8e942e258c`.
- **Licence:** Open Government Licence. Acknowledgement, as the two catalogue
  entries give it: "Contains NERC materials ©NERC 2024" (data.gov.uk) and
  "Contains data supplied by UKRI." (NGDC). Work by BGS under the GWatt
  project (NE/S004262/1), funded by UKRI NERC. Some underlying material, such
  as the image DLIS files, belongs to Geothermal Engineering Ltd.
- **Not copied:** the dataset's `readme.rtf` (33 MB, mostly an embedded
  figure). Its text is summarised below.

## Contents

12,031 features from 907 to 5,158 m MD (two header rows: names, then units).
Columns: Measured Depth (m), Hole Azimuth and Hole Deviation (deg), Average
Caliper (in), Dip Classification, True Dip Azimuth and True Dip Inclination
(deg), Breakout_Height_N (m) and Breakout_Width_N (deg); -9999 marks no value.

| Classification | Count |
|---|---|
| Conductive fracture | 8,600 |
| Resistive fracture | 3,353 |
| Breakout | 75 |
| Induced fracture | 2 |
| Fault | 1 |

For a breakout, the depth is the midpoint, the azimuth is the breakout
direction, the height is its length along the hole, and the width is its
angular width in degrees. Methods per the readme: planes were picked one by
one on amplitude and travel-time images. Fractures with wall damage on the
travel-time image are "conductive"; the rest are "resistive". Breakouts are
interpreted as in Kingdon et al. (2016), *Mar. Petrol. Geol.* 73, 1–20.
Washouts hide features, so faults (and presumably breakouts inside washouts)
are undercounted.

## Against Reinecker et al. (2021)

| | Reinecker et al. (2021) | This dataset |
|---|---|---|
| Breakouts, 900–4,000 m MD | 27, 139 m total | 24, 46.0 m total, between 2,210 and 3,962 m |
| Breakout widths there | not given | 30–63°, median 52° |
| Breakouts below 4,000 m MD | "extensive", excluded (inclination > 15°) | 51, 140.5 m total, widths 35–90°, inclination 14–36° |
| Drilling-induced tensile fractures | 13, 168 m, c. 2,500–3,700 m MD | 2 "induced fractures", at 2,665 and 3,667 m MD |
| SHmax azimuth | 134° ± 25° (breakouts) | breakout azimuths mostly 215–250°, i.e. SHmax ~125–160° |

The two interpretations agree on where the wall failed and on the stress
orientation, but not on the counts. The BGS picks fewer breakouts, with much
less total length, in the near-vertical section, and only 2 tensile
fractures, at the two ends of the paper's tensile-fracture interval. Most
likely the operator's interpretation (behind the paper) and the BGS one used
different picking criteria; tensile fractures are often axial traces that a
plane-picking workflow would not record. This dataset is the one with depths
and widths, so it is what the strength calibration (src/calibrate_ud1.py)
fits, and the paper's counts are a second opinion.

The deviated 8.5" section is not used in the calibration, because the
vertical-hole Kirsch solution does not apply at 14–36° inclination. It is
still relevant to the breakout limit used in the verdicts: UD-1 was drilled
to its target through breakouts of 35–90° there (its mud losses came from an
open fracture zone at 4,890 m MD, not from the wall). That is consistent with
the 90° criterion of Zoback (2007), and suggests that 63°, the widest breakout
in the near-vertical section, is a conservative limit. The widths are not
directly comparable with a vertical hole, so this supports the criterion
without calibrating it.

Note: 140.5 m below 4,000 m MD is close to the paper's 139 m above it. That
could be coincidence, or a sign the paper's total covers a different interval
than its text says. Not resolvable from these two sources.
