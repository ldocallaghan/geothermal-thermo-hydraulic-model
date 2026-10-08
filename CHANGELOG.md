# Changelog

## v1.2.2: the verdicts hold while circulating (October 2026)

Corrections from a review of v1.2.1. No verdict changes.

- **Trips.** v1.2.1 checked the breakout bound at static fluid weight, citing
  connections and trips, but with the circulating wall temperature. That
  pairing holds for a connection, minutes with the wall still cool, but not
  for a trip: at 11 to 13 km a bit change takes the bottom of the hole out of
  circulation for a day or more, and the wall loses most of its cooling. The
  README and the code now say that the cooled verdicts hold while circulating
  and that the uncooled sensitivity is the guide to a trip. Modelling the
  reheating over the drilling cycle, with the mitigations it allows, is left
  for v1.3.
- **Flow for GO** now requires every case to be within the frictional cap,
  not only the worst-verdict case. No result changes.
- **README.** The opening carries the rig-capacity caveat on higher flows.
- **Fig. 6.** The site names on the right panel no longer run together.

## v1.2.1: commercial pipe as the base case (October 2026)

Corrections from a review of v1.2.

- **The site verdicts now use pipe that can be bought.** v1.2 tried pipe types
  from best insulated to worst and picked the dual-wall pipe, which exists
  only in a modelling study. The base case is now the best-insulated pipe that
  is commercially available, TK-Drakōn coated pipe; dual-wall pipe and
  vacuum-insulated tubing are sensitivities.
- **Flow for GO.** For each site the table also gives the lowest flow at which
  the same pipe gives GO within the pump limit, with its standpipe pressure,
  flagged where it exceeds the 700 gal/min at which FORGE 16B was drilled.
- **Circulating pressure on the fracture bound.** The annular friction loss,
  as an equivalent density (0.02 to 0.04 SG at the sites), is added to the
  fluid weight on the upper bound of the window. The breakout bound stays at
  static fluid weight, which the wall sees whenever circulation stops. The
  static check on both bounds is a sensitivity.
- **Solver seeds.** A solve with a depth-varying exposure that fails from one
  starting profile is retried from others. Where it converges the answer is
  the same to within 0.01 °C; before, which FORGE samples failed, and so
  dropped out of a run mean, depended on library versions, and one slow test
  failed on another machine. `requirements-lock.txt` records the versions the
  reference results were made with.
- **README.** The models' accuracy is stated by model; Eavor's full 47 to
  75 °F estimate of the pipe's benefit is quoted; the Fig. 1 caption notes the
  TK-Drakōn case is above the tool limit at 30 kg/s; the Holtzman et al. (2023)
  reference is completed.

### Verdicts, v1.2 → v1.2.1 (target 400 °C)

| Site | v1.2 (dual-wall pipe) | v1.2.1 (TK-Drakōn) | Flow for GO |
|---|---|---|---|
| Upper Rhine / Soultz | GO | CONDITIONAL [GO..CONDITIONAL], 1.07 SG | 60 kg/s |
| United Downs | GO | CONDITIONAL [GO..CONDITIONAL], 1.07 SG | 80 kg/s |
| Larderello | GO | GO | |
| Pannonian Basin | GO | CONDITIONAL [GO..CONDITIONAL], 1.06 SG | 70 kg/s |
| Newberry | GO | GO | |

Each change follows from the pipe: TK-Drakōn needs 35 to 50 kg/s to keep the
bit below 200 °C, and the wall is then at 184 to 188 °C, against 60 to 65 °C
with dual-wall pipe. With circulating pressure on the fracture bound, United
Downs is NO-GO in its worst case with an uncooled wall or with the 63° limit.

## v1.2: circulation checked against field data (October 2026)

The circulation model (Model 1) was rebuilt on real well and string geometry,
benchmarked against a published simulator and validated against the FORGE
16B insulated drill pipe trial. The sites were rerun with it, Newberry was
added, and the README now places the work in its literature.

### Circulation model

- **Well and string geometry.** Model 1 runs along measured depth through a
  survey, casing and cement, and a string of segments (drill pipe, heavy-weight
  pipe, collars), each with its own pipe type. Rock temperature is taken at
  true vertical depth.
- **Named pipe types.** The pipe wall is a set of layers in series, with an
  optional fraction of bare tool joints. Conventional, internally and
  externally coated and dual-wall pipe follow Wu et al. (2025); NOV's
  TK-Drakōn coating (0.635 mm at 0.162 W/m·K) is the best-insulated pipe that is
  commercially available; vacuum-insulated tubing is a sensitivity. The 0.02 and 0.10 W/m·K
  walls are gone; the old 0.02 wall survives only as "legacy vacuum tubing"
  for the earlier results.
- **Exposure time by depth.** The rock term can take each depth's time since
  the bit passed, from a rate of penetration or a drilling record, at least
  one hour.
- **Hydraulics.** Friction loss in each piece of pipe and annulus, bit nozzle
  loss, buoyancy and standpipe pressure against a 51.7 MPa pump limit, and the
  lowest annular velocity against a hole-cleaning reference: the 1.245 m/s at
  which FORGE 16B was drilled and cleaned. No published minimum could be
  checked.
- **Friction heating**, found missing by the benchmark: the dual-wall pipe
  there loses about 70 MPa to friction and warms the mud by about 18 °C.
- **Numerics.** Changes of string or hole are ramped over 1 m in a single
  domain, an exposure history is smoothed, a solve can be seeded from an
  earlier one, friction heat uses liquid properties, and a solution outside the
  water table's range is reported as failed.

### Benchmark and validation

- **Wu et al. (2025)** (`benchmark_utaustin.py`): conventional and dual-wall
  pipe within 7 °C; coated pipe 30 to 60 °C hotter with their stated 1 mm
  coatings and 10 to 36 °C colder with their conductivity across the whole
  wall. Nothing fitted.
- **FORGE 16B** (`validate_forge16b.py`): an annulus film multiplier of 0.75
  fitted on the runs without insulated pipe, Eavor's pipe at an effective
  4.9 W/m·K fitted on BHA 11 (close to a 1 mm coating at 0.47 W/m·K), and BHA 12
  predicted blind to within 7 °F, or 4 to 6 °F with the film unfitted. The
  sites use the unfitted film, with 0.75 as a sensitivity.

### Sites

- **Flow and pipe.** Each site is drilled on Wu et al.'s well design. The flow
  is the larger of the survival and hole-cleaning flows within the pump limit,
  on the best-insulated pipe that allows it; this replaces the lowest
  survivable flow. Hole cleaning sets 30 kg/s everywhere, on dual-wall pipe.
- **Newberry** added, for 400 °C at 3.6 km, with the NWG 55-29 cross-check
  (`crosscheck_newberry.py`). The volcanic breakout width is matched at 25 to
  35 K of wall cooling, but the model also breaks out the granodiorite, where
  none were logged, so Newberry is speculative. Its friction coefficient
  (0.55, and 0.70) is the one its stresses were derived with; stress cases can
  now carry their own.
- **Commercial-pipe sensitivity.** The dual-wall pipe exists only in a
  modelling study, so each site also reports the verdict with TK-Drakōn.

### Verdicts, v1.1.1 → v1.2 (target 400 °C)

| Site | v1.1.1 | v1.2, dual-wall pipe | v1.2, commercial pipe | Why |
|---|---|---|---|---|
| Upper Rhine / Soultz | CONDITIONAL [GO..CONDITIONAL] | GO | CONDITIONAL [GO..CONDITIONAL] | wall 148 → 64 °C: hole cleaning sets 30 kg/s, on dual-wall pipe |
| United Downs | CONDITIONAL [GO..CONDITIONAL] | GO | CONDITIONAL [GO..CONDITIONAL] | wall 165 → 65 °C, the same reason |
| Larderello | GO | GO | GO | wall 55 → 48 °C |
| Pannonian Basin | CONDITIONAL [GO..CONDITIONAL] | GO | CONDITIONAL [GO..CONDITIONAL] | wall 139 → 62 °C, the same reason |
| Newberry | (new) | GO | GO | wall 49 °C; speculative |

The uncooled and 63° sensitivities don't depend on the circulation and are
unchanged: both well-characterised sites are NO-GO with an uncooled wall and
the 63° limit together.

### README

The introduction treats the loop as a way of drilling, with its heat a
by-product: a single closed loop is conduction-limited to a few MW, and the
case for superhot rock assumes open-loop flow. A related-work section, the
benchmark and validation, the Newberry cross-check and the commercial-pipe
results are new; the figures are regenerated; the data credits cover the
FORGE and NWG 55-29 datasets.

## v1.1.1: how far the verdicts depend on wall cooling (October 2026)

Corrections from a review of v1.1. The verdicts themselves are unchanged;
what changed is how honestly the output and README describe them.

- **The "no wall cooling" sensitivity was not uncooled.** It removed the
  thermal hoop stress but left the wall at the cooled temperature, so the
  rock kept the strength gained from cooling. It now puts the wall at rock
  temperature. Uncooled, no case is GO at either well-characterised site:
  Soultz needs at least 1.33 SG against a fracture limit of 1.37, United
  Downs 1.29 against 1.32.
- **New sensitivity: an uncooled wall with UD-1's 63° breakout limit.**
  Both well-characterised sites are NO-GO under it (Soultz 1.58 against
  1.37 SG, United Downs 1.54 against 1.32), and Pannonian at worst.
- **README.** The summary and the site assessment now say that the open
  mud-weight window at both sites comes from cooling the wall, with the
  uncooled and combined numbers. The limitations add that the cooling at the
  targets (235–252 K) is about six times the 0–40 K the UD-1 calibration
  covered, that applying the full thermal stress at the wall surface makes the
  cooled results an upper bound, and that the calibrated rock strength holds
  only for the stresses assumed with it. United Downs is described as measured
  or calibrated, not measured. UD-1's deviated section, drilled through
  breakouts up to 90°, is noted as consistent with the 90° criterion.
- **Pannonian.** The full strike-slip range spans GO to CONDITIONAL with the
  wall cooled and CONDITIONAL to NO-GO without, not "GO to NO-GO" as v1.1's
  notes said (that dated from an earlier rule). The equal-stress worst case is
  now described as the wall yielding all round, not as a 180° breakout.
- **Tests.** The water-property table is now loaded on first use, so the fast
  test suite no longer builds it on a fresh clone (seconds, not minutes).

## v1.1: stability calibrated against UD-1 (October 2026)

Models 1 to 4 are unchanged. What changed is the hole-stability calculation,
the site data it runs on, and how results are presented.
`evaluate(site, v10=True)` still reproduces the v1.0 table exactly, and the
test suite holds it to that (`tests/golden/v10_site_table.json`).

### Why

v1.0 ruled out United Downs as NO-GO on stability. UD-1, drilled there to
5,275 m MD (5,058 m TVD) with breakouts, shows the stability calculation was
the problem: v1.0's stress inputs were not the site's, pore pressure never
entered the breakout check, the verdict demanded breakouts be suppressed
entirely, and the cooled wall's thermal stress was ignored.

### Verdicts, v1.0 → v1.1 (target 400 °C)

| Site | v1.0 | v1.1 | Tier |
|---|---|---|---|
| Upper Rhine / Soultz | CONDITIONAL | CONDITIONAL [GO..CONDITIONAL], +0.5 MPa | evidence-based |
| United Downs | NO-GO (stability) | CONDITIONAL [GO..CONDITIONAL], +2.9 MPa | evidence-based |
| Larderello | GO | GO | speculative |
| Pannonian Basin | GO | CONDITIONAL [GO..CONDITIONAL], +3.1 MPa | speculative |

The README gives the reason for each change.

### Stability physics

- **Stress profiles with depth.** Each site supplies Sv, Shmin and SHmax as
  gradient plus intercept, a pore-pressure profile, its source, and the depth
  range the data cover. Results report how far the target lies below the data.
- **Frictional admissibility.** Effective S1/S3 is checked against
  ((1+μ²)^½+μ)², with μ = 0.8, the value the United Downs SHmax was derived
  with (0.6 and 0.85 as sensitivities). Below a profile's data, SHmax is capped
  at the frictional limit. An inadmissible stress state can't be GO.
- **Effective-stress breakout check.** Kirsch hoop stress, Mohr-Coulomb in
  effective stress, no filter-cake credit beyond Pw − Pp.
- **Per-site rock strength**, in place of a global 200 MPa.
- **Breakout width**, in closed form.
- **Mud window and width verdict.** GO if the breakout is within 90° at
  hydrostatic mud (Zoback 2007), CONDITIONAL if heavier mud below Shmin (less
  0.05 SG) gets it there, NO-GO if not. Tensile fractures are reported, with
  zero wall tensile strength, but don't bound the window: UD-1 logged them and
  still reached TD.
- **Thermal hoop stress** at the cooled wall, with the wall temperature from
  Model 1's circulating bottom-hole temperature. Strength, thermal stress,
  isotropic closure and Model 4's creep closure all use that one wall
  temperature; v1.0's 200 °C assumption is gone.
- **Cases.** A site can carry several stress cases (a measured range, or regime
  bounds) and several strength cases; every pair is run, the table row is the
  worst (worst verdict, then narrowest mud window), and the verdict shows the
  range.

### Calibration and checks

- **UD-1.** `src/calibrate_ud1.py` fits rock strength to the BGS image log of
  UD-1 (24 breakouts in the 12.25" section), with wall cooling while drilling
  bracketed at 0–40 K. Weak zones 118–147 MPa (median), intact rock at least
  172–203 MPa; United Downs carries these as four strength cases. The widest
  breakout in a trouble-free section, 63°, is reported as a site-calibrated
  sensitivity on the 90° limit.
- **Soultz.** `src/crosscheck_soultz.py`: breakouts at 5 km with an open window,
  and a breakout onset depth consistent with GPK4, with nothing fitted.

### Data

- **Data basis.** Every site input states its basis and source. Sites without
  measured or calibrated stress and strength, checked against a well, are
  reported separately as speculative: Larderello and Pannonian.
- **United Downs stresses** now follow Reinecker et al. (2021): Sv 25.275 MPa/km,
  Shmin 13.21 MPa/km + 3 MPa, SHmax 25.99 MPa/km + 5.9 MPa, pore pressure
  9.494 MPa/km below a 61 m fluid level. SHmax/Shmin falls from 2.55 to 1.97.
  v1.0's ratios (Shmin 0.55 Sv, SHmax 1.40 Sv) were attributed in the code to
  Pine & Batchelor and the Rosemanowes experiment, but were never traced to a
  value in a source. The published profile replaces them because it is the
  site's own, checked against the paper.
- **Soultz stresses** follow Valley & Evans (2007), valid 1.5–5.0 km, with
  SHmax as a measured range (0.90–1.05 Sv). **Soultz strength** is their lab UCS
  of unaltered granite, 100–130 MPa, replacing an unsourced 170 MPa.
- **Larderello and Pannonian** have no published stress magnitudes or rock
  strengths that we found; their stresses are bounded by the faulting regime
  (Liotta & Brogi; Bada et al. 2007). Pannonian's interior is strike-slip,
  locally transtensional today; v1.0 assumed extension.
- **Temperatures.** United Downs follows Reinecker et al.'s ~180 °C at 5 km
  (v1.0: 190 °C), which moves 400 °C from 12.5 to 12.9 km. Larderello follows
  the Venelle-2 well (507–517 °C at 2.9 km; Bertani et al. 2018), which puts
  400 °C at 2.4 km, where v1.0 had 2.85 km.
- Sources and transcriptions: `data/ud1/`, `data/sites/stress_sources.md`.

### Also

- A pytest suite (`tests/`), with the v1.0 anchors, a frozen v1.0 table, and
  tests for every v1.1 change.
- Figures regenerated: `comparative_sites.png` (breakout width and data depth by
  site), `site_dashboard.png`; new `ud1_calibration.png`. Model 1–3 figures are
  unchanged because those models are.

## v1.0

Four European sites scored for a 400 °C target with Models 1 to 5.
