# Changelog

## v1.3.2: the margins on the window (October 2026)

Additions from a review of v1.3.1. No verdict changes.

- **Shmin margin.** For each site, Shmin is reduced in every stress case, with
  SHmax held and with SHmax scaled with it, until the window over the cycle
  closes. The frictional cap on SHmax below the data couples the two: at the
  friction coefficient of 0.8 the worst cases at Soultz and United Downs lie
  within a few MPa of the cap, so a lower Shmin pulls SHmax down and the
  window stays open (none within 20% at Soultz; 19% at United Downs, in a
  stress state needing friction 1.55). With the cap at 1.0, SHmax can stay on
  its trend, and the window closes for a fall of 3.1% at Soultz and 2.2% at
  United Downs (friction at least 0.87 and 0.84), or 6.6% and 4.6% with SHmax
  scaled. The review's estimate of about 2% assumed SHmax held without the
  cap. Both sets are reported, with the friction each closure needs. The
  response is not monotonic, because at large reductions the cap binds even
  at 1.0, so the reductions are scanned and the first closure refined.
- **Flow ceiling.** For each site, the highest flow from the one used at which
  the verdict over the cycle holds: 70 kg/s at Soultz and 60 kg/s at United
  Downs, both within three 14-P-220 pumps; the pump limit elsewhere. The
  circulating table no longer gives a flow for GO, which over the cycle read
  as a target past the ceiling.
- **Drilling-state label.** When a connection sets a fluid heavier than water,
  the drilling state of that case is CONDITIONAL. The site verdict is
  unaffected.
- **README.** The opening states the margin and points to the trip fluid up
  the open hole (Fig. 8); the cycle table gives the flow range; Section 5
  discusses the Shmin margin.
- **The FORGE trip test.** The reviewer's environment runs Python 3.11
  against 3.14 in the lock file, which accounts for the difference in the
  no-node fit factor.

## v1.3.1: the wall is held by fluid weight (October 2026)

Corrections from a review of v1.3. No base verdict changes; the README is
reframed around the v1.3 result.

- **Connection on freshly drilled rock.** The rock the bit has just drilled is
  circulated against only until the pumps stop, then sits through the
  connection at static pressure. From the 10-second Pason record over the
  whole well (303 connections), that circulation is 20 s at the 10th
  percentile and 1 min at the median; the base case uses the 10th percentile,
  the median is a sensitivity. The connection durations now come from the same
  record: 2.8 min at the median and 6.3 min at the 90th percentile, against
  2.0 and 5.0 min from the 1-minute record over the trial week. The review
  found United Downs without a window here; that comparison added the
  circulating friction to the requirement and also took it off the limit. With
  it counted once, United Downs needs 1.26 SG against 1.29.
- **Drilling fluid against the limit while circulating.** When a connection
  sets the drilling fluid, the fluid is now also checked against the fracture
  limit less the circulating friction. At the highest flow within the pump
  limit, both well-characterised sites become NO-GO in their worst case, so a
  higher flow narrows the window over the cycle.
- **The one-hour formation allowance** now applies to the fresh-rock
  connection as well as to drilling, since both fall within that hour.
- **Full stress field in the first minutes.** The site output reports the
  failed depth from the full stress field after 0 to 10 minutes of cooling
  against the uncooled depth: the cooled skin pushes it from 1.9 to 2.1 cm at
  Soultz after a minute, and below the uncooled depth within 5 minutes.
- **When breakouts form.** Moore et al. (2011) imaged breakouts at the bit
  minutes after drilling, widening over 30 minutes to 3 days; Wenning et al.
  (2017) logged breakout growth over more than a year in crystalline rock in
  the COSC-1 borehole. The README cites both: the one-hour allowance is a
  sensitivity the data do not support, and later growth is unmodelled.
- **README.** Section 5 leads with the verdicts over the cycle; the
  assessment while circulating is kept for comparison. The opening, Section 1
  and Section 3.2 say that the wall is held by fluid weight and that cooling
  serves the tools and slows creep. The comparison with FORGE 16B's trips is
  described as too coarse to confirm or reject the reheating model. The narrow
  windows are stated. Fig. 8 adds the trip fluid needed up the open hole.
- **The FORGE trip test.** The review reported a no-node fit factor of 0.391
  against 0.357 pinned, with numpy 2.4.4 and scipy 1.17.1. Those versions
  reproduce the pinned values exactly here, so the difference has another
  cause, and the tolerance is unchanged.

### Verdicts, v1.3 → v1.3.1 (target 400 °C)

| Site | Verdict | Drilling fluid (limit) | Trip fluid (limit) |
|---|---|---|---|
| Upper Rhine / Soultz | CONDITIONAL, unchanged | 1.31 SG (1.34), unchanged | 1.30 SG (1.37) |
| United Downs | CONDITIONAL, unchanged | 1.25 → 1.26 SG (1.29) | 1.27 SG (1.32) |
| Larderello | GO | water | water |
| Pannonian Basin | CONDITIONAL, unchanged | 1.42 SG (2.48) | 1.39 SG (2.50) |
| Newberry | GO | water | water |

## v1.3: the wall through the drilling cycle (October 2026)

The verdicts now follow the borehole wall through a bit run and the trip that
ends it, with the rock behind the wall resolved in temperature and stress.

### Wall temperature through the cycle

- **Transient wall model** (`wall_thermal.py`). Implicit finite volumes for
  radial conduction on a log-spaced grid. While circulating, the wall exchanges
  heat with the annulus at Model 1's temperature and film coefficient. With
  the pumps off, the fluid in the hole is one well-mixed volume coupled to the
  wall by the laminar film (Nu 4.36); the case with no fluid is a sensitivity.
  The solver is checked against the constant-temperature cylinder (Carslaw and
  Jaeger, 1959) to within 1.2% and the recovery after constant heat flow to
  within 2%, by numerical Laplace inversion (Abate and Valkó, 2004), and for
  grid convergence and energy conservation. The line source (Bullard, 1947)
  differs from the cylinder by 6 to 13% over the times of interest and is
  reported for comparison only.
- **Cycle timings from FORGE 16B** (`drilling_cycle.py`), from the Pason
  1-minute record and the daily reports: trips at 1,694 ft/h out and 1,833 ft/h
  in, 3.9 h of routine surface time, connections of 2.0 min (median) and
  5.0 min (90th percentile), 94 ft stands. The first trip in the record began
  on 21 May.
- **Check against FORGE 16B's trips** (`validate_forge_pauses.py`). The return
  temperature over the first bottoms-up after each of the five trips is
  compared with a plug-flow column of the modelled fluid temperatures, scaled
  by one factor fitted on four trips and used to predict the fifth in turn:
  0.55, with a misfit of 13 °F against anomalies of −10 to +38 °F (16 °F with
  no fluid in the hole). The comparison was first planned at the bottom of the
  hole alone, which the returns cannot resolve, and on the pauses with the bit
  on bottom, most of which were surveys or reduced-flow circulation; both were
  dropped. Five trips through a plug-flow column are weak evidence.

### Failure behind the wall

- **Stresses in the rock** (`model5_convergence_confinement.py`): the Kirsch
  field with wellbore pressure (Jaeger et al., 2007) and the plane-strain
  thermal stresses of the temperature field (Timoshenko and Goodier, 1970),
  with Mohr-Coulomb at the local temperature.
- **Width judged at the wall, with the temperature of the rock behind it.**
  In elastic Mohr-Coulomb the failed region reaches wider just behind the wall
  than at the wall, so the extent of failure in the rock is not a breakout
  width. The width is judged at the wall, with the temperature of the rock at
  the depth the uncooled breakout reaches (about 3 cm at Soultz). This replaces
  the planned measure, the extent of failure at any radius.
- **Recalibration.** UD-1's weak-zone strength becomes 119 to 147 MPa (from
  118 to 147), with the intact bound unchanged. Soultz widths at 5 km rise by
  up to 3°, with every window still open. Newberry's volcanics are matched at
  25 to 40 K of cooling; the granodiorite is still not reproduced.

### Verdicts over the cycle

- **Drilling, connection and trip.** Each case of stress and strength is
  checked while drilling (circulating pressure), after a 90th-percentile
  connection and at the end of a trip (static pressure). Drilling is judged as
  the bit exposes the rock: elastic failure is immediate, and the rock at the
  breakout's depth is then still at formation temperature. The drilling fluid
  is the lightest static weight that holds the wall while drilling and through
  a connection; the safe pause is how long the wall holds at that weight with
  the pumps off; the trip fluid is what the wall needs at the end of the trip.
  A one-hour allowance for breakouts to form is a sensitivity.
- **Sensitivities:** half and double trip time, no fluid in the hole, the
  dual-wall pipe, the highest flow within the pump limit, and the trip
  verdict up the open hole.
- **Rig capacity.** The flow for GO is compared with three pumps of the rating
  of NOV's 14-P-220 triplex pump (1,980 hydraulic hp, up to 1,215 gal/min); all
  three sites' flows are within it.

### Verdicts, v1.2.2 → v1.3 (target 400 °C)

| Site | v1.2.2, while circulating | v1.3, over the cycle | Drilling fluid | Trip fluid (limit) |
|---|---|---|---|---|
| Upper Rhine / Soultz | CONDITIONAL, 1.07 SG | CONDITIONAL | 1.31 SG | 1.30 SG (1.37) |
| United Downs | CONDITIONAL, 1.07 SG | CONDITIONAL | 1.25 SG | 1.27 SG (1.32) |
| Larderello | GO | GO | water | water (1.33) |
| Pannonian Basin | CONDITIONAL, 1.06 SG | CONDITIONAL | 1.42 SG | 1.39 SG (2.50) |
| Newberry | GO | GO | water | water (1.36) |

No verdict changes, but the weights do. The drilling fluid rises because the
wall is judged as the bit exposes it, before circulation has cooled the rock
the breakout grows into, so it is close to the uncooled weight; with a
one-hour allowance it would be 1.14 SG at Soultz and 1.12 SG at United Downs.
At Soultz that weight holds the hole through a 44-hour trip; at United Downs
it holds for 23 hours, and the fluid must be raised to 1.27 SG for a 52-hour
trip. The v1.2.2 table, while circulating, is kept unchanged as a column of
the site table.

### Other

- New figures show the trip check (Fig. 3), the failure margin around the hole
  at Soultz (Fig. 6) and the wall at Soultz through a bit run and a trip
  (Fig. 8). The Soultz dashboard shows the fluid weights through the cycle.
- README: the reheating check, the stresses behind the wall, the cycle table,
  and limitations for thermo-poroelastic effects (Ghassemi et al., 2009), surge
  and swab, bit life, the single source of cycle timings and the time
  breakouts take to form.
- The v1.2.2 site table is pinned with the circulation model of v1.2; the
  v1.3 table is pinned alongside it.

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
