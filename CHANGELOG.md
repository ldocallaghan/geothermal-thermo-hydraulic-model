# Changelog

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
