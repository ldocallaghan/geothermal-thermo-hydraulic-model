# Site stress profiles (v1.1, C1 and D5)

Where each site's `stress_cases` in the code come from. United Downs has its own
file, `data/ud1/reinecker2021.md`.

| Site | Basis | Profile | Data cover |
|---|---|---|---|
| Upper Rhine / Soultz | measured range | Valley & Evans (2007) | 1.5–5.0 km |
| United Downs | measured | Reinecker et al. (2021) | to 2.0 km (Shmin) |
| Larderello | regime bounds | normal/strike-slip transition, SHmax ≈ Sv | none |
| Pannonian | regime bounds | strike-slip, locally transtensional, SHmax ≈ Sv | none |

## Upper Rhine Graben / Soultz-sous-Forêts

Valley, B. & Evans, K. F. (2007). "Stress state at Soultz-sous-Forêts to 5 km
depth from wellbore failure and hydraulic observations." *Proceedings,
32nd Workshop on Geothermal Reservoir Engineering*, Stanford, SGP-TR-183.
Read from the paper directly. Conclusions, verbatim, "valid between depths of
1.5-5.0 km (excluding the effects of local stress heterogeneity)":

- Sv = −1.30 + 25.50 z[km] MPa (±2 MPa)
- Shmin = −1.78 + 14.06 z[km] MPa, a fit to casing-shoe pressures in the large
  GPK1, GPK3 and GPK4 injections (GPK2 excluded as about 4 MPa low)
- −1.17 + 22.95 z[km] ≤ SHmax < −1.37 + 26.78 z[km] MPa, i.e. 0.90 Sv ≤ SHmax ≤
  1.05 Sv, from breakouts (lower bound) and the absence of tensile fractures in
  the lower half of GPK4 (upper bound)
- The upper bound from frictional strength alone, at μ = 1.0, is 1.21 Sv
  (Evans 2005), so the wellbore-failure bound is the tighter one.

The code carries three cases: the SHmax lower bound, its midpoint, and its upper
bound. Pore pressure is hydrostatic (the paper treats the wells as near-
hydrostatic). Shmin/Sv = 0.54 at 5 km, matching v1.0's ratio.

## United Downs

Data-coverage depth: Reinecker et al. (section 6.3) fit Shmin to refracture
pressures from the Rosemanowes hydrofracture tests of Pine, Ledingham &
Merrifield (1983b), "In-situ stress measurement in the Carnmenellis granite –
II. Hydrofracture tests at Rosemanowes quarry to depths of 2000 m", *Int. J.
Rock Mech. Min. Sci. Geomech. Abstr.* 20(2), 63–72. So the Shmin data end at
2.0 km. SHmax is derived from Shmin at μ = 0.8, so it covers no deeper. (The
Rosemanowes wells themselves reached 2.65 km, Reinecker et al. section 1.)
Checked against the published paper on 3 October.

## Larderello

Liotta, D. & Brogi, A. "Pliocene-Quaternary fault kinematics in the Larderello
geothermal area (Italy): insights for the interpretation of the present stress
field." Read from the manuscript as submitted to *Geothermics* (GEOT_2019_13);
**the published volume, pages and DOI are not yet confirmed.**

No stress magnitudes and no pore pressures. Regime, from fault-slip inversion
of Pliocene-Quaternary faults compared with focal mechanisms of local
earthquakes (M < 2, 2–6 km deep, bracketing the 2.85 km target):

- NW-trending faults are dominantly normal (σ1 near vertical, σ3 NE-SW).
  NE-trending faults are left-lateral oblique to strike-slip, and can also be
  reactivated as normal faults.
- Focal mechanisms are "normal to strike-slip"; Albarello et al. (2005) found a
  strike-slip component "nowadays dominant". Breakouts in the inner Northern
  Apennines give a NE-oriented Shmin with mainly normal solutions, plus a
  strike-slip component at Larderello (Mariucci et al. 1999; Montone et al.
  1999, 2004).
- The authors explain this as competition between crustal stretching and
  heat-flow-driven uplift, with the intermediate stress axis switching
  between vertical and horizontal.

σ2 switching between vertical and horizontal means Sv ≈ SHmax, which is the
`transitional_bounds()` construction.

Pore pressure: hydrostatic by default, and the paper gives no basis to change
it. It points both ways: the shallow reservoir is vapour-dominated, so likely
underpressured, while the paper reports over-pressured fluid injection and
hydrofracturing along active faults and non-double-couple focal mechanisms
attributed to fluid pressure. A sourced figure is still to collect.

## Pannonian Basin

No measured stress magnitudes found. Bada, G., Dövényi, P., Horváth, F.,
Szafián, P. & Windhoffer, G. (2007). "Present-day stress field in the Pannonian
Basin and the surrounding Alpine-Carpathian-Dinaric orogens." *Földtani Közlöny*
137(3). The basin interior (eastern Transdanubia, Great Hungarian Plain) is
strike-slip, locally transtensional. Thrusting and transpression are at the
margins. This contradicts v1.0's normal-faulting ratios (Shmin 0.60 Sv,
SHmax 0.90 Sv), which reflected the Miocene back-arc extension, no longer active.

Bounded as transtensional (decided 3 October): SHmax = Sv, with Shmin from the
frictional floor to Sv. The wider strike-slip range, with SHmax up to the
friction cap, spans GO to NO-GO at 9.7 km and was rejected as uninformative.

## Regime-bound construction

`transitional_bounds()` in `site_evaluation.py`: SHmax = Sv, and Shmin at 0, 25,
50, 75 and 100% of the way from the frictional floor Pp + (Sv − Pp)/R(μ) to Sv,
with hydrostatic Pp and μ = 0.8 (D1). Every case is admissible by construction.
Until step 5 adds range verdicts, the table row uses the case needing the most
mud (decided 3 October), which is the floor case.
