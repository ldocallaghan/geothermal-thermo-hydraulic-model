# Site stress profiles (v1.1, C1 and D5)

Where each site's `stress_cases` in the code come from. United Downs has its own
file, `data/ud1/reinecker2021.md`.

| Site | Basis | Profile | Data cover |
|---|---|---|---|
| Upper Rhine / Soultz | measured range | Valley & Evans (2007) | 1.5–5.0 km |
| United Downs | measured | Reinecker et al. (2021) | to 2.8 km (Shmin) |
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

Data-coverage depth: Reinecker et al. fit Shmin to Rosemanowes hydrofracture
tests (Pine et al. 1983). The Rosemanowes programme drilled boreholes to at most
2.8 km, so 2.8 km is an upper bound on where the Shmin data end. SHmax is
derived from Shmin at μ = 0.8, so it covers no deeper. The exact depths of the
Pine et al. (1983) tests are not confirmed; they may stop shallower.

## Larderello

No measured stress magnitudes found. Regime from structural field data and
focal mechanisms of local earthquakes in the Lago Basin: NW-striking normal
faults and NE-striking left-lateral oblique-to-strike-slip faults active
together, with the intermediate stress switching between vertical and
horizontal. That puts SHmax close to Sv. Attributed to Brogi and co-authors
(CNR / University of Bari repositories); **the full citation is not yet
confirmed**, because the repository pages refused automated access.

Pore pressure: hydrostatic by default. The field is vapour-dominated and
probably underpressured, which would raise effective stresses; a sourced
figure is still to collect.

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
