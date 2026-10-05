# Site data sources

Where each site's stress cases, rock strength, pore pressure and temperature in
the code come from, and what each rests on (the `data_basis` of each
`SiteProfile`). United Downs' main source has its own file,
`data/ud1/reinecker2021.md`.

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

Wellbore observations used by the Soultz cross-check (`src/crosscheck_soultz.py`),
same paper, GPK3 and GPK4 ultrasonic (UBI) logs run 12-15 hours after drilling:

- About 10% of the logged length broken out in each well; breakouts at 5 km.
- GPK4: no high-confidence breakouts above 3,000 m TVD, sparse to 3,670 m TVD,
  dense below. Reproducing the 3,670 m onset needs SHmax ≥ 0.9 Sv.
- Tensile fractures almost continuous to 2,180 m TVD in GPK4, sporadic below.
- Wall thermal stress at logging time −17.1 MPa (3,160 m TVD) to −31.3 MPa
  (2,235 m TVD), from MWD bit temperatures and the logging tool's temperature.
- UCS of ten samples of unaltered Soultz granite: **100–130 MPa**. The model
  uses these as Soultz's strength cases, in place of v1.0's unsourced 170 MPa.
- Annulus pressure kept near hydrostatic while drilling.

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

**Temperature.** Bertani, R. et al. (2018). "The first results of the
DESCRAMBLE project." *Proc. 43rd Workshop on Geothermal Reservoir
Engineering*, Stanford, SGP-TR-213. Venelle-2 was deepened from 2.2 km
(350 °C) to 2.9 km: static logs give ≥ 504 °C at 2,815 m and 507–517 °C at
~2.9 km (not fully re-equilibrated, so lower bounds). v1.0's profile read
400 °C at 2.85 km, over 100 K too cold there. v1.1 follows Venelle-2:
15 °C + 152 °C/km to 2.2 km, 250 °C/km to 2,815 m, 100 °C/km below, which puts
400 °C at ~2.4 km.

Drilling record, same paper: leak-off tests at 2,500, 2,585 and 2,616 m (no
values given), total circulation loss at 2,334 m, two differentially stuck
pipes at 2,695 and 2,709 m with 1.35-1.5 SG mud, then water to TD. No
breakout or image-log data are published, so Venelle-2 can't check the
stability model. A search for Larderello rock strength found only an
analogue study of micaschist exposed on Elba (CNR).

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

Stress magnitudes: none found. The newest compilation, Békési, Porkoláb &
Wéber (2023), "Stress field of the Pannonian region", *Földtani Közlöny*
153(4), DOI 10.23928/foldt.kozl.2023.153.4.mapB, maps SHmax orientations only.

Pore pressure: the deep regime of the Great Hungarian Plain is reported as
overpressured by 1–35 MPa above hydrostatic (Almási 2001, PhD thesis,
University of Alberta, seen only through a search summary; **not yet
verified**). The model keeps hydrostatic Pp and labels it "assumed". Those
data are from the Neogene basin fill; the 9.7 km target is in basement,
where nothing is known.

Rock strength and temperature: no basement strength data found; the
45–50 °C/km gradient is regional and is extrapolated to 9.7 km with no deep
well cited.

Bounded as transtensional, matching the description of the interior: SHmax =
Sv, with Shmin from the frictional floor to Sv. The wider strike-slip range,
with SHmax up to the friction cap, gives verdicts at 9.7 km from GO to
CONDITIONAL with the wall cooled, and from CONDITIONAL to NO-GO without
cooling.

## Regime-bound construction

`transitional_bounds()` in `site_evaluation.py`: SHmax = Sv, and Shmin at 0, 25,
50, 75 and 100% of the way from the frictional floor Pp + (Sv − Pp)/R(μ) to Sv,
with hydrostatic Pp and μ = 0.8. Every case is admissible by construction. The
verdict reports the range across the cases; the table row's numbers come from
the case with the worst verdict, then the narrowest mud window.
