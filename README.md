# Deep and superhot-rock geothermal: a first-principles model stack and site-evaluation tool

A self-contained set of coupled physical models, in Python, for the engineering questions that decide a deep or superhot-rock geothermal well: whether it can be drilled at all, whether the tools survive down there, whether the hole stays open, and whether it returns more energy than it costs. Alongside them is a tool that scores a real candidate site against all of it.

The models are built from first principles with real material data: IAPWS-95 water properties, with rock and in-situ-stress parameters from the literature. Every constant is sourced and adjustable, and every model states its own assumptions and limits. The aim was an analysis that survives being checked, including the results that undercut the original premise.

> These are order-of-magnitude models: 1-D or axisymmetric, mostly quasi-steady, with several coupling parameters that are uncertain and flagged as such. Their job is to locate the binding constraints and put numbers on the trade-offs. They do not replace a full reservoir or geomechanical simulator. The limitations are listed per model and again at the end.

> **Version 1.1** rebuilt the hole-stability calculation and checked it against two real wells: UD-1 at United Downs, the calibration well, and the GPK wells at Soultz, a cross-check. Every site input now states what it rests on, and sites without measured stresses and rock strength are reported separately, as speculative. All the 400 °C targets lie kilometres below the deepest stress data, so every verdict is an extrapolation; see [Site results](#site-results-target-400-c) and [Limitations](#limitations-read-before-trusting-any-number). What changed, and why, is in [CHANGELOG.md](CHANGELOG.md).

---

## The questions, and what the models found

| # | Model | Question | Headline result |
|---|-------|----------|-----------------|
| 1 | `model1_coupled.py` | Can a circulating coolant loop keep the bottom-hole assembly survivable while exporting heat at 10-20 km? | Yes, within a flow and insulation window. A single closed loop is conduction-limited at around 3-4 MW. Water stays dense and supercritical, with no flashing, at anything deeper than about 2 km. |
| 2 | `model2_spallation.py` | Can a cold-coolant quench fracture hot rock at the cutting face? | Tensile spallation is confinement-limited. It works only in low-K0 (extensional) crust below the brittle-ductile transition. Elsewhere the quench does not break the rock outright, though it does lower the effective cutting energy. |
| 3 | `model3_optimiser.py` | Where does coupling drilling to thermal balance actually win? | Quench raises the rate of penetration by 1.5 to 2.3 times in the assisted regime, and up to 7.5 times in a narrow spallation window. The model outputs a siting and operating map. |
| 4 | `model4_hole_stability.py` | Does cooling the hot ductile rock keep the hole open? | Cooling slows time-dependent creep closure by a factor of 10^5 to 10^6 at a 200 °C wall; in the v1.1 site runs, with Model 1's cooler walls (55 to 165 °C), by 1.7×10^6 to 1.3×10^7. Convective heat resupply to the cooled zone is negligible, with a Péclet number well below 1. |
| 5 | `model5_convergence_confinement.py` | Is breakout at depth a collapse or a manageable yielded zone? | Breakout produces a millimetre-scale yielded zone held by fluid pressure, well short of collapse. Cooling adds 3 to 6 km of stable depth. In v1.1 the site verdicts use an effective-stress, breakout-width check with the wall's thermal stress, calibrated against UD-1 (below). |

### The equations

The core of each model, in the notation the code uses.

**Model 1** solves two coupled advection balances for the down and up legs, with the rock wall as a fixed-temperature boundary:

$$\dot{m}c_p\frac{dT_d}{dz}=UA_i(T_u-T_d)$$

$$\dot{m}c_p\frac{dT_u}{dz}=UA_i(T_u-T_d)-UA_o\big(T_{\text{rock}}(z)-T_u\big)$$

with boundary conditions at the surface and at the bit:

$$T_d(0)=T_{\text{inj}},\qquad T_u(L)=T_d(L)+\frac{Q_{\text{face}}}{\dot{m}c_p}$$

Output per unit length is set by conduction into the rock, which is why diameter barely helps:

$$q'(z)=\frac{2\pi k_{\text{rock}}\Delta T}{\ln(r_\infty/r_w)},\qquad r_\infty\approx 2\sqrt{\alpha t}$$

**Model 2** puts the quenched face under constrained thermal stress and asks whether it beats the confining stress plus the rock's tensile strength:

$$\sigma_T=\frac{E\alpha\Delta T}{1-\nu},\qquad \sigma_T-K_0\rho g z>T_0$$

The transient temperature field is the convective half-space solution:

$$T(x,t)=T_{\text{rock}}+(T_{\text{cold}}-T_{\text{rock}})\big[\mathrm{erfc}(\eta)-e^{-\eta^2}\mathrm{erfcx}(\eta+\beta)\big],\qquad \eta=\frac{x}{2\sqrt{\alpha t}},\quad \beta=\frac{h\sqrt{\alpha t}}{k}$$

**Model 3** turns the quench into an effective cutting energy and a penetration rate, then feeds the face heat back into Model 1:

$$\text{MSE}_{\text{eff}}=\text{MSE}(1-D),\qquad \text{ROP}=\frac{P_{\text{mech}}}{\text{MSE}_{\text{eff}}\cdot A_{\text{bit}}}$$

$$Q_{\text{face}}=\underbrace{\text{MSE}\cdot A\cdot\text{ROP}}_{\text{cutting}}+\underbrace{\rho c\cdot A\cdot\text{ROP}\cdot\Delta T}_{\text{cuttings sensible heat}}+Q_{\text{cond}}$$

**Model 4** rates time-dependent creep closure with an Arrhenius power law, and checks that convective resupply stays negligible:

$$\dot\varepsilon=A\sigma^{n}\exp\left(-\frac{Q}{RT}\right),\qquad \mathrm{Pe}=\frac{vL}{\alpha}\ll 1$$

**Model 5** takes the plastic zone from the ground-reaction curve:

$$k=\frac{1+\sin\phi}{1-\sin\phi},\qquad p_{\text{cr}}=\frac{2p_0-\sigma_{cm}}{1+k}$$

**Hole stability in the site verdicts (v1.1).** The hoop stress at the wall of a vertical hole, at angle θ from the SHmax azimuth, includes the thermal stress from the wall being cooled from rock temperature to the circulating fluid's:

$$\sigma_\theta(\theta)=S_H+S_h-2(S_H-S_h)\cos 2\theta-P_w+\Delta\sigma_T,\qquad \Delta\sigma_T=-\frac{E\alpha}{1-\nu}\,(T_{\text{rock}}-T_{\text{wall}})$$

The wall breaks out where Mohr-Coulomb fails in effective stress, with pore pressure P_p and no filter-cake credit beyond P_w − P_p:

$$\sigma_\theta-P_p>\sigma_{cm}(T_{\text{wall}})+k\,(P_w-P_p)$$

That inequality holds over a lobe of width W around the Shmin azimuth, in closed form:

$$W=\arccos\!\left(-\frac{a}{d}\right),\qquad a=S_H+S_h-P_w+\Delta\sigma_T-P_p-\sigma_{cm}-k(P_w-P_p),\quad d=2(S_H-S_h)$$

Wells are drilled with breakouts, so the verdict is on width: **GO** if W ≤ 90° at hydrostatic mud, **CONDITIONAL** if a heavier mud below Shmin (less a 0.05 SG margin) brings it there, **NO-GO** if none does. Every stress state is also checked for frictional admissibility (μ = 0.8):

$$\frac{S_1-P_p}{S_3-P_p}\le\left(\sqrt{1+\mu^2}+\mu\right)^2$$

### Model figures

![Model 1 coupled counterflow temperature profiles](figures/model1_profiles.png)

Model 1, G = 35 K/km with 450 °C at 12.4 km. Baseline insulation leaves the bit at 318 °C. Vacuum tubing brings it to 171 °C. Raising flow to 10 kg/s brings it to 59 °C and the export to 1.90 MW.

![Model 2 quench spallation feasibility against depth](figures/model2_feasibility.png)

Model 2, available quench tension against the tension required at each confining ratio K0. Spallation needs the red curve above the dashed line and the rock still brittle, which is the shaded window at low K0.

![Model 3 coupled optimiser siting map](figures/model3_siting_map.png)

Model 3, ROP gain from quench across target temperature and stress ratio. The 7.5x band is pure spallation, at K0 below 0.5 and below the brittle-ductile transition. Everywhere else the gain runs 1.5 to 2.3x.

### Integrating tool

- `site_evaluation.py` takes a `SiteProfile`, which holds a site's temperature profile, its stress profiles with depth (one or more cases), its rock-strength cases, and the data basis of each input. It runs Models 1 to 5 against every stress and strength case and returns a verdict with the operating envelope. The table row is the worst case; the verdict shows the range across cases.
- `comparative_sites.py` scores four European provinces, split into evidence-based and speculative tiers.
- `calibrate_ud1.py` fits rock strength to the breakouts logged in UD-1; `crosscheck_soultz.py` checks the model against the Soultz wells.

![Site evaluation dashboard for Soultz-sous-Forets](figures/site_dashboard.png)

`site_evaluation.py` output for Soultz-sous-Forêts: thermal profile at production flow, the measured stress profiles (shaded where the data are), and the scorecard. Target 400 °C at 10.7 km, 5.7 km below the stress data. Tool survival OK at 148 °C, 4.1 MW_th. Breakouts 70 to 90° wide at hydrostatic mud across the nine stress and strength cases; the worst, 90.5°, needs +0.5 MPa of mud. Verdict CONDITIONAL [GO..CONDITIONAL].

### Site results (target 400 °C)

Sites are split by what their verdicts rest on. A site is **evidence-based** only if its stress magnitudes and rock strength are measured or calibrated there, and the stability model has been checked against a real well there. Otherwise it is **speculative**: its verdict rests on assumptions, and the tables are not comparable.

**Evidence-based**

| Site | Depth to 400 °C | Verdict (range across cases) | Mud | Target below data (stress / temperature) | Rests on |
|------|-----------------|---------|-----|-----|-----|
| Upper Rhine Graben / Soultz (FR) | 10.7 km | CONDITIONAL [GO..CONDITIONAL] | +0.5 MPa (1.02 SG) | 5.7 km / 5.7 km | stresses measured to 5 km (Valley & Evans 2007); lab UCS 100–130 MPa; checked on GPK3/GPK4 |
| United Downs (UK) | 12.9 km | CONDITIONAL [GO..CONDITIONAL] | +2.9 MPa (1.04 SG) | 10.9 km / 7.8 km | stresses measured (Reinecker et al. 2021; Shmin data to 2 km); strength calibrated on UD-1 |

**Speculative**

| Site | Depth to 400 °C | Verdict | Missing |
|------|-----------------|---------|---------|
| Larderello (IT) | 2.4 km | GO | stress magnitudes (regime only), rock strength (unsourced), any wellbore-failure data. Temperature is measured (Venelle-2). |
| Pannonian Basin (HU) | 9.7 km | CONDITIONAL [GO..CONDITIONAL] | stress magnitudes (regime only), rock strength (unsourced), any well check, deep temperature data; pore pressure likely overpressured, not modelled |

Sensitivities, printed beside every verdict by `comparative_sites.py`: without the wall's thermal stress ("no wall cooling"), every site but Larderello is CONDITIONAL. At UD-1's site-calibrated breakout limit of 63° in place of 90° (the widest breakout UD-1 logged in a section drilled without trouble), the same three are CONDITIONAL.

**Against v1.0.** The v1.0 verdicts, which `evaluate(site, v10=True)` still reproduces exactly, and why each changed:

| Site | v1.0 | v1.1 | Why |
|------|------|------|-----|
| Upper Rhine / Soultz | CONDITIONAL (+3 MPa to stop breakout entirely) | CONDITIONAL [GO..CONDITIONAL] (+0.5 MPa, one case of nine) | Measured SHmax range and the lab strength (100–130 MPa, not the unsourced 170). Breakouts are now judged on width, and wall cooling offsets the effective-stress penalty. Eight of nine cases are GO; the worst sits at 90.5°. |
| United Downs | NO-GO (stability) | CONDITIONAL [GO..CONDITIONAL] (+2.9 MPa) | v1.0's stress ratios (SHmax/Shmin 2.55) were not the site's; the published profile gives 1.97. Width replaces full suppression, cooling relieves about 120 MPa of hoop stress, and the strength is calibrated on UD-1. Three of four strength cases are GO. The target is 12.9 km, not 12.5: the site is ~180 °C at 5 km, not 190. |
| Larderello | GO | GO (speculative) | Same verdict, but now flagged as resting on assumed stress and strength. 400 °C is at 2.4 km, not 2.9: Venelle-2 measured 507–517 °C at 2.9 km. |
| Pannonian Basin | GO | CONDITIONAL [GO..CONDITIONAL] (speculative) | v1.0 assumed extension; the basin interior is strike-slip, locally transtensional today (Bada et al. 2007). Bounded by the regime, the equal-stress case breaks out all round at hydrostatic mud and needs +3 MPa. |

![Breakout width against depth to 400 C for each site, and how deep each site's data reach](figures/comparative_sites.png)

Left: breakout width at hydrostatic mud across each site's cases; the verdict turns on the 90° line. Right: where each site's stress and temperature data end against its target. Filled markers are evidence-based, hollow speculative.

### Calibration against UD-1

`calibrate_ud1.py` fits rock strength to the 24 breakouts the BGS image log records in UD-1's near-vertical 12.25" section (900–4,000 m MD), with the published stresses and pore pressure, 1.05 SG mud and the wall's thermal stress. Each breakout's width fixes the strength of the rock it formed in; the rock that didn't break out bounds intact strength from below. No mud temperatures were published, so wall cooling while drilling is bracketed at 0–40 K.

| | 0 K cooling | 40 K cooling |
|---|---|---|
| Weak zones (where breakouts formed), median UCS | 147 MPa | 118 MPa |
| Intact rock (no breakout to 4 km), UCS at least | 203 MPa | 172 MPa |

The fit predicts breakouts starting at about 2.9 km TVD (logged: one at 2,210 m, the rest from 3,095 m MD), widths of 53–57° in the lowest 350 m (logged 50–53°), and a hole drillable to TD. The 51 breakouts in the deviated 8.5" section, not fitted, fall between the intact and weak-zone curves.

![Breakout width with depth in UD-1, logged and predicted, and the strength each breakout implies](figures/ud1_calibration.png)

### Cross-check against Soultz

`crosscheck_soultz.py` runs the same model on the Soultz stresses with nothing fitted. It predicts breakouts at 5 km with an open mud window in 46 of 48 strength, cooling and stress cases, as in GPK3 and GPK4. UD-1's weak-zone strengths, carried over unchanged, put the onset of breakouts at 2.6–3.8 km; GPK4's first breakouts are at about 3.0 km, dense below 3.67 km. With Soultz's own lab strength and the lower SHmax bound, the onset lands at 3.64 km, matching Valley & Evans' own analysis.

---

## Key engineering conclusions

- Closed-loop power is limited by conduction and barely moves with diameter. Ten times the borehole diameter gives about 1.9 times the heat, a logarithmic gain, while cutting cost grows with the square of the radius. Output moves with contact length, so laterals, more wells, or more open fracture area. (`model1_diameter_scaling.py`)
- Depth raises the grade of the heat more than the power, and it is capped by the brittle-ductile transition near 400 °C, which sits between 5.5 and 19 km depending on gradient. Active cooling removes tool survival as the binding constraint, and the rock becomes the constraint instead. (`model1_depth_limits.py`)
- The mechanism is prior art. Cold-quench rock fracture appears in Holtzman et al. (GRC Transactions 47, 2023) and in the cryogenic and thermal-shock EGS stimulation literature. What this repository adds is the coupled model stack: real material properties, site scoring, and known trade-offs reproduced end to end.

---

## Running it

From the repository root:

```bash
pip install -r requirements.txt
python src/water_table.py        # one-time: build the cached IAPWS-95 property table (~2.5 min)
python src/site_evaluation.py    # full verdict for the Upper Rhine Graben / Soultz site
python src/comparative_sites.py  # four-site comparison, evidence-based and speculative tiers
python src/calibrate_ud1.py      # UD-1 strength calibration and figure (no water table needed)
python src/crosscheck_soultz.py  # Soultz cross-check (no water table needed)
```

Each model file runs standalone and prints its own analysis (`python src/model1_coupled.py`, and so on). The figure scripts (`src/*_figures.py`) regenerate the PNGs into `figures/`.

Tests:

```bash
pytest                           # v1.0 anchors and v1.1 stability tests (fast, no water table)
pytest -m slow                   # full site runs and the frozen v1.0 table, need the water table
python tests/report_v10_anchors.py   # print the anchors as a table
```

### Layout
```
README.md                 this file
CHANGELOG.md              what changed in v1.1, and why
requirements.txt
data/ud1/                 Reinecker et al. (2021) transcription, BGS UD-1 image log
data/sites/               stress, strength and temperature sources for every site
figures/                  generated PNGs (embedded above)
tests/                    pytest suite; golden/ holds the frozen v1.0 table
src/
  geo_constants.py        sourced constants (geotherm, rock, drilling, geometry)
  water_props.py          IAPWS-95 wrapper (point queries)
  water_table.py          cached, vectorised IAPWS-95 property table (T,P grid)
  model1_coupled.py       coupled 1-D counterflow borehole heat exchanger (BVP)
  model1_steadystate.py   single-depth bottom-hole heat balance
  model1_diameter_scaling.py, model1_depth_limits.py   scaling analyses
  model2_spallation.py    thermoelastic quench fracture vs confining stress
  model3_optimiser.py     coupled drilling and thermal optimiser, siting map
  model4_hole_stability.py            creep closure and convection (Péclet) check
  model5_convergence_confinement.py   ground-reaction-curve hole stability
  site_evaluation.py      SiteProfile, stress profiles, mud window, integrated verdict
  comparative_sites.py    four European provinces, by data tier
  calibrate_ud1.py        rock-strength calibration against UD-1
  crosscheck_soultz.py    stability cross-check against the Soultz wells
  *_figures.py            figure generators (write PNGs to ../figures)
```

---

## Limitations (read before trusting any number)

**The site verdicts are extrapolations.** Every 400 °C target lies far below the deepest stress data at its site: 5.7 km at Soultz, 10.9 km at United Downs, and no stress data at all for the two speculative sites. Below the data, stresses follow the measured linear trends, capped at frictional equilibrium. That is a defensible guess, not a measurement, and it is the largest uncertainty in every verdict. Temperatures are extrapolated too, below 5 km at Soultz and United Downs, at gradients (35 and 28 °C/km) carried over from v1.0 without a source.

**What is calibrated, and on what.** UD-1 at United Downs is the calibration well: its rock strength comes from the breakouts in the BGS interpretation of its image log. Soultz's model is checked against GPK3 and GPK4 but not fitted. Its strength is the published lab UCS, which is usually higher than in-situ rock-mass strength (UD-1's weak zones came out at 118–147 MPa). Larderello and Pannonian are calibrated on nothing: their strengths (140, 160 MPa) are v1.0 values with no source, and their stresses are bounds from the faulting regime.

**What remains uncalibrated or assumed:**

- Wall cooling while UD-1 was drilled is unknown; the 0–40 K bracket moves the fitted strength by about 30 MPa. At the 400 °C targets the wall temperature comes from Model 1 and is assumed to equal the circulating fluid's, the most cooling possible.
- The calibration uses the BGS interpretation (24 breakouts, 46 m in the 12.25" section). The operator's, summarised by Reinecker et al. (2021), has 27 breakouts and 139 m, and 13 tensile fractures to the BGS's 2.
- Tensile fractures: with zero wall tensile strength (chosen because a critically stressed wall fails from existing flaws), the model predicts them at every depth at both UD-1 and Soultz, where they were logged only in parts of the hole. It reports tensile fractures but cannot locate them.
- The breakout limit is 90° (Zoback 2007). UD-1's widest breakout in a trouble-free section is 63°; verdicts at 63° are printed as a sensitivity.
- The mud window's upper bound is Shmin less 0.05 SG, a typical drilling margin, not a site value.
- Pore pressure is hydrostatic wherever it isn't reported. The deep Great Hungarian Plain is reported overpressured, which the Pannonian verdict does not include.
- Rock stiffness and thermal conductivity are generic granite values at every site.
- The breakout check assumes a vertical hole. Model 5's isotropic closure check is still in total stress, and Model 4's creep uses its own lithostatic stress estimate rather than the site profiles.

**Model limitations, unchanged from v1.0:**

- The models are low-dimensional. Model 1 is 1-D, and Models 2 and 5 are axisymmetric or plane-strain.
- The rock-coupling and stability models are mostly quasi-steady, so they do not resolve time dependence. Multi-year reservoir thermal decline, which drives the production economics, is left out on purpose, since the drilling face meets fresh rock.
- Two couplings are uncertain and would need experimental calibration: the quench-to-effective-MSE parameters (`chi_macro`, `chi_micro`), and the high-temperature rock flow law. Both come from the literature, and results that depend on them are reported as sensitivities.
