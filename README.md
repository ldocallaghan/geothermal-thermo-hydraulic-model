# Deep and superhot-rock geothermal drilling: a first-principles model stack and site-evaluation tool

This repository contains a set of coupled physical models, written in Python, that examine whether a single circulating coolant loop could be used to drill into rock at around 400 °C, keep the drilling tools alive, keep the borehole open and return the heat to surface. The models are applied to four European candidate sites. The borehole-stability model, which decides most of the site verdicts, has been calibrated against the breakouts logged in the UD-1 well at United Downs, Cornwall, and checked against the GPK wells at Soultz-sous-Forêts. At the two sites where stress and rock strength have been measured or calibrated, the borehole appears to be stable with a modest increase in drilling-fluid weight, but only because the circulating fluid cools the wall. Without that cooling the range of workable fluid weights almost closes, and with the stricter breakout limit suggested by UD-1 it closes entirely. Both targets also lie 6 to 11 km below the deepest stress data, so the results are projections from much shallower measurements.

## 1. Introduction

Rock at temperatures of around 400 °C, often called superhot rock, could yield considerably more power per well than conventional geothermal systems, because water at those conditions is supercritical and carries much more energy. Reaching it is primarily a drilling problem. Electronics and seals in the bottom-hole assembly fail above roughly 200 °C, hot crystalline rock is slow to cut, and a hot, highly stressed borehole wall may break out or creep shut. The deepest borehole yet drilled, the Kola superdeep well, reached 12.3 km, where the rock proved hotter and more plastic than expected; in most continental crust, 400 °C lies deeper still. Where it is shallower, as at Larderello in Italy, the Venelle-2 well recorded temperatures of 507 to 517 °C at 2.9 km depth (Bertani et al., 2018).

The concept examined here relies on a single closed loop to address all three problems. Cold fluid is pumped down an insulated pipe within the drill string to the bit, where it keeps the tools below their survival limit and quenches the rock face. The quench may crack the rock directly (thermal spallation) or weaken it so that the bit cuts faster, and it cools the borehole wall, which strengthens the rock and slows creep. The fluid then returns up the annulus in contact with the hot rock, and the heat it gains on the way up is the energy product.

The models were built to test each part of this concept from first principles, using IAPWS-95 water properties and rock and stress parameters taken from the literature. They are order-of-magnitude models: one-dimensional or axisymmetric, and mostly steady-state. Their purpose is to identify the constraints that govern feasibility and to quantify the trade-offs between them, not to replace a full reservoir or geomechanical simulation.

## 2. Drilling and producing with a single loop

The thermal balance of the loop (Model 1, `model1_coupled.py`) shows that tool survival depends above all on insulating the downgoing pipe. For rock at 450 °C at 12.4 km depth, fluid reaches the bit at 318 °C with conventional insulation, 171 °C with vacuum-insulated tubing, and 59 °C when the flow is also raised to 10 kg/s (Fig. 1). Water in the loop remains dense and supercritical below about 2 km, with no flashing. The heat a single closed loop can export is limited by conduction through the surrounding rock to around 3 to 4 MW. Because that conduction depends only logarithmically on the borehole radius, a tenfold increase in diameter yields only about 1.9 times the heat; output scales instead with the length of hole in contact with hot rock.

![Model 1 coupled counterflow temperature profiles](figures/model1_profiles.png)

*Fig. 1. Temperature of the downgoing fluid, the returning fluid and the rock along a 12.4 km loop (Model 1), for conventional insulation, vacuum-insulated tubing, and vacuum tubing at higher flow.*

The quench was expected to make superhot drilling much faster, but this expectation is only partly borne out. Thermal spallation (Model 2, `model2_spallation.py`) requires the thermal tension at the quenched face to exceed both the rock's tensile strength and the confining stress, which grows with depth. It is therefore possible only in crust where the horizontal stress is low relative to the vertical (a ratio K0 below about 0.5, typical of extensional settings) and the rock is still brittle (Fig. 2). Where spallation occurs, the coupled drilling model (Model 3, `model3_optimiser.py`) predicts rates of penetration up to about 7.5 times those without the quench. Elsewhere the quench does not break the rock outright but lowers the energy needed to cut it, and the gain is 1.5 to 2.3 times (Fig. 3). Similar quench-assisted fracture has been reported by Holtzman et al. (2023) and in the literature on cryogenic and thermal-shock stimulation of enhanced geothermal systems.

![Model 2 quench spallation feasibility against depth](figures/model2_feasibility.png)

*Fig. 2. Thermal tension available from the quench (red) against the tension needed to crack the rock at each stress ratio K0 (dashed). Spallation is possible in the shaded window, at low K0 and below the brittle-ductile transition.*

![Model 3 coupled optimiser siting map](figures/model3_siting_map.png)

*Fig. 3. Gain in rate of penetration from the quench across target temperature and stress ratio K0. The 7.5× band is pure spallation; elsewhere the gain is 1.5 to 2.3×.*

Cooling also keeps the hole open. Hot rock creeps slowly into an open borehole, and the hotter the rock, the faster it creeps. Model 4 (`model4_hole_stability.py`) shows that cooling the wall to 200 °C slows this creep by a factor of 10^5 to 10^6, and that convection carries too little heat back into the cooled rock to undo the effect. Where the wall does break out, Model 5 (`model5_convergence_confinement.py`) shows that the damaged zone is only millimetres deep and is held in place by the fluid pressure, well short of collapse.

The approach meets a further limit at about 400 °C, where granite begins to deform by flow rather than by fracture: the brittle-ductile transition, reached at depths between about 5.5 and 19 km depending on the geothermal gradient. Cooling does not remove this limit for spallation. At the bit, the quench chills only a skin a few millimetres thick before it is cut away, while the rock that would have to crack lies behind it at full temperature and yields by flowing instead. Above about 400 °C, the quench can therefore only weaken the face, as it does in most rock. Cooling does still help keep the hole open, because the cooled shell around the wall, about 1.5 m thick after a week, keeps the hot rock away from the stress concentration at the wall (Model 4). The 400 °C targets used here sit at the transition; whether the approach can be pushed hotter rests on that hole-stability argument, which has not been tested in practice.

## 3. Borehole stability

The site verdicts depend mainly on whether the borehole wall can be kept stable. Drilling a vertical hole concentrates the in-situ stresses around its circumference, and where the resulting hoop stress exceeds the strength of the rock, the wall spalls in two lobes on opposite sides of the hole, in the direction of the minimum horizontal stress (Shmin). These features are known as breakouts. Wells are routinely drilled with breakouts present; problems arise when the breakouts become too wide. Following the empirical criterion of Zoback (2007), a site is rated GO if breakouts stay within 90° of the circumference with drilling fluid at hydrostatic weight. It is rated CONDITIONAL if a heavier fluid brings them within 90° while its pressure stays below Shmin; above Shmin the fluid would fracture the rock and be lost into it. It is rated NO-GO if any fluid heavy enough to hold breakouts within 90° would also exceed Shmin. The equations are given in the Appendix.

### 3.1 Calibration against UD-1

UD-1 was drilled into the Carnmenellis granite with drilling fluid below 1.05 SG (specific gravity relative to water) in its 12.25-inch section, and stress magnitudes at the site have been published by Reinecker et al. (2021). The British Geological Survey's interpretation of the UD-1 borehole images (BGS, 2024) records 24 breakouts in the near-vertical part of the well between 900 and 4,000 m measured depth, each with its depth and angular width. Given the stresses, the width of a breakout fixes the strength of the rock in which it formed, so the rock strength can be fitted directly (`calibrate_ud1.py`). The intervals that did not break out set a lower bound on the strength of intact rock.

The temperature of the drilling fluid at the wall was not published, and its cooling effect was therefore bracketed between 0 and 40 K. Over that range the median strength of the rock that broke out is 118 to 147 MPa, and intact rock must have a strength of at least 172 to 203 MPa. With these values the model predicts breakouts from about 2.9 km depth, compared with one logged at 2.2 km and the rest from 3.1 km (measured depth, which differs from vertical depth by less than 10 m in this section), and widths of 53 to 57° near the base of the section, compared with 50 to 53° logged (Fig. 4). The 51 breakouts recorded in the deeper, deviated part of the well were not used in the fit; they fall between the predictions for weak and intact rock. The interpretation summarised by Reinecker et al. (2021) contains more breakouts over the same interval (27, totalling 139 m against 46 m), and the calibration is accordingly specific to the BGS reading of the images.

![Breakout width with depth in UD-1, logged and predicted, and the strength each breakout implies](figures/ud1_calibration.png)

*Fig. 4. Left: breakout widths logged in UD-1 (filled, 12.25-inch section; hollow, deviated 8.5-inch section) and widths predicted for weak-zone and intact rock. Right: rock strength implied by each breakout over the 0 to 40 K range of wall cooling.*

### 3.2 Cross-check at Soultz-sous-Forêts

The Soultz wells GPK3 and GPK4 reached 5 km with breakouts, and their stress state has been characterised to that depth by Valley and Evans (2007). With nothing fitted to Soultz (`crosscheck_soultz.py`), the model predicts breakouts at 5 km with a drillable range of fluid weights in 46 of 48 combinations of rock strength, wall cooling and stress. Using the weak-zone strengths from UD-1 without adjustment, it places the first breakouts at 2.6 to 3.8 km depth; in GPK4 they appear at about 3.0 km and become dense below 3.7 km. With the laboratory strength of Soultz granite (100 to 130 MPa) and the lower bound on the maximum horizontal stress, the onset falls at 3.64 km, consistent with the analysis of Valley and Evans (2007). That a strength calibrated in one granite reproduces the failure pattern in another is encouraging, although two wells are a small basis for generalisation.

## 4. Site assessment

The four sites were assessed for a target rock temperature of 400 °C (`comparative_sites.py`). They are not equally well characterised, and the results are reported in two groups accordingly. At Soultz and United Downs the stress magnitudes and rock strength have been measured or calibrated at the site, and the stability model has been checked against wells drilled there. At Larderello and in the Pannonian Basin, no stress magnitudes or rock strengths were found in the literature; the stresses are bounded only by the type of faulting observed (Liotta and Brogi, manuscript; Bada et al., 2007), and the strengths have no source. The results for those two sites show where data would be worth collecting, and they can't be compared with the first two. For the Pannonian Basin the table covers the transtensional part of the strike-slip regime described for the basin interior, with the maximum horizontal stress close to the vertical; across the full strike-slip range, with the maximum horizontal stress up to the frictional limit, the verdict at 9.7 km spans GO to CONDITIONAL with the wall cooled, and CONDITIONAL to NO-GO without cooling.

| Site | Depth to 400 °C | Stability | Heat at 10 kg/s | Basis |
|------|-----------------|-----------|-----------------|-------|
| Soultz-sous-Forêts (FR) | 10.7 km | CONDITIONAL: GO in 8 of 9 cases; worst case needs 1.02 SG | 4.1 MW | measured stress to 5 km; laboratory strength; checked against GPK3/GPK4 |
| United Downs (UK) | 12.9 km | CONDITIONAL: GO in 3 of 4 cases; worst case needs 1.04 SG | 3.9 MW | measured stress (Shmin to 2 km); strength calibrated on UD-1 |
| Larderello (IT) | 2.4 km | GO | 1.2 MW | stress regime only; strength unsourced; temperature measured in Venelle-2 |
| Pannonian Basin (HU) | 9.7 km | CONDITIONAL: GO in 4 of 5 cases; worst case needs 1.05 SG | 3.6 MW | stress regime only; strength unsourced; no deep wells |

Where the data give a range of values (the measured range of maximum horizontal stress at Soultz, the four calibrated strengths at United Downs, the regime bounds elsewhere), every combination was evaluated, and the table reports the number of cases rated GO and the drilling-fluid weight required in the worst case. In the Pannonian worst case the two horizontal stresses are equal, and the wall then yields all the way round rather than forming breakouts; a slightly heavier fluid prevents this. The heat figures are early-life values for a single loop.

At both well-characterised sites the fluid weight required is close to that of water, but this depends on cooling the wall. With the wall at rock temperature, no case at either site is rated GO: Soultz then needs at least 1.33 SG against a fracture limit of 1.37, and United Downs at least 1.29 against 1.32. UD-1's widest breakout in the section drilled without difficulty was 63°, and with that stricter limit in place of 90°, the cooled sites remain CONDITIONAL, with windows of 1.28 to 1.37 SG at Soultz and 1.30 to 1.32 SG at United Downs. With both an uncooled wall and the 63° limit, both sites are NO-GO: Soultz would need 1.58 SG and United Downs 1.54, above their fracture limits. The 63° limit is probably conservative, since UD-1 was drilled to its target through breakouts of 35 to 90° in its deviated lower section; the widths there are not directly comparable with a vertical hole, but they are consistent with the 90° criterion. The open window at both sites is therefore a result of the cooling, which is the concept's central claim and the part of it most in need of testing.

![Breakout width against depth to 400 C for each site, and how deep each site's data reach](figures/comparative_sites.png)

*Fig. 5. Left: range of breakout widths at hydrostatic fluid weight across each site's cases, against depth to 400 °C; filled markers show the well-characterised sites. Right: depth reached by each site's stress and temperature data, against its target.*

The most important qualification is the depth of extrapolation (Fig. 5). The stress data at Soultz extend to 5 km and those at United Downs to 2 km, so the 400 °C targets lie 5.7 and 10.9 km below the deepest measurement. Below the data, the stresses are assumed to continue along their measured trends, limited by the frictional strength of the crust. This is a reasonable assumption but an untested one, and it is the largest single uncertainty in the verdicts. The full evaluation of the Soultz site is shown in Fig. 6.

![Site evaluation dashboard for Soultz-sous-Forets](figures/site_dashboard.png)

*Fig. 6. Evaluation of the Soultz site (`site_evaluation.py`): temperature along the loop at production flow, measured stress profiles (shaded where data exist), and a summary of tool survival, energy, drilling rate, creep and stability.*

## 5. Limitations and further work

Beyond the extrapolation of stress with depth, several assumptions limit the results. Temperatures below 5 km at Soultz and United Downs are extrapolated at gradients of 35 and 28 °C/km that have no published source. The wall temperature at the targets is assumed to equal that of the circulating fluid, which gives the greatest possible cooling and so the most favourable stability. That cooling is 235 K at United Downs and 252 K at Soultz, about six times the 0 to 40 K range covered by the UD-1 calibration. The full thermal stress is also applied at the wall surface, whereas breakouts form just behind the bit, where the wall has been cooled for minutes and only a thin skin is cold; the relief deeper in the rock, where a breakout grows, will be smaller. The stability results with cooling are therefore an upper bound on its benefit. The rock strength at UD-1 was calibrated with the stress profile taken as given, so the two cannot be checked independently with the data available: any error in the stresses at 2 to 4 km is carried into the fitted strength, and nothing ensures the error is the same at 12.9 km. The rock strength at Soultz comes from laboratory tests, which generally overestimate the strength of rock in place; at UD-1 the calibrated weak-zone strength is lower than the Soultz laboratory range. Pore pressure is taken as hydrostatic where it has not been reported, although the deep Pannonian Basin is reported to be overpressured. Rock stiffness and thermal conductivity are generic values for granite at every site, the stability check assumes a vertical hole, and the drilling-fluid margin below Shmin (0.05 SG) is a typical engineering allowance rather than a site value.

The model predicts drilling-induced tensile fractures at all depths in both calibration wells, whereas they were logged only over parts of each hole. This follows from setting the tensile strength of the wall to zero, which is appropriate for a naturally fractured, critically stressed rock mass but leaves the model unable to say where tensile fractures will occur. Two parameters in the drilling model, the reduction in cutting energy due to the quench and the high-temperature flow law of the rock, are taken from the literature and would need laboratory calibration. The models are mostly steady-state, and the long-term cooling of the rock around a producing well, which determines production economics, is deliberately left out, since the bit always meets fresh hot rock while drilling.

The most useful additional data would be measurements of how far drilling fluid cools the wall of a hot, deep hole, stress measurements below 5 km at either well-characterised site, the drilling-fluid temperatures recorded while UD-1 was drilled, and stress magnitudes and rock strengths at Larderello and in the Pannonian Basin. The sources for all site data are documented in `data/sites/stress_sources.md` and `data/ud1/`.

## 6. Reproducing the results

From the repository root:

```bash
pip install -r requirements.txt
python src/water_table.py        # one-time build of the IAPWS-95 property table (~2.5 min)
python src/site_evaluation.py    # full evaluation of the Soultz site
python src/comparative_sites.py  # all four sites
python src/calibrate_ud1.py      # UD-1 strength calibration and Fig. 4
python src/crosscheck_soultz.py  # Soultz cross-check
```

Each model file also runs on its own and prints its analysis, for example `python src/model1_coupled.py`. The scripts `src/*_figures.py` regenerate the figures in `figures/`. The test suite is run with `pytest`, which takes a few seconds and doesn't need the property table, or `pytest -m slow` to include the full site runs, which do (the table is built on first use if it isn't there). A set of reference results is held fixed by the tests, so that any change to the models shows up as a difference.

```
data/ud1/                 UD-1 sources: notes on Reinecker et al. (2021), BGS image-log interpretation
data/sites/               stress, strength and temperature sources for every site
figures/                  generated figures
tests/                    test suite; golden/ holds the fixed reference results
src/
  geo_constants.py        physical constants, each with its source
  water_props.py, water_table.py      IAPWS-95 water properties
  model1_coupled.py       thermal balance of the coolant loop
  model1_steadystate.py, model1_diameter_scaling.py, model1_depth_limits.py
  model2_spallation.py    quench fracture against the confining stress
  model3_optimiser.py     coupled drilling and thermal model
  model4_hole_stability.py            creep closure of the cooled hole
  model5_convergence_confinement.py   yielding and breakout at the wall
  site_evaluation.py      site data, stability verdict and integrated evaluation
  comparative_sites.py    the four European sites
  calibrate_ud1.py        rock-strength calibration against UD-1
  crosscheck_soultz.py    cross-check against the Soultz wells
  *_figures.py            figure scripts
```

## Appendix: equations

The equations are given in the notation used in the code.

The loop model (Model 1) solves coupled heat balances for the downgoing and returning fluid, with the rock wall at the local rock temperature:

$$\dot{m}c_p\frac{dT_d}{dz}=UA_i(T_u-T_d)$$

$$\dot{m}c_p\frac{dT_u}{dz}=UA_i(T_u-T_d)-UA_o\big(T_{\text{rock}}(z)-T_u\big)$$

$$T_d(0)=T_{\text{inj}},\qquad T_u(L)=T_d(L)+\frac{Q_{\text{face}}}{\dot{m}c_p}$$

Heat output per unit length is limited by conduction into the rock:

$$q'(z)=\frac{2\pi k_{\text{rock}}\Delta T}{\ln(r_\infty/r_w)},\qquad r_\infty\approx 2\sqrt{\alpha t}$$

Spallation (Model 2) requires the constrained thermal stress of the quenched face to exceed the confining stress and the tensile strength, with the face temperature given by the convective half-space solution:

$$\sigma_T=\frac{E\alpha\Delta T}{1-\nu},\qquad \sigma_T-K_0\rho g z>T_0$$

$$T(x,t)=T_{\text{rock}}+(T_{\text{cold}}-T_{\text{rock}})\big[\mathrm{erfc}(\eta)-e^{-\eta^2}\mathrm{erfcx}(\eta+\beta)\big],\qquad \eta=\frac{x}{2\sqrt{\alpha t}},\quad \beta=\frac{h\sqrt{\alpha t}}{k}$$

The drilling model (Model 3) converts the quench into an effective mechanical specific energy (MSE) and rate of penetration (ROP), and returns the heat generated at the face to Model 1:

$$\text{MSE}_{\text{eff}}=\text{MSE}(1-D),\qquad \text{ROP}=\frac{P_{\text{mech}}}{\text{MSE}_{\text{eff}}\cdot A_{\text{bit}}}$$

$$Q_{\text{face}}=\text{MSE}\cdot A\cdot\text{ROP}+\rho c\cdot A\cdot\text{ROP}\cdot\Delta T+Q_{\text{cond}}$$

Creep closure (Model 4) follows an Arrhenius power law, and convective heat resupply is checked through the Péclet number:

$$\dot\varepsilon=A\sigma^{n}\exp\left(-\frac{Q}{RT}\right),\qquad \mathrm{Pe}=\frac{vL}{\alpha}\ll 1$$

Yielding at the wall (Model 5) follows the ground-reaction curve:

$$k=\frac{1+\sin\phi}{1-\sin\phi},\qquad p_{\text{cr}}=\frac{2p_0-\sigma_{cm}}{1+k}$$

The stability verdict uses the hoop stress at the wall of a vertical hole, at angle θ from the direction of the maximum horizontal stress, including the thermal stress from cooling the wall:

$$\sigma_\theta(\theta)=S_H+S_h-2(S_H-S_h)\cos 2\theta-P_w+\Delta\sigma_T,\qquad \Delta\sigma_T=-\frac{E\alpha}{1-\nu}\,(T_{\text{rock}}-T_{\text{wall}})$$

The wall fails where the Mohr-Coulomb criterion is exceeded in effective stress, with pore pressure P_p and drilling-fluid pressure P_w:

$$\sigma_\theta-P_p>\sigma_{cm}(T_{\text{wall}})+k\,(P_w-P_p)$$

which holds over a breakout of width

$$W=\arccos\!\left(-\frac{a}{d}\right),\qquad a=S_H+S_h-P_w+\Delta\sigma_T-P_p-\sigma_{cm}-k(P_w-P_p),\quad d=2(S_H-S_h)$$

Each stress state is also required to be within the frictional strength of the crust, with a friction coefficient μ of 0.8 (Reinecker et al., 2021):

$$\frac{S_1-P_p}{S_3-P_p}\le\left(\sqrt{1+\mu^2}+\mu\right)^2$$

## References

Bada, G., Dövényi, P., Horváth, F., Szafián, P., Windhoffer, G., 2007. Present-day stress field in the Pannonian Basin and the surrounding Alpine-Carpathian-Dinaric orogens. Földtani Közlöny 137 (3).

Bertani, R., et al., 2018. The first results of the DESCRAMBLE project. Proceedings, 43rd Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-213.

British Geological Survey, 2024. United Downs 1 borehole image interpretation. NERC EDS National Geoscience Data Centre. https://doi.org/10.5285/a040a246-8691-428b-b146-13874e743b0a. Contains NERC materials ©NERC 2024.

Holtzman et al., 2023. Geothermal Resources Council Transactions 47.

Liotta, D., Brogi, A. Pliocene-Quaternary fault kinematics in the Larderello geothermal area (Italy): insights for the interpretation of the present stress field. Manuscript submitted to Geothermics.

Reinecker, J., Gutmanis, J., Foxford, A., Cotton, L., Dalby, C., Law, R., 2021. Geothermal exploration and reservoir modelling of the United Downs deep geothermal project, Cornwall (UK). Geothermics 97, 102226.

Valley, B., Evans, K.F., 2007. Stress state at Soultz-sous-Forêts to 5 km depth from wellbore failure and hydraulic observations. Proceedings, 32nd Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-183.

Zoback, M.D., 2007. Reservoir Geomechanics. Cambridge University Press.
