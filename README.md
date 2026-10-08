# Quenchwell: models of circulating temperature and borehole stability for drilling into superhot rock

Quenchwell is an open-source tool, written in Python, for screening sites for drilling into rock at around 400 °C. For a given site it takes the temperature profile, in-situ stresses, rock strength and well design, and selects the drill pipe and flow rate that keep the drilling tools below their temperature limit and carry the rock cuttings out of the hole. It then calculates the temperature and thermal stress of the borehole wall and the range of drilling-fluid weights that keeps the wall stable, and reports a verdict together with the basis for that verdict. The tool is applied here to five candidate sites.

The circulation model has been validated against temperatures measured in the insulated drill pipe trial in Utah FORGE well 16B. The stability model has been calibrated against breakouts in the UD-1 well at United Downs, Cornwall, and checked against wells at Soultz-sous-Forêts and Newberry.

At Soultz and United Downs, where stress and rock strength have been measured or calibrated, the borehole appears stable while circulating, with a commercially available insulated drill pipe and a drilling fluid of about 1.07 SG. A higher flow, if a rig can deliver it, or a better-insulated pipe that has so far only been modelled, would allow drilling with fluid at the weight of water. The stability depends on the circulating fluid cooling the wall. When circulation stops the wall reheats, and over the day or more needed to change a bit at these depths it approaches rock temperature, at which the range of workable fluid weights almost closes at Soultz and closes at United Downs. Whether the hole remains open through a trip has not yet been modelled.

The targets at both sites lie 6 to 11 km below the deepest stress measurements, so the results are projections from much shallower data.

The name is taken from Quenchwell, a hamlet in Feock, Cornwall, a few miles from United Downs.

## 1. Introduction

Rock at temperatures of around 400 °C, often called superhot rock, could yield considerably more power per well than conventional geothermal systems, because water at those conditions is supercritical and carries much more energy. Reaching it is primarily a drilling problem. Electronics and seals in the bottom-hole assembly fail above roughly 200 °C, hot crystalline rock is slow to cut, and a hot, highly stressed borehole wall may break out or creep shut. The deepest borehole yet drilled, the Kola superdeep well, reached 12.3 km, where the rock proved hotter and more plastic than expected; in most continental crust, 400 °C lies deeper still. Where it is shallower, as at Larderello in Italy, the Venelle-2 well recorded temperatures of 507 to 517 °C at 2.9 km depth (Bertani et al., 2018).

The concept examined here relies on a single circulating loop to address all three problems. Cold fluid is pumped down an insulated pipe within the drill string to the bit, where it keeps the tools below their survival limit and quenches the rock face. The quench may crack the rock directly (thermal spallation) or weaken it so that the bit cuts faster, and it cools the borehole wall, which strengthens the rock and slows creep. The fluid then returns up the annulus in contact with the hot rock and brings heat back to surface while the well is being drilled. As a production system, a single closed loop is limited by conduction through the rock to a few megawatts of heat, as Model 1 shows below and as the closed-loop field results and studies cited in Section 2 find. The fivefold to tenfold advantage per well claimed for superhot rock assumes open-loop flow through fractured rock (Clean Air Task Force, 2022), which is how Mazama Energy plans to develop Newberry. The models here are therefore concerned with the loop as a means of drilling the well, and the heat it returns is treated as a by-product.

The models were built to test each part of this concept from first principles, using IAPWS-95 water properties and rock and stress parameters taken from the literature. They are one-dimensional or axisymmetric, and mostly steady-state. The circulation model (Model 1) has been validated against field measurements to within a few degrees (Section 3.1). The models of spallation, drilling rate and creep (Models 2 to 4) rest on parameters taken from the literature and remain order-of-magnitude estimates. Their purpose is to identify the constraints that govern feasibility and to quantify the trade-offs between them. A full reservoir or geomechanical simulation would still be needed to design a well.

## 2. Related work

The use of insulated drill pipe to keep the bit cool in hot rock is well established. Sandia National Laboratories and Drill Cool Systems built insulated drill pipe in the 1990s and tested it in the laboratory and in a geothermal well, comparing circulating temperatures with those in conventional pipe (Finger et al., 2000). A later DOE-funded project to develop it for high-temperature drilling stopped after its first phase, because a pipe strong enough for the hole left too small a bore for the flow (Champness et al., 2008). Eavor ran internally and externally coated pipe at Utah FORGE in 2023 and has published its field experience and coating properties (Vetsak et al., 2024). Circulating temperatures in wells have been calculated analytically since Raymond (1969), Holmes and Swift (1970) and Kabir et al. (1996), and numerically with Sandia's GEOTEMP (Mitchell, 1981) and GEOTEMP2 (Mitchell et al., 1984), which Ujyo et al. (2026) have extended to supercritical water as GEOTEMPSC; they find that circulation cools the rock only a few metres from the wall. Model 1 is a steady-state model of the same kind. Wu et al. (2025) compared conventional, coated and dual-wall pipe in deep vertical and horizontal wells with a transient thermo-hydraulic model, and Pearce and Pink (2024, 2025) concluded that superhot wells can be drilled by combining existing technologies, with insulated pipe the main means of keeping tools below their 175 to 200 °C limits.

The mechanics of borehole stability used here are also established. The stress concentration around a borehole, breakout formation under a Mohr-Coulomb criterion and the use of breakout width to infer stress are standard (Zoback, 2007). Cooling a borehole wall adds a thermal hoop stress that suppresses breakouts and promotes tensile fractures, the combination studied by Hals and Berre (2012). In practice, the Japan Beyond-Brittle Project describes drilling at Kakkonda into rock above 500 °C with continuous circulation, and sets a target of cooling the hole below 160 °C (Muraoka et al., 2014), and Kruszewski and Wittig (2018) reviewed the failures in 20 wells drilled into supercritical conditions, many of them in casing and cement under thermal load. On the production side, Eavor's Geretsried project now reports about 8.5 MW of heat from its first loop, more than 2 MW per lateral pair falling to an expected 1.3 MW after five years (Eavor, 2026), and Beckers and Johnston (2022) estimated about 22 MW of heat and 2.2 MW of electricity from a 7.5 km loop with more than 90 km of laterals at 30 °C/km.

The contribution of this repository is to link these elements for named sites in a single open tool, whose results are fixed by its test suite. The circulating temperature from a validated model is carried through to the thermal stress at the wall, to the range of drilling-fluid weights that keeps the wall stable, and to the depth of the target below the site's data. The stability model has been calibrated in one well and checked in two others, and the circulation model has been compared with a published simulator and validated against field measurements.

## 3. Drilling with a single loop

The thermal balance of the loop (Model 1, `model1_coupled.py`) is solved for the fluid flowing down the drill string and back up the annulus. The well is represented by its casing and cement, a string of drill pipe, heavy-weight pipe and collars, and a pipe wall made up of its steel and insulation layers. The model includes the heat generated by friction in the pipe, annulus and bit nozzles, and the time for which each depth of the wall has been exposed to circulation. For rock at 450 °C at 12.4 km depth, at the 30 kg/s needed to lift the cuttings in an 8.67-inch hole, fluid reaches the bit at 390 °C with conventional pipe, 222 °C with pipe internally coated with NOV's TK-Drakōn, and 56 °C with the dual-wall pipe of Xiao et al. (2022) (Fig. 1). Tool survival depends above all on the insulation of the downgoing pipe. Water in the loop remains dense and supercritical below about 2 km, with no flashing. The heat returned is limited by conduction through the surrounding rock, to a few megawatts per loop once the rock has been circulated against for a year (Fig. 1); the higher figures in the site table are for rock freshly exposed while drilling. Because that conduction depends only logarithmically on the borehole radius, a tenfold increase in diameter yields only about 1.9 times the heat, and output scales instead with the length of hole in contact with hot rock.

![Model 1 coupled counterflow temperature profiles](figures/model1_profiles.png)

*Fig. 1. Temperature of the downgoing fluid, the returning fluid and the rock along a 12.4 km well to 450 °C rock (Model 1), at 30 kg/s, for conventional pipe, TK-Drakōn coated pipe and dual-wall pipe, with the rock exposed for a year. At this flow the TK-Drakōn case is above the 200 °C tool limit; a higher flow brings it below.*

### 3.1 Checking the circulation model

Model 1 was first compared with the transient simulator of Wu et al. (2025), with their published well, string, fluid and flow and nothing fitted (`benchmark_utaustin.py`). The comparison identified one omission in Model 1, the heat generated by friction. Their dual-wall pipe loses about 70 MPa to friction, which warms the mud by about 18 °C; with this heat included, Model 1 agrees with their results for conventional and dual-wall pipe to within 7 °C in the vertical well and 2 °C in the horizontal well. The results for coated pipe do not agree. With the 1 mm coatings described in their text, Model 1 is 30 to 60 °C hotter than their results, and with the single conductivity per pipe given in their table applied across the whole wall, it is 10 to 36 °C colder. Their results lie between the two, and the paper does not state which interpretation their model uses. The standpipe pressures from Model 1 are 10 to 25% higher than theirs, of which 1.4 MPa is the loss across a bit whose nozzles they do not specify.

| Base case, 600 gal/min, 24 °C/km | Wu et al. | Model 1 | Model 1, whole-wall reading |
|---|---|---|---|
| Vertical, conventional | 179 °C | 172 °C | 172 °C |
| Vertical, internally coated | 76 °C | 126 °C | 66 °C |
| Vertical, externally coated | 121 °C | 153 °C | 91 °C |
| Vertical, dual-wall | 59 °C | 61 °C | 61 °C |
| Horizontal, conventional | 206 °C | 207 °C | 207 °C |
| Horizontal, dual-wall | 65 °C | 67 °C | 66 °C |

The field comparison uses the trial of Eavor's insulated drill pipe in Utah FORGE well 16B in May 2023 (`validate_forge16b.py`). Four bits were run in 9.5-inch hole at 60 to 70° inclination: two on regular pipe, one on a full string of insulated pipe and one with about 70% insulated pipe above regular pipe. The well's survey, casing, string, rock properties and drilling record come from the public FORGE dataset, the flow and inlet temperature minute by minute from the rig's records, and the measured MWD temperatures from Eavor's trial report. The model is run at half-hourly intervals while each bit was drilling and compared with the fluid temperature in the pipe at the MWD sensor. Without fitting, the model overestimated the temperature on the runs with regular pipe by 4 to 12 °F. A multiplier of 0.75 on the annulus heat-transfer coefficient, fitted to those runs, brings them within 5 °F. One effective conductivity for Eavor's pipe, fitted on the full insulated run, gives 4.9 W/m·K across the pipe wall, which is within about 15% of the wall resistance of a 1 mm coating at 0.47 W/m·K on steel, the coating Eavor has described (Vetsak et al., 2024). The run with partial insulation was then predicted with no further fitting.

| Run | Measured MWD | Fitted model | Model with the film unfitted |
|---|---|---|---|
| BHA 10, regular pipe | 180 °F | 184 °F | 192 °F |
| BHA 11, insulated pipe (fitted) | 149 °F | 149 °F | 149 °F |
| BHA 12, partial, blind | 164 °F | 157 °F | 160 °F |
| BHA 12 after levelling, blind | 150 to 160 °F | 159 °F | 161 °F |
| BHA 13, regular pipe | 220 °F | 218 °F | 228 °F |

Eavor estimated from the data how each change between the first two runs moved the MWD temperature. Rebuilt one factor at a time, the model attributes −9 °F to the higher flow (Eavor −12 to −14), +17 °F to the warmer inlet (+23), +9 °F to the added depth (+7) and −61 °F to the insulated pipe. Eavor's breakdown gives −47 to −49 °F for the pipe, which they describe as conservative, and their overall estimate of its benefit, against the runs before and after, is 47 to 75 °F. The site evaluations use the unfitted coefficient, since the multiplier absorbs all the effects the model omits at FORGE and there is no evidence that it applies to a 12 km well; it is reported as a sensitivity. With the coefficient unfitted, BHA 12 is predicted to within 4 to 6 °F, and the temperatures on the runs with regular pipe are overestimated, so the model errs towards a warmer wall and therefore towards caution.

![Model 1 against the FORGE 16B trial](figures/validate_forge16b.png)

*Fig. 2. Left: modelled MWD temperature while drilling each run (points), against the measured run averages (bars) and BHA 12's levelled range (shaded). Right: the change between the first two runs, one factor at a time, from the model and from Eavor's estimate.*

### 3.2 Drilling rate, creep and the brittle-ductile transition

The quench was expected to make superhot drilling much faster, but this expectation is only partly borne out. Thermal spallation (Model 2, `model2_spallation.py`) requires the thermal tension at the quenched face to exceed both the rock's tensile strength and the confining stress, which grows with depth. It is therefore possible only in crust where the horizontal stress is low relative to the vertical (a ratio K0 below about 0.5, typical of extensional settings) and the rock is still brittle (Fig. 3). Where spallation occurs, the coupled drilling model (Model 3, `model3_optimiser.py`) predicts rates of penetration up to about 7.5 times those without the quench. Elsewhere the quench lowers the energy needed to cut the rock without breaking it outright, and the gain is 1.5 to 2.3 times (Fig. 4). Similar quench-assisted fracture has been reported by Holtzman et al. (2023) and in the literature on cryogenic and thermal-shock stimulation of enhanced geothermal systems.

![Model 2 quench spallation feasibility against depth](figures/model2_feasibility.png)

*Fig. 3. Thermal tension available from the quench (red) against the tension needed to crack the rock at each stress ratio K0 (dashed). Spallation is possible in the shaded window, at low K0 and below the brittle-ductile transition.*

![Model 3 coupled optimiser siting map](figures/model3_siting_map.png)

*Fig. 4. Gain in rate of penetration from the quench across target temperature and stress ratio K0. The 7.5× band is pure spallation; elsewhere the gain is 1.5 to 2.3×.*

Cooling also keeps the hole open. Hot rock creeps slowly into an open borehole, and the hotter the rock, the faster it creeps. Model 4 (`model4_hole_stability.py`) shows that cooling the wall to 200 °C slows this creep by a factor of 10^5 to 10^6, and that convection carries too little heat back into the cooled rock to undo the effect. Where the wall does break out, Model 5 (`model5_convergence_confinement.py`) shows that the damaged zone is only millimetres deep and is held in place by the fluid pressure, well short of collapse.

The approach meets a further limit at about 400 °C, where granite stops fracturing and begins to flow: the brittle-ductile transition, reached at depths between about 5.5 and 19 km depending on the geothermal gradient. Cooling does not remove this limit for spallation. At the bit, the quench chills only a skin a few millimetres thick before it is cut away, while the rock that would have to crack lies behind it at full temperature and yields by flowing. Above about 400 °C, the quench can therefore only weaken the face, as it does in most rock. Cooling does still help keep the hole open, because the cooled shell around the wall, about 1.5 m thick after a week, keeps the hot rock away from the stress concentration at the wall (Model 4). The 400 °C targets used here sit at the transition; whether the approach can be pushed hotter rests on that hole-stability argument, which has not been tested in practice.

## 4. Borehole stability

The site verdicts depend mainly on whether the borehole wall can be kept stable. Drilling a vertical hole concentrates the in-situ stresses around its circumference, and where the resulting hoop stress exceeds the strength of the rock, the wall spalls in two lobes on opposite sides of the hole, in the direction of the minimum horizontal stress (Shmin). These features are known as breakouts. Wells are routinely drilled with breakouts present; problems arise when the breakouts become too wide. Following the empirical criterion of Zoback (2007), a site is rated GO if breakouts stay within 90° of the circumference with drilling fluid at hydrostatic weight. It is rated CONDITIONAL if a heavier fluid brings them within 90° while its pressure stays below Shmin; above Shmin the fluid would fracture the rock and be lost into it. It is rated NO-GO if any fluid heavy enough to hold breakouts within 90° would also exceed Shmin. The equations are given in the Appendix.

### 4.1 Calibration against UD-1

UD-1 was drilled into the Carnmenellis granite with drilling fluid below 1.05 SG (specific gravity relative to water) in its 12.25-inch section, and stress magnitudes at the site have been published by Reinecker et al. (2021). The British Geological Survey's interpretation of the UD-1 borehole images (BGS, 2024) records 24 breakouts in the near-vertical part of the well between 900 and 4,000 m measured depth, each with its depth and angular width. Given the stresses, the width of a breakout fixes the strength of the rock in which it formed, so the rock strength can be fitted directly (`calibrate_ud1.py`). The intervals that did not break out set a lower bound on the strength of intact rock.

The temperature of the drilling fluid at the wall was not published, and its cooling effect was therefore bracketed between 0 and 40 K. Over that range the median strength of the rock that broke out is 118 to 147 MPa, and intact rock must have a strength of at least 172 to 203 MPa. With these values the model predicts breakouts from about 2.9 km depth, compared with one logged at 2.2 km and the rest from 3.1 km (measured depth, which differs from vertical depth by less than 10 m in this section), and widths of 53 to 57° near the base of the section, compared with 50 to 53° logged (Fig. 5). The 51 breakouts recorded in the deeper, deviated part of the well were not used in the fit; they fall between the predictions for weak and intact rock. The interpretation summarised by Reinecker et al. (2021) contains more breakouts over the same interval (27, totalling 139 m against 46 m), and the calibration is accordingly specific to the BGS reading of the images.

![Breakout width with depth in UD-1, logged and predicted, and the strength each breakout implies](figures/ud1_calibration.png)

*Fig. 5. Left: breakout widths logged in UD-1 (filled, 12.25-inch section; hollow, deviated 8.5-inch section) and widths predicted for weak-zone and intact rock. Right: rock strength implied by each breakout over the 0 to 40 K range of wall cooling.*

### 4.2 Cross-check at Soultz-sous-Forêts

The Soultz wells GPK3 and GPK4 reached 5 km with breakouts, and their stress state has been characterised to that depth by Valley and Evans (2007). With nothing fitted to Soultz (`crosscheck_soultz.py`), the model predicts breakouts at 5 km with a drillable range of fluid weights in 46 of 48 combinations of rock strength, wall cooling and stress. Using the weak-zone strengths from UD-1 without adjustment, it places the first breakouts at 2.6 to 3.8 km depth; in GPK4 they appear at about 3.0 km and become dense below 3.7 km. With the laboratory strength of Soultz granite (100 to 130 MPa) and the lower bound on the maximum horizontal stress, the onset falls at 3.64 km, consistent with the analysis of Valley and Evans (2007). That a strength calibrated in one granite reproduces the failure pattern in another is encouraging, although two wells are a small basis for generalisation.

### 4.3 Cross-check at Newberry

NWG 55-29 at Newberry, Oregon, was imaged with a borehole televiewer in 2008 during a programme of injecting cold water to cool the well, at up to 277 °C. Davatzes and Hickman (2011) found breakouts throughout the volcanic rocks above about 8,610 ft and none in the granodiorite below, and derived the stresses from them with a strength taken from the neutron porosity log. They set the thermal stress of the cooled wall to zero. With their stresses and porosity-derived strength, the well's logs from the Geothermal Data Repository (Cladouhos, 2012) and nothing fitted (`crosscheck_newberry.py`), the model matches the volcanics' modal breakout width of 36° when the wall is 25 to 35 K below its static temperature. This is the same order as the 20 K of cooling implied by the log's maximum temperature against the static temperature at its base, and it is the first test of the thermal term at about 290 °C. The model also predicts breakouts over most of the logged granodiorite, where none were seen, as does the authors' own criterion. The porosity relation gives the granodiorite 60 to 82 MPa against 51 to 74 MPa in the volcanics, too little contrast to separate them, so the pattern by rock type is not reproduced. Their friction coefficient (0.55, with 0.70 as an alternative) is used for this well and for the site, because their stresses were derived with it; it differs from the 0.8 used at United Downs. The Shmin they quote at 8,420 ft corresponds to their formula at 0.70, and both cases are run.

## 5. Site assessment

The five sites were assessed for a target rock temperature of 400 °C (`comparative_sites.py`). They are not equally well characterised, and the results are reported in two groups accordingly. At Soultz and United Downs the stress magnitudes and rock strength have been measured or calibrated at the site, and the stability model has been checked against wells drilled there. At Larderello and in the Pannonian Basin, no stress magnitudes or rock strengths were found in the literature; the stresses are bounded only by the type of faulting observed (Liotta and Brogi, manuscript; Bada et al., 2007), and the strengths have no source. Newberry has stresses from a breakout analysis and strengths derived from logs, but its well check did not reproduce the breakout pattern, so it is grouped with them. Mazama Energy is developing Newberry as an open-loop enhanced geothermal system, and the table assesses drilling there; production is outside its scope. For the Pannonian Basin the table covers the transtensional part of the strike-slip regime described for the basin interior, with the maximum horizontal stress close to the vertical.

In the model, each site is drilled to the well design of Wu et al. (2025), with casing to 4 km and 8.67-inch open hole below, and water as the circulating fluid. The drill pipe is the best-insulated pipe that is commercially available, which is pipe internally coated with NOV's TK-Drakōn. The flow is the larger of the lowest flow that keeps the bit below 200 °C and the lowest that lifts the rock cuttings to surface (hole cleaning), taken as the annular velocity at which FORGE 16B was drilled, provided that the standpipe pressure remains below 51.7 MPa. At the three deep sites the flow is set by the tool limit, at 35 to 50 kg/s, which places the circulating temperature at the bit at 184 to 188 °C; at Larderello and Newberry it is set by hole cleaning, at 30 kg/s. The wall temperature for the verdict is the circulating temperature at the bit, with the rock freshly exposed at the rate of penetration from Model 3. The lower bound of the fluid-weight window, which keeps breakouts within the limit, is checked at static fluid weight, which is what the wall sees when circulation stops for a connection, for minutes while it is still cool. The upper bound, which prevents the rock from fracturing, is checked at the circulating pressure, which adds the annular friction loss of 0.02 to 0.04 SG.

| Site | Depth to 400 °C | Pipe and flow | Wall | Stability | Flow for GO | Heat returned while drilling | Basis |
|------|-----------------|---------------|------|-----------|-------------|-----------------|-------|
| Soultz-sous-Forêts (FR) | 10.7 km | TK-Drakōn, 40 kg/s | 184 °C | CONDITIONAL: GO in 6 of 9 cases; worst case needs 1.07 SG | 60 kg/s (956 gal/min) | 5.5 MW | measured stress to 5 km; laboratory strength; checked against GPK3/GPK4 |
| United Downs (UK) | 12.9 km | TK-Drakōn, 50 kg/s | 188 °C | CONDITIONAL: GO in 2 of 4 cases; worst case needs 1.07 SG | 80 kg/s (1,275 gal/min) | 6.7 MW | measured stress (Shmin to 2 km); strength calibrated on UD-1 |
| Larderello (IT) | 2.4 km | TK-Drakōn, 30 kg/s | 71 °C | GO | | 2.5 MW | stress regime only; strength unsourced; temperature measured in Venelle-2 |
| Pannonian Basin (HU) | 9.7 km | TK-Drakōn, 35 kg/s | 188 °C | CONDITIONAL: GO in 3 of 5 cases; worst case needs 1.06 SG | 70 kg/s (1,116 gal/min) | 4.5 MW | stress regime only; strength unsourced; no deep wells |
| Newberry (US) | 3.6 km | TK-Drakōn, 30 kg/s | 90 °C | GO | | 2.8 MW | stress from breakouts; strength from logs; well check not reproduced |

Where the data give a range of values (the measured range of maximum horizontal stress at Soultz, the four calibrated strengths at United Downs, the regime bounds and friction coefficients elsewhere), every combination was evaluated, and the table reports the number of cases rated GO and the drilling-fluid weight required in the worst case. The flow for GO is the lowest flow at which the same pipe gives GO in every case within the pump limit. All three exceed the 600 to 700 gal/min at which FORGE 16B was drilled, and whether a rig can deliver them at these depths has not been established. With the dual-wall pipe of Xiao et al. (2022) the wall is cooled to 48 to 64 °C and every site is GO at 30 kg/s, and the same holds with vacuum-insulated tubing; neither has been made as drill pipe. With FORGE's fitted annulus coefficient, or with the upper bound checked at static fluid weight, the verdicts are unchanged. The heat figures are early-life values for a single loop at the flow used.

At both well-characterised sites the fluid weight required is close to that of water, but this depends on cooling the wall. With the wall at rock temperature, no case at either site is rated GO: Soultz then needs at least 1.33 SG against a fracture limit of 1.34 while circulating, and at United Downs the window closes in the worst case, at 1.29 SG. UD-1's widest breakout in the section drilled without difficulty was 63°, and with that stricter limit in place of 90°, Soultz remains CONDITIONAL and United Downs is NO-GO in its worst case. With both an uncooled wall and the 63° limit, both sites are NO-GO: Soultz would need 1.58 SG and United Downs 1.54, above their fracture limits. The 63° limit is probably conservative, since UD-1 was drilled to its target through breakouts of 35 to 90° in its deviated lower section; the widths there are not directly comparable with a vertical hole, but they are consistent with the 90° criterion. The open window at both sites is therefore a result of the cooling, which is the concept's central claim and the part of it most in need of testing. The cooled verdicts hold while the fluid is circulating. A trip to change the bit takes the bottom of the hole out of circulation for a day or more at these depths, and by line-source conduction the wall loses most of its cooling within a day, so the uncooled sensitivity is the better guide to the hole's condition during a trip.

![Breakout width against depth to 400 C for each site, and how deep each site's data reach](figures/comparative_sites.png)

*Fig. 6. Left: range of breakout widths at hydrostatic fluid weight across each site's cases, against depth to 400 °C; filled markers show the well-characterised sites. Right: depth reached by each site's stress and temperature data, against its target.*

The most important qualification is the depth of extrapolation (Fig. 6). The stress data at Soultz extend to 5 km and those at United Downs to 2 km, so the 400 °C targets lie 5.7 and 10.9 km below the deepest measurement. At Newberry the target lies 1.0 km below the stress data and 0.6 km below the temperature data. Below the data, the stresses are assumed to continue along their measured trends, limited by the frictional strength of the crust. This is a reasonable assumption but an untested one, and it is the largest single uncertainty in the verdicts. The full evaluation of the Soultz site is shown in Fig. 7.

![Site evaluation dashboard for Soultz-sous-Forets](figures/site_dashboard.png)

*Fig. 7. Evaluation of the Soultz site (`site_evaluation.py`): temperature along the loop for the pipe and flow behind the verdict, measured stress profiles (shaded where data exist), and a summary of tool survival, heat, drilling rate, creep and stability.*

## 6. Limitations and further work

Beyond the extrapolation of stress with depth, several assumptions limit the results. The flows at which the deep sites become GO exceed the rate at which FORGE 16B was drilled, and no source has been found for the flows a rig can deliver at these depths. The circulation model has been validated against a single well, FORGE 16B, in rock at about 220 °C, below the superhot range, and the conductivity fitted for Eavor's pipe combines the effects of the coating and the tool joints, as neither has been published separately. Water is assumed as the circulating fluid at every site; a drilling mud would raise the standpipe pressure and could rule out the dual-wall pipe, as it does in Wu et al. (2025). The hole-cleaning requirement is taken from the annular velocity in a single well, as no published minimum could be found. Temperatures below 5 km at Soultz and United Downs are extrapolated at gradients of 35 and 28 °C/km that have no published source.

The wall temperature at the targets is taken as that of the circulating fluid at the bit, which gives the greatest possible cooling and so the most favourable stability. That cooling is about 215 K with the commercial pipe at the flows used, and about 340 K with the dual-wall pipe, against 0 to 40 K in the UD-1 calibration and 25 to 35 K at Newberry. The full thermal stress is also applied at the wall surface, whereas breakouts form just behind the bit, where the wall has been cooled for minutes and only a thin skin is cold; the relief deeper in the rock, where a breakout grows, will be smaller. The reheating of the wall when circulation stops, at connections and above all during trips, is not modelled, and the verdicts with cooling apply while the fluid is circulating. Modelling of the reheating over the drilling cycle, and of failure behind the wall face, is planned; until then the stability results with cooling are an upper bound on its benefit. The rock strength at UD-1 was calibrated with the stress profile taken as given, so the two cannot be checked independently with the data available: any error in the stresses at 2 to 4 km is carried into the fitted strength, and nothing ensures the error is the same at 12.9 km. The rock strength at Soultz comes from laboratory tests, which generally overestimate the strength of rock in place, and at Newberry from a porosity relation that does not distinguish the rock in which breakouts formed from the rock in which they did not. Pore pressure is taken as hydrostatic where it has not been reported, although the deep Pannonian Basin is reported to be overpressured. The stability check assumes a vertical hole, and the drilling-fluid margin below Shmin (0.05 SG) is a typical engineering allowance, with no site-specific value behind it. The models do not cover drilling with total losses of circulation, in which no fluid returns up the annulus, or the behaviour of rock above the brittle-ductile transition beyond the creep model.

The model predicts drilling-induced tensile fractures at all depths in the calibration wells, whereas they were logged only over parts of each hole. This follows from setting the tensile strength of the wall to zero, which is appropriate for a naturally fractured, critically stressed rock mass but leaves the model unable to say where tensile fractures will occur. Two parameters in the drilling model, the reduction in cutting energy due to the quench and the high-temperature flow law of the rock, are taken from the literature and would need laboratory calibration. The long-term cooling of the rock around a producing well, which determines production economics, is deliberately left out, since the bit always meets fresh hot rock while drilling.

The most useful additional data would be downhole temperatures from a hot, deep well drilled with insulated pipe, measurements of how far the circulating fluid cools the wall, stress measurements below 5 km at either well-characterised site, the drilling-fluid temperatures recorded while UD-1 was drilled, and stress magnitudes and rock strengths at Larderello and in the Pannonian Basin. The sources for all site data are documented in `data/sites/stress_sources.md`, `data/ud1/`, `data/newberry/`, `data/forge16b/` and `data/materials.md`.

## 7. Reproducing the results

From the repository root:

```bash
pip install -r requirements.txt
python src/water_table.py         # one-time build of the IAPWS-95 property table (~2.5 min)
python src/site_evaluation.py     # full evaluation of the Soultz site
python src/comparative_sites.py   # all five sites
python src/calibrate_ud1.py       # UD-1 strength calibration and Fig. 5
python src/crosscheck_soultz.py   # Soultz cross-check
python src/crosscheck_newberry.py # Newberry cross-check
python src/benchmark_utaustin.py  # benchmark against Wu et al. (2025)
python src/validate_forge16b.py   # FORGE 16B validation and Fig. 2 (about 10 min)
```

Each model file also runs on its own and prints its analysis, for example `python src/model1_coupled.py`. The scripts `src/*_figures.py` regenerate the remaining figures in `figures/`. The test suite is run with `pytest`, which takes a few seconds and doesn't need the property table, or `pytest -m slow` to include the full runs, which do (the table is built on first use if it isn't there). A set of reference results is held fixed by the tests, so that any change to the models shows up as a difference. The FORGE and Newberry extracts can be rebuilt from the public datasets with the scripts in their data folders.

```
data/forge16b/            FORGE 16B survey, temperatures, string, drilling record and trial results
data/newberry/            NWG 55-29 logs and notes on Davatzes and Hickman (2011)
data/benchmarks/          inputs and results of Wu et al. (2025)
data/ud1/                 UD-1 sources: notes on Reinecker et al. (2021), BGS image-log interpretation
data/sites/               stress, strength and temperature sources for every site
data/materials.md         pipe, coating and cement properties, and which pipe can be bought
figures/                  generated figures
tests/                    test suite; golden/ holds the fixed reference results
src/
  geo_constants.py        physical constants, each with its source
  water_props.py, water_table.py      IAPWS-95 water properties
  well_geometry.py        well, casing, string and named pipe types
  model1_coupled.py       thermal balance and hydraulics of the circulating loop
  model1_steadystate.py, model1_diameter_scaling.py, model1_depth_limits.py
  model2_spallation.py    quench fracture against the confining stress
  model3_optimiser.py     coupled drilling and thermal model
  model4_hole_stability.py            creep closure of the cooled hole
  model5_convergence_confinement.py   yielding and breakout at the wall
  site_evaluation.py      site data, flow and pipe choice, stability verdict
  comparative_sites.py    the five sites
  calibrate_ud1.py        rock-strength calibration against UD-1
  crosscheck_soultz.py    cross-check against the Soultz wells
  crosscheck_newberry.py  cross-check against NWG 55-29
  benchmark_utaustin.py   benchmark against Wu et al. (2025)
  validate_forge16b.py    validation against the FORGE 16B trial
  *_figures.py            figure scripts
```

## Data credits

The FORGE data are from the Geothermal Data Repository, used under CC-BY 4.0: Utah FORGE: Well 16B(78)-32 Drilling Data (McLennan et al., 2023, https://doi.org/10.15121/1998591) and Utah FORGE: Deep Wells Temperature Surveys as of September 2022 (Jones, 2022, https://doi.org/10.15121/1893533). The NWG 55-29 logs are from Newberry EGS Demonstration: Well 55-29 Stimulation Data (Cladouhos, AltaRock Energy, https://doi.org/10.15121/1357916), also CC-BY 4.0. The UD-1 image interpretation contains NERC materials ©NERC 2024. The extracts in `data/` are derived from these datasets, and the README in each folder lists the files and how to rebuild the extracts.

## Licence

The code is released under the MIT licence (`LICENSE`). The files in `data/` keep the licences of the datasets they come from, given above and in each folder's README.

## Appendix: equations

The equations are given in the notation used in the code.

The loop model (Model 1) solves coupled heat balances for the downgoing and returning fluid along measured depth z, with the rock at the temperature of its true vertical depth, and with the heat dissipated by friction (pressure gradient dp/dz) added to each stream:

$$\dot{m}c_p\frac{dT_d}{dz}=UA_i(T_u-T_d)+\frac{\dot{m}}{\rho}\left|\frac{dp}{dz}\right|_d$$

$$\dot{m}c_p\frac{dT_u}{dz}=UA_i(T_u-T_d)-UA_o\big(T_{\text{rock}}(z)-T_u\big)-\frac{\dot{m}}{\rho}\left|\frac{dp}{dz}\right|_u$$

$$T_d(0)=T_{\text{inj}},\qquad T_u(L)=T_d(L)+\frac{Q_{\text{face}}}{\dot{m}c_p}+\frac{\Delta p_{\text{bit}}}{\rho c_p}$$

The conductances per unit length are resistances in series: the films on each side, each layer of the pipe wall, and on the rock side any casing and cement and the conduction into the rock, which limits the heat per unit length:

$$\frac{1}{UA_i}=\frac{1}{2\pi r_b h_d}+\sum_j\frac{\ln(r_{j+1}/r_j)}{2\pi k_j}+\frac{1}{2\pi r_o h_u},\qquad \frac{1}{UA_o}=\frac{1}{2\pi r_w h_u}+\sum_c\frac{\ln(r_{c+1}/r_c)}{2\pi k_c}+\frac{\ln(r_\infty/r_r)}{2\pi k_{\text{rock}}}$$

$$q'(z)=\frac{2\pi k_{\text{rock}}\Delta T}{\ln(r_\infty/r_r)},\qquad r_\infty=r_r+2\sqrt{\alpha t(z)}$$

where t(z) is the time the wall at z has been exposed to circulation, at least one hour. Where a pipe has bare tool joints over a fraction f of its length, UA_i is the length-weighted mean of the insulated body and bare steel. Friction follows the laminar and Blasius factors, and the bit nozzles lose

$$\Delta p_{\text{bit}}=\frac{\rho Q^2}{2C_d^2A_n^2},\qquad C_d=0.95$$

The standpipe pressure is the sum of the friction losses down the pipe and up the annulus, the bit loss, and the buoyancy of the hotter annulus, the integral of (ρ_u − ρ_d) g over true vertical depth.

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

Each stress state is also required to be within the frictional strength of the crust, with a friction coefficient μ of 0.8 (Reinecker et al., 2021), or at Newberry the 0.55 and 0.70 of Davatzes and Hickman (2011):

$$\frac{S_1-P_p}{S_3-P_p}\le\left(\sqrt{1+\mu^2}+\mu\right)^2$$

## References

Bada, G., Dövényi, P., Horváth, F., Szafián, P., Windhoffer, G., 2007. Present-day stress field in the Pannonian Basin and the surrounding Alpine-Carpathian-Dinaric orogens. Földtani Közlöny 137 (3).

Beckers, K.F., Johnston, H.E., 2022. Techno-economic performance of Eavor-Loop 2.0. Proceedings, 47th Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-223.

Bertani, R., et al., 2018. The first results of the DESCRAMBLE project. Proceedings, 43rd Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-213.

British Geological Survey, 2024. United Downs 1 borehole image interpretation. NERC EDS National Geoscience Data Centre. https://doi.org/10.5285/a040a246-8691-428b-b146-13874e743b0a. Contains NERC materials ©NERC 2024.

Champness, T., Worthen, T., Finger, J., 2008. Development and application of insulated drill pipe for high temperature, high pressure drilling. Drill Cool Systems Inc., final report to DOE/NETL.

Cladouhos, T.T., 2012. Newberry EGS Demonstration: Well 55-29 Stimulation Data. AltaRock Energy. Geothermal Data Repository, submission 271. https://doi.org/10.15121/1357916. CC-BY 4.0.

Clean Air Task Force, 2022. Superhot Rock Energy: A Vision for Firm, Global Zero-Carbon Energy.

Davatzes, N.C., Hickman, S.H., 2011. Preliminary analysis of stress in the Newberry EGS well NWG 55-29. Geothermal Resources Council Transactions 35.

Eavor, 2026. Technical update from Geretsried: what we built, what we learned and what comes next. 21 May 2026. https://eavor.com/blog/technical-update-from-geretsried-what-we-built-what-we-learned-and-what-comes-next/

Finger, J.T., Jacobson, R.D., Champness, A.T., 2000. Development and testing of insulated drillpipe. IADC/SPE Drilling Conference, New Orleans. SAND2000-0285C.

Hals, K.M.D., Berre, I., 2012. Thermal fracturing of geothermal wells and the effects of borehole orientation. arXiv:1212.2763.

Holmes, C.S., Swift, S.C., 1970. Calculation of circulating mud temperatures. Journal of Petroleum Technology 22 (6), 670–674.

Holtzman, B.K., Groebner, N., Mittal, T., 2023. Quench-spallation drilling: a novel drilling head design for routine heat mining above the brittle-ductile transition. Geothermal Resources Council Transactions 47, 2960–2978.

Jones, C., 2022. Utah FORGE: Deep Wells Temperature Surveys as of September 2022. Geothermal Data Repository, submission 1421. https://doi.org/10.15121/1893533. CC-BY 4.0.

Kabir, C.S., Hasan, A.R., Kouba, G.E., Ameen, M.M., 1996. Determining circulating fluid temperature in drilling, workover, and well-control operations. SPE Drilling & Completion 11 (2), 74–79.

Kruszewski, M., Wittig, V., 2018. Review of failure modes in supercritical geothermal drilling projects. Geothermal Energy 6, 28.

Liotta, D., Brogi, A. Pliocene-Quaternary fault kinematics in the Larderello geothermal area (Italy): insights for the interpretation of the present stress field. Manuscript submitted to Geothermics.

McLennan, J., Mock, B., Swearingen, L., Baldwin, R., Hodder, M., Vetsak, A., Kuhns, A.T., Breland, J., England, K., 2023. Utah FORGE: Well 16B(78)-32 Drilling Data. Geothermal Data Repository, submission 1516. https://doi.org/10.15121/1998591. CC-BY 4.0.

Mitchell, R.F., 1981. Downhole temperature prediction for drilling geothermal wells. Sandia National Laboratories, SAND-81-0036C.

Mitchell, R.F., Mondy, L.A., Duda, L.E., 1984. GEOTEMP2: advanced wellbore thermal simulator. Sandia National Laboratories.

Muraoka, H., Asanuma, H., Tsuchiya, N., Ito, T., Mogi, T., Ito, H., 2014. The Japan Beyond-Brittle Project. Scientific Drilling 17, 51–59.

NOV, 2024. TK-Drakōn insulating coating. Tuboscope product flyer D392006697-MKT-001.

Pearce, R., Pink, T., 2024. Drilling for superhot geothermal energy: a technology gap analysis. Cascade Institute, technical paper 2024-6.

Pearce, R., Pink, T., 2025. Drilling for superhot geothermal energy: a technology gap analysis. Proceedings, 50th Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-229.

Raymond, L.R., 1969. Temperature distribution in a circulating drilling fluid. Journal of Petroleum Technology 21 (3), 333–341.

Reinecker, J., Gutmanis, J., Foxford, A., Cotton, L., Dalby, C., Law, R., 2021. Geothermal exploration and reservoir modelling of the United Downs deep geothermal project, Cornwall (UK). Geothermics 97, 102226.

Ujyo, S., Hyodo, M., Okabe, T., Naganawa, S., Burnell, J., 2026. A wellbore hydrothermal simulation technology applicable to supercritical conditions. Proceedings, 51st Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-230.

Valley, B., Evans, K.F., 2007. Stress state at Soultz-sous-Forêts to 5 km depth from wellbore failure and hydraulic observations. Proceedings, 32nd Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-183.

Vetsak, A., Besoiu, C., Zatonski, V., Torre, A., Hodder, M., Toews, M., 2024. The insulated drill pipe: field experience and thermal model validation. SPE/IADC Drilling Conference. https://doi.org/10.2118/217753-MS

Wu, Y., Zhang, Y., Bettir, N., Ashok, P., van Oort, E., 2025. A comprehensive evaluation of drill pipe insulation for downhole temperature management using physics-based models. Proceedings, 50th Workshop on Geothermal Reservoir Engineering, Stanford University, SGP-TR-202.

Xiao, D., Hu, Y., Wang, Y., Deng, H., Zhang, J., Tang, B., Xi, J., Tang, S., Li, G., 2022. Wellbore cooling and heat energy utilization method for deep shale gas horizontal well drilling. Applied Thermal Engineering 213, 118684.

Zoback, M.D., 2007. Reservoir Geomechanics. Cambridge University Press.
