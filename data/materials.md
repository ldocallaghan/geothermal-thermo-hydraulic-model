# Material properties for the circulation model

| Property | Value | Source |
|---|---|---|
| Steel (drill pipe, casing) | 45 W/m·K | Wu et al. (2025), Tables 3 to 5 |
| Internally coated drill pipe | 1 mm coating, 0.47 W/m·K | Vetsak et al. (2024), as cited by Wu et al. (2025) |
| Externally coated drill pipe | 1 mm coating, 1.31 W/m·K | Vetsak et al. (2024), as cited by Wu et al. (2025) |
| Dual-wall drill pipe | phenolic resin fill, 80 mm bore, 0.092 W/m·K | Xiao et al. (2022), as cited by Wu et al. (2025) |
| Vacuum-insulated tubing | about 0.04 W/m·K averaged over body and connection; body as low as 0.001–0.005; industry range 0.02–0.08 | Dufrene et al. (2025) |
| Cement sheath | 1 W/m·K as a central value, sensitivity 0.1–10 | Dufrene et al. (2025), simulation choice |

Dufrene, C., Cambre, E., Sfeir, J., Kozikowski, M., Bois, A.-P. (2025).
"Enhancing Thermal Efficiency in Super-Hot Geothermal Systems: Optimizing
Premium Connections for Advanced Vacuum Insulated Tubulars." *Proceedings, 50th
Workshop on Geothermal Reservoir Engineering*, Stanford, SGP-TR-229. They find
70–80% of a vacuum-insulated string's heat loss is at the connections.

Vacuum-insulated tubing is a production tubular. None of the sources here
describes a drill pipe built that way, so for drilling it is a sensitivity.
The model's previous 0.02 W/m·K matches the best-case end of the vacuum-tubing
range, with connections; 0.10 W/m·K matches none of the named products.

## Hole cleaning

No citable minimum annular velocity could be checked. PetroWiki's "Hole
cleaning" page now redirects to OnePetro, the open review and thesis found
were not reachable, and Finger and Blankenship (2010), *Handbook of Best
Practices for Geothermal Drilling* (SAND2010-6048), calls for "high annular
velocity to lift the cuttings" without a number.

The model therefore takes its hole-cleaning reference from practice: the
annular velocity at which a comparable hole was actually drilled and cleaned.
FORGE 16B was drilled through the trial interval, at 60 to 70° inclination in
9-1/2-inch open hole, at 600 to 700 gal/min. The rig's daily reports give the
annular velocity around the 5-1/2-inch drill pipe as 245.1 ft/min at 600 gal/min
(1.245 m/s) and 286.0 ft/min at 700 gal/min (`data/forge16b/rig_hydraulics.csv`).
The reference is the lower, 1.245 m/s. It is a deviated-hole value, so it is
also used for vertical hole, where cuttings settle less readily and it is
conservative. No vertical value from practice is available: Wu et al. (2025)
run 600 gal/min in a 0.22 m hole with 139.7 mm pipe (1.67 m/s), but as a
modelling choice, not a field observation.

## Bit nozzles

Nozzle pressure drop is ρQ²/(2C_d²A²) with a discharge coefficient
C_d = 0.95. This reproduces the FORGE rig's own figure: 192 psi for eight 14/32-inch nozzles at 600 gal/min and
8.40 lb/gal, against 193 psi in the daily report. Where a bit's nozzles aren't
known, the model uses that FORGE bit's total flow area, 1.203 in² (776 mm²),
as a stated default.
