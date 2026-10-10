# Quenchwell: a summary for readers

*Draft for the author to edit. Figures are from v1.3.3, pinned in
`tests/golden/v133_site_table.json`.*

## What the tool does

Quenchwell is a set of screening models for drilling into rock at about
400 °C with a single circulating loop. It calculates the temperature of the
circulating fluid and of the rock around the hole, the stability of the
borehole wall through the drilling cycle (drilling, connections and trips),
and the fluid weight, flow and circulation schedule that the cycle requires,
together with the temperatures the downhole tools and the coated drill pipe
reach. The circulation model is checked against the FORGE 16B
insulated-drill-pipe trial and against the calculations of Wu et al. (2025);
the stability model is calibrated against breakouts logged in UD-1 and
checked against the Soultz wells GPK3 and GPK4 and against NWG 55-29 at
Newberry. Every input carries its source, and a site is reported as
evidence-based only where its stresses and rock strength were measured there
and the stability model has been checked against a well at the site.

Each evaluation ends in a status: GO; CONDITIONAL with its conditions; NO-GO
with the parts that fail; or INDETERMINATE where a prerequisite is unknown,
such as the behaviour of the pipe coating above its rated temperature.

## The reference case: Soultz-sous-Forêts

At Soultz the rock reaches 400 °C at 10.7 km. The stress magnitudes of
Valley and Evans (2007) cover 1.5 to 5 km and are extrapolated to the
target, and the rock strength is the laboratory range of 100 to 130 MPa.

The water in the hole, at the temperatures it has while circulating, exerts
0.94 SG at the target (equivalent static density). To hold the wall as the
bit exposes the rock, the worst case of stress and strength needs 1.31 SG,
which a barite-weighted mud of 1.36 SG at surface provides. Circulating that
mud adds 0.07 SG of annular friction, and the fracture limit while
circulating, 1.30 SG, then lies below what the wall needs. Through a 44-hour
trip the column heats and the trip fluid loses 3.1 MPa at the bottom, so a
fluid that holds the wall at the end of the trip exceeds the fracture limit
at its start. In both states the window is closed for the case of the
highest SHmax and the weakest rock, and the site is NO-GO; in the other
cases it is open. The tool stays below 200 °C while drilling (197 °C) and
through connections, but running back into the hole it needs circulation for
about 70% of the run in, and the coated drill pipe exceeds its 204 °C
rating on the way in.

Across the five sites, Soultz and United Downs are NO-GO, Pannonian is
INDETERMINATE (the coating, and no single fluid suits all its stress cases),
Larderello is CONDITIONAL on a weighted fluid that the shallower open hole
needs, and Newberry is GO with water.

## Claims and evidence

| Claim | Support | Status |
|---|---|---|
| Model 1 predicts the circulating temperature in an insulated string | FORGE 16B trial; benchmark against Wu et al. (2025) | Validated at 2.5 km and about 200 °C; extrapolated to the targets |
| The rock behind the wall reheats through a trip as the wall solver gives | FORGE 16B trips, from the surface returns | Checked; the record is too coarse to confirm or reject the reheating |
| Breakout widths follow from the stresses and a calibrated strength | BGS interpretation of the UD-1 image logs | Calibrated at UD-1 |
| The model places breakout onset at Soultz without fitting | GPK3 and GPK4 (Valley and Evans, 2007) | Cross-checked |
| The thermal stress of the cooled wall is of the right size | NWG 55-29 volcanics (Davatzes and Hickman, 2011) | Matched at 25 to 35 K of cooling; the granodiorite is not reproduced |
| The stresses at the targets | Measurements to 5 km at Soultz, and shallower at the other sites | Extrapolated; the largest single uncertainty |
| The rheology of the weighted mud | One measured barite-weighted mud (Anawe and Folayan, 2018) | A single pair, with a sensitivity |
| The coating's behaviour above 204 °C | The manufacturer's rating only | Unknown, so INDETERMINATE |

## Two assumptions

**Breakouts form as the bit exposes the rock.** Moore et al. (2011) found
breakouts already present in images logged at the bit, minutes after
drilling, in four boreholes in sediment. The model applies this to hot
granite and judges the drilling state with the rock at formation
temperature. If breakouts took an hour to form, circulation would
cool the rock in the meantime and the drilling fluid at Soultz would fall
from 1.31 to 1.18 SG; the trip would still close the window.

**The breakout's width is judged at the wall, with the temperature of the
rock at the depth the uncooled breakout reaches.** That depth is about
1.6 cm at Soultz. A check with the full stress field over the first minutes
of cooling finds the failed depth growing to 1.7 cm after one minute and
falling below the uncooled depth within five, so the heuristic does not
hide a worse case over that time.

## The open question about cooling

The original premise was that cooling the wall would hold it. In the model,
the rock that decides the breakout is still at formation temperature when
the bit exposes it, so the fluid's weight has to hold the wall, and cooling
serves the tools and slows creep. Whether this holds in hot crystalline rock
depends on how quickly breakouts form there, which has not been measured.

## What data would settle it

- Image logs at the bit and repeated hours to days later in hot crystalline
  rock, to show when breakouts form. The FORGE image logs are being assessed
  for this.
- Stress measurements below 5 km at Soultz and near the target depths
  elsewhere, which would replace the extrapolation.
- Downhole temperatures at the bit through circulation, connections and
  trips, to test the tool and coating results directly.
- The coating's qualification at 250 to 300 °C, and the rheology of weighted
  water-based muds above 150 °C, where conventional systems begin to break
  down.
