# Evaluating a new site

This guide takes a site from its inputs to a result in JSON. It assumes the
environment in the README's section on reproducing the results: Python
3.14.4 with the packages in `requirements-lock.txt`.

## 1. Copy the template

`examples/site_template.py` defines one site, as a `SiteProfile`, with every
field documented: its units, the range `validate_site` accepts, and where its
value comes from. The values in the template are those of Soultz-sous-Forets,
the reference case. Copy the file and replace each value with the site's own.

The fields fall into four groups.

- **Temperature.** A geotherm, any function from depth in metres to rock
  temperature in °C, and the target depth and temperature. For measured
  points, `site_evaluation.layered_geotherm([(z, T), ...], gradient_below)`
  interpolates between them and extends below the deepest at a stated
  gradient. The geotherm must rise with depth through the open hole and give
  the target temperature at the target depth within 5 °C.
- **Stress.** One or more stress cases, each a linear profile of Sv, Shmin and
  SHmax with depth, with the pore pressure, the depths the magnitude data
  cover and the source. A measured range is given as several cases (Soultz
  has three, across its range of SHmax). Below the data, SHmax is capped at
  the frictional limit at the site's friction coefficient. A case may state
  its faulting regime, and its stresses are then checked against that
  regime's order.
- **Rock.** Thermal conductivity, Young's modulus (which sets the thermal
  stress) and strength. Strength is given as one or more strength cases,
  each a uniaxial compressive strength; each is run against every stress
  case.
- **Sources.** `data_basis` records, for the stress, strength, pore pressure,
  temperature and the check against a real well, the basis of the value
  (measured, calibrated, extrapolated and so on, from `site_evaluation.BASES`)
  and its source. The site is evidence-based only if its stresses and
  strength are measured or calibrated there and the stability model has been
  checked against a well at the site; otherwise it is reported as
  speculative, in a separate table. Sources for the conductivity and the
  modulus go under `conductivity` and `stiffness`.

## 2. Validate the inputs

```
python -c "import sys; sys.path[:0] = ['src', 'examples']; \
import site_evaluation as se, site_template as t; print(se.validate_site(t.make_site()))"
```

`validate_site` raises an error listing every problem it finds: a value
outside its plausible range (most often a unit slip, kilometres for metres or
megapascals for pascals), temperature falling down the open hole, stresses
out of order, or pore pressure above Shmin. It returns, as warnings, any
input without a source and any temperature fall above the casing shoe.

## 3. Evaluate the site

`examples/new_site.py` builds the site from the template, validates it and,
with `--evaluate`, runs it:

```
python examples/new_site.py --evaluate --json result.json
```

An evaluation takes a few minutes for a site drilled with water and up to
about fifteen for one that needs a weighted fluid, because the fluid's weight
is iterated with the circulation it sets. With the template unchanged, the
script reports that the site is the repository's Soultz and compares its
row of the site table with the pinned one.

To run one of the repository's sites from the command line:

```
python src/site_evaluation.py --site soultz --json soultz.json
```

## 4. Read the result

The printed report and the JSON hold the same result. The parts, in order:

1. **Inputs and their sources**, with the site's tier.
2. **The operating programme**: the drill pipe, the flow, the drilling fluid
   and how many passes its weight took to converge, the temperature at the
   bit, the standpipe pressure and the circulating friction as an equivalent
   density.
3. **The cycle**: the verdict in each state (drilling, connection, the
   connection on freshly drilled rock, and the trip) for the case that
   decides it, with the drilling and trip fluids, their limits and the safe
   pause.
4. **The common windows**: the fluid's density at surface that suits every
   stress and strength case in every state down to the casing shoe, with the
   bound, state, case and depth that set each limit. Where no one fluid
   suits every case, the deciding pair of cases is named.
5. **The tool and the coating**: the tool's temperature in each state, the
   circulation schedule it needs where FORGE's practice does not keep it
   below 200 °C, and the coating against its rating.
6. **The margins**: the range of flow over which the verdict holds, and the
   fall in Shmin that would close the window.
7. **The status**: GO, CONDITIONAL with its conditions, NO-GO with the parts
   that fail, or INDETERMINATE where a solver did not converge or the coating
   is above a rating whose qualification is unknown.

The JSON also records the git commit and the package versions that produced
it.
