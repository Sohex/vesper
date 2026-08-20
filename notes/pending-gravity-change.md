# Pending: surface gravity to ~12.81 m/s2

**Applied** (see Status below); kept as the record of the decision, because
applying it was a larger job
than editing a config line, and because a decision that lives only in a
conversation is a decision that gets silently passed over.

## Why

Current declared gravity is 10.1989 m/s2 at 1.20 Earth radii, which implies:

    M = (g/g_earth) * R^2 = 1.4976 Earth masses
    mean density = 4.78 g/cm3

Earth is 5.51 g/cm3, and at 1.2 Earth radii self-compression alone should push an
Earth-composition body *above* that, toward 6. So 4.78 does not describe a
scaled-up Earth. It sits close to the iron-free MgSiO3 mass-radius relation, and
implies a planet with a negligible metallic core.

That has consequences the project was not intending to assert:

- **No dynamo**, so no intrinsic magnetic field. That bears on atmospheric
  retention against a K dwarf's wind, which this project has otherwise treated as
  unproblematic.
- **A different differentiation history**, and therefore a different
  basalt-to-granite production ratio. `pedology/` keys its whole texture chain off
  a granite/basalt quartz distinction calibrated on Earth. An iron-free mantle is
  not obviously entitled to that calibration.

## The target

    g = 12.81 m/s2  ->  M = 1.881009 Earth masses,  density 6.002 g/cm3

Gravity is the declared value and mass follows, because Orogen scales relief by
gravity and the terrain therefore cannot be separated from it. `derive()` refuses
a mass and gravity that disagree, so the pair cannot drift apart.

Earth-like life remains reasonable at that gravity, which is the constraint that
matters for this world.

## What a gravity change does not move

Orogen runs in model units and applies the 1/g relief scaling only at export,
verified upstream on two 30k-region planets: all nine hashes equal, peak
5.769 km against 4.593 km. So the terrain hash does not move -- which makes it
insufficient identity, and `lib/orogen.py` now reads `gravityMS2` and the
preflight checks it against config -- the basin catalogue keeps all 3,629 ids,
and existing carve verdicts replay; recomputing them is an accuracy choice
(relief compresses 20% and the water balance follows the climate), not a
forced restart. Erosion also ran in model units, so the landscape's shape is
identical and only the vertical scale changes.

## What it requires

Orogen applies the 1/g relief scaling only at export (the correction above),
so the gravity change itself moves no hash and no basin id. What it requires:

1. **A re-export at the new gravity** -- delivered here as a fresh base
   generation, because the request bundles the crust/fill zoning and the
   cover-chain fix, and those changes, not gravity, are what can move the
   terrain hash and the basin catalogue.
2. **Recomputed carve verdicts**, as an accuracy choice rather than a forced
   restart: relief compresses by about 20% and the water balance follows the
   climate.
3. **The full geography loop again** on the new base: hydrography, boundary
   conditions, climate, carve verdict. The methods and the tooling all
   survive.
4. **A re-derived flux.** Relief affects land albedo through lithology
   exposure and affects circulation through orography, so the baseline flux
   has to be re-measured rather than carried over.

The pedology and biosphere chains survive as method; their calibrations carry
their own sourced/declared labels.

## Status

Applied. `planet.yaml` carries 12.81 and 1.881009, the build is registered,
and the geography loop restarted from the new base. The pre-build
prediction -- relief scales as 1/g, so the peak should fall from 5.769 km to
about 4.6 km -- was confirmed at 4.593 km.
