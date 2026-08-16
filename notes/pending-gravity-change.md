# Pending: surface gravity to ~12.81 m/s2

**Decided, not yet applied.** Recorded here because applying it is a larger job
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

    density 6.0 g/cm3  ->  M = 1.880 Earth masses,  g = 12.81 m/s2

Earth-like life remains reasonable at that gravity, which is the constraint that
matters for this world.

## What it invalidates, which is nearly everything geographic

Orogen scales maximum terrain relief as 1/g. The geography in `source/` therefore
cannot be separated from the gravity it was generated with -- `CLAUDE.md` already
says this, in the context of gravity being canonical rather than derived. The
change requires:

1. **A fresh base terrain generation** at the new gravity. Every build in
   `source/` is superseded, including `carved-zoned-v5`.
2. **A new basin catalogue.** Basin ids are stable across carve iterations
   because detection runs on the pre-conditioning surface, but that surface moves
   when relief scales, so `basinCatalogue` will not be `2d1f8e57` any more. No
   existing carve verdict transfers.
3. **The full geography loop again**: hydrography, boundary conditions, climate,
   carve verdict. The methods and the tooling all survive; the numbers do not.
4. **A re-derived flux.** Relief affects land albedo through lithology exposure
   and affects circulation through orography, so the baseline flux has to be
   re-measured rather than carried over.

The pedology and biosphere chains survive as method. Their calibrations should be
revisited against a mantle that is not Earth's, per the differentiation point
above.

## Sequencing

Best done at a clean point rather than mid-iteration, because it restarts the
geography loop from the beginning. The work in flight when this was decided --
the iteration-2 climate on `carved-zoned-v5` -- is worth finishing only for what
it teaches about method, not for its numbers.

`derive()` in `run_exoplasim.py` already refuses a mass and gravity that
disagree, so the two cannot be changed independently by accident.
