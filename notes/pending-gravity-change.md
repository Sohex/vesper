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

    g = 12.81 m/s2  ->  M = 1.881009 Earth masses,  density 6.002 g/cm3

Gravity is the declared value and mass follows, because Orogen scales relief by
gravity and the terrain therefore cannot be separated from it. `derive()` refuses
a mass and gravity that disagree, so the pair cannot drift apart.

Earth-like life remains reasonable at that gravity, which is the constraint that
matters for this world.

## Correction: it invalidates far less than this note first claimed

Written before the generator was tested, and wrong in the direction that mattered.
Orogen runs its pipeline in model units and applies the 1/g relief scaling only
at the model-unit-to-km conversion on export. Two builds differing only in
gravity are bit-identical in every hash. Verified upstream on two 30k-region
planets: all nine hashes equal, peak elevation 5.769 km against 4.593 km.

So of the four things this note predicted:

- **The terrain hash does not move.** Which makes it insufficient identity: a
  build from another gravity would pass the allowlist while every vertical
  quantity was off by the ratio. `lib/orogen.py` now reads `gravityMS2` and the
  preflight checks it against config.
- **The basin catalogue does not move.** All 3,629 ids stay valid.
- **Existing carve verdicts transfer and replay.** They *should* still be
  recomputed, because relief compresses by 20% and the water balance follows the
  climate -- but that is a choice about accuracy, not a forced restart.
- **Erosion did not respond to gravity either**, having also run in model units.
  The landscape's shape is identical; only the vertical scale changes.

The relief prediction was right and the mechanism behind it was wrong, which is
the least useful way to be right.

## What it does invalidate

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


## Status

Base regeneration requested from World Orogen on 2026-08-16 at g = 12.81, no
carve list, with the crust/fill zoning and the cover-chain fix. `planet.yaml` is
deliberately NOT yet edited: the config and the artifacts must move together, and
the terrain is the thing that takes time.

When the build lands, in order:

1. Register the hash in `lib/orogen.py`; expect `basinCatalogue` to have moved,
   which invalidates every carve verdict rather than merely ageing it.
2. Edit `planet.yaml`: `gravity_m_s2` 12.81, `mass_earth` 1.881009,
   `source_build` to the new build.
3. `python scripts/check_consistency.py`, which will list everything stale.
4. Rebuild hydrography, boundary conditions, albedo, roughness.
5. Re-derive the baseline flux from scratch. Do not carry 0.945 across: relief
   compresses by about 20%, which moves lithology exposure and orography, and the
   flux was measured against neither.
6. Carve loop again from the pre-carve base.

Prediction made before the build exists, so it can be scored: relief scales as
1/g, so the highest point should fall from 5.769 km to about 4.6 km.
