# The re-baseline, in order, with the circularities named

Written before starting, because three of the surface fields a climate run
consumes cannot be built without a climatology, and discovering that one field at
a time costs a run each time.

## What needs a climatology before it can be built

    229  soil water capacity   pedology's soil map, which needs climate
    174-176  albedo            only for the lake compositing; the lithology part
                               needs no climate
    lakes                      surface_water.py solves them against a climatology

Everything else -- topography, land mask, roughness, forest fraction, the
lithology half of albedo -- is a pure function of the terrain and can be built
immediately.

So the first climate run on a new terrain is necessarily a **bootstrap**: it
exists to produce the climatology that the remaining fields need. Its own numbers
are not the baseline.

## Order

**0. Terrain.** `precarve-zoned-g1281`, registered. Note this is a *pre-carve*
base, so the first verdict of the loop is iteration 1 again, though the existing
verdict transfers if we choose to replay it rather than recompute.

**1. Config.** Gravity to 12.81, mass to 1.881009, `source_build` to the new
build. Then `scripts/check_consistency.py`, which will list everything stale.

**2. Terrain-only products.** Hydrography; boundary conditions; roughness;
albedo without lakes.

**3. Flux bracket.** Two caveats on how the slope gets quoted.

*Points must clear the extrapolated-offset criterion, not just the drift test.*
The first T21 point at 0.90 passed every drift criterion and still had +0.348 K
of approach left, with a fitted time constant matching the expected one. A point
that enters a bracket a third of a kelvin low bends the slope.

*The bracket varies two things.* Flux is moved by semimajor axis, which is
correct here -- it keeps the star fixed, so `k25v` stays valid, where moving
luminosity across a 10% flux range would shift the effective temperature by about
65 K and leave the 4900-5000 K window the spectrum was interpolated in. But the
year moves with it, 189.6 to 175.2 days across the bracket, so the measured slope
carries a small year-length effect through seasonality. Second-order for an
annual mean; state it rather than imply the bracket varied one thing.

The knobs swap at lock: semimajor axis while searching, luminosity afterwards,
once the calendar is fixed and LPJ-GUESS is compiled against it.

*Original text follows.* The flux must be re-derived, not carried: relief compresses
about 20%, which moves lithology exposure and orography, and mean land rock
albedo has separately moved from 0.3135 to 0.2810 through the cover-chain fix.
Two converged runs spanning the design band, then interpolate. Do **not**
extrapolate from one; the project has paid for that twice.

**4. Bootstrap climatology.** From the better of the two, with `--seasonal-output`
so snapshots exist. The spin-up runs without them, and building a climatology
afterwards fails on the missing snapshots -- which has already happened once.

**5. The fields that needed a climatology.** Lakes; pedology soil map; then 229,
and set `model.soil_water_source`. Expect this to matter: land runoff reads 2.8%
against Earth's ~35% under the uniform bucket, and that is the single largest
known distortion in the land hydrology.

**6. The real baseline run**, with lakes, soil water, roughness and both energy
diagnostic sets. This is the climatology every downstream product uses.

**7. Verdict, biosphere, and the loop.** Carve verdict on step 6's climatology;
LPJ-GUESS on the same; then back to the terrain if basins flip.

## Cost

Three converged climate runs before the first verdict: two for the flux bracket,
one for the real baseline. The bootstrap is one of the bracket runs rather than a
fourth, since any converged climatology will do for fields that are themselves
about to be recomputed.

## The trap this plan exists to avoid

Enabling `soil_water_source` after the baseline has run means the baseline was
computed with a uniform soil bucket, which changes evaporation, which changes
P - E, which is the numerator of the carve criterion. A verdict taken on that
climate is a verdict on the wrong evaporation. The same applies to lakes, which
darken 7% of the land. Both have to be in place *before* the run whose
climatology the verdict uses -- not merely before the verdict.
