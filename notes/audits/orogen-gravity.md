# Audit: how gravity enters Orogen

*Audited 2026-08-16, against `vendor/orogen` at `cf-fork` 3907c58. Tasks in
the `bd` issue tracker under the `grav` label.*

Two questions were put: whether gravity being applied as a post-hoc scaling
factor rather than entering generation warrants fixing, and whether the unscaled
bathymetric side is a bug. Neither is a defect. The one real gravity gap in the
generator is somewhere else entirely.

## What the code does

`planet-params.js` sets `reliefScale = REFERENCE_GRAVITY_MS2 / gravityMS2`,
exactly 1 at Earth gravity. It is applied at the model-unit-to-kilometres
conversion and nowhere else. Tectonics, erosion and basin detection all run in
dimensionless model units, so two builds differing only in gravity are
bit-identical in every hash.

The rule is: **positive heights scale, negative ones do not.**

## Bathymetry is correctly left unscaled

Ocean depth is an isostatic balance between the water-plus-oceanic-crust column
and the continental column. Write the mass balance out and g multiplies every
term, so it cancels: ocean depth is set by density contrasts and crustal
thicknesses, and a planet at 2 g has the same ocean depth as one at 1 g. Scaling
it would be wrong.

The defect was not the physics but its distribution. The same rule was written
out at three call sites -- `js/basins.js`, `js/data-export.js`,
`tools/export-maps.mjs` -- and only one carried the reasoning. The site computing
the exported `elevation_km` had no comment on the sign asymmetry at all, where a
bare `h > 0 ? h * reliefScale : h` reads like a sign bug.

Fixed: `scaledHeightKm()` in `js/color-map.js` is the single definition, carries
the argument, and all three sites call it. All 106 fork tests pass.

## The land scaling is exact, not approximate

The intuition that gravity should affect erosion and slopes is right in general
and does not apply here, for two reasons that have to be checked against the code
rather than reasoned about in the abstract.

**This mesh cannot represent a hillslope.** Cells are 15.19 km apart.
`terrain-post.js` records that a genuine 34-degree angle of repose would need a
model slope of about 4300, while the median downhill neighbour pair sits at
0.037 degrees, the 99th percentile at 0.37, and the steepest pair anywhere at
16.9. `talusSlope` is a near-inert floor already cleared by about 83% of slopes,
so the thermal step behaves as linear diffusion. The angle of repose is indeed
gravity-independent -- g cancels in both the driving stress `rho g h sin(theta)`
and the resisting `mu rho g h cos(theta)`, which is why repose angles on Mars
resemble Earth's -- but repose angles are two orders of magnitude below what this
grid resolves. Every slope the model carries is a channel or landscape slope.

**The erosion law has n = 1, and at n = 1 the post-hoc scaling is exact.** The
hydraulic step solves

    h_new = (h + factor * h_receiver) / (1 + factor),
    factor = K * flow^m * dt / dx

which is the implicit form of `dh/dt = -K A^m S` with slope exponent one. At
steady state `K A^m S = U` gives `S ~ 1/K`, and stream power puts `K ~ rho g`, so
`S ~ 1/g`. Scaling every height by `1/gamma` reproduces that exactly. **For the
fluvial landscape the shortcut and the correct calculation are the same
operation.**

How near to steady state, on this build's settings -- hydraulic slider 0.6, so
K = 0.00036, m = 0.5, dt = 1, 12 iterations:

| upstream cells | factor/iteration | residual relief after 12 iterations |
| --- | --- | --- |
| 1 | 0.181 | 0.136 |
| 10 | 0.573 | 0.004 |
| 100 | 1.812 | 0.000 |
| 1000+ | 5.7+ | 0.000 |

Anything with ten or more upstream cells is at grade. Running erosion at the
correct gravity instead of scaling afterwards moves single-cell headwaters from
13.6% residual relief to 7.8% and changes nothing else measurably: a **5.7% of
local relief** difference, confined to the finest headwaters.

The live uncertainty is not the method but the exponent. **n = 1 is a choice**,
and nothing in the generator states it. At n = 2 the correct scaling would be
`g^(-1/2)` = 0.875 rather than `g^(-1)` = 0.766, so the current behaviour would
be over-correcting relief by about 14%. That single number is the whole remaining
question. `GRAV-4`, `GRAV-5`.

A smaller point, unmodelled and unmodellable here: the 1/g argument is about
**strength-supported** relief. Isostatically compensated land -- a plateau on
thick crust -- is gravity-independent for the ocean's own reason. Separating the
two needs a crustal-thickness field the model does not carry, so uniform 1/g on
land slightly over-suppresses plateaus.

## The real gravity gap is glacial

Glacial erosion runs at slider 0.8, eight iterations, as
`deepening ~ iceFlow^0.6`, and **carries no gravity term at all**. This is the
one process where gravity enters strongly and non-linearly: Glen's flow law makes
ice velocity go as the cube of driving stress, `u ~ (rho g h sin a)^3`, which is
**2.23x** at this world's gravity. A high-gravity world should show markedly more
glacial erosion -- deeper troughs, more overdeepening, and by extension different
closed-basin geometry.

The catch is that `iceFlow` is a heuristic accumulation rather than a Glen's-law
velocity, so multiplying it by `g^3` would attach a rigorous factor to an
imprecise quantity and read as more trustworthy than it is. Fixing this properly
means rebuilding the glacial model, which is its own project.

**And it belongs downstream, not here.** A glacial valley is 1 to 5 km wide
against a 15.19 km mesh cell, so glacial erosion is a SUB-GRID process on this
mesh: 0.07 to 0.33 of a cell. Attaching a `g^3` velocity term to a process the
grid cannot resolve buys nothing, whatever the term's rigour. The place it can
be resolved is the downscaling pass, which needs sub-grid hypsometry persisted
anyway -- see the closing section of `notes/glacier-rough-pass.md`, which reached
the same conclusion from the other direction.

The scale reinforces it. This world's glaciers run a few times Earth's mountain
glaciation with no ice sheets at all, so the affected area is small as well as
unresolved. `GRAV-6` is therefore a declared gap deferred to downscaling, which
is Option B below rather than an outstanding task.

## SWOT: make gravity physical, or roll with what we have

### Option A -- address it properly

Give the erosion laws physical units: a real K containing rho and g, an explicit
n, and a Glen's-law ice velocity.

**Strengths.** Removes the only unmodelled gravity dependence, glacial erosion,
a 2.23x effect landing on the feature this world is most interested in. Forces n
to be chosen rather than inherited from whatever the implicit solver supported.
Makes the generator defensible at any gravity instead of calibrated near Earth's.

**Weaknesses.** Expensive and invalidating: every build regenerates, every carve
verdict replays against different terrain, and the basin catalogue stops being
bit-identical across a gravity change -- currently a documented feature that lets
verdicts survive one. The fluvial half of the work provably buys nothing. And
Orogen has no sediment budget, no timescale and no real uplift rate, so
"physical units" would be real K over heuristic everything else.

**Opportunities.** The glacial model is the weakest part of the generator
independent of gravity; this would be a reason to rebuild it, improving the world
at Earth gravity too. It would also let the project state a defensible n, which
nothing currently does.

**Threats.** High risk of precision theatre -- a rigorous factor multiplying a
heuristic reads as more trustworthy than it is, which is `docs/src/practice/failure-modes.md`
class 9 wearing different clothes. Recalibrating the sliders against a changed
erosion law is open-ended, and prior worlds would no longer reproduce, breaking
the planet code `01eshm059lt0b9mpgro2y83t` as an identity.

### Option B -- roll with what we have

Keep the post-hoc 1/g amplitude scaling, state n = 1 as the assumption it is, and
carry glacial erosion as a declared gap.

**Strengths.** Provably correct for the fluvial landscape, which is essentially
all of the resolved terrain. Costs nothing. Preserves bit-identical basin
catalogues across gravity, so carve verdicts replay. Every existing build stays
valid and datable.

**Weaknesses.** Glacial overdeepening is under-modelled at high gravity on a
world that has glaciers as a deliberate feature. n = 1 is unstated and
unexamined; if it should be 2, relief is over-suppressed by about 14%.

**Opportunities.** The n-uncertainty is a single number and can be bracketed
cheaply: generate one build at `g^(-1/2)` instead of `g^(-1)` and compare basin
statistics. A build is one pass, so this answers the only live question without
touching the generator.

**Threats.** The correctness argument depends on the network being at grade,
which depends on the hydraulic slider. At a much lower setting the landscape
would be transient and the equivalence would weaken. Anything that moves that
slider silently invalidates this conclusion, and nothing currently warns of it.

### Recommendation

**Option B, plus one line of code.** The fluvial case is settled and Option A
cannot improve it. The glacial gap is real but cannot be closed by adding a
gravity term to a heuristic; it needs the glacial model rebuilt, and that should
be justified on its own merits rather than smuggled in as a gravity fix.

The one line is to record `n = 1` at the erosion call site, because the entire
correctness argument rests on it and nothing says so.

## What is not affected

- Land/sea area split: scaling is monotonic and sea level sits at elevation 0, so
  the mask does not move.
- Ocean volume and sea level: bathymetry is unscaled and the basin is unchanged,
  so the water budget stays consistent.
- Basin ids and carve verdicts across a gravity change: the catalogue is
  bit-identical, so an existing verdict replays. Already in `CLAUDE.md`.
