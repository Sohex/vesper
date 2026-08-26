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

## The glacial altitude gate is gravity-invariant, and that is correct

The gate that decides which high ground glaciates is expressed in the model's
dimensionless elevation parameter. The obvious reading is that this makes the
snowline track MAXIMUM RELIEF, which goes as 1/g, rather than a freezing
altitude set by a lapse rate and a surface temperature. Write out both sides and
the reading does not survive:

- physical height is `reliefScale * f(elev)` with `reliefScale = g_ref/g`,
  because a crustal root fails at `sigma/(rho g)`. The ceiling goes as 1/g.
- the height at which a surface temperature reaches freezing is
  `(T_s - T_freeze)/Gamma`, and a dry adiabat is `Gamma = g/cp`. The freezing
  height ALSO goes as 1/g, at fixed surface temperature and composition.

The two carry one factor of `1/g` each, so their ratio is gravity-free and a
DIMENSIONLESS gate is the invariant form. Re-anchoring the gate to a fixed
number of kilometres would introduce an error rather than remove one. So the
gravity half of `PHYS-13`'s finding is not a defect; `vendor/orogen/js/glacial-ice.js`
carries the algebra and `tools/test-glacial.mjs` pins it, so a later gravity term
on the gate has to argue with it rather than be added quietly.

The cancellation is exact only for a dry adiabat at Earth's `cp` and Earth's
surface temperature, and none of those three holds here. The residual is a
fraction of the gate: this world's measured warm-season environmental lapse rate
is shallower than `g/cp` implies, so the true freezing height sits somewhat
above where the dimensionless gate puts it. The TEMPERATURE term the gate is
missing entirely is the whole of the effect beside that, which is why the fix
is an ice mask from a climatology rather than a gravity term. `PHYS-13`.

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
mesh: 0.07 to 0.33 of a cell. Four times the region count does not change the
verdict, only the fraction: the mean edge is 7.60 km on the build the pipeline
points at, which is 0.13 to 0.66 of a cell. Attaching a `g^3` velocity term to a process the
grid cannot resolve buys nothing, whatever the term's rigour. The place it can
be resolved is the downscaling pass, which needs sub-grid hypsometry persisted
anyway -- see the closing section of `notes/glacier-rough-pass.md`, which reached
the same conclusion from the other direction.

The scale reinforces it. This world's glaciers run a few times Earth's mountain
glaciation with no ice sheets at all, so the affected area is small as well as
unresolved. `GRAV-6` is therefore a declared gap deferred to downscaling, which
is Option B below rather than an outstanding task.

### The exponent is not one number, and 2.23x is the wrong one for this term

The `2.23x` above is `g^3` on ice VELOCITY at fixed thickness, and that is not
the quantity `iceFlow` stands for. Both halves have to be settled before a
factor can be written down, and neither is settled by Glen's law alone.

**Which quantity is held fixed.** Depth-averaged deformation velocity is
`u = 2A/(n+2) (rho_i g sin a)^n H^(n+1)` with `n = 3`, so `u ~ g^3` at fixed
`H`. But the flux is `q = u H ~ g^n H^(n+2)`, and `q` at a point is set by
upstream accumulation, which is a climate quantity and not a gravity one. Hold
`q` fixed instead, let `H` find its own value, and `H ~ g^(-n/(n+2))` and
`u ~ g^(n/(n+2)) = g^0.6`. `iceFlow` is an accumulation down the ice drainage
graph -- a supply, seeded from the glaciation index and summed downstream -- so
it is the FLUX reading that applies to it, not the fixed-thickness one. Orogen
carries no thickness and no mass balance, so it cannot decide this from the
inside.

**Which erosion law.** Erosion is not velocity. Abrasion goes as `u_b^l` with
`l` between 1 and 2 depending on whose law, and quarrying after Iverson (2012),
as iSOSIA implements it, goes as `K_q p_e^3 u_b (slope + 0.55)^2` where the
effective pressure `p_e ~ rho_i g H` carries its own gravity. Both are on the
SLIDING velocity, which is a Weertman-type law Orogen does not have either.

At `g = 1.30626` Earth the bracket that follows is:

| reading | exponent | factor |
| --- | --- | --- |
| E ~ u, flux fixed (the one `iceFlow` is) | g^0.60 | 1.17 |
| E ~ u^2, flux fixed | g^1.20 | 1.38 |
| quarrying, p_e^3 u, flux fixed | g^1.80 | 1.62 |
| E ~ u, thickness fixed (the `2.23x` quoted) | g^3.00 | 2.23 |
| E ~ u^2 or quarrying, thickness fixed | g^6.00 | 4.97 |

The bracket spans a factor of 4.2 and its width is dominated by a STRUCTURAL
choice, not by a coefficient: whether thickness or flux is the conserved
quantity. That is the finding. A `g^3` on `iceFlow^0.6` would be neither end of
it, and it would not even be dimensionally the same object, since `0.6` is an
exponent on an accumulated supply and `n` is an exponent on a stress.

This does not reopen the row and it does not change the recommendation. It
sharpens it: the reason not to bolt a gravity term onto `iceFlow` is not only
that the quantity is heuristic, it is that the model does not carry the state --
thickness, sliding, effective pressure -- that decides which exponent the term
would take. The bracket is what to carry until it does.

### ExoPlaSim's `dsnowz` is not the missing thickness, and the missing state is a solve

The state the bracket needs is not held anywhere else in the pipeline either.
`glaciermod.f90` carries `dsnowz`, the snow and ice column in metres of liquid
water equivalent, and it is prognostic and genuinely uncapped: `landmod`'s
`DSMAX` hard-clip is off, `config/planet.yaml` declaring
`surface.glaciers.max_snow_depth_m: -1.0`. The question is whether that column is
an admissible `H` for a Glen's-law erosion term. It is not, for four reasons,
and each of them is a property of the two codebases rather than of a run.

**It solves a different equation.** `dsnowz` integrates the local surface mass
balance and nothing else, `dH/dt = b`. A glacier's thickness solves continuity
with a flux divergence, `dH/dt = b - div q` with `q = u H` under Glen's law, and
the two agree only where `div q` vanishes. `glaciermod`'s own header says as
much: it is a coupling slot for an external gridpoint glacier model, which is
expected to SET the height. `exoplasim/notes/parameter-decisions.md` states the
consequence in the model's terms -- without flow the accumulation zone thickens
into a tower rather than spreading into a sheet.

**Where a glacier erodes, `dsnowz` is zero.** Quarrying and abrasion are carried
by the trunk and by the ablation zone, and ice is there only because it flowed in
from above. `dsnowz` is zero by construction wherever the local balance is
negative, and largest at the accumulation summit, where a real ice sheet is
nearest to stagnant and coldest at the bed. So it is not a noisy estimate of the
right field; it is anti-correlated with it in space.

**It cannot decide the exponent, which is the whole width of the bracket.** The
bracket above is wide because of a structural choice: thickness conserved gives
`u ~ g^3`, flux conserved gives `u ~ g^0.6`. Substituting a thickness that
carries no ice dynamics does not MEASURE that choice, it imposes the
thickness-fixed limb, because `dsnowz` has no mechanism by which faster flow can
thin the column. The bracket would appear to collapse to the erosion-rate
exponent while having actually been settled by an artefact of the coupling. That
is a worse position than the honest bracket, not a better one.

**It is a cell mean against a law that is quartic in thickness.**
`glaciermod.f90` records it on the line that flags a glacier: `dsnowz` there is
the mean over the cell, while the 30 m it is compared against is a column
property of ice. Glen's law makes `u ~ H^(n+1) = H^4` and the flux `~ H^5`, so
the mean of the fourth power is not the fourth power of the mean, and the gap is
set by the sub-grid distribution rather than being a constant. Orogen's mesh is
finer than a T42 cell by more than an order of magnitude, so one cell mean would
be evaluated inside a quartic for hundreds of mesh regions.

**AND THE INSTRUMENT IS NOWHERE NEAR THE EFFECT.** MEASURED 2026-08-26 from the
raw `newdsnow` dump `glacierstop` writes, across all six T21 runs on disk, every
one with the glacier module enabled and the snow cap lifted, at 25 to 50 orbits:

| quantity | across the six runs |
| --- | --- |
| deepest `dsnowz` | 0.122 to 0.211 m water equivalent |
| the same as ice at `rhoglac` 850 | 0.14 to 0.25 m |
| cells reaching `GLACELIM`, 2.0 m | 0 of 2,048, every run |
| cells reaching the module's 30 m flow threshold | 0 of 2,048, every run |
| cells with the glacier persistence clock complete | 0, every run |

Glen's law would be evaluated more than two orders of magnitude below the
thickness at which ice flows at all. The runs are short and none is a baseline,
so this is a floor on what exists rather than an equilibrium, but run length is
not the binding reason: `phys-13` measured the cause, which is that the
model evaluates its own high ground about 7.8 K too warm on average and 21.9 K
in the top tenth from the sub-grid orography excess, measured at T42, and no ice
survives that at any truncation the model can afford. These runs are T21, which
is coarser, so the excess there is larger rather than smaller.

**What would be admissible, and what it costs.** The erosion terms want a basal
sliding velocity and, for quarrying, an effective pressure. Both are OUTPUTS of
an ice-dynamics solve: continuity with `q = u H` under Glen's law, a
Weertman-type sliding law, and a subglacial hydrology assumption. There is no
field ExoPlaSim can hand across that substitutes for that solve, because the
missing state is the solve and not one of its variables.

What ExoPlaSim can legitimately supply is the other side of the same solve. A
shallow-ice model needs a surface mass balance `b(x)` as an INPUT, and that is a
climate quantity. This is the strong direction rather than the weak one: the
published reference implementation this row already names drives its ablation
from 2 m air temperature alone, with no shortwave, albedo or spectrum, which is
structurally blind to the term this world differs on most over ice, while
ExoPlaSim computes melt with all three under this star's measured spectrum. So
the half of a shallow-ice model that is the other ninety per cent is the half
this project is in a position to provide, and the half it would have to acquire
is the flow algebra, whose two routes `clim-62` costs.

**What this changes and what it does not.** The bracket stands at its full
width, and so does the sub-grid argument: a 1 to 5 km glacial valley against the
7.60 km mesh of the build the pipeline currently points at. Option B stands.
What changes is that the row's blocker is now named exactly rather than left as
"the model does not carry a thickness": a thickness field would not have
unblocked it, and the thing that would is a solver.

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
be justified on its own merits rather than smuggled in as a gravity fix. The
exponent bracket above is what makes that concrete: the model would have to grow
a thickness and a sliding law before it could say which power of `g` to use.

The one line is to record `n = 1` at the erosion call site, because the entire
correctness argument rests on it and nothing says so.

## What is not affected

- Land/sea area split: scaling is monotonic and sea level sits at elevation 0, so
  the mask does not move.
- Ocean volume and sea level: bathymetry is unscaled and the basin is unchanged,
  so the water budget stays consistent.
- Basin ids and carve verdicts across a gravity change: the catalogue is
  bit-identical, so an existing verdict replays. Already in `CLAUDE.md`.
