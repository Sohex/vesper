# Orogen has no time axis, and what to do about it

Written 2026-08-17 because this keeps being rediscovered and re-litigated, once
per component that needs a duration. It is settled. Read this instead of working
it out again.

## The fact

**World Orogen has no time axis at all.** Not a coarse one, not an implicit one.
Its own lithology module states it: "no ages, no stratigraphy and no
unconformities, because there are no timesteps to hang them on". The export
manifest carries no duration anywhere -- `params.hydraulicErosion`,
`thermalErosion` and `glacialErosion` are dimensionless intensity sliders, and
`vendor/orogen/tools/README.md` says of basin infill that the generator has "no
sediment budget or timescale, so modelling infill is left to whoever has one".

The consequence that matters: **a duration is UNDEFINED here, not unmeasured.**
Asking the generator how old a surface is, or how long a landscape has been
relaxing, asks it for a concept it does not have. There is no field to add and no
setting to read. Anything that wants a duration has to supply one, and this
document is about how.

`notes/audits/orogen-gravity.md` reaches the same conclusion from the other
direction: putting real physical units on the erosion law would give "real K over
heuristic everything else", because there is no sediment budget, no timescale and
no real uplift rate for the K to act through. That is `docs/src/practice/failure-modes.md`
class 9 -- a rigorous factor multiplying a heuristic reads as more trustworthy
than it is.

## What to do instead: an expected-value argument over a stationary population

The move that works, and that this project has now used three times, is to stop
asking "when" and ask **"what does a randomly chosen moment look like?"** A
landscape, a volcanic province and a lake population are all stationary
populations with a turnover; the terrain is one sample from such a population;
so the question has an answer even though the duration does not.

Earth is the calibration, because Earth is also one randomly chosen moment.

**Volcanic provinces, in `pedology/config/pedogenesis.yaml`.** Andic soils need
ongoing ash resupply, so they need to know whether a province is active. Rather
than ask for an eruption age: a flood basalt emplaces most of its volume in about
a million years and then persists for hundreds, so a randomly chosen moment finds
a given province erupting with probability of order one percent. Expected value
says these provinces are old, decisively, and `flood_basalt` and `oib` are
excluded from andic soils as the correct default rather than a hedge. The arc is
the exception because it is not a province but a boundary process, active for as
long as the plate is subducting -- which is the whole time the terrain exists.

**Deposit types, in `docs/src/reference/economic-minerals.md`.** Anything whose defining
control is an absolute age cannot be placed, only its tectonic setting can. The
concrete loss is komatiite-hosted nickel, which Naldrett puts at 2.7 to 1.9
Ga. That is stated as something the world cannot have rather than approximated.

**The carve verdict, in `hydrography/notes/retain-fraction.md`.** How much of a
basin's rim survives its overflow depends on the relaxation window. Rather than
declare one: match the DENSITY of standing through-flowing impounded basins
against Earth's, at the size this mesh can resolve. That fixes the incision
coefficient at 161 m per (m3/s)^0.5, with a bracket from the Poisson error on
Earth's small sample.

## The two ways this goes wrong

**Declaring a duration and then reasoning as though it were measured.** The
window is not 1e5 years or 1e6 years; there is no fact of the matter. A number
declared here should be labelled as calibrated to an outcome, not as a duration.

**Calibrating against Earth's survivors instead of Earth's density.** Survivors
are not a random sample of anything: Earth's large standing lakes are almost all
tectonically maintained rift basins, so the edge of that population measures
subsidence and youth rather than incision. The density of survivors per unit land
is the robust statistic; the identity of any individual survivor is not. This one
cost real work -- see `hydrography/notes/retain-fraction.md`, where reading the
survivor edge gave an answer 300 times away from the density answer, and the
size-class mismatch between Earth's small lakes and this mesh's large basins
flipped the sign of the conclusion once on top of that.
