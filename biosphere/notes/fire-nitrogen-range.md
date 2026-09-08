# The simulated fire operator on the accepted run: what it is worth, and whether its nitrogen flux is possible

Worldbuilding. Vesper is an invented planet; this note is about what its
simulated vegetation model's fire operator is doing on the one accepted run --
how much of the carbon cycle it moves, and whether the nitrogen it releases is
inside two conservation bounds. Nothing here is about fire on Earth except where
Earth is named as a distance to report.

Two questions in one measurement session on one run, and they are in this order
because the second only matters if the first says fire matters.

FIRE-8 exists because `assess_lpj_run.py` checks that values are finite,
physical and closing, and asks nothing about whether a flux is ORDINARY. A large
positive emission is finite, is positive, and CLOSES -- the nitrogen it removes
is nitrogen the model had -- so the row's own argument is that nothing in the
tree currently looks. This note is the looking, for the one axis that can be
settled without choosing a fire model.

## The criteria, and why they are criteria

Both are stated before the comparison and neither is fitted to the distribution
it judges. CLAUDE.md's rule is that a threshold chosen after the run it judges is
not a threshold; the tail below had already been measured when this was written,
so the bounds are derived from the INPUT side, where the observed distribution
cannot reach them.

**Bound 1, the stock bound.** In one simulated year a gridcell cannot lose more
nitrogen to fire than the nitrogen it holds. Fire moves nitrogen out of the
vegetation, litter and soil pools, so a year in which the fire flux exceeds the
cell's total nitrogen pool is not a large fire; it is a defect.

**Bound 2, the rate bound.** Over the retained record a gridcell is at
equilibrium, so its nitrogen losses equal its nitrogen inputs. Fire is one of at
least three loss pathways -- fire gases, soil gases and leaching -- so the
retained-record mean fire loss must be strictly below the retained-record mean
input. A cell above it is either not at equilibrium or is losing nitrogen it
never received.

Neither bound involves a fire model, a burned area or an Earth comparison. Both
have a right answer, which is what makes them tests rather than descriptions.

## What was measured

Measured on 2026-09-07, on `lpj_1e6a2b9ca51a4eff9592992cad96677b` -- the first
accepted LPJ run, 8600 spin-up, 1253 retained, 1617 gridcells,
`canonical-10m-carve2`, `npatch 5`. All 2,026,101 gridcell-years of the retained
record, no subsampling.

The fire flux is the sum of the four fire columns of `ngases.out`: `NH3_fire`,
`NOx_fire`, `N2O_fire` and `N2_fire`. The inputs are `NH4dep`, `NO3dep`, `fix`
and `fert` from `nflux.out`, which that file writes NEGATED. The stock is
`npool.out`'s `Total`.

**A unit trap sits between those three files and it is not documented anywhere a
reader would meet it.** `nflux.out` and `ngases.out` are multiplied by
`M2_PER_HA` at `commonoutput.cpp:1983` and `:2085` and are in kgN/ha;
`npool.out` at `:2057` is not, and is in kgN/m2. The comment at `:196` reading
"kgN/m2, converted by 1e4 to kgN/ha" is about the PRECISION the column is given
and not about the column being converted. Comparing the two without the factor
of 10^4 makes every cell-year fail bound 1 by about three orders of magnitude,
which is a tidy-looking wrong answer of exactly the kind CLAUDE.md's
check-the-instrument rule is about.

## Results

**Bound 1 holds, with room.**

    cell-years compared              2,026,101
    cell-years where fire > pool     0
    worst fire/pool ratio            0.0712
    at                               lon -123.75, lat -2.77, year 9118
    that year                        831.19 kgN/ha burned from a pool of 11,675 kgN/ha

**Bound 2 holds, and one cell approaches it.**

    gridcells compared               1,617
    cells with no nitrogen input     0
    cells whose mean fire loss exceeds mean input   0
    mean fire/input ratio            0.185
    worst ratio                      0.981  (fire 9.813, input 9.999 kgN/ha/yr)
    largest mean fire loss           10.494 kgN/ha/yr against an input of 12.809

## What this settles

**The tail is a transient, not a rate.** The 831 kgN/ha in the worst
gridcell-year is 7.1 per cent of that cell's nitrogen pool, released in one year
by one fire after a long accumulation. Over the retained record the same cells
lose, on average, 18.5 per cent of their nitrogen input to fire. Those two
numbers are consistent with a fire regime of long return intervals and large
accumulated fuel, and they are not consistent with a flux the model is
manufacturing.

So the measured tail is NOT out of physical range, and FIRE-8's
"reject impossible ranges" does not reject it. That is a real answer rather than
an absence of one: two conservation bounds were available, both were applied to
every gridcell-year, and both pass.

**A cross-check fell out and it confirms the emission partition.** The figure
recorded on FIRE-8 is 797.1 kgN/ha for the worst cell-year, against 831.19 here.
The ratio is 0.959, which is exactly `NOx_FIRERATIO + N2_FIRERATIO` --
`0.237 + 0.722` from `guess.cpp:60-65`, after Delmas et al. (1995). The earlier
figure summed two species of four; the four are in their declared proportions to
five digits, so the partition is intact and cannot be what produced the
magnitude.

## What this does not settle, and where it goes instead

Whether the fire REGIME is right. Both bounds constrain nitrogen and neither
constrains how often a cell burns or how much of it burns, and a model with too
few, too large fires satisfies both exactly as well as a correct one does. The
quantities that would settle that are the fire return interval and the burned
fraction, and there is no observed burned area on this world to compare them
against -- which is why FIRE-9's matched model-form comparison, not a range
check, is where that question belongs.

Nor does it say the emission is right in ABSOLUTE terms. Earth fire emission from
burned area runs of order 1 to 10 kgN/ha/yr, and this run's mean of 0.98 sits
inside that while its per-cell means run to 10.5. Under CLAUDE.md's rule Earth's
figure is a distance to report and never a target to solve onto, and the distance
here is unremarkable.

Both bounds are worth keeping as acceptance criteria precisely because they
passed: they cost one pass over two tables, they have right answers, and the
defect they would catch -- a fire operator drawing nitrogen the model never had --
is one that closure testing cannot see, because a flux that removes nitrogen the
model had closes whatever its size.

## What fire is worth in this world

Measured on the same run, over the same 2,026,101 gridcell-years, from
`cflux.out`. That table is in kgC/m2/yr and is written without the `M2_PER_HA`
factor the nitrogen flux tables carry. `Veg` is the NEGATED net primary
production (`commonoutput.cpp:1488` reports `-NPP`), `Soil` is heterotrophic
respiration plus organic leaching, and `Fire` is the fire carbon flux. At
equilibrium the net is zero, so production is balanced by soil respiration plus
fire, and fire's share of the carbon TURNOVER is `Fire / (Soil + Fire)`.

    NPP                        0.15309 kgC/m2/yr
    soil respiration + leach   0.12374
    fire                       0.01639
    fire / (soil + fire)       0.1169
    fire / NPP                 0.1070

**Fire moves about a tenth of this world's simulated carbon turnover**, and it
does so through the operator that is actually running: GLOBFIRM, on the
configuration `run_lpj_guess.py` writes.

It is not spread evenly, and the distribution is the more useful half:

    per-gridcell fire share of turnover, 1617 cells
      median                 0.0029
      p75                    0.0921
      p90                    0.1577
      p99                    0.3206
      maximum                0.3595
      cells above 10%          369   (22.8%)
      cells above 25%           52   (3.2%)
      cells with no fire       669   (41.4%)

So two fifths of the land never burns and roughly a quarter of it routes more
than a tenth of its carbon through fire, with the most fire-dominated gridcells
approaching a third.

**This is what makes the rest of the fire work worth doing, and it is worth
stating as a number rather than as an assumption.** A ten per cent term in the
carbon cycle is far above any noise floor, so every question about the fire
operator -- the deleted burn-probability floor, the Earth-fitted curve it runs,
the four defects registered against the effects layer it does not yet run -- is
a question about a material quantity rather than a rounding. It also means
FIRE-9's comparison is well posed: a no-fire arm would move about a tenth of the
carbon turnover, which no ensemble spread is going to hide.

It says nothing about whether a tenth is the RIGHT share. Earth's terrestrial
fire flux is of the same order as a fraction of NPP, and under CLAUDE.md's rule
that is a distance to report rather than a target to solve onto. What the number
establishes is the STAKE, not the answer.

## Where the deleted burn-probability floor would have bound

`biosphere/config/fire.yaml` records the deletion of GlobFIRM's
`fireprob = 0.001` floor and, until this measurement, said the number of
gridcells it bound on "cannot be stated here. It needs a run". It does not.

The floor and the reporting cap are THE SAME CONSTANT, which is what makes the
extent readable from a run that has no floor in it. `commonoutput.cpp` writes a
fire return interval of exactly 1000 years for any patch whose `fireprob` falls
below 0.001, and the deleted line raised any such `fireprob` to 0.001. So a
capped cell-year in `firert.out` is precisely a cell-year the floor would have
bound, and the accepted run reports one without carrying the other.

    cell-years                        2,026,101
    at the 1000-year cap              1,209,759   59.71%
    gridcells at the cap in EVERY year      791   48.9% of 1617
    gridcells at the cap in some year      1183   73.2%
    gridcells never at the cap              434   26.8%

    mean burned fraction where capped     1.65e-10
    mean burned fraction where not        1.90e-02

**The floor was not rounding a small number up.** Where it bound, the model's
own burned fraction is 1.65e-10 -- approximately zero, to eight orders of
magnitude below the cells that burn -- and the floor would have replaced it with
1e-3. It did that on 59.71 per cent of gridcell-years and in every year of the
record on 791 gridcells, which is about half this world's simulated land.

The consequence follows from `fire.yaml`'s own arithmetic once the spin-up is
the derived one. A 0.001 per model-year hazard of a stand-replacing fire over
`nyear_spinup 8600` is 8.6 expected fires per patch and a 99.98 per cent chance
of at least one. On the gridcells it reaches, the floor does not perturb the
cohort age structure; it replaces it.

**What this does not settle**, and it is the part that still needs the paired
run `comparison_arm` specifies: what that did to the vegetation. The extent is
readable because the deleted constant is also a diagnostic the run reports; the
standing biomass, cohort age structure and soil carbon of those gridcells under
a restored floor are not, because the accepted run has no floor in it to remove.

No LPJ-GUESS run was performed for this note. It reads the tables of a run that
already existed.
