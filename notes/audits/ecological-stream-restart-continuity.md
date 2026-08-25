# Is a segmented run the same experiment as the whole run?

**Recorded:** 2026-08-24
**Scope:** the Vesper climate model's new ecological output stream (EFOR-2) and,
through it, whether the model's own state reproduces across a segment boundary.
Measured, not reasoned: the runs below were executed.

Worldbuilding frame: this is a correctness check on a simulation of the planet
Vesper. Nothing here is about the real world.

## Result

Two findings, and they are independent.

**The ecological stream's structure is exact across a restart.** A run split in
two writes the same number of intervals as the same run taken whole, over
intervals that agree on their start, end and duration in absolute seconds, with
the same record headers. The interval that STRADDLES the segment boundary is
closed at the right step and covers the span it declares. That is the property
`naccueco` and the nineteen `aeco*` accumulators are serialized for, and it is
the property a replay protocol later depends on.

**The model's own state did not reproduce across a restart, and the divergence
was immediate rather than accumulated.** With the stream reporting one interval
per timestep, the first three intervals of a six-step run are bit-identical to
the same six steps taken as 3 + 3; the FOURTH interval, the first timestep after
the restart, already differs by 1.19 K in near-surface air temperature and
1341 Pa in surface pressure. The final restarts of the two runs differ in 60 of
222 records.

The second finding is not caused by the first. With the stream switched off the
two restarts differ by exactly the same pair of hashes, so it is a property of
the model rather than of anything the stream does. It was tracked as
`world-8yyh`, and the section below has its cause, which is one unguarded line
in `glacierini`, and the demonstration that guarding it makes a segmented run
bit-identical to the whole run in every restart record.

## What was run

T21, ten levels, four threads, on the surface fields and namelists of
`exoplasim/runs/run_2b20e3324bb0`, cold started from `TSST_EQ = 300.0` and
`TSST_POL = 273.0` in `icemod_nl`. Three preparations were needed and each is
worth recording, because the bed that exists in the tree does not run as it
stands:

- `plasim_namelist` in that run and in `exoplasim/bench/bed_check` carries
  `NGUIDBG`, and `carbonmod_namelist` carries `ZETA`. Neither is declared by the
  current sources, and gfortran stops the model on the first one it reaches.
- No restart in `exoplasim/runs/` can be resumed by the current model:
  `glaciermod.f90` reads `persistt` under the default `nexcheck`, and no
  restart in the tree carries it. Every continuation test therefore has to make
  its own restart first.
- A cold start refuses without a sea surface temperature climatology or a
  declared `tsst_eq`/`tsst_pol` pair (`icemod.f90:ice_cold_start`), and
  `bed_check` has neither.

The arms:

    whole      N_RUN_STEPS = 6,  NECO = 1, NECOSTEP = 1, cold
    segment 1  N_RUN_STEPS = 3,  same bed, cold
    segment 2  N_RUN_STEPS = 3,  segment 1's plasim_status as plasim_restart

and the comparison is
`exoplasim/scripts/compare_eco_streams.py whole/plasim_eco seg1/plasim_eco seg2/plasim_eco`,
which separates a structural failure, a bounds failure and a payload failure
because they have different owners.

The same pair was run at `NECOSTEP = 8` over 24 steps split 12 + 12, so that the
boundary falls four steps INSIDE an interval rather than on one. Structure and
bounds hold there too, which is the case that matters: an interval closed by
step count alone would have closed four steps late and declared a span it did
not cover.

That pair is also what caught the interval ORIGIN. `ecogp` first wrote the
bounds from `nstep`, which at that point in the loop is the index of the step
just integrated rather than of the next one, so every block was labelled one
timestep early. The blocks still tiled without gap or overlap and their
durations were right, so every closure check EFOR-7 would run passed on them.
The bounds now read [0, 21600], [21600, 43200], [43200, 64800] seconds for the
first three intervals of a cold run at eight steps of 2700 s, which is the
tiling from absolute zero that they should always have been.

## The measurements

Per-timestep divergence, 6 steps whole against 3 + 3, worst cell:

| interval | near-surface air temperature | surface pressure | total precipitation |
| --- | --- | --- | --- |
| 1 | 0 | 0 | 0 |
| 2 | 0 | 0 | 0 |
| 3 | 0 | 0 | 0 |
| 4 | 1.190 K | 1341 Pa | 0 |
| 5 | 0.692 K | 126 Pa | 0 |
| 6 | 2.654 K | 1052 Pa | 0 |

Intervals 1 to 3 are the pre-restart steps and interval 4 is the first
post-restart step. Precipitation stays identical throughout, which is a
consequence of a six-step run producing almost none rather than evidence that
the moist path is unaffected; the 24-step pair does show precipitation
diverging.

Restart records differing after 24 steps taken as 12 + 12, from
`exoplasim/scripts/diff_restarts.py`: 60 of 222, worst relative difference 1.0
on `dcc` (cloud cover), and the leapfrog previous-level arrays `sdm` and `spm`
differing in 4830 of 5060 and 483 of 506 elements.

Both runs reproduce themselves exactly when repeated, so this is not a threading
race and the determinism gate is not implicated.

## The cause, measured 2026-08-24

**`glacierini` overwrote the leapfrog minus level on the restart path.** It ends
each of its two branches with `call mpscsp(sp,spm,1)`, which scatters the
current surface pressure into `spm`, the t - dt level. On a cold start that is
right: a run beginning from rest has the same state at both levels. On a restart
it is not, and neither call was guarded. `surfini` is called from `prolog`
BELOW the restart read, and it calls `glacierini`, so the sequence was
`read_atmos_restart` restoring the true `spm` and `glacierini` replacing it with
`sp` a few hundred lines later. The equivalent scatter in `surfini` itself has
carried an `if (nrestart == 0)` guard all along; the two in `glaciermod.f90` did
not.

The instrument that named it is a restart round trip: take N steps, then restart
from that state and take NO steps. Nothing integrates between the two
`plasim_status` files, so every record that differs is state the restart does not
carry. There were two, and only one of them was live:

| record | what it is | relative difference | elements |
| --- | --- | ---: | ---: |
| `spm` | leapfrog log surface pressure, t - dt | 1.256 | 483/506 |
| `dglac` | glacier flag | 1.0 | 1019/2048 |

`spm` after the round trip is `sp` before it, exactly and element for element,
which is the signature of the scatter rather than of an integration. `dglac`
differs because `glacierini` recomputes the flag from the snow depth while
`glacierstep` maintains it, and it does not affect the comparison: both arms of
a segmented run take at least one step and end with the maintained value.

**The demonstration.** With the two calls guarded and nothing else changed, a run
split into two segments is bit-identical to the same run taken whole in ALL 222
restart records, at 2 = 1 + 1, at 6 = 3 + 3 and at 24 = 12 + 12, and the
round trip loses nothing. Unguarded, the same bed differs in 46 of 222 records by
the first post-restart step. Both `NGLACIER` settings were run, because the two
calls sit in the two branches of one switch and only one of them is compiled out
by a namelist.

**This was never a property of any tree.** The same round trip run on `main`, 238
commits below the branch that first measured the divergence, loses `spm` in the
same way and by the same amount, and the same guard fixes it there. What the
project had proved before this is a different property:
`notes/audits/model-reproducibility.md` establishes that the model gives the same
bytes twice from ONE restart, which a model can satisfy while integrating a
different weather either side of every boundary, and this one did.

**A cold-start comparison needs a declared `SEED`.** With `KICK > 0` and no
`SEED` in `plasim_namelist`, `initrandom` draws the initial perturbation from
`system_clock`, so two cold starts are two different experiments and any
comparison between them measures the clock. The bed used here declares one, and
`verify_restart_continuity.py` refuses a bed that does not.

## Two things the experiment also established

**The stream does not perturb the model.** The restarts written with `NECO = 1`
and with `NECO = 0` are byte-identical, at both split points. The stream reads
state that already exists and writes to a unit of its own, so it is a pure
diagnostic and turning it on cannot change a result.

**`restartmod`'s record limit was already below the model's record count.**
`restart_ini` stops the model when a restart file holds `nresdim` records, and
`nresdim` was 200 while the call sites could emit 215 distinct names; the
ecological stream's twenty-one took the total to 236 and turned the latent
defect into a hard stop, reported as "Too many variables in restart file" at the
START of the second segment, after the first has finished and written a file
that looks complete. The constant is now 512 and
`restart_schema.check_record_limit` holds it against the inventory, so the bound
is derived rather than carried.

## What this means for the forcing artifact

The ecological forcing artifact EFOR-3 builds from this stream can be assembled
across segment boundaries without a seam in its INTERVALS: no gap, no overlap,
no interval that covers less than it declares. That is what EFOR-7's acceptance
checks test and it holds.

The weather either side of a boundary now holds too: with the `glacierini` guard
in place a forcing block assembled from more than one model call is the block a
single call would have written, so EFOR-6's replay seam tests the seam it means
to test and nothing else. Every run, climatology and forcing block produced by
more than one model call BEFORE that guard is a different experiment from the one
it was taken to be, and carries a discontinuity at each call boundary whose size
is the divergence measured above.
