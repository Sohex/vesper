# The conversion decomposition reproduces: `zcnow` against the divergence advance

Worldbuilding frame: this is a correctness record for the Vesper climate model's
threaded build. Everything below is about a Fortran routine and the numbers it
prints; nothing here is about the simulated planet.

`nenergy = 2` adds `CONVDECOMP`, the control `world-0ov` is diagnosed from, and
two of its nine columns did not reproduce across identical runs. This says what
the source actually does, what was changed, and what the fixed instrument was
shown to do. `world-tqh4`.

## The mechanism, read off the source

`spectrala` runs inside the OpenMP parallel region, and the shared build gives
each thread a slice of the spectral state rather than a private copy.
`plasimmod.f90:assoc_spectral` points every partial at its own slice:

    lo = mypid * NSPP + 1
    hi = lo + NSPP - 1
    sdp => sd(lo:hi,:)

`NESP = NSPP * NPRO`, so the slices tile the divergence array exactly and no
byte of `sd` is owned by two threads. Three consequences follow, and the third
is the defect.

**`sd` is advanced by a write through that pointer.** Step 4.b of `spectrala`
executes `sdp = 2.0 * sdt - adm`, which lands in `sd(lo:hi,:)`. That, and not
`mpsyncsp`, is what carries the divergence from t to t+dt; `mpsyncsp` is an
`!$omp barrier` and nothing else (`mpimod_omp.f90`), so it publishes the write
that has already happened. The earlier reading of this defect had the advance
inside `mpsyncsp`, and that is not what the routine does.

**The control's copy is of the WHOLE array.** At the head of the routine,

    if (nenergy > 1 .or. nconvtime > 0) then
       allocate(zcnow(NESP,NLEV))
       zcnow(:,:) = sd(:,:)
    endif

`zcnow` is a local allocatable in a routine called from inside the parallel
region, so every thread has its own and every thread copies all `NESP` modes.

**Nothing separated the copy from the advance.** Between the copy and step 4.b
there is no barrier: steps 1, 2, 3, 3a and 3b are local arithmetic on `NSPP`
slices and the routine's only other barriers are the one above the copy, which
guards the previous phase, and `mpsyncsp` after 4.b. So a thread held up inside
the copy reads slices that other threads have already advanced, and its `zcnow`
is a mixture of the divergence at t and at t+dt, slice by slice, with the mix
set by where each thread was when the copy passed it.

Terms 3 and 4 of the decomposition -- `Ct` and `Dt` -- are the only two that
read `zcnow`. Terms 1 and 2 read `zcsdt`, which arrives through `mpgallsp` and
its two barriers; term 5 reads `zsd`, gathered the same way; term 6 reads `sd`
itself, after `mpsyncsp`. That is exactly the split the symptom showed: seven
columns bit-identical and two that moved.

**A thread's own slice was never at risk**, and this is why the model term is
not affected. The `nconvtime` block indexes `zcnow(jsg,jlev2)` with
`jsg = jsp + mypid * NSPP`, which is `lo:hi` and nothing else. That slice is
written only by this thread, and only at 4.b, which is after the copy in program
order.

## The change

One `!$omp barrier`, immediately after the copy and inside the same branch:

    if (nenergy > 1 .or. nconvtime > 0) then
       allocate(zcnow(NESP,NLEV))
       zcnow(:,:) = sd(:,:)
    !$omp barrier
    endif

It is conditional, and that is legal here because `nenergy` and `nconvtime` are
broadcast by `mpbci` before any of this runs, so the team enters the branch
together or not at all. A configuration that sets neither pays nothing. No array
was added, so the thread team's 32 MB working-set target is untouched. The cost
where the branch is taken is one barrier among the several the routine already
runs, which a 200-step wall clock on this host cannot resolve, so it is not
quoted as a number.

## What it could reach

`config/planet.yaml` declares `energy_diagnostics: true`, which
`run_exoplasim.py:energy_diagnostics_level` renders as `NENERGY = 1`, and it
does not declare `conversion_time_level`, so `NCONVTIME` is never written. A
configured run therefore does not allocate `zcnow` at all and no configured run
was ever exposed. Reaching this needed `energy_diagnostics: 2` or
`conversion_time_level: true` explicitly, which is the decomposition control and
the `world-0ov` model arm.

## The criterion

Fixed before any of the measurements below were taken.

Three runs of one configuration -- one binary, one bed, one restart, one thread
count -- must produce BYTE-IDENTICAL `CONVDECOMP` lines in `plasim_diag`, every
print, all nine columns. Nothing is decomposed differently between the repeats,
so every floating point operation is the same operation in the same order and
the only thing that can move a printed digit is a race; a tolerance would be a
way of not noticing one. The pre-fix binary must FAIL that criterion, and fail
it on `ct` and `dt` alone, or the arm is measuring nothing. `plasim_status` must
be one hash across every run of both binaries, because the decomposition is a
diagnostic and neither the defect nor the barrier may move what the model
integrates.

## The race is intermittent, so the control is forced

Measured 2026-08-26. The copy takes microseconds and a thread needs tens of
microseconds to get from the barrier to step 4.b, so on a quiet host the race
needs the copying thread to be descheduled. It is not rare enough to ignore and
not frequent enough to assume: three 200-step arms caught it in three runs of
three, and three 600-step arms over the same trajectory caught it in none.

An intermittent control is not a control, so the arms below oversubscribe --
every thread of the team on a fraction of the cores, unbound, with
`OMP_WAIT_POLICY=passive` -- which makes the delay routine instead of lucky.
Oversubscription changes no arithmetic: the same operations in the same order,
with only the timing moved. The two binaries are run ALTERNATELY, so a change in
host conditions between them cannot be read as the fix.

## T21, the rung the defect was found at

Bed: `make_profile_bed.py` from `run_14906cb7b914`, restart `7f4c1b6547824e1c`,
T21 l10, 30 minutes a step, dry adiabatic (`nswr`, `nlwr`, `nevap`, `nshfl`,
`nstress`, `nvdiff`, `nprl`, `nprc`, `ndca`, `nshallow` all zero), `nenergy = 2`,
the energy fixer off, `ndiag = 20`, 200 steps, seven prints. Binaries built by
`build_model.py --res T21 --ranks 8 --no-publish`: pre-fix
`37026652ffef1e59`, which is byte-identical to the registered
`most_plasim_t21_l10_p8.x`, and post-fix `09b0382c84a1be2e`. Eight threads on two
cores, five repeats each, alternated.

| | `Ct` at the last print, per run |
| --- | --- |
| pre-fix | 2.19992, 2.21728, 2.17264, 2.21490, 2.23887 |
| post-fix | 2.19130, 2.19130, 2.19130, 2.19130, 2.19130 |

Pre-fix: five runs, five answers, and `ct` and `dt` are the only columns that
move. Post-fix: byte-identical over all seven prints and all nine columns.
`plasim_status` is one hash, `a88cfedb26d5`, across all ten runs.

The post-fix value 2.19130 is also what the pre-fix binary printed at that print
in the unstressed 600-step runs, where the race did not fire. The barrier does
not produce a new number; it makes the number the code was always meant to
compute the only one obtainable.

On the same bed unstressed, three pre-fix runs gave `Cimp - Ct` of -0.00341,
-0.00123 and -0.00618 W/m2 at the last print and three post-fix runs gave
+0.00102 three times.

## T42, the rung the decomposition will be read at

Bed: `exoplasim/bench/bed_check`, a cold T42 start with `SEED = 88888888` and
the two keys the model no longer declares pruned, 45 minutes a step, dry
adiabatic, `nenergy = 2`, the fixer off, 200 steps, seven prints. Binaries built
by `build_model.py --res T42 --ranks 16 --no-publish`: pre-fix
`726288048303a236`, post-fix `ec3c16e74083bea4`. Sixteen threads on four cores,
three repeats each, alternated.

| | `Cimp - Ct` at the last print, per run |
| --- | --- |
| pre-fix | -0.0818, +0.0710, +1.0365 |
| post-fix | -0.1828, -0.1828, -0.1828 |

Pre-fix: three runs, three answers, `ct` and `dt` again the only columns that
move. Post-fix: byte-identical over every print and every column.
`plasim_status` is one hash, `f331b13fcc54`, across all six runs.

**The size of what the defect was worth at T42 is the whole point.** The
displacement `world-0ov` measures is -0.96 W/m2, and on this bed the pre-fix
scatter in the same quantity spans more than a watt. This bed is a cold start
and its magnitudes are not the spun-up dry arm's, so the number to carry forward
is not 1.04 but the fact that the scatter and the effect are the same size at
this rung. An arm taken through the pre-fix instrument could have returned any
sign for the displacement, including the right one.

## `nconvtime` was not affected, tested rather than argued

Same bed, `nenergy = 1` and `nconvtime = 1`, so `zcnow` exists only for the
model term and no decomposition reads it. Ten minutes a step, 200 steps, three
repeats of each binary. All six runs write one `plasim_status`,
`ca7ed8e1d2736373`. The own-slice indexing above says the term never read a
raced value, and the model term reads the same divergence with the barrier as
without it.

That also explains why this survived: `verify_shared_determinism.sh` compares
`plasim_status`, and at `nconvtime = 0` the defect never reached it. The racy
runs above agree with the clean ones on the restart, bit for bit.

## Where the `nconvtime` term runs at T21

Measured on the way, same bed, post-fix binary, 200 steps, `nenergy = 1`. The
model's own guard put the explicit gravity-wave limit at 18.7 minutes and
refused 30 and 20. Below the guard, `nconvtime = 1` trapped on a floating-point
exception at 15 minutes and ran to 200 steps at 12 and at 10; `nconvtime = 0`
ran to 200 steps at 15.

**The failure at 15 is the term and not the timestep**, and that much the bed
shows: the same steps, the same restart, one namelist key apart.

**What it does not show is where the term becomes usable.** Reaching 200 steps
is not stability, and the two things a 200-step bed can tell apart are narrower
than they look. `plasim.f90:stability_check` writes `Abort_Message` and executes
a bare Fortran `stop`, which exits 0, so a run that aborted and a run that
finished carry the same exit status; only a trap, at 136, is distinguishable
from either. And the growth this term produces is slow enough at 10 and 12
minutes that 200 steps do not reach the range that trips
`-ffpe-trap=overflow`, so what separated 15 from 12 on that bed was the bed's
LENGTH.

`world-bt3b` establishes the boundary by deriving it rather than probing for
it. The term takes the reference conversion's divergence half off `sdt` -- the
centred mean of t-dt and t+dt -- and puts it on the divergence at t, and
`sdt - sd` is the second time difference: `O(dt^2)` for a smooth mode and, for
the leapfrog computational mode which alternates sign every step, exactly
`-2 sd`. So the term feeds that mode, the Robert-Asselin filter is the only
thing damping it, and the boundary is a function of `PNU`; at `PNU = 0` there
is none at any step, which no gravity-wave CFL can express. At T21 with the
model's own `PNU`, hyperdiffusion and sponge the boundary is 8.27 minutes
against a guard that admitted 19.0, and the linearised step grows 1.0313 per
step at 10 minutes and 1.0635 at 12 -- amplifications of 471 and 2.2e5 over
this bed's 200 steps, against 4.5e8 at 15.

`exoplasim/analysis/conversion_time_stability_sweep.json` carries the table and
`exoplasim/scripts/conversion_time_stability.py` regenerates it without
building or running the model. **So "an arm that wants
`conversion_time_level` at T21 runs at 12 minutes or below" was wrong: 12
minutes grows, it grows slowly, and a 200-step bed cannot see it.** The model
now measures its own amplification at startup and refuses, so a caller does not
have to know this.

## The boundary, measured against a criterion fixed before the arms ran

T21 l10 p8, warm from `run_0d41aa82c287`'s restart `58a3d2b977af8c8f`, dry
adiabatic, `nenergy = 1`, the fixer off, `ndiag = 200`, 4000 steps, one arm per
(timestep, `nconvtime`). Termination is read from `Abort_Message` and from the
exit status separately, never from the exit status alone.

The criterion, fixed against the derivation before any arm ran: the three arms
at 15, 12 and 10 minutes with the term on must all terminate abnormally within
4000 steps and the arms at 8 and 6 must not, and the step each dies at must
scale with the predicted growth rate.

    dt    term   outcome                        last diagnostic print
    15    on     SIGFPE                         step  200
    12    on     SIGFPE                         step  400
    10    on     SIGFPE                         step 1000
     8    on     completed 4000 steps           step 4000
     6    on     completed 4000 steps           step 4000
    15    off    completed 4000 steps           step 4000
    10    off    completed 4000 steps           step 4000

**Twelve minutes is not a stable timestep; it is a slow one.** It dies inside
600 steps, three times the length of the bed that reported it as completing.
Ten minutes dies inside 1200. The two arms with the term off run to 4000 steps
at 15 and at 10, so the failure is the term and not the timestep, and the
boundary sits between 8 and 10 minutes -- where the derivation puts it, at 8.27.

**The three deaths are one exponential.** The prints are 200 steps apart, so
each arm's death is bracketed, and the predicted growth rates turn those
brackets into the amplification each arm reached: 4.4e8 to 2e17 at 15 minutes,
5.0e10 to 1.1e16 at 12, and 2.3e13 to 1.1e16 at 10. The three intervals overlap
in 2.3e13 to 1.1e16, so a single amplification -- the range from the initial
perturbation to the overflow the declared `-ffpe-trap` refuses -- accounts for
all three deaths at the predicted rates. That is the second half of the
criterion, and it is what says the rate is right and not only its sign.

No arm wrote an `Abort_Message`. `stability_check` tests the divergence,
vorticity and temperature at ONE gridpoint of the top level, and the mode that
kills these arms is at the truncation, so the trap reaches overflow first. An
arm that had died the other way would have exited 0.

**These arms were run on `most_plasim_t21_l10_p8.x` `55a375dbe2cf470f`, which
predates the guard, and a current binary will not reproduce them.** The guard
now measures its own amplification at startup and stops before integrating a
configuration that grows, so the three arms above the boundary are refused
rather than run -- which is the repair working. `conversion_time_arms.py`
records that as a third outcome, neither a completion nor a death, and it is
what a re-take of any comparison at this rung now has to work within: the
boundary can still be approached from below.

## The bed of the section above cannot be rebuilt

Attempting the longer arms found it. `run_14906cb7b914`'s restart was written
on 2026-08-25 and the model has gained the land-tile surface temperature record
since, so a binary built from the current source stops in its first step with
`Error in get_restart_array: Requested array {dlt_ts} was not found`. That is
the same class `lib/rungs.py` records for the earlier T21 restarts, which
`convert_restart` refuses for the partial-cell tile records.

`run_0d41aa82c287` is the replacement and is a better parent on its own terms:
T21 l10 at dt 45, 108 orbits, `equilibrated_for_worldbuilding`, restart
`58a3d2b977af8c8f`, and its restart carries the current record set. A stability
boundary is a property of the configuration rather than of the state, which is
the argument `filter_timestep_matrix.py` track B and `make_profile_bed.py` both
make, so nothing about the question moves with the parent.

## Rebuild

This is a `vendor/exoplasim` change, so every registered binary is stale by rule
4 and the model compiles one executable per (resolution, layers, ranks)
configuration. `exoplasim/scripts/rebuild_binaries.py`, then `--verify`. The
arms above were built with `--no-publish` and left the registry untouched.
