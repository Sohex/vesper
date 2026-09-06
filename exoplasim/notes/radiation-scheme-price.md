# What the radiation scheme costs, and what a band-resolved one would

*Worldbuilding frame: this is a COMPUTE measurement of the Vesper project's
climate model on this desktop, and a comparison against a candidate radiation
code held in `references/`. Nothing here is about the simulated planet. Measured
2026-08-25 against `radmod.f90` as it now stands, at NLEV = 10, the layer count
every registry binary is built with.*

CLIM-84 asks the scheme to be priced before the routines in it are optimised.
CLIM-61 asks a band-resolved shortwave and longwave to be priced against the
broadband scheme. They are one question and this note answers the cost half of
it, plus what has become of the accuracy half since both rows were written.

## The instrument is retired instructions, not wall clock

Every number below is a retired-instruction count from `perf stat -e
instructions:u`, or a static count of call sites in a compiled object. That
choice is not a convenience. An instruction count is a property of the code and
the input; it does not move with what else the machine is doing, and it
reproduces here to about one part in 10^8 across repeats. A wall clock on a
shared desktop does not, and a contended timing number is worse than no number
because it looks like a measurement.

WHAT THAT INSTRUMENT CANNOT SAY. Instructions are not cycles. A vector call is
more instruction-efficient than eight scalar ones and need not be eight times
faster, because the limit may be elsewhere. So this note establishes the SHAPE
and the ORDER of the two schemes' cost and refuses to convert either into a cost
per timestep. That conversion is what still needs a quiet machine.

## The present scheme is O(NLEV^2) in scalar `pow` behind masks

Two facts about `lwr` set the whole cost.

**Its absorptances are evaluated per LEVEL PAIR, not per level.** The
transmissivity loop is `do jlev = 1, NLEV` around `do jlev2 = jlev, NLEV`, which
is NLEV(NLEV+1)/2 = 55 pairs at NLEV = 10. A broadband absorptance is a
nonlinear fit to a whole path, so optical depths cannot be added and every pair
of levels has the fit re-evaluated on its own accumulated amount. The inner body
carries 19 transcendental evaluations: 4 non-integer `**`, 6 `EXP` and 9
`LOG`/`LOG10`, counting the CH4 and N2O band CLIM-42 added. That is 1045 per
cell per longwave call.

`swr` is O(NLEV): 31 transcendental sites inside its level loops, of which 5 are
non-integer `**`, for 310 per cell per shortwave call.

**Every non-integer `**` in the radiation compiles to a SCALAR `pow`.** All 22
source sites sit inside a `where` block, and a masked array assignment is not
vectorised by GCC, so libmvec's eight-wide `pow` is never reached. Compiling
`radmod.f90` with the declared flag line and reading the object confirms it:
`nm -u` lists `pow`, `exp`, `log`, `log10` and the vector `_ZGVeN8v_exp`,
`_ZGVeN8v_log`, `_ZGVeN8v_sin`, `_ZGVeN8v_cos` -- and no `_ZGVeN8vv_pow` at all.
The relocation count is 87 scalar `pow` call sites, 146 scalar `exp`, 98 scalar
`log` and 25 `log10` against 169 vector `exp` and 44 vector `log`.

That is why radiation's transcendentals are 9.74% of T170 in
`where-the-time-goes.md`, `pow` alone 6.69%.

## What each transcendental costs, measured

Retired instructions per evaluation, glibc, `-march=znver4`, loop overhead
subtracted by a control that does the same arithmetic without the call. The
scalar figures use the model's own exponents and argument ranges; the vector
figures are the same expression written so GCC can vectorise it.

| form | instructions per element |
| --- | ---: |
| scalar `pow`, masked, as `radmod` writes it | 116 |
| scalar `pow`, unmasked C loop, model exponents | 111 to 112 |
| scalar `exp`, masked | 45.3 |
| scalar `log` | 53.0 |
| scalar `exp(a*log(x))` for `x**a` | 93 |
| vector `pow`, 8-wide `_ZGVeN8vv_pow` | 15.2 |
| vector `exp`, 8-wide `_ZGVeN8v_exp` | 5.1 |
| vector `exp(a*log(x))` | 12.3 |

Two things fall out of that table before any comparison.

**The exp/log identity CLIM-84 lists as a fallback is worth 17%, not a factor.**
93 against 112 for the same result, and it is less accurate. `merge()` in place
of `where` does not recover the vectorisation either: measured at 107
instructions per element, indistinguishable from the masked form, because both
branches are then evaluated and GCC still declines the vector call.

**Removing the mask is worth 7.6x on `pow` and 8.9x on `exp` IN INSTRUCTIONS,
and 4.25% of the model's wall clock.** WORLD-43RK took the arm and measured it in
place at T21: unmasking `lwr`'s four absorptance selections removes 23.79% of the
whole model's retired instructions and 4.25% of its run, and IPC falls from 1.601
to 1.292. The masked scalar `pow` calls are independent per element with no
carried dependence, and an out-of-order window was absorbing most of their
latency; the vector `pow` that replaces eight of them has neither the latency nor
the throughput to be eight times faster. **So a per-element instruction ratio in
this table is an upper bound on a cost ratio and over-states it by about 5.6
here.** `masked-radiation-and-four-bytes.md` carries that measurement, and the
change remains a separate decision rather than a free win: the vector library's
ULP bound is looser than the scalar one, so it changes results.

## The two schemes side by side, transcendental work only

Per cell per radiation call, at NLEV = 10. The present scheme's `pow` is scalar
by measurement above; its `exp` and `log` are a mixture, so those are bracketed
between the vector and masked-scalar costs.

| | evaluations | instructions |
| --- | ---: | ---: |
| present `lwr`, 55 level pairs x 19 | 1045 | 30,200 to 66,700 |
| present `swr`, 10 levels x 31 | 310 | 7,200 to 18,200 |
| **present total** | **1355** | **37,400 to 84,900** |
| ga7 longwave, 88 to 166 solves x 10 layers x 1 to 2 `exp` | 880 to 3,320 | 4,500 to 16,900 |
| ga7 shortwave, 46 to 97 solves x 10 layers x 1 to 2 `exp` | 460 to 1,940 | 2,300 to 9,900 |
| **ga7 total** | **1,340 to 5,260** | **6,800 to 26,800** |

The solve counts are read out of the spectral files themselves rather than from
their summary headers. `references/socrates/data/spectra/ga7/sp_sw_ga7` block 5
holds 22 (band, gas) entries over 6 bands, 97 k-terms in total and 46 summed over
each band's largest gas; `sp_lw_ga7` holds 51 entries over 9 bands, 166 k-terms
in total and 88 summed over each band's largest. The low end assumes one gas per
band is treated by exponential-sum fitting and the rest by equivalent extinction,
which is how SOCRATES is normally configured; the high end assumes every k-term
gets its own solve. One or two `exp` per layer per solve is the layer
transmission plus, where the two-stream carries scattering, its eigenvalue.

**CLIM-61's premise does not survive this.** The row estimated "about an order of
magnitude of radiative work and ... roughly a twofold whole-model slowdown" from
the band counts alone: 6 shortwave bands with up to 12 k-terms and 9 longwave
bands against 2 and 1. The band counts are right and the conclusion does not
follow, because two things run the other way and together more than cancel it.
The broadband scheme is quadratic in levels where a k-term scheme is linear, and
the broadband scheme's kernel is a scalar `pow` behind a mask at 116
instructions where a k-term scheme's is a vectorisable `exp` at 5.1. On
transcendental work the candidate is CHEAPER at every corner of both brackets.

WHAT THIS COMPARISON EXCLUDES, and it is not small. The candidate's two-stream
and adding algebra, its k-table interpolation in pressure and temperature, its
per-band cloud overlap, its per-band Planck integral, and SOCRATES's own code
overhead are all uncounted; so is the present scheme's non-transcendental
algebra, which is also quadratic in levels. The candidate would also carry
roughly (bands x k-terms x levels) more state per column, which is a working-set
question against the 32 MB per-team target and is not answered here. So this is
evidence that the candidate is not an order of magnitude dearer, and it is not
evidence about which is faster.

## The accuracy half has moved, and it weakens the case for swapping

CLIM-61 was written on the claim that the broadband scheme "returns about 40
percent of a line-by-line answer for CO2 and the same for the trace-gas band".
**That claim was retracted.** `exoplasim/notes/trace-gas-band.md` traced it to
Earth greenhouse figures quoted from general knowledge, never sourced, and
probably not even the same quantity -- an attribution with band overlap shared
out, against a removal with water vapour left in place -- and reports no
identified defect anywhere in the chain.

What replaced it is a comparison that mostly confirms the present scheme.
`corrk-cross-check.md` puts the re-derived shortwave weights against HITRAN2020
correlated-k at the operating paths: `h2osww` agrees to 1.4%, `co2sww` to 0.03%,
the solar-weighted water vapour absorptance to 0.25%, and the Earth check passes
on both sides. Three things are left that a band-resolved scheme would fix
structurally:

- **The per-band attribution inside the CO2 total is wrong** and only survives by
  cancellation: at 2.7 um the derivation charges the band 0.000841 where
  correlated-k gives 0.000016, and at 2.0 um it undercharges by about the same
  amount with the opposite sign.
- **13% of the CO2 shortwave absorption falls outside every band Howard
  measured**, so the reconstruction cannot see it at all.
- **There is no water vapour continuum term**, which is why `h2o_sw_level` is
  carried as a bracket rather than a number.

- And the MAINTENANCE argument the row makes stands untouched: seven derived
  per-star corrections onto Lacis-Hansen and Sasamori, each separately derived
  and citable, exist because a fitted broadband absorptance cannot be re-weighted
  structurally and has to be re-derived per host. A band-resolved scheme performs
  that as configuration.

## The verdict on ordering, which is what CLIM-84 asked for

**Do not optimise the present scheme's transcendentals.** The optimisation would
be work thrown away if the scheme is swapped, which is the exact outcome CLIM-84
exists to prevent, and the two routes inside the 6.69% `pow` share are now
priced. The exp/log identity is worth 17% of a `pow`. Removing the `where` masks
is worth 23.79% of the model's retired instructions and 4.25% of its wall clock
at T21, measured under WORLD-43RK, falling as the rung rises, and it is not
bit-identical.

**And the two-scheme comparison above is not a cost comparison, which the
measurement below then proved the hard way.** It is a comparison of instruction
counts, and WORLD-43RK measured what an instruction count of this kind is worth
in place: the mask removal's 23.79% of instructions bought 4.25% of the clock.
The same correction applies to the candidate's vectorisable `exp`. What that
table establishes is that the candidate is not an order of magnitude dearer IN
TRANSCENDENTAL WORK, and it says in its own last paragraph what it excludes --
the two-stream and adding algebra, the k-table interpolation, the per-band cloud
overlap, the per-band Planck integral and SOCRATES's own overhead. Those excluded
terms turn out to be most of the cost, and the candidate is dearer by a factor of
8.6 in CPU time per column per call. The table is kept because the reasoning in
it is sound and the failure is instructive: a count of one kind of operation is
not a cost, however carefully the count is taken.

**The swap is not compelled by accuracy either.** The 40-percent figure that
motivated it is gone and the re-weighted broadband scheme reproduces
correlated-k to within a few percent on every quantity that has been compared.
What argues for it is structural: three known holes that re-weighting cannot
close, and per-star derivations a band-resolved scheme would perform as
configuration -- two of the ten this project actually runs, by the graph count
below.

## The buy criterion, fixed before the cost per timestep was measured

Declared 2026-09-05, before the candidate had been compiled or run. The
measurement that follows is judged against these and nothing else.

**What the price is compared against is a COMMISSIONING, not a run.** One
commissioning is the bootstrap plus the baseline plus their settling, and
`lib/run_lengths.py` puts its span at 37.8 to 77.4 orbits. At the ladder's
measured wall clock -- `notes/audits/resolution-ladder-wall-clock.md`, 13.45 s an
orbit at T21, 75.74 at T42, 491.96 at T85 -- that is 8 to 17 minutes at T21, 48
to 98 minutes at T42, and 5.2 to 10.6 hours at T85. The rung the criterion is
stated at is T85, because that is where the answer changes.

- **AFFORDABLE at a whole-model slowdown of 1.5x or less.** A T85 commissioning
  stays inside one overnight, so the cadence of loop A is still set by the
  physics rather than by the clock.
- **MARGINAL between 1.5x and 3x.** Buyable, but the commissioning span has to be
  re-derived against the new clock before it is, and the escalation route through
  the rungs is a separate decision.
- **PROHIBITIVE above 3x.** A T85 commissioning no longer fits a day, and the
  iteration cadence of loop A becomes the binding constraint on the world.

**And a second criterion that is not a time.** A thread team's working set on one
die targets 32 MB. The candidate carries roughly (bands x k-terms x levels) more
state per column than the broadband scheme, and at T170 on sixteen the model is
already about eight times over that target
(`docs/src/reference/environment.md`). So the increment is what is judged:
**a candidate that adds more per-die state than the whole present radiation
holds is a regression on the standing constraint whatever it costs in seconds**,
and the disposition then is to block the columns rather than to refuse the
scheme.

**Rule 7 removes the other half of the arithmetic.** No canonical climatology
lineage has been declared, so every build, run and climatology in the tree is
disposable. The price of this swap is what the NEXT cycle runs. There is no
invalidated output to weigh against it and none is priced in below.

**What would make the measurement not a measurement.** The model's own startup at
T21 is about 1.6 s (`where-the-time-goes.md`), and a 2.94 s bed against it
returned a confident number with the wrong SIGN. Both arms below are therefore
sized against their own startup and the size is stated.

## The maintenance half, counted from the graph rather than asserted

Read 2026-09-06 from `config/pipeline.yaml`. CLIM-61's case rests on seven manual
per-star derivations that a band-resolved scheme would perform as configuration,
and the count has always been the argument. What has never been established is
how much of this project's per-star work that actually is.

**A change to the stellar spectrum invalidates 77 steps**, taken as the
transitive closure of `stellar_spectrum` through every `needs` in the graph.
**Ten of those are per-star RADIATIVE re-derivations**: `shortwave_band_weights`,
`cloud_band_weight`, `ice_albedo`, `playa_albedo`, `rock_albedo_bands`,
`snow_albedo_grain`, `snow_albedo_zenith`, `soil_albedo_wetting`,
`surface_albedo` and `vegetation_albedo`.

**A band-resolved scheme retires TWO of the ten and leaves eight.** The seven
corrections are written by `shortwave_band_weights` -- which produces
`shortwave_band_weights.json` and `h2o_sw_level.json`, so `h2o_sw_weight`,
`h2o_sw_level`, `co2_sw_weight`, `ozone_scale` and `ozone_uv_weight` are its --
by `cloud_band_weight` for `cloud_absorption_scale`, and by
`lib/stellar.py:ozone_visible_weight`. The other eight are SURFACE albedos, and
the candidate does not touch them: its albedo argument is one broadband value per
column, so a band-resolved atmosphere does not resolve the surface and PHYS-14
stays a separate decision. That is a correction to this row's own framing, and it
is against the swap rather than for it.

**One of the two is inside the loop, and that is the part worth more than the
count.** `shortwave_band_weights` declares `needs: [stellar_spectrum,
baseline_climatology]`, so it is not a step a spectrum change simply re-runs: it
reads a converged climatology, and re-deriving it after a spectrum change needs a
commissioning first. `cloud_band_weight` needs only the spectrum. A band-resolved
scheme re-weights its k-tables against a stellar spectrum with no run at all --
`examples/trappist1/mk_ga_trappist` rewrites block 2 of a shortwave file and
carries every other block over, with no `corr_k` pass -- so what it retires is
one loop-coupled derivation and one standalone one.

**The surface of the seven, for scale.** Eleven scripts, twenty notes, twenty run
arm files and sixteen lines of `config/planet.yaml` name at least one of them.
That is what a per-star re-derivation has to stay consistent across, and
`check_consistency.py` is what holds it there.

## What the cost per timestep is measured against

`exoplasim/scripts/radiation_cost_per_column.py` -> `exoplasim/analysis/radiation_cost_per_column.json`.

**The model side.** A T21 executable built from the tree's own source at the
invocation, with `--no-publish` so it stays in its build directory and cannot
become the binary a run picks up, and its sha recorded beside the number. The bed
is cut by `make_profile_bed.py` from a settled T21 run, which PRUNES the keys the
current `namelist` statements no longer declare and FORCES what
`config/planet.yaml` says. That is not a convenience: a hand-copied bed carried
`tswr3` in `radmod_namelist`, the key WORLD-F9IG removed when `swr` stopped
reading three tuned coefficients, and the model aborted in `readnl` before its
first step. An aborted model is 0.04 s of wall clock and a radiation share of
ZERO, and a zero share divides into a whole-model slowdown of exactly one, which
is an ordinary-looking answer with nothing behind it. The harness now refuses a
non-zero exit and a long arm that costs no more than a short one.

**The source it is against.** `radmod.f90` is byte-identical to the tree's
canonical copy at the time of measurement, which is what the radiation number
requires. The rest of the model is four commits behind -- a canopy snow store
with interception and unloading, a glacier reference temperature, one specific
heat of ice for both columns, and a comment on an optical conversion -- and none
of them touches `swr`, `lwr` or `radstep`. They add land-surface work, so the
radiation SHARE measured here is a slight over-estimate of the canonical tree's:
the numerator is right and the denominator is a little small.

**The candidate side.** SOCRATES ga7 through the `runes` interface, built by
`exoplasim/scripts/socrates_cost_bench.sh` at `config/planet.yaml`'s own
production flag line less `-ffpe-trap`, so what is compared is the code and not
the compiler. The build is refused unless it reproduces SOCRATES's own
`examples/runes` gfortran reference; it reproduces it to a worst relative
1.14e-08 over 24 of about 400 printed values, which is a last digit at eight
decimal places. ga7 is 6 shortwave bands over 22 (band, gas) entries at 97
k-terms and 9 longwave bands over 51 entries at 166, counted out of the spectral
files themselves.

**The unit.** CPU-seconds per column per radiation call, from `perf stat`'s
`task-clock` with `OMP_WAIT_POLICY=passive`, because the model is threaded and
the candidate's driver is not and both schemes are embarrassingly parallel over
columns. The per-step cost is the DIFFERENCE between an 8,000-step and a
16,000-step arm, so the bed's own startup cancels rather than being assumed
small; `where-the-time-goes.md` records a confident T42 number on a bed shorter
than its startup that reversed sign when the bed was lengthened.

## The cost per timestep, measured 2026-09-06

Load through the run: minimum 2.98, median 6.46, maximum 10.62 over ten samples,
and the measuring job is itself sixteen threads on thirty-two logical cores, so
most of that is its own. T21, sixteen threads, three interleaved rounds, warm
from `run_0d41aa82c287`. The profile attributes 98.1 per cent of its samples to a
named object.

**THE BED IS 46 TIMES ITS OWN STARTUP.** 339,918 ms of CPU on the 8,000-step arm
against a startup of 7,378 ms, and 672,458 ms on the 16,000-step arm; the
difference over the extra 8,000 steps is 41.57 ms of CPU a step. That is the
sizing this measurement had to meet and it meets it emphatically, against the
2.94 s bed at 1.8 times its startup that `where-the-time-goes.md` records
returning a confident number with the wrong sign.

**RADIATION IS 31.99 PER CENT OF T21.** `swr` 10.47, `lwr` 7.51, `radstep`'s own
0.34, and libm 13.67. libmvec is 6.32 per cent and is reported beside rather than
inside, because every non-integer `**` in the radiation sits behind a `where` and
none of them reaches a vector call.

**THE PRICE, per column per radiation call, in CPU-seconds:**

| | present | candidate | ratio | slowdown T21 | slowdown T85 |
| --- | ---: | ---: | ---: | ---: | ---: |
| clear sky, every column lit | 6.49e-06 | 3.39e-05 | 5.22 | 2.35 | 2.16 |
| clear sky, diurnal spread | 6.49e-06 | 3.48e-05 | 5.36 | 2.40 | 2.20 |
| **cloudy, every column lit** | 6.49e-06 | 5.57e-05 | **8.57** | **3.42** | **3.08** |
| **cloudy, diurnal spread** | 6.49e-06 | 5.62e-05 | **8.66** | **3.45** | **3.10** |

**The cloudy rows are the representative ones.** The model runs `nswrcl = 1` and
always carries cloud, and its own 6.49e-06 is measured in situ with that cloud in
it, so cloudy against cloudy is the like-for-like comparison. Cloud costs the
candidate 61 per cent; the zenith treatment costs it 2.6 per cent, which settles
that arm -- neither scheme saves anything worth measuring on a dark column.

**The rung transfer, derived rather than assumed.** Radiation cost per orbit goes
as columns per step and the whole-model cost per orbit is measured, so from
`notes/audits/resolution-ladder-wall-clock.md` radiation grows 32.0 times from
T21 to T85 while the model grows 37.3, and the radiation share at T85 is 0.858
times its T21 value: 27.45 per cent. The independent check is the far end --
`where-the-time-goes.md` measured 24.63 per cent at T170 -- and a share falling
32.0, 27.5, 24.6 across T21, T85 and T170 is the trend the arithmetic predicts.

## The verdict against the criterion fixed before the measurement

**MARGINAL at the clear-sky end and PROHIBITIVE at the cloudy end, and the cloudy
end is the one this model runs.** The bar declared above was 1.5 or less
affordable, 1.5 to 3 marginal, above 3 prohibitive, stated at T85 because that is
where the answer changes. The representative arm lands at 3.08 to 3.10.

**In a commissioning that is 5.2 to 10.6 hours at T85 becoming 16 to 33.** Loop
A's iteration cadence would be set by the clock rather than by the physics, which
is exactly what the prohibitive band was drawn to mean.

**THIS REVERSES THE DIRECTION THE INSTRUCTION COUNT IMPLIED, and the earlier
section already said why it could.** That comparison found the candidate cheaper
at every corner on transcendental work, and listed what it excluded: the
two-stream and adding algebra, the k-table interpolation in pressure and
temperature, the per-band cloud overlap, the per-band Planck integral and
SOCRATES's own code overhead. Those are not a correction to the transcendental
count, they are most of the cost. A count of one kind of operation is not a cost,
and this is the second time on this row that a counted quantity has pointed the
wrong way.

**The working-set criterion PASSES, and it was the one with a disposition
attached.** The candidate carries 42.09 KB of marginal state per column with
10.56 MB fixed for the spectral tables, one copy a process. Eight threads on a
die inside the 32 MB target is 65 columns a thread, and the per-column cost is
flat across a sweep from 25 to 1,600 columns -- 2.59e-05 to 3.62e-05 with no
trend, which is scatter and not slope. So blocking the columns is free, and the
standing constraint is not what refuses this scheme.

**What the price does NOT settle.** It is a cost for `ga7` as shipped.
`sbin/Ccorr_k` can drop absorbers and k-terms this world has no use for, and ga7
carries eight shortwave absorbers including sulphur dioxide and carbonyl sulphide
and twelve longwave including four CFCs, an HCFC and an HFC. A trimmed spectral
file is the obvious lever and it is not priced here; what can be said is that the
gap is a factor of 8.6 and the trim would have to close most of it.

**So the decision is now unambiguous in shape.** The swap is not compelled by
accuracy -- `corrk-cross-check.md` has the re-weighted broadband scheme within a
few per cent of correlated-k on every quantity compared. It is not compelled by
maintenance -- the graph count above says a band-resolved scheme retires two of
ten per-star radiative re-derivations, the other eight being surface albedos it
does not touch. And it now costs a factor of three at the rung that matters. What
remains for it is three structural holes re-weighting cannot close: the per-band
attribution inside the CO2 total, wrong by two orders of magnitude at 2.7 um and
surviving by cancellation; 13 per cent of the CO2 shortwave absorption falling
outside every band Howard measured; and no water vapour continuum term, which is
why `h2o_sw_level` is a bracket. **Those three now have to be worth trebling a
commissioning, and the honest recommendation is that on this evidence they are
not -- but that is a judgement about what this world's numbers must be able to
bear, and it is the author's to take rather than this note's.**
