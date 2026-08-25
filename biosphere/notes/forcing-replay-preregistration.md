# Pre-registration: the LPJ forcing replay and uncertainty protocol

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf; its climate is ExoPlaSim and its biosphere is LPJ-GUESS. What
follows is about how the simulated weather is presented to the simulated
vegetation.

**Nothing here has been run.** Every threshold below is fixed now, with its
derivation, so that a later result can be scored rather than rationalised.
Re-deriving a threshold from a better argument is allowed; moving one after
seeing what it decides is not. A criterion chosen after the run it judges is not
a criterion.

The units are `biosphere/notes/time-base-unit-contract.md`'s; the fields and
their meanings are `biosphere/notes/ecological-forcing-field-contract.md`'s;
the artifact this protocol replays is EFOR-3's and its acceptance checks are
EFOR-7's. None of those is restated here.

## What this protocol decides

LPJ-GUESS spin-up is far longer than the post-equilibrium climate segment it is
economical to archive, so some forcing has to repeat. This protocol says WHICH
sequence the accepted biosphere result is computed on, and how much of that
result's spread is attributable to that choice rather than to the climate.

It is not a sensitivity study of the biosphere. It is the construction rule for
its input, plus the uncertainty that construction carries.

## Definitions, fixed here

**Accepted block.** A contiguous run of whole orbits that: comes from a single
ExoPlaSim model call; lies after the run's declared equilibration; carries the
ecological stream; and passes EFOR-7's acceptance checks.

The single-call requirement is not fastidiousness. `world-8yyh` measured the
model's state failing to reproduce across a restart, with the divergence
appearing in the FIRST timestep after the boundary and reaching 1.19 K in
near-surface air temperature; `notes/audits/ecological-stream-restart-continuity.md`
carries it. A block assembled across a model call therefore contains a
discontinuity of its own, and a protocol whose central test is a seam test
cannot use a sequence with an undeclared seam already in it. If `world-8yyh`
closes, this requirement is relaxed and this paragraph is what says why it
existed.

**Replay.** Presenting a block's intervals to LPJ-GUESS in their recorded order,
wrapping from the last interval of the block to the first.

**Seam.** The join between a block's last interval and its first on the next
pass.

**Phasing.** Which interval of a block the replay starts from. Rotating the
phase moves the seam without changing the block.

## The central case

Cyclic replay of ONE accepted block, whole and in order, with its full interval
chronology intact, beginning and ending at declared orbital boundaries.

For a stellar-cycle run the block is the whole cycle, replayed in its physical
order. The cycle is not shuffleable, its length is the block length, and the
only admissible phase rotation is by a whole number of cycles.

## The thresholds

Every bar below is an artefact-of-construction bar. The rule generating them,
fixed here: **an artefact of how the forcing was constructed must be at most one
fifth of the smallest effect this protocol exists to resolve.** One fifth
because an artefact at that size cannot change the sign or the ranking of the
effect, and a smaller ratio would set bars below what the model's own patch
sampling can resolve.

| quantity | bar | derivation |
| --- | --- | --- |
| land-mean NPP | 1.0% | one fifth of 4.8%, the concavity effect measured on this model over a three-phase forcing cycle (`biosphere/README.md`). That is the smallest effect the protocol has to resolve, and it is a measurement on this model rather than an assumption |
| tree cover, fraction of land with tree FPC > 0.1 | 1.0 percentage point | the only route by which tree cover re-enters the pipeline is surface albedo. Land-mean albedo runs from 0.276 bare to 0.179 vegetated, so one point of canopy is 0.001 of land-mean albedo, two orders of magnitude inside the albedo bracket's own width. An artefact this size cannot move loop C's verdict |
| between-block spread, land-mean NPP | 15% | half the half-width of the registered scored band, 0.7 to 1.4 times Earth's, in `biosphere/notes/productivity-prediction.md`. A spread wider than this would mean the score is decided by which block was drawn |

`world-kso7` owns replacing the tree-cover bar with one derived from a measured
canopy-to-verdict sensitivity rather than from the albedo endmembers.

## The four tests, and what each can fail

### 1. The seam test

Run the central case twice on the same block, phased half a block apart, so the
seam falls in a different place. Compare the equilibrium result.

**Passes** when land-mean NPP differs by at most the 1.0% bar and tree cover by
at most the 1.0 point bar.

**Fails** when either bar is exceeded. A failure says the vegetation is
responding to the join and not only to the weather, so cyclic replay at that
block length is reporting an artefact.

The seam is tested rather than assumed harmless because the model's own
smoothing does not reach it: `mtemp_min20` averages twenty simulation years and
therefore averages across many passes, but daily interception, snow, phenology
and the growth-efficiency mortality buffer do not, and those are where the
concavity effect lives.

### 2. The block-length ladder

Repeat the central case at block lengths of 1, 2, 4, 8 and 16 orbits, drawn from
the same accepted run, each block a prefix of the next where the archive allows
it. Fixed here: the ladder, and the stopping rule.

**Stop** at the first length where BOTH the seam test passes and the result
moves from the previous rung by at most the 1.0% and 1.0 point bars.

**If 16 orbits does not settle it**, retain a bracket over the lengths rather
than taking the longest as truth. This mirrors loop D's rule in
`config/pipeline.yaml` and for the same reason: the longest rung run is not
evidence that it converged.

### 3. Alternative blocks

Run the central case at the chosen length on at least three DISJOINT accepted
blocks. Their spread is the forcing-construction uncertainty, and it is reported
with the answer rather than folded into it.

**Fails** when the spread in land-mean NPP exceeds 15%. That failure is a
statement about the archive, not about the biosphere: the accepted climate is
too short to determine the vegetation, and no central value may be quoted from
it.

### 4. Season-preserving reordering

Permute the ORBITS within a block, each orbit kept whole so no seasonal cycle is
ever broken, and replay the permuted block. Three permutations, drawn with
`numpy.random.default_rng(20260824)` and taken in order, so that which
permutations are used is fixed here rather than chosen later.

**Passes** when the permutations move the result by no more than the
between-block spread from test 3.

**Fails** when they move it further. The answer then depends on the ORDER of the
orbits rather than on the climate's distribution over them, and a protocol built
on repeating one order is reporting that order.

## The weather generator, and the conditions on it

A Vesper-trained stochastic weather generator is admissible ONLY as a labelled
fallback, and only when the archived chronology is too short for a sensitivity
the protocol requires. Three prohibitions, fixed here and not negotiable by a
later result:

- **Never Earth GWGEN.** Its statistics are Earth observations. LPJ-GUESS's
  BLAZE fire model forces `weathergenerator "GWGEN"`, which is one of the two
  reasons the run harness selects GLOBFIRM instead; that substitution is not an
  opening for the generator by another route.
- **Never Earth bias correction.** ISIMIP-style adjustment maps a model onto an
  observed distribution. This world has no observed distribution, so the
  operation is undefined here rather than merely inadvisable.
- **Never merged with a replay result.** A generator result is reported under
  its own label alongside the replay it supplements, never averaged into it.

Before any generator is fitted, the accepted block is split into a training part
and a held-out part, and the split is declared. The generator must then
reproduce, on the HELD-OUT part alone:

| property | what must be reproduced |
| --- | --- |
| wet/dry spell length | the distribution of consecutive wet and consecutive dry interval runs |
| precipitation amount | the distribution of per-interval depth, including the upper tail |
| temperature extremes | the distribution of the interval maxima and minima |
| incident shortwave extremes | the distribution of the interval maxima |
| serial dependence | the lag-1 autocorrelation of each carried field |
| cross-variable dependence | the mean temperature and the mean incident shortwave conditioned separately on wet and on dry intervals |

**The bar, fixed here.** For each row, the generator's quantiles at 0.05, 0.25,
0.50, 0.75 and 0.95 must lie inside a 90% bootstrap interval computed from the
held-out sample itself, resampling whole orbits rather than intervals so the
resampling does not destroy the serial dependence being tested. Resampling whole
orbits also means the interval widens as the held-out part shortens, which is
the correct behaviour: a short hold-out cannot validate a generator, and this
bar says so by becoming easy to pass and reporting how wide it had to be.

Richardson's first-order wet/dry model and Hempel's independent-bias-correction
result are evidence that these properties matter. They are not a
parameterisation for this world, and citing them is not a licence to import one.

## What counts as the protocol failing

Distinct from the biosphere producing an unexpected number. The protocol has
failed, and must be replaced rather than reported, when any of these holds:

1. **The seam test fails at every rung of the block ladder.** Cyclic replay
   cannot then be made to look like continuous weather at any available length,
   and the answer has to come from a block long enough to need no repetition.
2. **The between-block spread exceeds its bar at the longest rung.** The archive
   does not determine the biosphere, and no central value may be quoted from it
   however tidy the central value looks.
3. **Season-preserving reordering moves the answer by more than the
   between-block spread.** The result is a property of one order rather than of
   the climate.
4. **A generator result is reported after failing held-out validation**, or
   merged into a replay result, or trained on anything but this world's own
   accepted stream.

The first three are outcomes the runs can produce. The fourth is a procedural
failure and is listed with them deliberately: it is the one most likely to
happen by drift rather than by decision.

## What is deliberately not registered here

- **The equilibrium criterion for an LPJ-GUESS run.** BIO-15 owns it, and every
  comparison above is between accepted equilibria by that criterion.
- **Which block lengths the archive will actually hold.** That depends on
  climate output that does not exist, and the ladder above is what will be run
  over whatever it does hold, shortened from the top if necessary and reported
  as shortened.
- **The spatial support.** SPAT-8's convergence protocol chooses it, and this
  protocol runs on whichever support that selects.
- **Any number that would be a result.** No central value, no expected spread,
  no expected block length. The point of registering thresholds before the runs
  is lost if an expectation is registered beside them.
