# Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software
versions and the terrain hash go into each run manifest and each analysis
report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any
run finished, and were then applied against a result that cleared one of them
by 0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. A bracket's DIRECTION is a claim like any other and needs
its own argument: saying which end is the bound, or that a fuller treatment
would move a correction one way, is a result to be shown rather than asserted
while stating the bracket. `analysis/vegetation_albedo.py` asserted that a
canopy treatment gives a smaller ratio than a leaf one and drew its bracket from
that; the assertion holds for one of CLM's three leaf classes and fails for the
other two, because canopy albedo depends on leaf transmittance as well as
reflectance and the ratio between them is not constant across classes. This is how albedo, evaporation and the carve verdict are
all handled.

Claims are checked against the artifact rather than the documentation.
Almost every class in [failure modes](failure-modes.md) was found that way.

## No tuned values

**A constant whose only justification is that it was fitted, calibrated or
adjusted until a comparison came out is not an acceptable input to this world.**
Replace it with a sourced value, derive it, declare it with a bracket that gets
swept, or record it as irreducible with the argument for why no sourced form
exists. Those four are the whole of the disposition space, and "leave it as it
is" is not in it.

The reason is not tidiness. A tuned value is the one kind of number whose
derivation cannot be inspected, cannot be checked and cannot be carried to a
different planet or star, because there is no derivation: the number IS the
residual of a fit, and the fit was to somewhere else. Every other unsatisfying
number here has a route back to something. An Earth measurement has provenance
for the wrong planet, which is a correction to make. A published parameterisation
has provenance in a paper whose domain can be read and compared against this
world's. A stability bound has an argument. A threshold fixed in advance is a
criterion, and it is allowed to be a judgement precisely because it was fixed
before the result. A tuning has none of that, so nothing downstream of it can
ever be checked either.

**Know which of the four you are looking at, because the sweep is worthless if
it reports every number in the tree.** `notes/audits/tuned-values.md` is the
enumeration and carries the evidence per value;
`notes/audits/opaque-constants.md` and the implicit-Earth audits are the
neighbouring classes and are deliberately separate from it.

**A tuning that names itself is still a tuning.** The model source and the
component configs are now largely honest about which of their constants were
fitted; several say so in their own declarations, and one says the value was
"picked out thin air". Labelling is what makes the sweep possible and it is not
the resolution. The register in `biosphere/config/ntransform.yaml` is the form to
copy: it names each unsourced constant, states what the paper does and does not
supply, and its gate refuses on them under `--strict`, so the label has teeth.

**A tuning that no run reaches is a trap, not a non-problem.** Record it with the
switch that arms it. `oceanmod`'s Earth radius sat dormant under a scheduled A/B
that would have come back wrong by a factor of 1.44 with nothing in the output to
say so; the ocean-diffusivity arms that had already run were relabelled after the
fact. Dormancy changes the priority and not the disposition.

**The cost of removing one is not the output it invalidates.** The canonical
climatology lineage does not exist, so nothing is commissioned, and a change that
"would require re-commissioning" costs nothing today: it is simply what the next
cycle implements. Counting the invalidated runs as a price invents a sunk cost
and biases every one of these decisions toward keeping the tuned value, which is
the outcome the rule exists to prevent. This is CLAUDE.md rule 7 applied to a
place where it is easy to forget.

Two worked examples, one of each verdict, so the line is visible.
`aeolian/config/dust.yaml` chose Kok (2014) over Marticorena-Bergametti because
K14 derives the emitted size distribution from fragmentation physics rather than
fitting it, "so it carries one fewer unconstrained knob. On a synthetic planet
with no aerosol optical depth observations to tune against, a globally tuned
constant is worth nothing." That is the refusal, taken at the point where a
scheme is chosen, which is the cheapest place to take it. Against it,
`canexch.h`'s `ALPHAA_NLIM` is an ecosystem-level scalar chosen to make simulated
carbon pools agree with published estimates of Earth's, and there is no
observation of this world that could replace it: that one is irreducible, it is
recorded as irreducible, and what is done about it is a bounded sensitivity
rather than a fix.

Fitting is not the offence; opacity is. `build_vesper_header.py` fits a solstice
offset by grid search, and it is not a tuned value: it fits against the
climatology's own declination, reports its rms and maximum residuals into the
header and the provenance JSON, and re-derives per planet. The test is whether a
reader can see what was fitted, to what, and how well.

## Documents and numbers

- **Do not write current values into prose.** `world_state.json` is generated
  from the artifacts and is the only place they belong. A number earns a place in
  a document only if it is a decision, a threshold, an identity, or if the
  magnitude carries an argument that fails without it: a convergence criterion, a
  terrain hash, "the two windows do not overlap", "exactly 180 degrees". Mean
  temperature, basin counts, land composition and the active build are none of
  those, and every one of them was wrong in these documents within a day of being
  written. `notes/` is the exception and the opposite: those are dated records of
  what was measured, so they keep their numbers and gain a "measured on".
- **Rewrite superseded content; do not mark it.** A reader grepping for a number
  lands on the number, not on the warning above it.
- Keep citations and table cells on one source line, even where that breaks
  column alignment. They get copied out.
- Claims about convergence and equilibration are stated with their exact criteria
  and are labelled honestly when they miss. One cold case is called
  "quasi-equilibrated" for missing its threshold by 0.004 W/m2. Preserve that
  standard rather than rounding results into passes.
- **An undocumented component is not complete.** Work here is picked up by
  someone with no memory of it -- assume a brick to the head between any two
  sessions, because a fresh session IS that. The test for done is not "does it
  work", it is "can someone who has never seen it find it and use it without
  reading the diff".

  What this does NOT mean is recording numbers. Anything that can be looked up
  or re-derived should be, and writing it into prose is the failure the first
  convention above exists to prevent. What it means is that the THING must be
  findable: a new module belongs in the `lib/` list, a new step in
  `config/pipeline.yaml`, a new component in the layout and in
  `docs/src/pipeline/components.md`, a
  model change as a commit under `vendor/exoplasim`, a new convention as a one-line directive in `CLAUDE.md`, argued here. A
  component that works
  and is invisible will be reimplemented beside itself, which is how this project
  came to have four copies of a path resolver and three of a grid convention.

  One test resolves both: **record what a thing IS and where it lives, never
  what it currently SAYS.**

- **A sentence earns its place by the future work it can inform.** The test
  for any passage is the decision or action a future reader could take
  differently because of it; a passage that changes nothing is byproduct,
  however true. The recurring genres: narration of how the document came to
  say what it says (git carries that); dates on durable arguments and
  definitions (a date belongs on a measurement that will drift, on a
  pre-registered threshold where the ordering is the point, or where a reader
  must judge staleness -- nowhere else); commentary about the document
  itself; and ceremony ("verified rather than asserted") where the evidence
  or the status label already carries the epistemic weight. What is NOT
  byproduct: the one-line cost history that calibrates a rule, the trap
  warning, the status label, and the argument behind a live decision.
  Failure-modes class 25.

## The issue tracker

Work is tracked in `bd` (beads): `bd ready` for what is unblocked, `bd show
<id>` for one issue, `bd list --label clim` for an area. These conventions
predate it and survived the move off a markdown table.

- **Every issue cites the document that justifies it.** An issue with no
  source is not yet a finding and probably needs one. Findings and issues are
  kept apart: a document under `notes/audits/` says what is true and carries
  its evidence, the issue says what to do about it, so a finding can be read
  without being re-litigated and an issue can be closed without editing the
  argument that produced it.
- **An issue that turns out to be wrong is closed with its reason, not
  deleted.** The reasoning is the point: a refutation records why something
  that looked like work turned out not to be, and several of the closed ones
  reverse an earlier conclusion. Deleting them loses the argument and invites
  the same issue being opened again.
- **Ids are never reused, and beads issues them.** A reused id makes every
  citation of it ambiguous forever. Hand-allocated ids did that twice: taking
  the max over six area tables and not the seventh issued a duplicate HYD-17,
  and two different tasks were both issued CLIM-19, which is why the second
  now sits under `clim-19-b`. Areas are labels, not id space.
- **An issue names the step it touches**, as a `step:<id>` label whose id comes
  from `config/pipeline.yaml`. It says WHERE the work lands and is annotation,
  not a verdict. The carve gate -- does anything outstanding still move what
  the verdict is computed from -- is a judgement made by READING the open
  issues, because the answer is what closing one would CHANGE and that is in
  its own prose. Computing it from these markers was tried and removed: a
  marker names a location, not an effect, so the two are not interchangeable
  at any count. `pipeline.py` reads no tracker, so nothing recorded in beads
  can change what the pipeline planner reports.
- **An issue also names the batch it must land in**, as a `batch:<n>` label.
  Where `step:<id>` says WHERE work lands, `batch:<n>` says WHEN, and the
  ordering is by BLAST RADIUS under rule 7 rather than by area or priority:
  batch 1 is a change that makes the build and everything below it worthless,
  batch 2 the compiled model, batch 3 the instruments that judge a run, then
  loop A, loop B, loop C and the escalation in the order `sequencing.md` runs
  them, with a last batch for what gates nothing. The batch is a JUDGEMENT
  taken by reading the issue, exactly as the carve gate is, and for the same
  reason: what a change makes worthless is a fact about its own prose. So
  nothing computes a VERDICT from it -- it does not decide the carve gate, and
  `pipeline.py` reads no tracker -- and the one gate that reads it,
  `smoke_test.py`'s batch check, asks only that every unclosed row HAS one and
  that it is one of the nine. Checking that an annotation is present is not the
  same as taking a decision from it, and the distinction is the one the step
  marker already draws. Maintain it by hand: it goes stale as rows close, and
  nothing regenerates it. A row blocked by a row in a LATER batch is a loop
  back-edge and is expected; a batch whose rows are mostly blocked forward is
  mis-cut.
- **A row that cannot move without the author carries `needs-decision`**, and
  `decision-blocks-loop-a` where the answer gates the current pass. That is a
  different claim from `needs-permission`, which says the design and fixture
  work may proceed and only the EXECUTION -- a model run, a coupled response
  case -- waits. Conflating them buries the handful of answers the project is
  actually stopped on inside the much larger set it is merely not authorised
  to run.
- **Regenerating a derived artifact is a STEP, not an issue.** It belongs in
  the ordering in `sequencing.md`, beside the run that consumes it. An issue
  that says "rebuild X before the next run" is tracking state, and state is
  what `check_consistency.py` and `world_state.json` are for. Archived HYD-16
  and BIO-1 were both such rows and should never have been.
- **An issue is code, physics, a decision or a measurement that settles a
  mechanism.** If closing it would produce only a number that the next
  iteration regenerates, it is not an issue; the durable half is whatever it
  teaches.
- **Blocking and related are different edges.** `bd dep add <issue>
  <blocker>` says the first cannot proceed until the second lands, and drives
  `bd ready`. `bd dep relate` says two issues inform each other -- a consumer,
  a shared boundary, a precedent -- without either waiting. Recording a
  consumer as a blocker makes `bd ready` lie in the direction that hides work.

## Prose registers

- **Write about the simulation in the simulation's terms.** Some of this
  project's subject matter shares surface vocabulary with sensitive real-world
  registers, and prose that borrows one of those registers gets misread by
  automated review at real cost, over content that is only ever a toy model of
  an invented planet. Three rules. Name the simulated subject explicitly: the
  model's snow albedo, the simulation's rainfall, a gridcell's vegetation --
  a sentence should not parse as being about the real world or about people.
  Keep established technical terms whose context already disambiguates them
  (extinction efficiency, albedo, flood seeding); a euphemism for a standard
  term is one quantity with two names, which is its own failure class. And any
  document that can be read standalone states the worldbuilding frame up
  front, as `biosphere/notes/fire-model-audit.md` does. THE RULE IS SUBJECT TO
  ITSELF: do not enumerate the phrases it exists to avoid, here or anywhere -- an enumeration in a
  file loaded into every session re-supplies, on every turn, exactly what it
  guards against. Commits f5de28f, 47a68b3 and 9957a0b hold the specifics for
  anyone who needs them.
