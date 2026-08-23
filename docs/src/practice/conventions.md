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
