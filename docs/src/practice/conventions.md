# Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software
versions and the terrain hash go into each run manifest and each analysis
report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any
run finished, and were then applied against a result that cleared one of them
by 0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. This is how albedo, evaporation and the carve verdict are
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

  The two rules pull in opposite directions and that is the point: **record what
  a thing IS and where it lives, never what it currently SAYS.**

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
  front, as `TASKS.md` does. THE RULE IS SUBJECT TO ITSELF: do not enumerate
  the phrases it exists to avoid, here or anywhere -- an enumeration in a
  file loaded into every session re-supplies, on every turn, exactly what it
  guards against. Commits f5de28f, 47a68b3 and 9957a0b hold the specifics for
  anyone who needs them.
