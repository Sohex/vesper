# Vesper

How geography, climate, water and life are computed for this planet, in what
order, and why the order is not a straight line.

Vesper is a super-Earth orbiting a mid-K dwarf: larger than Earth, higher
gravity, a longer day, and more obliquity. Every one of those is declared in
`config/planet.yaml`, which is the only place they are written down. The flux,
the year, the active build and the mean surface temperature move every
iteration; they are in `world_state.json`, which is generated. See
[chapter 5](pipeline/state.md).

**The pipeline chapters are CANONICAL for the pipeline's REASONING**: what each
stage is for, why the order is what it is, why the loops close, and the traps.
The GRAPH itself -- every step, what it writes, what must precede it -- is
`config/pipeline.yaml`, and the two do not overlap. These chapters reference
step ids and do not restate them, because one graph written in two places is
the duplication both exist to prevent.

The rule that follows, and it is a rule rather than an aspiration: **an
artifact that no step in `config/pipeline.yaml` generates does not exist.** If
a script writes something that is in no step, either the step is missing or the
product is one nobody should be reading. Both are defects, and the second is
worse, because an unregistered artifact has no derivable consumers and so no
way to know what it invalidates when it changes. `smoke_test.py` fails when a
generator is absent from the graph. `CLAUDE.md` rule 7 depends on all of this:
"what is now worthless" is only answerable from a graph that is complete.

## How the documentation is organized

Detail lives with the thing it describes. Each component directory has a
`README.md` saying what it does and how to run it. Dated findings live in
`notes/` with their evidence -- `notes/audits/` for audits, and the remaining
files there for measurements and pending changes that carry a date. `TASKS.md`
says what to do about a finding; `world_state.json` holds every current value.

This book holds what crosses component boundaries and does not change with the
current iteration:

- **The pipeline**: the reasoning behind the graph in `config/pipeline.yaml`,
  chapter numbers 0 through 6 preserved from the document this book was split
  from, so citations of the form "section 4" keep resolving.
- **Reference**: builds and identity, the vendored upstreams, the environment,
  the vocabulary, and the standing topic notes.
- **Working practice**: the project's conventions, its working agreements, and
  how it goes wrong ([failure modes](practice/failure-modes.md)).

`CLAUDE.md` at the project root is the map, not the manual: identity, pointers,
and the invariants that bite. Its rule numbers are stable and cited from code.
