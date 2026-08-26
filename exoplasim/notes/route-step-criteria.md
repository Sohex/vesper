# What it takes to raise a step on the escalation route

*Worldbuilding frame: this is about the Vesper project's climate model and the
order in which its configurations are run. Nothing here is about the simulated
planet.*

**Declared 2026-08-26, before the cells that judge it were measured.** The grid
being re-taken under WORLD-37TN can move the route's per-rung timestep, and a
threshold chosen after the run it judges is not a threshold. So the adoption
rule is written first and the numbers are read into it afterwards.

## What the route is choosing

`docs/src/pipeline/sequencing.md` section D fixes the shape: resolution and
timestep alternate, one variable at a time, and `convert_restart.py` refuses a
conversion across a change of step. **So the five entries carry three free
variables, not five.**

| route entry | what it is | free variable |
| --- | --- | --- |
| 1. T21 at dt A | commissioning | **A** |
| 2. T21 at dt B | settling, so the conversion happens at constant dt | -- |
| 3. T42 at dt B | commissioning | **B** |
| 4. T42 at dt C | settling | -- |
| 5. T85 at dt C | commissioning | **C** |

B is the step the T21 -> T42 conversion happens at, so it is BOTH T42's
commissioning step and the step T21 has to settle to. C is the same for
T42 -> T85. A step therefore cannot be maximised for one rung in isolation:
raising C asks T42 to endure C as well as T85.

The route today is A = 45, B = 30, C = 22.5.

## The rule

**C-ROUTE-1, adoption.** A step is adopted for a rung's COMMISSIONING block only
if both hold:

1. the rung's cell at that step and kappa 8 is `no_refusal_in_steps` in the
   grid re-taken on the current model source, with its executable sha, its
   declared damping and its staging shas recorded; and
2. `lib/rungs.py:COMMISSIONING_EVIDENCE` carries an `endured` row for that
   exact (rung, step) pair, naming the run and the orbits it reached.

A refusal-clean cell is NECESSARY AND NEVER SUFFICIENT. The reason is measured
rather than cautious: T42 at dt 45 is refusal-clean at every column any probe
has measured and blew up in its forty-seventh orbit, and again in its
sixty-eighth with the damping corrected. A 600-step probe is a tenth of an orbit
at dt 45. A route set from refusal cells buys a blow-up on the second day.

**C-ROUTE-2, the conversion step.** The step a conversion happens at must
additionally be `endured` at the DONOR rung, because the donor settles at that
step before it converts. So B needs endurance at T21 and at T42; C needs
endurance at T42 and at T85.

**C-ROUTE-3, no evidence, no raise.** Where the endurance row for a candidate
step is missing, the step is NOT adopted, whatever the refusal grid says. The
route keeps the step it has and the finding records exactly which run would
settle it: rung, step, orbits, and what it costs. A missing row is not a pass.

**C-ROUTE-4, the saving is a number.** The route's total is the sum over its
five blocks of (orbits bought) x (seconds per orbit at that rung and that step).
Orbits come from `lib/run_lengths.py` -- the commissioning window from the
MEMORY time, the settling block from the RELAXATION time, both bracketed --
and seconds per orbit from the probe's cost half. Both inputs are brackets, so
the total is reported as a bracket and never as a point.

**C-ROUTE-5, what a cost number has to clear.** The cost half differences two
wall times on one bed, which cancels startup but not contention, and this host
carries a floor of about five busy cores that no lock can clear. A per-orbit
cost is quoted only if the two passes agree to better than the saving being
argued for; otherwise the ratio between rungs is reported and the absolutes are
not. `docs/src/practice/failure-modes.md` class 34.

**C-ROUTE-6, which model an endurance row describes.** Every row in
`COMMISSIONING_EVIDENCE` today was taken on source that batch 2 has since
changed -- rainmod, fluxmod, radmod, seamod, landmod, plasim and outmod all
moved, and the binaries were rebuilt from the new source. Read strictly,
C-ROUTE-1 would then adopt nothing and freeze the route, which is not what it is
for. So: **an endurance row on superseded source is evidence, and it is
evidence about a model that no longer exists.** It supports keeping a step and
it does not condemn one. A `blew_up` row on superseded source makes its pair
UNMEASURED on the current model rather than known-bad, and the answer to an
unmeasured pair is to measure it. A row taken on the current source supersedes
any older row for the same pair outright.

This is stated here, before the arm that tests it ran, because it is the rule
that decides whether T42 at dt 45 is a ceiling to respect or a candidate to
measure. It is a candidate.

## What this rules in and out before any cell is read

Rule C-ROUTE-3 means the answer can be "no raise is supportable", and that is a
result rather than a failure. It also means the largest saving available is not
whichever step the refusal grid permits: it is whichever step has an endurance
row, and endurance rows are bought one commissioning-length run at a time.

## What the evidence supports, rung by rung

*Read against the rule above, 2026-08-26, on model source `a041a1e9`.*

### The endurance evidence, and what it is worth

`lib/rungs.py:COMMISSIONING_EVIDENCE` is the whole of it:

| pair | verdict | orbits | run | on disk? |
| --- | --- | ---: | --- | --- |
| T21 dt 45 | endured | 50 | `run_ec32946bec89` | yes |
| T21 dt 30 | endured | 35 | `run_14906cb7b914` | yes |
| T42 dt 45 | blew_up | 46 | `run_900548ae632e` | **no, and no stub** |
| T42 dt 30 | endured | 84 | `run_1d39fef9bfc2` | **no, and no stub** |
| T42 dt 22.5 | -- | -- | -- | -- |
| T85 any step | -- | -- | -- | -- |

**Three of the six rows a route decision needs are missing, and two of the four
that exist cannot be checked.** `run_1d39fef9bfc2` and `run_900548ae632e` are
absent from `exoplasim/runs/`, absent from the fifty-seven stubs under
`archive/runs/`, and absent from the one `INDEX_AT_DELETION.json` in the tree.
So the row backing the step the route converts at today is a table entry with no
artifact behind it. That is WORLD-WW6Z's failure mode and it has happened at
least twice on this one pair.

Every row was also taken on source batch 2 has since changed, so C-ROUTE-6
applies to all of them.

### A, the T21 commissioning step: stays 45

T21 is refusal-clean well above 45 -- the boundary is between dt 120 and dt 150
-- and has an endurance row at 45 and at 30 and at nothing coarser. C-ROUTE-1
part 2 is not met at 60, so **A stays 45**. What it would take is one T21
endurance arm at dt 60; T21 is the cheapest rung on the route and the smallest
of the three savings, so it is the last one worth buying.

### B, the step the T21 -> T42 conversion happens at: 30, unless dt 45 endures at T42

T21 endured 50 orbits at dt 45, so the donor can settle at 45. The target
cannot be shown to: T42's only row at 45 is a blow-up on superseded source,
which C-ROUTE-6 makes UNMEASURED rather than known-bad. **B stays 30 until a
T42 endurance arm at dt 45 on the current source says otherwise.**

### C, the step the T42 -> T85 conversion happens at: raise it to 30

This is the recommendation and it is the largest saving on the route.

**C = 22.5 rests on nothing about either rung.** T85 has no endurance row at any
step. T42 has no endurance row at 22.5 either -- the pair has never been run to
a commissioning length. So the step the route's most expensive block runs at is
supported by no endurance evidence on either side of its conversion.

**C = 30 rests on strictly more.** T42 at dt 30 has the longest endurance row in
the table, 84 orbits, which is exactly what C-ROUTE-2 asks of the donor. T85's
own row is missing either way, so raising C from 22.5 to 30 does not give up any
evidence: it moves from a step neither rung has been shown to endure to one the
donor has.

**C = 45 needs two things this session did not establish**: the T42 endurance arm
at dt 45 has to endure, and T85's refusal cell at dt 45 has to be clean.

### The saving, in the unit that needs no clock

An orbit is a fixed span of model time and the model does the same work every
step, so steps per orbit are exactly `4387.23 h / dt` and the saving from
raising a step is exact. Orbit counts are `lib/run_lengths.py`'s brackets: a
commissioning window from the memory time, a settling block from the relaxation
time, an approach of 70 orbits from cold and 15 after a conversion.

| route | T21 steps | T42 steps | T85 steps |
| --- | --- | --- | --- |
| A 45, B 30, C 22.5 (today) | 967k-1774k | 957k-2155k | 1158k-2616k |
| A 45, B 30, C 30 | 967k-1774k | 935k-2107k | **869k-1962k** |
| A 45, B 30, C 45 | 967k-1774k | 913k-2058k | **579k-1308k** |
| A 45, B 45, C 45 | 945k-1726k | **623k-1404k** | **579k-1308k** |

**C 22.5 -> 30 removes exactly 25 per cent of the T85 steps** and 2.3 per cent
of the T42 ones. **C 22.5 -> 45 removes exactly 50 per cent.** Taking B to 45 as
well removes a further 35 per cent of the T42 steps. T85 is the largest single
block in every row, so C is where the money is, exactly as expected.

These percentages are arithmetic and carry no instrument error. What steps
cannot give is the weight of a T85 step against a T21 step, so they do not sum
to one number: that needs the cost half, on a machine this session did not get.
