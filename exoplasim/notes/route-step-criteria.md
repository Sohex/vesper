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

## What this rules in and out before any cell is read

Rule C-ROUTE-3 means the answer can be "no raise is supportable", and that is a
result rather than a failure. It also means the largest saving available is not
whichever step the refusal grid permits: it is whichever step has an endurance
row, and endurance rows are bought one commissioning-length run at a time.
