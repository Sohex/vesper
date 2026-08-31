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

The route today is A = B = C = 45: `lib/rungs.py:ESCALATION_ROUTE` is
three entries at one step, so the settling blocks are gone and the free
variables collapse to one. The rows below were written when it was 45/30/22.5
and the rule they state is unchanged; what has moved is the evidence.

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

*Re-read 2026-08-30, after WORLD-YYX8 bought the T42 arms.*

| pair | verdict | orbits | run | on disk? |
| --- | --- | ---: | --- | --- |
| T21 dt 45 | endured, converged | 108 | `run_0d41aa82c287` | yes, this source |
| T21 dt 30 | endured | 35 | `run_14906cb7b914` | yes |
| **T42 dt 45** | **endured, converged** | **144** | `run_88ed6f9d34ae` | **yes, this source** |
| T42 dt 30 | endured | 84 | `run_1d39fef9bfc2` | no, identity only |
| T42 dt 22.5 | -- | -- | -- | -- |
| T85 any step | -- | -- | -- | -- |

**The two rows the route's first two rungs stand on are now taken on this
source, and both converged.** `run_88ed6f9d34ae` is a T42 cold start at dt 45
that ran 144 orbits and met all six convergence criteria, and its paired
converted arm `run_8d0ae7e2d02c` did the same at 89. Under C-ROUTE-6 that
supersedes `run_900548ae632e`'s blow-up in its forty-seventh orbit outright:
that run was on source batch 2 has replaced, and its record, its provenance and
its reproducer are all gone.

`run_1d39fef9bfc2` and `run_900548ae632e` remain absent from
`exoplasim/runs/`, from every `INDEX_ENTRY.json` stub under `archive/runs/` and
from the one `INDEX_AT_DELETION.json` in the tree. What `archive/runs/` now
holds for each is a `RECONSTRUCTED.json`, which is identity recovered from
surviving analysis artifacts and carries no orbit count and no status, so it
cannot re-read a row in either direction. That is WORLD-WW6Z.

The two T21 dt 30 and T42 dt 30 rows were taken on source batch 2 has since
changed, so C-ROUTE-6 still applies to them.

### A, the T21 commissioning step: stays 45

T21 is refusal-clean well above 45 -- the boundary is between dt 120 and dt 150
-- and has an endurance row at 45 and at 30 and at nothing coarser. The row at
45 is now `run_0d41aa82c287`, 108 orbits on this source and converged on all six
criteria. C-ROUTE-1 part 2 is still not met at 60, so **A stays 45**. What it
would take is one T21 endurance arm at dt 60; T21 is the cheapest rung on the
route and the smallest of the three savings, so it is the last one worth buying.

### B, the step the T21 -> T42 conversion happens at: 45, and it is measured

**B = 45, and the condition this section carried is met.** It read "B stays 30
until a T42 endurance arm at dt 45 on the current source says otherwise". That
arm is `run_88ed6f9d34ae`: T42, dt 45, cold, 144 orbits, all six convergence
criteria met, on this source. Its paired converted arm `run_8d0ae7e2d02c` ran 89
at the same pair and also converged. C-ROUTE-2 asks the donor to endure the
conversion step and T21 endures 45 over 108 converged orbits, so both sides of
the T21 -> T42 conversion are now carried by a run that exists.

### C, the step the T42 -> T85 conversion happens at: 45 on the donor's side

**C = 45 and the donor half is settled; the target half is not.** C-ROUTE-2
splits into two questions and only one of them is about T42. The donor side is
answered: T42 endures 45 for a commissioning span from both initial conditions,
which is what the conversion needs of the rung it leaves.

**T85 still has no endurance row at any step**, so C-ROUTE-3 applies to the
target side exactly as before: T85's refusal cell at 45 is clean, measured on
the grid built for it, and a refusal-clean cell is necessary and never
sufficient. What would settle it is one T85 commissioning-length arm at dt 45.
That is the most expensive single arm on the route and it is also the last one
the route needs, so nothing cheaper substitutes for it.

**What the T42 arms bought, in the unit the next section uses.** The route runs
45 throughout, which against the 22.5 it once carried at T85 is half the steps
at the most expensive rung, and both settling blocks are gone rather than
skipped. Before these arms that rested on a refusal cell at T42 and a blow-up on
source nobody has; it now rests on 144 converged orbits.

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

### The arm that settled B, and what it cost

`notes/audits/epilog-adenergy-use-after-free.md` blocked this arm outright:
`epilog` freed `adenergy` and then wrote it to the restart under the same
`nenergy > 0` guard, so every multi-orbit run this tree prepared died in its
first orbit. That is fixed -- the deallocation block moved below the restart
write and the ordering is structural now -- and the arm ran.

What it took, recorded beside the load because a wall-clock number on this host
is meaningless without it: T42 at dt 45 on eight threads pinned to one die runs
at about 94 s an orbit at loads of 11 to 28, and the two arms plus their T21
donor came to 341 orbits over about four and a half hours of held host lock, a
third of which was spent queued behind other agents. The estimate made before
the purchase was 100 s an orbit and about three hours of lock for the pair; the
extra came from the declared extension to the top of the commissioning bracket
and from the T21 donor, which was not in the estimate because every T21 restart
on disk turned out to predate the partial-cell tile records the current model
writes, so `convert_restart` refused all of them and a donor had to be cut
fresh.

**What is still unbought is the T85 arm.** Its refusal cell is clean at 45 and
its endurance row is empty, which is the same shape T42 was in this morning.
