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
| T21 dt 45 | endured, converged | 108 | `run_0d41aa82c287` | yes |
| T21 dt 30 | endured | 35 | `run_14906cb7b914` | yes |
| **T42 dt 45** | **endured, converged** | **144** | `run_88ed6f9d34ae` | **yes** |
| T42 dt 30 | endured | 84 | `run_1d39fef9bfc2` | no, identity only |
| T42 dt 22.5 | -- | -- | -- | -- |
| T85 any step | -- | -- | -- | -- |

**WHICH MODEL SOURCE EACH ROW DESCRIBES IS NOT A COLUMN HERE, and that is the
point of C-ROUTE-6 rather than an omission.** A rebuild moves every binary at
once and moves nothing that says so, so a currency column is a statement that
goes false without being edited -- which is what this table did through the
wet-soil merge, while every static gate stayed green.
`lib/rungs.py:evidence_source_currency` derives it per row from the run's own
manifest against `exoplasim/binary_manifest.json`, and `run_exoplasim.py`
prints the answer at launch beside the pair's other caveats. Ask it; do not
read it here.

`run_88ed6f9d34ae` is a T42 cold start at dt 45
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
45 is now `run_0d41aa82c287`, 108 orbits and converged on all six
criteria. C-ROUTE-1 part 2 is still not met at 60, so **A stays 45**. What it
would take is one T21 endurance arm at dt 60; T21 is the cheapest rung on the
route and the smallest of the three savings, so it is the last one worth buying.

### B, the step the T21 -> T42 conversion happens at: 45, and it is measured

**B = 45, and the condition this section carried is met.** It read "B stays 30
until a T42 endurance arm at dt 45 on the current source says otherwise". That
arm is `run_88ed6f9d34ae`: T42, dt 45, cold, 144 orbits, all six convergence
criteria met. Its paired converted arm `run_8d0ae7e2d02c` ran 89
at the same pair and also converged. C-ROUTE-2 asks the donor to endure the
conversion step and T21 endures 45 over 108 converged orbits, so both sides of
the T21 -> T42 conversion are now carried by a run that exists.

### Above T21 the route CANNOT take a cold start, and that is a property of the configuration rather than a preference

*Measured 2026-08-31, pre-flight on the T42 arms world-90y6 bought.*

`config/planet.yaml` sets `model.soil_water_source: pedology`, and it argues
the case: the land-mean capacity is well under the model's uniform 0.5 m
bucket, which changes evaporation, which changes P - E, which is the numerator
of the carve criterion. With it set, `run_exoplasim.py` REFUSES to start
without staged code 229.

**Code 229 exists at T21 and nowhere else, and it is the ONLY field in the
family in that position.** `exoplasim/inputs/t42` carries all thirteen of the
others.

**WHAT PINS IT IS NOT THE EXOPLASIM BUILDER.** That was the first answer and it
was only half of one: `build_surface_soil_water.py` took no `--grid`, read the
grid and the land mask off a climatology and called `require_configured_grid`
on it, so the one rung a run had reached was the only rung the field could be
cut at. It now takes a `--grid` like the rest of the family and takes the
ownership mask from surface code 0172, which `build_boundary_conditions.py`
writes at every exported rung and which a climatology's `lsm` is a copy of index
for index. A climatology passed to it is a cross-check and never a carrier.
Verified by reproducing both staged T21 files bit for bit through the new path.

**WHAT PINS IT IS `land_column_states_<res>.txt`, and that is a real dependency
rather than a carrier.** The states file is per build and per rung.
`pedology/scripts/land_column_properties.py` derives it from `soilmap_<res>.txt`
alone, but `pedology/scripts/build_soil.py` reads a climatology for temperature,
precipitation, runoff, evaporation and elevation -- content, not a grid -- and
refuses one whose shape is not its grid export's. Both now take their rung from
a `--grid` and the climatology declaration is per rung, so what is left is not a
carrier at all: a soil at a rung needs a CLIMATE at that rung, and no
climatology above T21 exists anywhere in this tree.

That leaves every carrier in the chain fixed and its content still at one rung,
and the route's consequence below unchanged until a T42 climate reaches
pedology. WORLD-QGB6 and WORLD-CCX6 carry the measurements.

#### The soil waits for a climate at the rung, and is not remapped onto it

*Decided 2026-08-31, WORLD-512R.*

The alternative was to remap the T21 climatology onto the target grid through
`lib/remap.py` and run `soil` there on the remapped fields, which would put a
land column and a staged code 0229 at a rung before any run had reached it.
**The soil waits.** Three things carry it, and the first is the one that would
still hold if the other two changed.

**A remapped climatology is not a less determined climate, it is another
world's.** `lib/paths.py:best_available_climatology` draws the line the
best-available agreement rests on, and draws it explicitly to separate itself
from the no-fallback rule: the bootstrap and the baseline are the SAME world at
two stages of determination -- same build, same terrain hash, same config,
checked by `require_build` at every call site -- so neither answer comes from
somewhere else. A T21 climatology is a climate that closed on the T21 land mask
and the T21 orography. Remapped, its precipitation and its evaporation balance
on a mask no run integrated, and `build_soil.py` takes runoff as P - E, which is
also the numerator of the carve criterion. That is a plausible number computed
on a different world, which is the failure shape the no-fallback rule names.
The best-available agreement is about the STAGE a step reads, and the rung is
not a stage.

**The saving it was worth is already spent.** The route has to buy a first arm
at each new rung regardless, because a restart template can only be cut from a
run of the target rung on the current staging. That arm's own climatology is
exactly the input `soil` needs, and it is the same `build_climatology.py`
minutes either way. Where the correct input and the approximate one cost the
same, there is nothing to trade.

**Nothing in the schedule would run earlier.** The first arm at a rung above
T21 cannot read a staged 0229 whatever exists, because it is the arm that must
run before the soil at that rung can be built at all; and the route does not
want a cold arm above T21 for its own sake. So a soil built on a remapped
climate would sit unread until the arm that supersedes its input had finished.

What the remap WOULD have to carry if it were ever taken up is stated here so
the question is not reopened cheaply: the product would have to be labelled
remapped through `lib/spatial_support.py`'s identity and semantics contract, at
every artifact it reached, so that no consumer could mistake it for a climate at
that rung. That is a cost on top of the remap and not a mitigation of it.

#### The order that unblocks the route

Once an arm at the target rung has a climatology, the soil chain is minutes and
the order is fixed:

1. `build_climatology.py` on an arm of the target rung whose staged surface
   family is the current one, and DECLARE it -- pointing a component at another
   component's output is deliberate, never defaulted. The declaration is PER
   RUNG: `config/planet.yaml`'s two climatology keys take a rung-to-path
   mapping beside the scalar form, and the scalar answers for `model.resolution`
   alone. `build_climatology.py --output` has to name a directory of its own at
   the new rung, because the product's filename is `<label>_regular_climatology
   .nc` and carries no rung.
2. `soil --grid source/<build>/exoplasim-<rung>` at that rung, then
   `land_column_properties --soil-map soilmap_<rung>.txt --states
   land_column_states_<rung>.txt`, which refuses a pair naming two rungs.
3. `surface_soil_water --grid source/<build>/exoplasim-<rung>`, which stages
   codes 0229 and 2290 there.
4. From then on every arm at that rung, cold or converted, carries the pedology
   capacity: a cold one reads the staged field, a converted one takes `dwmax`
   from a template cut at that rung after step 3.

**Step 2 needs nothing but step 1.** `build_soil.py --grid` is the whole rung
decision: the grid export it integrates the lithology onto, the soil map's
name, BIO-11's rootable fraction and which climatology is resolved all follow
it, and an export of a build other than the configured one is refused. The two
climatology keys are declared per rung, so a rung with nothing declared is
REFUSED naming what to build and declare rather than served the configured
rung's file. `--climatology` remains a cross-check pinned by sha to the file
declared for the grid's rung, which is what keeps BIO-19's mode pinning: a
rebuild of iteration 0 after a baseline exists still reproduces the bootstrap
soil. WORLD-CCX6 measured the carrier at T21 -- the soil map and the land
column states both reproduce bit for bit through the new path -- and
WORLD-QGB6 is the same fix one level down.

`config/pipeline.yaml`'s `surface_soil_water` needs `boundary_conditions` and
`land_column_properties`, which is what it reads. The DETERMINATION argument the
climatology edge stood for -- a staged surface field is an INPUT to the run
whose output would otherwise build it -- is enforced where the climate actually
enters, at `soil`, which names its own stage and refuses a pairing it did not
declare.

**The consequence reaches further than a cold start, because `dwmax` is
`STATIC_GRID, TARGET`.** `restart_schema.py` takes field capacity from the
TARGET TEMPLATE and never from the donor, which is correct -- a donor's
capacity belongs to another soil column -- and it means a CONVERTED arm's
capacity is its own rung's staging too. So the route cannot carry the pedology
field up the ladder by conversion either: every rung above T21 runs the uniform
0.5 m bucket, whichever initial condition it starts from, until 229 can be
staged at that rung. `config/planet.yaml` argues at the key itself that the
uniform bucket is materially wrong here -- the land-mean capacity is well under
it, which changes evaporation, which changes P - E.

So at every rung above T21 the two initial conditions are not
interchangeable, and the asymmetry is one-signed:

| arm | where `dwmax` comes from | possible today at T42 |
| --- | --- | --- |
| cold | the staged `.sra`, read by `landini` only when `nrestart == 0` | **no**, pedology has no land column states at this rung |
| converted | the restart, which FREEZES it | yes, at whatever capacity the template froze |

**This is the route working rather than failing.** The escalation route exists
so that no rung above T21 pays for a cold start, and a converted arm carries
the T21 field up with it. What is new is that the route is now the only
option there, so a cold arm at T42 or T85 is not a fallback the route can take
if a conversion misbehaves.

Two consequences to state before they are needed. A cold arm at a rung above
T21 is buyable only at `soil_water_source: uniform`, which is a DIFFERENT
configuration and has to be declared as one -- it is what the T42 arms of
WORLD-YYX8 ran, before the key was flipped. And the T85 endurance arm this
chapter asks for is a converted arm for this reason as well as for its price,
which C-ROUTE-2 already assumed and which is now also a requirement.

### A restart template is superseded by its own staging, and nothing said so

*Measured 2026-08-31, same pre-flight.*

`exoplasim/inputs/templates/T42_l10_p8.rest` was cut against a staged surface
family that has since gained codes 1743, 1751 and 1761.
`run_exoplasim.py:conversion_surface_reason` catches it and refuses: a
converted state's static surface records are the TEMPLATE's, so a template cut
against superseded staging would hand the run a surface it did not stage.

The refusal is the gate working. What is worth carrying is the ORDER it
imposes, because it is not obvious and it is the same at every rung: **a
template can only be cut from a run of the target rung on the current staging,
so the first arm at a new rung is a cold one and every converted arm at that
rung waits on it.** Where a cold arm is not available -- which is every rung
above T21, per the section above -- the template has to come from the first
run made at that rung on whatever configuration could start, and the
conversion inherits that configuration.

### C, the step the T42 -> T85 conversion happens at: 45 on the donor's side

**C = 45 and the donor half is settled; the target half is not.** C-ROUTE-2
splits into two questions and only one of them is about T42. The donor side is
answered: T42 endures 45 for a commissioning span from both initial conditions,
which is what the conversion needs of the rung it leaves.

**T85 still has no endurance row at any step**, so C-ROUTE-3 applies to the
target side exactly as before, and a refusal-clean cell is necessary and never
sufficient. What would settle it is one T85 commissioning-length arm at dt 45.
That is the most expensive single arm on the route and it is also the last one
the route needs, so nothing cheaper substitutes for it.

**AND THE ARM CANNOT BE BOUGHT TODAY: three of its four inputs do not exist.**
This is not the arm's price and it is not a reason to keep the step; it is the
work in front of the purchase, and it is enumerated so the purchase can be
sized rather than attempted.

1. **The T85 staged surface family is absent.** There is no
   `exoplasim/inputs/t85` and no run directory carrying `N128_surf_*.sra`. The
   T85 cells in `exoplasim/analysis/stability_probe.json` were staged from
   `run_2a50670d8f2b`, which is one of the runs in
   `archive/runs/RECORDLESS.json` and has no record, no payload and no build.
   `source/canonical-10m-carve2/exoplasim-T85` exists, so the family can be
   built; nothing has built it. It also has to carry codes 1743, 1751 and 1761,
   which `config/planet.yaml`'s `surface.soil_albedo_moisture.three_point`
   makes a refusal in `landini` rather than an omission.
2. **The T85 restart template is absent.** `exoplasim/inputs/templates` holds
   `T85_l10_p16_omp.rest.provenance.json` and no `.rest` beside it, and that
   provenance records the template as cut from `precarve-craton-10m` -- the
   pre-carve build, not the active one. `convert_restart.py` needs a target
   template, so the conversion the route makes at this hop has no target.
3. **T85's refusal cell at 45 is on superseded source.** Every T85 probe cell
   names executable `6722f7280e09`, and that is not the sha
   `exoplasim/binary_manifest.json` registers for
   `most_plasim_t85_l10_p8.x` now. C-ROUTE-1 part 1 asks for the cell in the
   grid re-taken on the current model source, so the cell has to be re-taken.
   It is the cheapest of the three and it depends on the first.
4. The T42 state to convert from exists as `run_88ed6f9d34ae`.

**And the arm's own price is a bracket wide enough that it is not yet a
purchase decision.** `lib/run_lengths.py` puts a converted commissioning span
at 53 to 92 orbits. A T85 orbit has never been timed on this host: the ratio to
a T42 orbit is between 4, which is the gridpoint count, and 8, which is what the
Legendre transform scales by, so 6.3 to 12.5 minutes an orbit against T42's
measured 94 s on eight pinned threads. That is **5.6 to 19.2 hours of held host
lock at p8**. A single arm is not a pair, so it may take all sixteen threads,
and the speedup there is unmeasured too.

**What closes the bracket costs minutes, and it is the second thing to buy after
the staging.** C-ROUTE-5's cost half differences two wall times on one bed:
`stability_probe.py` at T85 and dt 45 over 600 and 1200 steps, differenced,
gives seconds per step directly, and 5850 steps is one orbit at this step. Run
the same two passes at p8 and p16 and the thread question closes with it. Until
that number exists the arm is priced across a factor of three and the route
cannot say what raising C is worth in wall clock, which is the whole of
C-ROUTE-4.

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

**What is still unbought is the T85 arm, and it is further off than the T42 arm
was.** T42's endurance row was empty and everything in front of it existed: the
staged family, the restart template, the refusal cell on the source that ran it.
T85 has none of those, and section C above enumerates them. The endurance row is
the last thing bought at this rung, not the first.
