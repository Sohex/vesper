# The BVOC activation contract

Worldbuilding. Vesper is an invented planet and this document is about the
simulation of it: a vegetation model's dormant volatile-emission module, the
atmosphere a climate model prescribes around it, and the declaration file that
governs whether the two are allowed to be connected. Nothing here is an
observation of anything.

This is a CONTRACT, not a finding. What is true about the vendored source is in
[`bvoc-soa-atmospheric-coupling-audit.md`](bvoc-soa-atmospheric-coupling-audit.md),
which this rests on entirely; what to do about it is in the `bd` tracker under
the `bvoc` label. What is here is the third thing: what the project has DECIDED
must hold before the simulated biosphere's volatile organic source may be
switched on, expressed so that a machine refuses rather than a reader
remembers.

The enforcement is `biosphere/scripts/bvoc_gate.py` over
`biosphere/config/bvoc.yaml`. Every rule below is a named refusal in that gate,
and the names in `CAPITALS` are the codes it prints.

## The short answer

The source stays off. Turning it on requires the declarations in sections 2 to
5, none of which exists, and every one of them is somebody's open issue. That
is a complete answer rather than a deferral: the work below is not "check
again later", it is a list with an owner per line.

## 1. What fail-closed means here, and the one direction it is closed in

`ifbvoc 1` reads like one switch and is four models: the simulated plants'
production of speciated volatile carbon, its oxidation by radicals in an
atmosphere whose oxidant state this project has not derived, the formation and
removal of the organic aerosol that oxidation would make, and what that aerosol
would do to the model's radiation and to its clouds. Only the first exists in
the vendored source, and it is an Earth-calibrated diagnostic that is not closed
in the plant carbon budget.

Fail-closed therefore means: **the absent case refuses, it does not default.**
A default is the mechanism by which an Earth measurement becomes this world's
number without anyone deciding it should, and the vendored PFT file is already
full of exactly that -- several of its emission capacity entries say in the
source that they were copied from another PFT.

The gate is closed in one direction only, and the asymmetry is deliberate. A
biosphere run with the volatile source OFF is a correct run and must not be
blocked, so with no request the gate reports what is undeclared and exits 0. It
is a REQUEST to activate that is refused, and refused with every unmet
precondition named at once rather than one per attempt.

Three places carry that:

- `biosphere/config/bvoc.yaml` holds one field per precondition, and the
  sentinel `undeclared` is the only value any of them ships with. There is no
  default. The file's absence is itself a refusal,
  `BVOC-DECLARATION-MISSING`.
- `biosphere/scripts/bvoc_gate.py` turns each remaining sentinel into a
  `BVOC-UNDECLARED-*` refusal carrying the issue that supplies it, and adds the
  refusals in sections 3 to 5 that no sentinel can express.
- `biosphere/scripts/run_lpj_guess.py` asks the gate and writes `ifbvoc` into
  the generated instruction file from the answer. It is written rather than
  inherited from the imported PFT file, because an inherited zero cannot be told
  apart from nobody having decided.

## 2. What the source has to be able to emit (BVOC-1)

Activation is an output and acceptance problem before it is a chemistry
problem, because the vendored fork's output tables cannot express what the
module already computes. The emission module carries nine monoterpene compounds
internally; `commonoutput.cpp` declares output parameters for an endocyclic and
an other class and for nothing else. The gate's evidence block reports that
count against nine on every run, read from the vendored source rather than
asserted here.

The chemistry the compounds would be fed to is branch-dependent on exactly the
identity that collapse destroys: the SOA literature this project has read
distinguishes the highly oxidized molecule yields of alpha- and beta-pinene, and
reports isoprene suppressing monoterpene nucleation in some regimes. A total
terpene carbon flux cannot reconstruct either. Adding the tables is BVOC-11.

The declared preconditions, each refusing by name:

| Declaration | What it has to be | Owner |
| --- | --- | --- |
| `source.compound_resolved_output` | monthly flux per monoterpene compound | BVOC-11 |
| `source.patch_resolved_output` | flux carrying the patch axis, not only the PFT axis | BVOC-11 |
| `source.storage_restart_fixture` | a fixture writing, reloading and comparing the monoterpene storage pool across a restart | BVOC-3 |
| `source.carbon_closure` | emitted and stored volatile carbon removed from a named plant pool, and the synthesis energy charged exactly once | BVOC-3 |
| `source.environmental_response_port` | the rotation-rate, photoperiod, year-length, leaf-energy and hydraulic dependencies ported off Earth's constants | BVOC-4 |
| `source.missing_scope_bracket` | a named lower model plus a bracket for the volatile classes the retained module does not represent | BVOC-5 |

The storage row is the one that looks satisfied and is not. The storage pool IS
in the vendored individual serializer, so the restart record exists; the gate
reports that it is there. What does not exist is anything that exercises it, and
a serialized field nobody has round-tripped is an assumption with a data type.

### Acceptance, for output that does not exist yet

`--check-run` rejects an activated run's volatile output rather than accepting
it by default. Every check has an answer that is wrong rather than merely
different, which is what makes it a test:

- `BVOC-OUTPUT-MISSING` -- a compound's table is absent.
- `BVOC-OUTPUT-INCOMPLETE` -- the simulated cell set does not match the one the
  productivity output covers. A cell present there and absent here is a missing
  cell, not a cell with no emission.
- `BVOC-OUTPUT-IMPOSSIBLE` / `BVOC-OUTPUT-NOT-FINITE` -- a negative or
  non-finite emission.
- `BVOC-OUTPUT-CALENDAR` / `BVOC-CALENDAR-UNKNOWN` -- a monthly table without
  one column per month of the simulated year, or a check run with no generated
  header to say how many months that is. The simulated year is shorter than
  Earth's and its twelve months are correspondingly shorter, so a monthly
  emission column here is not an Earth month and the two must not be differenced
  without the conversion the run manifest records.
- `BVOC-PROVENANCE-*` -- a run manifest with no source build or no hash of the
  climate forcing. A volatile flux that cannot be attributed to a climate cannot
  enter a loop whose whole structure is that climate drives emission and
  emission feeds back on climate.

The units are carried rather than assumed: the module reports volatile fluxes
in milligrams of carbon per square metre and the plant carbon pools are in
kilograms per square metre. Six orders of magnitude between an emission and the
ledger it must be subtracted from is not a rounding difference, and BVOC-3 is
where it gets subtracted.

## 3. Putting the inherited Earth numbers on a declared footing (BVOC-2)

The vendored PFT file's emission capacities are Earth measurements taken under
standard laboratory conditions, converted at initialisation into electron
fractions at a reference CO2, temperature, photon flux and daylength. Those are
measurement conditions, not this world's plant traits, and the conversion also
carries the model's Earth photon conversion factor.

**Keeping them is allowed. Keeping them silently is not.** Vesper's biosphere is
a declared Earth-analogue and transplanted Earth biochemistry is its central
case; the audit's finding is not that the numbers are wrong but that a PFT mean
hides the dominant source uncertainty. The species-composition study this
project has read shows why the mean is especially weak here: emission capacity
is not one of the traits the PFTs were constructed by, so aggregation moves the
apparent capacity on its own, in opposite directions for isoprene and
monoterpenes.

The contract is therefore:

1. The capacities, emitter fractions, per-compound speciation, seasonality and
   storage strategies live in PCAR-5's covarying, provenance-stamped plant trait
   registry, with the current PFT means retained as ONE named transplanted-Earth
   baseline rather than as the answer.
2. Every one of them carries a bracket, and the brackets have a FLOOR. The
   floors are declared in `biosphere/config/bvoc.yaml` as multiplicative
   factors, read from the global emissions framework and the
   species-composition study recorded in `references/INDEX.md`. A declared
   bracket narrower than its floor is refused as `BVOC-BRACKET-TOO-NARROW`.
   The direction matters: a narrow bracket is a claim to know this world's
   plants better than the Earth studies know Earth's, and it arrives looking
   like progress.
3. Capacities are not copied between PFTs to fill a gap, and they are not
   retuned to obtain a desired aerosol forcing. A capacity that improves an
   agreement is not thereby a better capacity; physics is not a knob and neither
   is a trait.

The floor is a floor and nothing is computed from it. It can only ever refuse a
declaration, never supply one.

## 4. The oxidant environment those numbers would be used in (BVOC-6)

A volatile flux is not an aerosol. Between them sits an oxidant field this
project has not derived, and the yield that connects them is not a constant:
the SOA literature this project has read finds even the DIRECTION of several
nitrogen-oxide effects to be pathway-dependent.

**No central SOA yield is declared, and a fixed mass yield is a screening bound
that is labelled as one.** What is declared instead is a bracket, and the
contract says what a bracket has to be made of:

- `oxidants.interface` -- the reduced chemistry's signature. It is driven by
  this star's ultraviolet, the modelled ozone column, water vapour, temperature
  and whichever soil and fire nitrogen-oxide sources the project actually
  adopts, and it returns precursor loss AND the changes to ozone, oxidizing
  capacity, carbon monoxide and methane lifetime. A chemistry that returns only
  precursor loss has hidden the half of its effect that feeds back on climate.
- `oxidants.instrument` -- what derived the bracket ends, over what domain. A
  parametric lifetime lifted from another model supplies a SHAPE and not a
  bracket. The partitioned-lifetime methane box recorded in
  `notes/external-model-survey.md` is the clearest case: separate sink
  timescales with oxidant sensitivities to the precursors, a structure that
  transfers, and coefficients fitted to Earth emissions and Earth stratospheric
  chemistry that do not. The tropospheric gas-phase box model recorded as
  acquired in `notes/external-source-inventory.md` is an instrument rather than
  another fit, because it carries the photolysis explicitly and this star's
  ultraviolet is exactly what a parametric lifetime cannot follow. BVOC-6's
  record carries the third case and its refusal: a family of emulators fitted
  as polynomials in ABSOLUTE mole inventories is anchored to Earth's
  atmospheric mass and is refused whatever its fit quality.
- `oxidants.bracket` -- both ends, and they must differ. A bracket whose ends
  agree is a central value with two names.

### The inconsistency that must not be combined silently

`config/planet.yaml` prescribes the simulated atmosphere's methane and ozone
from a photochemical calculation that holds modern Earth's biogenic surface
fluxes FIXED. Driving a volatile source from this world's own simulated
vegetation while those stay fixed is not one atmospheric state: the volatile
flux consumes the radicals that set methane's lifetime, and it either forms or
destroys tropospheric ozone depending on the nitrogen oxides.

The gate refuses that combination as `BVOC-TRACE-GAS-STATE-NOT-CLOSED`, naming
the prescribed keys it found. The refusal clears only when
`oxidants.trace_gas_state_closed` is declared true, which is BVOC-6's and
WET-10's to declare and nobody else's -- WET-10 because the same prescribed
methane is the only thing standing between the dormant wetland source and the
same problem from the other side.

This refusal is the one that would otherwise never announce itself. Nothing
crashes when a biosphere-derived flux is combined with a fixed-Earth-flux
oxidant state; it just produces numbers.

## 5. The cloud arm (BVOC-9)

Mass and optical depth do not determine particle number, and the climate model
has no pathway from aerosol to droplet at all -- the gate's evidence block
reports the absence from the radiation source. That makes the cloud effect a
STRUCTURAL question rather than a parameter, and a structural question decided
after the results are in is not decided at all.

It is therefore pre-registered separately, in
[`bvoc-cloud-sensitivity-preregistration.md`](bvoc-cloud-sensitivity-preregistration.md),
which fixes the three arms, the decision rule and the threshold before any arm
has run. The declaration carries only which arm is in force, what it is carried
by, and the verdict once it has been scored. Two refusals sit here:
`BVOC-CLOUD-OPTICAL-PROXY`, for an arm carried on the prescribed aerosol
optical-depth field, and `BVOC-CLOUD-VERDICT-MISSING`, for building the third
arm before the verdict that authorises it exists.

## 6. What would have to exist to turn it on

In order, because each constrains the next:

1. Compound- and patch-resolved emission tables in the vendored fork
   (BVOC-11), so that identity survives the output layer.
2. The volatile carbon and its synthesis energy closed in the plant ledger, and
   the storage pool round-tripped through a restart (BVOC-3).
3. The environmental response ported off Earth's rotation, photoperiod, year
   and leaf-energy constants (BVOC-4), and the unrepresented volatile classes
   bracketed (BVOC-5).
4. The capacities and their brackets in PCAR-5's trait registry (BVOC-2).
5. The oxidant bracket, with both ends and a named instrument, and the trace-gas
   state either closed or the source held off (BVOC-6, WET-10).

Only then is `requested: true` a question the gate can answer with anything
other than a refusal. Nothing on that list is blocked on machine time; all of it
is blocked on declarations that have owners.

## Sources

The audit carries the literature and the file-level evidence; this document does
not restate either. `references/INDEX.md` records which of those sources have
been read.
