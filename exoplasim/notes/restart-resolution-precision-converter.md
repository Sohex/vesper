# Restart conversion across horizontal resolution and precision

*Worldbuilding frame: the restart is the saved state of the Vesper climate
model. Nothing here is about the real world.*

*This is the contract for converting ExoPlaSim restart states, and the argument
behind it. The converter is `exoplasim/scripts/convert_restart.py`, built out
of `restart_format.py` (framing), `restart_schema.py` (what every record is)
and `restart_transforms.py` (the two operators); `exoplasim/README.md` says
what each is for. Everything up to the model boundary is built and proven by
`convert_restart.py --self-test`. What is NOT built is the half that needs the
target executable: the template writer and the model-owned post-load fixup.*

## Decision

Replace `tools/pusca.f90` with a schema-aware, model-assisted converter. Retain
its triangular spectral mapping algorithm, but not its file traversal, record
list, or compiled-real representation.

The first supported contract is deliberately narrow:

- horizontal resolution changes where `NLON = 2*NLAT` and PlaSim's triangular
  truncation follows from that grid;
- four-byte to eight-byte and eight-byte to four-byte real conversion;
- identical atmospheric, soil, and ocean vertical levels;
- identical timestep, calendar, enabled physics, and physical configuration;
- a clean restart template written by the exact target build and configuration;
- a new initial condition, not a bitwise continuation of the donor experiment.

The design has no dependency on an MPI build or rank layout. The external
spectral representation should remain a canonical PlaSim layout even if SHTns
uses a different internal representation.

## Why the existing converter is not a base to extend in place

`vendor/exoplasim/exoplasim/tools/pusca.f90` reads the file sequentially and
expects an obsolete record sequence. The current file begins with `nstep`,
`naccuout`, `nlat`, `nlon`, `nlev`, and `nrsp`; `pusca` asks for `nlat`
immediately after `nstep`, so it stops on a current file before converting an
array.

It also:

- names obsolete fields such as `st2`, `st3`, `sr1`, and `sr2`;
- omits `naccuout`, the random seed, humidity under the current names, surface
  state, land, ocean, sea ice, forcing state, and diagnostics;
- emits only a few header and spectral records, not a loadable current restart;
- uses the executable's default `real` kind, so source and target precision
  cannot differ;
- assumes every expected record is the next record rather than indexing by
  name;
- performs no grid remapping, invariant checking, provenance, or target-model
  validation.

The salvageable part is `convert_restart_array`: for each level, unpack the two
real coefficients of every `(m,n)` mode into triangular work arrays, then pack
the target triangle. Shared modes survive, modes above a lower target
truncation disappear, and modes newly present at a higher truncation are zero.

## Current restart contract

The restart is a gfortran sequential unformatted file. Every logical quantity
is two physical records:

1. a 16-byte, space-padded ASCII name;
2. its data payload.

Every physical record has a four-byte leading length, its payload, and the same
four-byte trailing length. The files in this repository are little-endian.
There is no type tag, shape, version, index, or configuration fingerprint.

`restart_ini` builds an in-memory name table and `reseek` scans by name, so the
model itself does not require file order. Its current table is fixed at 200
entries and rejects entry 200, making 199 the practical maximum. A current
restart in this tree already contains 199 named quantities. Adding in-file
metadata therefore requires increasing or dynamically allocating that table
first.

Real payload width is the model's compiled default. Integers remain four bytes.
A data length divisible by eight is not evidence that the payload is a
double-precision real array: integer arrays and four-byte real arrays can have
the same length. Conversion must be driven by a schema.

The geometry is:

```text
NLON = 2 * NLAT
NTRU = (NLON - 1) / 3
NRSP = (NTRU + 1) * (NTRU + 2)
NCSP = NRSP / 2
NUGP = NLON * NLAT
```

The spectral packing is equivalent to `real coefficient(2,NCSP)`, in the order
`m = 0...NTRU`, `n = m...NTRU`, with the two real components adjacent for each
mode. Gridpoint payloads are global Fortran arrays with longitude as the
fastest-varying horizontal index.

## Architecture

```text
source restart + source manifest
                |
                v
       schema-aware reader
                |
       field-policy transforms <--- target restart template
                |
                v
      target restart + JSON report
                |
                v
   target-build read and one-step check
```

The target template is the central design choice. It supplies:

- the exact target record set and order;
- target precision and record lengths;
- enabled optional modules;
- target-resolution static surface and forcing fields;
- correct target compiler random-seed shape;
- clean accumulator reset values, including nonzero sentinels;
- a schema compatibility check against the exact target executable.

The converter starts with the template and overlays only state whose policy
explicitly permits transfer. It must never use an evolved target run as an
unexamined source of fallback state.

### Producing the template

The short-term route is a target build run for one step with output cadence one,
so the ordinary shutdown writer emits every enabled record after the model's
own accumulator reset. Validate that every accumulation counter is at its reset
value.

The cleaner long-term route is a target-build template mode that runs normal
initialization, invokes model-owned reset routines, and writes a complete
restart without advancing the climate. That avoids maintaining cold defaults
or sentinel values in Python.

## Components

### 1. Framing and record library

Create one reusable parser/writer rather than extending the three subtly
different restart readers currently under `exoplasim/scripts`.

It must:

- validate every leading and trailing marker;
- reject truncation, negative lengths, duplicate names, malformed name/data
  pairing, and trailing garbage;
- recognize names only in the name position, not by treating any printable
  16-byte data payload as a name;
- retain raw bytes as well as decoded arrays;
- preserve target-template order;
- write through a temporary file and rename only after validation;
- never overwrite the source or target template;
- support a no-write inspection mode.

Little-endian four-byte markers are sufficient for v1. Other marker sizes or
endianness should fail with a precise message rather than being guessed.

### 2. Machine-readable restart schema

The schema needs one entry for every restart name the supported build can emit,
including optional records. Each entry carries:

- logical type: `int32`, compiler seed, or model real;
- symbolic shape;
- subsystem and feature predicate;
- semantic class;
- conversion action;
- reset or target-template policy;
- validation bounds or conservation quantity where applicable.

Suggested semantic classes are:

```text
header
time-and-phase
prognostic-spectral
static-spectral
prognostic-grid
static-grid
derived-grid
accumulator
opaque-build-state
```

Generate the mechanical record inventory from the Fortran
`put_restart_*`/`mpputgp` call sites, then review and check in the semantic
policy. There are more potential call-site names than the 199 enabled in the
current configuration, so deriving the inventory reduces drift. Semantic
actions still require human review.

An enabled name absent from the schema is a hard error. A schema entry absent
from the exact target template is acceptable only when its feature predicate is
false.

### 3. Precision discovery and conversion

Read the integer headers first. Infer source real width using several invariant
records:

```text
len(sp) == NRSP * real_bytes
len(sz) == NRSP * NLEV * real_bytes
len(dls) == NLAT * NLON * real_bytes
```

All available checks must agree on four or eight bytes. A command-line override
may validate the result but should not suppress a disagreement.

Decode transformations into float64 work arrays and cast exactly once when
writing the target payload:

- four to eight bytes is exact but does not recover precision lost earlier;
- eight to four bytes is lossy and must report maximum absolute and relative
  casting error;
- reject overflow and newly introduced non-finite values;
- leave integer records as four-byte integers;
- copy the random seed only if its target payload length agrees, otherwise
  require an explicit replacement seed or retain the documented template seed.

### 4. Spectral resolution conversion

For every transferable spectral field and level:

1. Decode adjacent coefficient pairs into a semantic map keyed by
   `(component,m,n)`.
2. Allocate the target triangular map as zero.
3. Copy every mode present in both source and target.
4. Pack in the target ordering.

For an increase in resolution, every new high mode is exactly zero. For a
decrease, report the L2 norm and maximum magnitude of discarded coefficients
for every field and level.

Do not prefix-copy a packed array. The offset of a later `m` block changes with
truncation even though the shared `(m,n)` coefficients have the same meaning.

Field policy for the baseline configuration:

| records | action |
| --- | --- |
| `sz`, `sd`, `st`, `sq`, `sp` | project shared prognostic modes |
| `szm`, `sdm`, `stm`, `sqm`, `spm` | project the matching leapfrog-history modes |
| `so` | take target spectral orography |
| `sr` | take target restoration state/configuration |
| spectral output accumulators | retain clean target-template reset values |

When humidity is gridpoint/semi-Lagrangian rather than spectral, `sq` and `sqm`
are absent and the gridpoint `dq` policy applies instead.

The future SHTns implementation should import and export this canonical
triangular convention at the restart boundary. If it instead exposes native
SHTns ordering or normalization in the file, that build must supply an explicit
adapter and cross-build restart compatibility becomes coupled to the library.

### 5. Gaussian-grid remapping

Longitude is periodic and uniform. Construct latitude cell areas from the
cumulative Gaussian quadrature weights, then build reusable source-to-target
overlap weights. A separable sparse latitude/longitude representation is
sufficient; a heavyweight external regridding system is not required.

Field policies differ by physical meaning:

| class | examples | policy |
| --- | --- | --- |
| intensive continuous | surface, soil, and ocean temperature | monotone or conservative interpolation within the applicable surface class |
| reservoir per area | soil water, snow water | conservative remap and report global inventory |
| bounded fraction | sea-ice concentration | conservative remap followed by bounded correction |
| fraction plus thickness | sea ice | remap concentration and ice volume, then recover thickness |
| categorical/static | land, glacier, ocean masks | take target boundary initialization |
| static continuous | terrain, albedo climatology, roughness, field capacity | take target boundary initialization |
| geometric | grid-cell area | recompute on target grid |
| diagnostic/derived | relative humidity, active albedo, saturation caches | recompute in the target model |

Land and ocean must be remapped separately under the target mask. Newly created
land or ocean cells need a declared fallback, preferably the target template or
the nearest valid source cell of the same class. The report must quantify the
inventory introduced or removed by mask changes and by any bounds correction.

Monthly/climatological fields with 14 slices are remapped slice by slice when
they are prognostic. Static climatologies come from the target template.

### 6. Model-owned post-load fixup

The restart contains fields that are functions of other restart state. A target
template's cold relative humidity or albedo is not consistent with remapped
soil water, snow, and temperature. Reproducing those formulas in Python would
create a second implementation of the physics.

Add a target-model converted-restart path, controlled explicitly through the
namelist or wrapper, that runs after all component restart reads and:

- recomputes grid areas and derived surface fields;
- enforces land/ocean/glacier consistency;
- recomputes humidity, albedo, roughness, and related caches through existing
  model routines;
- performs documented bounds corrections;
- resets all model and component accumulation windows through their own reset
  routines.

The converter report should list which records are expected to change during
this fixup. A target read/write smoke test can then distinguish intended fixup
changes from unexplained ones.

## Record action summary

| record class | source | target template | output action |
| --- | --- | --- | --- |
| geometry headers | validate | authoritative | target value |
| vertical dimensions | validate | validate | require equality |
| absolute step | authoritative if timestep equal | ignored | copy source |
| forcing phase state | validate/cast | configuration check | source or recomputed target equivalent |
| prognostic spectral | authoritative | shape/schema | project modes |
| dynamic grid reservoirs | authoritative | fallback | remap |
| static boundary state | provenance only | authoritative | target value |
| derived state | inputs only | placeholder | target-model recompute |
| accumulators/counters | ignored | authoritative clean reset | target reset value |
| optional unknown state | unknown | unknown | fail |

`nstep` is important to more than output naming: calendar and stellar-cycle
phase use it. It must not be reset merely because the converted file starts a
new experiment lineage.

## Timestep changes are a separate feature

V1 must require identical source and target timesteps.

The restart stores `nstep`, not elapsed seconds, and the calendar and stellar
cycle derive phase from that count. It also stores the current and previous
leapfrog states. Copying both into a target with a different step duration
changes the represented derivative and can create a computational-mode shock.

Supporting timestep changes later requires all of:

- reconstructing elapsed physical time and a target step count;
- rebasing calendar and stellar-cycle phase exactly;
- rebuilding leapfrog history for the target step, or adding a model-owned
  startup integration path;
- validating output cadence in physical time;
- recording that the state is temporally reinitialized rather than continued.

Do not approximate this in the first converter by merely scaling `nstep`.

## Vertical resolution changes are a separate feature

V1 must also require equal `NLEV`, `NLSOIL`, and `NLEV_OCE`.

Atmospheric vertical conversion needs interpolation in pressure or hybrid
coordinates, conservation of dry mass and tracers, hydrostatic consistency,
and a treatment of below-ground target levels. Soil and ocean layers need
heat- and water-content-aware remapping. This is independent of horizontal
spectral projection and should have its own specification and tests.

## Interface

Proposed explicit command:

```bash
python exoplasim/scripts/convert_restart.py SOURCE OUTPUT \
  --source-manifest SOURCE_RUN/run_manifest.json \
  --target-template TARGET_TEMPLATE \
  --target-config config/target.yaml \
  --report OUTPUT.conversion.json
```

Useful inspection commands:

```bash
python exoplasim/scripts/convert_restart.py SOURCE --inspect
python exoplasim/scripts/convert_restart.py SOURCE --check-template TARGET
python exoplasim/scripts/convert_restart.py SOURCE OUTPUT ... --dry-run
```

Resolution and target precision should normally be inferred from the target
template and checked against the target config. Explicit values are assertions,
not competing sources of truth.

The converter must refuse:

- a missing source manifest for a resolution change unless every otherwise
  absent compatibility value is supplied explicitly;
- different timesteps or vertical dimensions;
- incompatible record sets or physics switches;
- changed physical configuration other than the approved static target-grid
  representation;
- an unclean target template;
- unknown records;
- overwrite of an existing output without an explicit force option.

## Conversion report and run provenance

Write provenance beside the output rather than adding restart records until the
199-record limit is fixed. The JSON report should contain:

- converter/schema version and source-code hash;
- input, target-template, and output SHA-256;
- source and target build/executable identities;
- source and target configuration hashes;
- source and target `NLAT`, `NLON`, `NTRU`, `NRSP`, vertical dimensions,
  timestep, and real width;
- every record's type, source/target byte lengths, and action;
- casting error per record;
- discarded spectral norms;
- grid conservation residuals, clipping, mask changes, and fallback cells;
- accumulator/reset treatment;
- target validation command and result.

`run_exoplasim.py` should ingest this report into `initial_state`, validate the
target surface hashes, and record the original donor as provenance. The
converted state begins a new run lineage. It is an initial condition whose
eventual equilibrium is set by the target model, not evidence that the source
and target are the same experiment.

The existing `--restart-from` surface guard must recognize a valid conversion
report. Today it correctly refuses a donor restart whose embedded roughness,
capacity, or albedo would supersede newly staged surface files. A converted
restart based on the exact target template resolves that problem only if the
report's target surface hashes match the current run.

## Verification

### Built and passing, none of them needing the model

`convert_restart.py --self-test`, six sections against the restart on disk.
Every positive claim carries a control that fails: a prefix copy of the packed
spectral array, a longitude axis without its wrap, a record of the wrong width
out-voted rather than refused.

- malformed, truncated, wrong-endian, and marker-mismatch files are rejected;
- parser/writer round-trip is byte-identical;
- integer, seed, scalar-real, spectral, and grid shapes are distinguished by
  schema rather than length heuristics;
- four/eight-byte precision inference agrees across invariant records;
- four-to-eight conversion is exact;
- eight-to-four conversion reports the exact cast error and rejects overflow;
- spectral packing matches the model's `(component,m,n)` order;
- up-conversion preserves shared modes and zeroes all new modes;
- down-conversion preserves shared modes and reports discarded norms;
- periodic-longitude and Gaussian-area remapping preserve constants and global
  integrals;
- sea-ice volume survives the concentration/thickness special case;
- unknown fields and configuration mismatches fail before output is created.

Every positive transformation test needs a negative control: a swapped
coefficient pair, wrong triangular offset, nonperiodic longitude edge, or
incorrect field policy that the test demonstrably rejects.

### End to end

The first five need no model and are in the self-test. The rest need the target
executable and a template it wrote, and are tracked separately.

1. Same resolution and precision, preservation mode: byte-identical output.
2. Same resolution, eight to four to eight bytes: only declared cast error.
3. T21 to T42 to T21: every shared spectral coefficient restored exactly at a
   fixed precision.
4. T42 to T21: discarded high-mode norms agree with an independent calculation.
5. Target static fields: byte-identical to the target template after conversion.
6. Target executable: converted restart loads and completes one step under
   bounds checks and floating-point traps.
7. Target read/write: every record name, type, and payload length matches the
   target schema after the model rewrites it.
8. Short integration: all fields remain finite and physical bounds hold.
9. Inventory: dry mass, atmospheric water, soil water, snow, ice volume, and
   ocean heat residuals stay inside declared tolerances or are fully accounted
   for by target-mask changes.
10. Settling run: the converted target approaches the same climate as a cold
    target control after an explicitly excluded adjustment interval.

The last test is statistical and physical, not a restart checksum. Resolution
conversion deliberately changes the represented state and the model is chaotic.

## What remains, and why it stops where it does

The converter, the template writer and the run integration are built. What is
left is one thing, and it is blocked on the resolution ladder rather than on
the converter.

**The template writer is settled and is not a model change.**
`build_restart_template.py` cuts a template from a short run of the target
build -- one T21 orbit costs 11.5 s of wall clock -- and normalises its
accumulation window through `reset_restart_accumulators.py`. The objection to a
Python-side writer was that it would mean maintaining the model's cold defaults
and sentinel values by hand; that is answered, because the clean values are
derived from `outreset` by `restart_schema.py` and held against it by
`scripts/smoke_test.py` rather than typed out. The one-step recipe was rejected
for a different reason than expected: a run does NOT end on an output boundary,
so its restart carries a partial window whatever it is asked to do. The
converter refuses a template that has one.

**`run_exoplasim.py` reads the report.** A converted restart carries its
conversion report beside it, the report's surface hashes are checked against
the run's staged `.sra` files, and the donor is recorded as provenance rather
than as lineage. A restart with neither a run manifest nor a report beside it
is now refused: that case used to skip the surface guard entirely, because the
guard reads the donor's manifest and there was none.

**The model-owned post-load fixup is REFUSED, on a measurement at a real rung
change.** The derived records -- `dalb`, `dsalb1`, `dsalb2`, `dz0`, `dqsat` --
arrive holding the template's values, and the model recomputes each of them
during its first timestep. Measured at T21 to T42 by running the same converted
state twice, once as the converter produces it and once with those five
replaced by the donor's own values remapped onto the target grid, which is the
perfect-fixup proxy because the donor's derived state is consistent with the
prognostics that were converted.

The template was a one-orbit T42 cold start carrying 188 times the donor's ice,
so its albedo was 0.104 out in the area-weighted mean and 0.63 at worst: 33.6
W/m2 of reflected shortwave instantaneously.

IT LASTS EXACTLY ONE TIMESTEP, and the arithmetic says so rather than a
comparison. The first output record averages 160 steps, so 33.6 W/m2 predicts
0.210 on that record; measured, 0.213. An error surviving into a second
timestep would double it.

A fixup therefore removes a one-timestep transient worth 0.21 W/m2 on one
record and 0.005 W/m2 over the orbit, against an energy fixer this project
carries at 0.42. It is not worth a model change and its rule-4 rebuild. What
does the work instead is `first_record_tainted` on the segment, which a
consumer already refuses.

The arms are `run_953ee807d32f` and `run_2c42c68fe9ca`, one T42 orbit each,
seeded from the same conversion of `run_2b20e3324bb0`'s restart onto a template
cut from `run_8102b89a08ac`. All three ran on executable `eae6b0e89357`, so the
comparison is paired on everything but the five records.

`dqsat` is identically zero in every restart on disk -- written, carried, never
populated -- so its rebuild is trivially satisfied.

**THE TEMPLATE'S OWN STATE IS WHAT MATTERS, and it is not the derived records.**
On the same conversion, 124 target land cells and 132 target ocean cells found
no source of their own class: the cells a moved coastline creates. They take
the template's value by declared fallback, and with a cold template that is
the template's ice -- 103% of the sea-ice volume residual is exactly that. The
derived records recover in a timestep; a prognostic reservoir arriving on a
coastline cell does not. Cut a template from a run in a state near the donor's,
and read `from_template_fallback` in the report before trusting an inventory.

## The decisions the contract rests on

1. **The exact target template is mandatory, precision-only conversions
   included.** A conversion without the target's own record set, precision and
   static fields is a guess at what the target expects, and the convenience of
   the simple path is not worth a second way of being wrong.
2. **A mismatched random-seed shape is refused, not resolved.** The seed's
   length belongs to the compiler, so neither the donor's nor the template's is
   inferable as the right answer: the caller names one, with `--seed` or
   `--keep-template-seed`.
3. **The fallback for a target cell with no same-class source overlap is the
   target template**, and every such cell is counted in the report. It is the
   only value on hand that belongs to the target grid.
4. **The tolerance on a configuration real is one part in a million**, fixed
   before any conversion was run. It is far looser than the cast a value may
   have crossed and far tighter than a configuration change, which moves these
   by percent. Conservation is not given a tolerance at all: the remap conserves
   a global integral to rounding, and where a mask makes that impossible the
   residual is reported rather than tested against a number.
5. **New metadata goes in a sidecar.** `restartmod.f90` fixes its name table at
   200 entries and rejects entry 200, a current restart holds 199, and raising
   the limit is a model change that makes every binary stale under rule 4.
6. **The template writer is not settled**, because it is the one decision that
   cannot be taken without running the model. Its two candidates are under
   "What remains" above.

