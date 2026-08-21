# Restart conversion across horizontal resolution and precision

*This is an implementation specification for converting ExoPlaSim restart
states. Written 2026-08-21 from the current model source and restart artifacts.
It is a design note; none of the converter described here has been implemented.*

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

### Unit tests

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

### End-to-end tests

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

## Delivery plan

### Phase 1: format and schema, 2-3 days

- consolidate the restart parser/writer;
- generate and review the current record inventory;
- implement source/target inspection and linting;
- define the target-template contract.

### Phase 2: precision conversion, 1-2 days

- typed decoding and encoding;
- precision inference and validation;
- cast/error reporting;
- same-resolution target read test.

### Phase 3: spectral projection, 2-3 days

- semantic triangular unpack/pack;
- prognostic/static field policy;
- projection and negative-control tests.

### Phase 4: grid and coupled-surface state, 5-8 days

- Gaussian conservative weights;
- land/ocean class handling;
- reservoir and sea-ice special cases;
- conservation reporting.

### Phase 5: model fixup and wrapper integration, 3-5 days

- converted-restart namelist path;
- model-owned derived-state rebuild and reset;
- template generation;
- `run_exoplasim.py` provenance and surface validation.

### Phase 6: hardening, 4-6 days

- full target-build matrix;
- FPE/bounds smoke runs;
- settling and invariant tests;
- operator documentation.

A useful atmosphere-only prototype is approximately one engineer-week. A
coupled, provenance-complete converter is approximately three to four
engineer-weeks, excluding timestep and vertical-resolution conversion.

## Decisions to make before implementation

1. Whether the exact target template is mandatory even for precision-only
   conversions. Requiring it is safer; omitting it makes the simple path more
   convenient.
2. Whether a mismatched random-seed shape keeps the template seed or requires a
   user-supplied deterministic seed.
3. The fallback for target land/ocean cells with no same-class source overlap.
4. Conservation tolerances for each precision and field class.
5. Whether to increase the restart record limit now or keep all new metadata in
   sidecars for the first release.
6. Whether the target-build template writer is a model mode or a tightly
   controlled one-step wrapper recipe.

