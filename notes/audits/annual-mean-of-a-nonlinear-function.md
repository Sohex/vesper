# A nonlinear function evaluated on an annual mean

An audit of one defect class across Vesper's components, run on 2026-08-27 after
finding it twice in hydrography.

## The class

A quantity is computed by reducing a climatology's time bins to an annual mean
FIRST and then applying a function that is not linear. Because the function is
not linear, `f(mean(x))` is not `mean(f(x))`, and the difference is one-signed
whenever the curvature is. The correct order is to evaluate per bin and reduce
afterwards.

**It does not announce itself.** Both arms produce a plausible number on the
right grid with the right units, the code reads as an ordinary annual-mean
read, and nothing downstream can tell. What makes a case findable is the
nonlinearity, so the useful search is for the OPERATION rather than for the
symptom: exponentials, powers other than one, products of two varying fields,
divisions by a varying field, `min`/`max` clamps -- which are strongly nonlinear
at the clamp and easy to miss -- thresholds, interpolation onto a curved table,
and iterative solves.

Measured in the two hydrography cases the audit started from, the term was
1.170x at the basin sinks and 1.569x over land. It is not a rounding.

## What is fixed

`hydrography/scripts/carve_verdict.py` and
`hydrography/scripts/export_carve_list.py` both evaluated the Penman
combination once on annual-mean air. Both now evaluate per bin and reduce
afterwards, and both carry the annual evaluation as the other end of a declared
bracket. `hydrography/notes/carve-verdict-interval.md` has the argument, the
ocean check that says neither end is universally right, and the numbers.
`hydrography/scripts/surface_water.py` was fixed earlier and is recorded in
`hydrography/notes/lake-balance-integration.md`.

**The correct order is already understood in this project and applied
unevenly.** `penman_open_water` integrates the DIURNAL cycle explicitly, for
exactly this reason, and then was handed annual-mean air.
`pedology/scripts/build_soil.py` reads `mint` and `maxt` per bin to build a
frost fraction, with a comment saying a bin can straddle freezing where the bin
mean never does, and computes its weathering intensity on annual means twenty
lines away. `aeolian/scripts/build_dust.py` runs its emission per bin and
reduces afterwards; `build_sea_salt.py` beside it does not.

## What is found and not fixed

None of the following is in the files this audit's author owns. Each is
recorded with the operation and the input, so the next person does not have to
re-find it. Ranked by how likely the effect is to be large.

| where | the nonlinear operation | on |
| --- | --- | --- |
| `aeolian/scripts/build_sea_salt.py` `main()` | `U10 ** 3.5` per mode, a cubic-with-clamp temperature weight, an iterative Charnock solve for `u10`, and a power law in precipitation | annual-mean wind, surface temperature, humidity and precipitation |
| `pedology/scripts/build_soil.py` `weathering_intensity()` | `exp(T) * Q ** beta`, then a two-sided clip | annual-mean temperature and runoff |
| `pedology/scripts/weathering_fluxes.py` `dessert_co2()` | `runoff * 323.44 * exp(0.0642 * T)` | annual-mean temperature and runoff |
| `pedology/scripts/thermostat_efficiency.py` `read_weathering()` | the same WHAK intensity, used as the WEIGHT in the efficiency number | annual-mean temperature and runoff |
| `minerals/scripts/build_downstream_prospectivity.py` `main()` | the same WHAK intensity | time means taken with an unweighted `.mean(axis=0)` |
| `aeolian/scripts/build_volcanic_sulfate.py` `main()` | hygroscopic growth factor and mass-extinction efficiency interpolated onto a curved table, times a varying burden | annual-mean relative humidity |
| `exoplasim/scripts/dust_forcing.py`, `build_surface_dust.py` | Planck weighting, which goes as `T ** 4` | a temperature collapsed to one annual AND global scalar |

Lower down, and listed because a reader should not have to re-derive that they
were considered: `pedology/scripts/solute_routing.py` and
`phosphorus_budget.py` clamp an annual-mean runoff at zero and do nothing else
nonlinear, and `exoplasim/scripts/cloud_optical_depth_bracket.py` and
`build_surface_albedo.py` apply a curved function to a mean in a DIAGNOSTIC
only, with the shipped number computed per bin.

**The runoff clamp is the one case where the annual order is argued rather than
assumed.** `hydrography/scripts/surface_water.py` and
`config/land_water_ledger.yaml` set out why the annual mean of `P - E` is
exactly the runoff a cell generated over a closed cycle, and why clamping per
bin would count the wet season's supply twice. Anything reading catchment
runoff inherits that argument; the exponential sitting beside it does not.

What would settle each row is the same measurement in every case: evaluate per
bin, reduce, and compare against the annual evaluation on the same inputs. The
ratio is the answer, and it is cheap.

## The declared intervals are honoured

`config/land_water_ledger.yaml` is the only file in the tree declaring an
`interval_floor`. Two terms declare `climatology_bin` -- `lake_precipitation`
and `open_water_evaporation` -- and both producers now evaluate per bin. Every
other term is `instantaneous`, where the producer is the model timestep and
nothing is reduced here, or `annual`. No term declares a floor its producer
does not meet.

## Two defects found on the way, of other classes

**`hydrography/scripts/build_groundwater.py` line 166 unpacks seven values from
a function that returns eight.** `surface_water.climate_fields` gained a staged
albedo as an eighth return; the two callers inside `surface_water.py` were
updated and this cross-module one was not, so it raises on unpack. The comment
in `surface_water.py` recording that the per-bin path "has raised since" for
this exact reason is still there, one caller short.

**Two aeolian scripts defeat the record-count weighting they call for.**
`build_sea_salt.py` and `build_volcanic_sulfate.py` call
`climatology.bin_weights(np.arange(nbin))`, passing an index array where the
bin CENTRES belong. Evenly spaced centres take the degenerate branch in
`lib/climatology.py`, which returns equal weights, so the weighting the module
exists to apply never happens. The real `time` variable is open in both files.
`build_dust.py` passes the real centres and is the comparison.
