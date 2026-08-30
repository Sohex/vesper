# A nonlinear function evaluated on an annual mean

An audit of one defect class across Vesper's components, run on 2026-08-27 after
finding it twice in hydrography, and re-verified against source on 2026-08-30
with the first measurements taken.

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

`minerals/scripts/build_downstream_prospectivity.py` took its time means with an
unweighted `.mean(axis=0)` and evaluated the WHAK intensity on them. The
weighting is now `lib/climatology.py`'s, and the evaluation interval is declared
as a bracket in `downstream_prospectivity_report.json` rather than chosen. The
measurement is in the next section but one.

`analysis/emissivity_contrast.py`, `maps/build_basemap.py` and
`scripts/world_state.py` each reduced a climatology's time axis with an
unweighted mean. The first two now weight by record count; the third's
reduction was unreachable after the annual mean above it and now refuses
instead of silently collapsing a level axis.

`lib/lapse.py` fits the annual-MEAN profile rather than averaging the per-bin
fits, and a least-squares slope is nonlinear in the profile. It is declared
rather than changed, and the measurement is with the others below.

`scripts/smoke_test.py` gained "no nonlinear function is evaluated on a naive
time mean": a lint over every producer that reduces a netCDF variable over
axis 0 with a bare `.mean()`, and a fixture on the real WHAK intensity that
requires the two arms to agree EXACTLY when the bins are identical and to
differ by exactly `cosh(dT/tau)` when they are not. The identical-bin half is
the control: without it the check could pass on two arms that were not two
evaluations of one function.

**The correct order is already understood in this project and applied
unevenly.** `penman_open_water` integrates the DIURNAL cycle explicitly, for
exactly this reason, and then was handed annual-mean air.
`pedology/scripts/build_soil.py` reads `mint` and `maxt` per bin to build a
frost fraction, with a comment saying a bin can straddle freezing where the bin
mean never does, and computes its weathering intensity on annual means five
lines below it. The proximity is misleading and worth stating, because it is
what decides whether a row is a fix or a bracket: the per-bin `tas`, `pr` and
`evap` are consumed inside the `with nc.Dataset(...)` block and discarded, and
only `bin_min` and `bin_max` survive it. A per-bin weathering intensity there
needs the arrays retained, not the two statements reordered.
`aeolian/scripts/build_dust.py` runs its emission per bin and reduces
afterwards; `build_sea_salt.py` beside it does not.

## What is found and not fixed

None of the following is in the files this audit's author owns. Each is
recorded with the operation and the input, so the next person does not have to
re-find it. Ranked by how likely the effect is to be large.

Every row below was re-verified against source on 2026-08-30, and the last
column is the one that decides the disposition: a site whose binned arrays are
still open can be fixed where it stands, and a site whose dataset has closed
needs the read restructured before the question can even be asked.

| where | the nonlinear operation | on | binned arrays still in scope |
| --- | --- | --- | --- |
| `aeolian/scripts/build_sea_salt.py` `main()` | `U10 ** 3.5` per mode, a cubic-with-clamp temperature weight, an iterative Charnock solve for `u10`, and a power law in precipitation | annual-mean wind, surface temperature, humidity and precipitation | yes, `solve()` closes over them |
| `pedology/scripts/build_soil.py` `weathering_intensity()` | `exp(T) * Q ** beta`, then a two-sided clip | annual-mean temperature and runoff | no, the dataset has closed |
| `pedology/scripts/weathering_fluxes.py` `dessert_co2()` | `runoff * 323.44 * exp(0.0642 * T)`, a product of two varying fields as well as an exponential | annual-mean temperature and runoff | no, the reads are inline temporaries |
| `pedology/scripts/thermostat_efficiency.py` `read_weathering()` | the same WHAK intensity, used as the WEIGHT in the efficiency number | annual-mean temperature and runoff | no |
| `minerals/scripts/build_downstream_prospectivity.py` `main()` | the same WHAK intensity | time means taken with an unweighted `.mean(axis=0)` | yes, and it already reads them for the warmest and coldest bin |
| `aeolian/scripts/build_volcanic_sulfate.py` `main()` | hygroscopic growth factor and mass-extinction efficiency interpolated onto a curved table, times a varying burden | annual-mean relative humidity | yes, `solve()` closes over them |
| `exoplasim/scripts/dust_forcing.py`, `build_surface_dust.py` | Planck weighting, which goes as `T ** 4` | a temperature collapsed to one annual AND global scalar | no |

**The minerals row was the only one taking a naive time mean.** The other six
already weight by record count through `lib/climatology.py`; what they get wrong
is the order, not the weighting. And in the two exoplasim files the GLOBAL
collapse is the larger of the two errors, because the pole-to-equator spread in
`ts` is far wider than the seasonal one, so a per-bin fix on its own would
address the smaller term.

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
ratio is the answer, and it is cheap. Two of them have now been taken and are in
the section below; both came out inside the uncertainty of the law being
evaluated, which is a result about this class rather than about those two rows,
and is the reason the remaining rows are worth measuring before they are worth
changing.

## The declared intervals are honoured

`config/land_water_ledger.yaml` is the only file in the tree declaring an
`interval_floor`. Two terms declare `climatology_bin` -- `lake_precipitation`
and `open_water_evaporation` -- and both producers now evaluate per bin. Every
other term is `instantaneous`, where the producer is the model timestep and
nothing is reduced here, or `annual`. No term declares a floor its producer
does not meet.

## The measurement, where it has been taken

The question each row asks is one number: evaluate per bin, reduce, divide by
the annual evaluation on the same inputs.

**The WHAK intensity, at the minerals consumer.** Measured on
`baseline_regular_climatology.nc` against `canonical-10m-carve2`, on the mesh,
over the cells each rule's climate gates actually select. Only `bauxite` and
`laterite_ni` scale by the intensity; `supergene_cu` and `placer_au` do not
read it.

| | gated cells | bin arm / annual arm | cells below 1 | seasonal T range, median |
| --- | --- | --- | --- | --- |
| bauxite | 14,158 | 1.0087 | 0 | 3.4 K |
| laterite_ni | 19,209 | 1.0080 | 0 | 4.1 K |
| all land, for comparison | 4,328,732 | 1.215 | 0 | 57.0 K |

**The sign came out exactly as convexity requires**, and that is the result
worth keeping: `exp()` is convex, so the bin arm can never fall below the
annual one, and across 33,367 gated cells not one did. Over all land the ratio
runs to 23.9x on a single cell.

**The size is another matter, and it is smaller than the instrument.** The
gates select wet and warm, which on this world is also where the seasonal swing
is smallest: 3 to 4 K against a land median of 57. So the emitted 0-1
prospectivity moves by 0.8 to 0.9% in the area mean and 0.12% at the best cell.
Against that, `pedogenesis.yaml` records Dunne's `S_y.x` as 0.13 log units, a
factor of 1.35 on the law itself, and the choice of temperature e-folding
between 13.7 and 8.91 K is a further 1.4 to 1.5x on the thermostat. A 0.9%
correction sits two orders inside the law's own uncertainty. It is reported as
a bracket rather than applied as a fix, and the bracket is in the report the
consumer already reads.

**Taking the runoff per bin as well is a much larger move -- 1.25x for bauxite
and 1.16x for laterite -- and it is refused rather than adopted.** That is the
one term the land water ledger has already settled: annual `P - E` is exactly
the water a cell drained over a closed cycle, and clamping it per bin counts
the wet season's supply twice. It is measured and reported beside the arms so
that the next reader does not have to re-derive that it was considered.

**The environmental lapse rate.** Fitting the annual-mean profile gives 6.799
K/km; the record-weighted mean of the twelve per-bin fits gives 6.763, a
difference of -0.036 K/km or 0.53%. It is one-signed in practice -- the
annual-mean profile fits steeper in 937 of 1,019 land cells -- but a slope is
not a convex functional of the profile, so the sign was not predictable in
advance the way the exponential's was. The per-bin rate varies by 1.21 K/km
within a cell across the twelve bins, 34 times the correction and 18% of the
rate, so the correction is well inside this instrument's own scatter and is
declared at the site rather than applied.

## Two defects found on the way, of other classes

Both were live when this audit was written and both are now closed. They are
kept because the SHAPE of each is what a future reader needs, and because the
first is the reason the audit's own line numbers should not be trusted without
re-reading the source.

**A cross-module caller unpacked a prefix of a tuple that had grown.**
`surface_water.climate_fields` gained a staged albedo as an eighth return and
`build_groundwater.py` unpacked seven. All three call sites now unpack the whole
tuple, so a ninth return raises rather than shifting every name one place along.
The same file's `measure_carve_effect` was separately calling
`per_basin_forcing` with nine arguments where the signature takes ten, which
would have raised the moment the step ran.

**A bin INDEX where the bin centres belong.** `build_sea_salt.py` and
`build_volcanic_sulfate.py` called `climatology.bin_weights(np.arange(nbin))`,
and evenly spaced integers assert that the bins are even instead of measuring
them. Both now pass the file's own `time` variable, and `lib/climatology.py`
no longer has a branch that could answer: `_checked_axis` refuses integers from
zero outright, and `infer_ntimes` refuses an evenly spaced axis whenever the
consistent record counts disagree about the weights. The route is closed at
both ends, and `scripts/smoke_test.py` holds the gate.
