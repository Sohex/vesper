# Audit: the grid convention across the export/model seam, and what runoff is averaged over

*Audited 2026-08-17, against `baseline_regular_climatology.nc` and
`coupling_exoplasim-T42.nc` on `precarve-craton`. Read-only audit; nothing here
was fixed. Every number is measured on those artifacts.*

Commissioned as a general sweep for gaps rather than against a specific
suspicion. It found one defect that invalidates the current carve list, one
missing invariant that would have caught it in a single line, one unbracketed
choice in the criterion's denominator, and one labelling mismatch that is the
reason the first defect looked like a fix.

Findings are tagged **[numeric]** where computed here, **[inspection]** where
read out of code or artifacts.

Findings 1, 2 and 4 are one defect and its consequences. Finding 3 is
independent and composes with it.

---

## 1. `basin_means` matches longitudes across the seam, and shifts every basin 180 degrees

**[numeric]** `hydrography/scripts/carve_verdict.py:basin_means` reconciles the
coupling matrix's columns with the climatology's by matching longitude labels:

    wrap = lambda a: (np.asarray(a) + 180.0) % 360.0 - 180.0
    remap = np.abs(wrap(field_lon)[None, :] - wrap(coupling_lon)[:, None]).argmin(axis=1)

`CLAUDE.md` rule 3 forbids exactly this: *"Never match by longitude across the
export/ExoPlaSim boundary. The two label the same grid differently. Map by
index."* The correct mapping is the identity, and the remap shifts by 64 of 128
columns, which is half a planet.

**The proof is the land mask, and it is exact.** The mask ExoPlaSim was given is
written by `lib/gridding.py:region_cells` as `col = (lon + 180)/360 * nlon`,
which is the same expression `build_hydrography.py:113` uses to place a mesh
region in a coupling column. So the two share an index convention by
construction, and the file the model read back is bit-identical to the mask it
was handed:

| `orogen_T42_surf_0172.sra` against the climatology's `lsm` | cells agreeing |
| --- | ---: |
| index for index | **1.0000** |
| shifted by 64 columns | 0.5955 |

The climatology's `lon` variable reads 0 to 360 because ExoPlaSim labels its own
grid and has never heard of Orogen. That label is not a coordinate
correspondence, which is what rule 3 is about.

**The physical check agrees and is the one to keep.** Endorheic catchments are
inland, so almost none of their area should fall on a cell the model calls ocean:

| mapping | catchment area on model-ocean cells |
| --- | ---: |
| index for index | **1.01%** |
| the remap in use | **49.77%** |

Half of every catchment integral is currently taken over open ocean.

**Size of the error on the verdict.** Recomputing with the identity mapping and
changing nothing else. The first row reproduces
`hydrography/analysis/carve_verdict.json` exactly, which is the check that the
recomputation is faithful:

| mapping | penman carves | wet carves | basins with no runoff |
| --- | ---: | ---: | ---: |
| the remap in use | 1539 | 1978 | 1643 |
| index for index | 1725 | 2559 | 1062 |

**1,172 of 3,621 basins get a different verdict, 32.4% of the catalogue.** 493
flip from carve to preserve and 679 from preserve to carve.

**The shift is independent of the evaporation scheme**, which is worth stating
because PHYS-4 landed between the first measurement of this and the one above.
Before the stability correction the same comparison read 1561 against 1758, with
1,187 basins differing; after it, 1539 against 1725 with 1,172. A change to how E
is estimated moves both mappings together and leaves the mapping error where it
was, so nothing downstream of Penman can absorb this.

**This is a recurrence, not a new defect.** `lib/orogen.py` records against the
`carved-zoned` hash that the build was superseded because *"every basin
integrated its antipode's rainfall"* and *"of its 1,522 carves, 850 are not
justified by the climate meant to justify them and cannot be un-cut"*. The remap
is the fix that was applied then. It was correct then:
`build_boundary_conditions.py` still went through the equirectangular PNG path,
whose convention differed from the coupling's. That script now integrates from
the native mesh on the coupling's own convention, and the correction outlived
the thing it corrected.

`basin_means` made `field_lon` a required argument so that omitting it could not
silently mean "do the wrong thing". The argument is supplied by every caller and
supplying it is the wrong thing.

**Reach, which is wider than the verdict.**

- `export_carve_list.py:299` calls `basin_means`, so `carve_list.txt` on disk
  carries the shift. That is the artifact that leaves the project.
- `surface_water.py:127` calls it too, and `region_grid_cells` at line 164
  repeats the remap for per-region sampling. Lake *areas* are therefore solved
  from antipodal climate. Lake *positions* are geometric and are unaffected.
- `surface_water.py:134` samples lake precipitation and evaporation at each
  basin sink with `col = mod(round((sink_lon - lon[0]) / dlon), nlon)`, which
  measures an Orogen longitude from the climatology's first label and lands in
  the same place, 65 columns out. This one SURVIVED the fix that closed
  HYD-13 and was closed under GRID-1 instead: 52.47% of basin sinks landed on
  cells the model calls ocean against 3.04% index for index. A third site,
  `maps/build_basemap.py:upsample`, was found by the same invariant at the same
  time. `notes/audits/carve-criterion-terms.md` finding 8.
- `exoplasim/inputs/t42/albedo_report.json` records lakes covering 11.28% of
  land and worth **-0.018024** on land-mean albedo, sourced from that
  `surface_water.nc`.

That last item is what makes this expensive. The baseline run's albedo fields
174 to 176 were built from lake extents computed on antipodal climate, so the
baseline climatology the verdict reads is itself contaminated. The remedy is
therefore a **re-run of the baseline** in the sense `WORKFLOW.md` section 0
defines: rebuild `surface_water.nc`, rebuild 174 to 176, run, rebuild the
climatology, re-take the verdict. `WORKFLOW.md` A2 requires the cycle run to be
centred on a *final* baseline, so this lands before CYC-1 and not after.

## 2. Nothing checks the coupling against the mask, only that `cell_lon` exists

**[inspection]** `scripts/check_consistency.py:223-233` opens each coupling
matrix, asserts the variable `cell_lon` is present, and stops. Its own docstring
names the antipodal bug as the reason the check exists:

> **Coupling convention.** `coupling_*.nc` must carry `cell_lon`. Without it the
> longitude convention cannot be checked [...] Every basin read its antipode for
> an entire iteration.

The convention is not checked. Only the field that would allow checking it is,
and that field is what the incorrect mapping is built from.

This is `CLAUDE.md`'s class 17 exactly, a check that cannot fail: presence of a
variable has no right answer to be wrong about. The table in finding 1 does have
one. The fraction of coupling area landing on cells the model calls ocean is
bounded above by coastal rounding, measured at 1.01% here, and 49.77% is not a
different opinion about a grid, it is a wrong answer.

## 3. Catchment runoff is an annual mean, which is the conserved quantity

**[numeric]** `carve_verdict.py` forms runoff as the catchment-area-weighted mean
of annual-mean `P - E`, clamped at zero after the aggregation. The clamp is
commented and deliberate: an unclamped negative drove the equilibrium lake area
to large negative values.

The **spatial** ordering is defensible. ExoPlaSim's `mkradv` advects runoff
downhill and a cell can re-evaporate water routed into it, so a negative local
`P - E` over land is physical, and the catchment integral is the net supply to
the sink provided the routing stays inside the catchment. 38.9% of land cells,
1,597 of 4,105, run annual-mean `P - E` below zero.

The **temporal** ordering is not. A catchment delivers water in its wet season
and delivers none in its dry one; it cannot deliver a negative amount. Averaging
over the orbit first charges the dry-season deficit against the wet-season
supply. Measured on the baseline's 12 bins, with the index mapping of finding 1
so the two compose:

| runoff treatment | basins with zero runoff | catchment-mean runoff |
| --- | ---: | ---: |
| annual mean, clamped after aggregation (in use) | 1062 | 161.2 mm/yr |
| clamped per cell, annual | 655 | 162.2 mm/yr |
| clamped per cell per bin | 6 | **266.1 mm/yr** |

The land-mean effect on runoff generation alone is 1.92x.

**Per-bin clamping is not an upper bound, and the annual mean is the conserved
answer.** Settled 2026-08-17 in `notes/audits/carve-criterion-terms.md` finding 4,
which supersedes the reading this table was given. Over one annual cycle at steady
state a land cell's storage returns to where it started, so `annual(P - E)` equals
the runoff that cell generated, exactly -- an identity, not an estimate. The
negative bins are the store being drawn down, and the water that refills it in the
wet season is water that did not run off, so adding the negative bins back counts
the wet season twice. The model's own soil water says so cell by cell: the ratio
of its seasonal range to the discarded deficit has a median of 0.988, quartiles
0.95 to 1.07.

So the criterion's denominator already uses the right quantity. The remaining
choice is where to clamp, and after aggregation is where the residual noise on a
non-negative quantity has averaged down furthest.

**The 1,062 basins with no runoff are dry catchments, not an averaging artifact.**
Under the same identity, a catchment whose evaporation consumes its precipitation
delivers nothing. They are auto-preserved by `index = inf` only in the RATIO form
of the test; `export_carve_list.py` decides them on the discharge form, which asks
whether the lake surface itself gains, and on the corrected evaporation 243 of
them do. What is worth chasing about that population is window length rather than
treatment: they concentrate at high latitude, where a five-orbit climatology is
least periodic.

## 4. The climate figures are labelled 180 degrees from the maps

**[inspection]** `exoplasim/scripts/analyze_climatology.py:57` relabels the
climatology's axis to -180..180 and reorders the field to match:

    display_lon = (lon + 180.0) % 360.0 - 180.0
    order = np.argsort(display_lon)
    return display_lon[order], field[..., order]

The figure is internally consistent, because `decorate` reorders the coastline
contour through the same function. But index 0 of the model grid is Orogen's
-180, so every longitude tick on `baseline_koppen_geiger.png` and its siblings
names a meridian 180 degrees from the same feature in `maps/`, which renders on
the export's own -180..180.

Nothing computed reads these figures, so no result is wrong because of it. It is
recorded here because it is the visible evidence for the false premise in finding
1: two products that share a grid appear to disagree about longitude, and
reconciling their labels is the wrong repair.

---

## What I would chase first

Finding 1, and specifically the invariant in finding 2 alongside it. The verdict
recomputation is cheap and the answer is already in this document; what is not
cheap is that the baseline itself has to be re-run, and what is not recoverable
is a carve list shipped to Orogen before that happens. Land finding 2 in the same
change, because this defect has now occurred twice with a fix in between, and the
distinguishing feature both times was that nothing tested the mapping against
something that could fail.

Finding 3 belongs in the same pass rather than a later one. It moves the same
denominator as `HYD-11`, `DUST-10` and `DUST-11`, and the baseline re-run that
finding 1 forces is the natural place to settle all of them at once rather than
re-running for each in turn. It was settled that way, and the answer was that the
denominator is already right; `notes/audits/carve-criterion-terms.md` carries it
along with the three items it was settled beside.

---

## Tasks

Tracked in `TASKS.md`, not restated here.

| finding | id |
| --- | --- |
| 1. the antipodal longitude remap | `HYD-13` |
| 2. the coupling convention is not checked | `CONS-2` |
| 3. seasonal rectification in catchment runoff | `HYD-14` |
| 4. climate figure longitude labels | `CLIM-8` |
