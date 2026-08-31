# The absent tile was undefined on every partial cell, and the seams are ordered

**Measured:** 2026-08-30, on `canonical-10m-carve2`, terrain hash
`f496ae9fd749`, at T21, from the restart of a one-orbit cold start
(`run_71b36de31f0c`) written by the carrier-only executable. The staged
fractional support is `exoplasim/inputs/t21/orogen_T21_surf_1720.sra`.

This is worldbuilding. Vesper is an invented super-Earth; the quantities below
are properties of that planet's modelled surface and of the Fortran that
advances it. Nothing here concerns the real world.

SPAT-5 selected a conservative tile model and landed its carrier: surface code
1720 reaches the model as separately owned, immutable `dlf`, `ylf` and `xlf`,
and the per-tile boundary and exchange bundles exist and are restart records.
The physics is still binary. This note records two things found while sizing
the next slice: the bundles were undefined on exactly the cells they exist
for, and the six remaining implementation seams have a dependency order that
was not stated.

## Only one surface owner advances on a cell, and neither initialises off-mask

`landmod` gates its physics on `dls(:) > 0.0` and `seamod` gates all of
`seastep` on `dls(:) < 0.5`. `landmod.f90:1022` binarises `dls` to exactly 0
or 1, so the two masks are complementary and exhaustive: exactly one owner
advances on any cell.

Neither owner initialises outside its own mask either. `landini` sets soil
state under `if (dls(jhor) > 0.0)`; `seaini` sets `dts`, `dsst`, `dmld`,
`dicec`, `diced` under `where(dls(:) < 0.5)`.

The two owners keep separate arrays -- `landmod::dts` and `seamod::dts` are
different variables, both `NHOR` -- so the storage for a tile model is already
there. What is not there is a value. `landmod::dts` and `landmod::dqs` are
declared `= 1.0e20`, a deliberate sentinel; `seamod::dts` and `seamod::dqs`
carry no initialiser at all. Every other archived boundary field -- `drhs`,
`dz0`, `dalb`, `dsalb` -- is a shared `plasimmod` array and is always defined.

## The completion covered the pure cells and only the pure cells

`spat5_complete_absent_tile` exists to give the tile that did not advance a
finite copy of the one that did, because the flux and radiation kernels
evaluate every lane before combination. It was keyed on the fraction's
endpoints, `dlf == 0.0` and `dlf == 1.0` -- which is the complement of the
population it is for.

Read out of `run_71b36de31f0c`'s restart at T21, against the staged code-1720
field:

  pure land cells      586
  pure ocean cells     409
  partial cells       1053

  `dlt_ts` == 1.0e20   620 cells, all partial, all ocean-owned by `dls`
  `dot_ts` == 0.0      433 cells, all partial, all land-owned by `dls`

620 + 433 = 1053, and no pure cell carried either value. The partition is
exact: the sentinel and the unassigned zero sat on precisely the cells whose
existence motivated the tile decision, and both were written into the restart
under a name that says they are tile state.

The values were inert. Nothing in the current physics reads the boundary
bundles, which is why the carrier slice's A/B came back byte-identical and why
this is not a correction to any published number. What it is, is the trap the
next slice steps on: `spat5_tile_turbulent` reads exactly these arrays, on
exactly these cells, and the build carries `-ffpe-trap=invalid,zero,overflow`,
so a 1.0e20 surface temperature entering the flux kernel is an overflow rather
than a wrong answer.

The repair keys the completion on `dls`, the execution mask, which is what
"the tile that did not advance" actually means, and runs it on continuations
as well as cold starts -- `fluxstep` and `radstep` precede the first
`surfstep`, so a restart written before the fix would otherwise carry the
sentinel live for one timestep. `analysis/partial_surface_gate.py` holds both
properties.

## What the repair moved, measured

The same one-orbit T21 cold start on the rebuilt executable
(`run_ef6fa5447420`), against `run_71b36de31f0c`:

  restart records byte-identical      236 of 244
  restart records differing             8, and all eight are tile bundles:
                                        `dlt_ts`, `dlt_qs`, `dlt_z0`,
                                        `dlt_alb`, `dlt_sa1`, `dlt_sa2`,
                                        `dot_ts`, `dot_qs`
  climate output variables differing    0 of 99, bit-identical

  sentinels in any tile bundle          0
  unassigned zeros in any tile bundle   0
  max |land tile - ocean tile|          0, on every cell and every field

No prognostic or diagnostic record moved. The last line is the invariant this
slice is allowed to hold and the next slice is required to break: while the
absent tile is a copy, the two bundles agree exactly, so any area weighting of
them returns the binary answer. When `spat5_tile_state` lands, that zero
becomes a real difference on the 1053 partial cells and the flux bracket in
`analysis/coastline_flux_bracket.json` becomes a prediction to test rather
than a bound.

## What the land tile would start from, read off the staged fields

`config/partial_surface.yaml` declares eight initial-state fields and the
decision gate now reads the ones that reach the model through a staged `.sra`
rather than trusting their `status`. The population is the cells the binary
ownership mask (code 172) gives to the ocean while the land fraction (code
1720) is positive: 620 of the 1053 partial cells, matched by index off two
fields on the same Gaussian grid.

  code 173 land roughness        0.3446..3.5809 m over the 620, 610 distinct
  code 229 field capacity        0.01654..0.26549 m over the 620, 605 distinct
  code 174 background albedo     0.06 on all 620, one distinct value
  codes 175/176 albedo bands     0.07028 and 0.05364 on all 620, one each
  code 212 vegetation cover      0 on all 620

against 0.003..0.7183, 0.01991..0.25679, 0.10381..0.26049, 0.072..0.197,
0.11774..0.29979 and 0..0.5 over the 586 pure-land cells. Roughness and field
capacity are integrated from the mesh and vary cell to cell; the albedo pair
and the canopy fraction are one value each, which is the open-water constant
the binary mask wrote. No staged field carries a sentinel and no cell with a
positive land fraction is unpopulated: the roughness builder's own land
population is 1639 cells, which is exactly 586 pure land plus 1053 partial.

**Those 620 cells hold 11.578% of the mesh land area binned to the grid.** That
is the size of what the two blocked fields cost, and it is why the seam order
refuses the model change rather than merely preferring the fields first: a land
tile switched on today would begin advancing an ninth of this world's land at
open-water albedo under no canopy, and the shortwave error would be attributed
to the tile model.

### The relief statistic thins on those cells, and the thin part is immaterial

`cell_moments` returns a region count beside the variance because the count is
what says whether a spread means anything, and `build_surface_roughness.py`
reports its median without gating on it. Over the 620 the count runs from 1,
with a fifth percentile of 4, against a minimum of 565 over pure land and 335
over the land-owned partial cells: the sparse support is exclusive to the cells
the land tile newly advances on.

  fewer than 30 land regions      112 cells, 0.0279% of all land area
  exactly one land region          10 cells, relief sigma exactly 0, so the
                                   roughness is that one region's cover class
                                   verbatim, 0.6251..2.9698 m
  z0 above the pure-land maximum  212 cells, 2.0547% of all land area, median
                                   land-region count 112

The sparse tail is real and is not worth a repair: 0.028% of land area cannot
move a surface-energy term against the 0.12 W m-2 state-storage tolerance the
tile model was selected on, which is the same materiality test the decision
gate already applies. The 212 rough cells are not the sparse ones -- their
median support is 112 regions -- and a coastal cell whose land runs from sea
level to a summit genuinely holds more sub-grid relief than an interior
plateau, so the high tail is the orographic term working rather than failing.

## The remaining seams are ordered

The decision gate listed six missing seams as a set. They are not a set.

While execution is binary, the absent tile's boundary is a COPY of the present
one -- that is what the repair above makes true, and it is the honest state of
the model rather than a defect. It follows that an exchange seam taken before
tile state would evaluate the same surface twice and area-weight two copies of
it. The result is arithmetically the binary answer, reported through per-tile
diagnostics that are not tiles, on a seam that the gate would then credit as
implemented. The order is therefore:

1. `spat5_tile_state` in `landmod` and `seamod` -- both tiles advance wherever
   their fraction is positive.
2. `spat5_tile_restart` -- that state is serialised.
3. `spat5_tile_turbulent` and `spat5_tile_radiative` -- dual evaluation.
4. `spat5_tile_diagnostics` -- per-tile output.

`analysis/partial_surface_decision_gate.py` now carries the order, names the
next actionable seam, and fails if a seam is present while a prerequisite is
absent.

## What the first seam actually needs

The gating cost of `spat5_tile_state` is not the gate expression; it is the
initial state behind it. A land tile on an ocean-owned partial cell has no
soil temperature, water or snow, and an ocean tile on a land-owned partial
cell has no sea surface temperature, mixed-layer depth or ice. Those are
boundary-condition questions -- what a 90%-land cell's ocean tile starts at --
and they belong with `build_boundary_conditions.py`, which already emits code
1720, rather than inside the model.

## Every site the first seam has to move

Read off the model source. `landmod` gates on `dls(:) > 0.0` at thirty-odd
sites, strictly positive, so its form is already "advance wherever my share is
positive" and only the binarisation makes that equivalent to a binary mask.
That asymmetry is what makes the model side cheaper than it looks, and it has
ONE exception.

  landmod.f90:2777  `wetalb`, `if (dls(jhor) > 0.5)`. The only hard threshold
                    in landmod, and it is the wet-soil albedo: the land tile's
                    shortwave. It has to become `dlf(jhor) > 0.0` with the
                    rest, or the tile advances a soil moisture whose albedo
                    never responds on any cell the ocean owns.
  landmod.f90:1494  `dglac(:) > 0.5 .and. dls(:) > 0.0`. The 0.5 is the
                    GLACIER fraction, not the ownership mask. It is a separate
                    question and stays.

`seamod` has seven sites, all `dls(:) < 0.5`, and every one of them has to
become `dlf(:) < 1.0`:

  seamod.f90:150   `seaini`, the transfer of cts, cicec, ciced, csnow, csst
                   and cmld into the puma arrays. THIS IS THE INITIAL-STATE
                   SITE, and it is why world-cp27 comes first.
  seamod.f90:186   the cold-start sea-ice albedo ramp and sea saturation
  seamod.f90:258   flux accumulation into the ice coupler
  seamod.f90:278   the accumulated fluxes divided by the accumulation count
  seamod.f90:315   `icestep`'s return into dsst, dts, dmld, dicec, diced,
                   dsnow, dsmelt, dsndch
  seamod.f90:330   the dependent surface variables: dt(:,NLEP), the ice albedo
                   ramp, dqs, drhs, the Charnock roughness, dalb, dsalb, dz0
  seamod.f90:391   `seaalbarchive`, the albedo diagnostic radmod actually used

THE BINARISATION IS THREE SITES, NOT ONE. Each module rounds its own copy of
the mask right after reading it, so removing one leaves the others rounding
silently:

  landmod.f90:1022  `dls`
  oceanmod.f90:349  `yls`
  icemod.f90:357    `xls`

`analysis/partial_surface_decision_gate.py` checks all three and names the ones
that remain.

