# External sources worth acquiring, and why each one

*Written 2026-08-22. Vesper is a fictional planet and this project simulates it.
This is a shopping list with reasons, not an acquisition record: what has
actually been fetched, pinned and read is `references/INDEX.md`, and what was
learned from it is `notes/external-model-survey.md`.*

## The rule for earning a place here

Two ways, and the second matters more.

**It serves an open row.** A row names a mechanism and an external
implementation of that mechanism exists. This is cheap to justify and easy to
act on.

**It probes a domain with no external comparison at all.** A rule keyed to open
rows cannot find a blind spot, because a blind spot is a thing nobody has
written a row about. Survey section 9 records that discovery and what correcting
for it turned up. So this list deliberately includes domains where the project
has open rows and no external eye, and one where it has almost no rows at all --
which is itself the signal.

Entries say which of the two they are. Licences below are what the projects are
generally understood to carry and **must be verified at acquisition**, not
taken from here; `references/INDEX.md` is where a verified licence gets
recorded.

## Tier 0: already on disk, never opened

**All five rows below were read on 2026-08-22 and are closed.** Survey sections
42, 44, 45, 46 and 47 carry the findings.

The cheapest reading available, and it should precede any download.
`references/climber-x/` was pinned for PALADYN and SEMI, and most of the tree
was never entered.

| path | what it is | serves |
| --- | --- | --- |
| `src/ice_sico/` | SICOPOLIS, 27k lines, behind the 1.2k-line `src/ice/` adapter that also serves Yelmo | CLIM-62 route C and its separability question. DONE 2026-08-22, survey section 45 |
| `src/geo/` | 21 files including `drainage_basins.f90`, `coast_cells.f90`, `connect_ocn.f90` | HYD's drainage and basins, OCN-11's coast and bathymetry contract, SURF-7's discharge mask. **The nearest peer to this project's hydrography that exists on this disk.** DONE 2026-08-22, survey section 47 -- and note the correction: section 7a had already triaged this directory at FILE level, so it was unread rather than unopened. `gia.f90`, `vilma.F90`, `fastearth.F90`, `sed.f90`, `q_geo.f90`, `topo_filter.f90` and `geo_grid.f90` were absent from that triage entirely |
| `src/ch4/` and `src/bnd/ch4*.f90` | the methane module and its radiative forcing path | WET-10, which has to close a methane flux, and BVOC-6. DONE 2026-08-22, survey section 42 |
| `src/lndvc/`, `src/bnd/` | land vegetation carbon; boundary and forcing construction | BIO, EFOR |
| `references/socrates/` correlated-k tooling | the k-table generator, as opposed to the radiation code already read | CLIM-61, which prices a band-resolved scheme. DONE 2026-08-22, survey section 44 |

## Tier 1: open rows, no external comparison in the domain

**CaMa-Flood** -- ACQUIRED 2026-08-22, `references/cama-flood/`. Global river routing with floodplain inundation from sub-grid
topography. *Row-driven and blind-spot-driven at once.* Serves WET-2's seasonal
inundation, ANUT-6's routing, SURF-7's discharge mask and GRID-2's sub-grid
hypsometry. The blind spot behind it is larger than any of those: nothing here
represents transient water storage in a channel network, so floodplain area --
which sets evaporation and methane-emitting area -- has no representation to be
wrong. Yamazaki et al.; freely available for research, licence to verify.

**CAABA/MECCA** -- ACQUIRED 2026-08-22, `references/caaba-mecca/`, release 4.6.0. A tropospheric chemistry box model. *Both.* BVOC-6 names
OH/O3/NO3 with no external reference and BVOC-10 closes a loop through them.
The blind spot: **methane lifetime on this world is unknown.** A K dwarf's
ultraviolet changes photolysis rates, so Earth's roughly nine years does not
transfer, and WET-10 has to close a methane budget without it. A box model is
small, portable, and answerable offline. GPL, part of the MESSy project.

**Note on the two cGENIE "lite" accelerators**, which an earlier draft of this list expected to be worth reading: they are not. Survey section 46c records that `ocnlite` is an empty template and `goldlite` an abandoned prototype that has never been run.

**ArcSDM** -- ACQUIRED 2026-08-22, `references/arcsdm/`. Mineral prospectivity
mapping. *Blind-spot-driven.* MIN has 1 open row of 6 issued, which reads as a
finished domain and may instead be an unexamined one. The tension worth
surfacing: prospectivity methods are SUPERVISED -- weights-of-evidence and fuzzy
logic are calibrated against known deposit occurrences -- and this world has no
occurrences to train on. Whether `minerals/` is doing something defensible
without a training set, or is a plausibility field wearing a method's clothes,
is not currently asked anywhere. ArcSDM is open source.

**Landlab** -- ACQUIRED 2026-08-22, `references/landlab/`, AND IT DOES NOT SERVE THE ROW IT WAS FETCHED FOR.
Of 61 components none is glacial, so GRAV-6's Glen's-law question has no source
in any tree held here and that gap stands open. Re-keyed on inspection to what
it does carry: `priority_flood_flow_router` as a direct comparison for
`hydrography/drainage.py`, `flexure` for the isostasy blind spot, `lithology`
for LITH, and fluvial erosion laws whose gravity dependence is implicit through
shear stress and is a different question from GRAV-6's. *Row-driven.*

**Read the no-time-axis rule first** (`docs/src/reference/no-time-axis.md`):
Orogen has no time axis and a landscape evolution model is explicitly
time-evolving, so what transfers is the FORM of a law, never a rate.

**GRAV-6 IS NOW AN UNSERVED ROW.** Its subject is glacial erosion under Glen's
law, and no tree held here carries a glacial erosion rule. SICOPOLIS supplies
the ice DYNAMICS and survey section 45a prices gravity through it, but abrasion
and quarrying laws are a separate literature. If that row is to be served, it
needs its own acquisition and this list does not currently contain one.

**ESMF** -- ACQUIRED 2026-08-22, `references/esmf/`. Conservative regridding between unstructured and
structured grids. *Row-driven.* SPAT has 11 open rows of 11 issued and no
external eye, and rule 3 exists because mesh-to-grid transfers have silently
matched zero cells three times. ESMF is the reference implementation of
conservative remapping and of what "conservative" is allowed to mean. Its
first-order conservative and second-order options would give `lib/gridding.py`
something to be checked against that can fail.

## Tier 2: blind-spot probes, with no row to justify them

These are on the list BECAUSE nothing names them.

**PySDM** -- ACQUIRED 2026-08-22, `references/pysdm/`. Super-droplet cloud microphysics, taken over a production two-moment scheme for readability. Aerosol here is offline dust,
sea salt and volcanic sulfate, and the clouds are diagnostic. **Nothing connects
condensation nuclei to droplet number to cloud albedo**, so the aerosol indirect
effect is not wrong here, it is absent, and no row says so. Large; scope before
committing.

**USGS Astrogeology planetary cartography standards, and the IAU nomenclature
conventions** -- *cartography and coordinates.* `maps/` is explicitly pictures
rather than analysis, which is a fair reason for it to have no rows. The reason
to read anyway is one directory over: **this project already has a documented
longitude hazard serious enough to be rule 3**, and planetary cartography has
settled conventions for exactly the questions behind it -- planetographic
against planetocentric latitude, positive-east against positive-west longitude,
where a prime meridian comes from on a body with no Greenwich, and how a
graticule is labelled when the body is not Earth. Whether this project's
conventions agree with any of them has never been asked. Documents rather than
source, and cheap.

**FATES** -- ACQUIRED 2026-08-22, `references/fates/`. Vegetation demography. *Blind-spot-driven.* DEMO has 6 open
rows of 6 issued and the only demographic model consulted is the one this
project runs. A cohort model's assumptions are hard to see from inside it, and
FATES is the modern independent implementation of the same idea.

**Global NEWS 2 and IMAGE-GNM** -- riverine nutrient export. NO PUBLIC SOURCE for either, checked 2026-08-22, so the papers are the acquisition rather than a tree: Beusen et al. (2015) `10.5194/gmd-8-4045-2015` for IMAGE-GNM and Mayorga et al. (2010) `10.1016/j.envsoft.2010.01.007` for Global NEWS 2. *Row-driven.* ANUT
has 10 open rows of 10 and no external comparison. These are the standard global
models for how nutrients actually reach an ocean, which is ANUT-6's subject.

**SPITFIRE** -- ACQUIRED 2026-08-22, and its host tree was then PROMOTED to a
vendored subtree at `vendor/lpjml/` because the work intended on LPJmL is
fork-shaped and it is planned into the pipeline, so the fire reading is now
incidental to why the tree is here. Process-based fire. *Row-driven.* FIRE has 10 open rows of 10
and the only scheme consulted is BLAZE, which arrived with LPJ-GUESS. SPITFIRE
is the independent scheme in the same lineage and reaches a different answer
about ignition and spread.

## Deliberately not on this list

**More exoplanet GCMs.** Survey sections 9 and 15 read ExoRT, ExoCAM and Isca
and closed the question. Another one is unlikely to change an answer.

**More EMICs.** CLIMBER-X, cGENIE and BIG-MITgcm cover the class, and section 6
records what that closed.

**A second ocean biogeochemistry model.** BIOGEM and MARBL were compared in
section 16 and OCN-4 has what it needs to choose.

**GPlates and plate reconstruction.** Orogen generates the geography and has no
time axis; a reconstruction engine has nothing to reconstruct here.

## Acquiring one

`references/INDEX.md` states the procedure and the reason for it: extract from a
PINNED tarball rather than clone, record the revision and the licence in the
table, and keep the bulk out of git. A tree acquired without its pinned revision
recorded is not reproducible once a scratchpad is cleared, which is the failure
that procedure exists to prevent.
