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

The cheapest reading available, and it should precede any download.
`references/climber-x/` was pinned for PALADYN and SEMI, and most of the tree
was never entered.

| path | what it is | serves |
| --- | --- | --- |
| `src/ice_sico/` | SICOPOLIS, 27k lines, behind the 1.2k-line `src/ice/` adapter that also serves Yelmo | CLIM-62 route C and its separability question. STARTED: the adapter measurement is already in CLIM-62 |
| `src/geo/` | 21 files including `drainage_basins.f90`, `coast_cells.f90`, `connect_ocn.f90` | HYD's drainage and basins, OCN-11's coast and bathymetry contract, SURF-7's discharge mask. **The nearest peer to this project's hydrography that exists on this disk** |
| `src/ch4/` and `src/bnd/ch4*.f90` | the methane module and its radiative forcing path | WET-10, which has to close a methane flux, and BVOC-6 |
| `src/lndvc/`, `src/bnd/` | land vegetation carbon; boundary and forcing construction | BIO, EFOR |
| `references/socrates/` correlated-k tooling | the k-table generator, as opposed to the radiation code already read | CLIM-61, which prices a band-resolved scheme |

## Tier 1: open rows, no external comparison in the domain

**CaMa-Flood** -- global river routing with floodplain inundation from sub-grid
topography. *Row-driven and blind-spot-driven at once.* Serves WET-2's seasonal
inundation, ANUT-6's routing, SURF-7's discharge mask and GRID-2's sub-grid
hypsometry. The blind spot behind it is larger than any of those: nothing here
represents transient water storage in a channel network, so floodplain area --
which sets evaporation and methane-emitting area -- has no representation to be
wrong. Yamazaki et al.; freely available for research, licence to verify.

**CAABA/MECCA** -- a tropospheric chemistry box model. *Both.* BVOC-6 names
OH/O3/NO3 with no external reference and BVOC-10 closes a loop through them.
The blind spot: **methane lifetime on this world is unknown.** A K dwarf's
ultraviolet changes photolysis rates, so Earth's roughly nine years does not
transfer, and WET-10 has to close a methane budget without it. A box model is
small, portable, and answerable offline. GPL, part of the MESSy project.

**ArcSDM, or the USGS three-part assessment method** -- mineral prospectivity
mapping. *Blind-spot-driven.* MIN has 1 open row of 6 issued, which reads as a
finished domain and may instead be an unexamined one. The tension worth
surfacing: prospectivity methods are SUPERVISED -- weights-of-evidence and fuzzy
logic are calibrated against known deposit occurrences -- and this world has no
occurrences to train on. Whether `minerals/` is doing something defensible
without a training set, or is a plausibility field wearing a method's clothes,
is not currently asked anywhere. ArcSDM is open source.

**Landlab** -- a modular landscape evolution framework, one component per
process. *Row-driven.* GRAV-6 records that glacial erosion carries no gravity
term while Glen's law makes ice velocity go as `g^3`. LITH and PHYS-13 sit
beside it. Landlab is the most readable LEM and separates its erosion laws
cleanly enough to read one without the rest. **Read the no-time-axis rule first**
(`docs/src/reference/no-time-axis.md`): Orogen has no time axis and a LEM is
explicitly time-evolving, so the transferable thing is the FORM of a
gravity-scaled erosion law, not a rate. MIT licence.

**ESMF/ESMPy, or SCRIP** -- conservative regridding between unstructured and
structured grids. *Row-driven.* SPAT has 11 open rows of 11 issued and no
external eye, and rule 3 exists because mesh-to-grid transfers have silently
matched zero cells three times. ESMF is the reference implementation of
conservative remapping and of what "conservative" is allowed to mean. Its
first-order conservative and second-order options would give `lib/gridding.py`
something to be checked against that can fail.

## Tier 2: blind-spot probes, with no row to justify them

These are on the list BECAUSE nothing names them.

**A two-moment cloud microphysics scheme** -- PySDM for readability, or
Morrison-Gettelman as the production reference. Aerosol here is offline dust,
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

**FATES, or ED2** -- vegetation demography. *Blind-spot-driven.* DEMO has 6 open
rows of 6 issued and the only demographic model consulted is the one this
project runs. A cohort model's assumptions are hard to see from inside it, and
FATES is the modern independent implementation of the same idea.

**Global NEWS 2, or IMAGE-GNM** -- riverine nutrient export. *Row-driven.* ANUT
has 10 open rows of 10 and no external comparison. These are the standard global
models for how nutrients actually reach an ocean, which is ANUT-6's subject.

**SPITFIRE** -- process-based fire. *Row-driven.* FIRE has 10 open rows of 10
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
