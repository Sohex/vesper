# Where external data comes from

A record of routes into data this project does not generate, so the next
person needing a dataset starts from a list rather than from a search engine
-- kept because a validation once stalled on an API quota that a bulk mirror
would have made irrelevant.

**Check the AWS Registry of Open Data first**, at `registry.opendata.aws`. It
indexes hundreds of public datasets held in S3 with no credentials and no quota:
ERA5 and other reanalyses, satellite archives, elevation, land cover, genomics.
Anonymous access is `fsspec.filesystem("s3", anon=True)` or a plain HTTPS GET
against `https://<bucket>.s3.amazonaws.com/<key>`, and listing is
`?list-type=2&prefix=...&delimiter=/`.

The reason to look there before an API is that an API is rate limited, its
answers depend on when you asked, and a validation whose result depends on which
points were throttled that minute looks like a measurement. A bucket is neither.

## What has been checked, and what it turned out to be

**Open-Meteo, `s3://openmeteo`.** Their own mirror of the archive their API
serves, in the `.om` format, readable by HTTP range request through the `omfiles`
package rather than by downloading a chunk. That matters: an ERA5 chunk is 95 MB
for 21 days, shaped `(721, 1440, 504)` as latitude by longitude by hour with a
chunk shape of `(1, 6, 504)`, so one row and six columns is a single fetch.

**It is a rolling window, not the archive.** Checked 2026-08-17:
`data/copernicus_era5/precipitation/` held 81 chunks numbered 904 to 984, which
is late 2021 to mid-2026, and some variables carry only 38. There is no
1940-onward history, so it cannot serve a 1991-2020 normal and the API remains
the only route to one. Do not re-derive this; check the chunk range if you need
a period.

**Open-Meteo's archive API**, `archive-api.open-meteo.com`. ERA5 by point, no
key, and an HOURLY request quota that a thirty-year daily series exhausts in
about 130 calls. It returns HTTP 429 with "Hourly API request limit exceeded",
which no amount of spacing inside the hour will fix. Cache every response and
make a re-run resume rather than restart;
`pedology/scripts/validate_against_earth.py` and
`hydrography/scripts/validate_lake_solver.py` both do.

**HydroSHEDS**, `data.hydrosheds.org`. HydroBASINS and HydroLAKES as plain HTTPS
downloads. It refuses `HEAD` outright and returns 403 without a browser
`User-Agent`, so a naive probe reports the data as unavailable when it is not.

**Copernicus DEM**, `s3://copernicus-dem-90m` and `copernicus-dem-30m`, no
credentials. One-degree COG tiles, 4.2 MB each at 90 m, named
`Copernicus_DSM_COG_30_N37_00_E095_00_DEM`. Read a tile straight out of S3 with
`fsspec` and `rasterio.io.MemoryFile` rather than downloading it. Used for HYD-7
to measure what representing a real basin at mesh scale costs its storage.

**GLHYMPS**, Gleeson et al. (2014), on Borealis (a Dataverse instance) at
`doi:10.5683/SP2/DLGXYO`. `GLHYMPS.zip` is 1.15 GB, CC-BY 4.0, and the file list
and sizes come from the Dataverse API without credentials:
`https://borealisdata.ca/api/datasets/:persistentId/?persistentId=<doi>`. The
2011 companion compilation is `doi:10.5683/SP2/TTJNIU`. This is the permeability
GW-3 would run an Earth calibration on.

**Berghuijs et al. (2022) global recharge**, Zenodo `10.5281/zenodo.7611675`,
CC-BY 4.0, `RechargeTotal.nc` at 6.2 GB in mm/yr with a companion recharge
fraction. RECHARGE rather than precipitation, which is the quantity a water
table model needs and the one that is hard to find. Its record states exclusions
by temperature and by aridity class, and WHICH side of the aridity threshold is
excluded has not been read off the paper; if it drops arid regions it is not
usable for this project, whose interesting basins are all arid. Check that
before downloading 6.2 GB. Moeck et al. (2020) at
`opendata.eawag.ch/dataset/globalscale_groundwater_moeck` is a point compilation
of measured recharge rates as CSV, which is small and is the cross-check rather
than the field.

**Fan et al. (2013) water table depth: THE HOST IS GONE.** Checked 2026-08-20.
`glowasis.deltares.nl` does not resolve at all -- DNS NXDOMAIN, while
`deltares.nl` itself resolves, so the GLOWASIS subdomain has been decommissioned
rather than moved within the site. That URL is what the paper names, and it is
also what every independent trail still points at: the HESS 2019 Amazon paper's
own data-availability section and Zeng et al. (2018) both cite the same dead
catalogue. The other routes were checked and are closed too: the Science
supplement carrying Databases S1 to S3 returns HTTP 403 on all three URL forms,
Fan's Rutgers page is 404, HydroShare returns no matching resource, Zenodo
carries nothing under the title, and the Borealis "Groundwaterscapes" dataset
Fan co-authors holds analysis code and a classification GeoTIFF but no water
table compilation.

This blocks GW-3, and it is worth being exact about why a substitute will not
do. The 1,603,781 well observations are the only EXTERNAL check available on
this project's water table; every other check it passes is an identity, a
conservation law or a reduction. Fan's SIMULATED equilibrium water table would
not serve even if a mirror turned up, because her model is the exponential
depth-decay formulation this project abandoned, so agreement or disagreement
would confound the formulation with the implementation. What would unblock it is
the compilation itself, from the authors -- Miguez-Macho is at Santiago de
Compostela -- or from wherever Deltares moved GLOWASIS.

**Zenodo**, for datasets published with a paper. The record API,
`https://zenodo.org/api/records/<id>`, lists files and sizes without
authentication, which is worth checking before starting a multi-gigabyte
download.

**LMD Generic PCM**, `https://web.lmd.jussieu.fr/~lmdz/planets/generic`. Not a
dataset but a second exoplanet GCM, fetched 2026-08-18 and living outside this
repo at `~/git/generic_pcm`. Source is `LMDZ.GENERIC` plus `LMDZ.COMMON`, 32 MB
and CeCILL licensed, so it is read for DESIGN and cited, and no code moves across.
The 338 MB `datagcm.tar.gz` unpacks to 1.2 GB and is the part worth having: it
holds `corrk_data`, `stellar_spectra`, `aerosol_properties`, `continuum` and
`surface_data`. The dated `datagcm_*.tar.gz` snapshots beside it are historical
copies of the same thing and were not fetched.

What it is good for here, and this is the point of writing it down:

- **`stellar_spectra/BT-Settl_stellar_spectra_grid/`** is the same grid `k25v` was
  interpolated from, 5,803 files at R = 100 across Teff, log g and metallicity.
  Teff 4900 and 5000 both exist and this star's derived log g
  (`build_stellar_spectrum.py`) snaps to the 4.5 grid point, so log g 4.5 at
  solar metallicity is the comparison point. That makes the band-1 flux fraction
  an independently checkable number rather than one this project computed alone.
- **`corrk_data/N2-0.000376CO2-H2Ovar_2026`** and its siblings up to 0.95 CO2 are
  correlated-k tables from modern line lists for an N2 atmosphere with variable
  water. This planet's atmosphere is N2-dominated with the CO2
  `config/planet.yaml` declares, so these bracket it, and a
  correlated-k treatment is an independent route to the shortwave CO2 and water
  vapour absorption that PHYS-1 and PHYS-6 derived from Howard's 1950s band data.
- **`aerosol_properties/optprop_dustvis_n50.dat`** and `optprop_dustir_n50.dat` for
  dust optics, against the DUST-12 choice.
- Its physics tree has NO wet deposition and no deposition velocity, checked by
  grep 2026-08-18. That is useful negatively: a second independent GCM has the same
  gap ExoPlaSim does, so DUST-7's scavenging is a common omission rather than a
  peculiarity, and there is nothing to borrow.

## The GW-3 Earth calibration: four inputs, all live

Checked 2026-08-20. The water table solver has passed only identities,
conservation laws and reductions, none of which can say it is right about a real
water table. Fan et al. (2013)'s well compilation is the external test, and its
own host is gone: `glowasis.deltares.nl` is NXDOMAIN while `deltares.nl`
resolves, so the subdomain is decommissioned rather than moved, and the Science
supplement 403s on every URL form. The exact Database S1 filenames are in
`references/INDEX.md` and return nothing anywhere indexed.

**What recovered the trail was the JRC AquaKnow mirror**, `aquaknow.jrc.ec.europa.eu/node/18916`,
which carries Fan's supplementary. That gives the schema to reproduce -- latitude,
longitude, land elevation in m above sea level, water table depth in m below land
surface -- and links a live THREDDS at `thredds-gfnl.usc.es/thredds/catalog/GLOBALWTDFTP`,
Miguez-Macho's own institution, holding continental annual and monthly means of
the SIMULATED field. Model output only; the observations are not under that tree.

So the compilation is rebuilt from the national archives it was made from.

**Recharge: Berghuijs et al. (2022), Zenodo, CC-BY 4.0, 6.2 GB in mm/yr.**
`10.1029/2022GL099010`. This was the gate, on the worry that it might exclude
arid regions, and the worry was backwards: 99% of its observation-based recharge
values come from regions with climate aridity above 0.75, and the vast majority
originate from Australia. What it excludes is PERMAFROST, "regions that can have
permafrost (i.e., mean temperature below -2C)", and what it is thin on is the
very wet -- Congo, Amazonia, southeastern Asia fall outside its observational
range. That is the opposite weighting from the one feared and the right one for
this world.

**United States wells: `api.waterdata.usgs.gov`, no key.** The old
`waterservices.usgs.gov/nwis/gwlevels` API was frozen in November 2025 and
decommissioned in February 2026; it 301s to a decommissioning notice. The
replacement is an OGC API collection:

    https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items?parameter_code=72019

**Parameter code 72019 is depth to water level in feet below land surface**,
which is Fan's quantity directly rather than a head needing conversion. Records
carry a monitoring location id, value and unit, time, vertical datum and approval
status.

**Australia: obtained, 75,321 sites.** The route is the Australian Groundwater
Explorer's per-state download, not an API, and the path to it took two wrong
turns worth recording.

    https://hosting.wsapi.cloud.bom.gov.au/arcgis/rest/services/groundwater/Bore_Water_Levels/FeatureServer/0/query

That FeatureServer holds **251,681 bores**, no key, and carries per bore exactly
what the schema needs plus two things Fan's four columns do not: `landelev` with
its datum and survey method, and `ftypeclass`, the bore's PURPOSE -- "Monitoring"
against production -- which is the anthropogenic-drawdown filter this whole
rebuild exists to make possible. It also carries `water_cnt`, `water_dmin` and
`water_dmax`, so the sampling of every bore is known.

**No public API serves the level VALUES.** The FeatureServer has one layer and
no related table; the siblings are Aquifer_Boundaries, Bore_Purpose,
Bore_Hydrochemistry and Bore_Salinity. `data.gov.au`'s NGIS package offers only
portal links. BOM's SOS2 endpoint returns an empty body to a GetCapabilities.
NGIS is documented as excluding time series.

**The levels ship BESIDE the geodatabase, not inside it.** Each state's zip
contains `gdb_<STATE>/ngis_<STATE>.gdb` AND, as siblings, `level_<state>.csv`,
`salinity_<state>.csv`, `hydrochem_*.csv` and `gw_state_README.txt`. The
geodatabase's nine layers are bores, logs, aquifers and management zones, and
`NGIS_Bore`'s `WaterLevel`, `WaterCount`, `WaterDateMin` and `WaterDateMax` are
flags and counts rather than values. Looking only inside the `.gdb` reads as the
data being absent when it is one directory level away, and that costs a download.
60,020,336 level readings across the eight states.

**Three datums, and they are not interchangeable.** `obs_point_datum` labels
every reading: `DTW` is depth below the ground surface and is Fan's quantity
exactly; `RSWL (mAHD)` is water level elevation on the Australian Height Datum;
`SWL` is depth below the bore's own reference point. DTW needs no elevation at
all, which matters because `LandElev` is populated for 10.8% of Western
Australian bores and none of Tasmania's.

**The casing stickup is `TsRefElev - LandElev`, and NOT `RefElev - LandElev`.**
`RefElev` is a copy of `LandElev` in this data, so using it applies no correction
at all and SWL reads 0.570 m too deep in the median on 91.8% of paired readings.
`TsRefElev` is the time-series reference and gives a median stickup of 0.633 m,
5th to 95th percentile 0.13 to 1.08 m. It is populated for about 4% of bores, so
SWL is dropped where it is absent rather than assumed zero.

**Both conversions are checked against an identity that can fail**, which is what
licenses them: where one bore reports two datums on one date, `DTW + RSWL` must
equal `LandElev` and `SWL - stickup` must equal `DTW`. Median residuals +0.000 m
over 238,185 pairs and -0.004 m over 156,037.

**The quality flag is reported, not filtered on.** Codes run A to F, are defined
neither in the shipped README nor on the Explorer metadata page, and are
per-agency: every Western Australian reading is `quality-A` and no Australian
Capital Territory or Tasmanian reading is. Filtering on A selects which STATE
survives, which is a selection bias wearing a quality label. Keeping everything
and reporting the flag turns it into a measurement: restricting to `quality-A`
moves a site's depth by a median of 0.0000 m with a 90th percentile of 0.075 m,
so it does not matter, and that is now known rather than assumed.

`hydrography/scripts/build_earth_wtd_sites.py` does the assembly and writes its
provenance. **BOM overwrites the Explorer extract about twice a year**, so the
result is tracked rather than ignored: it is a snapshot of a moving source and a
future download will not reproduce it.

One incidental finding worth keeping: `www.bom.gov.au` 403s a bare `curl` and
returns 200 to an ordinary browser User-Agent. It is a bot filter on the header,
not a restriction on the data, which BOM publishes for download under CC-BY.

**Take Australia first even though it is small by well count.** Fan's set is 88% United States and Canada, which are
humid-to-semiarid mid-latitudes, and this world is dominated by arid endorheic
basins: a calibration covering most of the wells and none of the regime is
testing the wrong thing. Australia is the arid analogue, and at 251,681 records
it is three times the 83,617 Australian sites Fan had.

**Canada is fragmented and is the lowest priority.** There is no single national
well database as there is for the United States; the holdings are provincial. It
is also the same humid-temperate regime the United States already covers, so it
adds well count rather than regime coverage.

Topography is Copernicus DEM 90 m from the AWS registry, the route HYD-7 already
used, and permeability is GLHYMPS at Borealis `10.5683/SP2/DLGXYO` under CC-BY.

**Two things decide whether the comparison is valid, and both are recorded
before it runs.** Compare the model AT OBSERVATION LOCATIONS, never
model-everywhere against observations-everywhere: Fan states the bias herself,
her model reading deeper than observations in arid regions because wells sit
where people live, which is valleys and oases. And filter production wells if the
metadata allows, because Fan's model neglects pumping, irrigation and drainage
while her observations do not -- this world has none of them, and rebuilding from
primary sources is what makes that filter possible at all, since NWIS carries
well-use metadata and Fan's four columns do not.

## A private gcc 16.1.1, from the pacman cache

Not a dataset but a build tool, and it is here because the route into it is the
same kind of thing: something this project needs, does not generate, and would
otherwise be re-discovered by search. Built 2026-08-18.

**Why it exists, and why nothing uses it.** A system update on 2026-08-18 took
`gcc-fortran` from 16.1.1 to 16.2.1 at about the same time as ExoPlaSim orbits
stopped completing. The two were read as cause and effect and this prefix was
built to get the old compiler back. That reading was WRONG. The orbits were
never slow, they were DEADLOCKED, on every compiler equally, because `nlowio` is
not broadcast to the non-root MPI tasks; `notes/audits/nlowio-collective-deadlock.md`
carries the evidence. One full T42 orbit takes 77.7 s built with 16.1.1 and
78.2 s built with 16.2.1, and the PRE-update binary deadlocks exactly like the
post-update one. Nothing in this project needs this prefix, and the model is
built with the system compiler.

The recipe is kept because a real toolchain regression will want it one day.
What follows is that recipe, verified on 2026-08-18. The version in the heading is the one that
happened to be cached then; check `/var/log/pacman.log` for what is actually
available rather than assuming.

**Where it comes from.** `/var/cache/pacman/pkg/`, which still holds the
packages pacman replaced. `/var/log/pacman.log` says which version was in place
when the fast binaries were built: the line to look for is the `upgraded
gcc-fortran` entry and the version on the left of the arrow. On 2026-08-18 that
was `16.1.1+r595+g171d15ac6959-1`, and several other 16.1.1 builds are cached
alongside it, so match the full version string rather than the `16.1.1` prefix.

**Four packages are needed, not three.** `gcc-libs` at this version is an empty
metapackage carrying only `.PKGINFO`, because this distribution splits the
runtime libraries into their own packages. The Fortran runtime comes from
`libgfortran` and the unwinder from `libgcc`, and without those two the prefix
links against the 16.2.1 copies in `/usr/lib`. Extract with `bsdtar`, never with
`pacman -U`:

    mkdir -p ~/toolchains/gcc-16.1.1
    cd /var/cache/pacman/pkg
    bsdtar -x -f gcc-16.1.1+r595+g171d15ac6959-1-x86_64_v4.pkg.tar.zst -C ~/toolchains/gcc-16.1.1 --strip-components=1 usr
    bsdtar -x -f gcc-fortran-16.1.1+r595+g171d15ac6959-1-x86_64_v4.pkg.tar.zst -C ~/toolchains/gcc-16.1.1 --strip-components=1 usr
    bsdtar -x -f libgcc-16.1.1+r595+g171d15ac6959-1-x86_64_v4.pkg.tar.zst -C ~/toolchains/gcc-16.1.1 --strip-components=1 usr
    bsdtar -x -f libgfortran-16.1.1+r595+g171d15ac6959-1-x86_64_v4.pkg.tar.zst -C ~/toolchains/gcc-16.1.1 --strip-components=1 usr

`--strip-components=1 usr` selects only the `usr/` members and drops that
component, so the result is a normal prefix with `bin/` and `lib/` at the top
and no package metadata. It comes to 300 MB.

**Nothing else is required to run it.** No `--sysroot`, no `-B`, no
`LD_LIBRARY_PATH`. GCC's driver computes its own exec prefix from `argv[0]`, so
a driver at `PREFIX/bin/gfortran` finds `PREFIX/lib/gcc/x86_64-pc-linux-gnu/16/f951`
by itself; `gfortran -print-search-dirs` shows every path pointing inside the
prefix. `as` and `ld` are deliberately NOT in the prefix and come from the
system binutils on `PATH`. Note that binutils DID move in the same transaction,
2.47-2 to 2.47-4, so "the assembler did not change" is not one of the things
this prefix holds fixed. It is the same upstream version at a new pkgrel, and
nothing has ever implicated it.

**To build with it**, override the compiler OpenMPI's wrappers call. ExoPlaSim's
`bld/compilerargs` invokes `mpif90`, `mpicc` and `mpicxx`, and the three
variables are read by the wrappers rather than baked into them:

    export OMPI_FC=$HOME/toolchains/gcc-16.1.1/bin/gfortran
    export OMPI_CC=$HOME/toolchains/gcc-16.1.1/bin/gcc
    export OMPI_CXX=$HOME/toolchains/gcc-16.1.1/bin/g++

Confirm the override took, because a wrapper that ignores it fails silently and
produces a working binary at the wrong speed: `OMPI_FC=... mpif90 --version`
must report 16.1.1, not 16.2.1.

**To check what actually built a finished binary**, read its `.comment` section
with `readelf -p .comment <binary>`. A 16.2.1 build of anything on this system
shows two strings, `16.1.1` and `16.2.1`, because glibc's `crt1.o` was compiled
by the older gcc and has not been rebuilt since; a genuine 16.1.1 build shows
only `16.1.1`. The presence of a `16.2.1` string is the tell, not the presence
of a `16.1.1` one.

**The runtime library needs no special handling, and this was tested rather
than assumed.** A binary built here records `libgfortran.so.5` as its only
Fortran `NEEDED` entry, carries no `RPATH` or `RUNPATH`, and therefore resolves
at run time to the system `/usr/lib/libgfortran.so.5`, which is now 16.2.1. That
is fine. Checked 2026-08-18: the two copies define the identical set of eleven
`GFORTRAN_*` version nodes and export 1,760 symbols each, with nothing present
in the 16.1.1 copy that is missing from the 16.2.1 one, and a test binary needs
only `GFORTRAN_8`. Link time already prefers the private copy, since the
prefix's own `lib/` precedes `/usr/lib` in the library search path.

`-static-libgfortran` is NOT available and is not needed. No `libgfortran.a` is
packaged anywhere on this distribution, in the prefix or in `/usr`, so the flag
fails at link with "cannot find -lgfortran". If a future change ever does make
the shared runtime a problem, the working lever is
`LD_LIBRARY_PATH=$HOME/toolchains/gcc-16.1.1/lib`, which was verified to move
both `libgfortran.so.5` and `libgcc_s.so.1` onto the private copies; an `-Wl,-rpath`
at link time would make that permanent. Neither is needed today.

C++ links against the system `libstdc++.so.6`, since the shared C++ runtime is
in a separate `libstdc++` package that is not extracted here. The same soname
argument applies and a C++ test program compiles and runs.

**Verified 2026-08-18** on this prefix: `gfortran --version` reports 16.1.1,
a Fortran hello-world compiles and runs, `mpif90` under `OMPI_FC` reports 16.1.1,
and an `MPI_Allreduce` program built through `mpif90` with ExoPlaSim's exact flag
set runs correctly on 4 ranks. Note that OpenMPI's `mpi.mod` is built by the
system 16.2.1 compiler and is read without complaint, because both compilers
write GFORTRAN module format version 16; that compatibility is the one thing
here most likely to break on a future GCC major bump, and it fails loudly at
compile time if it does.

**This prefix is not managed by pacman.** It receives no security updates, it is
invisible to `pacman -Qo` and to every dependency check, and nothing will ever
tell you it is out of date. It is a BUILD TOOL, not a general
compiler: do not put it on `PATH`, and do not use it for anything but a
deliberate compiler comparison. There is no regression to wait out, so nothing
here is pending; delete the directory whenever the disk is wanted. Re-creating it needs nothing but the
package files, so if the pacman cache is ever cleared with `paccache` or
`pacman -Sc`, copy those four files somewhere durable first.

## Where it goes once fetched

Bulk under `<component>/data/reference/` or `references/`, excluded by
`.gitignore`, with a README recording the exact URL, the DOI of the paper behind
it and which fields are used. The extract that an analysis actually reads is
tracked, because it is kilobytes and it makes the result reproducible without
the bulk. `hydrography/data/reference/README.md` is the worked example.
