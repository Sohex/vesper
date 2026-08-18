# Where external data comes from

A record of routes into data this project does not generate, so the next person
needing a dataset starts from a list rather than from a search engine. Written
2026-08-17 after a validation stalled on an API quota that a bulk mirror would
have made irrelevant.

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
  Teff 4900 and 5000 both exist and this star's log g is 4.57, so log g 4.5 at
  solar metallicity is the comparison point. That makes the band-1 flux fraction
  an independently checkable number rather than one this project computed alone.
- **`corrk_data/N2-0.000376CO2-H2Ovar_2026`** and its siblings up to 0.95 CO2 are
  correlated-k tables from modern line lists for an N2 atmosphere with variable
  water. This planet is N2-dominated at 450 ppm CO2, so these bracket it, and a
  correlated-k treatment is an independent route to the shortwave CO2 and water
  vapour absorption that PHYS-1 and PHYS-6 derived from Howard's 1950s band data.
- **`aerosol_properties/optprop_dustvis_n50.dat`** and `optprop_dustir_n50.dat` for
  dust optics, against the DUST-12 choice.
- Its physics tree has NO wet deposition and no deposition velocity, checked by
  grep 2026-08-18. That is useful negatively: a second independent GCM has the same
  gap ExoPlaSim does, so DUST-7's scavenging is a common omission rather than a
  peculiarity, and there is nothing to borrow.

## Where it goes once fetched

Bulk under `<component>/data/reference/` or `references/`, excluded by
`.gitignore`, with a README recording the exact URL, the DOI of the paper behind
it and which fields are used. The extract that an analysis actually reads is
tracked, because it is kilobytes and it makes the result reproducible without
the bulk. `hydrography/data/reference/README.md` is the worked example.
