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

**Zenodo**, for datasets published with a paper. The record API,
`https://zenodo.org/api/records/<id>`, lists files and sizes without
authentication, which is worth checking before starting a multi-gigabyte
download.

## Where it goes once fetched

Bulk under `<component>/data/reference/` or `references/`, excluded by
`.gitignore`, with a README recording the exact URL, the DOI of the paper behind
it and which fields are used. The extract that an analysis actually reads is
tracked, because it is kilobytes and it makes the result reproducible without
the bulk. `hydrography/data/reference/README.md` is the worked example.
