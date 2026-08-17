# Earth reference data for the lake solver validation

What HYD-4 is tested against. The bulk is excluded from git and the extract is
not: `endorheic_lakes.json` beside this file is the table
`validate_lake_solver.py` actually reads, and it is enough to re-run the test
without re-downloading six gigabytes.

Same reasoning as `pedology/data/reference/`, and as `references/INDEX.md`:
exclude the bulk, keep the record.

## What is here, and where it came from

| what | source | size |
| --- | --- | --- |
| `HydroLAKES_polys_v10_shp/` | HydroLAKES v1.0, Messager, Lehner, Grill, Nedeva, Schmitt (2016), *Estimating the volume and age of water stored in global lakes using a geo-statistical approach*, Nat. Commun. 7, 13603. `10.1038/ncomms13603`. Downloaded from `https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip` | 2.2 GB unpacked |
| `hybas_<region>_lev05_v1c.*` | HydroBASINS v1.c level 5, Lehner, Grill (2013), *Global river hydrography and network routing*, Hydrol. Process. 27(15), 2171-2186. `10.1002/hyp.9740`. Downloaded per region from `https://data.hydrosheds.org/file/hydrobasins/standard/hybas_<region>_lev05_v1c.zip`, regions af ar as au eu gr na sa si | 60 MB |
| `glev_evaporation_rate.csv` | GLEV, Zhao, Li, Zhang, Gao (2022), *Evaporative water loss of 1.42 million global lakes*, Nat. Commun. 13, 3686. `10.1038/s41467-022-31125-6`. Dataset at `https://doi.org/10.5281/zenodo.4646621`, file `0_evaporation_rate.csv` | 4.0 GB |

The HydroSHEDS host refuses `HEAD` and returns 403 without a browser
`User-Agent`; a plain `curl -L -A "Mozilla/5.0 ..."` works. Zenodo needs neither.

## Which fields the test uses, and why those

From HydroLAKES: `Lake_area`, `Wshd_area` and `Dis_avg`, which give the catchment
runoff depth `R = Dis_avg / (Wshd_area - Lake_area)` under one set of
definitions. That independence from the lake is the entire test, and it is what
the endorheic-lake literature could not supply -- see
`hydrography/notes/lake-solver-validation.md`.

`Dis_avg` is modelled by WaterGAP rather than gauged. It is independent of the
equation under test, which is what the test requires, but it is not an
observation, and the note says what that costs.

From HydroBASINS: `ENDO` and `MAIN_BAS`, to find each endorheic system and its
most downstream lake. From GLEV: monthly evaporation rates keyed by `Hylak_id`,
averaged to an annual rate.

Lake precipitation is not here. It comes from ERA5 through the Open-Meteo archive
at each lake's pour point, cached under `hydrography/analysis/cache_era5/`, which
is excluded for the same reason.

## exorheic_impounded_lakes.json

Earth's sill-impounded through-flowing basins: a lake that persists although a
river crosses its sill, which is the terrestrial analogue of a MARGINAL carve
verdict. Extracted from HydroLAKES joined to HydroBASINS by pour point, with the
selection declared in the file and the argument in
`hydrography/notes/retain-fraction.md`. It is what falsified `Q_full = 1 m3/s`,
which gives every one of these basins a retain of 0.
