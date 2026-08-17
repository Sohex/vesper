#!/usr/bin/env python3
"""Average a contiguous post-equilibrium ExoPlaSim segment into climatologies.

Three products, from one pass over the same orbits.

**The averaged climatology**, as before, unchanged. One 12-bin year and one
32-snapshot year, averaged across the segment. This is what the described
climatology and every downstream product has always used.

**Per-orbit climatologies**, one file per model year, written with `--per-year`.
The average is the right forcing for a fixed climate and the wrong one for a
variable star: `biosphere/` cycles through however many years it is given, and a
single repeating year cannot represent a stellar cycle at all. These are what
feed it.

**A climate series**, always written and small enough to track. One number per
time bin per orbit, globally and over land, for each quantity worth watching
over time rather than in a map, with annual means derived from the bins. It exists because three separate questions in this project turned out
to need it and none of them could be answered from an average:

  - whether a drift is still running, which the average hides by construction
  - the amplitude and phase lag of the climate against a stellar cycle, which is
    the whole point of running one
  - whether the budget residuals are constant in time. The 0.446 W/m2 that does
    not close between the top of the atmosphere and the surface was shown to be
    structural by comparing three separate runs; within one run it was invisible.

It also carries the coldest and warmest month per orbit, which no average can
reconstruct and which is what PFT survival thresholds actually see: 13.8% of
this world's land sits within one stellar cycle's swing of such a threshold.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS

from orbit import EARTH_CALENDAR_YEAR_DAYS


COORDINATES = {"time", "lat", "lon", "lev", "levp", "fourier", "modes"}
ORBITAL_GEOMETRY = {"nu", "lambda", "zdec", "rdist", "rasc"}


def identity_from_run(run_dir: Path) -> dict:
    """What world a climatology describes, read from the run that produced it.

    A climatology is the most widely shared artifact in this project --
    pedology, hydrography and the biosphere all read one -- and until now it
    carried orbit counts and source filenames but nothing that said which
    terrain, flux or spectrum it was. So a consumer could not tell a current
    climatology from a superseded one, and none of them checked. That is the
    root of a whole class of failure here: a stale climatology produces a
    plausible number from the wrong world rather than an error.
    """
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return {}
    m = json.loads(manifest_path.read_text(encoding="utf-8"))
    phys = m.get("physical") or {}
    out = {
        "vesper_run_id": m.get("run_id"),
        "vesper_source_build": m.get("source_build"),
        "vesper_config_sha256": m.get("config_sha256"),
        "vesper_geography": phys.get("geography"),
        "vesper_flux_ratio": phys.get("flux_ratio"),
        "vesper_stellar_spectrum": phys.get("stellar_spectrum"),
    }
    return {k: v for k, v in out.items() if v is not None}


def average_files(paths: list[Path], output: Path, product: str,
                  identity: dict | None = None) -> None:
    """Write a compressed, same-grid mean across corresponding model orbits."""
    if not paths:
        raise ValueError("No input files supplied")
    for path in paths:
        if not path.is_file():
            raise RuntimeError(f"Missing climatology input {path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(paths[0]) as template, Dataset(output, "w", format="NETCDF4") as dst:
        for name, dim in template.dimensions.items():
            dst.createDimension(name, len(dim))

        dst.setncatts({name: template.getncattr(name) for name in template.ncattrs()})
        dst.setncattr("climatology_product", product)
        dst.setncattr("climatology_orbit_count", len(paths))
        dst.setncattr("climatology_start_year_index", int(paths[0].name.split(".")[-2]))
        dst.setncattr("climatology_end_year_index", int(paths[-1].name.split(".")[-2]))
        dst.setncattr("climatology_created_utc", datetime.now(timezone.utc).isoformat())
        dst.setncattr("climatology_source_files", ",".join(path.name for path in paths))
        # Which world this describes. Without it a climatology is anonymous and
        # every consumer has to be told out of band which build it belongs to.
        for key, value in (identity or {}).items():
            dst.setncattr(key, value)

        for name, source in template.variables.items():
            fill = source.getncattr("_FillValue") if "_FillValue" in source.ncattrs() else None
            kwargs = {"zlib": True, "complevel": 4, "shuffle": True}
            if fill is not None:
                kwargs["fill_value"] = fill
            target = dst.createVariable(name, source.dtype, source.dimensions, **kwargs)
            target.setncatts(
                {
                    attr: source.getncattr(attr)
                    for attr in source.ncattrs()
                    if attr != "_FillValue"
                }
            )

            first = np.ma.asarray(source[:])
            should_average = "time" in source.dimensions and name not in (
                COORDINATES | ORBITAL_GEOMETRY
            )
            if should_average:
                total = np.ma.asarray(first, dtype=np.float64)
                for path in paths[1:]:
                    with Dataset(path) as nc:
                        if name not in nc.variables or nc[name].shape != source.shape:
                            raise RuntimeError(f"Variable mismatch for {name} in {path}")
                        total += np.ma.asarray(nc[name][:], dtype=np.float64)
                target[:] = total / len(paths)
            else:
                for path in paths[1:]:
                    with Dataset(path) as nc:
                        if name not in nc.variables or nc[name].shape != source.shape:
                            raise RuntimeError(f"Variable mismatch for {name} in {path}")
                        if not np.ma.allclose(first, nc[name][:], rtol=1e-6, atol=1e-7):
                            raise RuntimeError(
                                f"Coordinate/orbital variable {name} differs in {path}"
                            )
                target[:] = first


# Quantities worth one number per orbit. Fluxes stay in their native units and
# are area-weighted; nothing here is annualised, because the point is to compare
# orbits with each other rather than with Earth.
SERIES_FIELDS = {
    "ts": "surface temperature, K",
    "tas": "2 m air temperature, K",
    "pr": "precipitation, m/s",
    "evap": "evaporation, m/s (negative upward)",
    "mrro": "river-routed net water flux, m/s (NOT local runoff)",
    "rst": "TOA net shortwave, W/m2",
    "rlut": "TOA net longwave, W/m2",
    "ntr": "TOA net radiation, W/m2",
    "rss": "surface net shortwave, W/m2",
    "rls": "surface net longwave, W/m2",
    "hfss": "surface sensible heat flux, W/m2",
    "hfls": "surface latent heat flux, W/m2",
    "hfns": "surface net heat flux, W/m2",
    "clt": "cloud area fraction",
    "alb": "surface albedo",
    "sic": "sea ice cover",
    "snd": "snow thickness, m",
}


def climate_series(paths: list[Path], years: list[int]) -> dict:
    """Per orbit AND per time bin, globally and over land, plus the residuals.

    Seasonal rather than annual, because the two questions a stellar cycle raises
    are how far the climate swings and *which season* swings. Those are different
    questions and an annual mean answers only the first. A cycle that deepens
    winter without touching summer, or the reverse, is a different world for
    anything with a survival threshold on the coldest month, and this world has
    13.8% of its land within one cycle's swing of such a threshold.

    Annual means are derived from the bins here rather than accumulated
    separately, so the two cannot drift apart. Bins are equal-length time
    averages, so a plain mean over them is the annual mean.

    Deliberately small even so: about 17 fields by two masks by twelve bins per
    orbit, rounded, which is a few hundred KB for a full stellar cycle and
    tracked rather than regenerated.
    """
    seasonal: dict[str, dict[str, list]] = {"global": {}, "land": {}}
    extremes: dict[str, list] = {"land_coldest_month_c": [],
                                 "land_warmest_month_c": []}
    closure: dict[str, list] = {"toa_minus_surface_w_m2": [],
                                "land_p_minus_e_minus_mrro_mm_per_orbit": []}

    def r6(value):
        return float(f"{value:.6g}")

    for path in paths:
        with Dataset(path) as data:
            lat = np.asarray(data["lat"][:], dtype=float)
            nlon = len(data["lon"][:])
            land = np.asarray(data["lsm"][0], dtype=float) > 0.5
            weight = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
            weight = weight / weight.sum()
            land_weight_total = float(weight[land].sum())

            def per_bin(name: str, mask=None) -> list[float]:
                """Area-weighted mean of each time bin, so the season survives."""
                field = np.asarray(data[name][:], dtype=float)
                if mask is None:
                    return [r6(float((bin_field * weight).sum()))
                            for bin_field in field]
                if not land_weight_total:
                    return [0.0] * field.shape[0]
                return [r6(float((bin_field[mask] * weight[mask]).sum()
                                 / land_weight_total)) for bin_field in field]

            for name in SERIES_FIELDS:
                if name not in data.variables:
                    continue
                seasonal["global"].setdefault(name, []).append(per_bin(name))
                seasonal["land"].setdefault(name, []).append(per_bin(name, land))

            tas = np.asarray(data["tas"][:], dtype=float) - 273.15
            lw = weight[land]
            extremes["land_coldest_month_c"].append(
                r6(float(np.average(tas.min(axis=0)[land], weights=lw))))
            extremes["land_warmest_month_c"].append(
                r6(float(np.average(tas.max(axis=0)[land], weights=lw))))

            # Residuals per bin as well as per orbit. If the energy gap is
            # structural it should be flat across seasons too, which is a
            # stronger statement than flat across orbits and costs nothing.
            ntr = per_bin("ntr")
            hfns = per_bin("hfns")
            closure["toa_minus_surface_w_m2"].append(
                [r6(a - b) for a, b in zip(ntr, hfns)])
            seconds = 86400.0 * EARTH_CALENDAR_YEAR_DAYS
            p = per_bin("pr", land)
            e = per_bin("evap", land)
            q = per_bin("mrro", land)
            closure["land_p_minus_e_minus_mrro_mm_per_orbit"].append(
                [r6((pi + ei - qi) * seconds * 1000.0)
                 for pi, ei, qi in zip(p, e, q)])

    def annual(block: dict) -> dict:
        return {name: [r6(sum(bins) / len(bins)) for bins in per_year]
                for name, per_year in block.items()}

    return {
        "years": years,
        "bins_per_orbit": len(next(iter(seasonal["global"].values()))[0])
        if seasonal["global"] else 0,
        "fields": SERIES_FIELDS,
        "seasonal": seasonal,
        "annual": {"global": annual(seasonal["global"]),
                   "land": annual(seasonal["land"])},
        "extremes": extremes,
        "closure": closure,
        "layout_note": (
            "seasonal[mask][field] is [orbit][bin]; annual[mask][field] is the "
            "mean over bins and is derived from it, so the two cannot disagree. "
            "closure entries are [orbit][bin]."),
        "closure_note": (
            "toa_minus_surface should be zero in steady state and is not; it has "
            "run near -0.45 W/m2 across every run measured. "
            "land_p_minus_e_minus_mrro is expected to be large and positive: "
            "mrro is river-routed net water flux, not local runoff, so the "
            "residual is the water the rivers carried away."),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--start-year", type=int, required=True)
    parser.add_argument("--end-year", type=int, required=True)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "climatology")
    parser.add_argument("--label", default="baseline")
    parser.add_argument("--per-year", action="store_true",
                        help="also write one climatology per orbit, which is "
                             "what biosphere/ needs to see a variable star")
    args = parser.parse_args()
    if args.end_year < args.start_year:
        raise ValueError("--end-year must not precede --start-year")
    if not re.fullmatch(r"[a-z0-9_-]+", args.label):
        raise ValueError("--label may contain only lowercase letters, digits, _ and -")

    run_dir = args.run_dir.resolve()
    years = range(args.start_year, args.end_year + 1)
    regular = [run_dir / f"MOST.{year:05d}.nc" for year in years]
    snapshots = [run_dir / "snapshots" / f"MOST_SNAP.{year:05d}.nc" for year in years]
    output_dir = args.output.resolve()
    regular_output = output_dir / f"{args.label}_regular_climatology.nc"
    snapshot_output = output_dir / f"{args.label}_snapshot_climatology.nc"
    identity = identity_from_run(run_dir)
    average_files(regular, regular_output,
                  "time-bin means averaged across model orbits", identity)
    average_files(snapshots, snapshot_output,
                  "instantaneous orbital snapshots averaged across model orbits",
                  identity)

    per_year_outputs = []
    if args.per_year:
        for year, path in zip(years, regular):
            target = output_dir / f"{args.label}_year{year:05d}_regular_climatology.nc"
            average_files([path], target, f"single model orbit {year}", identity)
            per_year_outputs.append(str(target))

    series = climate_series(regular, list(years))
    series_path = output_dir / f"{args.label}_climate_series.json"
    series_path.write_text(json.dumps(series, indent=2) + "\n", encoding="utf-8")

    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest.setdefault("climatologies", {})[args.label] = {
            "start_year_index": args.start_year,
            "end_year_index": args.end_year,
            "orbit_count": args.end_year - args.start_year + 1,
            "regular": str(regular_output),
            "snapshots": str(snapshot_output),
            "per_year": per_year_outputs,
            "climate_series": str(series_path),
        }
        # Retain the original convenience key for existing tooling.
        if args.label == "baseline":
            manifest["climatology"] = manifest["climatologies"][args.label]
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {"regular": str(regular_output), "snapshots": str(snapshot_output),
             "per_year": per_year_outputs, "climate_series": str(series_path)},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
