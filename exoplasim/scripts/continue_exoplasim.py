#!/usr/bin/env python3
"""Resume a prepared ExoPlaSim run from its latest restart file."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

import exoplasim as exo
from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import CONFIG, INPUTS, RUNS  # noqa: E402
from run_exoplasim import (  # noqa: E402
    LANDMAP,
    TOPOMAP,
    surface_field_report,
    REGULAR_CODES,
    SNAPSHOT_CODES,
    derive,
    file_sha256,
    run_id,
)


REQUIRED = {
    "ts", "tas", "pr", "ua", "va", "snd", "sic", "sit", "mrso",
    "mrro", "evap", "ntr", "hfns", "lsm", "sg", "nu",
}


def output_years(run_dir: Path) -> list[int]:
    years = []
    for path in run_dir.glob("MOST.*.nc"):
        match = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if match:
            years.append(int(match.group(1)))
    return sorted(years)


def validate_year(path: Path, expected_times: int | None = None) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing model output {path}")
    with Dataset(path) as nc:
        if expected_times is not None and len(nc.dimensions["time"]) != expected_times:
            raise RuntimeError(
                f"{path} has {len(nc.dimensions['time'])} time samples; "
                f"expected {expected_times}"
            )
        missing = sorted(REQUIRED - set(nc.variables))
        if missing:
            raise RuntimeError(f"{path} is missing required variables: {missing}")
        bad = [name for name in sorted(REQUIRED) if not np.isfinite(nc[name][:]).all()]
        if bad:
            raise RuntimeError(f"{path} has non-finite variables: {bad}")


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def year_diagnostics(path: Path) -> dict:
    with Dataset(path) as nc:
        weights = leggauss(len(nc.dimensions["lat"]))[1][::-1]
        values = {}
        for name in ["ts", "ntr", "hfns", "sic", "pr"]:
            field = np.asarray(nc[name][:], dtype=float)
            mean = global_mean(field, weights).mean()
            if name == "pr":
                mean *= 86400.0 * 1000.0
            values[name] = float(mean)
    diag_path = path.with_name(path.name.replace("MOST.", "MOST_DIAG.").replace(".nc", ""))
    runtime = None
    if diag_path.is_file():
        match = re.search(
            r"Seconds per sim year:\s*([0-9.]+)",
            diag_path.read_text(encoding="ascii", errors="replace"),
        )
        if match:
            runtime = float(match.group(1))
    return {
        "year_index": int(path.name.split(".")[1]),
        "surface_temperature_k": values["ts"],
        "toa_net_radiation_w_m2": values["ntr"],
        "surface_downward_heat_flux_w_m2": values["hfns"],
        "planetary_sea_ice_fraction": values["sic"],
        "precipitation_mm_day": values["pr"],
        "native_runtime_seconds": runtime,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--flux-ratio", type=float, default=None)
    parser.add_argument("--orbits", type=int, default=5)
    parser.add_argument(
        "--seasonal-output", action="store_true",
        help="Also retain instantaneous seasonal snapshots for this segment",
    )
    args = parser.parse_args()
    if args.orbits < 1:
        raise ValueError("--orbits must be positive")

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    flux_ratio = float(
        config["orbit"]["baseline_flux_earth"]
        if args.flux_ratio is None else args.flux_ratio
    )
    derived = derive(config, flux_ratio)
    identifier = run_id(config, flux_ratio)
    run_dir = (RUNS / identifier).resolve()
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"No prepared run manifest at {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["config_sha256"] != file_sha256(config_path):
        raise RuntimeError("Current configuration differs from the run manifest; refusing to resume")

    years = output_years(run_dir)
    if not years or years != list(range(years[-1] + 1)):
        raise RuntimeError(f"Run outputs are absent or non-contiguous: {years}")
    start_year = years[-1] + 1
    restart = run_dir / f"MOST_REST.{years[-1]:05d}"
    if not restart.is_file():
        raise RuntimeError(f"Missing restart file {restart}")

    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    landmap = LANDMAP.resolve()
    topomap = TOPOMAP.resolve()
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
        inityear=start_year,
        workdir=str(run_dir),
        modelname=identifier,
        outputtype=model_cfg["output_type"],
        hyperthreading=False,
    )
    model.configure(
        flux=derived["stellar_flux_w_m2"],
        startemp=float(star["effective_temperature_k"]),
        starradius=derived["stellar_radius_solar"],
        pN2=float(atmosphere["pN2_bar"]),
        pO2=float(atmosphere["pO2_bar"]),
        pAr=float(atmosphere["pAr_bar"]),
        pCO2=float(atmosphere["pCO2_bar"]),
        rotationperiod=derived["rotation_days"],
        year=derived["orbital_year_earth_days"],
        gravity=derived["gravity_m_s2"],
        radius=float(planet["radius_earth"]),
        eccentricity=float(planet["eccentricity"]),
        obliquity=float(planet["obliquity_degrees"]),
        lonvernaleq=float(config["orbit"]["longitude_vernal_equinox_degrees"]),
        fixedorbit=True,
        keplerian=True,
        meananomaly0=0.0,
        seaice=bool(surface["sea_ice"]),
        ozone=bool(atmosphere["ozone"]),
        mldepth=float(surface["mixed_layer_depth_m"]),
        twobandalbedo=bool(config["radiation"]["two_band_albedo"]),
        timestep=float(model_cfg["timestep_minutes"]),
        physicsfilter=model_cfg["physics_filter"],
        landmap=str(landmap),
        topomap=str(topomap),
        restartfile=str(restart),
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=(int(derived["snapshot_interval_steps"]) if args.seasonal_output else 0),
        otherargs={
            "N_DAYS_PER_YEAR@plasim_namelist": str(
                derived["rotations_per_orbit_namelist"]
            )
        },
    )
    model._add_postcodes("example.nl", REGULAR_CODES)
    model.cfgpostprocessor(
        ftype="regular",
        extension=model_cfg["output_type"],
        variables=[str(code) for code in REGULAR_CODES],
        mode="grid",
        times=int(model_cfg["regular_output_bins_per_orbit"]),
        timeaverage=True,
        interpolatetimes=False,
    )
    if args.seasonal_output:
        model._add_postcodes("snapshot.nl", SNAPSHOT_CODES)
        model.cfgpostprocessor(
            ftype="snapshot",
            extension=model_cfg["output_type"],
            variables=[str(code) for code in SNAPSHOT_CODES],
            mode="grid",
            times=None,
            timeaverage=False,
            interpolatetimes=False,
        )

    surface_field_report(run_dir, config)
    started = datetime.now(timezone.utc).isoformat()
    try:
        model.run(years=args.orbits, crashifbroken=True, clean=True)
        new_diagnostics = []
        for year in range(start_year, start_year + args.orbits):
            output = run_dir / f"MOST.{year:05d}.nc"
            validate_year(
                output,
                expected_times=int(model_cfg["regular_output_bins_per_orbit"]),
            )
            if args.seasonal_output:
                validate_year(
                    run_dir / "snapshots" / f"MOST_SNAP.{year:05d}.nc",
                    expected_times=int(model_cfg["seasonal_samples_per_orbit"]),
                )
            new_diagnostics.append(year_diagnostics(output))
    except Exception:
        manifest["status"] = "failed_during_resume"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        raise

    is_post_equilibrium_climatology = bool(
        args.seasonal_output
        and "equilibrium_cutoff_year_index" in manifest
        and start_year > int(manifest["equilibrium_cutoff_year_index"])
    )
    manifest["status"] = (
        "climatology_complete"
        if is_post_equilibrium_climatology
        else "spinup_in_progress"
    )
    manifest["completed_orbits"] = start_year + args.orbits
    manifest.setdefault("segments", []).append(
        {
            "start_year_index": start_year,
            "end_year_index": start_year + args.orbits - 1,
            "seasonal_output": args.seasonal_output,
            "purpose": (
                "post_equilibrium_climatology"
                if is_post_equilibrium_climatology
                else "spinup"
            ),
            "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "input_restart_sha256": file_sha256(restart),
            "diagnostics": new_diagnostics,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(new_diagnostics, indent=2))


if __name__ == "__main__":
    main()
