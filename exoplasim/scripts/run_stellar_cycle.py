#!/usr/bin/env python3
"""Run or resume smooth sinusoidal stellar-flux experiments with ExoPlaSim."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import struct
import sys

import exoplasim as exo
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import CONFIG, INPUTS, PATCHES, RUNS  # noqa: E402
from run_exoplasim import (  # noqa: E402
    surface_sra,
    geography_tag,
    spectrum_tag,
    stage_surface_extras,
    surface_field_report,
    REGULAR_CODES,
    derive,
    file_sha256,
)
from continue_exoplasim import validate_year, year_diagnostics  # noqa: E402


def load_cases(config: dict) -> dict:
    """Cycle amplitudes, from config rather than hardcoded.

    They were literals here while the baseline was 0.90. The baseline is now
    0.96 and the amplitudes moved with it, so leaving them in the source would
    have silently applied the old range about the new centre.
    """
    cases = config.get("stellar_cycle", {}).get("cases")
    if not cases:
        raise RuntimeError("config/planet.yaml has no stellar_cycle.cases")
    return cases
JULIAN_YEAR_DAYS = 365.25


def restart_integer(path: Path, wanted: str) -> int:
    """Read an integer from ExoPlaSim's gfortran sequential restart format."""
    with path.open("rb") as handle:
        while marker := handle.read(4):
            name_size = struct.unpack("<i", marker)[0]
            name = handle.read(name_size).decode("ascii", errors="replace").strip()
            if struct.unpack("<i", handle.read(4))[0] != name_size:
                raise ValueError(f"Invalid record marker in {path}")
            value_size = struct.unpack("<i", handle.read(4))[0]
            value = handle.read(value_size)
            if struct.unpack("<i", handle.read(4))[0] != value_size:
                raise ValueError(f"Invalid value record marker in {path}")
            if name == wanted:
                if value_size != 4:
                    raise ValueError(f"{wanted} is not a 4-byte integer in {path}")
                return struct.unpack("<i", value)[0]
    raise KeyError(f"{wanted} was not found in {path}")


def completed_years(run_dir: Path) -> list[int]:
    years = []
    for path in run_dir.glob("MOST.*.nc"):
        match = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if match:
            years.append(int(match.group(1)))
    return sorted(years)


def cycle_run_id(config: dict, case: str, period_earth_years: float) -> str:
    cases = load_cases(config)
    p = config["planet"]
    a = config["atmosphere"]
    m = config["model"]
    lo = round(100 * cases[case]["minimum_flux_earth"])
    hi = round(100 * cases[case]["maximum_flux_earth"])
    identifier = (
        f"{str(m['resolution']).lower()}l{int(m['layers'])}p{int(m['ncpus'])}"
        f"_cycle{period_earth_years:g}ey_s{lo:03d}-{hi:03d}"
        f"_co2{round(1e6 * float(a['pCO2_bar'])):04d}ppm"
        f"_rot{float(p['rotation_hours']):g}h"
        f"_obl{float(p['obliquity_degrees']):g}"
        f"_e{round(1000 * float(p['eccentricity'])):03d}"
    )
    identifier += spectrum_tag(config) + f"_g{geography_tag(config)}"
    return identifier.replace(".", "p")


def make_model(
    config: dict,
    derived: dict,
    identifier: str,
    run_dir: Path,
    source_dir: Path,
    restart: Path,
    output_start_year: int,
    cycle: dict,
) -> exo.Model:
    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    landmap = surface_sra(config, 172).resolve()
    topomap = surface_sra(config, 129).resolve()
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
        inityear=output_start_year,
        workdir=str(run_dir),
        source=str(source_dir),
        modelname=identifier,
        outputtype=model_cfg["output_type"],
        hyperthreading=False,
    )
    model.configure(
        flux=cycle["mean_flux_w_m2"],
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
        glaciers={
            "toggle": bool(surface.get("glaciers", {}).get("enabled", False)),
            "mindepth": float(surface.get("glaciers", {}).get("min_snow_depth_m", 2.0)),
            "initialh": float(surface.get("glaciers", {}).get("initial_height_m", -1.0)),
        },
        ozone=bool(atmosphere["ozone"]),
        mldepth=float(surface["mixed_layer_depth_m"]),
        twobandalbedo=bool(config["radiation"]["two_band_albedo"]),
        timestep=float(model_cfg["timestep_minutes"]),
        physicsfilter=model_cfg["physics_filter"],
        landmap=str(landmap),
        topomap=str(topomap),
        restartfile=str(restart),
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=0,
        otherargs={
            "N_DAYS_PER_YEAR@plasim_namelist": str(
                derived["rotations_per_orbit_namelist"]
            ),
            # 6068 is divisible by 164, preventing an unflushed partial output
            # accumulator at the end of each modeled orbit.
            "NSTPW@plasim_namelist": "164",
            "NSOLCYCLE@radmod_namelist": "1",
            "GSOLSTART@radmod_namelist": str(cycle["start_model_step"]),
            "GSOLAMP@radmod_namelist": repr(cycle["semi_amplitude_w_m2"]),
            "GSOLPERIOD@radmod_namelist": repr(cycle["period_model_steps"]),
            "GSOLPHASE@radmod_namelist": "0.0",
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
    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--baseline", type=str, default=None,
                        help="run id or directory to centre the cycle on")
    parser.add_argument("--restart", type=str, default=None,
                        help="restart filename within the baseline run")
    parser.add_argument("--period-earth-years", type=float, default=8.0)
    parser.add_argument(
        "--cycles", type=float, default=4.0,
        help="Run at least this many complete stellar cycles in total",
    )
    parser.add_argument(
        "--orbits", type=int,
        help="Override the target output-orbit count (useful for validation)",
    )
    args = parser.parse_args()
    if args.period_earth_years <= 0 or args.cycles <= 0:
        raise ValueError("Cycle period and number of cycles must be positive")

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    cases = load_cases(config)
    if args.case not in cases:
        raise SystemExit(f"unknown case {args.case!r}; config defines {sorted(cases)}")
    mean_ratio = 0.5 * (
        cases[args.case]["minimum_flux_earth"]
        + cases[args.case]["maximum_flux_earth"]
    )
    if not math.isclose(mean_ratio, float(config["orbit"]["baseline_flux_earth"])):
        raise ValueError("Cycle mean must equal the configured baseline flux")
    derived = derive(config, mean_ratio)
    identifier = cycle_run_id(config, args.case, args.period_earth_years)
    run_dir = (RUNS / identifier).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    source_dir = (INPUTS / "exoplasim_cycle_t42").resolve()
    executable = source_dir / "most_plasim_t42_l10_p8.x"
    patch_file = (PATCHES / "exoplasim-3.4.2-star-cycle.patch").resolve()
    if not executable.is_file():
        raise RuntimeError("Build the patched cycle executable first")

    # The baseline a cycle varies about must be named, not hardcoded and not
    # chosen by sort order. A cycle is variance around a mean, so starting it
    # from whichever restart sorts last -- of whichever run someone once wired in
    # -- centres the whole experiment on an unknown state.
    if args.baseline is None:
        raise SystemExit(
            "--baseline is required: the run id or directory whose equilibrated "
            "restart this cycle starts from. "
            "`python exoplasim/scripts/index_runs.py` lists what is available.")
    cand = Path(args.baseline)
    baseline_dir = (cand if cand.is_dir() else RUNS / args.baseline).resolve()
    if not (baseline_dir / "run_manifest.json").is_file():
        raise SystemExit(f"{baseline_dir} has no run_manifest.json")
    if args.restart is None:
        raise SystemExit(
            "--restart is required: the restart file to begin from, e.g. "
            "MOST_REST.00059. Naming it is the point -- a cycle centred on the "
            "wrong state is not detectable from its own output.")
    initial_restart = baseline_dir / args.restart
    if not initial_restart.is_file():
        raise SystemExit(f"no restart at {initial_restart}")
    start_step = restart_integer(initial_restart, "nstep")
    timestep_seconds = float(config["model"]["timestep_minutes"]) * 60.0
    period_steps = args.period_earth_years * JULIAN_YEAR_DAYS * 86400.0 / timestep_seconds
    earth_constant = float(config["orbit"]["earth_solar_constant_w_m2"])
    semi_amplitude_ratio = 0.5 * (
        cases[args.case]["maximum_flux_earth"]
        - cases[args.case]["minimum_flux_earth"]
    )
    cycle = {
        "case": args.case,
        "minimum_flux_earth": cases[args.case]["minimum_flux_earth"],
        "maximum_flux_earth": cases[args.case]["maximum_flux_earth"],
        "mean_flux_earth": mean_ratio,
        "mean_flux_w_m2": mean_ratio * earth_constant,
        "semi_amplitude_w_m2": semi_amplitude_ratio * earth_constant,
        "period_earth_years": args.period_earth_years,
        "period_local_orbits": period_steps / derived["runsteps_per_orbit"],
        "period_model_steps": period_steps,
        "start_model_step": start_step,
        "phase_convention": "sinusoid begins at mean flux and rises",
        "stellar_spectrum_treatment": "fixed 4965 K blackbody spectral partition",
    }
    target_orbits = (
        args.orbits
        if args.orbits is not None
        else math.ceil(args.cycles * period_steps / derived["runsteps_per_orbit"])
    )
    if target_orbits < 1:
        raise ValueError("Target orbit count must be positive")

    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        immutable = {
            "config_sha256": file_sha256(config_path),
            "patched_executable_sha256": file_sha256(executable),
            "cycle": cycle,
        }
        for key, expected in immutable.items():
            if manifest.get(key) != expected:
                raise RuntimeError(f"Existing manifest has incompatible {key}")
    else:
        manifest = {
            "schema_version": 1,
            "run_id": identifier,
            "status": "prepared",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "config_path": str(config_path),
            "config_sha256": file_sha256(config_path),
            "source_config": config,
            "fixed_orbit_parameters": derived,
            "cycle": cycle,
            "initial_restart": str(initial_restart),
            "initial_restart_sha256": file_sha256(initial_restart),
            "patch_path": str(patch_file),
            "patch_sha256": file_sha256(patch_file),
            "patched_executable": str(executable),
            "patched_executable_sha256": file_sha256(executable),
            "target_orbits": target_orbits,
            "diagnostics": [],
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    years = completed_years(run_dir)
    if years and years != list(range(years[-1] + 1)):
        raise RuntimeError(f"Outputs are non-contiguous: {years}")
    output_start = years[-1] + 1 if years else 0
    if output_start >= target_orbits:
        print(json.dumps({"status": "already_complete", "orbits": output_start}))
        return
    restart = run_dir / f"MOST_REST.{years[-1]:05d}" if years else initial_restart
    if not restart.is_file():
        raise RuntimeError(f"Restart file is missing: {restart}")

    model = make_model(
        config, derived, identifier, run_dir, source_dir, restart,
        output_start, cycle,
    )
    manifest["status"] = "running"
    manifest["target_orbits"] = target_orbits
    manifest["last_started_utc"] = datetime.now(timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    stage_surface_extras(run_dir, config)
    surface_field_report(run_dir, config)
    try:
        for year in range(output_start, target_orbits):
            model.run(years=1, crashifbroken=True, clean=True)
            output = run_dir / f"MOST.{year:05d}.nc"
            validate_year(
                output,
                expected_times=int(config["model"]["regular_output_bins_per_orbit"]),
            )
            diagnostic = year_diagnostics(output)
            diagnostic["cycle_start_phase"] = float(
                ((year * derived["runsteps_per_orbit"]) / period_steps) % 1.0
            )
            diagnostic["cycle_end_phase"] = float(
                (((year + 1) * derived["runsteps_per_orbit"]) / period_steps) % 1.0
            )
            manifest["diagnostics"] = [
                d for d in manifest["diagnostics"] if d["year_index"] != year
            ] + [diagnostic]
            manifest["completed_orbits"] = year + 1
            manifest["completed_cycles"] = (
                (year + 1) * derived["runsteps_per_orbit"] / period_steps
            )
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
    except BaseException:
        manifest["status"] = "interrupted_or_failed"
        manifest["last_stopped_utc"] = datetime.now(timezone.utc).isoformat()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        raise

    manifest["status"] = "simulation_complete"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"run_dir": str(run_dir), "completed_orbits": target_orbits}, indent=2))


if __name__ == "__main__":
    main()
