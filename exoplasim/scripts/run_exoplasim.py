#!/usr/bin/env python3
"""Prepare, validate, and optionally run a reproducible ExoPlaSim experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess

import exoplasim as exo
from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, INPUTS, RUNS


EARTH_STANDARD_GRAVITY = 9.80665
EARTH_SIDEREAL_YEAR_DAYS = 365.2568983
SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0

# Climate and balance fields retained in compact annual products.
REGULAR_CODES = [
    50, 51, 52, 53, 54, 110, 129, 130, 131, 132, 133, 134, 135, 139,
    140, 141, 142, 143, 144, 146, 147, 151, 157, 160, 163, 164, 167,
    168, 170, 171, 172, 174, 175, 176, 177, 178, 179, 180,
    181, 182, 184, 203, 204, 205, 207, 208, 209, 210, 211, 218, 221,
    230, 232, 238, 259, 260, 261, 262, 263, 264, 267, 318, 320, 321,
]

# Instantaneous seasonal output omits redundant flux details and expensive
# three-dimensional humidity and vertical-motion fields.
SNAPSHOT_CODES = [
    50, 51, 52, 53, 54, 129, 130, 131, 132, 134, 139, 140, 141, 142,
    143, 144, 151, 160, 163, 164, 167, 168, 170, 171, 172,
    175, 180, 181, 182, 210, 211, 218, 230, 232, 259, 260, 261, 263,
    318, 320, 321,
]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_version(command: list[str]) -> str | None:
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return (result.stdout or result.stderr).splitlines()[0]


def derive(config: dict, flux_ratio: float) -> dict:
    planet = config["planet"]
    star = config["star"]
    orbit = config["orbit"]
    model = config["model"]
    rotation_days = float(planet["rotation_hours"]) / 24.0
    # Gravity is taken from the configuration rather than derived from mass and
    # radius: World Orogen scaled this terrain's relief as 1/g, so the value it
    # ran with is canonical. Mass is checked against it so the two cannot drift.
    gravity = float(planet["gravity_m_s2"])
    implied_mass = gravity * float(planet["radius_earth"]) ** 2 / EARTH_STANDARD_GRAVITY
    if abs(implied_mass - float(planet["mass_earth"])) > 5e-4:
        raise ValueError(
            f"planet.mass_earth {planet['mass_earth']} disagrees with "
            f"planet.gravity_m_s2 {gravity} at radius {planet['radius_earth']} "
            f"(implied mass {implied_mass:.6f})"
        )
    semimajor_axis_au = math.sqrt(float(star["luminosity_solar"]) / flux_ratio)
    orbital_year_days = (
        EARTH_SIDEREAL_YEAR_DAYS
        * math.sqrt(semimajor_axis_au**3 / float(star["mass_solar"]))
    )
    timestep_minutes = float(model["timestep_minutes"])
    runsteps = round(orbital_year_days * 1440.0 / timestep_minutes)
    rotations_per_orbit = round(orbital_year_days / rotation_days)
    stellar_radius_solar = (
        math.sqrt(float(star["luminosity_solar"]))
        / (float(star["effective_temperature_k"]) / SOLAR_EFFECTIVE_TEMPERATURE_K) ** 2
    )
    return {
        "gravity_m_s2": gravity,
        "gravity_earth": gravity / EARTH_STANDARD_GRAVITY,
        "rotation_days": rotation_days,
        "stellar_radius_solar": stellar_radius_solar,
        "stellar_flux_ratio_earth": flux_ratio,
        "stellar_flux_w_m2": flux_ratio * float(orbit["earth_solar_constant_w_m2"]),
        "semimajor_axis_au": semimajor_axis_au,
        "orbital_year_earth_days": orbital_year_days,
        "orbital_year_seconds": orbital_year_days * 86400.0,
        "rotations_per_orbit_namelist": rotations_per_orbit,
        "runsteps_per_orbit": runsteps,
        "actual_run_duration_earth_days": runsteps * timestep_minutes / 1440.0,
        "snapshot_interval_steps": max(
            1, runsteps // int(model["seasonal_samples_per_orbit"])
        ),
    }


# The two SRA files that carry the geography. Every script that locates a run
# directory has to agree on these, so they live here.
LANDMAP = INPUTS / "t42" / "orogen_T42_surf_0172.sra"
TOPOMAP = INPUTS / "t42" / "orogen_T42_surf_0129.sra"


def surface_input_paths(config: dict) -> list[Path]:
    """Every SRA file that defines this run's surface, in a stable order."""
    return [INPUTS / "t42" / f"orogen_T42_surf_{code:04d}.sra"
            for code in sorted(intended_surface_codes(config))]


def geography_tag(config: dict) -> str:
    """Short digest of every surface input, for the run directory name.

    Without this, two runs sharing a config land in the same directory.
    ExoPlaSim's finalize() then selects output as the last match of
    sorted(glob("MOST*")), so a shorter new run in a directory holding a longer
    old one copies out the previous run's final year under the new name,
    silently. The prepare guard refuses that case, but naming the directory after
    the inputs stops it arising at all.

    It covers albedo as well as mask and topography, because the planned
    experiment is exactly two runs differing only in albedo. Digesting the
    geography alone would put both endmembers in one directory.
    """
    digest = hashlib.sha256()
    for path in surface_input_paths(config):
        if not path.is_file():
            raise RuntimeError(
                f"{path} is missing; run build_boundary_conditions.py and "
                "build_surface_albedo.py before preparing a run"
            )
        digest.update(file_sha256(path).encode("ascii"))
    return digest.hexdigest()[:8]


# Surface fields we supply ourselves. Everything else falls back to a uniform
# namelist default, which for a world that is not Earth is mostly the right
# answer: Earth's roughness and vegetation maps are tied to Earth's continents.
# Albedo is the exception, because the substrate genuinely varies and we know how
# from the lithology, so it is supplied when land_albedo_source is not uniform.
BASE_SURFACE_CODES = {129, 172}
ALBEDO_SURFACE_CODES = {174, 175, 176}


def intended_surface_codes(config: dict) -> set[int]:
    if str(config["model"].get("land_albedo_source", "uniform")) == "uniform":
        return set(BASE_SURFACE_CODES)
    return BASE_SURFACE_CODES | ALBEDO_SURFACE_CODES


def stage_surface_extras(run_dir: Path, config: dict) -> list[int]:
    """Copy the surface fields we generate into the run directory.

    Has to run after `configure()`, which does `rm workdir/*.sra` whenever a
    landmap or topomap is given and then writes back only those two
    (`__init__.py:2938`). Anything staged before that call is silently deleted.
    """
    staged = []
    for code in sorted(intended_surface_codes(config) - BASE_SURFACE_CODES):
        src = INPUTS / "t42" / f"orogen_T42_surf_{code:04d}.sra"
        if not src.is_file():
            raise RuntimeError(
                f"{src} is missing. Run build_surface_albedo.py, or set "
                "model.land_albedo_source: uniform to accept ExoPlaSim's 0.22."
            )
        shutil.copyfile(src, run_dir / f"N{int(config['model']['latitudes']):03d}"
                                       f"_surf_{code:04d}.sra")
        staged.append(code)
    return staged

# The uniform values those fallbacks take, from plasim/src/landmod.f90 preset
# block and array declarations. Recorded in the manifest so a run says what its
# land surface actually was rather than leaving it to be inferred.
UNIFORM_LAND_DEFAULTS = {
    "dz0clim_roughness_m": 2.0,
    "dwmax_field_capacity_m": "wsmax (Earth value)",
    "dalbcl_background_albedo": 0.22,
    "dforest_fraction": 0.5,
    "dglac_glacier_mask": 0.0,
}


def surface_field_report(run_dir: Path, config: dict) -> dict:
    """Record which surface fields are set from file and which are uniform.

    `configure()` runs `rm workdir/*.sra` whenever a landmap or topomap is
    given (__init__.py:2938) and then writes only those two. So every other
    surface field is absent, and surfmod reads each one with
    `inquire(exist=)` and skips it silently when it is (surfmod.f90:125).

    That is not a failure. landmod presets the fields to uniform namelist
    values before reading, so the fallback is 2.0 m roughness, 0.22 albedo,
    Earth field capacity, 0.5 forest fraction and no glaciers, not zeros. For a
    custom world that is a defensible baseline and better than imprinting
    Earth's geography. It just has to be a declared choice rather than an
    accident, which is what `model.uniform_land_surface` is for.

    Because the wipe happens regardless of resolution, the shipped-data limit
    at T42 does not constrain us: a T85 run would default identically.
    """
    nlat = int(config["model"]["latitudes"])
    present = set()
    for f in run_dir.glob(f"N{nlat:03d}_surf_*.sra"):
        try:
            present.add(int(f.stem.rsplit("_", 1)[1]))
        except (IndexError, ValueError):
            continue

    missing = sorted(intended_surface_codes(config) - present)
    if missing:
        raise RuntimeError(
            f"{run_dir} is missing surface fields we supply ourselves: "
            f"{missing}. surfmod skips these silently, so the run would use "
            "flat terrain or an all-ocean mask without saying so."
        )

    declared = bool(config["model"].get("uniform_land_surface", False))
    defaulted = sorted(present - intended_surface_codes(config))
    if not declared:
        raise RuntimeError(
            "Every surface field other than topography and the land mask will "
            "fall back to a uniform namelist default, because configure() "
            "clears them. Set model.uniform_land_surface: true in "
            "config/planet.yaml to declare that this is intended."
        )
    return {
        "from_file": sorted(present),
        "unexpected_extra_files": defaulted,
        "uniform_defaults_declared": declared,
        "uniform_values": UNIFORM_LAND_DEFAULTS,
        "note": ("Fields absent from the run directory are preset to uniform "
                 "namelist values by landmod and never read from file."),
    }


def run_id(config: dict, flux_ratio: float) -> str:
    p = config["planet"]
    a = config["atmosphere"]
    m = config["model"]
    identifier = (
        f"{str(m['resolution']).lower()}l{int(m['layers'])}p{int(m['ncpus'])}"
        f"r{int(m['precision_bytes'])}"
        f"_s{round(100 * flux_ratio):03d}"
        f"_co2{round(1e6 * float(a['pCO2_bar'])):04d}ppm"
        f"_rot{float(p['rotation_hours']):g}h"
        f"_obl{float(p['obliquity_degrees']):g}"
        f"_e{round(1000 * float(p['eccentricity'])):03d}"
        # Physics switches that change the answer get a marker, so an on/off
        # comparison cannot land in one directory. Only non-defaults are named,
        # to keep the identifier readable.
        + ("_glac" if config["surface"].get("glaciers", {}).get("enabled") else "")
        + f"_g{geography_tag(config)}"
    )
    return identifier.replace(".", "p")


def validate_outputs(run_dir: Path) -> dict:
    regular = run_dir / "MOST.00000.nc"
    snapshot = run_dir / "snapshots" / "MOST_SNAP.00000.nc"
    required_regular = {
        "ts", "tas", "pr", "ua", "va", "snd", "sic",
        "sit", "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg", "nu",
    }
    required_snapshot = {
        "ts", "tas", "pr", "ua", "va", "snd", "sic",
        "sit", "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg", "nu",
    }
    result = {}
    for label, path, required in [
        ("regular", regular, required_regular),
        ("snapshot", snapshot, required_snapshot),
    ]:
        if not path.is_file():
            raise RuntimeError(f"Missing {label} output: {path}")
        with Dataset(path) as dataset:
            present = set(dataset.variables)
            missing = sorted(required - present)
            if missing:
                raise RuntimeError(f"{path} is missing required fields: {missing}")
            nonfinite = [
                name for name in sorted(required)
                if not np.isfinite(dataset.variables[name][:]).all()
            ]
            if nonfinite:
                raise RuntimeError(f"{path} has non-finite required fields: {nonfinite}")
            result[label] = {
                "path": str(path),
                "sha256": file_sha256(path),
                "time_samples": len(dataset.dimensions["time"]),
                "required_fields_present": sorted(required),
            }
    return result


def read_sra(path: Path, expected_code: int, nlat: int, nlon: int) -> np.ndarray:
    with path.open("r", encoding="ascii") as handle:
        header = np.fromstring(handle.readline(), sep=" ", dtype=np.int64)
        values = np.loadtxt(handle).reshape(-1)
    expected = np.array([expected_code, 0, 20260811, 0, nlon, nlat, 0, 0])
    if not np.array_equal(header, expected):
        raise ValueError(f"Unexpected SRA header in {path}: {header.tolist()}")
    if values.size != nlat * nlon or not np.isfinite(values).all():
        raise ValueError(f"Invalid SRA payload in {path}")
    return values.reshape(nlat, nlon)


def namelist_value(path: Path, key: str) -> str:
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip().upper()
        if stripped.startswith(key.upper() + " ") or stripped.startswith(key.upper() + "="):
            return line.split("=", 1)[1].strip()
    raise KeyError(f"{key} not found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--flux-ratio", type=float, default=None)
    parser.add_argument(
        "--run-years", type=int, default=0,
        help="Run this many model orbits after preparation",
    )
    parser.add_argument(
        "--force-prepare", action="store_true",
        help="Allow re-preparing an existing run with no climate outputs",
    )
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    flux_ratio = float(
        config["orbit"]["baseline_flux_earth"]
        if args.flux_ratio is None else args.flux_ratio
    )
    if flux_ratio <= 0:
        raise ValueError("Flux ratio must be positive")
    derived = derive(config, flux_ratio)
    identifier = run_id(config, flux_ratio)
    run_dir = (RUNS / identifier).resolve()
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    # Any of these left over from an earlier run makes finalize() and the
    # crash-tolerant rewind pick up the wrong world's state, so refuse the lot
    # rather than only the primary output, and do it regardless of --force-prepare.
    stale = []
    if run_dir.exists():
        for pattern in ("MOST.*", "MOST_REST.*", "MOST*DIAG*",
                        "snapshots/*", "highcadence/*"):
            stale.extend(sorted(run_dir.glob(pattern)))
    if stale:
        raise RuntimeError(
            f"{run_dir} already holds model artifacts ({len(stale)}, e.g. "
            f"{stale[0].name}); this command will not overwrite or implicitly "
            "resume them. Use a fresh directory or move the old run aside."
        )
    if run_dir.exists() and any(run_dir.iterdir()) and not args.force_prepare:
        raise RuntimeError(
            f"{run_dir} already exists and is nonempty; inspect it or pass --force-prepare"
        )

    nlat = int(config["model"]["latitudes"])
    nlon = int(config["model"]["longitudes"])
    landmap = LANDMAP.resolve()
    topomap = TOPOMAP.resolve()
    land = read_sra(landmap, 172, nlat, nlon)
    topo = read_sra(topomap, 129, nlat, nlon)
    if not np.all(np.isin(land, [0.0, 1.0])):
        raise ValueError("Land SRA is not binary")
    if np.any(topo[land == 0.0] != 0.0):
        raise ValueError("Ocean geopotential must be zero")

    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
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
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=int(derived["snapshot_interval_steps"]),
        otherargs={
            "N_DAYS_PER_YEAR@plasim_namelist": str(
                derived["rotations_per_orbit_namelist"]
            )
        },
    )
    # Shipped output lists omit orbital phase and several hydrology fields.
    # This private helper is stable in the pinned release and edits those
    # postprocessor-code lists only.
    model._add_postcodes("example.nl", REGULAR_CODES)
    model._add_postcodes("snapshot.nl", SNAPSHOT_CODES)
    model.cfgpostprocessor(
        ftype="regular",
        extension=model_cfg["output_type"],
        # ExoPlaSim 3.4.2 documents integer codes but silently drops derived
        # variables when integers are used. String codes take the correct path.
        variables=[str(code) for code in REGULAR_CODES],
        mode="grid",
        times=int(model_cfg["regular_output_bins_per_orbit"]),
        timeaverage=True,
        interpolatetimes=False,
    )
    model.cfgpostprocessor(
        ftype="snapshot",
        extension=model_cfg["output_type"],
        variables=[str(code) for code in SNAPSHOT_CODES],
        mode="grid",
        times=None,
        timeaverage=False,
        interpolatetimes=False,
    )
    staged = stage_surface_extras(run_dir, config)
    model.exportcfg(str(run_dir / f"{identifier}.cfg"))
    surface_report = surface_field_report(run_dir, config)

    checks = {
        "N_RUN_STEPS": namelist_value(run_dir / "plasim_namelist", "N_RUN_STEPS"),
        "N_DAYS_PER_YEAR": namelist_value(run_dir / "plasim_namelist", "N_DAYS_PER_YEAR"),
        "NSTPS": namelist_value(run_dir / "plasim_namelist", "NSTPS"),
        "SIDEREAL_YEAR": namelist_value(run_dir / "planet_namelist", "SIDEREAL_YEAR"),
        "ROTSPD": namelist_value(run_dir / "planet_namelist", "ROTSPD"),
        "GSOL0": namelist_value(run_dir / "planet_namelist", "GSOL0"),
        "ECCEN": namelist_value(run_dir / "planet_namelist", "ECCEN"),
        "OBLIQ": namelist_value(run_dir / "planet_namelist", "OBLIQ"),
        "STARBBTEMP": namelist_value(run_dir / "radmod_namelist", "STARBBTEMP"),
        "NSIMPLEALBEDO": namelist_value(run_dir / "radmod_namelist", "NSIMPLEALBEDO"),
    }
    # ExoPlaSim copies its whole run directory in, so every previously built
    # executable is present. Name the one this run will actually use rather than
    # taking the last glob match, which sorts p8 after p16.
    exe_path = run_dir / (
        f"most_plasim_t{int(str(model_cfg['resolution']).lstrip('Tt'))}"
        f"_l{int(model_cfg['layers'])}_p{int(model_cfg['ncpus'])}.x"
    )
    if not exe_path.is_file():
        raise RuntimeError(f"expected executable {exe_path} is not in the run directory")

    manifest = {
        "schema_version": 1,
        "run_id": identifier,
        "status": "prepared",
        "config_path": str(config_path),
        "config_sha256": file_sha256(config_path),
        "source_config": config,
        "derived_parameters": derived,
        "geography": {
            "landmap": str(landmap),
            "landmap_sha256": file_sha256(landmap),
            "topomap": str(topomap),
            "topomap_sha256": file_sha256(topomap),
        },
        # ExoPlaSim names the binary most_plasim_t<res>_l<layers>_p<ncpus>.x, with
        # no precision in the path, and reuses whatever is already there rather
        # than recompiling. So the compiled precision is invisible shared state.
        # Recording the binary's digest is what makes a silent swap auditable.
        "executable": {
            "path": str(exe_path) if exe_path else None,
            "sha256": file_sha256(exe_path) if exe_path and exe_path.is_file() else None,
            "note": ("Built precision is not encoded in the filename. If "
                     "model.precision_bytes changes, pass recompile=True or "
                     "remove the binary, or the old one is silently reused."),
        },
        "software": {
            "python": platform.python_version(),
            "exoplasim": getattr(exo, "__version__", "3.4.2"),
            "numpy": np.__version__,
            "gfortran": command_version(["gfortran", "--version"]),
        },
        "namelist_checks": checks,
        "surface_fields": surface_report,
        "surface_fields_staged": staged,
        "postprocessor": {
            "regular_codes": REGULAR_CODES,
            "snapshot_codes": SNAPSHOT_CODES,
        },
    }
    manifest_path = run_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"run_dir": str(run_dir), "derived": derived, "namelists": checks},
            indent=2,
        )
    )

    if args.run_years:
        try:
            model.run(years=args.run_years, crashifbroken=True, clean=True)
            manifest["output_validation"] = validate_outputs(run_dir)
            manifest["status"] = (
                "smoke_complete" if args.run_years == 1 else "run_complete"
            )
            manifest["completed_orbits"] = args.run_years
        except Exception:
            manifest["status"] = "failed"
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            raise
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
