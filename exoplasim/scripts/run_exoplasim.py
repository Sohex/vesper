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
import uuid

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

# High-cadence output is a gust distribution and nothing else. At one sample
# every fourth timestep a full field list runs to tens of gigabytes an orbit for
# variables nothing reads, so this is the near-surface wind alone: `spd` is what
# the Weibull is fitted to, and `ua`/`va` are kept because a direction is what
# distinguishes a real gust from a reversing mean. DUST-5.
HIGH_CADENCE_CODES = [131, 132, 259]


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


def surface_sra(config: dict, code: int) -> Path:
    """Path to a generated surface field, keyed by the configured resolution."""
    resolution = str(config["model"]["resolution"]).upper()
    return INPUTS / resolution.lower() / f"orogen_{resolution}_surf_{code:04d}.sra"


def stellar_spectrum_path(config: dict) -> str | None:
    """Resolve the configured stellar spectrum to a file ExoPlaSim will accept.

    `configure(starspec=...)` needs a real path and both the file and its
    `_hr.dat` companion; it does not search the package. Returning None falls
    back to the blackbody at `startemp`.

    `inputs/stellarspectra/` is searched before the ExoPlaSim package, because
    spectra we generate are tracked here and the package lives in an untracked
    `.venv` that any reinstall resets. It also lets a project spectrum shadow a
    package one of the same name rather than silently losing to it.
    """
    name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        return None
    roots = (
        INPUTS / "stellarspectra",
        Path(exo.__file__).resolve().parent / "stellarspectra",
    )
    for base in roots:
        path = base / f"{name}.dat"
        companion = base / f"{name}_hr.dat"
        if path.is_file() and companion.is_file():
            return str(path)
    available = sorted(
        {f.stem for base in roots for f in base.glob("*.dat")
         if not f.stem.endswith("_hr")}
    )
    raise RuntimeError(
        f"stellar spectrum {name!r} needs {name}.dat and {name}_hr.dat in one of "
        f"{[str(r) for r in roots]}, which do not both exist. Available: {available}"
    )


def stage_stellar_spectrum(model, run_dir: Path, spectrum: str | None) -> str | None:
    """Copy the spectrum into the run directory and shorten the namelist entry.

    `radmod.f90` declares `starfile` and `starfilehr` as `character(len=80)`,
    and `configure()` absolutises whatever path it is given. A venv path is
    easily longer than that: ours is 90 characters, and Fortran truncates it at
    80, emits a namelist *warning* rather than an error, then dies with an
    end-of-file inside `readdat`. The failure names a path that does not exist,
    which is not obviously a length problem.

    The model runs with the run directory as its working directory, so staging
    the two files there and referring to them by bare name sidesteps the limit
    entirely.
    """
    if not spectrum:
        return None
    src = Path(spectrum)
    hires = src.with_name(f"{src.stem}_hr.dat")
    for f in (src, hires):
        shutil.copyfile(f, run_dir / f.name)
    model._edit_namelist("radmod_namelist", "STARFILE", f"'{src.name}'")
    model._edit_namelist("radmod_namelist", "STARFILEHR", f"'{hires.name}'")
    return src.name


def surface_input_paths(config: dict) -> list[Path]:
    """Every SRA file that defines this run's surface, in a stable order."""
    return [surface_sra(config, code)
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
ALBEDO_SURFACE_CODES = {174, 175, 176, 212}
# dwmax, the soil water bucket whose overflow *is* ExoPlaSim's runoff. Supplied
# from pedology when asked for; otherwise the uniform namelist default stands.
SOIL_WATER_SURFACE_CODES = {229}
ROUGHNESS_SURFACE_CODES = {173}

# PlaSim's own 28-term energy decomposition, denergy(NHOR,28), written to these
# codes when nenergy > 0. The instrument for the constant -0.455 W/m2 that does
# not close between the top of the atmosphere and the surface; see
# notes/water-and-energy-closure.md.
#
# Term 15 needs patches/exoplasim-3.4.2-energy-diagnostics.patch applied and
# ExoPlaSim rebuilt, or it reports the latent heating with sublimation and
# vaporisation swapped.
ENERGY_DIAGNOSTIC_CODES = list(range(360, 388))

# The same 28 terms resolved per level, under nener3d. A column total says the
# atmosphere is losing energy; these say at which level, which is the difference
# between a number and a lead.
ENERGY_3D_CODES = list(range(460, 488))

# Which module each denergy term is accumulated in, from the assignment sites in
# the PlaSim source. Attribution only: the individual terms are not separately
# documented upstream, and naming them by physics would be inventing detail that
# has not been traced. Module plus index is enough to say where a residual is
# being introduced, which is the question being asked.
ENERGY_TERM_MODULE = {
    1: "plasim", 2: "plasim", 3: "plasim", 4: "plasim", 5: "plasim",
    6: "miscmod", 7: "fluxmod", 8: "fluxmod",
    9: "radmod", 10: "radmod",
    11: "rainmod", 12: "rainmod", 13: "rainmod", 14: "rainmod",
    15: "rainmod", 16: "rainmod",
    17: "radmod", 18: "radmod", 19: "radmod", 20: "radmod",
    21: "fluxmod", 22: "fluxmod",
    23: "plasim", 24: "plasim", 25: "plasim",
    26: "plasim", 27: "plasim", 28: "radmod",
}


def register_energy_diagnostic_codes() -> int:
    """Teach the postprocessor about codes 360-387 before it runs.

    pyburn resolves every requested code against its own `ilibrary`, which ships
    119 entries and none of these. Asking for them produces a run that completes,
    writes a restart, and then dies in postprocessing with "Going to stop here
    just in case", naming neither the code nor the reason. So enabling the
    diagnostics without this crashes at the last step of an otherwise good run.

    `ilibrary` is consulted at call time rather than captured at import, so
    inserting here is enough and the vendored tree stays untouched. That matters:
    `.venv` is untracked and any reinstall would silently drop a patch, whereas
    this travels with the script that depends on it.
    """
    from exoplasim import pyburn
    for i, code in enumerate(ENERGY_DIAGNOSTIC_CODES, start=1):
        module = ENERGY_TERM_MODULE.get(i, "unknown")
        pyburn.ilibrary.setdefault(str(code), [
            f"denergy{i:02d}",
            f"plasim_energy_budget_term_{i:02d}_{module}",
            "W m-2",
        ])
    for i, code in enumerate(ENERGY_3D_CODES, start=1):
        module = ENERGY_TERM_MODULE.get(i, "unknown")
        pyburn.ilibrary.setdefault(str(code), [
            f"dener3d{i:02d}",
            f"plasim_energy_budget_term_{i:02d}_{module}_per_level",
            "W m-2",
        ])
    return len(ENERGY_DIAGNOSTIC_CODES) + len(ENERGY_3D_CODES)


def intended_surface_codes(config: dict) -> set[int]:
    """Which surface fields this run supplies rather than leaving at defaults.

    `model.soil_water_source` defaults to `uniform` when the key is absent, so a
    config predating this option behaves exactly as it did.

    Setting it to `pedology` makes 229 mandatory here, and 229 comes from a soil
    weathered under a climatology. On a terrain that has none yet, the bootstrap
    run has to go out with the key at `uniform`; see WORKFLOW section 6A, which
    also says why the flip cannot happen mid-run.
    """
    codes = set(BASE_SURFACE_CODES)
    if str(config["model"].get("land_albedo_source", "uniform")) != "uniform":
        codes |= ALBEDO_SURFACE_CODES
    if str(config["model"].get("soil_water_source", "uniform")) != "uniform":
        codes |= SOIL_WATER_SURFACE_CODES
    if str(config["model"].get("roughness_source", "uniform")) != "uniform":
        codes |= ROUGHNESS_SURFACE_CODES
    return codes


def energy_diagnostics_enabled(config: dict) -> bool:
    """Whether to ask PlaSim for its 28-term energy decomposition.

    Defaults False when the key is absent, so a config predating this option
    behaves exactly as it did and `config_sha256` does not move. Turning it on
    adds 28 output fields and is meant for short diagnostic segments, not for
    production runs.
    """
    return bool(config["model"].get("energy_diagnostics", False))


def disable_low_io(model) -> None:
    """Turn off PlaSim's low-I/O output path, which writes a corrupt first bin.

    `NLOWIO = 1` is PlaSim's default (`plasimmod.f90:151`) and accumulates fields
    over each output interval, dividing in place at write time. The first record
    of every model call comes out wrong under it: the free troposphere is right
    to 5% and the boundary layer is missing entirely, so bottom-level wind reads
    7.5x the other bins and humidity 27% low. Scalars are within 2%.

    Proven to be this path rather than the postprocessor by running one orbit at
    `NLOWIO = 0` through the same `pyburn` averaging: the humidity ratio goes
    from 0.714 to 0.930, which is inside the ordinary seasonal spread. See
    `exoplasim/notes/first-output-bin.md`.

    The exact line is NOT pinned, and a patch is deliberately not written on the
    obvious candidate: `naccuout` persists across runs through the restart while
    the accumulators do not (`plasim.f90:767`), but a counter error scales a
    field uniformly and this one changes its vertical structure. So the setting
    is turned off rather than the bug fixed; CLIM-5 carries the patch.

    What it costs: about 2.4 GB per orbit against 96 MB, because output becomes
    182 instantaneous records per orbit, one roughly every 24 hours, instead of
    12 accumulated bins. `pyburn` still averages them to 12 for the .nc, so
    nothing downstream changes shape. Samples are also strictly more information
    than an accumulation, since an average can be recomputed from them and an
    accumulation cannot be undone -- which is what makes a proper mean wind
    speed recoverable at all.
    """
    model._edit_namelist("plasim_namelist", "NLOWIO", "0")


def enable_energy_diagnostics(model, config: dict) -> bool:
    """Set nenergy in plasim_nl, which the Python API does not expose.

    Same mechanism `stage_stellar_spectrum` uses for STARFILE: the namelist is
    edited directly because `configure()` has no parameter for it.
    """
    if not energy_diagnostics_enabled(config):
        return False
    model._edit_namelist("plasim_namelist", "NENERGY", "1")
    if config["model"].get("energy_diagnostics_3d", False):
        model._edit_namelist("plasim_namelist", "NENER3D", "1")
    return True


def stage_surface_extras(run_dir: Path, config: dict) -> list[int]:
    """Copy the surface fields we generate into the run directory.

    Has to run after `configure()`, which does `rm workdir/*.sra` whenever a
    landmap or topomap is given and then writes back only those two
    (`__init__.py:2938`). Anything staged before that call is silently deleted.
    """
    staged = []
    for code in sorted(intended_surface_codes(config) - BASE_SURFACE_CODES):
        src = surface_sra(config, code)
        if not src.is_file():
            builder = ("build_surface_soil_water.py"
                       if code in SOIL_WATER_SURFACE_CODES
                       else "build_surface_albedo.py")
            setting = ("model.soil_water_source"
                       if code in SOIL_WATER_SURFACE_CODES
                       else "model.land_albedo_source")
            raise RuntimeError(
                f"{src} is missing. Run {builder}, or set {setting}: uniform "
                f"to accept ExoPlaSim's namelist default."
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


def spectrum_tag(config: dict) -> str:
    """Marker naming the stellar spectrum in a run directory.

    The spectrum sets the weighting for every snow, ice and glacier albedo, so
    two runs differing only in it are different climates. Until the k2/K2.5V
    correction there was nothing here, and a re-baseline would have landed in the
    completed run's directory: the same silent-overwrite that `geography_tag`
    exists to prevent, reached by a different route.

    Directories written before this was added carry no marker and are all `k2`.
    Recomputing an id for one now yields a name that does not exist on disk, so a
    continuation of a pre-fix run fails loudly rather than resuming the wrong
    world. That is the intended direction to fail in.
    """
    name = config.get("radiation", {}).get("stellar_spectrum")
    return f"_{name}" if name else "_bb"


def flux_tag(flux_ratio: float) -> str:
    """Flux in the run id, at whatever precision the value actually needs.

    This was round(100 * flux), so 0.945 and 0.94 both produced `s094` and two
    different climates would have shared a directory. The geography digest
    happened to separate them the first time only because the terrain changed in
    the same step, which is luck rather than a guard.

    Thousandths, with a single trailing zero dropped, so every id written under
    the old rule is reproduced exactly: 0.90 -> 090, 0.96 -> 096, 1.00 -> 100,
    and 0.945 -> 0945. Existing run directories stay findable.
    """
    tag = f"{round(flux_ratio * 1000):04d}"
    return tag[:-1] if tag.endswith("0") else tag


def physical_fingerprint(config: dict, flux_ratio: float) -> dict:
    """Everything about a run that makes it a different climate.

    Recorded in the manifest and indexed by `index_runs.py` so a run is findable
    by what it is. It is deliberately NOT the run's identity -- see `run_id`.
    """
    p = config["planet"]
    a = config["atmosphere"]
    m = config["model"]
    return {
        "resolution": str(m["resolution"]),
        "layers": int(m["layers"]),
        "ranks": int(m["ncpus"]),
        "precision_bytes": int(m["precision_bytes"]),
        "flux_ratio": round(float(flux_ratio), 6),
        "co2_ppm": round(1e6 * float(a["pCO2_bar"]), 3),
        "rotation_hours": float(p["rotation_hours"]),
        "obliquity_degrees": float(p["obliquity_degrees"]),
        "eccentricity": float(p["eccentricity"]),
        "glaciers": bool(config["surface"].get("glaciers", {}).get("enabled")),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "geography": geography_tag(config),
    }


def run_id(config: dict, flux_ratio: float) -> str:
    """A UUID. Not derived from anything.

    This used to spell out the physical parameters --
    `t42l10p16r8_s0968_co20450ppm_rot30h_obl32_e020_glac_k25v_g83d1b976` -- so
    that two different climates could not share a directory, which matters
    because ExoPlaSim's `finalize()` picks output as the last glob match.

    **A derived identifier can only separate runs along the dimensions it
    encodes, and that is not a property you can maintain.** The ozone band-weight
    patch changed the physics and moved nothing in the encoded set, so the
    pre-patch and post-patch runs at the same flux computed the same identifier
    and landed in the same directory. Adding the patch stack to the name would
    have fixed that instance and not the next one; the encoded set is a list of
    everything someone has thought of so far.

    A UUID is unique unconditionally, so no run can ever collide with another
    regardless of what changed between them, including things nothing here
    models. What the run *was* belongs in the manifest, which is written anyway,
    and in `exoplasim/runs/INDEX.json`, which is generated from the manifests.

    The cost is that a continuation cannot recompute the directory name and must
    be given it. That is the safer direction: recomputation could silently
    resolve to a different run, and being handed an id cannot.
    """
    return "run_" + uuid.uuid4().hex[:12]


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
    parser.add_argument(
        "--restart-from", type=Path, default=None,
        help="Seed the initial state from an existing MOST_REST file instead of "
             "cold-starting. Only the spin-up path changes, not the equilibrium.",
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
    # Announce the id on its own line, in a machine-readable form. It is a UUID,
    # so a caller cannot reconstruct it, and the alternative -- picking the
    # newest or first run_* directory afterwards -- is selection by sort order,
    # which is the pattern this project has now been bitten by four times.
    print(f"RUN_ID={identifier}", flush=True)
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
    landmap = surface_sra(config, 172).resolve()
    topomap = surface_sra(config, 129).resolve()
    land = read_sra(landmap, 172, nlat, nlon)
    topo = read_sra(topomap, 129, nlat, nlon)
    if not np.all(np.isin(land, [0.0, 1.0])):
        raise ValueError("Land SRA is not binary")
    if np.any(topo[land == 0.0] != 0.0):
        raise ValueError("Ocean geopotential must be zero")

    # Seeding the initial state from a nearby equilibrium is worth about 20 of
    # the 70 orbits a cold start spends climbing from its 269 K initial state.
    # It is an initial condition and nothing else: the equilibrium a run settles
    # to is set by its forcing, so this changes the path and not the answer. The
    # source is recorded in the manifest because the path is no longer a function
    # of the configuration alone.
    restart_seed = None
    if args.restart_from is not None:
        restart_seed = args.restart_from.resolve()
        if not restart_seed.is_file():
            raise RuntimeError(f"--restart-from {restart_seed} does not exist")
        # A restart FREEZES every land surface boundary condition.
        #
        # landmod.f90's landini reads dwmax, dz0clim, dz0climo, dalbcl, dalbcl1
        # and dalbcl2 from the restart file when restart > 0, and from the
        # surface .sra files only on a cold start. So seeding from a restart
        # silently discards any surface field that has changed since -- soil
        # water, roughness and all three albedo bands.
        #
        # This was found by running exactly that: a baseline seeded from the
        # bootstrap, with a new soil-water map and lake-composited albedo, would
        # have reproduced the bootstrap and shown that lakes and soil water do
        # nothing. It only surfaced because it also happened to trap a SIGFPE in
        # landini. A silent null result is the failure mode this guard exists to
        # prevent.
        src_manifest = restart_seed.parent / "run_manifest.json"
        if src_manifest.is_file():
            src = json.loads(src_manifest.read_text(encoding="utf-8"))
            was = set((src.get("surface_fields") or {}).get("from_file") or [])
            now = intended_surface_codes(config)
            reason = None
            if was != now:
                reason = (f"codes differ: {sorted(was)} then, {sorted(now)} now")
            else:
                # Same codes, possibly different content.
                old_h = src.get("surface_field_sha256") or {}
                changed = [c for c in sorted(now)
                           if surface_sra(config, c).is_file()
                           and old_h.get(str(c)) not in
                           (None, file_sha256(surface_sra(config, c)))]
                if changed:
                    reason = f"content changed for code(s) {changed}"
                elif not old_h:
                    reason = ("the source run recorded no surface field hashes, "
                              "so content changes cannot be ruled out")
            if reason:
                raise RuntimeError(
                    "--restart-from refused: " + reason + ". A restart reads soil "
                    "water, roughness and albedo from ITSELF, not from the .sra "
                    "files -- landmod's landini takes dwmax, dz0clim and all "
                    "three dalbcl bands from the restart when restart > 0 -- so "
                    "those changes would be silently discarded and the run would "
                    "reproduce its parent. Cold-start instead.")

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
        restartfile=None if restart_seed is None else str(restart_seed),
        flux=derived["stellar_flux_w_m2"],
        startemp=float(star["effective_temperature_k"]),
        starspec=stellar_spectrum_path(config),
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
    regular_codes = list(REGULAR_CODES)
    if energy_diagnostics_enabled(config):
        regular_codes = regular_codes + ENERGY_DIAGNOSTIC_CODES
        if config["model"].get("energy_diagnostics_3d", False):
            regular_codes = regular_codes + ENERGY_3D_CODES
    model._add_postcodes("example.nl", regular_codes)
    model._add_postcodes("snapshot.nl", SNAPSHOT_CODES)
    model.cfgpostprocessor(
        ftype="regular",
        extension=model_cfg["output_type"],
        # ExoPlaSim 3.4.2 documents integer codes but silently drops derived
        # variables when integers are used. String codes take the correct path.
        variables=[str(code) for code in regular_codes],
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
    spectrum = stage_stellar_spectrum(model, run_dir, stellar_spectrum_path(config))
    # Ozone column scaling, set through the namelist because the Python API does
    # not expose o3scale. The model prescribes an Earth column and derives
    # nothing about it from the host star; Segura et al. (2003) measure 0.794 of
    # Earth's for a K2V host at 1 PAL O2.
    o3 = config["model"].get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        model._edit_namelist("radmod_namelist", "O3SCALE", f"{float(o3)}")
        print(f"ozone column scaled to {float(o3)} of Earth's (Segura et al. 2003)")
    for key, name in (("ozone_uv_weight", "O3UVW"), ("ozone_visible_weight", "O3VISW")):
        w = config["model"].get(key)
        if w is not None and float(w) != 1.0:
            model._edit_namelist("radmod_namelist", name, f"{float(w)}")

    disable_low_io(model)
    if enable_energy_diagnostics(model, config):
        n = register_energy_diagnostic_codes()
        print(f"energy diagnostics on: nenergy=1, {n} codes 360-387 registered "
              "with the postprocessor. Term 15 needs "
              "patches/exoplasim-3.4.2-energy-diagnostics.patch and a rebuild.")
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
        # What this run IS, as opposed to what it is called. The id is a UUID and
        # carries no meaning, so this is the only place the physics is written
        # down in a machine-readable form; exoplasim/scripts/index_runs.py builds
        # exoplasim/runs/INDEX.json from it.
        "physical": physical_fingerprint(config, flux_ratio),
        # Which terrain this run is of, by name, at the top level rather than
        # only inside the config copy below. Everything downstream verifies build
        # identity through lib/provenance.py, and build_climatology was looking
        # for it here and finding nothing, so every climatology went out
        # unstamped and every check on it went quiet.
        "source_build": config.get("source_build"),
        "status": "prepared",
        "config_path": str(config_path),
        "config_sha256": file_sha256(config_path),
        "source_config": config,
        "derived_parameters": derived,
        "initial_state": ({"cold_start": True} if restart_seed is None else {
            "cold_start": False,
            "restart_from": str(restart_seed),
            "restart_from_sha256": file_sha256(restart_seed),
            "restart_from_run": restart_seed.parent.name,
            "note": "Initial condition only; equilibrium is set by the forcing.",
        }),
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
        # Per-code hashes, so a surface field that CHANGED CONTENT under the same
        # code is detectable. The set of codes alone is not enough: the albedo
        # rebuilt with lakes composited is still code 174, and a restart would
        # have discarded it silently.
        "surface_field_sha256": {
            str(c): file_sha256(surface_sra(config, c))
            for c in sorted(intended_surface_codes(config))
            if surface_sra(config, c).is_file()},
        "stellar_spectrum": spectrum,
        "postprocessor": {
            "regular_codes": regular_codes,
            "energy_diagnostics": energy_diagnostics_enabled(config),
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
