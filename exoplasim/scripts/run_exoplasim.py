#!/usr/bin/env python3
"""Prepare, validate, and optionally run a reproducible ExoPlaSim experiment."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
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

from _paths import CONFIG, INPUTS, PROJECT_ROOT, RUNS
import reset_restart_accumulators


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

# High-cadence output is a gust distribution and nothing else, so this is the
# near-surface wind alone: `spd` is what the Weibull is fitted to, and `ua`/`va`
# are kept because a direction is what distinguishes a real gust from a reversing
# mean. DUST-5.
#
# This trims the POSTPROCESSED file. The raw one the model writes is 15 GB for a
# single T42 orbit at one sample every fourth timestep regardless of this list,
# because the model's output path takes no field list; it is cleaned up after
# pyburn runs, but the disk has to be free first.
HIGH_CADENCE_CODES = [131, 132, 259]

# The shortwave gas band weights, and the ONE list of them. Lacis and Hansen's
# absorptances are fractions of SOLAR flux, so a non-solar host needs every one
# re-weighted; radmod.f90 carries a key per term and this is where the config
# reaches them.
#
# Each is written only when it differs from the model's OWN default, so a run's
# namelist says what departs from the scheme as shipped -- and the defaults are
# not all 1.0. co2sww's is 0.0, because upstream has no shortwave CO2 term at
# all and zero is what reproduces upstream; comparing it against 1.0 would
# silently drop a weight of 1.0 and silently write one of 0.0, which is the
# inversion of what is meant.
#
# IT LIVES HERE BECAUSE continue_exoplasim.py IMPORTS IT. configure() rewrites
# the namelist on every continuation, so each key has to be reapplied per
# segment, and for a while that was a second copy of this tuple. PHYS-9 then
# added h2o_sw_level to this copy and not to that one, so H2OSWL applied for the
# orbits run_exoplasim.py prepared and silently reverted to the model default on
# every continuation after -- a physics change partway through a run, which is
# the failure mode both copies' comments called the hardest to notice in a long
# spin-up. One list cannot drift from itself.
SHORTWAVE_GAS_KEYS = (
    ("ozone_uv_weight", "O3UVW", 1.0),
    ("ozone_visible_weight", "O3VISW", 1.0),
    ("h2o_sw_weight", "H2OSWW", 1.0),
    # A LEVEL, not a weight, and separate from H2OSWW on purpose: that one is a
    # star-over-Sun ratio and an error in Eq. 21's absolute level divides out of
    # it. PHYS-9.
    ("h2o_sw_level", "H2OSWL", 1.0),
    ("co2_sw_weight", "CO2SWW", 0.0),
)

# The LONGWAVE trace gases, CLIM-42. Same one-list rule and the same reason:
# `continue_exoplasim.py` imports this, because a key applied at prepare and not
# reapplied per segment is a physics change partway through a run.
#
# These come from the `atmosphere` block rather than `model`, because they are
# composition and not a scheme weight -- pCO2_bar's neighbours. radmod_nl wants
# a VOLUME MIXING RATIO in ppmv, so the conversion is the partial pressure over
# the total, and the model's own default is 0.0, meaning absent.
TRACE_GAS_KEYS = (
    ("pCH4_bar", "CH4", 0.0),
    ("pN2O_bar", "N2O", 0.0),
)


def trace_gas_ppmv(config: dict) -> dict:
    """{namelist key: ppmv} for the longwave trace gases the config declares."""
    atmosphere = config.get("atmosphere", {})
    total = sum(float(v) for k, v in atmosphere.items()
                if k.startswith("p") and k.endswith("_bar"))
    if total <= 0.0:
        raise SystemExit("atmosphere block declares no partial pressures")
    out = {}
    for key, name, default in TRACE_GAS_KEYS:
        bar = atmosphere.get(key)
        if bar is None:
            continue
        ppmv = 1e6 * float(bar) / total
        if ppmv != default:
            out[name] = ppmv
    return out


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


def configure_otherargs(derived: dict) -> dict:
    """Namelist keys passed through `configure()`, for a prepare OR a resume.

    ONE definition, because constructing the model over an existing run
    directory re-copies the shipped namelists over the configured ones, so
    every one of these has to be reapplied per segment. That was a second copy
    in `continue_exoplasim.py` carrying only N_DAYS_PER_YEAR, so CLIM-17's
    TFREEZE applied to the orbits this script prepared and reverted to
    icemod.f90's compiled Earth value on every continuation after.

    At the declared salinity that particular revert is worth 5 mK and nothing
    else, but the salinity BRACKET arms are worth -0.04 to -0.36 K and are run
    as continuation segments: the key would have silently not applied, and an
    A/B designed to have a right answer would have measured zero and looked
    like it had confirmed one.
    """
    return {
        "N_DAYS_PER_YEAR@plasim_namelist": str(
            derived["rotations_per_orbit_namelist"]
        ),
        # Sea water's freezing point, from the declared salinity rather than
        # from icemod.f90's compiled-in Earth value. `configure()` has no
        # parameter for it and it is an ordinary `icemod_nl` key, so it goes
        # the same route N_DAYS_PER_YEAR does. CLIM-17.
        "TFREEZE@icemod_namelist": f"{derived['sea_water_freezing_point_k']:.4f}",
        # CLIM-16, oceanmod_nl. Written unconditionally, defaults included, so
        # the namelist in the run directory records what the arm actually ran
        # with instead of leaving it to the binary's compiled value.
        "NHDIFF@oceanmod_namelist": str(derived["ocean_horizontal_diffusion"]),
        "HDIFFK@oceanmod_namelist":
            f"{derived['ocean_horizontal_diffusivity_m2_s']:.6g}",
        # PHYS-11, radmod_nl. The compiled values are tswr3 0.0055 and
        # acl2 (0.05, 0.10, 0.20); at scale 1.0 these reproduce them exactly.
        "TSWR3@radmod_namelist":
            f"{0.0055 * derived['cloud_absorption_scale']:.6g}",
        "ACL2@radmod_namelist": ", ".join(
            f"{v * derived['cloud_absorption_scale']:.6g}"
            for v in (0.05, 0.10, 0.20)),
    }


def freezing_point_k(salinity_psu: float) -> float:
    """Freezing point of sea water at the surface, from its salinity.

    UNESCO (1983) / Millero's polynomial, the standard one, at zero gauge
    pressure. It exists because `icemod.f90` hardcodes `TFREEZE = 271.25` with
    the comment "at S=34.7", so the model carries Earth's ocean unless something
    tells it otherwise, and nothing in this project did until CLIM-17.

    The formula is checked against that constant rather than trusted: at
    S = 34.7 it must reproduce 271.25 K, which is a right answer the model
    already knows and therefore a test that can fail.
    """
    s = float(salinity_psu)
    if not 0.0 <= s <= 42.0:
        raise ValueError(f"salinity {s} psu is outside the range the UNESCO "
                         "polynomial is fitted over (0 to 42)")
    celsius = -0.0575 * s + 1.710523e-3 * s**1.5 - 2.154996e-4 * s**2
    return 273.15 + celsius


_REFERENCE_SALINITY = 34.7
_REFERENCE_TFREEZE = 271.25
# Tolerance is the precision the model states the constant to, two decimals, and
# not tighter: the polynomial gives 271.2449 and 271.25 is that rounded. A
# tighter bound fails on the rounding rather than on a disagreement, which it
# did when this was first written at 0.005.
if abs(freezing_point_k(_REFERENCE_SALINITY) - _REFERENCE_TFREEZE) > 0.01:
    raise RuntimeError(
        f"the freezing-point polynomial gives "
        f"{freezing_point_k(_REFERENCE_SALINITY):.4f} K at S = "
        f"{_REFERENCE_SALINITY}, where icemod.f90's TFREEZE says "
        f"{_REFERENCE_TFREEZE}. One of the two is wrong and this is not a "
        "difference to average over.")


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
        "ocean_salinity_psu": float(config["ocean"]["salinity_psu"]),
        "sea_water_freezing_point_k": freezing_point_k(
            config["ocean"]["salinity_psu"]),
        # CLIM-16. Ocean horizontal heat transport EXISTS; a constant
        # diffusivity is a BOUND on the missing transport rather than the
        # transport, and it is bracketed rather than tuned. NLEV_OCE is 1
        # (oceanmod.f90:15) so hdiffk is a one-element array and a scalar
        # assignment fills it.
        "ocean_horizontal_diffusion": int(bool(
            config["ocean"].get("horizontal_diffusion", False))),
        "ocean_horizontal_diffusivity_m2_s": float(
            config["ocean"].get("horizontal_diffusivity_m2_s", 1.0e3)),
        # PHYS-11. Scales the two ABSORPTION-like cloud keys only, tswr3 and
        # the acl2 triplet. The scattering keys are held and acllwr is a
        # thermal-band constant with no stellar dependence, so neither gets an
        # arm on this argument. 1.0 must reproduce the compiled values exactly,
        # which is the arm the prediction says has to come out bit-identical.
        "cloud_absorption_scale": float(
            config["model"].get("cloud_absorption_scale", 1.0)),
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


def verify_stellar_spectrum(model, config: dict) -> None:
    """Refuse to run if the configured spectrum never reached the namelist.

    `radmod.f90:813` takes the spectrum branch only when `NSTARFILE > 0`;
    otherwise `solarini` builds a Planck curve at `STARBBTEMP` and says so
    nowhere except one line of `MOST_DIAG` that nothing reads. Constructing an
    `exo.Earthlike` on an existing run directory re-copies the shipped
    namelists over the configured ones, so any driver that rebuilds the model
    and calls `configure()` without `starspec` reverts the star in silence.

    `continue_exoplasim.py` and `run_stellar_cycle.py` both did. Every run on
    this build integrated its first orbit against `k25v` and every orbit after
    it against a 4965 K blackbody: the model's own log reports
    `Energy fraction below 0.75 microns` as 0.384383 in `MOST_DIAG.00000` and
    0.418350 from `MOST_DIAG.00001` onward, and that fraction weights the snow
    and sea-ice albedo, which moved the broadband snow albedo the model uses
    from 0.5382 to 0.5623.
    """
    name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        return
    namelist = Path(model.workdir) / "radmod_namelist"
    entries = {}
    for line in namelist.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            entries[key.strip().upper()] = value.strip().strip(" ,").strip("'\"")
    if entries.get("NSTARFILE") != "1" or not entries.get("STARFILEHR"):
        raise RuntimeError(
            f"config names stellar spectrum {name!r}, but {namelist} carries "
            f"NSTARFILE={entries.get('NSTARFILE')!r} and "
            f"STARFILEHR={entries.get('STARFILEHR')!r}. The model would run on a "
            "blackbody at STARBBTEMP instead. Pass starspec= to configure() and "
            "call stage_stellar_spectrum() after it.")
    staged = Path(model.workdir) / entries["STARFILEHR"]
    if not staged.is_file():
        raise RuntimeError(
            f"radmod_namelist names {entries['STARFILEHR']} but it is not in "
            f"{model.workdir}; readdat would die at end of file.")
def stellar_spectrum_digest(config: dict) -> dict | None:
    """Content digests of the two spectrum FILES the model will actually read.

    The config names a spectrum; the model reads a file. `build_stellar_spectrum.py`
    writes `<name>.dat` and `<name>_hr.dat` from `star.spectral_type` and
    `star.effective_temperature_k`, IN PLACE and under the same name, so the
    whole radiative input can change while every recorded name stays `k25v`.
    That difference is worth +0.024 on broadband snow albedo at the size it has
    already been measured at, and it moves every derived surface albedo with it.
    A name is a label here; the sha256 is the identity.

    `stellar_spectrum_path` has already proved both files exist, so this does
    not have to handle a half-present pair. Returns None only when the config
    names no spectrum, which means a blackbody at `startemp` and nothing to
    hash.

    This is the run-manifest twin of `spectrum_sha256` in
    `biosphere/generated/vesper_provenance.json`, deliberately: the biosphere
    half of this was already safe, and one mechanism recorded in two places
    beats a second mechanism.
    """
    path = stellar_spectrum_path(config)
    if not path:
        return None
    src = Path(path)
    hires = src.with_name(f"{src.stem}_hr.dat")
    return {
        "name": config.get("radiation", {}).get("stellar_spectrum"),
        "file": src.name,
        "sha256": file_sha256(src),
        "hr_file": hires.name,
        "hr_sha256": file_sha256(hires),
    }


def require_stellar_spectrum(manifest: dict, config: dict) -> bool:
    """Refuse to continue a run whose spectrum FILE has changed under it.

    The resume guard compares parsed config values, and `star.spectral_type` is
    on its inert list. That entry is CORRECT about the config key -- no script
    passes the spectral type to the model -- and WRONG about the artifact,
    because the value reaches the radiation through `k25v.dat`, which the key
    only names. Comparing the file is what closes that route; loosening the
    config comparison would not, since the spectrum can also be regenerated with
    no config edit at all.

    Returns True when the manifest was BACKFILLED and the caller must write it
    out. A run prepared before this existed carries no digest, and refusing
    those outright would make every existing run unresumable for a defect they
    predate. So an unstamped run is stamped with what it is about to run on,
    said out loud, and guarded from the next resume onward -- the same direction
    `lib/provenance.py:require_build` takes with `allow_unstamped`. Segments
    record their own digest either way, so which orbits are covered stays
    readable rather than being inferred from the top-level value.
    """
    current = stellar_spectrum_digest(config)
    recorded = manifest.get("stellar_spectrum_digest")
    if recorded is None:
        manifest["stellar_spectrum_digest"] = current
        print("  warning: this run was prepared before its spectrum was recorded "
              "by content, so the orbits already in it cannot be checked against "
              f"{(current or {}).get('file', 'a blackbody')}. Stamping it now; "
              "the next resume is guarded.")
        return True
    if recorded == current:
        return False
    differing = sorted(
        key for key in ("name", "file", "sha256", "hr_file", "hr_sha256")
        if (recorded or {}).get(key) != (current or {}).get(key))
    raise RuntimeError(
        "The stellar spectrum this run was prepared on is not the one on disk; "
        "refusing to resume:\n  "
        + "\n  ".join(f"{k}: {(recorded or {}).get(k)!r} -> {(current or {}).get(k)!r}"
                      for k in differing)
        + "\nThe spectrum weights every snow, ice and glacier albedo, so orbits "
          "either side of this are different climates. Start a new run rather "
          "than extending this one."
    )


def surface_input_paths(config: dict) -> list[Path]:
    """Every SRA file that defines this run's surface, in a stable order."""
    return [surface_sra(config, code)
            for code in sorted(intended_surface_codes(config))]


def geography_tag(config: dict) -> str:
    """Short digest of every surface input, carried in the run's identity.

    NOT the directory name any more: names are UUIDs, per CLAUDE.md rule 6,
    because a parameter-built name separates runs only along the dimensions it
    encodes and both ExoPlaSim and LPJ-GUESS collided that way. What this value
    is now is `physical_fingerprint["geography"]`, and it is load-bearing there
    -- it is the field that lets `INDEX.json` tell two runs on different surfaces
    apart, which is how the albedo bracket's cases are distinguishable at all.

    What it used to prevent is now prevented by the prepare guard. Two runs
    sharing a config landing in one directory made ExoPlaSim's finalize() select
    output as the last match of sorted(glob("MOST*")), so a shorter new run in a
    directory holding a longer old one copied out the previous run's final year
    under the new name, silently. The guard refuses that case directly.

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
# ddustcol, the prescribed band-1 column dust optical depth. DUST-11. Only read
# by a binary carrying patches/exoplasim-3.4.2-prescribed-dust.patch; an
# unpatched one cannot parse NDUSTRAD and aborts in radini_, which is the loud
# failure this arrangement is designed to produce rather than avoid.
# CLIM-39 made ndustrad a COUNT of prescribed species rather than a switch, and
# the model now reads one surface field per species: species s is code 1810+s,
# up to NAERSP in radmod. Only 1811 is listed here because only 1811 has a
# generator -- `build_surface_dust.py`. Routing a second species in is CLIM-40,
# and it adds the code here in the same commit as the generator that writes it,
# per the rule that an artifact no step produces does not exist.
DUST_SURFACE_CODES = {1811}
# dsrcw, ddrage and dwpr, the source map for the INTERACTIVE emission scheme.
# DUST-3. A different thing from 1811 and not an alternative to it: 1811 is a
# prescribed column optical depth that nothing responds to, while these three
# are a source map that the model's own winds drive. Only read by a binary
# carrying patches/exoplasim-3.4.2-dust-emission.patch; an unpatched one cannot
# parse LDUSTEMIT and aborts in aero_ini, which is the loud failure this
# arrangement is designed to produce rather than avoid.
DUST_EMISSION_SURFACE_CODES = {1801, 1802, 1803}

# PlaSim's own 28-term energy decomposition, denergy(NHOR,28), written to these
# codes when nenergy > 0. Turned on to name the gap between the top of the
# atmosphere and the surface; it closed to 0.02 W/m2 and so ruled itself out, and
# is kept as the audit on the gridpoint physics. See
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


def load_conversion_report(restart: Path, resolved: Path) -> dict | None:
    """The report `convert_restart.py` wrote beside a converted restart.

    Looked for beside the file rather than passed as a flag, so a converted
    state carries its own provenance wherever it is moved to and cannot be
    seeded from without it. The report's own output hash must match the file:
    a report beside a DIFFERENT restart is worse than none, because it would
    vouch for a state it does not describe.
    """
    for candidate in (Path(str(restart) + ".conversion.json"),
                      restart.with_suffix(restart.suffix + ".conversion.json"),
                      restart.with_suffix(".conversion.json")):
        if candidate.is_file():
            report = json.loads(candidate.read_text(encoding="utf-8"))
            break
    else:
        return None
    got = file_sha256(resolved)
    want = (report.get("output") or {}).get("sha256")
    if want != got:
        raise RuntimeError(
            f"{candidate} describes a restart with sha256 {want} and "
            f"{resolved} has {got}. A conversion report beside the wrong "
            "restart vouches for a state it does not describe; regenerate it "
            "with convert_restart.py --report.")
    report["__path__"] = str(candidate)
    return report


def conversion_surface_reason(report: dict, config: dict) -> str | None:
    """Why this converted state's surface does not match what the run staged.

    A converted restart's static surface records are the TEMPLATE's, and the
    template was cut from a run that staged a particular set of `.sra` files.
    If this run stages different ones, landmod will read the template's from
    the restart and the staged files will be silently discarded -- which is the
    same failure the donor-run guard below exists for, reached by a different
    route.
    """
    template = report.get("target_template") or {}
    hashes = template.get("surface_field_sha256")
    if not hashes:
        return ("the conversion report carries no surface hashes for its "
                "target template, so nothing can say which surface the "
                "converted state froze")
    now = intended_surface_codes(config)
    transparent = DUST_SURFACE_CODES
    was = {int(c) for c in hashes}
    if (was - transparent) != (now - transparent):
        return (f"the template staged codes {sorted(was)} and this run stages "
                f"{sorted(now)}")
    changed = [c for c in sorted(now)
               if surface_sra(config, c).is_file()
               and hashes.get(str(c)) != file_sha256(surface_sra(config, c))]
    if changed:
        return f"content changed for code(s) {changed} since the template"
    return None


def intended_surface_codes(config: dict) -> set[int]:
    """Which surface fields this run supplies rather than leaving at defaults.

    `model.soil_water_source` defaults to `uniform` when the key is absent, so a
    config predating this option behaves exactly as it did.

    Setting it to `pedology` makes 229 mandatory here, and 229 comes from a soil
    weathered under a climatology. On a terrain that has none yet, the bootstrap
    run has to go out with the key at `uniform`; see `docs/src/pipeline/sequencing.md` loop A, which
    also says why the flip cannot happen mid-run.
    """
    codes = set(BASE_SURFACE_CODES)
    if str(config["model"].get("land_albedo_source", "uniform")) != "uniform":
        codes |= ALBEDO_SURFACE_CODES
    if str(config["model"].get("soil_water_source", "uniform")) != "uniform":
        codes |= SOIL_WATER_SURFACE_CODES
    if str(config["model"].get("roughness_source", "uniform")) != "uniform":
        codes |= ROUGHNESS_SURFACE_CODES
    if str(config["model"].get("dust_source", "none")) != "none":
        codes |= DUST_SURFACE_CODES
    if str(config["model"].get("dust_emission", "none")) != "none":
        codes |= DUST_EMISSION_SURFACE_CODES
    return codes


def energy_diagnostics_enabled(config: dict) -> bool:
    """Whether to ask PlaSim for its 28-term energy decomposition.

    Defaults False when the key is absent, so a config predating this option
    behaves exactly as it did and `config_sha256` does not move. Turning it on
    adds 28 output fields and is meant for short diagnostic segments, not for
    production runs.
    """
    return bool(config["model"].get("energy_diagnostics", False))


def energy_diagnostics_level(config: dict) -> int:
    """NENERGY: 1 for the 28-term decomposition, 2 to add the conversion control.

    `energy_diagnostics: 2` asks the model for the split of its adiabatic
    conversion across the semi-implicit scheme, which is a CONTROL for world-0ov
    and costs three extra spectral transforms a timestep. `true` and `1` both
    mean 1, so a config predating this reads exactly as it did.
    """
    level = config["model"].get("energy_diagnostics", False)
    if isinstance(level, bool):
        return 1 if level else 0
    return int(level)


def set_low_io(model, low_io: bool) -> None:
    """Turn off PlaSim's low-I/O accumulation, so output is instantaneous.

    `NLOWIO = 1` is PlaSim's default (`plasimmod.f90:151`) and accumulates fields
    over each output interval, dividing in place at write time. `NLOWIO = 0`
    writes instantaneous records instead, 182 an orbit rather than 12 bins, and
    `pyburn` averages them to the same 12 for the .nc, so nothing downstream
    changes shape.

    WHY THIS IS THE DEFAULT. An accumulation cannot be undone and a sample set
    can always be averaged, so samples are strictly more information. The
    accumulation is also not the mean you would compute yourself: binned `spd`
    under `NLOWIO = 1` sits between the speed of the time-mean vector and the
    mean of instantaneous speeds. Anything reading variance, extremes or single
    records -- DUST-5's gust distribution above all -- needs the samples, and a
    spin-up reading scalars does not, which is what `--low-io` in
    `continue_exoplasim.py` is for. See `exoplasim/README.md`.

    A SEPARATE DEFECT, NOW FIXED. `NLOWIO = 1` used to corrupt the first output
    record of every model call, because `naccuout` survived a restart while the
    accumulators did not: bottom-level wind read 7.5x the other bins and
    humidity 27% low, while scalars stayed within 2%.
    `exoplasim-3.4.2-lowio-first-record.patch` repairs that and was verified on
    2026-08-18 against a binary with it reversed, which cut the artifact by about
    four orders of magnitude. It does NOT change what is accumulated, so it does
    not make the two regimes interchangeable and it is not why this function
    exists. Runs written before it still carry the defect; see
    `exoplasim/notes/first-output-bin.md`.

    WHAT IT COSTS. Raw output is several times larger, and with the postprocessor
    fixed the clean regime is about 1.27x an orbit rather than the 3.7x once on
    record -- almost all of that gap was a quadratic reader in `pyburn` and not
    the I/O mode. Model time is identical either way. See
    `notes/audits/pyburn-postprocessing-cost.md`. Delete run directories once
    their climatologies are extracted.

    It has to be reapplied on every continuation, because `configure()` rewrites
    the namelist each time.

    WRITTEN EXPLICITLY IN BOTH DIRECTIONS, never left to the compiled default.
    `nlowio` defaults to 1 in `plasimmod.f90` and is read on NROOT only; the
    broadcast that makes it agree across ranks is ours, and a namelist that
    states the regime is the artifact that says which one ran. Leaving it unsaid
    is what made the deadlock in `notes/audits/nlowio-collective-deadlock.md`
    invisible for a day.
    """
    model._edit_namelist("plasim_namelist", "NLOWIO", "1" if low_io else "0")


def declare_cold_start_seed(model, config: dict, is_cold: bool) -> None:
    """Give a cold start a DECLARED initial kick instead of the clock's.

    `initrandom` at `plasim.f90:2068` reads the namelist `SEED` when `seed(1)`
    is non-zero and calls `system_clock` otherwise. Nothing here wrote `SEED`,
    so every cold-started run drew a different initial condition and no cold
    result could be checked by re-running it. The resolution-ladder stability
    boundary was measured entirely on cold starts, which makes each of its cells
    one draw rather than a verdict.

    Only for a cold start. A run seeded from a restart reads its state, so the
    kick never happens and writing a seed there would suggest a control this has
    no part in.
    """
    if not is_cold:
        return
    seed = int(config["model"].get("cold_start_seed", 0))
    if seed == 0:
        raise RuntimeError(
            "model.cold_start_seed is 0 or absent, which is the value that makes "
            "initrandom fall back to the system clock. Declare it; a cold run "
            "nobody can reproduce is a measurement nobody can check.")
    model._edit_namelist("plasim_namelist", "SEED", str(seed))
    print(f"cold start: SEED = {seed} (declared; without it initrandom takes "
          f"the system clock and the run is unreproducible)")


def read_applied_energy_fix(run_dir: Path) -> dict | None:
    """What the energy fixer actually had to put back, from the run's own diag.

    The fixer restores the energy the adiabatic conversion loses (`world-0ov`),
    and by doing so it HIDES that loss: with it on, `denergy26 - denergy27`
    reports a residual near zero and the size of the defect is the correction
    instead. So the correction has to be readable, or a known 0.85 W/m2 quietly
    becomes an unknown one that can grow. `world-mzy`.

    Taken from the model's own print rather than recomputed, so this cannot
    disagree with what was applied. The first samples are dropped: the first
    step out of a restart shows an imbalance of order 250 W/m2 and the rate
    limiter spends its first steps climbing out of it, which is startup and not
    the defect.
    """
    diags = sorted(run_dir.glob("MOST_DIAG.*"))
    if not diags:
        return None
    # Anchored on the header rather than counted from the end, because the line
    # carries a varying number of trailing diagnostics: after the label come the
    # step and then the applied W/m2, and that order is what the model writes.
    values = []
    for line in diags[-1].read_text(errors="replace").splitlines():
        if "ENERGY FIXER applied" not in line:
            continue
        tail = line.split("K/day", 1)[-1].split()
        if len(tail) < 2:
            continue
        try:
            values.append(float(tail[1]))
        except ValueError:
            continue
    if len(values) < 8:
        return None
    settled = values[len(values) // 2:]
    return {"applied_w_m2_mean": sum(settled) / len(settled),
            "applied_w_m2_min": min(settled),
            "applied_w_m2_max": max(settled),
            "samples": len(settled),
            "note": "the energy the fixer had to restore; a CORRECTION and not "
                    "physics. It is the size of the defect on world-0ov, which "
                    "the fixer otherwise hides. Second half of the run only."}


def declare_energy_fixer(model, config: dict) -> bool:
    """Switch on the global energy fixer, which is a CORRECTION and not physics.

    The adiabatic step must conserve total energy: the column enthalpy it gives
    up has to equal the kinetic energy it takes on. It does not. Measured on a
    dry adiabatic run it loses about 0.8 W/m2, and the whole of that is the
    reference conversion's two halves being taken at different time levels of
    the semi-implicit scheme: the advective half is explicit at t, the
    divergence half is applied on `sdt`, the centred mean of t-dt and t+dt, and
    the displacement between them is -0.96 W/m2. That defect is `world-0ov`.
    This does not fix it.

    What this does is put the missing energy back as a uniform warming, the way
    ECHAM, the IFS and CAM all do. Without it the surface silently supplies the
    shortfall -- across three filter settings the adiabatic term moved 2.08 W/m2
    while the closed budget moved 0.013, with sensible heat flux up 2.17 and
    latent up 1.48, about a tenth of the sensible flux on this planet -- and
    downstream components read that as physics.

    IT MASKS WHAT IT COMPENSATES, so the applied increment is printed by the
    model at `ndiag` cadence and belongs on the manifest as a health metric: a
    change in the underlying defect has to be able to show rather than be
    absorbed. `world-mzy`.

    The fixer reads `denergy26` and `denergy27`, so it needs the energy
    diagnostics. Rather than switch them on quietly underneath the caller, which
    would make a declared setting mean something it does not say, this raises.
    """
    fix = config["model"].get("energy_fixer")
    if fix is None:
        raise RuntimeError(
            "model.energy_fixer is absent. It is a correction for a known "
            "defect (world-0ov) and whether a run carries it changes what its "
            "surface fluxes mean, so it is declared and never defaulted.")
    if not fix:
        model._edit_namelist("plasim_namelist", "NENERGYFIX", "0")
        print("energy fixer: OFF (declared). The core's conversion loses about "
              "0.8 W/m2 and the surface supplies it; see world-0ov.")
        return False
    if not config["model"].get("energy_diagnostics"):
        raise RuntimeError(
            "model.energy_fixer needs model.energy_diagnostics, because the "
            "fixer is driven by denergy26 and denergy27. Declare both rather "
            "than having one switch the other on underneath you.")
    model._edit_namelist("plasim_namelist", "NENERGYFIX", "1")
    print("energy fixer: ON (declared). A CORRECTION, not physics: it restores "
          "the energy the adiabatic conversion loses without restoring where it "
          "went. world-mzy; the defect is world-0ov.")
    return True


def declare_robert_filter(model, config: dict) -> float | None:
    """PNU, the leapfrog time filter's coefficient. world-0ov.

    The leapfrog scheme carries a computational mode that alternates sign every
    step, and the Robert-Asselin filter is what damps it. It matters here beyond
    stability: the adiabatic sink is the reference conversion's two halves being
    taken at different time levels, the displacement between them is the second
    time difference of the divergence, and a mode that alternates sign every
    step contributes to that difference whatever the timestep is.

    IT IS IN `planet_nl` AND NOT `plasim_nl`, which is where a reader looks
    first: `pnu` is declared in `plasimmod.f90` at 0.0 and set to 0.1 by
    `p_earth.f90`'s `planet_ini`, so the value that runs comes from the planet
    module and the key that overrides it sits in the planet namelist. Written to
    the wrong one the model aborts with "Cannot match namelist object name pnu",
    which is the loud failure and not the silent one.

    Absent means the value the planet module sets, so a config predating this
    reads exactly as it did.
    """
    value = config["model"].get("robert_filter")
    if value is None:
        return None
    model._edit_namelist("planet_namelist", "PNU", repr(float(value)))
    print(f"Robert time filter: PNU = {float(value)}")
    return float(value)


def declare_parmode(config: dict) -> str:
    """Which compiled parallel mode this run uses. DECLARED, never defaulted.

    `mpi` distributes NLAT over ranks and launches through mpiexec; `omp`
    distributes it over threads of one process, which is what the SHTns
    transform path requires because SHTns's parallelism is threads. They are
    different binaries built from a different flag line, and `nshtns` is 1 by
    default in the threaded build and 0 in the MPI one -- so the parmode decides
    which TRANSFORM integrates the run, not just how the work is spread.

    It is declared rather than defaulted because upstream's API knew only `mpi`
    and named its executable without a parmode, so a project that had moved to
    the threaded build could not ask for it and silently kept running the other
    one. That is what happened here: 179 runs, none of them threaded, while the
    threaded binaries were built, verified and gated on. world-bdh.
    """
    mode = config["model"].get("parmode")
    if mode is None:
        raise RuntimeError(
            "model.parmode is absent. It selects which compiled binary and "
            "therefore which spectral transform integrates the run, so it is "
            "declared and never defaulted. world-bdh.")
    if mode not in ("mpi", "omp"):
        raise RuntimeError(f"model.parmode must be 'mpi' or 'omp', not {mode!r}")
    return mode


def declare_dealias_conversion(model, config: dict) -> bool:
    """Truncate V.grad(ln ps) to the retained modes before the products. world-ly5.

    Hoskins and Simmons (1975) section 2, the paper `spectrala` cites: the
    transform grid "is insufficient for removing aliased interactions" for "the
    triple correlation involved in the energy conversion term and in the
    vertical advection terms", which need `M_g >= 4M + 1`. This model runs
    `NLON = 3*NTRU + 1`, the grid that dealiases a product of two band-limited
    fields and not of three, so the conversion is aliased by construction.

    ON, `calcgp` projects `zvgpg` onto the retained modes before anything
    multiplies it, which makes every product downstream quadratic in band
    limited fields again. That is HS75's own option (ii) and the cheaper of the
    two they name; the other is a `4M + 1` transform grid.

    IT CHANGES WHAT THE MODEL INTEGRATES, and it costs two extra transform pairs
    a timestep. It also has no SHTns path -- the projection uses the Legendre
    decomposition's partials -- and the model refuses rather than transforming
    through a path whose partials mean something else.

    Absent means off, so a config predating this reads exactly as it did.
    """
    if not config["model"].get("dealias_conversion", False):
        return False
    model._edit_namelist("plasim_namelist", "NDEALIAS", "1")
    print("dealiasing: V.grad(ln ps) is truncated to the retained modes before "
          "the conversion and the vertical advection use it. A CHANGE TO THE "
          "DYNAMICS; world-ly5.")
    return True


def declare_conversion_time_level(model, config: dict) -> bool:
    """Put the reference conversion's two halves on one time level. world-0ov.

    `calcgp` carries the advective half of the adiabatic reference conversion at
    time t and `spectrala` applies its divergence half on `sdt`, the centred
    mean of t-dt and t+dt. The identity that makes the two cancel the momentum
    equations' reference pressure-gradient work is the global integral of a
    mass-flux divergence, which holds at one time level and not across two, and
    the reference geopotential weights the split by up to 3.8 times `t0`. On the
    dry adiabatic arm the displacement is 0.96 W/m2 against a sink of 0.79.

    ON, this puts the divergence half back at time t. IT CHANGES WHAT THE MODEL
    INTEGRATES, and it takes that half out of the semi-implicit treatment in the
    temperature equation while the divergence solve still treats the temperature
    implicitly, so the timestep it is stable at is its own question and is NOT
    the one `model.resolution_timestep_minutes` answers.

    Absent means off, so a config predating this reads exactly as it did.
    """
    if not config["model"].get("conversion_time_level", False):
        return False
    model._edit_namelist("plasim_namelist", "NCONVTIME", "1")
    print("conversion time level: the reference conversion's divergence half is "
          "taken at t, not on sdt. A CHANGE TO THE DYNAMICS; world-0ov.")
    return True


def declare_hyperdiffusion(model, config: dict) -> dict:
    """Write the derived horizontal diffusion, overriding the compiled branch.

    ExoPlaSim hard-codes these for T21 and T42 only (`plasim.f90:1388`), and
    every finer rung silently falls through to T21's -- which at T170 is four
    times too weak at the truncation and, worse, damps at 41% of the LOCAL
    CASCADE RATE at half the truncation, because `nhdiff` is an absolute
    wavenumber worth 71% of T21's spectrum and 8.8% of T170's.

    Setting them here is enough and needs no source change: `prolog` calls
    `readnl` then `initpm`, `readnl` applies the branch and THEN reads the
    namelist, and `initpm` builds the operator from whatever survived. Rule 4
    does not apply.

    `nhdiff` is an integer wavenumber, so the fraction is rounded; the rounding
    is reported rather than hidden, because at coarse truncations one wavenumber
    is a several-percent change in where the damping starts.

    exoplasim/notes/resolution-tuned-parameters.md has the derivation.
    """
    model_cfg = config["model"]
    hd = model_cfg.get("hyperdiffusion")
    if not hd:
        raise RuntimeError(
            "model.hyperdiffusion is absent. Without it the run inherits T21's "
            "damping at every rung above T42, which is the defect that block "
            "exists to fix; see exoplasim/notes/resolution-tuned-parameters.md.")
    rung = str(model_cfg["resolution"]).upper()
    table = hd["timescales_days"]
    if rung not in table:
        raise RuntimeError(
            f"model.hyperdiffusion.timescales_days has no entry for {rung}. The "
            f"damping for a rung is derived from the rule, so a rung with no "
            f"entry has no derived damping and must not fall back to T21's.")
    tau = table[rung]
    ntru = int(rung.lstrip("Tt"))
    nhdiff = int(round(float(hd["cutoff_fraction"]) * ntru))
    # EVERY LEVEL, WRITTEN OUT. `TDISS*` and `NDEL` are per-level arrays, and a
    # Fortran namelist assignment of a scalar to an array sets element 1 only:
    # levels 2 to NLEV then keep whatever readnl's compiled branch left there,
    # which at T21 is ExoPlaSim's own table in days and at T42 is that branch's
    # values already in SECONDS. `dayseccheck` discriminates on maxval, so a
    # mixed array converts nothing and level 1 is read in the wrong unit. The
    # `n*value` replication is what ExoPlaSim's own wrapper writes; world-720.
    layers = int(model_cfg["layers"])
    keys = {"NDEL": f"{layers}*{int(hd['order_alpha'])}",
            "NHDIFF": f"{nhdiff}",
            "TDISSD": f"{layers}*{float(tau['divergence'])}",
            "TDISSZ": f"{layers}*{float(tau['vorticity'])}",
            "TDISST": f"{layers}*{float(tau['temperature'])}",
            "TDISSQ": f"{layers}*{float(tau['humidity'])}"}
    for key, value in keys.items():
        model._edit_namelist("plasim_namelist", key, value)
    print(f"hyperdiffusion: {rung} alpha={hd['order_alpha']} "
          f"nhdiff={nhdiff} (n*/N={nhdiff/ntru:.3f}, asked {hd['cutoff_fraction']}) "
          f"tau_vorticity={tau['vorticity']} d")
    return {"rung": rung, "nhdiff": nhdiff, "cutoff_fraction_actual": nhdiff / ntru,
            "order_alpha": int(hd["order_alpha"]), "timescales_days": dict(tau),
            "eddy_wind_m_s": float(hd["eddy_wind_m_s"])}


def declare_dry_constants(model, config: dict) -> dict:
    """The dry thermodynamic constants, the cold-start profile and the sponge.

    Four things ExoPlaSim otherwise takes from `p_earth.f90` and its compiled
    defaults, grouped because they are all read by the same dry column and all
    were inherited rather than chosen.

    AKAP. `p_earth.f90:41` declares R/cp as 0.286 and `readnl` derives the
    specific heat of dry air from it as `acpd = gascon/akap`. GASCON reaches the
    model derived from the declared composition and this did not, so the model
    ran a cp inconsistent with its own gas constant. Both come off
    `lib/lapse.py:gas_properties` here, which is the one place in this tree that
    turns the composition into (R, cp), and the R it returns is checked against
    the GASCON the model was actually configured with rather than assumed to
    agree. world-cwu.

    T0. The semi-implicit reference temperature, written as a per-level list
    because it is an array and a namelist scalar would set element 1 only.
    world-bmf.

    TGR, ALR, DTROP. The cold start's initial profile, and TGR additionally sets
    the orographic reduction of surface pressure on EVERY start. world-wmw.

    TFRC. Rayleigh drag on the top levels, in seconds, from a declaration in
    rotations. Written for every level so the value is the config's and not
    `readnl`'s unconditional NLEV==10 branch. world-aee.
    """
    import lapse

    model_cfg = config["model"]
    layers = int(model_cfg["layers"])
    gas_constant, cp = lapse.gas_properties(config)
    akap = gas_constant / cp
    # A CHECK THAT CAN FAIL: the same composition through two implementations.
    # If ExoPlaSim's mean molecular weight and this module's molar masses ever
    # disagree, akap and gascon stop being two views of one gas and the model
    # gets a cp that belongs to neither.
    configured = float(model.gascon)
    if abs(configured - gas_constant) > 0.01:
        raise RuntimeError(
            f"lib/lapse.py derives R = {gas_constant:.4f} J/kg/K from the "
            f"declared composition and ExoPlaSim was configured with GASCON = "
            f"{configured:.4f}. akap = R/cp is only meaningful beside the "
            "gascon it was derived with; reconcile the two before running.")

    profile = model_cfg["cold_start_profile"]
    t0 = float(model_cfg["semi_implicit_reference_temperature_k"])
    sponge = [float(x) for x in model_cfg["rayleigh_sponge_rotations"]]
    if len(sponge) != layers:
        raise ValueError(
            f"model.rayleigh_sponge_rotations has {len(sponge)} entries and the "
            f"model has {layers} levels. TFRC is per level and a short list is "
            "a level whose drag nobody declared.")
    # TFRC is entered in SECONDS. One rotation is the SIDEREAL day, which is
    # also the model's own unit of time after world-rt1, so a sponge declared in
    # rotations is the same number of nondimensional units at any rotation rate.
    rotation_s = float(config["planet"]["rotation_hours"]) * 3600.0
    tfrc = ",".join(f"{x * rotation_s:.6g}" for x in sponge)

    model._edit_namelist("planet_namelist", "AKAP", f"{akap:.8g}")
    model._edit_namelist("planet_namelist", "ALR",
                         f"{float(profile['lapse_rate_k_per_m']):.8g}")
    model._edit_namelist("plasim_namelist", "TGR",
                         f"{float(profile['surface_temperature_k']):.8g}")
    model._edit_namelist("plasim_namelist", "DTROP",
                         f"{float(profile['tropopause_height_m']):.8g}")
    model._edit_namelist("plasim_namelist", "T0", f"{layers}*{t0:.8g}")
    model._edit_namelist("plasim_namelist", "TFRC", tfrc)
    print(f"dry constants: akap={akap:.6f} (cp={cp:.2f} J/kg/K beside gascon "
          f"{configured:.4f}), t0={t0:g} K, cold start "
          f"{profile['surface_temperature_k']:g} K / "
          f"{float(profile['lapse_rate_k_per_m']) * 1000:g} K/km / "
          f"{profile['tropopause_height_m']:g} m, "
          f"sponge {sponge[:2]} rotations")
    return {"akap": akap, "cp_j_kg_k": cp, "gascon_j_kg_k": configured,
            "t0_k": t0, "cold_start_profile": dict(profile),
            "rayleigh_sponge_rotations": sponge,
            "rayleigh_sponge_seconds": [x * rotation_s for x in sponge]}


def declare_dynamics_only(model, config: dict) -> bool:
    """Strip every diabatic source, leaving the dynamical core alone.

    A DIAGNOSTIC CONFIGURATION, not a climate. With radiation, vertical
    diffusion and surface fluxes off and precipitation, convection and dry
    convective adjustment switched out of `rainstep`, nothing heats or cools the
    atmosphere. Anything that then appears in the energy budget is the numerics,
    and -- the point of it -- there is no surface exchange left to absorb it.

    That matters because the adiabatic residual is otherwise INVISIBLE: measured
    across three filters, `26 - 27` moved 2.08 W/m2 while the closed budget moved
    0.013, because the surface silently supplied whatever the dynamics destroyed.
    Take the surface away and the same sink has nowhere to hide.

    THE DAMPING STAYS ON, deliberately. `26 - 27` is measured across `spectrala`
    and the hyperdiffusion is applied in `spectrald`, so the damping does not
    enter the identity; removing it would only make the run blow up sooner. The
    physics filter is left to the config, so it can be varied as its own arm.
    """
    if not bool(config["model"].get("dynamics_only", False)):
        return False
    # RADIATION STAYS SWITCHED ON AND IS EMPTIED INSTEAD. `nrad = 0` skips
    # `radstep` entirely, which leaves `eccf` at its initialised zero, and
    # `outsc_` then writes `1.0/sqrt(eccf)` -- a divide by zero that traps under
    # the project's floating-point traps at the first output record. Setting
    # `nswr` and `nlwr` to zero leaves the orbital geometry computed and removes
    # every radiative tendency, which is what this mode actually wants.
    # NOT `nflux`. It is accepted by `plasim_nl` and read by nothing -- the only
    # occurrence in the source is its own namelist declaration -- so setting it
    # is silent and does nothing. Measured: with `nflux = 0` the latent flux was
    # still 6.5 W/m2. The live switches are fluxmod's.
    for key in ("NEVAP", "NSHFL", "NSTRESS", "NVDIFF"):
        model._edit_namelist("fluxmod_namelist", key, "0")
    for key in ("NSWR", "NLWR"):
        model._edit_namelist("radmod_namelist", key, "0")
    for key in ("NPRL", "NPRC", "NDCA", "NSHALLOW"):
        model._edit_namelist("rainmod_namelist", key, "0")
    print("dynamics only: nswr=0 nlwr=0, no evaporation, sensible flux, "
          "surface stress or vertical diffusion, and large-scale, convective, "
          "shallow and dry-adjustment heating all off. Not a climate.")
    return True


def enable_energy_diagnostics(model, config: dict) -> bool:
    """Set nenergy in plasim_nl, which the Python API does not expose.

    Same mechanism `stage_stellar_spectrum` uses for STARFILE: the namelist is
    edited directly because `configure()` has no parameter for it.
    """
    if not energy_diagnostics_enabled(config):
        return False
    model._edit_namelist("plasim_namelist", "NENERGY",
                         str(energy_diagnostics_level(config)))
    if config["model"].get("energy_diagnostics_3d", False):
        model._edit_namelist("plasim_namelist", "NENER3D", "1")
    return True


def enable_prescribed_dust(model, run_dir: Path, config: dict) -> dict | None:
    """Switch on the prescribed dust the patched radiation reads. DUST-11.

    Every value comes from the provenance file `build_surface_dust.py` wrote
    beside the field, so the run cannot be given a longwave ratio, a scale height
    or a burden that disagrees with the field it is applied to. There is no
    default for `DUSTQLW` here and none in the Fortran either: a shortwave-only
    dust is a larger error than no dust, so `radini` aborts rather than assume
    one, and this raises rather than supply one.

    `AEROFILE` goes in `radmod_namelist` rather than `aero_namelist`, because the
    prescribed path never runs `aero_ini` -- that is only called when the
    semi-Lagrangian tracer transport is on, and nothing here is transported. The
    file is staged into the run directory and named bare, for the same
    `character(len=80)` reason `stage_stellar_spectrum` does it.
    """
    if str(config["model"].get("dust_source", "none")) == "none":
        return None
    field = surface_sra(config, sorted(DUST_SURFACE_CODES)[0])
    prov_path = field.with_name(field.stem + "_provenance.json")
    if not prov_path.is_file():
        raise RuntimeError(
            f"{prov_path} is missing. Run exoplasim/scripts/build_surface_dust.py; "
            "the namelist values live with the field, not in the config.")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    values = prov["namelist_values"]
    if float(values.get("DUSTQLW", 0.0)) <= 0.0:
        raise RuntimeError(
            f"{prov_path} carries no thermal-infrared absorption ratio. "
            "Shortwave-only dust is worse than no dust; see notes/dust.md.")
    aerofile = PROJECT_ROOT / prov["aerofile"]
    if not aerofile.is_file():
        raise RuntimeError(
            f"{aerofile} is missing. Run exoplasim/scripts/dust_aerofile.py.")
    shutil.copyfile(aerofile, run_dir / aerofile.name)
    model._edit_namelist("radmod_namelist", "AEROFILE", f"'{aerofile.name}'")
    scale = float(config["model"].get("dust_scale", values["DUSTSC"]))
    for key, value in (("NDUSTRAD", "1"),
                       ("DUSTSC", f"{scale}"),
                       ("DUSTHSC", f"{float(values['DUSTHSC'])}"),
                       ("DUSTQLW", f"{float(values['DUSTQLW'])}")):
        model._edit_namelist("radmod_namelist", key, value)
    return {"variant": prov["variant"], "dustsc": scale,
            "dusthsc": float(values["DUSTHSC"]),
            "dustqlw": float(values["DUSTQLW"]),
            "aerofile": aerofile.name,
            "field_sha256": prov["output_sha256"],
            "land_mean_band1_optical_depth":
                prov["field_statistics"]["land_mean_band1_optical_depth"]}


def enable_dust_emission(model, run_dir: Path, config: dict) -> dict | None:
    """Write the `aero_nl` group the in-model emission scheme reads. DUST-3.

    This is the half of the item that is easiest to skip and it is load-bearing:
    the emission patch and the deposition patch both add namelist keys with NO
    Fortran defaults, so until something writes this group none of them is
    reachable and `aero_ini` aborts rather than run. `aero_namelist` is written
    by ExoPlaSim's own Python API and this project had never touched it.

    Every value comes from the provenance file
    `aeolian/scripts/build_dust_source_fields.py` wrote beside the fields, so a
    run cannot be given a threshold, a Weibull shape or a roughness that
    disagrees with the source map it is applied to. Same arrangement as
    `enable_prescribed_dust` and for the same reason.

    THE REMOVAL SWITCHES RIDE HERE TOO. `ldepvel` and `lwetdep` come from
    exoplasim-3.4.2-aerosol-deposition.patch, they default off, and their
    coefficients live in the same `aeolian/config/dust.yaml`. They are read
    directly from that config rather than from the field provenance because they
    describe the atmosphere rather than the surface, and nothing about them
    depends on which terrain the fields were built on.

    Turning emission on also turns on the aerosol machinery it lives in:
    `L_AERO` gates `aero_main` and `l_source = 2` selects the dust case. Every
    run this project has made so far set `L_AERO = 0`, which is why `aerocore`
    had never executed and why seven defects survived in it.
    """
    if str(config["model"].get("dust_emission", "none")) == "none":
        return None
    field = surface_sra(config, min(DUST_EMISSION_SURFACE_CODES))
    prov_path = field.with_name(field.stem + "_provenance.json")
    if not prov_path.is_file():
        raise RuntimeError(
            f"{prov_path} is missing. Run "
            "aeolian/scripts/build_dust_source_fields.py; the namelist values "
            "live with the fields, not in the config.")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    values = dict(prov["namelist_values"])

    dust_cfg_path = PROJECT_ROOT / "aeolian" / "config" / "dust.yaml"
    dust_cfg = yaml.safe_load(dust_cfg_path.read_text(encoding="utf-8"))
    removal = dust_cfg["removal"]

    # The removal terms are opt-in per run, because each one is a separate A/B
    # arm and A3 wants them switchable from the namelist alone.
    if bool(config["model"].get("dust_dry_deposition", False)):
        vd = config["model"].get("dust_deposition_velocity_m_s")
        if vd is None:
            raise RuntimeError(
                "model.dust_dry_deposition is on but "
                "model.dust_deposition_velocity_m_s is unset. vdaero is the "
                "NON-gravitational part of dry deposition and there is no "
                "defensible default; see the removal block of "
                "aeolian/config/dust.yaml and the deposition patch header.")
        values["LDEPVEL"] = 1
        values["VDAERO"] = float(vd)
    if bool(config["model"].get("dust_wet_scavenging", False)):
        values["LWETDEP"] = 1
        values["SCAVA"] = float(removal["scavenging_a"])
        values["SCAVB"] = float(removal["scavenging_b"])

    model._edit_namelist("plasim_namelist", "L_AERO", "1")
    model._edit_namelist("aero_namelist", "l_source", "2")

    # THE SHIPPED aero_namelist CANNOT BE READ, and nothing has noticed because
    # `aero_ini` runs only when `L_AERO > 0` and every run this project has made
    # set it to 0. It carries `aerofile = 0`, an unquoted integer for a
    # `character(len=80)`, which gfortran rejects with iostat 5010 -- and
    # `aero_ini` reads without an iostat, so the first run to turn the aerosol
    # on aborts in the namelist read before any of this is reached. Quoting it
    # is the fix and it has to happen here, because this is the first thing that
    # sets L_AERO at all.
    aerofile = PROJECT_ROOT / "exoplasim" / "data" / "dust" / "vesper_dust_aerosol.dat"
    if not aerofile.is_file():
        raise RuntimeError(
            f"{aerofile} is missing. Run exoplasim/scripts/dust_aerofile.py.")
    shutil.copyfile(aerofile, run_dir / aerofile.name)
    model._edit_namelist("aero_namelist", "aerofile", f"'{aerofile.name}'")

    # RADIATIVELY INERT, and that is not a preference. `l_aerorad = 1` would put
    # the emitted dust into the shortwave through `radmod`'s own `apart`, which
    # `aero_ini` never populates from the namelist -- upstream defect 1 in
    # `aeolian/notes/in-model-dust.md`, still open because it lives in
    # `radmod.f90`. At this world's effective radius that path gives an optical
    # depth 1/385 of intent, and the longwave term is missing on top of it,
    # which `notes/dust.md` prices at several kelvin of spurious cooling. So the
    # only non-wrong setting available on this branch is off. The file is staged
    # and named anyway, so turning it on once item 5 lands is a namelist change.
    model._edit_namelist("aero_namelist", "l_aerorad", "0")

    for key, value in values.items():
        model._edit_namelist("aero_namelist", key, f"{value}")
    return {
        "z0_bracket_end": prov["z0_bracket_end"],
        "aeolian_z0_m": prov["aeolian_z0_m"],
        "terrain_hash": prov["terrain_hash"],
        "field_sha256": prov["output_sha256"],
        "namelist_values": values,
        "source_cells": prov["field_statistics"]["srcw_cells_nonzero"],
    }


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
            if code in SOIL_WATER_SURFACE_CODES:
                builder, setting = ("build_surface_soil_water.py",
                                    "model.soil_water_source")
            elif code in DUST_SURFACE_CODES:
                builder, setting = ("build_surface_dust.py", "model.dust_source")
            elif code in DUST_EMISSION_SURFACE_CODES:
                builder, setting = ("build_dust_source_fields.py",
                                    "model.dust_emission")
            else:
                builder, setting = ("build_surface_albedo.py",
                                    "model.land_albedo_source")
            fallback = ("none"
                        if code in DUST_SURFACE_CODES | DUST_EMISSION_SURFACE_CODES
                        else "uniform")
            raise RuntimeError(
                f"{src} is missing. Run {builder}, or set {setting}: {fallback} "
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
        "parmode": str(m.get("parmode", "mpi")),
        "precision_bytes": int(m["precision_bytes"]),
        "flux_ratio": round(float(flux_ratio), 6),
        "co2_ppm": round(1e6 * float(a["pCO2_bar"]), 3),
        "ch4_ppmv": round(1e6 * float(a.get("pCH4_bar", 0.0)), 4),
        "n2o_ppmv": round(1e6 * float(a.get("pN2O_bar", 0.0)), 4),
        "rotation_hours": float(p["rotation_hours"]),
        "obliquity_degrees": float(p["obliquity_degrees"]),
        "eccentricity": float(p["eccentricity"]),
        "glaciers": bool(config["surface"].get("glaciers", {}).get("enabled")),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "geography": geography_tag(config),
    }


def run_id(config: dict, flux_ratio: float) -> str:
    """A UUID. Not derived from anything.

    This used to spell out the physical parameters in the directory name, so
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


def expected_namelist_keys(config: dict) -> dict:
    """What the CONFIG says the run's namelists must contain, {file: {KEY: value}}.

    Derived from `config/planet.yaml` and from nothing else -- deliberately not
    from the staging code, because a check built out of the staging lists tests
    only that they agree with themselves. CONS-9.
    """
    m = config["model"]
    want: dict = {"radmod_namelist": {}, "icemod_namelist": {}, "plasim_namelist": {},
                  "planet_namelist": {}}
    for key, name, default in SHORTWAVE_GAS_KEYS:
        v = m.get(key)
        if v is not None and float(v) != default:
            want["radmod_namelist"][name] = float(v)
    o3 = m.get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        want["radmod_namelist"]["O3SCALE"] = float(o3)
    salinity = config.get("ocean", {}).get("salinity_psu")
    if salinity is not None:
        want["icemod_namelist"]["TFREEZE"] = round(freezing_point_k(salinity), 4)
    if energy_diagnostics_enabled(config):
        want["plasim_namelist"]["NENERGY"] = float(energy_diagnostics_level(config))
    if m.get("conversion_time_level", False):
        want["plasim_namelist"]["NCONVTIME"] = 1.0
    if m.get("dealias_conversion", False):
        want["plasim_namelist"]["NDEALIAS"] = 1.0
    # THE FIXER IS A DECLARED SWITCH AND BELONGS IN THE CHECK. It is written on
    # both branches -- 1 when on and 0 when off -- so its absence from a staged
    # namelist means a continuation dropped it, which is exactly what this check
    # exists to catch. failure-modes class 22.
    if m.get("energy_fixer") is not None:
        want["plasim_namelist"]["NENERGYFIX"] = 1.0 if m["energy_fixer"] else 0.0
    if m.get("robert_filter") is not None:
        want["planet_namelist"]["PNU"] = float(m["robert_filter"])
        if m.get("energy_diagnostics_3d", False):
            want["plasim_namelist"]["NENER3D"] = 1.0
    # The scalar half of `declare_dry_constants`. Each of these reverts to a
    # p_earth.f90 or plasimmod.f90 Earth default if a continuation drops it, and
    # TGR reverting moves the realised mean surface pressure, so the class-22
    # failure here is a change in the mean state rather than a lost switch. T0
    # and TFRC are per-level lists and cannot go through this float comparison.
    profile = m.get("cold_start_profile")
    if profile is not None:
        import lapse
        gas_constant, cp = lapse.gas_properties(config)
        # Formatted the way `declare_dry_constants` writes it, so the comparison
        # is against the value the namelist can hold and not against one more
        # digit than it carries.
        want["planet_namelist"]["AKAP"] = float(f"{gas_constant / cp:.8g}")
        want["planet_namelist"]["ALR"] = float(profile["lapse_rate_k_per_m"])
        want["plasim_namelist"]["TGR"] = float(profile["surface_temperature_k"])
        want["plasim_namelist"]["DTROP"] = float(profile["tropopause_height_m"])
    return {f: keys for f, keys in want.items() if keys}


def verify_staged_namelists(run_dir: Path, config: dict) -> dict:
    """Every config key that departs from the model default reached the namelist.

    CONS-9. Reads the RUN, which is the only artifact that records what was
    actually integrated -- `docs/src/practice/failure-modes.md` class 12 said so and class 22
    is what happens when nobody checks: `h2o_sw_level` was declared at 1.127,
    staged by the prepare script, dropped by every continuation, and no guard
    noticed for a whole run. `config_drift` compared the config against the
    manifest and both said 1.127; this compares the config against the FILE.

    Raises rather than warns. A run whose radiation does not match its own
    configuration is not a cheaper run, it is a different world.
    """
    checked, wrong = {}, []
    for fname, keys in expected_namelist_keys(config).items():
        path = run_dir / fname
        for key, want in keys.items():
            try:
                got = float(namelist_value(path, key))
            except (KeyError, FileNotFoundError):
                wrong.append(f"{key} absent from {fname}, config wants {want:g}")
                continue
            except ValueError:
                wrong.append(f"{key} in {fname} is not a number")
                continue
            if abs(got - want) > 1e-9:
                wrong.append(f"{key} in {fname} is {got:g}, config wants {want:g}")
            checked[f"{key}@{fname}"] = got
    if wrong:
        raise SystemExit(
            "the staged namelists do not match config/planet.yaml:\n  "
            + "\n  ".join(wrong)
            + "\nThe namelist in the run directory is what the model integrates. "
              "Fix the staging rather than the check; see docs/src/practice/failure-modes.md "
              "class 22.")
    return checked


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
    # THE DEFAULT IS THE CHEAP REGIME, because the orbits this script writes are
    # a spin-up: a run is prepared cold or off a seed and integrates toward
    # equilibrium, and nothing should be building a climatology from the first
    # block. The clean regime is bought per segment, at the end, for the orbits
    # an analysis actually reads -- `continue_exoplasim.py --purpose
    # post_equilibrium_climatology`, which refuses low I/O outright.
    #
    # This was `disable_low_io()` with no way to say otherwise, so every orbit
    # of every run paid for instantaneous samples and the spin-up ones threw
    # them away. Safe to flip because the bin-0 defect is fixed and verified
    # (`exoplasim/notes/first-output-bin.md`) and because `nlowio` is now
    # broadcast, without which NLOWIO = 0 deadlocks a patched model
    # (`notes/audits/nlowio-collective-deadlock.md`).
    #
    # What low I/O costs is not accuracy but INFORMATION: it writes interval
    # accumulations where the clean regime writes instantaneous samples, and an
    # accumulation cannot be undone. Anything reading variance, extremes or
    # single records needs --clean-io.
    # BENCHMARK OVERRIDES. Both change the machine's layout and neither changes
    # the physics, so they belong on the command line rather than in
    # config/planet.yaml -- a config edit would move config_sha256 and make
    # every existing run unresumable, for a question about core placement.
    # `exoplasim/notes/rank-layout-benchmark.md` is what they are for.
    parser.add_argument(
        "--ncpus", type=int, default=None,
        help="override model.ncpus for this run. Selects a different compiled "
             "binary, since ExoPlaSim builds one per (resolution, layers, "
             "ranks, parmode); NLAT must divide by it")
    parser.add_argument(
        "--mpi-opts", type=str, default=None,
        help="extra flags for mpiexec, e.g. "
             "'--map-by pe-list=0,1,2,3,4,5,6,7:ordered --bind-to core'. The "
             ":ordered qualifier is load-bearing -- without it the ranks share "
             "one pool instead of getting a core each")
    # DIAGNOSTIC. The model writes an output record every `nafter` timesteps and
    # defaults to one per EARTH day, 32 steps of 45 minutes. This planet's day is
    # 30 hours, which is exactly 40 steps, so the default samples the diurnal
    # cycle at a 32/40 = 4/5 ratio: consecutive records advance the phase by
    # four fifths of a day and land on FIVE PHASES, forever. The regular output's
    # annual mean therefore carries a fixed aliasing bias, while the ocean and
    # ice streams accumulate every step and do not. CLIM-11.
    #
    # Choose a value coprime with the 40-step day to spread the sampling over
    # all forty phases -- 37 or 39 -- and the bias should collapse if that is
    # what it is. Diagnostic only: it changes what is written, not what is
    # integrated, but it moves the record count and so the climatology.
    parser.add_argument(
        "--writes-per-day", type=int, default=None,
        help="NWPD, output records per day. The model derives its write "
             "interval as nafter = mtspd / nwpd, so this is the control and "
             "NAFTER is not a namelist key -- setting that one aborts the run "
             "with `Cannot match namelist object name nafter`. Raising it "
             "samples more diurnal phases")
    parser.add_argument(
        "--clean-io", dest="low_io", action="store_false", default=True,
        help="run this block at NLOWIO = 0, writing instantaneous samples "
             "rather than interval accumulations. The default is the cheap "
             "regime, because a prepared run is a spin-up; pass this when the "
             "block itself has to be read as data")
    parser.add_argument(
        "--run-years", type=int, default=0,
        help="Run this many model orbits after preparation",
    )
    parser.add_argument(
        "--force-prepare", action="store_true",
        help="Allow re-preparing an existing run with no climate outputs",
    )
    parser.add_argument(
        "--superseded-surface-ok", action="store_true",
        help="allow --restart-from when the staged surface fields have moved "
             "since the donor run. The new fields are DISCARDED: a restart "
             "reads albedo, roughness and soil water from itself. Only valid "
             "for a paired A/B, where every arm inherits the same superseded "
             "surface and the difference is what is measured. Stamped on the "
             "run manifest so the run is self-labelling.")
    parser.add_argument(
        "--restart-from", type=Path, default=None,
        help="Seed the initial state from an existing MOST_REST file instead of "
             "cold-starting. Only the spin-up path changes, not the equilibrium.",
    )
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    # Applied to the loaded config rather than threaded through, so the binary
    # name, the physical fingerprint and the model construction cannot disagree
    # about how many ranks this run has. `config_sha256` hashes the FILE, so
    # provenance still records the config as written, and the fingerprint
    # records the ranks that actually ran.
    if args.ncpus is not None:
        if args.ncpus < 1:
            raise ValueError("--ncpus must be positive")
        config["model"]["ncpus"] = int(args.ncpus)
    mpi_opts = args.mpi_opts
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
    superseded_surface = None
    conversion = None
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
        # CLIM-31. A restart carries the donor's PARTIAL ACCUMULATION and the
        # counters that normalise it, so a seeded run opens mid-window: it adds
        # its own steps to somebody else's sum and divides by its own count. The
        # damage is one output record and it is visible on a field that cannot
        # vary -- under low I/O the land mask reads 0.95341 in orbit 0 bin 0,
        # exactly 1.0 everywhere else, and a cold-started low-I/O run reads 1.0
        # throughout.
        #
        # Zeroed in a COPY, so the donor run is never modified: a restart is the
        # only record of where a run was, and editing one in place would make a
        # completed run unreproducible to fix a defect in the run seeded from
        # it. Every record it rewrites is one the model itself puts back to a
        # clean value at an interval boundary -- which is zero for most and a
        # large sentinel for the two running minima -- so the copy describes a
        # run sitting exactly at the start of an accumulation window.
        seeded_from = restart_seed

        # BEFORE the copy, and before the run directory exists. A guard that
        # fires after doing work leaves the work behind: this one used to run
        # after the seed was written and left a run directory holding nothing
        # but a MOST_REST.seed on every refusal.
        #
        # A CONVERTED restart is not a donor run's file and has no run
        # directory to carry a manifest. Its provenance is the conversion
        # report `convert_restart.py` wrote beside it, which names the template
        # its static surface records came from and the hash of every staged
        # surface file that template was cut against. That is exactly what the
        # guard below needs, so a converted state is checked rather than
        # exempted.
        conversion = load_conversion_report(seeded_from, restart_seed)
        src_manifest = seeded_from.parent / "run_manifest.json"
        if conversion is not None:
            reason = conversion_surface_reason(conversion, config)
            if reason and args.superseded_surface_ok:
                print(f"  OVERRIDDEN (--superseded-surface-ok): {reason}")
                superseded_surface = reason
                reason = None
            if reason:
                raise RuntimeError(
                    "--restart-from refused a converted restart: " + reason
                    + ". The converted state's static surface records come "
                    "from its target template, not from this run's staged "
                    ".sra files, so a mismatch means the run would integrate "
                    "a surface it did not stage. Rebuild the template against "
                    "the current surface with build_restart_template.py.")
            print(f"  conversion report: "
                  f"{conversion['source']['truncation']} -> "
                  f"{conversion['target_template']['truncation']}, donor "
                  f"{(conversion.get('source_manifest') or {}).get('run_id')}, "
                  "target surface matches this run's staged files")
        elif not src_manifest.is_file():
            raise RuntimeError(
                f"--restart-from refused: {seeded_from} has neither a "
                "run_manifest.json beside it nor a conversion report, so "
                "nothing says which surface it froze. A restart FREEZES soil "
                "water, roughness and albedo -- landmod reads them from the "
                "restart, not from the .sra -- and with no provenance that "
                "cannot be ruled out. Seed from a run directory, or convert "
                "with convert_restart.py, which writes the report this needs.")
        if src_manifest.is_file():
            src = json.loads(src_manifest.read_text(encoding="utf-8"))
            was = set((src.get("surface_fields") or {}).get("from_file") or [])
            now = intended_surface_codes(config)
            reason = None
            # The guard above is about fields landmod reads only on a cold
            # start. It does not apply to every surface field, and refusing a
            # restart for one it does not apply to costs a cold start for no
            # reason. `radini` reads the dust field with `mpsurfgp` outside any
            # `nrestart` test, where landmod's block is inside `if (nrestart ==
            # 0)`, so a restarted run picks it up. Adding or removing dust is
            # therefore a forcing change and an ordinary perturbation
            # experiment, which is exactly what seeding from an equilibrium is
            # for. The exemption is stated per code rather than assumed for the
            # class: anything else added here must be checked the same way,
            # in the Fortran and not from the name.
            transparent = DUST_SURFACE_CODES
            if (was - transparent) != (now - transparent):
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
            if reason and args.superseded_surface_ok:
                # Deliberately overridden. The guard protects the CANONICAL
                # chain: a spin-up seeded this way would silently discard the
                # newer surface and reproduce its parent. An A3 forcing A/B is
                # the one case where that does not matter, because every arm
                # inherits the SAME superseded surface and what is being
                # measured is the DIFFERENCE between arms, not the absolute
                # climate. The override is a flag rather than a deletion so it
                # cannot be taken by accident, and it is stamped on the manifest
                # below so the resulting run says what it is.
                print(f"  OVERRIDDEN (--superseded-surface-ok): {reason}. "
                      "The staged surface is DISCARDED and this run inherits "
                      "its parent's. Valid for a paired A/B, invalid for "
                      "anything whose absolute climate is read.")
                superseded_surface = reason
                reason = None
            if reason:
                raise RuntimeError(
                    "--restart-from refused: " + reason + ". A restart reads soil "
                    "water, roughness and albedo from ITSELF, not from the .sra "
                    "files -- landmod's landini takes dwmax, dz0clim and all "
                    "three dalbcl bands from the restart when restart > 0 -- so "
                    "those changes would be silently discarded and the run would "
                    "reproduce its parent. Cold-start instead.")

        # Every refusal above has now been passed, so the copy is made and the
        # run directory created. Zeroed in a COPY so the donor is never
        # modified: a restart is the only record of where a run was, and
        # editing one in place would make a completed run unreproducible to fix
        # a defect in the run seeded from it.
        restart_seed = run_dir / "MOST_REST.seed"
        run_dir.mkdir(parents=True, exist_ok=True)
        reset_info = reset_restart_accumulators.reset(seeded_from, restart_seed)
        print(f"seed {seeded_from.name} from {seeded_from.parent.name}: "
              f"zeroed {len(reset_info['zeroed'])} accumulator records, "
              f"{len(reset_info['sentinels'])} put back to a nonzero clean "
              f"value, {len(reset_info['never_reset_by_the_model'])} the model "
              "resets nowhere left as they were (CLIM-31)")

    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    parmode = declare_parmode(config)
    if parmode == "omp" and mpi_opts:
        raise RuntimeError(
            "--mpi-opts was given but model.parmode is 'omp', which launches one "
            "process and no mpiexec. The flags would be silently dropped.")
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
        workdir=str(run_dir),
        modelname=identifier,
        outputtype=model_cfg["output_type"],
        hyperthreading=False,
        mpi_opts=mpi_opts,
        parmode=parmode,
    )
    print(f"parmode: {parmode} "
          f"({'threads, SHTns transform' if parmode == 'omp' else 'ranks, legmod transform'})")
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
        filterkappa=float(model_cfg["filter_kappa"]),
        filterpower=int(model_cfg["filter_power"]),
        landmap=str(landmap),
        topomap=str(topomap),
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=int(derived["snapshot_interval_steps"]),
        otherargs=configure_otherargs(derived),
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
    verify_stellar_spectrum(model, config)
    # Ozone column scaling, set through the namelist because the Python API does
    # not expose o3scale. The model prescribes an Earth column and derives
    # nothing about it from the host star; Segura et al. (2003) measure 0.794 of
    # Earth's for a K2V host at 1 PAL O2.
    o3 = config["model"].get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        model._edit_namelist("radmod_namelist", "O3SCALE", f"{float(o3)}")
        print(f"ozone column scaled to {float(o3)} of Earth's (Segura et al. 2003)")
    # SHORTWAVE_GAS_KEYS is the one list; see its definition for why it is not
    # written out here.
    for key, name, default in SHORTWAVE_GAS_KEYS:
        w = config["model"].get(key)
        if w is not None and float(w) != default:
            model._edit_namelist("radmod_namelist", name, f"{float(w)}")
            print(f"{name} = {float(w)} (radmod.f90 default {default})")

    for name, ppmv in trace_gas_ppmv(config).items():
        model._edit_namelist("radmod_namelist", name, f"{ppmv:.6g}")
        print(f"{name} = {ppmv:.6g} ppmv (radmod.f90 default 0.0, meaning absent)")

    dust = enable_prescribed_dust(model, run_dir, config)
    if dust is not None:
        print(f"prescribed dust ON: {dust['variant']} burden, land-mean band-1 "
              f"optical depth {dust['land_mean_band1_optical_depth']:.4f} x "
              f"{dust['dustsc']}, DUSTQLW {dust['dustqlw']:.5f}. Needs "
              "patches/exoplasim-3.4.2-prescribed-dust.patch and a rebuild.")

    emission = enable_dust_emission(model, run_dir, config)
    if emission is not None:
        print(f"in-model dust emission ON: {emission['source_cells']} source "
              f"cells, aeolian z0 {emission['aeolian_z0_m']:g} m "
              f"({emission['z0_bracket_end']} end). Needs "
              "patches/exoplasim-3.4.2-dust-emission.patch, the two aerosol "
              "patches under it, and a rebuild.")

    if args.writes_per_day is not None:
        if args.writes_per_day < 1:
            raise ValueError("--writes-per-day must be positive")
        model._edit_namelist("plasim_namelist", "NWPD",
                             str(int(args.writes_per_day)))
        print(f"NWPD = {args.writes_per_day} writes/day "
              f"(default 1 samples 5 diurnal phases; see CLIM-11)")
    declare_cold_start_seed(model, config, args.restart_from is None)
    hyperdiffusion = declare_hyperdiffusion(model, config)
    dry_constants = declare_dry_constants(model, config)
    dynamics_only = declare_dynamics_only(model, config)
    robert_filter = declare_robert_filter(model, config)
    conversion_time_level = declare_conversion_time_level(model, config)
    dealias_conversion = declare_dealias_conversion(model, config)
    set_low_io(model, args.low_io)
    energy_fixer = declare_energy_fixer(model, config)
    if enable_energy_diagnostics(model, config):
        n = register_energy_diagnostic_codes()
        print(f"energy diagnostics on: nenergy=1, {n} codes 360-387 registered "
              "with the postprocessor. Term 15 needs "
              "patches/exoplasim-3.4.2-energy-diagnostics.patch and a rebuild.")
    model.exportcfg(str(run_dir / f"{identifier}.cfg"))
    surface_report = surface_field_report(run_dir, config)

    # CONS-9. AFTER every staging call and BEFORE the run, because the namelist
    # on disk is the only artifact that records what will actually be
    # integrated. Raises, so a run whose radiation does not match its own
    # configuration never starts.
    staged_namelists = verify_staged_namelists(run_dir, config)
    print(f"namelists verified {len(staged_namelists)} config-set keys present "
          f"with the declared values: {', '.join(sorted(staged_namelists))}")

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
    exe_suffix = "_omp" if str(model_cfg.get("parmode")) == "omp" else ""
    exe_path = run_dir / (
        f"most_plasim_t{int(str(model_cfg['resolution']).lstrip('Tt'))}"
        f"_l{int(model_cfg['layers'])}_p{int(model_cfg['ncpus'])}{exe_suffix}.x"
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
        # The DONOR is the provenance, not the copy this run reads. Recording
        # the copy would name a file that exists only inside this run directory
        # and whose sha nothing else can match, which is the opposite of what a
        # provenance field is for.
        "initial_state": ({"cold_start": True} if restart_seed is None else {
            "cold_start": False,
            "superseded_surface_override": superseded_surface,
            "restart_from": str(seeded_from),
            "restart_from_sha256": file_sha256(seeded_from),
            "restart_from_run": seeded_from.parent.name,
            "accumulators_reset": [n for n, _ in reset_info["zeroed"]],
            # A converted state BEGINS a lineage. The donor is provenance, not
            # continuation: resolution conversion changes the represented state
            # on purpose and the target model owns the equilibrium.
            "conversion": None if conversion is None else {
                "report": str(conversion["__path__"]),
                "converter_version": conversion.get("converter_version"),
                "from": conversion["source"]["truncation"],
                "to": conversion["target_template"]["truncation"],
                "real_bytes": [conversion["source"]["real_bytes"],
                               conversion["target_template"]["real_bytes"]],
                "donor_run": (conversion.get("source_manifest") or {}).get("run_id"),
                "template": conversion["target_template"].get("path"),
                "template_cut_from_run":
                    conversion["target_template"].get("cut_from_run"),
                "records_the_model_must_rebuild":
                    conversion.get("expected_to_change_in_model_fixup"),
                "whole_run_accumulators_restarted":
                    conversion.get("whole_run_accumulators_restarted"),
            },
            "seed_copy_sha256": file_sha256(restart_seed),
            "note": "Initial condition only; equilibrium is set by the forcing. "
                    "The accumulator records listed were zeroed in a copy before "
                    "the run read it, so this run does not inherit the donor's "
                    "partial accumulation window (CLIM-31); every other record "
                    "is byte-identical to the donor.",
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
            "path": str(exe_path),
            "sha256": file_sha256(exe_path),
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
        "namelist_checks": {**checks, **staged_namelists},
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
        # The spectrum by CONTENT, because the line above is a filename and the
        # file is regenerated in place. `continue_exoplasim.py` compares this on
        # every resume; without it a regenerated `k25v.dat` changes the star
        # mid-run and nothing says so. CONS-3.
        "stellar_spectrum_digest": stellar_spectrum_digest(config),
        # Null when the run has no prescribed dust, which is most of them. The
        # values are copied from the field's own provenance rather than the
        # config, so a run says what burden and what longwave ratio it actually
        # carried instead of what a setting asked for.
        "prescribed_dust": dust,
        # What damping this run actually got, including the rounded nhdiff, so a
        # result can be read against the operator that produced it rather than
        # against the config that was meant to.
        "hyperdiffusion": hyperdiffusion,
        # The dry thermodynamic constants, the cold-start profile and the top
        # sponge, all of which used to be ExoPlaSim's Earth defaults and none of
        # which appeared in a manifest.
        "dry_constants": dry_constants,
        "dynamics_only": dynamics_only,
        "conversion_time_level": conversion_time_level,
        "dealias_conversion": dealias_conversion,
        "robert_filter": robert_filter,
        # A CORRECTION and not physics; what it is correcting is world-0ov and
        # the fixer itself is world-mzy. On the manifest because whether a run
        # carries it changes what its surface fluxes mean.
        "energy_fixer": energy_fixer,
        # Null when the run has no interactive emission, which is most of them.
        # Like prescribed_dust, the values are copied from the fields' own
        # provenance rather than from the config, so the manifest says what the
        # run carried and not what a setting asked for.
        "dust_emission": emission,
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
            if energy_fixer:
                applied = read_applied_energy_fix(run_dir)
                if applied is None:
                    raise RuntimeError(
                        "the energy fixer was declared on and the run's diag "
                        "carries no readable correction. The fixer hides the "
                        "defect it compensates, so a run that cannot report "
                        "what it applied is a run whose fluxes nobody can "
                        "check. world-mzy.")
                manifest["energy_fixer_applied"] = applied
                print(f"energy fixer applied {applied['applied_w_m2_mean']:+.4f} "
                      f"W/m2 on average over the second half "
                      f"({applied['applied_w_m2_min']:+.3f} to "
                      f"{applied['applied_w_m2_max']:+.3f}); that is the size of "
                      f"the defect on world-0ov, not a physical flux")
            # REGISTER THE BLOCK AS A SEGMENT, like a continuation does.
            # Without this the orbits this script writes carry no purpose and
            # no `low_io`, and CLIM-9's two defaults then disagree about them:
            # unlabelled counts as production for convergence and as TAINTED
            # for I/O regime. That was harmless while the block was always
            # NLOWIO = 0 and merely conservative; with the cheap regime as the
            # default it would be a run whose first block really is
            # accumulations and says so nowhere.
            #
            # The purpose is `spinup` and is not a parameter: this script
            # prepares a run and integrates it from a cold start or a seed,
            # which is the definition of one. Orbits meant to be read as data
            # come from `continue_exoplasim.py --purpose
            # post_equilibrium_climatology`, on a run something has assessed.
            manifest.setdefault("segments", []).append({
                "start_year_index": 0,
                "end_year_index": args.run_years - 1,
                "seasonal_output": True,
                "low_io": bool(args.low_io),
                "high_cadence": False,
                "purpose": "spinup",
                # WHY A FIRST RECORD IS REFUSABLE, and it is no longer CLIM-31.
                #
                # CLIM-31 was the donor's partial accumulation window arriving
                # in the seed: under low I/O the land mask, which cannot vary,
                # read 0.95341 in the first record. That is FIXED. The
                # accumulator reset takes its record set and its clean values
                # from `restart_schema.py` rather than from a hand list that
                # covered 49 of 114, and a converted seeded run now reads
                # exactly 1.0 there.
                #
                # What remains is narrower and applies only to a CONVERTED
                # state. Its derived surface records -- albedo, roughness --
                # come from the target template and are not consistent with the
                # prognostics that were remapped onto it, so the model's first
                # timestep runs on the template's albedo before rebuilding it.
                # One timestep, and then it is gone: measured at T21 to T42
                # against a template 0.104 out in area-weighted albedo, the
                # first record carries 0.213 W/m2 of reflected shortwave, which
                # is 33.6 W/m2 spread over the 160 steps the record averages.
                # Under clean I/O the first record is an instantaneous sample
                # rather than an accumulation, so it does not pick this up.
                # world-eyb.
                "first_record_tainted": bool(args.low_io and conversion is not None),
                "stellar_spectrum_digest": stellar_spectrum_digest(config),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            })
        except Exception:
            manifest["status"] = "failed"
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            raise
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
