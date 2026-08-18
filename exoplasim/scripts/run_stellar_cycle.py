#!/usr/bin/env python3
"""Run or resume a superposed-sinusoid stellar-flux experiment with ExoPlaSim.

    python exoplasim/scripts/run_stellar_cycle.py \\
        --baseline run_xxxxxxxxxxxx --restart MOST_REST.00059

The forcing is two sinusoids about `orbit.baseline_flux_earth`, defined by
`stellar_cycle.components` in `config/planet.yaml` and applied at every
radiation timestep by the patched `radmod.f90`. Nothing about the forcing is
set here.

The cycle is bolometric and grey. That is a measured decision rather than a
simplification of convenience: a spot-driven swing of this size moves the
Lacis-Hansen band-1 fraction by 0.007, which is worth under 0.02 K on a world
this nearly ice-free. It would not be negligible on a much colder branch --
see exoplasim/notes/parameter-decisions.md for the calculation and the
condition under which it stops holding.
"""

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
    stage_surface_extras,
    surface_field_report,
    stage_stellar_spectrum,
    stellar_spectrum_path,
    verify_stellar_spectrum,
    physical_fingerprint,
    run_id,
    REGULAR_CODES,
    derive,
    file_sha256,
)
from continue_exoplasim import validate_year, year_diagnostics  # noqa: E402


JULIAN_YEAR_DAYS = 365.25

# radmod.f90 carries exactly two component slots after the dual-sinusoid patch:
# the unsuffixed triple and a `2`-suffixed one. Naming them here rather than
# accepting whatever the config lists means a third component is a loud
# KeyError instead of a silently dropped forcing.
COMPONENT_SLOTS = {"medium": "", "long": "2"}


def load_components(config: dict) -> dict:
    """The cycle's components, from config rather than hardcoded.

    The amplitudes were literals in this file while the baseline flux was 0.90.
    The baseline has moved twice since, so a literal here would have applied an
    old range about a new centre without anything noticing.

    Two superposed sinusoids, not one. The medium component is climatic and the
    long one is geomorphic, and their periods are deliberately non-commensurate
    so grand minima differ in depth rather than repeating -- see the block's own
    notes in `config/planet.yaml`.
    """
    block = config.get("stellar_cycle", {}) or {}
    components = block.get("components")
    if not components:
        raise RuntimeError("config/planet.yaml has no stellar_cycle.components")
    unknown = set(components) - set(COMPONENT_SLOTS)
    if unknown:
        raise RuntimeError(
            f"stellar_cycle.components has {sorted(unknown)}, but the patched "
            f"radmod carries only {sorted(COMPONENT_SLOTS)}. Adding a component "
            "needs a namelist slot in "
            "exoplasim/patches/exoplasim-3.4.2-star-cycle.patch first.")
    for name, spec in components.items():
        for key in ("period_earth_years", "amplitude_flux_peak_to_peak"):
            if spec.get(key) is None:
                raise RuntimeError(f"stellar_cycle.components.{name} has no {key}")
        if float(spec["period_earth_years"]) <= 0:
            raise RuntimeError(f"component {name} has a non-positive period")
    # The declared total is a cross-check on the parts, not an input. It exists
    # because prose elsewhere quotes it, and prose goes stale.
    declared = block.get("total_amplitude_flux_peak_to_peak")
    actual = sum(float(s["amplitude_flux_peak_to_peak"]) for s in components.values())
    if declared is not None and not math.isclose(declared, actual, rel_tol=1e-9):
        raise RuntimeError(
            f"stellar_cycle.total_amplitude_flux_peak_to_peak is {declared} but "
            f"the components sum to {actual}")
    return components


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


def completed_years(run_dir: Path, manifest: dict) -> list[int]:
    """Which orbits are finished, from the files, checked against the manifest.

    Enumerating the whole set and requiring it to be contiguous is not the same
    thing as choosing an artifact by sort order: nothing here is selected, and a
    gap is an error rather than a silently shorter run. The manifest's own count
    is then required to agree, so a run interrupted between writing an output and
    recording it cannot resume from a different place than it reports.
    """
    years = sorted(int(m.group(1)) for path in run_dir.glob("MOST.*.nc")
                   if (m := re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)))
    if years and years != list(range(years[-1] + 1)):
        raise RuntimeError(f"Outputs are non-contiguous: {years}")
    recorded = manifest.get("completed_orbits")
    if recorded is not None and recorded != len(years):
        raise RuntimeError(
            f"manifest records {recorded} completed orbits but {len(years)} are "
            f"on disk in {run_dir}. Resolve by hand: one of them is wrong, and "
            "guessing which would resume the cycle at the wrong phase.")
    return years


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
            **cycle["namelist"],
        },
    )
    # The cycle scales `gsol0` and leaves the band split alone, so the split has
    # to be the right one to begin with. Building the model on an existing run
    # directory re-copies the shipped namelists, so the spectrum is staged again
    # here; without it `solarini` falls back to a blackbody at `STARBBTEMP` and
    # the cycle would be centred on a differently-coloured star than the
    # baseline it is a variance about.
    stage_stellar_spectrum(model, run_dir, stellar_spectrum_path(config))
    verify_stellar_spectrum(model, config)
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
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--baseline", type=str, default=None,
                        help="run id or directory to centre the cycle on")
    parser.add_argument("--restart", type=str, default=None,
                        help="restart filename within the baseline run")
    parser.add_argument(
        "--cycles", type=float, default=2.0,
        help="Run at least this many periods of the LONGEST component",
    )
    parser.add_argument(
        "--orbits", type=int,
        help="Override the target output-orbit count (useful for validation)",
    )
    parser.add_argument("--run-id", type=str, default=None,
                        help="resume this existing run directory")
    args = parser.parse_args()
    if args.cycles <= 0:
        raise ValueError("Number of cycles must be positive")

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    components = load_components(config)
    # The cycle is variance about the baseline. There is no separate mean to
    # configure and no opportunity for the two to disagree, which is what the
    # superseded min/max case representation allowed.
    mean_ratio = float(config["orbit"]["baseline_flux_earth"])
    derived = derive(config, mean_ratio)
    identifier = args.run_id or run_id(config, mean_ratio)
    run_dir = (RUNS / identifier).resolve()
    if args.run_id and not (run_dir / "run_manifest.json").is_file():
        raise SystemExit(f"--run-id {identifier} has no run_manifest.json")
    run_dir.mkdir(parents=True, exist_ok=True)
    source_dir = (INPUTS / "exoplasim_cycle_t42").resolve()
    model_cfg = config["model"]
    # Named from the config, not hardcoded. This said p8 while
    # build_star_cycle_exoplasim.sh had moved to p16, so the cycle run would have
    # died on a missing executable -- or worse, found a stale p8 left over.
    executable = source_dir / (
        f"most_plasim_t{int(str(model_cfg['resolution']).lstrip('Tt'))}"
        f"_l{int(model_cfg['layers'])}_p{int(model_cfg['ncpus'])}.x")
    patch_file = (PATCHES / "exoplasim-3.4.2-star-cycle.patch").resolve()
    if not executable.is_file():
        raise RuntimeError(
            f"{executable} is not built. Run "
            "exoplasim/scripts/build_star_cycle_exoplasim.sh first.")

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
    earth_constant = float(config["orbit"]["earth_solar_constant_w_m2"])
    steps_per_orbit = derived["runsteps_per_orbit"]

    resolved, namelist = {}, {}
    for name, spec in components.items():
        suffix = COMPONENT_SLOTS[name]
        period_years = float(spec["period_earth_years"])
        period_steps = period_years * JULIAN_YEAR_DAYS * 86400.0 / timestep_seconds
        semi_ratio = 0.5 * float(spec["amplitude_flux_peak_to_peak"])
        phase = float(spec.get("phase_cycles", 0.0))
        resolved[name] = {
            "period_earth_years": period_years,
            "period_local_orbits": period_steps / steps_per_orbit,
            "period_model_steps": period_steps,
            "amplitude_flux_peak_to_peak": float(spec["amplitude_flux_peak_to_peak"]),
            "semi_amplitude_flux_earth": semi_ratio,
            "semi_amplitude_w_m2": semi_ratio * earth_constant,
            "phase_cycles": phase,
            "namelist_slot": f"GSOL*{suffix or '(unsuffixed)'}",
        }
        namelist[f"GSOLAMP{suffix}@radmod_namelist"] = repr(semi_ratio * earth_constant)
        namelist[f"GSOLPERIOD{suffix}@radmod_namelist"] = repr(period_steps)
        namelist[f"GSOLPHASE{suffix}@radmod_namelist"] = repr(phase)

    longest = max(resolved.values(), key=lambda c: c["period_earth_years"])
    total_ptp = sum(c["amplitude_flux_peak_to_peak"] for c in resolved.values())
    cycle = {
        "components": resolved,
        "mean_flux_earth": mean_ratio,
        "mean_flux_w_m2": mean_ratio * earth_constant,
        # The envelope, reached only when the components align. With
        # non-commensurate periods that is approached rather than attained, which
        # is the design: see stellar_cycle in config/planet.yaml.
        "aligned_minimum_flux_earth": mean_ratio - 0.5 * total_ptp,
        "aligned_maximum_flux_earth": mean_ratio + 0.5 * total_ptp,
        "total_amplitude_flux_peak_to_peak": total_ptp,
        "start_model_step": start_step,
        "phase_convention": "each sinusoid begins at mean flux and rises",
        "stellar_spectrum_treatment": (
            "fixed 4965 K blackbody spectral partition; the cycle is a grey "
            "multiplier. A spot-driven 6% bolometric swing shifts the band-1 "
            "fraction by 0.007, worth under 0.02 K at this ice cover -- see "
            "exoplasim/notes/parameter-decisions.md."),
        "namelist": namelist,
    }
    target_orbits = (
        args.orbits
        if args.orbits is not None
        else math.ceil(args.cycles * longest["period_local_orbits"])
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
            # What this run IS. The id is a UUID and carries no meaning, so
            # index_runs.py reads this; without it a cycle run indexes as
            # legacy and its forcing is invisible to everything downstream.
            # `forcing` is what distinguishes it from a static run at the
            # same flux, which is otherwise an identical fingerprint.
            "physical": {**physical_fingerprint(config, mean_ratio),
                         "forcing": "stellar_cycle",
                         "cycle_components": {
                             n: [c["period_earth_years"],
                                 c["amplitude_flux_peak_to_peak"]]
                             for n, c in resolved.items()}},
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

    years = completed_years(run_dir, manifest)
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
            # Phase per component, not one number. With two superposed
            # sinusoids a single phase is meaningless, and the whole point
            # of the design is that the two disagree.
            for edge, offset in (("start", 0), ("end", 1)):
                diagnostic[f"cycle_{edge}_phase"] = {
                    n: float(((year + offset) * steps_per_orbit
                              / c["period_model_steps"] + c["phase_cycles"]) % 1.0)
                    for n, c in resolved.items()}
            diagnostic["forcing_flux_earth"] = {
                edge: mean_ratio + sum(
                    c["semi_amplitude_flux_earth"]
                    * math.sin(2 * math.pi * diagnostic[f"cycle_{edge}_phase"][n])
                    for n, c in resolved.items())
                for edge in ("start", "end")}
            manifest["diagnostics"] = [
                d for d in manifest["diagnostics"] if d["year_index"] != year
            ] + [diagnostic]
            manifest["completed_orbits"] = year + 1
            manifest["completed_cycles"] = {
                n: (year + 1) * steps_per_orbit / c["period_model_steps"]
                for n, c in resolved.items()}
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
