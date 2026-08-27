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
from provenance import applied_removals, config_drift  # noqa: E402
from restart_surface import verify_restart_surface_fields  # noqa: E402
from segments import SEGMENT_PURPOSES  # noqa: E402
from run_exoplasim import (  # noqa: E402
    declare_dry_constants,
    declare_hyperdiffusion,
    prepare_thread_stack,
    declare_dynamics_only,
    declare_hyperdiffusion,
    declare_energy_fixer,
    declare_robert_filter,
    declare_conversion_time_level,
    declare_dealias_conversion,
    declare_storm_diagnostics,
    declare_ecological_stream,
    refuse_eco_codes,
    SHORTWAVE_GAS_KEYS,
    trace_gas_ppmv,
    verify_staged_namelists,
    configure_otherargs,
    surface_sra,
    stage_surface_extras,
    surface_field_report,
    stage_stellar_spectrum,
    stellar_spectrum_path,
    stellar_spectrum_digest,
    intended_surface_codes,
    require_stellar_spectrum,
    verify_stellar_spectrum,
    REGULAR_CODES,
    ENERGY_DIAGNOSTIC_CODES,
    ENERGY_3D_CODES,
    set_low_io,
    enable_energy_diagnostics,
    energy_diagnostics_enabled,
    register_energy_diagnostic_codes,
    HIGH_CADENCE_CODES,
    SNAPSHOT_CODES,
    derive,
    file_sha256,
    restore_run_executable,
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



# Config keys that no script passes to the model. They are declarations for
# readers, not inputs, so a change to one cannot alter a run and must not block
# resuming it.
#
# Deliberately a short, explicit allowlist rather than a rule. This guard already
# has a history: it once compared raw file bytes, so an edited comment blocked a
# legitimate resume, which was fixed by comparing parsed values. This is the same
# failure one level up -- a semantically real change to a physically inert key.
# The fix is to name the inert keys, not to loosen the comparison. Anything not
# listed here is assumed to reach the model.
#
# A key that has been REMOVED from the configuration is a THIRD case and does not
# belong here: an inert entry excuses a key that still exists and whose value
# may change again, while a removal is a one-time event with a value that can be
# recorded. `lib/provenance.py:REMOVED_CONFIG_KEYS` is where it is declared, and
# the guard checks both halves of that declaration before it excuses anything.
# world-51wj.
INERT_CONFIG_KEYS = {
    # PROSE CARRIED INSIDE THE YAML. A resume reads nothing out of a `note`;
    # `run_stellar_cycle.py` takes only `period_earth_years` and
    # `amplitude_flux_peak_to_peak` out of `stellar_cycle.components`. Both keys
    # were already traced and declared inert for the BIOSPHERE consumer, whose
    # entry says this check "used to fail on" a comment edit -- and the lesson
    # was applied to one consumer only, so editing the same two notes blocked a
    # resume today. `inert` is per consumer by design, so the repair is to trace
    # them here too rather than to share a list.
    "stellar_cycle.components.medium.note",
    "stellar_cycle.components.long.note",
    "star.spectral_type",     # a label; the model gets effective_temperature_k
                              # and the spectrum file, not this. THE ENTRY IS
                              # CORRECT ABOUT THE KEY AND WRONG ABOUT THE
                              # ARTIFACT: nothing passes the spectral type to
                              # the model, but `build_stellar_spectrum.py`
                              # writes `k25v.dat` from it, in place and under
                              # the same name, so the value reaches the
                              # radiation through a file the key only names.
                              # A config comparison cannot see that, and
                              # tightening this list would not help -- the
                              # spectrum can also be regenerated with no config
                              # edit at all. `require_stellar_spectrum` below
                              # compares the FILE, which is the route the value
                              # actually takes. CONS-3.
    "star.surface_uv_relative_to_earth",
                              # a design declaration; ExoPlaSim models no
                              # ultraviolet and no script reads this. Spelled in
                              # full because it was spelled `star.surface_uv`
                              # here, which is not a key in the config and so
                              # excused nothing; `provenance.unknown_inert_keys`
                              # is the check that now says so.
    "schema_version",         # bookkeeping
    # Which climatology DOWNSTREAM components read. BOTH keys, because there are
    # two climatologies: the bootstrap is the run on terrain-only surface fields
    # and the baseline is the run on those fields. Nothing on the run or resume
    # path touches either: only the surface-field builders do, and what they
    # produce is guarded by `surface_field_report`'s presence check on every
    # resume; content is not compared (the per-code sha comparison runs only on
    # the --restart-from prepare path). Leaving either here blocks a resume for
    # the entirely expected act of naming the climatology the run itself
    # produced, which is a false positive that trains people to reach for a
    # bypass.
    "baseline_climatology",
    "bootstrap_climatology",
    # The cartographic declaration. `lib/provenance.py` carries the trace and
    # the argument; repeated here because this guard is a deliberate explicit
    # allowlist rather than an import of someone else's, and because the run
    # and resume path is the one consumer that would strand a run in flight.
    # Nothing passes a spin direction or a zero meridian to the model.
    "planet.rotation_direction",
    "planet.longitude_positive",
    "planet.prime_meridian",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--flux-ratio", type=float, default=None)
    parser.add_argument("--run", type=str, default=None,
                        help="run id or directory to continue (required)")
    parser.add_argument("--orbits", type=int, default=5)
    # WHAT THESE ORBITS ARE FOR, declared rather than inferred. It used to be
    # inferred from `--seasonal-output` plus the run's equilibrium cutoff, which
    # labelled a three-orbit low-I/O verification segment as climatology input;
    # it had to be corrected by hand, and by then `build_climatology.py` was the
    # only thing standing between it and a production climatology. The flags say
    # what the model was asked to WRITE. Only the caller knows what the orbits
    # are FOR, and every rule over the flags gets that wrong in a new way as
    # soon as a flag is added. Required, and there is no default, because a
    # default is the inference again with fewer places to notice it.
    parser.add_argument(
        "--purpose", required=True, choices=SEGMENT_PURPOSES,
        help="what this segment is for. `spinup` integrates toward equilibrium; "
             "`post_equilibrium_climatology` is orbits meant to be read as this "
             "world's climate, and needs a run already assessed; `diagnostic` "
             "measures the model rather than the planet -- an I/O verification, "
             "a high-cadence wind sample, a block on a differently patched "
             "binary -- and is kept out of convergence windows and "
             "climatologies")
    # Snapshots carry orbital phase, which the regular output does not, and
    # analyze_climatology cannot map output bins onto orbital position without
    # them. Defaulting them OFF meant a 60-orbit baseline finished with no
    # snapshots and needed a 12-orbit top-up; that had happened before. They
    # cost one extra file per bin and the run is the expensive part, so the
    # default is now ON and the flag turns them off.
    parser.add_argument(
        "--no-seasonal-output", dest="seasonal_output",
        action="store_false", default=True,
        help="Skip instantaneous seasonal snapshots for this segment. They are "
             "written by default because analyze_climatology needs the orbital "
             "phase they carry and a run without them has to be extended.",
    )
    # DUST-5 wants a gust distribution rather than a climatology: 32 seasonal
    # snapshots resolve synoptic variance and not the diurnal and sub-daily
    # variance that actually lifts dust, so a Weibull fitted to them is too
    # narrow and the emission it implies is a lower bound. Sampling every fourth
    # timestep gives 1462 samples a cell over one orbit, which resolves a
    # 30-hour day. The postprocessed field list is trimmed to the winds by
    # HIGH_CADENCE_CODES; note that this does NOT shrink what the model writes.
    # The raw `MOST_HC.NNNNN` is 15 GB for one orbit at T42 whatever is asked
    # for, because the model's own output path does not take a field list. It is
    # deleted once pyburn has run, but the disk has to be there first.
    # A spin-up orbit and a climatology orbit want different things, and the
    # reason is what the two regimes MEAN rather than what they cost. Low I/O
    # writes the model's own accumulation over each output interval; clean I/O
    # writes instantaneous records that pyburn averages into the same twelve
    # bins. An accumulation cannot be undone and a sample set can always be
    # averaged, so anything reading variance, extremes or single records needs
    # the clean regime, while a spin-up reading scalars does not. The
    # accumulation is also not the mean you would compute: binned `spd` under
    # low I/O sits between the speed of the time-mean vector and the mean of
    # instantaneous speeds.
    #
    # The corrupt first record that used to ride along with low I/O is FIXED and
    # verified, see exoplasim/notes/first-output-bin.md, so it is no longer a
    # reason to avoid the cheap regime for spin-up.
    #
    # Model time is identical either way. The wall-clock gap was 104 s against
    # 387 s an orbit when measured 2026-08-17; almost all of that was a
    # quadratic reader in pyburn rather than the I/O mode, and with it fixed the
    # clean regime costs about 1.27x. See
    # notes/audits/pyburn-postprocessing-cost.md. Segments record which regime
    # they ran, and `build_climatology.py` refuses to mix them.
    io_mode = parser.add_mutually_exclusive_group()
    io_mode.add_argument(
        "--low-io", dest="low_io", action="store_true", default=None,
        help="force PlaSim's low-I/O accumulation ON for this segment: cheaper "
             "per orbit, and every orbit carries interval ACCUMULATIONS rather "
             "than instantaneous samples, which cannot be undone afterwards")
    io_mode.add_argument(
        "--clean-io", dest="low_io", action="store_false", default=None,
        help="force it OFF, writing instantaneous samples. Needed by anything "
             "reading variance, extremes or single records")
    # EFOR-2, and it is per SEGMENT for the same reason the I/O regime is:
    # `configure()` rewrites the namelist on every continuation, so a segment
    # that does not restate it does not write the stream. The stream cannot
    # change a result either way.
    parser.add_argument(
        "--ecological-stream", action="store_true",
        help="NECO = 1 for this segment, writing the biosphere's own output "
             "stream into MOST_ECO.NNNNN. Off by default; the restart is "
             "byte-identical either way")
    parser.add_argument(
        "--eco-interval-steps", type=int, default=None,
        help="NECOSTEP for this segment. Left unset the model uses mtspd, one "
             "absolute 24-hour day exactly")
    parser.add_argument(
        "--high-cadence", action="store_true",
        help="write near-surface wind every fourth timestep for this segment, "
             "into highcadence/MOST_HC.NNNNN.nc. For DUST-5")
    parser.add_argument(
        "--high-cadence-interval", type=int, default=4,
        help="timesteps between high-cadence samples (default 4)")
    # The prepare-time flag of the same name authorises ADOPTING a donor's
    # surface. It says nothing about the segments that follow, and
    # `stage_surface_extras` rewrites the current `.sra` into the run directory
    # on every resume, so such a run comes to assert a surface the model will
    # not read. Restating it here keeps that a claim about this segment, made on
    # purpose and stamped on it. See restart_surface.py.
    parser.add_argument(
        "--superseded-surface-ok", action="store_true",
        help="this run adopted its donor's surface at prepare time and the "
             "staged .sra no longer describes what the model reads. Continue "
             "anyway, and stamp the segment with what it is actually "
             "integrating. Valid for a paired A/B, never for the canonical "
             "chain")
    args = parser.parse_args()
    if args.orbits < 1:
        raise ValueError("--orbits must be positive")
    if args.high_cadence and args.high_cadence_interval < 1:
        raise ValueError("--high-cadence-interval must be positive")
    # THE I/O REGIME FOLLOWS THE PURPOSE, which is the whole reason the purpose
    # is declared rather than inferred. A spin-up integrates toward equilibrium
    # and nothing reads its output, so it gets the cheap regime; the other two
    # exist to be READ, so they get instantaneous samples. Either flag overrides.
    if args.low_io is None:
        args.low_io = args.purpose == "spinup"
    # These two combinations are not judgement calls. The reason for the first
    # has CHANGED and the refusal has not: it used to be that a low-I/O orbit
    # carried a corrupt first output record, and that defect is fixed and
    # verified. What remains is not a bug but an inequality -- low I/O writes
    # interval ACCUMULATIONS where the clean regime writes instantaneous
    # samples, a mean can be recovered from samples and an accumulation cannot
    # be undone -- so a climatology built on it is permanently poorer than one
    # that was not. An orbit with no seasonal snapshots carries no orbital
    # phase, which is the second.
    if args.purpose == "post_equilibrium_climatology":
        if args.low_io:
            raise SystemExit(
                "--purpose post_equilibrium_climatology with --low-io: low I/O "
                "writes interval accumulations, not instantaneous samples, and "
                "that cannot be undone afterwards -- variance, extremes and "
                "single records are gone from those orbits for good. Drop "
                "--low-io, or declare the segment --purpose diagnostic.")
        if not args.seasonal_output:
            raise SystemExit(
                "--purpose post_equilibrium_climatology with --no-seasonal-output: "
                "analyze_climatology needs the orbital phase only the snapshots "
                "carry, so such a segment has to be extended before it can "
                "produce a climatology. Use --purpose spinup, or drop "
                "--no-seasonal-output.")

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    # A run id is a UUID and cannot be recomputed, so a continuation must be told
    # which run to continue. That is the safer direction: the old behaviour
    # rebuilt the name from config and could silently resolve to a different run
    # than the one intended -- which it did, when a patch changed the physics
    # without moving anything the name encoded.
    if args.run is None:
        raise SystemExit(
            "--run is required: pass the run id or its directory. "
            "`python exoplasim/scripts/index_runs.py` lists what is on disk.")
    cand = Path(args.run)
    run_dir = (cand if cand.is_dir() else RUNS / args.run).resolve()
    identifier = run_dir.name
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"No prepared run manifest at {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    # Compare the parsed configuration against the copy the manifest already
    # stores, not the file's bytes. Hashing the raw file makes an edited comment
    # indistinguishable from an edited parameter, which blocks a legitimate
    # resume and says nothing about why. This reports the offending keys.
    # THE THREAD COUNT IS ADOPTED FROM THE MANIFEST, NOT COMPARED WITH IT.
    # `run_exoplasim.py --ncpus N` writes N into the loaded config before the
    # manifest is stamped, so the manifest carries N while the config FILE
    # still declares whatever it declares, and the comparison below then
    # refuses a run it should extend -- with no flag able to say otherwise,
    # which left such a run unresumable at all. It cost a re-run of a paired
    # soil-thermal experiment. world-q4gh.
    #
    # ADOPTION AND NOT AN EXEMPTION, because the count is not free to differ.
    # A continuation has to integrate on the same executable as the segment
    # before it -- the run directory holds exactly one, and ExoPlaSim compiles
    # one per (resolution, layers, threads) -- so a config declaring 16 against
    # a manifest of 8 must resolve the p8 binary. Merely dropping the key from
    # the comparison would let it resolve a p16 binary into a run whose earlier
    # orbits are p8, which is worse than the refusal it replaces. The value is
    # taken from the manifest INTO the config so that everything downstream,
    # the binary resolver included, reads one number.
    manifest_ncpus = (manifest.get("source_config", {})
                      .get("model", {}).get("ncpus"))
    if manifest_ncpus is not None and \
            config["model"].get("ncpus") != manifest_ncpus:
        print(f"  adopting ncpus = {manifest_ncpus} from the run manifest, "
              f"against config's {config['model'].get('ncpus')}: a "
              f"continuation integrates on the executable the run already "
              f"holds. world-q4gh.")
        config["model"]["ncpus"] = manifest_ncpus

    drift = config_drift(manifest["source_config"], config,
                         INERT_CONFIG_KEYS)
    if drift:
        raise RuntimeError(
            "Configuration differs from the run manifest; refusing to resume:\n  "
            + "\n  ".join(drift)
        )
    # A key REMOVED from the configuration reads to the comparison above exactly
    # as an edited parameter does, so it refuses both unless the removal was
    # declared in advance with the value the key held. That refusal killed a
    # commissioning at 28 orbits, and a declaration that nobody sees at the
    # moment it is used is a silence however well it is written down; say which
    # ones this resume rested on. lib/provenance.py:REMOVED_CONFIG_KEYS.
    for line in applied_removals(manifest["source_config"], config):
        print(f"  resuming across a declared removal: {line}")

    # The other half of that comparison, and the half a parsed-value diff cannot
    # do: the spectrum reaches the model as a FILE, and the file is regenerated
    # in place under a name that never moves. Checked here, before anything
    # expensive, because a guard that fires after the orbits are integrated has
    # cost exactly what it exists to save. CONS-3.
    if require_stellar_spectrum(manifest, config):
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    # The flux comes from the RUN, not from the config, and that is the whole
    # point of reading it here rather than above. A run prepared with
    # --flux-ratio carries a flux the config does not: `source_config` still
    # holds the baseline, so config_drift sees nothing to complain about, and a
    # continuation that took the config's value silently integrated a bracket
    # point at the baseline flux. That is exactly what happened to the first 0.91
    # point, which spent 45 orbits at 0.945 while its manifest said 0.91 --
    # detected only because two runs 0.035 apart converged to the same
    # temperature, which is not physical.
    manifest_flux = float(manifest["physical"]["flux_ratio"])
    if args.flux_ratio is not None and abs(args.flux_ratio - manifest_flux) > 1e-9:
        raise SystemExit(
            f"--flux-ratio {args.flux_ratio:g} does not match this run's "
            f"{manifest_flux:g}. A continuation cannot change the physics of the "
            "run it continues; start a new run instead."
        )
    flux_ratio = manifest_flux
    derived = derive(config, flux_ratio)

    years = output_years(run_dir)
    if not years or years != list(range(years[-1] + 1)):
        raise RuntimeError(f"Run outputs are absent or non-contiguous: {years}")
    start_year = years[-1] + 1
    restart = run_dir / f"MOST_REST.{years[-1]:05d}"
    if not restart.is_file():
        raise RuntimeError(f"Missing restart file {restart}")

    # The one part of the old inference worth keeping, now as a CHECK on the
    # declaration rather than a substitute for it. `assess_convergence.py` writes
    # `equilibrium_cutoff_year_index`, so its absence means the run has never
    # been assessed, and calling orbits post-equilibrium on a run nothing has
    # judged is a claim the manifest can refuse.
    cutoff = manifest.get("equilibrium_cutoff_year_index")
    # AN ARM RUN HAS NO LINEAGE TO CONTINUE INTO. The prepare stamps
    # `canonical_lineage_eligible` false when the caller named the binary, which
    # means the executable is one build_model.py refused to publish -- another
    # precision, another flag line, a patched model source. The climatology a
    # post-equilibrium segment exists to produce is what a canonical lineage is
    # made of, so it is the one purpose an arm may not be given; a spin-up or a
    # diagnostic segment on an arm is exactly what an arm is for. world-u5pf.
    if (args.purpose == "post_equilibrium_climatology"
            and manifest.get("canonical_lineage_eligible") is False):
        arm = (manifest.get("executable") or {}).get("arm") or {}
        raise SystemExit(
            f"{identifier} was integrated by the arm {arm.get('build_tag')}, "
            f"which is not the registry's binary and is absent from "
            f"binary_manifest.json by construction. A "
            f"post_equilibrium_climatology segment produces a climatology, and "
            f"a climatology from an unpublished arm cannot belong to the "
            f"canonical lineage. Use --purpose spinup or --purpose diagnostic, "
            f"and compare the arm against production climatologically rather "
            f"than by adopting its output.")
    if args.purpose == "post_equilibrium_climatology":
        if cutoff is None:
            raise SystemExit(
                f"{identifier} has no equilibrium_cutoff_year_index, so it has "
                "never been assessed and nothing may treat it as settled. Run "
                "assess_convergence.py first, or use --purpose spinup.")
        if start_year <= int(cutoff):
            raise SystemExit(
                f"orbit {start_year} is not past this run's assessed "
                f"equilibrium cutoff at {int(cutoff)}. Assess it again after "
                "the spin-up you are about to add, or use --purpose spinup.")

    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    landmap = surface_sra(config, 172).resolve()
    topomap = surface_sra(config, 129).resolve()
    # The same stack the prepare gives the threaded model. A continuation
    # launches the binary itself, so without this every segment after the first
    # would run on the 8 MB default -- which is exactly the shape of the defect
    # SHORTWAVE_GAS_KEYS' comment above describes for the namelist keys.
    thread_stack = prepare_thread_stack(config)
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
        inityear=start_year,
        workdir=str(run_dir),
        modelname=identifier,
        outputtype=model_cfg["output_type"],
    )
    # IMMEDIATELY AFTER THE CONSTRUCTION, because that is what overwrote it:
    # `exo.Model.__init__` copies the registry's executable into the working
    # directory every time one is made. For an arm run this puts the arm back;
    # for an ordinary run it is the prepare's binary_manifest.json check asked
    # again at the boundary where a rebuild lands. world-bdb5, world-u5pf.
    segment_exe, segment_exe_provenance = restore_run_executable(
        run_dir, model_cfg, manifest)
    model.configure(
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
        # DECLARED, not inherited. Earthlike.configure supplies vtype=4 and
        # modeltop=50.0 from its own signature, so NEQSIG and PTOP reached the
        # model from a library subclass default that no document here named.
        # Passed explicitly so the vertical grid and the model top are config,
        # and verified below like every other config-set key. world-a05.
        vtype=int(model_cfg["vertical_grid"]),
        modeltop=float(model_cfg["model_top_hpa"]),
        timestep=float(model_cfg["timestep_minutes"]),
        physicsfilter=model_cfg["physics_filter"],
        # RESTATED, and they have to be. `configure()` writes FILTERKAPPA and
        # NFILTEREXP unconditionally from its OWN defaults (8.0 and 8), so a
        # continuation that passes only `physicsfilter` rewrites the filter
        # exponent to 8 while the config declares 16. The filter is this model's
        # small-scale damping, so that is a change in the integrated dynamics
        # partway through a run, not a lost switch. world-8bs.
        filterkappa=float(model_cfg["filter_kappa"]),
        filterpower=int(model_cfg["filter_power"]),
        landmap=str(landmap),
        topomap=str(topomap),
        restartfile=str(restart),
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=(int(derived["snapshot_interval_steps"]) if args.seasonal_output else 0),
        highcadence=(
            {"toggle": 1, "start": 1,
             "end": int(derived["runsteps_per_orbit"]) * args.orbits,
             "interval": int(args.high_cadence_interval)}
            if args.high_cadence else
            {"toggle": 0, "start": 0, "end": 0, "interval": 4}),
        # IMPORTED, not restated. See configure_otherargs' docstring: this was
        # a second copy carrying only N_DAYS_PER_YEAR, and CLIM-17's TFREEZE
        # reverted to the compiled Earth value on every segment because of it.
        otherargs=configure_otherargs(derived),
    )
    # Constructing the model on an existing run directory re-copies the shipped
    # namelists over the configured ones, so the spectrum has to be staged again
    # here and not only at prepare time. Without it `NSTARFILE` stays 0 and
    # `solarini` falls back to a blackbody at `STARBBTEMP`, which is what every
    # continuation on this build did: 0.3844 of the flux below 0.75 um on the
    # first orbit and 0.4184 on every orbit after it.
    stage_stellar_spectrum(model, run_dir, stellar_spectrum_path(config))
    verify_stellar_spectrum(model, config)

    # Without this a continuation silently drops the diagnostics: nenergy is
    # namelist state that configure() rebuilds, and the codes are not in
    # REGULAR_CODES. The run would keep going and the terms would simply stop
    # appearing partway through, which is the failure mode that is hardest to
    # notice in a long spin-up.
    o3 = config["model"].get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        model._edit_namelist("radmod_namelist", "O3SCALE", f"{float(o3)}")
    # The shortwave gas band weights, IMPORTED rather than restated.
    # configure() rewrites the namelist on every continuation, so a weight not
    # reapplied here stops applying partway through a run -- and that is not
    # hypothetical: this was a second copy of the tuple, PHYS-9 added
    # h2o_sw_level to the other one, and H2OSWL silently reverted to the model
    # default on every segment of run_78c22fb1a1bd after the first.
    for key, name, default in SHORTWAVE_GAS_KEYS:
        w = config["model"].get(key)
        if w is not None and float(w) != default:
            model._edit_namelist("radmod_namelist", name, f"{float(w)}")

    # The longwave trace gases go the same way and for the same reason. Missing
    # this is exactly the failure the comment above records: CH4 and N2O default
    # to 0.0, meaning ABSENT, so a segment that does not reapply them integrates
    # a world with no methane in it while the manifest says otherwise. CLIM-42.
    for name, ppmv in trace_gas_ppmv(config).items():
        model._edit_namelist("radmod_namelist", name, f"{ppmv:.6g}")

    # Every continuation re-runs configure(), which rewrites the namelist, so
    # this has to be reapplied here and not only at prepare time. It is also
    # independent of the energy diagnostics, which it used to be nested inside.
    eco_stream = declare_ecological_stream(
        model, args.ecological_stream, args.eco_interval_steps)
    if eco_stream["enabled"]:
        print(f"  ecological stream ON for this segment: NECO = 1, NECOSTEP = "
              f"{eco_stream['necostep']} ({eco_stream['interval']}). Read "
              f"directly, not through pyburn.")
    set_low_io(model, args.low_io)
    if args.low_io:
        print("  NLOWIO = 1 for this segment: cheaper, and every orbit carries "
              "interval accumulations rather than instantaneous samples. "
              "build_climatology.py refuses these orbits without "
              "--allow-low-io, which is the guard rather than a nuisance.")
    else:
        print("  NLOWIO = 0 for this segment: instantaneous samples, for orbits "
              "something will read as data.")

    # CONS-9, and this is the path the defect was actually on: configure()
    # rewrites the namelists on every continuation, so a key that is not
    # reapplied here silently reverts partway through a run. Checked against the
    # config after re-staging and before the segment starts.
    # RE-APPLY BEFORE VERIFYING. `configure()` rewrites the namelists on every
    # continuation, so every declared switch has to be set again here -- and the
    # verification below is only meaningful once they have been. It ran first,
    # which made it fail on a key this script was about to write and, worse,
    # pass on the ones it never wrote at all: the energy fixer and the
    # dynamics-only switches were dropped on every continuation and nothing
    # said so. failure-modes class 22, which is the class this script's own
    # comment above cites.
    regular_codes = list(REGULAR_CODES)
    if energy_diagnostics_enabled(config):
        enable_energy_diagnostics(model, config)
        register_energy_diagnostic_codes()
        regular_codes = regular_codes + ENERGY_DIAGNOSTIC_CODES
        if config["model"].get("energy_diagnostics_3d", False):
            regular_codes = regular_codes + ENERGY_3D_CODES
    declare_dynamics_only(model, config)
    # THE THIRD INSTANCE of the defect this file already records twice, for the
    # stellar spectrum and for the shortwave gas weights: configure() re-copies
    # the shipped namelists over the configured ones on every continuation, so
    # anything not reapplied here stops applying partway through a run.
    #
    # Hyperdiffusion is the worst of the three to lose, because the fallback is
    # not "no damping" but the compiled branch at plasim.f90:1443 -- a DIFFERENT
    # OPERATOR. At T21 that is ndel 2 against the configured 4, grad^4 instead
    # of grad^8, with humidity damped 7.4x harder and vorticity 2.0x. The
    # baseline this project's climatology rests on integrated 84 of its 85
    # orbits that way. world-1nz.
    declare_hyperdiffusion(model, config)
    declare_energy_fixer(model, config)
    declare_robert_filter(model, config)
    declare_conversion_time_level(model, config)
    declare_dealias_conversion(model, config)
    # CLIM-54 and CLIM-55, on the same footing as every switch above: the
    # shipped hurricane_namelist that configure() re-copies happens to carry
    # both keys at 0, so this segment would integrate the same either way. It is
    # reapplied so that the value is this project's declaration rather than a
    # coincidence, and so verify_staged_namelists below has something to check.
    declare_storm_diagnostics(model)
    # AKAP, T0, TGR, DTROP, ALR and TFRC are on the same footing as the switches
    # above: configure() rewrites planet_namelist and plasim_namelist here, so a
    # continuation that did not re-apply them would revert to p_earth.f90's
    # constants partway through a run. TGR reaching the orographic surface
    # pressure reduction makes that a change in the mean state, not only in the
    # dry column.
    declare_dry_constants(model, config)
    # The hyperdiffusion is on exactly the same footing, and was the half this
    # script never wrote at all: `configure()` does not touch NDEL, NHDIFF or
    # TDISS*, so a continuation left them at `readnl`'s compiled T21/T42 branch
    # and integrated the derived damping on the prepare segment and ExoPlaSim's
    # on every segment after it. world-8bs.
    declare_hyperdiffusion(model, config)
    staged_namelists = verify_staged_namelists(run_dir, config)
    print(f"  namelists verified: {len(staged_namelists)} config-set keys "
          f"present with the declared values")
    refuse_eco_codes(regular_codes, "REGULAR_CODES")
    refuse_eco_codes(SNAPSHOT_CODES, "SNAPSHOT_CODES")
    model._add_postcodes("example.nl", regular_codes)
    model.cfgpostprocessor(
        ftype="regular",
        extension=model_cfg["output_type"],
        variables=[str(code) for code in regular_codes],
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

    if args.high_cadence:
        model._add_postcodes("highcadence.nl", HIGH_CADENCE_CODES)
        model.cfgpostprocessor(
            ftype="highcadence",
            extension=model_cfg["output_type"],
            variables=[str(code) for code in HIGH_CADENCE_CODES],
            mode="grid",
            times=None,
            timeaverage=False,
            interpolatetimes=False,
        )

    stage_surface_extras(run_dir, config)
    surface_field_report(run_dir, config)
    # Presence is not the property that failed. `surface_field_report` above
    # says the .sra files are in the run directory, which is the whole of what
    # a cold start needs; a resume takes every one of those fields out of the
    # restart instead, and nothing compared the two. CLIM-67 was a restart
    # branch that read the per-cell dwmax and then assigned a namelist scalar
    # over it, with every file, manifest and report still well-formed. This
    # compares CONTENT, on every code this project stages. CLIM-70.
    surface_restart = verify_restart_surface_fields(
        run_dir, restart, intended_surface_codes(config), manifest=manifest,
        allow_superseded=args.superseded_surface_ok)
    print(f"  restart surface verified against the "
          f"{surface_restart['reference']} reference: "
          f"{surface_restart['matched']} of {len(surface_restart['codes'])} "
          f"codes carry the field they were built from")
    started = datetime.now(timezone.utc).isoformat()
    try:
        # `clean` deletes the RAW outputs once pyburn has written the netCDF,
        # and for a high-cadence segment the raw file is the deliverable: the
        # netCDF pyburn writes is the ordinary twelve-bin average, and
        # `aeolian/scripts/extract_high_cadence_wind.py` reads the raw stream
        # because that is where the per-sample records are. So a high-cadence
        # segment keeps everything and drops the raws it does not need itself.
        model.run(years=args.orbits, crashifbroken=True,
                  clean=not args.high_cadence)
        if args.high_cadence:
            for year in range(start_year, start_year + args.orbits):
                for name in (f"MOST.{year:05d}", f"MOST_SNAP.{year:05d}"):
                    raw = run_dir / name
                    if raw.is_file():
                        raw.unlink()
            kept = sorted((run_dir / "highcadence").glob("MOST_HC.[0-9]*"))
            kept = [p for p in kept if p.suffix == ""]
            print(f"  high-cadence raw kept: {[p.name for p in kept]}")
            print("  extract it with aeolian/scripts/extract_high_cadence_wind.py")
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

    # A diagnostic segment says nothing about where the run is settling, so it
    # must not move the run's status either way: a run that was equilibrated
    # before a three-orbit I/O check is equilibrated after it.
    if args.purpose == "post_equilibrium_climatology":
        manifest["status"] = "climatology_complete"
    elif args.purpose == "spinup":
        manifest["status"] = "spinup_in_progress"
    manifest["completed_orbits"] = start_year + args.orbits

    # The binary that ran is the one sitting in the run directory: ExoPlaSim
    # copies it there and runs it in place, so this is what integrated the
    # orbits above rather than what some index says should have. Re-hashed here
    # rather than reusing what `restore_run_executable` settled, so the segment
    # records the file as it stands after the run and not as it was checked.
    segment_exe_sha = file_sha256(segment_exe) if segment_exe.is_file() else None
    prepared_exe_sha = (manifest.get("executable") or {}).get("sha256")
    if segment_exe_sha and prepared_exe_sha and segment_exe_sha != prepared_exe_sha:
        print(f"NOTE: this segment was integrated by {segment_exe_sha[:16]}, and the run "
              f"was prepared with {prepared_exe_sha[:16]}. The binary was rebuilt at some "
              "point in this run's life. That is recorded per segment and is legitimate; "
              "it is NOT legitimate across the low-I/O restart-layout change, which an "
              "older restart cannot survive.")

    manifest.setdefault("segments", []).append(
        {
            "start_year_index": start_year,
            "end_year_index": start_year + args.orbits - 1,
            "thread_stack": thread_stack,
            "seasonal_output": args.seasonal_output,
            "low_io": bool(args.low_io),
            "ecological_stream": eco_stream,
            "high_cadence": bool(args.high_cadence),
            "purpose": args.purpose,
            # What star these particular orbits were integrated against. The
            # top-level digest is what the run was PREPARED on and is what the
            # resume guard compares; this is per segment because the two came
            # apart once already, when every continuation reverted the namelist
            # to a blackbody and only the later segments ran on the file.
            "stellar_spectrum_digest": stellar_spectrum_digest(config),
            "started_utc": started,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "input_restart_sha256": file_sha256(restart),
            # WHICH BINARY integrated these particular orbits. The top-level
            # `executable` block is what the run was PREPARED with, and a
            # continuation is where the two come apart: CLAUDE.md rule 4 says
            # rebuild every binary after any change under vendor/exoplasim, so
            # a run that spans a rebuild has segments no single executable
            # produced. Recorded per segment for the same reason the stellar
            # spectrum digest above is, and it matters more here, because the
            # low-I/O change alters the restart layout and a resume across it
            # is not merely unattributed but wrong.
            "executable_sha256": segment_exe_sha,
            # AND WHOSE it was: the manifest entry it matched, or the arm it is.
            # A sha alone cannot say whether the binary was a registered one,
            # and a segment integrated by an unregistered binary is
            # unattributable in exactly the way world-qnue's stability probe
            # entries are.
            "executable_provenance": segment_exe_provenance,
            # WHICH SURFACE these particular orbits were integrated on, checked
            # inside the restart rather than taken from the run directory. A
            # run that adopted a donor's surface keeps integrating it while the
            # staged .sra beside it is refreshed on every resume, so this is
            # the only per-segment record of which of the two the model read.
            "surface_restart_check": surface_restart,
            "diagnostics": new_diagnostics,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(new_diagnostics, indent=2))


if __name__ == "__main__":
    main()
