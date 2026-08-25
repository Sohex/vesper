"""Check that an artifact from another component describes the world we are in.

Pointing one module at another's output is a decision, and it should look like
one. This is the check that makes it deliberate instead of assumed.

## Why

The components here hand each other files: the climate model's climatology
drives pedology, hydrography and the biosphere; hydrography's coupling matrix
drives the carve verdict; pedology's soil map drives LPJ-GUESS. Every one of
those pairings is only meaningful if both sides describe the SAME terrain under
the SAME climate, and nothing about a NetCDF file on disk makes that visible.

The failure is silent by construction. A superseded climatology has the same
grid, the same variables and the same units as a current one, so every consumer
reads it happily and produces a number that is simply about a different planet.
This project has hit that repeatedly -- a hardcoded climatology default that
outlived the terrain it was built on, a flat per-component directory holding
whichever build was active when it was last written, a coupling matrix paired
with another build's basin catalogue.

## The two mechanisms, and when to use which

**Namespacing** is the stronger one and should be preferred: write per-build,
under `<component>/data/<source_build>/`, and resolve with
`builds.component_data(..., strict=True)`. A mismatch is then impossible rather
than merely detectable, because the wrong file is not at the path at all.

**Stamping and checking** is for artifacts that cannot be namespaced -- usually
because they are large, shared, or named by something other than the build. Those
carry their identity as attributes or in a provenance sidecar, and the consumer
verifies it here. That is what this module is for.

Use `require_build` at the point of reading, not at the end. A check that runs
after the expensive part has already used the wrong input still wastes the run.

And do not put a silent default across a component boundary. A default that
resolves to another component's latest output is an assumption wearing the
costume of a convenience: it is right until the day it is not, and on that day
it produces a number rather than an error. Require the argument, or resolve it
per build and strictly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Stamped onto climatologies by exoplasim/scripts/build_climatology.py.
BUILD_ATTR = "vesper_source_build"
GEOGRAPHY_ATTR = "vesper_geography"


def active_build(config: dict | None = None) -> str:
    """The build every component is supposed to be working on."""
    if config is None:
        import yaml
        config = yaml.safe_load(
            (PROJECT_ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    build = config.get("source_build")
    if not build:
        raise SystemExit("config/planet.yaml names no source_build")
    return build


def artifact_build(path: Path) -> str | None:
    """Which build an artifact says it came from, or None if it does not say.

    Understands a NetCDF file with the stamped attributes, a JSON provenance
    document, and a data file sitting beside a `*_provenance.json`. Returns None
    rather than raising when an artifact carries no identity at all: that is a
    different problem from carrying the wrong one, and the caller decides how
    strict to be about it.
    """
    path = Path(path)
    if path.suffix == ".json":
        try:
            return json.loads(path.read_text(encoding="utf-8")).get("source_build")
        except (OSError, json.JSONDecodeError):
            return None

    sidecar = path.with_name(path.stem + "_provenance.json")
    if sidecar.is_file():
        try:
            return json.loads(sidecar.read_text(encoding="utf-8")).get("source_build")
        except (OSError, json.JSONDecodeError):
            return None

    if path.suffix == ".nc":
        try:
            import netCDF4 as nc
            with nc.Dataset(path) as data:
                if BUILD_ATTR in data.ncattrs():
                    return str(data.getncattr(BUILD_ATTR))
        except Exception:
            return None
    return None


def require_build(path: Path, what: str, config: dict | None = None,
                  allow_unstamped: bool = True) -> str | None:
    """Raise unless `path` describes the active build.

    `what` names the artifact in the error, because "does not match" is useless
    without saying which of several inputs is the wrong one.

    `allow_unstamped` governs artifacts produced before identity was stamped, or
    by a generator that does not stamp. The default permits them with a warning
    printed, because refusing outright would make every pre-existing artifact
    unusable; pass False where the pairing is expensive or changes the terrain,
    and an unidentifiable input should stop the run.
    """
    want = active_build(config)
    got = artifact_build(path)
    if got is None:
        if allow_unstamped:
            print(f"  warning: {what} ({Path(path).name}) carries no build "
                  f"identity, so it cannot be checked against {want}")
            return None
        raise SystemExit(
            f"{what} ({path}) carries no build identity and this consumer "
            f"requires one. Regenerate it, or pass the check explicitly.")
    if got != want:
        raise SystemExit(
            f"{what} ({path}) was built from {got!r}, but config/planet.yaml "
            f"names {want!r}. These describe different worlds; pairing them "
            f"would silently mix one terrain's rows with another's. Rebuild it "
            f"for {want!r}, or change source_build deliberately.")
    return got


def config_drift(recorded: dict, current: dict,
                 inert: frozenset[str] | set[str] = frozenset(),
                 path: str = "") -> list[str]:
    """Semantic differences between two parsed configurations, deepest first.

    Returns one `key: old -> new` line per differing leaf, and an empty list when
    the two configurations mean the same thing. Keys in `inert` are skipped; a
    block name skips the whole block, since the name is tested before recursing.

    Compare PARSED VALUES, never a hash of `config/planet.yaml`. A file hash
    cannot tell an edited comment from an edited parameter, so every guard built
    on one reports a stale artifact for a documentation change. That fired on the
    resume guard first and on `check_consistency.py`'s biosphere check second,
    which is why the function lives here instead of in either of them: one
    mechanism, so the answer to "has the config moved under this artifact" cannot
    differ between the two places that ask it.

    `inert` is per CONSUMER and is deliberately not shared, because reachability
    is a property of the consumer, not of the key. `baseline_climatology` cannot
    change a run in flight and is inert for a resume; it names the climatology
    `build_vesper_header.py` fits the solstice offset against, so it is not inert
    for the biosphere. A shared list would have to be the intersection, and the
    intersection is the one nobody checks.

    A key earns a place in an `inert` set only by being traced to nothing, and
    the trace belongs in a comment beside it. `unknown_inert_keys` is the check
    that the entry at least names a live key: `star.surface_uv` sat in the resume
    guard's list for as long as the guard existed, matching
    `star.surface_uv_relative_to_earth` never.
    """
    out = []
    for key in sorted(set(recorded) | set(current)):
        full = f"{path}{key}"
        if full in inert:
            continue
        a, b = recorded.get(key), current.get(key)
        if isinstance(a, dict) and isinstance(b, dict):
            out += config_drift(a, b, inert, f"{full}.")
        elif a != b:
            out.append(f"{full}: {a!r} -> {b!r}")
    return out


# The cartographic declaration: what the world's spin and its zero meridian ARE,
# recorded so the question has an answer, and read by nothing that makes an
# artifact. Traced the same way as everything else here, `grep -c` returning 0
# for all three names in `build_surface_albedo.py`, `build_surface_soil_water.py`,
# `build_boundary_conditions.py`, `build_surface_roughness.py`,
# `build_vesper_header.py`, `build_vesper_pfts.py`, `build_lpj_driver.py`,
# `run_exoplasim.py` and `continue_exoplasim.py`.
#
# Unioned into each consumer's set rather than shared as one, because
# reachability is a property of the consumer: if a generator ever learns to
# orient a field from these, it drops the group and states what it reads. The
# honesty guard in `check_consistency.py` re-runs the grep and fails the moment
# one of them does.
#
# `model` is the wrong home and so is `SURFACE_UNREAD_MODEL_KEYS`: these are
# `planet.` keys, and a set whose name promises one block should not carry
# another.
CARTOGRAPHIC_DECLARATION_KEYS = frozenset({
    "planet.rotation_direction",
    "planet.longitude_positive",
    "planet.prime_meridian",
})


# Moved here from `scripts/check_consistency.py` so that ONE registry answers
# "has the config moved under this artifact" for every consumer that asks.
# `scripts/pipeline.py` is the second asker and appeared after this set did.
BIOSPHERE_INERT_CONFIG_KEYS = {
    "schema_version",                     # bookkeeping; written into provenance
                                          # records, never read as an input
    "star.activity",                      # a design declaration; no script
                                          # reads it
    "star.surface_uv_relative_to_earth",  # a design declaration; ExoPlaSim
                                          # models no ultraviolet and no script
                                          # reads it
    # Prose carried inside the YAML rather than in a comment above it, and the
    # in-band twin of the comment edit this check used to fail on. Only
    # `period_earth_years` and `amplitude_flux_peak_to_peak` are read out of
    # `stellar_cycle.components`, by `run_stellar_cycle.py`.
    "stellar_cycle.components.medium.note",
    "stellar_cycle.components.long.note",
} | CARTOGRAPHIC_DECLARATION_KEYS


# Per-CONSUMER inert config keys for the STAGED SURFACE FIELDS, keyed by the
# `config/pipeline.yaml` step id that writes them. Read by both
# `scripts/check_consistency.py` and `scripts/pipeline.py`, so the answer to
# "has the config moved under this field" cannot differ between the two places
# that ask it -- the same requirement `config_drift`'s docstring states.
#
# Every entry below is TRACED, not assumed: the block or key does not appear in
# the builder at all, and the trace is the grep in the comment beside it. A key
# that is merely believed harmless does not belong here; over-flagging costs a
# generator that runs in seconds, and under-flagging costs a climate run.
#
# `model` is deliberately NOT blanket-inert for any of them. All three read
# `config["model"]` wholesale and then pull keys out of it, so marking the block
# inert would hide a real change; only the keys traced to the run scripts alone
# are listed.

SURFACE_UNREAD_MODEL_KEYS = {
    # DERIVED, not asserted: each key below is one whose NAME does not appear
    # anywhere in the generator's source, so the generator cannot be reading it.
    # `check_consistency.py` re-runs that grep and FAILS if any of these names
    # turns up in the file, which is what keeps the list from going quietly
    # wrong when a builder learns to read a new key.
    #
    # The limitation, stated because it is the way this can be wrong: the trace
    # is TEXTUAL and per file. A key reached indirectly through a helper in
    # `lib/` that takes the whole config would not appear here and would be
    # wrongly inert. Re-derive against the helper too if one starts doing that.
    #
    # `model.compile_flags` sits beside `model.precision_bytes` in all four sets and for
    # the same reason: both describe how the MODEL BINARY is compiled, and
    # these four builders make staged surface fields out of the export without
    # compiling anything. Traced the same way, `grep -c compile_flags`
    # returning 0 in each of the four generators. It is the whole block rather
    # than its leaves because config_drift tests a block name before recursing,
    # so naming the block skips everything under it.
    "surface_albedo": frozenset({
        # Read only by `build_surface_soil_water.py --lakes`, which writes the
        # bucket depth over the lake fraction. Traced: `grep -c lake_dwmax_m`
        # is 0 in this generator.
        "model.lake_dwmax_m",
        "model.co2_sw_weight", "model.energy_diagnostics",
        "model.energy_diagnostics_3d", "model.h2o_sw_level",
        "model.h2o_sw_level_bracket",
        "model.h2o_sw_weight", "model.layers", "model.ncpus",
        "model.output_type", "model.ozone_scale",
        "model.ozone_uv_weight", "model.ozone_visible_weight",
        "model.physics_filter", "model.cold_start_seed",
        "model.hyperdiffusion",
        "model.cold_start_profile",
        "model.semi_implicit_reference_temperature_k",
        "model.rayleigh_sponge_rotations",
        "model.resolution_timestep_minutes",
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
        "model.ozone_height_m", "model.ozone_spread_m",
        "model.energy_fixer",
        "model.conversion_time_level",
        "model.dealias_conversion",
        "model.robert_filter",
        "model.compile_flags",
        "model.precision_bytes",
        "model.regular_output_bins_per_orbit", "model.roughness_source",
        "model.seasonal_samples_per_orbit", "model.soil_water_source",
        "model.timestep_minutes", "model.uniform_land_surface",
        "model.vegetation_albedo_bracket",
        # world-9m5. tree_albedo and grass_albedo themselves are READ here and
        # are deliberately absent; only their brackets are unread.
        "model.tree_albedo_bracket", "model.grass_albedo_bracket",
        # world-36g. The band pair moved out of this generator when the split
        # became a per-material RATIO: the shapes come from
        # analysis/rock_albedo_bands.json and analysis/vegetation_albedo.json
        # and are applied to whatever level each material carries, so the
        # config pair is now read only by run_exoplasim, as ALBFOREST.
        "model.vegetation_albedo_bands",
    }),
    "surface_roughness": frozenset({
        # Read only by `build_surface_soil_water.py --lakes`, which writes the
        # bucket depth over the lake fraction. Traced: `grep -c lake_dwmax_m`
        # is 0 in this generator.
        "model.lake_dwmax_m",
        "model.co2_sw_weight", "model.energy_diagnostics",
        "model.energy_diagnostics_3d", "model.geography_land_threshold",
        "model.h2o_sw_level", "model.h2o_sw_weight",
        "model.h2o_sw_level_bracket",
        "model.latitudes", "model.layers",
        "model.lithology_albedo_overrides", "model.longitudes",
        "model.ncpus", "model.output_type", "model.ozone_scale",
        "model.ozone_uv_weight", "model.ozone_visible_weight",
        "model.physics_filter", "model.cold_start_seed",
        "model.hyperdiffusion",
        "model.cold_start_profile",
        "model.semi_implicit_reference_temperature_k",
        "model.rayleigh_sponge_rotations",
        "model.resolution_timestep_minutes",
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
        "model.ozone_height_m", "model.ozone_spread_m",
        "model.energy_fixer",
        "model.conversion_time_level",
        "model.dealias_conversion",
        "model.robert_filter",
        "model.compile_flags",
        "model.precision_bytes",
        "model.regular_output_bins_per_orbit",
        "model.seasonal_samples_per_orbit", "model.soil_water_source",
        "model.timestep_minutes", "model.uniform_land_surface",
        "model.vegetation_albedo", "model.vegetation_albedo_bands",
        "model.vegetation_albedo_bracket",
        # world-9m5. The two cover endmembers moved out of
        # build_surface_albedo.py's argparse defaults and into the config
        # beside model.vegetation_albedo, so they are inert wherever that is.
        "model.tree_albedo", "model.tree_albedo_bracket",
        "model.grass_albedo", "model.grass_albedo_bracket",
    }),
    "surface_soil_water": frozenset({
        # Not `model.` keys, and the first non-model entries here. Both
        # describe OCEAN HEAT TRANSPORT in the climate model (CLIM-16);
        # this builder makes a soil water field out of the export and
        # compiles nothing. Traced the same way, and the re-grep in
        # check_consistency.py now covers non-model prefixes so the
        # claim is checked rather than merely asserted.
        "ocean.horizontal_diffusion",
        "ocean.horizontal_diffusivity_m2_s",
        "model.barren_rock_classes", "model.co2_sw_weight",
        "model.energy_diagnostics", "model.energy_diagnostics_3d",
        "model.geography_land_threshold", "model.h2o_sw_level",
        "model.h2o_sw_level_bracket",
        "model.h2o_sw_weight", "model.land_albedo_source",
        "model.layers", "model.lithology_albedo_overrides",
        "model.ncpus", "model.output_type", "model.ozone_scale",
        "model.ozone_uv_weight", "model.ozone_visible_weight",
        "model.physics_filter", "model.cold_start_seed",
        "model.hyperdiffusion",
        "model.cold_start_profile",
        "model.semi_implicit_reference_temperature_k",
        "model.rayleigh_sponge_rotations",
        "model.resolution_timestep_minutes",
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
        "model.ozone_height_m", "model.ozone_spread_m",
        "model.energy_fixer",
        "model.conversion_time_level",
        "model.dealias_conversion",
        "model.robert_filter",
        "model.compile_flags",
        "model.precision_bytes",
        "model.regular_output_bins_per_orbit", "model.roughness_source",
        "model.seasonal_samples_per_orbit", "model.timestep_minutes",
        "model.vegetation_albedo", "model.vegetation_albedo_bands",
        "model.vegetation_albedo_bracket",
        # world-9m5. The two cover endmembers moved out of
        # build_surface_albedo.py's argparse defaults and into the config
        # beside model.vegetation_albedo, so they are inert wherever that is.
        "model.tree_albedo", "model.tree_albedo_bracket",
        "model.grass_albedo", "model.grass_albedo_bracket",
    }),
    "boundary_conditions": frozenset({
        # Read only by `build_surface_soil_water.py --lakes`, which writes the
        # bucket depth over the lake fraction. Traced: `grep -c lake_dwmax_m`
        # is 0 in this generator.
        "model.lake_dwmax_m",
        "model.barren_rock_classes", "model.co2_sw_weight",
        "model.energy_diagnostics", "model.energy_diagnostics_3d",
        "model.h2o_sw_level", "model.h2o_sw_weight",
        "model.h2o_sw_level_bracket",
        "model.land_albedo_source", "model.layers",
        "model.lithology_albedo_overrides", "model.ncpus",
        "model.output_type", "model.ozone_scale",
        "model.ozone_uv_weight", "model.ozone_visible_weight",
        "model.physics_filter", "model.cold_start_seed",
        "model.hyperdiffusion",
        "model.cold_start_profile",
        "model.semi_implicit_reference_temperature_k",
        "model.rayleigh_sponge_rotations",
        "model.resolution_timestep_minutes",
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
        "model.ozone_height_m", "model.ozone_spread_m",
        "model.energy_fixer",
        "model.conversion_time_level",
        "model.dealias_conversion",
        "model.robert_filter",
        "model.compile_flags",
        "model.precision_bytes",
        "model.regular_output_bins_per_orbit", "model.roughness_source",
        "model.seasonal_samples_per_orbit", "model.soil_water_source",
        "model.timestep_minutes", "model.uniform_land_surface",
        "model.vegetation_albedo", "model.vegetation_albedo_bands",
        "model.vegetation_albedo_bracket",
        # world-9m5. The two cover endmembers moved out of
        # build_surface_albedo.py's argparse defaults and into the config
        # beside model.vegetation_albedo, so they are inert wherever that is.
        "model.tree_albedo", "model.tree_albedo_bracket",
        "model.grass_albedo", "model.grass_albedo_bracket",
    }),
}

# Top-level config BLOCKS a staged-field builder never touches, traced by the
# same grep as the model keys above and listed by hand because a block name
# skips a whole subtree and deserves to be read rather than computed.
_SURFACE_BLOCKS = {
    # build_surface_albedo.py mentions `orbit` (baseline_flux_earth) and
    # `baseline_climatology`, so neither is inert for it.
    "surface_albedo": frozenset({
        "star", "atmosphere", "radiation", "surface", "ocean", "stellar_cycle",
    }),
    # build_surface_roughness.py mentions no block but `model`.
    "surface_roughness": frozenset({
        "planet", "star", "orbit", "atmosphere", "radiation", "surface",
        "ocean", "stellar_cycle", "baseline_climatology",
    }),
    # build_surface_soil_water.py mentions planet, surface, ocean and
    # baseline_climatology. `surface` and `ocean` are almost certainly substring
    # hits in prose rather than reads, and they are left OUT of the inert set
    # anyway: the textual trace over-detects a read, which costs a regeneration
    # that runs in seconds, and under-detecting one costs a climate run.
    "surface_soil_water": frozenset({
        "star", "orbit", "atmosphere", "radiation", "stellar_cycle",
    }),
    # build_boundary_conditions.py mentions `planet`, so it is not inert there.
    "boundary_conditions": frozenset({
        "star", "orbit", "atmosphere", "radiation", "surface", "ocean",
        "stellar_cycle", "baseline_climatology",
    }),
}


SURFACE_INERT_CONFIG_KEYS = {
    step: blocks | SURFACE_UNREAD_MODEL_KEYS[step] | CARTOGRAPHIC_DECLARATION_KEYS
    for step, blocks in _SURFACE_BLOCKS.items()
}


def config_stamp(config: dict, generator: str) -> dict:
    """The provenance every generated input carries, in one place.

    `CLAUDE.md` requires every generated product to record its provenance, and
    three staged surface fields did not, which is why config drift under them
    was unobservable rather than merely present. This is the shape the biosphere
    generators already write and the shape `config_drift` consumes.

    `source_config` is the PARSED configuration and is what a drift check reads.
    `config_sha256` is the file hash and is kept only so a record written before
    a generator learned to keep the parsed config can still be compared at all;
    it cannot separate an edited comment from an edited parameter, which is why
    nothing prefers it. `generator` is here so the remedy is read off the
    artifact rather than guessed -- one message naming one script for three
    artifacts was wrong for two of them.
    """
    import subprocess
    root = PROJECT_ROOT
    cfg_path = root / "config" / "planet.yaml"
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                                capture_output=True, text=True,
                                check=True).stdout.strip()
    except Exception:
        commit = None
    return {
        "generator": generator,
        "source_config": config,
        "config_sha256": hashlib.sha256(cfg_path.read_bytes()).hexdigest(),
        "git_commit": commit,
    }


def unknown_inert_keys(inert, config: dict) -> list[str]:
    """Entries of an `inert` set that name nothing in `config`.

    An allowlist entry that matches no key is not harmless: it reads as a
    decision that was made and it silently does nothing, so the key it was meant
    to excuse still blocks. Every entry must name a live key.
    """
    def present(dotted: str) -> bool:
        node = config
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return False
            node = node[part]
        return True
    return sorted(k for k in inert if not present(k))


# The one registry, keyed by the `config/pipeline.yaml` step id that writes the
# artifact. A step ABSENT from here is not checked for drift, and that is the
# safe default rather than an oversight: an inert set is a claim about what a
# consumer reads, `config_drift`'s docstring requires each entry to be traced,
# and an untraced empty set would flag every step on any parameter change and
# teach people to re-run generators to silence a message.
INERT_CONFIG_KEYS = dict(SURFACE_INERT_CONFIG_KEYS)
INERT_CONFIG_KEYS.update({
    "vesper_header": BIOSPHERE_INERT_CONFIG_KEYS,
    "vesper_pfts": BIOSPHERE_INERT_CONFIG_KEYS,
    "lpj_driver": BIOSPHERE_INERT_CONFIG_KEYS,
})


def artifact_drift(path, config: dict, inert) -> list[str] | None:
    """Config drift under one generated artifact, or None if it carries no stamp.

    Returns the `config_drift` lines, an empty list when nothing moved, and None
    when the artifact records no `source_config` at all -- which is a different
    answer and must not be collapsed into "fine". An artifact with no stamp is
    unobservable rather than current, and the caller decides whether that is a
    warning or a failure.
    """
    from pathlib import Path as _P
    path = _P(path)
    if path.suffix != ".json" or not path.is_file():
        return None
    try:
        rec = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(rec, dict):
        return None
    was = rec.get("source_config")
    if was is None:
        return None
    return config_drift(was, config, inert)
