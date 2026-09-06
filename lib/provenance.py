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

**A path keyed by something OTHER than the build** is the third case and the
staged surface fields are it. `exoplasim/inputs/<rung>/orogen_<RUNG>_surf_<code>.sra`
is keyed by the RUNG, and `surface_albedo` rewrites it per BUILD, so the
directory holds one build's field at a time and the filename cannot say which.
Namespacing is not available -- the model reads that path -- so
`staged_surface_field` is the one door: it checks the build at the read and
returns the record a consumer stamps into its own product. A DELIBERATE
cross-build read declares itself by NAMING the build it means, never by a flag
that turns the check off.

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

from paths import rel  # noqa: E402  same package directory

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


def build_stamp(build: str) -> dict:
    """What a product must record about the build it was made on.

    READING A REFUSED BUILD IS LEGAL AND HAS TO STAY LEGAL. A build is
    disposable until a climate run has consumed it, results already computed
    from a superseded terrain stay readable and datable, and
    `notes/audits/orogen-resolution.md` measured one on purpose. Registration
    and activation are separate for that reason, and this is not a refusal.

    What goes wrong is later, not now: a product made on a refused build lands
    in `<component>/data/<name>` beside a production one with nothing telling
    them apart, and the mistake is made months afterwards by whoever reads the
    artifact. So the ARTIFACT carries the verdict. `activation_refusal` is the
    registry's own words when the build may not be named as `source_build`,
    and None when it may; `registered` is False for a build the registry has
    never seen, which is a third state and not a refusal.
    """
    import orogen as _orogen
    entry = _orogen.registry_entry(name=build)
    return {
        "build": build,
        "terrain_hash": _orogen.terrain_hash_for_name(build),
        "registered": entry is not None,
        "activation_refusal": _orogen.activation_refusal(name=build),
    }


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
                 path: str = "", removed: dict | None = None) -> list[str]:
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
    is a property of the consumer, not of the key. Neither climatology key can
    change a run in flight, so both are inert for a resume; `bootstrap_climatology`
    names the climatology `build_vesper_header.py` fits the solstice offset
    against and `baseline_climatology` names the one `build_lpj_driver.py` reads,
    so neither is inert for the biosphere. A shared list would have to be the
    intersection, and the intersection is the one nobody checks.

    A key earns a place in an `inert` set only by being traced to nothing, and
    the trace belongs in a comment beside it. `unknown_inert_keys` is the check
    that the entry at least names a live key: `star.surface_uv` sat in the resume
    guard's list for as long as the guard existed, matching
    `star.surface_uv_relative_to_earth` never.

    `removed` is a different question and is answered separately. A key present
    in `recorded` and absent from `current` reads here as `{...} -> None`, which
    is what an edited parameter reads as too, so this cannot tell a REMOVAL from
    a value change and must refuse both. `REMOVED_CONFIG_KEYS` is where a
    removal is declared in advance with the value the key held, and both halves
    of that declaration are checked before it excuses anything -- see
    `_removal_verdict`. It is not a relaxation of the comparison and there is no
    flag: an undeclared removal is still drift, and so is a declared one under a
    value the declaration does not record.
    """
    if removed is None:
        removed = REMOVED_CONFIG_KEYS
    out = []
    for key in sorted(set(recorded) | set(current)):
        full = f"{path}{key}"
        if full in inert:
            continue
        a, b = recorded.get(key), current.get(key)
        if key in recorded and key not in current:
            out += _removal_verdict(full, a, removed)
            continue
        if isinstance(a, dict) and isinstance(b, dict):
            out += config_drift(a, b, inert, f"{full}.", removed)
        elif a != b:
            out.append(f"{full}: {a!r} -> {b!r}")
    return out


def _removal_verdict(full: str, was, removed: dict) -> list[str]:
    """What a key present in `recorded` and absent from `current` is worth.

    Both halves of the declaration are checked, and a failure of either is drift
    rather than a warning, because the direction that must stay loud is the one
    that refuses.
    """
    spec = removed.get(full)
    if spec is None:
        return [f"{full}: {was!r} -> removed, and no removal is declared for it. "
                f"A key that vanished from the configuration is indistinguishable "
                f"here from one whose value changed; declare the removal in "
                f"lib/provenance.py:REMOVED_CONFIG_KEYS, with the value it held, "
                f"or restore the key."]
    if not any(was == v for v in spec["values"]):
        return [f"{full}: {was!r} -> removed by {spec['removed_by']}, whose "
                f"declaration records the values this key held as "
                f"{spec['values']!r}. This artifact holds none of them, so the "
                f"key's VALUE moved under it before the key itself went. That is "
                f"drift, and the removal does not excuse it."]
    return []


def applied_removals(recorded: dict, current: dict, removed: dict | None = None,
                     path: str = "") -> list[str]:
    """Which declared removals a comparison of these two configurations rested on.

    `config_drift` returns only problems, so a resume that crossed a declared
    removal would otherwise say nothing at all -- and a declaration nobody sees
    at the moment it is used is the shape of a silence, whatever is written in
    the source. The caller prints these.
    """
    if removed is None:
        removed = REMOVED_CONFIG_KEYS
    out = []
    for key in sorted(set(recorded) | set(current)):
        full = f"{path}{key}"
        a, b = recorded.get(key), current.get(key)
        if key in recorded and key not in current:
            spec = removed.get(full)
            if spec is not None and any(a == v for v in spec["values"]):
                out.append(f"{full}, removed by {spec['removed_by']}; "
                           f"{spec['moved_to_note']}")
        elif isinstance(a, dict) and isinstance(b, dict):
            out += applied_removals(a, b, removed, f"{full}.")
    return out


# A KEY THAT NO LONGER EXISTS ANYWHERE, declared, so the guard can tell it from
# a key whose value changed.
#
# WHAT THIS IS FOR. `config_drift` compares an artifact's recorded configuration
# against the current one, and a key that was REMOVED reads there as
# `{...} -> None`, which is exactly what an edited parameter reads as. The guard
# cannot tell them apart and must therefore refuse both, which is the right way
# round: a false resume integrates orbits under a configuration nobody declared
# and is silent, a false refusal costs a restart and is loud. That refusal
# killed a commissioning in flight when world-f997 moved the ladder's timestep
# table into `lib/rungs.py`, and at the top of the escalation route the same
# event costs most of a day. world-51wj.
#
# WHAT IS NOT WANTED IS A LOOSER GUARD. Nothing below relaxes the comparison.
# A removal is stated in advance, in the source, with the value the key held
# copied verbatim, and the guard then checks BOTH halves -- the same shape the
# fork's `mainline_divergences` files use, which is the one that has held here:
#
#   the key must be GONE from the current configuration. A declaration for a key
#   that is still there is stale, and `removal_problems` fails on it.
#
#   the artifact's recorded value must be one of the values declared below,
#   verbatim. A run whose manifest holds anything else had the key's VALUE move
#   under it before the key went, which is real drift and is still refused.
#
# So a removal cannot be declared without reading what the configuration held,
# and it cannot be used to wave through a value that changed. There is no flag.
#
# THE RESIDUAL RISK, stated because it is the way this can be wrong. The
# declaration asserts that nothing on the run path read the key, so that its
# removal cannot change what the model integrates. `removal_problems` checks
# that assertion the only mechanical way available -- the key's name appears in
# no Python source in the tree outside the paths `moved_to` names -- and that is
# a TEXTUAL trace, the same limitation `SURFACE_UNREAD_MODEL_KEYS` records. A
# key reached through a helper that takes the whole configuration and looks the
# name up by construction would not be found. `settles` is where the human half
# of the argument goes, and it is not optional.
#
# Each entry carries:
#   removed_by      the issue that removed the key
#   owner           the issue that declared the removal
#   values          every value `config/planet.yaml` is known to have held for
#                   this key, verbatim. A run manifest matching none of them is
#                   still refused
#   moved_to        repo-relative paths where the meaning now lives, and the
#                   only paths exempt from the name trace. Empty means the
#                   quantity is gone rather than moved
#   moved_to_note   one line, printed at the resume that rests on this entry
#   settles         why the removal cannot change what a run integrated
REMOVED_CONFIG_KEYS = {
    "model.resolution_timestep_minutes": {
        "removed_by": "world-f997",
        "owner": "world-51wj",
        "values": ({"T21": 45.0, "T42": 45.0, "T85": 45.0,
                    "T127": 30.0, "T170": 22.5},),
        "moved_to": ("lib/rungs.py",),
        "moved_to_note": "the per-rung ceiling now lives in lib/rungs.py",
        "settles": (
            "A TABLE OF CEILINGS that nothing read at run time. It recorded the "
            "coarsest step each rung was measured to carry, one row per rung, "
            "while the step a run is actually configured at is "
            "`model.timestep_minutes` -- a different key, which world-f997 did "
            "not touch. `lib/rungs.py` now holds the ceiling, and its docstring "
            "is where the three quantities the phrase 'the timestep at rung X' "
            "has named in this tree are kept apart. Nothing on the run or the "
            "resume path ever asked the configuration for this key: its name "
            "appears in no Python source in the tree."),
    },
    "surface.land_water_column.layer_capacity_fraction": {
        "removed_by": "world-6nla",
        "owner": "world-6nla",
        "values": ([0.013333, 0.32, 0.666667], [0.333333, 0.666667], [1.0]),
        "moved_to": ("exoplasim/scripts/run_exoplasim.py",),
        "moved_to_note": "the fallback capacity shape is now derived from "
                         "`layer_thickness_m` in run_exoplasim.py",
        "settles": (
            "ONE LINE OF ARITHMETIC ON THE NEXT LINE. The key was the "
            "thickness-proportional share a uniform profile gives, each entry "
            "of `layer_thickness_m` over the column's total, and `landini` "
            "renormalises whatever it is given, so only the SHAPE reached the "
            "model. `run_exoplasim.py` now computes that shape from the "
            "thicknesses and writes the same `DSOILWF` at the same six "
            "significant digits, which reproduces every value recorded above "
            "exactly at the thicknesses that produced them. What the removal "
            "ends is the case that already happened: the surface cut moved the "
            "thicknesses and left the shape behind."),
    },
    "model.cold_start_profile.lapse_rate_k_per_m": {
        "removed_by": "world-6nla",
        "owner": "world-6nla",
        "values": (0.008489,),
        "moved_to": ("lib/lapse.py", "exoplasim/scripts/run_exoplasim.py"),
        "moved_to_note": "the cold start's lapse rate is now derived by "
                         "lib/lapse.py:cold_start_lapse_rate_k_per_m",
        "settles": (
            "A FRACTION OF THIS ATMOSPHERE'S OWN DRY ADIABAT, and the fraction "
            "is what transfers rather than the rate: Earth's standard "
            "6.5 K/km is a fixed share of Earth's g/cp, and the same share of "
            "this planet's g/cp is what the key held. Both sides of that are "
            "already declared here -- the composition and the gravity -- so "
            "the key was a second statement of a derivation the configuration "
            "already determines. `run_exoplasim.py` writes `ALR` from it and "
            "records what it wrote on the run manifest, so a reader still "
            "finds the number a run integrated; what is gone is the copy that "
            "could not learn the composition had moved."),
    },
    "model.cold_start_profile.tropopause_height_m": {
        "removed_by": "world-6nla",
        "owner": "world-6nla",
        "values": (9220.0,),
        "moved_to": ("lib/lapse.py", "exoplasim/scripts/run_exoplasim.py"),
        "moved_to_note": "the cold start's tropopause height is now derived by "
                         "lib/lapse.py:cold_start_tropopause_height_m",
        "settles": (
            "A FIXED NUMBER OF PRESSURE SCALE HEIGHTS, which is the invariant "
            "the place a profile turns isothermal has to hold across a change "
            "of atmosphere. The key was upstream's 12 km scaled by this "
            "atmosphere's R T / g over Earth's, so it is a function of the "
            "cold start's SURFACE TEMPERATURE and of the composition and "
            "gravity -- all three declared here, and the surface temperature "
            "one line above where the key used to sit. It went stale exactly "
            "that way: the value recorded above was computed at a surface "
            "temperature the configuration no longer carries."),
    },
    "model.ozone_height_m": {
        "removed_by": "world-6nla",
        "owner": "world-6nla",
        "values": (15311.8,),
        "moved_to": ("exoplasim/scripts/run_exoplasim.py",),
        "moved_to_note": "BO3 is now derived by "
                         "run_exoplasim.py:ozone_profile_lengths_m",
        "settles": (
            "UPSTREAM'S OWN COMPILED LENGTH, SCALED. `radmod.f90` declares "
            "`bo3` and `mko3` places the ozone profile against it as a "
            "geometric height fitted to Earth; what has to be held across a "
            "change of atmosphere is the PRESSURE that height corresponds to, "
            "so the key was `bo3` times this atmosphere's (R/g) over Earth's. "
            "Both halves are readable -- the length out of the model source, "
            "the ratio out of `lib/lapse.py` on the declared composition and "
            "gravity -- so the key was a transcription of a product of two "
            "things the tree already holds, and it could not learn that either "
            "had moved."),
    },
    "model.ozone_spread_m": {
        "removed_by": "world-6nla",
        "owner": "world-6nla",
        "values": (3828.0,),
        "moved_to": ("exoplasim/scripts/run_exoplasim.py",),
        "moved_to_note": "CO3 is now derived by "
                         "run_exoplasim.py:ozone_profile_lengths_m",
        "settles": (
            "The same statement as `model.ozone_height_m` about `radmod.f90`'s "
            "`co3`: the second of the two lengths `mko3` places the ozone "
            "profile against, scaled by the same ratio for the same reason, "
            "and derived beside it in one function so the two cannot part."),
    },
}


def removal_problems(config: dict, root=None, removed: dict | None = None) -> list[str]:
    """Every declared removal that the tree no longer bears out.

    A declaration excuses a refusal, so a stale one is worse than none: it is a
    decision that was made once and keeps applying to a tree that has moved. The
    same argument `unknown_inert_keys` makes about an allowlist entry naming
    nothing, in the opposite direction -- here the failure is an entry naming
    something that is BACK.

    Three ways an entry can be wrong, and each fails by name:

      the key is in `config/planet.yaml` again, so nothing was removed
      `moved_to` names a path that does not exist
      the key's name appears in a Python source outside `moved_to`, so something
      may read it and the removal's whole claim is contradicted
    """
    if removed is None:
        removed = REMOVED_CONFIG_KEYS
    root = PROJECT_ROOT if root is None else Path(root)
    dirs = [root / "scripts", root / "lib", root / "maps"] + [
        root / c / "scripts" for c in
        ("exoplasim", "hydrography", "pedology", "biosphere", "minerals", "aeolian")]

    def present(dotted: str) -> bool:
        node = config
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return False
            node = node[part]
        return True

    out = []
    for key, spec in sorted(removed.items()):
        if present(key):
            out.append(f"{key} is declared removed by {spec['removed_by']} and "
                       f"config/planet.yaml still has it")
        # This module holds the declaration, so it names every removed key by
        # construction. Exempting it here rather than by a path in `moved_to`
        # keeps `moved_to` meaning "where the quantity went".
        exempt = {Path(__file__).resolve()}
        for rel_path in spec["moved_to"]:
            path = root / rel_path
            if not path.exists():
                out.append(f"{key} declares it moved to {rel_path}, which is "
                           f"not there")
            exempt.add(path.resolve())
        leaf = key.rsplit(".", 1)[-1]
        for d in dirs:
            if not d.is_dir():
                continue
            for f in sorted(d.glob("*.py")):
                if f.resolve() in exempt:
                    continue
                if leaf in f.read_text(encoding="utf-8", errors="ignore"):
                    out.append(f"{key} is declared removed and {rel(f)} names "
                               f"it, so the claim that nothing reads it is "
                               f"contradicted by the source")
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
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        # world-et25. The three constants that select between a form the
        # model derives and a declared literal: rainmod_nl's GAMMA and
        # RCRITWIDTH and fluxmod_nl's VDIFF_LAMM. All three are run-script
        # keys for the same reason clwref above is, and they arrived in the
        # config as declarations rather than as new numbers -- the value the
        # model integrates has not moved. Traced the same way,
        # `grep -c precip_reevaporation_gamma`,
        # `grep -c cloud_fraction_subgrid_width` and
        # `grep -c asymptotic_mixing_length_m` each returning 0 in all four.
        "model.precip_reevaporation_gamma",
        "model.cloud_fraction_subgrid_width",
        "model.asymptotic_mixing_length_m",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
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
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        # world-et25. The three constants that select between a form the
        # model derives and a declared literal: rainmod_nl's GAMMA and
        # RCRITWIDTH and fluxmod_nl's VDIFF_LAMM. All three are run-script
        # keys for the same reason clwref above is, and they arrived in the
        # config as declarations rather than as new numbers -- the value the
        # model integrates has not moved. Traced the same way,
        # `grep -c precip_reevaporation_gamma`,
        # `grep -c cloud_fraction_subgrid_width` and
        # `grep -c asymptotic_mixing_length_m` each returning 0 in all four.
        "model.precip_reevaporation_gamma",
        "model.cloud_fraction_subgrid_width",
        "model.asymptotic_mixing_length_m",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
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
        # `model.geography_land_threshold` is NOT listed here, and it is the
        # one entry these sets lost to a change in what the generator reads.
        # WORLD-QGB6 moved this builder off a climatology and onto surface code
        # 0172, and the mask is built from the threshold, so the value is
        # load-bearing on what is staged even though the key is never pulled out
        # of `config["model"]` here: with `--lakes` the blend is written only
        # where the mask owns land, so a threshold change moves code 0229's
        # bytes; without it the threshold still decides the unmatched-land
        # refusal and the `land_cells` and `land_mean_capacity_m` the record
        # carries. It does NOT move code 2290, whose split is a property of the
        # soil column and takes `land` only into its report.
        # THE INPUT-HASH CHAIN DOES NOT COVER IT, which is why the key comes out
        # of the set rather than resting on `inputs`: the mask's hash moves only
        # once `build_boundary_conditions.py` has re-run, so between the config
        # edit and that rerun the 0229 record would report itself current on a
        # mask cut at the old threshold.
        "model.h2o_sw_level",
        "model.h2o_sw_level_bracket",
        "model.h2o_sw_weight", "model.land_albedo_source",
        # `model.layers` is NOT listed here, and it is the one entry in these
        # four sets whose textual trace cannot be made to hold. The trace is a
        # substring match on the leaf name, and this generator's whole subject
        # is the land water column's LAYERS: the word is ordinary English in
        # every second line of it, so a claim that the name is absent would be
        # true only until the next sentence. It is dropped rather than defended,
        # which costs a regeneration of a generator that runs in seconds
        # whenever the atmosphere's vertical layer count moves. That is the
        # cheap direction this list's own header names.
        "model.lithology_albedo_overrides",
        "model.ncpus", "model.output_type", "model.ozone_scale",
        "model.ozone_uv_weight", "model.ozone_visible_weight",
        "model.physics_filter", "model.cold_start_seed",
        "model.hyperdiffusion",
        "model.cold_start_profile",
        "model.semi_implicit_reference_temperature_k",
        "model.rayleigh_sponge_rotations",
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        # world-et25. The three constants that select between a form the
        # model derives and a declared literal: rainmod_nl's GAMMA and
        # RCRITWIDTH and fluxmod_nl's VDIFF_LAMM. All three are run-script
        # keys for the same reason clwref above is, and they arrived in the
        # config as declarations rather than as new numbers -- the value the
        # model integrates has not moved. Traced the same way,
        # `grep -c precip_reevaporation_gamma`,
        # `grep -c cloud_fraction_subgrid_width` and
        # `grep -c asymptotic_mixing_length_m` each returning 0 in all four.
        "model.precip_reevaporation_gamma",
        "model.cloud_fraction_subgrid_width",
        "model.asymptotic_mixing_length_m",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
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
        "model.filter_kappa",
        "model.filter_power", "model.cloud_absorption_scale",
        # world-n1nu. The CCM3 cloud water reference, a rainmod_nl key
        # the run scripts write; these builders make a staged surface
        # field out of the export and no cloud scheme runs in them.
        # Traced the same way, `grep -c cloud_water_reference_kg_m3`
        # returning 0 in each of the four.
        "model.cloud_water_reference_kg_m3",
        # world-et25. The three constants that select between a form the
        # model derives and a declared literal: rainmod_nl's GAMMA and
        # RCRITWIDTH and fluxmod_nl's VDIFF_LAMM. All three are run-script
        # keys for the same reason clwref above is, and they arrived in the
        # config as declarations rather than as new numbers -- the value the
        # model integrates has not moved. Traced the same way,
        # `grep -c precip_reevaporation_gamma`,
        # `grep -c cloud_fraction_subgrid_width` and
        # `grep -c asymptotic_mixing_length_m` each returning 0 in all four.
        "model.precip_reevaporation_gamma",
        "model.cloud_fraction_subgrid_width",
        "model.asymptotic_mixing_length_m",
        "model.land_longwave_emissivity", "model.sea_longwave_emissivity",
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
    # `baseline_climatology`, so neither is inert for it. `bootstrap_climatology`
    # does not appear in it at all: the script calls no resolver, takes the file
    # from --climatology so it can run before any climatology exists, and the
    # one it must be given is the baseline the LPJ-GUESS driver was built from.
    "surface_albedo": frozenset({
        "star", "atmosphere", "radiation", "surface", "ocean", "stellar_cycle",
        "bootstrap_climatology",
    }),
    # build_surface_roughness.py mentions no block but `model`.
    "surface_roughness": frozenset({
        "planet", "star", "orbit", "atmosphere", "radiation", "surface",
        "ocean", "stellar_cycle", "baseline_climatology",
        "bootstrap_climatology",
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
        "stellar_cycle", "baseline_climatology", "bootstrap_climatology",
    }),
}


SURFACE_INERT_CONFIG_KEYS = {
    step: blocks | SURFACE_UNREAD_MODEL_KEYS[step] | CARTOGRAPHIC_DECLARATION_KEYS
    for step, blocks in _SURFACE_BLOCKS.items()
}


def _staged_reports(rung_dir: Path):
    """Every provenance record beside a staged surface field, code by code.

    Yields `(code, report_path, record)`. The two naming conventions in the
    directory -- `*report*.json` from the albedo, roughness and boundary
    builders and `*_provenance.json` from the dust and soil-water builders --
    are BOTH scanned, because which one a generator chose is an accident of
    when it was written and a reader that knew only one would report an
    unstamped field for half the directory.
    """
    for path in sorted(list(rung_dir.glob("*report*.json"))
                       + list(rung_dir.glob("*_provenance.json"))):
        try:
            rec = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(rec, dict):
            continue
        codes = rec.get("codes")
        if codes is None and rec.get("code") is not None:
            codes = [rec["code"]]
        for code in codes or ():
            try:
                yield int(code), path, rec
            except (TypeError, ValueError):
                continue


def staged_surface_records(codes, config: dict | None = None, *,
                           root: Path | None = None) -> list[Path]:
    """The provenance records beside the staged fields for these codes.

    Deduplicated and sorted: one record can cover several codes, as
    `albedo_report.json` covers seven. Returns paths and reads nothing else, so
    a caller that wants to COPY the records rather than interpret them -- which
    is what `run_exoplasim.py` does, because the staged record is overwritten
    at the next staging while the run's copy of the field is not -- does not
    have to know which of the directory's two naming conventions a generator
    chose.
    """
    root = PROJECT_ROOT if root is None else Path(root)
    if config is None:
        import yaml
        config = yaml.safe_load(
            (root / "config" / "planet.yaml").read_text(encoding="utf-8"))
    rung_dir = (root / "exoplasim" / "inputs"
                / str(config["model"]["resolution"]).lower())
    wanted = {int(c) for c in codes}
    return sorted({path for code, path, _ in _staged_reports(rung_dir)
                   if code in wanted})


def _record_build_name(rec: dict) -> str | None:
    """The build a staged-field record NAMES, wherever it put the name.

    `build_surface_dust.py` writes `source_build` at the top level and the
    albedo, roughness and soil-water builders write it inside the
    `config_stamp` block as `source_config.source_build`. Both are the same
    assertion and `artifact_build` already treats either as the artifact's
    identity, so a resolver that recognised only one would disagree with the
    function beside it about what stamps an artifact.
    """
    for value in (rec.get("source_build"), rec.get("build"),
                  (rec.get("source_config") or {}).get("source_build")):
        if isinstance(value, str) and value:
            return value
    return None


def _record_terrain_hash(rec: dict, path: Path, code: int, root: Path):
    """The terrain hash a staged-field record asserts, by hash or by name.

    A RECORD THAT NAMES ITS BUILD IS STAMPED. The registry in `lib/orogen.py`
    is what turns a name into a terrain hash, and it is the same registry
    `build_stamp` writes the explicit hash from, so the two routes cannot
    disagree about a registered build. Accepting only the explicit hash made
    this door refuse `orogen_T21_surf_0229_provenance.json`, which carries the
    name and not the hash -- a correct-looking refusal on a field that was
    stamped, which pushes the caller into hand-rolling the check.

    THE CHECK THAT FIRES WHEN THEY DISAGREE is here rather than in the writer,
    because a record carrying both is the ordinary case and a record whose two
    halves name different builds is unusable in either direction.
    """
    import orogen as _orogen
    explicit = rec.get("terrain_hash")
    name = _record_build_name(rec)
    by_name = _orogen.terrain_hash_for_name(name) if name else None
    if explicit and by_name and explicit != by_name:
        raise SystemExit(
            f"{rel(path, root)} stamps code {code} with terrain hash "
            f"{str(explicit)[:16]} and names build {name!r}, which the registry "
            f"in lib/orogen.py resolves to {by_name[:16]}. A record cannot name "
            f"two builds; re-run the generator that wrote it.")
    return explicit or by_name


def run_surface_field(code: int, run: str, *, root: Path | None = None):
    """The staged `.sra` ONE RUN consumed, for one surface code.

    THE FAILURE THIS CLOSES. `exoplasim/inputs/<rung>/` holds the field the NEXT
    run will read, and a run's climatology is in equilibrium with the field
    THAT run read. The two are the same file only until the next restage, and a
    restage is an ordinary event: `build_surface_soil_water.py` rewrites code
    229 whenever pedology re-derives the column. An offline calculation that
    pairs a run's climatology with whatever sits in `exoplasim/inputs/` is then
    computing on two worlds, and nothing about either file says so. Measured on
    `run_67323a923013`: the staged code 229 and the field that run consumed
    differ on 1598 of 2048 cells, and the climatology's soil water exceeds the
    STAGED bucket on 55 land cells while fitting the run's own bucket exactly.

    The run's own `run_manifest.json` records `surface_field_sha256` per code
    and the run directory keeps its copy of every `.sra` it staged, so the
    answer is already written down and needs only a door. This is that door.
    The returned record names both hashes and `matches_staged`, so a product
    that stamps it says whether the tree has moved on since the run -- which is
    a fact to record rather than a refusal, because the run's field is the
    right answer either way.

    Raises `SystemExit` when the run, its manifest, its record of the code or
    the file itself is absent, and when the bytes on disk do not match what the
    manifest recorded.
    """
    root = PROJECT_ROOT if root is None else Path(root)
    run_dir = root / "exoplasim" / "runs" / str(run)
    if not run_dir.is_dir():
        raise SystemExit(
            f"no run directory at {rel(run_dir, root)}; "
            f"exoplasim/runs/INDEX.json is the record of what exists")
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(
            f"{rel(run_dir, root)} has no run_manifest.json, so which surface "
            f"fields it consumed is unrecorded and cannot be recovered from "
            f"the directory alone")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded = (manifest.get("surface_field_sha256") or {}).get(str(int(code)))
    if not recorded:
        from_file = (manifest.get("surface_fields") or {}).get("from_file") or []
        if int(code) not in {int(c) for c in from_file}:
            raise SystemExit(
                f"{rel(run_dir, root)} staged no surface code {code}: the run "
                f"read the namelist default for it, so there is no field to "
                f"pair with its climatology. run_manifest.json's "
                f"`surface_fields.uniform_values` says what the default was.")
        raise SystemExit(
            f"{rel(manifest_path, root)} lists code {code} among the fields "
            f"staged but records no sha256 for it, so the bytes in the run "
            f"directory cannot be checked against what the run read")

    nlat = int(manifest.get("source_config", {}).get("model", {})
               .get("latitudes") or 0)
    candidates = sorted(run_dir.glob(f"N*_surf_{int(code):04d}.sra"))
    if nlat:
        preferred = run_dir / f"N{nlat:03d}_surf_{int(code):04d}.sra"
        if preferred.is_file():
            candidates = [preferred]
    if not candidates:
        raise SystemExit(
            f"{rel(run_dir, root)} records a sha256 for surface code {code} but "
            f"holds no N*_surf_{int(code):04d}.sra; the run directory is the "
            f"only copy of what the run read and it is gone")
    path = candidates[0]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != recorded:
        raise SystemExit(
            f"{rel(path, root)} hashes to {digest[:16]} and "
            f"{rel(manifest_path, root)} records {recorded[:16]} for code "
            f"{code}. The run directory's copy has been rewritten since the run, "
            f"so it is not evidence of what the run read.")

    staged_path = (root / "exoplasim" / "inputs"
                   / str(manifest.get("physical", {}).get("resolution", "")).lower()
                   / f"orogen_{manifest.get('physical', {}).get('resolution', '')}"
                     f"_surf_{int(code):04d}.sra")
    staged_sha = (hashlib.sha256(staged_path.read_bytes()).hexdigest()
                  if staged_path.is_file() else None)
    # The generator's own record, if the run kept it. `run_exoplasim.py` copies
    # it in beside the field; a run taken before that did not, and the staged
    # original has been overwritten since, so None here means the derivation is
    # gone rather than that it can be looked up. world-hl06.
    kept = sorted({p for c, p, _ in _staged_reports(run_dir) if c == int(code)})
    return {
        "code": int(code),
        "run": str(run),
        "path": rel(path, root),
        "sha256": digest,
        "source_build": manifest.get("source_build"),
        "provenance_record": rel(kept[0], root) if kept else None,
        "staged_path": rel(staged_path, root) if staged_sha else None,
        "staged_sha256": staged_sha,
        "matches_staged": None if staged_sha is None else staged_sha == digest,
    }


def staged_surface_field(code: int, config: dict | None = None, *,
                         for_build: str | None = None, root: Path | None = None,
                         paired_with_run: str | None = None):
    """The staged `.sra` for one surface code, WITH the build it was staged from.

    `exoplasim/inputs/<rung>/orogen_<RUNG>_surf_<code>.sra` is keyed by the RUNG
    ALONE, and `surface_albedo` is a loop A step that rewrites it per build, so
    the directory holds exactly one build's field at a time and the path cannot
    say which. Five consumers built that path by hand from
    `model.resolution`. Rule 5 says pointing one component at another's output
    is deliberate and never defaulted; this is the door that makes it so.

    Returns the record a consumer must put in its own product:
    `path`, `sha256`, `terrain_hash`, `build`, `declared`. Stamping it is the
    other half -- the check makes a wrong read loud NOW, and the record makes
    the pairing auditable AFTERWARDS, which is what a verdict already written
    needs.

    **The deliberate cross-build read is declared by NAMING THE BUILD, not by
    turning the check off.** `for_build` is a build name from the registry in
    `lib/orogen.py`; the read then succeeds only against THAT build and fails
    against any other, including the active one. The carve overshoot
    measurement is a cross-build read by construction -- the staged albedo is
    the right one while the carved build is staged and the wrong one for
    anything re-run on the pre-carve build afterwards -- so it needs a way to
    say which build it means, and a boolean `--allow-any` would have let the
    same silence back in under a flag.

    **A record NAMES its build or carries its terrain hash, and either is a
    stamp.** `_record_terrain_hash` resolves a name through the registry in
    `lib/orogen.py` and refuses when a record carries both and they disagree.

    **`paired_with_run` is for a caller reading a RUN's climatology.** The
    staged field is what the next run will read; a run's climatology is in
    equilibrium with what that run read, and the two part company at every
    restage. Naming the run makes this door refuse the pairing rather than
    return a field from another iteration; `run_surface_field` is the door onto
    the field the run actually consumed.

    Raises `SystemExit` when the file is absent, when no provenance record
    beside it names the code, or when the build it was staged from is not the
    one asked for. An unstamped field is UNOBSERVABLE and is refused rather
    than trusted, because the whole failure this closes is silent.
    """
    import builds as _builds
    import orogen as _orogen
    from orogen import _KNOWN_TERRAIN_HASHES
    if config is None:
        import yaml
        config = yaml.safe_load(
            (PROJECT_ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    root = PROJECT_ROOT if root is None else Path(root)
    rung = str(config["model"]["resolution"]).upper()
    rung_dir = root / "exoplasim" / "inputs" / rung.lower()
    path = rung_dir / f"orogen_{rung}_surf_{int(code):04d}.sra"
    if not path.is_file():
        raise SystemExit(
            f"staged surface code {code} is not at {rel(path, root)}; run the "
            f"config/pipeline.yaml step that writes it")

    stamps = [(p, rec) for c, p, rec in _staged_reports(rung_dir) if c == int(code)]
    hashes = {_record_terrain_hash(rec, p, int(code), root) for p, rec in stamps}
    hashes.discard(None)
    if not hashes:
        raise SystemExit(
            f"{rel(path, root)} has no provenance record beside it naming code "
            f"{code}, so the build it was staged from is unknown. An unstamped "
            f"staged field is not the same as a current one; re-run its "
            f"generator so it writes its report.")
    if len(hashes) > 1:
        raise SystemExit(
            f"the records beside {rel(path, root)} disagree about which build "
            f"code {code} was staged from: {sorted(h[:16] for h in hashes)}")
    staged = hashes.pop()

    wanted = _builds.terrain_hash(config)
    declared = for_build is not None
    if declared:
        # orogen.py owns the registry and both its indexes. Rebuilding the
        # inverse here took the LAST of a duplicate pair, where the index in
        # orogen.py raises on one -- and a duplicate name is exactly the case
        # that cannot be resolved downstream, because a name addresses both a
        # source/ directory and a per-build data directory.
        wanted = _orogen.terrain_hash_for_name(for_build)
        if wanted is None:
            raise SystemExit(
                f"--for-build {for_build!r} is not a build in lib/orogen.py's "
                f"registry; a cross-build read is declared by NAMING the build "
                f"it means, so an unregistered name cannot declare anything")
    if staged != wanted:
        who = (_KNOWN_TERRAIN_HASHES.get(staged) or {}).get("name", "an "
                                                            "unregistered build")
        raise SystemExit(
            f"{rel(path, root)} was staged from {who} ({staged[:16]}), and this "
            f"read is against {for_build or active_build(config)} "
            f"({wanted[:16]}). The path is keyed by the rung alone, so the "
            f"directory holds one build's field at a time. Re-run the "
            f"surface_albedo step for the build you mean, or, if the "
            f"cross-build read is deliberate, declare it by naming the build.")

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if paired_with_run is not None:
        consumed = run_surface_field(int(code), paired_with_run, root=root)
        if consumed["sha256"] != digest:
            raise SystemExit(
                f"{rel(path, root)} is not the code {code} field "
                f"{paired_with_run} ran on: the run consumed "
                f"{consumed['sha256'][:16]} and the staged tree now holds "
                f"{digest[:16]}. A run's climatology is in equilibrium with the "
                f"field THAT run read, so pairing it with the staged one is a "
                f"read across two worlds. Call "
                f"lib/provenance.py:run_surface_field for the field the run "
                f"consumed, or drop `paired_with_run` if you mean the field the "
                f"NEXT run will read.")

    return {
        "code": int(code),
        "path": rel(path, root),
        "sha256": digest,
        "terrain_hash": staged,
        "build": (_KNOWN_TERRAIN_HASHES.get(staged) or {}).get("name"),
        "declared_cross_build": for_build,
        "paired_with_run": paired_with_run,
    }


INPUT_STAMP_KEY = "source_inputs"


def input_stamp(paths) -> dict:
    """`{INPUT_STAMP_KEY: {repo-relative path: sha256}}` for files a generator READ.

    THE CONFIG STAMP CANNOT SEE THESE AND THAT IS THE WHOLE REASON THIS EXISTS.
    `config_drift` compares parsed configuration values, which is the right
    check for a configuration key and is blind to a DERIVED FILE the generator
    also reads. A derivation re-run -- a new spectrum, a corrected proxy, a rock
    class added -- moves the file while every configuration key stays put, and
    the artifact built from the old file goes on describing the old shapes with
    nothing saying so. That is world-nvs2, seen first on the two band-shape
    files under `analysis/` that set staged surface codes 175 and 176.

    A file is hashed and not merely named, because the whole failure is that the
    NAME is stable across the change: `analysis/rock_albedo_bands.json` is the
    same path before and after its generator re-runs, exactly as
    `build_stellar_spectrum.py` rewrites `<name>.dat` in place.

    A path that does not exist is recorded as `None` rather than dropped. An
    absent input is a fact about the artifact, and dropping it would make a
    generator that later STOPS reading a file indistinguishable from one whose
    input has gone.

    Pass only files whose CONTENT the generator consumed. A mesh or grid export
    is covered by `terrain_hash` and does not belong here; a configuration file
    is covered by `source_config` and does not either.
    """
    out = {}
    for path in paths:
        path = Path(path)
        out[rel(path, PROJECT_ROOT)] = (
            hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file()
            else None)
    return {INPUT_STAMP_KEY: out}


def input_drift(recorded: dict, root: Path | None = None) -> list[str]:
    """Recorded inputs that have moved since the artifact was written.

    Takes the `{path: sha256}` mapping, not the whole record, and returns one
    line per input that has changed, appeared or gone -- empty when every
    recorded input is byte-identical to the file at its path. The caller decides
    what an empty record means: nothing to compare is not the same answer as
    nothing moved, and `artifact_input_drift` keeps the two apart by returning
    `None` for an artifact that carries no stamp at all.
    """
    root = PROJECT_ROOT if root is None else Path(root)
    lines = []
    for name, was in sorted((recorded or {}).items()):
        path = root / name
        now = (hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file()
               else None)
        if now == was:
            continue
        if was is None:
            lines.append(f"{name} did not exist when this was built and does now")
        elif now is None:
            lines.append(f"{name} was read at {was[:12]} and is now missing")
        else:
            lines.append(f"{name} has been rewritten since this was built "
                         f"({was[:12]} -> {now[:12]})")
    return lines


def artifact_input_drift(path, root: Path | None = None) -> list[str] | None:
    """Input drift under one generated artifact, or None if it carries no stamp.

    The same three-valued shape as `artifact_drift`, and for the same reason:
    an artifact that records no inputs is UNOBSERVABLE rather than current.
    """
    path = Path(path)
    if path.suffix != ".json" or not path.is_file():
        return None
    try:
        rec = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if not isinstance(rec, dict) or INPUT_STAMP_KEY not in rec:
        return None
    return input_drift(rec[INPUT_STAMP_KEY], root)


def config_stamp(config: dict, generator: str, inputs=None) -> dict:
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

    `inputs` are DERIVED FILES the generator read that no configuration key
    names. They are hashed by `input_stamp`, which carries the argument for why
    the config stamp cannot stand in for them.
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
    stamp = {
        "generator": generator,
        "source_config": config,
        "config_sha256": hashlib.sha256(cfg_path.read_bytes()).hexdigest(),
        "git_commit": commit,
    }
    # `inputs` is the half `source_config` cannot cover: derived files the
    # generator read, hashed so a re-derivation under an unchanged name is
    # visible. Absent means the generator reads no such file; an EMPTY sequence
    # means it does and none were present, and the two are written differently
    # so the check can tell them apart. `input_stamp` says why.
    if inputs is not None:
        stamp.update(input_stamp(inputs))
    return stamp


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


def _selftest() -> int:
    """The staged-field door's refusals, each with a right answer.

    Synthetic: a four-value `.sra` and hand-built records in a temporary
    directory, so nothing here depends on which build is staged or on a run
    existing. What it checks is the DISPOSITION of a record, which is the part
    that decided wrongly: the door read only an explicit `terrain_hash` and
    refused a record that named its build, which is the stamp four of the five
    records beside a staged field carry and the fifth does not.
    """
    import shutil
    import tempfile
    import orogen as _orogen

    problems: list[str] = []
    checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal checks
        checks += 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}"
              + (f"  ({detail})" if detail else ""))
        if not ok:
            problems.append(name)

    import yaml
    live = yaml.safe_load(
        (PROJECT_ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    build = active_build(live)
    known = _orogen.terrain_hash_for_name(build)
    if known is None:
        print(f"the configured source_build {build!r} is not in the registry, "
              f"so there is nothing to resolve a name against")
        return 1
    # The rung comes from the config rather than being written down here, for
    # the same reason no artifact path carries one: the fixture is about the
    # door's disposition and holds at whatever rung the tree is configured at.
    rung = str(live["model"]["resolution"]).upper()

    header = ("     229       0 20260811       0       4       1"
              "       0       0\n")
    root = Path(tempfile.mkdtemp())
    try:
        rung_dir = root / "exoplasim" / "inputs" / rung.lower()
        rung_dir.mkdir(parents=True)
        sra = rung_dir / f"orogen_{rung}_surf_0229.sra"
        sra.write_text(header + "  0.1  0.2  0.3  0.4\n")
        record = rung_dir / f"orogen_{rung}_surf_0229_provenance.json"
        config = {"model": {"resolution": rung}, "source_build": build}

        def door(rec: dict):
            record.write_text(json.dumps(rec))
            try:
                return staged_surface_field(229, config, root=root), None
            except SystemExit as refusal:
                return None, str(refusal)

        got, why = door({"code": 229, "source_config": {"source_build": build}})
        check("a record that NAMES its build is a stamp",
              got is not None and got["terrain_hash"] == known, why or "")

        got, why = door({"code": 229, "source_build": build})
        check("the name is read wherever the writer put it",
              got is not None and got["terrain_hash"] == known, why or "")

        got, why = door({"code": 229, "terrain_hash": known,
                         "source_config": {"source_build": build}})
        check("a record carrying both, agreeing, opens the door",
              got is not None and got["terrain_hash"] == known, why or "")

        got, why = door({"code": 229, "terrain_hash": "0" * 64,
                         "source_config": {"source_build": build}})
        check("a record whose hash and name are different builds is refused",
              got is None and "cannot name two builds" in (why or ""),
              (why or "")[:60])

        got, why = door({"code": 229, "field": "no stamp of any kind"})
        check("a record with no stamp at all is still refused",
              got is None and "no provenance record" in (why or ""),
              (why or "")[:60])

        runs = root / "exoplasim" / "runs" / "run_synthetic"
        runs.mkdir(parents=True)
        (runs / "run_manifest.json").write_text(json.dumps({
            "physical": {"resolution": rung},
            "source_config": {"model": {"latitudes": 4}},
            "source_build": build,
            "surface_fields": {"from_file": [229]},
            "surface_field_sha256": {
                "229": hashlib.sha256(sra.read_bytes()).hexdigest()},
        }))
        shutil.copy(sra, runs / "N004_surf_0229.sra")

        consumed = run_surface_field(229, "run_synthetic", root=root)
        check("a run's field resolves to the run's own copy, and says whether "
              "the staged tree has moved on",
              consumed["path"].endswith("N004_surf_0229.sra")
              and consumed["matches_staged"] is True, consumed["path"])
        check("a run that kept no provenance record says so rather than "
              "reaching for the staged one",
              consumed["provenance_record"] is None)

        shutil.copy(record, runs / record.name)
        consumed = run_surface_field(229, "run_synthetic", root=root)
        check("a run that kept the record names its own copy",
              consumed["provenance_record"] is not None
              and consumed["provenance_record"].endswith(record.name),
              str(consumed["provenance_record"]))

        try:
            run_surface_field(1720, "run_synthetic", root=root)
            check("a code the run never staged is refused as absent", False)
        except SystemExit as refusal:
            check("a code the run never staged is refused as absent",
                  "staged no surface code" in str(refusal))

        (runs / "N004_surf_0229.sra").write_text(header + "  0.9  0.9  0.9  0.9\n")
        try:
            run_surface_field(229, "run_synthetic", root=root)
            check("a rewritten copy in the run directory is refused", False)
        except SystemExit as refusal:
            check("a rewritten copy in the run directory is refused",
                  "has been rewritten since the run" in str(refusal))

        shutil.copy(sra, runs / "N004_surf_0229.sra")
        sra.write_text(header + "  0.5  0.5  0.5  0.5\n")
        record.write_text(json.dumps(
            {"code": 229, "source_config": {"source_build": build}}))
        try:
            staged_surface_field(229, config, root=root,
                                 paired_with_run="run_synthetic")
            check("a staged field that is not the one a named run read is refused",
                  False)
        except SystemExit as refusal:
            check("a staged field that is not the one a named run read is refused",
                  "two worlds" in str(refusal))
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print(f"\n{checks} checks, {len(problems)} failed")
    return 1 if problems else 0


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true",
                    help="run the staged-field door's refusals against their "
                         "right answers on synthetic input; needs no staged "
                         "field and no run on disk")
    args = ap.parse_args()
    if not args.selftest:
        ap.error("this module is a library; --selftest is the only thing to run")
    return _selftest()


if __name__ == "__main__":
    raise SystemExit(main())
