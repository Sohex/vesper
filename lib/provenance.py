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
