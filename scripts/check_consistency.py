#!/usr/bin/env python3
"""Preflight: does everything in tree describe the same world?

Run this before an expensive run and after changing `source_build`. It exists
because the same class of error kept arriving by different routes, and every
instance was silent: array shapes matched, fields looked plausible, and the
answer was wrong.

    python scripts/check_consistency.py

Exit status is 0 when everything agrees, 1 when something disagrees. Warnings
report things that are stale but harmless to the active configuration.

## What it checks, and which bug each one is for

**Terrain hash across artifacts.** Every product carrying a `terrain_hash` must
match the configured build. The carve verdict once ran with one terrain's basins
and another's coupling matrix and raised an IndexError purely because the counts
differed; had they matched it would have paired row *i* of one build with row *i*
of another and produced a number.

**The registry does not refuse the configured build.** `lib/orogen.py` registers
a build so that results computed from it stay readable, which is a different
question from whether the tree may point at it; entries that answer the second
one "no" carry a `refusal`. Config named one of those for as long as nobody
read the comment above its entry, and batch P0C measured the basin catalogue on
two builds because it could not tell which was active.

**Stale flat data directories.** Hydrography is per-build now. A leftover
`<component>/data/*.nc` from an earlier build is not an error, but it is what
`world_state.py` was reading when it reported 2,107 basins for a build with
2,540.

**Coupling convention.** The share of coupling catchment area landing on
model-ocean cells must stay under 10%, index for index, and the matrix numbers its columns
-180..180 while an ExoPlaSim climatology numbers its own 0..360. Every basin read
its antipode for an entire iteration.

**Surface inputs present and current.** Every code the config says it supplies
must exist, and must not predate the build it claims to describe.

**Carve list completeness.** A list sent to Orogen must cover the whole
catalogue with unique ids, and its carried-forward entries must all be retain 0,
since the loop is only monotone if an already-carved basin stays carved.

**Config internal consistency.** Gravity against mass, and the declared orbit
against the flux it is derived from.

**The configured soil saturation endpoints against the contract that derives
them.** `pedology/config/land_column_properties.yaml` derives the degree of
saturation at an empty and a full store as medians over the build's own land
cells, so they move with every soil rebuild. `run_exoplasim.py` already refuses
when `landmod.f90`'s compiled defaults drift from that contract; nothing refused
when `config/planet.yaml` did, and that is the copy the armed moisture-dependent
soil albedo reads. A key the config leaves unset is written at the compiled
default and is passed, because unset is one copy of the number instead of two.

**The configured flux against the DERIVED one.** `config/pipeline.yaml` makes
`baseline_run` need `design_flux`, but that edge only says the artifact exists:
adopting its winner into `config/planet.yaml` is a separate human act, because
the flux fixes the semi-major axis and that orbit is compiled into LPJ-GUESS. So
the edge alone leaves the failure available -- derive a flux, do not adopt it,
buy hours of baseline run at the old number. No `design_flux.json` at all is the
legitimate provisional state and warns; a winner that differs from the
configured value refuses.

**A design flux that was CHOSEN rather than derived.** The artifact says which
of the two it is and this refuses one that does not, because a flux the
derivation returned and a flux a person picked because the derivation refused
are different kinds of number. A chosen one must carry the ground the derivation
refused on, the finding it rests on and what would reopen it, and the ground has
to be one `derive_design_flux.py` itself declares choosable -- that list is read
from the generator, never restated here. The two laundering paths are checked
against a stubbed derivation rather than argued: a derivation that ANSWERS
cannot be overridden, and a refusal outside the declared set cannot be converted
into a choice.

**What the runs on disk actually integrated.** Three cases over run directories,
in `audit_runs`. The staged namelists against the config and the manifest each
run carries, the staged namelists against the hyperdiffusion block the manifest
declares, and the per-code `surface_field_sha256` against the `.sra` files that
are there now. Until these existed, `run_exoplasim.py`'s `expected_namelist_keys`
and `verify_staged_namelists` were the only guard on any of it and both fire only
at prepare and continue time, so a run already on disk was invisible: the T21 run
the project's numbers rest on records a hyperdiffusion block its namelist does not
carry, and a person reading run directories is what found it.

**The flux-to-kelvin slope against the runs it was measured on, and against the
build those runs are ON.** `lib/sensitivity.py` declares one sensitivity for the
whole project. Its predecessor was called FALLBACK_SLOPE, nothing ever fell back
to it, and it survived two terrain changes and a resolution change while three
incompatible values accumulated around it. This recomputes the declared slope
from the runs it names every time -- and, because recomputing agreed with itself
for a whole build after both named runs were deleted, it also asks whether those
runs are in the LIVE index and on the active `source_build`. The self-test beside
it drives the same checks with the superseded declaration and requires them to
fire, so the currency half is exercised rather than assumed.

**Generated biosphere inputs against the config values they were derived from.**
Parsed values, not a hash of `config/planet.yaml`, because a file hash reports an
edited comment as a stale artifact and the remedy it names -- re-run the
generators -- both works and teaches the wrong thing.

## Checking the checker

    python scripts/check_consistency.py --self-test

The audit above compares artifacts, so nothing in it says whether the comparison
itself can still tell a real change from a cosmetic one. `--self-test` asks that
directly: it edits the live `config/planet.yaml` IN MEMORY and asserts what the
comparison must say about each edit, and it builds a run directory that agrees
with itself and asserts what `audit_runs` must say about each way of breaking
it. A comment must not be drift, an allowlisted
key must not be drift, a parameter must be, a deleted key must be, and every
allowlist entry must name a key that exists. The parameter case is the one that
matters: a comparison that fired at nothing would pass every negative case here
and no positive one, which is why each negative case is paired with the positive
that proves the edit landed.

It also drives the activation guard with a build the registry is known to refuse.
That guard's live input is the configured build, so a case reading the config
would report whatever the tree happens to be in the middle of rather than
whether the guard works.

**The stellar band split.** `lib/stellar.py` must still reproduce
`radmod.f90:solarini`, the shipped spectrum file must still represent the
BT-Settl blend it was resampled from, and every consumer that stores the band-1
share must agree. The first is the identity at `radmod.f90:207`: a 5772 K
spectrum through this code gives 0.517. The second is the one with an outside
answer, and it is why the third is not enough on its own -- every stored copy
agreed with every other for as long as the resampler was point-sampling the
source, because they were all copies of the same wrong integral.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from netCDF4 import Dataset          # noqa: E402
import numpy as np                   # noqa: E402
import yaml                          # noqa: E402

import builds                        # noqa: E402
import rungs                         # noqa: E402
from paths import rel                # noqa: E402
from gridding import coupling_ocean_fraction   # noqa: E402
from provenance import applied_removals, artifact_drift, artifact_input_drift, BIOSPHERE_INERT_CONFIG_KEYS, INERT_CONFIG_KEYS, config_drift, REMOVED_CONFIG_KEYS, removal_problems, unknown_inert_keys   # noqa: E402


def land_sea_mask():
    """The climatology's land mask, for the coupling-alignment invariant."""
    import climatology as clim
    from paths import climatology_path
    with Dataset(climatology_path()) as ds:
        return clim.annual_mean_of(ds, "lsm")


FAIL, WARN, OK = "FAIL", "warn", "ok"


class Report:
    def __init__(self):
        self.rows: list[tuple[str, str, str]] = []

    def add(self, status, check, detail=""):
        self.rows.append((status, check, detail))

    @property
    def failed(self) -> bool:
        return any(s == FAIL for s, _, _ in self.rows)

    def show(self):
        width = max(len(c) for _, c, _ in self.rows) + 2
        for status, check, detail in self.rows:
            mark = {OK: "  ok  ", WARN: " warn ", FAIL: " FAIL "}[status]
            print(f"[{mark}] {check:<{width}} {detail}")
        n_fail = sum(1 for s, _, _ in self.rows if s == FAIL)
        n_warn = sum(1 for s, _, _ in self.rows if s == WARN)
        print(f"\n{len(self.rows)} checks, {n_fail} failed, {n_warn} warnings")


def nc_terrain(path: Path) -> str | None:
    try:
        with Dataset(path) as ds:
            return getattr(ds, "terrain_hash", None)
    except OSError:
        return None


def json_terrain(path: Path) -> str | None:
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return d.get("terrain_hash") or (d.get("source") or {}).get("terrain_hash")


def sha256_of(path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_runs(rep: "Report", runs: Path, runner) -> None:
    """What the runs already on disk actually integrated, read back off them.

    Three cases, and they exist because nothing looked. The guard that a run
    got the damping, the filter and the radiation its config asked for is
    `run_exoplasim.py`'s `expected_namelist_keys` plus
    `verify_staged_namelists`, and both fire only while a run is being prepared
    or resumed. A run already on disk -- prepared before the guard, seeded from
    a template, or edited -- was invisible to every check in the tree, and the
    runs here turned out not to agree with their own manifests.

    SEVERITY IS NOT A NEW JUDGEMENT. Each case applies the SAME predicate as
    the prepare-time guard it mirrors, at the same severity.
    `verify_staged_namelists` raises on a staged key whose value is not what
    the config asks for, so a wrong value here fails. The restart guard in
    `run_exoplasim.py` refuses to seed from a run whose recorded
    `surface_field_sha256` no longer matches the `.sra` on disk -- "content
    changed for code(s) ..." -- so a changed surface here fails.

    The one thing a retrospective check meets that a prepare-time one cannot
    is a run that PREDATES a key. `expected_namelist_keys` supplies the
    model's own default for a config that never declared the key, so an
    absent namelist entry on such a run means the model integrated its
    compiled value and NOTHING RECORDED WHICH. That is unobservable rather
    than wrong, and it warns -- the same answer this file already gives a
    climatology carrying no build identity. What separates the two is the
    run's own `namelist_checks` block: an entry there is the run asserting it
    staged that key and read it back, so the same key missing from the file
    is the manifest contradicted, and that fails.

    Takes the runs directory and the runner module rather than reading ROOT,
    so `--self-test` can drive it over a fixture whose right answer is known
    in advance instead of only over whatever happens to be on disk.
    """
    manifests = sorted(runs.glob("run_*/run_manifest.json"))

    # -- a run against the namelists it staged -------------------------------
    try:
        wrong, unrecorded, unstaged, seen = [], {}, [], 0
        for manifest in manifests:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            cfg = data.get("source_config")
            if not cfg:
                continue
            run = manifest.parent
            if not (run / "plasim_namelist").is_file():
                # Nothing was ever staged; the run died before `configure()`
                # wrote its namelists. There is no claim to contradict.
                unstaged.append(run.name)
                continue
            seen += 1
            claimed = set(data.get("namelist_checks") or {})
            bad = []
            for fname, keys in runner.expected_namelist_keys(cfg).items():
                path = run / fname
                for key, want in keys.items():
                    try:
                        raw = runner.namelist_value(path, key)
                    except (KeyError, FileNotFoundError):
                        if f"{key}@{fname}" in claimed or key in claimed:
                            bad.append(f"{key} claimed staged, absent from "
                                       f"{fname}")
                        else:
                            unrecorded.setdefault(run.name, set()).add(key)
                        continue
                    raw = raw.rstrip(",").strip()
                    # ONE comparison path, `run_exoplasim.namelist_elements`,
                    # which is what verify_staged_namelists uses at prepare
                    # time. It expands Fortran's `n*value` replication, so a
                    # replicated array and the same array written out element
                    # by element compare equal -- they are the same value to
                    # the model. LENGTH IS COMPARED FIRST, and that is what
                    # catches world-720: a scalar where NLEV elements are
                    # wanted sets element 1 and leaves levels 2..NLEV on the
                    # compiled values, and it fails here on its length rather
                    # than on its spelling.
                    wanted = ([float(v) for v in want]
                              if isinstance(want, (list, tuple))
                              else [float(want)])
                    try:
                        got = runner.namelist_elements(raw)
                    except ValueError:
                        bad.append(f"{key}={raw!r} is not a number")
                        continue
                    if len(got) != len(wanted):
                        bad.append(f"{key} has {len(got)} element(s), "
                                   f"config wants {len(wanted)}")
                        continue
                    for i, (g, w) in enumerate(zip(got, wanted)):
                        if abs(g - w) > 1e-9:
                            bad.append(f"{key}[{i + 1}]={g:g} not {w:g}")
                            break
            if bad:
                wrong.append(f"{run.name}: " + ", ".join(bad))
        if wrong:
            rep.add(FAIL, "runs vs the namelists they staged",
                    "; ".join(wrong[:4])
                    + (f" (+{len(wrong) - 4} more runs)" if len(wrong) > 4 else "")
                    + " -- the namelist in the run directory is what the "
                      "model integrated, so these runs did not integrate "
                      "what their own configuration asks for")
        else:
            rep.add(OK, "runs vs the namelists they staged",
                    f"{seen} runs carry namelists and every key their config "
                    f"declares is in them with the declared value")
        if unrecorded:
            lines = [f"{name}: {', '.join(sorted(keys))}"
                     for name, keys in sorted(unrecorded.items())]
            rep.add(WARN, "namelist keys the runs never recorded",
                    "; ".join(lines[:3])
                    + (f" (+{len(lines) - 3} more runs)" if len(lines) > 3 else "")
                    + " -- absent from the namelist and never claimed on the "
                      "manifest, so the model used its compiled value and no "
                      "artifact says which")
        if unstaged:
            rep.add(WARN, "runs that staged no namelists",
                    f"{', '.join(unstaged)} have a manifest and no "
                    f"plasim_namelist, so there is nothing to read back")
    except Exception as exc:
        rep.add(WARN, "runs vs the namelists they staged", f"not checked: {exc}")

    # -- a sentinel-selected constant against the number the model derived --
    #
    # world-et25. Three keys select a FORM rather than carrying a value:
    # `vdiff_lamm`, `gamma` and `rcritwidth` each derive when the namelist gives
    # a negative sentinel. So for those three the namelist records only that the
    # model chose, and a run whose record stops there cannot say what it
    # integrated -- two runs at different rotation rates or different rungs
    # carry the same `-1` and integrated different constants, and an arm on such
    # a constant would have to be a code fork.
    #
    # `run_exoplasim.py` reads each one back out of the model's own
    # initialisation print and stamps it under `derived_model_constants`. This
    # is the read-back: a run that staged one of the three and recorded no value
    # for it FAILS, and where the run staged a DECLARED literal the recorded
    # value must be that literal, which is a right answer the model already
    # knows. A run that never started has no print to read and is not counted.
    try:
        unrecorded, disagreed, seen = [], [], 0
        for manifest in manifests:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            run = manifest.parent
            if not any(run.glob("MOST_DIAG.*")):
                continue
            recorded = data.get("derived_model_constants") or {}
            for fname, key, _marker, name, _units in runner.DERIVED_CONSTANT_ECHOES:
                try:
                    staged = float(runner.namelist_value(run / fname, key)
                                   .rstrip(",").strip())
                except (KeyError, FileNotFoundError, ValueError, OSError):
                    # The run predates the route for this key entirely. That is
                    # the "namelist keys the runs never recorded" warning above
                    # and a different defect from this one.
                    continue
                seen += 1
                entry = recorded.get(name)
                if entry is None:
                    unrecorded.append(f"{run.name}: {key} = {staged:g}")
                    continue
                value = float(entry["value"])
                if staged < 0.0:
                    if value <= 0.0:
                        disagreed.append(
                            f"{run.name}: {key} staged at its sentinel and the "
                            f"recorded {name} is {value:g}, which is not a "
                            f"constant the model can have derived")
                elif abs(value - staged) > 1e-6 * max(1.0, abs(staged)):
                    disagreed.append(
                        f"{run.name}: {key} = {staged:g} was declared and the "
                        f"model echoed {value:g}. The model did not read the "
                        f"key the run staged")
        if unrecorded or disagreed:
            rep.add(FAIL, "runs vs the constants the model derived",
                    "; ".join((unrecorded + disagreed)[:4])
                    + (f" (+{len(unrecorded) + len(disagreed) - 4} more)"
                       if len(unrecorded) + len(disagreed) > 4 else "")
                    + " -- a sentinel says the model chose the number, so a run "
                      "that cannot state which number it chose is one nothing "
                      "can reconstruct or build an arm against")
        else:
            rep.add(OK, "runs vs the constants the model derived",
                    f"{seen} sentinel-selected keys across the runs that have "
                    f"started, each with the value the model printed for it")
    except Exception as exc:
        rep.add(WARN, "runs vs the constants the model derived",
                f"not checked: {exc}")

    # -- a run against the hyperdiffusion its manifest declares --------------
    #
    # world-e6m, and the narrowest of the three: it does not go through the
    # config at all. The manifest's `hyperdiffusion` block is the run's own
    # record of the damping it was prepared with, `plasim_namelist` is what the
    # model read, and both are inside the run directory, so nothing outside it
    # can move the answer and the right answer is equality. It is separate from
    # the case above because a manifest that carries the block has made the
    # claim even where the config route cannot see it: run_2b20e3324bb0 records
    # a full hyperdiffusion block and its namelist has no NDEL, NHDIFF or TDISS*
    # at all, so the model fell back to `readnl`'s compiled T21 branch -- a
    # grad^4 operator where the block says grad^8, and humidity damped 7.4x too
    # hard. world-1nz found it by hand.
    try:
        wrong, seen, silent = [], 0, []
        for manifest in manifests:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            hd = data.get("hyperdiffusion")
            nl = manifest.parent / "plasim_namelist"
            if not hd or not nl.is_file():
                continue
            layers = int((data.get("physical") or {}).get("layers")
                         or data["source_config"]["model"]["layers"])
            tau = hd["timescales_days"]
            # Per-level keys are wanted as a LIST of `layers` elements, and
            # compared through the same `namelist_elements` expansion the
            # namelist case uses. Length carries world-720: a scalar sets
            # element 1 and leaves levels 2..NLEV on readnl's compiled branch,
            # and it fails on having one element where `layers` are wanted.
            want = {"NHDIFF": [float(hd["nhdiff"])],
                    "NDEL": [float(int(hd["order_alpha"]))] * layers,
                    "TDISSD": [float(tau["divergence"])] * layers,
                    "TDISSZ": [float(tau["vorticity"])] * layers,
                    "TDISST": [float(tau["temperature"])] * layers,
                    "TDISSQ": [float(tau["humidity"])] * layers}
            seen += 1
            missing, differs = [], []
            for key, expect in want.items():
                try:
                    raw = runner.namelist_value(nl, key).rstrip(",").strip()
                except (KeyError, FileNotFoundError):
                    missing.append(key)
                    continue
                try:
                    got = runner.namelist_elements(raw)
                except ValueError:
                    differs.append(f"{key}={raw!r} is not a number")
                    continue
                if len(got) != len(expect):
                    differs.append(f"{key} has {len(got)} element(s), "
                                   f"the manifest declares {len(expect)}")
                    continue
                for i, (g, w) in enumerate(zip(got, expect)):
                    if abs(g - w) > 1e-9:
                        differs.append(f"{key}[{i + 1}]={g:g} not {w:g}")
                        break
            # One line per run: the same five keys go wrong together on
            # every run that has them wrong at all, and a per-key list
            # buries which runs are affected under how.
            if missing:
                silent.append(f"{manifest.parent.name} carries none of "
                              f"{', '.join(missing)}")
            if differs:
                wrong.append(f"{manifest.parent.name}: "
                             + ", ".join(differs))
        if wrong or silent:
            lines = silent + wrong
            rep.add(FAIL, "runs vs the hyperdiffusion they declare",
                    "; ".join(lines[:4])
                    + (f" (+{len(lines) - 4} more runs)" if len(lines) > 4 else "")
                    + " -- an absent key leaves the model on readnl's "
                      "compiled branch and a scalar leaves levels 2..NLEV "
                      "on it, so these runs damped differently from what "
                      "their manifests record")
        else:
            rep.add(OK, "runs vs the hyperdiffusion they declare",
                    f"{seen} runs declare a hyperdiffusion block and each "
                    f"carries it per level in its own namelist")
    except Exception as exc:
        rep.add(WARN, "runs vs the hyperdiffusion they declare",
                f"not checked: {exc}")

    # -- a run against the surface fields it staged --------------------------
    #
    # A `.sra` lives under one name per code per rung and is regenerated in
    # place, so an albedo rebuilt with lakes composited is still code 174. The
    # run that was made on the old one keeps reading its restart's copy while
    # the new file sits beside it in the run directory, which is why
    # `run_exoplasim.py` refuses to seed a new run from a donor whose recorded
    # `surface_field_sha256` no longer matches the file on disk. That refusal
    # only ever looked at the donor of the run being prepared. This asks it of
    # every run there is, so a regenerated input cannot orphan one silently.
    try:
        orphaned, unstamped, seen = [], [], 0
        for manifest in manifests:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            cfg = data.get("source_config")
            hashes = data.get("surface_field_sha256")
            if not cfg:
                continue
            if not hashes:
                unstamped.append(manifest.parent.name)
                continue
            seen += 1
            moved = []
            for code, was in sorted(hashes.items(), key=lambda kv: int(kv[0])):
                path = runner.surface_sra(cfg, int(code))
                if not path.is_file():
                    moved.append(f"{code} is gone")
                elif sha256_of(path) != was:
                    moved.append(f"{code} {was[:8]} -> {sha256_of(path)[:8]}")
            if moved:
                orphaned.append(f"{manifest.parent.name} "
                                f"({cfg['model']['resolution']}): "
                                + ", ".join(moved))
        if orphaned:
            rep.add(FAIL, "runs vs the surface fields they staged",
                    "; ".join(orphaned[:4])
                    + (f" (+{len(orphaned) - 4} more runs)"
                       if len(orphaned) > 4 else "")
                    + " -- the file under that code has been regenerated "
                      "since, so the run is on a surface the tree no longer "
                      "holds and nothing can be seeded from it")
        else:
            rep.add(OK, "runs vs the surface fields they staged",
                    f"{seen} runs record per-code hashes and every one still "
                    f"matches the .sra at that code")
        if unstamped:
            rep.add(WARN, "runs that record no surface hashes",
                    f"{', '.join(unstamped)} predate surface_field_sha256, "
                    f"so which surface they integrated cannot be settled")
    except Exception as exc:
        rep.add(WARN, "runs vs the surface fields they staged",
                f"not checked: {exc}")


def activation_row(build: str, terrain_hash: str | None) -> tuple[str, str, str]:
    """One report row: may `build` be named as `source_build`?

    A function rather than a paragraph inside `main()` because `--self-test`
    drives it with a build the registry is known to refuse, and with one it
    allows. A guard whose only input is the live config reports the tree rather
    than itself.
    """
    import orogen
    # By hash where there is an export to hash, because the hash is the identity;
    # by name otherwise, so the check still answers when the build is absent,
    # which is the state that makes the terrain hash optional here.
    refusal = orogen.activation_refusal(terrain_hash, name=build)
    by = "hash" if terrain_hash is not None else "name"
    if refusal:
        return (FAIL, "build activatable",
                f"lib/orogen.py refuses {build!r} as source_build: {refusal}")
    if orogen.registry_entry(terrain_hash, name=build) is None:
        return (WARN, "build activatable",
                f"{build!r} has no entry in lib/orogen.py; a build with no "
                f"registry entry cannot be refused and cannot be vouched for")
    return (OK, "build activatable", f"registry allows {build} (by {by})")


def self_test() -> int:
    """Assert what the config comparison must say about known edits.

    Each case is an edit whose correct verdict is known in advance, which is the
    only kind of case worth writing: see docs/src/practice/failure-modes.md class 17. They
    are made against the LIVE config, so they also fail if a key one of them
    names stops existing, rather than passing on a config nobody has.
    """
    import hashlib
    import re
    text = (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
    base = yaml.safe_load(text)
    inert = BIOSPHERE_INERT_CONFIG_KEYS

    def edited(pattern: str, replace) -> tuple[str, dict]:
        m = re.search(pattern, text, re.M)
        if m is None:
            raise AssertionError(f"no line matching {pattern!r} in config/planet.yaml")
        out = text[:m.start()] + replace(m) + text[m.end():]
        return out, yaml.safe_load(out)

    failures, cases = [], []

    def case(name: str, got, want, why: str):
        cases.append((name, got == want, why))
        if got != want:
            failures.append(f"{name}: expected {want!r}, got {got!r} ({why})")

    # A comment. The edit the old file-hash comparison could not distinguish
    # from a parameter change, and the reason this function exists.
    commented = "# self-test: a comment, which changes no value\n" + text
    case("a comment is not drift",
         config_drift(base, yaml.safe_load(commented), inert), [],
         "parsed values are identical")
    case("a comment does move the file hash",
         (hashlib.sha256(commented.encode()).hexdigest()
          == hashlib.sha256(text.encode()).hexdigest()), False,
         "the comparison this check replaced would have called it stale")

    # A parameter that reaches the artifact: the header's solstice offset is
    # fitted at this obliquity. THIS is the case that matters. Bumped by 1
    # rather than set to a literal, so the test carries no current value.
    _, cfg = edited(r"^(\s*obliquity_degrees:\s*)([0-9.]+)(.*)$",
                    lambda m: f"{m.group(1)}{float(m.group(2)) + 1.0}{m.group(3)}")
    case("an edited parameter is drift",
         [d.split(":")[0] for d in config_drift(base, cfg, inert)],
         ["planet.obliquity_degrees"],
         "it is read by build_vesper_header.py and sizes nothing else")

    # A key that vanishes. A generator reading a config missing a key it needs
    # raises; one reading a config that has GAINED a key does not, and both are
    # differences the artifact was not built from.
    dropped = "\n".join(line for line in text.splitlines()
                        if not re.match(r"^\s*eccentricity:", line))
    case("a deleted parameter is drift",
         [d.split(":")[0] for d in config_drift(base, yaml.safe_load(dropped), inert)],
         ["planet.eccentricity"],
         "recorded -> None is a difference like any other")

    # A DECLARED REMOVAL, both halves, over the same deletion. This is the one
    # relaxation the guard admits and it has to be shown to bite in the
    # direction that matters: a removal declared with the value the key held
    # stops being drift, and the same removal under any other value does not.
    # world-51wj. The fixture declares its own registry rather than reaching for
    # REMOVED_CONFIG_KEYS, so the case cannot pass because the real one happens
    # to be empty.
    without = yaml.safe_load(dropped)
    held = base["planet"]["eccentricity"]
    declared = {"planet.eccentricity": {
        "removed_by": "self-test", "owner": "self-test",
        "values": (held,), "moved_to": (), "moved_to_note": "self-test",
        "settles": "a fixture, not a decision about the world"}}
    case("a declared removal is not drift",
         config_drift(base, without, inert, "", declared), [],
         "the value the artifact holds is the value the declaration records")
    case("a declared removal under another value is drift",
         [d.split(":")[0] for d in config_drift(base, without, inert, "",
                                                {"planet.eccentricity":
                                                 dict(declared["planet.eccentricity"],
                                                      values=(held + 1.0,))})],
         ["planet.eccentricity"],
         "the key's value moved under the artifact before the key went")
    case("a declared removal is named where it applied",
         [line.split(",")[0] for line in
          applied_removals(base, without, declared)],
         ["planet.eccentricity"],
         "a declaration nobody sees at the moment it is used is a silence")

    # The declarations themselves, and a control for each way one can go stale.
    case("every declared removal still holds",
         removal_problems(base), [],
         "a removal for a key that is back, or that something names, excuses "
         "a refusal it has no right to")
    case("a removal for a key still in the config fails",
         bool(removal_problems(base, removed={"planet.eccentricity":
                                              declared["planet.eccentricity"]})),
         True,
         "otherwise the check above would pass on a declaration that removed nothing")
    case("a removal something still names fails",
         bool(removal_problems(without, removed={"planet.eccentricity":
                                                 declared["planet.eccentricity"]})),
         True,
         "eccentricity is named by run_exoplasim.py, so the claim is contradicted")

    # An allowlisted key, both ways round: inert only because it is allowlisted,
    # not because the edit failed to land.
    _, activity = edited(r"^(\s*activity:\s*)(\S+)(.*)$",
                         lambda m: f"{m.group(1)}{m.group(2)}-selftest{m.group(3)}")
    case("an allowlisted key is not drift",
         config_drift(base, activity, inert), [],
         "star.activity is a declaration no script reads")
    case("the same edit is drift without the allowlist",
         [d.split(":")[0] for d in config_drift(base, activity, frozenset())],
         ["star.activity"],
         "otherwise the case above would pass on an edit that never happened")

    # Every allowlist entry names a live key. `star.surface_uv` sat in the
    # resume guard's list excusing nothing, because the key is
    # `star.surface_uv_relative_to_earth`.
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    from continue_exoplasim import INERT_CONFIG_KEYS
    case("biosphere allowlist entries all exist",
         unknown_inert_keys(inert, base), [], "an entry that matches nothing excuses nothing")
    case("resume allowlist entries all exist",
         unknown_inert_keys(INERT_CONFIG_KEYS, base), [], "same, for the resume guard")

    # -- the retrospective run gates, over a fixture whose answer is known ----
    #
    # `audit_runs` reads run directories, and every run on disk today fails it.
    # That says nothing about whether it can tell a run that agrees with itself
    # from one that does not: a gate that fired at everything would look exactly
    # the same. So the fixture below is a run built to AGREE -- namelists
    # written from the same `expected_namelist_keys` the gate reads back, a
    # hyperdiffusion block derived from the same config, and surface files
    # hashed as staged -- and it is asserted GREEN first. Each case after it
    # breaks exactly one thing and names the verdict in advance.
    import tempfile
    from copy import deepcopy
    sys.path.append(str(ROOT / "exoplasim" / "scripts"))
    import run_exoplasim as runner

    NAMELISTS, HYPERDIFF, SURFACE = ("runs vs the namelists they staged",
                                     "runs vs the hyperdiffusion they declare",
                                     "runs vs the surface fields they staged")
    DERIVED = "runs vs the constants the model derived"

    def build_fixture(root: Path) -> Path:
        """A run directory that agrees with its own manifest in every way."""
        cfg = deepcopy(base)
        model = cfg["model"]
        rung = str(model["resolution"]).upper()
        run = root / "runs" / "run_selftest"
        run.mkdir(parents=True)
        want = runner.expected_namelist_keys(cfg)

        def as_fortran(value):
            """A wanted value as a namelist would carry it.

            A sequence is written out ELEMENT BY ELEMENT rather than as
            `n*value`, deliberately: the two are the same value to the model,
            so a fixture written the long way and a gate that expands
            replication must agree. If they ever stop agreeing, the green case
            below is what says so.
            """
            if isinstance(value, (list, tuple)):
                return ", ".join(f"{float(v):g}" for v in value)
            return value

        for fname, keys in want.items():
            (run / fname).write_text(
                "".join(f" {key} = {as_fortran(value)}\n"
                        for key, value in keys.items()),
                encoding="ascii")
        inputs = root / "inputs" / rung.lower()
        inputs.mkdir(parents=True)
        hashes = {}
        for code in sorted(runner.intended_surface_codes(cfg)):
            sra = inputs / f"orogen_{rung}_surf_{code:04d}.sra"
            sra.write_text(f"self-test surface field {code}\n", encoding="ascii")
            hashes[str(code)] = sha256_of(sra)
        hd = model["hyperdiffusion"]
        # A RUN THAT STARTED. The derived-constant gate reads the model's own
        # initialisation print, so a fixture with no diag would be skipped by it
        # and every case below would pass on a check that never ran. The stamp
        # is what `run_exoplasim.py` writes: for each of the three keys the
        # fixture staged at a sentinel, a positive constant the model derived.
        (run / "MOST_DIAG.00000").write_text("self-test diag\n", encoding="ascii")
        derived = {}
        for fname, key, _marker, name, units in runner.DERIVED_CONSTANT_ECHOES:
            staged = float(want[fname][key])
            derived[name] = {
                "value": 200.5476 if staged < 0 else staged,
                "branch": "derived" if staged < 0 else "declared",
                "namelist_key": f"{key}@{fname}",
                "namelist_value": staged, "units": units}
        (run / "run_manifest.json").write_text(json.dumps({
            "run_id": "run_selftest",
            "physical": {"resolution": rung, "layers": int(model["layers"])},
            "source_config": cfg,
            # Every key claimed, so a deleted one is the manifest contradicted
            # rather than a key the run predates. The case below removes a key
            # from here as well, to reach the other branch.
            "namelist_checks": {f"{k}@{f}": v
                                for f, ks in want.items() for k, v in ks.items()},
            "hyperdiffusion": {
                "rung": rung,
                "nhdiff": round(float(hd["cutoff_fraction"]) * int(rung.lstrip("T"))),
                "order_alpha": int(hd["order_alpha"]),
                "timescales_days": hd["timescales_days"][rung],
            },
            "surface_field_sha256": hashes,
            "derived_model_constants": derived,
        }, indent=1) + "\n", encoding="utf-8")
        return run

    def verdicts(edit=None) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run = build_fixture(root)
            if edit is not None:
                edit(run, root)
            saved = runner.INPUTS
            # `surface_sra` resolves against the runner's module-level INPUTS,
            # so the fixture's surface has to be reachable through it.
            runner.INPUTS = root / "inputs"
            try:
                r = Report()
                audit_runs(r, root / "runs", runner)
            finally:
                runner.INPUTS = saved
            return {check: (status, detail) for status, check, detail in r.rows}

    def drop(run: Path, fname: str, key: str) -> None:
        path = run / fname
        path.write_text("".join(
            line for line in path.read_text(encoding="ascii").splitlines(True)
            if not line.strip().upper().startswith(key + " ")), encoding="ascii")

    def unclaim(run: Path, *keys: str) -> None:
        path = run / "run_manifest.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["namelist_checks"] = {k: v for k, v in data["namelist_checks"].items()
                                   if k.split("@")[0] not in keys}
        path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")

    HD_KEYS = ("NHDIFF", "NDEL", "TDISSD", "TDISSZ", "TDISST", "TDISSQ")

    green = verdicts()
    case("a run that agrees with itself passes all four",
         {check: status for check, (status, _) in green.items()},
         {NAMELISTS: OK, HYPERDIFF: OK, SURFACE: OK, DERIVED: OK},
         "otherwise every case below would pass on a gate that only ever fails")
    case("and the fixture is a run the gates actually read",
         [green[check][1].split()[0] for check in (NAMELISTS, HYPERDIFF, SURFACE)]
         + [green[DERIVED][1].split()[0] != "0"],
         ["1", "1", "1", True],
         "an empty runs directory would report four OKs as well")

    case("a staged key with the wrong value fails",
         verdicts(lambda run, root: (run / "plasim_namelist").write_text(
             (run / "plasim_namelist").read_text(encoding="ascii")
             .replace("FILTERKAPPA = ", "FILTERKAPPA = 1"), encoding="ascii")
         )[NAMELISTS][0], FAIL,
         "the filter exponent run_2b20e3324bb0 got wrong is this shape")

    # The fixture writes arrays out element by element, so collapsing NDEL to
    # a scalar means keeping the FIRST element and dropping the rest -- which
    # is exactly what a Fortran scalar assignment to an array does, and what
    # the runs on disk actually carry.
    case("a per-level array staged as a scalar fails",
         verdicts(lambda run, root: (run / "plasim_namelist").write_text(
             re.sub(r"^ NDEL = ([^,\n]+).*$", r" NDEL = \1",
                    (run / "plasim_namelist").read_text(encoding="ascii"),
                    flags=re.MULTILINE),
             encoding="ascii"))[NAMELISTS][0], FAIL,
         "a namelist scalar sets element 1 and leaves 2..NLEV compiled; world-720")

    case("a claimed key absent from the file fails",
         verdicts(lambda run, root: drop(run, "plasim_namelist",
                                         "FILTERKAPPA"))[NAMELISTS][0], FAIL,
         "the manifest says it was staged and read back, and it is not there")

    def predates(run, root):
        drop(run, "plasim_namelist", "FILTERKAPPA")
        unclaim(run, "FILTERKAPPA")
    predated = verdicts(predates)
    case("a key the run never claimed does not fail",
         predated[NAMELISTS][0], OK,
         "the model used its compiled value; unobservable, not wrong")
    case("and it is reported rather than passed over",
         predated["namelist keys the runs never recorded"][0], WARN,
         "otherwise the case above would be indistinguishable from not looking")

    def unstamped_constant(run, root):
        path = run / "run_manifest.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["derived_model_constants"].pop("asymptotic_mixing_length_m")
        path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    case("a run that staged a sentinel and recorded no derived value fails",
         verdicts(unstamped_constant)[DERIVED][0], FAIL,
         "the namelist says the model chose the number and nothing says which")

    def wrong_constant(run, root):
        path = run / "run_manifest.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data["derived_model_constants"]["asymptotic_mixing_length_m"].update(
            {"branch": "declared", "namelist_value": 160.0, "value": 200.5476})
        path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
        (run / "fluxmod_namelist").write_text(
            (run / "fluxmod_namelist").read_text(encoding="ascii")
            .replace("VDIFF_LAMM = -1", "VDIFF_LAMM = 160"), encoding="ascii")
    case("a declared literal the model did not echo back fails",
         verdicts(wrong_constant)[DERIVED][0], FAIL,
         "the key was staged and the model integrated something else")

    def never_started(run, root):
        (run / "MOST_DIAG.00000").unlink()
    case("a run that never started is not judged on a print it has not made",
         verdicts(never_started)[DERIVED][1].split()[0], "0",
         "a prepared run has no initialisation print to read")

    def no_damping(run, root):
        for key in HD_KEYS:
            drop(run, "plasim_namelist", key)
        unclaim(run, *HD_KEYS)
    case("a declared hyperdiffusion block with no keys in the namelist fails",
         verdicts(no_damping)[HYPERDIFF][0], FAIL,
         "run_2b20e3324bb0 and run_8102b89a08ac are exactly this; world-1nz")

    case("a surface field regenerated under the same code fails",
         verdicts(lambda run, root: next(
             (root / "inputs").rglob("*.sra")).write_text("edited\n",
                                                          encoding="ascii")
         )[SURFACE][0], FAIL,
         "the name never moves, so only the digest can say the file changed")

    def unstamped(run, root):
        path = run / "run_manifest.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        data.pop("surface_field_sha256")
        path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    without = verdicts(unstamped)
    case("a run that records no surface hashes does not fail",
         without[SURFACE][0], OK, "there is no claim to contradict")
    case("and it is reported rather than passed over",
         without["runs that record no surface hashes"][0], WARN,
         "which surface it integrated cannot be settled either way")

    # The activation guard, against the registry as it stands. Driven by name
    # and not by the live config on purpose. Which build the config names moves,
    # so a case reading it asserts whichever state the tree is in on the day
    # rather than that the guard separates the two.
    import orogen
    refused = sorted(e["name"] for e in orogen._KNOWN_TERRAIN_HASHES.values()
                     if e.get("refusal"))
    allowed = sorted(e["name"] for e in orogen._KNOWN_TERRAIN_HASHES.values()
                     if not e.get("refusal"))
    # A withdrawn export and a carve verdict decided on the wrong climate are
    # permanently wrong terrains, so this list does not empty. If it ever does,
    # the positive case has no input and saying so is the only honest verdict --
    # a self-test that quietly skips its positive case is the thing it exists
    # against.
    case("the registry still refuses something",
         bool(refused), True,
         "otherwise the case below has nothing to drive it and can only pass")
    if refused:
        case("a refused build fails activation",
             activation_row(refused[0], None)[0], FAIL,
             f"{refused[0]} carries a refusal in lib/orogen.py")
    case("an allowed build passes activation",
         activation_row(allowed[0], None)[0], OK,
         f"{allowed[0]} carries none, so the case above is the refusal and "
         f"not the check failing at everything")
    case("a build in no entry is not vouched for",
         activation_row("no-such-build-selftest", None)[0], WARN,
         "unregistered is neither refused nor allowed")

    width = max(len(n) for n, _, _ in cases) + 2
    for name, ok, why in cases:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name:<{width}} {why}")
    print(f"\n{len(cases)} cases, {len(failures)} failed")
    for f in failures:
        print(f"  {f}")
    return 1 if failures else 0


# ---------------------------------------------------------------------------
# A DESIGN FLUX THAT WAS CHOSEN
#
# `derive_design_flux.py` refuses on this world: four land bands have a NEGATIVE
# warmest-bin response between the two measured flux points, so the per-band
# projection has nothing to scale by, and the guard is right. A refusal no run
# at any flux would lift leaves a design decision as the only thing that can
# settle the number, so the script records the choice INSTEAD of the derivation,
# on the one ground it declares as choosable, and only after running the
# derivation and catching the refusal.
#
# WHAT THIS GATE IS FOR, and what it deliberately is not. It cannot re-run the
# derivation -- that reads a climatology and ten orbits of model output, and it
# is the recorder's job to prove the refusal was real. What it can do is refuse
# a record that carries a number without an argument: the ground the derivation
# refused on, the finding the choice rests on, and what would reopen it are all
# required, and the ground has to be one the generator itself declares choosable
# rather than a phrase invented for the record. That list is READ from the
# generator rather than restated here, because two copies of it would be free to
# disagree and the copy in this file would be the one nobody edits.
def _check_chosen_flux(rep, design: Path, record: dict) -> None:
    """Refuse a chosen design flux that does not carry its argument."""
    label = "the design flux was CHOSEN, not derived"
    chosen = record.get("chosen")
    if not isinstance(chosen, dict):
        rep.add(FAIL, label,
                f"{rel(design)} says basis `chosen` and carries no `chosen` "
                f"block. The block is where the argument lives; without it the "
                f"artifact is a number asserting it was decided")
        return

    missing = [k for k in ("refused_at", "refusal", "reopens_when", "evidence")
               if not str(chosen.get(k) or "").strip()]
    if missing:
        rep.add(FAIL, label,
                f"{rel(design)} is missing {', '.join(missing)}. A chosen flux "
                f"records WHY the derivation refused, WHAT the choice rests on "
                f"and WHAT would reopen it; a record short of any of the three "
                f"cannot be audited by the reader who inherits the number")
        return

    scripts_dir = str(ROOT / "exoplasim" / "scripts")
    if scripts_dir not in sys.path:
        sys.path.append(scripts_dir)
    try:
        import derive_design_flux
        choosable = dict(derive_design_flux.CHOOSABLE_REFUSALS)
    except Exception as exc:
        rep.add(FAIL, label,
                f"exoplasim/scripts/derive_design_flux.py does not offer "
                f"CHOOSABLE_REFUSALS: {exc}. The grounds a chosen flux may "
                f"stand on are the generator's to declare and this check reads "
                f"them from it")
        return

    ground = str(chosen["refused_at"])
    if ground not in choosable:
        rep.add(FAIL, label,
                f"{rel(design)} says the derivation refused at {ground!r}, "
                f"which derive_design_flux.py does not declare as a ground a "
                f"choice may stand on. Its choosable refusals are "
                f"{sorted(choosable)}, and what puts a refusal in that set is "
                f"that no run at any flux would lift it")
        return

    evidence = ROOT / str(chosen["evidence"])
    if not evidence.is_file():
        rep.add(FAIL, label,
                f"{rel(design)} rests on {chosen['evidence']}, which is not on "
                f"disk. The finding a chosen flux cites has to be readable by "
                f"whoever inherits the number")
        return

    rep.add(OK, label,
            f"{record['design_flux']} chosen on {ground}, evidence "
            f"{chosen['evidence']}; reopens on {choosable[ground][:64]}...")


# THE RECORDER'S OWN REFUSALS, run against a stubbed derivation so the two ways
# a chosen flux could launder a derivation are checked rather than argued. The
# check is here for the reason the build-activation cases below are: a guard
# nothing exercises is a guard nobody knows still works, and both of these are
# cheap -- they replace the expensive half with a stub and never read a field.
def _check_chosen_flux_cannot_launder(rep) -> None:
    label = "a chosen flux cannot replace a derivation that would answer"
    scripts_dir = str(ROOT / "exoplasim" / "scripts")
    if scripts_dir not in sys.path:
        sys.path.append(scripts_dir)
    try:
        import derive_design_flux as ddf
    except Exception as exc:
        rep.add(FAIL, label, f"derive_design_flux.py does not import: {exc}")
        return

    evidence = ROOT / "notes" / "audits"
    candidates = sorted(evidence.glob("*.md"))
    if not candidates:
        rep.add(WARN, label, "no note under notes/audits/ to cite in the stub")
        return
    flux = ddf.DECLARED["anchor_flux"]

    def with_derive(stub):
        real = ddf.derive
        ddf.derive = stub
        try:
            ddf.chosen_record(flux, candidates[0], "run_stub")
        except SystemExit as exc:
            return str(exc)
        except Exception as exc:                      # pragma: no cover
            return f"raised {type(exc).__name__}: {exc}"
        finally:
            ddf.derive = real
        return None

    answered = with_derive(lambda *a, **k: {"design_flux": flux,
                                            "basis": "derived"})
    if answered is None:
        rep.add(FAIL, label,
                "chosen_record wrote a record for a derivation that ANSWERED. "
                "A flux is chosen when the instrument refuses, never when it "
                "disagrees, and this is the laundering path")
    else:
        rep.add(OK, label, "a derivation that answers is adopted, not overridden")

    def unlisted(*a, **k):
        raise ddf.Refusal("not-a-declared-ground", "stub refusal")

    unlisted_result = with_derive(unlisted)
    label2 = "a chosen flux cannot stand on an undeclared refusal"
    if unlisted_result is None:
        rep.add(FAIL, label2,
                "chosen_record wrote a record for a refusal that "
                "CHOOSABLE_REFUSALS does not list. Membership is the whole "
                "test that separates a refusal no run would lift from one that "
                "is fixed by supplying an input")
    else:
        rep.add(OK, label2,
                f"{len(ddf.CHOOSABLE_REFUSALS)} declared ground(s); anything "
                f"else refuses")


_MISSING = object()


def _config_node(config: dict, dotted: str):
    """The value at a dotted path in the configuration, or `_MISSING`."""
    node = config
    for step in dotted.split("."):
        if not isinstance(node, dict) or step not in node:
            return _MISSING
        node = node[step]
    return node


def _written_tolerance(value) -> float:
    """Half a unit in the last place the declaration actually wrote.

    Arithmetic on the written precision, not a number picked to make today's
    values pass: a declaration carrying three decimals cannot honestly differ
    from its generator by more than half of the third.

    Exponent form is read the same way rather than as an integer: `1.719635e+6`
    carries six decimals on a mantissa scaled by a million, so its last place is
    a whole unit and not a tenth. A YAML file that writes a large derived
    quantity that way would otherwise get a tolerance a million times tighter
    than what it wrote.
    """
    written = str(value)
    mantissa, _, exponent = written.lower().partition("e")
    places = len(mantissa.split(".")[1]) if "." in mantissa else 0
    shift = int(exponent) if exponent else 0
    return 0.5 * 10.0 ** (shift - places)


def check_derived_config_values(rep: "Report", config: dict) -> None:
    """Every DERIVED literal in `config/planet.yaml` against the artifact behind it.

    Each row is a value some generator computes, prints, and a human retypes
    into the configuration. That route has no check on it, and it has already
    drifted repeatedly -- the CO2 shortwave fit in the model source, the dust
    notes against the albedos the code now reads, and
    `model.land_longwave_emissivity` against a build the registry refuses.
    `notes/audits/frozen-derived-quantities.md` is the enumeration.

    Compared to the FULL-PRECISION value in the artifact at the tolerance the
    configuration's own rounding implies, element by element for the rows that
    are lists: a bracket is two numbers and either of them can go stale on its
    own.

    A row may also require the artifact to be STAMPED with the configured
    build. That is not decoration on those rows: `land_longwave_emissivity` is
    DEFINED as an area weighting over the build's land, so every terrain change
    moves it, and a row that compared only numbers would pass an artifact
    computed on a build this tree no longer points at -- which is exactly how
    that key came to carry another terrain's mean.

    The generator is named in each row so a failure says what to re-run.
    """
    try:
        analysis = ROOT / "analysis"
        sw = ROOT / "exoplasim" / "analysis" / "shortwave_band_weights.json"
        cloud = ROOT / "exoplasim" / "analysis" / "cloud_band_weight.json"
        # `h2o_sw_level` is read out of its OWN artifact and not out of
        # `shortwave_band_weights.json`, although that report carries the same
        # node. The full report cannot be written while `baseline_climatology`
        # is null -- its prediction blocks open a climatology and have no Earth
        # fallback by design -- so a row pointed at it would turn a stale
        # derived artifact into a failing gate. `--level` writes the artifact
        # below from two declared constants and needs no climatology, no
        # spectrum and no network, so this row can be regenerated whenever it
        # fails. It catches a retyping error into config and it CANNOT catch
        # the climatology moving under `CORRK_PATH_CM`, the water path the
        # ratio it is built from was evaluated at; `check_water_path_currency`
        # below is what asks that question.
        level = ROOT / "exoplasim" / "analysis" / "h2o_sw_level.json"
        veg = analysis / "vegetation_albedo.json"
        ice = analysis / "ice_properties.json"
        rows = [
            ("model.h2o_sw_weight", sw, ("h2o", "weight"), False,
             "exoplasim/scripts/shortwave_band_weights.py"),
            ("model.co2_sw_weight", sw, ("co2", "weight"), False,
             "exoplasim/scripts/shortwave_band_weights.py"),
            ("model.h2o_sw_level", level, ("value",), False,
             "exoplasim/scripts/shortwave_band_weights.py --level"),
            ("model.cloud_absorption_scale", cloud, ("weight", "central"), False,
             "exoplasim/scripts/cloud_band_weight.py"),
            ("model.vegetation_albedo", veg, ("vegetation_albedo",), False,
             "analysis/vegetation_albedo.py"),
            ("model.vegetation_albedo_bracket", veg,
             ("vegetation_albedo_bracket",), False, "analysis/vegetation_albedo.py"),
            ("model.vegetation_albedo_bands", veg,
             ("vegetation_albedo_bands",), False, "analysis/vegetation_albedo.py"),
            ("model.tree_albedo", veg, ("cover_albedo", "tree"), False,
             "analysis/vegetation_albedo.py"),
            ("model.tree_albedo_bracket", veg, ("cover_albedo_bracket", "tree"),
             False, "analysis/vegetation_albedo.py"),
            ("model.grass_albedo", veg, ("cover_albedo", "grass"), False,
             "analysis/vegetation_albedo.py"),
            ("model.grass_albedo_bracket", veg, ("cover_albedo_bracket", "grass"),
             False, "analysis/vegetation_albedo.py"),
            ("surface.cryosphere.snow_fusion_j_kg", ice,
             ("pure_ice_ih", "melting_enthalpy_j_kg"), False,
             "analysis/ice_properties.py"),
            ("model.lithology_albedo_overrides.playa_clastic.albedo",
             analysis / "playa_albedo.json", ("reconciled_albedo",), False,
             "analysis/playa_albedo.py"),
            # STAMPED, because this one is an area weighting over the build's
            # land and moves with every terrain.
            ("model.land_longwave_emissivity", analysis / "land_emissivity.json",
             ("land_mean",), True, "analysis/emissivity_contrast.py --land-mean"),
        ]
        import provenance
        build = provenance.active_build(config)
        want_hash = builds.terrain_hash(config)
        problems, checked = [], []
        for key, path, where, stamped, generator in rows:
            declared = _config_node(config, key)
            if declared is _MISSING:
                problems.append(f"{key} is not in config/planet.yaml")
                continue
            if not path.is_file():
                problems.append(f"{key}: {rel(path)} is missing; run {generator}")
                continue
            report = json.loads(path.read_text(encoding="utf-8"))
            # Walked with a guard rather than bare indexing: a missing node
            # would otherwise raise into this block's own except and downgrade
            # every OTHER comparison here to a warning, which is the one failure
            # mode a consistency check must not have.
            node, missing = report, False
            for step in where:
                if not isinstance(node, dict) or step not in node:
                    problems.append(
                        f"{key}: {rel(path)} has no {'.'.join(where)}; it "
                        f"predates the node, so re-run {generator}")
                    missing = True
                    break
                node = node[step]
            if missing:
                continue
            if stamped:
                got = report.get("terrain_hash")
                if got != want_hash:
                    problems.append(
                        f"{key}: {rel(path)} carries terrain "
                        f"{str(got)[:8]} and the configured build {build} is "
                        f"{str(want_hash)[:8]}. This key is an area weighting "
                        f"over the build's land, so the artifact describes "
                        f"another terrain: re-run {generator}")
                    continue
            declared_list = isinstance(declared, (list, tuple))
            node_list = isinstance(node, (list, tuple))
            if declared_list != node_list:
                problems.append(
                    f"{key} is {'a list' if declared_list else 'a scalar'} in "
                    f"config and {'a list' if node_list else 'a scalar'} in "
                    f"{rel(path)}; the two do not describe the same quantity")
                continue
            pairs = (list(zip(declared, node)) if declared_list
                     else [(declared, node)])
            if declared_list and len(declared) != len(node):
                problems.append(
                    f"{key} has {len(declared)} entries in config and "
                    f"{len(node)} in {rel(path)}: re-run {generator}")
                continue
            bad = []
            for index, (mine, theirs) in enumerate(pairs):
                tol = _written_tolerance(mine)
                if abs(float(mine) - float(theirs)) > tol:
                    where_at = f"[{index}]" if declared_list else ""
                    bad.append(f"{where_at}{float(mine):g} against "
                               f"{float(theirs):.6g} (tolerance {tol:g})")
            if bad:
                problems.append(
                    f"{key}: {'; '.join(bad)} in {rel(path)} "
                    f"({'.'.join(where)}), past what its own written precision "
                    f"allows: re-run {generator} or retype it")
            else:
                checked.append(key.split(".", 1)[-1])
        rep.add(FAIL if problems else OK,
                "derived config values match their artifacts",
                "; ".join(problems) if problems else
                f"{len(checked)} compared: {', '.join(checked)}")
    except Exception as exc:
        rep.add(WARN, "derived config values match their artifacts",
                f"not checked: {exc}")


def check_derived_config_rules(rep: "Report", config: dict) -> None:
    """Derived config literals whose producer is a RULE rather than an artifact.

    Same class as `check_derived_config_values` and the same two dispositions
    behind it; what differs is that nothing writes these to disk, so the
    re-derivation is the rule itself, stated once, here. Each is compared at the
    tolerance the configuration's own written precision implies.
    """
    import lapse
    import stellar

    model = config.get("model", {})
    problems, checked = [], []

    # `semi_implicit_reference_temperature_k`. The mass-weighted mean of the
    # `setzt` profile over the model's own sigma levels, weights dsigma, which
    # is what `plasim.f90` builds T0 from. The rule is what the config comment
    # states; running it here is what stops the core's adiabatic conversion
    # defect resting on a number nobody re-derives.
    try:
        declared = float(model["semi_implicit_reference_temperature_k"])
        derived = lapse.semi_implicit_reference_temperature_k(config)
        tol = _written_tolerance(model["semi_implicit_reference_temperature_k"])
        if abs(declared - derived) > tol:
            problems.append(
                f"semi_implicit_reference_temperature_k is {declared:g} in "
                f"config and the dsigma-weighted setzt profile gives "
                f"{derived:.4f}, past the {tol:g} its own precision allows")
        else:
            checked.append(f"t0 {declared:g} K")
    except Exception as exc:                       # noqa: BLE001 - reported
        problems.append(f"semi_implicit_reference_temperature_k not checked: {exc}")

    # `ozone_visible_weight`. The Chappuis band's share of this star's flux over
    # its share of the Sun's, on the spectrum file the model reads. The spectrum
    # is rewritten in place by `build_stellar_spectrum.py`, and this weight went
    # stale the last time it was: a resampling that moved the band-1 share moved
    # this with it and nothing said so.
    try:
        declared = float(model["ozone_visible_weight"])
        derived = stellar.ozone_visible_weight()
        tol = _written_tolerance(model["ozone_visible_weight"])
        if abs(declared - derived) > tol:
            problems.append(
                f"ozone_visible_weight is {declared:g} in config and the "
                f"configured spectrum gives {derived:.4f} over the Chappuis "
                f"band, past the {tol:g} its own precision allows: re-derive "
                f"it with lib/stellar.py:ozone_visible_weight")
        else:
            checked.append(f"o3visw {declared:g}")
    except Exception as exc:                       # noqa: BLE001 - reported
        problems.append(f"ozone_visible_weight not checked: {exc}")

    rep.add(FAIL if problems else OK,
            "derived config values match the rules that produce them",
            "; ".join(problems) if problems else
            f"{len(checked)} re-derived: {', '.join(checked)}")


def _written_token(path: Path, key: str) -> str | None:
    """The literal a YAML file wrote for `key`, so its precision can be read.

    `str(value)` cannot answer this: YAML parses `1.719635e+6` and `0.9350` into
    floats whose `repr` has lost the form they were written in, and the form is
    exactly what says how far a declaration may honestly sit from its generator.
    """
    import re
    match = re.search(rf"^\s*{re.escape(key)}\s*:\s*([-+0-9.eE]+)\s*$",
                      path.read_text(encoding="utf-8"), re.MULTILINE)
    return match.group(1) if match else None


def check_land_column_thermal_constants(rep: "Report") -> None:
    """The land column contract's four cryosphere thermal constants against
    `analysis/ice_properties.py`, which is what produces them.

    All four are properties of a DENSITY the model compiles -- the glacier pair
    of `rhoglac`, the snow pair of `rhosnow` -- and `ice_properties.py` already
    holds the compiled Fortran literals to what it computes. It never opens this
    yaml, so the contract's copy was a fourth restatement sitting outside that
    loop: a bracket on either density would have moved the model and the
    derivation together and left the contract behind.
    """
    contract = ROOT / "pedology" / "config" / "land_column_properties.yaml"
    artifact = ROOT / "analysis" / "ice_properties.json"
    rows = [
        ("glacier_ice_heat_capacity_j_m3_k", ("glacial_ice", "heat_capacity_j_m3_k")),
        ("glacier_ice_thermal_conductivity_w_m_k",
         ("glacial_ice", "conductivity_w_m_k")),
        ("snow_heat_capacity_j_m3_k",
         ("snow_conductivity_bracket", "declared_heat_capacity_j_m3_k")),
        ("snow_thermal_conductivity_w_m_k",
         ("snow_conductivity_bracket", "declared_conductivity_w_m_k")),
    ]
    label = "the land column contract's cryosphere thermal constants"
    try:
        if not artifact.is_file():
            rep.add(FAIL, label,
                    f"{rel(artifact)} is missing; run analysis/ice_properties.py")
            return
        declared = (yaml.safe_load(contract.read_text(encoding="utf-8"))
                    ["thermal"]["constants"])
        report = json.loads(artifact.read_text(encoding="utf-8"))
        problems, checked = [], []
        for key, where in rows:
            if key not in declared:
                problems.append(f"{key} is not in {rel(contract)}")
                continue
            node = report
            for step in where:
                if not isinstance(node, dict) or step not in node:
                    node = None
                    break
                node = node[step]
            if node is None:
                problems.append(
                    f"{key}: {rel(artifact)} has no {'.'.join(where)}; it "
                    f"predates the node, so re-run analysis/ice_properties.py")
                continue
            token = _written_token(contract, key) or str(declared[key])
            tol = _written_tolerance(token)
            if abs(float(declared[key]) - float(node)) > tol:
                problems.append(
                    f"{key} is {token} in {rel(contract)} and {float(node):.7g} "
                    f"in {rel(artifact)}, past the {tol:g} its own written "
                    f"precision allows: re-run analysis/ice_properties.py and "
                    f"retype it")
            else:
                checked.append(key)
        rep.add(FAIL if problems else OK, label,
                "; ".join(problems) if problems else
                f"{len(checked)} compared against {rel(artifact)}")
    except Exception as exc:                       # noqa: BLE001 - reported
        rep.add(WARN, label, f"not checked: {exc}")


# The run `model.cold_start_profile.surface_temperature_k` is measured on. Named
# here rather than in `config/planet.yaml` for the reason `lib/sensitivity.py`
# names its bracket runs in code: a run id inside a YAML comment is reachable by
# a reader and by nothing else, and what makes this a check is that something
# opens the run index.
COLD_START_RUN = "run_432e5e46adef"


def check_cold_start_currency(rep: "Report", config: dict) -> None:
    """The cold start's surface temperature against the run it is taken from.

    `model.cold_start_profile.surface_temperature_k` is this model's own
    converged global-mean surface temperature at the declared baseline flux, so
    it is a MEASUREMENT and it can go stale in two ways at once: the number can
    move, and the run it was taken from can leave the tree. Both have happened
    to it. The check therefore asks the same two questions
    `lib/sensitivity.py:verify` asks of the slope -- does the declaration match
    the run it names, and is that run still on the world this tree describes --
    and it names the run here rather than in the configuration because a run id
    in a YAML comment is not reachable by anything.
    """
    label = "the cold start's surface temperature is current"
    try:
        import sensitivity
        want_build = config["source_build"]
        want_flux = float(config["orbit"]["baseline_flux_earth"])
        declared = config["model"]["cold_start_profile"]["surface_temperature_k"]
        live = sensitivity._live_entries()
        archived = sensitivity._archived_entries()
        entry = live.get(COLD_START_RUN)
        problems = []
        if entry is None:
            where = ("survives only as archived identity, on build "
                     f"{archived[COLD_START_RUN].get('source_build')}"
                     if COLD_START_RUN in archived else "is in no run index at all")
            rep.add(FAIL, label,
                    f"{COLD_START_RUN} {where}; the temperature it carries "
                    f"cannot be recomputed, so re-derive this declaration "
                    f"against a live run on {want_build}")
            return
        on = entry.get("source_build")
        if on != want_build:
            problems.append(f"{COLD_START_RUN} was run on build {on}, not the "
                            f"configured {want_build}")
        flux = float(entry.get("physical", {}).get("flux_ratio", float("nan")))
        if not math.isclose(flux, want_flux, abs_tol=1e-9):
            problems.append(f"{COLD_START_RUN} is at flux ratio {flux}, not the "
                            f"declared orbit.baseline_flux_earth {want_flux}")
        report = sensitivity._report_metrics(f"{COLD_START_RUN}_convergence.json")
        if report is None:
            problems.append(f"{COLD_START_RUN} has no convergence report")
        else:
            if not report.get("sufficiently_equilibrated_for_worldbuilding"):
                problems.append(
                    f"{COLD_START_RUN} no longer passes every convergence "
                    f"criterion: {report.get('failed_criteria')}")
            asymptote = (report.get("metrics") or {}).get("temperature_asymptote_k")
            if asymptote is None:
                problems.append(f"{COLD_START_RUN} has no fitted asymptote")
            else:
                tol = _written_tolerance(declared)
                if abs(float(declared) - float(asymptote)) > tol:
                    problems.append(
                        f"the declared {float(declared):g} K is not "
                        f"{COLD_START_RUN}'s fitted asymptote {asymptote:.4f} K, "
                        f"past the {tol:g} its own precision allows")
        rep.add(FAIL if problems else OK, label,
                "; ".join(problems) if problems else
                f"{float(declared):g} K is {COLD_START_RUN}'s fitted asymptote, "
                f"on {want_build} at flux ratio {want_flux:g}")
    except Exception as exc:                       # noqa: BLE001 - reported
        rep.add(WARN, label, f"not checked: {exc}")


def check_eddy_wind(rep: "Report", config: dict) -> None:
    """`model.hyperdiffusion.eddy_wind_m_s` against every run that measured it.

    THE PRODUCER IS A SET AND NOT AN ARTIFACT, which is what makes this row
    different from the ones above. `measure_eddy_wind.py` writes one verdict per
    run under `exoplasim/analysis/eddy-wind/`, and the declared wind is an
    aggregate over them with a declared bracket around it. So the honest check
    is not equality with any one file: it is that the declaration still sits
    inside the spread of the measurements that exist, and that the set has not
    emptied out from under it.

    The bracket is the configuration's own and is not recomputed here. A
    tolerance derived from today's files would move whenever a run was added,
    which is a criterion chosen after the result it judges.

    THE BUILD IS REPORTED AND IS NOT A REFUSAL. Every measurement on disk is on
    a build that is not the configured one, and the configuration says in full
    why that is not a defect to fix by re-running: these runs were damped with
    the PREVIOUS value, so the sweep that closes the bracket is the next run's
    own re-measurement rather than a re-measurement of these.
    """
    label = "the declared eddy wind sits inside its own measurements"
    try:
        declared = float(config["model"]["hyperdiffusion"]["eddy_wind_m_s"])
        bracket = 0.05          # the +/- 5 per cent the configuration declares
        directory = ROOT / "exoplasim" / "analysis" / "eddy-wind"
        files = sorted(directory.glob("*.json")) if directory.is_dir() else []
        if not files:
            rep.add(FAIL, label,
                    f"{rel(directory)} holds no measurement at all, so the "
                    f"declared {declared:g} m/s has nothing behind it; run "
                    f"exoplasim/scripts/measure_eddy_wind.py")
            return
        measured, builds_seen, unreadable = {}, set(), []
        for path in files:
            try:
                report = json.loads(path.read_text(encoding="utf-8"))
                measured[path.stem] = float(
                    report["measured"]["eddy_wind_m_s"])
                builds_seen.add(report.get("source_build")
                                or (report.get("run") or {}).get("source_build"))
            except Exception:                      # noqa: BLE001 - reported
                unreadable.append(path.name)
        if not measured:
            rep.add(FAIL, label,
                    f"none of the {len(files)} file(s) in {rel(directory)} "
                    f"carries measured.eddy_wind_m_s: {', '.join(unreadable)}")
            return
        low, high = min(measured.values()), max(measured.values())
        problems = []
        if not (declared * (1 - bracket) <= low
                and high <= declared * (1 + bracket)):
            problems.append(
                f"the {len(measured)} measurement(s) span {low:.4g} to "
                f"{high:.4g} m/s and the declared {declared:g} carries a "
                f"+/-{bracket:.0%} bracket, which does not cover them: "
                f"re-derive the aggregate")
        if unreadable:
            problems.append(f"unreadable: {', '.join(unreadable)}")
        note = ""
        stale = {b for b in builds_seen if b and b != config["source_build"]}
        if stale:
            note = (f"; every measurement is on {', '.join(sorted(stale))} and "
                    f"the configured build is {config['source_build']}, which "
                    f"the declaration accounts for: the sweep that closes the "
                    f"bracket is the next run's own re-measurement")
        rep.add(FAIL if problems else OK, label,
                "; ".join(problems) if problems else
                f"{len(measured)} measurement(s) span {low:.4g} to {high:.4g} "
                f"m/s, inside {declared:g} +/-{bracket:.0%}{note}")
    except Exception as exc:                       # noqa: BLE001 - reported
        rep.add(WARN, label, f"not checked: {exc}")


def check_water_path_currency(rep: "Report", config: dict) -> None:
    """`CORRK_PATH_CM` against the water path the best available climatology has.

    THE SHORTWAVE GATE ABOVE CANNOT ASK THIS. `h2o_sw_level.json` is written by
    `shortwave_band_weights.py --level` from `CORRK_RATIO_TO_EQ21` and
    `H2O_CONTINUUM_FRACTION` and from nothing else, so the row that holds
    `model.h2o_sw_level` against it regenerates its own reference from the
    frozen value: it catches a retyping error into the configuration and it
    cannot catch the climatology moving under the PATH the ratio was evaluated
    at. A check that regenerates its own reference from the frozen value is not
    a check on that value.

    This is what asks. `CORRK_PATH_CM` is the model's own effective water path,
    magnified, and it was measured on a baseline climatology that no longer
    exists; the best available climatology in this tree is what can be measured
    now. Reported rather than refused, and the reason is measured rather than
    preferred: `h2o_sw_weight` and `co2_sw_weight` are quoted at the SAME path
    and their report cannot be regenerated while `baseline_climatology` is null,
    so moving this one alone would leave the three quoted at two different
    paths. What settles it is a baseline climatology on the configured build,
    at which point all three move together. world-wtt3.
    """
    label = "the shortwave water path against the best available climatology"
    try:
        sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
        import shortwave_band_weights as swbw
        from paths import best_available_climatology

        best = best_available_climatology()
        if not best.path.is_file():
            rep.add(WARN, label, f"{rel(best.path)} is not on disk")
            return
        with Dataset(best.path) as data:
            hus = np.asarray(data.variables["hus"][:])
            air_t = np.asarray(data.variables["ta"][:])
            surface_p = np.asarray(data.variables["ps"][:]) * 100.0
            sigma = np.asarray(data.variables["lev"][:])
            sigma_half = np.asarray(data.variables["levp"][:])
            lat = np.asarray(data.variables["lat"][:])
        gravity = float(config["planet"]["gravity_m_s2"])
        dsigma = np.diff(sigma_half)
        column = np.zeros_like(surface_p)
        for k in range(len(sigma)):
            column += (0.1 * dsigma[k] * hus[:, k] * surface_p / gravity
                       * np.sqrt(273.0 / air_t[:, k])
                       * sigma[k] * surface_p / 1.0e5)
        weights = np.cos(np.deg2rad(lat))[None, :, None] * np.ones_like(column)
        here = (float((column * weights).sum() / weights.sum())
                * swbw.WATER_MAGNIFICATION)
        declared = swbw.CORRK_PATH_CM
        rep.add(WARN if abs(here - declared) > 0.05 * declared else OK, label,
                f"declared {declared:g} cm against {here:.4f} on "
                f"{rel(best.path)} ({best.stage}). The declaration was measured "
                f"on a baseline that no longer exists; it is not moved alone "
                f"because h2o_sw_weight and co2_sw_weight are quoted at the "
                f"same path and cannot be regenerated without a baseline. "
                f"world-wtt3")
    except Exception as exc:                       # noqa: BLE001 - reported
        rep.add(WARN, label, f"not checked: {exc}")


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Verify that everything in tree describes the same world. "
                    "Exit 1 on disagreement. Takes no arguments; it audits the "
                    "artifacts as they are.")
    parser.add_argument(
        "--self-test", action="store_true",
        help="check the checker instead of the tree: assert what the config "
             "comparison says about a comment, an allowlisted key, an edited "
             "parameter and a deleted one. Touches nothing on disk.")
    if parser.parse_args().self_test:
        return self_test()
    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    rep = Report()
    build = str(config.get("source_build", ""))
    # A build being ABSENT is a legitimate state, not a crash. One is disposable
    # until a climate run has consumed it, so between a generator change and the
    # next generation there is deliberately nothing in source/. Report it and run
    # the checks that do not need a build, rather than exiting before the summary
    # and leaving the caller with a traceback instead of a work list.
    want = None
    try:
        want = builds.terrain_hash(config)
    except Exception as exc:
        rep.add(FAIL, "active build",
                f"{build!r} is not present: {exc}. Generate it -- the recipe is "
                f"in source/README.md -- then register its hash in lib/orogen.py")
    if want is not None:
        rep.add(OK, f"active build", f"{build}  {want[:16]}")

    # -- the registry does not refuse the build the config names -------------
    #
    # `lib/orogen.py` is the registry, and being registered there is not
    # permission to be active. Some entries are recorded precisely so that a
    # result already computed from them stays readable: a withdrawn export, three
    # superseded carve verdicts, and a 10M-region measurement artifact generated
    # for the resolution audit. Each of those carries a `refusal`, and a build
    # with one may not be named as `source_build`.
    #
    # Nothing enforced that, and the two files disagreed in silence with config
    # naming a build whose entry refused it. What that costs is not abstract.
    # Batch P0C could not tell which build was
    # active and measured the basin catalogue on both, and relaxing minCells adds
    # 4,463 basins on precarve-craton against 8 on precarve-craton-10m, so
    # "is the catalogue resolution-limited" came back TRUE on one build and FALSE
    # on the other.
    #
    # Here rather than in smoke_test.py, though the check itself is a static read
    # and would cost that tier nothing. The tier is set by the cost of the STATE,
    # not of the check. A refused build named in config is a legitimate transient
    # -- the export being regenerated, the config still pointing at the old one --
    # and a per-commit gate that fails on it stops every commit in the tree,
    # including the commits that resolve it. What must never happen is that the
    # disagreement survives into something that MEASURES, and this file is the
    # gate that runs before an expensive run and after a change to `source_build`.
    if build:
        rep.add(*activation_row(build, want))

    # -- both named climatologies belong to the active build -----------------
    #
    # A climatology is the most widely shared artifact here: pedology,
    # hydrography and the biosphere all read one. It has the same grid, the
    # same variables and the same units whichever world it describes, so a
    # superseded one produces a plausible number rather than an error. The
    # consumers now check at the point of reading; this is the same check
    # applied once, ahead of an expensive run.
    #
    # BOTH KEYS, because there are two climatologies and the graph asks for one
    # or the other per step. The bootstrap is the run on terrain-only surface
    # fields and is what a step's `needs` edge requires to EXIST; the baseline
    # is the run on those fields once they exist, and a state-dependent step
    # READS it once it is named, through
    # `lib/paths.py:best_available_climatology`. No list of steps here: the
    # register is `config/pipeline.yaml` and a copy of it drifts. A stale one of
    # either climatology is the same defect and neither is checked by the
    # other.
    sys.path.insert(0, str(ROOT / "lib"))
    from provenance import artifact_build
    for key, label in (("bootstrap_climatology", "bootstrap climatology"),
                       ("baseline_climatology", "baseline climatology")):
        declared = config.get(key)
        if not declared:
            rep.add(WARN, label,
                    f"none named in config; anything needing {key} will raise")
            continue
        clim = ROOT / declared
        if not clim.is_file():
            rep.add(FAIL, label,
                    f"config names {declared}, which does not exist")
            continue
        got = artifact_build(clim)
        if got is None:
            rep.add(WARN, label,
                    f"{Path(declared).name} carries no build identity; it "
                    f"predates the stamping in build_climatology.py")
        elif got != build:
            rep.add(FAIL, label,
                    f"{Path(declared).name} is on {got}, config names {build}")
        else:
            rep.add(OK, label, f"on {build}")

    # -- the configured soil saturation endpoints against the contract -------
    #
    # `pedology/config/land_column_properties.yaml` DERIVES the degree of
    # saturation at an empty and a full store as medians over the build's own
    # land cells, and `pedology/scripts/land_column_properties.py` refuses when
    # they stop being that statistic. A run reaches those numbers through copies:
    # `landmod.f90` carries them as compiled `landmod_nl` defaults, and
    # `config/planet.yaml` restates the surface layer's pair because the
    # moisture-dependent soil albedo is an ARM and an arm's endpoints have to be
    # declared where the arm is. `run_exoplasim.py` refuses at import when the
    # compiled defaults drift from the contract. Nothing refused when the CONFIG
    # did, and that is the copy an armed term reads.
    #
    # WHY THIS IS THE TIER FOR IT. The medians move with every soil rebuild,
    # which is every turn of loop A, so this is exactly a "do the artifacts
    # agree" question about a derived quantity, and what it protects is the run
    # that comes after it: a stale upper endpoint means the modelled albedo
    # mixing reaches its wet end at a saturation the soil no longer has.
    #
    # AN UNSET KEY IS NOT A FAILURE. `run_exoplasim.py` writes each of these
    # keys into `landmod_nl` whatever the config says; a key the config leaves
    # out is written at `landmod.f90`'s compiled default, which that script
    # already holds to the contract. Unset is one copy of the number instead of
    # two.
    contract_path = (ROOT / "pedology" / "config"
                     / "land_column_properties.yaml")
    label = "configured soil saturation vs the land column contract"
    if not contract_path.is_file():
        rep.add(FAIL, label, f"{rel(contract_path)} does not exist")
    else:
        contract = yaml.safe_load(contract_path.read_text(encoding="utf-8"))
        mapping = contract["thermal"]["saturation_mapping"]
        tol = float(mapping["tolerance"])
        surface_map = contract["surface_layer"]["saturation_mapping"]
        node = contract
        for part in str(surface_map["sr_at_field_capacity_from"]).split("."):
            node = node[part]
        wanted = {
            ("soil_albedo_moisture", "saturation_at_empty_layer"):
                (float(surface_map["sr_at_air_dry"]),
                 "surface_layer.saturation_mapping.sr_at_air_dry"),
            ("soil_albedo_moisture", "saturation_at_full_layer"):
                (float(node), str(surface_map["sr_at_field_capacity_from"])),
            ("soil_thermal", "saturation_at_empty_store"):
                (float(mapping["sr_at_wilting_point"]),
                 "thermal.saturation_mapping.sr_at_wilting_point"),
            ("soil_thermal", "saturation_at_full_store"):
                (float(mapping["sr_at_field_capacity"]),
                 "thermal.saturation_mapping.sr_at_field_capacity"),
        }
        surface_cfg = config.get("surface", {}) or {}
        agreed, unset = [], []
        for (block_name, key), (declared, source) in sorted(wanted.items()):
            stated = (surface_cfg.get(block_name, {}) or {}).get(key)
            if stated is None:
                unset.append(f"surface.{block_name}.{key}")
            elif abs(float(stated) - declared) > tol:
                rep.add(FAIL, label,
                        f"config/planet.yaml surface.{block_name}.{key} is "
                        f"{float(stated):g} and {rel(contract_path)} derives "
                        f"{declared:g} at {source}, outside the contract's own "
                        f"tolerance of {tol:g}. The endpoint is a median over "
                        f"this build's land cells and the soil has been rebuilt "
                        f"under it; re-derive with "
                        f"`pedology/scripts/land_column_properties.py "
                        f"--update-declaration` and carry the new value into "
                        f"every restatement it names")
            else:
                agreed.append(f"surface.{block_name}.{key}")
        if agreed:
            rep.add(OK, label, f"{len(agreed)} declared and agreeing: "
                               f"{', '.join(agreed)}")
        if unset:
            rep.add(OK, label, f"{len(unset)} unset, so written at "
                               f"landmod.f90's compiled default, which "
                               f"run_exoplasim.py holds to the same contract: "
                               f"{', '.join(unset)}")

    # -- the configured flux is the DERIVED flux -----------------------------
    #
    # `config/pipeline.yaml` makes `baseline_run` need `design_flux`, and that
    # edge only says the ARTIFACT exists. Adopting its winner into
    # `config/planet.yaml` is a separate human act, because the flux fixes the
    # semi-major axis and that orbit is compiled into LPJ-GUESS. So the edge
    # alone leaves the whole failure available: derive a flux, do not adopt it,
    # buy hours of baseline run at the old number, and every artifact below it
    # describes a world nothing intended.
    #
    # WHY THIS IS THE TIER FOR IT. It is a static read of two files, so it could
    # sit in `smoke_test.py`; it reads a GENERATED artifact, so it belongs where
    # the artifact comparisons are, and this is what runs before an expensive
    # run, which is the thing being protected.
    #
    # NO ARTIFACT IS NOT A FAILURE. A provisional flux with no design_flux.json
    # on disk is the legitimate current state of this tree: the derivation has
    # not been run on this terrain yet, and refusing it would refuse the
    # commissioning that produces the artifact. A design_flux.json whose winner
    # DIFFERS from the configured value is the state to refuse, and the window
    # between writing the artifact and adopting its number is one edit rather
    # than a resting state. An exploratory derivation writes elsewhere with
    # --output; the registered path carries the answer.
    design = ROOT / "exoplasim" / "analysis" / "design_flux.json"
    configured = float(config["orbit"]["baseline_flux_earth"])
    if not design.is_file():
        rep.add(WARN, "configured flux vs the derived one",
                f"{configured} is PROVISIONAL: {rel(design)} does not exist, so "
                f"no derivation has been run on this terrain. Loop C re-derives "
                f"the flux on every new terrain")
    else:
        try:
            record = json.loads(design.read_text(encoding="utf-8"))
            derived_flux = float(record["design_flux"])
        except Exception as exc:
            rep.add(FAIL, "configured flux vs the derived one",
                    f"{rel(design)} does not carry a readable `design_flux`: {exc}")
        else:
            # DERIVED OR CHOSEN, AND THE ARTIFACT HAS TO SAY WHICH. A flux the
            # derivation returned and a flux a person picked because the
            # derivation refused are two different kinds of number, and every
            # consumer that reads this file is entitled to know which one it is
            # holding without going and reading a note. An artifact with no
            # `basis` is one written before this distinction existed or by hand,
            # and both are refused rather than assumed to be the derived case.
            basis = record.get("basis")
            if basis not in ("derived", "chosen"):
                rep.add(FAIL, "configured flux vs the derived one",
                        f"{rel(design)} carries basis {basis!r}. Every artifact "
                        f"derive_design_flux.py writes says `derived` or "
                        f"`chosen`; a number without one of those is a flux "
                        f"whose standing nothing reading it can establish")
            elif basis == "chosen":
                _check_chosen_flux(rep, design, record)

            if abs(derived_flux - configured) < 1e-9:
                rep.add(OK, "configured flux vs the derived one",
                        f"{configured}, adopted, basis {basis}")
            else:
                rep.add(FAIL, "configured flux vs the derived one",
                        f"{rel(design)} carries {derived_flux} (basis {basis}) "
                        f"and config/planet.yaml runs at {configured}. Adopt "
                        f"it, or re-derive on this terrain; a baseline run "
                        f"bought at the unadopted number is hours spent on a "
                        f"world nothing intended, and the flux fixes a "
                        f"semi-major axis compiled into LPJ-GUESS")

    _check_chosen_flux_cannot_launder(rep)

    # -- the stellar band split, and everything carrying a copy of it -------
    #
    # `zsolar1` weights the two-band snow, sea-ice, glacier and ground albedos.
    # It had three values in three artifacts and a fourth in the model, and the
    # model's was a 4965 K blackbody because no run's namelist carried the
    # spectrum. The identity below is the only statement in the scheme with a
    # right answer rather than a plausible one, so it is what the reproduction
    # is checked against.
    try:
        import stellar
        identity = stellar.solar_partition_identity()
        rep.add(OK, "solarini reproduction",
                f"5772 K gives {identity:.6f}; radmod.f90:207 says "
                f"{stellar.SOLAR_PARTITION}")
    except Exception as exc:
        rep.add(FAIL, "solarini reproduction", str(exc))
    else:
        band1 = stellar.band_fractions()[0]
        # -- the file against the blend it was resampled from ---------------
        #
        # The comparison below this one asks whether every stored copy of the
        # band-1 share agrees with every other, which is a question that can
        # only ever be answered "they differ" or "they do not". This one has a
        # right answer: `<name>_provenance.json` carries the same two
        # quantities integrated from the BT-Settl blend at ITS OWN resolution,
        # written by build_stellar_spectrum.py at build time, and the shipped
        # 2048-point file is supposed to be a resampling of exactly that. A
        # resampler that does not conserve flux makes it not, which is what
        # notes/audits/stellar-spectrum-oracle.md found. The source itself is
        # 26 MB of BT-Settl outside the repository; only its two integrals are
        # needed here, so the check costs nothing and needs nothing.
        try:
            low, high = stellar.spectrum_paths()
            prov = low.with_name(low.stem + "_provenance.json")
            if not prov.is_file():
                rep.add(WARN, "spectrum against its source",
                        f"{rel(high)} has no provenance record beside it, so it "
                        "cannot be checked against its blend. Rebuilding with "
                        "build_stellar_spectrum.py writes the record but "
                        "REWRITES the spectrum in place, and existing runs "
                        "refuse to resume across that (CONS-3); rebuild only "
                        "alongside a re-baseline, not to clear this warning")
            else:
                record = json.loads(prov.read_text(encoding="utf-8"))
                products = record.get("products") or {}
                describes = [
                    f"{f.name} has changed since {rel(prov)} was written"
                    for f in (low, high)
                    if (products.get(f.name) or {}).get("sha256") != sha256_of(f)]
                source = record.get("source_resolution")
                if describes:
                    rep.add(FAIL, "spectrum against its source",
                            "; ".join(describes) + " -- the record describes a "
                            "different file, so its integrals prove nothing")
                elif not source:
                    rep.add(FAIL, "spectrum against its source",
                            f"{rel(prov)} carries no source_resolution block, so "
                            "the file has never been checked against the blend "
                            "it came from. Rebuilding REWRITES the spectrum in "
                            "place and existing runs refuse to resume across "
                            "that (CONS-3); rebuild only alongside a "
                            "re-baseline")
                else:
                    d_band1 = band1 - float(source["band1_fraction"])
                    ratio = stellar.cross_section_ratio()
                    d_ratio = ratio / float(source["cross_section_ratio_um4"]) - 1.0
                    bad = []
                    if abs(d_band1) > stellar.SOURCE_BAND1_TOLERANCE:
                        bad.append(f"band 1 is {d_band1:+.2e} from source "
                                   f"resolution, over "
                                   f"{stellar.SOURCE_BAND1_TOLERANCE:.0e}")
                    if abs(d_ratio) > stellar.SOURCE_CROSS_SECTION_TOLERANCE:
                        bad.append(f"zcross/z1 is {d_ratio:+.2e} relative, over "
                                   f"{stellar.SOURCE_CROSS_SECTION_TOLERANCE:.0e}")
                    rep.add(FAIL if bad else OK, "spectrum against its source",
                            "; ".join(bad) + "; the resampler is the first thing "
                            "to look at" if bad else
                            f"{high.name} reproduces its BT-Settl blend: band 1 "
                            f"{d_band1:+.2e}, zcross/z1 {d_ratio:+.2e} relative")
        except Exception as exc:
            rep.add(WARN, "spectrum against its source", f"not checked: {exc}")

        stale = []
        for path, key in (
                (ROOT / "analysis" / "dust_optics.json",
                 "stellar_flux_fraction_band1"),
                (ROOT / "exoplasim" / "data" / "dust"
                 / "vesper_dust_aerosol.provenance.json",
                 "stellar_flux_fraction_band1")):
            if not path.is_file():
                continue
            stored = json.loads(path.read_text(encoding="utf-8")).get(key)
            if stored is None or abs(float(stored) - band1) > 5.0e-4:
                stale.append(f"{rel(path)} carries {stored}")
        rep.add(FAIL if stale else OK, "band-1 share is one value",
                "; ".join(stale) + f"; canonical is {band1:.4f}" if stale
                else f"{band1:.6f} from {stellar.spectrum_paths()[1].name}, "
                     "and every stored copy agrees")

    if want is None:
        # Everything below compares artifacts against the build. With no build
        # there is nothing to compare against, and each of those checks would
        # report a second, derived failure for the same one cause. The stellar
        # checks are ABOVE this line because none of them touches the build,
        # and a fresh clone with no terrain payload is exactly the state in
        # which a wrong spectrum file would otherwise go unnoticed.
        rep.add(WARN, "artifact checks",
                "skipped: they compare against the active build, which is absent")
        rep.show()
        return 1 if rep.failed else 0

    # -- terrain hash across every artifact that records one -----------------
    data = builds.component_data("hydrography", config)
    checked = 0
    for path in sorted(data.glob("*.nc")):
        got = nc_terrain(path)
        checked += 1
        if got is None:
            rep.add(WARN, f"terrain {path.name}", "carries no terrain_hash")
        elif got != want:
            rep.add(FAIL, f"terrain {path.name}",
                    f"is {got[:16]}, build is {want[:16]}")
        else:
            rep.add(OK, f"terrain {path.name}", got[:16])
    if not checked:
        rep.add(FAIL, "hydrography products", f"none in {data.relative_to(ROOT)}")

    resolution = str(config["model"]["resolution"]).lower()
    for path in sorted((ROOT / "exoplasim" / "inputs" / resolution).glob("*report*.json")):
        got = json_terrain(path)
        if got is None:
            rep.add(WARN, f"terrain {path.name}", "records no terrain_hash")
        elif got != want:
            rep.add(FAIL, f"terrain {path.name}",
                    f"is {got[:16]}, build is {want[:16]}")
        else:
            rep.add(OK, f"terrain {path.name}", got[:16])
    for path in (ROOT / "pedology" / "analysis" / "soil_report.json",):
        got = json_terrain(path)
        if got is not None:
            rep.add(OK if got == want else FAIL, f"terrain {path.name}",
                    got[:16] if got == want else f"is {got[:16]}, build is {want[:16]}")

    # -- gravity, which the terrain hash cannot see --------------------------
    try:
        from orogen import Export as _Export
        ex = _Export(builds.mesh_export(config))
        declared = float(config["planet"]["gravity_m_s2"])
        if abs(ex.gravity_m_s2 - declared) > 1e-6:
            rep.add(FAIL, "gravity build vs config",
                    f"build was generated at {ex.gravity_m_s2} m/s2, config "
                    f"declares {declared}; every vertical km is off by "
                    f"{declared / ex.gravity_m_s2:.4f}x")
        else:
            rep.add(OK, "gravity build vs config", f"{declared} m/s2")
    except Exception as exc:
        rep.add(WARN, "gravity build vs config", f"not checked: {exc}")

    # -- stale flat data -----------------------------------------------------
    flat = ROOT / "hydrography" / "data"
    if data != flat:
        stale = [p for p in flat.glob("*.nc") if nc_terrain(p) not in (None, want)]
        if stale:
            rep.add(WARN, "stale flat hydrography/data",
                    f"{len(stale)} files from another build; scripts defaulting "
                    "there will read the wrong terrain")

    # -- coupling convention -------------------------------------------------
    couplings = sorted(data.glob("coupling_*.nc"))
    if not couplings:
        rep.add(FAIL, "coupling matrices", "none found")
    # The invariant needs a climatology, and during a commissioning there is
    # not one yet: hydrography runs before the bootstrap, so the couplings
    # exist while the land mask to check them against does not. Fetched ONCE
    # and outside the loop, because `climatology_path` raises SystemExit with
    # none named, which took the whole report down between the two steps --
    # exactly where rule 8 requires this script to run. Unavailable is
    # reported per coupling below and is NOT a pass.
    lsm = None
    if couplings:
        try:
            lsm = land_sea_mask()
        except SystemExit as exc:
            lsm_unavailable = str(exc)
        except Exception as exc:                                  # noqa: BLE001
            lsm_unavailable = f"{type(exc).__name__}: {exc}"

    for path in couplings:
        # An invariant, not an existence assertion. `cell_lon` being PRESENT
        # says nothing about whether the mapping built from it is right, and
        # the mapping was wrong twice while this check passed both times.
        # Endorheic catchments are inland, so almost none of their area may
        # land on a cell the model calls ocean: 1.01% index for index against
        # 49.77% when longitude labels are matched. See docs/src/practice/failure-modes.md
        # class 17 and notes/audits/grid-convention-and-runoff.md.
        if lsm is None:
            rep.add(WARN, f"convention {path.name}",
                    f"not checkable yet: {lsm_unavailable}")
            continue
        with Dataset(path) as ds:
            same_grid = int(ds.n_lon) == lsm.shape[1]
        if not same_grid:
            # A coupling for another resolution cannot be checked against this
            # climatology. Reported rather than skipped silently, and NOT a
            # pass: nothing has been verified about it.
            rep.add(WARN, f"convention {path.name}",
                    "built for another resolution; not checkable against the "
                    "active climatology")
            continue
        try:
            frac = coupling_ocean_fraction(path, lsm)
            ok = frac <= 0.10
            detail = (f"{frac:.2%} of catchment area on model-ocean cells"
                      + ("" if ok else " -- NOT index-aligned"))
        except Exception as exc:                                  # noqa: BLE001
            ok, detail = False, f"not checked: {exc}"
        rep.add(OK if ok else FAIL, f"convention {path.name}", detail)

    # -- surface inputs ------------------------------------------------------
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        import run_exoplasim as rx
        codes = sorted(rx.intended_surface_codes(config))
        missing = [c for c in codes if not rx.surface_sra(config, c).is_file()]
        rep.add(FAIL if missing else OK, "surface inputs",
                f"missing {missing}" if missing else
                f"{len(codes)} present: {codes}")
        if not missing:
            # Currency by CONTENT, via the reports written alongside the .sra
            # files, not by mtime. An mtime comparison against the build
            # manifest fires whenever a file is restored or touched without
            # changing: restoring the export from a backup gave the manifest a
            # new mtime and every surface input read as stale, which is the same
            # false positive the cycle-binary check had. The terrain hash the
            # reports already carry is the thing that actually matters, and it
            # is checked above; this adds the check that a report EXISTS for the
            # active resolution, since a missing report is how an .sra survives
            # a build change unnoticed.
            res = str(config["model"]["resolution"]).lower()
            reports = sorted((ROOT / "exoplasim" / "inputs" / res).glob("*report*.json"))
            stale = []
            for r in reports:
                try:
                    got = json.loads(r.read_text(encoding="utf-8")).get("terrain_hash")
                except (OSError, json.JSONDecodeError):
                    got = None
                if got != want:
                    stale.append(r.name)
            rep.add(FAIL if (stale or not reports) else OK, "surface inputs current",
                    f"reports describe another build: {stale}" if stale else
                    (f"no surface reports under exoplasim/inputs/{res}"
                     if not reports else
                     f"{len(reports)} reports all describe {want[:16]}"))
    except Exception as exc:
        rep.add(WARN, "surface inputs", f"not checked: {exc}")

    # -- carve list ----------------------------------------------------------
    carve = data / "carve_list.json"
    if carve.is_file():
        d = json.loads(carve.read_text(encoding="utf-8"))
        entries = d["basins"]
        ids = [e["id"] for e in entries]
        # Against the basins this build still has plus the ones an earlier pass
        # carved, NOT against the build manifest's preserved list: that is the
        # post-carve survivors, and a list sent to Orogen has to cover the base
        # catalogue, which is larger.
        carried_n = sum(1 for e in entries if e.get("carried_from_previous_pass"))
        try:
            with Dataset(data / "basins.nc") as ds:
                here = ds.dimensions["basin"].size
        except Exception:
            here = None
        dup = len(ids) - len(set(ids))
        rep.add(FAIL if dup else OK, "carve list ids",
                f"{dup} duplicates" if dup else f"{len(ids)} unique")
        bad = [e["id"] for e in entries
               if e.get("carried_from_previous_pass") and float(e["retain"]) != 0.0]
        rep.add(FAIL if bad else OK, "carve list monotone",
                f"{len(bad)} carried entries are not retain 0" if bad else
                "carried entries all retain 0")
        if here is not None:
            expect = here + carried_n
            rep.add(OK if len(ids) == expect else FAIL, "carve list accounts for every basin",
                    f"{len(ids)} entries = {here} decidable + {carried_n} carried"
                    if len(ids) == expect else
                    f"{len(ids)} entries, expected {here} + {carried_n} = {expect}")

    # -- config --------------------------------------------------------------
    try:
        derived = rx.derive(config, float(config["orbit"]["baseline_flux_earth"]))
        rep.add(OK, "config gravity vs mass", "consistent")
        rep.add(OK, "derived orbit",
                f"a = {derived['semimajor_axis_au']:.6f} AU, "
                f"year = {derived['orbital_year_earth_days']:.3f} d")
    except Exception as exc:
        rep.add(FAIL, "config", str(exc))

    # -- the active timestep against the ladder registry ---------------------
    #
    # `model.timestep_minutes` is the step ONE run is configured at, and it is
    # judged against two quantities that `lib/rungs.py` declares: the coarsest
    # step the rung is MEASURED to start clean at, and the step the escalation
    # route runs that rung at. Both, not one -- a step can be comfortably below
    # the ceiling and still not be the escalation the project declared, and a
    # route step on a rung nothing has probed is a step nothing has qualified.
    #
    # This checks the DECLARED configuration, so both verdicts are failures
    # here. A diagnostic arm deliberately runs off the route and reports rather
    # than refusing; `lib/rungs.py:timestep_problems` says why the two callers
    # differ. WORLD-F997, WORLD-J37.
    try:
        model = config["model"]
        active, rung = rungs.configured_timestep(config)
        problems = rungs.timestep_problems(rung, active)
        for problem in problems:
            rep.add(FAIL, "timestep vs the ladder registry", problem)
        if not problems:
            rep.add(OK, "timestep vs the ladder registry",
                    f"{rung} at {active} min is on the escalation route and at "
                    f"or below the measured ceiling "
                    f"{rungs.stability_ceiling(rung)}, kappa "
                    f"{model.get('filter_kappa')}")
    except Exception as exc:
        rep.add(FAIL, "timestep vs the ladder registry", f"not checked: {exc}")

    # -- the registry's own tables against what they restate -----------------
    #
    # The ceiling table is read out of the probe grid and the route table is
    # read out of the chapter that decides the route, so both can drift from
    # their source silently. They are checked HERE as well as in smoke_test
    # because this is what runs before an expensive run.
    for label, problems in (
            ("stability ceilings vs the probe grid",
             rungs.check_stability_ceilings(ROOT)),
            ("escalation route vs sequencing.md section D",
             rungs.check_timestep_restatements(ROOT))):
        for problem in problems:
            rep.add(FAIL, label, problem)
        if not problems:
            rep.add(OK, label, "agree")

    # -- the energy fixer is declared, coherent, and matches the runs that ran --
    #
    # `model.energy_fixer` is a CORRECTION for the conversion defect (world-0ov),
    # not physics, and whether a run carries it changes what its surface fluxes
    # mean. Two things can go wrong quietly. It can be switched on without the
    # diagnostics it is driven by, in which case the model runs with the switch
    # set and the fixer does nothing. And a run on disk can claim it on its
    # manifest while its namelist says otherwise, which is the claim the artifact
    # is supposed to settle.
    try:
        model = config["model"]
        fix = model.get("energy_fixer")
        if fix is None:
            rep.add(FAIL, "energy fixer declared",
                    "model.energy_fixer is absent; whether a run carries the "
                    "correction changes what its fluxes mean, so it is declared")
        elif fix and not model.get("energy_diagnostics"):
            rep.add(FAIL, "energy fixer declared",
                    "model.energy_fixer is on and model.energy_diagnostics is "
                    "off; the fixer is driven by denergy26 and denergy27 and "
                    "would be a no-op with the switch set")
        else:
            rep.add(OK, "energy fixer declared",
                    f"{'on' if fix else 'off'}, a correction for world-0ov")

        runs = ROOT / "exoplasim" / "runs"
        seen, bad = 0, []
        for manifest in sorted(runs.glob("run_*/run_manifest.json")):
            claim = json.loads(manifest.read_text()).get("energy_fixer")
            if claim is None:
                continue          # written before the fixer existed
            nl = manifest.parent / "plasim_namelist"
            if not nl.is_file():
                continue
            on = any(line.strip().upper().replace(" ", "").startswith("NENERGYFIX=1")
                     for line in nl.read_text().splitlines())
            seen += 1
            if bool(claim) != on:
                bad.append(f"{manifest.parent.name} claims {bool(claim)}, "
                           f"namelist says {on}")
        if bad:
            rep.add(FAIL, "runs vs the energy fixer they claim",
                    "; ".join(bad[:4]))
        else:
            rep.add(OK, "runs vs the energy fixer they claim",
                    f"{seen} runs carry the stamp and all agree with their namelist")

        # -- what each run INTEGRATED with, against what was declared -------
        #
        # This lived here as a second gate and is now one: `audit_runs` above
        # runs "runs vs the hyperdiffusion they declare" over every run against
        # its OWN manifest, through `run_exoplasim.namelist_elements`. That
        # matters because the keys are written in Fortran's `n*value`
        # replication -- runs on disk carry `TDISSD = 10*0.0737` -- and the
        # version that stood here parsed with a bare `float()`, so a correctly
        # replicated array raised, was skipped, and was reported MISSING. It
        # could not see world-720's real defect either, a scalar where NLEV
        # elements are wanted, because it never compared lengths.
        #
        # The finding that motivated it stands and is world-1nz: a continuation
        # re-copies the shipped namelists, so a setting written only at prepare
        # time reverts partway through, and the hyperdiffusion fallback is a
        # DIFFERENT OPERATOR, ndel 2 against 4 at T21.

        # -- and what the fixer actually had to put back --
        #
        # The fixer HIDES the defect it compensates: with it on, denergy26 minus
        # denergy27 reports a residual near zero and the size of the loss is the
        # correction instead. A correction nobody checks turns a known 0.85 W/m2
        # into an unknown one that can grow, so the recorded magnitude is
        # checked rather than merely stored. The bound is an order and a half
        # above the defect, which catches a change in kind without firing on the
        # controller's ordinary wander.
        loud, quiet = [], []
        for manifest in sorted(runs.glob("run_*/run_manifest.json")):
            data = json.loads(manifest.read_text())
            if not data.get("energy_fixer"):
                continue
            applied = data.get("energy_fixer_applied")
            if applied is None:
                quiet.append(manifest.parent.name)
                continue
            mean = float(applied["applied_w_m2_mean"])
            if abs(mean) > 10.0:
                loud.append(f"{manifest.parent.name} applied {mean:+.2f} W/m2")
        if loud:
            rep.add(FAIL, "what the energy fixer had to put back",
                    "; ".join(loud[:4]) + " -- the defect on world-0ov is about "
                    "0.6 W/m2, so this is a change in kind, not a wander")
        elif quiet:
            rep.add(WARN, "what the energy fixer had to put back",
                    f"{len(quiet)} fixer runs predate the manifest field and "
                    f"report nothing; the runner now refuses to finish without it")
        else:
            rep.add(OK, "what the energy fixer had to put back",
                    "every fixer run records its correction and none is a "
                    "change in kind")
    except Exception as exc:
        rep.add(WARN, "energy fixer declared", f"not checked: {exc}")

    # -- what the runs already on disk actually integrated -------------------
    #
    # Everything above this line reads the configuration and the artifacts a
    # generator wrote. `audit_runs` reads RUN DIRECTORIES, which is where the
    # model's own record of what it integrated lives. Its argument, and the
    # reasoning for each case and its severity, are on the function.
    try:
        # Appended, not inserted: `lib/` must keep priority, and nothing under
        # `exoplasim/scripts/` shares a module name with it today.
        scripts = str(ROOT / "exoplasim" / "scripts")
        if scripts not in sys.path:
            sys.path.append(scripts)
        import run_exoplasim as runner          # noqa: E402
    except Exception as exc:
        rep.add(WARN, "runs vs the namelists they staged",
                f"not checked: run_exoplasim.py did not import: {exc}")
    else:
        audit_runs(rep, ROOT / "exoplasim" / "runs", runner)

    # -- the derived diffusion table matches the rule it claims to come from --
    #
    # `model.hyperdiffusion.timescales_days` is DERIVED: vorticity damps on the
    # advective time at the smallest resolved scale, pi*radius/(NTRU*eddy_wind),
    # and the other three keep their ratios to it. A table that has drifted from
    # its own rule is worse than no table, because it reads as derived.
    try:
        model = config["model"]
        hd = model.get("hyperdiffusion") or {}
        if not hd:
            rep.add(WARN, "hyperdiffusion vs its rule", "no block declared")
        else:
            radius = float(config["planet"]["radius_earth"]) * 6371e3
            wind = float(hd["eddy_wind_m_s"])
            ratios = hd["ratios_to_vorticity"]
            bad = []
            for rung, tau in hd["timescales_days"].items():
                ntru = int(str(rung).lstrip("Tt"))
                base = math.pi * radius / (ntru * wind) / 86400.0
                want = {"vorticity": base,
                        "divergence": base * float(ratios["divergence"]),
                        "temperature": base * float(ratios["temperature"]),
                        "humidity": base * float(ratios["humidity"])}
                for name, value in want.items():
                    got = float(tau[name])
                    # 1% of the value: the table is quoted to four decimals and
                    # the ratios to three, so exact equality is not available.
                    if abs(got - value) > 0.01 * value:
                        bad.append(f"{rung}.{name}: table {got:g}, rule {value:.4f}")
            rep.add(FAIL if bad else OK, "hyperdiffusion vs its rule",
                    "; ".join(bad) if bad else
                    f"{len(hd['timescales_days'])} rungs derive from "
                    f"pi*a/(NTRU*{wind} m/s)")
    except Exception as exc:
        rep.add(WARN, "hyperdiffusion vs its rule", f"not checked: {exc}")

    # -- the filter is confined to the scales it is meant to damp -------------
    #
    # WHAT THE FILTER DOES, because the gate follows from where it is applied.
    # `legmod.f90` folds `f(n) = exp(-kappa (n/NTRU)^gamma)` into every per-mode
    # transform weight in both directions, and what passes through those weights
    # is the nonlinear TENDENCY: no filtered value is ever stored back into the
    # prognostic spectral state. So the filter attenuates the SOURCE that fills a
    # wavenumber by a fixed factor rather than damping the state at a rate, a
    # source scaled by `f^2` against an unchanged sink equilibrates at
    # `spec/law = f(x)^2 = exp(-2 kappa x^gamma)`, and the filter takes half of
    # what the flow would otherwise put at
    #
    #     x_f = (ln 2 / (2 kappa)) ** (1/gamma),   x = n/NTRU
    #
    # which carries no timestep and no rung. That REPLACES the rate-law crossover
    # `(dt / (2 kappa tau))**(1/(gamma-1))` this gate used to solve. Two timestep
    # pairs differing in MPSTEP alone, same binary by sha and same cold start,
    # move the measured spectrum by 0.1 and 0.2 sigma where that law predicts
    # several times the instrument's own three-sigma bar.
    # `exoplasim/notes/filter-spectral-price.md`.
    #
    # THE FLOOR IS UNCHANGED AT 0.60. What moved is the quantity, not the bar.
    # Re-choosing the bar to suit the corrected quantity would be choosing a
    # criterion after seeing the result it judges.
    #
    # THIS GATE DOES NOT PREDICT THE MEASURED BITE POINT, and the old one's
    # claim to was its second defect: it reported 0.625 at T21 where the run
    # measured 0.561 and 0.655 at T42 where the run measured 0.714, wrong in
    # both directions at once. A run's departure from its own inertial range
    # carries the derived hyperdiffusion and the spectrum's own steepening as
    # well, and the hyperdiffusion is the larger of the two levers -- reaching
    # every model level moves the T21 spectrum further than the whole of
    # `filter_kappa` is worth. That criterion is a property of a RUN and has an
    # instrument: `exoplasim/scripts/spectral_tail.py` fails a run whose bite
    # point sits below 0.6 of the truncation. A config gate has no spectrum to
    # read, and asserting one from arithmetic is what produced the disagreement.
    #
    # WHAT IT CAN SEE. `x_f` goes as `kappa**(-1/gamma)`, so at gamma 16 it is
    # nearly blind to kappa and is in practice a guard on GAMMA -- which is the
    # degree of freedom that moves it, and which was 8 until 2026-08-23. The
    # kappa the floor holds up to at the declared gamma is reported so the
    # margin is a number, and the filter's own depletion at 0.8 of the
    # truncation is reported beside it because that one IS linear in kappa.
    try:
        model = config["model"]
        kappa = float(model["filter_kappa"])
        gamma = int(model["filter_power"])
        floor = float(model.get("filter_confinement_floor", 0.60))
        reach = (math.log(2.0) / (2.0 * kappa)) ** (1.0 / gamma)
        kappa_ceiling = math.log(2.0) / (2.0 * floor ** gamma)
        depth_dex = 2.0 * kappa * 0.8 ** gamma / math.log(10.0)
        if not model.get("physics_filter"):
            rep.add(OK, "filter confinement", "no filter configured")
        elif reach < floor:
            rep.add(FAIL, "filter confinement",
                    f"kappa {kappa:g}, gamma {gamma}: the filter takes half of "
                    f"the flow's own amplitude by {reach:.3f} of the truncation, "
                    f"below the {floor:.2f} floor -- everything above that is "
                    f"resolution being damped away. At gamma {gamma} the floor "
                    f"holds to kappa {kappa_ceiling:.4g}")
        else:
            rep.add(OK, "filter confinement",
                    f"kappa {kappa:g}, gamma {gamma}: half-attenuation at "
                    f"{reach:.3f} of the truncation (floor {floor:.2f}, which at "
                    f"this gamma holds to kappa {kappa_ceiling:.4g}); the filter "
                    f"alone takes {depth_dex:.3f} dex at 0.8 of the truncation")
    except Exception as exc:
        rep.add(WARN, "filter confinement", f"not checked: {exc}")

    # -- config blocks declare their determination status --------------------
    #
    # A value that has never been decided must not be indistinguishable from one
    # that has. stellar_cycle carried an amplitude and a centre from before there
    # was a baseline flux to centre on, and both were later read back as though
    # they had been chosen.
    try:
        import re as _re
        text = (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
        lines = text.splitlines()
        unmarked = []
        for i, line in enumerate(lines):
            if not _re.match(r"^[a-z_]+:\s*$", line):
                continue
            # Scan the WHOLE contiguous comment block above the key, not a fixed
            # number of lines. A fixed window silently fails a block whose
            # rationale grew past it, which is backwards: the blocks with the
            # most to explain are exactly the ones that most need a status, and
            # an 8-line window flagged `orbit` as unmarked when its marker was
            # the first line of an 11-line header.
            j = i - 1
            while j >= 0 and (lines[j].startswith("#") or not lines[j].strip()):
                if not lines[j].strip() and j < i - 1:
                    break          # a blank line ends the block's own header
                j -= 1
            window = "\n".join(lines[j + 1:i])
            if not any(k in window for k in
                       ("DETERMINED", "DECLARED", "PROVISIONAL", "DERIVED",
                        "MIXED", "TRANSITIVE")):
                unmarked.append(line.split(":")[0])
        rep.add(FAIL if unmarked else OK, "config blocks declare their status",
                f"unmarked: {', '.join(unmarked)}" if unmarked else
                "every top-level block is marked")
    except Exception as exc:
        rep.add(WARN, "config blocks declare their status", f"not checked: {exc}")

    # -- generated biosphere inputs vs the config they were derived from -----
    #
    # biosphere/generated/ holds a C++ header, a driver and a PFT file, each
    # DERIVED from config/planet.yaml and each recording the config it was built
    # against. That design is right and it is why the year length is not written
    # down anywhere by hand. But nothing verified the derivation had actually
    # been re-run, so all three drifted -- and to three DIFFERENT config states,
    # which is worse than one stale artifact because they disagree with each
    # other as well as with the config.
    #
    # The year length sizes arrays in a compiled header, so it cannot be checked
    # at runtime. This is the check.
    #
    # It compares PARSED VALUES, not a hash of the file. A hash of
    # config/planet.yaml cannot tell an edited comment from an edited parameter,
    # so this fired on a two-line comment change and named re-running three
    # generators as the remedy. That is worse than not firing: the remedy sounds
    # right, it makes the message go away, and what it teaches is to re-run a
    # generator to silence a hash rather than because anything moved. The resume
    # guard had the same defect and the same fix, and `config_drift` is now the
    # one implementation of it -- see lib/provenance.py.
    try:
        gen = ROOT / "biosphere" / "generated"
        cur = sha256_of(ROOT / "config" / "planet.yaml")
        stale = []
        for prov in sorted(gen.glob("*_provenance.json")):
            rec = json.loads(prov.read_text(encoding="utf-8"))
            yl = rec.get("year_length_days")
            # Each artifact names its own generator, so the remediation is read
            # off the artifact rather than guessed. It used to name one script
            # for all three, which was wrong for two of them.
            who = rec.get("generator") or "(generator not recorded)"
            what = f"{prov.name}" + (f" (year_length_days={yl})" if yl else "")
            was = rec.get("source_config")
            if was is None:
                # A record from before the generators kept the parsed config.
                # Its bytes are all there is, so the file hash is all that can be
                # compared -- and it is reported as what it is, a comparison that
                # cannot separate a comment from a parameter, rather than as a
                # verdict on the artifact.
                if rec.get("config_sha256") != cur:
                    stale.append(f"{what} records no source_config and the "
                                 f"config file has changed, which may be only a "
                                 f"comment -> rerun {who} to settle it")
                continue
            drift = config_drift(was, config, BIOSPHERE_INERT_CONFIG_KEYS)
            if drift:
                stale.append(f"{what}: " + "; ".join(drift) + f" -> rerun {who}")
        if not list(gen.glob("*_provenance.json")):
            rep.add(WARN, "generated biosphere inputs", "none present")
        else:
            rep.add(FAIL if stale else OK, "generated biosphere inputs current",
                    "derived from an older config: " + "; ".join(stale)
                    if stale else "all derived from the current config")
    except Exception as exc:
        rep.add(WARN, "generated biosphere inputs", f"not checked: {exc}")

    # -- staged surface fields vs the config they were built from ------------
    #
    # The same question the biosphere check above asks, of the fields the model
    # reads at run start. It went unasked until 2026-08-19, when SPEC-5's
    # star-reweighted canopy sat unstaged through two commits: `surface inputs
    # current` above tests the TERRAIN hash, and the terrain had not moved.
    #
    # Two failures, not one, because they are different. An artifact whose
    # config has drifted is WORTHLESS and its generator has to re-run. An
    # artifact with no `source_config` at all is UNOBSERVABLE, which is not the
    # same as current, and the remedy is the same re-run either way.
    try:
        import pipeline as _pipeline
        graph = _pipeline.load()
        build = _pipeline.active_build()
        drifted, unstamped = [], []
        for step in graph["steps"]:
            inert = INERT_CONFIG_KEYS.get(step["id"])
            if inert is None:
                continue
            for w in step.get("writes", []):
                path = _pipeline.resolve(w, build)
                if not path.is_file() or path.suffix != ".json":
                    continue
                drift = artifact_drift(path, config, inert)
                if drift is None:
                    unstamped.append(f"{path.name} -> rerun {step.get('script')}")
                elif drift:
                    drifted.append(f"{path.name}: {'; '.join(drift)} -> rerun "
                                   f"{step.get('script')}")
        rep.add(FAIL if drifted else OK, "generated inputs vs their config",
                "built from an older config: " + "; ".join(drifted) if drifted
                else "every stamped artifact matches the current config")
        if unstamped:
            rep.add(WARN, "generated inputs without a config stamp",
                    "; ".join(unstamped))
        # The inert sets are claims about what a generator reads, and a claim
        # that has gone stale is worse than none. Re-run the trace: a key listed
        # inert whose NAME appears in the generator's source is a claim the
        # source contradicts.
        contradicted = []
        for step in graph["steps"]:
            inert = INERT_CONFIG_KEYS.get(step["id"])
            script = step.get("script")
            if inert is None or not script:
                continue
            src_path = ROOT / script
            if not src_path.is_file():
                continue
            src = src_path.read_text(encoding="utf-8")
            for key in sorted(inert):
                # Any block, not just `model.`. The sets were model-only when
                # this was written, so skipping other prefixes cost nothing;
                # once `ocean.horizontal_diffusion` joined one, skipping meant
                # the honesty guard silently stopped covering the entries most
                # likely to be wrong. Trace every key or the claim is unchecked.
                if "." not in key:
                    continue
                if key.split(".", 1)[1] in src:
                    contradicted.append(f"{step['id']} lists {key} inert and "
                                        f"{script} names it")
        rep.add(FAIL if contradicted else OK, "inert config sets vs their generators",
                "; ".join(contradicted) if contradicted
                else "no generator names a key its own inert set calls unread")
    except Exception as exc:
        rep.add(WARN, "generated inputs vs their config", f"not checked: {exc}")

    # The same honesty question asked of the removals, which are the other thing
    # that excuses a refusal. An inert entry goes stale by naming nothing; a
    # removal goes stale by naming something that is BACK, or that a script has
    # since learned to read. Either way it is a decision made once that keeps
    # applying to a tree that has moved, so it is checked here rather than only
    # under --self-test. world-51wj.
    stale = removal_problems(config, ROOT)
    rep.add(FAIL if stale else OK, "declared config removals vs the tree",
            "; ".join(stale) if stale
            else f"{len(REMOVED_CONFIG_KEYS)} declared, each gone from the "
                 f"config and named by no script outside where it moved")

    # -- staged fields vs the DERIVED FILES they were built from -------------
    #
    # The check above compares an artifact's CONFIG stamp against the current
    # config, which is the right question for a configuration key and cannot see
    # this one. A generator also reads derived files that no config key names --
    # the two band-shape derivations under `analysis/` that set staged surface
    # codes 175 and 176, a dust field, an optics table, a land-column contract.
    # Re-running one of those derivations moves the file while every config key
    # stays put, so the staged .sra goes on describing the previous shapes and
    # the config check reports everything current. world-nvs2.
    #
    # It compares HASHES and not names, because the name is what is stable
    # across the change: `analysis/rock_albedo_bands.json` is the same path
    # before and after its own generator re-runs.
    #
    # It walks the pipeline graph's own `writes` rather than `INERT_CONFIG_KEYS`
    # ON PURPOSE. An inert set is a claim about which config keys a generator
    # reads and has to be traced by hand, so steps whose trace has not been done
    # are deliberately absent from it -- `surface_dust` among them. An input
    # hash needs no such claim, so this covers every step that stamps one.
    try:
        import pipeline as _pipeline
        graph = _pipeline.load()
        build = _pipeline.active_build()
        moved, stamped = [], 0
        for step in graph["steps"]:
            for w in step.get("writes", []):
                path = _pipeline.resolve(w, build)
                lines = artifact_input_drift(path)
                if lines is None:
                    continue
                stamped += 1
                if lines:
                    moved.append(f"{path.name}: {'; '.join(lines)} -> rerun "
                                 f"{step.get('script')}")
        if stamped == 0:
            rep.add(WARN, "generated inputs vs the files they read",
                    "no artifact on disk records the files it read")
        else:
            rep.add(FAIL if moved else OK,
                    "generated inputs vs the files they read",
                    "; ".join(moved) if moved
                    else f"{stamped} stamped artifacts, every recorded input "
                         f"unchanged since it was read")
    except Exception as exc:
        rep.add(WARN, "generated inputs vs the files they read",
                f"not checked: {exc}")

    # -- runs vs the spectrum file they were integrated against --------------
    #
    # The run manifest records the spectrum by CONTENT as well as by name,
    # because `build_stellar_spectrum.py` rewrites `<name>.dat` in place from
    # `star.spectral_type` and `star.effective_temperature_k`. The name is
    # therefore stable across a change that moves every snow, ice and glacier
    # albedo, and the resume guard compares parsed config values, which cannot
    # see a file rewritten under an unchanged name. `continue_exoplasim.py`
    # refuses to resume across it; this says the same thing before an expensive
    # run rather than at the moment one is being extended.
    #
    # A finished run on another spectrum is a fact rather than a defect -- it is
    # a different climate and stays readable as one -- so a mismatch is reported
    # and not failed. What it means is that the run cannot be extended.
    try:
        import run_exoplasim as _rx
        current = _rx.stellar_spectrum_digest(config)
        runs = ROOT / "exoplasim" / "runs"
        same, other, unstamped = [], [], []
        for m in sorted(runs.glob("*/run_manifest.json")):
            try:
                rec = json.loads(m.read_text(encoding="utf-8")).get("stellar_spectrum_digest")
            except (OSError, json.JSONDecodeError):
                rec = None
            if rec is None:
                unstamped.append(m.parent.name)
            elif rec == current:
                same.append(m.parent.name)
            else:
                other.append(m.parent.name)
        if not (same or other or unstamped):
            rep.add(WARN, "runs vs stellar spectrum", "no run manifests present")
        elif other or unstamped:
            detail = []
            if other:
                detail.append(f"{len(other)} on another spectrum file and not "
                              f"resumable: {other}")
            if unstamped:
                detail.append(f"{len(unstamped)} record no spectrum digest and "
                              f"predate the check: {unstamped}")
            rep.add(WARN, "runs vs stellar spectrum", "; ".join(detail))
        else:
            name = (current or {}).get("file", "a blackbody")
            rep.add(OK, "runs vs stellar spectrum",
                    f"{len(same)} runs all on {name}")
    except Exception as exc:                                       # noqa: BLE001
        rep.add(WARN, "runs vs stellar spectrum", f"not checked: {exc}")

    # -- binaries vs the patches they should contain -------------------------
    #
    # ExoPlaSim builds one executable per (resolution, layers, ranks)
    # configuration, so
    # changing the source and rebuilding touches only the configuration in use
    # and leaves the rest silently stale. That is failure class 11, and it fired
    # three times in one day. Vendoring the model as a subtree does not fix it:
    # a stale binary is still a stale binary, so this check outlived the patch
    # stack it was originally written to police.
    # ASKED OF `rebuild_binaries.py:verify`, which is the operation rule 4 names,
    # rather than restated here. This block used to carry its own copy of the
    # comparison and the copy had already fallen behind the original in three
    # ways: it globbed `plasim/run` alone, so an executable in `plasim/bin` kept
    # a provenance nothing looked at; it could not see a MATRIX row with NO
    # binary, which is rule 4's other half and is invisible to anything that
    # globs the directory; and it never compared the toolchain, so a compiler or
    # a flag line that had moved under unchanged sources passed. `verify` builds
    # nothing and links nothing -- a manifest read and a sha per source. world-wkci.
    try:
        scripts = str(ROOT / "exoplasim" / "scripts")
        if scripts not in sys.path:
            sys.path.append(scripts)
        import rebuild_binaries as rb           # noqa: E402
        report = rb.verify()
    except Exception as exc:                                       # noqa: BLE001
        rep.add(WARN, "binaries", f"not checked: {exc}")
    else:
        if not report["manifest_present"]:
            rep.add(FAIL, "binary manifest",
                    "absent; run exoplasim/scripts/rebuild_binaries.py")
        elif not report["built"]:
            # A tree with no executables is UNBUILT, not inconsistent. That is
            # the resting state of a fresh worktree, where the build directories
            # are deliberately not linked in, and calling it a failure would
            # make this gate red everywhere before any work had gone wrong.
            rep.add(WARN, "binaries", "none built")
        else:
            bad = []
            if report["absent"]:
                bad.append("provenance unknown (not in the manifest, or in it "
                           f"under another sha): {', '.join(report['absent'])}")
            if report["stale"]:
                bad.append("built from source that has since moved: "
                           f"{', '.join(report['stale'])}")
            if report["missing"]:
                # Only once SOMETHING is built. A partial set is exactly rule 4's
                # shape -- a rebuild that refreshed the configuration in use and
                # left the rest -- and it is the case the old block could not see.
                bad.append("in the matrix and not built: "
                           f"{', '.join(report['missing'])}")
            if report["drift"]:
                bad.append("toolchain has moved since the build: "
                           + "; ".join(report["drift"]))
            rep.add(FAIL if bad else OK, "binaries carry current patches",
                    "; ".join(bad) if bad else
                    f"{len(report['built'])} executables match the manifest")

    # -- the cycle executable can parse the cycle the config asks for ---------
    #
    # The star-cycle patch is deliberately NOT resident, so this binary is
    # outside the matrix above and nothing else would notice it going stale.
    # That is the same gap failure class 11 came through: exactly one of six
    # executables was newer than the patched source, and the run died in
    # radini_ with "Cannot match namelist object name".
    #
    # Checking the namelist names are actually IN the binary is stronger than
    # checking a timestamp, and it is the specific thing that fails: adding a
    # component to config/planet.yaml without rebuilding leaves a binary whose
    # namelist has no slot for it, and Fortran rejects the whole group.
    try:
        m = config["model"]
        run_dir = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "run"
        # THE NAME IS THE CONFIGURATION AND NOTHING ELSE: resolution, layers and
        # the thread count, which is what `build_model.executable_name` composes
        # and what a run goes looking for. A check that names it any other way
        # asserts a binary that is not built and never looks at the one that is.
        # world-bdh, world-38b.
        exe = run_dir / (f"most_plasim_t{int(str(m['resolution']).lstrip('Tt'))}"
                         f"_l{int(m['layers'])}_p{int(m['ncpus'])}.x")
        components = (config.get("stellar_cycle") or {}).get("components") or {}
        slots = {"medium": "", "long": "2"}
        # There is no separate cycle executable any more. The star-cycle change
        # is resident in the subtree and gated by the namelist: nsolcycle
        # defaults to 0 and both amplitudes to 0.0, which reduces the guard to
        # `gsolinst = gsol0`, so the ordinary binaries are bit-exact identical
        # to an unpatched model until the cycle is configured on. What is worth
        # checking is therefore only that the binary the run will use actually
        # carries a namelist slot for every component the config declares.
        if not components:
            rep.add(WARN, "cycle capability", "no stellar_cycle.components")
        elif not exe.is_file():
            rep.add(FAIL, "cycle capability",
                    f"{exe.name} absent; run exoplasim/scripts/rebuild_binaries.py")
        else:
            blob = exe.read_bytes()
            missing = [f"gsolamp{slots[n]}" for n in components
                       if n in slots
                       and f"gsolamp{slots[n]}".encode() not in blob]
            unslotted = sorted(set(components) - set(slots))
            problems = []
            if unslotted:
                problems.append(f"no namelist slot for {unslotted}")
            if missing:
                problems.append(f"binary lacks {missing}")
            rep.add(FAIL if problems else OK,
                    "cycle capability is compiled into the run binary",
                    "; ".join(problems) if problems else
                    f"{len(components)} components, all present in {exe.name}")
    except Exception as exc:
        rep.add(WARN, "cycle executable", f"not checked: {exc}")

    # -- the flux-to-kelvin slope still matches the runs it was measured on, --
    #    and those runs are still on the world this tree describes
    try:
        import sensitivity
        sensitivity.test_identity()
        stale = sensitivity.verify()
        rep.add(FAIL if stale else OK, "flux-to-kelvin slope",
                "; ".join(stale) if stale else
                f"{sensitivity.SLOPE_K_PER_FLUX_RATIO} K per unit flux ratio, "
                f"reproduced from the runs lib/sensitivity.py names, both live "
                f"on {sensitivity.active_build()}")
    except Exception as exc:
        rep.add(WARN, "flux-to-kelvin slope", f"not checked: {exc}")

    # -- and the currency check can still tell a superseded bracket apart ----
    #
    # Separate row because it is a different question. The one above asks
    # whether the declaration is current; this asks whether the instrument that
    # answers that can return no. It failed silently for a whole build.
    try:
        import sensitivity
        sensitivity.test_currency_refuses_a_superseded_measurement()
        rep.add(OK, "the slope's currency check can fail",
                "the superseded bracket's two deleted runs are both refused")
    except AssertionError as exc:
        rep.add(FAIL, "the slope's currency check can fail", str(exc))
    except Exception as exc:
        rep.add(WARN, "the slope's currency check can fail",
                f"not checked: {exc}")

    check_derived_config_values(rep, config)
    check_derived_config_rules(rep, config)
    check_land_column_thermal_constants(rep)
    check_cold_start_currency(rep, config)
    check_eddy_wind(rep, config)
    check_water_path_currency(rep, config)
    rep.show()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
