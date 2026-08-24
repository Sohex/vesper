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

**What the runs on disk actually integrated.** Three cases over run directories,
in `audit_runs`. The staged namelists against the config and the manifest each
run carries, the staged namelists against the hyperdiffusion block the manifest
declares, and the per-code `surface_field_sha256` against the `.sra` files that
are there now. Until these existed, `run_exoplasim.py`'s `expected_namelist_keys`
and `verify_staged_namelists` were the only guard on any of it and both fire only
at prepare and continue time, so a run already on disk was invisible: the T21 run
the project's numbers rest on records a hyperdiffusion block its namelist does not
carry, and a person reading run directories is what found it.

**The flux-to-kelvin slope against the runs it was measured on.** `lib/sensitivity.py`
declares one sensitivity for the whole project. Its predecessor was called
FALLBACK_SLOPE, nothing ever fell back to it, and it survived two terrain changes
and a resolution change while three incompatible values accumulated around it.
This recomputes the declared slope from the run index every time.

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
from paths import rel                # noqa: E402
from gridding import coupling_ocean_fraction   # noqa: E402
from provenance import artifact_drift, BIOSPHERE_INERT_CONFIG_KEYS, INERT_CONFIG_KEYS, config_drift, unknown_inert_keys   # noqa: E402


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
    case("a run that agrees with itself passes all three",
         {check: status for check, (status, _) in green.items()},
         {NAMELISTS: OK, HYPERDIFF: OK, SURFACE: OK},
         "otherwise every case below would pass on a gate that only ever fails")
    case("and the fixture is a run the gates actually read",
         [green[check][1].split()[0] for check in (NAMELISTS, HYPERDIFF, SURFACE)],
         ["1", "1", "1"],
         "an empty runs directory would report three OKs as well")

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

    width = max(len(n) for n, _, _ in cases) + 2
    for name, ok, why in cases:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name:<{width}} {why}")
    print(f"\n{len(cases)} cases, {len(failures)} failed")
    for f in failures:
        print(f"  {f}")
    return 1 if failures else 0


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

    # -- the named baseline climatology belongs to the active build ----------
    #
    # A climatology is the most widely shared artifact here: pedology,
    # hydrography and the biosphere all read one. It has the same grid, the
    # same variables and the same units whichever world it describes, so a
    # superseded one produces a plausible number rather than an error. The
    # consumers now check at the point of reading; this is the same check
    # applied once, ahead of an expensive run.
    declared = config.get("baseline_climatology")
    if not declared:
        rep.add(WARN, "baseline climatology",
                "none named in config; anything needing a climate will raise")
    else:
        clim = ROOT / declared
        if not clim.is_file():
            rep.add(FAIL, "baseline climatology",
                    f"config names {declared}, which does not exist")
        else:
            sys.path.insert(0, str(ROOT / "lib"))
            from provenance import artifact_build
            got = artifact_build(clim)
            if got is None:
                rep.add(WARN, "baseline climatology",
                        f"{Path(declared).name} carries no build identity; it "
                        f"predates the stamping in build_climatology.py")
            elif got != build:
                rep.add(FAIL, "baseline climatology",
                        f"{Path(declared).name} is on {got}, config names {build}")
            else:
                rep.add(OK, "baseline climatology", f"on {build}")

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

    # -- the active timestep matches the per-rung table ----------------------
    #
    # The filter buys timestep, so the step a rung can take is a MEASURED
    # property of that rung and of `model.filter_kappa`, not a preference. It
    # lives in `model.resolution_timestep_minutes`; `model.timestep_minutes` is
    # the active value. Two places holding a timestep is two places for them to
    # disagree, and the disagreement would be silent -- the model runs happily
    # at a step nobody chose for the resolution it is running.
    try:
        model = config["model"]
        table = model.get("resolution_timestep_minutes") or {}
        rung = str(model.get("resolution", "")).upper()
        active = float(model["timestep_minutes"])
        if not table:
            rep.add(WARN, "timestep vs the per-rung table", "no table declared")
        elif rung not in table:
            rep.add(FAIL, "timestep vs the per-rung table",
                    f"resolution {rung} has no entry; the adopted step for a rung "
                    f"is measured, so a rung with no entry has no adopted step")
        elif abs(float(table[rung]) - active) > 1e-9:
            rep.add(FAIL, "timestep vs the per-rung table",
                    f"model.timestep_minutes is {active} and the table says "
                    f"{table[rung]} for {rung}")
        else:
            rep.add(OK, "timestep vs the per-rung table",
                    f"{rung} at {active} min, kappa {model.get('filter_kappa')}")
    except Exception as exc:
        rep.add(WARN, "timestep vs the per-rung table", f"not checked: {exc}")

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
    # The physics filter, not the hyperdiffusion, is this model's small-scale
    # damping: two to three orders of magnitude stronger at every scale. Where
    # it overtakes the flow's own cascade is therefore where the model stops
    # resolving, and at gamma 8 that was 0.41 of the truncation -- most of the
    # rung being paid for and thrown away. This is a DEFAULT guard rather than a
    # diagnostic to remember: a filter that reaches too far down fails here.
    #
    # Crossover of 2*kappa*x^gamma/dt against the cascade rate x/tau_vorticity:
    #   x = (dt / (2 kappa tau))**(1/(gamma-1))
    try:
        model = config["model"]
        hd = model.get("hyperdiffusion") or {}
        rung = str(model.get("resolution", "")).upper()
        kappa = float(model["filter_kappa"])
        gamma = int(model["filter_power"])
        dt = float(model["timestep_minutes"]) * 60.0
        tau = float(hd["timescales_days"][rung]["vorticity"]) * 86400.0
        reach = (dt / (2.0 * kappa * tau)) ** (1.0 / (gamma - 1))
        floor = float(model.get("filter_confinement_floor", 0.60))
        if not model.get("physics_filter"):
            rep.add(OK, "filter confinement", "no filter configured")
        elif reach < floor:
            rep.add(FAIL, "filter confinement",
                    f"{rung} at kappa {kappa:g}, gamma {gamma}, dt {dt/60:g} min "
                    f"overtakes the cascade at {reach:.3f} of the truncation, "
                    f"below the {floor:.2f} floor: everything above that is "
                    f"resolution being damped away")
        else:
            rep.add(OK, "filter confinement",
                    f"{rung} damps from {reach:.3f} of the truncation "
                    f"(floor {floor:.2f})")
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
    try:
        manifest = ROOT / "exoplasim" / "binary_manifest.json"
        run_dir = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "run"
        on_disk = sorted(run_dir.glob("most_plasim_*.x")) if run_dir.is_dir() else []
        if not manifest.is_file():
            rep.add(FAIL, "binary manifest",
                    "absent; run exoplasim/scripts/rebuild_binaries.py")
        elif not on_disk:
            rep.add(WARN, "binaries", "none built")
        else:
            mf = json.loads(manifest.read_text(encoding="utf-8"))
            known = mf.get("binaries", {})
            # Source keys are relative to the package root, not by basename,
            # because two files of the same name can live under
            # different directories.
            src = ROOT / "vendor" / "exoplasim" / "exoplasim"
            bad = []
            for exe in on_disk:
                rec = known.get(exe.name)
                if rec is None:
                    bad.append(f"{exe.name}: not in manifest")
                    continue
                if rec["sha256"] != sha256_of(exe):
                    bad.append(f"{exe.name}: sha differs from manifest")
                    continue
                for name, want in (rec.get("sources") or {}).items():
                    got = sha256_of(src / name) if (src / name).is_file() else None
                    if got != want:
                        bad.append(f"{exe.name}: built from an older {name}")
            rep.add(FAIL if bad else OK, "binaries carry current patches",
                    "; ".join(bad) if bad else
                    f"{len(on_disk)} executables match the manifest")
    except Exception as exc:
        rep.add(WARN, "binaries", f"not checked: {exc}")

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

    # -- the flux-to-kelvin slope still matches the runs it was measured on ---
    try:
        import sensitivity
        sensitivity.test_identity()
        stale = sensitivity.verify()
        rep.add(FAIL if stale else OK, "flux-to-kelvin slope",
                "; ".join(stale) if stale else
                f"{sensitivity.SLOPE_K_PER_FLUX_RATIO} K per unit flux ratio, "
                f"reproduced from the runs lib/sensitivity.py names")
    except Exception as exc:
        rep.add(WARN, "flux-to-kelvin slope", f"not checked: {exc}")

    # -- config's derived radiation scalings against the artifacts that made --
    #    them
    #
    # Each of these is DERIVED: a generator computes it, prints it, and a human
    # retypes it into `config/planet.yaml`. That route has no check on it, and
    # it has already drifted twice elsewhere -- the CO2 shortwave fit in the
    # model source, and the dust notes against the albedos the code now reads.
    #
    # Compared to the FULL-PRECISION value in the artifact, at the tolerance the
    # config's own rounding implies: config carries three decimals, so half a
    # unit in the last place is the largest disagreement that can be honest
    # rounding, and anything beyond it is drift. That bound is arithmetic on the
    # written precision, not a number picked to make today's values pass.
    #
    # The generator is named in each row so a failure says what to re-run.
    try:
        model = config["model"]
        sw = ROOT / "exoplasim" / "analysis" / "shortwave_band_weights.json"
        cloud = ROOT / "exoplasim" / "analysis" / "cloud_band_weight.json"
        rows = [
            ("h2o_sw_weight", sw, ("h2o", "weight"),
             "exoplasim/scripts/shortwave_band_weights.py"),
            ("co2_sw_weight", sw, ("co2", "weight"),
             "exoplasim/scripts/shortwave_band_weights.py"),
            ("cloud_absorption_scale", cloud, ("weight", "central"),
             "exoplasim/scripts/cloud_band_weight.py"),
        ]
        problems, checked = [], []
        for key, path, where, generator in rows:
            if key not in model:
                problems.append(f"{key} is not in config/planet.yaml")
                continue
            if not path.is_file():
                problems.append(f"{key}: {rel(path)} is missing; run {generator}")
                continue
            node = json.loads(path.read_text(encoding="utf-8"))
            for step in where:
                node = node[step]
            declared = float(model[key])
            # Half a unit in the last decimal place the config actually wrote.
            written = str(model[key])
            places = len(written.split(".")[1]) if "." in written else 0
            tol = 0.5 * 10.0 ** (-places)
            if abs(declared - float(node)) > tol:
                problems.append(
                    f"{key} is {declared:g} in config but {float(node):.6g} in "
                    f"{rel(path)} ({'.'.join(where)}), past the {tol:g} its own "
                    f"{places} decimals allow: re-run {generator} or retype it")
            else:
                checked.append(f"{key} {declared:g}")
        rep.add(FAIL if problems else OK,
                "derived radiation scalings match their artifacts",
                "; ".join(problems) if problems else
                f"{len(checked)} compared: {', '.join(checked)}")
    except Exception as exc:
        rep.add(WARN, "derived radiation scalings match their artifacts",
                f"not checked: {exc}")
    rep.show()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
