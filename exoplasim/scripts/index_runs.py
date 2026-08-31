#!/usr/bin/env python3
"""The ledger of every run that has EXISTED, since a run id says nothing.

    python exoplasim/scripts/index_runs.py           # rescan and merge
    python exoplasim/scripts/index_runs.py --find flux_ratio=0.968
    python exoplasim/scripts/index_runs.py --register exoplasim/runs/run_xxx
    python exoplasim/scripts/index_runs.py --recordless

Run ids are UUIDs. That is deliberate -- see `run_id` in `run_exoplasim.py` --
and it moves the entire burden of "what was this run" onto the manifest each run
already writes. This reads those manifests and produces
`exoplasim/runs/INDEX.json`, the human- and machine-readable answer.

A LEDGER, NOT A DIRECTORY LISTING, and that distinction is the whole of
WORLD-WW6Z. This once rebuilt `INDEX.json` from a scan of `exoplasim/runs/` and
wrote whatever the scan returned. The identity of a run lives in a manifest
inside its own directory, which is untracked and dies with the payload, so the
index recorded what EXISTS rather than what has EXISTED: a run created and
deleted between two scans left nothing tracked anywhere, and a rescan after any
deletion erased the row of a run that had been indexed. Forty-seven runs reached
that state -- among them both arms of the T42 ladder comparison an audit rests
on -- and none of them can be identified now, let alone re-measured. CLAUDE.md
rule 6 says this file is the only record of what each run was; it can only be
that if a row is written at CREATION and never removed by a scan.

So two things hold the invariant, and both are needed:

  REGISTER AT CREATION   `run_exoplasim.py` calls `register()` as soon as it has
                         written the manifest, before the model integrates a
                         single orbit, and again when the block ends however it
                         ends. `continue_exoplasim.py` re-registers after each
                         segment. A run therefore has a tracked row from the
                         moment it exists, and the window in which a deletion
                         could erase it entirely is closed rather than narrowed.
  MERGE, NEVER REPLACE   a rescan updates the rows it can see and KEEPS every
                         row it cannot. A directory that has gone keeps its last
                         known row and gains `payload_present: false`, the date
                         the disappearance was first noticed, and whether an
                         `archive/runs/` stub holds its extracted products. The
                         row is the identity; the payload was only ever the
                         evidence.

`scripts/archive_runs.py` is the deliberate deletion path and writes an
`INDEX_ENTRY.json` stub from the row it is deleting. That stub is still the
place a dead run's kept files live. What changed is that it is no longer the
ONLY thing standing between a deleted run and a dead string, because a hand
`rm -rf`, a crash during preparation and a purge all bypass it.

Runs written before the UUID change carry their parameters in their directory
name instead, and most carry no `physical` block. Those are indexed from their
manifest where possible and flagged `legacy` where not, rather than being parsed
back out of the name: a name is not a record.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import sys as _sys
if str(ROOT / "lib") not in _sys.path:
    _sys.path.insert(0, str(ROOT / "lib"))
from paths import rel  # noqa: E402
RUNS = ROOT / "exoplasim" / "runs"
INDEX = RUNS / "INDEX.json"
ARCHIVE = ROOT / "archive" / "runs"
RECORDLESS = ARCHIVE / "RECORDLESS.json"

FIELDS = ["resolution", "layers", "ranks", "flux_ratio", "co2_ppm",
          "rotation_hours", "obliquity_degrees", "eccentricity", "glaciers",
          "stellar_spectrum", "geography"]

# `run_` and twelve hex digits is exactly what `run_id()` in `run_exoplasim.py`
# emits, so a string of this shape anywhere in the tree names a run that once
# existed under `exoplasim/runs/`. `--recordless` looks for these.
RUN_ID = re.compile(r"run_[0-9a-f]{12}")

# The two that are not runs: a docstring's illustrative id and a fixture id in a
# gate. Both happen to be well-formed, and neither ever named a directory. Named
# literally rather than filtered by a rule like "cited only by a .py file",
# because such a rule would silently drop a REAL run whose only surviving
# citation is in a script, which is the failure this enumeration exists to
# measure. Adding to this set is therefore a deliberate edit with a reason.
NOT_RUNS = {
    "run_1a2b3c4d5e6f",   # exoplasim/scripts/derive_design_flux.py, an example
    "run_aaaaaaaaaaaa",   # scripts/smoke_test.py, a fixture
}


def high_cadence(m: dict) -> dict:
    """The high-cadence stream orbit by orbit, and whether it was substituted.

    THE SUBSTITUTION IS WHY THIS IS INDEXED AT ALL. ExoPlaSim answers a failed
    high-cadence conversion by calling `integritycheck`, which re-runs pyburn
    under `example.nl` -- the REGULAR variable set and the regular twelve-bin
    average -- and returns success. The segment then reports normally and the
    run directory holds a twelve-bin average under the high-cadence name, which
    nothing about its name or its location distinguishes from the real product.
    Orbit 127 of run_67323a923013 is that file, and it is the accepted
    baseline.

    So the index carries it. Run ids are UUIDs and `exoplasim/runs/` is not
    tracked, which makes this the only durable record of what a run's
    high-cadence stream is; `substituted_conversion` names the file that was
    reported, kept on disk beside the corrected one because it is the evidence
    of what a consumer would have read.

    Derived from the manifest and never written here, so re-running the
    conversion moves this with it. `continue_exoplasim.py
    --reconvert-high-cadence` is what fills those manifest blocks in.
    """
    raws = {int(r["year_index"]): r for r in (m.get("high_cadence_raw") or [])}
    convs = {int(c["year_index"]): c
             for c in (m.get("high_cadence_conversions") or [])}
    orbits = {}
    for year in sorted(set(raws) | set(convs)):
        raw, conv = raws.get(year, {}), convs.get(year, {})
        orbits[str(year)] = {
            "raw": raw.get("file"),
            "raw_bytes": raw.get("bytes"),
            "raw_sha256": raw.get("sha256"),
            "expected_samples": raw.get("expected_samples_per_orbit"),
            "converted_samples": conv.get("samples"),
            "converted_sha256": conv.get("sha256"),
            "substituted_conversion": conv.get("moved_aside"),
        }
    return orbits


def read(run_dir: Path) -> dict | None:
    mf = run_dir / "run_manifest.json"
    if not mf.is_file():
        return None
    try:
        m = json.loads(mf.read_text(encoding="utf-8"))
    except Exception:
        return None

    phys = m.get("physical")
    legacy = phys is None
    if legacy:
        # Recover what the manifest does hold. source_config is a full copy of
        # the configuration the run was prepared with, so most of it survives
        # even for runs that predate the fingerprint block.
        cfg = m.get("source_config") or {}
        p, a, mo = (cfg.get("planet", {}), cfg.get("atmosphere", {}),
                    cfg.get("model", {}))
        phys = {
            "resolution": mo.get("resolution"),
            "layers": mo.get("layers"),
            "ranks": mo.get("ncpus"),
            # The run's ACTUAL flux, from derived_parameters -- not
            # orbit.baseline_flux_earth, which is the config's baseline at
            # prepare time and is simply a different number whenever a run was
            # launched with --flux-ratio. Reading the baseline here reported the
            # 0.968 bootstrap as 0.945.
            "flux_ratio": (m.get("derived_parameters", {}) or {})
                          .get("stellar_flux_ratio_earth"),
            "co2_ppm": (round(1e6 * float(a["pCO2_bar"]), 3)
                        if a.get("pCO2_bar") is not None else None),
            "rotation_hours": p.get("rotation_hours"),
            "obliquity_degrees": p.get("obliquity_degrees"),
            "eccentricity": p.get("eccentricity"),
            "glaciers": ((cfg.get("surface", {}) or {})
                         .get("glaciers", {}) or {}).get("enabled"),
            "stellar_spectrum": (cfg.get("radiation", {}) or {})
                                .get("stellar_spectrum"),
            "geography": None,
        }

    orbits = len(list(run_dir.glob("MOST.*.nc")))
    conv = m.get("convergence_assessment") or {}
    return {
        "run_id": m.get("run_id", run_dir.name),
        "directory": run_dir.name,
        "status": m.get("status"),
        "legacy_name": legacy,
        "orbits_on_disk": orbits,
        "snapshots_on_disk": len(list((run_dir / "snapshots").glob("*.nc")))
                             if (run_dir / "snapshots").is_dir() else 0,
        "size_gb": round(sum(f.stat().st_size for f in run_dir.rglob("*")
                             if f.is_file()) / 1e9, 2),
        "executable_sha256": ((m.get("executable") or {}).get("sha256")
                              or m.get("executable_sha256")),
        # WHICH ARM, if the run was one. Run ids are UUIDs, so this index is the
        # only record of what each run WAS, and an arm run -- another precision,
        # another flag line, a patched model source -- is not the shipped model
        # and must not read as it. Null for every ordinary run. world-u5pf.
        "arm": ((m.get("executable") or {}).get("arm") or {}).get("build_tag"),
        "config_sha256": m.get("config_sha256"),
        "physical": phys,
        "source_build": (m.get("source_config") or {}).get("source_build"),
        # Every climatology this run produced, keyed by label, so a consumer can
        # resolve a climatology back to the run and build that made it WITHOUT
        # globbing the analysis directory and taking whatever sorts last.
        "climatologies": {
            label: {"regular": Path(e.get("regular") or "").name,
                    "snapshots": Path(e.get("snapshots") or "").name,
                    "orbit_count": e.get("orbit_count")}
            for label, e in (m.get("climatologies") or {}).items()},
        "converged": conv.get("sufficiently_equilibrated_for_worldbuilding",
                              conv.get("pass")),
        "convergence_metrics": conv.get("metrics") or {},
        "high_cadence": high_cadence(m),
    }

NOTE = ("Every run that has EXISTED, not every run directory that is on disk. "
        "Run ids are UUIDs and carry no meaning, so this is where the physics "
        "lives, and by CLAUDE.md rule 6 it is the only record of what each run "
        "was. Rows are written at CREATION by run_exoplasim.py and are NEVER "
        "removed by a rescan: a run whose payload has gone keeps its row and "
        "carries payload_present false. Generated by "
        "exoplasim/scripts/index_runs.py; do not edit.")


def load(index: Path = INDEX) -> list[dict]:
    """The rows the ledger already holds. Empty when there is no ledger yet."""
    if not index.is_file():
        return []
    try:
        return json.loads(index.read_text(encoding="utf-8")).get("runs") or []
    except (OSError, json.JSONDecodeError):
        return []


def _gone(row: dict, today: str) -> dict:
    """A row whose directory is no longer on disk, marked and left otherwise.

    The last known row IS the identity: orbit count, executable sha, config sha,
    source build, climatologies. None of that is re-derivable once the manifest
    has gone with the payload, so nothing here recomputes a field. It records
    that the evidence behind them is no longer readable, and whether
    `archive/runs/` holds the extracted products.

    `payload_gone_since` is set once and preserved, so the date means the day
    the disappearance was first noticed and not the day of the latest rescan.
    """
    row = dict(row)
    row["payload_present"] = False
    row.setdefault("payload_gone_since", today)
    row["archived"] = (ARCHIVE / row["directory"]).is_dir()
    return row


def merge(previous: list[dict], scanned: list[dict]) -> list[dict]:
    """The scan laid over the ledger, keeping every row the scan cannot see.

    Keyed by DIRECTORY and not by run id, because the two are not the same key:
    a crashed preparation is moved aside to `<directory>_crashed` and keeps the
    run id it was prepared under, so two rows legitimately share one id.
    """
    today = datetime.now(timezone.utc).date().isoformat()
    rows = {r["directory"]: r for r in previous}
    for row in scanned:
        row = dict(row)
        row["payload_present"] = True
        row["archived"] = (ARCHIVE / row["directory"]).is_dir()
        was = rows.get(row["directory"], {})
        if was.get("payload_gone_since") and not was.get("payload_present", True):
            # A directory that came back: a restore, or a purge that took the
            # payload and left the identity. Drop the gone marks rather than
            # carry a contradiction, and say once in the row that it happened.
            row["payload_returned"] = today
        rows[row["directory"]] = row
    seen = {r["directory"] for r in scanned}
    for name, row in rows.items():
        if name not in seen:
            rows[name] = _gone(row, today)
    return [rows[k] for k in sorted(rows)]


def write(rows: list[dict], index: Path = INDEX) -> None:
    """Write the ledger atomically, so a crash cannot truncate the only record."""
    payload = {
        "note": NOTE,
        "generated": datetime.now(timezone.utc).isoformat(),
        "runs": rows,
    }
    index.parent.mkdir(parents=True, exist_ok=True)
    tmp = index.with_name(index.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, index)


def register(run_dir: Path, index: Path = INDEX) -> dict | None:
    """Put this run in the ledger NOW, and return the row.

    CALLED AT CREATION, from `run_exoplasim.py` the moment the manifest exists
    and again when the block ends however it ends, and from
    `continue_exoplasim.py` after each segment. That is what makes this a record
    of what has existed: a scan can only ever see a run that is still on disk
    when someone happens to run it, and the runs this repository has lost were
    lost in exactly that gap.

    Merging rather than appending, so re-registering the same run updates its
    row and leaves every other row alone. Under `flock`, because two arms of a
    pair are two processes writing this one file and a lost update here is a
    lost run.
    """
    import fcntl

    index.parent.mkdir(parents=True, exist_ok=True)
    lock = index.with_name(index.name + ".lock")
    with open(lock, "w", encoding="utf-8") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        row = read(run_dir)
        if row is None:
            return None
        write(merge(load(index), [row]), index)
    return row


def recordless(root: Path = ROOT) -> list[dict]:
    """Every run id the tree cites that no record can identify.

    THE MEASURE OF THE DEFECT THIS FILE EXISTS TO CLOSE. A run id is twelve hex
    digits and means nothing on its own, so an id cited by a note, an audit or
    an analysis product is readable only through a record. This collects the ids
    that have none -- not in the ledger, not an `INDEX_ENTRY.json` stub, not an
    `INDEX_AT_DELETION.json` batch, not a `RECONSTRUCTED.json` -- and names the
    files that cite each, which is the whole of what survives of it.

    Generated rather than listed, so the enumeration shrinks by itself as
    identity is recovered and cannot go stale the way a hand-kept list does.
    """
    known: set[str] = set()
    for row in load():
        known |= {str(row.get("run_id") or ""), str(row.get("directory") or "")}
    if ARCHIVE.is_dir():
        for path in sorted(ARCHIVE.rglob("*.json")):
            # Not its own output. RECORDLESS.json lists the ids that have NO
            # record, so reading it as a record makes every id known and the
            # enumeration collapses to empty on its second run.
            if path == RECORDLESS:
                continue
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            entries = (loaded.get("runs", [loaded]) if isinstance(loaded, dict)
                       else loaded)
            for entry in entries if isinstance(entries, list) else []:
                if isinstance(entry, dict):
                    known |= {str(entry.get("run_id") or ""),
                              str(entry.get("directory") or "")}
        known |= {d.name for d in ARCHIVE.iterdir() if d.is_dir()}
    known.discard("")

    # `.beads` holds issue text, `runs` and `archive` are the records
    # themselves, and `vendor` is upstream source. None of them cites a run the
    # way a note or an analysis product does.
    skip = {".git", ".beads", "vendor", ".venv", "runs", "archive", ".claude"}
    cited: dict[str, list[str]] = {}
    for path in sorted(root.rglob("*")):
        rel_parts = path.relative_to(root).parts
        if set(rel_parts) & skip or path.is_symlink() or not path.is_file():
            continue
        if path.suffix not in {".json", ".md", ".py", ".yaml", ".txt"}:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for rid in set(RUN_ID.findall(text)):
            if rid not in known and rid not in NOT_RUNS:
                cited.setdefault(rid, []).append("/".join(rel_parts))
    return [{"run_id": rid, "cited_by": sorted(files)}
            for rid, files in sorted(cited.items())]


def write_recordless(rows: list[dict]) -> None:
    RECORDLESS.parent.mkdir(parents=True, exist_ok=True)
    RECORDLESS.write_text(json.dumps({
        "note": "Run ids this tree cites that NO record can identify: not in "
                "exoplasim/runs/INDEX.json, not an INDEX_ENTRY.json stub, not "
                "an INDEX_AT_DELETION.json batch, not a RECONSTRUCTED.json. "
                "Each was a directory under exoplasim/runs/ -- the id shape is "
                "what run_exoplasim.py emits and nothing else emits it -- and "
                "its configuration, executable sha, source build and orbit "
                "count went with its manifest. The files beside each id are the "
                "whole of what survives of it, and a claim resting on one of "
                "these runs rests on a run nobody can check. WORLD-WW6Z; "
                "notes/audits/resolution-ladder.md is the case that found it. "
                "Generated by exoplasim/scripts/index_runs.py --recordless.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "count": len(rows),
        "runs": rows,
    }, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--find", action="append", default=[],
                    help="filter as key=value, repeatable, e.g. resolution=T42")
    ap.add_argument("--register", type=Path, metavar="RUN_DIR",
                    help="write one run directory's row into the ledger and "
                         "exit. What run_exoplasim.py calls at creation")
    ap.add_argument("--recordless", action="store_true",
                    help="rewrite archive/runs/RECORDLESS.json: every run id "
                         "the tree cites that no record can identify")
    ap.add_argument("--output", type=Path, default=INDEX)
    args = ap.parse_args()

    if args.register:
        row = register(args.register.resolve(), args.output)
        if row is None:
            raise SystemExit(f"{rel(args.register)} has no run_manifest.json, "
                             "so there is nothing to register. A run is "
                             "registered once its manifest is written.")
        print(f"registered {row['directory']} in {rel(args.output)}")
        return

    if args.recordless:
        rows = recordless()
        write_recordless(rows)
        for row in rows:
            print(f"  {row['run_id']}  cited by {len(row['cited_by'])} file(s)")
        print(f"\n{len(rows)} run ids no record can identify; "
              f"wrote {rel(RECORDLESS)}")
        return

    scanned = [r for r in (read(d) for d in sorted(RUNS.iterdir()) if d.is_dir())
               if r is not None]
    rows = merge(load(args.output), scanned)

    shown = rows
    for f in args.find:
        k, _, v = f.partition("=")
        shown = [r for r in shown if str(r["physical"].get(k)).lower() == v.lower()]

    print(f"{'directory':66}{'res':5}{'flux':>7}{'orb':>5}{'GB':>7}  conv")
    for r in shown:
        p = r["physical"]
        flux = p.get("flux_ratio")
        gone = "" if r.get("payload_present", True) else "  PAYLOAD GONE"
        print(f"{r['directory'][:64]:66}{str(p.get('resolution') or '-'):5}"
              f"{(f'{flux:.3f}' if isinstance(flux,(int,float)) else '-'):>7}"
              f"{r['orbits_on_disk']:>5}{r['size_gb']:>7.2f}  "
              f"{'yes' if r['converged'] else ''}{gone}")
    here = [r for r in rows if r.get("payload_present", True)]
    print(f"\n{len(shown)} of {len(rows)} rows, {len(here)} with a payload on "
          f"disk, {sum(r['size_gb'] for r in here):.1f} GB")

    write(rows, args.output)
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
