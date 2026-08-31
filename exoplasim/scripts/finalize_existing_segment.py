#!/usr/bin/env python3
"""Validate and record a completed segment after post-run bookkeeping failed.

TWO THINGS FAIL TOGETHER when a run is killed mid-block, and this recovers both.
`run_exoplasim.py` writes the segment record AND the derived-constant stamp from
inside the same `if args.run_years:` block, so a run that does not reach the end
of that block has neither. The segment was recovered here from the start; the
stamp was not, and `check_consistency.py`'s "runs vs the constants the model
derived" is what noticed -- a run that staged `VDIFF_LAMM`, `RCRITWIDTH` or
`GAMMA` at its negative sentinel records only that the model chose, and without
the stamp nobody can say which number it integrated. WORLD-ET25.

The stamp is read from the run's own `MOST_DIAG`, which survives a kill, so this
is recovery and not reconstruction. `--constants-only` does that half alone, for
a run whose segment record is intact and whose stamp is missing.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import yaml

from _paths import CONFIG

from continue_exoplasim import validate_year, year_diagnostics
from run_exoplasim import read_derived_constants, stellar_spectrum_digest
from segments import SEGMENT_PURPOSES


def record_derived_constants(run_dir: Path) -> str:
    """Stamp `derived_model_constants` from the run's own diag. Idempotent.

    Returns a line saying what it did. Recovery, not reconstruction: every value
    is read back out of the model's own initialisation print, which is the same
    source `run_exoplasim.py` reads.
    """
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("derived_model_constants"):
        return f"{run_dir.name}: derived_model_constants already recorded"
    manifest["derived_model_constants"] = read_derived_constants(run_dir)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n",
                             encoding="utf-8")
    named = ", ".join(f"{k} {v['value']:g} ({v['branch']})" for k, v
                      in manifest["derived_model_constants"].items())
    return f"{run_dir.name}: recovered derived_model_constants -- {named}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    # `--constants-only` recovers the derived-constant stamp and nothing else,
    # so the two halves of a failed bookkeeping are separately recoverable: a
    # run whose segment IS recorded and whose stamp is not must not have a
    # second segment appended to reach the stamp.
    parser.add_argument("--constants-only", action="store_true",
                        help="record only what the model derived for the "
                             "sentinel-selected constants, from this run's own "
                             "MOST_DIAG, and leave the segments alone. For a "
                             "run whose segment record is intact")
    parser.add_argument("--start-year", type=int, default=None)
    parser.add_argument("--end-year", type=int, default=None)
    parser.add_argument("--seasonal-output", action="store_true")
    # Declared, not inferred, for the same reason continue_exoplasim.py declares
    # it: seasonal output plus a run past its cutoff was taken to mean
    # climatology input, and that mislabelled a low-I/O verification segment.
    parser.add_argument("--purpose", default=None, choices=SEGMENT_PURPOSES,
                        help="what the recovered segment was for; see "
                             "continue_exoplasim.py --purpose")
    # Tri-state on purpose. Omitting both leaves the key out, and an orbit with
    # no `low_io` key is treated as low-I/O by build_climatology.py, which is
    # the honest default for a segment nobody has said anything about. Asserting
    # it clean has to be an act.
    parser.add_argument("--low-io", dest="low_io", action="store_true",
                        default=None,
                        help="the segment ran with PlaSim's low-I/O accumulation on")
    parser.add_argument("--no-low-io", dest="low_io", action="store_false",
                        help="the segment ran with it off; without either flag "
                             "the orbits are treated as low-I/O")
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    if args.constants_only:
        raise SystemExit(record_derived_constants(run_dir))
    for name in ("start_year", "end_year", "purpose"):
        if getattr(args, name) is None:
            raise SystemExit(f"--{name.replace('_', '-')} is required unless "
                             "--constants-only is given")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    diagnostics = []
    for year in range(args.start_year, args.end_year + 1):
        validate_year(
            run_dir / f"MOST.{year:05d}.nc",
            expected_times=int(model["regular_output_bins_per_orbit"]),
        )
        if args.seasonal_output:
            validate_year(
                run_dir / "snapshots" / f"MOST_SNAP.{year:05d}.nc",
                expected_times=int(model["seasonal_samples_per_orbit"]),
            )
        diagnostics.append(year_diagnostics(run_dir / f"MOST.{year:05d}.nc"))

    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    existing = manifest.setdefault("segments", [])
    if any(
        segment["start_year_index"] == args.start_year
        and segment["end_year_index"] == args.end_year
        for segment in existing
    ):
        raise RuntimeError("Segment is already recorded in the manifest")
    if args.purpose == "post_equilibrium_climatology" and not args.seasonal_output:
        raise SystemExit(
            "--purpose post_equilibrium_climatology without --seasonal-output: "
            "the orbital phase a climatology needs is only in the snapshots.")
    record = {
        "start_year_index": args.start_year,
        "end_year_index": args.end_year,
        "seasonal_output": args.seasonal_output,
        "purpose": args.purpose,
        # The spectrum as it is on disk NOW. This is a recovery path, run right
        # after the segment it records, so that is the one the orbits ran on --
        # but it is a weaker statement than the same key written by
        # continue_exoplasim.py, which writes it before the model starts.
        "stellar_spectrum_digest": stellar_spectrum_digest(config),
        "recovered_after_bookkeeping_failure_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostics": diagnostics,
    }
    if args.low_io is not None:
        record["low_io"] = bool(args.low_io)
    existing.append(record)
    manifest["completed_orbits"] = args.end_year + 1
    # A diagnostic segment does not move the run's status; see
    # continue_exoplasim.py.
    if args.purpose == "post_equilibrium_climatology":
        manifest["status"] = "climatology_complete"
    elif args.purpose == "spinup":
        manifest["status"] = "spinup_in_progress"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    # THE OTHER HALF OF THE SAME FAILURE. See the module docstring: the runner
    # writes the segment and this stamp from one block, so a recovery that
    # restores only the segment leaves the run unable to say what it integrated.
    print(record_derived_constants(run_dir))
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
