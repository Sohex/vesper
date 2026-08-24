#!/usr/bin/env python3
"""Cut a restart TEMPLATE from a run of the target build.

Worldbuilding frame: the restart is the saved state of the Vesper climate
model. Nothing here is about the real world.

    python exoplasim/scripts/build_restart_template.py --from-run RUN_DIR

The output name is derived from the run's own geometry and thread count, so a
template cannot be filed under a rung it is not.

WHAT A TEMPLATE IS FOR. `convert_restart.py` starts from a template and
overlays only the state whose policy permits transfer, so the template is what
supplies everything the DONOR must not: the record set and its order, the real
width, the target-resolution static surface fields, the compiler's random-seed
shape, and a clean accumulation window. It is the reason a converted restart
cannot carry a donor's stale roughness or albedo into a run that staged new
ones.

WHY IT IS CUT FROM A RUN RATHER THAN WRITTEN HERE. Only the target executable
knows its own record set, and only a run of it against the target's own staged
`.sra` files has the target's static fields in them. Writing a template in
Python would mean keeping a second copy of the model's cold defaults in step
with the model by hand, and the four accumulators whose clean value is not zero
say how that ends.

THE ONE THING THIS ADDS to the run's own restart is a clean accumulation
window, through `reset_restart_accumulators.py`, whose values come from
`restart_schema.py` and are checked against `outreset` by
`scripts/smoke_test.py`. Doing it here rather than asking the run to stop on an
output boundary is what makes the recipe independent of where the run happened
to end.

WHAT IS NOT CLEAN, and is recorded rather than pretended away: the prognostic
state and the derived surface fields are the run's, at whatever point it
stopped. The converter overwrites every prognostic record from the donor, and
it names the derived ones for the model's own post-load fixup, so neither
reaches a converted run unexamined. What a template is NOT is a cold start.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import reset_restart_accumulators  # noqa: E402
import restart_format  # noqa: E402
import restart_schema  # noqa: E402
from _paths import COMPONENT_ROOT, MODEL_SRC, PROJECT_ROOT  # noqa: E402

MODEL_SOURCE = MODEL_SRC / "plasim" / "src"
TEMPLATES = COMPONENT_ROOT / "inputs" / "templates"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(path: Path) -> dict:
    """Everything about a template a caller has to be able to check.

    Returns the geometry, and the accumulator verdict per record. Raises if the
    file is not a restart at all; a template that is merely UNCLEAN comes back
    described rather than refused, because saying which record is dirty is more
    use than refusing to look.
    """
    records = restart_format.read(path)
    geometry, real_bytes = restart_schema.describe(records)
    inventory = restart_schema.inventory_from_source(MODEL_SOURCE)
    resets = restart_schema.model_resets_from_source(MODEL_SOURCE)

    unknown = sorted({r.name for r in records} - set(restart_schema.POLICY))
    dirty, clean, uncheckable = [], 0, []
    for rec in records:
        pol = restart_schema.POLICY.get(rec.name)
        if pol is None or pol.semantic != restart_schema.ACCUMULATOR:
            continue
        kind, value = restart_schema.derived_model_reset(rec.name, inventory, resets)
        if kind == "zero":
            want = b"\x00" * rec.nbytes
        elif kind == "sentinel":
            want = reset_restart_accumulators.clean_payload(rec, value, real_bytes)
        else:
            uncheckable.append(rec.name)
            continue
        if rec.payload == want:
            clean += 1
        else:
            dirty.append(rec.name)
    return {"geometry": geometry.label, "real_bytes": real_bytes,
            "records": len(records), "nesp": geometry.nesp,
            "nseedlen": geometry.nseedlen, "unknown_records": unknown,
            "accumulators_clean": clean, "accumulators_dirty": dirty,
            "accumulators_the_model_never_resets": uncheckable}


def provenance(run_dir: Path, source: Path, output: Path, report: dict) -> dict:
    """Who made this template, from what, and against which surface."""
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) \
        if manifest_path.is_file() else {}
    return {
        "schema_version": 1,
        "kind": "restart_template",
        "generated": datetime.now(timezone.utc).isoformat(),
        "generator": "exoplasim/scripts/build_restart_template.py",
        "generator_sha256": sha256(Path(__file__)),
        "schema_sha256": sha256(Path(restart_schema.__file__)),
        "output": {"path": str(output), "sha256": sha256(output)},
        "cut_from": {
            "run_dir": str(run_dir),
            "run_id": manifest.get("run_id"),
            "restart": str(source),
            "restart_sha256": sha256(source),
            "completed_orbits": manifest.get("completed_orbits"),
            "status": manifest.get("status"),
        },
        # The template's authority is exactly the executable that wrote it: its
        # record set, its real width and its seed shape are properties of this
        # binary and of no other.
        "executable": manifest.get("executable"),
        "physical": manifest.get("physical"),
        "source_build": manifest.get("source_build"),
        "config_sha256": manifest.get("config_sha256"),
        # What the target's static surface records ARE. A conversion is only
        # sound onto a run whose staged surface matches these.
        "surface_field_sha256": manifest.get("surface_field_sha256"),
        "surface_fields": manifest.get("surface_fields"),
        "template": report,
        "note": ("A clean accumulation window on the target executable's own "
                 "record set and static surface. The prognostic and derived "
                 "records are the donor run's and are overwritten or flagged "
                 "by convert_restart.py; this is not a cold start."),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-run", type=Path,
                    help="a run directory of the TARGET build, whose restart "
                         "the template is cut from")
    ap.add_argument("--restart", type=Path,
                    help="the restart inside it; defaults to plasim_restart")
    ap.add_argument("--output", type=Path,
                    help=f"where the template is written; defaults into "
                         f"{TEMPLATES.relative_to(PROJECT_ROOT)}/")
    ap.add_argument("--audit", type=Path,
                    help="report whether an existing template is clean, and "
                         "write nothing")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing template")
    args = ap.parse_args()

    if args.audit:
        report = audit(args.audit)
        print(json.dumps(report, indent=2))
        return 1 if report["accumulators_dirty"] or report["unknown_records"] else 0

    if args.from_run is None:
        ap.error("--from-run is required (or --audit)")
    source = args.restart or (args.from_run / "plasim_restart")
    if not source.is_file():
        raise SystemExit(f"{source} is not there; a template is cut from a "
                         "run that actually wrote a restart")

    before = restart_format.read(source)
    geometry, _ = restart_schema.describe(before)
    manifest_path = args.from_run / "run_manifest.json"
    physical = (json.loads(manifest_path.read_text(encoding="utf-8"))
                .get("physical", {}) if manifest_path.is_file() else {})
    # No parallel-mode word: world-38b left one parallel mode, so it named
    # nothing. The templates already on disk keep the name they were written
    # with; the .provenance.json beside each is what identifies it.
    default_name = (f"{geometry.label}_l{geometry.nlev}"
                    f"_p{physical.get('ranks', 'x')}.rest")
    output = args.output or (TEMPLATES / default_name)
    if output.exists() and not args.force:
        raise SystemExit(f"{output} exists; pass --force to replace it")
    output.parent.mkdir(parents=True, exist_ok=True)

    info = reset_restart_accumulators.reset(source, output)
    report = audit(output)
    if report["accumulators_dirty"]:
        raise SystemExit(
            "the template is not clean after the reset: "
            f"{', '.join(report['accumulators_dirty'])}. That is a defect in "
            "reset_restart_accumulators or in the schema's reset column, not "
            "something to write out and hope about.")
    if report["unknown_records"]:
        raise SystemExit(
            f"the restart holds records the schema does not name: "
            f"{', '.join(report['unknown_records'])}. Add them to "
            "restart_schema.POLICY first; see exoplasim/README.md.")

    side = output.with_suffix(output.suffix + ".provenance.json")
    side.write_text(json.dumps(
        provenance(args.from_run, source, output, report), indent=2) + "\n")

    print(f"wrote {output}")
    print(f"  {report['geometry']} at {report['real_bytes']}-byte reals, "
          f"{report['records']} records, NESP {report['nesp']}, "
          f"seed {report['nseedlen']} integers")
    print(f"  {report['accumulators_clean']} accumulators at their clean value "
          f"({len(info['sentinels'])} of them nonzero), "
          f"{len(report['accumulators_the_model_never_resets'])} the model "
          "resets nowhere")
    print(f"wrote {side}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
