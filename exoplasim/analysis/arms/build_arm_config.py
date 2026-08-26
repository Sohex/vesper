#!/usr/bin/env python
"""Build an ARM CONFIG from `config/planet.yaml` by overriding `model.*` keys.

WORLDBUILDING CONTEXT: Vesper is a fictional super-Earth and every arm below is
a configuration of the climate model that simulates it.

WHY TEXT AND NOT A YAML ROUND TRIP. Every value in `config/planet.yaml` carries
its derivation in the comment above it, and a round trip through `yaml.safe_load`
and `yaml.safe_dump` deletes all of it. The arm config is what
`run_exoplasim.py` hashes into `config_sha256` and what the run manifest points
at, so it is the record of what the arm integrated; a record stripped of every
derivation is not one. Overriding by line keeps the whole document and changes
only the values named.

`config/planet.yaml` IS NEVER WRITTEN. An arm config is a derived sibling, and
the keys it moves are the arm's own -- A3's condition that both arms differ by a
namelist key and nothing else is enforced by both arms being built from this one
base with a single override set each.
"""
import argparse
import hashlib
import json
import pathlib
import re
import sys

import yaml

HERE = pathlib.Path(__file__).resolve()
REPO = HERE.parents[3]
BASE = REPO / "config" / "planet.yaml"


def scalar(value) -> str:
    """Render one scalar for a YAML mapping value.

    NOT `yaml.safe_dump`: dumping a bare scalar emits a whole DOCUMENT, which
    carries a trailing `...` end marker, and splicing that into the middle of a
    mapping makes the rest of the file a second document that `safe_load`
    refuses. Caught by this module's own re-read check on the first arm built.
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (int, float)):
        return repr(value)
    return yaml.safe_dump(value, default_flow_style=True).strip().rstrip("...").strip()


def build(overrides: dict, out_path: pathlib.Path, note: str) -> tuple[pathlib.Path, str]:
    """Write an arm config with `model.*` overridden, and verify each landed."""
    lines = BASE.read_text(encoding="utf-8").splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines) if line.startswith("model:"))
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if re.match(r"^[a-z_]+:", lines[i]):
            end = i
            break

    remaining = dict(overrides)
    for i in range(start + 1, end):
        m = re.match(r"^  ([a-z_0-9]+):", lines[i])
        if m and m.group(1) in remaining:
            key = m.group(1)
            lines[i] = f"  {key}: {scalar(remaining.pop(key))}\n"

    inserted = [f"  # ARM OVERRIDE ({note}). Keys absent from the base config.\n"]
    for key, value in remaining.items():
        inserted.append(f"  {key}: {scalar(value)}\n")
    if len(inserted) > 1:
        lines[start + 1:start + 1] = inserted

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(lines), encoding="utf-8")

    # The check that can fail: re-read the file and assert every override landed
    # with the value asked for. A silent miss would make the two arms differ by
    # nothing and the A/B measure scatter.
    config = yaml.safe_load(out_path.read_text(encoding="utf-8"))
    for key, value in overrides.items():
        got = config["model"].get(key, "<absent>")
        if got != value:
            raise SystemExit(f"override {key} did not land: wanted {value!r}, got {got!r}")
    return out_path, hashlib.sha256(out_path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--note", required=True, help="what this arm is, for the header comment")
    parser.add_argument("--overrides", required=True, help="JSON object of model.* keys")
    args = parser.parse_args()
    path, sha = build(json.loads(args.overrides), args.out.resolve(), args.note)
    print(f"{path} {sha}")


if __name__ == "__main__":
    sys.exit(main())
