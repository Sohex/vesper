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

TWO WAYS IN, because a namelist key is not always a `model.*` key. `--overrides`
moves keys inside the `model` block, which is where most of them are declared;
`--set` takes DOTTED PATHS from the document root, for the namelist keys that
are declared elsewhere -- `surface.soil_thermal`, `surface.land_water_column`,
`surface.soil_albedo_moisture` are all groups of `landmod_nl` keys and all of
them are legitimate A3 arms. `--set` also INSERTS a path the base config does
not carry, which is what a block declaring nothing and running the model's
compiled default needs before an arm can move it.
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


def set_path(lines: list[str], path: str, value, note: str) -> None:
    """Set one DOTTED path in the config text, inserting the levels it lacks.

    Line-preserving for the reason the module docstring gives: the derivation of
    every value the arm does NOT move lives in the comment above it, and a YAML
    round trip deletes all of them.

    A path whose parent block is absent is inserted rather than refused. A block
    that declares nothing is a block running the model's compiled default, and
    that is exactly the state an arm has to be able to move off.
    """
    parts = path.split(".")
    start, end = 0, len(lines)
    for depth, key in enumerate(parts):
        indent = "  " * depth
        here = re.compile(rf"^{indent}{re.escape(key)}:")
        found = next((i for i in range(start, end) if here.match(lines[i])), None)
        if found is None:
            # The note marks the block, not every leaf under it: a comment
            # repeated above each of four values is noise, and the block it
            # sits on is what a reader needs told is an arm's and not the
            # base config's.
            block = ([f"{indent}# ARM OVERRIDE ({note}). Absent from the base config.\n"]
                     if depth < len(parts) - 1 else [])
            for d in range(depth, len(parts)):
                pad = "  " * d
                block.append(f"{pad}{parts[d]}:"
                             + (f" {scalar(value)}\n" if d == len(parts) - 1 else "\n"))
            lines[start:start] = block
            return
        if depth == len(parts) - 1:
            lines[found] = f"{indent}{key}: {scalar(value)}\n"
            return
        # Descend. The span ends at the next line indented no deeper than this
        # key, which is the first sibling or uncle -- anything past it belongs
        # to another block and a match there would be the wrong key.
        sibling = re.compile(rf"^ {{0,{len(indent)}}}[A-Za-z_0-9]+:")
        stop = next((i for i in range(found + 1, end) if sibling.match(lines[i])), end)
        start, end = found + 1, stop


def read_path(config: dict, path: str):
    node = config
    for part in path.split("."):
        if not isinstance(node, dict) or part not in node:
            return "<absent>"
        node = node[part]
    return node


def build(overrides: dict, out_path: pathlib.Path, note: str,
          paths: dict | None = None) -> tuple[pathlib.Path, str]:
    """Write an arm config with the named keys overridden, and verify each landed."""
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

    for path, value in (paths or {}).items():
        set_path(lines, path, value, note)

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
    for path, value in (paths or {}).items():
        got = read_path(config, path)
        if got != value:
            raise SystemExit(f"override {path} did not land: wanted {value!r}, got {got!r}")
    return out_path, hashlib.sha256(out_path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--note", required=True, help="what this arm is, for the header comment")
    parser.add_argument("--overrides", default="{}", help="JSON object of model.* keys")
    parser.add_argument("--set", dest="paths", default="{}",
                        help="JSON object of DOTTED paths from the document "
                             "root, for namelist keys declared outside `model`. "
                             "A path the base config lacks is inserted.")
    args = parser.parse_args()
    if json.loads(args.overrides) == {} and json.loads(args.paths) == {}:
        parser.error("nothing to override: pass --overrides, --set, or both. "
                     "An arm config identical to the base is the base.")
    path, sha = build(json.loads(args.overrides), args.out.resolve(), args.note,
                      json.loads(args.paths))
    print(f"{path} {sha}")


if __name__ == "__main__":
    sys.exit(main())
