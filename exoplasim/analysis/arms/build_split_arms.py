#!/usr/bin/env python
"""Build the two arms that split the batch-2 bundle's measured drift.

WORLDBUILDING CONTEXT: Vesper is a fictional super-Earth and both arms below are
configurations of the climate model that simulates it.

WHY THESE TWO. Every arm branched from `run_14906cb7b914` relaxes at least
12.5 K colder than that donor, against registered predictions summing to about
+1 K of warming. The drift confounds the compiled forcing terms with the only two
config changes the donor did not carry: the four `surface.cryosphere` constants,
and `surface.land_water_column` going from a 1-layer bucket to a 2-layer scheme.
Both are declared in `config/planet.yaml`, so both are arm-able against the same
control on the same restart, and that splits the drift into the part config
carries and the part the source carries. It is as far as the split can go without
a control binary, which was not preserved.

These two edit the `surface` block rather than `model`, which is why they are
here and not in `build_arm_config.py`: that one exists for A3 arms, where both
sides differ by a namelist key and nothing else.
"""
import pathlib
import re

import yaml

HERE = pathlib.Path(__file__).resolve().parent
CONTROL = HERE / "s8rv_10.yaml"


def block_bounds(lines: list[str], name: str) -> tuple[int, int]:
    start = next(i for i, l in enumerate(lines) if l.startswith(f"  {name}:"))
    end = next(i for i in range(start + 1, len(lines)) if re.match(r"^  [a-z_]+:", lines[i]))
    return start, end


def build_bucket(out: pathlib.Path) -> None:
    """The land water column back to the 1-layer bucket the donor ran."""
    lines = CONTROL.read_text(encoding="utf-8").splitlines(keepends=True)
    start, end = block_bounds(lines, "land_water_column")
    replacements = {
        "    scheme:": "    scheme: bucket\n",
        "    layers:": "    layers: 1\n",
        "    layer_thickness_m:": "    layer_thickness_m: [1.5]\n",
    }
    for i in range(start + 1, end):
        for key, value in replacements.items():
            if lines[i].startswith(key):
                lines[i] = value
    lines.insert(start + 1, "    # ARM OVERRIDE: the 1-layer bucket run_14906cb7b914 carried.\n")
    out.write_text("".join(lines), encoding="utf-8")
    got = yaml.safe_load(out.read_text(encoding="utf-8"))["surface"]["land_water_column"]
    assert got["scheme"] == "bucket" and got["layers"] == 1, got
    assert got["layer_thickness_m"] == [1.5], got
    print("bucket arm:", {k: got[k] for k in
                          ("scheme", "layers", "layer_thickness_m")})


def build_no_cryosphere(out: pathlib.Path) -> None:
    """The whole `cryosphere` block absent, which is what the donor's config had.

    The block goes rather than its four values, because `run_14906cb7b914`'s
    `source_config` carries no `surface.cryosphere` key at all and an emptied
    block parses as None, which is a third state and not the one being restored.
    """
    lines = CONTROL.read_text(encoding="utf-8").splitlines(keepends=True)
    start, end = block_bounds(lines, "cryosphere")
    kept = lines[:start] + lines[end:]
    out.write_text("".join(kept), encoding="utf-8")
    surface = yaml.safe_load(out.read_text(encoding="utf-8"))["surface"]
    assert "cryosphere" not in surface, "the cryosphere block survived the removal"
    print("no-cryosphere arm: surface keys", sorted(surface))


if __name__ == "__main__":
    build_bucket(HERE / "split_bucket.yaml")
    build_no_cryosphere(HERE / "split_no_cryosphere.yaml")
