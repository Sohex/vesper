#!/usr/bin/env python3
"""Compare two or more World Orogen exports: terrain, composition, basement.

    python scripts/compare_exports.py source/precarve-craton/exoplasim-T42 /tmp/probe/raw
    python scripts/compare_exports.py --glim A/exoplasim-<rung> B/exoplasim-<rung>

Takes export directories -- anything with a `manifest.json` and `raw/` beside it,
including the lean `--no-grid --only geometry,elevation,basins,lithology` probes
that a generator change is tested with.

## Why this is a tool

Because testing a generator change means generating a comparison build and
measuring it, and this session that measurement was rewritten from scratch five
times: for arc erodibility, arc cover thickness, the craton multiplier, the basin
factor, and plate count. Four of those five refuted a hypothesis that had looked
solid in argument. A build takes one to three minutes and settles what reasoning
about this generator's internals does not, so the friction of measuring should be
as close to zero as it can be.

## Surface against basement, and why both

`substrate_class` is what is exposed; `basement_rock` is what lies under the
cover. A rule can move basement substantially while the surface does not move at
all, because the cells it claims carry thick cover -- which is exactly what the
craton investigation found after three sweeps had been read as null results off
the surface figures alone. Print both, always.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
from orogen import LAND  # noqa: E402

# Hartmann & Moosdorf (2012) Table 1, first level, global, as percentages of
# land. Grouped to Orogen's classes; see notes/audits/orogen-lithology.md for
# what each grouping does and does not capture.
GLIM = {"plutonic": 6.8, "metamorphic": 13.0, "volcanic": 6.2,
        "siliciclastic": 30.9, "carbonate": 7.8, "evaporite": 0.3}
GROUPS = {
    "plutonic": ["granite", "granodiorite"],
    "metamorphic": ["gneiss", "schist", "quartzite", "melange"],
    "volcanic": ["morb", "oib", "flood_basalt", "arc_basalt", "arc_andesite",
                 "rift_bimodal"],
    "siliciclastic": ["shelf_clastic", "foreland_clastic", "continental_clastic",
                      "pelagic"],
    "carbonate": ["carbonate"],
    "evaporite": ["evaporite"],
}


def load(root: Path) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    raw = manifest["raw"]["fields"]
    names = raw if isinstance(raw, dict) else {f["name"]: f for f in raw}

    def field(name: str):
        spec = names[name]
        return np.fromfile(root / spec["path"], dtype=spec["dtype"])

    area = field("cell_area").astype(float)
    land = field("surface_class") == LAND
    elev = field("elevation_km").astype(float)
    classes = {c["id"]: c["code"]
               for c in manifest["lithology"]["rockClasses"]}
    total = area[land].sum()

    def share(values, code_field):
        return {code: 100 * area[land & (code_field == i)].sum() / total
                for i, code in classes.items()}

    # `substrate_class` was named `surface_rock` before 2026-08-16; accept both
    # so an archived or older export still compares.
    top = "substrate_class" if "substrate_class" in names else "surface_rock"
    return {
        "hash": manifest["hashes"]["finalElevation"][:16],
        "land_pct": 100 * total / area.sum(),
        "mean_km": float(np.average(elev[land], weights=area[land])),
        "max_km": float(elev[land].max()),
        "basins": (manifest.get("basins", {}).get("counts", {})
                   .get("preserved")),
        "surface": share(None, field(top).astype(int)),
        "basement": share(None, field("basement_rock").astype(int))
                    if "basement_rock" in names else {},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("exports", nargs="+", type=Path)
    ap.add_argument("--glim", action="store_true",
                    help="also print grouped ratios against Earth")
    args = ap.parse_args()

    rows = {}
    for path in args.exports:
        root = path if (path / "manifest.json").is_file() else path.parent
        if not (root / "manifest.json").is_file():
            raise SystemExit(f"no manifest.json at or above {path}")
        rows[path.name or str(path)] = load(root)

    width = max(len(k) for k in rows) + 2
    print(f"{'':22}" + "".join(f"{k:>{width}}" for k in rows))
    for label, key, fmt in (("terrain hash", "hash", "s"),
                            ("land %", "land_pct", ".3f"),
                            ("mean land km", "mean_km", ".4f"),
                            ("max land km", "max_km", ".4f"),
                            ("preserved basins", "basins", "d")):
        cells = []
        for r in rows.values():
            v = r[key]
            cells.append(f"{str(v):>{width}} " if fmt == "s" or v is None
                         else f"{v:{width}{fmt}}")
        print(f"  {label:20}" + "".join(cells))

    for which in ("surface", "basement"):
        if not any(rows[k][which] for k in rows):
            continue
        codes = sorted({c for r in rows.values() for c, v in r[which].items()
                        if v > 0.005},
                       key=lambda c: -max(r[which].get(c, 0) for r in rows.values()))
        print(f"\n{which.upper()} composition, % of land")
        print(f"{'':22}" + "".join(f"{k:>{width}}" for k in rows))
        for code in codes:
            print(f"  {code:20}"
                  + "".join(f"{r[which].get(code, 0):{width}.3f}"
                            for r in rows.values()))

    if args.glim:
        print(f"\nSURFACE grouped, ratio to GLiM Earth")
        print(f"{'':22}" + "".join(f"{k:>{width}}" for k in rows) + "   Earth %")
        for group, codes in GROUPS.items():
            vals = [sum(r["surface"].get(c, 0) for c in codes)
                    for r in rows.values()]
            print(f"  {group:20}"
                  + "".join(f"{v / GLIM[group]:{width}.2f}" for v in vals)
                  + f"   {GLIM[group]:.1f}")


if __name__ == "__main__":
    main()
