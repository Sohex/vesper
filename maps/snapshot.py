#!/usr/bin/env python3
"""Take the map's picture where the pipeline currently stands.

    python maps/snapshot.py --step surface_water
    python maps/snapshot.py --step hydrography --all
    python maps/snapshot.py --step orogen --width 6000

This is the driver `scripts/pipeline.py --plan` names after every step the map
is reachable from. It runs `build_basemap.py` when the raster on disk is not the
world the config now describes, then `render_projections.py`, and the frame lands
in `maps/data/<build>/<frame>/` with its row in `INDEX.json`.

**The point is a series.** One frame per state of the world, so a pass can be
watched changing -- the coastline the carve moved, the lakes the water balance
filled, the biomes the new climatology repainted -- rather than only its end
state being drawn. That only works if the frames are taken as the pass runs,
which is why the plan carries a line for each of them.

**Running it at every step is cheap because a frame is keyed on its inputs.**
The map is drawn from four things: the terrain, the classification, the
climatology and the lake solution. A step that moves none of them leaves the
picture identical, and this records that step against the frame that already
exists instead of rendering a second copy of it. `frames.py` carries the
argument. `--force` renders anyway, which is what to reach for after a change to
the RENDERING rather than to the world.

It generates no artifact of its own: `basemap` and `projections` are the steps,
and this decides which of them still has work to do.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import frames  # noqa: E402

BUILD = HERE / "build"
ROOT = HERE.parent


def basemap_is_current(climatology: Path | None,
                       basemap_width: int | None = None) -> tuple[bool, str]:
    """Does the raster on disk show the world the config now describes, at the
    width being asked for?

    Compares the base map's recorded fingerprint against one taken over the
    inputs as they resolve NOW. Cheap -- four hashes -- and it is the check that
    decides whether the 25-second rebuild is needed, so it has to be taken over
    the same identities the frame is keyed on rather than over file timestamps.
    A newer climatology written under the same name is exactly the case an mtime
    comparison gets wrong in the direction that matters.

    THE WIDTH IS PART OF THE QUESTION. `build_basemap` records its resolution in
    the provenance but outside the fingerprint, which is right -- the width does
    not change which world the picture shows -- yet the caller only forwards
    `--basemap-width` when a rebuild is already triggered. So asking for a
    different width against an unmoved fingerprint kept the old raster, and the
    option ran, reported and changed nothing. world-60x0.
    """
    prov = BUILD / "basemap_provenance.json"
    if not prov.is_file():
        return False, "no base map has been built"
    recorded = json.loads(prov.read_text(encoding="utf-8"))
    if "inputs" not in recorded:
        return False, "the base map predates the frame manifest"

    import build_basemap
    if climatology is not None:
        build_basemap._CLIMATOLOGY_OVERRIDE = climatology.resolve()
    src = build_basemap.builds.mesh_export()
    manifest = json.loads((src / "manifest.json").read_text(encoding="utf-8"))
    water = build_basemap.surface_water_path()
    now = frames.inputs(
        source_build=src.parent.name,
        terrain_hash=manifest["hashes"]["finalElevation"],
        climatology=build_basemap._climatology(),
        classification=build_basemap._classification(),
        surface_water=(water if water.exists() else None),
    )
    if frames.fingerprint(now) == recorded["fingerprint"]:
        if basemap_width is not None:
            built = (recorded.get("resolution") or [None])[0]
            if built != basemap_width:
                return False, (f"the base map on disk is {built} px wide and "
                               f"{basemap_width} was asked for")
        return True, "the base map is the world the config describes"
    moved = [k for k, v in now.items() if recorded["inputs"].get(k) != v
             and not k.endswith("_sha256")]
    moved += [k[:-7] for k, v in now.items() if k.endswith("_sha256")
              and recorded["inputs"].get(k) != v]
    return False, "moved since the base map was built: " + ", ".join(sorted(set(moved)))


def run(script: str, *argv: str) -> None:
    cmd = [sys.executable, str(HERE / script), *argv]
    print("+ " + " ".join(cmd[1:]))
    subprocess.run(cmd, check=True, cwd=str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", required=True,
                    help="the pipeline step just run; checked against the graph")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="REGULAR climatology to tint from, passed through to "
                         "build_basemap.py")
    ap.add_argument("--width", type=int, default=None,
                    help="render width for the projections")
    ap.add_argument("--basemap-width", type=int, default=None,
                    help="equirectangular width of the base map itself")
    ap.add_argument("--all", action="store_true",
                    help="every projection rather than the authagraph alone")
    ap.add_argument("--force", action="store_true",
                    help="rebuild and re-render even where nothing moved")
    args = ap.parse_args()
    frames.check_step(args.step)

    current, why = basemap_is_current(args.climatology, args.basemap_width)
    print(why)
    if args.force or not current:
        base_args = []
        if args.climatology is not None:
            base_args += ["--climatology", str(args.climatology)]
        if args.basemap_width is not None:
            base_args += ["--width", str(args.basemap_width)]
        run("build_basemap.py", *base_args)

    render_args = ["--step", args.step]
    if args.width is not None:
        render_args += ["--width", str(args.width)]
    if args.all:
        render_args.append("--all")
    if args.force:
        render_args.append("--force")
    run("render_projections.py", *render_args)


if __name__ == "__main__":
    main()
