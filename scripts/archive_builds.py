#!/usr/bin/env python3
"""Extract a superseded build's identity, then delete its payload.

    python scripts/archive_builds.py               # dry run
    python scripts/archive_builds.py --execute

Each Orogen export is about 2 GB of mesh, grids and a 66 MB manifest. Seven of
them is 13 GB, and six describe terrains this project has abandoned -- all the
pre-gravity-change builds were superseded wholesale when gravity moved to 12.81.

Deleting a payload is safe here in a way deleting a climate run's output is not,
because **the recipe survives**: Orogen regenerates terrain from the planet code
plus a carve list in one pass, and both are kept. A build is replaced wholesale
rather than edited, so an export is reproducible in a way a 60-orbit spin-up is
not.

What must survive is identity, not data. `lib/orogen.py` refuses a build it has
not been checked against, and results already computed from a superseded terrain
have to stay readable and datable. That needs the hashes and the headline
numbers, which is kilobytes, not the mesh.

The per-build hydrography products under `hydrography/data/<build>/` are NOT
touched: they are small, they are what downstream work actually cites, and they
are derived rather than regenerable in one pass.

Anything the registry in `lib/orogen.py` does not know is SKIPPED, as is any
build missing one of its five grid exports. `source_build` names the OLD build
for as long as it takes to generate a new one, so "not active" and "superseded"
are not the same claim, and only the registry can tell them apart.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source"
ARCHIVE = ROOT / "archive" / "builds"

sys.path.insert(0, str(ROOT / "lib"))
from orogen import _KNOWN_TERRAIN_HASHES  # noqa: E402

# A build is these five exports. Fewer means it is still being generated.
EXPECTED_GRIDS = ["exoplasim-T21", "exoplasim-T42", "exoplasim-T63",
                  "exoplasim-T85", "grid-512x256"]

# Everything needed to identify a build and to date a result computed from it.
IDENTITY_KEYS = ["hashes", "seed", "numRegions", "planet", "planetRadiusKm",
                 "landSeaMask", "code", "version", "generated"]


def stub(manifest: dict) -> dict:
    out = {k: manifest[k] for k in IDENTITY_KEYS if k in manifest}
    b = manifest.get("basins") or {}
    out["basins"] = {k: v for k, v in b.items()
                     if k in ("counts", "resolution", "selectionSource",
                              "thresholds")}
    lit = manifest.get("lithology") or {}
    out["lithology"] = {k: v for k, v in lit.items() if k != "perRegion"}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()

    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    active = cfg.get("source_build")

    freed = 0
    print(f"active build: {active}\n")
    for d in sorted(SOURCE.iterdir()):
        if not d.is_dir() or d.name == "maps":
            continue
        size = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
        if d.name == active:
            print(f"  LIVE     {d.name:26}{size/1e9:>6.2f} GB  untouched")
            continue
        man = d / "exoplasim-T42" / "manifest.json"
        if not man.is_file():
            print(f"  SKIP     {d.name:26}{size/1e9:>6.2f} GB  no T42 manifest")
            continue

        # Refuse anything the registry does not know. An unregistered build is
        # either mid-generation or unidentified, and deleting its payload is
        # wrong either way -- a build being written right now is not superseded,
        # it is the next active one. `source_build` still names the OLD build
        # while a new one is generated, so keying on that alone would archive the
        # build currently being produced.
        try:
            m = json.loads(man.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            print(f"  SKIP     {d.name:26}{size/1e9:>6.2f} GB  manifest unreadable "
                  f"(generating?)")
            continue
        if m.get("hashes", {}).get("finalElevation") not in _KNOWN_TERRAIN_HASHES:
            print(f"  SKIP     {d.name:26}{size/1e9:>6.2f} GB  NOT REGISTERED in "
                  f"lib/orogen.py -- register it or delete it by hand")
            continue
        missing = [g for g in EXPECTED_GRIDS if not (d / g / "manifest.json").is_file()]
        if missing:
            print(f"  SKIP     {d.name:26}{size/1e9:>6.2f} GB  incomplete, missing "
                  f"{', '.join(missing)}")
            continue

        freed += size
        print(f"  ARCHIVE  {d.name:26}{size/1e9:>6.2f} GB  -> identity stub")
        if args.execute:
            dest = ARCHIVE / d.name
            dest.mkdir(parents=True, exist_ok=True)
            (dest / "identity.json").write_text(
                json.dumps(stub(m), indent=2) + "\n", encoding="utf-8")
            for extra in ("README.txt", "SUPERSEDED.md"):
                for p in d.rglob(extra):
                    shutil.copy2(p, dest / f"{p.parent.name}_{extra}")
            shutil.rmtree(d)

    print(f"\n{'freed' if args.execute else 'would free'} {freed/1e9:.1f} GB")
    if not args.execute:
        print("dry run; pass --execute to do it")
        return
    (ARCHIVE / "README.md").write_text(
        "# Archived builds\n\n"
        "Identity of superseded World Orogen exports whose payload has been\n"
        "deleted. `identity.json` carries the hashes, seed, region count, planet\n"
        "parameters, land/sea mask summary and basin counts -- everything needed\n"
        "to recognise a build and to date a result computed from it.\n\n"
        "The payload is regenerable: Orogen rebuilds terrain from the planet code\n"
        "plus a carve list in one pass, and both are kept. That is why these can\n"
        "go and a climate run's output cannot.\n\n"
        "`lib/orogen.py` remains the registry and still refuses a build it has not\n"
        "been checked against. Per-build hydrography products under\n"
        "`hydrography/data/<build>/` were not touched.\n\n"
        f"Archived {datetime.now(timezone.utc).date().isoformat()}.\n",
        encoding="utf-8")
    print(f"wrote {ARCHIVE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
