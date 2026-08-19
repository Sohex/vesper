"""How much of this world's silicate weathering can reach the ocean.

The carbonate-silicate thermostat does not run on weathering. It runs on
weathering whose alkalinity is *delivered*: a closed basin weathers its
catchment, precipitates the carbonate on its own floor, and returns nothing to
the ocean, so its contribution never joins the feedback that stabilises CO2.
On a world where a large share of the land drains inward, the thermostat is
weaker than the total weathering rate suggests.

    thermostat efficiency = weathering over exorheic land / weathering over land

This measures it on the configured build, and separately on any other build in
`source/`, holding the climate fixed so the difference is basin geometry alone.
It also reports what a pending carve verdict would do, by reading the retain
fraction per basin from a carve list and treating a basin as delivering the
fraction of its weathering that its rim no longer holds back.

    python pedology/scripts/thermostat_efficiency.py
    python pedology/scripts/thermostat_efficiency.py --compare-builds \\
        --carve-list hydrography/data/carved-zoned-v4/carve_list.json

Two things to keep straight. The efficiency is *weathering*-weighted, not
area-weighted, and the two differ: endorheic land here is drier than exorheic
land and weathers about 20% less per unit area, so it decouples less alkalinity
than its area share implies. And a carve list is a request, not a terrain. Until
Orogen returns the carved export, the post-carve number is a projection onto the
current terrain's weathering field, not a measurement of a world that exists.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import (ANALYSIS, CONFIG, PEDOGENESIS, PROJECT_ROOT, SOURCE,
                    climatology_path)

import climatology  # noqa: E402  from lib/, put on sys.path by _paths
from paths import rel  # noqa: E402

from build_soil import EARTH_YEAR_DAYS, KELVIN, weathering_intensity


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_weathering(climatology: Path, params: dict) -> np.ndarray:
    """The same WHAK intensity `build_soil.py` computes, from the same fields."""
    with nc.Dataset(climatology) as data:
        temperature = climatology.annual_mean_of(data, "tas") - KELVIN
        scale = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        evaporation = -climatology.annual_mean_of(data, "evap")
        runoff = np.maximum(
            (climatology.annual_mean_of(data, "pr") - evaporation)
            * scale, 0.0)
    return weathering_intensity(runoff, temperature, params)


def read_terrain(build: str) -> dict:
    path = SOURCE / build / "exoplasim-T42" / "planet.nc"
    with nc.Dataset(path) as data:
        return {
            "path": path,
            "endorheic": np.asarray(data["is_endorheic"][:], dtype=float) > 0.5,
            "basin_index": np.asarray(data["basin_index"][:]).astype(int),
            "area": np.asarray(data["grid_cell_area"][:], dtype=float),
            # surface_class is the authoritative land definition; land_mask
            # floods the dry closed-basin floors, which is exactly the terrain
            # this measurement is about.
            "land": np.asarray(data["surface_class"][:]).astype(int) == 1,
        }


def efficiency(weathering: np.ndarray, terrain: dict,
               endorheic_fraction: np.ndarray | None = None) -> dict:
    """Weathering-weighted exorheic share, plus the area share for contrast."""
    land, area = terrain["land"], terrain["area"]
    closed = (terrain["endorheic"].astype(float) if endorheic_fraction is None
              else endorheic_fraction)
    weighted = weathering * area
    decoupled = (weighted * closed)[land].sum() / weighted[land].sum()
    endo = terrain["endorheic"] & land
    exo = (~terrain["endorheic"]) & land
    return {
        "endorheic_area_fraction": float((area * closed)[land].sum()
                                         / area[land].sum()),
        "decoupled_weathering_fraction": float(decoupled),
        "thermostat_efficiency": float(1.0 - decoupled),
        "mean_weathering_endorheic": float(
            np.average(weathering[endo], weights=area[endo])) if endo.any() else None,
        "mean_weathering_exorheic": float(
            np.average(weathering[exo], weights=area[exo])),
    }


def pending_carve(terrain: dict, build: str, carve_list: Path) -> tuple[np.ndarray, dict]:
    """Endorheic fraction per cell after applying a carve list's retain values.

    `retain` is 1 for a preserved basin, 0 for a carved one and in between for a
    marginal one whose rim is notched. Reading it as the fraction of the basin's
    weathering that stays decoupled treats a notched rim as partially delivering,
    which is the same linear reading `export_carve_list.py` puts on it.
    """
    carve_list = carve_list.resolve()
    verdict = json.loads(carve_list.read_text())
    by_id = {basin["id"]: basin for basin in verdict["basins"]}
    manifest = json.loads((SOURCE / build / "exoplasim-T42"
                           / "manifest.json").read_text())
    preserved = manifest["basins"]["preserved"]

    retain = np.ones(len(preserved))
    matched = 0
    for index, entry in enumerate(preserved):
        basin = by_id.get(entry.get("id"))
        if basin is not None:
            retain[index] = float(basin["retain"])
            matched += 1

    fraction = np.zeros(terrain["endorheic"].shape, dtype=float)
    inside = terrain["endorheic"] & (terrain["basin_index"] >= 0)
    fraction[inside] = retain[terrain["basin_index"][inside]]
    # Endorheic cells outside the catalogue are the smaller enclosed depressions
    # the basin thresholds never selected. No verdict carves them, so they stay
    # closed. See CLAUDE.md on why the naive union is wrong.
    fraction[terrain["endorheic"] & (terrain["basin_index"] < 0)] = 1.0
    provenance = {
        "carve_list": rel(carve_list),
        "carve_list_sha256": sha256(carve_list),
        "terrain_hash": verdict.get("terrain_hash"),
        "counts": verdict.get("counts"),
        "catalogue_basins": len(preserved),
        "basins_matched": matched,
    }
    return fraction, provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None)
    parser.add_argument("--carve-list", type=Path, default=None,
                        help="a carve_list.json whose verdict has been sent to "
                             "Orogen but whose terrain has not come back")
    parser.add_argument("--compare-builds", action="store_true",
                        help="measure every build in source/ against the same "
                             "climate, so the difference is basin geometry")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    pedo = yaml.safe_load(PEDOGENESIS.read_text())
    build = config["source_build"]
    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    from provenance import require_build
    require_build(climatology, "climatology")

    weathering = read_weathering(climatology, pedo["weathering"])
    terrain = read_terrain(build)
    current = efficiency(weathering, terrain)

    print(f"climatology  {rel(climatology)}")
    print(f"build        {build}\n")
    print(f"{'terrain':<32} {'endorheic area':>15s} {'decoupled W':>13s} "
          f"{'efficiency':>11s}")

    def show(name: str, result: dict) -> None:
        print(f"{name:32s} {result['endorheic_area_fraction']*100:14.1f}% "
              f"{result['decoupled_weathering_fraction']*100:12.1f}% "
              f"{result['thermostat_efficiency']:11.3f}")

    show(f"{build}, as exported", current)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatology": rel(climatology),
        "climatology_sha256": sha256(climatology),
        "source_build": build,
        "weathering": pedo["weathering"],
        "current": current,
        "definition": (
            "thermostat efficiency = weathering integrated over exorheic land "
            "divided by weathering over all land. Endorheic land weathers; its "
            "alkalinity never reaches the ocean, so it is decoupled from the "
            "carbonate-silicate feedback."),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }

    if args.carve_list:
        fraction, provenance = pending_carve(terrain, build, args.carve_list)
        projected = efficiency(weathering, terrain, fraction)
        show("after the pending carve", projected)
        report["pending_carve"] = {
            **provenance, **projected,
            "caveat": ("A projection onto the current terrain's weathering "
                       "field, not a measurement. The carved export does not "
                       "exist yet."),
        }

    if args.compare_builds:
        print()
        others = {}
        for other in sorted(p.name for p in SOURCE.iterdir() if p.is_dir()):
            path = SOURCE / other / "exoplasim-T42" / "planet.nc"
            if not path.is_file() or other == build:
                continue
            result = efficiency(weathering, read_terrain(other))
            others[other] = result
            show(f"{other}, same climate", result)
        report["other_builds"] = others
        report["comparison_note"] = (
            "Every build measured against one climatology, so the spread is "
            "basin geometry and not climate. Builds carved under superseded "
            "verdicts are not a trend: read them as separate hypotheses about "
            "which basins overflow.")

    endo = current["mean_weathering_endorheic"]
    exo = current["mean_weathering_exorheic"]
    print(f"\nmean weathering intensity   endorheic {endo:.3f}   exorheic {exo:.3f}")
    print("Endorheic land is drier, so it decouples less alkalinity than its "
          "area share implies.")

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    path = args.output or ANALYSIS / "thermostat_efficiency.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {rel(path)}")


if __name__ == "__main__":
    main()
