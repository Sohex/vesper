"""Score an LPJ-GUESS run against the pre-registered productivity prediction.

`notes/productivity-prediction.md` fixed ten numbered predictions, each with a
hit band and a clear-miss band, before the model had run a single Vesper
gridcell. This computes what the run actually did and reports hit or miss per
line. It does not adjust anything: a miss is the result.

    python biosphere/scripts/score_prediction.py biosphere/runs/<run_id>

The bands below are transcribed from that document and are not to be edited to
fit a result. If a band turns out to be badly posed, say so in the note and
re-register; do not quietly widen it here.

## Productivity is quoted over the nitrogen bracket, never from one arm

The Cleveland fixation slope `nfix_a` is a DECLARED BRACKET and not a value.
LPJ-GUESS states the range in `global.ins`, this project has no Vesper
observation of biological nitrogen fixation to narrow it with, and fixation
supplies an order more nitrogen than the declared deposition does. So the four
productivity lines are properties of a range, and a single run at the midpoint
does not measure them however precise it is.

    python biosphere/scripts/score_prediction.py <central run> \
        --nfix-arm <run at 0.102> --nfix-arm <run at 0.367>

Each arm is a run differing from the central one in `nfix_a` alone; the script
refuses an arm whose manifest says otherwise, and refuses an arm that does not
cover the cells being scored. With the ends supplied, lines 1 to 4 report the
span across them and a line is a HIT only when the WHOLE span lands in the hit
band. Without them the lines are reported as NOT QUOTED and are not scored:
collapsing an unresolved range onto its midpoint at the point of quotation is
the failure this refuses, and printing the midpoint with a verdict beside it is
exactly that collapse.

Lines 5 to 10 are structural and are scored from the central arm, which is what
they are registered against.
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

from _paths import CONFIG, PROJECT_ROOT, climatology_path
from gridding import gaussian_area_weights  # noqa: E402
from paths import rel  # noqa: E402

import climatology as climatology_lib
import orbit
from rootable import read_rootable
from lpj_output import reduce_table, require_lpj_acceptance
# ONE SOURCE for the declared fixation bracket: the script that writes it into
# the instruction file. Restating the ends here would be a second declaration of
# one quantity, and the two would drift the moment either moved.
from run_lpj_guess import NFIX_A_BRACKET
import lpj_pfts
import nc_geometry

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = COMPONENT_ROOT / "analysis"

# Earth reference, pinned in the prediction so the comparison cannot drift.
EARTH_LAND_AREA_KM2 = 149.0e6
EARTH_TOTAL_NPP_PGC = (55.0, 60.0)
EARTH_LAND_MEAN_NPP = (370.0, 400.0)

# Read from the file LPJ-GUESS itself parses, never restated here.
# Predictions 5, 8 and 9 are all tree-against-grass, so a missed type would
# move a scored line silently. lib/lpj_pfts.py carries the argument.
# Prediction 10 is "boreal NEEDLELEAF", which is the conjunction of two
# groups the PFT file declares and not a list of three codes: IBS is
# boreal too, and broadleaved, so an enumeration has to get the
# intersection right by hand every time a type is added.
# lib/lpj_pfts.py resolves it from the file the model reads.

# Cover above which a PFT counts as present in a cell. The registration said
# "tree FPC > 0.1" for prediction 5 and left 9 and 10 looser; the same threshold
# is applied to all three and stated here rather than chosen per line.
PRESENCE_FPC = 0.1


def read_table(path: Path, peers: list[Path] | None = None):
    """BIO-12 equilibrium mean per cell, plus its uncertainty report."""
    reduced = reduce_table(path, peers or ())
    return reduced.names, reduced.values, reduced.report


def band(value: float, hit: tuple[float, float],
         miss: tuple[float, float]) -> str:
    if hit[0] <= value <= hit[1]:
        return "HIT"
    if miss[0] <= value <= miss[1]:
        return "outside hit band, inside clear-miss band"
    return "CLEAR MISS"


def span_band(span: tuple[float, float], hit: tuple[float, float],
              miss: tuple[float, float]) -> str:
    """The verdict on a whole bracket, which is not the verdict on its middle.

    A prediction is hit when the range the model can produce lies inside the
    band, not when one arm of it does. A span straddling an edge is named as
    straddling it: that is the honest answer and it is a different result from
    either end's.
    """
    low, high = min(span), max(span)
    if hit[0] <= low and high <= hit[1]:
        return "HIT"
    if miss[0] <= low and high <= miss[1]:
        return "outside hit band, inside clear-miss band"
    if high < miss[0] or low > miss[1]:
        return "CLEAR MISS"
    return "STRADDLES the clear-miss edge"


def nfix_a_of(run: Path) -> float:
    """The fixation slope a run was made at, from its own manifest."""
    manifest = run / "run_manifest.json"
    if not manifest.is_file():
        raise SystemExit(f"{run} has no run_manifest.json, so the nfix_a it was "
                         "run at cannot be read and it cannot be an arm")
    physical = json.loads(manifest.read_text()).get("physical") or {}
    if "nfix_a" not in physical:
        raise SystemExit(f"{manifest} records no physical.nfix_a")
    return float(physical["nfix_a"])


def require_single_factor(central: Path, arm: Path) -> None:
    """Refuse an arm that differs from the central run in more than nfix_a.

    A sweep over one factor is only a sweep over that factor while everything
    else is held. Two runs that differ in their patch count as well as their
    fixation slope measure the sum of the two, and the bracket reported from
    them would be wider or narrower than the nitrogen bracket by an amount
    nothing records.
    """
    held = ("nyear", "npatch", "root_seed", "ntransform_profile", "ifbvoc",
            "run_peatland", "ifmethane", "nyear_spinup")
    a = json.loads((central / "run_manifest.json").read_text())
    b = json.loads((arm / "run_manifest.json").read_text())
    differ = [key for key in held
              if (a.get("physical") or {}).get(key) != (b.get("physical") or {}).get(key)]
    if a.get("source_build") != b.get("source_build"):
        differ.append("source_build")
    for name in ("driver", "soilmap", "pfts", "binary", "vesper_h"):
        was = ((a.get("inputs") or {}).get(name) or {}).get("sha256")
        now = ((b.get("inputs") or {}).get(name) or {}).get("sha256")
        if was != now:
            differ.append(name)
    if differ:
        raise SystemExit(
            f"{arm.name} differs from {central.name} in {', '.join(differ)} as "
            "well as nfix_a, so the span between them is not the nitrogen "
            "bracket. Re-run the arm changing --nfix-a alone.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="an LPJ-GUESS run directory")
    parser.add_argument("--climatology", type=Path, default=None,
                        help="the climatology the run was driven from, for the "
                             "forcing-derived lines and cell areas")
    parser.add_argument("--rootable", type=Path, default=None,
                        help="BIO-11 rootable-fraction artifact; default the "
                             "active build and rung")
    parser.add_argument("--equilibrium-peer", type=Path, action="append", default=[],
                        help="comparable LPJ run directory with another root seed "
                             "or patch count; repeat to measure stochastic spread")
    parser.add_argument("--nfix-arm", type=Path, action="append", default=[],
                        help=f"an LPJ run differing from the scored one in nfix_a "
                             f"alone. The declared bracket is "
                             f"{NFIX_A_BRACKET[0]} to {NFIX_A_BRACKET[1]} and both "
                             "ends are required before the productivity lines are "
                             "quoted at all; repeat the flag")
    args = parser.parse_args()

    run = args.run
    require_lpj_acceptance(run)
    manifest_path = run / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    config = yaml.safe_load(CONFIG.read_text())
    climatology = args.climatology or climatology_path()

    # Simulation years are this world's years. Everything registered is per Earth
    # year. Getting this wrong looks like missing every line low by 2.
    to_earth = manifest.get("annual_flux_per_orbit_to_per_earth_year")
    if to_earth is None:
        to_earth = 1.0 / orbit.earth_years_per_model_year(config)

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        # The bins are NOT equal length: pyburn splits the raw stream with an
        # integer linspace, so at 182 records in twelve bins two hold sixteen
        # records and the rest fifteen. A bare `.mean(axis=0)` asserts they are
        # equal; `annual_mean_of` measures them off the file's own bin centres.
        temperature = climatology_lib.annual_mean_of(data, "tas") - 273.15
        precip = (climatology_lib.annual_mean_of(data, "pr")
                  * 1000.0 * 86400.0 * orbit.EARTH_CALENDAR_YEAR_DAYS)

    radius_km = nc_geometry.planet_radius_m(config) / 1000.0
    cell_km2 = (gaussian_area_weights(lat, len(lon), what=str(climatology))
                * 4.0 * np.pi * radius_km ** 2)
    rootable, rootable_provenance = read_rootable(
        config, lat, lon, args.rootable, land=land)
    effective_area_km2 = cell_km2 * rootable
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)

    npp_names, npp, npp_window = read_table(
        run / "anpp.out", [peer / "anpp.out" for peer in args.equilibrium_peer])
    lai_names, lai, lai_window = read_table(
        run / "lai.out", [peer / "lai.out" for peer in args.equilibrium_peer])
    fpc_names, fpc, fpc_window = read_table(
        run / "fpc.out", [peer / "fpc.out" for peer in args.equilibrium_peer])

    # Map every simulated cell onto the grid so areas can weight it.
    simulated = np.zeros(land.shape, dtype=bool)
    npp_grid = np.zeros(land.shape)
    lai_grid = np.zeros(land.shape)
    tree_grid = np.zeros(land.shape)
    grass_grid = np.zeros(land.shape)
    c4_grid = np.zeros(land.shape)
    boreal_grid = np.zeros(land.shape)

    grass_pfts = lpj_pfts.grass()
    tree_idx = [i for i, n in enumerate(fpc_names)
                if n not in grass_pfts and n != "Total"]
    grass_idx = [i for i, n in enumerate(fpc_names) if n in grass_pfts]
    c4_idx = [i for i, n in enumerate(fpc_names) if n == "C4G"]
    boreal_pfts = lpj_pfts.members("boreal", "needleleaved")
    boreal_idx = [i for i, n in enumerate(fpc_names) if n in boreal_pfts]
    npp_total = npp_names.index("Total")
    lai_total = lai_names.index("Total")

    for j in range(len(lat)):
        for i in range(len(lon)):
            key = (round(float(lon_signed[i]), 2), round(float(lat[j]), 2))
            row = fpc.get(key)
            if row is None:
                continue
            simulated[j, i] = True
            npp_grid[j, i] = npp[key][npp_total]
            lai_grid[j, i] = lai[key][lai_total]
            tree_grid[j, i] = sum(row[k] for k in tree_idx)
            grass_grid[j, i] = sum(row[k] for k in grass_idx)
            c4_grid[j, i] = sum(row[k] for k in c4_idx)
            boreal_grid[j, i] = sum(row[k] for k in boreal_idx)

    covered = simulated & land & (rootable > 0.0)
    rootable_cells = land & (rootable > 0.0)
    n_land = int(rootable_cells.sum())
    n_done = int(covered.sum())
    if n_done == 0:
        raise SystemExit("no simulated cells matched the climatology grid")
    partial = n_done < n_land

    area = effective_area_km2[covered]
    land_area_km2 = float(cell_km2[land].sum())
    rootable_area_km2 = float(effective_area_km2[land].sum())

    def wmean(field: np.ndarray) -> float:
        return float(np.average(field[covered], weights=area))

    def wfrac(mask: np.ndarray) -> float:
        return float(np.average(mask[covered].astype(float), weights=area))

    # gC/m2 per Earth year, from kgC/m2 per simulation year.
    land_mean_npp = wmean(npp_grid) * 1000.0 * to_earth
    total_npp_pgc = land_mean_npp * rootable_area_km2 * 1e6 / 1e15
    earth_land_mean = sum(EARTH_LAND_MEAN_NPP) / 2.0
    earth_total = sum(EARTH_TOTAL_NPP_PGC) / 2.0

    # -- THE NITROGEN BRACKET, carried rather than collapsed --------------------
    #
    # The same reduction on the same cells with the same weights, once per arm,
    # so what separates the numbers is nfix_a and nothing else. The central run
    # is an arm like any other and is not privileged in the span.
    def npp_of(arm_run: Path) -> float:
        names, values, _ = read_table(arm_run / "anpp.out")
        total_column = names.index("Total")
        grid = np.zeros(land.shape)
        seen = np.zeros(land.shape, dtype=bool)
        for j in range(len(lat)):
            for i in range(len(lon)):
                row = values.get((round(float(lon_signed[i]), 2),
                                  round(float(lat[j]), 2)))
                if row is None:
                    continue
                grid[j, i] = row[total_column]
                seen[j, i] = True
        missing = int((covered & ~seen).sum())
        if missing:
            raise SystemExit(
                f"{arm_run.name} is missing {missing} of the {n_done} cells "
                "being scored, so its land mean is over a different land than "
                "the one it would be bracketing")
        return float(np.average(grid[covered], weights=area)) * 1000.0 * to_earth

    arms = {nfix_a_of(run): land_mean_npp}
    for arm_run in args.nfix_arm:
        require_lpj_acceptance(arm_run)
        require_single_factor(run, arm_run)
        value = nfix_a_of(arm_run)
        if value in arms:
            raise SystemExit(
                f"{arm_run.name} was run at nfix_a {value}, which {run.name} or "
                "another arm already covers")
        arms[value] = npp_of(arm_run)

    missing_ends = [end for end in NFIX_A_BRACKET
                    if not any(abs(end - value) < 1e-9 for value in arms)]
    bracket_carried = not missing_ends
    npp_span = (min(arms.values()), max(arms.values()))
    total_span = tuple(value * rootable_area_km2 * 1e6 / 1e15 for value in npp_span)

    # Prediction 7 is a property of the forcing, not of LPJ-GUESS: the
    # registration defined it through the Miami model, so it is recomputed the
    # same way rather than invented from model output.
    temp_limit = 3000 / (1 + np.exp(1.315 - 0.119 * temperature))
    prec_limit = 3000 * (1 - np.exp(-0.000664 * precip))
    water_limited = wfrac(prec_limit < temp_limit)
    nonrootable_area = float((cell_km2 * (1.0 - rootable))[land].sum())
    low_lai_rootable_area = float(
        effective_area_km2[covered & (lai_grid < 0.5)].sum())
    effective_barren_land_fraction = (
        (nonrootable_area + low_lai_rootable_area) / land_area_km2)

    productivity = [
        ("1  land-mean NPP relative to Earth",
         tuple(v / earth_land_mean for v in npp_span), (0.7, 1.4), (0.5, 2.0), "x"),
        ("2  land-mean NPP, gC/m2 per Earth year",
         npp_span, (280.0, 560.0), (200.0, 800.0), ""),
        ("3  total NPP relative to Earth",
         tuple(v / earth_total for v in total_span), (1.6, 2.8), (1.2, 3.5), "x"),
        ("4  total NPP, PgC per Earth year",
         total_span, (90.0, 175.0), (60.0, 240.0), ""),
    ]

    results = [
        ("5  tree cover, land fraction with tree FPC > 0.1",
         wfrac(tree_grid > PRESENCE_FPC), (0.40, 0.70), (0.25, 0.80), ""),
        ("6  effectively barren land, LAI < 0.5",
         effective_barren_land_fraction, (0.08, 0.25), (0.0, 0.40), ""),
        ("7  water-limited land (from the forcing, not the run)",
         water_limited, (0.50, 0.70), (0.35, 1.0), ""),
        ("8  grass cover exceeds tree cover",
         1.0 if wmean(grass_grid) > wmean(tree_grid) else 0.0,
         (1.0, 1.0), (1.0, 1.0), " (1 = yes)"),
        ("9  C4 grasses above 10% of land",
         wfrac(c4_grid > PRESENCE_FPC), (0.10, 1.0), (0.05, 1.0), ""),
        ("10 boreal needleleaf below 10% of land",
         wfrac(boreal_grid > PRESENCE_FPC), (0.0, 0.10), (0.0, 0.20), ""),
    ]

    print(f"run              {run.name}")
    print(f"cells scored     {n_done} of {n_land} rootable land cells"
          + ("   PARTIAL RUN, treat every line as indicative" if partial else ""))
    print(f"simulation year  {to_earth:.4f} Earth years")
    print(f"land area        {land_area_km2 / 1e6:.1f}e6 km2\n")
    print(f"rootable area    {rootable_area_km2 / 1e6:.1f}e6 km2\n")
    # The REPORTED span, which is the one the scored values were taken over and
    # the one the contract's drift bound certifies. The per-cell half's window is
    # a different and shorter span and is not what these numbers came from.
    print("equilibrium      "
          f"{fpc_window['reported']['complete_forcing_cycles']} complete forcing "
          f"cycles; seed uncertainty "
          f"{fpc_window['uncertainty']['seed']['status']}; patch uncertainty "
          f"{fpc_window['uncertainty']['patch_count']['status']}\n")
    arm_line = ", ".join(f"nfix_a {value:g}" for value in sorted(arms))
    print(f"nitrogen bracket {arm_line}"
          + ("" if bracket_carried else
             "   BRACKET NOT CARRIED, productivity is not quoted"))
    if not bracket_carried:
        print("  missing " + ", ".join(f"an arm at nfix_a {end:g}"
                                       for end in missing_ends)
              + ". Run each with --nfix-a and pass it as --nfix-arm.")
    print()

    scored = []
    for name, span, hit, miss, unit in productivity:
        low, high = min(span), max(span)
        shown = f"{low:.3f}" if low == high else f"{low:.3f}-{high:.3f}"
        verdict = (span_band(span, hit, miss) if bracket_carried else
                   "NOT QUOTED: the declared nfix_a bracket is not carried")
        print(f"  {name:52s} {shown:>17s}{unit:3s} {verdict}")
        scored.append({"prediction": name,
                       "span_over_arms_supplied": [low, high],
                       "hit_band": list(hit), "clear_miss_band": list(miss),
                       "verdict": verdict,
                       "quoted": bracket_carried})
    for name, value, hit, miss, unit in results:
        verdict = band(value, hit, miss)
        print(f"  {name:52s} {value:14.3f}{unit:6s} {verdict}")
        scored.append({"prediction": name, "value": value, "hit_band": list(hit),
                       "clear_miss_band": list(miss), "verdict": verdict,
                       "quoted": True})

    hits = sum(1 for r in scored if r["verdict"] == "HIT")
    quoted = sum(1 for r in scored if r["quoted"])
    print(f"\n  {hits} of {quoted} quoted predictions hit"
          + ("" if quoted == len(scored)
             else f"; {len(scored) - quoted} not quoted"))

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "run": str(run),
        "run_manifest": manifest,
        "climatology": rel(climatology),
        "cells_scored": n_done,
        "land_cells": n_land,
        "partial": partial,
        "annual_flux_per_orbit_to_per_earth_year": to_earth,
        "land_area_km2": land_area_km2,
        "rootable_area_km2": rootable_area_km2,
        "rootable_surface": rootable_provenance,
        "equilibrium_windows": {
            "anpp": npp_window, "lai": lai_window, "fpc": fpc_window},
        "prediction_area_basis": (
            "intensive vegetation quantities are weighted and normalised by "
            "BIO-11 rootable area; extensive totals multiply by that area. "
            "The effectively-barren diagnostic additionally counts the "
            "non-rootable share against total model land."),
        "nonrootable_area_km2": nonrootable_area,
        "land_mean_npp_gc_m2_earth_year": land_mean_npp,
        "total_npp_pgc_earth_year": total_npp_pgc,
        "nfix_a": {
            "declared_bracket": list(NFIX_A_BRACKET),
            "arms": {f"{value:g}": mean for value, mean in sorted(arms.items())},
            "bracket_carried": bracket_carried,
            "missing_ends": missing_ends,
            "land_mean_npp_span_gc_m2_earth_year": list(npp_span),
            "total_npp_span_pgc_earth_year": list(total_span),
            "note": ("productivity is a property of the declared fixation "
                     "bracket; a single arm is not a quotation of it"),
        },
        "earth_reference": {
            "land_area_km2": EARTH_LAND_AREA_KM2,
            "total_npp_pgc": list(EARTH_TOTAL_NPP_PGC),
            "land_mean_npp": list(EARTH_LAND_MEAN_NPP),
        },
        "predictions": scored,
        "hits": hits,
        "note": ("Bands are transcribed from notes/productivity-prediction.md and "
                 "are not to be edited to fit a result. A miss is the result. "
                 "Lines 1 to 4 are quoted over the declared nfix_a bracket and "
                 "are not scored until both ends are supplied."),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    out_dir = ANALYSIS / run.name
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "prediction_score.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {rel(path)}")


if __name__ == "__main__":
    main()
