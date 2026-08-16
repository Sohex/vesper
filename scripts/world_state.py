#!/usr/bin/env python3
"""Generate world_state.json: what this project currently knows about Vesper.

Generated, not hand-maintained. Every prose summary in this project has drifted
from its artifacts at least once, usually within a commit of the number changing,
so the state file reads the artifacts instead of restating them. Anything here
that cannot be derived is in `curated`, which lives in this file under version
control rather than being edited into the JSON.

    python scripts/world_state.py

Re-run it after anything that changes a build, a run, or a verdict.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 1

# --- things no artifact can tell us -----------------------------------------
# Judgement, not data. Kept in the generator so it is reviewed in diffs.
CURATED = {
    "open_questions": [
        "Does the climate response follow the albedo arithmetic? Land-mean albedo "
        "on carved-zoned is 0.1722 against the baseline's 0.2240, predicting about "
        "299.4 K. The chain is four linear steps across a response we have "
        "repeatedly declined to extrapolate.",
        "The habitable flux depends on the biosphere, which is assumed rather than "
        "modelled. Bare rock and vegetated reach 290-293 K at non-overlapping "
        "fluxes, so LPJ-GUESS determines not just the climate but where the planet "
        "should be placed.",
        "Iteration 1 systematically over-carved: carving removes bright evaporite, "
        "the world warms, and marginal basins would have stayed closed. The "
        "overshoot should be measured on iteration 2 by re-evaluating the carved "
        "set against the new climate.",
        "Roughly 228 basins with zero catchment runoff were preserved but have "
        "P > E over open water, so a lake there grows and spills. They should go "
        "through the marginal band with a retain fraction, not be flipped to 0.",
    ],
    "known_uncertainties": [
        {"item": "stellar spectrum", "value": "k25v", "range": ["k2", "k25v"],
         "worth_w_m2": 0.63, "worth_k": 0.59,
         "note": "Resolved going forward, inherited backwards. ExoPlaSim's k2.dat "
                 "is the star K2-18, an M2.5V at about 3450 K, not the spectral "
                 "type K2. Every run of the first three eras used it, so their "
                 "snow and ice are 0.10-0.17 too dark, worth 0.4-0.7 K on the warm "
                 "baseline and more wherever there is ice. k25v is now selected; "
                 "results predating it remain valid in kind with the direction of "
                 "the error known. See exoplasim/notes/stellar-spectrum-audit.md."},
        {"item": "mrro is river-routed, not local runoff", "value": None,
         "note": "Not a defect, a misreading, and recorded here because this "
                 "file previously called it one. landmod.f90's roffstep calls "
                 "mkradv, which advects runoff downhill and modifies its argument "
                 "in place, so the output is local generation minus river outflow "
                 "plus inflow. That explains the negative values, the nonzero "
                 "ocean values, the land shortfall correlating 0.60 with "
                 "precipitation, and the 99.1% global conservation. The CF name "
                 "surface_runoff is what misleads. Anything needing water that "
                 "drains through the local profile should use P - E; pedology "
                 "does. See exoplasim/notes/water-and-energy-closure.md."},
        {"item": "land runoff, and therefore the whole soil", "value": 0.192,
         "range": [0.192, 0.35],
         "note": "Much less under-determined than previously recorded here, once "
                 "the mrro defect above is set aside. The land runoff ratio is "
                 "18.8% against Earth's ~35%, not 2.8%: drier than Earth but not "
                 "extraordinarily so. Weathering intensity is 0.50 rather than "
                 "0.19, and the runoff-versus-precipitation bracket narrows from "
                 "14.6x to 5.6x. Land-mean clay 0.33, water capacity 263 mm. The "
                 "remaining spread is whether weathering follows drainage or "
                 "rainfall, which is a modelling choice, not a defect."},
        {"item": "nitrogen fixation, nfix_a", "value": 0.234,
         "range": [0.102, 0.367], "worth_npp_fraction": 0.182,
         "note": "The real nitrogen lever, and LPJ-GUESS documents the range "
                 "itself in global.ins as conservative/central/upper. Measured "
                 "over 18 cells: NPP 0.2356, 0.2642, 0.2785 kgC/m2 across the "
                 "three, a monotone 18.2% span. It is the Cleveland "
                 "fixation-versus-evapotranspiration fit, an Earth calibration "
                 "with no way to check it here. It does scale correctly with this "
                 "world's calendar, being a flux per unit AET."},
        {"item": "nitrogen deposition", "value": 0.5, "range": [0.1, 15.0],
         "note": "Declared with no basis, and measured to be second-order, which "
                 "is the opposite of what was recorded here before. At 0.5 "
                 "kgN/ha/yr deposition supplies 0.50 against 4.36 from biological "
                 "fixation and 25.4 from mineralisation, so 1.7% of the nitrogen "
                 "plants actually receive. Sweeping it over 150x, 0.1 to 15, moves "
                 "NPP by 12% non-monotonically, which at 18 cells and npatch 5 is "
                 "patch stochasticity. Recycling and fixation dominate; the "
                 "assumption stands but carries little."},
        {"item": "playa_clastic albedo", "value": 0.30, "range": [0.25, 0.33],
         "worth_w_m2": 1.99, "worth_k": 1.87,
         "note": "Dominant remaining albedo lever: 18.8% of land against the salt "
                 "crust's 2.1%, 2.48 W/m2 per 0.10. Itself a mixture of clay playa "
                 "and desert-varnished fan gravel. Recorded, not tuned."},
        {"item": "salt crust extent", "value": 0.106, "range": [0.10, 0.40],
         "note": "Degenerate with crust albedo in the effective value: 10.6% at "
                 "0.50 and 20% at 0.40 are the same planet from the climate side. "
                 "Only an independent measure of flooding frequency separates "
                 "them, and the lake solver is not validated enough to provide it."},
        {"item": "lake levels", "value": None,
         "note": "The equilibrium solver has never been validated against anything, "
                 "unlike Penman which was checked against the model's own ocean "
                 "cells to 1.7%. Least-checked numbers in the pipeline."},
        {"item": "energy budget closure", "value": 0.446,
         "note": "Net TOA runs 0.27-0.45 W/m2 more negative than the net surface "
                 "flux, which cannot both be true in steady state, and the "
                 "|mean TOA| < 0.5 criterion is applied inside a bias of its own "
                 "size. Better characterised now: a naive rss+rls+hfss+hfls sum "
                 "misses the 0.294 W/m2 consumed melting snow, which hfns books "
                 "correctly and which matches global snowmelt's latent heat of "
                 "fusion to three digits. So the surface terms are internally "
                 "consistent and the residual sits between the surface and the "
                 "top of the atmosphere. It is structural, not climatic: the gap "
                 "is -0.4529, -0.4656 and -0.4455 across three runs spanning "
                 "different fluxes, terrains and mean temperatures, a 2% spread. "
                 "That rules out sea ice, snow and fusion, all of which scale "
                 "with climate. It is not the spectral dissipation leak either, "
                 "since mkdheat returns friction and biharmonic diffusion heating "
                 "under ndheat, on by default. An atmosphere truly losing this "
                 "would cool 1.4 K per Earth year and these do not, so a heating "
                 "term is missing from the diagnostic sum or one diagnostic is "
                 "offset. PlaSim ships a 28-term decomposition on output codes "
                 "360-387 under nenergy, which would name it in a short run. "
                 "Note rainmod.f90:524 inverted als and alv in that diagnostic, now "
                 "fixed. Seasonal resolution reframes it further: per time bin "
                 "the gap swings about 14 W/m2 peak to peak, repeating closely "
                 "between orbits, so the annual residual is 3% of a large "
                 "seasonal storage term failing to cancel rather than a uniform "
                 "leak. See exoplasim/notes/water-and-energy-closure.md."},
        {"item": "retain fraction", "value": None,
         "note": "Measures distance from threshold, not incision capacity. A basin "
                 "that barely trickles over its sill is cut the same as one that "
                 "pours. Should become a function of overflow discharge."},
    ],
    "next_steps": [
        "Finish the carved-zoned baseline at 0.96 and compare against the 299.4 K "
        "prediction.",
        "Choose the flux from what it measures, not from the estimate.",
        "Supply equilibrium lake levels as --water-level so basin floors stop being "
        "exported as dry ground.",
        "Run LPJ-GUESS on the resulting climatology and replace the assumed "
        "biosphere, then score it against biosphere/notes/productivity-prediction.md.",
        "Carry the nfix_a bracket, 0.102-0.367, through to the productivity "
        "numbers rather than quoting the central value alone.",
        "Run a short segment with nenergy=1 to name the constant -0.455 W/m2 "
        "offset between the top-of-atmosphere and surface budgets, fixing the "
        "inverted als/alv in rainmod.f90:524 first.",
        "Enable model.soil_water_source: pedology and check the soil-runoff loop "
        "converges rather than oscillating; expect a weak response.",
        "Iteration-2 carve verdict, with retain as a function of discharge and the "
        "228 zero-runoff spillers handled.",
        "One T85 equilibrium once terrain and biosphere settle, for regional "
        "products.",
        "Stellar cycle last, 0.91-1.01 S-Earth over 8 Earth years.",
    ],
}


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def builds(active: str) -> dict:
    """The active build only.

    `world_state.json` is what the project currently knows, and a superseded
    build is not that -- it is history. Enumerating all seven put six wrong
    terrains next to the right one under identical keys, which is how a reader
    (or a script) picks the wrong basin count.

    The registry is a separate thing and stays complete: `lib/orogen.py` keeps
    every build it has been checked against, so results already computed from a
    superseded terrain stay readable and datable. Identity belongs in the
    registry; state belongs here.
    """
    out = {}
    for d in sorted((ROOT / "source").iterdir()):
        if d.name != active:
            continue
        man = d / "exoplasim-T42" / "manifest.json"
        if not man.is_file():
            continue
        m = json.loads(man.read_text(encoding="utf-8"))
        b, lit = m["basins"], m["lithology"]
        comp = {x["code"]: x["fraction"] for x in lit["compositionLand"]}
        out[d.name] = {
            "terrain_hash": m["hashes"]["finalElevation"],
            "catalogue_hash": m["hashes"]["basinCatalogue"],
            "seed": m["seed"],
            "regions": m["numRegions"],
            "selection_source": b.get("selectionSource"),
            "drainage_hypothesis": {k: v for k, v in
                                    (b.get("drainageHypothesis") or {}).items()
                                    if k != "entries"},
            "preserved_basins": b["counts"]["preserved"],
            "land_fraction": m.get("landSeaMask", {}).get("landFractionBySurfaceClass"),
            "endorheic_fraction_of_land":
                (b.get("drainageConsistency") or {}).get("fractionOfLand"),
            "lithology_land_fractions": {k: round(v, 5) for k, v in
                                         sorted(comp.items(), key=lambda kv: -kv[1])[:6]},
            "grids": sorted(p.name.replace("exoplasim-", "") for p in d.iterdir()
                            if p.is_dir() and p.name.startswith("exoplasim-")),
            "superseded": (d / "SUPERSEDED.md").is_file(),
        }
    return out


def climate_runs(active: str, index: list) -> list:
    """Runs on the active build only.

    Runs on a superseded terrain describe a world we no longer model. They stay
    on disk and stay in `exoplasim/runs/INDEX.json`, which is the record of what
    exists; this file is the record of what is true.
    """
    rows = []
    for row in index:
        if row.get("source_build") != active:
            continue
        run = ROOT / "exoplasim" / "runs" / row["directory"]
        d = json.loads((run / "run_manifest.json").read_text(encoding="utf-8"))
        conv = d.get("convergence_assessment") or {}
        # The stellar-cycle runs use fixed_orbit_parameters rather than
        # derived_parameters, and the earliest runs predate several keys, so read
        # everything defensively: a state file that crashes on an old artifact is
        # worse than one that reports a null.
        derived = d.get("derived_parameters") or d.get("fixed_orbit_parameters") or {}
        cfg = d.get("source_config") or {}
        model = cfg.get("model") or {}
        rows.append({
            "run_id": d.get("run_id", run.name),
            "status": d.get("status"),
            "orbits_complete": row.get("orbits_on_disk"),
            "flux_earth": derived.get("stellar_flux_ratio_earth"),
            "cycle": d.get("cycle", {}).get("case") if isinstance(d.get("cycle"), dict) else None,
            "resolution": model.get("resolution"),
            "land_albedo_source": model.get("land_albedo_source"),
            "source_build": cfg.get("source_build"),
            "glaciers": (cfg.get("surface") or {}).get("glaciers", {}).get("enabled"),
            "converged": conv.get("pass"),
            "mean_surface_temperature_k": (conv.get("metrics") or {}).get(
                "temperature_mean_k") or _mean_ts(run),
        })
    return rows


def _mean_ts(run: Path, last: int = 5) -> float | None:
    """Global mean surface temperature over the final orbits, from the output.

    Computed here rather than taken from a convergence assessment, because a run
    has output long before anything assesses it, and a state file reporting null
    for a finished run is worse than one that does the arithmetic.
    """
    # Enumerating one already-resolved run's own output files. This is not
    # artifact selection -- the run was chosen through the index -- but it is
    # still the last directory read in this file, and it exists only because a
    # run has output before anything assesses it.
    files = sorted(run.glob("MOST.*.nc"))
    if not files:
        return None
    try:
        from netCDF4 import Dataset
        from numpy.polynomial.legendre import leggauss
        vals = []
        for f in files[-last:]:
            with Dataset(f) as ds:
                w = leggauss(ds.dimensions["lat"].size)[1][::-1]
                a = np.asarray(ds["ts"][:])
                while a.ndim > 2:
                    a = a.mean(axis=0)
                vals.append(float((a.mean(axis=1) * w).sum() / w.sum()))
        return round(float(np.mean(vals)), 3)
    except Exception:
        return None


def hydrography_dir() -> Path:
    """Hydrography products for the configured build.

    These are per-build, because drainage is a property of the terrain. The flat
    `hydrography/data/` predates that and still holds whichever build was active
    when it was written, so preferring it would report another terrain's basins
    against the current one. That is how this file came to say 2,107 basins and
    1,522 carves while the active build had 2,540 and 1,089.
    """
    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    named = ROOT / "hydrography" / "data" / str(cfg.get("source_build", ""))
    return named if (named / "hydrography_report.json").is_file() else (
        ROOT / "hydrography" / "data")


def hydrography() -> dict:
    data = hydrography_dir()
    rep = read_json(data / "hydrography_report.json") or {}
    verdict = read_json(ROOT / "hydrography" / "analysis" / "carve_verdict.json") or {}
    carve = read_json(data / "carve_list.json") or {}
    # No fallback to another build's list. An earlier version globbed for any
    # carve_list.json and took the first alphabetically, which reported one
    # terrain's verdict against another's basins -- the same silent pairing this
    # file was fixed for once already. Absent means absent.
    return {
        "products_from": str(data.relative_to(ROOT)),
        # Two different quantities have both been called "the endorheic share"
        # and they differ by a factor of three or so, because a basin's catchment
        # is far larger than its floor. Named apart here before someone quotes
        # one as the other.
        "endorheic_note": {
            "drainage_share_of_land": "fraction of LAND AREA that drains to a "
                "closed basin; the hydrography figure, and the one the carve "
                "verdict is about",
            "basin_floor_share_of_land": "fraction of LAND AREA inside a "
                "preserved basin, i.e. is_endorheic cells; the lithology and "
                "thermostat figure",
        },
        "terrain_hash": (rep.get("source") or {}).get("terrain_hash"),
        "drainage": rep.get("drainage"),
        "basins": rep.get("basins"),
        "carve_verdict": {
            "bounds": verdict.get("bounds"),
            "agreement": verdict.get("agreement"),
            "penman_ocean_validation": verdict.get("penman_ocean_validation"),
        },
        "carve_list": {
            "counts": carve.get("counts"),
            "terrain_hash": carve.get("terrain_hash"),
            "climate": carve.get("climate"),
        } if carve else None,
    }


def cleared(field: str, depends_on: str, expected, found) -> dict:
    """A field whose dependency has moved, reported as absent rather than stale.

    A derived value computed from a superseded terrain is not a fact about this
    world, and carrying it with a prose caveat does not help: it stays readable,
    quotable and wrong. `current_climate` did exactly that -- it reported the
    pre-carve, wrong-spectrum climatology with a paragraph explaining not to
    quote it, which is an invitation rather than a guard.
    """
    return {
        "value": None,
        "invalidated": True,
        "depends_on": depends_on,
        "expected": expected,
        "found": found,
        "note": f"{field} was computed from a dependency that has since moved; "
                "regenerate it. Cleared rather than shown, because a stale "
                "derived value is indistinguishable from a current one.",
    }


def run_index() -> list:
    """The run set, from its manifest. Never from a directory listing.

    `exoplasim/runs/INDEX.json` is written by `index_runs.py`, which is the one
    place allowed to enumerate the runs directory. Everything downstream reads
    the index, so no consumer can pick an artifact by whatever sorts last -- the
    pattern behind ExoPlaSim's `finalize()` emitting the wrong world's output,
    behind this file previously reporting one terrain's verdict against
    another's basins, and behind `current_climate` reporting a superseded
    climatology because it happened to sort last.
    """
    idx = ROOT / "exoplasim" / "runs" / "INDEX.json"
    if not idx.is_file():
        raise SystemExit(
            "exoplasim/runs/INDEX.json is missing. Run "
            "`python exoplasim/scripts/index_runs.py` -- world_state resolves "
            "runs through that index and will not enumerate the directory.")
    return json.loads(idx.read_text(encoding="utf-8"))["runs"]


def climatology_build(regular_path: str | None, index: list) -> str | None:
    """Which build produced a climatology, resolved through the run index."""
    if not regular_path:
        return None
    name = Path(regular_path).name
    for row in index:
        for entry in (row.get("climatologies") or {}).values():
            if entry.get("regular") == name:
                return row.get("source_build")
    return None


def current_climate(active: str, index: list) -> dict | None:
    """The climate report belonging to the ACTIVE build, resolved by provenance.

    Candidate reports are named by the climatologies the run index says exist,
    not discovered by walking the analysis directory.
    """
    reports = []
    for row in index:
        for label in (row.get("climatologies") or {}):
            for cand in (ROOT / "exoplasim" / "analysis" / label /
                         "baseline_climate_report.json",
                         ROOT / "exoplasim" / "analysis" / "climatology" /
                         f"{label}_baseline_climate_report.json"):
                if cand.is_file():
                    reports.append(cand)
    if not reports:
        return cleared("current_climate", "source_build", active,
                       "no baseline_climate_report.json exists for any "
                       "climatology this build produced; run "
                       "exoplasim/scripts/analyze_climatology.py")
    # Pick by provenance, not by glob order. Taking the last match is how
    # ExoPlaSim's own finalize() emits the wrong world's result, and this
    # function had the same bug: it reported whichever report sorted last.
    latest = None
    for cand in reports:
        c = json.loads(cand.read_text(encoding="utf-8"))
        if climatology_build(c.get("source_regular"), index) == active:
            latest = cand
            break
    if latest is None:
        got = {str(c.relative_to(ROOT)): climatology_build(
                   json.loads(c.read_text(encoding="utf-8")).get("source_regular"),
                   index) for c in reports}
        return cleared("current_climate", "source_build", active, got)
    d = json.loads(latest.read_text(encoding="utf-8"))
    koppen = d.get("koppen_land_area_fractions") or {}
    biome = d.get("biome_land_area_fractions") or {}
    return {
        "source": str(latest.relative_to(ROOT)),
        "global_metrics": d.get("global_metrics"),
        "koppen_top": dict(sorted(koppen.items(), key=lambda kv: -kv[1])[:6]),
        "biome_top": dict(sorted(biome.items(), key=lambda kv: -kv[1])[:6]),
        "source_build": active,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=ROOT / "world_state.json")
    args = ap.parse_args()
    index = run_index()

    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    state = {
        "schema_version": SCHEMA_VERSION,
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "generator": "scripts/world_state.py",
        "note": "Generated from the artifacts. Do not edit; re-run the generator. "
                "Judgement that no artifact can supply is in `curated`, which "
                "lives in the generator under version control.",
        "planet": config["planet"],
        "star": config["star"],
        "orbit": config["orbit"],
        "atmosphere": config["atmosphere"],
        "active_configuration": {
            "source_build": config.get("source_build"),
            "model": config["model"],
            "surface": config["surface"],
            "radiation": config["radiation"],
            "stellar_cycle": config.get("stellar_cycle"),
        },
        "builds": builds(config.get("source_build")),
        "climate_runs": climate_runs(config.get("source_build"), index),
        "current_climate": current_climate(config.get("source_build"), index),
        "hydrography": hydrography(),
        "curated": CURATED,
    }
    args.output.write_text(json.dumps(state, indent=2, default=str) + "\n",
                           encoding="utf-8")

    b = state["builds"]
    print(f"active build  : {config.get('source_build')}")
    print(f"climate runs  : {len(state['climate_runs'])}, "
          f"{sum(1 for r in state['climate_runs'] if r['converged'])} converged")
    hy = state["hydrography"]
    if hy.get("basins"):
        print(f"basins        : {hy['basins']['count']}")
    if hy.get("carve_list") and hy["carve_list"].get("counts"):
        print(f"carve list    : {hy['carve_list']['counts']}")
    print(f"open questions: {len(CURATED['open_questions'])}, "
          f"uncertainties: {len(CURATED['known_uncertainties'])}")
    print(f"\nwrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
