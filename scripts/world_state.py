#!/usr/bin/env python3
"""Generate world_state.json: what this project currently knows about Vesper.

Generated, not hand-maintained. Every prose summary in this project has drifted
from its artifacts at least once, usually within a commit of the number changing,
so the state file reads the artifacts instead of restating them. Everything it
writes is derived from an artifact: what is judged rather than measured belongs
in a notes document, and what is to be done about it in the `bd` issue tracker.

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
import gridding  # noqa: E402  from lib/, via _paths: the one Gaussian quadrature

ROOT = Path(__file__).resolve().parents[1]
import sys as _sys
if str(ROOT / "lib") not in _sys.path:
    _sys.path.insert(0, str(ROOT / "lib"))
from paths import rel  # noqa: E402
# Aliased: there is a `builds()` reporter below, and the bare name would be
# rebound by it -- which `scripts/smoke_test.py` lints for.
import builds as build_registry  # noqa: E402
from provenance import build_stamp  # noqa: E402
import stellar  # noqa: E402
SCHEMA_VERSION = 1

# No `curated` block here. This file is generated state, and judgement is not
# state: every item that block carried was either already argued in a notes
# document, or an action, and in both cases the copy here went stale without
# anything failing. It still posed the carved-zoned albedo arithmetic as an open
# question two documents after it had been measured at 295.18 K, and it still
# reported the weathering bracket at a width pedology/README.md had already
# corrected. Findings live in `notes/` and `notes/audits/` with their evidence,
# what to do about them lives in the `bd` issue tracker, and both are read by
# than regenerated over.


def stellar_spectrum() -> dict:
    """The band split every consumer resolves through, and its evidence.

    Derived, not declared: `lib/stellar.py` reproduces `radmod.f90:solarini` on
    the configured spectrum, and the model computes the same number from the
    same file. It lives here rather than in `config/planet.yaml` because it is a
    property of the spectrum file, not a decision.

    `rcoeff_as_solarini_codes_it` is recorded beside the consistent value
    because they differ by 3.4x and only the second is physical; SPEC-2.
    """
    import hashlib
    low, high = stellar.spectrum_paths()
    band1, band2 = stellar.band_fractions(path=high)
    return {
        "spectrum": high.stem.removesuffix("_hr"),
        "hires_file": rel(high),
        "hires_sha256": hashlib.sha256(high.read_bytes()).hexdigest(),
        "generator": "lib/stellar.py, reproducing radmod.f90:solarini",
        "band_split_um": stellar.BAND_SPLIT_UM,
        "min_wavelength_nm": stellar.MIN_WAVELENGTH_NM,
        "flux_fraction_band1": band1,
        "flux_fraction_band2": band2,
        "solar_partition_identity": stellar.solar_partition_identity(),
        "rayleigh_coefficient": stellar.rayleigh_coefficient(),
        "rayleigh_coefficient_as_solarini_codes_it":
            stellar.rayleigh_coefficient(as_the_model_does=True),
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
        try:
            man = build_registry.mesh_export_of(d) / "manifest.json"
        except RuntimeError:
            continue          # a stub build with no payload
        if not man.is_file():
            continue
        m = json.loads(man.read_text(encoding="utf-8"))
        b, lit = m["basins"], m["lithology"]
        comp = {x["code"]: x["fraction"] for x in lit["compositionLand"]}
        out[d.name] = {
            # WORLD-8H15. world_state.json is the file everything quotes
            # the current build from, and it reported the build without
            # the registry's verdict on it. The registry can REFUSE the
            # build config names -- registration and activation are
            # separate on purpose -- so the state has to carry which it is
            # rather than leaving a reader to go and look.
            "build_verdict": build_stamp(d.name),
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
            # TWO endorheic quantities, and the key used to name neither. A
            # basin's CATCHMENT is far larger than its FLOOR, and the lithology
            # depends on the floor because Orogen assigns basin fill on
            # `is_endorheic`. Reporting one under an ambiguous name is how a
            # retention argument came to quote a catchment share as though it
            # were the fill. CONS-12 asked for this rename.
            "endorheic_catchment_share_of_land":
                (b.get("drainageConsistency") or {}).get("fractionOfLand"),
            "endorheic_basin_floor_share_of_land": basin_floor_share(d),
            # Every class, not the top six. Truncation hid `evaporite`, which
            # ranks tenth, so basin fill -- playa clastic plus evaporite, the
            # quantity a carve moves -- was not reconstructible from this file.
            "lithology_land_fractions": {k: round(v, 5) for k, v in
                                         sorted(comp.items(), key=lambda kv: -kv[1])},
            # Land elevation, area-weighted over surface_class land. Here
            # because it was being written into README prose instead, where it
            # survived a gravity correction that changed it by a quarter and
            # went on reading as current. reliefScale is applied at export, so
            # these move whenever gravity does even though no hash changes.
            "land_elevation_km": land_elevation(d),
            "grids": sorted(p.name.replace("exoplasim-", "") for p in d.iterdir()
                            if p.is_dir() and p.name.startswith("exoplasim-")),
            "superseded": (d / "SUPERSEDED.md").is_file(),
        }
    return out


def basin_floor_share(build_dir: Path) -> float | None:
    """Area-weighted share of land that is closed-basin FLOOR, not catchment.

    The lithology's basin fill is assigned on `is_endorheic`, so this is the
    quantity a fill fraction has to be read against. Land comes from
    `surface_class` per CLAUDE.md rule 1: `land_mask` drops the dry floors
    below sea level, which are exactly the cells this measures.

    A PRE-CARVE BUILD REPORTS A LIMIT, NOT A STATE. The carve list moves
    `is_endorheic` directly, so this falls as basins are opened.
    """
    try:
        root = build_registry.mesh_export_of(build_dir)
        man = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        raw = man["raw"]["fields"]
        names = raw if isinstance(raw, dict) else {f["name"]: f for f in raw}

        def field(name: str):
            spec = names[name]
            return np.fromfile(root / spec["path"], dtype=spec["dtype"])

        area = field("cell_area").astype(float)
        land = field("surface_class") == 1
        endo = field("is_endorheic").astype(bool)
        if not land.any():
            return None
        return round(float(area[land & endo].sum() / area[land].sum()), 6)
    except Exception:
        # Same as land_elevation: a missing payload is not a state error.
        return None


def land_elevation(build_dir: Path) -> dict | None:
    """Area-weighted land elevation from the native mesh, in physical km.

    From `elevation_km`, which carries the 1/g relief scaling, and masked by
    `surface_class` rather than `land_mask` -- the latter would drop the dry
    below-sea-level basin floors, which is most of what makes this world's
    hypsometry unusual.
    """
    try:
        root = build_registry.mesh_export_of(build_dir)
        man = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        raw = man["raw"]["fields"]
        names = raw if isinstance(raw, dict) else {f["name"]: f for f in raw}

        def field(name: str):
            spec = names[name]
            return np.fromfile(root / spec["path"], dtype=spec["dtype"])

        z = field("elevation_km").astype(float)
        area = field("cell_area").astype(float)
        land = field("surface_class") == 1
        if not land.any():
            return None
        return {
            "mean": round(float(np.average(z[land], weights=area[land])), 4),
            "max": round(float(z[land].max()), 4),
            "min": round(float(z[land].min()), 4),
        }
    except Exception:
        # The payload is deleted for every build but the active one, and a
        # missing mesh is not a state error.
        return None


def climate_runs(active: str, index: list) -> list:
    """Runs on the active build only.

    Runs on a superseded terrain describe a world we no longer model. They stay
    on disk and stay in `exoplasim/runs/INDEX.json`, which is the record of
    every run that has EXISTED; this file is the record of what is true.

    So a row is skipped unless its payload is still there. The index is a ledger
    (world-ww6z): it keeps the row of a run whose directory has gone, which is
    what makes a deleted run identifiable, and a row whose manifest is no longer
    readable cannot answer what this file asks of it. Its identity is in
    `archive/runs/`, which is where a reader who wants it should go.
    """
    rows = []
    for row in index:
        if row.get("source_build") != active:
            continue
        if not row.get("payload_present", True):
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
        import climatology as clim
        vals = []
        for f in files[-last:]:
            with Dataset(f) as ds:
                w = gridding.gaussian_row_weights(
                    np.asarray(ds["lat"][:], dtype=float),
                    what=f"{f}'s latitude axis")
                a = np.asarray(ds["ts"][:])
                if a.ndim > 2:      # leading axis is time, and its bins differ
                    a = clim.annual_mean(a, np.asarray(ds["time"][:]))
                if a.ndim != 2:
                    # `ts` is a surface field, so the time mean leaves (lat, lon)
                    # and there is nothing further to reduce. An extra axis here
                    # is a level axis, and collapsing one with an unweighted mean
                    # is a mass weighting nobody chose; refuse instead.
                    raise ValueError(
                        f"{f.name}: ts is {a.ndim}-dimensional after the annual "
                        "mean. The remaining axis is not time and this has no "
                        "weighting for it")
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
    return named


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
            # The last two are the old layout and the old doubled name,
            # kept so reports already on disk stay findable.
            for cand in (ROOT / "exoplasim" / "analysis" / "climatology" /
                         f"{label}_climate_report.json",
                         ROOT / "exoplasim" / "analysis" / label /
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
                "Only derived values live here. Findings are in `notes/`, with "
                "their evidence, and open work is in the `bd` issue tracker.",
        "planet": config["planet"],
        "star": config["star"],
        "orbit": config["orbit"],
        "atmosphere": config["atmosphere"],
        "stellar_spectrum": stellar_spectrum(),
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
    print(f"\nwrote {rel(args.output)}")


if __name__ == "__main__":
    main()
