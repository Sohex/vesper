#!/usr/bin/env python3
"""Write the carve verdict as a World Orogen `--preserve-basins` list.

Three outcomes per basin, and the middle one is the point:

  carved     Both evaporation estimates say the basin overflows. Retain 0, which
             Orogen treats as not preserved, so drainage enforcement carves the
             outlet as it would have without preservation at all.
  preserved  Neither estimate says it overflows. Retain 1, rim intact.
  marginal   The two estimates disagree. Retain between 0 and 1, which cuts a
             notch at the saddle and tapers it over the divide band, producing a
             through-flowing valley with a residual lake.

The retain fraction measures how close a basin is to its own threshold, in units
of the uncertainty we actually have.

A basin balances exactly at an evaporation of `E* = P + critical * runoff`. Under
the Penman estimate it evaporates `E_penman`. The fractional margin

    margin = (E_penman - E*) / E_penman
    retain = min(1, margin / TOLERANCE)

says how much open-water evaporation would have to fall before the basin starts
overflowing. A basin needing a 3% change is genuinely marginal and keeps little
of its rim; one needing 40% is comfortably closed and keeps all of it.

TOLERANCE is 0.25, set from the uncertainty that actually dominates: the biosphere
is assumed rather than modelled, and that assumption is worth 3.7 to 7.1 K, which
moves evaporation by considerably more than the 1.7% Penman itself was validated
to.

An earlier mapping interpolated the critical index between the two evaporation
estimates. It is retained in the sidecar as `retain_span` but is not used, because
it measures the wrong thing: the land-evaporation end is a deliberately weak lower
bound, so the span is enormous and the critical value lands near its low end
almost regardless of the basin. It put 70% of marginal basins above retain 0.9,
which is to say it called them marginal and then declined to cut them.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

import carve_verdict as cv
from _paths import CONFIG, DATA, PROJECT_ROOT
from orbit import orbital_year_days
from lake_balance import BasinSet


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climatology", type=Path, default=PROJECT_ROOT /
                    "exoplasim/analysis/climatology_s096/baseline_regular_climatology.nc")
    ap.add_argument("--coupling", type=Path, default=DATA / "coupling_exoplasim-T42.nc")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--out-list", type=Path, default=DATA / "carve_list.txt")
    ap.add_argument("--out-json", type=Path, default=DATA / "carve_list.json")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    basins = BasinSet()
    resolution = str(config["model"]["resolution"]).upper()

    with Dataset(args.climatology) as ds:
        am = cv.annual_mean
        pr, evap, mrro = am(ds, "pr"), -am(ds, "evap"), am(ds, "mrro")
        ts, tas = am(ds, "ts"), am(ds, "tas")
        ps_pa, rss, rls = am(ds, "ps") * 100.0, am(ds, "rss"), am(ds, "rls")
        q_air = np.asarray(ds["hus"][:]).mean(axis=0)[-1]
        wind = np.asarray(ds["spd"][:]).mean(axis=0)[-1]

    land_albedo = cv.read_sra_field(
        PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
        / f"orogen_{resolution}_surf_0174.sra", *ps_pa.shape)
    penman = np.maximum(cv.penman_open_water(
        ts, tas, q_air, wind, ps_pa, rss, rls, land_albedo,
        float(config["planet"]["gravity_m_s2"])), evap)

    means, _ = cv.basin_means(args.coupling,
                              {"pr": pr, "wet": evap, "pen": penman, "ro": mrro},
                              basins.n)
    # Orbital period varies with flux, so take it from the config rather than
    # hardcoding. The literal here was 189.6145 d, the 0.90-flux year, while this
    # baseline runs at 0.96 and 180.655 d. It cancels out of the aridity index
    # and the equilibrium lake area, both ratios, but it was making the reported
    # runoff depth 5% high.
    year_s = orbital_year_days(config) * 86400.0
    runoff = means["ro"] * year_s / 1000.0
    precip = means["pr"] * year_s / 1000.0
    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0

    def index(field):
        e = means[field] * year_s / 1000.0
        return np.where(runoff > 0, (e - precip) / np.where(runoff > 0, runoff, 1.0),
                        np.inf)

    idx_wet, idx_pen = index("wet"), index("pen")
    carve_wet, carve_pen = idx_wet <= crit, idx_pen <= crit
    carved = carve_pen                      # nests inside carve_wet
    preserved = ~carve_wet
    marginal = carve_wet & ~carve_pen
    assert int(carved.sum() + preserved.sum() + marginal.sum()) == basins.n

    TOLERANCE = 0.25
    e_pen = means["pen"] * year_s / 1000.0
    e_balance = precip + crit * runoff              # evaporation that exactly balances
    with np.errstate(divide="ignore", invalid="ignore"):
        margin = np.where(e_pen > 0, (e_pen - e_balance) / np.where(e_pen > 0, e_pen, 1.0), 1.0)
    retain = np.clip(margin / TOLERANCE, 0.0, 1.0)
    # A basin with no catchment runoff receives nothing and can never overflow,
    # whatever the evaporation is. The margin expression does not know that: with
    # runoff zero it reduces to a comparison of evaporation against precipitation
    # alone, which sent 228 dry salt pans to be carved wide open.
    retain = np.where(runoff > 0, retain, 1.0)
    retain = np.where(carved, 0.0, retain)

    # The superseded mapping, kept for comparison in the sidecar only.
    span = idx_pen - idx_wet
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(np.isfinite(span) & (span > 1e-9),
                     (crit - idx_wet) / np.where(span > 1e-9, span, 1.0), 0.5)
    retain_span = np.where(carved, 0.0, np.where(preserved, 1.0, np.clip(1.0 - f, 0.0, 1.0)))

    # Verdict now follows retain, since a basin can be preserved by the Penman
    # test yet sit close enough to its threshold to warrant a partial cut.
    verdict_name = np.where(carved, "carve",
                            np.where(retain >= 1.0, "preserve", "marginal"))
    n_carve = int(carved.sum())
    n_preserve = int((~carved & (retain >= 1.0)).sum())
    n_marginal = basins.n - n_carve - n_preserve

    with Dataset(DATA / "basins.nc") as ds:
        ids = [str(x) for x in ds["basin_id"][:]]

    header = f"""# Vesper carve verdict, first pass
#
# Produced from a converged ExoPlaSim climatology: {resolution}, 0.96 S-Earth,
# vegetated land surface, glaciers enabled, measured K2 spectrum, 292.88 K.
# Terrain {basins.terrain_hash[:16]}, basin catalogue unchanged.
#
# A basin overflows, and so should have its outlet carved, when
#     (E - P) / runoff  <=  catchment / area_at_spill - 1
# Open-water evaporation is the Penman combination equation with water's albedo
# and roughness, validated against the model over ocean cells to within 1.7%.
#
#   retain 1.0   {n_preserve:4d} basins  comfortably closed
#   retain 0<r<1 {n_marginal:4d} basins  within 25% of their threshold
#   retain 0.0   {n_carve:4d} basins  overflow under both estimates
#
# Carved basins are listed explicitly at retain 0 rather than omitted, so this
# file is the complete verdict rather than a subset of it. If your parser would
# rather they were absent, say so and we will drop them.
#
# id                                    retain
"""
    lines = [header]
    for i in np.argsort(-retain):
        lines.append(f"{ids[i]:<38s} {retain[i]:.4f}")
    args.out_list.parent.mkdir(parents=True, exist_ok=True)
    args.out_list.write_text("\n".join(lines) + "\n", encoding="ascii")

    sidecar = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "terrain_hash": basins.terrain_hash,
        "climatology": str(args.climatology),
        "climate": {
            "resolution": resolution, "flux_earth": 0.96,
            "land_surface": "vegetated", "mean_surface_temperature_k": 292.88,
            "note": "Converged on all six criteria. The biosphere is assumed, "
                    "not modelled, and that assumption is worth 3.7 to 7.1 K.",
        },
        "method": {
            "test": "(E - P) / runoff <= catchment / area_at_spill - 1",
            "open_water_evaporation": "Penman combination, water albedo and "
                                      "roughness, floored at the model's land rate",
            "penman_ocean_validation_ratio": 1.017,
            "retain_mapping": "min(1, ((E_penman - (P + critical*runoff)) / E_penman) / 0.25)",
            "retain_tolerance": TOLERANCE,
            "retain_tolerance_rationale": "fractional change in open-water "
                "evaporation that would flip the verdict; 0.25 is set by the "
                "assumed-biosphere uncertainty, which dominates Penman's own 1.7%",
        },
        "counts": {"carve": n_carve, "preserve": n_preserve, "marginal": n_marginal},
        "counts_by_penman_test_alone": {
            "carve": int(carved.sum()), "preserve": int(preserved.sum()),
            "disagreeing_with_land_evaporation": int(marginal.sum())},
        "basins": [
            {
                "id": ids[i],
                "verdict": str(verdict_name[i]),
                "retain": round(float(retain[i]), 4),
                "retain_span_superseded": round(float(retain_span[i]), 4),
                "evaporation_margin": None if not np.isfinite(margin[i])
                                      else round(float(margin[i]), 4),
                "critical_aridity_index": round(float(crit[i]), 4),
                "aridity_index_penman": None if not np.isfinite(idx_pen[i])
                                        else round(float(idx_pen[i]), 4),
                "aridity_index_land_evaporation": None if not np.isfinite(idx_wet[i])
                                                  else round(float(idx_wet[i]), 4),
                "catchment_km2": round(float(basins.catchment_km2[i]), 2),
                "area_at_spill_km2": round(float(basins.area_at_spill_km2[i]), 2),
                "runoff_km_per_year": round(float(runoff[i]), 8),
                "no_catchment_runoff": bool(runoff[i] <= 0),
            } for i in range(basins.n)
        ],
    }
    args.out_json.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    mid = (~carved) & (retain < 1.0)
    print(f"carve     {n_carve:5d}  retain 0.0")
    print(f"marginal  {n_marginal:5d}  retain {retain[mid].min():.3f}"
          f"-{retain[mid].max():.3f}, median {np.median(retain[mid]):.3f}")
    print(f"preserve  {n_preserve:5d}  retain 1.0")
    print(f"\nwrote {args.out_list}\n      {args.out_json}")


if __name__ == "__main__":
    main()
