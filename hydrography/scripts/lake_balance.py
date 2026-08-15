#!/usr/bin/env python3
"""Solve endorheic lake levels from a climate forcing.

The machinery here is climate-independent; the forcing is not. Nothing in this
file decides what the climate is. It takes runoff and evaporation as arguments
and returns the lake each basin settles at, which is the step that has to wait
for ExoPlaSim.

Run directly to exercise the solver against a uniform placeholder forcing. That
is a smoke test of the machinery and a sensitivity sweep, not a result about
this world.

The balance
-----------

For a basin with catchment area C and lake area A, at equilibrium the water
arriving equals the water leaving::

    r (C - A)  +  p A  =  e A

where r is runoff depth generated per unit of dry catchment, p is precipitation
falling directly on the lake, and e is evaporation from open water, all as
depths per unit time. Solving for area::

    A = r C / (e - p + r)

so the lake is set by the ratio of catchment supply to the net evaporative
demand over water. If ``e - p + r <= 0`` the lake can never evaporate what it
receives and the basin fills to its spill and overflows regardless of size.

Overflow cascades. On this terrain 2,640 of 3,629 basins spill into another
basin rather than to the ocean, so a full basin passes its surplus downstream
and the system has to be iterated to a fixed point.
"""

from __future__ import annotations

import argparse
import json

from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS, DATA

TERMINAL_OCEAN = -1


class BasinSet:
    """Basin hypsometry and topology, as built by `build_hydrography.py`."""

    def __init__(self, path=None):
        path = path or (DATA / "basins.nc")
        with Dataset(path) as ds:
            self.level_km = np.asarray(ds["level_km"][:])
            self.area_km2 = np.asarray(ds["flooded_area_km2"][:])
            self.volume_km3 = np.asarray(ds["volume_km3"][:])
            self.spill_km = np.asarray(ds["spill_km"][:])
            self.spill_target = np.asarray(ds["spill_target"][:])
            self.catchment_km2 = np.asarray(ds["catchment_km2"][:])
            self.capacity_km3 = np.asarray(ds["capacity_km3"][:])
            self.area_at_spill_km2 = np.asarray(ds["area_at_spill_km2"][:])
            self.terrain_hash = ds.terrain_hash
        self.n = self.level_km.shape[0]

    def level_for_area(self, basin: int, area_km2: float) -> float:
        """Invert the hypsometric curve. Area is monotone in level."""
        a, l = self.area_km2[basin], self.level_km[basin]
        return float(np.interp(area_km2, a, l))

    def volume_for_area(self, basin: int, area_km2: float) -> float:
        a, v = self.area_km2[basin], self.volume_km3[basin]
        return float(np.interp(area_km2, a, v))


def solve(
    basins: BasinSet,
    runoff_km_per_year: np.ndarray,
    lake_evaporation_km_per_year: np.ndarray,
    lake_precip_km_per_year: np.ndarray,
    *,
    max_iterations: int = 200,
    tolerance_km3: float = 1e-3,
):
    """Equilibrium lake area, level and volume for every basin.

    All three forcings are per-basin depths per year: `runoff_km_per_year` is
    generated over the dry catchment, the other two apply over open water.

    Returns a dict of arrays plus `converged`, which is False if the overflow
    cascade did not reach a fixed point. Do not quietly use a non-converged
    result; it means basins are exchanging overflow in a cycle.
    """
    n = basins.n
    supply = runoff_km_per_year * basins.catchment_km2      # km3/yr from land
    demand = lake_evaporation_km_per_year - lake_precip_km_per_year + runoff_km_per_year

    inflow_extra = np.zeros(n)      # overflow received from upstream basins
    area = np.zeros(n)
    overflow = np.zeros(n)
    converged = False

    for _ in range(max_iterations):
        total_supply = supply + inflow_extra
        with np.errstate(divide="ignore", invalid="ignore"):
            want = np.where(demand > 0, total_supply / demand, np.inf)
        area = np.minimum(np.nan_to_num(want, posinf=np.inf), basins.area_at_spill_km2)

        # A basin pinned at its spill passes on whatever it cannot evaporate.
        at_spill = area >= basins.area_at_spill_km2 - 1e-9
        overflow_new = np.where(
            at_spill,
            np.maximum(total_supply - basins.area_at_spill_km2 * demand, 0.0),
            0.0,
        )

        received = np.zeros(n)
        downstream = basins.spill_target
        routed = (downstream >= 0) & (overflow_new > 0)
        np.add.at(received, downstream[routed], overflow_new[routed])

        if np.max(np.abs(received - inflow_extra)) < tolerance_km3:
            inflow_extra = received
            overflow = overflow_new
            converged = True
            break
        inflow_extra = received
        overflow = overflow_new

    level = np.array([basins.level_for_area(b, area[b]) for b in range(n)])
    volume = np.array([basins.volume_for_area(b, area[b]) for b in range(n)])
    return {
        "area_km2": area,
        "level_km": level,
        "volume_km3": volume,
        "overflow_km3_per_year": overflow,
        "fills_to_spill": area >= basins.area_at_spill_km2 - 1e-9,
        "dry": area <= 0,
        "converged": converged,
    }


def _sweep(basins: BasinSet, out_path) -> dict:
    """Uniform-forcing sensitivity sweep. Machinery check, not a result."""
    planet_km2 = 734_492_839.55
    rows = []
    for runoff_mm in (10.0, 50.0, 200.0):
        for evap_mm in (400.0, 800.0, 1600.0):
            r = np.full(basins.n, runoff_mm * 1e-6)         # mm/yr -> km/yr
            e = np.full(basins.n, evap_mm * 1e-6)
            p = np.full(basins.n, runoff_mm * 1e-6)          # placeholder
            s = solve(basins, r, e, p)
            rows.append({
                "runoff_mm_per_year": runoff_mm,
                "lake_evaporation_mm_per_year": evap_mm,
                "converged": bool(s["converged"]),
                "lake_area_fraction_of_planet": float(s["area_km2"].sum() / planet_km2),
                "lake_volume_km3": float(s["volume_km3"].sum()),
                "basins_filling_to_spill": int(s["fills_to_spill"].sum()),
                "basins_dry": int(s["dry"].sum()),
            })
    result = {
        "note": ("Uniform placeholder forcing. This exercises the solver and shows "
                 "how lake extent scales; it is not a claim about this world's "
                 "climate. Replace with ExoPlaSim runoff and evaporation "
                 "integrated over each basin's catchment via coupling_*.nc."),
        "terrain_hash": basins.terrain_hash,
        "basins": basins.n,
        "capacity_km3": float(basins.capacity_km3.sum()),
        "sweep": rows,
    }
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basins", default=None)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    basins = BasinSet(args.basins)
    out = ANALYSIS / "lake_balance_sweep.json" if args.output is None else __import__(
        "pathlib").Path(args.output)
    res = _sweep(basins, out)
    print(f"{basins.n} basins, total capacity {res['capacity_km3']:,.0f} km3")
    print(f"{'runoff':>8} {'lake evap':>10} {'lake area':>11} {'volume km3':>13} "
          f"{'at spill':>9} {'dry':>6} {'conv':>5}")
    for r in res["sweep"]:
        print(f"{r['runoff_mm_per_year']:7.0f}m {r['lake_evaporation_mm_per_year']:9.0f}m "
              f"{r['lake_area_fraction_of_planet']*100:10.3f}% "
              f"{r['lake_volume_km3']:13,.0f} {r['basins_filling_to_spill']:9d} "
              f"{r['basins_dry']:6d} {str(r['converged']):>5}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()


def carve_verdict(basins: BasinSet, aridity_index: np.ndarray, land_area_km2: np.ndarray):
    """Which basins overflow persistently, and so should not be basins at all.

    A basin that overflows year on year incises its outlet. Over the 1e4 to 1e6
    years a landscape needs to relax, that cuts the sill down and drains the
    lake, and the depression stops existing. So `solve()` returning "pinned at
    spill, overflowing forever" is a transient, not a landscape state.

    The verdict is climate-free per basin up to a single number. From the
    balance, a basin overflows exactly when

        (E - P) / runoff  <=  catchment / area_at_spill - 1

    and the right-hand side is `basins.critical_aridity_index`, pure geometry.
    Pass the left-hand side per basin, from climate, and get back the basins
    whose outlets the terrain should have carved.

    Returns (carve, endorheic_land_km2). The area accounts for catchments that
    pass through a carved basin on their way somewhere else, so it is the land
    that still has nowhere to drain once the carving is done. Divide by total
    land yourself; passing only basin-draining land as `land_area_km2` and
    dividing by its own sum silently answers a different question.
    """
    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0
    carve = np.asarray(aridity_index) <= crit

    eff = np.arange(basins.n)
    for _ in range(basins.n):
        nxt = np.where(carve[eff] & (basins.spill_target[eff] >= 0),
                       basins.spill_target[eff], eff)
        nxt = np.where(carve[eff] & (basins.spill_target[eff] < 0), -1, nxt)
        nxt = np.where(eff < 0, -1, nxt)
        if np.array_equal(nxt, eff):
            break
        eff = nxt
    else:
        raise RuntimeError("carve cascade did not settle")

    return carve, float(np.sum(np.asarray(land_area_km2)[eff >= 0]))
