"""Route released solutes down the drainage network, per mesh region and per basin.

Meybeck's table 2C gives what a rock RELEASES per litre. Turning that into a
quantity a basin receives needs two more things, and both live here so that
every consumer gets the same ones:

- **the water**, which is local runoff generated at a region, not the discharge
  passing through it. Discharge is the accumulated sum down the tree, so
  weighting by it counts every upstream region again at every downstream one:
  a catchment's outlet region would carry the whole catchment's weight, and the
  chemistry of one cell near the sink would dominate a basin. Local generation
  is what each region contributes exactly once.

- **the join to the climate**, which is `gridding.climatology_cells` and nothing
  else. Longitude is index for index across the Orogen/ExoPlaSim seam. That has
  silently matched the wrong cells on three separate scripts, so this module
  does not reimplement it.

What "discharge weighting" then means, precisely: a basin's inflow chemistry is
the flux-weighted mean concentration, `sum(c_i q_i) / sum(q_i)` over the regions
of its catchment, with `q_i` the local runoff volume. Area weighting is the same
expression with `q_i` replaced by area, which is the special case of uniform
runoff -- and runoff on this world is anything but uniform, spanning four orders
of magnitude across a single large catchment.

The absolute flux `sum(c_i q_i)` falls out of the same sum and is the quantity
silica supply is about, so the two tasks share one routine rather than two.
"""

from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np

from _paths import PROJECT_ROOT  # noqa: F401  (adds lib/ to sys.path)

from climatology import annual_mean
from gridding import climatology_cells
from orogen import LAND, Export

# Meybeck's table 2C columns, in the order every consumer here expects them.
SPECIES = ["sio2_umol_l", "ca_ueq_l", "mg_ueq_l", "na_ueq_l", "k_ueq_l",
           "cl_ueq_l", "so4_ueq_l", "hco3_ueq_l"]

EARTH_YEAR_DAYS = 365.2422


def region_runoff_m3_yr(export: Export, grid_dir: Path, climatology: Path):
    """Local runoff volume generated at each mesh region, m3 per Earth year.

    `P - E`, clamped at zero, and NOT the model's `mrro`. `mrro` is river-routed
    net divergence rather than local generation, so it is negative in places and
    nonzero over ocean; using it understates land runoff by 6.6x. That correction
    was made independently in two components and failed to reach two others, which
    is why it is stated at the one place that computes it here.

    Clamped at zero because a catchment delivers zero or more, never less. A cell
    whose evaporation exceeds its precipitation contributes no water; it does not
    contribute negative water and it certainly does not contribute negative
    solute. `notes/failure-modes.md` class 7.

    Returns (runoff_m3_yr, runoff_mm_yr_per_region, row, col).
    """
    with nc.Dataset(climatology) as ds:
        clim_lat = np.asarray(ds["lat"][:], dtype=float)
        # Records per bin are unequal; lib/climatology.py weights by them.
        centres = np.asarray(ds["time"][:], dtype=float)
        pr = annual_mean(np.asarray(ds["pr"][:], dtype=float), centres)
        evap = -annual_mean(np.asarray(ds["evap"][:], dtype=float), centres)
    # m/s -> mm per Earth year.
    to_mm_yr = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
    runoff_mm = np.maximum((pr - evap) * to_mm_yr, 0.0)

    row, col = climatology_cells(export, grid_dir, clim_lat)
    per_region_mm = runoff_mm[row, col]
    area_km2 = export.cell_area.astype(np.float64)
    # mm -> m over km2 -> m2.
    volume = per_region_mm * 1e-3 * area_km2 * 1e6
    volume[export.surface_class != LAND] = 0.0
    return volume, per_region_mm, row, col


def region_release(export: Export, mapping: dict[str, str],
                   release_table: dict, columns: dict[str, int]) -> np.ndarray:
    """Per-region release concentration vector, shaped (len(SPECIES), n_regions).

    Keyed on `substrate_class`, which is the top of the cover/basement stack and
    is what water actually runs over. Ocean regions keep zero.
    """
    manifest = export.manifest
    classes = {c["id"]: c["code"] for c in manifest["lithology"]["rockClasses"]}
    unmapped = sorted({c for c in classes.values()
                       if c != "water" and c not in mapping})
    if unmapped:
        raise SystemExit(f"rock classes with no Meybeck mapping: {unmapped}")

    rock = export.substrate_class.astype(np.int32)
    release = np.zeros((len(SPECIES), export.n_regions))
    for cid, code in classes.items():
        if code == "water":
            continue
        sel = rock == cid
        if not sel.any():
            continue
        row = release_table[mapping[code]]
        for si, sp in enumerate(SPECIES):
            release[si, sel] = row[columns[sp]]
    return release


def basin_totals(terminal: np.ndarray, weight: np.ndarray,
                 release: np.ndarray, n_basins: int):
    """Weighted mean concentration and total solute per basin.

    `terminal` is `regions.nc`, the priority flood's answer: preserved basin
    index, -1 at the world ocean, -2 not land. Only regions with a terminal
    basin count, and each counts once.

    Returns (mean_concentration, total_solute, total_weight), the first shaped
    (species, basin) and carrying NaN where a basin received no weight at all.
    `total_solute` is in concentration-units times whatever `weight` is in, so
    passing runoff in m3/yr and Meybeck's umol/l gives umol/yr after the litre
    conversion the caller owns.
    """
    valid = terminal >= 0
    w = np.where(valid, weight, 0.0)
    total_weight = np.bincount(terminal[valid], weights=w[valid],
                               minlength=n_basins)
    total = np.zeros((release.shape[0], n_basins))
    mean = np.zeros_like(total)
    for si in range(release.shape[0]):
        total[si] = np.bincount(terminal[valid], weights=(w * release[si])[valid],
                                minlength=n_basins)
        mean[si] = np.where(total_weight > 0,
                            total[si] / np.maximum(total_weight, 1e-30), np.nan)
    return mean, total, total_weight
