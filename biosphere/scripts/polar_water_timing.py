#!/usr/bin/env python3
"""Does the simulated polar cap's water cross the surface when a plant could use it?

    python biosphere/scripts/polar_water_timing.py
    python biosphere/scripts/polar_water_timing.py --run lpj_<id>

This is a note about a simulated world: every quantity below is a modelled field
of Vesper, and every comparison to Earth is a distance to report.

`biosphere/notes/polar-cover-cold-filter-and-capture.md` bounds what an adapted
flora could capture at the cap by adding two terms of the run's own water
balance: the bare-soil evaporation a canopy would shade, and the drainage a
faster root system would intercept. That bound is a MASS argument. It never
asked WHEN either term happens, and a term that leaves while the ground is
frozen is recoverable by nothing, so the bound would be an overstatement rather
than a ceiling. This measures the timing.

WHERE THE FORCING COMES FROM. The driver LPJ-GUESS actually read, not the
climatology re-read on another calendar. The driver carries the model's own
cells, its own twelve intervals of the 183-day simulation year and its own air
temperature, so the growing-season window measured here is the one the run
integrated. Reconstructing it from the climatology would cross a boundary the
comparison does not need to cross, and the two label months differently.

THE CONTROLS, and they are the reason a null here would be believable. The two
monthly tables are summed over the year and compared against the annual figures
the note already published from `aaet.out` and `mevap.out`; a reduction, a unit
or a cell-set error shows up as a mismatch there before any timing statement is
reached. The growing-season length falls out as a third: 3.39 months above
freezing and 2.87 above 5 degC, derived here from the driver, are what the note
measured from the climatology.

WHAT IT CANNOT MEASURE. `tot_runoff.out` is annual, so the drainage term's month
is inferred from the surface flux rather than measured, and the script says so
rather than reporting an inferred number beside measured ones.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
from pathlib import Path

import numpy as np

from _paths import PROJECT_ROOT, RUNS       # noqa: F401  (RUNS documents the root)

import lpj_output

# The accepted run, and the one the note's figures are on.
DEFAULT_RUN = "lpj_1e6a2b9ca51a4eff9592992cad96677b"

# The driver's per-cell record: four doubles of soil state plus one per physical
# layer, then the four forcing fields. The layer count is not in the header, so
# it is solved from the file size and the solution must be exact.
FORCING_FIELDS = 4
POLAR_LAT = 75.0

# Two brackets rather than one, because "could hold leaf" has no single
# temperature: months above freezing, and months above the model's own GDD base.
LOOSE_C, STRICT_C = 0.0, 5.0

# The verdict rule, fixed before the numbers are read. On the loose bracket:
STANDS, FICTIONAL = 0.70, 0.30

# The note's annual polar figures, as positive controls.
CONTROL = {"mevap.out": 29.14, "maet.out": 3.66}
CONTROL_TOL = 0.02


def solve_layers(size: int, ncells: int, nint: int, nyears: int) -> int:
    """The soil-layer count the driver was written with, from its size."""
    head = 8 + 40 + 64 + nyears * nint * FORCING_FIELDS * 8
    per, rem = divmod(size - head, ncells)
    fixed = 16 + 8 + 4 * 8 + FORCING_FIELDS * nyears * nint * 8
    layers, layer_rem = divmod(per - fixed, 8)
    if rem or layer_rem or layers < 1:
        raise SystemExit(f"driver layout does not close on {size} bytes")
    return layers


def read_driver(path: Path):
    raw = path.read_bytes()
    ncells, nint, year_days, nyears, _sub, _pad = struct.unpack("<6i", raw[8:32])
    layers = solve_layers(len(raw), ncells, nint, nyears)
    off = 8 + 40 + 64
    span = nyears * nint * FORCING_FIELDS * 8
    intervals = np.frombuffer(raw[off:off + span],
                              dtype="<f8").reshape(nyears, nint, FORCING_FIELDS)
    off += span
    stride = 16 + 8 + (4 + layers) * 8 + FORCING_FIELDS * nyears * nint * 8
    lon = np.empty(ncells)
    lat = np.empty(ncells)
    n = nyears * nint * 8
    tas = np.empty((ncells, nyears * nint))
    pr = np.empty((ncells, nyears * nint))
    for c in range(ncells):
        b = off + c * stride
        lon[c], lat[c] = struct.unpack("<dd", raw[b:b + 16])
        f = b + 16 + 8 + (4 + layers) * 8
        tas[c] = np.frombuffer(raw[f:f + n], dtype="<f8")
        pr[c] = np.frombuffer(raw[f + n:f + 2 * n], dtype="<f8")
    # `build_lpj_driver.py` subtracts KELVIN when it assembles tas, so the file
    # is already in degC. Converting again here put the cap below absolute zero,
    # which is how the error announced itself.
    return lon, lat, tas, pr, intervals[0, :, 2] / 86400.0, year_days, layers


def monthly(run: Path, table: str, keys):
    """The equilibrium-window reduction of a monthly table, per cell."""
    reduced = lpj_output.reduce_table(run / table)
    names = [n for n in reduced.names if n != "Total"]
    index = [reduced.names.index(m) for m in names]
    out = np.full((len(keys), len(names)), np.nan)
    for i, key in enumerate(keys):
        value = reduced.values.get(key)
        if value is not None:
            out[i] = np.asarray(value)[index]
    if not np.isfinite(out).all():
        raise SystemExit(f"{table}: {int((~np.isfinite(out)).any(axis=1).sum())} "
                         "driver cells are absent from the table")
    return out, names


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=DEFAULT_RUN)
    args = parser.parse_args()
    run = RUNS / args.run

    manifest = json.loads((run / "run_manifest.json").read_text())
    driver = Path(manifest["inputs"]["driver"]["path"])
    digest = lpj_output.sha256(driver)
    if digest != manifest["inputs"]["driver"]["sha256"]:
        raise SystemExit(f"{driver} changed since the run: {digest}")
    print(f"driver verified against {run.name}'s manifest: {driver.name}")

    lon, lat, tas, pr_rate, days, year_days, layers = read_driver(driver)
    print(f"driver: {len(lat)} cells, {len(days)} intervals of "
          f"{days[0]:.3f} days, {year_days}-day year, {layers} soil layers")

    polar = np.abs(lat) > POLAR_LAT      # BOTH caps; the note's band is 234 cells
    print(f"cells poleward of {POLAR_LAT:g} degrees, both caps: {polar.sum()}")

    # The tables print two decimals, so the driver's coordinates are matched at
    # two and the match is checked for collisions rather than assumed.
    keys = [(round(float(a), 2), round(float(b), 2)) for a, b in zip(lon, lat)]
    if len(set(keys)) != len(keys):
        raise SystemExit("driver coordinates collide at the table's precision")

    evap, months = monthly(run, "mevap.out", keys)
    transp, _ = monthly(run, "maet.out", keys)

    for name, table, target in (("mevap.out", evap, CONTROL["mevap.out"]),
                                ("maet.out", transp, CONTROL["maet.out"])):
        annual = float(np.nansum(table[polar], axis=1).mean())
        off = abs(annual - target) / target
        print(f"  control {name}: twelve months sum to {annual:.3f} mm against "
              f"the note's {target:.2f} ({off * 100:.1f}% off)")
        if off > CONTROL_TOL:
            raise SystemExit("positive control failed; the reduction, the units "
                             "or the cell set is wrong and nothing below holds")

    t = tas[polar]
    p = pr_rate[polar] * days[None, :]
    e = evap[polar]
    x = transp[polar]

    print("\n== the simulated polar year, by cap ==")
    for cap, sel in (("north", lat[polar] > 0), ("south", lat[polar] < 0)):
        print(f"  {cap} cap, {sel.sum()} cells")
        print(f"  {'month':>6} {'tas degC':>9} {'precip mm':>10} {'frozen':>7} "
              f"{'soil evap':>10} {'transp':>8}")
        for m in range(len(months)):
            frozen = "yes" if t[sel, m].mean() < 0 else ""
            print(f"  {months[m]:>6} {t[sel, m].mean():9.2f} "
                  f"{p[sel, m].mean():10.3f} {frozen:>7} "
                  f"{e[sel, m].mean():10.3f} {x[sel, m].mean():8.3f}")

    for label, threshold in (("above 0 degC", LOOSE_C), ("above 5 degC", STRICT_C)):
        warm = t > threshold
        share = lambda a: float(np.nanmean(np.nansum(a * warm, axis=1)
                                           / np.nansum(a, axis=1)))
        print(f"\n== months {label}: {warm.sum(axis=1).mean():.2f} of "
              f"{len(months)} ==")
        print(f"  bare-soil evaporation inside the window : {share(e) * 100:.1f}%")
        print(f"  precipitation delivered inside it       : {share(p) * 100:.1f}%")
        print(f"  transpiration inside it                 : {share(x) * 100:.1f}%")
        if threshold == LOOSE_C:
            verdict = ("the bare-soil evaporation term STANDS" if share(e) >= STANDS
                       else "the term is LARGELY FICTIONAL" if share(e) <= FICTIONAL
                       else "PARTIAL")
            print(f"  -> {verdict} "
                  f"(declared rule: >={STANDS:.0%} stands, <={FICTIONAL:.0%} not)")

    print("\n== the melt against leaf-out ==")
    above0, above5 = t > LOOSE_C, t > STRICT_C
    first0 = np.where(above0.any(axis=1), above0.argmax(axis=1), np.nan)
    first5 = np.where(above5.any(axis=1), above5.argmax(axis=1), np.nan)
    both = np.isfinite(first0) & np.isfinite(first5)
    pack = np.where(~above0, p, 0.0).sum(axis=1)
    delivered_warm = np.nansum(p * above0, axis=1)
    print(f"  cells reaching 0 / 5 degC: {int(np.isfinite(first0).sum())} / "
          f"{int(np.isfinite(first5).sum())} of {polar.sum()}")
    print(f"  thaw to growth, paired over the {int(both.sum())} reaching both: "
          f"{np.mean(first5[both] - first0[both]):.3f} months = "
          f"{np.mean(first5[both] - first0[both]) * days[0]:.2f} days")
    index = np.arange(len(months))[None, :]
    before = float(np.nanmean(np.nansum(e * (index < first0[:, None]), axis=1)
                              / np.nansum(e, axis=1)))
    print(f"  bare-soil evaporation strictly before the thaw month: {before * 100:.1f}%")
    print(f"  precipitation accumulated below freezing: {pack.mean():.2f} mm of "
          f"{p.sum(axis=1).mean():.2f} ({pack.mean() / p.sum(axis=1).mean() * 100:.1f}%)")
    print(f"  water crossing the surface inside the warm window: "
          f"{(pack + delivered_warm).mean():.2f} mm over "
          f"{above0.sum(axis=1).mean():.2f} months "
          f"({pack.mean():.2f} as melt, {delivered_warm.mean():.2f} delivered warm)")

    print("\n== how concentrated the flux is ==")
    for name, table in (("bare-soil evaporation", e), ("precipitation", p),
                        ("transpiration", x)):
        total = np.nansum(table, axis=1)
        peak = np.nanmax(table, axis=1)
        held = total > 0
        print(f"  {name:<22}: {np.mean(peak[held] / total[held]) * 100:5.1f}% of "
              f"the annual total in its single largest month")
    print("\n  The drainage term's month is NOT measured: tot_runoff.out is "
          "annual. It is inferred to sit in the same month as the surface flux.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
