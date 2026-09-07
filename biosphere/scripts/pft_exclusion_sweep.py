#!/usr/bin/env python3
"""Where do the shipped plant functional types leave this world's ground empty?

    python biosphere/scripts/pft_exclusion_sweep.py
    python biosphere/scripts/pft_exclusion_sweep.py --run lpj_<id>

This is a note about a simulated world: every quantity below is a modelled field
of Vesper, and every comparison to Earth is a distance to report.

`world-orok` found one niche the twelve types do not reach, the polar cap, by
looking at it. This asks the question over the whole planet instead: for every
gridcell, which types are barred and by which declared limit, and is the cover
that remains running inside its own declared photosynthetic band.

TWO SIGNATURES, KEPT APART. A type barred from a cell is not evidence of
anything on its own: boreal types are correctly absent from the tropics and
tropical types from the caps, and the ground is full either way. What the polar
diagnosis actually rested on is the SECOND signature -- the cover that remains
runs outside the temperature band it declares -- and that is what generalises.
The count of barred types is reported beside it as context, never as the finding.

SURVIVAL AND ESTABLISHMENT ARE NOT THE SAME BAR. `tcmin_surv` kills, so a type
below it must carry exactly zero cover and that is this script's control.
`tcmin_est`, `tcmax_est`, `twmin_est` and `gdd5min_est` bar RECRUITMENT only, so
a type that established under conditions since departed persists and legitimately
carries cover under a violated limit. Conflating the two reported sixteen sound
cells as broken.

THE CALENDAR TOLERANCE. The driver's intervals are 15.083 days and LPJ-GUESS's
months are whole days, so the two cannot align: a day straddling two intervals
takes their duration-weighted mean and the model's coldest monthly mean comes
out slightly warmer than the driver's raw interval minimum. Every exception seen
was inside 0.31 degC and in that direction. The control therefore carries a
declared 0.5 degC margin, which a genuine misreading of a limit would miss by
far more.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

from _paths import COMPONENT_ROOT, GENERATED, RUNS

import lpj_output
import lpj_pfts

DEFAULT_RUN = "lpj_1e6a2b9ca51a4eff9592992cad96677b"
DRIVER = GENERATED / "vesper_driver.bin"
GDD_BASE = 5.0
CAL_TOL = 0.5

SURVIVAL = "tcmin_surv"
ESTABLISHMENT = ("tcmin_est", "tcmax_est", "twmin_est", "gdd5min_est")

BANDS = (("poleward of 75", 75.0, 90.0), ("60 to 75", 60.0, 75.0),
         ("45 to 60", 45.0, 60.0), ("30 to 45", 30.0, 45.0),
         ("15 to 30", 15.0, 30.0), ("equatorward of 15", 0.0, 15.0))

# The note's polar figures, as controls on the forcing read.
CONTROL = {"coldest": -68.75, "warmest": 31.11}
CONTROL_TOL = 0.02


def read_driver():
    """Borrow the reader that already states the driver's layout."""
    spec = importlib.util.spec_from_file_location(
        "polar_water_timing", COMPONENT_ROOT / "scripts" / "polar_water_timing.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_driver(DRIVER)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=DEFAULT_RUN)
    args = parser.parse_args()
    run = RUNS / args.run

    lon, lat, tas, pr_rate, days, year_days, _layers = read_driver()
    coldest, warmest = tas.min(axis=1), tas.max(axis=1)
    gdd5 = (np.maximum(0.0, tas - GDD_BASE) * days[None, :]).sum(axis=1)
    precip = (pr_rate * days[None, :]).sum(axis=1)
    polar = np.abs(lat) > 75.0

    for label, got in (("coldest", coldest[polar].mean()),
                       ("warmest", warmest[polar].mean())):
        want = CONTROL[label]
        off = abs(got - want) / abs(want)
        print(f"  control polar {label} month: {got:.2f} against the note's "
              f"{want:.2f} ({off * 100:.1f}% off)")
        if off > CONTROL_TOL:
            raise SystemExit("forcing control failed; nothing below holds")

    names = [n for n in lpj_pfts.names() if n != "Total"]
    declared = {n: lpj_pfts.parameters(n) for n in names}

    fpc = lpj_output.reduce_table(run / "fpc.out")
    keys = [(round(float(a), 2), round(float(b), 2)) for a, b in zip(lon, lat)]
    cover = np.array([[fpc.values[k][fpc.names.index(n)] for n in names]
                      for k in keys])

    killed = np.zeros(cover.shape, dtype=bool)
    barred = np.zeros(cover.shape, dtype=bool)
    for j, name in enumerate(names):
        p = declared[name]
        killed[:, j] = coldest < p.get(SURVIVAL, -1e3) - CAL_TOL
        barred[:, j] = ((coldest < p.get("tcmin_est", -1e3))
                        | (coldest > p.get("tcmax_est", 1e3))
                        | (warmest < p.get("twmin_est", -1e3))
                        | (gdd5 < p.get("gdd5min_est", 0.0)))

    violations = int((killed & (cover > 0)).sum())
    print(f"  control survival: types below {SURVIVAL} carrying cover: "
          f"{violations} (must be 0, margin {CAL_TOL} degC)")
    if violations:
        raise SystemExit("survival control failed")
    print(f"  types barred from ESTABLISHING yet carrying cover: "
          f"{int((barred & (cover > 0)).sum())} -- legitimate, they persist")

    # The signature that generalises: is the resident cover inside its own band?
    weights = np.broadcast_to(days[None, :], tas.shape)
    season = tas > GDD_BASE
    season_time = np.where(season, weights, 0.0).sum(axis=1)
    below = np.zeros(len(keys))
    above = np.zeros(len(keys))
    stopped = np.zeros(len(keys))
    held = np.zeros(len(keys))
    for j, name in enumerate(names):
        p = declared[name]
        low, high, top = (p.get("pstemp_low"), p.get("pstemp_high"),
                          p.get("pstemp_max"))
        if low is None or high is None or top is None:
            continue
        f = cover[:, j]
        def share(mask):
            return np.divide(np.where(season & mask, weights, 0.0).sum(axis=1),
                             season_time, out=np.zeros(len(keys)),
                             where=season_time > 0)
        below += f * share(tas < low)
        above += f * share(tas > high)
        stopped += f * share(tas > top)
        held += f
    ok = held > 0
    for arr in (below, above, stopped):
        arr[ok] /= held[ok]

    total = cover.sum(axis=1)
    print(f"\n== by band: what is barred, and whether what remains fits ==")
    print(f"{'band':>20} {'cells':>6} {'killed':>7} {'barred':>7} {'GDD5':>8} "
          f"{'precip':>7} {'FPC':>6} {'>high':>7} {'>max':>6}")
    for name, low, high in BANDS:
        sel = ok & (np.abs(lat) > low) & (np.abs(lat) <= high)
        print(f"{name:>20} {int(sel.sum()):6d} {killed[sel].sum(1).mean():7.2f} "
              f"{barred[sel].sum(1).mean():7.2f} {gdd5[sel].mean():8.1f} "
              f"{precip[sel].mean():7.1f} {total[sel].mean():6.3f} "
              f"{above[sel].mean() * 100:6.1f}% {stopped[sel].mean() * 100:5.1f}%")
    print(f"\n  '>high' and '>max' are the cover-weighted share of growing-season")
    print(f"  time the RESIDENT types spend above their own pstemp_high, where")
    print(f"  photosynthesis is already declining, and above pstemp_max, where it")
    print(f"  stops. A band with cover fitting its band reads near zero on both.")

    # WHO HOLDS THE GROUND THAT IS OCCUPIED. The exclusion table says what is
    # barred; this says what the survivors build, which is the other half of
    # whether a niche is underoccupied or merely differently occupied.
    grass = [names.index(n) for n in lpj_pfts.grass() if n in names]
    trees = [names.index(n) for n in lpj_pfts.trees() if n in names]
    g = cover[:, grass].sum(axis=1)
    t = cover[:, trees].sum(axis=1)
    vegetated = total > 0.05
    dominated = vegetated & (g > t)
    print(f"\n== who holds the ground ==")
    print(f"  grass carries {g.sum() / total.sum() * 100:.1f}% of the planet's cover; "
          f"{int(dominated.sum())} cells are grass-dominated against "
          f"{int((vegetated & (t >= g)).sum())} tree-dominated")
    print(f"{'band':>20} {'grass-dom':>10} {'% of band':>10} {'grass FPC':>10} "
          f"{'tree FPC':>9}")
    for name, low, high in BANDS:
        inband = (np.abs(lat) > low) & (np.abs(lat) <= high)
        sel = dominated & inband
        if not inband.any():
            continue
        print(f"{name:>20} {int(sel.sum()):10d} "
              f"{sel.sum() / inband.sum() * 100:9.1f}% "
              f"{g[sel].mean() if sel.any() else float('nan'):10.3f} "
              f"{t[sel].mean() if sel.any() else float('nan'):9.3f}")
    print(f"\n  total cover over all land cells: p25 {np.percentile(total, 25):.3f}, "
          f"median {np.percentile(total, 50):.3f}, p75 {np.percentile(total, 75):.3f}")

    hurt = ok & (stopped > 0.05)
    print(f"\n== cells losing more than 5% of the season to pstemp_max shutdown ==")
    print(f"  {int(hurt.sum())} cells, mean |lat| {np.abs(lat[hurt]).mean():.1f}, "
          f"FPC {total[hurt].mean():.3f}, warmest month {warmest[hurt].mean():.1f} degC")
    for name, low, high in BANDS:
        sel = hurt & (np.abs(lat) > low) & (np.abs(lat) <= high)
        if sel.any():
            print(f"    {name:>20}: {int(sel.sum()):4d} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())
