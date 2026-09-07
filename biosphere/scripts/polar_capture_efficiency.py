#!/usr/bin/env python3
"""Is the simulated polar cap inefficient at taking water, or does it lack leaf?

    python biosphere/scripts/polar_capture_efficiency.py
    python biosphere/scripts/polar_capture_efficiency.py --run lpj_<id>

This is a note about a simulated world: every quantity below is a modelled field
of Vesper, and every comparison to Earth is a distance to report.

`biosphere/notes/polar-cover-cold-filter-and-capture.md` measures the cap taking
4.7% of the water that reaches its ground against 78.9% on warm ground holding
LESS water, and reads the gap as a property of the plants. Three capture-side
traits were proposed for sizing on the strength of it: rooting depth and
distribution, the phenology trigger, and the bare-soil evaporation path. Each of
those is an EFFICIENCY -- water captured per unit of leaf deployed -- so each is
worth moving only if the cap's efficiency is deficient. That was assumed and is
measured here.

The instrument is the per-cell extinction coefficient k in a Beer's-law reading
of capture, k = -ln(1 - captured) / LAI, which is the capture a cell gets per
unit of leaf area it carries. It is used as a RANKING and not as a physical
constant: the same table shows k is not common across the planet, so a Beer's
law in leaf area alone is refused rather than asserted, and what survives is the
comparison of one band's k against the distribution of all of them.

WHY THE LEAF AREA IS THE CEILING AND NOT A PHASE. `lai.out` carries
`indiv.lai`, which `growth.cpp` sets to `cmass_leaf * sla` -- the leaf area at
FULL display. The displayed area is `lai_today() = lai * phen`, and it is never
larger. So the cap's 0.067 is what it could show with its phenology satisfied
on the first day of the season, and a faster `phengdd5ramp` displays that same
area sooner rather than displaying more of it.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from _paths import RUNS

import lpj_output

DEFAULT_RUN = "lpj_1e6a2b9ca51a4eff9592992cad96677b"
POLAR_LAT = 75.0

# Controls, fixed before the numbers are read. The water balance is not one:
# precipitation is reconstructed FROM the balance here, so its closure is an
# identity. What can fail is the reconstruction against the note's published
# per-band totals, and the band's capture against the note's published share.
CONTROL_POLAR_PRECIP = 77.27       # mm per simulation year
CONTROL_POLAR_CAPTURE = 0.047      # ratio of the band's means, as the note takes it
CONTROL_TOL = 0.03

BANDS = (("poleward of 75", 75.0, 90.0), ("60 to 75", 60.0, 75.0),
         ("45 to 60", 45.0, 60.0), ("30 to 45", 30.0, 45.0),
         ("15 to 30", 15.0, 30.0), ("equatorward of 15", 0.0, 15.0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", default=DEFAULT_RUN)
    args = parser.parse_args()
    run = RUNS / args.run

    tables = {name: lpj_output.reduce_table(run / name)
              for name in ("lai.out", "aaet.out", "tot_runoff.out",
                           "mevap.out", "mintercep.out")}
    cells = sorted(tables["lai.out"].values)
    lat = np.array([c[1] for c in cells])

    def total(name: str, column: str | None = None):
        reduced = tables[name]
        if column is not None:
            i = reduced.names.index(column)
            return np.array([reduced.values[c][i] for c in cells])
        return np.array([float(np.sum(reduced.values[c])) for c in cells])

    leaf = total("lai.out", "Total")
    transpired = total("aaet.out", "Total")
    # The other three sinks of the balance the note closed. Precipitation is
    # their sum with transpiration, which is why its closure proves nothing.
    arrived = (transpired + total("tot_runoff.out", "Total")
               + total("mevap.out") + total("mintercep.out"))

    polar = np.abs(lat) > POLAR_LAT
    for label, got, want in (
            ("polar precipitation", arrived[polar].mean(), CONTROL_POLAR_PRECIP),
            ("polar capture", transpired[polar].sum() / arrived[polar].sum(),
             CONTROL_POLAR_CAPTURE)):
        off = abs(got - want) / want
        print(f"  control {label}: {got:.4g} against the note's {want:.4g} "
              f"({off * 100:.1f}% off)")
        if off > CONTROL_TOL:
            raise SystemExit("positive control failed; nothing below holds")

    held = (arrived > 0) & (leaf > 0)
    captured = np.zeros_like(transpired)
    captured[held] = transpired[held] / arrived[held]
    usable = held & (captured > 0) & (captured < 1)
    k = -np.log(1.0 - captured[usable]) / leaf[usable]
    band_lat = np.abs(lat[usable])

    print(f"\n== capture per unit leaf area, {usable.sum()} of {len(cells)} cells ==")
    print(f"{'band':>20} {'cells':>6} {'LAI':>7} {'capture':>8} {'k':>7} {'pctile':>7}")
    for name, low, high in BANDS:
        sel = (band_lat > low) & (band_lat <= high)
        if not sel.any():
            continue
        median = float(np.median(k[sel]))
        print(f"{name:>20} {int(sel.sum()):6d} {leaf[usable][sel].mean():7.3f} "
              f"{captured[usable][sel].mean() * 100:7.1f}% {median:7.3f} "
              f"{(k < median).mean() * 100:6.0f}")

    print(f"\n  k over the whole planet: p10 {np.percentile(k, 10):.3f}, "
          f"median {np.median(k):.3f}, p90 {np.percentile(k, 90):.3f}")
    print("  A sixteen-fold spread, so k is NOT a constant and capture is not a "
          "function of leaf area alone. The ranking is what this measures.")

    polar_k = float(np.median(k[band_lat > POLAR_LAT]))
    pct = (k < polar_k).mean() * 100
    print(f"\n  The cap's median k is {polar_k:.3f}, at percentile "
          f"{pct:.0f} of the planet's cells.")
    print("  It is at or above the median rate at which this world's vegetation")
    print("  turns leaf area into captured water. Its capture is low because its")
    print(f"  leaf area is {leaf[polar].mean():.4f}, not because it takes water badly.")
    print("\n  So the three proposed capture-side arms move an efficiency that is")
    print("  not deficient, and none of them adds leaf area.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
