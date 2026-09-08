"""Two conservation bounds on the simulated fire nitrogen flux.

Worldbuilding. Vesper is an invented planet; this checks the nitrogen its
simulated vegetation model loses to fire against two quantities the same model
already reports, and nothing here is about fire on Earth.

WHY THIS EXISTS. `assess_lpj_run.py` checks that values are finite, physical and
closing, and asks nothing about whether a flux is ORDINARY. A large fire
emission is finite, is positive, and CLOSES -- the nitrogen it removes is
nitrogen the model had -- so no amount of closure testing can see it. FIRE-8 is
the row that adds "reject impossible ranges" and this is the part of it that can
be settled without first choosing a fire model.

THE TWO BOUNDS, and why each is a test rather than a description. CLAUDE.md's
rule is that a check needs a right answer: an identity, a conservation law, or a
quantity the other side already knows. Both of these are the second.

  STOCK   In one simulated year a gridcell cannot lose more nitrogen to fire
          than the nitrogen it holds. Fire moves nitrogen out of the vegetation,
          litter and soil pools, so a year whose fire flux exceeds the cell's
          total nitrogen pool is not a large fire, it is a defect.

  RATE    Over the retained record a gridcell is at equilibrium, so its nitrogen
          losses equal its nitrogen inputs. Fire is one of at least three loss
          pathways -- fire gases, soil gases, leaching -- so the retained-record
          mean fire loss must be strictly below the retained-record mean input.
          A cell above it is either not at equilibrium or is losing nitrogen it
          never received.

BOTH ARE DERIVED FROM THE INPUT SIDE, which is deliberate and is the reason they
may be applied to a run whose fire distribution has already been looked at.
CLAUDE.md: a criterion chosen after the run it judges is not a criterion. Neither
bound can be reached by the distribution it judges, because neither is written in
terms of it -- one is the cell's own stock and the other is the cell's own
supply.

THE UNIT TRAP THIS MODULE EXISTS TO GET RIGHT. `nflux.out` and `ngases.out` are
multiplied by `M2_PER_HA` at `commonoutput.cpp:1983` and `:2085` and are in
kgN/ha. `npool.out` at `:2057` is NOT, and is in kgN/m2. The comment at `:196`
reading "kgN/m2, converted by 1e4 to kgN/ha" is about the PRECISION the column is
given, not about the column being converted. Comparing the two without the factor
of 10^4 fails the stock bound on every gridcell-year by about three orders of
magnitude and looks entirely tidy doing it, which is the case CLAUDE.md's
check-the-instrument rule is about. `M2_PER_HA` below is read from the model's
own source rather than written here, so a fork that changed it moves this too.

    python biosphere/scripts/fire_nitrogen_bounds.py <run-dir>
    python biosphere/scripts/fire_nitrogen_bounds.py <run-dir> --strict

`--strict` exits non-zero when either bound is violated. Without it the module
reports and exits 0, because a violation is a finding about a RUN rather than
about the tree, and this is registered as a one-off rather than as a gate.

`biosphere/notes/fire-nitrogen-range.md` carries the first application of it and
what it settled.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from _paths import PROJECT_ROOT
from lpj_table import RowParseRequired, read as read_table  # noqa: E402

# The four fire species, in the order `commonoutput.cpp` writes them. Soil
# columns sit between them and are deliberately not summed: this bound is about
# the FIRE pathway, and including the soil operator's gases would let a soil
# defect be absorbed by a fire headroom or the reverse.
FIRE_SPECIES = ("NH3_fire", "NOx_fire", "N2O_fire", "N2_fire")

# `nflux.out` writes every supply term NEGATED (`-aNH4dep_gridcell` and its
# neighbours at commonoutput.cpp:1983-1986), so a positive supply is the
# negated sum. `fert` is included because it is a supply when it is nonzero and
# silently zero when it is not; leaving it out would understate the input on any
# run that ever turns land use on.
INPUT_FIELDS = ("NH4dep", "NO3dep", "fix", "fert")

# The nitrogen a cell holds, from `npool.out`.
POOL_FIELD = "Total"

GUESS_SOURCE = PROJECT_ROOT / "vendor" / "lpj-guess"


def m2_per_ha(root: Path = GUESS_SOURCE) -> float:
    """The model's own m2-per-hectare constant, read rather than written.

    The whole correctness of this module is one factor of 10^4 between two
    tables, so the factor is taken from the source that applies it. A fork that
    redefined it would move this check with it instead of leaving it quietly
    wrong.
    """
    for candidate in (root / "framework" / "guessmath.h",
                      root / "framework" / "guess.h"):
        if not candidate.is_file():
            continue
        text = candidate.read_text(encoding="utf-8", errors="replace")
        text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
        text = re.sub(r"//[^\n]*", " ", text)
        found = re.search(r"\bM2_PER_HA\s*=\s*([0-9.eE+-]+)", text)
        if found:
            return float(found.group(1))
    raise SystemExit(
        "M2_PER_HA is not declared anywhere this module knows to look under "
        f"{root}. Refusing rather than assuming 1e4: this check is one unit "
        "conversion and a guessed one would pass on every run.")


def _column(table, name: str) -> np.ndarray:
    if name not in table.fields:
        raise SystemExit(f"the table has no column {name!r}; it has "
                         f"{table.fields}")
    return table.values[:, table.fields.index(name)]


def _load(path: Path):
    try:
        return read_table(path)
    except RowParseRequired as exc:
        raise SystemExit(
            f"{path.name} could not be certified column-wise: {exc}. That is "
            f"lib/lpj_table.py's certify-or-decline contract declining, and a "
            f"malformed table is a finding in its own right -- see world-v5xp "
            f"for the shape it takes.") from exc


def bounds(run: Path) -> dict:
    """Both bounds over every gridcell-year of one run."""
    ngases = _load(run / "ngases.out")
    npool = _load(run / "npool.out")
    nflux = _load(run / "nflux.out")

    for other, name in ((npool, "npool.out"), (nflux, "nflux.out")):
        if other.rows != ngases.rows or not np.array_equal(other.key, ngases.key):
            raise SystemExit(
                f"{name} and ngases.out do not cover the same cell-years, so "
                f"no row-wise comparison between them is meaningful. Both are "
                f"key-ordered by lib/lpj_table.py, so this is a difference in "
                f"CONTENT and not in order.")

    factor = m2_per_ha()
    fire = sum(_column(ngases, s) for s in FIRE_SPECIES)      # kgN/ha/yr
    pool = _column(npool, POOL_FIELD) * factor                # kgN/m2 -> kgN/ha
    supply = -sum(_column(nflux, f) for f in INPUT_FIELDS)    # kgN/ha/yr

    # STOCK, per gridcell-year.
    positive = pool > 0.0
    stock_ratio = np.zeros_like(fire)
    stock_ratio[positive] = fire[positive] / pool[positive]
    stock_violations = int(np.count_nonzero(stock_ratio > 1.0))
    fire_without_pool = int(np.count_nonzero((~positive) & (fire > 0.0)))
    worst_stock = int(np.argmax(stock_ratio)) if fire.size else 0

    # RATE, per gridcell over its own retained record.
    starts = ngases.cell_starts()
    edges = np.append(starts, ngases.rows)
    counts = np.diff(edges)
    mean_fire = np.add.reduceat(fire, starts) / counts
    mean_supply = np.add.reduceat(supply, starts) / counts
    supplied = mean_supply > 0.0
    rate_ratio = np.full(mean_fire.shape, np.nan)
    rate_ratio[supplied] = mean_fire[supplied] / mean_supply[supplied]
    rate_violations = int(np.count_nonzero(rate_ratio > 1.0))
    worst_rate = int(np.nanargmax(rate_ratio)) if supplied.any() else 0
    cells = ngases.cells()

    return {
        "run": run.name,
        "m2_per_ha": factor,
        "cell_years": ngases.rows,
        "cells": len(cells),
        "years": len(ngases.years()),
        "stock": {
            "violations": stock_violations,
            "cells_with_fire_and_no_pool": fire_without_pool,
            "worst_ratio": float(stock_ratio[worst_stock]),
            "worst_at": {"lon": float(ngases.lon[worst_stock]),
                         "lat": float(ngases.lat[worst_stock]),
                         "year": int(ngases.year[worst_stock])},
            "worst_fire_kgN_ha": float(fire[worst_stock]),
            "worst_pool_kgN_ha": float(pool[worst_stock]),
            "max_fire_kgN_ha": float(fire.max()) if fire.size else 0.0,
        },
        "rate": {
            "violations": rate_violations,
            "cells_with_no_supply": int(np.count_nonzero(~supplied)),
            "mean_ratio": float(np.nanmean(rate_ratio)),
            "worst_ratio": float(rate_ratio[worst_rate]),
            "worst_at": {"lon": cells[worst_rate][0],
                         "lat": cells[worst_rate][1]},
            "worst_fire_kgN_ha_yr": float(mean_fire[worst_rate]),
            "worst_supply_kgN_ha_yr": float(mean_supply[worst_rate]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="an LPJ-GUESS run directory")
    parser.add_argument("--strict", action="store_true",
                        help="exit non-zero when either bound is violated")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    if not args.run.is_dir():
        raise SystemExit(f"{args.run} is not a directory")

    report = bounds(args.run)
    report["measured"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        s, r = report["stock"], report["rate"]
        print(f"Fire nitrogen against what the simulated cells hold and "
              f"receive.\n")
        print(f"  run          {report['run']}")
        print(f"  covered      {report['cell_years']:,} cell-years "
              f"({report['cells']} cells x {report['years']} years)")
        print(f"  kgN/m2 -> kgN/ha factor read from the model source: "
              f"{report['m2_per_ha']:g}\n")
        print(f"  STOCK  fire nitrogen against the nitrogen the cell holds")
        print(f"    cell-years where fire > pool   {s['violations']}")
        print(f"    worst ratio                    {s['worst_ratio']:.4f}  "
              f"at lon {s['worst_at']['lon']}, lat {s['worst_at']['lat']}, "
              f"year {s['worst_at']['year']}")
        print(f"      fire {s['worst_fire_kgN_ha']:.3f} kgN/ha from a pool of "
              f"{s['worst_pool_kgN_ha']:.1f} kgN/ha")
        if s["cells_with_fire_and_no_pool"]:
            print(f"    fire with no pool at all       "
                  f"{s['cells_with_fire_and_no_pool']}  <- impossible")
        print(f"\n  RATE   mean fire loss against mean nitrogen supply")
        print(f"    cells where fire > supply      {r['violations']} of "
              f"{report['cells']}")
        print(f"    mean ratio                     {r['mean_ratio']:.4f}")
        print(f"    worst ratio                    {r['worst_ratio']:.4f}  "
              f"at lon {r['worst_at']['lon']}, lat {r['worst_at']['lat']}")
        print(f"      fire {r['worst_fire_kgN_ha_yr']:.3f} against a supply of "
              f"{r['worst_supply_kgN_ha_yr']:.3f} kgN/ha/yr")
        violated = s["violations"] or r["violations"] or s["cells_with_fire_and_no_pool"]
        if violated:
            print(f"\n  A BOUND IS VIOLATED. Neither bound involves a fire "
                  f"model, a burned area or an Earth comparison, so a "
                  f"violation is a defect in the operator rather than a "
                  f"disagreement about fire regime.")
        else:
            print(f"\n  Both bounds hold. That does NOT say the fire regime is "
                  f"right: a model with too few, too large fires satisfies "
                  f"both exactly as well as a correct one, and return interval "
                  f"and burned fraction belong to fire-9.")

    if args.strict:
        s, r = report["stock"], report["rate"]
        if s["violations"] or r["violations"] or s["cells_with_fire_and_no_pool"]:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
