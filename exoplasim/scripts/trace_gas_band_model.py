#!/usr/bin/env python3
"""The CH4 and N2O band model, checked against the tables it comes from.

    python exoplasim/scripts/trace_gas_band_model.py
    python exoplasim/scripts/trace_gas_band_model.py --fortran

Worldbuilding. Vesper is an invented planet and this script is about the
simulation of it: the longwave band the climate model does not yet carry, and
the laboratory parameters that define it. The gas measurements are real; what
is invented is the atmosphere they will be applied to.

CLIM-42. `exoplasim/notes/trace-gas-band.md` is the design and the declared
tests; this is the first of those tests, and it is the one that can be run
without the model. It writes `exoplasim/analysis/trace_gas_band_model.json` and
`--fortran` prints the parameter block for `radmod.f90` so the constants reach
the model from here rather than being retyped.

## The band model

Cess and Ramanathan (1972) as modified by Ramanathan (1976), in the form Donner
and Ramanathan (1980) use it:

    A(U, beta) = 2 A0 ln[ 1 + U / sqrt(4 + U (1 + 1/beta)) ]     (1)
    U    = S W / A0                                              (2)
    beta = beta0 (P / P0)                                        (3)

`A` is the total band absorptance in cm-1, `W` the absorber amount in cm atm at
STP, `P` the broadening pressure against `P0` = 1 atm. `A0`, `beta0` and `S` are
the bandwidth, line shape and band intensity parameters.

## The three parameters are a MATCHED TRIPLE, and that decides the sourcing

`A0` and `beta0` are not independent measurements: Donner and Ramanathan
obtained them by fitting Eq. (1), at a particular `S`, to laboratory
absorptance. Substituting a different compilation's `S` into their `A0` and
`beta0` therefore breaks the fit rather than improving it, however modern the
substitute. So each band's `S` is taken from the source THAT PAPER used:

- **CH4 1306 cm-1** -- Cess and Chen, which the paper does not restate. It is
  recovered instead from the paper's own Table 2, eleven (P, W, absorptance)
  triples for this band: with `A0` and `beta0` fixed by Table 1, Eq. (1) has one
  free parameter and the table determines it.
- **N2O 1285 and 589 cm-1** -- McClatchey et al. (1973) Table 13, the band
  SYSTEM intensities, which is what the paper names. A system intensity is the
  right quantity: the paper's 1285 cm-1 analysis "includes the fundamental and
  the first hot band and, furthermore, accounts for contribution from four
  isotopes", and Table 13 sums a system rather than a single transition.

McClatchey quotes intensities per molecule cm-2; Donner's `W` is cm atm at STP,
so the conversion is Loschmidt's number and nothing else.

## What the checks can fail

**Table 2, and it is an identity against printed numbers.** Eq. (1) with the
recovered CH4 `S` must reproduce the paper's own "Eq. (1)" column. It does, on
ten of the eleven rows, to 1.6% relative. The eleventh is excluded as a
typesetting slip and the exclusion is declared here rather than discovered:
its entry reprints the number standing diagonally below it in the adjacent
column, and its neighbours fit to a fifth of a wavenumber.

**The unit conversion, against a band whose intensity is in BOTH sources.**
McClatchey Table 18 carries CH4's 1306 cm-1 intensity as well, so the same
arithmetic that converts the N2O numbers can be run on a band where an
independent answer already exists. It gives 158 against the 188 recovered from
Donner's Table 2 -- a 19% difference, which is the gap between two published
compilations of one band and NOT an arithmetic error, because any unit slip
here would be a factor of 2.7e19, 10 or 100 rather than 1.19. That is the whole
strength of this check and it is stated rather than overclaimed: it rules out a
unit error and it does not validate the value.

**Continuity and limits.** `A` must be zero at zero absorber, must rise
monotonically in `W`, and must approach the weak-line limit `A -> S W` as the
amount goes to zero, where the line shape drops out. That last one is what
makes `S` a band intensity rather than a fitted constant.
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from math import log, sqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402

OUTPUT = ANALYSIS / "trace_gas_band_model.json"

# Loschmidt's number: molecules per cm3 at STP, 273.15 K and 1 atm. This is the
# whole of the conversion from McClatchey's per-molecule intensities to Donner's
# per-cm-atm ones.
LOSCHMIDT = 2.6868e19

# Donner and Ramanathan (1980) Table 1. `A0` in cm-1 scales as (T/300)^0.5 and
# `beta0` as (300/T)^0.5; the values here are at 300 K.
BANDS = {
    "ch4_1306": dict(centre_cm1=1306.0, a0_300k=52.0, beta0_300k=0.17),
    "n2o_1285": dict(centre_cm1=1285.0, a0_300k=20.4, beta0_300k=1.12),
    "n2o_589": dict(centre_cm1=589.0, a0_300k=23.0, beta0_300k=1.08),
}

# McClatchey et al. (1973) Table 13, N2O band SYSTEM intensities, in units of
# 1e-20 cm-1 / (molecule cm-2). Table 18 carries the CH4 entry, which is used
# only as a check on the conversion and NOT as the model's CH4 intensity.
MCCLATCHEY = {
    "n2o_1285": 996e-20,
    "n2o_589": 118e-20,
    "ch4_1306_crosscheck": 5.87e-18,
}

# Donner and Ramanathan (1980) Table 2: total pressure atm, absorber amount
# cm atm STP, and the paper's own Eq. (1) absorptance in cm-1.
TABLE2 = [(1.0, 0.916, 51.6), (1.0, 0.505, 41.6), (1.0, 0.458, 36.5),
          (1.0, 0.231, 24.4), (0.3, 0.275, 19.1), (0.3, 0.137, 12.8),
          (0.3, 0.069, 8.2), (0.1, 0.092, 6.9), (0.1, 0.046, 4.5),
          (0.1, 0.023, 2.8), (0.03, 0.028, 2.1)]

# Declared BEFORE fitting: this row's Eq. (1) entry of 41.6 is the number
# printed diagonally below it in the adjacent column, and its neighbours fit to
# a fifth of a wavenumber. Excluded as a typesetting slip; the fit is reported
# both ways so the exclusion cannot hide behind the result.
TABLE2_SUSPECT = (1.0, 0.505, 41.6)

TOLERANCE_PCT = 2.0        # declared bar for the Table 2 identity


def absorptance(amount: float, pressure: float, a0: float, beta0: float,
                intensity: float) -> float:
    """Eq. (1): total band absorptance in cm-1."""
    if amount <= 0.0:
        return 0.0
    u = intensity * amount / a0
    beta = beta0 * pressure
    return 2.0 * a0 * log(1.0 + u / sqrt(4.0 + u * (1.0 + 1.0 / beta)))


def scaled(band: str, temperature: float) -> tuple[float, float]:
    """`A0` and `beta0` at a temperature, per Table 1's scaling."""
    b = BANDS[band]
    return (b["a0_300k"] * sqrt(temperature / 300.0),
            b["beta0_300k"] * sqrt(300.0 / temperature))


def fit_intensity(rows, a0: float, beta0: float) -> tuple[float, float, float]:
    """The one `S` that reproduces a set of (P, W, A) rows. Bisection on the
    residual's derivative is unnecessary: `A` rises monotonically in `S`, so a
    golden-section search on the rms is exact to the printed precision."""
    lo, hi = 1.0, 2000.0
    for _ in range(200):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        r1 = sum((absorptance(w, p, a0, beta0, m1) - a) ** 2 for p, w, a in rows)
        r2 = sum((absorptance(w, p, a0, beta0, m2) - a) ** 2 for p, w, a in rows)
        if r1 < r2:
            hi = m2
        else:
            lo = m1
    s = 0.5 * (lo + hi)
    diffs = [absorptance(w, p, a0, beta0, s) - a for p, w, a in rows]
    rms = sqrt(sum(d * d for d in diffs) / len(diffs))
    worst = max(abs(d) / a for (_, _, a), d in zip(rows, diffs))
    return s, rms, worst * 100.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fortran", action="store_true",
                    help="print the radmod.f90 parameter block and exit")
    args = ap.parse_args()

    a0, beta0 = BANDS["ch4_1306"]["a0_300k"], BANDS["ch4_1306"]["beta0_300k"]
    kept = [r for r in TABLE2 if r != TABLE2_SUSPECT]
    s_ch4, rms, worst_pct = fit_intensity(kept, a0, beta0)
    s_all, rms_all, _ = fit_intensity(TABLE2, a0, beta0)

    if worst_pct > TOLERANCE_PCT:
        raise SystemExit(f"the CH4 band intensity reproduces Donner Table 2 to "
                         f"only {worst_pct:.1f}%, against a declared {TOLERANCE_PCT}%")

    intensities = {"ch4_1306": s_ch4}
    for band in ("n2o_1285", "n2o_589"):
        intensities[band] = MCCLATCHEY[band] * LOSCHMIDT

    # The conversion check, on the one band both sources carry.
    s_ch4_mcc = MCCLATCHEY["ch4_1306_crosscheck"] * LOSCHMIDT
    ratio = s_ch4 / s_ch4_mcc
    if not 0.5 < ratio < 2.0:
        raise SystemExit(
            f"the per-molecule conversion puts McClatchey's CH4 intensity at "
            f"{s_ch4_mcc:.1f} against {s_ch4:.1f} recovered from Donner Table 2, "
            f"a factor of {ratio:.3g}. Two compilations of one band differ by "
            "tens of percent; a factor is a unit error.")

    # Limits, which the model must satisfy by construction.
    limits = {}
    for band, s in intensities.items():
        a0b, b0b = scaled(band, 300.0)
        limits[band] = {
            "zero_at_zero": absorptance(0.0, 1.0, a0b, b0b, s) == 0.0,
            "monotonic": all(absorptance(w, 1.0, a0b, b0b, s)
                             < absorptance(w * 1.05, 1.0, a0b, b0b, s)
                             for w in (1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0)),
            "weak_line_ratio_at_1e-6": absorptance(1e-6, 1.0, a0b, b0b, s) / (s * 1e-6),
        }
    for band, lim in limits.items():
        if not (lim["zero_at_zero"] and lim["monotonic"]
                and 0.99 < lim["weak_line_ratio_at_1e-6"] <= 1.0):
            raise SystemExit(f"{band} fails a limit: {lim}")

    if args.fortran:
        print("!     Donner and Ramanathan (1980) band model, Eq. (1)-(3).")
        print("!     Generated by exoplasim/scripts/trace_gas_band_model.py --fortran.")
        print("!     A0 scales as sqrt(T/300), beta0 as sqrt(300/T).")
        for band, s in intensities.items():
            b = BANDS[band]
            print(f"      parameter(za0_{band:9s} = {b['a0_300k']:8.3f})"
                  f"   ! bandwidth at 300 K, cm-1  ({b['centre_cm1']:.0f} cm-1)")
            print(f"      parameter(zbe_{band:9s} = {b['beta0_300k']:8.3f})"
                  f"   ! line shape at 300 K, 1 atm")
            print(f"      parameter(zsi_{band:9s} = {s:8.3f})"
                  f"   ! band intensity, cm-1 (cm atm)-1")
        return

    report = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
                             .replace(microsecond=0).isoformat(),
        "generator": "exoplasim/scripts/trace_gas_band_model.py",
        "task": "CLIM-42",
        "model": "Cess and Ramanathan (1972) / Ramanathan (1976), as used by "
                 "Donner and Ramanathan (1980) Eq. (1)-(3)",
        "loschmidt_per_cm3_stp": LOSCHMIDT,
        "bands": {band: dict(BANDS[band],
                             intensity_cm1_per_cm_atm=round(s, 3))
                  for band, s in intensities.items()},
        "intensity_sources": {
            "ch4_1306": "recovered from Donner and Ramanathan (1980) Table 2; "
                        "the paper takes it from Cess and Chen and does not "
                        "restate it",
            "n2o_1285": "McClatchey et al. (1973) Table 13, band system intensity",
            "n2o_589": "McClatchey et al. (1973) Table 13, band system intensity",
        },
        "table2_identity": {
            "rows_used": len(kept),
            "row_excluded": list(TABLE2_SUSPECT),
            "intensity": round(s_ch4, 2),
            "rms_cm1": round(rms, 4),
            "worst_relative_pct": round(worst_pct, 3),
            "declared_tolerance_pct": TOLERANCE_PCT,
            "intensity_including_excluded_row": round(s_all, 2),
            "rms_including_excluded_row": round(rms_all, 4),
        },
        "conversion_check": {
            "mcclatchey_ch4_1306": round(s_ch4_mcc, 2),
            "donner_table2_ch4_1306": round(s_ch4, 2),
            "ratio": round(ratio, 4),
            "reading": "rules out a unit error, which would be a factor of "
                       "1e19, 100 or 10; does not validate the value, because "
                       "two compilations of one band differ by tens of percent",
        },
        "limits": limits,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print("Donner and Ramanathan (1980) band model")
    print(f"  Table 2 identity: S(CH4 1306) = {s_ch4:.1f} cm-1 (cm atm)-1 "
          f"reproduces {len(kept)} rows to {worst_pct:.2f}% "
          f"(rms {rms:.3f} cm-1), bar {TOLERANCE_PCT}%")
    print(f"    including the suspect row: S = {s_all:.1f}, rms {rms_all:.3f}")
    print(f"  conversion check on CH4 1306: McClatchey {s_ch4_mcc:.1f} against "
          f"Donner {s_ch4:.1f}, ratio {ratio:.3f}")
    print("  band intensities, cm-1 (cm atm)-1:")
    for band, s in intensities.items():
        print(f"    {band:10s} A0 {BANDS[band]['a0_300k']:6.1f}  "
              f"beta0 {BANDS[band]['beta0_300k']:5.2f}  S {s:8.2f}")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
