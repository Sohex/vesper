"""Check a climate product against the ecological forcing field contract.

Worldbuilding frame: Vesper is a simulated super-Earth, its climate is
ExoPlaSim and its biosphere is LPJ-GUESS. Everything below is about the fields
those two models hand each other.

`biosphere/notes/ecological-forcing-field-contract.md` is the contract; this is
the half of it that can fail. It runs the identities in that note's "The
identities that can fail" section against an actual product, and it runs a set
of reduced fixtures built to be wrong in a named way on every invocation, so a
fixture that does not get the verdict it was built for is a defect in the
checker rather than in the product.

    python biosphere/scripts/check_forcing_contract.py                    # fixtures
    python biosphere/scripts/check_forcing_contract.py --climatology X.nc # and a product

Exit 0 when every fixture got its expected verdict and, where a product was
given, every REQUIRED check on it passed. A check can also report ABSENT, which
is what a field the producer does not write earns: it is a fact about the
product, not a pass and not a failure, and the contract forbids naming such a
field as though it were delivered.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from _paths import CONFIG, PROJECT_ROOT  # noqa: F401

PASS = "ok"
FAIL = "FAIL"
ABSENT = "absent"

# Every field the contract's rows resolve to on the ExoPlaSim side, with the
# expected sign of each. `+` means the values must be non-negative, `-` that
# they must be non-positive, `0` that no sign is asserted. The signs are not
# cosmetic: ExoPlaSim stores upward fluxes negative (radmod.f90 forms the
# shortwave upward flux as `zfu1 = -...` and then `dswfl = dfu + dfd`), so
# incident surface shortwave is `rss - ssru` and getting that sign wrong is a
# factor of two plus a flip.
SIGNS = {
    "rss": "+",     # surface net shortwave, positive downward
    "rst": "+",     # TOA net shortwave, positive downward
    "ssru": "-",    # surface shortwave upward
    "rsut": "-",    # TOA shortwave upward
    "stru": "-",    # surface longwave upward
    "hfls": "-",    # surface latent heat, upward
    "evap": "-",    # evaporation, upward
    "pr": "+",
    "prl": "+",
    "prc": "+",
    "prsn": "+",
    "czen": "+",
}

# Fields the contract needs and this producer does not write. Reported, never
# silently tolerated: a contract that names a field nothing produces has not
# closed anything. Each entry says where the field would have to come from.
NOT_PRODUCED = {
    "tasmax": "code 201 (atsama). Written by outmod.f90, absent from "
              "pyburn.ilibrary and from run_exoplasim.REGULAR_CODES.",
    "tasmin": "code 202 (atsami). Same.",
    "td2m": "code 168. In REGULAR_CODES and SNAPSHOT_CODES; outmod.f90 never "
            "writes it and pyburn does not derive it.",
    "uas": "code 165. Never written; there is no 10 m wind in this model's "
           "output.",
    "vas": "code 166. Never written.",
}

REQUIRED = ["tas", "ts", "pr", "prl", "prc", "prsn", "rss", "ssru", "rls",
            "stru", "ps", "hus", "maxt", "mint", "lsm"]


def _tol(field: np.ndarray) -> float:
    """Closure tolerance: relative to the field's own scale, not absolute.

    A precipitation rate is 1e-8 m/s and a radiative flux is 1e2 W/m2, so one
    absolute tolerance cannot serve both without either passing everything or
    failing everything.
    """
    scale = float(np.max(np.abs(field))) if field.size else 1.0
    return max(scale, 1.0) * 1e-9


def check_signs(data: dict) -> list[tuple[str, str, str]]:
    out = []
    for name, sign in SIGNS.items():
        if name not in data:
            continue
        v = np.asarray(data[name], dtype=float)
        bad = (v < -_tol(v)).sum() if sign == "+" else (
            (v > _tol(v)).sum() if sign == "-" else 0)
        out.append((
            f"{name} is {'non-negative' if sign == '+' else 'non-positive'}",
            PASS if bad == 0 else FAIL,
            "" if bad == 0 else f"{int(bad)} values on the wrong side of zero"))
    return out


def check_precipitation(data: dict) -> list[tuple[str, str, str]]:
    out = []
    if {"pr", "prl", "prc"} <= data.keys():
        pr = np.asarray(data["pr"], float)
        residual = float(np.max(np.abs(pr - np.asarray(data["prl"], float)
                                       - np.asarray(data["prc"], float))))
        ok = residual <= _tol(pr)
        out.append(("large-scale + convective closes on the total",
                    PASS if ok else FAIL, f"max residual {residual:.3e}"))
    if {"pr", "prsn"} <= data.keys():
        pr = np.asarray(data["pr"], float)
        prsn = np.asarray(data["prsn"], float)
        excess = float(np.max(prsn - pr))
        ok = excess <= _tol(pr)
        out.append((
            "snowfall is a SUBSET of the total, never a third addend",
            PASS if ok else FAIL,
            f"max prsn - pr = {excess:.3e}. rainmod.f90 forms dprs from the "
            f"same zprsc and zprsl that are already inside dprc and dprl, so "
            f"prl + prc + prsn double counts the snow."))
    return out


def check_shortwave(data: dict) -> list[tuple[str, str, str]]:
    if not {"rss", "ssru"} <= data.keys():
        return []
    rss = np.asarray(data["rss"], float)
    ssru = np.asarray(data["ssru"], float)
    incident = rss - ssru
    ok = bool(np.all(incident >= -_tol(rss)))
    out = [("incident shortwave = rss - ssru is non-negative",
            PASS if ok else FAIL,
            f"range {float(incident.min()):.3f} to {float(incident.max()):.3f} "
            f"W/m2")]
    lit = incident > 1.0
    if lit.any():
        implied = (-ssru[lit]) / incident[lit]
        ok = bool(np.all((implied >= -1e-9) & (implied <= 1.0 + 1e-9)))
        out.append(("the implied surface albedo lies in [0, 1]",
                    PASS if ok else FAIL,
                    f"{float(implied.min()):.3f} to {float(implied.max()):.3f}"))
    return out


def check_extrema(data: dict) -> list[tuple[str, str, str]]:
    """Which variable are maxt and mint extrema OF?

    An extremum over a window brackets the mean of the same variable over the
    same window, whatever the variable does inside it. So this identifies the
    variable rather than merely testing a plausible bound, and it is the check
    the wrong-variable defect was found by.
    """
    if not {"maxt", "mint"} <= data.keys():
        return []
    lo = np.asarray(data["mint"], float)
    hi = np.asarray(data["maxt"], float)
    out = [("maxt >= mint", PASS if bool(np.all(hi >= lo - 1e-6)) else FAIL,
            f"worst maxt - mint = {float(np.min(hi - lo)):.3f} K")]
    for name in ("ts", "tas"):
        if name not in data:
            continue
        v = np.asarray(data[name], float)
        outside = int(np.sum((v < lo - 1e-6) | (v > hi + 1e-6)))
        out.append((
            f"maxt and mint bracket {name}",
            PASS if outside == 0 else FAIL,
            "" if outside == 0 else
            f"{outside} of {v.size} outside, worst excursion "
            f"{float(np.max(np.maximum(lo - v, v - hi))):.3f} K"))
    return out


def check_declared_units(data: dict, units: dict) -> list[tuple[str, str, str]]:
    """A declared unit that the values contradict.

    pyburn declares relative humidity's units as `1` and computes it as a
    percentage clipped to [0, 100]. A consumer that trusts the attribute is out
    by a hundred, which is why the contract carries specific humidity instead.
    """
    out = []
    for name, unit in units.items():
        if name not in data or unit not in ("1", "nondimen"):
            continue
        v = np.asarray(data[name], float)
        top = float(np.max(np.abs(v))) if v.size else 0.0
        out.append((f"{name} declares units '{unit}' and stays within them",
                    PASS if top <= 1.5 else FAIL,
                    "" if top <= 1.5 else
                    f"reaches {top:.4g}, which is a percentage and not a "
                    f"fraction"))
    return out


def check_finite(data: dict) -> list[tuple[str, str, str]]:
    bad = sorted(n for n, v in data.items()
                 if not np.all(np.isfinite(np.asarray(v, float))))
    return [("every field is finite", PASS if not bad else FAIL,
             "" if not bad else f"non-finite in {bad}")]


def check_present(data: dict) -> list[tuple[str, str, str]]:
    out = []
    missing = [n for n in REQUIRED if n not in data]
    out.append(("every field the contract requires is present",
                PASS if not missing else FAIL,
                "" if not missing else f"missing {missing}"))
    for name, why in NOT_PRODUCED.items():
        if name in data:
            out.append((f"{name} is available", PASS, ""))
        else:
            out.append((f"{name} is available", ABSENT, why))
    return out


def run_checks(data: dict, units: dict | None = None
               ) -> list[tuple[str, str, str]]:
    """Every contract check that can be run on the fields given."""
    units = units or {}
    return (check_present(data) + check_finite(data) + check_signs(data)
            + check_precipitation(data) + check_shortwave(data)
            + check_extrema(data) + check_declared_units(data, units))


# ---------------------------------------------------------------------------
# Fixtures. Six of the seven are built to be wrong in a named way.
# ---------------------------------------------------------------------------

def _clean() -> dict:
    """A small artifact that satisfies every identity, by construction."""
    rng = np.random.default_rng(20260824)
    shape = (4, 3, 5)
    prl = rng.uniform(0.0, 2e-8, shape)
    prc = rng.uniform(0.0, 4e-8, shape)
    total = prl + prc
    prsn = total * rng.uniform(0.0, 1.0, shape)
    incident = rng.uniform(0.0, 400.0, shape)
    albedo = rng.uniform(0.05, 0.6, shape)
    ts = rng.uniform(240.0, 310.0, shape)
    swing = rng.uniform(0.5, 12.0, shape)
    return {
        "tas": ts - swing * rng.uniform(0.0, 0.5, shape),
        "ts": ts,
        "pr": total, "prl": prl, "prc": prc, "prsn": prsn,
        "rss": incident * (1.0 - albedo), "ssru": -incident * albedo,
        "rls": rng.uniform(-150.0, 60.0, shape),
        "stru": -rng.uniform(150.0, 600.0, shape),
        "ps": rng.uniform(800.0, 1050.0, shape),
        "hus": rng.uniform(0.0, 0.02, shape),
        "maxt": ts + swing, "mint": ts - swing,
        "lsm": (rng.uniform(0.0, 1.0, shape) > 0.5).astype(float),
        "czen": rng.uniform(0.0, 0.6, shape),
    }


def _fixtures() -> list[tuple[str, dict, dict, set]]:
    """Name, fields, declared units, and the EXACT set of checks that must fail.

    A set rather than one name, because a corruption legitimately breaks the
    checks derived from it as well: storing the upward shortwave positive breaks
    the incident term and the implied albedo along with the sign itself. Naming
    the whole set is what makes this a test rather than a smoke alarm, since a
    fixture that fails one MORE check than declared is as much a defect as one
    that fails one fewer.
    """
    out: list[tuple[str, dict, dict, set]] = [
        ("a clean artifact", _clean(), {}, set()),
    ]

    d = _clean()
    d["pr"] = d["pr"] + d["prsn"]
    out.append(("snowfall added as a third addend", d, {},
                {"large-scale + convective closes on the total"}))

    d = _clean()
    d["prsn"] = d["pr"] * 1.5
    out.append(("more snow than precipitation", d, {},
                {"snowfall is a SUBSET of the total, never a third addend"}))

    d = _clean()
    d["ssru"] = -d["ssru"]
    out.append(("upward shortwave stored positive", d, {},
                {"ssru is non-positive",
                 "incident shortwave = rss - ssru is non-negative",
                 "the implied surface albedo lies in [0, 1]"}))

    d = _clean()
    d["maxt"], d["mint"] = d["mint"], d["maxt"]
    out.append(("extrema swapped", d, {},
                {"maxt >= mint", "maxt and mint bracket ts",
                 "maxt and mint bracket tas"}))

    d = _clean()
    # Extrema of a DIFFERENT variable: the defect this checker was written for,
    # and the one that leaves every other identity intact.
    other = d["ts"] + 40.0
    d["maxt"] = other + 1.0
    d["mint"] = other - 1.0
    out.append(("extrema of a different variable", d, {},
                {"maxt and mint bracket ts", "maxt and mint bracket tas"}))

    d = _clean()
    d["hur"] = d["czen"] * 100.0
    out.append(("a percentage declared as a fraction", d, {"hur": "1"},
                {"hur declares units '1' and stays within them"}))
    return out


def run_fixtures(verbose: bool = False) -> int:
    """Every fixture must get the verdict it was built for. Returns failures."""
    failures = 0
    for name, data, units, expected in _fixtures():
        results = {label: (status, detail)
                   for label, status, detail in run_checks(data, units)}
        failed = {label for label, (status, _) in results.items()
                  if status == FAIL}
        ok = failed == expected
        why = "" if ok else (f"expected {sorted(expected)} to fail, "
                             f"got {sorted(failed)}")
        failures += 0 if ok else 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] fixture: {name}"
              + (f"\n           {why}" if why else ""))
        if verbose:
            for label in sorted(expected):
                if results.get(label, ("", ""))[1]:
                    print(f"           {label}: {results[label][1]}")
    return failures


def read_climatology(path: Path) -> tuple[dict, dict]:
    """Fields and their declared units, with level axes reduced to the lowest.

    The lowest MODEL LEVEL, not a 2 m diagnostic: `hus` and the winds carry a
    level axis and the contract's near-surface rows resolve to the bottom of it.
    """
    import netCDF4 as nc
    data: dict = {}
    units: dict = {}
    with nc.Dataset(path) as ds:
        for name in ds.variables:
            if name in ("lat", "lon", "lev", "levp", "time", "modes"):
                continue
            v = np.asarray(ds[name][:], dtype=float)
            if v.ndim == 4:
                v = v[:, -1]        # lowest model level
            if v.ndim != 3:
                continue
            data[name] = v
            units[name] = getattr(ds[name], "units", "")
    return data, units


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None,
                        help="a pyburn climatology or forcing artifact to "
                             "check. Without it only the fixtures run.")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    print("FIXTURES  (six of seven are built to be wrong in a named way)")
    failures = run_fixtures(args.verbose)

    if args.climatology is not None:
        print(f"\nPRODUCT   {args.climatology}")
        data, units = read_climatology(args.climatology)
        for label, status, detail in run_checks(data, units):
            mark = {PASS: "  ok  ", FAIL: " FAIL ", ABSENT: "absent"}[status]
            print(f"[{mark}] {label}" + (f"\n           {detail}" if detail
                                         else ""))
            if status == FAIL:
                failures += 1

    print(f"\n{failures} failed" if failures else "\nno failures")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
