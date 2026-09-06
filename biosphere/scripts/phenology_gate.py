"""The natural-vegetation phenology gate: where the seasonal landmarks come
from, and whether the model still reads them from there.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's leaf phenology, the days of its
simulation year the phenology keys on, and the source the model actually reads.

Summergreen phenology in LPJ-GUESS turns on two days of the year: the coldest,
where the growing-degree-day sum and the annual leaf-on sum are reset and
chilling detection is switched off, and the warmest, where chilling detection is
switched back on. The vendored source carried both as fixed ordinal dates on the
Earth calendar, one pair per hemisphere. This world's year is shorter than half
an Earth year, so the southern date named a day outside it and the southern
resets never fired at all, while northern chilling detection was switched off
and never switched back on.

They are now derived per gridcell from the temperature forcing itself:
`Climate::coldest_day` and `Climate::warmest_day`, read off a running day-of-year
mean of the air temperature. This module is the enforcement, and it can fail:

  ordinal    an Earth ordinal date still deciding a natural-vegetation
             phenology event, in any file the model compiles
  wiring     a reader of the landmarks not keyed on the derived pair, or the
             derivation no longer called where the forcing arrives
  bound      the chill-day count able to leave the lookup table it indexes
  restart    a landmark field missing from Climate::serialize, which would
             silently reset the seasonal cycle of every gridcell on resume

Ten fixtures run on every invocation, four of them built to be wrong in a named
way. A fixture that does not get the verdict it was built for is a defect in this
checker rather than in the model. They exercise the operator as the source
declares it over a northern, a southern, an equatorial, a flat and a
sixteen-orbit drifting forcing, and over a gridcell whose air temperature never
reaches the chilling base at all.

    python biosphere/scripts/phenology_gate.py            # status, exit 0

Every natural-vegetation reader of a seasonal landmark is covered by the wiring
list, the summergreen leaf litter release in `modules/somdynam.cpp` included: it
sheds over the month containing `Climate::coldest_day`, the same landmark the
degree-day and leaf-on sums reset on, rather than over a hemisphere's January or
July. A month is an Earth calendar artifact and this world's year is not twelve
of anything, so the month is taken from `Date::month_of` on this world's own
month lengths.

Nothing here is verified by execution. LPJ-GUESS builds and runs on this tree,
but no run it has made has passed acceptance (`world-qcse` has why), so the
source statements are made against the text of the source, and the fixtures run
a Python statement of the same operator rather than the compiled one.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT
from orbit import model_year_days  # noqa: E402
from paths import rel  # noqa: E402

REPORT = GENERATED / "phenology_gate_report.json"

FRAMEWORK = GUESS_SOURCE / "framework"
MODULES = GUESS_SOURCE / "modules"

# The Earth ordinal dates that used to decide these events, as they appeared in
# the source. Any of them back in a compiled file is a finding.
EARTH_ORDINALS = ["COLDEST_DAY_NHEMISPHERE", "COLDEST_DAY_SHEMISPHERE",
                  "WARMEST_DAY_NHEMISPHERE", "WARMEST_DAY_SHEMISPHERE"]


def read_source(path: Path) -> str:
    """Source text, whatever byte encoding the vendored file happens to be in.

    Two crop sources in this tree are ISO-8859 rather than UTF-8, so grep skips
    them by default and a reader that assumes UTF-8 raises. Both of them read
    the seasonal landmarks, so neither may be skipped here.
    """
    return path.read_text(encoding="latin-1")


def compiled_sources() -> list[Path]:
    """Every source the model builds, which is what a check has to cover."""
    paths = sorted(FRAMEWORK.glob("*.cpp")) + sorted(FRAMEWORK.glob("*.h"))
    paths += sorted(MODULES.glob("*.cpp")) + sorted(MODULES.glob("*.h"))
    return paths


# ---------------------------------------------------------------------------
# The operator, as the source declares it
# ---------------------------------------------------------------------------

def landmarks(cycle: list[float]) -> tuple[int, int]:
    """The coldest and warmest day of a seasonal cycle.

    A statement in Python of Climate::find_seasonal_landmarks. The extremum is
    taken of the cycle smoothed over one of this world's months, circularly,
    because a landmark is the centre of a season and not one day's weather. Where
    the cycle is flat to the last bit the two searches return the same day, and
    the warmest is then placed half a year from the coldest: the only choice that
    assumes nothing about which half of the year is which and still fires each
    reset once per year.
    """
    length = len(cycle)
    window = max(1, length // 12)

    cold = warm = 0
    coldest = warmest = 0.0
    for day in range(length):
        total = sum(cycle[(day + k - window // 2) % length] for k in range(window))
        if day == 0 or total < coldest:
            coldest, cold = total, day
        if day == 0 or total > warmest:
            warmest, warm = total, day

    if warm == cold:
        warm = (cold + length // 2) % length
    return cold, warm


def run_forcing(forcing, *, orbits: int, length: int, nyear_seasonal: int = 20,
                fixed_landmarks: tuple[int, int] | None = None,
                guard_le: bool = False, supplied: bool = True) -> dict:
    """Run the phenology counters over a forcing and report what happened.

    `forcing(year, day)` is the daily mean air temperature in degrees C.
    `fixed_landmarks` replaces the derived pair, which is how a fixture states
    what an Earth ordinal date does here. `guard_le` restores the chill-day guard
    the vendored source shipped. `supplied` is the Vesper input module, which
    hands the whole year over before it starts; False is a module that only
    delivers a day at a time, and reaches the same landmarks by accumulation.
    """
    cycle = [0.0] * length
    years_seen = 0
    cold, warm = 0, length // 2
    if fixed_landmarks is not None:
        cold, warm = fixed_landmarks

    gdd5 = 0.0
    chilldays = 0
    max_index = 0
    per_year = []

    for year in range(orbits):
        # The record is kept whatever decides the landmarks, so that a fixture
        # pinning them to a calendar date can still be told where the season
        # actually was.
        if supplied:
            weight = 1.0 / min(years_seen + 1, nyear_seasonal)
            cycle = [(1.0 - weight) * cycle[d] + weight * forcing(year, d)
                     for d in range(length)]
            if fixed_landmarks is None:
                cold, warm = landmarks(cycle)

        resets = {"gdd5": 0, "aphen": 0, "sensechill_on": 0}
        sensechill = True
        mtemp_last = None

        for day in range(length):
            temp = forcing(year, day)

            if not supplied:
                weight = 1.0 / min(years_seen + 1, nyear_seasonal)
                cycle[day] = (1.0 - weight) * cycle[day] + weight * temp

            if day == cold:
                gdd5 = 0.0
                sensechill = False
                resets["gdd5"] += 1
                resets["aphen"] += 1
            elif day == warm:
                sensechill = True
                resets["sensechill_on"] += 1

            gdd5 += max(0.0, temp - 5.0)
            if temp < 5.0 and (chilldays <= length if guard_le else chilldays < length):
                chilldays += 1
            max_index = max(max_index, chilldays)

            # The monthly-mean crossing that is the count's only other reset.
            month = min(11, day * 12 // length)
            window = max(1, length // 12)
            mtemp = sum(forcing(year, (day - k) % length)
                        for k in range(window)) / window
            if mtemp_last is not None and mtemp_last >= 5.0 and mtemp < 5.0 and sensechill:
                gdd5 = 0.0
                chilldays = 0
            mtemp_last = mtemp

        years_seen += 1
        if not supplied and fixed_landmarks is None:
            cold, warm = landmarks(cycle)
        per_year.append(resets)

    return {"coldest_day": cold, "warmest_day": warm,
            "max_chill_index": max_index, "per_year": per_year,
            "cycle": cycle}


def verdicts(result: dict, length: int) -> set[str]:
    """What went wrong in a run, by name. Empty is a correct run."""
    found = set()
    cold, warm = result["coldest_day"], result["warmest_day"]

    if not (0 <= cold < length and 0 <= warm < length):
        found.add("unreachable")
    if cold == warm:
        found.add("collapsed")
    if result["max_chill_index"] > length:
        found.add("overflow")

    for resets in result["per_year"]:
        for event, count in resets.items():
            if count > 1:
                found.add("multiple_resets")
            if count == 0:
                found.add("missing_reset")

    # A landmark is the centre of a season, and the search that finds it smooths
    # over one month, so a reset more than a month from the derived minimum is
    # not on the season it claims to be on. The criterion is the operator's own
    # window and is fixed here rather than after a run.
    cycle = result["cycle"]
    if any(value != cycle[0] for value in cycle) and 0 <= cold < length:
        derived_cold, _ = landmarks(cycle)
        separation = abs(cold - derived_cold)
        separation = min(separation, length - separation)
        if separation > max(1, length // 12):
            found.add("misplaced")
    return found


# ---------------------------------------------------------------------------
# Source checks
# ---------------------------------------------------------------------------

def check_source() -> list[dict]:
    findings = []

    for path in compiled_sources():
        text = read_source(path)
        for symbol in EARTH_ORDINALS:
            if symbol in text:
                findings.append({
                    "kind": "ordinal", "what": symbol,
                    "detail": f"an Earth ordinal date is back in {rel(path)}"})

    guess_h = read_source(FRAMEWORK / "guess.h")
    guess_cpp = read_source(FRAMEWORK / "guess.cpp")
    driver = read_source(MODULES / "driver.cpp")
    growth = read_source(MODULES / "growth.cpp")
    vesperinput = read_source(MODULES / "vesperinput.cpp")
    somdynam = read_source(MODULES / "somdynam.cpp")

    wiring = [
        ("the coldest-day reset of the summergreen degree-day sum",
         driver, r"date\.day\s*==\s*climate\.coldest_day"),
        ("the warmest-day return of chilling detection",
         driver, r"date\.day\s*==\s*climate\.warmest_day"),
        ("the daily record of the seasonal cycle",
         driver, r"climate\.accumulate_seasonal_cycle\(\)"),
        ("the yearly re-derivation of the landmarks",
         driver, r"climate\.find_seasonal_landmarks\(\)"),
        ("the coldest-day reset of the annual leaf-on sum",
         growth, r"date\.day\s*==\s*climate\.coldest_day"),
        ("the summergreen leaf litter release, on the coldest day's month",
         somdynam,
         r"date\.month_of\(patch\.get_climate\(\)\.coldest_day\)"),
        ("the month lengths that release is placed on",
         guess_h, r"int\s+month_of\(int\s+julian_day\)\s*const"),
        ("the Vesper forcing handing over the whole year",
         vesperinput, r"climate\.set_seasonal_cycle\(dtemp\)"),
        ("the seasonal cycle the landmarks are read off",
         guess_h, r"double\s+dtemp_seasonal\[Date::MAX_YEAR_LENGTH\]"),
        ("the window the seasonal cycle is averaged over",
         guess_h, r"const\s+int\s+NYEAR_SEASONAL\s*=\s*\d+"),
        ("the circular search for the extremum",
         guess_cpp, r"void\s+Climate::find_seasonal_landmarks"),
    ]
    for what, text, pattern in wiring:
        if not re.search(pattern, text):
            findings.append({"kind": "wiring", "what": what,
                             "detail": f"no match for /{pattern}/"})

    # The chill-day count indexes Pft::gdd0, which is one longer than the year.
    table = re.search(r"double\s+gdd0\[Date::MAX_YEAR_LENGTH\s*\+\s*(\d+)\]", guess_h)
    guard = re.search(r"climate\.chilldays\s*(<=?)\s*Date::MAX_YEAR_LENGTH", driver)
    if table is None:
        findings.append({"kind": "bound", "what": "Pft::gdd0",
                         "detail": "the lookup table is no longer sized from the year length"})
    elif guard is None:
        findings.append({"kind": "bound", "what": "the chill-day guard",
                         "detail": "no guard on the count that indexes Pft::gdd0"})
    else:
        highest = int(table.group(1)) - 1 if guard.group(1) == "<" else int(table.group(1))
        if highest > int(table.group(1)) - 1:
            findings.append({
                "kind": "bound", "what": "the chill-day guard",
                "detail": "the count may reach MAX_YEAR_LENGTH + "
                          f"{highest}, past the last entry of a table sized "
                          f"MAX_YEAR_LENGTH + {table.group(1)}"})

    serialized = re.search(r"void Climate::serialize.*?\n}", guess_cpp, re.S)
    if serialized is None:
        findings.append({"kind": "restart", "what": "Climate::serialize",
                         "detail": "not found"})
    else:
        for field in ["coldest_day", "warmest_day", "dtemp_seasonal",
                      "seasonal_cycle_years", "seasonal_cycle_supplied"]:
            if not re.search(rf"&\s*{field}\b", serialized.group(0)):
                findings.append({
                    "kind": "restart", "what": field,
                    "detail": "not in Climate::serialize, so a resumed run would "
                              "lose the seasonal cycle it derived"})
    return findings


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _fixtures(length: int) -> list[dict]:

    def seasonal(mean, amplitude, phase):
        return lambda year, day: mean - amplitude * math.cos(
            2.0 * math.pi * ((day / length) - phase))

    def drifting(year, day):
        # A star that is not constant: amplitude and phase both move, which is
        # the forcing a driver file of several years carries.
        amplitude = 12.0 + 3.0 * math.sin(2.0 * math.pi * year / 16.0)
        phase = 0.02 * math.sin(2.0 * math.pi * year / 16.0)
        return 5.0 - amplitude * math.cos(2.0 * math.pi * ((day / length) - phase))

    north = seasonal(5.0, 15.0, 0.0)
    south = seasonal(5.0, 15.0, 0.5)

    cases = [
        ("a northern cell, one repeating orbit",
         dict(forcing=north, orbits=3), None),
        ("a southern cell, the same cycle half a year out of phase",
         dict(forcing=south, orbits=3), None),
        ("an equatorial cell, a seasonal cycle of a fifth of a degree",
         dict(forcing=seasonal(26.0, 0.1, 0.3), orbits=3), None),
        ("a cell with no seasonal cycle at all",
         dict(forcing=lambda year, day: 18.0, orbits=3), None),
        ("sixteen orbits of a variable star",
         dict(forcing=drifting, orbits=16), None),
        ("a cell whose air temperature never reaches the chilling base",
         dict(forcing=lambda year, day: -20.0, orbits=5), None),
        ("a module that delivers the forcing one day at a time",
         dict(forcing=north, orbits=3, supplied=False), None),
        ("the Earth ordinal dates, on this world's calendar",
         dict(forcing=north, orbits=3, fixed_landmarks=(14, 195)), "unreachable"),
        ("a landmark pair that collapses onto one day",
         dict(forcing=north, orbits=3, fixed_landmarks=(14, 14)), "collapsed"),
        ("the chill-day guard the vendored source shipped",
         dict(forcing=lambda year, day: -20.0, orbits=5, guard_le=True), "overflow"),
        ("a southern cell given the northern hemisphere's landmark",
         dict(forcing=south, orbits=3, fixed_landmarks=(14, 105)), "misplaced"),
    ]

    results = []
    for label, kwargs, expect in cases:
        result = run_forcing(length=length, **kwargs)
        found = verdicts(result, length)
        ok = (not found) if expect is None else (expect in found)
        results.append({"fixture": label, "expected": expect or "clean",
                        "found": sorted(found), "pass": ok,
                        "coldest_day": result["coldest_day"],
                        "warmest_day": result["warmest_day"],
                        "max_chill_index": result["max_chill_index"]})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    length = model_year_days(config)

    findings = check_source()
    fixtures = _fixtures(length)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "year_length_days": length,
        "source": rel(GUESS_SOURCE),
        "findings": findings,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The natural-vegetation seasonal landmarks, checked against the source.\n")
        print(f"  simulation year: {length} days")
        broken = [case for case in fixtures if not case["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, "
                  f"found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  every landmark reader is on the derived pair, leaf litter")
            print("  included, the chill-day count cannot leave its table, and a")
            print("  resumed run keeps its cycle")
        print(f"\n  report: {rel(REPORT)}")

    if broken := [case for case in fixtures if not case["pass"]]:
        print("\nA fixture did not get the verdict it was built for. That is a defect in",
              file=sys.stderr)
        print("this checker, not in the model.", file=sys.stderr)
        return 2
    if findings:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
