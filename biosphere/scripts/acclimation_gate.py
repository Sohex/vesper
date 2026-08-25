"""The respiration acclimation gate: what the growth temperature is, and whether
the model still keeps one.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's autotrophic respiration, the memory its
basal rate is supposed to carry, and the source the model actually reads.

`respiration_acclimated()` replaces each simulated plant functional type's
`respcoeff` with `f(T_acc) = f_ref * 10^(-0.008 (T_acc - 10.15))`, where `T_acc`
is the GROWTH temperature the tissue has adjusted to. It was being handed the
current day's air temperature and the current 25 cm soil temperature, so the same
variable drove both the acute Lloyd and Taylor response and the acclimation
multiplier, which is no acclimation at all.

`Climate::tacc_air` and `Soil::tacc_root` carry the memory now: exponential
running means of the daily air and root-zone temperature with an e-folding time
of `acclim_resp_tau` ABSOLUTE days, seeded from the first temperature each
gridcell sees and serialized so a resumed run keeps them. This module is the
enforcement, and it can fail:

  wiring     the acclimation multiplier reading an instantaneous temperature
             again, or a growth temperature no longer updated once a day from
             the variable it is a mean of
  restart    a growth temperature missing from a serialize list, which would
             silently re-acclimate every gridcell from scratch on resume
  default    a default value appearing for acclim_resp_tau, or the refusal that
             requires it going missing. A default here is an undeclared
             physiological memory
  seeding    the growth temperature initialised to a constant rather than to
             the first temperature the gridcell sees

Eight fixtures run on every invocation, three of them built to be wrong in a
named way: the shipped behaviour with no memory at all, a resume that drops the
state, and a state seeded from zero rather than from the forcing. A fixture that
does not get the verdict it was built for is a defect in this checker rather than
in the model.

The gate also reports what the memory is worth on this world's seasonal cycle,
over the bracket the held literature supports, because the e-folding time is
BRACKETED and not measured. Gifford (2003) reports plant respiration acclimating
to a temperature change in as little as a week; QUINCY gives its lagged responses
a process-specific memory whose length is in a supplement this project does not
hold. So the bracket's fast end is sourced and its slow end is a convention, and
`--strict` refuses while the run's value is undeclared.

    python biosphere/scripts/acclimation_gate.py            # status, exit 0
    python biosphere/scripts/acclimation_gate.py --strict   # refuses while the
                                                            # option is unusable

Nothing here is verified by execution. LPJ-GUESS does not build on this tree, so
the source statements are made against the text of the source, and the fixtures
run a Python statement of the same update rule rather than the compiled one.
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

from _paths import CONFIG, GENERATED, GUESS_SOURCE
from orbit import model_year_days  # noqa: E402
from paths import rel  # noqa: E402

REPORT = GENERATED / "acclimation_gate_report.json"

FRAMEWORK = GUESS_SOURCE / "framework"
MODULES = GUESS_SOURCE / "modules"

# The e-folding time of the growth temperature, in ABSOLUTE days. The fast end is
# Gifford (2003), which reports respiration acclimating to a temperature change
# in as little as a week. The slow end is the month-long growth-temperature
# window the acclimating land-surface literature conventionally uses, and this
# project holds no source for it, so it is a convention and is labelled one.
TAU_BRACKET_DAYS = (7.0, 30.0)
TAU_BRACKET_SOURCES = {
    "fast": "Gifford (2003), references/gifford2003-plant-respiration.pdf: plant "
            "respiration acclimates to a temperature change in as little as a week",
    "slow": "convention, unsourced in this tree. QUINCY's process-specific memory "
            "lengths are in a supplement this project does not hold",
}

# The response of the basal rate to the growth temperature, from Sprugel et al.
# (1996) as the vendored source cites it.
RESP_ACC_EXPONENT = -0.008
RESP_ACC_REF_TEMP = 10.15


def read_source(path: Path) -> str:
    return path.read_text(encoding="latin-1")


# ---------------------------------------------------------------------------
# The update rule, as the source declares it
# ---------------------------------------------------------------------------

def acclimate(state: float | None, temp: float, tau: float) -> float:
    """One day of acclimation. A statement in Python of driver.cpp:acclimate."""
    if state is None:
        return temp
    return state + (1.0 - math.exp(-1.0 / tau)) * (temp - state)


def basal_rate(tacc: float) -> float:
    """The multiplier that replaces respcoeff, at a growth temperature."""
    return 10.0 ** (RESP_ACC_EXPONENT * (tacc - RESP_ACC_REF_TEMP))


def trajectory(forcing, *, days: int, tau: float, seed: float | None = None,
               restart_at: int | None = None, drop_state: bool = False,
               no_memory: bool = False) -> list[float]:
    """The growth temperature day by day.

    `seed` replaces the seeding from the first temperature, which is how a
    fixture states what initialising to a constant does. `restart_at` serializes
    and resumes on that day, and `drop_state` makes the resume lose the state.
    `no_memory` is the shipped behaviour, the instantaneous temperature.
    """
    state = seed
    out = []
    for day in range(days):
        temp = forcing(day)
        if no_memory:
            out.append(temp)
            continue
        if restart_at is not None and day == restart_at:
            state = None if drop_state else state
        state = acclimate(state, temp, tau)
        out.append(state)
    return out


def seasonal_response(tau: float, length: int, amplitude: float = 15.0,
                      mean: float = 5.0) -> dict:
    """What the memory is worth over one seasonal cycle of the growth temperature.

    Reported rather than asserted. A first-order lag driven at the seasonal
    frequency attenuates by 1/sqrt(1 + (omega tau)^2) and lags by atan(omega tau);
    the numbers below are measured off the trajectory rather than read off that
    result, and the two are printed together because a disagreement between them
    would mean the update rule is not the lag it is documented as.
    """
    def forcing(day):
        return mean - amplitude * math.cos(2.0 * math.pi * day / length)

    # Long enough for the transient from seeding to have decayed.
    days = length * 8
    series = trajectory(forcing, days=days, tau=tau)
    last = series[-length:]
    attenuation = (max(last) - min(last)) / (2.0 * amplitude)
    lag = (last.index(min(last)) - 0) % length
    if lag > length // 2:
        lag -= length

    omega = 2.0 * math.pi / length
    analytic_attenuation = 1.0 / math.sqrt(1.0 + (omega * tau) ** 2)
    analytic_lag = math.atan(omega * tau) / omega

    return {
        "tau_days": tau,
        "amplitude_kept": round(attenuation, 4),
        "amplitude_kept_first_order_lag": round(analytic_attenuation, 4),
        "lag_days": lag,
        "lag_days_first_order": round(analytic_lag, 2),
        "basal_rate_range_instantaneous": round(
            basal_rate(mean - amplitude) / basal_rate(mean + amplitude), 4),
        "basal_rate_range_acclimated": round(
            basal_rate(min(last)) / basal_rate(max(last)), 4),
    }


# ---------------------------------------------------------------------------
# Source checks
# ---------------------------------------------------------------------------

def check_source() -> list[dict]:
    findings = []

    guess_h = read_source(FRAMEWORK / "guess.h")
    guess_cpp = read_source(FRAMEWORK / "guess.cpp")
    params_cpp = read_source(FRAMEWORK / "parameters.cpp")
    driver = read_source(MODULES / "driver.cpp")
    canexch = read_source(MODULES / "canexch.cpp")
    soil_cpp = read_source(MODULES / "soil.cpp")

    wiring = [
        ("the aboveground growth temperature", guess_h, r"double\s+tacc_air\s*;"),
        ("the root-zone growth temperature", guess_h, r"double\s+tacc_root\s*;"),
        ("the daily acclimation of the air state", driver,
         r"acclimate\(climate\.tacc_air,\s*climate\.tacc_air_set,\s*climate\.temp\)"),
        ("the daily acclimation of the root-zone state", driver,
         r"acclimate\(soil\.tacc_root,\s*soil\.tacc_root_set,\s*soiltemp25\)"),
        ("the exponential update with the declared e-folding time", driver,
         r"1\.0\s*-\s*exp\(-1\.0\s*/\s*acclim_resp_tau\)"),
        ("the acclimation multiplier reading the growth temperature", canexch,
         r"pow\(10,\s*\(-0\.008\s*\*\s*\(tacc_air\s*-\s*10\.15\)\)\)"),
        ("the root acclimation multiplier reading the root-zone state", canexch,
         r"pow\(10,\s*\(-0\.008\s*\*\s*\(tacc_root\s*-\s*10\.15\)\)\)"),
        ("the caller passing the states and not the day's temperature", canexch,
         r"respiration_acclimated\(gtemp,\s*patch\.soil\.gtemp,\s*"
         r"climate\.tacc_air,\s*patch\.soil\.tacc_root"),
    ]
    for what, text, pattern in wiring:
        if not re.search(pattern, text):
            findings.append({"kind": "wiring", "what": what,
                             "detail": f"no match for /{pattern}/"})

    # The acute response and the acclimation multiplier must be different
    # variables. gtemp is the acute one and must still reach the function.
    if not re.search(r"resp_sap\s*=.*gtemp_air", canexch):
        findings.append({"kind": "wiring", "what": "the acute temperature response",
                         "detail": "sapwood respiration no longer carries gtemp_air"})

    for field, text, where in [("tacc_air", guess_cpp, "Climate::serialize"),
                               ("tacc_air_set", guess_cpp, "Climate::serialize"),
                               ("tacc_root", soil_cpp, "Soil::serialize"),
                               ("tacc_root_set", soil_cpp, "Soil::serialize")]:
        block = re.search(rf"void {where}.*?\n}}", text, re.S)
        if block is None or not re.search(rf"&\s*{field}\b", block.group(0)):
            findings.append({
                "kind": "restart", "what": field,
                "detail": f"not in {where}, so a resumed run would re-acclimate "
                          "every gridcell from scratch"})

    declared = re.search(r'declareitem\("acclim_resp_tau"[^;]*;', params_cpp)
    if declared is None:
        findings.append({"kind": "default", "what": "acclim_resp_tau",
                         "detail": "not declared in the instruction file grammar"})
    if not re.search(r"acclimated_respiration\s*&&\s*!itemparsed\(\"acclim_resp_tau\"\)",
                     params_cpp):
        findings.append({
            "kind": "default", "what": "the refusal",
            "detail": "acclimated_respiration 1 no longer requires acclim_resp_tau, "
                      "so a run could take an undeclared physiological memory"})
    if re.search(r"acclim_resp_tau\s*=\s*[0-9]", params_cpp):
        findings.append({"kind": "default", "what": "acclim_resp_tau",
                         "detail": "has been given a default value"})

    if not re.search(r"if\s*\(!is_set\)\s*{\s*state\s*=\s*temp", driver):
        findings.append({"kind": "seeding", "what": "the growth temperature",
                         "detail": "no longer seeded from the first temperature the "
                                   "gridcell sees"})
    return findings


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _fixtures(length: int) -> list[dict]:
    tau = sum(TAU_BRACKET_DAYS) / 2.0
    warm = 25.0

    def constant(day):
        return warm

    def step(day):
        return 5.0 if day < 100 else 25.0

    def weather(day):
        return 15.0 + 10.0 * math.sin(day) + 5.0 * math.sin(day / 3.0)

    cases = []

    # Constant. A state seeded from the forcing is at the answer on day 0 and
    # stays there, so there is no spin-in to wait out.
    series = trajectory(constant, days=50, tau=tau)
    cases.append(("a constant temperature holds the growth temperature at it",
                  all(abs(value - warm) < 1e-12 for value in series), "clean"))

    # Step. After tau days the state has covered 1 - 1/e of the step, which is
    # what an e-folding time means, and this is the check that the declared
    # constant is the constant the rule uses.
    series = trajectory(step, days=200, tau=tau)
    covered = (series[100 + int(tau) - 1] - 5.0) / 20.0
    cases.append((f"a step is {round(covered, 3)} of the way closed after one e-folding time",
                  abs(covered - (1.0 - math.exp(-1.0))) < 0.02, "clean"))

    # The step approaches, never overshoots, and never runs backwards.
    after = series[100:]
    cases.append(("the approach to a step is monotone and never overshoots",
                  all(a <= b + 1e-12 for a, b in zip(after, after[1:]))
                  and max(after) <= 25.0 + 1e-12, "clean"))

    # Restart. Serializing and resuming reproduces the uninterrupted trajectory
    # exactly, because the state is what carries the memory across the boundary.
    plain = trajectory(weather, days=300, tau=tau)
    resumed = trajectory(weather, days=300, tau=tau, restart_at=150)
    cases.append(("a resume that keeps the state reproduces the run exactly",
                  plain == resumed, "clean"))

    # Weather is damped: the growth temperature must vary less than the forcing
    # that drives it, which is the whole point of a memory.
    forced = [weather(day) for day in range(300)]
    cases.append(("synoptic variation is damped, not tracked",
                  (max(plain[100:]) - min(plain[100:]))
                  < 0.5 * (max(forced[100:]) - min(forced[100:])), "clean"))

    # The seasonal statement agrees with a first-order lag.
    response = seasonal_response(tau, length)
    cases.append(("the seasonal attenuation is the first-order lag it is documented as",
                  abs(response["amplitude_kept"]
                      - response["amplitude_kept_first_order_lag"]) < 0.02, "clean"))

    # Built to be wrong: the shipped behaviour, no memory at all.
    shipped = trajectory(weather, days=300, tau=tau, no_memory=True)
    cases.append(("the shipped behaviour, the instantaneous temperature",
                  shipped[100:] == forced[100:], "no_memory"))

    # Built to be wrong: a resume that drops the state re-acclimates from the
    # day it resumes on and does not reproduce the run.
    dropped = trajectory(weather, days=300, tau=tau, restart_at=150, drop_state=True)
    cases.append(("a resume that drops the state does not reproduce the run",
                  dropped != plain, "lost_memory"))

    # Built to be wrong: seeding from zero in a warm cell. The basal rate is
    # wrong by this factor on the first day and decays back over the memory.
    seeded = trajectory(constant, days=100, tau=tau, seed=0.0)
    error = basal_rate(seeded[0]) / basal_rate(warm)
    cases.append((f"seeding from zero starts a warm cell {round(error, 3)} times "
                  "off its basal rate", abs(error - 1.0) > 0.05, "spin_in"))

    return [{"fixture": label, "expected": expect, "pass": bool(ok)}
            for label, ok, expect in cases]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while the acclimated path has no declared memory")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    length = model_year_days(config)

    findings = check_source()
    fixtures = _fixtures(length)
    bracket = [seasonal_response(tau, length) for tau in TAU_BRACKET_DAYS]

    # What the run instruction file actually asks for.
    runner = (Path(__file__).resolve().parent / "run_lpj_guess.py").read_text()
    enabled = re.search(r"^acclimated_respiration (\d)", runner, re.M)
    declared_tau = re.search(r"^acclim_resp_tau ([0-9.]+)", runner, re.M)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "year_length_days": length,
        "source": rel(GUESS_SOURCE),
        "findings": findings,
        "fixtures": fixtures,
        "tau_bracket_days": list(TAU_BRACKET_DAYS),
        "tau_bracket_sources": TAU_BRACKET_SOURCES,
        "seasonal_response_over_bracket": bracket,
        "baseline_acclimated_respiration": int(enabled.group(1)) if enabled else None,
        "baseline_acclim_resp_tau": float(declared_tau.group(1)) if declared_tau else None,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The respiration acclimation memory, checked against the source.\n")
        broken = [case for case in fixtures if not case["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']} (expected {case['expected']})")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  the growth temperatures exist, are updated daily, survive a")
            print("  resume, and no default stands in for the memory's length")
        print(f"\n  what the memory is worth over one seasonal cycle:")
        for entry in bracket:
            print(f"    tau {entry['tau_days']:>4} d: keeps {entry['amplitude_kept']:.2f} "
                  f"of the cycle's amplitude, lags it {entry['lag_days']} days, "
                  f"basal rate ranges {entry['basal_rate_range_acclimated']:.2f} "
                  f"against {entry['basal_rate_range_instantaneous']:.2f} instantaneous")
        print(f"\n  baseline: acclimated_respiration "
              f"{report['baseline_acclimated_respiration']}, "
              f"acclim_resp_tau {report['baseline_acclim_resp_tau']}")
        print(f"\n  report: {rel(REPORT)}")

    if [case for case in fixtures if not case["pass"]]:
        print("\nA fixture did not get the verdict it was built for. That is a defect in",
              file=sys.stderr)
        print("this checker, not in the model.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and report["baseline_acclimated_respiration"] and \
            report["baseline_acclim_resp_tau"] is None:
        print("\n--strict: refused. acclimated_respiration is on and no acclim_resp_tau",
              file=sys.stderr)
        print("is declared, so the growth temperature has no memory length.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
