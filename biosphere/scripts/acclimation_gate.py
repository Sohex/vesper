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
  declaration
             `biosphere/config/respiration_acclimation.yaml` disagreeing with
             itself, or a run instruction naming a memory the declaration does
             not carry: outside the bracket, or set at all while the
             declaration says the memory is undeclared

Fixtures run on every invocation, several of them built to be wrong in a named
way: the shipped behaviour with no memory at all, a resume that drops the state,
a state seeded from zero rather than from the forcing, and four declarations the
declaration check must refuse. A fixture that does not get the verdict it was
built for is a defect in this checker rather than in the model.

The bracket, its sources and the one-factor sensitivity over it are DECLARED in
`biosphere/config/respiration_acclimation.yaml`; this module reads them rather
than carrying them, so the declaration is the thing a reader argues with. The
sensitivity is executed here on every invocation and its response goes into the
report: what fraction of this world's seasonal amplitude the growth temperature
keeps on each arm, how many days it lags, and what the basal multiplier's annual
range becomes against the instantaneous temperature the routine used to run on.
It carries no pass/fail bar, for the reason the declaration states. `--strict`
refuses while the run's value is undeclared.

    python biosphere/scripts/acclimation_gate.py            # status, exit 0
    python biosphere/scripts/acclimation_gate.py --strict   # refuses while the
                                                            # option is unusable

Nothing here is verified by execution. LPJ-GUESS does not build on this tree, so
the source statements are made against the text of the source, and the fixtures
run a Python statement of the same update rule rather than the compiled one.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, COMPONENT_ROOT, GENERATED, GUESS_SOURCE
from orbit import model_year_days  # noqa: E402
from paths import rel  # noqa: E402

REPORT = GENERATED / "acclimation_gate_report.json"
DECLARATION = COMPONENT_ROOT / "config" / "respiration_acclimation.yaml"

FRAMEWORK = GUESS_SOURCE / "framework"
MODULES = GUESS_SOURCE / "modules"

# The sentinel the declaration uses for a memory length nobody has chosen. It is
# not a placeholder for a default: the model refuses the acclimated path without
# a declared length, so this is what "the option is not in use" looks like.
UNDECLARED = "undeclared"

# The response of the basal rate to the growth temperature, from Sprugel et al.
# (1996) as the vendored source cites it.
RESP_ACC_EXPONENT = -0.008
RESP_ACC_REF_TEMP = 10.15


def read_source(path: Path) -> str:
    return path.read_text(encoding="latin-1")


def load_declaration(path: Path = DECLARATION) -> dict:
    """The declared memory bracket and the sensitivity registered over it."""
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def bracket_days(declaration: dict) -> tuple[float, float]:
    return tuple(float(v) for v in declaration["memory"]["bracket_days"])


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
# The declaration, checked against itself and against the run instruction
# ---------------------------------------------------------------------------

def check_declaration(declaration: dict, *, enabled: int | None,
                      run_tau: float | None) -> list[dict]:
    """The declared bracket agrees with itself and with what the run asks for.

    Split out from `main` and given only its inputs so the fixtures can drive it
    with declarations built to be wrong. Every case has a right answer: a
    bracket is or is not the interval its own two ends span, a run instruction
    does or does not name a memory the declaration carries.
    """
    findings = []
    memory = declaration.get("memory") or {}
    ends = memory.get("ends") or {}

    try:
        low, high = bracket_days(declaration)
    except Exception as exc:
        return [{"kind": "declaration", "what": "memory.bracket_days",
                 "detail": f"not a pair of numbers: {exc}"}]

    if not low < high:
        findings.append({"kind": "declaration", "what": "memory.bracket_days",
                         "detail": f"{low} is not below {high}, so the bracket "
                                   "is not an interval"})

    for name, expected in (("fast", low), ("slow", high)):
        end = ends.get(name) or {}
        value = end.get("value_days")
        if value is None or float(value) != expected:
            findings.append({
                "kind": "declaration", "what": f"memory.ends.{name}",
                "detail": f"value_days {value} is not the {name} end of the "
                          f"declared bracket, {expected}"})
        if not str(end.get("source") or "").strip():
            findings.append({
                "kind": "declaration", "what": f"memory.ends.{name}",
                "detail": "carries no source. An end may be a convention, but "
                          "it may not be silent about being one"})

    arms = (declaration.get("sensitivity") or {}).get("arms_days")
    if arms is None or [float(a) for a in arms] != [low, high]:
        findings.append({
            "kind": "declaration", "what": "sensitivity.arms_days",
            "detail": f"{arms} does not span the declared bracket [{low}, {high}], "
                      "so the registered sensitivity is over some other factor "
                      "range than the one the run is bracketed by"})

    declared = memory.get("run_value_days")
    if declared != UNDECLARED:
        try:
            declared = float(declared)
        except (TypeError, ValueError):
            findings.append({
                "kind": "declaration", "what": "memory.run_value_days",
                "detail": f"{declared!r} is neither a number nor {UNDECLARED!r}"})
            declared = None
        else:
            if not low <= declared <= high:
                findings.append({
                    "kind": "declaration", "what": "memory.run_value_days",
                    "detail": f"{declared} is outside the declared bracket "
                              f"[{low}, {high}]"})
            if not enabled:
                findings.append({
                    "kind": "declaration", "what": "memory.run_value_days",
                    "detail": f"names {declared} days while the run takes the "
                              "standard respiration path, so nothing is on that "
                              f"arm. {UNDECLARED!r} is what an unused memory says"})

    if run_tau is not None:
        if declared == UNDECLARED:
            findings.append({
                "kind": "declaration", "what": "acclim_resp_tau",
                "detail": f"the run instruction declares {run_tau} days and the "
                          "declaration says the memory is undeclared, so the run "
                          "carries a physiological memory nothing argued for"})
        elif declared is not None and float(run_tau) != declared:
            findings.append({
                "kind": "declaration", "what": "acclim_resp_tau",
                "detail": f"the run instruction declares {run_tau} days and the "
                          f"declaration says {declared}"})
        if not low <= float(run_tau) <= high:
            findings.append({
                "kind": "declaration", "what": "acclim_resp_tau",
                "detail": f"the run instruction declares {run_tau} days, outside "
                          f"the declared bracket [{low}, {high}]"})
    return findings


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _reference_declaration() -> dict:
    """A well-formed declaration, written here and not read from the tree.

    The declaration fixtures mutate this rather than the live file, so a fixture
    verdict is a statement about `check_declaration` and never about
    `biosphere/config/respiration_acclimation.yaml`. The live file is checked
    once, into `findings`.
    """
    return {
        "version": 1,
        "memory": {
            "bracket_days": [7.0, 30.0],
            "ends": {
                "fast": {"value_days": 7.0, "sourced": True, "source": "a paper"},
                "slow": {"value_days": 30.0, "sourced": True, "source": "a paper"},
            },
            "run_value_days": UNDECLARED,
        },
        "sensitivity": {"factor": "acclim_resp_tau", "arms_days": [7.0, 30.0],
                        "responds": [], "threshold": "none",
                        "propagated_as": "model-form uncertainty"},
    }


def _fixtures(length: int, declaration: dict) -> list[dict]:
    low, high = bracket_days(declaration)
    tau = (low + high) / 2.0
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

    # The declaration check, driven from the REFERENCE declaration rather than
    # the live one. A fixture that read the live file would report an edit to
    # biosphere/config/respiration_acclimation.yaml as a defect in this checker,
    # which is the wrong verdict on the wrong artifact: the live declaration is
    # what `findings` is for. Positive first, because a gate whose negatives all
    # fire and whose positives were never built proves only that it can say no.
    ref = _reference_declaration()
    ref_low, ref_high = bracket_days(ref)

    cases.append(("a declaration with no memory in use and the standard path is "
                  "accepted", not check_declaration(ref, enabled=0, run_tau=None),
                  "clean"))

    on_slow = copy.deepcopy(ref)
    on_slow["memory"]["run_value_days"] = ref_high
    cases.append(("a run on the slow arm, declared and asked for, is accepted",
                  not check_declaration(on_slow, enabled=1, run_tau=ref_high),
                  "clean"))

    # Built to be wrong: a memory outside the bracket the declaration argues.
    cases.append(("a run instruction asking for a memory past the slow end is "
                  "refused",
                  bool(check_declaration(on_slow, enabled=1, run_tau=ref_high * 2.0)),
                  "outside_bracket"))

    # Built to be wrong: a run carrying a memory the declaration never named.
    cases.append(("a run instruction with a memory the declaration leaves "
                  "undeclared is refused",
                  bool(check_declaration(ref, enabled=1, run_tau=ref_high)),
                  "undeclared_memory"))

    # Built to be wrong: a bracket that is not the interval its own ends span.
    widened = copy.deepcopy(ref)
    widened["memory"]["bracket_days"] = [ref_low, ref_high * 3.0]
    cases.append(("a bracket that is not the interval its own ends span is refused",
                  bool(check_declaration(widened, enabled=0, run_tau=None)),
                  "bracket_disagrees"))

    # Built to be wrong: an end that does not say where it came from.
    silent = copy.deepcopy(ref)
    silent["memory"]["ends"]["slow"]["source"] = ""
    cases.append(("an end that carries no source is refused",
                  bool(check_declaration(silent, enabled=0, run_tau=None)),
                  "silent_end"))

    # Built to be wrong: a sensitivity registered over some other range.
    narrowed = copy.deepcopy(ref)
    narrowed["sensitivity"]["arms_days"] = [ref_low, (ref_low + ref_high) / 2.0]
    cases.append(("a sensitivity registered over a range that is not the bracket "
                  "is refused",
                  bool(check_declaration(narrowed, enabled=0, run_tau=None)),
                  "sensitivity_off_bracket"))

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
    declaration = load_declaration()

    # What the run instruction file actually asks for.
    runner = (Path(__file__).resolve().parent / "run_lpj_guess.py").read_text()
    enabled = re.search(r"^acclimated_respiration (\d)", runner, re.M)
    declared_tau = re.search(r"^acclim_resp_tau ([0-9.]+)", runner, re.M)
    baseline_enabled = int(enabled.group(1)) if enabled else None
    baseline_tau = float(declared_tau.group(1)) if declared_tau else None

    findings = check_source() + check_declaration(
        declaration, enabled=baseline_enabled, run_tau=baseline_tau)
    fixtures = _fixtures(length, declaration)

    # The registered one-factor sensitivity, executed. The arms come from the
    # declaration, so what is reported here is what was registered there and not
    # a range this module chose.
    registered = declaration["sensitivity"]
    arms = [float(a) for a in registered["arms_days"]]
    bracket = [seasonal_response(tau, length) for tau in arms]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "year_length_days": length,
        "source": rel(GUESS_SOURCE),
        "declaration": rel(DECLARATION),
        "findings": findings,
        "fixtures": fixtures,
        "tau_bracket_days": list(bracket_days(declaration)),
        "tau_bracket_sources": {
            name: (end.get("source") or "").strip()
            for name, end in declaration["memory"]["ends"].items()},
        "tau_bracket_sourced": {
            name: bool(end.get("sourced"))
            for name, end in declaration["memory"]["ends"].items()},
        "registered_sensitivity": {
            "factor": registered["factor"],
            "arms_days": arms,
            "responds": registered["responds"],
            "threshold": registered["threshold"],
            "propagated_as": registered["propagated_as"],
            "response": bracket,
        },
        "seasonal_response_over_bracket": bracket,
        "declared_run_value_days": declaration["memory"]["run_value_days"],
        "baseline_acclimated_respiration": baseline_enabled,
        "baseline_acclim_resp_tau": baseline_tau,
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
        print(f"\n  the registered one-factor sensitivity over "
              f"{registered['factor']}, executed:")
        for entry in bracket:
            print(f"    tau {entry['tau_days']:>4} d: keeps {entry['amplitude_kept']:.2f} "
                  f"of the cycle's amplitude, lags it {entry['lag_days']} days, "
                  f"basal rate ranges {entry['basal_rate_range_acclimated']:.2f} "
                  f"against {entry['basal_rate_range_instantaneous']:.2f} instantaneous")
        print(f"    threshold: {registered['threshold']}, carried as "
              f"{registered['propagated_as']}")
        print(f"\n  declaration: {rel(DECLARATION)}, run_value_days "
              f"{report['declared_run_value_days']}")
        print(f"  baseline: acclimated_respiration "
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
