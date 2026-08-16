"""Rescale the PFT bioclimatic limits that are annual sums, and nothing else.

LPJ-GUESS's plant functional types are Earth's, and keeping them is a declared
choice: this is an Earth-analogue biosphere, not a prediction of alien
physiology. But some of their bioclimatic limits are calibrated per *year*, and
Vesper's year is half of Earth's, so those limits mean something different here.

Two kinds of parameter, and only one of them scales.

**Annual sums scale.** `gdd5min_est` is the minimum growing degree-days above
5 C a PFT needs *in a year* to establish. Degree-days accumulate per day at the
same rate on both worlds, so a 181-day year reaches roughly half the annual total
for identical temperatures. Left alone, Earth thresholds exclude nearly every
tree PFT for reasons that have nothing to do with the climate. This is not
hypothetical: the patched model run on Earth's own demo data collapses boreal
needleleaf and temperate broadleaf to grass.

**Everything else does not.** `tcmin_surv`, `tcmin_est`, `tcmax_est`,
`twmin_est` are temperatures, which are physiology and carry over unchanged.
`phengdd5ramp` is a *within-season* accumulation counted from the start of the
growing season, so it too is already in absolute time and must not be touched.
Rescaling it would be a real error, which is why it is named here explicitly
rather than merely omitted.

The scale factor is derived from the configured orbit, not written down, because
the year length moves with the stellar flux.

    python biosphere/scripts/build_vesper_pfts.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT

import orbit

EARTH_YEAR_DAYS = orbit.EARTH_SIDEREAL_YEAR_DAYS

# Two kinds of scaling, in opposite directions, and confusing them would be
# worse than doing neither.
#
# SCALED_DOWN are annual *sums*. A shorter year accumulates less, so a threshold
# has to come down by the year ratio to mean the same climate.
SCALED_DOWN = ("gdd5min_est", "gdd5min")

# SCALED_UP are *durations counted in years*. A simulation year is only 0.4946
# Earth years, so representing the same absolute span needs more of them: the
# reciprocal, 2.022.
#
# This is the same trap gdd5min was, one level up, and it is easy to miss because
# "500 years of spin-up" reads like an absolute statement and is not. At
# nyear_spinup 500 this world gets 247 Earth years of soil and vegetation
# development where Earth practice assumes 500.
#
#   nyear_spinup   time for vegetation and soil pools to reach steady state
#   distinterval   mean return time of generic patch-destroying disturbance
#   freenyears     time allowed to build an N pool before N limitation bites
SCALED_UP = ("nyear_spinup", "distinterval", "freenyears")

# Named so the decision not to scale them is explicit and reviewable, rather
# than an omission someone later reads as an oversight.
DELIBERATELY_UNSCALED = (
    "phengdd5ramp",   # within-season accumulation, already absolute time
    "tcmin_surv", "tcmin_est", "tcmax_est", "twmin_est", "twminusc",
    "gdd0_min", "gdd0_max",  # "no restriction" sentinels, 0 and 100000
    # Counted in growing seasons rather than in absolute time. One simulation
    # year is one seasonal cycle on this world just as on Earth, so these are
    # already in the right unit.
    "estinterval",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=GUESS_SOURCE / "data" / "ins" / "global.ins")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    orbital_days = orbit.orbital_year_days(config)
    factor = orbital_days / EARTH_YEAR_DAYS

    text = args.source.read_text()
    changes = []

    def rescale(match: re.Match) -> str:
        name, gap, value = match.group(1), match.group(2), float(match.group(3))
        # Sentinels meaning "no limit" stay sentinels; scaling 0 is a no-op but
        # scaling a large "no restriction" value would quietly become a limit.
        if value <= 0.0 or value >= 1e4:
            return match.group(0)
        new = value * factor
        changes.append({"parameter": name, "from": value, "to": round(new, 2),
                        "direction": "down, an annual sum"})
        return f"{name}{gap}{new:.1f}"

    def rescale_up(match: re.Match) -> str:
        name, gap, value = match.group(1), match.group(2), float(match.group(3))
        if value <= 0.0:
            return match.group(0)
        new = value / factor
        changes.append({"parameter": name, "from": value, "to": round(new, 2),
                        "direction": "up, a duration in years"})
        return f"{name}{gap}{new:.0f}"

    rescaled = re.compile(
        r"\b(" + "|".join(SCALED_DOWN) + r")(\s+)([0-9.]+)").sub(rescale, text)
    rescaled = re.compile(
        r"\b(" + "|".join(SCALED_UP) + r")(\s+)([0-9.]+)").sub(rescale_up, rescaled)

    rescaled_names = ", ".join(sorted({c["parameter"] for c in changes})) or "nothing"
    header = f"""!///////////////////////////////////////////////////////////////////////////////
!// GENERATED by biosphere/scripts/build_vesper_pfts.py. Do not edit.
!//
!// {args.source.name} with annual degree-day limits rescaled for Vesper's year.
!//
!// orbit        {orbital_days:.4f} Earth days at {config['orbit']['baseline_flux_earth']} S-Earth
!// Earth year   {EARTH_YEAR_DAYS:.4f} days
!// factor       {factor:.6f}
!// rescaled     {rescaled_names}
!// NOT rescaled {", ".join(DELIBERATELY_UNSCALED)}
!//
!// phengdd5ramp is deliberately untouched: it is a within-season accumulation,
!// already in absolute time, and scaling it would be a real error.
!//
!// generated {datetime.now(timezone.utc).isoformat(timespec="seconds")}
!///////////////////////////////////////////////////////////////////////////////

"""
    output = args.output or (GENERATED / "vesper_pfts.ins")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(header + rescaled)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": str(args.source),
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "orbital_year_earth_days": orbital_days,
        "earth_year_days": EARTH_YEAR_DAYS,
        "scale_factor": factor,
        "scaled_down_annual_sums": SCALED_DOWN,
        "scaled_up_year_counts": SCALED_UP,
        "deliberately_unscaled": DELIBERATELY_UNSCALED,
        "changes": changes,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    report["generator"] = "biosphere/scripts/build_vesper_pfts.py"
    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"orbit   {orbital_days:.4f} Earth days -> factor {factor:.6f}")
    print(f"changed {len(changes)} parameter values:")
    for change in changes:
        print(f"   {change['parameter']:14s} {change['from']:>8.1f} -> {change['to']:.1f}")
    print(f"\nwrote {output.relative_to(PROJECT_ROOT)}")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
