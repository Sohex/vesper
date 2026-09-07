#!/usr/bin/env python3
"""How many basins a hydrological-cycle suppression is worth. DUST-11, before the run.

    python aeolian/scripts/dust_runoff_sensitivity.py

This exists to make DUST-11's prediction FALSIFIABLE before the climate run that
tests it exists. It converts "dust lowers the simulated rainfall by x percent" into
the currency the carve is decided in, by perturbing the baseline climatology's
precipitation and evaporation and re-taking the verdict on the result.

## Why the existing equivalence needed measuring

`notes/dust.md` asserted that -10% on catchment runoff is worth about the same as
+10% on lake evaporation, which is 144 basins. That was an assumption about a
ratio, not a measurement, and it is wrong: the aridity index is

    (E_lake - P) / R

and its SIGN varies across the population -- the median basin on this world has
P > E and therefore a negative index. Scaling the denominator alone pushes the
positive half toward surviving and the negative half toward carving, and the two
largely cancel. Scaling precipitation moves the numerator as well, and in the
same direction for both halves, so it does not cancel.

Both curves are computed here, because the difference between them is the point:
one of them is what dust does and the other is what the shorthand assumed.

## What it is not

It is not a climate. The radiation is untouched, so this says nothing about the
lake limb DUST-10 measured, and the two must not be added. The verdict on a
world with dust in it is the one taken on a dust climatology, which is the run
`aeolian/notes/prescribed-dust-run.md` specifies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "lib") not in sys.path:
    sys.path.insert(0, str(ROOT / "lib"))
from paths import climatology_path, rel  # noqa: E402

VERDICT = ROOT / "hydrography" / "scripts" / "carve_verdict.py"
OUT = ROOT / "aeolian" / "analysis" / "dust_runoff_sensitivity.json"

# The precipitation reductions to span, and how much of each reduction
# evaporation takes with it. Both come from `notes/dust.md`: the fast adjustment
# is of order 10% here, and the three following fractions are the same three that
# note's own runoff table uses.
PRECIP_REDUCTIONS = (0.05, 0.10, 0.15)
EVAP_FOLLOWS = (1.00, 0.80, 0.50)
RUNOFF_SCALES = (0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10)


def run_verdict(work: Path, climatology: Path, runoff_scale: float = 1.0) -> dict:
    out = work / f"verdict_{abs(hash((str(climatology), runoff_scale))):x}.json"
    cmd = [sys.executable, str(VERDICT), "--climatology", str(climatology),
           "--output", str(out)]
    if runoff_scale != 1.0:
        cmd += ["--runoff-scale", str(runoff_scale)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        raise SystemExit(f"carve_verdict failed:\n{r.stdout}\n{r.stderr}")
    return json.loads(out.read_text(encoding="utf-8"))


def perturbed(src: Path, dst: Path, precip_scale: float, evap_scale: float) -> None:
    """Copy a climatology with pr and evap rescaled and nothing else touched."""
    import netCDF4 as nc
    shutil.copyfile(src, dst)
    with nc.Dataset(dst, "a") as ds:
        for name, scale in (("pr", precip_scale), ("evap", evap_scale)):
            v = ds[name]
            v[:] = np.array(v[:], dtype=float, copy=True) * scale


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()
    # PINNED TO THE BASELINE DELIBERATELY, and not through
    # `best_available_climatology` as the aerosol producers are. This is a
    # sensitivity study, so its reference has to be a NAMED climatology that
    # does not change stage underneath a result: a verdict perturbed against
    # the bootstrap on one pass and the baseline on the next is two studies
    # reported as one. The docstring's "the baseline climatology's
    # precipitation and evaporation" is that choice, and this is where it is
    # made. `aeolian/scripts/build_dust_source_fields.py` states the opposite
    # exception for the opposite reason.
    if args.climatology is None:
        args.climatology = climatology_path()

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        base = run_verdict(work, args.climatology)
        base_carve = base["bounds"]["penman"]["basins_carved"]
        base_runoff = base["runoff_source"]["catchment_mean_mm_per_year"]["p_minus_e"]
        print(f"baseline: {base_carve} basins carved, catchment runoff "
              f"{base_runoff:.2f} mm/yr\n")

        denominator = []
        print("  denominator only (runoff scaled, P and E left alone)")
        print(f"  {'dR%':>6} {'carved':>7} {'delta':>7}")
        for s in RUNOFF_SCALES:
            v = run_verdict(work, args.climatology, runoff_scale=s)
            c = v["bounds"]["penman"]["basins_carved"]
            denominator.append({"runoff_scale": s, "basins_carved": c,
                                "delta": c - base_carve})
            print(f"  {100 * (s - 1):+6.0f} {c:7d} {c - base_carve:+7d}")

        cycle = []
        print("\n  hydrological cycle suppressed (P and E both reduced)")
        print(f"  {'dP%':>6} {'E/P':>5} {'dR%':>7} {'carved':>7} {'delta':>7}")
        for dp in PRECIP_REDUCTIONS:
            for follow in EVAP_FOLLOWS:
                pert = work / f"clim_{dp}_{follow}.nc"
                perturbed(args.climatology, pert, 1.0 - dp, 1.0 - dp * follow)
                v = run_verdict(work, pert)
                c = v["bounds"]["penman"]["basins_carved"]
                r = v["runoff_source"]["catchment_mean_mm_per_year"]["p_minus_e"]
                cycle.append({
                    "precip_reduction": dp, "evaporation_follows": follow,
                    "catchment_runoff_mm_per_year": r,
                    "runoff_change": r / base_runoff - 1.0,
                    "basins_carved": c, "delta": c - base_carve})
                print(f"  {-100 * dp:+6.0f} {follow:5.2f} "
                      f"{100 * (r / base_runoff - 1):+7.1f} {c:7d} "
                      f"{c - base_carve:+7d}")
                pert.unlink()

    payload = {
        "note": "DUST-11's prediction, made before the climate run that tests it. "
                "Two sensitivities on the baseline climatology: scaling the "
                "criterion's DENOMINATOR alone, and suppressing the whole "
                "hydrological cycle. The radiation is untouched in both, so "
                "neither contains the lake limb DUST-10 measured and they must "
                "not be added to it. Generated by "
                "aeolian/scripts/dust_runoff_sensitivity.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatology": rel(args.climatology),
        "climatology_sha256": hashlib.sha256(
            args.climatology.read_bytes()).hexdigest(),
        "baseline_basins_carved": base_carve,
        "baseline_catchment_runoff_mm_per_year": base_runoff,
        "denominator_only": denominator,
        "hydrological_cycle": cycle,
        "why_the_two_differ":
            "The aridity index is (E_lake - P)/R and its sign varies across the "
            "population: the median basin here has P > E, so the index is "
            "negative. Scaling R alone makes a positive index larger and a "
            "negative one more negative, which moves the two halves toward "
            "opposite verdicts and largely cancels. Reducing P raises the "
            "numerator for every basin as well as cutting R, so it does not.",
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=ROOT).stdout.strip() or None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {rel(args.output)}")


if __name__ == "__main__":
    main()
