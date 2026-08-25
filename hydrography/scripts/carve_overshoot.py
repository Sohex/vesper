#!/usr/bin/env python3
"""How much of a carve the climate it produced would not have made.

    python hydrography/scripts/carve_overshoot.py \
        --applied      hydrography/data/<carved>/carve_list.json \
        --reevaluated  hydrography/analysis/carve_list_postcarve.json
    python hydrography/scripts/carve_overshoot.py --self-test

Worldbuilding. Vesper is an invented planet and this compares two verdicts on
the modelled drainage of its terrain.

## What the overshoot is, and why it is not the bracket

`docs/src/pipeline/loops.md`: the verdict map is ANTITONE. Carving removes the
bright closed-basin fill, the modelled land darkens, the modelled world warms,
open-water evaporation rises, and basins that were marginal would now stay
closed. So a larger carve set produces a smaller next verdict, and the set a
pass carves is an upper bound on the set that same terrain's own climate would
carve. The overshoot is that difference, counted in basins.

Loop A exits by carving the INTERSECTION of the verdicts taken at the two
bounding climates, and that intersection is a construction over ONE axis: the
vegetation state, which is unknown when the verdict is taken. It is not a bound
on the antitone feedback, and cannot be, as set algebra: both arms run on
PRE-CARVE terrain, so both are cold relative to the world their own carve
produces, and an overshooting basin is by definition one the intersection cut,
which puts it outside the bracketed set every time. The bracket and the
overshoot are two different uncertainties and the bracket measures one of them.

That is why this number is reported beside the bracket width rather than inside
it. A residual of a few basins against a bracket of hundreds says the exit
predicate is sound in practice; a residual of the same order as the bracket says
the loop was exited one pass early.

## The instrument, and the route that does not work

The obvious route is to re-run the verdict on the carved build. It cannot work.
A carved basin's rim has been breached, so on the finished terrain there is no
impoundment left to measure a catchment-to-spill-area ratio against, and it is
absent from that build's `basins.nc` entirely. `export_carve_list.py --previous`
says as much when it carries such an entry forward at retain 0: its verdict is
not re-decidable there.

The route that works holds the GEOMETRY still and moves only the CLIMATE. Take
the verdict again on the pre-carve build's own basins and coupling matrix, and
the carved build's baseline climatology:

    python hydrography/scripts/export_carve_list.py \
        --basins       hydrography/data/<precarve>/basins.nc \
        --coupling     hydrography/data/<precarve>/coupling_exoplasim-<rung>.nc \
        --climatology  <the carved build's baseline regular climatology> \
        --out-json     hydrography/analysis/carve_list_postcarve.json

which is a DELIBERATE cross-build read in the sense of CLAUDE.md rule 5, and is
the only combination that isolates the feedback: every term except the climate
is the one the applied verdict used. Two cautions on taking it. The endmember
arm has to be re-taken the same way if the applied verdict was an intersection,
or the two are not the same kind of object. And `export_carve_list.py` reads
background land albedo from `exoplasim/inputs/<rung>/`, which is staged per
rung and NOT per build, so it holds one build's albedo at a time: taking this
measurement is correct while the carved build is staged and silently wrong once
anything else is.

## The controls

`--self-test` runs four, on fixtures rather than on a real verdict, so the
comparison can fail on itself:

  IDENTITY     A verdict against itself is zero overshoot and zero undershoot.
               A join that lost or mispaired ids fails this.
  PERMUTATION  A verdict against the same verdict with its retains permuted
               must report flips. A comparison that never compares passes
               IDENTITY and fails this.
  COVERAGE     A re-evaluation missing one of the applied carves is refused.
               Scoring it as "did not flip" would report an overshoot of zero
               for the one case where the answer is unknown.
  SIGN         A re-evaluation that preserves everything makes the whole applied
               carve set overshoot; one that carves everything makes it zero.
"""
from __future__ import annotations

import argparse
import datetime
import json
import subprocess
from pathlib import Path

from _paths import ANALYSIS, PROJECT_ROOT  # noqa: F401


def verdicts(sidecar: dict) -> tuple[dict[str, float], set[str]]:
    """`{id: retain}` for the basins a pass DECIDED, plus the ids it carried.

    Retain is the source of truth rather than the `verdict` string, because
    `export_carve_list.py` derives the string from retain and not the reverse.
    A carried entry was carved by an earlier pass and re-decided by none, so it
    is separated here instead of being counted as a decision.
    """
    decided: dict[str, float] = {}
    carried: set[str] = set()
    for entry in sidecar["basins"]:
        if entry.get("carried_from_previous_pass"):
            carried.add(entry["id"])
            continue
        decided[entry["id"]] = float(entry["retain"])
    return decided, carried


def compare(applied: dict, reevaluated: dict,
            area_at_spill: dict[str, float] | None = None) -> dict:
    """The overshoot, its opposite, and what the bracket said."""
    applied_retain, applied_carried = verdicts(applied)
    new_retain, _ = verdicts(reevaluated)

    cut = {b for b, r in applied_retain.items() if r <= 0.0}
    kept = set(applied_retain) - cut

    # COVERAGE. A basin the applied pass cut and the re-evaluation does not
    # carry is a basin whose flip is UNKNOWN, and folding it in as "did not
    # flip" would report the wrong answer in exactly the case that matters.
    uncovered = sorted(cut - set(new_retain))
    if uncovered:
        raise SystemExit(
            f"{len(uncovered)} basins the applied verdict carved are absent "
            f"from the re-evaluation, first {uncovered[:3]}. Re-evaluate the "
            "PRE-CARVE build's basins against the post-carve climatology; the "
            "carved build no longer has the geometry to decide them.")

    still_cut = {b for b in cut if new_retain[b] <= 0.0}
    overshoot = sorted(cut - still_cut)
    undershoot = sorted(b for b in kept
                        if b in new_retain and new_retain[b] <= 0.0)

    bracket = applied.get("intersection") or {}
    bracketed = None if bracket.get("single_climate", True) else bracket.get("bracketed")

    out = {
        "applied": {
            "decided": len(applied_retain),
            "carried_from_earlier_passes": len(applied_carried),
            "carved": len(cut),
            "kept": len(kept),
        },
        "overshoot": {
            "basins": len(overshoot),
            "fraction_of_applied_carve": len(overshoot) / max(len(cut), 1),
            "ids": overshoot,
            "what": "carved by the applied verdict, and not carved by the "
                    "climate that carve produced. One-signed by the antitone "
                    "argument, and the honest measure of what the pass cost.",
        },
        "undershoot": {
            "basins": len(undershoot),
            "ids": undershoot,
            "what": "kept by the applied verdict and cut by the re-evaluation. "
                    "The antitone argument does not predict these, so a large "
                    "count is evidence that something other than the carve "
                    "feedback moved between the two climates.",
        },
        "against_the_bracket": {
            "bracketed_basins": bracketed,
            "overshoot_over_bracketed": (
                None if not bracketed else len(overshoot) / bracketed),
            "note": "The bracketed set is the two bounding climates' "
                    "disagreement, which is the VEGETATION uncertainty. An "
                    "overshooting basin was cut by both arms, so it is never in "
                    "the bracketed set: these are two uncertainties and the "
                    "exit predicate measures one. The ratio says whether the "
                    "unmeasured one is small beside the measured one.",
        },
        "not_re_decidable": {
            "basins": len(applied_carried),
            "why": "carved by a pass earlier than the applied one, so no "
                   "terrain in this lineage still holds their geometry.",
        },
    }
    if area_at_spill is not None:
        cut_area = sum(area_at_spill.get(b, 0.0) for b in cut)
        out["overshoot"]["area_at_spill_km2"] = sum(
            area_at_spill.get(b, 0.0) for b in overshoot)
        out["overshoot"]["area_share_of_applied_carve"] = (
            out["overshoot"]["area_at_spill_km2"] / cut_area if cut_area else None)
        out["overshoot"]["area_note"] = (
            "area_at_spill_km2 is the depression FOOTPRINT and is extensive, so "
            "it is summed. It is not the catchment, which is several times "
            "larger and is the quantity the criterion integrates over.")
    return out


def _fixture(retains: dict[str, float], bracketed: int | None = None) -> dict:
    return {
        "basins": [{"id": b, "retain": r,
                    "verdict": "carve" if r <= 0 else
                               ("preserve" if r >= 1 else "marginal")}
                   for b, r in retains.items()],
        "intersection": ({"single_climate": True} if bracketed is None
                         else {"single_climate": False, "bracketed": bracketed}),
    }


def self_test() -> int:
    base = {"a": 0.0, "b": 0.0, "c": 0.5, "d": 1.0}
    failures = []

    identity = compare(_fixture(base), _fixture(base))
    if identity["overshoot"]["basins"] or identity["undershoot"]["basins"]:
        failures.append("IDENTITY: a verdict against itself reported a flip")

    permuted = compare(_fixture(base), _fixture({"a": 1.0, "b": 0.5,
                                                 "c": 0.0, "d": 0.0}))
    if permuted["overshoot"]["basins"] != 2 or permuted["undershoot"]["basins"] != 2:
        failures.append(
            f"PERMUTATION: expected 2 and 2, got "
            f"{permuted['overshoot']['basins']} and "
            f"{permuted['undershoot']['basins']}")

    try:
        compare(_fixture(base), _fixture({"a": 0.0, "c": 0.5, "d": 1.0}))
    except SystemExit:
        pass
    else:
        failures.append("COVERAGE: a re-evaluation missing a carved basin was "
                        "accepted instead of refused")

    all_kept = compare(_fixture(base), _fixture({k: 1.0 for k in base}))
    if all_kept["overshoot"]["basins"] != 2:
        failures.append(
            f"SIGN: preserving everything should make the whole applied carve "
            f"set overshoot, got {all_kept['overshoot']['basins']} of 2")
    all_cut = compare(_fixture(base), _fixture({k: 0.0 for k in base}))
    if all_cut["overshoot"]["basins"] != 0 or all_cut["undershoot"]["basins"] != 2:
        failures.append(
            f"SIGN: carving everything should be zero overshoot and 2 "
            f"undershoot, got {all_cut['overshoot']['basins']} and "
            f"{all_cut['undershoot']['basins']}")

    for f in failures:
        print(f"FAIL  {f}")
    if failures:
        return 1
    print("ok  identity, permutation, coverage and sign")
    return 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--applied", type=Path, default=None,
                    help="carve_list.json of the pass whose carve was applied")
    ap.add_argument("--reevaluated", type=Path, default=None,
                    help="carve_list.json from re-taking that verdict on the "
                         "SAME basins under the climatology the carve produced")
    ap.add_argument("--basins", type=Path, default=None,
                    help="the pre-carve basins.nc, to price the overshoot in "
                         "depression footprint as well as in basins")
    ap.add_argument("--output", type=Path,
                    default=ANALYSIS / "carve_overshoot.json")
    ap.add_argument("--self-test", action="store_true",
                    help="run the four controls on fixtures and exit")
    args = ap.parse_args()

    if args.self_test:
        raise SystemExit(self_test())
    if args.applied is None or args.reevaluated is None:
        raise SystemExit("--applied and --reevaluated are both required; "
                         "--self-test needs neither")

    applied = json.loads(args.applied.read_text(encoding="utf-8"))
    reevaluated = json.loads(args.reevaluated.read_text(encoding="utf-8"))

    area = None
    if args.basins is not None:
        from netCDF4 import Dataset
        with Dataset(args.basins) as ds:
            area = {str(b): float(a) for b, a in
                    zip(ds["basin_id"][:], ds["area_at_spill_km2"][:],
                        strict=True)}

    result = compare(applied, reevaluated, area)
    result["inputs"] = {
        "applied": str(args.applied),
        "reevaluated": str(args.reevaluated),
        "basins": str(args.basins) if args.basins else None,
    }
    result["generated"] = datetime.datetime.now(datetime.UTC).isoformat()
    result["git_commit"] = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
        capture_output=True, text=True).stdout.strip() or None

    o = result["overshoot"]
    print(f"applied carve  {result['applied']['carved']:5d} basins")
    print(f"overshoot      {o['basins']:5d} "
          f"({o['fraction_of_applied_carve']:.1%} of it) would no longer carve")
    print(f"undershoot     {result['undershoot']['basins']:5d} kept basins "
          "would now carve")
    b = result["against_the_bracket"]["bracketed_basins"]
    if b:
        print(f"bracketed      {b:5d} basins, so the overshoot is "
              f"{o['basins'] / b:.2f} of the bracket")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
