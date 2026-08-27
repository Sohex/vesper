#!/usr/bin/env python3
"""The finalizer: is a sequence of individually converged loops jointly converged?

    python scripts/verify_joint_convergence.py
    python scripts/verify_joint_convergence.py --self-test

Worldbuilding. Vesper is an invented planet, and every artifact named here
describes the modelled state of it.

`docs/src/pipeline/loops.md`, final section. The pipeline's loops are nested and
loop D declares that it advances only after "replaying A, B and C" on the new
support. The ladder's cost argument does not buy that replay: the bulk of the
work sits at T21 and the minimum sufficient at the rungs above it, so the final
state carries terrain carved at T21, soil and vegetation converged at T21, and a
climate merely SETTLED at the operating support. Each loop exited. None of them
exited against the others' final state, and a sequence of individually converged
loops is not a jointly converged system.

This re-evaluates each loop's OWN exit predicate against the final state and
reports whether it still holds. It is a VERIFICATION and not another turn of any
loop: nothing iterates on its verdict, because re-entering A at the operating
support is a commissioning-scale purchase and that decision is the author's.

## The four predicates, each taken from the loop rather than invented

**A, the carve intersection.** Loop A's exit is a CONSTRUCTION rather than a
tolerance. The verdict map is antitone, so successive verdicts bracket instead
of converging, and the loop exits by carving the INTERSECTION of the verdicts
taken at the two bounding climates. Re-taken on the operating support's baseline
climatology, that intersection must be the set that was actually carved. Two
things are asked, and they fail differently: the re-taken list's own carved set
must BE the intersection of its own two arms, and that intersection must be the
applied set. The first catches a construction that stopped being an
intersection; the second catches a carve the operating support's climate would
not have made.

**B, the soil and biosphere loop.** The criteria in
`pedology/config/pedogenesis.yaml`, applied to the soil the final state carries
against the iteration before it: land-mean soil organic carbon under
`convergence.soil_carbon_relative_tolerance`, the share of land cells whose clay
moved under `convergence.texture_cells_moved_tolerance`, and the iteration count
within `convergence.maximum_iterations`.

**C, the vegetation-climate loop.** Its exit ALREADY is a re-take: the verdict on
MODELLED vegetation, and whether basins flip. Generalised here to the operating
support rather than the support it was first taken on. A flip in either
direction fails and names loop A, because that is what C's exit says to do
about one.

**D, the declared escalation.** The route's own invariants, asked of `lib/rungs.py`
rather than restated: the operating support is reached, and every change of rung
happened at CONSTANT dt. The first is read off the run chain that produced the
baseline climatology; the second is read off the same chain hop by hop, so it is
a statement about the runs that exist and not about the route that was declared.
`rungs._check_route()` covers the declared route and runs at import.

## NOT EVALUABLE is a verdict, and it is not a pass

Three statuses, and the third is why this file is worth having. On a tree with
no baseline climatology at the operating support there is nothing for A, B or C
to compare, and a finalizer that reports PASS when it could not test anything is
worse than one that refuses: it converts an absence of evidence into evidence.
Every refusal names what was missing and the command that would produce it.

The same rule covers the near miss, which is the one that would be silent: an
artifact at a COARSER support than the operating one is NOT a substitute for the
operating support's. A T21 baseline against a T85 operating support is a refusal
with the two rungs named, never a pass, because the whole point of the finalizer
is that a verdict taken on a coarser climate looks exactly like one taken on the
finer.

Every verdict carries WHAT IT COMPARED under `compared`: which climatology,
which support, which carve list, which soil, which runs. A pass with an empty
`compared` is a vacuous pass and a reader can see it as one.

## Exit status

0 when all four predicates hold, 1 when any of them fails, 2 when none failed
and at least one could not be evaluated. Two rather than zero, because a caller
that treats a refusal as a pass is the failure this file exists to prevent.

## Checking the checker

    python scripts/verify_joint_convergence.py --self-test

Each predicate is driven to FAIL on a fixture rather than asserted, and each
failing case is paired with the passing case that proves the fixture is what
moved it. The A case is the interesting one: a carved set that is NOT the
intersection of the two arms it claims to be built from, which the check must
catch, and an applied set that is not the re-taken intersection, which it must
also catch. A finalizer whose checks cannot fail is exactly the thing it exists
to prevent.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

import yaml                                  # noqa: E402

import rungs                                 # noqa: E402
from paths import rel                        # noqa: E402

PASS, FAIL, NOT_EVALUABLE = "pass", "FAIL", "not evaluable"

REPORT = ROOT / "analysis" / "joint_convergence.json"

# The share of a land cell's clay that counts as having MOVED between soil
# iterations. `pedology/config/pedogenesis.yaml` states it in the comment above
# `convergence.texture_cells_moved_tolerance` -- "fraction of land cells whose
# clay content moves by more than 0.02" -- and carries no key for it, so the
# tolerance has a key and the threshold it is a tolerance ON does not. Read from
# the config the moment a key exists; until then this is the config's own number
# and not a choice made here.
CLAY_MOVE_THRESHOLD = 0.02

# The two bounding climates of loop A's bracket, by `model.land_albedo_source`.
# `docs/src/pipeline/loops.md`: the arms differ in this and in nothing else, the
# cold arm is a BOUND rather than a world, and the warm arm is whichever
# vegetated surface the pass ran on.
COLD_ARM = "lithology"
WARM_ARMS = ("vegetated", "modelled")
# Loop C's re-take is on the MODELLED biosphere specifically. An assumed canopy
# is the arm loop A already used, so a verdict on one is not a re-take.
MODELLED_ARM = "modelled"


@dataclass
class Verdict:
    """One loop's exit predicate, re-evaluated. `compared` is not optional.

    A verdict with nothing under `compared` is a verdict about nothing, and the
    field exists so that a reader can tell a real pass from a vacuous one
    without going back to the tree.
    """
    loop: str
    predicate: str
    status: str
    detail: str
    compared: dict = field(default_factory=dict)

    def as_json(self) -> dict:
        return {"loop": self.loop, "predicate": self.predicate,
                "status": self.status, "detail": self.detail,
                "compared": self.compared}


# ---------------------------------------------------------------------------
# reading the tree
# ---------------------------------------------------------------------------

def read_json(path: Path):
    """The file's content, or None. Absence is a verdict here, not an error."""
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def climatology_rung(path: Path) -> str | None:
    """The rung a climatology is ON, from its own grid.

    Taken from the file's dimensions and never from a sidecar's `resolution`,
    which several products copy out of `config/planet.yaml` at the moment they
    run. `export_carve_list.py` does exactly that, so a re-take driven by an
    explicit `--climatology` records the CONFIGURED rung beside a climatology
    that may be on another one. The grid cannot lie about itself.
    """
    try:
        from netCDF4 import Dataset
    except ImportError:                                   # pragma: no cover
        return None
    try:
        with Dataset(path) as ds:
            nlat = len(ds.dimensions["lat"])
    except (OSError, KeyError):
        return None
    try:
        return rungs.rung_of_latitudes(nlat)
    except Exception:
        return None


def climatology_run_id(path: Path) -> str | None:
    """Which run a climatology was averaged from, off its own attributes."""
    try:
        from netCDF4 import Dataset
    except ImportError:                                   # pragma: no cover
        return None
    try:
        with Dataset(path) as ds:
            return getattr(ds, "vesper_run_id", None)
    except OSError:
        return None


def run_manifest(run_id: str, root: Path = ROOT) -> dict | None:
    """A run's manifest, live or archived.

    `scripts/archive_runs.py` keeps `run_manifest.json` under
    `archive/runs/<dir>/` when it deletes the output, so an archived run can
    still answer what it integrated. A run that is in neither place cannot, and
    that is a refusal rather than a guess.
    """
    for base in (root / "exoplasim" / "runs", root / "archive" / "runs"):
        got = read_json(base / run_id / "run_manifest.json")
        if got is not None:
            return got
    return None


def arm_of(climatology: Path, root: Path = ROOT) -> dict:
    """What a climatology IS: its rung, its run, and that run's land surface.

    The land surface is `model.land_albedo_source` off the run's own manifest,
    because it is the one thing that separates loop A's two bounding arms from
    each other and either of them from loop C's re-take, and no climatology
    carries it.
    """
    path = Path(climatology)
    out = {"path": rel(path), "exists": path.is_file(),
           "rung": None, "run_id": None, "land_albedo_source": None}
    if not out["exists"]:
        return out
    out["rung"] = climatology_rung(path)
    out["run_id"] = climatology_run_id(path)
    if out["run_id"]:
        manifest = run_manifest(out["run_id"], root)
        if manifest:
            model = (manifest.get("source_config") or {}).get("model") or {}
            out["land_albedo_source"] = model.get("land_albedo_source")
    return out


# The tree's top-level directories, for re-rooting a path a sidecar recorded as
# absolute in another checkout. Not every directory, only the ones an artifact
# path can start with.
_TOP_LEVEL = ("exoplasim", "hydrography", "pedology", "biosphere", "aeolian",
              "minerals", "analysis", "source", "maps", "config", "archive")


def resolve_under_root(path_text: str, root: Path = ROOT) -> Path:
    """A path recorded in a sidecar, resolved against THIS tree.

    Sidecars record whatever was on the command line, and several products
    record `str(args.climatology)`, which is absolute and belongs to whichever
    checkout ran them. A worktree reading a main checkout's report would
    otherwise resolve to the other tree's file and compare the wrong world. The
    path is re-rooted at its last top-level component, and returned untouched
    when it has none, so a genuinely external path is still reported as itself.
    """
    p = Path(path_text)
    if not p.is_absolute():
        return root / p
    parts = p.parts
    candidates = [root.joinpath(*parts[i:]) for i, part in enumerate(parts)
                  if part in _TOP_LEVEL]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    # Leftmost, not rightmost: `exoplasim/analysis/...` contains two top-level
    # names and only the first one starts the path this tree knows.
    return candidates[0] if candidates else p


def carved_ids(sidecar: dict) -> set[str]:
    """The basins a carve list actually carves.

    Retain, not the overflow test: retain is what Orogen cuts with, so a basin
    that overflows and keeps most of its rim is a marginal landform rather than
    a carve. `export_carve_list.py` says the same thing where it names the
    verdict.
    """
    out = set()
    for row in sidecar.get("basins") or []:
        if row.get("verdict") == "carve" or float(row.get("retain", 1.0)) <= 0.0:
            out.add(str(row["id"]))
    return out


def arm_carve_sets(sidecar: dict) -> tuple[set[str], set[str]] | None:
    """The two bounding arms' carve sets, per basin, off a carve list sidecar.

    None when the list was taken on one climate: `retain_warm_vegetated_arm` and
    `retain_cold_bare_rock_arm` are null on every row then, and
    `intersection.single_climate` says so at the top of the file.
    """
    warm, cold = set(), set()
    seen = False
    for row in sidecar.get("basins") or []:
        w, c = row.get("retain_warm_vegetated_arm"), row.get("retain_cold_bare_rock_arm")
        if w is None or c is None:
            continue
        seen = True
        if float(w) <= 0.0:
            warm.add(str(row["id"]))
        if float(c) <= 0.0:
            cold.add(str(row["id"]))
    return (warm, cold) if seen else None


def sidecar_ids(sidecar: dict) -> set[str]:
    return {str(row["id"]) for row in sidecar.get("basins") or []}


def read_soilmap(path: Path) -> dict | None:
    """A soilmap as {(lon, lat): {column: value}}, keyed by position.

    Columns are found BY NAME off the header rather than by position, so a
    column added to `build_soil.py`'s writer moves nothing here.
    """
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    if not lines:
        return None
    names = lines[0].split()
    rows = {}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(names):
            continue
        record = dict(zip(names, (float(p) for p in parts), strict=True))
        rows[(round(record["Lon"], 4), round(record["Lat"], 4))] = record
    return rows or None


# ---------------------------------------------------------------------------
# A -- the carve intersection, re-taken on the operating support
# ---------------------------------------------------------------------------

A_PREDICATE = ("carve the INTERSECTION of the verdicts taken at the two "
               "bounding climates; re-taken on the operating support's "
               "baseline climatology, that intersection is the set that was "
               "actually carved")


def loop_a(applied: dict | None, retaken: dict | None, support: str,
           warm_arm: dict | None, cold_arm: dict | None,
           applied_path: str = "", retaken_path: str = "") -> Verdict:
    """Loop A's exit, re-evaluated. Takes data so a fixture can drive it.

    Order matters. Everything that could make the comparison meaningless is
    asked first and refuses, and only a comparison that could have gone either
    way is allowed to pass.
    """
    compared = {"operating_support": support,
                "applied_carve_list": applied_path,
                "retaken_carve_list": retaken_path,
                "warm_vegetated_arm": warm_arm,
                "cold_bare_rock_arm": cold_arm}

    def refuse(why):
        return Verdict("A", A_PREDICATE, NOT_EVALUABLE, why, compared)

    if applied is None:
        return refuse(
            f"no applied carve list at {applied_path}. Loop A has not exited on "
            "this lineage, so there is no carved set to hold a re-take against.")
    if retaken is None:
        return refuse(
            f"no re-taken carve list at {retaken_path}. Produce one with "
            "hydrography/scripts/export_carve_list.py on the PRE-CARVE basins "
            "and coupling under the operating support's baseline climatology, "
            "passing --endmember-climatology for the cold arm; "
            "hydrography/scripts/carve_overshoot.py's docstring carries the "
            "invocation and the cross-build caution.")

    arms = arm_carve_sets(retaken)
    if arms is None or (retaken.get("intersection") or {}).get("single_climate"):
        return refuse(
            "the re-taken carve list was taken on ONE climate, and loop A's "
            "exit is the intersection of two. Re-run it with "
            "--endmember-climatology.")

    for label, arm, want in (("warm vegetated", warm_arm, WARM_ARMS),
                             ("cold bare-rock", cold_arm, (COLD_ARM,))):
        if arm is None or not arm.get("exists"):
            return refuse(f"the {label} arm's climatology is not in this tree, "
                          "so what the re-take ran on cannot be established")
        if arm.get("rung") != support:
            return refuse(
                f"the {label} arm's climatology is at {arm.get('rung')} and the "
                f"operating support is {support}. A coarser climatology is not "
                "a substitute for the operating support's; the whole failure "
                "this checks for is a verdict taken on a coarser climate "
                "looking exactly like one taken on the finer.")
        if arm.get("land_albedo_source") not in want:
            return refuse(
                f"the {label} arm ran with land_albedo_source "
                f"{arm.get('land_albedo_source')!r}, and loop A's bracket is "
                f"over {' or '.join(want)}. The two arms differ in that key and "
                "in nothing else.")

    a_hash = applied.get("terrain_hash")
    r_hash = retaken.get("terrain_hash")
    compared["terrain_hash"] = {"applied": a_hash, "retaken": r_hash}
    if a_hash and r_hash and a_hash != r_hash:
        return refuse(
            f"the applied list is on terrain {str(a_hash)[:16]} and the re-take "
            f"on {str(r_hash)[:16]}. The re-take holds the GEOMETRY still and "
            "moves only the climate, so a differing terrain means it measured "
            "the next pass instead of this one.")

    applied_set = carved_ids(applied)
    warm, cold = arms
    intersection = warm & cold
    present = sidecar_ids(retaken)
    compared["counts"] = {
        "applied_carved": len(applied_set),
        "cut_by_warm_vegetated_arm": len(warm),
        "cut_by_cold_bare_rock_arm": len(cold),
        "intersection": len(intersection),
        "retaken_carved": len(carved_ids(retaken)),
        "retaken_basins": len(present),
    }

    uncovered = sorted(applied_set - present)
    if uncovered:
        return refuse(
            f"{len(uncovered)} of the {len(applied_set)} applied carves are "
            "absent from the re-take, so their verdict under the operating "
            "support's climate is unknown. Scoring an unknown as agreement "
            f"would report a pass for the one case with no answer. First: "
            f"{', '.join(uncovered[:3])}")

    # The construction, checked against itself first. A carve list's own carved
    # set is built as the larger of the two arms' retains, which is "carve only
    # what both carve" generalised to a fractional rim. If that is no longer
    # what the file holds, then nothing below is a test of loop A's exit -- it
    # is a test of some other set that the file happens to call carved.
    own = carved_ids(retaken)
    not_intersection = sorted((own ^ intersection))
    if not_intersection:
        return Verdict(
            "A", A_PREDICATE, FAIL,
            f"the re-taken list's own carved set is not the intersection of the "
            f"two arms it records: {len(own - intersection)} carved by neither "
            f"or by only one arm, {len(intersection - own)} in the intersection "
            f"and not carved. Loop A's exit IS the intersection, so a list that "
            f"carves anything else did not exit loop A. First: "
            f"{', '.join(not_intersection[:3])}", compared)

    would_not_carve = sorted(applied_set - intersection)
    would_now_carve = sorted(intersection - applied_set)
    compared["counts"]["applied_but_not_in_retaken_intersection"] = len(would_not_carve)
    compared["counts"]["in_retaken_intersection_but_not_applied"] = len(would_now_carve)
    if would_not_carve or would_now_carve:
        return Verdict(
            "A", A_PREDICATE, FAIL,
            f"the intersection re-taken at {support} is not the set that was "
            f"carved: {len(would_not_carve)} carved basins fall outside it and "
            f"{len(would_now_carve)} basins inside it were not carved. Loop A "
            f"did not exit against the final state. Re-entering it at "
            f"{support} is a commissioning-scale purchase and that decision is "
            f"not this file's to take.", compared)

    return Verdict("A", A_PREDICATE, PASS,
                   f"the intersection re-taken at {support} over "
                   f"{len(present)} basins is exactly the {len(applied_set)} "
                   f"that were carved", compared)


# ---------------------------------------------------------------------------
# B -- pedogenesis.yaml's criteria against the soil the final state carries
# ---------------------------------------------------------------------------

B_PREDICATE = ("the criteria in pedology/config/pedogenesis.yaml: land-mean "
               "soil organic carbon within soil_carbon_relative_tolerance, the "
               "share of land cells whose clay moved within "
               "texture_cells_moved_tolerance, and the iteration count within "
               "maximum_iterations")


def loop_b(criteria: dict, final: dict | None, previous: dict | None,
           iteration, soil_rung: str | None, support: str,
           final_path: str = "", previous_path: str = "") -> Verdict:
    """Loop B's exit, re-evaluated against the soil the final state carries.

    `final` and `previous` are soilmaps as `read_soilmap` returns them, so a
    fixture drives this without writing a file.
    """
    compared = {"operating_support": support,
                "soil_support": soil_rung,
                "final_soil": final_path,
                "previous_soil": previous_path,
                "iteration": iteration,
                "criteria": {
                    "soil_carbon_relative_tolerance":
                        criteria.get("soil_carbon_relative_tolerance"),
                    "texture_cells_moved_tolerance":
                        criteria.get("texture_cells_moved_tolerance"),
                    "maximum_iterations": criteria.get("maximum_iterations"),
                    "clay_move_threshold": CLAY_MOVE_THRESHOLD,
                    "clay_move_threshold_source":
                        "the comment above texture_cells_moved_tolerance in "
                        "pedology/config/pedogenesis.yaml; it carries no key",
                }}

    def refuse(why):
        return Verdict("B", B_PREDICATE, NOT_EVALUABLE, why, compared)

    if final is None:
        return refuse(f"no soil at {final_path}. Loop B's criteria are about "
                      "the soil the final state carries and there is none.")
    if previous is None:
        return refuse(
            f"no previous iteration's soil at {previous_path or '(none named)'}. "
            "Loop B's criteria are differences BETWEEN iterations, and the "
            "`soil` step writes one file per build and rung which each "
            "iteration overwrites, so the tree keeps only the latest. Name one "
            "with --previous-soil.")
    if soil_rung is not None and soil_rung != support:
        return refuse(
            f"the soil was built at {soil_rung} and the operating support is "
            f"{support}. A soil converged on a coarser climate is not the soil "
            "the final state carries.")

    shared = sorted(set(final) & set(previous))
    if not shared:
        return refuse("the two soils share no land cell, so nothing can be "
                      "differenced. They are on different grids or different "
                      "terrains.")

    prev_c = [previous[k]["soilc"] for k in shared]
    final_c = [final[k]["soilc"] for k in shared]
    mean_prev = sum(prev_c) / len(prev_c)
    mean_final = sum(final_c) / len(final_c)
    moved = sum(1 for k in shared
                if abs(final[k]["clay"] - previous[k]["clay"]) > CLAY_MOVE_THRESHOLD)
    moved_share = moved / len(shared)

    compared["measured"] = {
        "land_cells_compared": len(shared),
        "land_mean_soil_carbon_kg_m2": {"previous": round(mean_prev, 6),
                                        "final": round(mean_final, 6)},
        "texture_cells_moved_fraction": round(moved_share, 6),
    }

    problems = []
    # A previous mean of zero is iteration 0's declared initial state, and a
    # relative change against it is undefined rather than infinite. Say so; do
    # not report a pass and do not report a division.
    if mean_prev == 0.0 and mean_final == 0.0:
        carbon_change = 0.0
        compared["measured"]["soil_carbon_relative_change"] = 0.0
    elif mean_prev == 0.0:
        return refuse(
            "the previous iteration's land-mean soil carbon is zero, which is "
            "`organic.initial_soil_carbon_kg_m2` and not a soil the biosphere "
            "produced. A relative change against it is undefined, so loop B's "
            "carbon criterion has nothing to be relative to.")
    else:
        carbon_change = abs(mean_final - mean_prev) / abs(mean_prev)
        compared["measured"]["soil_carbon_relative_change"] = round(carbon_change, 6)

    carbon_tol = float(criteria["soil_carbon_relative_tolerance"])
    texture_tol = float(criteria["texture_cells_moved_tolerance"])
    if carbon_change > carbon_tol:
        problems.append(f"land-mean soil carbon moved {carbon_change:.4f} "
                        f"against a tolerance of {carbon_tol}")
    if moved_share > texture_tol:
        problems.append(f"{moved_share:.4f} of land cells moved clay by more "
                        f"than {CLAY_MOVE_THRESHOLD} against a tolerance of "
                        f"{texture_tol}")
    cap = criteria.get("maximum_iterations")
    if cap is not None and iteration is not None and int(iteration) > int(cap):
        problems.append(f"the soil is at iteration {iteration} against a cap of "
                        f"{cap}")

    if problems:
        return Verdict("B", B_PREDICATE, FAIL,
                       "; ".join(problems) + f". Measured over {len(shared)} "
                       f"land cells at {soil_rung}.", compared)
    return Verdict("B", B_PREDICATE, PASS,
                   f"soil carbon moved {carbon_change:.4f} of "
                   f"{carbon_tol} and {moved_share:.4f} of land cells moved "
                   f"clay against {texture_tol}, over {len(shared)} cells",
                   compared)


# ---------------------------------------------------------------------------
# C -- the verdict on modelled vegetation, at the operating support
# ---------------------------------------------------------------------------

C_PREDICATE = ("re-take the verdict on MODELLED vegetation; if basins flip, "
               "re-enter A. Generalised to the operating support rather than "
               "the support the verdict was first taken on")


def loop_c(applied_verdict: dict | None, modelled_verdict: dict | None,
           applied_arm: dict | None, modelled_arm: dict | None, support: str,
           applied_path: str = "", modelled_path: str = "") -> Verdict:
    """Loop C's exit, re-evaluated at the operating support.

    Verdict against verdict, on `carve_list_penman` both sides. Not verdict
    against carve list: a carve list's carved set is retain, which is the
    overflow test AND the incision, so comparing one with the other would
    report every marginal landform as a flip.
    """
    compared = {"operating_support": support,
                "verdict_loop_a_used": applied_path,
                "verdict_on_modelled_vegetation": modelled_path,
                "assumed_biosphere_arm": applied_arm,
                "modelled_biosphere_arm": modelled_arm}

    def refuse(why):
        return Verdict("C", C_PREDICATE, NOT_EVALUABLE, why, compared)

    if applied_verdict is None:
        return refuse(f"no verdict at {applied_path} to hold the re-take "
                      "against, so nothing can be said to have flipped")
    if modelled_verdict is None:
        return refuse(
            f"no verdict on modelled vegetation at {modelled_path}. Produce one "
            "with hydrography/scripts/carve_verdict.py against a baseline "
            "climatology whose run staged surface albedo in `modelled` mode, "
            "which needs an LPJ-GUESS run at the operating support.")
    if modelled_arm is None or not modelled_arm.get("exists"):
        return refuse("the re-take's climatology is not in this tree, so what "
                      "it ran on cannot be established")
    if modelled_arm.get("rung") != support:
        return refuse(
            f"the re-take ran on a {modelled_arm.get('rung')} climatology and "
            f"the operating support is {support}. Loop C's exit generalised to "
            "the operating support is exactly the thing a coarser climatology "
            "cannot answer.")
    if modelled_arm.get("land_albedo_source") != MODELLED_ARM:
        return refuse(
            f"the re-take's run staged land_albedo_source "
            f"{modelled_arm.get('land_albedo_source')!r}, and loop C's exit is "
            f"the verdict on {MODELLED_ARM!r} vegetation. An assumed canopy is "
            "the arm loop A already used, so a verdict on one is not a re-take.")

    before = set(applied_verdict.get("carve_list_penman") or [])
    after = set(modelled_verdict.get("carve_list_penman") or [])
    if not before and not after:
        return refuse("neither verdict carries a `carve_list_penman`, so there "
                      "is no verdict to compare")

    # A verdict lists only the basins it CARVES, so the two must be over the same
    # catalogue or every basin the re-take never saw counts as a flip. This is
    # the coverage refusal `carve_overshoot.py` takes for the same reason: a
    # basin that was never asked did not fail to flip.
    n_before, n_after = applied_verdict.get("basins"), modelled_verdict.get("basins")
    compared["catalogue_basins"] = {"verdict_loop_a_used": n_before,
                                    "verdict_on_modelled_vegetation": n_after}
    if n_before is not None and n_after is not None and n_before != n_after:
        return refuse(
            f"the two verdicts are over catalogues of {n_before} and {n_after} "
            "basins, so a difference between their carve lists is a difference "
            "of catalogue and not a flip")

    flipped_off = sorted(before - after)
    flipped_on = sorted(after - before)
    compared["counts"] = {
        "carved_under_the_verdict_loop_a_used": len(before),
        "carved_under_modelled_vegetation": len(after),
        "no_longer_carve": len(flipped_off),
        "newly_carve": len(flipped_on),
    }
    if flipped_off or flipped_on:
        return Verdict(
            "C", C_PREDICATE, FAIL,
            f"{len(flipped_off) + len(flipped_on)} basins flip on modelled "
            f"vegetation at {support}: {len(flipped_off)} no longer carve and "
            f"{len(flipped_on)} newly carve. Loop C's exit says a flip re-enters "
            f"loop A, and that decision is the author's.", compared)
    return Verdict("C", C_PREDICATE, PASS,
                   f"no basin flips between the verdict loop A used and the "
                   f"verdict on modelled vegetation at {support}; "
                   f"{len(before)} carve under both", compared)


# ---------------------------------------------------------------------------
# D -- the route's own invariants, asked of lib/rungs.py
# ---------------------------------------------------------------------------

D_PREDICATE = ("the route's own invariants: the operating support is reached, "
               "and every change of rung happened at constant dt")


def loop_d(chain: list[dict], support: str, declared_route_error: str | None
           ) -> Verdict:
    """Loop D's exit, re-evaluated on the runs that exist.

    `chain` is the restart lineage of the run the baseline climatology was
    averaged from, oldest first, each entry `{run_id, rung, timestep_minutes,
    converted}`. Taking the chain as data rather than reading it here is what
    lets a fixture drive a route that changed the rung and the step at once,
    which no run on this tree has done.

    The DECLARED route is `lib/rungs.py`'s and is checked there: `_check_route`
    runs at import and is re-run by the caller, so this never restates the
    route. What this asks is whether the runs that exist walked it.
    """
    compared = {"operating_support": support,
                "declared_route": [list(entry) for entry in rungs.ESCALATION_ROUTE],
                "chain": chain}

    if declared_route_error:
        return Verdict("D", D_PREDICATE, FAIL,
                       f"lib/rungs.py's own route invariants do not hold: "
                       f"{declared_route_error}", compared)
    if not chain:
        return Verdict("D", D_PREDICATE, NOT_EVALUABLE,
                       "no run chain: there is no baseline climatology naming a "
                       "run, so nothing says which supports were walked",
                       compared)

    problems = []
    reached = chain[-1].get("rung")
    if reached != support:
        problems.append(f"the chain ends at {reached} and the operating support "
                        f"is {support}, so the route's terminal rung is not "
                        f"reached")

    conversions = []
    for parent, child in zip(chain, chain[1:]):
        if parent.get("rung") == child.get("rung"):
            continue
        conversions.append({
            "from": parent.get("rung"), "to": child.get("rung"),
            "from_run": parent.get("run_id"), "to_run": child.get("run_id"),
            "dt_from": parent.get("timestep_minutes"),
            "dt_to": child.get("timestep_minutes"),
            "recorded_as_conversion": bool(child.get("converted")),
        })
        if parent.get("timestep_minutes") != child.get("timestep_minutes"):
            problems.append(
                f"{child.get('run_id')} converts {parent.get('rung')} at dt "
                f"{parent.get('timestep_minutes')} to {child.get('rung')} at dt "
                f"{child.get('timestep_minutes')}. A conversion across a change "
                f"of step reinterprets the donor's stored leapfrog levels as "
                f"spanning the target's step and copies nstep, so the converted "
                f"run's calendar and stellar phase move by the ratio")
        elif not child.get("converted"):
            problems.append(
                f"{child.get('run_id')} is at {child.get('rung')} and its donor "
                f"{parent.get('run_id')} at {parent.get('rung')}, and the "
                f"manifest records no conversion. A rung change that nothing "
                f"converted read the donor's state on another grid")
    compared["rung_changes"] = conversions

    if problems:
        return Verdict("D", D_PREDICATE, FAIL, "; ".join(problems), compared)
    if not conversions:
        return Verdict(
            "D", D_PREDICATE, NOT_EVALUABLE,
            f"the chain reaches {support} in {len(chain)} runs with no change of "
            f"rung on it, so the constant-dt invariant has no conversion to "
            f"judge. A route walked in one rung is not evidence about "
            f"conversions", compared)
    return Verdict("D", D_PREDICATE, PASS,
                   f"the chain reaches {support} and all {len(conversions)} "
                   f"changes of rung happened at constant dt", compared)


# ---------------------------------------------------------------------------
# gathering the final state
# ---------------------------------------------------------------------------

def run_chain(baseline: Path, root: Path = ROOT) -> list[dict]:
    """The restart lineage behind a climatology, oldest first.

    Walks `initial_state.restart_from_run` back to the cold start, taking the
    rung and the step from each run's own `source_config.model`, which is what
    the run integrated rather than what the config says now.
    """
    run_id = climatology_run_id(baseline)
    chain: list[dict] = []
    seen: set[str] = set()
    while run_id and run_id not in seen:
        seen.add(run_id)
        manifest = run_manifest(run_id, root)
        if manifest is None:
            break
        model = (manifest.get("source_config") or {}).get("model") or {}
        state = manifest.get("initial_state") or {}
        chain.append({
            "run_id": run_id,
            "rung": (str(model["resolution"]).upper()
                     if model.get("resolution") else None),
            "timestep_minutes": model.get("timestep_minutes"),
            "converted": state.get("conversion") is not None,
        })
        run_id = state.get("restart_from_run")
        if run_id is None:
            got = state.get("restart_from")
            if got:
                match = re.search(r"(run_[0-9a-f]+)", str(got))
                run_id = match.group(1) if match else None
    chain.reverse()
    return chain


def gather(args, root: Path = ROOT) -> list[Verdict]:
    """Read the tree once, then evaluate the four predicates over what it holds."""
    config = yaml.safe_load(
        (root / "config" / "planet.yaml").read_text(encoding="utf-8"))
    build = str(config.get("source_build") or "")
    support = rungs.ESCALATION_ROUTE[-1][0]

    declared_route_error = None
    try:
        rungs._check_route()
    except RuntimeError as exc:
        declared_route_error = str(exc)

    declared_baseline = config.get("baseline_climatology")
    baseline = root / declared_baseline if declared_baseline else None

    # -- D first: everything else is a statement about the operating support,
    #    and D is what says whether it was reached.
    chain = run_chain(baseline, root) if baseline and baseline.is_file() else []
    verdicts = [loop_d(chain, support, declared_route_error)]

    # -- A
    applied_path = (Path(args.applied) if args.applied else
                    root / "hydrography" / "data" / build / "carve_list.json")
    retaken_path = (Path(args.retaken) if args.retaken else
                    root / "hydrography" / "analysis" / "carve_list_final.json")
    applied = read_json(applied_path)
    retaken = read_json(retaken_path)
    warm_arm = cold_arm = None
    if retaken:
        if retaken.get("climatology"):
            warm_arm = arm_of(resolve_under_root(retaken["climatology"], root), root)
        endmember = (retaken.get("intersection") or {}).get("endmember_climatology")
        if endmember:
            cold_arm = arm_of(resolve_under_root(endmember, root), root)
    verdicts.append(loop_a(applied, retaken, support, warm_arm, cold_arm,
                           rel(applied_path), rel(retaken_path)))

    # -- B
    soil_report = read_json(root / "pedology" / "analysis" / "soil_report.json")
    criteria = yaml.safe_load(
        (root / "pedology" / "config" / "pedogenesis.yaml").read_text(
            encoding="utf-8"))["convergence"]
    final_soil_path = (Path(args.soil) if args.soil else
                       root / "pedology" / "data" / build /
                       f"soilmap_{support}.txt")
    previous_soil_path = Path(args.previous_soil) if args.previous_soil else None
    soil_rung = None
    iteration = None
    if soil_report:
        iteration = soil_report.get("iteration")
        clim = soil_report.get("climatology")
        if clim:
            soil_rung = climatology_rung(resolve_under_root(clim, root))
    verdicts.append(loop_b(
        criteria, read_soilmap(final_soil_path),
        read_soilmap(previous_soil_path) if previous_soil_path else None,
        iteration, soil_rung, support, rel(final_soil_path),
        rel(previous_soil_path) if previous_soil_path else ""))

    # -- C
    applied_verdict_path = (Path(args.verdict_applied) if args.verdict_applied
                            else root / "hydrography" / "analysis" /
                            "carve_verdict.json")
    modelled_verdict_path = (Path(args.verdict_modelled) if args.verdict_modelled
                             else root / "hydrography" / "analysis" /
                             "carve_verdict_modelled.json")
    applied_verdict = read_json(applied_verdict_path)
    modelled_verdict = read_json(modelled_verdict_path)
    applied_arm = modelled_arm = None
    if applied_verdict and applied_verdict.get("climatology"):
        applied_arm = arm_of(
            resolve_under_root(applied_verdict["climatology"], root), root)
    if modelled_verdict and modelled_verdict.get("climatology"):
        modelled_arm = arm_of(
            resolve_under_root(modelled_verdict["climatology"], root), root)
    verdicts.append(loop_c(applied_verdict, modelled_verdict, applied_arm,
                           modelled_arm, support, rel(applied_verdict_path),
                           rel(modelled_verdict_path)))

    verdicts.sort(key=lambda v: v.loop)
    return verdicts


def show(verdicts: list[Verdict]) -> None:
    mark = {PASS: "  ok  ", FAIL: " FAIL ", NOT_EVALUABLE: " n/e  "}
    for v in verdicts:
        print(f"[{mark[v.status]}] loop {v.loop}  {v.detail}")
    n_fail = sum(1 for v in verdicts if v.status == FAIL)
    n_ne = sum(1 for v in verdicts if v.status == NOT_EVALUABLE)
    print(f"\n{len(verdicts)} loops, {n_fail} failed, {n_ne} not evaluable")
    if n_fail:
        print("A failure NAMES a loop. Re-entering one at the operating support "
              "is a commissioning-scale purchase and nothing here takes that "
              "decision.")


def write_report(verdicts: list[Verdict], path: Path, root: Path = ROOT) -> None:
    payload = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "what": "each pipeline loop's OWN exit predicate, re-evaluated against "
                "the final state. docs/src/pipeline/loops.md, the finalizer.",
        "not_a_loop": "nothing iterates on this verdict. A failure names a loop "
                      "and the decision to re-enter it is the author's.",
        "operating_support": rungs.ESCALATION_ROUTE[-1][0],
        "status": (FAIL if any(v.status == FAIL for v in verdicts) else
                   NOT_EVALUABLE if any(v.status == NOT_EVALUABLE for v in verdicts)
                   else PASS),
        "loops": [v.as_json() for v in verdicts],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {rel(path, root)}")


# ---------------------------------------------------------------------------
# checking the checker
# ---------------------------------------------------------------------------

def _carve_list(rows, terrain="t0", climatology="clim.nc",
                endmember="clim_cold.nc", single=False):
    """A carve list sidecar, in the shape `export_carve_list.py` writes.

    `rows` is (id, retain, warm_retain, cold_retain). The three retains are
    given SEPARATELY on purpose: that is what lets a fixture carry a carved set
    which is not the intersection of the arms it records, which is the case the
    A predicate exists to catch and which no correct writer produces.
    """
    return {
        "terrain_hash": terrain,
        "climatology": climatology,
        "intersection": ({"single_climate": True} if single else
                         {"single_climate": False,
                          "endmember_climatology": endmember}),
        "basins": [
            {"id": i, "retain": r,
             "verdict": "carve" if r <= 0.0 else ("preserve" if r >= 1.0
                                                  else "marginal"),
             "retain_warm_vegetated_arm": None if single else w,
             "retain_cold_bare_rock_arm": None if single else c}
            for i, r, w, c in rows],
    }


def _arm(rung, source, exists=True):
    return {"path": "fixture", "exists": exists, "rung": rung,
            "run_id": "run_fixture", "land_albedo_source": source}


def _soil(cells):
    """A soilmap as `read_soilmap` returns one. `cells` is (lon, lat, clay, soilc)."""
    return {(lon, lat): {"Lon": lon, "Lat": lat, "clay": clay, "soilc": soilc}
            for lon, lat, clay, soilc in cells}


def self_test() -> int:
    cases: list[tuple[str, bool, str]] = []
    failures: list[str] = []
    reached: dict[str, set[str]] = {}

    def case(name, got, want, why):
        ok = got == want
        cases.append((name, ok, f"{got} ({why})"))
        if not ok:
            failures.append(f"{name}: expected {want}, got {got}")
        # Which statuses each predicate has been driven to. Asserted at the end,
        # because a predicate that can only refuse would pass every negative
        # case here and no positive one.
        if want in (PASS, FAIL, NOT_EVALUABLE):
            reached.setdefault(name.split()[0], set()).add(want)

    support = "T85"
    warm = _arm(support, "vegetated")
    cold = _arm(support, COLD_ARM)

    # -- A ------------------------------------------------------------------
    #
    # The intersection is {b1}: b1 is cut by both arms, b2 by the warm arm only.
    agreeing = _carve_list([("b1", 0.0, 0.0, 0.0),
                            ("b2", 1.0, 0.0, 1.0),
                            ("b3", 1.0, 1.0, 1.0)])
    applied_ok = _carve_list([("b1", 0.0, 0.0, 0.0),
                              ("b2", 1.0, 0.0, 1.0),
                              ("b3", 1.0, 1.0, 1.0)])
    case("A passes when the carved set IS the re-taken intersection",
         loop_a(applied_ok, agreeing, support, warm, cold).status, PASS,
         "b1 in both arms and b1 carved")

    # THE CASE THE PREDICATE EXISTS FOR: a carved set that is not the
    # intersection of the two verdicts the same file records. b2 is carved and
    # only the warm arm cuts it.
    not_intersection = _carve_list([("b1", 0.0, 0.0, 0.0),
                                    ("b2", 0.0, 0.0, 1.0),
                                    ("b3", 1.0, 1.0, 1.0)])
    case("A fails a carved set that is not the intersection of its own arms",
         loop_a(applied_ok, not_intersection, support, warm, cold).status, FAIL,
         "b2 carved, cut by the warm arm alone")

    # And the other direction: the construction is sound and the APPLIED set is
    # not what it produces.
    applied_extra = _carve_list([("b1", 0.0, 0.0, 0.0),
                                 ("b2", 0.0, 0.0, 1.0),
                                 ("b3", 1.0, 1.0, 1.0)])
    case("A fails an applied set the re-taken intersection does not contain",
         loop_a(applied_extra, agreeing, support, warm, cold).status, FAIL,
         "b2 was carved and the re-take's intersection is b1 alone")

    case("A refuses a re-take at a coarser support",
         loop_a(applied_ok, agreeing, support, _arm("T21", "vegetated"),
                cold).status, NOT_EVALUABLE,
         "a T21 climatology is not a substitute for T85's")
    case("A refuses a cold arm that is not the bare-rock bound",
         loop_a(applied_ok, agreeing, support, warm,
                _arm(support, "vegetated")).status, NOT_EVALUABLE,
         "both arms vegetated is one arm run twice")
    case("A refuses a re-take on one climate",
         loop_a(applied_ok, _carve_list([("b1", 0.0, 0, 0)], single=True),
                support, warm, cold).status, NOT_EVALUABLE,
         "loop A's exit is the intersection of two")
    case("A refuses a re-take that does not cover every applied carve",
         loop_a(applied_extra, _carve_list([("b1", 0.0, 0.0, 0.0)]), support,
                warm, cold).status, NOT_EVALUABLE,
         "b2 was carved and the re-take has no verdict on it")
    case("A refuses a re-take on another terrain",
         loop_a(applied_ok, _carve_list([("b1", 0.0, 0.0, 0.0)], terrain="t1"),
                support, warm, cold).status, NOT_EVALUABLE,
         "the re-take holds geometry still and moves only the climate")
    case("A refuses when nothing was carved to compare against",
         loop_a(None, agreeing, support, warm, cold).status, NOT_EVALUABLE,
         "no applied carve list")

    # -- B ------------------------------------------------------------------
    criteria = {"soil_carbon_relative_tolerance": 0.02,
                "texture_cells_moved_tolerance": 0.05,
                "maximum_iterations": 6}
    settled_prev = _soil([(0.0, 0.0, 0.30, 5.00), (1.0, 0.0, 0.40, 6.00),
                          (2.0, 0.0, 0.20, 4.00), (3.0, 0.0, 0.25, 5.00)])
    settled_now = _soil([(0.0, 0.0, 0.305, 5.02), (1.0, 0.0, 0.405, 6.02),
                         (2.0, 0.0, 0.205, 4.02), (3.0, 0.0, 0.255, 5.02)])
    case("B passes soil inside both tolerances",
         loop_b(criteria, settled_now, settled_prev, 3, support, support).status,
         PASS, "carbon moved 0.004 and no cell moved clay past 0.02")

    carbon_moved = _soil([(0.0, 0.0, 0.305, 6.00), (1.0, 0.0, 0.405, 7.00),
                          (2.0, 0.0, 0.205, 5.00), (3.0, 0.0, 0.255, 6.00)])
    case("B fails soil carbon outside its tolerance",
         loop_b(criteria, carbon_moved, settled_prev, 3, support,
                support).status, FAIL,
         "land-mean carbon moved about 0.19 against 0.02")

    texture_moved = _soil([(0.0, 0.0, 0.40, 5.00), (1.0, 0.0, 0.50, 6.00),
                           (2.0, 0.0, 0.205, 4.02), (3.0, 0.0, 0.255, 5.02)])
    case("B fails texture outside its tolerance",
         loop_b(criteria, texture_moved, settled_prev, 3, support,
                support).status, FAIL,
         "half the cells moved clay past 0.02 against a tolerance of 0.05")

    case("B fails a soil past the iteration cap",
         loop_b(criteria, settled_now, settled_prev, 7, support,
                support).status, FAIL, "iteration 7 against a cap of 6")
    case("B refuses a soil built at a coarser support",
         loop_b(criteria, settled_now, settled_prev, 3, "T21",
                support).status, NOT_EVALUABLE,
         "a soil converged on a coarser climate is not the final state's")
    case("B refuses with no previous iteration to difference",
         loop_b(criteria, settled_now, None, 3, support, support).status,
         NOT_EVALUABLE, "the criteria are differences between iterations")
    case("B refuses a previous iteration carrying no biosphere carbon",
         loop_b(criteria, settled_now,
                _soil([(0.0, 0.0, 0.30, 0.0), (1.0, 0.0, 0.40, 0.0),
                       (2.0, 0.0, 0.20, 0.0), (3.0, 0.0, 0.25, 0.0)]),
                1, support, support).status, NOT_EVALUABLE,
         "a relative change against the declared initial state is undefined")

    # -- C ------------------------------------------------------------------
    modelled = _arm(support, MODELLED_ARM)
    verdict_before = {"climatology": "clim.nc",
                      "carve_list_penman": ["b1", "b2", "b4"]}
    verdict_same = {"climatology": "clim_modelled.nc",
                    "carve_list_penman": ["b4", "b1", "b2"]}
    case("C passes when no basin flips",
         loop_c(verdict_before, verdict_same, warm, modelled, support).status,
         PASS, "the same three basins carve, order aside")

    verdict_flipped = {"climatology": "clim_modelled.nc",
                       "carve_list_penman": ["b1", "b2", "b5"]}
    case("C fails when a basin flips",
         loop_c(verdict_before, verdict_flipped, warm, modelled,
                support).status, FAIL, "b4 drops out and b5 appears")

    case("C refuses a re-take at a coarser support",
         loop_c(verdict_before, verdict_same, warm, _arm("T21", MODELLED_ARM),
                support).status, NOT_EVALUABLE,
         "the exit generalised to the operating support needs that support")
    case("C refuses a re-take on an assumed canopy",
         loop_c(verdict_before, verdict_same, warm, _arm(support, "vegetated"),
                support).status, NOT_EVALUABLE,
         "an assumed canopy is the arm loop A already used")
    case("C refuses with no modelled verdict",
         loop_c(verdict_before, None, warm, None, support).status,
         NOT_EVALUABLE, "there is no re-take")

    # -- D ------------------------------------------------------------------
    def hop(run, rung, dt, converted):
        return {"run_id": run, "rung": rung, "timestep_minutes": dt,
                "converted": converted}

    walked = [hop("r0", "T21", 45.0, False), hop("r1", "T42", 45.0, True),
              hop("r2", support, 45.0, True)]
    case("D passes a chain that reaches the support at constant dt",
         loop_d(walked, support, None).status, PASS,
         "two conversions, neither moving the step")

    moved_step = [hop("r0", "T21", 45.0, False), hop("r1", "T42", 30.0, True),
                  hop("r2", support, 30.0, True)]
    case("D fails a conversion across a change of step",
         loop_d(moved_step, support, None).status, FAIL,
         "T21 at 45 converted to T42 at 30")

    short = [hop("r0", "T21", 45.0, False), hop("r1", "T42", 45.0, True)]
    case("D fails a chain that stops below the operating support",
         loop_d(short, support, None).status, FAIL,
         f"the chain ends at T42 and the operating support is {support}")

    unconverted = [hop("r0", "T21", 45.0, False),
                   hop("r1", support, 45.0, False)]
    case("D fails a rung change nothing converted",
         loop_d(unconverted, support, None).status, FAIL,
         "the manifest records no conversion across the rung change")

    case("D fails when lib/rungs.py's own route invariants do not hold",
         loop_d(walked, support, "a conversion at a changed step").status, FAIL,
         "the declared route is what the runs are judged against")
    case("D refuses a chain with no rung change on it",
         loop_d([hop("r0", support, 45.0, False),
                 hop("r1", support, 45.0, False)], support, None).status,
         NOT_EVALUABLE,
         "a route walked in one rung is not evidence about conversions")
    case("D refuses with no chain at all",
         loop_d([], support, None).status, NOT_EVALUABLE,
         "no baseline climatology names a run")

    # -- every predicate must reach all three statuses -----------------------
    #
    # A predicate that can only refuse passes every negative case above and no
    # positive one, which is the shape this whole file exists to prevent. Taken
    # as a snapshot first, because each of these cases adds to `reached` itself.
    coverage = {loop: sorted(reached.get(loop, set())) for loop in "ABCD"}
    want_all = sorted({PASS, FAIL, NOT_EVALUABLE})
    for loop in "ABCD":
        cases.append((f"{loop} was driven to pass, to fail and to refuse",
                      coverage[loop] == want_all, str(coverage[loop])))
        if coverage[loop] != want_all:
            failures.append(f"loop {loop} only reached {coverage[loop]}")

    width = max(len(n) for n, _, _ in cases) + 2
    for name, ok, why in cases:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name:<{width}} {why}")
    print(f"\n{len(cases)} cases, {len(failures)} failed")
    for f in failures:
        print(f"  {f}")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Re-evaluate each pipeline loop's own exit predicate "
                    "against the final state. A verification, not a loop: "
                    "nothing iterates on its verdict. Exit 0 when all four "
                    "hold, 1 when any fails, 2 when none failed and at least "
                    "one could not be evaluated.")
    parser.add_argument(
        "--self-test", action="store_true",
        help="check the checker instead of the tree: drive each predicate to "
             "pass, to fail and to refuse on fixtures. Touches nothing on disk.")
    parser.add_argument(
        "--applied", type=Path, default=None,
        help="loop A: the carve list that was actually carved. Defaults to "
             "hydrography/data/<source_build>/carve_list.json")
    parser.add_argument(
        "--retaken", type=Path, default=None,
        help="loop A: a carve list re-taken at the operating support with both "
             "bounding climates. Defaults to "
             "hydrography/analysis/carve_list_final.json")
    parser.add_argument(
        "--soil", type=Path, default=None,
        help="loop B: the soil the final state carries. Defaults to "
             "pedology/data/<source_build>/soilmap_<operating support>.txt")
    parser.add_argument(
        "--previous-soil", type=Path, default=None,
        help="loop B: the previous iteration's soil. There is no default: the "
             "soil step overwrites one file per build and rung, so the tree "
             "keeps only the latest")
    parser.add_argument(
        "--verdict-applied", type=Path, default=None,
        help="loop C: the carve verdict loop A exited on. Defaults to "
             "hydrography/analysis/carve_verdict.json")
    parser.add_argument(
        "--verdict-modelled", type=Path, default=None,
        help="loop C: the carve verdict re-taken on MODELLED vegetation at the "
             "operating support. Defaults to "
             "hydrography/analysis/carve_verdict_modelled.json")
    parser.add_argument(
        "--output", type=Path, default=REPORT,
        help=f"where the report goes. Default {rel(REPORT)}")
    args = parser.parse_args()
    if args.self_test:
        return self_test()

    verdicts = gather(args)
    show(verdicts)
    write_report(verdicts, args.output)
    if any(v.status == FAIL for v in verdicts):
        return 1
    if any(v.status == NOT_EVALUABLE for v in verdicts):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
