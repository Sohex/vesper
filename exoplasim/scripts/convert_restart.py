#!/usr/bin/env python3
"""Convert a model restart across horizontal resolution and real precision.

Worldbuilding frame: the restart is the saved state of the Vesper climate
model. Nothing here is about the real world.

    python exoplasim/scripts/convert_restart.py SOURCE OUTPUT \
        --target-template TARGET_TEMPLATE --report OUTPUT.conversion.json

WHAT THIS IS FOR. Advancing a spun-up state up the resolution ladder --
T21, T42, T85, T127, T170 -- instead of paying for a cold start at every rung.

WHAT THE RESULT IS. A NEW INITIAL CONDITION, and nothing more. It is not a
bitwise continuation of the donor experiment and it is not evidence that the
source and the target are the same experiment: resolution conversion changes
the represented state on purpose, and the model is chaotic. The converted state
begins a new run lineage and its equilibrium belongs to the target model.

THE TARGET TEMPLATE IS THE CENTRE OF THE DESIGN. It is a clean restart written
by the exact target executable, and it supplies the record set and order, the
target precision, the target-resolution static surface fields, the compiler's
seed shape, and the accumulators' own reset values. The converter starts from
the template and overlays only state whose policy permits transfer. It never
uses an evolved target run as a source of fallback state -- that is how a
donor's stale surface field supersedes a newly staged one, which is the whole
subject of `restart_surface.py`.

THE CONTRACT IS DELIBERATELY NARROW. Equal vertical dimensions, equal timestep,
equal enabled physics, and a grid where NLON is twice NLAT. Timestep and
vertical-resolution conversion are separate features with their own reasons:
the file stores `nstep` rather than elapsed time and both leapfrog levels, so
changing the step duration changes the represented derivative.

WHAT IS NOT HERE, and is `--report`ed rather than silently done: the model-owned
post-load fixup. Fields that are functions of other restart state -- albedo,
roughness, saturation humidity -- cannot be rebuilt here without writing a
second implementation of the model's physics, so they arrive holding the
template's values and the report names every one of them.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import _paths
import restart_format as rf
import restart_schema as rs
import restart_transforms as rt
from restart_schema import ConversionError, _int, infer_real_bytes

# The version of the conversion contract itself. A change to what a policy
# DOES belongs in this number, because a report is the only record of how a
# state was made.
CONVERTER_VERSION = 1

# Configuration reals that must agree between source and target. Declared here,
# before any conversion has been run, because a tolerance chosen after seeing
# the answer is not a tolerance. One part in a million is far tighter than the
# 8-to-4-byte cast these values may have crossed and far looser than any real
# configuration change, which moves them by percent.
CONFIG_RTOL = 1e-6

REAL = {4: np.dtype("<f4"), 8: np.dtype("<f8")}


# ---------------------------------------------------------------------------
# Reading a file into geometry
# ---------------------------------------------------------------------------

@dataclass
class RestartState:
    """One parsed restart: its records, its geometry and its precision."""

    path: Path
    records: list
    by_name: dict
    geometry: rs.Geometry
    real_bytes: int

    @property
    def names(self) -> list:
        return [r.name for r in self.records]

    def decode(self, name: str) -> np.ndarray:
        """A real record as float64, whatever width it was stored at."""
        rec = self.by_name[name]
        return np.frombuffer(rec.payload, dtype=REAL[self.real_bytes]).astype(np.float64)


def load(path: Path) -> RestartState:
    records = rf.read(path)
    try:
        geometry, real_bytes = rs.describe(records)
    except ConversionError as exc:
        raise ConversionError(f"{path}: {exc}") from None
    return RestartState(path=Path(path), records=records,
                        by_name=rf.index(records), geometry=geometry,
                        real_bytes=real_bytes)


# ---------------------------------------------------------------------------
# Refusals: everything checked before a byte of output exists
# ---------------------------------------------------------------------------

def check_compatible(src: RestartState, tgt: RestartState,
                     inventory: dict) -> None:
    """Refuse every conversion this contract does not cover, before writing."""
    problems = []

    for state, where in ((src, "source"), (tgt, "template")):
        for rec in state.records:
            if rec.name not in rs.POLICY or rec.name not in inventory:
                continue
            if rs.POLICY[rec.name].action == rs.SEED:
                want = [state.geometry.nseedlen * 4]
            elif inventory[rec.name].writer == "put_restart_integer":
                want = [4]
            else:
                want = [n * state.real_bytes for n in
                        rs.candidate_counts(rec.name, state.geometry, inventory)]
            if rec.nbytes not in want:
                problems.append(
                    f"the {where}'s '{rec.name}' is {rec.nbytes} bytes and a "
                    f"{state.geometry.label} grid at {state.real_bytes}-byte "
                    f"reals gives {want}. The file and its own headers "
                    "disagree, so nothing downstream can trust either.")

    unknown_src = sorted(set(src.names) - set(rs.POLICY))
    unknown_tgt = sorted(set(tgt.names) - set(rs.POLICY))
    for names, where in ((unknown_src, src.path), (unknown_tgt, tgt.path)):
        if names:
            problems.append(
                f"{where} holds records with no schema entry: "
                f"{', '.join(names)}. An enabled name absent from the schema "
                "is a hard error, not a record to pass through.")

    missing = sorted(set(src.names) - set(tgt.names))
    extra = sorted(set(tgt.names) - set(src.names))
    if missing or extra:
        detail = []
        if missing:
            detail.append(f"only in the source: {', '.join(missing)}")
        if extra:
            detail.append(f"only in the template: {', '.join(extra)}")
        problems.append(
            "the two record sets differ, so the two runs do not have the same "
            "physics enabled (" + "; ".join(detail) + ")")

    for dim in ("nlev", "nlsoil", "nlev_oce"):
        a = getattr(src.geometry, dim), getattr(tgt.geometry, dim)
        if a[0] != a[1]:
            problems.append(
                f"{dim} is {a[0]} in the source and {a[1]} in the template. "
                "Vertical conversion needs pressure-coordinate interpolation "
                "and its own conservation argument, and is a separate feature.")

    vegetation = sorted(n for n in src.names
                        if rs.POLICY.get(n) and rs.POLICY[n].vegetation_owned
                        and n in ("dforest", "dwmax")
                        and any(v in src.by_name
                                for v in ("dcveg", "dcsoil")))
    if vegetation:
        problems.append(
            "the source carries SIMBA's carbon pools, so coupled vegetation "
            f"owns {', '.join(vegetation)} and they are prognostic rather "
            "than the staged boundary fields this converter treats them as. "
            "That case needs its own policy.")

    for name, pol in rs.POLICY.items():
        if pol.action != rs.REQUIRE_EQUAL or name not in src.by_name:
            continue
        if src.by_name[name].nbytes == 4 and tgt.by_name[name].nbytes == 4:
            a, b = _int(src.by_name[name]), _int(tgt.by_name[name])
            same = a == b
        else:
            a, b = src.decode(name), tgt.decode(name)
            same = a.shape == b.shape and np.allclose(a, b, rtol=CONFIG_RTOL,
                                                      atol=0.0)
        if not same:
            problems.append(
                f"'{name}' is {a} in the source and {b} in the template; "
                f"{pol.why}, so this is a configuration change rather than a "
                "change of resolution")

    # A template's accumulators must already be at the value the model resets
    # them to. The converter takes them verbatim, so a template cut from a run
    # mid-window hands the converted run somebody else's partial accumulation
    # -- which is CLIM-31 again, by another route. `build_restart_template.py`
    # is what makes one clean.
    dirty = []
    for rec in tgt.records:
        pol = rs.POLICY.get(rec.name)
        if pol is None or pol.semantic != rs.ACCUMULATOR:
            continue
        if pol.model_reset == "zero":
            want = b"\x00" * rec.nbytes
        elif pol.model_reset == "sentinel":
            count = rec.nbytes // tgt.real_bytes
            want = np.full(count, pol.reset_value,
                           dtype=REAL[tgt.real_bytes]).tobytes()
        else:
            continue
        if rec.payload != want:
            dirty.append(rec.name)
    if dirty:
        shown = ", ".join(dirty[:8]) + (" ..." if len(dirty) > 8 else "")
        problems.append(
            f"{len(dirty)} of the template's accumulators are not at the "
            f"value the model resets them to ({shown}). It was cut from a run "
            "mid-window, so a conversion onto it would open the new run on "
            "somebody else's partial accumulation. Run "
            "build_restart_template.py on it first.")

    if problems:
        raise ConversionError("\n  - ".join(["this conversion is refused:"] + problems))


# ---------------------------------------------------------------------------
# The conversion
# ---------------------------------------------------------------------------

@dataclass
class RecordReport:
    name: str
    semantic: str
    action: str
    source_bytes: int
    target_bytes: int
    detail: dict = field(default_factory=dict)


def _domain_mask(state: RestartState, domain: str):
    """Source cells of one surface class, from the model's own land-sea mask.

    `dls` is the mask every other one is built from, so it is the single
    source rather than `xls` or `yls`: three masks that ought to agree are
    three chances for them not to.
    """
    if domain == "global":
        return None
    land = state.decode("dls") > 0.5
    return land if domain == "land" else ~land


def convert(src: RestartState, tgt: RestartState, *,
            seed_override: bytes | None = None):
    """Every output record, in the template's order, with a report for each."""
    same_grid = src.geometry.nlat == tgt.geometry.nlat
    weights = None if same_grid else rt.build_weights(src.geometry.nlat,
                                                      tgt.geometry.nlat)
    src_area = rt.cell_area(src.geometry.nlat)
    tgt_area = rt.cell_area(tgt.geometry.nlat)
    masks = {d: _domain_mask(src, d) for d in ("global", "land", "ocean")}
    out_dtype = REAL[tgt.real_bytes]
    fractions: dict[str, np.ndarray] = {}      # remapped covers, for their pairs
    records, reports = [], []

    def levels(name: str, state: RestartState, per_level: int):
        flat = state.decode(name)
        if flat.size % per_level:
            raise ConversionError(
                f"'{name}' holds {flat.size} reals, not a whole number of "
                f"{per_level}-element levels")
        return flat.reshape(-1, per_level)

    def emit(name, payload, semantic, action, detail=None):
        records.append(rf.Record(name=name, payload=payload, offset=-1))
        reports.append(RecordReport(
            name=name, semantic=semantic, action=action,
            source_bytes=src.by_name[name].nbytes if name in src.by_name else 0,
            target_bytes=len(payload), detail=detail or {}))

    def cast(values: np.ndarray) -> tuple:
        """Target-width bytes, plus the error the narrowing introduced."""
        narrowed = values.astype(out_dtype)
        widened = narrowed.astype(np.float64)
        if not np.all(np.isfinite(widened)):
            raise ConversionError(
                "the cast to four-byte reals produced a value that is not "
                "finite; the source holds magnitudes this precision cannot "
                "represent")
        err = np.abs(widened - values)
        worst = float(err.max()) if err.size else 0.0
        scale = np.abs(values)
        rel = float(np.max(np.where(scale > 0, err / np.where(scale > 0, scale, 1.0), 0.0))) \
            if err.size else 0.0
        detail = {} if worst == 0.0 else {"cast_abs_error": worst,
                                          "cast_rel_error": rel}
        return narrowed.tobytes(), detail

    for rec in tgt.records:
        name = rec.name
        pol = rs.POLICY[name]

        if pol.action in (rs.TARGET, rs.RESET, rs.RECOMPUTE, rs.REQUIRE_EQUAL):
            emit(name, rec.payload, pol.semantic, pol.action)

        elif pol.action == rs.COPY:
            emit(name, src.by_name[name].payload, pol.semantic, pol.action,
                 {"value": _int(src.by_name[name])})

        elif pol.action == rs.SEED:
            if seed_override is not None:
                emit(name, seed_override, pol.semantic, pol.action,
                     {"source": "explicit"})
            elif src.by_name[name].nbytes == rec.nbytes:
                emit(name, src.by_name[name].payload, pol.semantic, pol.action,
                     {"source": "donor"})
            else:
                raise ConversionError(
                    f"the donor's seed is {src.by_name[name].nbytes} bytes and "
                    f"the target executable's is {rec.nbytes}. The random "
                    "state belongs to the compiler, so there is no right "
                    "answer here to infer: pass --seed to name one, or "
                    "--keep-template-seed to accept the template's.")

        elif pol.action == rs.PROJECT:
            src_lv = levels(name, src, src.geometry.nrsp)
            out = np.zeros((src_lv.shape[0], tgt.geometry.nrsp))
            l2, mx = 0.0, 0.0
            for k, level in enumerate(src_lv):
                out[k], (a, b) = rt.project_spectral(
                    level, src.geometry.ntru, tgt.geometry.ntru)
                l2, mx = float(np.hypot(l2, a)), max(mx, b)
            payload, detail = cast(out.ravel())
            detail.update({"discarded_l2": l2, "discarded_max": mx,
                           "modes_source": src.geometry.nrsp // 2,
                           "modes_target": tgt.geometry.nrsp // 2})
            emit(name, payload, pol.semantic, pol.action, detail)

        elif pol.action == rs.REMAP:
            if same_grid:
                payload, detail = cast(src.decode(name))
                emit(name, payload, pol.semantic, "copy_same_grid", detail)
                continue
            src_lv = levels(name, src, src.geometry.nugp)
            fallback = np.frombuffer(rec.payload, dtype=REAL[tgt.real_bytes]
                                     ).astype(np.float64).reshape(-1, tgt.geometry.nugp)
            out = np.empty((src_lv.shape[0], tgt.geometry.nugp))
            detail = {"remap": pol.remap, "domain": pol.domain}
            unfilled = 0
            for k, level in enumerate(src_lv):
                if pol.remap == rs.THICKNESS:
                    cover = fractions.get(pol.partner)
                    if cover is None:
                        raise ConversionError(
                            f"'{name}' is remapped as a volume with "
                            f"'{pol.partner}', which has not been converted "
                            "yet; the template's record order has moved")
                    volume, missed = rt.remap(level * src.decode(pol.partner),
                                              weights, mask=masks[pol.domain])
                    with np.errstate(invalid="ignore", divide="ignore"):
                        value = np.where(cover > 0, volume / np.where(cover > 0, cover, 1.0), 0.0)
                    detail["conserved_as"] = "volume with " + pol.partner
                else:
                    value, missed = rt.remap(level, weights, mask=masks[pol.domain])
                gaps = missed | ~np.isfinite(value)
                unfilled += int(gaps.sum())
                value = np.where(gaps, fallback[k % fallback.shape[0]], value)
                if pol.bounds is not None:
                    lo, hi = pol.bounds
                    clipped = np.clip(value, lo, hi)
                    detail["clipped"] = detail.get("clipped", 0) + \
                        int(np.count_nonzero(clipped != value))
                    value = clipped
                if pol.conserve or pol.remap == rs.RESERVOIR:
                    before = float((level * src_area).sum())
                    after = float((value * tgt_area).sum())
                    detail.setdefault("inventory", []).append(
                        {"level": k, "source": before, "target": after,
                         "relative_residual": (abs(after - before) / abs(before))
                         if before else 0.0})
                out[k] = value
            if pol.remap == rs.FRACTION:
                fractions[name] = out[0]
            detail["fallback_cells"] = unfilled
            if unfilled:
                detail["fallback"] = ("the target template, for target cells "
                                      "with no source overlap of their class")
            if pol.conserve:
                detail["conserves"] = pol.conserve
            payload, cast_detail = cast(out.ravel())
            detail.update(cast_detail)
            emit(name, payload, pol.semantic, pol.action, detail)

        else:
            raise ConversionError(
                f"'{name}' has action '{pol.action}', which this converter "
                "does not implement")

    return records, reports


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _template_identity(template: Path) -> dict:
    """What a consumer has to check the converted state against."""
    prov = template_provenance(template)
    if prov is None:
        return {"provenance": None,
                "warning": ("no provenance sidecar beside this template, so "
                            "nothing here says which staged surface its static "
                            "records are or which executable wrote it")}
    return {"provenance": str(Path(str(template) + ".provenance.json")),
            "executable": prov.get("executable"),
            "source_build": prov.get("source_build"),
            "config_sha256": prov.get("config_sha256"),
            "surface_field_sha256": prov.get("surface_field_sha256"),
            "cut_from_run": (prov.get("cut_from") or {}).get("run_id")}


def _geometry_dict(state: RestartState) -> dict:
    g = state.geometry
    return {"truncation": g.label, "nlat": g.nlat, "nlon": g.nlon,
            "ntru": g.ntru, "nrsp": g.nrsp, "nesp": g.nesp, "nugp": g.nugp,
            "nlev": g.nlev, "nlsoil": g.nlsoil, "nlev_oce": g.nlev_oce,
            "nseedlen": g.nseedlen, "real_bytes": state.real_bytes}


def template_provenance(template: Path) -> dict | None:
    """The sidecar `build_restart_template.py` wrote beside a template.

    It carries which executable wrote the template and, crucially, WHICH STAGED
    SURFACE its static records are. A conversion is only sound onto a run whose
    staged surface matches those hashes, and the conversion report is the only
    thing a consumer has to check that against.
    """
    side = Path(str(template) + ".provenance.json")
    if not side.is_file():
        return None
    return json.loads(side.read_text(encoding="utf-8"))


def build_report(src, tgt, out_path, reports, source_manifest) -> dict:
    recompute = [r.name for r in reports if r.action == rs.RECOMPUTE]
    # An accumulator the model never resets spans the whole run rather than one
    # output window, so taking the template's value restarts a sum that was
    # never meant to restart. That is what a new lineage means, and it is named
    # here rather than left to be discovered in an output.
    never_reset = sorted(
        r.name for r in reports
        if rs.POLICY[r.name].semantic == rs.ACCUMULATOR
        and rs.POLICY[r.name].model_reset == "none")
    return {
        "schema_version": 1,
        "converter_version": CONVERTER_VERSION,
        "converter_sha256": sha256(Path(__file__)),
        "schema_sha256": sha256(Path(rs.__file__)),
        "transforms_sha256": sha256(Path(rt.__file__)),
        "source": {"path": str(src.path), "sha256": sha256(src.path),
                   **_geometry_dict(src)},
        "target_template": {"path": str(tgt.path), "sha256": sha256(tgt.path),
                            **_geometry_dict(tgt),
                            **_template_identity(tgt.path)},
        "output": {"path": str(out_path), "sha256": sha256(out_path)},
        "source_manifest": source_manifest,
        "config_rtol": CONFIG_RTOL,
        "records": [{"name": r.name, "semantic": r.semantic, "action": r.action,
                     "source_bytes": r.source_bytes,
                     "target_bytes": r.target_bytes, **r.detail}
                    for r in reports],
        "expected_to_change_in_model_fixup": recompute,
        "whole_run_accumulators_restarted": never_reset,
        "status": "initial_condition",
        "note": ("A new initial condition at the target support. Not a bitwise "
                 "continuation and not evidence that the source and target are "
                 "the same experiment: the records listed under "
                 "expected_to_change_in_model_fixup hold the template's values "
                 "and must be rebuilt by the target model before the state is "
                 "self-consistent."),
    }


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def inspect(path: Path) -> int:
    state = load(path)
    inventory = rs.inventory_from_source(_paths.MODEL_SRC / "plasim" / "src")
    g = state.geometry
    print(f"{path}")
    print(f"  {g.label}  nlat {g.nlat}  nlon {g.nlon}  nlev {g.nlev}  "
          f"nrsp {g.nrsp}  nesp {g.nesp}  reals {state.real_bytes} bytes")
    print(f"  {len(state.records)} records, of {len(inventory)} this build can write")
    by_class: dict = {}
    for rec in state.records:
        pol = rs.POLICY.get(rec.name)
        by_class.setdefault(pol.semantic if pol else "UNKNOWN", []).append(rec.name)
    for semantic in sorted(by_class):
        names = by_class[semantic]
        print(f"  {semantic:22s} {len(names):3d}  {' '.join(sorted(names)[:6])}"
              + (" ..." if len(names) > 6 else ""))
    unknown = by_class.get("UNKNOWN", [])
    return 1 if unknown else 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", nargs="?", type=Path,
                    help="the restart to convert from")
    ap.add_argument("output", nargs="?", type=Path,
                    help="where the converted restart is written")
    ap.add_argument("--target-template", type=Path,
                    help="a clean restart written by the EXACT target "
                         "executable and configuration. It supplies the record "
                         "set and order, the precision, the target-resolution "
                         "static fields, the seed shape and the accumulators' "
                         "reset values.")
    ap.add_argument("--source-manifest", type=Path,
                    help="the donor run's run_manifest.json, recorded in the "
                         "report as the converted state's provenance")
    ap.add_argument("--report", type=Path,
                    help="where the conversion report JSON is written")
    ap.add_argument("--seed", type=Path,
                    help="a file of four-byte integers to use as the random "
                         "seed, where the donor's shape does not fit the "
                         "target executable")
    ap.add_argument("--keep-template-seed", action="store_true",
                    help="accept the target template's own seed instead")
    ap.add_argument("--inspect", action="store_true",
                    help="report what a restart holds and stop")
    ap.add_argument("--check-template", action="store_true",
                    help="run every refusal against the source and template "
                         "and stop, writing nothing")
    ap.add_argument("--dry-run", action="store_true",
                    help="convert and report, but write no restart")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing output")
    ap.add_argument("--self-test", action="store_true",
                    help="prove the parser, the schema and both transforms "
                         "against fixtures, each with a negative control")
    args = ap.parse_args()

    if args.self_test:
        import restart_convert_selftest as st
        return st.run()
    if args.source is None:
        ap.error("a source restart is required")
    if args.inspect:
        return inspect(args.source)
    if args.target_template is None:
        ap.error("--target-template is required. A conversion without the "
                 "exact target's own record set, precision and static fields "
                 "is a guess at what the target expects.")

    src, tgt = load(args.source), load(args.target_template)
    inventory = rs.inventory_from_source(_paths.MODEL_SRC / "plasim" / "src")
    gaps = rs.check_policy_covers_source(_paths.MODEL_SRC / "plasim" / "src")
    if gaps:
        raise ConversionError(
            "the schema no longer covers the model source:\n  - "
            + "\n  - ".join(gaps))
    check_compatible(src, tgt, inventory)
    if args.check_template:
        print(f"template accepted: {src.geometry.label} {src.real_bytes}-byte "
              f"-> {tgt.geometry.label} {tgt.real_bytes}-byte, "
              f"{len(tgt.records)} records")
        return 0

    seed_override = None
    if args.seed is not None:
        seed_override = args.seed.read_bytes()
    elif args.keep_template_seed:
        seed_override = tgt.by_name["seed"].payload

    records, reports = convert(src, tgt, seed_override=seed_override)
    if args.output is None:
        ap.error("an output path is required")
    if args.dry_run:
        print(f"dry run: {len(records)} records, nothing written")
    else:
        rf.write(args.output, records, overwrite=args.force)
        print(f"wrote {args.output} ({src.geometry.label} -> "
              f"{tgt.geometry.label}, {len(records)} records)")
    if args.report and not args.dry_run:
        manifest = json.loads(args.source_manifest.read_text()) \
            if args.source_manifest else None
        report = build_report(src, tgt, args.output, reports,
                              {"path": str(args.source_manifest),
                               "run_id": manifest.get("run_id")}
                              if manifest else None)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(f"wrote {args.report}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ConversionError, rf.RestartFormatError) as exc:
        print(f"convert_restart: {exc}", file=sys.stderr)
        sys.exit(2)
