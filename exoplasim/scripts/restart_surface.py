#!/usr/bin/env python3
"""Check that a restart still carries the per-cell surface fields it was given.

    python exoplasim/scripts/restart_surface.py exoplasim/runs/run_xxxxxxxxxxxx
    python exoplasim/scripts/restart_surface.py RUN_DIR --restart MOST_REST.00007
    python exoplasim/scripts/restart_surface.py --all --json out.json
    python exoplasim/scripts/restart_surface.py --self-test

Worldbuilding. Vesper is an invented planet and this module is about the
simulation of it: a checkpoint file from a toy climate model and the surface
boundary fields staged for it. Nothing here is a claim about the world.

## What this is a test of

This project stages per-cell surface fields it generates itself -- soil water
capacity from pedology, background albedo from the lithology, roughness from the
aeolian chain -- as `.sra` files, one per ExoPlaSim surface code. The model reads
each of them once, on a cold start, and from then on carries its own copy in the
restart. Every resume takes the field from the restart, not from the file.

So there are two ways a staged field can stop being the field the model uses, and
only one of them was ever checked. `surface_field_report` checks that the `.sra`
files are PRESENT in the run directory, which catches the cold-start half.
Nothing checked the resume half: a restart branch that re-derives a field the
staged inputs already own leaves every `.sra` in place, every manifest
well-formed, and every downstream product the wrong answer.

That is not hypothetical. `landmod.f90` used to read the per-cell `dwmax` out of
the restart and then assign the scalar namelist `wsmax` over it on every resume,
so a soil water capacity that varies by two orders of magnitude across the
continents was flattened to one number after the first orbit. TASKS.md CLIM-67
removed the block; this is the check that would have named it.

IT IS A TEST AND NOT A DIAGNOSTIC, which is the reason it belongs on the resume
path rather than in a report. The staged `.sra` is a quantity the other side
already knows: the restart's copy of a boundary field has a right answer, stated
below per code, and any departure from it is wrong rather than merely different.

SCOPED BY THE PROPERTY, NOT BY THE FIELD. `dwmax` is the instance that fired.
Every surface code this project stages itself is checked here, on the same
argument, whether or not it has ever gone wrong -- `dwcl`'s collapse from 14
months to a scalar sat latent in the same block, and a check written for `dwmax`
alone would have missed it. The code list comes from
`run_exoplasim.intended_surface_codes`, so a code added there is checked here
without editing this file.

## Which field is the right answer, and why there are two

A resume has to hold two separate things true, and they have different
references. Both are stated here because collapsing them into one comparison
either refuses a decision the project already took or misses the defect this
module exists for.

1. THE RESTART HAS NOT CHANGED SINCE THE RUN ADOPTED ITS SURFACE. For a
   cold-started run the reference is the staged `.sra`; for a run seeded with
   `--restart-from` it is `MOST_REST.seed`, the accumulator-zeroed copy of the
   donor restart that `run_exoplasim.py` keeps in the run directory. A departure
   means the model re-derived a boundary field between two segments, which is
   CLIM-67's class exactly. This one is never overridable: nothing legitimately
   rewrites a boundary field on a resume.
2. THE STAGED `.sra` IN THE RUN DIRECTORY IS THE FIELD THE MODEL IS USING. True
   by construction after a cold start. Deliberately false for a run prepared with
   `--superseded-surface-ok`, which adopts a donor restart and DISCARDS the
   staged surface -- and `stage_surface_extras` then rewrites the current `.sra`
   into that run directory on every resume, so the directory and the manifest
   both assert a surface the model will not read. That is refused by default and
   restated per resume with `--superseded-surface-ok` on the resume, which stamps
   the segment rather than going quiet.

## The right answer, per code, fixed before any restart was read

The model's only operations between the staged `.sra` and the restart record are
a list-directed ASCII read into `real*8` (`-fdefault-real-8`, `configure.sh:73`),
an MPI scatter, and the matching gather at write time (`mpimod.f90:719`). None of
those is arithmetic. This module parses the same ASCII text with the same
correctly-rounded decimal-to-binary conversion, so the two sides are bit-equal
when nothing has substituted the field.

THE TOLERANCE IS THEREFORE EXACT EQUALITY, and it is exact because a tolerance
would be answering a question nobody asked: no round-off is available on this
path, so the smallest difference this check can see is already a substitution.
A tolerance here would only widen the target for the class of defect the check
exists to catch.

Three declared transforms sit on the expected side, not on the comparison:

* Code 172 `dls`, the land/sea mask, is binarised by `landmod.f90:349-353`
  (`> 0.5` to 1.0, everything else to 0.0) after `surfmod.f90:474-475` clamps it
  into [0, 1]. The right answer for the restart record is the binarised staged
  mask, compared exactly.
* Codes 174, 175 and 176 are 14-month arrays in the restart and single-record
  files on disk. `get_surf_array` copies a single record to all 14 months
  (`surfmod.f90:152-155`); a 12-record file would be expanded cyclically
  (`surfmod.f90:146-151`). The right answer is that expansion, and every month
  must match exactly.
* Anything with no restart record at all is re-read from its `.sra` at every
  start and cannot be substituted from a restart. Verified as an absence, not
  skipped: if a record by that name appears, the field has become restart-borne
  and the class of this check changed under it.

## What is deliberately NOT compared, and why

A field the model is entitled to change is not a field the staged file is the
right answer for, and comparing one would be a diagnostic wearing a test's
clothes. Each exclusion is reported by name, never silently dropped.

* Code 129 `doro`. The restart's orography is not the staged field: it is scaled
  by `oroscale`, spectrally fitted and truncated to `noromax`, and with
  `nglacier` on it is `groundoro + glacieroro` (`glaciermod.f90:358`). There is no
  right answer available from the `.sra` alone, so this check does not claim one.
  The remaining gap is real and tracked: nothing verifies that a restart's
  topography still derives from the staged 129.
* Codes 229 `dwmax` and 212 `dforest` under `NVEG = 2`. Coupled vegetation
  overwrites both prognostically every timestep (`simba.f90:484,487`), so under
  that setting the staged field stops being the right answer the moment the run
  starts. Under `NVEG` 0 or 1 nothing writes them and they are compared.
* Everything, when `NAQUA` or `NDESERT` is set: those replace the mask and the
  orography wholesale (`surfmod.f90:400-418`) and the run is not on this world.

`NEWSURF = 2` is refused rather than excluded. It assigns the namelist scalars
`wsmax`, `dz0land` and `albland` over `dwmax`, `dz0clim` and `dalbcl`
(`landmod.f90:492-499`) -- the identical substitution CLIM-67 removed, differing
only in being namelist-driven. A run that stages per-cell fields and then asks
for them to be flattened is stating two incompatible intents, and the check says
so instead of picking one.

## The mapping is derived, not duplicated

`SURFACE_RESTART_FIELDS` names a restart record per surface code. The code-to-name
binding is `surfmod.f90`'s own `surfcode` table, so `check_code_bindings` reads
that table out of the vendored source and refuses if any entry here disagrees.
A renumbering upstream fails loudly rather than quietly comparing the wrong
record.

The restart record walker is `diff_restarts.read_records`; the `.sra` reader is
`sra.read_sra`. Neither is reimplemented here.

## Index mapping

Both sides are flat arrays of NLAT*NLON in the same latitude-major order:
`write_sra` ravels (nlat, nlon) in C order, and `mpputgp` gathers to the same
global NUGP ordering that `mpsurfgp` scattered from. The comparison is
positional and no coordinate is reconstructed on either side, which is the only
form rule 3 permits across this boundary.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import MODEL_SRC, PROJECT_ROOT, RUNS  # noqa: E402
from diff_restarts import read_records  # noqa: E402
from reset_restart_accumulators import read_records as walk_records  # noqa: E402
from sra import read_sra  # noqa: E402


@dataclass(frozen=True)
class SurfaceRestartField:
    """One staged surface code and the restart record it must survive as.

    `record` is the restart record name, `months` its second dimension there,
    and `expand` how a single-record `.sra` becomes that many months. `binarise`
    marks the land/sea mask, whose right answer is the thresholded file rather
    than the file.
    """

    code: int
    record: str
    months: int = 1
    binarise: bool = False


# Every code `intended_surface_codes` can return, with the restart record it
# lands in. A code absent from here is reported as unmapped rather than passed:
# adding a generator without saying what happens to its field on a resume is the
# gap this whole module is about.
SURFACE_RESTART_FIELDS = {
    f.code: f
    for f in (
        SurfaceRestartField(172, "dls", binarise=True),
        SurfaceRestartField(173, "dz0clim"),
        SurfaceRestartField(174, "dalbcl", months=14),
        SurfaceRestartField(175, "dalbcl1", months=14),
        SurfaceRestartField(176, "dalbcl2", months=14),
        SurfaceRestartField(212, "dforest"),
        SurfaceRestartField(229, "dwmax"),
    )
}

# Codes whose field never reaches a restart: `radini` and `aero_ini` read them
# from the `.sra` at every start and abort if the file is missing
# (`radmod.f90:1285`, `aeromod.f90:326-328`). The check on these is that the
# restart still carries no record by that name, which is what makes them
# unsubstitutable.
REREAD_EVERY_START = {
    1811: "ddustcol",
    1801: "dsrcw",
    1802: "ddrage",
    1803: "dwpr",
}

# Fields SIMBA owns once coupled vegetation is on: `simba.f90:484,487` assigns
# both every timestep under `nveg == 2`, so the staged file stops being their
# right answer the moment such a run starts. Named as a set rather than tested
# by literal, because anything else simba starts writing belongs here beside
# them.
VEGETATION_OWNED = {212, 229}

# Not comparable, with the reason reported rather than the code dropped. See the
# module docstring.
NOT_COMPARABLE = {
    129: ("the restart's doro is scaled by oroscale, spectrally fitted and "
          "truncated to noromax, and with nglacier on it carries the ice-sheet "
          "orography as well (glaciermod.f90:358), so the staged file is not "
          "the right answer for it"),
}


def surfmod_bindings(source: Path | None = None) -> dict[str, int]:
    """{array name: surface code} as `surfmod.f90` itself registers them."""
    path = source or (MODEL_SRC / "plasim" / "src" / "surfmod.f90")
    text = path.read_text(encoding="latin-1")
    out = {}
    for code, name in re.findall(
            r"call\s+surfcode\s*\(\s*(\d+)\s*,\s*'([^']+)'\s*\)", text):
        out[name.strip()] = int(code)
    if not out:
        raise RuntimeError(f"{path} registered no surface codes; the parse is wrong")
    return out


def check_code_bindings(source: Path | None = None) -> None:
    """Refuse if this module's code-to-record mapping disagrees with the model.

    The vendored source is the only authority for which code fills which array.
    Comparing the wrong record would be a check that cannot fail, so a drift here
    is fatal rather than a warning.
    """
    bindings = surfmod_bindings(source)
    wrong = []
    for field in SURFACE_RESTART_FIELDS.values():
        actual = bindings.get(field.record)
        if actual != field.code:
            wrong.append(f"{field.record} is code {actual} in surfmod.f90, "
                         f"not {field.code}")
    for code, name in REREAD_EVERY_START.items():
        actual = bindings.get(name)
        if actual != code:
            wrong.append(f"{name} is code {actual} in surfmod.f90, not {code}")
    if wrong:
        raise RuntimeError(
            "restart_surface.py disagrees with surfmod.f90 about which surface "
            "code fills which array, so it would compare the wrong restart "
            "record: " + "; ".join(wrong))


def namelist_int(run_dir: Path, filename: str, key: str, default: int) -> int:
    """One integer out of a staged namelist, or `default` if it is unset.

    Read from the run directory rather than from the config, because
    `configure()` rewrites the namelists on every continuation and the file is
    what the binary will actually parse.
    """
    path = run_dir / filename
    if not path.is_file():
        return default
    match = re.search(rf"^\s*{key}\s*=\s*(-?\d+)", path.read_text(encoding="latin-1"),
                      re.IGNORECASE | re.MULTILINE)
    return int(match.group(1)) if match else default


def expected_field(field: SurfaceRestartField, staged: np.ndarray) -> np.ndarray:
    """The restart contents a staged `.sra` entitles us to expect, as (months, n).

    The declared transforms live here and nowhere else: binarisation for the
    land/sea mask, and `get_surf_array`'s expansion of a single record to the
    14 months an annual-cycle array holds.
    """
    flat = np.asarray(staged, dtype=np.float64).ravel(order="C")
    if field.binarise:
        flat = np.where(flat > 0.5, 1.0, 0.0)
    return np.repeat(flat[None, :], field.months, axis=0)


def restart_field(raw: bytes, months: int, ncells: int) -> np.ndarray:
    """A restart record as (months, ncells).

    `put_restart_array` writes `pa(1:NUGP,1:klev)` from a Fortran array, so the
    payload is month-major: NUGP values for month 1, then month 2, and so on.
    """
    if len(raw) != 8 * months * ncells:
        raise RuntimeError(
            f"restart record is {len(raw)} bytes, expected "
            f"{8 * months * ncells} for {months} x {ncells} real*8. The restart "
            "layout is not what this check assumes; read it before trusting it.")
    return np.frombuffer(raw, dtype="<f8").reshape(months, ncells)


def verify_restart_surface_fields(run_dir: Path, restart: Path, codes: set[int],
                                  manifest: dict | None = None,
                                  allow_superseded: bool = False) -> dict:
    """Compare a restart's copy of every staged surface field against the file.

    Raises RuntimeError on any mismatch, on an unmapped code, on a `NEWSURF = 2`
    that would flatten the staged fields, or on a restart record appearing for a
    field that is supposed to be re-read every start. Returns the per-code
    verdicts for the run manifest.

    `codes` is `run_exoplasim.intended_surface_codes(config)`: the caller
    resolves it, so nothing here decides on its own which fields a run supplies.
    `manifest` is the run's own `run_manifest.json`, which is where a run says
    whether it adopted a donor's surface; passing None treats the run as cold
    started. `allow_superseded` restates that adoption for this resume.
    """
    check_code_bindings()

    nveg = namelist_int(run_dir, "plasim_namelist", "NVEG", 0)
    naqua = namelist_int(run_dir, "plasim_namelist", "NAQUA", 0)
    ndesert = namelist_int(run_dir, "plasim_namelist", "NDESERT", 0)
    newsurf = namelist_int(run_dir, "landmod_namelist", "NEWSURF", 0)

    if newsurf == 2:
        raise RuntimeError(
            f"{run_dir} stages per-cell surface fields and its landmod namelist "
            "sets NEWSURF = 2, which assigns the scalars wsmax, dz0land and "
            "albland over dwmax, dz0clim and dalbcl at every start "
            "(landmod.f90:492-499). That is the substitution CLIM-67 removed, "
            "asked for by namelist. Stage the fields or set NEWSURF = 0; the run "
            "cannot mean both.")
    if naqua or ndesert:
        raise RuntimeError(
            f"{run_dir} sets NAQUA={naqua} NDESERT={ndesert}, which replaces the "
            "land/sea mask and the orography with an aqua or desert planet. The "
            "staged surface fields are not this run's surface, so nothing here "
            "can be checked against them.")

    records = read_records(restart)

    # Which surface this run is entitled to be carrying. A run prepared with
    # --superseded-surface-ok adopted its donor's surface and discarded the
    # staged one, so the seed copy -- not the .sra -- is what its restart must
    # still agree with. The declaration is read off the manifest the prepare
    # path stamped rather than inferred from the disagreement itself, which
    # would make every disagreement its own excuse.
    override = ((manifest or {}).get("initial_state") or {}).get(
        "superseded_surface_override")
    seed = run_dir / "MOST_REST.seed"
    if override and not seed.is_file():
        raise RuntimeError(
            f"{run_dir} declares superseded_surface_override "
            f"({override}), so its surface came from a donor restart and not "
            f"from the staged .sra -- and {seed.name} is gone, so nothing on "
            "disk says what that surface was. There is no reference to check "
            "this restart against.")
    seed_records = read_records(seed) if override else {}

    verdicts, failures = [], []
    superseded = []
    for code in sorted(codes):
        staged = staged_path(run_dir, code)
        if code in NOT_COMPARABLE:
            verdicts.append({"code": code, "verdict": "not_comparable",
                             "reason": NOT_COMPARABLE[code]})
            continue
        if code in REREAD_EVERY_START:
            name = REREAD_EVERY_START[code]
            if name in records:
                failures.append(
                    f"code {code} ({name}) has become a restart record. It was "
                    "re-read from its .sra at every start, which is why nothing "
                    "could substitute it; that is no longer true and this check "
                    "has to compare it instead of asserting its absence.")
                verdicts.append({"code": code, "verdict": "unexpected_restart_record",
                                 "record": name})
            else:
                verdicts.append({"code": code, "verdict": "reread_every_start",
                                 "record": name})
            continue
        field = SURFACE_RESTART_FIELDS.get(code)
        if field is None:
            failures.append(
                f"code {code} is staged by this project and this check does not "
                "know which restart record it lands in, so nothing verifies that "
                "a resume keeps it. Add it to SURFACE_RESTART_FIELDS, or to "
                "REREAD_EVERY_START if the model re-reads it at every start.")
            verdicts.append({"code": code, "verdict": "unmapped"})
            continue
        if nveg == 2 and code in VEGETATION_OWNED:
            verdicts.append({
                "code": code, "record": field.record, "verdict": "model_evolved",
                "reason": ("NVEG = 2, so simba.f90 overwrites this field every "
                           "timestep and the staged file is not its right answer")})
            continue
        if staged is None:
            failures.append(
                f"code {code} has no staged N???_surf_{code:04d}.sra in "
                f"{run_dir}, so there is nothing to compare the restart's "
                f"{field.record} against.")
            verdicts.append({"code": code, "verdict": "staged_file_missing"})
            continue
        if field.record not in records:
            failures.append(
                f"code {code} is staged but the restart carries no record named "
                f"{field.record}. Either the field is not being written, or the "
                "restart layout has changed and this check is reading for a name "
                "that no longer exists.")
            verdicts.append({"code": code, "verdict": "missing_restart_record",
                             "record": field.record})
            continue

        nlat, nlon = _grid_from(staged)
        ncells = nlat * nlon
        from_file = expected_field(field, read_sra(staged, nlat, nlon))
        got = restart_field(records[field.record], field.months, ncells)

        if override:
            if field.record not in seed_records:
                failures.append(
                    f"code {code}: this run adopted a donor surface and the "
                    f"seed copy carries no {field.record}, so the field it is "
                    "using came from nowhere this run can name.")
                verdicts.append({"code": code, "verdict": "seed_record_missing",
                                 "record": field.record})
                continue
            want = restart_field(seed_records[field.record], field.months, ncells)
            reference = "seed"
        else:
            want = from_file
            reference = "staged"

        differing = int((want != got).sum())
        row = {"code": code, "record": field.record, "months": field.months,
               "cells": ncells, "reference": reference, "differing": differing,
               "binarised": field.binarise}
        if differing:
            delta = np.abs(want - got)
            row.update({
                "verdict": "mismatch",
                "max_abs": float(delta.max()),
                "restart_unique_values": int(np.unique(got).size),
                "reference_unique_values": int(np.unique(want).size),
                "restart_range": [float(got.min()), float(got.max())],
                "reference_range": [float(want.min()), float(want.max())],
            })
            source = ("the staged field" if reference == "staged"
                      else "the surface this run adopted from its donor")
            failures.append(
                f"code {code}: the restart's {field.record} is not {source}. "
                f"{differing} of {want.size} values differ, worst by "
                f"{delta.max():.6g}. The restart holds "
                f"{np.unique(got).size} distinct values over "
                f"[{got.min():.6g}, {got.max():.6g}]; the reference holds "
                f"{np.unique(want).size} over "
                f"[{want.min():.6g}, {want.max():.6g}].")
        else:
            row["verdict"] = "match"

        # The second invariant, and it is separate on purpose: even when the
        # restart still agrees with what the run adopted, the .sra sitting in
        # the run directory may be a different field entirely, because
        # stage_surface_extras rewrites the current one there on every resume.
        # The run directory and the manifest then assert a surface the model
        # will not read, which is the appearance of correctness this module
        # exists to remove.
        if reference == "seed":
            adrift = int((from_file != got).sum())
            row["staged_file_differing"] = adrift
            if adrift:
                superseded.append(
                    f"code {code} ({field.record}): {adrift} of "
                    f"{from_file.size} values in "
                    f"{staged.name} are not what the restart carries. The "
                    f"staged file holds {np.unique(from_file).size} distinct "
                    f"values over [{from_file.min():.6g}, "
                    f"{from_file.max():.6g}]; the model will integrate "
                    f"{np.unique(got).size} over [{got.min():.6g}, "
                    f"{got.max():.6g}].")
        verdicts.append(row)

    report = {
        "restart": restart.name,
        "nveg": nveg,
        "newsurf": newsurf,
        "reference": "seed" if override else "staged",
        "superseded_surface_override": override,
        "superseded_surface_restated": bool(superseded) and allow_superseded,
        "codes": verdicts,
        "matched": sum(1 for v in verdicts if v["verdict"] == "match"),
    }
    if failures:
        raise RuntimeError(
            f"{restart} does not carry the surface fields this run was built "
            "from, so resuming it would integrate a surface nothing staged:\n  "
            + "\n  ".join(failures)
            + "\n\nA restart is only as good as the boundary fields inside it. "
            "See exoplasim/README.md, 'Which resume path is valid', and "
            "TASKS.md CLIM-67 and CLIM-70.")
    if superseded and not allow_superseded:
        raise RuntimeError(
            f"{run_dir} was prepared with --superseded-surface-ok "
            f"({override}), so it integrates its donor's surface and not the "
            "one staged in its own directory. The staged files have since been "
            "rewritten and no longer describe what the model reads:\n  "
            + "\n  ".join(superseded)
            + "\n\nResuming is still legitimate -- a paired A/B measures the "
            "DIFFERENCE between arms and every arm carries the same superseded "
            "surface -- but it is a claim about this segment and not one the "
            "prepare inherited. Pass --superseded-surface-ok to restate it and "
            "stamp it on the segment, or start a run on the current surface. "
            "It is NOT valid for anything feeding the canonical chain.")
    if superseded:
        print("  NOTE: this segment integrates a SUPERSEDED surface, restated "
              "with --superseded-surface-ok and stamped on the segment:")
        for line in superseded:
            print(f"    {line}")
    return report


def staged_path(run_dir: Path, code: int) -> Path | None:
    """The staged `.sra` for one code, found by code rather than by resolution.

    `stage_surface_extras` names these `N<NLAT>_surf_<code>.sra` and `surfmod`
    reopens them by the same rule (`code_surf_file`), so the code identifies the
    file and the NLAT in the name is the model's own, not something to be
    inferred from a config here.
    """
    found = set(run_dir.glob(f"N???_surf_{code:04d}.sra"))
    if not found:
        return None
    if len(found) > 1:
        raise RuntimeError(
            f"{run_dir} holds {len(found)} files for surface code {code}: "
            f"{sorted(p.name for p in found)}. They are different resolutions "
            "and only one of them is this run's, so which field the model read "
            "cannot be decided here.")
    return found.pop()


def _grid_from(staged: Path) -> tuple[int, int]:
    """(nlat, nlon) out of the `.sra` header, which is where the model reads it.

    Header fields 5 and 6 are NLON and NLAT (`sra.write_sra`,
    `surfmod.f90:135`). Taken from the file rather than from the config so a
    resolution mismatch reads as a mismatch instead of a reshape error.
    """
    header = staged.read_text(encoding="ascii").split("\n", 1)[0].split()
    return int(header[5]), int(header[4])


def overwrite_record(raw: bytes, name: str, payload: bytes) -> bytes:
    """One restart record's payload replaced, same length, everything else kept.

    Only the self-test uses this, and only on a copy in a temporary directory. A
    restart is the sole record of where a run was, so nothing here writes into
    `exoplasim/runs/`.
    """
    out = bytearray(raw)
    current = None
    for start, _end, body in walk_records(raw):
        if len(body) == 16 and body.strip() and all(32 <= b < 127 for b in body):
            current = body.decode("ascii").strip()
            continue
        if current == name:
            if len(payload) != len(body):
                raise ValueError(f"{name} is {len(body)} bytes, not {len(payload)}")
            out[start + 4:start + 4 + len(body)] = payload
            return bytes(out)
        current = None
    raise KeyError(name)


def self_test() -> int:
    """Prove the check can fail, on the substitution it was written for.

    CLIM-67's block ran only in the restart branch, so it left no mark on the
    restarts this project happens to have: every run on disk was seeded from a
    donor that never staged code 229, and the block assigned `wsmax` over a
    `dwmax` that was already `wsmax`. The defect is therefore invisible in the
    artifacts, and a check that only ever ran against them would be untested
    against the thing it exists for.

    So the restart the defect WOULD have produced is built here, out of a real
    restart and the real staged fields, and the check is required to accept the
    one and refuse the other. Both cases are named in advance:

      A  the fields as a cold-started run writes them -> must PASS
      B  dwmax flattened to the scalar wsmax, which is what the removed
         `limitwater` block did on every resume -> must FAIL, naming dwmax
      C  dalbcl flattened to the scalar albland, the same substitution on a
         14-month array, which is where dwcl's latent instance sits -> must
         FAIL, naming dalbcl

    Nothing is compiled and no orbit is integrated: this reads restarts and
    `.sra` files that are already on disk and writes its copies to a temporary
    directory.
    """
    import shutil
    import tempfile

    fixtures = sorted(d.name for d in RUNS.glob("run_*")
                      if staged_path(d, 229) is not None
                      and latest_restart(d) is not None)
    run_dir = RUNS / fixtures[0] if fixtures else None
    if run_dir is None:
        print("self-test: no run directory with a staged code 229 and a "
              "restart, so there is nothing to build a fixture from")
        return 1
    restart = latest_restart(run_dir)
    codes = {172, 173, 174, 175, 176, 212, 229}

    with tempfile.TemporaryDirectory(prefix="restart_surface_selftest_") as tmp:
        work = Path(tmp)
        for name in ("plasim_namelist", "landmod_namelist"):
            if (run_dir / name).is_file():
                shutil.copyfile(run_dir / name, work / name)
        raw = restart.read_bytes()
        for code in sorted(codes):
            field = SURFACE_RESTART_FIELDS[code]
            staged = staged_path(run_dir, code)
            shutil.copyfile(staged, work / staged.name)
            nlat, nlon = _grid_from(staged)
            want = expected_field(field, read_sra(staged, nlat, nlon))
            raw = overwrite_record(raw, field.record,
                                   want.astype("<f8").tobytes())
        cold = work / "MOST_REST.00000"
        cold.write_bytes(raw)

        failures = []
        try:
            report = verify_restart_surface_fields(work, cold, codes)
            print(f"  A cold-start restart: PASS "
                  f"({report['matched']} of {len(report['codes'])} codes match)")
        except RuntimeError as exc:
            failures.append(f"A should have passed and did not: {exc}")

        # WSMAX_EARTH, the scalar the removed block assigned. Taken from the
        # value the restarts on disk actually carry rather than from a literal,
        # so the fixture is the substitution this project would have suffered.
        cells = 64 * 128
        flat = np.full(cells, 0.5, dtype="<f8")
        b = work / "MOST_REST.00001"
        b.write_bytes(overwrite_record(raw, "dwmax", flat.tobytes()))
        try:
            verify_restart_surface_fields(work, b, codes)
            failures.append("B should have failed and did not: a dwmax "
                            "flattened to a scalar was accepted")
        except RuntimeError as exc:
            if "dwmax" not in str(exc):
                failures.append(f"B failed without naming dwmax: {exc}")
            else:
                print("  B dwmax flattened to wsmax: REFUSED, naming dwmax")

        c = work / "MOST_REST.00002"
        c.write_bytes(overwrite_record(
            raw, "dalbcl",
            np.full(14 * cells, 0.22, dtype="<f8").tobytes()))
        try:
            verify_restart_surface_fields(work, c, codes)
            failures.append("C should have failed and did not: a dalbcl "
                            "flattened to a scalar was accepted")
        except RuntimeError as exc:
            if "dalbcl" not in str(exc):
                failures.append(f"C failed without naming dalbcl: {exc}")
            else:
                print("  C dalbcl flattened to albland: REFUSED, naming dalbcl")

    if failures:
        for line in failures:
            print(f"  SELF-TEST FAILED: {line}")
        return 1
    print(f"  self-test passed, on a fixture built from {run_dir.name}/"
          f"{restart.name}")
    return 0


def latest_restart(run_dir: Path) -> Path | None:
    """The numbered restart a resume would start from.

    Chosen by the highest ORBIT INDEX the run wrote, not by whatever sorts last:
    the indices are enumerated and required to be contiguous from zero, the same
    property `continue_exoplasim.output_years` asserts before it names a
    restart. A gap means the run directory is not a run and nothing here should
    guess which file it meant.
    """
    orbits = {int(p.suffix[1:]) for p in run_dir.glob("MOST_REST.*")
              if p.suffix[1:].isdigit()}
    if not orbits:
        return None
    if orbits != set(range(max(orbits) + 1)):
        raise RuntimeError(
            f"{run_dir} has non-contiguous restarts {sorted(orbits)}; name one "
            "with --restart rather than having it inferred.")
    return run_dir / f"MOST_REST.{max(orbits):05d}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dir", type=Path, nargs="?", help="a run directory")
    ap.add_argument("--all", action="store_true",
                    help=f"every run under {RUNS}")
    ap.add_argument("--restart", type=Path,
                    help="a specific restart, default the highest-numbered one")
    ap.add_argument("--json", type=Path, help="write the full report here")
    ap.add_argument("--superseded-surface-ok", action="store_true",
                    help="restate that a run adopted its donor's surface, so a "
                         "staged .sra that no longer describes what the model "
                         "reads is reported rather than refused")
    ap.add_argument("--self-test", action="store_true",
                    help="build the restart CLIM-67 would have written and "
                         "require this check to refuse it")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import yaml
    from _paths import CONFIG
    from run_exoplasim import intended_surface_codes

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    codes = intended_surface_codes(config)

    if args.all:
        run_dirs = sorted(p for p in RUNS.glob("run_*") if p.is_dir())
    elif args.run_dir:
        run_dirs = [args.run_dir]
    else:
        ap.error("give a run directory or --all")

    reports, bad = [], 0
    for run_dir in run_dirs:
        # Anchor ONLY an explicit --restart, and only when it is a bare name.
        # `latest_restart` globs inside run_dir and already returns a path
        # anchored there, so re-joining it produced run_dir/run_dir/MOST_REST.*
        # and every run directory failed on a path that does not exist.
        if args.restart is not None:
            restart = args.restart
            if not restart.is_absolute() and not restart.is_file():
                restart = run_dir / restart
        else:
            restart = latest_restart(run_dir)
        if restart is None:
            print(f"{run_dir.name}: no restart")
            continue
        manifest_path = run_dir / "run_manifest.json"
        manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                    if manifest_path.is_file() else None)
        try:
            report = verify_restart_surface_fields(
                run_dir, restart, codes, manifest=manifest,
                allow_superseded=args.superseded_surface_ok)
            print(f"{run_dir.name} {restart.name}: "
                  f"{report['matched']} of {len(report['codes'])} codes match "
                  f"the {report['reference']} reference")
            reports.append({"run": run_dir.name, "ok": True, **report})
        except RuntimeError as exc:
            bad += 1
            print(f"{run_dir.name} {restart.name}: FAILED\n{exc}\n")
            reports.append({"run": run_dir.name, "ok": False, "error": str(exc)})

    if args.json:
        args.json.write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
        try:
            shown = args.json.relative_to(PROJECT_ROOT)
        except ValueError:
            shown = args.json
        print(f"wrote {shown}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
