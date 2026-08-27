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
continents was flattened to one number after the first orbit. CLIM-67
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
a scatter, and the matching gather at write time (`mpimod_omp.f90:684`). None of
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
* Code 129's record is `groundsg` rather than the `doro` that `surfmod.f90`'s
  `surfcode` table binds the code to. `SurfaceRestartField.binds` carries that
  name so `check_code_bindings` still reads the model's own table for it.

## What is deliberately NOT compared, and why

A field the model is entitled to change is not a field the staged file is the
right answer for, and comparing one would be a diagnostic wearing a test's
clothes. Each exclusion is reported by name, never silently dropped.

* Code 129 `doro`, whose restart record is NOT the one the staged file is the
  right answer for. `doro` is the orography the dynamics is built on and it is
  `groundoro + glacieroro`, so it legitimately moves as ice accumulates. The
  LITHOGRAPHIC half is a separate record, `groundsg`, and that one is the staged
  field exactly: `glacierini` reads the `.sra` into `doro` on a cold start and
  copies it to `groundoro` unchanged, and no arithmetic but one multiplication
  by `oroscale` ever touches it (`glaciermod.f90:181,182`). Nothing spectral does:
  the fit goes into the spectral array `so`, and the gridpoint round trip that
  would overwrite `doro` with it is inside a `npro == 1` print block the threaded
  build never enters.

  So code 129 is compared against `groundsg`, exactly, and `doro` is held to the
  identity `groundsg + dglacsg` instead of to a file. Together those say the
  whole of the resumed orography: the lithographic half is the field that was
  staged, and the rest is this run's own ice. This is the check CLIM-70 declared
  it could not make, and it needed the second record rather than the spectral fit
  reproduced. It is withdrawn, by name and with the reason, under the one setting
  that takes the right answer away: `OROSCALE` other than 1.0, which makes the
  record the staged field times the scale rather than the staged field.

  `NGLACIER` does not take it away and used to withdraw it. Both halves hold at
  every setting: `groundoro` is assigned from `doro` before the nglacier switch
  and nowhere else, and under NGLACIER = 0 `oroini` never runs, `glacieroro`
  stays zero and the identity is true term by term. At that setting the identity
  is uninformative rather than wrong, and the `groundsg` comparison carries the
  whole of the check on its own. world-zq1k.
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

The restart record walker is `restart_format`; the `.sra` reader is
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
from paths import rel  # noqa: E402  from lib/, put on sys.path by _paths
import restart_format  # noqa: E402
from sra import read_sra  # noqa: E402


@dataclass(frozen=True)
class SurfaceRestartField:
    """One staged surface code and the restart record it must survive as.

    `record` is the restart record name, `months` its second dimension there,
    and `expand` how a single-record `.sra` becomes that many months. `binarise`
    marks the land/sea mask, whose right answer is the thresholded file rather
    than the file.

    `binds` is the name `surfmod.f90`'s own `surfcode` table gives the code, for
    the codes whose restart record is named something else. It is empty when the
    two agree, which is every code but 129: the `.sra` fills `doro` and the
    lithographic half of it is written back as `groundsg`.
    """

    code: int
    record: str
    months: int = 1
    binarise: bool = False
    binds: str = ""


# Every code `intended_surface_codes` can return, with the restart record it
# lands in. A code absent from here is reported as unmapped rather than passed:
# adding a generator without saying what happens to its field on a resume is the
# gap this whole module is about.
SURFACE_RESTART_FIELDS = {
    f.code: f
    for f in (
        SurfaceRestartField(129, "groundsg", binds="doro"),
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
    # The saturated soil albedo pair, PHYS-15. `landmod.f90:833-835` reads these
    # through `mpsurfgp` on every start and says why in its own comment: the
    # saturated endmember is a BOUNDARY CONDITION and not a state, so it cannot
    # drift, a run cannot have changed it, and carrying it would be three more
    # restart records for a field the surface file already holds. The dry pair
    # beside it, 174/175/176, IS carried, which is the distinction this entry
    # records. `landmod.f90:836` aborts when the field is absent under
    # nwetsoil = 1, so a resume with nothing staged stops rather than mixing
    # toward the sentinel.
    1742: "dalbwet",
    1750: "dalbwet1",
    1760: "dalbwet2",
}

# Fields SIMBA owns once coupled vegetation is on: `simba.f90:484,487` assigns
# both every timestep under `nveg == 2`, so the staged file stops being their
# right answer the moment such a run starts. Named as a set rather than tested
# by literal, because anything else simba starts writing belongs here beside
# them.
VEGETATION_OWNED = {212, 229}

def not_comparable(oroscale: float) -> dict[int, str]:
    """Codes whose right answer this run's settings take away, with the reason.

    Reported per code rather than dropped, so an exclusion is a verdict a reader
    meets rather than an absence they have to notice. There is one, and it is a
    setting: under `oroscale = 1.0` code 129 has a right answer and is compared,
    and under anything else it does not.

    NGLACIER DOES NOT TAKE CODE 129'S RIGHT ANSWER AWAY, and it used to withdraw
    the code whenever it was not 1. world-zq1k. `glacierini`'s cold branch sets
    `groundoro(:) = doro(:)` BEFORE the nglacier switch (`glaciermod.f90:182`),
    and the only other write to `groundoro` in the whole model is the restart
    read itself, so `groundsg` is the staged field times `oroscale` under every
    setting of NGLACIER. The other half holds too: under NGLACIER = 0 `oroini`
    is never called at all, `glacieroro` stays at the zero it is declared with,
    `doro` stays the field `glacierini` read from the `.sra`, and
    `doro == groundsg + dglacsg` is true term by term. The identity is
    UNINFORMATIVE at that setting rather than wrong, and an uninformative check
    that passes is not a reason to stop making it.

    `nglacier` is still reported beside the verdicts, because which of the two
    halves is carrying the information depends on it.
    """
    out: dict[int, str] = {}
    if oroscale != 1.0:
        out[129] = (
            f"OROSCALE = {oroscale}, and glacierini scales the staged field by "
            "it where that field enters the model (glaciermod.f90:181), so the "
            "restart's groundsg is the staged field TIMES oroscale and the "
            "staged file is not its right answer. Scaled exactly once, as "
            "world-6qee settled; the record is well defined, it is the .sra "
            "that is not the thing to compare it with")
    return out


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
        name = field.binds or field.record
        actual = bindings.get(name)
        if actual != field.code:
            wrong.append(f"{name} is code {actual} in surfmod.f90, "
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


def namelist_float(run_dir: Path, filename: str, key: str, default: float) -> float:
    """One real out of a staged namelist, or `default` if it is unset.

    The default matters here in a way it does not for the integers: OROSCALE is
    a `planet_nl` key nothing in this project writes, so an absent key is the
    normal case and `p_earth.f90`'s compiled 1.0 is what the model will use.
    Passing that value in is what makes reading the file mean something.
    """
    path = run_dir / filename
    if not path.is_file():
        return default
    match = re.search(rf"^\s*{key}\s*=\s*(-?[\d.]+(?:[eEdD][-+]?\d+)?)",
                      path.read_text(encoding="latin-1"),
                      re.IGNORECASE | re.MULTILINE)
    return float(match.group(1).replace("d", "e").replace("D", "E")) if match else default


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
    # THE ONE SETTING CODE 129'S RIGHT ANSWER DEPENDS ON, read from the run for
    # the reason `namelist_int` gives: `oroscale` decides whether `groundsg`
    # still equals the field that was staged. `nglacier` is read to be REPORTED
    # beside the verdicts -- it says which half of the check is informative, not
    # whether either half has an answer. world-zq1k.
    nglacier = namelist_int(run_dir, "glacier_namelist", "NGLACIER", 0)
    oroscale = namelist_float(run_dir, "planet_namelist", "OROSCALE", 1.0)
    excluded = not_comparable(oroscale)

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

    records = restart_format.payloads(restart)

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
    seed_records = restart_format.payloads(seed) if override else {}

    verdicts, failures = [], []
    superseded = []
    for code in sorted(codes):
        staged = staged_path(run_dir, code)
        if code in excluded:
            verdicts.append({"code": code, "verdict": "not_comparable",
                             "reason": excluded[code]})
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

    # THE SECOND HALF OF CODE 129, and it is what makes checking `groundsg`
    # stand for the topography rather than for one record. `doro` is the
    # orography the dynamics is built on, and it is not the staged field: it is
    # the lithographic half plus the ice sheet. The identity `oroini` writes it
    # by is exact -- one addition in the same precision, `glaciermod.f90:473`,
    # and nothing scales it afterwards since world-6qee moved oroscale to the
    # cold read -- so `doro == groundsg + dglacsg` is a right answer with no
    # tolerance at any oroscale, and
    # holding it means the whole of `doro`'s departure from the staged file is
    # the model's own ice and nothing else. Without it a substituted `doro`
    # would pass on a `groundsg` that was never read.
    #
    # Held at every NGLACIER. Under 0 the identity is the trivial one -- `oroini`
    # never runs, `glacieroro` stays at its declared zero and `doro` stays the
    # field `glacierini` read -- and a right answer that is trivially true is
    # still the answer. It goes red exactly where it should: on a `doro` that
    # came from somewhere neither record accounts for.
    glacier_identity = None
    if 129 in codes and 129 not in excluded:
        missing = [n for n in ("doro", "groundsg", "dglacsg") if n not in records]
        if missing:
            failures.append(
                "the restart carries groundsg but not " + ", ".join(missing) +
                ", so the orography it will resume on cannot be decomposed into "
                "the staged field and the model's ice.")
        else:
            ncells = len(records["groundsg"]) // 8
            oro = restart_field(records["doro"], 1, ncells)
            ground = restart_field(records["groundsg"], 1, ncells)
            ice = restart_field(records["dglacsg"], 1, ncells)
            differing = int((oro != ground + ice).sum())
            glacier_identity = {"record": "doro", "cells": ncells,
                                "differing": differing,
                                "verdict": "match" if not differing else "mismatch"}
            if differing:
                delta = np.abs(oro - (ground + ice))
                glacier_identity["max_abs"] = float(delta.max())
                failures.append(
                    f"the restart's doro is not groundsg + dglacsg: "
                    f"{differing} of {oro.size} values differ, worst by "
                    f"{delta.max():.6g}. The orography the model will resume on "
                    "is therefore neither the staged code 129 nor that field "
                    "plus this run's ice, so nothing on disk says where it came "
                    "from.")

    report = {
        "restart": restart.name,
        "nveg": nveg,
        "newsurf": newsurf,
        "nglacier": nglacier,
        "oroscale": oroscale,
        "glacier_orography_identity": glacier_identity,
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
            "CLIM-67 and CLIM-70.")
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
    records = restart_format.decode(raw)
    for i, rec in enumerate(records):
        if rec.name != name:
            continue
        if len(payload) != rec.nbytes:
            raise ValueError(f"{name} is {rec.nbytes} bytes, not {len(payload)}")
        records[i] = restart_format.Record(name=name, payload=payload,
                                           offset=rec.offset)
        return restart_format.encode(records)
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
      D  a donor surface whose adoption is not restated -> must FAIL
      E  the same, restated -> must PASS against the seed
      F  groundsg flattened to its own mean, a restart that re-derived its
         topography instead of carrying the staged 129 -> must FAIL, naming
         groundsg
      G  doro moved off groundsg + dglacsg, an orography substituted where every
         per-code comparison still passes -> must FAIL, naming doro

    Each case names the surface code it needs and reports itself skipped when
    the fixture does not stage it, rather than the fixture deciding whether the
    function runs at all.

    Nothing is compiled and no orbit is integrated: this reads restarts and
    `.sra` files that are already on disk and writes its copies to a temporary
    directory.
    """
    import shutil
    import tempfile

    # THE FIXTURE IS CHOSEN BY WHAT IT STAGES, not by one code. Requiring a
    # staged 229 made the whole self-test return early on a tree whose only runs
    # stage the seven codes this project currently generates, 229 not among
    # them: every case below was skipped and the function reported a failure
    # that was about the fixture rather than about the check. Each case now
    # states the code it needs and says so when the fixture cannot supply it, so
    # what runs is reported and what cannot is named.
    def _staged_codes(d: Path) -> set[int]:
        return {c for c in SURFACE_RESTART_FIELDS
                if staged_path(d, c) is not None}

    candidates = [d for d in sorted(RUNS.glob("run_*"))
                  if d.is_dir() and latest_restart(d) is not None
                  and _staged_codes(d)]
    if not candidates:
        print("self-test: no run directory with a restart and a staged surface "
              "code, so there is nothing to build a fixture from")
        return 1
    run_dir = max(candidates, key=lambda d: (len(_staged_codes(d)),
                                             len(list(d.glob("MOST_REST.0*"))),
                                             d.name))
    restart = latest_restart(run_dir)
    # 129 is left out here and owned by F/G below: its record is one of three
    # that have to stay mutually consistent, so substituting it in the shared
    # fixture would break the identity the other cases are not about.
    codes = _staged_codes(run_dir) - {129}

    with tempfile.TemporaryDirectory(prefix="restart_surface_selftest_") as tmp:
        work = Path(tmp)
        for name in ("plasim_namelist", "landmod_namelist",
                     "glacier_namelist", "planet_namelist"):
            if (run_dir / name).is_file():
                shutil.copyfile(run_dir / name, work / name)
        raw = restart.read_bytes()
        cells, seed_field, seed_want = None, None, None
        for code in sorted(codes):
            field = SURFACE_RESTART_FIELDS[code]
            staged = staged_path(run_dir, code)
            shutil.copyfile(staged, work / staged.name)
            nlat, nlon = _grid_from(staged)
            cells = nlat * nlon
            want = expected_field(field, read_sra(staged, nlat, nlon))
            # ANY COMPARED FIELD WILL DO for the donor case: what D and E need
            # is one record where the seed and the staged file disagree, not a
            # particular one. It was pinned to dwmax, so a fixture that does not
            # stage code 229 skipped the whole donor branch.
            if seed_field is None or field.record == "dwmax":
                seed_field, seed_want = field, want
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
        # The CELL COUNT comes from the fixture's own staged file, as every
        # other length in this function does. It was `64 * 128`, so a donor at
        # any rung but T21 made `overwrite_record` raise outside the try and
        # crash the self-test instead of reporting a verdict.
        flat = np.full(cells, 0.5, dtype="<f8")
        b = work / "MOST_REST.00001"
        if 229 not in codes:
            print("  B skipped: the fixture stages no code 229")
        else:
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
        if 174 not in codes:
            print("  C skipped: the fixture stages no code 174")
        else:
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
                    print("  C dalbcl flattened to albland: REFUSED, naming "
                          "dalbcl")

        # D AND E: THE DONOR-SURFACE BRANCH, which is the half of this function
        # the production call reaches through `manifest=` and
        # `allow_superseded=`. The self-test passed neither, so `override` was
        # always None and the code that decides whether the staged .sra or the
        # donor's seed copy is the authority was never executed by it -- a check
        # built from a shorter argument list than the call it stands for.
        # world-60x0.
        #
        # The seed differs from the staged file in exactly one record and the
        # restart carries the seed, which is the state a run prepared with
        # --superseded-surface-ok is in. D is that state UNRESTATED and must be
        # refused; E is the same state restated and must pass AGAINST THE SEED,
        # which is the assertion that says the branch ran at all.
        if seed_field is None:
            failures.append("D/E could not be built: the fixture stages no "
                            "comparable code, so there is no record to make the "
                            "seed differ in")
        else:
            seed_raw = overwrite_record(
                raw, seed_field.record,
                (np.asarray(seed_want, dtype="<f8") + 0.25).tobytes())
            (work / "MOST_REST.seed").write_bytes(seed_raw)
            d = work / "MOST_REST.00003"
            d.write_bytes(seed_raw)
            donor = {"initial_state":
                     {"superseded_surface_override": "run_selftest_donor"}}
            try:
                verify_restart_surface_fields(work, d, codes, manifest=donor)
                failures.append(
                    "D should have failed and did not: a run whose staged .sra "
                    "no longer describes the donor surface it integrates was "
                    "accepted without the adoption restated")
            except RuntimeError as exc:
                if "superseded-surface-ok" not in str(exc):
                    failures.append(f"D failed for the wrong reason: {exc}")
                else:
                    print("  D donor surface, adoption not restated: REFUSED")
            try:
                rep = verify_restart_surface_fields(
                    work, d, codes, manifest=donor, allow_superseded=True)
                if rep["reference"] != "seed":
                    failures.append(
                        f"E was checked against the {rep['reference']!r} "
                        "reference rather than the seed, so the donor-surface "
                        "branch did not run and D proved nothing about it")
                else:
                    print("  E donor surface, adoption restated: PASS against "
                          "the seed")
            except RuntimeError as exc:
                failures.append(f"E should have passed and did not: {exc}")

        # F AND G: THE OROGRAPHY, the pair CLIM-72 added. They are built from the
        # fixture run's own restart rather than from the synthetic one above,
        # because the orography records are three and the identity between them
        # has to survive intact for F to mean that groundsg alone was refused.
        #
        #   F  groundsg flattened to its own mean, which is a restart that
        #      re-derived its topography instead of carrying the staged field
        #      -> must FAIL, naming groundsg
        #   G  doro moved off groundsg + dglacsg, which is a substituted
        #      orography that leaves both other records untouched and would pass
        #      every per-code comparison -> must FAIL, naming doro
        oro_staged = staged_path(run_dir, 129)
        oro_records = restart_format.payloads(restart)
        if oro_staged is None or "groundsg" not in oro_records:
            print("  F/G skipped: the fixture stages no code 129 or its restart "
                  "carries no groundsg")
        else:
            shutil.copyfile(oro_staged, work / oro_staged.name)
            nlat, nlon = _grid_from(oro_staged)
            oro_cells = nlat * nlon
            oro_codes = codes | {129}
            base = restart.read_bytes()
            for name in ("groundsg", "dglacsg", "doro"):
                base = overwrite_record(base, name, oro_records[name])
            ok = work / "MOST_REST.00004"
            ok.write_bytes(base)
            try:
                rep = verify_restart_surface_fields(work, ok, {129})
                if rep["glacier_orography_identity"] is None:
                    failures.append(
                        "F's baseline reported no orography identity, so the "
                        "doro check did not run and G proves nothing")
            except RuntimeError as exc:
                failures.append(
                    f"F's baseline should have passed and did not: {exc}")

            ground = np.frombuffer(oro_records["groundsg"], dtype="<f8")
            f = work / "MOST_REST.00005"
            f.write_bytes(overwrite_record(
                base, "groundsg",
                np.full(oro_cells, float(ground.mean()), dtype="<f8").tobytes()))
            try:
                verify_restart_surface_fields(work, f, oro_codes)
                failures.append("F should have failed and did not: a groundsg "
                                "flattened to its mean was accepted")
            except RuntimeError as exc:
                if "groundsg" not in str(exc):
                    failures.append(f"F failed without naming groundsg: {exc}")
                else:
                    print("  F groundsg flattened to its mean: REFUSED, naming "
                          "groundsg")

            oro = np.frombuffer(oro_records["doro"], dtype="<f8")
            g = work / "MOST_REST.00006"
            g.write_bytes(overwrite_record(
                base, "doro", (oro + 1.0).tobytes()))
            try:
                verify_restart_surface_fields(work, g, oro_codes)
                failures.append("G should have failed and did not: a doro that "
                                "is not groundsg + dglacsg was accepted")
            except RuntimeError as exc:
                if "doro" not in str(exc):
                    failures.append(f"G failed without naming doro: {exc}")
                else:
                    print("  G doro off the glacier identity: REFUSED, naming "
                          "doro")

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
        print(f"wrote {rel(args.json)}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
