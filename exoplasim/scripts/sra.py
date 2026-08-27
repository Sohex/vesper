"""Write an ExoPlaSim surface field in its `.sra` text format.
## The format

Eight integers of header, then the field as rows of eight values, latitude-major
from the north. The field size must divide by eight because the format has no
way to express a partial row; ExoPlaSim reads it back with a fixed-width Fortran
format and a short row shifts everything after it.

The header's third field is a date stamp. ExoPlaSim does not interpret it, and
it is fixed rather than taken from the clock so that regenerating an unchanged
surface produces a byte-identical file. That matters here: the boundary
conditions are hashed into run manifests, and a timestamp would make every
rebuild look like a change.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Fixed, not the current date. See the module docstring: a clock here would make
# every regenerated surface field hash differently for no physical reason.
SRA_DATE_STAMP = 20260811



def refuse_a_write_through_a_symlink(path: Path) -> None:
    """Refuse to write a staged field that is a link into another checkout.

    `scripts/link_worktree.py` links `exoplasim/inputs/<rung>/` ENTRY BY ENTRY,
    because a worktree has to READ the staged fields to prepare a run. The
    directory is also regenerable output, and the four builders that write it
    resolve their default output to that same directory, so a generator run in
    a worktree opens a symlink and writes STRAIGHT THROUGH into the main
    checkout -- replacing the fields a commissioning run there is reading, with
    nothing recording which worktree did it. `check_consistency.py` would then
    report that run as being on a surface the tree no longer holds, which is
    the true statement and not the useful one.

    A SKIP IN `link_worktree.py` IS THE WRONG REPAIR, because the read is
    legitimate and a worktree that cannot stage a run cannot do anything. The
    repair is a refusal at the point of writing, which is here: `write_sra` is
    the single door all four builders go through, so one guard covers
    `build_surface_albedo.py`, `build_surface_soil_water.py`,
    `build_boundary_conditions.py` and `build_surface_roughness.py` at once.

    It fires on the path itself or on any symlinked ancestor, and it says what
    to pass instead. In the main checkout nothing on these paths is a link, so
    it never fires there. world-4o59.
    """
    path = Path(path)
    offender, target = None, None
    for candidate in (path, *path.parents):
        if candidate.is_symlink():
            offender, target = candidate, candidate.resolve()
            break
        if candidate == candidate.parent:
            break
    if offender is None:
        return
    raise SystemExit(
        f"{path} is reached through the symlink {offender}, which points at "
        f"{target}.\n"
        "  Writing here would write THROUGH the link into another checkout and "
        "replace a staged surface field a run there may be reading, and nothing "
        "would record which worktree did it.\n"
        "  A worktree links these entries so it can READ them to stage a run; "
        "generating into them is what is refused. Pass --output to a path "
        "inside this worktree, or run the generator in the main checkout. "
        "world-4o59."
    )

def write_sra(path: Path, code: int, field: np.ndarray) -> None:
    """Write one field to `path` under ExoPlaSim's numeric `code`.

    `field` is (nlat, nlon), or (nlev, nlat, nlon) for a code the model reads
    with more than one level. `surfmod`'s `get_surf_array` reads `klot`
    consecutive header-and-field records out of one file, which is how the
    14-month climatologies arrive; a levelled field is written the same way,
    every record under the same code. A 2-D field writes exactly one record, so
    this is what it always was for every existing caller.
    """
    refuse_a_write_through_a_symlink(path)
    levels = np.asarray(field, dtype=np.float64)
    if levels.ndim == 2:
        levels = levels[None, :, :]
    if levels.ndim != 3:
        raise ValueError(
            f"SRA field must be (nlat, nlon) or (nlev, nlat, nlon); got shape "
            f"{np.shape(field)}")
    nlat, nlon = levels.shape[1], levels.shape[2]
    if (nlat * nlon) % 8:
        raise ValueError(
            f"SRA field size must be divisible by 8; {nlat}x{nlon} is "
            f"{nlat * nlon}, which would write a short final row and shift "
            f"everything ExoPlaSim reads after it")
    header = [code, 0, SRA_DATE_STAMP, 0, nlon, nlat, 0, 0]
    with path.open("w", encoding="ascii") as handle:
        for level in levels:
            handle.write("".join(f" {value:11d}" for value in header) + "\n")
            for row in level.ravel(order="C").reshape(-1, 8):
                handle.write("".join(f" {value:12.5f}" for value in row) + "\n")


def read_sra(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """Read one field back out of an `.sra`, as (nlat, nlon).

    The inverse of `write_sra` and the reader for surface fields this component
    generated earlier in the chain -- background albedo, roughness, soil water --
    when a later step needs the field it actually gave the model rather than a
    reconstruction of it.

    The header is one line of eight integers; everything after it is the field.
    The size is asserted rather than inferred, because a wrong resolution reads
    as a reshape error only when the total happens not to divide.
    """
    lines = path.read_text(encoding="ascii").splitlines()
    values = np.array(" ".join(lines[1:]).split(), dtype=np.float64)
    if values.size != nlat * nlon:
        raise ValueError(
            f"{path}: {values.size} values, expected {nlat * nlon} for "
            f"{nlat}x{nlon}")
    return values.reshape(nlat, nlon)


def read_sra_records(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """Every record in a multi-record `.sra`, as (records, nlat, nlon).

    `read_sra` above asserts a single record, which is right for the fields
    that have one and is why a multi-record file fails it loudly rather than
    reshaping into nonsense. The layered soil-water capacity split is written
    one record per water layer, so it needs this.

    The record length is derived from the grid and the file's own line count
    rather than from a declared layer count, so a file with a different number
    of layers than the caller expected is a size error here and not a silent
    truncation.
    """
    lines = path.read_text(encoding="ascii").splitlines()
    per_record = 1 + -(-(nlat * nlon) // 8)
    if len(lines) % per_record:
        raise ValueError(
            f"{path}: {len(lines)} lines is not a whole number of "
            f"{per_record}-line records for {nlat}x{nlon}")
    records = len(lines) // per_record
    out = np.empty((records, nlat, nlon))
    for r in range(records):
        block = lines[r * per_record + 1:(r + 1) * per_record]
        out[r] = np.array(" ".join(block).split(),
                          dtype=np.float64).reshape(nlat, nlon)
    return out
