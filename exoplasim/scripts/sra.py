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


def write_sra(path: Path, code: int, field: np.ndarray) -> None:
    """Write one field to `path` under ExoPlaSim's numeric `code`."""
    nlat, nlon = field.shape
    flat = np.asarray(field, dtype=np.float64).ravel(order="C")
    if flat.size % 8:
        raise ValueError(
            f"SRA field size must be divisible by 8; {nlat}x{nlon} is "
            f"{flat.size}, which would write a short final row and shift "
            f"everything ExoPlaSim reads after it")
    header = [code, 0, SRA_DATE_STAMP, 0, nlon, nlat, 0, 0]
    with path.open("w", encoding="ascii") as handle:
        handle.write("".join(f" {value:11d}" for value in header) + "\n")
        for row in flat.reshape(-1, 8):
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
