"""Write an ExoPlaSim surface field in its `.sra` text format.

This lived in `convert_orogen.py`, which is superseded and now deleted. Four
current scripts imported it from there, so a live utility sat inside a module
whose docstring said "do not use" -- which reads, correctly, as current code
depending on dead code.

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
