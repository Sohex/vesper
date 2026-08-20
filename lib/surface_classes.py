"""Read the derived surface classes by NAME, never by integer code.

`pedology/scripts/build_surface_classes.py` writes `surface_cover` and
`duricrust` as integer codes and, being CF-conventional about it, writes the
legend alongside as `flag_values` and `flag_meanings`. So the mapping from name
to code is IN THE FILE and no consumer has to carry a copy of it.

That matters because the alternative is failure-modes class 1 -- one quantity,
several consumers, and a fix that reaches some of them. A consumer holding
`PAVEMENT = 3` keeps working silently if the class list ever gains a member
ahead of it, and it keeps working while pointing at the wrong class. Asking the
file removes the possibility rather than documenting it.
"""

from __future__ import annotations

from pathlib import Path

from netCDF4 import Dataset
import numpy as np


def cover_legend(path: Path, variable: str = "surface_cover") -> dict:
    """{name: code} for a classified field, from its own flag attributes."""
    with Dataset(path) as ds:
        if variable not in ds.variables:
            raise KeyError(f"{path} has no variable {variable}")
        v = ds[variable]
        try:
            names = str(v.flag_meanings).split()
            values = [int(x) for x in np.asarray(v.flag_values).ravel()]
        except AttributeError as exc:
            raise KeyError(
                f"{variable} in {path} carries no flag_values/flag_meanings "
                "legend; it predates the CF attributes and its codes cannot be "
                "resolved by name") from exc
    if len(names) != len(values):
        raise ValueError(
            f"{variable} in {path} has {len(values)} flag_values and "
            f"{len(names)} flag_meanings; the legend is not usable")
    return dict(zip(names, values))


def cover_mask(path: Path, name: str, variable: str = "surface_cover"):
    """Boolean mask over regions where the classified field is `name`.

    Raises on an unknown name rather than returning an empty mask, because a
    misspelled class silently selecting nothing is the shape of error that gets
    read as "this world has no pavement".
    """
    legend = cover_legend(path, variable)
    if name not in legend:
        raise KeyError(
            f"{variable} in {path} has no class {name!r}; it carries "
            f"{sorted(legend)}")
    with Dataset(path) as ds:
        field = np.asarray(ds[variable][:]).astype(np.int64)
    return field == legend[name]
