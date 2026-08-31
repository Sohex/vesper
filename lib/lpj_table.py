"""One column-wise reader for an LPJ-GUESS annual output table.

An LPJ-GUESS `.out` table is whitespace-delimited text whose first three columns
are Lon, Lat and Year and whose remaining columns are the simulated per-field
values, one row per gridcell and simulated year. A run BIO-14 assesses carries
millions of such rows in each of its outputs, and the acceptance gate reads
every row once for the merged table and once more in the rank piece it came
from, so the same data is parsed about seventeen times over BY DESIGN: checking
every rank is what the gate is for.

This module parses each column in one pass and states every check its callers
make as a vectorised test over those columns. It certifies ONLY the ordinary
case. Anything it cannot certify -- a header that is not Lon Lat Year, a ragged
row, a nonnumeric or nonfinite field, a nonintegral year, a duplicated
cell-year, or a coordinate close enough to a half-cent boundary that the
column-wise rounding could disagree with Python's -- raises `RowParseRequired`,
and the caller falls back to its own row-at-a-time parser. That is what keeps
the refusal messages exactly what they were: the fast path never diagnoses a
defect, it only declines to certify one.

The rounding is why the near-tie test exists. Cell identity is `round(lon, 2)`,
which is correctly-rounded decimal, while a column-wise multiply by a hundred
and rint is binary arithmetic that can land on the other side of a half-cent
tie. The test refuses the whole table when any coordinate sits within a
thousandth of a cent of a tie, so on a table this module does certify, the two
agree row for row.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Cell identity is the coordinate rounded to two decimal places, so a table is
# keyed on hundredths of a degree held as integers. The spans bound a whole
# sphere; the year span bounds a signed 32-bit year. Both are FIXED rather than
# taken from the table in hand, because two tables' keys are compared against
# each other and a packing that depended on either table's range would not be
# the same packing.
CENTS = 100
LON_CENTS = 180 * CENTS
LAT_CENTS = 90 * CENTS
LAT_SPAN = 2 * LAT_CENTS + 1
YEAR_OFFSET = 1 << 31
YEAR_SPAN = 1 << 32
# How far from a half-cent tie every coordinate must sit for the column-wise
# rounding to be certainly what Python's round would give. A coordinate written
# with two decimals is many orders of magnitude clear of this.
TIE_MARGIN = 1.0e-3


class RowParseRequired(Exception):
    """The column-wise parse cannot certify this table; parse it a row at a time."""


@dataclass(frozen=True)
class Table:
    """One certified output table, its rows in (lon, lat, year) order.

    `key` is strictly increasing, so the rows are unique and the table is
    already in the order a merge check compares two tables in.
    """

    fields: list[str]
    lon: np.ndarray
    lat: np.ndarray
    year: np.ndarray
    values: np.ndarray
    key: np.ndarray

    @property
    def rows(self) -> int:
        return int(self.values.shape[0])

    def cell_starts(self) -> np.ndarray:
        """The row index at which each distinct cell's block begins."""
        if self.rows == 0:
            return np.zeros(0, dtype=np.int64)
        change = np.flatnonzero((self.lon[1:] != self.lon[:-1])
                                | (self.lat[1:] != self.lat[:-1])) + 1
        return np.concatenate((np.zeros(1, dtype=change.dtype), change))

    def cells(self) -> list[tuple[float, float]]:
        """The distinct cells, ascending in (lon, lat)."""
        start = self.cell_starts()
        return [(float(lon), float(lat))
                for lon, lat in zip(self.lon[start], self.lat[start])]

    def years(self) -> list[int]:
        """The distinct years, ascending."""
        return [int(year) for year in np.unique(self.year)]

    def grid(self, ncells: int, nyears: int) -> np.ndarray:
        """The values as (cell, year, field), for a complete cell-year product.

        The caller establishes completeness first; `coverage` is what does that
        in the acceptance gate. Sorted order already groups a cell's rows and
        orders them by year, so this is a reshape and never a scatter.
        """
        if self.rows != ncells * nyears:
            raise ValueError("table is not a complete cell-year product")
        return self.values.reshape(ncells, nyears, len(self.fields))


def composite_key(lon_cents: np.ndarray, lat_cents: np.ndarray,
                  year: np.ndarray) -> np.ndarray:
    """One int64 per row, ordering rows exactly as (lon, lat, year) does.

    A single integer key is what lets a table be sorted, deduplicated and
    compared against another table in one pass each, rather than by a
    three-column lexical sort per operation.
    """
    cell = (lon_cents + LON_CENTS) * LAT_SPAN + (lat_cents + LAT_CENTS)
    return cell * YEAR_SPAN + (year + YEAR_OFFSET)


def header_fields(line: str) -> list[str] | None:
    """The value field names, or None if this is not an LPJ output header."""
    header = line.split()
    if len(header) < 4 or header[:3] != ["Lon", "Lat", "Year"]:
        return None
    return header[3:]


def _cents(column: np.ndarray, limit: int) -> np.ndarray:
    """A coordinate column in hundredths of a degree, or a refusal to certify."""
    scaled = column * CENTS
    cents = np.rint(scaled)
    if not np.all(np.abs(scaled - cents) < 0.5 - TIE_MARGIN):
        raise RowParseRequired("a coordinate sits too close to a half-cent tie")
    if not np.all(np.abs(cents) <= limit):
        raise RowParseRequired("a coordinate leaves the sphere")
    return cents.astype(np.int64)


def _rounded(cents: np.ndarray) -> np.ndarray:
    """Hundredths of a degree back to the float `round(x, 2)` returns.

    Python's round is correctly-rounded decimal, so the float for a given
    hundredth is taken from Python once per distinct coordinate and broadcast,
    never recomputed in binary over the column.
    """
    if cents.size == 0:
        return cents.astype(np.float64)
    low = int(cents.min())
    offsets = (cents - low).astype(np.intp)
    lookup = np.zeros(int(offsets.max()) + 1, dtype=np.float64)
    for offset in np.flatnonzero(np.bincount(offsets)):
        lookup[offset] = round((int(offset) + low) / CENTS, 2)
    return lookup[offsets]


def read(path: Path) -> Table:
    """Parse and certify one output table column-wise.

    Raises `RowParseRequired` for anything the column-wise checks cannot
    certify, including an empty or headerless file, so that the caller's own
    row-at-a-time parser is what produces every refusal message.
    """
    with path.open("rb") as handle:
        first = handle.readline()
    if not first:
        raise RowParseRequired("empty file")
    try:
        fields = header_fields(first.decode("utf-8"))
    except UnicodeDecodeError as exc:
        raise RowParseRequired("undecodable header") from exc
    if fields is None:
        raise RowParseRequired("header is not Lon Lat Year")
    width = len(fields) + 3
    try:
        with warnings.catch_warnings():
            # A table with a header and no rows is one this reader declines and
            # the caller then refuses in its own words, so numpy's notice of it
            # is not the gate's verdict and does not belong in the gate's output.
            warnings.filterwarnings("ignore", message="loadtxt: input contained no data")
            block = np.loadtxt(path, skiprows=1, dtype=np.float64,
                               comments=None, ndmin=2)
    except (ValueError, UnicodeDecodeError) as exc:
        raise RowParseRequired("column-wise parse failed") from exc
    if block.shape[0] == 0 or block.shape[1] != width:
        raise RowParseRequired("no data rows, or a width the header disagrees with")

    values = block[:, 3:]
    if not np.isfinite(block[:, :3]).all() or not np.isfinite(values).all():
        raise RowParseRequired("a nonfinite field")
    year_float = block[:, 2]
    if not np.all(np.abs(year_float) < YEAR_OFFSET):
        raise RowParseRequired("a year outside the packed range")
    year = year_float.astype(np.int64)
    if not np.array_equal(year.astype(np.float64), year_float):
        raise RowParseRequired("a nonintegral year")
    lon_cents = _cents(block[:, 0], LON_CENTS)
    lat_cents = _cents(block[:, 1], LAT_CENTS)

    key = composite_key(lon_cents, lat_cents, year)
    order = np.argsort(key, kind="stable")
    key = key[order]
    if key.size > 1 and np.any(key[1:] == key[:-1]):
        raise RowParseRequired("a duplicated cell-year")

    return Table(fields=fields, lon=_rounded(lon_cents[order]),
                 lat=_rounded(lat_cents[order]), year=year[order],
                 values=np.ascontiguousarray(values[order]), key=key)
