"""Framing for an ExoPlaSim restart: the one parser and writer.

Worldbuilding frame: the restart is the saved state of the Vesper climate
model. Nothing here is about the real world.

WHY THIS EXISTS. Five scripts under `exoplasim/scripts/` had grown their own
copy of this loop -- `compare_restarts.py`, `diff_restarts.py`,
`scale_restart.py`, `reset_restart_accumulators.py` and `restart_surface.py` --
and they differ in what they validate and in what they call a name. A converter
that rewrites a file rather than reading it cannot be built on a reader that
guesses, so the framing is settled once, here, and the rest import it.

THE FILE. A gfortran sequential unformatted stream. Every logical quantity is
two physical records: a 16-byte space-padded ASCII name, then its payload.
Every physical record is a four-byte little-endian length, the payload, and the
same length again. There is no type tag, no shape, no version and no
configuration fingerprint anywhere in the file, so nothing here can tell an
integer array from a four-byte real array of the same length: that is the
schema's job, not the parser's.

THE RULE THE OTHER READERS BREAK. They decide a record is a name by testing
whether it is 16 bytes of printable ASCII. Names and payloads strictly
alternate -- `put_restart_*` writes them in pairs and nothing else writes to
the unit -- so the position settles it and no payload can ever be mistaken for
a name. `zsolars` is a 16-byte payload and is one bit pattern away from being
read as the name of the record after it.
"""
from __future__ import annotations

import os
import struct
from dataclasses import dataclass
from pathlib import Path

# The name field is `character (len=16)` in `restartmod.f90`, space padded.
NAME_BYTES = 16
# Sequential unformatted record markers. gfortran writes four-byte markers on
# this host; another width or byte order is a file this parser has never seen
# and must not guess at.
MARKER = struct.Struct("<i")


class RestartFormatError(Exception):
    """The bytes are not a restart written by this model."""


@dataclass(frozen=True)
class Record:
    """One named quantity: its name, its payload, and where it came from."""

    name: str
    payload: bytes
    # Byte offset of the NAME record's leading marker, for error messages.
    offset: int

    @property
    def nbytes(self) -> int:
        return len(self.payload)


def _read_physical(raw: bytes, pos: int, path) -> tuple[bytes, int]:
    """One physical record at `pos`; returns its payload and the next offset."""
    if pos + MARKER.size > len(raw):
        raise RestartFormatError(
            f"{path}: truncated at byte {pos}: a record marker needs "
            f"{MARKER.size} bytes and {len(raw) - pos} remain")
    (n,) = MARKER.unpack_from(raw, pos)
    if n < 0:
        raise RestartFormatError(
            f"{path}: record at byte {pos} declares a negative length {n}. "
            "A four-byte marker read from a file written with eight-byte "
            "markers, or with the other byte order, looks like this.")
    end = pos + MARKER.size + n
    if end + MARKER.size > len(raw):
        raise RestartFormatError(
            f"{path}: record at byte {pos} declares {n} bytes but only "
            f"{len(raw) - pos - MARKER.size} remain; the file is truncated")
    (m,) = MARKER.unpack_from(raw, end)
    if m != n:
        raise RestartFormatError(
            f"{path}: record markers disagree at byte {pos} ({n} then {m}). "
            "This is not a sequential unformatted file written by this model.")
    return raw[pos + MARKER.size:end], end + MARKER.size


def read(path: Path) -> list[Record]:
    """Every named quantity in file order.

    Raises rather than repairing. A restart that does not parse is not a
    restart to convert from: the caller wanted a state and there is no state
    here to be had.
    """
    path = Path(path)
    raw = path.read_bytes()
    if not raw:
        raise RestartFormatError(f"{path}: empty file")

    records: list[Record] = []
    seen: set[str] = set()
    pos = 0
    while pos < len(raw):
        start = pos
        name_bytes, pos = _read_physical(raw, pos, path)
        if len(name_bytes) != NAME_BYTES:
            raise RestartFormatError(
                f"{path}: expected a {NAME_BYTES}-byte name at byte {start} "
                f"and found a {len(name_bytes)}-byte record. Names and "
                "payloads alternate; an odd count means a record was written "
                "without its name or the file is a splice.")
        try:
            name = name_bytes.decode("ascii").strip()
        except UnicodeDecodeError as exc:
            raise RestartFormatError(
                f"{path}: the name record at byte {start} is not ASCII "
                f"({exc}). Names are written as `character (len=16)`.") from exc
        if not name:
            raise RestartFormatError(
                f"{path}: the name record at byte {start} is blank")
        if pos >= len(raw):
            raise RestartFormatError(
                f"{path}: '{name}' at byte {start} has a name and no payload; "
                "the file ends mid-pair")
        payload, pos = _read_physical(raw, pos, path)
        if name in seen:
            raise RestartFormatError(
                f"{path}: '{name}' appears twice. `reseek` returns the FIRST "
                "match, so a duplicate is a state the model would silently "
                "half-read.")
        seen.add(name)
        records.append(Record(name=name, payload=payload, offset=start))

    if not records:
        raise RestartFormatError(f"{path}: no records")
    return records


def encode(records) -> bytes:
    """The file bytes for a sequence of records."""
    out = bytearray()
    for rec in records:
        name = rec.name.encode("ascii")
        if len(name) > NAME_BYTES:
            raise RestartFormatError(
                f"'{rec.name}' is {len(name)} bytes; the model's name field "
                f"is {NAME_BYTES} and would truncate it into a collision")
        for payload in (name.ljust(NAME_BYTES, b" "), bytes(rec.payload)):
            marker = MARKER.pack(len(payload))
            out += marker + payload + marker
    return bytes(out)


def write(path: Path, records, *, overwrite: bool = False) -> None:
    """Write records to `path` through a temporary file, then rename.

    A converter that fails half way through must not leave something that
    parses. The rename is atomic on this filesystem, so the output either does
    not exist or is complete.
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise RestartFormatError(
            f"{path} exists. Pass an explicit force option rather than "
            "overwriting a state something may have been run from.")
    raw = encode(records)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_bytes(raw)
    # Read it back through the parser before it takes the output's name: the
    # writer and the reader disagreeing is exactly the failure a rename cannot
    # protect against.
    read(tmp)
    os.replace(tmp, path)


def index(records) -> dict:
    """Records by name. `read` has already refused duplicates."""
    return {rec.name: rec for rec in records}
