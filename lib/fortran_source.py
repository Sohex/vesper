"""Reading a Fortran source the way a check has to read it.

WORLDBUILDING CONTEXT: Vesper is an invented planet and the sources these
functions read are the climate and ocean models that simulate it. Nothing here
is about the real world.

WHAT THIS IS FOR. Two gates in this tree state, in a declaration file, what a
compiled model's source has to still say -- `exoplasim/config/ocean_tier.yaml`
for the adopted atmosphere's slab, `config/cgenie_calibration.yaml` for the
adopted offline ocean's calibration set. Both have to answer the same four
questions of a Fortran file: what does this line say once the comments and the
continuations are gone, what does its initialiser evaluate to, and is this
symbol in a namelist block or is it not. Those answers are here so the two
gates cannot drift into two dialects of the same reader; a symbol that one gate
calls namelist-reachable and the other does not is a disagreement about the
model, and it must not be reachable from a difference in two regexes.

Fortran is case-insensitive, so everything here lowercases: `TMELT` and `tmelt`
are one symbol and a reader that told them apart would be checking the typing
rather than the model.
"""

from __future__ import annotations

import re

# A namelist block's body: names separated by commas and nothing else. The
# commas are what ends it -- the statement after the block begins with a token
# that is not preceded by one, and a looser pattern swallowed that token and
# lost the last key in the block with it.
NAME_LIST = r"[a-z_][a-z0-9_]*(?:\s*,\s*[a-z_][a-z0-9_]*)*"


def strip_comments(text: str) -> str:
    """Fortran source with `!` comments removed, line by line.

    Column-1 `!` and trailing `!` are the same comment marker in free and fixed
    form alike in the sources this reads. Where a `!` can appear inside a
    character literal, the caller matches the literal BEFORE this runs -- which
    is what both gates' refusal checks do.
    """
    out = []
    for line in text.split("\n"):
        cut = line.find("!")
        out.append(line if cut < 0 else line[:cut])
    return "\n".join(out)


def strip_fixed_form_comments(text: str) -> str:
    """The same, plus fixed-form comments: `c`, `C` or `*` in column 1.

    GOLDSTEIN and its sea ice are fixed-form `.F` and their commentary is almost
    all column-1 `c`. Without this a commented-out assignment reads as a live
    one, which is the exact failure a declaration of what the source RUNS has to
    avoid: `initialise_goldstein.F` carries `c     sodaylen = 86400.0` three
    lines above the print that reports the value it actually holds.
    """
    out = []
    for line in strip_comments(text).split("\n"):
        if line[:1] in ("c", "C", "*"):
            continue
        out.append(line)
    return "\n".join(out)


def join_continuations(text: str) -> str:
    """Free-form continuation lines joined into the statement they belong to."""
    return re.sub(r"&\s*\n\s*&?", " ", text)


def join_fixed_continuations(text: str) -> str:
    """Fixed-form continuations joined: any non-blank in column 6 continues.

    Fixed-form Fortran continues a statement by putting a character in column 6
    of the next line, and GOLDSTEIN uses digits, `$`, `+`, `:` and `&` for it
    interchangeably. A reader that only knew `&` would see the wind-stress
    scaling's four assignments as eight fragments.
    """
    out: list[str] = []
    for line in text.split("\n"):
        body = line[6:] if len(line) > 6 else ""
        marker = line[5:6] if len(line) > 5 else ""
        if out and marker.strip() and not line[:5].strip():
            out[-1] = out[-1] + " " + body
        else:
            out.append(line)
    return "\n".join(out)


def squash(text: str) -> str:
    """Comments gone, continuations joined, whitespace collapsed, lowercased."""
    return re.sub(r"\s+", " ", join_continuations(strip_comments(text))).lower()


def tight(text: str) -> str:
    """`squash` with every space removed as well.

    Used wherever the thing being matched is a Fortran STATEMENT, where
    whitespace carries nothing: a declaration, an assignment, an expression.
    `squash` is kept for messages, where the spaces are inside a character
    literal and are part of what is being matched.
    """
    return squash(text).replace(" ", "")


def squash_fixed(text: str) -> str:
    """`squash` for fixed-form source: column-1 comments and column-6 continuations."""
    joined = join_fixed_continuations(strip_fixed_form_comments(text))
    return re.sub(r"\s+", " ", join_continuations(joined)).lower()


def tight_fixed(text: str) -> str:
    """`squash_fixed` with every space removed as well."""
    return squash_fixed(text).replace(" ", "")


def initialiser(line: str) -> float | None:
    """The value a declaration line's own initialiser evaluates to.

    Both shapes the sources use: `parameter(CLFSN = 3.337E5)` and
    `real :: CPS = 3990.34`. The arithmetic is the source's own, so a value
    changed in the model and a value changed in a declaration file fail the
    same check.
    """
    text = squash(line).strip()
    m = re.search(r"parameter\s*\(\s*[a-z_][a-z0-9_]*\s*=\s*(.+?)\s*\)\s*$", text)
    if not m:
        m = re.search(r"::\s*[a-z_][a-z0-9_]*\s*(\([^)]*\))?\s*=\s*(.+?)\s*$", text)
        if not m:
            return None
        expression = m.group(2)
    else:
        expression = m.group(1)
    expression = expression.replace("d", "e").replace("D", "e")
    if not re.fullmatch(r"[0-9e+\-.]+", expression):
        return None
    try:
        return float(expression)
    except ValueError:
        return None


def namelist_symbols(text: str, *, fixed: bool = False) -> set[str]:
    """Every symbol in every `namelist/<name>/` block of one source file.

    A key is namelist-reachable exactly when it is in one of these, so this is
    what a declared reachability is checked against rather than a list written
    down twice.
    """
    names: set[str] = set()
    for block in namelist_blocks(text, fixed=fixed).values():
        names |= block
    return names


def namelist_blocks(text: str, *, fixed: bool = False) -> dict[str, set[str]]:
    """The same, per namelist name, so a key can be checked against ITS block.

    `fixed` picks the reading, and it is a choice rather than a union because
    the two readings are not both safe on one file. Column 1 means a comment in
    fixed form and means an identifier in free form, so reading a `.f90` as
    fixed would drop every declaration whose name begins with c -- and reading a
    `.F` as free would leave its namelist split across the continuation lines it
    is stated over. Taking both would let a mangled reading ADD a name, and an
    added name reports a compile-time constant as settable, which is the one
    error this column exists to prevent.
    """
    flat = squash_fixed(text) if fixed else squash(text)
    out: dict[str, set[str]] = {}
    for m in re.finditer(
            r"namelist\s*/\s*([a-z_][a-z0-9_]*)\s*/\s*(.+?)(?=namelist\s*/|$)", flat):
        run = re.match(NAME_LIST, m.group(2))
        if run:
            out.setdefault(m.group(1), set()).update(
                n.strip() for n in run.group(0).split(",") if n.strip())
    return out
