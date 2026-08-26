#!/usr/bin/env python3
"""Automatic locals whose ONLY definition is inside a WHERE, and are read anyway.

    python exoplasim/scripts/lint_masked_locals.py          # the unclassified
    python exoplasim/scripts/lint_masked_locals.py --all     # the whole population
    python exoplasim/scripts/lint_masked_locals.py --tsv     # machine-readable

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. A procedure declares an array local with no initialiser, so it
is automatic storage and its contents on entry are whatever the stack held. It
then writes that array only inside `where (mask)` blocks, and reads it. A
`where` selects which elements the ASSIGNMENT stores; it does not restrict which
elements of a right-hand side the compiler EVALUATES, and a vectorising
compilation evaluates all of them. So on every lane the mask excluded, the read
is a read of indeterminate memory -- and the declared flag line carries
`-ffpe-trap=invalid,zero,overflow` (`config/planet.yaml`), which turns whatever
the stack held into SIGFPE at the instruction that touched it.

That is not hypothetical. `tands` wrote nine soil locals only under
`where (dls(:) > 0.0)` and then divided by three of them, so on every sea lane
the divisor was a stack word; `mktsoil` had the same shape for four locals and
`mkdca` for two. world-d016 fixed all three, and
`notes/audits/masked-where-blocks.md` carries the reading.

The second face of the same defect does not need a trap at all. Where the read
is under a DIFFERENT mask than the write, or under no mask, the indeterminate
lane is not discarded -- it is stored, and the model carries a stack word
forward as a value.

WHY THIS IS NOT A RULE INSIDE `lint_implicit_save.py`. The two are adjacent and
they are exact complements, which is precisely why one pass cannot hold both.
`lint_implicit_save` fires on the PRESENCE of an initialiser in a procedure-body
declaration: an initialiser confers SAVE, SAVE is static storage, and under
`-fopenmp` static storage is one copy for the whole thread team. This pass fires
on the ABSENCE of one. Merging them would mean a single pass that reports a
declaration both for having an initialiser and for not having it, and the
remedies point in opposite directions: `lint_implicit_save`'s remedy is to
delete the initialiser and preset the array in the body, which is exactly the
precondition this class needs. A merged lint would flag its own remedy.

They also read different things. `lint_implicit_save` is a declaration-only
pass and never looks at an executable statement; this one is decided entirely
in the body, by where the writes sit relative to the `where` stack, which is
`lint_masked_domains.py`'s machinery. Neither half of this pass is a rule the
other lint could carry without becoming the other lint.

The three passes divide the same source three ways and the boundaries are
crisp: `lint_implicit_save` asks what STORAGE CLASS a declaration has,
`lint_masked_domains` asks whether an EXPRESSION inside a mask has a domain the
mask was standing in for, and this one asks whether a local has a DEFINED VALUE
on the lanes a mask discarded.

THE PASS OVER-REPORTS, DELIBERATELY, and that direction is what makes a clean
run evidence. Five places it errs towards reporting a local that is fine, and
none towards clearing one that is not:

  1. An element write inside an `if` does not clear the local. The lane axis
     is dimension 1 (NHOR) and a `where` masks that axis and nothing else, so
     `z(:) = ...` covers every lane and `z(jhor) = ...` covers one. An element
     write whose enclosing `do` is unconditional is taken to cover the axis;
     one under an `if` is not, because the condition selects lanes. That is the
     shape `if (dls(jhor) > 0.0) z(jhor) = ...` takes, and it is how `tands`'s
     `zsnowz` stayed reachable.
  2. A `where` / `elsewhere` pair that writes the same local in both arms does
     define every lane, and this pass does not pair the arms. Both writes are
     masked writes, so the local is reported and the pairing is the argument in
     its CLASSIFIED row.
  3. Every read counts, wherever it sits. A read that only ever reaches a
     discarded lane of a further masked assignment is still a read.
  4. The where-stack is popped only by an explicit `end where` or a procedure
     boundary, as in `lint_masked_domains`. A construct the parser fails to
     close makes more writes look masked, never fewer.
  5. An `allocatable` local is treated as an ordinary automatic one. Its
     contents after `allocate` are equally indeterminate.

WHAT IT CANNOT SEE, stated because a direction claim is worthless without it.
A whole-axis write at `where`-depth zero clears the local for good, so a write
inside an `if` block that may not be taken, or one covering only some elements
of a trailing dimension, reads as a definition here. Those are a different
class -- an unconditional read of a conditionally written local -- and
`notes/audits/uninitialised-reads-and-implicit-save.md` carries it. A local
defined through a `call` that writes an `intent(out)` dummy is not tracked
either. This pass answers one question: is every definition of this local
inside a mask.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import MODEL_SRC  # noqa: E402

SRC = MODEL_SRC / "plasim" / "src"

_TYPES = (
    r"double\s+precision",
    r"double\s+complex",
    r"character",
    r"integer",
    r"logical",
    r"complex",
    r"real",
    r"type\s*\(",
    r"class\s*\(",
)
_DECL_RE = re.compile(r"^\s*(" + "|".join(_TYPES) + r")\b", re.I)
_DIMENSION_STMT_RE = re.compile(r"^\s*dimension\b", re.I)
_PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:real|integer|logical|complex|double\s+precision|character)"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+([a-z_]\w*)",
    re.I,
)
_PROC_END_RE = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
_WHERE_OPEN_RE = re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?where\s*\(", re.I)
_END_WHERE_RE = re.compile(r"^\s*end\s*where\b", re.I)
_ELSEWHERE_RE = re.compile(r"^\s*else\s*where\b", re.I)
_ONE_LINE_IF_RE = re.compile(r"^\s*if\s*\(", re.I)

# Statement keywords that never open an assignment. `forall` and `associate`
# are not in the model source; if one appears, its body statements are read as
# ordinary assignments, which is the over-reporting direction only when the
# construct is a mask, and `forall` is not one.
_NOT_ASSIGNMENTS = (
    "do", "if", "else", "end", "enddo", "endif", "call", "write", "read",
    "print", "goto", "go", "return", "allocate", "deallocate", "nullify",
    "use", "implicit", "parameter", "common", "data", "save", "namelist",
    "format", "continue", "stop", "cycle", "exit", "select", "case",
    "interface", "contains", "entry", "external", "intrinsic", "equivalence",
    "open", "close", "rewind", "backspace", "inquire", "where", "elsewhere",
    "module", "subroutine", "function", "program", "then", "while",
    "dimension", "allocatable", "pointer", "target", "intent", "optional",
    "public", "private", "sequence", "type", "class", "block", "forall",
    "associate", "include", "pause", "assign",
)

_ID_RE = re.compile(r"[a-z_]\w*", re.I)


# ---------------------------------------------------------------------------
# THE CLASSIFIED POPULATION, one row per local with the reason it is safe.
#
# The key is (file, procedure, variable) and the value is the argument. A row is
# a claim; a row whose reason does not hold is a defect in this table. A local
# that is RENAMED or moved to another procedure drops out of the table and is
# reported again, which is the direction that matters.
#
# `--all` prints the whole population with its verdicts; the default run prints
# only what is NOT here, and exits 1 if anything is.
# Evidence: notes/audits/masked-where-blocks.md.
# ---------------------------------------------------------------------------

CLASSIFIED: dict[tuple[str, str, str], str] = {
    ("landmod.f90", "mkradv", "zrop"):
        "the where/elsewhere pair that writes it has an unconditional "
        "elsewhere, so the two arms partition every lane and zrop is defined "
        "everywhere before mpgagp reads it. This pass does not pair the arms "
        "of a construct, which is over-report 2",
    ("radmod.f90", "lwr", "zaco2"):
        "the where/elsewhere pair that writes it has an unconditional "
        "elsewhere -- zsumco2 <= 1.0 and its complement -- so the two arms "
        "partition every lane and zaco2 is defined everywhere before the "
        "clear-sky transmissivity reads it. Over-report 2, the same shape as "
        "mkradv's zrop. The sibling zah2o is not reported at all because the "
        "continuum adds an unmasked definition below the pair",
    ("radmod.f90", "lwr", "zao3"):
        "the same shape: a where(zsumo3 <= 0.01) / elsewhere pair with an "
        "unconditional elsewhere, so every lane is written before "
        "ztaucs reads it. Over-report 2",
    ("radmod.f90", "lwr", "zth2o"):
        "the same shape: a where(zsumwv <= 2.) / elsewhere pair with an "
        "unconditional elsewhere, so every lane is written before the CO2 "
        "overlap and the CH4 and N2O bands multiply through it. Over-report 2",
    ("oceanmod.f90", "hdiffo", "zdtx"):
        "the three index ranges of zdtx(0:NLON,:) are written between them on "
        "every lane: 1:NLON-1 by a where/elsewhere pair with an unconditional "
        "elsewhere, 0 by a second such pair, and NLON by an unmasked copy of "
        "0. Every read is of an index in 0:NLON, so none is indeterminate. "
        "Over-report 1 and 2 together -- a section write is not a whole-axis "
        "write and the arms are not paired",
}


# ---------------------------------------------------------------------------
# DEFERRED FILES. A file here is REPORTED and does not gate, with the issue
# that carries it. This is not a classification: it is a statement that the
# sites have not been read, and the count is printed separately so it cannot be
# mistaken for a clean one.
#
# The bar for a row is the same one lint_masked_domains.py's radmod remainder
# met: the file is being edited on another branch, so a per-local row written
# now is a claim about a version that may not survive.
#
# The table is EMPTY. radmod.f90 was the one entry: swr's 53 locals are now
# preset on every lane and lwr's three carry rows above, world-px61.
# ---------------------------------------------------------------------------

DEFERRED_FILES: dict[str, str] = {}


def squeeze(text: str) -> str:
    return " ".join(text.split())


def strip_comment(line: str) -> str:
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def logical_lines(text: str):
    """Yield (first_physical_lineno, joined_text) with continuations folded."""
    buf = ""
    start = None
    for n, line in enumerate(text.splitlines(), 1):
        code = strip_comment(line).rstrip()
        if not code.strip():
            continue
        if buf:
            code = code.lstrip()
            if code.startswith("&"):
                code = code[1:]
        else:
            start = n
        if code.rstrip().endswith("&"):
            buf += code.rstrip()[:-1]
            continue
        buf += code
        yield start, buf
        buf = ""
    if buf:
        yield start, buf


def match_paren(text: str, open_idx: int) -> int:
    depth = 0
    quote = None
    for i in range(open_idx, len(text)):
        ch = text[i]
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
    return -1


def split_top(text: str, sep: str = ",") -> list[str]:
    parts = []
    depth = 0
    cur = ""
    quote = None
    for ch in text:
        if quote:
            cur += ch
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            cur += ch
        elif ch in "([":
            depth += 1
            cur += ch
        elif ch in ")]":
            depth -= 1
            cur += ch
        elif ch == sep and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def dummy_args(head: str) -> set[str]:
    """The dummy argument names on a procedure heading, lowercased."""
    if "(" not in head:
        return set()
    close = match_paren(head, head.index("("))
    if close < 0:
        return set()
    inner = head[head.index("(") + 1 : close]
    return {a.strip().lower() for a in split_top(inner) if a.strip()}


def declared_arrays(text: str) -> list[tuple[str, bool]]:
    """(name, is_array) for every entity a type declaration introduces.

    Entities that carry an initialiser are dropped: an initialiser confers SAVE
    and static storage, so the variable is defined on entry and is
    `lint_implicit_save.py`'s class rather than this one.
    """
    head, sep, rhs = text.partition("::")
    if not sep:
        m = _DECL_RE.match(text)
        if not m:
            return []
        head = text[: m.end()]
        # A type specifier may carry its own parenthesised kind or length.
        rest = text[m.end() :]
        if rest.lstrip().startswith("("):
            close = match_paren(rest, rest.index("("))
            if close < 0:
                return []
            head += rest[: close + 1]
            rest = rest[close + 1 :]
        rhs = rest
    attrs = head.lower()
    if "parameter" in attrs or "intent" in attrs or "external" in attrs:
        return []
    shaped_by_attr = "dimension" in attrs
    out = []
    for entity in split_top(rhs):
        entity = entity.strip()
        if not entity or "=" in entity:
            continue
        name = entity.split("(")[0].split("*")[0].strip()
        if not re.fullmatch(r"[a-z_]\w*", name, re.I):
            continue
        out.append((name.lower(), shaped_by_attr or "(" in entity))
    return out


def assignment_of(stmt: str) -> tuple[str, str] | None:
    """(lhs, rhs) if this statement is an assignment, else None."""
    first = _ID_RE.match(stmt.strip())
    if first and first.group(0).lower() in _NOT_ASSIGNMENTS:
        return None
    depth = 0
    quote = None
    for i, ch in enumerate(stmt):
        if quote:
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
        elif ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        elif ch == "=" and depth == 0:
            if stmt[i + 1 : i + 2] == "=":
                return None
            if stmt[i - 1 : i] in ("=", "/", "<", ">"):
                return None
            return stmt[:i], stmt[i + 1 :]
    return None


def covers_lane_axis(lhs: str) -> bool:
    """Does this assignment target define every element of the masked axis?

    The masked axis is dimension 1, the horizontal NHOR index every `where` in
    this model tests. A bare name is the whole array; a subscript whose FIRST
    component is `:` covers that axis; anything else -- a scalar index, a
    partial section -- does not, and the local stays undefined here.
    """
    lhs = lhs.strip()
    if "(" not in lhs:
        return True
    close = match_paren(lhs, lhs.index("("))
    if close < 0:
        return False
    first = split_top(lhs[lhs.index("(") + 1 : close])[0].strip()
    return first == ":"


def strip_prefixes(text: str) -> tuple[str, bool, bool]:
    """Peel a one-line `where (mask)` or `if (cond)` prefix off a statement.

    Returns the statement it governs, whether a `where` prefix was peeled, and
    whether an `if` prefix was. A one-line `if` does not mask lanes -- its
    condition is a scalar -- but it does make the write conditional, which is
    what decides whether an element write covers the lane axis.
    """
    masked = False
    conditional = False
    while True:
        stripped = text.lstrip()
        if _WHERE_OPEN_RE.match(stripped):
            open_idx = stripped.lower().index("(", stripped.lower().index("where"))
            close = match_paren(stripped, open_idx)
            rest = stripped[close + 1 :].strip() if close > 0 else ""
            if not rest:
                return text, masked, conditional
            masked = True
            text = rest
            continue
        if _ONE_LINE_IF_RE.match(stripped):
            close = match_paren(stripped, stripped.index("("))
            rest = stripped[close + 1 :].strip() if close > 0 else ""
            if not rest or rest.lower().startswith("then"):
                return text, masked, conditional
            conditional = True
            text = rest
            continue
        return text, masked, conditional


def scan(filename: str, text: str):
    """Every local in `text` whose definitions are all inside a `where`."""
    sites = []
    proc = None
    args: set[str] = set()
    decls: dict[str, int] = {}
    where_write: dict[str, int] = {}
    clear_write: set[str] = set()
    read: dict[str, int] = {}
    masks: dict[str, str] = {}
    stack: list[str] = []
    ifdepth = 0
    in_module = False

    def flush():
        if proc is None:
            return
        for name, line in decls.items():
            if name in clear_write or name not in where_write:
                continue
            if name not in read:
                continue
            sites.append(
                dict(
                    file=filename,
                    line=line,
                    proc=proc,
                    name=name,
                    write=where_write[name],
                    readline=read[name],
                    mask=masks.get(name, "?")[:80],
                )
            )

    for lineno, raw in logical_lines(text):
        low = raw.lower()
        body = re.sub(r"^\s*\d+\s+", "      ", raw)

        if re.match(r"^\s*module\s+[a-z_]\w*\s*$", low):
            in_module = True
            continue
        if re.match(r"^\s*end\s*module\b", low):
            in_module = False
            continue
        if _PROC_END_RE.match(low):
            flush()
            proc = None
            args = set()
            decls, where_write, clear_write, read, masks = {}, {}, set(), {}, {}
            stack = []
            ifdepth = 0
            continue
        pm = _PROC_START_RE.match(low)
        if pm:
            flush()
            proc = pm.group(2)
            args = dummy_args(raw[pm.end(2) - len(pm.group(2)) :])
            decls, where_write, clear_write, read, masks = {}, {}, set(), {}, {}
            stack = []
            ifdepth = 0
            continue
        if proc is None or in_module:
            continue

        if _END_WHERE_RE.match(low):
            if stack:
                stack.pop()
            continue
        if _ELSEWHERE_RE.match(low):
            m = re.match(r"^\s*else\s*where\s*\(", low)
            if m:
                close = match_paren(body, m.end() - 1)
                if stack:
                    stack[-1] = ".not. " + body[m.end() : close]
            elif stack:
                stack[-1] = ".not. " + stack[-1]
            continue

        if _DECL_RE.match(body) or _DIMENSION_STMT_RE.match(body):
            src = body
            if _DIMENSION_STMT_RE.match(body):
                src = "real :: " + body.strip()[len("dimension") :]
            for name, is_array in declared_arrays(src):
                if name in args or not is_array:
                    continue
                decls.setdefault(name, lineno)
            continue

        if re.match(r"^\s*(save|data|common|equivalence|namelist)\b", low):
            # An explicit SAVE, a `data` initialisation or storage association
            # all give the name a defined value on entry. Drop it from the
            # population rather than reporting it: the storage class is
            # lint_implicit_save.py's question, not this one.
            for name in _ID_RE.findall(body):
                clear_write.add(name.lower())
            continue

        # An `if` block, and a `select case`, make the statements inside them
        # conditional. That does not select LANES the way a `where` does, but it
        # is what decides whether an element write -- `z(jhor) = ...` -- covers
        # the whole lane axis: an unconditional `do jhor` loop does, and the
        # same loop with `if (dls(jhor) > 0.0)` inside it does not.
        if re.match(r"^\s*end\s*(if|select)\b", low) or re.match(
            r"^\s*end(if|select)\b", low
        ):
            ifdepth = max(0, ifdepth - 1)
            continue

        stmt, one_line_mask, one_line_if = strip_prefixes(body)
        opened = None
        if _WHERE_OPEN_RE.match(stmt.lstrip()):
            s = stmt.lstrip()
            open_idx = s.lower().index("(", s.lower().index("where"))
            close = match_paren(s, open_idx)
            opened = s[open_idx + 1 : close] if close > 0 else "?"

        active = list(stack) + ([opened] if opened is not None else [])
        if one_line_mask:
            m = re.match(r"^\s*(?:[a-z_]\w*\s*:\s*)?where\s*\(", body.lstrip(), re.I)
            if m:
                s = body.lstrip()
                close = match_paren(s, m.end() - 1)
                if close > 0:
                    active = list(stack) + [s[m.end() : close]]

        asg = assignment_of(stmt) if opened is None else None
        if asg:
            lhs, rhs = asg
            m = _ID_RE.match(lhs.strip())
            if m:
                lhs_name = m.group(0).lower()
                if lhs_name in decls:
                    if active:
                        where_write.setdefault(lhs_name, lineno)
                        masks.setdefault(lhs_name, squeeze(" .and. ".join(active)))
                    elif covers_lane_axis(lhs) or (ifdepth == 0 and not one_line_if):
                        clear_write.add(lhs_name)
            # Everything but the name being assigned is a read: the peeled
            # prefix carries the mask, and a subscript on the left is an index
            # into the target rather than a definition of the name inside it.
            prefix = body[: len(body) - len(stmt)]
            lhs_sub = lhs[lhs.index("(") :] if "(" in lhs else ""
            scan_text = f"{prefix} {lhs_sub} {rhs}"
        else:
            scan_text = body

        for name in _ID_RE.findall(scan_text):
            low_name = name.lower()
            if low_name in decls:
                read.setdefault(low_name, lineno)

        if opened is not None:
            stack.append(opened)
        if re.search(r"\bthen\s*$", low) or re.match(
            r"^\s*select\s*case\b", low
        ):
            if not re.match(r"^\s*else\s*if\b", low) and not re.match(
                r"^\s*elseif\b", low
            ):
                ifdepth += 1

    flush()
    return sites


# ---------------------------------------------------------------------------
# THE FIXTURES. Reduced sources with the right answer known, run on every
# invocation, because a gate whose only evidence is that the tree passes cannot
# distinguish "no masked-only local" from "no pass". Each entry is
# (name, source, expected number of reported locals).
# ---------------------------------------------------------------------------

FIXTURES: tuple[tuple[str, str, int], ...] = (
    ("a local written only in a where and read is reported", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = 1.0 / z(:)
      endwhere
      end subroutine leaf
""", 1),
    ("a preset on every lane clears it", """
      subroutine leaf
      real z(NHOR)
      z(:) = 1.0
      where (dls(:) > 0.0)
       z(:) = 2.0
       y(:) = 1.0 / z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("a masked-only local nothing reads is not reported", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
      endwhere
      end subroutine leaf
""", 0),
    ("a read under a DIFFERENT mask is still a read", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
      endwhere
      where (dsnow(:) > 0.0)
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 1),
    ("a dummy argument is not a local", """
      subroutine leaf(z)
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("an initialised local is SAVEd, so it is the other lint's class", """
      subroutine leaf
      real :: z(NHOR) = 0.0
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("a scalar local is not of this class, a where cannot assign one", """
      subroutine leaf
      real z
      where (dls(:) > 0.0)
       y(:) = z
      endwhere
      end subroutine leaf
""", 0),
    ("an element write does NOT clear it, the lane axis is not covered", """
      subroutine leaf
      real z(NHOR)
      do jhor = 1, NHOR
       if (dls(jhor) > 0.0) z(jhor) = 1.0
      enddo
      where (dls(:) > 0.0)
       z(:) = 2.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 1),
    ("a whole-lane write in a trailing-dimension loop clears it", """
      subroutine leaf
      real z(NHOR,NLEV)
      do jlev = 1, NLEV
       z(:,jlev) = 1.0
      enddo
      where (dls(:) > 0.0)
       y(:) = z(:,1)
      endwhere
      end subroutine leaf
""", 0),
    ("a one-line where is a mask", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0) z(:) = 1.0
      y(:) = z(:)
      end subroutine leaf
""", 1),
    ("a one-line if is not a mask", """
      subroutine leaf
      real z(NHOR)
      if (nstep > 0) z(:) = 1.0
      where (dls(:) > 0.0)
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("continuations are folded, so a split write is one write", """
      subroutine leaf
      real z(NHOR)
      z(:) = 1.0                                                        &
     &     + 2.0
      where (dls(:) > 0.0)
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("the where stack does not leak past its end", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0)
       y(:) = 1.0
      endwhere
      z(:) = 1.0
      w(:) = z(:)
      end subroutine leaf
""", 0),
    ("a local used only as a subscript is still read", """
      subroutine leaf
      integer k(NHOR)
      where (dls(:) > 0.0)
       k(:) = 1
      endwhere
      y(k(1)) = 2.0
      end subroutine leaf
""", 1),
    ("a module variable is not a procedure local", """
      module m
      real z(NHOR)
      end module m
      subroutine leaf
      use m
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("each procedure is judged on its own body", """
      subroutine one
      real z(NHOR)
      z(:) = 1.0
      end subroutine one
      subroutine two
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine two
""", 1),
    ("two locals in one declaration are two sites", """
      subroutine leaf
      real z(NHOR), w(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
       w(:) = 2.0
       y(:) = z(:) + w(:)
      endwhere
      end subroutine leaf
""", 2),
    ("an elsewhere arm is still a masked write", """
      subroutine leaf
      real z(NHOR)
      where (dls(:) > 0.0)
       z(:) = 1.0
      elsewhere
       z(:) = 2.0
      endwhere
      y(:) = z(:)
      end subroutine leaf
""", 1),
    ("a data statement gives it a value on entry", """
      subroutine leaf
      real z(NHOR)
      data z /NHOR*0.0/
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 0),
    ("an allocatable is treated as automatic", """
      subroutine leaf
      real, allocatable :: z(:)
      allocate(z(NHOR))
      where (dls(:) > 0.0)
       z(:) = 1.0
       y(:) = z(:)
      endwhere
      end subroutine leaf
""", 1),
)


def selftest() -> list[str]:
    bad = []
    for name, source, want in FIXTURES:
        got = len(scan("fixture.f90", source))
        if got != want:
            bad.append(f"{name}: expected {want} reported, got {got}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--all", action="store_true", help="print the whole population with verdicts"
    )
    ap.add_argument("--tsv", action="store_true", help="one tab-separated row per site")
    ap.add_argument(
        "--file", default=None, help="restrict the pass to one source file name"
    )
    args = ap.parse_args()

    broken = selftest()
    if broken:
        print("the pass does not answer its own fixtures, so its verdict on the "
              "tree means nothing:")
        for b in broken:
            print(f"    {b}")
        return 1

    sites = []
    for path in sorted(SRC.glob("*.f90")):
        if args.file and path.name != args.file:
            continue
        sites.extend(scan(path.name, path.read_text(errors="replace")))
    for s in sites:
        s["verdict"] = CLASSIFIED.get((s["file"], s["proc"], s["name"]))

    deferred = [s for s in sites if s["file"] in DEFERRED_FILES]
    if args.all:
        shown = sites
    else:
        shown = [
            s for s in sites
            if not s["verdict"] and s["file"] not in DEFERRED_FILES
        ]

    if args.tsv:
        for s in shown:
            print(
                "\t".join(
                    (
                        s["file"],
                        str(s["line"]),
                        s["proc"],
                        s["name"],
                        str(s["write"]),
                        str(s["readline"]),
                        s["mask"],
                        s["verdict"] or "UNCLASSIFIED",
                    )
                )
            )
    else:
        by_file: dict[str, int] = {}
        for s in shown:
            by_file[s["file"]] = by_file.get(s["file"], 0) + 1
            print(f"{s['file']}:{s['line']}  {s['proc']}  {s['name']}")
            print(f"    written {s['write']} [mask: {s['mask']}], "
                  f"read {s['readline']}")
            if s["verdict"]:
                print(f"    verdict: {s['verdict']}")
        print()
        for f, n in sorted(by_file.items(), key=lambda kv: -kv[1]):
            print(f"{n:5d}  {f}")
        print(f"{len(shown):5d}  {'POPULATION' if args.all else 'UNCLASSIFIED'}")
        if deferred and not args.all:
            by_deferred: dict[str, int] = {}
            for s in deferred:
                by_deferred[s["file"]] = by_deferred.get(s["file"], 0) + 1
            print()
            for f, n in sorted(by_deferred.items()):
                print(f"{n:5d}  {f}  DEFERRED, not read: {DEFERRED_FILES[f]}")
    if args.all:
        return 0
    return 1 if shown else 0


if __name__ == "__main__":
    raise SystemExit(main())
