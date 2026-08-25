#!/usr/bin/env python3
"""Procedure-body declarations that acquire SAVE, which under threads is sharing.

    python exoplasim/scripts/lint_implicit_save.py          # the unclassified
    python exoplasim/scripts/lint_implicit_save.py --all    # the whole population
    python exoplasim/scripts/lint_implicit_save.py --tsv    # machine-readable

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. A local declared with an initialiser inside a procedure
acquires the SAVE attribute and persists across calls. Under `-fopenmp` a SAVEd
local is SHARED by the whole thread team, and `plasim.f90` opens one parallel
region around `mpstart` through `mpstop`, so every procedure in the model's call
tree runs on NPRO threads. `notes/audits/uninitialised-reads-and-implicit-save.md`
measured the race: eight threads writing their own id into a SAVEd scratch array
and reading it back, seven of the eight reading a value another thread wrote.

Under the MPI arm each rank was a process and its SAVEd locals were private by
construction, so the same source was correct. world-38b removed that arm. There
is no longer a configuration in which this class is harmless, and gfortran
reports nothing about it at any warning level tried -- so there is no flag that
shortens this and it is reading.

THE PASS OVER-REPORTS, DELIBERATELY. Four places it errs towards reporting a
declaration that is fine, and none towards missing one:

  1. A declaration it cannot attribute to an enclosing procedure is reported,
     not dropped.
  2. An explicit `save`, and a `data` statement, are reported alongside the
     implicit case. They are the identical sharing written out.
  3. A variable named on an `!$omp threadprivate` directive counts as answered
     only when the directive is in the SAME procedure. A directive in a module
     specification part privatises a module variable, a different construct,
     and does not clear a local of the same name.
  4. Continuation lines are joined, so a declaration split across lines is one
     site and cannot hide half of itself from the initialiser test.

The one shape it does NOT look for is a declaration without `::` that carries an
initialiser, because Fortran has no such shape: initialisation in a type
declaration requires the double colon. The old-style `character(len=18) datch`
declarations in `plasim.f90` are therefore not of this class, and testing them
for `=` finds only the `len=` inside the type specifier.

WHAT IT DOES NOT KNOW. Whether a site is a defect is not decidable from the
declaration. A local that is written before it is read on every call never
exposes its storage class whatever the thread count; a counter meant to persist
is a defect the moment two threads reach it. Both look identical here. So the
count is a POPULATION and the verdict is in the CLASSIFIED table below, one row
per site with the reason it is safe.
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
_DECL_RE = re.compile(r"^\s*(" + "|".join(_TYPES) + r")", re.I)
_PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:real|integer|logical|complex|double\s+precision|character)"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+([a-z_]\w*)",
    re.I,
)
_PROC_END_RE = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
_MODULE_RE = re.compile(r"^\s*module\s+([a-z_]\w*)\s*$", re.I)
_PROGRAM_RE = re.compile(r"^\s*program\s+([a-z_]\w*)", re.I)
_SAVE_STMT_RE = re.compile(r"^\s*save\b", re.I)
_DATA_STMT_RE = re.compile(r"^\s*data\s+", re.I)
_THREADPRIVATE_RE = re.compile(r"^\s*!\$omp\s*&?\s*threadprivate\s*\(", re.I)


# ---------------------------------------------------------------------------
# THE CLASSIFIED POPULATION, one row per site with the reason it is safe.
#
# clim-51 read every site the pass reports and recorded the verdict here rather
# than in a document, so that a declaration which CHANGES drops out of the table
# and is reported again. The key is (file, procedure, variable); the value is
# the argument. A row is a claim, and a row whose reason does not hold is a
# defect in this table.
#
# `--all` prints the whole population with its verdicts; the default run prints
# only what is NOT here, and exits 1 if anything is. Evidence:
# notes/audits/implicit-save-under-threads.md.
# ---------------------------------------------------------------------------

CLASSIFIED: dict[tuple[str, str, str], str] = {
    ("buildice.f90", "buildice", "dsnowz"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("buildice.f90", "buildice", "lsm"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnow1"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnow2"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnow3"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnow4"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnow5"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "dsnowz"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "xsnowz"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "lsm"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "persist"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "deltyrs"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("newsnow.f90", "newsnow", "delsnow"):
        "a main program's variables have static storage whatever their "
        "declaration, and a standalone offline program is one thread",
    ("calmod.f90", "ntodat", "mona"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("fft991mod.f90", "set99", "LFAX"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "obamp"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "obrate"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "obphas"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "ecamp"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "ecrate"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "ecphas"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "mvamp"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "mvrate"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("radmod.f90", "orb_params", "mvphas"):
        "read-only after its data statement, so one shared copy is what "
        "every thread wants",
    ("trc_routines.f90", "filns", "first"):
        "root-only by caller: tracer_main calls tpcore under if (mypid == "
        "NROOT), so one thread reaches this first-call flag. The ONLY site "
        "in the population whose safety is a property of the caller rather "
        "than of the code in front of it",
    ("trc_routines.f90", "filns", "cap1"):
        "root-only by caller, as `first` above; the cached value is a "
        "function of the grid, which is fixed for a run",
}


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


def logical_lines(path: Path):
    raw = path.read_text(errors="replace").splitlines()
    buf = ""
    start = None
    for n, line in enumerate(raw, 1):
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


def split_top(text: str, sep: str = ",") -> list[str]:
    """Split on `sep` at parenthesis depth zero."""
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


def statement_names(text: str, keyword: str) -> list[str]:
    """The variable names a `data` or `save` statement gives static storage to.

    `data` is `data name-list /values/ [, name-list /values/]...`, so the names
    are what precedes the first slash. `save` is a bare name list. Both are
    keyed on the NAME rather than on the statement text, so a coefficient edit
    does not silently change the key and drop the site out of the table.
    """
    body = text.strip()[len(keyword) :]
    if keyword == "data":
        body = body.split("/", 1)[0]
    names = []
    for part in split_top(body):
        name = part.split("(")[0].strip()
        if name and re.match(r"^[a-z_]\w*$", name, re.I):
            names.append(name)
    return names or [squeeze(text)[:60]]


def initialised_names(rhs: str) -> list[str]:
    """The names in an entity list that carry an `=` initialiser."""
    names = []
    for entity in split_top(rhs):
        if "=" not in entity:
            continue
        name = entity.split("=")[0]
        name = name.split("(")[0].split("*")[0].strip()
        if name:
            names.append(name)
    return names


def threadprivate_by_proc(path: Path) -> dict[str, set[str]]:
    """Names an `!$omp threadprivate` directive privatises, per procedure.

    A SAVEd local named on such a directive inside the SAME procedure is one
    copy per thread and not one copy per team, which is the whole of the
    hazard. A directive in a module specification part privatises a module
    variable, a different construct, so only same-procedure directives count.
    """
    out: dict[str, set[str]] = {}
    proc = None
    buf = ""
    for line in path.read_text(errors="replace").splitlines():
        low = line.lower()
        if re.match(r"^\s*end\s*(subroutine|function)\b", low):
            proc = None
        pm = _PROC_START_RE.match(strip_comment(line).lower())
        if pm:
            proc = pm.group(2)
        if _THREADPRIVATE_RE.match(low) or (buf and low.lstrip().startswith("!$omp")):
            buf += low
            if buf.rstrip().endswith("&"):
                continue
            names = re.findall(r"[a-z_]\w*", buf.split("(", 1)[1])
            if proc:
                out.setdefault(proc, set()).update(names)
            buf = ""
    return out


def scan(path: Path):
    proc = None
    in_module = False
    sites = []
    private = threadprivate_by_proc(path)
    for lineno, text in logical_lines(path):
        low = text.lower()

        if _MODULE_RE.match(low) and "procedure" not in low:
            in_module = True
            continue
        if re.match(r"^\s*end\s*module\b", low):
            in_module = False
            continue
        if _PROC_END_RE.match(low):
            proc = None
            continue
        pm = _PROC_START_RE.match(low)
        if pm:
            proc = pm.group(2)
            continue
        gm = _PROGRAM_RE.match(low)
        if gm:
            proc = gm.group(1)
            continue

        # Module-level and file-level static storage is intended; only a
        # PROCEDURE BODY acquires SAVE by surprise.
        if proc is None:
            continue

        kind = None
        names: list[str] = []
        if _SAVE_STMT_RE.match(low):
            kind = "save statement"
            names = statement_names(text, "save")
        elif _DATA_STMT_RE.match(low):
            kind = "data statement"
            names = statement_names(text, "data")
        elif _DECL_RE.match(low) and "::" in text:
            head, _, rhs = text.partition("::")
            attrs = head.lower()
            if "parameter" in attrs:
                continue
            got = initialised_names(rhs)
            if "save" in attrs:
                kind = "save attribute"
                names = got or [
                    n.split("(")[0].split("*")[0].strip()
                    for n in split_top(rhs)
                    if n.strip()
                ]
            elif got:
                kind = "initialiser"
                names = got

        if not kind:
            continue
        for name in names:
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind=kind,
                    name=name,
                    text=squeeze(text)[:110],
                    module=in_module,
                    private=name.lower() in private.get(proc, ()),
                )
            )
    return sites


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--all", action="store_true", help="print the whole population with verdicts"
    )
    ap.add_argument("--tsv", action="store_true", help="one tab-separated row per site")
    args = ap.parse_args()

    sites = []
    for path in sorted(SRC.glob("*.f90")):
        sites.extend(scan(path))
    for s in sites:
        s["verdict"] = CLASSIFIED.get((s["file"], s["proc"], s["name"]))

    shown = (
        sites
        if args.all
        else [s for s in sites if not s["verdict"] and not s["private"]]
    )

    if args.tsv:
        for s in shown:
            print(
                "\t".join(
                    (
                        s["file"],
                        str(s["line"]),
                        s["proc"],
                        s["kind"],
                        s["name"],
                        "threadprivate"
                        if s["private"]
                        else (s["verdict"] or "UNCLASSIFIED"),
                    )
                )
            )
    else:
        by_file: dict[str, int] = {}
        for s in shown:
            by_file[s["file"]] = by_file.get(s["file"], 0) + 1
            print(f"{s['file']}:{s['line']}  {s['proc']}  {s['name']}  ({s['kind']})")
            print(f"    {s['text']}")
            if s["private"]:
                print("    verdict: threadprivate in the same procedure")
            elif s["verdict"]:
                print(f"    verdict: {s['verdict']}")
        print()
        for f, n in sorted(by_file.items(), key=lambda kv: -kv[1]):
            print(f"{n:5d}  {f}")
        print(f"{len(shown):5d}  {'POPULATION' if args.all else 'UNCLASSIFIED'}")
    if args.all:
        return 0
    return 1 if shown else 0


if __name__ == "__main__":
    raise SystemExit(main())
