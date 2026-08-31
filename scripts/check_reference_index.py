#!/usr/bin/env python3
"""Every filename references/INDEX.md asserts is a file that is there.

INDEX.md's first column is either a filename in backticks or prose. A backticked
filename is an ASSERTION that the PDF is on disk, and rows marked **read** carry
numbers transcribed out of it. When the file goes and the row does not, the row
still reads as a citation whose artifact a reader can open, and the project's own
standard -- claims are checked against the artifact rather than the documentation
-- has quietly become unavailable for exactly the rows that most need it.

The index already carries the shape for a source that is cited but not held: a
prose first cell saying so, with what was tried. This check exists to make the
absence of a file force that rewrite rather than sit unnoticed.

The library is untracked, so a checkout without it has nothing to check and the
check reports that instead of failing. A worktree linked by scripts/link_worktree.py
sees the main checkout's library and is checked normally.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "references"
INDEX = REFERENCES / "INDEX.md"

# The first cell of a table row, when it is a single backticked token.
_ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|")


def asserted_files(text: str) -> list[tuple[int, str]]:
    """Every (line number, path) the index's first column claims is on disk.

    A first cell that is prose, `*(not held)*`, or a directory (trailing `/`) is
    not a file assertion. Only the FIRST cell counts: citation cells quote other
    filenames -- an OCR variant, a companion volume -- and those are prose about
    a file, not the row's own claim to hold one.
    """
    out = []
    for number, line in enumerate(text.splitlines(), start=1):
        match = _ROW.match(line)
        if match is None:
            continue
        cell = match.group(1).strip()
        if cell.endswith("/"):
            continue
        if not cell.lower().endswith(".pdf"):
            continue
        out.append((number, cell))
    return out


def check_reference_index() -> list[str]:
    if not INDEX.exists():
        return [f"{INDEX} does not exist"]

    rows = asserted_files(INDEX.read_text())
    if not rows:
        return [f"{INDEX} asserts no filenames at all; the row pattern has changed"]

    held = {p.name for p in REFERENCES.rglob("*.pdf")}
    if not held:
        print(f"    references/ holds no PDFs in this checkout; "
              f"{len(rows)} asserted filenames not checked")
        return []

    problems = []
    for number, cell in rows:
        path = REFERENCES / cell
        if path.exists():
            continue
        # Named at the top level but filed in a subdirectory is a path error in
        # the row, not a missing paper; say which it is.
        elsewhere = sorted(p.relative_to(REFERENCES).as_posix()
                           for p in REFERENCES.rglob(Path(cell).name))
        if elsewhere:
            problems.append(f"INDEX.md:{number} names `{cell}`, which is at "
                            f"{', '.join(elsewhere)}")
        else:
            problems.append(f"INDEX.md:{number} names `{cell}`, which is not on disk. "
                            f"Re-fetch it, or rewrite the row to the shape the index "
                            f"uses for a source it cites and does not hold")
    return problems


def main() -> None:
    argparse.ArgumentParser(
        description="Check that every filename references/INDEX.md asserts is on "
                    "disk. Takes no arguments.").parse_args()
    problems = check_reference_index()
    for problem in problems:
        print(f"  {problem}")
    if problems:
        print(f"\n{len(problems)} row(s) assert a file that is not there")
        sys.exit(1)
    print("every filename references/INDEX.md asserts is on disk")


if __name__ == "__main__":
    main()
