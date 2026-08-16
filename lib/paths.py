"""Path display helpers shared across components.

One function, and it exists because of a specific recurring bug.

`Path.relative_to` RAISES when the path is not under the given root. Every script
here prints "wrote <path>" relative to the project root, and that print happens
*after* the artifact has been written. So passing `--output` to somewhere outside
the tree -- /tmp, a scratch directory, an absolute path typed by hand -- writes
the file correctly and then dies with a ValueError, leaving an artifact on disk
and a traceback on the terminal that looks like the run failed.

`pedology/scripts/build_soil.py` hit this and left a soil map with no provenance
record. It was then fixed at one of its two call sites and not the other, which
is the more instructive half of the story: the fix has to be applied everywhere
the pattern appears, not just where it was first observed. See
`notes/failure-modes.md`.

Provenance records have the same problem for a different reason: a path recorded
as absolute is not portable, and a path that raises is worse than either.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def rel(path: Path | str, root: Path | None = None) -> str:
    """Path relative to the project root, or absolute if it lies outside.

    Never raises. Use for anything printed to a human or written into a
    provenance record.
    """
    p = Path(path)
    base = Path(root) if root is not None else PROJECT_ROOT
    try:
        return str(p.relative_to(base))
    except ValueError:
        return str(p)
