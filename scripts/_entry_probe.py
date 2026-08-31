#!/usr/bin/env python3
"""One entry point, started in this interpreter, with a classified verdict.

    python scripts/_entry_probe.py --mode import <script>
    python scripts/_entry_probe.py --mode help   <script>

`scripts/verify_entry_points.py` spawns one of these per entry point and reads
the single JSON line it writes to stdout. It is `_`-prefixed so the sweep does
not try to start the prober with itself.

WHY A SEPARATE PROCESS SAYS WHAT HAPPENED RATHER THAN THE PARENT GUESSING.
The parent used to decide "this script does not start" from the child's exit
status alone, which equates two unrelated events: an import graph that no longer
resolves, and a gate that started perfectly well and refused on purpose. Both
leave a non-zero status. The distinction is not recoverable from outside the
child -- scraping a traceback out of stderr guesses at it, and a per-script
allowlist would go stale the moment a script changed its mind -- but INSIDE the
child it is an exception type, so this runs the script under `except` clauses
that name it and reports the class it caught.

THE TWO MODES. `import` executes the module body with a `__name__` that is not
`__main__`, so every module-level statement runs, every import resolves or does
not, and the `if __name__ == "__main__"` block does not fire. `help` is
`python <script> --help`, which additionally builds the argparse parser. The
parent picks per script: `help` where the source constructs an
`ArgumentParser`, because there the parser short-circuits before the work, and
`import` where it does not, because there nothing short-circuits and starting
the script RUNS it.

FAITHFULNESS. `python <script>` puts the script's own directory at the head of
`sys.path`, and many scripts here rely on that to reach a sibling; the prober
therefore REPLACES its own directory with the script's rather than adding to it,
so a script that only imports its neighbour because the prober is in `scripts/`
fails here as it would in the wild. `sys.argv` is set to what the mode implies.
The working directory is the parent's, which is the repository root.

THE VERDICT CHANNEL. Started scripts write to stdout and stderr, and `--help`
writes a screenful, so fd 1 and fd 2 are pointed at a scratch file for the
duration and the verdict goes to a saved duplicate of fd 1. Redirecting at the
FILE DESCRIPTOR rather than swapping `sys.stdout` is what makes it hold for a
C extension or a grandchild process, which write to the descriptor directly.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import runpy
import sys
import tempfile

# The four verdicts. `ok` and `refusal` are outcomes of a script that started;
# `import` and `start` are outcomes of one that did not.
OK = "ok"
IMPORT = "import"      # ImportError / ModuleNotFoundError: the graph is broken
START = "start"        # any other exception: a statement or the parser raised
REFUSAL = "refusal"    # SystemExit with a non-zero code: the script chose this


def _last_line(text: str) -> str:
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    return lines[-1] if lines else "(no output)"


def probe(script: Path, mode: str) -> dict:
    """Start `script` and name what happened. Never raises."""
    sys.path[0] = str(script.parent)
    sys.argv = [str(script)] + (["--help"] if mode == "help" else [])
    # `__main__` for `help`, because that is what `python <script> --help` is;
    # anything else for `import`, because the guard must not fire.
    run_name = "__main__" if mode == "help" else f"_probe_{script.stem}"

    saved_out, saved_err = os.dup(1), os.dup(2)
    scratch = tempfile.TemporaryFile(mode="w+", encoding="utf-8",
                                     errors="replace")
    kind, detail = OK, ""
    try:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(scratch.fileno(), 1)
        os.dup2(scratch.fileno(), 2)
        try:
            runpy.run_path(str(script), run_name=run_name)
        except ImportError as e:
            kind, detail = IMPORT, f"{type(e).__name__}: {e}"
        except SystemExit as e:
            code = e.code
            if code not in (0, None):
                kind = REFUSAL
                detail = code if isinstance(code, str) else f"exit {code}"
        except BaseException as e:  # noqa: BLE001 -- classifying, not handling
            kind = START
            detail = f"{type(e).__name__}: {e}".strip().splitlines()[0]
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(saved_out, 1)
            os.dup2(saved_err, 2)
    finally:
        os.close(saved_out)
        os.close(saved_err)

    scratch.seek(0)
    written = scratch.read()
    scratch.close()
    if kind in (START, REFUSAL) and written.strip():
        detail = f"{detail}: {_last_line(written)}" if detail else _last_line(written)
    return {"kind": kind, "detail": detail, "mode": mode}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("script", type=Path)
    ap.add_argument("--mode", choices=("import", "help"), required=True)
    a = ap.parse_args()
    verdict = probe(a.script.resolve(), a.mode)
    sys.stdout.write(json.dumps(verdict) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
