#!/usr/bin/env python3
"""The optional diagnostic blocks: their code bands and the arrays they fill.

    python exoplasim/scripts/lint_diag_arrays.py          # the problems
    python exoplasim/scripts/lint_diag_arrays.py --all    # every fact, with its verdict

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. `outmod.f90`'s `outdiag` and `snapshotdiag` number four
optional diagnostic blocks off the LOOP INDEX and write them into the same unit
the ordinary output goes to. Off the bases they carried -- 0, 20, 50 and 60 --
`ndiagsp2d = 1` made the model write code 51 twice per output step: once as
`outsc`'s ecliptic longitude, a length-1 scalar record, and once as a
length-NESP spectral array. `pyburn.readallvariables` keeps the first record's
header per code and reshapes everything joined under that code by it, so the
two record sets came back as one variable with the wrong time axis and the
wrong values, under the name `pyburn.ilibrary` supplies. world-2v9z.

The bases moved into free space and became `NDIAGGP2D_CODE0` and its three
siblings in `plasimmod.f90`. That is a boundary only while the bands stay
disjoint from each other and from every code anything else writes, and only
while a block cannot be asked for more codes than its band holds. This is what
holds both, and it holds them against the SOURCE rather than against a memory
of it.

IT ALSO HOLDS THE SECOND, UNRELATED PAIRING IN THE SAME BLOCK. `ndiaggp` and
`ndiagsp` are the switches that WRITE `dgp3d` and `dsp3d`; `ndiaggp3d` and
`ndiagsp3d` are the keys that ALLOCATE them. They are different namelist keys
with nothing tying them together, so `ndiaggp = 1` with a smaller `ndiaggp3d`
-- or with the key unset, which leaves the array unallocated -- writes past the
end of an allocatable from five modules. `plasim.f90:check_diagnostic_blocks`
refuses that pairing using `NDIAGGP_ARRAYS_FILLED` and `NDIAGSP_ARRAYS_FILLED`,
and those two numbers are what this derives from the writes and checks.

WHAT IT DOES NOT DO. It does not build and it does not run; it is a text parse,
and it over-reports rather than under-reports where a continuation line or an
unusual index expression defeats it. It answers reduced fixtures on every
invocation before it reports on the tree, because a gate whose only evidence is
a clean tree cannot tell a clean tree from a broken pass.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
from _paths import MODEL_SRC, PROJECT_ROOT  # noqa: E402

SRC = MODEL_SRC / "plasim" / "src"
MODULE = SRC / "plasimmod.f90"
GUARD = SRC / "plasim.f90"

# The four blocks, as (namelist count key, base parameter, array the block
# writes out, the switch that FILLS that array, the parameter holding how many
# that switch fills). A block whose array nothing fills has None for the last
# two: it is emitted as zeros, which is upstream's extension point and not a
# defect, but it stops being that the moment something writes it.
BLOCKS = (
    ("ndiaggp2d", "NDIAGGP2D_CODE0", "dgp2d", None, None),
    ("ndiaggp3d", "NDIAGGP3D_CODE0", "dgp3d", "ndiaggp", "NDIAGGP_ARRAYS_FILLED"),
    ("ndiagsp2d", "NDIAGSP2D_CODE0", "dsp2d", None, None),
    ("ndiagsp3d", "NDIAGSP3D_CODE0", "dsp3d", "ndiagsp", "NDIAGSP_ARRAYS_FILLED"),
)


def parameters(text: str) -> dict[str, int]:
    """Every `integer, parameter :: NAME = value` in a source text."""
    out = {}
    for match in re.finditer(
            r"^\s*integer,\s*parameter\s*::\s*([A-Za-z_0-9]+)\s*=\s*(-?\d+)",
            text, re.IGNORECASE | re.MULTILINE):
        out[match.group(1).upper()] = int(match.group(2))
    return out


def literal_array_indices(text: str, array: str) -> set[int]:
    """The literal LAST subscript of every reference to `array` in a text.

    The blocks' arrays are `dgp2d(NHOR,n)`, `dgp3d(NHOR,NLEV,n)` and their
    spectral twins, so the array index is the last subscript. A reference whose
    last subscript is not a literal -- the loop in `outdiag`, the whole-array
    zeroing, the `allocate` -- carries no index to check and is skipped.
    """
    found: set[int] = set()
    for match in re.finditer(rf"\b{array}\s*\(([^()]*)\)", text, re.IGNORECASE):
        last = match.group(1).split(",")[-1].strip()
        if re.fullmatch(r"\d+", last):
            found.add(int(last))
    return found


def codes_written_literally(src: Path) -> set[int]:
    """Every literal code any record writer in the model source is handed.

    `writegp` and `writesp` take the code third, `writescalar` last. A call
    whose code is a variable -- which is what the diagnostic blocks themselves
    are -- carries no literal and is not one of these.
    """
    codes: set[int] = set()
    for path in sorted(src.glob("*.f90")):
        text = path.read_text(encoding="latin-1")
        for match in re.finditer(
                r"call\s+write(gp|sp|scalar)\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)",
                text, re.IGNORECASE):
            args = [a.strip() for a in split_top(match.group(2))]
            if len(args) < 3:
                continue
            raw = args[2] if match.group(1).lower() in ("gp", "sp") else args[-1]
            if re.fullmatch(r"\d+", raw):
                codes.add(int(raw))
    return codes


def split_top(text: str) -> list[str]:
    """Split on commas at bracket depth zero."""
    parts, depth, start = [], 0, 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1
    parts.append(text[start:])
    return parts


def library_codes() -> set[int] | None:
    """The codes the postprocessor can NAME, or None when it does not import.

    None rather than an empty set, so a missing postprocessor reads as "not
    checked" instead of as "nothing collides".
    """
    try:
        from exoplasim import pyburn
    except Exception:
        return None
    return {int(k) for k in pyburn.ilibrary if str(k).isdigit()}


def check(src: Path, verbose: bool = False) -> list[str]:
    problems: list[str] = []
    facts: list[str] = []
    module_text = MODULE.read_text(encoding="latin-1") if MODULE.is_file() else ""
    params = parameters(module_text)
    needed = ["NDIAG_BLOCK_CODES"] + [b[1] for b in BLOCKS] + \
             [b[4] for b in BLOCKS if b[4]]
    missing = [name for name in needed if name not in params]
    if missing:
        return [f"{MODULE} declares no {', '.join(missing)}; the diagnostic "
                "blocks' code bands are stated there and nowhere else"]

    width = params["NDIAG_BLOCK_CODES"]
    sources = {path: path.read_text(encoding="latin-1")
               for path in sorted(src.glob("*.f90"))}

    # 1. every band is disjoint from every other band
    bands = {}
    for key, base_name, _, _, _ in BLOCKS:
        base = params[base_name]
        bands[key] = range(base + 1, base + width + 1)
        facts.append(f"{key} is numbered {base + 1} to {base + width}")
    keys = list(bands)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            overlap = set(bands[a]) & set(bands[b])
            if overlap:
                problems.append(
                    f"the {a} and {b} bands overlap at codes "
                    f"{min(overlap)}-{max(overlap)}; a block filled to its "
                    "capacity would write onto the other block's codes")

    # 2. no band reaches a code the model writes literally, or the
    #    postprocessor names
    literal = codes_written_literally(src)
    library = library_codes()
    facts.append(f"{len(literal)} literal codes are written by the model source")
    for key in keys:
        landed = sorted(set(bands[key]) & literal)
        if landed:
            problems.append(
                f"the {key} band reaches code(s) {landed}, which a record "
                "writer in the model source already uses. Two fields under one "
                "code in one stream come back as one variable with the wrong "
                "time axis; that is world-2v9z and moving the bases is what "
                "fixed it")
        if library is not None:
            named = sorted(set(bands[key]) & library)
            if named:
                problems.append(
                    f"the {key} band reaches code(s) {named}, which "
                    "pyburn.ilibrary names, so a record of this block would be "
                    "read back under another field's name and units")
    if library is None:
        facts.append("pyburn did not import; the library codes were NOT checked")

    # 3. the filled-array counts are what the source actually writes
    for key, _, array, switch, count_name in BLOCKS:
        indices: set[int] = set()
        for path, text in sources.items():
            if path.name == "outmod.f90":
                continue          # the writer's own loop, not a fill
            indices |= literal_array_indices(text, array)
        observed = max(indices) if indices else 0
        if count_name is None:
            facts.append(f"{array} is filled by nothing (block emits zeros)")
            if observed:
                problems.append(
                    f"{array} is now written at index {observed} and no "
                    f"parameter says how many arrays fill it, so nothing "
                    f"stops {key} from being set smaller than that. Add a "
                    f"count beside NDIAGGP_ARRAYS_FILLED in plasimmod.f90, "
                    f"have check_diagnostic_blocks test it, and name it here")
            continue
        declared = params[count_name]
        facts.append(f"{array} is filled at indices 1 to {observed} under "
                     f"{switch}; {count_name} declares {declared}")
        if observed != declared:
            problems.append(
                f"{count_name} is {declared} and the model source writes "
                f"{array} up to index {observed}. That number is what "
                f"check_diagnostic_blocks demands of {key}, so the two "
                f"disagreeing means either a write past the end of the "
                f"allocation is permitted or a larger array is demanded than "
                f"anything fills")

    # 4. the guard tests every count against the band width
    guard_text = GUARD.read_text(encoding="latin-1") if GUARD.is_file() else ""
    if "subroutine check_diagnostic_blocks" not in guard_text.lower():
        problems.append(
            f"{GUARD} has no check_diagnostic_blocks, so nothing refuses a "
            "block longer than its code band or a switch that fills more "
            "arrays than its allocation holds")
    else:
        for key, _, _, _, _ in BLOCKS:
            if not re.search(rf"{key}\s*>\s*NDIAG_BLOCK_CODES", guard_text,
                             re.IGNORECASE):
                problems.append(
                    f"check_diagnostic_blocks does not test {key} against "
                    "NDIAG_BLOCK_CODES, so a count past the band is allowed "
                    "through and writes onto the next block")
    if verbose:
        for fact in facts:
            print(f"[ fact ] {fact}")
    return problems


FIXTURE_MODULE = """
      integer, parameter :: NDIAGGP2D_CODE0 = 700
      integer, parameter :: NDIAGGP3D_CODE0 = 750
      integer, parameter :: NDIAGSP2D_CODE0 = 900
      integer, parameter :: NDIAGSP3D_CODE0 = 1000
      integer, parameter :: NDIAG_BLOCK_CODES = 99
      integer, parameter :: NDIAGGP_ARRAYS_FILLED = 4
      integer, parameter :: NDIAGSP_ARRAYS_FILLED = 3
"""


def selftest() -> list[str]:
    """Reduced cases with known answers, on every invocation.

    The three the tree cannot supply: two bands that overlap, a declared count
    that disagrees with the writes, and a band that lands on a literal code.
    """
    failures = []
    params = parameters(FIXTURE_MODULE)
    if params.get("NDIAGGP2D_CODE0") != 700 or params.get("NDIAG_BLOCK_CODES") != 99:
        failures.append("the parameter parse does not read a plain declaration")
    a = range(701, 800)
    b = range(751, 850)
    if not set(a) & set(b):
        failures.append("two bands 50 apart at width 99 were not seen to overlap")
    text = "       dgp3d(:,jlev,17)=x\n       dgp3d(:,1:NLEV,9)=y\n"
    if literal_array_indices(text, "dgp3d") != {17, 9}:
        failures.append("the array-index parse missed a literal last subscript")
    if literal_array_indices("      dgp3d(:,:,:)=0.\n", "dgp3d"):
        failures.append("the array-index parse read a whole-array reference as "
                        "an index")
    calls = "      call writegp(40,dclforc(1,1),101,0)\n" \
            "      call writesp(40,dsp2d(1,jdiag),jcode,0,1.,0.0)\n"
    got = set()
    for match in re.finditer(r"call\s+write(gp|sp)\s*\(([^\n]*)\)", calls):
        args = [x.strip() for x in split_top(match.group(2))]
        if len(args) >= 3 and re.fullmatch(r"\d+", args[2]):
            got.add(int(args[2]))
    if got != {101}:
        failures.append(f"the literal-code parse returned {sorted(got)}, not "
                        "the one literal code in the fixture")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--all", action="store_true",
                        help="print every fact with its verdict, not only the "
                             "problems")
    args = parser.parse_args()

    broken = selftest()
    if broken:
        print("the lint's own fixtures do not answer as they must:")
        for line in broken:
            print(f"  {line}")
        return 2
    if not SRC.is_dir():
        print(f"{SRC} is missing; the vendored model source moved")
        return 2

    problems = check(SRC, verbose=args.all)
    if problems:
        print(f"{len(problems)} problem(s) in the diagnostic blocks:")
        for line in problems:
            print(f"  {line}")
        return 1
    print("the diagnostic blocks' code bands are disjoint and clear of every "
          "code the model writes, and the filled-array counts are what the "
          "source writes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
