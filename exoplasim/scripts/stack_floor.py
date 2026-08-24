#!/usr/bin/env python3
"""The per-thread stack a rung needs, as a FLOOR read off the model source.

    python exoplasim/scripts/stack_floor.py                 # the whole ladder
    python exoplasim/scripts/stack_floor.py --rung T170 --threads 16 --detail

Worldbuilding frame: an instrument on the Vesper climate model's build. Nothing
here is about the simulated planet.

WHAT IT COMPUTES AND WHAT IT DOES NOT. `plasim.f90` opens ONE parallel region
around the whole run, so every thread's stack has to hold the deepest call
chain of the entire model, not of one phase. `-fopenmp` implies `-frecursive`,
which moves the local arrays gfortran would otherwise place in static storage
onto that stack. This walks the model's own call graph from the region's entry
points and reports the heaviest chain: the sum, over the routines on it, of the
declared local arrays whose dimensions are constant at that rung.

That is a FLOOR and is labelled one. It counts declared locals only. It does
not count the temporaries gfortran materialises for whole-array and `where`
expressions -- the physics is written in array notation, so those exist and are
not small -- nor saved registers, alignment, or library frames. Nothing in the
source fixes their size; only `-fstack-usage` on a real compile does, which is
what `world-p4b` asks for and what this does not replace.

WHY IT IS TRUSTWORTHY AS A FLOOR. The same parse, summed over every routine
instead of over one chain, must reproduce what a compiled binary actually holds
in static storage, because those are the same declarations. `--validate` does
exactly that against a built executable's BSS symbols; the two agreements on
record are in `exoplasim/notes/thread-stack-floor.md`.
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

from _paths import MODEL_SRC  # noqa: F401  (also puts lib/ on sys.path)

import rungs

SRC = MODEL_SRC / "plasim" / "src"

# The files CMakeLists compiles, in its order. resmod is generated and carries
# only the three parameters this script is given as arguments; the C file and
# the deleted parallel variants are not Fortran this walks.
SOURCES = """shtnsmod mpimod_omp fftmod guimod_stub rainmod simba p_earth
carbonmod hurricanemod utilities_omp plasim plasimmod specblock calmod gaussmod
legmod outmod miscmod fluxmod radmod surfmod landmod glaciermod seamod icemod
oceanmod restartmod tracermod tpcore trc_routines aeromod aerocore lsgmod
cpl_stub""".split()

# `plasim.f90`'s parallel region, in the order it calls them.
ENTRY_POINTS = ["mpstart", "setfilenames", "opendiag", "mrdimensions",
                "allocate_arrays", "prolog", "master", "epilog", "mpstop"]

DECL = re.compile(
    r"^\s*(real|integer|logical|complex|character|double\s+precision)"
    r"\s*(\([^)]*\))?\s*((?:,\s*[a-z]+\s*(?:\([^)]*\))?)*)\s*(::)?\s*(.*)$", re.I)
ROUTINE = re.compile(
    r"^\s*(?:recursive\s+)?(subroutine|(?:real|integer|logical|"
    r"double\s+precision)?\s*function)\s+([a-z_]\w*)", re.I)
CALL = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?call\s+([a-z_]\w*)", re.I)
# On the stack only if it is neither heap nor static nor the caller's.
NOT_A_FRAME = ("allocatable", "pointer", "save", "parameter", "intent",
               "external")


def parameters(nlat: int, nlev: int, npro: int) -> dict[str, int]:
    """`plasimmod.f90`'s derived parameter block, evaluated here.

    Restated rather than parsed because these are the definitions every array
    bound in the model is written against, and a wrong one here is silent.
    """
    p = {"NLAT": nlat, "NLEV": nlev, "NPRO": npro, "NTRACE": 1, "NAERO": 1,
         "NROOT": 0, "NLEV_OCE": 1, "NPRHOR": 1}
    p["NLON"] = 2 * nlat
    p["NTRU"] = (p["NLON"] - 1) // 3
    p["NLPP"] = nlat // npro
    p["NLHP"] = p["NLPP"] // 2
    p["NHOR"] = p["NLON"] * p["NLPP"]
    p["NUGP"] = p["NLON"] * nlat
    p["NPGP"] = p["NUGP"] // 2
    p["NLEM"], p["NLEP"], p["NLSQ"] = nlev - 1, nlev + 1, nlev * nlev
    p["NTP1"] = p["NTRU"] + 1
    p["NRSP"] = (p["NTRU"] + 1) * (p["NTRU"] + 2)
    p["NCSP"] = p["NRSP"] // 2
    p["NSPP"] = (p["NRSP"] + npro - 1) // npro
    p["NESP"] = p["NSPP"] * npro
    p["NVCT"] = 2 * (nlev + 1)
    p["NZOM"] = 2 * p["NTP1"]
    return p


def _value(expr: str, p: dict[str, int]) -> int | None:
    text = expr.strip()
    if not text:
        return None
    text = re.sub(r"\b([A-Za-z_]\w*)\b",
                  lambda m: str(p[m.group(1).upper()])
                  if m.group(1).upper() in p else m.group(1), text)
    try:
        out = eval(text, {"__builtins__": {}}, {})  # noqa: S307 - arithmetic only
    except Exception:
        return None
    return int(out) if isinstance(out, (int, float)) else None


def _top_level_split(text: str) -> list[str]:
    depth, current, parts = 0, "", []
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += ch
    parts.append(current)
    return parts


def _array_bytes(dims: str, element: int, p: dict[str, int]) -> int | None:
    count = 1
    for part in _top_level_split(dims):
        part = part.strip()
        if ":" in part:
            lo, _, hi = part.partition(":")
            lo, hi = _value(lo, p), _value(hi, p)
            if lo is None or hi is None:
                return None            # assumed shape: the caller owns it
            count *= hi - lo + 1
        else:
            n = _value(part, p)
            if n is None:
                return None
            count *= n
    return count * element


def parse(precision: int, p: dict[str, int]):
    """{routine: bytes of declared local arrays}, {routine: {called routines}}."""
    element = {"real": precision, "integer": 4, "logical": 4,
               "complex": 2 * precision, "character": 1, "double": 8}
    frames: dict[str, int] = defaultdict(int)
    detail: dict[str, list] = defaultdict(list)
    calls: dict[str, set] = defaultdict(set)
    defined: set[str] = set()

    for stem in SOURCES:
        path = SRC / f"{stem}.f90"
        if not path.is_file():
            raise SystemExit(f"{path} is not there; the source list is stale")
        routine, dummies = None, set()
        for raw in path.read_text(errors="replace").splitlines():
            line = raw.split("!")[0].rstrip()
            lowered = line.lower().strip()
            start = ROUTINE.match(line)
            if start and not lowered.startswith("end"):
                routine = start.group(2).lower()
                defined.add(routine)
                head = line[line.find("("):] if "(" in line else ""
                dummies = {w.lower() for w in re.findall(r"[a-z_]\w*", head)}
                continue
            if re.match(r"^\s*end\s+(subroutine|function)\b", lowered) \
               or lowered == "end":
                routine = None
                continue
            if routine is None:
                continue
            called = CALL.match(line)
            if called:
                calls[routine].add(called.group(1).lower())
                continue
            decl = DECL.match(line)
            if not decl:
                continue
            typename = decl.group(1).lower().split()[0]
            kind, attributes, sep, rest = (decl.group(2), (decl.group(3) or "").lower(),
                                           decl.group(4), decl.group(5))
            if not sep and not re.match(r"^\s*[a-z_]\w*\s*\(", rest, re.I):
                continue               # an old-style scalar, or not a declaration
            if any(a in attributes for a in NOT_A_FRAME):
                continue
            width = element[typename]
            if kind:
                explicit = _value(kind.strip("()").replace("kind=", "")
                                  .replace("len=", ""), p)
                if explicit:
                    width = explicit
            for one in _top_level_split(rest):
                one = one.split("=")[0].strip()
                named = re.match(r"^([a-z_]\w*)\s*\((.*)\)\s*$", one, re.I | re.S)
                if not named or named.group(1).lower() in dummies:
                    continue
                size = _array_bytes(named.group(2), width, p)
                if size:
                    frames[routine] += size
                    detail[routine].append((named.group(1).lower(), size))
    return frames, detail, calls, defined


def heaviest_chain(frames, calls, defined):
    """The deepest CONCURRENT set: the heaviest path out of the parallel region."""
    memo, active, recursive = {}, set(), set()

    def walk(routine):
        if routine in memo:
            return memo[routine]
        if routine in active:
            recursive.add(routine)
            return 0, []
        active.add(routine)
        best = (0, [])
        for callee in sorted(calls.get(routine, ())):
            if callee in defined:
                deep = walk(callee)
                if deep[0] > best[0]:
                    best = deep
        active.discard(routine)
        memo[routine] = (frames.get(routine, 0) + best[0], [routine] + best[1])
        return memo[routine]

    reachable = [walk(e) for e in ENTRY_POINTS if e in defined]
    return max(reachable, key=lambda t: t[0]), sorted(recursive)


def validate(binary: Path, frames) -> None:
    """The same declarations, summed a different way, against a real binary.

    gfortran without -frecursive puts these locals in static storage, one BSS
    symbol each named `<local>.<n>`. Their total is the total this parse sees,
    so the two must agree -- a check with a right answer rather than a
    comparison that can only differ.
    """
    import subprocess
    out = subprocess.run(["nm", "-S", str(binary)], check=True, text=True,
                         capture_output=True).stdout
    measured, symbols = 0, 0
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2] in ("b", "B") \
           and re.search(r"\.\d+$", parts[3]) and "_MOD_" not in parts[3]:
            measured += int(parts[1], 16)
            symbols += 1
    parsed = sum(frames.values())
    print(f"  binary  {measured / 1e6:9.1f} MB over {symbols} BSS symbols")
    print(f"  parsed  {parsed / 1e6:9.1f} MB over {len(frames)} routines")
    if measured:
        print(f"  agree to {100 * abs(parsed - measured) / measured:.1f}%")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rung", default=None, help="one rung; default is the ladder")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--levels", type=int, default=10)
    ap.add_argument("--precision", type=int, default=8, help="bytes per real")
    ap.add_argument("--detail", action="store_true",
                    help="name the routines on the heaviest chain")
    ap.add_argument("--validate", type=Path, default=None,
                    help="a built executable to check the parse against")
    args = ap.parse_args()

    ladder = [args.rung] if args.rung else list(rungs.RUNGS)
    print(f"{args.levels} levels, {args.threads} threads, "
          f"{args.precision}-byte reals")
    print(f"{'rung':>6}  {'all locals':>12}  {'heaviest chain':>15}")
    for rung in ladder:
        nlat, _, _ = rungs.geometry(rung)
        p = parameters(nlat, args.levels, args.threads)
        frames, detail, calls, defined = parse(args.precision, p)
        (total, chain), recursive = heaviest_chain(frames, calls, defined)
        print(f"{rung:>6}  {sum(frames.values()) / 1e6:9.1f} MB  "
              f"{total / 1e6:12.1f} MB")
        if args.detail:
            for routine in chain:
                print(f"          {frames.get(routine, 0) / 1e6:9.3f} MB  {routine}")
                for name, size in sorted(detail.get(routine, []),
                                         key=lambda kv: -kv[1])[:5]:
                    print(f"            {size / 1e6:9.3f} MB  {name}")
            if recursive:
                print(f"          recursion broken at: {', '.join(recursive)}")
        if args.validate:
            validate(args.validate, frames)
    print("\nThese are FLOORS: declared local arrays only. Compiler temporaries "
          "for\nwhole-array and `where` expressions are on the same stack and "
          "are not\nvisible here. See exoplasim/notes/thread-stack-floor.md.")


if __name__ == "__main__":
    main()
