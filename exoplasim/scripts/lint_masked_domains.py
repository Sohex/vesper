#!/usr/bin/env python3
"""Domain-sensitive intrinsics and divisions that sit inside a masked WHERE block.

    python exoplasim/scripts/lint_masked_domains.py            # the survivors
    python exoplasim/scripts/lint_masked_domains.py --all      # every site, guarded or not
    python exoplasim/scripts/lint_masked_domains.py --tsv      # machine-readable

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. A Fortran `where` block does NOT protect its right-hand sides.
The mask selects which elements the ASSIGNMENT stores; the compiler is free to
evaluate the expression on every element of the array, and a vectorising one
does exactly that. The build profile carries
`-ffpe-trap=invalid,zero,overflow` (`config/planet.yaml`), so a lane that was
going to be discarded raises SIGFPE instead of being discarded. That is not
hypothetical: a T42 run died eight orbits in and the fault moved to a new site
each time one was clamped (world-bhs, world-5a0).

`-O2` rather than `-O3` removed the loop transforms that were firing. It did
not remove the class, because `-O2` vectorises too. So the class needs a source
answer, and a source answer needs the population enumerated.

THE PASS OVER-REPORTS, DELIBERATELY, and that direction is the whole of its
evidential value. Three places it errs towards reporting a site that is safe,
and none towards clearing one that is not:

  1. The WHERE stack is popped only by an explicit `end where` / `endwhere`,
     or by a procedure boundary. A construct this parser fails to close keeps
     reporting sites after it, never fewer.
  2. A guard is recognised only from a small, literal set of shapes (an
     argument already wrapped in `max`/`min`/`abs`, or a divisor already
     wrapped in `max`/`sign`). An unrecognised but correct guard is reported.
  3. Division is reported as its own class even though a divisor is only a
     hazard when it can reach zero, which the pass cannot know.

So a clean run is evidence and a dirty one is a work list. A pass whose
direction is unstated is not evidence at all, which is why this paragraph is
here.

WHAT IT DOES NOT KNOW. It cannot tell whether a mask's complement actually
holds an out-of-domain value -- that needs the physics, and it is why the
output is read in context rather than acted on mechanically. It also has no
opinion on `merge(a,b,mask)`, which has the identical hazard for the same
reason, so those sites are reported under their own class.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import MODEL_SRC  # noqa: E402

SRC = MODEL_SRC / "plasim" / "src"

# The intrinsics whose argument has a restricted domain, in every spelling the
# source uses. Longest first so `log10` is not matched as `log`.
DOMAIN_INTRINSICS = (
    "alog10",
    "log10",
    "alog",
    "log",
    "sqrt",
    "exp",
    "asin",
    "acos",
    "acosh",
    "atanh",
)
_INTRINSIC_RE = re.compile(
    r"(?<![a-z0-9_])(" + "|".join(DOMAIN_INTRINSICS) + r")\s*\(", re.I
)
_MERGE_RE = re.compile(r"(?<![a-z0-9_])merge\s*\(", re.I)

# Block openers/closers. `where (mask)` with nothing after the closing paren
# opens a construct; with a statement after it, the mask covers that one
# assignment and the construct is not open.
_WHERE_OPEN_RE = re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?where\s*\(", re.I)
_END_WHERE_RE = re.compile(r"^\s*end\s*where\b", re.I)
_ELSEWHERE_RE = re.compile(r"^\s*else\s*where\b", re.I)
_PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:real|integer|logical|complex|double\s+precision|character)"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function|program|module)\s+[a-z_]\w*",
    re.I,
)
_PROC_END_RE = re.compile(r"^\s*end\s*(subroutine|function|program|module)\b", re.I)

# Guard shapes already in the tree. An argument wrapped in one of these has had
# its domain bounded by hand; a divisor wrapped in max/sign has had its pole
# floored. Recognised literally, so an unrecognised guard is REPORTED.
_ARG_GUARDS = ("max", "min", "amax1", "amin1", "abs", "dim")
_DIV_GUARDS = ("max", "amax1", "sign", "abs")


# ---------------------------------------------------------------------------
# THE CLASSIFIED SURVIVORS, and why each is not a defect.
#
# world-5a0 read every site the pass reported and made the arguments safe where
# a floor could be argued free. The sites below are the ones where NO guard was
# added, each with the argument that makes it safe. They are keyed on the file,
# the procedure, the intrinsic and the ARGUMENT TEXT rather than on a line
# number, so the table survives edits above them and a changed argument drops
# out of the table and is reported again -- which is the direction that matters.
#
# `--all` prints them; the default run prints only what is NOT here, so a new
# masked intrinsic is visible the moment it is added. Adding a row is a claim
# with a reason, and a row whose reason does not hold is a defect in this table.
#
# Evidence: notes/audits/masked-where-blocks.md.
# ---------------------------------------------------------------------------

CLASSIFIED: dict[tuple[str, str, str, str], str] = {
    ("radmod.f90", "subroutine swr", "alog10", "1.5+max(0.,zlwp(:))"):
        "argument is at least 1.5 by the floor beside it",
    ("radmod.f90", "subroutine swr", "alog", "3.+0.1*ztau(:)"):
        "ztau follows the floored alog10 above and is preset to 1, so the "
        "argument is at least 3",
    ("radmod.f90", "subroutine swr", "sqrt", "273./dt(:,jlev)"):
        "the mask is losun and the argument is the gridpoint temperature, "
        "which the mask has no bearing on; masked and kept lanes are equally "
        "safe and this is not a masked-domain site",
    ("radmod.f90", "subroutine swr", "exp", "zaertf1(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf1(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf1s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf1s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf2(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf2(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf2s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf2s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "log", "1.+zcb1*zco2(:)"):
        "the column amount is preset to 0 unconditionally and accumulates "
        "non-negative terms, so the argument is at least 1",
    ("radmod.f90", "subroutine swr", "log", "1.+zcb2*zco2(:)"):
        "the column amount is preset to 0 unconditionally and accumulates "
        "non-negative terms, so the argument is at least 1",
    ("radmod.f90", "subroutine swr", "log", "1.0-0.144"):
        "constant argument",
    ("radmod.f90", "subroutine swr", "exp", "zscf(:)*log(1.0-0.144)"):
        "constant logarithm times a non-negative scale factor",
    ("radmod.f90", "subroutine swr", "log",
     "1.0-(0.219/(1.+0.816*max(0.,zmu0(:))))"):
        "the floor holds the argument at or above 0.781 for any zmu0",
    ("radmod.f90", "subroutine swr", "exp",
     "zscf(:)*log(1.0-(0.219/(1.+0.816*max(0.,zmu0(:)))))"):
        "the logarithm inside it is bounded, and zscf is non-negative",
    ("radmod.f90", "subroutine lwr", "alog", "ztau0(:)"):
        "clamped to [zero, 1-zero] on every lane by the AMIN1/MAX two lines "
        "above the where, unconditionally",
    ("seamod.f90", "subroutine seaini", "exp",
     "ra2*(dt(:,NLEP)-TMELT) /ra4d(dt(:,NLEP),ra4)"):
        "ra4d floors the pole; the quotient tends to ra2 for large dt and to "
        "a large negative for small dt, which underflows rather than traps",
    ("seamod.f90", "subroutine seastep", "exp",
     "ra2*(dt(:,NLEP)-TMELT) /ra4d(dt(:,NLEP),ra4)"):
        "ra4d floors the pole; the quotient tends to ra2 for large dt and to "
        "a large negative for small dt, which underflows rather than traps",
    ("seamod.f90", "subroutine seastep", "sqrt", "dtaux(:)**2+dtauy(:)**2"):
        "a sum of squares cannot be negative on any lane",
}


def squeeze(text: str) -> str:
    """Collapse runs of whitespace, so a continuation's padding is not a key."""
    return " ".join(text.split())


def strip_comment(line: str) -> str:
    """Drop a trailing `!` comment, respecting quoted strings."""
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
    """Yield (first_physical_lineno, joined_text) with continuations folded."""
    raw = path.read_text(errors="replace").splitlines()
    buf = ""
    start = None
    for n, line in enumerate(raw, 1):
        code = strip_comment(line).rstrip()
        if not code.strip():
            if buf:
                continue
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
    """Index of the `)` closing the `(` at open_idx, or -1."""
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


def outer_call(text: str, name_start: int) -> str:
    """The token immediately enclosing the call that starts at name_start."""
    i = name_start - 1
    while i >= 0 and text[i] in " \t":
        i -= 1
    if i < 0 or text[i] != "(":
        return ""
    j = i - 1
    while j >= 0 and text[j] in " \t":
        j -= 1
    end = j + 1
    while j >= 0 and (text[j].isalnum() or text[j] == "_"):
        j -= 1
    return text[j + 1 : end].lower()


def arg_of(text: str, name_end: int) -> str:
    close = match_paren(text, name_end - 1)
    return text[name_end:close] if close > 0 else text[name_end:]


def leading_call(expr: str) -> str:
    """The function name a parenthesised expression opens with, if any."""
    m = re.match(r"\s*([a-z_]\w*)\s*\(", expr, re.I)
    return m.group(1).lower() if m else ""


def divisor_sites(text: str):
    """Yield (divisor_text,) for each `/` whose right operand is parenthesised.

    Only the parenthesised form is reported. A `/ x` by a bare name is a
    hazard too, but reporting every one of those buries the signal, and the
    docstring says so rather than pretending the pass is complete.
    """
    for m in re.finditer(r"/\s*\(", text):
        open_idx = m.end() - 1
        close = match_paren(text, open_idx)
        if close < 0:
            continue
        yield text[open_idx + 1 : close], open_idx


def scan(path: Path):
    stack: list[str] = []
    proc = "(file scope)"
    sites = []
    for lineno, text in logical_lines(path):
        low = text.lower()

        if _PROC_END_RE.match(low):
            # A procedure boundary clears the stack. Over-reporting means the
            # stack is never cleared EARLIER than the truth, so this is the one
            # place it is cleared without an explicit `end where`.
            stack.clear()
            proc = "(file scope)"
        pm = _PROC_START_RE.match(low)
        if pm:
            stack.clear()
            proc = text.strip().split("!")[0].strip()

        if _END_WHERE_RE.match(low):
            if stack:
                stack.pop()
            continue

        single_mask = None
        if _ELSEWHERE_RE.match(low):
            m = re.match(r"^\s*else\s*where\s*\(", low)
            if m:
                close = match_paren(text, m.end() - 1)
                if stack:
                    stack[-1] = ".not. " + text[m.end() : close]
            elif stack:
                stack[-1] = ".not. " + stack[-1]
        elif _WHERE_OPEN_RE.match(low):
            open_idx = low.index("(", low.index("where"))
            close = match_paren(text, open_idx)
            mask = text[open_idx + 1 : close] if close > 0 else "?"
            rest = text[close + 1 :].strip() if close > 0 else ""
            if rest:
                single_mask = mask
            else:
                stack.append(mask)
                continue

        masks = list(stack) + ([single_mask] if single_mask else [])
        if not masks:
            continue

        for m in _INTRINSIC_RE.finditer(text):
            name = m.group(1).lower()
            arg = arg_of(text, m.end())
            guarded = leading_call(arg) in _ARG_GUARDS or outer_call(
                text, m.start()
            ) in _ARG_GUARDS
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind=name,
                    guarded=guarded,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail=squeeze(arg)[:90],
                )
            )
        for m in _MERGE_RE.finditer(text):
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind="merge",
                    guarded=False,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail="",
                )
            )
        for div, _idx in divisor_sites(text):
            guarded = leading_call(div) in _DIV_GUARDS
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind="divide",
                    guarded=guarded,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail=squeeze(div)[:90],
                )
            )
    return sites


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--all",
        action="store_true",
        help="include the guarded sites and the classified survivors",
    )
    ap.add_argument("--tsv", action="store_true", help="one tab-separated row per site")
    ap.add_argument(
        "--kind",
        default="intrinsic",
        choices=("intrinsic", "divide", "merge", "any"),
        help="which class to report (default: the domain-sensitive intrinsics)",
    )
    args = ap.parse_args()

    sites = []
    for path in sorted(SRC.glob("*.f90")):
        sites.extend(scan(path))

    def wanted(s):
        if args.kind == "intrinsic":
            return s["kind"] in DOMAIN_INTRINSICS
        if args.kind == "any":
            return True
        return s["kind"] == args.kind

    sites = [s for s in sites if wanted(s)]
    for s in sites:
        s["classified"] = CLASSIFIED.get(
            (s["file"], s["proc"], s["kind"], s["detail"])
        )
    if not args.all:
        sites = [s for s in sites if not s["guarded"] and not s["classified"]]

    if args.tsv:
        for s in sites:
            print(
                "\t".join(
                    (
                        s["file"],
                        str(s["line"]),
                        s["kind"],
                        "guarded"
                        if s["guarded"]
                        else ("classified" if s["classified"] else "bare"),
                        s["proc"],
                        s["mask"],
                        s["detail"],
                    )
                )
            )
    else:
        by_file: dict[str, int] = {}
        for s in sites:
            by_file[s["file"]] = by_file.get(s["file"], 0) + 1
            print(f"{s['file']}:{s['line']}  {s['kind']}  [mask: {s['mask'][:70]}]")
            print(f"    {s['detail']}")
            if s["classified"]:
                print(f"    classified: {s['classified']}")
        print()
        for f, n in sorted(by_file.items(), key=lambda kv: -kv[1]):
            print(f"{n:5d}  {f}")
        print(f"{len(sites):5d}  TOTAL")
    if args.all:
        return 0
    return 1 if sites else 0


if __name__ == "__main__":
    raise SystemExit(main())
