#!/usr/bin/env python3
"""Writes to the shared diagnostics unit that a thread other than 0 can reach.

    python exoplasim/scripts/lint_diag_writes.py           # the unanswered sites
    python exoplasim/scripts/lint_diag_writes.py --all     # every site with its verdict
    python exoplasim/scripts/lint_diag_writes.py --tsv     # machine-readable
    python exoplasim/scripts/lint_diag_writes.py --procs   # the call-graph verdict per procedure

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. `opendiag` in `plasim.f90` opens `nud` under
`if (mypid == NROOT)`, and `nud` is 6. A Fortran unit belongs to the PROCESS,
so that one `open` makes unit 6 the diagnostics file for the whole thread team,
and `write(*,...)` lands there too because `*` is unit 6 as well. `plasim.f90`
then opens `!$omp parallel num_threads(NPRO)` around `mpstart` through
`mpstop`, so every statement in the model's call tree is executed by NPRO
threads unless something on the path tests `mypid == NROOT`.

Under the deleted MPI arm a non-root rank's `write(nud,...)` went to that
rank's own stdout and out of the record, so an unguarded site was invisible.
Under threads it lands in `plasim_diag` beside root's, NPRO interleaved copies,
in an order libgfortran's per-statement unit lock does not fix. That makes the
diagnostics file non-deterministic run to run, and
`reproducibility_matrix.py` names a `plasim_diag`-only disagreement as this
candidate first.

world-0ihs read all 928 sites by hand, found thirteen holes and closed them.
This is the same two passes written down, so the next unguarded site fails a
gate instead of waiting for the next reader.
`notes/audits/shared-diagnostics-unit-under-threads.md` is the finding.

WHAT IT PERMITS. A `write(nud,...)` or `write(*,...)` in `plasim/src` passes on
any one of four answers, and on nothing else:

  1. An enclosing `if` branch tests `mypid == NROOT`. One thread reaches the
     statement, which is the answer for a quantity that is the same on every
     thread.
  2. An enclosing `if` branch tests `npro == 1`. `NPRO` is the thread count, so
     that admits exactly the single-thread case and there is no second writer.
  3. The statement is inside an `!$omp critical` region. The answer for a
     quantity that is the THREAD's own chunk, where a root guard would silence
     exactly the cells worth seeing, and for a leaf module with no `pumamod`
     dependency and so no `mypid` to test.
  4. The enclosing procedure is reached only through answers 1 to 3. That is
     the call-graph pass: `print_planet`, `readnl`, `wrspam` and the rest write
     unguarded because every `call` that reaches them is already under a root
     guard.

WHAT IT FORBIDS is everything else, which is one thing: a write that NPRO
threads reach with nothing between them and the unit.

THE PASSES ARE LEXICAL, AND OVER-REPORT DELIBERATELY. Five places they err
towards reporting a site that is answered, and none towards clearing one that
is not:

  1. The block stack is popped only by the `end` form that matches its top. An
     `end` this parser cannot match to the construct it closes CLEARS the whole
     stack, which drops a guard the source still has and never adds one.
  2. A condition containing `.or.` never counts as a guard, even when one arm
     of it is the root test, because the other arm may admit every thread.
  3. Only `mypid == NROOT` and `mypid .eq. NROOT` are read as the root test,
     in any spacing and either case. A test written some other way is reported.
  4. A procedure with no `call` site anywhere in `plasim/src` is NOT root-only.
     That covers the program unit, the entry points the parallel region opens
     on, and anything reached from outside this directory.
  5. `interface` bodies and derived-type definitions are walked like any other
     block rather than skipped, so a construct inside one cannot hide a site.

The one direction it can err the other way is a call CYCLE that no seed
reaches: mutually recursive procedures called from nowhere else keep the
benefit of the doubt. That is unreachable code by construction, and this tree
has none.

WHAT IT DOES NOT KNOW. Whether a site that IS answered is answered CORRECTLY --
whether the quantity really is the same on every thread, or really is the
thread's own chunk. That is the reading world-0ihs did, and the audit note
carries the verdict per site. This gate holds the line the reading established.
"""

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import MODEL_SRC  # noqa: E402

SRC = MODEL_SRC / "plasim" / "src"

# The shared unit, in every spelling that reaches it. `nud` is 6 and `*` is 6,
# so both land in plasim_diag. A `print` statement would too, and the tree has
# none; it is matched anyway so that adding one does not slip past.
_WRITE_RE = re.compile(r"^\s*write\s*\(\s*(nud|\*)\s*[,)]", re.I)
_PRINT_RE = re.compile(r"^\s*print\s*[*'\"]", re.I)

# The root test, in the two spellings and any spacing. `/=` and `.ne.` do not
# match, which is the direction that matters: a negative test is not a guard.
_ROOT_RE = re.compile(r"mypid\s*(?:==|\.eq\.)\s*nroot", re.I)
# NPRO is the thread count, so `npro == 1` admits exactly the serial case.
_SERIAL_RE = re.compile(r"npro\s*(?:==|\.eq\.)\s*1(?![0-9.])", re.I)
_OR_RE = re.compile(r"\.or\.", re.I)

_IF_THEN_RE = re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?if\s*\(", re.I)
_ELSE_IF_RE = re.compile(r"^\s*else\s*if\s*\(", re.I)
_ELSE_RE = re.compile(r"^\s*else\b", re.I)
_PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:real|integer|logical|complex|double\s+precision|character)"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function|program)\s+([a-z_]\w*)",
    re.I,
)
_PROC_END_RE = re.compile(r"^\s*end\s*(subroutine|function|program)\b", re.I)
_MODULE_RE = re.compile(r"^\s*module\s+([a-z_]\w*)\s*$", re.I)
_CALL_RE = re.compile(r"(?<![a-z0-9_])call\s+([a-z_]\w*)", re.I)
_OMP_CRITICAL_RE = re.compile(r"^\s*!\$omp\s+critical\b", re.I)
_OMP_END_CRITICAL_RE = re.compile(r"^\s*!\$omp\s+end\s+critical\b", re.I)

# Block openers that are not `if`, each with the `end` word that closes it.
# They carry no guard; they are tracked only so the stack stays in step and an
# `end if` is never matched against something that is not an `if`.
_OPENERS = (
    (re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?do\b", re.I), "do"),
    (re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?select\s*(?:case|type)\b", re.I), "select"),
    (re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?forall\s*\(", re.I), "forall"),
    (re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?associate\s*\(", re.I), "associate"),
    (re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?block\s*$", re.I), "block"),
    (re.compile(r"^\s*interface\b", re.I), "interface"),
    (re.compile(r"^\s*type\s*(?:,[^:]*)?::\s*[a-z_]\w*\s*$", re.I), "type"),
    (re.compile(r"^\s*type\s+[a-z_]\w*\s*$", re.I), "type"),
)
_END_RE = re.compile(
    r"^\s*end\s*(if|do|select|where|forall|associate|block|interface|type)\b", re.I
)
_BARE_END_RE = re.compile(r"^\s*end\s*$", re.I)
# A `where` construct: `where (mask)` with nothing after the closing paren.
_WHERE_RE = re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?where\s*\(", re.I)


# ---------------------------------------------------------------------------
# THE CLASSIFIED SURVIVORS, and why each is not a defect.
#
# A site belongs here only when the guard on it is REAL and this pass cannot
# see it, never to grandfather a site that is genuinely unguarded. Every write
# in the tree today is answered by one of the four rules above, so the table is
# empty; it exists because the first site whose answer the parser cannot
# express should be argued in one line here rather than by loosening a rule
# that then stops catching the next one.
#
# The key is (file, procedure, the statement text squeezed to one line), so a
# changed statement drops out of the table and is reported again.
# ---------------------------------------------------------------------------

CLASSIFIED: dict[tuple[str, str, str], str] = {}


def squeeze(text: str) -> str:
    return " ".join(text.split())


def strip_comment(line: str) -> str:
    """Drop a trailing `!` comment, respecting quoted strings.

    An `!$omp` directive is a comment to the compiler and code to us, so it is
    handed back whole rather than stripped.
    """
    if re.match(r"^\s*!\$omp\b", line, re.I):
        return line.rstrip()
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


def logical_lines(source: str):
    """Yield (first_physical_lineno, joined_text) with continuations folded.

    A continued statement is one site and cannot hide half of itself, and a
    guard split across two lines is still read as one condition.
    """
    raw = source.splitlines()
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


def split_top(text: str, sep: str) -> list[str]:
    """Split on `sep` at parenthesis depth zero, case-insensitively."""
    parts = []
    depth = 0
    cur = ""
    quote = None
    i = 0
    n = len(sep)
    while i < len(text):
        ch = text[i]
        if quote:
            cur += ch
            if ch == quote:
                quote = None
            i += 1
            continue
        if ch in "'\"":
            quote = ch
            cur += ch
            i += 1
            continue
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if depth == 0 and text[i : i + n].lower() == sep:
            parts.append(cur)
            cur = ""
            i += n
            continue
        cur += ch
        i += 1
    parts.append(cur)
    return parts


def unwrap(text: str) -> str:
    """Strip one layer of enclosing parentheses, repeatedly."""
    t = text.strip()
    while t.startswith("(") and match_paren(t, 0) == len(t) - 1:
        t = t[1:-1].strip()
    return t


def _is_root_term(term: str) -> bool:
    t = unwrap(term)
    return bool(re.fullmatch(r"mypid\s*(?:==|\.eq\.)\s*nroot", t, re.I))


def _is_serial_term(term: str) -> bool:
    t = unwrap(term)
    return bool(re.fullmatch(r"npro\s*(?:==|\.eq\.)\s*1", t, re.I))


def guard_of(cond: str) -> str:
    """Which of the recognised guards a condition is, or the empty string.

    Read from the boolean structure at depth zero rather than by searching the
    text, because `.and.` binds tighter than `.or.` in Fortran and the two need
    opposite treatment. A condition guards when EVERY top-level `.or.` disjunct
    carries the test as one of its top-level `.and.` conjuncts: then the test
    holds whenever the condition does, however the rest of it reads.

    The conjunct has to BE the test, not merely contain it, so `.not. (mypid ==
    NROOT)` and `f(mypid == NROOT)` are not guards. A single disjunct with no
    `.or.` is the ordinary case and falls out of the same rule.
    """
    disjuncts = split_top(cond, ".or.")
    for check, name in ((_is_root_term, "mypid == NROOT"),
                        (_is_serial_term, "npro == 1")):
        if all(any(check(t) for t in split_top(d, ".and.")) for d in disjuncts):
            return name
    return ""


def literal_constraints(conds: list[str], dummies: list[str]) -> dict[int, int]:
    """Dummy-argument positions a guard chain pins to an integer literal.

    A branch taken on `if (kpid == 0)`, where `kpid` is a dummy argument, is
    reachable only from the call sites that pass 0 there. That is the shape
    `calini`'s summary block has: `prolog` passes -1 and `readnl` passes 0, so
    the block is reached from `readnl` alone, and `readnl` is root-only.

    Only a bare `name == integer` conjunct at depth zero counts, and only when
    `name` is a dummy of the enclosing procedure. Anything else contributes no
    constraint, which leaves the site resting on the call-graph pass alone.
    """
    out: dict[int, int] = {}
    lower = [d.lower() for d in dummies]
    for cond in conds:
        for d in split_top(cond, ".or."):
            if len(split_top(cond, ".or.")) > 1:
                break
            for term in split_top(d, ".and."):
                m = re.fullmatch(
                    r"\s*([a-z_]\w*)\s*(?:==|\.eq\.)\s*([+-]?\d+)\s*",
                    unwrap(term), re.I,
                )
                if m and m.group(1).lower() in lower:
                    out[lower.index(m.group(1).lower())] = int(m.group(2))
    return out


def if_parts(text: str):
    """Split `if (cond) tail` into (cond, tail, opens_block).

    `opens_block` is true when the statement ends in `then`, so the guard
    covers a construct rather than the one statement after it.
    """
    m = re.search(r"if\s*\(", text, re.I)
    if not m:
        return "", "", False
    open_idx = m.end() - 1
    close = match_paren(text, open_idx)
    if close < 0:
        return "", "", False
    cond = text[open_idx + 1 : close]
    tail = text[close + 1 :].strip()
    return cond, tail, bool(re.match(r"^then\b", tail, re.I))


def dummy_args(header: str) -> list[str]:
    """The dummy argument names on a `subroutine` / `function` header."""
    m = re.search(r"\(", header)
    if not m:
        return []
    close = match_paren(header, m.start())
    if close < 0:
        return []
    return [a.strip() for a in split_top(header[m.start() + 1 : close], ",") if a.strip()]


def actual_args(text: str, call_end: int) -> list[str]:
    """The actual argument expressions on a `call name(...)` statement."""
    rest = text[call_end:]
    m = re.match(r"\s*\(", rest)
    if not m:
        return []
    open_idx = call_end + m.end() - 1
    close = match_paren(text, open_idx)
    if close < 0:
        return []
    return [a.strip() for a in split_top(text[open_idx + 1 : close], ",")]


def scan(name: str, source: str):
    """Pass 1: every write site and call site, with the guards enclosing it.

    Takes the source TEXT rather than a path so the fixtures below can run the
    same parser the tree does.
    """
    stack: list[tuple[str, str, str]] = []  # (kind, guard, condition)
    proc = None
    dummies: list[str] = []
    writes = []
    calls = []

    def guards() -> list[str]:
        return [g for _k, g, _c in stack if g]

    def conds() -> list[str]:
        return [c for k, _g, c in stack if k == "if" and c]

    for lineno, text in logical_lines(source):
        low = text.lower()

        if _OMP_END_CRITICAL_RE.match(low):
            if stack and stack[-1][0] == "critical":
                stack.pop()
            else:
                stack.clear()
            continue
        if _OMP_CRITICAL_RE.match(low):
            stack.append(("critical", "!$omp critical", ""))
            continue
        if low.lstrip().startswith("!$omp"):
            continue

        if _PROC_END_RE.match(low) or _BARE_END_RE.match(low):
            stack.clear()
            proc = None
            continue
        pm = _PROC_START_RE.match(low)
        if pm and "end" not in low.split()[0]:
            stack.clear()
            proc = pm.group(2).lower()
            dummies = dummy_args(text[pm.end(2) :])
            continue
        if _MODULE_RE.match(low):
            stack.clear()
            proc = None
            continue

        em = _END_RE.match(low)
        if em:
            want = em.group(1).lower()
            if stack and stack[-1][0] == want:
                stack.pop()
            else:
                # Cannot match this `end` to its opener. Clearing the stack
                # drops guards the source still has and never invents one,
                # which is the safe direction.
                stack.clear()
            continue

        # `else` / `else if` replace the top frame's guard rather than nesting.
        if _ELSE_IF_RE.match(low):
            cond, _tail, _opens = if_parts(text)
            if stack and stack[-1][0] == "if":
                stack[-1] = ("if", guard_of(cond), cond)
            else:
                stack.clear()
            continue
        if _ELSE_RE.match(low) and not re.match(r"^\s*elsewhere\b", low):
            if stack and stack[-1][0] == "if":
                stack[-1] = ("if", "", "")
            else:
                stack.clear()
            continue
        if re.match(r"^\s*else\s*where\b", low):
            continue

        # A `case` label inside a `select` opens no frame and carries no guard.
        if re.match(r"^\s*(case|class\s+is|class\s+default|type\s+is)\b", low):
            continue

        # `if (cond) then` opens a frame; `if (cond) stmt` guards that one
        # statement, so the guard is applied to the tail below rather than
        # pushed.
        inline_guard = ""
        inline_cond = ""
        if _IF_THEN_RE.match(low):
            cond, tail, opens = if_parts(text)
            if opens:
                stack.append(("if", guard_of(cond), cond))
                continue
            inline_guard = guard_of(cond)
            inline_cond = cond
            text, low = tail, tail.lower()

        opened = False
        for rx, kind in _OPENERS:
            if rx.match(low):
                stack.append((kind, "", ""))
                opened = True
                break
        if opened:
            continue
        if _WHERE_RE.match(low):
            m = re.search(r"where\s*\(", low)
            close = match_paren(text, m.end() - 1)
            if close > 0 and not text[close + 1 :].strip():
                stack.append(("where", "", ""))
                continue

        here = guards() + ([inline_guard] if inline_guard else [])
        if proc is None:
            # A statement outside any procedure has no thread of its own.
            continue

        if _WRITE_RE.match(low) or _PRINT_RE.match(low):
            writes.append(
                dict(
                    file=name,
                    line=lineno,
                    proc=proc,
                    guards=here,
                    conds=conds() + ([inline_cond] if inline_cond else []),
                    dummies=list(dummies),
                    text=squeeze(text)[:110],
                )
            )
        for m in _CALL_RE.finditer(text):
            calls.append(
                dict(
                    caller=proc,
                    callee=m.group(1).lower(),
                    answered=bool(here),
                    args=actual_args(text, m.end()),
                )
            )
    return writes, calls


def single_writer_procs(calls) -> set[str]:
    """Pass 2: the procedures only one thread, or one at a time, can be in.

    Fixed point of "every call site is answered, or is in a procedure that is
    itself answered". Computed as a forward propagation of the NEGATIVE, from
    two kinds of seed: a procedure with no call site anywhere in `plasim/src`,
    and a procedure with an unanswered call site in a procedure already known
    not to qualify. Propagating the negative rather than assuming it means the
    answer does not depend on iteration order.
    """
    callers_of: dict[str, list[tuple[str, bool]]] = defaultdict(list)
    callees_of: dict[str, set[str]] = defaultdict(set)
    procs = set()
    for c in calls:
        callers_of[c["callee"]].append((c["caller"], c["answered"]))
        callees_of[c["caller"]].add(c["callee"])
        procs.add(c["caller"])
        procs.add(c["callee"])

    # A procedure this directory never calls is an entry point as far as this
    # pass can tell, so it gets no benefit of the doubt.
    bad = {p for p in procs if not callers_of[p]}
    work = deque(bad)
    while work:
        p = work.popleft()
        for callee in callees_of[p]:
            if callee in bad:
                continue
            if any(caller == p and not answered for caller, answered in callers_of[callee]):
                bad.add(callee)
                work.append(callee)
    return procs - bad


def verdicts(writes, calls):
    """Annotate every write site with the answer that clears it, or "".

    Pass 2 lives here rather than in `main` so that the fixtures below run the
    same code the tree does. Returns the set of root-only procedures alongside,
    which `--procs` prints.
    """
    root_only = single_writer_procs(calls)
    sites_of: dict[str, list[dict]] = defaultdict(list)
    for c in calls:
        sites_of[c["callee"]].append(c)

    def reached_only_by_guarded_calls(w) -> str:
        """Whether the call sites this branch is reachable FROM are all answered.

        Without a literal constraint this is the plain call-graph question and
        `root_only` has already answered it. With one, the population narrows
        to the call sites that pass the pinned literal, and a branch reachable
        from none of them is reachable at all only from outside `plasim/src`,
        which this pass does not clear.
        """
        pinned = literal_constraints(w["conds"], w["dummies"])
        if not pinned:
            return ""
        reaching = []
        for c in sites_of[w["proc"]]:
            for idx, want in pinned.items():
                if idx < len(c["args"]):
                    actual = c["args"][idx].strip()
                    if re.fullmatch(r"[+-]?\d+", actual) and int(actual) != want:
                        break
            else:
                reaching.append(c)
        if not reaching:
            return ""
        if all(c["answered"] or c["caller"] in root_only for c in reaching):
            names = sorted({c["caller"] for c in reaching})
            pins = ", ".join(f"{w['dummies'][i]} == {v}" for i, v in sorted(pinned.items()))
            return f"reached under {pins} only from {', '.join(names)}, all guarded"
        return ""

    for s in writes:
        s["classified"] = CLASSIFIED.get((s["file"], s["proc"], s["text"]))
        if s["guards"]:
            s["verdict"] = s["guards"][-1]
        elif s["proc"] in root_only:
            s["verdict"] = "call graph: every path is guarded"
        elif reached_only_by_guarded_calls(s):
            s["verdict"] = reached_only_by_guarded_calls(s)
        elif s["classified"]:
            s["verdict"] = s["classified"]
        else:
            s["verdict"] = ""
    return root_only


# ---------------------------------------------------------------------------
# THE FIXTURES. Reduced sources run on every invocation, each built to be right
# or wrong in a NAMED way, because a gate whose only evidence is that the tree
# passes cannot distinguish "no unguarded write" from "no pass".
#
# Each entry is (name, source, expected number of unguarded sites).
# ---------------------------------------------------------------------------

FIXTURES: tuple[tuple[str, str, int], ...] = (
    ("a bare write is reported", """
      subroutine leaf
      write(nud,*) 'x'
      end subroutine leaf
""", 1),
    ("a root guard clears it", """
      subroutine leaf
      if (mypid == NROOT) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("the .eq. spelling clears it", """
      subroutine leaf
      if (mypid .eq. NROOT) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("a NEGATED root test does not clear it", """
      subroutine leaf
      if (mypid /= NROOT) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 1),
    ("the else arm of a root guard is not guarded", """
      subroutine leaf
      if (mypid == NROOT) then
      write(nud,*) 'x'
      else
      write(nud,*) 'y'
      endif
      end subroutine leaf
""", 1),
    ("a root test conjoined with anything still clears it", """
      subroutine leaf
      if (mypid == NROOT .and. n == 1 .and. (a == 0 .or. b == 0)) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("a root test DISJOINED with anything does not", """
      subroutine leaf
      if (n == 1 .or. mypid == NROOT) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 1),
    ("the guard does not survive the end of its if", """
      subroutine leaf
      if (mypid == NROOT) then
      n = 1
      endif
      write(nud,*) 'x'
      end subroutine leaf
""", 1),
    ("npro == 1 clears it", """
      subroutine leaf
      if (npro == 1) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("an omp critical clears it", """
      subroutine leaf
!$omp critical (nudwrite)
      write(nud,*) 'x'
!$omp end critical (nudwrite)
      end subroutine leaf
""", 0),
    ("the critical does not survive its end", """
      subroutine leaf
!$omp critical (nudwrite)
      n = 1
!$omp end critical (nudwrite)
      write(nud,*) 'x'
      end subroutine leaf
""", 1),
    ("a one-line if guards the statement after it", """
      subroutine leaf
      if (mypid == NROOT) write(nud,*) 'x'
      write(nud,*) 'y'
      end subroutine leaf
""", 1),
    ("write(*,...) and print are the same unit", """
      subroutine leaf
      write(*,*) 'x'
      print *,'y'
      end subroutine leaf
""", 2),
    ("a guard split across a continuation is still one condition", """
      subroutine leaf
      if (mypid == NROOT .and.                                          &
     &    n == 1) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("root-only by caller clears the callee", """
      subroutine top
      if (mypid == NROOT) then
      call leaf
      endif
      end subroutine top
      subroutine leaf
      write(nud,*) 'x'
      end subroutine leaf
""", 0),
    ("one unguarded caller is enough to report the callee", """
      subroutine top
      if (mypid == NROOT) then
      call leaf
      endif
      end subroutine top
      subroutine other
      call leaf
      end subroutine other
      subroutine leaf
      write(nud,*) 'x'
      end subroutine leaf
""", 1),
    ("a dummy pinned to a literal narrows the callers", """
      subroutine top
      if (mypid == NROOT) then
      call leaf(0)
      endif
      end subroutine top
      subroutine other
      call leaf(-1)
      end subroutine other
      subroutine leaf(kpid)
      if (kpid == 0) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 0),
    ("the same pin the other way round is reported", """
      subroutine top
      if (mypid == NROOT) then
      call leaf(0)
      endif
      end subroutine top
      subroutine other
      call leaf(-1)
      end subroutine other
      subroutine leaf(kpid)
      if (kpid == -1) then
      write(nud,*) 'x'
      endif
      end subroutine leaf
""", 1),
    ("a nested do does not swallow the enclosing guard", """
      subroutine leaf
      if (mypid == NROOT) then
      do j = 1, 10
      write(nud,*) 'x'
      enddo
      endif
      end subroutine leaf
""", 0),
    ("a select case inside a guard keeps it", """
      subroutine leaf
      if (mypid == NROOT) then
      select case (n)
      case (1)
      write(nud,*) 'x'
      end select
      endif
      end subroutine leaf
""", 0),
)


def selftest() -> list[str]:
    bad = []
    for name, source, want in FIXTURES:
        writes, calls = scan("fixture.f90", source)
        verdicts(writes, calls)
        got = len([w for w in writes if not w["verdict"]])
        if got != want:
            bad.append(f"{name}: expected {want} unguarded, got {got}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--all", action="store_true", help="print every site with its verdict"
    )
    ap.add_argument("--tsv", action="store_true", help="one tab-separated row per site")
    ap.add_argument(
        "--procs",
        action="store_true",
        help="print the call-graph verdict for every procedure that writes",
    )
    args = ap.parse_args()

    broken = selftest()
    if broken:
        print("the pass does not answer its own fixtures, so its verdict on the "
              "tree means nothing:")
        for b in broken:
            print(f"    {b}")
        return 1

    writes = []
    calls = []
    for path in sorted(SRC.glob("*.f90")):
        w, c = scan(path.name, path.read_text(errors="replace"))
        writes.extend(w)
        calls.extend(c)
    root_only = verdicts(writes, calls)

    if args.procs:
        by_proc: dict[str, int] = defaultdict(int)
        for s in writes:
            by_proc[s["proc"]] += 1
        for p in sorted(by_proc):
            mark = "root-only" if p in root_only else "every thread"
            print(f"{by_proc[p]:5d}  {p:28s}  {mark}")
        return 0

    shown = writes if args.all else [s for s in writes if not s["verdict"]]

    if args.tsv:
        for s in shown:
            print(
                "\t".join(
                    (s["file"], str(s["line"]), s["proc"], s["verdict"] or "UNGUARDED",
                     s["text"])
                )
            )
    else:
        by_file: dict[str, int] = {}
        for s in shown:
            by_file[s["file"]] = by_file.get(s["file"], 0) + 1
            print(f"{s['file']}:{s['line']}  {s['proc']}")
            print(f"    {s['text']}")
            print(f"    verdict: {s['verdict'] or 'UNGUARDED -- every thread reaches this'}")
        print()
        for f, n in sorted(by_file.items(), key=lambda kv: -kv[1]):
            print(f"{n:5d}  {f}")
        print(f"{len(shown):5d}  {'POPULATION' if args.all else 'UNGUARDED'}")
    if args.all:
        return 0
    return 1 if shown else 0


if __name__ == "__main__":
    raise SystemExit(main())
