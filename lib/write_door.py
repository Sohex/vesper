#!/usr/bin/env python3
"""The one refusal a generator of a linked artifact puts at its write door.

    python lib/write_door.py --self-test   # prove it both ways

`scripts/link_worktree.py` gives a worktree the ignored payload it cannot see
otherwise. Where a directory holds tracked content beside ignored payload it
cannot be one directory symlink, so the linker makes a real directory of one
link PER EXISTING FILE, and each of those links points INTO THE MAIN CHECKOUT.
A worktree that REGENERATES such a file therefore writes straight through, and
every artifact the main checkout built from those bytes is invalidated at a
moment nobody chose.

Measured on `vendor/lpj-guess/framework/vesper.h`: a worktree regenerated it,
correctly, and the write landed in the main checkout, whose LPJ-GUESS binary had
been built against the previous bytes. `build_lpj_guess.py --verify` began
refusing a binary nobody had touched. `notes/audits/worktree-write-through.md`
enumerates what else has that shape.

WHY A REFUSAL AND NOT A PERMISSION. A symlink carries no permissions of its own,
so the only mode that decides whether a write lands is the TARGET's, in the main
checkout, which is shared with every other worktree; a read-only link cannot be
built. The enforceable equivalent is a refusal at the point of writing, and it
needs one call per writer.

WHY THIS IS A MODULE AND NOT A FUNCTION EACH COMPONENT KEEPS. Seven groups of
regenerable artifact reach a linked path, in five components. Seven copies of
this is seven places for the walk, the message and the disposition to drift, and
a reader who finds one has no reason to think the others agree with it.

WHICH WRITERS CALL IT. Those that REGENERATE an artifact: the same path
rewritten from inputs, where the write replaces bytes another tree built from.
Not those that ACCUMULATE one -- `exoplasim/runs/` is linked whole on purpose so
that a climate run started in a worktree survives it, and writes there land in
the main checkout by arrangement. That line is the caller's to draw, which is
why the walk below reaches ancestors and has no opt-out.

WHY IT IS SEPARATE FROM `lib/paths.py`, whose `require_*` guards it otherwise
sits beside. Those ask whether an ARTIFACT is fit to read and open it to find
out; this asks where a PATH goes and touches no bytes. Keeping it dependency
free is what lets the cheapest writers in the tree call it -- a `.sra` text
writer and a C++ header generator both do, and neither should acquire netCDF4
to learn that its output is a link.

THIS IS NOT THE WHOLE GUARD, and the two halves are complements. A door
PREVENTS a write-through at the writer that would make it;
`scripts/check_worktree_links.py` REPORTS one that already landed, from a ledger
`link_worktree.py` records, and it reaches the groups that have no door because
they have no single writer to put one on. world-pkt5.
"""

from __future__ import annotations

from pathlib import Path

DEFAULT_INSTEAD = (
    "Write to a path inside this worktree, or run the generator in the main "
    "checkout.")


def write_through_link(path: Path | str) -> tuple[Path, Path] | None:
    """The link `path` is reached through and what it resolves to, or `None`.

    THE WALK REACHES ANCESTORS, not just the path. A per-file link makes the
    path itself a symlink, which is the shape every group in
    `notes/audits/worktree-write-through.md` has today; but whether a directory
    is linked per file or linked WHOLE is decided by whether any tracked file
    happens to sit in it, and that can change under a group without the group
    changing. Testing the path alone would make the guard depend on a fact
    nobody maintains, so it walks.

    THE COST OF WALKING IS THAT THIS CANNOT BE CALLED ON EVERY WRITE, and the
    caller contract is the other half of the design. A wholly-ignored directory
    is linked whole ON PURPOSE -- `exoplasim/runs/` is, so that a climate run
    started in a worktree survives the worktree -- and writes there land in the
    main checkout by arrangement. A writer of that kind does not call this door.
    The door is for a REGENERABLE artifact, where the same path is rewritten
    from inputs and the write replaces bytes another tree built from; a run
    directory is accumulated rather than regenerated, and that is the line.
    """
    path = Path(path)
    for candidate in (path, *path.parents):
        if candidate.is_symlink():
            return candidate, candidate.resolve()
        if candidate == candidate.parent:
            break
    return None


def refuse_a_write_through_a_symlink(path: Path | str, *,
                                     what: str = "a file other agents may be reading",
                                     instead: str = DEFAULT_INSTEAD) -> None:
    """Refuse to write `path` when it is a link into another checkout.

    `what` names the artifact in the message and `instead` says what the caller
    should do about it, because the disposition differs per writer: one takes an
    `--output`, another has a flag that skips the install step. The mechanism
    and the reason are the same everywhere, which is why they are here and not
    restated at each call.

    THERE IS NO OPT-OUT, deliberately. A write through the link is never the
    thing the caller wanted: it mutates a tree other agents are using, and
    nothing records which worktree did it. A caller that genuinely wants the
    main checkout's copy changed regenerates it there.
    """
    found = write_through_link(path)
    if found is None:
        return
    offender, target = found
    raise SystemExit(
        f"{path} is reached through the symlink {offender}, which points at "
        f"{target}.\n"
        f"  Writing here would write THROUGH the link and replace {what}, "
        f"invalidating whatever was built from those bytes, and nothing would "
        f"record which worktree did it.\n"
        f"  A worktree links this entry so it can READ it; generating into it "
        f"is what is refused. {instead}\n"
        f"  notes/audits/worktree-write-through.md.")


def _self_test() -> int:
    """Prove the refusal fires where it must and nowhere else.

    Every case has an answer that is known before it runs, and the negative
    cases are the point: a guard that refused every path would pass a suite of
    positive cases alone and would stop every generator in the main checkout.
    """
    import tempfile

    failures, ran = [], 0

    def case(name: str, got, want, why: str) -> None:
        nonlocal ran
        ran += 1
        ok = got == want
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name:<58} {why}")
        if not ok:
            failures.append(f"{name}: expected {want!r}, got {got!r}")

    def refused(path: Path) -> bool:
        try:
            refuse_a_write_through_a_symlink(path)
            return False
        except SystemExit:
            return True

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        main = root / "main"
        (main / "framework").mkdir(parents=True)
        real = main / "framework" / "generated.h"
        real.write_text("main checkout bytes\n", encoding="ascii")
        (main / "payload").mkdir()
        (main / "payload" / "big.nc").write_text("payload\n", encoding="ascii")

        wt = root / "worktree"
        (wt / "framework").mkdir(parents=True)
        # A PER-FILE link, which is the shape the linker makes in a directory
        # holding tracked content beside ignored payload.
        (wt / "framework" / "generated.h").symlink_to(real)
        # A WHOLLY-IGNORED directory, linked as one symlink. A write inside it
        # also lands in the main checkout. A writer that ACCUMULATES into such
        # a directory does not call this door at all; a writer that REGENERATES
        # a path under one is refused, because whether a directory is linked
        # whole or per file turns on whether a tracked file happens to sit in
        # it, which can change under a group without the group changing.
        (wt / "payload").symlink_to(main / "payload")
        (wt / "framework" / "own.h").write_text("mine\n", encoding="ascii")

        case("a per-file link is refused", refused(wt / "framework" / "generated.h"),
             True, "the measured vesper.h case")
        case("a real file beside it is not", refused(wt / "framework" / "own.h"),
             False, "otherwise the guard stops every generator everywhere")
        case("a file that does not exist yet is not",
             refused(wt / "framework" / "new.h"), False,
             "a first write into a real directory is the ordinary case")
        case("a path under a linked directory is refused",
             refused(wt / "payload" / "big.nc"), True,
             "the walk has to reach ancestors, not only the path")
        case("and so is one that does not exist yet",
             refused(wt / "payload" / "new.nc"), True,
             "the link is crossed whether or not the leaf is there")
        case("the main checkout's own copy is not refused", refused(real), False,
             "nothing on these paths is a link there, so it never fires")

        case("the offender named is the link and not the target",
             write_through_link(wt / "payload" / "big.nc"),
             (wt / "payload", main / "payload"),
             "a reader has to be told which link to stop crossing")

        # The refusal must SAY the two things a reader needs: where the write
        # would land, and what to do instead. A guard that refuses without
        # either sends the reader to read this module.
        try:
            refuse_a_write_through_a_symlink(
                wt / "framework" / "generated.h", what="a generated header",
                instead="Pass --no-install.")
            message = ""
        except SystemExit as refusal:
            message = str(refusal)
        case("the message carries the target and the caller's disposition",
             [str(real) in message, "a generated header" in message,
              "--no-install" in message], [True, True, True],
             "the disposition differs per writer and the mechanism does not")

        # The write the guard prevented must still be POSSIBLE, or the repair
        # is a component nobody can run from a worktree.
        elsewhere = wt / "framework" / "own.h"
        refuse_a_write_through_a_symlink(elsewhere)
        elsewhere.write_text("regenerated\n", encoding="ascii")
        case("a refused generator can still write inside the worktree",
             [elsewhere.read_text(encoding="ascii"),
              real.read_text(encoding="ascii")],
             ["regenerated\n", "main checkout bytes\n"],
             "the main checkout's bytes are untouched and the work is not lost")

    print(f"\n{ran} cases, {len(failures)} failed")
    for line in failures:
        print(f"  {line}")
    return 1 if failures else 0


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--self-test", action="store_true",
                        help="prove the refusal fires where it must and "
                             "nowhere else")
    args = parser.parse_args()
    if not args.self_test:
        parser.error("this module is a guard other writers call; "
                     "--self-test is the only thing it does on its own")
    sys.exit(_self_test())
