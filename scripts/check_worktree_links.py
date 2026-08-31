#!/usr/bin/env python3
"""Has this worktree written through one of its links into the main checkout?

    python scripts/check_worktree_links.py            # exit 1 if a target moved
    python scripts/check_worktree_links.py --quiet     # say nothing when clean
    python scripts/check_worktree_links.py --self-test # prove it both ways

`scripts/link_worktree.py` links a worktree's ignored payload back to the main
checkout, and where a directory holds tracked content beside ignored payload it
cannot make one directory symlink -- so it links PER FILE, and each of those
links is a symlink INTO THE MAIN CHECKOUT. A worktree that REGENERATES such a
file therefore writes straight through, and every artifact the main checkout
built from those bytes is invalidated at a moment nobody chose.

Measured 2026-08-31 on `vendor/lpj-guess/framework/vesper.h`, which is
generated, ignored, and sits beside tracked source. An agent regenerated it in
a worktree -- correctly, to take a generation timestamp out of a header whose
hash gates the binary -- and the write landed in the main checkout, whose
LPJ-GUESS binary had been built against the previous bytes.
`build_lpj_guess.py --verify` began refusing a binary nobody had touched.

WHY THIS EXISTS SEPARATELY FROM `link_worktree.py --check`, which asks the same
question. The row's acceptance is that the write reach THE AGENT THAT MADE IT,
and an agent regenerating a file does not then run the linker. What it does do
is commit, and `smoke_test.py` runs before every commit. This is a static read
of one JSON ledger plus a hash of the link targets under the ledger's budget --
74 MB of the 11.6 GB, by construction -- so it belongs in that tier.

`--self-test` builds a throwaway repository of the same shape -- a generated
ignored file beside tracked source, and a wholly-ignored directory -- adds a
REAL git worktree to it, and exercises the ledger against a genuine write
through the link and against a legitimate regeneration, which must remain
possible. Temporary directories only; it touches no worktree of this project.

IT IS A NO-OP OUTSIDE A WORKTREE, and a no-op in a worktree that has no ledger
yet. The second matters: worktrees standing when the ledger was introduced have
none, and failing them all at their next commit would punish agents for the
absence of a record they could not have written. They get a line telling them to
run the linker, which is what establishes one.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from link_worktree import (  # noqa: E402
    LEDGER_HASH_BUDGET, changed_since_link, read_ledger)

FIX = "python scripts/link_worktree.py"


def worktree_and_main() -> tuple[Path, Path] | None:
    """(worktree, main checkout), or None if this is not a worktree at all."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--path-format=absolute", "--show-toplevel",
             "--git-common-dir"],
            cwd=Path(__file__).resolve().parent, check=True,
            capture_output=True, text=True).stdout.split()
    except (subprocess.CalledProcessError, OSError):
        return None
    if len(out) != 2:
        return None
    top, common = Path(out[0]), Path(out[1])
    main = common.parent
    return None if top == main else (top, main)


def check_worktree_links() -> list[str]:
    """Failures, in `smoke_test.py`'s vocabulary. Empty everywhere else."""
    located = worktree_and_main()
    if located is None:
        return []
    wt, main = located
    if read_ledger(wt) is None:
        return []
    changed, linked_at = changed_since_link(wt, main)
    if not changed:
        return []
    return [
        f"{len(changed)} linked target(s) in {main} have changed since this "
        f"worktree was linked ({linked_at}). A per-file link points INTO the "
        "main checkout, so either something here regenerated one and wrote "
        "through, invalidating whatever the main checkout built from it, or "
        "the main checkout regenerated it and this worktree's results came "
        "from bytes that are gone. Settle which, rebuild what depended on "
        f"those bytes, then `{FIX}` to re-baseline."
    ] + [f"    {c}" for c in changed]


def self_test() -> list[str]:
    """Build a real worktree of a throwaway repo and prove the ledger both ways.

    Class 17: every case has a right answer. The two that matter are paired --
    a write through the link MUST be caught, and the regeneration that caused
    it MUST still have been possible, because a guard that made regenerating a
    generated file impossible would have replaced one defect with another.
    """
    import tempfile

    here = Path(__file__).resolve().parent
    bad: list[str] = []

    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp) / "repo"
        (repo / "scripts").mkdir(parents=True)
        (repo / "pkg").mkdir()
        (repo / "payload").mkdir()
        for name in ("link_worktree.py", "check_worktree_links.py"):
            (repo / "scripts" / name).write_text((here / name).read_text(),
                                                 encoding="utf-8")
        (repo / "pkg/source.c").write_text("int main(void){return 0;}\n")
        (repo / "pkg/gen.h").write_text("#define GEN 1\n")       # the vesper.h shape
        (repo / "pkg/small.h").write_text("v1\n")                # under the budget
        (repo / "pkg/big.bin").write_bytes(                      # over the budget
            b"x" * (LEDGER_HASH_BUDGET + 1))
        (repo / "payload/out.dat").write_text("run output\n")
        # Without a trailing slash, so the pattern covers the symlink that
        # stands in for the directory as well as the directory itself.
        (repo / ".gitignore").write_text(
            "pkg/gen.h\npkg/small.h\npkg/big.bin\npayload\nscripts/__pycache__\n")

        def git(*args, cwd=repo):
            subprocess.run(["git", *args], cwd=cwd, check=True,
                           capture_output=True, text=True)

        git("init", "-q", ".")
        git("config", "user.email", "self-test@invalid")
        git("config", "user.name", "self-test")
        git("add", ".gitignore", "pkg/source.c", "scripts/link_worktree.py",
            "scripts/check_worktree_links.py")
        git("commit", "-qm", "init")
        wt = repo / "wt"
        git("worktree", "add", "-q", str(wt), "-b", "self-test")

        def link(*extra):
            return subprocess.run(
                [sys.executable, str(repo / "scripts/link_worktree.py"),
                 "--worktree", str(wt), *extra],
                cwd=repo, capture_output=True, text=True)

        def check(tree=wt):
            return subprocess.run(
                [sys.executable, str(tree / "scripts/check_worktree_links.py")],
                cwd=repo, capture_output=True, text=True)

        link()
        ledger = repo / ".git/worktrees/wt/link-ledger.json"
        if not ledger.exists():
            return ["linking wrote no ledger, so nothing below can be tested"]
        if not (wt / "pkg/gen.h").is_symlink():
            return ["the probe did not reproduce the shape: pkg/gen.h is not a "
                    "per-file link in the worktree"]
        if not (wt / "payload").is_symlink():
            return ["the probe did not reproduce the shape: payload is not a "
                    "directory link in the worktree"]
        if check().returncode != 0:
            bad.append("a freshly linked worktree was reported as changed")

        # A GENUINE WRITE-THROUGH.
        (wt / "pkg/gen.h").write_text("#define GEN 2\n")
        if (repo / "pkg/gen.h").read_text() != "#define GEN 2\n":
            return ["the probe did not reproduce the hazard: the worktree's "
                    "write did not reach the main checkout"]
        r = link("--check")
        if r.returncode == 0 or "pkg/gen.h: CONTENT CHANGED" not in r.stdout:
            bad.append("--check missed a write through a per-file link")
        r = check()
        if r.returncode == 0 or "pkg/gen.h: CONTENT CHANGED" not in r.stdout:
            bad.append("the commit-time check missed a write through a link")

        # AND THE LEGITIMATE REGENERATION SURVIVED IT.
        if (repo / "pkg/gen.h").read_text() != "#define GEN 2\n":
            bad.append("the regenerated content did not survive the check")
        before = ledger.read_text(encoding="utf-8")
        link()
        if ledger.read_text(encoding="utf-8") == before:
            bad.append("re-linking left the ledger untouched, so a legitimate "
                       "regeneration can never be accepted")
        if check().returncode != 0:
            bad.append("re-linking did not re-baseline the ledger, so a "
                       "legitimate regeneration can never be accepted")

        # A rewrite with IDENTICAL bytes, under the budget, is not a finding.
        (repo / "pkg/small.h").write_text("v1\n")
        if check().returncode != 0:
            bad.append("an identical-bytes rewrite under the hash budget was "
                       "reported; the digest comparison is not being used")

        # Above the budget there is no digest, so a rewrite reports as one.
        (repo / "pkg/big.bin").touch()
        r = check()
        if r.returncode == 0 or "pkg/big.bin: REWRITTEN" not in r.stdout:
            bad.append("a rewrite above the hash budget was not reported")

        # A target that is gone.
        link()
        (repo / "pkg/small.h").unlink()
        r = check()
        if r.returncode == 0 or "pkg/small.h: GONE" not in r.stdout:
            bad.append("a deleted link target was not reported")
        (repo / "pkg/small.h").write_text("v1\n")

        # A write inside a DIRECTORY link is the arrangement, not the failure.
        (wt / "payload/new.dat").write_text("from the worktree\n")
        if not (repo / "payload/new.dat").exists():
            bad.append("the probe did not reproduce the directory link")
        link()
        if check().returncode != 0:
            bad.append("a write inside a directory link was reported; the "
                       "ledger's boundary is wrong")

        # A main checkout has no links to write through.
        if check(tree=repo).returncode != 0:
            bad.append("the check was not a no-op in a main checkout")

        # A worktree linked before the ledger existed must not be failed for it.
        ledger.unlink(missing_ok=True)
        if check().returncode != 0:
            bad.append("a worktree with no ledger was failed rather than told "
                       "how to get one")

        git("worktree", "remove", "--force", str(wt))
    return bad


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quiet", action="store_true",
                    help="print nothing when there is nothing to report")
    ap.add_argument("--self-test", action="store_true",
                    help="prove the ledger against a throwaway worktree and exit")
    args = ap.parse_args()

    if args.self_test:
        failures = self_test()
        for line in failures:
            print(f"  {line}")
        print("self-test: FAILED" if failures else "self-test: every case passes")
        sys.exit(1 if failures else 0)

    located = worktree_and_main()
    if located is None:
        if not args.quiet:
            print("not a worktree: a main checkout holds no links to write through.")
        return
    wt, main = located
    if read_ledger(wt) is None:
        if not args.quiet:
            print(f"no link ledger in {wt}.\n"
                  f"  This worktree was linked before the ledger existed, or by hand.\n"
                  f"  Run `{FIX}` to establish one; until then a write through a\n"
                  "  link here cannot be detected.")
        return

    failures = check_worktree_links()
    if not failures:
        if not args.quiet:
            print("no linked target has moved since this worktree was linked.")
        return
    for line in failures:
        print(line)
    sys.exit(1)


if __name__ == "__main__":
    main()
