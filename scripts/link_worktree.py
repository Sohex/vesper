#!/usr/bin/env python3
"""Link a worktree's ignored payloads back to the main checkout.

    python scripts/link_worktree.py                    # from inside a worktree
    python scripts/link_worktree.py --worktree PATH    # from the main checkout
    python scripts/link_worktree.py --dry-run
    python scripts/link_worktree.py --check            # report only, exit 1 if incomplete

A git worktree carries the tracked tree and nothing else. Everything this
project deliberately keeps out of history -- the Orogen export payloads, the
reference PDFs and bulk datasets, the climate run output, the Earth validation
caches, `node_modules`, `.venv` -- is simply absent, so a script that runs in the
main checkout dies on a missing file in a worktree. The fix has always been to
symlink them in by hand, and by hand it drifted: of the three worktrees standing
when this was written, one had the T42 export and no interpreter, one had every
export and an interpreter, and one had nothing. The one with everything also had
`exoplasim/runs/runs -> exoplasim/runs`, a link one level too deep, so every run
id under it resolved to nothing.

**The set of things to link is derived, not listed.** `git ls-files --others
--ignored --exclude-standard --directory` in the main checkout is exactly the
"present here, invisible there" set, and it is re-derived on every run, so a new
build or a new cache is picked up without editing this file. It also collapses a
wholly-ignored directory into one entry while listing files individually where
the parent holds tracked content -- which is what stops `exoplasim/runs/` from
being linked as a unit. `runs/INDEX.json` is tracked, so git lists the run
directories one by one and the link lands at `runs/<id>` where a reader expects
it. That is the deep-link bug above, fixed by construction rather than by care.

Three classes are held back, each because linking it is worse than not having
it:

- `.claude/` holds the worktrees themselves, and a worktree containing a link to
  the directory it lives in is a loop.
- `__pycache__`, `docs/book/`, `maps/build/` and `biosphere/generated/` are
  cheap regenerable output. Linking them means a worktree's build writes over
  the main checkout's. `biosphere/generated/` joined the list the moment it
  gained an ignore rule: ten gate scripts write it, several batches run them at
  once, and a shared link would have them overwriting each other's reports in
  the main checkout. Nothing reads it across a component boundary -- the record
  each gate argues lives in `biosphere/notes/`, which is tracked.
- Everything compiled from tracked source that a worktree may have edited:
  `vendor/exoplasim/` and `vendor/lpj-guess/build/`. This is rule 4's failure
  mode with the safety off. A linked binary directory means the worktree runs
  the model the MAIN checkout compiled, so an edit under `vendor/exoplasim`
  looks like a no-op; and a rebuild inside the worktree writes its executables
  over the main checkout's, which is the arm of any A/B this project is running.
  `--model-binaries` links them anyway, for a worktree that does not touch the
  model. `vendor/orogen/node_modules` is NOT in this class -- it is an install,
  not a build of tracked source -- and is linked.

`.venv` IS linked, and it carries a caveat this script prints rather than
solves: ExoPlaSim is installed editable from the main checkout's
`vendor/exoplasim`, so `import exoplasim` in a worktree reads the main
checkout's model source no matter which tree the interpreter was invoked from.
That is a property of the editable install, not of the link.

Idempotent, and safe to re-run after a build is added or archived: a correct
link is left alone, a link that points somewhere stale is repaired, a link into
the main checkout whose target no longer exists is removed, and a real file or
directory sitting where a link should go is reported and never touched. Links
that point at a different path in the main checkout than the one they occupy --
the deep-link bug above -- are reported too, and left for a human.

It also reports IGNORED PAYLOAD THAT LIVES ONLY IN THE WORKTREE, and `--check`
fails on it. Per-file linking has a consequence that catches people: a directory
holding tracked content beside ignored payload cannot be one symlink, so what
lands there is a real directory of per-file links, and a file CREATED there
afterwards is real, in this worktree alone, and ignored -- so it is never
committed and it is gone when the worktree is removed. `references/` lost
thirteen PDFs that way while their `INDEX.md` rows, being tracked, survived and
outlived their own artifacts; `exoplasim/runs/` and `source/` have the same
shape, where the thing that dies is a climate run or an export payload. Reported
and never moved: whether it belongs in the main checkout is a judgement.

IT ALSO KEEPS A LINK LEDGER, and reports the OTHER direction of that same shape.
A per-file link is a symlink INTO the main checkout, so a worktree that
REGENERATES a linked file writes straight through, and every artifact the main
checkout built from those bytes is invalidated at a moment nobody chose. That is
worse than losing a file: the damage lands on a compiled binary or a staged
field in a tree other agents are using, rather than on the file itself.
`vendor/lpj-guess/framework/vesper.h` did exactly that on 2026-08-31. So every
run records what each per-file link pointed at, into the worktree's own
`.git/worktrees/<name>/link-ledger.json`, and reports any target whose bytes
have moved since. It names the FACT and not the culprit, because from inside a
worktree a write from here and a regeneration in the main checkout are
indistinguishable and both matter. `scripts/check_worktree_links.py` asks the
same question at the moment an agent commits.

The last thing it does is check `git status` in the worktree. A symlink is a
file, not a directory, so an ignore rule ending in `/` does not match the link
that stands in for the directory it named -- which is why `.gitignore` carries a
second set of patterns for exactly these paths. If a link shows up as untracked,
this prints the pattern that would cover it and exits non-zero, because the
alternative is a worktree committing absolute paths that mean nothing anywhere
else.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

# Held back because the worktree lives inside it.
SKIP_SELF = (".claude/",)

# Held back because they are regenerable output, and a link makes the worktree
# write over the main checkout's copy.
SKIP_REGENERABLE = ("docs/book/", "maps/build/", "biosphere/generated/")

# Held back because they are compiled from tracked source the worktree may have
# edited. See the module docstring and CLAUDE.md rule 4.
SKIP_COMPILED = ("vendor/exoplasim/", "vendor/lpj-guess/build/")

# LINKED, and deliberately left out of the link ledger below. `.beads/` is the
# issue tracker's shared state and every `bd` command in every tree writes it
# through the link on purpose, so ledgering it would report the arrangement as
# a fault on every single check and drown the finding the ledger exists to
# make.
LEDGER_SKIP = (".beads/",)

# The ledger's hash budget: content is compared EXACTLY below it and by size
# and mtime above it. Measured 2026-08-31 against the standing payload -- 516
# of the 594 per-file links are under this and come to 74 MB altogether, which
# hashes in well under a second, while the other 78 are 11.5 GB and hashing
# them on every commit would cost more than the answer is worth. Above the
# budget a rewrite that happened to produce identical bytes still reports, and
# that is the instrument's honest resolution rather than a false positive: what
# is being detected is the WRITE THROUGH THE LINK, not the difference in the
# bytes.
LEDGER_HASH_BUDGET = 4 * 1024 * 1024

LEDGER_NAME = "link-ledger.json"


def run(args: list[str], cwd: Path) -> str:
    return subprocess.run(args, cwd=cwd, check=True, capture_output=True,
                          text=True).stdout


def locate(explicit: Path | None) -> tuple[Path, Path]:
    """Return (main checkout, worktree), refusing anything that is not both."""
    here = Path.cwd()
    common = Path(run(["git", "rev-parse", "--path-format=absolute",
                       "--git-common-dir"], here).strip())
    main = common.parent
    if explicit is not None:
        wt = explicit.resolve()
        if not (wt / ".git").exists():
            raise SystemExit(f"{wt} is not a git worktree (no .git)")
    else:
        wt = Path(run(["git", "rev-parse", "--path-format=absolute",
                       "--show-toplevel"], here).strip())
    if wt == main:
        raise SystemExit(
            f"{wt} IS the main checkout, and linking it to itself is not a\n"
            "thing to do. Run this from inside a worktree, or pass --worktree.")
    # A worktree of a DIFFERENT repository would link the wrong payloads in.
    other = Path(run(["git", "rev-parse", "--path-format=absolute",
                      "--git-common-dir"], wt).strip()).parent
    if other != main:
        raise SystemExit(f"{wt} is a worktree of {other}, not of {main}")
    return main, wt


def ignored_entries(main: Path, model_binaries: bool) -> list[str]:
    """The present-here, invisible-there set, minus the three held-back classes.

    Entries ending in "/" are wholly-ignored directories; the rest are files git
    listed individually because their parent holds tracked content.
    """
    out = run(["git", "ls-files", "--others", "--ignored", "--exclude-standard",
               "--directory"], main)
    skip = SKIP_SELF + SKIP_REGENERABLE
    if not model_binaries:
        skip += SKIP_COMPILED

    kept = []
    for line in out.splitlines():
        entry = line.strip()
        if not entry:
            continue
        if any(entry.startswith(p) for p in skip):
            continue
        if "__pycache__" in Path(entry).parts:
            continue
        kept.append(entry)

    # git lists both a directory and, in some ignore-rule shapes, its contents.
    # Dropping the contained entries matters for more than tidiness: linking the
    # directory first would put every later dst INSIDE the main checkout, so the
    # worktree would write there.
    dirs = sorted(e for e in kept if e.endswith("/"))
    return sorted(e for e in kept
                  if not any(e != d and e.startswith(d) for d in dirs))


def ancestor_is_link(dst: Path, wt: Path) -> Path | None:
    """The first symlinked ancestor of dst below the worktree root, if any.

    Walked from the root downwards and stopped at the root, so a symlink on the
    path to the worktree itself is not mistaken for one inside it.
    """
    rel = dst.relative_to(wt)
    walk = wt
    for part in rel.parts[:-1]:
        walk = walk / part
        if walk.is_symlink():
            return walk
    return None


def link_one(entry: str, main: Path, wt: Path, dry: bool) -> tuple[str, str]:
    """Return (state, detail) for one entry. States are the report's vocabulary."""
    rel = entry.rstrip("/")
    src = main / rel
    dst = wt / rel

    if not src.exists():
        return "vanished", f"{rel} is gone from the main checkout"

    caught = ancestor_is_link(dst, wt)
    if caught is not None:
        return "shadowed", (f"{rel} sits under the link "
                            f"{caught.relative_to(wt)}; writing it would write "
                            "into the main checkout")

    if dst.is_symlink():
        target = Path(os.readlink(dst))
        if target == src:
            return "ok", rel
        if not dry:
            dst.unlink()
            dst.symlink_to(src)
        return "repaired", f"{rel} pointed at {target}"

    if dst.exists():
        return "conflict", f"{rel} exists as a real {'directory' if dst.is_dir() else 'file'}"

    if not dry:
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.symlink_to(src)
    return "linked", rel


def sweep(wt: Path, main: Path, dry: bool) -> tuple[list[str], list[str]]:
    """Find the links into the main checkout that are wrong, two ways.

    DANGLING, and removed: the target is gone. An archived build leaves one of
    these at every path it exported to, and a dangling link reads as a
    present-but-broken artifact rather than an absent one, which is the more
    expensive of the two to diagnose.

    MISDIRECTED, and only reported: the link is at `X` in the worktree and
    points at `Y` in the main checkout, where `X != Y`. Every link this script
    makes mirrors its own path, so a mismatch is a hand-made link that landed
    somewhere else -- `exoplasim/runs/runs -> exoplasim/runs`, one level too
    deep, is the one that was standing when this was written, and it made every
    run id under it resolve to nothing. Reported rather than removed, because
    the target exists and a deliberate link elsewhere is somebody's decision.
    """
    removed, misdirected = [], []
    # os.walk with followlinks=False, and symlinked directories taken out of the
    # descent by hand: a linked export is gigabytes of the main checkout, and
    # walking into it costs minutes to find nothing. They are still INSPECTED --
    # a dangling one is exactly what this is looking for -- just not entered.
    for root, dirs, files in os.walk(wt, followlinks=False):
        rootp = Path(root)
        dirs[:] = [d for d in dirs if d not in (".git", ".claude")]
        linked_dirs = [d for d in dirs if (rootp / d).is_symlink()]
        dirs[:] = [d for d in dirs if d not in linked_dirs]
        for name in linked_dirs + files:
            path = rootp / name
            if not path.is_symlink():
                continue
            target = Path(os.readlink(path))
            if not target.is_absolute() or main not in target.parents:
                continue
            rel = path.relative_to(wt)
            if not target.exists():
                removed.append(str(rel))
                if not dry:
                    path.unlink()
            elif target.relative_to(main) != rel:
                misdirected.append(f"{rel} -> {target.relative_to(main)}")
    return removed, misdirected


def stranded(wt: Path, main: Path, model_binaries: bool) -> list[str]:
    """Ignored payload that exists ONLY in the worktree, and dies with it.

    A wholly-ignored directory is linked as one symlink, so a write inside it
    lands in the main checkout. A directory holding tracked content beside
    ignored payload CANNOT be: `references/INDEX.md` and
    `exoplasim/runs/INDEX.json` are tracked, so git lists their payload file by
    file and this script makes a REAL directory holding one link per EXISTING
    file. A file created there afterwards is a real file in the worktree alone.
    It is ignored, so it is never committed, and when the worktree is removed it
    is gone -- while a tracked row describing it survives and outlives its own
    artifact. That is rule 7's hazard, not a bookkeeping one: a reference PDF, a
    climate run under `exoplasim/runs/`, an export payload under `source/`.

    Reported and never moved. Where it belongs is a judgement -- the main
    checkout, or nowhere -- and this script does not make it.

    The three held-back classes are excluded, because a worktree building its
    own `vendor/exoplasim` is the arrangement, not the failure. A directory the
    main checkout does not have at all is reported whole and not descended: a
    run directory is thousands of files and one line is the finding.
    """
    skip = SKIP_SELF + SKIP_REGENERABLE
    if not model_binaries:
        skip += SKIP_COMPILED

    candidates: list[str] = []
    for root, dirs, files in os.walk(wt, followlinks=False):
        rel_root = Path(root).relative_to(wt)
        dirs[:] = sorted(d for d in dirs
                         if d not in (".git", ".claude", "__pycache__")
                         and not (Path(root) / d).is_symlink())
        for name in list(dirs):
            rel = (rel_root / name).as_posix()
            if any(f"{rel}/".startswith(p) for p in skip):
                dirs.remove(name)
            elif not (main / rel).exists():
                candidates.append(f"{rel}/")
                dirs.remove(name)
        for name in files:
            path = Path(root) / name
            if path.is_symlink():
                continue
            rel = (rel_root / name).as_posix()
            if any(rel.startswith(p) for p in skip):
                continue
            if not (main / rel).exists():
                candidates.append(rel)

    if not candidates:
        return []
    # One call, and check-ignore is the authority: a candidate git would happily
    # commit is ordinary new work and not this failure at all.
    out = subprocess.run(["git", "check-ignore", "--stdin"], cwd=wt, text=True,
                         input="\n".join(candidates), capture_output=True)
    return sorted(line.strip() for line in out.stdout.splitlines() if line.strip())


def ledger_path(wt: Path) -> Path:
    """Where the link ledger lives: the worktree's PRIVATE git directory.

    `.git/worktrees/<name>/` belongs to one worktree, is not a path git will
    ever offer to commit, and `git worktree remove` deletes it along with the
    tree it describes. So the ledger needs no ignore rule and cannot outlive
    its subject -- which matters, because a ledger that survived its worktree
    would go on to accuse the next one.
    """
    git_dir = run(["git", "rev-parse", "--path-format=absolute", "--git-dir"], wt)
    return Path(git_dir.strip()) / LEDGER_NAME


def stamp(src: Path) -> dict:
    """What the ledger records for one link target: size, mtime, and content.

    The digest is present only below `LEDGER_HASH_BUDGET`. Its absence is what
    tells the comparison it may not claim more than "this was rewritten".
    """
    st = src.stat()
    entry = {"size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if st.st_size <= LEDGER_HASH_BUDGET:
        digest = hashlib.blake2b(digest_size=16)
        with src.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        entry["digest"] = digest.hexdigest()
    return entry


def ledgered(rels: list[str], model_binaries: bool) -> list[str]:
    """The per-file links the ledger covers.

    DIRECTORY links are out of it, and the boundary is the point rather than an
    omission. A wholly-ignored directory is linked as ONE symlink precisely so
    that a worktree's write inside it lands in the main checkout and survives
    the worktree -- that is the arrangement, not the failure, and the payload
    under those links comes to hundreds of gigabytes. What the ledger covers is
    the other shape: an EXISTING file, linked individually because its parent
    holds tracked content, which a worktree regenerates in place and thereby
    replaces under everything the main checkout built from it.
    """
    skip = LEDGER_SKIP + SKIP_SELF + SKIP_REGENERABLE
    if not model_binaries:
        skip += SKIP_COMPILED
    return [r for r in rels if not any(r.startswith(pref) for pref in skip)]


def read_ledger(wt: Path) -> dict | None:
    """The ledger, or None if this worktree has never had one written."""
    path = ledger_path(wt)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def write_ledger(wt: Path, main: Path, rels: list[str]) -> None:
    """Record what every per-file link pointed at, as of now."""
    entries: dict[str, dict] = {}
    for rel in rels:
        src = main / rel
        try:
            if src.is_file() and not src.is_symlink():
                entries[rel] = stamp(src)
        except OSError:
            continue
    payload = {
        "main": str(main),
        "linked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "hash_budget": LEDGER_HASH_BUDGET,
        "entries": entries,
    }
    ledger_path(wt).write_text(json.dumps(payload, indent=1, sort_keys=True),
                              encoding="utf-8")


def changed_since_link(wt: Path, main: Path) -> tuple[list[str], str | None]:
    """Link targets in the main checkout whose bytes have moved since linking.

    THE WRITE-THROUGH DIRECTION. `stranded()` above catches a worktree's file
    dying with the worktree. This catches the opposite and sharper one: a
    per-file link is a symlink INTO THE MAIN CHECKOUT, so a worktree that
    REGENERATES one writes through, and every artifact the main checkout built
    from it is invalidated at a moment nobody chose.
    `vendor/lpj-guess/framework/vesper.h` did exactly that on 2026-08-31 -- an
    agent moved a generation timestamp out of it, correctly, and the write
    landed in the main checkout whose LPJ binary had been built against the
    previous bytes.

    IT REPORTS THE FACT AND NOT THE CULPRIT, because from inside a worktree the
    two readings are indistinguishable and BOTH matter. Either this worktree
    wrote through -- and whatever the main checkout built from those bytes is
    now stale -- or the main checkout regenerated legitimately, and this
    worktree's own results were computed from bytes that no longer exist. The
    disposition is the same either way: settle which it was, rebuild what
    depended on it, and re-run this script to re-baseline.

    Returns the findings and the timestamp they are measured against.
    """
    ledger = read_ledger(wt)
    if ledger is None:
        return [], None

    findings = []
    for rel, was in sorted(ledger.get("entries", {}).items()):
        src = main / rel
        if not src.exists():
            findings.append(f"{rel}: GONE from the main checkout")
            continue
        try:
            now = src.stat()
        except OSError as exc:
            findings.append(f"{rel}: unreadable in the main checkout ({exc})")
            continue
        if now.st_size != was["size"]:
            findings.append(f"{rel}: CONTENT CHANGED, "
                            f"{was['size']} -> {now.st_size} bytes")
        elif "digest" in was and now.st_size <= LEDGER_HASH_BUDGET:
            if stamp(src)["digest"] != was["digest"]:
                findings.append(f"{rel}: CONTENT CHANGED, same size")
        elif now.st_mtime_ns != was["mtime_ns"]:
            budget = LEDGER_HASH_BUDGET // (1 << 20)
            findings.append(f"{rel}: REWRITTEN, same size, content not compared "
                            f"(over the {budget} MB hash budget)")
    return findings, ledger.get("linked_at")


def check_git_clean(wt: Path, created: list[str]) -> list[str]:
    """Links that git can see. Each one is a missing .gitignore pattern.

    An ignore rule ending in "/" matches a directory and not the symlink that
    stands in for it, so the rules that hide a payload in the main checkout do
    not necessarily hide its link here.
    """
    out = run(["git", "status", "--porcelain", "--untracked-files=all"], wt)
    visible = {line[3:].strip().strip('"') for line in out.splitlines()
               if line[:2] in ("??", "A ", "AM")}
    return sorted(c for c in created if c in visible)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Link a worktree's ignored payloads back to the main checkout.")
    ap.add_argument("--worktree", type=Path,
                    help="the worktree to link; default is the one you are in")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change and touch nothing")
    ap.add_argument("--check", action="store_true",
                    help="like --dry-run, but exit 1 if anything is missing")
    ap.add_argument("--model-binaries", action="store_true",
                    help="also link vendor/exoplasim and the LPJ-GUESS build "
                         "(see the warning this prints)")
    ap.add_argument("--verbose", "-v", action="store_true",
                    help="name every path instead of counting them")
    args = ap.parse_args()

    dry = args.dry_run or args.check
    main_root, wt = locate(args.worktree)
    print(f"main checkout: {main_root}")
    print(f"worktree:      {wt}")
    if dry:
        print("DRY RUN -- nothing is written\n")
    else:
        print()

    states: dict[str, list[str]] = {}
    per_file: list[str] = []
    for entry in ignored_entries(main_root, args.model_binaries):
        state, detail = link_one(entry, main_root, wt, dry)
        states.setdefault(state, []).append(detail)
        # The ledger's subject is the per-file links, and only the ones that
        # are actually standing: a "would link" in a dry run is not a link.
        if not entry.endswith("/") and state in ("ok", "linked", "repaired"):
            per_file.append(entry)

    for state in ("linked", "repaired", "ok", "conflict", "shadowed", "vanished"):
        items = states.get(state, [])
        if not items:
            continue
        # A tree at rest is almost all "ok"; the states that need a human are the
        # ones worth spelling out whether or not --verbose was asked for.
        spell_out = args.verbose or state in ("conflict", "shadowed", "vanished")
        head = {"linked": "would link" if dry else "linked",
                "repaired": "would repair" if dry else "repaired",
                "ok": "already linked",
                "conflict": "NOT LINKED, something real is in the way",
                "shadowed": "NOT LINKED, the parent is already a link",
                "vanished": "SKIPPED, no longer in the main checkout"}[state]
        print(f"{head}: {len(items)}")
        if spell_out:
            for item in items:
                print(f"    {item}")

    removed, misdirected = sweep(wt, main_root, dry)
    if removed:
        print(f"{'would remove' if dry else 'removed'} {len(removed)} dangling "
              f"link(s) into the main checkout:")
        for r in removed:
            print(f"    {r}")
    if misdirected:
        print(f"\n{len(misdirected)} link(s) point at a DIFFERENT path in the main\n"
              "checkout than the one they sit at, so a reader following them lands\n"
              "somewhere other than where the layout says. Left alone; delete each by\n"
              "hand unless it was deliberate:")
        for m in misdirected:
            print(f"    {m}")

    orphans = stranded(wt, main_root, args.model_binaries)
    if orphans:
        print(f"\n{len(orphans)} ignored path(s) exist ONLY IN THIS WORKTREE and will be\n"
              "gone when it is removed. A directory that holds tracked content beside\n"
              "ignored payload is linked FILE BY FILE, so anything written there\n"
              "afterwards is a real file here and nowhere else -- and being ignored, it\n"
              "is never committed either. Copy each to the main checkout if it is worth\n"
              "keeping, or delete it:")
        for o in orphans:
            print(f"    {o}")

    ledgered_rels = ledgered(per_file, args.model_binaries)
    changed, linked_at = changed_since_link(wt, main_root)
    if changed:
        print(f"\n{len(changed)} linked target(s) in the main checkout have CHANGED since\n"
              f"this worktree was linked ({linked_at}). A per-file link is a symlink INTO\n"
              "the main checkout, so either something here regenerated one and wrote\n"
              "straight through -- invalidating whatever the main checkout built from it,\n"
              "at a moment nobody chose -- or the main checkout regenerated it and this\n"
              "worktree's own results were computed from bytes that are gone. Settle which,\n"
              "rebuild what depended on those bytes, and re-run this script to re-baseline:")
        for c in changed:
            print(f"    {c}")
    elif linked_at is None and not dry:
        print("\nno link ledger yet; writing one now, so the next run can say whether a\n"
              "linked target was written through.")

    if not dry:
        write_ledger(wt, main_root, ledgered_rels)

    if not args.model_binaries:
        print("\nvendor/exoplasim and vendor/lpj-guess/build were NOT linked: they are\n"
              "compiled from tracked source this worktree may have edited, so a link\n"
              "would both hide the edit and let a rebuild here overwrite the main\n"
              "checkout's binaries. Build them in the worktree, or pass\n"
              "--model-binaries if this worktree does not touch the model.")
    else:
        print("\n--model-binaries: vendor/exoplasim IS linked. A rebuild in this\n"
              "worktree now writes over the main checkout's executables, and an edit\n"
              "to the model source here does not change what runs. Do neither.")

    if (wt / ".venv").is_symlink():
        print("\n.venv is linked, and ExoPlaSim is installed editable from the MAIN\n"
              "checkout's vendor/exoplasim. `import exoplasim` in this worktree reads\n"
              "the main checkout's model source whichever tree you run from.")

    created = [c for c in states.get("linked", []) + states.get("ok", [])]
    if not dry:
        visible = check_git_clean(wt, created)
        if visible:
            print(f"\n{len(visible)} link(s) are VISIBLE TO GIT in the worktree, so an\n"
                  "`git add -A` here would commit absolute paths. Add a pattern for\n"
                  "each to .gitignore in the main checkout -- without a trailing\n"
                  "slash, which matches the directory and not the link:")
            for v in visible:
                print(f"    {v}")
            sys.exit(1)
        print("\ngit status in the worktree is clean of these links.")

    if args.check:
        # MISDIRECTED COUNTS. A link at `X` pointing at `Y` is the deep-link bug
        # in the module docstring -- `exoplasim/runs/runs -> exoplasim/runs`, one
        # level too deep, which made every run id under it resolve to nothing --
        # and it is the failure this script was written for. It is not "missing":
        # no entry is absent, so the payload count came out complete and --check
        # printed the warning and then exited 0 beside it. A worktree carrying
        # one is not correctly linked, whatever the count says.
        #
        # STRANDED PAYLOAD COUNTS TOO, and for the same reason: nothing is
        # absent, so every count above comes out complete while the worktree
        # holds artifacts that die with it. That is the failure this whole
        # arrangement produces silently, and --check is where a caller asks
        # whether the worktree is in a state it can be thrown away from.
        missing = len(states.get("linked", [])) + len(states.get("repaired", []))
        blocked = len(states.get("conflict", [])) + len(states.get("shadowed", []))
        if missing or blocked or misdirected or orphans or changed:
            print(f"\nINCOMPLETE: {missing} to link or repair, {blocked} blocked, "
                  f"{len(misdirected)} pointing somewhere else, "
                  f"{len(orphans)} living only here, "
                  f"{len(changed)} changed since linking.")
            sys.exit(1)
        print("\ncomplete: every ignored payload is linked, none of it lives only "
              "here, and no linked target has moved.")


if __name__ == "__main__":
    main()
