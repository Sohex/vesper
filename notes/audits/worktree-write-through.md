# A worktree regenerating a linked file writes into the main checkout

Audited 2026-08-31. Vesper is a simulated world; the tree under audit is this
project's own tooling, and the checkouts, binaries and agents below are real.

`notes/audits/worktree-stranded-payload.md` audited one direction of
`scripts/link_worktree.py`'s per-file linking: a file a worktree CREATES in a
mixed-tracking directory lives only there and dies with the worktree. This
audits the other direction. A per-file link is a symlink INTO the main checkout,
so a worktree that REGENERATES a linked file writes straight through, and every
artifact the main checkout built from those bytes is invalidated.

The question: what can be written through, what already stops it, and is the
existing backstop the whole answer?

## What is true

**The mechanism is the same per-file linking, and the two directions differ in
which file the worktree touches.** A directory holding tracked content beside
ignored payload cannot be one directory symlink, so the linker makes a real
directory of one link per EXISTING file. Write a NEW file there and it is real
in the worktree alone; write an EXISTING one and the write lands in the main
checkout.

**The write-through direction is the more expensive of the two.** A file that
dies with its worktree costs the work that made it, and the loss is confined to
the tree that took it. A file that writes through mutates a tree other agents
are using, at a moment nobody chose, and the damage lands on an artifact -- a
compiled binary, a staged surface field -- rather than on the file itself. Five
agents were live on this host on 2026-08-31.

**Measured 2026-08-31 on `vendor/lpj-guess/framework/vesper.h`.** It is
generated, gitignored, and sits beside tracked source, so a worktree gets a
per-file link to the main checkout's copy. An agent regenerated it in a worktree
-- correctly, to move a generation timestamp out of a header whose hash gates
the binary -- and the write landed in the main checkout, whose LPJ-GUESS binary
had been built against the previous bytes. `build_lpj_guess.py --verify` began
refusing a binary nobody had touched. One rebuild recovered it, and only because
the agent reported the write; the refusal alone names the symptom and not the
cause.

## What else has this shape

Every per-file link, grouped, measured 2026-08-31 against the standing payload.
The held-back classes are excluded, so this is what a worktree actually gets.

| group | files | bytes | can a worktree regenerate it? |
| --- | --- | --- | --- |
| `vendor/cgenie/**` `.o` `.mod` `.dep` `.a` | 405 | 31 MB | YES, and it is build output of tracked source. See below |
| `hydrography/data/reference/`, `earth_validation/` | 65 | 4324 MB | no: fetched external data |
| `source/<build>/<grid>/` regrids | 56 | 6772 MB | YES, and `source/` is read-only by rule 7 |
| `exoplasim/inputs/<rung>/` staged `.sra` | 34 | 35 MB | YES, and a door already refuses it |
| `hydrography/data/<build>/` products | 12 | 249 MB | YES |
| `exoplasim/analysis/climatology/` | 8 | 88 MB | YES, and downstream components read it |
| `exoplasim/analysis/` ladder and prerebuild records | 6 | 45 MB | YES |
| `.beads/` tracker state | 4 | small | YES, and every `bd` command in every tree does, by design |
| `pedology/data/reference/`, `references/` | 2 | 52 MB | no: fetched external data |
| `aeolian/analysis/dust_baseline.nc` | 1 | small | YES |
| `vendor/lpj-guess/framework/vesper.h` | 1 | small | YES. The measured case |

594 files and 11.6 GB in total, against 321 GB behind the 125 wholly-ignored
directory links. Writes inside a directory link also land in the main checkout,
and that is the ARRANGEMENT rather than the failure: a wholly-ignored directory
is one symlink precisely so a worktree's new climate run survives the worktree.

**`vendor/cgenie` was being linked and should never have been.** Every ignored
path under it is build output -- 234 `.o`, 81 `.dep`, 74 `.mod`, 16 `.a`, and no
ignored directory at all -- so a worktree that built the ocean model wrote its
objects and libraries straight over the main checkout's, and an edit to cgenie
source in a worktree read as a no-op behind the main checkout's objects. That is
exactly CLAUDE.md rule 4's failure mode, and exactly what `SKIP_COMPILED`
already held back for `vendor/exoplasim` and the LPJ-GUESS build. It is in that
tuple as of this audit, which is one prefix because there is no ignored
directory under the subtree to strand.

## What already protects, and what it does not reach

**Write doors refuse.** `lib/write_door.py` holds the one refusal and
`exoplasim/scripts/sra.py:write_sra` was the first to call it, covering
`exoplasim/inputs/<rung>/` for all four staged-field builders at once. That is
the right shape of repair: a skip in the linker would be wrong there, because
the READ is legitimate and a worktree that cannot read the staged fields cannot
stage a run. Which groups have a door and which need none is settled at the end
of this audit.

**One artifact class has a backstop.** `build_lpj_guess.py --verify` compares
the binary against the hash of the source it was built from, and it is what
caught `vesper.h`. It reaches exactly one artifact. Nothing hashes a staged
field, a climatology, a `source/` regrid or a cgenie object into a refusal of
the same kind, and the backstop fires at a LATER agent's build rather than at
the agent that wrote.

**`link_worktree.py --check` did not reach it.** It reported payload living only
in the worktree, which is the other direction, and a write-through leaves every
count complete: nothing is missing, nothing is stranded, nothing is misdirected.

## Why the two repairs the row proposed are not available

**A read-only link cannot be built.** A symlink carries no permissions of its
own; the only mode that decides whether a write lands is the TARGET's, in the
main checkout. Setting it would be a change to shared state -- it would refuse
the main checkout's own legitimate regeneration and every other worktree's read
path at the same time -- and it is per-file rather than per-worktree, so no
worktree can hold it without holding it for all of them. The enforceable
equivalent is a refusal at the write door, which is what `write_sra` already
does, and it needs one door per writer.

**A copy-on-write local shadow costs about a terabyte.** The per-file payload is
11.6 GB, and 84 worktrees were standing on 2026-08-30. It is also the wrong
answer on its own terms twice over: a copy is pinned to the moment it was taken,
so a worktree reads the earliest available input rather than the best, which is
the failure `docs/src/pipeline/loops.md` argues against; and it converts a
write-through into the direction `world-bga1` already closed, silently, with the
worktree's correct regeneration dying at removal after having looked right.

**Doing nothing is not enough, though the backstop did work.** The binary gate
covers one artifact class of the eight in the table above that a worktree can
regenerate, and it reports the symptom to whoever builds next rather than the
cause to whoever wrote. The convention that fills the gap is agents remembering
it: two declined to restage `.sra` fields from a worktree on 2026-08-31 having
reasoned it out unprompted, which is the thing a guard replaces.

## The repair

**A link ledger, and a check at the moment an agent commits.**
`scripts/link_worktree.py` records what every per-file link pointed at -- size
and mtime, plus a content digest below a hash budget -- into the worktree's own
`.git/worktrees/<name>/link-ledger.json`, which no ignore rule has to cover and
which `git worktree remove` deletes along with the tree it describes. Any target
whose bytes have moved since is reported on every run and fails `--check`.
`scripts/check_worktree_links.py` asks the same question and is where the write
reaches the agent that made it: an agent regenerating a file does not then run
the linker, but it does commit.

**It names the fact and not the culprit, and that is the honest limit.** From
inside a worktree, a write from here and a regeneration in the main checkout are
indistinguishable, and both matter: the first invalidates whatever the main
checkout built from those bytes, the second means this worktree's results came
from bytes that are gone. The disposition is the same either way, and re-running
the linker re-baselines -- which is how a legitimate regeneration is accepted
rather than refused.

**The hash budget is a cost line, not a threshold on the finding.** 516 of the
594 per-file links are under 4 MB and come to 74 MB altogether, which hashes in
a fraction of a second; the other 78 are 11.5 GB. Below the budget the content
comparison is exact, so a rewrite that produced identical bytes is correctly not
a finding. Above it a rewrite reports on mtime alone, and that is the
instrument's resolution rather than a false positive: what is detected is the
write through the link, not the difference in the bytes. The whole check runs in
0.13 s on the standing payload, which is what makes it a per-commit cost.

**`.beads/` is linked and deliberately outside the ledger**, because every `bd`
command in every tree writes it through on purpose; ledgering it would report
the arrangement as a fault on every check and drown the finding.

**Exercised both ways.** `check_worktree_links.py --self-test` builds a
throwaway repository of the same shape -- a generated ignored file beside
tracked source, a wholly-ignored directory, and files either side of the hash
budget -- adds a REAL git worktree, and answers a case for each claim
below. A genuine write
through the link must be caught by both doors; the regeneration that caused it
must have succeeded and must survive; re-linking must re-baseline, so that an
accepted regeneration is possible; an identical-bytes rewrite under the budget
must NOT report; a rewrite above it must; a deleted target must; a write inside
a directory link must not; a main checkout and a worktree with no ledger must
both be no-ops. Two mutations of the implementation were confirmed to fail it:
blinding the comparison, and stopping the linker from re-baselining.

## Which groups have a door, and which need none

`lib/write_door.py` is the one refusal, and `exoplasim/scripts/sra.py` states
only what is particular to a staged field on top of it. A door is worth putting
on a group when it has a single Python writer, the write is a regeneration a
worktree plausibly runs, and nothing else already refuses. Three groups meet
that and have one; three need no door, for three different reasons; three want
one and it has to be installed by the component that owns the writer.

| group | disposition |
| --- | --- |
| `exoplasim/inputs/<rung>/` staged `.sra` | door, on `write_sra`, covering all four builders |
| `vendor/lpj-guess/framework/vesper.h` | door, on `build_vesper_header.py`'s install step |
| `aeolian/analysis/dust_baseline.nc` | door, on `build_dust.py`, refusing the report and the field together |
| `vendor/cgenie/**` build objects | NO DOOR: nothing to refuse |
| `source/<build>/<grid>/` regrids | NO SYMLINK DOOR: the wrong instrument |
| `exoplasim/analysis/` prerebuild record | NO DOOR: there is no generator |
| `exoplasim/analysis/climatology/` | door wanted, two writers, world-cgct |
| `exoplasim/analysis/` ladder | door wanted, two writers, world-cgct |
| `hydrography/data/<build>/` | door wanted, six writers, world-pcnd |

**`vendor/cgenie` needs no door because the link is gone.** The subtree is in
`link_worktree.py`'s `SKIP_COMPILED`, so a worktree linked after this audit gets
no link there to write through, and a build in place stays in the worktree.
Deleting the exposure is a better repair than guarding it, and it is the only
one available here anyway: the writer is `make`, and both cgenie drivers already
build from a `git archive` export outside the repository for reasons of their
own. The worktrees standing before the change keep their links, and the ledger
is what reaches those.

**`source/<build>/<grid>` needs no symlink door because a symlink is not what is
wrong with that write.** `source/` is read-only by rule 7: overwriting an
existing build is a violation in the main checkout exactly as much as through a
link, so a guard that fires only in a worktree would license the same write in
the tree where it does most harm. The writer is also the vendored Node exporter,
`vendor/orogen/tools/export-planet.mjs`, which takes its output directory from
`--out` and has no Python door to hang a refusal on. The refusal that belongs
here is by EXISTENCE and not by topology, and it is world-811b.

**The prerebuild record needs no door because nothing generates it.** The
`.x` files under `exoplasim/analysis/prerebuild_binaries_<date>/` were copied by
hand before a rebuild, to keep an A/B reference arm; no script writes them and
`rebuild_binaries.py` never touches `analysis/`. A door is a guard at a write
door, and there is no door. The ledger reports a write through those links, and
for a group whose whole purpose is to hold bytes nobody regenerates, that is the
complete answer rather than a gap.

## What is still carried by convention

The three groups with a door prevent the write; the rest are reported by
`scripts/check_worktree_links.py` after it lands, which names the fact and not
the culprit. That is the honest limit: from inside a worktree, a write from here
and a regeneration in the main checkout are indistinguishable, and the
disposition is the same either way.
