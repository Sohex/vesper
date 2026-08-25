# Who copies an executable, and does any of them establish which one ran

*Worldbuilding frame: a BUILD audit of the Vesper project's climate model on
this desktop. Nothing here is about the simulated planet. Swept 2026-08-25 at
eb8a18d9, statically, over every `most_plasim_*.x` reference in the tree.*

CLAUDE.md rule 4 is about the REGISTRY: a rebuild refreshes only the
configuration you ran, so the rest go stale in silence. This sweep is about the
other side of it. `rebuild_binaries.py --verify` can truthfully report every
installed binary current while a consumer integrates a COPY the manifest does
not track, because a copy carries no provenance of its own: it is a file with a
name, and the name is the only thing tying it back to a source set and a flag
line. `--verify` answers "are the installed binaries current"; it cannot answer
"is this the one that ran".

That is not hypothetical. On 2026-08-24 `rebuild_binaries.py` rebuilt every
configuration from patched radmod source and `--verify` reported all five
current; `stability_probe.py` then integrated a run directory's older copy and
printed the pre-repair surface albedos, which read as the patch having failed to
take.

## The sweep

Every site that selects an executable and then runs it, or copies one somewhere
that will:

| consumer | where its executable comes from | does it establish which one ran |
| --- | --- | --- |
| `exoplasim/scripts/stability_probe.py` | a TEMPLATE RUN DIRECTORY, by glob | no, and its output recorded nothing. THE DEFECT: world-anl |
| `exoplasim/scripts/make_profile_bed.py` | `--binary-dir`, the registry by default | yes: sha against `binary_manifest.json`, refuses, `--allow-unregistered` to proceed |
| `exoplasim/scripts/reproducibility_matrix.py` | the registry, copied into a bed | yes: every bed executable against the manifest before a cell runs |
| `exoplasim/scripts/run_exoplasim.py` | the run directory, which ExoPlaSim filled from the registry at prepare | records the sha in `run_manifest.json`; does not compare it |
| `exoplasim/scripts/continue_exoplasim.py` | the copy sitting in the run directory | compares against the sha the PREPARE recorded and prints a note; never against the manifest |
| `exoplasim/scripts/run_stellar_cycle.py` | the registry, by composed name | records the sha; does not compare it |
| `thread_count_sweep.sh`, `verify_shared_determinism.sh`, `verify_weight_factorisation_model.sh` | the registry, `cp` into a scratch bed | prints the sha; does not compare it |
| `exoplasim/scripts/_bed_guard.sh` | -- | checks a bed's binary against the SURFACE resolution beside it, which is a different question |
| `exoplasim/scripts/sweep_compiler_flags.py` | builds its own arms | n/a: an arm is not a registry entry |
| `exoplasim/scripts/stack_floor.py` | -- | parses the name only |

## Two classes, and only one of them is severe

**A copy of unbounded age.** One consumer took its executable from a run
directory: `stability_probe.py`. A run directory holds the copy that was placed
in it when that run was staged, so its age is the run's age and has no upper
bound. `find_template` picked the most recently MODIFIED matching run directory,
so the copy was not even stable between two invocations of the same command.

**Registry-fresh but unverified.** Everything else resolves against the model's
`plasim/run/`, so its executable is as current as the registry was at the moment
it was taken. `scripts/check_consistency.py` calls `rebuild_binaries.verify()`
and CLAUDE.md rule 8 requires it before an expensive run, so the registry is
gated -- but by a gate a caller has to remember rather than one it cannot skip,
and a run prepared before a rebuild and continued after it is caught only by
`continue_exoplasim.py`'s per-segment note.

## Which recorded probe results came from which binary: they cannot be attributed

`exoplasim/analysis/stability_probe.json` carries no attribution for any entry
written before this sweep, and none can be recovered:

- a probe entry records the rung, the timestep, the filter strength and the
  wall clock, and no executable name, sha or template path;
- there is no per-entry timestamp. The file carries one `generated` field and
  it is rewritten on every write, so it dates the last probe and nothing else;
- `find_template` selected by modification time, so knowing WHEN a probe ran
  would still not name the directory it drew from;
- run directories are untracked, so what any of them held at that time is not
  recoverable. `exoplasim/runs/INDEX.json` records an `executable_sha256` per
  run, but it is a regenerated snapshot of what is on disk now.

The file's git history brackets an entry to the commit that added it, which
dates when it was COMMITTED and not which executable produced it. So every
number in that file is a measurement of an unnamed model, and the refusal
boundaries and per-step costs it holds are re-derivable and not checkable.

## What is checked now

`rebuild_binaries.py` gained `registered(name)` and `unregistered(exe)`: the
manifest lookup a consumer needs, kept beside the code that WRITES the manifest
so a second statement of the rule cannot drift from it. `stability_probe.py`
refuses when the executable it is about to run is not what the manifest
registers under that name, prints both shas, and records the name, the sha and
the build profile in every result. Restaging from the registry is available and
is `--restage`, never automatic: a silent substitution would leave "which
binary" exactly as unanswerable as it was.

`make_profile_bed.py` and `reproducibility_matrix.py` each still carry their own
statement of the same comparison. They are correct, and converging them onto
`rebuild_binaries.unregistered` is bookkeeping rather than a repair.
