# Working on a large blob: batch, chunk, checkpoint, report

A standing rule, written 2026-08-17 after a 15.3 GB postprocess ran for 73
minutes, reached 27.6 GB resident, wrote nothing, and had to be killed. Nothing
about that was surprising in hindsight and all of it was avoidable.

**Any step whose input is larger than a few GB gets all four of these**, and a
step that has none of them is not finished, however correct its arithmetic.

## Batch and chunk

Size the working set deliberately and declare the number. Do not let it be
whatever the input happens to be.

The failure to avoid is a library that looks like it streams and does not. The
case above: `pyburn.postprocess` was asked for three variables out of a raw file
with roughly 280 field-levels per sample, and its memory came out at 18.9 MB per
sample, which is the whole decoded record. **Restricting the output bought
nothing, because the cost was in the decode.** Measure MB per unit of input
before assuming a filter helps.

Chunk on a boundary the format actually has. Sequential Fortran output is
delimited by record headers and by timestep; a byte offset chosen without
reference to those is a corrupt chunk.

## Checkpoint

Write results as they are produced, record what is done, and make re-running the
same command resume rather than restart. The test is blunt: **kill it at 80% and
run it again.** If that costs more than the last chunk, it is not checkpointed.

Keep the checkpoint next to the output and validate it against the input before
trusting it -- a state file that has lost track of which input it belongs to is
worse than none, because it will silently skip work that was never done.

Do not clean up the input until the output is complete and verified. The 15.3 GB
raw above survived the kill only because ExoPlaSim deletes it on success, and it
is the reason 20 minutes of model time did not have to be repeated.

## Report progress

Print per chunk: which chunk of how many, elapsed, an estimate of what is left,
and resident memory. Flush it.

This is not decoration. The killed job had produced no output and no log line for
73 minutes, so there was no way to tell a slow job from a hung one from a job
that was going to exhaust memory -- and the answer to those three is different.
Resident memory in particular is what turns "it is still going" into "it will not
finish", and it costs one line to read from `/proc/self/status`.

## Before starting, not after

Estimate the footprint against a measured per-unit cost and compare it with what
the machine has. The arithmetic is short and it is the difference between a
scheduling decision and an OOM: the same 18.9 MB per sample predicted that a
regular ExoPlaSim orbit's postprocess needs about 7 GB, which is what decided
whether a second run could start alongside.

Whatever the estimate says, leave the headroom. Linux picks its OOM victim by
resident size, so the process that dies is the biggest one, which is usually the
one with the most work invested in it and rarely the one that caused the
pressure.
