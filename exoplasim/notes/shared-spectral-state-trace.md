# Tracing the shared-state race

*Appendix to `shared-spectral-state.md`. In progress, 2026-08-21.*

## What tracing established

A `trace_stamp` at the four phase boundaries of the timestep -- `pre-gpa`,
`post-gpa`, `post-spa`, `post-gpd`, `post-spd` -- one file and one unit per
thread, run twice and diffed.

**First run: `pre-gpa` and `post-gpa` identical, `post-spa` the first
difference.** So `gridpointa` was reproducible and `spectrala` was not.

**Second run, tracing the tendency instead of the state: the first eight
boundaries identical, `post-gpd` of the second iteration the first
difference.**

**The divergence point MOVES between runs.** That is characteristic of a race
rather than of a deterministic mistake, and it is also why tracing does not
converge cheaply: each build tells you where it went wrong that time, not where
it always goes wrong.

## What is ruled out

- The shared-state logic. One thread is deterministic and agrees with the
  one-process build to 1.9e-13.
- The pointers. Correctly associated, checked per thread.
- The publish. After `mpsyncsp` both threads read identical sums, so the
  barrier does what it should.
- The phase boundaries as I understood them. `gridpointa`'s reads are separated
  from `spectrala`'s writes by a barrier, and `gridpointd`'s from
  `spectrald`'s, and a barrier is not passed until every thread arrives -- so
  neither phase can overlap the next.

## What that leaves

Something reads a shared array, or writes one, outside the four phases the
barriers separate. The readers found so far are all accounted for:
`sp3fc` and `dv2uv` in `gridpointa`, `outaccu` after `spectrald`'s publish,
and the root-only output paths. One of them is not where it appears to be, or
there is a reader not yet found.

## The honest position

This is the plan's stated stopping condition approaching: *"the threaded build
cannot be made deterministic. A race on shared state that cannot be found is a
reason to abandon, not to tune."* It is not reached -- the search is bounded and
the detector works -- but the cost per cycle is high and the remaining
hypotheses are thin.

What would change that is a tool that models OpenMP synchronisation.
ThreadSanitizer would be exactly right and is unusable here because libgomp is
not instrumented; an OpenMP-aware race detector, or a libgomp built with TSan
annotations, would turn this from a search into a lookup. That is the thing to
get before spending more cycles on hypothesis-and-rebuild.
