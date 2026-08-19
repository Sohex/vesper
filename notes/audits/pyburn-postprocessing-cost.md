# Postprocessing was quadratic in record count, and the profiler hid it

Measured 2026-08-18, T42, 10 layers, one orbit at `NLOWIO = 0`: 2.4 GB of raw
model output, 120 requested codes, 281 variables written.

## Result

Driven with the exact kwargs `Model.run()` forwards for a configured
postprocessor -- `namelist` None and the project's 120-code `variables` list:

| | wall | np.append calls |
| --- | --- | --- |
| pyburn as shipped | 304.0 s | 3,062,749 |
| pyburn fixed | 29.4 s | 1 |

10.3x, with all 123 output variables **bitwise identical**. Both timings are
unprofiled, on the same raw file.

**It reproduces the pipeline's own recorded history.** `continue_exoplasim.py`
records 387 s of wall clock an orbit at `NLOWIO = 0`, measured 2026-08-17. This
bench gives 81 s of model plus 304 s of pyburn, 385 s. After the fix the same
orbit is 81 + 29 = 110 s, so an orbit costs about a third of what it did.

The model takes 77 to 93 s, so postprocessing went from roughly four times the
cost of the science to about a third of it.

**A caution about which invocation is measured.** pyburn ignores `variables`
whenever a `namelist` is also given (`pyburn.py`, "Scrape namelist"). The
pipeline never passes one -- `run()` calls `postprocess(dataname, None, ...)`
and `cfgpostprocessor` stores `namelist=None` -- so production processes exactly
the 120 codes it asks for and never reads `example.nl`. A hand-written harness
that passes `example.nl` instead processes its 216 codes, including `wa` and the
streamfunction that production never requests, and measures 363.5 s against
46.1 s. Same defect, heavier workload; quote the 304 to 29.4 figures.

## The defect that mattered

`readallvariables` accumulated every variable one record at a time:

    variables[kcode] = np.append(variables[kcode],field)

`np.append` reallocates and copies the entire accumulated array on every call,
so reading N records costs O(N^2) in copied bytes. One orbit is about 81,000
records, and this single line was **246 s of the 363 s**. Records are now
collected in a per-code list and joined once with `np.concatenate`, which is the
identical result because `np.append` flattens anyway.

`readvariablecode` carries the same pattern and is fixed the same way. It also
could not have worked as written: `if not variable:` raises on an ndarray of
more than one element. Nothing in this project calls it.

## The defect that did not matter

The omega/`wap`, `wa` and streamfunction branches each ran a triple loop over
(time, lat, lon), calling `np.append` twice and `scipy.integrate.cumulative_trapezoid`
once **per grid cell** -- 1.49e6 iterations an orbit each. That is genuinely bad
code and it is now integrated on the level axis for the whole array at once, but
it was worth only about 81 s of the 363 s.

## How this was nearly misdiagnosed, twice

The first profile said `np.append` was 246 s of 412 s with 6,044,637 calls, and
the only place with millions of `np.append` calls is the triple loops. That
reading was wrong in both directions.

**cProfile charges per call event.** Six million tiny appends inside the loops
cost almost nothing in numpy terms; what the profile was measuring was mostly
its own instrumentation. Vectorising them moved an unprofiled orbit from 363.5 s
to 282.8 s -- real, but not the 5x the profile implied.

**The same 246 s belonged to a different caller.** Re-profiling the vectorised
build showed `np.append` still at 245.8 s but now with **80,861 calls** instead
of 6,044,637. Same cost, one seventy-fifth the calls, which means the expensive
appends were always the large ones in the reader and never the small ones in the
loops. `readallvariables` held 271.5 s of 297.1 s.

The lesson is specific and worth keeping: **a profile attributes time to a
FUNCTION, and `np.append` had two callers with call counts three orders of
magnitude apart and wildly different per-call costs.** The tell was the ratio.
Sorting by `tottime` alone points at the wrong one, and comparing a profiled
baseline against an unprofiled treatment then flatters the wrong fix. Both
numbers here are unprofiled, and the call counts are reported beside them.

## What this changes elsewhere

Two earlier conclusions move.

`continue_exoplasim.py` records that `NLOWIO = 0` costs 104 s against 387 s of
wall clock per orbit "because turning the accumulation off multiplies the raw
volume pyburn has to chew". The volume claim is right and the price is now much
smaller, because the chewing was quadratic rather than proportional.

The measurement that two 8-rank runs beat one 16-rank run by 1.09x was taken
when postprocessing dominated and was single-threaded, so overlapping it was
worth something. With postprocessing at 46 s against a 77 s model, **run at 16
ranks**; the concurrency trick was a workaround for this defect.
