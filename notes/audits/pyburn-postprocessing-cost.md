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

# And then every value in the file was boxed into a Python float

Measured 2026-08-19, on a SYNTHETIC raw file, because no raw file survives
postprocessing and there was none left to re-measure the real one against. The
synthetic matches the geometry above: T42, 10 levels, 77 single-level codes and
36 ten-level codes over 185 output steps, 80,845 records of 128x64 float32,
2.65 GB. It reads correctly through unmodified pyburn and refactors to the
expected shapes, which is what makes it usable as a bench.

`exoplasim/scripts/bench_pyburn_read.py` is that bench: `--make` writes the file,
`--bench` times the reader on it, and `--verify --against <old pyburn.py>`
requires every variable and every `readfile()` key to match in value, shape and
dtype. Everything below was produced with it and is reproducible from it.

## Result

| | wall |
| --- | --- |
| `readallvariables`, quadratic accumulation fixed | 30.53 s |
| `readallvariables`, payload read with `np.frombuffer` | 0.96 s |

32x. Each timing is its own process, so neither carries the other's resident
memory; run head to head in one process both inflate, to 36.6 s against 0.64 s,
and that is the pair the `--verify` output shows.

**All 116 variables and all 118 refactored `readfile()` keys are identical in
VALUE, SHAPE AND DTYPE**, at full scale. Dtype is part of the check on purpose:
the shipped reader produced float64 by handing `np.asarray` a tuple of Python
floats, and a `np.frombuffer` reader that returned the file's native float32
instead would be silently cheaper and silently different downstream, where
`wap`, the streamfunction and the pressure gradients are derived in whatever
precision they are handed.

**The read was essentially the whole of what remained.** The section above
measured 29.4 s for the entire postprocess of a 2.4 GB orbit; this measures
30.5 s for the read alone of a 2.65 GB one, which is the same number scaled.
Postprocessing an orbit is now a second or two against 77 to 93 s of model,
rather than the third of it that the fix above left.

## The defect

With the quadratic accumulation gone, what was left in `readrecord` was

    data = struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl])

which builds a tuple of Python float OBJECTS -- 8192 of them for one T42 grid
record, of order 6.6e8 per orbit -- and hands it straight to `np.asarray`, which
immediately unboxes them again. It also slices a fresh 32 KB bytes object per
record and builds an 8193-character format string to do it.

The payload is now decoded with `np.frombuffer` in a new `_decoderecord`, which
creates no Python objects at all. The views it returns are promoted to float64
once per code at the concatenate, and float64 is exactly what `np.asarray` gave
the old tuple, so both values and dtype are unchanged. Word length is still
derived from the ratio of the record's length in bytes to its length in words,
as `_getknownwordlength` derived it, but without re-reading the header to do so
-- the shipped `readrecord` unpacked every header twice.

## Numba was the wrong tool, and this is why

The question that prompted this was whether the remaining cost could be taken
off with numba decorators. It could not, for two reasons that are worth keeping
because they generalise.

**The hot call is not compilable.** `struct.unpack` over a `bytes` object is not
supported in numba's nopython mode. Reaching it would mean first rewriting the
loop to work on a numpy view of the buffer -- and that rewrite IS the fix, so
the JIT would then have nothing left to compile.

**Nothing downstream is a numba shape either.** In `mode='grid'`, which is what
the pipeline uses, a grid variable passes through `_transformvar` untouched
(`gridvar = variable`, no arithmetic), and a spectral one goes to
`pyfft.sp2gp`, which is already f2py'd Fortran. The per-cell triple loops that
would have been a genuine numba candidate were vectorised in the fix above. The
writer is I/O and compression bound.

The general form: **before reaching for a JIT, check whether the interpreter is
doing arithmetic or just manufacturing objects.** Here it was manufacturing
objects, and the cure for that is to stop, not to compile it faster. No
dependency was added.

## Two functions that could never have run

`readvariablecode` and `_gettimevar` both carried the same reader pattern and
both call `readrecord(fbuffer,n,en,ml)` without its `mf` argument, so both raise
TypeError on their first line of work; `readvariablecode` additionally tested
`if not variable:` on an ndarray. Nothing in this project calls either, and the
section above noted the second defect without noticing that the first made it
moot. Both are fixed and now agree with `readallvariables` record for record.

## What this changes elsewhere

The conclusion above to **run at 16 ranks** rather than overlapping two 8-rank
runs is strengthened, not moved: the postprocessing that concurrency trick was
hiding is now a rounding error against the model.

The per-orbit wall clock quoted in `docs/src/pipeline/costs.md` should fall by
roughly the postprocessing share. That has NOT been measured end to end, because
it needs a run made after this change and none has been; the figure there is
marked accordingly.
