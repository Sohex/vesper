# There are no stored Legendre tables to choose, and the tuner cannot repeat itself

*Worldbuilding frame: a COMPUTE reading of the transform library the Vesper
project's climate model calls. Nothing here is about the simulated planet. Read
from source and probed at SHTns 3.7.5, revision 4e69ceb, the revision
`build_shtns.sh` pins, on 2026-08-22.*

CLIM-74 was opened to measure what `SHT_QUICK_INIT` costs against STORED
LEGENDRE TABLES. `nm` on the model binary shows only `_fly` kernels, which was
read as the library recomputing per transform what it could have streamed, and
`shtns-viability.md` carried a table of what those tables would occupy against
the 32 MB target -- 0.23 MB at T42 rising to 14.36 MB at T170 -- with the
prediction that stored tables win at T42, wash at T85 and lose at T170.

**There is nothing on the other side of that comparison.** In this revision the
matrix-based algorithms are absent, not unselected.

## What the source says

`init_sht_array_func` fills `sht_func` with `SHT_ODD`, `SHT_SV` and the
`SHT_FLY*` and `SHT_OMP*` slots. It never writes an `SHT_MEM` entry. `SHT_MEM`
survives in the `sht_algos` enum and in exactly one other place:

    static void choose_best_sht(...)
        const int on_the_fly_only = 1;      // only on-the-fly.
        ...
        i0 = SHT_MEM;
        if (on_the_fly_only) i0 = SHT_SV;   // only on-the-fly

and `shtns_set_grid_auto` opens with `const int on_the_fly = 1;  // only
on-the-fly algos are available`. The memory arithmetic in `shtns-viability.md`
is therefore arithmetic about a configuration the library cannot produce, and
the prediction it carried cannot be tested.

## So what the two modes actually differ in

Both are Gauss grids; `sht_quick_init` is `sht_gauss` with the tuning pass
skipped. The choice is between ON-THE-FLY VARIANTS, which differ in unroll depth
(NWAY):

- `sht_quick_init` pins `SHT_FLY2` through `set_sht_fly`, with no timing.
- `sht_gauss` runs `choose_best_sht`, which times the available variants and
  keeps the winner.

The choice is made SEVEN TIMES, once per transform type -- `syn`, `ana`, `vsy`,
`van`, `gsp`, `v3s`, `v3a` -- so a configuration is a vector of seven picks, not
one setting. How many variants are available is per-resolution: `alg_lim` is
capped by `nlat_2 / VSIZE2`, and the model calls `shtns_use_threads(1)`, so
`alg_end = SHT_OMP1` and the OpenMP variants are never candidates. Probed with
`probe_shtns_algo.c`, the candidate set is `s+v` and `fly1` through `fly4`.

## The tuner cannot repeat its own answer

Five timed initialisations at each of two rungs, pick vector printed each time,
`SHTNS_VERBOSE=2`:

| run | T42 (syn ana vsy van gsp v3s v3a) | T170 |
| ---: | --- | --- |
| 1 | fly4 fly2 fly1 fly2 fly2 fly1 fly2 | fly2 fly2 fly2 fly2 fly1 fly3 fly1 |
| 2 | fly4 fly4 fly1 fly2 fly2 fly1 fly4 | fly2 fly2 fly2 fly4 fly2 fly3 fly1 |
| 3 | fly4 fly2 fly2 fly2 fly2 fly1 fly1 | fly4 fly4 fly2 fly4 fly2 fly2 fly1 |
| 4 | fly4 fly2 fly4 fly3 fly1 fly1 fly4 | fly2 fly3 fly2 fly3 fly2 fly3 fly1 |
| 5 | fly2 fly3 fly2 fly2 fly2 fly2 fly2 | fly4 fly3 fly2 fly2 fly2 fly3 fly2 |

**Ten runs, ten distinct vectors.** This is the mechanism behind archive
CLIM-44, where eight runs of one binary gave four distinct restart hashes, and
it is now located rather than inferred: a different transform function per run
is a different rounding order, and the model's own restart hash follows it.

The reason it cannot repeat is visible in the timings it decides on. At T21 the
`syn` race is `t(fly1) = 7.11e+03` against `t(fly2) = 7.12e+03`, and the winner
is declared on a tenth of a percent from a single unrepeated sample. Most of the
seven races at every rung are decided inside that margin. `SHT_GAUSS` is not a
setting that is faster or slower than `SHT_QUICK_INIT`; it is a lottery over
variants that are mostly indistinguishable, and it is disqualified for
production on determinism alone whatever a bench says.

One race is not in the margin. At T170 `syn` the tuner reports `fly4` at
1.60e+06 against `fly2` at 1.85e+06, a 13% gap on the transform type the model
performs most, repeated across runs. That is the one place a pinned alternative
to `fly2` might be worth having, and it is a claim to be settled by measuring
the MODEL rather than by trusting the library's single-shot tuner --
`docs/src/practice/failure-modes.md` class 34.

## The deterministic route exists, and it is a file

`SHT_LOAD_SAVE_CFG` makes `config_load` read `shtns_cfg` from the working
directory before any timing. Line 1529 of `sht_init.c` runs `choose_best_sht`
only `if ((quick_init == 0) && (!cfg_loaded))`, so a loaded config suppresses the
tuning pass in either mode.

The file is plain text and AUTHORABLE. One line per configuration: version,
SIMD id, `lmax`, `mmax`, `mres`, `nphi`, `nlat`, `grid`, `nthreads`,
`req_flags`, `nlorder`, `-1`, then the algorithm NAMES per variant and type as
`fprint_ftable` writes them. `config_load` matches on the grid parameters, the
request flags and `strcmp(simd, _SHTNS_ID_)`, and resolves each name through
`sht_func`, accepting only non-null pointers.

So a per-rung pick can be DECLARED rather than raced for: a file, per
resolution, whose contents are a measured choice and whose sha is its
provenance. `SHT_QUICK_INIT + SHT_LOAD_SAVE_CFG` loads that file when it is
there and falls back to the deterministic `fly2` when it is not -- which is a
fallback that must be made to announce itself rather than silently change the
transform.

`config_save` appends rather than truncates, and only when the tuning pass ran.
A run that is allowed to save is a run that raced.

## Pilot: the authored config works, and the tuner's one big claim was noise

`probe_shtns_variant_cost.c`, T170, single-threaded, 2001 repetitions per type,
each variant forced through an authored `shtns_cfg` and CONFIRMED by reading
back what the library says it loaded rather than assuming the request took.
Median microseconds per call, with the interquartile range beside it:

| type | fly1 | fly2 | fly3 | fly4 | IQR band |
| --- | ---: | ---: | ---: | ---: | ---: |
| syn `SH_to_spat` | 282.5 | 286.5 | 298.3 | **277.2** | 6 to 28 |
| ana `spat_to_SH` | 328.7 | 324.0 | 328.3 | **322.4** | 7 to 15 |
| vsy `SHsphtor_to_spat` | **512.9** | 519.5 | 547.7 | 581.8 | 11 to 18 |
| van `spat_to_SHsphtor` | 606.0 | **601.3** | 607.9 | 604.3 | 13 to 22 |

**The 13% the tuner reported for T170 `syn` is not there.** Against fly2, fly4
is 3.2% faster on `syn`, not 13%, and the difference the tuner acted on came
from single unrepeated samples. The direction of the effect is real; its size
was noise-inflated by an order of magnitude, which is what an unrepeated
instrument does.

**No uniform variant wins.** fly4 takes `syn` by 3.2% and loses `vsy` by 12.0%,
and `vsy` is the more expensive call. A uniform-fly4 configuration is a net loss
at T170 even though it holds the single best number in the table.

So the per-rung pick is a per-TYPE vector at each rung, which is what the
library's own structure already said, and the pilot's purpose was to establish
that a vector can be pinned at all and to measure the scatter that any adoption
threshold has to clear. Both are now known: the route works, and the scatter is
2 to 4% of a median.

WHAT THIS IS WORTH AT MODEL LEVEL, stated before the full sweep rather than
after it. `shtns-viability.md` measures SHTns's own kernels at 0.00% of T170
model runtime with the wrappers at 0.74% together. A few percent of a few
percent is a model-level effect of order 0.1%, against a bench self-scatter
floor of 5%. **No model-level A/B can resolve this choice**, and one run at a
declared tolerance can only confirm that the restart still agrees. The
transform-level measurement above is the only instrument that reaches the
effect, which is `docs/src/practice/failure-modes.md` class 34 applied before
the measurement instead of after.
