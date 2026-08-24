# Every 16-rank binary was built single precision

Found 2026-08-18, while migrating ExoPlaSim to a subtree.

## What was true

`rebuild_binaries.py` invoked the model build as

    ./compile.sh -n <ranks> -p <ranks> -r <res> -v <layers>

`-n` is ncpus and is right. **`-p` is PRECISION IN BYTES**, and it was being
given the rank count.

`compile.sh` matches `-p` against `4`, `8`, `single` and `double`, and its
default is `prec=4`. So:

| build | flag | precision |
| --- | --- | --- |
| 8 ranks | `-p 8` | 8, correct by accident |
| **16 ranks** | `-p 16` | **4 -- matches no case, falls back to the default** |

`config/planet.yaml` declares `precision_bytes: 8`. Every 16-rank binary the
script produced was therefore single precision, including the T42 p16 that every
production run uses.

## Why nothing caught it

ExoPlaSim's `Model(precision=8)` rebuilt the executable at `-r8` on first use.
So the runs really were double precision -- the physics is not in question --
but the binary the manifest recorded had been replaced before a single orbit was
integrated. `binary_manifest.json` described an executable that never ran, and
`--verify` compared runs against a binary that no longer existed.

That is why the migration surfaced it. Vendoring made the build directory
stable and watched, so the rebuild-on-first-run showed up as a sha changing
under a clean tree instead of being invisible inside `.venv`.

## The evidence

Rebuilding with the precision read from the config changed **exactly the three
16-rank binaries and neither 8-rank one**, which is what the bug predicts:

| binary | before | after |
| --- | --- | --- |
| t21_l10_p8 | 72506ec46a872860 | 72506ec46a872860 |
| t42_l10_p8 | 7fe367ade133c662 | 7fe367ade133c662 |
| t21_l10_p16 | 5e14e23449538e06 | 972b979c8a47ccf3 |
| t42_l10_p16 | ae671f6b8eb0c1fa | **d1ef8bb078f241e3** |
| t85_l10_p16 | 9c8c75895b5a06cc | fedd8430dff34ebc |

`d1ef8bb078f241e3` is the binary ExoPlaSim had itself compiled at `-r8` during a
run before the fix, which is the confirmation that the corrected build produces
what the model wanted all along. After the fix a run recompiles nothing and
`--verify` stays clean across it.

## What it does not change

No result moves. Every orbit on disk was integrated by a double-precision
binary, because ExoPlaSim rebuilt before running. What was wrong was the
PROVENANCE: the manifest, and so `check_consistency.py`, were describing an
artifact that was replaced before use.

Precision is read from `config/planet.yaml` rather than passed positionally, so
the declared value and the compiled one cannot drift again.

Two of this record's mechanisms are gone rather than fixed, and both remove the
class rather than the instance. `compile.sh`'s positional parse was replaced by
CMake and `build_model.py` on 2026-08-22, which check every argument against a
list and exit non-zero on an unrecognised value; precision arrives as
`PLASIM_PRECISION`. And the rebuild-on-first-use path that made the runs
correct while the manifest was wrong is removed: a `Model` that finds no
executable raises and names `build_model.py`, because building one silently is
how the registry stopped describing the binaries that existed.
`notes/audits/model-build-driver.md` has both.
