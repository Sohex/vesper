# Declared checks that build their comparison arm from a shorter argument list

This is worldbuilding. Vesper is a fictional super-Earth, and the checks below
are gates in the code that generates it: they are guards on a simulation, not on
anything real.

Swept 2026-08-24 over `scripts/`, `lib/`, `hydrography/scripts/`,
`exoplasim/scripts/`, `biosphere/`, `pedology/`, `minerals/`, `aeolian/`,
`analysis/` and `maps/`. Every `--check`/`--verify`/`--test`/`--self-test` flag,
every round trip, and every function that builds a comparison arm.

## The shape

A check constructs a second run -- a re-solve, a reference arm, a round trip, a
control, an independent recomputation -- and builds it from a SHORTER argument
list than the thing it is checking. It passes for as long as the omitted terms
are off or absent, and it silently stops testing what its name claims the moment
one of them lands. It cannot fail in the direction that matters, so by
`docs/src/practice/conventions.md` it is not a test.

Three variants turned up often enough to be worth naming separately.

- **The short arm.** The instance that prompted the sweep: a re-solve missing
  two terms the primary solve carries.
- **The fixture that agrees with itself.** The comparison is built by the same
  code it checks, so it proves only that that code is self-consistent. Named as
  CONS-9 in `run_exoplasim.py:expected_namelist_keys`. It is not automatically a
  defect: `scripts/check_consistency.py:self_test` builds its run fixture from
  `expected_namelist_keys` deliberately, asserts the fixture GREEN first, and
  then breaks exactly one thing per case with the verdict named in advance.
  That is the correct mitigation and it is commented as such.
- **The verdict nobody acts on.** A bar is declared, evaluated, printed and
  written to a report, and then dropped: the process exits 0. This turned out to
  be the most common form by a wide margin.

## What was found, by component

Verdicts: DEFECT was fixed by forwarding the argument unless a row says
otherwise; JUSTIFIED means the difference is deliberate and the justification
now sits beside the code.

### `scripts/`

| check | checked call | checking call | difference | verdict |
| --- | --- | --- | --- | --- |
| `smoke_test.check_build_scoped_defaults` | a script's per-build data directory | `component_data(c, cfg).name` vs `cfg["source_build"]` | the second is the string the first builds the path out of; the component list was the literal `("hydrography", "pedology")` | DEFECT. Rewritten to read the component list and the expected directory from `config/pipeline.yaml`, so the graph and the resolver are two declarations that can disagree |
| `smoke_test.check_purge_never_reaches_the_terrain` | `pipeline.cmd_purge` -> `downstream(target, by_id)` | `downstream(seed, by_id)` over every seed | none; `map_affecting` belongs to `cmd_plan`, not to purge | JUSTIFIED |
| `smoke_test.check_slope_fit` | `Export.local_slope_deg` | the same function on a synthetic icosphere | no ET-shaped omission; the reference is analytic and the tolerance sits between two real implementations | JUSTIFIED |
| `smoke_test.check_gate_filter_matches_config` | the drivers' `filterkappa`/`nfilterexp` | `config/planet.yaml` | only those two keys; another hardcoded model parameter would drift the same way unseen | JUSTIFIED as far as it goes, and narrower than its name |
| `smoke_test` file list | the tree | `d.glob("*.py")` per `SCRIPT_DIRS`, plus `*.sh` for the two rung lints only | subdirectories are invisible to every lint | JUSTIFIED, argued at the two rung lints |
| `check_consistency.self_test` | `audit_runs(rep, runs, runner)` | `audit_runs(r, root/"runs", runner)` | none but `runner.INPUTS`, redirected for the fixture and commented | JUSTIFIED |
| `check_consistency.self_test` fixture | the staging lists | built from `runner.expected_namelist_keys` | the fixture-self-agreement shape, mitigated by asserting green first and then breaking one thing per case | JUSTIFIED, commented |
| `link_worktree.py --check` | `ignored_entries(main, model_binaries)` | the same call | none for the entry set, but the exit ignored `misdirected` | DEFECT. A link at `X` pointing at `Y` is the deep-link bug the module docstring names as its motivating incident, and it is not "missing", so `--check` printed the warning and exited 0 beside it |

### `hydrography/scripts/`

| check | checked call | checking call | difference | verdict |
| --- | --- | --- | --- | --- |
| `build_groundwater --uniqueness-check` | `gw.solve(..., et_max_m_s, et_lambda_m, fixed_head_m, aquifer_base_m, min_saturated_m)` | the same, plus `start_all_free` | none as it now stands; the two lists were written out separately, which is how GW-15's sink and GW-17's baselevels came to be in one and not the other | DEFECT of construction. One `solve_kwargs` is built and both calls take it whole |
| `groundwater.confined_test` | `solve` | `solve` without `et_max_m_s`, `aquifer_base_m`, `min_saturated_m` | the analytic head is the solution with no sink; the aquifer base selects the other form | JUSTIFIED, now said beside the call |
| `groundwater.dupuit_test` | `solve` | `solve` without `et_max_m_s` | the Dupuit parabola is the solution with no sink | JUSTIFIED, now said beside the call |
| `carve_verdict` humidity bracket | `penman_open_water(...)` | the same plus `column_relative_humidity` | a strict superset | NOT-THIS-CLASS |
| `earth_calibration.stage_diagnostics` | `stage_solve(..., river_km2)` | a literal `1e4` | `--river-km2` ran and changed nothing for this stage | DEFECT. Forwarded; the default is unchanged, so every table the notes cite reads as it did |
| `earth_calibration` result key | `--et-lambda` selects which solution is scored | the key had `region/confinement/drain/surface/edge` | two lambdas wrote over each other under one name | DEFECT. Appended when set, so existing keys keep their names |
| `earth_calibration.stage_benchmark` | `solution_name(drain, et_lambda)` | `earth_solution.npz` by name | nothing here reads the head or depth; the file is opened for masks no variant moves, and a variant run fails loudly on a missing file | JUSTIFIED, now said beside the load |
| `earth_calibration --surface` | -- | -- | already fixed and commented; the per-cell aggregation read `sol["depth"]` and the option ran, reported and changed nothing | CLOSED, recorded in `earth-calibration-criterion.md` |
| `validate_lake_solver` | `lake_balance.solve` | the analytic relation, not the solver | the hypsometric inversion and the overflow cascade are out | JUSTIFIED, stated in the module docstring |

### `exoplasim/scripts/`

| check | checked call | checking call | difference | verdict |
| --- | --- | --- | --- | --- |
| `corrk_cross_check.run_checks` | `flux_fractions(t376, sun)` and `flux_fractions(t376, star)` | `run_checks(t376, t1000, sun, ...)` | `star` was never a parameter, so the stellar partition was unchecked -- and `h2osww` and `co2sww` are both a stellar broadband over a solar one | DEFECT. `star` forwarded; the partition runs over both spectra and both pass |
| `run_exoplasim.verify_stellar_spectrum` | `stage_stellar_spectrum` writes `STARFILE` and `STARFILEHR` | `NSTARFILE != "1" or not STARFILEHR` | `name` was bound from the config and used only in error messages; `STARFILE` was never read | DEFECT. Both entries must name the configured spectrum and both files must exist |
| `rebuild_binaries --verify` | the rebuild clears `RUN` and `BIN`; `MATRIX` declares what must exist | globbed `RUN` only, and dropped `missing` from the exit | a stale executable in `BIN` was invisible; "IN THE MATRIX AND NOT BUILT" exited 0 | DEFECT, both halves |
| `restart_surface --self-test` | `verify_restart_surface_fields(run, restart, codes, manifest=, allow_superseded=)` | the same without either keyword | the donor-surface branch was never executed by the self-test; and `cells = 64 * 128` crashed it on any rung but T21 | DEFECT. Cases D and E added; the grid comes from the fixture |
| `convert_restart --check-template` | "run every refusal" | returned after `check_compatible` | the seed-size, whole-level, non-finite and partner-ordering refusals were all past the return | DEFECT. `seed_override` is computed first and the conversion runs with the records discarded |
| `convert_restart` explicit seed | the donor path refuses a size mismatch | `--seed` had no length bar at all | a wrong-width seed wrote a malformed restart in silence | DEFECT |
| `restart_convert_selftest` | `cv.convert(..., seed_override=)` | never called with it | `--seed` and `--keep-template-seed` had no case | DEFECT. Both directions added |
| `restart_convert_selftest` length control | the length heuristic the schema replaces | `_require(len(ints) == len(reals))` on two payloads it built | an arithmetic identity, with no heuristic run on either side | DEFECT. It now finds the byte lengths this donor shares between different `put_restart_*` writers, and fails if there are none |
| `transform_exactness` round trip | `STORAGE`, declared, "must reach" | printed a NOTE and returned 0 | the one row with a bar could miss it invisibly | DEFECT |
| `verify_legendre_parity` | the parity premise | printed FAIL, exited 0 | -- | DEFECT |
| `shortwave_band_weights` Earth CO2 | `EARTH_CO2_SHORTWAVE_W_M2`, declared | printed `inside`/`OUTSIDE`, raised nothing | every weight the file quotes was produced under a bar that could not stop it | DEFECT. Raises after the report is written |
| `shortwave_band_weights` LH74 Eq. 21 | -- | ratio min/max/median | no declared bound at all | OPEN, `world-zvk2`. A bound chosen now would be chosen after the run it judges |
| `stack_floor --validate` | a binary at one (rung, levels, threads) | a parse at `--rung`/`--levels`/`--threads` | nothing reconciled them; with no `--rung` it printed one percentage per rung | DEFECT for the configuration, fixed from the binary's registry name. The missing bar is `world-f1u3` |
| `verify_weight_factorisation.compare` | `legmod.f90:legini` | a Python transcription of it, against a second restatement of the same nine quantities | every one of the eight is the identical product reassociated; the negative control cancels a transcription error in either arm | DEFECT, OPEN as `world-g5va`. The fix needs `verify_inverse_transform`'s `extract`, which compiles |
| `predict_ocean_terms.validate_operator` | the operator | `allsea`, `k=1.0` | the control the identity requires | JUSTIFIED |
| `verify_*.sh`, `close_*_energy.py`, `restart_schema`, `dust_forcing` | -- | -- | each compiles the model's own modules or carries an independent right answer, several citing class 17 | NOT-THIS-CLASS |
| `bench_*`, `sweep_*`, `score_*`, `reproducibility_matrix`, `compare_*` | -- | -- | timing and factorial arms; no reference recomputation | NOT-THIS-CLASS |

### `lib/`, `pedology/`, `aeolian/`, `maps/`, `biosphere/`, `minerals/`, `analysis/`

| check | checked call | checking call | difference | verdict |
| --- | --- | --- | --- | --- |
| `lib/stellar.solar_partition_identity` | a real spectrum: `read_hires` + `_bands` + `_partition` | `_blackbody_grids` + `_planck` + `_partition` | the two shared `_partition` alone, so `assert_model_grid`, the row split and the `minwavel` deletion were outside "the one statement with a right answer" | DEFECT. The same Planck curve now also goes through `_bands`, at the tolerance already declared; the two arms differ by 4.8e-05 because the grid's first row sits a float below the cut and is deleted |
| `pedology/validate_against_earth` | `build_soil.weather_texture` | an inlined transcription | no `sand_to_silt_loss_ratio`, neither cap on the donor pools, and `clay_yield` read through a `.get(..., 1.0)` default the model does not have | DEFECT. The caps bound the model BELOW the arm, so the arm could only overstate |
| `pedology/brine_paths.single_lithology_check` | both routings are published | `args.weighting`'s alone, and `passes` never gated | half the numbers were never held against their own right answer, and the report went out with a failed definitional identity in it | DEFECT, both halves |
| `pedology/build_surface_classes.audit_rules` | -- | reports `never_fires`/`always_fires` without gating | argued for in `surface_classes.yaml` | JUSTIFIED |
| `pedology/thermostat_efficiency --compare-builds` | -- | fixed climate across builds | the fixed climate is the stated point | JUSTIFIED |
| `aeolian/build_dust_source_fields --self-test` | `build_dust.emission_over_weibull`, survival-space quadrature | the same rule, in `aerocore.f90` and in the mirror of it | none. The model's quadrature moved to survival space and the identity holds at 2.3e-16 against a declared 1e-10 | SETTLED, `world-2lsg`. The q-space rule's node ceiling did not move with the threshold, so it returned exactly zero wherever `u*t` cleared it; the derivation named the in-model half as the wrong one and that is the half that changed |
| `aeolian/build_dust_source_fields` drag | per-cell `z0_aeolian` from the erodible mosaic | the bracket-end plane of `source_fractions`' roughness stack, and a log-uniform `z0` fixture | none on this axis. `ln(zref/z0)` still takes a scalar in-model and cannot be folded into a static field, which is `world-h24h` | SETTLED, `world-2lsg`. A constant fixture cannot satisfy the identity any more |
| `aeolian/sea_salt_source._earth_check` | "the source functions and the mass integration" | the agreement of two ratios with each other | a common-mode error moves both together; and `_integrate_mass` re-implements the trapezoid rather than calling `mass_flux` | OPEN, `world-fgak`. Both halves need a threshold fixed in advance |
| `aeolian/build_sea_salt --bracket` | `solve(cfg)` | `solve(deepcopy(cfg))` | full re-run | NOT-THIS-CLASS |
| `aeolian/dust_runoff_sensitivity.perturbed` | -- | scales `pr` and `evap` only | the docstring covers it; the payload note overstates what the penman arm sees | JUSTIFIED-but-narrow |
| `maps/snapshot.basemap_is_current` | `build_basemap` records its resolution outside the fingerprint | the fingerprint alone | `--basemap-width` was forwarded only when a rebuild was already triggered, so a different width against an unmoved fingerprint kept the old raster | DEFECT |
| `biosphere/abiotic_nutrient_ledger.run_fixtures` | -- | declared expectations per fixture | the fixture-self-agreement shape done right | NOT-THIS-CLASS |
| `biosphere/build_vesper_pfts.self_check`, `bvoc_gate.check_run` | -- | -- | one tautological line carried by three real checks; `bvoc_gate` exits 1 | NOT-THIS-CLASS |
| `lib/sensitivity.verify`, `lib/rungs`, `lib/gridding`, `lib/paths` | -- | -- | genuine identities or external files; `sensitivity` declares itself arithmetic | NOT-THIS-CLASS |
| `minerals` ceilings, `analysis/trace_gas_forcing.check_fits`, `analysis/orogen_resolution_controls` | -- | -- | one-sided by design with `normalise` raising; external right answers; the control export recipe is `source/README.md`'s own `$COMMON` | NOT-THIS-CLASS |

## What the recorded results now mean

- **`hydrography/notes/water-table-convergence.md`, the 0.000e+00 uniqueness
  identity.** Measured 2026-08-20, before GW-15's sink and GW-17's baselevels
  existed, so both trajectories carried neither. It is an honest identity for
  the equation of that date and is not a result for the equation the component
  now solves. The note says so and names the command that would produce one.
- **`pedology/README.md`'s clay bias.** `+0.062` mean bias, correlation 0.840,
  and a mafic-felsic divergence of 0.231 against an observed 0.179 were measured
  against a scoring arm that lacked the model's own caps. Through
  `weather_texture` the same 16 localities give a mean bias of `+0.004`, a
  correlation of 0.845, and a divergence of 0.152. The verdict inverts: the
  model understates the spread rather than overstating it by 29%.
- **`corrk_cross_check`'s "all three hold".** Every run before this sweep
  checked the solar partition only. Re-run: both partitions close, 0.97204
  solar and 0.99217 stellar over the table span, with no gaps.
- **Anything quoted from `shortwave_band_weights`, `verify_legendre_parity` or
  `transform_exactness`.** Each carried a declared bar that could not change the
  exit code, so a past green result is evidence that the script ran, not that
  the bar was met. The reports on disk carry the numbers; read them rather than
  the exit.

## What is not settled

Three checks could not be repaired here and each has an issue: `world-g5va`
(the weight factorisation comparison, which needs a compiled extraction),
`world-2lsg` (the dust quadrature divergence, which needs every binary rebuilt),
and `world-fgak` plus `world-zvk2` and `world-f1u3` (three bars that must be
fixed before the run that judges them, not after).

Two checks could not be RUN in the worktree this was swept from, and their
fixes are static: `restart_surface --self-test` needs a run directory with a
staged code 229, and `brine_paths` and `shortwave_band_weights` both need a
`baseline_climatology`, which `config/planet.yaml` does not yet declare.
