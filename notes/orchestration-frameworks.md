# Adopting an external orchestration or provenance framework

*Written 2026-08-22. Claims about this tree are from reading it; claims about
each tool are from its own documentation and from its configuration files, cited
inline. Where a tool's behaviour was checked against an artifact here, the
artifact is named.*

Vesper is a fictional planet and this project simulates it. What follows compares
this project's own pipeline machinery against four external frameworks, so that
the comparison is not run a second time. It is not a decision to adopt anything;
the one decision it reaches is which framework is worth a spike, and what the
spike has to answer.

## The subject of the comparison

Three separable things get called "the pipeline" here, and only the first two are
what an external framework would replace.

**A declared graph.** `config/pipeline.yaml` states every step, what it writes,
what must precede it, its cost band, its loop membership and its gate. Steps are
declared whether or not any of them has ever run. `scripts/pipeline.py` queries
that declaration: `--status` for presence and config drift, `--plan` for the
ordered closure to a target, `--register` for the artifact table, and `--purge`
for reverse reachability with the traversal cut at `orogen`. It runs nothing.

**A record of what happened.** `exoplasim/runs/INDEX.json` is the only record of
what a UUID-named run was. Artifacts carry `vesper_source_build`,
`vesper_geography` and `vesper_config_sha256`, and `lib/provenance.py` verifies
them at the point of reading.

**Model configuration and execution.** `run_exoplasim.py` and
`continue_exoplasim.py` stage namelists, verify the stellar spectrum, match the
binary against the sha of the source it was built from, hash the config, assign
the run UUID and register it. Segments chain through restart files.

## ESM-Tools

**Verdict: no. The mismatch is the machine, not the models.**

`configs/components/` holds 28 entries, and not all of them are models: the
model configs `echam`, `fesom`, `icon`, `nemo`, `oifs`, `pism`, `jsbach` and
`lpj_guess` sit beside couplers and I/O libraries (`oasis3mct`, `xios`, `yac`,
`yaxt`). The machine library is named clusters: `albedo`, `aleph`, `blogin`,
`glogin`, `juwels`, `levante`, `mistral`, `nesh`, `olaf`, `ollie`, and a
`generic.yaml`.

`configs/machines/generic.yaml` is not a workstation mode. It sets
`batch_system: "slurm"` and carries an explicit guard: running an actual
simulation from it is an error, and it exists for local testing and
`esm_runscripts -c` dry checks only. ESM-Tools' unit of value is a scheduler
submission on a named cluster, and this project runs `mpiexec` on one host.

The design overlap is nonetheless real and worth recording, because it is the
shape of `continue_exoplasim.py`. The ESM-Tools workflow manager runs
`newrun -> prepcompute -> compute -> observe_compute -> tidy`, with a resubmit
that starts the next run, custom subjobs placed by `run_after`/`run_before`, a
`trigger_next_run` flag, and chunk controls (`run_only: first_run_in_chunk`,
`skip_chunk_number`). That is segment chaining with the resubmit loop, factored
better than ours.

The `lpj_guess` config is a git clone of the DKRZ ec-earth 4.1.1 fork plus a
cmake line, which is the pattern in miniature: what ESM-Tools supplies is a
curated build recipe per model, version and machine, and all three models here
are forks it does not carry.

## ESMValTool

**Verdict: no, on the same argument that already refuses CDO.**

ESMValTool is evaluation rather than orchestration, so it competes with
`analyze_climatology.py`, `error_budget.py`, `close_*_energy.py` and
`assess_convergence.py`, not with `pipeline.py`.

Checked against `exoplasim/runs/run_4182235e9781/MOST.00003.nc` and
`exoplasim/analysis/climatology/bootstrap_regular_climatology.nc`:

- `time` carries `units = "timesteps"`, `standard_name = "timestep_of_year"`,
  values 255, 735, 1215, 1695. It has no epoch and no calendar attribute, so it
  is not a CF time coordinate and iris does not load it as one.
- `lat` and `lon` carry `units = "deg"` rather than `degrees_north` and
  `degrees_east`, and no coordinate system is declared.

The geometry half is fixable. iris `area_weights` takes the radius from the
cube's coordinate system when it is present and spherical, and only otherwise
falls back to its Earth default, so a declared spheroid would be honoured. That
fix is worth making independently of this comparison and is CONV-5.

The time half is not fixable, and it is where the preprocessor's value is. No CF
calendar expresses this orbital period; the seasons here are orbital longitude
quadrants rather than months, so `seasonal_statistics` has nothing to bind to;
and the twelve bins are uneven by construction, which is why `lib/climatology.py`
weights by record count, where `climate_statistics` would weight by its own month
lengths.

The diagnostic library is the other half of the offering and it is comparison
against Earth observations and CMIP ensembles, which has no counterpart here.
What this project checks instead are conservation identities.

`docs/src/reference/environment.md` already refuses CDO on this argument:
calendar climatologies are unusable on a `units = timesteps` axis and Earth
remapping is a trap. ESMValTool is that argument applied to a larger tool.

## Cylc

**Verdict: not assessed.** Examination was declined at the point the question was
asked, so nothing here supports a verdict either way. The one observation
carried forward is that Cylc cycles over cycle points, which are dates or
integers, while every loop in `config/pipeline.yaml` exits on a predicate over
results.

## AiiDA

**Verdict: the candidate. Not adopted, and what stands between is a spike rather
than an argument.**

AiiDA covers more of this project's machinery than its materials-science plugin
ecosystem suggests, and four properties are worth stating because each one
removes a reason to dismiss it.

**It does not require a database server or a broker.** `verdi presto` creates a
profile on SQLite storage in place of PostgreSQL, and falls back to a built-in
ZeroMQ broker when RabbitMQ is absent. A broker-less profile is polling-based,
cannot run the daemon or `submit()` asynchronously, and cannot pause, play or
kill a process from another terminal, because those go over RabbitMQ's RPC.
Processes launched with `run()` execute in the calling interpreter and block it.
For runs in the `hours` cost band the daemon is wanted, so adoption implies
installing `rabbitmq-server`.

**It does not require a plugin per step.** `aiida-shell`'s `launch_shell_job`
runs any command-line program, writing the CalcJob and parser on the fly. Every
step in `config/pipeline.yaml` is a script with a command-line interface.

**It has an invalidation traversal, and a better-founded one than `--purge`.**
`verdi node delete` and `delete_nodes()` walk forward over the provenance graph
under configurable rules -- `create_forward` on by default, plus
`call_calc_forward` and `call_work_forward` -- and are a dry run until told
otherwise. The graph is an instance graph rather than a graph of step types, so
loop A unrolls into successive nodes and no equivalent of the `orogen` cut is
needed. That cut exists in `scripts/pipeline.py` because a step-typed graph
cannot distinguish one iteration from the next.

**It declares work before running it.** A `WorkChain` states its sequence in
`spec.outline()`. `aiida-workgraph` builds the graph as an object with `If`,
`While` and `Map` zones and renders it with `wg.to_html()` before execution, and
its documented convergence pattern is a `While` zone whose predicate compares the
two latest results against a threshold. A `CalcJob` run with
`metadata.dry_run=True` prepares the full submission folder without submitting
and reports the path in `dry_run_info`.

### Extras carry whatever this project needs a node to say

Node **attributes** are immutable once the node is stored, because they define
the node's content and changing them would invalidate the provenance of its
descendants. **Extras** are the mutable half: user-defined key-value pairs, set
after storing, queryable alongside attributes. Arbitrary project metadata on a
node -- a terrain hash, a cost band, the text of a gate -- is a supported use
rather than an extension.

Ancestry is queried directly. `QueryBuilder` joins with `with_ancestors` and
`with_descendants`, so "which build was this climatology computed on" is a link
query over recorded provenance and is correct by construction, where this
project answers it today by stamping `vesper_source_build` and checking it in
`lib/provenance.py:require_build`. For a value that enters from outside the
graph, the identity goes in an extra and is queried the same way.

### Gates are expressible, and better than the form they take here

A `gate:` in `config/pipeline.yaml` is a condition that must be answered by a
person before a step is meaningful. AiiDA supports that pattern directly: a
`WorkChain` runs a `while_` loop that pauses, poses its question by setting an
extra on its own process node, waits for a person to answer by writing the
reply into that extra, and resumes. `verdi process pause` and
`verdi process play` drive it by hand.

That is stronger than what this project does now, in one specific way: the
human's answer lands inside the provenance graph rather than in a commit
message or in nobody's record at all.

### What large artifacts cost

`use_symlinks` does not answer the question it appears to answer. It applies to
input staging: the contents of a `RemoteData` input are copied into the working
directory by default, and the option symlinks them instead.

Outputs work differently. Files named in the retrieve list -- which is what
`aiida-shell`'s `outputs` populates -- are copied into a `FolderData` node in the
repository. Files not named there stay in the working directory, and every
`CalcJob` emits a `remote_folder` output of type `RemoteData` that points at it,
which a later calculation reads through `remote_copy_list` without transferring
it.

So a multi-GB climatology can be kept out of the repository. What follows is an
arrangement to make rather than a guarantee lost: the content-addressed store
hashes and deduplicates what it holds, so an artifact left in the working
directory is covered by neither. Integrity on it is then something the wrapped
command supplies -- emit a checksum beside the artifact and retrieve the
checksum, which is small -- so the question for the spike is what that
arrangement costs in practice, not whether the property is available.

### The shape adoption would take

AiiDA would replace the record of what happened: `INDEX.json` as the run
registry, `downstream()` and `--purge`, and part of the stamping in
`lib/provenance.py`.

Whether `config/pipeline.yaml` survives is a judgement about legibility and not
a capability question. A WorkGraph carries the same graph, and extras carry the
cost bands and gates. What it does not carry as well is a reader: the YAML file
states 56 steps, their costs and their gates to someone who has never run the
pipeline and is not running it now, where the equivalent WorkGraph is a program
that has to be executed or rendered before it says anything. That is a real
consideration and it is a preference, so it should be decided as one.

Adoption belongs at a generation boundary. Rule 7 makes everything below a new
generation worthless anyway, so the existing tree is purged once by the existing
`--purge` and nothing is ingested. That removes the migration entirely, and it is
the reason the cost of trying is low.

### Before recording another limitation, check it

Four claims of the form "AiiDA does not do this" were made about this project's
requirements and all four were false: that it needs PostgreSQL and RabbitMQ,
that it needs a plugin per step, that it has no invalidation traversal, and that
it only records what has already run. Three further ones -- that it cannot ask
the semantic build question, that gates have to live outside it, and that an
un-retrieved artifact loses integrity checking irrecoverably -- were false in
the same way.

So the decision procedure for anything that looks like a gap: check it against
node extras, a custom parser, `aiida-shell` and `aiida-workgraph` before
recording it. The gaps that are real will be architectural consequences of
immutable stored nodes and of provenance meaning recorded links, not missing
features.

### What the spike has to answer

Every question below is about effort or fit. None is about whether AiiDA can do
the thing.

1. What does keeping a multi-GB climatology out of the store actually cost --
   does the `remote_folder` route work as documented, and does a retrieved
   checksum sidecar give back what the content-addressed store would have given?
2. Does `aiida-shell` wrap the step scripts as they stand, or do they have to be
   restructured to expose their inputs as nodes?
3. Can an Orogen generation sit inside the graph? Loop A exits by carving the
   intersection of the verdicts taken at the two bounding climates, which the
   `while_`-and-pause pattern expresses; the open part is the generation itself,
   which runs a vendored fork and writes the read-only `source/` tree.
4. What does `verdi node delete --dry-run` report against `pipeline.py --purge`
   on the same seed? This is a differential test with a right answer: the two
   should agree on everything within one pass, and the instance graph should be
   correct where the `orogen` cut is an approximation.
5. Is a WorkGraph program an acceptable replacement for reading
   `config/pipeline.yaml`? This one is a preference and is settled by looking at
   both, not by measuring.

CONV-6 carries the spike.
