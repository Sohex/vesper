# Tasks

Atomic, trackable work items. **Findings documents do not carry their own to-do
lists** -- an audit says what is true, this file says what to do about it, and
the two are kept apart so a finding can be read without being re-litigated and a
task can be closed without editing the argument that produced it.

Every task cites the document that justifies it. If a task has no source, it is
not yet a finding and probably needs one.

Status: `open` | `doing` | `blocked` | `done` | `wontfix` (with a reason).

## Conventions

- IDs are stable and never reused. Prefix by area: `LITH` lithology, `GRAV`
  gravity, `VOLC` volcanism and weathering fluxes, `SURF` derived surface
  classes, `CLIM` climate, `HYD` hydrography, `BIO` biosphere, `REF` references
  and provenance.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point.
- A closed task moves to `archive/tasks.md` rather than staying here, so this
  file shows only what is left. Ids are never reused, so a cited id that is not
  here is there.

---

## Open

| id | task | source | status |
| --- | --- | --- | --- |
| VOLC-2 | Give Orogen an eruption-age or resurfacing-rate field, so `flood_basalt` and `oib` can be resolved as andisol-bearing or not. They are reported as UNDETERMINED rather than counted | `pedology/config/pedogenesis.yaml` | open |
| VOLC-3 | Place silcrete and diatomite against per-basin silica supply. A large share of silica delivery goes to closed basins, which is the setting for both, but per-basin placement needs hydrography run on the arc-bearing build | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline |
| VOLC-5 | Apply a glass-dissolution enhancement to vitric terrain. Deferred deliberately: field rates on andesitic ash run ~20x laboratory, but only 0.45% of land is andic and 0.60% vitric, so the correction is small and needs VOLC-2 to place | `references/INDEX.md` | open |
| BIO-1 | Regenerate the LPJ-GUESS driver: `vesper_driver_provenance.json` is pinned to a superseded config and `check_consistency.py` has been failing on it. Its provenance also records no generator, so the check cannot name the remedy | `scripts/check_consistency.py` | blocked on the re-baseline -- the driver is built from a climatology, so rebuilding it now would only have to be redone |
| LITH-5 | Explain or correct the composition divergence from GLiM | `notes/audits/orogen-lithology.md` | open -- RE-MEASURED on the arc build and mostly resolved. The volcanic shortfall was the unreachable arc rules (0.49x -> 1.28x) and the metamorphic excess halved (2.06x -> 1.42x). Two new divergences appeared and are one cause: plutonic doubled because the unroofed arc root fires, and carbonate fell because the coastal forearc overwrites shelf. What remains genuinely unexplained is evaporite at 9.5x |
| LITH-6 | Decide whether `surface_rock` should be renamed; it carries consolidated lithology, and GLiM's largest land class (unconsolidated, 24.6%) has no counterpart | `notes/audits/orogen-lithology.md` | open |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here). Needs the glacial model rebuilt, not a factor bolted on; justify on its own merits | `notes/audits/orogen-gravity.md` | open |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline -- there is no climatology at all now, and `baseline_climatology` is null until one exists |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | blocked on the re-baseline, same reason as SURF-1 |
| SURF-3 | Loess: needs a dust emission scheme and a transport path, neither of which exists | `pedology/notes/derived-surface-classes.md` | blocked |


## Closed

Moved to `archive/tasks.md`, with the reasoning that closed each one. Ids are
never reused, so look there when a commit cites an id this file does not have.
