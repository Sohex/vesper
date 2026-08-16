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
  gravity, `SURF` derived surface classes, `CLIM` climate, `HYD` hydrography,
  `BIO` biosphere, `REF` references and provenance.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point.

---

## Open

| id | task | source | status |
| --- | --- | --- | --- |
| LITH-1 | Assign evaporite mineralogy downstream from the chemical divide, not from Orogen's single geometric `evaporite` class | `notes/audits/orogen-lithology.md` | open |
| LITH-2 | Make surface albedo mineralogy-dependent for closed-basin fill; the single 0.50 halite value is worth ~1.7 W/m2 per 0.10 of error | `notes/audits/orogen-lithology.md` | open |
| LITH-3 | Ground the 60 uncited rock-class constants (erodibility, albedo, density) or mark each explicitly as a stylistic choice | `notes/audits/orogen-lithology.md` | open |
| LITH-4 | Get correct DOIs and fetch bare-rock / evaporite albedo literature; three title-only queries matched nothing | `notes/audits/orogen-lithology.md` | open |
| LITH-5 | Explain or correct the composition divergence from GLiM: metamorphic 2.1x Earth, evaporite ~10x, volcanic 0.5x | `notes/audits/orogen-lithology.md` | open |
| LITH-6 | Decide whether `surface_rock` should be renamed; it carries consolidated lithology, and GLiM's largest land class (unconsolidated, 24.6%) has no counterpart | `notes/audits/orogen-lithology.md` | open |
| GRAV-4 | State `n = 1` explicitly at the erosion call site; the whole gravity correctness argument rests on it and nothing says so | `notes/audits/orogen-gravity.md` | open |
| GRAV-5 | Bracket the slope-exponent uncertainty: generate one build at `g^(-1/2)` instead of `g^(-1)` and compare basin statistics. One pass, and it settles the only live question | `notes/audits/orogen-gravity.md` | open |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here). Needs the glacial model rebuilt, not a factor bolted on; justify on its own merits | `notes/audits/orogen-gravity.md` | open |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-3 | Loess: needs a dust emission scheme and a transport path, neither of which exists | `pedology/notes/derived-surface-classes.md` | blocked |
| REF-1 | Hardie & Eugster (1970), the origin of the chemical divide, could not be fetched and is not held; cited through Eugster & Jones 1979 and Deocampo & Jones 2014 | `references/INDEX.md` | open |

## Done

| id | task | source | closed |
| --- | --- | --- | --- |
| LITH-0 | Determine whether the missing gypsum class is an Orogen oversight | `notes/audits/orogen-lithology.md` | wontfix -- it is not. Orogen assigns evaporite geometrically and has no basis for mineralogy; the decision belongs downstream. Superseded by LITH-1 |
| GRAV-0 | Determine whether unscaled bathymetry is a bug | `notes/audits/orogen-gravity.md` | wontfix -- g cancels in the isostatic balance, so scaling ocean depth would be wrong. The rule is now stated once in `scaledHeightKm` instead of being triplicated |
| GRAV-1 | Decide whether to give Orogen's erosion physical units so gravity enters generation | `notes/audits/orogen-gravity.md` | wontfix -- the erosion law is n = 1, at which post-hoc 1/g scaling and correct-gravity erosion are the SAME operation, and the network is at grade nearly everywhere. Difference is 5.7% of local relief on single-cell headwaters and nothing elsewhere. Superseded by GRAV-4/5/6 |
| GRAV-2 | Re-examine basin statistics for the gravity texture error | `notes/audits/orogen-gravity.md` | wontfix -- rested on the superseded hillslope reasoning. The fluvial scaling is correct, so any basin-size effect must come from glacial overdeepening instead: GRAV-6 |
| GRAV-3 | Revisit the glacier rough pass for the slope error | `notes/glacier-rough-pass.md` | wontfix -- there is no slope error to inherit. The pass keys on relief amplitude, which the 1/g scaling gets right |
