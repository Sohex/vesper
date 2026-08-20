# The docs-book restructure audit

Measured 2026-08-19, on the docs-book branch, immediately after the split of
WORKFLOW.md and CLAUDE.md into the book under `docs/src/`. Six independent
fresh-context reviewers, one per chapter group plus one repo-wide reference
sweep, each verifying claims against the code and configs rather than against
other prose. This is the natural moment for such a pass: every chapter was
moved, so every chapter was read. This document records what was found and
what was done about it; the fixes landed in the commit that adds this file.

## Method

Scopes: pipeline chapters 0-3; pipeline chapters 4-6; introduction, SUMMARY,
vocabulary, builds, design-intent; vendored-upstreams, environment, and the
five topic chapters; CLAUDE.md, README and the practice chapters; and a
repo-wide sweep of every documentation reference outside `docs/src/`. Each
reviewer verified script paths and flags against argparse, config keys against
the YAML, and cross-references against the tree. `smoke_test.py` passed (15
checks) before and after the fixes.

## Fixed: content that was wrong before the restructure

The split did not create these; reading everything surfaced them.

- The components chapter's data-flow diagrams named eleven generator scripts
  at `<component>/<script>` where every one lives at
  `<component>/scripts/<script>`; cited `carve_verdict.py --dust`, where the
  flag is `--dust-forcing`; omitted the T63 export and placed `maps` inside
  the per-build directory; and drew the dust optical path straight into the
  model, skipping the `surface_dust` step the register routes it through.
- The loops chapter presented DUST-3 as future work; it landed 2026-08-18,
  and what remains is enabling emission in a run (DUST-13).
- The sequencing chapter said the cycle-before-verdict ordering "is now CYC-1
  as well"; CYC-1's own row disclaims that, and the ordering is enforced by
  the graph (`carve_verdict` needs `stellar_cycle_run`). The same chapter
  called over-carving "the irreversible direction" while arguing elsewhere
  that a wrong verdict is recoverable; it now calls it the expensive
  direction. It also cited `aeolian/analysis/dust_runoff_sensitivity.json`,
  which no step generates and which is not on disk; it now cites the
  registered one-off script.
- `config-rationale.md` quoted `soil_water_source: uniform` where the config
  says `pedology`; carried the superseded single-sinusoid stellar-cycle era
  stacked above its two-component replacement, against the
  rewrite-superseded-content convention; used a status (CANDIDATE) the
  four-status system does not define, now PROVISIONAL with the nuance stated;
  called the baseline flux DETERMINED two hundred lines after recording its
  demotion to PROVISIONAL; and pointed three notes at `notes/` paths that
  live under `exoplasim/notes/`.
- `economic-minerals.md` said `r_t_craton` has no downstream reader
  (`build_prospectivity.py` reads it), keyed VMS on `backArcDist` (absent
  from minerals/; the hosts are melange, arc basalt, morb), and keyed
  magmatic Ni-Cu on `lipV` and on `flood_basalt` "alone" (lipV is expressed
  through the flood-basalt class, and oib carries a small weight).
- `no-time-axis.md` stated the incision coefficient as a fixed 161; it is
  solved by calibration on every run (`hydrography/notes/retain-fraction.md`),
  and the last recalibration left the old figure's bracket.
- `external-data.md` copied the star's log g and the atmosphere's CO2 into
  prose against rule 2; both now point at the deriving script and the config.
- `pedology/config/pedogenesis.yaml` and `build_soil.py` pointed at
  `notes/model.md`, which has never existed in git history.
- `smoke_test.py`'s generator-check docstring said the check reads the
  workflow document; it reads `config/pipeline.yaml`.
- Supersession markers that re-supplied the superseded value ("not 2:1 as an
  earlier revision said"; "Corrected 2026-08-19; it used to say three") were
  rewritten to state only the current claim.
- The run-count paragraph claimed "three of the six pair naturally and one
  cannot" while listing pairings that shared a member and an arithmetic that
  partitioned nothing. The design intent was never ambiguous -- the
  dependency edges are all documented -- so the paragraph now states the
  edges (bracket points need only the build; each arm's baseline needs its
  own bootstrap climatology; the cycle follows the final vegetated baseline,
  the bare arm gating nothing as a bound) and lets the schedule fall out.
  The seeding sentence now permits any converged vegetated restart, which
  the same paragraph's own --restart-from argument already licensed.
- Current values in prose, removed or repointed: the per-orbit wall-clock
  figures (the pyburn audit note keeps them), the model configuration
  transcribed from `config/planet.yaml`, the obliquity figure, the playa land
  share (which had already drifted from the config's own comment), the open
  task count, and the dated derivation story in section 5b.
- The early failure-mode classes (1-7) carried em dashes and non-ASCII math
  glyphs from before the punctuation convention; normalized.

## Fixed: seams the restructure itself made

- Deictic references that moved out from under their antecedents: "rule 7
  above" in the vocabulary; "the rule below" three times in
  working-agreements, now naming the defect-as-decision rule; "section 7"
  citations for rules that live in `docs/src/practice/conventions.md`;
  bare "section 4" cross-file mentions, now links.
- The compressed "re-baseline" ban in CLAUDE.md inverted the full
  definition's meaning (means-neither vs ambiguous-between); corrected.
- Three conventions lost their one-line directive in CLAUDE.md
  (thresholds-fixed, estimates-bracketed, claims-checked); restored.
- The simulation's-terms rule was the one convention still argued in place;
  the argument moved to the conventions chapter (citing the three commits by
  hash), the directive stayed.
- Assembly duplications in `builds.md` and `vendored-upstreams.md` trimmed;
  the disposable-until-consumed argument rule 7 points at was missing from
  `builds.md` and was added.
- One mechanical rewrite mapped a reference to the wrong chapter
  (`config/pipeline.yaml`'s design-flux gate pointed section 6 material at
  the state chapter); corrected to the sequencing chapter.

## Accepted, with reasons

- Magnitude-bearing measurements stay where the argument fails without them:
  the equilibrium-line excursion over the cycle, the T21 bracket-widening
  measurement, the dwmax-against-uniform contrast, the factor-of-2.5 damping
  range (consistent with CYC-1's 0.24-0.6), the mesh cell scale (stated once,
  approximate).
- `builds.md`'s illustrative pair of numbers ("binding on 73%", "1,938 basins
  carve") stays: the sentence exists to say the second kind does not survive,
  so a grepping reader lands on its own warning.
- The pointer-table row for `scripts/pipeline.py` in CLAUDE.md and the
  register section of the components chapter deliberately overlap: map row
  and canonical home.
- SUMMARY numbering shows 0, 3, 4, 5, 6 with the components chapter
  unnumbered (it holds old sections 1, 2 and 2b); the introduction says so.

## Open

- The ExoPlaSim fork carries a remote branch `cf-fork` identical to master --
  exactly the ownerless copy `vendored-upstreams.md` argues against. Deleting
  a remote branch is the repo owner's call, not a docs fix.
- `check_registered_in_workflow` in `smoke_test.py` keeps its pre-split name;
  rename whenever the function is next touched.
