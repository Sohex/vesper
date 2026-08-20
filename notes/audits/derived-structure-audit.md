# The derived-structure prose audit

Measured 2026-08-19, immediately after the docs-book restructure and its
audit. The trigger was one paragraph: the run-count text in
`docs/src/pipeline/sequencing.md` claimed "three of the six pair naturally and
one cannot" while listing pairings that shared a member, and the fix was to
state the dependency edges and let the schedule fall out. This audit swept the
whole documentation corpus for the same failure class.

## The class

Derived-structure prose: a passage asserts a derived summary -- a count, a
partition, a pairing, an exhaustive enumeration, an exception claim, a
self-referential qualification -- where (a) the summary is internally
incoherent (members overlap, counts do not add, the exception contradicts the
rule), (b) it contradicts primitive facts stated elsewhere or checkable in
code, or (c) it stands in place of the primitive facts, so a reader cannot
reconstruct it. The fix pattern is always the same: state the primitives, keep
a derived claim only when it checks out and adds decision value.

## Method

Six fresh-context reviewers, calibrated on the run-count exemplar, each over
one scope: pipeline chapters; reference chapters; CLAUDE.md, README and the
practice chapters; component READMEs; exoplasim notes plus the patches README;
the remaining component notes. Every count recomputed, every partition's
members listed, every exception traced; primitive facts verified against the
code, the configs, the patch files, git history, and in one case the netCDF
itself. `TASKS.md` and `notes/audits/` were out of scope as tracker rows and
dated records. Fifty-one findings; the fixes landed with this file. Fixes were
applied only where the correct statement is derivable from the same document's
own primitives or from verified code; where the truth needed re-deriving a
measurement, the finding was first recorded open and then settled from
surviving artifacts, git history, or by dated amendment -- the section below
says how, item by item.

## Fixed

Counts and partitions that did not hold:
- "Three quantities each depend on the other two" headed a section that states
  edges for two of the three pairs and denies the third; it now names the
  four pairwise couplings and the absent transpiration channel (loops).
- "Three of the surface fields... everything else is a pure function of the
  terrain" overlapped two members and omitted the conditional dust fields;
  "seven surface fields are supplied" now separates the seven unconditional
  codes from the three key-gated ones (sequencing, steps).
- "Both prospectivity fields" sat in the worthless-below-the-climatology list;
  the tectonic half is terrain-only by the register (sequencing).
- "Two branches... would turn one picture into four" undercounted the undrawn
  branches; the register carries dust, sea salt, volcanic sulfate and minerals
  (components).
- The durable set called itself "small and complete" and left out the build
  recipe it requires two lines later; the recipe is now a member (CLAUDE.md
  rule 7, builds chapter). "Rules 5 and 7 both cite modules from this list"
  survived rule 7's rewrite; only rule 5 does.
- failure-modes: "None announced itself" against three of its own loud
  classes; "nine times wider" mixing a half-width with a full span (4.5x);
  "the last two were wrong the same way" pairing two different lessons; the
  Penman-over-ocean comparison offered as a productive cross-check where its
  own class 17 rules it a non-check; the endorheic-share section quoting two
  stale percentages and a ratio against world_state; "cost more wall clock
  than any defect in this file" against three costlier classes in the file.
- The four-status system: both `config/planet.yaml`'s header and the
  config-rationale chapter claimed exactly four statuses while both files use
  MIXED, PARTLY DETERMINED, TRANSITIVE and DECLARED ABSENT; the compound
  labels are now defined in both places.
- no-time-axis counted the expected-value move "used three times" where one of
  the three is the move's identified limit (an absolute-age control has no
  stationary population); and stated the incision coefficient as a fixed 161
  where it is re-solved by calibration on every run.
- economic-minerals told a one-key-removal story git contradicts (borate's
  host set was replaced wholesale, three hosts added), repeated in the
  config's own comment; both now say so.
- README.md claimed mass, year and flux are "derived rather than stated";
  mass and flux are declared in the config, only the year is derived.
- Component READMEs: the withdrawn 13% endorheic comparator re-quoted after
  its own correction (about one-fifth); a pointer at a
  `probe_runoff_response.py` that does not exist in hydrography; "derived
  surface classes are designed but not built" about a component the same file
  documents as running; 18.8% for a ratio the file's own table gives as
  17.8%; a 14.6x "bracket" pedology had already corrected to a 3.41x
  sensitivity; "one paper closed three gaps" where the config attributes the
  relief bound to a second paper; the dust wind-tail correction quoted at
  8.8x against the 40x its own primitives give; a stale opt-in account of the
  I/O regime contradicting the code and the paragraph beside it; "105 fields"
  where the file's own rule says read counts from the manifest.
- Component notes: 1,549 for a both-estimates carve set the same note proves
  is 1,522; "splits in two" over classes that overlap by about a hundred
  rain-fed zero-runoff lakes; the gravity note's invalidation list and status
  step standing unrewritten below the correction that refutes them; the
  runoff "currency" attributed to the denominator against the file's own
  20x-smaller denominator-only measurement; "one was dropped" about a surface
  class that ships; a 10% caveat attached to a number (0.427) the note never
  adopted (0.396 is adopted); in-model-dust's two stacked patch-list eras (4
  patches vs 3, "none touches radmod.f90" against its own table) merged into
  one five-patch account with touched files verified from the patch files,
  and its ordering list that ran 1,2,3,4,5,4,5,6 merged to one entry per
  item at its current state; "of order nine kelvin" computed on the
  sensitivity BUDG-4 rules out (4 to 5 K on the canonical slope); "three
  levers" over four; a fill-share ranking ("second, ahead of 18.9%") that its
  own numbers place third; "root-sum-square is about 2.9 K" where the RSS is
  4.07 and 2.9 is the RMS; "there is no T42 8-rank point" against the
  recorded old-world sweep; "nine of them carry headers" where the count is
  twelve of eighteen with six carrying none; "three of the four channels
  close" over a table that lists three channels and never enumerates the
  four.

## Open, then settled the same day

All eight items below were initially recorded open; each was then taken to a
conclusion on 2026-08-19, by recomputation from surviving artifacts, by git
archaeology, or by amendment where numbers are pre-registered.

- The two biosphere cost figures were BOTH real measurements in
  `notes/lpj-guess-porting-audit.md`: 24 s/cell is the shipped demo
  configuration, 5.46 s/cell the measured `npatch 5` Vesper configuration
  (times 4,106 cells = the table's 6.2 CPU-h). The README's Cost section had
  quoted the demo rate for the pipeline's own runs; it now cites the measured
  one and names the difference.
- The "three of six pass" convergence tally: the six T21 bracket runs predate
  UUID naming and their assessments were never committed, so the tally is
  unreconstructable. The prose now says so, and states what the tabulated
  climatology means support (one of six over the |mean TOA| < 0.5 criterion
  then in force).
- water-and-energy-closure: the melt-booking lines were computed on clean
  orbit 66 while the table beside them is the climatology -- confirmed by
  arithmetic (-0.3539 - 0.4222 = -0.7761) -- and are now labelled, with the
  climatology's own identities (-0.4097, -0.9327) stated from its table. The
  ice/no-ice storage column is now stated as indicative (its estimator spread
  exceeds its disagreement with the total; the run is archived, so it cannot
  be recomputed). "What is left" now matches the re-measured weighting and
  the Status section: nothing of the -0.292 is attributed, and the spread
  argument is restated against 0.29 rather than the refuted remainder.
- forcing-bundle-predictions: amended in place with a dated note -- the rows'
  brackets sum to +0.35 to +1.37, and the headline's -0.3 endpoint carries an
  implicit -0.65 K PHYS-9 downside no row states. The registered numbers
  stand; the A/B scores against them.
- The economic-minerals "9.62% -> 3.72%" sentence transcribed the three
  columns of one report row (% land, % above 0.5, mean -- the format
  `build_downstream_prospectivity.py` prints); it now names them.
- The lake-solver quartile table was recomputed from the registered JSON:
  equal-count quartiles are 37/36/36/36 of 145, first row identical to the
  published one, three rows' statistics shifted slightly.
- The `surface_classes` path: the artifact is provenance-stamped
  (`vesper_source_build`), so the register's `pedology/analysis/` location
  needs no per-build namespacing; the script's default output now matches the
  register, and the pedology README follows.
- The dust fill-share chain resolves by naming builds: 26.6% on the uncarved
  base, 16.5% on `carved-zoned-v4` (where the lake solution exists;
  `parameter-decisions.md` dates that figure), 12.4% on `carved-zoned-v5`;
  the note now attributes each number.

## Accepted

failure-modes' "those are the four that work" reads as exhaustive over a table
one row of which fits only loosely; the row can be read under "a quantity the
other side already knows", so it stands. Deliberately marked tensions (rules
that pull in opposite directions), assertive emphasis, and coherent long
arguments were out of scope by design and generated no findings.
