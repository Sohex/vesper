# The docs-book restructure audit

Six fresh-context reviews of the restructured book, 2026-08-19, each verifying
its chapters against the code and configs; the findings and fixes are in git
(commits 5537fb3 through 67d99ec). What remains binding are the acceptances --
deliberate states a later reader might otherwise re-litigate:

- Magnitude-bearing measurements stay in prose where the argument fails
  without them: the equilibrium-line excursion over the cycle, the T21
  bracket-widening measurement, the dwmax-against-uniform contrast, the
  damping-factor range, the mesh cell scale.
- `docs/src/reference/builds.md` keeps its illustrative pair of numbers: the
  sentence exists to say the second kind does not survive a regeneration, so
  a grepping reader lands on its own warning.
- The pointer-table row for `scripts/pipeline.py` in CLAUDE.md and the
  register section of the components chapter deliberately overlap: map row
  and canonical home.
- SUMMARY numbering shows 0, 3, 4, 5, 6 with the components chapter
  unnumbered (it holds old sections 1, 2 and 2b); the introduction says so.
