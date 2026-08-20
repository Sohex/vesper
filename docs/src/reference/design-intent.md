# Design intent

- **The orbit and the biosphere are one choice, not two.** The flux windows that
  put this world in its design temperature range do not overlap between a
  vegetated surface and a bare-rock one, and the endmember spread near the target
  is wider than the target band. No single flux is robust to the vegetation
  question. Chosen: vegetated. The band is narrow, so re-derive after anything
  that moves land albedo.
- **The pipeline is a loop, not a line.** Drainage depends on climate, climate on
  drainage, and both on the biosphere. `docs/src/pipeline/loops.md` says why, and
  which loop is deliberately left open.
- **Do not reuse a sensitivity measured in one regime in another.** Bracket
  between two converged points that span the target rather than extrapolating
  from one. Doing the latter across the ice transition predicted 291.9 K for a
  run that converged at 287.47 K.
- **A correction's size depends on how much of the surface it acts on.** Estimate
  it against the state you are in, not the state it was first measured in.
- **A judgment made while a class is absent is not a judgment.** Values that
  nothing exposes go unchecked; three of them surfaced at once when a single
  tectonic bug was fixed. When a rule starts firing for the first time, audit
  everything it controls.

