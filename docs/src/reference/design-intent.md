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
- **An external model becomes a maintained fork, and that is the expected end
  state rather than a cost to weigh against adopting one.** These models were
  built for other people's priorities: Earth as the subject, the authors'
  hardware, and publication rather than a reproducible pipeline. This project's
  constraints differ in kind, so the question to ask about a candidate component
  is not whether it can be adopted unmodified, but whether its physics is worth
  the fork its runtime and its Earth content will require. Cost the fork, then
  decide on the physics; do not treat the fork as the objection.
  `docs/src/reference/vendored-upstreams.md` carries the mechanics, and the
  subtree layout exists so that maintaining one is ordinary work.
- **A judgment made while a class is absent is not a judgment.** Values that
  nothing exposes go unchecked; three of them surfaced at once when a single
  tectonic bug was fixed. When a rule starts firing for the first time, audit
  everything it controls.

- **This world's cloud droplet effective radius is Earth's stratiform relation,
  and that is declared rather than absent.** The shortwave cloud optics is
  Stephens (1978), whose fitted optical depths carry the effective radius
  implicitly rather than taking it as an input: his Eq. (7) is
  `tau = 1.5 W / r_e` and his p. 2125 states that the dependence has been
  "inherently parameterized" into the fits, over eight terrestrial stratiform
  cloud models. So there is no missing input to supply. What the model carries
  is a fixed `r_e(W)` measured on Earth's clouds, and nothing in this project
  determines a cloud condensation nucleus population that would give another.
  Supplying an effective radius means replacing the scheme with one that takes
  it, which is a scheme decision and not a constant.
  `exoplasim/notes/cloud-water-reference.md`.

- **The broadband radiation scheme stays, and a band-resolved replacement is
  refused on price.** Measured, not estimated: the correlated-k candidate costs
  8.6 times the present scheme per column per radiation call on the corner this
  model actually runs -- cloudy, which is every column, since `nswrcl=1` and the
  present scheme's own cost is measured in situ with that cloud in it. Radiation
  is 32 per cent of T21 wall clock and the share falls to 0.858 of that at T85,
  so the swap trebles a commissioning: 5.2-10.6 hours becomes 16-33. The
  criterion was fixed before the measurement and the cloudy corner lands past its
  prohibitive end. What the swap would have bought is three holes re-weighting
  cannot close -- the per-band attribution inside the CO2 total, wrong by two
  orders of magnitude at 2.7 um and surviving by cancellation; 13 per cent of the
  CO2 shortwave absorption falling outside every band Howard measured; and no
  water vapour continuum term, which is why `h2o_sw_level` is a bracket. Those
  three are not worth trebling every commissioning this world runs. The present
  scheme is within a few per cent of correlated-k on the comparison
  `exoplasim/notes/corrk-cross-check.md` records, which is what makes the refusal
  affordable rather than merely cheap.
  THE THREE HOLES REMAIN OPEN AND ARE NOT CLOSED BY THIS DECISION: the route to
  them is targeted gap fill inside the present scheme, priced one hole at a time,
  and NOT a scheme replacement re-argued from the same three. A future proposal
  to swap has to beat this measurement rather than restate the holes.
  `exoplasim/notes/radiation-scheme-price.md` carries the price, its bracket and
  the load it was taken under.
