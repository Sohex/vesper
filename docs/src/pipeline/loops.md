# 4. Why this is not a straight line

Three quantities each depend on the other two.

**Drainage depends on climate.** Which basins survive is a water balance.

**Climate depends on drainage.** Closed-basin fill is a large minority of this
planet's land and the brightest thing on it -- salt crust and playa clastics
against a much darker land mean. Carve the basins and the world gets darker.
This is also the channel through which a lithology bug reached the climate,
twice, so measure it on the surface the model actually sees rather than on the
rock table.

**Climate depends on the biosphere, and the biosphere on climate.** Bare rock
and a vegetated surface differ in land albedo by enough to be worth several
kelvin, and the two reach any given design mean at *non-overlapping* stellar
fluxes. That is the load-bearing fact: no single flux is robust to the
vegetation question, so the orbit and the biosphere are one choice, not two.

The two interact rather than adding. Vegetation paints everything that can
carry a canopy at a single value, so it masks bare-rock variation but not the
barren classes; a lithology change confined to closed-basin fill therefore
moves the *vegetated* albedo roughly twice as far as it moves the bare one.
Quote the vegetated figure when the question is what the climate will do.

**That coupling is radiative and aerodynamic, and it is not hydrological.**
Albedo, roughness and forest fraction carry the vegetation state into the
model; **what is absent is transpiration**, because the land surface is a
single bucket with no stomatal control, LAI dependence or rooting depth. The
absence is a limitation rather than a footnote: the channel would have been
worth the same order as this world's entire land runoff, runoff is the
denominator of the carve criterion, and the sign is not obvious either way.
So the biosphere-climate coupling is priced in kelvin and unpriced in
millimetres, and the carve is decided in millimetres. Evidence and the
model-source reading that establishes the absence:
`notes/audits/missed-couplings.md`, finding 4.

**Soil depends on the biosphere, and the biosphere on soil.** Texture, pH and
regolith depth are weathering products of lithology under a climate, but the
organic fraction is what the vegetation leaves behind, and it changes the bulk
density and water-holding capacity the vegetation then grows in. `pedology/`
therefore iterates against `biosphere/` rather than running once before it.

The loop is therefore: assume, compute, feed back, repeat.

**There is a fourth loop, and it is currently cut rather than closed: dust.**
Emission reads the wind and the surface, the burden sets an optical depth, the
optical depth changes the radiation, and the radiation changes the wind. The
offline chain in `aeolian/` runs that path exactly once per iteration and
stops, which is what makes it PRESCRIBED rather than interactive: the dust the
model sees is the dust the previous climate produced. That is defensible one
iteration deep and not for a converged answer -- the reopening test in
`notes/dust.md` fired. DUST-3 closes it by putting emission in the model.
Until then, `build_dust.py` sits BELOW the climatology in the register and
surface code 1811 sits ABOVE the next run, which reads as a contradiction in
the diagram and is really one loop drawn across two iterations.

**One loop is deliberately left open: the carbon cycle.** `config/planet.yaml`
fixes CO2, and nothing in this project solves the carbonate-silicate balance
that would set it -- a defensible choice for a snapshot climate, to be read as
an assumption rather than a result. The pipeline measures what the assumption
costs (a weathering rate is also a statement about the outgassing the world
needs to hold its atmosphere steady), and what exists is a bound on the size
of the assumption, not its resolution; closing the loop needs a weathering law
in pCO2 and a climate response to it, a coupled calculation across pedology
and exoplasim rather than an analysis in either.

Two things to hold on to, and one caution. The requirement is driven by **land
area, not weathering intensity** -- a big-land planet is a high-outgassing
planet or it is a cold one. Only the exorheic share joins the marine carbonate
feedback that stabilises CO2; `pedology/scripts/thermostat_efficiency.py`
measures the share. And the budget is a sum over lithologies in which a tiny
class with extreme solute chemistry can decide the total, so
`weathering_fluxes.py` reports per-class contributions: check which class is
on top before quoting the total.

**The verdict map is antitone, not monotone, and that changes what the loop
does.** Carving removes closed-basin fill, the brightest lithology, so the
land darkens, the world warms, open-water evaporation rises, and basins that
were marginal would now stay closed. A larger carve set produces a *smaller*
next verdict. An antitone map does not approach a fixed point from one side;
it oscillates, and successive verdicts bracket the answer rather than
converging onto it.

That is a better procedure than a one-sided approach, because a bracket is
measurable. Take the verdict at both bounding climates -- the cold, bright,
bare-rock end and the warm, dark, vegetated end -- and carve only the
intersection. Everything between the two is the BRACKETED set *by
construction* rather than by a tolerance chosen after the fact, and the width
of the bracket is the honest uncertainty on the carve.

**`bracketed`, not `marginal`.** `marginal` names a per-basin landform in the
carve list; `bracketed` names our uncertainty. The two were briefly one word
and it leaked into an exporter; the definitions and the third neighbour,
`disputed`, are in [the vocabulary](../reference/vocabulary.md).

**Both arms run at ONE flux, and it is the design flux.** Decided 2026-08-19.
The alternative -- each endmember at whatever flux keeps IT in the design
range -- collapses the bracket: those two worlds sit at the same global mean
by construction, so the intersection lands close to either verdict alone and
the reported width understates the real uncertainty, which something
downstream will quote. It is also the reading that matches what is actually
unknown: the orbit is chosen, not uncertain, while the vegetation state is
genuinely unknown when the verdict is taken, because the carve happens before
LPJ-GUESS has ever run and the verdict inherits an assumed biosphere.
Bracket the thing you do not know.

**The flux bracket itself runs on the VEGETATED branch only.** That is the
functional world, the one this project has chosen and the one the design flux
is defined against. Running it on both branches would be deriving two orbits
for a world that has one.

**The bare-rock arm is a BOUND, not a world.** Nobody claims Vesper sits at
the design flux without a biosphere; that run exists to produce the cold-end
verdict and nothing else. Say so wherever its numbers appear, because a reader
meeting its climatology in this repository will otherwise take it for a
description.

**Expect the bracket to come out wide, and do not read width as failure.** The
endmember spread widens as the world cools and sea ice grows back to amplify
it -- the T21 work measured it nearly doubling across a tenth of a unit of
flux -- and the bare-rock arm sits at the cold end of exactly that behaviour.
A large bracketed set is the correct answer to a genuinely uncertain question,
and if the bare arm grows glaciers where the vegetated one does not, that is a
finding about the bound rather than a reason to move the arm.

Over-carving is a budget item, not a lost landscape. A build is regenerated
from the planet code plus a verdict rather than edited, so an over-carve costs
a terrain, hydrography and boundary-condition rebuild and nothing else. That
is how `carved-zoned` was abandoned once its verdict turned out to have been
computed on antipodal climate: its carves could not be un-cut *within that
build*, and the build was replaced wholesale.

On iteration 2, the already-carved set should be re-evaluated against the new
climate and the number that would no longer have carved reported. That number
is the overshoot, and it is the honest measure of how much the first pass
cost.
