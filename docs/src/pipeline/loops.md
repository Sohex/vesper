# 4. Why this is not a straight line

## What the loop is FOR, which decides what each pass may read

**Bootstrapping and looping are the process of moving from the least determined
state to the most determined one.** That is not a description of the diagram; it
is the invariant every step is held to. A step reads the MOST DETERMINED input
available to it at the point it runs, and not the one its first pass happened to
have.

The corollary is what catches defects: an input pinned to an earlier stage's
artifact, when a later and better determined one exists, is holding the loop
back. It does not converge more slowly -- it converges to the wrong place, and
it does so silently, because the artifact it reads is a real artifact of a real
world and nothing about it looks wrong.

**The first pass is where this is confused with an ordering constraint, and the
two are different.** A derived field is an INPUT to the baseline run, so on the
first pass through a build there is nothing better than the bootstrap and the
bootstrap is correctly what it reads. That is a statement about what EXISTS, not
about what the step wants. On the second pass a baseline exists, and a step
still reading the bootstrap is no longer choosing the best available -- it is
choosing the first available.

**Where the invariant does NOT bite**, said so it is not applied by reflex: a
quantity that depends on the model calendar or on geometry rather than on the
climate STATE is equally determined at either stage.
`biosphere/scripts/build_vesper_header.py` fits a solstice offset against solar
declination, and the phase is a property of the calendar, so a baseline
climatology would tell it nothing the bootstrap does not. Reading the earlier
artifact there is not a violation; it is the same answer.

**The case that established this** is in `notes/audits/design-flux-two-point-response.md`.
`build_soil.py` took its runoff from the bootstrap climatology, which is
terrain-only by definition and therefore carries no lakes on any iteration,
while the same step took its vegetation from `lpj_run`, which reaches the
baseline climatology through `lpj_driver`. One step, two stages of one world,
and the split was a leftover rather than a decision. What it cost is the arid
tail of the soil: at most 8.2 per cent of this planet is inland open water, it
sits under the 74.4 per cent of land that drains internally, and that is exactly
where the thin soil is.

**How the invariant is enforced.** `lib/paths.py:best_available_climatology` is
the resolver a state-dependent step calls: it returns the baseline once
`config/planet.yaml` names one and the bootstrap before that, and it returns
WHICH, so the stage reaches the product's provenance as `climatology_stage` and
a first-pass artifact is tellable from a later one without re-deriving it. The
graph edge does not move with it. A step still declares `needs:
bootstrap_climatology` in `config/pipeline.yaml`, because that edge says which
artifact must EXIST before the step can run and the bootstrap is the one that
must: these steps run on a first pass, when no baseline exists at all. The edge
is the ordering constraint and the resolver is what the step READS, and
`scripts/smoke_test.py:check_climatology_needs_match_call_sites` holds each step
to the resolvers its declaration permits.

**This is not the silent fallback the no-fallback rule forbids**, and the
difference is not a matter of degree. That rule exists because a fallback
returns a plausible number computed from a DIFFERENT WORLD instead of an error:
the resolver default it replaced named pre-carve terrain under a superseded
spectrum. The bootstrap and the baseline are the SAME world at two stages of
determination, on one build and one terrain hash, checked by `require_build` at
every call site. Three things keep it that way and all three are load-bearing:
the choice is made on what config DECLARES rather than on what happens to be on
disk, so a named baseline that is missing raises instead of degrading quietly;
the choice is stamped on the product; and with neither key named the resolver
raises rather than reaching for a third option.


Four quantities are pairwise coupled: drainage with climate, climate with
the biosphere, and soil with the biosphere. The biosphere does not reach
drainage directly -- the channel that would couple them, transpiration, is
absent from the model, and the absence is priced below.

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
`notes/dust.md` fired. The loop stays cut until a run enables emission
(DUST-13). Until then, `build_dust.py` sits BELOW
the climatology in the register and
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
construction* rather than by a tolerance chosen after the fact.

**That bracket is over the vegetation state, and it does not cover the
antitone feedback.** The two arms differ in `land_albedo_source` and in
nothing else; both run on pre-carve terrain, so both are cold relative to the
world their own carve produces, and both cut more than that world would. The
set algebra makes it exact: an overshooting basin is one the intersection cut,
so both arms cut it, so it is never in the set the two arms disagree about.
The bracket width therefore carries no information about the overshoot, and
quoting the width alone as the uncertainty on the carve understates it. The
intersection is still the right set to carve, being the smallest defensible
one and therefore the one with the smallest feedback; what changes is that the
honest uncertainty is the width AND the overshoot, reported side by side in
basins. `notes/audits/carve-overshoot.md` has the argument, the route that
isolates the feedback, and the three ways of taking the number wrong.

**`bracketed`, not `marginal`.** `marginal` names a per-basin landform in the
carve list; `bracketed` names our uncertainty. The definitions and the third
neighbour, `disputed`, are in [the vocabulary](../reference/vocabulary.md).

**Both arms run at ONE flux, and it is the design flux.** Decided.
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

Once the carved build has its own baseline, the already-carved set is
re-evaluated against the climate the carve produced and the number that would
no longer have carved is reported. That number is the overshoot, and it is the
honest measure of how much the pass cost.

It is taken by holding the GEOMETRY still and moving only the CLIMATE:
re-take the verdict on the PRE-CARVE build's own basins and coupling matrix,
under the CARVED build's baseline climatology. Re-running it on the carved
build measures the next pass instead, because a carved basin's rim has been
breached and the finished terrain no longer holds the impoundment the
criterion is about. `hydrography/scripts/carve_overshoot.py` is the
comparison; it refuses a re-evaluation that does not cover every applied carve
rather than scoring an unknown as a zero.

## The finalizer: individually converged is not jointly converged

**A sequence of loops that each reached their own exit is not the same thing as
a system at a joint fixed point, and nothing in the ladder establishes the
second.** The loops are nested and D declares that it advances only after
"replaying A, B and C" on the new support, which is what would make the final
state a joint fixed point. In practice that replay is not bought: the cost
argument the ladder is run under puts the bulk of the work at T21 and the
minimum sufficient at the rungs above it, so the final state has terrain carved
at T21, soil and vegetation converged at T21, and a climate merely SETTLED at
T85. Each loop exited. None of them exited against the others' final state.

That is a defensible trade and it is not a defect. What is a defect is leaving
it unmeasured, because the failure is silent: every artifact is a real artifact
of a real world, and a verdict taken on a coarser climate looks exactly like one
taken on the finer.

**So the finalizer is a VERIFICATION rather than another iteration.** At the
final state, re-evaluate each loop's own exit predicate and record whether it
still holds:

| loop | what the finalizer re-evaluates |
| --- | --- |
| A | the carve intersection. Re-take the verdict at the two bounding climates on the operating support's baseline climatology, and require the intersection to be the set that was actually carved |
| B | `pedology/config/pedogenesis.yaml`'s criteria, against the soil the final state carries |
| C | its exit already IS a re-take -- the verdict on modelled vegetation, and whether basins flip. The finalizer generalises it to the operating support rather than the support it was first taken on |
| D | the route's own invariants: the operating support is reached, and every rung change happened at constant dt |

**Why verification and not one more turn of each loop.** Another iteration pays
for the three couplings that had already closed and does not say which one had
not. Verification localises the failure to a named loop, and for B, C and D it
is cheap: no new climate run, because the operating support's baseline already
exists.

**A IS NOT CHEAP AND THIS SECTION FIRST CLAIMED IT WAS.** Loop A's exit is the
intersection of the verdicts taken at the TWO BOUNDING CLIMATES, and the second
of those is `endmember_bootstrap_run` plus `endmember_baseline_run`, both
`cost: hours`, neither of which exists at the operating support. So re-taking
A's predicate there costs a pair of runs rather than minutes. The loop whose
non-closure is the most expensive to REPAIR is also the only one expensive to
DETECT, and that is the opposite of the ordering a cheap verification would
want. `scripts/verify_joint_convergence.py` reports A as NOT EVALUABLE until
those runs exist rather than substituting the single climate it does have,
because a bounding pair with one arm missing is not the construction loop A
exits on.

**What it can and cannot do.** It DETECTS non-closure and cannot repair it. A
failure on A at the operating support means re-entering loop A there, which is
the expensive thing the ladder's cost argument exists to avoid. The finalizer's
value is that it turns that cost into a measured decision rather than a
discovery, and that a passing run is EVIDENCE of a joint fixed point instead of
an assumption inherited from the nesting.

**It is `scripts/verify_joint_convergence.py`**, registered in
`config/pipeline.yaml`'s `checks` block. It reports NOT EVALUABLE as a third
verdict distinct from pass and fail, exits 2 when nothing failed and something
could not be tested, and never substitutes a coarser artifact for the operating
support's: a finalizer that passes on what it could not test is the failure it
exists to prevent.

**It runs once, after the ladder, and it is not a loop itself.** Nothing
iterates on its verdict automatically: a failure names a loop and the decision
to re-enter is the author's, because re-entering A at the operating support is a
commissioning-scale purchase.
