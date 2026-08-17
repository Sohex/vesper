# Aeolian

Mineral dust: emission, transport and deposition, computed offline from a
climatology rather than inside the climate model. `notes/dust.md` carries the
decision and why it is not in the GCM; this component is the implementation.

```bash
python aeolian/scripts/build_dust.py                          # the baseline source map
python aeolian/scripts/build_dust.py --variant arid_bare_ground   # the bracket on vegetation
```

Products are `analysis/dust_<variant>.json` and a matching `.nc` carrying annual
mean optical depth, deposition flux, emission and the erodible fraction.

## STATUS: THE MACHINERY WORKS AND THE NUMBERS DO NOT

Read this before quoting anything from this component.

The model runs end to end, every constant is grounded in a primary source, and
the provenance is stamped. **Its emission is three to four orders of magnitude
above Earth's**, which is about 2000 Tg per year: this build produces 3.4e6 Tg
per Earth year at the central roughness. The optical depths that follow are
physically impossible, around 50 rather than around 0.05.

So the reopening test in `notes/dust.md` is reported as CROSSED and that result
must not be used. A number that is four orders out does not answer a question
whose threshold is a factor of two.

Two causes are known and one is not.

**Known and fixed: the subgrid wind shape was declared and is measurable.**
Emission goes as roughly u* cubed above a threshold, and the mean u* over source
cells on this world sits BELOW the threshold, so every gram of dust comes from
the tail of the distribution. The shape parameter was declared at k = 2 with a
bracket of 1.5 to 3.0. Measured from the snapshot climatology's 32 instantaneous
samples it is **3.98**, outside that bracket, and worth two orders of magnitude
on the answer by itself. `weibull_shape_from_snapshots` now measures it rather
than taking the declared value. That is a lower bound on emission, because
samples 5.7 days apart resolve synoptic but not sub-daily variance.

**Known and not fixed: transport does not reach steady state.** The advection
solver reports `transport_converged: false` at 4000 iterations for every case,
so the loads are still growing when it stops and the optical depths are upper
bounds of an unconverged integration.

**Not diagnosed: the emission magnitude itself.** With the measured wind shape
the duty cycle above threshold is a few percent, and a hand calculation of Kok
equation 18 at representative values gives an annual mean around 4e-8 kg/m2/s,
which over this world's source area would be of order 1e5 Tg per year. The model
returns fifty times more than that. The discrepancy is somewhere in the u* chain
or in the flux assembly and it was not found.

## The structural finding, which matters more than the bug

**A 12-bin climatology has averaged away exactly the variance that drives dust
emission.** Because the mean friction velocity over the source cells lies below
the threshold, the emission is entirely a property of the wind distribution's
tail, and the tail is not in the driver. Measuring it from snapshots recovers
synoptic variance and still misses the diurnal and sub-daily part.

`notes/dust.md` framed the reopening test around feedback magnitude: run it in
the GCM if the forcing is large enough that a prescribed field stops being
defensible. This is a different failure mode and the note does not anticipate
it. The offline route is compromised here not because the feedback is missing
but because **the driver has lost the information the emission depends on**. An
in-model scheme sees every timestep's wind and does not have this problem at
all.

That is a decision for the pipeline rather than for this component, and it is
recorded here rather than resolved.

## What the component does get right

**The source map, which is the thing the built-in GCM scheme gets wrong.**
`fcoeff * land_mask` emits as much from forest as from salt pan. This one starts
from `substrate_class`, which is consolidated lithology plus closed-basin fill,
so the only unconsolidated material is the fill; weights the two barren classes
separately, because a cemented salt crust is not a silicate soil and Kok's
fragmentation theory does not describe halite cement; removes standing water
from the solved lake extent; and removes snow. It gives **16.4% of land** as
bare erodible ground against a playa fraction of 23.9%, and the difference is
lakes and crust rather than an assumption.

**The physics is grounded rather than recited.** Every constant traces to a
fetched primary source, listed in `config/dust.yaml` with what it is worth:
Kok et al. (2014) for emission, Kok (2011) for the emitted size distribution,
Fecan et al. (1999) for the soil-moisture threshold, Marticorena and Bergametti
(1995) for the drag partition, Sportisse (2007) for below-cloud scavenging.

**Settling uses this world's gravity.** At 12.81 m/s2 a given particle falls
1.31x faster than terrestrial intuition, so dust lifetime here is shorter than
Earth analogues suggest. That is in `settling_velocity` rather than in a comment.

## What is declared rather than derived, and where the brackets are

`config/dust.yaml` carries all of it. The widest are, in order:

1. **The subgrid wind shape**, now measured but only to an upper bound.
2. **The aeolian roughness of the erodible surface**, bracketed 3e-6 to 1e-3 m,
   worth a factor of 40 in emission across that range. The grid-cell roughness
   field is deliberately NOT used: its median over source cells is 0.49 m, and
   feeding a 15 km orographic variance to a scheme built for centimetre-scale
   roughness elements returns zero emission everywhere. Sheltering of a patch by
   the terrain around it is therefore not represented, which biases emission up.
3. **The evaporite erodible weight**, 0.1 with a bracket of 0.0 to 0.3, for
   crust cementation. Declared suppression, not measured efficiency.
4. **Vegetation cover**, which is not modelled at all. Non-barren land is
   assumed to carry a canopy and not emit, following the project's existing
   declared position, which is generous on a world whose median land runoff is a
   few mm per Earth year. `--variant arid_bare_ground` brackets it.

## What it does not do

No dust-climate feedback: the climatology is an input and does not respond.
No vertical structure: a well-mixed column of declared scale height advected by
a single steering wind. No inter-bin microphysics, which is correct for mineral
dust because it neither coagulates nor grows appreciably.
