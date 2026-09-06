# Offline ocean circulation

`ocean/` is the project-owned boundary around the adopted cGENIE/GOLDSTEIN
circulation host under `vendor/cgenie/`. The coupling is offline full-flux:
ExoPlaSim writes one chronological surface-forcing artifact, GOLDSTEIN advances
to an accepted equilibrium, and the ocean returns heat-flux convergence and
surface velocity to the next ExoPlaSim baseline. ExoPlaSim remains authoritative
for sea ice.

The machine-readable decision is `config/transport_loop.yaml`. It names the
handover, the controls that remain inert or bracketed, the complete provenance
chain, and every exit criterion. Run its static gate from any directory:

```bash
python ocean/scripts/transport_loop_gate.py
```

The report is `analysis/transport_loop_gate_report.json`. The gate runs no
model; it refuses when the declaration, canonical pipeline graph, planet config,
or the cGENIE source seams on which the architecture rests move apart.

The component's pipeline scripts are:

- `ocean_column_response` (`analysis/ocean_column_response.py`): OCN-1's exact
  no-run seasonal response and negative adoption verdict for the dormant
  multilayer slab path.
- `ocean_shortwave_penetration`
  (`analysis/ocean_shortwave_penetration.py`): OCN-6's energy-closing
  clear-water deposition bound and negative implementation verdict.
- `carbon_feedback_gate.py`: OCN-16's fail-closed prescribed-CO2 boundary and
  conservative reopening contract; it activates no carbon cycle.
- `transport_loop_gate.py`: OCN-5's declaration and source/graph gate.
- `build_spatial_support.py`: OCN-11's wet mask, volume, connectivity,
  bathymetry and atmosphere/ocean crossing.
- `build_forcing.py`: OCN-10's chronological heat, momentum, freshwater, salt
  and ice forcing artifact.
- `run_cgenie.py`: one UUID-named offline circulation integration.
- `assess_convergence.py`: the ocean integration's own equilibrium and budget
  verdict.
- `build_transport_return.py`: the heat-convergence and surface-velocity return
  fields, conservatively mapped onto the atmosphere support.
- `assess_transport_loop.py`: the across-pass OCN-5 exit predicate.

The GOLDSTEIN barotropic direct solve is intentionally retained. OCN-20's
pre-registered replacement threshold was a serial share above 20%; profiling
measured 1.83% at 36 x 36 x 16 and 2.57% at 72 x 72 x 16. The decision is
recorded in `notes/audits/cgenie-parallelism-and-coupling-support.md`; it is
reopened only for an accepted support above 144 x 144 or a fresh profile above
that threshold.

`vendor/cgenie` carries a second, unreachable surflux family that is kept
deliberately: `surf_ocn_sic.F`, `outm_surf_ocn_sic.f`, `gold_ocnsic_avg.F` and
`surf_ocn_sic_wrapper` compile into every executable and no line of `genie.F`
calls any of them. They are the coupling to an atmosphere this subtree no longer
contains, they compute surface fluxes from state on the ocean grid rather than
taking them, and they are NOT the driven path. That path is `flag_fluxatmos`,
which calls `surflux_goldstein_seaice` and refuses in `initialise_genie.F` while
the supply half is missing. `notes/audits/cgenie-unreachable-code-dispositions.md`
is the verdict, the evidence, and the rule that decides whether unreachable code
in this maintained fork is deleted or kept.

`exoplasim/scripts/verify_ocean_flux_channel.py` remains the independent OCN-2
instrument for surface code 903. Its published invocation for the pipeline row
uses `check ... --output exoplasim/analysis/ocean_flux_channel_report.json`.

`analysis/ocean_remap.py` owns the exact Gaussian-to-GOLDSTEIN crossing. Its
JSON report embeds the SPAT-1 contracts and SPAT-10 assessment, while its NPZ
weights carry the three contract identities. Before OCN-11 supplies wet volume
and topology, the GOLDSTEIN side is deliberately a candidate comparison
support rather than an accepted ocean support.

## Wiring a cgenie component means declaring its calibration state first

cGENIE carries machinery whose PURPOSE is to rescale modelled fields until they
match Earth observations, and every part of it is compiled into the executable
whether or not the recipe runs it. So before any cgenie component is wired to
anything here, the decision states which calibration switches are off, which
factors are at their identity, and which reference files are absent -- rather
than the answer being discovered from a result.

The declarations made so far are in `notes/audits/tuned-values.md` section 19,
which also carries what they found. Two of them are the reason this is a rule
and not a formality: `genie-rokgem`'s four 2-D weathering factors ship OFF their
identity behind a single logical, and `genie-embm`'s Atlantic-Pacific freshwater
adjustment ships NONZERO and is armed by `flag_ebatmos`, which reads as a choice
of atmosphere. The adopted path is EMBM-free and reaches neither, which is why
both are preconditions rather than defects in the world's numbers.

Generated data are namespaced under `data/<source_build>/`. Ocean runs are UUID
named under `runs/`; `runs/INDEX.json` is their identity record. No script in
this component silently selects another component's latest artifact.
