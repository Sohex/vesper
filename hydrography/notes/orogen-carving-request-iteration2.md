# Iteration-2 carve request for World Orogen

Ready to send. The Orogen session closed before this could be relayed, so it is
recorded here rather than lost.

## The list

    hydrography/data/carved-zoned-v4/carve_list.txt   (JSON sidecar alongside)

Same format and same catalogue as iteration 1: 3,629 entries, one line per basin,
`<basin_id>  <retain>`, retain in [0,1], carved basins listed explicitly at 0.

| retain | basins | |
| --- | ---: | --- |
| 0.0 | 1,838 | 1,089 carried forward from the pass that produced v4, 749 new |
| 0 < r < 1 | 47 | marginal |
| 1.0 | 1,744 | comfortably closed |

The 1,089 carried entries are basins v4 already carved. They are not re-decided:
the terrain that would justify re-examining them no longer exists, so they are
held at 0 to keep the loop monotone. Verified 3,629 unique ids, no duplicates,
every carried entry exactly 0.

Start from the pre-carve terrain, as with v4. Orogen regenerates from the planet
code plus the list in one pass, so a build is replaced rather than edited.

## Why v4 is under-carved

The verdict integrates catchment runoff and was reading ExoPlaSim's `mrro`. That
is not local runoff generation: `roffstep` calls `mkradv`, which advects runoff
downhill and modifies its argument in place, so `mrro` is local generation minus
river outflow plus river inflow. A net divergence, on the model's own grid and
its own downhill directions, which know nothing about our basins.

The right quantity is precipitation minus evaporation over the catchment.
Measured on v4's 2,540 remaining basins against one climatology, decided both
ways:

| | `mrro` | P - E |
| --- | ---: | ---: |
| carve | 176 | 749 |
| marginal | 90 | 47 |
| preserve | 2,274 | 1,744 |

591 basins that `mrro` preserves, P - E carves. Verdicts agree on 71.3% and
retain fractions correlate at 0.281, about as decorrelated as the longitude bug
was at 0.267. `mrro` reported no catchment runoff at all for 68.4% of basins,
which traces to a uniform soil-bucket depth that almost never overflows: a land
runoff ratio of 2.8% where Earth manages about 35%.

## The failure worth naming

Three consumers of catchment runoff. `pedology` switched to P - E and documented
why in `pedogenesis.yaml`. `surface_water.py` switched to P - E and documented
why in its own docstring. `export_carve_list.py` was never updated, and it is the
only one of the three whose output leaves the project and changes the terrain.

Two independent corrections of the same misunderstanding, neither of which
propagated to the consumer that mattered most. Worth a check across consumers
whenever a shared quantity is redefined, rather than fixing it where it was
noticed.

## What to ask for back

The usual export set with `raw/` on T42, and `manifest.hashes.finalElevation` so
the build can be registered in `lib/orogen.py`. Also worth asking: what
closed-basin fill comes out as. It went 20.9% of land pre-carve to 16.5% on v4,
749 more basins lose their rims here, and that number drives land albedo directly.
