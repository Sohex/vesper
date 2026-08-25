# Four externally raised mechanisms, checked against the record

Reviewed 2026-08-20. Worldbuilding frame: Vesper is a fictional planet and
every mechanism below concerns its simulation -- a toy GCM's boundary layer, a
slab ocean, an offline dust chain, a weathering budget. Nothing here is about
the real world. Each point was checked against artifacts and code rather than
answered from the argument alone.

## 1. Boundary layer depth under 12.81 m/s2, and basin low-level jets

Two halves, and they come out differently.

**The scale-height compression is carried by the model automatically.** The
dynamics run in sigma coordinates, so a layer at fixed sigma sits at a physical
height proportional to H = kT/(mu g); at 1.31 g every level is ~24% lower in
metres than the Earth-tuned intuition expects, with no code involved. The
offline aerosol chain was separately measured insensitive to its declared
scale heights: `notes/audits/unpriced-terms.md` finds the wet-removal lifetime
carries almost all of the dust bracket "and the scale height almost none", and
each config declares a bracket (dust 2000-5000 m) that spans the compression.

**The nocturnal low-level jet is a real, one-signed, unresolvable gap.** The
Bodele mechanism -- a sharp nocturnal inversion decoupling the surface, the
jet above it mixing down after sunrise -- is exactly what
`washington2005-bodele-low-level-jet.pdf` (held) documents, and its vertical
structure is decameters to a few hundred metres.

*Measured against the model's own grid, 2026-08-24.* The bottom two sigma
levels in the baseline climatology are 0.98282 and 0.92191. Through the
hypsometric relation at this world's gravity and a 280 K near-surface
temperature that is a lowest full level near 110 m and a second near 510 m, so
the model carries two levels below a kilometre and ONE below half a kilometre.
Both the inversion and the jet core would sit inside the single lowest layer,
whose wind is that layer's mean. The obstacle is therefore not that the jet is
weakly resolved; it is that the model has no interior to the layer the jet
lives in.

Gravity neither causes this nor relieves it. Sigma levels sit at heights
proportional to the scale height, so at 1.31 g every level is about 24% lower
in metres AND the nocturnal boundary layer compresses by the same factor; the
count of levels inside it does not move. What this world's gravity does is
concentrate the same surface cooling into less air, which strengthens the
inversion and makes the regime more likely to matter, not more likely to be
resolved.

The project's existing instrument for wind-tail problems, high-cadence gust
sampling (DUST-5, DUST-15), measures the tail THE MODEL PRODUCES; it cannot
recover a jet the model cannot represent. So basin deflation is understated
one-sidedly. It is a FLOOR rather than a bracket: no end of the declared
roughness or vegetation brackets moves it, because the missing wind is absent
from every arm.

DUST-16 declares the gap; the honest revisit is vertical resolution, not a
parameterisation fitted to nothing. What would settle whether L20 is that
revisit is one number nobody here has: an L20 binary's own sigma array. Any
level rule that doubles the count while keeping its shape roughly halves the
near-surface layer thickness in log-sigma, which puts the lowest level
somewhere in 30-60 m and the second somewhere in 130-250 m; that is a bracket
and not a result, and reading the array off a built L20 executable replaces it.
Even at the optimistic end it buys an interior to the jet layer rather than a
resolved inversion.

## 2. The 50 m slab at 32 degrees obliquity and a half-length year

The mechanism is real and the project has already priced its frame: the config
note on `mixed_layer_depth_m` says the half-length year damps seasonality
about twice as hard as Earth's ocean, and CLIM-33 (closed 2026-08-20) measured
the arms. Two things blunt the concern as raised:

- **The direction, as measured.** The claim is that a 50 m slab keeps marginal
  ice frozen through a short intense summer, i.e. a shallower ocean would hold
  LESS ice. The 25 m arm grew MORE ice, +0.53% of the planet, because the
  sharpened winter swing crosses TFREEZE before the sharpened summer clears
  anything. On this world's margins the winter side wins. Two orbits is a
  transient, so the equilibrium geography is not settled, but the sign as
  measured is opposite to the one the concern needs.
- **The bracket is common-mode.** Both endmember climates carry the same slab,
  so a seasonal-amplitude bias moves both arms together and largely divides
  out of the bracket, which is a difference. What survives is second order at
  the margins.

No action beyond what CLIM-33 already records. If a future measurement shows
marginal-ice geography moving the carve verdict, the lever is priced at
1/depth and one namelist key.

## 3. NIR-shifted flux, dust absorption, and convective shutdown over basins

Already the design of DUST-11 rather than a gap in it. The prescribed-dust
prediction (`aeolian/notes/prescribed-dust-run.md`) prices the fast adjustment
at -4% to -9% on net global precipitation, carries the column-heating overlap
explicitly, and declares the falsifiers before the run. The spectral half of
the concern is in the optics already: the aerofile's per-band Q values are
built on the k25v band split (DUST-12 machinery), so the NIR shift of the
incident flux is in the numbers, and the measured single-scattering albedos
(0.964-0.975, `analysis/dust_optics.json`) are this star's, not the Sun's.
Whether the elevated heating layer shuts convection down harder than the
declared range is precisely what the run measures; the mechanism adds a
diagnostic worth reading (the vertical heating profile over arid basins), now
listed in the note, and changes no threshold after the fact.

## 4. Endorheic carbonate burial is never recycled

Correct, and the half the project can state is now stated where the outgassing
number lives. `weathering_fluxes.py` already reports the implied outgassing
(~2x Earth, from 2.1x land area) and says nothing solves the carbonate-
silicate balance. The addition: on a world burying a large share of its
carbonate IN CLOSED BASINS ON CRATON, that outgassing is not resupplied by
subducted carbonate the way Earth's arcs are fed, so it draws one-way on the
interior. As a snapshot this changes no number -- steady state balances fluxes
wherever burial happens. As secular evolution it is a duration, and
`docs/src/reference/no-time-axis.md` refuses it; the sentence exists so that
if a long-term carbon balance is ever built, the recycling asymmetry is the
first trap on the list rather than a rediscovery.
