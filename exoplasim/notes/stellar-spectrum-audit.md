# The k2 stellar spectrum is not a K2 star

Every completed run in this project has used `exoplasim/stellarspectra/k2.dat` on
the understanding, recorded in `config/planet.yaml` and `CLAUDE.md`, that it is a
measured K2 dwarf spectrum matching this world's declared K2.5V host. It is not.
It is an early-to-mid M dwarf, around 3000 to 3450 K, and on the evidence below
it is almost certainly K2-18, an M2.5V star. The filename means the star K2-18,
not the spectral type K2.

This was found while deriving a photosynthetically-active-radiation fraction for
LPJ-GUESS, which is how a biosphere question turned into a climate one.

## How it was established

**The file's units are not in doubt.** `exoplasim/makestellarspec.py` builds
these files from PHOENIX BT-Settl spectra, converts to W/m2/micron, resamples
onto two log-spaced blocks (0.2-0.75 um and 0.75-100 um in the hi-res file), and
divides by the speed of light. Column two is a spectral flux density. `radmod.f90`
uses it by assigning it straight into `bb3`, the same array its blackbody branch
fills with a Planck function, so a density is what the model expects and what it
gets.

**The model's own log confirms the reading.** `MOST_DIAG.00000` in the 0.96
baseline run reports "Energy fraction below 0.75 microns: 0.11588". Integrating
the file independently gives 0.1155. For reference, `radmod.f90:198` states that
its default partitioning of 0.517 is what a 5772 K solar spectrum produces. A
4965 K blackbody, this world's declared effective temperature, gives 0.4125.

So the star this project has been simulating puts 11.6% of its shortwave below
0.75 um where the declared star would put about 41%.

**Fitting a blackbody to each shipped spectrum** places all four in the M dwarf
range, and identifies the set:

| file | best-fit T | peak | almost certainly |
| --- | --- | --- | --- |
| `trap.dat` | 2545 K | 1.16 um | TRAPPIST-1, M8V, 2566 K |
| `wolf.dat` | 3000 K | 0.82 um | Wolf 1061, M3V |
| `gj667.dat` | 3015 K | 0.75 um | GJ 667 C, M1.5V |
| `k2.dat` | 3085 K | 0.75 um | K2-18, M2.5V |

ExoPlaSim ships no K dwarf spectrum at all.

**Provenance agrees.** The four files entered ExoPlaSim in a single commit,
415d119, "Added stellar spectra and haze data for study", by Maureen Cohen in
April 2023, alongside haze constants for Wolf 1061. That is a study of tidally
locked M dwarf planets, and TRAPPIST-1 e, Wolf 1061 c, GJ 667 Cc and K2-18 b are
its targets.

**Two independent files of the same star agree.** ExoPlaSim separately ships
`plasim/src/specs/k2-18.dat`, produced by a different author through a different
path. It fits 3055 K against `k2.dat`'s 3085 K, and the two agree to 6.1% in
shape across 0.4 to 5 um. They are not byte-identical, so this is not a copy; it
is two renderings of the same star.

The one thing not proven is the exact identification. `k2.dat` sits 4.4% from
`gj667.dat` and 4.7% from `wolf.dat` against 6.3% from `k2-18.dat`, so on shape
alone it could be any early-to-mid M dwarf. What is proven, and is what matters,
is that it is not a K2V.

## What it changes

Total insolation is unaffected: `radmod` normalises, so the flux is whatever
`baseline_flux_earth` says. Only the spectral weighting moves, and it moves the
things that are wavelength-selective.

Recomputing `radmod`'s own albedo integrals over its own reflectance data, once
with the k2 file and once with a 4965 K blackbody. The k2 column reproduces the
run log to better than 0.5% on every line, which is what makes the other column
trustworthy:

| surface | as run (k2) | declared 4965 K | change |
| --- | --- | --- | --- |
| snow, mean | 0.3956 | 0.5683 | **+0.173** |
| snow, maximum | 0.5674 | 0.7664 | +0.199 |
| snow, minimum | 0.2658 | 0.3793 | +0.114 |
| sea ice, maximum | 0.4606 | 0.6628 | +0.202 |
| sea ice, minimum | 0.3308 | 0.4738 | +0.143 |
| glacier, minimum | 0.3957 | 0.5683 | +0.173 |
| bare ground | 0.2238 | 0.2132 | -0.011 |
| ocean | 0.0637 | 0.0689 | +0.005 |

Every frozen surface in every completed run is 0.11 to 0.20 too dark. Snow and
ice are bright in the visible and dark in the near infrared, so shifting the
weighting into the infrared makes them absorb, and the more so the redder the
star. Ground and ocean are nearly spectrally flat over this range and barely
move, and the ground term moves the other way.

The land albedo this project computes from lithology is supplied directly as
codes 174, 175 and 176, so it is not affected. The error is confined to snow,
ice and glacier, which `radmod` derives from the spectrum itself.

**On the current baseline the effect is modest.** Integrating over the 0.96
climatology: the annual mean frozen-surface fraction is 5.1% of the planet, and
it sits where insolation averages only 70 W/m2, so the shortwave not reflected
comes to 0.41 to 0.73 W/m2 planetary mean, roughly 0.4 to 0.7 K at this
project's own 0.94 K per W/m2. That is a first-order estimate at fixed climate
and excludes the feedback it would set off.

So the 292.97 K baseline is probably good to well under a kelvin, which is
inside the 3 K design band. It is not, however, small against the 0.05 K per
orbit convergence criteria, and it is not small against the question the stellar
cycle experiment exists to ask.

**Where it is not modest is the cold branch.** The whole point of enabling
glaciers, of supplying forest fraction rather than leaving it uniform, and of
running a 0.91 to 1.01 stellar cycle, is albedo feedback. The static cold
endpoint sits near 282.2 K with far more ice than the baseline has, and the
error scales with ice extent. An ice-albedo feedback computed with snow 0.17 too
dark is a weaker feedback than the real one, biased toward not glaciating. Any
bistability or cycle result computed this way understates the cold branch.

**Rayleigh scattering moves too, and in the same direction.** `radmod` derives
its Rayleigh coefficient from the spectrum, weighted as lambda^-4, so a redder
star scatters far less. I could not reproduce the model's logged 0.21080
absolutely (my transcription gives 0.42802, a factor of 2.03 out, cause not
found), so only the ratio is quoted here, computed identically across spectra:
the declared 4965 K star gives **1.20x** the k2 file's Rayleigh coefficient, and
the Sun 1.40x. More scattering means a higher planetary albedo than the 0.152
currently reported, again a small cooling.

**For the biosphere**, this settles the question it came from. The 400-700 nm
fraction of shortwave is 0.085 under the k2 file against 0.313 under a 4965 K
blackbody, a factor of 3.7 on gross primary productivity. LPJ-GUESS should use
the latter until a real spectrum exists. See
`biosphere/notes/lpj-guess-porting-audit.md`.

## The fix, applied

A correct spectrum exists as `inputs/stellarspectra/k25v.dat` and
`k25v_hr.dat`, built by `scripts/build_stellar_spectrum.py`, and is now
selected (see "Confirmed in the model" below).

Source is the BT-Settl (CIFIST2011) grid, Allard and Homeier 2012, served by the
SVO Theoretical Spectra Server because phoenix.ens-lyon.fr did not respond. The
grid steps 100 K in this range, so the declared 4965 K is interpolated between
the 4900 K and 5000 K models at log g 4.5, [M/H] 0, rather than rounded to
either. The blend is log-linear in flux, which is right for a Planck-like
function over 100 K; a linear blend would bias the blue end low. Interpolating
rather than rounding is worth doing here: the two endpoints differ by 4.7% in
the band-1 fraction.

Conversion is done by ExoPlaSim's own `makestellarspec.convert` rather than a
reimplementation, so the output is format-compatible by construction. Verified:
965 and 2048 rows, wavelength grids identical to the shipped files, hi-res
blocks 0.200-0.749 and 0.750-100 um, all flux finite and positive. The absolute
flux scale differs from the shipped files because these are unnormalised
surface fluxes of different stars; that is harmless, since every quantity
`radmod` derives from the spectrum is a ratio.

`makestellarspec.py` needed two fixes first, recorded in
`patches/exoplasim-3.4.2-makestellarspec.patch`:

1. Lines 221 and 307, `np.loadtxt(Path(__file__).parent.resolve()+"/wvref.txt")`.
   Adding a `str` to a `Path` raises `TypeError`.
2. `np.trapz`, nine call sites, removed in NumPy 2.0. This environment is on
   2.5.2. Aliased to `np.trapezoid` rather than replaced outright so the tool
   still works on older NumPy.

The patch is applied to the vendored copy in `.venv`, which any reinstall
resets; `build_stellar_spectrum.py` checks for it and prints the reapply command
if it is missing. `run_exoplasim.py` now searches `inputs/stellarspectra/` before
the package, so a tracked spectrum wins over an untracked one.

### What the new spectrum gives

| quantity | k2.dat, as run | **k25v, new** | blackbody 4965 K | Sun 5772 K |
| --- | --- | --- | --- | --- |
| energy fraction < 0.75 um | 0.1155 | **0.3823** | 0.4125 | 0.5074 |
| 400-700 nm fraction | 0.0783 | **0.3093** | 0.3221 | 0.3903 |
| snow, mean | 0.3956 | **0.5428** | 0.5683 | 0.6122 |
| snow, maximum | 0.5674 | **0.7373** | 0.7664 | 0.8162 |
| snow, minimum | 0.2658 | **0.3625** | 0.3793 | 0.4081 |
| sea ice, maximum | 0.4606 | **0.6330** | 0.6628 | 0.7142 |
| sea ice, minimum | 0.3308 | **0.4527** | 0.4738 | 0.5102 |
| bare ground | 0.2238 | **0.2163** | 0.2132 | 0.2072 |
| ocean | 0.0637 | **0.0682** | 0.0689 | 0.0704 |

The new spectrum sits slightly redder than a 4965 K blackbody, which is the
expected sign: line blanketing removes blue flux that a blackbody keeps. So the
blackbody estimates used in the audit above were close but a little optimistic,
and the corrections are 5 to 15% smaller than they implied. The direction and
the magnitude both stand.

One further inconsistency worth knowing before regenerating anything:
`radmod.f90:224-227` zeroes all flux below 316 nm, but only on the star-file
path. The blackbody branch keeps its ultraviolet. So a star file and a blackbody
of the same temperature are not compared on equal terms, and the 4965 K column
above is therefore a slight overestimate of the blue end. It does not change any
conclusion here.

## What has to be redone

Nothing is invalidated in kind, and the direction of every error is known, but
these carry a spectrum-dependent bias that should be stated wherever they are
quoted:

- Every completed run under `exoplasim/runs/`, all three eras.
- The albedo bracket, though least of all: it turns on ground albedo, which
  moves by only -0.011, so the vegetated-versus-bare conclusion stands.
- The 292.97 K baseline and the `climatology_s096` products, by roughly 0.4 to
  0.7 K at fixed climate and more once feedback is allowed.
- The stellar cycle results most of all, because they are an albedo-feedback
  experiment run with the feedback too weak.

The carve verdict depends on climate through evaporation over catchments, so it
inherits the same bias, but at a magnitude far below the endorheic-fraction
spread it already reports.

## Confirmed in the model

The run check the audit deferred has now been done. Under `k25v`,
`MOST_DIAG.00000` reports an energy fraction below 0.75 microns of 0.38438,
against 0.11588 under `k2` and 0.3823 from integrating the new file
independently. The prediction was 0.382. The spectrum the model uses is
therefore the one intended, and the audit's albedo table can be read as applying
to the runs from here.

Verified independently before the switch: blackbody fits to the four shipped
spectra give 3117 K (`k2`), 3046 K (`gj667`), 3033 K (`wolf`) and 2572 K
(`trap`), none of them within 1800 K of this star. `k25v` fits 4702 K, below its
nominal 4965 K in the direction line blanketing predicts.

## One thing the correction exposed

`run_id` encoded resolution, flux, CO2, rotation, obliquity, eccentricity,
glaciers and a digest of every surface input, but not the spectrum. So the
corrected baseline at 0.96 would have been written into the completed `k2` run's
directory, where `finalize()` takes `sorted(glob("MOST*"))[-1]` and would have
copied out the wrong world's output under the new name. This is the same failure
`geography_tag` was written to prevent, reached by a different route: any
physical input not in the directory name is a collision waiting for the run that
changes it.

Both `run_id` and `cycle_run_id` now carry a spectrum marker. Directories written
before it carry none and are all `k2`; recomputing an id for one yields a name
that does not exist, so continuing a pre-fix run fails loudly instead of resuming
the wrong world.
