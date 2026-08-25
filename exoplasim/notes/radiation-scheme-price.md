# What the radiation scheme costs, and what a band-resolved one would

*Worldbuilding frame: this is a COMPUTE measurement of the Vesper project's
climate model on this desktop, and a comparison against a candidate radiation
code held in `references/`. Nothing here is about the simulated planet. Measured
2026-08-25 against `radmod.f90` as it now stands, at NLEV = 10, the layer count
every registry binary is built with.*

CLIM-84 asks the scheme to be priced before the routines in it are optimised.
CLIM-61 asks a band-resolved shortwave and longwave to be priced against the
broadband scheme. They are one question and this note answers the cost half of
it, plus what has become of the accuracy half since both rows were written.

## The instrument is retired instructions, not wall clock

Every number below is a retired-instruction count from `perf stat -e
instructions:u`, or a static count of call sites in a compiled object. That
choice is not a convenience. An instruction count is a property of the code and
the input; it does not move with what else the machine is doing, and it
reproduces here to about one part in 10^8 across repeats. A wall clock on a
shared desktop does not, and a contended timing number is worse than no number
because it looks like a measurement.

WHAT THAT INSTRUMENT CANNOT SAY. Instructions are not cycles. A vector call is
more instruction-efficient than eight scalar ones and need not be eight times
faster, because the limit may be elsewhere. So this note establishes the SHAPE
and the ORDER of the two schemes' cost and refuses to convert either into a cost
per timestep. That conversion is what still needs a quiet machine.

## The present scheme is O(NLEV^2) in scalar `pow` behind masks

Two facts about `lwr` set the whole cost.

**Its absorptances are evaluated per LEVEL PAIR, not per level.** The
transmissivity loop is `do jlev = 1, NLEV` around `do jlev2 = jlev, NLEV`, which
is NLEV(NLEV+1)/2 = 55 pairs at NLEV = 10. A broadband absorptance is a
nonlinear fit to a whole path, so optical depths cannot be added and every pair
of levels has the fit re-evaluated on its own accumulated amount. The inner body
carries 19 transcendental evaluations: 4 non-integer `**`, 6 `EXP` and 9
`LOG`/`LOG10`, counting the CH4 and N2O band CLIM-42 added. That is 1045 per
cell per longwave call.

`swr` is O(NLEV): 31 transcendental sites inside its level loops, of which 5 are
non-integer `**`, for 310 per cell per shortwave call.

**Every non-integer `**` in the radiation compiles to a SCALAR `pow`.** All 22
source sites sit inside a `where` block, and a masked array assignment is not
vectorised by GCC, so libmvec's eight-wide `pow` is never reached. Compiling
`radmod.f90` with the declared flag line and reading the object confirms it:
`nm -u` lists `pow`, `exp`, `log`, `log10` and the vector `_ZGVeN8v_exp`,
`_ZGVeN8v_log`, `_ZGVeN8v_sin`, `_ZGVeN8v_cos` -- and no `_ZGVeN8vv_pow` at all.
The relocation count is 87 scalar `pow` call sites, 146 scalar `exp`, 98 scalar
`log` and 25 `log10` against 169 vector `exp` and 44 vector `log`.

That is why radiation's transcendentals are 9.74% of T170 in
`where-the-time-goes.md`, `pow` alone 6.69%.

## What each transcendental costs, measured

Retired instructions per evaluation, glibc, `-march=znver4`, loop overhead
subtracted by a control that does the same arithmetic without the call. The
scalar figures use the model's own exponents and argument ranges; the vector
figures are the same expression written so GCC can vectorise it.

| form | instructions per element |
| --- | ---: |
| scalar `pow`, masked, as `radmod` writes it | 116 |
| scalar `pow`, unmasked C loop, model exponents | 111 to 112 |
| scalar `exp`, masked | 45.3 |
| scalar `log` | 53.0 |
| scalar `exp(a*log(x))` for `x**a` | 93 |
| vector `pow`, 8-wide `_ZGVeN8vv_pow` | 15.2 |
| vector `exp`, 8-wide `_ZGVeN8v_exp` | 5.1 |
| vector `exp(a*log(x))` | 12.3 |

Two things fall out of that table before any comparison.

**The exp/log identity CLIM-84 lists as a fallback is worth 17%, not a factor.**
93 against 112 for the same result, and it is less accurate. `merge()` in place
of `where` does not recover the vectorisation either: measured at 107
instructions per element, indistinguishable from the masked form, because both
branches are then evaluated and GCC still declines the vector call.

**Removing the mask is worth 7.6x on `pow` and 8.9x on `exp`.** That is a much
larger return than anything the row listed, and it is a separate decision rather
than a free win: the vector library's ULP bound is looser than the scalar one,
so it changes results, and the masks are not decoration -- they select between
two absorptance formulas by path amount. It is tracked separately.

## The two schemes side by side, transcendental work only

Per cell per radiation call, at NLEV = 10. The present scheme's `pow` is scalar
by measurement above; its `exp` and `log` are a mixture, so those are bracketed
between the vector and masked-scalar costs.

| | evaluations | instructions |
| --- | ---: | ---: |
| present `lwr`, 55 level pairs x 19 | 1045 | 30,200 to 66,700 |
| present `swr`, 10 levels x 31 | 310 | 7,200 to 18,200 |
| **present total** | **1355** | **37,400 to 84,900** |
| ga7 longwave, 88 to 166 solves x 10 layers x 1 to 2 `exp` | 880 to 3,320 | 4,500 to 16,900 |
| ga7 shortwave, 46 to 97 solves x 10 layers x 1 to 2 `exp` | 460 to 1,940 | 2,300 to 9,900 |
| **ga7 total** | **1,340 to 5,260** | **6,800 to 26,800** |

The solve counts are read out of the spectral files themselves rather than from
their summary headers. `references/socrates/data/spectra/ga7/sp_sw_ga7` block 5
holds 22 (band, gas) entries over 6 bands, 97 k-terms in total and 46 summed over
each band's largest gas; `sp_lw_ga7` holds 51 entries over 9 bands, 166 k-terms
in total and 88 summed over each band's largest. The low end assumes one gas per
band is treated by exponential-sum fitting and the rest by equivalent extinction,
which is how SOCRATES is normally configured; the high end assumes every k-term
gets its own solve. One or two `exp` per layer per solve is the layer
transmission plus, where the two-stream carries scattering, its eigenvalue.

**CLIM-61's premise does not survive this.** The row estimated "about an order of
magnitude of radiative work and ... roughly a twofold whole-model slowdown" from
the band counts alone: 6 shortwave bands with up to 12 k-terms and 9 longwave
bands against 2 and 1. The band counts are right and the conclusion does not
follow, because two things run the other way and together more than cancel it.
The broadband scheme is quadratic in levels where a k-term scheme is linear, and
the broadband scheme's kernel is a scalar `pow` behind a mask at 116
instructions where a k-term scheme's is a vectorisable `exp` at 5.1. On
transcendental work the candidate is CHEAPER at every corner of both brackets.

WHAT THIS COMPARISON EXCLUDES, and it is not small. The candidate's two-stream
and adding algebra, its k-table interpolation in pressure and temperature, its
per-band cloud overlap, its per-band Planck integral, and SOCRATES's own code
overhead are all uncounted; so is the present scheme's non-transcendental
algebra, which is also quadratic in levels. The candidate would also carry
roughly (bands x k-terms x levels) more state per column, which is a working-set
question against the 32 MB per-team target and is not answered here. So this is
evidence that the candidate is not an order of magnitude dearer, and it is not
evidence about which is faster.

## The accuracy half has moved, and it weakens the case for swapping

CLIM-61 was written on the claim that the broadband scheme "returns about 40
percent of a line-by-line answer for CO2 and the same for the trace-gas band".
**That claim was retracted.** `exoplasim/notes/trace-gas-band.md` traced it to
Earth greenhouse figures quoted from general knowledge, never sourced, and
probably not even the same quantity -- an attribution with band overlap shared
out, against a removal with water vapour left in place -- and reports no
identified defect anywhere in the chain.

What replaced it is a comparison that mostly confirms the present scheme.
`corrk-cross-check.md` puts the re-derived shortwave weights against HITRAN2020
correlated-k at the operating paths: `h2osww` agrees to 1.4%, `co2sww` to 0.03%,
the solar-weighted water vapour absorptance to 0.25%, and the Earth check passes
on both sides. Three things are left that a band-resolved scheme would fix
structurally:

- **The per-band attribution inside the CO2 total is wrong** and only survives by
  cancellation: at 2.7 um the derivation charges the band 0.000841 where
  correlated-k gives 0.000016, and at 2.0 um it undercharges by about the same
  amount with the opposite sign.
- **13% of the CO2 shortwave absorption falls outside every band Howard
  measured**, so the reconstruction cannot see it at all.
- **There is no water vapour continuum term**, which is why `h2o_sw_level` is
  carried as a bracket rather than a number.

- And the MAINTENANCE argument the row makes stands untouched: seven derived
  per-star corrections onto Lacis-Hansen and Sasamori, each separately derived
  and citable, exist because a fitted broadband absorptance cannot be re-weighted
  structurally and has to be re-derived per host. A band-resolved scheme performs
  that as configuration.

## The verdict on ordering, which is what CLIM-84 asked for

**Do not optimise the present scheme's transcendentals.** The swap is not refused
on cost -- the measurement above says the candidate's transcendental work is at
worst comparable and probably several times smaller -- so the optimisation would
be work thrown away, which is the exact outcome CLIM-84 exists to prevent. The
6.69% `pow` share remains the ceiling on what that work could return, and the two
routes inside it are now priced: 17% from the exp/log identity, and 7.6x on the
kernel from removing the `where` masks, which is not bit-identical and is its own
decision.

**The swap is not compelled by accuracy either.** The 40-percent figure that
motivated it is gone and the re-weighted broadband scheme reproduces
correlated-k to within a few percent on every quantity that has been compared.
What argues for it is structural: three known holes that re-weighting cannot
close, and seven manual per-star derivations that a band-resolved scheme would
retire.

**So the decision is a maintenance-and-capability one, and it needs one more
number this host cannot give today:** a cost per timestep for both schemes,
measured rather than counted. The command is a wrapped `radstep` on
`bench/bed_t170cold` under `perf`, with the candidate stood up far enough to run
one column, on a machine with nothing else on it.
