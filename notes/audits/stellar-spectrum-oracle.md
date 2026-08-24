# The k25v spectrum against an independent copy of the grid it came from

Measured on 2026-08-18.

`flux_fraction_band1` is the one number that carries this star's colour into the
surface energy balance, and until now the project computed it alone.
`check_consistency.py` verified that every stored copy agreed with every other,
which is internal consistency and not a test: nothing outside this repository
knew the answer. Something does now. The LMD Generic PCM ships a grid of 5,803
BT-Settl spectra degraded to R = 100, including the 4900 K and 5000 K points
`k25v` is interpolated between, and it was produced independently of anything
here.

Two results, and they point in different directions.

**The source data are validated.** The LMD grid and the SVO grid points pinned in
`exoplasim/scripts/build_stellar_spectrum.py` are the same BT-Settl release, and
integrating them under the same band definition agrees to 0.000325 in the band-1
fraction against a tolerance of 0.0010 argued below. The interpolation between the
two grid temperatures, which is load-bearing because the 100 K step is worth 4.8%
in the fraction, is treated the same way on both sides.

**The shipped file is not the source.** `exoplasim.makestellarspec.convert`
resamples onto the model's 2048-point grid with `np.interp`, a point sample rather
than a flux-conserving rebin. Point sampling steps over absorption-line cores,
line density is highest in the blue, and band 1 is therefore biased high by
0.00196. That is 6x the oracle's agreement with the source and 2x the tolerance,
and it is demonstrable without any external data at all.

## What was compared, and what would have meant wrong

### Provenance of both sides

| | project | oracle |
| --- | --- | --- |
| grid | BT-Settl (CIFIST2011), Allard and Homeier 2012 | Phoenix BT-Settl, per `README_stellar_files` |
| served by | SVO Theoretical Spectra Server, `newov2/ssap.php`, fids 3697 and 3850 | SVO Theoretical Spectra Server, `svo2.cab.inta-csic.es/theory/newov2` |
| grid points | Teff 4900 and 5000 K, log g 4.5, [M/H] 0, alpha 0 | same, `bt-settl_Teff{4900,5000}_logg4.5_met+0.0_R00100.dat` |
| resolution | native, R about 130,000 in the optical, 395,448 rows | R = 100, 571 rows, 0.1000 to 29.8873 um |
| normalisation | surface flux, erg/cm2/s/A | unit integral over its own 0.1 to 29.8873 um window |
| local copy | `/tmp/vesper-btsettl-cache`, sha256 matching `k25v_provenance.json` for both endpoints | `~/git/generic_pcm/LMDZ.GENERIC/datagcm/stellar_spectra/BT-Settl_stellar_spectra_grid/` |

Release identity was established rather than assumed, because BT-Settl has several
releases that differ and a version difference reported as an error would be worse
than not checking. The project's own full-resolution SVO endpoint was degraded to
the LMD wavelength grid by flux-conserving rebin and divided by the LMD file.
The quotient is a single constant to a median of 0.04%:

| band | Teff 4900 residual | Teff 5000 residual |
| --- | --- | --- |
| 0.15 to 0.30 um | +1.97% mean, 4.63% rms | +2.26% mean, 4.82% rms |
| 0.30 to 0.50 um | +0.11% mean, 3.66% rms | -0.00% mean, 3.59% rms |
| 0.50 to 0.75 um | +0.36% mean, 0.69% rms | +0.32% mean, 0.64% rms |
| 0.75 to 1.50 um | +0.17% mean, 0.34% rms | +0.17% mean, 0.31% rms |
| 1.50 to 5.00 um | +0.03% mean, 0.19% rms | +0.03% mean, 0.19% rms |
| 5.00 to 29.0 um | -0.01% mean, 0.11% rms | -0.01% mean, 0.11% rms |

The residual is confined to the near ultraviolet, where line blanketing is
deepest and any two degradation kernels differ most, and it vanishes into the
infrared continuum. Two unrelated model atmospheres do not agree to 0.01% over
five to twenty nine microns. Two further coincidences confirm it. The SVO
5000/4900 flux ratio over the LMD window is 1.08355 against 1.08417 for
(5000/4900)^4, which is 0.06%; and the LMD/SVO scale factor changes between the
two temperatures by exactly that ratio, which is what a per-file renormalisation
to unit integral produces and which the LMD integrals confirm at 1.000000 for
every Teff sampled.

The renormalisation is harmless here. A log-linear blend of two spectra each
carrying its own constant scale is the true blend times a constant, and a band
fraction divides the constant out. Measured: the LMD blend gives 0.382061 at its
native normalisation and 0.382061 rescaled to the SVO absolute scale.

### The tolerance, and what it has to be smaller than

Set before the oracle value was computed, from terms that do not depend on it.

| term | bound | how bounded | allowed |
| --- | --- | --- | --- |
| resolution, R = 100 against the model's 2048-point grid | 0.000011 | the project's own source integrated at full resolution and again after flux-conserving degradation onto the LMD wavelength grid; no LMD flux involved | 0.00005 |
| degradation kernel, boxcar against Gaussian of FWHM lambda/R | 0.000016 | same spectrum, both kernels | included above |
| truncation, LMD stops at 29.8873 um | 0.000017 | a 4965 K Planck curve puts 4.51e-05 of its flux beyond that, so the denominator moves by at most that and the fraction by b1 times it; measured effect 0.000014 | 0.00003 |
| blend treatment, log-linear against linear at fraction 0.65 | 0.000069 | both blends of the same two endpoints | 0.00010 |
| the model's discretisation against a clean split at 0.75 um | 0.000007 | `lib/stellar.py` against a direct integral of the same file | 0.00002 |
| two independent renderings of the same grid point, not explained by kernel choice | 0.00032 to 0.00038 | per endpoint, SVO degraded against the LMD file | 0.00060 |

Sum 0.0008, rounded up to **0.0010 absolute on the band-1 fraction**, or 0.26%
relative. A crude worst-case bound on the resolution term, mis-assigning half a
resolution element of flux at the 0.75 um edge, would have been 0.0036, but that
bound is 300x the measured term and using it would have made the check unable to
fail.

The tolerance has to be well below the errors it exists to catch, and it is:

| error the check must detect | cost in band-1 |
| --- | --- |
| rounding to the 5000 K grid point instead of interpolating | +0.0063 |
| rounding to the 4900 K grid point instead of interpolating | -0.0114 |
| log g 4.0 instead of 4.5 | +0.0041 |
| [M/H] -0.5 instead of 0.0 | +0.0114 |
| [M/H] +0.5 instead of 0.0 | -0.0084 |

Every one of those is 4 to 11 times the tolerance. **Falsification, stated in
advance: any disagreement above 0.0010 means one side is wrong about which star
it is describing, about how the 100 K step is interpolated, or about how the band
is defined.**

### The numbers

All under one definition, the model's: flux below 0.75 um over total flux, with
everything below `minwavel` = 316.036116751 nm discarded from band 1 as
`radmod.f90:500` does, and both sides truncated at 29.8873 um so the supports
match. `zcross/z1` is the star-dependent factor of the Rayleigh coefficient,
`solarini`'s lambda^-4-weighted flux integral over the band-1 integral.

| spectrum at 4965 K | band-1 | zcross/z1 |
| --- | --- | --- |
| SVO full resolution, blended (the project's own source) | 0.382386 | 15.646050 |
| SVO degraded to the LMD grid, flux-conserving | 0.382375 | 15.647981 |
| **LMD R = 100 as shipped (the oracle)** | **0.382061** | **15.643039** |
| SVO on the model's 2048-point grid, flux-conserving | 0.382383 | 15.645261 |
| SVO on the model's 2048-point grid, `np.interp` as built | 0.384340 | 15.485561 |
| `k25v_hr.dat` as shipped | 0.384404 | 15.485500 |

Truncation costs 0.000014 on both sides, so the untruncated figures are the same
to five decimals; `lib/stellar.py` reports 0.384383 for the shipped file under
`solarini`'s own discretisation, against 0.384390 for a direct integral of it.

**Oracle against the project's source: -0.000325, inside the tolerance at 0.33x.
The source is validated.**

**Oracle against the project's canonical value: -0.002322, outside the tolerance
at 2.3x.** The gap is not between the two grids. It is between the project's
source and the project's own file.

### The mechanism, demonstrated without the oracle

`makestellarspec.py` resampled with `f2 = np.interp(w2, w, f)`. The model grid
`w2` is 1024 log-spaced points from 0.2 to 0.75 um and 1024 from 0.75 to 100 um,
so R is about 775 in band 1 and about 209 in band 2, and the spectrum being
sampled is R about 130,000. A point sample of a line-blanketed spectrum lands in
the continuum more often than in a line core, so it overestimates; line density
is highest in the blue, so it overestimates band 1 more than band 2.

Replacing the point sample with a flux-conserving rebin onto **the same 2048-point
grid** gives 0.382383 against the point sample's 0.384340, and that replacement
is what `makestellarspec.py` now does: `_rebin_conserving` integrates the source
across each output bin and divides by the bin width, so the output's own
trapezoidal integral reproduces the source's, with the bin edges taken as the
ARITHMETIC midpoints of the output grid because that is the choice under which
`np.trapz` over the output equals the sum of the bin integrals exactly. The
defect and the argument are recorded at the top of the file. The flux-conserving
version reproduces the full-resolution integral to 0.000003, so the grid is not
the problem and the sampling is. Per-interval, point sample over flux-conserving:

| interval | ratio |
| --- | --- |
| 0.200 to 0.316 um | 1.03128 |
| 0.316 to 0.500 um | 0.97965 |
| 0.500 to 0.750 um | 1.01094 |
| 0.750 to 1.500 um | 0.99708 |
| 1.500 to 5.000 um | 0.98626 |
| 5.000 to 29.89 um | 0.98939 |

The shipped file is confirmed to be exactly that point sample of exactly these two
SVO endpoints: `k25v_hr.dat` times the speed of light reproduces `np.interp` of
the blended source to a median relative deviation of 7.1e-08.

The same mechanism is worth -1.03% in `zcross/z1`, and the lambda^-4 weight is why
it is larger there than in the band fraction. `lib/stellar.py`'s Rayleigh
coefficient of 0.712517 would be about 0.7198 on a flux-conserving resample.
`lib/stellar.py` itself is not in question: SPEC-2 verified it against the patched
Fortran to seven significant figures, so what moves here is the file it reads, not
the code that reads it.

Also biased, by the same mechanism and not separately measured here: the 965-point
`k25v.dat`, which is interpolated from the same point-sampled array and is what
`solarini` reads for the snow, sea ice, glacier and ground albedo integrals.

## What else the grid constrains

**The 100 K step.** `docs/src/pipeline/sequencing.md` loop C and `exoplasim/notes/stellar-spectrum-audit.md`
both state 4.7% as the band-1 difference between the two endpoints, which is the
argument for interpolating rather than rounding. Independently: 4.77% from the LMD
files, 4.75% from the SVO endpoints at full resolution.

**The surface gravity.** Derived here rather than taken: with `star.mass_solar`
0.80, `star.luminosity_solar` and `star.effective_temperature_k` from
`config/planet.yaml`, R/Rsun = sqrt(L) (Tsun/T)^2 = 0.76880 and log g = 4.5696
cgs. The BT-Settl grid steps 0.5 dex, so 4.5 is the nearest point and is what
`build_stellar_spectrum.py` pins. Interpolating to 4.57 instead would move band-1
by -0.00053, below the tolerance and opposite in sign to the resampling bias.

**Metallicity is the more consequential pin, and it is not in the config.**
`config/planet.yaml` is canonical for the star and declares no metallicity;
`METALLICITY = 0.0` lives only in `build_stellar_spectrum.py` and in
`k25v_provenance.json`. The grid says 0.5 dex is worth +0.0114 or -0.0084 in
band-1, which is 5x the resampling bias and 20x the log g rounding. Solar
metallicity is a reasonable choice for a K2.5V star, and the point is not that it
is wrong but that the star's most powerful undeclared parameter is declared
nowhere the config can see.

**PAR.** The 400 to 700 nm share of the same total, which
`exoplasim/notes/stellar-spectrum-audit.md` records as 0.3093 for `k25v` and which
sets LPJ-GUESS gross primary productivity: 0.311991 from the SVO source at full
resolution, 0.311460 from the oracle, 0.313002 from the shipped file. Same sign,
same size, same cause.

## What the discrepancy moves downstream

Small, and worth stating precisely so it is not re-litigated.

The correct band-1 fraction is **lower** than the stored one by 0.00232, so every
two-band frozen-surface albedo the model builds is slightly **too bright**. The
local slope of broadband snow albedo against band-1 fraction, read off the audit's
own table between the `k25v` and 4965 K blackbody columns, is 0.55 to 0.84, so
snow albedo is high by 0.0013 to 0.0020. Against the 0.10 to 0.17 the wrong star
was worth, that is two orders of magnitude smaller.

At fixed climate: annual mean frozen-surface fraction 5.1% of the planet, sitting
where insolation averages 70 W/m2, gives 0.005 to 0.007 W/m2 of shortwave, or
about 0.005 K at this project's 0.94 K per W/m2. Run-to-run spread on a converged
pair is 0.23 K, so it is 40x below the noise the project already accepts, and far
below the 0.5 to 2 K that dust is priced at.

The Rayleigh term is 1.0% on `rcoeff`. That is negligible beside the 3.4x that
SPEC-2 closed, and it does not reopen it.

For the biosphere the PAR fraction moves by 0.0015 on 0.3115, or 0.5% on gross
primary productivity, against the 3.7x the wrong star was worth.

**So nothing computed is invalidated and no run needs redoing for this.** What is
worth knowing is the sequencing. Per PHYS-2 no completed run integrated the
spectrum at all beyond its first orbit, and SPEC-1's baseline re-run on the
declared spectrum has not happened. Regenerating the spectrum before that run is
therefore nearly free, and after it is not.

Regenerating it in place is already safe, and it is worth saying so because the
opposite is the obvious guess. `run_id` carries the spectrum as a NAME, `k25v`,
not as a content hash, so the directory name alone would not separate a corrected
file from the one it replaced. It does not have to: CONS-3 stamps
`stellar_spectrum_digest`, both sha256 and hr_sha256, into the run manifest, and
`require_stellar_spectrum` raises rather than resuming across a change; a fresh
run into a nonempty directory is refused separately. `smoke_test.py` exercises
both directions of that guard.

## How this becomes a standing check

Not by depending on the external bundle. It is 1.2 GB, it is not in this
repository, and `check_consistency.py` cannot require it.

The check that can fail is implemented (SPEC-3): `build_stellar_spectrum.py`
integrates the blend at source resolution and writes the band-1 fraction and
`zcross/z1` into `k25v_provenance.json`, and `check_consistency.py` compares
`lib/stellar.py`'s value against that stored number instead of only against
other copies of itself. That is a comparison with a right answer: the
spectrum file is supposed to represent the source, and a resampler that does not
conserve flux makes it not. It would have fired on the day the file was built, it
needs nothing outside the repository at check time, and its tolerance is set by
the flux-conserving rebin rather than by the gap it is meant to detect.

The external oracle stays what it is: a one-off, repeatable from this document,
whose value was establishing that the source and the interpolation are right so
that the residual could be attributed to the resampler rather than argued about.
