"""The star's shortwave band split, computed the way the model computes it.

ExoPlaSim's shortwave scheme is two bands divided at 0.75 um, and `zsolar1` --
the fraction of stellar flux below the split -- weights the two-band snow, sea
ice, glacier and ground albedos and divides the ozone transmissivity. It is the
one number that carries the host star's colour into the surface energy balance,
and getting it wrong is what previously made snow and ice 0.10 to 0.17 too dark
and cost this project a baseline re-run.

This module exists because three different values for it were in circulation at
once, from three different integrations of two different files, and none of them
was the one the model printed. `radmod.f90:solarini` is the authority on what
the split means, so this reproduces `solarini` rather than integrating the
spectrum in whatever way looks reasonable. The differences are not cosmetic:

- **The band edge belongs to band 1.** `solarini` integrates rows 1-1024 of the
  hi-res file as band 1 and rows 1025-2048 as band 2, then adds the interval
  BETWEEN them -- the gap from the last band-1 sample to the first band-2 one --
  to band 1. A naive `lam < 0.75` split drops that interval entirely.
- **Everything below `minwavel` is discarded.** `radmod.f90:235` zeroes the
  band-1 flux below 316.036116751 nm, so the ultraviolet the file carries down
  to 0.2 um is not in the model's normalisation. That is worth -0.0023 here.
- **The two halves of the file are separate grids**, log-spaced from 0.2 to
  0.75 um and from 0.75 to 100 um. The row index IS the band assignment, so a
  spectrum written at any resolution other than 2048 points would be silently
  mis-split. `assert_model_grid` refuses that case.

`k25v.dat`, the low-resolution companion, is NOT the file to use: it spans only
0.34 to 14.01 um, and `solarini` reads it only for the albedo integrals, taking
its energy fractions from the hi-res file. Integrating it gives a visibly
different share, because the ultraviolet below 0.34 um and everything past
14 um are simply absent from the denominator. That truncation artifact was once
carried into the dust optics.

## The two checks that can fail

`radmod.f90:207` states that its default partitioning of 0.517 is what a 5772 K
solar spectrum produces through this code. `solar_partition_identity` reproduces
0.517000 from `blackbody_band_fractions(5772.0)`, which is an identity with a
right answer rather than two formulations being compared. If the reproduction
drifts, this module is wrong about `solarini`, not merely different from it.

That checks the code. The FILE is checked separately, because a correct
integration of a wrong file is still wrong: the spectrum is supposed to
represent the BT-Settl blend it was resampled from, so the blend integrated at
its own resolution is a right answer for what this module should report.
`build_stellar_spectrum.py` writes that integral into `<name>_provenance.json`
and `check_consistency.py` compares it against `band1_fraction` and
`cross_section_ratio` here, within `SOURCE_BAND1_TOLERANCE` and
`SOURCE_CROSS_SECTION_TOLERANCE`. It caught a resampler that point-sampled an
R = 130,000 spectrum onto 2,048 points; `notes/audits/stellar-spectrum-oracle.md`
is the finding.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SPECTRA_DIR = PROJECT_ROOT / "exoplasim" / "inputs" / "stellarspectra"

BAND_SPLIT_UM = 0.75
"""ExoPlaSim's shortwave band edge. `radmod.f90:210`, `hinge = log10(7.5e-7)`."""

MIN_WAVELENGTH_NM = 316.036116751
"""`radmod.f90:86`, `minwavel`. Flux below this is removed, not just unresolved."""

MAX_WAVELENGTH_UM = 100.0
"""`radmod.f90:211`, `dl2 = (-4 - hinge)/1024`, i.e. the grid ends at 1e-4 m."""

HALF_ROWS = 1024
"""Rows per band in the hi-res file. `radmod.f90:229-233` splits by index."""

SOLAR_PARTITION = 0.517
"""What `radmod.f90:207` says a 5772 K spectrum gives through this code."""

SOURCE_BAND1_TOLERANCE = 5.0e-5
"""How far `band1_fraction` may sit from the same blend integrated at source resolution.

Set from the resampler's own accuracy, NOT from the size of any error it is
meant to catch. Two terms, both measured on the same spectrum so neither
depends on the star: the flux-conserving rebin onto the 2048-point grid
reproduces the full-resolution integral to about 3e-06, and `solarini`'s
discretisation -- the trapezoid that straddles `minwavel` and the band-edge
interval booked to band 1 -- differs from a clean split at 0.75 um by about
7e-06. Sum 1e-05, rounded up to 5e-05 for the arithmetic of a different blend.

That leaves the check able to fail at 40x below the +0.00196 the point-sampling
resampler was worth, and 200x below the 0.0114 that 0.5 dex of metallicity is
worth. A tolerance set from either of those would have passed the defect.
"""

SOURCE_CROSS_SECTION_TOLERANCE = 5.0e-4
"""The same, RELATIVE, for `cross_section_ratio`.

The lambda^-4 weight puts more of the integral in the blue, where the source is
most structured, so the rebin is less accurate here than on the band fraction:
about 5e-05 relative against the full-resolution integral. Rounded up an order
of magnitude, and still 20x below the 1.03% the point sample was worth.
"""

PLANCK_CONST_HC_OVER_K = 0.0143877735383
"""`radmod.f90`'s `const`, hc/k in metre kelvin."""


def spectrum_paths(name: str | None = None) -> tuple[Path, Path]:
    """The (low-resolution, hi-resolution) pair for a named spectrum.

    `name` defaults to `config/planet.yaml`'s `radiation.stellar_spectrum`.
    Project spectra shadow the ExoPlaSim package's, because the package lives in
    an untracked `.venv` that a reinstall resets while ours are tracked.
    """
    if name is None:
        import yaml
        config = yaml.safe_load(
            (PROJECT_ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
        name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        raise SystemExit("no stellar spectrum configured; "
                         "set radiation.stellar_spectrum in config/planet.yaml")
    stem = name[:-4] if name.endswith(".dat") else name
    roots = [SPECTRA_DIR]
    try:
        import exoplasim as _exo
        roots.append(Path(_exo.__file__).resolve().parent / "stellarspectra")
    except Exception:
        pass
    for root in roots:
        low, high = root / f"{stem}.dat", root / f"{stem}_hr.dat"
        if low.is_file() and high.is_file():
            return low, high
    raise SystemExit(
        f"stellar spectrum {stem} not found; looked for {stem}.dat and "
        f"{stem}_hr.dat under " + ", ".join(str(r) for r in roots))


def read_hires(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Wavelength in metres and spectral flux density, as `readdat` reads it.

    `readdat` skips one header line and then reads exactly `nitems` rows, so a
    file with more rows than the model asks for is truncated in silence.
    """
    data = np.loadtxt(path, skiprows=1)
    return data[:, 0] * 1.0e-6, data[:, 1]


def assert_model_grid(wavelength_m: np.ndarray) -> None:
    """Refuse a spectrum whose rows do not mean what `solarini` assumes.

    The model does not look at the wavelengths to decide which band a sample is
    in; it takes the first 1024 rows as band 1 and the rest as band 2. A file
    written at a different resolution, or with the split at a different row,
    would produce a partition that is wrong without producing an error.
    """
    if wavelength_m.size != 2 * HALF_ROWS:
        raise SystemExit(
            f"stellar spectrum has {wavelength_m.size} rows; radmod.f90 reads "
            f"exactly {2 * HALF_ROWS} and splits at row {HALF_ROWS}. Rebuild it "
            "with build_stellar_spectrum.py, which passes numwavelengths=2048.")
    split_um = BAND_SPLIT_UM * 1.0e-6
    if not (wavelength_m[HALF_ROWS - 1] < split_um <= wavelength_m[HALF_ROWS]):
        raise SystemExit(
            "stellar spectrum does not change band at row "
            f"{HALF_ROWS}: rows {HALF_ROWS} and {HALF_ROWS + 1} are "
            f"{wavelength_m[HALF_ROWS - 1] * 1e6:.6f} and "
            f"{wavelength_m[HALF_ROWS] * 1e6:.6f} um, and the model would put "
            "them on the wrong sides of the 0.75 um split.")


def _bands(wavelength_m: np.ndarray, flux: np.ndarray):
    """Split a hi-res file into `solarini`'s two bands, `minwavel` removed.

    The one place that knows the row split IS the band assignment and that the
    ultraviolet below `minwavel` is deleted rather than merely unresolved.
    """
    assert_model_grid(wavelength_m)
    wv1, bb1 = wavelength_m[:HALF_ROWS], flux[:HALF_ROWS].copy()
    wv2, bb2 = wavelength_m[HALF_ROWS:], flux[HALF_ROWS:]
    bb1[wv1 < MIN_WAVELENGTH_NM * 1.0e-9] = 0.0
    return wv1, bb1, wv2, bb2


def _partition(wv1, bb1, wv2, bb2) -> tuple[float, float]:
    """`solarini`'s trapezoidal integration, band edge included in band 1."""
    z1 = float(np.trapezoid(bb1, wv1))
    z2 = float(np.trapezoid(bb2, wv2))
    z1 += 0.5 * (bb1[-1] + bb2[0]) * (wv2[0] - wv1[-1])
    total = z1 + z2
    return z1 / total, z2 / total


def band_fractions(path: Path | None = None,
                   name: str | None = None) -> tuple[float, float]:
    """(below, above) 0.75 um for a spectrum file, as the model computes it.

    Pass either an explicit hi-res `path` or a spectrum `name`; with neither,
    the configured spectrum is used.
    """
    if path is None:
        path = spectrum_paths(name)[1]
    wavelength, flux = read_hires(Path(path))
    return _partition(*_bands(wavelength, flux))


def _blackbody_grids() -> tuple[np.ndarray, np.ndarray]:
    """`solarini`'s own logarithmic grids, `radmod.f90:210-217`."""
    lwmin = np.log10(MIN_WAVELENGTH_NM * 1.0e-9)
    hinge = np.log10(BAND_SPLIT_UM * 1.0e-6)
    lmax = np.log10(MAX_WAVELENGTH_UM * 1.0e-6)
    k = np.arange(HALF_ROWS)
    return (10.0 ** (lwmin + k * (hinge - lwmin) / HALF_ROWS),
            10.0 ** (hinge + k * (lmax - hinge) / HALF_ROWS))


def _planck(wavelength_m: np.ndarray, temperature_k: float) -> np.ndarray:
    """`solarini`'s Planck function, in its own arbitrary units."""
    c2 = PLANCK_CONST_HC_OVER_K / float(temperature_k)
    return 1.0 / (1.0e6 * wavelength_m) ** 5 / (np.exp(c2 / wavelength_m) - 1.0)


def blackbody_band_fractions(temperature_k: float) -> tuple[float, float]:
    """The same partition for `solarini`'s blackbody branch, on its own grid.

    This is what the model falls back to when no spectrum file reaches the
    namelist, and it is NOT the same answer: a photospheric model has line
    blanketing that a Planck curve does not, and the blanketing is in the blue.
    Kept here so the difference can be stated rather than guessed at.
    """
    wv1, wv2 = _blackbody_grids()
    return _partition(wv1, _planck(wv1, temperature_k),
                      wv2, _planck(wv2, temperature_k))


def solar_partition_identity(tolerance: float = 5.0e-4) -> float:
    """Reproduce `radmod.f90:207`'s stated 0.517 for the Sun, or raise.

    The one statement in the scheme with a right answer rather than a plausible
    one. Everything else here is checked against it.

    THROUGH BOTH ROUTES, at the same declared tolerance. The blackbody route is
    `_blackbody_grids` + `_planck` + `_partition`; a real spectrum takes
    `read_hires` + `_bands` + `_partition`, and the two shared only
    `_partition`, so `assert_model_grid`, the row split and the `minwavel`
    deletion -- everything that decides which flux lands in which band for the
    spectrum the model actually runs -- were outside the one statement with a
    right answer. The second arm below pushes the same Planck curve through
    `_bands`, which is that structure. The cut is NOT inert on this grid: the
    first row sits a float below `MIN_WAVELENGTH_NM` and is deleted, worth
    about 5e-5 in band 1, so the two arms are genuinely different arithmetic
    and both have to land on `SOLAR_PARTITION`. `read_hires` is the only part
    of the real-spectrum route this still does not reach. world-60x0.
    """
    band1, band2 = blackbody_band_fractions(5772.0)
    if abs(band1 + band2 - 1.0) > 1.0e-12:
        raise SystemExit(f"band fractions do not sum to 1: {band1} + {band2}")
    if abs(band1 - SOLAR_PARTITION) > tolerance:
        raise SystemExit(
            f"5772 K through solarini's grid gives {band1:.6f}, and "
            f"radmod.f90:207 says {SOLAR_PARTITION}. This module no longer "
            "reproduces the model's integration.")
    grid = np.concatenate(_blackbody_grids())
    through_bands, _ = _partition(*_bands(grid, _planck(grid, 5772.0)))
    if abs(through_bands - SOLAR_PARTITION) > tolerance:
        raise SystemExit(
            f"the same Planck curve through `_bands`, the split a real "
            f"spectrum takes, gives {through_bands:.6f} against "
            f"{SOLAR_PARTITION}. The row split, the grid assertion or the "
            "minwavel deletion no longer agrees with the model's.")
    return band1


def band_reflectances(wavelength_um: np.ndarray, reflectance: np.ndarray,
                      name: str | None = None, path: Path | None = None,
                      temperature_k: float | None = None) -> dict:
    """A surface reflectance integrated against the star, in the model's bands.

    The two-band shortwave scheme carries a star's spectral shape by moving
    flux BETWEEN the bands, so a surface whose two band constants are equal
    receives none of it. This is the integral that produces constants which are
    not equal, and it is the same operation for every surface: snow, glacier
    ice, sea ice, ocean, rock.

    `wavelength_um` and `reflectance` are the surface's own grid, ascending,
    reflectance as a fraction. With `temperature_k` the weighting is a Planck
    on `solarini`'s own grids, which is how the model's shipped constants were
    made; otherwise it is the named or configured stellar spectrum, split at
    the row boundary `_bands` applies rather than at a nominal 0.75 um.

    Returns band1, band2 and broadband reflectance, the band-1 flux fraction
    the recombination uses, the residual of that recombination against the
    broadband value, and `flux_outside_measured`: the share of stellar flux
    falling where the surface spectrum does not reach and its endpoint value is
    held instead. That last number is reported rather than hidden because it is
    the one assumption this integral makes.
    """
    if temperature_k is not None:
        wv1, wv2 = _blackbody_grids()
        bb1, bb2 = _planck(wv1, temperature_k), _planck(wv2, temperature_k)
    else:
        wavelength, flux = read_hires(path if path is not None
                                      else spectrum_paths(name)[1])
        wv1, bb1, wv2, bb2 = _bands(wavelength, flux)

    lo, hi = float(wavelength_um[0]), float(wavelength_um[-1])
    outside = 0.0
    total = 0.0
    means = []
    for wv, bb in ((wv1, bb1), (wv2, bb2)):
        um = wv * 1.0e6
        r = np.interp(um, wavelength_um, reflectance)  # endpoint-held outside
        means.append(float(np.trapezoid(r * bb, wv) / np.trapezoid(bb, wv)))
        total += float(np.trapezoid(bb, wv))
        beyond = (um < lo) | (um > hi)
        if beyond.any():
            outside += float(np.trapezoid(np.where(beyond, bb, 0.0), wv))

    band1, band2 = means
    f1, _ = _partition(wv1, bb1, wv2, bb2)
    # The same integral over both bands at once, which is what band1 and band2
    # must recombine to. Any residual is the band-edge interval `_partition`
    # books to band 1 while the two band means split it at the row boundary.
    broadband = f1 * band1 + (1.0 - f1) * band2
    um_all = np.concatenate([wv1, wv2]) * 1.0e6
    bb_all = np.concatenate([bb1, bb2])
    r_all = np.interp(um_all, wavelength_um, reflectance)
    direct = float(np.trapezoid(r_all * bb_all, np.concatenate([wv1, wv2]))
                   / np.trapezoid(bb_all, np.concatenate([wv1, wv2])))
    return {"band1": band1, "band2": band2, "broadband": broadband,
            "band1_flux_fraction": f1,
            "recombination_residual": broadband - direct,
            "flux_outside_measured": outside / total}


def _cross_section(wavelength_m: np.ndarray, flux: np.ndarray) -> float:
    """`solarini`'s `zcross`: the lambda^-4-weighted flux integral."""
    return float(np.trapezoid(flux / (wavelength_m * 1.0e6) ** 4, wavelength_m))


def _band_edge_cross_section(lo_w, lo_f, hi_w, hi_f) -> float:
    """The band-edge interval of `zcross`, which `solarini` books to band 1."""
    edge = 0.5 * (lo_f[-1] / (lo_w[-1] * 1e6) ** 4
                  + hi_f[0] / (hi_w[0] * 1e6) ** 4)
    return edge * (hi_w[0] - lo_w[-1])


def rayleigh_coefficient(name: str | None = None,
                         temperature_k: float | None = None,
                         as_the_model_does: bool = False) -> float:
    """`radmod.f90:311`'s `rcoeff`, the Rayleigh cross-section scaling.

    `rcoeff` multiplies the Rayleigh optical depth at `radmod.f90:1928`,
    normalised so a 5772 K spectrum gives exactly 1. A redder star scatters
    less, so it should fall below 1 and it does.

    ## The defect this reproduces

    `solarini` tabulates the 5772 K reference `bbg1`/`bbg2` on its OWN
    logarithmic grid at lines 224-226, and then, if a spectrum file was given,
    overwrites `wv1`/`wv2` with the file's wavelengths at lines 230-233. The
    reference integrals at lines 287-299 are evaluated afterwards, so they pair
    the reference's Planck VALUES with the file's WAVELENGTHS. The two grids do
    not span the same range -- the model's starts at `minwavel` and the file's
    at 0.2 um -- so the normalisation `zchi` is wrong whenever a spectrum file
    is used, and only then. The blackbody branch is self-consistent because
    nothing overwrites the grid.

    For `k25v` this is worth a factor of about 3.4, and `__main__` prints both
    numbers: `as_the_model_does=True` against the same star integrated on one
    grid. Passing the spectrum through to a run without fixing `solarini` would
    therefore weaken Rayleigh scattering 3.4x relative to the 0.862015 the runs
    have used, which is the 4965 K blackbody value and is self-consistent for
    the star it describes but not for this one. The FACTOR is a property of the
    defect and does not move when the spectrum file is rebuilt; the two values
    it was measured between do, so they are dated in SPEC-2
    beside the seven significant figures the patched Fortran agreed to.
    """
    if temperature_k is not None:
        gw1, gw2 = _blackbody_grids()
        wv1, wv2 = gw1, gw2
        bb1, bb2 = _planck(gw1, temperature_k), _planck(gw2, temperature_k)
    else:
        wavelength, flux = read_hires(spectrum_paths(name)[1])
        wv1, bb1, wv2, bb2 = _bands(wavelength, flux)
        gw1, gw2 = _blackbody_grids()
    ref1, ref2 = _planck(gw1, 5772.0), _planck(gw2, 5772.0)
    # The reference grid the model actually integrates the reference over.
    rw1, rw2 = (wv1, wv2) if as_the_model_does else (gw1, gw2)

    bridged = _band_edge_cross_section
    z1 = float(np.trapezoid(bb1, wv1)) + 0.5 * (bb1[-1] + bb2[0]) * (wv2[0] - wv1[-1])
    zcross = (_cross_section(wv1, bb1) + _cross_section(wv2, bb2)
              + bridged(wv1, bb1, wv2, bb2))
    zg = (float(np.trapezoid(ref1, rw1)) + float(np.trapezoid(ref2, rw2))
          + 0.5 * (ref1[-1] + ref2[0]) * (rw2[0] - rw1[-1]))
    zgcross = (_cross_section(rw1, ref1) + _cross_section(rw2, ref2)
               + bridged(rw1, ref1, rw2, ref2))
    return zcross * SOLAR_PARTITION / z1 / (zgcross / zg)


def cross_section_ratio(name: str | None = None,
                        path: Path | None = None) -> float:
    """`solarini`'s `zcross/z1` for a spectrum file, in um^-4.

    The star-dependent half of `rcoeff`: the lambda^-4-weighted flux integral
    over the band-1 flux integral, both over the whole file. `rcoeff` itself
    also carries the 5772 K reference normalisation, which the BT-Settl blend
    the file was built from knows nothing about, so this is the part that a
    source-resolution integral of that blend can be compared against.
    `check_consistency.py` does exactly that, against
    `SOURCE_CROSS_SECTION_TOLERANCE`.
    """
    if path is None:
        path = spectrum_paths(name)[1]
    wavelength, flux = read_hires(Path(path))
    wv1, bb1, wv2, bb2 = _bands(wavelength, flux)
    z1 = float(np.trapezoid(bb1, wv1)) + 0.5 * (bb1[-1] + bb2[0]) * (wv2[0] - wv1[-1])
    zcross = (_cross_section(wv1, bb1) + _cross_section(wv2, bb2)
              + _band_edge_cross_section(wv1, bb1, wv2, bb2))
    return zcross / z1


SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0
"""IAU 2015 Resolution B3 nominal solar effective temperature."""

OZONE_CHAPPUIS_BAND_UM = (0.44, 0.75)
"""Where Lacis and Hansen's visible ozone absorptance takes its flux from.

Their Eq. 8 gives ozone's visible absorption as a fraction of TOTAL INCIDENT
SOLAR flux, so the term carries the Sun's share of flux in the Chappuis band
inside it. On a non-solar host that share is wrong and the functional form,
which encodes ozone's cross-section shape, is not.
"""


def _band_share(wavelength_m: np.ndarray, flux: np.ndarray,
                lo_um: float, hi_um: float) -> float:
    """A band's share of a spectrum's whole integrated flux, edges interpolated.

    Interpolated rather than snapped to the nearest sample: the file's grid is
    logarithmic, so a snapped edge moves the band by a fraction of a percent and
    that is the size of the quantity this is used for.
    """
    um = np.asarray(wavelength_m, dtype=float) * 1.0e6
    inside = um[(um > lo_um) & (um < hi_um)]
    edges = np.concatenate([[lo_um], inside, [hi_um]])
    return float(np.trapezoid(np.interp(edges, um, flux), edges)
                 / np.trapezoid(flux, um))


def ozone_visible_weight(name: str | None = None) -> float:
    """`o3visw`: the Chappuis band's flux share here over the Sun's.

    DERIVED FROM THE SPECTRUM FILE THE MODEL READS, which is the whole reason it
    is a function and not a number: `build_stellar_spectrum.py` rewrites that
    file in place, and a resampling that moves the band-1 share moves this with
    it. The declared weight went stale exactly that way once, when the converter
    replaced a point sample with a flux-conserving rebin.

    The solar reference is a Planck curve at the nominal solar effective
    temperature, integrated on the SAME grid, so the comparison is a difference
    between two spectra rather than between two integrations.
    """
    wavelength, flux = read_hires(spectrum_paths(name)[1])
    lo, hi = OZONE_CHAPPUIS_BAND_UM
    star = _band_share(wavelength, flux, lo, hi)
    sun = _band_share(wavelength, _planck(wavelength, SOLAR_EFFECTIVE_TEMPERATURE_K),
                      lo, hi)
    return star / sun


def band1_fraction(name: str | None = None) -> float:
    """The canonical band-1 share for this world's star.

    The single value. Read it; do not integrate a spectrum yourself and do not
    copy the result into a script, a note or a namelist. `analysis/`,
    `world_state.json` and the dust optics all resolve through here, and the
    model computes the same number from the same file at run time.
    """
    solar_partition_identity()
    return band_fractions(name=name)[0]


if __name__ == "__main__":
    low, high = spectrum_paths()
    b1, b2 = band_fractions(path=high)
    bb1, bb2 = blackbody_band_fractions(4965.0)
    print(f"spectrum          : {high.name}")
    print(f"solar identity    : {solar_partition_identity():.6f} "
          f"(radmod.f90:207 says {SOLAR_PARTITION})")
    print(f"band 1 (< 0.75 um): {b1:.6f}")
    print(f"band 2 (> 0.75 um): {b2:.6f}")
    print(f"4965 K blackbody  : {bb1:.6f} / {bb2:.6f}  "
          "(what the model uses when no spectrum reaches the namelist)")
    print()
    print("Rayleigh coefficient, radmod.f90:311")
    print(f"  5772 K, the identity : "
          f"{rayleigh_coefficient(temperature_k=5772.0):.6f}  (must be 1)")
    print(f"  4965 K blackbody     : "
          f"{rayleigh_coefficient(temperature_k=4965.0):.6f}  "
          "(what every run after its first orbit used)")
    print(f"  this star, one grid  : {rayleigh_coefficient():.6f}")
    print(f"  zcross/z1, this star : {cross_section_ratio():.6f}  "
          "(the star-dependent half, what the provenance record checks)")
    print(f"  this star, as coded  : "
          f"{rayleigh_coefficient(as_the_model_does=True):.6f}  "
          "(the mismatched reference grid; 3.4x too weak)")
