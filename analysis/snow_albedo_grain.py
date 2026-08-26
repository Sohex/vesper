#!/usr/bin/env python3
"""A two-band snow albedo on a grain-radius axis, for this star, from the spectra the model already ships.

    python analysis/snow_albedo_grain.py

WORLDBUILDING CONTEXT, stated first because this file borrows vocabulary from a
real discipline: Vesper is an invented planet around a K2.5V host, and
everything below is about the climate model that simulates it -- the
reflectance of its modelled snow, the spectrum of its star, and the two
shortwave bands the model's radiation carries. Nothing here is a measurement of
Earth's cryosphere; the reflectance library is an Earth laboratory measurement
of a MATERIAL, and what this script does is re-weight it for a different star.

## What is here and what is missing

`specblock.f90` ships `fsnowalb`, `msnowalb` and `csnowalb`: fine, medium and
coarse granular snow reflectance over the same 965 wavelengths every other
endmember uses, from the JHU becknic library, at three stated effective grain
sizes. NO CODE PATH SELECTS AMONG THEM. `radmod.f90` weights the `iceblend`
family instead, and that family is not a grain axis: `combinedspec-snow_ice.ipynb`
built each blend as a mixture of five component spectra with the CLEAR ICE
fraction solved by bisection against a declared broadband albedo target under a
5772 K blackbody. So `landmod`'s `albsmin`-to-`albsmax` ramp interpolates
between two Earth-target mixtures on surface temperature, and the predictor is
temperature, not grain size.

The three grain-resolved spectra plus `lib/stellar.py`'s band integral are
between them enough to build the axis that is missing, for this star, from data
already in the tree. That is what this produces.

## Why the star matters here, and it is the row where it does

Snow reflectance is near unity through the visible and falls steeply through the
near infrared. A mid-K host puts far more of its output past the model's 0.75 um
band boundary than a G dwarf does, so the SAME snow is darker in broadband under
this star -- not because the material changed but because the flux moved into
the half of the spectrum where the material absorbs. The size of that shift is
this script's first output, and it is calculable rather than guessable.

## The functional form, and why it is fitted rather than imported

`references/climber-x/src/smb/smb_surface_par.f90:snow_albedo_dang` implements
Dang, Brandt and Warren (2015): a quadratic in `rn = log10(r/r0)` with `r0` at
100 um, one coefficient triple per band and per illumination. Its coefficients
are a fit to Earth's solar spectrum. Importing them onto this star would carry
the Sun's band weighting into a K dwarf's radiation, which is the whole of what
this row exists to avoid, so the FORM is taken and the COEFFICIENTS are fitted
here to the three star-weighted points. Three points and three coefficients: the
fit is an exact interpolation through them and is not a regression. What it buys
is a smooth, differentiable, extrapolable statement of the same three numbers in
the form a dust darkening term can scale against, which is what `dust-14` needs.

## The checks, all fixed before any of them was run

1. REPRODUCTION. `stellar.band_reflectances` is the integral used here, and it
   has to be the same arithmetic `radmod` runs. It is checked against the eight
   band pairs `radini` printed under "Finalized Albedos" for the blend family,
   which `analysis/ice_albedo.py` reproduces to 1e-9 by re-implementing
   `radmod`'s own loops. The two paths differ in a known way -- `radmod`
   integrates the numerator on the 965-point blend grid and this one interpolates
   the blend onto the 2048-point stellar grid -- so they agree to a grid
   resolution error rather than exactly. `REPRODUCTION_TOLERANCE` is one per cent
   of a reflectance that runs zero to one.
2. AGREEMENT OF THE THREE COPIES. The reflectances are read from
   `specblock.f90`, which is what the model compiles, and checked against
   `basespecs/*.npy`, which is what `surfacespecs.py` builds from, and against
   the JHU library files under `references/ecospeclib-all/`, which is where both
   came from. The grain radii are PARSED from the library headers rather than
   written down here, so a wrong spectrum cannot be paired with a right radius.
3. MONOTONICITY. Snow albedo falls with grain radius in both bands, because a
   larger grain gives a photon a longer path in ice per scattering event. A
   weighting that returned a rise would mean the spectra or the pairing were
   wrong, and this is the cheapest statement of that which can fail.

## What this does NOT establish

The grain radius is an axis, not a prediction. Nothing in this model evolves a
snow grain, so a run still has to be told which point on the axis its snow sits
at; how the radius EVOLVES is a mass balance's question and arrives with one.
The three points are also an Earth laboratory's fine, medium and coarse, which
is a bracket over the seasonal-snow range and not over firn or over a grain size
this world's own snow metamorphism would reach.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import stellar  # noqa: E402
from lib.paths import PROJECT_ROOT, rel  # noqa: E402

ROOT = PROJECT_ROOT
SPECBLOCK = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "specblock.f90"
BASESPECS = ROOT / "vendor" / "exoplasim" / "exoplasim" / "basespecs"
SURFACESPECS = ROOT / "vendor" / "exoplasim" / "exoplasim" / "surfacespecs.py"
ECOSPECLIB = ROOT / "references" / "ecospeclib-all"
OUTPUT = ROOT / "analysis" / "snow_albedo_grain.json"

# The three grain-resolved snow spectra, each named in all three places it
# appears: the array the model compiles, the file `surfacespecs.py` loads, and
# the library entry both were made from. The radius is NOT here -- it is parsed
# from the library header, so this table cannot pair a spectrum with a radius
# that belongs to another one.
GRAINS = [
    ("fsnowalb", "finesnow",
     "water.snow.finegranular.fine.all.fine_snw_.jhu.becknic.spectrum.txt"),
    ("msnowalb", "mediumsnow",
     "water.snow.mediumgranular.medium.all.medgran_snw_.jhu.becknic.spectrum.txt"),
    ("csnowalb", "coarsesnow",
     "water.snow.coarsegranular.coarse.all.coarse_snw_.jhu.becknic.spectrum.txt"),
]

R0_UM = 100.0
"""Dang's reference radius. The form is a quadratic in `log10(r/r0)`, so `r0`
sets only where the constant term is read off; it is carried at the published
value so a fitted triple and a published triple are comparable term by term."""

REPRODUCTION_TOLERANCE = 0.01
"""How far `stellar.band_reflectances` may sit from what `radini` printed.

FIXED BEFORE THE COMPARISON WAS RUN. This is not the 1e-9 `ice_albedo.py`
holds: that number checks a re-implementation of `radmod`'s own loops against
`radmod`, operation for operation, and either reproduces to floating point or
means nothing. This one compares two DIFFERENT quadratures of the same integral
over the same interval -- `radmod` integrates on the 965-point blend grid,
`band_reflectances` interpolates the blend onto the 2048-point stellar grid --
so the residual is a grid resolution error and has no reason to be zero. One
per cent of a reflectance that runs zero to one is the bar; anything larger
would mean the two are not the same integral.
"""

SPECTRUM_AGREEMENT = 1.0e-4
"""How far the three copies of one reflectance spectrum may differ, in percent
reflectance. The library ships four decimal places and both derived copies are
transcriptions of it, so the only permitted difference is the last printed
digit."""

SOLAR_TEMPERATURE_K = 5772.0
"""The blackbody `combinedspec-snow_ice.ipynb` solved its mixtures under, and
what `analysis/ice_albedo.py` and `lib/stellar.py` already use for the solar
arm. Carried so the star's effect is read off against a stated alternative
rather than against a remembered one."""

# Dang, Brandt and Warren (2015), as `snow_albedo_dang` implements them:
# constant, linear and quadratic coefficient in `rn = log10(r/r0)`, per band and
# per illumination. Held here to be COMPARED WITH, never used: the fitted
# coefficients below are what this star's arithmetic returns.
DANG_2015 = {
    "vis_diffuse": (0.9856, -0.0202, -0.0125),
    "nir_diffuse": (0.7493, -0.1820, -0.0388),
    "vis_direct": (0.9849, -0.0215, -0.0132),
    "nir_direct": (0.6596, -0.1927, -0.0229),
}
DANG_ZENITH_GRAIN_SCALING = {"vis": 0.781, "nir": 0.791}
"""`snow_albedo_dang`'s effective-radius correction for illumination angle,
`r_eff = r * (1 + k*(mu0 - 0.65)**2)`. Needed because the library spectra were
measured at a stated angle rather than hemispherically, so a comparison against
the published direct form has to be made at the geometry the library measured
in."""

DANG_BAND_EDGES_UM = {"vis": (0.2, 0.7), "nir": (0.7, 5.0)}
"""The band boundaries Dang's coefficients were fitted over. The model's are
0.34 to the row boundary near 0.75 and that boundary to 100 um, so the two are
NOT the same bands and the published triple is a reference point rather than an
answer. `band_edge_mismatch` below prices each end of the difference."""

# The eight blends `radini` weights, with what it printed for each. Copied from
# `analysis/ice_albedo.py`'s RECORDED block, which is checked against the model
# to 1e-9 there; here it is the fixed point check 1 measures against.
RECORDED = {
    "iceblend": (0.76382640121202905, 0.40611515565148087),
    "iceblendmax": (0.98278401164203333, 0.58529460781319775),
    "iceblendmin": (0.50782484289747332, 0.27261180703024718),
    "seaicemax": (0.89180978408407729, 0.47289790510791729),
    "seaicemin": (0.63575483791955212, 0.33949058411682098),
    "glacalbmin": (0.76378241785467849, 0.40619422519279852),
    "groundblend": (0.19354492153968911, 0.23002255831102814),
    "oceanblend": (0.075569967429263696, 0.063668098143514609),
}


def fortran_array(source: str, name: str, n: int = 965) -> np.ndarray:
    """One `real :: name(n) = (/ ... /)` initialiser out of the model's source.

    Read from the Fortran rather than from the `.npy` copy because the Fortran
    is what compiles: a spectrum edited in one and not the other would be a
    second declaration of one quantity, and check 2 exists to catch exactly
    that.
    """
    match = re.search(r"real\s*::\s*" + name + r"\(" + str(n) + r"\)\s*=\s*\(/(.*?)/\)",
                      source, re.S)
    if match is None:
        raise SystemExit(
            f"{rel(SPECBLOCK)} no longer declares {name}({n}) as an initialised "
            f"array. This script reads the spectra the model compiles, so a "
            f"rename or a reshape has to be settled here rather than worked "
            f"around.")
    body = re.sub(r"&\s*\n\s*&", "", match.group(1)).replace("&", "").replace("\n", "")
    values = np.array([float(t) for t in body.split(",") if t.strip()])
    if values.size != n:
        raise SystemExit(f"{name} parsed to {values.size} values, not {n}")
    return values


def library_entry(filename: str) -> tuple[float, np.ndarray, np.ndarray, str]:
    """A JHU becknic spectrum with its stated effective grain size.

    Returns the size in micrometres, the wavelength grid, the reflectance in
    percent, and the description line the size came from. The size is parsed
    rather than written down so it cannot drift away from the spectrum it
    belongs to.
    """
    path = ECOSPECLIB / filename
    text = path.read_text(errors="replace")
    header, _, body = text.partition("\n\n")
    described = [ln for ln in header.splitlines() if ln.startswith("Description:")]
    if not described:
        raise SystemExit(f"{rel(path)} carries no Description line to read a grain size from")
    size = re.search(r"([0-9.]+)\s*micrometers effective size", described[0])
    if size is None:
        raise SystemExit(
            f"{rel(path)} no longer states an effective grain size in its "
            f"Description. The grain axis is that number and there is no "
            f"second source for it in this tree.")
    rows = np.array([[float(v) for v in ln.split()]
                     for ln in body.splitlines() if ln.strip()])
    return float(size.group(1)), rows[:, 0], rows[:, 1], described[0].strip()


def dang_direct(band: str, radius_um: float, mu0: float) -> float:
    """Dang's published DIRECT albedo at a stated illumination angle.

    The library spectra are directional-hemispherical at 10 degrees from
    normal, so the published form has to be evaluated at that geometry for the
    comparison to be between two statements of one thing. `snow_albedo_dang`
    carries the angle in an effective radius rather than in the coefficients.
    """
    k = DANG_ZENITH_GRAIN_SCALING[band]
    r_eff = radius_um * (1.0 + k * (mu0 - 0.65) ** 2)
    rn = np.log10(r_eff / R0_UM)
    a0, a1, a2 = DANG_2015[f"{band}_direct"]
    return float(a0 + a1 * rn + a2 * rn * rn)


def fit_dang_form(radii_um, values) -> tuple[float, float, float]:
    """The quadratic in `log10(r/r0)` through three points. Exact, not a regression."""
    rn = np.log10(np.asarray(radii_um, dtype=float) / R0_UM)
    a2, a1, a0 = np.polyfit(rn, np.asarray(values, dtype=float), 2)
    return float(a0), float(a1), float(a2)


def _weighting(name: str | None, temperature_k: float | None):
    """One wavelength grid and one weighting curve, spanning the whole shortwave.

    `lib/stellar.py` hands its spectra out already split at the band boundary,
    because that is what the model's two-band scheme needs. Re-integrating a
    reflectance over SOMEBODY ELSE'S band edges needs them joined back up
    first, which is all this does.
    """
    if temperature_k is not None:
        wv1, wv2 = stellar._blackbody_grids()
        bb1 = stellar._planck(wv1, temperature_k)
        bb2 = stellar._planck(wv2, temperature_k)
    else:
        wavelength, flux = stellar.read_hires(stellar.spectrum_paths(name)[1])
        wv1, bb1, wv2, bb2 = stellar._bands(wavelength, flux)
    return np.concatenate([wv1, wv2]), np.concatenate([bb1, bb2])


def reflectance_over(lo_um: float, hi_um: float, wavelength_um, reflectance,
                     name: str | None = None,
                     temperature_k: float | None = None) -> float:
    """A reflectance flux-weighted over an ARBITRARY interval, not the model's bands.

    The comparison against Dang's published coefficients is only a comparison
    if both sides average over the same wavelengths, and they do not: his bands
    are 0.2-0.7 and 0.7-5.0 um and the model's are 0.34 to the row boundary and
    that boundary to 100. This is what lets the difference between the two be
    SPLIT into the part that is the band edges and the part that is the data,
    rather than attributed to one of them.
    """
    grid, curve = _weighting(name, temperature_k)
    um = grid * 1.0e-6
    inside = (grid >= lo_um * 1.0e-6) & (grid <= hi_um * 1.0e-6)
    r = np.interp(grid * 1.0e6, wavelength_um, reflectance)
    return float(np.trapezoid(r[inside] * curve[inside], grid[inside])
                 / np.trapezoid(curve[inside], grid[inside]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spectrum", default=None,
                    help="stellar spectrum name; the configured one by default")
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    source = SPECBLOCK.read_text()
    model_wavelength = fortran_array(source, "wavelengths")

    spec = importlib.util.spec_from_file_location("surfacespecs", SURFACESPECS)
    surfacespecs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(surfacespecs)
    if not np.allclose(surfacespecs.wvl, model_wavelength, atol=1.0e-9):
        raise SystemExit(
            "specblock.f90 and surfacespecs.py no longer carry the same "
            "wavelength grid, so the model's radiation and this project's "
            "analysis of it are indexing different wavelengths.")

    # ---- check 2: the three copies of each spectrum, and the grain radii ----
    grains, disagreements = [], []
    for fortran_name, npy_name, library_name in GRAINS:
        compiled = fortran_array(source, fortran_name)
        stored = np.load(BASESPECS / f"{npy_name}.npy")
        radius_um, library_wavelength, library_percent, description = library_entry(library_name)
        for label, other in (("basespecs", stored), ("ecospeclib", library_percent)):
            if other.shape != compiled.shape:
                disagreements.append(f"  {fortran_name} vs {label}: "
                                     f"{other.shape} against {compiled.shape}")
                continue
            worst = float(np.max(np.abs(other - compiled)))
            if worst > SPECTRUM_AGREEMENT:
                disagreements.append(f"  {fortran_name} vs {label}: worst "
                                     f"{worst:.6f} percent reflectance")
        if not np.allclose(library_wavelength, model_wavelength, atol=1.0e-9):
            disagreements.append(f"  {fortran_name}: the library entry is on a "
                                 f"different wavelength grid from the model's")
        grains.append((fortran_name, npy_name, library_name, radius_um,
                       np.clip(compiled / 100.0, 0.0, 1.0), description))
    if disagreements:
        raise SystemExit(
            "the model's compiled snow spectra, surfacespecs.py's copies and "
            "the library they came from no longer agree, so a grain radius "
            "read from the library no longer describes the spectrum the model "
            "runs:\n" + "\n".join(disagreements))

    # ---- check 1: is this the same integral radmod runs? ----
    failures = []
    for blend, printed in RECORDED.items():
        percent = getattr(surfacespecs, blend)
        got = stellar.band_reflectances(model_wavelength,
                                        np.clip(percent / 100.0, 0.0, 1.0),
                                        name=args.spectrum)
        for i, key in enumerate(("band1", "band2")):
            if abs(got[key] - printed[i]) > REPRODUCTION_TOLERANCE:
                failures.append(f"  {blend} band {i + 1}: {got[key]:.6f} "
                                f"against the printed {printed[i]:.6f}")
    if failures:
        raise SystemExit(
            "lib/stellar.py's band integral no longer agrees with what radmod "
            "printed for the blend family, so it is not the arithmetic the "
            "model runs and the grain axis below would be a second opinion "
            "rather than an extension:\n" + "\n".join(failures))
    reproduction_worst = max(
        abs(stellar.band_reflectances(model_wavelength,
                                      np.clip(getattr(surfacespecs, b) / 100.0, 0.0, 1.0),
                                      name=args.spectrum)[k] - p[i])
        for b, p in RECORDED.items() for i, k in enumerate(("band1", "band2")))

    # ---- the axis ----
    mu0 = float(np.cos(np.radians(10.0)))  # the library's measurement geometry
    rows, radii = [], []
    for fortran_name, npy_name, library_name, radius_um, reflectance, description in grains:
        star = stellar.band_reflectances(model_wavelength, reflectance, name=args.spectrum)
        sun = stellar.band_reflectances(model_wavelength, reflectance,
                                        temperature_k=SOLAR_TEMPERATURE_K)
        radii.append(radius_um)
        rows.append({
            "array": fortran_name,
            "basespec": npy_name,
            "library_entry": library_name,
            "library_description": description,
            "effective_grain_size_um": radius_um,
            "star": {"band1": round(star["band1"], 6),
                     "band2": round(star["band2"], 6),
                     "broadband": round(star["broadband"], 6)},
            "solar": {"band1": round(sun["band1"], 6),
                      "band2": round(sun["band2"], 6),
                      "broadband": round(sun["broadband"], 6)},
            "star_minus_solar_broadband": round(star["broadband"] - sun["broadband"], 6),
            "flux_outside_measured": round(star["flux_outside_measured"], 6),
            "dang_2015_direct_at_library_geometry": {
                "vis": round(dang_direct("vis", radius_um, mu0), 6),
                "nir": round(dang_direct("nir", radius_um, mu0), 6)},
        })

    # ---- check 3: monotonic in the radius, in both bands, on both arms ----
    order = np.argsort(radii)
    breaks = []
    for arm in ("star", "solar"):
        for band in ("band1", "band2", "broadband"):
            series = [rows[i][arm][band] for i in order]
            if any(b >= a for a, b in zip(series, series[1:])):
                breaks.append(f"  {arm} {band}: {series} over radii "
                              f"{[radii[i] for i in order]} um")
    if breaks:
        raise SystemExit(
            "the star-weighted snow albedo does not fall monotonically with "
            "grain radius, which it must: a larger grain is a longer path in "
            "ice per scattering event. The spectra and the radii are paired "
            "wrongly, or one of them is not what its name says:\n"
            + "\n".join(breaks))

    # ---- the fitted form, and what it is against the published one ----
    coefficients = {}
    for arm in ("star", "solar"):
        for band in ("band1", "band2", "broadband"):
            a0, a1, a2 = fit_dang_form(radii, [r[arm][band] for r in rows])
            coefficients[f"{arm}_{band}"] = [round(a0, 6), round(a1, 6), round(a2, 6)]

    # WHY THE SOLAR ARM AND THE PUBLISHED COEFFICIENTS DIFFER, SPLIT RATHER
    # THAN ASSERTED. Neither end of either band is Dang's, so a raw difference
    # between the two mixes a band-definition effect with a data one. Each
    # grain is re-integrated over Dang's OWN band edges under the same
    # blackbody, which separates them: what remains after that is a difference
    # between the JHU library's snow and Dang's two-stream model of it, and it
    # is not this row's to resolve.
    lo = float(model_wavelength[0])
    decomposition = []
    for row, (_, _, _, radius_um, reflectance, _) in zip(rows, grains):
        for band, model_band, model_key in (("vis", (lo, float(stellar.BAND_SPLIT_UM)), "band1"),
                                            ("nir", (float(stellar.BAND_SPLIT_UM), 100.0), "band2")):
            dlo, dhi = DANG_BAND_EDGES_UM[band]
            on_model = row["solar"][model_key]
            on_dang = reflectance_over(dlo, dhi, model_wavelength, reflectance,
                                       temperature_k=SOLAR_TEMPERATURE_K)
            published = row["dang_2015_direct_at_library_geometry"][band]
            decomposition.append({
                "effective_grain_size_um": radius_um,
                "band": band,
                "model_band_um": [round(model_band[0], 4), round(model_band[1], 4)],
                "dang_band_um": [dlo, dhi],
                "this_spectrum_on_model_band": round(on_model, 6),
                "this_spectrum_on_dang_band": round(on_dang, 6),
                "dang_published": round(published, 6),
                "band_definition_worth": round(on_dang - on_model, 6),
                "data_difference": round(published - on_dang, 6),
            })
    band_edge_mismatch = {
        "per_grain": decomposition,
        "note": "`band_definition_worth` is what re-integrating THE SAME "
                "spectrum over Dang's band edges is worth; `data_difference` "
                "is what is left, and it is a difference between the JHU "
                "library's snow and Dang's two-stream model of pure snow. In "
                "the near infrared the band definition is worth a few "
                "hundredths -- almost all of it the 0.70 to 0.75 um slice, "
                "which Dang counts as near infrared and the model counts as "
                "band 1, and where snow is still bright -- and the data "
                "difference is several times larger and near-constant in grain "
                "size. THAT IS THE REASON THIS ROW FITS THE FORM RATHER THAN "
                "IMPORTING THE COEFFICIENTS: the published constant term does "
                "not transfer onto these spectra even under the Sun, while the "
                "grain-size SLOPE does, which is the part a dust term scales "
                "against.",
    }

    lowres, hires = stellar.spectrum_paths(args.spectrum)
    band1_fraction, _ = stellar.band_fractions(name=args.spectrum)
    solar_band1_fraction, _ = stellar.blackbody_band_fractions(SOLAR_TEMPERATURE_K)

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/snow_albedo_grain.py",
        "finding": "the three grain-resolved snow spectra the model compiles "
                   "reach no code path, and the ramp that runs instead "
                   "interpolates two mixtures fitted to a solar broadband "
                   "target. Weighted against this star they give a two-band "
                   "snow albedo on a grain-radius axis, and the same snow is "
                   "darker in broadband here than under the Sun because the "
                   "flux moved past the band boundary rather than because the "
                   "material changed.",
        "spectrum_lowres": rel(lowres),
        "spectrum_hires": rel(hires),
        "spectrum_sha256": hashlib.sha256(hires.read_bytes()).hexdigest(),
        "specblock": rel(SPECBLOCK),
        "specblock_sha256": hashlib.sha256(SPECBLOCK.read_bytes()).hexdigest(),
        "band_split_um": float(stellar.BAND_SPLIT_UM),
        "band1_flux_fraction_star": round(band1_fraction, 6),
        "band1_flux_fraction_solar": round(solar_band1_fraction, 6),
        "solar_temperature_k": SOLAR_TEMPERATURE_K,
        "grain_reference_radius_um": R0_UM,
        "reproduction_tolerance": REPRODUCTION_TOLERANCE,
        "reproduction_worst": round(float(reproduction_worst), 8),
        "spectrum_agreement_percent": SPECTRUM_AGREEMENT,
        "grains": rows,
        "dang_2015_published": {k: list(v) for k, v in DANG_2015.items()},
        "fitted_coefficients": coefficients,
        "coefficient_order": "constant, linear and quadratic in log10(r/100um)",
        "band_edge_mismatch": band_edge_mismatch,
        "grain_size_reading": "the library states an 'effective size' in "
                              "micrometres and does not say radius or "
                              "diameter. It is read as a RADIUS, which is what "
                              "Dang's r0 = 100 um is and what makes 24, 82 and "
                              "178 um an ordinary seasonal-snow range; read as "
                              "diameters they would be finer than fresh snow. "
                              "A radius-for-diameter error would shift every "
                              "fitted linear coefficient by log10(2) times "
                              "itself and leave the three albedos untouched, "
                              "so it is a statement about the AXIS and not "
                              "about the values on it.",
        "what_this_does_not_establish": "the radius is an axis and not a "
                                        "prediction: nothing in this model "
                                        "evolves a snow grain, so a run is "
                                        "still told where on the axis its snow "
                                        "sits. The three points are a "
                                        "laboratory's fine, medium and coarse "
                                        "and bracket seasonal snow, not firn.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    print(f"grain axis, {rel(hires)} against a {SOLAR_TEMPERATURE_K:.0f} K blackbody\n")
    print(f"{'r (um)':>7}  {'star b1':>8}  {'star b2':>8}  {'star bb':>8}  "
          f"{'sun b1':>8}  {'sun b2':>8}  {'sun bb':>8}  {'d(bb)':>8}")
    for row in rows:
        s, u = row["star"], row["solar"]
        print(f"{row['effective_grain_size_um']:7.0f}  {s['band1']:8.4f}  "
              f"{s['band2']:8.4f}  {s['broadband']:8.4f}  {u['band1']:8.4f}  "
              f"{u['band2']:8.4f}  {u['broadband']:8.4f}  "
              f"{row['star_minus_solar_broadband']:+8.4f}")
    print(f"\nband 1 flux share: {band1_fraction:.4f} at this star, "
          f"{solar_band1_fraction:.4f} under the blackbody")
    print("\nfitted, constant/linear/quadratic in log10(r/100um):")
    for key in ("star_band1", "star_band2", "star_broadband",
                "solar_band1", "solar_band2", "solar_broadband"):
        a0, a1, a2 = coefficients[key]
        print(f"  {key:<16} {a0:+9.5f} {a1:+9.5f} {a2:+9.5f}")
    for label, key in (("vis_direct", "solar_band1"), ("nir_direct", "solar_band2")):
        a0, a1, a2 = DANG_2015[label]
        print(f"  {'Dang ' + label:<16} {a0:+9.5f} {a1:+9.5f} {a2:+9.5f}   "
              f"(different band edges; see band_edge_mismatch)")
    print(f"\nreproduced radmod's eight printed band pairs to "
          f"{reproduction_worst:.2e}, inside {REPRODUCTION_TOLERANCE}")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
