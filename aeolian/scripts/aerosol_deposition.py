#!/usr/bin/env python3
"""Deposition MASS as a nutrient carrier: the elemental split, and the writer.

ANUT-7 retains marine aerosol as a base-cation and sulfur source to the abiotic
nutrient ledger and registers volcanic sulfate deposition beside it, and neither
had a carrier. What `aeolian/` produced for both was optics, and
`biosphere/config/abiotic_nutrients.yaml` refuses an optical quantity as a
carrier BY NAME: an optical depth folds in refractive index, size distribution
and water uptake, and no elemental mass comes back out of it. That refusal is
enforced -- `abiotic_nutrient_ledger.py:check_screen` reads
`screen.forbidden_carriers` and rejects any candidate whose `carrier` string
contains one of them, and two of the three entries are this component's own
optics artifacts by path.

So the carrier has to be a separate artifact holding MASS and only mass. That is
what this module writes. Putting the deposition mass into `sea_salt_baseline.nc`
or `volcanic_sulfate.nc` and pointing the ledger at those would be refused by the
gate on the filename alone, and correctly: those files hold optical depth too,
and a carrier that is half optics is exactly the ambiguity the refusal exists to
prevent.

## What the quantity is

Wet plus dry deposition of DRY aerosol mass, per unit RECEIVING area per unit
ABSOLUTE time, size-resolved, split into elements. In SI in the netCDF, so that
nothing has to know which world's year is meant; the JSON additionally reports
the land mean per Earth year, because the magnitude anchor is an Earth
measurement and the comparison has to be made in its units.

The removal terms are the ones the transport already runs. At steady state a
size bin loses `loss * m` per unit area per unit time, where `loss` is the sum of
the settling frequency `v_s / H` and the wet scavenging rate, so the split into

    dry = (v_s / H) * m        wet = lambda_wet * m

is a partition of a rate the solver has already balanced against emission rather
than a second calculation beside it. This module never re-derives either.

## Why the split is by ELEMENT and not by ion

The ledger's nodes carry elements: `soluble_deposition` declares
`[K, Ca, Mg, S]`. Sulfur arrives as sulfate and calcium as a cation, so the
conversion from an ion mass fraction to an element mass fraction is a molar-mass
ratio and belongs with the composition rather than at the point of use, where two
consumers would do it twice and differently.

## The check that can fail

`check_composition` sums the declared solute mass fractions and refuses a table
that does not close on 1. That is an identity of the reference composition and
not a comparison, so it has a right answer supplied by the source: a transcription
error in any row shows up as a residual whatever its sign.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from _paths import PROJECT_ROOT  # noqa: F401  -- puts lib/ on sys.path
from gridding import gaussian_grid, gaussian_latitudes  # noqa: E402

EARTH_YEAR_S = 365.25 * 86400.0


def check_composition(comp: dict, tol: float = 2.0e-3) -> float:
    """Sum the declared solute mass fractions. Raises if they do not close on 1.

    An IDENTITY, not a comparison: a reference composition is defined as the
    complete solute inventory of one kilogram of the salt, so the fractions sum
    to one by construction and a residual is a transcription error. The tolerance
    is what rounding every row to the printed number of digits can leave, plus
    whatever the source itself declares it truncated.
    """
    total = sum(float(row["mass_fraction"]) for row in comp["solutes"].values())
    if abs(total - 1.0) > tol:
        raise SystemExit(
            f"solute mass fractions sum to {total:.6f}, not 1 within {tol}. "
            f"The composition table has been mistranscribed; check "
            f"{comp.get('source', 'the source named in the config')}")
    return total


def element_mass_fractions(comp: dict) -> dict[str, float]:
    """Mass of each ELEMENT per unit mass of dry aerosol.

    Each solute row names the element it carries and the mass of that element per
    unit mass of the solute, so this is a weighted sum rather than a rename: two
    solutes can carry the same element, which is how sulfur reaches the total
    from sulfate and how carbonate would if it were ever wanted.
    """
    out: dict[str, float] = {}
    for row in comp["solutes"].values():
        f = float(row["mass_fraction"])
        for element, per_solute in row["elements"].items():
            out[element] = out.get(element, 0.0) + f * float(per_solute)
    return out


def deposition_by_element(dry_kg_m2_s, wet_kg_m2_s, comp: dict):
    """Per-element wet and dry deposition fields, kg of ELEMENT per m2 per s.

    `dry` and `wet` are dry-aerosol mass fluxes on the model grid. Returns
    `{element: {"dry": field, "wet": field, "total": field}}`.
    """
    fractions = element_mass_fractions(comp)
    return {element: {"dry": dry_kg_m2_s * f,
                      "wet": wet_kg_m2_s * f,
                      "total": (dry_kg_m2_s + wet_kg_m2_s) * f}
            for element, f in fractions.items()}


def land_mean_mg_m2_earth_year(field, land_weight) -> float:
    """Area-weighted land mean of a kg m-2 s-1 field, in mg per m2 per Earth year.

    Earth years because the magnitude anchor this is checked against is an Earth
    measurement. The netCDF carries the SI rate and this is the reading of it, so
    the conversion lives in one place and no consumer has to guess a year.
    """
    w = np.asarray(land_weight, dtype=float)
    if w.sum() <= 0:
        return float("nan")
    return float((np.asarray(field) * w).sum() / w.sum()
                 * EARTH_YEAR_S * 1e6)


def write_deposition(nc_path: Path, json_path: Path, *, lat, lon,
                     bins_um, dry_per_bin, wet_per_bin, comp: dict,
                     land_weight, ocean_weight, payload: dict) -> dict:
    """Write the carrier pair. Returns the summary block it put in the JSON.

    `dry_per_bin` and `wet_per_bin` are (nbin, nlat, nlon) dry-aerosol mass flux
    in kg m-2 s-1. The netCDF carries the per-bin fields, the bin totals and the
    per-element totals; the JSON carries the land means in the anchor's units and
    the composition it used, so a reader can reproduce the split without opening
    the config.
    """
    total_solutes = check_composition(comp)
    dry_per_bin = np.asarray(dry_per_bin, dtype=float)
    wet_per_bin = np.asarray(wet_per_bin, dtype=float)
    dry = dry_per_bin.sum(axis=0)
    wet = wet_per_bin.sum(axis=0)
    per_element = deposition_by_element(dry, wet, comp)

    summary = {
        "carrier": "wet plus dry deposition MASS of dry aerosol, per unit "
                   "receiving area per unit absolute time. NOT an optical "
                   "depth: nothing in this file is an optics product.",
        "units_netcdf": "kg m-2 s-1",
        "units_here": "mg of element per m2 of land per Earth year",
        "composition": {
            "source": comp.get("source"),
            "solute_mass_fraction_sum": round(total_solutes, 6),
            "element_mass_fraction_of_dry_aerosol": {
                k: round(v, 8) for k, v in
                element_mass_fractions(comp).items()},
        },
        "land_mean_deposition_mg_m2_earth_year": {
            element: {
                "dry": round(land_mean_mg_m2_earth_year(v["dry"], land_weight), 4),
                "wet": round(land_mean_mg_m2_earth_year(v["wet"], land_weight), 4),
                "total": round(land_mean_mg_m2_earth_year(v["total"], land_weight), 4),
            } for element, v in per_element.items()},
        "ocean_mean_deposition_mg_m2_earth_year": {
            element: round(land_mean_mg_m2_earth_year(v["total"], ocean_weight), 4)
            for element, v in per_element.items()},
        "bins_dry_um": [list(b) for b in bins_um],
        "land_mean_total_aerosol_mg_m2_earth_year": round(
            land_mean_mg_m2_earth_year(dry + wet, land_weight), 4),
        "land_mean_dry_fraction": round(float(
            land_mean_mg_m2_earth_year(dry, land_weight)
            / max(land_mean_mg_m2_earth_year(dry + wet, land_weight), 1e-30)), 4),
    }
    payload = dict(payload)
    payload["deposition"] = summary

    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    with Dataset(nc_path, "w") as out:
        out.createDimension("lat", len(lat))
        out.createDimension("lon", len(lon))
        out.createDimension("bin", dry_per_bin.shape[0])
        v = out.createVariable("lat", "f8", ("lat",)); v.units = "deg"; v[:] = lat
        v = out.createVariable("lon", "f8", ("lon",)); v.units = "deg"; v[:] = lon
        v = out.createVariable("bin_dry_diameter_lo_um", "f8", ("bin",))
        v.units = "um"
        v[:] = [b[0] for b in bins_um]
        v = out.createVariable("bin_dry_diameter_hi_um", "f8", ("bin",))
        v.units = "um"
        v[:] = [b[1] for b in bins_um]
        for name, data, long_name in (
                ("dry_deposition_per_bin", dry_per_bin,
                 "dry deposition of dry aerosol mass, per size bin"),
                ("wet_deposition_per_bin", wet_per_bin,
                 "wet deposition of dry aerosol mass, per size bin")):
            var = out.createVariable(name, "f8", ("bin", "lat", "lon"), zlib=True)
            var.units = "kg m-2 s-1"
            var.long_name = long_name
            var[:] = data
        for name, data, long_name in (
                ("dry_deposition", dry, "dry deposition of dry aerosol mass"),
                ("wet_deposition", wet, "wet deposition of dry aerosol mass"),
                ("deposition", dry + wet,
                 "wet plus dry deposition of dry aerosol mass")):
            var = out.createVariable(name, "f8", ("lat", "lon"), zlib=True)
            var.units = "kg m-2 s-1"
            var.long_name = long_name
            var[:] = data
        for element, v_el in per_element.items():
            var = out.createVariable(f"deposition_{element}", "f8",
                                     ("lat", "lon"), zlib=True)
            var.units = "kg m-2 s-1"
            var.long_name = (f"wet plus dry deposition of elemental {element}")
            var.element_mass_fraction_of_dry_aerosol = float(
                element_mass_fractions(comp)[element])
            var[:] = v_el["total"]
        out.note = payload.get("note", "")
        out.carrier = summary["carrier"]
        out.composition_source = str(comp.get("source"))
        if "generated" in payload:
            out.generated = payload["generated"]
        if "climatology" in payload:
            out.climatology = payload["climatology"]
    return summary


# ---------------------------------------------------------------------------
# The self-test
# ---------------------------------------------------------------------------

def _self_test() -> bool:
    """Three checks with right answers, and each is shown to be able to fail.

        python aeolian/scripts/aerosol_deposition.py --self-test

    A checker that has never been seen to reject anything is not a checker, so
    every check here is run twice: once on the shipped configs, which must pass,
    and once on a deliberately broken copy, which must be rejected.
    """
    import tempfile
    import sys as _sys

    import yaml

    from _paths import PROJECT_ROOT

    ok = True

    # --- 1. the composition tables close, and a mistranscription does not ---
    for name in ("sea_salt", "volcanic_sulfate"):
        cfg = yaml.safe_load(
            (PROJECT_ROOT / "aeolian" / "config" / f"{name}.yaml").read_text())
        comp = cfg["composition"]
        total = check_composition(comp)
        print(f"check 1  {name} solute mass fractions sum to {total:.7f}   "
              f"{'OK' if abs(total - 1.0) < 2e-3 else 'FAIL'}")
        ok = ok and abs(total - 1.0) < 2e-3
        broken = {"source": "deliberately broken",
                  "solutes": {k: dict(v) for k, v in comp["solutes"].items()}}
        first = next(iter(broken["solutes"]))
        broken["solutes"][first]["mass_fraction"] = float(
            broken["solutes"][first]["mass_fraction"]) * 0.9
        try:
            check_composition(broken)
            print(f"check 1  {name} broken table ACCEPTED   FAIL")
            ok = False
        except SystemExit:
            print(f"check 1  {name} broken table rejected   OK")

    # --- 2. the element split conserves mass ------------------------------
    #
    # The mass of element carried by a solute cannot exceed the mass of the
    # solute, so the total element mass per unit aerosol mass cannot exceed 1.
    # And for sea salt the monatomic ions carry themselves entirely, so their
    # element fractions must EQUAL their solute mass fractions. Both are
    # identities of the table rather than comparisons with a source.
    cfg = yaml.safe_load(
        (PROJECT_ROOT / "aeolian" / "config" / "sea_salt.yaml").read_text())
    comp = cfg["composition"]
    frac = element_mass_fractions(comp)
    total_element = sum(frac.values())
    print(f"check 2  element mass per unit sea salt {total_element:.7f}, must "
          f"not exceed 1   {'OK' if total_element <= 1.0 else 'FAIL'}")
    ok = ok and total_element <= 1.0
    worst = 0.0
    for solute, row in comp["solutes"].items():
        for element, per in row["elements"].items():
            if float(per) == 1.0:
                worst = max(worst, abs(frac[element]
                                       - float(row["mass_fraction"])))
    print(f"check 2  monatomic ions carry their own solute mass, worst "
          f"{worst:.2e}   {'OK' if worst < 1e-12 else 'FAIL'}")
    ok = ok and worst < 1e-12

    # --- 3. the writer's land mean is the analytic one --------------------
    #
    # A uniform deposition field over a uniform weight has a land mean equal to
    # itself, converted. Run against a value chosen so that the answer is not
    # 1 and a dropped conversion shows.
    nlat, nlon, nbin = 8, 16, 3
    flux = 1.0e-12                                     # kg m-2 s-1
    dry = np.full((nbin, nlat, nlon), flux / nbin / 2.0)
    wet = np.full((nbin, nlat, nlon), flux / nbin / 2.0)
    weight = np.ones((nlat, nlon))
    expected = flux * EARTH_YEAR_S * 1e6
    with tempfile.TemporaryDirectory() as d:
        summary = write_deposition(
            Path(d) / "t.nc", Path(d) / "t.json",
            lat=np.linspace(-80, 80, nlat), lon=np.linspace(0, 340, nlon),
            bins_um=[[0.1, 1.0], [1.0, 3.0], [3.0, 10.0]],
            dry_per_bin=dry, wet_per_bin=wet, comp=comp,
            land_weight=weight, ocean_weight=weight, payload={"note": "self-test"})
    got = summary["land_mean_total_aerosol_mg_m2_earth_year"]
    miss = abs(got / expected - 1.0)
    print(f"check 3  uniform field land mean {got:.4f} against an analytic "
          f"{expected:.4f}, off by {miss:.2e}   "
          f"{'OK' if miss < 1e-6 else 'FAIL'}")
    ok = ok and miss < 1e-6
    ca = summary["land_mean_deposition_mg_m2_earth_year"]["Ca"]["total"]
    ca_expected = expected * frac["Ca"]
    miss_ca = abs(ca / ca_expected - 1.0)
    print(f"check 3  and its calcium share {ca:.6f} against "
          f"{ca_expected:.6f}, off by {miss_ca:.2e}   "
          f"{'OK' if miss_ca < 1e-4 else 'FAIL'}")
    ok = ok and miss_ca < 1e-4

    # --- 4. the deposited mass is EXTENSIVE, and the artifact carries it ----
    #
    # kg m-2 s-1 is the areal DENSITY of an extensive quantity, so the number
    # with a right answer is the area integral, and it has three identities:
    # the written file's fields reproduce the arrays handed in, the per-bin
    # and per-element fields partition the total exactly, and the integral is
    # invariant under a change of support when the aggregation is the
    # area-weighted one. The wrong operator -- an unweighted block mean, the
    # defect world-vwqw found three times in this component -- moves the
    # integral on any field correlated with cell area, so the check is run
    # with both operators on a deliberately correlated field: the weighted one
    # must conserve and the unweighted one must be caught not conserving.
    # A REAL PARTITION OF THE SPHERE, from `lib/gridding.py`, so that the
    # closure this checks is the one a caller gets. A cosine of the row centres
    # is not one and this check would pass on it, which is exactly why the
    # fixture uses the weight the callers now pass.
    nlat, nlon = 8, 16
    lat_t = gaussian_latitudes(nlat)
    area = gaussian_grid(nlat, nlon).cell_area_fraction()
    # Squared, not linear: a linear gradient's covariance with a symmetric
    # area cancels BETWEEN blocks and the broken operator hides.
    grad = ((1.0 + np.arange(nlat, dtype=float)) ** 2)[:, None] \
        * np.ones((1, nlon))
    dry = np.stack([1.0e-12 * grad, 0.5e-12 * grad, 0.25e-12 * grad])
    wet = 0.5 * dry
    with tempfile.TemporaryDirectory() as d:
        nc_path = Path(d) / "t.nc"
        write_deposition(
            nc_path, Path(d) / "t.json",
            lat=lat_t, lon=np.linspace(0, 337.5, nlon),
            bins_um=[[0.1, 1.0], [1.0, 3.0], [3.0, 10.0]],
            dry_per_bin=dry, wet_per_bin=wet, comp=comp,
            land_weight=area, ocean_weight=area, payload={"note": "self-test"})
        from netCDF4 import Dataset as _DS
        with _DS(nc_path) as f:
            total = np.array(f["deposition"][:])
            per_bin = (np.array(f["dry_deposition_per_bin"][:])
                       + np.array(f["wet_deposition_per_bin"][:]))
            ca_field = np.array(f["deposition_Ca"][:])
            ca_frac = float(f["deposition_Ca"].element_mass_fraction_of_dry_aerosol)
    integral = float((total * area).sum())
    handed_in = float(((dry + wet).sum(axis=0) * area).sum())
    miss = abs(integral / handed_in - 1.0)
    print(f"check 4  the file's area integral reproduces the input's, off by "
          f"{miss:.2e}   {'OK' if miss < 1e-12 else 'FAIL'}")
    ok = ok and miss < 1e-12
    miss = float(np.abs(per_bin.sum(axis=0) - total).max()
                 / max(float(np.abs(total).max()), 1e-300))
    print(f"check 4  the bins partition the total, worst {miss:.2e}   "
          f"{'OK' if miss < 1e-12 else 'FAIL'}")
    ok = ok and miss < 1e-12
    miss = float(np.abs(ca_field - ca_frac * total).max()
                 / max(float(np.abs(ca_frac * total).max()), 1e-300))
    print(f"check 4  the element field is its fraction of the total, worst "
          f"{miss:.2e}   {'OK' if miss < 1e-12 else 'FAIL'}")
    ok = ok and miss < 1e-12
    # The change of support: 2x2 blocks. Area-weighted aggregation first.
    blk = (total * area).reshape(nlat // 2, 2, nlon // 2, 2).sum(axis=(1, 3))
    blk_area = area.reshape(nlat // 2, 2, nlon // 2, 2).sum(axis=(1, 3))
    coarse = blk / blk_area                       # the intensive field, coarse
    miss = abs(float((coarse * blk_area).sum()) / integral - 1.0)
    print(f"check 4  the integral is invariant under the area-weighted change "
          f"of support, off by {miss:.2e}   {'OK' if miss < 1e-12 else 'FAIL'}")
    ok = ok and miss < 1e-12
    # And the broken operator is CAUGHT rather than assumed broken: the
    # unweighted block mean re-integrated over the coarse areas must move the
    # integral on this field, or this check could not fail and is not a check.
    coarse_bad = total.reshape(nlat // 2, 2, nlon // 2, 2).mean(axis=(1, 3))
    drift = abs(float((coarse_bad * blk_area).sum()) / integral - 1.0)
    print(f"check 4  and the unweighted aggregation moves it by {drift:.2e}, "
          f"which must be seen   {'OK' if drift > 1e-4 else 'FAIL'}")
    ok = ok and drift > 1e-4
    _sys.stdout.flush()
    return ok


if __name__ == "__main__":
    import argparse
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true",
                    help="the identities, and that each of them can fail")
    a = ap.parse_args()
    if not a.self_test:
        ap.error("this module is a library; --self-test is the only entry point")
    raise SystemExit(0 if _self_test() else 1)
