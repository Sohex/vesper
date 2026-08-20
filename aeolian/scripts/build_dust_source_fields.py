#!/usr/bin/env python3
"""Write the three dust-source boundary fields the in-model emission reads. DUST-3.

    python aeolian/scripts/build_dust_source_fields.py
    python aeolian/scripts/build_dust_source_fields.py --z0 low   # shelter bracket
    python aeolian/scripts/build_dust_source_fields.py --self-test

`exoplasim/patches/exoplasim-3.4.2-dust-emission.patch` puts Kok et al. (2014)
equation 18 inside ExoPlaSim and reads its per-cell source map as three surface
boundary fields. This writes them, and writes beside them the thirteen namelist
constants the scheme needs and nothing else knows.

## Why the map stays outside the model

Which ground can emit is a question about lithology, the lake solution and the
soil. Answering it needs the mesh export, the solved lake extent and the pedology
clay field, none of which ExoPlaSim has ever seen, and `substrate_class` in
particular is the distinction that makes dust interesting here: it separates
loose closed-basin fill from consolidated bedrock, so the model does not emit as
much from forest as from salt pan. That is the standing decision for this
component. The model gets prepared fields and computes everything that responds
to the weather.

## Three fields, and the split is not arbitrary

Two of the offline chain's terms are LINEAR prefactors on the flux and can be
multiplied together out here; the third enters inside the nonlinearity and
cannot.

| code | name | quantity |
| ---: | --- | --- |
| 1801 | `dsrcw` | erodible fraction x clipped clay fraction |
| 1802 | `ddrage` | Marticorena-Bergametti (1995) drag partition of the bed |
| 1803 | `dwpr` | Fecan et al. (1999) residual moisture w', percent |

`ddrage` is separate because it scales `u*` itself, so it sits inside the
threshold comparison AND inside the fragmentation exponent; folding it into 1801
would be wrong rather than merely untidy. `dwpr` is separate because the moisture
gate also needs the model's own soil water, which is already in memory.

## What the self-test proves, and it can fail

`--self-test` reconstructs the offline flux from ONLY the three fields and the
thirteen scalars, using a Python transcription of the Fortran the patch adds, and
requires it to reproduce `build_dust.emission_over_weibull` to 1e-10. That is the
design claim -- three fields and no fewer, none of them dropped and none of them
in the wrong place -- stated as an identity that has a right answer.

It then re-runs the same comparison with four deliberate MUTATIONS: the drag
partition removed, the gravity scaling on the threshold removed, the Weibull
collapsed to its mean, and the Fecan residual zeroed. Every one of them must
BREAK the identity. A test that cannot be made to fail is not evidence, and this
one carries its own demonstration that it can be.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from math import gamma
from pathlib import Path

import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, DUST_CONFIG, PROJECT_ROOT  # noqa: E402
from builds import component_data, grid_export, resolution_of, soilmap  # noqa: E402
from paths import climatology_path, rel  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))
from sra import write_sra  # noqa: E402

import build_dust as bd  # noqa: E402

# The codes the patch registers in `surface_ini`. 1801-1803 were reserved for
# this by the prescribed-dust patch, which took 1811 and left them alone.
CODE_SRCW, CODE_DRAGE, CODE_WPR = 1801, 1802, 1803

BASELINE = ANALYSIS / "dust_baseline.json"

# The exponent cap the Fortran applies to (u*/u*t)**alpha. It is a numerical
# guard against -ffpe-trap=overflow and not physics, so the self-test and the
# real-field check both assert that it never binds.
EXP_CAP = 50.0


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# -- the namelist the model needs --------------------------------------------

def namelist_values(config: dict, cfg: dict, weibull_shape: float,
                    z0_aeolian: float) -> dict:
    """Every `aero_nl` constant the emission scheme needs, from `dust.yaml`.

    None of these has a Fortran default and `aero_ini` aborts without them, on
    purpose: `aeolian/config/dust.yaml` owns them and a plausible-looking number
    in the model would be a fourth place for each one to go stale.

    The two thresholds leave here ALREADY SCALED for this planet's gravity, by
    the fourth root of the gravity ratio. The scaling is applied on this side
    because this is the side that reads `config/planet.yaml`. Both take the same
    factor, so `(u*st - u*st0)/u*st0` is unchanged and Cd and alpha are
    invariant under it.
    """
    em, mc, sw, src = (cfg["emission"], cfg["moisture"],
                       cfg["subgrid_wind"], cfg["source"])
    gscale = bd.gravity_threshold_scaling(em, float(config["planet"]["gravity_m_s2"]))
    return {
        "LDUSTEMIT": 1,
        "DUSTZ0": float(z0_aeolian),
        "DUSTUST0": float(em["u_star_st0_m_s"]) * gscale,
        "DUSTUSTT": float(em["u_star_st_typical_m_s"]) * gscale,
        "DUSTRA0": float(em["rho_a0_kg_m3"]),
        "DUSTCD0": float(em["cd0"]),
        "DUSTCE": float(em["ce"]),
        "DUSTCA": float(em["c_alpha"]),
        "DUSTWK": float(weibull_shape),
        "DUSTNQ": int(sw["quadrature_points"]),
        "DUSTFA": float(mc["a"]),
        "DUSTFB": float(mc["b"]),
        "DUSTSND": float(src["snow_suppression_depth_m"]),
        "DUSTWCV": 1.0e5 / (float(mc["soil_depth_m"])
                            * float(mc["bulk_density_kg_m3"])),
    }


def measured_weibull_shape(report: dict) -> float:
    """The subgrid wind shape, from the offline run that MEASURED it.

    NOT `dust.yaml`'s declared 2.0. The declared value is a placeholder and the
    config says so; the measurement is `subgrid_wind.weibull_shape_used` in
    `aeolian/analysis/dust_baseline.json`, fitted to high-cadence winds under
    DUST-5. There is no silent fallback here for the reason `build_dust.py`
    removed its own: omitting the measurement once took the offline emission by
    a factor of 40 with nothing in the output saying so.

    Taking it from the same file the in-model total will be COMPARED against is
    also what makes that comparison mean anything: both arms then carry the same
    subgrid distribution and differ only in what the model resolves itself.
    """
    shape = report.get("subgrid_wind", {}).get("weibull_shape_used")
    if shape is None:
        raise SystemExit(
            f"{rel(BASELINE)} carries no measured Weibull shape. Re-run "
            "aeolian/scripts/build_dust.py with --gust-samples; a shape fitted "
            "to snapshots instead of high-cadence winds is wrong by a factor "
            "of 40 in emission and this refuses to guess one.")
    return float(shape)


# -- the Fortran, transcribed -------------------------------------------------

def emission_from_fields(srcw, drage, wpr, spd, temp, rho, wsoil, snow,
                         nl_values, gascon, gravity, sigma_bottom):
    """What `dustsrc` in the patch computes, in Python, from the same inputs.

    A TRANSCRIPTION and not a second implementation: every line here mirrors one
    there, so if the two ever disagree one of them has been edited alone. It
    exists so the identity in `self_test` can be checked without compiling the
    model, which is the only check available before a rebuild.

    Returns the emission in kg/m2/s per cell, already multiplied by `srcw`.
    """
    v = nl_values
    zsig = min(max(sigma_bottom, 0.5), 0.999)
    zref = np.maximum((gascon * temp / gravity) * np.log(1.0 / zsig), 2.0)
    ust = drage * bd.VON_KARMAN * spd / np.log(zref / v["DUSTZ0"])

    w = wsoil * v["DUSTWCV"]
    excess = np.maximum(w - wpr, 0.0)
    fm = np.where(excess > 0.0,
                  np.sqrt(1.0 + v["DUSTFA"] * np.where(excess > 0.0, excess, 1.0)
                          ** v["DUSTFB"]),
                  1.0)

    ustst = v["DUSTUSTT"] * fm
    ut = ustst * np.sqrt(v["DUSTRA0"] / rho)
    rel_ = (ustst - v["DUSTUST0"]) / v["DUSTUST0"]
    cd = v["DUSTCD0"] * np.exp(-v["DUSTCE"] * rel_)
    al = v["DUSTCA"] * rel_

    k, n = v["DUSTWK"], int(v["DUSTNQ"])
    q = (np.arange(n) + 0.5) / n
    zq = (-np.log(1.0 - q)) ** (1.0 / k)
    scale = ust / gamma(1.0 + 1.0 / k)

    total = np.zeros_like(srcw, dtype=float)
    capped = False
    for f in zq:
        u = scale * f
        active = (u > ut) & (ust > 0.0)
        if not np.any(active):
            continue
        ratio = np.where(active, u / np.maximum(ut, 1e-30), 1.0)
        arg = al * np.log(ratio)
        capped = capped or bool(np.any(active & (arg > EXP_CAP)))
        power = np.exp(np.minimum(arg, EXP_CAP))
        flux = cd * rho * (u * u - ut * ut) / ustst * power
        total += np.where(active, flux, 0.0)
    total /= n

    emit = np.where(snow > v["DUSTSND"], 0.0, total * srcw)
    return np.where(srcw > 0.0, emit, 0.0), capped


# -- the fields ---------------------------------------------------------------

def build_fields(config: dict, cfg: dict, lat, lon, z0_aeolian: float):
    """The three boundary fields, from the map `build_dust.py` already computes."""
    lakes = component_data("hydrography", config, strict=True) / "surface_water.nc"
    erodible, land_fraction, per_class, class_detail, terrain = bd.source_fractions(
        config, cfg, lakes)
    clay = bd.soil_clay_grid(soilmap(config), lat, lon)
    clay = np.nan_to_num(clay, nan=0.0)
    clay_pct = clay * 100.0

    f_clay = np.clip(clay, 0.0, float(cfg["emission"]["f_clay_max"]))
    srcw = erodible * f_clay
    drage = bd.drag_efficiency(np.full_like(srcw, z0_aeolian), cfg)
    wpr = bd.fecan_residual_moisture(clay_pct, cfg)
    return srcw, drage, wpr, land_fraction, per_class, class_detail, terrain, lakes


# -- the test that can fail ---------------------------------------------------

def _reference_flux(cfg, srcw, drage, spd, temp, rho, wsoil, clay_pct, snow,
                    gravity, gascon, sigma_bottom, values):
    """The same emission through `build_dust`'s own functions, not the mirror."""
    zsig = min(max(sigma_bottom, 0.5), 0.999)
    zref = np.maximum((gascon * temp / gravity) * np.log(1.0 / zsig), 2.0)
    ust = drage * bd.VON_KARMAN * spd / np.log(zref / values["DUSTZ0"])
    grav_pct = wsoil * values["DUSTWCV"]
    fm = bd.moisture_threshold_factor(grav_pct, clay_pct, cfg)
    u_st = values["DUSTUSTT"] * fm
    u_t = u_st * np.sqrt(values["DUSTRA0"] / rho)
    u_star_st = u_t * np.sqrt(rho / values["DUSTRA0"])
    flux = bd.emission_over_weibull(np.maximum(ust, 1e-6), u_t, rho,
                                    np.ones_like(srcw), u_star_st, cfg)
    return np.where(snow > values["DUSTSND"], 0.0, flux * srcw)


def self_test(config: dict, cfg: dict, tol: float = 1e-10) -> None:
    """Reconstruct the offline flux from the fields alone, then break it."""
    rng = np.random.default_rng(20260818)
    n = 4096
    clay_pct = rng.uniform(0.0, 45.0, n)
    clay = clay_pct / 100.0
    srcw = (rng.uniform(0.0, 0.6, n)
            * np.clip(clay, 0.0, float(cfg["emission"]["f_clay_max"])))
    spd = rng.uniform(0.1, 22.0, n)
    temp = rng.uniform(215.0, 320.0, n)
    rho = rng.uniform(0.7, 1.4, n)
    wsoil = rng.uniform(0.0, 0.30, n)
    snow = np.where(rng.uniform(0.0, 1.0, n) < 0.1, 0.5, 0.0)

    gravity = float(config["planet"]["gravity_m_s2"])
    gascon = 287.05
    sigma_bottom = 0.976

    shape = 2.012
    values = namelist_values(config, cfg, shape,
                             float(cfg["drag_partition"]["aeolian_z0_m"]))
    drage = np.full(n, bd.drag_efficiency(
        np.array([values["DUSTZ0"]]), cfg)[0])
    wpr = bd.fecan_residual_moisture(clay_pct, cfg)

    # `emission_over_weibull` reads the gravity factor off the config dict, the
    # way main() sets it. Same factor namelist_values applied to the thresholds.
    cfg = json.loads(json.dumps(cfg))
    cfg["emission"]["_gravity_scaling"] = bd.gravity_threshold_scaling(
        cfg["emission"], gravity)
    cfg["subgrid_wind"]["weibull_shape"] = shape

    # The REFERENCE side is always the correct one. A mutation changes only what
    # the model is given, which is the thing the fields and the namelist are
    # supposed to carry; mutating both sides together would test nothing, which
    # is exactly what the first draft of this did.
    ref = _reference_flux(cfg, srcw, drage, spd, temp, rho, wsoil,
                          clay_pct, snow, gravity, gascon, sigma_bottom, values)
    ref_scale = max(float(np.abs(ref).max()), 1e-300)

    def compare(label, vals=None, drage_=None, wpr_=None, srcw_=None):
        mirror, capped = emission_from_fields(
            srcw if srcw_ is None else srcw_,
            drage if drage_ is None else drage_,
            wpr if wpr_ is None else wpr_,
            spd, temp, rho, wsoil, snow,
            values if vals is None else vals, gascon, gravity, sigma_bottom)
        return float(np.abs(mirror - ref).max() / ref_scale), capped, label

    err, capped, _ = compare("identity")
    print(f"identity              max relative error {err:.3e}   (tolerance {tol:g})")
    if capped:
        raise SystemExit(
            "the e**50 exponent cap BOUND during the self-test. It is a guard "
            "against -ffpe-trap=overflow and not physics, so if it binds the "
            "Fortran and this transcription have stopped agreeing and the "
            "in-model flux is no longer Kok 18.")
    if not err < tol:
        raise SystemExit(
            f"the three boundary fields plus the namelist do NOT reproduce "
            f"build_dust.emission_over_weibull: max relative error {err:.3e}. "
            "Something is either missing from the fields or applied in the "
            "wrong place; the design claim that three fields suffice is what "
            "just failed.")

    # Now prove the test has teeth. Each mutation must break the identity.
    unscaled = dict(values)
    unscaled["DUSTUST0"] = float(cfg["emission"]["u_star_st0_m_s"])
    unscaled["DUSTUSTT"] = float(cfg["emission"]["u_star_st_typical_m_s"])
    delta = dict(values)
    delta["DUSTWK"] = 200.0

    broken = [
        compare("drag partition removed", drage_=np.ones_like(drage)),
        compare("gravity scaling removed", vals=unscaled),
        compare("Weibull collapsed to its mean", vals=delta),
        compare("Fecan residual zeroed", wpr_=np.zeros_like(wpr)),
    ]

    for err_b, _capped, label in broken:
        status = "BREAKS" if err_b > tol else "DOES NOT BREAK"
        print(f"  {label:<32} {err_b:.3e}  {status}")
    if any(e <= tol for e, _c, _l in broken):
        raise SystemExit(
            "a deliberate mutation did NOT break the identity, so the check "
            "above is not evidence of anything. Fix the test before trusting "
            "the fields.")
    print("self-test passed: the identity holds and every mutation breaks it")


# -- main ---------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--dust-config", type=Path, default=DUST_CONFIG)
    ap.add_argument("--baseline-report", type=Path, default=BASELINE)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="supplies the grid; defaults to baseline_climatology")
    ap.add_argument("--z0", default="central", choices=("low", "central", "high"),
                    help="which end of the aeolian roughness bracket to write")
    ap.add_argument("--self-test", action="store_true",
                    help="check the identity and the mutations; writes nothing")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cfg = yaml.safe_load(args.dust_config.read_text(encoding="utf-8"))

    if args.self_test:
        self_test(config, cfg)
        return

    from netCDF4 import Dataset

    if args.climatology is None:
        args.climatology = climatology_path()
    model = config["model"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    resolution = resolution_of(grid_export(config))

    with Dataset(args.climatology) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
    if (lat.size, lon.size) != (nlat, nlon):
        raise SystemExit(
            f"climatology grid is {lat.size}x{lon.size}, config says {nlat}x{nlon}")

    dp = cfg["drag_partition"]
    z0 = {"central": float(dp["aeolian_z0_m"]),
          "low": float(dp["aeolian_z0_bracket_m"][0]),
          "high": float(dp["aeolian_z0_bracket_m"][1])}[args.z0]

    report = json.loads(args.baseline_report.read_text(encoding="utf-8"))
    shape = measured_weibull_shape(report)
    values = namelist_values(config, cfg, shape, z0)

    (srcw, drage, wpr, land_fraction, per_class, class_detail,
     terrain, lakes) = build_fields(config, cfg, lat, lon, z0)

    for name, field in (("dsrcw", srcw), ("ddrage", drage), ("dwpr", wpr)):
        if not np.isfinite(field).all() or field.min() < 0.0:
            raise SystemExit(f"{name} has negatives or NaNs; refusing to write")

    out_dir = PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for code, field in ((CODE_SRCW, srcw), (CODE_DRAGE, drage), (CODE_WPR, wpr)):
        path = out_dir / f"orogen_{resolution}_surf_{code:04d}.sra"
        write_sra(path, code, field)
        outputs[code] = path

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
    provenance = {
        "note": "DUST-3 in-model dust emission. Surface codes 1801 dsrcw, 1802 "
                "ddrage and 1803 dwpr, plus the aero_nl group they need. "
                "Generated by aeolian/scripts/build_dust_source_fields.py; do "
                "not edit.",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "codes": {
            "1801": "dsrcw, erodible fraction x clipped clay fraction, dimensionless",
            "1802": "ddrage, MB95 drag partition of the erodible bed, dimensionless",
            "1803": "dwpr, Fecan residual soil moisture w', percent",
        },
        "z0_bracket_end": args.z0,
        "aeolian_z0_m": z0,
        "source_build": str(config["source_build"]),
        "terrain_hash": terrain,
        "lake_solution": rel(lakes),
        "soilmap": rel(soilmap(config)),
        "climatology": rel(args.climatology),
        "config_sha256": sha256_of(args.config),
        "dust_config_sha256": sha256_of(args.dust_config),
        "baseline_report": rel(args.baseline_report),
        "baseline_report_sha256": sha256_of(args.baseline_report),
        "weibull_shape_source":
            "subgrid_wind.weibull_shape_used in the baseline report, which is "
            "the MEASURED shape from DUST-5 and not dust.yaml's declared "
            "placeholder. The in-model total is compared against that report, "
            "so both arms have to carry the same subgrid distribution.",
        "erodible_weights": per_class,
        "namelist_values": values,
        "namelist_group": "aero_nl",
        "requires": "exoplasim/patches/exoplasim-3.4.2-dust-emission.patch, "
                    "applied and every binary rebuilt. A binary without it "
                    "cannot parse LDUSTEMIT and aborts in aero_ini, which is "
                    "the loud failure and the one to want.",
        "field_statistics": {
            "srcw_land_area_weighted_mean": float(np.average(srcw, weights=weights)),
            "srcw_max": float(srcw.max()),
            "srcw_cells_nonzero": int((srcw > 0.0).sum()),
            "drage_min": float(drage.min()),
            "drage_max": float(drage.max()),
            "wpr_max_percent": float(wpr.max()),
            "land_fraction_mean": float(np.average(land_fraction, weights=weights)),
        },
        "outputs": {str(code): rel(p) for code, p in outputs.items()},
        "output_sha256": {str(code): sha256_of(p) for code, p in outputs.items()},
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    report_path = outputs[CODE_SRCW].with_name(
        outputs[CODE_SRCW].stem + "_provenance.json")
    report_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    print(f"aeolian z0         {z0:g} m ({args.z0} end of the bracket)")
    print(f"Weibull shape      {shape:.5f}  (measured)")
    print(f"quadrature points  {values['DUSTNQ']}")
    print(f"thresholds         u*st0 {values['DUSTUST0']:.5f}  "
          f"u*st {values['DUSTUSTT']:.5f} m/s, gravity-scaled")
    print(f"drag partition     {drage.min():.5f} to {drage.max():.5f}")
    print(f"source cells       {int((srcw > 0.0).sum())} nonzero, max "
          f"dsrcw {srcw.max():.5f}")
    for code, path in outputs.items():
        print(f"wrote {rel(path)}")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
