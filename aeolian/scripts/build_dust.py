#!/usr/bin/env python3
"""Emit, transport and deposit mineral dust, offline, from a climatology.

    python aeolian/scripts/build_dust.py

DUST-1. `notes/dust.md` carries the decision and the reasoning; this is the
implementation of "what to build instead". Products are a per-cell deposition
flux, which loess (SURF-3), the phosphorus return leg and the minerals overlay
all want, and an aerosol optical depth map, which DUST-2 prices.

## Why this is not in the GCM

ExoPlaSim's own source term is `fcoeff * land_mask`: a fixed mixing ratio in the
bottom layer of every land gridbox. It emits as much from forest as from salt
pan, which is exactly the distinction that makes dust interesting on a world
where a quarter of the land is closed-basin fill. Enabling it would add a bias
rather than remove one.

Everything this needs is outside the model anyway. The lithology split between
crust and loose clastics, the solved lake extent, the derived ephemeral crust
and the pedology clay field are all artifacts of other components, and none of
them exists inside ExoPlaSim.

## What it computes, in order

1. **Friction velocity** from the climatology's lowest-level wind through a
   logarithmic profile on the roughness field `build_surface_roughness.py`
   writes. The reference height comes from the model's own lowest sigma level
   through the hypsometric relation, not from an assumed 10 m.

2. **A subgrid wind distribution**, Weibull, integrated over. This is not
   optional and the reason is arithmetic: emission goes as roughly u* cubed
   above a threshold, so a T42 gridbox mean emits almost nothing while the same
   mean with realistic variance emits a great deal. Cakmur et al. (2004).

3. **A threshold friction velocity**, gated on soil moisture through Fecan et
   al. (1999) with the pedology clay field, and partitioned for roughness
   through Marticorena and Bergametti (1995).

4. **Emission**, Kok et al. (2014) equation 18, over the erodible fraction of
   each cell. The erodible fraction is the part this project can do and the GCM
   cannot: lithology weights from `config/dust.yaml`, minus standing water from
   the lake solution, minus snow, minus vegetation.

5. **Transport**, as a column of dust advected on the mass-weighted wind of the
   lower atmosphere, solved to steady state per time bin.

6. **Removal**, gravitational settling with this world's gravity plus
   below-cloud scavenging on `prc` and `prl`. Without the wet term the fine
   mode's lifetime is out by an order of magnitude.

7. **Optical depth**, column mass times the mass extinction efficiency already
   computed in `analysis/dust_optics.json` against this star's spectrum.

## What it does not do

No dust-climate feedback: the climatology is an input and does not respond. That
is the whole reason this is defensible offline and it is also its limit. If the
answer crosses the reopening test in `notes/dust.md` -- land-mean optical depth
above 0.10 or 1.5 W/m2 of global forcing -- then the feedback is too large for a
prescribed field and the emission scheme belongs in the model.

No vertical structure: dust is a well-mixed column of a declared scale height.
No inter-bin microphysics, which is right for mineral dust because it neither
coagulates nor grows appreciably.

Every constant is declared in `aeolian/config/dust.yaml` with its source, and
anything that could not be taken from a source carries a bracket that this
script reports rather than collapsing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

from _paths import ANALYSIS, CONFIG, DUST_CONFIG, PROJECT_ROOT  # noqa: E402
from builds import component_data, grid_export, mesh_export, resolution_of, soilmap
from gridding import land_fraction_of_class, region_cells
from orogen import LAND, Export
from paths import climatology_path, rel
from provenance import require_build

VON_KARMAN = 0.4
R_DRY = 287.05          # J/kg/K
EARTH_YEAR_S = 365.25 * 86400.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:                                   # noqa: BLE001
        return None


# -- size distribution --------------------------------------------------------

def emitted_mass_fractions(cfg: dict) -> tuple[np.ndarray, np.ndarray]:
    """Mass fraction and mass-median diameter per bin, from Kok (2011) eq 6.

    The brittle-fragmentation volume distribution, integrated over each bin.
    Returns fractions summing to one over the binned range, and the
    volume-weighted mean diameter of each bin, which is what settles.
    """
    from math import erf
    sd = cfg["size_distribution"]
    ds, sg, lam = (sd["median_diameter_um"], sd["geometric_sigma"],
                   sd["crack_length_um"])
    edges = np.asarray(sd["bin_edges_um"], dtype=float)

    def dvdlnd(d):
        return (d * (1.0 + np.array([erf(x) for x in np.atleast_1d(
            np.log(d / ds) / (np.sqrt(2.0) * np.log(sg)))]))
            * np.exp(-((d / lam) ** 3)))

    frac, dbar = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        d = np.exp(np.linspace(np.log(lo), np.log(hi), 200))
        w = dvdlnd(d)
        m = np.trapezoid(w, np.log(d))
        frac.append(m)
        dbar.append(np.trapezoid(w * d, np.log(d)) / max(m, 1e-30))
    frac = np.array(frac)
    return frac / frac.sum(), np.array(dbar)


# -- physics ------------------------------------------------------------------

def settling_velocity(diameter_um: float, rho_p: float, rho_a: np.ndarray,
                      temperature_k: np.ndarray, gravity: float) -> np.ndarray:
    """Stokes settling with the Cunningham slip correction, at this world's g.

    Gravity is not Earth's and this is where that bites: at 12.81 m/s2 a given
    particle settles 1.31x faster than terrestrial intuition, so dust lifetime
    here is shorter and the transported fraction smaller.
    """
    d = diameter_um * 1e-6
    # Sutherland's law for dynamic viscosity, so the cold high-latitude columns
    # are not given warm-air viscosity.
    mu = 1.458e-6 * temperature_k ** 1.5 / (temperature_k + 110.4)
    # Mean free path, and the slip correction that matters below about a micron.
    lam = mu / rho_a * np.sqrt(np.pi / (2.0 * R_DRY * temperature_k))
    kn = 2.0 * lam / d
    cc = 1.0 + kn * (1.257 + 0.4 * np.exp(-1.1 / kn))
    return (rho_p - rho_a) * gravity * d ** 2 * cc / (18.0 * mu)


def drag_efficiency(z0_m: np.ndarray, cfg: dict) -> np.ndarray:
    """Fraction of the wind stress that reaches the erodible bed. MB95."""
    dp = cfg["drag_partition"]
    z0_cm = np.maximum(z0_m * 100.0, dp["z0s_cm"] * 1.0001)
    z0s = dp["z0s_cm"]
    denom = np.log(dp["a"] * (dp["x_cm"] / z0s) ** dp["b"])
    return np.clip(1.0 - np.log(z0_cm / z0s) / denom, 0.0, 1.0)


def flag_anomalous_bins(spd: np.ndarray, factor: float = 2.5) -> list[int]:
    """Time bins whose near-surface wind is wildly out of line with the rest.

    THIS EXISTS BECAUSE OF A REAL DEFECT, not as defensive habit. Every orbit of
    run_b014469b8091 carries a corrupted first output bin: `spd`, `ua` and `va`
    are inflated in the lower troposphere, worsening downward from 1.02x at the
    model top to 7.5x at the bottom level. It is systematic rather than a
    restart shock -- orbits 86, 87, 88 and 90 all show 7.50 to 7.58 -- so it is
    a property of how the first output interval of each year is accumulated.

    Consuming it silently cost this component a factor of 1700 in emission: bin
    0 alone was 100.00% of the annual total, because emission goes as roughly
    u* cubed above a threshold and that bin's wind is 7.5x too large.

    The same field is read by `carve_verdict.py`, `surface_water.py` and
    `export_carve_list.py`, all of which take `spd.mean(axis=0)[-1]` for the
    Penman wind and therefore inherit a 1.55x inflation of open-water
    evaporation. That is not this component's to fix and it is reported instead.
    """
    bottom = spd[:, -1, :, :].mean(axis=(1, 2))
    median = np.median(bottom)
    return [int(i) for i in np.where(bottom > factor * median)[0]]


def weibull_shape_from_snapshots(snapshot: Path, mask) -> float | None:
    """Weibull shape implied by the instantaneous winds, over erodible cells.

    The regular climatology is a 12-bin mean and has averaged away exactly the
    variance that drives dust emission, so the shape parameter cannot be read
    from it. The snapshot climatology keeps 32 instantaneous samples per orbit
    and can, through the coefficient of variation:

        cv = sqrt( Gamma(1 + 2/k) / Gamma(1 + 1/k)^2 - 1 )

    THIS IS AN UPPER BOUND ON k, and therefore a LOWER bound on emission. The
    samples are about 5.7 days apart, so they resolve synoptic variance and not
    the diurnal and sub-daily variance that produces real dust events. A
    distribution built from them is narrower than the true one.
    """
    from math import gamma as gfn
    if not snapshot.is_file():
        return None
    with Dataset(snapshot) as ds:
        if "spd" not in ds.variables:
            return None
        spd = np.asarray(ds["spd"][:, -1, :, :], dtype=float)
    mean, std = spd.mean(axis=0), spd.std(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        cv = np.where(mean > 0.1, std / np.maximum(mean, 1e-6), np.nan)
    ks = np.linspace(1.0, 12.0, 400)
    cvs = np.array([np.sqrt(gfn(1 + 2 / k) / gfn(1 + 1 / k) ** 2 - 1) for k in ks])
    sel = mask & np.isfinite(cv)
    if not np.any(sel):
        return None
    return float(np.median(np.interp(cv[sel], cvs[::-1], ks[::-1])))


def moisture_threshold_factor(gravimetric_pct: np.ndarray,
                              clay_pct: np.ndarray, cfg: dict) -> np.ndarray:
    """Fecan et al. (1999) equations 14 and 15."""
    mc = cfg["moisture"]
    w_res = mc["clay_quadratic"] * clay_pct ** 2 + mc["clay_linear"] * clay_pct
    excess = np.maximum(gravimetric_pct - w_res, 0.0)
    return np.sqrt(1.0 + mc["a"] * excess ** mc["b"])


def emission_over_weibull(u_star_mean: np.ndarray, u_star_t: np.ndarray,
                          rho_a: np.ndarray, f_clay: np.ndarray,
                          u_star_st: np.ndarray, cfg: dict) -> np.ndarray:
    """Kok (2014) eq 18, integrated over a Weibull distribution of gridbox wind.

    Emission is evaluated at each quadrature point of the subgrid distribution
    and mass-weighted, rather than evaluated once at the gridbox mean. The
    difference is not a refinement: at u* cubed above a threshold, the mean of
    the flux and the flux of the mean differ by orders of magnitude.
    """
    em, sw = cfg["emission"], cfg["subgrid_wind"]
    k = sw["weibull_shape"]
    n = int(sw["quadrature_points"])
    from math import gamma
    # Weibull scale from the mean: mean = c Gamma(1 + 1/k).
    scale = u_star_mean / gamma(1.0 + 1.0 / k)
    # Quadrature on the cumulative distribution, which puts points where the
    # mass is and keeps the steep tail resolved.
    q = (np.arange(n) + 0.5) / n
    # u at each quantile: u = c (-ln(1-q))^(1/k)
    factor = (-np.log(1.0 - q)) ** (1.0 / k)

    st0 = em["u_star_st0_m_s"]
    cd = em["cd0"] * np.exp(-em["ce"] * (u_star_st - st0) / st0)
    alpha = em["c_alpha"] * (u_star_st - st0) / st0

    total = np.zeros_like(u_star_mean)
    for f in factor:
        u = scale * f
        active = u > u_star_t
        if not np.any(active):
            continue
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            flux = (cd * f_clay * rho_a * (u ** 2 - u_star_t ** 2)
                    / np.maximum(u_star_st, 1e-9)
                    * (u / np.maximum(u_star_t, 1e-9)) ** alpha)
        total += np.where(active & np.isfinite(flux), np.maximum(flux, 0.0), 0.0)
    return total / n


def advect_to_steady_state(emission, u, v, loss_rate, lat, lon, cfg):
    """Column dust mass at steady state: emission in, advection and loss out.

    Upwind advection on the regular longitude / Gaussian latitude grid, stepped
    explicitly until the field stops moving. A steady state exists because the
    loss term is everywhere positive.
    """
    tr = cfg["transport"]
    nlat, nlon = emission.shape
    radius = 6.371e6 * 1.2                       # Vesper, 1.2 Earth radii
    dlon = np.deg2rad(360.0 / nlon)
    dphi = np.abs(np.gradient(np.deg2rad(lat)))
    coslat = np.maximum(np.cos(np.deg2rad(lat)),
                        tr.get("polar_coslat_floor", 1e-3))
    dx = (radius * coslat * dlon)[:, None] * np.ones((1, nlon))
    dy = (radius * dphi)[:, None] * np.ones((1, nlon))

    speed = np.maximum(np.abs(u) / dx + np.abs(v) / dy, 1e-12)
    kdiff = tr["diffusion_m2_s"]
    dt = 0.4 / np.max(speed + 2 * kdiff * (1 / dx ** 2 + 1 / dy ** 2) + loss_rate.max())
    m = np.zeros_like(emission)
    for step in range(int(tr["max_iterations"])):
        # Upwind fluxes, periodic in longitude and closed at the poles.
        mp_e, mp_w = np.roll(m, -1, axis=1), np.roll(m, 1, axis=1)
        fx = np.where(u > 0, u * m - u * mp_w, u * mp_e - u * m) / dx
        mp_n = np.vstack([m[:1], m[:-1]])
        mp_s = np.vstack([m[1:], m[-1:]])
        fy = np.where(v > 0, v * m - v * mp_n, v * mp_s - v * m) / dy
        lap = ((mp_e - 2 * m + mp_w) / dx ** 2 + (mp_n - 2 * m + mp_s) / dy ** 2)
        new = m + dt * (emission - loss_rate * m - fx - fy + kdiff * lap)
        new = np.maximum(new, 0.0)
        change = np.max(np.abs(new - m)) / max(np.max(new), 1e-30)
        m = new
        if change < tr["convergence_fraction"]:
            return m, step + 1, True
    return m, int(tr["max_iterations"]), False


# -- inputs -------------------------------------------------------------------

def read_sra_grid(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """One record of an ExoPlaSim .sra surface field."""
    values = []
    for line in path.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if len(parts) == 8 and all(p.lstrip("-").isdigit() for p in parts[:2]):
            continue                                    # header record
        values.extend(float(p) for p in parts)
    arr = np.asarray(values, dtype=float)
    return arr[:nlat * nlon].reshape(nlat, nlon)


def soil_clay_grid(path: Path, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Clay fraction per grid cell from the pedology soil map."""
    lines = path.read_text(encoding="utf-8").splitlines()
    head = lines[0].split()
    ilon, ilat, iclay = head.index("Lon"), head.index("Lat"), head.index("clay")
    out = np.full((lat.size, lon.size), np.nan)
    lon_signed = np.where(lon > 180.0, lon - 360.0, lon)
    for line in lines[1:]:
        p = line.split()
        j = int(np.argmin(np.abs(lat - float(p[ilat]))))
        i = int(np.argmin(np.abs(lon_signed - float(p[ilon]))))
        out[j, i] = float(p[iclay])
    return out


def source_fractions(config: dict, cfg: dict, lakes: Path):
    """Erodible land fraction per cell, and the pieces it is built from.

    This is the function the GCM cannot have. `substrate_class` is consolidated
    lithology plus closed-basin fill, so the only unconsolidated material the
    export carries is the fill; bedrock has no loose fines to lift. Standing
    water and the ephemeral salt crust come from the lake solution.
    """
    export = Export(mesh_export(config))
    grid_dir = grid_export(config)
    lit = json.loads((export.root / "manifest.json").read_text(
        encoding="utf-8"))["lithology"]
    ids = {r["code"]: int(r["id"]) for r in lit["rockClasses"]}
    substrate = export.field("substrate_class")

    weights = cfg["source"]["erodible_weight"]
    erodible = np.zeros(substrate.shape[0])
    per_class = {}
    for code, w in weights.items():
        mask = substrate == ids[code]
        per_class[code] = float(w)
        erodible[mask] = w
    # Standing water removes a region from the source entirely: a flooded playa
    # emits nothing, which is the seasonal half of why fill is not all source.
    with Dataset(lakes) as ds:
        wet = np.asarray(ds["lake"][:]).astype(bool)
    erodible[wet] = 0.0

    is_land = export.surface_class == LAND
    cell, nlat, nlon = region_cells(export, grid_dir)
    area = export.cell_area.astype(np.float64)
    ncell = nlat * nlon
    land_area = np.zeros(ncell)
    erod_area = np.zeros(ncell)
    np.add.at(land_area, cell[is_land], area[is_land])
    np.add.at(erod_area, cell[is_land], (area * erodible)[is_land])
    total_area = np.zeros(ncell)
    np.add.at(total_area, cell, area)
    with np.errstate(invalid="ignore", divide="ignore"):
        erodible_of_cell = np.where(total_area > 0, erod_area / total_area, 0.0)
    land_fraction = np.where(total_area > 0, land_area / total_area, 0.0)

    detail = {}
    for code in weights:
        detail[code] = land_fraction_of_class(
            export, grid_dir, substrate == ids[code]).tolist()
    return (erodible_of_cell.reshape(nlat, nlon),
            land_fraction.reshape(nlat, nlon), per_class, detail,
            export.terrain_hash)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--dust-config", type=Path, default=DUST_CONFIG)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="defaults to the configured baseline_climatology")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--variant", default="baseline",
                    choices=("baseline", "arid_bare_ground"),
                    help="`arid_bare_ground` lets cells drier than the declared "
                         "threshold emit regardless of rock class, which "
                         "brackets the assumption that all non-barren land "
                         "carries a canopy")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cfg = yaml.safe_load(args.dust_config.read_text(encoding="utf-8"))
    clim_path = args.climatology or climatology_path()
    output = args.output or (ANALYSIS / f"dust_{args.variant}.json")

    # Cross-component reads are checked, not assumed. lib/provenance.py says why.
    require_build(clim_path, "climatology", config)
    lakes = component_data("hydrography", config, strict=True) / "surface_water.nc"
    require_build(lakes, "lake solution", config)

    gravity = float(config["planet"]["gravity_m_s2"])
    resolution = resolution_of(grid_export(config))

    with Dataset(clim_path) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        lev = np.asarray(ds["lev"][:], dtype=float)
        spd = np.asarray(ds["spd"][:], dtype=float)
        ua = np.asarray(ds["ua"][:], dtype=float)
        va = np.asarray(ds["va"][:], dtype=float)
        mrso = np.asarray(ds["mrso"][:], dtype=float)
        prc = np.asarray(ds["prc"][:], dtype=float)
        prl = np.asarray(ds["prl"][:], dtype=float)
        pr = np.asarray(ds["pr"][:], dtype=float)
        tas = np.asarray(ds["tas"][:], dtype=float)
        ps = np.asarray(ds["ps"][:], dtype=float) * 100.0     # hPa to Pa
        snd = np.asarray(ds["snd"][:], dtype=float)
    nbin, nlat, nlon = tas.shape

    # A corrupted time bin is not a rounding error here: emission goes as u*
    # cubed above a threshold, so one bin with 7.5x the wind is the whole answer.
    bad_bins = flag_anomalous_bins(spd)
    good = [t for t in range(nbin) if t not in bad_bins]
    if bad_bins:
        print(f"  EXCLUDING time bins {bad_bins}: near-surface wind is out of "
              f"line with the other bins by more than 2.5x. See "
              f"flag_anomalous_bins for what this defect is and who else "
              f"inherits it.")
    if len(good) < 2:
        raise SystemExit("too few usable time bins; refusing to guess a climate")

    z0 = read_sra_grid(
        PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
        / f"orogen_{resolution}_surf_0173.sra", nlat, nlon)
    clay = soil_clay_grid(soilmap(config), lat, lon)
    erodible, land_fraction, per_class, class_detail, terrain = source_fractions(
        config, cfg, lakes)

    # The lowest model level's height, from the model's own sigma coordinate
    # through the hypsometric relation. Assuming 10 m would misstate u* by the
    # ratio of the logarithms, which on a rough surface is not small.
    sigma_bottom = float(lev[-1]) if lev[-1] > lev[0] else float(lev[0])
    if sigma_bottom > 1.5:                       # levels given in hPa, not sigma
        sigma_bottom = sigma_bottom * 100.0 / float(np.mean(ps))
    z_ref = (R_DRY * tas / gravity) * np.log(1.0 / min(max(sigma_bottom, 0.5), 0.999))

    rho_a = ps / (R_DRY * tas)
    u_bottom = spd[:, -1, :, :]
    ua_col = np.average(ua[:, lev >= cfg["transport"]["steering_sigma"], :, :]
                        if np.any(lev >= cfg["transport"]["steering_sigma"])
                        else ua[:, -3:, :, :], axis=1)
    va_col = np.average(va[:, lev >= cfg["transport"]["steering_sigma"], :, :]
                        if np.any(lev >= cfg["transport"]["steering_sigma"])
                        else va[:, -3:, :, :], axis=1)

    f_eff = drag_efficiency(z0, cfg)
    clay_pct = np.nan_to_num(clay, nan=0.0) * 100.0
    f_clay = np.clip(np.nan_to_num(clay, nan=0.0), 0.0, 0.2)   # K14 caps fclay at 0.2

    # Gravimetric soil moisture, in percent. mrso is a depth of water; the
    # pedology bulk density converts it against the regolith it sits in.
    bulk_density = 1500.0
    depth_m = 0.5
    grav_pct = np.clip(mrso / max(depth_m, 1e-6) * 1000.0 / bulk_density * 100.0,
                       0.0, 100.0)

    frac_bin, d_bin = emitted_mass_fractions(cfg)
    rho_p = cfg["removal"]["particle_density_kg_m3"]
    em_cfg, src = cfg["emission"], cfg["source"]

    mode = cfg["drag_partition"].get("mode", "bracket")

    # The bracket is over the aeolian roughness of the erodible surface, which
    # is the parameter that actually carries the uncertainty. See the long note
    # in aeolian/config/dust.yaml: the grid-cell roughness is NOT used here.
    dp = cfg["drag_partition"]
    # The wind-tail shape is measured from the snapshot climatology rather than
    # declared, because the declared value turned out to be wrong by enough to
    # move the answer two orders of magnitude. See the config note.
    snap = clim_path.with_name(clim_path.name.replace("_regular_", "_snapshot_"))
    k_measured = weibull_shape_from_snapshots(snap, erodible > 0.05)
    if k_measured is not None:
        cfg["subgrid_wind"]["weibull_shape"] = k_measured
    ends = [("low", dp["aeolian_z0_bracket_m"][0]),
            ("central", dp["aeolian_z0_m"]),
            ("high", dp["aeolian_z0_bracket_m"][1])]

    w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, nlon))
    is_land = land_fraction > 0.5

    def gmean(field, mask=None):
        ww = w if mask is None else w * mask
        return float((field * ww).sum() / max(ww.sum(), 1e-30))

    optics = json.loads((PROJECT_ROOT / "analysis" / "dust_optics.json").read_text(
        encoding="utf-8"))
    band1 = next(r for r in optics["results"]
                 if r["indices"] == cfg["optics"]["band1_indices"]
                 and r["band"] == "band 1")
    band2 = next(r for r in optics["results"]
                 if r["indices"] == cfg["optics"]["band2_indices"]
                 and r["band"] == "band 2")
    f1 = float(optics["stellar_flux_fraction_band1"])
    # m2/g to m2/kg, flux-weighted across the two bands.
    mee = 1000.0 * (f1 * band1["mass_extinction_efficiency_m2_g"]
                    + (1 - f1) * band2["mass_extinction_efficiency_m2_g"])

    radius = 6.371e6 * 1.2
    planet_area = 4 * np.pi * radius ** 2
    outcomes, fields = {}, {}
    for shelter, z0a in ends:
        emission, load, converged = run_one(
            good, z0a, cfg, nbin, nlat, nlon, frac_bin, d_bin, rho_p, u_bottom,
            z_ref, rho_a, tas, grav_pct, clay_pct, f_clay, erodible,
            land_fraction, snd, pr, prc, prl, ua_col, va_col, lat, lon, gravity,
            args.variant)

        column = load.sum(axis=0)
        aod_annual = (column * mee)[good].mean(axis=0)
        deposition = np.zeros((nlat, nlon))
        for t in good:
            precip_mm_hr = (prc[t] + prl[t]) * 1000.0 * 3600.0
            wet = (cfg["removal"]["scavenging_a"]
                   * np.maximum(precip_mm_hr, 0.0) ** cfg["removal"]["scavenging_b"])
            for b, db in enumerate(d_bin):
                vs = settling_velocity(db, rho_p, rho_a[t], tas[t], gravity)
                deposition += load[b, t] * (
                    vs / cfg["transport"]["dust_scale_height_m"] + wet) / len(good)
        dep = deposition * EARTH_YEAR_S * 1000.0            # g/m2 per Earth year
        emit_mean = gmean(emission[good].mean(axis=0))
        outcomes[shelter] = {
            "aeolian_z0_m": z0a,
            "global_emission_Tg_per_earth_year":
                round(emit_mean * EARTH_YEAR_S * planet_area / 1e9, 3),
            "land_mean_aod": round(gmean(aod_annual, is_land), 5),
            "global_mean_aod": round(gmean(aod_annual), 5),
            "max_aod": round(float(aod_annual.max()), 4),
            "land_mean_deposition_g_m2_per_earth_year": round(gmean(dep, is_land), 4),
            "max_deposition_g_m2_per_earth_year": round(float(dep.max()), 3),
            "transport_converged": bool(converged),
        }
        fields[shelter] = {"aod": aod_annual, "deposition": dep,
                           "emission": emission[good].mean(axis=0)}

    lo = min(o["land_mean_aod"] for o in outcomes.values())
    hi = max(o["land_mean_aod"] for o in outcomes.values())
    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "variant": args.variant,
        "note": "Offline dust: emission, transport and deposition from a "
                "climatology rather than inside the climate model. notes/dust.md "
                "says why; aeolian/README.md says what the numbers mean. The two "
                "shelter ends bracket how much the surrounding terrain shields "
                "an erodible patch, which this component does not claim to know.",
        "source_build": config.get("source_build"),
        "terrain_hash": terrain,
        "climatology": rel(clim_path),
        "climatology_sha256": sha256(clim_path),
        "config_sha256": sha256(args.config),
        "dust_config_sha256": sha256(args.dust_config),
        "lake_solution": rel(lakes),
        "excluded_time_bins": bad_bins,
        "excluded_because": "near-surface wind out of line with the other bins "
                            "by more than 2.5x; see flag_anomalous_bins",
        "grid": {"resolution": resolution, "latitudes": nlat,
                 "longitudes": nlon, "time_bins": nbin},
        "size_bins": [
            {"edges_um": cfg["size_distribution"]["bin_edges_um"][i:i + 2],
             "mass_fraction_emitted": round(float(frac_bin[i]), 4),
             "mean_diameter_um": round(float(d_bin[i]), 3)}
            for i in range(len(frac_bin))
        ],
        "source_map": {
            "erodible_weights": per_class,
            "erodible_fraction_of_land": round(gmean(erodible, is_land), 5),
            "land_fraction": round(gmean(land_fraction), 4),
        },
        "mass_extinction_efficiency_m2_kg": round(mee, 2),
        "subgrid_wind": {
            "weibull_shape_used": cfg["subgrid_wind"]["weibull_shape"],
            "measured_from": rel(snap) if k_measured is not None else None,
            "note": "Measured from the 32 instantaneous samples of the snapshot "
                    "climatology, not declared. It is an UPPER bound on the "
                    "shape and so a LOWER bound on emission, because samples "
                    "5.7 days apart resolve synoptic but not sub-daily "
                    "variance. Emission spans two orders of magnitude across "
                    "the plausible range of this one parameter, which is the "
                    "headline uncertainty of this component.",
        },
        "shelter_bracket": outcomes,
        "reopening_test": {
            "from": "notes/dust.md",
            "land_mean_aod_threshold": 0.10,
            "forcing_threshold_w_m2": 1.5,
            "land_mean_aod_range": [lo, hi],
            "crosses": bool(hi > 0.10),
            "forcing_is_DUST_2": "this script does not price the forcing",
        },
        "caveats": [
            "No dust-climate feedback: the climatology is an input and does not "
            "respond. That is what makes this defensible offline and it is also "
            "its limit.",
            "Vegetation cover is not modelled. Non-barren land is assumed to "
            "carry a canopy and not emit, which is the project's declared "
            "position and is generous on a world whose median land runoff is a "
            "few mm per Earth year. The `arid_bare_ground` variant brackets it.",
            "Dust is a well-mixed column of declared scale height advected by a "
            "single steering wind, not a three-dimensional tracer.",
            "The evaporite erodible weight is a declared suppression for crust "
            "cementation, not a measured efficiency, and it is the widest "
            "bracket in aeolian/config/dust.yaml.",
            "The shelter bracket spans a factor that is not a small correction: "
            "feeding MB95 the cell-mean roughness shuts emission off almost "
            "everywhere, and treating the patch as exposed removes shelter that "
            "is physically real.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    field_path = output.with_suffix(".nc")
    with Dataset(field_path, "w") as ds:
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        ds.title = "Offline dust: annual mean optical depth and deposition"
        ds.terrain_hash = terrain
        ds.source_build = str(config.get("source_build"))
        ds.climatology = rel(clim_path)
        ds.variant = args.variant
        for name, data in (("lat", lat), ("lon", lon)):
            v = ds.createVariable(name, "f8", (name,))
            v[:] = np.asarray(data)
        v = ds.createVariable("erodible_fraction", "f8", ("lat", "lon"), zlib=True)
        v.units, v.long_name = "1", "fraction of the cell that is bare erodible soil"
        v[:] = np.asarray(erodible)
        for shelter, f in fields.items():
            for name, data, units, note in [
                ("aod", f["aod"], "1", "annual mean dust optical depth"),
                ("deposition", f["deposition"], "g m-2 yr-1",
                 "annual mean dust deposition, per Earth year"),
                ("emission", f["emission"], "kg m-2 s-1",
                 "annual mean vertical dust flux"),
            ]:
                v = ds.createVariable(f"{name}_{shelter}", "f8", ("lat", "lon"),
                                      zlib=True)
                v.units, v.long_name = units, f"{note}, {shelter} end"
                v[:] = np.asarray(data)

    print(f"variant            {args.variant}")
    print(f"erodible land      "
          f"{100*report['source_map']['erodible_fraction_of_land']:.2f}% of land")
    for shelter, o in outcomes.items():
        print(f"  [{shelter:9s}] emission {o['global_emission_Tg_per_earth_year']:9.2f} Tg/yr"
              f"   land AOD {o['land_mean_aod']:.4f}   max {o['max_aod']:.3f}"
              f"   dep {o['land_mean_deposition_g_m2_per_earth_year']:.3f} g/m2/yr"
              f"   converged {o['transport_converged']}")
    print(f"reopening test     land-mean AOD {lo:.4f} to {hi:.4f} against 0.10 -> "
          f"{'CROSSES' if hi > 0.10 else 'below'}")
    print(f"\nwrote {rel(output)}\n      {rel(field_path)}")


def run_one(good, z0_aeolian, cfg, nbin, nlat, nlon, frac_bin, d_bin, rho_p, u_bottom,
            z_ref, rho_a, tas, grav_pct, clay_pct, f_clay, erodible,
            land_fraction, snd, pr, prc, prl, ua_col, va_col, lat, lon, gravity,
            variant):
    """One point in the roughness bracket: emit, transport, return the load."""
    em_cfg, src = cfg["emission"], cfg["source"]
    total_emission = np.zeros((nbin, nlat, nlon))
    load = np.zeros((len(frac_bin), nbin, nlat, nlon))
    converged = True
    for t in good:
        # The erodible patch is its own surface: the reference wind drives it
        # through ITS roughness, not the grid cell's. Shelter by surrounding
        # terrain is not represented and biases this high; see the config.
        u_star_patch = VON_KARMAN * u_bottom[t] / np.log(
            np.maximum(z_ref[t], 2.0) / z0_aeolian)
        # MB95 still partitions stress between the bed and whatever roughness
        # the patch itself carries, which is the ratio of its own z0 to the
        # smooth-bed value. A patch at the smooth limit keeps all of it.
        u_star_soil = u_star_patch * drag_efficiency(
            np.full_like(u_star_patch, z0_aeolian), cfg)

        u_st = em_cfg["u_star_st_typical_m_s"]
        u_t_dry = u_st * np.sqrt(em_cfg["rho_a0_kg_m3"] / rho_a[t])
        u_t = u_t_dry * moisture_threshold_factor(grav_pct[t], clay_pct, cfg)
        u_star_st = u_t * np.sqrt(rho_a[t] / em_cfg["rho_a0_kg_m3"])

        bare = erodible.copy()
        if variant == "arid_bare_ground":
            ab = src["arid_bare_ground"]
            dry = (pr[t] * EARTH_YEAR_S * 1000.0
                   < ab["precipitation_mm_per_earth_year"])
            bare = np.maximum(bare, np.where(dry, land_fraction * ab["weight"], 0.0))
        bare = np.where(snd[t] > src["snow_suppression_depth_m"], 0.0, bare)

        flux = emission_over_weibull(np.maximum(u_star_soil, 1e-6), u_t,
                                     rho_a[t], f_clay, u_star_st, cfg)
        total_emission[t] = flux * bare

        precip_mm_hr = (prc[t] + prl[t]) * 1000.0 * 3600.0
        wet_rate = (cfg["removal"]["scavenging_a"]
                    * np.maximum(precip_mm_hr, 0.0) ** cfg["removal"]["scavenging_b"])
        for b, (fb, db) in enumerate(zip(frac_bin, d_bin)):
            vs = settling_velocity(db, rho_p, rho_a[t], tas[t], gravity)
            loss = vs / cfg["transport"]["dust_scale_height_m"] + wet_rate
            m, _steps, ok = advect_to_steady_state(
                total_emission[t] * fb, ua_col[t], va_col[t], loss, lat, lon, cfg)
            load[b, t] = m
            converged &= ok
    return total_emission, load, converged


if __name__ == "__main__":
    main()
