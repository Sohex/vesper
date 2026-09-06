#!/usr/bin/env python3
"""Which lever moves this world's dust intensity, and by how much.

    python aeolian/scripts/dust_intensity_levers.py

world-03x. Menut et al. (2013) qualifies the per-lithology roughness this
component adopted in world-c8m in two ways, and both are claims about levers
rather than about constants:

  1. A TABULATED roughness biases dust flux high, because tabulated values are
     "more discrete and thus less variable and realistic" than a satellite
     retrieval. This component tabulates z0 by lithology, which is the
     configuration that finding is about.
  2. For a given soil data set two roughness data sets differ mainly in the
     SPATIAL DISTRIBUTION of the flux, while for a given roughness two soil
     texture data sets differ mainly in its INTENSITY.

Neither claim can be imported. They were measured on Earth, between two real
roughness maps and two real texture maps, and what this world has instead is one
map of each with a declared bracket on the level. So this script measures the
same decomposition here, on this project's own artifacts, and the answers are
different in both cases.

WHAT IT MEASURES, and every one of them is a ratio rather than a level, so the
declared constants that are common to both sides cancel:

  pattern_vs_intensity        does the roughness bracket move WHERE the dust is
                              or HOW MUCH there is
  tabulation_bias             emission with z0 distributed WITHIN a lithology,
                              against emission at the one tabulated value
  class_mixture_collapse      the same question for the between-class mosaic,
                              which `source_fractions` already collapses to one
                              geometric mean per cell
  roughness_bracket_decomposition  which lithology's z0 carries the bracket
  texture_sensitivity         emission under this project's own Earth-calibrated
                              clay field, and under its validation scatter
  drag_partition_worth        what using a drag partition at all is worth here,
                              which is what to quote beside a model that omits
                              one
  export_geometry_scale       whether the export's geometry can supply a
                              within-class roughness at all, which is the route
                              world-03x proposed and this rules out
  in_model_scalar_z0          what the scalar DUSTZ0 costs the in-model arm
                              against the offline mosaic (world-h24h)
  vegetation_bracket          what `--variant arid_bare_ground` is worth, which
                              is what the in-model arm cannot express
                              (world-4qem)

TRANSPORT IS NOT RUN. Every number here is a ratio of emission totals, and the
transport that turns emission into optical depth is linear in the source with
the climatology's winds held fixed, so a ratio of emissions is a ratio of
optical depths to the accuracy the `emission_to_optical_depth_transfer` block
measures on the baseline artifact's own three arms. That check is what licenses
the omission; it is reported rather than assumed.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from scipy.stats import norm

from _paths import ANALYSIS, CONFIG, DUST_CONFIG, PROJECT_ROOT  # noqa: E402

import climatology                                           # noqa: E402  from lib/, via _paths
import build_dust as bd                                      # noqa: E402
from builds import component_data, grid_export, mesh_export, soilmap
from gridding import gaussian_area_weights, region_cells
from orogen import LAND, Export
from paths import rel

BASELINE_JSON = ANALYSIS / "dust_baseline.json"
BASELINE_NC = ANALYSIS / "dust_baseline.nc"
EARTH_VALIDATION = PROJECT_ROOT / "pedology" / "analysis" / "earth_validation.json"

# THE CRITERIA, FIXED BEFORE ANY OF THE RATIOS BELOW WERE COMPUTED. A criterion
# chosen after the run it judges is not a criterion.
#
# The question is whether the intensity this component reports is an artefact of
# tabulating z0 by lithology. Two numbers already on the record set the scale it
# has to reach to be one: world-c8m moved the central land-mean optical depth by
# a factor of 4.4 when it replaced one scalar with two tabulated values, and the
# roughness bracket spans a factor of 10.6 on emission. So:
#
#   a tabulation bias under 1.25 in either direction cannot be the explanation,
#   because it is small against both;
#   one reaching 4.4 in the LOW direction would account for the whole world-c8m
#   jump and the intensity would be an artefact;
#   anything between is a partial account and is reported as a fraction of the
#   4.4 in the logarithm.
TABULATION_NEGLIGIBLE = 1.25
TABULATION_EXPLAINS = 4.4


def sha256_of(path: Path) -> str:
    return bd.sha256(path)


# -- the within-class roughness distribution ----------------------------------

def within_class_sigma(cfg: dict) -> dict:
    """The geometric spread of z0 WITHIN one arid surface class.

    Declared in `aeolian/config/dust.yaml` under `drag_partition.within_class_z0`
    and measured there from Prigent et al. (2005) Table 1, which is the only
    source in this repository that reports several in-situ z0 values for the
    SAME area. Read the config for the derivation; this reads the numbers.
    """
    wc = cfg["drag_partition"]["within_class_z0"]
    return {"central": float(wc["sigma_g"]),
            "bracket": [float(x) for x in wc["sigma_g_bracket"]]}


def lognormal_nodes(sigma_g: float, n: int) -> np.ndarray:
    """Multipliers on a class's central z0 at `n` equal-probability nodes."""
    return sigma_g ** norm.ppf((np.arange(n) + 0.5) / n)


def uniform_in_log_nodes(lo: float, hi: float, n: int) -> np.ndarray:
    """Multipliers of the class's own log-width, uniform in ln z0.

    The bounded alternative to the lognormal, and it is the one to prefer where
    the two disagree: a lognormal has no smooth end, so the integral keeps
    creeping up as the quadrature reaches further into a tail no measurement
    supports, while this one is bounded by the width the source actually spans.

    CENTRED ON THE ARM, not on the class's bracket midpoint. The arm already
    says where the class's central roughness sits; what is borrowed from the
    bracket here is its WIDTH in the logarithm and nothing else. Re-centring
    would move the level as well as spread it, and the level is what the arm is.
    """
    half = 0.5 * np.log(hi / lo)
    u = (np.arange(n) + 0.5) / n
    return np.exp((2.0 * u - 1.0) * half)


# -- the emission integral ----------------------------------------------------

class Emitter:
    """Total emission for a chosen roughness, texture and vegetation variant.

    One object holds the climatology and the source map so that a sweep costs a
    loop rather than a re-read. `populations` is a list of (z0_for_the_log_term,
    z0_for_the_drag_term, emitting_fraction) triples that are summed, which is
    what lets a cell carry a DISTRIBUTION of roughness rather than one value.
    """

    def __init__(self, config: dict, cfg: dict, clim_path: Path):
        self.config, self.cfg = config, cfg
        self.gravity = float(config["planet"]["gravity_m_s2"])
        with Dataset(clim_path) as ds:
            bin_centres = np.asarray(ds["time"][:], dtype=float)
            self.lat = np.asarray(ds["lat"][:], dtype=float)
            self.lon = np.asarray(ds["lon"][:], dtype=float)
            lev = np.asarray(ds["lev"][:], dtype=float)
            spd = np.asarray(ds["spd"][:], dtype=float)
            mrso = np.asarray(ds["mrso"][:], dtype=float)
            tas = np.asarray(ds["tas"][:], dtype=float)
            ps = np.asarray(ds["ps"][:], dtype=float) * 100.0
            self.snd = np.asarray(ds["snd"][:], dtype=float)
            self.pr = np.asarray(ds["pr"][:], dtype=float)
        self.nbin, self.nlat, self.nlon = tas.shape
        # The same corrupted-bin exclusion the baseline run applies. One bin
        # with 7.5x the wind is the whole annual total at u* cubed.
        self.bad = bd.flag_anomalous_bins(spd)
        self.good = [t for t in range(self.nbin) if t not in self.bad]
        # A lever is a RATIO of annual means, so both arms have to be annual
        # means: the bins hold different numbers of raw records and averaging
        # them alike is the climatology's weighting thrown away. Renormalised
        # over `good`, the way build_dust.py weights the same sum.
        w = climatology.bin_weights(bin_centres)
        self.bin_weight = np.zeros(self.nbin)
        self.bin_weight[self.good] = w[self.good] / w[self.good].sum()
        sigma_bottom = float(lev[-1]) if lev[-1] > lev[0] else float(lev[0])
        if sigma_bottom > 1.5:
            # AREA-weighted, like build_dust.py's copy of this branch: a plain
            # mean over a Gaussian grid counts a polar row and an equatorial
            # row alike.
            sigma_bottom = sigma_bottom * 100.0 / bd.area_mean(ps, self.lat)
        self.z_ref = (bd.R_DRY * tas / self.gravity) * np.log(
            1.0 / min(max(sigma_bottom, 0.5), 0.999))
        self.rho_a = ps / (bd.R_DRY * tas)
        self.u_bottom = spd[:, -1, :, :]
        self.grav_pct = bd.gravimetric_moisture_pct(mrso, cfg)
        # The soil map covers the cells it has soil for and leaves the rest NaN,
        # which the emission integral reads as zero clay and therefore zero
        # flux. The mask is kept apart from the filled field so that a texture
        # perturbation moves the clay the soil map ACTUALLY GAVE and does not
        # quietly turn cells it never covered into emitting ones.
        raw = bd.soil_clay_grid(soilmap(config), self.lat, self.lon)
        self.clay_known = np.isfinite(raw)
        self.clay = np.nan_to_num(raw, nan=0.0)
        cfg["emission"]["_gravity_scaling"] = bd.gravity_threshold_scaling(
            cfg["emission"], self.gravity)
        # The MEASURED subgrid wind shape the baseline ran at, not the config's
        # declared placeholder. The two differ by a factor of 40 in emission and
        # every ratio here has to be taken at the one the artifact carries.
        report = json.loads(BASELINE_JSON.read_text(encoding="utf-8"))
        cfg["subgrid_wind"]["weibull_shape"] = float(
            report["subgrid_wind"]["weibull_shape_used"])
        self.report = report
        self.weight = gaussian_area_weights(self.lat, self.nlon,
                                            what="the lever grid")

        lakes = component_data("hydrography", config, strict=True) / "surface_water.nc"
        (self.erodible, self.land_fraction, self.per_class, self.detail,
         self.terrain, self.z0_cell) = bd.source_fractions(config, cfg, lakes, None)
        self.lakes = lakes
        self.by_class = self._erodible_by_class(config, cfg, lakes)

    @staticmethod
    def _erodible_by_class(config: dict, cfg: dict, lakes: Path) -> dict:
        """Erodible area fraction per cell, per lithology, not summed.

        `source_fractions` returns the sum and one geometric-mean roughness.
        Resolving the mixture needs the classes apart, and this rebuilds them
        the same way: from `substrate_class` on `surface_class == LAND`, minus
        the solved lake extent. Land comes from `surface_class` and never from
        `land_mask`; the two disagree over exactly the dry closed-basin floor
        this component emits from.
        """
        export = Export(mesh_export(config))
        grid_dir = grid_export(config)
        lit = json.loads((export.root / "manifest.json").read_text(
            encoding="utf-8"))["lithology"]
        ids = {r["code"]: int(r["id"]) for r in lit["rockClasses"]}
        substrate = export.field("substrate_class")
        with Dataset(lakes) as ds:
            wet = np.asarray(ds["lake"][:]).astype(bool)
        is_land = export.surface_class == LAND
        cell, nlat, nlon = region_cells(export, grid_dir)
        area = export.cell_area.astype(np.float64)
        total = np.zeros(nlat * nlon)
        np.add.at(total, cell, area)
        out = {}
        for code, w in cfg["source"]["erodible_weight"].items():
            mask = (substrate == ids[code]) & is_land & ~wet
            acc = np.zeros(nlat * nlon)
            np.add.at(acc, cell[mask], area[mask] * float(w))
            out[code] = (acc / np.maximum(total, 1e-30)).reshape(nlat, nlon)
        return out

    def class_z0(self, code: str, end: str) -> float:
        by_class = self.cfg["drag_partition"]["aeolian_z0_by_class_m"]
        entry = by_class[code]
        return float({"low": entry["bracket"][0], "central": entry["z0"],
                      "high": entry["bracket"][1]}[end])

    def mosaic(self, end: str) -> np.ndarray:
        """The per-cell geometric-mean roughness of the class mosaic."""
        return np.asarray(self.z0_cell)[{"low": 0, "central": 1, "high": 2}[end]]

    def total(self, populations, clay=None, variant="baseline") -> float:
        """Area-weighted global emission, summed over the sub-populations."""
        # The arid-bare-ground term is a FLOOR on the emitting fraction of the
        # whole cell, not of one sub-population, so it cannot be applied inside
        # a sum over sub-populations without being counted once per term. It is
        # only ever asked for on the single-population case, and this refuses
        # the combination rather than returning a plausible wrong number.
        if variant != "baseline" and len(populations) != 1:
            raise ValueError(
                f"variant {variant!r} floors the cell's emitting fraction and "
                f"cannot be resolved over {len(populations)} sub-populations")
        cl = self.clay if clay is None else clay
        clay_pct = cl * 100.0
        f_clay = np.clip(cl, 0.0, float(self.cfg["emission"]["f_clay_max"]))
        em, src = self.cfg["emission"], self.cfg["source"]
        out = 0.0
        for t in self.good:
            acc = np.zeros((self.nlat, self.nlon))
            u_st = em["u_star_st_typical_m_s"] * em["_gravity_scaling"]
            u_t_dry = u_st * np.sqrt(em["rho_a0_kg_m3"] / self.rho_a[t])
            u_t = u_t_dry * bd.moisture_threshold_factor(
                self.grav_pct[t], clay_pct, self.cfg)
            u_star_st = u_t * np.sqrt(self.rho_a[t] / em["rho_a0_kg_m3"])
            bare_base = np.ones((self.nlat, self.nlon))
            if variant == "arid_bare_ground":
                ab = src["arid_bare_ground"]
                dry = (self.pr[t] * bd.EARTH_YEAR_S * 1000.0
                       < ab["precipitation_mm_per_earth_year"])
                extra = np.where(dry, self.land_fraction * ab["weight"], 0.0)
            for z0_log, z0_drag, frac in populations:
                u_star = bd.VON_KARMAN * self.u_bottom[t] / np.log(
                    np.maximum(self.z_ref[t], 2.0) / z0_log)
                u_star = u_star * bd.drag_efficiency(
                    np.broadcast_to(np.asarray(z0_drag, dtype=float),
                                    u_star.shape), self.cfg)
                flux = bd.emission_over_weibull(
                    np.maximum(u_star, 1e-6), u_t, self.rho_a[t], f_clay,
                    u_star_st, self.cfg)
                bare = frac * bare_base
                if variant == "arid_bare_ground":
                    bare = np.maximum(bare, extra)
                bare = np.where(self.snd[t] > src["snow_suppression_depth_m"],
                                0.0, bare)
                acc += flux * bare
            out += self.bin_weight[t] * float((acc * self.weight).sum())
        return out

    def point(self, end: str, **kw) -> float:
        """The treatment the component ships: one geometric-mean z0 per cell."""
        m = self.mosaic(end)
        return self.total([(m, m, self.erodible)], **kw)

    def class_resolved(self, end: str, **kw) -> float:
        """The class mosaic resolved instead of collapsed to its mean."""
        return self.total([(self.class_z0(c, end), self.class_z0(c, end),
                            self.by_class[c]) for c in self.by_class], **kw)

    def distributed(self, end: str, nodes_for, n: int, **kw) -> float:
        """z0 distributed WITHIN each class as well as resolved between them."""
        pops = []
        for c in self.by_class:
            centre = self.class_z0(c, end)
            for mult in nodes_for(c):
                z0 = centre * mult
                pops.append((z0, z0, self.by_class[c] / n))
        return self.total(pops, **kw)


# -- the blocks ---------------------------------------------------------------

def pattern_vs_intensity() -> dict:
    """Does the roughness bracket move WHERE the dust is, or HOW MUCH.

    Menut's decomposition puts roughness on the spatial side. That was measured
    between two roughness DATA SETS whose maps differ; this project's bracket is
    a level uncertainty applied to one map, so there is no reason to expect it
    to carry the same meaning, and it does not.
    """
    with Dataset(BASELINE_NC) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        nlon = ds.dimensions["lon"].size
        emission = {e: np.asarray(ds[f"emission_{e}"][:])
                    for e in ("low", "central", "high")}
        aod = {e: np.asarray(ds[f"aod_{e}"][:])
               for e in ("low", "central", "high")}
    w = gaussian_area_weights(lat, nlon, what=str(BASELINE_NC))
    total = {e: float((v * w).sum()) for e, v in emission.items()}
    out = {"emission_bracket_factor": total["low"] / total["high"],
           "pairs": {}}
    for a, b in (("low", "central"), ("central", "high"), ("low", "high")):
        pa, pb = emission[a] / total[a], emission[b] / total[b]
        m = (pa > 0) | (pb > 0)
        out["pairs"][f"{a}_vs_{b}"] = {
            "pattern_correlation": float(np.corrcoef(pa[m], pb[m])[0, 1]),
            "displaced_fraction_of_emission":
                float(0.5 * np.abs((pa - pb) * w).sum()),
        }
    out["reading"] = (
        "Across the whole roughness bracket the normalised emission pattern is "
        "nearly fixed while the total moves by the bracket factor above, so on "
        "this world the aeolian roughness is an INTENSITY lever and not a "
        "spatial one. Menut's decomposition compared two roughness maps that "
        "differ in space; a level bracket on one map cannot behave that way, "
        "and the finding does not transfer.")
    # The transfer that licenses running emission without transport.
    aod_land = {e: float(np.nanmean(v[emission[e] > 0])) for e, v in aod.items()}
    transfer = {}
    for a, b in (("low", "central"), ("central", "high")):
        transfer[f"{a}_over_{b}"] = {
            "emission_ratio": total[a] / total[b],
            "optical_depth_ratio": aod_land[a] / aod_land[b],
        }
    out["emission_to_optical_depth_transfer"] = {
        "ratios": transfer,
        "worst_disagreement": max(
            abs(v["emission_ratio"] / v["optical_depth_ratio"] - 1.0)
            for v in transfer.values()),
        "note": "Transport is linear in the source at fixed winds, so a ratio "
                "of emissions is a ratio of optical depths. Measured here on "
                "the baseline artifact's own three arms so that the ratios "
                "elsewhere in this file can be read as optical depth without "
                "re-running the transport.",
    }
    return out


def tabulation_bias(em: Emitter, sigma: dict, quadrature: int) -> dict:
    """Emission with z0 distributed within a class, against the tabulated value.

    This is the quantitative content of Menut's caution. Tabulating a class is
    exactly the collapse of a within-class distribution to a point, and emission
    is not linear in z0, so the collapse has a sign and a size.
    """
    out = {"sigma_g": sigma, "quadrature_points": quadrature, "ends": {}}
    zc = em.cfg["drag_partition"]["aeolian_z0_by_class_m"]
    for end in ("low", "central", "high"):
        base = em.class_resolved(end)
        arms = {}
        for label, sg in (("central", sigma["central"]),
                          ("smallest_spread", sigma["bracket"][0]),
                          ("largest_spread", sigma["bracket"][1])):
            nodes = lognormal_nodes(sg, quadrature)
            arms[label] = em.distributed(
                end, lambda _c, nodes=nodes: nodes, quadrature) / base
        # The bounded alternative, as wide in the logarithm as the class's own
        # measured range and centred on the arm.
        def bounded(code, end=end):
            e = zc[code]
            return uniform_in_log_nodes(float(e["bracket"][0]),
                                        float(e["bracket"][1]), quadrature)
        arms["bounded_uniform_in_log"] = em.distributed(
            end, bounded, quadrature) / base
        out["ends"][end] = arms
    lo = min(out["ends"]["central"].values())
    hi = max(out["ends"]["central"].values())
    out["central_end_bracket"] = [lo, hi]
    out["direction"] = "raises emission" if lo > 1.0 else (
        "lowers emission" if hi < 1.0 else "either sign")
    if lo >= 1.0:
        # The caution being tested says the tabulated answer is too HIGH, so a
        # ratio above one is the wrong sign for it whatever its size.
        out["verdict"] = ("the wrong sign: resolving the tabulation RAISES the "
                          "emission, so the tabulation cannot be why the "
                          "intensity is high")
    elif 1.0 / lo >= TABULATION_EXPLAINS:
        out["verdict"] = "accounts for the whole world-c8m jump"
    elif max(hi, 1.0 / lo) < TABULATION_NEGLIGIBLE:
        out["verdict"] = "negligible against the roughness bracket"
    else:
        out["verdict"] = "a partial account, in the direction stated"
    out["criteria"] = {
        "negligible_below": TABULATION_NEGLIGIBLE,
        "explains_the_jump_at": TABULATION_EXPLAINS,
        "fixed_before": "the ratios in this block were computed",
    }
    return out


def class_mixture_collapse(em: Emitter) -> dict:
    """The between-class collapse `source_fractions` already performs.

    It replaces the two lithologies in a cell with one erodible-area weighted
    geometric mean. That is the same species of collapse as the tabulation, and
    it needs no declared input to undo, so it is worth knowing whether it should
    be undone.
    """
    out = {}
    for end in ("low", "central", "high"):
        out[end] = em.class_resolved(end) / em.point(end)
    out["reading"] = (
        "Resolving the two lithologies rather than collapsing them to one "
        "geometric mean per cell changes the emission by the fractions above. "
        "The collapse stays: a correction this size does not earn a second "
        "code path through the emission integral, and the number is here so "
        "that the choice is priced rather than assumed.")
    return out


def roughness_bracket_decomposition(em: Emitter) -> dict:
    """Which lithology's roughness carries the bracket."""
    full = {e: em.class_resolved(e) for e in ("low", "central", "high")}
    out = {"full_bracket_factor": full["low"] / full["high"], "held": {}}
    for held, swept in (("evaporite", "playa_clastic"),
                        ("playa_clastic", "evaporite")):
        vals = {}
        for end in ("low", "central", "high"):
            pops = []
            for code in em.by_class:
                z0 = em.class_z0(code, "central" if code == held else end)
                pops.append((z0, z0, em.by_class[code]))
            vals[end] = em.total(pops)
        out["held"][f"{held}_at_its_centre"] = {
            "bracket_factor": vals["low"] / vals["high"],
            "swept": swept,
        }
    out["reading"] = (
        "The bracket is almost entirely the roughness of the class carrying "
        "erodible weight 1.0. Narrowing the intensity therefore means "
        "narrowing that one class's z0; the other class's range can be "
        "reinterpreted without moving the answer.")
    return out


def texture_sensitivity(em: Emitter) -> dict:
    """What the soil texture is worth on the intensity, on this world.

    Menut puts the intensity on the texture side. Here it cannot be, and the
    reason is structural rather than a matter of degree: this world's erodible
    substrate is closed-basin fill and evaporite crust, both clay-rich, so
    Kok's clay fraction is at its declared cap over most of the emission and the
    flux is pinned there. The cap is the standard treatment and is not in
    question; what it implies for the lever is.
    """
    validation = json.loads(EARTH_VALIDATION.read_text(encoding="utf-8"))
    sites = validation["sites"]
    predicted = np.array([s["clay_predicted"] for s in sites])
    observed = np.array([s["clay_observed_30_60cm"] for s in sites])
    rmse = float(np.sqrt(((predicted - observed) ** 2).mean()))
    slope = float(validation["verdict"]["regression_slope_observed_on_predicted"])
    intercept = float(validation["verdict"]["regression_intercept"])

    base = em.point("central")
    cap = float(em.cfg["emission"]["f_clay_max"])

    def perturb(f):
        """Move the clay the soil map gave, and leave the cells it did not."""
        return np.where(em.clay_known, np.clip(f(em.clay), 0.0, 1.0), em.clay)

    arms = {
        "earth_regression_calibrated": perturb(lambda c: slope * c + intercept),
        "plus_validation_rmse": perturb(lambda c: c + rmse),
        "minus_validation_rmse": perturb(lambda c: c - rmse),
        "halved": perturb(lambda c: c * 0.5),
        "doubled": perturb(lambda c: c * 2.0),
    }
    out = {
        "pedology_validation": {
            "n_sites": len(sites), "clay_rmse": rmse,
            "regression_slope_observed_on_predicted": slope,
            "regression_intercept": intercept,
            "source": rel(EARTH_VALIDATION),
        },
        "emission_ratio": {k: em.point("central", clay=v) / base
                           for k, v in arms.items()},
    }
    with Dataset(BASELINE_NC) as ds:
        emission = np.asarray(ds["emission_central"][:])
        erodible = np.asarray(ds["erodible_fraction"][:])
    live = erodible > 0
    out["clay_cap"] = {
        "f_clay_max": cap,
        "emission_weighted_fraction_at_or_above_cap":
            float((emission[live] * (em.clay[live] >= cap)).sum()
                  / max(float(emission[live].sum()), 1e-300)),
        "erodible_area_fraction_at_or_above_cap":
            float((erodible[live] * (em.clay[live] >= cap)).sum()
                  / max(float(erodible[live].sum()), 1e-300)),
        "note": "Kok equation 3 takes the flux as linear in the clay fraction "
                "up to a declared cap, above which the sandblasting relation "
                "it inherits rises without bound. Where the cap binds the "
                "texture cannot move the flux through that term at all, and "
                "only the Fecan residual is left -- which runs the other way, "
                "because more clay holds more water before the threshold "
                "moves.",
    }
    return out


def drag_partition_worth(em: Emitter) -> dict:
    """What using a drag partition at all is worth, on this world.

    Menut et al. (2013) records that using one lowers fluxes by a factor of 2 to
    3 against not using one. That is an Earth number over Earth's roughness
    distribution, and it is only useful here as the thing to quote when these
    numbers are set beside a model that omits the partition, so it is measured
    on this world's own roughness rather than borrowed.

    The comparison is against a drag efficiency of exactly one, which is a
    surface whose only roughness is its own grains. It is not an alternative
    treatment: MB95's partition is the physics, and this is the size of the term.
    """
    out = {"ends": {}}
    for end in ("low", "central", "high"):
        mosaic = em.mosaic(end)
        with_partition = em.total([(mosaic, mosaic, em.erodible)])
        # `drag_efficiency` reads its argument as a roughness; the no-shelter
        # limit is the smooth-bed value itself, where the function returns 1.
        smooth_bed = float(em.cfg["drag_partition"]["z0s_cm"]) / 100.0
        without = em.total([(mosaic, smooth_bed, em.erodible)])
        out["ends"][end] = {"with_over_without": with_partition / without}
    out["reading"] = (
        "The partition lowers the emission by these factors against a bed that "
        "keeps all of its stress. Quote it whenever a number from this "
        "component is set beside a model that does not partition the stress; "
        "it is not a bracket, because the partition is the physics rather than "
        "an option.")
    return out


def export_geometry_scale(config: dict, cfg: dict, lakes: Path) -> dict:
    """Whether the export's geometry can supply a within-class roughness at all.

    world-03x proposed building the within-class distribution from this
    project's own geometry, on the grounds that a playa has a smooth interior
    and rougher margins. It cannot, and this is why, measured rather than
    argued. Two numbers settle it.

    A connected component of erodible substrate is NOT a playa. Taken on
    `substrate_class` restricted to `surface_class == LAND` -- never
    `land_mask`, and never `land_mask | is_endorheic`, which misses exactly the
    depressions too small for the basin catalogue that this question is about --
    the components run to basin-fill provinces, and the height of a region above
    its component's own floor is therefore not a height above a playa floor.

    And the scale gap is six orders of magnitude. The mesh edge is kilometres
    and Orogen designs terrain down to about twenty; the drag partition
    partitions stress between a bed and centimetre-scale roughness ELEMENTS.
    That is the category error `aeolian/config/dust.yaml` already records for
    the grid-cell roughness field, restated at region scale.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    export = Export(mesh_export(config))
    lit = json.loads((export.root / "manifest.json").read_text(
        encoding="utf-8"))["lithology"]
    ids = {r["code"]: int(r["id"]) for r in lit["rockClasses"]}
    substrate = export.field("substrate_class")
    with Dataset(lakes) as ds:
        wet = np.asarray(ds["lake"][:]).astype(bool)
    erodible = np.zeros(export.n_regions, dtype=bool)
    for code in cfg["source"]["erodible_weight"]:
        erodible |= substrate == ids[code]
    erodible &= export.surface_class == LAND
    erodible &= ~wet

    off, lst = export.adjacency
    src = np.repeat(np.arange(export.n_regions, dtype=np.int64), np.diff(off))
    dst = lst.astype(np.int64)
    keep = erodible[src] & erodible[dst]
    graph = csr_matrix((np.ones(int(keep.sum()), dtype=np.int8),
                        (src[keep], dst[keep])),
                       shape=(export.n_regions, export.n_regions))
    _n, label = connected_components(graph, directed=False)
    _uniq, label = np.unique(label[erodible], return_inverse=True)
    npatch = int(label.max()) + 1

    area = export.cell_area.astype(np.float64)[erodible]
    elevation = export.field("elevation_km").astype(np.float64)[erodible]
    floor = np.full(npatch, np.inf)
    np.minimum.at(floor, label, elevation)
    patch_area = np.zeros(npatch)
    np.add.at(patch_area, label, area)
    height_m = (elevation - floor[label]) * 1000.0

    def weighted_quantile(values, weights, q):
        order = np.argsort(values)
        v, w = values[order], weights[order]
        return float(np.interp(q, np.cumsum(w) / w.sum(), v))

    return {
        "mesh_mean_edge_km": float(
            export.manifest["basins"]["resolution"]["avgEdgeKm"]),
        "erodible_regions": int(erodible.sum()),
        "connected_patches": npatch,
        "patch_area_km2_area_weighted": {
            "median": weighted_quantile(patch_area[label], area, 0.5),
            "p99": weighted_quantile(patch_area[label], area, 0.99),
        },
        "height_above_patch_floor_m_area_weighted_median":
            weighted_quantile(height_m, area, 0.5),
        "erodible_area_within_10m_of_its_patch_floor":
            float(area[height_m <= 10.0].sum() / area.sum()),
        "reading": "The components are basin-fill provinces rather than "
                   "playas, so height above one's floor does not measure "
                   "position within a playa. And a mesh edge of kilometres "
                   "cannot carry the centimetre-scale roughness elements the "
                   "drag partition wants, which is why the tabulation bias "
                   "above depends on the SPREAD of the within-class "
                   "distribution and not on how the geometry arranges it: a "
                   "cell of this grid already contains the whole population.",
    }


def in_model_scalar_z0(em: Emitter) -> dict:
    """What the scalar DUSTZ0 costs the in-model arm. world-h24h.

    Boundary field 1802 carries the per-cell drag partition, so the mosaic
    reaches the model there. The friction velocity's `ln(zref/z0)` term does not:
    `dustsrc` divides by the scalar namelist `dustz0`. Two questions follow and
    they have very different answers, which is why they are measured apart:
    whether the scalar is the RIGHT scalar, and whether a scalar can do at all.
    """
    dp = em.cfg["drag_partition"]
    class_blind = {"low": float(dp["aeolian_z0_bracket_m"][0]),
                   "central": float(dp["aeolian_z0_m"]),
                   "high": float(dp["aeolian_z0_bracket_m"][1])}
    out = {"ends": {}}
    for end in ("low", "central", "high"):
        mosaic = em.mosaic(end)
        offline = em.total([(mosaic, mosaic, em.erodible)])
        scalar = bd.mosaic_scalar_z0(mosaic, em.erodible, em.lat)
        out["ends"][end] = {
            "mosaic_area_weighted_geometric_mean_m": scalar,
            "class_blind_scalar_m": class_blind[end],
            "in_model_over_offline_with_class_blind_scalar":
                em.total([(class_blind[end], mosaic, em.erodible)]) / offline,
            "in_model_over_offline_with_mosaic_mean":
                em.total([(scalar, mosaic, em.erodible)]) / offline,
        }
    out["reading"] = (
        "The class-blind scalar is not the mosaic's mean and never was: it is "
        "the pre-world-c8m bracket, whose ends sit outside the per-lithology "
        "ones. Feeding it to the log term puts the in-model arm off the offline "
        "one by the first ratio above, which is a comparison defeated before it "
        "starts. Setting DUSTZ0 to the mosaic's own area-weighted geometric "
        "mean leaves the second ratio, and that residual is the whole of what a "
        "fourth boundary field carrying z0 per cell would recover.")
    return out


def vegetation_bracket(em: Emitter) -> dict:
    """What `--variant arid_bare_ground` is worth. world-4qem.

    A static boundary field cannot carry a per-timestep precipitation test, so
    the in-model arm expresses the baseline variant only. What that gap costs is
    the size of the bracket end it cannot reach.
    """
    out = {"ends": {}}
    for end in ("low", "central", "high"):
        mosaic = em.mosaic(end)
        pops = [(mosaic, mosaic, em.erodible)]
        out["ends"][end] = (em.total(pops, variant="arid_bare_ground")
                            / em.total(pops))
    ab = em.cfg["source"]["arid_bare_ground"]
    out["threshold_mm_per_earth_year"] = float(
        ab["precipitation_mm_per_earth_year"])
    out["bare_bedrock_weight"] = float(ab["weight"])
    out["reading"] = (
        "The variant lets any cell drier than the declared threshold emit "
        "regardless of lithology, at the declared bare-bedrock weight. The "
        "ratios above are what the in-model arm cannot express, and they are "
        "one-signed: the in-model total is a floor by that much.")
    return out


# -- main ---------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--dust-config", type=Path, default=DUST_CONFIG)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="defaults to the climatology the baseline report was "
                         "built on, because every ratio here is taken against "
                         "that artifact's numbers")
    ap.add_argument("--quadrature", type=int, default=33,
                    help="nodes across the within-class roughness distribution")
    ap.add_argument("--output", type=Path,
                    default=ANALYSIS / "dust_intensity_levers.json")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cfg = yaml.safe_load(args.dust_config.read_text(encoding="utf-8"))
    cfg["_planet_radius_earth"] = config["planet"]["radius_earth"]

    report = json.loads(BASELINE_JSON.read_text(encoding="utf-8"))
    if report["source_build"] != str(config["source_build"]):
        raise SystemExit(
            f"{rel(BASELINE_JSON)} was built on {report['source_build']} and "
            f"config/planet.yaml says {config['source_build']}. Every ratio "
            "here is taken against that artifact; re-run build_dust.py first.")
    clim_path = args.climatology or (PROJECT_ROOT / report["climatology"])
    if not clim_path.exists():
        raise SystemExit(f"{rel(clim_path)} does not exist")

    em = Emitter(config, cfg, clim_path)
    sigma = within_class_sigma(cfg)

    out = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": bd.git_commit(),
        "note": "world-03x. Which lever moves this world's dust intensity. "
                "Every figure is a RATIO of emission totals; transport is not "
                "run and the emission_to_optical_depth_transfer block is what "
                "licenses that.",
        "source_build": str(config["source_build"]),
        "terrain_hash": em.terrain,
        "climatology": rel(clim_path),
        "baseline_report": rel(BASELINE_JSON),
        "baseline_report_sha256": sha256_of(BASELINE_JSON),
        "dust_config_sha256": sha256_of(args.dust_config),
        "excluded_time_bins": em.bad,
        "weibull_shape_used": cfg["subgrid_wind"]["weibull_shape"],
        "pattern_vs_intensity": pattern_vs_intensity(),
        "tabulation_bias": tabulation_bias(em, sigma, args.quadrature),
        "class_mixture_collapse": class_mixture_collapse(em),
        "roughness_bracket_decomposition": roughness_bracket_decomposition(em),
        "texture_sensitivity": texture_sensitivity(em),
        "drag_partition_worth": drag_partition_worth(em),
        "export_geometry_scale": export_geometry_scale(config, cfg, em.lakes),
        "in_model_scalar_z0": in_model_scalar_z0(em),
        "vegetation_bracket": vegetation_bracket(em),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")
    tb = out["tabulation_bias"]
    print(f"  tabulation bias, central end   "
          f"{tb['central_end_bracket'][0]:.3f} to "
          f"{tb['central_end_bracket'][1]:.3f}  ({tb['verdict']})")
    print(f"  roughness bracket on emission  "
          f"{out['pattern_vs_intensity']['emission_bracket_factor']:.2f}x")
    ts = out["texture_sensitivity"]["emission_ratio"]
    print(f"  texture, Earth-calibrated clay {ts['earth_regression_calibrated']:.3f}")
    print(f"  vegetation bracket, central    "
          f"{out['vegetation_bracket']['ends']['central']:.3f}")


if __name__ == "__main__":
    main()
