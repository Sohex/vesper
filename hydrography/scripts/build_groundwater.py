"""The steady-state water table under a real climate, and what it moves.

`groundwater.py` owns the discretisation and the solve; this supplies the
forcing, writes the products, and takes the measurement the carve criterion is
waiting for. One product, one measurement:

    data/<build>/water_table.nc     per region: head, depth, seepage, transmissivity
                                    per basin: net groundwater exchange

    python hydrography/scripts/build_groundwater.py

RECHARGE IS P-E, the same field `surface_water.py` runs the lakes on, read
through the same coupling convention. At steady state that is the right
quantity and not an approximation to one: whatever falls on land and does not
evaporate has to leave, and this model decides which of the two doors it leaves
by. A cell whose water table reaches the surface returns its recharge as
seepage, which IS the surface runoff the existing balance already routes; a cell
whose water table is below the surface passes some of it sideways instead. That
is the whole difference between this and the surface-only balance, and it is why
zero permeability has to reproduce that balance exactly.

WHAT THIS MEASURES AND DOES NOT DECIDE. `Qg`, the net groundwater exchange per
basin, is the term Fan (2019) writes into the catchment budget and every water
balance in this project currently drops. It is reported here, per basin, with
the count of carve verdicts it would move. It is NOT wired into
`carve_verdict.py`: that is a decision about loop A and it is taken by a person
reading this measurement, not by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402

import carve_verdict as cv  # noqa: E402
import groundwater as gw  # noqa: E402
import lake_balance as lb  # noqa: E402
import surface_water as sw  # noqa: E402
from orogen import LAND, Export  # noqa: E402
from paths import best_available_climatology, rel  # noqa: E402

import builds  # noqa: E402
import orbit  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "config" / "groundwater.yaml"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def rock_permeability(export: Export, cfg: dict, sigma: float = 0.0):
    """Intrinsic permeability per region, m2, and the unassigned mask.

    Read through the export's own `rockClasses` legend rather than an integer
    table held here, for the reason `lib/surface_classes.py` gives: a class
    inserted upstream would silently shift every code below it, and a consumer
    holding its own copy keeps working while pointing at the wrong rock.

    `sigma` shifts every class by that many of ITS OWN standard deviations, so
    +1 and -1 are the bracket Gleeson's spread implies rather than a uniform
    factor. The classes do not share a sigma and must not be shifted as if they
    did: carbonate's is 1.5 and siliciclastic's is 2.5.
    """
    legend = {rc["code"]: int(rc["id"])
              for rc in export.manifest["lithology"]["rockClasses"]}
    hydro = cfg["hydrolithology"]
    mapping = cfg["rock_class_hydrolithology"]

    unknown = set(mapping) - set(legend)
    if unknown:
        raise SystemExit(
            f"groundwater.yaml names rock classes this export does not have: "
            f"{sorted(unknown)}")
    missing = set(legend) - set(mapping)
    if missing:
        raise SystemExit(
            f"this export has rock classes groundwater.yaml does not map: "
            f"{sorted(missing)}. Every class needs a hydrolithology or an "
            f"explicit null; there is no default.")

    n_class = max(legend.values()) + 1
    log_k = np.full(n_class, np.nan)
    for code, name in mapping.items():
        if name is None:
            continue
        h = hydro[name]
        log_k[legend[code]] = float(h["log10_k_m2"]) + sigma * float(h["log10_k_sigma"])

    rock = export.substrate_class.astype(np.int64)
    k = np.where(np.isnan(log_k[rock]), 0.0, 10.0 ** np.nan_to_num(log_k[rock]))
    unassigned = np.isnan(log_k[rock])
    return k, unassigned, legend


def aquifer_thickness(export: Export, cfg: dict, land: np.ndarray):
    """Saturated thickness per region, metres, and how it was obtained. GW-18.

    `constant` is Gleeson's "on the order of 100 m" everywhere, which is what
    this component has always run and is a multiplicative constant on every
    transmissivity.

    `cover_thickness` reads the export's own surviving cover thickness and
    floors it at the constant. The floor is not a fudge: the constant is the
    depth Gleeson's lithology maps describe, that depth exists wherever there is
    rock, and what the cover adds is basin FILL on top of a shallow subsurface
    already there. Exposed basement keeps the floor rather than dropping out of
    the conductive network on a zero.

    AND IT IS A THICKNESS, NOT A HISTORY. `docs/src/reference/no-time-axis.md`
    forbids asking this terrain for an age or an accumulation rate. This asks it
    for neither: `cover_thickness` is the thickness that is there now.
    """
    source = str(cfg["aquifer"].get("thickness_source", "constant"))
    constant = float(cfg["aquifer"]["thickness_m"])
    if source == "constant":
        return np.full(export.n_regions, constant), {
            "source": "constant", "constant_m": constant}
    if source != "cover_thickness":
        raise SystemExit(
            f"aquifer.thickness_source {source!r} is not one this solver has; "
            "it is 'constant' or 'cover_thickness'")

    floor = float(cfg["aquifer"].get("thickness_floor_m", constant))
    cover_m = export.field("cover_thickness").astype(np.float64) * 1000.0
    if not np.all(np.isfinite(cover_m)):
        raise SystemExit(
            "the export's cover_thickness carries non-finite values; a "
            "thickness that is not a number cannot become a transmissivity")
    thick = np.maximum(cover_m, floor)
    on_land = thick[land]
    return thick, {
        "source": "cover_thickness",
        "floor_m": floor,
        "at_floor_land_fraction": float((cover_m[land] <= floor).mean()),
        "land_median_m": float(np.median(on_land)),
        "land_p95_m": float(np.percentile(on_land, 95)),
        "land_max_m": float(on_land.max()),
    }


def recharge_field(export: Export, config: dict, clim_path: Path):
    """P-E per region, m/s, on the mesh. The same field the lakes are solved on.

    Goes through `surface_water.climate_fields` and
    `surface_water.region_grid_cells` rather than repeating either. The second
    of those is the export/ExoPlaSim grid join, which `CLAUDE.md` rule 3 governs
    and which this project has got wrong on three separate scripts; there is one
    copy of it in `lib/gridding.py` and this is not a fourth.

    **THE SINK'S `ET_max` IS THE BIN MEAN OF THE PER-BIN PENMAN, not one Penman
    on annual-mean air.** Penman is nonlinear in everything it reads, so the two
    are different numbers, and this path used to take the second: it called
    `climate_fields` with no `bin_index` where `carve_verdict` and
    `surface_water` had both already moved. `land_water_ledger.yaml` holds
    `open_water_evaporation` at `interval_floor: climatology_bin`, so the annual
    evaluation was against a standing decision rather than a simplification.

    THE BIN MEAN IS THE RIGHT COMBINATION HERE AND NOT ONE END OF A BRACKET,
    which is the question `carve_verdict._INTERVAL_BRACKET` makes every caller
    of Penman answer. Two reasons, and they are separate. The sink is
    `E(d) = ET_max exp(-d / lambda)`, nonlinear in the DEPTH and linear in
    `ET_max`; a steady-state solve holds `d` fixed through the year, so the
    annual mean of `E` is `exp(-d / lambda)` times the annual mean of `ET_max`
    and that mean is the only thing the solver can be given. And the body doing
    the evaporating is the ground, not a lake: the seasonal heat storage that
    makes the annual evaluation the defensible end for a deep water body is a
    ground heat flux of a few W/m2 against a net radiation cycle an order
    larger, so the no-storage end is the near one rather than a bound.

    Both arms are returned. The annual one is no longer used for the solve; it
    is measured against the bin mean and the spread goes on the report, because
    the SIZE of a correction belongs on the artifact rather than in a note.
    """
    sw._CLIM_FILE = clim_path
    # Eight values, and the staged albedo is the one this path used to discard.
    # It is needed now: the per-bin Penman below reads the same staged field,
    # and taking it from here is what keeps the two arms on ONE albedo. Unpack
    # the whole tuple rather than a prefix of it, so a ninth return raises here
    # instead of silently shifting every name one place along.
    (lat, lon, runoff, precip, evap_annual, model_runoff, lsm,
     staged_albedo) = sw.climate_fields(config)
    # Through `carve_verdict.bin_mean_open_water`, which is the one place that
    # owns the per-bin Penman loop and its weights. A second loop written here
    # would be a fourth copy of a convention this component has already had to
    # reconcile twice.
    land_albedo = cv.read_sra_field(PROJECT_ROOT / staged_albedo["path"],
                                    *lsm.shape)
    # `clip_at_zero` because the consumer is a SINK: a water table cannot gain
    # water from a bin whose Penman is negative, so the floor belongs inside the
    # mean rather than after it. `max(x, 0)` is a rectifier, which is the shape
    # that makes an interval error large rather than second order.
    evaporation = cv.bin_mean_open_water(
        clim_path, land_albedo, float(config["planet"]["gravity_m_s2"]),
        cfg=config, clip_at_zero=True)
    row, col = sw.region_grid_cells(export, lat)
    return runoff[row, col], (lat, lon, runoff, precip, evaporation, lsm,
                              evap_annual)


def basin_exchange(terminal, seepage_m3_s, supply_m3_s, n_basins,
                   et_m3_s=None):
    """Net groundwater import per basin, m3/s. Fan (2019)'s `Qg`, signed.

    A basin's catchment receives its own recharge and returns some of it as
    seepage. The difference is what crossed the catchment boundary underground:
    positive where the basin gains water the surface balance never sees,
    negative where it loses water the surface balance credits it with.

    **THE EVAPOTRANSPIRATION HAS TO BE ADDED BACK, and this is not a
    correction, it is the identity.** Per cell the solver enforces

        seepage = supply - et + (net lateral inflow)

    so summed over a catchment, `sum(seepage) - sum(supply) + sum(et)` is the
    net lateral inflow across its boundary, which is what `Qg` means. Drop the
    ET term and every drop the sink evaporated is charged to the basin as
    groundwater EXPORT, which is a different physical claim about a different
    destination: exported water arrives somewhere else and evaporated water does
    not.

    Leaving it out inverted the answer rather than shading it. It put 2,853
    basins in deficit against 546 in surplus and a median absolute shift of 90%
    of a basin's own recharge, where the exchange is a redistribution that sums
    to nearly nothing. `et_m3_s` is None only when the sink is off, where the
    term is zero and this is the arithmetic it always was.
    """
    gain = np.zeros(n_basins)
    give = np.zeros(n_basins)
    evap = np.zeros(n_basins)
    sel = terminal >= 0
    np.add.at(gain, terminal[sel], seepage_m3_s[sel])
    np.add.at(give, terminal[sel], supply_m3_s[sel])
    if et_m3_s is not None:
        np.add.at(evap, terminal[sel], np.asarray(et_m3_s)[sel])
    return gain - give + evap, gain, give


def artifact_path(args, data: Path) -> Path:
    """Where the water table goes. A VARIANT NEVER TAKES THE PRODUCT'S NAME.

    `water_table.nc` is what `builds.component_data` resolves for every
    downstream consumer, so a GW-18 or GW-24 arm writing it would be read as
    the product rather than as the experiment it is. `--output` still wins,
    for a caller naming an artifact deliberately.
    """
    if args.output is not None:
        return args.output
    tag = ""
    if args.aquifer_thickness_source not in (None, "constant"):
        tag += f"_{args.aquifer_thickness_source}"
    if args.unconfined:
        tag += "_unconfined"
    return data / f"water_table{tag}.nc"


def report_path(args, out: Path) -> Path:
    """Where this run's report goes. One artifact, one report, no collisions.

    The artifact refuses to overwrite an existing file; the report used to have
    no such protection and a fixed name, so every variant run -- a sigma arm, an
    operator-noise member, a groundwater-ET arm -- wrote over the central
    `groundwater_report.json` on its way past. The unconverged branch was worse:
    it wrote the fixed name whatever the run was, so a variant that failed to
    converge replaced the record of the run that had.

    So the name follows the ARTIFACT when one is named, and otherwise carries
    the variant in it. The names below are unchanged for the default run and for
    the sigma and noise arms, because `notes/` and the closed issues cite them.
    """
    if args.output is not None:
        return ANALYSIS / f"{out.stem}_report.json"
    if args.reduction_test:
        return ANALYSIS / "groundwater_reduction_test.json"
    if args.divide_test:
        return ANALYSIS / "groundwater_divide_test.json"
    if args.operator_noise:
        return ANALYSIS / (f"groundwater_report_noise{args.operator_noise:g}"
                           f"_seed{args.noise_seed}.json")
    if args.sigma:
        return ANALYSIS / f"groundwater_report_sigma{args.sigma:+g}.json"
    if args.et_lambda is not None:
        return ANALYSIS / f"groundwater_report_et{args.et_lambda:g}.json"
    # GW-18 and GW-24 are variants by the same argument as the four above: each
    # changes the transmissivity and so the field, and neither was in this list
    # when it arrived, so a sourced-thickness or unconfined run without
    # `--output` wrote over the record of the run every note cites.
    tag = ""
    if args.aquifer_thickness_source not in (None, "constant"):
        tag += f"_{args.aquifer_thickness_source}"
    if args.unconfined:
        tag += "_unconfined"
    return ANALYSIS / f"groundwater_report{tag}.json"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, default=None,
                    help="hydrography products for this build; defaults to "
                         "data/<source_build>/")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="regular climatology to take the recharge from; "
                         "defaults to the configured bootstrap_climatology")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--sigma", type=float, default=0.0,
                    help="shift every hydrolithology by this many of its own "
                         "Gleeson standard deviations; -1 and +1 are the bracket")
    ap.add_argument("--reduction-test", action="store_true",
                    help="run at zero permeability and check the surface-only "
                         "identity instead of solving the real case")
    ap.add_argument("--divide-test", action="store_true",
                    help="run at uniform permeability with a terrain-following "
                         "table and check the catchments against the surface ones")
    ap.add_argument("--et-lambda", type=float, default=None,
                    help="override the e-folding depth of GW-15's ET sink, "
                         "metres. Default is config `evapotranspiration.lambda_m`. "
                         "ET_max is the Penman field the lakes and the carve "
                         "verdict take")
    ap.add_argument("--no-groundwater-et", action="store_true",
                    help="turn GW-15's sink OFF. It is on by default because it "
                         "is not optional physics: without it the table pins at "
                         "the surface over most of the land. This exists for the "
                         "reduction identity and for reproducing pre-GW-15 runs")
    ap.add_argument("--no-baselevels", action="store_true",
                    help="turn GW-17's river and lake fixed heads OFF, leaving "
                         "the ocean as the only boundary. On by default for the "
                         "same reason: a river IS the water table where it sits")
    ap.add_argument("--aquifer-thickness-source", default=None,
                    choices=["constant", "cover_thickness"],
                    help="GW-18: override config's aquifer.thickness_source. "
                         "`cover_thickness` takes the saturated thickness from "
                         "the export's own surviving cover, floored at the "
                         "constant, so the range is sourced rather than assumed")
    ap.add_argument("--unconfined", action="store_true",
                    help="GW-24: T = K (h - z_bottom) rather than T = K D, with "
                         "the aquifer base at surface minus the thickness. The "
                         "transmissivity then depends on the head, so the "
                         "problem stops being a LINEAR complementarity problem "
                         "and --uniqueness-check stops being an identity. Read "
                         "hydrography/notes/subgrid-water-table.md first")
    ap.add_argument("--max-outer", type=int, default=60)
    ap.add_argument("--operator-noise", type=float, default=0.0,
                    help="multiply every face's w/l by lognormal noise of this "
                         "relative width, to ask what the operator's own "
                         "truncation error (GW-8, about 0.12) does to the "
                         "answer. Not a model parameter")
    ap.add_argument("--noise-seed", type=int, default=0)
    ap.add_argument("--uniqueness-check", action="store_true",
                    help="re-solve from the opposite initial active set and "
                         "require the identical head field. The problem is a "
                         "linear complementarity problem with a symmetric "
                         "positive definite matrix, so its solution is unique "
                         "and this is an identity rather than a comparison")
    args = ap.parse_args()

    config = yaml.safe_load((PROJECT_ROOT / "config/planet.yaml").read_text())
    cfg = yaml.safe_load(CONFIG.read_text())

    data = args.data if args.data is not None else builds.component_data(
        "hydrography", config, strict=True)
    # THE BEST AVAILABLE, which is the baseline once one is named and the
    # bootstrap before that. Recharge is P-E off the climatology, so the water
    # table is a function of the climate STATE and the bootstrap is the best
    # answer only on the pass where it is the only one. The
    # `needs: bootstrap_climatology` edge in config/pipeline.yaml is unchanged
    # and states what must EXIST for a first pass to run. Resolved through
    # `lib/paths.py` rather than out of the config here: that is the one copy,
    # and it carries `require_configured_grid`, which a hand-rolled read of the
    # config key does not.
    clim, clim_stage = best_available_climatology(args.climatology)
    # NAMED BY WHAT WAS DECLARED, READ BY WHERE THE BYTES ARE, and the two are
    # not the same path in a worktree. `scripts/link_worktree.py` links the
    # climatology payload in from the main checkout, so `resolve()` follows the
    # link OUT of the project root and `rel()` can then only hand back an
    # absolute path. That path is what gets stamped on `water_table.nc` and
    # compared against `surface_water.nc`'s repo-relative `forcing` stamp
    # below, so a resolved name refuses every worktree run and writes an
    # artifact naming a file that does not exist in any other tree.
    clim_name = rel(args.climatology if args.climatology is not None else clim)
    clim = clim.resolve()
    if not clim.is_file():
        raise SystemExit(
            f"config/planet.yaml names {clim_name} as the {clim_stage} "
            "climatology and it does not exist")
    # THE TWO HALVES OF ONE WATER BALANCE, checked rather than assumed. This
    # takes the same recharge field `surface_water.py` runs the lakes on, so
    # the surface and the subsurface have to be forced by the same run; the two
    # scripts resolve it separately, and once one of them can prefer a baseline
    # they can differ. `surface_water.nc` names its own forcing on the
    # artifact, which is what makes this an identity to test rather than a
    # convention to remember.
    sw_nc = data / "surface_water.nc"
    if sw_nc.is_file():
        with Dataset(sw_nc) as _ds:
            sw_forcing = getattr(_ds, "forcing", None)
        # COMPARED AS FILES AND NOT AS SPELLINGS. The stamp is repo-relative,
        # so it is re-anchored and resolved before the test: one file reachable
        # under two names is one forcing, and two names for one file must not
        # read as a disagreement about the climate.
        if sw_forcing:
            sw_path = Path(sw_forcing)
            if not sw_path.is_absolute():
                sw_path = PROJECT_ROOT / sw_path
            if sw_path.resolve() != clim:
                raise SystemExit(
                    f"{rel(sw_nc)} was forced by {sw_forcing} and this run "
                    f"resolved {clim_name} as the {clim_stage} climatology. "
                    "The lakes and the water table are the surface and "
                    "subsurface halves of one balance and cannot come from "
                    "two climates: re-run `surface_water` on this "
                    "climatology, or pass --climatology to force both onto "
                    "the same one.")

    build = builds.build_root(config)
    export = Export(builds.mesh_export(config))
    sw._DATA = data
    land = export.surface_class == LAND

    print(f"build {build.name}, {export.n_regions:,} regions")
    print("reconstructing the Voronoi dual")
    geom = gw.Geometry(export)
    closes = geom.closes_on_sphere()
    print(f"  {geom.src.size:,} faces, area closes on the sphere to {closes:.8f}")
    if args.operator_noise:
        # GW-8 SENSITIVITY, not a model knob. The mesh operator sits about 12%
        # above its analytic eigenvalue past l = 1, which is first-order
        # truncation on an irregular mesh. Injecting that as multiplicative
        # noise on the face coefficients asks what it does to the answer -- in
        # particular whether the carve flip COUNT survives it, or only the
        # direction of the flips.
        rng = np.random.default_rng(args.noise_seed)
        sig = float(args.operator_noise)
        geom.geom = geom.geom * rng.lognormal(
            -0.5 * np.log1p(sig * sig), np.sqrt(np.log1p(sig * sig)),
            size=geom.geom.size)
        print(f"  operator noise: w/l perturbed by {sig:.0%} relative, "
              f"seed {args.noise_seed}")
    if abs(closes - 1.0) > 1e-6:
        raise SystemExit(
            f"the Voronoi dual does not tile the sphere ({closes:.8f}); the "
            "geometry is wrong and nothing below it is worth computing")

    # -- subsurface --------------------------------------------------------
    k_m2, unassigned, legend = rock_permeability(export, cfg, args.sigma)
    fluid = cfg["fluid"]
    gravity = float(config["planet"]["gravity_m_s2"])
    k0 = gw.conductivity(k_m2, float(fluid["density_kg_m3"]), gravity,
                         float(fluid["dynamic_viscosity_pa_s"]))
    if args.aquifer_thickness_source is not None:
        cfg["aquifer"]["thickness_source"] = args.aquifer_thickness_source
    thickness_m, thickness_note = aquifer_thickness(export, cfg, land)
    # GW-24. The unconfined form reads a per-region aquifer BASE, and the base is
    # the same field the confined form reads as a multiplier: `surface - D`. One
    # thickness field, two readings, and neither invents a second geometry.
    unconfined = (str(cfg["aquifer"].get("transmissivity", "confined"))
                  == "unconfined") or args.unconfined
    min_saturated = float(cfg["aquifer"].get("min_saturated_thickness_m", 1.0))

    policy = str(cfg["unassigned"]["policy"])
    if policy not in ("exclude", "impermeable"):
        raise SystemExit(f"unassigned.policy {policy!r} is not one this solver has")
    conductive = land & ~unassigned if policy == "exclude" else land.copy()
    if policy == "impermeable":
        k0 = np.where(unassigned, 0.0, k0)

    print(f"  gravity {gravity} m/s2 -> conductivity is "
          f"{gravity / 9.80665:.3f}x Earth's for the same permeability")
    print(f"  unassigned lithology on {int((unassigned & land).sum()):,} land "
          f"regions ({(unassigned & land).sum() / land.sum() * 100:.2f}%), "
          f"policy {policy!r}")
    live = k0[land][k0[land] > 0]
    if thickness_note["source"] == "constant":
        print(f"  aquifer thickness {thickness_m[0]:.0f} m, constant")
    else:
        print(f"  aquifer thickness from {thickness_note['source']}, floored at "
              f"{thickness_note['floor_m']:.0f} m: land median "
              f"{thickness_note['land_median_m']:.0f} m, 95th "
              f"{thickness_note['land_p95_m']:.0f} m, max "
              f"{thickness_note['land_max_m']:.0f} m, "
              f"{thickness_note['at_floor_land_fraction']:.1%} of land on the floor")
    print(f"  transmissivity {'UNCONFINED, T = K (h - z_bottom)' if unconfined else 'confined, T = K D'}"
          f"; conductivity spans {np.log10(live.max() / live.min()):.1f} orders "
          f"of magnitude across lithologies")

    # -- forcing -----------------------------------------------------------
    print(f"reading the climatology {clim_name}")
    recharge, fields = recharge_field(export, config, clim)
    lat, lon, runoff_grid, precip_grid, evap_grid, lsm, evap_annual_grid = fields
    recharge = np.where(land, recharge, 0.0)
    # None on every arm that does not run the sink, so the report says "no
    # interval was measured" rather than carrying a stale one from a branch
    # that never evaluated it.
    evaporation_interval = None

    surface_m = export.elevation_km.astype(np.float64) * 1000.0

    if args.reduction_test:
        print("\nREDUCTION IDENTITY: zero permeability must reproduce the "
              "surface-only balance exactly")
        k0 = np.zeros_like(k0)
        # THE TERMS THE SURFACE-ONLY BALANCE DOES NOT HAVE ARE PART OF THE
        # CASE, not a configuration the caller may vary. `surface_water.nc` is
        # solved with no groundwater sink and no imposed water-table heads, so a
        # reduction carrying either is not reducing to it: with the sink on a
        # cell returns its recharge less what the sink took, and with local
        # baselevels on the river cells leave the network entirely. The identity
        # went silent that way once already -- world-qq10 is the same class in
        # the uniqueness arm -- so the case sets them here and says so rather
        # than depending on the caller remembering two flags.
        if not (args.no_groundwater_et and args.no_baselevels):
            print("  GW-15's sink and GW-17's baselevels are OFF for this arm: "
                  "the balance\n  being reduced to has neither, so an identity "
                  "carrying one is not one")
        args.no_groundwater_et = True
        args.no_baselevels = True
    if args.divide_test:
        # IMPOSED, not solved. The declared test is what a terrain-following
        # table does, so the table is set to the terrain rather than solved for
        # and then compared to the terrain it did not quite reproduce. Uniform
        # permeability, so nothing but the topography can steer the flow, and
        # the answer is fixed in advance: the groundwater divides must be the
        # topographic ones.
        # ON THE FILLED SURFACE, which is the whole correction. `terminal` is
        # a priority flood on the filled surface, so a groundwater trace on the
        # raw one is describing different drainage and cannot agree with it.
        print("\nDIVIDE TEST: uniform permeability, table following the "
              "FILLED surface")
        with Dataset(data / "regions.nc") as ds:
            filled_m = np.asarray(ds["filled_km"][:]).astype(np.float64) * 1000.0
        k0 = np.where(land, 1e-4, 0.0)
        res = gw.terrain_following(export, geom, k0_m_s=k0,
                                   thickness_m=thickness_m,
                                   surface_m=filled_m, conductive=conductive)
    else:
        # GW-15's sink, off unless asked for. ET_max is the Penman field already
        # read above for the lakes, so one evaporation rule governs the lakes,
        # the carve verdict and the water table. `region_grid_cells` is the
        # export/ExoPlaSim join CLAUDE.md rule 3 governs, and it is reused here
        # rather than repeated.
        et_cfg = cfg.get("evapotranspiration", {})
        bl_cfg = cfg.get("baselevels", {})
        et_on = et_cfg.get("enabled", False) and not args.no_groundwater_et
        et_lambda = (args.et_lambda if args.et_lambda is not None
                     else float(et_cfg.get("lambda_m", 1.0)))
        et_max = None
        if et_on:
            row, col = sw.region_grid_cells(export, lat)
            et_max = np.where(land, np.clip(evap_grid[row, col], 0.0, None), 0.0)
            # THIS WORLD'S YEAR, from lib/orbit. The rate reported here was per
            # 365.25 days, which is Earth's, so the only number a reader could
            # compare against a published potential-ET was in the wrong units.
            mm_per_year = orbit.orbital_year_days(config) * 86400.0 * 1000.0
            print(f"groundwater ET on, lambda {et_lambda} m, ET_max from Penman: "
                  f"land median "
                  f"{np.median(et_max[land]) * mm_per_year:.0f} mm/yr")
            # WHAT THE INTERVAL IS WORTH, measured on the run that uses it. The
            # solve takes the bin mean; the annual evaluation is the number this
            # path used to take, and it is kept here as the size of the change
            # rather than as a second answer. In the local-sink regime the depth
            # is `lambda ln(ET_max A / supply)` exactly, so a ratio `f` between
            # the two arms is a depth shift of `lambda ln f` and needs no second
            # solve to state. `sink_fraction` below says where that regime holds.
            et_annual = np.where(
                land, np.clip(evap_annual_grid[row, col], 0.0, None), 0.0)
            # ON THE CELLS THE LOCAL BALANCE GIVES A POSITIVE DEPTH ON, under
            # BOTH arms. Where `ET_max` cannot take the cell's own recharge the
            # local balance puts the table at the surface, so a depth shift
            # quoted there is a shift on a depth of zero. The cells that CROSS
            # that line are the categorical part of the change and are counted
            # separately.
            _pos = land & (recharge > 0)
            _l = _pos & (et_annual > recharge) & (et_max > recharge)
            _f = et_max[_l] / et_annual[_l]
            _unpinned = int((_pos & (et_annual <= recharge)
                             & (et_max > recharge)).sum())
            evaporation_interval = {
                "used": "the bin mean of the per-bin Penman",
                "other_end": "one Penman evaluation on annual-mean air",
                "why_the_bin_mean_is_used": (
                    "the sink is linear in ET_max and nonlinear in the depth, "
                    "and a steady-state solve holds the depth fixed through the "
                    "year, so the annual mean of ET_max is the only input the "
                    "solver can take. The evaporating body is the ground, whose "
                    "seasonal heat storage is a ground heat flux of a few W/m2 "
                    "against a net radiation cycle an order larger, so the "
                    "no-storage end of cv._INTERVAL_BRACKET is the near one "
                    "rather than a bound"),
                "et_max_land_median_mm_per_year_bin_mean": float(
                    np.median(et_max[_l]) * mm_per_year),
                "et_max_land_median_mm_per_year_annual": float(
                    np.median(et_annual[_l]) * mm_per_year),
                "capacity_ratio_bin_over_annual": float(
                    et_max[_l].sum() / max(et_annual[_l].sum(), 1e-30)),
                "ratio_percentiles": {str(p): float(np.percentile(_f, p))
                                      for p in (5, 25, 50, 75, 95)},
                "land_fraction_bin_below_annual": float((_f < 1.0).mean()),
                "regions_scored": int(_l.sum()),
                "regions_unpinned_by_the_bin_mean": _unpinned,
                "regions_unpinned_share_of_recharging_land": float(
                    _unpinned / max(int(_pos.sum()), 1)),
                "unpinned_note": (
                    "the local balance cannot take the cell's recharge under "
                    "the annual evaluation and can under the bin mean, so the "
                    "table leaves the surface. That moves at_surface and the "
                    "seepage a pinned cell returns, which is a different kind "
                    "of change from a depth shift"),
                "sign_is_not_one_signed": (
                    "Penman is not convex in temperature alone: the "
                    "Delta/(Delta+gamma) energy weighting saturates and the "
                    "aerodynamic term is a product of wind with a deficit, so "
                    "the curvature of the composite changes sign. The bin mean "
                    "falls BELOW the annual evaluation on the warm cells with "
                    "the smallest seasonal swing"),
                "depth_shift_m_at_lambda": {
                    str(lam): {
                        "median": float(lam * np.log(np.median(_f))),
                        "p95": float(lam * np.log(np.percentile(_f, 95))),
                    } for lam in (0.5, float(et_lambda), 2.0)},
                "depth_shift_note": (
                    "d = lambda ln(ET_max A / supply) in the local-sink regime, "
                    "so the shift is lambda ln f exactly. Read it against the "
                    "lambda bracket in groundwater.yaml, which multiplies the "
                    "same depth by four"),
            }
            print(f"  interval: bin mean is "
                  f"{evaporation_interval['capacity_ratio_bin_over_annual']:.3f}x "
                  f"the annual evaluation over land, below it on "
                  f"{evaporation_interval['land_fraction_bin_below_annual']:.1%} "
                  f"of land;\n  worth "
                  f"{evaporation_interval['depth_shift_m_at_lambda'][str(float(et_lambda))]['median']:+.3f} m "
                  f"on the median local-sink depth and "
                  f"{evaporation_interval['depth_shift_m_at_lambda'][str(float(et_lambda))]['p95']:+.3f} m "
                  f"at its 95th")
        else:
            et_lambda = None
            print("groundwater ET OFF: the table will pin at the surface")

        # GW-17. A river or lake surface IS the water table there, so it is a
        # fixed head. `surface_water.nc` is the only thing that knows where they
        # are, which is why `groundwater` needs `surface_water` in the graph.
        fixed_head = None
        if bl_cfg.get("enabled", False) and not args.no_baselevels:
            sw_path = data / "surface_water.nc"
            if not sw_path.exists():
                raise SystemExit(
                    f"{rel(sw_path)} is missing and local baselevels are on. "
                    "Run `surface_water` first, or pass --no-baselevels and "
                    "accept the ocean as the only boundary.")
            with Dataset(sw_path) as ds:
                lake = np.asarray(ds["lake"][:]) > 0
                disch = np.asarray(ds["discharge_m3_s"][:])
            wet = land & ((disch >= float(bl_cfg.get("discharge_min_m3_s", 10.0)))
                          | (lake if bl_cfg.get("include_lakes", True) else False))
            fixed_head = np.full(export.n_regions, np.nan)
            fixed_head[wet] = surface_m[wet]
            print(f"local baselevels: {int(wet.sum()):,} river and lake cells "
                  f"({wet.sum() / max(int(land.sum()), 1):.1%} of land) at fixed head")
        else:
            print("local baselevels OFF: the ocean is the only fixed head")

        # GW-24. The aquifer BASE, where it is asked for, is the surface less
        # the same thickness the confined form multiplies by. It is NaN off
        # land, which is what tells the solver a face has no column on that
        # side and must take its neighbour's.
        base_m = np.where(land, surface_m - thickness_m, np.nan) if unconfined else None

        # ONE ARGUMENT LIST, BUILT ONCE. The uniqueness check below re-solves
        # the same problem from the opposite initial active set, and while the
        # two calls were written out separately the check quietly stopped being
        # an identity every time a term landed here and not there: GW-15's sink
        # and GW-17's baselevels were both missing from the re-solve for as long
        # as they existed. A second trajectory through a DIFFERENT model is not
        # a comparison, so the two calls now differ in exactly one keyword and
        # a new term cannot reach one without reaching the other. world-60x0.
        solve_kwargs = dict(
            k0_m_s=k0, thickness_m=thickness_m, recharge_m_s=recharge,
            surface_m=surface_m, conductive=conductive,
            max_outer=args.max_outer, et_max_m_s=et_max,
            et_lambda_m=et_lambda, fixed_head_m=fixed_head,
            aquifer_base_m=base_m, min_saturated_m=min_saturated)

        print("solving the water table")
        res = gw.solve(export, geom, **solve_kwargs)

    area_m2 = geom.volume_area_m2
    supply = recharge * area_m2
    # Closure is a statement about a solved water balance. The divide test
    # imposes a table rather than solving one, so there is no balance to close
    # and reporting a residual for it would invent a failure.
    clos = None
    if not args.divide_test:
        clos = gw.closure(res, recharge, area_m2, export)
        print(f"\nCLOSURE  recharge {clos['recharge_m3_s']:,.3f} m3/s   "
              f"discharge {clos['discharge_m3_s']:,.3f} m3/s")
        print(f"         relative residual {clos['relative_residual']:.3e} "
              f"against {clos['tolerance']:.0e}   "
              f"{'pass' if clos['passes'] else 'FAIL'}")
        print(f"         seepage {clos['seepage_m3_s']:,.1f}   "
              f"to the ocean underground {clos['to_ocean_m3_s']:,.1f}")

    depth = res["depth_m"]
    print(f"\nwater table depth on land: median {np.median(depth[land]):.2f} m, "
          f"at the surface on {(depth[land] <= 0.01).mean() * 100:.1f}% of land")
    drained = np.clip(depth[land] / np.maximum(thickness_m[land], 1e-9), 0.0, 1.0)
    print(f"  saturated column drained, median {np.median(drained):.3%}, 95th "
          f"{np.percentile(drained, 95):.3%}: what the confined transmissivity "
          f"is wrong by,\n  and therefore what GW-24's unconfined form is worth "
          f"here. `--unconfined` is decided on this.")

    with Dataset(data / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:])
    basins = lb.BasinSet(data / "basins.nc")
    qg, gain, give = basin_exchange(terminal, res["seepage_m3_s"], supply,
                                    basins.n, res.get("et_m3_s"))

    report = {
        "source_build": build.name,
        "terrain_hash": export.terrain_hash,
        "forcing": clim_name,
        # WHICH STAGE the recharge came from. lib/paths.py.
        "climatology_stage": clim_stage,
        "forcing_sha256": sha256(clim),
        "config_sha256": sha256(CONFIG),
        "sigma": args.sigma,
        "operator_noise": args.operator_noise,
        "noise_seed": args.noise_seed,
        "unassigned_policy": policy,
        # WHICH INTERVAL THE SINK'S PENMAN WAS EVALUATED OVER, and what the
        # other end would have been. None when the sink is off.
        "evaporation_interval": evaporation_interval,
        "gravity_m_s2": gravity,
        "conductivity_vs_earth": gravity / 9.80665,
        "geometry": {
            "faces": int(geom.src.size),
            "voronoi_area_over_sphere": closes,
            "export_cell_area_over_sphere": float(
                geom.volume_area_m2.sum() / (4 * np.pi * geom.radius_m ** 2)),
            "note": ("flux uses the Voronoi area its own faces bound; water "
                     "volumes use the export's cell_area so the reduction "
                     "identity against surface_water.py is exact"),
        },
        "laplace_beltrami": {
            "criterion": gw.LAPLACE_TOLERANCE,
            "relative_rms": {str(k): v for k, v in
                             gw.laplace_beltrami_error(geom).items()},
        },
        "closure": clos,
        "converged": bool(res.get("converged", True)),
        "convergence": {
            "criterion_residual": gw.RESIDUAL_TOLERANCE,
            "criterion_flip_fraction": gw.FLIP_TOLERANCE,
            "final_residual": res.get("final_residual"),
            "criterion_note": ("the water balance residual on free cells over "
                               "total land recharge, and a still active set. "
                               "Not the head step: a bar on the step can be "
                               "met by a shrinking step rather than by "
                               "convergence"),
            "trace_columns": ["pass", "free_cells", "cells_flipped",
                              "residual", "infeasible_pinned_cells"],
            "seepage_clipped_m3_s": res.get("seepage_clipped_m3_s"),
            "trace": res.get("residual_trace", []),
        },
        "outer_iterations": res.get("outer_iterations"),
        "water_table": {
            "median_depth_m_land": float(np.median(depth[land])),
            "at_surface_fraction_land": float((depth[land] <= 0.01).mean()),
            "unassigned_land_fraction": float(
                (unassigned & land).sum() / max(int(land.sum()), 1)),
            "aquifer_thickness": thickness_note,
            # GW-24, MEASURED FROM THE CONFINED RUN ITSELF. Under `T = K D` the
            # transmissivity is wrong by exactly the fraction of the saturated
            # column the water table has drained, `depth / D`, because that is
            # the column the unconfined form would have removed. So the confined
            # solve reports how much the term it omits is worth without anyone
            # having to solve for it, and the decision to enable `--unconfined`
            # rests on this number rather than on the argument. Under the
            # unconfined form it is the correction ALREADY applied.
            "saturated_column_drained": {
                str(p): float(np.percentile(
                    np.clip(depth[land] / np.maximum(thickness_m[land], 1e-9),
                            0.0, 1.0), p))
                for p in (50, 75, 95, 99)},
            "unconfined": bool(unconfined),
            "at_transmissivity_floor": int(res["at_transmissivity_floor"].sum()),
            "transmissivity_note": (
                "T = K (h - z_bottom): the saturated column is solved for, so a "
                "cell at the minimum saturated thickness carries a LOWER BOUND "
                "on its depth rather than a value, and at_transmissivity_floor "
                "counts them"
                if unconfined else
                "T = K D: the transmissivity does not depend on the head, so "
                "there is no depth floor and no cell has a transmissivity that "
                "underflows"),
        },
        "git_commit": subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
    }

    # -- the declared tests ------------------------------------------------
    if args.reduction_test:
        seep = res["seepage_m3_s"]
        exact = np.array_equal(seep[land], supply[land])
        worst = float(np.abs(seep[land] - supply[land]).max())
        qg_max = float(np.abs(qg).max())
        print(f"\n  seepage equals recharge on every land cell: {exact}")
        print(f"  worst absolute difference {worst:.6e} m3/s")
        print(f"  largest |Qg| over {basins.n:,} basins {qg_max:.6e} m3/s")
        report["reduction_identity"] = {
            "criterion": "exact: seepage == recharge per cell, Qg == 0 per basin",
            "seepage_equals_recharge_exactly": bool(exact),
            "worst_absolute_difference_m3_s": worst,
            "largest_abs_basin_qg_m3_s": qg_max,
            "passes": bool(exact and qg_max == 0.0),
        }
        print("  " + ("PASS" if report["reduction_identity"]["passes"] else "FAIL"))

    if args.divide_test:
        rec = gw.groundwater_receiver(export.n_regions, res["src"], res["dst"],
                                      res["face_flux_m3_s"], land)
        gwterm = gw.trace_terminals(rec, terminal, land)
        flat = land & (filled_m > surface_m + 1e-6)
        unfilled = land & ~flat
        same = gwterm == terminal
        a_all = float(area_m2[land][same[land]].sum() / area_m2[land].sum())
        a_unf = float(area_m2[unfilled][same[unfilled]].sum()
                      / area_m2[unfilled].sum())
        flat_share = float(area_m2[flat].sum() / area_m2[land].sum())
        print(f"\n  outside the pits the flood filled: {a_unf * 100:.2f}% of "
              f"land area agrees (criterion "
              f"{gw.DIVIDE_AGREEMENT_UNFILLED * 100:.0f}%)")
        print(f"  those pits are {flat_share * 100:.2f}% of land and carry no "
              f"gradient to trace; over all land {a_all * 100:.2f}%")
        report["divide_test"] = {
            "criterion_unfilled": gw.DIVIDE_AGREEMENT_UNFILLED,
            "surface": "filled_km, the same surface regions.nc:terminal uses",
            "agreement_outside_filled_pits": a_unf,
            "agreement_all_land_area": a_all,
            "filled_pit_share_of_land": flat_share,
            "passes": bool(a_unf >= gw.DIVIDE_AGREEMENT_UNFILLED),
        }
        print("  " + ("PASS" if report["divide_test"]["passes"] else "MISS"))

        # THE TRACE IDENTITY. Synthesise a flux field whose only outgoing flux
        # at each land cell is the face to that cell's own surface `receiver`.
        # The trace must then reproduce `terminal` exactly, because it is being
        # handed the surface network itself. Nothing about groundwater enters.
        with Dataset(data / "regions.nc") as ds:
            surf_recv = np.asarray(ds["receiver"][:]).astype(np.int64)
        ssrc, sdst = res["src"], res["dst"]
        synth = np.zeros(ssrc.size)
        # positive flux is dst into src, so an outgoing flux from src is -1
        synth[surf_recv[ssrc] == sdst] = -1.0
        synth[surf_recv[sdst] == ssrc] = +1.0
        srec = gw.groundwater_receiver(export.n_regions, ssrc, sdst, synth, land)
        sterm = gw.trace_terminals(srec, terminal, land)
        agree = float(area_m2[land][(sterm == terminal)[land]].sum()
                      / area_m2[land].sum())
        rec_agree = float((srec[land] == surf_recv[land]).mean())
        ok = bool(agree >= gw.TRACE_IDENTITY_EXACT)
        report["trace_identity"] = {
            "criterion": gw.TRACE_IDENTITY_EXACT,
            "what": ("the trace machinery handed the surface network's own "
                     "receiver as a flux field; nothing about groundwater "
                     "enters, so it must reproduce terminal exactly"),
            "receiver_recovered_fraction": rec_agree,
            "terminal_agreement_land_area": agree,
            "passes": ok,
        }
        print(f"\n  TRACE IDENTITY: receiver recovered on {rec_agree * 100:.4f}% "
              f"of land regions, terminal on {agree * 100:.4f}% of land area")
        print("  " + ("PASS" if ok else "MISS"))

    # -- GW-4: what the term would move ------------------------------------
    if not (args.reduction_test or args.divide_test):
        report["basin_exchange"] = measure_carve_effect(
            export, data, basins, terminal, qg, give, gain, config, clim,
            fields, area_m2)

    # -- products ----------------------------------------------------------
    # UNIQUENESS, which replaces the cross-scheme check the nonlinear model
    # needed. With transmissivity independent of the head there is no second
    # algebra to compare against -- Picard and Kirchhoff both collapse onto this
    # same linear system. What linearity buys instead is stronger: the matrix is
    # symmetric positive definite, so the complementarity problem has exactly
    # ONE solution, and any two active-set trajectories must land on the
    # identical head field. Starting from every cell free rather than every cell
    # pinned walks a completely different path to it. That is an identity, it
    # can genuinely fail, and it is what would catch a bug in the release and
    # pin logic that a single run cannot see.
    #
    # THE SECOND TRAJECTORY MUST SOLVE THE SAME PROBLEM. Every argument the
    # sink and the baselevels are given for is a term in the equation, so a
    # re-solve that omits one is comparing two models rather than two paths to
    # the same solution, and cannot fail in the direction that matters. That is
    # why `solve_kwargs` is assembled once above and both calls take it whole:
    # the two differ in `start_all_free` and in nothing else, by construction.
    #
    # AND UNDER THE UNCONFINED FORM IT IS NOT AN IDENTITY. GW-24's
    # transmissivity depends on the head, so the matrix is not fixed and the
    # symmetric-positive-definite argument above does not apply. The bar moves
    # from bit-comparison to `UNCONFINED_HEAD_RELATIVE`, declared in
    # `groundwater.py` before any unconfined run.
    if (args.uniqueness_check and not args.divide_test
            and res.get("converged", True) is not False):
        print("\nUNIQUENESS: re-solving from the opposite initial active set")
        alt = gw.solve(export, geom, **solve_kwargs, start_all_free=True)
        d = np.abs(alt["head_m"] - res["head_m"])[conductive]
        scale = max(float(np.abs(res["head_m"][conductive]).max()), 1.0)
        rel_head = float(d.max() / scale)
        bar_head = (gw.UNCONFINED_HEAD_RELATIVE if unconfined
                    else gw.SCHEME_HEAD_RELATIVE)
        ok = bool(alt.get("converged") and rel_head < bar_head)
        report["uniqueness"] = {
            "criterion_relative_head": bar_head,
            "is_identity": not unconfined,
            "max_absolute_head_difference_m": float(d.max()),
            "relative_head_difference": rel_head,
            "alternate_converged": bool(alt.get("converged")),
            "alternate_passes": alt.get("outer_iterations"),
            "passes": ok,
        }
        print(f"  max |h_altstart - h| {d.max():.3e} m, relative {rel_head:.3e} "
              f"against {bar_head:.0e}"
              + ("" if not unconfined else
                 " (a MEASUREMENT, not an identity: the unconfined "
                 "transmissivity depends on the head)"))
        print("  " + ("PASS" if ok else "MISS"))

    if not res.get("converged", True):
        # The report goes out; the water table does not. An unconverged head
        # field written to `data/` would be read downstream as a result, and
        # nothing about the file would say it was not one.
        ANALYSIS.mkdir(parents=True, exist_ok=True)
        rp = report_path(args, artifact_path(args, data))
        rp.write_text(json.dumps(report, indent=2) + "\n")
        print(f"\nwrote {rp}")
        print("NOT writing water_table.nc: the solve did not converge, and an "
              "unconverged head field is not a result.")
        return 1

    # WHICH TERM SET THE DEPTH, per cell. GW-22 scored this model on two
    # continents and found it has a regime: where the sink takes the recharge
    # locally the depth is a recharge map (Australia, model spread 1.91 m
    # against an observed 18.36), and where lateral flow carries it the solve is
    # doing physics (United States, 34.02 m against 34.26, and Pearson +0.26 on
    # bores the USGS labels unconfined). The two are separable here.
    et_cell = res.get("et_m3_s")
    with np.errstate(divide="ignore", invalid="ignore"):
        sink_fraction = np.where(
            supply > 0,
            (et_cell if et_cell is not None else np.zeros_like(supply)) / supply,
            np.nan).astype(np.float32)
    sink_fraction = np.clip(np.nan_to_num(sink_fraction, nan=np.nan), 0.0, 1.0)
    lp = np.isfinite(sink_fraction) & land
    if lp.any():
        print(f"\n  sink fraction on land: median {np.nanmedian(sink_fraction[lp]):.3f}, "
              f"above 0.9 on {np.nansum(sink_fraction[lp] > 0.9) / lp.sum():.1%} of land")
        print("  (near 1 the depth is a recharge map; GW-22 measured that regime)")

    out = artifact_path(args, data)
    if out.exists():
        raise SystemExit(
            f"{out} exists. This step writes a new artifact and never "
            "overwrites one; move it aside or pass --output.")
    with Dataset(out, "w") as ds:
        ds.createDimension("region", export.n_regions)
        ds.createDimension("basin", basins.n)
        ds.title = f"Steady-state water table under {clim.stem}"
        ds.terrain_hash = export.terrain_hash
        ds.setncattr("vesper_source_build", build.name)
        ds.forcing = clim_name
        ds.climatology_stage = clim_stage
        ds.sigma = args.sigma
        ds.caveat = (
            "An equilibrium water table, not an aquifer with a history. "
            "Permeability carries 1.5 to 2.5 orders of magnitude of Gleeson "
            "spread, so depth is a bracket; run --sigma -1 and +1 for it. "
            "Valley-to-ridge convergence is sub-grid at this mesh and absent. "
            "READ sink_fraction BEFORE USING depth_m: where it approaches 1 the "
            "evapotranspiration sink took the recharge locally and the depth is "
            "a function of recharge rather than a flow solution. GW-22 scored "
            "both regimes against real bores."
        )
        for name, dat, dtype, dim, units, note in [
            ("head_m", res["head_m"], "f4", "region", "m",
             "water table elevation above sea level"),
            ("depth_m", depth, "f4", "region", "m",
             "water table below the land surface"),
            ("seepage_m3_s", res["seepage_m3_s"], "f4", "region", "m3 s-1",
             "groundwater returning to the surface"),
            ("transmissivity_m2_s", res["transmissivity_m2_s"], "f4", "region",
             "m2 s-1",
             ("T = K (h - z_bottom), the unconfined saturated column, GW-24. "
              "NOT Fan's exponential decay: exp(h/f) is convex, so a cell mean "
              "carries exp(sigma^2/2f^2) and has no value at this spacing, "
              "while this form is LINEAR in the head and so has an exact one"
              if unconfined else
              "T = K D at the aquifer thickness. NOT Fan's exponential decay, "
              "which GW-9 removed: exp(h/f) is convex, so a cell mean carries "
              "exp(sigma^2/2f^2) and has no value at this spacing")),
            ("at_surface", res["pinned"].astype(np.int8), "i1", "region", "1",
             "water table pinned at the land surface"),
            ("sink_fraction", sink_fraction, "f4", "region", "1",
             "share of a cell's recharge removed by groundwater ET rather than "
             "carried laterally. Near 1 the depth is the local balance "
             "lambda ln(et_max A / supply), a recharge map in a water table's "
             "units. Below 1 the sink did NOT set the depth, which is weaker "
             "than the flow solve setting it: a cell pinned at the surface sheds "
             "its recharge as seepage and also scores low, so read at_surface "
             "alongside. GW-22 measured the regimes against real bores"),
            ("excluded", res["excluded"].astype(np.int8), "i1", "region", "1",
             "lithology has no assigned permeability"),
            ("basin_groundwater_m3_s", qg, "f8", "basin", "m3 s-1",
             "net groundwater import; Fan (2019) Qg, positive into the basin"),
            ("basin_seepage_m3_s", gain, "f8", "basin", "m3 s-1",
             "seepage within the catchment"),
            ("basin_recharge_m3_s", give, "f8", "basin", "m3 s-1",
             "recharge falling on the catchment"),
        ]:
            v = ds.createVariable(name, dtype, (dim,), zlib=True)
            v.units = units
            if note:
                v.long_name = note
            v[:] = dat

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    rp = report_path(args, out)
    rp.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}\nwrote {rp}")
    return 0


def measure_carve_effect(export, data, basins, terminal, qg, give, gain,
                         config, clim, fields, area_m2) -> dict:
    """GW-4. How many carve verdicts would the `Qg` term move?

    A MEASUREMENT. The criterion is `carve_verdict.py`'s and is not restated
    here: a basin overflows when its aridity index falls at or below the
    geometric `critical_aridity_index` that `basins.nc` already carries. What
    changes is the denominator. The surface balance divides the lake's
    evaporative deficit by the catchment's own recharge depth; with groundwater
    the basin is fed by what actually seeps inside its catchment, which is that
    recharge plus whatever crossed the divide underground.

    Nothing here is written back into `carve_verdict.py`. Whether this term
    belongs in the criterion is a loop A decision and it is taken by reading
    this, not by running it.
    """
    # The BIN MEAN of the per-bin Penman, which is what `carve_verdict.py`
    # decides its primary bound on. This function used to read the annual
    # evaluation, so it measured the groundwater term against the OTHER end of
    # `cv._INTERVAL_BRACKET` from the verdict it claims to be moving. The annual
    # arm is unpacked and unused here: this reports a count of flips, and
    # reporting it twice would be two numbers for one question.
    lat, lon, runoff_grid, precip_grid, evap_grid, lsm, _ = fields
    sinks = np.array([b.sink for b in export.basins])
    # Only the lake fluxes are taken from here. The catchment runoff depth is
    # recomputed on the mesh below, for BOTH arms, because the comparison has to
    # isolate the groundwater term: `per_basin_forcing` bins the catchment over
    # whole grid cells, and pairing that binned depth against a mesh-resolved
    # one would fold hydrography's known binning approximation into a number
    # that is supposed to measure groundwater and nothing else.
    _, lake_precip, lake_evap = sw.per_basin_forcing(
        basins.n, lat, lon, runoff_grid, precip_grid, evap_grid, sinks,
        export, lsm, config)

    year_s = orbit.orbital_year_days(config) * sw.SECONDS_PER_DAY
    to_km_yr = year_s / 1000.0

    catch_m2 = np.zeros(basins.n)
    sel = terminal >= 0
    np.add.at(catch_m2, terminal[sel], area_m2[sel])
    ok = catch_m2 > 0

    # Depths, m/s, over the catchment. The baseline is the catchment's own
    # recharge; the groundwater number adds the EXCHANGE and nothing else.
    #
    # It used to be `gain`, the seepage. Without a sink that is the same number,
    # exactly: `qg = gain - give`, so `give + qg == gain` identically, and this
    # change is inert on every run made before GW-15. WITH a sink they part
    # company, because seepage is recharge minus what evaporated plus what was
    # exchanged, and attributing that difference to groundwater EXCHANGE credits
    # the sink's evaporation to a mechanism that did not do it. The question
    # this function asks is what the exchange does to the verdict, so the
    # exchange is what varies and everything else is held.
    #
    # It also matches `carve_verdict.py`'s own supply term under GW-14,
    # `r_eff = (surface_runoff * C + Qg) / C`, which is the point: one quantity
    # measured here and applied there rather than two that drift.
    surface_depth = np.where(ok, give / np.maximum(catch_m2, 1e-30), 0.0)
    seeped_depth = np.where(ok, (give + qg) / np.maximum(catch_m2, 1e-30), 0.0)

    with Dataset(data / "basins.nc") as ds:
        crit = np.asarray(ds["critical_aridity_index"][:])

    deficit = (lake_evap - lake_precip) * to_km_yr        # km/yr
    def verdict(depth_m_s):
        r = depth_m_s * to_km_yr
        with np.errstate(divide="ignore", invalid="ignore"):
            idx = np.where(r > 0, deficit / np.maximum(r, 1e-30), np.inf)
        return idx, (idx <= crit)

    idx_s, carve_s = verdict(surface_depth)
    idx_g, carve_g = verdict(seeped_depth)
    scored = ok & np.isfinite(crit)
    flipped = scored & (carve_s != carve_g)

    net = qg / np.maximum(give, 1e-30)
    big = scored & (np.abs(net) > 0.10)
    return {
        "what_this_is": ("a measurement, not a change: Qg is reported per basin "
                         "and is NOT wired into carve_verdict.py"),
        "both_arms_use_mesh_resolved_catchment_depth": True,
        "why": ("the surface arm is recomputed on the mesh rather than taken "
                "from the coupling matrix, so the only difference between the "
                "two arms is the groundwater term"),
        "basins_scored": int(scored.sum()),
        "total_abs_exchange_m3_s": float(np.abs(qg).sum()),
        "total_recharge_m3_s": float(give.sum()),
        "exchange_over_recharge": float(np.abs(qg).sum() / max(give.sum(), 1e-30)),
        "basins_gaining": int((scored & (qg > 0)).sum()),
        "basins_losing": int((scored & (qg < 0)).sum()),
        "basins_shifted_over_10_percent": int(big.sum()),
        "median_abs_fractional_shift": float(np.median(np.abs(net[scored]))),
        "carve_verdicts": {
            "criterion": ("carve_verdict.py's own: overflow when the aridity "
                          "index is at or below basins.nc's "
                          "critical_aridity_index"),
            "would_carve_surface_only": int(carve_s[scored].sum()),
            "would_carve_with_groundwater": int(carve_g[scored].sum()),
            "flipped": int(flipped.sum()),
            "flipped_to_carve": int((flipped & carve_g).sum()),
            "flipped_to_hold": int((flipped & ~carve_g).sum()),
        },
    }


if __name__ == "__main__":
    raise SystemExit(main())
