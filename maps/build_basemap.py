"""Render a natural-colour equirectangular base map of Vesper.

Terrain comes from the native mesh of the active build (nearest region per
pixel, so coastlines and closed-basin floors stay sharp). Colour
comes from the baseline ExoPlaSim climatology: land is tinted by its biome
class, ocean by depth, ice and snow by the warmest-month surface temperature
corrected from the climatology's own orography to the high-resolution one at
the environmental lapse rate `lib/lapse.py` measures from that same file.

The colour is illustrative. The biome field is on the climate grid of whichever
climatology is read, which is coarser than the mesh by orders of magnitude, and
was computed on the pre-carve terrain; so it is a plausible tint over the
geography being drawn rather than a result about it. The provenance sidecar
names the grid, read from the climatology rather than assumed.

Lakes and rivers are the exception: those come from the water balance solved in
hydrography/, on this terrain, and are a result rather than a tint.

Writes maps/build/basemap.npy (H x W x 3 uint8) and a provenance sidecar.
"""

import json
import subprocess
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# lib/ itself as well as the package: `gridding` imports `orogen` as a
# top-level module, the way every other component's `_paths.py` sets it up.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from lib import builds  # noqa: E402
import gridding  # noqa: E402
import lapse  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import frames  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = Path(__file__).resolve().parent / "build"
# Set by `--climatology`; None means take the configured one.
_CLIMATOLOGY_OVERRIDE: Path | None = None


# Resolved at call time, not import time, and from the config rather than a
# hardcoded directory.
def _climatology() -> Path:
    """The climatology FILE. It used to return the directory and callers
    appended `baseline_regular_climatology.nc`, which finds nothing once a
    product is labelled anything else -- a bootstrap climatology is
    `bootstrap_regular_climatology.nc`.

    `--climatology` overrides it, which is the only way to draw a map off a
    climatology `config/planet.yaml` deliberately does not name. Early in a
    cycle the one that exists is the bootstrap's, and the config key is held at
    null so that nothing reads it by DEFAULT; naming it per invocation is a
    caller's declaration rather than a fallback, and the tint is illustrative
    either way. `paths.climatology_path` names this flag in its own error.
    """
    if _CLIMATOLOGY_OVERRIDE is not None:
        return _CLIMATOLOGY_OVERRIDE
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "lib"))
    from paths import climatology_path
    return climatology_path()


def _classification() -> Path:
    """The classification product beside the configured climatology.

    Derived from the climatology's own filename rather than hardcoded, so it
    follows whatever label that product carries. Three call sites used to divide
    the climatology DIRECTORY by a fixed `baseline_classification.nc`, which
    broke the moment `_climatology()` started returning a file.
    """
    clim = _climatology()
    return clim.with_name(clim.name.replace("_regular_climatology.nc",
                                            "_classification.nc"))


def _grid_label() -> str:
    """The climate grid of the product actually read, as `LONxLAT`.

    Read from the file rather than named, for the reason the label above is:
    this string is emitted into `basemap_provenance.json` and inherited by every
    frame, so a hardcoded rung is a wrong number written INTO an artifact rather
    than merely into a comment. The rung the maps are drawn from moves with the
    escalation route, and nothing here would have learned it.
    """
    with nc.Dataset(_climatology()) as ds:
        return f"{ds.dimensions['lon'].size}x{ds.dimensions['lat'].size}"


def _caveat() -> str:
    """What a reader must know about the tint, taken from the product read."""
    label = _climatology().name.replace("_regular_climatology.nc", "")
    caveat = (f"biome and temperature fields are on the {_grid_label()} climate "
              f"grid and come from the "
              f"`{label}` climatology; the tint is illustrative, not a result")
    if label.startswith("bootstrap"):
        caveat += (". A bootstrap climatology is the FIRST run on the build, on "
                   "terrain-only surface fields, so its biomes are not the "
                   "baseline's and no number here is a baseline number")
    return caveat


# Render grid. `--width` moves it, because the sampling that is enough for one
# mesh is not enough for the next: 5760x2880 is about six pixels per region over
# a 2.5M-region export and under two over the ~10M-region one. Height follows
# width, since the grid is equirectangular. The nearest-region lookup is cached
# per size, so a change here costs one re-query and not one per map.
WIDTH, HEIGHT = 5760, 2880

# Biome palette, indexed by biome_index from baseline_classification.nc.
# 0 is "not land" and never drawn.
BIOME_RGB = np.array(
    [
        (0, 0, 0),  # 0  (unused)
        (26, 84, 44),  # 1  everwet tropical forest
        (44, 112, 54),  # 2  monsoon tropical forest
        (139, 160, 74),  # 3  seasonal tropical woodland / savanna
        (214, 176, 111),  # 4  desert
        (188, 166, 113),  # 5  steppe / semidesert
        (150, 160, 89),  # 6  mediterranean woodland / shrubland
        (74, 133, 63),  # 7  temperate / subtropical forest
        (94, 141, 100),  # 8  cool oceanic forest / heath
        (61, 118, 69),  # 9  continental temperate / mixed forest
        (52, 94, 70),  # 10 boreal forest / woodland
        (150, 154, 133),  # 11 tundra / alpine
    ],
    dtype=np.float64,
)

# Ocean depth ramp, (depth_km, rgb). Interpolated on sqrt(depth).
OCEAN_RAMP = [
    (0.0, (134, 190, 216)),
    (0.3, (98, 162, 200)),
    (1.5, (60, 124, 176)),
    (4.0, (38, 90, 140)),
    (9.0, (20, 55, 96)),
]

# Inland water, kept a shade greener and lighter than the sea so a lake reads as
# a lake rather than as a bay that lost its connection.
LAKE_RAMP = [
    (0.0, (108, 172, 178)),
    (0.05, (72, 140, 158)),
    (0.4, (44, 102, 126)),
    (2.0, (28, 72, 96)),
]
RIVER_RGB = np.array([70, 136, 158.0])
RIVER_MIN_M3_S = 120.0     # below this nothing is drawn
RIVER_FULL_M3_S = 8000.0   # at this the line is at full weight

ROCK_EVAPORITE, ROCK_PLAYA = 18, 19
SALT_RGB = np.array([236, 232, 222.0])
PLAYA_RGB = np.array([203, 187, 155.0])
ROCK_RGB = np.array([138, 126, 112.0])
ICE_RGB = np.array([238, 242, 246.0])
SEAICE_RGB = np.array([226, 234, 240.0])


def pixel_directions():
    """Unit vectors for the centre of every pixel of the equirectangular grid."""
    lon = (np.arange(WIDTH) + 0.5) / WIDTH * 360.0 - 180.0
    lat = 90.0 - (np.arange(HEIGHT) + 0.5) / HEIGHT * 180.0
    lo = np.deg2rad(lon)[None, :]
    la = np.deg2rad(lat)[:, None]
    coslat = np.cos(la)
    # Export convention: lat = asin(y), lon = atan2(x, z).
    x = coslat * np.sin(lo) * np.ones_like(la)
    y = np.sin(la) * np.ones_like(lo)
    z = coslat * np.cos(lo) * np.ones_like(la)
    return lat, lon, np.stack([x, y, z], axis=-1)


def region_index(src):
    """Nearest mesh region for each pixel. Cached; the query is the slow part.

    THE CACHE KEY IS THE MESH AND NOT ONLY THE RESOLUTION. It was
    `region_index_{WIDTH}x{HEIGHT}.npy` and returned unconditionally when the
    file existed, so a build change silently reused another build's answer: the
    first map of `canonical-10m-carve1` was drawn through an index built six
    days earlier on `precarve-craton`, addressing 2,500,001 regions of a
    10,000,005-region mesh. Every pixel read the wrong region. It did not look
    like noise -- both meshes cover the same sphere from the same seed, so the
    result was a plausible world with a radial smear at the antimeridian, which
    is the kind of wrong that gets shipped.

    The region COUNT keys the file, because that is what makes two of these
    incompatible, and the terrain hash goes in a sidecar and is refused on
    mismatch: two builds at the same region count are different meshes if
    anything upstream of the mesh moved. Belt and braces on purpose -- this is
    an artifact whose wrongness is invisible in its own output.
    """
    manifest = json.loads((src / "manifest.json").read_text(encoding="utf-8"))
    regions = int(manifest["numRegions"])
    terrain = str(manifest.get("hashes", {}).get("finalElevation", ""))
    cache = BUILD_DIR / f"region_index_{WIDTH}x{HEIGHT}_n{regions}.npy"
    stamp = cache.with_suffix(".json")
    if cache.exists() and stamp.exists():
        recorded = json.loads(stamp.read_text(encoding="utf-8"))
        if (recorded.get("numRegions") == regions
                and recorded.get("finalElevation") == terrain):
            idx = np.load(cache)
            if int(idx.max()) < regions:
                return idx
        print(f"  cached index is for another mesh; rebuilding")
    xyz = np.stack(
        [np.fromfile(src / f"raw/{c}.bin", dtype=np.float32) for c in "xyz"], axis=1
    ).astype(np.float64)
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)
    print(f"  building tree over {len(xyz):,} regions")
    tree = cKDTree(xyz)
    _, _, dirs = pixel_directions()
    print(f"  querying {WIDTH * HEIGHT:,} pixels")
    _, idx = tree.query(dirs.reshape(-1, 3), workers=-1)
    idx = idx.reshape(HEIGHT, WIDTH).astype(np.int32)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    np.save(cache, idx)
    stamp.write_text(json.dumps(
        {"numRegions": regions, "finalElevation": terrain,
         "width": WIDTH, "height": HEIGHT}, indent=2) + "\n", encoding="utf-8")
    return idx


def raw(src, name, dtype):
    return np.fromfile(src / f"raw/{name}.bin", dtype=dtype)


def upsample(field, smooth=True):
    """A climate field (lat x lon, south-to-north) onto the render grid.

    Bilinear, then a Gaussian of about half a climate cell. Without the blur the
    later clipped blends (ice, sea ice) plateau along straight bilinear
    contours and the map grows visible 128x64 polygons.
    """
    from scipy.ndimage import gaussian_filter, map_coordinates

    ds = nc.Dataset(_classification())
    clat = np.asarray(ds["lat"][:], dtype=float)
    clon = np.asarray(ds["lon"][:], dtype=float)
    ds.close()

    field = np.asarray(field, dtype=float)
    if clat[0] < clat[-1]:  # want north-first to match the render grid
        field = field[::-1]
        clat = clat[::-1]

    lat, lon, _ = pixel_directions()
    # Wrap longitude by padding a column at each end.
    padded = np.concatenate([field[:, -1:], field, field[:, :1]], axis=1)

    # Columns come from `lib/gridding.py`, which owns the convention, and the
    # climatology's own `lon` variable is deliberately not consulted for them.
    # This interpolated against those LABELS, which run 0..360 while the render
    # grid runs -180..180, and so drew every climate-derived layer on the
    # basemap 180 degrees from the terrain beneath it: the biome tint, the ice,
    # the sea ice. The mapping is the identity index for index -- the model's
    # own land mask reads back 1.0000 against the mesh that way and 0.5955
    # shifted by half the grid. See notes/audits/grid-convention-and-runoff.md.
    # Latitude IS interpolated by value, because the render rows are uniform and
    # the model's are Gaussian, so there is a genuine resample to do there.
    row = np.interp(lat, clat[::-1], np.arange(len(clat))[::-1])
    col = gridding.column_fraction(lon, len(clon)) + 1.0
    rr, cc = np.meshgrid(row, col, indexing="ij")
    out = map_coordinates(padded, [rr, cc], order=1, mode="nearest")
    if smooth:
        sigma = 0.5 * HEIGHT / len(clat)
        out = gaussian_filter(out, sigma=(sigma, sigma), mode=("nearest", "wrap"))
    return out


def polar_smooth(field, lat_deg):
    """Smooth a climate field along longitude, harder towards the poles.

    128 cells run round a latitude circle whatever its size, so near the pole
    they hold detail finer than the model can carry. Reprojected onto anything
    that shows the pole as a point, that detail becomes rays. Filtering with a
    width that grows as 1/cos(lat) removes it before it is drawn.
    """
    from scipy.ndimage import gaussian_filter1d

    out = np.array(field, dtype=float)
    for j, lat in enumerate(lat_deg):
        sigma = min(0.6 * (1.0 / max(np.cos(np.deg2rad(lat)), 1e-3) - 1.0), 40.0)
        if sigma > 0.3:
            out[j] = gaussian_filter1d(out[j], sigma, mode="wrap")
    return out


def fill_ocean_gaps(field, land):
    """Extend land values across the sea by nearest land cell.

    Coastal pixels of the render grid sit inside climate cells that are mostly
    ocean, so interpolating a land field against its raw zeros over water puts
    an ice sheet on every coast. Fill first, interpolate second.
    """
    from scipy.ndimage import distance_transform_edt

    field = np.asarray(field, dtype=float)
    tiled = np.tile(~np.asarray(land, dtype=bool), (1, 3))
    _, (rr, cc) = distance_transform_edt(tiled, return_indices=True)
    n = field.shape[1]
    out = np.tile(field, (1, 3))[rr, cc][:, n : 2 * n]
    return out


def surface_water_path() -> Path:
    """The lake and river solution for the build being drawn.

    Per-build. This read was the flat hydrography/data/, which holds whatever
    build was current when it was last written -- so a map could be drawn with
    one terrain's coastlines and another terrain's lakes, and would say nothing
    about it. Named here rather than inline because the provenance block has to
    hash the same file the compositing reads.
    """
    return builds.component_data("hydrography") / "surface_water.nc"


def load_surface_water(idx):
    """Lakes and rivers from the hydrography water balance, per pixel.

    Optional on purpose: the maps are built from `source/` and a climatology,
    and this is the one input that comes from another component's solver. If it
    has not been run, the map is drawn without standing water rather than
    failing, and says so.
    """
    path = surface_water_path()
    if not path.exists():
        print(f"  no surface_water.nc for this build; drawing without lakes or rivers")
        return None
    ds = nc.Dataset(path)
    lake = np.asarray(ds["lake"][:]).astype(bool)[idx]
    depth = np.asarray(ds["lake_depth_km"][:])[idx]
    discharge = np.asarray(ds["discharge_m3_s"][:])[idx]
    ds.close()
    return lake, depth, discharge


def lake_colour(depth_km):
    d = np.sqrt(np.clip(depth_km, 0, None))
    xs = np.array([np.sqrt(p[0]) for p in LAKE_RAMP])
    out = np.zeros(depth_km.shape + (3,))
    for c in range(3):
        out[..., c] = np.interp(d, xs, [p[1][c] for p in LAKE_RAMP])
    return out


def ocean_colour(depth_km):
    d = np.sqrt(np.clip(depth_km, 0, None))
    xs = np.array([np.sqrt(p[0]) for p in OCEAN_RAMP])
    out = np.zeros(depth_km.shape + (3,))
    for c in range(3):
        out[..., c] = np.interp(d, xs, [p[1][c] for p in OCEAN_RAMP])
    return out


def hillshade(elev_km, lat_deg, radius_km, exaggeration=14.0):
    """Standard hillshade on the equirectangular grid, sun from the north-west."""
    dy_km = np.pi * radius_km / HEIGHT
    coslat = np.clip(np.cos(np.deg2rad(lat_deg)), 0.08, None)[:, None]
    dx_km = 2 * np.pi * radius_km * coslat / WIDTH

    e = elev_km * exaggeration
    dzdx = (np.roll(e, -1, axis=1) - np.roll(e, 1, axis=1)) / (2 * dx_km)
    dzdy = np.empty_like(e)
    dzdy[1:-1] = (e[:-2] - e[2:]) / (2 * dy_km)
    dzdy[0] = dzdy[1]
    dzdy[-1] = dzdy[-2]

    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(dzdy, -dzdx)
    az, alt = np.deg2rad(315.0), np.deg2rad(45.0)
    shade = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(shade, 0, 1)


def main():
    import argparse
    global WIDTH, HEIGHT
    ap = argparse.ArgumentParser(description=__doc__ or "Render the basemap")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="REGULAR climatology file to tint from, overriding "
                         "config/planet.yaml's baseline_climatology")
    ap.add_argument("--width", type=int, default=WIDTH,
                    help=f"equirectangular width in pixels, height is half of "
                         f"it (default {WIDTH})")
    args = ap.parse_args()
    if args.width != WIDTH:
        if args.width % 2:
            raise SystemExit("--width must be even: the height is half of it")
        WIDTH, HEIGHT = args.width, args.width // 2
    if args.climatology is not None:
        global _CLIMATOLOGY_OVERRIDE
        _CLIMATOLOGY_OVERRIDE = args.climatology.resolve()
        # The classification name is derived from this one, so a file that does
        # not carry the suffix would silently resolve the tint to the wrong
        # product instead of failing.
        if not _CLIMATOLOGY_OVERRIDE.name.endswith("_regular_climatology.nc"):
            raise SystemExit(f"--climatology must name a REGULAR climatology "
                             f"(*_regular_climatology.nc), got "
                             f"{_CLIMATOLOGY_OVERRIDE.name}")
        for pth in (_CLIMATOLOGY_OVERRIDE, _classification()):
            if not pth.exists():
                raise SystemExit(f"missing {pth}")

    src = builds.mesh_export()
    manifest = json.loads((src / "manifest.json").read_text())
    radius_km = manifest["planet"]["radiusKm"]
    print(f"build {src.relative_to(ROOT)}  finalElevation {manifest['hashes']['finalElevation'][:8]}")

    print("resolving mesh regions per pixel")
    idx = region_index(src)

    elev = raw(src, "elevation_km", np.float32)[idx].astype(np.float64)
    sclass = raw(src, "surface_class", np.uint8)[idx]
    rock = raw(src, "substrate_class", np.uint8)[idx]
    endo = raw(src, "is_endorheic", np.uint8)[idx].astype(bool)
    land = sclass == 1

    lat_deg, _, _ = pixel_directions()

    print("reading climatology")
    cls = nc.Dataset(_classification())
    biome = np.asarray(cls["biome_index"][:])
    clim_elev_m = np.asarray(cls["surface_elevation"][:])
    cls.close()

    clm = nc.Dataset(_climatology())
    tas = np.asarray(clm["tas"][:])  # (12, 64, 128)
    sic = np.asarray(clm["sic"][:]).mean(axis=0)
    lsm = np.asarray(clm["lsm"][:]).mean(axis=0)
    clm.close()
    warmest = tas.max(axis=0)

    clat = np.asarray(nc.Dataset(_classification())["lat"][:])
    clim_land = lsm > 0.5

    def prepare(field, mask=None):
        """Fill over sea, damp the polar rays, then put it on the render grid."""
        if mask is not None:
            field = fill_ocean_gaps(field, mask)
        return upsample(polar_smooth(field, clat))

    # Biome colours first, then spread over the sea so coastal pixels do not
    # interpolate against an ocean cell.
    biome_rgb = BIOME_RGB[np.clip(biome, 0, len(BIOME_RGB) - 1)]
    land_rgb = np.stack(
        [prepare(biome_rgb[..., c], clim_land & (biome > 0)) for c in range(3)], axis=-1
    )
    warmest_hi = prepare(warmest, clim_land)
    clim_elev_hi = prepare(clim_elev_m / 1000.0, clim_land)
    sic_hi = prepare(sic)

    print("compositing")
    rgb = ocean_colour(-elev)

    # Sea ice over open water, from the annual-mean fraction.
    a = np.clip((sic_hi - 0.12) / 0.45, 0, 1)[..., None] * 0.92
    rgb = rgb * (1 - a) + SEAICE_RGB * a

    # Land.
    lrgb = land_rgb.copy()

    # Bare rock takes over on steep high ground; the biome tint is a lowland
    # statement and the climate grid cannot see a mountain.
    rock_frac = np.clip((elev - 1.8) / 2.2, 0, 1)[..., None]
    lrgb = lrgb * (1 - 0.75 * rock_frac) + ROCK_RGB * (0.75 * rock_frac)

    # Evaporite crust and playa fill, which is where the endorheic drainage put it.
    salt = ((rock == ROCK_EVAPORITE) * 0.85 + (rock == ROCK_PLAYA) * 0.55)[..., None]
    salt_rgb = np.where((rock == ROCK_EVAPORITE)[..., None], SALT_RGB, PLAYA_RGB)
    lrgb = lrgb * (1 - salt) + salt_rgb * salt
    # Closed-basin floors read drier still.
    dry = (endo & land & (elev < 0.2))[..., None] * 0.18
    lrgb = lrgb * (1 - dry) + PLAYA_RGB * dry

    # Permanent snow and ice: warmest month below freezing after a lapse-rate
    # correction from the CLIMATOLOGY'S OWN orography, `clim_elev_hi`, to this
    # raster's. Neither rung is named here: the climatology carries its
    # `surface_elevation` and the rung it was run at is whichever one produced
    # the file `--climatology` points at, so a truncation written into this
    # comment goes stale the first time the ladder moves. The rate is measured,
    # not Earth's 6.5 (PHYS-12), and it is the warm-season one because the field
    # being extrapolated is the warmest month.
    #
    # THIS IS A RENDER AND NOT A REDUCTION. It is the same freezing-height
    # criterion the ice mask uses, evaluated per PIXEL against the mesh region
    # under it, and it stays here rather than reading
    # `hydrography/data/<build>/support_<grid>.nc`: that artifact holds the
    # distribution inside a MODEL cell, and what this needs is the elevation at
    # a point. GRID-2 counted this as one of three implementations of one
    # criterion; the two that had to converge are the mask and the statistic,
    # and they now both come from the support artifact.
    # The file is passed rather than left to default: lapse.py resolves its own
    # from config/planet.yaml, which would ignore --climatology and take the
    # rate from a different climatology than the tint.
    lapse_k_per_km = lapse.environmental_lapse_k_per_km(
        season="warmest", climatology_file=_climatology())
    t_surface = warmest_hi - lapse_k_per_km * (elev - clim_elev_hi)
    ice = np.clip((273.15 - t_surface) / 4.0, 0, 1)[..., None]
    lrgb = lrgb * (1 - ice) + ICE_RGB * ice

    rgb = np.where(land[..., None], lrgb, rgb)

    # Lakes and rivers, which are the one part of the colour that is a result
    # rather than a tint: a water balance solved in hydrography/, not a guess.
    water = np.zeros(elev.shape, bool)
    water_path = surface_water_path()
    surface_water = load_surface_water(idx)
    if surface_water is not None:
        lake, lake_depth, discharge = surface_water

        lrgb = lake_colour(lake_depth)
        # A lake in a place cold enough to hold permanent snow is frozen.
        lrgb = lrgb * (1 - ice) + ICE_RGB * ice
        rgb = np.where((lake & land)[..., None], lrgb, rgb)

        # Rivers are drawn at the resolution the mesh actually has: one region
        # across, about 15 km, so their width is not information. Weight is,
        # so the blend follows discharge rather than the line getting fatter.
        weight = np.clip(
            (np.log10(np.maximum(discharge, 1.0)) - np.log10(RIVER_MIN_M3_S))
            / (np.log10(RIVER_FULL_M3_S) - np.log10(RIVER_MIN_M3_S)),
            0, 1,
        )
        # Any river drawn at all is drawn solidly. A line one region wide at 20%
        # opacity is not a river, it is a smudge; discharge sets the last half.
        river = np.where(weight > 0, 0.55 + 0.45 * weight, 0.0) * land * ~lake
        rgb = rgb * (1 - river[..., None]) + RIVER_RGB * river[..., None]
        water = (lake & land) | (river > 0.35)

    # Relief. Standing water is flat, so it takes the muted treatment the sea
    # gets rather than the land hillshade, which would emboss a lake surface.
    shade = hillshade(elev, lat_deg, radius_km)
    lit = 0.62 + 0.76 * shade
    lit = np.where(land & ~water, lit, 0.90 + 0.18 * shade)
    rgb *= lit[..., None]

    rgb = np.clip(rgb, 0, 255).astype(np.uint8)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    np.save(BUILD_DIR / "basemap.npy", rgb)
    # The projection orientations are chosen against this, so that interrupts
    # and cut meridians fall in the sea where the geography allows it.
    np.save(BUILD_DIR / "land.npy", land)

    from PIL import Image

    Image.fromarray(rgb).save(BUILD_DIR / "basemap_equirectangular.png")

    # The identity of everything the raster was drawn from, and the fingerprint
    # over it. `render_projections.py` copies this into the frame it writes
    # rather than re-deriving it from the config, so a frame records the world
    # THIS raster holds and not the one the config currently names.
    inputs = frames.inputs(
        source_build=src.parent.name,
        terrain_hash=manifest["hashes"]["finalElevation"],
        climatology=_climatology(),
        classification=_classification(),
        surface_water=(water_path if water_path.exists() else None),
    )

    prov = {
        "source_build": src.parent.name,
        "source_export": str(src.relative_to(ROOT)),
        "final_elevation_hash": manifest["hashes"]["finalElevation"],
        "climatology": str(_climatology().relative_to(ROOT)),
        # Derived from the product actually read, not fixed: the caveat a reader
        # needs is which climatology this is, and a bootstrap's biomes are not a
        # baseline's.
        "climatology_caveat": _caveat(),
        "inputs": inputs,
        "fingerprint": frames.fingerprint(inputs),
        "resolution": [WIDTH, HEIGHT],
        "planet_radius_km": radius_km,
        "numpy": np.__version__,
        "git_commit": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True
        ).stdout.strip(),
    }
    (BUILD_DIR / "basemap_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {BUILD_DIR / 'basemap.npy'}  {rgb.shape}")


if __name__ == "__main__":
    main()
