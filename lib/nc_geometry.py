"""The geometry a NetCDF product declares about itself, and the reader that refuses a wrong one.

WHAT THIS IS FOR. `ncwa`, `iris`, `xarray` and everything built on them assume
Earth when a file does not say otherwise, and each assumption returns an
ordinary-looking number rather than an error. A file with no coordinate system
gets a 6371 km sphere and so a cell area wrong by the square of the radius
ratio. A file with a `lat` axis and no quadrature weight gets an UNWEIGHTED mean
over Gaussian rows, which is wrong by about 1.8 per cent in the polar row and
converges to that rather than to zero. A file with `units = "deg"` rather than
the CF axis names gets its axes guessed from position.
`docs/src/reference/environment.md` enumerates those three traps; this module is
what a writer calls so its own file is not subject to them.

WHAT IT DOES NOT DO. It is not a remapper and not a second grid convention.
Every edge, every node, every label and every weight here is CONSTRUCTED BY
`lib/gridding.py`, which owns the one grid convention; this module writes what
that module returns into a file and reads it back. It does not invent a
metadata vocabulary either: `lib/spatial_support.py` owns the versioned identity
and semantics contract, and `geometry_sha256` below is that module's
`canonical_digest` over that module's geometry payload, so a file's stamp and an
artifact's contract carry the same identity for the same grid.

THE TWO LONGITUDE CONVENTIONS, WHICH IS WHY A NAME AND NOT A NUMBER. Orogen
labels the columns -180..180 and ExoPlaSim labels the SAME columns 0..360 from
index 0. Every derived number a reader could check agrees between them: the
column count, the column width, the width of a bounds pair, the sum of the area
weights, the pole-to-pole span. What differs is the origin, and an export centre
offered as a model label lands exactly half a column out, which matches every
cell half a planet away instead of matching none. So the declaration carries the
convention BY NAME, `read_grid` reconstructs the named convention's axis from
`lib/gridding.py` and refuses a file whose axis is not it, and a consumer holding
two files compares two names rather than two axes. CLAUDE.md rule 3.

THE PLANET IS READ, NEVER COPIED. `planet_radius_m` is the one expression that
turns `config/planet.yaml`'s `planet.radius_earth` into metres and
`sphere_area_m2` the one that forms the area of the sphere from it, against the
several spellings of the conversion each script used to carry. CLAUDE.md
rule 2.

TWO SUPPORTS, BECAUSE THIS PROJECT WRITES TWO. `declare_grid` is for a product on
a model grid, which is small enough that the per-cell weight and area belong in
the file. `declare_region_mesh` is for a product on the Orogen mesh, whose ten
million region areas live in the durable export and are POINTED AT rather than
copied: an eighty-megabyte area array duplicated into every derived product is
the derived-quantity-written-down-somewhere defect at scale, and what a mesh
product's reader needs is to be told that its `region` axis is not a grid and
that an unweighted reduction over it is not an area mean.
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import numpy as np
import yaml

import gridding
import spatial_support

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PLANET_CONFIG = PROJECT_ROOT / "config" / "planet.yaml"

EARTH_RADIUS_M = 6_371_000.0
"""Earth's mean radius, the unit `planet.radius_earth` is expressed in.

Here because `config/planet.yaml` states the radius as a RATIO and something has
to hold the denominator. It is a unit conversion and NOT a planetary parameter:
nothing in this project may use it as a radius. It had been spelled `6371e3`,
`6.371e6` and `6_371_000.0` in different scripts, several of which formed the
area of the sphere from it inline, which is what `sphere_area_m2` is for."""

DECLARATION_VERSION = "vesper-nc-geometry/1"
"""Stamped as `vesper_geometry_declaration`. A reader keys on this, not on the
presence of any one attribute, so a file written before the declaration existed
is distinguishable from one written after it and not merely missing a field."""

MODEL_LABELS = "model_labels"
"""ExoPlaSim's own labels: column `k` is `k * 360 / nlon` degrees east from
index 0, so column 0 is labelled 0. Every climatology carries this axis, and so
does every product built by copying a climatology's `lon` array."""

EXPORT_CENTRES = "export_centres"
"""The export's own centres: column `k` spans `[-180 + k*dlon, -180 + (k+1)*dlon)`
and is labelled at its midpoint, so column 0 is labelled just east of -180. Every
product binned from the mesh through `lib/gridding.py:column` is on this axis,
because that expression is what put each region in its column."""

CONVENTIONS = (MODEL_LABELS, EXPORT_CENTRES)

CRS_VARIABLE = "crs"
WEIGHT_VARIABLE = "cell_weight"
AREA_VARIABLE = "cell_area"

# A label read off a netCDF float32 axis agrees with the constructed axis to a
# few parts in a million and no closer, which is the bar `lib/gridding.py`'s own
# axis doors carry for the same reason. It is four orders tighter than the half
# column that separates the two conventions at the coarsest rung on the ladder
# (2.8 degrees at T21), which is the thing it exists to catch.
_AXIS_TOLERANCE_DEG = 1.0e-3


class DeclaredGrid(NamedTuple):
    """What a declared file says its grid is, once every claim in it has held.

    `spec` is a `gridding.GridSpec` CONSTRUCTED from the declared row and column
    counts, not assembled from the axes on disk: `GridSpec.source` refuses a spec
    built out of labels read off a file, and reconstructing a coordinate from
    labels is what rule 3 forbids. The file's axes are checked AGAINST it.
    """

    spec: gridding.GridSpec
    convention: str
    radius_m: float
    geometry_sha256: str

    @property
    def nlat(self) -> int:
        return self.spec.nlat

    @property
    def nlon(self) -> int:
        return self.spec.nlon


def planet_config(path: Path = PLANET_CONFIG) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def planet_radius_m(config: dict | None = None) -> float:
    """This planet's radius in metres, from `config/planet.yaml` and nowhere else.

    The one expression. A caller that has already loaded the config passes it;
    one that has not gets it loaded here. CLAUDE.md rule 2: the value is read at
    the point of use and never copied into a script or into prose.
    """
    cfg = config if config is not None else planet_config()
    return float(cfg["planet"]["radius_earth"]) * EARTH_RADIUS_M


def sphere_area_m2(config: dict | None = None) -> float:
    """Area of the whole sphere, in square metres.

    Here so that no caller forms `4 * pi * r**2` from a radius of its own. Two
    already had, from a hardcoded radius.
    """
    return 4.0 * np.pi * planet_radius_m(config) ** 2


def longitude_axis(convention: str, nlon: int) -> np.ndarray:
    """The longitude axis of the named convention, constructed by `lib/gridding.py`.

    The ONE place the two conventions are named beside each other, and neither
    is derived here: `model_longitude_labels` and `GridSpec.cell_centres` are
    both `lib/gridding.py`'s. A convention that is not one of the two is refused
    rather than defaulted, because defaulting is how a half-column rotation gets
    written into a file that then looks right.
    """
    if convention == MODEL_LABELS:
        return gridding.model_longitude_labels(int(nlon))
    if convention == EXPORT_CENTRES:
        return gridding.gaussian_grid(1, int(nlon)).cell_centres()[1]
    raise ValueError(
        f"{convention!r} is not a longitude convention. It is {MODEL_LABELS!r} "
        f"for a model's own labels or {EXPORT_CENTRES!r} for the export's "
        "centres, and the two are different names for the same columns rather "
        "than a coordinate transform. CLAUDE.md rule 3.")


def longitude_bounds(convention: str, spec: gridding.GridSpec) -> np.ndarray:
    """Cell boundaries in the named convention's frame, shaped `(nlon, 2)`.

    Both conventions bound a column at half a column width either side of its
    own label, so the WIDTHS are identical and only the origin differs. That is
    exactly why the file also carries the convention by name: a reader comparing
    bounds alone cannot tell the two apart.
    """
    labels = longitude_axis(convention, spec.nlon)
    half = 0.5 * spec.dlon
    return np.stack([labels - half, labels + half], axis=1)


def latitude_bounds(spec: gridding.GridSpec) -> np.ndarray:
    """Row boundaries in degrees, shaped `(nlat, 2)`, in the grid's own row order.

    Taken from the spec's sine-of-latitude edges, which for a Gaussian grid ARE
    the quadrature weights laid end to end. A bounds pair in degrees is for a
    reader; the weight variable beside it is what a weighted mean must use,
    because the rows are equal in SINE and not in degrees.
    """
    edges = np.rad2deg(np.arcsin(np.clip(spec.sin_edges, -1.0, 1.0)))
    return np.stack([edges[:-1], edges[1:]], axis=1)


def geometry_digest(spec: gridding.GridSpec, convention: str) -> str:
    """The identity of this grid, through `lib/spatial_support.py`'s digest.

    The payload carries the convention as well as the edges, so the two axes
    that agree on every number a reader could compare come out with DIFFERENT
    identities. Without the convention in the payload the digest would be one
    more number both conventions produce.
    """
    return spatial_support.canonical_digest({
        "kind": "atmosphere_gaussian_grid",
        "shape": [int(spec.nlat), int(spec.nlon)],
        "constructor": str(spec.source),
        "longitude_convention": convention,
        "longitude_edges": list(map(float, spec.lon_edges)),
        "sine_latitude_edges": list(map(float, spec.sin_edges)),
    })


def _existing_axis(ds, name):
    if name not in ds.variables:
        return None
    return np.asarray(ds.variables[name][:], dtype=np.float64)


def _detect_convention(lon_axis: np.ndarray, nlon: int, what: str) -> str:
    """Which convention this axis IS, or a refusal naming both.

    Detection rather than assumption, and it refuses an axis that is neither
    instead of picking the closer one. An axis a third of a column out of true
    is not a slightly wrong model axis; it is an axis nothing in this project
    constructs, and snapping it to the nearest convention is the failure this
    whole module exists to make impossible.
    """
    offsets = {}
    for convention in CONVENTIONS:
        want = longitude_axis(convention, nlon)
        offsets[convention] = float(np.abs(lon_axis - want).max())
    best = min(offsets, key=offsets.get)
    if offsets[best] <= _AXIS_TOLERANCE_DEG:
        return best
    raise SystemExit(
        f"{what}: the {nlon}-column longitude axis on this file is neither of "
        f"the project's two conventions. It is {offsets[MODEL_LABELS]:.4g} "
        f"degrees from the model's labels and {offsets[EXPORT_CENTRES]:.4g} "
        "degrees from the export's centres. Do not adopt the nearer one: "
        "construct the axis from lib/gridding.py in the convention the writer "
        "actually used. CLAUDE.md rule 3.")


def declare_grid(ds, *, nlat: int | None = None, nlon: int | None = None,
                 convention: str | None = None, lat_dim: str = "lat",
                 lon_dim: str = "lon", config: dict | None = None,
                 what: str = "this product") -> DeclaredGrid:
    """Write this file's geometry into it. Call LAST, inside the open Dataset.

    Last, because it attaches `grid_mapping` and `cell_measures` to every
    variable that spans the two grid dimensions, and a variable created after
    this returns does not get them.

    It creates the `lat` and `lon` coordinate variables when the file has none
    -- several products carry the dimensions and no axis at all, which is the
    worst of the three cases, since a reader then has nothing to be wrong about
    -- and when the file already has them it CHECKS them rather than replacing
    them. Either way `convention` is required for a file with no longitude axis
    and verified for a file with one; there is no default, because a default is
    how a half-column rotation becomes a stated fact.

    Returns the `DeclaredGrid` it wrote, so a caller can put the identity in its
    own provenance record without re-deriving it.
    """
    nlat = int(len(ds.dimensions[lat_dim]) if nlat is None else nlat)
    nlon = int(len(ds.dimensions[lon_dim]) if nlon is None else nlon)
    spec = gridding.gaussian_grid(nlat, nlon)

    lat_axis = _existing_axis(ds, lat_dim)
    if lat_axis is None:
        lat_axis = gridding.gaussian_latitudes(nlat)
        ds.createVariable(lat_dim, "f8", (lat_dim,))[:] = lat_axis
    else:
        # Through `lib/gridding.py`'s own door, which refuses an equally spaced
        # axis. That refusal is the second trap caught at the write door rather
        # than in a consumer: a file whose rows are not Gaussian must not be
        # stamped with Gaussian quadrature weights.
        gridding.gaussian_row_weights(lat_axis, what=f"{what}: the {lat_dim} axis")

    lon_axis = _existing_axis(ds, lon_dim)
    if lon_axis is None:
        if convention is None:
            raise ValueError(
                f"{what} has a {lon_dim} dimension and no {lon_dim} variable, so "
                "the convention its columns are on cannot be detected and must "
                f"be declared. Pass convention={MODEL_LABELS!r} for a product "
                "carrying a model's own labels or "
                f"convention={EXPORT_CENTRES!r} for one binned from the mesh "
                "through lib/gridding.py:column.")
        lon_axis = longitude_axis(convention, nlon)
        ds.createVariable(lon_dim, "f8", (lon_dim,))[:] = lon_axis
    else:
        found = _detect_convention(lon_axis, nlon, f"{what}: the {lon_dim} axis")
        if convention is not None and convention != found:
            raise SystemExit(
                f"{what} declares its columns as {convention!r} and carries "
                f"{found!r}. These are different names for the same columns, so "
                "the file would agree with itself on every number and disagree "
                "with the world by half a grid. CLAUDE.md rule 3.")
        convention = found

    lat_var, lon_var = ds.variables[lat_dim], ds.variables[lon_dim]
    lat_var.units = "degrees_north"
    lat_var.standard_name = "latitude"
    lat_var.long_name = "latitude of the Gaussian node"
    lat_var.axis = "Y"
    lat_var.bounds = f"{lat_dim}_bnds"
    lon_var.units = "degrees_east"
    lon_var.standard_name = "longitude"
    lon_var.axis = "X"
    lon_var.bounds = f"{lon_dim}_bnds"
    lon_var.long_name = (
        "longitude on ExoPlaSim's own label axis, degrees east from column 0"
        if convention == MODEL_LABELS else
        "longitude at the export's own column centres, -180 to 180")
    lon_var.setncattr("vesper_longitude_convention", convention)

    if "bnds" not in ds.dimensions:
        ds.createDimension("bnds", 2)
    for name, values in ((f"{lat_dim}_bnds", latitude_bounds(spec)),
                         (f"{lon_dim}_bnds", longitude_bounds(convention, spec))):
        dim = lat_dim if name.startswith(lat_dim) else lon_dim
        if name not in ds.variables:
            ds.createVariable(name, "f8", (dim, "bnds"))
        ds.variables[name][:] = values

    radius = planet_radius_m(config)
    if CRS_VARIABLE not in ds.variables:
        ds.createVariable(CRS_VARIABLE, "i4", ())
    crs = ds.variables[CRS_VARIABLE]
    crs[...] = 0
    crs.grid_mapping_name = "latitude_longitude"
    # `earth_radius` is the attribute iris reads to build a spherical GeogCS,
    # and a GeogCS is what makes its `area_weights` use this radius instead of
    # falling back to Earth's default. The two axis attributes are the same
    # sphere spelled the way a projection library expects to find it.
    crs.earth_radius = radius
    crs.semi_major_axis = radius
    crs.semi_minor_axis = radius
    crs.inverse_flattening = 0.0
    crs.longitude_of_prime_meridian = 0.0
    crs.long_name = ("the sphere this world is, from config/planet.yaml's "
                     "planet.radius_earth. It is NOT Earth's, and a reader that "
                     "substitutes Earth's is wrong by the square of the ratio")

    weights = spec.cell_area_fraction()
    if WEIGHT_VARIABLE not in ds.variables:
        ds.createVariable(WEIGHT_VARIABLE, "f8", (lat_dim, lon_dim))
    weight = ds.variables[WEIGHT_VARIABLE]
    weight[:] = weights
    weight.units = "1"
    weight.long_name = ("share of the sphere in each cell: the Gauss-Legendre "
                        "quadrature weight of the row divided by the columns. "
                        "Sums to one over the grid. THIS IS WHAT A GLOBAL MEAN "
                        "IS WEIGHTED BY -- ncwa -w cell_weight -a lat,lon -- "
                        "and cos(lat) is not it")

    if AREA_VARIABLE not in ds.variables:
        ds.createVariable(AREA_VARIABLE, "f8", (lat_dim, lon_dim))
    area = ds.variables[AREA_VARIABLE]
    area[:] = spec.cell_area(radius)
    area.units = "m2"
    area.standard_name = "cell_area"
    area.grid_mapping = CRS_VARIABLE
    area.long_name = ("area of each cell on this planet's sphere; sums to "
                      "4 pi R^2 for the R the crs variable declares")

    digest = geometry_digest(spec, convention)
    ds.setncattr("Conventions", "CF-1.10")
    ds.setncattr("vesper_geometry_declaration", DECLARATION_VERSION)
    ds.setncattr("vesper_support_kind", "atmosphere_gaussian_grid")
    ds.setncattr("vesper_grid_latitudes", nlat)
    ds.setncattr("vesper_grid_longitudes", nlon)
    ds.setncattr("vesper_longitude_convention", convention)
    ds.setncattr("vesper_geometry_sha256", digest)
    ds.setncattr("vesper_planet_radius_m", radius)
    ds.setncattr("vesper_grid_constructor", str(spec.source))

    reserved = {lat_dim, lon_dim, f"{lat_dim}_bnds", f"{lon_dim}_bnds",
                CRS_VARIABLE, WEIGHT_VARIABLE, AREA_VARIABLE}
    for name, var in ds.variables.items():
        if name in reserved or not {lat_dim, lon_dim} <= set(var.dimensions):
            continue
        var.grid_mapping = CRS_VARIABLE
        var.cell_measures = f"area: {AREA_VARIABLE}"

    return DeclaredGrid(spec=spec, convention=convention, radius_m=radius,
                        geometry_sha256=digest)


def declare_region_mesh(ds, *, n_regions: int, terrain_hash: str,
                        region_dim: str = "region", config: dict | None = None,
                        area_variable: str = "cell_area",
                        measure: str | None = None) -> str:
    """Declare a product on the Orogen mesh: what its axis is and what it is not.

    A mesh product's `region` axis is not a grid and its regions are not equal
    in area, so an unweighted reduction over it is not an area mean and no tool
    can tell without being told. The areas themselves are NOT copied in: they
    are ten million numbers that already exist in the durable export, and a copy
    in every derived product is a derived quantity written down somewhere,
    frozen the moment the terrain moves. So the file names where they are and
    the reader fetches them from the build the file already stamps.

    `measure` overrides the sentence that says where the axis's measure lives,
    for a file whose axis is not the mesh's own -- a basin catalogue indexes the
    mesh but a basin's area is its catchment, which the file itself carries.

    Returns the geometry digest, for the caller's provenance record.
    """
    radius = planet_radius_m(config)
    if CRS_VARIABLE not in ds.variables:
        ds.createVariable(CRS_VARIABLE, "i4", ())
    crs = ds.variables[CRS_VARIABLE]
    crs[...] = 0
    crs.grid_mapping_name = "latitude_longitude"
    crs.earth_radius = radius
    crs.semi_major_axis = radius
    crs.semi_minor_axis = radius
    crs.inverse_flattening = 0.0
    crs.longitude_of_prime_meridian = 0.0
    crs.long_name = ("the sphere this world is, from config/planet.yaml's "
                     "planet.radius_earth; any length or area derived from "
                     "this file is on it and not on Earth's")

    digest = spatial_support.canonical_digest({
        "kind": "orogen_region_mesh",
        "shape": [int(n_regions)],
        "constructor": "World Orogen export mesh, region order as exported",
        "terrain_hash": str(terrain_hash),
    })
    ds.setncattr("Conventions", "CF-1.10")
    ds.setncattr("vesper_geometry_declaration", DECLARATION_VERSION)
    ds.setncattr("vesper_support_kind", "orogen_region_mesh")
    ds.setncattr("vesper_region_count", int(n_regions))
    ds.setncattr("vesper_geometry_sha256", digest)
    ds.setncattr("vesper_planet_radius_m", radius)
    ds.setncattr("vesper_region_measure", measure or (
        f"NOT IN THIS FILE. Region areas are the export's {area_variable}, on "
        f"the mesh of terrain_hash {terrain_hash}. The regions are UNEQUAL in "
        f"area, so a reduction over the {region_dim} axis without them -- "
        f"ncwa -a {region_dim}, a bare mean, a histogram of region counts -- is "
        "a mean over regions and not over the world."))
    return digest


def read_grid(ds, *, lat_dim: str = "lat", lon_dim: str = "lon",
              config: dict | None = None,
              what: str = "this file") -> DeclaredGrid:
    """The grid a file declares, with every claim in it checked. Refuses, never guesses.

    What it refuses, and each of these has a right answer rather than a
    difference: a file with no declaration; a convention that is not one of the
    two; an axis that is not the named convention's; bounds that are not the
    constructed spec's; a declared radius that is not this planet's; a digest
    that is not the digest of what the file actually carries.
    """
    attrs = set(ds.ncattrs())
    if "vesper_geometry_declaration" not in attrs:
        raise SystemExit(
            f"{what} carries no geometry declaration, so nothing in it says what "
            "sphere, what rows or which longitude convention it is on. A reader "
            "would supply Earth's. Regenerate it with a writer that calls "
            "lib/nc_geometry.py:declare_grid.")
    version = ds.getncattr("vesper_geometry_declaration")
    if version != DECLARATION_VERSION:
        raise SystemExit(
            f"{what} declares geometry version {version!r} and this reader is "
            f"{DECLARATION_VERSION!r}. The declaration is versioned so that a "
            "change to it is a refusal rather than a field read under the wrong "
            "meaning.")
    if ds.getncattr("vesper_support_kind") != "atmosphere_gaussian_grid":
        raise SystemExit(
            f"{what} is on the {ds.getncattr('vesper_support_kind')!r} support "
            "and not a grid. There is no grid in it to read.")

    nlat = int(ds.getncattr("vesper_grid_latitudes"))
    nlon = int(ds.getncattr("vesper_grid_longitudes"))
    for dim, count in ((lat_dim, nlat), (lon_dim, nlon)):
        if dim in ds.dimensions and len(ds.dimensions[dim]) != count:
            raise SystemExit(
                f"{what} declares {count} for {dim} and its {dim} dimension is "
                f"{len(ds.dimensions[dim])}. The file disagrees with itself.")
    spec = gridding.gaussian_grid(nlat, nlon)

    convention = ds.getncattr("vesper_longitude_convention")
    if convention not in CONVENTIONS:
        raise SystemExit(
            f"{what} declares longitude convention {convention!r}, which is "
            f"neither {MODEL_LABELS!r} nor {EXPORT_CENTRES!r}.")

    lon_axis = _existing_axis(ds, lon_dim)
    if lon_axis is None:
        raise SystemExit(f"{what} declares a grid and carries no {lon_dim} variable.")
    want = longitude_axis(convention, nlon)
    off = float(np.abs(lon_axis - want).max())
    if off > _AXIS_TOLERANCE_DEG:
        raise SystemExit(
            f"{what} declares its columns as {convention!r} and its {lon_dim} "
            f"axis is up to {off:.4g} degrees from that axis. Half a column out "
            "is the export's centres offered as the model's labels, which "
            "matches every cell half a planet away rather than none of them. "
            "CLAUDE.md rule 3.")
    lat_axis = _existing_axis(ds, lat_dim)
    if lat_axis is None:
        raise SystemExit(f"{what} declares a grid and carries no {lat_dim} variable.")
    gridding.gaussian_row_weights(lat_axis, what=f"{what}: the {lat_dim} axis")

    for name, expected in ((f"{lat_dim}_bnds", latitude_bounds(spec)),
                           (f"{lon_dim}_bnds", longitude_bounds(convention, spec))):
        if name not in ds.variables:
            raise SystemExit(
                f"{what} declares a grid and carries no {name}. Without bounds "
                "a reader has centres only, and centres do not say whether the "
                "rows are Gaussian or equally spaced.")
        got = np.asarray(ds.variables[name][:], dtype=np.float64)
        if got.shape != expected.shape or not np.allclose(
                got, expected, rtol=0.0, atol=1e-9):
            raise SystemExit(
                f"{what}: {name} is not the boundary of the {nlat}x{nlon} grid "
                "it declares. The file's bounds and its shape describe two "
                "different grids.")

    radius = planet_radius_m(config)
    declared_radius = float(ds.getncattr("vesper_planet_radius_m"))
    if abs(declared_radius - radius) > 1.0:
        raise SystemExit(
            f"{what} declares a sphere of {declared_radius:.1f} m and "
            f"config/planet.yaml's planet is {radius:.1f} m. Every area in the "
            "file is off by the square of the ratio; it is a product of another "
            "planet or of a superseded config.")

    digest = geometry_digest(spec, convention)
    if ds.getncattr("vesper_geometry_sha256") != digest:
        raise SystemExit(
            f"{what} carries a geometry identity that is not the identity of "
            "the grid it describes. Something rewrote one of the two.")
    return DeclaredGrid(spec=spec, convention=convention, radius_m=radius,
                        geometry_sha256=digest)


def require_grid(ds, nlat: int, nlon: int, *, lat_dim: str = "lat",
                 lon_dim: str = "lon", convention: str | None = None,
                 config: dict | None = None,
                 what: str = "this file") -> DeclaredGrid:
    """The declared grid, refused unless it is the one the caller asked for.

    THE DOOR A CONSUMER USES. A product built on one rung read as another, or a
    model-labelled file read as an export-centred one, is a wrong answer that
    survives every check made of the numbers themselves. Both are refused here,
    by comparing what the file DECLARES against what the caller NAMES, which is
    a comparison with a right answer rather than a difference.
    """
    grid = read_grid(ds, lat_dim=lat_dim, lon_dim=lon_dim, config=config, what=what)
    if (grid.nlat, grid.nlon) != (int(nlat), int(nlon)):
        raise SystemExit(
            f"{what} is a {grid.nlat}x{grid.nlon} product and it is being read "
            f"as {int(nlat)}x{int(nlon)}. These are different rungs of the "
            "ladder in lib/rungs.py and no mapping between them is defined here.")
    if convention is not None and grid.convention != convention:
        raise SystemExit(
            f"{what} is on {grid.convention!r} and is being read as "
            f"{convention!r}. The two label the same columns differently; map "
            "by index and share one coordinate source. CLAUDE.md rule 3.")
    return grid


# ---------------------------------------------------------------------------
# The selftest. Every check below has an answer fixed before it runs: an
# identity the round trip must reproduce, a closure against config/planet.yaml,
# or a refusal that must fire. The negative ones carry the weight -- a reader
# that accepted everything would pass every positive check here.
# ---------------------------------------------------------------------------

_CHECKS = 16


def _selftest() -> int:
    import tempfile
    from netCDF4 import Dataset

    problems: list[str] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    def refuses(name: str, fn) -> None:
        try:
            fn()
        except (SystemExit, ValueError) as exc:
            print(f"[  ok  ] {name}")
            print(f"         {str(exc).splitlines()[0][:96]}")
            return
        check(name, False, "it was accepted")

    radius = planet_radius_m()
    total_area = sphere_area_m2()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # A fine-rung product carrying a model-labelled axis, written the way a
        # climatology-derived product is: the axis already on the file.
        fine = tmp / "model_labelled.nc"
        with Dataset(fine, "w") as ds:
            ds.createDimension("lat", 64)
            ds.createDimension("lon", 128)
            ds.createVariable("lat", "f8", ("lat",))[:] = gridding.gaussian_latitudes(64)
            ds.createVariable("lon", "f8", ("lon",))[:] = gridding.model_longitude_labels(128)
            ds.createVariable("field", "f4", ("lat", "lon"))[:] = 1.0
            written = declare_grid(ds, what="the fine fixture")

        check("a model-labelled axis is detected as model_labels",
              written.convention == MODEL_LABELS, written.convention)

        with Dataset(fine) as ds:
            read = read_grid(ds, what="the fine fixture")
            weights = np.asarray(ds.variables[WEIGHT_VARIABLE][:], dtype=np.float64)
            areas = np.asarray(ds.variables[AREA_VARIABLE][:], dtype=np.float64)
            crs_radius = float(ds.variables[CRS_VARIABLE].earth_radius)
            field_attrs = set(ds.variables["field"].ncattrs())
            lon_units = ds.variables["lon"].units
            lat_units = ds.variables["lat"].units

        # THE ROUND TRIP. A file read back must reproduce the grid it was
        # written from, exactly: same edges, same identity, same convention.
        check("the round trip reproduces the longitude edges exactly",
              bool(np.array_equal(read.spec.lon_edges, written.spec.lon_edges)))
        check("the round trip reproduces the sine-latitude edges exactly",
              bool(np.array_equal(read.spec.sin_edges, written.spec.sin_edges)))
        check("the round trip reproduces the geometry identity",
              read.geometry_sha256 == written.geometry_sha256)
        check("the round trip reproduces the convention",
              read.convention == written.convention)

        # THE CLOSURE. The declared areas sum to the sphere config/planet.yaml
        # describes, and the weights to one. Both are exact statements about a
        # quadrature, not tolerances chosen to pass.
        residual = abs(areas.sum() - total_area) / total_area
        check("the declared cell areas sum to 4 pi R^2 from config/planet.yaml",
              residual <= 1e-12, f"relative residual {residual:.3g}")
        check("the declared cell weights sum to one",
              abs(weights.sum() - 1.0) <= 1e-12,
              f"{weights.sum() - 1.0:.3g} away")
        check("the declared weights are lib/gridding.py's area weights",
              bool(np.array_equal(weights, gridding.gaussian_area_weights(
                  gridding.gaussian_latitudes(64), 128))))
        check("the crs declares this planet's radius and not Earth's",
              abs(crs_radius - radius) < 1.0 and abs(crs_radius - EARTH_RADIUS_M) > 1.0,
              f"{crs_radius:.1f} m")

        # THE CF NAMES, which are what stops an axis being guessed from position.
        check("the axes carry the CF names and not units = deg",
              (lat_units, lon_units) == ("degrees_north", "degrees_east"),
              f"{lat_units!r}, {lon_units!r}")
        check("a data variable is bound to the crs and the area measure",
              {"grid_mapping", "cell_measures"} <= field_attrs,
              str(sorted(field_attrs)))

        # A coarse-rung product with NO axes on the file, written the way a mesh-binned
        # product is: the convention declared and the axis constructed.
        coarse = tmp / "export_centred.nc"
        with Dataset(coarse, "w") as ds:
            ds.createDimension("lat", 32)
            ds.createDimension("lon", 64)
            ds.createVariable("field", "f4", ("lat", "lon"))[:] = 0.0
            declare_grid(ds, convention=EXPORT_CENTRES, what="the coarse fixture")
        with Dataset(coarse) as ds:
            export_grid = read_grid(ds, what="the coarse fixture")
        check("a file with no axes gets the declared convention's axis",
              export_grid.convention == EXPORT_CENTRES
              and (export_grid.nlat, export_grid.nlon) == (32, 64))

        # THE NEGATIVE CONTROLS.
        with Dataset(coarse) as ds:
            refuses("a file read as another rung is refused",
                    lambda: require_grid(ds, 64, 128, what="the coarse fixture"))
            refuses("an export-centred file read as model-labelled is refused",
                    lambda: require_grid(ds, 32, 64, convention=MODEL_LABELS,
                                         what="the coarse fixture"))

        # A file whose axis is the export's centres while its declaration says
        # the model's labels: every number in it is plausible and it is rotated
        # by half the planet. This is the case the declaration exists for.
        lying = tmp / "half_a_grid_out.nc"
        with Dataset(lying, "w") as ds:
            ds.createDimension("lat", 32)
            ds.createDimension("lon", 64)
            ds.createVariable("field", "f4", ("lat", "lon"))[:] = 0.0
            declare_grid(ds, convention=EXPORT_CENTRES, what="the rotated fixture")
        with Dataset(lying, "a") as ds:
            ds.setncattr("vesper_longitude_convention", MODEL_LABELS)
        with Dataset(lying) as ds:
            refuses("an axis half a column from the convention it declares is refused",
                    lambda: read_grid(ds, what="the rotated fixture"))

        # An equally spaced latitude axis is not a Gaussian one, and must not be
        # stamped with Gaussian quadrature weights.
        equal = tmp / "equally_spaced.nc"
        with Dataset(equal, "w") as ds:
            ds.createDimension("lat", 32)
            ds.createDimension("lon", 64)
            ds.createVariable("lat", "f8", ("lat",))[:] = np.linspace(87.0, -87.0, 32)
            ds.createVariable("lon", "f8", ("lon",))[:] = gridding.model_longitude_labels(64)
            refuses("an equally spaced latitude axis is refused at the write door",
                    lambda: declare_grid(ds, what="the equally spaced fixture"))

        # An undeclared file is refused rather than read with Earth substituted.
        bare = tmp / "undeclared.nc"
        with Dataset(bare, "w") as ds:
            ds.createDimension("lat", 32)
            ds.createDimension("lon", 64)
        with Dataset(bare) as ds:
            refuses("a file with no declaration is refused rather than assumed",
                    lambda: read_grid(ds, what="the undeclared fixture"))

        # The mesh support: declared, and not readable as a grid.
        mesh = tmp / "mesh.nc"
        with Dataset(mesh, "w") as ds:
            ds.createDimension("region", 11)
            ds.createVariable("terminal", "i4", ("region",))[:] = 0
            declare_region_mesh(ds, n_regions=11, terrain_hash="0" * 64)
        with Dataset(mesh) as ds:
            refuses("a mesh product read as a grid is refused",
                    lambda: read_grid(ds, what="the mesh fixture"))

    print(f"\n{_CHECKS} checks, {len(problems)} failed")
    return 1 if problems else 0


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true",
                    help="write, read back and refuse synthetic products; "
                         "needs no build and no run on disk")
    args = ap.parse_args()
    if not args.selftest:
        ap.error("this module is a library; --selftest is the only thing to run")
    return _selftest()


if __name__ == "__main__":
    raise SystemExit(main())
