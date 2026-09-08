World Orogen data export
========================

Seed: 16236323
Regions: 10000005
Grid: 64 x 32 equirectangular

Everything is described in manifest.json: for each field you get its path,
dtype, shape, units and a description. All binaries are flat little-endian
arrays with no header.

Reading a field in Python:

    import json, numpy as np
    m = json.load(open('manifest.json'))
    f = next(x for x in m['grid']['fields'] if x['name'] == 'elevation')
    elev = np.fromfile(f['path'], dtype=f['dtype']).reshape(f['shape'])
    # elev[0] is the northernmost row, elev[:, 0] is longitude -180.

Coordinates: lat = asin(y), lon = atan2(x, z); longitude increases eastward from +z
Elevation: km (elevation_km) / dimensionless (elevation), sea level at 0.
The `elevation` field is the model's internal shaping parameter, not kilometres: land maps to height through a nonlinear hypsometric curve where 0.5 is ~1.1 km and 1.0 is 6 km, and ocean is linear at 10 km per unit. The parameter is not bounded by 1; above it the curve is linear at 6 km per unit, which is where its shape function's domain ends. Use `elevation_km` for physical orography. The land/ocean split is elevation > 0 in either field.

raw/  — one value per mesh region (irregular Delaunay sphere mesh). Region
        positions are in raw/lat.bin, raw/lon.bin and raw/x|y|z.bin; cell areas
        in raw/cell_area.bin. Mesh topology is included so the dual mesh can be
        rebuilt and neighbours walked.
grid/ — the same fields resampled to a regular lat/lon grid. Continuous fields
        are area-weighted means of the regions falling in each cell; categorical
        fields (plate index, boundary type, Köppen class, masks) take the value
        of the region containing the cell centre.
