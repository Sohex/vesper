/**
 * Minimal NetCDF classic writer (CDF-1 and CDF-2 / 64-bit offset).
 *
 * Written by hand rather than pulled from a dependency: the format is small and
 * fully specified, and the export path needs to work from a checkout with no
 * native build step. Only fixed-size (non-record) variables are emitted, which
 * removes the record-striding rules entirely.
 *
 * Reference: https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
 * All header and data values are big-endian; everything is padded to 4 bytes.
 */

const NC_BYTE = 1, NC_CHAR = 2, NC_SHORT = 3, NC_INT = 4, NC_FLOAT = 5, NC_DOUBLE = 6;

// Every integer field is promoted to NC_INT rather than mapped to its narrowest
// NetCDF type. NC_BYTE is signed in NetCDF-3 with no portable way to mark it
// unsigned, and byte/short variables need trailing padding when the element
// count is odd — which some readers mishandle. Widening costs a few bytes on
// mask fields and removes both problems.
const TYPE_INFO = {
    int8:    { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v, false) },
    uint8:   { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v, false) },
    int16:   { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v, false) },
    uint16:  { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v, false) },
    int32:   { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v, false) },
    uint32:  { nc: NC_INT,    size: 4, write: (dv, o, v) => dv.setInt32(o, v | 0, false) },
    float32: { nc: NC_FLOAT,  size: 4, write: (dv, o, v) => dv.setFloat32(o, v, false) },
    float64: { nc: NC_DOUBLE, size: 8, write: (dv, o, v) => dv.setFloat64(o, v, false) },
};

const pad4 = (n) => (n + 3) & ~3;

class Writer {
    constructor() { this.parts = []; this.len = 0; }
    push(bytes) { this.parts.push(bytes); this.len += bytes.length; }
    u32(v) { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0, v, false); this.push(b); }
    u64(v) { const b = new Uint8Array(8); new DataView(b.buffer).setBigUint64(0, BigInt(v), false); this.push(b); }
    name(s) {
        const bytes = new TextEncoder().encode(s);
        this.u32(bytes.length);
        const padded = new Uint8Array(pad4(bytes.length));
        padded.set(bytes);
        this.push(padded);
    }
    concat() {
        const out = new Uint8Array(this.len);
        let p = 0;
        for (const part of this.parts) { out.set(part, p); p += part.length; }
        return out;
    }
}

function attrBytes(w, name, value) {
    w.name(name);
    if (typeof value === 'string') {
        const bytes = new TextEncoder().encode(value);
        w.u32(NC_CHAR);
        w.u32(bytes.length);
        const padded = new Uint8Array(pad4(bytes.length));
        padded.set(bytes);
        w.push(padded);
    } else {
        const arr = Array.isArray(value) ? value : [value];
        w.u32(NC_DOUBLE);
        w.u32(arr.length);
        const b = new Uint8Array(pad4(arr.length * 8));
        const dv = new DataView(b.buffer);
        arr.forEach((v, i) => dv.setFloat64(i * 8, v, false));
        w.push(b);
    }
}

function attrListBytes(w, attrs) {
    const entries = Object.entries(attrs || {}).filter(([, v]) => v !== undefined && v !== null);
    if (!entries.length) { w.u32(0); w.u32(0); return; }
    w.u32(0x0C);                 // NC_ATTRIBUTE
    w.u32(entries.length);
    for (const [k, v] of entries) attrBytes(w, k, v);
}

/**
 * @param {object} spec
 *   dimensions: { name: length, ... }
 *   variables:  [{ name, dtype, dimensions: [dimName...], data: TypedArray, attributes: {} }]
 *   attributes: global attributes
 * @returns {Uint8Array}
 */
export function writeNetCDF(spec) {
    const dimNames = Object.keys(spec.dimensions);
    const dimIndex = new Map(dimNames.map((n, i) => [n, i]));

    // Total data size decides CDF-1 vs CDF-2.
    let dataTotal = 0;
    const varSizes = spec.variables.map(v => {
        const info = TYPE_INFO[v.dtype];
        if (!info) throw new Error(`Unsupported dtype for NetCDF: ${v.dtype}`);
        const count = v.dimensions.reduce((a, d) => a * spec.dimensions[d], 1);
        const size = pad4(count * info.size);
        dataTotal += size;
        return { count, size, info };
    });
    // CDF-2 (64-bit offset) once the data no longer fits classic 32-bit offsets.
    // spec.version forces a version, which is how the 64-bit path gets tested
    // without generating a 2 GB file.
    const version = spec.version ?? (dataTotal > 0x7FFFFFF0 ? 2 : 1);
    const offsetBytes = version === 2 ? 8 : 4;

    const header = new Writer();
    header.push(new Uint8Array([0x43, 0x44, 0x46, version]));  // "CDF" + version
    header.u32(0);                                             // numrecs = 0

    // Dimensions
    if (dimNames.length) {
        header.u32(0x0A);            // NC_DIMENSION
        header.u32(dimNames.length);
        for (const n of dimNames) { header.name(n); header.u32(spec.dimensions[n]); }
    } else { header.u32(0); header.u32(0); }

    // Global attributes
    attrListBytes(header, spec.attributes);

    // Variables. Data offsets aren't known until the header length is known,
    // and the header length depends on the offset field width — which is fixed
    // per version, so a single sizing pass with placeholder offsets is exact.
    const varHeaderStart = header.len;
    const emitVars = (offsets) => {
        const w = new Writer();
        w.u32(0x0B);                 // NC_VARIABLE
        w.u32(spec.variables.length);
        spec.variables.forEach((v, i) => {
            w.name(v.name);
            w.u32(v.dimensions.length);
            for (const d of v.dimensions) {
                const idx = dimIndex.get(d);
                if (idx === undefined) throw new Error(`Variable ${v.name} uses undeclared dimension ${d}`);
                w.u32(idx);
            }
            attrListBytes(w, v.attributes);
            w.u32(varSizes[i].info.nc);
            w.u32(varSizes[i].size);
            if (version === 2) w.u64(offsets[i]); else w.u32(offsets[i]);
        });
        return w.concat();
    };

    const placeholder = emitVars(spec.variables.map(() => 0));
    const headerLen = varHeaderStart + placeholder.length;

    const offsets = [];
    let cursor = headerLen;
    for (const vs of varSizes) { offsets.push(cursor); cursor += vs.size; }

    const varBytes = emitVars(offsets);
    if (varBytes.length !== placeholder.length) {
        throw new Error('NetCDF header sizing mismatch');
    }
    header.push(varBytes);

    const out = new Uint8Array(headerLen + dataTotal);
    out.set(header.concat(), 0);

    // Write the data.
    //
    // NetCDF is big-endian, and the obvious implementation — one DataView
    // setter per value — costs a call per element. On a 2048x1024 grid with a
    // hundred variables that is a quarter of a billion calls and dominates the
    // whole export. Instead, fill a native-endian typed-array view of the
    // output buffer and byte-swap the region in bulk, which Node does in native
    // code. Falls back to the per-element path when alignment forbids a view.
    const dv = new DataView(out.buffer);
    spec.variables.forEach((v, i) => {
        const { info, count } = varSizes[i];
        const base = offsets[i];
        const data = v.data;

        if (info.size === 4 && base % 4 === 0) {
            const view = info.nc === NC_FLOAT
                ? new Float32Array(out.buffer, base, count)
                : new Int32Array(out.buffer, base, count);
            for (let k = 0; k < count; k++) view[k] = data[k];
            Buffer.from(out.buffer, base, count * 4).swap32();
            return;
        }
        if (info.size === 8 && base % 8 === 0) {
            const view = new Float64Array(out.buffer, base, count);
            for (let k = 0; k < count; k++) view[k] = data[k];
            Buffer.from(out.buffer, base, count * 8).swap64();
            return;
        }
        for (let k = 0; k < count; k++) info.write(dv, base + k * info.size, data[k]);
    });

    return out;
}

/** dtype that survives a NetCDF round-trip for a given source dtype. */
export function netcdfDtype(dtype) {
    return TYPE_INFO[dtype] ? dtype : 'float32';
}
