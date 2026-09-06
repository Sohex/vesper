#!/usr/bin/env python
"""Does the shallow-ice diffusion converge on this project's mesh, and what does it cost?

    python analysis/shallow_ice_route_a.py --regions 80000 --steps 20

WORLDBUILDING FRAME. Vesper is a simulated super-Earth. Everything below is a
property of this project's region mesh and of a numerical scheme on it; nothing
here is about the real world.

CLIM-62 asks for the cost of each route to a shallow-ice flow solver, and
declares the prior gate: if the doubly nonlinear diffusion will not converge on
this mesh then Route A is unavailable at any cost and the comparison is over.
This is that gate. The criteria it is judged against, and the date they were
fixed, are in `notes/audits/shallow-ice-solver-cost.md`.

## The identity, because a test needs a right answer

The Halfar (1981) similarity solution for isothermal shallow ice at zero mass
balance on a flat bed. For Glen exponent n = 3,

    H(r, t) = H0 (t0/t)^(1/9) [1 - ((t0/t)^(1/18) r/R0)^(4/3)]^(3/7)
    R(t)    = R0 (t/t0)^(1/18)
    t0      = (1/18) (7/4)^3 R0^4 / (Gamma H0^7),   Gamma = (2/5) A (rho g)^3

Three things about it can FAIL rather than merely differ. Volume is conserved
exactly, because the mass balance is zero. The margin radius is closed form. And
the profile itself is closed form everywhere. A dome that merely looks plausible
is not evidence of anything.

Halfar is planar and this mesh is spherical, so the comparison carries a
contamination of order `(R0/R_planet)^2 / 6`. It is reported, and the dome is
sized so that it sits far below GW-8's measured 12 per cent operator noise floor
rather than beside it.

## The scheme, and the one choice that is not free

Two-point flux finite volume on the mesh's own Voronoi dual, which is the
discretisation `hydrography/scripts/groundwater.py` already owns and whose face
widths and areas that module reconstructs exactly. The flux across the face
between cells i and j is `D_ij g_ij (H_j - H_i)` with `g_ij = face_width /
generator_separation`, and the step is implicit in H with D lagged, which is a
Picard iteration.

THE CHOICE IS THE SURFACE SLOPE. `D` depends on `|grad z_S|^(n-1)`, and a
two-point flux scheme sees only the NORMAL component of that gradient across
each face. The transverse component is invisible to it. Both are computed here:
`--slope normal` uses the two-point difference alone, `--slope full`
reconstructs a cell-centred gradient by least squares over each cell's
neighbours and averages its magnitude to the face. The two are reported side by
side because the difference between them is a cost of Route A rather than a
detail of it -- the normal-only form is what the existing operator gives for
free, and the full form is a per-cell least-squares solve this project does not
have.

## The free boundary

`H >= 0` is an obstacle, and the ice margin is where it binds. The active set
here is the same shape as the one `groundwater.py:solve` runs for `h <= z`, with
the sense reversed: cells driven negative by a solve are pinned at zero, the
reduced system is re-solved on the rest, and a pinned cell whose balance turns
positive is released. What it does NOT inherit is that solver's termination
argument, which rests on the matrix being fixed and SPD; here `D` depends on the
unknown, and the report says what that costs by measuring it rather than by
asserting it.

## What is measured

  - the Picard residual per outer iteration, and whether it falls monotonically
    to the declared tolerance or stalls the way GW-9's exponential
    transmissivity stalled;
  - the margin radius against `R0 (t/t0)^(1/18)`, in cells;
  - the profile's relative RMS against Halfar, against GW-8's 12 per cent floor;
  - volume conservation, which is an identity here;
  - wall clock per step and its scaling in the region count, with the host load
    it was taken under recorded beside it.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "hydrography" / "scripts"))

SECONDS_PER_YEAR = 3.15576e7

# Glen's law, from references/big-mitgcm/MITgcmIS.py, which declares
# `Aglen = 0.6e-15` per year with `nglen = 3`. Neither constant changes any
# verdict here: both enter only through Gamma, which sets the timescale t0 and
# cancels out of every relative error and every iteration count.
A_GLEN_PER_YEAR = 0.6e-15
N_GLEN = 3

# Declared before the run, and mirrored from the note.
PICARD_TOLERANCE = 1e-8      # relative L2 of the nonlinear residual
PICARD_MAX = 60
OPERATOR_NOISE_FLOOR = 0.12  # GW-8, relative RMS above l = 1


def host_load() -> dict:
    one, five, fifteen = os.getloadavg()
    return dict(load_1min=round(one, 2), load_5min=round(five, 2),
                load_15min=round(fifteen, 2))


def halfar(r, t, h0, r0, gamma):
    """The similarity solution, and the margin radius that goes with it."""
    t0 = (1.0 / 18.0) * (7.0 / 4.0) ** 3 * r0 ** 4 / (gamma * h0 ** 7)
    s = (t0 / t) ** (1.0 / 18.0)
    inner = 1.0 - np.minimum(1.0, (s * r / r0) ** (4.0 / 3.0))
    h = h0 * (t0 / t) ** (1.0 / 9.0) * np.maximum(0.0, inner) ** (3.0 / 7.0)
    return h, t0, r0 / s


def cell_gradient_operator(geom, points):
    """Least-squares tangent-plane gradient at each cell, from its neighbours.

    Returns the two per-face coefficient arrays needed to evaluate
    `|grad H|` at a cell: for cell i the gradient is `sum_j c_ij (H_j - H_i)`
    with `c_ij` a tangent-plane vector, precomputed here because it does not
    move with H. This is the piece a two-point flux scheme does not have and
    the reason `--slope full` is a separate arm.
    """
    n = points.shape[0]
    src, dst = geom.src, geom.dst
    deg = np.bincount(np.concatenate([src, dst]), minlength=n)
    order = np.argsort(np.concatenate([src, dst]), kind="stable")
    nbr = np.concatenate([dst, src])[order]
    own = np.concatenate([src, dst])[order]
    offs = np.concatenate([[0], np.cumsum(deg)])

    # A tangent frame per cell.
    p = points / np.linalg.norm(points, axis=1, keepdims=True)
    ref = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    bad = np.abs(p[:, 2]) > 0.9
    ref[bad] = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(p, ref)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(p, e1)

    radius = geom.radius_m
    d = (p[nbr] - p[own]) * radius
    u = np.einsum("ij,ij->i", d, e1[own])
    v = np.einsum("ij,ij->i", d, e2[own])

    # Normal equations per cell, 2x2, assembled by segment sums.
    suu = np.add.reduceat(u * u, offs[:-1])
    svv = np.add.reduceat(v * v, offs[:-1])
    suv = np.add.reduceat(u * v, offs[:-1])
    det = suu * svv - suv * suv
    det = np.where(np.abs(det) < 1e-30, 1e-30, det)
    return dict(nbr=nbr, own=own, offs=offs, u=u, v=v,
                inv=(svv / det, -suv / det, suu / det), n=n)


def grad_magnitude(op, h):
    """`|grad H|` at each cell, from the precomputed least-squares operator."""
    dh = h[op["nbr"]] - h[op["own"]]
    su = np.add.reduceat(op["u"] * dh, op["offs"][:-1])
    sv = np.add.reduceat(op["v"] * dh, op["offs"][:-1])
    a, b, c = op["inv"]
    gu = a * su + b * sv
    gv = b * su + c * sv
    return np.sqrt(gu * gu + gv * gv)


def face_diffusivity(h, geom, gamma, slope_mode, gradop):
    """`D = Gamma H_face^(n+2) |grad z_S|^(n-1)`, at each face. Flat bed."""
    src, dst = geom.src, geom.dst
    h_face = 0.5 * (h[src] + h[dst])
    if slope_mode == "normal":
        slope = np.abs(h[dst] - h[src]) / geom.length_m
    else:
        gm = grad_magnitude(gradop, h)
        slope = 0.5 * (gm[src] + gm[dst])
    return gamma * h_face ** (N_GLEN + 2) * slope ** (N_GLEN - 1)


def assemble(d_face, geom, area, dt, n):
    """`(A/dt) - div(D grad .)`, a symmetric M-matrix for D >= 0."""
    w = d_face * geom.geom
    src, dst = geom.src, geom.dst
    rows = np.concatenate([src, dst, src, dst])
    cols = np.concatenate([dst, src, src, dst])
    vals = np.concatenate([-w, -w, w, w])
    m = csr_matrix((vals, (rows, cols)), shape=(n, n))
    return m + diags(area / dt)


def residual(h, h_old, geom, area, dt, gamma, slope_mode, gradop):
    """The NONLINEAR residual, which is what convergence is declared on.

    GW-9's exponential transmissivity held a residual flat while its head step
    fell thirtyfold, so the head step is not the bar and is not reported as one.
    """
    d = face_diffusivity(h, geom, gamma, slope_mode, gradop)
    w = d * geom.geom
    flux = w * (h[geom.dst] - h[geom.src])
    div = np.zeros_like(h)
    np.add.at(div, geom.src, flux)
    np.add.at(div, geom.dst, -flux)
    return area * (h - h_old) / dt - div


def step(h_old, geom, area, dt, gamma, slope_mode, gradop, tol, itmax,
         damping=1.0):
    """One implicit step: Picard on D, with an active set on `H >= 0`.

    `damping` under-relaxes the update, `H <- (1-w) H + w H_new`. Undamped is
    `damping = 1`. It is a parameter rather than a fixed choice because how much
    damping a step needs, and what that costs in linear solves, IS the price of
    Route A: each Picard pass is one sparse factorisation of the whole mesh.
    """
    n = h_old.size
    h = h_old.copy()
    trace = []
    pinned = np.zeros(n, dtype=bool)
    solves = 0
    for _ in range(itmax):
        d = face_diffusivity(h, geom, gamma, slope_mode, gradop)
        m = assemble(d, geom, area, dt, n)
        rhs = area / dt * h_old
        free = ~pinned
        if free.sum() == 0:
            break
        sub = m[free][:, free].tocsc()
        new = np.zeros(n)
        new[free] = spsolve(sub, rhs[free])
        solves += 1
        if damping < 1.0:
            new = (1.0 - damping) * h + damping * new
        neg = new < 0.0
        if neg.any():
            pinned |= neg
            new[neg] = 0.0
        h = new
        r = residual(h, h_old, geom, area, dt, gamma, slope_mode, gradop)
        scale = max(1e-30, float(np.linalg.norm(area * h_old / dt)))
        rel = float(np.linalg.norm(r)) / scale
        trace.append(rel)
        rel_set = pinned & (r < 0.0)
        if rel_set.any():
            pinned &= ~rel_set
            continue
        if rel < tol:
            return h, trace, True, solves
    return h, trace, False, solves


def check_gradient(n_regions, cfg, seed=0):
    """The least-squares gradient against an identity, before it is quoted.

    On a sphere of radius R the field `H = R cos(theta)` has tangential gradient
    magnitude exactly `sin(theta)`, everywhere except the poles where the
    identity is degenerate. That is a right answer the reconstruction cannot
    see, which is what makes it a test rather than a comparison. The two-point
    operator's own noise floor is GW-8's 12 per cent, so a reconstruction worse
    than that would make the `full` arm meaningless.
    """
    export, geom, pts, radius_km, _ = build_mesh(n_regions, cfg, seed)
    op = cell_gradient_operator(geom, pts)
    z = pts[:, 2]
    got = grad_magnitude(op, radius_km * 1000.0 * z)
    exact = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    away = exact > 0.2
    rel = np.abs(got[away] - exact[away]) / exact[away]
    return dict(field="H = R cos(theta), |grad H| = sin(theta)",
                cells=int(away.sum()),
                median_relative_error=float(np.median(rel)),
                p90_relative_error=float(np.percentile(rel, 90)),
                max_relative_error=float(rel.max()),
                against=OPERATOR_NOISE_FLOOR)


def build_mesh(n_regions, cfg, seed=0):
    import groundwater as gw
    import orogen as og
    radius_km = float(cfg["planet"]["radius_earth"]) * 6371.0
    jitter = 0.0
    try:
        reg = og.REGISTRY[cfg["source_build"]] if hasattr(og, "REGISTRY") else None
        jitter = float(reg.get("jitter", 0.0)) if reg else 0.0
    except Exception:
        jitter = 0.0
    if jitter <= 0.0:
        jitter = 0.5   # sphere-mesh.js's own default, and the value the
                       # refinement sweep in groundwater.py runs at
    export = gw.orogen_mesh(n_regions, jitter=jitter, radius_km=radius_km,
                            seed=seed)
    geom = gw.Geometry(export)
    pts = np.stack([export.x, export.y, export.z], axis=1).astype(np.float64)
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    return export, geom, pts, radius_km, jitter


def run(n_regions, steps, h0_m, r0_km, slope_mode, cfg, dt_factor=20.0,
        damping=1.0, seed=0):
    export, geom, pts, radius_km, jitter = build_mesh(n_regions, cfg, seed)
    n = export.n_regions
    area = geom.flux_area_m2
    radius_m = radius_km * 1000.0

    ga = float(cfg["planet"]["gravity_m_s2"])
    rho = float(cfg["surface"]["glacial_ice"]["density_kg_m3"]) \
        if "glacial_ice" in cfg.get("surface", {}) else 917.0
    a_glen = A_GLEN_PER_YEAR / SECONDS_PER_YEAR
    gamma = 2.0 / (N_GLEN + 2) * a_glen * (rho * ga) ** N_GLEN

    # Great-circle distance from a pole to every cell; the dome is centred there.
    centre = np.array([0.0, 0.0, 1.0])
    cosang = np.clip(pts @ centre, -1.0, 1.0)
    r = radius_m * np.arccos(cosang)

    r0 = r0_km * 1000.0
    _, t0, _ = halfar(np.zeros(1), 1.0, h0_m, r0, gamma)
    t_start = t0
    h = halfar(r, t_start, h0_m, r0, gamma)[0]

    cell_m = math.sqrt(4.0 * math.pi * radius_m ** 2 / n)
    cells_across_radius = r0 / cell_m

    # The step is sized against the explicit diffusion limit so that the
    # implicit scheme is being asked something a scheme could plausibly want,
    # rather than a step so small the nonlinearity never bites.
    d0 = face_diffusivity(h, geom, gamma, "normal", None)
    w = d0 * geom.geom
    rowsum = np.zeros(n)
    np.add.at(rowsum, geom.src, w)
    np.add.at(rowsum, geom.dst, w)
    dt_explicit = float(np.min(area[rowsum > 0] / rowsum[rowsum > 0]))
    dt = dt_factor * dt_explicit

    gradop = cell_gradient_operator(geom, pts) if slope_mode == "full" else None

    vol0 = float((area * h).sum())
    traces, converged, wall, nsolve = [], [], [], []
    t = t_start
    for _ in range(steps):
        t += dt
        c0 = time.perf_counter()
        h, tr, ok, ns = step(h, geom, area, dt, gamma, slope_mode, gradop,
                             PICARD_TOLERANCE, PICARD_MAX, damping)
        wall.append(time.perf_counter() - c0)
        traces.append(tr)
        converged.append(bool(ok))
        nsolve.append(ns)

    exact, _, margin_exact = halfar(r, t, h0_m, r0, gamma)
    inside = exact > 0.0
    rms = float(np.sqrt(((h[inside] - exact[inside]) ** 2).mean())
                / max(1e-30, np.sqrt((exact[inside] ** 2).mean())))
    have_ice = h > 1.0
    margin_model = float(r[have_ice].max()) if have_ice.any() else 0.0
    vol1 = float((area * h).sum())

    # Did the residual FALL, or stall the way GW-9's did?
    stalls = []
    for tr in traces:
        if len(tr) >= 10:
            stalls.append(float(tr[-1] / max(tr[len(tr) // 2], 1e-300)))
    return dict(
        regions=n, jitter=jitter, slope=slope_mode, dt_factor=dt_factor,
        damping=damping,
        mean_cell_km=round(cell_m / 1000.0, 3),
        cells_across_dome_radius=round(cells_across_radius, 2),
        sphericity_contamination=round((r0 / radius_m) ** 2 / 6.0, 6),
        dt_explicit_years=round(dt_explicit / SECONDS_PER_YEAR, 4),
        dt_years=round(dt / SECONDS_PER_YEAR, 4),
        steps=steps,
        every_step_converged=bool(all(converged)),
        steps_converged=int(sum(converged)),
        picard_iterations_median=float(np.median([len(t_) for t_ in traces])),
        picard_iterations_max=int(max(len(t_) for t_ in traces)),
        final_relative_residual_max=float(max(t_[-1] for t_ in traces)),
        residual_second_half_ratio_max=(max(stalls) if stalls else None),
        margin_radius_km_model=round(margin_model / 1000.0, 2),
        margin_radius_km_exact=round(margin_exact / 1000.0, 2),
        margin_error_cells=round((margin_model - margin_exact) / cell_m, 3),
        profile_relative_rms=round(rms, 5),
        operator_noise_floor=OPERATOR_NOISE_FLOOR,
        volume_relative_change=round((vol1 - vol0) / vol0, 8),
        linear_solves_per_step_median=float(np.median(nsolve)),
        residual_first_trace=[float(f"{x:.4g}") for x in traces[0][:12]],
        wall_s_per_step_median=round(float(np.median(wall)), 4),
        wall_s_per_linear_solve=round(float(np.median(wall))
                                      / max(1.0, float(np.median(nsolve))), 4),
        wall_s_total=round(float(sum(wall)), 3),
    )


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", type=int, nargs="+", default=[80000])
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--dome-height-m", type=float, default=3000.0)
    ap.add_argument("--dome-radius-km", type=float, default=750.0)
    ap.add_argument("--dt-factor", type=float, nargs="+", default=[20.0],
                    help="implicit step as a multiple of the explicit diffusion limit")
    ap.add_argument("--damping", type=float, nargs="+", default=[1.0],
                    help="Picard under-relaxation weight; 1.0 is undamped")
    ap.add_argument("--slope", nargs="+", default=["normal", "full"],
                    choices=["normal", "full"])
    ap.add_argument("--out", type=Path,
                    default=ROOT / "analysis/shallow_ice_route_a.json")
    args = ap.parse_args()

    cfg = yaml.safe_load((ROOT / "config/planet.yaml").read_text())
    load_before = host_load()
    gradient_check = check_gradient(min(args.regions), cfg) \
        if "full" in args.slope else None
    arms = []
    for n in args.regions:
        for mode in args.slope:
            for f in args.dt_factor:
                for w in args.damping:
                    arms.append(run(n, args.steps, args.dome_height_m,
                                    args.dome_radius_km, mode, cfg,
                                    dt_factor=f, damping=w))
                    print(json.dumps(arms[-1]))
    report = dict(
        measured_on="2026-09-05",
        host_load_before=load_before,
        host_load_after=host_load(),
        dome_height_m=args.dome_height_m,
        dome_radius_km=args.dome_radius_km,
        picard_tolerance=PICARD_TOLERANCE,
        picard_max=PICARD_MAX,
        glen_exponent=N_GLEN,
        glen_rate_factor_per_year=A_GLEN_PER_YEAR,
        gradient_reconstruction_check=gradient_check,
        arms=arms,
    )
    args.out.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
