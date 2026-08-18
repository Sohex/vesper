#!/usr/bin/env python3
"""Write the `aerofile` ExoPlaSim's radiation reads, from this world's dust optics.

    python exoplasim/scripts/dust_aerofile.py

DUST-9. `radmod.f90` reads eight numbers through
`readdat(aerofile,1,8,aeroqs)`: one header line, then Qext, Qsca, Qback and g for
band 1, then the same four for band 2. This derives them from
`analysis/dust_optics.json`, which already band-averages Mie calculations on the
k25v spectrum across exactly ExoPlaSim's 0.75 um split.

Generated rather than hand-written, for the same reason `dust_optics.py` exists:
a number that lives only in a scratchpad is a number the repository cannot check.

## The conversion, which is where this goes wrong silently

`dust_optics.json` reports a MASS extinction efficiency in m2/g, integrated over
a lognormal size distribution. ExoPlaSim wants a dimensionless EXTINCTION
EFFICIENCY against the geometric cross-section of one sphere of radius `apart`,
because `radmod.f90:1844` builds optical depth as

    aod = nrho * PI * apart**2 * Qext * dz

with `nrho` coming from `mmr2n`, which converts a mass mixing ratio using the
same `apart` and `rhop`:

    nrho = mmr * rho_air / ((4/3) * PI * apart**3 * rhop)

Substituting one into the other, the geometry collapses:

    aod = column_mass * 3 * Qext / (4 * apart * rhop)

so ExoPlaSim's effective mass extinction efficiency is `3 Qext / (4 apart rhop)`
and the number this file must carry is

    Qext = 4 * apart * rhop * MEE / 3

**Only the ratio `Qext / apart` affects optical depth.** `apart` cancels between
the number density and the cross-section. It does NOT cancel in the
gravitational settling inside `aerocore`, which is the whole reason `apart` still
has to be chosen physically rather than for optical convenience.

## What `apart` is set to here, and why it is not the optical radius

That identity is what lets the two radii be decoupled, and DUST-8 decided they
should be. Setting `apart` from the SETTLING requirement and building `Qext` at
that same radius leaves the optical depth, the single-scattering albedo, the
backscatter ratio and the band-2 ratio all unchanged, because every one of them
is a ratio the rescale preserves. The LMD Generic PCM decouples the same pair
explicitly: `radius(iq)` drives sedimentation and `reffrad` drives the optics.

So `apart` here is the BURDEN-MATCHED radius, the one whose settling velocity is
`1 / sum(f_b / v_b)` over the Kok (2011) emitted distribution in
`aeolian/config/dust.yaml`, evaluated with `build_dust.settling_velocity` so
that the in-model mode and the offline chain agree on what settling means. The
optically weighted radius, `r_mod exp(2.5 ln^2 sigma)` of the distribution
`dust_optics.py` integrates over, keeps the dust up several times too long for
the burden, and that is a real error in the AOD the in-model chain exists to
produce.

**The Mie plausibility check therefore applies to the UNSCALED number.** `Qext`
at the optical radius is a dimensionless efficiency and can be sanity-checked
against Mie theory for a micron particle; `Qext` at the burden-matched radius is
larger by the ratio of the two radii and is NOT an efficiency at all. Both are
printed, and only the first means anything as a physical cross-section.

**What the single mode still costs is deposition, not optical depth.** In steady
state total deposition equals total emission whatever the settling velocity, so
the cost is WHERE it lands: one lifetime is one travel distance, and a single
burden-matched mode carries the coarse fraction too far and the sub-micron
fraction not far enough. The offline chain therefore keeps producing the
size-resolved deposition field that pedology and the phosphorus budget read, and
the in-model chain is not asked for it. DUST-8.

## The upstream defects this file depends on

Both are patched and neither is resident yet; they are in `PENDING_PATCHES` in
`exoplasim/scripts/rebuild_binaries.py` and land at the next rebuild.
`exoplasim-3.4.2-aerosol-apart.patch` gives `radmod` the namelist's `apart`,
which it never received, and without which the optical depth is
`(50e-9/apart)**2` of intent. `exoplasim-3.4.2-aerocore-defects.patch` fixes the
integer `(4/3)` in `mmr2n` that made `nrho` 33% high. This file assumes both,
and the run refuses to start without the first because `apart` is what the
aerofile is declared valid for.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "lib") not in sys.path:
    sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dust_indices import selection  # noqa: E402
from paths import rel  # noqa: E402


def _aeolian_build_dust():
    """Import `aeolian/scripts/build_dust.py`, `_paths` collision and all.

    The settling velocity and the emitted size distribution belong to the
    aeolian component, and the whole point of matching this mode to the burden
    is that the in-model particle and the offline chain's bins settle by the
    SAME formula. Copying it here would be a second statement of it.

    Every component keeps its own `_paths.py` and only one can hold that module
    name, so the aeolian directory goes to the front of `sys.path` for the
    duration of this import and whatever held `_paths` is put back after.
    """
    import importlib
    aeolian = str(ROOT / "aeolian" / "scripts")
    saved = sys.modules.pop("_paths", None)
    sys.path.insert(0, aeolian)
    try:
        return importlib.import_module("build_dust")
    finally:
        sys.path.remove(aeolian)
        sys.modules.pop("_paths", None)
        if saved is not None:
            sys.modules["_paths"] = saved


_build_dust = _aeolian_build_dust()
emitted_mass_fractions = _build_dust.emitted_mass_fractions
settling_velocity = _build_dust.settling_velocity

OPTICS = ROOT / "analysis" / "dust_optics.json"
PLANET_CFG = ROOT / "config" / "planet.yaml"
DUST_CFG = ROOT / "aeolian" / "config" / "dust.yaml"
OUT = ROOT / "exoplasim" / "data" / "dust" / "vesper_dust_aerosol.dat"

# The lognormal `dust_optics.py` integrates over. Kept here rather than imported
# so this file states what its own numbers are valid for.
R_MOD_UM, SIGMA_G, RHO_P_KG_M3 = 0.295, 2.0, 2600.0

# Reference state the two radii are COMPARED at. It is not a claim about the
# atmosphere: both the per-bin velocities and the inversion back to a radius use
# it, so it cancels to first order, and this world's gravity cancels exactly
# because Stokes velocity is linear in g on both sides.
REF_TEMPERATURE_K, REF_AIR_DENSITY_KG_M3 = 288.0, 1.15

# radmod's un-populated default, for the penalty report only.
RADMOD_DEFAULT_APART_M = 50e-9


def effective_radius_m(r_mod_um: float, sigma_g: float) -> float:
    """Area-weighted mean radius of a lognormal, <r^3>/<r^2>."""
    return r_mod_um * 1e-6 * np.exp(2.5 * np.log(sigma_g) ** 2)


def burden_matched_radius_m(cfg: dict, gravity: float) -> tuple[float, float]:
    """Radius whose settling velocity reproduces the emitted burden's.

    A single mode has one lifetime, and which lifetime it should have depends on
    what the mode is FOR. This one is for the optical depth and the
    precipitation response, so it is matched on burden: at fixed emission the
    steady-state burden goes as the harmonic mean `1 / sum(f_b / v_b)` of the
    per-bin settling velocities, and the radius that gives that velocity is the
    one mode that carries the right amount of mass. Matching on the mass-weighted
    velocity instead would give a much larger radius and a much smaller burden.

    Returns the radius in metres and the velocity it was matched to, in m/s.
    """
    frac, dbar_um = emitted_mass_fractions(cfg)
    temperature = np.array([REF_TEMPERATURE_K])
    density = np.array([REF_AIR_DENSITY_KG_M3])
    # One particle, one density. The aerofile's Qext and the settling that picks
    # its radius have to describe the same grain, so this refuses rather than
    # quietly matching a radius for a particle the optics do not describe.
    rho_p = float(cfg["removal"]["particle_density_kg_m3"])
    if abs(rho_p - RHO_P_KG_M3) > 1e-9:
        raise SystemExit(
            f"aeolian/config/dust.yaml has particle_density_kg_m3 = {rho_p} but "
            f"this file builds Qext at {RHO_P_KG_M3}. They are the same grain; "
            "fix one.")

    def velocity(diameter_um: float) -> float:
        return float(settling_velocity(diameter_um, RHO_P_KG_M3, density,
                                       temperature, gravity)[0])

    target = 1.0 / float(np.sum(frac / np.array([velocity(d) for d in dbar_um])))
    # Monotonic in diameter over any range this distribution reaches, so a
    # bisection is exact rather than a fit.
    lo, hi = 1e-3, 1e3
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if velocity(mid) < target:
            lo = mid
        else:
            hi = mid
    diameter_um = np.sqrt(lo * hi)
    return diameter_um * 0.5e-6, target


def qext_from_mee(mee_m2_kg: float, apart_m: float, rho_p: float) -> float:
    """Dimensionless extinction efficiency reproducing a mass efficiency."""
    return 4.0 * apart_m * rho_p * mee_m2_kg / 3.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--optics", type=Path, default=OPTICS)
    ap.add_argument("--output", type=Path, default=OUT)
    ap.add_argument("--dust-config", type=Path, default=DUST_CFG)
    args = ap.parse_args()

    import yaml
    optics = json.loads(args.optics.read_text(encoding="utf-8"))
    cfg = yaml.safe_load(args.dust_config.read_text(encoding="utf-8"))

    # aeolian/config/dust.yaml is the one file that decides which refractive
    # indices this world's dust has, and every consumer reads it from there.
    # DUST-12.
    sel = selection(cfg)
    band1 = next(r for r in optics["results"]
                 if r["indices"] == sel["band1"] and r["band"] == "band 1")
    band2 = next(r for r in optics["results"]
                 if r["indices"] == sel["band2"] and r["band"] == "band 2")
    f1 = float(optics["stellar_flux_fraction_band1"])

    planet = yaml.safe_load(PLANET_CFG.read_text(encoding="utf-8"))
    gravity = float(planet["planet"]["gravity_m_s2"])
    r_optical = effective_radius_m(R_MOD_UM, SIGMA_G)
    apart, v_burden = burden_matched_radius_m(cfg, gravity)
    scale = apart / r_optical
    rows = []
    for band in (band1, band2):
        mee = band["mass_extinction_efficiency_m2_g"] * 1000.0     # m2/g -> m2/kg
        qext = qext_from_mee(mee, apart, RHO_P_KG_M3)
        qsca = qext * band["single_scattering_albedo"]
        qback = qsca * band["backscatter_fraction"]
        rows.append({"band": band["band"], "indices": band["indices"],
                     "wavelength_um": band["wavelength_um"],
                     "mee_m2_kg": mee, "qext": qext, "qsca": qsca,
                     "qback": qback, "g": band["asymmetry_parameter"],
                     # The Mie-plausible number, at the OPTICAL radius. The Qext
                     # written to the file is this times apart/r_optical and is
                     # not an efficiency; see the module docstring.
                     "qext_at_optical_radius":
                         qext_from_mee(mee, r_optical, RHO_P_KG_M3)})

    print(f"apart = burden-matched settling radius       = {apart * 1e6:.4f} um")
    print(f"        optically weighted radius            = {r_optical * 1e6:.4f} um")
    print(f"        Q scaling apart/r_optical            = {scale:.4f}")
    print(f"        burden-matched settling velocity     = {v_burden * 100:.4f} cm/s")
    print(f"rhop  = {RHO_P_KG_M3:.0f} kg/m3\n")
    for r in rows:
        print(f"{r['band']}  {r['wavelength_um'][0]}-{r['wavelength_um'][1]} um  "
              f"({r['indices']})")
        print(f"   MEE {r['mee_m2_kg']:7.2f} m2/kg -> Qext {r['qext']:.4f}   "
              f"Qsca {r['qsca']:.4f}   Qback {r['qback']:.4f}   g {r['g']:.4f}")
        print(f"   Mie plausibility, Qext at the optical radius "
              f"{r['qext_at_optical_radius']:.4f}")

    # -- validation: round-trip the optical depth the offline chain reports -----
    mee_chain = 1000.0 * (f1 * band1["mass_extinction_efficiency_m2_g"]
                          + (1 - f1) * band2["mass_extinction_efficiency_m2_g"])
    mee_back = f1 * 3 * rows[0]["qext"] / (4 * apart * RHO_P_KG_M3) \
        + (1 - f1) * 3 * rows[1]["qext"] / (4 * apart * RHO_P_KG_M3)
    print(f"\nround trip: chain MEE {mee_chain:.2f} m2/kg, "
          f"reconstructed from Qext {mee_back:.2f} m2/kg, "
          f"error {abs(mee_back / mee_chain - 1) * 100:.4f}%")
    if abs(mee_back / mee_chain - 1) > 1e-6:
        raise SystemExit("the Qext conversion does not round-trip; do not use this file")

    dust = json.loads((ROOT / "aeolian" / "analysis"
                       / "dust_baseline.json").read_text(encoding="utf-8"))
    aod_land = dust["shelter_bracket"]["central"]["land_mean_aod"]
    column = aod_land / dust["mass_extinction_efficiency_m2_kg"]
    aod_check = column * mee_back
    print(f"            at the DUST-1 burden, column {column * 1e3:.3f} g/m2 "
          f"-> AOD {aod_check:.4f} against the chain's {aod_land:.4f}")

    # -- the AOD invariance DUST-8 rests on, checked rather than asserted ------
    #
    # Decoupling the optical radius from the settling radius is only free if the
    # optical depth is untouched by the swap. It is, because aod goes as
    # Qext/apart and both scale together -- but "because the algebra says so" is
    # exactly the kind of claim that survives a sign error, so it is evaluated.
    mee_optical = f1 * 3 * rows[0]["qext_at_optical_radius"] \
        / (4 * r_optical * RHO_P_KG_M3) \
        + (1 - f1) * 3 * rows[1]["qext_at_optical_radius"] \
        / (4 * r_optical * RHO_P_KG_M3)
    if abs(mee_back / mee_optical - 1) > 1e-9:
        raise SystemExit(
            "rescaling apart changed the mass extinction efficiency by "
            f"{abs(mee_back / mee_optical - 1) * 100:.6f}%; it must not change "
            "it at all. Do not use this file.")
    print(f"            invariance: the same MEE at {r_optical * 1e6:.4f} um and "
          f"at {apart * 1e6:.4f} um, so optical depth, SSA and backscatter "
          "are unmoved by the rescale")

    penalty = (RADMOD_DEFAULT_APART_M / apart) ** 2
    print(f"\nIF exoplasim-3.4.2-aerosol-apart.patch is not resident, radmod uses "
          f"50 nm and the optical depth comes out {penalty:.6f}x this, i.e. "
          f"1/{1 / penalty:.0f}. Check rebuild_binaries.py --verify.")

    header = ("# Mineral dust optical constants for Vesper, k25v-weighted. "
              "Qext Qsca Qback g, band 1 then band 2. "
              f"VALID ONLY FOR apart={apart:.4e} m, rhop={RHO_P_KG_M3:.0f} kg/m3. "
              "Generated by exoplasim/scripts/dust_aerofile.py; do not edit.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as fh:
        fh.write(header + "\n")
        for r in rows:
            for key in ("qext", "qsca", "qback", "g"):
                fh.write(f"{r[key]:.6e}\n")

    sidecar = args.output.with_suffix(".provenance.json")
    sidecar.write_text(json.dumps({
        "note": "Provenance for the aerofile ExoPlaSim's radmod.f90 reads. "
                "Qext is a dimensionless efficiency against the geometric cross "
                "section of one sphere of radius apart, NOT a mass efficiency; "
                "see exoplasim/scripts/dust_aerofile.py for the conversion.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "apart_m": apart, "rhop_kg_m3": RHO_P_KG_M3,
        "apart_note": "BURDEN-MATCHED settling radius: the single mode whose "
                      "Stokes velocity is the harmonic mass mean of the Kok "
                      "(2011) emitted bins. Chosen over the optically weighted "
                      "radius because the optics do not constrain apart at all "
                      "-- only Qext/apart enters the optical depth -- while "
                      "settling does. DUST-8.",
        "optical_radius_m": r_optical,
        "q_scale_from_optical_radius": scale,
        "burden_matched_settling_velocity_m_s": v_burden,
        "settling_reference_state": {
            "temperature_k": REF_TEMPERATURE_K,
            "air_density_kg_m3": REF_AIR_DENSITY_KG_M3,
            "gravity_m_s2": gravity,
            "note": "The state the two radii are compared at, not a claim about "
                    "the atmosphere. Gravity cancels exactly and the reference "
                    "state cancels to first order, because the per-bin "
                    "velocities and the inversion use the same one.",
        },
        # The namelist values live WITH the aerofile rather than in a config, so
        # a run cannot be given a radius that disagrees with the optics it is
        # applied to. That is the pattern build_surface_dust.py already uses for
        # DUSTQLW, and DUST-3 item 4 reads these when it writes aero_namelist.
        "namelist_values": {"APART": apart, "RHOP": RHO_P_KG_M3},
        "source_optics": rel(args.optics),
        "source_dust_config": rel(args.dust_config),
        "stellar_flux_fraction_band1": f1,
        "bands": rows,
        "effective_mass_extinction_efficiency_m2_kg": mee_back,
        "upstream_patches_assumed_applied": [
            "exoplasim-3.4.2-aerosol-apart.patch, without which radmod keeps its "
            f"50e-9 default and optical depth is {penalty:.6f}x this file's "
            "intent.",
            "exoplasim-3.4.2-aerocore-defects.patch, without which the integer "
            "(4/3) in mmr2n makes number density 33% high and optical depth "
            "with it.",
        ],
    }, indent=1) + "\n", encoding="utf-8")

    print(f"\nwrote {rel(args.output)}")
    print(f"      {rel(sidecar)}")


if __name__ == "__main__":
    main()
