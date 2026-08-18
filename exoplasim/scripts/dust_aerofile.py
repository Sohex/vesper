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

## What `apart` is set to here, and what that costs

The effective radius of the Balkanski distribution `dust_optics.py` integrates
over: `r_eff = r_mod * exp(2.5 ln^2 sigma)`, which for r_mod = 0.295 um and
sigma = 2.0 is 0.98 um. That is the optically weighted radius, so the optical
depth is right by construction.

It is NOT the mass-weighted radius. The emitted distribution has a volume-median
diameter of 3.4 um, so a single mode at 0.98 um settles far too slowly for the
coarse fraction and this file therefore overstates how long the mass stays up.
That is DUST-8, and it is the argument for per-bin sizes rather than one
effective radius; this file is what a single mode can do, with the cost stated.

## Two upstream defects this file depends on

`radmod`'s own `apart` is never populated from the namelist -- `aero_ini`
use-associates only `l_aerorad` and `aerofile` -- so it stays at its 50e-9 haze
default while `aerocore` uses the namelist value. Until that is patched, the
optical depth is wrong by `(50e-9 / apart)**2`, about 1/400 here. And
`aerocore.f90:1156` computes sphere volume with an integer `(4/3)`, so `nrho` is
33% high. Both are reported upstream; this file assumes both are fixed, and the
`--check` output says what it would be if they are not.
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
from paths import rel  # noqa: E402

OPTICS = ROOT / "analysis" / "dust_optics.json"
DUST_CFG = ROOT / "aeolian" / "config" / "dust.yaml"
OUT = ROOT / "exoplasim" / "data" / "dust" / "vesper_dust_aerosol.dat"

# The lognormal `dust_optics.py` integrates over. Kept here rather than imported
# so this file states what its own numbers are valid for.
R_MOD_UM, SIGMA_G, RHO_P_KG_M3 = 0.295, 2.0, 2600.0

# radmod's un-populated default, for the `--check` warning only.
RADMOD_DEFAULT_APART_M = 50e-9


def effective_radius_m(r_mod_um: float, sigma_g: float) -> float:
    """Area-weighted mean radius of a lognormal, <r^3>/<r^2>."""
    return r_mod_um * 1e-6 * np.exp(2.5 * np.log(sigma_g) ** 2)


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

    band1 = next(r for r in optics["results"]
                 if r["indices"] == cfg["optics"]["band1_indices"]
                 and r["band"] == "band 1")
    band2 = next(r for r in optics["results"]
                 if r["indices"] == cfg["optics"]["band2_indices"]
                 and r["band"] == "band 2")
    f1 = float(optics["stellar_flux_fraction_band1"])

    apart = effective_radius_m(R_MOD_UM, SIGMA_G)
    rows = []
    for band in (band1, band2):
        mee = band["mass_extinction_efficiency_m2_g"] * 1000.0     # m2/g -> m2/kg
        qext = qext_from_mee(mee, apart, RHO_P_KG_M3)
        qsca = qext * band["single_scattering_albedo"]
        qback = qsca * band["backscatter_fraction"]
        rows.append({"band": band["band"], "indices": band["indices"],
                     "wavelength_um": band["wavelength_um"],
                     "mee_m2_kg": mee, "qext": qext, "qsca": qsca,
                     "qback": qback, "g": band["asymmetry_parameter"]})

    print(f"apart = effective radius of the optics distribution = {apart * 1e6:.4f} um")
    print(f"rhop  = {RHO_P_KG_M3:.0f} kg/m3\n")
    for r in rows:
        print(f"{r['band']}  {r['wavelength_um'][0]}-{r['wavelength_um'][1]} um  "
              f"({r['indices']})")
        print(f"   MEE {r['mee_m2_kg']:7.2f} m2/kg -> Qext {r['qext']:.4f}   "
              f"Qsca {r['qsca']:.4f}   Qback {r['qback']:.4f}   g {r['g']:.4f}")

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

    penalty = (RADMOD_DEFAULT_APART_M / apart) ** 2
    print(f"\nIF the upstream `apart` defect is unpatched, radmod uses 50 nm and the "
          f"optical depth comes out {penalty:.5f}x this, i.e. 1/{1 / penalty:.0f}.")

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
        "apart_note": "Effective (area-weighted) radius of the optics lognormal, "
                      "r_mod exp(2.5 ln^2 sigma). Optically right and "
                      "mass-weighted WRONG: the emitted volume-median diameter is "
                      "3.4 um, so a single mode here settles too slowly. DUST-8.",
        "source_optics": rel(args.optics),
        "stellar_flux_fraction_band1": f1,
        "bands": rows,
        "effective_mass_extinction_efficiency_m2_kg": mee_back,
        "upstream_defects_assumed_fixed": [
            "radmod's apart is never populated from the namelist and stays 50e-9, "
            f"which would make optical depth {penalty:.5f}x this file's intent.",
            "aerocore.f90:1156 computes sphere volume with an integer (4/3), so "
            "number density is 33% high and optical depth with it.",
        ],
    }, indent=1) + "\n", encoding="utf-8")

    print(f"\nwrote {rel(args.output)}")
    print(f"      {rel(sidecar)}")


if __name__ == "__main__":
    main()
