#!/usr/bin/env python3
"""Fold per-rank perf samples into a grid<->spectral share, with its spread.

    python exoplasim/scripts/score_transform_profile.py <perf_dir> [<perf_dir> ...]

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
desktop. Nothing here is about the simulated planet.

Reads `perf.rank??.data` written by `profile_transforms.sh`, buckets symbols by
the table below, and reports each bucket's share of samples as a median over
ranks with the rank-to-rank spread beside it. The spread is not decoration: a
bulk-synchronous model has every rank doing the same transform work, so ranks
that disagree by much mean the profile caught something other than the
transform.

THE BUCKETS ARE DECLARED HERE, AHEAD OF THE MEASUREMENT, and the two boundaries
that decide the answer are drawn deliberately:

  `spectral_step` is NOT part of the transform. `spectrala`, `spectrald`,
  `makebm`, `minvers` and `hdiffo` are arithmetic on spectral coefficients --
  the semi-implicit solve and the hyperdiffusion. A different transform library
  leaves every one of them exactly where it is, so counting them would inflate
  the addressable share.

  `transform_mpi` IS part of it. `mpsumsc` and its siblings exist only because
  the Legendre sum is split across ranks; they are the decomposition's cost, and
  any scheme that changes who owns the latitude axis changes them too. Counting
  them separately is what keeps a transform speedup from being claimed over
  collectives it does not touch.

Unbucketed symbols are reported rather than swallowed, heaviest first, so a
bucket that is wrong shows up as a fat `unbucketed` line instead of a quiet
under-count.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics as st
import subprocess
from collections import defaultdict
from pathlib import Path

# Symbol -> bucket. Exact matches on the gfortran symbol (trailing underscore),
# taken from `nm` on most_plasim_t42_l10_p16.x.
BUCKETS: dict[str, tuple[str, ...]] = {
    # The Legendre transform proper, including the four routines that fuse
    # physics into it (mktend, qtend, dv2uv, uv2dv).
    "legendre": (
        "fc2sp_", "sp2fc_", "sp2fcdmu_", "sp2fl_",
        "dv2uv_", "uv2dv_", "mktend_", "qtend_", "invlega_", "invlegd_",
    ),
    # The longitudinal FFT. fftini_ builds the trig table and is startup, but it
    # reallocates on a resolution change, so it sits here where it is visible.
    "fft": (
        "gp2fc_", "fc2gp_", "fftini_", "fftini_.part.0",
        "dfft2_", "dfft3_", "dfft4_", "dfft8_",
        "ifft2_", "ifft3_", "ifft4_", "ifft8_",
    ),
    # The Fortran collective wrappers. Their own time is marshalling; the wait
    # is inside the MPI library and is bucketed by DSO below.
    "mpi": (
        "mpsumsc_", "mpsum_", "mpsumr_", "mpsumbcr_",
        "mpgallsp_", "mpgasp_", "mpscsp_", "mpgacs_",
        "mpgagp_", "mpscgp_", "mpbcr_", "mpbci_",
    ),
    # Accumulated every timestep whether or not anything is written
    # (plasim.f90:688 is unguarded), so it is per-step model work that a
    # production run pays too -- NOT part of the output path it looks like.
    "diagnostics": ("outaccu_", "outreset_", "outdiag_", "diag_", "energy_", "entropy_s_"),
    # The rest of the column physics. Present so that `unbucketed` means
    # "a symbol the table does not know" rather than "physics".
    "physics": (
        "fluxstep_", "rainstep_", "surfstep_", "miscstep_", "icestep_",
        "oceanstep_", "landstep_", "seastep_", "vegstep_", "soilstep_",
        "glacierstep_", "glacierprep_", "carbonstep_", "hurricanestep_",
        "roffstep_", "clsgstep_", "lsg_", "planet_step_",
        "kuo_", "mkrain_", "mklsp_", "mkdca_", "mkshallow_", "mkcflux_",
        "vdiff_", "vdiffo_", "mkevap_", "mkdqtgp_", "surflx_", "getflx_",
        "getflxco_", "mkiflux_", "mkiflx_", "mkstress_", "mkshfl_",
        "mksst_", "mkice_", "mkicec_", "mktsoil_", "skintemp_",
        "iceget_", "oceanget_", "gsettle_", "subsnow_", "newtonraphson_",
        "cape_", "stability_check_", "getshear_", "ql_", "qld_", "lv_", "rv_",
        "ev_", "es_cc_", "e_plcl_", "density_", "trho_", "gettcll_", "poti_",
        "tands_", "momint_", "vdiffo_", "mknudge_", "mkflcor_", "mkflukoi_",
    ),
    # THE SHTns PATH, declared before the first profile was read. It is the
    # counterpart of `legendre` AND `fft` together, not of `legendre` alone:
    # SHTns does the Legendre sum and the longitudinal FFT in one call, so a
    # comparison that set it against the Legendre bucket by itself would
    # understate what it replaced.
    #
    # Split in two so the model's own marshalling is visible separately from
    # the library's arithmetic. `shtns_wrap` is the conversion, the packing and
    # the full-globe temporaries -- the part this project can still change --
    # and `shtns_lib` is SHTns itself. If the wrappers are a large share of the
    # pair, the temporaries are worth removing; if they are noise, they are not.
    "shtns_wrap": (
        "__shtnsmod_MOD_sh_sp2gp", "__shtnsmod_MOD_sh_dv2uv",
        "__shtnsmod_MOD_sh_sp2grad", "__shtnsmod_MOD_sh_gp2sp",
        "__shtnsmod_MOD_sh_uv2dv", "__shtnsmod_MOD_sh_advtend",
        "__shtnsmod_MOD_sh_dztend", "__shtnsmod_MOD_sh_slice",
        "__shtnsmod_MOD_analyse_uv", "__shtnsmod_MOD_unpack_sp",
        "__shtnsmod_MOD_unpack_dv", "__shtnsmod_MOD_shtns_setup",
    ),
    # Spectral-space arithmetic. Untouched by any transform swap; see above.
    "spectral_step": ("spectrala_", "spectrald_", "makebm_", "minvers_", "hdiffo_"),
    "radiation": (
        "radstep_", "swr_", "lwr_", "mko3_", "mkclouds_", "mkradv_",
        "solang_", "getalb_", "mkdheat_", "gen_orb_decl_", "orb_decl_",
    ),
    "dyn_grid": (
        "calcgp_", "gridpointa_", "gridpointd_", "prepare_uvps_",
        "gpot_", "energy_", "absvorticity_", "addfc_", "addfci_",
    ),
    "aerosol_tracer": (
        "aerocore_", "aero_main_", "aeroprof_", "aero_surf_", "dustprof_",
        "dustsrc_", "tracer_main_", "tpcore_", "fct3d_", "xtp_", "ytp_",
        "fxppm_", "fzppm_", "xadv_", "xmist_", "ymist_", "qckxyz_",
    ),
    "io": (
        "outgp_", "outsp_", "outsc_",
        "writegp_", "writesp_", "writescalar_", "write_short_",
        "snapshotgp_", "snapshotsp_", "snapshotsc_", "snapshotdiag_",
        "hcadencegp_", "hcadencesp_", "hcadencesc_", "hcadencediag_",
        "mpwritegph_",
    ),
    "startup": (
        "legini_", "inigau_", "inilat_", "inilat_.part.0", "initpm_", "initsi_",
        "initfd_", "prolog_", "readnl_", "readdat_", "readarray_",
        "restart_ini_", "read_atmos_restart_", "read_ice_surface_",
        "surface_ini_", "get_surf_array_", "get_restart_array_",
        "get_restart_integer_", "get_restart_real_", "get_restart_seed_",
        "oroini_", "radini_", "rainini_", "fluxini_", "surfini_", "seaini_",
        "iceini_", "oceanini_", "landini_", "soilini_", "vegini_", "miscini_",
        "glacierini_", "carbonini_", "hurricaneini_", "solarini_", "clsgini_",
        "aero_ini_", "tracer_ini_", "tracer_ini0_", "outini_", "calini_",
        "planet_ini_", "print_planet_", "mpstart_", "allocate_arrays_",
    ),
}

# DSO -> bucket, applied when the symbol is unknown. MPI library internals are
# the other half of transform_mpi: mpsumsc is a thin wrapper and most of the
# time lands inside Open MPI.
DSO_BUCKETS = (
    (re.compile(r"libmpi|libopen-pal|libopen-rte|libpmix|mca_"), "mpi"),
    (re.compile(r"\[kernel|\[unknown"), "kernel"),
    # The spin loop reads the clock through the vDSO, so this is MPI wait, not
    # kernel work. Bucketed before the kernel rule for that reason.
    (re.compile(r"\[vdso"), "mpi"),
)

# Symbols that are unambiguously Open MPI WAITING rather than moving bytes.
# Reported as a subtotal of `mpi` so the wait can be named as such.
SPIN_SYMBOLS = re.compile(r"opal_progress|mca_btl_sm_poll|__vdso_gettimeofday|"
                          r"opal_timer|sched_yield|poll_handle_frag")

SYMBOL_BUCKET = {sym: bucket for bucket, syms in BUCKETS.items() for sym in syms}

# SHTns's own arithmetic, matched by PREFIX because the library names one symbol
# per specialised kernel -- SHsphtor_to_spat_fly2_m0l and dozens like it -- and
# an exact list would go stale on a library upgrade without anyone noticing.
# Declared before the first profile was read, with the wrapper bucket above.
#
# The `_fly` in those names is worth knowing: it is SHTns's ON-THE-FLY
# algorithm, which recomputes the Legendre functions with SIMD instead of
# streaming stored tables. shtns_setup asks for SHT_QUICK_INIT, which picks by a
# fixed heuristic rather than by timing, because the timing made the model
# irreproducible. If these dominate, that choice is worth revisiting against a
# deterministic way of getting the stored-table path.
SYMBOL_PREFIX_BUCKETS = (
    ("SH_to_spat", "shtns_lib"),
    ("spat_to_SH", "shtns_lib"),
    ("SHsph", "shtns_lib"),
    ("SHtor", "shtns_lib"),
    ("SHqst", "shtns_lib"),
    ("shtns_", "shtns_lib"),
    ("fftw", "shtns_lib"),
)

# What the profile is being read FOR. The addressable share is what a different
# transform could touch at all; everything outside it bounds the answer by
# Amdahl no matter how fast the transform becomes.
#
# `mpi` is NOT in it, and the first arm is why. Open MPI polls a shared-memory
# queue while a rank waits at a collective, so a rank that arrives early burns
# CPU in `opal_progress` and perf counts it as work. Measured on arm 1, the
# same bucket ran 12.5% on the busiest rank and 44.5% on the idlest -- ranks
# doing identical transform work, so the variation is waiting, not transform.
# Waiting is slack: making the Legendre transform faster shortens the critical
# path and lets the idle ranks wait longer, it does not save their spin.
ADDRESSABLE = ("legendre", "fft")
# Excluded from the compute normalisation entirely: neither is model arithmetic.
NOT_COMPUTE = ("mpi", "kernel")

# A bulk-synchronous model gives every rank the same transform work. Ranks that
# disagree by more than this on the addressable share mean the profile is
# reading something else -- an unquiet machine, or a rank doing I/O for the
# others -- and the number is not to be quoted.
RANK_SPREAD_FLOOR_PCT = 20.0
# Buckets that should be empty in a bed built from run_4182235e9781: L_AERO=0
# and every output stream is off. A fat one means the bed is wrong, not that the
# model is slow.
EXPECT_EMPTY = ("aerosol_tracer", "io")
EXPECT_EMPTY_PCT = 2.0


def lost_samples(data: Path) -> int:
    """Samples perf dropped. A ring buffer that overflowed is a biased profile:
    what is lost is whatever ran while perf was behind, which is not random."""
    out = subprocess.run(["perf", "report", "--stdio", "-i", str(data)],
                         capture_output=True, text=True).stdout
    m = re.search(r"Total Lost Samples:\s*(\d+)", out)
    return int(m.group(1)) if m else 0


def rank_shares(data: Path) -> dict[str, float]:
    """Bucket -> percent of this rank's samples."""
    out = subprocess.run(
        ["perf", "report", "--stdio", "--quiet", "--no-children",
         "--percent-limit", "0", "-F", "overhead,dso,symbol", "-i", str(data)],
        capture_output=True, text=True, check=True,
    ).stdout

    shares: dict[str, float] = defaultdict(float)
    unbucketed: dict[str, float] = defaultdict(float)
    # Per-symbol shares, kept so a bucket can be opened up. A bucket total says
    # "physics, 12.9%", which is not a thing anyone can act on; the symbols under
    # it are.
    symbols: dict[str, float] = {}
    sym_dso: dict[str, str] = {}
    sym_bucket: dict[str, str] = {}
    for line in out.splitlines():
        m = re.match(r"\s*(\d+\.\d+)%\s+(\S+)\s+\[[.k]\]\s+(.+?)\s*$", line)
        if not m:
            continue
        pct, dso, sym = float(m.group(1)), m.group(2), m.group(3)
        symbols[sym] = symbols.get(sym, 0.0) + pct
        sym_dso[sym] = dso
        bucket = SYMBOL_BUCKET.get(sym)
        if bucket is None:
            for prefix, b in SYMBOL_PREFIX_BUCKETS:
                if sym.startswith(prefix):
                    bucket = b
                    break
        if bucket is None:
            for pattern, b in DSO_BUCKETS:
                if pattern.search(dso):
                    bucket = b
                    break
        if bucket is None:
            bucket = "unbucketed"
            unbucketed[sym] += pct
        shares[bucket] += pct
        sym_bucket[sym] = bucket
        if bucket == "mpi" and SPIN_SYMBOLS.search(sym):
            shares["_mpi_spin_named"] += pct
    shares["_unbucketed_top"] = unbucketed  # type: ignore[assignment]
    shares["_symbols"] = symbols            # type: ignore[assignment]
    shares["_sym_bucket"] = sym_bucket      # type: ignore[assignment]
    shares["_sym_dso"] = sym_dso            # type: ignore[assignment]
    return shares


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("perf_dirs", type=Path, nargs="+")
    ap.add_argument("--per-symbol", type=int, default=0, metavar="N",
                    help="also list the N heaviest symbols in each bucket, "
                         "which is what makes a bucket actionable")
    ap.add_argument("--exclude-startup", action="store_true", default=True,
                    help="renormalise over the integration, dropping the "
                         "startup BUCKET (default: on; a 600-step bed carries "
                         "~1.6 s of it). Note the hole: it drops symbols the "
                         "table calls startup, and startup work that lands in "
                         "a stripped library falls in `unbucketed` instead and "
                         "survives. oroini_ and roffini_ between them zero 26 "
                         "GiB through libc before step one, which at 25 steps "
                         "is most of the libc share and at 300 is noise -- so "
                         "a short profile overstates it and this flag does not "
                         "save you. Profile the length you mean to run.")
    args = ap.parse_args()

    ranks: list[dict[str, float]] = []   # one normalised bucket table per rank
    unbucketed_total: dict[str, float] = defaultdict(float)
    symbol_total: dict[str, float] = defaultdict(float)
    symbol_bucket: dict[str, str] = {}
    symbol_dso: dict[str, str] = {}
    passes = []

    for d in args.perf_dirs:
        files = sorted(d.glob("perf.rank??.data"))
        if not files:
            print(f"  (no perf.rank??.data in {d}, skipped)")
            continue
        meta = d / "pass.json"
        if meta.is_file():
            passes.append(json.loads(meta.read_text()))
        for f in files:
            lost = lost_samples(f)
            if lost:
                print(f"  ** {f.name} lost {lost} samples: the ring buffer overflowed, "
                      f"and what it dropped is whatever ran while perf was behind. "
                      f"Re-record with a larger PERF_MMAP_PAGES; do not score this.")
            s = rank_shares(f)
            for sym, pct in s.pop("_symbols").items():      # type: ignore[union-attr]
                symbol_total[sym] += pct
            symbol_bucket.update(s.pop("_sym_bucket"))      # type: ignore[arg-type]
            symbol_dso.update(s.pop("_sym_dso"))            # type: ignore[arg-type]
            top = s.pop("_unbucketed_top")
            for sym, pct in top.items():  # type: ignore[union-attr]
                unbucketed_total[sym] += pct
            spin = s.pop("_mpi_spin_named", 0.0)
            if args.exclude_startup:
                s.pop("startup", None)
            total = sum(s.values())
            if not total:
                continue
            table = {b: 100.0 * p / total for b, p in s.items()}
            # Shares of COMPUTE, which is where the transform question lives:
            # MPI wait and kernel time are not model arithmetic and a faster
            # transform does not remove them.
            compute = sum(p for b, p in table.items() if b not in NOT_COMPUTE)
            table["_compute"] = compute
            table["_mpi_spin_named"] = 100.0 * spin / total
            ranks.append(table)

    if not ranks:
        raise SystemExit("no samples read")

    nranks = len(ranks)
    names = sorted({b for table in ranks for b in table if not b.startswith("_")})
    per_rank = {b: [table.get(b, 0.0) for table in ranks] for b in names}

    if passes:
        walls = [p["wall_s"] for p in passes]
        shas = {p.get("status_sha256_16") for p in passes}
        print(f"passes: {len(passes)}  wall {min(walls):.1f}-{max(walls):.1f} s"
              f"  restart sha {'AGREES' if len(shas) == 1 else 'DISAGREES ' + str(shas)}")
    print(f"ranks read: {nranks}"
          f"{'  (startup excluded, renormalised over the integration)' if args.exclude_startup else ''}")
    print()

    print(f"{'bucket':16} {'median %':>9} {'min':>7} {'max':>7}")
    order = sorted(per_rank, key=lambda b: -st.median(per_rank[b]))
    for bucket in order:
        v = per_rank[bucket]
        print(f"{bucket:16} {st.median(v):9.2f} {min(v):7.2f} {max(v):7.2f}")

    # Transform as a share of each rank's COMPUTE, not of its samples.
    share_of_compute = [
        100.0 * sum(t.get(b, 0.0) for b in ADDRESSABLE) / t["_compute"]
        for t in ranks if t.get("_compute")
    ]
    mpi = [t.get("mpi", 0.0) for t in ranks]
    spin = [t.get("_mpi_spin_named", 0.0) for t in ranks]

    print()
    print(f"MPI (wait + transfer):   {st.median(mpi):5.1f}% of samples "
          f"[{min(mpi):.1f}, {max(mpi):.1f}]; of that, {st.median(spin):.1f}% is "
          f"NAMED spin (opal_progress, sm poll, vdso clock)")
    print(f"  The spread IS the imbalance: identical ranks, so the rank at "
          f"{min(mpi):.1f}% is closest to the critical path and the one at "
          f"{max(mpi):.1f}% spent the difference waiting.")

    med = st.median(share_of_compute)
    spread = 100.0 * (max(share_of_compute) - min(share_of_compute)) / med if med else 0.0
    print()
    print(f"TRANSFORM SHARE OF COMPUTE (legendre + fft): {med:.2f}% "
          f"[{min(share_of_compute):.2f}, {max(share_of_compute):.2f}], "
          f"rank spread {spread:.1f}%")
    if not med:
        print("  ** zero. Either these samples are not of the model, or the buckets "
              "no longer match its symbols; check the unbucketed list below.")
    elif spread > RANK_SPREAD_FLOOR_PCT:
        print(f"  ** rank spread exceeds the {RANK_SPREAD_FLOOR_PCT}% floor declared "
              f"before the run: DO NOT QUOTE this number, find out what the ranks "
              f"were doing differently.")

    # Wall-clock ceiling. The critical path is the busiest rank: its compute plus
    # the transfer nobody can skip. Slack on the other ranks is free, so the
    # saving a faster transform buys is its share of the BUSIEST rank's compute.
    busiest = min(range(len(ranks)), key=lambda i: ranks[i].get("mpi", 0.0))
    b_compute = ranks[busiest]["_compute"]
    b_share = 100.0 * sum(ranks[busiest].get(b, 0.0) for b in ADDRESSABLE) / b_compute
    ceiling = b_compute * b_share / 100.0
    print(f"  Critical-path rank spends {b_compute:.1f}% of its samples computing, "
          f"{b_share:.1f}% of which is transform.")
    print(f"  CEILING ON WALL TIME: an infinitely fast transform saves {ceiling:.1f}%; "
          f"halving the transform saves {ceiling / 2:.1f}%.")
    print(f"  Neither figure touches the collectives: a faster Legendre transform "
          f"moves the same bytes through the same reduce_scatter.")

    for bucket in EXPECT_EMPTY:
        if bucket in per_rank and st.median(per_rank[bucket]) > EXPECT_EMPTY_PCT:
            print(f"  ** {bucket} is {st.median(per_rank[bucket]):.1f}%, over the "
                  f"{EXPECT_EMPTY_PCT}% declared for a bed with it switched off: the bed "
                  f"is misconfigured, not the model slow.")

    if args.per_symbol:
        # A bucket total is a category, not a target. This opens each one so the
        # answer to "where do the cycles go" is a routine rather than a noun.
        print(f"\nheaviest symbols in each bucket (share of all samples)")
        by_bucket: dict[str, list[tuple[str, float]]] = defaultdict(list)
        for sym, pct in symbol_total.items():
            by_bucket[symbol_bucket.get(sym, "unbucketed")].append((sym, pct))
        order = sorted(by_bucket, key=lambda b: -sum(p for _, p in by_bucket[b]))
        for b in order:
            entries = sorted(by_bucket[b], key=lambda kv: -kv[1])
            tot = sum(p for _, p in entries) / max(nranks, 1)
            print(f"\n  {b}  ({tot:.2f}%)")
            for sym, pct in entries[:args.per_symbol]:
                share = pct / max(nranks, 1)
                if share < 0.01:
                    break
                dso = symbol_dso.get(sym, "")
                where = "" if dso.startswith("probe") or not dso else f"   [{dso}]"
                print(f"    {share:6.2f}%  {sym}{where}")

    if unbucketed_total:
        print(f"\nheaviest unbucketed symbols (sum over ranks, uncalibrated):")
        for sym, pct in sorted(unbucketed_total.items(), key=lambda kv: -kv[1])[:15]:
            print(f"  {pct / max(nranks, 1):6.2f}%  {sym}")


if __name__ == "__main__":
    main()
