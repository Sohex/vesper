#!/bin/bash
# Build the SOCRATES cost bench for CLIM-61, outside the tree.
#
# Worldbuilding frame: a COMPUTE instrument. It builds a candidate radiation
# code and a driver for it, so that one radiation call's cost per column can be
# compared with the Vesper climate model's own. Nothing it produces is about the
# simulated planet.
#
#     exoplasim/scripts/socrates_cost_bench.sh <workdir>
#
# then point the cost harness at what it built:
#
#     SOCRATES_BENCH=<workdir>/bench/socrates_bench \
#       scripts/lock_and_run -m "CLIM-61 radiation cost" \
#       python exoplasim/scripts/radiation_cost_per_column.py
#
# WHY OUTSIDE THE TREE. `references/socrates` is read-only reference, and in a
# worktree it is a directory symlink into the main checkout, so a build there
# writes its objects straight over the main checkout's -- the same hazard
# `docs/src/reference/environment.md` records for `vendor/cgenie`. The source is
# therefore copied out and built in a work directory the caller names. `data/`
# is linked rather than copied because it is 222 MB and is only read.
#
# THE FLAGS ARE THE MODEL'S, from `config/planet.yaml`'s production line, so
# what is compared is the code and not the compiler. `-ffpe-trap` is the one
# omission: SOCRATES is not written to run under it, and
# `notes/audits/aocl-and-model-build-flags.md` prices it under 1 per cent.
#
# THE BUILD IS CHECKED against SOCRATES's own key output before anything is
# timed. `examples/runes` ships a reference for gfortran, and this flag line
# reproduces it to a worst relative 1.14e-08 over 24 of about 400 printed
# values, which is a last digit at eight decimal places and is what
# -O2 -march=znver4 -funroll-loops buys over the reference's plain -O.
set -e
WORK="${1:?usage: socrates_cost_bench.sh <workdir>}"
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
SOC="$HERE/references/socrates"

mkdir -p "$WORK"
rm -rf "$WORK/socrates"
mkdir -p "$WORK/socrates"
cp -r "$SOC/src" "$SOC/make" "$SOC/sbin" "$SOC/build_code" "$WORK/socrates/"
ln -sfn "$SOC/data" "$WORK/socrates/data"
chmod +x "$WORK/socrates/build_code"

python3 - "$HERE" "$WORK" <<'PY'
import sys, pathlib, yaml
here, work = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
cfg = yaml.safe_load((here / "config/planet.yaml").read_text())
flags = [f for f in cfg["model"]["compile_flags"]["f90_opts"]
         if not f.startswith("-ffpe-trap")]
line = " ".join(flags)
# The LINK line drops the two flags that act on a COMPILE and not on a link,
# so the two lines say what each step actually does. gfortran accepts both at
# link and ignores them, which is exactly why they would go unnoticed.
link = " ".join(f for f in flags
                if f != "-cpp" and not f.startswith("-ffpe-summary"))
(work / "socrates/make/Mk_cmd").write_text(f"""#
# Vesper CLIM-61 cost bench. The flag line is config/planet.yaml's production
# f90_opts less -ffpe-trap, so this measures the code and not the compiler.
#
FORTCOMP        = gfortran {line} -c
LINK            = gfortran {link}
LIBLINK         = ar rvu
INCCDF_PATH     = /usr/include
LIBCDF_PATH     = /usr/lib
LIBCDF_NAME     = netcdff
OMPARG          = -fopenmp
FLAGS_ODEPACK   = -w -std=legacy -fallow-argument-mismatch

LIBSUFFIX       = a

.SUFFIXES: $(SUFFIXES) .f90
""")
import shutil, socket
shutil.copy(work / "socrates/make/Mk_cmd",
            work / f"socrates/make/Mk_cmd_{socket.gethostname()}")
print(line)
PY

( cd "$WORK/socrates" && ./build_code )

mkdir -p "$WORK/bench"
FL="$(python3 -c "
import yaml,sys
c=yaml.safe_load(open('$HERE/config/planet.yaml'))
print(' '.join(f for f in c['model']['compile_flags']['f90_opts']
                if not f.startswith('-ffpe-trap')))")"
( cd "$WORK/socrates/bin" \
  && gfortran $FL -c -I. "$HERE/exoplasim/scripts/socrates_cost_bench.f90" \
       -o socrates_bench.o \
  && gfortran $FL -o "$WORK/bench/socrates_bench" socrates_bench.o radlib.a )

cp "$HERE/references/socrates/data/spectra/ga7/sp_sw_ga7"* "$WORK/bench/"
cp "$HERE/references/socrates/data/spectra/ga7/sp_lw_ga7"* "$WORK/bench/"

# The check that can fail, before anything is timed.
cp "$HERE/references/socrates/data/spectra/ga9/sp_sw_ga9"* "$WORK/bench/"
cp "$HERE/references/socrates/data/spectra/ga9/sp_lw_ga9"* "$WORK/bench/"
( cd "$WORK/bench" && "$WORK/socrates/bin/runes_driver" > kgo_here.txt )
python3 - "$WORK/bench/kgo_here.txt" \
         "$HERE/references/socrates/examples/runes/gfortran_8_1_0.txt" <<'PY'
import sys, pathlib
got = pathlib.Path(sys.argv[1]).read_text().splitlines()
ref = pathlib.Path(sys.argv[2]).read_text().splitlines()
if len(got) != len(ref):
    sys.exit(f"runes key output has {len(got)} lines against {len(ref)}")
worst, differing = 0.0, 0
for a, b in zip(got, ref):
    fa, fb = a.split(), b.split()
    if len(fa) != len(fb):
        sys.exit(f"shape differs: {a!r} against {b!r}")
    for x, y in zip(fa, fb):
        try:
            x, y = float(x), float(y)
        except ValueError:
            continue
        if x != y:
            differing += 1
            worst = max(worst, abs(x - y) / max(abs(y), 1e-30))
print(f"runes key output: {differing} values differ, worst relative {worst:.2e}")
# The bar is 1e-7 and the reason is measured rather than chosen: this flag line
# reproduces the reference to 1.14e-08 at its worst, over 24 of about 400
# printed values, and the reference is written to eight decimal places so a
# last-digit disagreement IS about 1e-8. A bar at 1e-8 would refuse the build
# that passes; one at 1e-6 would let a real numerical difference through.
if worst > 1e-7:
    sys.exit("the build does not reproduce SOCRATES's own key output")
PY
echo "socrates cost bench built at $WORK/bench/socrates_bench"
