#!/bin/bash
# Where does the OpenMP barrier time go, per thread and per die?
#
# The barrier is one of the two costs a stripped libc and libgomp hide, and as
# a single number it says nothing about whether it is barrier machinery, a
# serial region, or threads waiting on a slower peer. This attributes it per
# thread and maps each thread to the core and die it ran on, which is what
# separates "the team is too big" from "the work is unevenly spread".
#
# Needs a FRAME-POINTER build: DWARF fails on 94% of samples here and frame
# pointers on none, at a measured cost of -1.17%. Build one with compile.sh -g
# and pass it as PROBE.
#
#   PROBE=/path/to/most_plasim_t170_l10_p16_omp_fp.x \
#     exoplasim/scripts/attribute_barrier_wait.sh [bed] [steps]
#
# A sample is classified by its LEAF: kernel if the leaf address is in the
# kernel half of the address space, barrier if the leaf is in a stripped
# library and no model frame appears anywhere in the stack, model otherwise.
# --sample-cpu is what makes the die column a measurement rather than an
# inference from thread creation order.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SCRATCH="${SCRATCH:-${TMPDIR:-/tmp}/attribute_barrier_wait}"
mkdir -p "$SCRATCH"
BED="${1:-bed_t170cold}"
STEPS="${2:-300}"
: "${PROBE:?set PROBE to a frame-pointer build, compile.sh -g}"

d="$SCRATCH/run"; rm -rf "$d"; mkdir -p "$d"
cp -a "$REPO/exoplasim/bench/$BED/." "$d/"
rm -f "$d"/*.x "$d/plasim_status" "$d/Abort_Message"; rm -rf "$d"/perf_*
sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $STEPS /" "$d/plasim_namelist"
sed -i "/^ *NSHTNS *=/d" "$d/plasim_namelist"
sed -i "2i\\ NSHTNS      =     1" "$d/plasim_namelist"
cp "$PROBE" "$d/probe.x"

( cd "$d" || exit 1
  ulimit -s unlimited
  OMP_NUM_THREADS=16 OMP_PROC_BIND=close OMP_PLACES=cores OMP_STACKSIZE=512M \
    perf record --quiet -F 97 -g --call-graph fp --sample-cpu --mmap-pages 32 \
        -o md.data -- ./probe.x > run.log 2>&1 )

lost=$(perf report -i "$d/md.data" --stdio 2>/dev/null | awk '/Total Lost Samples/{print $5}')
echo "lost samples: ${lost:-?}"

cd "$d" || exit 1
perf script -i md.data --comms probe.x --max-stack 40 -F tid,cpu,ip,sym,dso 2>/dev/null > s.txt
python3 - s.txt <<'PY'
import sys, re, collections

# Anchor the header pattern. An unanchored "digits [digits]" also matches
# callchain lines and silently triples the barrier share.
HDR   = re.compile(r"^\s*(\d+)\s+\[(\d+)\]\s*$")
FRAME = re.compile(r"^\s+([0-9a-f]+) (.*?) \((.*)\)\s*$")

def parse(path):
    tid = cpu = None; frames = []; out = []
    for line in open(path):
        if not line.strip():
            if tid is not None: out.append((tid, cpu, frames))
            tid = cpu = None; frames = []; continue
        m = FRAME.match(line)
        if m:
            frames.append((m.group(1), m.group(2))); continue
        h = HDR.match(line)
        if h:
            if tid is not None: out.append((tid, cpu, frames))
            tid, cpu, frames = h.group(1), int(h.group(2)), []
    if tid is not None: out.append((tid, cpu, frames))
    return out

def classify(frames):
    if not frames: return "empty"
    ip, _ = frames[0]
    if len(ip) == 16 and ip.startswith("ffff"): return "kernel"
    for _, sym in frames:
        if sym != "[unknown]": return "model"
    return "barrier"

st = parse(sys.argv[1])
overall = collections.Counter(classify(f) for _, _, f in st)
n = len(st)
print(f"\n{n} stacks")
for k, v in overall.most_common():
    print(f"    {k:>8} {v:7d}  {100*v/n:5.1f}%")

# Both SMT siblings of a core are the same core, and that is what
# OMP_PLACES=cores binds to; without the fold a bound thread looks unpinned.
nsib = (max(c for _, c, _ in st) + 1) // 2
tot = collections.Counter(); wait = collections.Counter()
cores = collections.defaultdict(collections.Counter)
for tid, cpu, fr in st:
    k = classify(fr)
    if k == "kernel": continue
    tot[tid] += 1; cores[tid][cpu % nsib] += 1
    if k == "barrier": wait[tid] += 1

rows = []
for t, m in tot.items():
    core, cn = cores[t].most_common(1)[0]
    rows.append((core, 100*cn/m, m, 100*wait[t]/m))
rows.sort()
half = nsib // 2
print(f"\n{'core':>5} {'on it':>7} {'samples':>8} {'barrier':>8}  die")
for core, pin, m, w in rows:
    die = f"CCD0/{'96MB'}" if core < half else f"CCD1/{'32MB'}"
    print(f"{core:>5} {pin:6.1f}% {m:8d} {w:7.1f}%  {die}")

d0 = [w for c, _, _, w in rows if c < half]
d1 = [w for c, _, _, w in rows if c >= half]
if d0 and d1:
    m0, m1 = sum(d0)/len(d0), sum(d1)/len(d1)
    lo = min(w for *_, w in rows); sp = max(w for *_, w in rows) - lo
    print(f"\nCCD0 mean {m0:5.2f}%   CCD1 mean {m1:5.2f}%   gap {m0-m1:+5.2f} points")
    print(f"floor {lo:5.2f}% is the critical path; the {sp:.2f}-point spread above "
          f"it is imbalance,")
    print(f"of which the die explains {100*abs(m0-m1)/sp:.0f}%.")
PY
