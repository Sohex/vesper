#!/bin/bash
# Where the model's data comes from: L2, this die's L3, the OTHER die, or DRAM.
#
#   exoplasim/scripts/profile_memory.sh <bed_dir> <binary> <threads> [nshtns]
#
# Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
# desktop. Nothing here is about the simulated planet.
#
# WHY THIS EXISTS. The cpu-clock profile says where cycles go and cannot say
# why. `docs/src/reference/environment.md` sets a 32 MB per-die working set as a
# standing target on the assumption that spilling past a die's L3 costs, and
# nothing in this repo has ever measured a fill. This does: it counts where
# every demand fill was served from, so the target can be checked rather than
# believed.
#
# THE EVENTS ARE ZEN-SPECIFIC AND CHOSEN FOR THIS PROCESSOR. Generic LLC-loads
# and LLC-load-misses report "not supported" here -- AMD does not map them --
# and there is no uncore PMU exposed, so there is no direct DRAM bandwidth
# counter either. `ls_any_fills_from_sys` is the core-side answer and its
# breakdown is exactly the question this machine poses:
#
#   local_l2      served by this core's own L2
#   local_ccx     served by this CCX's L3 -- the 32 MB on CCD1, 96 MB on CCD0
#   far_cache     served by the OTHER die's cache, which is a cross-die hop
#   dram_io_all   went to memory
#
# far_cache is the one worth watching on a 7950X3D. A thread team spread across
# both dies pays for every line it takes from the other one, and libgomp does
# not bind by default -- which is why the launcher pins.
#
# TWO PASSES, NOT ONE, and the reason is the counter file. Zen exposes six
# programmable core counters; asking for more multiplexes them, and a
# multiplexed count is an extrapolation from a fraction of the run rather than a
# measurement. Each pass here fits, and `cycles` appears in both so the two can
# be put on one denominator.
set -euo pipefail

bed="$(cd "${1:?usage: profile_memory.sh <bed_dir> <binary> <threads> [nshtns]}" && pwd)"
binary="${2:?binary name inside the bed}"
threads="${3:?thread count}"
nshtns="${4:-}"

[ -x "$bed/$binary" ] || { echo "no executable $bed/$binary" >&2; exit 1; }
[ -f "$bed/bed_manifest.json" ] || { echo "$bed is not a bed" >&2; exit 1; }

# USER TIME ONLY, pinned with :u on every event rather than left to the
# ambient setting. perf_event_paranoid decides the default -- at 2 it can only
# count user, at 0 it counts kernel too -- so an unpinned event set silently
# measures a different thing before and after someone runs sysctl. That
# happened here, between two arms of one comparison.
PASS_A="cycles:u,instructions:u,L1-dcache-loads:u,L1-dcache-load-misses:u"
PASS_B="cycles:u,ls_any_fills_from_sys.local_l2:u,ls_any_fills_from_sys.local_ccx:u,ls_any_fills_from_sys.far_cache:u,ls_any_fills_from_sys.dram_io_all:u"

cd "$bed"
if [ -n "$nshtns" ]; then
    sed -i "/^ *NSHTNS *=/d" plasim_namelist
    sed -i "2i\\ NSHTNS      =     $nshtns" plasim_namelist
fi
rm -f plasim_status Abort_Message

run () {
    local events="$1" out="$2"
    OMP_NUM_THREADS="$threads" OMP_PROC_BIND=close OMP_PLACES=cores \
    OMP_STACKSIZE=512M \
        bash -c "ulimit -s unlimited; exec perf stat -x, -e $events -o $out -- ./$binary" \
        >/dev/null 2>&1
}

echo "== warm-up, uncounted =="
OMP_NUM_THREADS="$threads" OMP_PROC_BIND=close OMP_PLACES=cores OMP_STACKSIZE=512M \
    bash -c "ulimit -s unlimited; exec ./$binary" >/dev/null 2>&1
[ -e Abort_Message ] && { echo "model aborted during warm-up" >&2; exit 1; }

echo "== pass A: instructions and L1 =="
run "$PASS_A" perf_stat_a.csv
echo "== pass B: where the fills came from =="
run "$PASS_B" perf_stat_b.csv

python3 - "$bed" "${nshtns:-default}" <<'PY'
import sys, csv, pathlib

bed, sw = pathlib.Path(sys.argv[1]), sys.argv[2]
vals, scaled = {}, []
for name in ("perf_stat_a.csv", "perf_stat_b.csv"):
    for row in csv.reader((bed / name).open()):
        if not row or row[0].startswith("#") or len(row) < 3:
            continue
        try:
            count = float(row[0])
        except ValueError:
            continue          # "<not supported>" and friends
        ev = row[2].split(":")[0]
        vals[ev] = count
        # perf appends the fraction of the run an event was actually counted
        # for when it multiplexes. Anything below 100 means extrapolation.
        if len(row) > 5 and row[5]:
            try:
                if float(row[5]) < 99.0:
                    scaled.append((ev, float(row[5])))
            except ValueError:
                pass

g = vals.get
cyc, ins = g("cycles", 0), g("instructions", 0)
l1l, l1m = g("L1-dcache-loads", 0), g("L1-dcache-load-misses", 0)
l2 = g("ls_any_fills_from_sys.local_l2", 0)
ccx = g("ls_any_fills_from_sys.local_ccx", 0)
far = g("ls_any_fills_from_sys.far_cache", 0)
dram = g("ls_any_fills_from_sys.dram_io_all", 0)
fills = l2 + ccx + far + dram

print()
print(f"NSHTNS={sw}")
print(f"  instructions per cycle      {ins / cyc:8.3f}" if cyc else "  no cycles counted")
if l1l:
    print(f"  L1 dcache load miss rate    {100 * l1m / l1l:8.3f}%   "
          f"({l1m:,.0f} of {l1l:,.0f})")
if fills:
    print()
    print("  where demand fills were served from")
    for label, n in (("local L2", l2), ("this die's L3", ccx),
                     ("the OTHER die", far), ("DRAM", dram)):
        print(f"    {label:<16} {100 * n / fills:7.2f}%   {n:>16,.0f}")
    # A fill is a 64-byte line, so this is the traffic the core actually pulled.
    print(f"  bytes from DRAM             {dram * 64 / 1024**3:8.2f} GiB")
    print(f"  bytes from the other die    {far * 64 / 1024**3:8.2f} GiB")
if scaled:
    print()
    print("  WARNING: these were multiplexed and are extrapolations, not counts:")
    for ev, pct in scaled:
        print(f"    {ev} counted for {pct:.1f}% of the run")
PY
