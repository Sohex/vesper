#!/bin/bash
# Does anything uninitialised reach a stored value, at the optimisation level
# this project ships?
#
#   exoplasim/scripts/verify_uninitialised_reads.sh <bed> <res> <n>
#
# Worldbuilding frame: a correctness check on the Vesper climate model's build.
# Nothing here is about the simulated planet.
#
# WHY THIS EXISTS AND WHY IT IS NOT THE TRAPPING FORM. `-finit-real=snan` with
# `-ffpe-trap=invalid` is the obvious check: poison every local real on entry
# and let the first read of one fault. world-5rs measured its positive control
# and it does not arm where it matters -- a deliberate uninitialised read faults
# at `-Og` and exits 0 at `-O2` and `-O3`, because gfortran folds the signalling
# NaN to a quiet one before the arithmetic and both binaries carry the folded
# constant in `.rodata` beside the signalling one. Production runs `-O2`. So the
# trapping form lives in the `poisoned` profile, which also carries `-Og`, and
# it certifies a build nothing integrates.
#
# `-finit-real=zero` is not the alternative. It was removed because it cost
# 26.03% of T170 and the restart was IDENTICAL without it, so it HID whatever it
# was protecting against rather than reporting it
# (`exoplasim/notes/the-zeroing-is-an-init-flag.md`). A check that can fail beats
# a tax that cannot.
#
# THE PROPAGATING FORM IS WHAT THAT TRADE ACTUALLY RESTS ON, and it works at
# `-O2` because it needs no trap to arm. Build the shipped flag line with
# `-finit-real=snan` added and FE_INVALID MASKED, run the same bed, and compare
# the restart bit for bit against the shipped build's. An uninitialised read
# then yields a NaN instead of whatever the stack held, and a NaN that reaches a
# stored field is just as visible in the file as a trap would have been. The
# fold to a quiet NaN does not defeat it: a quiet NaN in a restart record is
# still not the number the other arm wrote.
#
# WHAT A PASS DOES AND DOES NOT COVER. It covers reads whose value REACHES the
# restart within the run length. A local read before it is written, whose value
# is then discarded or overwritten before anything stored depends on it, is
# invisible here and is also harmless to the answer. It says nothing about
# integer or logical locals, which `-finit-real` does not touch, and nothing
# about a branch this bed does not take.
#
# THE MASKING IS THE PART TO GET RIGHT. `-ffpe-trap=invalid,zero,overflow` is in
# the declared flag line, so the poisoned arm DROPS that flag and adds back
# `-ffpe-trap=zero,overflow`. Dropping nothing and adding a second `-ffpe-trap`
# would leave the last one to win at the mercy of argument order, and adding
# `-ffpe-trap=none` would unmask the divide and the overflow as well, which are
# not what is under test here. `build_model.py --drop-flag` errors if the flag it
# is asked to remove is not in the declared line, so a change to
# `model.compile_flags.f90_opts` breaks this loudly rather than silently leaving
# the trap armed.
#
# TWO SWITCH SETTINGS, AND BOTH ARE A VERDICT.
#
#   NSHTNS=1 is what production runs, so it is the arm the question is asked
#   about.
#   NSHTNS=0 selects legmod and isolates the MODEL from the library. SHTns reads
#   a lane of its own scratch that it has not written -- located in
#   `the-zeroing-is-an-init-flag.md` down to the faulting instruction and the
#   address -- and under this arm's poison that lane can hold a signalling NaN
#   left by a Fortran frame below it. So a difference at NSHTNS=1 alone is
#   attributable to the library, and a difference at both is the model's own.
#
# Both must be bit identical, and the failure message says which case it is in.
# Neither is excused: a NaN reaching a stored value is a defect wherever it was
# born, and `the-zeroing-is-an-init-flag.md` records the library's over-read as
# safe only while the stack garbage under it is not a NaN.
#
# THE CONTROL is a deliberate uninitialised read that reaches a stored field:
# an external subroutine whose only local is never written, called once per
# `gridpointd` from inside the parallel region, folding that local into the
# temperature tendency through a factor that is zero at run time. In the shipped
# arm the local holds ordinary stack leftovers and the product is zero, so that
# arm integrates its usual trajectory; in the poisoned arm it is a NaN, zero
# times a NaN is a NaN, and the restart fills. External and not inline so that
# `-O2` cannot fold the read away, and scaled by a run-time zero rather than a
# literal so that neither arm's arithmetic is disturbed by the size of whatever
# it read. The two arms MUST differ. A gate nobody has seen fail is not a gate.
set -euo pipefail

bed="$(cd "${1:?usage: verify_uninitialised_reads.sh <bed> <res> <n>}" && pwd)"
res="${2:?}"
n="${3:?}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"
PKG="$REPO/vendor/exoplasim/exoplasim"
SRC="$PKG/plasim/src"
BUILD="$REPO/.venv/bin/python $REPO/exoplasim/scripts/build_model.py"
# exoplasim/bench/_uninit, or wherever UNINIT_WORK points. A git worktree
# reaches exoplasim/bench through a symlink into the MAIN checkout, which every
# worktree shares, so a worktree that runs this without setting the variable
# deletes and rewrites work another one is using -- this script opens with
# `rm -rf "$WORK"`.
WORK="${UNINIT_WORK:-$REPO/exoplasim/bench/_uninit}"

# DECLARED BEFORE THE ARMS RUN.
STEPS=60          # the length the 2026-08-22 measurement was taken at
# The verdict is BIT IDENTITY of plasim_status. There is no tolerance and no
# advisory record list, because the two arms are one source through two flag
# lines that change no arithmetic: `-finit-real` writes locals the model is
# meant to write itself, and `-ffpe-trap` sets MXCSR. Equality is the right
# answer, so any difference at all is the finding rather than a threshold
# question, and `verify_shtns_model.sh`'s advisory policy has nothing to
# forgive here. When they do differ, compare_restarts --norm names the record.

require_settled_bed "$bed"
require_bed_grid "$bed" "$res"

dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/bin"
restore() { ( cd "$REPO" && git checkout -- vendor/exoplasim/exoplasim/plasim/src ); }
trap restore EXIT

# The two flag lines. `shipped` is the declared line untouched; `poisoned` is
# that line with FE_INVALID masked and the local-real poison added.
# `--drop-flag=-x` and not `--drop-flag -x`: every one of these values starts
# with a dash, and argparse reads a dashed value as the next option.
POISON_FLAGS=(--drop-flag=-ffpe-trap=invalid,zero,overflow
              --extra-flag=-ffpe-trap=zero,overflow
              --extra-flag=-finit-real=snan)

build_arm() {
    local arm="$1"
    local stamp="$WORK/.stamp"
    local flags=()
    case "$arm" in poisoned*) flags=("${POISON_FLAGS[@]}") ;; esac
    : > "$stamp"
    # --no-publish: an ARM must not land in the model run directory under the
    # registry's naming, or the next run picks up a binary nobody asked for and
    # check_consistency reports unknown provenance. world-v3d.
    built=$( $BUILD --res "$res" --ranks "$n" --no-publish --print-path \
        "${flags[@]+"${flags[@]}"}" 2>"$WORK/build_$arm.log" ) || {
        echo "build failed: $arm (see $WORK/build_$arm.log)" >&2; exit 1; }
    [ -f "$built" ] && [ "$built" -nt "$stamp" ] || {
        echo "build produced nothing newer than the stamp: $arm "\
             "(see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$built" "$WORK/bin/$arm.x"
    echo "  built $arm  $(sha256sum "$WORK/bin/$arm.x" | cut -c1-16)"
}

# THE CONTROL PATCH. Appends an external subroutine with one never-written local
# and calls it once per gridpointd, inside the parallel region, on the
# temperature tendency. `real(nqspec-1,8)` is zero for every configuration this
# project runs and the compiler cannot fold it, so the shipped arm's trajectory
# is untouched and the poisoned arm's fills with NaN.
patch_control() {
    local before after
    before=$(grep -c '^      call ctl_uninit$' "$SRC/plasim.f90" || true)
    [ "$before" = 0 ] || {
        echo "control patch: plasim.f90 already carries the control" >&2; exit 1; }
    grep -q '^      mmrt(:,:)=0\.$' "$SRC/plasim.f90" || {
        echo "control patch: the gridpointd tendency-zeroing block has moved;" >&2
        echo "  this control inserts its call after 'mmrt(:,:)=0.' and that" >&2
        echo "  line is not there. Re-anchor it before trusting this gate." >&2
        exit 1; }
    sed -i 's/^      mmrt(:,:)=0\.$/      mmrt(:,:)=0.\n      call ctl_uninit/' \
        "$SRC/plasim.f90"
    cat >> "$SRC/plasim.f90" <<'FORTRAN'

!     CONTROL: a deliberate uninitialised read, inserted by
!     exoplasim/scripts/verify_uninitialised_reads.sh. External so -O2 cannot
!     fold it away; scaled by a run-time zero so the shipped arm is undisturbed.
      subroutine ctl_uninit
      use pumamod
      real :: zctl(NHOR)
      gtdt(:,:) = gtdt(:,:) + zctl(1) * real(nqspec-1,8)
      return
      end subroutine ctl_uninit
FORTRAN
    # A marker the smoke test refuses to see committed and build_model.py
    # refuses to publish under. A control patch is a deliberate corruption of
    # the model source, and one was once committed by a `git add -A` that ran
    # while a check was still working.
    sed -i '1i\!     CONTROL PATCH IN PROGRESS -- must not be committed' "$SRC/plasim.f90"
    after=$(grep -c '^      call ctl_uninit$' "$SRC/plasim.f90" || true)
    [ "$after" = 1 ] || {
        echo "control patch: expected one call site, found $after" >&2; exit 1; }
}

run_arm() {
    local arm="$1" sw="$2" tag="$3"
    local d="$WORK/run_$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp -a "$bed"/. "$d"/
    ( cd "$d"
      rm -f ./*.x plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $STEPS /" plasim_namelist
      sed -i "/^ *NSHTNS *=/d" plasim_namelist
      sed -i "2i\\ NSHTNS      =     $sw" plasim_namelist
      cp -f "$WORK/bin/$arm.x" ./probe.x
      export OMP_NUM_THREADS="$n" OMP_PROC_BIND=close OMP_PLACES=cores
      export OMP_STACKSIZE=512M
      ulimit -s unlimited
      ./probe.x >run.log 2>&1 ) >/dev/null 2>&1 || true
    [ -f "$d/plasim_status" ]
}

statsha() { sha256sum "$WORK/run_$1/plasim_status" | cut -c1-24; }

norm() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol 0 --norm 2>/dev/null || true
}

echo "$res, $n thread(s), $STEPS steps, bed $bed"
echo "declared before the arms ran: the poisoned arm's plasim_status must be"
echo "  BIT IDENTICAL to the shipped arm's at NSHTNS=1 and at NSHTNS=0, and the"
echo "  control must break that. No tolerance: the two flag lines change no"
echo "  arithmetic, so equality is the right answer."
echo "poisoned flag line: shipped ${POISON_FLAGS[*]}"
echo

rc=0
build_arm shipped
build_arm poisoned

for sw in 1 0; do
    case "$sw" in
      1) what="SHTns, what production runs" ;;
      0) what="legmod, the model isolated from the library" ;;
    esac
    echo
    echo "==== NSHTNS=$sw ($what): must be bit identical ===="
    if run_arm shipped "$sw" "shipped$sw" && run_arm poisoned "$sw" "poisoned$sw"; then
        a="$(statsha "shipped$sw")"; b="$(statsha "poisoned$sw")"
        if [ "$a" = "$b" ]; then
            echo "  [  ok  ] $a"
        else
            echo "  [ FAIL ] shipped $a, poisoned $b"
            echo "           worst record: $(norm "shipped$sw" "poisoned$sw")"
            echo "           Something uninitialised reached a stored value."
            if [ "$sw" = 1 ]; then
                echo "           If the NSHTNS=0 arm below agrees, the read is in"
                echo "           SHTns's own scratch rather than in the model."
            fi
            rc=1
        fi
    else
        echo "  [ FAIL ] an arm produced no restart (see $WORK/run_*/run.log)"
        rc=1
    fi
done

echo
echo "==== the control, a deliberate uninitialised read: must NOT agree ===="
patch_control
build_arm shipped_ctl
build_arm poisoned_ctl
if run_arm shipped_ctl 0 "cshipped" && run_arm poisoned_ctl 0 "cpoisoned"; then
    a="$(statsha cshipped)"; b="$(statsha cpoisoned)"
    if [ "$a" = "$b" ]; then
        echo "  [ FAIL ] the control agreed at $a, so this check cannot see an"
        echo "           uninitialised read reaching a stored field -- the defect"
        echo "           it exists to catch."
        rc=1
    else
        echo "  [  ok  ] control rejected: shipped $a, poisoned $b"
    fi
else
    # A control that cannot produce a restart has not been shown to differ for
    # the right reason, so it is not evidence either way.
    echo "  [ FAIL ] the control produced no restart, so it tests nothing"
    rc=1
fi
restore

echo
if [ "$rc" = 0 ]; then
    echo "PASS: nothing uninitialised reaches a stored value on this bed."
else
    echo "FAIL: see above."
fi
exit "$rc"
