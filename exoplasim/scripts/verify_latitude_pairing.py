#!/usr/bin/env python
"""Does `ilatperm` actually pair mirror latitudes on a process?

The paired latitude decomposition rests entirely on one integer function in
`plasimmod.f90`. Every grid-space transfer calls it, and if it is wrong the model
still runs -- it just distributes latitudes somewhere else and then applies a
north-south symmetry that no longer holds, which shows up as a plausible wrong
climate rather than as an error.

So this checks the REAL function text, not a reimplementation of it: the source
is lifted verbatim out of `plasimmod.f90`, compiled inside a parameter module
that supplies one (NLAT, NPRO), run, and its map checked against properties
that a wrong indexing cannot satisfy:

  1. it is a bijection of the slots onto the latitudes;
  2. slot `l` and slot `NLPP+1-l` of the same process are MIRROR latitudes,
     `g` and `NLAT+1-g`, which is the property the symmetric transform needs;
  3. the first NLHP slots of process `r` are the contiguous northern block
     `r*NLHP+1 .. (r+1)*NLHP`, so latitude 1 is northernmost as `inigau`
     orders it and the symmetric branch's reference latitude is the northern
     one, as its sign convention assumes;
  4. at NPRO == 1 it is the identity, which is what makes the single-process
     executable the trivial case rather than a special one;
  5. where NPRO does not divide NLAT/2 it is the identity too, the documented
     fallback to the stock contiguous layout.

And a NEGATIVE CONTROL: the same checks are run against a deliberately broken
variant, which must fail. A test that cannot fail is not a test -- the filter
fold shipped with one of those, and this is the same class of change.
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT  # noqa: E402

SRC = "vendor/exoplasim/exoplasim/plasim/src/plasimmod.f90"

# (NLAT, NPRO) covering the ladder and both sides of the divisibility fallback.
CASES = [(32, 1), (32, 8), (32, 16), (32, 32),
         (64, 1), (64, 8), (64, 16), (64, 32),
         (128, 8), (128, 16), (128, 32),
         (192, 8), (192, 16), (192, 32),
         (256, 8), (256, 16), (256, 32)]

PARAMOD = """
      module pumamod
      integer, parameter :: NLAT = {nlat}
      integer, parameter :: NPRO = {npro}
      integer, parameter :: NLPP = NLAT / NPRO
      integer, parameter :: NLHP = NLPP / 2
      logical, parameter :: LPAIRLAT = (mod(NLAT,2*NPRO) == 0)
      contains
{func}
      end module pumamod

      program dump
      use pumamod
      integer :: j
      do j = 1 , NLAT
         write(*,'(i6)') ilatperm(j)
      enddo
      end program dump
"""


def extract(text: str) -> str:
    """The ilatperm function, verbatim, from mpimod.f90."""
    m = re.search(r"^      integer function ilatperm.*?^      end function ilatperm",
                  text, re.M | re.S)
    if m is None:
        raise SystemExit("ilatperm not found in " + SRC)
    return m.group(0)


def run_case(func: str, nlat: int, npro: int, workdir: Path) -> list[int]:
    f = workdir / f"perm_{nlat}_{npro}.f90"
    f.write_text(PARAMOD.format(nlat=nlat, npro=npro, func=func))
    exe = workdir / f"perm_{nlat}_{npro}.x"
    # -J and cwd both, because gfortran writes the .mod beside the WORKING
    # directory otherwise and leaves it in the repository.
    subprocess.run(["gfortran", "-ffixed-line-length-132", "-fcheck=all",
                    "-J", str(workdir), "-o", str(exe), str(f)], check=True,
                   capture_output=True, text=True, cwd=workdir)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
    return [int(x) for x in out.stdout.split()]


def check(perm: list[int], nlat: int, npro: int) -> list[str]:
    """Every way this map can be wrong, stated so it can fail."""
    bad = []
    nlpp = nlat // npro
    nlhp = nlpp // 2
    paired = (nlat % (2 * npro)) == 0

    if sorted(perm) != list(range(1, nlat + 1)):
        bad.append("not a bijection of 1..NLAT")
        return bad                      # every other check reads garbage

    if not paired:
        if perm != list(range(1, nlat + 1)):
            bad.append("NPRO does not divide NLAT/2, so the map must be the "
                       "identity fallback and is not")
        return bad

    if npro == 1 and perm != list(range(1, nlat + 1)):
        bad.append("one process, so the map must be the identity and is not")

    for r in range(npro):
        block = perm[r * nlpp:(r + 1) * nlpp]
        for l in range(1, nlhp + 1):
            g, gm = block[l - 1], block[nlpp - l]
            if gm != nlat + 1 - g:
                bad.append(f"process {r}: slots {l} and {nlpp + 1 - l} hold "
                           f"latitudes {g} and {gm}, which are not mirrors")
        want = list(range(r * nlhp + 1, (r + 1) * nlhp + 1))
        if block[:nlhp] != want:
            bad.append(f"process {r}: northern half is {block[:nlhp]}, "
                       f"expected the contiguous block {want}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()

    text = (PROJECT_ROOT / SRC).read_text()
    good = extract(text)

    # Two controls, because the two ways this function goes wrong fail
    # different checks and a test that only catches one of them is half a test.
    #
    #   stride  the northern block advances by NLPP instead of NLHP. Latitudes
    #           collide, so it is not even a bijection.
    #   shift   the southern blocks are handed round the processes by one. It
    #           IS a bijection and every process still holds NLPP latitudes;
    #           they are simply not each other's mirrors. This is the subtler
    #           mistake and the one a bijection check alone would wave through.
    controls = {
        "stride": good.replace("ilatperm = ir * NLHP + il",
                               "ilatperm = ir * NLPP + il"),
        "shift": good.replace("ilatperm = NLAT - ir * NLHP - NLPP + il",
                              "ilatperm = NLAT - mod(ir+1,NPRO) * NLHP - NLPP + il"),
    }
    for nm, txt in controls.items():
        if txt == good:
            raise SystemExit(f"could not build the {nm} control: the line it "
                             f"perturbs has moved")

    failures = 0
    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)
        print("the function as it stands")
        for nlat, npro in CASES:
            bad = check(run_case(good, nlat, npro, workdir), nlat, npro)
            tag = f"NLAT {nlat:4d}  NPRO {npro:3d}"
            if bad:
                failures += 1
                print(f"[ FAIL ] {tag}")
                for b in bad:
                    print(f"         {b}")
            else:
                print(f"[  ok  ] {tag}")

        print()
        print("the negative controls, which must fail")
        for nm, txt in controls.items():
            caught = 0
            checked = 0
            for nlat, npro in CASES:
                if npro == 1 or (nlat % (2 * npro)) != 0:
                    continue        # identity either way; a control cannot show there
                checked += 1
                if check(run_case(txt, nlat, npro, workdir), nlat, npro):
                    caught += 1
            if caught == checked:
                print(f"[  ok  ] all {checked} paired cases reject the {nm} control")
            else:
                failures += 1
                print(f"[ FAIL ] {checked - caught} of {checked} paired cases "
                      f"ACCEPT the {nm} control, so these checks prove nothing")

    print()
    print(f"{len(CASES)} cases, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
