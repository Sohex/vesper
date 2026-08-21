#!/usr/bin/env python
"""Do the symmetric inverse transforms produce the right field at the mirror latitude?

`sp2fc` and `sp2fcdmu` now compute a mirror pair of latitudes from one pass over
the spectral modes, by summing the symmetric and antisymmetric modes into
separate accumulators and combining them twice. The mirror's weight matrix is
never read. That is only valid if

    qi(w,-mu) =  s(w) * qi(w,mu)        s(w) = (-1)**(m+n),  qi carries P
    qj(w,-mu) = -s(w) * qj(w,mu)                             qj carries dP/dmu

and getting the SIGN of the second one backwards produces a model that runs and
is wrong. So this is the test with a right answer: it is not a comparison
against another run.

HOW. The two transform loops are lifted verbatim out of `legmod.f90` and
compiled against a stub module whose `qi` and `qj` this script supplies -- an
associated Legendre table computed by scipy at Gauss-Legendre nodes, which is an
implementation the model shares nothing with. A single spectral mode is set and
the loops are run. The transform is by definition the matrix-vector product with
that table, so the answer at EVERY latitude is known in advance, including the
mirrors the paired branch never reads.

Supplying the table also makes the check free of any normalisation or phase
convention: whatever scipy's differs from legini's, the model is being asked to
reproduce the table it was given.

Both parities of mode are tested, and both Fourier components, at every mode of
the truncation.

And a NEGATIVE CONTROL: the mirror's combination is rebuilt with the two parity
sums swapped, which is the exact mistake the sign argument above exists to
prevent. It must fail.
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.special import assoc_legendre_p

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT  # noqa: E402

SRC = "vendor/exoplasim/exoplasim/plasim/src/legmod.f90"

STUB = """
      module legmod
      integer, parameter :: NTRU = {ntru}
      integer, parameter :: NTP1 = NTRU + 1
      integer, parameter :: NCSP = NTP1 * (NTP1 + 1) / 2
      integer, parameter :: NLAT = {nlat}
      integer, parameter :: NLON = NLAT + NLAT
      integer, parameter :: NLPP = NLAT
      integer, parameter :: NLHP = NLPP / 2
      logical, parameter :: LPAIRLAT = .true.
      real :: qi(NCSP,NLPP)
      real :: qj(NCSP,NLPP)
      end module legmod

      program drive
      use legmod
      implicit none
      real :: sp(2,NCSP)
      real :: fc(2,NLON/2,NLPP)
      integer :: w, m, l
      open(20,file='qi.dat',status='old')
      read(20,*) qi
      close(20)
      open(21,file='qj.dat',status='old')
      read(21,*) qj
      close(21)
      open(30,file='fc.dat',status='replace')
      do w = 1 , NCSP
         sp(:,:) = 0.0
         sp(1,w) = 1.0
         sp(2,w) = 2.0
         call sp2fc(sp,fc)
         do l = 1 , NLPP
            do m = 1 , NTP1
               write(30,'(2e26.17)') fc(1,m,l), fc(2,m,l)
            enddo
         enddo
         call sp2fcdmu(sp,fc)
         do l = 1 , NLPP
            do m = 1 , NTP1
               write(30,'(2e26.17)') fc(1,m,l), fc(2,m,l)
            enddo
         enddo
      enddo
      close(30)
      end program drive
"""


def extract(text: str, name: str) -> str:
    m = re.search(r"^subroutine %s\(.*?^end\n" % name, text, re.M | re.S)
    if m is None:
        raise SystemExit(f"{name} not found in {SRC}")
    return m.group(0)


def legendre_tables(ntru: int, nlat: int):
    """qi and qj at Gauss-Legendre nodes, north to south, from scipy.

    Ordered north to south because that is how inigau hands the model its
    latitudes, and it is what makes local l and NLPP+1-l a mirror pair.
    """
    mu, _ = np.polynomial.legendre.leggauss(nlat)
    mu = mu[::-1]                       # north first, as inigau orders them
    ncsp = (ntru + 1) * (ntru + 2) // 2
    qi = np.zeros((ncsp, nlat))
    qj = np.zeros((ncsp, nlat))
    modes = []
    lm = 0
    for m in range(ntru + 1):
        for n in range(m, ntru + 1):
            vals = assoc_legendre_p(n, m, mu, diff_n=1, norm=True)
            qi[lm, :] = vals[0]
            qj[lm, :] = vals[1]
            modes.append((m, n))
            lm += 1
    return qi, qj, modes


def run(func_text: str, ntru: int, nlat: int, qi, qj, workdir: Path):
    src = workdir / "drive.f90"
    src.write_text(STUB.format(ntru=ntru, nlat=nlat) + "\n" + func_text + "\n")
    exe = workdir / "drive.x"
    # -J and cwd both, because gfortran writes the .mod beside the WORKING
    # directory otherwise and leaves it in the repository.
    subprocess.run(["gfortran", "-fdefault-real-8", "-J", str(workdir),
                    "-o", str(exe), str(src)],
                   check=True, capture_output=True, text=True, cwd=workdir)
    np.savetxt(workdir / "qi.dat", qi.flatten(order="F"))
    np.savetxt(workdir / "qj.dat", qj.flatten(order="F"))
    subprocess.run([str(exe)], check=True, cwd=workdir, capture_output=True)
    out = np.loadtxt(workdir / "fc.dat")
    ncsp, ntp1 = qi.shape[0], ntru + 1
    # per mode: sp2fc block then sp2fcdmu block, each (NLPP, NTP1, 2)
    return out.reshape(ncsp, 2, nlat, ntp1, 2)


def check(got, qi, qj, modes, ntru: int, nlat: int, tol: float):
    """The transform IS the matrix-vector product with the table it was given."""
    bad = []
    ntp1 = ntru + 1
    for w, (m0, n0) in enumerate(modes):
        for which, table in ((0, qi), (1, qj)):
            for l in range(nlat):
                for mi in range(ntp1):
                    want1 = table[w, l] * 1.0 if mi == m0 else 0.0
                    want2 = table[w, l] * 2.0 if mi == m0 else 0.0
                    g1, g2 = got[w, which, l, mi, 0], got[w, which, l, mi, 1]
                    scale = max(abs(table[w, l]), 1.0)
                    if abs(g1 - want1) > tol * scale or abs(g2 - want2) > tol * scale:
                        bad.append((("sp2fc", "sp2fcdmu")[which], m0, n0, l, mi,
                                    g1, want1))
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ntru", type=int, default=10)
    ap.add_argument("--nlat", type=int, default=16)
    ap.add_argument("--tol", type=float, default=1e-12,
                    help="absolute agreement, scaled by the weight's own magnitude. "
                         "The transform is a sum of at most NTP1 products of order "
                         "one numbers, so anything above a few ulp is a wrong answer "
                         "rather than a rounding difference.")
    args = ap.parse_args()

    text = (PROJECT_ROOT / SRC).read_text()
    good = extract(text, "sp2fc") + "\n" + extract(text, "sp2fcdmu")

    # The controls: swap the two parity sums in the mirror's combination, once
    # in each routine. That is what a sign error in the parity relation looks
    # like, and the two routines carry OPPOSITE signs, so a control in only one
    # of them leaves the other untested.
    controls = {
        "sp2fc": good.replace("      fc(1,m,k) = ze1 - zo1\n      fc(2,m,k) = ze2 - zo2",
                              "      fc(1,m,k) = zo1 - ze1\n      fc(2,m,k) = zo2 - ze2"),
        "sp2fcdmu": good.replace("      fc(1,m,k) = zo1 - ze1\n      fc(2,m,k) = zo2 - ze2",
                                 "      fc(1,m,k) = ze1 - zo1\n      fc(2,m,k) = ze2 - zo2"),
    }
    for nm, txt in controls.items():
        if txt == good:
            raise SystemExit(f"could not build the {nm} control: the lines it "
                             f"perturbs have moved")

    qi, qj, modes = legendre_tables(args.ntru, args.nlat)
    print(f"NTRU {args.ntru}, NLAT {args.nlat}, {len(modes)} modes, "
          f"tolerance {args.tol:.0e}")

    failures = 0
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        bad = check(run(good, args.ntru, args.nlat, qi, qj, wd),
                    qi, qj, modes, args.ntru, args.nlat, args.tol)
        if bad:
            failures += 1
            print(f"[ FAIL ] {len(bad)} disagreements, first five:")
            for r in bad[:5]:
                print(f"         {r[0]} mode (m={r[1]},n={r[2]}) latitude {r[3]+1} "
                      f"wavenumber {r[4]+1}: got {r[5]:.17e} want {r[6]:.17e}")
        else:
            print(f"[  ok  ] every mode reproduces its own Legendre weight at "
                  f"every one of the {args.nlat} latitudes, both components")

        print()
        print("the negative controls, which must fail")
        for nm, txt in controls.items():
            bad = check(run(txt, args.ntru, args.nlat, qi, qj, wd),
                        qi, qj, modes, args.ntru, args.nlat, args.tol)
            if bad:
                print(f"[  ok  ] swapping the parity sums at the mirror in {nm} "
                      f"is rejected ({len(bad)} disagreements)")
            else:
                failures += 1
                print(f"[ FAIL ] the {nm} control PASSES, so this check proves "
                      f"nothing about it")

    print()
    print(f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
