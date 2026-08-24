#!/usr/bin/env python
"""Do legmod's three inverse transforms apply the table and the per-mode factor they are given?

Worldbuilding frame: a correctness check on the Vesper project's climate model.
Nothing here is about the simulated planet.

WHAT IS UNDER TEST. `sp2fc`, `sp2fcdmu` and `dv2uv` are the synthesis half of
the spectral transform, spectral coefficients to Fourier coefficients on the
latitudes a thread holds. Each is a matrix-vector product with `pmat` or `qmat`
scaled by a factor that depends only on the MODE, and that factor is the fork's
own change to these routines: `legini` stores the table and the factor apart --
`fsp` for the two scalar transforms, `fmu` and `fmv` for `dv2uv` -- and the
routines fold the factor into the FIELD once instead of into the weight once per
latitude. A factor applied to the wrong mode, dropped, or paired with the wrong
matrix is a model that runs and is wrong.

WHAT THE RIGHT ANSWER IS. The three loops are lifted verbatim out of
`legmod.f90` and compiled against a stub module whose `pmat`, `qmat`, `fsp`,
`fmu`, `fmv` and `plavor` this script supplies. The transform is BY DEFINITION
the matrix-vector product with what it was given, so the answer at every
latitude and every zonal wavenumber is known before the driver runs, and it is
computed here in numpy from the same tables. The tables themselves are an
associated Legendre table from scipy at Gauss-Legendre nodes, which keeps the
magnitudes those of a real transform; whether `legini` BUILDS that table
correctly is `verify_banded_transform.sh`'s question and is not asked here.

`dv2uv` is the case this exists for. It mixes `pmat*fmu` with `qmat*fmv` across
both Fourier components of both input fields, with four sign patterns, so
pattern-matching `sp2fc` onto it gives an answer that looks right. Each of its
four outputs is checked against its own analytic value at every latitude, every
wavenumber and every level.

PLANETARY VORTICITY. `dv2uv` takes ABSOLUTE vorticity and removes the planetary
part from its RESULT rather than from `pz` in place, because `pz` aliases the
shared spectral state and a subtract-restore around the loop is a write race.
The property that makes those two the same thing is an identity independent of
how the removal is coded:

    dv2uv(pd, pz) with plavor = P  ==  dv2uv(pd, pz - P at mode w=2) with P = 0

so the reference below simply subtracts P from the first Fourier component of
mode w=2 of the vorticity and takes the ordinary matrix-vector product. That
pins the mode, the component, the two outputs it touches, the two signs and the
two factors, with no appeal to the model's own formula. It does NOT pin the term
under the threaded grid bands, which is world-siq.

THE NEGATIVE CONTROLS. Seven, each the source with one line perturbed, each of
which must be rejected:

    sp2fc/wrongtable    pmat read where the routine's own matrix is pmat -> qmat
    sp2fc/nofactor      the per-mode factor dropped from the hoist
    sp2fcdmu/wrongtable qmat -> pmat
    sp2fcdmu/nofactor   the per-mode factor dropped from the hoist
    dv2uv/swapfactor    fmu and fmv exchanged on the vorticity scaling
    dv2uv/signflip      the sign joining dv2uv's two terms in pv component 2
    dv2uv/plavorsign    the sign of the planetary vorticity removal in pu

Every one of them is built by a textual substitution that must actually apply;
if a perturbation is a no-op because the line has moved, this exits rather than
running a control-free comparison. A gate that cannot fail is worse than one
that does not run.
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

# The planetary vorticity coefficient the dv2uv driver runs at. Any non-zero
# value of order one does; it is fixed here so the bound below is fixed too.
PLAVOR = 1.7

STUB = """
      module legmod
      integer, parameter :: NTRU = {ntru}
      integer, parameter :: NTP1 = NTRU + 1
      integer, parameter :: NCSP = NTP1 * (NTP1 + 1) / 2
      integer, parameter :: NLAT = {nlat}
      integer, parameter :: NLON = NLAT + NLAT
      integer, parameter :: NLPP = NLAT
      real :: pmat(NCSP,NLPP)
      real :: qmat(NCSP,NLPP)
      real :: fsp(NCSP)
      end module legmod

      program drive
      use legmod
      implicit none
      real :: sp(2,NCSP)
      real :: fc(2,NLON/2,NLPP)
      integer :: w, m, l
!     pmat, qmat and fsp all come from the reference, so the loops below are
!     being asked to reproduce the product of the two things they were handed.
!     fsp is NOT unity here: a per-mode factor of one hides a factor applied to
!     the wrong mode, and that is one of the two mistakes the hoist invites.
      open(20,file='pmat.dat',status='old')
      read(20,*) pmat
      close(20)
      open(21,file='qmat.dat',status='old')
      read(21,*) qmat
      close(21)
      open(22,file='fsp.dat',status='old')
      read(22,*) fsp
      close(22)
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


STUB_DV = """
      module legmod
      integer, parameter :: NTRU = {ntru}
      integer, parameter :: NTP1 = NTRU + 1
      integer, parameter :: NCSP = NTP1 * (NTP1 + 1) / 2
      integer, parameter :: NESP = NCSP + NCSP
      integer, parameter :: NLAT = {nlat}
      integer, parameter :: NLON = NLAT + NLAT
      integer, parameter :: NLEV = {nlev}
      integer, parameter :: NLPP = NLAT
      real :: plavor
      real :: pmat(NCSP,NLPP)
      real :: qmat(NCSP,NLPP)
      real :: fmu(NCSP)
      real :: fmv(NCSP)
      end module legmod

      program drive
      use legmod
      implicit none
      real :: pd(2,NESP/2,NLEV)
      real :: pz(2,NESP/2,NLEV)
      real :: pu(2,NLON/2,NLPP,NLEV)
      real :: pv(2,NLON/2,NLPP,NLEV)
      real :: za(4,NLEV)
      integer :: w, m, l, v
!     dv2uv's two weights are pmat*fmu and qmat*fmv, so the table and the
!     per-mode factors are supplied apart, exactly as legini stores them.
!     NLEV is more than one and the levels carry different fields, so a level
!     index pinned to 1 in the scaling loop is a wrong answer here.
      open(20,file='pmat.dat',status='old')
      read(20,*) pmat
      close(20)
      open(21,file='qmat.dat',status='old')
      read(21,*) qmat
      close(21)
      open(22,file='fmu.dat',status='old')
      read(22,*) fmu
      close(22)
      open(23,file='fmv.dat',status='old')
      read(23,*) fmv
      close(23)
      open(24,file='amp.dat',status='old')
      read(24,*) za
      close(24)
      plavor = {plavor}
      open(30,file='uv.dat',status='replace')
      do w = 1 , NCSP
         pz(:,:,:) = 0.0
         pd(:,:,:) = 0.0
         do v = 1 , NLEV
            pz(1,w,v) = za(1,v)
            pz(2,w,v) = za(2,v)
            pd(1,w,v) = za(3,v)
            pd(2,w,v) = za(4,v)
         enddo
         call dv2uv(pd,pz,pu,pv)
         do v = 1 , NLEV
            do l = 1 , NLPP
               do m = 1 , NTP1
                  write(30,'(4e26.17)') pu(1,m,l,v), pu(2,m,l,v),        &
     &                                  pv(1,m,l,v), pv(2,m,l,v)
               enddo
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


def perturb(text: str, name: str, old: str, new: str) -> str:
    """One substitution that must apply exactly once, or the gate refuses.

    A control built by a substitution that silently did nothing leaves the
    check comparing the source against itself, which passes and means nothing.
    """
    if text.count(old) != 1:
        raise SystemExit(f"could not build the {name} control: the line it "
                         f"perturbs occurs {text.count(old)} times in {SRC}, "
                         f"not once")
    return text.replace(old, new)


def legendre_tables(ntru: int, nlat: int):
    """pmat, qmat and the per-mode factors, north to south.

    Ordered north to south because that is how inigau hands the model its
    latitudes.  The tables are scipy's; the factors follow legini's own
    definitions -- fsp is the physics filter, fmu is m/(n(n+1)) times it and
    fmv is 1/(n(n+1)) times it -- with a filter that is not unity so that a
    factor carried to the wrong mode cannot hide.
    """
    mu, _ = np.polynomial.legendre.leggauss(nlat)
    mu = mu[::-1]                       # north first, as inigau orders them
    ncsp = (ntru + 1) * (ntru + 2) // 2
    pmat = np.zeros((ncsp, nlat))
    qmat = np.zeros((ncsp, nlat))
    modes = []
    lm = 0
    for m in range(ntru + 1):
        for n in range(m, ntru + 1):
            vals = assoc_legendre_p(n, m, mu, diff_n=1, norm=True)
            pmat[lm, :] = vals[0]
            qmat[lm, :] = vals[1]
            modes.append((m, n))
            lm += 1
    fsp = np.zeros(ncsp)
    fmu = np.zeros(ncsp)
    fmv = np.zeros(ncsp)
    for lm, (m, n) in enumerate(modes):
        # A filter shaped like the model's: indexed by TOTAL wavenumber, order
        # one, and different at every n. The shape does not matter; that it is
        # not constant is the whole point.
        filt = 1.0 / (1.0 + 0.37 * n)
        znn1 = 0.0 if n == 0 else 1.0 / (n * (n + 1))
        fsp[lm] = filt
        fmu[lm] = znn1 * m * filt
        fmv[lm] = znn1 * filt
    return pmat, qmat, fsp, fmu, fmv, modes


def build_and_run(stub: str, func_text: str, workdir: Path, tag: str,
                  files: dict, outname: str, **fmt):
    src = workdir / f"{tag}.f90"
    src.write_text(stub.format(**fmt) + "\n" + func_text + "\n")
    exe = workdir / f"{tag}.x"
    # -J and cwd both, because gfortran writes the .mod beside the WORKING
    # directory otherwise and leaves it in the repository.
    subprocess.run(["gfortran", "-fdefault-real-8", "-ffree-line-length-none",
                    "-J", str(workdir), "-o", str(exe), str(src)],
                   check=True, capture_output=True, text=True, cwd=workdir)
    for name, arr in files.items():
        np.savetxt(workdir / name, arr)
    subprocess.run([str(exe)], check=True, cwd=workdir, capture_output=True)
    return np.loadtxt(workdir / outname)


def run_sp(func_text: str, ntru, nlat, pmat, qmat, fsp, workdir, tag):
    files = {"pmat.dat": pmat.flatten(order="F"),
             "qmat.dat": qmat.flatten(order="F"),
             "fsp.dat": fsp}
    out = build_and_run(STUB, func_text, workdir, tag, files, "fc.dat",
                        ntru=ntru, nlat=nlat)
    ncsp = pmat.shape[0]
    # per mode: the sp2fc block then the sp2fcdmu block, each (NLPP, NTP1, 2)
    return out.reshape(ncsp, 2, nlat, ntru + 1, 2)


def run_dv(func_text: str, ntru, nlat, nlev, pmat, qmat, fmu, fmv, amp,
           workdir, tag):
    files = {"pmat.dat": pmat.flatten(order="F"),
             "qmat.dat": qmat.flatten(order="F"),
             "fmu.dat": fmu, "fmv.dat": fmv,
             "amp.dat": amp.flatten(order="F")}
    out = build_and_run(STUB_DV, func_text, workdir, tag, files, "uv.dat",
                        ntru=ntru, nlat=nlat, nlev=nlev,
                        plavor=repr(float(PLAVOR)))
    ncsp = pmat.shape[0]
    return out.reshape(ncsp, nlev, nlat, ntru + 1, 4)


def reference_sp(pmat, qmat, fsp, modes, ntru, nlat):
    """The transform IS the matrix-vector product with the table it was given.

    Only one mode is driven at a time, so every wavenumber row but that mode's
    own is zero -- which is half of what is being checked, because a weight
    written to the wrong row is exactly what a mis-stepped `w` produces.
    """
    ncsp = pmat.shape[0]
    want = np.zeros((ncsp, 2, nlat, ntru + 1, 2))
    for w, (m0, _n0) in enumerate(modes):
        for which, table in ((0, pmat), (1, qmat)):
            want[w, which, :, m0, 0] = table[w, :] * fsp[w] * 1.0
            want[w, which, :, m0, 1] = table[w, :] * fsp[w] * 2.0
    return want


def reference_dv(pmat, qmat, fmu, fmv, amp, modes, ntru, nlat, nlev):
    """dv2uv's four outputs, from the same tables, with planetary vorticity
    taken off the vorticity coefficient rather than off the result."""
    ncsp = pmat.shape[0]
    qu = pmat * fmu[:, None]
    qv = qmat * fmv[:, None]
    want = np.zeros((ncsp, nlev, nlat, ntru + 1, 4))
    for w, (m0, _n0) in enumerate(modes):
        for v in range(nlev):
            vor = np.zeros((2, ncsp))
            div = np.zeros((2, ncsp))
            vor[0, w], vor[1, w] = amp[0, v], amp[1, v]
            div[0, w], div[1, w] = amp[2, v], amp[3, v]
            # w = 2 in Fortran is mode index 1 here: zonal wavenumber 0,
            # total wavenumber 1, the only mode planetary vorticity occupies.
            vor[0, 1] -= PLAVOR
            for j in range(ncsp):
                if vor[0, j] == 0.0 and vor[1, j] == 0.0 \
                        and div[0, j] == 0.0 and div[1, j] == 0.0:
                    continue
                mj = modes[j][0]
                u, vq = qu[j, :], qv[j, :]
                want[w, v, :, mj, 0] += vq * vor[0, j] + u * div[1, j]
                want[w, v, :, mj, 1] += vq * vor[1, j] - u * div[0, j]
                want[w, v, :, mj, 2] += u * vor[1, j] - vq * div[0, j]
                want[w, v, :, mj, 3] += -u * vor[0, j] - vq * div[1, j]
    return want


def compare(got, want, scale, tol):
    """Absolute agreement, scaled by the magnitude the answer itself carries."""
    bad = np.abs(got - want) > tol * scale
    return int(bad.sum()), np.abs(got - want).max()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ntru", type=int, default=10)
    ap.add_argument("--nlat", type=int, default=16)
    ap.add_argument("--nlev", type=int, default=2)
    ap.add_argument("--tol", type=float, default=1e-12,
                    help="absolute agreement, scaled by the answer's own "
                         "magnitude. The transform is a sum of at most NTP1 "
                         "products of order one numbers, so anything above a "
                         "few ulp is a wrong answer rather than a rounding "
                         "difference.")
    args = ap.parse_args()

    text = (PROJECT_ROOT / SRC).read_text()
    sp2fc = extract(text, "sp2fc")
    sp2fcdmu = extract(text, "sp2fcdmu")
    dv2uv = extract(text, "dv2uv")

    hoist = ("   zsp(1,j) = fsp(j) * sp(1,j)\n"
             "   zsp(2,j) = fsp(j) * sp(2,j)")
    nohoist = ("   zsp(1,j) = sp(1,j)\n"
               "   zsp(2,j) = sp(2,j)")

    def sp_pair(a, b):
        return a + "\n" + b

    good_sp = sp_pair(sp2fc, sp2fcdmu)
    sp_controls = {
        "sp2fc/wrongtable": sp_pair(perturb(
            sp2fc, "sp2fc/wrongtable",
            "         fc(1,m,l) = fc(1,m,l) + pmat(w,l) * zsp(1,w)\n"
            "         fc(2,m,l) = fc(2,m,l) + pmat(w,l) * zsp(2,w)",
            "         fc(1,m,l) = fc(1,m,l) + qmat(w,l) * zsp(1,w)\n"
            "         fc(2,m,l) = fc(2,m,l) + qmat(w,l) * zsp(2,w)"), sp2fcdmu),
        "sp2fc/nofactor": sp_pair(
            perturb(sp2fc, "sp2fc/nofactor", hoist, nohoist), sp2fcdmu),
        "sp2fcdmu/wrongtable": sp_pair(sp2fc, perturb(
            sp2fcdmu, "sp2fcdmu/wrongtable",
            "         fc(1,m,l) = fc(1,m,l) + qmat(w,l) * zsp(1,w)\n"
            "         fc(2,m,l) = fc(2,m,l) + qmat(w,l) * zsp(2,w)",
            "         fc(1,m,l) = fc(1,m,l) + pmat(w,l) * zsp(1,w)\n"
            "         fc(2,m,l) = fc(2,m,l) + pmat(w,l) * zsp(2,w)")),
        "sp2fcdmu/nofactor": sp_pair(
            sp2fc, perturb(sp2fcdmu, "sp2fcdmu/nofactor", hoist, nohoist)),
    }

    dv_controls = {
        "dv2uv/swapfactor": perturb(
            dv2uv, "dv2uv/swapfactor",
            "     zzu(1,jm) = fmu(jm)*pz(1,jm,v) ; zzu(2,jm) = fmu(jm)*pz(2,jm,v)",
            "     zzu(1,jm) = fmv(jm)*pz(1,jm,v) ; zzu(2,jm) = fmv(jm)*pz(2,jm,v)"),
        "dv2uv/signflip": perturb(
            dv2uv, "dv2uv/signflip",
            "        pv(2,m,l,v)=pv(2,m,l,v)-pmat(w,l)*zzu(1,w)-qmat(w,l)*zdv(2,w)",
            "        pv(2,m,l,v)=pv(2,m,l,v)+pmat(w,l)*zzu(1,w)-qmat(w,l)*zdv(2,w)"),
        "dv2uv/plavorsign": perturb(
            dv2uv, "dv2uv/plavorsign",
            "    pu(1,1,l,v) = pu(1,1,l,v) - qmat(2,l) * fmv(2) * plavor",
            "    pu(1,1,l,v) = pu(1,1,l,v) + qmat(2,l) * fmv(2) * plavor"),
    }

    pmat, qmat, fsp, fmu, fmv, modes = legendre_tables(args.ntru, args.nlat)
    # The four amplitudes per level: vorticity components 1 and 2 then
    # divergence components 1 and 2. Different at every level, and no two of
    # the eight equal, so a component or a level read from the wrong place is
    # a wrong number rather than the same number.
    amp = np.array([[1.0, 7.0], [2.0, 11.0], [3.0, 13.0], [5.0, 17.0]])
    amp = amp[:, :args.nlev]

    want_sp = reference_sp(pmat, qmat, fsp, modes, args.ntru, args.nlat)
    want_dv = reference_dv(pmat, qmat, fmu, fmv, amp, modes,
                           args.ntru, args.nlat, args.nlev)
    scale_sp = np.maximum(np.abs(want_sp), 1.0)
    scale_dv = np.maximum(np.abs(want_dv), 1.0)

    print(f"NTRU {args.ntru}, NLAT {args.nlat}, NLEV {args.nlev}, "
          f"{len(modes)} modes, plavor {PLAVOR}, tolerance {args.tol:.0e}")

    failures = 0
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)

        n, worst = compare(run_sp(good_sp, args.ntru, args.nlat, pmat, qmat,
                                  fsp, wd, "sp"), want_sp, scale_sp, args.tol)
        if n:
            failures += 1
            print(f"[ FAIL ] sp2fc and sp2fcdmu: {n} disagreements, "
                  f"worst {worst:.3e}")
        else:
            print(f"[  ok  ] sp2fc and sp2fcdmu reproduce the table times the "
                  f"per-mode factor at every one of the {args.nlat} latitudes, "
                  f"both components (worst {worst:.3e})")

        n, worst = compare(run_dv(dv2uv, args.ntru, args.nlat, args.nlev, pmat,
                                  qmat, fmu, fmv, amp, wd, "dv"),
                           want_dv, scale_dv, args.tol)
        if n:
            failures += 1
            print(f"[ FAIL ] dv2uv: {n} disagreements, worst {worst:.3e}")
        else:
            print(f"[  ok  ] dv2uv reproduces all four outputs analytically at "
                  f"every latitude and level, planetary vorticity included "
                  f"(worst {worst:.3e})")

        print()
        print("the negative controls, which must fail")
        for i, (nm, txt) in enumerate(sorted(sp_controls.items())):
            n, worst = compare(run_sp(txt, args.ntru, args.nlat, pmat, qmat,
                                      fsp, wd, f"spc{i}"),
                               want_sp, scale_sp, args.tol)
            if n:
                print(f"[  ok  ] {nm} is rejected ({n} disagreements, "
                      f"worst {worst:.3e})")
            else:
                failures += 1
                print(f"[ FAIL ] the {nm} control PASSES, so this check proves "
                      f"nothing about it")

        for i, (nm, txt) in enumerate(sorted(dv_controls.items())):
            n, worst = compare(run_dv(txt, args.ntru, args.nlat, args.nlev,
                                      pmat, qmat, fmu, fmv, amp, wd, f"dvc{i}"),
                               want_dv, scale_dv, args.tol)
            if n:
                print(f"[  ok  ] {nm} is rejected ({n} disagreements, "
                      f"worst {worst:.3e})")
            else:
                failures += 1
                print(f"[ FAIL ] the {nm} control PASSES, so this check proves "
                      f"nothing about it")

    print()
    print(f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
