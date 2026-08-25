#!/usr/bin/env python
"""Are the eight weights the transform sites need pmat/qmat times legini's own factors?

Worldbuilding frame: a correctness check on the Vesper project's climate model.
Nothing here is about the simulated planet.

WHAT CLIM-48 CLAIMED. `legini` used to build eight `NCSP x NLPP` matrices, and
the model touched every one each timestep, which is 48.4 MB a die at T127 and
114.9 MB at T170 against CCD1's 32 MB of L3 -- the measured cause of the load
imbalance in `notes/rank-imbalance-and-weight-traffic.md`. If all eight are one
of two arrays scaled by a per-MODE factor and a per-LATITUDE factor, then two
matrices plus a handful of vectors hold the same information and the resident
footprint falls fourfold. `legini` now stores exactly that: `pmat`, `qmat`,
`gwdc` and the six vectors `fsp`, `fgp`, `fmu`, `fmv`, `fmm`, `fmq`.

WHAT THIS CHECKED BEFORE, AND WHY IT WAS NOT A CHECK. It reproduced legini's
recurrence in Python, built the eight matrices from that transcription, rebuilt
each from the same transcription's P and Q, and compared. Both sides were the
identical product reassociated, so it proved float reassociation and nothing
about legini; a transcription error anywhere in the Python -- in the recurrence
that fed BOTH sides most of all -- cancelled exactly, and the negative control
was caught only because it perturbed one side of a self-comparison. world-g5va.

WHAT IT CHECKS NOW. Two sides that share no code.

THE SUBJECT is `legini` itself, lifted verbatim out of `legmod.f90` by
`extract` and compiled against a stub module that supplies its inputs. The
tables and the six factor vectors compared here are the ones the model's own
initialisation produces, at the declared precision, from the source in the tree.

THE REFERENCE is scipy, through `banded_transform_reference.legendre_table` --
`sph_legendre_p_all`, which shares no recurrence with `legini` -- and the
quadrature from `gauss_weight_reference.gauss_legendre` in extended precision.
The six per-mode factors are written here from the filter's declared formula.
`exoplasim/notes/banded-transform-reference.md` states what that reference is
independent of and what it is not; read it before reading a result off this.

WHAT THE RIGHT ANSWER IS. Each of the eight weights is a stated product of a
table and factors, so its value at every mode and every latitude is known
before the driver runs. A factor carried to the wrong mode, paired with the
wrong table, or built from the wrong filter is a different number here, not the
same number in a different order.

WHICH TABLE PAIRS WITH WHICH FACTOR IS READ OUT OF THE SOURCE, not asserted
here. `pairings` extracts each transform routine from `legmod.f90` and requires
the pairing this file assumes to appear in it, so a routine that starts
multiplying `pmat` by `fmv` makes this gate refuse rather than pass on a
description that has gone stale.

WHAT THIS DOES NOT COVER. `inigau`: the nodes and weights are supplied to
`legini` rather than computed by it, and `verify_gauss_weights.sh` is where that
is checked. The transform LOOPS: `verify_inverse_transform.py` drives those
against the matrix-vector product with the tables they are handed, and
`verify_banded_transform.sh` drives the analysis direction across thread bands.
And the model: whether the separated form integrates the same climate is
`verify_weight_factorisation_model.sh`.

THE NEGATIVE CONTROLS are textual perturbations of the extracted Fortran, each
of which must be caught. Every one must apply exactly once or the gate refuses,
so a control that has become a no-op cannot leave the comparison control-free.

ONE MISTAKE IS INVISIBLE HERE AND IS NAMED RATHER THAN LEFT OUT. Crossing the
two filters -- carrying `skgpsp` where `skspgp` belongs -- cannot be detected by
anything, because `legini` builds both from ONE shape selected by `nfilter` and
turns each on with `ngptfilter` and `nspvfilter`, so with both on the two arrays
are identical element for element. A control that swaps them is a no-op, and a
gate that reported it as caught would be reporting the shape it shares rather
than the direction it carries. What IS checked is that the filter is read at the
TOTAL wavenumber, which is the indexing mistake the same block invites.
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import banded_transform_reference as btr  # noqa: E402
from _paths import PROJECT_ROOT  # noqa: E402
from gauss_weight_reference import gauss_legendre  # noqa: E402
import rungs  # noqa: E402  -- the ladder registry; _paths put lib on the path

SRC = PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src" / "legmod.f90"

# The eight weights, each as (which table, which per-mode factor, which
# per-latitude factor). This is the table legmod.f90's own declarations carry
# beside `fsp` and its siblings; `pairings` below requires the source to agree.
WEIGHTS = {
    "qi": ("pmat", "fsp", None),
    "qj": ("qmat", "fsp", None),
    "qc": ("pmat", "fgp", "gwd"),
    "qe": ("qmat", "fgp", "gwdc"),
    "qu": ("pmat", "fmu", None),
    "qv": ("qmat", "fmv", None),
    "qm": ("pmat", "fmm", "gwdc"),
    "qq": ("pmat", "fmq", "gwdc"),
}
NAMES = tuple(WEIGHTS)
FACTORS = ("fsp", "fgp", "fmu", "fmv", "fmm", "fmq")

# What each transform routine must be seen doing, so that the table above is
# read out of the source rather than remembered. The hoisted forms scale the
# FIELD by the factor and then multiply by the bare table, which is the whole
# point of storing them apart, so the requirement is on both halves.
PAIRING_EVIDENCE = {
    "fc2sp": ["pmat(w,l)*fgp(w)*gwd(l)"],
    "sp2fc": ["zsp(1,j) = fsp(j) * sp(1,j)", "pmat(w,l) * zsp(1,w)"],
    "sp2fcdmu": ["zsp(1,j) = fsp(j) * sp(1,j)", "qmat(w,l) * zsp(1,w)"],
    "dv2uv": ["zzu(1,jm) = fmu(jm)*pz(1,jm,v)", "zzv(1,jm) = fmv(jm)*pz(1,jm,v)",
              "pmat(w,l)*zzu(2,w)", "qmat(w,l)*zzv(1,w)",
              "qmat(2,l) * fmv(2) * plavor", "pmat(2,l) * fmu(2) * plavor"],
    "mktend": ["qmat(w,l)*fgp(w)*gwdc(l)", "pmat(w,l)*fmm(w)*gwdc(l)",
               "pmat(w,l)*fmq(w)*gwdc(l)"],
}

# Each control is (name, what it breaks, the line, the line it becomes). The
# substitution must apply exactly once to the extracted legini.
CONTROLS = [
    ("qm_mode", "fmm loses the zonal wavenumber it carries",
     "      fmm(lm) = m * skgpsp(n+1)",
     "      fmm(lm) = n * skgpsp(n+1)"),
    ("qu_mode", "fmu loses the zonal wavenumber, so qu becomes qv",
     "      fmu(lm) = znn1 * m * skspgp(n+1)",
     "      fmu(lm) = znn1 * skspgp(n+1)"),
    ("filter_index", "the filter is read at the zonal wavenumber, not the total",
     "      fgp(lm) = skgpsp(n+1)",
     "      fgp(lm) = skgpsp(m+1)"),
    ("recurrence", "one sign in the P recurrence",
     "         z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)",
     "         z2 = zsin * zpli(lm-1) + z1 * zpli(lm-2)"),
    ("qq_gain", "fmq's n(n+1)/2 becomes n(n+1)",
     "      fmq(lm) = n * (n+1) * 0.5_8 * skgpsp(n+1)",
     "      fmq(lm) = n * (n+1) * skgpsp(n+1)"),
]

STUB = """
      module legmod
      integer, parameter :: NTRU = {ntru}
      integer, parameter :: NTP1 = NTRU + 1
      integer, parameter :: NCSP = NTP1 * (NTP1 + 1) / 2
      integer, parameter :: NLAT = {nlat}
      integer, parameter :: NLPP = NLAT
      integer, parameter :: NROOT = 0
!     Not NROOT, so legini's diagnostic print of the filter stays out of the
!     way. It writes to nud and says nothing this reads.
      integer :: mypid = 1
      integer :: nud = 6
      integer :: nfilter = {nfilter}
      integer :: ngptfilter = 1
      integer :: nspvfilter = 1
      integer :: nfilterexp = {nfilterexp}
      real :: filterkappa = {filterkappa}
      real :: landhoskn0 = {landhoskn0}
!     inigau's outputs, supplied rather than computed: whether the model builds
!     these correctly is verify_gauss_weights.sh's question.
      real :: gwd(NLPP)
      real :: sid(NLPP)
      real :: pmat(NCSP,NLPP)
      real :: qmat(NCSP,NLPP)
      real :: fsp(NCSP), fgp(NCSP), fmu(NCSP), fmv(NCSP), fmm(NCSP), fmq(NCSP)
      real :: gwdc(NLPP)
      real :: skgpsp(NTP1), skspgp(NTP1)
      end module legmod

      program drive
      use legmod
      implicit none
      open(20,file='sid.dat',status='old')
      read(20,*) sid
      close(20)
      open(21,file='gwd.dat',status='old')
      read(21,*) gwd
      close(21)
      call legini
      open(30,file='pmat.out',status='replace')
      write(30,'(es26.17e3)') pmat
      close(30)
      open(31,file='qmat.out',status='replace')
      write(31,'(es26.17e3)') qmat
      close(31)
      open(32,file='fac.out',status='replace')
      write(32,'(es26.17e3)') fsp, fgp, fmu, fmv, fmm, fmq, gwdc
      close(32)
      end program drive
"""


def extract(text: str, name: str) -> str:
    """One subroutine, verbatim, from `legmod.f90`.

    `legini` takes no arguments, so the pattern is the bare statement rather
    than the parenthesised one `verify_inverse_transform.extract` matches.
    """
    m = re.search(r"^subroutine %s\b.*?^end$" % name, text, re.M | re.S)
    if m is None:
        raise SystemExit(f"{name} not found in {SRC}")
    return m.group(0)


def perturb(text: str, name: str, old: str, new: str) -> str:
    """One substitution that must apply exactly once, or the gate refuses.

    A control built by a substitution that silently did nothing leaves the
    check comparing the source against itself, which passes and means nothing.
    """
    if text.count(old) != 1:
        raise SystemExit(
            f"could not build the {name} control: the line it perturbs occurs "
            f"{text.count(old)} times in legini, not once. The control has to "
            f"be re-derived against the source as it now stands.")
    return text.replace(old, new)


def pairings(text: str) -> list[str]:
    """Require each transform routine to pair the tables and factors WEIGHTS says.

    Read out of the source rather than remembered, so a routine that starts
    multiplying `pmat` by `fmv` makes this refuse instead of passing on a
    description that has gone stale.
    """
    missing = []
    for routine, evidence in PAIRING_EVIDENCE.items():
        body = extract(text, routine)
        for line in evidence:
            if line not in body:
                missing.append(f"{routine}: {line}")
    return missing


def run_legini(legini_text: str, ntru: int, nlat: int, filt: dict,
               nodes, weights, workdir: Path, tag: str) -> dict:
    """Compile the extracted legini and return everything it built."""
    src = workdir / f"{tag}.f90"
    src.write_text(STUB.format(ntru=ntru, nlat=nlat, **filt) + "\n"
                   + legini_text + "\n")
    exe = workdir / f"{tag}.x"
    # -J and cwd both, because gfortran writes the .mod beside the WORKING
    # directory otherwise and leaves it in the repository.
    subprocess.run(["gfortran", "-fdefault-real-8", "-ffree-line-length-none",
                    "-J", str(workdir), "-o", str(exe), str(src)],
                   check=True, capture_output=True, text=True, cwd=workdir)
    np.savetxt(workdir / "sid.dat", np.asarray(nodes, dtype=np.float64))
    np.savetxt(workdir / "gwd.dat", np.asarray(weights, dtype=np.float64))
    subprocess.run([str(exe)], check=True, cwd=workdir, capture_output=True)
    ncsp = (ntru + 1) * (ntru + 2) // 2
    out = {
        "pmat": np.loadtxt(workdir / "pmat.out").reshape(nlat, ncsp).T,
        "qmat": np.loadtxt(workdir / "qmat.out").reshape(nlat, ncsp).T,
    }
    fac = np.loadtxt(workdir / "fac.out")
    for i, name in enumerate(FACTORS):
        out[name] = fac[i * ncsp:(i + 1) * ncsp]
    out["gwdc"] = fac[len(FACTORS) * ncsp:]
    return out


def filter_curve(ntru: int, filt: dict) -> np.ndarray:
    """skgpsp and skspgp, from the declared formula rather than from legini.

    Indexed by n from 1 to NTP1, which is the model's own indexing: legmod.f90
    counts n from zero in the mode loop and reads `skgpsp(n+1)`. Only the
    filters legini builds are here; anything else makes this refuse rather than
    fall back to unity, because a filter that is silently one hides every
    factor that carries it.
    """
    ntp1 = ntru + 1
    n = np.arange(1, ntp1 + 1, dtype=np.float64)
    nf = filt["nfilter"]
    if nf == 0:
        return np.ones(ntp1)
    if nf == 1:                                   # Cesaro
        return 1.0 - n / ntp1
    if nf == 2:                                   # exponential
        return np.exp(-filt["filterkappa"] * (n / ntru) ** filt["nfilterexp"])
    if nf == 3:                                   # Lander-Hoskins
        n0 = filt["landhoskn0"] * ntru / 21.0
        return np.exp(-((n * (n + 1.0)) / (n0 * (n0 + 1.0))) ** 2)
    if nf == 4:                                   # Riesz-2
        return (1.0 - n / ntp1) ** 2
    raise SystemExit(f"nfilter {nf} is not one of the filters legini builds")


def reference(ntru: int, nlat: int, filt: dict, nodes, weights) -> dict:
    """The tables from scipy and the six factors from the declared formula."""
    pmat, qmat = btr.legendre_table(ntru, nodes)
    modes = btr.triangular_modes(ntru)
    sk = filter_curve(ntru, filt)
    ms = np.array([m for m, _ in modes], dtype=np.float64)
    ns = np.array([n for _, n in modes], dtype=np.float64)
    # legmod.f90 counts n from zero here and reads skgpsp(n+1), so mode (m,n)
    # takes the curve's n-th entry counting from zero.
    f = sk[ns.astype(int)]
    znn1 = np.where(ns > 0, 1.0 / np.maximum(ns * (ns + 1.0), 1.0), 0.0)
    csq = 1.0 - np.asarray(nodes, dtype=np.float64) ** 2
    return {
        "pmat": np.asarray(pmat, dtype=np.float64),
        "qmat": np.asarray(qmat, dtype=np.float64),
        "fsp": f, "fgp": f,
        "fmu": znn1 * ms * f, "fmv": znn1 * f,
        "fmm": ms * f, "fmq": ns * (ns + 1.0) * 0.5 * f,
        "gwdc": np.asarray(weights, dtype=np.float64) / csq,
    }


def weight(tables: dict, name: str, gwd: np.ndarray) -> np.ndarray:
    """One of the eight, as (NCSP, NLAT)."""
    table, factor, per_lat = WEIGHTS[name]
    out = tables[table] * tables[factor][:, None]
    if per_lat == "gwd":
        out = out * gwd[None, :]
    elif per_lat == "gwdc":
        out = out * tables["gwdc"][None, :]
    return out


def compare(got: dict, ref: dict, gwd: np.ndarray) -> dict:
    """Worst difference per weight, scaled by the reference's own largest entry."""
    worst = {}
    for name in NAMES:
        a = weight(got, name, gwd)
        b = weight(ref, name, gwd)
        scale = max(float(np.abs(b).max()), 1.0)
        worst[name] = float(np.abs(a - b).max() / scale)
    return worst


def declared_filter() -> dict:
    """The filter the model is configured to run, from config/planet.yaml.

    Not a filter chosen here: the factors under test carry it, so running one
    the model does not would leave the configured one unchecked.
    """
    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "planet.yaml").read_text())
    model = cfg["model"]
    spec = str(model["physics_filter"])
    by_name = {"none": 0, "cesaro": 1, "exp": 2, "hoskins": 3, "riesz": 4}
    nfilter = 0
    for token in spec.split("|"):
        if token in by_name:
            nfilter = by_name[token]
    return {"nfilter": nfilter,
            "nfilterexp": int(model.get("filter_gamma", 8)),
            "filterkappa": float(model.get("filter_kappa", 8.0)),
            "landhoskn0": float(model.get("landhoskn0", 15.0))}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--res", default="T21",
                    help="the rung to run at; its truncation and latitude count "
                         "come from lib/rungs.py, which is where the ladder is "
                         "written down")
    ap.add_argument("--cesaro", action="store_true",
                    help="run the Cesaro filter instead of the configured one. "
                         "It carries an EXACT zero at n = NTP1, so a "
                         "factorisation that recovers P by dividing qi by the "
                         "filter fails here rather than in production.")
    args = ap.parse_args()

    nlat, _nlon, ntru = rungs.geometry(args.res)
    filt = {"nfilter": 1, "nfilterexp": 8, "filterkappa": 8.0,
            "landhoskn0": 15.0} if args.cesaro else declared_filter()

    # THE BOUND, DECLARED BEFORE ANY OF IT RUNS. legini walks an NTRU-long
    # chain in m and, inside it, an NTRU-long chain in n, each step a
    # subtract-and-scale contributing eps to values of order one, so the
    # compounded error is bounded by the product of the two chain lengths.
    # 4*NTRU^2*eps is the envelope verify_banded_transform.f90 derives and
    # measures for the tables alone; the eight weights carry up to two further
    # multiplies, so this doubles it. It is absolute, scaled by the largest
    # entry of the weight it is applied to, because the filter can be exactly
    # zero and a relative bar there has no meaning.
    eps = float(np.finfo(np.float64).eps)
    tol = 8 * ntru * ntru * eps

    nodes, weights = gauss_legendre(nlat)
    gwd = np.asarray(weights, dtype=np.float64)

    shape = ("Cesaro, carrying an exact zero at n = NTP1" if args.cesaro
             else "as config/planet.yaml declares it")
    print(f"{args.res}: NTRU {ntru}, {nlat} latitudes, "
          f"nfilter {filt['nfilter']} ({shape})")
    print(f"declared before the arms ran: 8*NTRU^2*eps = {tol:.3e}, absolute, "
          f"scaled by each weight's own largest entry")
    print("subject:   legini, compiled verbatim out of legmod.f90")
    print("reference: scipy's sph_legendre_p_all, and the filter's declared formula")
    print()

    text = SRC.read_text()
    legini = extract(text, "legini")

    failures = 0

    missing = pairings(text)
    total = sum(len(v) for v in PAIRING_EVIDENCE.values())
    if missing:
        failures += 1
        print("  [ FAIL ] the source does not pair the tables and factors this "
              "gate assumes:")
        for line in missing:
            print(f"           {line}")
    else:
        print(f"  [  ok  ] all {total} pairings read back out of legmod.f90's "
              f"transform routines")
    print()

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        got = run_legini(legini, ntru, nlat, filt, nodes, weights, work, "shipped")
        ref = reference(ntru, nlat, filt, nodes, weights)

        print("the six per-mode factors, against the declared formula")
        for name in FACTORS:
            scale = max(float(np.abs(ref[name]).max()), 1.0)
            d = float(np.abs(got[name] - ref[name]).max() / scale)
            state = "  ok  " if d <= tol else " FAIL "
            if d > tol:
                failures += 1
            print(f"  [{state}] {name}  {d:.3e}")
        print()

        print("the eight weights, each formed from both sides")
        worst = compare(got, ref, gwd)
        for name in NAMES:
            table, factor, per_lat = WEIGHTS[name]
            d = worst[name]
            state = "  ok  " if d <= tol else " FAIL "
            if d > tol:
                failures += 1
            print(f"  [{state}] {name} = {table} * {factor} * {per_lat or '1'}"
                  f"   {d:.3e}")
        print()

        print("the negative controls, every one of which must be caught")
        for name, what, old, new in CONTROLS:
            broken = perturb(legini, name, old, new)
            arm = run_legini(broken, ntru, nlat, filt, nodes, weights,
                             work, f"ctl_{name}")
            bad = compare(arm, ref, gwd)
            if max(bad.values()) > tol:
                hit = [k for k in NAMES if bad[k] > tol]
                print(f"  [  ok  ] {name}: {what} -- caught on {' '.join(hit)}")
            else:
                failures += 1
                print(f"  [ FAIL ] {name}: {what} -- NOT caught, so this gate "
                      f"cannot see it")

    print()
    print("stored: 2 matrices of NCSP x NLPP, 6 vectors of NCSP, 2 scalars a latitude")
    print("against 8 matrices of NCSP x NLPP -- a quarter of the resident traffic")
    print(f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
