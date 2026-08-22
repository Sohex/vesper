#!/usr/bin/env python
"""What the model's hot loops cycle through one die's L3, by component.

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model. Nothing
it produces is about the simulated planet.

WHY THIS IS A SCRIPT AND NOT A TABLE IN PROSE. The per-die working set decides
whether a change is worth taking -- `docs/src/reference/environment.md` sets the
target at 32 MB, which is CCD1 here and is what a part without stacked cache is
likely to have or less. A number that has to be recomputed by hand every time a
resolution or a thread count moves is a number that will be quoted stale, and
one component of it will be quoted as the whole. Both of those happened while
the rule was being written.

WHAT IS COUNTED. The arrays the dynamical core cycles every timestep: the
Legendre weights, the Fourier fields the direct transform holds live, the
spectral state, the tendency partials, and mkdheat's scratch. Each is counted
with its STORAGE CLASS, because that is where the arithmetic goes wrong:
threadprivate means one copy a thread and so `threads_per_die` copies on a die,
while shared means one copy however many threads read it.

WHAT IS NOT COUNTED, and why the total is a floor rather than an estimate: the
physics modules' own grid arrays and their several hundred compiler temporaries.
They are real and they are transient, they do not persist across a timestep the
way these do, and counting them properly means counting what the compiler
invents. This number is what the core cannot avoid touching.
"""
import argparse
import sys

# (NTRU, NLAT, NLON) per spectral truncation, as compile.sh resolves them.
RESOLUTIONS = {
    "T21":  (21,   32,   64),
    "T31":  (31,   48,   96),
    "T42":  (42,   64,  128),
    "T85":  (85,  128,  256),
    "T127": (127, 192,  384),
    "T170": (170, 256,  512),
    "T213": (213, 320,  640),
}

CEILING_MB = 32.0   # docs/src/reference/environment.md


def budget(res: str, nthreads: int, per_die: int, nlev: int = 10):
    ntru, nlat, nlon = RESOLUTIONS[res]
    ncsp = (ntru + 1) * (ntru + 2) // 2
    nrsp = (ntru + 1) * (ntru + 2)
    nlpp = nlat // nthreads
    nspp = -(-nrsp // nthreads)
    nesp = nspp * nthreads
    nhor = nlon * nlpp
    B = 8  # -fdefault-real-8

    # name, elements, storage class. Three classes and the third is the one that
    # is easy to get wrong:
    #   private  one copy a thread, so per_die copies on a die
    #   shared   one copy, and every thread reads all of it
    #   sliced   one copy, but a thread only ever touches its own NSPP rows of
    #            it, so a die touches per_die/nthreads of the array. The
    #            reduction partials are this: thread t writes slot t and reads
    #            rows [t*NSPP,(t+1)*NSPP) of EVERY slot, never the rest.
    rows = [
        ("Legendre weights pmat, qmat",      2 * ncsp * nlpp,            "private"),
        ("per-mode factors fsp..fmq",        6 * ncsp,                   "shared"),
        ("mktend Fourier fields, six",       6 * 2 * (nlon // 2) * nlpp * nlev, "private"),
        ("dv2uv pre-scaled fields, four",    4 * 2 * ncsp,               "private"),
        ("spectral state sd,st,sz,sq,sr",    5 * nesp * nlev,            "shared"),
        ("spectral state sp",                nesp,                       "shared"),
        ("tendency partials zpsd..zpsq",     4 * nesp * nlev * nthreads, "sliced"),
        ("mkdheat gathered zhd,zhz,zhq,zhe", 4 * nesp * nlev,            "shared"),
        ("mkdheat partials zhf1,zhf2,zhef",  3 * nesp * nlev * nthreads, "sliced"),
        ("mkdheat grid scratch, seven",      7 * nhor * nlev,            "private"),
    ]

    out, total = [], 0.0
    for name, elems, cls in rows:
        if cls == "private":
            copies, frac = per_die, 1.0
        elif cls == "sliced":
            copies, frac = 1, per_die / nthreads
        else:
            copies, frac = 1, 1.0
        mb = elems * B * copies * frac / 1e6
        total += mb
        out.append((name, cls, copies, mb))
    return out, total, dict(NCSP=ncsp, NESP=nesp, NLPP=nlpp, NSPP=nspp, NHOR=nhor)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--res", default="T170", choices=sorted(RESOLUTIONS))
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--per-die", type=int, default=8,
                    help="threads sharing one L3. 8 on a 7950X3D CCD.")
    ap.add_argument("--levels", type=int, default=10)
    ap.add_argument("--ceiling", type=float, default=CEILING_MB)
    a = ap.parse_args()

    if a.threads <= 0 or RESOLUTIONS[a.res][1] % a.threads:
        print(f"{a.res} has NLAT={RESOLUTIONS[a.res][1]}, which {a.threads} does not divide",
              file=sys.stderr)
        return 2

    rows, total, dims = budget(a.res, a.threads, a.per_die, a.levels)
    print(f"{a.res}, {a.threads} threads, {a.per_die} sharing a die, {a.levels} levels")
    print("  " + "  ".join(f"{k}={v}" for k, v in dims.items()))
    print()
    print(f"  {'component':36s} {'class':9s} {'copies':>6s} {'MB/die':>8s}")
    for name, cls, copies, mb in sorted(rows, key=lambda r: -r[3]):
        print(f"  {name:36s} {cls:9s} {copies:6d} {mb:8.2f}")
    print(f"  {'':36s} {'':9s} {'':>6s} {'-'*8}")
    print(f"  {'total the core cannot avoid':36s} {'':9s} {'':>6s} {total:8.2f}")
    print()
    print(f"  target {a.ceiling:.0f} MB a die: {'MET' if total <= a.ceiling else 'OVER by %.2f MB (%.1fx)' % (total - a.ceiling, total / a.ceiling)}")
    print()
    print("  Physics grid arrays and compiler temporaries are NOT counted, so this")
    print("  is a floor. A component's own number is never the budget: check the")
    print("  total, and check the copies -- threadprivate means one a thread.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
