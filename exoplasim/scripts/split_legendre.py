"""Split the legendre bucket by routine: how much of it is a plain transform.

    python exoplasim/scripts/split_legendre.py

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
desktop. Nothing here is about the simulated planet.

The share alone says how big the prize is. This says what shape the work is,
which is the difference between a swap and a rewrite. Three classes, and the
boundary that matters is the second one, not the first:

  plain   scalar analysis and synthesis. One library call each.
  vector  dv2uv and uv2dv are (divergence, vorticity) <-> (u, v), which is a
          spheroidal/toroidal transform and a standard primitive -- NOT a
          fusion problem, though PlaSim folds a cos(phi) scaling into the
          weights that would have to move into a separate pass.
  fused   mktend and qtend genuinely combine physics with the Legendre sum,
          mixing four weight matrices (P*gw, Q*gw/cos2, P*gw*m/cos2 and
          P*gw*n(n+1)/2/cos2) to make three spectral tendencies at once. Each
          decomposes into two or three library calls plus rescaling passes.

`sp2fl` and `invlega`/`invlegd` read as 0.00% because they are thin loops whose
samples land in the `sp2fc` they call.
"""
import re
import subprocess
import sys
from pathlib import Path

FUSED = {"mktend_", "qtend_"}
VECTOR = {"dv2uv_", "uv2dv_"}
PLAIN = {"fc2sp_", "sp2fc_", "sp2fcdmu_", "sp2fl_", "sp3fc_", "invlega_", "invlegd_"}
FFT = {"gp2fc_", "fc2gp_", "dfft2_", "dfft3_", "dfft4_", "dfft8_",
       "ifft2_", "ifft3_", "ifft4_", "ifft8_"}

BEDS = sys.argv[1:] or sorted(
    d.name for d in Path("exoplasim/bench").iterdir()
    if d.is_dir() and d.name.startswith("bed_") and any(d.glob("perf_*"))
)

for bed in BEDS:
    tot: dict[str, float] = {}
    n = 0
    for f in sorted(Path("exoplasim/bench", bed).glob("perf_*/perf.rank??.data")):
        out = subprocess.run(
            ["perf", "report", "--stdio", "--quiet", "--no-children",
             "--percent-limit", "0", "-F", "overhead,symbol", "-i", str(f)],
            capture_output=True, text=True).stdout
        n += 1
        for line in out.splitlines():
            m = re.match(r"\s*(\d+\.\d+)%\s+\[[.k]\]\s+(.+?)\s*$", line)
            if m and m.group(2) in FUSED | VECTOR | PLAIN | FFT:
                tot[m.group(2)] = tot.get(m.group(2), 0.0) + float(m.group(1))

    def kind_of(s):
        return ("fused" if s in FUSED else "vector" if s in VECTOR
                else "fft" if s in FFT else "plain")

    print(f"\n=== {bed}  ({n} rank-passes)  % of all samples ===")
    for s, v in sorted(tot.items(), key=lambda kv: -kv[1]):
        if v / n >= 0.005:
            print(f"  {v / n:6.2f}%  {s:12} {kind_of(s)}")
    by = {k: sum(v for s, v in tot.items() if kind_of(s) == k) / n
          for k in ("plain", "vector", "fused", "fft")}
    grand = sum(by.values())
    print(f"  ---- plain {by['plain']:.2f}%  vector {by['vector']:.2f}%  "
          f"fused {by['fused']:.2f}%  fft {by['fft']:.2f}%   total {grand:.2f}%")
    print(f"  maps to a library primitive directly: "
          f"{100 * (by['plain'] + by['vector'] + by['fft']) / grand:.0f}% of it; "
          f"needs decomposing: {100 * by['fused'] / grand:.0f}%")
