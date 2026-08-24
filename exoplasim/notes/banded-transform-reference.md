# The independent second side for the threaded grid bands

Worldbuilding frame: measurements on the Vesper project's climate model. The
fields and tables here are model arrays; nothing is about the real world.

`verify_threaded_numerics.sh` was written as the threaded build against an
independent build. world-38b left one build, and the gate refused rather than
comparing the thing under test against itself. world-d5l took the decision to
build a standalone independent driver instead of keeping only the control arm
or retiring the gate. This records what that driver measures, what it is
independent of, what it is not, and the numbers it returns across the ladder.

Measured on 2026-08-24, gfortran 16.2.1, the `checked` profile of
`config/planet.yaml`, on this host.

## What is independent of the model, and what is not

The subject is the model's own analysis chain, compiled out of the model tree
at NPRO threads with OMPSHARED: `inigau`, the `mpscdn` scatter, `legini`,
`assoc_grid`, `fc2sp`, and `mpsumsc` followed by `mpgallsp`. The reference is a
file written before any of that runs.

INDEPENDENT. The associated Legendre functions come from scipy's
`sph_legendre_p_all`, which shares no code and no recurrence with `legini`. The
mapping between the two normalisations is derived and not fitted:
`pmat = sqrt(2*pi) * Pbar_n^m` and `qmat = (1 - mu^2) d/dmu` of the same, the
first because the model puts P(0,0) at sqrt(1/2) where the fully normalised
harmonic puts it at 1/sqrt(4*pi), the second because the model's wind variables
carry a factor of cos(phi).

NOT INDEPENDENT, and this is the judgement call the design turns on. The
Gaussian nodes and weights come from `gauss_weight_reference.py`, which runs the
SAME textbook construction `inigau` runs -- Newton on P_n by the three-term
recurrence, then w = 2/((1-x^2) P'^2) -- in x86 longdouble rather than float64.
That adjudicates rounding by eight orders and does not adjudicate method: a
shared misreading of what Gauss-Legendre quadrature is would be invisible to
both sides together.

What closes that hole is not a second implementation but two identities, checked
at the rung the reference is being built for, before anything is written:

| identity | exact value | measured, T21 | measured, T170 |
| --- | --- | --- | --- |
| sum of the weights | 2 | off by 2.168e-19 | off by 2.168e-19 |
| sum of w P_k, 1 <= k <= 2*NLAT-1 | 0 | worst 3.015e-19 | worst 2.202e-19 |

A rule on NLAT nodes satisfying both is the Gauss-Legendre rule, because that
rule is unique. So the method is pinned by an identity even though the two
implementations are related, and the residuals are at longdouble eps.

## Why the right answer is a right answer

The driver is handed Fourier coefficients the REFERENCE synthesised from a known
spectral vector, and the model's analysis of them has to return that vector,
because

    sum over l of gwd(l) * pmat(w,l) * pmat(w',l) = delta(w,w')

for two modes sharing an m. That is orthonormality under the quadrature: a
property of the functions rather than of either implementation.

Handing the model coefficients IT synthesised would make this a round trip, and
a round trip holds on a wrong weight matrix because synthesis and analysis drop
the same factor together. That is exactly the shape of the weight pre-scaling
defect the retired MPI arm once caught, so the synthesis is done on the
reference side and nowhere else.

The reference's own error on that identity, which is the instrument's error bar:

| rung | orthonormality residual of the reference table |
| --- | --- |
| T21 | 2.109e-15 |
| T170 | 1.497e-14 |

Both are two to three orders inside the bound the analysis arm runs at, so the
comparison is limited by the model and not by the reference.

## The bounds, and where each comes from

None was chosen after seeing what the model produced.

- Band placement has NO bound. Which global latitude a thread holds is an
  integer, so the check is which reference node each of the thread's rows
  matches most closely and the answer has to be the row the decomposition says
  it owns.
- The weights are compared at `verify_gauss_weights.sh`'s own bar, 1e-12
  relative, which that gate derives as an order above the 2*n*eps floor of
  evaluating P_n by the three-term recurrence. It is transplanted rather than
  re-derived on purpose: two gates with two bars on one quantity is how a
  quantity ends up with two answers. The node VALUES are left to that gate
  entirely.
- gwd/cos^2 carries that bar plus eps/(1-mu^2) at the pole-most node, because
  `legini` forms cos^2 as 1 - sid*sid and at the top of the ladder sid is
  0.99997, where that subtraction loses four digits.
- The Legendre tables are bounded at 4*NTRU^2*eps. The FORM is the recurrence:
  an NTRU-long chain in m with an NTRU-long chain in n inside it, each step a
  subtract-and-scale contributing eps to values of order one. The COEFFICIENT is
  an envelope over the declared ladder rather than a fit, measured below.
- The analysis is bounded at 8*NLAT*eps*(1+Pmax) with Pmax the largest |pmat| in
  the reference table. The sum is NLAT terms each at most Pmax * gwd(l) * max|fc|;
  the reference normalises max|fc| to one and the weights sum to two, so the
  summands total at most 2*Pmax and NLAT roundings of that is 2*NLAT*eps*Pmax.
  The remaining term is the quadrature's own relative error on a result of order
  one, and the whole is doubled twice.
- The repeat arm has no bound either: four repeats in one process must agree bit
  for bit. A reduction summing in arrival order gives a different answer run to
  run rather than a wrong one, which a tolerance cannot see.

## `legini`'s tables against the reference, across the ladder

`legini`'s own recurrence reproduced in Python and compared with the scipy
table, so this is a property of the algorithm and not of the compiled build.
`dP` is absolute; `dQ` is absolute and Q grows with NTRU, which is why the gate
scales it by the table's own largest magnitude.

| rung | NTRU | worst dP | dP / (NTRU^2 eps) | worst dQ |
| --- | --- | --- | --- | --- |
| T21 | 21 | 1.554e-14 | 0.16 | 9.193e-14 |
| T31 | 31 | 4.352e-14 | 0.20 | 1.457e-13 |
| T42 | 42 | 7.416e-14 | 0.19 | 3.140e-13 |
| T63 | 63 | 1.861e-13 | 0.21 | 6.839e-13 |
| T85 | 85 | 2.585e-13 | 0.16 | 1.178e-12 |
| T106 | 106 | 2.922e-13 | 0.12 | 2.172e-12 |
| T127 | 127 | 1.219e-12 | 0.34 | 2.658e-12 |
| T170 | 170 | 2.220e-12 | 0.35 | 4.833e-12 |

The ratio to NTRU^2*eps stays between 0.12 and 0.35 over the whole ladder, so
4*NTRU^2*eps clears it by between eleven and thirty-three times. The ratio is
still climbing at T170, so the envelope holds over `lib/rungs.py` and is not
argued beyond it: a rung above T170 needs it measured again before the bound is
used there.

## What the compiled gate returns

Every rung below passes in all three driving cases. T21 was run on four threads
and the rest on eight; T42 was run on 2, 4, 8 and 16 and returned the same table
differences to the last digit at each, so the tables do not depend on the thread
count.

| rung | weights | gw/cos2 | P | Q, scaled | analysis, worst of the three cases | analysis bound |
| --- | --- | --- | --- | --- | --- | --- |
| T21 | 2.60e-15 | 2.59e-15 | 1.11e-14 | 5.31e-15 | 2.53e-16 | 2.17e-13 |
| T31 | 2.379e-14 | 2.387e-14 | 3.775e-14 | 5.555e-15 | 6.453e-16 | 3.789e-13 |
| T42 | 2.785e-14 | 2.785e-14 | 1.449e-13 | 1.477e-14 | 7.633e-16 | 5.652e-13 |
| T63 | 6.69e-14 | 6.70e-14 | 1.02e-13 | 1.76e-14 | 7.77e-16 | 9.99e-13 |
| T85 | 2.81e-13 | 2.81e-13 | 1.60e-13 | 1.79e-14 | 1.17e-15 | 1.50e-12 |
| T106 | 2.75e-13 | 2.75e-13 | 1.32e-13 | 2.56e-14 | 8.47e-16 | 2.06e-12 |
| T127 | 3.84e-13 | 3.84e-13 | 5.30e-13 | 3.65e-14 | 9.84e-16 | 2.68e-12 |
| T170 | 2.53e-13 | 2.53e-13 | 8.07e-13 | 3.93e-14 | 1.47e-15 | 4.06e-12 |

The tightest margin anywhere is the weights at T127, 2.6 times inside a bar
`verify_gauss_weights.sh` owns, so that gate fails before this one does. The
analysis arm sits three orders inside its bound at every rung.

## The controls, and what each proves

| control | what it breaks | rejected by | measured at T21 |
| --- | --- | --- | --- |
| noweight | the quadrature weight dropped from `fc2sp` | the analysis arm | 5.06 against a 2.17e-13 bound |
| slideband | every band slid one latitude row, so neighbours share one | the analysis arm | 0.265 against the same bound |
| samelats | the scatter hands every thread the master's latitudes | the placement arm | every non-master row wrong |

`noweight` and `slideband` leave the table arm passing and `samelats` leaves the
analysis arm failing as well, which is the two arms behaving as their split
predicts: a thread's weight matrices come from the scatter, and its Fourier
rows come from the band pointer, and only a control that moves the scatter moves
both.

## What a pass does not cover

- Anything that needs a timestep. The driver does not integrate the model, so
  it is silent on `mkdheat` and on every physics term.
- `dv2uv` under the bands. It is a synthesis routine with no partial sum and no
  reduction, and its only band dependence is the `pmat` and `qmat` rows for the
  thread's own latitudes, which the table arm checks directly. Its loop and its
  mirror combination are checked by `verify_symmetric_transform.py` against a
  scipy table. What is not checked anywhere is `dv2uv`'s planetary vorticity
  term under the bands, and pinning that needs the vector convention
  `probe_shtns_vector_conventions.f90` was written to discover.
- The growth of a last-bit difference with run length. That needs two model arms
  differing at last-bit scale and there is one build. `verify_shtns_model.sh`
  carries the same curve on NSHTNS=0 against NSHTNS=1, which is still a live
  pair.
