# The radiation's masks and four-byte reals compound, and multiplicatively

*Worldbuilding frame: this is a COMPUTE measurement of the Vesper project's
climate model on this desktop. Nothing here is about the simulated planet.
Measured 2026-08-25 at T21, NLEV 10, sixteen threads, against `radmod.f90` as it
now stands.*

CLIM-84 established that every non-integer `**` in the radiation sits inside a
`where` block, so none of them reaches libmvec's eight-wide `pow`;
`radiation-scheme-price.md` carries that finding and what each transcendental
costs in instructions. CLIM-59 established that four-byte reals are worth tens of
per cent of wall clock, and that the mechanism is SIMD LANES rather than memory
bandwidth: this host is znver4 with AVX-512 and the declared flag line carries
`-march=znver4`, so a vector register holds sixteen floats against eight doubles
and everything that vectorises does twice the work per instruction.

Neither measurement was taken with the other in place, and the two are not
independent: removing a mask moves work out of the bucket four bytes cannot help
and into the bucket where four bytes doubles the lanes. **WORLD-43RK is whether
they compound.** They do. Their gains multiply, and the mask removal is worth
four times less in TIME than its instruction count says.

## Four arms, because two would measure one effect at an arbitrary setting

The arms are the two factors crossed: the declared eight-byte precision and four
bytes, against `lwr`'s absorptances masked and unmasked. A two-arm comparison
answers "what is the mask worth" only for whichever precision the arms happened
to be built at, and that is the question rather than a detail of it.

The four arms differ by the executable and by nothing else. One bed, one
namelist, one thread count, one launcher; the schedule is interleaved and the
starting arm rotates by round, so over any four consecutive rounds each arm
occupies each position once. `exoplasim/scripts/bench_ab.py` carries the
schedule and the analysis.

The bed is a cold start from the declared `cold_start_seed` at T21, output
streams off, built by `make_profile_bed.py` from a T21 run's inputs. A cold bed
rather than a warm one because a four-byte build cannot read an eight-byte
restart, and an identical cold start is the comparison
`notes/audits/single-precision-spin-up.md` already specifies for a precision
arm.

## The arm: what the unmasked form is, and why it is admissible

`lwr`'s transmissivity loop runs over NLEV(NLEV+1)/2 level PAIRS, and four of
its absorptances are written as a `where`/`elsewhere` that selects between a
fractional-power fit below a path amount and a logarithmic fit above it: the
water vapour 6.3 um band, the CO2 band, the Boer water vapour transmissivity at
the CO2 overlap, and ozone. The arm replaces each with one unmasked array
assignment, the branch chosen arithmetically by `0.5+SIGN(0.5,threshold-path)`.
`exoplasim/scripts/patch_unmasked_absorptance.py` applies and reverts it, and
carries both forms verbatim so neither can drift.

**Both branches are then evaluated on every lane, and that is safe here.** A
`where` selects the assignment and not the evaluation, and the declared flag
line carries `-ffpe-trap=invalid,zero,overflow`, so an unmasked form is
admissible only if every lane is in domain. Every argument in these eight
expressions already carries an unconditional floor INSIDE it -- `max(0.,.)`
before each fractional power, so the base is never negative and `pow(0.,a)` is
`0.` for positive `a`; `max(1.E-30,.)` before each logarithm, so the argument is
never zero and never negative. That was true before the arm existed and is what
makes the arm available. `notes/audits/masked-where-blocks.md` argues the class.

**The restructuring itself is exact.** The selector is exactly `1.` or exactly
`0.`, so one term is the chosen branch unchanged and the other is `0.` times a
finite number. `exoplasim/scripts/probe_unmasked_absorptance.f90` evaluates all
four both ways over path amounts spanning fourteen decades plus the exact
thresholds and exact zero, under the model's own trap, and with vectorisation
disabled the two forms agree on every element of 65536, bit for bit, in both
precisions. No lane raised invalid, divide-by-zero or overflow in any run of it.

**What is not bit identical is the library, and only the library.** With
vectorisation enabled the same probe reports a maximum absolute difference of
2.2e-16 in eight-byte reals and 1.2e-7 in four-byte, on absorptances that are
dimensionless fractions the routine bounds into [0,1] two lines later. That is
libmvec's looser ULP bound and nothing else. It still means every result moves,
and rule 4 still applies: it is a model source change, so every published binary
goes stale.

**The compiled object confirms the vectorisation.** Reading `radmod.f90.o` out
of the two builds, by relocation count:

| symbol | masked | unmasked |
| --- | ---: | ---: |
| scalar `pow` | 87 | 64 |
| `_ZGVeN8vv_pow` | 0 | 16 |
| scalar `log10` | 25 | 9 |
| `_ZGVeN8v_log10` | 0 | 12 |
| scalar `log` | 98 | 91 |
| `_ZGVeN8v_log` | 44 | 48 |

The masked column reproduces CLIM-84's counts exactly. The scalar `pow` left in
the unmasked column is the peel and epilogue of the vectorised loops plus the
eighteen `**` sites this arm does not touch.

## Three instruments, because on this host one of them cannot answer

Every arm was run three times over: retired instructions, unhalted cycles, and
the wall clock, all through the same interleaved schedule. The reason is that
they say different things and only one of them survives a shared desktop.

| instrument | worst self-scatter | what it says |
| --- | ---: | --- |
| `instructions:u` | 0.0017% | how much work there is |
| `cycles:u` | 4.09% | what that work costs, at whatever clock |
| wall clock, quiet host | 1.87% | what the run took |
| wall clock, contended host | 11.27% | what the host was doing as well |

The counters run under `OMP_WAIT_POLICY=passive`, so an idle thread sleeps
rather than retiring instructions spinning. The instruction count reproduced to
about one part in 10^5 across four rounds while the host's one-minute load fell
from 30.4 to 8.1 underneath it; that is the property it was chosen for, and it
is why `radiation-scheme-price.md` uses it.

**THE WALL CLOCK IS A MEASUREMENT OF THE HOST AS WELL, AND THIS ONE TOOK THREE
ATTEMPTS.** The bed was sized on the reasoning that the smallest effect to be
resolved is the mask removal alone, whose ceiling was the radiation's 6.69 per
cent `pow` share of T170, so a run of about 30 seconds -- forty times what a
200-step bed took -- should put the noise well under it. Sizing was not the
problem.

The first schedule, eight rounds on a 16000-step bed while four other worktrees
were working, gave the eight-byte masked arm 29.92, 32.29, 35.02, 30.28 and
38.24 seconds in its first five rounds, put a single arm at 219.26 seconds
against 22 in every other, and then put all four arms above 200 seconds while
another worktree compiled. It was abandoned. **Before the surviving rounds were
seen and before any gain was computed**, the rule was fixed that a round in which
any arm exceeds three times that arm's own median is a machine event rather than
a measurement, and is dropped WHOLE -- every arm of it together, so the pairing
survives. `bench_ab.py --outlier-factor` carries it and always prints what it
removed. The second schedule, on an 8000-step bed, dropped two of its eight
rounds under that rule and still scattered 11.3 per cent in its noisiest arm,
which is wider than the effect; `bench_ab`'s own floor refused that arm, and
correctly.

The third schedule is the one reported below. It was taken under a claim on a
host at load 1.49 with nothing else running, and its worst arm scattered 1.87
per cent against a smallest effect of 4.25. Nothing was dropped, and every arm
was faster than the base in every round. **The difference between the second
schedule and the third is entirely the machine**: same binaries, same bed, same
schedule, 11.3 per cent scatter against 1.9.

WHAT THE COUNTERS CANNOT SAY: neither is a duration, and cycles are taken at
whatever frequency the part chose. They agree with the quiet clock on every arm
to within its round-to-round range, which is the check that all three are
measuring the same thing.

## The four arms

T21 NLEV 10 on sixteen threads. Counters over 4000 timesteps, four rounds
(instructions) and six rounds (cycles), host load 8 to 30 with one cGENIE cost
probe and another worktree's builds resident throughout. Wall clock over 16000
timesteps, six rounds, host load 1.49 at the start and nothing else running. No
round dropped in any of the three.

| arm | instructions | cycles | IPC | seconds |
| --- | ---: | ---: | ---: | ---: |
| eight-byte, masked | 6.2194e11 | 3.8842e11 | 1.601 | 30.23 |
| eight-byte, unmasked | 4.7398e11 | 3.6681e11 | 1.292 | 28.99 |
| four-byte, masked | 5.0153e11 | 2.8328e11 | 1.770 | 22.50 |
| four-byte, unmasked | 4.0337e11 | 2.7263e11 | 1.480 | 21.85 |

Paired per-round gain against the eight-byte masked arm, with the range over
rounds:

| arm | instructions | cycles | seconds |
| --- | ---: | ---: | ---: |
| eight-byte, unmasked | -23.79% | -5.70% [-4.59, -6.32] | -4.25% [-3.68, -4.30] |
| four-byte, masked | -19.36% | -26.93% [-24.03, -28.36] | -25.53% [-24.67, -26.17] |
| four-byte, unmasked | -35.14% | -29.75% [-29.16, -30.84] | -27.63% [-27.33, -28.99] |

Per-arm self-scatter on the quiet clock schedule was 0.85, 0.37, 1.31 and 1.87
per cent, in the order of the table. Every arm was faster than the base in every
round of every instrument. Each arm's restart sha differs from every other,
which is the numerics change stated separately from the cost and never folded
into it.

**The four-byte arm reproduces CLIM-59 independently.** 25.53 per cent on a
quiet clock and 26.93 in cycles, against the 27.5 per cent that row measured at
T21 on a different bed. Different bed, different day, different instrument, and
the three agree to two points.

## The mask removal is worth four times less in time than in instructions

Nearly a quarter of the whole model's retired instructions at this rung were in
four lines of `lwr`, and removing them bought 5.70 per cent of its cycles and
4.25 per cent of its wall clock.
The arithmetic on the instruction side corroborates the size of the work: the
base retires 75,900 instructions per gridpoint per timestep, `radstep` runs on
every timestep, and four masked absorptances over 55 level pairs at 116
instructions for a masked scalar `pow` is 25,500 per gridpoint per longwave call
before the logarithms are counted.

What did not follow is the cost. **IPC falls from 1.601 to 1.292 when the masks
go**, and by the same proportion at four bytes. The masked scalar `pow` calls
were largely hidden -- they are independent per element, the loop around them has
no carried dependence, and an out-of-order core with a wide window absorbs a
great deal of that. The vector `pow` that replaces eight of them retires a
quarter of the instructions and has neither the latency nor the throughput to be
eight times faster.

**This corrects an inference, not a measurement.** CLIM-84's per-element
instruction counts stand -- 116 for a masked scalar `pow` against 15.2 for a
vector one is what those two forms retire. What does not stand is reading a
factor of 7.6 in instructions as a factor of 7.6 in cost. Measured in place, the
whole-model instruction saving over-states the cycle saving by 4.2 and the wall
clock by 5.6.

## They compound, and the interaction is not resolved against multiplying

| | instructions | cycles | seconds |
| --- | ---: | ---: | ---: |
| unmasking, at eight bytes | -23.79% | -5.70% | -4.25% |
| unmasking, at four bytes | -19.57% | -4.12% | -3.23% |
| four bytes, masked | -19.36% | -26.93% | -25.53% |
| four bytes, unmasked | -14.90% | -25.70% | -24.63% |
| both together | -35.14% | -29.75% | -27.63% |
| multiplying the two single-factor gains predicts | -38.54% | -30.89% | -28.53% |
| **excess over that prediction** | **+3.40** | **+1.14 [-1.23, +2.22]** | **+0.34 [-0.20, +1.65]** |

Two things follow, and the first is the answer to the row.

**Each change is worth having with the other in place.** Nothing collapses:
unmasking is worth 3.23 per cent of a four-byte model's run and four bytes is
worth 24.63 per cent of an unmasked one, and together they retire a third of the
instructions and take 27.63 per cent off the clock.

**Whether they interfere is settled in instructions and NOT settled in cost.**
The instruction excess is 3.40 points with a spread of 0.0013 across rounds, so
they are sub-multiplicative there by three orders of magnitude more than the
noise, and the reason is plain: both act on the same work, because unmasking's
whole effect is to hand the radiation's `pow` to the vectoriser, where four
bytes has already halved what a lane costs. In cycles the excess is 1.14 points
against a round-to-round spread of 3.45, and on the quiet clock 0.34 against
1.85. **In cost the two are consistent with multiplying, and a departure from it
is not resolved even by the quiet bed**: the interaction is a third of a point
and the bed resolves about two. Resolving it would need a bed longer by roughly
the square of that ratio, and nothing downstream turns on the answer, so it is
recorded as unresolved rather than pursued.

The hypothesis this row was filed on -- that unmasking moves work into the bucket
four bytes doubles, so the two should EXCEED either alone -- is confirmed in its
conclusion and wrong in its arithmetic. They overlap rather than stack.

## What this does not settle

- **The rung.** One bed at one rung. The mix shifts with resolution: the
  spectral transforms take a larger share at T85 against a physics column count
  that grows more slowly, so the radiation's share, and therefore unmasking's,
  falls as the rung rises, while CLIM-59 measured four bytes moving the other
  way, 27.5 per cent at T21 against 39.2 at T85. Both single-factor numbers are
  T21 numbers and so is their interaction.
- **The other eighteen `**` sites.** This arm touches only the four
  formula-selector masks inside `lwr`'s level-pair loop. `swr`'s `where(losun)`
  blocks are work-elision masks rather than formula selectors: unmasking those
  computes the shortwave on the night side, which is more work and not less, and
  one of them raises `zmu0**1.7` on a negative cosine, which the declared trap
  would catch. They are a different question and were left alone.
- **The state the branch mix comes from.** The bed is a cold start, so the
  masked arm's cost depends on how many lanes take the power branch, and a
  spun-up atmosphere holds a different mix. The unmasked arm evaluates both
  branches and does not care, so the measured gain is an upper bound on what a
  spun-up bed would give.
- **Whether four bytes is SOUND**, which is CLIM-59's own question and is
  untouched here. This note prices four bytes; it says nothing about whether a
  four-byte spin-up lands in the right place.

## The verdict, which does not change the sequencing

**CLIM-84's ordering stands and is not softened.** The change is not bit
identical, `merge()` was measured and does not recover the vectorisation at 107
instructions per element, and CLIM-61 is open: if a band-resolved scheme is
adopted, `lwr` and `swr` are deleted outright and this work is thrown away. The
measurement makes the case for taking it WEAKER rather than stronger -- 4.25 per
cent of wall clock at T21, falling as the rung rises, for a change that moves
every result. **Do not adopt it, and the reason is now a number rather than an
ordering argument alone.**

What HAS changed is an input to CLIM-61's own decision, and it runs against
deciding that row on cost. `radiation-scheme-price.md` priced the candidate as at
worst comparable and probably several times cheaper on transcendental work, and
rested that on the broadband kernel being a masked scalar `pow` at 116
instructions against a k-term scheme's vectorisable `exp` at 5.1. Both halves of
that ratio are instruction counts, and this note measures what an instruction
count is worth in place: a 23.79 per cent instruction saving bought 4.25 per cent
of wall clock, an over-statement of 5.6. **So the scheme comparison in that note
should not be read as a cost comparison in either direction, and CLIM-61 should
be settled on its structural argument -- three known holes in the re-weighted
broadband scheme and seven manual per-star derivations a band-resolved one
performs as configuration -- rather than on a price neither scheme has been
measured at.**

## The evidence

`exoplasim/analysis/radiation_mask_arms_instructions_t21.json`,
`..._cycles_t21.json` and `..._wallclock_t21.json` hold every round of the three
schedules reported above, with each arm's executable, launcher, namelist
override, restart sha, per-round times and the factorial breakdown.
`..._wallclock_contended_t21.json` is the second, refused schedule, kept because
what a contended clock looks like is worth being able to recognise.

## Reproduce

    python scripts/machine.py --claim "world-43rk" --minutes 45 --who <id>
    python exoplasim/scripts/make_profile_bed.py --from-run <a T21 run> \
        --dest exoplasim/bench/43rk/bed_t21 --cold --steps 4000 \
        --binary most_plasim_t21_l10_p16.x --binary-dir <that run>

    # masked arms
    python exoplasim/scripts/patch_unmasked_absorptance.py --revert
    python exoplasim/scripts/build_model.py --res T21 --ranks 16 --precision 8 --no-publish --print-path
    python exoplasim/scripts/build_model.py --res T21 --ranks 16 --precision 4 --no-publish --print-path
    # unmasked arms. Copy each executable out before the next build: the build
    # tag names every input that changes a byte of the executable EXCEPT the
    # state of the model source, so these four share two directories.
    python exoplasim/scripts/patch_unmasked_absorptance.py --apply
    python exoplasim/scripts/build_model.py --res T21 --ranks 16 --precision 8 --no-publish --print-path
    python exoplasim/scripts/build_model.py --res T21 --ranks 16 --precision 4 --no-publish --print-path
    python exoplasim/scripts/patch_unmasked_absorptance.py --revert

    python exoplasim/scripts/bench_ab.py --bed exoplasim/bench/43rk/bed_t21 \
        --threads 16 --rounds 6 --factorial --counter cycles:u \
        --nl N_RUN_STEPS=4000 \
        --arm "fp64 masked=<x>" --arm "fp64 vector=<x>" \
        --arm "fp32 masked=<x>" --arm "fp32 vector=<x>"

`--counter instructions:u` for the work, and no `--counter` for the clock. The
equivalence half needs no bed, no binary and no claim:

    gfortran -O2 -cpp -ffpe-trap=invalid,zero,overflow -ffpe-summary=none \
        -march=znver4 -funroll-loops -g -fdefault-real-8 -fno-tree-vectorize \
        exoplasim/scripts/probe_unmasked_absorptance.f90 -o probe && ./probe
