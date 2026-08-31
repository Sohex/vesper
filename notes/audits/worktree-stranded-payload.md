# Ignored payload written in a worktree lives only there

Audited 2026-08-30. Vesper is a simulated world; the literature indexed here is
real, and the tree under audit is this project's own tooling.

The question: `references/INDEX.md` cited thirteen PDFs by filename that were not
on disk, eleven of them marked **read** and several carrying numbers transcribed
out of them. Where did they go, and what else goes the same way?

## What is true

**The thirteen rows were real and the files were absent.** Measured 2026-08-30
against `/home/cfutro/docs/world/references/`, which then held 383 PDFs at the
top level and 450 counting subdirectories. None of the thirteen filenames matched
any of them, at the top level or below.

**The mechanism is `scripts/link_worktree.py`'s per-file linking, and it is a
property of the directory's shape.** The set of things to link is derived from
`git ls-files --others --ignored --exclude-standard --directory` in the main
checkout. That listing collapses a WHOLLY-ignored directory into one entry and
lists files INDIVIDUALLY where the parent holds tracked content. So:

| directory | tracked content | what the worktree gets | a file written there afterwards |
| --- | --- | --- | --- |
| `biosphere/runs/` | none | one directory symlink | lands in the main checkout |
| `references/` | `INDEX.md` | a real directory of one link per existing PDF | is real HERE and nowhere else |
| `exoplasim/runs/` | `INDEX.json` | a real directory of one link per existing run | is real HERE and nowhere else |
| `source/` | `README.md` | a real directory of one link per existing build | is real HERE and nowhere else |

The second column is the whole of the difference. `references/` cannot be a
directory symlink, because `INDEX.md` is tracked and a worktree must be able to
edit its own copy.

**The consequence has three steps and no warning at any of them.** A PDF fetched
into a worktree's `references/` is a real file there. It matches `.gitignore`, so
it is never committed. The worktree is removed and it is gone -- while the
`INDEX.md` row written for it in the same session, being tracked, is committed
and survives. The row outlives its own artifact, and what it outlives it by is
the entire distance between **read** and **held**: a row marked read with numbers
in it, whose file is absent, cannot be checked against the artifact, which is the
standard the index exists to hold.

**The same shape puts a climate run and an export payload at the same risk**, and
those are rule 7's durable set rather than a bookkeeping matter. A run generated
under a worktree's `exoplasim/runs/<uuid>/` is real only there; `runs/INDEX.json`
is tracked and would carry the row for it.

**A sweep of the 84 standing worktrees on 2026-08-30 found 21 files in this
state.** None of them belonged to this audit's worktree.

## What is not recoverable

Which sessions wrote which of the thirteen rows, and whether every one of them
was lost this way rather than some other. The payload is untracked, so git
carries no record of a file that was never committed, and the mechanism above is
established from the tooling's behaviour rather than from an observation of the
loss. What is established is that the mechanism exists, that it produces exactly
this outcome, and that it was live on 2026-08-30 in 21 places.

## The thirteen, and how they came back

All thirteen are held again as of 2026-08-30, under the filenames `INDEX.md`
asserts. Eight returned through `paperfetch` from open-access routes; five were
supplied directly after the automated routes refused them, and those five are
recorded below with what refused.

| filename stem | route | quoted numbers re-checked against the file |
| --- | --- | --- |
| `blackadar1962` | supplied | eq. 24 `l = kz/(1 + kz/lambda)`, eq. 25 `lambda = 0.00027 G/f`, the 32-degree Brookhaven deflection at z0 = 1 m and G = 10 m/s, the "ruled out as a factor affecting characteristics of the free atmosphere" argument for G/f, the "a priori just as acceptable" concession to u*/f, and 8 pages: all hold |
| `goll2012` | `unpaywall:publisher` | row quotes none |
| `jackson_1996` | `scholar` | row quotes none; the cumulative form `Y = 1 - beta^d` and the per-biome beta fits are there, attributed within the paper to Gale and Grigal (1987) |
| `king_2005` | supplied | `Db` proportional to `H^1.5` and structural mass to `H^4` above 8 m: both hold |
| `mcgroddy2004` | supplied | CV of foliar C:P 79 per cent against C:N 59 per cent, over Table 1's n = 59 and n = 55 senesced-litter stands, on a molar basis: all hold |
| `mcmahon_1973` | supplied | critical height as `(E/(rho*g))^(1/3) * d^(2/3)` with rho*g written as weight per unit volume, and the 2/3 exponent surviving hollow, conical, paraboloid and tapered columns up to the constant: all hold |
| `papastefanou_2024` | supplied | the named missing processes -- ground water uptake, hydraulic redistribution, dynamic rooting, lag effects, sub-daily aggregation: all hold |
| `parton2010-forcent` | `unpaywall:repository` | the row's claim that it contains no phosphorus at all: holds, zero occurrences of the stem in the full text |
| `pilegaard_2013` | `europepmc` | the 60 per cent WFPS crossover, NO:N2O near 1 there, and the 15 to 65 per cent per-soil optimum range: all hold |
| `reich2004-leaf` | `europepmc` | **13.8 IS NOT IN THIS PAPER.** See below |
| `yang2011-hedley` | `unpaywall:publisher` | "labile P is generally much higher than vegetation demand, even in highly weathered soils commonly considered P limited": holds |
| `yang2013-global-soil-phosphorus` | `unpaywall:publisher` | labile, organic, occluded and secondary at 3.6 +/- 3, 8.6 +/- 6, 12.2 +/- 8 and 3.2 +/- 2 PgP against a 40.6 +/- 18 total, and the Hedley-labile definition as inorganic plus organic: all hold |
| `yuan2011-fine-root` | `unpaywall:publisher` | 631 live-root samples, p = 0.271 against Reich and Oleksyn and p = 0.120 against Wright et al., live-root N:P 16.0 against leaf 13.8 and 18.2: all hold |

**The one quoted number that did not verify.** The `reich2004-leaf` row named
Reich and Oleksyn (2004) as the source of the 13.8 end of `PFRAC_LEAFTOROOT`'s
bracket. The paper is the compilation -- 5,087 observations of leaf N and P for
1,280 species at 452 sites -- and it presents leaf N:P as fits and binned figures
against mean annual temperature and latitude, printing no grand mean anywhere.
13.8 is Yuan et al. (2011)'s summary statistic OF that compilation, stated in
their live-root comparison. The number is sound and its provenance was one
citation short; `references/INDEX.md` now sends a reader chasing it to the paper
that prints it.

## Which of the thirteen are read by code or config

Eight of the thirteen are cited from something that runs, not only from prose.
That is what makes an absent file a live exposure rather than a filing matter:
the constant is in the tree and its derivation was not.

| paper | where it is read |
| --- | --- |
| `blackadar1962` | `vendor/exoplasim/exoplasim/plasim/src/fluxmod.f90`, `config/planet.yaml`, `exoplasim/scripts/run_exoplasim.py`, and every arm file under `exoplasim/analysis/arms/` |
| `pilegaard_2013` | `vendor/lpj-guess/modules/ntransform.cpp`, `biosphere/config/ntransform.yaml` |
| `mcgroddy2004` | `vendor/lpj-guess/framework/guess.h` (`PFRAC_MINTOMAX`), `vendor/lpj-guess/modules/somdynam.cpp`, `biosphere/config/somdynam.yaml` |
| `yuan2011-fine-root` | `vendor/lpj-guess/framework/guess.h` (`PFRAC_LEAFTOROOT`) |
| `reich2004-leaf` | `vendor/lpj-guess/framework/guess.h`, as the bracket's lower end |
| `yang2013-global-soil-phosphorus` | `vendor/lpj-guess/modules/somdynam.cpp`, as the definition of the Hedley-labile pool `PMASS_SAT` is converted into |
| `parton2010-forcent` | `vendor/lpj-guess/framework/guess.h`, `modules/somdynam.cpp`, `modules/driver.cpp` -- cited there for the CENTURY C:P ratios it does not contain |
| `jackson_1996` | `vendor/lpj-guess/data/ins/*.ins`, as `rootdistribution "jackson"` and the per-PFT `root_beta` values |

The remaining five -- `goll2012`, `yang2011-hedley`, `mcmahon_1973`, `king_2005`,
`papastefanou_2024` -- are cited from prose only: `biosphere/notes/`,
`notes/audits/` and this index.

## What refused the five that were supplied by hand

Recorded because it is what a future fetch will hit, not as a narrative.

| paper | what refused |
| --- | --- |
| Blackadar (1962), `10.1029/JZ067i008p03095` | not open access; `doi.org` 403, no repository copy, Scholar found no candidate, last resort 403 |
| King (2005), `10.1016/j.ecolmodel.2004.11.017` | not open access; publisher landing page carried no reachable PDF, no repository copy, last resort 403 |
| McGroddy et al. (2004), `10.1890/03-0351` | not open access; publisher 403, the HAL record `hal-04928275` carries metadata and no file, Scholar captcha, last resort 403 |
| McMahon (1973), `10.1126/science.179.4079.1201` | not open access; publisher 403, Scholar captcha, last resort 403 |
| Papastefanou et al. (2024), `10.1088/1748-9326/ad8f48` | GOLD OPEN ACCESS and blocked by bot management at every route: the IOP article PDF URL serves a captcha to `paperfetch` and to a browser-impersonating curl; the Birmingham repository file returns 403 even with the Cloudflare cookie its own landing page sets; the VU and DOAJ records carry no file |

## The two checks

**`scripts/check_reference_index.py`** reads `references/INDEX.md`'s first column
and fails on a backticked filename that is not on disk. It reproduced all
thirteen and nothing else, including no false positive on the OCR variant quoted
inside the `meybeck1987` row's citation cell -- only the FIRST cell is a claim to
hold a file. It reports rather than fails in a checkout without the untracked
library.

**`scripts/link_worktree.py --check`** now fails on ignored payload that exists
only in the worktree. It is the more general of the two: it catches the same
failure in `exoplasim/runs/` and `source/`, where the artifact that dies is a
climate run or an export payload rather than a citation. Exercised both ways --
a clean worktree passes, and a probe PDF under `references/` plus a probe run
directory under `exoplasim/runs/` are both reported, the run directory as one
line rather than as its contents.

The first check catches the symptom in the one place it has already been paid
for. The second catches the cause everywhere it can occur.

## The other direction

This audit is about the file a worktree CREATES, which lives only there. The
same per-file linking has a second and sharper direction: an EXISTING linked
file is a symlink into the main checkout, so a worktree that regenerates one
writes through and invalidates whatever the main checkout built from those
bytes. `notes/audits/worktree-write-through.md` audits it.
