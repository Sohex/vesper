#!/usr/bin/env python3
"""Apply or revert the WORLD-43RK arm: lwr's absorptances without the masks.

    python exoplasim/scripts/patch_unmasked_absorptance.py --apply
    python exoplasim/scripts/patch_unmasked_absorptance.py --check
    python exoplasim/scripts/patch_unmasked_absorptance.py --revert

Worldbuilding frame: this switches the Vesper climate model's longwave source
between two forms of the same four absorptance fits, for a COMPUTE measurement.
Nothing here is about the simulated planet.

WHY THIS IS A SCRIPT AND NOT A COMMIT. The unmasked form is measured and NOT
adopted: `exoplasim/notes/masked-radiation-and-four-bytes.md` carries the
numbers, and CLIM-84's sequencing verdict stands -- the change is not bit
identical, and it must not be taken before CLIM-61 settles whether the
broadband scheme survives at all. So the tree holds the masked form, this file
holds the arm, and either can be reproduced exactly.

IT IS A TEXT SUBSTITUTION, NOT A PATCH FILE. `exoplasim/patches/` is the
authored record of changes that WERE taken and is applied to nothing; an arm
that was not taken does not belong there. The two forms are carried here
verbatim, so --check can say which one the tree holds and neither can drift.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from _paths import MODEL_SRC  # noqa: E402

RADMOD = MODEL_SRC / "plasim" / "src" / "radmod.f90"

DECL_MASKED = """\
      real zth2o(NHOR)          ! water vapor - co2 overlap transmissivity
      real zbdl(NHOR)           ! layer evective downward rad."""

DECL_UNMASKED = """\
      real zth2o(NHOR)          ! water vapor - co2 overlap transmissivity
      real zsel(NHOR)           ! branch selector, exactly 1. or exactly 0.
      real zbdl(NHOR)           ! layer evective downward rad."""

BODY_MASKED = """\
        where(zsumwv(:) <= 0.01)
         zah2o(:)=0.846*max(0.,zsumwv(:)+3.59E-5)**0.243-zh2o0a
        elsewhere
         zah2o(:)=0.24*ALOG10(max(1.E-30,zsumwv(:)+0.01))+zah2oc
        endwhere
!
!     b) continuum
!
        if(th2oc > 0.) then
         zah2o(:)=AMIN1(zah2o(:)+(1.-exp(-th2oc*zsumwv(:))),1.)
        endif
!
!     co2 absorption:
!
        where(zsumco2(:) <= 1.0)
         zaco2(:)=0.0676*max(0.,zsumco2(:)+0.01022)**0.421-zco20
        elsewhere
         zaco2(:)=0.0546*ALOG10(max(1.E-30,zsumco2(:)))+zaco2c
        endwhere
!
!     Boer et al. (1984) scheme for t(h2o) at co2 overlapp
!
        where(zsumwv(:)<= 2.)
         zth2o(:)=1.-(0.832*max(0.,zsumwv(:)+0.0286)**0.26-zh2o0)
        elsewhere
         zth2o(:)=max(0.,zth2oc-0.1196*log(max(1.E-30,zsumwv(:)-0.6931)))
        endwhere
!
!     o3 absorption:
!
        where(zsumo3(:) <= 0.01)
         zao3(:)= 0.209*max(0.,zsumo3(:)+7.E-5)**0.436 - zao30
        elsewhere
         zao3(:)= 0.0212*log10(max(1.E-30,zsumo3(:)))+zao3c
        endwhere
!
"""

BODY_UNMASKED = """\
!     WORLD-43RK, THE UNMASKED FORM. Each of these four absorptances selects
!     between a fractional-power fit below a path amount and a logarithmic fit
!     above it. Written as `where`/`elsewhere` the assignment is masked, and GCC
!     does not vectorise a masked array assignment, so every one of these
!     compiled to a SCALAR `pow` and libmvec's eight-wide `pow` was never
!     reached. Written as one unmasked array assignment per absorptance, with
!     the branch chosen arithmetically, the same expressions vectorise.
!
!     WHY EVALUATING BOTH BRANCHES ON EVERY LANE IS SAFE. A `where` selects the
!     assignment and not the evaluation, so an unmasked form is only admissible
!     if every lane is in domain -- `config/planet.yaml` compiles with
!     -ffpe-trap=invalid,zero,overflow and an excluded lane out of domain raises
!     at the instruction that evaluated it. Every argument here already carries
!     an unconditional floor INSIDE it: max(0.,.) before each fractional power,
!     so the base is never negative and pow(0.,a)=0. for a>0.; and
!     max(1.E-30,.) before each logarithm, so the argument is never zero and
!     never negative. Both branches are therefore finite on every lane whatever
!     the path amount is, and neither can raise the invalid, the divide-by-zero
!     or the overflow. That was true before this change; it is what makes the
!     change available. notes/audits/masked-where-blocks.md.
!
!     WHY THE BLEND IS EXACT. zsel is 0.5+SIGN(0.5,.), which is exactly 1. or
!     exactly 0. and nothing between, so one term is the branch value unchanged
!     and the other is 0. times a FINITE number, which is a signed zero and adds
!     nothing. The restructuring moves no result. What does move the result is
!     libmvec: its vector routines carry a looser ULP bound than the scalar ones
!     it replaces, so this is not bit identical to the masked form and the
!     difference lives entirely in the library, not in the algebra.
!
!     h2o 6.3mu: below 0.01 the Sasamori power fit, above it the log fit.
!
        zsel(:)=0.5+SIGN(0.5,0.01-zsumwv(:))
        zah2o(:)=zsel(:)*(0.846*max(0.,zsumwv(:)+3.59E-5)**0.243-zh2o0a)      &
     &          +(1.-zsel(:))*(0.24*ALOG10(max(1.E-30,zsumwv(:)+0.01))+zah2oc)
!
!     b) continuum
!
        if(th2oc > 0.) then
         zah2o(:)=AMIN1(zah2o(:)+(1.-exp(-th2oc*zsumwv(:))),1.)
        endif
!
!     co2 absorption:
!
        zsel(:)=0.5+SIGN(0.5,1.0-zsumco2(:))
        zaco2(:)=zsel(:)*(0.0676*max(0.,zsumco2(:)+0.01022)**0.421-zco20)     &
     &          +(1.-zsel(:))*(0.0546*ALOG10(max(1.E-30,zsumco2(:)))+zaco2c)
!
!     Boer et al. (1984) scheme for t(h2o) at co2 overlapp
!
        zsel(:)=0.5+SIGN(0.5,2.-zsumwv(:))
        zth2o(:)=zsel(:)*(1.-(0.832*max(0.,zsumwv(:)+0.0286)**0.26-zh2o0))    &
     &          +(1.-zsel(:))                                                 &
     &           *max(0.,zth2oc-0.1196*log(max(1.E-30,zsumwv(:)-0.6931)))
!
!     o3 absorption:
!
        zsel(:)=0.5+SIGN(0.5,0.01-zsumo3(:))
        zao3(:)=zsel(:)*(0.209*max(0.,zsumo3(:)+7.E-5)**0.436-zao30)          &
     &         +(1.-zsel(:))*(0.0212*log10(max(1.E-30,zsumo3(:)))+zao3c)
!
"""


PAIRS = ((DECL_MASKED, DECL_UNMASKED), (BODY_MASKED, BODY_UNMASKED))


def state(text: str) -> str:
    """Which form the file holds, or an error naming what it holds instead."""
    masked = all(text.count(a) == 1 for a, _ in PAIRS)
    unmasked = all(text.count(b) == 1 for _, b in PAIRS)
    if masked and not unmasked:
        return "masked"
    if unmasked and not masked:
        return "unmasked"
    return "unrecognised"


def main() -> None:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--apply", action="store_true",
                   help="replace the four masked absorptance selections with "
                        "the unmasked form")
    g.add_argument("--revert", action="store_true",
                   help="put the masked form back")
    g.add_argument("--check", action="store_true",
                   help="say which form radmod.f90 holds; exit 0 either way")
    args = ap.parse_args()

    text = RADMOD.read_text(encoding="utf-8")
    now = state(text)
    if now == "unrecognised":
        raise SystemExit(
            f"{RADMOD} holds neither form of lwr's absorptance block. Either it "
            f"has been edited since this arm was written, or the arm has already "
            f"been adopted. Read the file rather than forcing either form onto "
            f"it; exoplasim/notes/masked-radiation-and-four-bytes.md says what "
            f"the two forms are.")
    if args.check:
        print(f"{RADMOD.name}: {now}")
        return
    want = "unmasked" if args.apply else "masked"
    if now == want:
        print(f"{RADMOD.name} is already {want}; nothing to do")
        return
    for a, b in PAIRS:
        text = text.replace(a, b) if args.apply else text.replace(b, a)
    RADMOD.write_text(text, encoding="utf-8")
    print(f"{RADMOD.name}: {now} -> {want}")
    if args.apply:
        print("The tree now holds a measurement arm that is NOT adopted. Revert "
              "it before committing anything under vendor/exoplasim.")


if __name__ == "__main__":
    sys.exit(main())
