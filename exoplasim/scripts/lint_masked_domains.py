#!/usr/bin/env python3
"""Domain-sensitive intrinsics and divisions that sit inside a masked WHERE block.

    python exoplasim/scripts/lint_masked_domains.py            # the survivors
    python exoplasim/scripts/lint_masked_domains.py --all      # every site, guarded or not
    python exoplasim/scripts/lint_masked_domains.py --tsv      # machine-readable
    python exoplasim/scripts/lint_masked_domains.py --kind divide   # one class alone

The default is `--kind any` and that is the gate: the intrinsics, the divisions
and the merges together. It was `--kind intrinsic` while radmod.f90's 45 divisor
keys were unread, and world-px61 read them.

Worldbuilding frame: this reads the Vesper climate model's Fortran source.
Nothing here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. A Fortran `where` block does NOT protect its right-hand sides.
The mask selects which elements the ASSIGNMENT stores; the compiler is free to
evaluate the expression on every element of the array, and a vectorising one
does exactly that. The build profile carries
`-ffpe-trap=invalid,zero,overflow` (`config/planet.yaml`), so a lane that was
going to be discarded raises SIGFPE instead of being discarded. That is not
hypothetical: a T42 run died eight orbits in and the fault moved to a new site
each time one was clamped (world-bhs, world-5a0).

`-O2` rather than `-O3` removed the loop transforms that were firing. It did
not remove the class, because `-O2` vectorises too. So the class needs a source
answer, and a source answer needs the population enumerated.

THE PASS OVER-REPORTS, DELIBERATELY, and that direction is the whole of its
evidential value. Three places it errs towards reporting a site that is safe,
and none towards clearing one that is not:

  1. The WHERE stack is popped only by an explicit `end where` / `endwhere`,
     or by a procedure boundary. A construct this parser fails to close keeps
     reporting sites after it, never fewer.
  2. A guard is recognised only from a small, literal set of shapes (an
     argument already wrapped in `max`/`min`/`abs`, or a divisor already
     wrapped in `max`/`sign`). An unrecognised but correct guard is reported.
  3. Division is reported as its own class even though a divisor is only a
     hazard when it can reach zero, which the pass cannot know. All three
     spellings are reported -- a parenthesised expression, a bare name, and a
     name with a subscript, section or argument list -- because the trap does
     not care how the divisor was written.

So a clean run is evidence and a dirty one is a work list. A pass whose
direction is unstated is not evidence at all, which is why this paragraph is
here.

WHAT IT DOES NOT KNOW. It cannot tell whether a mask's complement actually
holds an out-of-domain value -- that needs the physics, and it is why the
output is read in context rather than acted on mechanically. It also has no
opinion on `merge(a,b,mask)`, which has the identical hazard for the same
reason, so those sites are reported under their own class.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import MODEL_SRC  # noqa: E402

SRC = MODEL_SRC / "plasim" / "src"

# The intrinsics whose argument has a restricted domain, in every spelling the
# source uses. Longest first so `log10` is not matched as `log`.
DOMAIN_INTRINSICS = (
    "alog10",
    "log10",
    "alog",
    "log",
    "sqrt",
    "exp",
    "asin",
    "acos",
    "acosh",
    "atanh",
)
_INTRINSIC_RE = re.compile(
    r"(?<![a-z0-9_])(" + "|".join(DOMAIN_INTRINSICS) + r")\s*\(", re.I
)
_MERGE_RE = re.compile(r"(?<![a-z0-9_])merge\s*\(", re.I)

# Block openers/closers. `where (mask)` with nothing after the closing paren
# opens a construct; with a statement after it, the mask covers that one
# assignment and the construct is not open.
_WHERE_OPEN_RE = re.compile(r"^\s*(?:[a-z_]\w*\s*:\s*)?where\s*\(", re.I)
_END_WHERE_RE = re.compile(r"^\s*end\s*where\b", re.I)
_ELSEWHERE_RE = re.compile(r"^\s*else\s*where\b", re.I)
_PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:(?:real|integer|logical|complex|double\s+precision|character)"
    r"(?:\s*\([^)]*\))?\s+)?"
    r"(subroutine|function|program|module)\s+[a-z_]\w*",
    re.I,
)
_PROC_END_RE = re.compile(r"^\s*end\s*(subroutine|function|program|module)\b", re.I)

# Guard shapes already in the tree. An argument wrapped in one of these has had
# its domain bounded by hand; a divisor wrapped in max/sign has had its pole
# floored. Recognised literally, so an unrecognised guard is REPORTED.
_ARG_GUARDS = ("max", "min", "amax1", "amin1", "abs", "dim")
_DIV_GUARDS = ("max", "amax1", "sign", "abs")


# ---------------------------------------------------------------------------
# THE CLASSIFIED SURVIVORS, and why each is not a defect.
#
# world-5a0 read every site the pass reported and made the arguments safe where
# a floor could be argued free. The sites below are the ones where NO guard was
# added, each with the argument that makes it safe. They are keyed on the file,
# the procedure, the intrinsic and the ARGUMENT TEXT rather than on a line
# number, so the table survives edits above them and a changed argument drops
# out of the table and is reported again -- which is the direction that matters.
#
# `--all` prints them; the default run prints only what is NOT here, so a new
# masked intrinsic is visible the moment it is added. Adding a row is a claim
# with a reason, and a row whose reason does not hold is a defect in this table.
#
# Evidence: notes/audits/masked-where-blocks.md.
# ---------------------------------------------------------------------------

CLASSIFIED: dict[tuple[str, str, str, str], str] = {
    ("radmod.f90", "subroutine swr", "alog", "3.+0.1*ztau2(:)"):
        "ztau2 is preset to 1 and is otherwise a product of non-negative "
        "factors, so the argument is at least 3",
    ("radmod.f90", "subroutine swr", "sqrt", "273./dt(:,jlev)"):
        "the mask is losun and the argument is the gridpoint temperature, "
        "which the mask has no bearing on; masked and kept lanes are equally "
        "safe and this is not a masked-domain site",
    ("radmod.f90", "subroutine swr", "exp", "zaertf1(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf1(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf1s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf1s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf2(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf2(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "zaertf2s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "exp", "-zaertf2s(:,jlev)"):
        "preset to 0 and, where evaluated, MIN(25.,nonneg)",
    ("radmod.f90", "subroutine swr", "log", "1.+zcb1*zco2(:)"):
        "the column amount is preset to 0 unconditionally and accumulates "
        "non-negative terms, so the argument is at least 1",
    ("radmod.f90", "subroutine swr", "log", "1.+zcb2*zco2(:)"):
        "the column amount is preset to 0 unconditionally and accumulates "
        "non-negative terms, so the argument is at least 1",
    ("radmod.f90", "subroutine swr", "log", "1.0-0.144"):
        "constant argument",
    ("radmod.f90", "subroutine swr", "exp", "zscf(:)*log(1.0-0.144)"):
        "constant logarithm times a non-negative scale factor",
    ("radmod.f90", "subroutine swr", "log",
     "1.0-(0.219/(1.+0.816*max(0.,zmu0(:))))"):
        "the floor holds the argument at or above 0.781 for any zmu0",
    ("radmod.f90", "subroutine swr", "exp",
     "zscf(:)*log(1.0-(0.219/(1.+0.816*max(0.,zmu0(:)))))"):
        "the logarithm inside it is bounded, and zscf is non-negative",
    ("radmod.f90", "subroutine lwr", "alog", "ztau0(:)"):
        "clamped to [zero, 1-zero] on every lane by the AMIN1/MAX two lines "
        "above the where, unconditionally",
    ("seamod.f90", "subroutine seaini", "exp",
     "ra2*(dt(:,NLEP)-TMELT) /ra4d(dt(:,NLEP),ra4)"):
        "ra4d floors the pole; the quotient tends to ra2 for large dt and to "
        "a large negative for small dt, which underflows rather than traps",
    ("seamod.f90", "subroutine seastep", "exp",
     "ra2*(dt(:,NLEP)-TMELT) /ra4d(dt(:,NLEP),ra4)"):
        "ra4d floors the pole; the quotient tends to ra2 for large dt and to "
        "a large negative for small dt, which underflows rather than traps",
    ("seamod.f90", "subroutine seastep", "sqrt", "dtaux(:)**2+dtauy(:)**2"):
        "a sum of squares cannot be negative on any lane",

    # -----------------------------------------------------------------------
    # THE CLASSIFIED DIVISIONS.
    #
    # world-d016 re-derived the population after the pass was extended to
    # bare-name divisors -- 133 parenthesised sites became 345 across three
    # spellings -- and read every one outside radmod.f90. Each row below is a
    # DIVISOR that cannot vanish, with the argument that says why. They fall
    # into five shapes:
    #
    #   * a scalar or a named constant, where the mask has no lane bearing;
    #   * a field positive on every gridpoint, sea and land alike;
    #   * an expression bounded away from zero by the physics of its terms;
    #   * a divisor already floored, here or under world-5a0;
    #   * a local now PRESET on every lane, so the masked evaluation has a
    #     defined divisor rather than whatever the stack held.
    #
    # A row saying "the argument is in the comment at the site" is one where a
    # floor or a preset was added, and the argument for its being a no-op on
    # the lanes the mask keeps belongs beside the code it changed.
    #
    # radmod.f90's 45 keys were the last of them, and they are below. They were
    # left unread while that file was edited on another branch, because a row
    # keyed on argument text from a version that no longer exists is a claim
    # about nothing; world-px61 read them once it settled. With them in, every
    # class this pass reports is classified, so `--kind any` is the default and
    # the gate covers the divisions and the merges as well as the intrinsics.
    # -----------------------------------------------------------------------

    ("fluxmod.f90", "subroutine mkevap", "divide", "deltsec"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the timestep in seconds, set from the calendar in plasim.f90 and "
        "positive for every run",
    ("fluxmod.f90", "subroutine mkevap", "divide", "dp(:)"):
        "the surface pressure, positive on every gridpoint whatever the mask "
        "selects",
    ("fluxmod.f90", "subroutine mkevap", "divide", "ra4d(dt(:,NLEP),ra4)"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("fluxmod.f90", "subroutine mkevap", "divide", "ra4d(dt(:,NLEP),ra4i)"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("fluxmod.f90", "subroutine mkevap", "divide", "zkonst2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "dsigma(NLEV)/deltsec2/ga, a product and quotient of positive "
        "constants",
    ("icemod.f90", "subroutine iceout", "divide", "crhosn"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the snow density, an icemod scalar at 330 kg/m3",
    ("icemod.f90", "subroutine icestep", "divide", "CRHOI"):
        "parameter(CRHOI = 920.) in icemod's header, so the divisor is "
        "nonzero by declaration",
    ("icemod.f90", "subroutine icestep", "divide", "real(naccuo)"):
        "the accumulation counter, incremented unconditionally on the "
        "statement above the block that divides by it, so it is at least 1",
    ("icemod.f90", "subroutine icestep", "divide", "xdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the coupled timestep, solar_day divided by the timesteps per day",
    ("icemod.f90", "subroutine icestep", "divide", "zrhoilfdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "CRHOI*CLFI/xdt, a product and quotient of positive constants",
    ("icemod.f90", "subroutine make_ice_thickness", "divide", "cmaxn-cminn"):
        "parameter(cminn=0.1) and parameter(cmaxn=0.9) in the same routine, "
        "so the divisor is 0.8 by declaration",
    ("icemod.f90", "subroutine mkcflux", "divide", "CKAPI"):
        "parameter(CKAPI = 2.03) in icemod's header, so the divisor is "
        "nonzero by declaration",
    ("icemod.f90", "subroutine mkcflux", "divide", "CKAPSN"):
        "parameter(CKAPSN = 0.31) in icemod's header, so the divisor is "
        "nonzero by declaration",
    ("icemod.f90", "subroutine mkflukoi", "divide", "taunc"):
        "the enclosing `if (taunc > 0.)` is a SCALAR branch and not a mask, "
        "so the division is not reached when taunc is zero; the `where` "
        "inside it selects lanes and not whether the branch runs",
    ("icemod.f90", "subroutine mkflukoi", "divide", "xdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the coupled timestep, solar_day divided by the timesteps per day",
    ("icemod.f90", "subroutine mkice", "divide", "zrhoilfdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "CRHOI*CLFI/xdt, a product and quotient of positive constants",
    ("icemod.f90", "subroutine mkicec", "divide", "2.*max(picedo(:),1.0e-30)"):
        "floored under world-d016; the argument is in the comment at the site",
    ("icemod.f90", "subroutine mkicec", "divide", "hlead"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the lead-closing growth scale, an icemod scalar at 0.5 m",
    ("icemod.f90", "subroutine subsnow", "divide", "CRHOS"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the sea water density, an icemod scalar at 1030 kg/m3",
    ("icemod.f90", "subroutine subsnow", "divide", "xdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the coupled timestep, solar_day divided by the timesteps per day",
    ("landmod.f90", "subroutine landstep", "divide", "tmelt-263.16"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "tmelt is the melting point and the difference with the ramp's lower "
        "end is a constant",
    ("landmod.f90", "subroutine mkradv", "divide", "zarea(1,:)"):
        "a gridcell area, positive on every cell of the Gaussian grid",
    ("landmod.f90", "subroutine mkradv", "divide", "zarea(:,:)"):
        "a gridcell area, positive on every cell of the Gaussian grid",
    ("landmod.f90", "subroutine mkradv", "divide", "zarea(NLON,:)"):
        "a gridcell area, positive on every cell of the Gaussian grid",
    ("landmod.f90", "subroutine mkradv", "divide", "zarea(jlon+1,:)"):
        "a gridcell area, positive on every cell of the Gaussian grid",
    ("landmod.f90", "subroutine mkradv", "divide", "zarea(jlon,:)"):
        "a gridcell area, positive on every cell of the Gaussian grid",
    ("landmod.f90", "subroutine tands", "divide", "(ALS-ALV)*1000."):
        "the latent heat of fusion, the difference of two plasimmod "
        "constants, times a positive literal",
    ("landmod.f90", "subroutine tands", "divide", "deltsec"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the timestep in seconds, set from the calendar in plasim.f90 and "
        "positive for every run",
    ("landmod.f90", "subroutine tands", "divide", "rhosnow"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the snow density, a landmod scalar",
    ("landmod.f90", "subroutine tands", "divide",
     "snowdiff*zsoilz(:,1)+zdiff(:,1)*zsnowz(:)"):
        "the world-d016 presets in tands give every lane, sea included, "
        "zsoilz(:,1) = dsoilz(1) = 0.4 m, zdiff1 = soildifsat, zctop = soilcapsat "
        "and zsnowz = 0, and each is overwritten with a positive value on the "
        "lanes the mask keeps; the divisor is therefore at least "
        "snowdiff*dsoilz(1)",
    ("landmod.f90", "subroutine tands", "divide", "zctop(:)"):
        "the world-d016 presets in tands give every lane, sea included, "
        "zsoilz(:,1) = dsoilz(1) = 0.4 m, zdiff1 = soildifsat, zctop = soilcapsat "
        "and zsnowz = 0, and each is overwritten with a positive value on the "
        "lanes the mask keeps",
    ("landmod.f90", "subroutine tands", "divide",
     "zctop(:)*zztop(:)/deltsec+2.*zdiff1(:)/zsoilz1(:)"):
        "the world-d016 presets in tands give every lane, sea included, "
        "zsoilz(:,1) = dsoilz(1) = 0.4 m, zdiff1 = soildifsat, zctop = soilcapsat "
        "and zsnowz = 0, and each is overwritten with a positive value on the "
        "lanes the mask keeps; every term is then positive",
    ("landmod.f90", "subroutine tands", "divide", "zsoilz1(:)"):
        "the world-d016 presets in tands give every lane, sea included, "
        "zsoilz(:,1) = dsoilz(1) = 0.4 m, zdiff1 = soildifsat, zctop = soilcapsat "
        "and zsnowz = 0, and each is overwritten with a positive value on the "
        "lanes the mask keeps",
    ("landmod.f90", "subroutine tands", "divide", "zztop(:)"):
        "dztop, assigned to every lane in the preset block above the mask",
    ("oceanmod.f90", "subroutine addfc", "divide", "ymld(:,1)"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth",
    ("oceanmod.f90", "subroutine addfc", "divide", "zcpsdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "dtmix/(CRHOS*CPS), a quotient of positive constants",
    ("oceanmod.f90", "subroutine hdiffo", "divide", "dlam"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "2*pi/NLON, the ocean grid's longitude spacing",
    ("oceanmod.f90", "subroutine hdiffo", "divide", "dphi(jlat)"):
        "the latitude spacing of the Gaussian grid, positive for every jlat",
    ("oceanmod.f90", "subroutine hdiffo", "divide", "dtmix"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "solar_day divided by the timesteps per day, positive for every run",
    ("oceanmod.f90", "subroutine mkfc", "divide", "dtmix"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "solar_day divided by the timesteps per day, positive for every run",
    ("oceanmod.f90", "subroutine mkfc", "divide", "taunc"):
        "the enclosing `if (taunc > 0.)` is a SCALAR branch and not a mask, "
        "so the division is not reached when taunc is zero; the `where` "
        "inside it selects lanes and not whether the branch runs",
    ("oceanmod.f90", "subroutine mkiflux", "divide", "zcpsdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "dtmix/(CRHOS*CPS), a quotient of positive constants",
    ("oceanmod.f90", "subroutine mksst", "divide", "ymld(:,1)"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth",
    ("oceanmod.f90", "subroutine mksst", "divide", "zcpsdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "dtmix/(CRHOS*CPS), a quotient of positive constants",
    ("oceanmod.f90", "subroutine vdiffo", "divide", "ymld(:,1)+zk(:,1)"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth, and zk is a non-negative diffusion coefficient",
    ("oceanmod.f90", "subroutine vdiffo", "divide",
     "ymld(:,NLEV_OCE)+zk(:,nlem_oce) *(1.-zebs(:,nlem_oce))"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth, zk is non-negative, and zebs is below 1 for the same reason",
    ("oceanmod.f90", "subroutine vdiffo", "divide",
     "ymld(:,jlev)+zk(:,jlev) +zk(:,jlem)*(1.-zebs(:,jlem))"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth, zk is non-negative, and zebs is zk/(ymld+zk) from the level "
        "above, which is below 1",
    ("oceanmod.f90", "subroutine vdiffo", "divide",
     "ymld(:,jlev+1)+ymld(:,jlev)"):
        "ymld is assigned dlayer(jlev) for EVERY gridpoint in the ocean "
        "initialisation, land and sea alike, and dlayer is a positive layer "
        "depth; a sum of two of them cannot vanish",
    ("plasim.f90", "subroutine initpm", "divide",
     "TWOPI * max(restim,1.0e-30)"):
        "floored under world-d016; the argument is in the comment at the site",
    ("plasim.f90", "subroutine initpm", "divide", "TWOPI * max(tfrc,1.0e-30)"):
        "floored under world-d016; the argument is in the comment at the site",
    ("rainmod.f90", "subroutine kuo", "divide", "1.-rhbeta"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "rhbeta is the namelist relative-humidity threshold of the Kuo beta "
        "closure",
    ("rainmod.f90", "subroutine kuo", "divide", "deltsec2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "twice the timestep in seconds, positive for every run",
    ("rainmod.f90", "subroutine kuo", "divide", "ga"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the surface gravity, which config/planet.yaml declares positive",
    ("rainmod.f90", "subroutine kuo", "divide", "zalcpe(:,jlev)"):
        "ALV/(acpd*(1.+ADV*dq)), assigned for every lane at the head of kuo; "
        "ALV and acpd are positive and 1+ADV*q cannot vanish for a mass "
        "fraction",
    ("rainmod.f90", "subroutine mkclouds", "divide", "real(max(1,icctot(:)))"):
        "floored at 1 under world-5a0",
    ("rainmod.f90", "subroutine mkclouds", "divide", "sigma(jlev)*dp(:)"):
        "a sigma level times the surface pressure, both positive on every "
        "lane",
    ("rainmod.f90", "subroutine mkdca", "divide", "zsum1(:)"):
        "preset to 1 at the head of the iteration under world-d016; the "
        "argument is in the comment at the site",
    ("rainmod.f90", "subroutine mkdca", "divide", "zsumq1(:)"):
        "preset to 1 at the head of the iteration under world-d016; the "
        "argument is in the comment at the site",
    ("rainmod.f90", "subroutine mklsp", "divide", "1.-(1./rdbrv-1.)*zqsat(:)"):
        "1/rdbrv - 1 is about 0.608, so the divisor vanishes only at a "
        "specific humidity of about 1.645 kg/kg, which is above the 1 kg/kg a "
        "mass fraction cannot exceed",
    ("rainmod.f90", "subroutine mklsp", "divide",
     "1.0+zlcpe(:)*ra2s(zt(:))*(TMELT-ra4s(zt(:))) *zqsat(:)*zcor(:)/ra4d(zt(:),ra4s(zt(:)))**2"):
        "one plus a product of non-negative factors: zlcpe is a latent heat "
        "over acpd*(1+ADV*q), ra2s is positive, TMELT minus either saturation "
        "pole is positive, zqsat and zcor are non-negative and ra4d is at "
        "least 1, so the divisor is at least 1",
    ("rainmod.f90", "subroutine mklsp", "divide",
     "1.0+zlcpe(:)*ra2s(ztn(:))*(TMELT-ra4s(ztn(:))) *zqsat(:)*zcor(:)/ra4d(ztn(:),ra4s(ztn(:)))"):
        "one plus a product of non-negative factors: zlcpe is a latent heat "
        "over acpd*(1+ADV*q), ra2s is positive, TMELT minus either saturation "
        "pole is positive, zqsat and zcor are non-negative and ra4d is at "
        "least 1, so the divisor is at least 1",
    ("rainmod.f90", "subroutine mklsp", "divide", "deltsec2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "twice the timestep in seconds, positive for every run",
    ("rainmod.f90", "subroutine mklsp", "divide", "dp(:)*sigma(jlev)"):
        "the surface pressure times a sigma level, both positive on every "
        "lane",
    ("rainmod.f90", "subroutine mklsp", "divide", "ra4d(zt(:),ra4s(zt(:)))"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("rainmod.f90", "subroutine mklsp", "divide", "ra4d(ztn(:),ra4s(ztn(:)))"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("rainmod.f90", "subroutine mklsp", "divide", "rdbrv"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the ratio of the dry and vapour gas constants, about 0.622",
    ("rainmod.f90", "subroutine mkrain", "divide",
     "1.0+zlcpe(:)*ra2s(zt(:))*(TMELT-ra4s(zt(:))) *zqsat(:)*zcor(:)/ra4d(zt(:),ra4s(zt(:)))**2"):
        "one plus a product of non-negative factors: zlcpe is a latent heat "
        "over acpd*(1+ADV*q), ra2s is positive, TMELT minus either saturation "
        "pole is positive, zqsat and zcor are non-negative and ra4d is at "
        "least 1, so the divisor is at least 1",
    ("rainmod.f90", "subroutine mkrain", "divide", "acpd*(1.+ADV*dq(:,jlev))"):
        "the dry specific heat times 1+ADV*q; ADV is 0.608 and q is a mass "
        "fraction, so the second factor is between 1 and 1.7",
    ("rainmod.f90", "subroutine mkrain", "divide", "deltsec2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "twice the timestep in seconds, positive for every run",
    ("rainmod.f90", "subroutine mkrain", "divide", "dsigma(jlev)"):
        "a sigma layer thickness, positive by construction of the vertical "
        "grid",
    ("rainmod.f90", "subroutine mkrain", "divide", "ra4d(zt(:),ra4s(zt(:)))"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("rainmod.f90", "subroutine mkshallow", "divide", "1.+adv*dq(:,jlev)"):
        "ADV is 0.608 and q is a mass fraction, so this is between 1 and 1.7 "
        "on every lane",
    ("rainmod.f90", "subroutine mkshallow", "divide", "acpd"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the dry-air specific heat, a positive plasimmod constant",
    ("rainmod.f90", "subroutine mkshallow", "divide", "deltsec2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "twice the timestep in seconds, positive for every run",
    ("rainmod.f90", "subroutine mkshallow", "divide", "dsigma(1)+zkdiff(:,1)"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide",
     "dsigma(1)+zkdiff(:,1)/zskap(1)"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide",
     "dsigma(NLEV)+zkdiff(:,NLEM)*(1.-zebs(:,NLEM))"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide",
     "dsigma(NLEV)+zkdiff(:,NLEM)/zskap(NLEV) *(1.-zebs(:,NLEM)/zskap(NLEM))"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide",
     "dsigma(jlev)+(zkdiff(:,jlev) +zkdiff(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem))) /zskap(jlev)"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide",
     "dsigma(jlev)+zkdiff(:,jlev) +zkdiff(:,jlem)*(1.-zebs(:,jlem))"):
        "zkdiff is preset to zero for every lane and written only where "
        "kshallow > 0, so on a lane the mask discards the divisor is exactly "
        "the dsigma term, a positive layer thickness; on a lane it keeps, "
        "zkdiff is non-negative and the zebs factor is below 1, so the "
        "divisor is at least that term there too",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(1)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(NLEM)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(NLEV)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(jlem)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(jlep)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("rainmod.f90", "subroutine mkshallow", "divide", "zskap(jlev)"):
        "sigma(jlev)**akap for a positive sigma, computed for every level at "
        "the head of the routine and independent of the mask",
    ("seamod.f90", "subroutine seaini", "divide", "1.-(1./rdbrv-1.)*dqs(:)"):
        "1/rdbrv - 1 is about 0.608, so the divisor vanishes only at a "
        "specific humidity of about 1.645 kg/kg, which is above the 1 kg/kg a "
        "mass fraction cannot exceed",
    # -----------------------------------------------------------------------
    # radmod.f90, 45 keys over 103 sites, world-px61.
    #
    # THE TWO REMAINDERS WERE ONE READING. Several of these divisors ARE the
    # locals lint_masked_locals.py reports in the same file: swr computed its
    # whole two-stream into automatics written only under where(losun(:)), so
    # on a night lane the divisor was whatever the stack held. All 53 are now
    # PRESET on every lane at the transparent atmosphere over a black surface
    # -- zero absorber amount, unit transmissivity, zero reflectivity, zero
    # flux -- so every quotient below divides by exactly 1 on a lane the mask
    # discards, including 1.-r*r, which becomes 1.-0.*0. Each row therefore
    # says what holds on a lane losun KEEPS, which is the only lane whose value
    # is stored.
    #
    # THE BEAM COSINE IS NON-NEGATIVE, and ten rows rest on it: solang zeroes
    # gmu0 and overwrites it only where the cosine exceeds sin(dawn), and
    # radini now refuses a negative dawn, so zmu0 is at or above zero on every
    # lane whatever the mask selects.
    #
    # NINE ROWS CITE world-2223 RATHER THAN A BOUND. The clear-sky
    # transmissivities and the adding method's reflectivities are not confined
    # by anything swr does, and the mask never confined them either -- the
    # expression is the same on a kept lane as on a discarded one, so it is not
    # a masked-domain site. Whether either can reach zero is that issue.
    # -----------------------------------------------------------------------

    ("radmod.f90", "subroutine lwr", "divide", "1.-ztau0(:)"):
        "ztau0 is clamped to [zero, 1-zero] on every lane by the AMIN1/MAX one "
        "line above the where, unconditionally, so the divisor is at least "
        "1.E-6",
    ("radmod.f90", "subroutine lwr", "divide", "ALOG(ztau0(:))"):
        "the same unconditional clamp holds ztau0 at or below 1-1.E-6, so the "
        "logarithm is at most about -1.E-6 and cannot be zero; it is in domain "
        "for the same reason, which is the alog row above",
    ("radmod.f90", "subroutine lwr", "divide", "sigma(NLEV)-sigma(NLEM)"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the gap between the two lowest full-level sigmas, which the vertical "
        "coordinate makes strictly positive",

    ("radmod.f90", "subroutine swr", "divide", "(1.+138.6*zo3(:))**0.805"):
        "zo3 is a column ozone amount: preset to 0 on every lane and otherwise "
        "a sum of non-negative layer amounts, so the base is at least 1 and so "
        "is the power",
    ("radmod.f90", "subroutine swr", "divide",
     "(1.+141.5*zwv(:))**0.635+5.925*zwv(:)"):
        "zwv is a column water vapour amount, preset to 0 and otherwise a sum "
        "of non-negative layer amounts, so the first term is at least 1 and "
        "the second is non-negative",
    ("radmod.f90", "subroutine swr", "divide", "1.+(103.6*zo3(:))**3"):
        "zo3 is non-negative for the reason above, so its cube is too and the "
        "divisor is at least 1",
    ("radmod.f90", "subroutine swr", "divide",
     "1.+0.042*zo3(:)+0.000323*zo3(:)**2"):
        "both terms are non-negative for a non-negative zo3, so the divisor is "
        "at least 1",
    ("radmod.f90", "subroutine swr", "divide", "1.+0.816*max(0.,zmu0(:))"):
        "the floor holds the second term at or above 0, so the divisor is at "
        "least 1",
    ("radmod.f90", "subroutine swr", "divide", "1.+bb*zmu0(:)"):
        "bb is a positive parameter and zmu0 is the beam cosine, non-negative "
        "on every lane, so the divisor is at least 1",
    ("radmod.f90", "subroutine swr", "divide",
     "1.+zb1(:)*ztau1(:)/max(1.E-30,zmu0(:))"):
        "zb1 is tswr1*SQRT(max(0.,zmu0)) and non-negative, ztau1 is preset to "
        "0 and otherwise a product of non-negative factors, and the inner "
        "divisor is floored, so the divisor is at least 1",
    ("radmod.f90", "subroutine swr", "divide", "1.+zb3*ztau1(:)"):
        "zb3 is tswr1*SQRT(zmu00)/zmu00, a non-negative scalar set above the "
        "level loop, and ztau1 is non-negative, so the divisor is at least 1",
    ("radmod.f90", "subroutine swr", "divide",
     "1.-(o3visw*0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2) "
     "+o3uvw*1.082*zo3(:)/((1.+13"):
        "the clear-sky band-1 transmissivity for the upward beam, on the "
        "pass's 90-character key. On a lane losun discards it is exactly 1: "
        "zo3 is preset to 0 and zo3t and zxo3t are zero on every lane by the "
        "block above the absorber loop. On a lane losun keeps it is the same "
        "expression the downward beam stores as a VALUE, so the mask does not "
        "decide it and this is not a masked-domain site. world-2223 carries "
        "whether swr's clear-sky transmissivities can be driven to zero",
    ("radmod.f90", "subroutine swr", "divide",
     "1.-co2sww*(zca1*LOG(1.+zcb1*zco2(:)) "
     "+zca2*LOG(1.+zcb2*zco2(:)))/zsolar2"):
        "the clear-sky CO2 transmissivity for the upward beam, exactly the "
        "ozone row's shape: 1 on a discarded lane because zco2, zco2t and "
        "zyco2t are all zero there, and the same expression as the stored "
        "downward value on a kept one. world-2223",
    ("radmod.f90", "subroutine swr", "divide",
     "1.-h2osww*h2oswl*2.9*zwv(:) "
     "/((1.+141.5*zwv(:))**0.635+5.925*zwv(:)) /zsolar2"):
        "the clear-sky water vapour transmissivity for the upward beam, "
        "exactly the ozone row's shape: 1 on a discarded lane because zwv, "
        "zwvt and zywvt are all zero there, and the same expression as the "
        "stored downward value on a kept one. world-2223",
    ("radmod.f90", "subroutine swr", "divide",
     "1.-zr1s(:,jlev)*zrl1s(:,jlev)"):
        "the adding method's 1 - R_above*R_below, the geometric series of "
        "round trips between two layers. Both factors are locals now preset to "
        "0 on every lane, so a discarded lane divides by exactly 1; on a kept "
        "lane both are layer reflectivities and the mask has no bearing on "
        "either. world-2223 carries whether the scheme bounds them",
    ("radmod.f90", "subroutine swr", "divide",
     "1.-zr2s(:,jlev)*zrl2s(:,jlev)"):
        "the band-2 twin of the row above: both factors preset to 0, so a "
        "discarded lane divides by exactly 1. world-2223",
    ("radmod.f90", "subroutine swr", "divide", "1.-zra1s(:)*zrb1s(:,jlev)"):
        "the same adding-method denominator in the downward and upward loops, "
        "with the combined layer above in place of the interface value: both "
        "factors preset to 0, so a discarded lane divides by exactly 1. "
        "world-2223",
    ("radmod.f90", "subroutine swr", "divide", "1.-zra2s(:)*zrb2s(:,jlev)"):
        "the band-2 twin of the row above: both factors preset to 0, so a "
        "discarded lane divides by exactly 1. world-2223",
    ("radmod.f90", "subroutine swr", "divide", "1.0+ztcon(:)"):
        "ztcon is b*tau/mu in the conservative-scattering limit: a "
        "non-negative mixture backscatter ratio times an optical depth the "
        "mask makes positive, over a beam cosine that is floored at 1.E-6 or "
        "is the 0.5 diffusivity cosine. It is preset to 0 besides, so the "
        "divisor is at least 1 on every lane",
    ("radmod.f90", "subroutine swr", "divide", "ALOG(3.+0.1*ztau2(:))"):
        "ztau2 is preset to 1 and is otherwise a product of non-negative "
        "factors, so the argument is at least 3 and the logarithm at least "
        "ln 3",
    ("radmod.f90", "subroutine swr", "divide", "EXP(zaertf1(:,jlev))"):
        "zaertf1 is preset to 0 and, where evaluated, is MIN(25., a "
        "non-negative), so the exponential lies in [1, e**25]",
    ("radmod.f90", "subroutine swr", "divide", "EXP(zaertf2(:,jlev))"):
        "zaertf2 is preset to 0 and, where evaluated, is MIN(25., a "
        "non-negative), so the exponential lies in [1, e**25]",
    ("radmod.f90", "subroutine swr", "divide", "dt(:,jlev)"):
        "the mask is losun and the divisor is the gridpoint temperature, which "
        "the mask has no bearing on and which is positive on every lane; this "
        "is the sqrt row above, as a division",
    ("radmod.f90", "subroutine swr", "divide", "ga"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "gravity, which config/planet.yaml declares positive",
    ("radmod.f90", "subroutine swr", "divide", "qex1(jaer)"):
        "radini aborts on a species whose band-1 extinction efficiency is not "
        "positive, and the loop that divides by it runs over exactly the "
        "species radini validated",
    ("radmod.f90", "subroutine swr", "divide", "zaerd1(:,jlev)"):
        "the aerosol two-stream denominator (u+1)**2 e**t - (u-1)**2 e**-t, "
        "with u floored at 1 and t at or above 0 by the guards argued at its "
        "site, so it is at least 4u; preset to 1 besides",
    ("radmod.f90", "subroutine swr", "divide", "zaerd1s(:,jlev)"):
        "the diffuse-beam twin of the row above, floored and preset the same "
        "way",
    ("radmod.f90", "subroutine swr", "divide", "zaerd2(:,jlev)"):
        "the band-2 twin of zaerd1, floored and preset the same way",
    ("radmod.f90", "subroutine swr", "divide", "zaerd2s(:,jlev)"):
        "the band-2 diffuse twin of zaerd1, floored and preset the same way",
    ("radmod.f90", "subroutine swr", "divide", "zexp(:)"):
        "exp of MIN(25., a non-negative exponent), because the optical depth, "
        "the square root beside it and both beam cosines are non-negative, so "
        "zexp is at least 1; preset to 1 besides",
    ("radmod.f90", "subroutine swr", "divide", "zmu0(:)**1.7+0.065"):
        "zmu0 is the beam cosine, non-negative on every lane, so the divisor "
        "is at least 0.065 and the power is in domain",
    ("radmod.f90", "subroutine swr", "divide", "zmu0(:)+0.15"):
        "zmu0 is the beam cosine, non-negative on every lane, so the divisor "
        "is at least 0.15",
    ("radmod.f90", "subroutine swr", "divide", "zmu0+zero"):
        "the floor is written into the expression: zero is a parameter at "
        "1.E-6 and zmu0 is the non-negative beam cosine, so the divisor is at "
        "least 1.E-6",
    ("radmod.f90", "subroutine swr", "divide", "zmu00"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the diffusivity cosine for the scattered beam, set to 0.5 above the "
        "level loop",
    ("radmod.f90", "subroutine swr", "divide", "zr(:)"):
        "the cloud two-stream denominator (u+1)**2 zexp - (u-1)**2/zexp. zexp "
        "is at least 1 and u at or above 0, so zr is at least 4u and vanishes "
        "only where u does, which is where zuz is at or below 0. zuz is "
        "zun + 2*zb2*zom0 = 1 - zom0*(1-2*zb2) with zom0 capped at 0.9999 and "
        "zb2 non-negative, so it is at least 1.E-4 whenever 2*zb2 is at most "
        "1; zb2 is at most tswr2/ALOG(3.), which at the 0.065 this project "
        "leaves tswr2 at is 0.059. Preset to 1 besides",
    ("radmod.f90", "subroutine swr", "divide", "zsolar1"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the fraction of the stellar flux below 0.75 um, an integral of the "
        "Planck function over a band divided by the total and therefore "
        "strictly between 0 and 1",
    ("radmod.f90", "subroutine swr", "divide", "zsolar2"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the fraction of the stellar flux above 0.75 um, strictly between 0 "
        "and 1 for the reason zsolar1 is",
    ("radmod.f90", "subroutine swr", "divide", "ztco2t(:)"):
        "a running clear-sky CO2 transmissivity, preset to 1 on every lane, so "
        "a discarded lane divides by exactly 1. On a kept lane it is "
        "1 - A(u)/zsolar2 for the column above the level, which the mask does "
        "not decide. world-2223",
    ("radmod.f90", "subroutine swr", "divide", "ztco2u(:)"):
        "the upward-beam twin of the row above, preset to 1 the same way. "
        "world-2223",
    ("radmod.f90", "subroutine swr", "divide", "zto3t(:)"):
        "a running clear-sky ozone transmissivity, preset to 1 on every lane, "
        "so a discarded lane divides by exactly 1. On a kept lane it is "
        "1 - A(u)/zsolar1 for the column above the level, which the mask does "
        "not decide. world-2223",
    ("radmod.f90", "subroutine swr", "divide", "zto3u(:)"):
        "the upward-beam twin of the row above, preset to 1 the same way. "
        "world-2223",
    ("radmod.f90", "subroutine swr", "divide", "ztwvt(:)"):
        "a running clear-sky water vapour transmissivity, preset to 1 on every "
        "lane, so a discarded lane divides by exactly 1. On a kept lane it is "
        "1 - A(u)/zsolar2 for the column above the level, which the mask does "
        "not decide. world-2223",
    ("radmod.f90", "subroutine swr", "divide", "ztwvu(:)"):
        "the upward-beam twin of the row above, preset to 1 the same way. "
        "world-2223",
    ("radmod.f90", "subroutine swr", "divide", "zun(:)"):
        "1 - zom0 with zom0 capped at 0.9999 by the MIN one line above it, so "
        "the divisor is at least 1.E-4; preset to 1 besides",
    ("radmod.f90", "subroutine swr", "divide", "zwfit"):
        "a named constant, `real, parameter :: zwfit = 10.0`, the bottom of "
        "Stephens's fitted range",

    ("seamod.f90", "subroutine seaini", "divide", "dicealbdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the width of the sea-ice albedo ramp, a seamod scalar at 10 K",
    ("seamod.f90", "subroutine seaini", "divide", "psurf"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the reference surface pressure, which config/planet.yaml declares "
        "positive",
    ("seamod.f90", "subroutine seaini", "divide", "ra4d(dt(:,NLEP),ra4)"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("seamod.f90", "subroutine seaini", "divide", "rdbrv"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the ratio of the dry and vapour gas constants, about 0.622",
    ("seamod.f90", "subroutine seastep", "divide", "1.-(1./rdbrv-1.)*dqs(:)"):
        "1/rdbrv - 1 is about 0.608, so the divisor vanishes only at a "
        "specific humidity of about 1.645 kg/kg, which is above the 1 kg/kg a "
        "mass fraction cannot exceed",
    ("seamod.f90", "subroutine seastep", "divide", "dicealbdt"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the width of the sea-ice albedo ramp, a seamod scalar at 10 K",
    ("seamod.f90", "subroutine seastep", "divide", "dp(:)"):
        "the surface pressure, positive on every gridpoint whatever the mask "
        "selects",
    ("seamod.f90", "subroutine seastep", "divide", "ga*dp(:)"):
        "gravity times the surface pressure, both positive on every gridpoint",
    ("seamod.f90", "subroutine seastep", "divide", "ra4d(dt(:,NLEP),ra4)"):
        "ra4d is `max(pt - ppole, 1.0)` in plasimmod, so the divisor is at "
        "least 1 on every lane",
    ("seamod.f90", "subroutine seastep", "divide", "rdbrv"):
        "a scalar with no lane dependence, so the mask has no bearing on it: "
        "the ratio of the dry and vapour gas constants, about 0.622",
    ("seamod.f90", "subroutine seastep", "divide", "real(naccua)"):
        "the accumulation counter, incremented unconditionally on the "
        "statement above the block that divides by it, so it is at least 1",
}


def squeeze(text: str) -> str:
    """Collapse runs of whitespace, so a continuation's padding is not a key."""
    return " ".join(text.split())


def strip_comment(line: str) -> str:
    """Drop a trailing `!` comment, respecting quoted strings."""
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def logical_lines(path: Path):
    """Yield (first_physical_lineno, joined_text) with continuations folded."""
    raw = path.read_text(errors="replace").splitlines()
    buf = ""
    start = None
    for n, line in enumerate(raw, 1):
        code = strip_comment(line).rstrip()
        if not code.strip():
            if buf:
                continue
            continue
        if buf:
            code = code.lstrip()
            if code.startswith("&"):
                code = code[1:]
        else:
            start = n
        if code.rstrip().endswith("&"):
            buf += code.rstrip()[:-1]
            continue
        buf += code
        yield start, buf
        buf = ""
    if buf:
        yield start, buf


def match_paren(text: str, open_idx: int) -> int:
    """Index of the `)` closing the `(` at open_idx, or -1."""
    depth = 0
    quote = None
    for i in range(open_idx, len(text)):
        ch = text[i]
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
    return -1


def outer_call(text: str, name_start: int) -> str:
    """The token immediately enclosing the call that starts at name_start."""
    i = name_start - 1
    while i >= 0 and text[i] in " \t":
        i -= 1
    if i < 0 or text[i] != "(":
        return ""
    j = i - 1
    while j >= 0 and text[j] in " \t":
        j -= 1
    end = j + 1
    while j >= 0 and (text[j].isalnum() or text[j] == "_"):
        j -= 1
    return text[j + 1 : end].lower()


def arg_of(text: str, name_end: int) -> str:
    close = match_paren(text, name_end - 1)
    return text[name_end:close] if close > 0 else text[name_end:]


def leading_call(expr: str) -> str:
    """The function name a parenthesised expression opens with, if any."""
    m = re.match(r"\s*([a-z_]\w*)\s*\(", expr, re.I)
    return m.group(1).lower() if m else ""


_BARE_DIV_RE = re.compile(r"(?<!/)/\s*([a-z_]\w*)\s*(\()?", re.I)


def divisor_sites(text: str):
    """Yield (divisor_text, index) for every `/` whose right operand is a value.

    Three forms, and the population is the union of them:

      * `/ (expr)`, the parenthesised divisor;
      * `/ name`, a scalar or a whole array;
      * `/ name(...)`, an element, a section, or a function result.

    world-d016 re-derived the population after the parenthesised form alone had
    been read. A division by a bare name is the identical hazard -- the mask
    does not stop the lane being evaluated, and `-ffpe-trap=zero` does not care
    how the divisor was spelled -- so leaving it out made the count look like a
    measurement of the class when it was a measurement of the parser.

    `//` is string concatenation and `(/ ... /)` is an array constructor;
    neither is a division and neither is yielded. A numeric literal divisor is
    not yielded either: it cannot vary by lane, so the mask has no bearing on
    it.
    """
    for m in re.finditer(r"/\s*\(", text):
        open_idx = m.end() - 1
        close = match_paren(text, open_idx)
        if close < 0:
            continue
        yield text[open_idx + 1 : close], open_idx
    for m in _BARE_DIV_RE.finditer(text):
        if text[max(0, m.start() - 1) : m.start() + 1] == "(/":
            continue
        if m.group(2):
            close = match_paren(text, m.end() - 1)
            if close < 0:
                continue
            yield text[m.start() + 1 : close + 1].strip(), m.start()
        else:
            yield m.group(1), m.start()


def scan(path: Path):
    stack: list[str] = []
    proc = "(file scope)"
    sites = []
    for lineno, text in logical_lines(path):
        low = text.lower()

        if _PROC_END_RE.match(low):
            # A procedure boundary clears the stack. Over-reporting means the
            # stack is never cleared EARLIER than the truth, so this is the one
            # place it is cleared without an explicit `end where`.
            stack.clear()
            proc = "(file scope)"
        pm = _PROC_START_RE.match(low)
        if pm:
            stack.clear()
            # The KIND and the NAME, without the dummy argument list. A key is
            # read by a person and an argument list runs to two hundred
            # characters on some of these; it also moves when an argument is
            # added, which would drop a classified site out of the table for a
            # change that has nothing to do with its divisor.
            head = text.strip().split("!")[0].strip()
            name = re.match(
                r".*?\b(subroutine|function|program|module)\s+([a-z_]\w*)",
                head, re.I,
            )
            proc = f"{name.group(1).lower()} {name.group(2)}" if name else head

        if _END_WHERE_RE.match(low):
            if stack:
                stack.pop()
            continue

        single_mask = None
        if _ELSEWHERE_RE.match(low):
            m = re.match(r"^\s*else\s*where\s*\(", low)
            if m:
                close = match_paren(text, m.end() - 1)
                if stack:
                    stack[-1] = ".not. " + text[m.end() : close]
            elif stack:
                stack[-1] = ".not. " + stack[-1]
        elif _WHERE_OPEN_RE.match(low):
            open_idx = low.index("(", low.index("where"))
            close = match_paren(text, open_idx)
            mask = text[open_idx + 1 : close] if close > 0 else "?"
            rest = text[close + 1 :].strip() if close > 0 else ""
            if rest:
                single_mask = mask
            else:
                stack.append(mask)
                continue

        masks = list(stack) + ([single_mask] if single_mask else [])
        if not masks:
            continue

        for m in _INTRINSIC_RE.finditer(text):
            name = m.group(1).lower()
            arg = arg_of(text, m.end())
            guarded = leading_call(arg) in _ARG_GUARDS or outer_call(
                text, m.start()
            ) in _ARG_GUARDS
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind=name,
                    guarded=guarded,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail=squeeze(arg)[:90],
                )
            )
        for m in _MERGE_RE.finditer(text):
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind="merge",
                    guarded=False,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail="",
                )
            )
        for div, _idx in divisor_sites(text):
            guarded = leading_call(div) in _DIV_GUARDS
            sites.append(
                dict(
                    file=path.name,
                    line=lineno,
                    proc=proc,
                    kind="divide",
                    guarded=guarded,
                    mask=" .and. ".join(masks),
                    text=text.strip(),
                    detail=squeeze(div)[:90],
                )
            )
    return sites


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--all",
        action="store_true",
        help="include the guarded sites and the classified survivors",
    )
    ap.add_argument("--tsv", action="store_true", help="one tab-separated row per site")
    ap.add_argument(
        "--kind",
        default="any",
        choices=("intrinsic", "divide", "merge", "any"),
        help="which class to report (default: all three, which is the gate)",
    )
    args = ap.parse_args()

    sites = []
    for path in sorted(SRC.glob("*.f90")):
        sites.extend(scan(path))

    def wanted(s):
        if args.kind == "intrinsic":
            return s["kind"] in DOMAIN_INTRINSICS
        if args.kind == "any":
            return True
        return s["kind"] == args.kind

    sites = [s for s in sites if wanted(s)]
    for s in sites:
        s["classified"] = CLASSIFIED.get(
            (s["file"], s["proc"], s["kind"], s["detail"])
        )
    if not args.all:
        sites = [s for s in sites if not s["guarded"] and not s["classified"]]

    if args.tsv:
        for s in sites:
            print(
                "\t".join(
                    (
                        s["file"],
                        str(s["line"]),
                        s["kind"],
                        "guarded"
                        if s["guarded"]
                        else ("classified" if s["classified"] else "bare"),
                        s["proc"],
                        s["mask"],
                        s["detail"],
                    )
                )
            )
    else:
        by_file: dict[str, int] = {}
        for s in sites:
            by_file[s["file"]] = by_file.get(s["file"], 0) + 1
            print(f"{s['file']}:{s['line']}  {s['kind']}  [mask: {s['mask'][:70]}]")
            print(f"    {s['detail']}")
            if s["classified"]:
                print(f"    classified: {s['classified']}")
        print()
        for f, n in sorted(by_file.items(), key=lambda kv: -kv[1]):
            print(f"{n:5d}  {f}")
        print(f"{len(sites):5d}  TOTAL")
    if args.all:
        return 0
    return 1 if sites else 0


if __name__ == "__main__":
    raise SystemExit(main())
