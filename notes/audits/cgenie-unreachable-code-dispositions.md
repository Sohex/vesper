# Unreachable code in the vendored ocean host, and what decides whether it goes

WORLDBUILDING CONTEXT, stated first because this document borrows vocabulary
from a real discipline: **Vesper is a fictional planet and this is engineering
work on the simulation of it.** Ocean, sea ice, surface flux and atmosphere name
modelled quantities and numerical properties of the adopted offline ocean host,
not observations of anything.

Audited 2026-09-05 against `vendor/cgenie` at `e462b20e4`, and against upstream
`derpycode/cgenie.muffin` master read through the GitHub commit API on the same
date.

## The question, and why it does not have one answer

`notes/audits/cgenie-embm-free-path-and-threading.md` left three things in the
vendored ocean host that no configuration reaches, and the working agreement
asks of each whether it should exist before anything is built or repaired.

The answer is not the same for all three, and the reason is the fork. This
subtree is a MAINTAINED fork: upstream is pulled into `Sohex/cgenie.muffin`
first and `vendor/cgenie` pulls from there, so an upstream commit still arrives.
`notes/audits/dead-code-and-unreachable-paths.md` could ask "should this file
exist" of every file in the climate model because that fork was declared HARD on
2026-08-21 and upstream compatibility stopped being a constraint. No such
declaration covers this subtree, so "upstream might want it" still defends
something here, and the question has to be asked twice rather than once.

## The test, and its measurable half

A thing goes when BOTH hold.

**Nothing here reaches it.** Established by a whole-tree sweep rather than a grep
for the identifier: the makefile object lists, `genie_loop_wrappers.f90` and
every `call` in `genie.F`, `src/xml-config/xml/definition.xml` and the XSLT that
turns a configuration into a namelist, `genie.job` and the shell around it, and
the COMMON block headers that carry a name into every file that includes one.

**And upstream maintenance of it cannot deliver value here.** This is the half
that decides, and it is measurable rather than a judgement: the upstream commit
series for the exact path says how often the file is touched and what the
touches are. Three outcomes recur.

- Maintained, and the maintenance has delivered a repair this project wants.
  KEEP. A deletion forfeits the next such repair and turns it into a subtree
  pull conflict.
- Maintained, but every commit in the series is to a file that cannot run. The
  maintenance is then evidence of NON-use, not of use. DELETE.
- Made dead by a change of this project's own. Upstream never had the deadness
  and has nothing to maintain. DELETE.

## 1. `config2xml.py` is deleted, and the maintenance is what condemned it

`genie-main/config2xml.py` and its only caller `genie-main/configs/config2xml.sh`
are removed. Both were mainline cGENIE, unmodified in this fork: the sole commit
touching either here is the subtree merge `0811a413c`.

**It has never been runnable in the history that exists.** `ast.parse` refuses it
at line 27, a missing comma between two entries of the `configflags` dict, which
is a syntax error under Python 2 and Python 3 alike. ALL THIRTEEN upstream
revisions of the file fail at that same line, from `e143897c` on 2018-05-14, the
import to GitHub, through `72508d72` on 2021-08-18, which is master today. The
twelve commits after the import are the tracer-table series -- vivianite,
greenalite, FeS2, FeOOH, the Os tracers -- so the file has been edited a dozen
times over seven years and executed by nobody, because it cannot be. The
maintenance is evidence of non-use.

**Repairing the comma would not have produced a usable tool.** The body is Python
2, with statement `print` at lines 864 and 963, and this host has no Python 2.
Its `configflags` map is a hand-copy of the model registry that
`src/xml-config/xml/definition.xml` owns authoritatively, and it has drifted:
`ma_flag_plasimatmos` is absent, and so is `ma_flag_fluxatmos`, the atmosphere
flag world-e3gp added for the driven path. A repaired converter would silently
drop both of the atmosphere flags that matter here.

**Deleting it closes no route.** The live old-style-config to XML converter is a
different file, `genie-main/translate_config.py`, wired into `genie.job` at line
33 as `TRANSLATE_CONFIG` and reached by `genie.job -O`.

It is also Python 2 and also does not run on this host, and the consequence
reaches further than the deleted file, so it is stated here rather than left to
be rediscovered. `runmuffin.sh:395` runs `./genie_example.job -O -f <config>`,
so the shipped route from any of the 635 flat `.config` base configurations to a
run goes through `translate_config.py` and is unavailable here for want of a
Python 2 interpreter. **Nothing in this project wants that route.**
`analysis/cgenie_cost.py` composes the XML directly, `<param name="world">`
included, so an arbitrary shipped topography and an arbitrary grid are both
reachable without a converter, and `analysis/cgenie_profile.py` and
`analysis/cgenie_omp.py` build on it. The converters are a legacy input format's
on-ramp, not the model's.

`config2xml.sh` goes with it because it is the only caller, and because the file
list it reads, `genie-main/configs/files`, is not in the tree either.

**What is kept.** The five XML configurations whose `<job author>` attribute
records that `config2xml.py` produced them once. Two of them,
`eb_go_gs_test.xml` and `eb_go_gs_ac_bg_test.xml`, are what
`analysis/cgenie_cost.py`, `analysis/cgenie_profile.py` and
`analysis/cgenie_omp.py` seed every generated configuration from, so the output
of the converter is load-bearing where the converter is not. The author string is
provenance text inside an attribute, not a reference to a file.

**The cost of the deletion, measured.** Upstream touched `config2xml.py` twelve
times in seven years, about 1.7 commits a year, each a mechanical tracer-table
addition. A subtree pull carrying one of those will raise a delete/modify
conflict, and the resolution is to keep it deleted. `config2xml.sh` has one
upstream commit ever, the import.

## 2. `dzu` is out of `/ocn_vars/`, because this project is what made it dead

`dzu(2,maxk)` is removed from `genie-goldstein/src/fortran/ocean.cmn`, from both
the type declaration and the `/ocn_vars/` COMMON list, and its zeroing is removed
from `initialise_goldstein.F`. Two commented-out lines carried the name as well
and go with it: the `crma` continuation of the pre-generalised-grid type
declaration, and the commented `/ocn_vars/` continuation beside the live one.
The two remaining `crma` lines in the header carry no `dzu` and stay. The
identifier now appears nowhere in `genie-goldstein`, so a later grep for it
lands on nothing rather than on a comment.

`velc` was its only reader and writer anywhere in the tree, and
`initialise_goldstein.F` only zeroed it. Threading `velc` over `j` made it the
local `dzu_col`, because one array in a COMMON block is one array for the whole
process and would have been shared between the threads. The COMMON member has
been written by nothing and read by nothing since. This is the third outcome of
the test: upstream never had a dead `dzu`, so there is no upstream maintenance of
the deadness to forfeit.

**Why the layout change is safe.** `/ocn_vars/` is declared in exactly one place,
`ocean.cmn:75`, and the 39 files that include the header all take the new layout
together. No `EQUIVALENCE` statement appears anywhere in `genie-goldstein` or
`genie-embm`, and no routine writes the block as a unit: every restart and
netCDF write in the component names its arrays individually.

**Acceptance, and it is the bar the threading was held to. It passes.** Every
float variable of both shipped regression cases is BIT-FOR-BIT what the tree at
`e462b20e4` writes, at `worst_relative` 0.0 with no variable differing: all 28
variables of `eb_go_gs`'s GOLDSTEIN year-20 annual average at one thread AND at
sixteen, and all 88 of `eb_go_gs_ac_bg`'s BIOGEM three-dimensional fields.
`analysis/cgenie_omp.py` is the driver, the reference arm is the `base` arm built
from `e462b20e4` rather than the shipped netCDF, and both arms are built in an
exported tree outside the repository.

Against the shipped `genie-knowngood/` netCDF the removal reproduces the two
differences the unedited tree already has and no others: two cells of `uvel` at
6.86e-27 of the field range, and `bio_fpart_CaCO3_13C` and `phys_u` at
2.99e-45. Those are the same two the threading note records, so the comparison
against the shipped reference is unchanged as well.

The bar is bit-for-bit, so it does not move with what else the host is doing,
and no wall clock from these arms is kept.

**Also checked, and deliberately left in place.** `genie-embm/src/fortran/embm.cmn`
declares its own `dzu(2,maxk)` in `/embm_vars/`, a different COMMON block, and it
is dead in the same way: the only assignments to it, `initialise_embm.F:587-588`,
are commented out, and nothing reads it. That one is mainline deadness rather
than deadness this project created, EMBM is not threaded so the shared-storage
trap does not apply to it, and `embm.cmn` is one of the more actively maintained
headers upstream at ten commits with the last on 2024-03-04. It fails the second
half of the test and stays. Recorded so it is not re-found and re-argued.

## 3. `surf_ocn_sic` and its two companions are KEPT

`genie-goldstein/src/fortran/surf_ocn_sic.F` (822 lines),
`outm_surf_ocn_sic.f` (254) and `gold_ocnsic_avg.F` (105) stay, with
`surf_ocn_sic_wrapper` in `genie-main/genie_loop_wrappers.f90` and their three
rows in `genie-goldstein/src/fortran/makefile`. All of it is mainline cGENIE,
unmodified in this fork.

**Nothing reaches it, and the sweep is wider than the earlier one.** No line of
`genie.F` calls `surf_ocn_sic_wrapper`; the wrapper's only two mentions anywhere
are its own definition and its own `end subroutine`. The atmosphere it was
written for is absent: there is no `genie-igcm` module in the tree and
`MODULE_NAMES` in `genie-main/makefile` cannot name one. `outm_surf_ocn_sic` is
called only from `surf_ocn_sic:748`.

### 3a. `gold_ocnsic_avg` is called by nothing at all, including `surf_ocn_sic`

Sharper than the earlier reading, which had it existing to serve `surf_ocn_sic`.
The only three occurrences of the name in the whole subtree are two comment lines
and its own `subroutine` statement. It is a compiled orphan one level deeper than
the routine it is named after, and upstream has one commit on it ever, the 2018
import.

### 3b. What keeping it costs, measured against what would notice

The static storage is the only cost that could bite, because every field cGENIE
holds is in a named COMMON block sized from the grid macros, and past roughly a
gigabyte of it the default `-mcmodel=small` cannot reach the data and the link
fails.

Twenty arrays in `genie-main/genie_global.f90` are referenced by the declaration
and by `surf_ocn_sic_wrapper` and by nothing else: `rough_ocn`, `albavg_ocn`, the
four `ocean_lowestlu2/v2/u3/v3_ocn`, the four `ocean_lowestlt/lq/lp/lh_ocn`,
`ocean_atm_netsolar_ocn`, `ocean_atm_netlong_ocn`, `ocean_netsolar_ocn`,
`ocean_netlong_ocn`, `ocean_latent_ocn`, `ocean_sensible_ocn`, `ocean_evap_ocn`,
`atmos_latent_ocn`, `atmos_sensible_ocn` and `atmos_evap_ocn`. One integer,
`istep_gsurf`, goes with them. Each array is `ilon1_ocn` by `ilat1_ocn`, which is
`GOLDSTEINNLONS` by `GOLDSTEINNLATS`, at eight bytes an element because the
gfortran block compiles `FLAGR8=-fdefault-real-8`.

That is 207 kB at 36 x 36, 829 kB at 72 x 72 and 3.3 MB at 144 x 144. Against the
ceiling it would have to move, 3.3 MB is three parts in a thousand. The effect is
smaller than the instrument, so there is no memory argument for the deletion and
none is claimed. The remaining cost of keeping it is one object file per build.

### 3c. What deleting it would have cost, and this is what decides

**Upstream maintains it in exactly the dimension this project cares about.**
`surf_ocn_sic.F` has two upstream commits, and the second, `72927373` on
2024-01-15, replaces four occurrences of a hardcoded `3600.0*24.0` with
`sodaylen`. That is Earth's day length coming out of a scale factor and a
namelist quantity going in: the implicit-Earth class
`notes/audits/ocean-tier-implicit-earth.md` and OCN-12 exist to find, arriving
free from upstream. A deletion forfeits the next one and converts it into a pull
conflict.

**Three of this project's documents cite the file as evidence.**
`notes/external-model-survey.md` section 27 cites it for the finding that cGENIE
has no sea-surface temperature or salinity restoring, which is half of OCN-5's
requirement satisfied by construction. `notes/audits/ocean-tier-implicit-earth.md`
cites `surf_ocn_sic.F:475-477` beside
`genie-goldsteinseaice/src/fortran/surflux_goldstein_seaice.F:159` for
GOLDSTEIN's salinity-dependent freezing point being computed correctly.
`notes/audits/cgenie-embm-free-path-and-threading.md` section 1c is where the
routine was first read. Each finding survives on its other citation, but the
routine is part of the evidence base rather than inert weight, and the deletion
would have made three references unresolvable without archaeology.

### 3d. The hazard the routine actually poses, and where it is answered

The reason to want it gone is that its name makes it look like the answer to
"how do I drive GOLDSTEIN without an atmosphere". It is not, and the reason is
settled: it takes atmospheric STATE and computes the turbulent fluxes from bulk
formulae on the OCEAN grid, and derives the wind stress from the lowest-level
winds itself, which is the coarse-grid flux evaluation the offline exchange
exists to avoid. A surface flux is nonlinear in the state.

The answer to that question is `flag_fluxatmos`, which now exists, is registered
in `definition.xml`, calls `surflux_goldstein_seaice`, and REFUSES at
`initialise_genie.F` with a message naming the missing supply rather than
integrating an ocean forced by zeros. A reader who reaches for `surf_ocn_sic`
finds nothing calling it; a reader who reaches for the flag finds a refusal that
says what is missing. The hazard is a naming hazard and it is answered by a
recorded verdict, which is cheaper than a permanent divergence from a community
model.
