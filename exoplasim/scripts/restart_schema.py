"""What every restart record IS: shape, type, meaning and conversion policy.

Worldbuilding frame: these are the saved fields of the Vesper climate model.
Nothing here is about the real world.

WHY A SCHEMA AT ALL. The restart file carries no type tag, no shape and no
version. A payload length divisible by eight is not evidence of a
double-precision array -- an integer array and a four-byte real array of the
same element count have the same length -- so nothing can be converted by
inspecting the bytes. Every decision this module makes is driven by a name.

TWO HALVES, DELIBERATELY SEPARATE.

`inventory_from_source()` reads the `put_restart_*`, `mpputgp` and `mpputsp`
call sites out of the model source and returns the MECHANICAL facts: which
names a build can emit, from which module, with which symbolic shape. It is
derived rather than typed out, so it cannot drift from the model.

`POLICY` is the SEMANTIC half, and it is checked in because no parser can read
it out of Fortran: whether a field is prognostic or derived, whether remapping
it conserves an inventory or averages an intensity, and whether it survives a
resolution change at all. It was reviewed once against the declarations and the
comments beside them.

`check_policy_covers_source()` holds the two halves together, and it can fail
in the direction that matters: add a record to the model and forget its policy,
and the check goes red rather than the converter guessing. It also holds the
reset column against `outreset` and its per-module equivalents, which is what
stops "an accumulator's clean value is zero" from quietly becoming false -- it
already is for four of them.

THE RESUME COLUMN is the third thing POLICY carries, and it is a DECISION rather
than a description: what should happen when a restart does not carry the record.
`read_sites_from_source()` is its mechanical half -- which readers touch a name,
whether any read is inside a lowered `nexcheck`, whether anything asks for it
with `has_restart_array` -- and `check_resume_policy()` holds the two against
each other in both directions, so a decision the code does not implement and a
handling nobody decided both go red. world-t16e;
`notes/audits/absent-restart-records.md` argues it per record.

IF YOU ADD OR REMOVE A RESTART RECORD, THIS FILE IS THE ONE TO EDIT, and
`scripts/smoke_test.py` is what fails until you have. `exoplasim/README.md`
carries the contract under "If you add or remove a restart record"; the tools
that read this table are `convert_restart.py`, `reset_restart_accumulators.py`
and `restart_surface.py`.
"""
from __future__ import annotations

import functools
import re
import struct
from dataclasses import dataclass, replace as _replace
from pathlib import Path

import _paths

# The vendored model this module reads shapes and constants out of, by default.
# Every source-reading function here still takes an explicit `src_dir`, so a
# caller working against another tree passes one; this is only what
# `describe()` uses when it is handed a file and nothing else.
MODEL_SRC = _paths.MODEL_SRC / "plasim" / "src"

# ---------------------------------------------------------------------------
# The mechanical half: what the model source says it writes.
# ---------------------------------------------------------------------------

# Which sources `plasim.x` is built from is read out of the build itself rather
# than listed here, so a file that sits in `plasim/src` and is in no source list
# cannot contribute a record this executable emits, and no hand list has to keep
# saying which those are. The configurable slots are expanded to every value the
# build file offers, because a record or a reset must be the same across the
# configurations rule 4 counts binaries over.
CMAKE = "CMakeLists.txt"
_SOURCES_BLOCK = re.compile(r"set\(_sources(?P<body>.*?)\)", re.DOTALL)
_CMAKE_SET = re.compile(r"^\s*set\(\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s+"
                        r"(?P<value>[A-Za-z0-9_.\"]+)", re.MULTILINE)
# `set()` gives a cache slot ONE value, its default. `require_one_of` is where
# the admitted values are, and a slot whose alternatives are only there would
# otherwise be read as having none: PLASIM_FFT defaults to fftmod and
# `lib/rungs.py` picks fft991mod for the longitude counts fftmod cannot factor.
_CMAKE_REQUIRE = re.compile(r"^\s*require_one_of\(\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s+"
                            r"\S+\s+(?P<values>[^)]*)\)", re.MULTILINE)


def compiled_modules(plasim_dir: Path, suffixes: tuple = (".f90",)) -> tuple:
    """Every source the model executable is built from, in build order.

    Parsed from `plasim/CMakeLists.txt`, which is what `build_model.py` drives.
    The build this replaced staged a second copy of every source under
    `plasim/bld/` and drifted from `plasim/src/`; that directory is gone, and
    reading the build file rather than any staged copy is what keeps it gone.

    `suffixes` widens it past Fortran. The restart scan wants the `.f90` set
    alone; `rebuild_binaries.model_sources` wants the C stub too, because it is
    asking what a binary was compiled FROM rather than what can emit a record.
    """
    text = (Path(plasim_dir) / CMAKE).read_text(encoding="utf-8")
    block = _SOURCES_BLOCK.search(text)
    if block is None:
        raise ValueError(f"{plasim_dir / CMAKE} has no _sources list; the "
                         "build file has been restructured")
    choices: dict = {}
    for m in _CMAKE_SET.finditer(text):
        choices.setdefault(m["name"], set()).add(m["value"].strip('"'))
    for m in _CMAKE_REQUIRE.finditer(text):
        choices.setdefault(m["name"], set()).update(m["values"].split())
    out = []
    for token in block["body"].split():
        token = token.strip('"')
        if not token.endswith(suffixes):
            continue
        stem = token.rsplit(".", 1)[0]
        if "/" in stem:
            # resmod.f90 is generated into the build directory from the
            # requested geometry. It holds no restart call and no reset.
            continue
        var = re.fullmatch(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}", stem)
        if var is None:
            out.append(token)
            continue
        suffix = token[len(stem):]
        for value in sorted(v for v in choices.get(var[1], ()) if v):
            if value and not value.startswith("$"):
                out.append(f"{value}{suffix}")
    return tuple(dict.fromkeys(out))

# `put_restart_array(yn,pa,k1,k2,k3)` writes `pa(1:k1,1:k3)`, so the WRITTEN
# length is k1*k3 and the declared second dimension k2 never reaches the file.
# `mpputgp(yn,p,NHOR,klev)` gathers to NUGP*klev; `mpputsp(yn,p,NSPP,klev)`
# gathers to NRSP*klev.
_CALL = re.compile(
    r"""^(?P<lead>[^!\n]*?)
        call\s+(?P<fn>put_restart_integer|put_restart_real|put_restart_seed
                    |put_restart_array|mpputgp|mpputsp)
        \s*\(\s*'(?P<name>[A-Za-z0-9_]+)'\s*,(?P<rest>.*)$""",
    re.IGNORECASE | re.VERBOSE)
# A same-line guard: `if (nqspec == 1) call mpputsp(...)`.
_GUARD = re.compile(r"if\s*\((?P<cond>[^)]*)\)\s*$", re.IGNORECASE)


def _split_args(rest: str) -> list[str]:
    """The remaining actual arguments of a call, by depth-aware comma split."""
    args, depth, cur = [], 0, []
    for ch in rest:
        if ch == "(":
            depth += 1
        elif ch == ")":
            if depth == 0:
                break
            depth -= 1
        if ch == "," and depth == 0:
            args.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    args.append("".join(cur).strip())
    return [a for a in args if a]


def _dim(token: str):
    """A dimension as an int where it is a literal, else its symbol.

    A PRODUCT STAYS ONE SYMBOL, `NLEV*28`, which is how `adener3d(NHOR,NLEV,28)`
    reaches `mpputgp` as one gridpoint record. `Geometry.resolve` splits it and
    multiplies the factors. It is kept as a string rather than expanded into a
    tuple because a shape is a sequence of symbols that gets sorted and printed,
    and a tuple inside one is neither orderable against the rest nor readable.
    """
    token = token.strip()
    return int(token) if token.isdigit() else token.upper()


@dataclass(frozen=True)
class SourceRecord:
    """One `put_*` call site, as the source states it."""

    name: str
    module: str
    line: int
    writer: str
    # The actual argument the call site passes: `aasqsp` is written from
    # `aasqout` and the six spectral accumulators are all named this way, so
    # anything looking for the variable cannot look for the record name.
    variable: str
    # Every shape the source can write this name with: ('NUGP', 'NLEV'),
    # ('NRSP', 1), ('NESP', 'NLEV'), (2, 1), or () for a scalar. More than one
    # where a switch chooses between them -- `dq` is the surface level under
    # spectral humidity and the whole column under semi-Lagrangian.
    shapes: tuple
    guard: str | None


def inventory_from_source(src_dir: Path) -> dict[str, SourceRecord]:
    """Every restart name `plasim.x` can emit, from the call sites themselves."""
    src_dir = Path(src_dir)
    found: dict[str, SourceRecord] = {}
    for module in compiled_modules(src_dir.parent):
        path = src_dir / module
        if not path.is_file():
            raise FileNotFoundError(f"{path} is missing; the model source moved")
        for lineno, line in enumerate(path.read_text(encoding="utf-8",
                                                     errors="ignore").splitlines(), 1):
            m = _CALL.match(line)
            if m is None:
                continue
            fn, name, args = m["fn"].lower(), m["name"], _split_args(m["rest"])
            variable = re.split(r"[(\s]", args[0])[0] if args else name
            if fn == "put_restart_integer":
                shape = ()
            elif fn == "put_restart_real":
                shape = ()
            elif fn == "put_restart_seed":
                shape = ("NSEEDLEN",)
            elif fn == "put_restart_array":
                shape = (_dim(args[1]), _dim(args[3]))
            elif fn == "mpputgp":
                shape = ("NUGP", _dim(args[2]))
            else:                                    # mpputsp
                shape = ("NRSP", _dim(args[2]))
            g = _GUARD.search(m["lead"].strip())
            if name in found:
                # A name written from more than one branch keeps every shape
                # it can take; the first call site keeps the attribution.
                prior = found[name]
                if shape not in prior.shapes:
                    found[name] = SourceRecord(
                        name=name, module=prior.module, line=prior.line,
                        writer=prior.writer, variable=prior.variable,
                        shapes=prior.shapes + (shape,), guard=prior.guard)
                continue
            found[name] = SourceRecord(
                name=name, module=module, line=lineno, writer=fn,
                variable=variable, shapes=(shape,),
                guard=g["cond"].strip() if g else None)
    return found


# ---------------------------------------------------------------------------
# The semantic half: what the converter does with each record.
# ---------------------------------------------------------------------------

# Semantic classes, as the specification names them.
HEADER = "header"
TIME_AND_PHASE = "time-and-phase"
PROGNOSTIC_SPECTRAL = "prognostic-spectral"
STATIC_SPECTRAL = "static-spectral"
PROGNOSTIC_GRID = "prognostic-grid"
STATIC_GRID = "static-grid"
DERIVED_GRID = "derived-grid"
ACCUMULATOR = "accumulator"
OPAQUE = "opaque-build-state"
CONTROLLER = "controller-state"    # a correction the model is tracking, not a
                                   # state variable and not an accumulator

# Conversion actions.
TARGET = "target"                # take the target template's value
COPY = "copy"                    # carry the source value across unchanged
REQUIRE_EQUAL = "require_equal"  # source and target must agree, or refuse
PROJECT = "project"              # triangular spectral projection
REMAP = "remap"                  # conservative Gaussian remap
RESET = "reset"                  # the target template's clean reset value
RECOMPUTE = "recompute"          # target placeholder; the model rebuilds it
SEED = "seed"                    # the compiler's random-number state

# Gridpoint remap policies, from the physical meaning of the field.
INTENSIVE = "intensive"          # area-weighted mean; a temperature, a rate
RESERVOIR = "reservoir"          # per-area stock; the same operator, but the
                                 # global inventory is reported and must close
FRACTION = "fraction"            # bounded 0..1 cover; remap then clip
THICKNESS = "thickness"          # remap the VOLUME with its partner fraction
GEOMETRIC = "geometric"          # recomputed from the target grid


@dataclass(frozen=True)
class Policy:
    semantic: str
    action: str
    remap: str | None = None
    # Land, ocean or global. Land and ocean are remapped separately under the
    # target mask so a coastline change cannot average a soil temperature into
    # a sea surface temperature.
    domain: str = "global"
    bounds: tuple | None = None
    # The partner record for a fraction/thickness pair, and the inventory a
    # reservoir has to close against.
    partner: str | None = None
    conserve: str | None = None
    # The fraction this quantity is measured PER, where that is not the whole
    # cell. An ice thickness is a depth over the ice-covered part, so its
    # inventory is thickness x cover x area, and comparing thickness x area
    # before and after says a volume-conserving remap moved the volume by 634%.
    measured_per: str | None = None
    # True where SIMBA takes ownership of the field under coupled vegetation.
    vegetation_owned: bool = False
    # True where the model REBUILDS this record from other state before it uses
    # it again, so whatever a converted restart carries survives one timestep.
    # It is not a licence to put anything there: between the template's value
    # and the donor's own, the donor's is consistent with the prognostics that
    # were converted and the template's is consistent with nothing in the file.
    # So these are REMAPPED like any other field and merely NAMED here, which
    # is what the report's expected_to_change_in_model_fixup list is built from.
    rebuilt_by_model: bool = False
    # For an accumulator, what `outreset` and its equivalents do to it at an
    # interval boundary: "zero", "sentinel" for one whose clean value is not
    # zero, or "none" for one the model deliberately never resets. Anything
    # that writes a clean accumulator has to know this: zeroing a running
    # minimum makes its minimum zero forever.
    model_reset: str | None = None
    reset_value: float | None = None
    # WHAT HAPPENS WHEN A RESTART DOES NOT CARRY THIS RECORD, and it is a
    # decision per record rather than a property of the reader. `required` means
    # the model stops with the record's name printed, which is the right answer
    # wherever the alternative is resuming from a state that is not the state --
    # rule 7 says every run is disposable until the canonical climatology
    # lineage is declared, so refusing an old restart costs a re-commissioning
    # and nothing else. `optional` means the model resumes without it, from
    # `resume_cold`, and is only defensible where what is lost is a partial
    # accumulation window. `resume_marker` names the version record the read
    # sits behind where one does. world-t16e; the enumeration and the
    # per-record argument are in notes/audits/absent-restart-records.md.
    resume: str = "required"
    resume_marker: str | None = None
    resume_cold: str = ""
    why: str = ""


def _acc(names):
    return {n: Policy(ACCUMULATOR, RESET, model_reset="zero",
                      why="accumulated over the output window; the window "
                          "restarts at the target") for n in names}


POLICY: dict[str, Policy] = {}

# --- headers and geometry --------------------------------------------------
POLICY.update({
    "nlat": Policy(HEADER, TARGET, why="the target grid's own"),
    "nlon": Policy(HEADER, TARGET, why="the target grid's own"),
    "nlev": Policy(HEADER, REQUIRE_EQUAL,
                   why="v1 requires equal vertical dimensions"),
    "nrsp": Policy(HEADER, TARGET, why="follows from the target truncation"),
    "nlsoil": Policy(HEADER, REQUIRE_EQUAL,
                     why="soil layers are a separate conversion feature"),
    "nlev_oce": Policy(HEADER, REQUIRE_EQUAL,
                       why="ocean layers are a separate conversion feature"),
    "nstep": Policy(
        TIME_AND_PHASE, COPY,
        why="not only output naming: the calendar and the stellar cycle take "
            "their phase from this count, so resetting it would move the "
            "planet in its orbit. Sound only because v1 requires equal "
            "timesteps."),
})

# --- opaque build and configuration state ----------------------------------
POLICY.update({
    "seed": Policy(OPAQUE, SEED,
                   why="the compiler's random-number state; its length is a "
                       "property of the executable, not of the grid"),
    "nicec2d": Policy(OPAQUE, REQUIRE_EQUAL,
                      why="a physics switch: ice thickness from cover"),
    "accuvers": Policy(OPAQUE, TARGET,
                       why="marks a restart whose accumulator set is complete"),
    "fixedlon": Policy(OPAQUE, REQUIRE_EQUAL,
                       why="the fixed substellar longitude is configuration"),
    "zsolars": Policy(OPAQUE, TARGET,
                      why="the two-band split of the stellar constant, written "
                          "by radstop as a configuration fingerprint. NOTHING "
                          "READS IT BACK: radstart recomputes the split from "
                          "the namelist every time, so the stored value can "
                          "only ever record what the donor run was configured "
                          "with. TARGET and not REQUIRE_EQUAL for exactly that "
                          "reason -- a converter that refused over it would be "
                          "refusing on a value the model is about to discard, "
                          "and converting a restart onto a new stellar "
                          "configuration is a thing this project does"),
})

# --- controller state ------------------------------------------------------
POLICY.update({
    "denergyfix": Policy(
        CONTROLLER, TARGET,
        why="the energy fixer's integrated correction, carried across a "
            "restart so a segment does not begin uncorrected (world-fsr). "
            "TARGET rather than COPY across a CONVERSION: it is a heating rate "
            "the controller reached against one truncation's own conversion "
            "defect, and that defect is resolution-dependent, so a T21 "
            "correction is not the T85 one. The target starts from its "
            "template's value and re-converges in one window."),
})

# --- accumulation counters -------------------------------------------------
POLICY.update(_acc(["naccuout", "naccuice", "naccuo", "naccuoce", "naccua"]))

# --- spectral state --------------------------------------------------------
for _n, _what in (("sz", "vorticity"), ("sd", "divergence"),
                  ("st", "temperature"), ("sq", "specific humidity"),
                  ("sp", "log surface pressure")):
    POLICY[_n] = Policy(PROGNOSTIC_SPECTRAL, PROJECT,
                        why=f"spectral {_what}, current leapfrog level")
for _n, _what in (("szm", "vorticity"), ("sdm", "divergence"),
                  ("stm", "temperature"), ("sqm", "specific humidity"),
                  ("spm", "log surface pressure")):
    POLICY[_n] = Policy(PROGNOSTIC_SPECTRAL, PROJECT,
                        why=f"spectral {_what}, previous leapfrog level; "
                            "projected with its partner or the represented "
                            "time derivative changes")
POLICY["so"] = Policy(STATIC_SPECTRAL, TARGET,
                      why="spectral orography belongs to the target's own "
                          "terrain, scaled and fitted at its truncation")
POLICY["sr"] = Policy(STATIC_SPECTRAL, TARGET,
                      why="the restoration state is configuration")

# --- spectral accumulators -------------------------------------------------
POLICY.update(_acc(["aasosp", "aaspsp", "aastsp", "aasqsp", "aasdsp", "aaszsp"]))
# --- accumulated orbital scalars -------------------------------------------
POLICY.update(_acc(["aorbnu", "alambm", "azdecl", "ardist", "arasc"]))

# --- gridpoint accumulators ------------------------------------------------
POLICY.update(_acc([
    "aprl", "aprc", "aprs", "aevap", "ashfl", "alhfl", "aroff", "asmelt",
    "asndch", "acc", "assol", "asthr", "atsol", "atthr", "ataux", "atauy",
    "atsolu", "assolu", "asthru", "aqvi", "atsa", "ats0", "atsama", "atsami",
    # WORLD-3QFZ: the surface downward solar flux by shortwave band. Written
    # every time and read only behind accuvers >= 2.0, so a restart that
    # predates them takes outreset rather than a diluted first window.
    "afdsw1", "afdsw2",
    # world-5qy: PlaSim's 28-term energy decomposition accumulated over the
    # output window, and the same set resolved per level under nener3d. They
    # are written only where nenergy and nener3d allocate them, so their
    # presence is asked for by name rather than promised by accuvers.
    "adenergy", "adener3d",
    "azmuz", "asigrain", "tempmax", "tempmin",
    "agpi", "aventi", "alaav", "ampoti", "avrmpi", "acapen", "alnb", "achim",
    "aadq", "aammr", "aanrho", "aadmld", "aadt", "aadwatc", "aadsnow", "aadql",
    "aadust3", "aadcc", "aadtd5", "aadls", "aadz0", "aadalb", "aadsalb1",
    "aadsalb2", "aadtsoil", "aadtd2", "aadtd3", "aadtd4", "aadicec", "aadiced",
    "aadforest", "aadwmax", "aadglac", "aadqo3", "aagroundoro", "aaglacieroro",
    "xflxicea", "xheata", "xofluxa", "xqmelta", "xcfluxa", "xcfluxra",
    "xcfluxna", "xsmelta", "ximelta", "xtsfluxa", "xfluxca", "xscflxa",
    "xcpmea", "xcroffa", "xstoia",
    "yheata", "yfssta", "yifluxa", "ydssta", "yfldoa", "yqhda",
    "cheata", "cpmea", "cprsa", "croffa", "ctauxa", "ctauya", "cust3a",
    "cshfla", "cshdta", "clhfla", "clhdta", "cswfla", "clwfla",
]))

# --- static boundary state: the target's own -------------------------------
# These are the restart's copies of staged surface fields. `restart_surface.py`
# is the module that says which staged code lands in which record, and its
# whole subject is that a donor's copy must never supersede a newly staged one.
POLICY.update({
    "dls": Policy(STATIC_GRID, TARGET, why="land-sea mask, staged code 172"),
    "dz0clim": Policy(STATIC_GRID, TARGET, why="roughness, staged code 173"),
    "dz0climo": Policy(STATIC_GRID, TARGET,
                       why="the orographic part of roughness, from the "
                           "target's own terrain"),
    "dalbcl": Policy(STATIC_GRID, TARGET, why="background albedo, code 174"),
    "dalbcl1": Policy(STATIC_GRID, TARGET, why="albedo below 0.75 um, code 175"),
    "dalbcl2": Policy(STATIC_GRID, TARGET, why="albedo above 0.75 um, code 176"),
    "dforest": Policy(STATIC_GRID, TARGET, vegetation_owned=True,
                      why="forest cover, staged code 212"),
    "dwmax": Policy(STATIC_GRID, TARGET, vegetation_owned=True,
                    why="field capacity, staged code 229"),
    "dsoilwfc": Policy(STATIC_GRID, TARGET,
                       why="WORLD-VJBZ. The capacity of each land water layer "
                           "as a fraction of dwmax, per cell, staged code "
                           "2290. A STATIC_GRID taken from the TARGET on the "
                           "same terms dwmax is: both are the same integral of "
                           "the land column property contract over the same "
                           "profile, one summed and one cut, so a conversion "
                           "that took dwmax from the target template and this "
                           "from the source would install a split belonging to "
                           "another soil. Written every step whatever "
                           "nlandwcol is; a restart that lacks the record "
                           "falls back to the namelist dsoilwf, which is the "
                           "split exactly at one layer and is refused above "
                           "it"),
    "doro": Policy(STATIC_GRID, TARGET,
                   why="oroscale-scaled and spectrally fitted at the target "
                       "truncation, and under nglacier it carries the ice "
                       "sheet's orography too"),
    "groundsg": Policy(STATIC_GRID, TARGET, why="bare-ground orography"),
    "dglacsg": Policy(STATIC_GRID, TARGET, why="glacier orography"),
    "dglac": Policy(STATIC_GRID, TARGET, why="glacier mask"),
    "darea": Policy(STATIC_GRID, TARGET, remap=GEOMETRIC,
                    why="grid-cell area weights; a property of the target grid"),
    "dtcl": Policy(STATIC_GRID, TARGET, why="climatological surface temperature"),
    "dwcl": Policy(STATIC_GRID, TARGET, why="climatological soil wetness"),
    "xls": Policy(STATIC_GRID, TARGET, domain="ocean", why="the ice model's mask"),
    "yls": Policy(STATIC_GRID, TARGET, domain="ocean", why="the ocean's mask"),
    "xclsst": Policy(STATIC_GRID, TARGET, domain="ocean",
                     why="climatological sea surface temperature"),
    "xclicec": Policy(STATIC_GRID, TARGET, domain="ocean",
                      why="climatological ice cover"),
    "xcliced": Policy(STATIC_GRID, TARGET, domain="ocean",
                      why="climatological ice thickness"),
    "yclsst": Policy(STATIC_GRID, TARGET, domain="ocean",
                     why="the ocean's climatological sea surface temperature"),
})

# --- state the model rebuilds before it uses it again -------------------
# Remapped like anything else, and named so the report can say the target
# model will overwrite them. Taking the template's values instead put a
# one-orbit cold start's cloud field into a converted state: 0.39 in
# column cloud cover away from the donor's own. world-eyb.
POLICY.update({
    "dalb": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                  rebuilt_by_model=True,
                  why="surface albedo, rebuilt each step from snow, ice, vegetation and the background climatology"),
    "dsalb1": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                  rebuilt_by_model=True,
                  why="the same, below 0.75 um"),
    "dsalb2": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                  rebuilt_by_model=True,
                  why="the same, above 0.75 um"),
    "dz0": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                  rebuilt_by_model=True,
                  why="roughness, rebuilt from the climatology, snow and ice"),
    "dqsat": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                  rebuilt_by_model=True,
                  why="saturation humidity, a function of temperature and pressure. Identically zero in every restart on disk: written, carried, never populated"),
})


# --- prognostic gridpoint state --------------------------------------------
POLICY.update({
    "drhs": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                   why="surface wetness factor"),
    **{n: Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE,
                 why="SPAT-5 restart-required tile boundary temperature or humidity")
       for n in ("dlt_ts", "dlt_qs", "dot_ts", "dot_qs")},
    **{n: Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                 why="SPAT-5 restart-required tile wetness or albedo")
       for n in ("dlt_rhs", "dlt_alb", "dlt_sa1", "dlt_sa2",
                 "dot_rhs", "dot_alb", "dot_sa1", "dot_sa2")},
    **{n: Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                 why="SPAT-5 restart-required tile roughness")
       for n in ("dlt_z0", "dot_z0")},
    "dicec": Policy(PROGNOSTIC_GRID, REMAP, FRACTION, domain="ocean",
                    bounds=(0.0, 1.0), partner="diced", why="ice cover"),
    "diced": Policy(PROGNOSTIC_GRID, REMAP, THICKNESS, domain="ocean",
                    partner="dicec", conserve="sea_ice_volume",
                    measured_per="dicec",
                    why="ice thickness; the VOLUME is what survives a "
                        "coastline change, not the thickness"),
    "dwatc": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                    conserve="soil_water", bounds=(0.0, None),
                    why="soil water, metres of water"),
    "drunoff": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                      why="surface runoff rate"),
    "dwatcl": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="soil_water", bounds=(0.0, None),
                     why="LSHY-3. Soil water by layer, metres of water, for "
                         "the selectable layered land column. A RESERVOIR on "
                         "the same terms as dwatc, which remains the column "
                         "total: the two are one inventory carried at two "
                         "resolutions, and remapping either as an intensity "
                         "would not conserve it. Written on every step "
                         "whatever nlandwcol is, so a run can switch schemes "
                         "at a restart boundary without a cold start; a "
                         "restart written before the layers existed is "
                         "rebuilt from dwatc on the declared layer shape, "
                         "which is exact at one layer"),
    "dsoili": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="soil_water", bounds=(0.0, None),
                     why="LSHY-5. Soil ice by layer, metres of WATER "
                         "equivalent, so it is the same substance as dwatcl "
                         "and conserves against the same inventory: freeze "
                         "and thaw move mass between the two and neither is a "
                         "source. A RESERVOIR for that reason, and remapping "
                         "it as an intensity would break the pair. Zero under "
                         "the default nlandwphase = 0, where the land column "
                         "has no soil ice at all"),
    "ddrain": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="LSHY-3. Drainage out of the base of the land "
                         "column, metres per second. Identically zero under "
                         "the default impermeable lower boundary, and the "
                         "term the land water ledger's `drainage` crossing "
                         "reads once a lower boundary is selected"),
    "adrain": Policy(ACCUMULATOR, RESET, model_reset="zero",
                     why="WORLD-P9QQ. The output accumulator for ddrain, "
                         "declared in landmod beside the flux it accumulates "
                         "rather than beside aroff in plasimmod, and reset by "
                         "outmod's outreset like every other output "
                         "accumulator. An ACCUMULATOR and not a reservoir: it "
                         "is a rate summed over the output window and divided "
                         "by the counter, so a restart that lacks the record "
                         "loses a partial window and nothing more, which is "
                         "why landmod reads it under a lowered nexcheck"),
    "dust3": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE,
                    why="friction velocity cubed, the coupling quantity"),
    "dcc": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                  rebuilt_by_model=True,
                  why="cloud cover per level, and it is DIAGNOSTIC rather than "
                      "integrated: `rainmod.f90:1918` zeroes it and "
                      "`mkclouds` rebuilds the whole field from humidity and "
                      "temperature every timestep under NCLOUDS = 1, which is "
                      "what the runs set. Its restart value is used once, by "
                      "the radiation in `fluxstep`, before `rainstep` "
                      "overwrites it"),
    "dql": Policy(DERIVED_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                  rebuilt_by_model=True,
                  why="cloud liquid water. Zeroed on the line above dcc at "
                      "`rainmod.f90:1917` and rebuilt in the same routine, so "
                      "it is diagnostic on the same terms"),
    "dt": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE,
                 why="the surface level of the temperature array; the levels "
                     "above it are spectral and travel in st"),
    "dq": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                 why="surface humidity where humidity is spectral, and the "
                     "whole column where it is semi-Lagrangian"),
    "mmr": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                  why="aerosol mass mixing ratio, semi-Lagrangian humidity only"),
    "nrho": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                   why="aerosol number density, semi-Lagrangian humidity only"),
    "dtsl": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                   why="the land module's surface temperature"),
    "dtsm": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                   why="the land module's mean surface temperature"),
    "dqs": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                  bounds=(0.0, None), why="land surface humidity"),
    "driver": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="river_water", bounds=(0.0, None),
                     why="river runoff water in store, metres"),
    "duroff": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="zonal runoff velocity"),
    "dvroff": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="meridional runoff velocity"),
    "dsnowt": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="snow temperature"),
    "dsnowz": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="land_snow", bounds=(0.0, None),
                     why="snow depth in metres of water equivalent"),
    "dcansn": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="land_snow", bounds=(0.0, None),
                     why="canopy snow in metres of water equivalent, the "
                         "reservoir cansnowmod holds between the snowfall and "
                         "the ground pack. It is the SAME SUBSTANCE as dsnowz "
                         "and conserves with it, because what a remap takes "
                         "out of one it must not lose: the canopy is where "
                         "some of the cell's snow is standing, not a different "
                         "kind of snow"),
    "persistt": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                       bounds=(0.0, None),
                       why="seconds of continuous snow cover at or above "
                           "glacelim; the glacier persistence clock. An "
                           "elapsed TIME per cell, so it remaps as an "
                           "intensity and not as an inventory"),
    "dsoilt": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="soil temperature per layer"),
    "xts": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                  why="the ice model's surface temperature"),
    "xicec": Policy(PROGNOSTIC_GRID, REMAP, FRACTION, domain="ocean",
                    bounds=(0.0, 1.0), partner="xiced", why="ice cover"),
    "xiced": Policy(PROGNOSTIC_GRID, REMAP, THICKNESS, domain="ocean",
                    partner="xicec", conserve="ice_model_volume",
                    measured_per="xicec",
                    why="ice thickness; remapped as volume with its cover"),
    "xsnow": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="ocean",
                    conserve="ice_snow", bounds=(0.0, None),
                    measured_per="xicec",
                    why="snow on sea ice, water equivalent. A DEPTH ON THE ICE "
                        "COLUMN rather than over the cell: `icemod.f90:591` "
                        "adds it to xiced in the same length units after a "
                        "density conversion, so its inventory carries the ice "
                        "cover the way xiced's does"),
    "xicecc": Policy(PROGNOSTIC_GRID, REMAP, FRACTION, domain="ocean",
                     bounds=(0.0, 1.0),
                     why="prognostically computed ice cover"),
    "ysst": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                   conserve="ocean_heat", why="ocean temperature per layer"),
    "yiflux": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                     why="residual heat flux into ice"),
    "yfldo": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                    why="flux from the deep ocean"),
    "dts": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                  why="the sea module's surface temperature"),
})
# Instantaneous coupler exchange fields: one model call's worth of flux, held
# so the other component can read it. Intensive, and on the ocean side.
POLICY.update({
    n: Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
              why="coupler exchange field, one model call of flux")
    for n in ("cheat", "cpme", "croff", "ctaux", "ctauy", "cust3", "csnow")
})

# --- SIMBA, present only under coupled vegetation --------------------------
POLICY.update({
    "dcveg": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                    conserve="vegetation_carbon", bounds=(0.0, None),
                    why="vegetation carbon"),
    "dcsoil": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                     conserve="soil_carbon", bounds=(0.0, None),
                     why="soil carbon"),
})
POLICY.update(_acc(["agpp", "agppl", "agppw", "alitter", "anogrow", "anpp",
                    "aresh", "adcsoil", "adcveg", "adlai"]))


# The three accumulators `outmod.f90:outreset` does not simply zero. Every tool
# that writes a clean accumulator has to carry these, and the general rule --
# "an accumulator's clean value is zero" -- is wrong for all three.
# `check_policy_covers_source` holds them against the source, which is how
# `atsami` was found: reading the list by eye had missed it.
POLICY["tempmin"] = Policy(
    ACCUMULATOR, RESET, model_reset="sentinel", reset_value=1.0e3,
    why="a running MINIMUM over the window. `outmod.f90:2469` resets it to "
        "1.0e3, not to zero, and a zeroed one reports a minimum of 0 K for "
        "the rest of the run")
POLICY["atsami"] = Policy(
    ACCUMULATOR, RESET, model_reset="sentinel", reset_value=1.0e10,
    why="the running MINIMUM surface air temperature over the window. "
        "`outmod.f90:2463` resets it to 1.0e10 while its maximum partner "
        "atsama resets to zero, and a zeroed one reports a minimum of 0 K "
        "for the rest of the run")
POLICY["asndch"] = Policy(
    ACCUMULATOR, RESET, model_reset="none",
    why="net accumulated snow depth change. `outmod.f90:2459` has its reset "
        "commented out on purpose -- 'Let the net snow change keep "
        "accumulating' -- so it spans the whole run rather than the window")


# THE ECOLOGICAL STREAM'S ACCUMULATORS, EFOR-2. A second output stream with its
# own interval and its own counter, written for the biosphere rather than for a
# climate diagnostic; `outmod.f90:ecoaccu` fills them and `ecoreset` clears them.
#
# They are here rather than under `_acc` because the four extrema do not reset
# to zero and because a converted restart has to know that these are means over
# a PARTIAL interval whose length is `naccueco`, not over `naccuout`. Divide by
# the wrong counter and every field in the first block after a conversion is
# scaled by the ratio of the two.
POLICY.update(_acc([
    "aecotas", "aecots", "aecops", "aecohus", "aecowind", "aecoswd", "aecoswu",
    "aecoswn", "aecolwn", "aecolwu", "aecoczen", "aecopr", "aecoprsn",
    "aecoprc", "aecoevap"]))
POLICY.update({
    "naccueco": Policy(
        ACCUMULATOR, RESET, model_reset="zero",
        why="how many timesteps of the ecological interval the accumulators "
            "above already hold. Separate from naccuout because the two "
            "streams have different intervals, and a counter that outlives "
            "what it counts is the defect `accuvers` exists for"),
    "ecovers": Policy(
        OPAQUE, TARGET,
        why="marks a restart that carries the ecological stream's partial "
            "interval, on the same terms as accuvers. mpgetgp scatters an "
            "uninitialised buffer when a record is absent, so the reader "
            "needs a marker rather than a lowered nexcheck; a restart without "
            "it starts the interval clean"),
})
for _n, _v in (("aecotasmx", -1.0e3), ("aecotsmx", -1.0e3),
               ("aecotasmn", 1.0e3), ("aecotsmn", 1.0e3)):
    POLICY[_n] = Policy(
        ACCUMULATOR, RESET, model_reset="sentinel", reset_value=_v,
        why="a running extremum of the ecological interval. `ecoreset` gives "
            "it a sentinel and not zero, in BOTH directions: unlike the "
            "regular stream's tempmax, a maximum zeroed here would report 0 K "
            "for the rest of the run on any cell that never reaches it")


# ---------------------------------------------------------------------------
# WHICH RECORDS A RESTART IS ALLOWED NOT TO CARRY. world-t16e.
# ---------------------------------------------------------------------------
#
# `required` is the default and it is the answer for all but the rows below. A
# record the file does not carry stops the model with its name printed, and that
# is the right failure wherever resuming would mean integrating from a state
# that is not the state. It is affordable because of rule 7: no canonical
# climatology lineage is declared, so every run and every climatology is
# disposable and the cost of a refused restart is a re-commissioning.
#
# `optional` is only ever justified by what is LOST. An accumulator over an
# output window loses a partial window: the first record after the resume is a
# mean over fewer timesteps than it declares, and `outreset` makes even that go
# away by discarding the whole interval. A PROGNOSTIC loses state the cell
# integrates from, and no reset makes that go away -- which is why `persistt`
# is on this list as a decision and not as an omission.
_OPTIONAL = {
    # The ecological stream, EFOR-2. Nineteen accumulators and their counter
    # behind one marker; without it the stream's interval starts clean, which
    # is what a run that predates the stream would have done anyway.
    **{n: ("ecovers", "ecoreset: the interval starts clean and says so")
       for n in ("aecotas", "aecots", "aecops", "aecohus", "aecowind",
                 "aecoswd", "aecoswu", "aecoswn", "aecolwn", "aecolwu",
                 "aecoczen", "aecopr", "aecoprsn", "aecoprc", "aecoevap",
                 "aecotasmx", "aecotsmx", "aecotasmn", "aecotsmn")},
    "naccueco": ("ecovers", "0, with the accumulators it counts"),
    # WORLD-3QFZ's band split, behind the accumulator-set marker.
    "afdsw1": ("accuvers", "outreset: the whole output interval is discarded"),
    "afdsw2": ("accuvers", "outreset: the whole output interval is discarded"),
    # world-5qy. Not behind accuvers, because accuvers promises a COMPLETE set
    # and these two are written only where nenergy and nener3d allocate them.
    # `read_atmos_restart` asks for them with has_restart_array instead, and a
    # run that wants them and does not find them discards the interval.
    "adenergy": (None, "outreset: the whole output interval is discarded"),
    "adener3d": (None, "outreset: the whole output interval is discarded"),
    # WORLD-VJBZ's per-cell capacity split, through mpgetgp_found. A restart
    # written before surface code 2290 existed has no per-layer split, and
    # landwfrac falls back to the namelist shape -- which is the split EXACTLY
    # at one water layer, and is refused above it wherever a land cell's dwmax
    # came from a soil, so the fallback cannot quietly stand in for a field it
    # is not equal to.
    "dsoilwfc": (None, "landwfrac falls back to the namelist dsoilwf, which is "
                       "the split exactly at one water layer and is refused "
                       "above it wherever dwmax came from a soil"),
    # world-onw8's four, through mpgetgp_found. The layered store and the
    # drainage rebuild from the single-layer fields a pre-LSHY-3 restart has.
    "dwatcl": (None, "landini rebuilds the layered store from dwatc"),
    # BIO-32's canopy snow store, through the same mpgetgp_found path. A
    # restart written before the store existed carries no canopy, and an empty
    # canopy is exactly the state that restart's world was in.
    "dcansn": (None, "zero: a restart written before the canopy snow store "
                     "existed held no canopy snow"),
    "dsoili": (None, "landini rebuilds the layered store from dwatc"),
    "ddrain": (None, "zero: no drainage is in flight at the resume"),
    "adrain": (None, "zero: the drainage output window starts clean"),
    # The six spectral accumulators and the five accumulated orbital scalars,
    # read inside the accuvers branch and under a lowered nexcheck. CLIM-2 and
    # CLIM-5 added them and a restart written before that has no complete set.
    **{n: ("accuvers", "outreset: the whole output interval is discarded")
       for n in ("aasosp", "aaspsp", "aastsp", "aasqsp", "aasdsp", "aaszsp",
                 "aorbnu", "alambm", "azdecl", "ardist", "arasc")},
    # The stamped geometry, world-4yz. Read under a lowered nexcheck into
    # locals pre-set to -1, and the mismatch test is written so that an
    # UNSTAMPED restart passes: an unstamped restart is not a mismatched one.
    **{n: (None, "-1, and the geometry check treats an unstamped restart as "
                 "unstamped rather than as wrong")
       for n in ("nlat", "nlon", "nlev", "nrsp")},
    # The two version markers themselves, and the fixer's correction, all read
    # under a lowered nexcheck into a scalar the caller pre-set.
    "accuvers": (None, "-1.0, which is what makes the marker work"),
    "ecovers": (None, "-1.0, which is what makes the marker work"),
    "denergyfix": (None, "0.0, the declared value: that segment starts "
                         "uncorrected and reconverges in one window"),
}
#
# THE DECISION THAT IS NOT ON THE LIST. `persistt` is glaciermod's snowpack
# persistence clock and it stays REQUIRED. It is a prognostic, not an
# accumulator: a cell whose snowpack has stood above `glacelim` for most of
# `glacpersec` is one step from becoming glacier, and a clock read as zero puts
# that cell back to the start of the threshold. The loss is bounded by the
# threshold and is silent, which is the worst pair -- an accumulator's loss is
# bounded by one output window and shows up in one record. So a restart written
# before the clock existed is refused, loudly, with the record's name printed,
# and the answer to it is to re-commission rather than to resume.
for _n, (_marker, _cold) in _OPTIONAL.items():
    POLICY[_n] = _replace(POLICY[_n], resume="optional",
                          resume_marker=_marker, resume_cold=_cold)
del _n, _marker, _cold


# `nexcheck = 0` and `nexcheck = 1` bracket a region where an absent record
# leaves the caller's variable alone instead of stopping the model.
_NEXCHECK = re.compile(r"^\s*nexcheck\s*=\s*(?P<v>[01])\s*$", re.IGNORECASE)
_READ = re.compile(
    r"""^[^!\n]*?
        call\s+(?P<fn>get_restart_integer|get_restart_real|get_restart_seed
                    |get_restart_array|has_restart_array
                    |mpgetgp_found|mpgetgp|mpgetsp)
        \s*\(\s*'(?P<name>[A-Za-z0-9_]+)'""",
    re.IGNORECASE | re.VERBOSE)


def read_sites_from_source(src_dir: Path) -> dict:
    """{record: how the model READS it}, from the call sites themselves.

    Three facts per name, and each of them is what decides whether an absent
    record stops the model: which readers touch it, whether any read is inside a
    lowered-`nexcheck` region, and whether anything ever ASKS for it with
    `has_restart_array`. The write side is `inventory_from_source`; this is the
    other half, and `check_policy_covers_source` holds the decision table
    against it.
    """
    src_dir = Path(src_dir)
    out: dict = {}
    for module in compiled_modules(src_dir.parent):
        path = src_dir / module
        if not path.is_file():
            continue
        checked = True
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            bare = line.split("!")[0]
            m = _NEXCHECK.match(bare)
            if m is not None:
                checked = m["v"] == "1"
                continue
            m = _READ.match(line)
            if m is None:
                continue
            fn, name = m["fn"].lower(), m["name"]
            site = out.setdefault(name, {"readers": set(), "unchecked": False,
                                         "asked": False})
            if fn == "has_restart_array":
                site["asked"] = True
                continue
            site["readers"].add(fn)
            if not checked and fn.startswith("get_restart_"):
                site["unchecked"] = True
    return out


def check_resume_policy(src_dir: Path) -> list[str]:
    """The resume decision per record, held against how the model reads it.

    A CHECK WITH A RIGHT ANSWER IN BOTH DIRECTIONS. A record recorded as
    `optional` that nothing in the model makes optional is a decision the code
    does not implement, and a record recorded as `required` that the model
    already resumes without is a decision nobody took. Either way the table and
    the source disagree about what happens when a restart is short a record,
    and that question is only ever answered at the start of the next segment.
    """
    inventory = inventory_from_source(src_dir)
    reads = read_sites_from_source(src_dir)
    problems = []

    def handled(name: str) -> bool:
        site = reads.get(name)
        if site is None:
            return False
        if "mpgetgp_found" in site["readers"] or site["unchecked"] or site["asked"]:
            return True
        marker = POLICY[name].resume_marker if name in POLICY else None
        return bool(marker) and marker in POLICY and \
            POLICY[marker].resume == "optional" and handled(marker)

    for name in sorted(set(POLICY) & set(inventory)):
        pol = POLICY[name]
        if pol.resume not in ("required", "optional"):
            problems.append(
                f"'{name}' records resume={pol.resume!r}, which is neither "
                "'required' nor 'optional'")
            continue
        if pol.resume_marker and pol.resume_marker not in POLICY:
            problems.append(
                f"'{name}' is read behind the marker '{pol.resume_marker}', "
                "which is not a record this schema knows")
            continue
        if name not in reads:
            if pol.resume == "optional":
                problems.append(
                    f"'{name}' is recorded as optional on resume and no call "
                    "site reads it, so there is nothing for the record's "
                    "absence to be optional about")
            continue
        if pol.resume == "optional" and not handled(name):
            problems.append(
                f"'{name}' is recorded as optional on resume and the model "
                f"reads it with {sorted(reads[name]['readers'])} at nexcheck 1 "
                "behind no marker this schema knows, so an absent record stops "
                "the model. Read it through mpgetgp_found, ask for it with "
                "has_restart_array, put it behind a marker, or record it as "
                "required.")
        if pol.resume == "required" and handled(name):
            problems.append(
                f"'{name}' is recorded as required on resume and the model "
                "already resumes without it. Say what it starts from and "
                "record it as optional, or take the handling back out.")
    return problems


def check_record_limit(src_dir: Path, emittable: int) -> list[str]:
    """`restartmod.nresdim` is above the number of names the model can emit.

    A CHECK WITH A RIGHT ANSWER, and one that has already been wrong.
    `restart_ini` walks the restart file naming each record and stops the model
    the moment it reaches `nresdim`. The limit stood at 200 while the call sites
    could emit 215 names, so a configuration that wrote enough of the optional
    records produced a restart it could not read back -- and the failure lands
    at the START of the following segment, after the segment that wrote it has
    finished and reported success.

    The bound is the emittable count and not the count in any one file, because
    which optional records a configuration writes is a namelist question and
    this constant is compiled in.
    """
    text = (Path(src_dir) / "restartmod.f90").read_text(encoding="utf-8")
    m = re.search(r"nresdim\s*=\s*(\d+)", text)
    if m is None:
        return ["restartmod.f90 declares no nresdim; restart_ini's record "
                "limit cannot be checked"]
    limit = int(m.group(1))
    if limit <= emittable:
        return [f"restartmod.f90 sets nresdim = {limit} and the model can emit "
                f"{emittable} distinct restart records. restart_ini stops the "
                f"model at the limit, so a run that writes them all produces a "
                f"restart it cannot read back."]
    return []


def check_policy_covers_source(src_dir: Path) -> list[str]:
    """Every name the model can write has a policy, and no policy is orphaned.

    The failure this catches is a real one and it has a right answer: a record
    added to the model with no policy would otherwise reach the converter as an
    unknown name at the moment somebody needed the conversion.
    """
    inventory = inventory_from_source(src_dir)
    problems = []
    problems += check_record_limit(src_dir, len(inventory))
    problems += check_resume_policy(src_dir)
    for name in sorted(set(inventory) - set(POLICY)):
        rec = inventory[name]
        problems.append(
            f"{rec.module}:{rec.line} writes '{name}' and restart_schema.POLICY "
            "has no entry for it")
    for name in sorted(set(POLICY) - set(inventory)):
        problems.append(
            f"POLICY names '{name}' and no call site in any source "
            f"{CMAKE} compiles writes it")

    # The reset column, held against the source that owns it. This is what
    # keeps "an accumulator's clean value is zero" from becoming a rule the
    # model has quietly stopped obeying.
    resets = model_resets_from_source(src_dir)
    for name in sorted(set(POLICY) & set(inventory)):
        pol = POLICY[name]
        if pol.semantic != ACCUMULATOR:
            continue
        kind, value = derived_model_reset(name, inventory, resets)
        if kind != pol.model_reset:
            problems.append(
                f"'{name}' is recorded as model_reset={pol.model_reset!r} and "
                f"the source resets {inventory[name].variable} as {kind!r}")
        elif kind == "sentinel" and value != pol.reset_value:
            problems.append(
                f"'{name}' is recorded as resetting to {pol.reset_value} and "
                f"the source resets it to {value}")
    return problems


# ---------------------------------------------------------------------------
# Geometry: the symbols a shape is written in, resolved for one file.
# ---------------------------------------------------------------------------

# `integer, parameter :: NLSOILWX = 8` and the fixed-form `parameter(NLSOIL=5)`.
# Only a plain integer literal is taken: `NLAT = NLAT_ATM` and
# `NGEN = max(NUGP, NESP*NLEV)` are derived from symbols this module already
# owns, and a parser that guessed at them would answer with the wrong build's
# geometry.
_PARAMETER = re.compile(
    r"^\s*(?:integer\s*(?:\([^)]*\))?\s*,\s*parameter\s*::\s*"
    r"(?P<n1>[a-z][a-z0-9_]*)\s*=\s*(?P<v1>[-+]?\d+)"
    r"|parameter\s*\(\s*(?P<n2>[a-z][a-z0-9_]*)\s*=\s*(?P<v2>[-+]?\d+)\s*\))\s*$",
    re.IGNORECASE)


@functools.lru_cache(maxsize=8)
def _parameter_dimensions(src_dir: str) -> tuple:
    return tuple(sorted(_scan_parameters(Path(src_dir)).items()))


def _scan_parameters(src_dir: Path) -> dict:
    out: dict = {}
    for module in compiled_modules(src_dir.parent):
        path = src_dir / module
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.lstrip().startswith("!"):
                continue
            m = _PARAMETER.match(line.split("!")[0])
            if m is None:
                continue
            name = (m["n1"] or m["n2"]).upper()
            value = int(m["v1"] if m["v1"] is not None else m["v2"])
            out.setdefault(name, value)
    return out


def parameter_dimensions_from_source(src_dir=MODEL_SRC) -> dict:
    """{SYMBOL: value} for every compile-time integer parameter in the model.

    A restart record's shape is written in whatever symbols its call site uses,
    and not all of them are derivable from NLAT. `dwatcl(NHOR,NLSOILWX)` is
    shaped by a parameter declared in `landcolumn.f90`, so the length of that
    record is a property of the SOURCE the writing executable was built from
    and of nothing in the file. Reading it here is what keeps the schema from
    carrying a second copy of a number the model already states, which is the
    same reason `inventory_from_source` reads the call sites rather than
    listing them.

    Consulted only for symbols `Geometry.resolve` does not already own, so a
    source `NLAT` local to a helper module can never shadow the file's own.
    """
    return dict(_parameter_dimensions(str(Path(src_dir))))


@dataclass(frozen=True)
class Geometry:
    """The dimensions a restart's record lengths are built out of.

    `plasimmod.f90` derives most of them from NLAT alone, and this reproduces
    that derivation rather than reading lengths back out of the file.

    NESP and NSEEDLEN are the exceptions and they are carried rather than
    derived. NESP is NRSP rounded UP to a multiple of the process count, so a
    restart's spectral-accumulator length depends on how many threads the run
    held -- and it does not run backwards: 506 rounds to 512 under both eight
    and sixteen, so a file cannot say which it was. NSEEDLEN belongs to the
    compiler. Both are properties of the executable that wrote the file, which
    is why the target template is where the target's values come from.
    """

    nlat: int
    nlev: int
    nlsoil: int
    nlev_oce: int
    nesp: int
    nseedlen: int
    # A COMPILE-TIME parameter, not a header field: the model writes `dwatcl`
    # and `dsoili` at NLSOILWX columns whatever `nlsoilw` the namelist asks
    # for, zeroing the layers above it. It reaches a geometry from the source
    # the writing executable was built from -- `parameter_dimensions_from_
    # source` -- and is zero only where a caller built a Geometry by hand and
    # did not supply it, which `resolve` refuses rather than silently sizing a
    # record at nothing.
    nlsoilwx: int = 0

    @property
    def nlon(self) -> int:
        return 2 * self.nlat

    @property
    def ntru(self) -> int:
        return (self.nlon - 1) // 3

    @property
    def nrsp(self) -> int:
        return (self.ntru + 1) * (self.ntru + 2)

    @property
    def ncsp(self) -> int:
        return self.nrsp // 2

    @property
    def nugp(self) -> int:
        return self.nlon * self.nlat

    @property
    def nlep(self) -> int:
        return self.nlev + 1

    def __post_init__(self):
        if self.nesp < self.nrsp:
            raise ValueError(
                f"NESP {self.nesp} is below NRSP {self.nrsp}; the spectral "
                "accumulators cannot hold fewer modes than the state does")

    def resolve(self, token) -> int:
        if isinstance(token, int):
            return token
        if "*" in token:
            n = 1
            for part in token.split("*"):
                part = part.strip()
                n *= int(part) if part.isdigit() else self.resolve(part)
            return n
        try:
            value = {"NUGP": self.nugp, "NRSP": self.nrsp, "NESP": self.nesp,
                     "NLEV": self.nlev, "NLEP": self.nlep, "NLSOIL": self.nlsoil,
                     "NLEV_OCE": self.nlev_oce, "NSEEDLEN": self.nseedlen,
                     "NLAT": self.nlat, "NLON": self.nlon,
                     "NLSOILWX": self.nlsoilwx}[token]
        except KeyError:
            raise KeyError(
                f"'{token}' is a dimension no geometry here knows. It came out "
                "of a call site, so the model has grown a symbol this module "
                "has not been taught. Compile-time parameters are read by "
                "parameter_dimensions_from_source(); add the symbol to "
                "Geometry if the model derives it from NLAT instead."
            ) from None
        # Only the compile-time parameters are guarded here. NLSOIL,
        # NLEV_OCE and NSEEDLEN are HEADER fields and legitimately zero on a
        # file that carries no soil, no ocean or no seed; a shape using one of
        # those on such a file gives a candidate count of zero, which matches
        # no payload and is refused where the count is compared. NLSOILWX has
        # no header to be absent from, so a zero there is a caller that did not
        # supply it and would silently size every land-water record at nothing.
        if token == "NLSOILWX" and value <= 0:
            raise ValueError(
                "NLSOILWX is a compile-time parameter of the model source and "
                "this Geometry was built without one, so `dwatcl` and `dsoili` "
                "would be sized at nothing. Build the Geometry through "
                "describe(), or pass nlsoilwx from "
                "parameter_dimensions_from_source().")
        return value

    def elements(self, shape) -> int:
        n = 1
        for token in shape:
            n *= self.resolve(token)
        return n

    @property
    def label(self) -> str:
        return f"T{self.ntru}"


def element_count(name: str, geometry: Geometry, inventory) -> int:
    """How many elements `name` holds in a file with this geometry.

    Where the source can write a name with more than one shape, every
    alternative that a geometry can express is a candidate and the caller
    settles it against the payload it actually has.
    """
    return geometry.elements(inventory[name].shapes[0])


def candidate_counts(name: str, geometry: Geometry, inventory) -> list[int]:
    return [geometry.elements(shape) for shape in inventory[name].shapes]


def check_every_shape_resolves(src_dir=MODEL_SRC) -> list[str]:
    """Every dimension symbol a call site writes a shape in is one we can size.

    `check_policy_covers_source` compares NAMES and stays green while a shape
    names a symbol no geometry knows, which is what `dwatcl(NHOR,NLSOILWX)`
    did: the policy covered `dwatcl`, and every attempt to predict its length
    raised instead. This asks the resolver, and it is source-only -- no restart
    file, no run output -- so it answers on a fresh checkout and on the day a
    new dimension lands rather than the next time somebody converts a state.

    The probe geometry is deliberately not a rung: nothing here depends on the
    SIZES, only on whether each symbol is one `Geometry.resolve` handles.
    """
    params = parameter_dimensions_from_source(src_dir)
    probe = Geometry(nlat=2, nlev=1, nlsoil=1, nlev_oce=1, nesp=6, nseedlen=1,
                     nlsoilwx=params.get("NLSOILWX", 0))
    problems = []
    for name, rec in sorted(inventory_from_source(src_dir).items()):
        for shape in rec.shapes:
            for token in shape:
                if isinstance(token, int):
                    continue
                try:
                    probe.resolve(token)
                except (KeyError, ValueError) as exc:
                    problems.append(
                        f"'{name}' is written at {rec.module}:{rec.line} with "
                        f"shape {shape} and {exc}")
    return sorted(set(problems))


# ---------------------------------------------------------------------------
# What the model does to an accumulator at an interval boundary
# ---------------------------------------------------------------------------

# An assignment of a bare numeric literal, outside a declaration and outside a
# comment. `aast(:,j) = 0.`, `tempmin(:) = 1.0e3` and `naccuout=0` all match.
_RESET = re.compile(
    r"^\s*(?P<var>[a-z][a-z0-9_]*)\s*(\([^)]*\))?\s*="
    r"\s*(?P<value>[-+]?(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?)\s*$",
    re.IGNORECASE)


def model_resets_from_source(src_dir: Path) -> dict:
    """{variable: the constant the model assigns it}, over the compiled modules.

    An accumulator's clean value is whatever `outreset` and its per-module
    equivalents put there, and that is a fact about the source rather than a
    convention: `tempmin` goes to 1.0e3 because it is a running minimum, and a
    tool that assumed zero would report a minimum of 0 K for the rest of a run.
    Deriving it is what makes the `model_reset` column falsifiable instead of
    asserted.
    """
    src_dir = Path(src_dir)
    resets: dict = {}
    for module in compiled_modules(src_dir.parent):
        for line in (src_dir / module).read_text(encoding="utf-8",
                                                 errors="ignore").splitlines():
            if "::" in line or line.lstrip().startswith("!"):
                continue
            m = _RESET.match(line.split("!")[0])
            if m is None:
                continue
            resets.setdefault(m["var"].lower(), set()).add(float(m["value"]))
    return resets


def derived_model_reset(name: str, inventory: dict, resets: dict):
    """('zero' | 'sentinel' | 'none', value) for one record, from the source."""
    values = resets.get(inventory[name].variable.lower())
    if not values:
        return "none", None
    if values == {0.0}:
        return "zero", None
    nonzero = sorted(values - {0.0})
    return "sentinel", nonzero[-1]


class ConversionError(Exception):
    """A restart cannot be worked on, and working on it partly would be worse."""


def _int(record) -> int:
    if record.nbytes != 4:
        raise ConversionError(
            f"'{record.name}' should be a four-byte integer and is "
            f"{record.nbytes} bytes")
    return struct.unpack("<i", record.payload)[0]


def infer_real_bytes(by_name: dict, nlat: int, nlev: int, nrsp: int) -> int:
    """Four or eight, agreed by every record whose element count is known.

    The element counts come from the integer headers, which are always four
    bytes, so this does not assume the answer it is looking for. Every check
    must agree: a length divisible by eight is not evidence on its own, since
    an integer array and a four-byte real array of the same count are the same
    size, and a disagreement means the file is not what the headers say.
    """
    invariants = {"sp": nrsp, "sz": nrsp * nlev, "dls": nlat * 2 * nlat}
    votes = {}
    for name, count in invariants.items():
        if name not in by_name:
            continue
        nbytes = by_name[name].nbytes
        if nbytes % count:
            raise ConversionError(
                f"'{name}' is {nbytes} bytes and holds {count} elements, "
                "which is not a whole number of bytes each")
        votes[name] = nbytes // count
    if not votes:
        raise ConversionError(
            "none of sp, sz or dls is present, so the real width cannot be "
            "established from anything the headers already pin down")
    widths = set(votes.values())
    if len(widths) > 1:
        detail = ", ".join(f"{n} says {w}" for n, w in sorted(votes.items()))
        raise ConversionError(
            f"the records disagree about the real width ({detail}). The file "
            "does not match its own headers.")
    width = widths.pop()
    if width not in (4, 8):
        raise ConversionError(
            f"a real width of {width} bytes; this model is built at four or "
            "eight and nothing here can guess at another")
    return width


def describe(records, src_dir=MODEL_SRC) -> tuple:
    """(Geometry, real width) for a parsed restart, from its own headers.

    `src_dir` supplies the compile-time parameters the headers cannot: a
    restart states NLAT, NLEV and NLSOIL, but nothing in it says how many
    layers `dwatcl` was written at. That is NLSOILWX in `landcolumn.f90` and it
    comes from the source, so a caller reading a file written by another tree
    passes that tree's source rather than trusting this one's.
    """
    by_name = {rec.name: rec for rec in records}
    for needed in ("nlat", "nlon", "nlev", "nrsp"):
        if needed not in by_name:
            raise ConversionError(
                f"no '{needed}' record. A restart this model wrote always "
                "opens with its geometry.")
    nlat, nlon = _int(by_name["nlat"]), _int(by_name["nlon"])
    nlev, nrsp = _int(by_name["nlev"]), _int(by_name["nrsp"])
    if nlon != 2 * nlat:
        raise ConversionError(
            f"NLON {nlon} is not twice NLAT {nlat}; no grid contract here "
            "covers that")
    real_bytes = infer_real_bytes(by_name, nlat, nlev, nrsp)
    geometry = Geometry(
        nlat=nlat, nlev=nlev,
        nlsoil=_int(by_name["nlsoil"]) if "nlsoil" in by_name else 0,
        nlev_oce=_int(by_name["nlev_oce"]) if "nlev_oce" in by_name else 0,
        nesp=(by_name["aasosp"].nbytes // real_bytes
              if "aasosp" in by_name else nrsp),
        nseedlen=by_name["seed"].nbytes // 4 if "seed" in by_name else 0,
        nlsoilwx=parameter_dimensions_from_source(src_dir).get("NLSOILWX", 0))
    if geometry.nrsp != nrsp:
        raise ConversionError(
            f"the file says NRSP is {nrsp} and a T{geometry.ntru} grid gives "
            f"{geometry.nrsp}; the header and the grid disagree")
    return geometry, real_bytes


# A declaration with an initial value: `real :: atsami(NHOR) = 1.E10`.
_DECLARED = re.compile(
    r"::\s*(?P<var>[a-z][a-z0-9_]*)\s*(\([^)]*\))?\s*="
    r"\s*(?P<value>[-+]?(\d+\.?\d*|\.\d+)([eEdDqQ][-+]?\d+)?)\s*(!.*)?$",
    re.IGNORECASE)


def _fortran_real(token: str) -> float:
    return float(re.sub("[dDqQ]", "e", token))


def declared_initials_from_source(src_dir: Path) -> dict:
    """{variable: its declared initial value}, over the compiled modules."""
    src_dir = Path(src_dir)
    declared: dict = {}
    for module in compiled_modules(src_dir.parent):
        for line in (src_dir / module).read_text(encoding="utf-8",
                                                 errors="ignore").splitlines():
            if line.lstrip().startswith("!"):
                continue
            m = _DECLARED.search(line.split("!")[0])
            if m is not None:
                declared.setdefault(m["var"].lower(),
                                    _fortran_real(m["value"]))
    return declared


def check_first_window_matches_the_rest(src_dir: Path) -> list[str]:
    """An accumulator starts a cold run at the value its reset would give it.

    A right answer rather than a comparison. `nhcstp` starts at 1 and the model
    writes output at `mod(nhcstp,nafter) == 0`, so a cold-started run's FIRST
    output window accumulates from the declared initial value with no reset
    before it -- every later window starts from `outreset`. If the two differ,
    the first record of every cold start is a different quantity from the rest
    of the run, and nothing says so.

    `atsami` is why this exists: the running minimum surface air temperature
    was declared at 0.0 and reset to 1.0e10, so `AMIN1(0.0, anything)` held it
    at zero and the first output record of every cold start reported a minimum
    of 0 K. Its maximum partner `atsama` is declared and reset at 0.0 alike,
    which is what the pair should look like.
    """
    inventory = inventory_from_source(src_dir)
    resets = model_resets_from_source(src_dir)
    declared = declared_initials_from_source(src_dir)
    problems = []
    for name in sorted(POLICY):
        pol = POLICY[name]
        if pol.semantic != ACCUMULATOR or name not in inventory:
            continue
        var = inventory[name].variable.lower()
        kind, value = derived_model_reset(name, inventory, resets)
        if kind == "none" or var not in declared:
            continue
        want = 0.0 if kind == "zero" else value
        if declared[var] != want:
            problems.append(
                f"'{name}' ({var}) is declared at {declared[var]:g} and reset "
                f"to {want:g}, so a cold start's first output window is a "
                "different quantity from every window after it")
    return problems
