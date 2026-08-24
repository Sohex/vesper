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
and the check goes red rather than the converter guessing.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# The mechanical half: what the model source says it writes.
# ---------------------------------------------------------------------------

# The modules `plasim/bld/make_plasim` links into `plasim.x` that write restart
# records. `plasim_dummy.f90` and `icemod_template.f90` are NOT in that list --
# the first is a separate program and the second is the template `icemod.f90`
# was generated from -- so their names are not records this executable can
# emit, and treating them as optional would weaken the unknown-record error
# into nothing.
WRITING_MODULES = (
    "plasim.f90", "landmod.f90", "glaciermod.f90", "icemod.f90",
    "oceanmod.f90", "seamod.f90", "radmod.f90", "simba.f90",
)

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
    """A dimension as an int where it is a literal, else its symbol."""
    token = token.strip()
    return int(token) if token.isdigit() else token.upper()


@dataclass(frozen=True)
class SourceRecord:
    """One `put_*` call site, as the source states it."""

    name: str
    module: str
    line: int
    writer: str
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
    for module in WRITING_MODULES:
        path = src_dir / module
        if not path.is_file():
            raise FileNotFoundError(f"{path} is missing; the model source moved")
        for lineno, line in enumerate(path.read_text(encoding="utf-8",
                                                     errors="ignore").splitlines(), 1):
            m = _CALL.match(line)
            if m is None:
                continue
            fn, name, args = m["fn"].lower(), m["name"], _split_args(m["rest"])
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
                        writer=prior.writer, shapes=prior.shapes + (shape,),
                        guard=prior.guard)
                continue
            found[name] = SourceRecord(
                name=name, module=module, line=lineno, writer=fn,
                shapes=(shape,), guard=g["cond"].strip() if g else None)
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
    # True where SIMBA takes ownership of the field under coupled vegetation.
    vegetation_owned: bool = False
    why: str = ""


def _acc(names, *, klev=None):
    return {n: Policy(ACCUMULATOR, RESET,
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
    "zsolars": Policy(OPAQUE, REQUIRE_EQUAL,
                      why="the two-band split of the stellar constant is "
                          "derived from configuration, not from the grid"),
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

# --- derived state: the target model rebuilds it ---------------------------
POLICY.update({
    "dalb": Policy(DERIVED_GRID, RECOMPUTE,
                   why="albedo is a function of snow, ice, vegetation and the "
                       "background climatology; reproducing the formula here "
                       "would be a second implementation of the physics"),
    "dsalb1": Policy(DERIVED_GRID, RECOMPUTE, why="the same, below 0.75 um"),
    "dsalb2": Policy(DERIVED_GRID, RECOMPUTE, why="the same, above 0.75 um"),
    "dz0": Policy(DERIVED_GRID, RECOMPUTE,
                  why="roughness follows from the climatology, snow and ice"),
    "dqsat": Policy(DERIVED_GRID, RECOMPUTE,
                    why="saturation humidity is a function of temperature and "
                        "pressure at the target's own levels"),
})

# --- prognostic gridpoint state --------------------------------------------
POLICY.update({
    "drhs": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                   why="surface wetness factor"),
    "dicec": Policy(PROGNOSTIC_GRID, REMAP, FRACTION, domain="ocean",
                    bounds=(0.0, 1.0), partner="diced", why="ice cover"),
    "diced": Policy(PROGNOSTIC_GRID, REMAP, THICKNESS, domain="ocean",
                    partner="dicec", conserve="sea_ice_volume",
                    why="ice thickness; the VOLUME is what survives a "
                        "coastline change, not the thickness"),
    "dwatc": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="land",
                    conserve="soil_water", bounds=(0.0, None),
                    why="soil water, metres of water"),
    "drunoff": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                      why="surface runoff rate"),
    "dust3": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE,
                    why="friction velocity cubed, the coupling quantity"),
    "dcc": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, 1.0),
                  why="cloud cover per level"),
    "dql": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, bounds=(0.0, None),
                  why="cloud liquid water, a mixing ratio"),
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
    "dsoilt": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="land",
                     why="soil temperature per layer"),
    "xts": Policy(PROGNOSTIC_GRID, REMAP, INTENSIVE, domain="ocean",
                  why="the ice model's surface temperature"),
    "xicec": Policy(PROGNOSTIC_GRID, REMAP, FRACTION, domain="ocean",
                    bounds=(0.0, 1.0), partner="xiced", why="ice cover"),
    "xiced": Policy(PROGNOSTIC_GRID, REMAP, THICKNESS, domain="ocean",
                    partner="xicec", conserve="ice_model_volume",
                    why="ice thickness; remapped as volume with its cover"),
    "xsnow": Policy(PROGNOSTIC_GRID, REMAP, RESERVOIR, domain="ocean",
                    conserve="ice_snow", bounds=(0.0, None),
                    why="snow on sea ice, water equivalent"),
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


def check_policy_covers_source(src_dir: Path) -> list[str]:
    """Every name the model can write has a policy, and no policy is orphaned.

    The failure this catches is a real one and it has a right answer: a record
    added to the model with no policy would otherwise reach the converter as an
    unknown name at the moment somebody needed the conversion.
    """
    inventory = inventory_from_source(src_dir)
    problems = []
    for name in sorted(set(inventory) - set(POLICY)):
        rec = inventory[name]
        problems.append(
            f"{rec.module}:{rec.line} writes '{name}' and restart_schema.POLICY "
            "has no entry for it")
    for name in sorted(set(POLICY) - set(inventory)):
        problems.append(
            f"POLICY names '{name}' and no call site in "
            f"{', '.join(WRITING_MODULES)} writes it")
    return problems


# ---------------------------------------------------------------------------
# Geometry: the symbols a shape is written in, resolved for one file.
# ---------------------------------------------------------------------------

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
        try:
            return {"NUGP": self.nugp, "NRSP": self.nrsp, "NESP": self.nesp,
                    "NLEV": self.nlev, "NLEP": self.nlep, "NLSOIL": self.nlsoil,
                    "NLEV_OCE": self.nlev_oce, "NSEEDLEN": self.nseedlen,
                    "NLAT": self.nlat, "NLON": self.nlon}[token]
        except KeyError:
            raise KeyError(
                f"'{token}' is a dimension no geometry here knows. It came out "
                "of a call site, so the model has grown a symbol this module "
                "has not been taught.") from None

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
