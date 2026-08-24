"""The ExoPlaSim resolution ladder: the one place it is written down.

Worldbuilding frame: the ladder is the set of Gaussian grids the Vesper
climate model can be solved on. Nothing here is about the real world.

WHY THIS IS ITS OWN MODULE and not part of `lib/gridding.py`: the build script
and two shell probes need it, and neither should have to import numpy and the
World Orogen export reader to find out how many latitudes T85 has.
`lib/gridding.py` re-exports everything below, so there is one definition and
two doors rather than two definitions.

WHAT IT REPLACED. `config/planet.yaml` stated `resolution`, `latitudes` and
`longitudes` independently, so a valid truncation could be paired with another
grid's dimensions and nothing said so; `build_model.py`, the SHTns variant
sweep, `verify_shtns_equivalence.sh` and `run_shtns_probe.sh` each carried
their own copy of the mapping, and two of the four were missing rungs the
others had. SPAT-2.
"""
from __future__ import annotations


# THE TABLE IS SELF-CHECKING rather than merely agreed: PlaSim derives the
# truncation from the grid as `(NLON-1)/3` with `NLON = 2*NLAT`, so a rung's
# name and its latitude count are not two facts but one, and `_check_rungs()`
# below rejects any row where they disagree. That is a right answer, not a
# convention.
RUNGS = {"T21": 32, "T31": 48, "T42": 64, "T63": 96,
         "T85": 128, "T106": 160, "T127": 192, "T170": 256}


def _check_rungs() -> None:
    for rung, nlat in RUNGS.items():
        if (2 * nlat - 1) // 3 != int(rung[1:]):
            raise RuntimeError(
                f"the rung table is inconsistent: {rung} claims {nlat} "
                f"latitudes, and a {2 * nlat}-column Gaussian grid truncates "
                f"at T{(2 * nlat - 1) // 3}")


_check_rungs()


def geometry(rung: str) -> tuple[int, int, int]:
    """(latitudes, longitudes, truncation) for a ladder rung.

    The only accepted spellings are the table's. A rung outside it must be an
    ERROR rather than a default: `-r 170` once built T21 and said nothing.
    """
    key = str(rung).upper()
    if key not in RUNGS:
        raise RuntimeError(
            f"{rung!r} is not a ladder rung. The ladder is "
            f"{', '.join(RUNGS)}; anything else needs a row here and a grid "
            "in the export before it can be run.")
    nlat = RUNGS[key]
    return nlat, 2 * nlat, int(key[1:])


def rung_of_latitudes(nlat: int) -> str:
    """The rung with this many latitudes, or an error naming what is close."""
    for rung, n in RUNGS.items():
        if n == int(nlat):
            return rung
    raise RuntimeError(
        f"no ladder rung has {nlat} latitudes; the ladder is "
        + ", ".join(f"{r} ({n})" for r, n in RUNGS.items()))


def model_grid(config: dict) -> tuple[str, int, int]:
    """The configured (rung, latitudes, longitudes), refusing a stale pairing.

    `config/planet.yaml` carries all three because `read_sra` validates every
    staged surface file against `latitudes` and `longitudes`, so a resolution
    changed without them refuses its own inputs. That makes the three a single
    fact written three times, and this is where it is checked rather than
    trusted -- the same shape as gravity and mass, which `run_exoplasim.py`
    refuses when they disagree.
    """
    model = config["model"]
    rung = str(model["resolution"]).upper()
    nlat, nlon, _ = geometry(rung)
    declared = int(model["latitudes"]), int(model["longitudes"])
    if declared != (nlat, nlon):
        raise RuntimeError(
            f"config/planet.yaml sets model.resolution {rung} with "
            f"latitudes {declared[0]} and longitudes {declared[1]}, and {rung} "
            f"is {nlat} by {nlon}. The three are one fact: fix whichever is "
            "stale rather than leaving a truncation paired with another "
            "grid's dimensions.")
    return rung, nlat, nlon


def fft_module(rung: str) -> str:
    """The Fortran FFT module that can transform this rung's longitudes.

    `fftmod` is the radix 8-4-3-2 transform: `gp2fc` does one pass of radix 8,
    then radix 4 while four or more remain, then a single radix 3 or radix 2
    tail. It therefore transforms `8 * 4**k * r` with `r` in 1, 2, 3, and
    NOTHING ELSE -- its `nallowed` table is that set, not an independent fact.
    `fft991mod` is the FFT991 package, whose `set99` factorises with 8, 6, 5,
    4, 3 and 2 and covers the lengths `fftmod` cannot.

    DERIVED RATHER THAN LISTED, because a table is what went wrong: the
    postprocessor named the two longitude counts that need `fft991mod` and
    then tested them against a LATITUDE, so T63 and T106 both selected the
    module that cannot transform them and the extension's bare Fortran `stop`
    killed the interpreter with no traceback. world-i38.

    `vendor/exoplasim/exoplasim/pyburn.py:_fftmodule` states the same rule for
    the postprocessor's f2py extensions. That package stays importable on its
    own and cannot import this one, so the rule is written twice on purpose;
    each names the other.
    """
    _, nlon, _ = geometry(rung)
    return "fftmod" if _fftmod_can_transform(nlon) else "fft991mod"


def _fftmod_can_transform(nlon: int) -> bool:
    """Whether `fftmod`'s radix loop covers `nlon` columns. See `fft_module`."""
    n = int(nlon)
    if n % 8:
        return False
    r = n // 8
    while r % 4 == 0:
        r //= 4
    return r in (1, 2, 3)


# THE TWO RESTATEMENTS THAT CANNOT BE REMOVED. The model's build system and the
# vendored ExoPlaSim package each need the ladder before this module is
# reachable: CMake has no Python, and `vendor/exoplasim/exoplasim/__init__.py`
# has to stay importable as an installed package that knows nothing about this
# repository. So the ladder is written three times on purpose, and the other two
# are CHECKED against this one rather than trusted. T31 is why: CMake's list
# omitted 48 while this table declared T31, so `build_model.py --res T31`
# resolved cleanly and then died at cmake configure, and the package's dispatch
# chain accepted seven rungs where the ladder declares eight. world-ajx.
RESTATEMENTS = (
    ("vendor/exoplasim/exoplasim/plasim/CMakeLists.txt",
     r'require_one_of\(PLASIM_NLAT\s+"\$\{PLASIM_NLAT\}"\s+([0-9 ]+)\)'),
    ("vendor/exoplasim/exoplasim/__init__.py",
     r'RESOLUTIONS = (\{[^}]*\})'),
)


def check_restatements(root) -> list[str]:
    """Every place the ladder is restated, against this table. Empty when they agree.

    A CHECK WITH A RIGHT ANSWER rather than a convention: each restatement has
    to name exactly the latitude counts in `RUNGS`, and any disagreement is
    reported as the two sets. Takes the repository root rather than resolving
    one, so this module keeps knowing nothing but the ladder.
    """
    import ast
    import re
    from pathlib import Path

    want = set(RUNGS.values())
    problems = []
    for rel, pattern in RESTATEMENTS:
        path = Path(root) / rel
        if not path.is_file():
            problems.append(f"{rel} is not there, and it restates the ladder")
            continue
        m = re.search(pattern, path.read_text(encoding="utf-8"))
        if m is None:
            problems.append(
                f"{rel} no longer carries the ladder in the shape this check "
                f"reads ({pattern!r}); it restates the ladder and cannot go "
                f"unchecked")
            continue
        body = m.group(1)
        if body.lstrip().startswith("{"):
            got = set(ast.literal_eval(body).values())
        else:
            got = {int(t) for t in body.split()}
        if got != want:
            problems.append(
                f"{rel} allows latitudes {sorted(got)} and the ladder is "
                f"{sorted(want)}: missing {sorted(want - got)}, extra "
                f"{sorted(got - want)}")
    return problems
