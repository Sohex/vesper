"""Check a climate product against the ecological forcing field contract.

Worldbuilding frame: Vesper is a simulated super-Earth, its climate is
ExoPlaSim and its biosphere is LPJ-GUESS. Everything below is about the fields
those two models hand each other.

`biosphere/config/ecological_forcing_contract.yaml` is the versioned contract
and `biosphere/notes/ecological-forcing-field-contract.md` is its argument; this
is the half of both that can fail. It does two things on every invocation.

THE DECLARATION, against the producer that has to satisfy it. Every field row
names a service-format code, and `outmod.f90:ecogp` is where those codes are
written, so the two can be checked against each other: a declared code the model
does not write is a contract naming a field nobody produces, and a code the
model writes that no row declares is a field crossing the seam with no agreed
meaning. It also refuses a row missing any column of the declared schema, a
process whose required field is not a declared field, a cadence that is neither
a declared class nor the `undeclared` sentinel with an owner, and an identity
naming a check this module does not implement.

THE PRODUCT, against the identities. The identities in the contract's own list
run against an actual product, and a set of reduced fixtures built to be wrong
in a named way runs on every invocation, so a fixture that does not get the
verdict it was built for is a defect in the checker rather than in the product.

    python biosphere/scripts/check_forcing_contract.py                    # declaration and fixtures
    python biosphere/scripts/check_forcing_contract.py --climatology X.nc # and a product

THE CONTRACT IS NOT THE TRANSPORT. Nothing here says how the artifact is
stored, and nothing here may be relaxed because a container would find a row
inconvenient to carry. EFOR-3 builds the artifact and EFOR-4 the reader; both
are held to the declaration.

Exit 0 when every fixture got its expected verdict and, where a product was
given, every REQUIRED check on it passed. A check can also report ABSENT, which
is what a field the producer does not write earns: it is a fact about the
product, not a pass and not a failure, and the contract forbids naming such a
field as though it were delivered.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, PROJECT_ROOT  # noqa: F401

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "ecological_forcing_contract.yaml"

PASS = "ok"
FAIL = "FAIL"
ABSENT = "absent"

# Every field the contract's rows resolve to on the ExoPlaSim side, with the
# expected sign of each. `+` means the values must be non-negative, `-` that
# they must be non-positive, `0` that no sign is asserted. The signs are not
# cosmetic: ExoPlaSim stores upward fluxes negative (radmod.f90 forms the
# shortwave upward flux as `zfu1 = -...` and then `dswfl = dfu + dfd`), so
# incident surface shortwave is `rss - ssru` and getting that sign wrong is a
# factor of two plus a flip.
SIGNS = {
    "rss": "+",     # surface net shortwave, positive downward
    "rst": "+",     # TOA net shortwave, positive downward
    "ssru": "-",    # surface shortwave upward
    "rsut": "-",    # TOA shortwave upward
    "stru": "-",    # surface longwave upward
    "hfls": "-",    # surface latent heat, upward
    "evap": "-",    # evaporation, upward
    "pr": "+",
    "prl": "+",
    "prc": "+",
    "prsn": "+",
    "czen": "+",
}

# Fields the contract needs whose presence in a given product is not
# guaranteed. Reported, never silently tolerated: a contract that names a field
# a product does not carry has not closed anything. Each entry says WHY it can
# be missing, and the two reasons are different -- one is a field nothing in the
# model writes, the other a field a product predates.
MAY_BE_ABSENT = {
    "tasmax": "code 201 (atsama), written by outmod.f90 all along. In "
              "pyburn.ilibrary and in run_exoplasim.REGULAR_CODES since "
              "world-j0az, so a product postprocessed before that does not "
              "carry it and one postprocessed after does.",
    "tasmin": "code 202 (atsami). Same.",
    "uas": "code 165. Never written; there is no 10 m wind in this model's "
           "output. The mean near-surface wind SPEED is the ecological "
           "stream's `ecowind`, code 614, and not a component pair.",
    "vas": "code 166. Never written.",
}

REQUIRED = ["tas", "ts", "pr", "prl", "prc", "prsn", "rss", "ssru", "rls",
            "stru", "ps", "hus", "maxt", "mint", "lsm"]


def _tol(field: np.ndarray) -> float:
    """Closure tolerance: relative to the field's own scale, not absolute.

    A precipitation rate is 1e-8 m/s and a radiative flux is 1e2 W/m2, so one
    absolute tolerance cannot serve both without either passing everything or
    failing everything.
    """
    scale = float(np.max(np.abs(field))) if field.size else 1.0
    return max(scale, 1.0) * 1e-9


def check_signs(data: dict) -> list[tuple[str, str, str]]:
    out = []
    for name, sign in SIGNS.items():
        if name not in data:
            continue
        v = np.asarray(data[name], dtype=float)
        bad = (v < -_tol(v)).sum() if sign == "+" else (
            (v > _tol(v)).sum() if sign == "-" else 0)
        out.append((
            f"{name} is {'non-negative' if sign == '+' else 'non-positive'}",
            PASS if bad == 0 else FAIL,
            "" if bad == 0 else f"{int(bad)} values on the wrong side of zero"))
    return out


def check_precipitation(data: dict) -> list[tuple[str, str, str]]:
    out = []
    if {"pr", "prl", "prc"} <= data.keys():
        pr = np.asarray(data["pr"], float)
        residual = float(np.max(np.abs(pr - np.asarray(data["prl"], float)
                                       - np.asarray(data["prc"], float))))
        ok = residual <= _tol(pr)
        out.append(("large-scale + convective closes on the total",
                    PASS if ok else FAIL, f"max residual {residual:.3e}"))
    if {"pr", "prsn"} <= data.keys():
        pr = np.asarray(data["pr"], float)
        prsn = np.asarray(data["prsn"], float)
        excess = float(np.max(prsn - pr))
        ok = excess <= _tol(pr)
        out.append((
            "snowfall is a SUBSET of the total, never a third addend",
            PASS if ok else FAIL,
            f"max prsn - pr = {excess:.3e}. rainmod.f90 forms dprs from the "
            f"same zprsc and zprsl that are already inside dprc and dprl, so "
            f"prl + prc + prsn double counts the snow."))
    return out


def check_shortwave(data: dict) -> list[tuple[str, str, str]]:
    if not {"rss", "ssru"} <= data.keys():
        return []
    rss = np.asarray(data["rss"], float)
    ssru = np.asarray(data["ssru"], float)
    incident = rss - ssru
    ok = bool(np.all(incident >= -_tol(rss)))
    out = [("incident shortwave = rss - ssru is non-negative",
            PASS if ok else FAIL,
            f"range {float(incident.min()):.3f} to {float(incident.max()):.3f} "
            f"W/m2")]
    lit = incident > 1.0
    if lit.any():
        implied = (-ssru[lit]) / incident[lit]
        ok = bool(np.all((implied >= -1e-9) & (implied <= 1.0 + 1e-9)))
        out.append(("the implied surface albedo lies in [0, 1]",
                    PASS if ok else FAIL,
                    f"{float(implied.min()):.3f} to {float(implied.max()):.3f}"))
    return out


# Each extremum pair with the variable it is an extremum OF. maxt/mint are
# extrema of the SURFACE temperature and tasmax/tasmin of the near-surface AIR
# temperature; the contract names them apart for exactly this reason, and each
# is held only against its own variable. Requiring either pair to bracket the
# other temperature would fail on correct data: the two variables differ by the
# dry-adiabatic reduction fluxmod applies, and on the bootstrap climatology tas
# lies outside [mint, maxt] in 15,561 of 24,576 cell-bins by up to 28.3 K.
EXTREMA_PAIRS = (("maxt", "mint", "ts"),
                 ("tasmax", "tasmin", "tas"))


def check_extrema(data: dict) -> list[tuple[str, str, str]]:
    """Is each extremum pair an extremum of the variable it is named for?

    An extremum over a window brackets the mean of the same variable over the
    same window, whatever the variable does inside it. So this identifies the
    variable rather than merely testing a plausible bound, and it is the check
    the wrong-variable defect was found by.
    """
    out = []
    for himax, lomin, owns in EXTREMA_PAIRS:
        if not {himax, lomin} <= data.keys():
            continue
        lo = np.asarray(data[lomin], float)
        hi = np.asarray(data[himax], float)
        out.append((f"{himax} >= {lomin}",
                    PASS if bool(np.all(hi >= lo - 1e-6)) else FAIL,
                    f"worst {himax} - {lomin} = {float(np.min(hi - lo)):.3f} K"))
        if owns not in data:
            continue
        v = np.asarray(data[owns], float)
        outside = int(np.sum((v < lo - 1e-6) | (v > hi + 1e-6)))
        out.append((
            f"{himax} and {lomin} bracket {owns}",
            PASS if outside == 0 else FAIL,
            "" if outside == 0 else
            f"{outside} of {v.size} outside, worst excursion "
            f"{float(np.max(np.maximum(lo - v, v - hi))):.3f} K"))
    return out


def check_declared_units(data: dict, units: dict) -> list[tuple[str, str, str]]:
    """A declared unit that the values contradict.

    pyburn declares relative humidity's units as `1` and computes it as a
    percentage clipped to [0, 100]. A consumer that trusts the attribute is out
    by a hundred, which is why the contract carries specific humidity instead.
    """
    out = []
    for name, unit in units.items():
        if name not in data or unit not in ("1", "nondimen"):
            continue
        v = np.asarray(data[name], float)
        top = float(np.max(np.abs(v))) if v.size else 0.0
        out.append((f"{name} declares units '{unit}' and stays within them",
                    PASS if top <= 1.5 else FAIL,
                    "" if top <= 1.5 else
                    f"reaches {top:.4g}, which is a percentage and not a "
                    f"fraction"))
    return out


def check_finite(data: dict) -> list[tuple[str, str, str]]:
    bad = sorted(n for n, v in data.items()
                 if not np.all(np.isfinite(np.asarray(v, float))))
    return [("every field is finite", PASS if not bad else FAIL,
             "" if not bad else f"non-finite in {bad}")]


def check_present(data: dict) -> list[tuple[str, str, str]]:
    out = []
    missing = [n for n in REQUIRED if n not in data]
    out.append(("every field the contract requires is present",
                PASS if not missing else FAIL,
                "" if not missing else f"missing {missing}"))
    for name, why in MAY_BE_ABSENT.items():
        if name in data:
            out.append((f"{name} is available", PASS, ""))
        else:
            out.append((f"{name} is available", ABSENT, why))
    return out


# ---------------------------------------------------------------------------
# The declaration, against the producer that has to satisfy it.
# ---------------------------------------------------------------------------

# Every column a field row has to carry. The schema block in the declaration
# argues each one; this is the list the rows are held to, and a row missing any
# of them is refused rather than defaulted, because every default here would be
# a meaning nobody agreed.
ROW_COLUMNS = ("name", "producer_code", "producer_expression", "contract_unit",
               "time_base", "kind", "semantics", "area_basis", "sign",
               "converts")

KINDS = ("mean", "extremum", "instantaneous")
SEMANTICS = ("intensive", "extensive", "extremal")

UNDECLARED = "undeclared"


def _split_args(text: str) -> list[str]:
    """Split a Fortran argument list on top-level commas.

    Nested calls are the reason: `writescalar(143,real(nstep+1-naccueco)*deltsec,600)`
    has a comma-free inner call, but a naive split on the first `)` truncates
    the argument and loses the code entirely.
    """
    parts, depth, current = [], 0, []
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        if char == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    parts.append("".join(current).strip())
    return parts


def _ecogp_codes(source_text: str) -> set[int]:
    """Every service-format code the ecological stream writes.

    Parsed from the emitter itself rather than from a list kept beside it: a
    list would be a second declaration of the same fact, and the failure this
    check exists for is the two disagreeing.

    The code sits in a different argument position in the two writers.
    `writegp(unit, field, code, level)` puts it third and a LEVEL fourth, so
    taking the last numeric argument would read every field as level zero;
    `writescalar(unit, value, code)` puts it last.
    """
    start = source_text.find("subroutine ecogp")
    if start < 0:
        return set()
    end = source_text.find("subroutine ecoreset", start)
    body = source_text[start:end if end > 0 else len(source_text)]
    codes = set()
    for match in re.finditer(r"call\s+write(gp|scalar)\s*\((.*)\)\s*$",
                             body, flags=re.M):
        args = _split_args(match.group(2))
        if len(args) < 3 or args[0].strip() != "143":
            continue
        raw = args[2] if match.group(1) == "gp" else args[-1]
        if re.fullmatch(r"\d+", raw.strip()):
            codes.add(int(raw.strip()))
    return codes


def check_declaration(declaration: dict, source_text: str
                      ) -> list[tuple[str, str, str]]:
    """The contract against the producer, and against itself."""
    out: list[tuple[str, str, str]] = []

    def verdict(label: str, bad: list[str]) -> None:
        out.append((label, PASS if not bad else FAIL,
                    "" if not bad else "; ".join(bad)))

    version = declaration.get("version")
    verdict("the contract declares an integer version",
            [] if isinstance(version, int) and version >= 1 else
            [f"version is {version!r}"])

    rows = declaration.get("fields", [])
    interval = declaration.get("interval_record", [])

    bad = [f"{row.get('name', '?')} has no {column}"
           for row in rows for column in ROW_COLUMNS
           if column not in row or row[column] in (None, "")]
    verdict("every field row carries every column of the declared schema", bad)

    bad = [f"{row['name']} declares kind {row.get('kind')!r}"
           for row in rows if row.get("kind") not in KINDS]
    bad += [f"{row['name']} declares semantics {row.get('semantics')!r}"
            for row in rows if row.get("semantics") not in SEMANTICS]
    verdict("every field row declares a known kind and semantics", bad)

    declared_extensive = [row["name"] for row in rows
                          if row.get("semantics") == "extensive"]
    claim = declaration.get("extensive_fields_carried")
    bad = []
    if claim == "none" and declared_extensive:
        bad = [f"declares no extensive field and {declared_extensive} are extensive"]
    elif claim != "none" and not declared_extensive:
        bad = [f"declares extensive_fields_carried {claim!r} and no row is extensive"]
    verdict("the extensive-field claim matches the rows", bad)

    seen_names: set[str] = set()
    seen_codes: dict[int, str] = {}
    bad = []
    for row in rows:
        if row["name"] in seen_names:
            bad.append(f"two rows named {row['name']}")
        seen_names.add(row["name"])
        code = row.get("producer_code")
        if code in seen_codes:
            bad.append(f"{row['name']} and {seen_codes[code]} both claim code {code}")
        elif code is not None:
            seen_codes[code] = row["name"]
    verdict("no two rows share a name or a producer code", bad)

    written = _ecogp_codes(source_text)
    verdict("the emitter was found and writes codes",
            [] if written else ["outmod.f90:ecogp parsed to no codes at all"])

    declared_codes = {row["producer_code"] for row in rows
                      if row.get("producer_code") is not None}
    declared_codes |= {row["producer_code"] for row in interval
                       if row.get("producer_code") is not None}
    if written:
        missing = sorted(declared_codes - written)
        verdict("every declared code is one the model writes",
                [] if not missing else
                [f"codes {missing} are declared and outmod.f90:ecogp writes none of them"])
        undeclared = sorted(written - declared_codes)
        verdict("every code the stream writes is declared",
                [] if not undeclared else
                [f"codes {undeclared} cross the seam with no row and so with no agreed meaning"])

    bad = [f"{row.get('name', '?')} is {row.get('status')!r}"
           for row in interval
           if row.get("status") not in ("produced", "derived", "not_produced")]
    bad += [f"{row['name']} is not produced and names no owner"
            for row in interval
            if row.get("status") == "not_produced" and not row.get("owner")]
    verdict("every interval-record row declares a status, and an absent one an owner",
            bad)

    classes = set(declaration.get("cadence_classes", {}))
    bad = []
    for process in declaration.get("processes", []):
        pid = process.get("id", "?")
        for field in process.get("fields", []):
            if field not in seen_names:
                bad.append(f"{pid} requires {field}, which is not a declared field")
        cadence = process.get("cadence")
        if cadence == UNDECLARED:
            if not process.get("owner"):
                bad.append(f"{pid} leaves its cadence undeclared and names no owner")
            if not process.get("why_undeclared"):
                bad.append(f"{pid} leaves its cadence undeclared and does not say why")
        elif cadence not in classes:
            bad.append(f"{pid} declares cadence {cadence!r}, which is not a declared class")
        elif not process.get("source"):
            bad.append(f"{pid} declares a cadence and names no source for it")
    verdict("every process names declared fields and a cadence it can account for", bad)

    implemented = {name for name in globals() if name.startswith("check_")}
    bad = []
    for identity in declaration.get("identities", []):
        target = identity.get("check")
        if target == UNDECLARED:
            if not identity.get("owner"):
                bad.append(f"{identity.get('id')} has no check and names no owner")
        elif target not in implemented:
            bad.append(f"{identity.get('id')} names {target!r}, which this module does not implement")
        if not identity.get("statement"):
            bad.append(f"{identity.get('id')} states no identity")
    verdict("every identity names a check this module implements, or an owner", bad)

    bad = [f"{entry.get('what', '?')} is not declared disposable"
           for entry in declaration.get("legacy", [])
           if entry.get("disposable") is not True or not entry.get("why")]
    verdict("every legacy transport is declared disposable, with a reason", bad)

    known = set(REQUIRED) | set(MAY_BE_ABSENT) | set(SIGNS)
    bad = [f"{row['name']} maps to product field {row['product_name']!r}, "
           "which this module's product checks do not know"
           for row in rows
           if row.get("product_name") and row["product_name"] not in known]
    verdict("every row's product name is one the product checks know", bad)

    return out


def run_declaration(declaration: dict | None = None,
                    source_text: str | None = None
                    ) -> list[tuple[str, str, str]]:
    if declaration is None:
        declaration = yaml.safe_load(DECLARATION.read_text())
    if source_text is None:
        source_text = (PROJECT_ROOT / declaration["producer"]["source_file"]
                       ).read_text()
    return check_declaration(declaration, source_text)


# Declaration fixtures. The first is the declaration as it stands and has to
# come back clean; every other is built to be wrong in a named way.
#
# The expectation is a SET of check names and not one, on the same terms as the
# product fixtures: a corruption legitimately breaks the checks derived from it
# as well. Moving one row's code onto another's takes a code out of the declared
# set, so the model-writes-it check and the every-code-declared check both fire,
# and removing the wind row takes a field two processes require with it. Naming
# the whole set is what makes this a test, since a fixture that breaks one MORE
# check than declared is as much a defect as one that breaks one fewer.
def _declaration_fixtures() -> list[dict]:
    import copy

    declaration = yaml.safe_load(DECLARATION.read_text())
    source_text = (PROJECT_ROOT / declaration["producer"]["source_file"]).read_text()

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def row(d, name):
        for entry in d["fields"]:
            if entry["name"] == name:
                return entry
        raise KeyError(name)

    def process(d, pid):
        for entry in d["processes"]:
            if entry["id"] == pid:
                return entry
        raise KeyError(pid)

    def unowned_gap(d):
        """An interval row declared absent from the producer with no owner.

        The wrongness is BUILT rather than borrowed from whichever row happens
        to be unproduced today. When world-wq8i made the orbital position a
        produced record, no row was `not_produced` any more and a fixture that
        only removed an owner stopped being wrong at all -- which is a fixture
        that passes by accident, the one thing this block exists to refuse.
        """
        entry = d["interval_record"][-1]
        entry["status"] = "not_produced"
        entry.pop("owner", None)

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a field row missing a column of the schema",
         mutate(lambda d: row(d, "wind_speed").pop("converts")),
         "every field row carries every column of the declared schema"),
        ("a row declaring a kind the schema does not have",
         mutate(lambda d: row(d, "wind_speed").__setitem__("kind", "total")),
         "every field row declares a known kind and semantics"),
        ("a row declared extensive against the no-extensive-field claim",
         mutate(lambda d: row(d, "total_precipitation").__setitem__("semantics", "extensive")),
         "the extensive-field claim matches the rows"),
        ("two rows claiming one producer code",
         mutate(lambda d: row(d, "wind_speed").__setitem__("producer_code", 613)),
         {"no two rows share a name or a producer code",
          "every code the stream writes is declared"}),
        ("a declared code the model does not write",
         mutate(lambda d: row(d, "wind_speed").__setitem__("producer_code", 999)),
         {"every declared code is one the model writes",
          "every code the stream writes is declared"}),
        ("a code the stream writes that no row declares",
         mutate(lambda d: d["fields"].remove(row(d, "wind_speed"))),
         {"every code the stream writes is declared",
          "every process names declared fields and a cadence it can account for"}),
        ("an interval row absent from the producer and naming no owner",
         mutate(unowned_gap),
         "every interval-record row declares a status, and an absent one an owner"),
        ("a process requiring a field the contract does not declare",
         mutate(lambda d: process(d, "bio-23")["fields"].append("relative_humidity")),
         "every process names declared fields and a cadence it can account for"),
        ("a process on a cadence class the contract does not define",
         mutate(lambda d: process(d, "bio-23").__setitem__("cadence", "hourly")),
         "every process names declared fields and a cadence it can account for"),
        ("a cadence left undeclared with no owner",
         mutate(lambda d: process(d, "bio-13").pop("owner")),
         "every process names declared fields and a cadence it can account for"),
        ("an identity naming a check this module does not implement",
         mutate(lambda d: d["identities"][0].__setitem__("check", "check_nothing")),
         "every identity names a check this module implements, or an owner"),
        ("a legacy transport not declared disposable",
         mutate(lambda d: d["legacy"][0].__setitem__("disposable", False)),
         "every legacy transport is declared disposable, with a reason"),
        ("a version that is not an integer",
         mutate(lambda d: d.__setitem__("version", "1.0")),
         "the contract declares an integer version"),
    ]

    results = []
    for label, candidate, expect in cases:
        failed = {name for name, status, _ in check_declaration(candidate, source_text)
                  if status == FAIL}
        wanted = set() if expect is None else (
            expect if isinstance(expect, set) else {expect})
        results.append({"fixture": label, "expected": sorted(wanted) or "clean",
                        "found": sorted(failed), "pass": failed == wanted})
    return results


def run_checks(data: dict, units: dict | None = None
               ) -> list[tuple[str, str, str]]:
    """Every contract check that can be run on the fields given."""
    units = units or {}
    return (check_present(data) + check_finite(data) + check_signs(data)
            + check_precipitation(data) + check_shortwave(data)
            + check_extrema(data) + check_declared_units(data, units))


# ---------------------------------------------------------------------------
# Fixtures. Seven of the eight are built to be wrong in a named way.
# ---------------------------------------------------------------------------

def _clean() -> dict:
    """A small artifact that satisfies every identity, by construction."""
    rng = np.random.default_rng(20260824)
    shape = (4, 3, 5)
    prl = rng.uniform(0.0, 2e-8, shape)
    prc = rng.uniform(0.0, 4e-8, shape)
    total = prl + prc
    prsn = total * rng.uniform(0.0, 1.0, shape)
    incident = rng.uniform(0.0, 400.0, shape)
    albedo = rng.uniform(0.05, 0.6, shape)
    ts = rng.uniform(240.0, 310.0, shape)
    swing = rng.uniform(0.5, 12.0, shape)
    # The air temperature and its own extrema, built so the air pair brackets
    # tas and does NOT bracket ts: that is what the two pairs being different
    # variables means, and a fixture in which both pairs bracket both
    # temperatures could not tell them apart.
    tas = ts - 20.0 - swing * rng.uniform(0.0, 0.5, shape)
    aswing = rng.uniform(0.5, 8.0, shape)
    return {
        "tas": tas,
        "tasmax": tas + aswing, "tasmin": tas - aswing,
        "ts": ts,
        "pr": total, "prl": prl, "prc": prc, "prsn": prsn,
        "rss": incident * (1.0 - albedo), "ssru": -incident * albedo,
        "rls": rng.uniform(-150.0, 60.0, shape),
        "stru": -rng.uniform(150.0, 600.0, shape),
        "ps": rng.uniform(800.0, 1050.0, shape),
        "hus": rng.uniform(0.0, 0.02, shape),
        "maxt": ts + swing, "mint": ts - swing,
        "lsm": (rng.uniform(0.0, 1.0, shape) > 0.5).astype(float),
        "czen": rng.uniform(0.0, 0.6, shape),
    }


def _fixtures() -> list[tuple[str, dict, dict, set]]:
    """Name, fields, declared units, and the EXACT set of checks that must fail.

    A set rather than one name, because a corruption legitimately breaks the
    checks derived from it as well: storing the upward shortwave positive breaks
    the incident term and the implied albedo along with the sign itself. Naming
    the whole set is what makes this a test rather than a smoke alarm, since a
    fixture that fails one MORE check than declared is as much a defect as one
    that fails one fewer.
    """
    out: list[tuple[str, dict, dict, set]] = [
        ("a clean artifact", _clean(), {}, set()),
    ]

    d = _clean()
    d["pr"] = d["pr"] + d["prsn"]
    out.append(("snowfall added as a third addend", d, {},
                {"large-scale + convective closes on the total"}))

    d = _clean()
    d["prsn"] = d["pr"] * 1.5
    out.append(("more snow than precipitation", d, {},
                {"snowfall is a SUBSET of the total, never a third addend"}))

    d = _clean()
    d["ssru"] = -d["ssru"]
    out.append(("upward shortwave stored positive", d, {},
                {"ssru is non-positive",
                 "incident shortwave = rss - ssru is non-negative",
                 "the implied surface albedo lies in [0, 1]"}))

    d = _clean()
    d["maxt"], d["mint"] = d["mint"], d["maxt"]
    out.append(("extrema swapped", d, {},
                {"maxt >= mint", "maxt and mint bracket ts"}))

    d = _clean()
    # Extrema of a DIFFERENT variable: the defect this checker was written for,
    # and the one that leaves every other identity intact.
    other = d["ts"] + 40.0
    d["maxt"] = other + 1.0
    d["mint"] = other - 1.0
    out.append(("extrema of a different variable", d, {},
                {"maxt and mint bracket ts"}))

    d = _clean()
    # The AIR pair built from the SURFACE temperature: the same defect on the
    # other pair, and the one that would silently restore the wrong diurnal
    # range to climate.dtr now that a product carries both pairs.
    d["tasmax"] = d["ts"] + 1.0
    d["tasmin"] = d["ts"] - 1.0
    out.append(("air extrema of the surface temperature", d, {},
                {"tasmax and tasmin bracket tas"}))

    d = _clean()
    d["hur"] = d["czen"] * 100.0
    out.append(("a percentage declared as a fraction", d, {"hur": "1"},
                {"hur declares units '1' and stays within them"}))
    return out


def run_fixtures(verbose: bool = False) -> int:
    """Every fixture must get the verdict it was built for. Returns failures."""
    failures = 0
    for name, data, units, expected in _fixtures():
        results = {label: (status, detail)
                   for label, status, detail in run_checks(data, units)}
        failed = {label for label, (status, _) in results.items()
                  if status == FAIL}
        ok = failed == expected
        why = "" if ok else (f"expected {sorted(expected)} to fail, "
                             f"got {sorted(failed)}")
        failures += 0 if ok else 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] fixture: {name}"
              + (f"\n           {why}" if why else ""))
        if verbose:
            for label in sorted(expected):
                if results.get(label, ("", ""))[1]:
                    print(f"           {label}: {results[label][1]}")
    return failures


def read_climatology(path: Path) -> tuple[dict, dict]:
    """Fields and their declared units, with level axes reduced to the lowest.

    The lowest MODEL LEVEL, not a 2 m diagnostic: `hus` and the winds carry a
    level axis and the contract's near-surface rows resolve to the bottom of it.
    """
    import netCDF4 as nc
    data: dict = {}
    units: dict = {}
    with nc.Dataset(path) as ds:
        for name in ds.variables:
            if name in ("lat", "lon", "lev", "levp", "time", "modes"):
                continue
            v = np.asarray(ds[name][:], dtype=float)
            if v.ndim == 4:
                v = v[:, -1]        # lowest model level
            if v.ndim != 3:
                continue
            data[name] = v
            units[name] = getattr(ds[name], "units", "")
    return data, units


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None,
                        help="a pyburn climatology or forcing artifact to "
                             "check. Without it only the fixtures run.")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    failures = 0

    declaration = yaml.safe_load(DECLARATION.read_text())
    print(f"DECLARATION  {DECLARATION.relative_to(PROJECT_ROOT)}  "
          f"version {declaration['version']}")
    for label, status, detail in run_declaration(declaration):
        mark = {PASS: "  ok  ", FAIL: " FAIL ", ABSENT: "absent"}[status]
        print(f"[{mark}] {label}" + (f"\n           {detail}" if detail else ""))
        if status == FAIL:
            failures += 1

    dfix = _declaration_fixtures()
    broken = [f for f in dfix if not f["pass"]]
    print(f"\nDECLARATION FIXTURES  ({len(dfix) - 1} of {len(dfix)} are built "
          "to be wrong in a named way)")
    print(f"[{'  ok  ' if not broken else ' FAIL '}] "
          f"{len(dfix) - len(broken)} of {len(dfix)} got their verdict")
    for case in broken:
        print(f"           BROKEN: {case['fixture']}: expected "
              f"{case['expected']!r}, found {case['found']}")
    failures += len(broken)

    print(f"\nPRODUCT FIXTURES  ({len(_fixtures()) - 1} of {len(_fixtures())} are built "
          "to be wrong in a named way)")
    failures += run_fixtures(args.verbose)

    if args.climatology is not None:
        print(f"\nPRODUCT   {args.climatology}")
        data, units = read_climatology(args.climatology)
        for label, status, detail in run_checks(data, units):
            mark = {PASS: "  ok  ", FAIL: " FAIL ", ABSENT: "absent"}[status]
            print(f"[{mark}] {label}" + (f"\n           {detail}" if detail
                                         else ""))
            if status == FAIL:
                failures += 1

    print(f"\n{failures} failed" if failures else "\nno failures")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
