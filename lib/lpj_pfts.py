"""Read the LPJ-GUESS plant functional types by NAME and LIFEFORM, from the
file the model itself reads.

`biosphere/generated/vesper_pfts.ins` is the instruction file LPJ-GUESS parses
at startup, so it is the only statement of which types exist and which of them
are grass that cannot disagree with the run. Every consumer that needs to split
cover into tree and grass asks here rather than carrying a copy.

That matters because the alternative is failure-modes class 1 -- one quantity,
several consumers, and a fix that reaches some of them. The split was stated in
four places at once: `build_surface_albedo.py`, which mixes tree and grass
cover into the surface albedo at different albedos; `score_prediction.py`,
whose structural predictions 5, 8 and 9 are all tree-against-grass;
`derive_cover_tolerance.py`, which prices the tolerance the acceptance contract
carries; and `biosphere/config/equilibrium_window.yaml`, which is that
contract. A fifth type of grass added to this world would have had to be found
in all four, and missing one splits cover silently rather than failing.

WHY ONLY THE GRASSES ARE LISTED ANYWHERE. Every consumer takes trees as the
COMPLEMENT -- everything that is not grass and not the `Total` column -- and
`equilibrium_window.yaml` says why in its own comment: "a list of a dozen tree
codes silently misses a newly added one". `trees()` is here for a caller that
wants the set stated positively, but a caller splitting a table should keep
using the exclusion, because the table's columns are the run's and this file's
types are the configuration's, and a column with no matching type must land
somewhere rather than vanish.

GROUPS AND PFTS SHARE A NAMESPACE and are resolved separately. The file carries
both `group "C3G"` and `pft "C3G"`, so a single mapping keyed on the name loses
the group -- and since the lifeform keyword sits in the GROUP, every grass then
comes back with no lifeform at all. Keeping one dict is how the first version
of this reader returned an empty grass set on a file that plainly declares two.
"""

from __future__ import annotations

import re
from pathlib import Path

import builds

DEFAULT_PFT_FILE = builds.PROJECT_ROOT / "biosphere" / "generated" / "vesper_pfts.ins"

# The lifeform keywords LPJ-GUESS accepts in a group or PFT block. Bare tokens,
# not `key value` pairs, which is why they are matched as whole lines.
LIFEFORMS = ("tree", "grass")


def _parse(text: str) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """(groups, pfts), each name to the statements in its block.

    Takes TEXT and not a path, so the parse is testable without a filesystem
    and this module writes nothing. Comments are stripped on `!`, which is the
    file's comment marker in both the `!//` banner form and the trailing form
    the generator writes onto every rescaled line.
    """
    groups: dict[str, list[str]] = {}
    pfts: dict[str, list[str]] = {}
    body: list[str] | None = None
    name = kind = ""
    for raw in text.splitlines():
        line = raw.split("!")[0].rstrip()
        opened = re.match(r'\s*(group|pft)\s+"([\w.]+)"', line)
        if opened:
            kind, name, body = opened.group(1), opened.group(2), []
            continue
        if body is not None and line.strip() == ")":
            (groups if kind == "group" else pfts)[name] = body
            body = None
            continue
        if body is not None and line.strip():
            body.append(line.strip())
    return groups, pfts


def _blocks(path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """`_parse` on a file, refusing an absent one by naming its generator."""
    path = Path(path)
    if not path.is_file():
        raise SystemExit(
            f"{path} does not exist; generate it with "
            "biosphere/scripts/build_vesper_pfts.py")
    return _parse(path.read_text(encoding="utf-8"))


# A numeric statement in a block: `key value`, with the comment already stripped.
_NUMERIC = re.compile(r"^([A-Za-z_]\w*)\s+(-?[\d.]+(?:[eE][-+]?\d+)?)$")


def parameters(name: str, path: Path | None = None) -> dict[str, float]:
    """Every numeric parameter one type resolves to, group references expanded.

    A PFT block states some values itself and inherits the rest by naming a
    group, which may name further groups; a later statement overrides an
    earlier one, which is the order LPJ-GUESS itself reads them in. So a
    consumer that greps the block for `tcmin_surv` finds nothing for the ten
    types that inherit it, and one that greps the whole file finds whichever
    copy comes last.

    A GROUP AND A PFT MAY SHARE A NAME, and `pft "C3G" ( C3G )` is exactly that
    -- the type carries no statements of its own and refers to the group beside
    it. So the lookup takes which of the two a reference means rather than
    preferring one, and preferring the PFT resolved C3G to an empty block that
    reported no cold limit for the type holding most of the polar cover.
    """
    groups, pfts = _blocks(Path(path or DEFAULT_PFT_FILE))
    if name not in pfts and name not in groups:
        raise SystemExit(f"no plant functional type or group named {name!r}")
    return _parameters(name, groups, pfts, as_group=False, seen=frozenset())


def _parameters(name, groups, pfts, as_group, seen) -> dict[str, float]:
    lines = groups.get(name, []) if as_group else pfts.get(name, groups.get(name, []))
    resolved: dict[str, float] = {}
    for line in lines:
        statement = line.split("!")[0].strip()
        if not statement:
            continue
        if statement in groups and statement not in seen:
            resolved.update(_parameters(statement, groups, pfts, True,
                                        seen | {statement}))
            continue
        matched = _NUMERIC.match(statement)
        if matched:
            resolved[matched.group(1)] = float(matched.group(2))
    return resolved


def lifeforms(path: Path | None = None) -> dict[str, str]:
    """{plant functional type: "tree" or "grass"} for every declared type.

    RAISES on a type whose lifeform cannot be resolved, rather than dropping it.
    A type missing from this mapping would be silently counted as a tree by
    every consumer, since they all take trees as the complement of grass.
    """
    path = Path(path or DEFAULT_PFT_FILE)
    return _lifeforms(*_blocks(path), where=path)


def _lifeforms(groups: dict[str, list[str]], pfts: dict[str, list[str]],
               where: object) -> dict[str, str]:
    """The lifeform of every parsed type, refusing one that does not resolve."""

    def resolve(block: list[str], seen: frozenset[str] = frozenset()) -> str | None:
        for keyword in LIFEFORMS:
            if keyword in block:
                return keyword
        # Later includes win, which is the file's own precedence: a block
        # states its general group first and narrows afterwards.
        for statement in reversed(block):
            if statement in groups and statement not in seen:
                found = resolve(groups[statement], seen | {statement})
                if found:
                    return found
        return None

    resolved = {name: resolve(block) for name, block in pfts.items()}
    missing = sorted(name for name, form in resolved.items() if form is None)
    if missing:
        raise SystemExit(
            f"{where} declares {', '.join(missing)} with no resolvable lifeform. "
            "Every consumer takes trees as the complement of grass, so a type "
            "with no lifeform would be counted as a tree without saying so.")
    return resolved


def grass(path: Path | None = None) -> tuple[str, ...]:
    """The grass types, in the file's own declaration order."""
    return tuple(name for name, form in lifeforms(path).items() if form == "grass")


def trees(path: Path | None = None) -> tuple[str, ...]:
    """The tree types, in the file's own declaration order.

    For a caller that wants the set stated positively. A caller splitting a
    TABLE should take the complement of `grass()` instead; see the module note.
    """
    return tuple(name for name, form in lifeforms(path).items() if form == "tree")


def groups(path: Path | None = None) -> tuple[str, ...]:
    """Every group the file declares, in declaration order.

    The file carries a whole taxonomy beside the lifeform -- `boreal`,
    `temperate`, `tropical`, `needleleaved`, `broadleaved`, `evergreen`,
    `summergreen` and the three shade-tolerance classes -- and a consumer
    wanting any of those should ask rather than list the members.
    """
    return tuple(_blocks(Path(path or DEFAULT_PFT_FILE))[0])


def members(*group_names: str, path: Path | None = None) -> tuple[str, ...]:
    """The types belonging to ALL of the named groups, in declaration order.

    Intersection rather than union, because the groupings this file declares
    are orthogonal axes and the useful sets are their conjunctions: "boreal
    needleleaf" is `members("boreal", "needleleaved")`, and stating it as a
    list of three codes is a restatement that a fourth boreal needleleaf type
    would silently fall out of.

    REFUSES a group the file does not declare. A misspelling would otherwise
    return the empty set, and an empty set flows into a caller as "no type is
    boreal" rather than as an error.
    """
    path = Path(path or DEFAULT_PFT_FILE)
    return _members(*_blocks(path), group_names=group_names, where=path)


def _members(group_blocks: dict[str, list[str]], pfts: dict[str, list[str]],
             *, group_names: tuple[str, ...], where: object) -> tuple[str, ...]:
    unknown = [name for name in group_names if name not in group_blocks]
    if unknown:
        raise SystemExit(
            f"{where} declares no group named {', '.join(unknown)}. It declares "
            f"{', '.join(group_blocks)}.")

    def includes(block: list[str], seen: frozenset[str] = frozenset()) -> set[str]:
        reached = set()
        for statement in block:
            if statement in group_blocks and statement not in seen:
                reached.add(statement)
                reached |= includes(group_blocks[statement], seen | {statement})
        return reached

    return tuple(name for name, block in pfts.items()
                 if set(group_names) <= includes(block))


def names(path: Path | None = None) -> tuple[str, ...]:
    """Every declared type, in the file's own declaration order."""
    return tuple(lifeforms(path))


def _selftest() -> int:
    """Exercised on strings, so this module writes nothing anywhere.

    The parse takes text for exactly this reason: a reader that needed a
    filesystem to be tested would have to write fixtures, and a module that
    writes is a module `smoke_test.py` has to ask whether `config/pipeline.yaml`
    declares as a generator. This one declares no artifact because it makes none.
    """
    checks: list[tuple[str, bool]] = []
    live = lifeforms()
    checks.append(("the live file resolves every declared type",
                   all(live.values()) and len(live) >= 2))
    checks.append(("grass and trees partition the declared types",
                   set(grass()) | set(trees()) == set(names())
                   and not set(grass()) & set(trees())))
    checks.append(("the grasses are found through the group of the same name",
                   set(grass()) == {"C3G", "C4G"}))
    checks.append(("boreal needleleaf is the conjunction of two groups",
                   members("boreal", "needleleaved") == ("BNE", "BINE", "BNS")))
    checks.append(("the same conjunction excludes the boreal BROADLEAF",
                   "IBS" in members("boreal")
                   and "IBS" not in members("boreal", "needleleaved")))

    def forms(text: str) -> dict[str, str]:
        return _lifeforms(*_parse(text), where="<fixture>")

    def group_of(text: str, *wanted: str) -> tuple[str, ...]:
        return _members(*_parse(text), group_names=wanted, where="<fixture>")

    def refuses(call) -> bool:
        try:
            call()
        except SystemExit:
            return True
        return False

    # A pft whose lifeform is reachable ONLY through a group sharing its name,
    # which is the arrangement the live file uses for both grasses.
    checks.append(("a group and pft sharing a name do not collide",
                   forms('group "G" (\n\tgrass\n)\n\npft "G" (\n\tG\n)\n') == {"G": "grass"}))
    checks.append(("a type with no resolvable lifeform is refused",
                   refuses(lambda: forms('pft "X" (\n\tinclude 1\n)\n'))))
    # The value every consumer of a bioclimatic limit reads. C3G states none of
    # these itself, so a resolver that stopped at the PFT block would report the
    # type holding most of the polar cover as having no cold limit at all.
    c3g = parameters("C3G")
    checks.append(("an inherited limit resolves through the same-named group",
                   c3g.get("tcmin_surv") == -1000.0 and c3g.get("pstemp_high") == 30.0))
    checks.append(("a type stating its own value overrides the group's",
                   parameters("BNS").get("tcmin_surv") == -1000.0
                   and parameters("BNE").get("tcmin_surv") == -31.0))
    checks.append(("an unknown type is refused rather than resolving to nothing",
                   refuses(lambda: parameters("NOSUCHPFT"))))
    checks.append(("a lifeform is followed through nested groups",
                   forms('group "common" (\n\ttree\n)\n\n'
                         'group "boreal" (\n\tcommon\n)\n\n'
                         'pft "BNE" (\n\tboreal\n)\n') == {"BNE": "tree"}))
    checks.append(("group membership is followed transitively",
                   group_of('group "a" (\n\ttree\n)\n\n'
                            'group "b" (\n\ta\n)\n\n'
                            'pft "P" (\n\tb\n)\n', "a") == ("P",)))
    checks.append(("a misspelled group is refused, not answered empty",
                   refuses(lambda: members("borreal"))))
    checks.append(("an absent file is refused with the generator named",
                   refuses(lambda: lifeforms(Path("/nonexistent/vesper_pfts.ins")))))

    for label, ok in checks:
        print(f"[{' ok ' if ok else 'FAIL'}] {label}")
    failures = sum(not ok for _, ok in checks)
    print(f"\n{failures} failures")
    return failures


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--pft-file", type=Path, default=None)
    args = parser.parse_args()
    if args.self_test:
        return 1 if _selftest() else 0
    for name, form in lifeforms(args.pft_file).items():
        print(f"{name:8s} {form}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
