#!/usr/bin/env python3
"""Preflight: does everything in tree describe the same world?

Run this before an expensive run and after changing `source_build`. It exists
because the same class of error kept arriving by different routes, and every
instance was silent: array shapes matched, fields looked plausible, and the
answer was wrong.

    python scripts/check_consistency.py

Exit status is 0 when everything agrees, 1 when something disagrees. Warnings
report things that are stale but harmless to the active configuration.

## What it checks, and which bug each one is for

**Terrain hash across artifacts.** Every product carrying a `terrain_hash` must
match the configured build. The carve verdict once ran with one terrain's basins
and another's coupling matrix and raised an IndexError purely because the counts
differed; had they matched it would have paired row *i* of one build with row *i*
of another and produced a number.

**Stale flat data directories.** Hydrography is per-build now. A leftover
`<component>/data/*.nc` from an earlier build is not an error, but it is what
`world_state.py` was reading when it reported 2,107 basins for a build with
2,540.

**Coupling convention.** `coupling_*.nc` must carry `cell_lon`. Without it the
longitude convention cannot be checked, and the matrix numbers its columns
-180..180 while an ExoPlaSim climatology numbers its own 0..360. Every basin read
its antipode for an entire iteration.

**Surface inputs present and current.** Every code the config says it supplies
must exist, and must not predate the build it claims to describe.

**Carve list completeness.** A list sent to Orogen must cover the whole
catalogue with unique ids, and its carried-forward entries must all be retain 0,
since the loop is only monotone if an already-carved basin stays carved.

**Config internal consistency.** Gravity against mass, and the declared orbit
against the flux it is derived from.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

from netCDF4 import Dataset          # noqa: E402
import numpy as np                   # noqa: E402
import yaml                          # noqa: E402

import builds                        # noqa: E402

FAIL, WARN, OK = "FAIL", "warn", "ok"


class Report:
    def __init__(self):
        self.rows: list[tuple[str, str, str]] = []

    def add(self, status, check, detail=""):
        self.rows.append((status, check, detail))

    @property
    def failed(self) -> bool:
        return any(s == FAIL for s, _, _ in self.rows)

    def show(self):
        width = max(len(c) for _, c, _ in self.rows) + 2
        for status, check, detail in self.rows:
            mark = {OK: "  ok  ", WARN: " warn ", FAIL: " FAIL "}[status]
            print(f"[{mark}] {check:<{width}} {detail}")
        n_fail = sum(1 for s, _, _ in self.rows if s == FAIL)
        n_warn = sum(1 for s, _, _ in self.rows if s == WARN)
        print(f"\n{len(self.rows)} checks, {n_fail} failed, {n_warn} warnings")


def nc_terrain(path: Path) -> str | None:
    try:
        with Dataset(path) as ds:
            return getattr(ds, "terrain_hash", None)
    except OSError:
        return None


def json_terrain(path: Path) -> str | None:
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return d.get("terrain_hash") or (d.get("source") or {}).get("terrain_hash")


def sha256_of(path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    import argparse
    argparse.ArgumentParser(
        description="Verify that everything in tree describes the same world. "
                    "Exit 1 on disagreement. Takes no arguments; it audits the "
                    "artifacts as they are.").parse_args()
    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    rep = Report()
    build = str(config.get("source_build", ""))
    try:
        want = builds.terrain_hash(config)
    except Exception as exc:
        print(f"cannot resolve the configured build {build!r}: {exc}")
        return 1
    rep.add(OK, f"active build", f"{build}  {want[:16]}")

    # -- terrain hash across every artifact that records one -----------------
    data = builds.component_data("hydrography", config)
    checked = 0
    for path in sorted(data.glob("*.nc")):
        got = nc_terrain(path)
        checked += 1
        if got is None:
            rep.add(WARN, f"terrain {path.name}", "carries no terrain_hash")
        elif got != want:
            rep.add(FAIL, f"terrain {path.name}",
                    f"is {got[:16]}, build is {want[:16]}")
        else:
            rep.add(OK, f"terrain {path.name}", got[:16])
    if not checked:
        rep.add(FAIL, "hydrography products", f"none in {data.relative_to(ROOT)}")

    resolution = str(config["model"]["resolution"]).lower()
    for path in sorted((ROOT / "exoplasim" / "inputs" / resolution).glob("*report*.json")):
        got = json_terrain(path)
        if got is None:
            rep.add(WARN, f"terrain {path.name}", "records no terrain_hash")
        elif got != want:
            rep.add(FAIL, f"terrain {path.name}",
                    f"is {got[:16]}, build is {want[:16]}")
        else:
            rep.add(OK, f"terrain {path.name}", got[:16])
    for path in (ROOT / "pedology" / "analysis" / "soil_report.json",):
        got = json_terrain(path)
        if got is not None:
            rep.add(OK if got == want else FAIL, f"terrain {path.name}",
                    got[:16] if got == want else f"is {got[:16]}, build is {want[:16]}")

    # -- gravity, which the terrain hash cannot see --------------------------
    try:
        from orogen import Export as _Export
        ex = _Export(builds.mesh_export(config))
        declared = float(config["planet"]["gravity_m_s2"])
        if abs(ex.gravity_m_s2 - declared) > 1e-6:
            rep.add(FAIL, "gravity build vs config",
                    f"build was generated at {ex.gravity_m_s2} m/s2, config "
                    f"declares {declared}; every vertical km is off by "
                    f"{declared / ex.gravity_m_s2:.4f}x")
        else:
            rep.add(OK, "gravity build vs config", f"{declared} m/s2")
    except Exception as exc:
        rep.add(WARN, "gravity build vs config", f"not checked: {exc}")

    # -- stale flat data -----------------------------------------------------
    flat = ROOT / "hydrography" / "data"
    if data != flat:
        stale = [p for p in flat.glob("*.nc") if nc_terrain(p) not in (None, want)]
        if stale:
            rep.add(WARN, "stale flat hydrography/data",
                    f"{len(stale)} files from another build; scripts defaulting "
                    "there will read the wrong terrain")

    # -- coupling convention -------------------------------------------------
    couplings = sorted(data.glob("coupling_*.nc"))
    if not couplings:
        rep.add(FAIL, "coupling matrices", "none found")
    for path in couplings:
        try:
            with Dataset(path) as ds:
                has = "cell_lon" in ds.variables
        except OSError:
            has = False
        rep.add(OK if has else FAIL, f"convention {path.name}",
                "declares cell_lon" if has else
                "predates cell_lon; its longitude convention cannot be checked")

    # -- surface inputs ------------------------------------------------------
    sys.path.insert(0, str(ROOT / "exoplasim" / "scripts"))
    try:
        import run_exoplasim as rx
        codes = sorted(rx.intended_surface_codes(config))
        missing = [c for c in codes if not rx.surface_sra(config, c).is_file()]
        rep.add(FAIL if missing else OK, "surface inputs",
                f"missing {missing}" if missing else
                f"{len(codes)} present: {codes}")
        if not missing:
            # Currency by CONTENT, via the reports written alongside the .sra
            # files, not by mtime. An mtime comparison against the build
            # manifest fires whenever a file is restored or touched without
            # changing: restoring the export from a backup gave the manifest a
            # new mtime and every surface input read as stale, which is the same
            # false positive the cycle-binary check had. The terrain hash the
            # reports already carry is the thing that actually matters, and it
            # is checked above; this adds the check that a report EXISTS for the
            # active resolution, since a missing report is how an .sra survives
            # a build change unnoticed.
            res = str(config["model"]["resolution"]).lower()
            reports = sorted((ROOT / "exoplasim" / "inputs" / res).glob("*report*.json"))
            stale = []
            for r in reports:
                try:
                    got = json.loads(r.read_text(encoding="utf-8")).get("terrain_hash")
                except (OSError, json.JSONDecodeError):
                    got = None
                if got != want:
                    stale.append(r.name)
            rep.add(FAIL if (stale or not reports) else OK, "surface inputs current",
                    f"reports describe another build: {stale}" if stale else
                    (f"no surface reports under exoplasim/inputs/{res}"
                     if not reports else
                     f"{len(reports)} reports all describe {want[:16]}"))
    except Exception as exc:
        rep.add(WARN, "surface inputs", f"not checked: {exc}")

    # -- carve list ----------------------------------------------------------
    carve = data / "carve_list.json"
    if carve.is_file():
        d = json.loads(carve.read_text(encoding="utf-8"))
        entries = d["basins"]
        ids = [e["id"] for e in entries]
        # Against the basins this build still has plus the ones an earlier pass
        # carved, NOT against the build manifest's preserved list: that is the
        # post-carve survivors, and a list sent to Orogen has to cover the base
        # catalogue, which is larger.
        carried_n = sum(1 for e in entries if e.get("carried_from_previous_pass"))
        try:
            with Dataset(data / "basins.nc") as ds:
                here = ds.dimensions["basin"].size
        except Exception:
            here = None
        dup = len(ids) - len(set(ids))
        rep.add(FAIL if dup else OK, "carve list ids",
                f"{dup} duplicates" if dup else f"{len(ids)} unique")
        bad = [e["id"] for e in entries
               if e.get("carried_from_previous_pass") and float(e["retain"]) != 0.0]
        rep.add(FAIL if bad else OK, "carve list monotone",
                f"{len(bad)} carried entries are not retain 0" if bad else
                "carried entries all retain 0")
        if here is not None:
            expect = here + carried_n
            rep.add(OK if len(ids) == expect else FAIL, "carve list accounts for every basin",
                    f"{len(ids)} entries = {here} decidable + {carried_n} carried"
                    if len(ids) == expect else
                    f"{len(ids)} entries, expected {here} + {carried_n} = {expect}")

    # -- config --------------------------------------------------------------
    try:
        derived = rx.derive(config, float(config["orbit"]["baseline_flux_earth"]))
        rep.add(OK, "config gravity vs mass", "consistent")
        rep.add(OK, "derived orbit",
                f"a = {derived['semimajor_axis_au']:.6f} AU, "
                f"year = {derived['orbital_year_earth_days']:.3f} d")
    except Exception as exc:
        rep.add(FAIL, "config", str(exc))

    # -- config blocks declare their determination status --------------------
    #
    # A value that has never been decided must not be indistinguishable from one
    # that has. stellar_cycle carried an amplitude and a centre from before there
    # was a baseline flux to centre on, and both were later read back as though
    # they had been chosen.
    try:
        import re as _re
        text = (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
        lines = text.splitlines()
        unmarked = []
        for i, line in enumerate(lines):
            if not _re.match(r"^[a-z_]+:\s*$", line):
                continue
            # Scan the WHOLE contiguous comment block above the key, not a fixed
            # number of lines. A fixed window silently fails a block whose
            # rationale grew past it, which is backwards: the blocks with the
            # most to explain are exactly the ones that most need a status, and
            # an 8-line window flagged `orbit` as unmarked when its marker was
            # the first line of an 11-line header.
            j = i - 1
            while j >= 0 and (lines[j].startswith("#") or not lines[j].strip()):
                if not lines[j].strip() and j < i - 1:
                    break          # a blank line ends the block's own header
                j -= 1
            window = "\n".join(lines[j + 1:i])
            if not any(k in window for k in
                       ("DETERMINED", "PROVISIONAL", "UNDETERMINED", "TRANSITIVE",
                        "CANDIDATE")):
                unmarked.append(line.split(":")[0])
        rep.add(FAIL if unmarked else OK, "config blocks declare their status",
                f"unmarked: {', '.join(unmarked)}" if unmarked else
                "every top-level block is marked")
    except Exception as exc:
        rep.add(WARN, "config blocks declare their status", f"not checked: {exc}")

    # -- generated biosphere inputs vs the config they were derived from -----
    #
    # biosphere/generated/ holds a C++ header, a driver and a PFT file, each
    # DERIVED from config/planet.yaml and each recording the config_sha256 it was
    # built against. That design is right and it is why the year length is not
    # written down anywhere by hand. But nothing verified the derivation had
    # actually been re-run, so all three drifted -- and to three DIFFERENT config
    # states, which is worse than one stale artifact because they disagree with
    # each other as well as with the config.
    #
    # The year length sizes arrays in a compiled header, so it cannot be checked
    # at runtime. This is the check.
    try:
        gen = ROOT / "biosphere" / "generated"
        cur = sha256_of(ROOT / "config" / "planet.yaml")
        stale = []
        for prov in sorted(gen.glob("*_provenance.json")):
            rec = json.loads(prov.read_text(encoding="utf-8"))
            was = rec.get("config_sha256")
            if was != cur:
                yl = rec.get("year_length_days")
                # Each artifact names its own generator, so the remediation is
                # read off the artifact rather than guessed. It used to name one
                # script for all three, which was wrong for two of them.
                who = rec.get("generator") or "(generator not recorded)"
                stale.append(f"{prov.name}"
                             + (f" (year_length_days={yl})" if yl else "")
                             + f" -> rerun {who}")
        if not list(gen.glob("*_provenance.json")):
            rep.add(WARN, "generated biosphere inputs", "none present")
        else:
            rep.add(FAIL if stale else OK, "generated biosphere inputs current",
                    "derived from an older config: " + "; ".join(stale)
                    if stale else "all derived from the current config")
    except Exception as exc:
        rep.add(WARN, "generated biosphere inputs", f"not checked: {exc}")

    # -- binaries vs the patches they should contain -------------------------
    #
    # ExoPlaSim builds one executable per (resolution, layers, ranks) triple, so
    # patching the source and running rebuilds only the configuration in use and
    # leaves the rest silently stale. That is failure class 11, and it fired
    # three times in one day. .venv is untracked and reinstallable, and a
    # reinstall discards every patch without warning, so this is also the check
    # that says a reinstall has happened.
    try:
        manifest = ROOT / "exoplasim" / "patches" / "binary_manifest.json"
        run_dir = (ROOT / ".venv" / "lib" / "python3.12" / "site-packages"
                   / "exoplasim" / "plasim" / "run")
        on_disk = sorted(run_dir.glob("most_plasim_*.x")) if run_dir.is_dir() else []
        if not manifest.is_file():
            rep.add(FAIL, "binary manifest",
                    "absent; run exoplasim/scripts/rebuild_binaries.py")
        elif not on_disk:
            rep.add(WARN, "binaries", "none built")
        else:
            mf = json.loads(manifest.read_text(encoding="utf-8"))
            known = mf.get("binaries", {})
            src = (ROOT / ".venv" / "lib" / "python3.12" / "site-packages"
                   / "exoplasim" / "plasim" / "src")
            bad = []
            for exe in on_disk:
                rec = known.get(exe.name)
                if rec is None:
                    bad.append(f"{exe.name}: not in manifest")
                    continue
                if rec["sha256"] != sha256_of(exe):
                    bad.append(f"{exe.name}: sha differs from manifest")
                    continue
                for name, want in (rec.get("sources") or {}).items():
                    got = sha256_of(src / name) if (src / name).is_file() else None
                    if got != want:
                        bad.append(f"{exe.name}: built from an older {name}")
            rep.add(FAIL if bad else OK, "binaries carry current patches",
                    "; ".join(bad) if bad else
                    f"{len(on_disk)} executables match the manifest")
    except Exception as exc:
        rep.add(WARN, "binaries", f"not checked: {exc}")

    # -- the cycle executable can parse the cycle the config asks for ---------
    #
    # The star-cycle patch is deliberately NOT resident, so this binary is
    # outside the matrix above and nothing else would notice it going stale.
    # That is the same gap failure class 11 came through: exactly one of six
    # executables was newer than the patched source, and the run died in
    # radini_ with "Cannot match namelist object name".
    #
    # Checking the namelist names are actually IN the binary is stronger than
    # checking a timestamp, and it is the specific thing that fails: adding a
    # component to config/planet.yaml without rebuilding leaves a binary whose
    # namelist has no slot for it, and Fortran rejects the whole group.
    try:
        cyc = ROOT / "exoplasim" / "inputs" / "exoplasim_cycle_t42"
        m = config["model"]
        exe = cyc / (f"most_plasim_t{int(str(m['resolution']).lstrip('Tt'))}"
                     f"_l{int(m['layers'])}_p{int(m['ncpus'])}.x")
        components = (config.get("stellar_cycle") or {}).get("components") or {}
        slots = {"medium": "", "long": "2"}
        if not components:
            rep.add(WARN, "cycle executable", "no stellar_cycle.components")
        elif not exe.is_file():
            rep.add(FAIL, "cycle executable",
                    f"{exe.name} absent; run "
                    "exoplasim/scripts/build_star_cycle_exoplasim.sh")
        else:
            blob = exe.read_bytes()
            missing = [f"gsolamp{slots[n]}" for n in components
                       if n in slots
                       and f"gsolamp{slots[n]}".encode() not in blob]
            unslotted = sorted(set(components) - set(slots))
            patch = (ROOT / "exoplasim" / "patches"
                     / "exoplasim-3.4.2-star-cycle.patch")
            problems = []
            if unslotted:
                problems.append(f"no namelist slot for {unslotted}")
            if missing:
                problems.append(f"binary lacks {missing}")
            # Compare by CONTENT, not mtime. An mtime test called a correct
            # binary stale as soon as a git stash/pop rewrote the patch file
            # without changing a byte of it, and it would equally have missed a
            # patch edited in place with a preserved timestamp.
            mf = (ROOT / "exoplasim" / "patches" / "cycle_binary_manifest.json")
            if not mf.is_file():
                problems.append("no cycle_binary_manifest.json; rebuild to record one")
            else:
                rec = json.loads(mf.read_text(encoding="utf-8"))
                if rec.get("executable_sha256") != sha256_of(exe):
                    problems.append("executable differs from its build manifest")
                elif patch.is_file() and rec.get("patch_sha256") != sha256_of(patch):
                    problems.append("built from an older star-cycle patch")
            rep.add(FAIL if problems else OK,
                    "cycle executable matches the configured cycle",
                    "; ".join(problems) if problems else
                    f"{len(components)} components, all present in {exe.name}")
    except Exception as exc:
        rep.add(WARN, "cycle executable", f"not checked: {exc}")

    rep.show()
    return 1 if rep.failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
