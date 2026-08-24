"""Proof for the restart converter: every claim, with something that can fail.

Worldbuilding frame: fixtures for the Vesper climate model's saved state.
Nothing here is about the real world.

Run through `convert_restart.py --self-test`.

WHAT IS AND IS NOT PROVEN HERE. Nothing in this file runs the model, so
nothing here can say a converted state LOADS, settles, or approaches the same
climate as a cold target control. Those are the tests that need a target
executable and a target template it wrote, and they are tracked separately.
What is proven is everything upstream of that: the framing, the schema against
the model source, both transforms against identities that have right answers,
and the whole conversion against a fixture whose answer is known in advance.

EVERY POSITIVE TEST HAS A NEGATIVE CONTROL, because a transformation test that
cannot fail measures nothing. The controls are the specific mistakes: a prefix
copy of a packed spectral array, a longitude axis treated as though it did not
wrap, a record-length heuristic in place of a schema.
"""
from __future__ import annotations

import struct
import tempfile
from pathlib import Path

import numpy as np

import _paths
import convert_restart as cv
import reset_restart_accumulators as ra
import restart_format as rf
import restart_schema as rs
import restart_transforms as rt
import rungs

MODEL_SRC = _paths.MODEL_SRC / "plasim" / "src"

# Declared before any of it was run. A same-width conversion is exact, so the
# only tolerance a fixture needs is for an integral accumulated over cells.
INTEGRAL_RTOL = 1e-12


class Failed(Exception):
    pass


def _require(condition, what: str) -> None:
    if not condition:
        raise Failed(what)


def _refuses(fn, what: str, naming: str | None = None) -> None:
    """`fn` must raise. A test whose negative case passes silently is not one."""
    try:
        fn()
    except (cv.ConversionError, rf.RestartFormatError, ValueError) as exc:
        if naming and naming not in str(exc):
            raise Failed(f"{what}: refused without naming '{naming}': {exc}")
        return
    raise Failed(f"{what}: was accepted and should not have been")


def _next_rung(rung: str) -> str:
    """The rung above `rung` on the ladder, for a fixture that changes support.

    Taken from the registry rather than written as T42, so the fixture follows
    the ladder if it ever gains a rung between these two. SPAT-2.
    """
    ladder = [r for r, _ in sorted(rungs.RUNGS.items(), key=lambda kv: kv[1])]
    i = ladder.index(rung)
    if i + 1 >= len(ladder):
        raise Failed(f"{rung} is the top of the ladder; this fixture needs a "
                     "rung above the donor's")
    return ladder[i + 1]


def _donor() -> Path:
    """The one real restart in the tree, or nothing to test against."""
    runs = sorted(p for p in (_paths.RUNS).glob("run_*/plasim_restart"))
    if not runs:
        raise Failed("no run directory holds a plasim_restart to test against")
    return runs[0]


# ---------------------------------------------------------------------------
# A synthetic target template. FIXTURES ONLY.
# ---------------------------------------------------------------------------

def _synthetic_template(src: cv.RestartState, nlat: int, real_bytes: int,
                        path: Path) -> Path:
    """A template at another resolution, built from the schema rather than run.

    THIS IS NOT HOW A TEMPLATE IS MADE. A real one is written by the exact
    target executable, which is the only thing that knows its own record set,
    its accumulators' reset values and its compiler's seed shape. This builds
    one from the schema so the conversion machinery can be exercised without a
    model run, and it is deliberately not reachable from the command line.

    Its land-sea mask is the source's, remapped and binarised, so the fixture
    has a coastline that is consistent between the two resolutions rather than
    an invented one.
    """
    _, _, ntru = rungs.geometry(rungs.rung_of_latitudes(nlat))
    g = rs.Geometry(nlat=nlat, nlev=src.geometry.nlev,
                    nlsoil=src.geometry.nlsoil, nlev_oce=src.geometry.nlev_oce,
                    nesp=_round_up(ntru, src), nseedlen=src.geometry.nseedlen)
    inventory = rs.inventory_from_source(MODEL_SRC)
    dtype = cv.REAL[real_bytes]
    weights = rt.build_weights(src.geometry.nlat, nlat)
    land, _ = rt.remap(src.decode("dls"), weights)
    land = (land > 0.5).astype(np.float64)

    records = []
    for rec in src.records:
        name = rec.name
        pol = rs.POLICY[name]
        clean_zero = (pol.semantic == rs.ACCUMULATOR and pol.model_reset == "zero")
        if name == "nlat":
            payload = struct.pack("<i", g.nlat)
        elif name == "nlon":
            payload = struct.pack("<i", g.nlon)
        elif name == "nrsp":
            payload = struct.pack("<i", g.nrsp)
        elif rec.nbytes == 4 and inventory[name].writer == "put_restart_integer":
            # Counters arrive at the start of a window; nlev and nlsoil are
            # geometry the target must agree with.
            payload = b"\x00" * 4 if clean_zero else rec.payload
        elif name == "seed":
            payload = rec.payload
        elif rs.POLICY[name].action == rs.REQUIRE_EQUAL:
            # Configuration the target must agree with rather than re-derive.
            # A fixture that zeroed these would be refused, and rightly.
            payload = rec.payload
        else:
            counts = rs.candidate_counts(name, g, inventory)
            src_counts = rs.candidate_counts(name, src.geometry, inventory)
            which = src_counts.index(rec.nbytes // src.real_bytes)
            n = counts[which]
            if name == "dls":
                values = land
            elif pol.semantic == rs.STATIC_GRID and n % g.nugp == 0:
                # A recognisable value, so a fallback that reaches for the
                # template shows up in the report as coming from here.
                values = np.full(n, -7.0)
            elif pol.model_reset == "sentinel":
                # A template has to arrive CLEAN, and clean is not zero for a
                # running minimum. The converter refuses a dirty one, so a
                # fixture that zeroed these would be refused and rightly.
                values = np.full(n, pol.reset_value)
            else:
                values = np.zeros(n)
            payload = values.astype(dtype).tobytes()
        records.append(rf.Record(name=name, payload=payload, offset=-1))
    rf.write(path, records, overwrite=True)
    return path


def _round_up(ntru: int, src: cv.RestartState) -> int:
    """NESP for a target truncation, at the source's own thread count.

    A real template carries its own; a fixture has to pick one, and picking the
    donor's keeps the fixture's spectral accumulators the shape a run at the
    same thread count would give.
    """
    nrsp = (ntru + 1) * (ntru + 2)
    npro = next(c for c in (1, 2, 4, 8, 16, 32, 64)
                if -(-src.geometry.nrsp // c) * c == src.geometry.nesp)
    return -(-nrsp // npro) * npro


# ---------------------------------------------------------------------------
# The tests
# ---------------------------------------------------------------------------

def test_framing(tmp: Path, donor: Path) -> list[str]:
    said = []
    original = donor.read_bytes()
    records = rf.read(donor)
    _require(rf.encode(records) == original,
             "parser and writer do not round-trip the donor byte for byte")
    said.append(f"parser round-trips {len(records)} records byte for byte")

    _refuses(lambda: rf.read(_write(tmp / "truncated", original[:-40])),
             "a truncated file", naming="truncated")
    bad = bytearray(original)
    bad[4 + 16:4 + 16 + 4] = struct.pack("<i", 99)      # trailing marker of rec 1
    _refuses(lambda: rf.read(_write(tmp / "marker", bytes(bad))),
             "disagreeing record markers", naming="markers disagree")
    _refuses(lambda: rf.read(_write(tmp / "dup", original + original[:4 + 16 + 4] +
                                    original[4 + 16 + 4:4 + 16 + 4 + 12])),
             "a duplicate record name", naming="twice")
    name_only = original[:4 + 16 + 4]
    _refuses(lambda: rf.read(_write(tmp / "odd", name_only)),
             "a name with no payload", naming="ends mid-pair")
    said.append("refuses truncation, marker disagreement, duplicates and a "
                "half-written pair, each naming what it found")

    # The positional rule: a payload that happens to be printable is a payload.
    printable = rf.Record(name="zsolars", payload=b"ABCDEFGHIJKLMNOP", offset=-1)
    header = rf.Record(name="nstep", payload=struct.pack("<i", 7), offset=-1)
    rf.write(tmp / "printable", [header, printable], overwrite=True)
    back = rf.read(tmp / "printable")
    _require([r.name for r in back] == ["nstep", "zsolars"],
             "a 16-byte printable PAYLOAD was read as a record name")
    _require(back[1].payload == b"ABCDEFGHIJKLMNOP",
             "the printable payload did not survive")
    said.append("a 16-byte printable payload stays a payload: names are found "
                "by position, not by looking like text")
    return said


def _write(path: Path, data: bytes) -> Path:
    path.write_bytes(data)
    return path


def test_schema(donor: Path) -> list[str]:
    said = []
    gaps = rs.check_policy_covers_source(MODEL_SRC)
    _require(not gaps, "the schema does not cover the model source: " + "; ".join(gaps))
    inventory = rs.inventory_from_source(MODEL_SRC)
    said.append(f"policy covers all {len(inventory)} names the build can write, "
                "with none orphaned")

    state = cv.load(donor)
    for rec in state.records:
        if rs.POLICY[rec.name].action == rs.SEED:
            want = [state.geometry.nseedlen * 4]
        elif inventory[rec.name].writer == "put_restart_integer":
            want = [4]
        else:
            want = [n * state.real_bytes
                    for n in rs.candidate_counts(rec.name, state.geometry, inventory)]
        _require(rec.nbytes in want,
                 f"'{rec.name}' is {rec.nbytes} bytes and the schema predicts {want}")
    said.append(f"the schema predicts every one of {len(state.records)} record "
                f"lengths in a {state.geometry.label} restart")

    # Negative control: the length heuristic the schema replaces cannot tell an
    # eight-byte real array from an integer array of twice the count.
    ints = np.arange(state.geometry.nugp * 2, dtype="<i4").tobytes()
    reals = np.zeros(state.geometry.nugp, dtype="<f8").tobytes()
    _require(len(ints) == len(reals),
             "the length-heuristic control is not set up: the two payloads "
             "should be the same size")
    said.append("negative control: an int32 array and a real8 array of the "
                "same byte length are indistinguishable without the schema")
    return said


def test_spectral() -> list[str]:
    said = []
    rng = np.random.default_rng(20260824)
    v21 = rng.standard_normal(2 * len(rt.triangular_modes(21)))
    up, lost = rt.project_spectral(v21, 21, 42)
    _require(lost == (0.0, 0.0), "an increase in truncation discarded something")
    modes21, modes42 = rt.triangular_modes(21), rt.triangular_modes(42)
    at42 = {m: k for k, m in enumerate(modes42)}
    for k, mode in enumerate(modes21):
        j = at42[mode]
        _require(np.array_equal(up[2 * j:2 * j + 2], v21[2 * k:2 * k + 2]),
                 f"mode {mode} did not survive the increase in truncation")
    new = [k for k, m in enumerate(modes42) if m not in set(modes21)]
    _require(all(up[2 * k] == 0.0 and up[2 * k + 1] == 0.0 for k in new),
             "a mode the source never resolved came up nonzero")
    said.append("T21 to T42 preserves every shared mode and zeroes all "
                f"{len(new)} new ones")

    back, _ = rt.project_spectral(up, 42, 21)
    _require(np.array_equal(back, v21), "T21 to T42 to T21 is not the identity")
    said.append("T21 to T42 to T21 restores every coefficient bit for bit")

    v42 = rng.standard_normal(2 * len(modes42))
    down, (l2, mx) = rt.project_spectral(v42, 42, 21)
    kept = set(modes21)
    lost_values = np.concatenate([v42[2 * k:2 * k + 2]
                                  for k, m in enumerate(modes42) if m not in kept])
    _require(np.isclose(l2, np.linalg.norm(lost_values), rtol=0, atol=1e-12),
             "the reported discarded L2 does not match an independent sum")
    _require(np.isclose(mx, np.abs(lost_values).max()),
             "the reported discarded maximum does not match")
    said.append(f"T42 to T21 reports the discarded norm exactly ({l2:.4f} over "
                f"{len(lost_values) // 2} modes)")

    # Negative control: the prefix copy the specification warns against. It
    # agrees for the first m block and diverges after it, so a test that could
    # not tell them apart would not be testing the packing at all.
    prefix = np.zeros_like(up)
    prefix[:v21.size] = v21
    _require(not np.array_equal(prefix, up),
             "a prefix copy and a mode-wise projection agree, so this test "
             "cannot detect the packing mistake it exists for")
    differing = np.count_nonzero(prefix != up)
    said.append(f"negative control: a prefix copy differs from the projection "
                f"in {differing} coefficients and would be caught")
    return said


def test_grid() -> list[str]:
    said = []
    rng = np.random.default_rng(9)
    for nlat_s, nlat_t in ((32, 64), (64, 32), (32, 128), (128, 64)):
        w = rt.build_weights(nlat_s, nlat_t)
        area_s, area_t = rt.cell_area(nlat_s), rt.cell_area(nlat_t)
        flat, missed = rt.remap(np.full(w.nugp_src, 3.5), w)
        _require(not missed.any(), f"{nlat_s} to {nlat_t} left a target cell unfilled")
        _require(np.allclose(flat, 3.5, rtol=0, atol=1e-13),
                 f"{nlat_s} to {nlat_t} did not preserve a constant")
        f = rng.standard_normal(w.nugp_src)
        out, _ = rt.remap(f, w)
        before, after = float((f * area_s).sum()), float((out * area_t).sum())
        _require(abs(after - before) <= INTEGRAL_RTOL * abs(before),
                 f"{nlat_s} to {nlat_t} did not conserve the global integral: "
                 f"{before} became {after}")
    said.append("constants and global integrals survive every direction "
                "between nlat 32, 64 and 128, to 1e-12 relative")

    w = rt.build_weights(32, 32)
    f = rng.standard_normal(w.nugp_src)
    _require(np.array_equal(rt.remap(f, w)[0], f),
             "a same-resolution remap is not the identity")
    said.append("a same-resolution remap is bit-identical, so a conversion "
                "that only changes precision cannot move the state")

    mask = np.zeros(w.nugp_src, dtype=bool)
    mask[:1024] = True
    out, missed = rt.remap(np.ones(w.nugp_src), w, mask=mask)
    _require(missed.any(), "a masked remap filled cells with no eligible source")
    _require(np.allclose(out[~missed], 1.0),
             "a masked remap leaked values from outside the class")
    said.append(f"a masked remap fills only from its own class and names the "
                f"{int(missed.sum())} cells it could not fill")

    # Negative control: longitude without its wrap. The seam cells are the ones
    # that lose their partner, so conservation fails where the tiling is gone.
    w2 = rt.build_weights(32, 64)
    keep = w2.src % 64 != 63
    maimed = rt.RemapWeights(tgt=w2.tgt[keep], src=w2.src[keep],
                             area=w2.area[keep], nugp_src=w2.nugp_src,
                             nugp_tgt=w2.nugp_tgt, shape_src=w2.shape_src,
                             shape_tgt=w2.shape_tgt)
    g = rng.standard_normal(w2.nugp_src)
    broken, _ = rt.remap(g, maimed)
    a_s, a_t = rt.cell_area(32), rt.cell_area(64)
    residual = abs(float((np.nan_to_num(broken) * a_t).sum())
                   - float((g * a_s).sum())) / abs(float((g * a_s).sum()))
    _require(residual > INTEGRAL_RTOL,
             "dropping a longitude column still conserved the integral, so "
             "this test cannot detect a broken longitude axis")
    said.append(f"negative control: dropping one longitude column breaks "
                f"conservation by {residual:.2e} relative and would be caught")
    return said


def test_precision(tmp: Path, donor: Path) -> list[str]:
    said = []
    state = cv.load(donor)
    g = state.geometry
    _require(state.real_bytes == 8, "the donor is not an eight-byte restart")
    _require(cv.infer_real_bytes(state.by_name, g.nlat, g.nlev, g.nrsp) == 8,
             "precision inference disagrees with the file it just read")
    said.append("precision inference agrees across sp, sz and dls, from "
                "element counts the four-byte headers already pin down")

    # Negative control: make one invariant record the wrong width and the
    # inference must refuse rather than take a majority.
    faked = dict(state.by_name)
    faked["dls"] = rf.Record(name="dls", offset=-1,
                             payload=np.zeros(g.nugp, dtype="<f4").tobytes())
    _refuses(lambda: cv.infer_real_bytes(faked, g.nlat, g.nlev, g.nrsp),
             "invariant records disagreeing about the real width",
             naming="disagree")
    said.append("negative control: one record of the other width is refused, "
                "not out-voted")

    values = state.decode("st")
    narrowed = values.astype("<f4").astype(np.float64)
    widened_back = narrowed.astype("<f4").astype(np.float64)
    _require(np.array_equal(narrowed, widened_back),
             "four to eight bytes is not exact")
    worst = float(np.abs(narrowed - values).max())
    _require(worst > 0.0,
             "narrowing this field lost nothing, so the cast-error report "
             "cannot be checked against it")
    said.append(f"four to eight bytes is exact; eight to four costs at most "
                f"{worst:.3e} on the temperature field and is reported")
    return said


def test_end_to_end(tmp: Path, donor: Path) -> list[str]:
    said = []
    inventory = rs.inventory_from_source(MODEL_SRC)
    src = cv.load(donor)

    # 1. The donor onto a template cut from itself: same grid, same precision.
    # The template has to be CLEAN, which the donor is not -- a run does not end
    # on an output boundary -- so it is normalised first, exactly as
    # build_restart_template.py does it. The right answer is then the cleaned
    # donor, byte for byte: nothing but the accumulation window may move.
    # Named from the donor's own geometry rather than from a literal, so the
    # fixture cannot be filed under a rung it is not.
    same_rung = tmp / f"template_{src.geometry.label.lower()}"
    ra.reset(donor, same_rung)
    clean = cv.load(same_rung)
    cv.check_compatible(src, clean, inventory)
    records, reports = cv.convert(src, clean)
    rf.write(tmp / "identity", records, overwrite=True)
    _require((tmp / "identity").read_bytes() == same_rung.read_bytes(),
             "a conversion onto the donor's own grid and precision moved "
             "something other than the accumulation window")
    moved = sum(1 for a, b in zip(rf.read(tmp / "identity"), src.records)
                if a.payload != b.payload)
    _require(moved > 0, "the cleaned template is identical to the donor, so "
                        "this test cannot tell a conversion from a copy")
    said.append(f"the donor converted onto a template cut from itself "
                f"reproduces it byte for byte across {len(records)} records, "
                f"with only the {moved} accumulation records moved")

    # 2. Up one rung and back down. The spectral state has a right answer here:
    # projection keeps every shared (m,n) exactly, so the round trip is the
    # identity on the prognostic spectral records however far apart the rungs
    # are. The rung above is taken from the registry rather than fixed at T42,
    # which also makes this a non-doubling ratio and works the remap harder.
    # One rung up, named from the fixture's own geometry.
    finer = rungs.RUNGS[_next_rung(src.geometry.label)]
    up_path = tmp / f"template_{_next_rung(src.geometry.label).lower()}"
    up_state = cv.load(_synthetic_template(src, finer, 8, up_path))
    cv.check_compatible(src, up_state, inventory)
    recs42, rep42 = cv.convert(src, up_state)
    stepped = tmp / f"stepped_{_next_rung(src.geometry.label).lower()}"
    rf.write(stepped, recs42, overwrite=True)
    mid = cv.load(stepped)
    _require(mid.geometry.label == _next_rung(src.geometry.label),
             "the fixture did not reach the next rung up")

    # Back down onto the same clean T21 template. A dirty one is refused, and
    # the donor itself is dirty: a run does not end on an output boundary.
    cv.check_compatible(mid, clean, inventory)
    recs21, rep21 = cv.convert(mid, clean)
    returned = tmp / f"returned_{src.geometry.label.lower()}"
    rf.write(returned, recs21, overwrite=True)
    back = cv.load(returned)
    spectral = [n for n, p in rs.POLICY.items()
                if p.action == rs.PROJECT and n in src.by_name]
    for name in spectral:
        _require(back.by_name[name].payload == src.by_name[name].payload,
                 f"'{name}' did not survive T21 to T42 to T21 unchanged")
    said.append(f"{src.geometry.label} to {mid.geometry.label} and back "
                f"returns all {len(spectral)} spectral prognostic records "
                "bit for bit")

    # 3. What the report has to carry for the fields that cannot round-trip.
    remapped = [r for r in rep42 if r.action == rs.REMAP]
    _require(remapped, "no record was remapped going up a resolution")
    reservoirs = [r for r in remapped if "inventory" in r.detail]
    _require(reservoirs, "no reservoir reported its inventory")
    # Conservation through the WHOLE pipeline, not just the operator: a field
    # remapped over the sphere with no mask to interrupt it must arrive at T42
    # with the same integral it left T21 with, cast included. Every reservoir
    # in this configuration is masked, so the unmasked case is checked on a
    # field that carries no `conserve` label rather than left to an empty loop.
    unmasked = [n for n, pol in rs.POLICY.items()
                if pol.action == rs.REMAP and pol.domain == "global"
                and n in src.by_name]
    _require(unmasked, "no field is remapped over the whole sphere, so "
                       "end-to-end conservation cannot be checked here")
    area21, area42 = rt.cell_area(src.geometry.nlat), rt.cell_area(mid.geometry.nlat)
    for name in unmasked:
        before = src.decode(name).reshape(-1, src.geometry.nugp)
        after = mid.decode(name).reshape(-1, mid.geometry.nugp)
        for k, (a, b) in enumerate(zip(before, after)):
            total = float((a * area21).sum())
            if abs(total) < 1e-30:
                continue
            moved = abs(float((b * area42).sum()) - total) / abs(total)
            _require(moved <= INTEGRAL_RTOL,
                     f"'{name}' level {k} is remapped over the whole sphere "
                     f"and its integral moved by {moved:.2e}")
    said.append(f"{len(unmasked)} unmasked fields keep their global integral "
                "through the whole conversion, to 1e-12 relative")
    masked = [r for r in reservoirs if rs.POLICY[r.name].domain != "global"]
    for r in masked:
        _require("fallback_cells" in r.detail,
                 f"'{r.name}' is remapped under a mask and does not report "
                 "how many target cells it could not fill")
    said.append(f"{len(reservoirs)} reservoirs report their inventory; the "
                f"{len(masked)} remapped under a mask also report the cells "
                "the mask left for the template to fill")

    static = [n for n, pol in rs.POLICY.items()
              if pol.action == rs.TARGET and n in src.by_name]
    for name in static:
        _require(mid.by_name[name].payload == up_state.by_name[name].payload,
                 f"'{name}' is the target's own and did not arrive from the "
                 "template byte for byte")
    said.append(f"all {len(static)} static and target-owned records arrive "
                "from the template byte for byte, none of them the donor's")

    # Eight bytes down to four and back. The only thing that may have moved is
    # the narrowing, so the round trip is checked against the cast error the
    # report declared rather than against a tolerance invented here.
    narrow = _synthetic_template(src, src.geometry.nlat, 4, tmp / "template_fp32")  # noqa: E501
    fp32_state = cv.load(narrow)
    _require(fp32_state.real_bytes == 4, "the four-byte fixture is not four-byte")
    recs32, rep32 = cv.convert(src, fp32_state)
    rf.write(tmp / "fp32", recs32, overwrite=True)
    wide = cv.load(tmp / "fp32")
    recs_back, _ = cv.convert(wide, clean)
    rf.write(tmp / "fp64_again", recs_back, overwrite=True)
    widened = cv.load(tmp / "fp64_again")
    declared = {r.name: r.detail.get("cast_abs_error", 0.0) for r in rep32}
    moved = 0
    for name, pol in rs.POLICY.items():
        if pol.action not in (rs.PROJECT, rs.REMAP) or name not in src.by_name:
            continue
        before, after = src.decode(name), widened.decode(name)
        worst = float(np.abs(after - before).max())
        _require(worst <= declared[name] * (1.0 + 1e-12),
                 f"'{name}' moved by {worst:.3e} through eight to four to "
                 f"eight bytes and the report declared {declared[name]:.3e}")
        moved += worst > 0.0
    said.append(f"eight to four to eight bytes moves {moved} records, none of "
                "them by more than the cast error the report declared")

    recompute = [r.name for r in rep42 if r.action == rs.RECOMPUTE]
    _require(sorted(recompute) == sorted(
        n for n, p in rs.POLICY.items()
        if p.action == rs.RECOMPUTE and n in src.by_name),
        "the report does not name every derived record the model must rebuild")
    said.append("the report names every derived record the target model has "
                f"to rebuild before the state is self-consistent: {', '.join(sorted(recompute))}")

    # 4. Refusals, each before an output exists.
    fabricated = list(src.records) + [rf.Record(name="notarecord",
                                                payload=b"\x00" * 8, offset=-1)]
    rf.write(tmp / "unknown", fabricated, overwrite=True)
    _refuses(lambda: cv.check_compatible(cv.load(tmp / "unknown"), src, inventory),
             "a record with no schema entry", naming="notarecord")
    thinner = [r for r in src.records if r.name != "dust3"]
    rf.write(tmp / "thinner", thinner, overwrite=True)
    _refuses(lambda: cv.check_compatible(src, cv.load(tmp / "thinner"), inventory),
             "a template with a different record set", naming="dust3")
    _refuses(lambda: rf.write(tmp / "identity", records),
             "overwriting an existing output without being told to",
             naming="exists")
    # The donor is itself the control here: a run does not end on an output
    # boundary, so its own restart is a template with a partial window in it.
    _refuses(lambda: cv.check_compatible(src, src, inventory),
             "a template whose accumulation window is partial",
             naming="not at the value the model resets them to")
    said.append("refuses an unknown record, a template whose record set "
                "differs, a template cut from a run mid-window, and an "
                "overwrite it was not asked for")
    return said


def run() -> int:
    donor = _donor()
    sections = []
    with tempfile.TemporaryDirectory(prefix="convert_restart_selftest_") as d:
        tmp = Path(d)
        for title, fn in (("framing", lambda: test_framing(tmp, donor)),
                          ("schema", lambda: test_schema(donor)),
                          ("spectral projection", test_spectral),
                          ("Gaussian remap", test_grid),
                          ("precision", lambda: test_precision(tmp, donor)),
                          ("end to end", lambda: test_end_to_end(tmp, donor))):
            try:
                sections.append((title, fn(), None))
            except Failed as exc:
                sections.append((title, [], str(exc)))

    failed = [t for t, _, err in sections if err]
    for title, lines, err in sections:
        print(f"{title}:")
        for line in lines:
            print(f"  ok   {line}")
        if err:
            print(f"  FAIL {err}")
    print()
    if failed:
        print(f"self-test FAILED in {len(failed)} of {len(sections)} sections: "
              + ", ".join(failed))
        return 1
    print(f"self-test passed, {len(sections)} sections, against {donor}")
    return 0
