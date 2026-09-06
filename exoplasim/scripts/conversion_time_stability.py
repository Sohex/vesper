#!/usr/bin/env python3
"""What timestep the NCONVTIME term is stable at, and it is not the guard's.

    python exoplasim/scripts/conversion_time_stability.py
    python exoplasim/scripts/conversion_time_stability.py --rung T21 --dt 15
    python exoplasim/scripts/conversion_time_stability.py --against <run>/plasim_diag

Worldbuilding frame: a numerical experiment on the Vesper project's climate
model. Nothing here is about the simulated planet.

WHAT `NCONVTIME` DOES AND WHY ITS TIMESTEP IS ITS OWN QUESTION.
`model.conversion_time_level` takes the adiabatic reference conversion's
divergence half off `sdt` -- the centred mean of t-dt and t+dt that
`spectrala`'s semi-implicit solve produces -- and puts it on the divergence at
t, so the two halves of the conversion meet at one time level. It is world-0ov's
first repair route and is refuted; `exoplasim/notes/resolution-tuned-parameters.md`
carries the refutation. It remains as a CONTROL.

WHY THE GUARD IT CARRIED CANNOT BE THE BOUNDARY. The model refuses the setting
above `a / (c sqrt(N(N+1)))` with `c = sqrt(R T0 / (1 - kappa))`, the explicit
gravity-wave limit. That speed is right -- taken on the model's own fastest
external mode instead, the largest eigenvalue of the semi-implicit vertical
structure matrix, the limit moves by 1% at every rung -- and it is still not
where the term becomes usable, because the mode the term destabilises is not a
gravity wave. `sdt - sd` is the second time difference:
O(dt^2) for a smooth mode and, for the LEAPFROG COMPUTATIONAL MODE, which
alternates sign every step, exactly `-2 sd`. The term therefore feeds the
computational mode, the only thing that damps that mode is the Robert-Asselin
filter, and the boundary is a function of PNU. At PNU = 0 there is no stable
timestep at all, which no gravity-wave CFL can express. world-bt3b.

WHAT THIS COMPUTES. The one-step amplification of the model's own linearised
adiabatic step: `spectrala`'s divergence solve with the nonlinear tendencies
zero, step 3a's term, the leapfrog advance, the time filter in the two halves
the model splits it into, and `spectrald`'s implicit hyperdiffusion, Rayleigh
drag and Newtonian cooling. Every matrix is rebuilt from `config/planet.yaml`
through the same arithmetic `initpm`, `initsi` and `makebm` use, so this is the
model's map and not a model of it.

TWO ANSWERS, AND THEY ARE FOR DIFFERENT QUESTIONS. The BOUNDARY reported here
is the exact spectral radius of the one-step map, solved rather than iterated,
which is what a caller needs before buying a run. The MODEL cannot solve an
eigenvalue problem at startup, so its guard iterates -- and a finite iteration
on a neutral map is biased high by about 1.8e-4 per step, because the neutral
modes sit exactly on the unit circle and the map is not normal. No threshold
can be set under that bias, so the guard runs a CONTROL ARM instead: the same
iteration with the term off, whose right answer is one, and whose measured
distance from one is the instrument's own error on that configuration. It
refuses when the effect exceeds that error. `verdict` reproduces exactly that,
and `--against` compares it with what the model printed, arm for arm, so a
disagreement is a defect in one implementation rather than a question about
either.

VALIDATION AGAINST THE MODEL'S OWN BEHAVIOUR is in
`exoplasim/notes/convdecomp-reproducibility.md`, which carries the T21 arms and
the criterion they were judged against.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
from _paths import ANALYSIS, CONFIG, PROJECT_ROOT  # noqa: E402
import lapse  # noqa: E402  from lib/
import rungs  # noqa: E402  the one rung-to-dimension mapping

OUT = ANALYSIS / "conversion_time_stability.json"
# The sweep writes BESIDE the ladder table, never over it: one is every rung
# at its route step and the other is one rung at many steps, and a path shared
# between them is one --sweep away from deleting the ladder. The RUNG IS A
# FIELD AND NOT IN THE NAME, because a resolution literal in an artifact path
# is a rung baked into a filename; smoke_test refuses one.
OUT_SWEEP = ANALYSIS / "conversion_time_stability_sweep.json"

# The model's own reference day for the write interval and the timestep count,
# `plasim.f90`'s `day_24hr`. Not the planet's rotation.
SECONDS_PER_ABSOLUTE_DAY = 86400.0
# `p_earth.f90:82` declares PNU and `config/planet.yaml` may override it with
# `model.robert_filter`. Absent there, this is what the model runs.
PNU_DEFAULT = 0.1
# The measured span, matching `conversion_time_amplification` in the model.
DISCARD, MEASURE = 1000, 3000
# The margin an EXACT spectral radius has to clear to be called growth. This is
# for `boundary_minutes` and the ladder table, which solve the eigenvalue
# problem rather than iterating it, so the only error is double-precision
# roundoff on a 42x42 solve. The model's own guard cannot use a number like
# this and does not: it iterates, and a finite iteration on a neutral map
# carries a bias of about 1.8e-4 per step that no threshold can be set under.
# `verdict` below is what the model does instead.
NEUTRAL = 1.0e-10


class Column:
    """The model's vertical structure at one configuration, built its way."""

    def __init__(self, cfg: dict, rung: str, pnu: float):
        self.rung = rung
        self.ntru = int(rung.lstrip("Tt"))
        self.nlev = int(cfg["model"]["layers"])
        self.pnu = pnu
        gas_constant, cp = lapse.gas_properties(cfg)
        self.akap = gas_constant / cp
        self.gascon = gas_constant
        self.t0_k = float(cfg["model"]["semi_implicit_reference_temperature_k"])
        self.plarad = float(cfg["planet"]["radius_earth"]) * 6371000.0
        self.sidereal = float(cfg["planet"]["rotation_hours"]) * 3600.0
        self.ww = 2.0 * math.pi / self.sidereal
        ptop = float(cfg["model"]["model_top_hpa"]) * 100.0
        psurf = 100000.0
        if int(cfg["model"]["vertical_grid"]) != 4:
            raise SystemExit(
                f"model.vertical_grid is {cfg['model']['vertical_grid']} and "
                "this builds initpm's neqsig == 4 quartic only. The sigma grid "
                "is what tau and g are built from, so another one is another "
                "column and this would silently answer for the wrong model.")
        n = self.nlev
        sigmah = np.array([0.75 * (j / n) + 1.75 * (j / n) ** 3
                           - 1.5 * (j / n) ** 4 for j in range(1, n + 1)])
        top = ptop / psurf
        sigmah = sigmah - sigmah[0]
        sigmah = sigmah / sigmah[-1]
        sigmah = sigmah * (1.0 - top) + top
        self.sigmah = sigmah
        self.dsigma = np.empty(n)
        self.dsigma[0] = sigmah[0]
        self.dsigma[1:] = sigmah[1:] - sigmah[:-1]
        rdsig = 0.5 / self.dsigma
        cv = self.plarad * self.ww
        self.ct = cv * cv / gas_constant
        self.t0 = np.full(n, self.t0_k / self.ct)
        tkp = self.akap * self.t0
        t01s2 = np.zeros(n)
        t01s2[:n - 1] = self.t0[1:] - self.t0[:-1]
        zalp = np.zeros(n)
        zalp[1:] = np.log(sigmah[1:]) - np.log(sigmah[:-1])
        g = np.zeros((n, n))
        g[0, 0] = 1.0
        for j in range(2, n + 1):
            g[j - 1, j - 1] = 1.0 - zalp[j - 1] * sigmah[j - 2] / self.dsigma[j - 1]
            g[j - 1, 0:j - 1] = zalp[j - 1]
        cm = np.zeros((n, n))
        for j in range(n):
            cm[j, :] = g[:, j] * (self.dsigma[j] / self.dsigma[:])
        tau = np.zeros((n, n))
        zt01s2, zsig = t01s2[0], sigmah[0]
        tau[0, 0] = 0.5 * zt01s2 * (zsig - 1.0) + tkp[0]
        tau[1:, 0] = 0.5 * zt01s2 * self.dsigma[1:]
        for jl in range(2, n + 1):
            zttm, zsigm = zt01s2, zsig
            zt01s2, zsig = t01s2[jl - 1], sigmah[jl - 1]
            for j2 in range(1, n + 1):
                ztm = 1.0 if j2 <= jl else 0.0
                ztmm = 1.0 if j2 < jl else 0.0
                ztau = zttm * (zsigm - ztmm)
                if jl < n:
                    ztau += zt01s2 * (zsig - ztm)
                ztau = ztau * rdsig[jl - 1] * self.dsigma[j2 - 1]
                if j2 <= jl:
                    ztau += tkp[jl - 1] * cm[j2 - 1, jl - 1]
                tau[j2 - 1, jl - 1] = ztau
        self.g, self.c, self.tau, self.tkp = g, cm, tau, tkp
        # The half nconvtime moves off sdt, `tkp(jlev) * c(jlev2,jlev)`.
        tkc = np.zeros((n, n))
        for jl in range(n):
            for j2 in range(jl + 1):
                tkc[j2, jl] = tkp[jl] * cm[j2, jl]
        self.tkc = tkc
        # The damping, converted as `initpm` converts it: the namelist value in
        # days becomes seconds through the 24-hour day, then a rate in the
        # model's own time unit through the sidereal day.
        hd = cfg["model"]["hyperdiffusion"]
        if rung not in hd["timescales_days"]:
            raise SystemExit(
                f"model.hyperdiffusion.timescales_days has no {rung}, so this "
                "rung has no declared damping and the map cannot be built for "
                "it. The damping is what bounds the growth.")
        table = hd["timescales_days"][rung]
        self.ndel = int(hd["order_alpha"])
        self.nhdiff = int(round(float(hd["cutoff_fraction"]) * self.ntru))
        rate = lambda days: self.sidereal / (2.0 * math.pi
                                             * float(days) * SECONDS_PER_ABSOLUTE_DAY)
        self.tdissd = rate(table["divergence"])
        self.tdisst = rate(table["temperature"])
        sponge = [float(x) for x in cfg["model"]["rayleigh_sponge_rotations"]]
        if len(sponge) != n:
            raise SystemExit(
                f"model.rayleigh_sponge_rotations has {len(sponge)} entries "
                f"and the model has {n} levels")
        self.tfrc = np.array([
            self.sidereal / (2.0 * math.pi * x * self.sidereal) if x > 0 else 0.0
            for x in sponge])

    def sak(self, jn: int) -> float:
        """`readnl`'s hyperdiffusion shape at total wavenumber `jn`."""
        if jn < self.nhdiff:
            return 0.0
        zakk = 1.0 / float(self.ntru - self.nhdiff) ** self.ndel
        return zakk * float(jn - self.nhdiff) ** self.ndel

    def gravity_wave_limit_minutes(self) -> float:
        """The limit `plasim.f90`'s NCONVTIME guard used to test against."""
        cgw = math.sqrt(self.gascon * self.t0_k / (1.0 - self.akap))
        return self.plarad / (cgw * math.sqrt(self.ntru * (self.ntru + 1.0))) / 60.0

    def structure_gravity_wave_limit_minutes(self) -> float:
        """The same limit on the model's OWN fastest external mode.

        A CHECK THAT CAN FAIL, and it is what says the guard's estimate was not
        the problem. `sqrt(R T0 / (1 - kappa))` is a single-level estimate of
        the external mode's phase speed; the model's actual gravity-wave speeds
        are `sqrt(cn * lambda)` over the eigenvalues of the vertical structure
        matrix `makebm` builds from `t0`, `dsigma`, `g` and `tau`, and the
        fastest is the largest of them. If this and the estimate disagreed by
        much, the guard would have been mis-derived rather than aimed at the
        wrong mode, and the repair would be a different one.
        """
        n = self.nlev
        m = np.zeros((n, n))
        for j1 in range(n):
            for j2 in range(n):
                m[j2, j1] = (self.t0[j1] * self.dsigma[j2]
                             + np.dot(self.g[:, j1], self.tau[j2, :]))
        lam = max(abs(np.linalg.eigvals(m)))
        cv = self.plarad * self.ww
        return self.plarad / (math.sqrt(lam) * cv
                              * math.sqrt(self.ntru * (self.ntru + 1.0))) / 60.0

    def growth(self, dt_minutes: float, jn: int, nconvtime: bool = True) -> float:
        """Amplification per step of the linearised step at wavenumber `jn`."""
        n = self.nlev
        delt = dt_minutes * 60.0 * self.ww
        delt2 = 2.0 * delt
        cn = jn * (jn + 1.0)
        mf = np.zeros((n, n))
        for j1 in range(n):
            for j2 in range(n):
                mf[j2, j1] = delt * delt * (self.t0[j1] * self.dsigma[j2]
                                            + np.dot(self.g[:, j1], self.tau[j2, :]))
        mf += np.eye(n) / cn
        bm1 = np.linalg.inv(mf)
        sak = self.sak(jn)
        fd = 1.0 / (1.0 + delt2 * (self.tdissd * sak + self.tfrc))
        ft = 1.0 / (1.0 + delt2 * (self.tdisst * sak))
        pnu21 = 1.0 - 2.0 * self.pnu
        # The one-step map as a matrix; its spectral radius is the answer, and
        # an eigenvalue solve says it exactly where the model's iteration
        # approaches it.
        size = 2 * (2 * n + 1)
        a = np.zeros((size, size))
        half = 2 * n + 1
        for k in range(size):
            e = np.zeros(size)
            e[k] = 1.0
            zd, zt, zp = e[0:n], e[n:2 * n], e[2 * n]
            adm, atm, apm = e[half:half + n], e[half + n:half + 2 * n], e[half + 2 * n]
            zz = adm / cn + delt * (self.g.T @ atm + self.t0 * apm)
            sdt = zz @ bm1
            spt = self.dsigma @ sdt
            stt = -(self.tau.T @ sdt)
            if nconvtime:
                stt = stt + (self.tkc.T @ (sdt - zd))
            sdm = pnu21 * zd + self.pnu * adm
            stm = pnu21 * zt + self.pnu * atm
            spm = pnu21 * zp + self.pnu * apm
            dp = (2.0 * sdt - adm) * fd
            tp = (delt2 * stt + atm) * ft
            pp = apm - delt2 * spt
            sdm = sdm + self.pnu * dp
            stm = stm + self.pnu * tp
            spm = spm + self.pnu * pp
            out = np.zeros(size)
            out[0:n], out[n:2 * n], out[2 * n] = dp, tp, pp
            out[half:half + n] = sdm
            out[half + n:half + 2 * n] = stm
            out[half + 2 * n] = spm
            a[:, k] = out
        return float(max(abs(np.linalg.eigvals(a))))

    def iterated_growth(self, dt_minutes: float, jn: int,
                        nconvtime: bool = True) -> float:
        """The model's OWN estimator: the linearised step, iterated.

        `plasim.f90:conversion_time_amplification` in Python, down to the start
        vector and the span, so `--against` compares two implementations of one
        procedure rather than an iteration against an eigenvalue solve. The two
        differ on a NEUTRAL configuration -- the iteration is biased high by
        about 1.8e-4 per step, because the neutral modes sit exactly on the unit
        circle and the map is not normal -- and that difference is the whole
        reason the model runs a control arm instead of testing a threshold.
        """
        n = self.nlev
        delt = dt_minutes * 60.0 * self.ww
        delt2 = 2.0 * delt
        cn = jn * (jn + 1.0)
        mf = np.zeros((n, n))
        for j1 in range(n):
            for j2 in range(n):
                mf[j2, j1] = delt * delt * (self.t0[j1] * self.dsigma[j2]
                                            + np.dot(self.g[:, j1], self.tau[j2, :]))
        mf += np.eye(n) / cn
        bm1 = np.linalg.inv(mf)
        sak = self.sak(jn)
        fd = 1.0 / (1.0 + delt2 * (self.tdissd * sak + self.tfrc))
        ft = 1.0 / (1.0 + delt2 * (self.tdisst * sak))
        pnu21, pnu = 1.0 - 2.0 * self.pnu, self.pnu
        jl = np.arange(1, n + 1, dtype=float)
        zd, zt = np.sin(1.0 * jl), np.cos(2.0 * jl)
        zdm, ztm = np.cos(0.7 * jl), np.sin(1.3 * jl)
        zp, zpm, zlog = 0.25, -0.125, 0.0
        for step in range(1, DISCARD + MEASURE + 1):
            zz = zdm / cn + delt * (self.g.T @ ztm + self.t0 * zpm)
            sdt = zz @ bm1
            spt = self.dsigma @ sdt
            stt = -(self.tau.T @ sdt)
            if nconvtime:
                stt = stt + (self.tkc.T @ (sdt - zd))
            zdmn = pnu21 * zd + pnu * zdm
            ztmn = pnu21 * zt + pnu * ztm
            zpmn = pnu21 * zp + pnu * zpm
            zdn = (2.0 * sdt - zdm) * fd
            ztn = (delt2 * stt + ztm) * ft
            zpn = zpm - delt2 * spt
            zdmn = zdmn + pnu * zdn
            ztmn = ztmn + pnu * ztn
            zpmn = zpmn + pnu * zpn
            zd, zt, zp, zdm, ztm, zpm = zdn, ztn, zpn, zdmn, ztmn, zpmn
            norm = math.sqrt(zd @ zd + zt @ zt + zdm @ zdm + ztm @ ztm
                             + zp * zp + zpm * zpm)
            if step > DISCARD:
                zlog += math.log(norm)
            zd, zt, zdm, ztm = zd / norm, zt / norm, zdm / norm, ztm / norm
            zp, zpm = zp / norm, zpm / norm
        return math.exp(zlog / MEASURE)

    def verdict(self, dt_minutes: float) -> dict:
        """The model's own two-arm test, recomputed here.

        The measured arm against a CONTROL whose right answer is known -- the
        unmodified semi-implicit scheme is neutrally stable, which is what every
        production run rests on -- so the control's distance from one is what
        this instrument can resolve on this configuration, measured rather than
        assumed. An effect smaller than that is noise however tidy it looks.
        """
        on, jn_on, off, jn_off = 0.0, 0, 0.0, 0
        for jn in range(1, self.ntru + 1):
            g = self.iterated_growth(dt_minutes, jn, True)
            if g > on:
                on, jn_on = g, jn
            g = self.iterated_growth(dt_minutes, jn, False)
            if g > off:
                off, jn_off = g, jn
        return {"on": on, "on_wavenumber": jn_on, "off": off,
                "off_wavenumber": jn_off, "effect": on - off,
                "resolution": abs(off - 1.0), "refused": (on - off) > abs(off - 1.0)}

    def worst_growth(self, dt_minutes: float, nconvtime: bool = True):
        """`(growth, wavenumber)` over every total wavenumber the rung resolves."""
        best, at = 0.0, 0
        for jn in range(1, self.ntru + 1):
            g = self.growth(dt_minutes, jn, nconvtime)
            if g > best:
                best, at = g, jn
        return best, at

    def boundary_minutes(self, lo: float = 0.05, hi: float = 240.0):
        """The coarsest step whose worst mode does not grow, or None.

        None means every step down to `lo` grows, which is what PNU = 0 gives:
        the modification has no stable timestep at all there, and that is a
        statement about the scheme and not about how small a step one can
        afford.
        """
        if self.worst_growth(lo)[0] > 1.0 + NEUTRAL:
            return None
        for _ in range(44):
            mid = 0.5 * (lo + hi)
            if self.worst_growth(mid)[0] <= 1.0 + NEUTRAL:
                lo = mid
            else:
                hi = mid
        return lo


def _real(raw: str) -> float:
    return float(raw.replace("D", "E").replace("d", "e"))


def model_amplification(diag: Path):
    """What `plasim.f90`'s guard printed: both arms, the effect and the step."""
    text = diag.read_text(encoding="latin-1", errors="replace")
    on = off = dt = None
    for m in re.finditer(r"NCONVTIME: amplification\s+([0-9.EeDd+-]+)\s+per "
                         r"step at total wavenumber\s+(\d+)", text):
        on = (_real(m.group(1)), int(m.group(2)))
    for m in re.finditer(r"NCONVTIME: control arm\s+([0-9.EeDd+-]+)\s+per "
                         r"step at total wavenumber\s+(\d+)", text):
        off = (_real(m.group(1)), int(m.group(2)))
    for m in re.finditer(r"NCONVTIME: effect .*?runs at\s+([0-9.]+)\s+min",
                         text, re.S):
        dt = float(m.group(1))
    if on is None or off is None or dt is None:
        return None
    return {"on": on[0], "on_wavenumber": on[1], "off": off[0],
            "off_wavenumber": off[1], "timestep_minutes": dt}


def provenance(cfg_path: Path) -> dict:
    try:
        commit = subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        commit = None
    src = PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
    digest = hashlib.sha256()
    for path in sorted(src.glob("*.f90")):
        digest.update(path.read_bytes())
    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(cfg_path.relative_to(PROJECT_ROOT)),
        "config_sha256": hashlib.sha256(cfg_path.read_bytes()).hexdigest(),
        "model_source_sha256": digest.hexdigest(),
        "commit": commit,
        "numpy": np.__version__,
        "discard_steps": DISCARD,
        "measured_steps": MEASURE,
        "neutral_tolerance": NEUTRAL,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--rung", default=None,
                        help="one rung; default is every rung the escalation "
                             "route visits plus every rung with a measured "
                             "stability ceiling")
    parser.add_argument("--dt", type=float, default=None,
                        help="one timestep in minutes; default is the route's "
                             "step for the rung")
    parser.add_argument("--pnu", type=float, default=None,
                        help="Robert-Asselin coefficient; default is "
                             "model.robert_filter, or the model's own 0.1")
    parser.add_argument("--against", type=Path, default=None,
                        help="a run's plasim_diag; compare the amplification "
                             "the model printed with the one computed here")
    parser.add_argument("--sweep", default=None,
                        help="comma-separated timesteps in minutes; with "
                             "--rung, record the growth at each of them so a "
                             "note can cite the artifact instead of restating "
                             "the numbers")
    parser.add_argument("--output", type=Path, default=None,
                        help="default is the ladder table, or the sweep table "
                             "beside it when --sweep is given")
    args = parser.parse_args()
    if args.output is None:
        args.output = OUT_SWEEP if args.sweep else OUT

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    pnu = args.pnu if args.pnu is not None else float(
        cfg["model"].get("robert_filter", PNU_DEFAULT))

    if args.against is not None:
        printed = model_amplification(args.against)
        if printed is None:
            print(f"{args.against} carries no NCONVTIME report; the run did "
                  "not enable conversion_time_level, or it was made by a binary "
                  "built before the guard")
            return 2
        dt = printed["timestep_minutes"]
        rung = args.rung or str(cfg["model"]["resolution"]).upper()
        column = Column(cfg, rung, pnu)
        here = column.verdict(dt)
        # A CHECK WITH A RIGHT ANSWER, and BOTH ARMS ARE IN IT. The same
        # iteration from the same start over the same span must return the same
        # numbers here and in the model, arm for arm and wavenumber for
        # wavenumber; a difference is a defect in one implementation and not a
        # question about either. The tolerance is what two orderings of the
        # same double-precision arithmetic can differ by over the span, and
        # nothing about the physics.
        problems = []
        for key, label in (("on", "measured arm"), ("off", "control arm")):
            if abs(printed[key] - here[key]) > 1e-9 * max(1.0, here[key]):
                problems.append(f"the {label} differs: model {printed[key]!r}, "
                                f"here {here[key]!r}")
            if printed[key + "_wavenumber"] != here[key + "_wavenumber"]:
                problems.append(
                    f"the {label}'s worst wavenumber differs: model "
                    f"{printed[key + '_wavenumber']}, here "
                    f"{here[key + '_wavenumber']}")
        print(f"{rung} at dt {dt} min, PNU {pnu}")
        print(f"  measured arm  model {printed['on']:.12f} at n="
              f"{printed['on_wavenumber']}   here {here['on']:.12f} at n="
              f"{here['on_wavenumber']}")
        print(f"  control arm   model {printed['off']:.12f} at n="
              f"{printed['off_wavenumber']}   here {here['off']:.12f} at n="
              f"{here['off_wavenumber']}")
        print(f"  effect {here['effect']:.3e} against a resolution of "
              f"{here['resolution']:.3e}: "
              f"{'refused' if here['refused'] else 'passed'}")
        if problems:
            print("  THEY DISAGREE, and the guard's refusals rest on the "
                  "model's:")
            for line in problems:
                print("    " + line)
            return 1
        print("  they agree, arm for arm")
        return 0

    if args.rung and args.sweep:
        wanted = [(args.rung.upper(), float(x))
                  for x in args.sweep.replace(",", " ").split()]
    elif args.rung:
        wanted = [(args.rung.upper(), args.dt)]
    else:
        wanted = [(r, dt) for r, dt in rungs.ESCALATION_ROUTE]
        for r in sorted(rungs.STABILITY_CEILING_MINUTES):
            if not any(r == w for w, _ in wanted):
                wanted.append((r, rungs.stability_ceiling(r)))

    rows = []
    for rung, dt in wanted:
        if rung not in cfg["model"]["hyperdiffusion"]["timescales_days"]:
            continue
        column = Column(cfg, rung, pnu)
        guard = column.gravity_wave_limit_minutes()
        boundary = column.boundary_minutes()
        row = {"rung": rung, "timestep_minutes": dt, "pnu": pnu,
               "gravity_wave_limit_minutes": round(guard, 4),
               "gravity_wave_limit_from_structure_matrix_minutes":
                   round(column.structure_gravity_wave_limit_minutes(), 4),
               "stability_boundary_minutes":
                   None if boundary is None else round(boundary, 4),
               "guard_over_boundary":
                   None if boundary is None else round(guard / boundary, 4)}
        if dt is not None:
            grow, jn = column.worst_growth(dt)
            control, _ = column.worst_growth(dt, nconvtime=False)
            row["growth_per_step"] = round(grow, 8)
            row["worst_total_wavenumber"] = jn
            row["growth_per_step_nconvtime_off"] = round(control, 8)
            # THE MODEL'S OWN TWO-ARM TEST, so the artifact carries both the
            # exact answer and what the guard would actually decide.
            v = column.verdict(dt)
            row["guard"] = {
                "measured_arm": v["on"], "measured_arm_wavenumber": v["on_wavenumber"],
                "control_arm": v["off"], "control_arm_wavenumber": v["off_wavenumber"],
                "effect": v["effect"], "resolution": v["resolution"],
                "refused": v["refused"]}
        rows.append(row)
        print(f"{rung:<5} dt {str(dt):>6}  guard {guard:7.2f} min  "
              f"boundary "
              f"{'none at any step' if boundary is None else f'{boundary:7.2f} min'}"
              + (f"  growth {row['growth_per_step']:.6f} at n="
                 f"{row['worst_total_wavenumber']}" if dt is not None else ""))

    result = {
        "what": ("the one-step amplification of the model's linearised "
                 "adiabatic step with NCONVTIME on, and the coarsest timestep "
                 "at which it does not grow"),
        "provenance": provenance(args.config),
        "rows": rows,
        "note": ("guard_over_boundary is how much coarser a step the model's "
                 "explicit gravity-wave limit admits than the scheme is stable "
                 "at. A null boundary means the modification grows at every "
                 "step down to 0.05 min, which is what PNU = 0 gives: the "
                 "Robert-Asselin filter is the only thing damping the leapfrog "
                 "computational mode the term feeds. "
                 "gravity_wave_limit_from_structure_matrix_minutes is the same "
                 "limit taken on the model's own fastest external mode rather "
                 "than on a single-level estimate of it, and the two agreeing "
                 "is what says the guard was aimed at the wrong mode rather "
                 "than mis-derived. world-bt3b."),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
