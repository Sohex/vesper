#!/usr/bin/env python3
"""Who has the machine, and is it quiet enough to time something on.

Wall clock, throughput and per-step cost are measurements of a MACHINE STATE as
much as of a model. This project runs many agents at once on one host, so a
timing taken while someone else is integrating is not a slow number, it is a
number of a different experiment. Three separate measurements were contaminated
that way in one session before this file existed.

    python scripts/machine.py --check                 # may I start heavy work?
    python scripts/machine.py --claim "spat-11 wall arms" --minutes 30
    python scripts/machine.py --release

A CLAIM IS ADVISORY AND THAT IS DELIBERATE. Nothing here can stop a process, and
a hard lock would strand the machine whenever an agent died holding it. What it
does is make the state VISIBLE and give a claim an expiry, so the failure mode is
a stale claim someone can see and override rather than an invisible collision.

Claim before a timing measurement. Check before anything heavy. A claim says
"my numbers are worthless if you start now"; it does not say "do not use the
machine" -- work that is not timing-sensitive may contend freely and should say
so when it starts.
"""
import argparse, json, os, time
from pathlib import Path

# OUTSIDE THE REPOSITORY, DELIBERATELY. Every fan-out agent works in its own git
# worktree, so a claim file under the project root would be a DIFFERENT file for
# each of them and would coordinate nothing -- which is the one failure mode this
# exists to prevent. One host, one path, shared by every worktree and every
# session.
CLAIM = Path("/tmp/vesper-machine-claim.json")
# Above this, a timing measurement is not worth taking: the host has 32 logical
# cores and the model runs 16 threads, so one integration is load ~16 on its own.
QUIET_LOAD = 4.0


def _read():
    if not CLAIM.is_file():
        return None
    try:
        c = json.loads(CLAIM.read_text())
    except Exception:
        return None
    if time.time() > c.get("expires_at", 0):
        return None                      # expired claims are not claims
    return c


# What a heavy job looks like in `ps`. Named rather than inferred from CPU alone,
# because a compiler at 98 per cent is load and a model at 1595 per cent is a
# different kind of problem, and the reader needs to know which.
HEAVY = ("most_plasim", "genie.exe", "run_exoplasim", "continue_exoplasim",
         "cgenie_cost", "stability_probe", "make_profile_bed", "build_model")


def _running():
    """Heavy jobs on this host, whoever started them.

    A CLAIM IS NOT ENOUGH ON ITS OWN. Work that is insensitive to load still
    GENERATES load, so it never claims -- it has no reason to. A checker that
    only read the claim file would see nothing, conclude the host was quiet, and
    be wrong. That is not hypothetical: an eighty-five-orbit precision pair and a
    serial cGENIE cost probe contended for an hour, each invisible to the other,
    with no claim held by either. So this reads the process table, which needs no
    cooperation from anyone.
    """
    try:
        import subprocess
        out = subprocess.run(["ps", "-eo", "pcpu,args"], capture_output=True,
                             text=True, timeout=10).stdout
    except Exception:
        return []
    found = []
    for line in out.splitlines()[1:]:
        for name in HEAVY:
            if name in line and "machine.py" not in line:
                try:
                    cpu = float(line.split(None, 1)[0])
                except (ValueError, IndexError):
                    cpu = 0.0
                if cpu >= 20.0:
                    found.append((cpu, line.split(None, 1)[1][:88]))
                break
    return sorted(found, reverse=True)


def check(verbose=True):
    load1 = os.getloadavg()[0]
    c = _read()
    busy = _running()
    if verbose:
        print(f"load {load1:.2f} over {os.cpu_count()} logical cores")
        if c:
            left = (c["expires_at"] - time.time()) / 60.0
            print(f"CLAIMED by {c['who']} for {c['purpose']!r}, {left:.0f} min left")
            print("  A timing measurement is in progress. Starting heavy work now")
            print("  invalidates it. Wait, or say in your report that you contended.")
        else:
            print("unclaimed")
        if busy:
            print(f"  {len(busy)} heavy job(s) running, claimed or not:")
            for cpu, cmd in busy[:5]:
                print(f"    {cpu:7.1f}% {cmd}")
        if load1 > QUIET_LOAD:
            print(f"  NOT QUIET: load {load1:.2f} is above {QUIET_LOAD}. Do not take a")
            print("  wall-clock number here. Retired instructions under")
            print("  OMP_WAIT_POLICY=passive survive this; wall clock does not.")
        else:
            print("  quiet enough to time")
    return 1 if (c or busy or load1 > QUIET_LOAD) else 0


def workers(cap=8):
    """How many worker processes a load-INSENSITIVE job may take right now.

    A gate that spawns a process per unit is not timing-sensitive and never
    claims, but it still generates load, and load is what turns somebody else's
    claimed measurement into a different experiment. So this reads the same two
    sources `check` does -- a live claim, and heavy jobs in the process table --
    and drops to one worker while either says the host is in use. A job that can
    wait finishing slowly is a cost to nobody; the same job taking eight cores
    across a claimed measurement costs that measurement.

    Otherwise it is `cap` or a quarter of the logical cores, whichever is
    smaller, so a gate on a shared host is never the heavy job someone else has
    to wait for. Nothing here is a permission system: a caller may pass its own
    count, and this is the number it should use when it has no reason to.
    """
    cores = os.cpu_count() or 1
    if _read() or _running():
        return 1
    return max(1, min(cap, cores // 4))


def claim(who, purpose, minutes):
    c = _read()
    if c:
        print(f"REFUSED: already claimed by {c['who']} for {c['purpose']!r}")
        return 1
    CLAIM.write_text(json.dumps({
        "who": who, "purpose": purpose,
        "claimed_at": time.time(), "expires_at": time.time() + minutes * 60,
        "load_at_claim": os.getloadavg()[0],
    }, indent=1))
    print(f"claimed for {minutes} min at load {os.getloadavg()[0]:.2f}: {purpose}")
    print("RECORD THE LOAD WITH YOUR NUMBERS. A timing without the machine state")
    print("it was taken under cannot be compared against a later one.")
    return 0


def release():
    if CLAIM.is_file():
        CLAIM.unlink()
        print("released")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true")
    p.add_argument("--claim", metavar="PURPOSE")
    p.add_argument("--who", default=os.environ.get("CLAUDE_AGENT", "unnamed"))
    p.add_argument("--minutes", type=int, default=30)
    p.add_argument("--release", action="store_true")
    a = p.parse_args()
    if a.claim:
        raise SystemExit(claim(a.who, a.claim, a.minutes))
    if a.release:
        raise SystemExit(release())
    raise SystemExit(check())
