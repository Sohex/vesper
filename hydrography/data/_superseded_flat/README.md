# Superseded flat products

These are the hydrography products that used to live directly in
`hydrography/data/`, from before drainage products were namespaced by build.

`carve_list.json` declares terrain hash `821aa71b...`, which is
**precarve-unzoned** -- the original pre-carve build. `basins.nc` here holds
**2,107** basins against the active build's 3,629.

They are kept rather than deleted because they are not duplicates of anything in
a per-build directory and they date results that were computed from them. They
are moved out of the flat path because every script that used to default there
would silently pair one terrain's rows with another's columns, and one of them
-- `carve_verdict.py` -- was doing exactly that: reading basin ids from this
directory while computing verdicts from the active build, joined by a `zip` that
truncated to the shorter without error.

Nothing should read from this directory. If you need a product for a specific
build, it is in `hydrography/data/<build>/`.
