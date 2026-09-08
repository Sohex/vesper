# Archived builds

Identity of superseded World Orogen exports whose payload has been
deleted. `identity.json` carries the hashes, seed, region count, planet
parameters, land/sea mask summary and basin counts -- everything needed
to recognise a build and to date a result computed from it.

The payload is regenerable: Orogen rebuilds terrain from the planet code
plus a carve list in one pass, and both are kept. That is why these can
go and a climate run's output cannot.

`lib/orogen.py` remains the registry and still refuses a build it has not
been checked against. Per-build hydrography products under
`hydrography/data/<build>/` were not touched.

Archived 2026-09-08.
