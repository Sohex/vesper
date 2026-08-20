# Conventions this project holds to

Provenance travels with every artifact: config, input hashes, software
versions and the terrain hash go into each run manifest and each analysis
report.

Thresholds are fixed before results are seen. The albedo bracket's convergence
and marginality criteria were written into the script's docstring before any
run finished, and were then applied against a result that cleared one of them
by 0.05 K.

Estimates that cannot be verified are bracketed rather than guessed, and the
bracket is reported. This is how albedo, evaporation and the carve verdict are
all handled.

Claims are checked against the artifact rather than the documentation. Several
findings in this project's history came from comparing two products that were
supposed to agree and finding they did not.
