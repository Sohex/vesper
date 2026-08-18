#!/usr/bin/env bash
set -euo pipefail

component_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
project_dir="$(cd "$component_dir/.." && pwd)"
package_dir="$project_dir/.venv/lib/python3.12/site-packages/exoplasim"
source_file="$package_dir/plasim/src/radmod.f90"
run_dir="$package_dir/plasim/run"
bin_dir="$package_dir/plasim/bin"
target_dir="$component_dir/inputs/exoplasim_cycle_t42"
patch_file="$component_dir/patches/exoplasim-3.4.2-star-cycle.patch"
executable="most_plasim_t42_l10_p16.x"
# The base this patch applies ON TOP OF, not pristine ExoPlaSim 3.4.2.
#
# radmod.f90 now carries three resident patches besides this one, and the pin has
# to name the stack rather than a version:
#
#   pristine 3.4.2                eb8e9e1c0127940e607828899ff5dea6653c9835c...
#   + ozone-band-weights
#   + prescribed-dust             (DUST-11; radmod and surfmod)
#   + h2o-shortwave-weight        (PHYS-1; radmod)
#                                 -> the sha below
#
# Regenerated 2026-08-17 against that base, which is what the comment this
# replaces said would be necessary and it was right: hunk 2 FAILED outright,
# because the prescribed-dust patch had appended its own keys to the same
# radmod_nl continuation the cycle keys attach to. The other four hunks applied
# at offsets. A patch that fails one hunk and offsets four is a patch that has to
# be regenerated rather than re-checked, and the cycle keys now sit after the
# dust keys instead of after minwavel.
#
# If a FOURTH patch lands on radmod.f90 this sha moves again and this patch has
# to be regenerated against the new base. That is not a nuisance to be worked
# around -- it is the check that stops us building a cycle binary on a source we
# have not looked at.
base_sha="5dd9ddfb765d4391e31ba5af89d64bf0dde0db528f9622cb0f6b7faf73e89083"

if [[ ! -f "$source_file" || ! -f "$run_dir/$executable" || ! -d "$bin_dir" ]]; then
  echo "ExoPlaSim 3.4.2 source or baseline executable is missing" >&2
  exit 1
fi
if [[ "$(sha256sum "$source_file" | cut -d' ' -f1)" != "$base_sha" ]]; then
  echo "radmod.f90 does not match the pinned ExoPlaSim 3.4.2 source" >&2
  exit 1
fi

scratch_dir="$(mktemp -d)"
cp "$run_dir/$executable" "$scratch_dir/$executable"
restore_vendor_tree() {
  patch --silent --reverse --strip=1 --directory="$package_dir" < "$patch_file" || true
  cp "$scratch_dir/$executable" "$run_dir/$executable"
  cp "$scratch_dir/$executable" "$bin_dir/$executable"
  rm -rf "$scratch_dir"
}
trap restore_vendor_tree EXIT

patch --forward --strip=1 --directory="$package_dir" < "$patch_file"
# Rank count must match model.ncpus in config/planet.yaml: the run_id encodes
# it, and a cycle run is long enough that 16 ranks against 8 is hours.
(cd "$package_dir" && ./compile.sh -n 16 -p 16 -r T42 -v 10)
mkdir -p "$target_dir"
cp -a "$run_dir/." "$target_dir/"
# Keep only the executable this script actually built with the cycle patch.
# cp -a brings the whole run directory, so every OTHER resolution and rank count
# arrives too -- steady binaries sitting in a directory named "cycle". They fail
# loudly rather than silently, because a steady binary cannot parse nsolcycle in
# the namelist, but a directory should hold what its name claims.
find "$target_dir" -maxdepth 1 -name 'most_plasim_*.x' ! -name "$executable" -delete
sha256sum "$target_dir/$executable"

# Record WHAT this binary was built from, by content rather than by timestamp.
# check_consistency.py used to compare mtimes, which reports a perfectly good
# binary as stale the moment anything rewrites the patch file without changing
# it -- a git stash/pop does exactly that. A hash cannot produce that false
# alarm, and unlike an mtime it also catches a patch edited in place.
manifest="$component_dir/patches/cycle_binary_manifest.json"
cat > "$manifest" <<JSON
{
  "note": "Written by build_star_cycle_exoplasim.sh. Identifies the cycle executable by content. Do not edit.",
  "built_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "executable": "$executable",
  "executable_sha256": "$(sha256sum "$target_dir/$executable" | cut -d" " -f1)",
  "patch": "$(basename "$patch_file")",
  "patch_sha256": "$(sha256sum "$patch_file" | cut -d" " -f1)",
  "base_radmod_sha256": "$base_sha"
}
JSON
echo "wrote $manifest"
