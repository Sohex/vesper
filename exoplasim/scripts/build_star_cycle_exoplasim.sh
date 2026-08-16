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
# radmod.f90 now carries the ozone band-weight patch as well, so pinning the
# pristine 3.4.2 sha made this script refuse to run at all. The guard is still
# worth having -- it is what stops us patching a source we have not checked --
# but its premise is "3.4.2 + ozone band weights", and it has to say so.
#
#   pristine 3.4.2                eb8e9e1c0127940e607828899ff5dea6653c9835c...
#   + exoplasim-3.4.2-ozone-band-weights.patch   -> the sha below
#
# Verified 2026-08-16: the star-cycle patch applies to that base with all five
# hunks clean at offsets of 10 to 12 lines, which is exactly the shift the ozone
# patch introduces. The two patches coexist. If a THIRD patch lands on
# radmod.f90, this sha moves again and the comment above needs another line.
base_sha="7fd39458a87a0bc042d17b0b93c2e45cbfc4fc04a560965c974efa734019c0ba"

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
