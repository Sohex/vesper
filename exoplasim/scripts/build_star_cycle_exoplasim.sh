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
executable="most_plasim_t42_l10_p8.x"
base_sha="eb8e9e1c0127940e607828899ff5dea6653c9835cf9d90b79215f7fc2d0f5261"

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
(cd "$package_dir" && ./compile.sh -n 8 -p 8 -r T42 -v 10)
mkdir -p "$target_dir"
cp -a "$run_dir/." "$target_dir/"
sha256sum "$target_dir/$executable"
