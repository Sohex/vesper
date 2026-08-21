#!/bin/bash
# Do mpimod_omp's 39 routines compute what mpimod's compute?
#
#   exoplasim/scripts/verify_omp_collectives.sh [nthreads] [nlat]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# parallel layer. Nothing here is about the simulated planet.
#
# Each collective has an answer that can be written down in advance, so this
# checks against those answers rather than against another run of the model.
#
# THE VALUES DEPEND ON THE MODE INDEX, and that is the point. A first version
# had every thread contribute a constant, which cannot fail on the mistake this
# code invites -- summing or scattering the wrong slice -- because a uniform
# contribution makes any index confusion cancel. The same flaw shipped once
# already in this project, in the filter fold's indexing test.
set -euo pipefail
nthreads="${1:-16}"
nlat="${2:-32}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/vendor/exoplasim/exoplasim/plasim/src"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

cat > "$W/resmod.f90" <<RES
      module resmod
      parameter(NLAT_ATM = $nlat)
      parameter(NLEV_ATM = 10)
      parameter(NPRO_ATM = $nthreads)
      end module resmod
RES

# restartmod's entry points, stubbed: this exercises the collectives, not the
# restart file format.
cat > "$W/stubs.f90" <<'STUB'
      subroutine get_restart_array(yn,pa,k1,k2,k3)
      character (len=*) :: yn
      integer :: k1,k2,k3
      real :: pa(k2,k3)
      return
      end
      subroutine put_restart_array(yn,pa,k1,k2,k3)
      character (len=*) :: yn
      integer :: k1,k2,k3
      real :: pa(k2,k3)
      return
      end
      subroutine get_surf_array(yn,pa,k1,k2,kread)
      character (len=*) :: yn
      integer :: k1,k2,kread(1)
      real :: pa(k1,k2)
      kread(1) = 0
      return
      end
STUB

cp "$SRC/plasimmod.f90" "$SRC/mpimod_omp.f90" "$W/"
cp "$ROOT/exoplasim/scripts/verify_omp_collectives.f90" "$W/"
( cd "$W" && gfortran -fopenmp -O2 -cpp -ffixed-line-length-132 -fdefault-real-8 \
     -J "$W" -o coll.x resmod.f90 plasimmod.f90 mpimod_omp.f90 stubs.f90 \
     verify_omp_collectives.f90 )
echo "NLAT $nlat, $nthreads threads"
OMP_STACKSIZE=512M "$W/coll.x"
