! =================
! SUBROUTINE INIGAU
! =================

! Gaussian abscissas and weights, to machine precision.
!
! WHY THIS WAS REWRITTEN. The previous implementation evaluated the Legendre
! polynomial as a TRIGONOMETRIC SERIES, sum of z1*cos(j*acos(p)), and took the
! weight as z4*(1-z^2)/z5^2 -- which squares that series' error into the weight.
! Measured against an extended-precision reference at NLAT 192, its weights were
! wrong by 9.9e-11 at the pole-most latitude and by a few times 1e-12 near it,
! while SHTns's own weights were exact to 2.7e-16 and bit-identical over most of
! the hemisphere. numpy's leggauss is also wrong there, by 4.2e-11, so three
! independent double-precision implementations disagreed and only one was right.
!
! WHAT IT COST. The weights enter the FORWARD transform and not the inverse --
! `fc2sp` carries gwd, `uv2dv` and `mktend` carry gwdc = gwd/cos^2, and nothing
! in the synthesis direction touches them. So every grid-to-spectral transform
! this model has ever done carried that error, and it showed up as a floor of
! about 2e-11 on the analysis arms of `verify_shtns_equivalence.sh` -- field
! blind, band-limit blind, and present in the scalar arms as much as the vector
! ones. CLIM-61 and `exoplasim/notes/shtns-viability.md` have the measurements.
!
! THE METHOD IS THE TEXTBOOK ONE, and that is the point: Newton on P_n by the
! three-term recurrence from Tricomi's asymptotic start, then
! w = 2/((1-x^2) P'^2). `exoplasim/scripts/gauss_weight_reference.py` is the
! same algorithm in extended precision and is what this is checked against.
!
! WHAT LIMITS IT NOW, since it is no longer the algorithm. Evaluating P_n by
! the recurrence accumulates about n*eps, and the weight squares the derivative,
! so the floor is roughly 2*n*eps -- about 1.1e-13 at NLAT 256, which is what is
! measured. Converging the ANGLE instead of the node was tried, to avoid forming
! 1-z^2 near the pole where z is 0.99993: it changed the weights by nothing
! outside noise, because the recurrence and not the cancellation is the limit,
! and it made the nodes worse near the equator. Beating 1e-13 here needs a
! different algorithm, not a rearrangement of this one.

subroutine inigau(klat,pz0,pzw)        ! pz0 & pzw are (kind=8) reals !!!
implicit none
integer                  :: klat       ! Number of Gaussian latitudes
real (kind=8)            :: pz0(klat)  ! Gaussian abscissas, descending from +1
real (kind=8)            :: pzw(klat)  ! Gaussian weights
integer                  :: jlat       ! Latitudinal loop index
integer                  :: jiter      ! Iteration loop index
integer                  :: j          ! Recurrence loop index
integer      , parameter :: NITER = 100 ! Maximum # of iterations
real (kind=8), parameter :: PI    =  3.141592653589793_8
real (kind=8), parameter :: ZEPS  =  1.0e-15 ! Convergence criterion, on the STEP
real (kind=8) :: z      ! the node being converged
real (kind=8) :: zp0    ! P_{n-1}(z)
real (kind=8) :: zp1    ! P_n(z)
real (kind=8) :: zt     ! recurrence scratch
real (kind=8) :: zdp    ! P'_n(z)
real (kind=8) :: zstep  ! Newton step

do jlat = 1 , klat/2

!  Tricomi's asymptotic root, good to about 1e-5, which Newton squares down in
!  three or four passes.

   z = cos(PI * (4.0_8*jlat - 1.0_8) / (4.0_8*klat + 2.0_8))

   do jiter = 1 , NITER
      zp0 = 1.0_8
      zp1 = z
      do j = 2 , klat
         zt  = zp1
         zp1 = ((2*j - 1) * z * zp1 - (j - 1) * zp0) / j
         zp0 = zt
      enddo ! j
      zdp   = klat * (z * zp1 - zp0) / (z * z - 1.0_8)
      zstep = zp1 / zdp
      z     = z - zstep
      if (abs(zstep) < ZEPS) exit ! converged
   enddo ! jiter

!  The derivative again at the converged node: the one from the last iteration
!  belongs to the node before the final step, and the weight goes as its square.

   zp0 = 1.0_8
   zp1 = z
   do j = 2 , klat
      zt  = zp1
      zp1 = ((2*j - 1) * z * zp1 - (j - 1) * zp0) / j
      zp0 = zt
   enddo ! j
   zdp = klat * (z * zp1 - zp0) / (z * z - 1.0_8)

   pz0(jlat) = z
   pzw(jlat) = 2.0_8 / ((1.0_8 - z * z) * zdp * zdp)
   pz0(klat-jlat+1) = -z
   pzw(klat-jlat+1) = pzw(jlat)
enddo ! jlat

return
end subroutine inigau
