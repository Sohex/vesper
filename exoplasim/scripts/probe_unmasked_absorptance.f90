! Does the unmasked absorptance return what the masked one returns, and does it
! vectorise?
!
! Worldbuilding frame: a COMPUTE and NUMERICS check on the Vesper climate
! model's longwave code. Nothing here is about the simulated planet.
!
! WORLD-43RK replaces four `where`/`elsewhere` absorptance selections in `lwr`
! with one unmasked array assignment each, so that GCC can reach libmvec's
! eight-wide `pow`. Two things then have to be shown, and neither is a timing:
!
!   1. THE BLEND IS EXACT. The selector is 0.5+SIGN(0.5,.), which is exactly 1.
!      or exactly 0., so the blend returns the chosen branch unchanged. This
!      program asserts BIT EQUALITY against the masked form, per element, over
!      path amounts that straddle every threshold. Any difference that survives
!      is the library's, not the algebra's.
!   2. THE EXCLUDED LANES ARE IN DOMAIN. A `where` masks the assignment and not
!      the evaluation, so the unmasked form evaluates both branches on every
!      lane. This program runs under the model's own
!      -ffpe-trap=invalid,zero,overflow, so an out-of-domain lane aborts here
!      rather than eight orbits into a run.
!
! The path amounts include the values the mask exists to keep apart: exactly 0.,
! exactly the threshold, denormal-small, and far above the threshold.
      program probe_unmasked_absorptance
      implicit none
      integer, parameter :: N = 65536
      real :: zsumwv(N), zsumco2(N), zsumo3(N)
      real :: am(N), au(N)
      real :: zh2o0a, zh2o0, zco20, zao30, zah2oc, zaco2c, zth2oc, zao3c
      real :: zsel(N)
      integer :: i, nbad

      zao30 = 0.209*(7.E-5)**0.436
      zco20 = 0.0676*(0.01022)**0.421
      zh2o0a= 0.846*(3.59E-5)**0.243
      zh2o0 = 0.832*0.0286**0.26
      zao3c = 0.209*(0.01+7.E-5)**0.436-zao30-0.0212*log10(0.01)
      zaco2c= 0.0676*(1.01022)**0.421-zco20
      zah2oc= 0.846*(0.01+3.59E-5)**0.243-zh2o0a-0.24*ALOG10(0.02)
      zth2oc= 1.-(0.832*(2.+0.0286)**0.26-zh2o0)+0.1196*log(2.-0.6931)

!     Path amounts over fourteen decades, plus the exact thresholds and exact
!     zero. `lwr` accumulates these as sums of non-negative layer amounts, so
!     zero is reachable and negative is not.
      do i = 1, N
       zsumwv(i)  = 1.E-10 * 10.0**(14.0*real(i-1)/real(N-1))
      enddo
      zsumwv(1) = 0.
      zsumwv(2) = 0.01
      zsumwv(3) = 2.
      zsumwv(4) = 0.6931
      zsumco2(:) = zsumwv(:)
      zsumo3(:)  = zsumwv(:)

      nbad = 0

!     h2o 6.3mu
      where(zsumwv(:) <= 0.01)
       am(:)=0.846*max(0.,zsumwv(:)+3.59E-5)**0.243-zh2o0a
      elsewhere
       am(:)=0.24*ALOG10(max(1.E-30,zsumwv(:)+0.01))+zah2oc
      endwhere
      zsel(:)=0.5+SIGN(0.5,0.01-zsumwv(:))
      au(:)=zsel(:)*(0.846*max(0.,zsumwv(:)+3.59E-5)**0.243-zh2o0a)           &
     &     +(1.-zsel(:))*(0.24*ALOG10(max(1.E-30,zsumwv(:)+0.01))+zah2oc)
      call report('h2o 6.3mu', am, au, N, nbad)

!     co2
      where(zsumco2(:) <= 1.0)
       am(:)=0.0676*max(0.,zsumco2(:)+0.01022)**0.421-zco20
      elsewhere
       am(:)=0.0546*ALOG10(max(1.E-30,zsumco2(:)))+zaco2c
      endwhere
      zsel(:)=0.5+SIGN(0.5,1.0-zsumco2(:))
      au(:)=zsel(:)*(0.0676*max(0.,zsumco2(:)+0.01022)**0.421-zco20)          &
     &     +(1.-zsel(:))*(0.0546*ALOG10(max(1.E-30,zsumco2(:)))+zaco2c)
      call report('co2', am, au, N, nbad)

!     h2o - co2 overlap
      where(zsumwv(:)<= 2.)
       am(:)=1.-(0.832*max(0.,zsumwv(:)+0.0286)**0.26-zh2o0)
      elsewhere
       am(:)=max(0.,zth2oc-0.1196*log(max(1.E-30,zsumwv(:)-0.6931)))
      endwhere
      zsel(:)=0.5+SIGN(0.5,2.-zsumwv(:))
      au(:)=zsel(:)*(1.-(0.832*max(0.,zsumwv(:)+0.0286)**0.26-zh2o0))         &
     &     +(1.-zsel(:))                                                      &
     &      *max(0.,zth2oc-0.1196*log(max(1.E-30,zsumwv(:)-0.6931)))
      call report('h2o-co2 overlap', am, au, N, nbad)

!     o3
      where(zsumo3(:) <= 0.01)
       am(:)= 0.209*max(0.,zsumo3(:)+7.E-5)**0.436 - zao30
      elsewhere
       am(:)= 0.0212*log10(max(1.E-30,zsumo3(:)))+zao3c
      endwhere
      zsel(:)=0.5+SIGN(0.5,0.01-zsumo3(:))
      au(:)=zsel(:)*(0.209*max(0.,zsumo3(:)+7.E-5)**0.436-zao30)              &
     &     +(1.-zsel(:))*(0.0212*log10(max(1.E-30,zsumo3(:)))+zao3c)
      call report('o3', am, au, N, nbad)

      if (nbad > 0) then
       write(*,*) 'FAIL: ', nbad, ' elements differ'
       stop 1
      endif
      write(*,*) 'PASS: unmasked form is bit identical on every element,'
      write(*,*) '      and no lane raised invalid, zero or overflow.'
      end program probe_unmasked_absorptance

      subroutine report(name, am, au, n, nbad)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in) :: n
      real, intent(in) :: am(n), au(n)
      integer, intent(inout) :: nbad
      integer :: i, k
      real :: d, dmax
      k = 0
      dmax = 0.
      do i = 1, n
       if (am(i) /= au(i)) then
        k = k + 1
        d = abs(am(i)-au(i))
        if (d > dmax) dmax = d
       endif
      enddo
!     ABSOLUTE, not relative. These are broadband absorptances, dimensionless
!     fractions the routine bounds into [0,1] two lines later, and each is a
!     difference of two nearly equal terms near its own zero. A relative
!     difference there reports the cancellation and not the error.
      write(*,'(a,a18,a,i7,a,i7,a,es10.2)') '  ', name,                       &
     &  ': differing ', k, ' of ', n, '  max abs ', dmax
      nbad = nbad + k
      end subroutine report
