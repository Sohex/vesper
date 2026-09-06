!     ==================================================================
!     Does the canopy snow store conserve the snow that falls into it, and
!     does it reduce to no store at all where the model has no canopy?
!
!     cansnowmod puts a reservoir between a vegetated gridcell's snowfall and
!     its ground snowpack. Two things about it have right answers written down
!     in advance rather than taken from another run:
!
!     THE REDUCTION. With no plant area index -- which is what this world
!     reports until the vegetation component reports one -- there is no canopy
!     to hold anything, so the throughfall has to be the snowfall bit for bit
!     and the interception weight has to be zero. The ground pack then receives
!     exactly what it received before this routine existed.
!
!     THE MASS BALANCE. Over any driven sequence of steps, the snow that fell
!     equals the throughfall that reached the ground plus the snow still held in
!     the canopy. That is an identity and not a comparison: it is accumulated
!     OUTSIDE the routine, from its inputs and outputs alone, so an arithmetic
!     path that loses or invents mass fails here rather than showing up as a
!     surface water balance that does not close.
!
!     And four limits: the store never exceeds its capacity, never goes
!     negative, empties on the step the air goes above freezing, and approaches
!     the load its snowfall and its timescale imply.
!
!     Built and run by verify_canopy_snow.sh.
!
!     Worldbuilding frame: a correctness check on the snow held in the canopy of
!     a simulated planet's vegetated land.
!     ==================================================================
      program verify_canopy_snow
      use cansnowmod
      implicit none

!     Essery (2013) Eq. (3) and the sentences under it, as landmod declares
!     them. `tau` is a DURATION and is in seconds: the paper's ten days is a
!     physical unloading time and not a count of this world's days, which are
!     30 hours long.
      real, parameter :: ceff  = 0.25          ! interception efficiency
      real, parameter :: tau   = 8.64e5        ! unloading timescale (s)
      real, parameter :: capai = 2.0e-4        ! capacity per unit plant area
                                               ! index (m water equivalent)
      real, parameter :: tmelt = 273.16        ! the model's melting point (K)
      real, parameter :: cold  = 263.0         ! an air temperature below it
      real, parameter :: dt    = 1800.0        ! a step of the model's order
      real, parameter :: off   = -1.0          ! the inert plant area index

      integer :: nfail
      nfail = 0

      call reduction(nfail)
      call conserves(nfail)
      call capacity(nfail)
      call flushes(nfail)
      call approaches(nfail)

      write(*,*)
      if (nfail == 0) then
         write(*,*) "PASS: the store reduces to no store and conserves what "// &
     &              "falls into it"
      else
         write(*,'(a,i4)') " FAIL: checks failed: ", nfail
         stop 1
      endif

      contains

!     THE REDUCTION. No plant area index, so no capacity: the throughfall is
!     the snowfall exactly, at every cover and every snowfall rate, and the
!     interception weight is zero. Bitwise, because that is what makes the
!     ground pack's tendency the expression it already was.
      subroutine reduction(kfail)
      integer, intent(inout) :: kfail
      integer :: j,k
      real :: zcov,zsf,zcan,zthr,zfcan,zfint
      do j=0,20
         zcov = real(j)/20.0
         do k=0,10
            zsf = real(k)*1.0e-8
            call canopysnowstep(zcov,off,capai,ceff,tau,zsf,0.0,cold,   &
     &                          tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
            if (zthr /= zsf .or. zcan /= 0.0 .or. zfint /= 0.0) then
               write(*,'(a,2e16.8)') " FAIL reduction: snowfall, through = ", &
     &                               zsf,zthr
               kfail = kfail + 1
            endif
         enddo
      enddo
      write(*,*) "reduction at no plant area index: throughfall is the snowfall"
      end subroutine reduction

!     THE MASS BALANCE. Drive the store through a season with a snowfall that
!     starts, stops and restarts, and an air temperature that crosses the
!     melting point in the middle of it. Accumulate what fell and what reached
!     the ground here, from the routine's own arguments, and compare.
      subroutine conserves(kfail)
      integer, intent(inout) :: kfail
      integer :: j
      real :: zfell,zreached,zsf,ztair,zcan,zthr,zfcan,zfint,zres,zscale
      zfell    = 0.0
      zreached = 0.0
      zcan     = 0.0
      do j=1,4000
         zsf = 0.0
         if (j <= 1200 .or. (j > 2000 .and. j <= 2600)) zsf = 2.0e-8
         if (j > 3000 .and. j <= 3200) zsf = 4.0e-7
         ztair = cold
         if (j > 1500 .and. j <= 1560) ztair = tmelt + 2.0
         if (j > 3500) ztair = tmelt + 5.0
         call canopysnowstep(0.8,3.1,capai,ceff,tau,zsf,zcan,ztair,     &
     &                       tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
         zfell    = zfell    + zsf*dt
         zreached = zreached + zthr*dt
         if (zthr < 0.0) then
            write(*,'(a,i6,e16.8)') " FAIL negative throughfall at step ",j,zthr
            kfail = kfail + 1
         endif
         if (zcan < 0.0) then
            write(*,'(a,i6,e16.8)') " FAIL negative store at step ",j,zcan
            kfail = kfail + 1
         endif
         if (zfcan < 0.0 .or. zfcan > 1.0) then
            write(*,'(a,i6,e16.8)') " FAIL cover fraction out of range ",j,zfcan
            kfail = kfail + 1
         endif
      enddo
      zres   = zfell - zreached - zcan
      zscale = zfell
      write(*,'(a,e14.6,a,e14.6)') " fell ",zfell," m w.e., reached the "// &
     &                             "ground ",zreached
      write(*,'(a,e12.4)') " mass residual over 4000 steps, relative: ",  &
     &                     abs(zres)/zscale
      if (abs(zres) > 1.0e-12*zscale) then
         write(*,*) " FAIL: the store does not conserve what falls into it"
         kfail = kfail + 1
      endif
      end subroutine conserves

!     THE CAPACITY. Snow hard enough for long enough and the load has to stop
!     at cover*pai*capai and the rest has to reach the ground.
      subroutine capacity(kfail)
      integer, intent(inout) :: kfail
      integer :: j
      real :: zcan,zthr,zfcan,zfint,zcap
      zcap = 0.8*3.1*capai
      zcan = 0.0
      do j=1,2000
         call canopysnowstep(0.8,3.1,capai,ceff,tau,1.0e-6,zcan,cold,   &
     &                       tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
      enddo
      write(*,'(a,2e14.6)') " load after a long heavy snowfall, capacity: ",  &
     &                      zcan,zcap
      if (zcan > zcap*(1.0+1.0e-12) .or. zfcan /= 1.0) then
         write(*,*) " FAIL: the load left its capacity"
         kfail = kfail + 1
      endif
      end subroutine capacity

!     ABOVE FREEZING the canopy sheds its load on that step, and it reaches the
!     ground rather than disappearing.
      subroutine flushes(kfail)
      integer, intent(inout) :: kfail
      real :: zcan,zthr,zfcan,zfint,zheld
      integer :: j
      zcan = 0.0
      do j=1,200
         call canopysnowstep(0.8,3.1,capai,ceff,tau,1.0e-7,zcan,cold,   &
     &                       tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
      enddo
      zheld = zcan
      call canopysnowstep(0.8,3.1,capai,ceff,tau,0.0,zcan,tmelt+0.1,    &
     &                    tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
      write(*,'(a,2e14.6)') " load before the thaw, throughfall over it: ",  &
     &                      zheld,zthr*dt
      if (zcan /= 0.0 .or. abs(zthr*dt - zheld) > 1.0e-12*zheld) then
         write(*,*) " FAIL: the thaw did not put the load on the ground"
         kfail = kfail + 1
      endif
      end subroutine flushes

!     THE STEADY LOAD. Under a constant snowfall below the capacity the
!     relaxation has an equilibrium written down in advance: c*S_f*tau.
      subroutine approaches(kfail)
      integer, intent(inout) :: kfail
      integer :: j
      real :: zcan,zthr,zfcan,zfint,zsf,zeq
      zsf = 1.0e-9
      zeq = ceff*0.8*zsf*tau
      zcan = 0.0
      do j=1,20000
         call canopysnowstep(0.8,3.1,capai,ceff,tau,zsf,zcan,cold,      &
     &                       tmelt,dt,1.0,zcan,zthr,zfcan,zfint)
      enddo
      write(*,'(a,2e14.6)') " steady load, and c*S_f*tau: ",zcan,zeq
      if (abs(zcan - zeq) > 1.0e-12*zeq) then
         write(*,*) " FAIL: the load does not approach its own equilibrium"
         kfail = kfail + 1
      endif
!     And at equilibrium the ground receives every flake that falls.
      if (abs(zthr - zsf) > 1.0e-12*zsf) then
         write(*,*) " FAIL: a steady canopy does not pass on a steady snowfall"
         kfail = kfail + 1
      endif
      end subroutine approaches

      end program verify_canopy_snow
