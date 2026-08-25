!     The land column reduction check.
!
!     LSHY-3 replaces the scalar soil water bucket behind a SELECTABLE
!     reduction. That word carries the whole requirement: the bucket must stay
!     reachable, and it must be the EXACT reduction of the column that replaces
!     it, not a similar answer under simple conditions.
!
!     This program is the arm that can fail. It drives `bucket_step` and
!     `column_step` from `landcolumn.f90` -- the model's own kernels, not a
!     lifted copy -- with the same flux sequence over the same store, and
!     requires BITWISE equality of the store and the runoff at every step.
!     Not agreement to a tolerance: a reduction right to 1e-12 cannot tell a
!     refactor from a physics change, and this check exists precisely to tell
!     them apart.
!
!     THE NEGATIVE CONTROLS are the point of the program. Three configurations
!     that MUST break the equality are run beside the reduction, because a
!     comparison that cannot fail proves only that both sides were run:
!
!       control 1  two layers instead of one. The column now backs water up
!                  through a second store, so a step that would have overflowed
!                  the bucket does not always overflow the column.
!       control 2  free drainage instead of an impermeable base. Water that the
!                  bucket returned as surface runoff leaves downward instead,
!                  so the runoff differs and the store does not refill.
!       control 3  the wetness limiter's general form at the same knee with an
!                  exponent of two. The default path is a branch and not
!                  algebra, so this must move the answer.
!
!     A control that does NOT break the equality is a defect in this program,
!     and the program exits non-zero for it just as it does for a failed
!     reduction.
!
!     The flux sequence is deterministic and is built to visit every branch:
!     overfilling, emptying past zero, exact-capacity arrivals and long dry
!     spells. It is a fixed sequence and not a random one, so a failure is
!     reproducible from the seed printed in the header.

      program verify_land_column_reduction
      use landcolumn
      implicit none

      integer, parameter :: NSTEP = 20000
      real,    parameter :: DT    = 1800.0      ! s, a plausible model timestep
      real,    parameter :: WMAX  = 0.5         ! m, ExoPlaSim's uniform default

      real    :: zwb, zroffb, zwbn
      real    :: zwc(NLSOILWX), zcap(NLSOILWX), zwcn(NLSOILWX)
      real    :: zice(NLSOILWX)
      real    :: zroffc, zdrnc
      real    :: zflux
      integer :: jstep, jfail, jlay
      integer :: icontrol
      logical :: lreduction_ok
      logical :: lcontrol_broke(3)
      real    :: zwatres, zeneres, zwrongres
      character(len=48) :: yname(3)

      yname(1) = 'two layers instead of one'
      yname(2) = 'free drainage instead of impermeable'
      yname(3) = 'wetness limiter exponent two'

      write(*,*) 'land column reduction check'
      write(*,*) '  steps        ', NSTEP
      write(*,*) '  timestep (s) ', DT
      write(*,*) '  capacity (m) ', WMAX
      write(*,*) ''

!     ------------------------------------------------------------------
!     THE REDUCTION. One layer, impermeable base, layer capacity = WMAX.
!     ------------------------------------------------------------------

      zwb      = 0.5 * WMAX
      zwc(:)   = 0.0
      zice(:)  = 0.0
      zwc(1)   = zwb
      zcap(:)  = 0.0
      zcap(1)  = WMAX
      jfail    = 0

      do jstep = 1, NSTEP
       zflux = flux_at(jstep)
       call bucket_step(zwb, WMAX, zflux, DT, zwbn, zroffb)
       call column_step(1, zwc, zice, zcap, zflux, DT, LOWER_IMPERMEABLE,     &
     &                  zwcn, zroffc, zdrnc)
       if (zwbn /= zwcn(1) .or. zroffb /= zroffc .or. zdrnc /= 0.0) then
        if (jfail < 5) then
         write(*,'(a,i8)') '  MISMATCH at step ', jstep
         write(*,'(a,3e24.16)') '    flux, bucket w, column w  ',            &
     &                          zflux, zwbn, zwcn(1)
         write(*,'(a,3e24.16)') '    bucket roff, column roff, drain ',      &
     &                          zroffb, zroffc, zdrnc
        endif
        jfail = jfail + 1
       endif
       zwb    = zwbn
       zwc(:) = zwcn(:)
      enddo

      lreduction_ok = (jfail == 0)
      if (lreduction_ok) then
       write(*,*) '  reduction    EXACT over every step, bitwise'
      else
       write(*,'(a,i8,a)') '   reduction    FAILED on ', jfail, ' steps'
      endif

!     ------------------------------------------------------------------
!     THE NEGATIVE CONTROLS. Each must break the equality.
!     ------------------------------------------------------------------

      lcontrol_broke(:) = .false.

      do icontrol = 1, 2
       zwb     = 0.5 * WMAX
       zwc(:)  = 0.0
       zice(:) = 0.0
       zcap(:) = 0.0
       if (icontrol == 1) then
!       Two layers of half the capacity each. Same column capacity.
        zcap(1) = 0.5 * WMAX
        zcap(2) = 0.5 * WMAX
        zwc(1)  = 0.5 * WMAX
        zwc(2)  = 0.0
        zwb     = zwc(1) + zwc(2)
       else
        zcap(1) = WMAX
        zwc(1)  = zwb
       endif

       do jstep = 1, NSTEP
        zflux = flux_at(jstep)
        call bucket_step(zwb, WMAX, zflux, DT, zwbn, zroffb)
        if (icontrol == 1) then
         call column_step(2, zwc, zice, zcap, zflux, DT, LOWER_IMPERMEABLE,   &
     &                    zwcn, zroffc, zdrnc)
        else
         call column_step(1, zwc, zice, zcap, zflux, DT, LOWER_FREEDRAIN,     &
     &                    zwcn, zroffc, zdrnc)
        endif
        if (zroffb /= zroffc) lcontrol_broke(icontrol) = .true.
        zwb    = zwbn
        zwc(:) = zwcn(:)
       enddo
      enddo

!     Control 3 is on the wetness limiter rather than the store.

      lcontrol_broke(3) = .false.
      do jstep = 1, NSTEP
       zflux = 0.4 * WMAX * real(mod(jstep, 100)) / 100.0
       if (land_wetness(zflux, WMAX, 0.4, 0.0, 1)                            &
     &     /= land_wetness(zflux, WMAX, 0.4, 0.0, 2)) then
        lcontrol_broke(3) = .true.
       endif
      enddo

      write(*,*) ''
      do icontrol = 1, 3
       if (lcontrol_broke(icontrol)) then
        write(*,'(a,i1,a,a)') '   control ', icontrol,                        &
     &        ' broke the equality as required: ', trim(yname(icontrol))
       else
        write(*,'(a,i1,a,a)') '   control ', icontrol,                        &
     &        ' DID NOT break the equality: ', trim(yname(icontrol))
       endif
      enddo

!     ------------------------------------------------------------------
!     The wetness limiter's own reduction: the default parameters must
!     reproduce the active expression bitwise at every store.
!     ------------------------------------------------------------------

      jfail = 0
      do jstep = 0, NSTEP
       zflux = WMAX * real(jstep) / real(NSTEP)
       if (land_wetness(zflux, WMAX, 0.4, 0.0, 1)                            &
     &     /= AMIN1(1., zflux / (0.4 * WMAX))) then
        jfail = jfail + 1
       endif
      enddo
      write(*,*) ''
      if (jfail == 0) then
       write(*,*) '  wetness      default form reproduces the active '//     &
     &            'expression bitwise'
      else
       write(*,'(a,i8,a)') '   wetness      FAILED on ', jfail, ' stores'
       lreduction_ok = .false.
      endif

!     ------------------------------------------------------------------
!     LSHY-5. Freeze and thaw: water and energy must both close, and a
!     control that books the wrong latent heat must break the energy
!     identity while leaving the water identity untouched.
!     ------------------------------------------------------------------
!
!     THE TOLERANCE IS DERIVED, NOT CHOSEN. A number picked after seeing the
!     residual is not a criterion, and a fixed constant here would be one:
!     this program runs at two precisions and the same physics leaves residuals
!     three orders apart between them.
!
!     So each case carries its OWN bound, computed from the arithmetic. Water
!     is a sum of two stores and rounds at a few machine epsilons. Energy is
!     read through a temperature DIFFERENCE near 273 K, so its relative
!     precision is degraded by the ratio of the absolute temperature to the
!     change -- a case that moves the layer by a millikelvin cannot be checked
!     to better than eps times 273000. The routine returns the worst ratio of
!     a residual to its own derived bound, and the criterion is that no case
!     exceeds its bound. That is one number at both precisions.
!
!     Control 4 is reported the same way and must exceed its bound by orders
!     of magnitude, because a wrong latent heat is a factor of seven and a
!     half and not a rounding.

      write(*,*) ''
      call phase_fixtures(zwatres, zeneres, zwrongres)
      write(*,'(a,e12.4)') '   phase water  worst residual over bound  ',     &
     &                     zwatres
      write(*,'(a,e12.4)') '   phase energy worst residual over bound  ',     &
     &                     zeneres
      write(*,'(a,e12.4)') '   control 4    wrong latent heat over bound ',   &
     &                     zwrongres

      if (zwatres > 1.0) then
       write(*,*) '  phase water  FAILED: freeze and thaw do not conserve'
       lreduction_ok = .false.
      else
       write(*,*) '  phase water  closes to the arithmetic'
      endif
      if (zeneres > 1.0) then
       write(*,*) '  phase energy FAILED: the latent heat and the '//         &
     &            'temperature change disagree'
       lreduction_ok = .false.
      else
       write(*,*) '  phase energy closes on the declared latent heat'
      endif
      if (zwrongres <= 1.0) then
       write(*,*) '  control 4    DID NOT break the energy identity: '//      &
     &            'the check cannot see a wrong latent heat'
       lreduction_ok = .false.
      else
       write(*,*) '  control 4    broke the energy identity as required'
      endif

      if (.not. lreduction_ok) stop 1
      do icontrol = 1, 3
       if (.not. lcontrol_broke(icontrol)) stop 1
      enddo
      write(*,*) ''
      write(*,*) 'PASS'

      contains

!     ==========================
!     SUBROUTINE PHASE_FIXTURES
!     ==========================
!
!     Drives `phase_step` through a full freeze-thaw cycle in one layer and
!     returns the worst relative residual on each of the two identities, plus
!     the residual a deliberately wrong latent heat produces.
!
!     WATER:  liquid + ice is unchanged by the exchange.
!     ENERGY: the temperature change times the layer's heat capacity equals
!             the latent heat of the mass that changed phase. Those are the
!             same number read two ways, so an implementation that gets the
!             water right and the energy wrong is exactly what this separates.
!
!     Control 4 calls the same kernel with the latent heat of VAPORISATION
!     where fusion belongs, which is the mistake `landmod.f90`'s snowmelt
!     avoids by using (ALS - ALV) rather than either constant alone. It must
!     leave the water identity intact and break the energy one, because that
!     is precisely what a wrong latent heat does.

      subroutine phase_fixtures(pwater, penergy, pwrong)
      real, intent(out) :: pwater    ! worst water residual over its own bound
      real, intent(out) :: penergy   ! worst energy residual over its own bound
      real, intent(out) :: pwrong    ! the same, under the wrong latent heat
      real, parameter :: TMELT = 273.16
      real, parameter :: LFUS  = 2.8345e6 - 2.5008e6   ! als - alv, as landmod uses
      real, parameter :: LVAP  = 2.5008e6
      real, parameter :: RHOW  = 1000.0
      real, parameter :: CAP   = 2.4e6                 ! landmod soilcap
      real, parameter :: DZ    = 1.5                   ! the contract's column
      real :: zliq, zice0, ztem, zliqn, zicen, ztemn
      real :: zbefore, zafter, zmoved, zlat, zdte, zres, zeps, zbound
      integer :: jt

!     A few roundings on a two-store transfer. Eight is the count of floating
!     point operations between the input and the residual, not a fitted factor.
      zeps = EPSILON(1.0)

      pwater  = 0.
      penergy = 0.
      pwrong  = 0.

!     A sweep across temperature and across how much water is present, so the
!     energy-limited branch and the water-limited branch are both visited in
!     both directions.
      do jt = -300, 300
       ztem  = TMELT + real(jt) * 0.1
       zliq  = 0.05 + 0.002 * real(mod(jt + 300, 41))
       zice0 = 0.04 + 0.002 * real(mod(jt + 137, 37))

       call phase_step(zliq, zice0, ztem, CAP, DZ, TMELT, LFUS, RHOW,          &
     &                 zliqn, zicen, ztemn)

       zbefore = zliq + zice0
       zafter  = zliqn + zicen
       zres    = ABS(zafter - zbefore) / AMAX1(zbefore, 1.0e-30)
       pwater  = AMAX1(pwater, zres / (8.0 * zeps))

       zmoved = zicen - zice0
       zlat   = RHOW * ABS(zmoved) * LFUS
       zdte   = ABS(ztemn - ztem) * CAP * DZ
       if (zlat > 0.0 .and. ztemn /= ztem) then
!       The bound: the temperature difference is formed from two numbers near
!       273 K, so it carries eps of the ABSOLUTE temperature as an error and
!       that error is a fraction eps*T/|dT| of the difference itself.
        zbound  = 8.0 * zeps * ztem / ABS(ztemn - ztem)
        zres    = ABS(zdte - zlat) / zlat
        penergy = AMAX1(penergy, zres / zbound)
       endif

!      Control 4: the same exchange booked at the latent heat of vaporisation.
       call phase_step(zliq, zice0, ztem, CAP, DZ, TMELT, LVAP, RHOW,          &
     &                 zliqn, zicen, ztemn)
       zbefore = zliq + zice0
       zafter  = zliqn + zicen
       zres    = ABS(zafter - zbefore) / AMAX1(zbefore, 1.0e-30)
!      A wrong latent heat must NOT break the water identity. If it does, the
!      water check is not measuring what it claims to, so this residual is
!      folded into the water verdict rather than ignored.
       pwater = AMAX1(pwater, zres / (8.0 * zeps))
       zmoved = zicen - zice0
       zlat   = RHOW * ABS(zmoved) * LFUS
       zdte   = ABS(ztemn - ztem) * CAP * DZ
       if (zlat > 0.0 .and. ztemn /= ztem) then
        zbound = 8.0 * zeps * ztem / ABS(ztemn - ztem)
        pwrong = AMAX1(pwrong, ABS(zdte - zlat) / zlat / zbound)
       endif
      enddo

      return
      end subroutine phase_fixtures

!     The flux sequence. Deterministic, and built to visit every branch of
!     both kernels: sustained wetting past capacity, sustained drying past
!     empty, and arrivals that land exactly on a boundary.
      pure function flux_at(kstep) result(pflux)
      integer, intent(in) :: kstep
      real :: pflux
      integer :: iphase

      iphase = mod(kstep, 400)
      if (iphase < 120) then
!      wetting, hard enough to overflow
       pflux =  3.0e-6 * real(1 + mod(kstep, 7))
      else if (iphase < 200) then
!      exactly the capacity in one step, which is the boundary AMIN1 decides
       pflux = WMAX / DT
      else if (iphase < 360) then
!      drying, hard enough to take the store past empty
       pflux = -2.0e-6 * real(1 + mod(kstep, 5))
      else
!      a long dry spell at a rate that leaves a remainder
       pflux = -1.0e-9
      endif

      return
      end function flux_at

      end program verify_land_column_reduction
