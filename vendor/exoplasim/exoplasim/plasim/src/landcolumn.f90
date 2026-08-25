!     ===================
!     MODULE LANDCOLUMN
!     ===================
!
!     The land liquid-water kernels. Two schemes behind one interface, and the
!     scalar bucket is the EXACT reduction of the layered column at one layer
!     with an impermeable base.
!
!     WHY THIS IS A SEPARATE FILE. It depends on nothing -- no plasimmod, no
!     resmod, no grid -- so it compiles on its own and the reduction can be
!     verified by a standalone program rather than asserted in a comment.
!     `exoplasim/scripts/verify_land_column_reduction.f90` does exactly that:
!     it drives both kernels with the same flux sequence and requires bitwise
!     equality. A reduction that is only claimed is not a reduction.
!
!     BIT-EXACTNESS IS A DESIGN CONSTRAINT AND NOT AN ASPIRATION. The two
!     kernels below perform the same operations in the same order on the same
!     values at the reduction, which is why `column_step` uses
!     AMIN1(cap, w) to take a layer to its capacity rather than the
!     algebraically identical w - AMAX1(0., w - cap). The second is not the
!     same number in floating point, and a reduction that is right to a
!     tolerance cannot tell a refactor from a physics change.
!
!     WHAT THE COLUMN IS. A parsimonious multilayer store: the surface flux
!     enters the top layer, each layer fills to its capacity and passes its
!     excess to the one below, and whatever the base cannot take backs up
!     through the column and leaves as surface runoff. With free drainage the
!     base excess leaves downward instead, which is the lower boundary
!     `nlandwdrain` selects. It is a tipping-bucket cascade and NOT a
!     gradient-driven flow: there is no unsaturated conductivity, no matric
!     potential and no upward flux, and pedology/config/land_column_properties.yaml
!     carries all three as undeclared for exactly that reason. Registering it
!     as one hypothesis of three is the point; it is not a claim that it is
!     the right one.

      module landcolumn
      implicit none

!     Compile-time ceiling on the number of liquid water layers. Eight is
!     above PALADYN's five and below anything that would make the threadprivate
!     copy of the state significant against the working-set target: at NHOR
!     for one thread's latitude slice this is the same shape of cost as the
!     five soil temperature layers already carried.
      integer, parameter :: NLSOILWX = 8

!     Lower boundary codes for `column_step`.
      integer, parameter :: LOWER_IMPERMEABLE = 0
      integer, parameter :: LOWER_FREEDRAIN   = 1

      contains

!     =====================
!     SUBROUTINE BUCKET_STEP
!     =====================
!
!     The active scheme. One store, one capacity, and runoff is the store
!     overflowing. Transcribed unchanged from landmod's WANDR, which is what
!     makes it the thing the column is checked against rather than a second
!     opinion about what WANDR did.
!
!     pwnew comes back clipped at zero. That clip CREATES WATER when the net
!     flux takes the store past empty, while the latent heat that removed it
!     has already been paid to the atmosphere. It is preserved here because
!     preserving it is the whole point of a reduction, and it is booked in
!     hydrography/config/land_water_ledger.yaml as `bucket_floor_creation`
!     from the `unowned_source` boundary. WORLD-HSSB owns measuring it.

      pure subroutine bucket_step(pw, pwmax, pflux, pdt, pwnew, proff)
      real, intent(in)  :: pw        ! store before the step (m)
      real, intent(in)  :: pwmax     ! capacity (m)
      real, intent(in)  :: pflux     ! net water flux to the surface (m/s)
      real, intent(in)  :: pdt       ! timestep (s)
      real, intent(out) :: pwnew     ! store after the step (m)
      real, intent(out) :: proff     ! surface runoff (m/s)
      real :: zw

      zw    = pw + pdt * pflux
      proff = AMAX1(0., zw - pwmax) / pdt
      pwnew = AMAX1(AMIN1(pwmax, zw), 0.)

      return
      end subroutine bucket_step

!     =====================
!     SUBROUTINE COLUMN_STEP
!     =====================
!
!     The layered column. `pcap` is the capacity of each layer, in metres of
!     water, and its sum is the column capacity. `klay` layers, `klower` the
!     lower boundary.
!
!     THE REDUCTION. At klay = 1, klower = LOWER_IMPERMEABLE and
!     pcap(1) = pwmax, this executes
!
!         zexc = AMAX1(0., zw - pcap(1))
!         zw   = AMIN1(pcap(1), zw)
!         proff = zexc / pdt
!         pwnew(1) = AMAX1(zw, 0.)
!
!     which is `bucket_step` operation for operation. Both loops are empty at
!     one layer, so nothing else runs and nothing else can differ.

      pure subroutine column_step(klay, pw, pcap, pflux, pdt, klower,          &
     &                            pwnew, proff, pdrn)
      integer, intent(in)  :: klay              ! number of layers in use
      real,    intent(in)  :: pw(NLSOILWX)      ! layer store before (m)
      real,    intent(in)  :: pcap(NLSOILWX)    ! layer capacity (m)
      real,    intent(in)  :: pflux             ! net water flux to the surface (m/s)
      real,    intent(in)  :: pdt               ! timestep (s)
      integer, intent(in)  :: klower            ! lower boundary code
      real,    intent(out) :: pwnew(NLSOILWX)   ! layer store after (m)
      real,    intent(out) :: proff             ! surface runoff (m/s)
      real,    intent(out) :: pdrn              ! drainage out of the base (m/s)
      real    :: zw(NLSOILWX)
      real    :: zexc
      integer :: jlay

      zw(:) = pw(:)

!     The surface flux enters the top layer. There is no infiltration capacity
!     here, so this scheme generates saturation-excess runoff only, exactly as
!     the bucket does. Infiltration excess is `flow.infiltration_capacity` in
!     the property contract and is undeclared.

      zw(1) = zw(1) + pdt * pflux

!     A withdrawal larger than the top layer holds draws on the layers below
!     it, because that is where the water is. The bucket cannot do this -- it
!     has one store and floors it -- so this is the one place the column is
!     strictly better rather than merely different, and it costs the reduction
!     nothing: at one layer the loop does not run and the floor below fires
!     exactly as the bucket's does.

      do jlay = 1, klay - 1
       if (zw(jlay) < 0.) then
        zw(jlay+1) = zw(jlay+1) + zw(jlay)
        zw(jlay)   = 0.
       endif
      enddo

!     Downward cascade: each layer fills to capacity and passes the rest on.

      do jlay = 1, klay - 1
       zexc      = AMAX1(0., zw(jlay) - pcap(jlay))
       zw(jlay)  = AMIN1(pcap(jlay), zw(jlay))
       zw(jlay+1) = zw(jlay+1) + zexc
      enddo

!     The base.

      zexc     = AMAX1(0., zw(klay) - pcap(klay))
      zw(klay) = AMIN1(pcap(klay), zw(klay))
      if (klower == LOWER_FREEDRAIN) then
       pdrn = zexc / pdt
       zexc = 0.
      else
       pdrn = 0.
      endif

!     Whatever the base could not take backs up through the column. At one
!     layer this loop does not run, which is what keeps the reduction exact.

      do jlay = klay - 1, 1, -1
       zw(jlay) = zw(jlay) + zexc
       zexc     = AMAX1(0., zw(jlay) - pcap(jlay))
       zw(jlay) = AMIN1(pcap(jlay), zw(jlay))
      enddo

      proff = zexc / pdt

!     The same floor at zero the bucket applies, layer by layer, and the same
!     unowned crossing with it.

      pwnew(:) = 0.
      do jlay = 1, klay
       pwnew(jlay) = AMAX1(zw(jlay), 0.)
      enddo

      return
      end subroutine column_step

!     ======================
!     FUNCTION LAND_WETNESS
!     ======================
!
!     The surface wetness factor the evaporation is throttled by. ONE
!     definition, where there were three copies of one expression in landmod.
!
!     THE GENERAL FORM AND ITS BRACKET. Both active limiters in the two models
!     this project has read sit on the same three axes,
!
!         beta = min(1, max((theta - theta_low)/(theta_high - theta_low), 0))**c
!
!     with theta the store as a fraction of capacity. This model is c = 1 with
!     the knee at drhsfull = 0.4 of capacity, so it evaporates at the potential
!     rate anywhere above two fifths full. cGENIE's ENTS is c = 4 with the knee
!     at a full store. At equal fractional fill the two differ by a factor of
!     39 at 0.4, 16 at 0.5, 8 at 0.6 and 2 at 0.8, agreeing only when full.
!     Neither is right and this is a BRACKET rather than an adoption; carrying
!     the general form makes it a runtime selection instead of a code fork.
!     Compare the SHAPE only: the capacities differ too, this model's being a
!     uniform scalar against ENTS's carbon-dependent one.
!
!     theta_low is the residual water fraction and both models pin it at zero,
!     so each approaches zero flux only as its store approaches empty and water
!     held below the wilting point is still available to be removed. That is
!     `states.residual` in the property contract, and it is undeclared.
!
!     DO NOT CITE Egea et al. (2011) for this form. That paper's beta
!     multiplies PHOTOSYNTHESIS, is applied to assimilation and deliberately
!     not re-applied in stomatal conductance, and is evaluated on root-zone
!     averaged moisture weighted by a root distribution. This beta multiplies
!     the surface evaporation flux over vegetated and bare ground alike, from a
!     store with no depth. Same algebra, different flux.
!
!     THE DEFAULT PATH IS THE ORIGINAL EXPRESSION, BY BRANCH AND NOT BY
!     ALGEBRA. `x**1.0` with a real exponent is not bitwise `x`, so the general
!     form cannot be used to reproduce the current one; the branch is what
!     keeps the default bit-identical.

      pure function land_wetness(pw, pwmax, pfull, plow, kexp) result(prhs)
      real,    intent(in) :: pw      ! store (m)
      real,    intent(in) :: pwmax   ! capacity (m)
      real,    intent(in) :: pfull   ! theta_high, fraction of capacity
      real,    intent(in) :: plow    ! theta_low, fraction of capacity
      integer, intent(in) :: kexp    ! the exponent c, an integer so **1 is exact
      real :: prhs
      real :: zx

      if (pwmax <= 0.0) then
       prhs = 1.0
       return
      endif

      if (plow == 0.0 .and. kexp == 1) then
!      The active form, unchanged.
       prhs = AMIN1(1., pw / (pfull * pwmax))
      else
       zx   = (pw / pwmax - plow) / (pfull - plow)
       prhs = AMIN1(1., AMAX1(0., zx)) ** kexp
      endif

      return
      end function land_wetness

      end module landcolumn
