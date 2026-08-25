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

      pure subroutine column_step(klay, pw, pice, pcap, pflux, pdt, klower,    &
     &                            pwnew, proff, pdrn)
      integer, intent(in)  :: klay              ! number of layers in use
      real,    intent(in)  :: pw(NLSOILWX)      ! layer liquid store before (m)
      real,    intent(in)  :: pice(NLSOILWX)    ! layer ice, m of water equivalent
      real,    intent(in)  :: pcap(NLSOILWX)    ! layer capacity (m)
      real,    intent(in)  :: pflux             ! net water flux to the surface (m/s)
      real,    intent(in)  :: pdt               ! timestep (s)
      integer, intent(in)  :: klower            ! lower boundary code
      real,    intent(out) :: pwnew(NLSOILWX)   ! layer liquid store after (m)
      real,    intent(out) :: proff             ! surface runoff (m/s)
      real,    intent(out) :: pdrn              ! drainage out of the base (m/s)
      real    :: zw(NLSOILWX)
      real    :: zcap(NLSOILWX)
      real    :: zexc
      integer :: jlay

      zw(:) = pw(:)

!     FROZEN PORE SPACE IS OCCUPIED PORE SPACE. Ice in a layer is not available
!     to hold liquid, so it comes off the capacity and a frozen layer fills and
!     overflows sooner. That is how ice impedes infiltration here: not by a
!     fitted impedance factor on a conductivity this column does not have, but
!     by taking up the room. `frozen_impedance` below is the conductivity form,
!     declared for the gradient-driven hypothesis that would need it.
!
!     At zero ice this is `pcap(j) - 0.0`, which is bitwise `pcap(j)`, so the
!     reduction to the bucket is untouched.
      do jlay = 1, klay
       zcap(jlay) = AMAX1(0., pcap(jlay) - pice(jlay))
      enddo

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
       zexc      = AMAX1(0., zw(jlay) - zcap(jlay))
       zw(jlay)  = AMIN1(zcap(jlay), zw(jlay))
       zw(jlay+1) = zw(jlay+1) + zexc
      enddo

!     The base.

      zexc     = AMAX1(0., zw(klay) - zcap(klay))
      zw(klay) = AMIN1(zcap(klay), zw(klay))
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
       zexc     = AMAX1(0., zw(jlay) - zcap(jlay))
       zw(jlay) = AMIN1(zcap(jlay), zw(jlay))
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

!     ====================
!     SUBROUTINE PHASE_STEP
!     ====================
!
!     Freeze and thaw in one soil layer, LSHY-5. Water and energy both close,
!     and they close for the same reason: the mass moved between the two stores
!     and the temperature change are the SAME number read two ways, so an
!     implementation that gets one right cannot get the other wrong.
!
!     WHAT DRIVES IT is the layer's cold or heat content relative to melting,
!     not a rate and not a threshold on air temperature. A layer at ptmelt - dT
!     can freeze at most the water whose latent heat of fusion would raise it
!     back to ptmelt, and a layer above melting can thaw at most the ice whose
!     fusion would lower it back. Both are limited by the water actually there.
!     That is what makes the exchange conservative rather than parameterised,
!     and it is the property the fixtures check.
!
!     THE MODEL HAS NO SOIL ICE TODAY. Its five soil temperature layers evolve
!     on a fixed heat capacity and conductivity, no water in the bucket is ever
!     frozen, and melt water therefore always infiltrates whatever the soil
!     temperature. LPJ-GUESS carries `Frac_ice` per layer and reduces available
!     liquid under freezing, so the two columns disagree about whether water
!     that reached the ground is liquid. This kernel is the climate side of
!     closing that; `hydrography/config/land_water_ledger.yaml` books it as
!     `soil_freezing` and `soil_thaw`.
!
!     LATENT HEAT: plfus, and never plv or pls. The model declares a latent
!     heat of vaporisation and one of sublimation and derives fusion as the
!     difference, so the caller passes the difference. A phase term carrying the
!     wrong latent heat conserves water and loses energy by a factor of seven
!     and a half, which is what the wrong-latent-heat fixture demonstrates.

      pure subroutine phase_step(pliq, pice, ptem, pcap, pdz, ptmelt, plfus,   &
     &                           prhow, pliqn, picen, ptemn)
      real, intent(in)  :: pliq    ! liquid water in the layer (m)
      real, intent(in)  :: pice    ! ice in the layer, m of water equivalent
      real, intent(in)  :: ptem    ! layer temperature (K)
      real, intent(in)  :: pcap    ! volumetric heat capacity (J/m3/K)
      real, intent(in)  :: pdz     ! layer thickness (m)
      real, intent(in)  :: ptmelt  ! melting point (K)
      real, intent(in)  :: plfus   ! latent heat of fusion (J/kg)
      real, intent(in)  :: prhow   ! water density (kg/m3)
      real, intent(out) :: pliqn   ! liquid after (m)
      real, intent(out) :: picen   ! ice after (m water equivalent)
      real, intent(out) :: ptemn   ! temperature after (K)
      real :: zheat, zcapdz, zmove, zlat

      pliqn  = pliq
      picen  = pice
      ptemn  = ptem
      zcapdz = pcap * pdz
      if (zcapdz <= 0.0 .or. plfus <= 0.0) return

!     Energy per unit area between the layer and melting. Positive is cold
!     content available to freeze liquid; negative is heat available to thaw.
      zheat = zcapdz * (ptmelt - ptem)

      if (zheat > 0.0 .and. pliq > 0.0) then
       zmove = AMIN1(pliq, zheat / (prhow * plfus))
       pliqn = pliq - zmove
       picen = pice + zmove
       zlat  = prhow * zmove * plfus
       ptemn = ptem + zlat / zcapdz
      else if (zheat < 0.0 .and. pice > 0.0) then
       zmove = AMIN1(pice, (-zheat) / (prhow * plfus))
       picen = pice - zmove
       pliqn = pliq + zmove
       zlat  = prhow * zmove * plfus
       ptemn = ptem - zlat / zcapdz
      endif

      return
      end subroutine phase_step

!     ==========================
!     FUNCTION FROZEN_IMPEDANCE
!     ==========================
!
!     The factor ice in the pore space multiplies a hydraulic conductivity by.
!     DECLARED HERE AND USED BY NOTHING YET, because the scheme that would use
!     it is the gradient-driven hypothesis and this column has no conductivity:
!     `pedology/config/land_column_properties.yaml` carries
!     `flow.saturated_conductivity` as undeclared and says why. The layered
!     column impedes flow the other way, by taking the ice out of the capacity
!     in `column_step`, which needs no conductivity at all.
!
!     The form is (1 - f_ice) ** kexp on the ice fraction of the pore space,
!     which is one at no ice and zero at a fully frozen pore, so the reduction
!     to the no-ice case is exact rather than asymptotic. kexp is an axis and
!     not a constant, on the same grounds as the evaporation limiter's: this is
!     a place two models would differ and the difference should be selectable.

      pure function frozen_impedance(pice, ppore, kexp) result(pfac)
      real,    intent(in) :: pice   ! ice, m of water equivalent
      real,    intent(in) :: ppore  ! pore volume, m of water equivalent
      integer, intent(in) :: kexp
      real :: pfac
      real :: zf

      if (ppore <= 0.0) then
       pfac = 0.0
       return
      endif
      zf   = AMIN1(1., AMAX1(0., pice / ppore))
      pfac = (1. - zf) ** kexp

      return
      end function frozen_impedance

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
