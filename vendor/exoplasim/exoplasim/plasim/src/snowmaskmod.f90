!     ==================================================================
!     THE CANOPY/SNOW MASK, K*
!
!     How much of a snow-covered gridcell's ground the vegetation hides from
!     the shortwave, and how bright the thing doing the hiding is.
!
!     Upstream masks snow with one scalar. The land module blends a forested
!     snow endmember against the exposed one at the cell's canopy cover
!     fraction, and the endmember itself was built once at initialisation from
!     a fixed masking fraction. That makes the hidden fraction
!     dforest*forcov, a constant times the cover, with no dependence on how
!     deep the snow is against the canopy, on how much plant area the canopy
!     actually carries, or on whether the canopy is itself holding snow.
!
!     This routine is the assessed operator. It returns the same number under
!     its defaults -- exactly, to the last bit, which is the check that can
!     fail -- and generalises it along three axes that the vegetation
!     component can drive once it reports them:
!
!     BURIAL. Snow that reaches into the canopy hides the canopy rather than
!     the other way round. With the plant area distributed uniformly in the
!     vertical between the ground and the canopy top, the share still standing
!     above the snow surface is 1 - z_snow/h_canopy, and that is the exposed
!     fraction below. z_snow is metres OF SNOW, not of water equivalent; the
!     caller converts. A canopy height at or below zero means no burial, which
!     is the state of a world whose vegetation does not report a height.
!
!     PLANT AREA. A canopy hides the ground in proportion to how little of it
!     can be seen through. The gap fraction of a canopy carrying plant area
!     index L is exp(-k L) for an extinction coefficient k, so the masking
!     efficiency is 1 - exp(-k L). A plant area index at or below zero keeps
!     the two masking efficiencies the caller passes, which is what preserves
!     the one-scalar reduction.
!
!     This is the gap-fraction form Essery (2013) evaluates as the middle of
!     three, and the one that paper found gives large-scale answers close to a
!     two-stream canopy when both are given realistic cover and parameters. A
!     two-stream scheme is therefore a bracket on this rather than an assumed
!     improvement, which is why this is the operator and not that one. What is
!     taken from that source is the FORM and the extinction coefficient, never
!     an albedo: every albedo here is the caller's, computed under this world's
!     star.
!
!     Two efficiencies and not one, because the caller has two: masking is
!     stronger over bright snow than over aged snow, since light escaping the
!     snowpack upward can be intercepted on the way out. Their ratio is
!     carried across to the plant-area law rather than restated, so whatever
!     the caller means by that contrast is what this means by it.
!
!     INTERCEPTED SNOW. Snow held in the canopy is not snow on the ground: it
!     brightens the thing that is doing the masking instead of being masked.
!     The interception weight returned here moves the forested endmember
!     toward the exposed-snow endmember, so a canopy completely loaded with
!     snow gives the albedo of open snow whatever its cover is, which is the
!     limit that has to hold.
!
!     The endmember it moves toward is the caller's own exposed snow rather
!     than a fixed canopy-snow albedo, which is what the schemes Essery (2013)
!     surveys use. That keeps this world's star out of a number measured under
!     another one, and it makes the weight an EFFECTIVE optical fraction rather
!     than a geometric one: canopy snow is patchier and more shaded than the
!     snowpack, so the fraction that brightens the canopy like open snow is
!     smaller than the fraction of the canopy holding snow.
!
!     UNLOADING IS NOT REPRESENTED: the weight is a standing fraction the
!     caller declares, not a canopy snow store with a mass balance, and
!     predicting it needs one. See BIO-32.
!
!     The routine deliberately reads no model state. Everything it needs is an
!     argument, so exoplasim/scripts/verify_snow_mask.sh can compile it alone
!     and check it against answers written down in advance.
!
!     Worldbuilding frame: this is the surface albedo of a simulated planet's
!     snow-covered vegetated land, in a climate model of it.
!     ==================================================================

      module snowmaskmod
      implicit none

      contains

      pure subroutine snowcanopymask(pcover,psnowm,phcan,ppai,pext,     &
     &                               pcovmx,pcovmn,pint,                &
     &                               pfcovmx,pfcovmn,pfint,pkmx,pkmn)

      real, intent(in)  :: pcover  ! canopy cover fraction of the cell (0-1)
      real, intent(in)  :: psnowm  ! ground snow depth (m OF SNOW)
      real, intent(in)  :: phcan   ! canopy height (m); <= 0 means no burial
      real, intent(in)  :: ppai    ! plant area index; <= 0 keeps pcovmx/pcovmn
      real, intent(in)  :: pext    ! extinction coefficient of the gap fraction
      real, intent(in)  :: pcovmx  ! masking efficiency, bright-snow end (0-1)
      real, intent(in)  :: pcovmn  ! masking efficiency, aged-snow end (0-1)
      real, intent(in)  :: pint    ! intercepted-snow fraction of the canopy

      real, intent(out) :: pfcovmx ! cover to blend the bright endmember at
      real, intent(out) :: pfcovmn ! cover to blend the aged endmember at
      real, intent(out) :: pfint   ! interception weight, clamped to 0-1
      real, intent(out) :: pkmx    ! K*, the masked ground fraction, bright end
      real, intent(out) :: pkmn    ! K*, the masked ground fraction, aged end

      real :: zexp     ! share of the canopy still above the snow surface
      real :: zeffmx   ! masking efficiency of that share, bright end
      real :: zeffmn   ! masking efficiency of that share, aged end
      real :: zratmx   ! what the plant-area law does to the caller's efficiency
      real :: zratmn

!     Burial. Exactly 1.0 with no canopy height, which is what makes the
!     default reduction bit-for-bit and not merely close.

      zexp = 1.0
      if (phcan > 0.0) then
         zexp = 1.0 - psnowm/phcan
         zexp = max(0.0,min(1.0,zexp))
      endif

!     Masking efficiency. Same, exactly the caller's numbers with no plant
!     area index, so the ratios below are exactly one.

      zeffmx = pcovmx
      zeffmn = pcovmn
      if (ppai > 0.0) then
         zeffmn = 1.0 - exp(-pext*ppai)
         zeffmx = zeffmn
         if (pcovmn > 0.0) zeffmx = min(1.0,zeffmn*pcovmx/pcovmn)
      endif

      zratmx = 1.0
      zratmn = 1.0
      if (pcovmx > 0.0) zratmx = zeffmx/pcovmx
      if (pcovmn > 0.0) zratmn = zeffmn/pcovmn

!     The cover the land module blends its forested endmember at. It carries
!     the plant-area ratio because the endmember was built at the caller's
!     efficiency: blending at cover*ratio hides cover*ratio*pcovmx of the
!     ground, which is cover*zeffmx, the fraction the canopy really hides.

      pfcovmx = pcover*zexp*zratmx
      pfcovmn = pcover*zexp*zratmn

!     The canopy cannot hide more ground than the cell has.

      if (pcovmx > 0.0) pfcovmx = min(pfcovmx,1.0/pcovmx)
      if (pcovmn > 0.0) pfcovmn = min(pfcovmn,1.0/pcovmn)

      pfcovmx = max(0.0,pfcovmx)
      pfcovmn = max(0.0,pfcovmn)

      pfint = max(0.0,min(1.0,pint))

!     K* itself, reported rather than used: it is the hidden fraction, and the
!     caller works in covers because that is the blend it already has.

      pkmx = pfcovmx*pcovmx
      pkmn = pfcovmn*pcovmn

      return
      end subroutine snowcanopymask

      end module snowmaskmod
