!     ==================================================================
!     THE CANOPY SNOW STORE
!
!     Snow held in a vegetated gridcell's canopy: how much of the falling
!     snow it catches, how fast it loses it again, and how much of the ground
!     snowpack's supply therefore arrives late.
!
!     WHAT IT REPLACES. snowmaskmod already distinguishes intercepted canopy
!     snow from ground snow -- its interception weight moves the forested snow
!     endmember toward the caller's exposed-snow endmember, so a fully loaded
!     canopy gives the albedo of open snow whatever its cover. What it did not
!     do was PREDICT that weight: the weight was a standing fraction the caller
!     declared, so interception and unloading were not represented and the
!     seasonal course of canopy snow was absent. This routine is the store that
!     predicts it, and it is a mass balance rather than a diagnostic, so the
!     canopy is a reservoir the falling snow passes through instead of a
!     surface property.
!
!     THE SCHEME. Essery (2013), Geophys. Res. Lett. 40, 5521-5525,
!     10.1002/grl.51008, Sect. 3.2, Eq. (3):
!
!         dS_c/dt = c S_f - S_c / tau
!
!     for a canopy snow load S_c and a snowfall rate S_f, with the load limited
!     to a capacity that scales with the plant area index, all of it removed
!     when the air is above freezing, and the canopy snow cover fraction taken
!     as the load over the capacity. That is the same source snowmaskmod takes
!     its gap-fraction form and its extinction coefficient from, and it is the
!     canopy-snow half of the scheme that paper evaluates against a two-stream
!     canopy and against observed seasonal albedo cycles.
!
!     The step below is the EXACT solution of that equation over the interval,
!     with the snowfall held constant across it, rather than a forward
!     difference: the removal is a relaxation with its own timescale, and an
!     explicit step would carry a stability limit the caller would have to
!     respect and would lose mass at the limit rather than refusing.
!
!     WHAT IT CONSERVES, AND WHAT THAT LETS A CHECK DO. Every path out of the
!     falling snow is accounted: what the canopy catches, what it releases as
!     the load relaxes, what it cannot hold above its capacity, and what it
!     sheds when the air is above freezing. Throughfall is the sum of the ones
!     that reach the ground, so over any sequence of steps the snow that fell
!     equals the throughfall plus the change in the store, exactly. That is an
!     identity with a right answer, and exoplasim/scripts/verify_canopy_snow.sh
!     asserts it over a driven sequence rather than inside one step.
!
!     WHAT IT DOES NOT CARRY. Canopy sublimation. Intercepted snow sits in the
!     air rather than on the ground and a real canopy loses part of its load to
!     the atmosphere rather than to the soil beneath it, but that needs a canopy
!     energy balance this model does not have. So the store conserves mass and
!     the snow it hands the ground pack is an UPPER BOUND: the modelled ground
!     pack may receive less than this and cannot receive more. The store also
!     carries no temperature of its own; it is released to the ground pack,
!     whose own thermodynamics then take it.
!
!     THE INERT REDUCTION. With no plant area index, or no declared capacity,
!     there is no canopy to hold anything: the throughfall is the snowfall to
!     the last bit, the store stays empty, and the interception weight is the
!     one the caller declares. That is the state this world is in until the
!     vegetation component reports a plant area index per gridcell, and it is
!     what exoplasim/scripts/verify_canopy_snow.sh checks bitwise.
!
!     THE WEIGHT IS NOT THE COVER. What this returns as `pfcan` is a GEOMETRIC
!     fraction: the share of the canopy carrying snow. snowmaskmod's
!     interception weight is an EFFECTIVE OPTICAL fraction against the caller's
!     own exposed-snow endmember, and canopy snow is patchier and more shaded
!     than a snowpack, so the optical fraction is the smaller of the two. The
!     conversion between them is the caller's `pintc`, declared rather than
!     assumed, and at zero the store runs and drives nothing optical -- which is
!     what keeps the albedo path a separate decision from the water path.
!
!     The routine deliberately reads no model state. Everything it needs is an
!     argument, so exoplasim/scripts/verify_canopy_snow.sh can compile it alone
!     and check it against answers written down in advance.
!
!     Worldbuilding frame: this is the snow held in the canopy of a simulated
!     planet's vegetated land, in a climate model of it.
!     ==================================================================

      module cansnowmod
      implicit none

      contains

      pure subroutine canopysnowstep(pcover,ppai,pcapai,pceff,ptau,     &
     &                               psnowf,pcansn,ptair,ptmelt,pdt,    &
     &                               pintc,                             &
     &                               pcansnew,pthrough,pfcan,pfint)

      real, intent(in)  :: pcover  ! canopy cover fraction of the cell (0-1)
      real, intent(in)  :: ppai    ! plant area index; <= 0 means no store
      real, intent(in)  :: pcapai  ! canopy snow capacity per unit plant area
                                   ! index (m water equivalent)
      real, intent(in)  :: pceff   ! interception efficiency (0-1)
      real, intent(in)  :: ptau    ! unloading timescale (s)
      real, intent(in)  :: psnowf  ! snowfall rate onto the cell (m w.e./s)
      real, intent(in)  :: pcansn  ! canopy snow at the start of the step
                                   ! (m water equivalent, cell mean)
      real, intent(in)  :: ptair   ! air temperature (K)
      real, intent(in)  :: ptmelt  ! melting point (K)
      real, intent(in)  :: pdt     ! length of the step (s)
      real, intent(in)  :: pintc   ! effective optical fraction per unit
                                   ! geometric canopy snow cover

      real, intent(out) :: pcansnew ! canopy snow at the end of the step
      real, intent(out) :: pthrough ! throughfall to the ground pack (m w.e./s)
      real, intent(out) :: pfcan    ! geometric canopy snow cover fraction
      real, intent(out) :: pfint    ! interception weight for snowcanopymask

      real :: zcov     ! cover, held inside its own range
      real :: zeff     ! efficiency, likewise: a fraction outside 0-1 would
                       ! make the throughfall negative rather than refusing
      real :: zcap     ! capacity of this cell's canopy (m w.e., cell mean)
      real :: zeq      ! the load this snowfall and this timescale hold (m w.e.)
      real :: znet     ! what the canopy takes from the falling snow, net of
                       ! what it gives back over the same step (m w.e.)
      real :: zover    ! what it cannot hold (m w.e.)
      real :: zflush   ! what it sheds above freezing (m w.e.)
      real :: zstore   ! the load through the step (m w.e.)

      zcov = max(0.0,min(1.0,pcover))
      zeff = max(0.0,min(1.0,pceff))
      zcap = 0.0
      if (ppai > 0.0 .and. pcapai > 0.0) zcap = zcov*ppai*pcapai

!     No canopy to hold anything. The store empties into the ground pack and
!     the throughfall is the snowfall itself, which with an empty store is the
!     snowfall to the last bit. This is the state the model ships in.

      if (zcap <= 0.0 .or. pdt <= 0.0 .or. ptau <= 0.0) then
         pcansnew = 0.0
         pthrough = psnowf
         if (pdt > 0.0) pthrough = psnowf + pcansn/pdt
         pfcan    = 0.0
         pfint    = 0.0
         return
      endif

!     Interception and unloading, together, as the exact solution of
!     dS/dt = c S_f - S/tau over the step with the snowfall held across it.
!     They are one expression and not two because they act at once: a load
!     relaxes while it is being added to, and splitting them into a full
!     interception followed by a decay of the OLD load leaves the steady state
!     high by half the ratio of the step to the timescale. The equilibrium
!     below is what the routine has to approach and is what
!     verify_canopy_snow.sh checks it against.

      zeq  = zeff*zcov*psnowf*ptau
      znet = (zeq - pcansn)*(1.0 - exp(-pdt/ptau))

      zstore = pcansn + znet

!     What the canopy cannot hold falls through.

      zover  = max(0.0,zstore - zcap)
      zstore = zstore - zover

!     Above freezing the canopy sheds its load: it melts and drips, or it
!     slides. Both reach the ground, which is why this is throughfall and not
!     a loss.

      zflush = 0.0
      if (ptair > ptmelt) then
         zflush = zstore
         zstore = 0.0
      endif

      pcansnew = zstore

!     Throughfall: what fell and the canopy did not keep, plus what it could
!     not hold and what it shed. Those are the only paths out of the falling
!     snow, so the snow that fell equals this plus the change in the store,
!     exactly and at every step.
!
!     It cannot go negative: the canopy's net gain over a step is at most
!     c*cover*S_f*tau*(1 - exp(-dt/tau)), which is at most S_f*dt for any
!     efficiency and cover inside their own range.

      pthrough = (psnowf*pdt - znet + zover + zflush)/pdt

      pfcan = min(1.0,zstore/zcap)
      pfint = max(0.0,min(1.0,pintc*pfcan))

      return
      end subroutine canopysnowstep

      end module cansnowmod
