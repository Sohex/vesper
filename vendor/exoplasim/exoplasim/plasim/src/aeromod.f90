!
!     **********************************************************
!     *  Parameters and subroutines for aerosol transport       *
!     **********************************************************

      module aeromod

!     **********************************************************
!     * This module contains the parameters, arrays and        *
!     * subroutines that are needed for transporting aerosols   *
!     * using the Flux-Form Semi-Lagrangian (FFSL) algorithm   *
!     * developed by S.-J. Lin (now at GFDL).                  *
!     **********************************************************
!     * The original transport code (in F77) was written by    *
!     * S.-J. Lin. Adaptation for the Planet Simulator was     *
!     * done by Hui Wan (MPI-M).
!     **********************************************************

      use pumamod

      logical,parameter :: aero_debug  = .FALSE.
      logical,parameter :: aero_zcross = .TRUE.
      logical,parameter :: aero_deform = .FALSE.

      logical,parameter :: aero_fill = .FALSE.
      logical,parameter :: aero_mfct = .FALSE.

      integer,parameter :: aero_iord = 2
      integer,parameter :: aero_jord = 2
      integer,parameter :: aero_kord = 3
      
      integer :: l_source = 1 ! 1 = photochemical haze (source at top level)
                                      ! 2 = dust (source at bottom level)
      integer :: l_bulk = 1 ! 1 = N2 atmosphere, 
                            ! 2 = H2 atmosphere

      integer,parameter :: aero_cnst = 1   ! 1 = constant preserving
                                           ! 2 = mass conserving

!      integer,parameter :: aero_j1  = 2  ! 1st lat. outside polar cap
!      integer,parameter :: aero_j2  = NLAT + 1 - aero_j1 
                                         ! last lat. outside polar cap 
      real :: apart = 50e-9 ! Radius of aerosol particle in m - DECLARED IN RADMOD AS WELL
      real :: rhop = 1000 ! Density of aerosol particle in kg/m3
      real :: fcoeff = 10e-13 ! Haze mass mixing ratio in kg/kg

!     Removal. Both terms are OFF by default, so one executable can run both
!     arms of an A/B and differ only by a namelist key. That is WORKFLOW.md
!     A3 and it is not a style preference: the low-I/O patch changes the
!     restart layout, so two arms built either side of a rebuild cannot share
!     a restart at all, and a term that cannot be switched off cannot be
!     tested.
!
!     ldepvel  0 = the legacy bottom-level sink, mmr = mmr*0.01 every step.
!                  It is a stability hack rather than a parameterisation: the
!                  implied deposition velocity is two orders of magnitude
!                  above any measured one, and because it does not scale with
!                  the timestep the implied velocity is inversely
!                  proportional to it. A rate that depends on the integration
!                  step is not a rate.
!              1 = a deposition velocity. The bottom layer keeps
!                  exp(-vdaero*dt/dz) per step, with dz built from the
!                  layer's own pressure thickness and gas density, so the
!                  rate is a property of the atmosphere and not of the step.
!     vdaero   the NON-gravitational dry deposition velocity, m/s: turbulent
!              transfer and impaction. Sedimentation to the surface is
!              already carried by the settling flux leaving the bottom layer,
!              so putting a settling term here as well would count it twice.
!              There is no default: aero_ini refuses ldepvel = 1 with vdaero
!              at or below zero rather than invent one, because the value
!              belongs with the rest of the removal budget in
!              aeolian/config/dust.yaml.
!     lwetdep  0 = no wet removal, which is what ExoPlaSim has and what the
!                  LMD Generic PCM has as well.
!              1 = below-cloud scavenging at Lambda = scava * p**scavb, with
!                  p the total precipitation rate in mm/h. That is the
!                  Sportisse (2007) form the offline chain in aeolian/
!                  already uses, so the two chains agree in form and not only
!                  in magnitude.
!     scava    that A, in s-1 per (mm/h)**B. No default, same reasoning as
!              vdaero; aeolian/config/dust.yaml carries it and its bracket.
!     scavb    that B, dimensionless. No default.
      integer :: ldepvel = 0
      integer :: lwetdep = 0
      real :: vdaero = 0.0
      real :: scava  = 0.0
      real :: scavb  = 0.0

!     Wind-driven dust emission, Kok et al. (2014) equation 18. OFF by default,
!     for the same WORKFLOW.md A3 reason the removal switches are: both arms of
!     an A/B have to come from one executable.
!
!     WHAT IS IN THE MODEL AND WHAT IS NOT. The per-cell source map stays
!     OUTSIDE: which ground can emit is a question about lithology, the lake
!     solution and the soil, and answering it needs the mesh export that the
!     climate model has never seen. Three surface boundary fields carry the
!     answer in, and everything time-varying -- friction velocity, soil water,
!     snow -- comes from the model's own state. That is the same boundary this
!     project already draws for albedo, roughness and soil-water capacity.
!
!       1801  dsrcw   erodible fraction x clipped clay fraction. Both are
!                     LINEAR prefactors on the flux, so their product is
!                     sufficient and a second field would be redundant.
!       1802  ddrage  the Marticorena-Bergametti (1995) drag partition, the
!                     fraction of the stress that reaches the erodible bed.
!                     Separate because it scales u* itself, so it sits inside
!                     the threshold comparison and inside the exponent and
!                     cannot be commuted out to a prefactor.
!       1803  dwpr    the Fecan et al. (1999) residual moisture w', percent.
!                     Separate because the moisture gate also needs the
!                     model's own soil water, which is already in memory.
!
!     ldustemit  0 = the legacy source, mmr(bottom) = fcoeff*land, which emits
!                    as much from forest as from salt pan and SETS rather than
!                    adds. This is the default and it reproduces the unpatched
!                    model exactly.
!                1 = Kok 18 evaluated over a Weibull distribution of subgrid
!                    wind and ADDED as a flux: d(mmr) = F g dt / dp.
!
!     Everything below is a calibration and every one of them lives in
!     aeolian/config/dust.yaml. None has a Fortran default, for the reason
!     vdaero has none: a plausible-looking number here would be a fourth place
!     for a constant that one file already owns, and it would go stale in
!     silence. aero_ini ABORTS on ldustemit = 1 with any of them unset.
!
!     dustz0   aeolian roughness of the erodible bed itself, m. NOT the grid
!              cell roughness: MB95 partitions stress between a bed and
!              centimetre-scale roughness ELEMENTS, and a 15 km orographic
!              variance is not a roughness element. dust.yaml carries the
!              argument and the bracket.
!     dustust0 standardized threshold friction velocity of an optimally
!              erodible bed, m/s, ALREADY SCALED for this planet's gravity.
!              The scaling is the fourth root of the gravity ratio and it is
!              applied outside, where config/planet.yaml's gravity lives.
!     dustustt the same for a typical erodible soil, m/s, likewise scaled.
!     dustra0  the standardization air density, kg/m3. Kok equation 6.
!     dustcd0  Cd0 of equation 18b, dimensionless.
!     dustce   Ce of equation 18b, dimensionless.
!     dustca   C_alpha of equation 18a, dimensionless.
!     dustwk   Weibull shape of the subgrid wind distribution. NOT optional:
!              emission goes as roughly u* cubed above a threshold, so the
!              flux of the mean and the mean of the flux differ by orders of
!              magnitude. Measured rather than declared; see dust.yaml.
!     dustnq   quadrature points across that distribution.
!     dustfa   A of Fecan equation 14.
!     dustfb   b of Fecan equation 14.
!     dustsnd  snow depth in m above which a cell emits nothing.
!     dustwcv  converts the model's soil water, a DEPTH in m, into gravimetric
!              soil moisture in percent. It is 100000/(soil depth x bulk
!              density) and both of those are declared in dust.yaml.
      integer :: ldustemit = 0
      integer :: dustnq  = 0
      real :: dustz0   = 0.0
      real :: dustust0 = 0.0
      real :: dustustt = 0.0
      real :: dustra0  = 0.0
      real :: dustcd0  = 0.0
      real :: dustce   = 0.0
      real :: dustca   = 0.0
      real :: dustwk   = 0.0
      real :: dustfa   = 0.0
      real :: dustfb   = 0.0
      real :: dustsnd  = 0.0
      real :: dustwcv  = 0.0

!     The three boundary fields, distributed as landmod holds its own surface
!     fields, and gathered ONCE into the global arrays aerocore works on.
!     aerocore runs serial on NROOT with the whole global grid, so a per-step
!     gather of a field that never changes would be pure cost.
      real :: dsrcw(NHOR)  = 0.0   ! code 1801, distributed
      real :: ddrage(NHOR) = 0.0   ! code 1802, distributed
      real :: dwpr(NHOR)   = 0.0   ! code 1803, distributed
      real :: gsrcw(NLON,NLAT)  = 0.0  ! the same, gathered and flipped
      real :: gdrage(NLON,NLAT) = 0.0
      real :: gwpr(NLON,NLAT)   = 0.0

      end module aeromod

!     ==================
!     SUBROUTINE AERO_INI
!     ==================

      subroutine aero_ini
      use aeromod
      use radmod, only: l_aerorad, aerofile, rad_apart => apart
      
      namelist/aero_nl/l_source,l_bulk,apart,rhop,fcoeff,l_aerorad,aerofile  &
     &                ,ldepvel,vdaero,lwetdep,scava,scavb                    &
     &                ,ldustemit,dustz0,dustust0,dustustt,dustra0            &
     &                ,dustcd0,dustce,dustca,dustwk,dustnq                   &
     &                ,dustfa,dustfb,dustsnd,dustwcv

      if (mypid==NROOT) then
         open(11,file=aero_namelist)
         read(11,aero_nl)
         close(11)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * AEROMOD ",a34)')
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist AERO_NL from <aero_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,aero_nl)
!
!        Refuse a switch that is on without the coefficient it needs, rather
!        than carry a plausible-looking default that nothing in this project
!        decided. Both values live in aeolian/config/dust.yaml.
!
         if (ldepvel == 1 .and. vdaero <= 0.0) then
            write(nud,*) '* ldepvel = 1 needs a positive vdaero in m/s.'
            write(nud,*) '* It is the non-gravitational dry deposition'
            write(nud,*) '* velocity; settling is already carried by the'
            write(nud,*) '* sedimentation flux. Take it from the removal'
            write(nud,*) '* block of aeolian/config/dust.yaml.'
            call mpabort('aero_nl: ldepvel = 1 without vdaero')
         endif
         if (lwetdep == 1 .and. (scava <= 0.0 .or. scavb <= 0.0)) then
            write(nud,*) '* lwetdep = 1 needs positive scava and scavb.'
            write(nud,*) '* Lambda = scava * p**scavb with p in mm/h.'
            write(nud,*) '* Take both from the removal block of'
            write(nud,*) '* aeolian/config/dust.yaml, which also carries the'
            write(nud,*) '* bracket on scava.'
            call mpabort('aero_nl: lwetdep = 1 without scava and scavb')
         endif
!
!        The emission scheme's calibration has no Fortran defaults either, and
!        for the same reason: aeolian/config/dust.yaml owns every one of these
!        and a number invented here would be a fourth place for it to go stale.
!        The abort is the loud failure and it is the one to want -- a silently
!        zeroed threshold would emit from the whole planet at every wind.
!
         if (ldustemit == 1) then
            if (dustz0   <= 0.0 .or. dustust0 <= 0.0 .or.                    &
     &          dustustt <= 0.0 .or. dustra0  <= 0.0 .or.                    &
     &          dustcd0  <= 0.0 .or. dustce   <= 0.0 .or.                    &
     &          dustca   <= 0.0 .or. dustwk   <= 0.0 .or.                    &
     &          dustnq   <= 0   .or. dustfa   <= 0.0 .or.                    &
     &          dustfb   <= 0.0 .or. dustsnd  <= 0.0 .or.                    &
     &          dustwcv  <= 0.0) then
               write(nud,*) '* ldustemit = 1 needs the whole Kok (2014)'
               write(nud,*) '* calibration: dustz0, dustust0, dustustt,'
               write(nud,*) '* dustra0, dustcd0, dustce, dustca, dustwk,'
               write(nud,*) '* dustnq, dustfa, dustfb, dustsnd, dustwcv.'
               write(nud,*) '* All of them live in aeolian/config/dust.yaml'
               write(nud,*) '* and are written into aero_namelist by'
               write(nud,*) '* exoplasim/scripts/run_exoplasim.py from the'
               write(nud,*) '* provenance file beside the 1801-1803 fields.'
               call mpabort('aero_nl: ldustemit = 1 without its calibration')
            endif
            if (dustustt < dustust0) then
               write(nud,*) '* dustustt is below dustust0. u*st0 is the'
               write(nud,*) '* OPTIMALLY erodible bed and is the minimum of'
               write(nud,*) '* the threshold curve, so nothing can sit under'
               write(nud,*) '* it. Kok equation 18b would return Cd above Cd0'
               write(nud,*) '* and the exponent would change sign.'
               call mpabort('aero_nl: dustustt < dustust0')
            endif
         endif
!
!        UPSTREAM DEFECT. radmod declares its own `apart` and this namelist
!        sets aeromod's, so the transport used the radius that was asked for
!        while the radiation kept the 50 nm photochemical-haze default. At
!        fixed number density the optical depth goes as apart squared, so the
!        shortwave aerosol came out (50e-9/apart)**2 of intent -- 1/385 at the
!        optical effective radius of this world's dust and 1/1948 at the
!        burden-matched one. Hand the value across; radini broadcasts it,
!        which is why this can sit inside the NROOT block. The two variables
!        stay separate because radmod cannot use aeromod: aeromod already uses
!        radmod, and make_plasim compiles it second.
!
         rad_apart = apart
      endif
      
      return
      end subroutine aero_ini

!     ===================
!     SUBROUTINE AERO_SURF
!     ===================

      subroutine aero_surf

!     Read the three dust-source boundary fields and put them where aerocore
!     can see them. Called from plasim.f90 by EVERY rank, which is the whole
!     reason this is not folded into aero_ini: aero_ini runs inside an
!     `if (mypid == NROOT)` block, and mpsurfgp is collective -- it broadcasts
!     the read flag and scatters the field, so a rank that does not enter it
!     hangs.
!
!     Once, not per step. aerocore runs serial on NROOT over the whole global
!     grid, so it needs these gathered; they are pure functions of terrain,
!     lithology, the lake solution and the soil and never change, so gathering
!     them every timestep would be cost for nothing.
!
!     THE LATITUDE FLIP IS NOT OPTIONAL. mpgagp returns a field in the MODEL's
!     latitude order, while plasim.f90 hands aerocore its tracer array flipped
!     south-to-north (`daeros(:,NLAT+1-jlat,:,1) = zmmr(:,jlat,:)`). Upstream
!     omitted the flip on the land mask and the solar zenith angle and drove
!     the aerosol source in the mirror hemisphere; that is defect 7 of
!     exoplasim-3.4.2-aerocore-defects.patch and this is the same convention,
!     written the same way round, on purpose.

      use aeromod
      implicit none

      real :: zgath(NLON,NLAT,1)
      real :: zmax
      integer :: j

!     ldustemit is read by aero_ini on NROOT only. Everything downstream of it
!     runs on NROOT too, but THIS routine does not: mpsurfgp is collective, so
!     every rank has to agree about whether it is called at all.

      call mpbci(ldustemit)
      if (ldustemit /= 1) return

      dsrcw(:)  = 0.0
      ddrage(:) = 0.0
      dwpr(:)   = 0.0
      call mpsurfgp('dsrcw' ,dsrcw ,NHOR,1)
      call mpsurfgp('ddrage',ddrage,NHOR,1)
      call mpsurfgp('dwpr'  ,dwpr  ,NHOR,1)

      call mpmaxval(dsrcw,NHOR,1,zmax)
      if (zmax <= 0.0) then
         if (mypid == NROOT) then
            write(nud,*) 'DUST EMISSION: no source field was read.'
            write(nud,*) 'Expected surface codes 1801, 1802 and 1803 in the'
            write(nud,*) 'run directory, written by'
            write(nud,*) 'aeolian/scripts/build_dust_source_fields.py.'
            write(nud,*) 'surfmod skips a missing surface file in silence, so'
            write(nud,*) 'the alternative to this abort is a run that emits'
            write(nud,*) 'nothing and says nothing.'
         endif
         call mpabort('ldustemit=1 but surface code 1801 is absent or zero')
      endif

      call mpgagp(zgath,dsrcw,1)
      do j = 1 , NLAT
         gsrcw(:,NLAT+1-j) = zgath(:,j,1)
      end do
      call mpgagp(zgath,ddrage,1)
      do j = 1 , NLAT
         gdrage(:,NLAT+1-j) = zgath(:,j,1)
      end do
      call mpgagp(zgath,dwpr,1)
      do j = 1 , NLAT
         gwpr(:,NLAT+1-j) = zgath(:,j,1)
      end do

      if (mypid == NROOT) then
         write(nud,'(/," *********************************************")')
         write(nud,'(" * DUST EMISSION (codes 1801-1803) is ON     *")')
         write(nud,'(" *********************************************")')
         write(nud,*) 'max erodible x clay prefactor  ',maxval(gsrcw)
         write(nud,*) 'drag partition, min and max    ',minval(gdrage),       &
     &                                                  maxval(gdrage)
         write(nud,*) 'Fecan w (percent), max         ',maxval(gwpr)
         write(nud,*) 'aeolian roughness (m)          ',dustz0
         write(nud,*) 'u*st0, u*st typical (m/s)      ',dustust0,dustustt
         write(nud,*) 'Weibull shape, quadrature      ',dustwk,dustnq
      endif

      return
      end subroutine aero_surf
  
  
!     ======================
!     SUBROUTINE AERO_MAIN
!     ======================

      subroutine aero_main

      use pumamod, only: du,dv,dp,du0,dv0,dp0,daeros,numrhos, &
                         NLON,NLAT,NLEV,NAERO,NHOR,   &
                         mypid,NROOT,sigmah,dt,dls,dswfl,dprl,dprc, &
                         dsnow,dwatc
      use tracermod
      use aeromod
      use radmod, only: gmu0, l_aerorad ! Use cosine of solar zenith angle from radmod;

      implicit none

      real :: zu   (NLON,NLAT,NLEV)
      real :: zv   (NLON,NLAT,NLEV)
      real :: zps0 (NLON,NLAT)
      real :: zps1 (NLON,NLAT)

      real :: x (NLON+1,NLAT,NLEV,NAERO)  ! for GUI output
      real :: y (NLON+1,NLAT,NLEV)         ! for GUI output
      real ::   angle(NLON,NLAT) ! Array for cosine of solar zenith angle
      real ::   aerosw(NLON,NLAT,NLEV) ! Array for SW flux 
      real ::   land(NLON,NLAT) ! Array for binary land mask

!     Gather buffer. mpgagp returns a field in the MODEL's latitude order,
!     while daeros and numrhos are handed to aerocore flipped south-to-north
!     (see plasim.f90, `daeros(:,NLAT+1-jlat,:,1) = zmmr(:,jlat,:)`). Every
!     field gathered here therefore has to be flipped the same way before it
!     is used against them. Upstream did not, so the land mask and the solar
!     zenith angle drove the aerosol source in the wrong hemisphere.
      real ::   zgath(NLON,NLAT,NLEV)

      real ::   prec(NLON,NLAT) ! Total precipitation rate (m/s), for lwetdep
      real ::   zprec(NHOR)     ! the same before gathering

      real ::   snow(NLON,NLAT)  ! Snow depth (m), for the emission gate
      real ::   wsoil(NLON,NLAT) ! Soil water (m), for the Fecan threshold

      integer :: j,jc

      character(len=9) :: aero_name

!     --- 

      call prepare_uvps( zu,zv,zps0,zps1,      & ! output
                         du0,dv0,dp0,du,dv,dp)   ! input

      if (l_aerorad == 0) then ! No radiative transfer
       select case (l_source) ! Choose your aerosol source
       case(1) ! Case 1: photochemical haze
         call solang ! Use subroutine from radmod to calculate solar zenith angle
         call mpgagp(zgath,gmu0,1) ! Gather from nodes
         do j=1,NLAT
            angle(:,NLAT+1-j) = zgath(:,j,1)
         end do
       case(2) ! Case 2: dust
         call mpgagp(zgath,dls,1) ! Import land-sea mask from landmod and reshape to match grid size
         do j=1,NLAT
            land(:,NLAT+1-j) = zgath(:,j,1)
         end do
       end select
      end if
      
      if (l_aerorad == 1) then ! Include radiative transfer
       select case (l_source) ! Choose aerosol source
       case(1) ! Case 1: photochemical haze     
        call mpgagp(zgath,dswfl,NLEV) ! Gather SW flux from nodes
        do j=1,NLAT
           aerosw(:,NLAT+1-j,:) = zgath(:,j,:)
        end do
       case(2) ! Case 2: dust
        call mpgagp(zgath,dls,1) ! Import land-sea mask from landmod and reshape to match grid size
        do j=1,NLAT
           land(:,NLAT+1-j) = zgath(:,j,1)
        end do
       end select
      end if 

!     Precipitation for the wet-scavenging term, large scale plus convective,
!     in m/s. Gathered only when the term is on, so that lwetdep = 0 costs
!     nothing and reproduces the unpatched model exactly. Flipped in latitude
!     like every other field gathered here.

      prec(:,:) = 0.0
      if (lwetdep == 1) then
         zprec(:) = dprl(:) + dprc(:)
         call mpgagp(zgath,zprec,1)
         do j=1,NLAT
            prec(:,NLAT+1-j) = zgath(:,j,1)
         end do
      end if

!     Snow depth and soil water for the emission scheme, gathered only when it
!     is on so that ldustemit = 0 costs nothing and reproduces the unpatched
!     model exactly. These two are the whole of what the source term takes from
!     the model's evolving state that is not already in aerocore: the wind is
!     zu and zv, which prepare_uvps has already gathered AND flipped, and the
!     air density and temperature are rhog and temp inside aerocore.
!
!     Flipped in latitude like every other field gathered here. See aero_surf
!     for why that is not optional and what it cost upstream.

      snow(:,:)  = 0.0
      wsoil(:,:) = 0.0
      if (ldustemit == 1) then
         call mpgagp(zgath,dsnow,1)
         do j=1,NLAT
            snow(:,NLAT+1-j) = zgath(:,j,1)
         end do
         call mpgagp(zgath,dwatc,1)
         do j=1,NLAT
            wsoil(:,NLAT+1-j) = zgath(:,j,1)
         end do
      end if

      if (mypid == NROOT .and. aero_debug) then
         write(nud,'(a,f11.2)') '* max aero u   =',maxval(abs(zu))
         write(nud,'(a,f11.2)') '* max v   =',maxval(abs(zv))
         write(nud,'(a,f11.2)') '* max ps0 =',maxval(zps0)
         write(nud,'(a,f11.2)') '* max ps1 =',maxval(zps1)
       ! write(nud,*)
       ! write(nud,*) 'ps0 NP'
       ! write(nud,*)  dp0(1:NLON)
       ! write(nud,*) 'zps0 NP'
       ! write(nud,*) zps0(1:NLON,NLAT)
       ! write(nud,*) 'ps0 SP'
       ! write(nud,*)  dp0(NLON*(NLAT-1)+1:)
       ! write(nud,*) 'zps0 SP'
       ! write(nud,*) zps0(1:NLON,1)
      end if
    
      if (mypid == NROOT) then

         call aerocore(daeros,numrhos,l_source,sigmah,dt,         &
                      zps0,zps1,zu,zv,                    &
                      dtoa,dtdx,apart,rhop,fcoeff,        &
                      aero_iord,aero_jord,aero_kord,      & 
                      NAERO,NLON,NLAT,NLEV,dap,dbk,       &
                      iml,ffsl_j1,ffsl_j2,js0,jn0,        &
                      colae,colad,rcolad,dlat,rcap,       &
                      aero_cnst,aero_deform,aero_zcross,  &
                      aero_fill,aero_mfct,aero_debug,nud, &
                      angle,land,aerosw,l_aerorad,prec, &
                      snow,wsoil)

!        preparation for the GUI output: 
!        invert the meridional direction and add the 360 deg. longitude

         do j=1,NLAT
            x(1:NLON,j,:,:) = daeros(:,NLAT+1-j,:,:)
            x(NLON+1,j,:,:) = daeros(1,NLAT+1-j,:,:)
         end do

!        send all tracer fields to output

         do jc=1,NAERO
            write(aero_name,'(a,i2.2)') 'DAEROS',jc
            call guiput(aero_name // char(0), x(1,1,1,jc), NLON+1,NLAT,NLEV)
         enddo

!        check the correlation between tracers 3 and 4

!        y(:,:,:) = x(:,:,:,3)+x(:,:,:,4)
!        call guiput('TRC03+04' // char(0), y, NLON+1,NLAT,NLEV)

      end if

      return
      end subroutine aero_main
