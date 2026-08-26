      module icemod

      use resmod
!
!     version identifier (date)
!
      character(len=80) :: version = '06.03.2013 by Larry'
!
!
!     Parameter
!
      parameter(NPRO = NPRO_ATM)        ! Number of processes (resmod)
      parameter(NLAT = NLAT_ATM)        ! Number of latitudes (resmod)
      parameter(NLON = NLAT + NLAT)     ! Number of longitudes
      parameter(NLPP = NLAT / NPRO)     ! Latitudes per process
      parameter(NHOR = NLON * NLPP)     ! Horizontal part
      parameter(NROOT = 0)              ! Master node
!
!
!     WHAT THE MODELLED SEA ICE AND ITS SNOW ARE MADE OF. All four below were
!     compile-time parameters, reachable only by a source edit, and each was an
!     Earth measurement standing where nothing said it had been chosen for this
!     world. They are icemod_nl keys now, at unchanged values, and
!     config/planet.yaml declares them under `surface.cryosphere` with a source
!     apiece. WORLD-04OK; notes/audits/cryosphere-material-properties.md has the
!     argument and analysis/ice_properties.py the arithmetic.
!
!     THE THREE SEA-ICE PROPERTIES ARE DECLARED, NOT DERIVED, and for the reason
!     CLFI below is already declared for: the density, the specific heat and the
!     conductivity of sea ice are all functions of its brine volume, hence of
!     the ice's own salinity and temperature, and THIS MODEL CARRIES NEITHER AS
!     A VARIABLE. It has one number per cell, a thickness. So each of the three
!     is a stated position about what this world's sea ice is, and IAPWS-06's
!     pure ice Ih is the bound it sits against rather than a substitute for it.
!
!     GRAVITY CANCELS OUT OF ALL THREE. It reaches a material property of ice
!     only through overburden pressure, and under a column far thicker than this
!     model carries the density response is parts per million at 12.81 m/s2 as
!     at 9.81. Where the surface gravity does reach this set is the SNOW's
!     density, which is set by compaction -- landmod's rhosnow, GRAV-8 -- and
!     from there the snow conductivity that now follows it.
!
      real :: CRHOI  = 920.    ! density of the modelled sea ice (kg/m**3)
      real :: CPI    = 2070.   ! specific heat of the modelled sea ice (J/(kg*K))
      real :: CKAPI  = 2.03    ! heat conductivity in the modelled sea ice (W/(m*K))
!
!     CLFSN IS THE ONE OF THE FOUR THAT IS A PURE SUBSTANCE. It is the melting
!     enthalpy of the modelled SNOW, and snow is ice Ih plus air with no brine
!     in it, so unlike the three above it is derivable: IAPWS-06 gives 333444.87
!     J/kg at the triple point against IAPWS-95's liquid water. The compiled
!     value below is 0.08 per cent above that and is kept, because run_exoplasim
!     writes the sourced number over it and the compiled value is only what
!     stands when nothing does.
      real :: CLFSN  = 3.337E5 ! melting enthalpy of the modelled snow (J/kg)
!
!     THE SNOW CONDUCTIVITY IS LANDMOD'S, NOT A SECOND COPY. WORLD-A9S5. This
!     was `parameter(CKAPSN = 0.31)`, the same number landmod declared as
!     `snowdiff` and used for the same purpose -- the conductive resistance of a
!     snow layer, `thickness / k` -- so the snow on the modelled sea ice and the
!     snow on the modelled soil were two statements of one material property
!     that nothing compared. landmod now DERIVES its value from the snow density
!     through Fourteau et al. (2021), and `iceini` takes it, so the two cannot
!     hold different snow and a bracket on the density moves both. The value
!     below is only what a run with no sea points never reads.
      real :: ckapsn = 0.3170  ! heat conductivity in the modelled snow (W/(m*K))
!
!     CPSN IS GONE. It was `parameter(CPSN = 2090.)`, the specific heat of snow,
!     and it was READ NOWHERE: a grep of the whole vendored tree finds the
!     declaration and no use. Deleting it is the fix, because exposing a
!     constant that reaches no arithmetic would put a knob in icemod_nl that
!     silently does nothing, and landmod's CPSNOW is the live declaration of the
!     same quantity. WORLD-04OK.
!
!     THE SNOW DENSITY IS LANDMOD'S, NOT A SECOND COPY. GRAV-8. 330 was declared
!     twice, here as a hardcoded parameter and in landmod as the namelist key
!     rhosnow, and both were used the same way: water equivalent times 1000 over
!     the density, to reach a physical thickness. The parameter here could not
!     be bracketed where it stood. `iceini` now takes landmod's value, so the
!     declared 330 below is only what a run with no sea points never reads.
      real :: crhosn = 330.             ! DENSITY OF SNOW (kg/m**3)
!
!     THE MELTING POINT IS PUMAMOD'S, NOT A SECOND COPY. This was
!     `parameter(TMELT=273.16)`, a compile-time constant, while `tmelt` is a
!     `planet_nl` key that `p_earth.f90` declares as "the freezing point every
!     soil, snow, sea and ice routine tests against". Both were live: seamod's
!     sea-ice albedo ramp anchors on pumamod's, icestep's skin-temperature melt
!     anchored on this one, and a configuration that moved the namelist key
!     would have split the two silently, with no line carrying the melting
!     point at all. Same class and same remedy as the snow density above:
!     `iceini` takes the one declaration, so the value below is only what
!     stands until it is called.
      real :: tmelt = 273.16            ! melting temp. for snow (0 deg C)
!
!     namelist parameters
!
      integer :: nice       = 1    ! compute ice yes/no (1/0)
      integer :: nseaice    = 1    ! Toggle sea ice (1/0)
      integer :: newsurf    = 0    ! 1: read surface data after restart
      integer :: nsnow      = 1    ! allow snow on ice yes/no (1/0)
      integer :: ntskin     = 1    ! compute skin temperature (0=clim.)
      integer :: ntspd      = 32   ! number of time steps per day
      integer :: noutput    = 1    ! master switch for output
      integer :: nout       = 32   ! output each nout timesteps
      integer :: nfluko     = 0    ! switch for flux correction
                                   ! (0 = none, 1 = heat-budget, 2 = newtonian)
      integer :: ncpl_ice_ocean = 1! coupling intervall ice - ocean

      integer :: nperpetual_ice = 0! perpetual climate conditions (day)
      integer :: nprint = 0        ! debug print out
      integer :: nprhor = 0        ! gp for debug printout
      integer :: ngui   = 0        ! switch for gui
      integer :: naout  = 0        ! no additional output fields 
!
!
!     THE OCEAN'S SALINITY REACHES THE MODEL THROUGH FOUR NUMBERS, NOT ONE.
!     All four are icemod_nl keys, and icemod passes the last three to
!     oceanini so that the two modules cannot hold different sea water.
!     Set them together: a bracket that moves only the freezing point moves
!     one of the four ways salinity acts.
!
!     TFREEZE  the freezing point, which sets where ice forms at all.
!     CRHOS    sea water density. It is the mixed-layer heat capacity with
!              CPS, and it is the snow-ice flooding threshold in subsnow as
!              the DIFFERENCE CRHOS-CRHOI, where a one per cent density
!              error is a ten per cent threshold error.
!     CPS      sea water specific heat. The compiled default is sea water's
!              at S=34.7 and its freezing point, from the UNESCO (1983)
!              polynomial. It was 4180, which is FRESH water at about 25 C.
!     CLFI     the heat of fusion of sea ice, depressed below pure ice's
!              3.337e5 by brine. DECLARED rather than derived: it is a
!              function of the ice's own salinity and temperature and this
!              model carries neither as a variable.
!
      real :: TFREEZE   =  271.25  ! freezing temp. for sea ice at S=34.7
      real :: CRHOS     = 1030.    ! density of sea water (kg/m**3)
      real :: CPS       = 3990.34  ! specific heat of sea water (J/(kg*K))
      real :: CLFI      = 3.28E5   ! heat of fusion of sea ice (J/kg)
      
!
!     THE DECLARED COLD START. This model has no sea surface temperature or
!     sea-ice climatology to begin from unless one is supplied as surface
!     codes 169, 210 and 211, and a world that has none is not a world whose
!     ocean can be guessed: an SST field is what this model PRODUCES. So a
!     cold start with no climatology begins from a stated profile instead of
!     a constructed one, and refuses to run if nothing states it.
!
!     The profile is hemispherically symmetric by construction --
!     tsst_pol + (tsst_eq-tsst_pol)*cos(lat)**2 -- because an initial
!     condition that is asymmetric between the hemispheres puts a difference
!     into the answer that nothing in the world put there. What it replaced
!     did exactly that: with code 169 absent xclsst held its -999 sentinel,
!     every ocean cell tested below the freezing point, and make_ice_thickness
!     shaped the resulting cover with an Earth Arctic-Antarctic table.
!
!     hice_ini is the thickness the cold start puts on cells whose declared
!     SST is at or below the freezing point. Zero is an ice-free cold start
!     and is the default, which is the side of the hysteresis a flux sweep
!     cannot otherwise approach from.
!
      real :: tsst_eq       = -999.! cold-start SST at the equator (K)
      real :: tsst_pol      = -999.! cold-start SST at the poles (K)
      real :: hice_ini      =  0.  ! cold-start sea-ice thickness (m)
!
      real :: taunc         =  0.  ! time scale for newtonian cooling
      real :: xmind         = 0.1  ! minimal ice thickness (m)
      real :: xmaxd         = 9.0  ! maximal ice thickness (m; neg. = no limit)
      real :: thicec        = 0.5  ! threshold to obtain make mask from comp. 
!
!     THE SECOND COMPACTNESS THRESHOLD, and it is a second one rather than a
!     restatement of thicec. thicec decides whether a cell is MASKED as iced,
!     at icestep and at mkicec; cicemin decides whether falling snow lands on
!     the ice or into the water, at subsnow's two sites. Both test xicec and
!     both stood at 0.5, but only thicec was reachable, so a bracket that moved
!     the mask left the snow partition on the compiled value. Exposed rather
!     than merged: the model uses them for different decisions and nothing in
!     the source says they are one number. Same remedy hlead got.
      real :: cicemin       = 0.5  ! minimum compactness to be ice
!
!     THE LEAD-CLOSING SCALE, and the whole of this model's lead
!     parameterisation. mkicec closes a cell's compactness with an e-folding
!     of hlead metres of new ice growth, and icestep then thresholds the
!     result at thicec into a hard mask, so hlead alone sets how much ice has
!     to grow before a cell counts as iced for albedo and roughness. It was a
!     local variable credited to Hippler 1979, unreachable from any namelist,
!     on a world whose year is half Earth's and which therefore grows a
!     different amount of ice per season.
!
      real :: hlead         = 0.5  ! lead-closing growth scale (m)
!
!     global integer
!
      integer :: nud        = 6    ! unit for messages
      integer :: nicec2d    = 0    ! 1: compute thickness from cover
      integer :: nstep      = 0    ! time step
      integer :: nrestart   = 0    ! restart switch
!
!     global real
!
      real :: xdt           = 0.   ! timestep (sec.)
      real :: solar_day  = 86400.0 ! length of day [sec]
!
!     global arrays
!
      real :: xshfl(NHOR)     = 0.   ! surface sensible heat flx
      real :: xshdt(NHOR)     = 0.   ! derivative of shfl w.r.t. temp
      real :: xlhfl(NHOR)     = 0.   ! surface latent heat flx
      real :: xlhdt(NHOR)     = 0.   ! derivative of slfl w.r.t. temp
      real :: xswfl(NHOR)     = 0.   ! surface short wave radiation
      real :: xlwfl(NHOR)     = 0.   ! surface long wave radiation
!
      real :: xts(NHOR)       = 0.   ! surface temperature (K)
      real :: xsst(NHOR)      = 0.   ! sea surface temperature (K)
      real :: xmld(NHOR)      = 0.   ! mixed layer depth (m) (from ocean)
      real :: xls(NHOR)       = 0.   ! land sea mask (1/0)
      real :: xiced(NHOR)     = 0.   ! ice thickness (m)
      real :: xicec(NHOR)     = 0.   ! ice cover (1/0)
      real :: xsnow(NHOR)     = 0.   ! snow depth (h2o equiv. m)
      real :: xsmelt(NHOR)    = 0.   ! snow melt  (m/s water equiv.)
      real :: xsmflx(NHOR)    = 0.   ! flux for snow melt  (w/m2)
      real :: ximelt(NHOR)    = 0.   ! flux for ice melt/freeze  (w/m2)
      real :: xsndch(NHOR)    = 0.   ! snow depth change (m/s water equiv.)
      real :: xqmelt(NHOR)    = 0.   ! res. heat flux for ice (W/m^2)
      real :: xstoi(NHOR)     = 0.   ! snow converted to ice (m h2o)
!
      real :: xicecc(NHOR)    = 0.  ! ice cover computed prognostically (frac.)
      real :: xaheat(NHOR)    = 0.  ! heat flux from atmosphere (W/m^2)
      real :: xheat(NHOR)     = 0.  ! heat flux from atm. (modified; W/m^2)
      real :: xcflux(NHOR)    = 0.  ! conductive heatflux through ice (W/m^2)
      real :: xcfluxf(NHOR)   = 0.  ! delta c-heatflux wrt tfreeze (W/m^2)
      real :: xcfluxr(NHOR)   = 0.  ! res. c-heatflux due to xmaxd (W/m^2)
      real :: xcfluxn(NHOR)   = 0.  ! res. c-heatflux due neg. ice (W/m^2)
      real :: xprs(NHOR)      = 0.  ! snow precipitation from atmosphere (m/s)
      real :: xpme(NHOR)      = 0.  ! fresh water flux (P-E only; m/s)
      real :: xroff(NHOR)     = 0.  ! runoff (m/s)
      real :: xtaux(NHOR)     = 0.  ! zonal wind stress (pa)
      real :: xtauy(NHOR)     = 0.  ! meridional wind stress (pa)
      real :: xust3(NHOR)     = 0.  ! ustar**3 (m**3/s**3)
      real :: xoheat(NHOR)    = 0.  ! heat flux input from ocean (w/m2)
      real :: xoflux(NHOR)    = 0.  ! heat flux from ocean (modified; w/m2)
      real :: xtsflux(NHOR)   = 0.  ! flux to warm/cool ice/snow (w/2)
      real :: xfluxc(NHOR)    = 0.  ! cond. heatflux (w/m2)
      real :: xscflx(NHOR)    = 0.  ! flux from snow -> ice conversion (w/m2)
      real :: xgw(NHOR)       = 0.  ! gaussian weights
      real :: xcoldsst(NHOR)  = 0.  ! declared cold-start SST profile (K)
!
!     Climatological fields
!
      real :: xclsst(NHOR,0:13) =-999.! climatological sst
      real :: xclicec(NHOR,0:13)=-999.! climatological ice cover
      real :: xcliced(NHOR,0:13)=-999.! climatological ice thickness
      real :: xflxice(NHOR,0:13)= 0.! flux correction (W/m^2)
      real :: xclsst2(NHOR)   = 0.  ! climatological sst
      real :: xclssto(NHOR)   = 0.  ! climatological sst (t-1)
      real :: xclicec2(NHOR)  = 0.  ! climatological ice cover
      real :: xcliced2(NHOR)  = 0.  ! climatological ice thickness
      real :: xflxice2(NHOR)  = 0.  ! flux correction (W/m^2)
!
!     fluxes for coupling to the ocean (accumulated)
!
      integer :: naccuo       = 0   ! counter for accumulation
      real :: cheat(NHOR)     = 0.  ! heat flux to the ocean (w/m2)
      real :: cpme(NHOR)      = 0.  ! fresh water flux (p-e; m/s)
      real :: croff(NHOR)     = 0.  ! runoff ( m/s)
      real :: ctaux(NHOR)     = 0.  ! zonal wind stress (pa)
      real :: ctauy(NHOR)     = 0.  ! meridional wind stress (pa)
      real :: cust3(NHOR)     = 0.  ! ustar**3 (m**3/s**)
      real :: csnow(NHOR)     = 0.  ! snow depth (m h2o eqv.)
!
!     accumulated diagnostics
!
      integer :: naccuout     = 0   ! counter for accumulation
      real :: xflxicea(NHOR)  = 0.  ! flux correction (w/m2)
      real :: xheata(NHOR)    = 0.  ! flux from the atmosphere (w/m2)
      real :: xofluxa(NHOR)   = 0.  ! flux from the ocean (w/m2)
      real :: xqmelta(NHOR)   = 0.  ! res flux to ice (w/m2)
      real :: xcfluxa(NHOR)   = 0.  ! flux to the ocean (w/m2)
      real :: xsmelta(NHOR)   = 0.  ! flux for snow melt (w/m2)
      real :: ximelta(NHOR)   = 0.  ! flux used for ice melt (w/m2)
      real :: xtsfluxa(NHOR)  = 0.  ! flux to warm/cool ice/snow (w/2)
      real :: xfluxca(NHOR)   = 0.  ! cond. heatflux (w/m2)
      real :: xcpmea(NHOR)    = 0.  ! fresh water (p-e; m/s) for lsg
      real :: xcroffa(NHOR)   = 0.  ! runoff  (m/s) for lsg
      real :: xstoia(NHOR)    = 0.  ! snow converted to ice (m h2o)
      real :: xscflxa(NHOR)   = 0.  ! flux from snow -> ice conversion (w/m2)
      real :: xcfluxra(NHOR)  = 0.  ! flux for limiting ice to xmaxd (w/m2)
      real :: xcfluxna(NHOR)  = 0.  ! flux due to neg. ice (w/m2)
!
!     entropy diagnostics
!
!
!     additional fields
!
      real,allocatable :: xaout(:,:)

      real :: deglat(NLPP) = 0.0    ! latitude in degrees
!
!     Parallel Stuff
!
      integer :: mpinfo  = 0
      integer :: mypid   = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO
!

!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(cheat,cicemin,ckapi,ckapsn,clfi,clfsn,cpi,cpme,cps,crhoi,crhos,crhosn,&
!$omp&  croff,csnow,ctaux,ctauy,cust3,deglat,&
!$omp&  mpinfo,mypid,myworld,&
!$omp&  naccuo,naccuout,naout,ncpl_ice_ocean,newsurf,nfluko,ngui,nice,nicec2d,nout,noutput,&
!$omp&  nperpetual_ice,nprhor,nprint,nproc,nrestart,nseaice,nsnow,nstep,ntskin,ntspd,nud,solar_day,&
!$omp&  taunc,tfreeze,thicec,tmelt,version,xaheat,xaout,xcflux,xcfluxa,xcfluxf,xcfluxn,xcfluxna,xcfluxr,&
!$omp&  xcfluxra,xclicec,xclicec2,xcliced,xcliced2,xclsst,xclsst2,xclssto,xcpmea,xcroffa,xdt,&
!$omp&  xfluxc,xfluxca,xflxice,xflxice2,xflxicea,xgw,xheat,xheata,xicec,xicecc,xiced,ximelt,ximelta,&
!$omp&  xlhdt,xlhfl,xls,xlwfl,xmaxd,xmind,xmld,xoflux,xofluxa,xoheat,xpme,xprs,xqmelt,xqmelta,xroff,&
!$omp&  xscflx,xscflxa,xshdt,xshfl,xsmelt,xsmelta,xsmflx,xsndch,xsnow,xsst,xstoi,xstoia,xswfl,xtaux,&
!$omp&  xtauy,xts,xtsflux,xtsfluxa,xust3,xcoldsst,tsst_eq,tsst_pol,hice_ini,hlead)

      end module icemod


!     ===========================
!     SUBROUTINE READ_ICE_SURFACE
!     ===========================

      subroutine read_ice_surface
      use icemod

!     THE -999 SENTINEL IS THE RECORD OF WHAT THIS WORLD HAS. mpsurfgp leaves
!     its argument untouched when the file is absent, so a field that stays at
!     -999 was not read, and every branch below that would prescribe the
!     surface from a climatology is refused on that test rather than fed a
!     construction. The sentinel is written to the restart with the field, so
!     the test answers the same on a continuation as on a cold start. Restore
!     it before each read: read_ice_surface is called again with newsurf = 1
!     after the restart arrays have been loaded.

      xclsst(:,:)  = -999.
      xclicec(:,:) = -999.
      xcliced(:,:) = -999.

      call mpsurfgp('xls',xls,NHOR,1)
      call mpsurfgp('xclsst' ,xclsst ,NHOR,14)
      if (nice .ge. 0) then
      call mpsurfgp('xclicec',xclicec,NHOR,14)
      call mpsurfgp('xcliced',xcliced,NHOR,14)
      endif
      
      if (nseaice == 0) then
         xclicec(:,:) = 0.
         xcliced(:,:) = 0.
      endif
      
!     make sure, that land sea mask values are 0 or 1

      where (xls(:) > 0.5)
         xls(:) = 1.0
      elsewhere
         xls(:) = 0.0
      endwhere

      call mpmaxval(xclicec,NHOR,14,zmax)
      if (nice .ge. 0) then
      if (zmax > 5.0) then
         xclicec(:,:) = xclicec(:,:) * 0.01
         if (mypid == NROOT) &
         write(nud,*) 'ice cover {xclicec} converted from % to fraction'
      endif

!     NO ICE COVER IS CONSTRUCTED FROM THE SST FIELD. What stood here tested
!     xclsst against the freezing point, and with xclsst at its sentinel that
!     put full cover on every ocean cell of every world that ships no code
!     210. The cover stays at the sentinel, iceget clamps it to zero for the
!     model, and the cold start builds a DECLARED initial state instead.

      if (zmax < 0.0 .and. mypid == NROOT) &
         write(nud,*) 'no ice cover {xclicec}: none read, and none constructed'

      call mpmaxval(xcliced,NHOR,14,zmaxd)
      if (zmax >= 0.0 .and. zmaxd < 0.0) then ! cover read, thickness not
         nicec2d = 1
         call make_ice_thickness
         if (mypid == NROOT) then
         write(nud,*) 'ice thickness {xcliced} computed from ice cover'
         write(nud,*) 'WARNING: by the CCM3 cover-to-thickness relation, which'
         write(nud,*) 'is fitted to Earth and is asymmetric between the'
         write(nud,*) 'hemispheres. Supply code 211 to avoid it.'
         endif
      endif
      endif

      if (nseaice == 0) then
         xclicec(:,:) = 0.0
         xcliced(:,:) = 0.0
      endif
!     correct climatological ice with land-sea mask
!     A sentinel is not a value to correct: zeroing it over land would leave
!     the array part -999 and part 0, and the test above reads a global max.

      do jm = 0 , 13
         where (xls(:) >= 1.0 .and. xclicec(:,jm) >= 0.0)
            xclicec(:,jm) = 0.0
         endwhere
         where (xls(:) >= 1.0 .and. xcliced(:,jm) >= 0.0)
            xcliced(:,jm) = 0.0
         endwhere
      enddo
      return
      end subroutine read_ice_surface


!     =================
!     SUBROUTINE ICEINI
!     =================

      subroutine iceini(kstep,krestart,koutput,kdpy,kgui,pts,psst,pmld  &
     &                 ,picec,piced,psnow,ktspd,psolday,pdeglat         &
     &                 ,prhosnow,psnowdiff,ptmelt                       &
                       ,icemod_namelist,oceanmod_namelist,ice_output    &
                       ,ocean_output)
      use icemod
      character (*) :: icemod_namelist
      character (*) :: ice_output
      character (*) :: oceanmod_namelist
      character (*) :: ocean_output
!
      real :: pts(NHOR)
      real :: psst(NHOR)
      real :: pmld(NHOR)
      real :: picec(NHOR)
      real :: piced(NHOR)
      real :: psnow(NHOR)
      real :: pdeglat(NLPP)
      real :: prhosnow
      real :: psnowdiff             ! landmod's snowdiff, derived from prhosnow
      real :: ptmelt                ! pumamod's melting point, a planet_nl key
      real (kind=8) :: zsi(NLAT)
      real (kind=8) :: zgw(NLAT)
      real :: zgw2(NLON,NLAT)

      logical :: lxsnow
!
      namelist/icemod_nl/nout,nfluko,nperpetual_ice,ntspd,nprint,nprhor &
     &               ,nice,nseaice,nsnow,ntskin,ncpl_ice_ocean,taunc   &
     &               ,xmind,xmaxd,thicec,cicemin,TFREEZE,CRHOS,CPS,CLFI  &
     &               ,CRHOI,CPI,CKAPI,CLFSN                              &
     &               ,tsst_eq,tsst_pol,hice_ini,hlead,newsurf,naout
!
!     copy input parameter to icemod
!
      nstep     = kstep
      nrestart  = krestart
      noutput   = koutput
      ngui      = kgui
      ntspd     = ktspd
      solar_day = psolday
      deglat(:) = pdeglat(:)
!     landmod's rhosnow, the one declaration of the snow density. GRAV-8.
      crhosn    = prhosnow
!     landmod's snowdiff, the one declaration of the snow conductivity, which
!     landini derived from that same density. WORLD-A9S5. Assigned BEFORE the
!     namelist is read, like the density above, because it is not an icemod_nl
!     key: a run that wants different snow moves rhosnow and gets both.
      ckapsn    = psnowdiff
!     pumamod's tmelt, the one declaration of the melting point. planet_nl sets
!     it and plasim.f90 broadcasts it before surfini, so every thread has the
!     configured value here.
      tmelt     = ptmelt

!     compute grids properties
!
      if(mypid == NROOT) then
       call inigau(NLAT,zsi,zgw)
       do jlat=1,NLAT
        zgw2(:,jlat)=zgw(jlat)
       enddo
      endif
      call mpscgp(zgw2,xgw,1)
!
!     get process id
!
      call ompi_info(nproc,mypid)
!
!     print version number and read namelist
!
      if (mypid == NROOT) then
         open(11,file=icemod_namelist)
         read(11,icemod_nl)
         close(11)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * ICEMOD ",a34," *")') trim(version)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist ICEMOD_NL from <",a15,"> *")') trim(icemod_namelist)
         write(nud,'(" *********************************************")')
         write(nud,icemod_nl)
         
         if (nseaice == 0) TFREEZE = 0.0
         
      endif

      call mpbci(nout)
      call mpbci(nfluko)
      call mpbci(nice)
      call mpbci(nseaice)
      call mpbci(newsurf)
      call mpbci(nsnow)
      call mpbci(ntskin)
      call mpbci(ntspd)
      call mpbci(nperpetual_ice)
      call mpbci(ncpl_ice_ocean)
      call mpbci(nprint)
      call mpbci(nprhor)
      call mpbci(naout)
      call mpbcr(taunc)
      call mpbcr(xmind)
      call mpbcr(xmaxd)
      call mpbcr(thicec)
      call mpbcr(cicemin)
      call mpbcr(TFREEZE)
      call mpbcr(CRHOS)
      call mpbcr(CPS)
      call mpbcr(CLFI)
!     The four material properties of the modelled ice and its snow, WORLD-04OK.
!     ckapsn is NOT here: it is landmod's, handed in above and already the same
!     on every thread because landini derived it on every thread.
      call mpbcr(CRHOI)
      call mpbcr(CPI)
      call mpbcr(CKAPI)
      call mpbcr(CLFSN)
      call mpbcr(tsst_eq)
      call mpbcr(tsst_pol)
      call mpbcr(hice_ini)
      call mpbcr(hlead)
!
!     set time step
!
      xdt   = solar_day / real(ntspd)
      taunc = solar_day * taunc

      if (nrestart == 0) then ! read start file
       
         call read_ice_surface
         call ice_cold_start

      else ! (nrestart /= 0)
!
!        restart from restart file
!
         if (mypid == NROOT) then
            call get_restart_integer('nstep',nstep)
            call get_restart_integer('naccuice',naccuout)
            call get_restart_integer('naccuo',naccuo)
            call get_restart_integer('nicec2d',nicec2d)
         endif
  
         call mpbci(nstep)
         call mpbci(naccuout)
         call mpbci(naccuo)
         call mpbci(nicec2d)
  
         call mpgetgp('xls'     ,xls     ,NHOR, 1)
         call mpgetgp('xts'     ,xts     ,NHOR, 1)
         call mpgetgp('xicec'   ,xicec   ,NHOR, 1)
         call mpgetgp('xiced'   ,xiced   ,NHOR, 1)
         call mpgetgp('xsnow'   ,xsnow   ,NHOR, 1)
         call mpgetgp('cheat'   ,cheat   ,NHOR, 1)
         call mpgetgp('cpme'    ,cpme    ,NHOR, 1)
         call mpgetgp('croff'   ,croff   ,NHOR, 1)
         call mpgetgp('ctaux'   ,ctaux   ,NHOR, 1)
         call mpgetgp('ctauy'   ,ctauy   ,NHOR, 1)
         call mpgetgp('cust3'   ,cust3   ,NHOR, 1)
         call mpgetgp('csnow'   ,csnow   ,NHOR, 1)
         call mpgetgp('xflxicea',xflxicea,NHOR, 1)
         call mpgetgp('xheata'  ,xheata  ,NHOR, 1)
         call mpgetgp('xofluxa' ,xofluxa ,NHOR, 1)
         call mpgetgp('xqmelta' ,xqmelta ,NHOR, 1)
         call mpgetgp('xcfluxa' ,xcfluxa ,NHOR, 1)
         call mpgetgp('xcfluxra',xcfluxra,NHOR, 1)
         call mpgetgp('xcfluxna',xcfluxna,NHOR, 1)
         call mpgetgp('xsmelta' ,xsmelta ,NHOR, 1)
         call mpgetgp('ximelta' ,ximelta ,NHOR, 1)
         call mpgetgp('xtsfluxa',xtsfluxa,NHOR, 1)
         call mpgetgp('xfluxca' ,xfluxca ,NHOR, 1)
         call mpgetgp('xscflxa' ,xscflxa ,NHOR, 1)
         call mpgetgp('xcpmea'  ,xcpmea  ,NHOR, 1)
         call mpgetgp('xcroffa' ,xcroffa, NHOR, 1)
         call mpgetgp('xstoia'  ,xstoia  ,NHOR, 1)
         call mpgetgp('xicecc'  ,xicecc  ,NHOR, 1)

         if (newsurf == 1) then 
            call read_ice_surface
         else
            call mpgetgp('xclicec' ,xclicec ,NHOR,14)
            call mpgetgp('xcliced' ,xcliced ,NHOR,14)
            call mpgetgp('xclsst'  ,xclsst  ,NHOR,14)
         endif ! (newsurf == 1)
        
!>> +++AYP   
         if (mypid==NROOT) inquire(file='restart_xsnow',exist=lxsnow)
         call mpbcl(lxsnow)
         if (lxsnow) call readarray(xsnow,'restart_xsnow')
!>> ---AYP         
                  
      endif ! (nrestart == 0)
!
!     WHAT NEEDS A CLIMATOLOGY, AND WHAT HAPPENS WHEN THERE IS NONE.
!
!     Three branches of this model do not integrate a surface, they prescribe
!     one from an observed climatology, and each is refused rather than run
!     against a construction. The test is the -999 sentinel in the field
!     itself, which survives the restart, so a continuation that switches one
!     of these on is refused the same way a cold start is.
!
!       nice == 0    prescribes ice cover and thickness from codes 210, 211
!       ntskin == 0  takes the skin temperature from the SST climatology, 169
!       nfluko /= 0  RELAXES the modelled state toward the climatology. That
!                    is the model's q-flux, and both of its branches read
!                    Earth artefacts here: nfluko = 1 reads code 709, and
!                    nfluko = 2 relaxes ice toward xcliced2 and compares SST
!                    against xclssto. Relaxing this world toward a field that
!                    was constructed from a sentinel is the whole of the
!                    hazard, and it is refused at the point of switching on.
!
      call mpmaxval(xclsst ,NHOR,14,zclsst)
      call mpmaxval(xcliced,NHOR,14,zcliced)
!
      if (nfluko /= 0 .and. zclsst < 0.) then
         call mpabort('icemod: nfluko needs a sea surface temperature '      &
     &              //'climatology (code 169) to relax toward, and this '    &
     &              //'run has none')
      endif
      if (nfluko == 2 .and. nseaice > 0 .and. zcliced < 0.) then
         call mpabort('icemod: nfluko = 2 relaxes sea ice toward a '         &
     &              //'thickness climatology (code 211), and this run has '  &
     &              //'none')
      endif
      if (nice == 0 .and. nseaice > 0 .and. zcliced < 0.) then
         call mpabort('icemod: nice = 0 prescribes sea ice from a '          &
     &              //'climatology (codes 210 and 211), and this run has '   &
     &              //'none')
      endif
      if (ntskin == 0 .and. zclsst < 0.) then
         call mpabort('icemod: ntskin = 0 takes the skin temperature from '  &
     &              //'a sea surface temperature climatology (code 169), '   &
     &              //'and this run has none')
      endif
!
!     read flux correction
!
      if (nfluko == 1) then
         xflxice(:,:) = -999.
         call mpsurfgp('xflxice',xflxice,NHOR,14)
         call mpmaxval(xflxice,NHOR,14,zflxice)
         if (zflxice < -900.) then
            call mpabort('icemod: nfluko = 1 needs the ice flux correction ' &
     &                 //'field (code 709), and this run has none')
         endif
      endif
!
!     open output file
!
      if (mypid == NROOT .and. noutput > 0) then
         open(71,file=ice_output,form='unformatted')
      endif ! (mypid == NROOT)
!
!     initialize ocean
!
      call oceanini(nstep,nrestart,noutput,kdpy,ngui,xsst,xmld,xoheat   &
     &             ,ntspd,solar_day,oceanmod_namelist,ocean_output      &
     &             ,TFREEZE,CRHOS,CPS,CLFI,CRHOI,xcoldsst)
!
      xoflux(:)=xoheat(:)
!
!     initialize skintemperature
!
      if (nrestart == 0) xts(:)=xsst(:)
!
!     copy output from icemod
!
      pts(:)=xts(:)
      psst(:)=xsst(:)
      pmld(:)=xmld(:)
      picec(:)=xicec(:)
      piced(:)=xiced(:)
      psnow(:)=xsnow(:)
!
!
      if(naout > 0) then
       allocate(xaout(NHOR,naout))
      endif
!
      return
      end subroutine iceini

!     =========================
!     SUBROUTINE ICE_COLD_START
!     =========================

      subroutine ice_cold_start
      use icemod
!
!     The initial ice and sea surface temperature of a run that starts from no
!     restart. Two cases, and the model must be told which it is in rather
!     than guessing: either a sea surface temperature climatology was read, in
!     which case the cold start is that climatology as it always was, or none
!     was, in which case it is the DECLARED profile in icemod_nl.
!
!     xicecc is set on BOTH paths. It is the prognostic compactness, mkicec
!     only grows it from whatever it holds, and icestep copies it into xicec
!     every step -- so leaving it at zero beside a non-zero thickness gave a
!     first output bin with ice thickness everywhere and ice cover identically
!     nowhere. Thick ice carrying no albedo is not a state this model can
!     start from.
!
      real :: zsst(NHOR)
      real :: zclsst
      real :: zpi
      integer :: jlat, jhor1, jhor2
!
      zpi = 4.*ATAN(1.)
!
      call mpmaxval(xclsst,NHOR,14,zclsst)
!
      if (zclsst > 0.) then
!
!        A climatology was read. Start from it.
!
         call iceget
         xiced(:)  = xcliced2(:)
         xicecc(:) = xclicec2(:)
!        Thickness without compactness is the same non-state from the other
!        direction, and it is what a run supplying code 211 and not 210 gets.
         where (xiced(:) > 0. .and. xicecc(:) <= 0.)
            xicecc(:) = 1.
         endwhere
      else
!
!        None was. Start from the declared profile.
!
         if (tsst_eq < 0. .or. tsst_pol < 0.) then
            call mpabort('icemod: a cold start with no sea surface '        &
     &                 //'temperature climatology (code 169) needs '        &
     &                 //'tsst_eq and tsst_pol declared in icemod_nl')
         endif
         do jlat = 1 , NLPP
            jhor1 = (jlat-1)*NLON + 1
            jhor2 = jlat*NLON
            zsst(jhor1:jhor2) = tsst_pol                                    &
     &          + (tsst_eq - tsst_pol)*COS(deglat(jlat)*zpi/180.)**2
         enddo
         xcoldsst(:) = zsst(:)
!
         xiced(:)  = 0.
         xicecc(:) = 0.
         where (xls(:) < 0.5 .and. zsst(:) <= TFREEZE)
            xiced(:)  = hice_ini
         endwhere
         where (xiced(:) > 0.)
            xicecc(:) = 1.
         endwhere
!
         if (mypid == NROOT) then
            write(nud,*) '* cold start from the declared SST profile:'
            write(nud,*) '*   equator ',tsst_eq,' K, pole ',tsst_pol,' K'
            write(nud,*) '*   sea ice ',hice_ini,' m below ',TFREEZE,' K'
         endif
      endif
!
!     the ice mask icestep applies to compactness, applied to the same field
!
      where (xicecc(:) >= thicec)
         xicec(:) = 1.
      elsewhere
         xicec(:) = 0.
      endwhere
!
      return
      end subroutine ice_cold_start

!     =====================================================================
!     SUBROUTINE icestep
!     =====================================================================

      subroutine icestep(pheat,pshfl,pshdt,plhfl,plhdt,plwfl,pswfl      &
     &                  ,ppme,proff,pprs,ptaux,ptauy,pust3              &
     &                  ,pts,picec,piced,psnow,psmelt,psndch,psst,pmld)
      use icemod
!
      real :: pheat(NHOR),pshfl(NHOR),pshdt(NHOR),plhfl(NHOR),plhdt(NHOR)
      real :: plwfl(NHOR),pswfl(NHOR),ppme(NHOR),proff(NHOR),pprs(NHOR)
      real :: ptaux(NHOR),ptauy(NHOR),pust3(NHOR),pts(NHOR),picec(NHOR)
      real :: piced(NHOR),psnow(NHOR),psmelt(NHOR),psndch(NHOR),psst(NHOR)
      real :: pmld(NHOR)
!
      real :: zsnowold(NHOR),zcflux(NHOR)
      real :: zicedold(NHOR)
!
!     debug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
      real,allocatable :: zprf5(:),zprf6(:),zprf7(:)
!
!     set some helpful bits
!
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     reset arrays
!
      zicedold(:) = 0.0
      ximelt(:)=0.
      xsmelt(:)=0.
      xsmflx(:)=0.
      xtsflux(:)=0.
      xfluxc(:)=0.
      xscflx(:)=0.
      xcfluxn(:)=0.
      xsndch(:)=0.
      xcfluxr(:)=0.
!
!     copy input to icemod
!
      xaheat(:)=pheat(:)
      xheat(:)=xaheat(:)
      xshfl(:)=pshfl(:)
      xshdt(:)=pshdt(:)
      xlhfl(:)=plhfl(:)
      xlhdt(:)=plhdt(:)
      xlwfl(:)=plwfl(:)
      xswfl(:)=pswfl(:)
      xpme(:)=ppme(:)
      xroff(:)=proff(:)
      xprs(:)=pprs(:)
      xtaux(:)=ptaux(:)
      xtauy(:)=ptauy(:)
      xust3(:)=pust3(:)
!
!     get climatology
!
      call iceget
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       call mpgagp(zprf1,xheat,1)
       call mpgagp(zprf7,xoflux,1)
       call mpgagp(zprf2,xclicec2,1)
       call mpgagp(zprf3,xcliced2,1)
       call mpgagp(zprf4,xsst,1)
       call mpgagp(zprf5,xicec,1)
       call mpgagp(zprf6,xiced,1)
       if(mypid==NROOT) then
        write(nud,*)'In icestep: nstep= ',nstep
        write(nud,*)'sst (old time step): ',zprf4(nprhor)
        write(nud,*)'heatflx from atm: ',zprf1(nprhor)
        write(nud,*)'heatflx from oce: ',zprf7(nprhor)
        write(nud,*)'clim. ice c and d: ',zprf2(nprhor),zprf3(nprhor)
        write(nud,*)'actual ice c and d: ',zprf5(nprhor),zprf6(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
      endif
!
!     zsnowold IS PRESET ON EVERY LANE. It is an automatic local with no
!     initialiser, was written only inside this mask, and is read below under the
!     same one -- `xsndch(:)=xsndch(:)+(xsnow(:)-zsnowold(:))`. A `where` masks
!     the ASSIGNMENT and not the evaluation, so on a land lane that difference
!     was formed against whatever the stack held, and the declared -ffpe-trap
!     turns an overflow or a signalling word there into SIGFPE. world-d016's
!     class. zicedold beside it already carries its own preset above.
!
!     THE PRESET IS INERT ON THE LANES THE MASK KEEPS. A sea lane assigns
!     zsnowold here, before anything reads it, and the only read is under the
!     identical mask. Copying xsnow makes the discarded lane's difference
!     exactly zero, which is what a lane with no ice column has changed.
      zsnowold(:)=xsnow(:)
      where(xls(:) < 0.5)
       zicedold(:)=xiced(:)
      endwhere
!
!     add new snow on ice
!
      if(nsnow==1) then
       where(xicec(:) >= cicemin) 
        xsnow(:)=xsnow(:)+xprs(:)*xdt
        xprs(:)=0.
       endwhere
      endif
!
!     compute surface temperature
!
      if(ntskin==1) then
       call skintemp
      else
       where(xiced(:)+1.E3/crhosn*xsnow(:) >= 0.1)
        xts(:)=xclsst2(:)
       elsewhere
        xts(:)=xsst(:)
       endwhere
      endif
!
!     make conductive heatflux
!
      call mkcflux
!
!     make new snow (alter conductive heatflux if neccessary)
!
      if(nsnow==1) then
       call subsnow
      else
       xsnow(:)=0.
      endif
!
!     make new ice
!
!     a) thickness
!
      call mkice
!
!     b) compactness (diagnostics for nice=0)
!
      call mkicec(zicedold,xiced,xicecc)
!
      if(nice == 0) then
!
!     climatological ice and 
!     flux correction diagnostics
!
       where(xiced(:) > 0. .or. xcliced2(:) > 0.)
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilfdt-xcflux(:)
        xcflux(:)=0.                    
        xiced(:)=xcliced2(:)
        xicec(:)=xclicec2(:)
       elsewhere
        xflxice2(:)=0.
        xicec(:)=0.
       endwhere
!
!     depug print out if needed
!
       if (nprint==2) then
        allocate(zprf1(NLON*NLAT))
        allocate(zprf2(NLON*NLAT))
        call mpgagp(zprf1,xflxice2,1)
        call mpgagp(zprf2,xiced,1)
        if(mypid==NROOT) then
         write(nud,*)'diagnose flux correction:'
         write(nud,*)'heat for correction: ',zprf1(nprhor)
         write(nud,*)'new iced: ',zprf2(nprhor)
        endif
        deallocate(zprf1)
        deallocate(zprf2)
       endif
       
       if (nseaice == 0) then
         xiced(:) = 0.
         xicec(:) = 0.
         xflxice2(:) = 0.
         xcflux(:) = 0.
       endif
       
      else if (nice < 0) then
        xiced(:) = 0.
        xcflux(:) = 0.
        xicec(:) = 0.
        xflxice2(:) = 0.
!
      else
!
!     a) set compactness
!
       xicec(:)=xicecc(:)         
!
!      b) flux correction (if switched on)
!
       if(nfluko == 0) then
        xflxice2(:)=0.
       elseif(nfluko == 1) then
        call getflx
        call addfci
        where(xiced(:) > 0. .and. xcliced2(:) > 0)
         xicec(:)=xclicec2(:)
        endwhere    
        where((xcliced2(:) > 0. .and. xicecc(:) > 0.999)                 &
     &    .or.(xcliced2(:)>0. .and. xmaxd>0. .and. xiced(:)>0.9*xmaxd))
         xicec(:)=1.
        endwhere
       elseif(nfluko == 2) then
        call mkflukoi
        call addfci
       endif
!
      endif
      
      if (nseaice == 0) then
        xicec(:) = 0.
        xiced(:) = 0.
        xflxice2(:) = 0.
        ximelt(:) = 0.
      endif
      
      
!
!     correct sea ice to a maximum of xmaxd and update ximelt
!
!     THE THICKNESS LIMIT PAYS FOR ITSELF, LOCALLY. Melting the excess takes
!     latent heat, and getiflx used to take it as a GLOBAL area-weighted sum
!     spread over every other cell with ice below xmaxd: a heat transport with
!     no physical carrier, and one that silently dropped the remainder
!     whenever the demand exceeded the capacity. The heat now comes out of the
!     same cell's conductive flux to the ocean, which is local, conserving,
!     and the only sink the cell has.
!
!     Set xmaxd negative to switch the limit off entirely, which is what a
!     world with no Earth Arctic to calibrate against should do: xmaxd is not
!     only a clamp, it also zeroes the conductive flux in mkcflux and skintemp
!     once ice reaches it, so a positive value stops basal growth outright.
!
      zcflux(:)=0.
      if(xmaxd >= 0.) then 
       where(xiced(:) > xmaxd)
        zcflux(:)=(xiced(:)-xmaxd)*zrhoilfdt
        xiced(:)=xmaxd
        xcfluxr(:)=xcfluxr(:)+zcflux(:)
        ximelt(:)=ximelt(:)+zcflux(:)
        xcflux(:)=xcflux(:)-zcflux(:)
!
!       diagnose the lost ice as accumulated snow 
!       (to make the budged from the atm. output) 
!
        xsndch(:)=xsndch(:)+zcflux(:)*1000./CRHOI/zrhoilfdt!/xdt
       end where
      endif
!
!     depug print out if needed
!
      if (nprint==2) then 
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,zcflux,1)
       call mpgagp(zprf2,xcfluxr,1)
       call mpgagp(zprf3,xcflux,1)
       call mpgagp(zprf4,xiced,1)
       if(mypid==NROOT) then 
        write(nud,*)'limit ice thickness (if > xmaxd):'
        write(nud,*)'heat for local correction: ',zprf1(nprhor)
        write(nud,*)'res. heat flux (after global adj.): ',zprf2(nprhor)
        write(nud,*)'new iced: ',zprf4(nprhor)
        write(nud,*)'new cflux: ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
!     make ice mask form compactness
!
      where(xicec(:) >= thicec)
       xicec(:)=1.
      elsewhere
       xicec(:)=0.
      end where
      
      if (nseaice == 0) then
        xicec(:) = 0.
        xiced(:) = 0.
        xflxice2(:) = 0.
        ximelt(:) = 0.
        xsnow(:) = 0.
      endif
      
!
!     correct snow with new ice (ice was melted below snow)
!     and melt snow (modify heat flux)
!
      zcflux(:)=0.
      where(xiced(:) <= 0. .and. xsnow(:) > 0.)
       zcflux(:)=xsnow(:)*1000.*CLFSN/xdt
       xsmelt(:)=xsmelt(:)+xsnow(:)/xdt
       xsmflx(:)=xsmflx(:)+zcflux(:)
       xcflux(:)=xcflux(:)-zcflux(:)
       xsnow(:)=0.
      end where
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,zcflux,1)
       call mpgagp(zprf2,xcflux,1)
       call mpgagp(zprf3,xsmflx,1)
       if(mypid==NROOT) then
        write(nud,*)'correct snow with new ice mask:'
        write(nud,*)'heat for correction: ',zprf1(nprhor)
        write(nud,*)'new flux for snow melt: ',zprf3(nprhor)
        write(nud,*)'new flux into ocean: ',zprf2(nprhor)
       endif
       zcflux(:)=xaheat(:)+xoheat(:)-xtsflux(:)-xsmflx(:)-ximelt(:)     &
     &          -xcflux(:)+xcfluxr(:)
       call mpgagp(zprf1,zcflux,1)  
       call mpgagp(zprf2,xflxice2,1) 
       if(mypid==NROOT) then
        write(nud,*)'final check for balance :'
        write(nud,*)'uncoupled: ',zprf1(nprhor)
        write(nud,*)'coupled without flux correction: ',zprf1(nprhor)
        zzz=zprf1(nprhor)+zprf2(nprhor)
        write(nud,*)'coupled with flux correction: ',zzz
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif

!
!     correct wind stress and ustar**3 according to ice
!
!!      where(xicec >= cicemin)
!!       xtaux(:)=0.
!!       xtauy(:)=0.
!!       xust3(:)=0.
!!      end where
!
!     possible input for other ocean models? (correct for ice/snow)
!
!     correct fresh water flux into ocean due to ice/snow changes
!
!     where(xls(:) < 0.5) 
!      xpme(:)=xpme(:)                                                  &
!    &        -(xsnow(:)-zsnowold(:))/xdt                               &
!    &        -(xiced(:)-zicedold(:))/xdt*CRHOI/1.E3
!     endwhere
!
!     output
!
!     accumulate
!
      xflxicea(:)=xflxicea(:)+xflxice2(:)
      xheata(:)=xheata(:)+xaheat(:)
      xofluxa(:)=xofluxa(:)+xoheat(:)
      xqmelta(:)=xqmelta(:)+xqmelt(:)
      ximelta(:)=ximelta(:)+ximelt(:)
      xsmelta(:)=xsmelta(:)+xsmflx(:)
      xcfluxa(:)=xcfluxa(:)+xcflux(:)
      xcfluxra(:)=xcfluxra(:)+xcfluxr(:)
      xcfluxna(:)=xcfluxna(:)+xcfluxn(:)
      xtsfluxa(:)=xtsfluxa(:)+xtsflux(:)
      xfluxca(:)=xfluxca(:)+xfluxc(:)
      xscflxa(:)=xscflxa(:)+xscflx(:)
      xcpmea(:)=xcpmea(:)+xpme(:)
      xcroffa(:)=xcroffa(:)+xroff(:)
      xstoia(:)=xstoia(:)+xstoi(:)
      naccuout=naccuout+1
!
!     write out
!
      if(nout > 0) then
      if(mod(nstep,nout) == 0) then
       xflxicea(:)=xflxicea(:)/REAL(naccuout)
       xheata(:)=xheata(:)/REAL(naccuout)
       xofluxa(:)=xofluxa(:)/REAL(naccuout)
       xqmelta(:)=xqmelta(:)/REAL(naccuout)
       xcfluxa(:)=xcfluxa(:)/REAL(naccuout)
       xcfluxra(:)=xcfluxra(:)/REAL(naccuout)
       xcfluxna(:)=xcfluxna(:)/REAL(naccuout)
       ximelta(:)=ximelta(:)/REAL(naccuout)
       xsmelta(:)=xsmelta(:)/REAL(naccuout)
       xtsfluxa(:)=xtsfluxa(:)/REAL(naccuout)
       xfluxca(:)=xfluxca(:)/REAL(naccuout)
       xscflxa(:)=xscflxa(:)/REAL(naccuout)
       xcpmea(:)=xcpmea(:)/REAL(naccuout)
       xcroffa(:)=xcroffa(:)/REAL(naccuout)
       xstoia(:)=xstoia(:)/REAL(naccuout)
       if(noutput > 0) call iceout
       xflxicea(:)=0.
       xheata(:)=0.
       xofluxa(:)=0.
       xqmelta(:)=0.
       xcfluxa(:)=0.
       xcfluxra(:)=0.
       xcfluxna(:)=0.
       ximelta(:)=0.
       xsmelta(:)=0.
       xtsfluxa(:)=0.
       xfluxca(:)=0.
       xscflxa(:)=0.
       xcpmea(:)=0.
       xcroffa(:)=0.
       xstoia(:)=0.
       naccuout=0
      endif
      endif
!
!     ocean coupling
!
!     accumulate fluxes
!
      where(xls(:) < 0.5)
       cheat(:)=cheat(:)+xcflux(:)
       cpme(:)=cpme(:)+xpme(:)
       croff(:)=croff(:)+xroff(:)
       ctaux(:)=ctaux(:)+xtaux(:)
       ctauy(:)=ctauy(:)+xtauy(:)
       cust3(:)=cust3(:)+xust3(:)
       csnow(:)=csnow(:)+xsnow(:)
      end where
      naccuo=naccuo+1
      if(mod(nstep,ncpl_ice_ocean) == 0) then
       where(xls(:) < 0.5)
        cheat(:)=cheat(:)/real(naccuo)
        cpme(:)=cpme(:)/real(naccuo)
        croff(:)=croff(:)/real(naccuo)
        ctaux(:)=ctaux(:)/real(naccuo)
        ctauy(:)=ctauy(:)/real(naccuo)
        cust3(:)=cust3(:)/real(naccuo)
        csnow(:)=csnow(:)/real(naccuo)
       endwhere

       call oceanstep(xicec,xiced,cheat,cpme,croff,ctaux,ctauy,cust3    &
     &               ,csnow,xsst,xmld,xoheat,xcliced2)
       xoflux(:)=xoheat(:)
!
       where(xls(:) < 0.5)
        cheat(:)=0.
        cpme(:)=0.
        croff(:)=0.
        ctaux(:)=0.
        ctauy(:)=0.
        cust3(:)=0.
        csnow(:)=0.
       end where
       naccuo=0
!
      endif
!
!     snow diagnostics
!
      where(xls(:) < 0.5) xsndch(:)=xsndch(:)+(xsnow(:)-zsnowold(:))!/xdt
!
!     copy output from icemod
!
      pts(:)=xts(:)
      picec(:)=xicec(:)
      piced(:)=xiced(:)
      psnow(:)=xsnow(:)
      psmelt(:)=xsmelt(:)
      psndch(:)=xsndch(:)
      psst(:)=xsst(:)
      pmld(:)=xmld(:)
!
!     advance time step
!
      nstep=nstep+1
!
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,xoflux,1)
       call mpgagp(zprf2,xicec,1)
       call mpgagp(zprf3,xiced,1)
       call mpgagp(zprf4,xsst,1)
       if(mypid==NROOT) then
        write(nud,*)'final ice c and d to atm and oce: ',zprf2(nprhor),zprf3(nprhor)
        write(nud,*)'heatflx and sst from ocean: ',zprf1(nprhor),zprf4(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
      return
      end subroutine icestep

!     ===================================
!     SUBROUTINE ICESTOP
!     ===================================

      subroutine icestop
      use icemod
!
!     close output file
!
      if (mypid == NROOT) then
       close(71)
      endif
!
!     write restart file
!
      if (mypid == NROOT) then
         call put_restart_integer('naccuice',naccuout)
         call put_restart_integer('naccuo',naccuo)
         call put_restart_integer('nicec2d',nicec2d)
      endif

      call mpputgp('xls'     ,xls     ,NHOR, 1)
      call mpputgp('xts'     ,xts     ,NHOR, 1)
      call mpputgp('xicec'   ,xicec   ,NHOR, 1)
      call mpputgp('xiced'   ,xiced   ,NHOR, 1)
      call mpputgp('xsnow'   ,xsnow   ,NHOR, 1)
      call mpputgp('xclicec' ,xclicec ,NHOR,14)
      call mpputgp('xcliced' ,xcliced ,NHOR,14)
      call mpputgp('xclsst'  ,xclsst  ,NHOR,14)
      call mpputgp('cheat'   ,cheat   ,NHOR, 1)
      call mpputgp('cpme'    ,cpme    ,NHOR, 1)
      call mpputgp('croff'   ,croff   ,NHOR, 1)
      call mpputgp('ctaux'   ,ctaux   ,NHOR, 1)
      call mpputgp('ctauy'   ,ctauy   ,NHOR, 1)
      call mpputgp('cust3'   ,cust3   ,NHOR, 1)
      call mpputgp('csnow'   ,csnow   ,NHOR, 1)
      call mpputgp('xflxicea',xflxicea,NHOR, 1)
      call mpputgp('xheata'  ,xheata  ,NHOR, 1)
      call mpputgp('xofluxa' ,xofluxa ,NHOR, 1)
      call mpputgp('xqmelta' ,xqmelta ,NHOR, 1)
      call mpputgp('xcfluxa' ,xcfluxa ,NHOR, 1)
      call mpputgp('xcfluxra',xcfluxra,NHOR, 1)
      call mpputgp('xcfluxna',xcfluxna,NHOR, 1)
      call mpputgp('xsmelta' ,xsmelta ,NHOR, 1)
      call mpputgp('ximelta' ,ximelta ,NHOR, 1)
      call mpputgp('xtsfluxa',xtsfluxa,NHOR, 1)
      call mpputgp('xfluxca' ,xfluxca ,NHOR, 1)
      call mpputgp('xscflxa' ,xscflxa ,NHOR, 1)
      call mpputgp('xcpmea'  ,xcpmea  ,NHOR, 1)
      call mpputgp('xcroffa' ,xcroffa ,NHOR, 1)
      call mpputgp('xstoia'  ,xstoia  ,NHOR, 1)
      call mpputgp('xicecc'  ,xicecc  ,NHOR, 1)
!
!     finalize ocean
!
      call oceanstop
!
!
      if(naout > 0) then
       deallocate(xaout)
      endif
!
      call finishup(xsnow,'newxsnow')
!       
      return
      end subroutine icestop

!     =====================================================================
!     SUBROUTINE mkice
!     =====================================================================

      subroutine mkice
      use icemod
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:)
!
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,xoflux,1)
       if(mypid==NROOT) then
        write(nud,*)'In mkice:'
        write(nud,*)'conductive hf and oce fl: ',zprf1(nprhor),zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
!     add flux from ocean
!
      where(xls(:) < 0.5) 
       xcflux(:)=xcflux(:)+xoflux(:)
      endwhere
!
!     make new ice
!
      where(xsst(:) <= TFREEZE .or. xiced(:) > 0.)               
!
!      add flux due to snow conversion
!
       ximelt(:)=xcflux(:)+xscflx(:)-xcfluxf(:)
!
!      new ice thickness 
!
       xiced(:)=xiced(:)-ximelt(:)/zrhoilfdt
!
!     reset heat flux 
!     (if ice > 0 and SST > TFREEZE (eg different climatologies)
!     use part of the conductive heat flux to cool the ocean)
!
       xcflux(:)=xcfluxf(:)
!
      end where
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,ximelt,1)
       call mpgagp(zprf3,xiced,1)
       if(mypid==NROOT) then
        write(nud,*)'new ice and flux for ice: ',zprf3(nprhor),zprf2(nprhor)
        write(nud,*)'new (residual) conductive hf: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
!     correct negative ice thickness (warm ocean)
!
      where(xiced(:) <= 0.)
       xcflux(:)=xcflux(:)-xiced(:)*zrhoilfdt
       ximelt(:)=ximelt(:)+xiced(:)*zrhoilfdt
       xcfluxn(:)=xcfluxn(:)-xiced(:)*zrhoilfdt
       xiced(:)=0.
      end where
!
!     set infinitisimal sea ice to zero
!
      where(ABS(xiced(:)) < 1.E-9) xiced(:)=0.
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,ximelt,1)
       call mpgagp(zprf3,xiced,1)
       if(mypid==NROOT) then
        write(nud,*)'final ice thickness: ',zprf3(nprhor)
        write(nud,*)'final flux for ice: ',zprf2(nprhor)
        write(nud,*)'final conductive hf: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      return
      end subroutine mkice

!     =====================================================================
!     SUBROUTINE mkicec
!     =====================================================================

      subroutine mkicec(picedo,picedn,picec)
      use icemod
      real :: picedo(NHOR)  ! old thickness (input)
      real :: picedn(NHOR)  ! new thickness (input)
      real :: picec(NHOR)   ! old and new compactness (input & output)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
!
!     compute new compactness (following Hippler '79; diagnostics)
!
      where(picedn(:) > picedo(:))
       picec(:)=picec(:)+(1.-picec(:))*(picedn(:)-picedo(:))/hlead
       picec(:)=AMIN1(picec(:),1.)
      endwhere
!     THE DIVISOR IS FLOORED BECAUSE THE MASK DOES NOT PROTECT IT. A `where`
!     selects which lanes the ASSIGNMENT stores and leaves the compiler free to
!     evaluate the right-hand side on every lane; the declared -ffpe-trap=zero
!     turns a lane that was going to be discarded into SIGFPE. world-d016, the
!     same mechanism as world-bhs and world-5a0.
!     Here the divisor IS the mask, negated: a lane the mask keeps has
!     picedo(:) > picedn(:) >= 0 and so picedo(:) > 0, and a lane it discards
!     has picedo(:) >= picedn(:), zero included. The floor is inert on every
!     kept lane -- 1.0e-30 m of ice is thirty decades below any thickness this
!     model carries -- and keeps the discarded lane's quotient finite, which a
!     floor at tiny() would not: -ffpe-trap=overflow is declared beside the
!     zero trap.
      where(picedn(:) < picedo(:))
       picec(:)=picec(:)+picec(:)*(picedn(:)-picedo(:))                 &
     &         /(2.*max(picedo(:),1.0e-30))
       picec(:)=AMAX1(picec(:),0.)
      endwhere
      where(picedn(:) <= 0.)
       picec(:)=0.
      endwhere
!
!     debug printout
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,picedo,1)
       call mpgagp(zprf2,picedn,1)
       call mpgagp(zprf3,picec,1)
       if(mypid==NROOT) then
        write(nud,*)'In mkicec: '
        write(nud,*)'old ice thickness: ',zprf1(nprhor)
        write(nud,*)'new ice thickness: ',zprf2(nprhor)
        write(nud,*)'new ice compactness: ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
      end subroutine mkicec

!     =====================================================================
!     SUBROUTINE subsnow
!     =====================================================================

      subroutine subsnow
      use icemod
!
      real :: zhice(NHOR) = 0. ! new ice thickness
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zhice)
      real :: zdice(NHOR) = 0. ! ice thickness change due to snow conversion
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zdice)
      real :: zdsnow(NHOR)= 0. ! snow depth change due to snow conversion
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zdsnow)
      real :: zqmelt(NHOR)= 0. ! residual qmelt going into ice
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zqmelt)
!
!     dbug arrays
!
      real,allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
!
!     a) melt snow falling into water (cool water, energy balance)
!
      where(xicec(:) < cicemin .and. xprs(:) > 0.)
       xcflux(:)=xcflux(:)-xprs(:)*1000.*CLFSN
       xsmelt(:)=xprs(:)+xsmelt(:)
       xsmflx(:)=xsmflx(:)+xprs(:)*1000.*CLFSN
       xprs(:)=0.
      endwhere
!
!     b) snow melt on ice
!
      zqmelt(:)=xqmelt(:)
      do jhor=1,NHOR
       if(xiced(jhor) > 0.) then
        if(xqmelt(jhor) > xsnow(jhor)*1000.*CLFSN/xdt) then
         xsmelt(jhor)=xsmelt(jhor)+xsnow(jhor)/xdt
         xsmflx(jhor)=xsmflx(jhor)+xsnow(jhor)*1000.*CLFSN/xdt
         zqmelt(jhor)=xqmelt(jhor)-xsnow(jhor)*1000.*CLFSN/xdt
         xsnow(jhor)=0.
        elseif(xqmelt(jhor) > 0. .and. xsnow(jhor) > 0.) then
         xsmelt(jhor)=xsmelt(jhor)+xqmelt(jhor)/CLFSN/1000.
         xsnow(jhor)=xsnow(jhor)-xqmelt(jhor)*xdt/CLFSN/1000.
         xsmflx(jhor)=xsmflx(jhor)+xqmelt(jhor)
         zqmelt(jhor)=0.
        endif
!
!     add residual flux to conductive heat flux 
!
        xcflux(jhor)=xcflux(jhor)+zqmelt(jhor)
!
       else
!
!       snow falls into water
!
        xsmelt(jhor)=xsmelt(jhor)+xsnow(jhor)/xdt
        xcflux(jhor)=xcflux(jhor)-xsnow(jhor)*1000.*CLFSN/xdt
        xsmflx(jhor)=xsmflx(jhor)+xsnow(jhor)*1000.*CLFSN/xdt
        xsnow(jhor)=0.
       endif
      enddo
!
!     dbug print out
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,xsmelt,1)
       call mpgagp(zprf2,xsnow,1)
       call mpgagp(zprf3,xsmflx,1)
       call mpgagp(zprf4,xcflux,1)
       if(mypid==NROOT) then
        write(nud,*)'in subsnow:'
        write(nud,*)'snow melt, and flx used for melting: ',zprf1(nprhor)    &
     &        ,zprf3(nprhor)
        write(nud,*)'new snow (m_snow): ',zprf2(nprhor)*1.E3/crhosn
        write(nud,*)'new conductive heat flux: ',zprf4(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
!     c) convert snow to sea ice if snow/ice interface is below sea level:
!     note: xsnow = h2o equiv.)
!     note: if CLFSN\=CLFI, this conversion does not conserve mass and
!     energy. Here, both is conserved by an additional flux 
!     (i.e. by using zdsnow*1.E3*CLFSN instead of zdice*CRHOI*CLFI)
!
      xstoi(:)=0.
      zdice(:)=0.
      zdsnow(:)=0.
      where(xls(:) < 0.5 .and. 1.E3*xsnow(:) > (CRHOS-CRHOI)*xiced(:))
       zhice(:)=(1.E3*xsnow(:)+CRHOI*xiced(:))/CRHOS
       zdice(:)=zhice(:)-xiced(:)
       zdsnow(:)=-zdice(:)*CRHOI/1.E3
       xsnow(:)=xsnow(:)+zdsnow(:)
       xscflx(:)=zdsnow*1.E3*CLFSN/xdt
!
!     diagnose ice thickness change from snow conversion
!     (e.g. to correct P-E in LSG-coupling)
!
       xstoi(:)=zdice(:)*CRHOI/1000./xdt
!
!     diagnose the converted snow as snowmelt 
!
       xsmelt(:)=xsmelt(:)-zdsnow/xdt
       xsmflx(:)=xsmflx(:)-zdsnow*1.E3*CLFSN/xdt
      endwhere
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,zdice,1)
       call mpgagp(zprf2,zdsnow,1)
       call mpgagp(zprf3,xsnow,1)
       call mpgagp(zprf4,xscflx,1)
       if(mypid==NROOT) then
        zzfli=zprf1(nprhor)*CRHOI*CLFI/xdt
        zzfls=zprf2(nprhor)*1000.*CLFSN/xdt
        zdfl=zzfli+zzfls
        write(nud,*)'snow change by snow -> ice conv. (m/s,W): '        &
     &              ,zprf2(nprhor),zzfls
        write(nud,*)'ice change by snow -> ice conv. (m/s,W): '         &
     &              ,zprf1(nprhor),zzfli
        write(nud,*)'new snow (m_snow): ',zprf3(nprhor)*1.E3/crhosn
        write(nud,*)'flux to build ice (total, residual): '             &
     &              ,zprf4(nprhor),zdfl
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
      return
      end subroutine subsnow

!     =====================================================================
!     SUBROUTINE MKCFLUX
!     =====================================================================

      subroutine mkcflux
      use icemod
!
      real :: zckap(NHOR) = 0.
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zckap)
      real :: zhsnow(NHOR)= 0.
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zhsnow)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
!
!     preset
!
      xfluxc(:)=0.
!
!     compute conductive heat flux
!
      where(xiced(:) < xmind) 
       xcflux(:)=xheat(:)
       xcfluxf(:)=0.
      end where
!     THE DIVISOR IS FLOORED BECAUSE THE MASK DOES NOT PROTECT IT. A `where`
!     selects which lanes the ASSIGNMENT stores and leaves the compiler free to
!     evaluate the right-hand side on every lane; the declared -ffpe-trap=zero
!     turns a lane that was going to be discarded into SIGFPE. world-d016, the
!     same mechanism as world-bhs and world-5a0.
!     A lane the mask keeps has xiced(:) >= xmind, so the divisor is at least
!     xmind/CKAPI; a lane it discards can have xiced(:) and xsnow(:) both zero,
!     which is open water and is the common case, and the divisor is then
!     EXACTLY zero. zhsnow is threadprivate with a zero initialiser and was
!     written only inside this mask, so a discarded lane carried a stale value
!     as well. It is computed on every lane now, which costs one multiply and
!     makes the value a discarded lane holds a defined one.
!
!     The floor is inert on every kept lane by twenty-eight decades and keeps
!     the discarded lane's quotient finite: the numerator there is zhsnow+xiced,
!     which is the same near-zero quantity, so the quotient is of order one.
      zhsnow(:)=1.E3/crhosn*xsnow(:)
      where(xiced(:) >= xmind)
       zckap(:)=(zhsnow(:)+xiced(:))                                    &
     &         /max(zhsnow(:)/CKAPSN+xiced(:)/CKAPI,1.0e-30)
      endwhere
!
!     limit ice thickness 
!
      if(xmaxd >= 0.) then
       where(xiced(:) >= xmaxd)
        zckap(:)=0.
       endwhere
      endif
!
!     FLOORED FOR THE REASON GIVEN AT zckap ABOVE: a lane the mask discards can
!     be open water with zhsnow(:) and xiced(:) both zero, and the divisor is
!     then exactly zero. A kept lane has xiced(:) >= xmind, so the floor is
!     inert there, and on a discarded lane zckap(:) is zero, so the quotient is
!     zero rather than an overflow. world-d016.
      where(xiced(:) >= xmind) 
       xfluxc(:)=zckap(:)*(xts(:)-xsst(:))                              &
     &          /max(zhsnow(:)+xiced(:),1.0e-30)
       xcflux(:)=xfluxc(:)
       xcfluxf(:)=zckap(:)*(TFREEZE-xsst(:))                            &
     &           /max(zhsnow(:)+xiced(:),1.0e-30)
      end where
!
      if(nfluko > 0. .or. nice == 0) xcfluxf(:)=0.
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,zhsnow,1)
       call mpgagp(zprf2,zckap,1)
       call mpgagp(zprf3,xcflux,1)
       call mpgagp(zprf4,xcfluxf,1)
       if(mypid==NROOT) then
        write(nud,*)'in mkcflux: '
        write(nud,*)'snow depth (in m snow): ',zprf1(nprhor)
        write(nud,*)'kappa and conductive hf: ',zprf2(nprhor),zprf3(nprhor)
        write(nud,*)'flux to adjust water temp. below ice: ',zprf4(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
      return
      end subroutine mkcflux

!     =====================================================================
!     SUBROUTINE iceout
!     =====================================================================

      subroutine iceout
      use icemod
!
      integer :: ih(8)
!
      real :: zsnow(NHOR) = 0.
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zsnow)
!
      real,allocatable :: zprf1(:),zprf2(:)
!
      where (xls(:) < 0.5) zsnow(:) = 1000./crhosn *xsnow(:)
!
      call ntomin(nstep,nmin,nhour,nday,nmonth,nyear)
!
      ih(2) = 0
      ih(3) = nyear*10000 +nmonth*100 +nday
      ih(4) = nhour*100 +nmin
      ih(5) = NLON
      ih(6) = NLAT
      ih(7) = 0
      ih(8) = 0
!
      ih(1) = 701
      call mpwritegph(71,xheata,NHOR,1,ih)
      ih(1) = 702
      call mpwritegph(71,xofluxa,NHOR,1,ih)
      ih(1) = 703
      call mpwritegph(71,xtsfluxa,NHOR,1,ih)
      ih(1) = 704
      call mpwritegph(71,xsmelta,NHOR,1,ih)
      ih(1) = 705
      call mpwritegph(71,ximelta,NHOR,1,ih)
      ih(1) = 706
      call mpwritegph(71,xcfluxa,NHOR,1,ih)
      ih(1) = 707
      call mpwritegph(71,xfluxca,NHOR,1,ih)
      ih(1) = 708
      call mpwritegph(71,xqmelta,NHOR,1,ih)
      ih(1) = 709
      call mpwritegph(71,xflxicea,NHOR,1,ih)
      ih(1) = 710
      call mpwritegph(71,xicec,NHOR,1,ih)
      ih(1) = 711
      call mpwritegph(71,xiced,NHOR,1,ih)
      ih(1) = 712
      call mpwritegph(71,xscflxa,NHOR,1,ih)
      ih(1) = 713
      call mpwritegph(71,xcfluxra,NHOR,1,ih)
      ih(1) = 714
      call mpwritegph(71,xcfluxna,NHOR,1,ih)
      ih(1) = 739
      call mpwritegph(71,xts,NHOR,1,ih)
      ih(1) = 741
      call mpwritegph(71,zsnow,NHOR,1,ih)
      ih(1) = 769
      call mpwritegph(71,xsst,NHOR,1,ih)
      ih(1) = 772
      call mpwritegph(71,xls,NHOR,1,ih)
      ih(1) = 790
      call mpwritegph(71,xclicec2,NHOR,1,ih)
      ih(1) = 791
      call mpwritegph(71,xcliced2,NHOR,1,ih)
      ih(1) = 792
      call mpwritegph(71,xicecc,NHOR,1,ih)    
      ih(1) = 794
      call mpwritegph(71,xcpmea,NHOR,1,ih)
      ih(1) = 795
      call mpwritegph(71,xcroffa,NHOR,1,ih)
      ih(1) = 796
      call mpwritegph(71,xstoia,NHOR,1,ih)
      if(naout > 0) then
       do ja=1,naout
        ih(1)=750+ja
        call mpwritegph(71,xaout(1,ja),NHOR,1,ih)
       enddo
      endif
!
!     diagnostics
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,xicec,1)
       call mpgagp(zprf2,xiced,1)
       if(mypid==NROOT) then
        write(nud,*)'In iceout:'
        write(nud,*)'ice compactness: ',zprf1(nprhor)
        write(nud,*)'ice thickness:   ',zprf2(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
      return
      end subroutine iceout

!====================================================================
!     SUBROUTINE MKFLUKOI
!====================================================================

      subroutine mkflukoi
      use icemod
!
      zrhoilf=CRHOI*CLFI
!
      xflxice2(:)=0.
      if(taunc > 0.) then
       where(xls(:) < 0.5)
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilf/taunc
       end where
      else
       where(xls(:) < 0.5)
        xflxice2(:)=(xiced(:)-xcliced2(:))*zrhoilf/xdt                  
       end where
      endif
!
      return
      end subroutine mkflukoi

!     =====================================================================
!     SUBROUTINE addfci
!     =====================================================================

      subroutine addfci
      use icemod
!
      real :: zmelt(NHOR) = 0.
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zmelt)
      real :: zflr(NHOR) = 0.
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zflr)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:)
!
      zrhoilfdt=CRHOI*CLFI/xdt
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       call mpgagp(zprf1,xflxice2,1)
       call mpgagp(zprf2,xclsst2,1)
       call mpgagp(zprf3,xcliced2,1)
       if(mypid==NROOT) then
        write(nud,*)'In addfci:'
        write(nud,*)'fluko: ',zprf1(nprhor)
        write(nud,*)'clsst: ',zprf2(nprhor)
        write(nud,*)'cliced: ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
      endif
!
!     add flux correction (ensure energy conservation)
!
      zflr(:)=0.
      do jhor=1,NHOR
       if(xls(jhor) < 1.) then
!
!     if ice or (clim. ice and t <= tclim):
!     make new ice (from fluko) and substract fluko from conductive hfl
!
        if((xsst(jhor) <= TFREEZE) .or. (xiced(jhor) > 0.)              &
     &    .or. (xcliced2(jhor) > 0. .and. xsst(jhor) <= xclssto(jhor))) &
     &  then
         zmelt(jhor)=xcflux(jhor)+xflxice2(jhor)
         xiced(jhor)=xiced(jhor)-zmelt(jhor)/zrhoilfdt
         ximelt(jhor)=ximelt(jhor)+zmelt(jhor)
         xcflux(jhor)=0.
        else
!
!     else: give flux correction to the ocean  
!
         zflr(jhor)=xflxice2(jhor)
        endif
       endif
      enddo
!
      where(xls(:) < 1.)
       xcflux(:)=xcflux(:)+zflr(:)
      endwhere
!
!!    zsum(1)=SUM(zflr(:)*xgw(:),MASK=(xls(:) < 1.))
!!    zsum(2)=SUM(xgw(:),MASK=(xls(:) < 1.))
!!    call mpsumbcr(zsum,2)
!!    if(zsum(1) /= 0.) then
!!     where(xls(:) < 1.)
!!      xcflux(:)=xcflux+zsum(1)/zsum(2)
!!     end where
!!    endif
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       call mpgagp(zprf1,zmelt,1)
       call mpgagp(zprf2,xiced,1)
       call mpgagp(zprf3,xcflux,1)
       call mpgagp(zprf4,zflr,1)
       if(mypid==NROOT) then
        write(nud,*)'new ice and flux for ice: ',zprf2(nprhor),zprf1(nprhor)
        write(nud,*)'residual flux correction: ',zprf4(nprhor)
        write(nud,*)'new flux into ocean: ',zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
      endif
!
!     correct negative ice thickness (warm ocean)
!
      where(xiced(:) <= 0.)
       xcflux(:)=xcflux(:)-xiced(:)*zrhoilfdt
       ximelt(:)=ximelt(:)+xiced(:)*zrhoilfdt
       xcfluxn(:)=xcfluxn(:)-xiced(:)*zrhoilfdt
       xiced(:)=0.
      end where
!
!     set infinitisimal sea ice to zero
!
      where(ABS(xiced(:)) < 1.E-9) xiced(:)=0.
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       call mpgagp(zprf1,xcflux,1)
       call mpgagp(zprf2,xiced,1)
       if(mypid==NROOT) then
        write(nud,*)'final ice: ',zprf2(nprhor)
        write(nud,*)'final conductive hf: ',zprf1(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
      endif
!
      return
      end subroutine addfci

!====================================================================
!     SUBROUTINE SKINTEMP
!====================================================================

      subroutine skintemp
      use icemod
!
      parameter(stb=5.67E-8) ! stefan boltzmann constant
!
      real :: ztso(NHOR)         ! old ts for diagnostics
      real :: zhsnow(NHOR)       ! snow depth (m snow)
      real :: zckap_mean(NHOR)   ! kappa (ice+snow layer)
      real :: zcflux(NHOR)       ! conductive heat flux (diagnostics)
!
!     debug arrays
!
      real, allocatable :: zprf1(:),zprf2(:),zprf3(:),zprf4(:),zprf5(:)
      real, allocatable :: zprf6(:),zprf7(:),zprf8(:),zprf9(:),zprf10(:)
      real, allocatable :: zprf11(:)
!
!     F.LUNKEIT   UNIHH     FEB-07
!
!     PURPOSE.
!     --------
!     CALCULATE ICE SKIN-TEMPERATURE AS A PROGNOSTIC VARIABLE
!
!     
      zrcpl=xmind*CRHOI*CPI
!
!     preset fluxes
!
      zcflux(:)=xheat(:)
      xtsflux(:)=0.
      xqmelt(:)=0.
!
!     save old ts of ice (max of tmelt)
!
      ztso(:)=AMIN1(TMELT,xts(:))
!
!     start calculations
!
      do jhor=1,NHOR
       zhsnow(jhor)=1.E3/crhosn*xsnow(jhor)
       if(xiced(jhor) >= xmind) then
        zcpdt=zrcpl/xdt
        zckap_mean(jhor)=(zhsnow(jhor)+xiced(jhor))                     &
     &                  /(zhsnow(jhor)/CKAPSN+xiced(jhor)/CKAPI)
!
!     new skin temperature (implicit w.r.t. heat conduction)
!
        zkapz=zckap_mean(jhor)/(xiced(jhor)+zhsnow(jhor))
!
!     set ice flux to 0 if ice >= xmaxd
!
        if(xiced(jhor) >= xmaxd .and. xmaxd >= 0.) zkapz=0.
!
        zflx=xheat(jhor)+zkapz*xsst(jhor)
        xts(jhor)=(zcpdt*ztso(jhor)+zflx)/(zcpdt+zkapz)
!
!     residual flux going into ice/snow melt
!
        if(xts(jhor) > TMELT) then
         xqmelt(jhor)=xheat(jhor)+zkapz*(xsst(jhor)-TMELT)              &
     &               -zcpdt*(TMELT-ztso(jhor))
         xts(jhor)=TMELT
        endif
!
!     diagnose fluxes
!
        zcflux(jhor)=zkapz*(xts(jhor)-xsst(jhor))
        xtsflux(jhor)=zcpdt*(xts(jhor)-ztso(jhor))
!
       else
        xts(jhor)=xsst(jhor)
        if(xiced(jhor) > 0.) then
         xqmelt(jhor)=xheat(jhor)
        endif
       endif
      enddo
!
!     depug print out if needed
!
      if (nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       allocate(zprf8(NLON*NLAT))
       allocate(zprf9(NLON*NLAT))
       allocate(zprf10(NLON*NLAT))
       allocate(zprf11(NLON*NLAT))
       call mpgagp(zprf1,zhsnow,1)
       call mpgagp(zprf2,zckap_mean,1)
       call mpgagp(zprf5,ztso,1)
       call mpgagp(zprf6,xts,1)
       call mpgagp(zprf7,xheat,1)
       call mpgagp(zprf8,xtsflux,1)
       call mpgagp(zprf9,xqmelt,1)
       call mpgagp(zprf10,zcflux,1)
       call mpgagp(zprf11,xoflux,1)
       if(mypid==NROOT) then
        write(nud,*)'in skintemp: '
        write(nud,*)'modified atm. heat flux: ',zprf7(nprhor)
        write(nud,*)'modified oce. heat flux: ',zprf11(nprhor)
        write(nud,*)'snow depth (in m snow): ',zprf1(nprhor)
        write(nud,*)'kappa: ',zprf2(nprhor)
        write(nud,*)'old and new Ts: ',zprf5(nprhor),zprf6(nprhor)
        write(nud,*)'global min Ts: ',MINVAL(zprf6)
        write(nud,*)'conductive heat flux (w. new ts): ',zprf10(nprhor)
        write(nud,*)'heat flux used to warm/cool ice: ',zprf8(nprhor)
        write(nud,*)'res. flux to melt snow/ice (T>TM): ',zprf9(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
       deallocate(zprf8)
       deallocate(zprf9)
       deallocate(zprf10)
       deallocate(zprf11)
      endif
!
!     entropy diagnostics
!
!
      return
      end subroutine skintemp


!     =================
!     SUBROUTINE ICEGET
!     =================

      subroutine iceget
      use icemod
!
!     get sea ice climatology for the actual time step 
!     (using linear interpolation from monthly means)
!
      call momint(nperpetual_ice,nstep,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      xclssto(:) = zgw1 * xclsst(:,jm1) + zgw2 * xclsst(:,jm2) ! SST (t-1)
      call momint(nperpetual_ice,nstep+1,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      xclsst2(:) = zgw1 * xclsst(:,jm1) + zgw2 * xclsst(:,jm2) ! SST
      xcliced2(:) = zgw1 * xcliced(:,jm1) + zgw2 * xcliced(:,jm2)
      xclicec2(:) = zgw1 * xclicec(:,jm1) + zgw2 * xclicec(:,jm2)
      xclicec2(:)=AMAX1(xclicec2(:),0.)
      xcliced2(:)=AMAX1(xcliced2(:),0.)
!
!     no compactness for hice=0
!
      where (xcliced2(:) <= 0.)
       xclicec2(:) = 0.
      endwhere
!
!     debug output
!
      if (nprint == 2) then
       if (mypid == NROOT) then
        write(nud,*)'In iceget:'
        write(nud,*)' jm1,jm2= ',jm1,jm2,' gw1,gw2= ',zgw1,zgw2
       endif
      endif  
!
      return
      end subroutine iceget

!     =================
!     SUBROUTINE GETFLX
!     =================

      subroutine getflx
      use icemod

      call momint(nperpetual_ice,nstep+1,jm1,jm2,zgw)
      xflxice2(:) = (1.0 - zgw) * xflxice(:,jm1) + zgw * xflxice(:,jm2)

      return
      end subroutine getflx

!     =============================
!     SUBROUTINE MAKE_ICE_THICKNESS
!     =============================

      subroutine make_ice_thickness
      use icemod

      parameter(cminn=0.1 ,cmins=0.25)
      parameter(cmaxn=0.9 ,cmaxs=1.0 )
      parameter(hminn=0.25,hmins=0.25)
      parameter(hmaxn=3.0 ,hmaxs=0.50)

      real :: zc(NLON,NLAT,0:13)
      real :: zd(NLON,NLAT,0:13)

      real, parameter :: zhfac(0:13)=(/0.912,0.942,1.,1.058,1.124,1.161,1.175 &
                        ,1.058,0.931,0.883,0.88,0.876,0.912,0.942/)

!     convert ice compactness to thickness (see CCM3 report pp 127-129)

      call mpgagp(zc,xclicec,14)
      zd(:,:,:) = 0.0
      
      if (mypid == NROOT) then
         do jm = 0 , 13
   
      !     northern hemisphere
      
            do jlat = 1 , NLAT/2
               where (zc(:,jlat,jm) >= cmaxn)
                  zd(:,jlat,jm)=hmaxn*zhfac(jm)
               elsewhere (zc(:,jlat,jm) >= cminn)
                  zd(:,jlat,jm)=(hminn+(hmaxn-hminn)*(zc(:,jlat,jm)-cminn) &
                               /(cmaxn-cminn))*zhfac(jm)
               elsewhere (zc(:,jlat,jm) > 0.0)
                  zd(:,jlat,jm)=hminn*zhfac(jm)
               endwhere
            enddo
      
      !     southern hemisphere
      
            do jlat = NLAT/2+1 , NLAT
               where (zc(:,jlat,jm) >= cmins)
                  zd(:,jlat,jm)=hmaxs
               elsewhere (zc(:,jlat,jm) > 0.0)
                  zd(:,jlat,jm)=hmins+zc(:,jlat,jm)
               endwhere
            enddo
   
         enddo ! jm
      endif ! (mypid == NROOT)

      call mpscgp(zd,xcliced,14)

      return
      end subroutine make_ice_thickness

     
