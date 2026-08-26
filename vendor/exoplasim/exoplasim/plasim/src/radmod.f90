      module radmod
!
!     radiation module for PUMA
!
!**   1) make global PUMA variables available
!
      use pumamod
!
!**   2) define global parameters for *subs* included in radmod
!
!*    2.0) version identifier (date)
!
      character(len=80) :: rversion = '23.09.2019 by Adiv'
!
!*    2.1)  constant parameters
!

      parameter(SBK = 5.67E-8)  ! Stefan-Bolzman Const.
!
!     Number of aerosol SPECIES the radiation can carry at once. CLIM-39.
!
!     The scheme used to carry exactly one, because ssa, the backscatter ratio
!     and the band-2 extinction ratio were scalars: one set of optical
!     properties for the whole planet. That is not a resolution question but a
!     species question, so this is a plain parameter and not a resolution one.
!     Four covers what aeolian/ produces and can price -- mineral dust, sea
!     salt, volcanic sulfate -- with one slot spare for the carbonaceous
!     aerosol that CLIM-29 is waiting on a biosphere run for. Raising it costs
!     memory in aodsp and nothing else.
!
!     NOT pumamod's NAERO, and the two must not be conflated. NAERO is how many
!     TRANSPORTED tracers aerocore carries and sizes daeros and numrhos; NAERSP
!     is how many RADIATIVE species radmod mixes, and most of those are
!     prescribed columns that aerocore never sees. One transported tracer can
!     occupy one radiative species, which is what makes the prescribed and
!     interactive paths coexist, but the counts are independent.
!
      parameter(NAERSP = 4)      ! max # of aerosol species carried at once
!       parameter(zsolar1=0.517)  
!       parameter(zsolar2=0.483)  

!
!*    2.2) namelist parameters (see *sub* radini)
!

      real    :: starbbtemp = 5772.0 ! Star's blackbody surface temperature (K)
      logical :: lstarfile = .false.
      integer :: nstarfile = 0      ! integer version of the logical
      integer :: l_aerorad = 0 ! Aerosol scattering and absorption off (0) or on (1)
      character(len=80) :: starfile = " " !Name of input stellar spectrum file
      character(len=80) :: starfilehr = " " !Name of hi-res version of input spectrum
      character(len=128) :: aerofile = " " ! Name/path to file constaining aerosol optical data
      
      real    :: gsol0   = 1367.0 ! solar constant (set in planet module)
      integer :: nsolcycle = 0    ! switch for sinusoidal stellar-flux cycle
      integer :: gsolstart = 0    ! absolute model step at cycle phase zero
      real    :: gsolamp = 0.0    ! cycle semi-amplitude (W/m2)
      real    :: gsolperiod = 1.0 ! cycle period (model timesteps)
      real    :: gsolphase = 0.0  ! phase offset (cycles; zero=mean, rising)
!     Second, independent component. A real activity cycle is not one
!     sinusoid -- epsilon Eridani carries a short and a long period at once --
!     and this world superposes a medium climatic cycle on a long geomorphic
!     one. The two periods are deliberately non-commensurate, so the deepest
!     minima drift instead of repeating on a fixed beat. Amplitude zero
!     disables it and a one-component run behaves exactly as before.
      real    :: gsolamp2 = 0.0    ! second semi-amplitude (W/m2)
      real    :: gsolperiod2 = 1.0 ! second period (model timesteps)
      real    :: gsolphase2 = 0.0  ! second phase offset (cycles)
      real    :: solclat = 1.0    ! cos of lat of insolation if ncstsol=1
      real    :: solcdec = 1.0    ! cos of dec of insolation if ncstsol=1
      real    :: clgray  = -1.0   ! cloud grayness (-1 = computed)
      real    :: th2oc   = 0.024  ! absorption coefficient h2o continuum (lwr)
      real    :: tswr1   = 0.077  ! tuning of cloud albedo range1
      real    :: tswr2   = 0.065  ! tuning of cloud back scattering c. range2
      real    :: tswr3   = 0.0055 ! tuning of cloud s. scattering alb. range2
      real    :: tpofmt  = 1.00   ! tuning of point of mean transmittance
      real    :: acllwr  = 0.100  ! mass absorption coefficient for clouds (lwr)
!
!     SURFACE LONGWAVE EMISSIVITY, split by the land-sea mask. Upstream wrote
!     zeps = dls + 0.98*(1-dls) as a literal in lwr, so the modelled land was a
!     perfect blackbody by construction and everything the mask does not call
!     land -- ocean and sea ice -- was 0.98, with no comment, no unit and no
!     source anywhere. Named here so a run records what it emitted with; the
!     defaults reproduce that literal exactly and no run changes.
!
!     Both are broadband thermal emissivities, dimensionless, 0 to 1. The
!     compiled 1.0 over land is a blackbody and is not a measurement; what the
!     configuration passes is derived, per rock class, from measured
!     directional-hemispherical reflectance spectra area-weighted over this
!     world's land, and sits near 0.94 with a preparation bracket of 0.92 to
!     0.95. Sea water is 0.985 to 0.99, so 0.98 is slightly low.
!
!     ONE SCALAR PER SURFACE IS ENOUGH, MEASURED AND NOT ASSUMED. A per-cell
!     land field out of the same lithology map buys, over a scalar set at that
!     field's own land mean, under half a W/m2 of surface net longwave at the
!     resolutions this project runs, against a criterion of 1.4 -- the top of
!     the model's own dry adiabatic energy sink, which the surface fluxes pay
!     for. So there is no NHOR array here, deliberately, and the reason is
!     recorded rather than left to be re-derived.
      real    :: elwland = 1.0    ! surface lw emissivity, land
      real    :: elwsea  = 0.98   ! surface lw emissivity, ocean and sea ice
      real    :: a0o3    = 0.25   ! parameter to define o3 profile
      real    :: a1o3    = 0.11   ! parameter to define o3 profile
      real    :: aco3    = 0.08   ! parameter to define o3 profile
      real    :: bo3     = 20000. ! parameter to define o3 profile
      real    :: co3     = 5000.  ! parameter to define o3 profile
      real    :: toffo3  = 0.25    ! parameter to define o3 profile
      real    :: o3scale = 1.0    ! scale o3 concentration
!     Spectral re-weighting of the Lacis & Hansen (1974) ozone absorptances for a
!     non-solar host. Their three terms give absorptance as a fraction of TOTAL
!     incident SOLAR flux, so each carries the Sun's share of flux in the band it
!     represents. Under a differently-shaped spectrum those shares change, and
!     the unmodified code applies solar band weights to a non-solar star. Both
!     default to 1.0, reproducing Lacis & Hansen exactly, so a solar-host run is
!     bit-identical.
      real    :: o3uvw   = 1.0    ! weight, Hartley-Huggins UV terms
      real    :: o3visw  = 1.0    ! weight, Chappuis visible term
!     The same correction in the larger term. Lacis & Hansen's water vapour
!     absorptance is their Eq. 21, a fit to Yamamoto (1962), and Yamamoto states
!     the definition outright: the ratio to the SOLAR CONSTANT of the energy
!     absorbed by the whole air column. It is a fraction of total incident flux
!     and the Sun's spectrum is inside it, because Yamamoto built it by weighting
!     laboratory band absorptivities with the solar flux and summing them.
!
!     Dividing by zsolar2 below converts it to a fraction of band-2 flux and then
!     multiplies it back by band-2 flux, so the two cancel and the ABSORBED FLUX
!     is Lacis & Hansen's solar value for any star at all. A redder host puts
!     more of its flux in the near infrared where water vapour absorbs, and gets
!     the Sun's absorption anyway. h2osww is that band re-weighting.
!
!     Default 1.0 reproduces Lacis & Hansen exactly, so a solar-host run is
!     bit-identical.
      real    :: h2osww  = 1.0    ! weight, near-infrared H2O bands
!
!     h2oswl is a SECOND and independent correction to the same term, and the
!     two must not be folded into one key. h2osww re-weights Eq. 21 for a
!     non-solar host and is a star-over-Sun RATIO, so any error in Eq. 21's
!     absolute level divides straight out of it. h2oswl is that level.
!
!     Eq. 21 is a fit to Yamamoto (1962), and a modern line list absorbs 12 to
!     13% more at the same absorber amount before any continuum: correlated-k
!     over Eq. 21 runs 1.10 to 1.15 across the range, median 1.134, and 1.127 at
!     this planet's operating path. The cause is a pressure treatment, Yamamoto's
!     bands being measured at a Curtis-Godson effective pressure while Eq. 21 is
!     refitted as though they held at standard pressure. Neither side carries the
!     MT_CKD-convention continuum and the continuum absorbs in the WINDOWS, so it
!     adds to the line-by-line side and 1.127 is a floor: with it the level is
!     1.163 over a bracket of 1.129 to 1.206, and the bracket ends are arms to
!     run rather than an error bar.
!
!     Default 1.0 reproduces Lacis & Hansen exactly, so a solar-host run is
!     bit-identical and a rebuilt binary reproduces every run that exists until
!     the namelist turns this on. exoplasim/notes/corrk-cross-check.md is the
!     measurement; TASKS.md PHYS-9 is the decision to carry it as a key rather
!     than as a line in the error budget.
      real    :: h2oswl  = 1.0    ! level, near-infrared H2O absorptance
!     Shortwave CO2, which this scheme does not have at all. swr carries ozone in
!     band 1 and water vapour in band 2 and nothing else; CO2 appears only in
!     lwr, from Sasamori (1968). Lacis & Hansen did not parameterise the
!     near-infrared CO2 bands either, so the port is faithful and the ABSORBER is
!     simply missing -- and a redder host puts about 1.5x the Sun's share of its
!     flux into those bands.
!
!     The absorptance below is built the way Yamamoto built the water vapour one:
!     Howard, Burch & Williams (1956) band absorptions weighted by a SOLAR
!     spectrum and charged only with what water vapour leaves them. So it is a
!     fraction of TOTAL incident flux and is divided by zsolar2 for the same
!     reason Eq. 21 is, and co2sww is the band re-weighting for a non-solar host,
!     exactly as h2osww is.
!
!     ZERO IS THE DEFAULT AND MEANS THE TERM IS ABSENT, not that it is
!     solar-weighted. Upstream has no shortwave CO2, so zero is what reproduces
!     it and a rebuilt binary is bit-identical until the namelist turns this on.
!     1.0 is the solar-weighted term; the value for a given host is the ratio of
!     its CO2 absorptance to the Sun's, which is what
!     exoplasim/scripts/shortwave_band_weights.py derives.
      real    :: co2sww  = 0.0    ! weight, near-infrared CO2 bands; 0 = absent

!     CH4 AND N2O IN THE LONGWAVE. CLIM-42.
!
!     Sasamori (1968) fits water vapour, CO2 and ozone and nothing else, so
!     these gases had no coefficient to set and no key to turn: adding either
!     means adding a band. The band model is Donner and Ramanathan (1980),
!     which is the paper that did exactly this job -- a band ABSORPTANCE, the
!     currency this scheme is already written in, rather than a
!     top-of-atmosphere forcing fit.
!
!     ZERO IS ABSENT and is the default, so a rebuilt binary reproduces the
!     old answer bit for bit until the namelist turns these on. That is the
!     reduction identity the change is tested against.
!
!     Both are scalars rather than fields because the photochemistry these
!     values come from puts both gases well mixed through this world's
!     troposphere -- CH4 falls 1.601 to 1.583 ppmv over the lowest 20 km and
!     N2O is flat to 15 km. A field would assert a structure that is not there.
!     exoplasim/notes/trace-gas-band.md has the argument and the sourcing.
      real    :: ch4     = 0.0    ! CH4 volume mixing ratio (ppmv); 0 = absent
      real    :: n2o     = 0.0    ! N2O volume mixing ratio (ppmv); 0 = absent
      integer :: no3     = 1      ! switch for ozon (0=no,1=yes,2=datafile)
      integer :: nsol    = 1      ! switch for solang (1/0=yes/no)
      integer :: nswr    = 1      ! switch for swr (1/0=yes/no)
      integer :: nlwr    = 1      ! switch for lwr (1/0=yes/no)
      integer :: necham  = 1      ! switch for using ECHAM-3 solar zenith angle 
                                  ! dependence for ocean albedo (1/0=yes/no)
      integer :: necham6  = 0     ! switch for using ECHAM-6 solar zenith angle 
                                  ! dependence for ocean albedo (overrides necham) (1/0=yes/no)
      integer :: nclouds = 1      ! switch for cloud sw effects (1/0=yes/no)
      integer :: nswrcl  = 1      ! switch for computed cloud props.(1/0=y/n)
      integer :: nrscat  = 1      ! switch for rayleigh scat. (1/0=yes/no)
      integer :: newrsc  = 0      ! switch for layer-by-layer rayleigh scat. (1/0=yes/no)
      integer :: nradice = 1      ! Whether to include sea ice reflectance (1/0=yes/no)
      integer :: ndcycle = 1      ! switch for daily cycle of insolation
                                  !  0 = daily mean insolation)
      integer :: ncstsol = 0      ! switch to set constant insolation
                                  ! on the whole planet (0/1)=(off/on)
      integer :: iyrbp   = -50    ! Year before present (1950 AD)
                                  ! default = 2000 AD
                                  
      integer :: npbroaden = 1    ! Should pressure broadening depend on surface pressure (1/0)
      integer :: nfixed  = 0      ! Switch for fixed zenith angle (0/1=no/yes)
      real    :: slowdown = 1.0   ! Factor by which to change diurnal insolation cycle
      real    :: desync = 0.0     ! Degrees per minute by which substellar point drifts (+/-)
      
      
      real    :: minwavel = 316.036116751 ! Minimum wavelength to use when computing spectra [nm]
      
      integer :: nstartemp = 0    ! Switch for using the star's bb temp to determine sw (0/1)
      integer :: nsimplealbedo = 1  ! Compute broadband albedo and use it for both bands
      
      real :: rcl1(3)=(/0.15,0.30,0.60/) ! cloud albedos spectral range 1
      real :: rcl2(3)=(/0.15,0.30,0.60/) ! cloud albedos spectral range 2
      real :: acl2(3)=(/0.05,0.10,0.20/) ! cloud absorptivities spectral range 2
      
!
!     PER-SPECIES optical constants. These were scalars, which is what limited
!     the scheme to one aerosol: dust ABSORBS and sea salt's single-scattering
!     albedo is 1 to within 1e-5, so no one value describes both, and the two
!     are co-located over the ocean wherever dust has been transported off the
!     land. aeolian/notes/multi-species-aerosol.md section 2.
!
      real :: ssa1(NAERSP) = 0. ! Single scattering albedo band 1
      real :: ssa2(NAERSP) = 0. ! Single scattering albedo band 2
      real :: qex1(NAERSP) = 0. ! Extinction efficiency band 1
      real :: qex2(NAERSP) = 0. ! Extinction efficiency band 2
      real :: bscat1(NAERSP) = 0. ! Backscattering ratio band 1
      real :: bscat2(NAERSP) = 0. ! Backscattering ratio band 2
!
!     The aerofile's own species axis, which already existed set to 1:
!     readdat(filename,ndim,nitems,kdata) fills kdata(nitems,ndim), so ndim IS
!     the species count and the file gains a COLUMN per species rather than
!     changing shape. Eight rows: Qext, Qsca, Qback, g for band 1, then the
!     same four for band 2.
!
      real :: aeroqs(8,NAERSP) = 0.  ! Array to read in aerosol optical constants
!
      integer :: naerosp = 0    ! number of active species; prescribed first,
                                ! then the transported one if it is on
      real :: aqlw(NAERSP) = 0.  ! thermal-IR absorption ratio per active
                                ! species, gathered from dustqlw and aeroqlw so
                                ! the longwave loop does not branch on path
      real :: apart = 50e-09 ! Aerosol particle radius. AEROMOD DECLARES ITS
                             ! OWN; aero_ini copies that one into this and
                             ! radini broadcasts it. The default is the
                             ! photochemical haze and is not this world's dust.
!
!*    2.2b) PRESCRIBED DUST (DUST-11)
!
!     A supplied, non-interacting dust field: the radiation sees a column
!     optical depth read as a surface boundary field and nothing transports,
!     emits or removes anything. That is enough to answer what dust does to
!     precipitation and runoff one iteration deep, and it deliberately does NOT
!     touch aerocore, so the bottom-level sink, the settling term and the
!     interactive number density are all out of the path.
!
!     ddustcol is the BAND 1 (0.34-0.75 um) column extinction optical depth.
!     Band 2 follows from the aerofile's own ratio of extinction efficiencies
!     and the thermal infrared from dustqlw, so one field carries all three and
!     the spectral ratios stay where the optics are.
!
!     ndustrad is now a COUNT of prescribed species rather than a switch, and
!     0 and 1 mean exactly what they meant before. Species s reads surface code
!     1810+s: 1811 is the dust field that already exists, 1812 onward are the
!     species added since. Each carries its own scale height, because a sea
!     salt layer sits in the boundary layer and dust does not, and a shared one
!     would put the second species at the first one's height while the column
!     total still looked right.
!
      integer :: ndustrad = 0     ! number of PRESCRIBED aerosol species (0 = off)
      real    :: dustsc(NAERSP)  = 1.0    ! multiplier on the prescribed column optical depth
      real    :: dusthsc(NAERSP) = 3000.0 ! aerosol scale height (m), concentration e-folding
      real    :: dustqlw(NAERSP) = 0.0    ! thermal-IR ABSORPTION optical depth per unit
                                  ! band-1 extinction optical depth. There is no
                                  ! defensible default: a shortwave-only dust is
                                  ! worse than no dust, so radini ABORTS if this
                                  ! is left at zero for any active species.
      logical :: ldustchk = .true.  ! report the column normalisation once
!
!*    2.2c) INTERACTIVE AEROSOL, AND ITS LONGWAVE (DUST-3 item 5)
!
!     The aerosol aerocore transports acted in the two SHORTWAVE bands and
!     nowhere else, so the model could cool with dust and could not warm with
!     it. On this world that is not a refinement: analysis/dust_forcing.json
!     prices shortwave-only dust against the full calculation at several W/m2
!     in the global mean, which is second in the whole error budget.
!
!     The physics is the prescribed path's, unchanged: a grey ABSORBER at the
!     same 1.66 diffusivity the cloud term uses, multiplied into the total
!     layer transmissivity so that overlap with water vapour and CO2 is handled
!     by construction. What is new here is only WHERE the optical depth comes
!     from -- the model's own number density instead of a boundary field.
!
!     aeroqlw is that ratio for the interactive aerosol and dustqlw is the one
!     for the prescribed field. They are kept apart because the two paths carry
!     different particles: the prescribed field is this world's dust at the
!     offline chain's size distribution, while the interactive tracer is
!     whatever aero_nl's apart and rhop describe. There is no defensible
!     default for either, so radini ABORTS on an interactive aerosol with
!     aeroqlw left at zero, exactly as it does for ndustrad with dustqlw.
!
!     The two paths were mutually exclusive and radini refused both at once,
!     because a prescribed column and a transported one are two aerosols and
!     they wrote into ONE slot, so adding their optical depths counted one of
!     them twice. With a species array they no longer share a slot: the
!     transported tracer is one species of it, at index ndustrad+1, and the
!     refusal is gone rather than preserved. CLIM-39.
!
!     That is what keeps DUST-13 open. While the exclusion stood, choosing
!     interactive dust put sea salt out of the radiation permanently, because
!     ndustrad and iaerint could not both be on.
!
      real    :: aeroqlw = 0.0    ! thermal-IR ABSORPTION optical depth per unit
                                  ! band-1 extinction optical depth, INTERACTIVE
                                  ! aerosol. No defensible default; radini
                                  ! aborts if it is zero with the aerosol on.
      integer :: iaerint = 0      ! 1 where the transported aerosol acts on the
                                  ! radiation, else 0. Set once in radini from
                                  ! l_aero and l_aerorad, so swr, lwr and
                                  ! radstep all ask the same question.
!
!*    2.3) arrays
!

      real :: gmu0(NHOR)                   ! cosine of solar zenit angle
      real :: gmu1(NHOR)                   ! cosine of solar zenit angle
      real :: ddustcol(NHOR,NAERSP)     = 0. ! prescribed band-1 column optical depth
      real :: ddustod(NHOR,NLEV,NAERSP) = 0. ! band-1 optical depth per layer
      real :: daerod(NHOR,NLEV)  = 0.      ! interactive band-1 optical depth per layer
!
!     The per-species band-1 extinction optical depth per layer, gathered from
!     the prescribed columns and the transported tracer once per radiation
!     step. swr and lwr both read THIS, for the reason aeroprof already gives
!     for daerod: a second copy of the arithmetic in either is a copy that can
!     drift, and the layer thicknesses depend on the temperature profile.
!
      real :: aodsp(NHOR,NLEV,NAERSP) = 0.  ! band-1 optical depth per layer per species
!       real :: dtdtlwr(NHOR,NLEV)           ! lwr temperature tendencies (now in pumamod)
!       real :: dtdtswr(NHOR,NLEV)           ! swr temperature tendencies (now in pumamod)

      real, allocatable :: dqo3cl(:,:,:)   ! climatological O3 (used if NO3=2)

      real :: zsolars(2) = 0.0             ! Container for storing solar constants
      
!
!*    2.4) scalars
!

      real :: gdist2 = 1.        ! Earth-sun distance factor ( i.e. (1/r)**2 )
      real :: time4rad = 0.      ! CPU time for radiation
      real :: time4swr = 0.      ! CPU time for short wave radiation
      real :: time4lwr = 0.      ! CPU time for long wave radiation
      
      real :: zsolar1 = 0.517    ! spectral partitioning 1 (wl < 0.75mue)
      real :: zsolar2 = 0.483    ! spectral partitioning 2 (wl > 0.75mue)
      real :: rcoeff = 1.0       ! Rayleigh scattering coefficient for cross section dependence
      
!
!     2.5 orbital parameters
!
      integer, parameter :: ORB_UNDEF_INT  = 2000000000  
      real :: obliqr   ! Earth's obliquity in radians
      real :: meananom0r = 0.0 ! Initial mean anomaly in radians
      real :: lambm0   ! Mean longitude of perihelion at the
                       ! vernal equinox (radians)
      real :: mvelpp   ! Earth's moving vernal equinox longitude
                       ! of perihelion plus pi (radians)
      real :: eccf=0.  ! Earth-sun distance factor ( i.e. (1/r)**2 )
      real :: orbnu=0. ! Earth true anomaly in radians.
      real :: lambm=0. ! Solar ecliptic longitude in radians
      real :: rasc=0.  ! Solar right ascension in radians
      real :: zcdayf=0. ! Fractional day
      real :: zdeclf=0. !Declination angle
      integer :: iyrad ! Year AD to calculate orbit for
      logical, parameter :: log_print = .true.
                       ! Flag to print-out status information or not.
                       ! (This turns off ALL status printing including)
                       ! (error messages.)
!
!     2.6 extended entropy/energy diagnostics
!
      real :: dftde1(NHOR,NLEP),dftde2(NHOR,NLEP)
      real :: dftue1(NHOR,NLEP),dftue2(NHOR,NLEP)
      real :: dftu0(NHOR,NLEP),dftd0(NHOR,NLEP)
!
!     auxiliary variables for solar zenit angle calculations
!
      real :: solclatcdec      ! cos(lat)*cos(decl) 
      real :: solslat          ! sin(lat)
      real :: solsdec          ! sin(decl)
      real :: solslatsdec      ! sin(lat)*sin(decl) 
      real :: zmuz             ! temporary zenit angle   
!

!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(a0o3,a1o3,acl2,acllwr,aco3,aerofile,aeroqlw,aeroqs,aodsp,apart,aqlw,bo3,bscat1,&
!$omp&  bscat2,ch4,clgray,co2sww,co3,daerod,ddustcol,ddustod,desync,dftd0,dftde1,dftde2,dftu0,&
!$omp&  dftue1,dftue2,dqo3cl,dusthsc,dustqlw,dustsc,eccf,elwland,elwsea,gdist2,gmu0,gmu1,gsol0,&
!$omp&  gsolamp,gsolamp2,&
!$omp&  gsolperiod,gsolperiod2,gsolphase,gsolphase2,gsolstart,h2oswl,h2osww,iaerint,iyrad,iyrbp,&
!$omp&  l_aerorad,lambm,lambm0,ldustchk,lstarfile,meananom0r,minwavel,mvelpp,n2o,naerosp,nclouds,&
!$omp&  ncstsol,ndcycle,ndustrad,necham,necham6,newrsc,nfixed,nlwr,no3,npbroaden,nradice,nrscat,&
!$omp&  nsimplealbedo,nsol,nsolcycle,nstarfile,nstartemp,nswr,nswrcl,o3scale,o3uvw,o3visw,obliqr,&
!$omp&  orbnu,qex1,qex2,rasc,rcl1,rcl2,rcoeff,rversion,slowdown,solcdec,solclat,solclatcdec,solsdec,&
!$omp&  solslat,solslatsdec,ssa1,ssa2,starbbtemp,starfile,starfilehr,th2oc,time4lwr,time4rad,&
!$omp&  time4swr,toffo3,tpofmt,tswr1,tswr2,tswr3,zcdayf,zdeclf,zmuz,zsolar1,zsolar2,zsolars)

      end module radmod

!
!     radiation subroutines
!

!     ===================
!     SUBROUTINE SOLARINI
!     ===================

      subroutine solarini
      use radmod
      use specblock
      
!       parameter(planckh = 6.62607004e-34)
!       parameter(boltzk = 1.38064852e-23 )
!       parameter(cc = 299792458.0        )
      parameter(const = 0.0143877735383)    !hc/k
      !parameter(chig0 = 11.234333860319996) !spectrum-weighted optical depth coefficient for 5772K
      
      real :: wv1(1024) !Wavelengths in meters up to 0.75 microns
      real :: wv2(1024) !Wavelength in meters starting at 0.75 microns
!     The 5772 K Rayleigh reference is tabulated on the grid BELOW, and wv1/wv2
!     are overwritten with the spectrum file's wavelengths further down when one
!     is given. Keeping the reference's own wavelengths here is what makes the
!     normalisation integrals self-consistent; before this they paired the
!     reference's values with the file's abscissae. SPEC-2.
      real :: wvg1(1024) !Reference wavelengths below 0.75 microns
      real :: wvg2(1024) !Reference wavelengths above 0.75 microns
      real :: wvm1(1024) !Wavelengths in microns up to 0.75 microns
      real :: wvm2(1024) !Wavelength in microns starting at 0.75 microns
      real :: bb1(1024) !Planck function for x<0.75 microns
      real :: bb2(1024) !Planck function for x>0.75 microns
      real :: bbg1(1024) !Planck function for x<0.75 microns
      real :: bbg2(1024) !Planck function for x>0.75 microns
      real :: bb3(965) !Planck function for albedo wavelengths
      real :: kdata(2048,2)
      real :: kdata2(965,2)
      
      
      real dl1,dl2,hinge,const1,const2,z1,z2,znet,wmin,lwmin,w1,w2,f1,f2,x
      real zout1,zout2 !Band flux lying outside the albedo grid
      integer k,nw,j
     
      if (mypid == NROOT) then
        
        constg = const/5772.0 !G star
        
        !wmin = const/(starbbtemp*36.841361) !Wavelength where exponential term is <=1.0e-16
        wmin = minwavel ! Set minimum wavelength to 316 nm; we don't include UV. 
                          ! This produces zsolar1=0.517 at Teff=5772 K.
        lwmin = log10(wmin)
        
        hinge = log10(7.5e-7) !We care about amounts above and below 0.75 microns
        dl1 = (hinge-lwmin)/1024.0
        dl2 = (-4-hinge)/1024.0
        
        do k=1,1024
          wv1(k) = 10**(lwmin+(k-1)*dl1)
          wv2(k) = 10**(hinge+(k-1)*dl2)
        enddo
        do k=1,1024
          wvm1(k) = (1.0e6 * wv1(k))**5
          wvm2(k) = (1.0e6 * wv2(k))**5
        enddo
        
        do k=1,1024
           wvg1(k) = wv1(k)
           wvg2(k) = wv2(k)
           bbg1(k) = 1.0/wvm1(k) * 1.0/(exp(constg/wv1(k))-1)
           bbg2(k) = 1.0/wvm2(k) * 1.0/(exp(constg/wv2(k))-1)
        enddo
        
        if (lstarfile) then ! Specific input spectrum was given
           call readdat(starfilehr,2,2048,kdata) !We keep the hi-res stuff for energy fractions
           wv1(:) = kdata(1:1024,1)*1.0e-6
           bb1(:) = kdata(1:1024,2)
           wv2(:) = kdata(1025:2048,1)*1.0e-6
           bb2(:) = kdata(1025:2048,2)
           do k=1,1024
              if (wv1(k) .lt. minwavel) bb1(k)=0. !Remove flux at wavelengths below 316 nm.
           enddo
           
           ! Scan through high-res wavelengths and re-sample to bb3 wavelengths
           call readdat(starfile,2,965,kdata2)
           bb3(:) = kdata2(:,2)
            
        else   ! Use blackbody spectrum
              
           !snowalbedos(:) = 0.25*(fsnowalb(:)+2.0*msnowalb(:)+csnowalb(:)) !assume mostly med-grain
           
!            const1 = 2*planckh*(cc**2)
           const2 = const/starbbtemp
           
           do k=1,1024 !Compute the Planck function
             bb1(k) = 1.0/wvm1(k) * 1.0/(exp(const2/wv1(k))-1) !const1/wv1(k)**5
             bb2(k) = 1.0/wvm2(k) * 1.0/(exp(const2/wv2(k))-1)
!              write(nud,*) wv1(k),bb1(k),wv2(k),bb2(k)
           enddo      !The scaling and units don't actually matter, because we're going to normalize
           
           do k=1,965 !Compute the Planck function for the wavelengths at which we have albedo data
             bb3(k) = 1.0/(wavelengths(k))**5 * 1.0/(exp(1.0e6*const2/wavelengths(k))-1)
           enddo
           
        endif
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
!           a1 = a1 + 0.5*(bb3(k)*fsnowalb(k)+bb3(k+1)*fsnowalb(k+1))* &
!      &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
            a1 = a1 + 0.5*(bb3(k)*iceblend(k)+bb3(k+1)*iceblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
!           a2 = a2 + 0.5*(bb3(k)*fsnowalb(k)+bb3(k+1)*fsnowalb(k+1)))* &
!      &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
            a2 = a2 + 0.5*(bb3(k)*iceblend(k)+bb3(k+1)*iceblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        
        z1 = 0.0
        z2 = 0.0
        
        zg1 = 0.0
        zg2 = 0.0
        
        zcross1 = 0.0
        zcross2 = 0.0
        
        zgcross1 = 0.0
        zgcross2 = 0.0
        
        do k=1,1023    !Do a trapezoidal integration above and below 0.75 microns
          z1 = z1 + 0.5*(bb1(k)+bb1(k+1))*(wv1(k+1)-wv1(k))
          z2 = z2 + 0.5*(bb2(k)+bb2(k+1))*(wv2(k+1)-wv2(k))
          zg1 = zg1 + 0.5*(bbg1(k)+bbg1(k+1))*(wvg1(k+1)-wvg1(k))
          zg2 = zg2 + 0.5*(bbg2(k)+bbg2(k+1))*(wvg2(k+1)-wvg2(k))
          zcross1 = zcross1 + 0.5*(bb1(k)/((wv1(k)*1.0e6)**4)+bb1(k+1)/((wv1(k+1)*1.0e6)**4)) &
     &                         *(wv1(k+1)-wv1(k))
          zcross2 = zcross2 + 0.5*(bb2(k)/((wv2(k)*1.0e6)**4)+bb2(k+1)/((wv2(k+1)*1.0e6)**4)) &
     &                         *(wv2(k+1)-wv2(k))
          zgcross1 = zgcross1+0.5*(bbg1(k)/((wvg1(k)*1.0e6)**4)+bbg1(k+1)/((wvg1(k+1)*1.0e6)**4)) &
     &                         *(wvg1(k+1)-wvg1(k))
          zgcross2 = zgcross2+0.5*(bbg2(k)/((wvg2(k)*1.0e6)**4)+bbg2(k+1)/((wvg2(k+1)*1.0e6)**4)) &
     &                         *(wvg2(k+1)-wvg2(k))
        enddo
        z1 = z1 + 0.5*(bb1(1024)+bb2(1))*(wv2(1)-wv1(1024))

!       The albedo grid does not span the model's bands. `wavelengths` runs
!       0.34 to 14.01 microns, while band 1 begins at minwavel and band 2 runs
!       to 100 microns. a1 and a2 below are integrated on the ALBEDO grid and
!       normalised by z1 and z2, which are integrated on the SPECTRUM grid, so
!       without the two terms added to every a1/a2 pair the quotient implicitly
!       assigns ZERO reflectance to the flux outside the albedo grid. For a K
!       dwarf that is 1.48 percent of band 1 and 0.06 percent of band 2, and it
!       is one-signed dark on every surface. Each blend is held at its endpoint
!       value across the gap, so numerator and denominator cover one interval.
        zout1 = 0.0
        zout2 = 0.0
        do k=1,1023
          if (wv1(k+1) .le. 1.0e-6*wavelengths(1)) then
            zout1 = zout1 + 0.5*(bb1(k)+bb1(k+1))*(wv1(k+1)-wv1(k))
          endif
          if (wv2(k) .ge. 1.0e-6*wavelengths(965)) then
            zout2 = zout2 + 0.5*(bb2(k)+bb2(k+1))*(wv2(k+1)-wv2(k))
          endif
        enddo
        zcross1 = zcross1+0.5*(bb1(1024)/((wv1(1024)*1.0e6)**4)+bb2(1)/((wv2(1)*1.0e6)**4)) &
     &                         *(wv2(1)-wv1(1024))
        zg1 = zg1 + 0.5*(bbg1(1024)+bbg2(1))*(wvg2(1)-wvg1(1024))
        zgcross1 = zgcross1+0.5*(bbg1(1024)/((wvg1(1024)*1.0e6)**4)+bbg2(1)/((wvg2(1)*1.0e6)**4)) &
     &                         *(wvg2(1)-wvg1(1024))
        
        zg = zg1+zg2
        zgcross = zgcross1 + zgcross2
        zchi = zgcross / zg !spectrum-weighted cross section for 5772 K
        rcoeff = (zcross1 + zcross2) * zsolar1 / z1 / zchi !Using default zsolar=0.517 here
        
        ! effective optical depth is the spectral average of the cross-section, normalized to 
        ! 5772 K input blackbody. There's already a spectral dependence due to z1/z2 partitioning,
        ! so we compute the true weighting and normalize to the partitioning and solar result
!         
!         We want tau = <sigma>/<sigma_g>*tau_g, where 
!         
!                        int_0^inf[F(w) w^-4 dw] 
!             <sigma> = -------------------------
!                          int_0^inf[F(w) dw]    
!                          
!         so:
!         
!                int_0^inf[F(w) w^-4 dw]        int_0^inf[F_g(w) dw]
!         tau = ------------------------- x --------------------------- x tau_g
!                  int_0^inf[F(w) dw]        int_0^inf[F_g(w) w^-4 dw]
!         
!         We need to somehow account for the fact that we have two bands, especially because that
!         will impart a Z1/Z1_g scaling all on its own. We can do this by multiplying by 1:
!         
!                int_0^w2[F(w) dw]     int_0^inf[F(w) w^-4 dw]        int_0^inf[F_g(w) dw]
!         tau = ------------------- x ------------------------- x -------------------------- x tau_g
!                int_0^w2[F(w) dw]       int_0^inf[F(w) dw]        int_0^inf[F_g(w) w^-4 dw]
!                
!         when we rearrange:
!          
!                int_0^w2[F(w) dw]      int_0^inf[F(w) w^-4 dw]        int_0^inf[F_g(w) dw]
!         tau = -------------------- x ------------------------ x -------------------------- x tau_g
!                int_0^inf[F(w) dw]        int_0^w2[F(w) dw]        int_0^inf[F_g(w) w^-4 dw]       
!         
!         This new first term out front is equal to Z1, the partitioning fraction. So 
!          
!                      int_0^inf[F(w) w^-4 dw]        int_0^inf[F_g(w) dw]
!         tau =  Z1 x ------------------------- x --------------------------- x tau_g
!                         int_0^w2[F(w) dw]        int_0^inf[F_g(w) w^-4 dw]       
!                
!         We also know that PlaSim's energy partitioning scheme will impart a factor of Z1/Z1_g, so
!         if we know what we really have is
!         
!                    Z1
!         tau = R x ---- x tau_g
!                   Z1_g
!         
!         then we can solve for R:
!         
!                      int_0^inf[F(w) w^-4 dw]        int_0^inf[F_g(w) dw]
!         R =  Z1_g x ------------------------- x --------------------------- 
!                         int_0^w2[F(w) dw]        int_0^inf[F_g(w) w^-4 dw]  
!
!                        zcross1 + zcross2          zg1 + zg2
!           = zsolar1 x ------------------- x ---------------------
!                               z1             zgcross1 + zgcross2
!                       
        zdenom1 = 0.01/z1
        zdenom2 = 0.01/z2
        
        a1 = a1 + iceblend(1)*zout1   !Flux below the albedo grid
        a2 = a2 + iceblend(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        znet = z1+z2
        
        z1 = z1/znet
        z2 = 1.0-z1
        
        zsolar1 = z1
        zsolar2 = z2
        
        write(nud,*) "Energy fraction below 0.75 microns:",zsolar1
        write(nud,*) "Energy fraction above 0.75 microns:",zsolar2
        write(nud,*) "Rayleigh scattering coefficient:",rcoeff
        
        zsolars(1) = zsolar1
        zsolars(2) = zsolar2
        
        dsnowalb(1) = a1
        dsnowalb(2) = a2
        
        write(nud,*) "Snow albedo below 0.75 microns:",dsnowalb(1)
        write(nud,*) "Snow albedo above 0.75 microns:",dsnowalb(2)
        write(nud,*) "Overall snow albedo:",z1*dsnowalb(1)+z2*dsnowalb(2)
        
        if (nsimplealbedo>0.5) dsnowalb(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*iceblendmin(k)+bb3(k+1)*iceblendmin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*iceblendmin(k)+bb3(k+1)*iceblendmin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + iceblendmin(1)*zout1   !Flux below the albedo grid
        a2 = a2 + iceblendmin(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dsnowalbmn(1) = a1
        dsnowalbmn(2) = a2
        
        write(nud,*) "Minimum snow albedo below 0.75 microns:",dsnowalbmn(1)
        write(nud,*) "Minimum snow albedo above 0.75 microns:",dsnowalbmn(2)
        write(nud,*) "Overall minimum snow albedo:",z1*dsnowalbmn(1)+z2*dsnowalbmn(2)
        
        if (nsimplealbedo>0.5) dsnowalbmn(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*iceblendmax(k)+bb3(k+1)*iceblendmax(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*iceblendmax(k)+bb3(k+1)*iceblendmax(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + iceblendmax(1)*zout1   !Flux below the albedo grid
        a2 = a2 + iceblendmax(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dsnowalbmx(1) = a1
        dsnowalbmx(2) = a2
        
        write(nud,*) "Maximum snow albedo below 0.75 microns:",dsnowalbmx(1)
        write(nud,*) "Maximum snow albedo above 0.75 microns:",dsnowalbmx(2)
        write(nud,*) "Overall maximum snow albedo:",z1*dsnowalbmx(1)+z2*dsnowalbmx(2)
                
        if (nsimplealbedo>0.5) dsnowalbmx(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*seaicemin(k)+bb3(k+1)*seaicemin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*seaicemin(k)+bb3(k+1)*seaicemin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + seaicemin(1)*zout1   !Flux below the albedo grid
        a2 = a2 + seaicemin(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dicealbmn(1) = a1
        dicealbmn(2) = a2
        
        write(nud,*) "Minimum sea ice albedo below 0.75 microns:",dicealbmn(1)
        write(nud,*) "Minimum sea ice albedo above 0.75 microns:",dicealbmn(2)
        write(nud,*) "Overall minimum sea ice albedo:",z1*dicealbmn(1)+z2*dicealbmn(2)
        
        if (nsimplealbedo>0.5) dicealbmn(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*seaicemax(k)+bb3(k+1)*seaicemax(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*seaicemax(k)+bb3(k+1)*seaicemax(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + seaicemax(1)*zout1   !Flux below the albedo grid
        a2 = a2 + seaicemax(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dicealbmx(1) = a1
        dicealbmx(2) = a2
        
        write(nud,*) "Maximum sea ice albedo below 0.75 microns:",dicealbmx(1)
        write(nud,*) "Maximum sea ice albedo above 0.75 microns:",dicealbmx(2)
        write(nud,*) "Overall maximum sea ice albedo:",z1*dicealbmx(1)+z2*dicealbmx(2)
        
        if (nsimplealbedo>0.5) dicealbmx(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*glacalbmin(k)+bb3(k+1)*glacalbmin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*glacalbmin(k)+bb3(k+1)*glacalbmin(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + glacalbmin(1)*zout1   !Flux below the albedo grid
        a2 = a2 + glacalbmin(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dglacalbmn(1) = a1
        dglacalbmn(2) = a2
        
        write(nud,*) "Minimum glacier albedo below 0.75 microns:",dglacalbmn(1)
        write(nud,*) "Minimum glacier albedo above 0.75 microns:",dglacalbmn(2)
        write(nud,*) "Overall minimum glacier albedo:",z1*dglacalbmn(1)+z2*dglacalbmn(2)
        
        if (nsimplealbedo>0.5) dglacalbmn(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*groundblend(k)+bb3(k+1)*groundblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*groundblend(k)+bb3(k+1)*groundblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + groundblend(1)*zout1   !Flux below the albedo grid
        a2 = a2 + groundblend(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        dgroundalb(1) = a1
        dgroundalb(2) = a2
        
        write(nud,*) "Ground albedo below 0.75 microns:",dgroundalb(1)
        write(nud,*) "Ground albedo above 0.75 microns:",dgroundalb(2)
        write(nud,*) "Overall ground albedo:",z1*dgroundalb(1)+z2*dgroundalb(2)
        
        if (nsimplealbedo>0.5) dgroundalb(:) = z1*a1 + z2*a2
        
        a1 = 0.0
        a2 = 0.0
        do k=1,41 !Compute insolation-weighted albedo below 0.75 microns
            a1 = a1 + 0.5*(bb3(k)*oceanblend(k)+bb3(k+1)*oceanblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        do k=42,964 !Compute insolation-weighted albedo above 0.75 microns
            a2 = a2 + 0.5*(bb3(k)*oceanblend(k)+bb3(k+1)*oceanblend(k+1))* &
       &              1.0e-6*(wavelengths(k+1)-wavelengths(k))
        enddo
        a1 = a1 + oceanblend(1)*zout1   !Flux below the albedo grid
        a2 = a2 + oceanblend(965)*zout2 !Flux above the albedo grid
        a1 = zdenom1*a1 !Percent -> Decimal; normalization
        a2 = zdenom2*a2
        
        doceanalb(1) = a1
        doceanalb(2) = a2
        
        write(nud,*) "Ocean albedo below 0.75 microns:",doceanalb(1)
        write(nud,*) "Ocean albedo above 0.75 microns:",doceanalb(2)
        write(nud,*) "Overall ocean albedo:",z1*doceanalb(1)+z2*doceanalb(2)
        
        if (nsimplealbedo>0.5) doceanalb(:) = z1*a1 + z2*a2
        
        
!     The nine put_restart_array calls that stood here wrote to nwriunit (34)
!     from solarini, which runs at INITIALISATION when no restart file is open
!     on that unit, so Fortran connected fort.34 and put them there. Nothing
!     ever read that file. They were the orphaned write half of a round trip
!     whose read half is commented out above (the get_restart_array block at
!     nstarfile > 0), where solarini was made unconditional: these arrays are
!     recomputed from the namelist at every start and are not checkpointed.
!     CLIM-38, notes/audits/zsolars-restart-overread.md.
                               
        
      endif
      
      
      
      call mpbcrn(zsolars,2)
      call mpbcr(zsolar1)
      call mpbcr(zsolar2)
      call mpbcr(rcoeff)
      call mpbcrn(dsnowalb,2)
      call mpbcrn(dsnowalbmn,2)
      call mpbcrn(dsnowalbmx,2)
      call mpbcrn(dglacalbmn,2)
      call mpbcrn(dicealbmn,2)
      call mpbcrn(dicealbmx,2)
      call mpbcrn(dgroundalb,2)
      call mpbcrn(doceanalb,2)
      
!       call mpputgp('zsolars',zsolars,2,1)
!       call mpputgp('dsnowalb',dsnowalb,2,1)
!       call mpputgp('dsnowalbmn',dsnowalbmn,2,1)
!       call mpputgp('dsnowalbmx',dsnowalbmx,2,1)
!       call mpputgp('dicealbmn',dicealbmn,2,1)
!       call mpputgp('dicealbmx',dicealbmx,2,1)
!       call mpputgp('dglacalbmn',dglacalbmn,2,1)
!       call mpputgp('dgroundalb',dgroundalb,2,1)
!       call mpputgp('doceanalb',doceanalb,2,1)
      
      end subroutine solarini
      
!     =================
!     SUBROUTINE RADINI
!     =================

      subroutine radini
      use radmod
!
      logical :: lexaero        ! does the aerosol optics file exist
      integer :: jaer           ! aerosol species index
      character (len=16) :: ysurf  ! surface field name for a species
!
!     initialize radiation
!     this *sub* is called by PUMA (PUMA-interface)
!
!     this *sub* reads the radiation namelist *radmod_nl*
!     and broadcasts the parameters
!
!     the following PUMA *subs* are used:
!
!     mpbci  : broadcasts 1 integer
!     mpbcr  : broadcasts 1 real
!
!     the following PUMA variables are used:
!
!     mypid  : process id (used for mpp)
!     nroot  : id of root process (used for mpp)
!
!**   0) define namelist
!
      namelist/radmod_nl/ndcycle,ncstsol,solclat,solcdec,no3,co2,ch4,n2o &
     &               ,iyrbp,nswr,nlwr,nfixed,slowdown,nradice,npbroaden,desync    &
     &               ,o3uvw,o3visw,h2osww,h2oswl,co2sww   &
     &               ,a0o3,a1o3,aco3,bo3,co3,toffo3,o3scale,newrsc,necham,necham6   &
     &               ,nsol,nclouds,nswrcl,nrscat,rcl1,rcl2,acl2,clgray,tpofmt   &
     &               ,acllwr,tswr1,tswr2,tswr3,th2oc,dawn,starbbtemp,nstartemp  &
     &               ,elwland,elwsea                                            &
     &               ,nsimplealbedo,nstarfile,starfile,starfilehr,minwavel      &
     &               ,ndustrad,dustsc,dusthsc,dustqlw,aerofile,aeroqlw          &
     &               ,nsolcycle,gsolstart,gsolamp,gsolperiod,gsolphase   &
     &               ,gsolamp2,gsolperiod2,gsolphase2
!
!     namelist parameter:
!
!     ndcycle : switch for daily cycle 1=on/0=off
!     ncstsol : switch to set constant insolation 
!     solclat : constant cosine of latitude of insolation 
!     solcdec : constant solar declination
!     no3     : switch for ozon 1=on/0=off
!     co2     : co2 concentration (ppmv)
!     iyrbp   : Year before present (1950 AD); default = 2000 AD
!     nswr    : switch for short wave radiation (dbug) 1=on/0=off
!     nlwr    : switch for long wave radiation (dbug) 1=on/0=off
!     nsol    : switch for solar insolation (dbug) 1=on/0=off
!     nswrcl  : switch for computed or prescribed cloud props. 1=com/0=pres
!     nrscat  : switch for rayleigh scattering (dbug) 1=on/0=off
!     o3scale : factor for scaling o3
!     rcl1(3) : cloud albedos spectral range 1
!     rcl2(3) : cloud albedos spectral range 2
!     acl2(3) : cloud absorptivities spectral range 2
!     clgray  : cloud grayness
!     tpofmt  ! tuning of point of mean (lwr) transmissivity in layer
!     acllwr  ! mass absorption coefficient for clouds (lwr)
!     tswr1   ! tuning of cloud albedo range1
!     tswr2   ! tuning of cloud back scattering c. range2
!     tswr3   ! tuning of cloud s. scattering alb. range2
!     th2oc   ! absorption coefficient for h2o continuum
!     elwland : surface longwave emissivity over land (1)
!     elwsea  : surface longwave emissivity over ocean and sea ice (1)
!     dawn    : zenith angle threshhold for night
!
!     aerosol namelist parameters:
!     l_aerorad  : turns aerosol scattering/absorption effects on/off
!     aerofile   : dat file containing optical constants for aerosols
!
!     following parameters are read from the planet module
!
!     gsol0   : solar constant (w/m2)
!
!     NO TRUNCATION IS A SPECIAL CASE, and the table that made three of them
!     one never fired. Upstream carried a per-(NTRU, NLEV) shortwave tuning
!     here -- tswr1, tswr2, tswr3 and th2oc, selected through a flag named
!     jtune, with branches for T21/T1, T31 and T42 -- and every branch of it
!     was unreachable at every truncation. Each one reached its coefficients
!     only when ndcycle was not 1; ndcycle's compiled default is 1 (:205) and
!     read(11,radmod_nl) is below this point, so ndcycle held its compiled
!     value whenever the test was made and the flag was always 0. The block
!     ended by announcing its own failure, on every run at every truncation:
!     'No radiation setup for this resolution ... you may need to tune'.
!
!     What the modelled radiation uses now is what it has always used: the
!     module defaults at :72-75, and whatever the namelist sets over them.
!     This project sets TSWR3 from model.cloud_absorption_scale and leaves
!     tswr1, tswr2 and th2oc at those defaults.
!     world-ys9, world-677x; exoplasim/notes/resolution-tuned-parameters.md.
!
!**   1) read and print version & namelist parameters
!
!     A LIVE REMNANT OF AN EARTH CALENDAR, and inert. iyrbp is a year before
!     1950 AD, and it reaches nothing here: it feeds orb_params, which is
!     Berger's Milankovitch series for EARTH's orbital elements, and radini
!     reaches that only at nfixorb == 0. run_exoplasim.py passes fixedorbit=True
!     at every call, so nfixorb is 1 and the series is unreachable. Throwing
!     that switch would compute this world's eccentricity, obliquity and
!     longitude of perihelion from Earth's polynomial fits at a year AD.
!     world-9d1.
      iyrbp = 1950 - n_start_year

      if (mypid==NROOT) then
         open(11,file=radmod_namelist)
         read(11,radmod_nl)
         close(11)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * RADMOD ",a34," *")') trim(rversion)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist RADMOD_NL from <radmod_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,radmod_nl)
         
         minwavel = minwavel*1.0e-9
         
         oldfixedlon = fixedlon
         if (nrestart > 0.) call get_restart_real('fixedlon',fixedlon)
         if ((fixedlon .ne. oldfixedlon) .and. (desync .eq. 0.0)) fixedlon=oldfixedlon
         
         if ((necham.eq.1).and.(necham6.eq.1)) necham=0 !necham6 overrides necham
      endif ! (mypid==NROOT)
!
!     broadcast namelist parameter
!
      call mpbci(ndcycle)
      call mpbci(ncstsol)
      call mpbci(no3)
      call mpbci(nfixed)
      call mpbcr(fixedlon)
      call mpbcr(desync)
      call mpbcr(slowdown)
      call mpbcr(a0o3)
      call mpbcr(a1o3)
      call mpbcr(aco3)
      call mpbcr(bo3)
      call mpbcr(co3)
      call mpbcr(toffo3)
      call mpbcr(o3uvw)
      call mpbcr(o3visw)
      call mpbcr(h2osww)
      call mpbcr(h2oswl)
      call mpbcr(co2sww)
      call mpbcr(elwland)
      call mpbcr(elwsea)
      call mpbcr(o3scale)
      call mpbcr(co2)
!
!     CH4 AND N2O MUST BE BROADCAST, and forgetting it is invisible.
!
!     radmod_nl is read on NROOT only, so a namelist variable that is not
!     broadcast keeps its default on every other rank. For these two the
!     default is 0.0, meaning absent -- so the band simply does not run on
!     ranks 1 and up, and at T42 on 8 ranks that left the term acting on the
!     eight polar latitude rows and nowhere else. It cost a day of looking for
!     a physics explanation for a term that came out ninety times too weak.
!
!     The same class is already recorded twice in this project:
!     notes/audits/nlowio-collective-deadlock.md is a collective placed behind
!     an unbroadcast nlowio, and PHYS-9 is a namelist key applied at prepare
!     and not per segment. CLIM-42.
      call mpbcr(ch4)
      call mpbcr(n2o)
      call mpbcr(gsol0)
      call mpbci(nsolcycle)
      call mpbci(gsolstart)
      call mpbcr(gsolamp)
      call mpbcr(gsolperiod)
      call mpbcr(gsolphase)
      call mpbcr(gsolamp2)
      call mpbcr(gsolperiod2)
      call mpbcr(gsolphase2)
      call mpbcr(solclat)
      call mpbcr(solcdec)
      call mpbcr(clgray)
      call mpbcr(dawn)
      call mpbcr(th2oc)
      call mpbcr(tpofmt)
      call mpbcr(acllwr)
      call mpbcr(tswr1)
      call mpbcr(tswr2)
      call mpbcr(tswr3)
      call mpbcrn(rcl1,3)
      call mpbcrn(rcl2,3)
      call mpbcrn(acl2,3)
      call mpbci(iyrbp)
      call mpbci(nswr)
      call mpbci(nlwr)
      call mpbci(nsol)
      call mpbci(nrscat)
      call mpbci(nswrcl)
      call mpbci(nclouds)
      call mpbci(necham)
      call mpbci(necham6)
      call mpbci(nradice)
      call mpbci(npbroaden)

      call mpbcr(starbbtemp)
      call mpbci(nstartemp)
      call mpbci(nsimplealbedo)
      call mpbci(nstarfile)
      call mpbcr(minwavel)
      
      call mpbci(l_aerorad)
!
!     aero_ini has already run, on NROOT only, and has copied aero_nl's
!     particle radius into radmod's own. This is the broadcast that was
!     missing: without it every rank kept the 50 nm default and the
!     shortwave aerosol optical depth was (50e-9/apart)**2 of intent.
!     Harmless where aero_ini never ran, because both copies are then the
!     same default.
!
      call mpbcr(apart)

      call mpbci(ndustrad)
      call mpbcrn(dustsc,NAERSP)
      call mpbcrn(dusthsc,NAERSP)
      call mpbcrn(dustqlw,NAERSP)
      call mpbcr(aeroqlw)

!      
!     determine stellar parameters      
!
      if (nstarfile > 0) then
!         if (nrestart > 0.) then
!           if (mypid == NROOT) then
!             call get_restart_array("zsolars",zsolars,2,2,1)
!             call get_restart_array('dsnowalb',dsnowalb,2,2,1)
!             call get_restart_array('dsnowalbmn',dsnowalbmn,2,2,1)
!             call get_restart_array('dsnowalbmx',dsnowalbmx,2,2,1)
!             call get_restart_array('dicealbmn',dicealbmn,2,2,1)
!             call get_restart_array('dicealbmx',dicealbmx,2,2,1)
!             call get_restart_array('dglacalbmn',dglacalbmn,2,2,1)
!             call get_restart_array('dgroundalb',dgroundalb,2,2,1)
!             call get_restart_array('doceanalb',doceanalb,2,2,1)
!             zsolar1 = zsolars(1)
!             zsolar2 = zsolars(2)
!             write(nud,*) "Read zsolar1 from restart: ",zsolar1
!             write(nud,*) "Read zsolar2 from restart: ",zsolar2
!             write(nud,*) "Read snow albedo <0.75 um from restart: ",dsnowalb(1)
!             write(nud,*) "Read snow albedo >0.75 um from restart: ",dsnowalb(2)
!             write(nud,*) "Read snow min albedo <0.75 um from restart: ",dsnowalbmn(1)
!             write(nud,*) "Read snow min albedo >0.75 um from restart: ",dsnowalbmn(2)
!             write(nud,*) "Read snow max albedo <0.75 um from restart: ",dsnowalbmx(1)
!             write(nud,*) "Read snow max albedo >0.75 um from restart: ",dsnowalbmx(2)
!             write(nud,*) "Read glacier min albedo <0.75 um from restart: ",dglacalbmn(1)
!             write(nud,*) "Read glacier min albedo >0.75 um from restart: ",dglacalbmn(2)
!             write(nud,*) "Read sea ice min albedo <0.75 um from restart: ",dicealbmn(1)
!             write(nud,*) "Read sea ice min albedo >0.75 um from restart: ",dicealbmn(2)
!             write(nud,*) "Read sea ice max albedo <0.75 um from restart: ",dicealbmx(1)
!             write(nud,*) "Read sea ice max albedo >0.75 um from restart: ",dicealbmx(2)
!             write(nud,*) "Read ground albedo <0.75 um from restart: ",dgroundalb(1)
!             write(nud,*) "Read ground albedo >0.75 um from restart: ",dgroundalb(2)
!             write(nud,*) "Read ocean albedo <0.75 um from restart: ",doceanalb(1)
!             write(nud,*) "Read ocean albedo >0.75 um from restart: ",doceanalb(2)
!           endif
!           call mpbcr(zsolar1)
!           call mpbcr(zsolar2)
!           call mpbcrn(dsnowalb,2)
!           call mpbcrn(dsnowalbmn,2)
!           call mpbcrn(dsnowalbmx,2)
!           call mpbcrn(dglacalbmn,2)
!           call mpbcrn(dicealbmn,2)
!           call mpbcrn(dicealbmx,2)
!           call mpbcrn(dgroundalb,2)
!           call mpbcrn(doceanalb,2)
! !           call mpputgp(
!         else
!           call solarini 
!         endif
        lstarfile = .true.
        call solarini
        nstartemp = 1
        call mpbci(nstartemp)
      else if (nstartemp > 0) then
        call solarini
      else
        call mpbcr(zsolar1)
        call mpbcr(zsolar2)
        call mpbcr(rcoeff)
        call mpbcrn(dsnowalb,2)
        call mpbcrn(dsnowalbmn,2)
        call mpbcrn(dsnowalbmx,2)
        call mpbcrn(dglacalbmn,2)
        call mpbcrn(dicealbmn,2)
        call mpbcrn(dicealbmx,2)
        call mpbcrn(dgroundalb,2)
        call mpbcrn(doceanalb,2)
      endif
      
!
!     determine orbital parameters
!

      iyrad = 1950 - iyrbp
      if (nfixorb == 1) then ! fixed orbital params (default AMIP II)
         iyrad = ORB_UNDEF_INT
      endif
      call orb_params(iyrad, eccen, obliq, meananom0, mvelp                          &
     &               ,obliqr, meananom0r, lambm0, mvelpp, log_print, ngenkeplerian &
     &               ,mypid, nroot,nud)
     
     
     call mpbcr(meananom0r)
     call mpbci(ngenkeplerian)

!
!     read climatological ozone
!
      if (no3 == 2) then
         allocate(dqo3cl(NHOR,NLEV,0:13))
         dqo3cl(:,:,:) = 0.0
         call mpsurfgp('dqo3cl',dqo3cl,NHOR,NLEV*14)
      endif
!
!     set co2 3d-field (enable external co2 by if statement)
!
!
      if(co2 > 0.) then
       dqco2(:,:)=co2
      endif
      
!
!     How many aerosol species this run carries, asked ONCE and before the
!     optics are read, because the aerofile's column count is this number.
!     Prescribed species occupy 1..ndustrad and the transported tracer, if it
!     is on, is ndustrad+1.
!
      iaerint = 0
      if (l_aero > 0 .and. l_aerorad == 1) iaerint = 1
      if (ndustrad < 0) call mpabort('ndustrad must not be negative')
      naerosp = ndustrad + iaerint
      if (naerosp > NAERSP) then
       if (mypid == NROOT) then
        write(nud,*) 'aerosol species requested: ',naerosp
        write(nud,*) 'NAERSP in radmod is        : ',NAERSP
        write(nud,*) 'Raise NAERSP and rebuild every binary (CLAUDE.md rule 4).'
       endif
       call mpabort('more aerosol species than NAERSP')
      endif
!
      if (naerosp > 0) then
       if (mypid == NROOT) then
!
!     readdat opens with the default status, so a missing aerofile is CREATED
!     empty and the failure arrives as an end-of-file inside a utility rather
!     than as a statement about the aerosol. Say what is wrong instead.
!
!     `aerofile` is settable from radmod_nl as well as aero_nl, because the
!     prescribed-dust path does not run aero_ini: that is only called when the
!     semi-Lagrangian tracer transport is on, and nothing here is transported.
!     aero_ini runs BEFORE radini, so radmod_nl wins when both name a file.
!
        inquire(file=aerofile,exist=lexaero)
        if (.not. lexaero) then
         write(nud,*) 'aerosol optics file not found: ',trim(aerofile)
         call mpabort('aerofile is missing')
        endif
!
!     One COLUMN per species. With naerosp = 1 this reads exactly the file the
!     single-species path always read, so an existing aerofile stays valid.
!
        call readdat(aerofile,naerosp,8,aeroqs) ! Get Qextinction, Qscattering, Qbackscatter, g for band 1 & 2

        do jaer = 1,naerosp
         if (aeroqs(1,jaer) <= 0. .or. aeroqs(5,jaer) <= 0.) then
          write(nud,*) 'aerofile species ',jaer,' has non-positive Qext'
          write(nud,*) 'band 1: ',aeroqs(1,jaer),' band 2: ',aeroqs(5,jaer)
          write(nud,*) 'Expected ',naerosp,' column(s) in ',trim(aerofile)
          call mpabort('aerofile has too few columns or a zero Qext')
         endif
         ssa1(jaer) = aeroqs(2,jaer)/aeroqs(1,jaer) ! Single scattering albedo band 1 (qscat/qext)
         ssa2(jaer) = aeroqs(6,jaer)/aeroqs(5,jaer) ! Single scattering albedo band 2
         bscat1(jaer) = aeroqs(3,jaer)/aeroqs(2,jaer) ! Backscatter ratio band 1
         bscat2(jaer) = aeroqs(7,jaer)/aeroqs(6,jaer) ! Backscatter ratio band 2
         qex1(jaer) = aeroqs(1,jaer) ! Extinction efficiency band 1
         qex2(jaer) = aeroqs(5,jaer) ! Extinction efficiency band 2
        enddo
       endif

        call mpbcrn(ssa1,NAERSP) ! Broadcast optical constants
        call mpbcrn(ssa2,NAERSP)
        call mpbcrn(bscat1,NAERSP)
        call mpbcrn(bscat2,NAERSP)
        call mpbcrn(qex1,NAERSP)
        call mpbcrn(qex2,NAERSP)
      endif
!
!     prescribed dust: the column optical depth field, and the two things that
!     make enabling it without them a silent wrong answer rather than a loud one
!
      ddustcol(:,:) = 0.
      do jaer = 1,ndustrad
       if (dustqlw(jaer) <= 0.) then
        if (mypid == NROOT) then
         write(nud,*) 'PRESCRIBED AEROSOL species ',jaer,': dustqlw is ',dustqlw(jaer)
         write(nud,*) 'The aerosol acts in the two SHORTWAVE bands only unless'
         write(nud,*) 'a thermal-infrared absorption ratio is supplied, and a'
         write(nud,*) 'shortwave-only aerosol is a larger error than none.'
        endif
        call mpabort('every prescribed aerosol species requires dustqlw > 0')
       endif
       write(ysurf,'("ddustcol",i0)') jaer
       if (jaer == 1) ysurf = 'ddustcol'   ! code 1811 keeps the name it has
       call mpsurfgp(trim(ysurf),ddustcol(1,jaer),NHOR,1)
       call mpmaxval(ddustcol(1,jaer),NHOR,1,zdustmx)
       if (zdustmx <= 0.) then
        if (mypid == NROOT) then
         write(nud,*) 'PRESCRIBED AEROSOL: no field was read for species ',jaer
         write(nud,*) 'Expected surface code ',1810+jaer,' in the run directory.'
        endif
        call mpabort('a prescribed aerosol species has no surface field')
       endif
       if (mypid == NROOT) then
        write(nud,'(/," *********************************************")')
        write(nud,'(" * PRESCRIBED AEROSOL species ",i1," is ON         *")') jaer
        write(nud,'(" *********************************************")')
        write(nud,*) 'surface code                    ',1810+jaer
        write(nud,*) 'max band-1 column optical depth ',zdustmx
        write(nud,*) 'scale factor                    ',dustsc(jaer)
        write(nud,*) 'scale height (m)                ',dusthsc(jaer)
        write(nud,*) 'thermal-IR absorption ratio     ',dustqlw(jaer)
        write(nud,*) 'single scattering albedo band 1 ',ssa1(jaer)
        write(nud,*) 'band 2 / band 1 extinction      ',qex2(jaer)/qex1(jaer)
       endif
      enddo
!
!     interactive aerosol: the one way of enabling it that is silently wrong
!     rather than loudly wrong. It is no longer exclusive with a prescribed
!     column: it is species ndustrad+1 of the same array.
!
      if (iaerint == 1 .and. aeroqlw <= 0.) then
       if (mypid == NROOT) then
        write(nud,*) 'INTERACTIVE AEROSOL: aeroqlw is ',aeroqlw
        write(nud,*) 'The aerosol acts in the two SHORTWAVE bands only unless'
        write(nud,*) 'a thermal-infrared absorption ratio is supplied, and a'
        write(nud,*) 'shortwave-only aerosol is a larger error than none.'
       endif
       call mpabort('l_aerorad=1 requires aeroqlw > 0')
      endif
      if (iaerint == 1 .and. mypid == NROOT) then
       write(nud,'(/," *********************************************")')
       write(nud,'(" * INTERACTIVE AEROSOL RADIATION is ON       *")')
       write(nud,'(" *********************************************")')
       write(nud,*) 'species index                   ',ndustrad+1
       write(nud,*) 'particle radius (m)             ',apart
       write(nud,*) 'thermal-IR absorption ratio     ',aeroqlw
       write(nud,*) 'band 2 / band 1 extinction      ',qex2(ndustrad+1)/qex1(ndustrad+1)
      endif
!
!     Gather the thermal-IR absorption ratios into one per-species array, so
!     the longwave sums over species instead of branching on which path each
!     came from. dustqlw and aeroqlw stay separate in the namelist because the
!     two paths genuinely carry different particles: the prescribed field is
!     this world's aerosol at the offline chain's size distribution, and the
!     transported tracer is whatever aero_nl's apart and rhop describe.
!
      aqlw(:) = 0.
      do jaer = 1,ndustrad
       aqlw(jaer) = dustqlw(jaer)
      enddo
      if (iaerint == 1) aqlw(ndustrad+1) = aeroqlw
!
      if (naerosp > 0 .and. mypid == NROOT) then
       write(nud,*) 'AEROSOL SPECIES CARRIED: ',naerosp,' of NAERSP ',NAERSP
      endif
!
      return
      end subroutine radini

!     ==================
!     SUBROUTINE RADSTEP
!     ==================

      subroutine radstep
      use radmod
!
!     do the radiation calculations
!     this *sub* is called by PUMA (PUMA-interface)
!
!     no PUMA *subs* are used
!
!     the following PUMA variables are used/modified:
!
!     ga               : gravity accelleration (m/s2) (used)
!     acpd             : specific heat of dry air (J/kgK) (used)
!     ADV              : ACPV/acpd - 1  (used)
!     sigma(NLEV)      : sigma of T-levels (used)
!     dp(NHOR)         : surface pressure (Pa) (used)
!     dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
!     dtdt(NHOR,NLEP)  : temperature tendencies (K/s) (modified)
!     dswfl(NHOR,NLEP) : short wave radiation (W/m2)  (modified)
!     dlwfl(NHOR,NLEP) : long wave radiation (W/m2)   (modified)
!     dfu(NHOR,NLEP)   : short wave radiation upward (W/m2) (modified)
!     dfd(NHOR,NLEP)   : short wave radiation downward (W/m2) (modified)
!     dftu(NHOR,NLEP)  : long wave radiation upward (W/m2) (modified)
!     dftd(NHOR,NLEP)  : long wave radiation downward (W/m2) (modified)
!     dflux(NHOR,NLEP) : total radiation (W/m2) (modified)
!
!     the following radiation *subs* are called:
!
!     solang           : calc. cosine of solar zenit angle
!     mko3             : calc. ozon distribution
!     swr              : calc. short wave radiation fluxes
!     lwr              : calc. long wave radiation fluxes
!
!
!**   0) define local arrays
!

      real zdtdt(NHOR,NLEV)    ! temperature tendency due to rad (K/s)
      real zdh(NHOR,NLEV)      ! Thickness of an atmospheric layer (m)
      real zfice(NHOR)         ! Temporary backup sea ice array
!
!     allocatable arrays for diagnostic
!

      real, allocatable :: zprf1(:,:)
      real, allocatable :: zprf2(:,:)
      real, allocatable :: zprf3(:,:)
      real, allocatable :: zprf4(:,:)
      real, allocatable :: zprf5(:,:)
      real, allocatable :: zprf6(:,:)
      real, allocatable :: zprf7(:,:)
      real, allocatable :: zprf8(:,:)
      real, allocatable :: zprf9(:,:)
      real, allocatable :: zprf10(:,:)
      real, allocatable :: zprf11(:,:)
      real, allocatable :: zprf12(:,:)
      real, allocatable :: zcc(:,:)
      real, allocatable :: zalb1(:)
      real, allocatable :: zalb2(:)
      real, allocatable :: zdtdte(:,:)
!
!     cpu time estimates
!
      if(ntime == 1) call mksecond(zsec,0.)
!
!**   1) set all fluxes to zero
!

      dfu(:,:)   = 0.0         ! short wave radiation upward
      dfd(:,:)   = 0.0         ! short wave radiation downward
      dfdsw1(:)  = 0.0         ! surface short wave downward, band 1
      dfdsw2(:)  = 0.0         ! surface short wave downward, band 2
      dftu(:,:)  = 0.0         ! long wave radiation upward
      dftd(:,:)  = 0.0         ! long wave radiation downward
      dswfl(:,:) = 0.0         ! total short wave radiation
      dlwfl(:,:) = 0.0         ! total long wave radiation
      dftue1(:,:)= 0.0         ! entropy
      dftue2(:,:)= 0.0         ! entropy
      
      if (nradice==0) then
        zfice(:) = dicec(:)
        dicec(:) = 0.0
      endif
!
!**   2) compute cosine of solar zenit angle for each gridpoint
!
      if(nsol==1) call solang
!
!**   3) compute ozon distribution
!
      if(no3<3) call mko3
!
!**   3b) distribute the prescribed dust column over the layers
!
      if(ndustrad >= 1) call dustprof
!
!**   3c) build the interactive aerosol's optical depth per layer
!
      if(iaerint == 1) call aeroprof
!
!**   3d) gather every species into one per-layer array for swr and lwr
!
      if(naerosp > 0) call aerogather
!
!**   4) short wave radiation
!
!     a) if clear sky diagnostic is switched on:
!

      if(ndiagcf > 0) then
       allocate(zcc(NHOR,NLEP))
       allocate(zalb1(NHOR))
       allocate(zalb2(NHOR))
       zcc(:,:)=dcc(:,:)
       zalb1(:) = dsalb(1,:)
       zalb2(:) = dsalb(2,:)
       dcc(:,:)=0.
       if(nswr==1) call swr
       dclforc(:,1)=dswfl(:,NLEP)
       dclforc(:,3)=dswfl(:,1)
       dclforc(:,5)=dfu(:,1)
       dclforc(:,6)=dfu(:,NLEP)
       dcc(:,:)=zcc(:,:)
       dsalb(1,:) = zalb1(:)
       dsalb(2,:) = zalb2(:)
       deallocate(zalb1)
       deallocate(zalb2)
      end if

!
!     b) normal computation
!

      if(ntime == 1) call mksecond(zsec1,0.)
      if(nswr==1) call swr
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4swr=time4swr+zsec1
      endif

!
!**   5) long wave radiation
!
!
!     a) if clear sky diagnostic is switched on:
!

      if(ndiagcf > 0) then
       zcc(:,:)=dcc(:,:)
       dcc(:,:)=0.
       if(nlwr==1) call lwr
       dclforc(:,2)=dlwfl(:,NLEP)
       dclforc(:,4)=dlwfl(:,1)
       dclforc(:,7)=dftu(:,NLEP)
       dcc(:,:)=zcc(:,:)
       deallocate(zcc)
      end if

!
!     b) normal computation
!

      if(ntime == 1) call mksecond(zsec1,0.)
      if(nlwr==1) call lwr
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4lwr=time4lwr+zsec1
      endif

!
!**   6) Total flux
!

      dflux(:,:)=dlwfl(:,:)+dswfl(:,:)

!
!**   6a) Get altitudes
!

      do jlev=NLEV,2,-1
       zdh(:,jlev)=-dt(:,jlev)*gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))
      enddo
      zdh(:,1)=-dt(:,1)*gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5
      
!
!**   7) compute tendencies and add them to PUMA dtdt
!

      do jlev = 1 , NLEV
       jlep=jlev+1
       zdtdt(:,jlev)=-ga*(dflux(:,jlep)-dflux(:,jlev))                  &
     &              /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
       dtdt(:,jlev)=dtdt(:,jlev)+zdtdt(:,jlev)
       dtdtswr(:,jlev)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))                &
     &              /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
       dtdtlwr(:,jlev)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))                &
     &              /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
       dconv(:,jlev) = (dflux(:,jlev)-dflux(:,jlep))/(zdh(:,jlev)+1.0e-9) !add small in case of zero zdh
      enddo
      
!
!**   7a) Restore sea ice distribution if applicable
!
      if (nradice==0) dicec(:) = zfice(:)
        

!
!**   8) dbug printout if nprint=2 (see pumamod)
!

      if (nprint==2) then
       allocate(zprf1(NLON*NLAT,NLEP))
       allocate(zprf2(NLON*NLAT,NLEP))
       allocate(zprf3(NLON*NLAT,NLEP))
       allocate(zprf4(NLON*NLAT,NLEP))
       allocate(zprf5(NLON*NLAT,NLEP))
       allocate(zprf6(NLON*NLAT,NLEP))
       allocate(zprf7(NLON*NLAT,NLEP))
       allocate(zprf8(NLON*NLAT,NLEV))
       allocate(zprf9(NHOR,NLEV))
       allocate(zprf10(NHOR,NLEV))
       allocate(zprf11(NLON*NLAT,NLEV))
       allocate(zprf12(NLON*NLAT,NLEV))
       call mpgagp(zprf1,dfd,NLEP)
       call mpgagp(zprf2,dfu,NLEP)
       call mpgagp(zprf3,dswfl,NLEP)
       call mpgagp(zprf4,dftd,NLEP)
       call mpgagp(zprf5,dftu,NLEP)
       call mpgagp(zprf6,dlwfl,NLEP)
       call mpgagp(zprf7,dflux,NLEP)
       call mpgagp(zprf8,zdtdt,NLEV)
       do jlev = 1 , NLEV
        jlep=jlev+1
        zprf9(:,jlev)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))                 &
     &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
        zprf10(:,jlev)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))                &
     &                /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
       enddo
       call mpgagp(zprf11,zprf9,NLEV)
       call mpgagp(zprf12,zprf10,NLEV)
       if(mypid==NROOT) then
        do jlev=1,NLEP
         write(nud,*)'L= ',jlev,' swd= ',zprf1(nprhor,jlev)                  &
     &                    ,' swu= ',zprf2(nprhor,jlev)                  &
     &                    ,' swt= ',zprf3(nprhor,jlev)
         write(nud,*)'L= ',jlev,' lwd= ',zprf4(nprhor,jlev)                  &
     &                    ,' lwu= ',zprf5(nprhor,jlev)                  &
     &                    ,' lwt= ',zprf6(nprhor,jlev)
         write(nud,*)'L= ',jlev,' totalflux= ',zprf7(nprhor,jlev)
        enddo
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' dtdt= ',zprf8(nprhor,jlev)                 &
     &                    ,' dtsw= ',zprf11(nprhor,jlev)                &
     &                    ,' dtlw= ',zprf12(nprhor,jlev)
        enddo
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
       deallocate(zprf12)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       do jlev = 1 , NLEV
        jlep=jlev+1
        dgp3d(:,jlev,5)=-ga*(dlwfl(:,jlep)-dlwfl(:,jlev))               &
     &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
        dgp3d(:,jlev,6)=-ga*(dswfl(:,jlep)-dswfl(:,jlev))               &
     &               /(dsigma(jlev)*dp(:)*acpd*(1.+ADV*dq(:,jlev)))
       enddo
      end if
!
!     energy diagnostics
!
      if(nenergy > 0) then
       allocate(zdtdte(NHOR,NLEV))
       denergy(:,9)=0.
       denergy(:,10)=0.
       denergy(:,17)=0.
       denergy(:,18)=0.
       denergy(:,19)=0.
       denergy(:,20)=0.
       denergy(:,28)=0.
       do jlev=1,NLEV
        jlep=jlev+1
        denergy(:,9)=denergy(:,9)-(dlwfl(:,jlep)-dlwfl(:,jlev))  
        denergy(:,10)=denergy(:,10)-(dswfl(:,jlep)-dswfl(:,jlev))
        denergy(:,17)=denergy(:,17)-(dftd(:,jlep)-dftd(:,jlev)) 
        denergy(:,18)=denergy(:,18)-(dftu(:,jlep)-dftu(:,jlev))
        denergy(:,19)=denergy(:,19)-(dftue1(:,jlep)-dftue1(:,jlev))
        denergy(:,20)=denergy(:,20)-(dftue2(:,jlep)-dftue2(:,jlev)) 
        denergy(:,28)=denergy(:,28)+dt(:,jlev)*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,9)=-(dlwfl(:,jlep)-dlwfl(:,jlev))
         dener3d(:,jlev,10)=-(dswfl(:,jlep)-dswfl(:,jlev))
         dener3d(:,jlev,17)=-(dftd(:,jlep)-dftd(:,jlev))
         dener3d(:,jlev,18)=-(dftu(:,jlep)-dftu(:,jlev))
         dener3d(:,jlev,19)=-(dftue1(:,jlep)-dftue1(:,jlev))
         dener3d(:,jlev,20)=-(dftue2(:,jlep)-dftue2(:,jlev))
         dener3d(:,jlev,28)=dt(:,jlev)*dsigma(jlev)
        endif
       enddo
       deallocate(zdtdte)
      endif

      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4rad=time4rad+zsec
      endif

      return
      end subroutine radstep

!     ==================
!     SUBROUTINE RADSTOP
!     ==================

      subroutine radstop
      use radmod
!
!     finalizes radiation
!     this *sub* is called by PUMA (PUMA-interface)
!
!     for the inclosed parameterizations only a dummy *sub* is needed
!
!     no PUMA *subs* are used
!
!     no PUMA variables are used
!
      if (mypid==NROOT) call put_restart_real('fixedlon',fixedlon)
!     zsolars is a GLOBAL PAIR, not a distributed gridpoint field. mpputgp
!     gathers NHOR per rank into a local z(NUGP,klev) and writes all of it, so
!     this read 510 elements past the end of zsolars(2) on every rank and wrote
!     8190 elements of buffer into the restart. solarini already uses the right
!     idiom for this array. failure-modes.md class 18, CLIM-37.
!
!     THIS WRITE IS A CONFIGURATION FINGERPRINT AND NOT CHECKPOINTED STATE.
!     Nothing reads it back: the get_restart_array block in radini is commented
!     out, because CLIM-38 made solarini unconditional and the pair is rebuilt
!     from the namelist and the spectrum at every start. It is kept because it
!     is the only place a restart, handed on without its run directory, records
!     which two-band split of the stellar constant the orbits were integrated
!     with -- and the pair is live model state, not a diagnostic: swr forms the
!     band-weighted surface albedo from it at nstartemp = 1. A restart whose
!     zsolars disagrees with the configuration it is being resumed under is a
!     different star, and the fingerprint is what makes that visible.
!     exoplasim/scripts/restart_schema.py carries the conversion policy and
!     says the same thing; world-5rq.
      if (mypid == NROOT) call put_restart_array('zsolars',zsolars,2,2,1)
      

      if(mypid == NROOT .and. ntime == 1) then
       write(nud,*)'******************************************'
       write(nud,*)' CPU usage in RADSTEP (ROOT process only):  '
       write(nud,*)'    All routines : ',time4rad,' s'
       write(nud,*)'    Short wave   : ',time4swr,' s'
       write(nud,*)'    Long  wave   : ',time4lwr,' s'
       write(nud,*)'******************************************'
      endif
!
      return
      end subroutine radstop

!     =================
!     SUBROUTINE SOLANG
!     =================

      subroutine solang
      use radmod
!
!     compute cosine of zenit angle including daily cycle
!
!     the following PUMA variables are used:
!
!     PI         : pi=3.14...
!     nstep      : PUMA time step
!     sid(NLPP)  : sines of gaussian latitudes
!     csq(NLPP)  : cosine**2 of gaussian latitudes
!     cola(NLPP) : cosine of latitude
!
!
!**   1) compute day of the year and hour of the day
!
      if (nperpetual > 0) then
         zcday = nperpetual / real(m_days_per_year)
      else
         zcday = mod(nstep,n_steps_per_year) / real(n_steps_per_year) !Actual fractional progression
      endif

      call ntomin(nstep,imin,ihou,iday,imon,iyea)
      
      istp = mod(nstep,int(ntspd*slowdown+0.5))
      imin = (istp * mpstep*ntspd) / int(ntspd*slowdown+0.5)
      ihou = imin / 60
      imin = mod(imin,60)      
      
!
!**   2) compute declination [radians]
!
      if (ngenkeplerian == 0) then
          call orb_decl(zcday, eccen, mvelpp, lambm0, obliqr, orbnu, lambm, rasc, zdecl, eccf)
      else
!           write(6,*) meananom0r
          call gen_orb_decl(zcday, eccen, obliqr, mvelpp, orbnu, lambm, rasc, zdecl, eccf)
      endif
      zcdayf = zcday
      zdeclf = zdecl
!
!**   3) compute zenith angle
!
      gmu0(:) = 0.0
      zmuz    = 0.0
      zdawn = sin(dawn * PI / 180.0) ! compute dawn/dusk angle 
      zrlon = TWOPI / NLON           ! scale lambda to radians
      zrtim = rotspd * TWOPI / 1440.0         ! scale time   to radians
      zmins = ihou * 60 + imin
      
      if (nfixed==1) then
        if (mypid==NROOT) fixedlon = fixedlon + desync*mpstep
        call mpbcr(fixedlon)
        zrtim = TWOPI
        zmins = 1.0 - (fixedlon/360.)  !Think about how to fix this: there's a dep
        zdecl = obliqr                 !on rotspd. Maybe zrtim = TWOPI/1440.0?
      endif
      jhor = 0
      if (ncstsol==0) then
       do jlat = 1 , NLPP
        do jlon = 0 , NLON-1
         jhor = jhor + 1
         zhangle = zmins * zrtim + jlon * zrlon - PI
         if (ngenkeplerian==1) zhangle = zhangle - rasc
         if (zhangle < -PI) zhangle = zhangle + TWOPI
         if (zhangle > PI) zhangle = zhangle - TWOPI
         
         if (nfixed==1) zhangle = zhangle + PI
         
         zmuz=sin(zdecl)*sid(jlat)+cola(jlat)*cos(zdecl)*cos(zhangle)
         if (zmuz > zdawn) gmu0(jhor) = zmuz
        enddo
       enddo
      else
       solclatcdec=solclat*solcdec
       solslat=sqrt(1-solclat*solclat)
       solsdec=sqrt(1-solcdec*solcdec)
       solslatsdec=solslat*solsdec
       do jlat = 1 , NLPP
        do jlon = 0 , NLON-1
         jhor = jhor + 1
         if (ndcycle == 1) then 
          zhangle = zmins * zrtim - PI
          zmuz=solslatsdec+solclatcdec*cos(zhangle)
         else
          zmuz=solslatsdec+solclatcdec/PI
         endif
         if (zmuz > zdawn) gmu0(jhor) = zmuz
        enddo
       enddo
      endif
!
!**  4) copy earth-sun distance (1/r**2) to radmod
!
      gdist2=eccf

      return
      end subroutine solang


!     ===============
!     SUBROUTINE MKO3
!     ===============

      subroutine mko3
      use radmod
      implicit none
!
!     compute ozon distribution
!
!     the following variables from pumamod are used:
!
!     TWOPI    : 2*PI
!     nstep    : PLASIM time step
!     sid(NLAT): double precision sines of gaussian latitudes
!
!     local parameter and arrays
!
      integer :: jlat   ! latitude index
      integer :: jlev   ! level    index
      integer :: jh1
      integer :: jh2
      integer :: imonth ! current month
      integer :: jm     ! index to next or previous month

      real, parameter :: zt0  = 255.
      real, parameter :: zro3 =   2.14
      real, parameter :: zfo3 = 100.0 / zro3

      real :: zcday       ! current day
      real :: zconst
      real :: zw          ! interpolation weight

      real :: za(NHOR)
      real :: zh(NHOR)
      real :: zo3t(NHOR)
      real :: zo3(NHOR)

      if (no3 == 1) then ! compute synthetic ozone distribution
         zcday = mod(nstep,n_steps_per_year) / real(n_steps_per_year)
         do jlat = 1 , NLPP
            jh2 = jlat * NLON     ! horizonatl index for end   of latitude
            jh1 = jh2  - NLON + 1 ! horizontal index for start of latitude
            za(jh1:jh2)=a0o3+a1o3*ABS(sid(jlat))                        &
               +aco3*sid(jlat)*cos(TWOPI*(zcday-toffo3)) 
         enddo ! jlat

         zconst  = exp(-bo3/co3)
         zo3t(:) = za(:)
         zh(:)   = 0.0

         do jlev=NLEV,2,-1
            zh(:)=zh(:)-dt(:,jlev)*GASCON/ga*alog(sigmah(jlev-1)/sigmah(jlev))
            zo3(:)=-(za(:)+za(:)*zconst)/(1.+exp((zh(:)-bo3)/co3))+zo3t(:)
            dqo3(:,jlev)=zo3(:)*ga/(zfo3*dsigma(jlev)*dp(:))
            zo3t(:)=zo3t(:)-zo3(:)
         enddo
         dqo3(:,1) = zo3t(:) * ga / (zfo3 * dsigma(1) * dp(:))
      elseif (no3 == 2) then ! interpolate from climatological ozone
         call momint(nperpetual,nstep+1,imonth,jm,zw)
         dqo3(:,:) = (1.0 - zw) * dqo3cl(:,:,imonth) + zw * dqo3cl(:,:,jm)
      endif ! no3

      if (o3scale /= 1.0) dqo3(:,:) = o3scale * dqo3(:,:)

      return
      end subroutine mko3

!     ===================
!     SUBROUTINE DUSTPROF
!     ===================

      subroutine dustprof
      use radmod
!
!     Spread the PRESCRIBED column dust optical depth over the model layers.
!
!     The offline chain this field comes from carries dust as a well-mixed
!     column of declared scale height, so the vertical shape here is the same
!     one: the concentration decays with `dusthsc` and the layer burden is that
!     concentration times the layer thickness.
!
!     The weights are NORMALISED, so the column optical depth is the prescribed
!     one to roundoff whatever the temperature profile does to the thicknesses.
!     That is the point: it makes the model's own column an IDENTITY against the
!     boundary field rather than a number that has to be believed, and the check
!     below is written to fail if the normalisation ever stops holding.
!
      real :: zdz(NHOR,NLEV)   ! layer thickness (m)
      real :: zzc(NHOR,NLEV)   ! height of the layer centre above ground (m)
      real :: zw(NHOR,NLEV)    ! unnormalised layer dust burden
      real :: zsum(NHOR)       ! column normalisation
      real :: zcol(NHOR)       ! column optical depth, for the check
      real :: zres             ! worst column residual, for the check
      integer :: jaer          ! aerosol species index
!
      ddustod(:,:,:) = 0.
      if (ndustrad < 1) return
!
!     layer thickness in m, built exactly as the shortwave aerosol block builds
!     it, so the two cannot drift apart
!
      do jlev = NLEV,2,-1
       zdz(:,jlev) = -dt(:,jlev)*gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))
      enddo
      zdz(:,1) = -dt(:,1)*gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5
!
!     height of each layer centre above the ground
!
      zzc(:,NLEV) = 0.5*zdz(:,NLEV)
      do jlev = NLEV-1,1,-1
       zzc(:,jlev) = zzc(:,jlev+1)+0.5*(zdz(:,jlev+1)+zdz(:,jlev))
      enddo
!
!     burden per layer, then normalise onto the prescribed column. Per species,
!     because each carries its own scale height: a sea salt layer sits in the
!     boundary layer and dust does not, and one shared height would put the
!     second species at the first one's while the column total still looked
!     right.
!
      do jaer = 1,ndustrad
       zsum(:) = 0.
       do jlev = 1,NLEV
        zw(:,jlev) = EXP(-zzc(:,jlev)/dusthsc(jaer))*zdz(:,jlev)
        zsum(:) = zsum(:)+zw(:,jlev)
       enddo
       do jlev = 1,NLEV
        ddustod(:,jlev,jaer) = dustsc(jaer)*ddustcol(:,jaer)*zw(:,jlev)/MAX(zsum(:),1.E-30)
       enddo
!
!     the identity, reported once PER SPECIES. A nonzero residual means the
!     vertical distribution is not conserving the column and every optical
!     depth below is wrong by that much. Per species and not over the sum, or
!     a species given the wrong scale height hides inside a correct total.
!
       if (ldustchk) then
        zcol(:) = 0.
        do jlev = 1,NLEV
         zcol(:) = zcol(:)+ddustod(:,jlev,jaer)
        enddo
        zres = MAXVAL(ABS(zcol(:)-dustsc(jaer)*ddustcol(:,jaer)))
        if (mypid == NROOT) then
         write(nud,*) 'PRESCRIBED AEROSOL species ',jaer
         write(nud,*) '  max |column - prescribed| = ',zres
         write(nud,*) '  scale height (m)          = ',dusthsc(jaer)
         write(nud,*) '  layer 1 (top) mass share  = ',                     &
     &                MAXVAL(zw(:,1)/MAX(zsum(:),1.E-30))
         write(nud,*) '  layer NLEV mass share     = ',                     &
     &                MAXVAL(zw(:,NLEV)/MAX(zsum(:),1.E-30))
        endif
       endif
      enddo
      ldustchk = .false.
!
      return
      end subroutine dustprof

!     ==================
!     SUBROUTINE AEROPROF
!     ==================

      subroutine aeroprof
      use radmod
!
!     Band-1 extinction optical depth per layer for the INTERACTIVE aerosol,
!     built once per radiation step from the transported number density.
!
!     It exists so that the shortwave and the longwave read ONE field. They
!     used to be unable to disagree only because the longwave had no aerosol
!     term at all; now that it has one, a second copy of this arithmetic in
!     swr would be a copy that can drift, and the layer thicknesses depend on
!     the temperature profile, so it is not a constant that could be built once
!     and kept.
!
!     nrho is floored at one particle per cubic metre, which is what swr did
!     before this and is kept: the two-stream factors below divide by the
!     optical depth.
!
      real :: zdz(NHOR,NLEV)   ! layer thickness (m)
      integer :: jlev
!
      daerod(:,:) = 0.
      if (iaerint /= 1) return
!
      nrho(:,:) = max(1.0,nrho(:,:))
!
      do jlev = NLEV,2,-1
       zdz(:,jlev) = -dt(:,jlev)*gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))
      enddo
      zdz(:,1) = -dt(:,1)*gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5
!
      do jlev = 1,NLEV
       daerod(:,jlev) = nrho(:,jlev)*PI*(apart**2)*qex1(ndustrad+1)*zdz(:,jlev)
      enddo
!
      return
      end subroutine aeroprof

!     ==================
!     SUBROUTINE AEROGATHER
!     ==================

      subroutine aerogather
      use radmod
!
!     Gather every active species' band-1 optical depth per layer into one
!     array, so swr and lwr read ONE field. Prescribed species occupy
!     1..ndustrad and the transported tracer, if it is on, is ndustrad+1.
!
!     This is where the exclusion between the two paths used to live. They no
!     longer share a slot, so there is nothing to refuse: their optical depths
!     sit in different species and are mixed rather than added into one.
!
      integer :: jaer, jlev
!
      aodsp(:,:,:) = 0.
      do jaer = 1,ndustrad
       do jlev = 1,NLEV
        aodsp(:,jlev,jaer) = ddustod(:,jlev,jaer)
       enddo
      enddo
      if (iaerint == 1) then
       do jlev = 1,NLEV
        aodsp(:,jlev,ndustrad+1) = daerod(:,jlev)
       enddo
      endif
!
      return
      end subroutine aerogather

!     ==============
!     SUBROUTINE SWR
!     ==============

      subroutine swr
      use radmod
!
!     calculate short wave radiation fluxes
!
!     this parameterization of sw-radiation bases on transmissivities
!     from Lacis & Hansen (1974) for clear sky (H2O,O3,Rayleigh)
!     and Stephens (1978) + Stephens et al. (1984) for clouds.
!     for the verical integration, the adding method is used.
!     some aspects of the realisation are taken from
!     'a simple radiation parameterization for use in mesoscale models'
!     by S. Bakan (Max-Planck Institut fuer Meteorologie)
!     (unfortunately neither published nor finished)
!
!     no PUMA *subs* are used
!
!     the following PUMA variables are used/modified:
!
!     ga               : gravity acceleration (m/s2) (used)
!     sigma(NLEV)      : full level sigma  (used)
!     sigmah(NLEV)     : half level sigma  (used)
!     dsigma(NLEV)     : delta sigma (half level)  (used)
!     dp(NHOR)         : surface pressure (Pa) (used)
!     dalb(NHOR)       : surface albedo (used)
!     dsalb(2,NHOR)    : band-specific surface albedo (used)
!     dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
!     dql(NHOR,NLEP)   : cloud liquid water content (kg/kg) (used)
!     dcc(NHOR,NLEP)   : cloud cover (frac.) (used)
!     nrho(NHOR,NLEP)  : number density of aerosol (particles/m3) (used if l_aerorad = 1) (used)
!     dswfl(NHOR,NLEP) : short wave radiation (W/m2)  (modified)
!     dfu(NHOR,NLEP)   : short wave radiation upward (W/m2) (modified)
!     dfd(NHOR,NLEP)   : short wave radiation downward (W/m2) (modified)

!
!     0) define local parameters and arrays
!
      parameter(zero=1.E-6)     ! if insolation < zero : fluxes=0.
      parameter(zbetta=1.66)    ! magnification factor water vapour
      parameter(zmbar=1.9)      ! magnification factor ozon
      parameter(zro3=2.14)      ! ozon density (kg/m**3 STP)
      parameter(zfo3=100./zro3) ! transfere o3 to cm STP
!     CO2, the same constants lwr already uses, so the shortwave column and the
!     longwave one are the same quantity. zpv2pm is the molecular weight ratio
!     that turns a VOLUME mixing ratio into a mass one: hydrostatic balance
!     converts total pressure into mass, not partial pressure, and leaving it out
!     understates the column by a factor of 1.52.
      parameter(zmmair=0.0289644)  ! molecular weight air (kg/mol)
      parameter(zmmco2=0.0440098)  ! molecular weight co2 (kg/mol)
      parameter(zpv2pm=zmmco2/zmmair) ! co2 ppvol to ppmass
      parameter(zrco2=1.9635)      ! co2 density (kg/m**3 STP)
      parameter(zfco2=100./zrco2)  ! transfere co2 to cm STP
!     The solar-weighted CO2 absorptance, as a fraction of TOTAL incident flux
!     against the absorber amount u in atmos-cm:
!
!         A(u) = zca1 ln(1 + zcb1 u) + zca2 ln(1 + zcb2 u)
!
!     Two logarithms rather than Lacis & Hansen's Eq. 21 form because that form
!     fits this curve worse and wants a negative coefficient in its denominator,
!     which can go singular on a column nothing here forbids. Quoted over 1 to
!     1e4 atmos-cm, which the model never leaves: the thinnest sigma layer
!     carries a few percent of the column and the smallest magnification is
!     zbetta.
!
!     FITTED TO A LINE LIST, not to Howard's 1956 band set. The level comes from
!     HITRAN2020 through the Generic PCM correlated-k tables, integrated against
!     a 5772 K spectrum and multiplied per band by the water vapour transmission
!     so the two gases do not both claim the same photons, by
!     exoplasim/scripts/corrk_cross_check.py --fit. It replaces a fit to
!     Howard's bands that is 7.2% stronger at this planet's path.
!
!     REFITTED once since, unchanged in method, when world-olt closed the hole
!     the correlated-k bundle's join used to leave between its IR and VI band
!     sets: 1974.95 to 2000 cm-1 was inside the table span and inside no band,
!     and every band-weighted total priced it as transparent. Closing it adds
!     1.7e-4 of the flux at 5 um and moves the absorptance at this planet's path
!     by +0.5%, worth +0.01 W/m2 of insolation -- an eighth of the fit's own rms
!     residual, so the coefficients follow their derivation rather than the
!     result being worth chasing.
!
!     WHY THE WHOLE CURVE AND NOT THE TWO BANDS THAT WERE WRONG. The defect
!     found was in the per-band water overlap: Howard's band-mean absorptance,
!     spread uniformly across a 1000 cm-1 interval, smears saturation out of the
!     band cores, putting the 2.7 um clear fraction at 0.171 where the line list
!     says 0.003 and the 2.0 um one at 0.479 where it says 0.718. Correcting
!     only those makes the TOTAL worse: summed over Howard's eight intervals the
!     line list gives 0.004520 against the derivation's 0.005536, 18.4% low,
!     while over all bands it gives 0.005123, 7.5% low. About an eighth of the
!     CO2 shortwave absorption falls outside every interval Howard measured, and
!     the overlap error was partly standing in for it. Taking the level from the
!     line list across the range needs no band bookkeeping and cannot cancel two
!     errors by accident.
!
!     The form is a worse fit to this curve than to Howard's, and the cost sits
!     where the model does not go: 8.8% of itself at u = 1e4, but 4.5% over
!     100 to 1000 atmos-cm where a T42 column sits, against 4.0% for the fit it
!     replaces. co2sww is unaffected -- it is a star-over-Sun RATIO and the two
!     derivations agree on it to 0.03%.
      parameter(zca1=3.1020E-4)
      parameter(zcb1=19.857)
      parameter(zca2=3.8291E-3)
      parameter(zcb2=3.9587E-3)
      parameter(aa=0.2542857142857143)
      parameter(bb=0.8229693877551021)
      parameter(c0=0.14997959183673468)
!
      real zt1(NHOR,NLEP),zt2(NHOR,NLEP)    ! transmissivities 1-l
      real zr1s(NHOR,NLEP),zr2s(NHOR,NLEP)  ! reflexivities l-1 (scattered)
      real zrl1(NHOR,NLEP),zrl2(NHOR,NLEP)  ! reflexivities l-NL (direct)
      real zrl1s(NHOR,NLEP),zrl2s(NHOR,NLEP)! reflexivities l-NL (scattered)
!
      real ztb1(NHOR,NLEV),ztb2(NHOR,NLEV)   ! layer transmissivity (down)
      real ztb1u(NHOR,NLEV),ztb2u(NHOR,NLEV) ! layer transmissivity (up)
      real zrb1(NHOR,NLEV),zrb2(NHOR,NLEV)   ! layer reflexivity (direct)
      real zrb1s(NHOR,NLEV),zrb2s(NHOR,NLEV) ! layer reflexibity (scattered)
!
      real zo3l(NHOR,NLEV)   ! ozon amount (top-l)
      real zxo3l(NHOR,NLEV)  ! effective ozon amount (top-l)
      real zwvl(NHOR,NLEV)   ! water vapor amount (top-l)
      real zywvl(NHOR,NLEV)  ! effective water vapor amount (top-l)
      real zco2l(NHOR,NLEV)  ! co2 amount (top-l)
      real zyco2l(NHOR,NLEV) ! effective co2 amount (top-l)
      real zrcs(NHOR,NLEV)   ! clear sky reflexivity (downward beam)
      real zrcsu(NHOR,NLEV)  ! clear sky reflexivity (upward beam)
!
      real zftop1(NHOR),zftop2(NHOR) ! top solar radiation
      real gsolinst,zcyclephase,zcyclephase2 ! instantaneous cyclic stellar flux
      real zfu1(NHOR),zfu2(NHOR)     ! upward fluxes
      real zfd1(NHOR),zfd2(NHOR)     ! downward fluxes
!
      real zmu0(NHOR)              ! zenit angle
      real zmu1(NHOR)              ! zenit angle
      real zmu00                   ! Diffusivity factor for scattered/diffuse light
      real zcs(NHOR)               ! clear sky part
      real zscf(NHOR)              ! Pressure scale factor
      real zm(NHOR)                ! magnification factor
      real zo3(NHOR)               ! ozon amount
      real zo3t(NHOR)              ! total ozon amount (top-sfc)
      real zxo3t(NHOR)             ! effective total ozon amount (top-sfc)
      real zto3(NHOR),zto3u(NHOR)  ! ozon transmissivity (downward/upward beam)
      real zto3t(NHOR),zto3tu(NHOR)! total ozon transmissivities (d/u)
      real zwv(NHOR)               ! water vapor amount
      real zwvt(NHOR)              ! total water vapor amount (top-sfc)
      real zywvt(NHOR)             ! total effective water vapor amount (top-sfc)
      real ztwv(NHOR),ztwvu(NHOR)  ! water vapor trasmissivity (d/u)
      real ztwvt(NHOR),ztwvtu(NHOR)! total water vapor transmissivities (d/u)
      real zco2(NHOR)              ! co2 amount
      real zco2t(NHOR)             ! total co2 amount (top-sfc)
      real zyco2t(NHOR)            ! total effective co2 amount (top-sfc)
      real ztco2(NHOR),ztco2u(NHOR)! co2 transmissivity (d/u)
      real ztco2t(NHOR),ztco2tu(NHOR)! total co2 transmissivities (d/u)
      real zzco2                   ! co2 mass mixing ratio (kg/kg)
!
      real zra1(NHOR),zra2(NHOR)   ! reflexivities combined layer (direct)
      real zra1s(NHOR),zra2s(NHOR) ! reflexivities combined layer (scatterd)
      real zta1(NHOR),zta2(NHOR)   ! transmissivities combined layer (di)
      real zta1s(NHOR),zta2s(NHOR) ! transmissivities combined layer (sc)
      real z1mrabr(NHOR)           ! 1/(1.-rb*ra(*))
!
      real zrcl1(NHOR,NLEV),zrcl2(NHOR,NLEV)  ! cloud reflexivities (direct)
      real zrcl1s(NHOR,NLEV),zrcl2s(NHOR,NLEV)! cloud reflexivities (scattered)
      real ztcl2(NHOR,NLEV),ztcl2s(NHOR,NLEV) ! cloud transmissivities
      
      real zaert1(NHOR,NLEV),zaert2(NHOR,NLEV) ! aerosol transmissivities (direct)
      real zaerr1(NHOR,NLEV),zaerr2(NHOR,NLEV) ! aerosol reflectivities (direct)
      real zaert1s(NHOR,NLEV),zaert2s(NHOR,NLEV) ! aerosol transmissivities (scattered)
      real zaerr1s(NHOR,NLEV),zaerr2s(NHOR,NLEV) ! aerosol reflectivities (scattered)

    ! Local intermediate arrays for aerosol calculations
!
!     The u-factors were scalars because the mixture's optical properties were.
!     They are now per cell, computed inside the level loop, because an
!     external mixture's single-scattering albedo and backscatter ratio depend
!     on which species are present in THAT layer and in what proportion.
!     CLIM-39, aeolian/notes/multi-species-aerosol.md section 2.
!
      real :: zaeru1(NHOR),zaeru2(NHOR) ! U-factors, per cell
      real :: ztemp1(NHOR),ztemp2(NHOR) ! per cell
      real :: zssa1(NHOR),zssa2(NHOR)   ! mixture single-scattering albedo
      real :: zbs1(NHOR),zbs2(NHOR)     ! mixture backscatter ratio
      real :: zext1(NHOR),zext2(NHOR)   ! sum of extinction optical depth
      real :: zsca1(NHOR),zsca2(NHOR)   ! sum of scattering optical depth
      real :: zbsc1(NHOR),zbsc2(NHOR)   ! sum of backscattered optical depth
      logical :: lcons1(NHOR),lcons2(NHOR) ! conservative-scattering cells
      real :: ztcon(NHOR)               ! b*tau/mu in the conservative limit
      real :: zepsc                     ! conservative-scattering threshold
      integer :: knz(NHOR)              ! species contributing in this cell
      integer :: klast(NHOR)            ! index of the last one that did
      integer :: jaer                   ! aerosol species index
      real :: zaertf1(NHOR,NLEV),zaertf2(NHOR,NLEV) ! Effective optical depth (direct light)
      real :: zaertf1s(NHOR,NLEV),zaertf2s(NHOR,NLEV) ! (scattered light)
      real :: zaerd1(NHOR,NLEV),zaerd2(NHOR,NLEV) ! Denominator (direct light)
      real :: zaerd1s(NHOR,NLEV),zaerd2s(NHOR,NLEV) ! (scattered light)
      real :: aod1(NHOR,NLEV),aod2(NHOR,NLEV) ! Aerosol optical depth
      integer :: iaeron ! 1 where an aerosol acts on the shortwave, else 0
!
!     arrays for diagnostic cloud properties
!
      real :: zlwp(NHOR)
      real :: zwl(NHOR)
      real :: ztau1(NHOR)
      real :: ztau2(NHOR)
      real :: zlog(NHOR)
      real zb2(NHOR),zom0(NHOR),zuz(NHOR),zun(NHOR),zr(NHOR)
      real zexp(NHOR),zu(NHOR),zb1(NHOR)
!
!     ONE BROADBAND OPTICAL DEPTH PER SHORTWAVE BAND, because Stephens (1978)
!     p. 2124 fits one per band and radmod splits its bands where he splits
!     his. Both fits are least squares to Mie calculations on the eight
!     standard cloud models of Part 1, in the form
!
!       log10(tau_N) = a + b ln(log10 W)      W the liquid water path, g m-2
!
!     with (a,b) = (0.2633, 1.7095) over 0.3 to 0.75 um, Eq. (10a), and
!     (0.3492, 1.6518) over 0.75 to 4.0 um, Eq. (10b). That is a power law in
!     log10 W: 10**(a + b ln x) is 10**a times x**(b ln 10), so the pairs below
!     are (10**a, b*ln 10) for each band. Eq. (10b) lies above Eq. (10a) at
!     every path, so a single optical depth serving both bands makes the
!     visible band's cloud too bright. world-jgen.
!
      real, parameter :: ztaua1 = 1.8336  ! 10**0.2633,  Eq. (10a) prefactor
      real, parameter :: ztaup1 = 3.9363  ! 1.7095*ln10, Eq. (10a) exponent
      real, parameter :: ztaua2 = 2.2346  ! 10**0.3492,  Eq. (10b) prefactor
      real, parameter :: ztaup2 = 3.8034  ! 1.6518*ln10, Eq. (10b) exponent
!
!     THE FIT'S FLOOR IS 10 g m-2 AND THIS MODEL'S HIGH CLOUD LIVES BELOW IT.
!     Figs. 1a and 1b span 10 to 10000 g m-2 and the fitted form is singular at
!     W = 1 and undefined below it, so below 10 g m-2 every value is an
!     extrapolation and the choice of continuation is a physical statement
!     rather than a guard against a domain error. Extending the fit itself is
!     the wrong statement: its elasticity d ln tau / d ln W is
!     b*ln10 / (ln10 * log10 W), which is already 1.7 at 10 g m-2 and grows
!     without bound as W falls towards 1, so the extrapolated optical depth
!     collapses far faster than linearly. Read back through Stephens Eq. (7),
!     tau_N = 1.5 W / r_e, that is a droplet effective radius growing past
!     100 um as the cloud thins, which no cloud does.
!
!     The continuation used below is Eq. (7) at fixed r_e: tau is linear in W
!     under ZWFIT and matches the fit's own value at ZWFIT, so the effective
!     radius is held at what the fit implies at the bottom of its range --
!     8.2 um in band 1, 6.7 um in band 2 -- rather than allowed to run away. It
!     is continuous at ZWFIT, it introduces no constant the fit does not
!     already contain, and it returns zero optical depth at zero cloud water,
!     which the inherited `1.5 +` offset did not: that offset left 0.055 of
!     optical depth in a layer holding no cloud water at all. Stephens,
!     Ackerman and Smith (1984) p. 690 state that reflection below an optical
!     depth of about 2 needs a parameterization of its own; this is the
!     weakest continuation that does not assert something false about the
!     droplets. world-jgen.
!
      real, parameter :: zwfit = 10.0     ! g m-2, bottom of the fitted range
!
      logical losun(NHOR)         ! flag for gridpoints with insolation
!
!     cosine of zenith angle
!
      if (ndcycle == 0) then ! compute zonal means
         js = 1
         je = NLON
         do jlat = 1 , NLPP
            icnt = count(gmu0(js:je) > 0.0)
            if (icnt > 0) then
               zsum = sum(gmu0(js:je))
               zmu0(js:je) = zsum / icnt ! used for clouds
               zmu1(js:je) = zsum / NLON ! used for insolation
            else
               zmu0(js:je) = 0.0
               zmu1(js:je) = 0.0
            endif
            js = js + NLON
            je = je + NLON
         enddo ! jlat
      else
         zmu0(:) = gmu0(:)
         zmu1(:) = gmu0(:)
      endif ! (ndcycle == 0)
!
!     top solar radiation downward
!
      gsolinst = gsol0
      if (nsolcycle > 0 .and. gsolperiod > 0.0) then
         zcyclephase = TWOPI * (real(nstep-gsolstart) / gsolperiod + gsolphase)
         gsolinst = gsol0 + gsolamp * sin(zcyclephase)
         if (gsolperiod2 > 0.0 .and. gsolamp2 /= 0.0) then
            zcyclephase2 = TWOPI * (real(nstep-gsolstart) / gsolperiod2   &
     &                              + gsolphase2)
            gsolinst = gsolinst + gsolamp2 * sin(zcyclephase2)
         endif
      endif
      zftop1(:) = zsolar1 * gsolinst * gdist2 * zmu1(:) !Adjust down here for redder spectrum. --AYP
      zftop2(:) = zsolar2 * gsolinst * gdist2 * zmu1(:)

!     from this point on, all computations are made only for
!     points with solar insolation > zero
!
      losun(:) = (zftop1(:) + zftop2(:) > zero)
!
!     cloud properites
!
      zcs(:) = 1.0 ! Clear sky fraction (1.0 = clear sky)
      zmu00  = 0.5
      zb3    = tswr1 * SQRT(zmu00) / zmu00
      zb4    = tswr2 * SQRT(zmu00)
      zb5    = tswr3 * zmu00 * zmu00
!
!     prescribed
!
!     THE CLOUD OPTICS ARE PRESET TO THE CLEAR-SKY VALUES, above the branch and
!     not inside it, because at NCLOUDS = 0 nothing below writes them at all and
!     every read of them further down is guarded by a multiplication rather than
!     by a branch: `zrcl1*dcc*nclouds` in the range-1 reflectivity, `zrcl2` and
!     `zrcl2s` in the range-2 pair, `(1.-ztcl2)*dcc*nclouds` in the
!     transmissivity. Zero times a stale non-finite word is not zero, it is the
!     invalid the declared flag line traps on. Reflectivity 0 and transmissivity
!     1 are what no cloud means, so the preset is also the right answer and not
!     merely a defined one. Bit-identical at NCLOUDS = 1: the NSWRCL = 0 arm
!     writes every element of all six, and the NSWRCL = 1 arm made this same
!     preset itself before its masked writes. world-35en, and the same class as
!     world-5a0 and world-bhs below.
!
      zrcl1(:,:)=0.0
      zrcl2(:,:)=0.0
      ztcl2(:,:)=1.0
      zrcl1s(:,:)=0.0
      zrcl2s(:,:)=0.0
      ztcl2s(:,:)=1.0
!
      if (nclouds==1) then
      if (nswrcl == 0) then
       do jlev=1,NLEV
        if(sigma(jlev) <= 1./3.) then
         zrcl1s(:,jlev)=rcl1(1)/(rcl1(1)+zmu00)
         zrcl1(:,jlev)=zcs(:)*rcl1(1)/(rcl1(1)+zmu0(:))                 &
     &                +(1.-zcs(:))*zrcl1s(:,jlev)
         zrcl2s(:,jlev)=AMIN1(1.-acl2(1),rcl2(1)/(rcl2(1)+zmu00))
         zrcl2(:,jlev)=AMIN1(1.-acl2(1),zcs(:)*rcl2(1)/(rcl2(1)+zmu0(:))&
     &                                +(1.-zcs(:))*zrcl2s(:,jlev))
         ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(1)
         ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(1)
        elseif(sigma(jlev) > 1./3. .and. sigma(jlev) <= 2./3.) then
         zrcl1s(:,jlev)=rcl1(2)/(rcl1(2)+zmu00)
         zrcl1(:,jlev)=zcs(:)*rcl1(2)/(rcl1(2)+zmu0(:))                 &
     &                +(1.-zcs(:))*zrcl1s(:,jlev)
         zrcl2s(:,jlev)=AMIN1(1.-acl2(2),rcl2(2)/(rcl2(2)+zmu00))
         zrcl2(:,jlev)=AMIN1(1.-acl2(2),zcs(:)*rcl2(2)/(rcl2(2)+zmu0(:))&
     &                                 +(1.-zcs(:))*zrcl2s(:,jlev))
         ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(2)
         ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(2)
        else
         zrcl1s(:,jlev)=rcl1(3)/(rcl1(3)+zmu00)
         zrcl1(:,jlev)=zcs(:)*rcl1(3)/(rcl1(3)+zmu0(:))                 &
     &                +(1.-zcs(:))*zrcl1s(:,jlev)
         zrcl2s(:,jlev)=AMIN1(1.-acl2(3),rcl2(3)/(rcl2(3)+zmu00))
         zrcl2(:,jlev)=AMIN1(1.-acl2(3),zcs(:)*rcl2(3)/(rcl2(3)+zmu0(:))&
     &                                 +(1.-zcs(:))*zrcl2s(:,jlev))
         ztcl2s(:,jlev)=1.-zrcl2s(:,jlev)-acl2(3)
         ztcl2(:,jlev)=1.-zrcl2(:,jlev)-acl2(3)
        endif
        zcs(:)=zcs(:)*(1.-dcc(:,jlev))
       enddo
      else
!
!     THE CHAIN'S SCRATCH IS PRESET, for the reason the magnification factor
!     further down is set unconditionally. Every one of these is written inside
!     the where below and read by the next line of the same block, and a `where`
!     masks the ASSIGNMENT rather than the evaluation: a lane with no cloud at
!     this level is never stored to, so the read takes whatever the stack held.
!     A lane the mask KEEPS is stored before every read of it, so presetting
!     cannot move a result the model uses. ztau2 is preset to 1 and not to 0
!     because 1000/ztau2 is formed on every lane, and zwl to 1 because it is
!     raised to a non-integer power on every lane. world-5a0.
!
       zlwp(:) = 0.0
       zwl(:)  = 1.0
       ztau1(:)= 0.0
       ztau2(:)= 1.0
       zlog(:) = 0.0
       zb1(:)  = 0.0
       zb2(:)  = 0.0
       zom0(:) = 0.0
       zun(:)  = 1.0
       zuz(:)  = 1.0
       zu(:)   = 1.0
       zexp(:) = 1.0
       zr(:)   = 1.0
       do jlev=1,NLEV
        where(losun(:) .and. (dcc(:,jlev) > 0.))
!
!     THE LIQUID WATER PATH IS FLOORED AT ZERO INSIDE THE LOGARITHM, and that
!     floor is what makes the rest of this chain safe. dql is set to 0 and then
!     to MAX(dql,1.E-9) in rainmod, so zlwp cannot be negative on any lane the
!     model has computed and the floor cannot bind on one. What it removes is a
!     stale lane, whose word can be anything at all. world-5a0.
!
!     THE TWO BRANCHES OF tau(W) ARE ONE EXPRESSION, and it is safe for every
!     zlwp from zero up. Above ZWFIT, min(1,zlwp/ZWFIT) is 1 and zwl is
!     log10(zlwp), so each band evaluates its own fit. Below ZWFIT, zwl is
!     log10(ZWFIT) = 1 exactly, so zwl**p is 1 and what is left is
!     ztaua_b * zlwp/ZWFIT: linear in the water path, matching the fit at
!     ZWFIT. The base of the power is never below 1 and never negative, so
!     there is no domain error to guard, and no offset is added to the water
!     path. world-jgen.
!
         zlwp(:) = min(1000.0,1000.*dql(:,jlev)*dp(:)/ga*dsigma(jlev))
         zwl(:)  = ALOG10(max(zwfit,max(0.,zlwp(:))))
         ztau1(:)= ztaua1*zwl(:)**ztaup1*min(1.0,max(0.,zlwp(:))/zwfit)
         ztau2(:)= ztaua2*zwl(:)**ztaup2*min(1.0,max(0.,zlwp(:))/zwfit)
!
!     ztau2 IS FLOORED ONLY WHERE IT IS DIVIDED INTO, and the floor is inert.
!     The optical depth now reaches zero on a layer holding no cloud water, so
!     1000/ztau2 needs a divisor that cannot be zero. At the floor the layer is
!     already transparent -- zexp is 1, so the band-2 reflectivity carries the
!     factor (zexp - 1/zexp) = 0 and the transmissivity is 4u/4u = 1 -- and the
!     only thing zlog still reaches is zom0, which multiplies nothing that
!     survives. So the value of the floor cannot move a result; it only keeps
!     the divide defined. world-jgen.
!
         zlog(:) = log(max(1.E-30,1000.0 / max(1.E-10,ztau2(:))))
         zb2(:)  = zb4 / ALOG(3.+0.1*ztau2(:))
         zom0(:) = min(0.9999,1.0 - zb5 * zlog(:))
         zun(:)  = 1.0 - zom0(:)
         zuz(:)  = zun(:) + 2.0 * zb2(:) * zom0(:)
         zu(:)   = SQRT(max(0.,zuz(:)/zun(:)))
         zexp(:) = exp(min(25.0,ztau2(:)*SQRT(max(0.,zuz(:)*zun(:)))/zmu00))
         zr(:)   = (zu(:)+1.)*(zu(:)+1.)*zexp(:)                      &
     &           - (zu(:)-1.)*(zu(:)-1.)/zexp(:)
         zrcl1s(:,jlev)=1.-1./(1.+zb3*ztau1(:))
         ztcl2s(:,jlev)=4.*zu(:)/zr(:)
         zrcl2s(:,jlev)=(zu(:)*zu(:)-1.)/zr(:)*(zexp(:)-1./zexp(:))

         zb1(:)  = tswr1*SQRT(max(0.,zmu0(:)))
         zb2(:)  = tswr2*SQRT(max(0.,zmu0(:)))/ALOG(3.+0.1*ztau2(:))
         zom0(:) = min(0.9999,1.-tswr3*zmu0(:)*zmu0(:)*zlog(:))
         zun(:)  = 1.0 - zom0(:)
         zuz(:)  = zun(:) + 2.0 * zb2(:) * zom0(:)
         zu(:)   = SQRT(max(0.,zuz(:)/zun(:)))
         zexp(:) = exp(min(25.0,ztau2(:)*SQRT(max(0.,zuz(:)*zun(:)))/max(1.E-30,zmu0(:))))
         zr(:)   = (zu(:)+1.)*(zu(:)+1.)*zexp(:)                      &
     &           - (zu(:)-1.)*(zu(:)-1.)/zexp(:)
         zrcl1(:,jlev)=1.-1./(1.+zb1(:)*ztau1(:)/max(1.E-30,zmu0(:)))
         ztcl2(:,jlev)=4.*zu(:)/zr(:)
         zrcl2(:,jlev)=(zu(:)*zu(:)-1.)/zr(:)*(zexp(:)-1./zexp(:))
         zrcl1(:,jlev)=zcs(:)*zrcl1(:,jlev)+(1.-zcs(:))*zrcl1s(:,jlev)
         ztcl2(:,jlev)=zcs(:)*ztcl2(:,jlev)+(1.-zcs(:))*ztcl2s(:,jlev)
         zrcl2(:,jlev)=zcs(:)*zrcl2(:,jlev)+(1.-zcs(:))*zrcl2s(:,jlev)
        endwhere
        zcs(:)=zcs(:)*(1.-dcc(:,jlev))
       enddo ! jlev
      endif ! (nswrcl == 0)
      endif ! (nclouds == 1)
!
!     magnification factor
!
!     UNCONDITIONAL, not under where(losun). These are the running column
!     integrals the loop below accumulates into, and setting them only in
!     daylight leaves every night lane holding whatever the stack already had.
!     The loop then reads those lanes -- zo3t(:)=zo3t(:)+... is a read -- so a
!     stale non-finite word propagates out of a lane the mask was supposed to
!     have excluded, and the FPE trap catches it inside the radiation. Setting
!     them for every lane changes nothing where the mask keeps the result, and
!     zm's own formula is safe for ANY zmu0: the square root's argument is
!     1 + 1224*zmu0^2, which cannot be negative. world-bhs.
      zm(:)=35./SQRT(1.+1224.*zmu0(:)*zmu0(:))
!
!     absorber amount and clear sky fraction
!
      zcs(:)=1.
      zo3t(:)=0.
      zxo3t(:)=0.
      zwvt(:)=0.
      zywvt(:)=0.
      zco2t(:)=0.
      zyco2t(:)=0.
!
!     THE WHOLE TWO-STREAM IS PRESET ON EVERY LANE, for the reason the
!     magnification factor above is set unconditionally and the cloud-optics and
!     aerosol chains already are. Every array below has ALL of its definitions
!     inside a where(losun(:)) block and is read inside one, and a `where` masks
!     the ASSIGNMENT rather than the evaluation: on a night lane the store is
!     discarded but the right-hand side may still be computed, so the read is of
!     whatever the stack held. That reaches the transmissivity quotients
!     (1./zto3t, 1./ztwvu, 1./ztco2t and their five siblings), the adding
!     method's 1./(1.-r*r) at four sites, and EXP and LOG of stale words -- all
!     under -ffpe-trap=invalid,zero,overflow. world-px61, and the same class as
!     world-5a0 and world-bhs.
!
!     THE VALUES ARE THE TRANSPARENT ATMOSPHERE OVER A BLACK SURFACE, which is
!     what the arithmetic here produces on a lane carrying no absorber, no
!     scatterer and no insolation: zero absorber amount, unit transmissivity,
!     zero reflectivity, zero flux. Every quotient then divides by exactly 1 --
!     zto3t and its siblings are 1, and 1.-r*r is 1.-0.*0. -- and the two
!     absorptance denominators at the upward beam are 1.-A(0)/zsolar, which is
!     also exactly 1, because zo3t, zxo3t, zwvt, zywvt, zco2t and zyco2t are
!     already zero on every lane by the block above and the per-level copies
!     below are now zero too.
!
!     NOTHING STORED MOVES. On every lane losun keeps, each of these is assigned
!     inside the same masked block, on the same pass, before any read of it: the
!     four combined-layer pairs at the preset above the downward loop, the
!     per-level R and T at the head of each loop body, the column copies at the
!     absorber loop, and z1mrabr, zscf, zfd* and zfu* one statement before their
!     use. The preset is visible only on lanes whose stores are discarded.
!
      zo3(:)     = 0.
      zwv(:)     = 0.
      zco2(:)    = 0.
      zo3l(:,:)  = 0.
      zxo3l(:,:) = 0.
      zwvl(:,:)  = 0.
      zywvl(:,:) = 0.
      zco2l(:,:) = 0.
      zyco2l(:,:)= 0.
!
      zto3(:)    = 1.
      zto3u(:)   = 1.
      zto3t(:)   = 1.
      zto3tu(:)  = 1.
      ztwv(:)    = 1.
      ztwvu(:)   = 1.
      ztwvt(:)   = 1.
      ztwvtu(:)  = 1.
      ztco2(:)   = 1.
      ztco2u(:)  = 1.
      ztco2t(:)  = 1.
      ztco2tu(:) = 1.
!
      zt1(:,:)   = 1.
      zt2(:,:)   = 1.
      ztb1(:,:)  = 1.
      ztb2(:,:)  = 1.
      ztb1u(:,:) = 1.
      ztb2u(:,:) = 1.
      zta1(:)    = 1.
      zta2(:)    = 1.
      zta1s(:)   = 1.
      zta2s(:)   = 1.
!
      zr1s(:,:)  = 0.
      zr2s(:,:)  = 0.
      zrl1(:,:)  = 0.
      zrl2(:,:)  = 0.
      zrl1s(:,:) = 0.
      zrl2s(:,:) = 0.
      zrb1(:,:)  = 0.
      zrb2(:,:)  = 0.
      zrb1s(:,:) = 0.
      zrb2s(:,:) = 0.
      zrcs(:,:)  = 0.
      zrcsu(:,:) = 0.
      zra1(:)    = 0.
      zra2(:)    = 0.
      zra1s(:)   = 0.
      zra2s(:)   = 0.
!
      z1mrabr(:) = 1.
      zscf(:)    = 0.
      zfd1(:)    = 0.
      zfd2(:)    = 0.
      zfu1(:)    = 0.
      zfu2(:)    = 0.
!
!     CO2 is well mixed, so its mass mixing ratio is one scalar for the column.
!     co2 is the namelist volume mixing ratio in ppmv, the same quantity lwr
!     reads through dqco2.
!
      zzco2=zpv2pm*1.E-6*co2
      do jlev=1,NLEV
       where(losun(:))
        zo3(:)=zfo3*dsigma(jlev)*dp(:)*dqo3(:,jlev)/ga
        zo3t(:)=zo3t(:)+zo3(:)
        zxo3t(:)=zcs(:)*(zxo3t(:)+zm(:)*zo3(:))                         &
     &          +(1.-zcs(:))*(zxo3t(:)+zmbar*zo3(:))
        zo3l(:,jlev)=zo3t(:)
        zxo3l(:,jlev)=zxo3t(:)
        zwv(:)=0.1*dsigma(jlev)*dq(:,jlev)*dp(:)/ga                     &
     &        *SQRT(273./dt(:,jlev))*sigma(jlev)*dp(:)/100000.
        zwvt(:)=zwvt(:)+zwv(:)
        zywvt(:)=zcs(:)*(zywvt(:)+zm(:)*zwv(:))                         &
     &          +(1.-zcs(:))*(zywvt(:)+zbetta*zwv(:))
        zwvl(:,jlev)=zwvt(:)
        zywvl(:,jlev)=zywvt(:)
!
!     CO2 amount, reduced to standard pressure the same way the water vapour
!     amount above it is and lwr's own CO2 amount is (radmod.f90 zqco2): the fit
!     is stated at standard pressure and the scheme reaches it by scaling the
!     amount by sigma*ps/p0 rather than by scaling the pressure. No temperature
!     factor, because Howard's CO2 constants carry none and lwr applies none.
!
        zco2(:)=zfco2*zzco2*dsigma(jlev)*dp(:)/ga                        &
     &         *sigma(jlev)*dp(:)/100000.
        zco2t(:)=zco2t(:)+zco2(:)
        zyco2t(:)=zcs(:)*(zyco2t(:)+zm(:)*zco2(:))                       &
     &           +(1.-zcs(:))*(zyco2t(:)+zbetta*zco2(:))
        zco2l(:,jlev)=zco2t(:)
        zyco2l(:,jlev)=zyco2t(:)
        zcs(:)=zcs(:)*(1.-dcc(:,jlev)*nclouds)
        zrcs(:,jlev) = (aa/(1.+bb*zmu0(:))*zcs+c0*(1.-zcs(:)-dcc(:,NLEV)*nclouds))   &
     &                  *dsigma(jlev)*(dp(:)/101100.0)*(9.80665/ga)*nrscat*newrsc
        zrcsu(:,jlev)= c0*(1.-dcc(:,NLEV))*dsigma(jlev)*(dp(:)/101100.0)*(9.80665/ga)*nrscat*newrsc
        
       endwhere
      end do
      
!     aerosol blockaerosol block
!     Initialise aerosol transmissivity and reflectivity arrays to 1.0 and 0.0, respectively
!     so the aerosol block has no effect if l_aerorad==0.

      zaert1(:,:) = 1.0
      zaert2(:,:) = 1.0 
      zaerr1(:,:) = 0.0
      zaerr2(:,:) = 0.0 
      zaert1s(:,:) = 1.0
      zaert2s(:,:) = 1.0
      zaerr1s(:,:) = 0.0
      zaerr2s(:,:) = 0.0

      iaeron = 0
      if (naerosp > 0) iaeron = 1

      if (iaeron == 1) then

!     THE INTERMEDIATES ARE PRESET TOO, for the reason the magnification factor
!     above is set unconditionally. Every one of them is written inside a
!     where(aod1 > 0.) and read inside the same block, and a `where` masks the
!     ASSIGNMENT and not the evaluation: a lane with no aerosol at this level is
!     never stored to, so the read takes whatever the stack held, and EXP of a
!     stale non-finite word raises the invalid the production profile traps on.
!     A lane the mask KEEPS is stored before it is read at every one of these,
!     so presetting them cannot move a result the model uses. The values are the
!     no-aerosol ones the conservative branch below already writes: zero
!     effective optical depth, unit denominator, unit u-factor. world-5a0.
!
!     These stay inside the iaeron test, unlike the transmissivities above.
!     zaert1 and its siblings are read at every level whether or not an aerosol
!     acts; these are read only from inside this block.

      zaertf1(:,:) = 0.0
      zaertf2(:,:) = 0.0
      zaertf1s(:,:) = 0.0
      zaertf2s(:,:) = 0.0
      zaerd1(:,:) = 1.0
      zaerd2(:,:) = 1.0
      zaerd1s(:,:) = 1.0
      zaerd2s(:,:) = 1.0
      zaeru1(:) = 1.0
      zaeru2(:) = 1.0
      ztemp1(:) = 0.0
      ztemp2(:) = 0.0
      ztcon(:) = 0.0

      ! Aerosol two-stream multiscattering radiative transfer approximation from
      ! Stephens (1978), with data read from outside the model instead of a parameterization
      ! for the effective optical depht and single scattering albedo
      !
      ! The species live in aodsp, gathered by aerogather from the prescribed
      ! columns dustprof built and the transported tracer aeroprof built. Band
      ! 2 follows from each species' OWN ratio of extinction efficiencies, so
      ! the band split stays with the optics rather than being restated here.
      !
      ! The threshold below is where the exact expression stops being the more
      ! accurate one. Its error goes as machine epsilon over (1-ssa), because
      ! that difference is formed by subtraction; the conservative limit's
      ! error goes as (1-ssa) itself. The two cross at SQRT(epsilon), which in
      ! 8-byte arithmetic is about 1.5e-8 and in 4-byte about 3.4e-4. Sea salt
      ! at 1-1e-5 therefore takes the exact branch in a double build and the
      ! conservative one in a single build, which is the right answer in both.

       zepsc = SQRT(EPSILON(1.0))

       do jlev=1,NLEV

      ! Mix the species present in this layer. tau adds; the single-scattering
      ! albedo is the extinction-weighted mean and the backscatter ratio the
      ! scattering-weighted mean, per band, because band 2's optical depth is a
      ! different split of the same species.

        zext1(:) = 0.
        zsca1(:) = 0.
        zbsc1(:) = 0.
        zext2(:) = 0.
        zsca2(:) = 0.
        zbsc2(:) = 0.
        knz(:)   = 0
        klast(:) = 1
        do jaer = 1,naerosp
         where (aodsp(:,jlev,jaer) > 0.)
          zext1(:) = zext1(:) + aodsp(:,jlev,jaer)
          zsca1(:) = zsca1(:) + ssa1(jaer)*aodsp(:,jlev,jaer)
          zbsc1(:) = zbsc1(:) + bscat1(jaer)*ssa1(jaer)*aodsp(:,jlev,jaer)
          zext2(:) = zext2(:) + aodsp(:,jlev,jaer)*qex2(jaer)/qex1(jaer)
          zsca2(:) = zsca2(:) + ssa2(jaer)*aodsp(:,jlev,jaer)*qex2(jaer)/qex1(jaer)
          zbsc2(:) = zbsc2(:) + bscat2(jaer)*ssa2(jaer)*aodsp(:,jlev,jaer)*qex2(jaer)/qex1(jaer)
          knz(:)   = knz(:) + 1
          klast(:) = jaer
         endwhere
        enddo

        aod1(:,jlev) = zext1(:)
        aod2(:,jlev) = zext2(:)

      ! Where exactly one species is present the mixture IS that species, and
      ! taking its constants directly rather than forming (s*tau)/tau keeps the
      ! single-species answer bit-for-bit. The division is not exact in IEEE
      ! arithmetic and would otherwise move the existing dust-only answer by an
      ! ulp for no physical reason.

      ! zext1 is floored the way the three quotients beside it already are. A
      ! lane with no species present has zext1 exactly 0 and knz 0, so the
      ! elsewhere excludes it -- and a `where` masks the ASSIGNMENT, not the
      ! evaluation, so the quotient is still formed and the production profile
      ! traps on it. world-5a0.

        zssa1(:) = 0.
        zbs1(:)  = 0.
        zssa2(:) = 0.
        zbs2(:)  = 0.
        where (knz(:) == 1)
         zssa1(:) = ssa1(klast(:))
         zbs1(:)  = bscat1(klast(:))
         zssa2(:) = ssa2(klast(:))
         zbs2(:)  = bscat2(klast(:))
        elsewhere (knz(:) > 1)
         zssa1(:) = zsca1(:)/MAX(zext1(:),TINY(1.0))
         zbs1(:)  = zbsc1(:)/MAX(zsca1(:),TINY(1.0))
         zssa2(:) = zsca2(:)/MAX(zext2(:),TINY(1.0))
         zbs2(:)  = zbsc2(:)/MAX(zsca2(:),TINY(1.0))
        endwhere

        lcons1(:) = (1.0-zssa1(:)) <= zepsc
        lcons2(:) = (1.0-zssa2(:)) <= zepsc

      ! THE TWO FLOORS BELOW ARE THE lcons TEST, MADE STRUCTURAL. The mask
      ! excludes the conservative lanes because the exact expression is
      ! singular there, and a `where` does not: it masks the ASSIGNMENT and
      ! leaves the compiler free to evaluate the right-hand side on every lane,
      ! which at 1-ssa = 0 is a division by zero and at 1-ssa < 0 a square root
      ! of a negative. Both are trapped by the production profile. On a lane the
      ! mask KEEPS, 1-ssa exceeds zepsc and the u-factor is 1 + 2*b*ssa/(1-ssa),
      ! which is at least 1, so neither floor can bind on a result the model
      ! uses. Flooring the u-factor at 1 rather than at 0 is what keeps the
      ! denominator below it away from zero: it is (u+1)^2 e^t - (u-1)^2 e^-t,
      ! which for u >= 1 and t >= 0 is at least 4u. world-5a0.

        where(losun(:) .and. aod1(:,jlev) > 0. .and. .not. lcons1(:))
         zaeru1(:) = SQRT(MAX(1.0,(1.0-zssa1(:)+2*zbs1(:)*zssa1(:))/MAX(1.0-zssa1(:),zepsc))) ! u-factor band 1
         ztemp1(:) = SQRT(MAX(0.0,(1.0-zssa1(:))*(1.0-zssa1(:)+2*zbs1(:)*zssa1(:))))
         zaertf1(:,jlev) = MIN(25.,(ztemp1(:)*aod1(:,jlev))/(zmu0+zero))  ! effective t band 1
         zaerd1(:,jlev) = (((zaeru1(:)+1.0)**2.0)*EXP(zaertf1(:,jlev)) - ((zaeru1(:)-1.0)**2.0)/EXP(zaertf1(:,jlev))) ! denominator band 1
         zaert1(:,jlev) = (4.0*zaeru1(:))/zaerd1(:,jlev) ! transmission band 1
         zaerr1(:,jlev) = (zaeru1(:) + 1.0)*(zaeru1(:) - 1.0)*(EXP(zaertf1(:,jlev))-EXP(-zaertf1(:,jlev)))/zaerd1(:,jlev) ! reflection band 1

         ! Next do scattered light. The DIFFUSE beam takes zmu00, so every
         ! quantity below it must be the s one: zaertf1s and zaerd1s were
         ! computed here and then never read, and the transmission and
         ! reflection were built from the direct-beam zaertf1 and zaerd1, which
         ! made the diffuse stream bit-identical to the direct one and dropped
         ! the factor-of-2 diffusivity the scattered beam carries. zaeru1 is a
         ! function of the single-scattering albedo and the backscatter ratio
         ! only, so it is shared by both beams and is correct as it stands.
         zaertf1s(:,jlev) = MIN(25.,(ztemp1(:)*aod1(:,jlev))/zmu00)  ! effective t band 1 using zmu00 not zmu0!
         zaerd1s(:,jlev) = (((zaeru1(:)+1.0)**2.0)*EXP(zaertf1s(:,jlev)) - ((zaeru1(:)-1.0)**2.0)*EXP(-zaertf1s(:,jlev))) ! denominator band 1
         zaert1s(:,jlev) = (4.0*zaeru1(:))/zaerd1s(:,jlev) ! transmission band 1
         zaerr1s(:,jlev) = (zaeru1(:) + 1.0)*(zaeru1(:) - 1.0)*(EXP(zaertf1s(:,jlev))-EXP(-zaertf1s(:,jlev)))/zaerd1s(:,jlev) ! reflection band 1
        endwhere

      ! CONSERVATIVE SCATTERING, band 1. As ssa goes to 1 the u-factor diverges
      ! and ztemp goes to zero while their product stays finite: u*t tends to
      ! 2*b*tau/mu, so the two-stream collapses to T = 1/(1+b*tau/mu) and
      ! R = (b*tau/mu)/(1+b*tau/mu). That is the limit of the expressions
      ! above, not a different scheme.
      !
      ! The diffuse beam takes the same limit at mu = zmu00, so it gets its own
      ! b*tau/mu and not a copy of the direct answer. Copying was what this
      ! branch did when it was written, which made the defect above it read as
      ! deliberate rather than as the transcription it was.

        where(losun(:) .and. aod1(:,jlev) > 0. .and. lcons1(:))
         ztcon(:) = zbs1(:)*aod1(:,jlev)/(zmu0+zero)
         zaertf1(:,jlev) = 0.
         zaerd1(:,jlev) = 1.0
         zaert1(:,jlev) = 1.0/(1.0+ztcon(:))
         zaerr1(:,jlev) = ztcon(:)/(1.0+ztcon(:))
         ztcon(:) = zbs1(:)*aod1(:,jlev)/zmu00
         zaertf1s(:,jlev) = 0.
         zaerd1s(:,jlev) = 1.0
         zaert1s(:,jlev) = 1.0/(1.0+ztcon(:))
         zaerr1s(:,jlev) = ztcon(:)/(1.0+ztcon(:))
        endwhere

        where(losun(:) .and. aod1(:,jlev) > 0. .and. .not. lcons2(:))
         zaeru2(:) = SQRT(MAX(1.0,(1.0-zssa2(:)+2*zbs2(:)*zssa2(:))/MAX(1.0-zssa2(:),zepsc))) ! u-factor band 2
         ztemp2(:) = SQRT(MAX(0.0,(1.0-zssa2(:))*(1.0-zssa2(:)+2*zbs2(:)*zssa2(:))))
         zaertf2(:,jlev) = MIN(25.,(ztemp2(:)*aod2(:,jlev))/(zmu0+zero)) ! effective t band 2
         zaerd2(:,jlev) = (((zaeru2(:)+1.0)**2.0)*EXP(zaertf2(:,jlev)) - ((zaeru2(:)-1.0)**2.0)/EXP(zaertf2(:,jlev))) ! denominator band 2
         zaert2(:,jlev) = (4.0*zaeru2(:))/zaerd2(:,jlev) ! transmission band 2
         zaerr2(:,jlev) = (zaeru2(:) + 1.0)*(zaeru2(:) - 1.0)*(EXP(zaertf2(:,jlev))-EXP(-zaertf2(:,jlev)))/zaerd2(:,jlev) ! reflection band 2

         zaertf2s(:,jlev) = MIN(25.,(ztemp2(:)*aod2(:,jlev))/zmu00) ! effective t band 2
         zaerd2s(:,jlev) = (((zaeru2(:)+1.0)**2.0)*EXP(zaertf2s(:,jlev)) - ((zaeru2(:)-1.0)**2.0)*EXP(-zaertf2s(:,jlev))) ! denominator band 2
         zaert2s(:,jlev) = (4.0*zaeru2(:))/zaerd2s(:,jlev) ! transmission band 2
         zaerr2s(:,jlev) = (zaeru2(:) + 1.0)*(zaeru2(:) - 1.0)*(EXP(zaertf2s(:,jlev))-EXP(-zaertf2s(:,jlev)))/zaerd2s(:,jlev) ! reflection band 2
        endwhere

        where(losun(:) .and. aod1(:,jlev) > 0. .and. lcons2(:))
         ztcon(:) = zbs2(:)*aod2(:,jlev)/(zmu0+zero)
         zaertf2(:,jlev) = 0.
         zaerd2(:,jlev) = 1.0
         zaert2(:,jlev) = 1.0/(1.0+ztcon(:))
         zaerr2(:,jlev) = ztcon(:)/(1.0+ztcon(:))
         ztcon(:) = zbs2(:)*aod2(:,jlev)/zmu00
         zaertf2s(:,jlev) = 0.
         zaerd2s(:,jlev) = 1.0
         zaert2s(:,jlev) = 1.0/(1.0+ztcon(:))
         zaerr2s(:,jlev) = ztcon(:)/(1.0+ztcon(:))
        endwhere

       enddo ! levels loop
       ! if (mypid == NROOT) then
        ! write(nud,*) "Aerosol effective optical depth 1:",zaertf1
        ! write(nud,*) "Aerosol denominator 1:",zaerd1
        ! write(nud,*) "Aerosol direct transmission 1:",zaert1
        ! write(nud,*) "Aerosol direct reflection 1:",zaerr1
        ! write(nud,*) "Aerosol diffuse transmission 1:",zaert1s
        ! write(nud,*) "Aerosol diffuse reflection 1:",zaerr1s
       ! endif
      endif
!
!     compute optical properties
!
!     downward loop
!
!     preset
!
      where(losun(:))
       zta1(:)=1.
       zta1s(:)=1.
       zra1(:)=0.
       zra1s(:)=0.
       zta2(:)=1.
       zta2s(:)=1.
       zra2(:)=0.
       zra2s(:)=0.
!
       zto3t(:)=1.
       zo3(:)=zxo3t(:)+zmbar*zo3t(:)
       zto3tu(:)=1.                                                     &
     &          -(o3visw*0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)   &
     &           +o3uvw*1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)               &
     &           +o3uvw*0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1
       ztwvt(:)=1.
       zwv(:)=zywvt(:)+zbetta*zwvt(:)
       ztwvtu(:)=1.-h2osww*h2oswl*2.9*zwv(:)                            &
     &            /((1.+141.5*zwv(:))**0.635+5.925*zwv(:))              &
     &            /zsolar2
       ztco2t(:)=1.
       zco2(:)=zyco2t(:)+zbetta*zco2t(:)
       ztco2tu(:)=1.-co2sww*(zca1*LOG(1.+zcb1*zco2(:))                  &
     &                      +zca2*LOG(1.+zcb2*zco2(:)))/zsolar2
!
!     clear sky scattering (Rayleigh scatterin lower most level only)
!
       zrcs(:,NLEV) = zrcs(:,NLEV)*newrsc
       zrcsu(:,NLEV) = zrcsu(:,NLEV)*newrsc
!
!      R = 1 - e^((ps/p0)*ln(T0))
!
       zscf(:) = rcoeff*dp(:)/101100.0*9.80665/ga
       zrcsu(:,NLEV)=zrcsu(:,NLEV) + (1.0-exp(zscf(:)*log(1.0-0.144))) &
     &                                *(1-newrsc)*nrscat*(1-dcc(:,NLEV)*nclouds)
       zrcs(:,NLEV)= zrcs(:,NLEV) + (1.0-exp(zscf(:)*log(1.0-(0.219/(1.+0.816*max(0.,zmu0(:)))))))&
     &                              * zcs(:)*(1-newrsc)*nrscat                    &
     &                            + (1.0-exp(zscf(:)*log(1.0-0.144)))*(1.-zcs(:)-dcc(:,NLEV))&
     &                              * nrscat*(1-newrsc)
       
      endwhere
!
      do jlev=1,NLEV
       where(losun(:))
        zt1(:,jlev)=zta1(:)
        zt2(:,jlev)=zta2(:)
        zr1s(:,jlev)=zra1s(:)
        zr2s(:,jlev)=zra2s(:)
!
!     set single layer R and T:
!
!     1. spectral range 1:
!
!     a) R
!     clear part: rayleigh scattering (only lowermost level)
!     cloudy part: cloud albedo
!     aerosols: reflected direct light (zaerr1) and reflected scattered light (zaerr1s)
!     in clear sky portion only (i.e. (1-dcc))
!
        zrb1(:,jlev)=zrcs(:,jlev)+zrcl1(:,jlev)*dcc(:,jlev)*nclouds+zaerr1(:,jlev)*(1.-dcc(:,jlev))*iaeron
        !zrb1(:,jlev) = zta1*zrcs(:,jlev)+(1-zta1)*zrcsu(:,jlev)+zrcl1(:,jlev)*dcc(:,jlev)
        zrb1s(:,jlev)=zrcsu(:,jlev)+zrcl1s(:,jlev)*dcc(:,jlev)*nclouds+zaerr1s(:,jlev)*(1.-dcc(:,jlev))*iaeron
!
!     b) T
!
!     ozon absorption
!
!     downward beam
!
        zo3(:)=zxo3l(:,jlev)
        zto3(:)=(1.                                                     &
     &          -(o3visw*0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)   &
     &           +o3uvw*1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)               &
     &           +o3uvw*0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1)        &
     &         /zto3t(:)
        zto3t(:)=zto3t(:)*zto3(:)
!
!     upward scattered beam
!
        zo3(:)=zxo3t(:)+zmbar*(zo3t(:)-zo3l(:,jlev))
        zto3u(:)=zto3tu(:)                                              &
     &         /(1.-(o3visw*0.02118*zo3(:)/(1.+0.042*zo3(:)+0.000323*zo3(:)**2)&
     &              +o3uvw*1.082*zo3(:)/((1.+138.6*zo3(:))**0.805)            &
     &              +o3uvw*0.0658*zo3(:)/(1.+(103.6*zo3(:))**3))/zsolar1)
        zto3tu(:)=zto3tu(:)/zto3u(:)
!
!     total T = 1-(A(ozon)+R(rayl.))*(1-dcc)-R(cloud)*dcc
!
!
!     Band 1 aerosol ABSORPTION. zaert1 and zaert1s were computed above and then
!     never used: upstream applies the aerosol to band 2's transmission and to
!     band 1's reflection only, so band 1 scattered but could not absorb. Band 1
!     carries a large share of this star's flux and mineral dust absorbs there,
!     so leaving it out biases atmospheric shortwave absorption low, which is the
!     term the whole precipitation response runs through.
!
!     Booked to band 1's own convention rather than band 2's. Here T is defined
!     as 1 - A - R with the reflection already subtracted through zrb1, so what
!     is removed is the absorption alone, 1 - T_aer - R_aer. Band 2 subtracts
!     1 - T_aer and adds R_aer back through zrb2. Both conserve; they differ
!     only in where the reflected part is booked.
!
        ztb1(:,jlev)=1.-(1.-zto3(:))*(1.-dcc(:,jlev))-zrb1(:,jlev)              &
     &              -(1.-zaert1(:,jlev)-zaerr1(:,jlev))*(1.-dcc(:,jlev))*iaeron
        ztb1u(:,jlev)=1.-(1.-zto3u(:))*(1.-dcc(:,jlev))-zrb1s(:,jlev)           &
     &               -(1.-zaert1s(:,jlev)-zaerr1s(:,jlev))*(1.-dcc(:,jlev))*iaeron
!
!     make combined layer R_ab, R_abs, T_ab and T_abs
!
        z1mrabr(:)=1./(1.-zra1s(:)*zrb1s(:,jlev))
        zra1(:)=zra1(:)+zta1(:)*zrb1(:,jlev)*zta1s(:)*z1mrabr(:)
        zta1(:)=zta1(:)*ztb1(:,jlev)*z1mrabr(:)
        zra1s(:)=zrb1s(:,jlev)+ztb1u(:,jlev)*zra1s(:)*ztb1(:,jlev)      &
     &                        *z1mrabr(:)
        zta1s(:)=ztb1u(:,jlev)*zta1s(:)*z1mrabr(:)
!
!     2. spectral range 2:
!
!     a) R
!
!     cloud albedo
!     aerosol scattering from clear sky part
!
        zrb2(:,jlev)=zrcl2(:,jlev)*dcc(:,jlev)*nclouds+zaerr2(:,jlev)*(1.-dcc(:,jlev))*iaeron
        zrb2s(:,jlev)=zrcl2s(:,jlev)*dcc(:,jlev)*nclouds+zaerr2s(:,jlev)*(1.-dcc(:,jlev))*iaeron
!
!     b) T
!
!     water vapor absorption
!
!     downward beam
!
       zwv(:)=zywvl(:,jlev)
       ztwv(:)=(1.-h2osww*h2oswl*2.9*zwv(:)                             &
     &            /((1.+141.5*zwv(:))**0.635+5.925*zwv(:))              &
     &            /zsolar2)                                             &
     &        /ztwvt(:)
       ztwvt(:)=ztwvt(:)*ztwv(:)
!
!     CO2 absorption, downward beam
!
       zco2(:)=zyco2l(:,jlev)
       ztco2(:)=(1.-co2sww*(zca1*LOG(1.+zcb1*zco2(:))                   &
     &                     +zca2*LOG(1.+zcb2*zco2(:)))/zsolar2)         &
     &         /ztco2t(:)
       ztco2t(:)=ztco2t(:)*ztco2(:)
!
!     upward scattered beam
!
       zwv(:)=zywvt(:)+zbetta*(zwvt(:)-zwvl(:,jlev))
       ztwvu(:)=ztwvtu(:)                                               &
     &         /(1.-h2osww*h2oswl*2.9*zwv(:)                            &
     &            /((1.+141.5*zwv(:))**0.635+5.925*zwv(:))              &
     &            /zsolar2)
       ztwvtu(:)=ztwvtu(:)/ztwvu(:)
!
!     CO2 absorption, upward scattered beam
!
       zco2(:)=zyco2t(:)+zbetta*(zco2t(:)-zco2l(:,jlev))
       ztco2u(:)=ztco2tu(:)                                             &
     &          /(1.-co2sww*(zca1*LOG(1.+zcb1*zco2(:))                  &
     &                      +zca2*LOG(1.+zcb2*zco2(:)))/zsolar2)
       ztco2tu(:)=ztco2tu(:)/ztco2u(:)
!
!     total T = 1-A(water vapor)*(1.-dcc)-(A(cloud)+R(cloud))*dcc
!
        ztb2(:,jlev)=1.-(1.-ztwv(:))*(1.-dcc(:,jlev)*nclouds)                   &
     &              -(1.-ztco2(:))*(1.-dcc(:,jlev)*nclouds)                     &
     &              -(1.-ztcl2(:,jlev))*dcc(:,jlev)*nclouds                     &
                    -(1.-zaert2(:,jlev))*(1.-dcc(:,jlev))*iaeron
        ztb2u(:,jlev)=1.-(1.-ztwvu(:))*(1.-dcc(:,jlev)*nclouds)                 &
     &               -(1.-ztco2u(:))*(1.-dcc(:,jlev)*nclouds)                   &
     &               -(1.-ztcl2s(:,jlev))*dcc(:,jlev)*nclouds                   &
                     -(1.-zaert2s(:,jlev))*(1.-dcc(:,jlev))*iaeron
!
!     make combined layer R_ab, R_abs, T_ab and T_abs
!
        z1mrabr(:)=1./(1.-zra2s(:)*zrb2s(:,jlev))
        zra2(:)=zra2(:)+zta2(:)*zrb2(:,jlev)*zta2s(:)*z1mrabr(:)
        zta2(:)=zta2(:)*ztb2(:,jlev)*z1mrabr(:)
        zra2s(:)=zrb2s(:,jlev)+ztb2u(:,jlev)*zra2s(:)*ztb2(:,jlev)      &
     &                       *z1mrabr(:)
        zta2s(:)=ztb2u(:,jlev)*zta2s(:)*z1mrabr(:)
       endwhere
      enddo
      where(losun(:))
       zt1(:,NLEP)=zta1(:)
       zt2(:,NLEP)=zta2(:)
       zr1s(:,NLEP)=zra1s(:)
       zr2s(:,NLEP)=zra2s(:)
!
!     upward loop
!
!     make upward R
!

! Currently: we use the same albedo for both spectral ranges.

       zra1s(:)=dalb(:)*(1-nstartemp) + dsalb(1,:)*nstartemp
       zra2s(:)=dalb(:)*(1-nstartemp) + dsalb(2,:)*nstartemp
       
!
!      set albedo for the direct beam (for ocean use ECHAM3 param unless necham=0)
       dsalb(1,:)=dls(:)*dsalb(1,:)   +   (1.-dls(:)) * dicec(:)*dsalb(1,:)              &
     &           + (1.-dls(:)) * (1.-dicec(:)) * AMIN1(0.05/(zmu0(:)+0.15),0.15)*necham*(1-necham6) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1-necham)*necham6 &
     &  *(0.026/(zmu0(:)**1.7+0.065)+0.15*(zmu0(:)-1)*(zmu0(:)-0.5)*(zmu0(:)-0.1)+0.0082) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1.-necham)*(1-necham6)*dsalb(1,:)
       dsalb(2,:)=dls(:)*dsalb(2,:)   +   (1.-dls(:)) * dicec(:)*dsalb(2,:)              &
     &           + (1.-dls(:)) * (1.-dicec(:)) * AMIN1(0.05/(zmu0(:)+0.15),0.15)*necham*(1-necham6) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1-necham)*necham6 &
     &  *(0.026/(zmu0(:)**1.7+0.065)+0.15*(zmu0(:)-1)*(zmu0(:)-0.5)*(zmu0(:)-0.1)+0.0082) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1.-necham)*(1-necham6)*dsalb(2,:)
       
       dalb(:) = (zsolars(1)*dsalb(1,:) + zsolars(2)*dsalb(2,:))*nstartemp  &
     &           + (dls(:)*dalb(:)      + (1.-dls(:)) * dicec(:)*dalb(:)              &
     &           + (1.-dls(:)) * (1.-dicec(:)) * AMIN1(0.05/(zmu0(:)+0.15),0.15)*necham*(1-necham6) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1-necham)*necham6 &
     &  *(0.026/(zmu0(:)**1.7+0.065)+0.15*(zmu0(:)-1)*(zmu0(:)-0.5)*(zmu0(:)-0.1)+0.0082) &
     &           + (1.-dls(:)) * (1.-dicec(:)) * (1.-necham)*(1-necham6)*dalb(:))*(1-nstartemp)
       
     
       zra1(:)=dsalb(1,:)*nstartemp + dalb(:)*(1-nstartemp)
       zra2(:)=dsalb(2,:)*nstartemp + dalb(:)*(1-nstartemp)
         
! Ice-free ocean albedo is min(0.05/(phi+0.15), 0.15)--reflection and scattering is higher at low phi
       
      endwhere
      do jlev=NLEV,1,-1
       where(losun(:))
        zrl1(:,jlev+1)=zra1(:)
        zrl2(:,jlev+1)=zra2(:)
        zrl1s(:,jlev+1)=zra1s(:)
        zrl2s(:,jlev+1)=zra2s(:)
        zra1(:)=zrb1(:,jlev)+ztb1(:,jlev)*zra1(:)*ztb1u(:,jlev)         &
     &                      /(1.-zra1s(:)*zrb1s(:,jlev))
        zra1s(:)=zrb1s(:,jlev)+ztb1u(:,jlev)*zra1s(:)*ztb1u(:,jlev)     &
     &                        /(1.-zra1s(:)*zrb1s(:,jlev))
        zra2(:)=zrb2(:,jlev)+ztb2(:,jlev)*zra2(:)*ztb2u(:,jlev)         &
     &                      /(1.-zra2s(:)*zrb2s(:,jlev))
        zra2s(:)=zrb2s(:,jlev)+ztb2u(:,jlev)*zra2s(:)*ztb2u(:,jlev)     &
     &                        /(1.-zra2s(:)*zrb2s(:,jlev))
       endwhere
      enddo
      where(losun(:))
       zrl1(:,1)=zra1(:)
       zrl2(:,1)=zra2(:)
       zrl1s(:,1)=zra1s(:)
       zrl2s(:,1)=zra2s(:)
      endwhere
!
!     fluxes at layer interfaces
!
      do jlev=1,NLEP
       where(losun(:))
        z1mrabr(:)=1./(1.-zr1s(:,jlev)*zrl1s(:,jlev))
        zfd1(:)=zt1(:,jlev)*z1mrabr(:)
        zfu1(:)=-zt1(:,jlev)*zrl1(:,jlev)*z1mrabr(:)
        z1mrabr(:)=1./(1.-zr2s(:,jlev)*zrl2s(:,jlev))
        zfd2(:)=zt2(:,jlev)*z1mrabr(:)
        zfu2(:)=-zt2(:,jlev)*zrl2(:,jlev)*z1mrabr(:)
        dfu(:,jlev)=zfu1(:)*zftop1(:)+zfu2(:)*zftop2(:)
        dfd(:,jlev)=zfd1(:)*zftop1(:)+zfd2(:)*zftop2(:)
!       THE BAND SPLIT, KEPT RATHER THAN SUMMED AWAY. WORLD-3QFZ. Each band's
!       contribution to `dfd` above, stored at every level so that what survives
!       the loop is the last one, NLEP, the surface. `dfd` keeps its own
!       statement unchanged so the sum is contracted exactly as it was.
        dfdsw1(:)=zfd1(:)*zftop1(:)
        dfdsw2(:)=zfd2(:)*zftop2(:)
        dswfl(:,jlev)=dfu(:,jlev)+dfd(:,jlev)
       endwhere
      enddo
!
      return
      end subroutine swr

!     ==============
!     SUBROUTINE LWR
!     ==============

      subroutine lwr
      use radmod
!
!     compute long wave radiation
!
!     o3-, co2-, h2o- and cloud-absorption is considered
!
!     clear sky absorptivities from Sasamori 1968
!     (J. Applied Meteorology, 7, 721-729)
!
!     no PUMA *subs* are used
!
!     the following PUMA variables are used/modified:
!
!     ga               : gravity accelleration (m/s2) (used)
!     sigma(NLEV)      : full level sigma  (used)
!     sigmah(NLEV)     : half level sigma  (used)
!     dsigma(NLEV)     : delta sigma (half level)  (used)
!     dp(NHOR)         : surface pressure (Pa) (used)
!     dq(NHOR,NLEP)    : specific humidity (kg/kg) (used)
!     dt(NHOR,NLEP)    : temperature (K) (used)
!     dcc(NHOR,NLEP)   : cloud cover (frac.) (used)
!     dlwfl(NHOR,NLEP) : long wave radiation (W/m2)  (modified)
!     dftu(NHOR,NLEP)  : long wave radiation upward (W/m2) (modified)
!     dftd(NHOR,NLEP)  : long wave radiation downward (W/m2) (modified)
!
!     0) define local parameters
!
      parameter(zmmair=0.0289644)      ! molecular weight air (kg/mol)
      parameter(zmmco2=0.0440098)      ! molecular weight co2 (kg/mol)
      parameter(zpv2pm=zmmco2/zmmair)  ! transfere co2 ppvol to ppmass
      parameter(zrco2=1.9635)          ! co2 density (kg/m3 stp)
      parameter(zro3=2.14)             ! o3 density (kg/m3 stp)
      parameter(zttop=0.)              ! t at top of atmosphere
!
!     scaling factors for transmissivities:
!
!     uh2o in g/cm**2 (zfh2o=0.1=1000/(100*100) *kg/m**2)
!     uco2 in cm-STP  (zfco2=100./rco2 * kg/m**2)
!     uo3  in cm-STP  (zfo3=100./ro3 * kg/m**2)
!
      parameter(zfh2o=0.1)
      parameter(zfco2=100./zrco2)
      parameter(zfo3=100./zro3)
      parameter(zt0=295.)
!
!**   local arrays
!
      real zbu(NHOR,0:NLEP)     ! effective SBK*T**4 for upward radiation
      real zbd(NHOR,0:NLEP)     ! effective SBK*T**4 for downward radiation
      real zst4h(NHOR,NLEP)     ! SBK*T**4  on half levels
      real zst4(NHOR,NLEP)      ! SBK*T**4  on ull levels
      real ztau(NHOR,NLEV)      ! total transmissivity
      real zq(NHOR,NLEV)        ! modified water vapour
      real zqo3(NHOR,NLEV)      ! modified ozon
      real zqco2(NHOR,NLEV)     ! modified co2
      real ztausf(NHOR,NLEV)    ! total transmissivity to surface
      real ztaucs(NHOR,NLEV)    ! clear sky transmissivity
      real ztaucc0(NHOR,NLEV)   ! layer transmissivity cloud
      real ztaucc(NHOR)         ! cloud transmissivity
      real ztaudu0(NHOR,NLEV)   ! layer transmissivity dust
      real ztaudu(NHOR)         ! dust transmissivity
!
!     The aerosol's thermal-IR absorption is additive IN THE EXPONENT, so N
!     species cost a sum and no restructuring: the grey absorbers overlap by
!     the same random-overlap assumption the cloud term already makes.
!
      real zqsum(NHOR)          ! sum of qlw*od over species
      integer :: jaer           ! aerosol species index
!
!     CH4 AND N2O. CLIM-42; exoplasim/notes/trace-gas-band.md is the argument
!     and exoplasim/analysis/trace_gas_band_model.json carries these constants
!     with the checks that produced them. Do not retype them: the generator is
!     `exoplasim/scripts/trace_gas_band_model.py --fortran`.
!
!     Donner and Ramanathan (1980), the Cess and Ramanathan (1972) band model
!     as modified by Ramanathan (1976):
!
!         A(U,beta) = 2 A0 ln[ 1 + U / sqrt(4 + U (1 + 1/beta)) ]
!         U = S W / A0            beta = beta0 (P / P0)
!
!     A0, beta0 and S are a MATCHED TRIPLE -- A0 and beta0 were fitted through
!     that equation at a particular S -- so no one of them is swapped for a
!     newer value on its own. CH4's S is recovered from the paper's own Table
!     2, which pins it; N2O's two come from McClatchey et al. (1973), which is
!     the compilation the paper names.
      parameter(zch4a0 =  52.000)   ! CH4 1306 cm-1 bandwidth at 300 K, cm-1
      parameter(zch4be =   0.170)   ! CH4 1306 cm-1 line shape at 300 K, 1 atm
      parameter(zch4si = 187.690)   ! CH4 1306 cm-1 intensity, cm-1 (cm atm)-1
      parameter(zch4wn =1306.000)   ! band centre, cm-1
      parameter(zn2aa0 =  20.400)   ! N2O 1285 cm-1 bandwidth at 300 K
      parameter(zn2abe =   1.120)   ! N2O 1285 cm-1 line shape at 300 K
      parameter(zn2asi = 267.605)   ! N2O 1285 cm-1 intensity
      parameter(zn2awn =1285.000)
      parameter(zn2ba0 =  23.000)   ! N2O 589 cm-1 bandwidth at 300 K
      parameter(zn2bbe =   1.080)   ! N2O 589 cm-1 line shape at 300 K
      parameter(zn2bsi =  31.704)   ! N2O 589 cm-1 intensity
      parameter(zn2bwn = 589.000)
!
!     Volume mixing ratio in ppmv to absorber amount in cm at STP, per kg/m2 of
!     AIR. The gas's molecular weight CANCELS -- it appears once converting a
!     volume ratio to a mass ratio and once again in the STP density -- so this
!     one factor serves every gas, and the CO2 path above computes the same
!     number the long way round as zfco2*zzf2.
      parameter(zvol2cm = 7.7384E-5)
!
!     THE CO2 OVERLAP AT 589 cm-1, which is not optional. That band sits 78
!     cm-1 from CO2's 667 cm-1 fundamental, well inside the 15 um band, so
!     carrying it without the overlap credits N2O with absorption CO2 already
!     provides -- an OVERSTATEMENT, and not the conservative error that
!     dropping the band would be.
!
!     ln(tau) is a fit to a correlated-k band mean over 546-630 cm-1 from the
!     LMD Generic PCM's HITRAN 2020 tables, corrected onto the 566-612 cm-1
!     window this band occupies, in ln(u), ln(P/atm) and 1/T. Worst error on
!     the resulting multiplier is 0.061 and that is beyond this world's whole
!     column; over the column it costs at most 0.05 W/m2 on a term whose own
!     bracket is 0.42. exoplasim/analysis/co2_overlap_589.json has the fit,
!     the cross-check against Ramanathan (1976) Appendix A, and the reason the
!     two are compared as absorptances rather than as transmissivities.
      real, parameter :: zovl(7) = (/ -1.2876101, 0.61550971, -0.033606505,   &
     &                    0.18164327, 0.016230565, -740.96407, 46.117656 /)
!
      real zqch4(NHOR,NLEV)     ! CH4 amount per layer, cm STP
      real zqn2o(NHOR,NLEV)     ! N2O amount per layer, cm STP
      real zqair(NHOR,NLEV)     ! layer air mass, kg/m2
      real zplay(NHOR,NLEV)     ! layer pressure, atm
      real zsch4(NHOR)          ! accumulated CH4 amount
      real zsn2o(NHOR)          ! accumulated N2O amount
      real zsup(NHOR)           ! accumulated amount * pressure, for the mean
      real zsut(NHOR)           ! accumulated amount * temperature
      real zsco2c(NHOR)         ! accumulated CO2 amount, cm STP, UNWEIGHTED
      real zpeff(NHOR)          ! path mean pressure, atm
      real zteff(NHOR)          ! path mean temperature, K
      real zach4(NHOR)          ! CH4 absorptivity
      real zan2oa(NHOR)         ! N2O 1285 absorptivity
      real zan2ob(NHOR)         ! N2O 589 absorptivity
      real ztco2b(NHOR)         ! CO2 transmissivity in the 589 region
      real zsair(NHOR)          ! accumulated air mass along the path
      real zta0(NHOR)           ! bandwidth parameter at the path temperature
      real ztbe(NHOR)           ! line shape parameter at the path pressure
      real zuu(NHOR)            ! dimensionless optical pathlength
      real zaa(NHOR)            ! band absorptance, cm-1
      real zbb(NHOR)            ! Planck radiance at the band centre
      real zlu(NHOR)            ! ln of the CO2 amount, for the overlap fit
      real zlp(NHOR)            ! ln of the path pressure
      real ztauov(NHOR)         ! CO2 optical depth in the 589 region
      real ztau0(NHOR)          ! approx. layer transmissivity
      real zsumwv(NHOR)         ! effective water vapor amount
      real zsumo3(NHOR)         ! effective o3 amount
      real zsumco2(NHOR)        ! effective co2 amount
      real zsfac(NHOR)          ! scaling factor
      real zah2o(NHOR)          ! water vapor absorptivity
      real zaco2(NHOR)          ! co2 absorptivity
      real zao3(NHOR)           ! o3 absorptivity
      real zth2o(NHOR)          ! water vapor - co2 overlap transmissivity
      real zbdl(NHOR)           ! layer evective downward rad.
      real zeps(NHOR)           ! surface emissivity
      real zps2(NHOR)           ! ps**2
      real zsigh2(NLEP)         ! sigmah**2
!
!     entropy
!
      real zbue1(NHOR,0:NLEP),zbue2(NHOR,0:NLEP)

!
!**   1) set some necessary (helpful) bits
!

      zero=1.E-6                     ! small number

      zao30= 0.209*(7.E-5)**0.436    ! to get a(o3)=0 for o3=0
      zco20=0.0676*(0.01022)**0.421  ! to get a(co2)=0 for co2=0
      zh2o0a=0.846*(3.59E-5)**0.243  ! to get a(h2o)=0 for h2o=0
      zh2o0=0.832*0.0286**0.26       ! to get t(h2o)=1 for h2o=0
!
!     to make a(o3) continues at 0.01cm: 
!
      zao3c=0.209*(0.01+7.E-5)**0.436-zao30-0.0212*log10(0.01) 
!
!     to make a(co2) continues at 1cm:
!
      zaco2c=0.0676*(1.01022)**0.421-zco20
!
!     to make a(h2o) continues at 0.01gm:
!
      zah2oc=0.846*(0.01+3.59E-5)**0.243-zh2o0a-0.24*ALOG10(0.02)
!
!     to make t(h2o) continues at 2gm :
!
      zth2oc=1.-(0.832*(2.+0.0286)**0.26-zh2o0)+0.1196*log(2.-0.6931)

      zsigh2(1)=0.
      zsigh2(2:NLEP)=sigmah(1:NLEV)**2

      if (npbroaden .gt. 0.5) then
         zps2(:)=dp(:)*dp(:)
      else
         zps2(:)=101100.0*101100.0 !If no pressure broadening, assume as much broadening as for 1 bar
      endif
!
!**   2) calc. stb*t**4 + preset fluxes
!
      dftu(:,:)=0.
      dftd(:,:)=0.
!
!*    stb*t**4 on full and half levels
!
!     full levels and surface
!
      zst4(:,1:NLEP)=SBK*dt(:,1:NLEP)**4
!
!     half level (incl. toa and near surface)
!
      zst4h(:,1)=zst4(:,1)                                              &
     &          -(zst4(:,1)-zst4(:,2))*sigma(1)/(sigma(1)-sigma(2))
      do jlev=2,NLEV
       jlem=jlev-1
       zst4h(:,jlev)=(zst4(:,jlev)*(sigma(jlem)-sigmah(jlem))           &
     &               +zst4(:,jlem)*(sigmah(jlem)-sigma(jlev)))          &
     &              /(sigma(jlem)-sigma(jlev))
      enddo
      where((zst4(:,NLEV)-zst4h(:,NLEV))                                &
     &     *(zst4(:,NLEV)-zst4(:,NLEP)) > 0.)
       zst4h(:,NLEP)=zst4(:,NLEP)
      elsewhere
       zst4h(:,NLEP)=zst4(:,NLEV)                                       &
     &              +(zst4(:,NLEV)-zst4(:,NLEM))*(1.-sigma(NLEV))       &
     &              /(sigma(NLEV)-sigma(NLEM))
      endwhere
!
!*    top downward flux, surface grayness and surface upward flux
!
      zbd(:,0)=SBK*zttop**4 !We could add IR flux from M dwarf host star here? --AYP
      zeps(:)=elwland*dls(:)+elwsea*(1.-dls(:))
      zbu(:,NLEP)=zeps(:)*zst4(:,NLEP)
      zbue1(:,NLEP)=0.
      zbue2(:,NLEP)=zbu(:,NLEP)
!
!**   3) vertical loops
!
!
!     a) modified absorber amounts and cloud transmissivities
!
      do jlev=1,NLEV
       jlep=jlev+1
       zzf1=sigma(jlev)*dsigma(jlev)/ga/100000.
       zzf2=zpv2pm*1.E-6                        !get co2 in pp mass (kg/kg-stp)
       zzf3=-1.66*acllwr*1000.*dsigma(jlev)/ga
       zsfac(:)=zzf1*zps2(:)
       zq(:,jlev)=zfh2o*zsfac(:)*dq(:,jlev)
       zqo3(:,jlev)=zfo3*zsfac(:)*dqo3(:,jlev)
       zqco2(:,jlev)=zfco2*zzf2*zsfac(:)*dqco2(:,jlev) 
!
!     The three amounts above are PRESSURE-WEIGHTED: zsfac carries zps2, so
!     what Sasamori's fits receive is an amount times a pressure. The band
!     model below wants the two SEPARATELY, because its beta = beta0 P/P0 is
!     where pressure enters, so these are the plain amounts and the plain
!     layer pressure beside them. CLIM-42.
!
       zqair(:,jlev)=dsigma(jlev)*dp(:)/ga
       zplay(:,jlev)=sigma(jlev)*dp(:)/101325.
       zqch4(:,jlev)=zvol2cm*ch4*zqair(:,jlev)
       zqn2o(:,jlev)=zvol2cm*n2o*zqair(:,jlev)
       if(clgray > 0) then
        ztaucc0(:,jlev)=1.-dcc(:,jlev)*clgray
       else
        ztaucc0(:,jlev)=1.-dcc(:,jlev)*(1.-exp(zzf3*dql(:,jlev)*dp(:)))
       endif
!
!     PRESCRIBED DUST in the longwave.
!
!     ExoPlaSim's aerosol acts in the two shortwave bands only and this solver
!     had no aerosol term of any kind, so the model could cool with dust and
!     could not warm with it. For this world that is not a refinement: the
!     shortwave-only error is of the same order as the whole quantity.
!
!     Dust enters exactly where cloud does, as an additional layer
!     transmissivity multiplied into the total. It is a grey ABSORBER -- no
!     longwave scattering, which for a micron-scale particle in the thermal
!     infrared is the standard approximation and is stated rather than assumed
!     -- at the same 1.66 diffusivity the cloud term uses.
!
!     Because it multiplies rather than adds, the overlap with water vapour and
!     CO2 is taken care of by construction: where the gas is already opaque the
!     dust adds nothing. That is the part the offline estimate in
!     analysis/dust_forcing.json has to approximate with a declared window
!     transmittance, and it is why the offline longwave is an upper bound.
!
!     The scheme is broadband, so multiplying a grey dust transmissivity into a
!     broadband gaseous absorptivity is a random-overlap assumption. It is the
!     same assumption already made for cloud one line above.
!
!     Every species enters the same term with its OWN absorption ratio and its
!     own per-layer optical depth, and the absorption adds in the exponent. N
!     species therefore cost a sum over aqlw(jaer)*aodsp(:,jlev,jaer) and
!     nothing about the scheme restructures. The prescribed columns and the
!     transported tracer are species alike here: neither replaces the other,
!     and both contribute whenever both are active. CLIM-39.
!
       if (naerosp > 0) then
        zqsum(:) = 0.
        do jaer = 1,naerosp
         zqsum(:) = zqsum(:) + aqlw(jaer)*aodsp(:,jlev,jaer)
        enddo
        ztaudu0(:,jlev)=exp(-1.66*zqsum(:))
       else
        ztaudu0(:,jlev)=1.
       endif
      enddo
!
!     b) transmissivities, effective radiations and fluxes
!
      do jlev=1,NLEV
       jlem=jlev-1
       ztaucc(:)=1.
       ztaudu(:)=1.
       zsumwv(:)=0.
       zsumo3(:)=0.
       zsumco2(:)=0.
       zsch4(:)=0.
       zsn2o(:)=0.
       zsco2c(:)=0.
       zsair(:)=0.
       zsup(:)=0.
       zsut(:)=0.
!
!     transmissivities
!
       do jlev2=jlev,NLEV
        jlep2=jlev2+1
        zsumwv(:)=zsumwv(:)+zq(:,jlev2)
        zsumo3(:)=zsumo3(:)+zqo3(:,jlev2)
        zsumco2(:)=zsumco2(:)+zqco2(:,jlev2)
!
!     clear sky transmisivity
!
!     h2o absorption:
!
!     a) 6.3mu
!
        where(zsumwv(:) <= 0.01)
         zah2o(:)=0.846*max(0.,zsumwv(:)+3.59E-5)**0.243-zh2o0a
        elsewhere
         zah2o(:)=0.24*ALOG10(max(1.E-30,zsumwv(:)+0.01))+zah2oc
        endwhere
!
!     b) continuum
!
        if(th2oc > 0.) then
         zah2o(:)=AMIN1(zah2o(:)+(1.-exp(-th2oc*zsumwv(:))),1.)
        endif
!
!     co2 absorption:
!
        where(zsumco2(:) <= 1.0)
         zaco2(:)=0.0676*max(0.,zsumco2(:)+0.01022)**0.421-zco20
        elsewhere
         zaco2(:)=0.0546*ALOG10(max(1.E-30,zsumco2(:)))+zaco2c
        endwhere
!
!     Boer et al. (1984) scheme for t(h2o) at co2 overlapp
!
        where(zsumwv(:)<= 2.)
         zth2o(:)=1.-(0.832*max(0.,zsumwv(:)+0.0286)**0.26-zh2o0)
        elsewhere
         zth2o(:)=max(0.,zth2oc-0.1196*log(max(1.E-30,zsumwv(:)-0.6931)))
        endwhere
!
!     o3 absorption:
!
        where(zsumo3(:) <= 0.01)
         zao3(:)= 0.209*max(0.,zsumo3(:)+7.E-5)**0.436 - zao30
        elsewhere
         zao3(:)= 0.0212*log10(max(1.E-30,zsumo3(:)))+zao3c
        endwhere
!
!     CH4 and N2O. CLIM-42.
!
!     Each band's absorptance comes out of Donner and Ramanathan's Eq. (1) in
!     cm-1 and has to reach Sasamori's currency, which is a FRACTION of the
!     broadband flux. The conversion is the share of the Planck function the
!     band sits on, pi B(v,T) / (sigma T**4), exact for a narrow band, and it
!     is where this term's temperature dependence enters.
!
!     The path's pressure and temperature are amount-weighted means, which is
!     the Curtis-Godson choice. The existing three absorbers fold pressure in
!     through zps2 instead; the two conventions are kept apart rather than
!     mixed, which is why the amounts above are carried unweighted.
!
!     Water vapour overlaps CH4 1306 and N2O 1285 and is handled by
!     multiplying through zth2o, which is what Donner and Ramanathan do and
!     line for line what this routine already does to CO2. N2O 589 is
!     overlapped by CO2 instead, and that is the fit below.
!
        if (ch4 > 0. .or. n2o > 0.) then
         zsch4(:)=zsch4(:)+zqch4(:,jlev2)
         zsn2o(:)=zsn2o(:)+zqn2o(:,jlev2)
         zsco2c(:)=zsco2c(:)+zvol2cm*dqco2(:,jlev2)*zqair(:,jlev2)
         zsair(:)=zsair(:)+zqair(:,jlev2)
         zsup(:)=zsup(:)+zqair(:,jlev2)*zplay(:,jlev2)
         zsut(:)=zsut(:)+zqair(:,jlev2)*dt(:,jlev2)
         zpeff(:)=zsup(:)/MAX(zsair(:),1.E-30)
         zteff(:)=MAX(zsut(:)/MAX(zsair(:),1.E-30),100.)
!
         if (ch4 > 0.) then
          zta0(:)=zch4a0*SQRT(zteff(:)/300.)
          ztbe(:)=MAX(zch4be*SQRT(300./zteff(:))*zpeff(:),1.E-12)
          zuu(:)=zch4si*zsch4(:)/zta0(:)
          zaa(:)=2.*zta0(:)                                                    &
     &          *LOG(1.+zuu(:)/SQRT(4.+zuu(:)*(1.+1./ztbe(:))))
          zbb(:)=1.191042E-8*zch4wn**3/(EXP(1.4387769*zch4wn/zteff(:))-1.)
          zach4(:)=zaa(:)*PI*zbb(:)/(SBK*zteff(:)**4)*zth2o(:)
         else
          zach4(:)=0.
         endif
!
         if (n2o > 0.) then
          zta0(:)=zn2aa0*SQRT(zteff(:)/300.)
          ztbe(:)=MAX(zn2abe*SQRT(300./zteff(:))*zpeff(:),1.E-12)
          zuu(:)=zn2asi*zsn2o(:)/zta0(:)
          zaa(:)=2.*zta0(:)                                                    &
     &          *LOG(1.+zuu(:)/SQRT(4.+zuu(:)*(1.+1./ztbe(:))))
          zbb(:)=1.191042E-8*zn2awn**3/(EXP(1.4387769*zn2awn/zteff(:))-1.)
          zan2oa(:)=zaa(:)*PI*zbb(:)/(SBK*zteff(:)**4)*zth2o(:)
!
          zlu(:)=LOG(MAX(zsco2c(:),1.E-10))
          zlp(:)=LOG(MAX(zpeff(:),1.E-10))
          ztauov(:)=EXP(zovl(1)+zovl(2)*zlu(:)+zovl(3)*zlu(:)*zlu(:)             &
     &           +zovl(4)*zlp(:)+zovl(5)*zlu(:)*zlp(:)                         &
     &           +zovl(6)/zteff(:)+zovl(7)*zlu(:)/zteff(:))
          ztco2b(:)=EXP(-MIN(ztauov(:),50.))
!
          zta0(:)=zn2ba0*SQRT(zteff(:)/300.)
          ztbe(:)=MAX(zn2bbe*SQRT(300./zteff(:))*zpeff(:),1.E-12)
          zuu(:)=zn2bsi*zsn2o(:)/zta0(:)
          zaa(:)=2.*zta0(:)                                                    &
     &          *LOG(1.+zuu(:)/SQRT(4.+zuu(:)*(1.+1./ztbe(:))))
          zbb(:)=1.191042E-8*zn2bwn**3/(EXP(1.4387769*zn2bwn/zteff(:))-1.)
          zan2ob(:)=zaa(:)*PI*zbb(:)/(SBK*zteff(:)**4)*ztco2b(:)
         else
          zan2oa(:)=0.
          zan2ob(:)=0.
         endif
        else
         zach4(:)=0.
         zan2oa(:)=0.
         zan2ob(:)=0.
        endif
!
!     total clear sky transmissivity
!
        ztaucs(:,jlev2)=1.-zah2o(:)-zao3(:)-zaco2(:)*zth2o(:)                 &
     &                   -zach4(:)-zan2oa(:)-zan2ob(:)
!
!     bound transmissivity:
!
        ztaucs(:,jlev2)=AMIN1(1.-zero,AMAX1(zero,ztaucs(:,jlev2)))
!
!     cloud transmisivity assuming random overlap
!
        ztaucc(:)=ztaucc(:)*ztaucc0(:,jlev2)
!
!     dust transmisivity, accumulated the same way
!
        ztaudu(:)=ztaudu(:)*ztaudu0(:,jlev2)
!
!     total transmissivity
!
        ztau(:,jlev2)=ztaucs(:,jlev2)*ztaucc(:)*ztaudu(:)
       enddo
!
!     upward and downward effective SBK*T**4
!
       if(jlev == 1) then
        ztau0(:)=1.-zero
        do jlev2=1,NLEV
         jlep2=jlev2+1
         ztau0(:)=ztaucs(:,jlev2)/ztau0(:)*ztaucc0(:,jlev2)*ztaudu0(:,jlev2)*tpofmt
         ztau0(:)=AMIN1(1.-zero,MAX(zero,ztau0(:)))
         where((zst4(:,jlev2)-zst4h(:,jlev2))                           &
     &        *(zst4(:,jlev2)-zst4h(:,jlep2)) > 0.)
          zbd(:,jlev2)=0.5*zst4(:,jlev2)                                &
     &                +0.25*(zst4h(:,jlev2)+zst4h(:,jlep2))
          zbu(:,jlev2)=zbd(:,jlev2)
          zbue1(:,jlev2)=zbu(:,jlev2)
          zbue2(:,jlev2)=0.
         elsewhere
          zbd(:,jlev2)=(zst4h(:,jlep2)-ztau0(:)*zst4h(:,jlev2))         &
     &                /(1.-ztau0(:))                                    &
     &                -(zst4h(:,jlev2)-zst4h(:,jlep2))/ALOG(ztau0(:))
          zbu(:,jlev2)=zst4h(:,jlev2)+zst4h(:,jlep2)-zbd(:,jlev2)
          zbue1(:,jlev2)=zbu(:,jlev2)
          zbue2(:,jlev2)=0.
         endwhere
         ztau0(:)=ztaucs(:,jlev2)
        enddo
       endif
!
!     fluxes
!
       dftu(:,jlev)=dftu(:,jlev)-zbu(:,jlev)
       dftd(:,jlev)=dftd(:,jlev)+zbd(:,jlem)
       dftue1(:,jlev)=dftue1(:,jlev)-zbue1(:,jlev)
       dftue2(:,jlev)=dftue2(:,jlev)-zbue2(:,jlev)
       zbdl(:)=zbd(:,jlem)-zbd(:,jlev)
       do jlev2=jlev,NLEV
        jlep2=jlev2+1
        dftu(:,jlev)=dftu(:,jlev)                                       &
     &              -(zbu(:,jlep2)-zbu(:,jlev2))*ztau(:,jlev2)
        dftue1(:,jlev)=dftue1(:,jlev)                                   &
     &              -(zbue1(:,jlep2)-zbue1(:,jlev2))*ztau(:,jlev2)
        dftue2(:,jlev)=dftue2(:,jlev)                                   &
     &              -(zbue2(:,jlep2)-zbue2(:,jlev2))*ztau(:,jlev2)
        dftd(:,jlep2)=dftd(:,jlep2)                                     &
     &               +zbdl(:)*ztau(:,jlev2)
       enddo
       if(jlev == 1) then
        dftu0(:,1)=-zbu(:,1)*(1.-ztau(:,1))
        do jlev2=2,NLEV
         jlem=jlev2-1
         dftu0(:,jlev2)=-zbu(:,jlev2)*(ztau(:,jlem)-ztau(:,jlev2))
        enddo
       endif
!
!     collect transmissivity to surface
!
       ztausf(:,jlev)=ztau(:,NLEV)
      enddo
      do jlev=1,NLEV-1
       jlep=jlev+1
       dftd0(:,jlev)=-zbd(:,jlev)*(ztausf(:,jlep)-ztausf(:,jlev))
      enddo
      dftd0(:,NLEV)=-zbd(:,jlev)*(1.-ztausf(:,NLEV))
!
!     complite surface lwr
!
      dftu(:,NLEP)=dftu(:,NLEP)-zbu(:,NLEP)
      dftd(:,NLEP)=dftd(:,NLEP)+zbd(:,NLEV)
      dftue1(:,NLEP)=dftue1(:,NLEP)-zbue1(:,NLEP)
      dftue2(:,NLEP)=dftue2(:,NLEP)-zbue2(:,NLEP)
!
!     correct for non black surface
!
!     THE REFLECTED FLUX HAS TO LEAVE THE SURFACE LEVEL AS WELL AS PASS THE
!     ATMOSPHERE. `ztausf` is dimensioned (NHOR,NLEV) so the loop below cannot
!     reach NLEP, and without the line before it the reflected part of the
!     downward longwave was propagated up through every atmospheric level and
!     never debited from the surface. `dlwfl(:,NLEP)` is what `landmod` and
!     `seamod` settle the surface energy budget with and what output code 177
!     reports, so the surface was absorbing (1-eps)*LWdown that it had just
!     reflected, while the lowest atmospheric layer paid for it. The total
!     column conserved and the partition did not.
!
!     At eps = 1 the term is identically zero, which is why the omission was
!     invisible for land while ELWLAND was 1.0. It was NOT invisible over
!     water: ELWSEA is 0.98, so the ocean and sea-ice surface has been carrying
!     a spurious 0.02*LWdown of absorbed longwave. The transmissivity from the
!     surface to itself is one, so the surface term is the whole reflection.
!
      zeps(:)=(1.-zeps(:))*dftd(:,NLEP)
      dftu(:,NLEP)=dftu(:,NLEP)-zeps(:)
      dftue2(:,NLEP)=dftue2(:,NLEP)-zeps(:)
      do jlev=1,NLEV
       dftu(:,jlev)=dftu(:,jlev)                                        &
     &             -ztausf(:,jlev)*zeps(:)
       dftue2(:,jlev)=dftue2(:,jlev)                                    &
     &             -ztausf(:,jlev)*zeps(:)
      enddo
!
!     total longwave radiation
!
      dlwfl(:,:)=dftu(:,:)+dftd(:,:)
!
      return
      end subroutine lwr

!     **********************
!     Generic Orbit Routines
!     **********************

!     =====================
!     SUBROUTINE GEN_ORB_DECL
!     =====================

!     Given a mean anomaly, eccentricity, obliquity, and moving longitude of vernal equinox,
!     compute the declination and true anomaly. This uses a Newton-Raphson iterator and
!     will work with reasonable accuracy for any bound orbit (eccen<1.0).

      subroutine gen_orb_decl(yearfraction, eccen, obliqr, mvelpp, trueanomaly, lamb, rasc, zdecl, eccf)
      use radmod, only : TWOPI, meananom0r, nfixed
      !Inputs
      real :: eccen        ! Eccentricity
      real :: yearfraction ! Elapsed fraction of the year 
      real :: obliqr       ! Obliquity in radians
      real :: mvelpp       ! Earth's moving vernal equinox longitude
!                          ! of perihelion plus pi (radians)
      !Internal
      real :: meananomaly
      real :: eccenanomaly ! Eccentric anomaly
      real thyng
      real anomarg
      real invrho
      !Outputs
      real :: trueanomaly  ! True anomaly in radians
      real :: zdecl        ! Solar declination in radians
      real :: eccf         ! Eccentricity factor for insolation
      real :: lamb         ! True anomaly - longitude of vernal equinox
      real :: rasc         ! Right ascension
      
      if (nfixed > 0) then
          trueanomaly = 0.
          eccf = 1.
      else
!           write(6,*) yearfraction
!           write(6,*) TWOPI
!           write(6,*) meananom0r
          meananomaly = yearfraction*TWOPI + meananom0r
!           write(6,*) meananomaly
          do while (meananomaly > TWOPI)
              meananomaly = meananomaly - TWOPI
          enddo
          do while (meananomaly < 0)
              meananomaly = meananomaly + TWOPI
          enddo
!           write(6,*) meananomaly
          
          if (eccen > 0.) then
              call newtonraphson(meananomaly,eccen,eccenanomaly)
              
              thyng = tan(eccenanomaly*0.5)
              anomarg = sqrt((1+eccen)/(1-eccen) * thyng*thyng)
              
              if (thyng .lt. 0.) trueanomaly = 2*atan(0.0-anomarg)
              if (thyng .ge. 0.) trueanomaly = 2*atan(anomarg)
              
              if (trueanomaly .lt. 0) trueanomaly = trueanomaly + TWOPI
              
              trueanomaly = mod(trueanomaly,TWOPI)
              
              invrho = 1./(1 - eccen*cos(eccenanomaly))
              eccf = invrho*invrho
          else  !For a circular orbit we don't need to do all that calculation
              trueanomaly = meananomaly
              eccf = 1.
          endif
      endif
      lamb = MOD(mvelpp+trueanomaly, TWOPI)
      zdecl  = asin(sin(obliqr)*sin(lamb))
      rasc = atan2(cos(obliqr)*sin(lamb),cos(lamb))
      
      return
      end subroutine gen_orb_decl
      
      
!     ========================
!     SUBROUTINE NEWTONRAPHSON
!     ========================

      subroutine newtonraphson(meananom,eccen,ee)
      use radmod, only: PI
      
      real meananom
      real eccen
      real ee
      real e0
      real thyng
      integer ict
      logical thresh
      
      if (eccen .lt. 0.5) then
        ee = meananom
      else
        ee = PI !prevents crazy excursions due to divide-by-zero
      endif
      
!       write(6,*) ee
      
      ict = 0
      thresh = .false.
      
      do while (thresh .neqv. .true.)
        e0 = ee
!         write(6,*) ee
        thyng = 1-eccen*cos(ee)
        if (thyng .lt. 1.0e-15) thyng=1.0e-15
        ee = ee - (ee-(meananom+eccen*sin(ee)))/thyng
        if (abs(ee-e0) .le. 1.0e-14) thresh = .true.
        ict = ict + 1
        if (ict .gt. 100.0) thresh = .true.
      enddo 
      
      return
      end subroutine newtonraphson
      
      
      
!     ====================
!     Earth Orbit Routines (form CCM3)
!     ====================
!
!     Routines for calculation of oribital parameters and solar
!     declination angle.
!
!     Based on f77 routines in NCAR CCM3:
!
!     Subroutines contained
!
!       orb_params --- Calculate the orbital parameters for a given
!                               situation/year.
!       orb_decl ----- Calculate the solar declination angle and
!                               Earth/Sun distance factor for a given
!                               time of the year.
!                               to use.
!
!     Code history
!
!         Original version:  Erik Kluzek
!         Date:              Oct/1997
!
!
!     Version information:
!
!         CVS: $Id: orb.F,v 1.8.8.1 1998/12/02 17:29:07 erik Exp $
!         CVS: $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/orb.F,v $
!         CVS: $Name: ccm3_6_16_brnchT_amip2_9 $
!

       module orbconst
!
!         parameters for orbital calculations
!
          implicit none
          real, parameter :: ORB_ECCEN_MIN  =   0.0                ! minimum value for eccen
          real, parameter :: ORB_ECCEN_MAX  =   0.1                ! maximum value for eccen
          real, parameter :: ORB_OBLIQ_MIN  = -180.0                ! minimum value for obliq
          real, parameter :: ORB_OBLIQ_MAX  = +180.0                ! maximum value for obliq
          real, parameter :: ORB_MVELP_MIN  =   0.0                ! minimum value for mvelp
          real, parameter :: ORB_MVELP_MAX  = 360.0                ! maximum value for mvelp
          real, parameter :: ORB_UNDEF_REAL = 1.e36                ! undefined/unset/invalid value
          real, parameter :: ORB_DEFAULT    = ORB_UNDEF_REAL       ! flag to use default orbit
          integer, parameter :: ORB_UNDEF_INT  = 2000000000        ! undefined/unset/invalid value
          integer, parameter :: ORB_NOT_YEAR_BASED = ORB_UNDEF_INT ! flag to not use input year
       end module orbconst

!     =====================
!     SUBROUTINE ORB_PARAMS
!     =====================

      subroutine orb_params(iyear_AD, eccen, obliq, meananom0, mvelp,   &
     &                      obliqr, meananom0r, lambm0, mvelpp, log_print,    &
     &                      ngenkeplerian, mypid, nroot,nud)
!
!      Calculate earth's orbital parameters using Dave Threshers
!      formula which came from Berger, Andre.  1978
!      "A Simple Algorithm to Compute Long-Term Variations
!      of Daily Insolation".  Contribution 18, Institute of Astronomy and
!      Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve,
!      Belgium.
!
!      Original Author: Erik Kluzek
!      Date:            Oct/97
!
      use orbconst
      implicit none

!     Input Arguments
!     ---------------
      real :: eccen        ! Earth's orbital eccentricity
      real :: obliq        ! Earth's obliquity in degree's
      real :: meananom0    ! Initial mean anomaly
      real :: mvelp        ! Earth's moving vernal equinox longitude
      integer :: iyear_AD  ! Year to calculate orbit for..
      logical :: log_print ! Flag to print-out status information or not.
                           ! (This turns off ALL status printing including)
                           ! (error messages.)
      integer :: mypid     ! process id (PUMA MPI)
      integer :: nroot     ! process id of root (PUMA MPI)
      integer :: nud       ! write unit for diagnostic messages
      integer :: ngenkeplerian

!     Output Arguments
!     ----------------
      real :: obliqr  ! Earth's obliquity in radians
      real :: meananom0r ! Initial mean anomaly in radians
      real :: lambm0  ! Mean longitude of perihelion at the
!                     ! vernal equinox (radians)
      real :: mvelpp  ! Earth's moving vernal equinox longitude
!                     ! of perihelion plus pi (radians)
!
! Parameters for calculating earth's orbital characteristics
! ----------
      integer, parameter :: poblen = 47  ! number of elements in the series to calc obliquity
      integer, parameter :: pecclen = 19 ! number of elements in the series to calc eccentricity
      integer, parameter :: pmvelen = 78 ! number of elements in the series to calc vernal equinox
      real :: degrad          ! degrees to radians conversion factor
      real :: obamp(poblen)   ! amplitudes for obliquity cosine series
      real :: obrate(poblen)  ! rates for obliquity cosine series
      real :: obphas(poblen)  ! phases for obliquity cosine series
      real :: ecamp(pecclen)  ! amplitudes for eccentricity/fvelp cosine/sine series
      real :: ecrate(pecclen) ! rates for eccentricity/fvelp cosine/sine series
      real :: ecphas(pecclen) ! phases for eccentricity/fvelp cosine/sine series
      real :: mvamp(pmvelen)  ! amplitudes for mvelp sine series
      real :: mvrate(pmvelen) ! rates for mvelp sine series
      real :: mvphas(pmvelen) ! phases for mvelp sine series
      real :: yb4_1950AD      ! number of years before 1950 AD
!
      real, parameter :: psecdeg = 1./3600. ! arc seconds to degrees conversion
!
!  Cosine series data for computation of obliquity:
!  amplitude (arc seconds), rate (arc seconds/year), phase (degrees).
!
      data obamp /-2462.2214466D0, -857.3232075D0, -629.3231835D0,      &
     &             -414.2804924D0, -311.7632587D0,  308.9408604D0,      &
     &             -162.5533601D0, -116.1077911D0,  101.1189923D0,      &
     &              -67.6856209D0,   24.9079067D0,   22.5811241D0,      &
     &              -21.1648355D0,  -15.6549876D0,   15.3936813D0,      &
     &               14.6660938D0,  -11.7273029D0,   10.2742696D0,      &
     &                6.4914588D0,    5.8539148D0,   -5.4872205D0,      &
     &               -5.4290191D0,    5.1609570D0,    5.0786314D0,      &
     &               -4.0735782D0,    3.7227167D0,    3.3971932D0,      &
     &               -2.8347004D0,   -2.6550721D0,   -2.5717867D0,      &
     &               -2.4712188D0,    2.4625410D0,    2.2464112D0,      &
     &               -2.0755511D0,   -1.9713669D0,   -1.8813061D0,      &
     &               -1.8468785D0,    1.8186742D0,    1.7601888D0,      &
     &               -1.5428851D0,    1.4738838D0,   -1.4593669D0,      &
     &                1.4192259D0,   -1.1818980D0,    1.1756474D0,      &
     &               -1.1316126D0,    1.0896928D0/
!
      data obrate /31.609974D0, 32.620504D0, 24.172203D0,               &
     &             31.983787D0, 44.828336D0, 30.973257D0,               &
     &             43.668246D0, 32.246691D0, 30.599444D0,               &
     &             42.681324D0, 43.836462D0, 47.439436D0,               &
     &             63.219948D0, 64.230478D0,  1.010530D0,               &
     &              7.437771D0, 55.782177D0,  0.373813D0,               &
     &             13.218362D0, 62.583231D0, 63.593761D0,               &
     &             76.438310D0, 45.815258D0,  8.448301D0,               &
     &             56.792707D0, 49.747842D0, 12.058272D0,               &
     &             75.278220D0, 65.241008D0, 64.604291D0,               &
     &              1.647247D0,  7.811584D0, 12.207832D0,               &
     &             63.856665D0, 56.155990D0, 77.448840D0,               &
     &              6.801054D0, 62.209418D0, 20.656133D0,               &
     &             48.344406D0, 55.145460D0, 69.000539D0,               &
     &             11.071350D0, 74.291298D0, 11.047742D0,               &
     &              0.636717D0, 12.844549D0/
!
      data obphas /251.9025D0, 280.8325D0, 128.3057D0,                  &
     &             292.7252D0,  15.3747D0, 263.7951D0,                  &
     &             308.4258D0, 240.0099D0, 222.9725D0,                  &
     &             268.7809D0, 316.7998D0, 319.6024D0,                  &
     &             143.8050D0, 172.7351D0,  28.9300D0,                  &
     &             123.5968D0,  20.2082D0,  40.8226D0,                  &
     &             123.4722D0, 155.6977D0, 184.6277D0,                  &
     &             267.2772D0,  55.0196D0, 152.5268D0,                  &
     &              49.1382D0, 204.6609D0,  56.5233D0,                  &
     &             200.3284D0, 201.6651D0, 213.5577D0,                  &
     &              17.0374D0, 164.4194D0,  94.5422D0,                  &
     &             131.9124D0,  61.0309D0, 296.2073D0,                  &
     &             135.4894D0, 114.8750D0, 247.0691D0,                  &
     &             256.6114D0,  32.1008D0, 143.6804D0,                  &
     &              16.8784D0, 160.6835D0,  27.5932D0,                  &
     &             348.1074D0,  82.6496D0/
!
!  Cosine/sine series data for computation of eccentricity and
!  fixed vernal equinox longitude of perihelion (fvelp):
!  amplitude, rate (arc seconds/year), phase (degrees).
!
      data ecamp /0.01860798D0,  0.01627522D0, -0.01300660D0,           &
     &            0.00988829D0, -0.00336700D0,  0.00333077D0,           &
     &           -0.00235400D0,  0.00140015D0,  0.00100700D0,           &
     &            0.00085700D0,  0.00064990D0,  0.00059900D0,           &
     &            0.00037800D0, -0.00033700D0,  0.00027600D0,           &
     &            0.00018200D0, -0.00017400D0, -0.00012400D0,           &
     &            0.00001250D0/
!
      data ecrate /4.2072050D0,  7.3460910D0, 17.8572630D0,             &
     &            17.2205460D0, 16.8467330D0,  5.1990790D0,             &
     &            18.2310760D0, 26.2167580D0,  6.3591690D0,             &
     &            16.2100160D0,  3.0651810D0, 16.5838290D0,             &
     &            18.4939800D0,  6.1909530D0, 18.8677930D0,             &
     &            17.4255670D0,  6.1860010D0, 18.4174410D0,             &
     &             0.6678630D0/
!
      data ecphas /28.620089D0, 193.788772D0, 308.307024D0,             &
     &            320.199637D0, 279.376984D0,  87.195000D0,             &
     &            349.129677D0, 128.443387D0, 154.143880D0,             &
     &            291.269597D0, 114.860583D0, 332.092251D0,             &
     &            296.414411D0, 145.769910D0, 337.237063D0,             &
     &            152.092288D0, 126.839891D0, 210.667199D0,             &
     &             72.108838D0/
!
!  Sine series data for computation of moving vernal equinox
!  longitude of perihelion:
!  amplitude (arc seconds), rate (arc seconds/year), phase (degrees).
!
      data mvamp /7391.0225890D0, 2555.1526947D0, 2022.7629188D0,       &
     &           -1973.6517951D0, 1240.2321818D0,  953.8679112D0,       &
     &            -931.7537108D0,  872.3795383D0,  606.3544732D0,       &
     &            -496.0274038D0,  456.9608039D0,  346.9462320D0,       &
     &            -305.8412902D0,  249.6173246D0, -199.1027200D0,       &
     &             191.0560889D0, -175.2936572D0,  165.9068833D0,       &
     &             161.1285917D0,  139.7878093D0, -133.5228399D0,       &
     &             117.0673811D0,  104.6907281D0,   95.3227476D0,       &
     &              86.7824524D0,   86.0857729D0,   70.5893698D0,       &
     &             -69.9719343D0,  -62.5817473D0,   61.5450059D0,       &
     &             -57.9364011D0,   57.1899832D0,  -57.0236109D0,       &
     &             -54.2119253D0,   53.2834147D0,   52.1223575D0,       &
     &             -49.0059908D0,  -48.3118757D0,  -45.4191685D0,       &
     &             -42.2357920D0,  -34.7971099D0,   34.4623613D0,       &
     &             -33.8356643D0,   33.6689362D0,  -31.2521586D0,       &
     &             -30.8798701D0,   28.4640769D0,  -27.1960802D0,       &
     &              27.0860736D0,  -26.3437456D0,   24.7253740D0,       &
     &              24.6732126D0,   24.4272733D0,   24.0127327D0,       &
     &              21.7150294D0,  -21.5375347D0,   18.1148363D0,       &
     &             -16.9603104D0,  -16.1765215D0,   15.5567653D0,       &
     &              15.4846529D0,   15.2150632D0,   14.5047426D0,       &
     &             -14.3873316D0,   13.1351419D0,   12.8776311D0,       &
     &              11.9867234D0,   11.9385578D0,   11.7030822D0,       &
     &              11.6018181D0,  -11.2617293D0,  -10.4664199D0,       &
     &              10.4333970D0,  -10.2377466D0,   10.1934446D0,       &
     &             -10.1280191D0,   10.0289441D0,  -10.0034259D0/
!
      data mvrate /31.609974D0, 32.620504D0, 24.172203D0,               &
     &              0.636717D0, 31.983787D0,  3.138886D0,               &
     &             30.973257D0, 44.828336D0,  0.991874D0,               &
     &              0.373813D0, 43.668246D0, 32.246691D0,               &
     &             30.599444D0,  2.147012D0, 10.511172D0,               &
     &             42.681324D0, 13.650058D0,  0.986922D0,               &
     &              9.874455D0, 13.013341D0,  0.262904D0,               &
     &              0.004952D0,  1.142024D0, 63.219948D0,               &
     &              0.205021D0,  2.151964D0, 64.230478D0,               &
     &             43.836462D0, 47.439436D0,  1.384343D0,               &
     &              7.437771D0, 18.829299D0,  9.500642D0,               &
     &              0.431696D0,  1.160090D0, 55.782177D0,               &
     &             12.639528D0,  1.155138D0,  0.168216D0,               &
     &              1.647247D0, 10.884985D0,  5.610937D0,               &
     &             12.658184D0,  1.010530D0,  1.983748D0,               &
     &             14.023871D0,  0.560178D0,  1.273434D0,               &
     &             12.021467D0, 62.583231D0, 63.593761D0,               &
     &             76.438310D0,  4.280910D0, 13.218362D0,               &
     &             17.818769D0,  8.359495D0, 56.792707D0,               &
     &              8.448301D0,  1.978796D0,  8.863925D0,               &
     &              0.186365D0,  8.996212D0,  6.771027D0,               &
     &             45.815258D0, 12.002811D0, 75.278220D0,               &
     &             65.241008D0, 18.870667D0, 22.009553D0,               &
     &             64.604291D0, 11.498094D0,  0.578834D0,               &
     &              9.237738D0, 49.747842D0,  2.147012D0,               &
     &              1.196895D0,  2.133898D0,  0.173168D0/
!
      data mvphas /251.9025D0, 280.8325D0, 128.3057D0,                  &
     &             348.1074D0, 292.7252D0, 165.1686D0,                  &
     &             263.7951D0,  15.3747D0,  58.5749D0,                  &
     &              40.8226D0, 308.4258D0, 240.0099D0,                  &
     &             222.9725D0, 106.5937D0, 114.5182D0,                  &
     &             268.7809D0, 279.6869D0,  39.6448D0,                  &
     &             126.4108D0, 291.5795D0, 307.2848D0,                  &
     &              18.9300D0, 273.7596D0, 143.8050D0,                  &
     &             191.8927D0, 125.5237D0, 172.7351D0,                  &
     &             316.7998D0, 319.6024D0,  69.7526D0,                  &
     &             123.5968D0, 217.6432D0,  85.5882D0,                  &
     &             156.2147D0,  66.9489D0,  20.2082D0,                  &
     &             250.7568D0,  48.0188D0,   8.3739D0,                  &
     &              17.0374D0, 155.3409D0,  94.1709D0,                  &
     &             221.1120D0,  28.9300D0, 117.1498D0,                  &
     &             320.5095D0, 262.3602D0, 336.2148D0,                  &
     &             233.0046D0, 155.6977D0, 184.6277D0,                  &
     &             267.2772D0,  78.9281D0, 123.4722D0,                  &
     &             188.7132D0, 180.1364D0,  49.1382D0,                  &
     &             152.5268D0,  98.2198D0,  97.4808D0,                  &
     &             221.5376D0, 168.2438D0, 161.1199D0,                  &
     &              55.0196D0, 262.6495D0, 200.3284D0,                  &
     &             201.6651D0, 294.6547D0,  99.8233D0,                  &
     &             213.5577D0, 154.1631D0, 232.7153D0,                  &
     &             138.3034D0, 204.6609D0, 106.5938D0,                  &
     &             250.4676D0, 332.3345D0,  27.3039D0/
!
!     Local variables
!     ---------------
      integer i        ! Index for series summations
      real :: obsum    ! Obliquity series summation
      real :: cossum   ! Cosine series summation for eccentricity/fvelp
      real :: sinsum   ! Sine series summation for eccentricity/fvelp
      real :: fvelp    ! Fixed vernal equinox longitude of perihelion
      real :: mvsum    ! mvelp series summation
      real :: beta     ! Intermediate argument for lambm0
      real :: years    ! Years to time of interest (negative = past;
!                      ! positive = future)
      real :: eccen2   ! eccentricity squared
      real :: eccen3   ! eccentricity cubed
      real :: pi       ! pi
!
! radinp and algorithms below will need a degrees to radians conversion
! factor.
!
      pi     =  4.*atan(1.)
      degrad = pi/180.
!
! Check for flag to use input orbit parameters
!
      if ( iyear_AD .eq. ORB_NOT_YEAR_BASED ) then
!
! Check input obliq, eccen, and mvelp to ensure reasonable
!
         if( obliq .eq. ORB_UNDEF_REAL )then
          if ( log_print ) then
           if(mypid==nroot) then
            write(nud,*)'(orb_params) Have to specify orbital parameters:'
            write(nud,*) 'Either set: '                                   &
     &                ,'iyear_AD, OR [obliq, eccen, and mvelp]:'
            write(nud,*)'iyear_AD is the year to simulate the orbit for ' &
     &                ,'(ie. 1950): '
            write(nud,*)'obliq, eccen, mvelp specify the orbit directly:'
            write(nud,*)'The AMIP II settings (for a 1995 orbit) are: '
            write(nud,*)' obliq = 23.4441'
            write(nud,*)' eccen = 0.016715'
            write(nud,*)' mvelp = 102.7'
           end if
          end if
          stop 999
        else if ( log_print ) then
          if(mypid==nroot) then
           write(nud,*)'(orb_params) Use input orbital parameters: '
          end if
         end if
         if( (obliq.lt.ORB_OBLIQ_MIN).or.(obliq.gt.ORB_OBLIQ_MAX) ) then
          if ( log_print ) then
           if(mypid==nroot) then
             write(nud,*) '(orb_params): Input obliquity unreasonable: '  &
     &                  ,obliq
           end if
          end if
          stop 999
         end if
         if( ((eccen.lt.ORB_ECCEN_MIN).or.(eccen.gt.ORB_ECCEN_MAX)).and.(ngenkeplerian==0) ) then
          if ( log_print ) then
           if(mypid==nroot) then
            write(nud,*) '(orb_params): Input eccentricity unreasonable: '&
     &                 ,eccen
           end if
          end if
          stop 999
         else if ((ngenkeplerian==1) .and. ((eccen.lt.0).or.(eccen.ge.1))) then
          if ( log_print ) then
           if(mypid==nroot) then
            write(nud,*) '(orb_params): Input eccentricity unreasonable: '&
     &                 ,eccen
           end if
          end if
          stop 999
         end if
         if( (mvelp.lt.ORB_MVELP_MIN).or.(mvelp.gt.ORB_MVELP_MAX) ) then
          if ( log_print ) then
           if(mypid==nroot) then
             write(nud,*)'(orb_params): Input mvelp unreasonable: ', mvelp
           endif
          end if
          stop 999
         end if
        eccen2 = eccen*eccen
        eccen3 = eccen2*eccen
      else
!
! Otherwise calculate based on years before present
!
        yb4_1950AD = 1950.0 - float(iyear_AD)
        if ( abs(yb4_1950AD) .gt. 1000000.0 )then
          if ( log_print ) then
           if(mypid==nroot) then
            write(nud,*)'(orb_params) orbit only valid for years+-1000000'
            write(nud,*)'(orb_params) Relative to 1950 AD'
            write(nud,*)'(orb_params) # of years before 1950: ',yb4_1950AD
            write(nud,*)'(orb_params) Year to simulate was  : ',iyear_AD
           end if
          end if
          stop 999
        end if
!
!
! The following calculates the earth's obliquity, orbital eccentricity
! (and various powers of it) and vernal equinox mean longitude of
! perihelion for years in the past (future = negative of years past),
! using constants (see parameter section) given in the program of:
!
! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
! of Daily Insolation.  Contribution 18, Institute of Astronomy and
! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
!
! and formulas given in the paper (where less precise constants are also
! given):
!
! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
!
! The algorithm is valid only to 1,000,000 years past or hence.
! For a solution valid to 5-10 million years past see the above author.
! Algorithm below is better for years closer to present than is the
! 5-10 million year solution.
!
! Years to time of interest must be negative of years before present
! (1950) in formulas that follow.
!
        years = - yb4_1950AD
!
! In the summations below, cosine or sine arguments, which end up in
! degrees, must be converted to radians via multiplication by degrad.
!
! Summation of cosine series for obliquity (epsilon in Berger 1978) in
! degrees. Convert the amplitudes and rates, which are in arc seconds, into
! degrees via multiplication by psecdeg (arc seconds to degrees conversion
! factor).  For obliq, first term is Berger 1978's epsilon star; second
! term is series summation in degrees.
!
        obsum = 0.0
        do i = 1, poblen
          obsum = obsum +                                               &
     &           obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years +        &
     &                                obphas(i))*degrad)
        end do
        obliq = 23.320556 + obsum
!
! Summation of cosine and sine series for computation of eccentricity
! (eccen; e in Berger 1978) and fixed vernal equinox longitude of perihelion
! (fvelp; pi in Berger 1978), which is used for computation of moving vernal
! equinox longitude of perihelion.  Convert the rates, which are in arc
! seconds, into degrees via multiplication by psecdeg.
!
        cossum = 0.0
        do i = 1, pecclen
          cossum = cossum +                                             &
     &            ecamp(i)*cos((ecrate(i)*psecdeg*years +               &
     &                          ecphas(i))*degrad)
        end do
!
        sinsum = 0.0
        do i = 1, pecclen
          sinsum = sinsum +                                             &
     &            ecamp(i)*sin((ecrate(i)*psecdeg*years +               &
     &                          ecphas(i))*degrad)
        end do
!
! Use summations to calculate eccentricity
!
        eccen2 = cossum*cossum + sinsum*sinsum
        eccen = sqrt(eccen2)
        eccen3 = eccen2*eccen
!
! A series of cases for fvelp, which is in radians.
!
        if (abs(cossum) .le. 1.0E-8) then
          if (sinsum .eq. 0.0) then
            fvelp = 0.0
          else if (sinsum .lt. 0.0) then
            fvelp = 1.5*pi
          else if (sinsum .gt. 0.0) then
            fvelp = .5*pi
          endif
        else if (cossum .lt. 0.0) then
          fvelp = atan(sinsum/cossum) + pi
        else if (cossum .gt. 0.0) then
          if (sinsum .lt. 0.0) then
            fvelp = atan(sinsum/cossum) + 2.0*pi
          else
            fvelp = atan(sinsum/cossum)
          endif
        endif
!
! Summation of sine series for computation of moving vernal equinox longitude
! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
! first term is fvelp in degrees; second term is Berger 1978's psi bar times
! years and in degrees; third term is Berger 1978's zeta; fourth term is
! series summation in degrees.  Convert the amplitudes and rates, which are
! in arc seconds, into degrees via multiplication by psecdeg.  Series summation
! plus second and third terms constitute Berger 1978's psi, which is the
! general precession.
!
        mvsum = 0.0
        do i = 1, pmvelen
          mvsum = mvsum +                                               &
     &           mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years +        &
     &                                mvphas(i))*degrad)
        end do
        mvelp = fvelp/degrad + 50.439273*psecdeg*years + 3.392506       &
     &  + mvsum
!
! Cases to make sure mvelp is between 0 and 360.
!
        do while (mvelp .lt. 0.0)
          mvelp = mvelp + 360.0
        end do
        do while (mvelp .ge. 360.0)
          mvelp = mvelp - 360.0
        end do
      end if  ! end of test on whether to calculate or use input orbital params
!
! Orbit needs the obliquity in radians
!
      obliqr = obliq*degrad
      meananom0r = meananom0*degrad
!
! 180 degrees must be added to mvelp since observations are made from the
! earth and the sun is considered (wrongly for the algorithm) to go around
! the earth. For a more graphic explanation see Appendix B in:
!
! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth's Orbital
! Periods.  J. of Geophysical Research 98:10,341-10,362.
!
! Additionally, orbit will need this value in radians. So mvelp becomes
! mvelpp (mvelp plus pi)
!
      mvelpp = (mvelp + 180.)*degrad
!
! Set up an argument used several times in lambm0 calculation ahead.
!
      beta = sqrt(1. - eccen2)
!
! The mean longitude at the vernal equinox (lambda m nought in Berger
! 1978; in radians) is calculated from the following formula given in
! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
! 1978) is 0.
!
      lambm0 = 2.*((.5*eccen + .125*eccen3)*(1. + beta)*sin(mvelpp)     &
     &            - .25*eccen2*(.5 + beta)*sin(2.*mvelpp)               &
     &            + .125*eccen3*(1./3. + beta)*sin(3.*mvelpp))
!
      if ( log_print ) then
       if(mypid==nroot) then
        write(nud,'(/," *****************************************")')
        write(nud,'(" *     Computed Orbital Parameters       *")')
        write(nud,'(" *****************************************")')
        write(nud,'(" * Year AD           =  ",i16  ," *")') iyear_AD
        write(nud,'(" * Eccentricity      =  ",f16.6," *")') eccen
        write(nud,'(" * Obliquity (deg)   =  ",f16.6," *")') obliq
        write(nud,'(" * Obliquity (rad)   =  ",f16.6," *")') obliqr
        write(nud,'(" * Long of perh(deg) =  ",f16.6," *")') mvelp
        write(nud,'(" * Long of perh(rad) =  ",f16.6," *")') mvelpp
        write(nud,'(" * Long at v.e.(rad) =  ",f16.6," *")') lambm0
        write(nud,'(" *****************************************")')
       end if
      end if
!
!
      return
      end subroutine orb_params


!     ===================
!     SUBROUTINE ORB_DECL
!     ===================

      subroutine orb_decl(calday,eccen,mvelpp,lambm0,obliqr,tnu,lamb,rasc,delta,eccf)
      use pumamod, only: mcal_days_per_year,ndatim
!
!     Compute earth/orbit parameters using formula suggested by
!     Duane Thresher.
!
!     Original version:  Erik Kluzek
!     Date:              Oct/1997
!
!     Modification: 22-Feb-2006 (ek) - get days/yr from pumamod
!
      implicit none
!
!     Input arguments
!     ---------------
      real :: calday     ! Calendar day, including fraction
      real :: eccen      ! Eccentricity
      real :: obliqr     ! Earth's obliquity in radians
      real :: lambm0     ! Mean longitude of perihelion at the
!                        ! vernal equinox (radians)
      real :: mvelpp     ! Earth's moving vernal equinox longitude
!                        ! of perihelion plus pi (radians)
!
!     Output arguments
!     ----------------
      real :: delta      ! Solar declination angle in radians
      real :: eccf       ! Earth-sun distance factor ( i.e. (1/r)**2 )
      real :: lamb    ! Lambda, the earth's longitude of perihelion
      real :: tnu     ! Earth's true anomaly in radians
!
!     Local variables
!     ---------------
      real, parameter :: ve = 80.5 ! Calday of vernal equinox
!                                  ! correct for Jan 1 = calday 1
      real, parameter :: pie = 3.141592653589793D0
!
      real :: lambm   ! Lambda m, earth's mean longitude of perihelion (radians)
      real :: lmm     ! Intermediate argument involving lambm
      real :: invrho  ! Inverse normalized sun/earth distance
      real :: sinl    ! Sine of lmm
      real :: rasc         ! Right ascension
!
! Compute eccentricity factor and solar declination using
! day value where a round day (such as 213.0) refers to 0z at
! Greenwich longitude.
!
! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
! 35:2362-2367.
!
! To get the earth's true longitude (position in orbit; lambda in Berger 1978),
! which is necessary to find the eccentricity factor and declination,
! must first calculate the mean longitude (lambda m in Berger 1978) at
! the present day.  This is done by adding to lambm0 (the mean longitude
! at the vernal equinox, set as March 21 at noon, when lambda = 0; in radians)
! an increment (delta lambda m in Berger 1978) that is the number of
! days past or before (a negative increment) the vernal equinox divided by
! the days in a model year times the 2*pi radians in a complete orbit.
!

!     EARTH'S CALENDAR, HARDCODED, and it does not scale. ve is Earth calendar
!     day 80.5 of a 365-day year, so ve/365 = 0.2205 is Earth's phase from
!     1 January to its vernal equinox; calday is a FRACTION of the year in this
!     fork, so the offset is a constant that belongs to another planet's
!     calendar. Dormant on ONE KEYWORD: run_exoplasim.py passes keplerian=True,
!     which routes solang to gen_orb_decl, whose phase comes from mvelpp and
!     meananom0r and is correct and config-reachable. Setting keplerian=False
!     lands here. world-9d1.
      lambm  = lambm0 + (calday - ve/365.)*2.*pie                            !& Moving to more robust system
            ! / (mcal_days_per_year + ndatim(7)) ! ndatim(7) = leap year
      lmm    = lambm  - mvelpp
!
! The earth's true longitude, in radians, is then found from
! the formula in Berger 1978:
!
      sinl   = sin(lmm)
      lamb   = lambm  + eccen*(2.*sinl                                  &
     &         + eccen*(1.25*sin(2.*lmm)                                &
     &         + eccen*((13.0/12.0)*sin(3.*lmm) - 0.25*sinl)))
!
! Using the obliquity, eccentricity, moving vernal equinox longitude of
! perihelion (plus), and earth's true longitude, the declination (delta)
! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
! rho will be used), and thus the eccentricity factor (eccf), can be calculated
! from formulas given in Berger 1978.
!
      invrho = (1. + eccen*cos(lamb - mvelpp))                          &
     &         / (1. - eccen*eccen)
!
! Set solar declination and eccentricity factor
!
      delta  = asin(sin(obliqr)*sin(lamb))
      eccf   = invrho*invrho
      tnu = MOD(lamb - mvelpp, 2*pie)
      rasc = atan2(cos(obliqr)*sin(lamb),cos(lamb))
!
      return
      end subroutine orb_decl
