!
!     ***********************
!     * Planet Simulator 17 *
!     ***********************

!     *****************
!     * Frank Lunkeit *       University of Hamburg
!     * Edilbert Kirk *      Meteorological Institute
!     *****************   Dept.: Theoretical Meteorology
!                             Head: Klaus Fraedrich

!     **************
!     * Name rules *
!     **************

!     i - local integer
!     j - loop index
!     k - integer dummy parameter
!     N - integer constants

!     g - real gridpoint arrays
!     p - real dummy parameter
!     s - real spectral arrays
!     z - local real

      module pumamod

!     ****************************************************************
!     * The module resmod defines all resolution parameters          *
!     * NPRO_ATM, NLAT_ATM, NLEV_ATM                                 *
!     * If MPI (Message Passing Interface) is to be used it will     *
!     * also include a 'use mpi' statement                           *
!     ****************************************************************

      use resmod

!     ****************************************************************
!     * The number of processes for processing on parallel machines  *
!     * It must be a power of 2 e.g. (1, 2, 4, 8, ... NLAT)          *
!     ****************************************************************

      parameter(NPRO = NPRO_ATM)          ! Number of processes (CPUs)

!     ****************************************************************
!     * Set the horizontal resolution by specifying the number of    *
!     * latitudes. Use only values, listed below (FFT restrictions)  *
!     ****************************************************************

      parameter(NLAT = NLAT_ATM)                !  Number of latitudes

!     currently allowed values:  2,  4,  32,  48,  64, 128,  192,  256
!     corresponding truncation: T1, T2, T21, T31, T42, T85, T127, T170

!     ****************************************************************
!     * Set the vertical resolution, a minimum of 5 is recommended   *
!     ****************************************************************

      parameter(NLEV = NLEV_ATM)                    ! Number of levels

!      *********************
!      * filenames & units *
!      *********************

      integer :: nud = 6 ! plasim_diag 
      integer :: nut = 7 ! reusable temporary unit

      character (256) :: plasim_namelist     = "plasim_namelist"
      character (256) :: radmod_namelist     = "radmod_namelist"
      character (256) :: miscmod_namelist    = "miscmod_namelist"
      character (256) :: fluxmod_namelist    = "fluxmod_namelist"
      character (256) :: rainmod_namelist    = "rainmod_namelist"
      character (256) :: surfmod_namelist    = "surfmod_namelist"
      character (256) :: plasim_output       = "plasim_output"
      character (256) :: plasim_snapshot     = "plasim_snapshot"
      character (256) :: plasim_hcadence     = "plasim_hcadence"
      character (256) :: plasim_diag         = "plasim_diag"
      character (256) :: plasim_restart      = "plasim_restart"
      character (256) :: plasim_status       = "plasim_status"
      character (256) :: planet_namelist     = "planet_namelist"
      character (256) :: efficiency_dat      = "efficiency.dat"
      character (256) :: icemod_namelist     = "icemod_namelist"
      character (256) :: ice_output          = "ice_output"
      character (256) :: oceanmod_namelist   = "oceanmod_namelist"
      character (256) :: ocean_output        = "ocean_output"
      character (256) :: landmod_namelist    = "landmod_namelist"
      character (256) :: vegmod_namelist     = "vegmod_namelist"
      character (256) :: seamod_namelist     = "seamod_namelist"
      character (256) :: aero_namelist       = "aero_namelist"

!     ****************************************************************
!     * Don't touch the following parameter definitions !            *
!     ****************************************************************

!     ********************
!     * Global Constants *
!     ********************
!
      parameter(NTRACE = 1)                ! # of tracers 1st. reserved for q
      parameter(NAERO = 1)                 ! # of aerosols (tracers with extra gravitational settling term)
      parameter(NLON = NLAT + NLAT)        ! Number of longitudes
      parameter(NTRU = (NLON-1) / 3)       ! Triangular truncation
      parameter(NLPP = NLAT / NPRO)        ! Latitudes per process
      parameter(NLHP = NLPP / 2)           ! Half of them, one hemisphere
      parameter(NHOR = NLON * NLPP)        ! Horizontal part
      parameter(NUGP = NLON * NLAT)        ! Number of gridpoints
      parameter(NPGP = NLON * NLAT / 2)    ! Dimension of packed fields
      parameter(NLEM = NLEV - 1)           ! Levels - 1
      parameter(NLEP = NLEV + 1)           ! Levels + 1
      parameter(NLSQ = NLEV * NLEV)        ! Levels squared
      parameter(NTP1 = NTRU + 1)           ! Truncation + 1
      parameter(NRSP =(NTRU+1)*(NTRU+2))   ! No of real global    modes
      parameter(NCSP = NRSP / 2)           ! No of complex global modes
      parameter(NSPP = (NRSP+NPRO-1)/NPRO) ! Modes per process
      parameter(NESP = NSPP * NPRO)        ! Dim of spectral fields
      parameter(NVCT = 2 * (NLEV+1))       ! Dim of Vert. Coord. Tab
      parameter(NZOM  = 2 * NTP1)          ! Dim for zonal mean diagnostics
      parameter(NROOT = 0)                 ! Master node

!     ****************************************************************
!     * PAIRED LATITUDE DECOMPOSITION                                *
!     *                                                              *
!     * A latitude and its mirror carry the same Legendre magnitudes *
!     * -- P(-mu) = (-1)**(m+n) P(mu) -- so a transform that holds   *
!     * both on one process does half the multiplies and reads the   *
!     * weight matrix once instead of twice. The stock decomposition *
!     * hands each process a CONTIGUOUS block of latitudes, which    *
!     * puts a latitude and its mirror on different processes and    *
!     * makes that saving unreachable on anything but one process.   *
!     *                                                              *
!     * With LPAIRLAT the scatter is permuted instead, so local      *
!     * latitude l and NLPP+1-l are a mirror pair on every process.  *
!     * It needs NPRO to divide NLAT/2, and falls back to the stock  *
!     * contiguous layout when it does not. At NPRO == 1 the         *
!     * permutation is the identity, which is why the single-process *
!     * executable is the trivial case here and not a special one.   *
!     *                                                              *
!     * The permutation lives entirely in the five grid-space        *
!     * transfers of mpimod; every latitude-dependent quantity is    *
!     * derived from the scattered sid/gwd/csq/rcs and follows it.   *
!     ****************************************************************

!     THE PAIRED DECOMPOSITION IS RETIRED. It existed to let legmod fold a
!     latitude and its mirror together, and SHTns replaces legmod -- and cannot
!     use the permuted layout anyway, since it needs the grid in latitude order.
!     Carrying two grid layouts through the rest of that conversion is a cost
!     this project has paid once. The symbol stays only so the paired branches
!     can be deleted in a change that is provably inert on its own.
      logical,parameter :: LPAIRLAT = .false.

      parameter(EZ     = 1.63299310207D0)  ! ez = 1 / sqrt(3/8)
      parameter(PI     = 3.14159265359D0)  ! Pi
      parameter(TWOPI  = PI + PI)          ! 2 Pi
      parameter(RV     = 461.51)           ! Gas constant for water vapour
      parameter(ACPV   = 1870.)            ! Specific heat for water vapour
      parameter(TMELT_CO2 = 148.0)         ! Melting point (CO2) - for Mars

!     ***************
!     * Date & Time *
!     ***************

      integer :: nstep           =       0 ! current timestep
      integer :: nhcstp            =     1 ! Timestep relative to run start
      integer :: nstep1          =       0 ! start timestep for this run
      integer :: mstep           =       0 ! timestep # in current month
      integer :: mocd            =       0 ! month countdown
      integer :: n_start_year    =    0001 ! start year
      integer :: n_start_month   =      01 ! start month
      integer :: n_start_step    =       0 ! start step since 1-Jan-0000
      integer :: n_days_per_month=      30 ! needed for time interpolation
      integer :: m_days_per_month=      30 ! standard days
      integer :: n_days_per_year =     360 ! set to 365 for real calendar
      integer :: m_days_per_year =     360 ! standard days
      integer :: mcal_days_per_year=   360 ! Days in a year in calmod (set by calini)
      integer :: n_steps_per_year=    11520! Number of timesteps per year
      integer :: n_run_years     =       0 ! years to run
      integer :: n_run_months    =       0 ! months to run
      integer :: n_run_days      =       0 ! days  to run (debugging)
      integer :: n_run_steps     =       0 ! steps to run (debugging)
      real :: mpstep             =       0.0 ! minutes/timestep = 1day/ntspd
      integer :: ntspd           =       0 ! number of timesteps per day
      integer :: mtspd           =       0 ! number of timesteps per standard day
      integer :: nwpd            =       1 ! number of writes per day
      integer :: nlowio          =       1 ! Low I/O mode (0/1)
      integer :: nstpw           =       0 ! Timesteps between writes (0=use nwpd)
      integer :: nstps           =       0 ! Steps per snapshot (0=use ntspd)
      integer :: ndatim(7)       =      -1 ! date & time array
      real    :: tmstart         =     0.0 ! start of run

      

!     **************************
!     * Global Integer Scalars *
!     **************************

      integer :: kick     =  1  ! add noise for kick > 0
      integer :: mars     =  0  ! global switch for planet mars
      integer :: noutput  =  1  ! master switch for output: 0=no output
      integer :: nsnapshot = 0  ! switch for snapshot output
      integer :: nhcadence = 0  ! Switch for high-cadence snapshot output
      integer :: hcstartstep = -1 ! Timestep to start high-cadence output
      integer :: hcendstep = -1 ! Timestep on which to end high-cadence output (exclusive)
      integer :: hcinterval = 1 ! Number of timesteps per high-cadence output
      integer :: nafter   =  0  ! write data interval: 0 = once per day
      integer :: naqua    =  0  ! 1: switch to aqua planet mode
      integer :: ndesert  =  0  ! 1: switch to desert planet mode
      integer :: nveg     =  0  ! 0: off; 1: diagnostic vegetation; 2: coupled vegetation
      integer :: ncoeff   =  0  ! number of modes to print
      integer :: ndiag    =  0  ! write diagnostics interval 0 = every 10th. day
      integer :: ngui     =  0  ! 1: run with GUI
      integer :: sellon   =  1  ! index of longitude for column mode
      integer :: nkits    =  3  ! number of initial timesteps
      integer :: nrestart =  0  ! 1 for true, 0 for false
      integer :: nrad     =  1  ! switches radiation off/on  1/0
      integer :: nflux    =  1  ! vertical diffusion 1/0
      integer :: nadv     =  1  ! advection 1/0=(y/n)
      integer :: nhordif  =  1  ! horizontal diffusion 1/0=(y/n)
      integer :: neqsig   =  0  ! equidistant sigma levels (1/0)=(y/n) !2=log-equidistant; 3=pseudolog; 4=lin-equidistant
      integer :: nprint   =  0  ! comprehensive print out (only for checks!)
      integer :: nprhor   =  0  ! grid point for print out (only for checks!)
      integer :: npacksp  =  0  ! pack spectral fields on output
      integer :: npackgp  =  0  ! pack gridpoint fields on output
      integer :: naccuout =  0  ! accumulation counter for diagnistics
      integer :: ndiaggp  =  0  ! switch for frank's gp-diagnostic arrays
      integer :: ndiagsp  =  0  ! switch for frank's sp-diagnostic arrays
      integer :: ndiagcf  =  0  ! switch for cloud forcing diagnostic
      integer :: ndiaggp2d=  0  ! number of additional 2-d gp-diagnostic arrays
      integer :: ndiaggp3d=  0  ! number of additional 3-d gp-diagnostic arrays
      integer :: ndiagsp2d=  0  ! number of additional 2-d sp-diagnostic arrays
      integer :: ndiagsp3d=  0  ! number of additional 3-d sp-diagnostic arrays
      integer :: ndivdamp =  0  ! divergence damping countdown
      integer :: nhdiff   = 15  ! critical wavenumber for horizontal diffusion
      integer :: ntime    =  0  ! switch for time use diagnostics
      integer :: nperpetual = 0 ! radiation day for perpetual integration
      integer :: n_sea_points=0 ! number of sea points on grid
      integer :: nentropy = 0   ! switch for entropy diagnostics
      integer :: nentro3d = 0   ! switch for 3d entropy diagnostics
      integer :: nenergy  = 0   ! switch for energy diagnostics
      integer :: nener3d  = 0   ! switch for 3d energy diagnostics
      integer :: ndheat   = 1   ! switch for heating due to momentum dissipation
      integer :: nseedlen = 0   ! length of random seed (set by lib call)
      integer :: nsela    = 1   ! enable (1) or disable (0) Semi Lagrangian Advection
      integer :: nspinit  = 0   ! switch for LnPs initialization
      integer :: nsponge  = 0   ! switch for top sponge layer
      integer :: nstratosponge = 0 ! Switch for Newtonian cooling in hybrid stratosphere
      integer :: nqspec   = 1   ! 1: spectral q   0: gridpoint q (semi-Langrangian)
!     SHTns is opt-in while it proves itself, so legmod stays the default and
!     the model is unchanged at 0. It needs every latitude in one address space
!     and the grid in latitude order, so it is refused on anything but the
!     threaded build without LPAIRLAT -- shtns_setup checks rather than assumes.
!     SHTns needs every latitude in one address space, so it is the default
!     where that holds and unavailable where it does not. legmod stays as the
!     reference verify_shtns_model.sh compares against, reachable with NSHTNS=0.
#ifdef OMPSHARED
      integer :: nshtns   = 1   ! 1: SHTns transforms   0: legmod's own
#else
      integer :: nshtns   = 0   ! the MPI build cannot use SHTns
#endif
      integer :: nrdrag   = 0   ! 1: Apply Rayleigh fraction to 20-layer atmosphere
      integer :: l_aero    = 1   ! 1: Aerosols on; this also enables the semi-Lagrangian advection tracer grid initialisations
!>>> AYP -- NEEDED AS PART OF GLACIERMOD      
      integer :: nglspec = 0
!>>> AYP      
      integer :: nhurricane = 0 ! 1/0=On/Off switch for hurricane metrics and monitoring

!     ***********************
!     * Global Real Scalars *
!     ***********************

      real :: als   = 2.8345E6! Latent heat of sublimation
      real :: alv   = 2.5008E6! Latent heat of vaporization
      real :: plavor=     EZ  ! planetary vorticity
      real :: dawn  =     0.0 ! angle threshold for solar radiation
      real :: deltsec  =  0.0 ! timestep [sec]
      real :: deltsec2 =  0.0 ! timestep [sec] * 2
      real :: delt            ! deltsec * Omega (ww)
      real :: delt2           ! 2 * delt
      real :: dtep  =     0.0
      real :: dtns  =     0.0
      real :: dttl  =     0.0  ! Tidally-locked substellar-antistellar temperature diff. [K]
      real :: dtrop = 12000.0
      real :: dttrp =     2.0
      real :: tgr   =   288.0  ! Temperature ground in mean profile
      real :: psurf =101100.0  ! global mean surface pressure
      real :: ptop  =  7500.0  ! Upper pressure to use to anchor upper levels (Pa)
                               ! (actual top pressure is ~0.5 this value)
      real :: ptop2 =     1.0  ! TOA pressure (Pa) for hybrid scheme with stratosphere                         
      real :: pfac  =       1  ! pN2 / 1.0e5
      real :: fixedlon = 0.0   ! Longitude of fixed solar zenith
      real :: time0 =     0.0  ! start time (for performance estimates)
      real :: co2   =   360.0  ! atm. co2 concentration (ppmv)
      real :: umax  =     0.0  ! diagnostic U max
      real :: t2mean=     0.0  ! diagnostic T2m mean
      real :: tmelt = 273.16   ! Melting point (H2O)
      real :: precip=     0.0  ! diagnostic precipitation mean
      real :: evap  =     0.0  ! diagnostic evaporation
      real :: olr   =     0.0  ! Outgoing longwave radiation
      real :: dampsp=     0.0  ! damping time (days) for sponge layer
      real :: taucool=   10.0  ! Cooling timescale for stratosphere (days)
      real :: frcmod=     1.0  ! modifier for rayleigh drag timescale
      
      real :: gpimax   =  0.0  ! Maximum Genesis Potential Index
      real :: ventimin =  0.0  ! Maximum Ventilation Index
      real :: laavmax  =  0.0  ! Maximum Lower Atmospheric Absolute Vorticity
      real :: mpotimax =  0.0  ! Global maximum of max potential intensity
      real :: vrmpimax =  0.0  ! Global maximum of ventilation-reduced max potential intensity
      integer :: nwritehurricane = 0 ! Do we write hurricane output on this timestep?

!     **************************
!     * Global Spectral Arrays *
!     **************************

      real, target ::  sd(NESP,NLEV) = 0.0 ! Spectral Divergence
      real, target ::  st(NESP,NLEV) = 0.0 ! Spectral Temperature
      real, target ::  sz(NESP,NLEV) = 0.0 ! Spectral Vorticity
      real, target ::  sq(NESP,NLEV) = 0.0 ! Spectral Specific Humidity
      real, target ::  sp(NESP)      = 0.0 ! Spectral Pressure (ln Ps)
      real ::  so(NESP)      = 0.0 ! NOT shared: fc2sp writes all of it, per thread ! Spectral Orography
      real, target ::  sr(NESP,NLEV) = 0.0 ! Spectral Restoration Temperature
      
      real :: sdipolep(NSPP) = 0.0 ! Spectral tidally-locked temperature dipole
      real :: sdipole(NESP) = 0.0 ! Spectral tidally-locked temperature dipole

#ifdef OMPSHARED
!     THE PARTIALS ARE NOT STORAGE. Each is this thread's slice of the full
!     array above, so writing sdp IS writing sd and the gather that used to
!     assemble sd out of every thread's sdp becomes a barrier that copies
!     nothing. Associated once a thread by assoc_spectral, called from mpstart.
!
!     Shared build only. Under MPI a rank's partial has to be its own array:
!     mpi_allgather forbids a send buffer aliasing its receive buffer, and the
!     MPI build is the reference this one is verified against.
      real, pointer :: sdp(:,:) => NULL() ! Spectral Divergence  Partial
      real, pointer :: stp(:,:) => NULL() ! Spectral Temperature Partial
      real, pointer :: szp(:,:) => NULL() ! Spectral Vorticity   Partial
      real, pointer :: sqp(:,:) => NULL() ! Spectral S.Humidity  Partial
      real, pointer :: spp(:)   => NULL() ! Spectral Pressure    Partial
      real, pointer :: srp(:,:) => NULL() ! Spectral Restoration Partial
#else
      real :: sdp(NSPP,NLEV) = 0.0 ! Spectral Divergence  Partial
      real :: stp(NSPP,NLEV) = 0.0 ! Spectral Temperature Partial
      real :: szp(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Partial
      real :: sqp(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Partial
      real :: spp(NSPP)      = 0.0 ! Spectral Pressure    Partial
      real :: srp(NSPP,NLEV) = 0.0 ! Spectral Restoration Partial
#endif
      real :: sop(NSPP)      = 0.0 ! Spectral Orography   Partial, NOT a slice

      real :: sdt(NSPP,NLEV) = 0.0 ! Spectral Divergence  Tendency
      real :: stt(NSPP,NLEV) = 0.0 ! Spectral Temperature Tendency
      real :: szt(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Tendency
      real :: sqt(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Tendency
      real :: spt(NSPP)      = 0.0 ! Spectral Pressure    Tendency

!     THE TENDENCY PARTIALS, one slot per process, and NOT storage private to
!     the routine that fills them. mktend, qtend and fc2sp write a process's
!     contribution straight into the array the reduction reads, so mpsumscp
!     sums where the numbers already are instead of staging every partial
!     through a buffer first. That staging was 68 GB over a 300-step T127 run.
!
!     ONE SLOT UNDER MPI, where a rank's address space already separates it
!     from every other rank's, so the slot IS the rank's own partial and
!     mpsumscp is the reduce-scatter it always was. NPRO slots under threads,
!     which share an address space and would otherwise overwrite each other.
!     Sequence association makes the two cases the same actual argument.
!
!     Four are live at once in gridpointa -- mktend fills three and fc2sp the
!     fourth -- so they cannot share one buffer between them the way the
!     collectives' scratch does. gridpointa and gridpointd do not overlap, so
!     they do share these.
#ifdef OMPSHARED
      integer, parameter :: NPART = NPRO
#else
      integer, parameter :: NPART = 1
#endif
      real :: zpsd(NESP,NLEV,0:NPART-1) = 0.0 ! Divergence  partial, by process
      real :: zpst(NESP,NLEV,0:NPART-1) = 0.0 ! Temperature partial, by process
      real :: zpsz(NESP,NLEV,0:NPART-1) = 0.0 ! Vorticity   partial, by process
      real :: zpsq(NESP,NLEV,0:NPART-1) = 0.0 ! S.Humidity  partial, by process
      real :: zpsp(NESP,     0:NPART-1) = 0.0 ! Pressure    partial, by process

!     MKDHEAT'S SCRATCH. The frictional heating term runs every timestep and
!     held thirteen arrays on the stack, 18.7 MB a thread at T170 and 299 MB
!     across sixteen, of which six were FULL spectral arrays -- one complete
!     copy of the global field per thread. These four replace them and there is
!     one of each: a module array is shared between threads unless it is
!     threadprivate, and every rank has its own, so the same declaration says
!     the right thing on both builds.
      real :: zhd(NESP,NLEV) = 0.0 ! divergence,  gathered
      real :: zhz(NESP,NLEV) = 0.0 ! vorticity,   gathered
      real :: zhq(NESP,NLEV) = 0.0 ! humidity,    gathered
      real :: zhe(NESP,NLEV) = 0.0 ! kinetic energy loss, gathered after reduction

!     and its three reduction partials, one slot per process, exactly as the
!     tendency partials above.
      real :: zhf1(NESP,NLEV,0:NPART-1) = 0.0 ! heating from the first wind pair
      real :: zhf2(NESP,NLEV,0:NPART-1) = 0.0 ! heating from the second
      real :: zhef(NESP,NLEV,0:NPART-1) = 0.0 ! kinetic energy loss, partial
      
      real :: sdm(NSPP,NLEV) = 0.0 ! Spectral Divergence  Minus
      real :: stm(NSPP,NLEV) = 0.0 ! Spectral Temperature Minus
      real :: szm(NSPP,NLEV) = 0.0 ! Spectral Vorticity   Minus
      real :: sqm(NSPP,NLEV) = 0.0 ! Spectral S.Humidity  Minus
      real :: spm(NSPP)      = 0.0 ! Spectral Pressure    Minus
      real :: srm(NSPP,NLEV) = 0.0 ! Spectral Restoration Minus

      real :: sak(NESP,NLEV)   = 0.0 ! horizontal diffusion
      real :: sakpp(NSPP,NLEV) = 0.0 ! horizontal diffusion partial
      real :: sqout(NESP,NLEV) = 0.0 ! specific humidity for output
      real :: spnorm(NESP)     = 0.0 ! Factors for output normalization

      integer :: nindex(NESP) = NTRU ! Holds wavenumber
      integer :: nscatsp(NPRO)= NSPP ! Used for reduce_scatter op
      integer :: ndel(NLEV)   =    2 ! ndel for horizontal diffusion
      
      integer :: nfilter = 0 ! What kind of physics filter to use 
                             ! (0/1/2/3/4 = None/Cesaro/Exp/Lander-Hoskins/Riesz-2)
      integer :: ngptfilter = 0 ! Whether or not to filter GP->SP
    ! Filtering GP->SP may mitigate effects of truncation-scale features in physical tendencies
      integer :: nspvfilter = 0 ! Whether of not to filter SP->GP
    ! Filtering SP->GP may mitigate effects of truncation-scale features in model spectral fields
      real :: landhoskn0 = 15.0 ! Land-Hoskins filter critical wavenumber (default gives 0.1 at n=N for T21)
      integer :: nfilterexp = 8 ! Exponential filter strength
      real :: filterkappa = 8.0 ! Exponential filter 

      real, allocatable :: sdd(:,:) ! Difference between instances
      real, allocatable :: std(:,:) ! Difference between instances
      real, allocatable :: szd(:,:) ! Difference between instances
      real, allocatable :: spd(:)   ! Difference between instances

!     ************************************************
!     * Global Gridpoint Arrays (un-dimensionalized) *
!     ************************************************

#ifdef OMPSHARED
!     THE GRID FIELDS THAT CROSS THE TRANSFORM are one shared globe, and a
!     thread's name for one is a POINTER to its own band of it. SHTns wants the
!     grid whole and contiguous; the physics wants a thread to address only its
!     own latitudes. A band gives both, and it costs no memory: NPRO times NHOR
!     IS NUGP, so one shared globe is exactly the bands it replaces.
!
!     The physics does not change. NHOR is still the band from the thread's
!     side, so a whole-array statement over one of these still covers this
!     thread's latitudes and nothing else -- a whole-array expression is not a
!     call boundary and has nothing to copy. What DOES copy is one of these
!     passed wholesale to an explicit-shape dummy, because the stride between
!     levels is the globe and not the band; `probe_grid_contiguity.f90` measures
!     that, and the copy-free forms are assumed-shape, one level at a time, and
!     the base-address idiom this model already uses.
!
!     Associated once a thread by assoc_grid, called from mpstart.
      real, target :: gd_g(NUGP,NLEV) = 0. ! divergence, whole globe
      real, target :: gt_g(NUGP,NLEV) = 0. ! temperature (-t0), whole globe
      real, target :: gz_g(NUGP,NLEV) = 0. ! absolut vorticity, whole globe
      real, target :: gq_g(NUGP,NLEV) = 0. ! spec. humidity, whole globe
      real, target :: gu_g(NUGP,NLEV) = 0. ! zonal wind (*cos(phi)), whole globe
      real, target :: gv_g(NUGP,NLEV) = 0. ! meridional wind (*cos(phi)), whole globe
      real, target :: gtdt_g(NUGP,NLEV) = 0. ! t-tendency, whole globe
      real, target :: gqdt_g(NUGP,NLEV) = 0. ! q-tendency, whole globe
      real, target :: gudt_g(NUGP,NLEV) = 0. ! u-tendency, whole globe
      real, target :: gvdt_g(NUGP,NLEV) = 0. ! v-tendency, whole globe
      real, target :: gp_g(NUGP) = 0. ! surface pressure or ln(ps), whole globe

!     The nonlinear terms gridpointa builds and then transforms. They were
!     locals of that routine, contiguous per thread, and they move here for the
!     same reason the prognostic fields did: SHTns analyses the GLOBE in one
!     call and cannot be handed a band. calcgp takes them by use association
!     rather than as arguments now, because a band of a full-globe array is not
!     a contiguous (NHOR,NLEV) block and passing one to an explicit-shape dummy
!     copies it in and out again on every call, on both transform paths.
      real, target :: gtn_g(NUGP,NLEV) = 0. ! t nonlinear term, whole globe
      real, target :: gqn_g(NUGP,NLEV) = 0. ! q nonlinear term, whole globe
      real, target :: gut_g(NUGP,NLEV) = 0. ! u*t, whole globe
      real, target :: gvt_g(NUGP,NLEV) = 0. ! v*t, whole globe
      real, target :: guz_g(NUGP,NLEV) = 0. ! u forcing, whole globe
      real, target :: gvz_g(NUGP,NLEV) = 0. ! v forcing, whole globe
      real, target :: gke_g(NUGP,NLEV) = 0. ! kinetic energy, whole globe
      real, target :: guq_g(NUGP,NLEV) = 0. ! u*q, whole globe
      real, target :: gvq_g(NUGP,NLEV) = 0. ! v*q, whole globe
      real, target :: gvpp_g(NUGP) = 0. ! vertical integral of div, whole globe

      real, pointer :: gd(:,:) => NULL() ! divergence
      real, pointer :: gt(:,:) => NULL() ! temperature (-t0)
      real, pointer :: gz(:,:) => NULL() ! absolut vorticity
      real, pointer :: gq(:,:) => NULL() ! spec. humidity
      real, pointer :: gu(:,:) => NULL() ! zonal wind (*cos(phi))
      real, pointer :: gv(:,:) => NULL() ! meridional wind (*cos(phi))
      real, pointer :: gtdt(:,:) => NULL() ! t-tendency
      real, pointer :: gqdt(:,:) => NULL() ! q-tendency
      real, pointer :: gudt(:,:) => NULL() ! u-tendency
      real, pointer :: gvdt(:,:) => NULL() ! v-tendency
      real, pointer :: gp(:) => NULL() ! surface pressure or ln(ps)
      real, pointer :: gtn(:,:) => NULL() ! t nonlinear term
      real, pointer :: gqn(:,:) => NULL() ! q nonlinear term
      real, pointer :: gut(:,:) => NULL() ! u*t
      real, pointer :: gvt(:,:) => NULL() ! v*t
      real, pointer :: guz(:,:) => NULL() ! u forcing
      real, pointer :: gvz(:,:) => NULL() ! v forcing
      real, pointer :: gke(:,:) => NULL() ! kinetic energy
      real, pointer :: guq(:,:) => NULL() ! u*q
      real, pointer :: gvq(:,:) => NULL() ! v*q
      real, pointer :: gvpp(:) => NULL() ! vertical integral of divergence
#else
      real :: gtn(NHOR,NLEV)  = 0. ! t nonlinear term
      real :: gqn(NHOR,NLEV)  = 0. ! q nonlinear term
      real :: gut(NHOR,NLEV)  = 0. ! u*t
      real :: gvt(NHOR,NLEV)  = 0. ! v*t
      real :: guz(NHOR,NLEV)  = 0. ! u forcing
      real :: gvz(NHOR,NLEV)  = 0. ! v forcing
      real :: gke(NHOR,NLEV)  = 0. ! kinetic energy
      real :: guq(NHOR,NLEV)  = 0. ! u*q
      real :: gvq(NHOR,NLEV)  = 0. ! v*q
      real :: gvpp(NHOR)      = 0. ! vertical integral of divergence
      real :: gd(NHOR,NLEV)   = 0. ! divergence
      real :: gt(NHOR,NLEV)   = 0. ! temperature (-t0)
      real :: gz(NHOR,NLEV)   = 0. ! absolut vorticity
      real :: gq(NHOR,NLEV)   = 0. ! spec. humidity
      real :: gu(NHOR,NLEV)   = 0. ! zonal wind (*cos(phi))
      real :: gv(NHOR,NLEV)   = 0. ! meridional wind (*cos(phi))
      real :: gtdt(NHOR,NLEV) = 0. ! t-tendency
      real :: gqdt(NHOR,NLEV) = 0. ! q-tendency
      real :: gudt(NHOR,NLEV) = 0. ! u-tendency
      real :: gvdt(NHOR,NLEV) = 0. ! v-tendency
      real :: gp(NHOR)        = 0. ! surface pressure or ln(ps)
#endif
      real :: gpj(NHOR)       = 0. ! dln(ps)/dphi

      real :: rcsq(NHOR)      = 0. ! 1/cos(phi)**2

!     *********************************************
!     * Global Gridpoint Arrays (dimensionalized) *
!     *********************************************

      real :: dt(NHOR,NLEP)   = 0.     ! temperature 
      real :: dq(NHOR,NLEP)   = 0.     ! spec. humidity
      real :: mmr(NHOR,NLEP) = 0.      ! Aerosol array (kg/kg)
      real :: nrho(NHOR,NLEP) = 0.     ! Aerosol array (particles/m3)
      real :: du(NHOR,NLEP)   = 0.     ! zonal wind [m/s]
      real :: dv(NHOR,NLEP)   = 0.     ! meridional wind [m/s]
      real :: dp(NHOR)        = 0.     ! surface pressure
      real :: dqsat(NHOR,NLEP)= 0.     ! saturation humidity
      real :: dqt(NHOR,NLEP)  = 0.     ! adiabatic q-tendencies (for eg kuo)
      real :: mmrt(NHOR,NLEP) = 0.     ! mmr tendency array
      real :: dcc(NHOR,NLEP)  = 0.     ! cloud cover
      real :: dql(NHOR,NLEP)  = 0.     ! Liquid water content
      real :: dqo3(NHOR,NLEV) = 0.     ! ozon concentration (kg/kg)
      real :: dqco2(NHOR,NLEV)= 0.     ! co2 concentration (ppmv)
      real :: dw(NHOR,NLEV)   = 0.     ! vertical velocity (dp/dt)
      real :: dtdt(NHOR,NLEP) = 0.     ! t-tendency
      real :: dqdt(NHOR,NLEP) = 0.     ! q-tendency
      real :: dudt(NHOR,NLEP) = 0.     ! u-tendency
      real :: dvdt(NHOR,NLEP) = 0.     ! v-tendency
      real :: dp0(NHOR)       = 0.     ! surface pressure at time t
      real :: du0(NHOR,NLEP)  = 0.     ! zonal wind at time t 
      real :: dv0(NHOR,NLEP)  = 0.     ! meridional wind at time t 
      real :: dtrace(NLON,NLAT,NLEV,NTRACE) = 1.0 ! Trace array
      real :: daeros(NLON,NLAT,NLEV,NAERO) = 0.0 ! Aerosol array (kg/kg) - for aerocore
      real :: numrhos(NLON,NLAT,NLEV,NAERO) = 0.0 ! Aerosol array (particles/m3) - for radmod
      
      real :: mint(NHOR) = 0.0 !Minimum troposphere temperature 
      
      ! Hurricane metrics
      real :: gpi(NHOR)    =    0.0  ! Genesis Potential Index
      real :: venti(NHOR)  =    1.0  ! Ventilation Index
      real :: laav(NHOR)   =    0.0  ! Lower Atmospheric Absolute Vorticity
      real :: mpoti(NHOR)  =    0.0  ! Max Potential Intensity
      real :: vrmpi(NHOR)  =    0.0  ! Ventilation-reduced maximum potential intensity
      real :: capen(NHOR)  =   -1.0  ! Convective Available Potential Energy
      real :: lnb(NHOR)    =   -1.0  ! Level of neutral buoyancy (hPa)
      real :: chim(NHOR)   =   -1.0  ! Tropospheric entropy deficit
      real :: agpi(NHOR)    =    0.0  ! Acc. Genesis Potential Index
      real :: aventi(NHOR)  =    0.0  ! Acc. Ventilation Index
      real :: alaav(NHOR)   =    0.0  ! Acc. Lower Atmospheric Absolute Vorticity
      real :: ampoti(NHOR)  =    0.0  ! Acc. Max Potential Intensity
      real :: avrmpi(NHOR)  =    0.0  ! Acc. Ventilation-reduced maximum potential intensity
      real :: acapen(NHOR)  =    0.0  ! Acc. Convective Available Potential Energy
      real :: alnb(NHOR)    =    0.0  ! Acc. Level of neutral buoyancy (hPa)
      real :: achim(NHOR)   =    0.0  ! Acc. Tropospheric entropy deficit

!     *************
!     * Radiation *
!     *************

      real :: dalb(NHOR)               ! albedo
      real :: dsalb(2,NHOR)               ! spectral weighted albedo
      real :: dsnowalb(2)           = 0.6  ! spectral weighted snow albedo
      real :: dgroundalb(2)         = 0.2  ! spectral weighted ground albedo
      real :: doceanalb(2)          = 0.069  ! spectral weighted ocean albedo
      real :: dsnowalbmx(2)         = 0.8
      real :: dsnowalbmn(2)         = 0.4
      real :: dicealbmx(2)          = 0.7
      real :: dicealbmn(2)          = 0.5
      real :: dglacalbmn(2)         = 0.6
      real :: dswfl(NHOR,NLEP)         ! net solar radiation
      real :: dlwfl(NHOR,NLEP)         ! net thermal radiation
      real :: dflux(NHOR,NLEP)         ! net radiation (SW + LW)
      real :: dfu(NHOR,NLEP)           ! solar radiation upward
      real :: dfd(NHOR,NLEP)           ! solar radiation downward
      real :: dftu(NHOR,NLEP)          ! thermal radiation upward
      real :: dftd(NHOR,NLEP)          ! thermal radiation downward
      real :: dtdtlwr(NHOR,NLEV)       ! lwr temperature tendencies
      real :: dtdtswr(NHOR,NLEV)       ! swr temperature tendencies
      real :: dconv(NHOR,NLEV)         ! flux convergence (-dF/dz)

!     ***********
!     * SURFACE *
!     ***********

      real :: drhs(NHOR)  = 0.  ! surface wetness
      real :: dls(NHOR)   = 1.  ! land(1)/sea(0) mask
      real :: dz0(NHOR)   = 0.  ! rougthness length
      real :: diced(NHOR) = 0.  ! ice thickness
      real :: dicec(NHOR) = 0.  ! ice cover
      real :: dtaux(NHOR) = 0.  ! x-surface wind stress
      real :: dtauy(NHOR) = 0.  ! y-surface wind stress
      real :: dust3(NHOR) = 0.  ! u-star**3 (needed eg. for coupling)
      real :: dshfl(NHOR) = 0.  ! surface sensible heat flx
      real :: dlhfl(NHOR) = 0.  ! surface latent heat flx
      real :: devap(NHOR) = 0.  ! surface evaporation
      real :: dtsa(NHOR)  = 0.  ! surface air temperature
      real :: dmld(NHOR)  = 0.  ! mixed-layer depth (output from ocean)
!
      real dshdt(NHOR),dlhdt(NHOR)
      
      real :: tdipolep(NHOR) = 0.  ! Tidally-locked temp dipole (partial) [K]
      real :: tdipole(NUGP) = 0.  ! Tidally-locked temperature dipole [K]
      

!     *********
!     * WATER *
!     *********

      real dprc(NHOR)       ! Convective Precip   (m/s)
      real dprl(NHOR)       ! Large Scale Precip  (m/s)
      real dprs(NHOR)       ! Snow Fall           (m/s)
      real dqvi(NHOR)       ! vertical integrated specific humidity (kg/m**2)

!     *********
!     * BIOME *
!     *********

      real :: dforest(NHOR) = 0.5  ! forest cover (fract.)
      real :: dwmax(NHOR)   = 0.0  ! field capacity (m)

!     ********
!     * SOIL *
!     ********

      real :: dwatc(NHOR)   = 0.  ! soil wetness (m)
      real :: drunoff(NHOR) = 0.  ! surface runoff (m/s)
      real :: dsnow(NHOR)   = 0.  ! snow depth (m)
      real :: dsmelt(NHOR)  = 0.  ! snow melt (m/s water eq.)
      real :: dsndch(NHOR)  = 0.  ! snow depth change (m/s water eq.)
      real :: dtsoil(NHOR)  = 0.  ! soil temperature uppermost level (K)
      real :: dtd2(NHOR)    = 0.  ! soil temperature level 2 (K)
      real :: dtd3(NHOR)    = 0.  ! soil temperature level 3 (K)
      real :: dtd4(NHOR)    = 0.  ! soil temperature level 4 (K)
      real :: dtd5(NHOR)    = 0.  ! soil temperature lowermost level (K)
      real :: dglac(NHOR)   = 0.  ! glacier mask (0.,1.)

!     *********************
!     * Diagnostic Arrays *
!     *********************

      integer :: ndl(NLEV) = 0

      real :: csu(NLAT,NLEV),csv(NLAT,NLEV),cst(NLAT,NLEV),csm(NLAT,NLEV)
      real :: ccc(NLAT,NLEV)
      real :: span(NESP)

      real, allocatable :: dgp2d(:,:),dsp2d(:,:)     ! 2-d diagnostics
      real, allocatable :: dgp3d(:,:,:),dsp3d(:,:,:) ! 3-d diagnostics
      real, allocatable :: dclforc(:,:)   ! cloud forcing diagnostics
      real, allocatable :: dentropy(:,:)  ! entropy diagnostics
      real, allocatable :: dentro3d(:,:,:)! entropy diagnostics 3d
      real, allocatable :: denergy(:,:)   ! energy diagnostics
      real, allocatable :: dener3d(:,:,:) ! energy diagnostics 3d
      real, allocatable :: adenergy(:,:)   ! accumulated energy diagnostics
      real, allocatable :: adener3d(:,:,:) ! accumulated energy diagnostics 3d
      real, allocatable :: dentrop(:)     ! ps for entropy diagnostics
      real, allocatable :: dentrot(:,:)   ! t for entropy diagnostics
      real, allocatable :: dentroq(:,:)   ! q for entropy diagnostics
      real, allocatable :: dentro(:)      ! 2d entropy for diagnostics

!
!     accumulated output
!
      real :: aorbnu = 0. !Accumulated true anomaly
      real :: alambm = 0. !Accumulated solar ecliptic longitude
      real :: arasc  = 0. !Accumualted solar right ascension
      real :: azdecl = 0. !Accumualted solar declination angle
      real :: ardist = 0. !Accumulated solar distance modulus

      real :: aevap(NHOR) = 0. ! acculumated evaporation
      real :: aprl(NHOR)  = 0. ! acculumated lage scale precip.
      real :: aprc(NHOR)  = 0. ! acculumated convective precip.
      real :: aprs(NHOR)  = 0. ! acculumated snow fall
      real :: ashfl(NHOR) = 0. ! acculumated sensible heat flux
      real :: alhfl(NHOR) = 0. ! acculumated latent heat flux
      real :: aroff(NHOR) = 0. ! acculumated surface runoff
      real :: asmelt(NHOR)= 0. ! acculumated snow melt
      real :: asndch(NHOR)= 0. ! acculumated snow depth change
      real :: acc(NHOR)   = 0. ! acculumated total cloud cover
      real :: assol(NHOR) = 0. ! acculumated surface solar radiation
      real :: asthr(NHOR) = 0. ! acculumated surface thermal radiation
      real :: atsol(NHOR) = 0. ! acculumated top solar radiation
      real :: atthr(NHOR) = 0. ! acculumated top thermal radiation
      real :: assolu(NHOR)= 0. ! acculumated surface solar radiation upward
      real :: asthru(NHOR)= 0. ! acculumated surface thermal radiation upward
      real :: atsolu(NHOR)= 0. ! acculumated top solar radiation upward
      real :: ataux(NHOR) = 0. ! acculumated zonal wind stress
      real :: atauy(NHOR) = 0. ! acculumated meridional wind stress
      real :: aqvi(NHOR)  = 0. ! acculumated vertical integrated q
      real :: atsa(NHOR)  = 0. ! accumulated surface air temperature
      real :: atsama(NHOR)= 0. ! maximum surface air temperature
      real :: atsami(NHOR)= 0. ! minimum surface air temperature
      real :: ats0(NHOR)  = 0. ! accumulated surface temperature
      real :: azmuz(NHOR) = 0. ! mean cosine of solar zenith angle
      
      real :: sigrain(NHOR) = 0. !instantaneous weathering-significant precipitation [mm/day]
      real :: asigrain(NHOR) = 0. !accumulated weathering-significant precipitation [mm/day]
      real :: tempmax(NHOR) = 0. !accumulated maximum temperature
      real :: tempmin(NHOR) = 1.0e3 !accumulated minimum temperature
      
      real :: aaso(NESP)          = 0. !Accumulated quantities
      real :: aasp(NESP)          = 0.
      real :: aast(NESP,NLEV)     = 0.
      real :: aasqout(NESP,NLEV)  = 0.
      real :: aasd(NESP,NLEV)     = 0.
      real :: aasz(NESP,NLEV)     = 0.
      real :: aadq(NHOR,NLEP)     = 0.
      real :: aammr(NHOR,NLEP)    = 0. ! Accumulated aerosol mmr (kg/kg)
      real :: aanrho(NHOR,NLEP)   = 0. ! Accumulated aerosol number density (particles/m3)
      real :: aadmld(NHOR)        = 0.
      real :: aadt(NHOR,NLEP)     = 0.
      real :: aadwatc(NHOR)       = 0.
      real :: aadsnow(NHOR)       = 0.
      real :: aadql(NHOR,NLEP)    = 0.
      real :: aadust3(NHOR)       = 0.
      real :: aadcc(NHOR,NLEP)    = 0.
      real :: aadtd5(NHOR)        = 0.
      real :: aadls(NHOR)         = 0.
      real :: aadz0(NHOR)         = 0.
      real :: aadalb(NHOR)        = 0.
      real :: aadsalb1(NHOR)        = 0.
      real :: aadsalb2(NHOR)        = 0.
      real :: aadtsoil(NHOR)      = 0.
      real :: aadtd2(NHOR)        = 0.
      real :: aadtd3(NHOR)        = 0.
      real :: aadtd4(NHOR)        = 0.
      real :: aadicec(NHOR)       = 0.
      real :: aadiced(NHOR)       = 0.
      real :: aadforest(NHOR)     = 0.
      real :: aadwmax(NHOR)       = 0.
      real :: aadglac(NHOR)       = 0.
      real :: aadqo3(NHOR,NLEV)   = 0.
      real :: aagroundoro(NHOR)   = 0.
      real :: aaglacieroro(NHOR)  = 0.
      
      
!     *******************
!     * Latitude Arrays *
!     *******************

      character(len=3) chlat(NLAT)
      real (kind=8) :: sid(NLAT)     ! sin(phi)
      real (kind=8) :: gwd(NLAT)     ! Gaussian weights
      real :: csq(NLAT)              ! cos(phi)**2
      real :: cola(NLAT)             ! cos(phi)
      real :: rcs(NLAT)              ! 1 / cos(phi)
      real :: deglat(NLPP)           ! latitude in degrees

!     ****************
!     * Level Arrays *
!     ****************

      real :: tdissd(NLEV)  =  0.20 ! diffusion time scale for divergence [days]
      real :: tdissz(NLEV)  =  1.10 ! diffusion time scale for vorticity [days]
      real :: tdisst(NLEV)  =  5.60 ! diffusion time scale for temperature [days]
      real :: tdissq(NLEV)  =  0.1  ! diffusion time scale for sp. humidity [days]
      real :: restim(NLEV)  =   0.0
      real :: t0(NLEV)      = 250.0
      real :: tfrc(NLEV)    =   0.0
      real :: sigh(NLEV)    =   0.0
      real damp(NLEV)
      real dsigma(NLEV)
      real rdsig(NLEV)
      real sigma(NLEV)
      real sigmah(NLEV)
      real t01s2(NLEV)
      real tkp(NLEV)                                
      real c(NLEV,NLEV)
      real g(NLEV,NLEV)
      real tau(NLEV,NLEV)
      real bm1(NLEV,NLEV,NTRU)

!     ******************
!     * Parallel Stuff *
!     ******************

      integer :: mpinfo  = 0
      integer :: mypid   = 0
!     Which slot of the tendency partials above is this process's own. It is
!     mypid where the processes share an address space and 0 where they do not,
!     so the arrays are indexed the same way in both builds and neither the
!     routines that fill them nor the reduction has to know which build it is.
!     Set by mpstart, whichever mpimod supplies it.
      integer :: mypart  = 0
      integer :: myworld = 0
      integer :: nproc   = NPRO
      character (80),allocatable :: ympname(:)

      ! ***************
      ! * Random seed *
      ! ***************
      
      integer              :: seed(8) = 0 ! settable in namelist
      integer, allocatable :: meed(:)     ! machine dependent seed
      
!     **********************
!     * Multirun variables *
!     **********************

      integer :: mrworld =  0   ! MPI communication
      integer :: mrinfo  =  0   ! MPI info
      integer :: mrpid   = -1   ! MPI instance id
      integer :: mrnum   =  0   ! MPI number of instances
      integer :: mintru  =  0   ! Lowest resolution of all instances
      integer :: mrdim   =  0   ! Exchange dimension  (min. NRSP)
      integer :: nsync   =  0   ! Synchronization on or off
      integer, allocatable :: mrtru(:) ! Truncations of members
      
      real    :: syncstr  =  0.0 ! Coupling strength (0 .. 1)
      real    :: synctime =  0.0 ! Coupling time [days]

!     **********************************************
!     * version identifier (to be set in puma.f90) *
!     * will be printed in puma sub *prolog*       *
!     **********************************************

      character(len=80) :: plasimversion

!     ******************************
!     * Planet dependent variables *
!     ******************************


      character(len=80) :: yplanet=" " ! Planet name

      integer :: nfixorb = 0           ! Global switch to fix orbit
      integer :: ngenkeplerian = 0     ! If 1, compute a general Keplerian orbit 

      real :: akap   = 0.0             ! Kappa
      real :: alr    = 0.0             ! Lapse rate
      real :: ga     = 0.0             ! Gravity
      real :: gascon = 0.0             ! Gas constant for dry air
      real :: plarad = 0.0             ! Planet radius
      real :: pnu    = 0.0             ! Time filter
      real :: sidereal_day  = 0.0      ! Length of sidereal day [sec]
      real :: solar_day     = 0.0      ! Length of solar day [sec]
      real :: day_24hr      = 86400.0  ! Seconds in a standard day
      real :: sidereal_year = 0.0      ! Length of sidereal year [sec]
      real :: tropical_year = 0.0      ! Length of tropical year [sec]
      real :: ww     = 0.0             ! Omega used for scaling
      real :: oroscale = 1.0           ! Orography scaling
      real :: ra1    = 0.0             !
      real :: ra2    = 0.0             !
      real :: ra4    = 0.0             !
      real :: acpd   = 0.0             ! acpd = gascon / akap ! Specific heat for dry air
      real :: adv    = 0.0             ! acpv / acpd - 1.0
      real :: cv     = 0.0             ! cv = plarad * ww
      real :: ct     = 0.0             ! ct = CV * CV / gascon
      real :: pnu21  = 0.0             ! pnu21 = 1.0 - 2.0 * pnu ! Time filter 2
      real :: rdbrv  = 0.0             ! rdbrv = gascon / RV  ! Rd/Rv
      real :: rotspd = 1.0             ! rotation speed (factor)
      real :: eccen  = 0.0             ! Eccentricity of Orbit
      real :: obliq  = 0.0             ! Obliquity of Orbit
      real :: mvelp  = 0.0             ! Longitude of moving vernal equinox
      real :: meananom0 = 0.0          ! Initial mean anomaly in degrees

!     ******************************************
!     * GUI (Graphical User Interface for X11) *
!     ******************************************

      integer, parameter :: PUMA   = 0
      integer, parameter :: SAM    = 1
      integer, parameter :: PLASIM = 2
      integer, parameter :: PLALSG = 3
      parameter (NPARCS = 5)          ! Number of GUI parameters
      character(6) :: yguinam(NPARCS) ! Variable names for GUI display
      integer(kind=4) :: nguidbg   = 0        ! 1: GUI debug printout
      integer(kind=4) :: model     = PLASIM
      integer :: nshutdown = 0        ! Flag for shutdown request
      integer :: ntpal     = 2        ! Color pallette for temperature
      real(kind=4) :: parc(NPARCS)            ! Values of GUI parameters
      real(kind=4) :: crap(NPARCS)            ! Backup of parc(NPARCS)
      real(kind=4) :: guimin(NPARCS)          ! lower limit
      real(kind=4) :: guimax(NPARCS)          ! upper limit
      real(kind=4) :: guiinc(NPARCS)          ! increment
      real, allocatable :: sr1(:,:)   ! dummy array for PUMA guimod compatibility
      real, allocatable :: sr2(:,:)   ! dummy array for PUMA guimod compatibility
      logical :: ldisp   = .FALSE.    ! DISP changed by GUI
      logical :: ldtep   = .FALSE.    ! DTEP changed by GUI
      logical :: ldtns   = .FALSE.    ! DTNS changed by GUI
      logical :: lrotspd = .FALSE.    ! rotspd changed by GUI


!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(aadalb,aadcc,aadforest,aadglac,aadicec,aadiced,aadls,aadmld,aadq,aadql,aadqo3,&
!$omp&  aadsalb1,aadsalb2,aadsnow,aadt,aadtd2,aadtd3,aadtd4,aadtd5,aadtsoil,aadust3,aadwatc,aadwmax,&
!$omp&  aadz0,aaglacieroro,aagroundoro,aammr,aanrho,aasd,aaso,aasp,aasqout,aast,aasz,acapen,acc,&
!$omp&  achim,acpd,adener3d,adenergy,adv,aero_namelist,aevap,agpi,akap,alaav,alambm,alhfl,alnb,alr,&
!$omp&  als,alv,ampoti,aorbnu,aprc,aprl,aprs,aqvi,arasc,ardist,aroff,ashfl,asigrain,asmelt,asndch,&
!$omp&  assol,assolu,asthr,asthru,ataux,atauy,ats0,atsa,atsama,atsami,atsol,atsolu,atthr,aventi,&
!$omp&  avrmpi,azdecl,azmuz,bm1,c,capen,ccc,chim,chlat,co2,cola,crap,csm,csq,cst,csu,csv,ct,cv,&
!$omp&  daeros,dalb,damp,dampsp,dawn,day_24hr,dcc,dclforc,dconv,deglat,delt,delt2,deltsec,deltsec2,&
!$omp&  dener3d,denergy,dentro,dentro3d,dentrop,dentropy,dentroq,dentrot,devap,dfd,dflux,dforest,&
!$omp&  dftd,dftu,dfu,dglac,dglacalbmn,dgp2d,dgp3d,dgroundalb,dicealbmn,dicealbmx,dicec,diced,dlhdt,&
!$omp&  dlhfl,dls,dlwfl,dmld,doceanalb,dp,dp0,dprc,dprl,dprs,dq,dqco2,dqdt,dql,dqo3,dqsat,dqt,dqvi,&
!$omp&  drhs,drunoff,dsalb,dshdt,dshfl,dsigma,dsmelt,dsndch,dsnow,dsnowalb,dsnowalbmn,dsnowalbmx,&
!$omp&  dsp2d,dsp3d,dswfl,dt,dtaux,dtauy,dtd2,dtd3,dtd4,dtd5,dtdt,dtdtlwr,dtdtswr,dtep,dtns,dtrace,&
!$omp&  dtrop,dtsa,dtsoil,dttl,dttrp,du,du0,dudt,dust3,dv,dv0,dvdt,dw,dwatc,dwmax,dz0,eccen,&
!$omp&  efficiency_dat,evap,filterkappa,fixedlon,fluxmod_namelist,frcmod,g,ga,gascon,gd,gp,gpi,&
!$omp&  gpimax,gpj,gq,gqdt,gqn,gtn,gut,gvt,guz,gvz,gke,guq,gvq,gvpp,&
!$omp&  gt,gtdt,gu,gudt,guiinc,guimax,guimin,gv,gvdt,gwd,gz,&
!$omp&  hcendstep,hcinterval,&
!$omp&  hcstartstep,ice_output,icemod_namelist,kick,l_aero,laav,laavmax,landhoskn0,landmod_namelist,&
!$omp&  ldisp,ldtep,ldtns,lnb,lrotspd,m_days_per_month,m_days_per_year,mars,mcal_days_per_year,&
!$omp&  meananom0,meed,mint,mintru,miscmod_namelist,mmr,mmrt,mocd,model,mpinfo,mpoti,mpotimax,&
!$omp&  mpstep,mrdim,mrinfo,mrnum,mrpid,mrtru,mrworld,mstep,mtspd,mvelp,mypart,mypid,myworld,&
!$omp&  n_days_per_month,n_days_per_year,n_run_days,n_run_months,n_run_steps,n_run_years,&
!$omp&  n_sea_points,n_start_month,n_start_step,n_start_year,n_steps_per_year,naccuout,nadv,nafter,&
!$omp&  naqua,ncoeff,ndatim,ndel,ndesert,ndheat,ndiag,ndiagcf,ndiaggp,ndiaggp2d,ndiaggp3d,ndiagsp,&
!$omp&  ndiagsp2d,ndiagsp3d,ndivdamp,ndl,nener3d,nenergy,nentro3d,nentropy,neqsig,nfilter,&
!$omp&  nfilterexp,nfixorb,nflux,ngenkeplerian,nglspec,ngptfilter,ngui,nguidbg,nhcadence,nhcstp,&
!$omp&  nhdiff,nhordif,nhurricane,nindex,nkits,nlowio,noutput,npackgp,npacksp,nperpetual,nprhor,&
!$omp&  nprint,nproc,nqspec,nrad,nrdrag,nrestart,nrho,nscatsp,nseedlen,nsela,&
!$omp&  nshtns,nshutdown,nsnapshot,&
!$omp&  nspinit,nsponge,nspvfilter,nstep,nstep1,nstps,nstpw,nstratosponge,nsync,ntime,ntpal,ntspd,&
!$omp&  nud,numrhos,nut,nveg,nwpd,nwritehurricane,obliq,ocean_output,oceanmod_namelist,olr,oroscale,&
!$omp&  parc,pfac,planet_namelist,plarad,plasim_diag,plasim_hcadence,plasim_namelist,plasim_output,&
!$omp&  plasim_restart,plasim_snapshot,plasim_status,plasimversion,plavor,pnu,pnu21,precip,psurf,&
!$omp&  ptop,ptop2,ra1,ra2,ra4,radmod_namelist,rainmod_namelist,rcs,rcsq,rdbrv,rdsig,restim,rotspd,&
!$omp&  sak,sakpp,sdd,sdipole,sdipolep,sdm,sdp,sdt,seamod_namelist,seed,sellon,sid,sidereal_day,&
!$omp&  sidereal_year,sigh,sigma,sigmah,sigrain,so,solar_day,sop,span,spd,spm,spnorm,spp,spt,sqm,&
!$omp&  sqout,sqp,sqt,sr1,sr2,srm,srp,std,stm,stp,stt,surfmod_namelist,syncstr,synctime,szd,szm,szp,&
!$omp&  szt,t0,t01s2,t2mean,tau,taucool,tdipole,tdipolep,tdissd,tdissq,tdisst,tdissz,tempmax,&
!$omp&  tempmin,tfrc,tgr,time0,tkp,tmelt,tmstart,tropical_year,umax,vegmod_namelist,venti,ventimin,&
!$omp&  vrmpi,vrmpimax,ww,yguinam,ympname,yplanet)

      contains

!     ==========================
!     SUBROUTINE ASSOC_SPECTRAL
!     ==========================

!     Point this thread's partials at its own slice of the shared arrays.
!
!     Called once a thread, from mpstart, after mypid is known. Under any build
!     but the shared one the partials are ordinary storage and this does
!     nothing -- the body compiles away with the pointers it refers to.

      subroutine assoc_spectral
#ifdef OMPSHARED
      integer :: lo, hi

      lo = mypid * NSPP + 1
      hi = lo + NSPP - 1

      sdp => sd(lo:hi,:)
      stp => st(lo:hi,:)
      szp => sz(lo:hi,:)
      sqp => sq(lo:hi,:)
      srp => sr(lo:hi,:)
      spp => sp(lo:hi)
#endif
      return
      end subroutine assoc_spectral


      subroutine assoc_grid
#ifdef OMPSHARED
      integer :: lo, hi

!     A thread's band of the shared globe. The band is the same NHOR rows the
!     thread's private array used to be, so this is a relabelling of storage and
!     nothing else -- the model computes the same numbers in the same order.
!
!     WHAT THE ROW ORDER IS, and it is the trap. Row mypid*NHOR + i is thread
!     mypid's LOCAL point i. That equals GLOBAL latitude order only when the
!     scatter is contiguous. Under LPAIRLAT it is not: the permutation gives a
!     thread a latitude and its mirror rather than a block, so the shared array
!     is in thread order and not in latitude order.
!
!     That does not affect this routine or the physics, both of which only ever
!     address a thread's own band. It affects anything that reads the array AS
!     A GLOBE -- SHTns above all, which needs latitudes in its own Gauss order.
!     SHTns and LPAIRLAT are therefore exclusive, and that is not a loss:
!     LPAIRLAT exists to let legmod fold a mirror pair together, and SHTns
!     replaces legmod.
      lo = mypid * NHOR + 1
      hi = lo + NHOR - 1

      gd   => gd_g(lo:hi,:)
      gt   => gt_g(lo:hi,:)
      gz   => gz_g(lo:hi,:)
      gq   => gq_g(lo:hi,:)
      gu   => gu_g(lo:hi,:)
      gv   => gv_g(lo:hi,:)
      gtdt => gtdt_g(lo:hi,:)
      gqdt => gqdt_g(lo:hi,:)
      gudt => gudt_g(lo:hi,:)
      gvdt => gvdt_g(lo:hi,:)
      gp   => gp_g(lo:hi)
      gtn  => gtn_g(lo:hi,:)
      gqn  => gqn_g(lo:hi,:)
      gut  => gut_g(lo:hi,:)
      gvt  => gvt_g(lo:hi,:)
      guz  => guz_g(lo:hi,:)
      gvz  => gvz_g(lo:hi,:)
      gke  => gke_g(lo:hi,:)
      guq  => guq_g(lo:hi,:)
      gvq  => gvq_g(lo:hi,:)
      gvpp => gvpp_g(lo:hi)
#endif
      return
      end subroutine assoc_grid


!     ==================
!     FUNCTION ILATPERM
!     ==================

!     Permuted slot -> global latitude, the whole of the paired decomposition.
!
!     Slot kp is the kp-th latitude of the scatter buffer, so it lands on
!     process ir = (kp-1)/NLPP at local index il = kp - ir*NLPP. The first
!     NLHP local latitudes of a process are a northern block; the rest are
!     those same latitudes' mirrors, in reverse, which is what makes local
!     il and NLPP+1-il a mirror pair.
!
!     Without LPAIRLAT this is the identity and the layout is the stock
!     contiguous one. At NPRO == 1 it is ALSO the identity, since ir is 0 and
!     NLHP is NLAT/2: slot il > NLHP maps to NLAT - NLAT + il = il.

      integer function ilatperm(kp)
      integer :: kp
      integer :: ir
      integer :: il

      if (.not. LPAIRLAT) then
         ilatperm = kp
         return
      endif

      ir = (kp - 1) / NLPP    ! process holding the slot
      il =  kp - ir * NLPP    ! its local latitude there

      if (il <= NLHP) then
         ilatperm = ir * NLHP + il                 ! northern block
      else
         ilatperm = NLAT - ir * NLHP - NLPP + il   ! the mirror of it
      endif

      return
      end function ilatperm

      end module pumamod
