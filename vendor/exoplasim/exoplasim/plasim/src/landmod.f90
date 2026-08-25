

      module landmod
      use pumamod
      use landcolumn
!
!     version identifier (date)
!
      character(len=80) :: lversion = '19.09.2019 by Adiv'
!
!     parameters
!
      parameter(NLSOIL=5)
      parameter(WSMAX_EARTH = 0.5) ! Initial value vor Earth
!
!     namelist parameters
!
!     THE SUBGRID SNOW-COVER SCALE, and it was a bare 0.01 at seven sites. The
!     snow-covered fraction of a cell is taken as dsnow/(dsnow + snowcovz), so
!     snowcovz is the snow depth at which half the cell is covered: it stands in
!     for the distribution of snow depth WITHIN the cell and therefore depends
!     on how big the cell is. 1 cm is the canonical value and carries no NLAT
!     term here or upstream. Anchored to T21; named and put in landmod_nl so a
!     rung study can move it, rather than left as a literal in seven places.
!     world-khn.
      real    :: snowcovz = 0.01  ! snow depth at half cell cover (m water eq.)
      integer :: nlandt   = 1     ! switch for land model (1/0 : prog./clim)
      integer :: nlandw   = 1     ! switch for soil model (1/0 : prog./clim)
      integer :: newsurf  = 0     ! (dtcl,dwcl) 1: update from file, 2:reset 
      integer :: nwatcini = 0     ! (0/1) initialize water content of soil
      integer :: nwetsoil = 0     ! (0/1) Soil albedo responds to water content
      real    :: alblandnl  = 0.2   ! albedo for land
      real    :: albland  = 0.22
      
      real  :: dgroundalbnl(2)  = 0.22
                   ! Non-spectral versions
      real    :: albsmin  = 0.4   ! min. albedo for snow
      real    :: albsmax  = 0.8   ! max. albedo for snow
      real    :: albsminf = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin  = 0.6   ! min. albedo for glaciers
      real    :: albgmax  = 0.8   ! max. albedo for glaciers
      real    :: alblandmax = 0.40  ! max land albedo; i.e. desert sand
                   ! lambda < 0.75 microns
      real    :: albsmin1  = 0.4   ! min. albedo for snow
      real    :: albsmax1  = 0.8   ! max. albedo for snow
      real    :: albsminf1 = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf1 = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin1  = 0.6   ! min. albedo for glaciers
      real    :: albgmax1  = 0.8   ! max. albedo for glaciers
                   ! lambda > 0.75 microns
      real    :: albsmin2  = 0.4   ! min. albedo for snow
      real    :: albsmax2  = 0.8   ! max. albedo for snow
      real    :: albsminf2 = 0.3   ! min. albedo for snow (with forest)
      real    :: albsmaxf2 = 0.4   ! max. albedo for snow (with forest)
      real    :: albgmin2  = 0.6   ! min. albedo for glaciers
      real    :: albgmax2  = 0.8   ! max. albedo for glaciers

!     Snow under a canopy, per band.
!
!     A forested snow endmember is a mixture of the exposed snow and the canopy
!     that hides the rest of it, so the only spectrum-free quantity in it is the
!     masked FRACTION. Upstream carries the mixture already evaluated, as one
!     broadband ratio (albsmaxf/albsmax = 0.5, albsminf/albsmaxf = 0.75), and a
!     ratio evaluated under one spectrum cannot be reapplied either side of
!     0.75 micron: a canopy is dark in band 1 and bright in band 2 while snow
!     runs the other way, so one factor is wrong in both bands at once.
!
!     albforest is the canopy albedo per band and forcov* the masked fraction:
!     albsmaxf_b = forcovmx*albforest(b) + (1-forcovmx)*albsmax_b, and the same
!     for the minimum with forcovmn. The defaults are what upstream's two ratios
!     imply when solved against its Earth-Sun canopy endmember 0.15:
!     0.6153846*0.15 + 0.3846154*0.8 = 0.4 = albsmaxf and 0.4*0.15 + 0.6*0.4 =
!     0.3 = albsminf, exactly. An Earth configuration therefore reproduces
!     upstream to the last digit, which is a check that can fail.
!
!     The two fractions differ because upstream's two pairs differ: multiple
!     scattering between a bright snowpack and the canopy masks more of it than
!     the same canopy masks of aged snow. They are separate keys rather than one
!     averaged fraction so that neither pair has to be discarded.
      real    :: albforest(2) = 0.15   ! canopy albedo, per band
      real    :: forcovmx = 0.6153846  ! canopy fraction masking max-alb snow
      real    :: forcovmn = 0.4        ! canopy fraction masking min-alb snow

!     The structural axes of the canopy/snow mask, all inert at these values.
!
!     forcovmx and forcovmn above make the hidden fraction dforest*forcov: a
!     constant times the cover, blind to how deep the snow is against the
!     canopy, to how much plant area the canopy carries, and to whether the
!     canopy is itself holding snow. snowmaskmod generalises that along three
!     axes and returns the same number to the last bit at the values below,
!     which is the check that can fail.
!
!     They are SCALARS standing in for fields. The canopy cover they multiply,
!     dforest, arrives per gridcell as surface code 212; a canopy height and a
!     plant area index vary the same way and will arrive the same way once the
!     vegetation component reports them, at which point these become the
!     fallback for a cell the field does not cover. BIO-33 is that wiring.
!
!     They are inert here because this world does not report their inputs yet.
!     forhgt needs a canopy height, which is GRAV-7's; forpai needs a plant or
!     stem area index per gridcell, which is what the vegetation component will
!     report under BIO-17; forint needs a canopy snow store with a mass balance
!     rather than a declared standing fraction, which is BIO-32. Turning any of
!     them on moves the surface energy balance of every snow-covered forested
!     cell, so it is a declared decision and not a default. BIO-30.
!
!     forext is the extinction coefficient of the canopy gap fraction,
!     exp(-forext*forpai), read only when forpai is positive so it changes
!     nothing at these values. One is the coefficient in the gap fraction
!     Essery (2013) reports for the scheme it evaluates, with the plant area
!     index counting leaves and stems. A leaf-angle-resolved reading would put
!     the geometric factor in the coefficient instead, near a half for a
!     randomly oriented canopy under diffuse light, and would want twice the
!     plant area index for the same masking; which convention the vegetation
!     component reports its index in therefore has to travel with it.
      real    :: forhgt   = -1.0  ! canopy height for snow burial (m; <=0 off)
      real    :: forpai   = -1.0  ! canopy plant area index (<=0 off)
      real    :: forext   =  1.0  ! extinction coefficient of the gap fraction
      real    :: forint   =  0.0  ! intercepted-snow fraction of the canopy

!     River routing, used by roffini here and by oroini in glaciermod.
!
!     u = zcvel/zdx * |grad(zoro)|**roffexp, and zoro is GEOPOTENTIAL, so the
!     slope term already carries ga**roffexp while an open-channel velocity
!     scales as sqrt(ga). Upstream's zcvel = 4.2 absorbs neither. The
!     coefficient multiplying the geopotential slope is roffvel*ga**(0.5-roffexp),
!     which leaves the whole expression proportional to sqrt(ga) once the
!     geopotential slope is divided back to a topographic one. At ga = 9.80665
!     that product is 4.200, upstream's fitted value, so an Earth configuration
!     is unchanged.
!
!     roffpit is the elevation added to a local minimum per pit-filling pass.
!     Upstream adds 1 m2/s2, which is a HEIGHT only after division by ga; this
!     is the height, and 0.101972 m is upstream's increment at Earth gravity.
      real    :: roffvel  = 2.022845 ! routing velocity coefficient (per sqrt(g))
      real    :: roffexp  = 0.18     ! routing velocity slope exponent
      real    :: roffpit  = 0.101972 ! pit-filling elevation increment (m)

      real    :: dz0land  = 2.0   ! roughness length land
      real    :: drhsland = 0.25  ! wetness factor land
      real    :: drhsfull = 0.4   ! threshold above which drhs=1 [frac. of wsmax]
      real    :: dzglac   = -1.   ! threshold of orography to be glacier (-1=none)
      real    :: dztop    = 0.20  ! thickness of the uppermost soil layer (m)
      real    :: dsmax    = 5.00  ! maximum snow depth (m-h20; -1 = no limit)

      
      real    :: wsmax    = WSMAX_EARTH ! max field capacity of soil water (m)
      real    :: dwatcini = 0           ! water content of soil (m)

!     THE LAND LIQUID WATER SCHEME, LSHY-3. A SELECTION, and the default is
!     the scheme that has always run.
!
!       nlandwcol = 0  the scalar bucket. One store, one capacity, runoff is
!                      the store overflowing. Bit-identical to what this model
!                      did before the selection existed, and that is checked
!                      rather than asserted: exoplasim/scripts/
!                      verify_land_column_reduction.sh drives both kernels with
!                      one flux sequence and requires bitwise equality.
!       nlandwcol = 1  the layered column, nlsoilw layers deep. At nlsoilw = 1
!                      with an impermeable base it IS the bucket, operation for
!                      operation, which is what makes it a reduction and not a
!                      resemblance.
!
!     The three registered hypotheses of LSHY-3 are the bucket, this
!     parsimonious multilayer column, and gradient-driven flow. The third is
!     NOT here: it needs an unsaturated conductivity and a matric potential,
!     and pedology/config/land_column_properties.yaml carries both as
!     undeclared. Registering two and naming the third is the honest state.
      integer :: nlandwcol   = 0   ! 0: scalar bucket, 1: layered column
      integer :: nlsoilw     = 1   ! liquid water layers when nlandwcol = 1
      integer :: nlandwdrain = 0   ! lower boundary: 0 impermeable, 1 free drain

!     Layer capacity as a fraction of dwmax. Only the first nlsoilw entries are
!     read and they are renormalised to sum to one, so a partial list is a
!     shape rather than an error. The default puts everything in one layer,
!     which is the reduction.
      real    :: dsoilwf(NLSOILWX) = (/1.0, 0.0, 0.0, 0.0,                    &
     &                                 0.0, 0.0, 0.0, 0.0/)

!     THE EVAPORATION LIMITER, as three axes rather than one hardcoded shape.
!     beta = min(1, max((theta - drhslow)/(drhsfull - drhslow), 0)) ** nrhsexp
!     with theta the store as a fraction of capacity. The defaults below are
!     this model's active form exactly, and the default path is taken by BRANCH
!     and not by algebra, because x**1.0 with a real exponent is not bitwise x.
!     nrhsexp is an integer for the same reason.
!
!     The bracket the other end of this sits at: cGENIE's ENTS uses nrhsexp = 4
!     with the knee at a full store, which at equal fractional fill differs
!     from the form below by a factor of 39 at 0.4 of capacity and 2 at 0.8,
!     agreeing only when full. Neither is right. Making it selectable is what
!     turns that into a runtime bracket instead of a code fork.
      real    :: drhslow  = 0.0   ! theta_low, fraction of capacity
      integer :: nrhsexp  = 1     ! the exponent c

!     SOIL PHASE, LSHY-5. Off by default, which is what this model does today:
!     its five soil temperature layers carry no water and no phase, so melt
!     water always infiltrates whatever the soil temperature and freeze/thaw
!     neither absorbs nor releases latent heat in the ground. LPJ-GUESS carries
!     an ice fraction per layer and reduces available liquid under freezing, so
!     the two columns presently disagree about whether water that reached the
!     ground is liquid.
!
!     nlandwphase = 1 needs nlandwcol = 1: phase is a property of a LAYER, and
!     the scalar bucket has no layer to freeze. Ice occupies pore space, so it
!     comes off the layer's capacity in `column_step` and a frozen layer
!     overflows sooner; that is the infiltration impedance, and it needs no
!     conductivity, which this column does not have.
      integer :: nlandwphase = 0  ! 0: no soil ice, 1: freeze and thaw

!     THE WATER COLUMN'S THICKNESSES, in metres, which the model has never had.
!     `dwmax` is a CAPACITY in metres of water and says nothing about how deep
!     the column is, and `dsoilwf` is a share of that capacity, not a thickness.
!     Phase needs a depth, because the temperature that decides it lives on the
!     SOIL TEMPERATURE layers, `dsoilz`, and those share no boundary with the
!     water layers. The mapping is by midpoint: a water layer takes the
!     temperature of whichever temperature layer contains its centre.
!
!     The default is one layer of 1.5 m, which is
!     pedology/config/land_column_properties.yaml's declared column base and
!     LPJ-GUESS's physical profile. Inert while nlandwphase = 0.
      real    :: dsoilwz(NLSOILWX) = (/1.5, 0.0, 0.0, 0.0,                    &
     &                                 0.0, 0.0, 0.0, 0.0/)

!     SIMBA - fixed parameters

      real    :: rlue     =  3.4E-10 ! Recommended by Pablo Paiewonsky
      real    :: co2conv  =  8.3E-4  ! Recommended by Pablo Paiewonsky
      real    :: tau_veg  = 10.0  ! [years] - in landini scaled to seconds
      real    :: tau_soil = 42.0  ! [years] - in landini scaled to seconds
      real    :: rinifor  =  0.5
      real    :: rnbiocats=  0.0
!
!     Surface thermal scalars (snow similar to the sea ice module).
!
!     All seven are namelist keys, in landmod_nl below. They were compiled-in
!     constants, so a run could not say what thermal inertia its land carried
!     and could not vary it; soildiff in particular was reachable from nowhere
!     at all while soilcap was reachable only through cpsoil.
!
!     They are SCALARS and the model has no per-cell field for any of them, so
!     one value covers every lithology. soildiff and soilcap are moist mineral
!     soil: thermal inertia sqrt(k*rho*c) = 2078 J/m2/K/s**0.5, against roughly
!     625 for a dry playa or salt crust, so a surface dominated by evaporite
!     and playa clastics is damped by about 3.3x too much. Giving those classes
!     their own inertia needs a field, not a different scalar.
!
!     rhosnow is a settled snow density and converts water equivalent to the
!     physical snow thickness that insulates the soil column below; rhoglac in
!     glaciermod does the same for ice thickness in the orography. Both are set
!     by overburden compaction, which scales with gravity, so both are low for a
!     planet with stronger surface gravity than the one they were measured on.
!
      real :: rhosnow  = 330.    ! snow density (kg/m**3)
      real :: soildiff = 1.8     ! heat diffusivity of the soil (W/m/K)
      real :: sicediff = 2.03    ! heat diffusivity of ice      (W/m/K)
      real :: snowdiff = 0.31    ! heat diffusivity of snow     (W/m/K)
      real :: soilcap  = 2.4E6   ! heat capacity of the soil  (J/m**3/K)
      real :: sicecap  = 2.07E6  ! heat capacity of ice       (J/m**3/K)
      real :: snowcap  = 0.6897E6! heat capacity of snow      (J/m**3/K)
!
!     global arrays
!
!     a) surface definitions
!
      real :: doro(NHOR) = 0.0     ! orography (m2/s2)
      real :: dts(NHOR)  = 1.0e20  ! surface temperature (K)
      real :: dtsm(NHOR) = 1.0e20  ! surface temperature (K)
      real :: dqs(NHOR)  = 1.0e20  ! surface humidity    (kg/kg)
!
!     b) runoff
!
      real :: driver(NHOR) = 0.0   ! surface water (river runoff) (m)
      real :: duroff(NHOR) = 0.0   ! zonal runoff-velocity  (1/s)
      real :: dvroff(NHOR) = 0.0   ! meridional runoff-velocity  (1/s)
      real :: darea(NHOR)  = 1.0   ! area weights
!
!     c) soil
!
      real :: dsoilz(NLSOIL)=(/0.4,0.8,1.6,3.2,6.4/)  ! soil layer thickness (m)
      real :: dsoilt(NHOR,NLSOIL) = 0.0    ! soil temperatur (K)
      real :: dsnowt(NHOR)        = 0.0    ! snow temperatur (K)
      real :: dtclsoil(NHOR)      = 0.0    ! clim soil temp. (initilization) (K)
      real :: dsnowz(NHOR)        = 0.0    ! snow depth (m water equivalent)
      real :: dwater(NHOR)        = 0.0    ! surface water for soil (m/s)
!
!     The layered liquid store, and the drainage out of its base. Both are
!     inert while nlandwcol = 0: dwatcl is initialised from dwatc and carried
!     through the restart so that switching the scheme on does not need a cold
!     start, and ddrain is identically zero because an impermeable base has no
!     drainage. dwatc remains the COLUMN TOTAL under either scheme, because
!     fluxmod's evaporation limiter, simba's water stress, aeromod and outmod
!     all read it and none of them knows about layers.
      real :: dwatcl(NHOR,NLSOILWX) = 0.0  ! liquid water by layer (m)
      real :: dsoili(NHOR,NLSOILWX) = 0.0  ! soil ice by layer (m water equiv.)
      real :: ddrain(NHOR)          = 0.0  ! drainage out of the column base (m/s)
!
!     The drainage's output accumulator, WORLD-P9QQ. `outmod` fills it, divides
!     it by the output counter and resets it on the pattern `aroff` uses for
!     the surface runoff, and it lives HERE rather than beside `aroff` in
!     plasimmod because the land column owns the flux it accumulates. Without
!     it the ledger's `drainage` crossing names a producer the postprocessor
!     cannot deliver: `ddrain` reaches the restart and nothing else, so a run
!     on nlandwcol = 1 with nlandwdrain = 1 computes a flux no product carries.
!     Identically zero under the default impermeable base, which is why the
!     scheme could land without it.
      real :: adrain(NHOR)          = 0.0  ! accumulated drainage (m/s)
!
!     e) climatological surface
!
      real :: dtcl(NHOR,0:13)   =  0.0  ! climatological surface temperature
      real :: dwcl(NHOR,0:13)   = -1.0  ! climatological soil wetness
      real :: dalbcl(NHOR,0:13) =  0.22  ! climatological background albedo
      real :: dalbcl1(NHOR,0:13) = 0.22  ! climatological background albedo (<.75 um)
      real :: dalbcl2(NHOR,0:13) = 0.22  ! climatological background albedo (>.75 um)
      real :: dtclim(NHOR)      =  0.0  ! climatological surface temperature
      real :: dwclim(NHOR)      =  0.0  ! climatological soil wetness
      real :: dz0clim(NHOR)     =  2.0  ! climatological z0  (total)
      real :: dz0climo(NHOR)    =  0.0  ! climatological z0  (from topograhpy only)
      real :: dalbclim(NHOR)    =  2.0  ! climatological background albedo
      real :: dalbclim1(NHOR)   =  2.0  ! climatological background albedo (<.75 um)
      real :: dalbclim2(NHOR)   =  2.0  ! climatological background albedo (>.75 um)
!

!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(albforest,albgmax,albgmax1,albgmax2,albgmin,albgmin1,albgmin2,albland,alblandmax,&
!$omp&  alblandnl,albsmax,albsmax1,albsmax2,albsmaxf,albsmaxf1,albsmaxf2,albsmin,albsmin1,albsmin2,&
!$omp&  albsminf,albsminf1,albsminf2,co2conv,dalbcl,dalbcl1,dalbcl2,dalbclim,dalbclim1,dalbclim2,&
!$omp&  darea,dgroundalbnl,doro,dqs,drhsfull,drhsland,driver,dsmax,dsnowt,dsnowz,dsoilt,dsoilz,dtcl,&
!$omp&  dtclim,dtclsoil,dts,dtsm,duroff,dvroff,dwatcini,dwater,dwcl,dwclim,dz0clim,dz0climo,dz0land,&
!$omp&  dwatcl,dsoili,ddrain,adrain,dsoilwf,dsoilwz,drhslow,nlandwcol,nlsoilw,nlandwdrain,nrhsexp,nlandwphase,dzglac,dztop,&
!$omp&  forcovmn,forcovmx,lversion,newsurf,nlandt,nlandw,nwatcini,nwetsoil,rhosnow,forext,forhgt,forint,forpai,&
!$omp&  snowcovz,&
!$omp&  rinifor,rlue,rnbiocats,roffexp,roffpit,roffvel,&
!$omp&  sicecap,sicediff,snowcap,snowdiff,soilcap,soildiff,tau_soil,tau_veg,wsmax)

      end module landmod

!     ==================
!     FUNCTION NCOUNTSEA
!     ==================

      function ncountsea(plsm)
      use pumamod
      real :: plsm(NHOR)

      call mpsumval(plsm,NHOR,1,zsum)
!     THE FRACTION OF THE SPHERE, not the fraction of the cells. world-mt5.
!     zsum/NUGP is a cell count, and on a Gaussian grid a polar cell covers a
!     small fraction of the area of an equatorial one, so the counted percentage
!     is not the land fraction of the world and differs from it by an amount
!     that is a function of NLAT. The counts stay because n_sea_points is a
!     count; the percentage reported beside them is the area.
      call gpareamean(plsm,zlfrac)
      ilpo    = nint(zsum)
      ispo    = NUGP - ilpo
      ilperc  = nint(100.0 * zlfrac)
      isperc  = 100 - ilperc
  
      if (mypid == NROOT) then
         write(nud,'(a,i6,a,i6,a,i3,a)') &
              ' Land:',ilpo,' from',NUGP,' = ',ilperc,'% of area'
         write(nud,'(a,i6,a,i6,a,i3,a)') &
              ' Sea: ',ispo,' from',NUGP,' = ',isperc,'% of area'
      endif
      ncountsea = ispo

      return
      end

!     ==================
!     SUBROUTINE LANDINI
!     ==================

      subroutine landini
      use landmod
      use radmod
      use snowmaskmod
      use restartmod, only: nexcheck
!
!     initialize land surface
!
      namelist/landmod_nl/nlandt,nlandw,albland,dz0land,drhsland        &
     &                ,nlandwcol,nlsoilw,nlandwdrain,dsoilwf              &
     &                ,drhslow,nrhsexp,nlandwphase,dsoilwz                &
     &                ,dsnowalbmn,dsnowalbmx,dglacalbmn,dsnowalb        &
     &                ,dsmax,wsmax,drhsfull,dzglac,dztop,dsoilz         &
     &                ,rlue,co2conv,tau_veg,tau_soil                    &
     &                ,rnbiocats,nwetsoil,soilcap                       &
     &                ,albforest,forcovmx,forcovmn                      &
     &                ,forhgt,forpai,forext,forint                       &
     &                ,soildiff,sicediff,snowdiff,sicecap,snowcap       &
     &                ,rhosnow,roffvel,roffexp,roffpit                  &
     &                ,newsurf,rinifor,nwatcini,dwatcini,dgroundalb     &
     &                ,snowcovz
!
      dtclsoil(:) = tmelt
      dsoilt(:,:) = tmelt
      dsnowt(:)   = tmelt
      dtcl(:,:)   = tmelt
      dtclim(:)   = tmelt
      

      if (ndesert == 1 .and. nrestart == 0) then
         dwatcini = 0.0
         dwater(:) = 0.0
      endif
      
!       if (mypid==NROOT) then
!         
!        albsmax1 = dsnowalbmx(1)
!        albgmax1 = dsnowalbmx(1)
!        albsmax2 = dsnowalbmx(2)
!        albgmax2 = dsnowalbmx(2)
!        albsmaxf1 = 0.5*albsmax1
!        albsmaxf2 = 0.5*albsmax2
!        albsmin1 = dsnowalbmn(1)
!        albsmin2 = dsnowalbmn(2)
!        albgmin1 = dglacalbmn(1)
!        albgmin2 = dglacalbmn(2)
!        albsminf1 = 0.75*albsmaxf1
!        albsminf2 = 0.75*albsmaxf2
!        
!       endif
      
      if (mypid == NROOT) then
      open(12,file=landmod_namelist)
      read(12,landmod_nl)
      write(nud,'(/,"***********************************************")')
      write(nud,'("* LANDMOD ",a35," *")') trim(lversion)
      write(nud,'("***********************************************")')
      write(nud,'("* Namelist LANDMOD_NL from <",a17,"> *")') &
            landmod_namelist
      write(nud,'("***********************************************")')
      write(nud,landmod_nl)
      close(12)
      
      albsmax1 = dsnowalbmx(1)
      albgmax1 = dsnowalbmx(1)
      albsmax2 = dsnowalbmx(2)
      albgmax2 = dsnowalbmx(2)
      albsmin1 = dsnowalbmn(1)
      albsmin2 = dsnowalbmn(2)
      albgmin1 = dglacalbmn(1)
      albgmin2 = dglacalbmn(2)

!     Snow seen through a canopy: mix the band's exposed snow with the band's
!     canopy albedo at the masked fraction, rather than rescaling the band by a
!     ratio measured under a different spectrum. See the albforest block in the
!     module header for the derivation and for the Earth identity.

      albsmaxf1 = forcovmx*albforest(1) + (1.-forcovmx)*albsmax1
      albsmaxf2 = forcovmx*albforest(2) + (1.-forcovmx)*albsmax2
      albsminf1 = forcovmn*albforest(1) + (1.-forcovmn)*albsmin1
      albsminf2 = forcovmn*albforest(2) + (1.-forcovmn)*albsmin2

      write(nud,*) "Forested snow albedo, band 1: max",albsmaxf1, &
     &             " min",albsminf1
      write(nud,*) "Forested snow albedo, band 2: max",albsmaxf2, &
     &             " min",albsminf2

      endif

      if (wsmax < 0.0) wsmax = 0.0 ! Catch user error

      
      call mpbci(nlandt)
      call mpbci(nlandw)
      call mpbci(newsurf)
      call mpbci(nwetsoil)
      call mpbci(nwatcini)
      call mpbcr(albland)
      
      call mpbcr(albsmin)
      call mpbcr(albsmax)
      call mpbcr(albsminf)
      call mpbcr(albsmaxf)
      call mpbcr(albgmin)
      call mpbcr(albgmax)
      
      call mpbcr(albsmin1)
      call mpbcr(albsmax1)
      call mpbcr(albsminf1)
      call mpbcr(albsmaxf1)
      call mpbcr(albgmin1)
      call mpbcr(albgmax1)
      
      call mpbcr(albsmin2)
      call mpbcr(albsmax2)
      call mpbcr(albsminf2)
      call mpbcr(albsmaxf2)
      call mpbcr(albgmin2)
      call mpbcr(albgmax2)
      
      call mpbcr(dz0land)
      call mpbcr(drhsland)
      call mpbcr(drhsfull)
      call mpbci(nlandwcol)
      call mpbci(nlsoilw)
      call mpbci(nlandwdrain)
      call mpbci(nrhsexp)
      call mpbcr(drhslow)
      call mpbci(nlandwphase)
      call mpbcrn(dsoilwf,NLSOILWX)
      call mpbcrn(dsoilwz,NLSOILWX)

!     LSHY-3. Validate the scheme selection and normalise the layer shape, once
!     and here, so nothing downstream has to. A shape that does not sum to one
!     would silently change the column capacity away from dwmax, which is the
!     field the whole pedology loop feeds, so it is renormalised rather than
!     rejected -- and ONLY when the column is selected, so the default path
!     never touches dsoilwf and cannot be moved by rounding.
      if (nlandwcol == 1) then
       if (nlsoilw < 1 .or. nlsoilw > NLSOILWX) then
        if (mypid == NROOT) then
         write(nud,*)'*** nlsoilw = ',nlsoilw,' is outside 1 to ',NLSOILWX
        endif
        stop
       endif
       zsum = 0.
       do jlay=1,nlsoilw
        if (dsoilwf(jlay) < 0.0) then
         if (mypid == NROOT) write(nud,*)'*** dsoilwf must be non-negative'
         stop
        endif
        zsum = zsum + dsoilwf(jlay)
       enddo
       if (zsum <= 0.0) then
        if (mypid == NROOT) write(nud,*)'*** dsoilwf sums to zero over nlsoilw'
        stop
       endif
       do jlay=1,nlsoilw
        dsoilwf(jlay) = dsoilwf(jlay) / zsum
       enddo
       do jlay=nlsoilw+1,NLSOILWX
        dsoilwf(jlay) = 0.
       enddo
       if (mypid == NROOT) then
        write(nud,*)' *** LSHY-3: land liquid water on the layered column,'
        write(nud,*)' *** nlsoilw = ',nlsoilw,', lower boundary ',nlandwdrain
        write(nud,*)' *** layer shape ',(dsoilwf(jlay),jlay=1,nlsoilw)
       endif
      endif

!     LSHY-5. Phase is a property of a layer, so it needs the layered column.
      if (nlandwphase == 1) then
       if (nlandwcol /= 1) then
        if (mypid == NROOT) then
         write(nud,*)'*** nlandwphase = 1 needs nlandwcol = 1: the scalar'
         write(nud,*)'*** bucket has no layer to freeze'
        endif
        stop
       endif
       do jlay=1,nlsoilw
        if (dsoilwz(jlay) <= 0.0) then
         if (mypid == NROOT) then
          write(nud,*)'*** dsoilwz must be positive for every water layer in'
          write(nud,*)'*** use: phase needs a depth to find a temperature at'
         endif
         stop
        endif
       enddo
       if (mypid == NROOT) then
        write(nud,*)' *** LSHY-5: soil phase active on the water column,'
        write(nud,*)' *** layer thicknesses ',(dsoilwz(jlay),jlay=1,nlsoilw)
       endif
      endif
      call mpbcr(wsmax)
      call mpbcr(dwatcini)
      call mpbcr(dzglac)
      call mpbcr(dztop)
      call mpbcr(dsmax)
      call mpbcr(snowcovz)
      call mpbcr(rlue)
      call mpbcr(co2conv)
      call mpbcr(tau_veg)
      call mpbcr(tau_soil)
      call mpbcr(rinifor)
      call mpbcr(rnbiocats)
      call mpbcrn(dsoilz,NLSOIL)

      call mpbcrn(albforest,2)
      call mpbcr(forhgt)
      call mpbcr(forpai)
      call mpbcr(forext)
      call mpbcr(forint)
      call mpbcr(forcovmx)
      call mpbcr(forcovmn)
      call mpbcr(soildiff)
      call mpbcr(sicediff)
      call mpbcr(snowdiff)
      call mpbcr(soilcap)
      call mpbcr(sicecap)
      call mpbcr(snowcap)
      call mpbcr(rhosnow)
      call mpbcr(roffvel)
      call mpbcr(roffexp)
      call mpbcr(roffpit)

      call mpbcrn(dsnowalbmn,2)
      call mpbcrn(dsnowalbmx,2)
      call mpbcrn(dglacalbmn,2)
      call mpbcrn(dsnowalb,2)
      call mpbcrn(dgroundalb,2)
      

!     scale taus from years to seconds

      tau_veg  = tau_veg  * m_days_per_year * day_24hr
      tau_soil = tau_soil * m_days_per_year * day_24hr

!     ROOT PRINTS, EVERY THREAD STOPS -- the pattern the LSHY-3 and LSHY-5
!     refusals below already use. Both taus are broadcast above. world-0ihs.
      if (tau_veg < 1.0 .or. tau_soil < 1.0) then
         if (mypid == NROOT) then
          write(nud,*)' *** error: tau_veg = ',tau_veg,'  tau_soil = ',tau_soil
         endif
         stop 1
       endif

!     preset

      if (nrestart == 0) then
         dforest(:)= rinifor
      endif

      if (nrestart == 0) then
!
!*       preset some fields
!
         dz0clim(:)  = dz0land
         dwmax(:)    = wsmax
         dalbcl(:,:) = albland
         dalbcl1(:,:) = dgroundalb(1)
         dalbcl2(:,:) = dgroundalb(2)
!
!*       read surface parameters
!
         call mpsurfgp('doro'    ,doro    ,NHOR,1)
         call mpsurfgp('dls'     ,dls     ,NHOR,1)
         call mpsurfgp('dz0clim' ,dz0clim ,NHOR,1)
         call mpsurfgp('dz0climo',dz0climo,NHOR,1)
         call mpsurfgp('dglac'   ,dglac   ,NHOR,1)
         call mpsurfgp('dforest' ,dforest ,NHOR,1)
         call mpsurfgp('dwmax'   ,dwmax   ,NHOR,1)
         call mpsurfgp('dtclsoil',dtclsoil,NHOR,1)
   
         call mpsurfgp('dtcl',dtcl,NHOR,14)
         call mpsurfgp('dwcl',dwcl,NHOR,14)
         call mpsurfgp('dalbcl',dalbcl,NHOR,14)
         call mpsurfgp('dalbcl1',dalbcl1,NHOR,14)
         call mpsurfgp('dalbcl2',dalbcl2,NHOR,14)

!        make sure, that dwmax is positive

         where (dwmax(:) < 0.0) dwmax(:) = 0.0

         if (dwcl(1,1) < 0.0) then ! not in file
            do jm = 0 , 13
               dwcl(:,jm) = dwmax(:) * drhsfull * drhsland
            enddo
         endif

!        make sure, that land sea mask values are 0 or 1

         where (dls(:) > 0.5)
            dls(:) = 1.0
         elsewhere
            dls(:) = 0.0
         endwhere
         n_sea_points = ncountsea(dls)
!
!*       modify glacier mask according to dzglac
!

         if (dzglac > 0.0) then
            where (doro(:) > dzglac * ga) dglac(:)=1.0
         endif
!
!*       convert fractional glacier mask to binary mask
!
         where (dglac(:) > 0.5)
            dglac(:) = 1.0
         elsewhere
            dglac(:) = 0.0
         endwhere
!
!*    initialize soil
!
       call soilini
!
!*    initialize runoff (except for land only simulations)
!
       if (n_sea_points > 0) call roffini
!
!*    get new background albedo 
!
       call getalb
!
!*    set other surface variables
!
       do jhor=1,NHOR
        if(dls(jhor) > 0.0) then
         dtsm(jhor)=dts(jhor)
         dqs(jhor)=rdbrv*ra1s(dts(jhor))*EXP(ra2s(dts(jhor))*(dts(jhor)-tmelt)/ra4d(dts(jhor),ra4s(dts(jhor)))) &
     &            /psurf
         dqs(jhor)=dqs(jhor)/(1.-(1./rdbrv-1.)*dqs(jhor))
         dsnow(jhor)=dsnowz(jhor)
         if(dsnow(jhor) > 0.) then
          ! The canopy/snow mask. snowmaskmod returns the cover to blend each
          ! forested endmember at and the weight that moves it toward exposed
          ! snow, and reduces to dforest and zero at the inert defaults, so the
          ! four lines below are the same arithmetic they were.
          call snowcanopymask(dforest(jhor),dsnow(jhor)*1000./rhosnow,        &
          &                   forhgt,forpai,forext,forcovmx,forcovmn,forint, &
          &                   zfcovmx,zfcovmn,zfint,zkmx,zkmn)
          zsfmax =albsmaxf +zfint*(albsmax -albsmaxf)
          zsfmin =albsminf +zfint*(albsmin -albsminf)
          zsfmax1=albsmaxf1+zfint*(albsmax1-albsmaxf1)
          zsfmin1=albsminf1+zfint*(albsmin1-albsminf1)
          zsfmax2=albsmaxf2+zfint*(albsmax2-albsmaxf2)
          zsfmin2=albsminf2+zfint*(albsmin2-albsminf2)
          zalbmax=zfcovmx*zsfmax+(1.-zfcovmx)*albsmax
          zalbmin=zfcovmn*zsfmin+(1.-zfcovmn)*albsmin
          zdalb=(zalbmax-zalbmin)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow=MAX(zalbmin,MIN(zalbmax,zalbmax-zdalb))
          zalbmax1=zfcovmx*zsfmax1+(1.-zfcovmx)*albsmax1
          zalbmin1=zfcovmn*zsfmin1+(1.-zfcovmn)*albsmin1
          zdalb1=(zalbmax1-zalbmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow1=MAX(zalbmin1,MIN(zalbmax1,zalbmax1-zdalb1))
          zalbmax2=zfcovmx*zsfmax2+(1.-zfcovmx)*albsmax2
          zalbmin2=zfcovmn*zsfmin2+(1.-zfcovmn)*albsmin2
          zdalb2=(zalbmax2-zalbmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
          zalbsnow2=MAX(zalbmin2,MIN(zalbmax2,zalbmax2-zdalb2))
          dalb(jhor)=dalbclim(jhor)                                     &
     &        +(zalbsnow-dalbclim(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
          dsalb(1,jhor) = dalbclim1(jhor)                               &
     &        +(zalbsnow1-dalbclim1(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
          dsalb(2,jhor) = dalbclim2(jhor)                               &
     &        +(zalbsnow2-dalbclim2(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
          drhs(jhor)=1.
         else
          dalb(jhor)=dalbclim(jhor)
          dsalb(1,jhor)=dalbclim1(jhor)
          dsalb(2,jhor)=dalbclim2(jhor)
          if (dwmax(jhor) > 0.0)                                        &
     &    drhs(jhor)=land_wetness(dwatc(jhor),dwmax(jhor),drhsfull,      &
     &                            drhslow,nrhsexp)
         endif
         dz0(jhor)=dz0clim(jhor)
!
!        diagnostics: soil temps.
!
         dtsoil(jhor)=dsoilt(jhor,1)
         if(NLSOIL > 2) dtd2(jhor)=dsoilt(jhor,2)
         if(NLSOIL > 3) dtd3(jhor)=dsoilt(jhor,3)
         if(NLSOIL > 4) dtd4(jhor)=dsoilt(jhor,4)
         dtd5(jhor)=dsoilt(jhor,NLSOIL)

!
!*       modifications according to glacier mask
!

!        A cell that starts as glacier starts at the snow CAP, so with the cap
!        lifted (dsmax <= 0) it starts at zero here and glacierini, which runs
!        immediately after landini, raises it to glacelim. That is the glacier
!        module's own initial depth and the intended one; the line below is not
!        the only thing setting it.

         if(dglac(jhor) > 0.5) then
          dsnowz(jhor)=AMAX1(dsmax,0.)
          dsnow(jhor)=dsnowz(jhor)
          zdalb=(albgmax-albgmin)*(dts(jhor)-263.16)/(tmelt-263.16)
          zdalb1=(albgmax1-albgmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
          zdalb2=(albgmax2-albgmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
          dalb(jhor)=MAX(albgmin,MIN(albgmax,albgmax-zdalb))
          dsalb(1,jhor)=MAX(albgmin1,MIN(albgmax1,albgmax1-zdalb1))
          dsalb(2,jhor)=MAX(albgmin2,MIN(albgmax2,albgmax2-zdalb2))
          drhs(jhor)=1.0
         end if

!
!*    set PUMA temperature and moisture:
!

         dt(jhor,NLEP)=dts(jhor)
         dq(jhor,NLEP)=dqs(jhor)

        endif
       enddo

!      Same archival correction landstep makes; see the block at the end of
!      landstep for why dalb has to be rebuilt from dsalb rather than from the
!      non-spectral scalars.

       if (nstartemp == 1) then
        where(dls(:) > 0.0)                                             &
     &   dalb(:)=zsolars(1)*dsalb(1,:)+zsolars(2)*dsalb(2,:)
       endif

      else ! restart > 0

       call mpgetgp('dtsl'    ,dts     ,NHOR,     1)
       call mpgetgp('dtsm'    ,dtsm    ,NHOR,     1)
       call mpgetgp('dqs'     ,dqs     ,NHOR,     1)
       call mpgetgp('driver'  ,driver  ,NHOR,     1)
       call mpgetgp('duroff'  ,duroff  ,NHOR,     1)
       call mpgetgp('dvroff'  ,dvroff  ,NHOR,     1)
       call mpgetgp('darea'   ,darea   ,NHOR,     1)
       call mpgetgp('dwmax'   ,dwmax   ,NHOR,     1)
       call mpgetgp('dtcl'    ,dtcl    ,NHOR,    14)
       call mpgetgp('dwcl'    ,dwcl    ,NHOR,    14)
       call mpgetgp('dsnowt'  ,dsnowt  ,NHOR,     1)
       call mpgetgp('dsnowz'  ,dsnowz  ,NHOR,     1)
       call mpgetgp('dsoilt'  ,dsoilt  ,NHOR,NLSOIL)
       call mpgetgp('dz0clim' ,dz0clim ,NHOR,     1)
       call mpgetgp('dz0climo',dz0climo,NHOR,     1)
       call mpgetgp('dalbcl'  ,dalbcl  ,NHOR,    14)
       call mpgetgp('dalbcl1' ,dalbcl1 ,NHOR,    14)
       call mpgetgp('dalbcl2' ,dalbcl2 ,NHOR,    14)

!      The layered store and the drainage, under a LOWERED nexcheck: a restart
!      written before LSHY-3 existed carries neither record, and a run that
!      resumes on the default scheme must not stop for a field the default
!      scheme never reads. What it gets instead is the block below, which
!      rebuilds the layers from the scalar store the atmospheric restart has
!      already restored -- the same construction soilini uses on a cold start,
!      so the two paths agree.
       nexcheck = 0
       dwatcl(:,:) = -1.0
       dsoili(:,:) = 0.
       adrain(:) = 0.
       call mpgetgp('dwatcl'  ,dwatcl  ,NHOR,NLSOILWX)
       call mpgetgp('dsoili'  ,dsoili  ,NHOR,NLSOILWX)
       call mpgetgp('ddrain'  ,ddrain  ,NHOR,     1)
!      `adrain` needs no rebuild if the record is absent. It is an
!      accumulator over the output window, so a restart written before it
!      existed loses a partial window and nothing else; zeroed above, it
!      starts the window the target run is going to finish anyway.
       call mpgetgp('adrain'  ,adrain  ,NHOR,     1)
       nexcheck = 1
       if (ALL(dwatcl(:,:) < 0.0)) then
        dwatcl(:,:) = 0.
        dsoili(:,:) = 0.
        do jlay=1,nlsoilw
         where(dls(:) > 0.0) dwatcl(:,jlay)=dwatc(:)*dsoilwf(jlay)
        enddo
        ddrain(:) = 0.
        if (mypid == NROOT) then
         write(nud,*)' *** LSHY-3: no dwatcl in this restart; the layered'
         write(nud,*)' *** store was rebuilt from dwatc on the declared'
         write(nud,*)' *** layer shape. At nlsoilw = 1 that is exact.'
        endif
       endif

       n_sea_points = ncountsea(dls)
       
!     dwmax and dwcl are restored above; drhs and dwatc are restored by
!     read_atmos_restart, which plasim.f90 calls before surfini. Nothing on
!     this path recomputes any of them. Resetting the soil water capacity from
!     the namelist on an existing restart is what newsurf=2 below is for: it is
!     namelist-driven, and it leaves the per-cell dwmax and the 14-month dwcl
!     as fields rather than flattening them to scalars.

      endif ! if (restart == 0)

      if (newsurf == 1) then
         call mpsurfgp('dtcl',dtcl,NHOR,14)
         call mpsurfgp('dwcl',dwcl,NHOR,14)
      endif

      if (newsurf == 2) then ! preset some fields
         dwmax(:)    = wsmax
         dz0clim(:)  = dz0land
         dalbcl(:,:) = albland
         dalbcl1(:,:) = dgroundalb(1)
         dalbcl2(:,:) = dgroundalb(2)
         dwcl(:,:)   = wsmax * drhsfull * drhsland
      endif

      return
      end subroutine landini


!     ===================
!     SUBROUTINE LANDSTEP
!     ===================

      subroutine landstep
      use landmod
      use radmod, only: nstartemp, zsolars
      use snowmaskmod

!
!     get climatological values if t and/or w are non interactive
!

      if(nlandt==0 .or. nlandw==0 ) call gettcll

!
!     soil
!

      call soilstep

!
!     runoff
!

      call roffstep
!
!     get new background albedo
!
      call getalb
!
!     set surface variables
!
      do jhor=1,NHOR
       if(dls(jhor) > 0.0) then
        dtsm(jhor)=dts(jhor)
        dqs(jhor)=rdbrv*ra1s(dts(jhor))*EXP(ra2s(dts(jhor))*(dts(jhor)-tmelt)/ra4d(dts(jhor),ra4s(dts(jhor))))  &
     &           /dp(jhor)
        dqs(jhor)=dqs(jhor)/(1.-(1./rdbrv-1.)*dqs(jhor))
        dsnow(jhor)=dsnowz(jhor)
        if(dsnow(jhor) > 0.) then
         ! The canopy/snow mask. snowmaskmod returns the cover to blend each
         ! forested endmember at and the weight that moves it toward exposed
         ! snow, and reduces to dforest and zero at the inert defaults, so the
         ! four lines below are the same arithmetic they were.
         call snowcanopymask(dforest(jhor),dsnow(jhor)*1000./rhosnow,        &
         &                   forhgt,forpai,forext,forcovmx,forcovmn,forint, &
         &                   zfcovmx,zfcovmn,zfint,zkmx,zkmn)
         zsfmax =albsmaxf +zfint*(albsmax -albsmaxf)
         zsfmin =albsminf +zfint*(albsmin -albsminf)
         zsfmax1=albsmaxf1+zfint*(albsmax1-albsmaxf1)
         zsfmin1=albsminf1+zfint*(albsmin1-albsminf1)
         zsfmax2=albsmaxf2+zfint*(albsmax2-albsmaxf2)
         zsfmin2=albsminf2+zfint*(albsmin2-albsminf2)
         zalbmax=zfcovmx*zsfmax+(1.-zfcovmx)*albsmax
         zalbmin=zfcovmn*zsfmin+(1.-zfcovmn)*albsmin
         zdalb=(zalbmax-zalbmin)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow=MAX(zalbmin,MIN(zalbmax,zalbmax-zdalb))
         zalbmax1=zfcovmx*zsfmax1+(1.-zfcovmx)*albsmax1
         zalbmin1=zfcovmn*zsfmin1+(1.-zfcovmn)*albsmin1
         zdalb1=(zalbmax1-zalbmin1)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow1=MAX(zalbmin1,MIN(zalbmax1,zalbmax1-zdalb1))
         zalbmax2=zfcovmx*zsfmax2+(1.-zfcovmx)*albsmax2
         zalbmin2=zfcovmn*zsfmin2+(1.-zfcovmn)*albsmin2
         zdalb2=(zalbmax2-zalbmin2)*(dts(jhor)-263.16)/(tmelt-263.16)
         zalbsnow2=MAX(zalbmin2,MIN(zalbmax2,zalbmax2-zdalb2))
         dalb(jhor)=dalbclim(jhor)                                     &
     &       +(zalbsnow-dalbclim(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
         dsalb(1,jhor) = dalbclim1(jhor)                               &
     &       +(zalbsnow1-dalbclim1(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
         dsalb(2,jhor) = dalbclim2(jhor)                               &
     &       +(zalbsnow2-dalbclim2(jhor))*dsnow(jhor)/(dsnow(jhor)+snowcovz)
         drhs(jhor)=1.
        else
         dalb(jhor)=dalbclim(jhor)
         dsalb(1,jhor)=dalbclim1(jhor)
         dsalb(2,jhor)=dalbclim2(jhor)
         if (dwmax(jhor) > 0.0)                                        &
     &   drhs(jhor)=land_wetness(dwatc(jhor),dwmax(jhor),drhsfull,       &
     &                           drhslow,nrhsexp)
        endif

!
!     diagnostics: soil temps.
!

        dtsoil(jhor)=dsoilt(jhor,1)
        if(NLSOIL > 2) dtd2(jhor)=dsoilt(jhor,2)
        if(NLSOIL > 3) dtd3(jhor)=dsoilt(jhor,3)
        if(NLSOIL > 4) dtd4(jhor)=dsoilt(jhor,4)
        dtd5(jhor)=dsoilt(jhor,NLSOIL)

!
!*    set PUMA temperature and moisture:
!

        dt(jhor,NLEP)=dts(jhor)
        dq(jhor,NLEP)=dqs(jhor)

       end if
      enddo

!
!*    vegetation
!

      if (nveg > 0) call vegstep

!
!*     modifications according to glacier mask
!

      where(dglac(:) > 0.5 .and. dls(:) > 0.0)
       dalb(:)=MAX(albgmin,MIN(albgmax                                  &
     &  ,albgmax-(albgmax-albgmin)*(dts(:)-263.16)/(tmelt-263.16)))
       dsalb(1,:)=MAX(albgmin1,MIN(albgmax1                                  &
     &  ,albgmax1-(albgmax1-albgmin1)*(dts(:)-263.16)/(tmelt-263.16)))
       dsalb(2,:)=MAX(albgmin2,MIN(albgmax2                                  &
     &  ,albgmax2-(albgmax2-albgmin2)*(dts(:)-263.16)/(tmelt-263.16)))
       drhs(:)=1.0
      end where

!
!*    archive the albedo the radiation actually used
!
!     At nstartemp = 1 the shortwave reads dsalb and nothing else, and radmod
!     writes the flux-weighted combination of the two bands into dalb as the
!     diagnostic of what it used. radstep runs BEFORE surfstep, so everything
!     above has just overwritten that diagnostic with a value built from the
!     six non-spectral scalars, which landini never re-derives from the
!     spectrum and which are therefore still Earth-Sun broadband. dalb is what
!     outmod accumulates and writes as code 175, so 175 disagreed with the
!     radiation over every snow, forest-snow and glacier cell.
!
!     Over land radmod leaves dsalb untouched (its ocean direct-beam branch is
!     gated on 1-dls), so this reproduces radmod's own dalb exactly rather than
!     approximating it. At nstartemp = 0 the shortwave reads dalb itself and
!     the block above is the value it will use, so nothing happens here.

      if (nstartemp == 1) then
       where(dls(:) > 0.0)                                              &
     &  dalb(:)=zsolars(1)*dsalb(1,:)+zsolars(2)*dsalb(2,:)
      endif

      return
      end subroutine landstep

!     ===================
!     SUBROUTINE LANDSTOP
!     ===================

      subroutine landstop
      use landmod

      if (mypid == NROOT) then
         call put_restart_integer('nlsoil',NLSOIL)
      endif

      call mpputgp('dtsl'    ,dts     ,NHOR, 1)
      call mpputgp('dtsm'    ,dtsm    ,NHOR, 1)
      call mpputgp('dqs'     ,dqs     ,NHOR, 1)
      call mpputgp('driver'  ,driver  ,NHOR, 1)
      call mpputgp('duroff'  ,duroff  ,NHOR, 1)
      call mpputgp('dvroff'  ,dvroff  ,NHOR, 1)
      call mpputgp('darea'   ,darea   ,NHOR, 1)
      call mpputgp('dwmax'   ,dwmax   ,NHOR, 1)
      call mpputgp('dtcl'    ,dtcl    ,NHOR,14)
      call mpputgp('dwcl'    ,dwcl    ,NHOR,14)
      call mpputgp('dsnowt'  ,dsnowt  ,NHOR, 1)
      call mpputgp('dsnowz'  ,dsnowz  ,NHOR, 1)
      call mpputgp('dsoilt'  ,dsoilt  ,NHOR,NLSOIL)
      call mpputgp('dwatcl'  ,dwatcl  ,NHOR,NLSOILWX)
      call mpputgp('dsoili'  ,dsoili  ,NHOR,NLSOILWX)
      call mpputgp('ddrain'  ,ddrain  ,NHOR, 1)
      call mpputgp('adrain'  ,adrain  ,NHOR, 1)
      call mpputgp('dz0clim' ,dz0clim ,NHOR, 1)
      call mpputgp('dz0climo',dz0climo,NHOR, 1)
      call mpputgp('dalbcl'  ,dalbcl  ,NHOR,14)
      call mpputgp('dalbcl1' ,dalbcl1 ,NHOR,14)
      call mpputgp('dalbcl2' ,dalbcl2 ,NHOR,14)
      return
      end subroutine landstop

!     ================
!     SUBROUTINE TANDS
!     ================

      subroutine tands
      use landmod
!
      parameter(zsnowmax=1.)
      parameter(ztop=0.1)
!
      real zsnowz(NHOR)       ! new snow depth
      real zhfls(NHOR)        ! heatflux from soil
      real zhfla(NHOR)        ! heatflux from atmosphere
      real zhflm(NHOR)        ! heatflux used by snowmelt
      real zdsnowz(NHOR)      ! snow depth tendency
      real zcap(NHOR,NLSOIL)  ! heat capacity  of soil layers
      real zdiff(NHOR,NLSOIL) ! thermal conductivity of soil layers
      real zsoilz(NHOR,NLSOIL)! soil layer thicknesses
      real zcap1(NHOR)        ! heat capacity (upper soil layer)
      real zdiff1(NHOR)       ! thermal conductivity (upper soil layer)
      real zsoilz1(NHOR)      ! layer thickness (upper soil layer)
      real zctop(NHOR)        ! heat capacity (top soil/snow layer)
      real zsntop(NHOR)       ! snow depth (top soil/snow layer)
      real zztop(NHOR)        ! depth of the uppermost layer
!
!     debug
!
      integer kmaxl(1),kminl(1)
      real,allocatable :: zfpr1(:)
      real,allocatable :: zfpr2(:)
      real,allocatable :: zfpr3(:)
      real,allocatable :: zfpr4(:)
      real,allocatable :: zfpr5(:)
      real,allocatable :: zfpr6(:)
      real,allocatable :: zfpr7(:)
      real,allocatable :: zfpr8(:)
      real,allocatable :: zfpr9(:)
      real,allocatable :: zfprl(:,:)
!
!     preset
!
      dwater(:)=0.
      dsmelt(:)=0.
      dsndch(:)=dsnowz(:)
      zdsnowz(:)=0.
      zhfls(:)=0.
      zhflm(:)=0.
      zztop(:)=dztop
!
      do jlev=1,NLSOIL
       where(dls(:) > 0.0)
        zcap(:,jlev)=sicecap*dglac(:)+soilcap*(1.-dglac(:))
        zdiff(:,jlev)=sicediff*dglac(:)+soildiff*(1.-dglac(:))
        zsoilz(:,jlev)=dsoilz(jlev)
       endwhere
      enddo
!
      where(dsnowz(:) == 0.0) dsnowt(:)=tmelt
!
!     copy heatflux
!
      zhfla(:)=dshfl(:)+dlhfl(:)+dflux(:,NLEP)
!
!     precip destribution
!
      do jhor=1,NHOR

!
!     surface water
!

       if(dls(jhor) > 0.0) then
        if(dprs(jhor) > 0.) zdsnowz(jhor)=dprs(jhor)
        if(dsnowz(jhor) > 0.) zdsnowz(jhor)=zdsnowz(jhor)+devap(jhor)
        zsnowz(jhor)=AMAX1(0.,dsnowz(jhor)+zdsnowz(jhor)*deltsec)
        zdsnowz(jhor)=(zsnowz(jhor)-dsnowz(jhor))/deltsec
        dwater(jhor)=devap(jhor)+dprl(jhor)+dprc(jhor)-zdsnowz(jhor)
        dsnowz(jhor)=zsnowz(jhor)
       end if
      enddo
!
!     surface temperature and heat flux into the soil
!
      where(dls(:) > 0.0)
!
!     new properties of the top dztop meters (mixed soil/snow layer):
!     1: snow (m of ztop) 2: heat capacity
!
       zsntop(:)=AMIN1(dsnowz(:)*1000./rhosnow,zztop(:))
       zctop(:)=(snowcap*zsntop(:)+zcap(:,1)*(zztop(:)-zsntop(:)))      &
     &         /zztop(:)
!
!     new properties of the uppermost soil layer (mixed soil/snow)
!     1: snow depth (m) 2: total thickness 3: conductivity
!
       zsnowz(:)=AMIN1(dsnowz(:)*1000./rhosnow,zsnowmax)-zsntop(:)
       zsoilz1(:)=zsoilz(:,1)+zsnowz(:)
       zdiff1(:)=zdiff(:,1)*snowdiff*zsoilz1(:)                         &
     &          /(snowdiff*zsoilz(:,1)+zdiff(:,1)*zsnowz(:))
       zcap1(:)=(snowcap*zsnowz(:)+zcap(:,1)*zsoilz(:,1))               &
     &         /zsoilz1(:)
!
!     new surface temp. implicit w.r.t. conductive heat flux
!
       dts(:)=(zctop(:)*zztop(:)*dtsm(:)/deltsec+zhfla(:)               &
              +2.*zdiff1(:)*dsoilt(:,1)/zsoilz1(:))                     &
             /(zctop(:)*zztop(:)/deltsec+2.*zdiff1(:)/zsoilz1(:)) 
!
!     heat flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
!
      endwhere
!
!     snow melt
!
      where(dls(:) > 0.0 .and. dsnowz(:) > 0. .and. dts(:) > tmelt )
       dts(:)=tmelt
       zhflm(:)=AMAX1(0.,zhfla(:)                                       &
     &                  -zctop(:)*zztop(:)*(dts(:)-dtsm(:))/deltsec     &
     &                  +2.*zdiff1(:)*(dsoilt(:,1)-dts(:))/zsoilz1(:))  
       zsnowz(:)=AMAX1(0.,dsnowz(:)-zhflm(:)*deltsec/((ALS-ALV)*1000.))
       zdsnowz(:)=(zsnowz(:)-dsnowz(:))/deltsec
       zhflm(:)=-zdsnowz(:)*1000.*(ALS-ALV)
!
!     new snow depth (h2o equiv. and snow melt diagnostic)
!
       dsnowz(:)=zsnowz(:)
       dsmelt(:)=-zdsnowz(:)
!
!     heat flux and water flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
       dwater(:)=dwater(:)-zdsnowz(:)
!
!     new properties of the top dztop meters (mixed soil/snow layer):
!     1: snow (m of dztop) 2: heat capacity
!
       zsntop(:)=AMIN1(dsnowz(:)*1000./rhosnow,zztop(:))
       zctop(:)=(snowcap*zsntop(:)+zcap(:,1)*(zztop(:)-zsntop(:)))      &
     &         /zztop(:)
!
!     new properties of the uppermost soil layer (mixed soil/snow)
!     1: snow depth (m) 2: total thickness 3: conductivity 4: heat capacity
!
       zsnowz(:)=AMIN1(dsnowz(:)*1000./rhosnow,zsnowmax)-zsntop(:)
       zsoilz1(:)=zsoilz(:,1)+zsnowz(:)
       zdiff1(:)=zdiff(:,1)*snowdiff*zsoilz1(:)                         &
     &          /(snowdiff*zsoilz(:,1)+zdiff(:,1)*zsnowz(:))
       zcap1(:)=(snowcap*zsnowz(:)+zcap(:,1)*zsoilz(:,1))               &
     &         /zsoilz1(:)
!
      endwhere
!     write(nud,*) 'landmod 881 zsnowz = ',sum(zsnowz(:))
!
!     correct surface temp and flux where snow is totally melted
!
      where(dls(:) > 0.0 .and. dsnowz(:) == 0. .and. zhflm(:) > 0.)
!
!     new surface (as above)
!
       dts(:)=(zctop(:)*zztop(:)*dtsm(:)/deltsec+zhfla(:)-zhflm(:)      &
     &        +2.*zdiff1(:)*dsoilt(:,1)/zsoilz1(:))                     &
     &       /(zctop(:)*zztop(:)/deltsec+2.*zdiff1(:)/zsoilz1(:)) 
!
!     heat flux into the soil
!
       zhfls(:)=2.*zdiff1(:)/zsoilz1(:)*(dts(:)-dsoilt(:,1))
!
      endwhere
!
!     set new soil properties for the uppermost blended soil/snow layer
!

      where(dls(:) > 0.0)
       zsoilz(:,1)=zsoilz1(:)
       zcap(:,1)=zcap1(:)
       zdiff(:,1)=zdiff1(:)
      endwhere

!
!     deep soil temperatures
!

      call mktsoil(zhfls,zsoilz,zcap,zdiff)

!
!     snow temperature (at the moment set to ts)
!

      where(dls(:) > 0.0 .and. dsnowz(:) > 0.) dsnowt(:)=dts(:)
!
!     limit snow depth to maximum (if switched on)
!     (heat s conserved by artifical cooling of ts (dsnowt)
!      water is conserved)
!
      if(dsmax > 0.) then
       where(dls(:) > 0.0 .and. dsnowz(:) > dsmax)
        zdsnowz(:)=(dsmax-dsnowz(:))/deltsec
!
!     new snow depth (h2o equiv. and snow melt diagnostic)
!
        dsnowz(:)=dsmax
!
!     water flux into the soil
!
        dwater(:)=dwater(:)-zdsnowz(:)
!
!     cool ts (dsnowt)
!
        dts(:)=dts(:)+zdsnowz(:)*1000.*(ALS-ALV)*deltsec/zctop(:)/zztop(:)
!
!       diagnose the lost snow as snow melt
!
        dsmelt(:)=dsmelt(:)-zdsnowz(:)
       endwhere
      endif
!
!     diagnostic: snow depth change
!
      where(dls(:) > 0.) dsndch(:)=(dsnowz(:)-dsndch(:))!/deltsec
!
!     dbug output
!

      if(nprint == 1) then
       allocate(zfpr1(NLON*NLAT))
       allocate(zfpr2(NLON*NLAT))
       allocate(zfpr3(NLON*NLAT))
       allocate(zfpr4(NLON*NLAT))
       allocate(zfprl(NLON*NLAT,NLSOIL))
       call mpgagp(zfpr1,dprs,1)
       call mpgagp(zfpr2,dts,1)
       call mpgagp(zfpr3,dsnowz,1)
       call mpgagp(zfpr4,dls,1)
       call mpgagp(zfprl,dsoilt,NLSOIL)
       if(mypid == NROOT) then
        write(nud,*)'Land Surface Global Diagnostic:'
        zzmax=MAXVAL(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        zzmin=MINVAL(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        kminl=MINLOC(zfpr2(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX TS = ',zzmax,' NHOR= ',kmaxl(1)
        write(nud,*)'MIN TS = ',zzmin,' NHOR= ',kminl(1)
        zzmax=MAXVAL(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        zzmin=MINVAL(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        kminl=MINLOC(zfpr3(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX ZSNOW = ',zzmax,' NHOR= ',kmaxl(1)
        write(nud,*)'MIN ZSNOW = ',zzmin,' NHOR= ',kminl(1)
        do jlev=1,NLSOIL
         zzmax=MAXVAL(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         kmaxl=MAXLOC(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         zzmin=MINVAL(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         kminl=MINLOC(zfprl(:,jlev),MASK=(zfpr4(:) > 0.5))
         write(nud,*)'MAX TSOIL L= ',jlev,' = ',zzmax,' NHOR= ',kmaxl(1)
         write(nud,*)'MIN TSOIL L= ',jlev,' = ',zzmin,' NHOR= ',kminl(1)
        enddo
        zzmax=MAXVAL(zfpr1(:),MASK=(zfpr4(:) > 0.5))
        kmaxl=MAXLOC(zfpr1(:),MASK=(zfpr4(:) > 0.5))
        write(nud,*)'MAX PRS = ',zzmax,' NHOR= ',kmaxl(1)
       endif
       deallocate(zfpr1)
       deallocate(zfpr2)
       deallocate(zfpr3)
       deallocate(zfpr4)
       deallocate(zfprl)
      endif
      if(nprint == 2) then
       allocate(zfpr1(NLON*NLAT))
       allocate(zfpr2(NLON*NLAT))
       allocate(zfpr3(NLON*NLAT))
       allocate(zfpr4(NLON*NLAT))
       allocate(zfpr5(NLON*NLAT))
       allocate(zfpr6(NLON*NLAT))
       allocate(zfpr7(NLON*NLAT))
       allocate(zfpr8(NLON*NLAT))
       allocate(zfpr9(NLON*NLAT))
       allocate(zfprl(NLON*NLAT,NLSOIL))
       call mpgagp(zfpr1,dprs,1)
       call mpgagp(zfpr2,dts,1)
       call mpgagp(zfpr3,dsnowz,1)
       call mpgagp(zfpr4,zhfla,1)
       call mpgagp(zfpr5,zhflm,1)
       call mpgagp(zfpr6,zhfls,1)
       call mpgagp(zfpr8,dwater,1)
       call mpgagp(zfpr9,dsmelt,1)
       call mpgagp(zfprl,dsoilt,NLSOIL)
       if(mypid == NROOT) then
        write(nud,*)'Land Surface Local Diagnostic:'
        do jlev=1,NLSOIL
         write(nud,*)'TSOIL L= ',jlev,' = ',zfprl(NPRHOR,jlev)
        enddo
        write(nud,*)'TS             = ',zfpr2(NPRHOR)
        write(nud,*)'ZSNOW          = ',zfpr3(NPRHOR)
        write(nud,*)'HFLA           = ',zfpr4(NPRHOR)
        write(nud,*)'HFLX MELTING   = ',zfpr5(NPRHOR)
        write(nud,*)'HFLX TO SOIL   = ',zfpr6(NPRHOR)
        write(nud,*)'PRS (mm/d)     = ',zfpr1(NPRHOR)*1000.*deltsec*mtspd
        write(nud,*)'DWATER (mm/d)  = ',zfpr8(NPRHOR)*1000.*deltsec*mtspd
        write(nud,*)'SNOWMELT (mm/d)= ',zfpr9(NPRHOR)*1000.*deltsec*mtspd
       endif
       deallocate(zfpr1)
       deallocate(zfpr2)
       deallocate(zfpr3)
       deallocate(zfpr4)
       deallocate(zfpr5)
       deallocate(zfpr6)
       deallocate(zfpr7)
       deallocate(zfpr8)
       deallocate(zfpr9)
       deallocate(zfprl)
      endif
!
!     entropy diagnostics
!
!
      return
      end subroutine tands

!     ==================
!     SUBROUTINE MKTSOIL
!     ==================

      subroutine mktsoil(pftop,psoilz,pcap,pdiff)
      use landmod
!
!     input
!
      real pftop(NHOR)         ! surface heat flux
      real psoilz(NHOR,NLSOIL) ! layer thickness
      real pcap(NHOR,NLSOIL)   ! heat capacities
      real pdiff(NHOR,NLSOIL)  ! heat conductivities
!
!     local
!
      real ztn(NHOR,NLSOIL)    ! new temperatures
      real zebs(NHOR,NLSOIL)
      real zcap(NHOR,NLSOIL)
      real zdiff(NHOR,NLSOIL-1)
      real ztold(NHOR,NLSOIL)
!
!
!     implicit scheme for soiltemp
!
!     a) deep layer elimination (zero flux at bottom)
!

      jlev=NLSOIL
      jlem=NLSOIL-1
      where(dls(:) > 0.0)
       zdiff(:,jlem)=2.*pdiff(:,jlev)*pdiff(:,jlem)                     &
     &              /(pdiff(:,jlev)*psoilz(:,jlem)                      &
     &               +pdiff(:,jlem)*psoilz(:,jlev))
       zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
       zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlem))
       ztn(:,jlev)=zcap(:,jlev)*dsoilt(:,jlev)*zebs(:,jlev)
      endwhere

!
!     middle layer elimination
!

      do jlev=NLSOIL-1,2,-1
       jlep=jlev+1
       jlem=jlev-1
       where(dls(:) > 0.0)
        zdiff(:,jlem)=2.*pdiff(:,jlev)*pdiff(:,jlem)                    &
     &                /(pdiff(:,jlev)*psoilz(:,jlem)                    &
     &                 +pdiff(:,jlem)*psoilz(:,jlev))
        zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
        zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlem)                     &
     &                  +zdiff(:,jlev)*(1.-zdiff(:,jlev)*zebs(:,jlep)))
        ztn(:,jlev)=zebs(:,jlev)*(zcap(:,jlev)*dsoilt(:,jlev)           &
     &                           +zdiff(:,jlev)*ztn(:,jlep))
       endwhere
      enddo

!
!     top layer elimination
!

       jlev=1
       jlep=2
       where(dls(:) > 0.0)
        zcap(:,jlev)=pcap(:,jlev)*psoilz(:,jlev)/deltsec
        zebs(:,jlev)=1./(zcap(:,jlev)+zdiff(:,jlev)                     &
     &                               *(1.-zdiff(:,jlev)*zebs(:,jlep)))
        ztn(:,jlev)=zebs(:,jlev)*(zcap(:,jlev)*dsoilt(:,jlev)           &
     &                           +zdiff(:,jlev)*ztn(:,jlep)+pftop(:))
       endwhere

!
!     back-substitution
!

      do jlev=2,NLSOIL
       jlem=jlev-1
       where(dls(:) > 0.0)
        ztn(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*zdiff(:,jlem)*ztn(:,jlem)
       endwhere
      enddo

!
!     new temperatures
!

      do jlev=1,NLSOIL
       where(dls(:) > 0.0)
        dsoilt(:,jlev)=ztn(:,jlev)
       endwhere

!
!      do not allow warming of permanet glaciers
!

       where(dls(:) > 0.0 .and. dglac(:) > 0.5)
        dsoilt(:,jlev)=AMIN1(dsoilt(:,jlev),tmelt)
       endwhere
      enddo
!
!     entropy diagnostics
!
!
      return
      end subroutine mktsoil

!     ================
!     SUBROUTINE WANDR
!     ================

      subroutine wandr
      use landmod
!
!     calculate soil water and runoff
!
!     THE SCALAR BUCKET, and the default. The three lines this used to be are
!     now `bucket_step` in landcolumn.f90, unchanged in value and in order, so
!     that the layered column can be checked against them by a standalone
!     program rather than by inspection. The `where` becomes a loop because a
!     kernel takes scalars; the operations per cell are the same operations on
!     the same values, so the result is the same bits.
!
      drunoff(:)=0.
      ddrain(:)=0.
      do jhor=1,NHOR
       if (dls(jhor) > 0.0) then
        call bucket_step(dwatc(jhor),dwmax(jhor),dwater(jhor),deltsec,   &
     &                   zwnew,zroff)
        dwatc(jhor)   = zwnew
        drunoff(jhor) = zroff
       endif
      enddo
!
      return
      end subroutine wandr

!     ===================
!     SUBROUTINE WANDRCOL
!     ===================

      subroutine wandrcol
      use landmod
!
!     The layered liquid column, LSHY-3's second registered hypothesis.
!
!     `dwatc` stays the column total, because fluxmod's evaporation limiter,
!     simba's water stress factor, aeromod and outmod all read it and none of
!     them knows about layers. At nlsoilw = 1 the sum is over one element and
!     is exact, which is the other half of the reduction.
!
      real :: zcap(NLSOILWX)
      real :: zwl(NLSOILWX)
      real :: zil(NLSOILWX)
      real :: zwn(NLSOILWX)
!
      drunoff(:)=0.
      ddrain(:)=0.
      do jhor=1,NHOR
       if (dls(jhor) > 0.0) then
        zcap(:) = 0.
        zwl(:)  = 0.
        zil(:)  = 0.
        do jlay=1,nlsoilw
         zcap(jlay) = dwmax(jhor) * dsoilwf(jlay)
         zwl(jlay)  = dwatcl(jhor,jlay)
         zil(jlay)  = dsoili(jhor,jlay)
        enddo
        call column_step(nlsoilw,zwl,zil,zcap,dwater(jhor),deltsec,      &
     &                   nlandwdrain,zwn,zroff,zdrn)
        zsum = 0.
        do jlay=1,nlsoilw
         dwatcl(jhor,jlay) = zwn(jlay)
         zsum = zsum + zwn(jlay)
        enddo
        dwatc(jhor)   = zsum
        drunoff(jhor) = zroff
        ddrain(jhor)  = zdrn
       endif
      enddo
!
      return
      end subroutine wandrcol

!     ===================
!     SUBROUTINE LANDPHASE
!     ===================

      subroutine landphase
      use landmod
!
!     LSHY-5. Freeze and thaw the liquid in each water layer against the soil
!     temperature at its depth, conserving water and energy.
!
!     THE TWO COLUMNS SHARE NO BOUNDARY. The soil temperature layers are
!     `dsoilz`, five of them reaching 12.4 m, and the water layers are
!     `dsoilwz`, which the model did not have until this row. The mapping is
!     declared and is by MIDPOINT: a water layer takes the temperature of
!     whichever temperature layer contains its centre. That is a choice and it
!     is written down rather than emergent; the alternative -- making one
!     column a sub-partition of the other -- changes the soil heat solver's
!     geometry and is a larger decision than this row.
!
!     OPERATOR SPLIT. `tands` has already solved the soil temperatures for this
!     step and `wandrcol` has already moved the liquid. This adjusts both for
!     the phase change afterwards, which is what makes the exchange exactly
!     conservative: the latent heat leaves one store and arrives in the other
!     as the same number, with no solver between them to round it.
!
      real :: zwl(NLSOILWX)
      real :: zil(NLSOILWX)
      real :: zztop(NLSOILWX)
      real :: ztop, zmid
      integer :: itlay(NLSOILWX)
!
!     Which temperature layer each water layer's midpoint falls in. A property
!     of the two declared geometries and not of the cell, so it is worked out
!     once rather than per cell.
!
      ztop = 0.
      do jlay=1,nlsoilw
       zmid = ztop + 0.5*dsoilwz(jlay)
       zztop(jlay) = ztop
       ztop = ztop + dsoilwz(jlay)
       zbot = 0.
       itlay(jlay) = NLSOIL
       do jt=1,NLSOIL
        zbot = zbot + dsoilz(jt)
        if (zmid <= zbot) then
         itlay(jlay) = jt
         exit
        endif
       enddo
      enddo
!
      do jhor=1,NHOR
       if (dls(jhor) > 0.0) then
        zcapv = sicecap*dglac(jhor) + soilcap*(1.-dglac(jhor))
        do jlay=1,nlsoilw
         zwl(jlay) = dwatcl(jhor,jlay)
         zil(jlay) = dsoili(jhor,jlay)
        enddo
        zsum = 0.
        do jlay=1,nlsoilw
         it = itlay(jlay)
         call phase_step(zwl(jlay),zil(jlay),dsoilt(jhor,it),zcapv,      &
     &                   dsoilwz(jlay),tmelt,als-alv,1000.,             &
     &                   zliqn,zicen,ztemn)
         dwatcl(jhor,jlay) = zliqn
         dsoili(jhor,jlay) = zicen
         dsoilt(jhor,it)   = ztemn
         zsum = zsum + zliqn
        enddo
!
!       dwatc is the LIQUID total, not the water total. Frozen pore water is
!       not available to evaporate, and dwatc is what fluxmod's evaporation
!       limiter and simba's water stress read.
!
        dwatc(jhor) = zsum
       endif
      enddo
!
      return
      end subroutine landphase

!     ==================
!     SUBROUTINE SOILINI
!     ==================

      subroutine soilini
      use landmod
!
      if(nrestart == 0) then
!
!     initialize snow and temperatures
!
       where(dls(:) > 0.0)
        dsnowz(:)=0.
        dsnowt(:)=tmelt
        dts(:)=dtclsoil(:)
       endwhere
       do jlev=1,NLSOIL
        where(dls(:) > 0.0) dsoilt(:,jlev)=dts(:)
       enddo
!
!     initialize soil water and runoff
!
       if (nwatcini > 0) then
        where(dls(:) > 0.0)
         dwatc(:)=dwatcini
        endwhere
       else
        where(dls(:) > 0.0)
         dwatc(:)=drhsland*drhsfull*dwmax(:)
        endwhere
       endif
       drunoff(:)=0.
       ddrain(:)=0.
!
!      The layered store starts holding the same water as the scalar one,
!      distributed by the declared layer shape. At nlsoilw = 1 that puts all
!      of it in layer 1 and dwatcl(:,1) is dwatc exactly, so a cold start and
!      a scheme switch agree.
!
       dwatcl(:,:)=0.
       dsoili(:,:)=0.
       do jlay=1,nlsoilw
        where(dls(:) > 0.0) dwatcl(:,jlay)=dwatc(:)*dsoilwf(jlay)
       enddo
!
      endif
!
      return
      end subroutine soilini

!     ===================
!     SUBROUTINE SOILSTEP
!     ===================

      subroutine soilstep
      use landmod
!
!     calculate snow and temperatures
!
      if(nlandt==1) then
       call tands
      else
       where(dls(:) > 0.0)
        dwater(:)=devap(:)+dprl(:)+dprc(:)
        dts(:)=dtclim(:)
       endwhere
      end if
!
!     calculate soil water and runoff
!
      if(nlandw==1) then
       if(nlandwcol==1) then
        call wandrcol
        if(nlandwphase==1) call landphase
       else
        call wandr
       endif
      else
       drunoff(:)=0.
       ddrain(:)=0.
       where(dls(:) > 0.)
        drunoff(:)=AMAX1(0.,dwater(:))
        dwatc(:)=dwclim(:)
       endwhere
      endif
!
      return
      end subroutine soilstep

!     ==================
!     SUBROUTINE ROFFINI
!     ==================

      subroutine roffini
      use landmod
!
!     zcvel and zoroinc were parameter(4.2) and a literal 1.0 m2/s2. Both are
!     now derived from namelist keys so that they carry this planet's gravity:
!     see the roffvel block in the module header. glaciermod's oroini does the
!     same, and is the copy that survives, because glacierini runs after
!     landini and recomputes duroff and dvroff.
!
      real :: zcvel
      real :: zcexp
      real :: zoroinc
!
      real zuroff(NLON,NLAT)
      real zvroff(NLON,NLAT)
      real zoro(NLON,NLAT)
      real zoron(NLON,NLAT)
      real zlsm(NLON,NLAT)
      real zsi(NLON,NLAT)
      real zsir(NLON,NLPP)

!
      zcexp   = roffexp
      zcvel   = roffvel*ga**(0.5-roffexp)
      zoroinc = roffpit*ga
!
      ilat = NLAT ! using ilat suppresses compiler warnings for T1
!
      if(NRESTART==0) then

       do jlat=1,NLPP
        do jlon=1,NLON
         jhor=(jlat-1)*NLON+jlon
         darea(jhor)=gwd(jlat)
         zsir(jlon,jlat)=sid(jlat)
        enddo
       enddo

       call mpgagp(zsi,zsir,1)
       call mpgagp(zoro,doro,1)
       call mpgagp(zlsm,dls,1)

       if(mypid==NROOT) then

        zoro(:,:)=MAX(zoro(:,:),0.)
        where(zlsm(:,:) < 1.) zoro(:,:)=zlsm(:,:)-1.

!
!     iterate to remove local minima
!

 1000   continue
        jconv=0
        do jlat=2,ilat-1
         do jlon=2,NLON-1
          if(zlsm(jlon,jlat) > 0.                                       &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat-1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon,jlat+1)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon+1,jlat)                  &
     &      .and. zoro(jlon,jlat) <= zoro(jlon-1,jlat)) then
           zoron(jlon,jlat)=zoroinc+MIN(zoro(jlon+1,jlat),zoro(jlon-1,jlat)  &
     &                            ,zoro(jlon,jlat+1),zoro(jlon,jlat-1))
           jconv=jconv+1
          else
           zoron(jlon,jlat)=zoro(jlon,jlat)
          endif
         enddo
         if(zlsm(1,jlat) > 0.                                           &
     &     .and. zoro(1,jlat) <= zoro(1,jlat-1)                         &
     &     .and. zoro(1,jlat) <= zoro(1,jlat+1)                         &
     &     .and. zoro(1,jlat) <= zoro(2,jlat)                           &
     &     .and. zoro(1,jlat) <= zoro(NLON,jlat)) then
          zoron(1,jlat)=zoroinc+MIN(zoro(2,jlat),zoro(NLON,jlat)             &
     &                        ,zoro(1,jlat+1),zoro(1,jlat-1))
          jconv=jconv+1
         else
          zoron(1,jlat)=zoro(1,jlat)
         endif
         if(zlsm(NLON,jlat) > 0.                                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat-1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(NLON,jlat+1)                   &
     &     .and. zoro(NLON,jlat) <= zoro(1,jlat)                        &
     &     .and. zoro(NLON,jlat) <= zoro(NLON-1,jlat)) then
          zoron(NLON,jlat)=zoroinc+MIN(zoro(1,jlat),zoro(NLON-1,jlat)        &
     &                           ,zoro(NLON,jlat+1),zoro(NLON,jlat-1))
          jconv=jconv+1
         else
          zoron(NLON,jlat)=zoro(NLON,jlat)
         endif
        enddo
        do jlon=2,NLON-1
         if(zlsm(jlon,1) > 0.                                           &
     &    .and. zoro(jlon,1) <= zoro(jlon,2)                            &
     &    .and. zoro(jlon,1) <= zoro(jlon+1,1)                          &
     &    .and. zoro(jlon,1) <= zoro(jlon-1,1)) then
          zoron(jlon,1)=zoroinc+MIN(zoro(jlon+1,1),zoro(jlon-1,1)            &
     &                        ,zoro(jlon,2))
          jconv=jconv+1
         else
          zoron(jlon,1)=zoro(jlon,1)
         endif
         if(zlsm(jlon,NLAT) > 0.                                        &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon,NLAT-1)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon+1,NLAT)                    &
     &    .and. zoro(jlon,NLAT) <= zoro(jlon-1,NLAT)) then
          zoron(jlon,NLAT)=zoroinc+MIN(zoro(jlon+1,NLAT),zoro(jlon-1,NLAT)   &
     &                           ,zoro(jlon,NLAT-1))
          jconv=jconv+1
         else
          zoron(jlon,NLAT)=zoro(jlon,NLAT)
         endif
        enddo
        if(zlsm(1,1) > 0.                                               &
     &    .and. zoro(1,1) <= zoro(1,2)                                  &
     &    .and. zoro(1,1) <= zoro(2,1)                                  &
     &    .and. zoro(1,1) <= zoro(NLON,1)) then
         zoron(1,1)=zoroinc+MIN(zoro(2,1),zoro(NLON,1)                       &
     &                    ,zoro(1,2))
         jconv=jconv+1
        else
         zoron(1,1)=zoro(1,1)
        endif
        if(zlsm(NLON,NLAT) > 0.                                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON,NLAT-1)                    &
     &    .and. zoro(NLON,NLAT) <= zoro(1,NLAT)                         &
     &    .and. zoro(NLON,NLAT) <= zoro(NLON-1,NLAT)) then
         zoron(NLON,NLAT)=zoroinc+MIN(zoro(1,NLAT),zoro(NLON-1,NLAT)         &
     &                          ,zoro(NLON,NLAT-1))
         jconv=jconv+1
        else
         zoron(NLON,NLAT)=zoro(NLON,NLAT)
        endif
        if(zlsm(NLON,1) > 0.                                            &
     &    .and. zoro(NLON,1) <= zoro(NLON,2)                            &
     &    .and. zoro(NLON,1) <= zoro(1,1)                               &
     &    .and. zoro(NLON,1) <= zoro(NLON-1,1)) then
         zoron(NLON,1)=zoroinc+MIN(zoro(1,1),zoro(NLON-1,1)                  &
     &                       ,zoro(NLON,2))
         jconv=jconv+1
        else
         zoron(NLON,1)=zoro(NLON,1)
        endif
        if(zlsm(1,NLAT) > 0.                                            &
     &    .and. zoro(1,NLAT) <= zoro(1,NLAT-1)                          &
     &    .and. zoro(1,NLAT) <= zoro(2,NLAT)                            &
     &    .and. zoro(1,NLAT) <= zoro(NLON,NLAT)) then
         zoron(1,NLAT)=zoroinc+MIN(zoro(NLON,NLAT),zoro(2,NLAT)              &
     &                       ,zoro(1,NLAT-1))
         jconv=jconv+1
        else
         zoron(1,NLAT)=zoro(1,NLAT)
        endif
        zoro(:,:)=zoron(:,:)
        if(jconv > 0 ) goto 1000

        do jlat=1,NLAT
         do jlon=1,NLON-1
          zdx=TWOPI*cos(ASIN(zsi(jlon,jlat)))*plarad/real(NLON)
          zdh=(zoro(jlon+1,jlat)-zoro(jlon,jlat))/zdx
          zfac=1.
          if(zdh > 0.) zfac=-1.
          zuroff(jlon,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
         enddo
         zdx=TWOPI*cos(ASIN(zsi(NLON,jlat)))*plarad/real(NLON)
         zdh=(zoro(1,jlat)-zoro(NLON,jlat))/zdx
         zfac=1.
         if(zdh > 0.) zfac=-1.
         zuroff(NLON,jlat)=zfac*zcvel/zdx*ABS(zdh)**zcexp
        enddo

        do jlat=1,NLAT-1
         do jlon=1,NLON
          zdy=(ASIN(zsi(jlon,jlat+1))-ASIN(zsi(jlon,jlat)))*plarad
          zdh=(zoro(jlon,jlat+1)-zoro(jlon,jlat))/zdy
          zfac=1.
          if(zdh < 0.) zfac=-1.
          zvroff(jlon,jlat)=zfac*zcvel/zdy*ABS(zdh)**zcexp
         enddo
        enddo
        zvroff(1:NLON,NLAT)=0.
!!FL        
!        open(62,file='roffuv.srv',form='unformatted')
!        write(62) 131,0,0,0,NLON,NLAT,0,0
!        write(62) zuroff
!        write(62) 132,0,0,0,NLON,NLAT,0,0
!        write(62) zvroff
!        write(62) 129,0,0,0,NLON,NLAT,0,0
!        write(62) zoro
!        write(62) 172,0,0,0,NLON,NLAT,0,0
!        write(62) zlsm
!        close(62)
!!
       endif

       call mpscgp(zuroff,duroff,1)
       call mpscgp(zvroff,dvroff,1)

       driver(:)=0.
       drunoff(:)=0.

      endif

      return
      end subroutine roffini

!     ===================
!     SUBROUTINE ROFFSTEP
!     ===================

      subroutine roffstep
      use landmod

!
!     make advection
!

      call mkradv(driver,drunoff,duroff,dvroff,darea,dls)

!
!     new water
!

      where(dls(:) > 0.) driver(:)=driver(:)+drunoff(:)*deltsec

      return
      end subroutine roffstep

!     =================
!     SUBROUTINE MKRADV
!     =================

      subroutine mkradv(zriver,zroff,zuroff,zvroff,zarea,zls)
      use landmod

      real zarea(NLON,NLPP)
      real zuroff(NLON,NLPP)
      real zvroff(NLON,NLPP)
      real zriver(NLON,NLPP)
      real zrivin(NLON,NLPP)
      real zroff(NLON,NLPP)
      real zls(NLON,NLPP)

      real zrip(NLON,NLPP)
      real zrop(NLON,NLPP)

      real zrigl(NLON,NLAT),zriglp(NLON,NLAT)
      real zrogl(NLON,NLAT),zroglp(NLON,NLAT)
!

      zrivin(:,:)=zriver(:,:)*zarea(:,:)*zls(:,:)

!
!     make advection in zonal direction (no partitioning in mpp)
!

      do jlon=1,NLON-1
       where(zuroff(jlon,:) > 0.)
        zroff(jlon,:)=zroff(jlon,:)                                     &
     &               -zuroff(jlon,:)*zrivin(jlon,:)/zarea(jlon,:)
        zroff(jlon+1,:)=zroff(jlon+1,:)                                 &
     &                 +zuroff(jlon,:)*zrivin(jlon,:)/zarea(jlon+1,:)
       elsewhere
        zroff(jlon,:)=zroff(jlon,:)                                     &
     &               -zuroff(jlon,:)*zrivin(jlon+1,:)/zarea(jlon,:)
        zroff(jlon+1,:)=zroff(jlon+1,:)                                 &
     &                 +zuroff(jlon,:)*zrivin(jlon+1,:)/zarea(jlon+1,:)
       endwhere
      enddo
      where(zuroff(NLON,:) > 0.)
       zroff(NLON,:)=zroff(NLON,:)                                      &
     &              -zuroff(NLON,:)*zrivin(NLON,:)/zarea(NLON,:)
       zroff(1,:)=zroff(1,:)                                            &
     &           +zuroff(NLON,:)*zrivin(NLON,:)/zarea(1,:)
      elsewhere
       zroff(NLON,:)=zroff(NLON,:)                                      &
     &              -zuroff(NLON,:)*zrivin(1,:)/zarea(NLON,:)
       zroff(1,:)=zroff(1,:)                                            &
     &           +zuroff(NLON,:)*zrivin(1,:)/zarea(1,:)
      endwhere

!
!     make advection in meridional direction (partitioning in mpp)
!     simple schema (exchange all gps)
!

      call mpgagp(zrigl,zrivin,1)
      if(mypid==NROOT) then
       zriglp(:,1:NLAT-1)=zrigl(:,2:NLAT)
       zriglp(:,NLAT)=0.
      endif
      call mpscgp(zriglp,zrip,1)
      where(zvroff(:,:) < 0.)
       zroff(:,:)=zroff(:,:)+zvroff(:,:)*zrivin(:,:)/zarea(:,:)
       zrop(:,:)=-zvroff(:,:)*zrivin(:,:)
      elsewhere
       zroff(:,:)=zroff(:,:)+zvroff(:,:)*zrip(:,:)/zarea(:,:)
       zrop(:,:)=-zvroff(:,:)*zrip(:,:)
      endwhere
      call mpgagp(zroglp,zrop,1)
      if(mypid==NROOT) then
       zrogl(:,2:NLAT)=zroglp(:,1:NLAT-1)
       zrogl(:,1)=0.
      endif
      call mpscgp(zrogl,zrop,1)

      zroff(:,:)=zroff(:,:)+zrop(:,:)/zarea(:,:)

      return
      end subroutine mkradv

!     =================
!     SUBROUTINE GETTCL
!     =================

      subroutine gettcll
      use landmod
!
!     get surface temperature annual cycle
!
      call momint(nperpetual,nstep+1,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      dtclim(:)=zgw1*dtcl(:,jm1)+zgw2*dtcl(:,jm2)
      dwclim(:)=zgw1*dwcl(:,jm1)+zgw2*dwcl(:,jm2)
      return
      end subroutine gettcll

!     =================
!     SUBROUTINE GETALB
!     =================

      subroutine getalb
      use landmod
      
      real :: aa = 5.2
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(aa)
      real :: yy = 4.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(yy)
      real :: bf = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(bf)
      real :: al = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(al)
!
!     get surface background albedo from  annual cycle
!
      call momint(nperpetual,nstep+1,jm1,jm2,zgw2)
      zgw1 = 1.0 - zgw2
      dalbclim(:)=zgw1*dalbcl(:,jm1)+zgw2*dalbcl(:,jm2)
      dalbclim1(:)=zgw1*dalbcl1(:,jm1)+zgw2*dalbcl1(:,jm2)
      dalbclim2(:)=zgw1*dalbcl2(:,jm1)+zgw2*dalbcl2(:,jm2)
      
!      
!     Modify the surface background albedo according to soil water capacity
!
      if (nwetsoil > 0.5) then
        do jhor = 1,NHOR
          if (dls(jhor)>0.5) then
            if (dts(jhor)>273.15) then !We have wet soil. The wetter, the darker.
              bf = dwater(jhor)/wsmax
              al = 1.0/aa*(max(bf,1.0e-3)**((1.0-yy)/yy) - 1)**(1.0/yy)
              dalbclim(jhor) = alblandmax * exp(-(bf*25)**6)+ &
&                              al * (1 - exp(-(bf*25)**6) - exp(-((1-bf)*30)**9))+ &
&                              doceanalb(1) * exp(-((1-bf)*30)**9)
              dalbclim1(jhor) = dalbclim(jhor)
              dalbclim2(jhor) = alblandmax * exp(-(bf*25)**6)+ &
&                              al * (1 - exp(-(bf*25)**6) - exp(-((1-bf)*30)**9))+ &
&                              doceanalb(2) * exp(-((1-bf)*30)**9)
            else !Frozen soil.
              dalbclim(jhor) = albland
              dalbclim1(jhor) = dgroundalb(1)
              dalbclim2(jhor) = dgroundalb(2)
            endif
          endif
        enddo
      endif
      
      return
      end subroutine getalb

