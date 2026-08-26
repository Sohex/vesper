!	This is a rudimentary glacier model for capturing the effects of large ice sheets
!	on atmospheric circulation and thermodynamics. It's not intended to be used as a 
!	real-time coupling, merely as a coupling between an external gridpoint glacier
!	model that may run on timescales of years to decades to centuries, and PlaSim, running
!	on timescales of minutes to hours. The external model sets the glacier/snowpack height,
!	and this translates that into a change in orography through the surface geopotential
!	height. This does this by keeping track of two orography fields; the lithographic orography
!	and the glacier orography. The surface orography is the sum of the two. Every timestep, the
!	model checks to see if existing snow on a gridpoint has melted away (default reaching a
!       depth below 2.0 meters liquid water equivalent, GLACELIM). A gridpoint that holds snow
!	above that depth continuously for GLACPERSIST orbits becomes a glacier; if its snow depth
!	later falls below GLACELIM the clock restarts, and below 0.1 m the glacier flag is cleared.
!	The elapsed time is carried in the restart, so the criterion does not depend on how the
!	run was cut into segments.
!
!	The orography is recomputed by oroini, which runs at INITIALISATION only, so the ice
!	sheet's effect on the circulation appears at segment boundaries rather than continuously.
!	That is the coupling this module was written for and not a defect.
!
!	Due to the ground orography implementation being introduced, this module also introduces
!	the possibility of coupling to an orogeny model.

      module glaciermod
      use landmod
!
!     version identifier (date)
!
      character(len=80) :: gversion = '07.07.2017 by Adiv'
!
!     Parameter   
!

!
!     namelist parameters
!
      character (256) :: glacier_namelist = "glacier_namelist"
      
      integer :: nglacier = 1 ! 1 = implement glaciation, 0 = ignore
      real :: glacelim = 2.0 ! Minimum snow depth in meters liquid water equivalent that has to be maintained
                                 ! year-round to convert the gridpoint to a glacier.
      real :: icesheeth = -1.0 !Initial snow depth

!     How long the snow has to stay above glacelim, in ORBITS.
!
!     The module header says "continuously for an entire year". What the code
!     did was continuously for one INVOCATION of the model: persistflag was set
!     .TRUE. in the declaration, glacierstep only ever cleared it, and
!     glacierstop converted whatever survived. So the duration was the segment
!     length, it was not reachable from any namelist, and it silently changed
!     meaning with the caller's choice of segment. persistt below counts the
!     time instead, and it is carried in the restart, so the criterion is the
!     same however a run is cut into segments.
!
!     1.0 is one orbit, which is what the header asks for and what a
!     one-orbit-per-invocation caller was already getting.
      real :: glacpersist = 1.0
      real :: glacpersec  = 0.0 ! glacpersist in seconds; set in glacierprep
!
!     global arrays
!
      real :: groundoro(NHOR)  = 0.0
      real :: glacieroro(NHOR) = 0.0
      real :: netoro(NHOR) = 0.0
      real :: persistt(NHOR) = 0.0 ! seconds of continuous cover above glacelim
      
!
!     global scalars
!
!     rhoglac converts water equivalent to ice thickness for the orography, as
!     rhosnow in landmod does for the snow that insulates the soil column. Both
!     are settled densities set by overburden compaction, so both are measured
!     under the gravity of the planet they were measured on and both are low
!     for a planet whose surface gravity is higher. INHERITED rather than
!     rescaled: firn densification is not a closed-form function of gravity
!     that can be evaluated here, and the effect is bracketed at roughly 10 to
!     20 per cent on snow insulation thickness and on glacier relief rather
!     than guessed at. Both are namelist keys so the bracket can be run.
      real :: rhoglac  = 850.    ! glacial ice density (kg/m**3)

!     WHAT THE MODELLED GLACIAL ICE IS MADE OF, DERIVED FROM THAT DENSITY.
!     WORLD-FG8W.
!
!     `landmod`'s `sicecap` and `sicediff` are the heat capacity per unit volume
!     and the thermal conductivity of the ice `glaciermod` grows, and they weight
!     the soil column's thermal properties by `dglac`. Both followed nothing: the
!     capacity factorised exactly as one thousand times ice's specific heat --
!     LIQUID WATER's density, not the ice's -- so a bracket that moved `rhoglac`
!     moved the ice orography and left the thermal mass of the ice behind, which
!     is the defect GRAV-8 removed from `snowcap` and WORLD-A9S5 from `snowdiff`.
!
!     THEY ARE DERIVED HERE AND NOT IN `landini`, because this is where the
!     density they follow is declared. `glacierprep` runs before `landini` in
!     `surfini` and this module already uses `landmod`, so each thread writes its
!     own copy of both before anything reads them; `landini` no longer touches
!     either and they have left `landmod_nl`. Deriving them in `landini` instead
!     would need `landmod` to use `glaciermod`, which is a cycle: this module
!     uses that one.
!
!     THE CAPACITY is `rhoglac * CPGLAC`. `CPGLAC` is the specific heat of ice Ih
!     from IAPWS-06 at `TGLACREF`, computed by `analysis/ice_properties.py` from
!     the release's own Gibbs function and checked against every entry of its
!     numerical check table before it is reported. That script holds this literal
!     to the value it computes, so the two cannot drift.
!
!     THE CONDUCTIVITY is pure ice's, reduced for the air the density implies.
!     Two steps, both Yen (1981), CRREL Report 81-10:
!
!     - Pure ice, his Eq. (33), `lambda = a exp(b T)`, the whole-range arm of his
!       Table 3. He recommends it for practical use because it has the highest
!       correlation coefficient of his three arms, 0.9313 against 0.5962 for the
!       arm fitted above the 150-195 K data gap. IAPWS-06 CANNOT supply this: a
!       Gibbs function carries density, specific heat and compressibility and no
!       transport property at all.
!     - Bubbles, his Eq. (37), `2 rho / (3 rhoice - rho)`, which is Schwerdtfeger's
!       reduction of Maxwell's effective-medium result for randomly distributed
!       spherical air inclusions once the conductivity of air is dropped against
!       the ice's. Yen states the same relation twice, as Eq. (36)/(37) for dense
!       snow and as Eq. (70) for the bubbly ice inside sea ice, and the two agree
!       to two parts in a thousand at this density. Glacial ice IS bubbly, and
!       his Figure 22 is the same equation drawn.
!
!     `RHOICE_YEN1981` is the pure-ice density Yen's air-fraction relations are
!     written against and is part of THEM. It is not `icemod`'s `CRHOI`, which is
!     the density of the modelled SEA ice, and it is not `landmod`'s
!     `RHOICE_F2021`, which is the normalising density of a snow conductivity fit;
!     the three must not be deduplicated into one another because each belongs to
!     a different relation.
!
!     NOT `icemod`'s `CKAPI`, WHICH IS THE SAME NUMBER TODAY BY COINCIDENCE. That
!     is the conductivity of the modelled SEA ice, which is brine-bearing, and
!     Yen's Eq. (71) subtracts a brine term from exactly the bubbly ice computed
!     here -- so the two sit on opposite sides of pure ice's value for different
!     reasons and deduplicating them would assert that glacial ice is salty.
!
!     THE DECLARED TEMPERATURE IS THE LARGER UNCERTAINTY, not the density. Over
!     233.15 K to the melting point the conductivity moves by about a quarter and
!     the specific heat by about a sixth, against a fifteenth for the difference
!     between Yen's two conductivity arms. `TGLACREF` is the temperature
!     `landmod` already evaluates the snow conductivity at, so the two cryosphere
!     materials are stated at one temperature rather than two.
!     `notes/audits/cryosphere-material-properties.md` carries the sweep.
      real, parameter :: TGLACREF       = 263.15  ! declared ice temperature (K)
      real, parameter :: CPGLAC         = 2023.10 ! ice Ih specific heat (J/kg/K)
      real, parameter :: RHOICE_YEN1981 = 917.    ! Yen's pure ice density (kg/m3)
      real, parameter :: YENICE_A       = 9.828   ! Eq. (33) prefactor (W/m/K)
      real, parameter :: YENICE_B       = -0.0057 ! Eq. (33) exponent (1/K)

!

!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(glacelim,glacier_namelist,glacieroro,glacpersec,glacpersist,groundoro,&
!$omp&  gversion,icesheeth,netoro,nglacier,persistt,rhoglac)

      end module glaciermod
          
!      
!     =============================
!

      subroutine glacierprep
      use glaciermod

      namelist/glacier_nl/nglacier,glacelim,icesheeth,glacpersist,rhoglac
      
      if (mypid==NROOT) then
         open(23,file=glacier_namelist)
         read(23,glacier_nl)
         close(23)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * GLACIERMOD ",a34," *")') trim(gversion)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist GLACIER_NL from <glacier_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,glacier_nl)
      endif
      
      call mpbci(nglacier)
      call mpbcr(glacelim)
      call mpbcr(icesheeth)
      call mpbcr(glacpersist)
      call mpbcr(rhoglac)

!     The modelled glacial ice's thermal properties, from the density every
!     thread has just been given, so the two cannot drift apart. WORLD-FG8W; the
!     argument is above the declarations. `landini` runs after this and no longer
!     writes either.
      sicecap  = rhoglac * CPGLAC
      sicediff = YENICE_A * exp(YENICE_B * TGLACREF)                      &
     &         * 2.0 * rhoglac / (3.0 * RHOICE_YEN1981 - rhoglac)

      if (mypid==NROOT) then
         write(nud,'(" * Glacial ice at ",f7.2," K: cap ",e12.4,          &
     &        " J/m3/K, diff ",f7.4," W/m/K *")') TGLACREF,sicecap,sicediff
      endif

!     Orbits to seconds. m_days_per_year * day_24hr is the orbital period, the
!     same conversion landini uses for tau_veg and tau_soil, and initpm has
!     already set both from n_days_per_year and the rotation rate.

      glacpersec = glacpersist * m_days_per_year * day_24hr

      if (mypid==NROOT) then
         write(nud,'(" * Glacier persistence: ",f8.3," orbits =",e12.4," s *")') &
     &        glacpersist,glacpersec
      endif

!     ROOT PRINTS, EVERY THREAD STOPS -- glacpersist is broadcast, so the
!     branch is taken on every thread and nud is one shared unit. world-0ihs.
      if (glacpersec <= 0.0) then
         if (mypid == NROOT) then
            write(nud,*)' *** error: glacpersist = ',glacpersist
         endif
         stop 1
      endif
      
      if (mypid==NROOT) nglspec = nglacier
      
      call mpbci(nglspec)
      
      end subroutine glacierprep
      
!      
!     =============================
!

      subroutine glacierini(noromax)
      use glaciermod
      
      logical :: ldsnow
      real :: zoro(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zoro)
      real :: foro(NHOR) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(foro)
      real :: zsnow(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(zsnow)
      
      integer noromax
      
      if (nrestart > 0.) then
        call mpgetgp('dglac',dglac,NHOR,1)
        call mpgetgp('groundsg',groundoro ,NHOR,1)
        call mpgetgp('dglacsg' ,glacieroro,NHOR,1)
        call mpgetgp('doro'    ,doro      ,NHOR,1)
        call mpgetgp('persistt',persistt  ,NHOR,1)
      else
        call mpsurfgp('doro'    ,doro    ,NHOR,1)
!       THE ONE PLACE OROSCALE IS APPLIED, and it is on the cold path because
!       this is where the staged field first enters the model. ExoPlaSim
!       invokes the executable once per orbit, so anything scaled outside this
!       branch is scaled again on every invocation: groundoro comes back out of
!       the restart record groundsg above and doro out of 'doro', both already
!       carrying whatever scaling the cold start applied. oroini used to scale
!       groundoro and glacierini used to scale doro a second time on top of it,
!       so a non-unit oroscale raised the staged orography to the power of the
!       number of invocations. Inert at the compiled default of 1.0, which is
!       why it stood; multiplying by exactly 1.0 is the identity, so moving it
!       here leaves groundsg bit-equal to the staged field. world-6qee.
        doro(:) = doro(:) * oroscale
        groundoro(:) = doro(:)
        glacieroro(:) = 0.0
        if (icesheeth .ge. 0.0) then
          dglac(:) = 1.0
          dsnowz(:) = icesheeth
          do jhor=1,NHOR
            if (groundoro(jhor) == 0.0) then
              dsnowz(jhor) = 0.0
              dglac(jhor) = 0.0
            endif
          enddo
        else
          do jhor=1,NHOR
            if ((dglac(jhor) .gt. 0.0).and.(dsnowz(jhor) .eq. 0.0)) dsnowz(jhor)=glacelim
          enddo
        endif
        
          
        call mpgagp(zsnow,dsnowz,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"Initial Icesheet Heights")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zsnow))
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zsnow))
          write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zsnow))
        endif 
      endif
      
      if (nglacier .eq. 1) then
      
        if (nrestart > 0.) then     
          call mpgagp(zoro,doro,1)
        
          if (mypid==NROOT) then
            write(nud,'(/,"Topography before glaciers")')
            write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
            write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
            write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
          endif        
        endif
        
        if (mypid==NROOT) inquire(file='restart_dsnow',exist=ldsnow)
        call mpbcl(ldsnow)
        
        if (ldsnow .and. (nrestart > 0.)) then
        
           call readarray(dsnowz,'restart_dsnow')
           
        endif   
        
        call oroini       
           
!          Compute spectral orography
        so(:) = 0.
!       doro is groundoro + glacieroro as oroini just built it, and groundoro
!       was scaled where the staged field entered the model. Scaling it again
!       here applied oroscale twice to the ground part within one call and once
!       more on every later invocation. world-6qee.
        foro(:) = doro(:)
        
        
        call gp2fc(foro,NLON,NLPP)
        call fc2sp(foro,so)
        call mpsum(so,1)
        
        if (npro == 1) then ! print only in single core runs
           call sp2fc(so,doro)
           call fc2gp(doro,nlon,nlpp)
           write(nud,'(/,"Topography after spectral fitting")')
           write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
           write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
           write(nud,'("Mean:    ",f10.2," [m]")') ugpmean(doro) / ga
        endif
        
        if (mypid == NROOT) then
           so(:) = so(:) / (cv*cv)
           if (noromax < NTRU) then
            jr=-1
            do jm=0,NTRU
             do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if(jn > noromax) then
               so(jr)=0.
               so(ji)=0.
              endif
             enddo
            enddo
           endif ! (noromax < NTRU)
        
!          Initialize surface pressure
        
          if (nspinit > 0) then
             sp(:) = -so(:)*cv*cv / (gascon * tgr)
          endif
        endif ! (mypid == NROOT)
!       ON THE COLD PATH ONLY. This scatters the surface pressure into the
!       leapfrog MINUS level, which is what a run starting from rest wants:
!       the two levels are the same state at t = 0. On the restart path
!       read_atmos_restart has already restored spm as the t - dt level, and
!       glacierini runs AFTER it -- surfini is called from prolog below the
!       restart read -- so an unguarded scatter overwrote the minus level with
!       the current one. That made the first step after every segment boundary
!       a different step from the one the uninterrupted run took, and it was
!       the whole of the restart discontinuity: with this guard a run split
!       into two segments reproduces the same run taken whole in every restart
!       record, and without it the two differ in 46 of them by the first step.
!       world-8yyh; notes/audits/ecological-stream-restart-continuity.md.
        if (nrestart == 0) call mpscsp(sp,spm,1)
        
!       A CELL-MEAN DEPTH AGAINST A COLUMN THRESHOLD. 30 m is the minimum ice
!       thickness for a sheet to flow, which is a property of ice and not of the
!       grid, but dsnowz here is the mean over the cell: a coarse cell that is
!       half covered to 40 m does not flag while a fine cell fully covered to
!       31 m does. Anchored to T21 in that sense. world-khn.
        where (dsnowz(:) > 30.0) dglac = 1.0 !If we have more than 30 m of lq H2O equivalent in snow/ice, it's a glacier
                                             !30 meters of ice is the minimum thickness for an ice sheet to flow
        where (dls(:) < 0.5) persistt = 0.0
        
        call mpgagp(zoro,doro,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"New Topography after glaciers")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
          write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
        endif
      
      else    !nglacier == 0
       
        if (nrestart > 0.) then     
          call mpgagp(zoro,doro,1)
        
          if (mypid==NROOT) then
            write(nud,'(/,"Topography before smoothing")')
            write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
            write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
            write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
          endif    
        endif
          
        call oroini 
                
!        Compute spectral orography
        so(:) = 0.
!       doro is groundoro + glacieroro as oroini just built it, and groundoro
!       was scaled where the staged field entered the model. Scaling it again
!       here applied oroscale twice to the ground part within one call and once
!       more on every later invocation. world-6qee.
        foro(:) = doro(:)
        
        
        call gp2fc(foro,NLON,NLPP)
        call fc2sp(foro,so)
        call mpsum(so,1)
        
        if (npro == 1) then ! print only in single core runs
           call sp2fc(so,doro)
           call fc2gp(doro,nlon,nlpp)
           write(nud,'(/,"Topography after spectral fitting")')
           write(nud,'("Maximum: ",f10.2," [m]")') maxval(doro) / ga
           write(nud,'("Minimum: ",f10.2," [m]")') minval(doro) / ga
           write(nud,'("Mean:    ",f10.2," [m]")') ugpmean(doro) / ga
        endif
        
        if (mypid == NROOT) then
           so(:) = so(:) / (cv*cv)
           if (noromax < NTRU) then
            jr=-1
            do jm=0,NTRU
             do jn=jm,NTRU
              jr=jr+2
              ji=jr+1
              if(jn > noromax) then
               so(jr)=0.
               so(ji)=0.
              endif
             enddo
            enddo
           endif ! (noromax < NTRU)
        
!          Initialize surface pressure
        
          if (nspinit > 0) then
             sp(:) = -so(:)*cv*cv / (gascon * tgr)
          endif
        endif ! (mypid == NROOT)
!       ON THE COLD PATH ONLY. This scatters the surface pressure into the
!       leapfrog MINUS level, which is what a run starting from rest wants:
!       the two levels are the same state at t = 0. On the restart path
!       read_atmos_restart has already restored spm as the t - dt level, and
!       glacierini runs AFTER it -- surfini is called from prolog below the
!       restart read -- so an unguarded scatter overwrote the minus level with
!       the current one. That made the first step after every segment boundary
!       a different step from the one the uninterrupted run took, and it was
!       the whole of the restart discontinuity: with this guard a run split
!       into two segments reproduces the same run taken whole in every restart
!       record, and without it the two differ in 46 of them by the first step.
!       world-8yyh; notes/audits/ecological-stream-restart-continuity.md.
        if (nrestart == 0) call mpscsp(sp,spm,1)
        
        call mpgagp(zoro,doro,1)
        
        if (mypid==NROOT) then
          write(nud,'(/,"New Topography after glaciers")')
          write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
          write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
          write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
        endif
      
      endif !nglacier switch
      
      end subroutine glacierini

     
!     ==================
!     SUBROUTINE oroini
!     ==================

      subroutine oroini
      use glaciermod
!
!     zcvel and zoroinc were parameter(4.2) and a literal 1.0 m2/s2. Both now
!     come from landmod's roffvel, roffexp and roffpit so that they carry this
!     planet's gravity; the derivation and the Earth identity are in the
!     roffvel block in landmod's module header. This is the copy that survives:
!     glacierini runs after landini, so oroini overwrites whatever roffini put
!     in duroff and dvroff.
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
      
      real thing
      real dhsnow
      real radius
      
!
      zcexp   = roffexp
      zcvel   = roffvel*ga**(0.5-roffexp)
      zoroinc = roffpit*ga
!
      ilat = NLAT ! using ilat suppresses compiler warnings for T1
!
      do jlat=1,NLPP
       do jlon=1,NLON
        jhor=(jlat-1)*NLON+jlon
        darea(jhor)=gwd(jlat)
        zsir(jlon,jlat)=sid(jlat)
       enddo
      enddo 
      
      call mpgagp(zoro,doro,1)

      if(mypid==NROOT) then

       write(nud,'(/,"Topography before glaciers and before smoothing")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      
      if (nglacier .gt. 0.5) then
      
!     groundoro arrives already scaled -- from the cold start's own read, or
!     from the restart record that start wrote. It is NOT scaled here: oroini
!     runs on every model invocation and this compounded. world-6qee.
      
!     Glacier thickness enters the orography as geopotential, and every reader
!     of doro divides by ga to get metres back, so the gravity used here has to
!     be the model's. zgm is the GM the model's own surface gravity implies on
!     its own radius, which keeps the height dependence below while making
!     grav = ga at the surface.

      radius = plarad
      zgm    = ga*radius**2

      do i=1,NHOR ! Add elevation of snowpack/ice sheet
         dz = radius**2*groundoro(i)/(zgm-radius*groundoro(i))
         dhsnow = dsnowz(i)/(rhoglac/1000.0) !Convert from meters of lq H20 equivalent to actual snow depth
         grav = zgm/(radius+dz)**2
         glacieroro(i) = grav*dhsnow - grav/(dz+radius)*dhsnow**2
         doro(i) = groundoro(i) + glacieroro(i)
!          if(mypid==NROOT) then
!          write(nud,'("Gr,Sn: ",f10.2," ",f10.2," ",f10.2," [m]")') (groundoro(i) / ga,glacieroro(i)/ga)
!          write(nud,'("Net: ",f10.2," [m]")') (doro(i)/ga)
!          endif 

      enddo
      
      call mpgagp(zoro,glacieroro,1)
      
      if(mypid==NROOT) then

       write(nud,'(/,"Glacial Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      
      call mpgagp(zoro,groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Ground Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      call mpgagp(zoro,doro-groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Net-Ground Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      call mpgagp(zoro,doro-glacieroro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Net-Glacier Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      call mpgagp(zoro,glacieroro+groundoro,1)
      if(mypid==NROOT) then

       write(nud,'(/,"Ground + Glacier Topography")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
      endif
      
      endif
      
      call mpgagp(zsi,zsir,1)
      call mpgagp(zoro,doro,1)
      call mpgagp(zlsm,dls,1)

      if(mypid==NROOT) then

       write(nud,'(/,"Topography after glaciers and before smoothing")')
       write(nud,'("Maximum: ",f10.2," [m]")') (maxval(zoro) / ga)
       write(nud,'("Minimum: ",f10.2," [m]")') (minval(zoro) / ga)
       write(nud,'("Mean:    ",f10.2," [m]")') (ugpmean(zoro) / ga)
      
       zoro(:,:)=MAX(zoro(:,:),0.)
       where(zlsm(:,:) < 1.) zoro(:,:)=zlsm(:,:)-1.

!
!     iterate to remove local minima, but only if we're not already flat.
!
       if (ndesert < 0.5 .and. naqua < 0.5 .and. ((maxval(zoro)-minval(zoro))/ga>0.01)) then

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
       
       endif

      endif

      call mpscgp(zuroff,duroff,1)
      call mpscgp(zvroff,dvroff,1)

      driver(:)=0.
      drunoff(:)=0.

      return
      end subroutine oroini           
      
!      
!     =============================
!
      subroutine glacierstep
      use glaciermod
      
      if (nglacier .eq. 1) then
      
      do jhor = 1,NHOR
         if (dls(jhor) > 0.5) then
            if (dsnowz(jhor) < glacelim) then
               persistt(jhor) = 0.0                 !Snowpack below the persistence threshold: the clock restarts
            else
               persistt(jhor) = persistt(jhor) + deltsec
               if (persistt(jhor) >= glacpersec) dglac(jhor) = 1.0
            endif
            if (dsnowz(jhor) < 0.1) dglac(jhor) = 0.0 !If the snow/ice is basically gone, so is the glacier
         else
            persistt(jhor) = 0.0
         endif
      enddo
      
      endif
      
      return
      end subroutine glacierstep
      
!     
!     ============================
!

      subroutine glacierstop
      use glaciermod

      real :: rpersist(NHOR) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(rpersist)
      
      if (nglacier .eq. 1) then
      
      call finishup(dsnowz,'newdsnow')
      call finishup(asndch,'restart_snow')

!     The conversion is glacierstep's now: it fires the moment the accumulated
!     cover reaches glacpersec, so it no longer depends on where the caller cut
!     the run into segments. What is left here is the diagnostic field, which
!     is set rather than only raised so that a cell whose clock has restarted
!     reports zero.

      rpersist(:) = 0.0
      where (persistt(:) >= glacpersec) rpersist(:) = 1.0
      
      endif
      
      call finishup(dls,'lsm')
      call finishup(rpersist,'persist')
      
      call mpputgp('groundsg',groundoro ,NHOR,1)
      call mpputgp('dglacsg' ,glacieroro,NHOR,1)
      call mpputgp('doro'    ,doro      ,NHOR,1)
      call mpputgp('persistt',persistt  ,NHOR,1)
      
!       endif
      
      call mpputgp('dglac'   ,dglac     ,NHOR,1)

      return
      end subroutine glacierstop
      
!
!     ===============================


