!
!     ***********************
!     * Planet Simulator 17 *
!     ***********************

!     ********************
!     * Frank Lunkeit    *       University of Hamburg
!     * Edilbert Kirk    *      Meteorological Institute
!     * KLaus Fraedrich  *
!     * Valerio Lucarini *
!     ********************

!     **************
!     * Name rules *
!     **************

!     i - local integer
!     j - loop index
!     k - integer dummy parameter
!     l - logical
!     n - global integer

!     g - real gridpoint arrays
!     p - real dummy parameter
!     s - real spectral arrays
!     y - character
!     z - local real

      program plasim_main
      use pumamod

!     based on PUMA (Portable University Model of the Atmosphere)
!     which is based on SGCM (Simple Global Circulation Model)
!     by Ian James and James Dodd (April 1993)
!     Dept of Meteorology, University of Reading, UK.

!     There is a linear drag which can again vary with level,
!     the time scale is entered for each level in days in tfrc(NLEV).
!     the diffusion is del-ndel with a time scale for diffusion of
!     tdiss days on every variable at the truncation wavenumber.
!
!     UPDATE VERSION IDENTIFIER AFTER EACH CODE CHANGE!

plasimversion = "https://github.com/Edilbert/PLASIM/ : 15-Dec-2015"

!     THE WHOLE RUN IS ONE PARALLEL REGION, and it has to be: the model's
!     state is module data declared !$omp threadprivate, so a thread's copy
!     lives exactly as long as the region does. Opening a region per phase
!     would hand each phase a fresh, empty model.
!
!     Inert without -fopenmp, so the MPI build and the serial build read
!     this as five comment lines and are unchanged. Under -fopenmp the team
!     is fixed at NPRO -- the binary is compiled for a thread count the way
!     the MPI binary is compiled for a rank count -- and mpstart checks it.
!
!     No copyin clause. gfortran's TLS initialisation image carries a
!     module variable's declaration initialiser to every thread, which
!     mpstart verifies on entry rather than trusting.
!$omp parallel num_threads(NPRO) default(shared)
      call mpstart(-1)       ! -1: Start MPI   >=0 arg = MPI_COMM_WORLD
      call setfilenames
      call opendiag
      if (mrnum == 2) then
         call mrdimensions
      endif
      call allocate_arrays
      call prolog
      call master
      call epilog
      call mpstop
!$omp end parallel

      stop
      end

      ! ***************************
      ! * SUBROUTINE SETFILENAMES *
      ! ***************************
      
      subroutine setfilenames
      use pumamod
      
      character (3) :: mrext
      
      if (mrpid <  0) return ! no multirun
      
      write(mrext,'("_",i2.2)') mrpid
      
      plasim_namelist     = trim(plasim_namelist    ) // mrext
      radmod_namelist     = trim(radmod_namelist    ) // mrext
      miscmod_namelist    = trim(miscmod_namelist   ) // mrext
      fluxmod_namelist    = trim(fluxmod_namelist   ) // mrext
      rainmod_namelist    = trim(rainmod_namelist   ) // mrext
      surfmod_namelist    = trim(surfmod_namelist   ) // mrext
      plasim_output       = trim(plasim_output      ) // mrext
      plasim_diag         = trim(plasim_diag        ) // mrext
      plasim_restart      = trim(plasim_restart     ) // mrext
      plasim_status       = trim(plasim_status      ) // mrext
      planet_namelist     = trim(planet_namelist    ) // mrext
      efficiency_dat      = trim(efficiency_dat     ) // mrext
      icemod_namelist     = trim(icemod_namelist    ) // mrext
      ice_output          = trim(ice_output         ) // mrext
      oceanmod_namelist   = trim(oceanmod_namelist  ) // mrext
      ocean_output        = trim(ocean_output       ) // mrext
      landmod_namelist    = trim(landmod_namelist   ) // mrext
      vegmod_namelist     = trim(vegmod_namelist    ) // mrext
      seamod_namelist     = trim(seamod_namelist    ) // mrext
      aero_namelist       = trim(aero_namelist      ) // mrext
      
      return
      end

!     ***********************
!     * SUBROUTINE OPENDIAG *
!     ***********************

      subroutine opendiag
      use pumamod
      
      if (mypid == NROOT) then
         open(nud,file=plasim_diag)
      endif
      
      return
      end


!     ******************************
!     * SUBROUTINE ALLOCATE_ARRAYS *
!     ******************************

      subroutine allocate_arrays
      use pumamod
      
      if (mrnum == 2) then
         allocate(sdd(nesp,nlev))   ; sdd(:,:)  = 0.0
         allocate(std(nesp,nlev))   ; std(:,:)  = 0.0
         allocate(szd(nesp,nlev))   ; szd(:,:)  = 0.0
         allocate(spd(nesp     ))   ; spd(:  )  = 0.0
      endif
      
      return
      end subroutine allocate_arrays


!     =================
!     SUBROUTINE PROLOG
!     =================

      subroutine prolog
      use pumamod
      use shtnsmod, only: shtns_setup

      logical :: lrestart

      real (kind=8) :: zsid(NLAT)     ! sid, gwd, csq and rcs reordered for
      real (kind=8) :: zgwd(NLAT)     ! the scatter that follows
      real :: zcsq(NLAT)
      real :: zrcs(NLAT)

!     ************************************************************
!     * Initializations that cannot be run on parallel processes *
!     ************************************************************

      if (mypid == NROOT) then
         call cpu_time(tmstart)
         write(nud,'(54("*"))')
         write(nud,'("* ",17X,"PLANET SIMULATOR",17X," *")')
         write(nud,'("* ",a50," *")') plasimversion
         write(nud,'(54("*"))')
         if (mrnum > 1) then
            write(nud,'("* Instance ",i3," of ",i3,32x,"*")') &
                  mrpid+1, mrnum
            write(nud,'("* My    truncation  :",i5,27x,"*")') NTRU
            write(nud,'("* Other truncation  :",i5,27x,"*")') mrtru(2-mrpid)
            write(nud,'(54("*"))')
         endif
         write(nud,'("* Truncation   NTRU :",i5,27x,"*")') NTRU
         write(nud,'("* Levels       NLEV :",i5,27x,"*")') NLEV
         write(nud,'("* Latitudes    NLAT :",i5,27x,"*")') NLAT
         write(nud,'("* Longitudes   NLON :",i5,27x,"*")') NLON
         write(nud,'(54("*"))')
         if (NPRO > 1) then
            write(nud,'(54("*"))')
            do jpro = 1 , NPRO
              write(nud,'("* CPU",i4,1x,a42," *")') jpro-1,ympname(jpro)
            enddo
            write(nud,'(54("*"))')
         endif
      endif ! (mypid == NROOT)

      call planet_ini               ! Define planet
      call mpbcr(day_24hr)
      
      if (mypid == NROOT) then
         call restart_ini(lrestart,plasim_restart)
         if (lrestart) then
            nrestart = 1
            nkits    = 0
            ndivdamp = 0
         endif
         call surface_ini           ! Read boundary and other data
         call inigau(NLAT,sid,gwd)  ! Gaussian abscissas and weights
         call inilat                ! Set latitudinal arrays
         call readnl                ! Open and read <plasim_namelist>
!        AFTER readnl, NOT BEFORE IT. world-1o4. `print_planet` reports
!        sidereal_day and the orbit period derived from it, and readnl is where
!        sidereal_day is recomputed from rotspd. Called first, the table
!        reported a 23.93 h rotation and a 183.30-day orbit for a model that
!        then integrated 30 h and 146 days. It writes and reads nothing else, so
!        the only thing this ordering changes is that the numbers are true.
         call print_planet
         call initpm                ! Several initializations
         call initsi                ! Initialize semi implicit scheme
         call guistart              ! Initialize GUI
         if (nsela > 0) call tracer_ini0 ! initialize tracer data
         if (nsela > 0 .and. l_aero > 0) call aero_ini 
      endif ! (mypid == NROOT)

!     ***********************
!     * broadcast & scatter *
!     ***********************

!     These four are the whole of what a process knows about WHERE its
!     latitudes are, and everything else derives from them: deglat, cola and
!     rcsq just below, the Legendre weight matrices in legini, the zenith
!     angle in radmod. So reordering them here is what carries the paired
!     decomposition into every latitude-dependent quantity in the model.
!
!     inilat has already run and filled the global arrays, and nothing reads
!     them globally after this point -- tracer_ini0 is the last to do so and
!     it is above, inside the root-only block.


      call mpscdn(sid ,NLPP)  ! sine of latitude (kind=8)
      call mpscdn(gwd ,NLPP)  ! gaussian weights (kind=8)
      call mpscrn(csq ,NLPP)  ! cosine squared of latitude
      call mpscrn(rcs ,NLPP)  ! 1.0 / cos(lat)

      do jlat = 1 , NLPP
         deglat(jlat) = 180.0 / PI * asin(sid(jlat))
         cola(jlat) = sqrt(csq(jlat))
         rcsq(1+(jlat-1)*NLON:jlat*NLON) = 1.0 / csq(jlat)
      enddo

      call mpbci(nfixorb ) ! Global switch to fix orbit
      call mpbci(ngenkeplerian) ! Global switch for general Keplerian orbits
      call mpbci(ntpal   ) ! color pallet for T

      call mpbci(kick    ) ! add noise for kick > 0
      call mpbci(naqua   ) ! aqua planet switch
      call mpbci(ndesert ) ! desert planet switch
      call mpbci(nveg    ) ! vegetation switch
      call mpbci(noutput ) ! write data switch
      call mpbci(nafter  ) ! write data interval
!
!     nlowio and nstpw are read from plasim_nl on NROOT only, exactly like every
!     key broadcast around here, but unlike them they were never broadcast at
!     all. Every non-root task therefore kept the COMPILED DEFAULT nlowio = 1
!     (plasimmod.f90) while NROOT held whatever the namelist said. That is
!     invisible while the namelist also says 1, and it DEADLOCKS the model when
!     it says 0: the `if (nlowio .eq. 0)` blocks in outmod.f90 call writegp and
!     writesp, and writegp opens with the COLLECTIVE mpgagp. NROOT enters that
!     gather and no other task does, so the tasks desynchronise and the run sits
!     in mismatched collectives -- 100% CPU, no system time, no output past the
!     40-byte header, forever.
!
!     nafter is broadcast on the line above and is derived from both of these on
!     NROOT, so the output CADENCE was always consistent across tasks. Only the
!     branch that decides WHICH fields to write was not, which is why this hid
!     for as long as the namelist agreed with the default.
!
      call mpbci(nlowio  ) ! low-I/O accumulation mode (0/1)
      call mpbci(nstpw   ) ! timesteps between writes
      call mpbci(nwpd    ) ! number of writes per day
      call mpbci(nsnapshot) ! Switch for writing snapshots
      call mpbci(nstps   ) ! number of steps per snapshot
      call mpbci(nhcadence) ! Switch for a burst of high-cadence output
      call mpbci(hcstartstep) ! Timestep to start high-cadence output
      call mpbci(hcendstep)   ! Timestep to end high-cadence output (exclusive)
      call mpbci(hcinterval)  ! Number of timesteps per high-cadence write
      call mpbci(neco    ) ! Switch for the ecological output stream, EFOR-2
      call mpbci(necostep) ! Timesteps per ecological interval
      call mpbci(ncoeff  ) ! number of modes to print
      call mpbci(ndiag   ) ! write diagnostics interval
      call mpbci(ndivdamp) ! divergence damping countdown
      call mpbci(ngui    ) ! GUI on (1) or off (0)
      call mpbci(sellon )  ! index of longitude for column mode
      call mpbci(nkits   ) ! number of initial timesteps
      call mpbci(nrestart) ! 1: read restart file 0: initial run
      call mpbci(nqspec  ) ! 1: spectral q 0: grodpoint q
      call mpbci(nshtns  ) ! 1: SHTns transforms 0: legmod's own
      call mpbci(nsela   ) ! 1: semi lagrangian advection enabled
      call mpbci(l_aero  ) ! 1: aerosols enabled

!     The dust emission scheme's three boundary fields, codes 1801-1803.
!     This is called by EVERY rank and aero_ini is not: aero_ini sits inside
!     the mypid == NROOT block above, and mpsurfgp is collective, so the read
!     cannot live there. It is a no-op unless aero_nl sets ldustemit = 1.

      if (nsela > 0 .and. l_aero > 0) call aero_surf

      call mpbci(nstep   ) ! current timestep
      call mpbci(mstep   ) ! current timestep in month
      call mpbci(ntspd   ) ! number of timesteps per day
      call mpbci(mtspd   ) ! number of timesteps per standard day

      call mpbcr(mpstep)   ! minutes per timestep
      call mpbci(n_days_per_month)
      call mpbci(n_days_per_year)
      call mpbci(m_days_per_month)
      call mpbci(m_days_per_year)
      call mpbci(mcal_days_per_year)
      call mpbci(n_steps_per_year)
      call mpbci(n_run_steps)
      call mpbci(n_run_days)
      call mpbci(n_run_months)
      call mpbci(n_run_years)
      call mpbci(n_start_step)
      call mpbci(n_start_year)
      call mpbci(n_start_month)

      call mpbci(nadv    ) !
      call mpbci(nhordif ) !
      call mpbci(nrad    ) !
      call mpbci(neqsig  ) ! switch for equidistant sigma levels
      call mpbci(nhdiff  ) ! critical wavenumber for horizonal diffusion

      call mpbci(nprint    ) ! switch for extensive diagnostic pintout (dbug)
      call mpbci(nprhor    ) ! grid point to be printed (dbug)
      call mpbci(ndiaggp   ) ! switch for franks grid point diagnostics
      call mpbci(ndiagsp   ) ! switch for franks spectral diagnostics
      call mpbci(ndiagcf   ) ! switch for cloud forcing diagnostics
      call mpbci(nenergy  )  ! switch for energy diagnostics
      call mpbci(nener3d  )  ! switch for 3d energy diagnostics
!     WITHOUT THIS THE FIXER DEADLOCKS. The namelist is read on NROOT only, so
!     an unbroadcast switch is 0 everywhere else -- and this one guards a block
!     containing `mpsumbcr`, a COLLECTIVE. NROOT enters and waits for ranks that
!     skipped it. Exactly the failure `notes/audits/nlowio-collective-deadlock.md`
!     records for `nlowio`, which is why that one is broadcast fifty lines above.
      call mpbci(nenergyfix)  ! switch for the energy fixer, world-mzy
      call mpbci(nconvtime)   ! the conversion's time level, world-0ov
      call mpbci(ndealias)    ! the conversion's dealiasing, world-ly5
      call mpbci(ndiaggp3d ) ! no of 3d gp diagnostic arrays
      call mpbci(ndiaggp2d ) ! no of 2d gp diagnostic arrays
      call mpbci(ndiagsp3d ) ! no of 3d sp diagnostic arrays
      call mpbci(ndiagsp2d ) ! no of 2d sp diagnostic arrays
      call mpbci(ntime     ) ! switch to activate time consuming estimate
      call mpbci(nperpetual) ! day of perpetual integration
      call mpbci(ndheat)     ! switch for heating due to momentum dissipation
      call mpbci(nsponge)    ! switch for top sponge layer
      call mpbci(nstratosponge) ! Switch for stratospheric newtonian cooling

      call mpbcr(acpd    )   ! Specific heat for dry air
      call mpbcr(adv     )
      call mpbcr(akap    )
      call mpbcr(alr     )
!     als, alv and tmelt are threadprivate and planet_nl now sets them
!     (p_earth.f90), and only NROOT reads that namelist. Without these three
!     every thread but the root would run the plasimmod default while the root
!     ran the configured value. world-58v.
      call mpbcr(als     )
      call mpbcr(alv     )
      call mpbcr(tmelt   )
      call mpbcr(cv      )
      call mpbcr(ct      )
      call mpbcr(dtep    )
      call mpbcr(dtns    )
      call mpbcr(dtrop   )
      call mpbcr(dttrp   )
      call mpbcr(ga      )
      call mpbcr(gascon  )
      call mpbcr(tgr     )
      call mpbcr(plarad  )
      call mpbcr(pnu     )
      call mpbcr(pnu21   )
      call mpbcr(psurf   )
      call mpbcr(ptop    )
      call mpbcr(ptop2   )
      call mpbcr(ra1     )
      call mpbcr(ra2     )
      call mpbcr(ra4     )
      call mpbcr(ra1i    ) ! the over-ice set, world-ako
      call mpbcr(ra2i    )
      call mpbcr(ra4i    )
      call mpbcr(rdbrv   )
      call mpbcr(ww      )
      call mpbcr(solar_day)
      call mpbcr(sidereal_day)
      call mpbcr(tropical_year)
      call mpbcr(sidereal_year)
      call mpbcr(rotspd)
      call mpbcr(eccen)
      call mpbcr(obliq)
      call mpbcr(mvelp)
      call mpbcr(meananom0)
      call mpbcr(plavor)
      call mpbcr(dampsp)
      call mpbcr(taucool)
      
      call mpbcr(fixedlon)
      call mpbcr(dttl)
      

      call mpbcin(ndel  ,NLEV) ! ndel
      call mpbcin(ndl   ,NLEV)

      call mpbcrn(tdissd,NLEV)
      call mpbcrn(tdissz,NLEV)
      call mpbcrn(tdisst,NLEV)
      call mpbcrn(tdissq,NLEV)
      call mpbcrn(damp  ,NLEV)
      call mpbcrn(dsigma,NLEV)
      call mpbcrn(rdsig ,NLEV)
      call mpbcrn(restim,NLEV)
      call mpbcrn(sigma ,NLEV)
      call mpbcrn(sigmah,NLEV)
      call mpbcrn(t0    ,NLEV)
      call mpbcrn(t01s2 ,NLEV)
      call mpbcrn(tfrc  ,NLEV)
      call mpbcrn(tkp   ,NLEV)

      call mpbcrn(c     ,NLSQ)
      call mpbcrn(g     ,NLSQ)
      call mpbcrn(tau   ,NLSQ)

      call mpscin(nindex,NSPP)
      call mpscsp(sak,sakpp,NLEV)
      call mpbcrn(sigh  ,NLEV)
      
      call mpbci(nfilter)
      call mpbci(ngptfilter)
      call mpbci(nspvfilter)
      call mpbcr(landhoskn0)
      call mpbci(nfilterexp)
      call mpbcr(filterkappa)
      

!     Copy some calendar variables to calmod

      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd &
                  ,solar_day,-1,mpstep,mcal_days_per_year &
                  ,m_days_per_year,m_days_per_month,mtspd)
!                  ,day_24hr,-1)

      call mpbci(mcal_days_per_year)
      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
      call updatim(nstep)  ! set date & time array ndatim
!
!     allocate additional diagnostic arrays, if switched on
!

      if(ndiaggp2d > 0) then
       allocate(dgp2d(NHOR,ndiaggp2d))
       dgp2d(:,:)=0.
      end if
      if(ndiaggp3d > 0) then
       allocate(dgp3d(NHOR,NLEV,ndiaggp3d))
       dgp3d(:,:,:)=0.
      end if
      if(ndiagsp2d > 0) then
       allocate(dsp2d(NESP,ndiagsp2d))
       dsp2d(:,:)=0.
      end if
      if(ndiagsp3d > 0) then
       allocate(dsp3d(NESP,NLEV,ndiagsp3d))
       dsp3d(:,:,:)=0.
      end if
      if(ndiagcf > 0) then
       allocate(dclforc(NHOR,7))
       dclforc(:,:)=0.
      end if
      if(nenergy > 0) then
       allocate(denergy(NHOR,28))
       denergy(:,:)=0.
       allocate(adenergy(NHOR,28))
       adenergy(:,:)=0.
      end if
      if(nener3d > 0) then
       allocate(dener3d(NHOR,NLEV,28))
       dener3d(:,:,:)=0.
       allocate(adener3d(NHOR,NLEV,28))
       adener3d(:,:,:)=0.
      end if

      call legini
!     After legini, because shtns_setup reads NTRU and the grid legini has just
!     set up, and because legmod stays the transform until nshtns says otherwise.
      if (nshtns == 1) call shtns_setup

      if (nrestart > 0) then
         call read_atmos_restart
      else
         call initfd
      endif

      if (mypid == NROOT) then
         if (noutput > 0) call outini    ! Open output file <plasim_output>
      endif
!
!     THE ECOLOGICAL STREAM, EFOR-2. Independent of noutput, because it is a
!     different consumer's product and not a subset of the climate one.
!
!     mtspd is the number of timesteps in one absolute 24-hour day exactly:
!     prolog derives it as nint(day_24hr)/nint(mpstep*60), makes it even, and
!     then recomputes mpstep so that mtspd*mpstep*60 = day_24hr. So the default
!     interval is 86400 s to the model's own arithmetic, with no rounding and no
!     drift, which is what makes the interval bounds ecogp writes exact.
!
!     24 h is the ecological step LPJ-GUESS integrates on, not a claim that this
!     world's rotation is 24 h; it is 30 h, so the local solar phase advances by
!     0.8 of a rotation every interval and stays visible in the bounds.
!     biosphere/notes/time-base-unit-contract.md settles the absolute day.
      if (neco > 0) then
         if (necostep < 1) necostep = mtspd
         if (mypid == NROOT) then
            call ecoini
!           deltsec is not set until master, so the interval is reported in
!           timesteps here and in seconds by the stream itself: ecogp writes
!           the bounds, the duration and the orbital position as its first four
!           records.
            write(nud,*) 'Ecological stream on: ',necostep,                 &
     &                   ' timesteps per interval (mtspd = ',mtspd,')'
         endif
      endif
!
!*    initialize miscellaneous additional parameterization
!     which are included in *miscmod*
!

      call miscini

!
!*    initialize surface fluxes and vertical diffusion
!

      call fluxini

!
!*    initialize radiation
!

      call radini

!
!*    initialize carbon weathering
!

      call carbonini
      
!
!*    initialize convective and large scale rain and clouds
!

      call rainini

!
!*    initialize surface parameterizations or models
!

      call surfini
      

      if (mypid==NROOT) then
         write(nud,*) "==========Finalized Albedos=========="
         write(nud,*) "-----For lambda < 0.75 microns------ "
         write(nud,*) "Ground:",dgroundalb(1)
         write(nud,*) "Ocean<;",doceanalb(1) 
         write(nud,*) "Snow:",dsnowalb(1)
         write(nud,*) "Snow max:",dsnowalbmx(1)
         write(nud,*) "Snow min:",dsnowalbmn(1)
         write(nud,*) "Sea ice max:",dicealbmx(1) 
         write(nud,*) "Sea ice min:",dicealbmn(1) 
         write(nud,*) "Glacier min:",dglacalbmn(1)
         write(nud,*) "-----For lambda > 0.75 microns------ "
         write(nud,*) "Ground:",dgroundalb(2)
         write(nud,*) "Ocean<<<;",doceanalb(2) 
         write(nud,*) "Snow:",dsnowalb(2)
         write(nud,*) "Snow max:",dsnowalbmx(2)
         write(nud,*) "Snow min:",dsnowalbmn(2)
         write(nud,*) "Sea ice max:",dicealbmx(2) 
         write(nud,*) "Sea ice min:",dicealbmn(2) 
         write(nud,*) "Glacier min:",dglacalbmn(2)
      endif      
  
!
!*    initialize hurricane/storm diagnostics
!
     
      call hurricaneini(gascon)

!
!*    reset psurf according to orography
!

      if (mypid == NROOT) then
         write(nud,'(" Surface pressure with no topography = ",f10.2," [hPa]")') psurf * 0.01
         zmeanoro = so(1) * cv * cv / sqrt(2.0)
         zdlnp    = zmeanoro / gascon / tgr
         psurf    = EXP(LOG(psurf)-zdlnp)
         write(nud,'(" Mean of topographic height          = ",f10.2," [m]")') zmeanoro / ga
         write(nud,'(" Mean of surface pressure            = ",f10.2," [hPa]")') psurf * 0.01
         write(nud,'(" Initial step number : ",i6," ")') nstep
      end if
      call mpbcr(psurf)

!
!*    broadcast and scatter
!

#ifdef OMPSHARED
!     Nothing to broadcast and nothing to scatter. The spectral state is one
!     shared array: root has already written the words the others are about to
!     read, and each partial already IS the slice it would have been filled
!     from. A broadcast here would have every thread write the same shared
!     words -- a race on identical values, so harmless in outcome, reported by
!     ThreadSanitizer, and pure traffic; a scatter would only make the compiler
!     copy a non-contiguous section in and out. One barrier does for both.
!$omp barrier
#else
      call mpbcrn(sp,NESP)
      call mpbcrn(sd,NESP*NLEV)
      call mpbcrn(st,NESP*NLEV)
      call mpbcrn(sz,NESP*NLEV)
      call mpbcrn(sq,NESP*NLEV)

      call mpscsp(sd,sdp,NLEV)
      call mpscsp(st,stp,NLEV)
      call mpscsp(sz,szp,NLEV)
      call mpscsp(sq,sqp,NLEV)
      call mpscsp(sr,srp,NLEV)
      call mpscsp(sp,spp,1)
#endif
      call mpscsp(so,sop,1)

!
!*    close the namelist file
!

      if(mypid==NROOT) close(11)

!
!*    open efficiency diagnostic file
!

      if(ndheat > 1 .and. mypid == NROOT) then
       open(9,file=efficiency_dat,form='formatted')
      endif

!
!*    start time consuming calculations
!

      if(ntime==1) then
       call mksecond(zsec,0.)
       time0=zsec
      endif

      return
      end


!     =================
!     SUBROUTINE MASTER
!     =================

      subroutine master
      use pumamod

!     ***************************
!     * short initial timesteps *
!     ***************************

      if (mypid == NROOT) write(nud,*) "Preparing to do initial timesteps: ",nkits

      ikits = nkits
      do jkits=1,ikits
         deltsec  = (day_24hr / mtspd) / (2**nkits) !Timestep
         deltsec2 = deltsec + deltsec
         delt     = deltsec * ww       !Timestep in the model's own time unit
         delt2    = delt + delt
         if (mypid == NROOT) then
            write(nud,*) 'Initial timestep ',jkits,'   deltsec = ',deltsec
            write(nud,*) 'Initial timestep ',jkits,'   delt    = ',delt
         endif
         call gridpointa
         call makebm
         call spectrala
         call gridpointd
         call spectrald
         nkits = nkits - 1
      enddo
      
      if (mypid == NROOT) write(nud,*) "Dynamical core is now spun-up."

!     ****************************************************************
!     * The scaling factor "ww" is derived from the rotation "omega" *
!     * with 1 planetary rotation per sidereal day (2 Pi) of Earth   *
!     ****************************************************************

!     THE NONDIMENSIONAL TIMESTEP IS deltsec*ww AND NOT TWOPI/ntspd. world-r8o.
!
!     `ntspd` counts timesteps per SOLAR day and `ww` is built from the SIDEREAL
!     one, so `TWOPI/ntspd` nondimensionalises the timestep against a clock the
!     rest of the model does not use. The two differ by
!     n_days_per_year/(n_days_per_year-1), which is 0.69% here.
!
!     It did not show, because `ntspd = nint(solar_day)/nint(mpstep*60)` is an
!     INTEGER division: at every timestep on the rung table it truncates to
!     exactly `sidereal_day/deltsec` and the error is annulled. The truncation
!     stops annulling it above n_days_per_year-1 steps per day, which is below
!     dt = 12.4 min -- inside the range the T127 and T170 rungs are headed for.
!     Written this way the identity holds at every timestep instead of at the
!     ones where an integer division happens to agree, which is also what the
!     declaration of `delt` in plasimmod.f90 says it is.
!
      deltsec  = day_24hr / mtspd   ! timestep in seconds
      deltsec2 = deltsec + deltsec   ! timestep in seconds * 2
      delt     = deltsec * ww        ! timestep in the model's own time unit
      delt2    = delt + delt
      call makebm
!
!     NCONVTIME TAKES THE TEMPERATURE EQUATION'S REFERENCE CONVERSION OUT OF THE
!     SEMI-IMPLICIT TREATMENT, so the timestep it is stable at is the EXPLICIT
!     gravity-wave one and not the rung table's. For a spectral model the fastest
!     resolved external mode gives
!
!         dt  <  a / (c sqrt(N (N+1))),    c = sqrt(R T0 / (1 - kappa))
!
!     which at T42 on this planet is 9.5 minutes against a configured 22.5. Run
!     above it and the model blows up inside ten model days, which is what
!     happened. This refuses rather than letting a declared setting integrate
!     something that is not a solution. world-0ov.
!
!     NDEALIAS USES THE LEGENDRE PATH'S DECOMPOSITION -- a per-process partial
!     from fc2sp_t, reduce-scattered and gathered back -- and SHTns integrates
!     the whole globe and returns a finished field instead. The two are not
!     interchangeable, so this refuses rather than transforming through a path
!     whose partials mean something else. world-ly5.
      if (ndealias > 0 .and. nshtns == 1) then
         if (mypid == NROOT) write(nud,*)                                &
     &      'NDEALIAS has no SHTns path; run the MPI build. world-ly5'
         stop 'ndealias not implemented on the SHTns transform path'
      endif
      if (nconvtime > 0) then
         zcgw  = sqrt(gascon * t0(NLEV) * ct / (1.0 - akap))
         zcgwd = plarad / (zcgw * sqrt(real(NTRU) * real(NTRU+1)))
         if (mypid == NROOT) then
            write(nud,'(A,F8.1,A,F8.1,A)')                              &
     &         ' NCONVTIME: explicit gravity-wave timestep limit ',      &
     &         zcgwd/60.0,' min, this run runs at ',deltsec/60.0,' min'
         endif
         if (deltsec > zcgwd) then
            if (mypid == NROOT) write(nud,*)                            &
     &         'NCONVTIME needs a timestep at or below the explicit ',   &
     &         'gravity-wave limit; see world-0ov'
            stop 'nconvtime above the explicit gravity-wave timestep'
         endif
      endif

      if (mypid == NROOT .and. nsela > 0) then
         call tracer_ini
      endif

!     Use either month countdown (n_run_years * 12 + n_run_months)
!     or step-countdown (for debugging purposes)
      
      mocd = n_run_months                     ! month countdown
      nscd = n_run_days * mtspd + n_run_steps ! step countdown

      if (nrestart == 0) nstep = n_start_step ! timestep since 01-01-0001
      call updatim(nstep)  ! set date & time array ndatim

      nstep1 = nstep ! Remember start step for timing stats
      
      nhcstp = 1
      
!       if (mypid == NROOT) write(nud,'("*  nstep1: ",i6,"   *")') nstep1

      do while (mocd > 0 .or. nscd > 0)  ! main loop

 
!        ************************************************************
!        * calculation of non-linear quantities in grid point space *
!        ************************************************************

         call gridpointa

!        ******************************
!        * adiabatic part of timestep *
!        ******************************

         call spectrala

!        *****************************
!        * diabatic part of timestep *
!        *****************************

         call gridpointd

!          call outaccu
!          naccuout = naccuout + 1
         if (mypid == NROOT) then
            if (mod(nhcstp,nafter) == 0 .and. noutput > 0 ) then
               call outsp
            endif
            if (mod(nhcstp,nstps) == 0 .and. nsnapshot > 0) call snapshotsp
            if (nhcadence>0 .and. hcstartstep>0 .and. hcendstep>0) then
            !We use this triple condition to be really sure that this doesn't get
            !turned on by accident--this has the potential to not only create huge
            !amounts of output if misused, but also to actually damage computing infrastructure.
              if ((nhcstp .ge. hcstartstep) .and. (nhcstp<hcendstep)) then
                if (mod(nhcstp-hcstartstep,hcinterval)==0) call hcadencesp(141)
              endif
            endif
            if (nwritehurricane>0 .and. mod(nhcstp,hcinterval)==0) call hcadencesp(142)
            
            if (mod(nstep,ndiag) == 0 ) then
               call diag
            elseif (ngui > 0) then
               if (mod(nstep,ngui) == 0 ) call diag
            endif
         endif

         if (ngui > 0) call guistep_plasim
         call spectrald
         call outaccu
         if (neco > 0) call ecoaccu

         if (mod(nhcstp,nstps) == 0 .and. nsnapshot > 0) then
           call snapshotsc
           call snapshotgp
           koutdiag=ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf     &
     &             +nenergy
           if(koutdiag > 0) call snapshotdiag
         endif
         if (nhcadence>0 .and. hcstartstep>0 .and. hcendstep>0) then
            !We use this triple condition to be really sure that this doesn't get
            !turned on by accident--this has the potential to not only create huge
            !amounts of output if misused, but also to actually damage computing infrastructure.
           if ((nhcstp .ge. hcstartstep) .and. (nhcstp<hcendstep)) then
             if (mypid == NROOT) write(nud,*) "HC OUTPUT step",nhcstp
             if (mod(nhcstp-hcstartstep,hcinterval)==0) then
!               ONCE, NOT TWICE. Upstream calls this twice in a row -- it is in
!               the first squashed import of the subtree, so it has been there as
!               long as the fork has. Every gridpoint high-cadence record was
!               written twice, and the spectral half beside it once, so pyburn
!               derives the sample count from one and the array from the other
!               and cannot reshape: "cannot reshape array of size 7402780 into
!               shape (2926,10,506)", and 7402780 is exactly half of 2926*10*506.
!               The storm-capture block below writes one hcadencesp and one
!               hcadencegp, which is the intended pattern.
                call hcadencegp(141)
!               NO SCALAR AND NO DIAGNOSTIC RECORD, and that is the decision
!               rather than a gap. The regular stream calls outsc beside outgp
!               and the snapshot stream calls snapshotsc beside snapshotgp,
!               because both are read as climate fields and want the orbital
!               phase alongside them. This stream is not: it exists for the gust
!               distribution DUST-5 needs, its postprocessed field list is the
!               winds, its time axis comes from the code 139 record that
!               hcadencegp already writes, and one orbit of it is 15 GB before
!               anything is added. `hcadencesc` and `hcadencediag` are gone with
!               the commented-out call that used to stand here; a stream that
!               wants the orbital scalars can call outsc's codes itself.
!               world-4vp.
             endif
           endif
         endif
         if (nwritehurricane>0 .and. mod(nhcstp,hcinterval)==0) then
            if (mypid == NROOT) write(nud,*) "HC STORM CAPTURE step",nhcstp
            call hcadencesp(142)
            call hcadencegp(142)
         endif
         if (neco > 0) then
            if (mod(naccueco,necostep) == 0) then
               call ecogp
               call ecoreset
            endif
         endif
         if (mod(nhcstp,nafter) == 0) then
          if(noutput > 0) then
!            write(nud,*) "High-cadence step",nhcstp
           call outsc
           call outgp
           koutdiag=ndiaggp3d+ndiaggp2d+ndiagsp3d+ndiagsp2d+ndiagcf     &
     &             +nenergy
           if(koutdiag > 0) call outdiag
          endif
          call outreset
!           write(nud,*)"Called outreset"
         endif

         iyea  = ndatim(1)    ! current year
         imon  = ndatim(2)    ! current month
         nstep = nstep + 1
         nhcstp = nhcstp + 1
         mstep = mstep + 1
         call updatim(nstep)  ! set date & time array ndatim
         if (imon /= ndatim(2)) then
            mocd = mocd - 1 ! next month
            if (mypid == NROOT) then
               write(nud,"('Completed month ',I2.2,'-',I4.4)") imon,iyea  
            endif
            mstep = 0
         endif
         call stability_check ! aborts if model tends to explode
         if (nshutdown   > 0) return
         if (nscd  > 0) then
            nscd = nscd - 1
            if (nscd == 0) return
         endif
      enddo ! month countdown

      return
      end


!     =================
!     SUBROUTINE EPILOG
!     =================

      subroutine epilog
      use pumamod
      real    (kind=8) :: zut,zst
      integer (kind=8) :: imem,ipr,ipf,isw,idr,idw
!
!     close output file
!
      if (mypid == NROOT) close(40)
      
!
!     close snapshot file
!
      if (nsnapshot > 0 .and. mypid == NROOT) close(140)
      if (nhcadence > 0 .and. mypid == NROOT) close(141)
      if (neco > 0 .and. mypid == NROOT) close(143)
!
!     close efficiency diagnostic file
!
      if(ndheat > 1 .and. mypid == NROOT) close(9)
!
!     write restart file
!
      if (mypid == NROOT) then
         call restart_prepare(plasim_status)
         
         if (mod(naccuout,nstpw)==0) naccuout = 0
         
         call put_restart_integer('nstep'   ,nstep   )
         call put_restart_integer('naccuout',naccuout)
         call put_restart_integer('naccueco',naccueco)
         call put_restart_integer('nlat'    ,NLAT    )
         call put_restart_integer('nlon'    ,NLON    )
         call put_restart_integer('nlev'    ,NLEV    )
         call put_restart_integer('nrsp'    ,NRSP    )

!        Save current random number generator seed

         call random_seed(get=meed)
         call put_restart_seed('seed',meed,nseedlen)

         call put_restart_array('sz',sz,NRSP,NESP,NLEV)
         call put_restart_array('sd',sd,NRSP,NESP,NLEV)
         call put_restart_array('st',st,NRSP,NESP,NLEV)
         if (nqspec == 1) call put_restart_array('sq',sq,NRSP,NESP,NLEV)
         call put_restart_array('sr',sr,NRSP,NESP,NLEV)
         call put_restart_array('sp',sp,NRSP,NESP,   1)
         call put_restart_array('so',so,NRSP,NESP,   1)
      endif

      call mpputsp('szm',szm,NSPP,NLEV)
      call mpputsp('sdm',sdm,NSPP,NLEV)
      call mpputsp('stm',stm,NSPP,NLEV)
      if (nqspec == 1) call mpputsp('sqm',sqm,NSPP,NLEV)
      call mpputsp('spm',spm,NSPP,   1)
!
!     gridpoint restart
!
      call mpputgp('dls'    ,dls    ,NHOR,1)
      call mpputgp('drhs'   ,drhs   ,NHOR,1)
      call mpputgp('dalb'   ,dalb   ,NHOR,1)
      call mpputgp('dsalb1',dsalb(1,:),NHOR,1)
      call mpputgp('dsalb2',dsalb(2,:),NHOR,1)
      call mpputgp('dz0'    ,dz0    ,NHOR,1)
      call mpputgp('dicec'  ,dicec  ,NHOR,1)
      call mpputgp('diced'  ,diced  ,NHOR,1)
      call mpputgp('dwatc'  ,dwatc  ,NHOR,1)
      call mpputgp('drunoff',drunoff,NHOR,1)
      call mpputgp('dforest',dforest,NHOR,1)
      call mpputgp('dust3'  ,dust3  ,NHOR,1)
      call mpputgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpputgp('dql'    ,dql    ,NHOR,NLEV)
      call mpputgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpputgp('dt'     ,dt(1,NLEP),NHOR,1)
      if (nqspec == 1) then ! spectral: save only soil humidity
         call mpputgp('dq'  ,dq(1,NLEP),NHOR,1)
      else                  ! semi-langrange: save complete humidity array
         call mpputgp('dq'  ,dq,NHOR,NLEP)
      if (l_aero > 0) then
         call mpputgp('mmr' ,mmr,NHOR,NLEP)
         call mpputgp('nrho',nrho,NHOR,NLEP)
      endif
      endif
!
!     accumulated diagnostics
!
      call mpputgp('aprl'  ,aprl  ,NHOR,1)
      call mpputgp('aprc'  ,aprc  ,NHOR,1)
      call mpputgp('aprs'  ,aprs  ,NHOR,1)
      call mpputgp('aevap' ,aevap ,NHOR,1)
      call mpputgp('ashfl' ,ashfl ,NHOR,1)
      call mpputgp('alhfl' ,alhfl ,NHOR,1)
      call mpputgp('aroff' ,aroff ,NHOR,1)
      call mpputgp('asmelt',asmelt,NHOR,1)
      call mpputgp('asndch',asndch,NHOR,1)
      call mpputgp('acc'   ,acc   ,NHOR,1)
      call mpputgp('assol' ,assol ,NHOR,1)
      call mpputgp('asthr' ,asthr ,NHOR,1)
      call mpputgp('atsol' ,atsol ,NHOR,1)
      call mpputgp('atthr' ,atthr ,NHOR,1)
      call mpputgp('ataux' ,ataux ,NHOR,1)
      call mpputgp('atauy' ,atauy ,NHOR,1)
      call mpputgp('atsolu',atsolu,NHOR,1)
      call mpputgp('assolu',assolu,NHOR,1)
!     The surface downward solar flux by band, WORLD-3QFZ. Written every time,
!     read only behind the accumulator marker below, which is why that marker
!     had to move: a restart written before these existed carries a counter
!     these two cannot match, and a diluted first output window is exactly the
!     defect `accuvers` exists to refuse.
      call mpputgp('afdsw1',afdsw1,NHOR,1)
      call mpputgp('afdsw2',afdsw2,NHOR,1)
      call mpputgp('asthru',asthru,NHOR,1)
      call mpputgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpputgp('atsa'  ,atsa  ,NHOR,1)
      call mpputgp('ats0'  ,ats0  ,NHOR,1)
      call mpputgp('atsama',atsama,NHOR,1)
      call mpputgp('atsami',atsami,NHOR,1)
      
      call mpputgp('azmuz'       ,azmuz   ,NHOR,1) 
                                           
      call mpputgp('asigrain'    ,asigrain,NHOR,1)
      call mpputgp('tempmax'     ,tempmax ,NHOR,1) 
      
!     THE ECOLOGICAL STREAM'S PARTIAL INTERVAL, EFOR-2.
!
!     Serialized unconditionally, whether or not the stream is on. An
!     accumulator saved only while its switch is set is a switch that changes
!     the restart's contents, and restart_schema.py's inventory is over the
!     model's call sites and not over a namelist.
!
!     WHY IT HAS TO BE SAVED AT ALL: without it a continuation restarts the
!     interval it was in the middle of, so the first block after a model call
!     covers a shorter span than it declares. That is not a small error in one
!     record; it puts a false weather boundary into the sequence, at exactly the
!     place a replay protocol later tests for a seam. The counter goes with the
!     accumulators, because a counter that outlives what it counts is the defect
!     the regular stream's `accuvers` marker exists for.
      call mpputgp('aecotas'    ,aecotas    ,NHOR,1)
      call mpputgp('aecots'     ,aecots     ,NHOR,1)
      call mpputgp('aecops'     ,aecops     ,NHOR,1)
      call mpputgp('aecohus'    ,aecohus    ,NHOR,1)
      call mpputgp('aecowind'   ,aecowind   ,NHOR,1)
      call mpputgp('aecoswd'    ,aecoswd    ,NHOR,1)
      call mpputgp('aecoswu'    ,aecoswu    ,NHOR,1)
      call mpputgp('aecoswn'    ,aecoswn    ,NHOR,1)
      call mpputgp('aecolwn'    ,aecolwn    ,NHOR,1)
      call mpputgp('aecolwu'    ,aecolwu    ,NHOR,1)
      call mpputgp('aecoczen'   ,aecoczen   ,NHOR,1)
      call mpputgp('aecopr'     ,aecopr     ,NHOR,1)
      call mpputgp('aecoprsn'   ,aecoprsn   ,NHOR,1)
      call mpputgp('aecoprc'    ,aecoprc    ,NHOR,1)
      call mpputgp('aecoevap'   ,aecoevap   ,NHOR,1)
      call mpputgp('aecotasmx'  ,aecotasmx  ,NHOR,1)
      call mpputgp('aecotasmn'  ,aecotasmn  ,NHOR,1)
      call mpputgp('aecotsmx'   ,aecotsmx   ,NHOR,1)
      call mpputgp('aecotsmn'   ,aecotsmn   ,NHOR,1)
      call mpputgp('tempmin'     ,tempmin ,NHOR,1) 
                                           
!     The six accumulators below are SPECTRAL (NESP,NLEV) and are meaningful
!     only on NROOT, which is the only task that writes them out. They used to
!     go through mpputgp/mpgetgp, which move NHOR elements per level and so
!     transferred NHOR*NLEV of NESP*NLEV values: levels beyond the first
!     NHOR*NLEV/NESP came back as zero while naccuout came back whole, and the
!     first record of every model call was scaled by (nstpw-1)/naccuout on
!     those levels. Save them the way sz/sd/st/sp/so are saved, and under new
!     names so that a pre-patch restart is recognised rather than misread.
      if (mypid == NROOT) then
         call put_restart_array('aasosp' ,aaso   ,NESP,NESP,   1)
         call put_restart_array('aaspsp' ,aasp   ,NESP,NESP,   1)
         call put_restart_array('aastsp' ,aast   ,NESP,NESP,NLEV)
         call put_restart_array('aasqsp' ,aasqout,NESP,NESP,NLEV)
         call put_restart_array('aasdsp' ,aasd   ,NESP,NESP,NLEV)
         call put_restart_array('aaszsp' ,aasz   ,NESP,NESP,NLEV)
!        Accumulated orbital scalars: not saved at all before this patch.
         call put_restart_real('aorbnu'  ,aorbnu)
         call put_restart_real('alambm'  ,alambm)
         call put_restart_real('azdecl'  ,azdecl)
         call put_restart_real('ardist'  ,ardist)
         call put_restart_real('arasc'   ,arasc )
!        Marks a restart whose accumulator set is complete. 2.0 adds the two
!        band-resolved surface downward solar accumulators, WORLD-3QFZ.
         call put_restart_real('accuvers',2.0)
!        Marks a restart that carries the ecological stream's partial interval.
         call put_restart_real('ecovers',1.0)
!        THE ENERGY FIXER'S INTEGRATED CORRECTION. world-fsr.
!
!        `denergyfix` is a controller state, not a diagnostic: it is the uniform
!        heating currently being applied, reached in one window and then tracked.
!        Left out of the restart it returned to zero at every segment boundary,
!        so the first 2*ntspd steps of every segment ran with NO correction and
!        then stepped to the full value -- and any run shorter than one window,
!        which includes the 20-step transform gates, never applied one at all.
!
!        The three window accumulators are deliberately NOT saved beside it. The
!        first step out of a restart carries a start-up transient of order 250
!        W/m2 and the design discards the first window for exactly that reason;
!        splicing a partial window across a segment boundary would feed that
!        transient into the previous segment's average. A segment starts on the
!        previous segment's converged correction and measures its own first
!        clean window before changing it, which is what the design asks for.
         call put_restart_real('denergyfix',denergyfix)
      endif
!     THE ENERGY DIAGNOSTICS' PARTIAL WINDOW, world-5qy. `adenergy` and
!     `adener3d` are extended by `outaccu` every timestep and divided by
!     `naccuout` in `outgp`, exactly as every other accumulator is, and they were
!     the only two the restart did not carry. They are carried now: every restart
!     on disk is written mid-window, so the first output record after each resume
!     was dividing what had accumulated since the resume by a whole window.
!
!     WRITTEN ONLY WHERE THEY ARE ALLOCATED, which `nenergy` and `nener3d`
!     decide, so the record set follows the switch. That is why the `accuvers`
!     marker cannot stand for them: it promises a complete set, and a set whose
!     membership is a namelist question cannot be complete. Their presence is
!     ASKED on the read side instead, and a resume that wants them and does not
!     find them discards the partial interval the way a pre-marker restart does.
      if (nenergy > 0) call mpputgp('adenergy',adenergy,NHOR,28)
      if (nener3d > 0) call mpputgp('adener3d',adener3d,NHOR,NLEV*28)
!     Accumulated hurricane indices are gridpoint fields divided by naccuout in
!     outgp; they were not saved either.
      call mpputgp('agpi'         ,agpi        ,NHOR,1)
      call mpputgp('aventi'       ,aventi      ,NHOR,1)
      call mpputgp('alaav'        ,alaav       ,NHOR,1)
      call mpputgp('ampoti'       ,ampoti      ,NHOR,1)
      call mpputgp('avrmpi'       ,avrmpi      ,NHOR,1)
      call mpputgp('acapen'       ,acapen      ,NHOR,1)
      call mpputgp('alnb'         ,alnb        ,NHOR,1)
      call mpputgp('achim'        ,achim       ,NHOR,1)
      call mpputgp('aadq'        ,aadq    ,NHOR,NLEP) 
      call mpputgp('aammr'       ,aammr   ,NHOR,NLEP)
      call mpputgp('aanrho'      ,aanrho  ,NHOR,NLEP)
      call mpputgp('aadmld'      ,aadmld  ,NHOR,1)      
      call mpputgp('aadt'        ,aadt    ,NHOR,NLEP)   
      call mpputgp('aadwatc'     ,aadwatc ,NHOR,1)     
      call mpputgp('aadsnow'     ,aadsnow ,NHOR,1)     
      call mpputgp('aadql'       ,aadql   ,NHOR,NLEP)  
      call mpputgp('aadust3'     ,aadust3 ,NHOR,1)     
      call mpputgp('aadcc'       ,aadcc   ,NHOR,NLEP)  
      call mpputgp('aadtd5'      ,aadtd5  ,NHOR,1)      
      call mpputgp('aadls'       ,aadls   ,NHOR,1)       
      call mpputgp('aadz0'       ,aadz0   ,NHOR,1)       
      call mpputgp('aadalb'      ,aadalb  ,NHOR,1)      
      call mpputgp('aadsalb1'    ,aadsalb1,NHOR,1)    
      call mpputgp('aadsalb2'    ,aadsalb2,NHOR,1)    
      call mpputgp('aadtsoil'    ,aadtsoil,NHOR,1)    
      call mpputgp('aadtd2'      ,aadtd2  ,NHOR,1)      
      call mpputgp('aadtd3'      ,aadtd3  ,NHOR,1)      
      call mpputgp('aadtd4'      ,aadtd4  ,NHOR,1)      
      call mpputgp('aadicec'     ,aadicec ,NHOR,1)     
      call mpputgp('aadiced'     ,aadiced ,NHOR,1)     
      call mpputgp('aadforest'   ,aadforest   ,NHOR,1)   
      call mpputgp('aadwmax'     ,aadwmax     ,NHOR,1)     
      call mpputgp('aadglac'     ,aadglac     ,NHOR,1)     
      call mpputgp('aadqo3'      ,aadqo3      ,NHOR,NLEV) 
      call mpputgp('aagroundoro' ,aagroundoro ,NHOR,1) 
      call mpputgp('aaglacieroro',aaglacieroro,NHOR,1)

!
!*    finish graphical user interface
!

      call guistop

!
!*    finish miscellaneous additional parameterizations
!

      call miscstop

!
!*    finish surface fluxes and vertical diffusion
!

      call fluxstop

!
!*    finish radiation
!

      call radstop

!
!*    finish large scale and convective rain and clouds
!

      call rainstop

!
!*    finish surface parameterizations or models
!

      call surfstop

!
!*    finish carbon-silicate weathering parameterizations
!

      call carbonstop
      
      if (mypid == NROOT) then
         call restart_stop
      endif

!
!     DEALLOCATE THE DIAGNOSTIC ARRAYS, AND NOT ONE LINE EARLIER. This block sat
!     at the TOP of epilog until WORLD-0OV's arms found what that costs: since
!     `adenergy` and `adener3d` became restart records (world-5qy), freeing them
!     here and writing them at `mpputgp` below is a use-after-free, and it fires
!     at EVERY orbit boundary rather than once at the end -- so the standing
!     production configuration, which sets `energy_fixer` and therefore requires
!     `energy_diagnostics`, could not complete a single orbit.
!
!     The whole block moved rather than the two accumulators, and that is the
!     point: the order is now STRUCTURAL. Nothing epilog frees can precede the
!     restart, so the next array promoted to a restart record cannot re-create
!     this by being added to a list that looks unrelated. Freeing after the write
!     costs nothing here -- epilog is the end of the run.
!
      if(ndiaggp2d > 0) deallocate(dgp2d)
      if(ndiagsp2d > 0) deallocate(dsp2d)
      if(ndiaggp3d > 0) deallocate(dgp3d)
      if(ndiagsp3d > 0) deallocate(dsp3d)
      if(ndiagcf   > 0) deallocate(dclforc)
      if(nenergy   > 0) deallocate(denergy)
      if(nenergy   > 0) deallocate(adenergy)
      if(nener3d   > 0) deallocate(dener3d)
      if(nener3d   > 0) deallocate(adener3d)

      call hurricanestop
!
!     time consumption
!
      if (mypid == NROOT) then 
!        Get resource stats from function resources in file pumax.c
         ires = nresources(zut,zst,imem,ipr,ipf,isw,idr,idw)
         call cpu_time(tmstop)
         tmrun = tmstop - tmstart
         if (nstep > nstep1) then 
            zspy = tmrun * m_days_per_year * real(mtspd) / (nstep-nstep1)
            zypd = (24.0 * 3600.0 / zspy)                         ! siy / day
            write(nud,'(/,"****************************************")')
            if (zut > 0.0) &
            write(nud,  '("* User   time         : ", f10.3," sec *")') zut
            if (zst > 0.0) &
            write(nud,  '("* System time         : ", f10.3," sec *")') zst
            if (zut + zst > 0.0) tmrun = zut + zst
            write(nud,  '("* Total CPU time      : ", f10.3," sec *")') tmrun
            if (imem > 0) then
               zmem = imem * 0.000001
               if (zmem < 1.0) zmem = 1000.0 * zmem ! Could be KB or MB
            write(nud,  '("* Memory usage        : ", f10.3," MB  *")') zmem
            endif
            if (ipr > 0) &
            write(nud,  '("* Page reclaims       :", i7," pages   *")') ipr
            if (ipf > 0) &
            write(nud,  '("* Page faults         :", i7," pages   *")') ipf
            if (isw > 0) &
            write(nud,  '("* Page swaps          :", i7," pages   *")') isw
            if (idr > 0) &
            write(nud,  '("* Disk read           :", i7," blocks  *")') idr
            if (idw > 0) &
            write(nud,  '("* Disk write          :", i7," blocks  *")') idw
            write(nud,'("****************************************")')
            if (zspy < 600.0) then
               write(nud,'("* Seconds per sim year: ",i6,9x,"*")') nint(zspy)
            else if (zspy < 900000.0) then
               write(nud,'("* Minutes per sim year  ",i6,9x,"*")') nint(zspy/60.0)
            else    
               write(nud,'("* Days per sim year:    ",i6,5x,"*")') nint(zspy/86400.0)
            endif
               write(nud,'("* Sim years per day   :",i7,9x,"*")') nint(zypd)
            write(nud,'("****************************************")')
         endif   
      endif      
      return
      end


!     =============================
!     SUBROUTINE READ_ATMOS_RESTART
!     =============================

      subroutine read_atmos_restart
      use pumamod
      use restartmod, only: nexcheck

!     THE STAMPED GEOMETRY IS READ BEFORE ANY ARRAY IS. world-4yz.
!
!     `epilog` writes nlat, nlon, nlev and nrsp into every restart and nothing
!     read them back. `get_restart_array` is `read (nreaunit) pa(1:k1,:)` with
!     no iostat, and an unformatted sequential read that consumes FEWER values
!     than the record holds is legal and advances: a donor written at a higher
!     truncation is silently reinterpreted, the target taking the leading k1*k3
!     values of a longer record, which runs out of the donor's level one into
!     its tail rather than truncating the field. A lower-truncation donor hits
!     end of record and dies with no iostat to say why.
!
!     `exoplasim/scripts/restart_schema.py` validates all four, but that is the
!     CONVERTER: an unconverted restart from another rung passed every gate the
!     run path had. This is the run path's own gate, and it is placed above the
!     first array read because a mis-sized read has already done its damage by
!     the time anything downstream could notice.
!
!     A restart written before this project stamped them is possible in
!     principle, so nexcheck is lowered for the four reads and a value that did
!     not arrive is left at its sentinel and passes. An UNSTAMPED restart is not
!     a MISMATCHED one.

      nresbad = 0
      if (mypid == NROOT) then
         jrlat = -1
         jrlon = -1
         jrlev = -1
         jrrsp = -1
         nexcheck = 0
         call get_restart_integer('nlat',jrlat)
         call get_restart_integer('nlon',jrlon)
         call get_restart_integer('nlev',jrlev)
         call get_restart_integer('nrsp',jrrsp)
         nexcheck = 1
         if ((jrlat > 0 .and. jrlat /= NLAT) .or.                        &
     &       (jrlon > 0 .and. jrlon /= NLON) .or.                        &
     &       (jrlev > 0 .and. jrlev /= NLEV) .or.                        &
     &       (jrrsp > 0 .and. jrrsp /= NRSP)) then
            nresbad = 1
            write(nud,*) '*** RESTART GEOMETRY MISMATCH ***'
            write(nud,*) 'restart NLAT NLON NLEV NRSP: ',                &
     &                   jrlat,jrlon,jrlev,jrrsp
            write(nud,*) 'this binary NLAT NLON NLEV NRSP: ',            &
     &                   NLAT,NLON,NLEV,NRSP
            write(nud,*) 'the spectral records are sized by the writing'
            write(nud,*) 'truncation and would be resliced, not projected.'
            write(nud,*) 'Convert it: exoplasim/scripts/convert_restart.py'
            write(nud,*) 'see world-4yz'
         endif
      endif
!     EVERY TASK STOPS, not just NROOT. The reads above are root-only and the
!     next collective is the broadcast below, so a root-only stop would leave
!     the others in a collective nobody is going to enter.
      call mpbci(nresbad)
      if (nresbad > 0) stop 'restart written at a different resolution'

!     read scalars and full spectral arrays

      if (mypid == NROOT) then
         call random_seed(size=nseedlen)
         allocate(meed(nseedlen))
         call get_restart_integer('nstep'   ,nstep)
         call get_restart_integer('naccuout',naccuout)
         call get_restart_seed('seed',meed,nseedlen)
         call get_restart_array('sz',sz,NRSP,NESP,NLEV)
         call get_restart_array('sd',sd,NRSP,NESP,NLEV)
         call get_restart_array('st',st,NRSP,NESP,NLEV)
         if (nqspec == 1) call get_restart_array('sq',sq,NRSP,NESP,NLEV)
         call get_restart_array('sr',sr,NRSP,NESP,NLEV)
         call get_restart_array('sp',sp,NRSP,NESP,   1)
         call get_restart_array('so',so,NRSP,NESP,   1)
         call random_seed(put=meed)
      endif

      call mpbci(nstep)     ! broadcast current timestep
      call mpbci(naccuout)  ! broadcast accumulation timestep for diagnostics

!     read and scatter spectral arrays

      call mpgetsp('szm',szm,NSPP,NLEV)
      call mpgetsp('sdm',sdm,NSPP,NLEV)
      call mpgetsp('stm',stm,NSPP,NLEV)
      if (nqspec == 1) call mpgetsp('sqm',sqm,NSPP,NLEV)
      call mpgetsp('spm',spm,NSPP,   1)

!     read and scatter surface grids

      call mpgetgp('dls'    ,dls    ,NHOR,   1)
      call mpgetgp('drhs'   ,drhs   ,NHOR,   1)
      call mpgetgp('dalb'   ,dalb   ,NHOR,   1)
      call mpgetgp('dsalb1' ,dsalb(1,:),NHOR,1)
      call mpgetgp('dsalb2' ,dsalb(2,:),NHOR,1)
      call mpgetgp('dz0'    ,dz0    ,NHOR,   1)
      call mpgetgp('dicec'  ,dicec  ,NHOR,   1)
      call mpgetgp('diced'  ,diced  ,NHOR,   1)
      call mpgetgp('dwatc'  ,dwatc  ,NHOR,   1)
      call mpgetgp('drunoff',drunoff,NHOR,   1)
      call mpgetgp('dforest',dforest,NHOR,   1)
      call mpgetgp('dust3'  ,dust3  ,NHOR,   1)
      call mpgetgp('dcc'    ,dcc    ,NHOR,NLEV)
      call mpgetgp('dql'    ,dql    ,NHOR,NLEV)
      call mpgetgp('dqsat'  ,dqsat  ,NHOR,NLEV)
      call mpgetgp('dt'     ,dt(1,NLEP),NHOR,1)
      if (nqspec == 1) then ! spectral: read only soil humidity
         call mpgetgp('dq',dq(1,NLEP),NHOR,1)
      else                  ! semi-langrange: read complete humidity array
         call mpgetgp('dq',dq,NHOR,NLEP)
      if (l_aero > 0) then
         call mpgetgp('mmr',mmr,NHOR,NLEP)
         call mpgetgp('nrho',nrho,NHOR,NLEP)
      endif
      endif

!     read and scatter accumulated diagnostics

      call mpgetgp('aprl'  ,aprl  ,NHOR,1)
      call mpgetgp('aprc'  ,aprc  ,NHOR,1)
      call mpgetgp('aprs'  ,aprs  ,NHOR,1)
      call mpgetgp('aevap' ,aevap ,NHOR,1)
      call mpgetgp('ashfl' ,ashfl ,NHOR,1)
      call mpgetgp('alhfl' ,alhfl ,NHOR,1)
      call mpgetgp('aroff' ,aroff ,NHOR,1)
      call mpgetgp('asmelt',asmelt,NHOR,1)
      call mpgetgp('asndch',asndch,NHOR,1)
      call mpgetgp('acc'   ,acc   ,NHOR,1)
      call mpgetgp('assol' ,assol ,NHOR,1)
      call mpgetgp('asthr' ,asthr ,NHOR,1)
      call mpgetgp('atsol' ,atsol ,NHOR,1)
      call mpgetgp('atthr' ,atthr ,NHOR,1)
      call mpgetgp('ataux' ,ataux ,NHOR,1)
      call mpgetgp('atauy' ,atauy ,NHOR,1)
      call mpgetgp('atsolu',atsolu,NHOR,1)
      call mpgetgp('assolu',assolu,NHOR,1)
      call mpgetgp('asthru',asthru,NHOR,1)
      call mpgetgp('aqvi'  ,aqvi  ,NHOR,1)
      call mpgetgp('atsa'  ,atsa  ,NHOR,1)
      call mpgetgp('ats0'  ,ats0  ,NHOR,1)
      call mpgetgp('atsama',atsama,NHOR,1)
      call mpgetgp('atsami',atsami,NHOR,1)
       
      call mpgetgp('azmuz'       ,azmuz   ,NHOR,1) 
                                           
      call mpgetgp('asigrain'    ,asigrain,NHOR,1)
      call mpgetgp('tempmax'     ,tempmax ,NHOR,1) 
      call mpgetgp('tempmin'     ,tempmin ,NHOR,1) 
      
!     THE ECOLOGICAL STREAM'S PARTIAL INTERVAL, EFOR-2, behind its own version
!     marker for the same reason `accuvers` exists. mpgetgp passes an
!     UNINITIALISED buffer to get_restart_array and scatters it whatever
!     happens, so reading a record that is not there under a lowered nexcheck
!     would scatter garbage rather than leave the array alone. A restart written
!     before this stream existed therefore starts the interval clean and says
!     so, which is what it would have done anyway.
      zecovers = -1.0
      if (mypid == NROOT) then
         nexcheck = 0
         call get_restart_real('ecovers',zecovers)
         call get_restart_integer('naccueco',naccueco)
         nexcheck = 1
      endif
      call mpbcr(zecovers)
      call mpbci(naccueco)
      if (zecovers < 0.0) then
         call ecoreset
         if (mypid == NROOT) write(nud,*)                                    &
     &      'Restart predates the ecological stream: its interval starts clean'
      else
         call mpgetgp('aecotas'    ,aecotas    ,NHOR,1)
         call mpgetgp('aecots'     ,aecots     ,NHOR,1)
         call mpgetgp('aecops'     ,aecops     ,NHOR,1)
         call mpgetgp('aecohus'    ,aecohus    ,NHOR,1)
         call mpgetgp('aecowind'   ,aecowind   ,NHOR,1)
         call mpgetgp('aecoswd'    ,aecoswd    ,NHOR,1)
         call mpgetgp('aecoswu'    ,aecoswu    ,NHOR,1)
         call mpgetgp('aecoswn'    ,aecoswn    ,NHOR,1)
         call mpgetgp('aecolwn'    ,aecolwn    ,NHOR,1)
         call mpgetgp('aecolwu'    ,aecolwu    ,NHOR,1)
         call mpgetgp('aecoczen'   ,aecoczen   ,NHOR,1)
         call mpgetgp('aecopr'     ,aecopr     ,NHOR,1)
         call mpgetgp('aecoprsn'   ,aecoprsn   ,NHOR,1)
         call mpgetgp('aecoprc'    ,aecoprc    ,NHOR,1)
         call mpgetgp('aecoevap'   ,aecoevap   ,NHOR,1)
         call mpgetgp('aecotasmx'  ,aecotasmx  ,NHOR,1)
         call mpgetgp('aecotasmn'  ,aecotasmn  ,NHOR,1)
         call mpgetgp('aecotsmx'   ,aecotsmx   ,NHOR,1)
         call mpgetgp('aecotsmn'   ,aecotsmn   ,NHOR,1)
      endif
                                           
!     Spectral accumulators, NROOT only -- see the matching comment in epilog.
!     zaccuvers stays negative for a restart written before this patch: such a
!     file has no complete accumulator set, so the partial output interval it
!     carries is discarded and naccuout is zeroed to match. Without that the
!     counter would again outlive the accumulators it counts.
      zaccuvers = -1.0
      ienerok = 1
      if (mypid == NROOT) then
         nexcheck = 0
         call get_restart_real('accuvers',zaccuvers)
!        THE TWO ENERGY-DIAGNOSTIC ACCUMULATORS ARE ASKED FOR BY NAME, world-5qy.
!        `epilog` writes them only where they are allocated, so the marker cannot
!        promise them and their presence has to be a question. `has_restart_array`
!        reads nothing and moves no read position, so asking is free and safe at
!        any nexcheck; the answer travels to the other tasks through mpbci below,
!        because yresnam is threadprivate and only NROOT has opened the file.
         if (nenergy > 0) then
            call has_restart_array('adenergy',ifound)
            if (ifound == 0) ienerok = 0
         endif
         if (nener3d > 0) then
            call has_restart_array('adener3d',ifound)
            if (ifound == 0) ienerok = 0
         endif
         if (zaccuvers >= 2.0) then
            call get_restart_array('aasosp' ,aaso   ,NESP,NESP,   1)
            call get_restart_array('aaspsp' ,aasp   ,NESP,NESP,   1)
            call get_restart_array('aastsp' ,aast   ,NESP,NESP,NLEV)
            call get_restart_array('aasqsp' ,aasqout,NESP,NESP,NLEV)
            call get_restart_array('aasdsp' ,aasd   ,NESP,NESP,NLEV)
            call get_restart_array('aaszsp' ,aasz   ,NESP,NESP,NLEV)
            call get_restart_real('aorbnu'  ,aorbnu)
            call get_restart_real('alambm'  ,alambm)
            call get_restart_real('azdecl'  ,azdecl)
            call get_restart_real('ardist'  ,ardist)
            call get_restart_real('arasc'   ,arasc )
         endif
!        The fixer's integrated correction, world-fsr. Read under the same
!        lowered nexcheck: a restart written before it was saved leaves
!        denergyfix at its declared zero, which is what that segment would have
!        started from anyway. Only NROOT ever reads or applies it.
         call get_restart_real('denergyfix',denergyfix)
         nexcheck = 1
      endif
      call mpbcr(zaccuvers)
      call mpbci(ienerok)
!     ONE DECISION, AND NOTHING IS READ BEFORE IT IS TAKEN. The aa* set below
!     used to be read AFTER the `outreset` a discarded interval calls, so at
!     nlowio > 0 -- the only setting at which those arrays hold anything, since
!     `outaccu` and `outreset` both gate them on it -- naccuout went to zero and
!     the partial sums were then loaded back over the reset. That is the counter
!     outliving what it counts, by the route `accuvers` was built to close.
!     Everything the output interval consists of is now inside one branch.
      if (zaccuvers < 2.0 .or. ienerok == 0) then
         call outreset          ! no complete accumulator set: start clean
         if (mypid == NROOT) then
            if (zaccuvers < 2.0) then
               write(nud,*)                                                 &
     &         'Restart predates the accumulator set: partial output interval', &
     &         ' discarded and naccuout reset to 0'
            else
               write(nud,*)                                                 &
     &         'Restart predates the energy-diagnostic accumulators and this', &
     &         ' run has them on: partial output interval discarded and',    &
     &         ' naccuout reset to 0'
            endif
         endif
      else
         call mpgetgp('afdsw1'       ,afdsw1      ,NHOR,1)
         call mpgetgp('afdsw2'       ,afdsw2      ,NHOR,1)
         call mpgetgp('agpi'         ,agpi        ,NHOR,1)
         call mpgetgp('aventi'       ,aventi      ,NHOR,1)
         call mpgetgp('alaav'        ,alaav       ,NHOR,1)
         call mpgetgp('ampoti'       ,ampoti      ,NHOR,1)
         call mpgetgp('avrmpi'       ,avrmpi      ,NHOR,1)
         call mpgetgp('acapen'       ,acapen      ,NHOR,1)
         call mpgetgp('alnb'         ,alnb        ,NHOR,1)
         call mpgetgp('achim'        ,achim       ,NHOR,1)
         call mpgetgp('aadq'        ,aadq    ,NHOR,NLEP)
         call mpgetgp('aammr'       ,aammr   ,NHOR,NLEP)
         call mpgetgp('aanrho'      ,aanrho  ,NHOR,NLEP)
         call mpgetgp('aadmld'      ,aadmld  ,NHOR,1)
         call mpgetgp('aadt'        ,aadt    ,NHOR,NLEP)
         call mpgetgp('aadwatc'     ,aadwatc ,NHOR,1)
         call mpgetgp('aadsnow'     ,aadsnow ,NHOR,1)
         call mpgetgp('aadql'       ,aadql   ,NHOR,NLEP)
         call mpgetgp('aadust3'     ,aadust3 ,NHOR,1)
         call mpgetgp('aadcc'       ,aadcc   ,NHOR,NLEP)
         call mpgetgp('aadtd5'      ,aadtd5  ,NHOR,1)
         call mpgetgp('aadls'       ,aadls   ,NHOR,1)
         call mpgetgp('aadz0'       ,aadz0   ,NHOR,1)
         call mpgetgp('aadalb'      ,aadalb  ,NHOR,1)
         call mpgetgp('aadsalb1'    ,aadsalb1,NHOR,1)
         call mpgetgp('aadsalb2'    ,aadsalb2,NHOR,1)
         call mpgetgp('aadtsoil'    ,aadtsoil,NHOR,1)
         call mpgetgp('aadtd2'      ,aadtd2  ,NHOR,1)
         call mpgetgp('aadtd3'      ,aadtd3  ,NHOR,1)
         call mpgetgp('aadtd4'      ,aadtd4  ,NHOR,1)
         call mpgetgp('aadicec'     ,aadicec ,NHOR,1)
         call mpgetgp('aadiced'     ,aadiced ,NHOR,1)
         call mpgetgp('aadforest'   ,aadforest   ,NHOR,1)
         call mpgetgp('aadwmax'     ,aadwmax     ,NHOR,1)
         call mpgetgp('aadglac'     ,aadglac     ,NHOR,1)
         call mpgetgp('aadqo3'      ,aadqo3      ,NHOR,NLEV)
         call mpgetgp('aagroundoro' ,aagroundoro ,NHOR,1)
         call mpgetgp('aaglacieroro',aaglacieroro,NHOR,1)
!        The energy diagnostics' own partial window, world-5qy. Present, because
!        the branch this is in was chosen by asking for them.
         if (nenergy > 0) call mpgetgp('adenergy',adenergy,NHOR,28)
         if (nener3d > 0) call mpgetgp('adener3d',adener3d,NHOR,NLEV*28)
      endif


      return
      end subroutine read_atmos_restart


!     ==========================
!     SUBROUTINE STABILITY_CHECK
!     ==========================

      subroutine stability_check
      use pumamod

!     Some operating systems ignore floating point exceptions
!     and continue to run the program after an explosion,
!     e.g. after wrong settings for timestep or other conditions.
!     This subroutine checks some variables for valid ranges
!     and aborts the program if it's obviously broken.

      if (mypid == NROOT) then
         if (gd(1,1) > 1000.0 .or. gd(1,1) < -1000.0 .or. &
             gz(1,1) > 1000.0 .or. gz(1,1) < -1000.0 .or. &
             (gt(1,1)*ct) > 1000.0 .or. gt(1,1) < -t0(1)) then
            open(44,file='Abort_Message')
            write(44,*) 'Planet Simulator aborted'
            write(44,*) 'timestep = ',nhcstp
            write(44,*) 'gd(1,1) = ',gd(1,1)
            write(44,*) 'gz(1,1) = ',gz(1,1)
            write(44,*) 'gt(1,1) = ',gt(1,1)
            close(44)
   
            write(nud,*) 'Planet Simulator aborted'
            write(nud,*) 'timestep = ',nhcstp
            write(nud,*) 'gd(1,1) = ',gd(1,1)
            write(nud,*) 'gz(1,1) = ',gz(1,1)
            write(nud,*) 'gt(1,1) = ',gt(1,1)
   
            stop
         endif
      endif
      end subroutine stability_check
            

!     =================
!     SUBROUTINE INITFD
!     =================

      subroutine initfd
      use pumamod

      if (nkits < 1) nkits = 1
!     ============================================================
!     next subroutine call to set inital temperature.
!     model started from rest with stratification given by
!     trs calculated in setzt(ex). a perturbation is added to
!     constant log(sp) by subroutine noise, called from setzt(ex).
!     ============================================================
!
       
      if (mypid == NROOT) then
      
       call setzt
      endif
      call mpscsp(sp,spm,1)
      if (abs(dttl)>0.0) then    !setzt created a stellar-antistellar dipole
        call mpscsp(sr,srm,NLEV)
        if (mypid == NROOT) then
           st(:,:) = sr(:,:)
          stm(:,:) =srm(:,:)
           sz(3,:) = plavor
          szm(3,:) = plavor
        endif
!         call mpscsp(sr,stm,NLEV)
      else
        if (mypid == NROOT) then
           st(1,:) = sr(1,:)
          stm(1,:) = sr(1,:)
           sz(3,:) = plavor
          szm(3,:) = plavor
        endif
      endif
      return
      end

!     =================
!     SUBROUTINE READNL
!     =================

      subroutine readnl
      use pumamod

      namelist /plasim_nl/ &
                     kick    , mpstep  , nadv    , naqua   , ncoeff     &
                   , ndel    , ndheat  , ndiag   , ndiagcf , ndiaggp    &
                   , ndiaggp2d , ndiaggp3d , ndesert                    &
                   , ndiagsp   , ndiagsp2d , ndiagsp3d, dttl            &
                   , ndl     , neqsig                                   &
                   , ngui    , nhdiff  , nhordif , nkits                &
                   , noutput , nlowio  , nstpw   , nsnapshot, nstps     &
                   , nperpetual        , nprhor                         &
                   , nprint  , nqspec  , nrad    , nsela   , nshtns     &
                   , nsync                                             &
                   , ntime   , ntspd   , nveg    , nwpd    &
                   , n_start_year , n_start_month, n_run_steps          &
                   , n_run_years , n_run_months  , n_run_days           &
                   , n_days_per_month, n_days_per_year, fixedlon        &
                   , nhcadence, hcstartstep, hcendstep, hcinterval      &
                   , neco    , necostep                            &
                   , seed    , nfilter , ngptfilter, nspvfilter          &
                   , landhoskn0, nfilterexp, filterkappa                &
                   , syncstr , synctime, frcmod                        &
                   , dtep    , dtns    , dtrop   , dttrp                &
                   , tdissd  , tdissz  , tdisst  , tdissq  , tgr        &
                   , psurf   , ptop    , ptop2   , taucool              &
                   , restim  , t0      , tfrc    , nstratosponge        &
                   , sigh    , nenergy , nener3d , nsponge , dampsp     &
                   , nenergyfix, nconvtime, ndealias                    &
                   , l_aero
!
!     preset namelist parameter according to model set up
!
!     A LAYER-COUNT DEFAULT, AND ITS POSITION IS WHAT MAKES IT ONE. `tfrc` is a
!     per-level array, so a Rayleigh drag profile is a statement about the top
!     two of however many layers there are and cannot be written without a layer
!     count; that is what separates this from the NTRU==42 preset world-677x
!     deleted, where nhdiff was an ABSOLUTE wavenumber whose meaning changed
!     under the truncation it was keyed to. This runs BEFORE read(11,plasim_nl),
!     so a caller that declares TFRC wins and a caller that declares none gets a
!     sponge rather than none. This project declares all ten levels from
!     model.rayleigh_sponge_rotations on every prepare and every continuation
!     (run_exoplasim.py:declare_dry_constants, world-aee), so the preset is
!     overwritten on every run here.
!
!     THE TWO NLEV==20 BLOCKS THAT STOOD AFTER THE READ ARE GONE, world-helo,
!     and `nrdrag` with them: the switch had nothing left to switch. They set
!     tfrc when nrdrag was 1, and a preset placed after the namelist read is not
!     a default -- it is a compiled constant wearing a namelist's clothes, and it
!     would have overwritten a declared TFRC with nothing on disk to show for it,
!     which verify_staged_namelists cannot catch because the namelist FILE still
!     holds the declared value. Deleting them is a no-op here (nrdrag's compiled
!     default was 0 and nothing in this tree set it, nor does any namelist on
!     disk name the key) and is NOT a no-op for a twenty-layer caller, which is
!     the same footing as the NSHALLOW deletion: that caller writes the profile
!     into TFRC in plasim_nl, which is where a per-level friction profile
!     belongs. The values are recorded in
!     exoplasim/notes/resolution-tuned-parameters.md section 3.
      if (NLEV==10) then
         tfrc(1)      =  20.0 * day_24hr * frcmod !day_24hr
         tfrc(2)      = 100.0 * day_24hr * frcmod !day_24hr
         tfrc(3:NLEV) =   0.0 * day_24hr * frcmod !day_24hr
      endif
!

!     NO TRUNCATION IS A SPECIAL CASE. Upstream set nhdiff, ndel and the four
!     tdiss* here for NTRU 42 alone, so one truncation carried its own damping
!     and every other one fell through to the module defaults at
!     plasimmod.f90:206, :447 and :918-921. The values were not the problem;
!     privileging a truncation was, because nhdiff is an ABSOLUTE wavenumber
!     and the same 16 confines a different fraction of the spectrum at each.
!
!     The T42 numbers survive in config/planet.yaml's model.hyperdiffusion,
!     recorded there as inherited from T42 and applied as a FRACTION of the
!     truncation so the confinement is the same at every rung.
!     run_exoplasim.py:declare_hyperdiffusion writes NHDIFF, NDEL and all four
!     TDISS* into the namelist read below, from the prepare path and the
!     continuation path alike, and refuses a rung its table does not name.
!     A caller that writes none of those keys gets the module defaults at every
!     truncation, in [days], which dayseccheck converts.
!     world-677x; exoplasim/notes/resolution-tuned-parameters.md.
!
!     read namelist
!

      open(11,file=plasim_namelist,form='formatted')
      read (11,plasim_nl)

!     NO DEFAULT IS SET BELOW THIS LINE. A preset that wants to be a default goes
!     ABOVE the read, where the NLEV==10 block is, so a declared value beats it.
!     What follows resolves declarations that CONTRADICT each other -- aqua and
!     desert asked for together, vegetation on an aqua planet -- where there is
!     no consistent state to give the caller both. That is a different thing from
!     a compiled value replacing one the caller wrote and could have had.
!     world-helo.

      if ((ndesert == 1) .and. (naqua == 1)) then !If both toggled, turn off both
        naqua = 0
        ndesert = 0
      endif


!     aqua planet settings
      if (naqua == 1) then
         nveg = 0 ! switch off vegetation
      endif
      
!     We won't turn off vegetation for the desert planet, figuring that interested parties
!     might appreciate having the vegetation naturally respond to low moisture.
      solar_day = day_24hr
      sidereal_day = solar_day * (n_days_per_year-1) / n_days_per_year

!     set rotation dependent variables

      if (rotspd /= 1.0) then
!          if (n_days_per_year == 365) n_days_per_year = 360
!          solar_day = solar_day / rotspd
!          solar_day = (n_days_per_year-1.0)*sidereal_day/n_days_per_year
         sidereal_day = day_24hr / rotspd
!          n_days_per_year = n_days_per_year * rotspd
!          sidereal_day =(n_days_per_year*day_24hr)/(n_days_per_year+1.0)
         if (n_days_per_year /= 1) then
           solar_day = sidereal_day * n_days_per_year/(n_days_per_year-1)
         else
           solar_day = sidereal_day !In this case the solar day is infinite, so we set it to 1 year
         endif
!          day_24hr = 86400.0 ! WHY DOESN'T THIS WORK???
!        If there's one day per year, then solar_day is zero. But really, it's infinite
      endif
      
!     THE CALENDAR'S THREE LENGTHS, and all three are set here because calmod
!     copies them rather than deriving its own. m_days_per_month had no setter
!     at all and stayed at Earth's 30 while the year became 183, so twelve
!     months did not span a year. The division ROUNDS UP, so twelve months
!     always cover the orbit and the last one is the short one; rounding down
!     gives a thirteenth month. The floors keep step2cal30's mod and divide off
!     zero for an orbit shorter than a day or a year shorter than twelve.
!     world-x1k.
      m_days_per_year = max(1,nint(n_days_per_year * sidereal_day / day_24hr)) !24-hour days per year
      m_days_per_month = max(1,(m_days_per_year + 11) / 12) !24-hour days per month
      n_days_per_month = m_days_per_month
      
      
!       day_24hr IS NOT PURELY A UNIT CONVERSION, and the upstream comment that
!       said so is what hid world-rt1. Two distinct roles meet in this symbol:
!
!         - the seconds in a 24-hour day, which is the unit tfrc, restim,
!           tdiss*, dampsp and taucool are ENTERED in and the unit
!           config/planet.yaml derives its timescales_days in. That role is
!           genuinely a conversion and day_24hr is right for it.
!         - the model's unit of time, which is 1/ww = sidereal_day/TWOPI. That
!           role belongs to sidereal_day, and using day_24hr for it applied
!           every damping timescale over 1/rotspd times its namelist value.
!
!       The two coincide on Earth, where rotspd is 1, and nothing in the source
!       assigns day_24hr anywhere: it holds 86400.0 for the whole run. The two
!       roles are separated by name now, so neither depends on the other.

      ww    = TWOPI / sidereal_day ! Omega (scaling)
      acpd  = gascon / akap        ! Specific heat for dry air
      adv   = ACPV / acpd -1.0     ! Often used
      cv    = plarad * ww          ! cv
      ct    = cv * cv / gascon     ! ct
      pnu21 = 1.0 - 2.0 * pnu      ! Time filter 2
      rdbrv = gascon / RV          ! rd / rv
!
!     calendar and time control
!     set simulation length by using the following parameters:
!     "n_run_years" and "n_run_months"

!     mpstep <= 0 and ntspd <= 0 triggers automatic
!
!     A DEFAULT, AND IT STANDS, world-helo. It fires only for a caller that
!     supplied neither MPSTEP nor NTSPD, so it never overrides a declared step;
!     this project writes MPSTEP from model.timestep_minutes on every prepare and
!     every continuation, so nothing here reaches it.
!
!     The two hand-picked arms are kept rather than folded into the formula,
!     which is the change that would look like tidying. At nlat 32 the formula
!     gives 60 minutes, and 60 is exactly the coarsest step T21 has been MEASURED
!     to start clean at (lib/rungs.py STABILITY_CEILING_MINUTES): folding the arm
!     in would hand a caller who declared no timestep the rung's measured ceiling
!     as its default, with no margin. 45 sits below it. Nothing here says what
!     the T31 arm's 36 is worth, because T31 is not a rung this project probes.
      if (mpstep <= 0 .and. ntspd <= 0) then
         if (nlat <= 32) then       ! T21
            mpstep = 45
         else if (nlat <= 48) then  ! T31
            mpstep = 36
         else                       ! T42 and more
            mpstep = (30 * 64) / nlat
         endif
      endif

!     Make sure that (mpstep * 60) * ntspd = day_24hr

      if (mpstep > 0) then             ! timestep given in [min]
         mtspd = nint(day_24hr) / nint(mpstep * 60)
         mtspd = mtspd + mod(mtspd,2)  ! make even
      endif
      
      mpstep = day_24hr  / real(mtspd * 60)
      ntspd = nint(solar_day) / nint(mpstep * 60)
      nafter = mtspd
      if (nwpd > 0 .and. nwpd <= mtspd) then
         nafter = mtspd / nwpd
      endif
      
      n_steps_per_year = nint(sidereal_year / (mpstep*60.0))
      
      if (nlowio > 0 .and. nstpw > 0) nafter = nstpw
      
      if (ndiag < 1) ndiag = 10 * mtspd

      if (nstps == 0) nstps = mtspd
      
!
!     for column runs set horizontal diffusion coefficients to 0
!
      if (nhordif==0) then
       restim = 0.0
       tfrc   = 0.0
       tdissd = 0.0
       tdissz = 0.0
       tdisst = 0.0
       tdissq = 0.0
      endif

      if (synctime > 0.0) syncstr = 1.0 / (TWOPI * synctime)

      write(nud,'(/,"****************************************")')
      write(nud,'("* plasim_nl from    <",a16,"> *")') plasim_namelist
      write(nud,'("****************************************")')
      write(nud,plasim_nl)

!     Convert start date to timesteps since 1-Jan-0000

      call calini(n_days_per_month,n_days_per_year,n_start_step,ntspd &
                   ,solar_day,0,mpstep,mcal_days_per_year &
                   ,m_days_per_year,m_days_per_month,mtspd)
!                  ,day_24hr,0)
      
      call cal2step(n_start_step,mtspd,n_start_year,n_start_month,1,0,0)

!     Compute simulation time in [months]

      n_run_months = n_run_years * 12 + n_run_months
      iyea = n_run_months / 12
      imon = mod(n_run_months,12)

!     Print some values

      write(nud,'(/,"*************************************")')
      write(nud,'("* Solar    day      :",f10.1," [s] *")') solar_day
      write(nud,'("* Sidereal day      :",f10.1," [s] *")') sidereal_day
      write(nud,'("* Omega             :",f8.2," [s-6] *")') ww * 1.0e6
      write(nud,'("* Rotation Speed    :",f10.8,"     *")') rotspd
      write(nud,'("* Days / Year       :",i8,"       *")') n_days_per_year
      write(nud,'("* Days / Month      :",i8,"       *")') n_days_per_month
      write(nud,'("* Timestep          :",f8.3," [min] *")') mpstep
      write(nud,'("* Timesteps / write :",i8,"       *")') nafter
      write(nud,'("* Timesteps / day   :",i8,"       *")') ntspd
      if (iyea  > 1 .and. imon == 0) then
         write(nud,'("* Simulation time:  ",i7,"  years *")') iyea
      else if (iyea == 1 .and. imon == 0) then
         write(nud,'("* Simulation time:     one year     *")')
      else if (n_run_months > 1) then
         write(nud,'("* Simulation time:  ",i7," months *")') n_run_months
      else if (n_run_months == 1) then
         write(nud,'("* Simulation time:    one month     *")')
      else if (n_run_days  > 1) then
         write(nud,'("* Simulation time:  ",i7,"   days *")') n_run_days
      else if (n_run_days == 1) then
         write(nud,'("* Simulation time:      one day     *")')
      else if (n_run_steps  > 1) then
         write(nud,'("* Simulation time:  ",i7,"  steps *")') n_run_steps
      else if (n_run_steps == 1) then
         write(nud,'("* Simulation time:   single step    *")')
      endif
      write(nud,'("*************************************")')

!     set sponge layer time scale
!     The [days]-to-[sec] guard keeps day_24hr, which is the unit the value is
!     ENTERED in; the nondimensionalisation takes sidereal_day, which is the
!     unit the model INTEGRATES in. world-rt1.

      if(dampsp > 0.) then
       if(dampsp < (day_24hr/mtspd)) dampsp=dampsp*day_24hr
       dampsp=sidereal_day/(TWOPI*dampsp)
      endif

!     set franks diagnostics

      if(ndiaggp==1) then
       ndiaggp3d=21+ndiaggp3d
      end if
      if(ndiagsp==1) then
       ndiagsp3d=3+ndiagsp3d
      end if

      return
      end

      subroutine dayseccheck(pf,yn)
      use pumamod
      real :: pf(NLEV)
      character (len=*) :: yn

!     THE WHOLE ARRAY HAS TO BE IN ONE UNIT. world-720.
!
!     This routine decides [days] against [sec] from maxval alone and then
!     converts every level. A namelist that sets element 1 only -- which is what
!     a Fortran scalar assignment to an array key does -- leaves levels 2..NLEV
!     at their compiled defaults, so a scalar written in one unit beside
!     defaults in the other gives a mixed-unit array. maxval then picks a single
!     unit for the whole of it, one level is converted the wrong way, and that
!     level carries orders of magnitude too much or too little damping,
!     silently. Write the key as NLEV*value and the question does not arise.
!
!     A right answer rather than a comparison: an array whose positive entries
!     straddle the discriminator cannot be in one unit, whatever the units are.

      zmax = maxval(pf(:))
      zmin = day_24hr
      do jlev = 1 , NLEV
         if (pf(jlev) > 0.0) zmin = min(zmin,pf(jlev))
      enddo
      if (zmax >= (day_24hr / mtspd) .and. zmin < (day_24hr / mtspd)      &
     &    .and. zmax > 0.0) then
         write(nud,*) 'MIXED UNITS in ',trim(yn),': min ',zmin,' max ',zmax
         write(nud,*) 'the timestep is ',day_24hr/mtspd,' [sec], so these'
         write(nud,*) 'cannot all be [days] or all be [sec]. A namelist key'
         write(nud,*) 'written as a scalar sets element 1 only; write it as'
         write(nud,*) 'NLEV*value instead. see world-720'
         stop 'mixed units in a per-level timescale'
      endif
      if (zmax < (day_24hr / mtspd) .and. zmax > 0.0) then
         write(nud,*) 'old maxval(',trim(yn),') = ',zmax
         write(nud,*) 'assuming [days] - converting to [sec]'
         pf(:) = pf(:) * day_24hr
         write(nud,*) 'new maxval(',trim(yn),') = ',maxval(pf(:))
      endif
      return
      end

!     =================
!     SUBROUTINE INITPM
!     =================

      subroutine initpm
      use pumamod

      real (kind=8) radea,zakk

!     *************************************************************
!     * carries out all initialisation of model prior to running. *
!     * major sections identified with comments.                  *
!     * this s/r sets the model parameters and all resolution     *
!     * dependent quantities.                                     *
!     *************************************************************

      radea = plarad

!     *********************
!     * set vertical grid *
!     *********************

      if(neqsig==-1) then
       sigmah(:)=sigh(:)
       
      elseif(neqsig==1) then
       do jlev = 1 , NLEV
        sigmah(jlev) = real(jlev) / NLEV
       enddo
       
      elseif(neqsig==2) then
      
       zsk0 = log10(ptop/psurf) + 1
       zskf = 0.99
       dzsk = (zskf - zsk0) / REAL(NLEV-1)
       do jlev = 2 , NLEV
        zsk = zsk0 + (jlev-1)*dzsk - 1
        sigma(jlev) = 10**zsk
       enddo
       sigma(1) = 0.5*sigma(2)
       sigmah(NLEV) = 1.0
       sigmah(1:NLEV-1) = 0.5*(sigma(1:NLEV-1)+sigma(2:NLEV))
       
      elseif(neqsig==3) then
      
       zsk0 = log10(ptop/psurf) + 1
       zskf = 0.99
       dzsk = (zskf - zsk0) / REAL(NLEV-1)
       do jlev = 2 , NLEV
        zfk = real(jlev-1)/real(NLEV-1)
        zffk = (real(jlev-1)/real(NLEV))**0.25
        zsk = zsk0 + (jlev-1)*dzsk - 1
        sigma(jlev) = (1.0-exp(-((1-zfk)**2)/0.05))*10**zsk + exp(-((1-zfk)**2)/0.05)*zffk
       enddo
       sigma(1) = 0.5*sigma(2)
       sigmah(NLEV) = 1.0
       sigmah(1:NLEV-1) = 0.5*(sigma(1:NLEV-1)+sigma(2:NLEV))
       
      elseif(neqsig==4) then
      
!      THE QUARTIC, and what fixes its coefficients. sigmah(zsk) with
!      zsk = jlev/NLEV is a quartic with no quadratic term, so three free
!      coefficients, and two of them are determined rather than fitted:
!
!        sigmah(1) = 1                 the last half-level IS the surface
!        d sigmah/d zsk = 0 at zsk = 1 half-levels bunch at the surface, which
!                                      is what puts resolution in the boundary
!                                      layer
!
!      a1 + a3 + a4 = 1 and a1 + 3 a3 + 4 a4 = 0, so choosing a1 leaves no
!      freedom: a4 = -(1 - 4 a1)/... solves to a3 = 1.75, a4 = -1.5 at
!      a1 = 0.75. The ONE fitted number is a1, the slope at the model top, and
!      it sets how much of the column the upper half spans. It is upstream's and
!      carries no derivation there.
!
!      The three lines below are the whole of what neqsig==4 adds over the
!      neqsig==0 fallback: shift the top half-level to zero, normalise, and map
!      affinely onto [ptop/psurf, 1]. Without them the model top sits wherever
!      the polynomial puts it at this NLEV -- sigma 0.0766, or 7660 Pa at ten
!      layers -- and `ptop` is ignored. world-a05.
       do jlev=1,NLEV
        zsk=REAL(jlev)/REAL(NLEV)
        sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       zsk = ptop/psurf
       sigmah(:) = sigmah(:) - sigmah(1)
       sigmah(:) = sigmah(:)/sigmah(NLEV)
       sigmah(:) = sigmah(:)*(1.0-zsk) + zsk
       
      elseif ((neqsig==5) .and. (NLEV .gt. 10)) then
       
       !Bottom atmosphere (PlaSim's normal domain)
       do jlev=max(NLEV-10,1),NLEV  !max() so the bound is provably >= 1; the
         !branch already requires NLEV > 10, so this never changes the range
         zsk=REAL(jlev-(NLEV-9+1))/10.0  !As if it was a 10-layer atmosphere
         sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       zsk = ptop/psurf
       sigmah(NLEV-9:) = sigmah(NLEV-9:) - sigmah(NLEV-9)
       sigmah(NLEV-9:) = sigmah(NLEV-9:)/sigmah(NLEV)
       sigmah(NLEV-9:) = sigmah(NLEV-9:)*(1.0-zsk) + zsk
       
       !Upper atmosphere (stratosphere)
       zsk0 = log10(ptop2/psurf) + 1 !Top at 1 hPa
       zskf = sigmah(max(NLEV-10,1)) !Have to add this max statement so it'll compile with debug flag
       dzsk = (zskf - zsk0) / REAL(NLEV-8)
       do jlev = 2 , NLEV-8
          zsk = zsk0 + (jlev-1)*dzsk - 1
          sigma(jlev) = 10**zsk
       enddo
       sigma(NLEV-8) = 0.5*(sigmah(NLEV-9)+sigmah(NLEV-8))
       sigma(1) = 0.5*sigma(2)
       sigmah(1:NLEV-9) = 0.5*(sigma(1:NLEV-9)+sigma(2:NLEV-8))
         
      else
       do jlev=1,NLEV
        zsk=REAL(jlev)/REAL(NLEV)
        sigmah(jlev)=0.75*zsk+1.75*zsk**3-1.5*zsk**4
       enddo
       
      end if

      dsigma(1     ) = sigmah(1)
      dsigma(2:NLEV) = sigmah(2:NLEV) - sigmah(1:NLEV-1)

      rdsig = 0.5 / dsigma
 
      sigma(1     ) = 0.5 * sigmah(1)
      sigma(2:NLEV) = 0.5 * (sigmah(1:NLEV-1) + sigmah(2:NLEV))
        
!     DIMENSIONLESS DAMPING RATES. THE UNIT OF TIME IS 1/ww, NOT A 24-HOUR DAY.
!     world-rt1.
!
!     `ww = TWOPI/sidereal_day` is the model's unit of frequency, so a
!     dimensional timescale T in seconds becomes the rate `1/(T*ww)` =
!     `sidereal_day/(TWOPI*T)`. These lines divided by `day_24hr` instead. The
!     ratio is `day_24hr/sidereal_day` = `rotspd`, which is 1 on Earth and
!     hides there; off Earth every one of restim, tfrc, tdiss* and dampsp was
!     applied over 1/rotspd times its namelist value.
!
!     TWO DIFFERENT DAYS MEET HERE AND THEY ARE NOT INTERCHANGEABLE. The unit a
!     namelist value is ENTERED in stays the 24-hour day, because that is what
!     `config/planet.yaml`'s `hyperdiffusion.timescales_days` is derived in and
!     what `scripts/check_consistency.py` recomputes it in; `dayseccheck` below
!     therefore still multiplies by `day_24hr`. The unit the model INTEGRATES in
!     is the sidereal day. Changing either one changes what the model does, so
!     neither is a free conversion factor.
!
!     dayseccheck assumes units [days] if values < timestep
!     and converts values to [sec] (compatibilty routine)

      call dayseccheck(restim,"restim")
      call dayseccheck(tfrc  ,"tfrc"  )
      call dayseccheck(tdissd,"tdissd")
      call dayseccheck(tdissz,"tdissz")
      call dayseccheck(tdisst,"tdisst")
      call dayseccheck(tdissq,"tdissq")

!     THE DIVISOR IS FLOORED BECAUSE THE MASK DOES NOT PROTECT IT. A `where`
!     selects which lanes the ASSIGNMENT stores and leaves the compiler free to
!     evaluate the right-hand side on every lane, and the declared
!     -ffpe-trap=zero turns a discarded lane into SIGFPE. restim and tfrc are
!     both `= 0.0` in their plasimmod declarations, so on a run that sets
!     neither the divisor is zero on EVERY lane and the only reason this has
!     not fired is that -O2 has not vectorised initpm. world-d016, and the same
!     mechanism as world-bhs.
!
!     THE FLOOR IS A NO-OP ON THE LANES THE MASK KEEPS. dayseccheck above has
!     already converted both arrays to SECONDS, so a lane the mask keeps holds
!     a relaxation time of order the timestep or longer; 1.0e-30 s is twenty-five
!     orders of magnitude below the shortest timestep this model can take. It is
!     also large enough that the quotient stays finite, which a floor at
!     tiny() would not be: -ffpe-trap=overflow is declared beside the zero trap.
      where (restim > 0.0)
         damp = sidereal_day / (TWOPI * max(restim,1.0e-30))
      elsewhere
         damp = 0.0
      endwhere

      where (tfrc > 0.0)
          tfrc = sidereal_day / (TWOPI * max(tfrc,1.0e-30))
      elsewhere
          tfrc = 0.0
      endwhere

!     compute internal diffusion parameter (LAUERSON)

      do jlev=1,NLEV
       jdel = ndel(jlev)
       if (tdissd(jlev) > 0.0) then
        tdissd(jlev) = sidereal_day/(TWOPI*tdissd(jlev))
       else
        tdissd(jlev)=0.
       endif
       if (tdissz(jlev) > 0.0) then
        tdissz(jlev) = sidereal_day/(TWOPI*tdissz(jlev))
       else
        tdissz(jlev)=0.
       endif
       if (tdisst(jlev) > 0.0) then
        tdisst(jlev) = sidereal_day/(TWOPI*tdisst(jlev))
       else
        tdisst(jlev) = 0.
       endif
       if (tdissq(jlev) > 0.0) then
        tdissq(jlev) = sidereal_day/(TWOPI*tdissq(jlev))
       else
        tdissq(jlev)=0.
       endif
       zakk=1./(real(NTRU-nhdiff)**jdel)
       jr=-1
       do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            zsq = (jn - nhdiff)
            if(jn >= nhdiff) then
             sak(jr,jlev) = zakk*zsq**jdel
            else
             sak(jr,jlev) = 0.
            endif
            sak(ji,jlev) = sak(jr,jlev)
         enddo
       enddo
      enddo

!     set coefficients which depend on wavenumber

      zrsq2 = 1.0 / sqrt(2.0)

      jr=-1
      do jm=0,NTRU
         do jn=jm,NTRU
            jr=jr+2
            ji=jr+1
            nindex(jr)=jn
            nindex(ji)=jn
            spnorm(jr)=zrsq2
            spnorm(ji)=zrsq2
         enddo
         zrsq2=-zrsq2
      enddo

! finally make temperatures dimensionless

      dtns  = dtns  / ct
      dtep  = dtep  / ct
      dttrp = dttrp / ct
      t0    = t0    / ct ! (NLEV)

!     print out

      zakk=tdisst(NLEV)/(real(NTRU-nhdiff)**ndel(NLEV))
      write(nud,'(/," *************************************************")')
      if (zakk == 0.0) then
      write (nud,'(" * No lateral dissipation *")')
      else
      write(nud,'(" * Lateral dissipation",13x,"NDEL(",i3,") =",i2," *")')&
         NLEV,ndel(NLEV)
      write(nud,'(" * Diffusion coefficient = ",g14.4," [m**",i1,"] *")')&
         zakk*ww*radea**ndel(NLEV),ndel(NLEV)
      write(nud,'(" * e-folding time for smallest scale =",f5.1," days *")')&
         1.0/(TWOPI*tdisst(NLEV))
      endif
      write(nud,'(" *************************************************")')
      return
      end


!     =================
!     SUBROUTINE MAKEBM
!     =================

      subroutine makebm
      use pumamod

      zdeltsq = delt * delt

      do jlev1 = 1 , NLEV
         do jlev2 = 1 , NLEV
            zaq = zdeltsq * (t0(jlev1) * dsigma(jlev2)                  &
     &          + dot_product(g(:,jlev1),tau(jlev2,:)))
            bm1(jlev2,jlev1,:) = zaq
         enddo
      enddo

      do jn=1,NTRU
         do jlev = 1 , NLEV
            bm1(jlev,jlev,jn) = bm1(jlev,jlev,jn) + 1.0 / (jn*(jn+1))
         enddo
         call minvers(bm1(1,1,jn),NLEV)
      enddo
      return
      end

!     =================
!     SUBROUTINE INITSI
!     =================

      subroutine initsi
      use pumamod
!===========================================================
! carries out all initialisation of model prior to running.
! major sections identified with comments.
! this s/r sets the variables and arrays associated with
! the semi-implicit scheme.
!===========================================================
!
      dimension zalp(NLEV),zh(NLEV)
!
! this value, used in setting alpha(1), is irrelevant in the
! angular momentum conserving ecmwf scheme

      tkp = akap * t0
      t01s2(1:NLEM) = t0(2:NLEV) - t0(1:NLEM)
      t01s2(  NLEV) = 0.0

      zalp(2:NLEV) = log(sigmah(2:NLEV)) - log(sigmah(1:NLEM))

      g      = 0.0
      g(1,1) = 1.0
      do jlev = 2 , NLEV
         g(jlev,jlev) = 1.0 - zalp(jlev)*sigmah(jlev-1)/dsigma(jlev)
         g(jlev,1:jlev-1) = zalp(jlev)
      enddo

      do jlev = 1 , NLEV
         c(jlev,:) = g(:,jlev) * (dsigma(jlev) / dsigma(:))
      enddo

      zt01s2   = t01s2(1)
      zsig     = sigmah(1)
      tau(1,1) = 0.5 * zt01s2 * (zsig - 1.0) + tkp(1)
      tau(2:NLEV,1) = 0.5 * zt01s2 * dsigma(2:NLEV)

      do 1410 jlev=2,NLEV
        zttm=zt01s2
        zsigm=zsig
        zt01s2=t01s2(jlev)
        zsig=sigmah(jlev)
        do 1420 jlev2=1,NLEV
          ztm=0.
          ztmm=0.
          if(jlev2.le.jlev) ztm=1
          if(jlev2.lt.jlev) ztmm=1
          ztau=zttm*(zsigm-ztmm)
          if(jlev.lt.NLEV) ztau=ztau+zt01s2*(zsig-ztm)
          ztau=ztau*rdsig(jlev)*dsigma(jlev2)
          if(jlev2.le.jlev) ztau=ztau+tkp(jlev)*c(jlev2,jlev)
          tau(jlev2,jlev)=ztau
 1420   continue
 1410 continue
!
      zfctr=0.001*CV*CV/ga
      do 1500 jlev = 1 , NLEV
         zh(jlev) = dot_product(g(:,jlev),t0) * zfctr
 1500 continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      write(nud,9001)
      write(nud,9003)
      write(nud,9002)
      do jlev = 1 , NLEV
        write(nud,9004) jlev,sigma(jlev),t0(jlev),zh(jlev)
      enddo
      write(nud,9000)

      if (nprint > 0) then
         write(nud,9012)
         write(nud,9013) (jlev,jlev = 1 , 5)
         write(nud,9012)
         do jlev = 1 , NLEV
           write(nud,9014) jlev,(c(i,jlev),i=1,5)
         enddo
         write(nud,9012)
      endif
      return
 9000 format(1x,33('*'),/)
 9001 format(/,1x,33('*'))
 9002 format(1x,33('*'))
 9003 format(' * Lv *    Sigma Basic-T  Height *')
 9004 format(' *',i3,' * ',3f16.8,' *')
 9012 format(1x,69('*'))
 9013 format(' * Lv * C',i11,4i12,' *')
 9014 format(' *',i3,' * ',5f12.8,' *')
      end

!     ==================
!     SUBROUTINE MINVERS
!     ==================

      subroutine minvers(a,n)
      dimension a(n,n),b(n,n),indx(n)

      b = 0.0
      do j = 1 , n
         b(j,j) = 1.0
      enddo
      call ludcmp(a,n,indx)
      do j = 1 , n
         call lubksb(a,n,indx,b(1,j))
      enddo
      a = b
      return
      end

!     =================
!     SUBROUTINE LUBKSB
!     =================

      subroutine lubksb(a,n,indx,b)
      dimension a(n,n),b(n),indx(n)
      k = 0
      do i = 1 , n
         l    = indx(i)
         sum  = b(l)
         b(l) = b(i)
         if (k > 0) then
            do j = k , i-1
               sum = sum - a(i,j) * b(j)
            enddo
         else if (sum /= 0.0) then
            k = i
         endif
         b(i) = sum
      enddo

      do i = n , 1 , -1
         sum = b(i)
         do j = i+1 , n
            sum = sum - a(i,j) * b(j)
         enddo
         b(i) = sum / a(i,i)
      enddo
      return
      end

!     =================
!     SUBROUTINE LUDCMP
!     =================

      subroutine ludcmp(a,n,indx)
      dimension a(n,n),indx(n),vv(n)

      d = 1.0
      vv = 1.0 / maxval(abs(a),2)

      do 19 j = 1 , n
         do i = 2 , j-1
            a(i,j) = a(i,j) - dot_product(a(i,1:i-1),a(1:i-1,j))
         enddo
         aamax = 0.0
         do i = j , n
            if (j > 1)                                                  &
     &      a(i,j) = a(i,j) - dot_product(a(i,1:j-1),a(1:j-1,j))
            dum = vv(i) * abs(a(i,j))
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
         enddo
         if (j .ne. imax) then
            do 17 k = 1 , n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
   17       continue
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j) == 0.0) a(j,j) = tiny(a(j,j))
         if (j < n) a(j+1:n,j) = a(j+1:n,j) / a(j,j)
   19 continue
      return
      end

!     =================
!     SUBROUTINE INILAT
!     =================

      subroutine inilat
      use pumamod
      do jlat = 1 , NLAT
         csq(jlat)  = 1.0 - sid(jlat) * sid(jlat)
         rcs(jlat)  = 1.0 / sqrt(csq(jlat))
      enddo
      do jlat = 1 , NLAT/2
         ideg = nint(180.0/PI * asin(sid(jlat)))
         write(chlat(jlat),'(i2,a1)') ideg,'N'
         write(chlat(NLAT+1-jlat),'(i2,a1)') ideg,'S'
      enddo
      return
      end


!     =====================
!     SUBROUTINE INITRANDOM
!     =====================

      subroutine initrandom
      use pumamod
      integer :: i, clock

!     Set random number generator seed

      call random_seed(size=nseedlen)
      allocate(meed(nseedlen))

!     Take seed from namelist parameter 'SEED' ?

      if (seed(1) /= 0) then
         meed(:) = 0
         i = nseedlen
         if (i > 8) i = 8
         meed(1:i) = seed(1:i)
      else
         call system_clock(count=clock)
         meed(:) = clock + 37 * (/(i,i=1,nseedlen)/)
      endif
      call random_seed(put=meed)
      return
      end

!     ================
!     SUBROUTINE NOISE
!     ================

      subroutine noise
      use pumamod

!     if kick is set to 1 or 2
!     adds white noise perturbation to ln(surface pressure)
!     balanced initial state at t=0.
!     for kick=2, the white noise pertubation is
!     symmetric to the equator
!     eps sets magnitude of the noise

      itp1 = NTP1 ! Suppress compiler warnings for T1
      if (itp1 <= 2 .and. kick > 2) kick = 0 ! for T1
      zeps=1.e-4
      zscale=zeps/sqrt(2.0)

      write(nud,'(/," *****************************************")')
      if (kick == 1) then
         jsp1=2*NTP1+1
         do jsp=jsp1,NRSP
            call random_number(zrand)
            if (mrpid > 0) zrand = zrand + mrpid * 0.01
            sp(jsp)=sp(jsp)+zscale*(zrand-0.5)
         enddo
         write(nud,'(" *     White noise added (KICK = 1)      *")')
      elseif (kick == 2) then
         jr=2*NTP1-1
         do jm=1,NTRU
            do jn=jm,NTRU
               jr=jr+2
               ji=jr+1
               if (mod(jn+jm,2) == 0) then
                  call random_number(zrand)
                  if (mrpid > 0) zrand = zrand + mrpid * 0.01
                  sp(jr)=sp(jr)+zscale*(zrand-0.5)
                  sp(ji)=sp(ji)+zscale*(zrand-0.5)
               endif
            enddo
         enddo
         write(nud,'(" * Symmetric white noise added (KICK=2) *")')
      elseif (kick == 3) then
         sp(2*itp1+3) = zscale
         sp(2*itp1+4) = zscale * 0.5
         write(nud,'(" *  Mode sp(1,1) disturbed (KICK = 3)    *")')
      endif
      write(nud,'(" *****************************************")')
      return
      end

!     ================
!     SUBROUTINE SETZT
!     ================
      subroutine setzt
      use pumamod
!
      dimension ztrs(NLEV)
      dimension ztdprs(NUGP,NLEV)
      dimension zfac(NLEV)
      dimension zdpfac(NUGP)
      dimension srlev(NESP)
!
!*********************************************************************
!  this s/r sets up restoration temp field.
! the temperature at sigma = 1 is tgr, entered in kelvin.
! a lapse rate of ALR k/m is assumed under the tropopause and zero
! above. the actual profile tends to this away from then
! tropopause, with smooth interpolation depending on dttrp
! at the model tropopause.the height of
! the tropopause is given as dtrop m.
!*********************************************************************
!
      sr(:,:) = 0.0 ! NESP,NLEV
!
      zdttrp=dttrp*ct
!
      zsigprev=1.
      ztprev=tgr
      zzprev=0.
      do 1100 jlev=NLEV,1,-1
        zzp=zzprev+(gascon*ztprev/ga)*log(zsigprev/sigma(jlev))
        ztp=tgr-dtrop*ALR
        ztp=ztp+sqrt((.5*ALR*(zzp-dtrop))**2+zdttrp**2)
        ztp=ztp-.5*ALR*(zzp-dtrop)
        ztpm=.5*(ztprev+ztp)
        zzpp=zzprev+(gascon*ztpm/ga)*log(zsigprev/sigma(jlev))
        ztpp=tgr-dtrop*ALR
        ztpp=ztpp+sqrt((.5*ALR*(zzpp-dtrop))**2+zdttrp**2)
        ztpp=ztpp-.5*ALR*(zzpp-dtrop)
        ztrs(jlev)=ztpp
        zzprev=zzprev                                                   &
     &        +(.5*(ztpp+ztprev)*gascon/ga)*log(zsigprev/sigma(jlev))
        ztprev=ztpp
        zsigprev=sigma(jlev)
1100  continue
!
!     **********************************
!     * write out vertical information *
!     **********************************
!
      if (nprint > 0) then
      write(nud,9001)
      write(nud,9003)
      write(nud,9002)
      endif
!     
      do 1200 jlev = 1 , NLEV
         if (nprint > 0) write(nud,9004) jlev,sigma(jlev),ztrs(jlev)
         ztrs(jlev)=ztrs(jlev)/ct
 1200 continue
!     
      if (nprint > 0) write(nud,9002)
!
!******************************************************************
! loop to set array zfac - this controls temperature gradients as a
! function of sigma in tres. it is a sine wave from one at
! sigma = 1 to zero at stps (sigma at the tropopause) .
!******************************************************************
! first find sigma at dtrop
!
      zttrop=tgr-dtrop*ALR
      ztps=(zttrop/tgr)**(ga/(ALR*gascon))
!
! now the latitudinal variation in tres is set up ( this being in terms
! of a deviation from t0 which is usually constant with height)
!
      zsqrt2=sqrt(2.)
      zsqrt04=sqrt(0.4)
      zsqrt6=sqrt(6.)
      
      zsubst1 = - 1.0/3.463843223261725/ct * dttl * cos(fixedlon*PI/180.0)
      zsubst2 =   1.0/3.463843223261725/ct * dttl * sin(fixedlon*PI/180.0)
      ksas = NTP1+1
      ksas1 = 2*ksas-1
      ksas2 = 2*ksas
      
      
      do 2100 jlev = 1 , NLEV
        zfac(jlev)=sin(0.5*PI*(sigma(jlev)-ztps)/(1.-ztps))
        if (zfac(jlev).lt.0.0) zfac(jlev)=0.0
        sr(1,jlev)=zsqrt2*(ztrs(jlev)-t0(jlev))
        sr(3,jlev)=(1./zsqrt6)*dtns*zfac(jlev)
        sr(5,jlev)=-2./3.*zsqrt04*dtep*zfac(jlev)
        sr(ksas1,jlev) = zsubst1*zfac(jlev)
        sr(ksas2,jlev) = zsubst2*zfac(jlev)
 2100 continue
!
      call initrandom
      call printseed
      call noise
!
      return
 9001 format(/,1x,26('*'))
 9002 format(1x,26('*'))
 9003 format(' * Lv *    Sigma Restor-T *')
 9004 format(' *',i3,' * ',f8.3,f9.3,' *')
      end


!     ====================
!     SUBROUTINE PRINTSEED
!     ====================

      subroutine printseed
      use pumamod
      integer :: i

      write (nud,9020)
      write (nud,9010)
      do i = 1 , nseedlen
         write (nud,9000) i,meed(i)
      enddo
      write (nud,9010)
      write (nud,9020)
      return
 9000 format('* seed(',i1,') = ',i10,' *')
 9010 format('************************')
 9020 format(/)
      end


!     ===============
!     SUBROUTINE DIAG
!     ===============

      subroutine diag
      use pumamod
      if (mod(nstep,ndiag) == 0) then
         if (ncoeff .gt. 0) call prisp
         call xsect
      endif
      call energy
      return
      end

!     ================
!     SUBROUTINE PRISP
!     ================

      subroutine prisp
      use pumamod

      character(len=30) title

      scale = 100.0
      title = 'Vorticity [10-2]'
      do 100 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sz(1,jlev),jlev,title,scale)
  100 continue

      title = 'Divergence [10-2]'
      do 200 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(sd(1,jlev),jlev,title,scale)
  200 continue

      scale = 1000.0
      title = 'Temperature [10-3]'
      do 300 jlev = 1 , NLEV
         if (ndl(jlev).ne.0) call wrspam(st(1,jlev),jlev,title,scale)
  300 continue

      if (nqspec == 1) then
         scale = 1000.0
         title = 'Specific Humidity [10-3]'
         do jlev = 1 , NLEV
            if (ndl(jlev).ne.0) call wrspam(sq(1,jlev),jlev,title,scale)
         enddo
      endif

      title = 'Pressure [10-3]'
      call wrspam(sp,0,title,scale)

      return
      end

!     =====================
!     SUBROUTINE GPAREAMEAN
!     =====================

      subroutine gpareamean(pf,pmean)
!
!     The AREA-weighted global mean of a distributed gridpoint field.
!
!     The one place in this model that turns a gridpoint field into a global
!     mean, because sum/NUGP is not that mean: the Gaussian latitudes are
!     unequally spaced, so the arithmetic mean over the cells over-weights the
!     poles, and by an amount that is a function of NLAT. Two rungs then report
!     different global means for the same simulated field, from the quadrature
!     alone. world-mt5.
!
!     gwd holds this rank's Gaussian weights, indexed over its own NLPP
!     latitudes, and sums to 2 over the globe. The normalisation is taken from
!     the summed weight rather than assumed, so a rank holding no latitudes
!     costs nothing.
!
!     COLLECTIVE. mpsumbcr carries an OpenMP barrier, so every rank has to reach
!     this; call it OUTSIDE any `mypid == NROOT` guard and print the result
!     inside one. That is the shape every caller here uses.
!
      use pumamod
      real, intent(in)  :: pf(NHOR)
      real, intent(out) :: pmean
      real :: zgw(NHOR)
      real :: zs(2)
      integer :: jlat, jlon, jhor

      jhor = 0
      do jlat = 1 , NLPP
       do jlon = 1 , NLON
        jhor = jhor + 1
        zgw(jhor) = gwd(jlat)
       enddo
      enddo
      zs(1) = dot_product(pf,zgw)
      zs(2) = sum(zgw)
      call mpsumbcr(zs,2)
      pmean = zs(1) / zs(2)
      return
      end

!     ================
!     FUNCTION UGPMEAN
!     ================

      function ugpmean(pf)
!
!     The AREA-weighted global mean of a GATHERED gridpoint field, pf(NUGP).
!
!     The companion to gpareamean above, for the diagnostic prints that already
!     hold the whole globe on one rank. It computes the Gaussian weights from
!     inigau rather than from gwd, because gwd is scattered and a gathered
!     field is indexed by GLOBAL latitude. That also makes it collective-free,
!     so it is safe inside a `mypid == NROOT` guard, which is where every one of
!     those prints lives. NLAT is small and these are startup prints, so
!     recomputing the nodes costs nothing worth avoiding.
!
      use pumamod
      real :: pf(NUGP)
      real (kind=8) :: zsi(NLAT), zgw(NLAT)
      integer :: jlat

      call inigau(NLAT,zsi,zgw)
      zs = 0.0
      zw = 0.0
      do jlat = 1 , NLAT
       zs = zs + zgw(jlat) * sum(pf((jlat-1)*NLON+1:jlat*NLON))
       zw = zw + zgw(jlat) * NLON
      enddo
      ugpmean = zs / zw
      return
      end

!     ==============
!     FUNCTION RMSSP
!     ==============

      function rmssp(pf)
      use pumamod
      real pf(NESP,NLEV)

      zsum = 0.0
      do jlev = 1 , NLEV
         zsum = zsum + dsigma(jlev)                                     &
     &        * (dot_product(pf(1:NZOM,jlev),pf(1:NZOM,jlev)) * 0.5     &
     &        +  dot_product(pf(NZOM+1:NRSP,jlev),pf(NZOM+1:NRSP,jlev)))
      enddo
      rmssp = zsum
      return
      end

!     =================
!     SUBROUTINE ENERGY
!     =================

      subroutine energy
      use pumamod

      parameter (idim=6) ! Number of scalars for GUI timeseries
      real (kind=4) ziso(idim)

      ziso(1) = umax   ! maximum value of csu (zonal mean cross section)
      ziso(2) = t2mean - TMELT  ! mean of 2m temperature
      ziso(3) = precip * 1.0e9  ! mean precipitation
      ziso(4) = evap   * 1.0e9  ! mean evaporation
      ziso(5) = olr             ! OLR
      ziso(6) = minval(dt(:,NLEP)) ! Minimum of surface temperature
!     ziso(6) = 1000.0 * sum(dq(:,NLEV))/2048.0 ! Mean of surface wetness

      call guiput("SCALAR" // char(0) ,ziso,idim,1,1)

      return
      end

!     ==================
!     SUBROUTINE UPDATIM
!     ==================

      subroutine updatim(kstep)
      use pumamod
       
      if (n_days_per_year == 365) then
         call step2cal(kstep,ntspd,ndatim)
      else
         call step2cal30(kstep,ndatim)
      endif
      return
      end

!     =================
!     SUBROUTINE WRSPAM
!     =================

      subroutine wrspam(ps,klev,title,scale)
      use pumamod
!
      dimension ps(NRSP)
      character(len=30) title
      character(len=18) datch


      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title,klev
      write(nud,20000)
      write(nud,20020) (i,i=0,9)
      write(nud,20000)
      write(nud,20100) (cab(i),i=1,10)
      write(nud,20200) (cab(i),i=NTRU+2,NTRU+10)
      write(nud,20300) (cab(i),i=2*NTRU+2,2*NTRU+9)
      write(nud,20400) (cab(i),i=3*NTRU+1,3*NTRU+7)
      write(nud,20000)
      write(nud,'(1x)')

      return

20000 format(1x,78('*'))
20020 format(' * n * ',10i7,' *')
20030 format(' *   * ',a18,2x,a30,'  Level ',i2,11x,'*')
20100 format(' * 0 *',f8.2,9f7.2,' *')
20200 format(' * 1 *',8x,9f7.2,' *')
20300 format(' * 2 *',15x,8f7.2,' *')
20400 format(' * 3 *',22x,7f7.2,' *')
      contains
      function cab(i)
      cab=real(scale*sqrt(ps(i+i-1)*ps(i+i-1)+ps(i+i)*ps(i+i)))
      end function cab
      end

!     ===============
!     SUBROUTINE WRZS
!     ===============

      subroutine wrzs(zs,title,scale)
      use pumamod
!
      dimension zs(NLAT,NLEV)
      character(len=30) title
      character(len=18) datch

      ip = NLAT / 16
      ia = ip/2
      ib = ia + 7 * ip
      id = NLAT + 1 - ia
      ic = id - 7 * ip

      call ntodat(nstep,datch)
      write(nud,'(1x)')
      write(nud,20000)
      write(nud,20030) datch,title
      write(nud,20000)
      write(nud,20020) (chlat(i),i=ia,ib,ip),(chlat(j),j=ic,id,ip)
      write(nud,20000)
      do 200 jlev = 1 , NLEV
         write(nud,20100) jlev,((int(zs(i,jlev)*scale)),i=ia,ib,ip),      &
     &                       ((int(zs(j,jlev)*scale)),j=ic,id,ip),jlev
  200 continue
      write(nud,20000)
      write(nud,'(1x)')

20000 format(1x,78('*'))
20020 format(' * Lv * ',16(1x,a3),' * Lv *')
20030 format(' *    * ',a18,2x,a30,20x,'*')
20100 format(' * ',i2,' * ',16i4,' * ',i2,' *')
      end

!     ================
!     SUBROUTINE WRORB
!     ================

      subroutine wrorb(zf,title)
      use pumamod
      
      real zf
      character(len=30) title
      character(len=18) datch
      
      call ntodat(nstep,datch)
      write(nud,20030) datch,title,zf
      
20030 format('>>>   * ',a18,2x,a30,' = ',f7.3,10x,'*') 
      end
      
!     ================
!     SUBROUTINE XSECT
!     ================

      subroutine xsect
      use pumamod
      use radmod
      character(len=30) title

      title = 'Fractional Day'
      call wrorb(zcdayf,title)
      title = 'True Anomaly [deg]'
      call wrorb(orbnu*180./PI,title)
      title = 'Solar Declination [deg]'
      call wrorb(zdeclf*180./PI,title)
      title = 'Ecliptic Longitude [deg]'
      call wrorb(lambm*180./PI,title)
      title = 'Right Ascension'
      call wrorb(rasc*180./PI,title)
      title = 'Distance Modulus'
      call wrorb(eccf,title)
      scale = 10.0
      title = 'Zonal Wind [0.1 m/s]'
      call wrzs(csu,title,scale)
      title = 'Meridional Wind [0.1 m/s]'
      call wrzs(csv,title,scale)
      scale = 1.0
      title = 'Temperature [C]'
      call wrzs(cst,title,scale)
      scale = 10000.0
      title = 'specific humidity [0.1g/Kg]'
      call wrzs(csm,title,scale)
      scale = 100.0
      title = 'cloud cover [%]'
      call wrzs(ccc,title,scale)
      return
      end


!     =====================
!     SUBROUTINE GRIDPOINTA
!     =====================

      subroutine gridpointa
      use pumamod
#ifdef OMPSHARED
      use shtnsmod, only: sh_sp2gp, sh_dv2uv, sh_sp2grad, shgdmu, shgdlam,    &
     &                    sh_gp2sp, sh_dztend, sh_advtend, sh_slice
#endif
!
!*    Adiabatic Gridpoint Calculations
!
!     gtn, gqn, gut, gvt, guz, gvz, gke, guq, gvq and gvpp live in pumamod now,
!     as bands of full-globe arrays, because SHTns analyses the globe in one
!     call. gphi and gpmt stay local: neither is transformed.
      real gphi(NHOR,NLEV)
      real gpmt(NLON,NLPP)
!     The tendency partials are zpsd, zpst, zpsz, zpsq and zpsp in pumamod,
!     one slot per process, written in place and reduced where they lie. They
!     were locals here and were copied into the reduction's buffer on the way
!     past; that copy was 68 GB over a 300-step T127 run.

      real zgp(NLON,NLAT)

!     The Gaussian weight spread over this rank's gridpoints, and the reduction
!     buffer for the area-weighted global means below. world-mt5.
      real zgw(NHOR)
      real zmean(5)

      real (kind=4) zcs(NLAT,NLEV)
      real (kind=4) zsp(NESP)

!
!     inverse Legendre transformation (spectral to fourier domain)
!     st -> gt  sq -> gq  sd -> gd  sz -> gz
!     (sd,sz) -> (gu,gv)
!     sp -> gp  sp -> gpj (dlnps/dphi)
!

#ifdef OMPSHARED
      if (nshtns == 1) then
!        A BARRIER BEFORE, AND IT IS NOT SYMMETRY WITH THE ONE AFTER. legmod
!        writes only the calling thread's band, so nothing it does can disturb
!        another thread; these wrappers write the WHOLE GLOBE of every field,
!        including bands other threads are still reading from earlier in the
!        timestep. Without this, a thread that arrives early overwrites gu for a
!        thread still in the physics behind it. The symptom is not a wrong
!        answer, it is a last-bit difference that appears in perhaps one run in
!        three -- which is why two runs agreeing is not evidence of anything.
!$omp barrier
!        SHTns does the Legendre transform and the FFT in ONE call, so the
!        fields land in GRID space here and the fc2gp block below is skipped.
!        Everything between the two has to be read with that in mind, which is
!        the whole reason this is a branch and not a call swap.
         call sh_dv2uv(sd, sz, gu_g, gv_g, NLEV)
         call sh_sp2gp(sd, gd_g, NLEV)
         call sh_sp2gp(st, gt_g, NLEV)
         call sh_sp2gp(sz, gz_g, NLEV)
         if (nqspec == 1) call sh_sp2gp(sq, gq_g, NLEV)
         call sh_sp2gp(sp, gp_g, 1)
         call sh_sp2grad(sp, shgdmu, shgdlam, 1)
!        The wrappers return without synchronising, on purpose: see the note in
!        shtnsmod. Nothing above this line may be read until here.
!$omp barrier
         gpj(:) = shgdmu(mypid*NHOR+1:mypid*NHOR+NHOR)
         gpmt   = reshape(shgdlam(mypid*NHOR+1:mypid*NHOR+NHOR), [NLON,NLPP])
      else
         call invlega
      endif
#else
      call invlega
#endif

      if (ngui > 0 .or. mod(nstep,ndiag) == 0) then
        do jlev = 1 , NLEV
          do jlat = 1 , NLPP
            sec = CV / sqrt(csq(jlat))
            j1=(jlat-1)*NLON+1
            j2=jlat*NLON
!           THE ZONAL MEAN, and where it comes from depends on the transform.
!           In Fourier space the first element of a latitude row IS the m=0
!           coefficient, which is that mean; in grid space it is one longitude
!           and means nothing. So under SHTns the mean is taken over the row,
!           the way ccc below has always taken it.
            if (nshtns == 1) then
               csu(jlat,jlev) = SUM(gu(j1:j2,jlev))/real(NLON) * sec
               csv(jlat,jlev) = SUM(gv(j1:j2,jlev))/real(NLON) * sec
               cst(jlat,jlev) =(SUM(gt(j1:j2,jlev))/real(NLON)               &
     &                          + t0(jlev))*ct-TMELT
            else
               csu(jlat,jlev) =  gu(1+(jlat-1)*NLON,jlev) * sec
               csv(jlat,jlev) =  gv(1+(jlat-1)*NLON,jlev) * sec
               cst(jlat,jlev) =(gt(1+(jlat-1)*NLON,jlev) + t0(jlev))*ct-TMELT
            endif
            ccc(jlat,jlev) = SUM(dcc(j1:j2,jlev))/real(NLON)
            if (nqspec == 1) then
               if (nshtns == 1) then
                  csm(jlat,jlev) = SUM(gq(j1:j2,jlev))/real(NLON)
               else
                  csm(jlat,jlev) = (gq(1+(jlat-1)*NLON,jlev))
               endif
            else
               csm(jlat,jlev) = sum(dq(j1:j2,jlev))/real(NLON)
            endif
          enddo
        enddo
        umax = maxval(csu)
      endif

!     The zonal derivative of ln(ps), built from gp's FOURIER coefficients by
!     multiplying each by i*m. It cannot survive a transform that lands in grid
!     space, and under SHTns it does not have to: sh_sp2grad returned it above,
!     from the same call that returned the meridional one.
      if (nshtns /= 1) then
         do jlat = 1 , NLPP
            do jlon = 1 , NLON , 2
              gpmt(jlon  ,jlat) = -gp(jlon+1+(jlat-1)*NLON) * ((jlon-1)/2)
              gpmt(jlon+1,jlat) =  gp(jlon  +(jlat-1)*NLON) * ((jlon-1)/2)
            enddo
         enddo

         call fc2gp(gu  ,NLON,NLPP*NLEV)
         call fc2gp(gv  ,NLON,NLPP*NLEV)
         call fc2gp(gt  ,NLON,NLPP*NLEV)
         call fc2gp(gd  ,NLON,NLPP*NLEV)
         call fc2gp(gz  ,NLON,NLPP*NLEV)
         call fc2gp(gpj ,NLON,NLPP)
         call fc2gp(gpmt,NLON,NLPP)
         call fc2gp(gp  ,NLON,NLPP)
         if (nqspec == 1) call fc2gp(gq  ,NLON,NLPP*NLEV)
      endif
      gp = exp(gp)



      call calcgp(gpmt,gphi)

      gut = gu * gt
      gvt = gv * gt
      gke = gu * gu + gv * gv
      if (nqspec == 1) then
         guq = gu * gq
         gvq = gv * gq
      endif
!
!     add non linear geopotential terms
!
      gke(:,:) = gke(:,:) + gphi(:,:)

!
!     fft
!

#ifdef OMPSHARED
      if (nshtns == 1) then
!        THE REDUCTIONS ARE GONE, not merely faster. legmod integrates a
!        thread's own latitudes and leaves a partial for mpsumscp to sum;
!        SHTns integrates the globe and returns the finished field, so each
!        wrapper writes sdt, stt, szt and spt directly. That is the structural
!        half of what the forward conversion buys.
!
!        The first barrier is for the grid arrays, which every thread has just
!        filled a band of; the second is for the spectral tendencies, which the
!        wrappers fill by level and the caller reads whole.
!        THE TENDENCIES ARE HELD AS A SLICE, not as a whole field: sdt is
!        (NSPP,NLEV), the thread's modes at every level, because that is what
!        the semi-implicit step in spectrala works on. The wrappers are parallel
!        over LEVELS and produce every mode of the levels they own. The two
!        decompositions are orthogonal, so the whole field lands in slot 0 of
!        the partial scratch -- already shared, already this shape, and idle on
!        this path -- and each thread then takes its slice. Giving the whole
!        field its own arrays would cost 11.8 MB a die at T170 against a 32 MB
!        target, to save a copy of 82 KB a thread.
!$omp barrier
         call sh_gp2sp(gvpp_g, zpsp(1,0), 1)
         call sh_dztend(gvz_g, guz_g, gke_g, zpsd(1,1,0), zpsz(1,1,0), NLEV)
         call sh_advtend(gtn_g, gut_g, gvt_g, zpst(1,1,0), NLEV)
         if (nqspec == 1)                                               &
     &      call sh_advtend(gqn_g, guq_g, gvq_g, zpsq(1,1,0), NLEV)
!$omp barrier
         call sh_slice(zpsp(1,0), spt, 1)
         call sh_slice(zpsd(1,1,0), sdt, NLEV)
         call sh_slice(zpsz(1,1,0), szt, NLEV)
         call sh_slice(zpst(1,1,0), stt, NLEV)
         if (nqspec == 1) call sh_slice(zpsq(1,1,0), sqt, NLEV)
      else
#endif
      call gp2fc(gtn ,NLON,NLPP*NLEV)
      call gp2fc(gqn ,NLON,NLPP*NLEV)

      call gp2fc(gut ,NLON,NLPP*NLEV)
      call gp2fc(gvt ,NLON,NLPP*NLEV)
      call gp2fc(guz ,NLON,NLPP*NLEV)
      call gp2fc(gvz ,NLON,NLPP*NLEV)
      call gp2fc(gke ,NLON,NLPP*NLEV)
      call gp2fc(gvpp,NLON,NLPP     )
      if (nqspec == 1) then
         call gp2fc(guq ,NLON,NLPP*NLEV)
         call gp2fc(gvq ,NLON,NLPP*NLEV)
      endif
!
!     direct Legendre transformation (fourier domain to spectral domain)
!
      call fc2sp(gvpp,zpsp(1,mypart))
      call mktend(zpsd(1,1,mypart),zpst(1,1,mypart),zpsz(1,1,mypart),   &
     &            gtn,gvz,guz,gke,gut,gvt)
      call mpsumscp(zpsp,spt,1)
      call mpsumscp(zpst,stt,NLEV)
      call mpsumscp(zpsd,sdt,NLEV)
      call mpsumscp(zpsz,szt,NLEV)
      if (nqspec == 1) then
         call qtend(zpsq(1,1,mypart),gqn,guq,gvq)
         call mpsumscp(zpsq,sqt,NLEV)
      endif
#ifdef OMPSHARED
      endif
#endif
!
!     compute entropy
!
!
!     save u, v, and ps (at time t) for tracer transport
!
       dp0(:)=psurf*gp(:)
       do jlev=1,NLEV
        du0(:,jlev)=cv*gu(:,jlev)*SQRT(rcsq(:))
        dv0(:,jlev)=cv*gv(:,jlev)*SQRT(rcsq(:))
       enddo
!
!     save u v t q and ps (at time t) for further diagnostics
!
      if(nenergy > 0) then
       dp(:)  =dp0(:)
       du(:,:)=du0(:,:)
       dv(:,:)=dv0(:,:)
       do jlev=1,NLEV
        dt(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
        if (nqspec == 1) dq(:,jlev)=gq(:,jlev)*psurf/dp(:)
       enddo
      endif

!     Diagnostic output and GUI calls

      if (ngui>0 .or. mod(nstep,ndiag)==0 .or. mod(nstep,nafter)==0) then
         call mpgagp(zgp,gp,1);
         call guips(zgp)                ! Send Ps for Hovmoeller
         call guihor("DQVI"//char(0),dqvi,1,1.0,0.0)! Vertically integrated q
         call guigv("GU"  // char(0),gu)            ! Send u to GUI
         call guigv("GV"  // char(0),gv)            ! Send v to GUI
!        AREA-WEIGHTED GLOBAL MEANS, world-mt5. These were sum/(NLON*NLAT),
!        which is the arithmetic mean of the cells and not the mean over the
!        sphere: the Gaussian latitudes are unequally spaced, so an unweighted
!        sum over-weights the poles, and by an amount that is a function of
!        NLAT. The same simulated climate then reports a different global mean
!        T2m, precipitation, evaporation and OLR at two rungs from the
!        quadrature alone, which is the one thing a resolution comparison must
!        not do.
!
!        gwd is this rank's Gaussian weights, indexed over its own NLPP
!        latitudes, and sums to 2 over the globe. The normalisation is taken
!        from the summed weight rather than assumed, the way the conversion
!        diagnostic below does, so it stays right if a rank holds no latitudes.
!        Summing locally and reducing also drops four full-globe gathers per
!        diagnostic step, which the gathered form needed and this does not.
         jhor = 0
         do jlat = 1 , NLPP
          do jlon = 1 , NLON
           jhor = jhor + 1
           zgw(jhor) = gwd(jlat)
          enddo
         enddo
         zmean(1) = dot_product(dtsa(:),zgw(:))
         zmean(2) = dot_product(dprc(:),zgw(:)) + dot_product(dprl(:),zgw(:))
         zmean(3) = dot_product(devap(:),zgw(:))
         zmean(4) = -dot_product(dftu(:,1),zgw(:)) ! make positive upwards
         zmean(5) = sum(zgw(:))
         call mpsumbcr(zmean,5)
         t2mean = zmean(1) / zmean(5)   ! Mean of 2m temperature
         precip = zmean(2) / zmean(5)   ! Convective plus large scale precip
         evap   = zmean(3) / zmean(5)   ! Evaporation
         olr    = zmean(4) / zmean(5)   ! OLR = dftu level 1
         gp(:) = gp(:) - 1.0
         call gp2fc(gp,NLON,NLPP)
         call fc2sp(gp,span)
         call mpsum(span,1)
 
         call mpgacs(csu)
         call mpgacs(csv)
         call mpgacs(cst)
         call mpgacs(csm)
         call mpgacs(ccc)
 
         if (mypid == NROOT) then
            zcs(:,:) = csu(:,:)
            call guiput("CSU"  // char(0) ,zcs ,NLAT, NLEV,1)
            zcs(:,:) = csv(:,:)
            call guiput("CSV"  // char(0) ,zcs ,NLAT, NLEV,1)
            zcs(:,:) = cst(:,:)
            call guiput("CST"  // char(0) ,zcs ,NLAT, NLEV,1)
            zsp(:) = span(:)
            call guiput("SPAN" // char(0) ,zsp ,NCSP,-NTP1,1)
         endif

         ! send fields to GUI for column visualization
         call guigvcol("GUCOL"  // char(0),gu,sellon) ! u
         call guigvcol("GVCOL"  // char(0),gv,sellon) ! v 
         call guigtcol(dt,sellon) ! t 
         call guid3dcol("DCCCOL" // char(0),dcc,sellon,NLEP,100.0,0.0) !cl-cov 
         call guid3dcol("DQCOL" // char(0),dq,sellon,NLEP,1000.0,0.0)  !dq
         call guid3dcol("DTDTCOL" // char(0),dtdt,sellon,NLEP,         &
                        day_24hr,0.0)            ! t-tendency
         call guid3dcol("DQDTCOL" // char(0),dqdt,sellon,NLEP,         &
                        1000.0*day_24hr,0.0)     ! q-tendency
      endif
      if (nqspec == 0) then
         do jlev = 1 , NLEV
           ! dq(:,jlev) = dq(:,jlev) + gqn(:,jlev) * psurf / dp(:)
         enddo
      endif
      return
      end

!     =================
!     SUBROUTINE CALCGP
!     =================

      subroutine dealias_gp(pgp,klev,premoved)
!     Project a gridpoint field onto the retained spectral modes and synthesise
!     it back, so what enters the products below is band limited at NTRU.
!     world-ly5, and it is Hoskins and Simmons (1975) section 2 option (ii).
!
!     WHY. Their sentence: "There are terms for which this grid is insufficient
!     for removing aliased interactions. These terms are the triple correlation
!     involved in the energy conversion term and in the vertical advection
!     terms. These require in theory M_g >= 4M + 1." This model runs
!     NLON = 3*NTRU + 1, which dealiases a product of TWO band-limited fields
!     and not of three. Truncating `zvgpg` before it multiplies anything makes
!     every product below quadratic in band-limited fields again, which the grid
!     does dealias -- which is the whole of the remedy and is why it goes here,
!     at the source, rather than on each term.
!
!     `premoved` is the fraction of the field's own norm the projection removes,
!     globally and mass-unweighted. It is the measurement that decides the
!     hypothesis on its own: if `zvgpg` carries nothing above NTRU then this is
!     a no-op and the aliasing cannot be what the sink is made of.
!
      use pumamod
      integer, intent(in)    :: klev
      real,    intent(inout) :: pgp(NHOR,klev)
      real,    intent(out)   :: premoved
      real, allocatable :: zfc(:,:), zpart(:,:), zslice(:,:), zfull(:,:)
      real :: zw(NHOR)
      real :: zs(2)

      allocate(zfc(NHOR,klev))
      allocate(zpart(NESP,klev))
      allocate(zslice(NSPP,klev))
      allocate(zfull(NESP,klev))

      zfc(:,:) = pgp(:,:)
      call gp2fc(zfc,NLON,NLPP*klev)
      do jlev = 1 , klev
         call fc2sp_t(zfc(1,jlev),zpart(1,jlev))
      enddo
      call mpsumsc(zpart,zslice,klev)
      call mpgallsp(zfull,zslice,klev)
      do jlev = 1 , klev
         call sp2fc_t(zfull(1,jlev),zfc(1,jlev))
      enddo
      call fc2gp(zfc,NLON,NLPP*klev)

      jhor = 0
      do jlat = 1 , NLPP
       do jlon = 1 , NLON
        jhor = jhor + 1
        zw(jhor) = gwd(jlat)
       enddo
      enddo
      zs(:) = 0.0
      do jlev = 1 , klev
       zs(1) = zs(1)                                                    &
     &       + dot_product((pgp(:,jlev)-zfc(:,jlev))**2,zw)
       zs(2) = zs(2) + dot_product(pgp(:,jlev)**2,zw)
      enddo
      call mpsumbcr(zs,2)
      premoved = 0.0
      if (zs(2) > 0.0) premoved = sqrt(zs(1)/zs(2))

      pgp(:,:) = zfc(:,:)

      deallocate(zfc)
      deallocate(zpart)
      deallocate(zslice)
      deallocate(zfull)
      return
      end

      subroutine calcgp(gpm,gphi)
!     gtn, gqn, guz, gvz and gvpp come from pumamod rather than the argument
!     list. They are bands of full-globe arrays, and a band is not a contiguous
!     (NHOR,NLEV) block, so passing one to an explicit-shape dummy would copy it
!     in and out on every call -- on both transform paths, not just the new one.

!     *****************************************************
!     * computes nonlinear tendencies in grid point space *
!     *****************************************************

      use pumamod

      real gphi(NHOR,NLEV)
      real gpm(NHOR)
      real zsdotp(NHOR,NLEM),zsumd(NHOR)
      real ztpta(NHOR),ztptb(NHOR)
      real zvgpg(NHOR,NLEV)
      real gtd(NHOR,NLEM)
      real gqm(NHOR,NLEM)                                                       !NEU

      real gud(NHOR,NLEM)
      real gvd(NHOR,NLEM)

      real ztv1(NHOR,NLEV),ztv2(NHOR,NLEV),zq(NHOR)
      real zdealr

      do jlev = 1 , NLEV
         zvgpg(:,jlev) = rcsq  * (gu(:,jlev) * gpm + gv(:,jlev) * gpj)
!
!     set pseudo temperatures to use virtual temperatures
!     (note: gq=ps*q)
!
         if (nqspec == 1) then
            zq(:)=AMAX1(gq(:,jlev)/gp(:),0.)
         else
            zq(:)=max(dq(:,jlev),0.0)
         endif
         ztv1(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               -t0(jlev)
         ztv2(:,jlev)=(gt(:,jlev)+t0(jlev))*(1.+(1./rdbrv-1.)*zq(:))    &
     &               /(1.+ADV*zq(:))                                    &
     &               -t0(jlev)

      enddo
!
!     THE DEALIASING TRUNCATION. world-ly5, default off. See dealias_gp above.
!     It goes here, after zvgpg is formed and before anything multiplies it, so
!     that the conversion, the vertical advection and the surface pressure
!     tendency all see the same band-limited field.
      if (ndealias > 0) then
         call dealias_gp(zvgpg,NLEV,zdealr)
         if (mypid == NROOT) then
            ddealias(1) = ddealias(1) + zdealr
            ddealias(2) = ddealias(2) + 1.0
            if (mod(nstep,ndiag) == 0 .and. ddealias(2) > 0.0) then
               write(nud,'(A,I9,2E15.6)') ' DEALIAS removed fraction ',  &
     &            nstep, zdealr, ddealias(1)/ddealias(2)
            endif
         endif
      endif

!     *******
!     * gvpp *
!     *******

      zsumd = dsigma(1) * gd(:,1)
      gvpp   = dsigma(1) * zvgpg(:,1)
      zsdotp(:,1) = zsumd + gvpp

      do jlev = 2 , NLEV-1
         zsumd = zsumd + dsigma(jlev) * gd(:,jlev)
         gvpp   = gvpp   + dsigma(jlev) * zvgpg(:,jlev)
         zsdotp(:,jlev) = zsumd + gvpp
      enddo

      zsumd = zsumd + dsigma(NLEV) * gd(:,NLEV)
      gvpp   = gvpp   + dsigma(NLEV) * zvgpg(:,NLEV)

!     **************
!     * loop  400: *
!     **************

      do jlev = 1 , NLEM
         zsdotp(:,jlev) = (sigmah(jlev) * (zsumd+gvpp) - zsdotp(:,jlev))
         gtd(:,jlev) = zsdotp(:,jlev) * (gt(:,jlev+1) - gt(:,jlev))
         gqm(:,jlev) = zsdotp(:,jlev) * (gq(:,jlev+1) + gq(:,jlev))
         gud(:,jlev) = zsdotp(:,jlev) * (gu(:,jlev+1) - gu(:,jlev))
         gvd(:,jlev) = zsdotp(:,jlev) * (gv(:,jlev+1) - gv(:,jlev))
      enddo

!     *************
!     * top level *
!     *************

      zsumd = zvgpg(:,1) * dsigma(1)

      gtn(:,1) = gt(:,1) * gd(:,1) - akap * ztv2(:,1) * gd(:,1)         &
     &   - rdsig(1)*(gtd(:,1) + t01s2(1) * (sigmah(1)*gvpp-zsumd))

      gqn(:,1) = - rdsig(1) * gqm(:,1)

      guz(:,1) =-gu(:,1) * gz(:,1) - gpj * ztv1(:,1) - rdsig(1)*gvd(:,1)
      gvz(:,1) = gv(:,1) * gz(:,1) - gpm * ztv1(:,1) - rdsig(1)*gud(:,1)

      dw(:,1)=gd(:,1)*gp(:)*psurf*ww

!     ****************
!     * inner levels *
!     ****************

      do jlev = 2 , NLEV-1
         ztpta = c(1,jlev) *  zvgpg(:,1)
         ztptb = c(1,jlev) * (zvgpg(:,1) + gd(:,1))

         do jlej = 2 , jlev
            ztpta = ztpta + c(jlej,jlev) *  zvgpg(:,jlej)
            ztptb = ztptb + c(jlej,jlev) * (zvgpg(:,jlej) + gd(:,jlej))
         enddo

         zsumd = zsumd + zvgpg(:,jlev) * dsigma(jlev)

         gtn(:,jlev) = gt(:,jlev) * gd(:,jlev)                          &
     &       + akap * ztv2(:,jlev) * (zvgpg(:,jlev) - ztptb)            &
     &       + tkp(jlev) * (zvgpg(:,jlev) - ztpta)                      &
     &       - rdsig(jlev) * (gtd(:,jlev) + gtd(:,jlev-1)               &
     &                       +gvpp*(t01s2(jlev)*sigmah(jlev)             &
     &                            +t01s2(jlev-1)*sigmah(jlev-1))        &
     &                       -zsumd*(t01s2(jlev-1)+t01s2(jlev))         &
     &                       +zvgpg(:,jlev)*dsigma(jlev)*t01s2(jlev-1))

         gqn(:,jlev) = - rdsig(jlev)*(gqm(:,jlev) - gqm(:,jlev-1))

         guz(:,jlev) = - gu(:,jlev) * gz(:,jlev) - gpj * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gvd(:,jlev) + gvd(:,jlev-1))

         gvz(:,jlev) =   gv(:,jlev) * gz(:,jlev) - gpm * ztv1(:,jlev)   &
     &       - rdsig(jlev)*(gud(:,jlev) + gud(:,jlev-1))

         dw(:,jlev)=(zvgpg(:,jlev)-ztptb)*gp(:)*psurf*ww

      enddo

!     ****************
!     * bottom level *
!     ****************

      ztpta = c(1,NLEV) *  zvgpg(:,1)
      ztptb = c(1,NLEV) * (zvgpg(:,1) + gd(:,1))

      do jlej = 2 , NLEV
         ztpta = ztpta + c(jlej,NLEV) *  zvgpg(:,jlej)
         ztptb = ztptb + c(jlej,NLEV) * (zvgpg(:,jlej) + gd(:,jlej))
      enddo

      gtn(:,NLEV) = gt(:,NLEV) * gd(:,NLEV)                             &
     &   + akap*ztv2(:,NLEV)*(zvgpg(:,NLEV)-ztptb)                      &
     &   + tkp(NLEV)*(zvgpg(:,NLEV)-ztpta)                              &
     &   - rdsig(NLEV)*(gtd(:,NLEM)                                     &
     &                 +t01s2(NLEV-1)*(sigmah(NLEV-1)*gvpp-zsumd))

      gqn(:,NLEV) = rdsig(NLEV) * gqm(:,NLEM)

      guz(:,NLEV) = -gu(:,NLEV)*gz(:,NLEV) - gpj*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gvd(:,NLEM)
      gvz(:,NLEV) =  gv(:,NLEV)*gz(:,NLEV) - gpm*ztv1(:,NLEV)           &
     &    - rdsig(NLEV) * gud(:,NLEM)

      dw(:,NLEV)=(zvgpg(:,NLEV)-ztptb)*gp(:)*psurf*ww

!
!     compute non linear geopotential terms
!

      ztv1(:,:)=ztv1(:,:)-gt(:,:)
      do jlev=1,NLEV
      do jhor=1,NHOR
       gphi(jhor,jlev)=dot_product(g(:,jlev),ztv1(jhor,:))              &
     &                *2./rcsq(jhor)
      enddo
      enddo

      return
      end

!     ====================
!     SUBROUTINE CONVWEIGHT
!     ====================

      subroutine convweight(psd,pgq,pgp,pval)
!     The reference conversion's divergence half on ONE divergence, mass
!     weighted into W/m2 by denergy02's own formula and reduced over the globe.
!     A CONTROL AND NOT A MODEL TERM. world-0ov, world-pkf.
!
!     `spectrala` prints this on the four divergences the adiabatic step has to
!     hand, and the sink is the gap between two of them. What it could not say
!     is WHICH operation opens the gap, because everything between the adiabatic
!     t+dt and the state at t+dt is inside `spectrald`. This is the same
!     quantity, callable, so `spectrald` can print it either side of each write
!     to `sdp` and name the operation instead of bracketing it.
!
      use pumamod
      real, intent(in)  :: psd(NESP,NLEV)
      real, intent(in)  :: pgq(NHOR,NLEV)
      real, intent(in)  :: pgp(NHOR)
      real, intent(out) :: pval
      real :: zwrk(NESP,NLEV)
      real :: zgp(NHOR,NLEV)
      real :: zw(NHOR)
      real :: zs(2)

      do jlev = 1 , NLEV
         zwrk(:,jlev) = 0.0
         do jlev2 = 1 , jlev
            zwrk(:,jlev) = zwrk(:,jlev)                                 &
     &                   - tkp(jlev) * c(jlev2,jlev) * psd(:,jlev2)
         enddo
      enddo
      zwrk(:,:) = zwrk(:,:) * ct * ww
      call sp2fl(zwrk,zgp,NLEV)
      call fc2gp(zgp,NLON,NLPP*NLEV)

      jhor = 0
      do jlat = 1 , NLPP
       do jlon = 1 , NLON
        jhor = jhor + 1
        zw(jhor) = gwd(jlat)
       enddo
      enddo

      zs(:) = 0.0
      do jlev = 1 , NLEV
       zs(1) = zs(1)                                                    &
     &  + dot_product(zgp(:,jlev)*acpd*(1.+adv*pgq(:,jlev))             &
     &                *pgp(:)/ga*dsigma(jlev),zw)
      enddo
      zs(2) = sum(zw)
      call mpsumbcr(zs,2)
      pval = zs(1) / zs(2)

      return
      end

!     ====================
!     SUBROUTINE SPECTRALA
!     ====================

      subroutine spectrala
      use pumamod
#ifdef OMPSHARED
      use shtnsmod, only: sh_sp2gp, sh_dv2uv
#endif
!
!*    Add adiabatic and diabatic tendencies
!
!     the adiabatic tendencies are added using the semi implicit scheme
!     described in
!     Hoskins and Simmons 1975 (Q.J.R.Meteorol.Soc.,101,637-655) (HS75)
!     To compare the code directly with HS75 the following notes might be
!     helpful (in addition to the comments below):
!
!     - *x* means model variable x
!     - eq.x referres to equation x of HS75
!     - the script T of HS75 (e.g. 1st term rhs of eq.9) is stored in *stt*
!     - the script D of HS75 (e.g. 1st term rhs of eq.8) is stored in *sdt*
!     - the script P of HS75 (e.g. 1st term rhs of eq.10) is stored in -*spt*
!     - the surface geopotential is stored in *so*
!
!     - the vertical scheme has being changed to the ECMWF scheme
!       (see e.g. Simmons and Burridge 1981, Mon.Wea.Rev.,109,758-766).
!       in this scheme,  matrix g differs from that in HS75.
!
!     - in addition to the dry HS75 version also the tendencys of specific
!       humidity are processed
!
      real azm(NSPP,NLEV)
      real atm(NSPP,NLEV)
      real aqm(NSPP,NLEV)
      real adt(NSPP,NLEV)
      real adm(NSPP,NLEV)
      real zgt(NSPP,NLEV)
      real zgm(NSPP,NLEV)
      real apm(NSPP)
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
      real, allocatable :: zsd(:,:),zsz(:,:),zsq(:,:),zsp(:),zst(:,:)
      real, allocatable :: zekin(:,:),zepot(:,:)
!
!     zttgp, ztgp, zqgp, zqmgp, zugp, zvgp, zpgp and zpmgp are in pumamod now,
!     as bands of full-globe arrays, for the reason gtn and hdu are: SHTns
!     writes the GLOBE in one call and cannot be handed a band. zekin and zepot
!     stay here because nothing transforms into them.
!
!     the energy fixer's scratch: the imbalance, the column heat capacity it is
!     spread over, and the weight sum that turns the pair into a per-area rate
      real :: zfix(3)
      real :: zfixr, zfixd, zfixc
      real :: zfixw(NHOR)
!
!     THE CONVERSION DECOMPOSITION, a CONTROL and not a model term (nenergy > 1).
!     world-0ov. The reference conversion is split across the semi-implicit
!     scheme: its advective half, tkp*(zvgpg-ztpta), is explicit in calcgp and
!     its divergence half is the tkp*c part of tau, applied implicitly in step 3
!     below. An enthalpy budget that reads calcgp's gtn and takes everything else
!     BY DIFFERENCE therefore books that second half to advection, which is what
!     the first attribution did. These arrays measure it directly, in the same
!     arithmetic and the same units as denergy02, so the two halves can be added.
!     zcnow holds the divergence at time t. sd is still the state gridpointa
!     read: what advances it to t+dt is step 4.b below, where `sdp = 2 sdt -
!     adm` writes THROUGH THE POINTER into this thread's slice of sd, and
!     mpsyncsp afterwards only publishes that write. So the same term can be
!     evaluated at the explicit half's time level and the semi-implicit
!     displacement read off as the difference. Every thread copies the WHOLE
!     array while owning one slice of it, which is why the copy has to be
!     closed against 4.b by a barrier. See the copy.
      real, allocatable :: zcnow(:,:), zcsdt(:,:), zcwrk(:,:), zcgp(:,:)
      real :: zcw(NHOR)
      real :: zcs(10)
!
!*    0. save prognostic variables at (t-dt)
!        and the non-linear divergence tendency terms
!
!     The shared spectral arrays are read in FULL by the phase before this one
!     and written by SLICE in this one, and nothing else separates the two.
!     Without this a thread arriving early overwrites what another is still
!     reading. Inert without -fopenmp.
!$omp barrier

      apm(:)   = spm(:)   ! log surface pressure
      adm(:,:) = sdm(:,:) ! divergence
      azm(:,:) = szm(:,:) ! (absolut) vorticity
      atm(:,:) = stm(:,:) ! temperature
      adt(:,:) = sdt(:,:) ! divergence tendency
      if (nqspec == 1) aqm(:,:) = sqm(:,:) ! spec.  humidity
!
!     The control's copy of the divergence at time t, taken BEFORE the solve
!     overwrites sdt and BEFORE step 4.b advances sd. See the declaration.
!
!     Every thread copies the FULL array, and each owns only the slice
!     sd(mypid*NSPP+1 : mypid*NSPP+NSPP,:). A thread's own slice is safe --
!     nothing but that thread writes it, and not before 4.b -- but the rest of
!     sd belongs to the other threads, and they reach 4.b on their own
!     schedule. Without the barrier a thread's zcnow is a MIXTURE of the
!     divergence at t and at t+dt, slice by slice, and which slices are which
!     depends on where each thread happened to be when the copy passed them.
!     That is what made terms 3 and 4 of the decomposition -- the only two that
!     read zcnow -- fail to reproduce across identical runs while the other
!     seven columns stayed bit-identical, and it put the sign of Cimp - Ct
!     inside the scatter. world-tqh4.
!
!     The barrier is conditional and that is safe: nenergy and nconvtime are
!     broadcast by mpbci, so the team enters this branch together or not at
!     all, and the configurations that set neither pay nothing.
      if (nenergy > 1 .or. nconvtime > 0) then
         allocate(zcnow(NESP,NLEV))
         zcnow(:,:) = sd(:,:)
!$omp barrier
      endif
!
!*    do the advective time step
!
      if(nadv > 0) then
!
!*    1. calculate divergence on timelevel t (solving eq.17),
!        which will replace the divergence tendency sdt
!
!     1.a precompute (g*script T; *zgt*) and (phi-phi*=g*T, eq.11; *zgm*)
!         needed for the rhs of eq.17.
!         (note that phi is needed in eq.17 and, therefor,
!         the surface geopotential phi* is added later (loop 1.b))
!
       do jlev = 1 , NLEV
         do jsp=1,NSPP
            zgt(jsp,jlev) = dot_product(g(:,jlev),stt(jsp,:))
            zgm(jsp,jlev) = dot_product(g(:,jlev),atm(jsp,:))
         enddo
       enddo
!
!     1.b compute divergence at time t (i.e. solve eq.17)
!         and overwrite *sdt* with the result
!         for comparison with HS75 note:
!
!         - *spt* contains -(!)script P
!         - *spm* contains log(ps)(t-dt)
!         - surface geopotential (*so*) needs to be added (see 1.a)
!         - bm1 is the invers of matrix (1/cn I+B dt**2) (lhs eq.17)
!         - zz is set to the rhs of eq.17
!         - zsum holds the result of the dot product of rhs(eq17) and bm1
!           (therefor jlev2-loop)
!
       do jlev = 1 , NLEV
        do jsp = 1 , NSPP
          jn = nindex(jsp)
          zsum = 0.0
          if (jn > 0) then
            zq = 1.0 / (jn * (jn+1))
            do jlev2 = 1 , NLEV
              z0 = t0(jlev2)
              zt = zgt(jsp,jlev2) - z0 * spt(jsp)
              zm = zgm(jsp,jlev2) + z0 * spm(jsp)
              za = adt(jsp,jlev2) * zq + sop(jsp)
              zb = adm(jsp,jlev2) * zq
              zz = zb + delt * (zm + za + delt * zt)
              zsum = zsum + zz * bm1(jlev2,jlev,jn)
            enddo
          endif
          sdt(jsp,jlev) = zsum
        enddo
       enddo
       if (mypid == NROOT) then
        sdt(1:2,:) = 0.0 ! first mode should be zero (errors due to numerics)
       endif
!
!*    2. calculate (-log) surface pressure tendency from eq.15 (pi=*dsigma*)
!
       do jlev = 1 , NLEV
         spt = spt + dsigma(jlev) * sdt(:,jlev)
       enddo
!
!*    3. calculate temperature tendency from eq.14
!
       do jlev = 1 , NLEV
        do jsp = 1 , NSPP
         stt(jsp,jlev)=stt(jsp,jlev)-dot_product(tau(:,jlev),sdt(jsp,:))
        enddo
       enddo
!
!     3a. THE REFERENCE CONVERSION'S TWO HALVES, PUT ON ONE TIME LEVEL. world-0ov.
!
!     The reference conversion is applied in two places: `calcgp` carries the
!     advective half `tkp*(zvgpg-ztpta)` at time t, and its divergence half is
!     the `tkp*c` part of `tau` just applied, on `sdt` -- the centred mean of
!     t-dt and t+dt, since `sdp = 2 sdt - adm`. What makes the two halves cancel
!     the momentum equations' reference pressure-gradient work is
!     `<ps (V.grad ln ps + D)> = 0`, the global integral of a mass-flux
!     divergence, and that identity holds at ONE time level and not across two.
!     The reference geopotential weights the split by up to 3.8 times t0, so the
!     displacement is 0.96 W/m2 against a sink of 0.79.
!
!     This puts the divergence half back on the divergence at t, which is what
!     `zcnow` holds. It is a CHANGE TO WHAT THE MODEL INTEGRATES: the semi-
!     implicit scheme no longer treats that half implicitly in the temperature
!     equation, while the divergence solve above still treats the temperature
!     implicitly, so the timestep this is stable at is its own question.
!
      if (nconvtime > 0) then
       do jlev = 1 , NLEV
        do jsp = 1 , NSPP
         jsg = jsp + mypid * NSPP
         zsum = 0.0
         do jlev2 = 1 , jlev
          zsum = zsum + tkp(jlev) * c(jlev2,jlev)                       &
     &                * (sdt(jsp,jlev2) - zcnow(jsg,jlev2))
         enddo
         stt(jsp,jlev) = stt(jsp,jlev) + zsum
        enddo
       enddo
      endif
!
      endif
!
!*    3b. THE ENERGY FIXER. THIS IS A CORRECTION AND IT IS NOT PHYSICS.
!
!     The adiabatic step is supposed to conserve total energy: the column
!     enthalpy it gives up must equal the kinetic energy it takes on. In this
!     model it does not. Measured on a dry adiabatic run it loses about 0.8
!     W/m2, and the whole of that is the REFERENCE conversion's two halves being
!     taken at different time levels. `calcgp` carries the advective half,
!     `tkp*(zvgpg-ztpta)`, at t; the divergence half is the `tkp*c` part of
!     `tau` in step 3 above, applied on `sdt`, which is the centred mean of t-dt
!     and t+dt. The identity that makes the two halves cancel the momentum
!     equations' reference pressure-gradient work is `<ps (V.grad ln ps + D)> =
!     0`, and it is split across those two time levels; the displacement is
!     -0.96 W/m2. That defect is world-0ov and it is NOT what this code fixes;
!     `nenergy = 2` is the control that measures it.
!
!     What this does is put the missing energy back, as a uniform warming, so
!     the atmosphere is not left short. It is the same device ECHAM, the IFS and
!     CAM all carry, and it has the same cost: it restores the TOTAL without
!     restoring where the energy went, so it MASKS the defect it compensates.
!     That is why the applied increment is reported rather than absorbed
!     silently -- a change in the underlying defect has to be able to show. This
!     fixer and its reporting are world-mzy.
!
!     The increment applied here was computed from the PREVIOUS step, at the end
!     of the diagnostic block below. Applying it to `stt` rather than to `stp`
!     is deliberate: the diagnostics further down read `stt`, so denergy01, 02
!     and 26 all see the correction and cannot disagree with the state.
!
      if (nenergyfix > 0 .and. mypid == NROOT) then
         stt(1,:) = stt(1,:) + denergyfix
      endif

!     3a. Coupling for synchronization runs

      if (mrnum == 2 .and. nsync > 0) then
         call mrdiff(stp,std,NESP,NLEV)
         call mrdiff(sdp,sdd,NESP,NLEV)
         call mrdiff(szp,szd,NESP,NLEV)
         call mrdiff(spp,spd,NESP,   1)
         stp(:,:) = stp(:,:) + syncstr * std(:,:)
         sdp(:,:) = sdp(:,:) + syncstr * sdd(:,:)
         szp(:,:) = szp(:,:) + syncstr * szd(:,:)
         spp(:  ) = spp(:  ) + syncstr * spd(:  )
      endif
!
!*    4. time stepping and time filtering (1st part)
!        the time filtering is splitted into two parts:
!
!        1st part: xf(prel)=x+eps*(xf(-)-2x)
!        2nd part: xf=xf(prel)+eps*x(+)
!
!        together the complete filter is xf=x+eps(xf(-)-2x+x(+))
!
!        with: x        = the prognostic variable to be filtered
!              xf       = the filtered variable (current time step)
!              xf(-)    = the filtered variable (old time step)
!              xf(prel) = a preliminary filterd variable
!              x(+)     = x at time (t+2dt)
!              eps      = filter constant = *pnu*
!
!     4.a 1st part of time filter
!
      if (nkits == 0) then
         spm = pnu21 * spp + pnu * apm
         sdm = pnu21 * sdp + pnu * adm
         szm = pnu21 * szp + pnu * azm
         stm = pnu21 * stp + pnu * atm
         if (nqspec == 1) sqm = pnu21 * sqp + pnu * aqm
      endif
!
!     4.b add tendencies
!
      if(nadv > 0 ) then
       spp = apm - delt2 * spt   ! log surface pressure (negative tendencies)
       sdp =   2.0 * sdt - adm   ! note that sdt=adm+delt*tendencies (see 1.)
       stp = delt2 * stt + atm   ! temperature
       szp = delt2 * szt + azm   ! vorticity
       if (nqspec == 1) sqp = delt2 * sqt + aqm   ! spec. humidity
      else
       spp = apm   ! log surface pressure (negative tendencies)
       sdp = adm   ! note that sdt=adm+delt*tendencies (see 1.)
       stp = atm   ! temperature
       szp = azm   ! vorticity
       if (nqspec == 1) sqp = aqm   ! spec. humidity
      endif
!
!     conserve p
!
      if (mypid == NROOT) then
       spp(1) = 0.0
       spp(2) = 0.0
      endif

!*    now the adiabatic time step is finished beside the time filtering.
!     2nd part of time filtering is done in subroutine spectrald
!
!     finaly the partial arrays are gathered from all processors (mpi)
!
!     These are the whole of the gather traffic. Under the shared build the
!     partial IS the slice, so mpgathersp is a barrier and moves nothing; under
!     MPI it is the allgather this used to be.
      call mpsyncsp
!
!     franks diagnostic
!
      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,3)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
!
!     energy diagnostics
!
      if(nenergy > 0) then
       allocate(ztt(NESP,NLEV))
       allocate(zst(NESP,NLEV))
       if (nqspec == 1) allocate(zsq(NESP,NLEV))
       allocate(zsp(NESP))
       call mpgallsp(ztt,stt,NLEV)
       call mpgallsp(zst,atm,NLEV)
       if (nqspec == 1) call mpgallsp(zsq,aqm,NLEV)
       call mpgallsp(zsp,apm,1)
       ztt(:,:)=ztt(:,:)*ct*ww
#ifdef OMPSHARED
       if (nshtns == 1) then
!        SHTns lands in GRID space, so sp2fl and its fc2gp both go. `sq` and
!        `sp` are the state at t+dt here -- the partials alias them -- so this
!        is not gridpointa's transform repeated, it is a different time level.
!
!        THE LEADING BARRIER keeps a thread's write to the whole globe clear of
!        a thread still reading last step's values out of the same scratch; the
!        trailing one is the wrappers' nowait contract. Six calls between one
!        pair, because the sources are not written in between and the
!        destinations are disjoint.
!$omp barrier
         call sh_sp2gp(ztt, zttgp_g, NLEV)
         if (nqspec == 1) call sh_sp2gp(sq, zqgp_g, NLEV)
         if (nqspec == 1) call sh_sp2gp(zsq, zqmgp_g, NLEV)
         call sh_sp2gp(zst, ztgp_g, NLEV)
         call sh_sp2gp(zsp, zpmgp_g, 1)
         call sh_sp2gp(sp, zpgp_g, 1)
!$omp barrier
         if (nqspec /= 1) then
            zqgp(:,:) = 0.0
            zqmgp(:,:) = 0.0
         endif
       else
       call sp2fl(ztt,zttgp,NLEV)
       if (nqspec == 1) call sp2fl(sq,zqgp,NLEV)
       if (nqspec == 1) call sp2fl(zsq,zqmgp,NLEV)
       call sp2fl(zst,ztgp,NLEV)
       call sp2fl(zsp,zpmgp,1)
       call sp2fl(sp,zpgp,1)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       if (nqspec == 1) then
          call fc2gp(zqgp,NLON,NLPP*NLEV)
          call fc2gp(zqmgp,NLON,NLPP*NLEV)
       else
          zqgp(:,:) = 0.0
          zqmgp(:,:) = 0.0
       endif
       call fc2gp(ztgp,NLON,NLPP*NLEV)
       call fc2gp(zpgp,NLON,NLPP)
       call fc2gp(zpmgp,NLON,NLPP)
       endif
#else
       call sp2fl(ztt,zttgp,NLEV)
       if (nqspec == 1) call sp2fl(sq,zqgp,NLEV)
       if (nqspec == 1) call sp2fl(zsq,zqmgp,NLEV)
       call sp2fl(zst,ztgp,NLEV)
       call sp2fl(zsp,zpmgp,1)
       call sp2fl(sp,zpgp,1)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       if (nqspec == 1) then
          call fc2gp(zqgp,NLON,NLPP*NLEV)
          call fc2gp(zqmgp,NLON,NLPP*NLEV)
       else
          zqgp(:,:) = 0.0
          zqmgp(:,:) = 0.0
       endif
       call fc2gp(ztgp,NLON,NLPP*NLEV)
       call fc2gp(zpgp,NLON,NLPP)
       call fc2gp(zpmgp,NLON,NLPP)
#endif
       zpmgp(:)=psurf*exp(zpmgp(:))
       zpgp(:)=psurf*exp(zpgp(:))
       do jlev=1,NLEV
        ztgp(:,jlev)=ct*(ztgp(:,jlev)+t0(jlev))
        if (nqspec == 1) zqmgp(:,jlev)=zqmgp(:,jlev)*psurf/zpmgp(:)
        if (nqspec == 1) zqgp(:,jlev)=zqgp(:,jlev)*psurf/zpgp(:)
       enddo
       deallocate(ztt)
       deallocate(zst)
       if (nqspec == 1) deallocate(zsq)
       deallocate(zsp)
      endif       
      if(nenergy > 0) then
       allocate(zsd(NESP,NLEV))
       allocate(zsz(NESP,NLEV))
       allocate(zekin(NHOR,NLEV))
       allocate(zepot(NHOR,NLEV))
       call mpgallsp(zsd,adm,NLEV)
       call mpgallsp(zsz,azm,NLEV)
#ifdef OMPSHARED
       if (nshtns == 1) then
!        One call for the pair, and it lands in grid space. `zsz` is ABSOLUTE
!        vorticity, which is what both transforms take: legmod removes
!        plavor's mode from its result and sh_dv2uv takes it off the
!        coefficient.
!$omp barrier
         call sh_dv2uv(zsd, zsz, zugp_g, zvgp_g, NLEV)
!$omp barrier
       else
       call dv2uv(zsd,zsz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
       endif
#else
       call dv2uv(zsd,zsz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
#endif
       do jlev=1,NLEV
        zekin(:,jlev)=0.5*(zugp(:,jlev)*zugp(:,jlev)                    &
     &                    +zvgp(:,jlev)*zvgp(:,jlev))*cv*cv*rcsq(:)     &
     &               *zpmgp(:)   
        zepot(:,jlev)=ztgp(:,jlev)*acpd*(1.+adv*zqmgp(:,jlev))*zpmgp(:) 
       enddo
#ifdef OMPSHARED
       if (nshtns == 1) then
!        THE LEADING BARRIER IS NOT SYMMETRY. The pair above is written to the
!        SAME globe, and the loop just above reads it: without this a thread
!        that arrives early overwrites the wind at t-dt for a thread still
!        forming zekin from it. `sd` and `sz` here are the state at t+dt --
!        `sdp` and `szp` alias them, so the tendencies applied above have
!        already advanced them -- and the difference of the two kinetic
!        energies is what denergy27 is.
!$omp barrier
         call sh_dv2uv(sd, sz, zugp_g, zvgp_g, NLEV)
!$omp barrier
       else
       call dv2uv(sd,sz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
       endif
#else
       call dv2uv(sd,sz,zugp,zvgp)
       call fc2gp(zugp,NLON,NLPP*NLEV)
       call fc2gp(zvgp,NLON,NLPP*NLEV)
#endif
       denergy(:,27)=0.
       denergy(:,26)=0.
       denergy(:,2)=0.
       denergy(:,1)=0.
       do jlev=1,NLEV
        zekin(:,jlev)=0.5*(zugp(:,jlev)*zugp(:,jlev)                    &
     &                    +zvgp(:,jlev)*zvgp(:,jlev))*cv*cv*rcsq(:)     &
     &               *zpgp(:)                                           &
     &               -zekin(:,jlev)
        denergy(:,27)=denergy(:,27)                                     &
     &               -zekin(:,jlev)/deltsec2/ga*dsigma(jlev)
        denergy(:,1)=denergy(:,1)                                       &
     &              +(ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)              &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)*dsigma(jlev)/ga
        denergy(:,2)=denergy(:,2)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)/ga*dsigma(jlev)
        denergy(:,26)=denergy(:,26)                                     &
     &               +((ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)            &
     &                 *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)              &
     &                -zepot(:,jlev))/deltsec2/ga*dsigma(jlev)
        if(nener3d > 0) then 
         dener3d(:,jlev,27)=-zekin(:,jlev)/deltsec2/ga*dsigma(jlev)
         dener3d(:,jlev,1)=(ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)        &
     &               *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)*dsigma(jlev)/ga
         dener3d(:,jlev,2)=zttgp(:,jlev)                                &
     &              *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)/ga*dsigma(jlev)
         dener3d(:,jlev,26)=((ztgp(:,jlev)+zttgp(:,jlev)*deltsec2)      &
     &                     *acpd*(1.+adv*zqgp(:,jlev))*zpgp(:)          &
     &                     -zepot(:,jlev))/deltsec2/ga*dsigma(jlev)
        endif  
       enddo
!
!      THE ENERGY FIXER, second half: what to apply next step. See 3b above.
!
!      `denergy26 - denergy27` is the enthalpy the adiabatic step gave up less
!      the kinetic energy it took on, which must be zero and is not. The
!      correction is spread as a uniform temperature increment over the column
!      heat capacity, so the energy it returns is proportional to the local mass.
!
!      This is an integral controller with unit gain, not a one-shot correction.
!      The increment is already inside the imbalance measured above, so
!      subtracting the residual leaves exactly minus the raw imbalance: it
!      settles in one step and then tracks. Reading `denergy26 - denergy27` after
!      a run with the fixer on therefore reports the RESIDUAL, near zero, and the
!      size of the defect is `denergyfix` itself.
!
       if(nenergyfix > 0) then
        zfix(:) = 0.0
        jhor = 0
        do jlat = 1 , NLPP
         do jlon = 1 , NLON
          jhor = jhor + 1
          zfixw(jhor) = gwd(jlat)
         enddo
        enddo
!       THE TARGET IS THE SUM THAT MUST VANISH, not the adiabatic step alone.
!       Total enthalpy change is denergy26 (spectrala) plus 24, 23 and 25
!       (spectrald); total kinetic change is -denergy27 less the friction that
!       23 and 25 book back as heat. The friction cancels by construction --
!       mkdheat returns exactly what it removed -- so what is left to vanish is
!       `26 - 27 + 24`. Targeting `26 - 27` alone zeroes the conversion defect
!       and leaves the hyperdiffusion's +0.24 W/m2 as a net GAIN, which is
!       measured: the trend went from -0.66 to +0.22.
        zfix(1) = dot_product(denergy(:,26)-denergy(:,27),zfixw)         &
     &          + denergyd24*sum(zfixw)
        do jlev = 1 , NLEV
         zfix(2) = zfix(2) + dot_product(acpd*(1.+adv*zqgp(:,jlev))     &
     &                       *zpgp(:)/ga*dsigma(jlev),zfixw)
        enddo
        zfix(3) = sum(zfixw)
        call mpsumbcr(zfix,3)
        if (mypid == NROOT) then
!        WINDOWED, over one model day, and this is the point of the design.
!        The PER-STEP imbalance swings by about 250 W/m2 either way -- that is
!        the leapfrog's computational mode, damped at pnu = 0.1 and not gone --
!        while the thing being corrected is a systematic loss of order half a
!        watt. A controller that reacts step by step converges in the mean and
!        wanders across -0.2 to +1.8 getting there, which is what it did.
!        Averaging over ntspd steps takes the mode out and leaves the defect.
!
!        Only NROOT touches these accumulators, and only after the reduction,
!        which is what keeps them free of a race under the threaded build.
         denergyacc(1) = denergyacc(1) + zfix(1)
         denergyacc(2) = denergyacc(2) + zfix(2)
         denergyacc(3) = denergyacc(3) + zfix(3)
         nenergyacc = nenergyacc + 1
         if (nenergyacc >= ntspd) then
!           THE FIRST WINDOW IS DISCARDED. It carries the start-up transient --
!           the first step out of a restart shows an imbalance of order 250
!           W/m2 -- and one such sample still moves a 64-step mean by four.
            if (nenergywin > 0) then
!              A TENDENCY, because stt is one. Writing an INCREMENT here makes
!              only delt2 of it land, so the controller never sees its own
!              correction arrive and winds up without bound.
               zfixd = -denergyacc(1)/denergyacc(2)/(ct*ww)
!              Unit gain on a clean window average converges in ONE window; the
!              limit is a backstop against a pathological one, not a brake.
               zfixc = zfixd*ct*ww*denergyacc(2)/denergyacc(3)
               if (abs(zfixc) > 2.0) zfixd = zfixd*2.0/abs(zfixc)
               denergyfix = denergyfix + zfixd
            endif
            nenergywin = nenergywin + 1
            denergyacc(:) = 0.0
            nenergyacc = 0
         endif
         zfixr = denergyfix*ct*ww*zfix(2)/zfix(3)
!        A runaway has to announce itself rather than appear as a blow-up in the
!        dynamics. The defect being corrected is of order 1 W/m2; this bound is
!        two orders above it and catches divergence without firing on a start-up
!        transient.
         if (abs(zfixr) > 100.0) then
          write(nud,*) 'ENERGY FIXER DIVERGED: applying ',zfixr,' W/m2 at step ',nstep
          write(nud,*) 'the correction is not converging; see world-mzy'
          stop 'energy fixer diverged'
         endif
         if (mod(nstep,ndiag) == 0 .or. nstep < nstep1+40) then
          write(nud,'(A,I8,3E16.7)') ' ENERGY FIXER applied W/m2, K/day ', &
     &      nstep, zfixr, denergyfix*ct*ww*86400.0, zfix(1)/zfix(3)
         endif
        endif
       endif
!
!      THE CONVERSION DECOMPOSITION. A CONTROL AND NOT A MODEL TERM. world-0ov.
!
!      Everything here is measured in denergy02's own arithmetic and units, so
!      the pieces add to the terms the budget is already written in. What is new
!      is the IMPLICIT half of the reference conversion: `calcgp` carries only
!      `tkp*(zvgpg-ztpta)`, the advective half, and the divergence half is the
!      `tkp*c` part of `tau` applied in step 3 above. A budget that reads gtn and
!      takes the remainder by difference books that half to the advection, and
!      the two halves are of opposite sign and comparable size, so the sign of
!      the attribution depends on keeping them together.
!
!      Term 6 is the same conversion evaluated on the divergence at time t
!      rather than on `sdt`. `sdt` is exactly the centred mean of t-dt and t+dt,
!      the explicit half is at t, and 6 minus 4 is therefore the semi-implicit
!      displacement between the conversion's two ends on its own.
!
       if (nenergy > 1) then
        allocate(zcsdt(NESP,NLEV))
        allocate(zcwrk(NESP,NLEV))
        allocate(zcgp(NHOR,NLEV))
        call mpgallsp(zcsdt,sdt,NLEV)
        jhor = 0
        do jlat = 1 , NLPP
         do jlon = 1 , NLON
          jhor = jhor + 1
          zcw(jhor) = gwd(jlat)
         enddo
        enddo
        zcs(:) = 0.0
        zcs(1) = dot_product(denergy(:,2),zcw)
        zcs(2) = dot_product(denergy(:,26),zcw)
        zcs(3) = dot_product(denergy(:,27),zcw)
        zcs(10) = sum(zcw)
        do jterm = 1 , 6
         do jlev = 1 , NLEV
          zcwrk(:,jlev) = 0.0
          do jlev2 = 1 , NLEV
           if (jterm == 1) then
            if (jlev2 <= jlev) zcwrk(:,jlev) = zcwrk(:,jlev)            &
     &         - tkp(jlev) * c(jlev2,jlev) * zcsdt(:,jlev2)
           else if (jterm == 2) then
            zcwrk(:,jlev) = zcwrk(:,jlev)                               &
     &         - tau(jlev2,jlev) * zcsdt(:,jlev2)
            if (jlev2 <= jlev) zcwrk(:,jlev) = zcwrk(:,jlev)            &
     &         + tkp(jlev) * c(jlev2,jlev) * zcsdt(:,jlev2)
           else if (jterm == 3) then
            if (jlev2 <= jlev) zcwrk(:,jlev) = zcwrk(:,jlev)            &
     &         - tkp(jlev) * c(jlev2,jlev) * zcnow(:,jlev2)
           else if (jterm == 5) then
            if (jlev2 <= jlev) zcwrk(:,jlev) = zcwrk(:,jlev)            &
     &         - tkp(jlev) * c(jlev2,jlev) * zsd(:,jlev2)
           else if (jterm == 6) then
            if (jlev2 <= jlev) zcwrk(:,jlev) = zcwrk(:,jlev)            &
     &         - tkp(jlev) * c(jlev2,jlev) * sd(:,jlev2)
           endif
          enddo
          if (jterm == 4) zcwrk(:,jlev) = -tkp(jlev) * zcnow(:,jlev)
         enddo
         zcwrk(:,:) = zcwrk(:,:) * ct * ww
         call sp2fl(zcwrk,zcgp,NLEV)
         call fc2gp(zcgp,NLON,NLPP*NLEV)
         do jlev = 1 , NLEV
          zcs(3+jterm) = zcs(3+jterm)                                   &
     &     + dot_product(zcgp(:,jlev)*acpd*(1.+adv*zqgp(:,jlev))        &
     &                   *zpgp(:)/ga*dsigma(jlev),zcw)
         enddo
        enddo
        call mpsumbcr(zcs,10)
        if (mypid == NROOT) then
!        THE FIRST MODEL DAY IS DROPPED, for the reason the fixer's window drops
!        its first: the first step out of a restart carries an imbalance of order
!        250 W/m2, and this is a mean over an orbit.
         if (nstep > nstep1 + ntspd) then
          dconvacc(1:9) = dconvacc(1:9) + zcs(1:9) / zcs(10)
          dconvacc(10) = dconvacc(10) + 1.0
          nconvacc = nconvacc + 1
         endif
         if (mod(nstep,ndiag) == 0 .and. nconvacc > 0) then
          write(nud,'(A,I9,I9,9E15.6)') ' CONVDECOMP d02 d26 d27 cimp '//&
     &      'cvadv ct dt ctm ctp ', nstep, nconvacc,                    &
     &      dconvacc(1:9) / dconvacc(10)
         endif
        endif
        deallocate(zcsdt)
        deallocate(zcwrk)
        deallocate(zcgp)
       endif
       deallocate(zsd)
       deallocate(zsz)
       deallocate(zekin)
       deallocate(zepot)
      endif
      if (nenergy > 1 .or. nconvtime > 0) deallocate(zcnow)
!
      return
      end

!-Diabatic-subroutines  (Frank (Larry) 26-Nov-99)

!     =====================
!     SUBROUTINE GRIDPOINTD
!     =====================

      subroutine gridpointd
      use pumamod
#ifdef OMPSHARED
      use shtnsmod, only: sh_sp2gp, sh_dv2uv, sh_gp2sp, sh_uv2dv, sh_slice
#endif
!
!     The same tendency partials gridpointa uses, from pumamod. The two
!     routines do not overlap -- master runs one and then the other -- so one
!     set of slots serves both.

!
!*    Diabatic Gridpoint Calculations
!
      gudt(:,:)=0.
      gvdt(:,:)=0.
      gtdt(:,:)=0.
      gqdt(:,:)=0.
      dudt(:,:)=0.
      dvdt(:,:)=0.
      dtdt(:,:)=0.
      dqdt(:,:)=0.
      mmrt(:,:)=0.

!
!     transform to gridpoint domain
!

#ifdef OMPSHARED
      if (nshtns == 1) then
!        The leading barrier is required for the same reason as in gridpointa:
!        these wrappers write the whole globe, not this thread's band.
!$omp barrier
!        As in gridpointa: one call per field, landing in GRID space, so there
!        is no Fourier intermediate and no fc2gp. This site is the easy one --
!        nothing between the transform and the dimensionalising below reads a
!        Fourier coefficient. invlegd leaves gd and gz alone and so does this.
         call sh_dv2uv(sd, sz, gu_g, gv_g, NLEV)
         call sh_sp2gp(st, gt_g, NLEV)
         call sh_sp2gp(sp, gp_g, 1)
         if (nqspec == 1) call sh_sp2gp(sq, gq_g, NLEV)
!        The wrappers do not synchronise; the caller does.
!$omp barrier
      else
         call invlegd
         call fc2gp(gu  ,NLON,NLPP*NLEV)
         call fc2gp(gv  ,NLON,NLPP*NLEV)
         call fc2gp(gt  ,NLON,NLPP*NLEV)
         call fc2gp(gp  ,NLON,NLPP)
         if (nqspec == 1) call fc2gp(gq  ,NLON,NLPP*NLEV)
      endif
#else
      call invlegd

      call fc2gp(gu  ,NLON,NLPP*NLEV)
      call fc2gp(gv  ,NLON,NLPP*NLEV)
      call fc2gp(gt  ,NLON,NLPP*NLEV)
      call fc2gp(gp  ,NLON,NLPP)
      if (nqspec == 1) then
         call fc2gp(gq  ,NLON,NLPP*NLEV)
      endif
#endif

!
!     dimensionalize
!

      dp(:)=psurf*exp(gp(:))
      do jlev = 1 , NLEV
       du(:,jlev)=cv*gu(:,jlev)*SQRT(rcsq(:))
       dv(:,jlev)=cv*gv(:,jlev)*SQRT(rcsq(:))
       dt(:,jlev)=ct*(gt(:,jlev)+t0(jlev))
       if (nqspec == 1) dq(:,jlev)=gq(:,jlev)*psurf/dp(:)
      enddo
!
!     tracer transport
!
      if (nsela > 0 .and. nkits == 0 .and. nqspec == 0) then	  
       if (nqspec == 0) then
!       ROOT ONLY, and the flush with it. nud is one Fortran unit shared by
!       the whole thread team, so an unguarded line here is NPRO copies per
!       timestep in whatever order the threads reach it, and NPRO concurrent
!       flushes of one unit. The message reports a global switch, so root's
!       copy is the whole of it. world-0ihs.
        if (mypid == NROOT) then
         write(nud,*) 'Semi-Lagrangian q running'
         flush(nud)
        endif
        dqt(:,:) = dq(:,:) ! Save old value of q
        call mpgagp(zgq,dq,NLEV)
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          dtrace(:,NLAT+1-jlat,:,1) = zgq(:,jlat,:)
         enddo ! jlat
        endif ! mypid
       endif ! nqspec
       call tracer_main
       if (nqspec == 0) then
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          zgq(:,jlat,:) = dtrace(:,NLAT+1-jlat,:,1)
         enddo
        endif ! mypid
        call mpscgp(zgq,dq,NLEV)
        dqt(:,:) = (dq(:,:) - dqt(:,:)) / deltsec !  q advection term
       endif ! nqspec	  
      endif ! nkits
      call guihor("DQ" // char(0),dq,NLEV,1000.0,0.0)

!     aerosol transport
  
      if (nsela == 1 .and. nkits == 0 .and. l_aero > 0) then
        mmrt(:,:) = mmr(:,:) ! Save old value of mmr
        call mpgagp(zmmr,mmr,NLEV)
        call mpgagp(znrho,nrho,NLEV)
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          daeros(:,NLAT+1-jlat,:,1) = zmmr(:,jlat,:)
          numrhos(:,NLAT+1-jlat,:,1) = znrho(:,jlat,:)
         enddo ! jlat
        endif ! mypid
        call aero_main
        if (mypid == NROOT) then
         do jlat = 1 , NLAT
          zmmr(:,jlat,:) = daeros(:,NLAT+1-jlat,:,1)
          znrho(:,jlat,:) = numrhos(:,NLAT+1-jlat,:,1)
         enddo
        endif ! mypid
        call mpscgp(zmmr,mmr,NLEV)
        call mpscgp(znrho,nrho,NLEV)
        mmrt(:,:) = (mmr(:,:) - mmrt(:,:)) / deltsec !  q advection term       endif !  
      endif ! nkits
!
!     compute output specific humidity
!

!     Under low I/O, outaccu sums sqout on EVERY timestep, so sqout has to be
!     current on every timestep. Recomputing it only at the output step made
!     the accumulated humidity a sum over one stale value per interval, and
!     zero for the steps of a model call that precede the first update -- the
!     phase of which walks from call to call. Only the mod() test changes; the
!     transform itself is unchanged.
      if ((nlowio > 0 .or. mod(nstep,nafter)==0) .and. nqspec == 1) then
       do jlev=1,NLEV
        zqout(:,jlev)=gq(:,jlev)/exp(gp(:))
       enddo
#ifdef OMPSHARED
       if (nshtns == 1) then
!         THE REDUCTION GOES, and the copy replacing it is not optional. sqout
!         is threadprivate and the legmod branch below fills each copy with the
!         partial from that thread's own latitudes, which mpsum turns into the
!         whole field. sh_gp2sp returns the whole field already, so mpsum would
!         multiply it by the thread count; and because the wrapper is parallel
!         over LEVELS, writing straight into a threadprivate array would leave
!         each copy holding only the levels its own thread owned.
!$omp barrier
          call sh_gp2sp(zqout_g, sqout_g, NLEV)
!$omp barrier
          sqout(:,:) = sqout_g(:,:)
       else
       do jlev=1,NLEV
        call gp2fc(zqout(:,jlev),NLON,NLPP)
        call fc2sp(zqout(:,jlev),sqout(1,jlev))
       enddo
       call mpsum(sqout,NLEV)
       endif
#else
       do jlev=1,NLEV
        call gp2fc(zqout(:,jlev),NLON,NLPP)
        call fc2sp(zqout(:,jlev),sqout(1,jlev))
       enddo
       call mpsum(sqout,NLEV)
#endif
      endif

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
          write(nud,*)'------------ SUBROUTINE GRIDPOINTD -------------'
          write(nud,*)'NSTEP= ',nstep
          write(nud,*)'YY,MM,DD,HH,MIN= ',ndatim(1:5)
          write(nud,*) 'gp= ',MAXVAL(gp(1:NHOR)),MINVAL(gp(1:NHOR))
          write(nud,*) 'dp= ',MAXVAL(dp(1:NHOR)),MINVAL(dp(1:NHOR))
       endif
       call prdbug1
      endif

!
!     PARAMETERIZATION ROUTINES TO BE INCLUDED BELOW
!
!     a) miscellaneous parameterizations (like q-fixer)
!

      call miscstep

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*) 'After MISC:'
       endif
       call prdbug2
      endif

!
!     b) surface fluxes and vertical diffusion
!

      call fluxstep
     
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After FLUX:'
       endif
       call prdbug2
      endif

!
!     c) radiation
!

      if(nrad > 0) call radstep
     
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After RAD:'
       endif
       call prdbug2
      endif

!
!     d) convective and large scale rain and clouds
!

      call rainstep

!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'After RAIN:'
       endif
       call prdbug2
      endif

!
!     e) surface parmeterizations or models
!

      call surfstep

!
!     f) carbon-silicate weathering parameterization
!

      call carbonstep
      
      
!
!     g) newtonian cooling of stratosphere
!
      
      if (nstratosponge > 0) then
         mint(:) = minval(dt(:,:NLEV-12),2)
         do k=1,NLEV-12
            dtdt(:,k) = dtdt(:,k) + (mint(:) - (dt(:,k)+dtdt(:,k)))/(taucool*day_24hr)
         enddo
      endif
      
!
!     h) hurricane/storm diagnostics
!
      
      call hurricanestep   

!
!     END OF PARAMETERISATION ROUTINES
!
!     dbug print out
!

      if(nprint > 0) then
       if(mypid==NROOT) then
        write(nud,*)'FINAL:'
       endif
       call prdbug2
      endif

!
!     franks diagnostic
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,1)=dtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,11)=dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000
       enddo
      endif
!
!     energy diagnostics
!
      if(nenergy > 0) then
       denergy(:,4)=0.
       do jlev=1,NLEV
        denergy(:,4)=denergy(:,4)                                       &
     &              +dtdt(:,jlev)                                       &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,4)=dtdt(:,jlev)                                 &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo
      endif
!
!     de-dimensionalize
!

      do jlev = 1 , NLEV
       gudt(:,jlev)=dudt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gudt(:,jlev)
       gvdt(:,jlev)=dvdt(:,jlev)/SQRT(rcsq(:))/cv/ww                    &
     &             +gvdt(:,jlev)
       gtdt(:,jlev)=dtdt(:,jlev)/ct/ww                                  &
     &             +gtdt(:,jlev)
       gqdt(:,jlev)=dqdt(:,jlev)*dp(:)/ww/psurf                         &
     &             +gqdt(:,jlev)
      enddo

!
!     transform to spectral space
!

#ifdef OMPSHARED
      if (nshtns == 1) then
!        The last transform site. As in gridpointa: the finished field lands in
!        slot 0 of the partial scratch and each thread takes its slice, and the
!        four reductions go away because there is nothing partial to reduce.
!$omp barrier
         call sh_gp2sp(gtdt_g, zpst(1,1,0), NLEV)
         call sh_uv2dv(gudt_g, gvdt_g, zpsd(1,1,0), zpsz(1,1,0), NLEV)
         if (nqspec == 1) call sh_gp2sp(gqdt_g, zpsq(1,1,0), NLEV)
!$omp barrier
         call sh_slice(zpst(1,1,0), stt, NLEV)
         call sh_slice(zpsd(1,1,0), sdt, NLEV)
         call sh_slice(zpsz(1,1,0), szt, NLEV)
         if (nqspec == 1) call sh_slice(zpsq(1,1,0), sqt, NLEV)
      else
#endif
      call gp2fc(gtdt,NLON,NLPP*NLEV)
      call gp2fc(gudt,NLON,NLPP*NLEV)
      call gp2fc(gvdt,NLON,NLPP*NLEV)
      if (nqspec == 1) call gp2fc(gqdt,NLON,NLPP*NLEV)

!
!     direct Legendre transformation (fourier domain to spectral domain)
!
      do jlev = 1 , NLEV
         call fc2sp(gtdt(:,jlev),zpst(1,jlev,mypart))
         if (nqspec == 1) call fc2sp(gqdt(:,jlev),zpsq(1,jlev,mypart))
      enddo

      call uv2dv(gudt,gvdt,zpsd(1,1,mypart),zpsz(1,1,mypart))

      call mpsumscp(zpst,stt,NLEV)
      call mpsumscp(zpsd,sdt,NLEV)
      call mpsumscp(zpsz,szt,NLEV)
      if (nqspec == 1) call mpsumscp(zpsq,sqt,NLEV)
#ifdef OMPSHARED
      endif
#endif
      if (nqspec == 0) dq(:,:) = dq(:,:) + dqdt(:,:) * deltsec
      if (nsela == 1 .and. l_aero > 0) mmr(:,:) = mmr(:,:) + mmrt(:,:) * deltsec

      return
      end

!     ==================
!     SUBROUTINE PRDBUG1
!     ==================

      subroutine prdbug1
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)
      real zfpr5(NLON*NLAT,NLEV)
      real zfpr6(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dt,NLEV)
      call mpgagp(zfpr2,dq,NLEV)
      call mpgagp(zfpr3,du,NLEV)
      call mpgagp(zfpr4,dv,NLEV)
      call mpgagp(zfpr5,dcc,NLEV)
      call mpgagp(zfpr6,dql,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         write(nud,*)'L= ',jlev,' MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN T= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         write(nud,*)'L= ',jlev,' MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN Q= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         write(nud,*)'L= ',jlev,' MAX U= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN U= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         write(nud,*)'L= ',jlev,' MAX V= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN V= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr5(:,jlev))
         kmaxl=MAXLOC(zfpr5(:,jlev))
         zminv=MINVAL(zfpr5(:,jlev))
         kminl=MINLOC(zfpr5(:,jlev))
         write(nud,*)'L= ',jlev,' MAX C= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN C= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr6(:,jlev))
         kmaxl=MAXLOC(zfpr6(:,jlev))
         zminv=MINVAL(zfpr6(:,jlev))
         kminl=MINLOC(zfpr6(:,jlev))
         write(nud,*)'L= ',jlev,' MAX L= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN L= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' T= ',zfpr1(nprhor,jlev)                    &
     &                    ,' Q= ',zfpr2(nprhor,jlev)
         write(nud,*)'L= ',jlev,' U= ',zfpr3(nprhor,jlev)                    &
     &                    ,' V= ',zfpr4(nprhor,jlev)
         write(nud,*)'L= ',jlev,' C= ',zfpr5(nprhor,jlev)                    &
     &                    ,' L= ',zfpr6(nprhor,jlev)
        enddo
       endif
      endif

      call mpgagp(zfpr1(1,1),dls,1)
      call mpgagp(zfpr2(1,1),dp,1)
      call mpgagp(zfpr3(1,1),drhs,1)
      call mpgagp(zfpr4(1,1),dalb,1)
      call mpgagp(zfpr5(1,1),dt(1,NLEP),1)
      call mpgagp(zfpr6(1,1),dq(1,NLEP),1)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostic:'
        zmaxv=MAXVAL(zfpr5(:,1))
        kmaxl=MAXLOC(zfpr5(:,1))
        zminv=MINVAL(zfpr5(:,1))
        kminl=MINLOC(zfpr5(:,1))
        write(nud,*)'SURF MAX T= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN T= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr6(:,1))
        kmaxl=MAXLOC(zfpr6(:,1))
        zminv=MINVAL(zfpr6(:,1))
        kminl=MINLOC(zfpr6(:,1))
        write(nud,*)'SURF MAX Q= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN Q= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr2(:,1))
        kmaxl=MAXLOC(zfpr2(:,1))
        zminv=MINVAL(zfpr2(:,1))
        kminl=MINLOC(zfpr2(:,1))
        write(nud,*)'SURF MAX P= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN P= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr3(:,1))
        kmaxl=MAXLOC(zfpr3(:,1))
        zminv=MINVAL(zfpr3(:,1))
        kminl=MINLOC(zfpr3(:,1))
        write(nud,*)'SURF MAX RHS= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN RHS= ',zminv,' NHOR= ',kminl(1)
        zmaxv=MAXVAL(zfpr4(:,1))
        kmaxl=MAXLOC(zfpr4(:,1))
        zminv=MINVAL(zfpr4(:,1))
        kminl=MINLOC(zfpr4(:,1))
        write(nud,*)'SURF MAX ALB= ',zmaxv,' NHOR= ',kmaxl(1)
        write(nud,*)'SURF MIN ALB= ',zminv,' NHOR= ',kminl(1)
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        write(nud,*)'SURF: LS= ',zfpr1(nprhor,1)                             &
     &             ,' PS= ',zfpr2(nprhor,1)
        write(nud,*)'SURF: RHS= ',zfpr3(nprhor,1)                            &
     &             ,' ALB= ',zfpr4(nprhor,1)
        write(nud,*)'SURF: TS= ',zfpr5(nprhor,1)                             &
     &             ,' QS= ',zfpr6(nprhor,1)
       endif
      endif

      if(nprint==2) then
       call mpgagp(zfpr1(1,1),dicec,1)

       if(mypid==NROOT) then
        write(nud,*)'SURF: ICEC= ',zfpr1(nprhor,1)
       endif
      endif

      return
      end subroutine prdbug1

!     ==================
!     SUBROUTINE PRDBUG2
!     ==================

      subroutine prdbug2
      use pumamod

      integer kmaxl(1),kminl(1)

      real zfpr1(NLON*NLAT,NLEV)
      real zfpr2(NLON*NLAT,NLEV)
      real zfpr3(NLON*NLAT,NLEV)
      real zfpr4(NLON*NLAT,NLEV)

      call mpgagp(zfpr1,dtdt,NLEV)
      call mpgagp(zfpr2,dqdt,NLEV)
      call mpgagp(zfpr3,dudt,NLEV)
      call mpgagp(zfpr4,dvdt,NLEV)

      if(mypid==NROOT) then
       if(nprint==1) then
        write(nud,*)'Global diagnostics: '
        do jlev=1,NLEV
         zmaxv=MAXVAL(zfpr1(:,jlev))
         kmaxl=MAXLOC(zfpr1(:,jlev))
         zminv=MINVAL(zfpr1(:,jlev))
         kminl=MINLOC(zfpr1(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DT= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DT= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr2(:,jlev))
         kmaxl=MAXLOC(zfpr2(:,jlev))
         zminv=MINVAL(zfpr2(:,jlev))
         kminl=MINLOC(zfpr2(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DQ= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DQ= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr3(:,jlev))
         kmaxl=MAXLOC(zfpr3(:,jlev))
         zminv=MINVAL(zfpr3(:,jlev))
         kminl=MINLOC(zfpr3(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DU= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DU= ',zminv,' NHOR= ',kminl(1)
         zmaxv=MAXVAL(zfpr4(:,jlev))
         kmaxl=MAXLOC(zfpr4(:,jlev))
         zminv=MINVAL(zfpr4(:,jlev))
         kminl=MINLOC(zfpr4(:,jlev))
         write(nud,*)'L= ',jlev,' MAX DV= ',zmaxv,' NHOR= ',kmaxl(1)
         write(nud,*)'L= ',jlev,' MIN DV= ',zminv,' NHOR= ',kminl(1)
        enddo
       elseif(nprint==2) then
        write(nud,*)'Local diagnostics at nhor= ',nprhor,': '
        do jlev=1,NLEV
         write(nud,*)'L= ',jlev,' DT= ',zfpr1(nprhor,jlev)                   &
     &                    ,' DQ= ',zfpr2(nprhor,jlev)
         write(nud,*)'L= ',jlev,' DU= ',zfpr3(nprhor,jlev)                   &
     &                    ,' DV= ',zfpr4(nprhor,jlev)
        enddo
       endif
      endif

      return
      end subroutine prdbug2

!     ====================
!     SUBROUTINE SPECTRALD
!     ====================

      subroutine spectrald
      use pumamod
#ifdef OMPSHARED
      use shtnsmod, only: sh_sp2gp
#endif
!
      real :: zsdt1(NSPP,NLEV),zsdt2(NSPP,NLEV)
      real :: zszt1(NSPP,NLEV),zszt2(NSPP,NLEV)
      real :: zsum3(1) ! must be an array because of the NAG compiler
      real :: zd24(2) ! the fixer's carry of the temperature diffusion heating
!
!     franks diagnostics
!
      real, allocatable :: ztt(:,:)
      real, allocatable :: zst(:,:),zsq(:,:),zstt(:,:),zstt2(:,:)
      real, allocatable :: zsttd(:,:)
      real, allocatable :: zdtgp(:,:)
      real, allocatable :: zgw(:),zsum1(:)
!
!     zttgp, ztgp and zqgp are in pumamod, as bands of full-globe arrays --
!     spectrala's twins, and this routine runs after it, so one set of slots
!     serves both. zdtgp stays here: the ndheat > 1 efficiency block that
!     fills it is not converted, ndheat defaulting to 1 and nothing in this
!     project raising it.
!
!     prepare diagnostics of efficiency
!
!     The shared spectral arrays are read in FULL by the phase before this one
!     and written by SLICE in this one, and nothing else separates the two.
!     Without this a thread arriving early overwrites what another is still
!     reading. Inert without -fopenmp.
!$omp barrier

      if(ndheat > 1) then
       allocate(zst(NESP,NLEV))      
       if (nqspec == 1) allocate(zsq(NESP,NLEV))      
       allocate(zstt(NESP,NLEV)) 
       allocate(zstt2(NESP,NLEV)) 
       allocate(zdtgp(NHOR,NLEV)) 
       allocate(zsum1(4))
       allocate(zgw(NHOR))
       call mpgallsp(zst,stp,NLEV)
       call mpgallsp(zstt,stt,NLEV)
       if (nqspec == 1) call mpgallsp(zsq,sqp,NLEV)
       jhor=0
       do jlat=1,NLPP
        do jlon=1,NLON
         jhor=jhor+1
         zgw(jhor)=gwd(jlat)
        enddo
       enddo
      endif

!
!     add tendencies from diabatic parameterizations
!
!     THE CONVERSION EITHER SIDE OF EACH WRITE TO sdp. A CONTROL, world-pkf.
!     `spectrald` moves the reference conversion's divergence half by about two
!     watts every step and that displacement is the whole of the adiabatic sink;
!     these three prints say which of the two writes does it. `dsdiv` allocated
!     only under the control, so this is a no-op below nenergy = 2.
      if (nenergy > 1) then
         allocate(dsdiv(NESP,NLEV))
         call mpgallsp(dsdiv,sdp,NLEV)
         call convweight(dsdiv,dq,dp,dconvspd(1))
      endif

      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sdp = sdp + delt2 * sdt
      if (nqspec == 1) sqp = sqp + delt2 * sqt

      if (nenergy > 1) then
         call mpgallsp(dsdiv,sdp,NLEV)
         call convweight(dsdiv,dq,dp,dconvspd(2))
      endif

!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,1)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
      if(ndiaggp==1) then
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,20)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif
!
!     energy diagnostics
!
      if(nenergy > 0) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
#ifdef OMPSHARED
       if (nshtns == 1) then
!        SHTns lands in GRID space, so sp2fl and its fc2gp both go. The
!        barriers are the wrappers' contract: they write the whole globe and
!        return without synchronising.
!$omp barrier
         call sh_sp2gp(ztt, zttgp_g, NLEV)
!$omp barrier
       else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       endif
#else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
#endif
       deallocate(ztt)
       denergy(:,3)=0.
       do jlev=1,NLEV
        denergy(:,3)=denergy(:,3)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,3)=zttgp(:,jlev)                                &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
      endif

!     calculates spectral tendencies from restoration (if included)
!     and biharmonic diffusion (both implicitely).
!     add newtonian cooling and drag.

      do jlev = 1 , NLEV
       stt(:,jlev)=((-damp(jlev)-tdisst(jlev)*sakpp(1:NSPP,jlev))       &
     &              *stp(:,jlev)+damp(jlev)*srp(:,jlev))                &
     &          /(1.+delt2*(damp(jlev)+tdisst(jlev)*sakpp(1:NSPP,jlev)))
       if (nqspec == 1)                                                 &
     & sqt(:,jlev)=-tdissq(jlev)*sakpp(1:NSPP,jlev)*sqp(:,jlev)         &
     &            /(1.+delt2*tdissq(jlev)*sakpp(1:NSPP,jlev))
       zsdt1(:,jlev)=-tfrc(jlev)*sdp(:,jlev)                            &
     &              /(1.+delt2*tfrc(jlev))
       zsdt2(:,jlev)=-tdissd(jlev)*sakpp(1:NSPP,jlev)*sdp(:,jlev)       &
     &              /(1.+delt2*tdissd(jlev)*sakpp(1:NSPP,jlev))
       sdt(:,jlev)=zsdt1(:,jlev)+zsdt2(:,jlev)
       zszt1(:,jlev)=-tfrc(jlev)*szp(:,jlev)                            &
     &              /(1.+delt2*tfrc(jlev))
       zszt2(:,jlev)=-tdissz(jlev)*sakpp(1:NSPP,jlev)*szp(:,jlev)       &
     &              /(1.+delt2*tdissz(jlev)*sakpp(1:NSPP,jlev))
       szt(:,jlev)=zszt1(:,jlev)+zszt2(:,jlev)
      enddo
!
!     energy/entropy diagnostics
!
      if(nenergy > 0) then
       allocate(zsttd(NSPP,NLEV))
       do jlev=1,NLEV
        zsttd(:,jlev)=-tdisst(jlev)*sakpp(1:NSPP,jlev)*stp(:,jlev)      &
     &               /(1.+delt2*tdisst(jlev)*sakpp(1:NSPP,jlev))
       enddo
      endif

!     no friction on planetary vorticity
!     no diffusion on plavor (planetary vorticity)

      if (mypid == NROOT) then
         spp(1) = 0.0
         spp(2) = 0.0
         zszt1(3,:)=zszt1(3,:)+plavor*tfrc(:)/(1.+delt2*tfrc(:))
         zszt2(3,:)=zszt2(3,:)+plavor*tdissz(:)*sakpp(3,:)              &
     &                               /(1.+delt2*tdissz(:)*sakpp(3,:))
         szt(3,:)=szt(3,:)+plavor*(tfrc(:)/(1.+delt2*tfrc(:))           &
     &                            +tdissz(:)*sakpp(3,:)                 &
     &                            /(1.+delt2*tdissz(:)*sakpp(3,:)))
      endif
!
!     sponge layer ad the top
!     (relax t(1) to t(2))
!
      if(nsponge > 0) then
       stt(:,1)=dampsp*(stp(:,2)-stp(:,1))+stt(:,1)
      endif
!
!     add temp-tendencies due to momentum diffusion/dissipation
!     (if switched on)
!
      if(ndheat > 0) call mkdheat(zszt1,zszt2,zsdt1,zsdt2)
!
!     diagnostics of efficiency
!
      if(ndheat > 1) then
       call mpgallsp(zstt2,stt,NLEV)
       zstt(:,:)=zstt(:,:)+zstt2(:,:)
       if (nqspec == 1) call sp2fl(zsq,zqgp,NLEV)
       call sp2fl(zst,ztgp,NLEV)
       call sp2fl(zstt,zdtgp,NLEV)
       if (nqspec == 1) then
          call fc2gp(zqgp,NLON,NLPP*NLEV)
       else
          zqgp(:,:) = 0.0
       endif
       call fc2gp(ztgp,NLON,NLPP*NLEV)
       call fc2gp(zdtgp,NLON,NLPP*NLEV)
       zsum1(:)=0.
       do jlev=1,NLEV
        if (nqspec == 1) zqgp(:,jlev)=zqgp(:,jlev)*psurf/dp(:)
        ztgp(:,jlev)=ct*(ztgp(:,jlev)+t0(jlev))
        zdtgp(:,jlev)=ct*ww*zdtgp(:,jlev)
        zsum1(1)=zsum1(1)+SUM(zdtgp(:,jlev)*zgw(:)                      &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &                ,mask=(zdtgp(:,jlev) >= 0.))
        zsum1(2)=zsum1(2)+SUM(zdtgp(:,jlev)*zgw(:)                      &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) < 0.))
        zsum1(3)=zsum1(3)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         &
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) >= 0.))
        zsum1(4)=zsum1(4)+SUM(zdtgp(:,jlev)/ztgp(:,jlev)*zgw(:)         & 
     &                *acpd*(1.+adv*zqgp(:,jlev))*dp(:)/ga*dsigma(jlev) &
     &               ,mask=(zdtgp(:,jlev) < 0.))
       enddo
       zsum3(1)=SUM(zgw(:))
       call mpsumbcr(zsum1,4)
       call mpsumbcr(zsum3,1)
       zsum1(:)=zsum1(:)/zsum3(1)
       if(mypid == NROOT) then
        ztp=zsum1(1)/zsum1(3)
        ztm=zsum1(2)/zsum1(4)
        write(9,*) zsum1(:),zsum1(1)/zsum1(3),zsum1(2)/zsum1(4)         &
     &            ,(ztp-ztm)/ztp
       endif
       deallocate(zst)
       if (nqspec == 1) deallocate(zsq)
       deallocate(zstt)
       deallocate(zstt2)
       deallocate(zdtgp)
       deallocate(zsum1)
       deallocate(zgw)
      endif
!
!     add diabatic tendencies
!
      szp = szp + delt2 * szt
      stp = stp + delt2 * stt
      sdp = sdp + delt2 * sdt
      if (nqspec == 1) sqp = sqp + delt2 * sqt

      if (nenergy > 1) then
         call mpgallsp(dsdiv,sdp,NLEV)
         call convweight(dsdiv,dq,dp,dconvspd(3))
         deallocate(dsdiv)
         if (mypid == NROOT) then
            if (nstep > nstep1 + ntspd) then
               dconvspa(1:3) = dconvspa(1:3) + dconvspd(1:3)
               dconvspa(4) = dconvspa(4) + 1.0
            endif
            if (mod(nstep,ndiag) == 0 .and. dconvspa(4) > 0.0) then
               write(nud,'(A,I9,3E15.6)') ' CONVSPD in mid out ',       &
     &            nstep, dconvspa(1:3) / dconvspa(4)
            endif
         endif
      endif

!     initial divergence damping for a smooth start of the model in case
!     of steep orography (Mars) or unusual initial conditions

      if (ndivdamp > 0) then
         zdd = 1.0 / (1.0 + ndivdamp * 0.02)
         sdp(:,:) = sdp(:,:) * zdd
!        zdd is the same number on every thread, so root's line is the record
!        and the other NPRO-1 copies were only interleaving. world-0ihs.
         if (mypid == NROOT) write(nud,*) '### Damping with ',zdd
         ndivdamp = ndivdamp - 1        
      endif

      if (nkits == 0) then
         szm = szm + pnu * szp
         stm = stm + pnu * stp
         sdm = sdm + pnu * sdp
         spm = spm + pnu * spp
         if (nqspec == 1) sqm = sqm + pnu * sqp
      endif

!     These are the whole of the gather traffic. Under the shared build the
!     partial IS the slice, so mpgathersp is a barrier and moves nothing; under
!     MPI it is the allgather this used to be.
      call mpsyncsp
!
!     franks diagnostic
!

      if(ndiagsp==1) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       dsp3d(:,1:NLEV,2)=ztt(:,1:NLEV)*ct*ww
       deallocate(ztt)
      endif
      if(ndiaggp==1) then
       call mkdqtgp
       do jlev=1,NLEV
        dgp3d(:,jlev,21)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      endif
!
!     energy diagnostics
!
      if(nenergy > 0) then
       allocate(ztt(NESP,NLEV))
       call mpgallsp(ztt,stt,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
#ifdef OMPSHARED
       if (nshtns == 1) then
!$omp barrier
         call sh_sp2gp(ztt, zttgp_g, NLEV)
!$omp barrier
       else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       endif
#else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
#endif
       denergy(:,5)=0.
       do jlev=1,NLEV
        denergy(:,5)=denergy(:,5)                                       &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,5)=zttgp(:,jlev)                                &
     &                  *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
       call mpgallsp(ztt,zsttd,NLEV)
       ztt(:,:)=ztt(:,:)*ct*ww
#ifdef OMPSHARED
       if (nshtns == 1) then
!        THE LEADING BARRIER IS NOT SYMMETRY: the denergy05 loop above reads
!        the same globe this call overwrites, and a thread that arrives early
!        would take the diffusion heating away from a thread still forming 05.
!$omp barrier
         call sh_sp2gp(ztt, zttgp_g, NLEV)
!$omp barrier
       else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
       endif
#else
       call sp2fl(ztt,zttgp,NLEV)
       call fc2gp(zttgp,NLON,NLPP*NLEV)
#endif
       denergy(:,24)=0.
       do jlev=1,NLEV
        denergy(:,24)=denergy(:,24)                                     &
     &              +zttgp(:,jlev)                                      &
     &              *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,24)=zttgp(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif 
       enddo
!      Carry the temperature diffusion's unbalanced heating to the fixer, which
!      runs in the next spectrala and cannot see this routine's terms otherwise.
       if(nenergyfix > 0) then
        zd24(1) = 0.0
        zd24(2) = 0.0
        jhor = 0
        do jlat = 1 , NLPP
         do jlon = 1 , NLON
          jhor = jhor + 1
          zd24(1) = zd24(1) + denergy(jhor,24)*gwd(jlat)
          zd24(2) = zd24(2) + gwd(jlat)
         enddo
        enddo
        call mpsumbcr(zd24,2)
        denergyd24 = zd24(1)/zd24(2)
       endif
       deallocate(ztt)
      endif
      if(nenergy > 0) then
       deallocate(zsttd)
      endif
!
      return
      end subroutine spectrald

!     ===============
!     SUBROUTINE MKDQ
!     ===============

!      subroutine mkdq(psq,pgq,KLEV)
!      use pumamod
!
!      real psq(NESP,klev),pgq(NHOR,klev)
!
!      real zps(NHOR)
!
!      call sp2fl(psq,pgq,KLEV)
!      call fc2gp(pgq,NLON,NLPP*KLEV)
!      call sp2fl(sp,zps,1)
!      call fc2gp(zps,NLON,NLPP)
!
!      zps(:)=psurf*EXP(zps(:))
!      do jlev=1,KLEV
!       pgq(:,jlev)=pgq(:,jlev)*psurf/zps(:)
!      enddo
!
!      return
!      end subroutine mkdq


!     ===================
!     SUBROUTINE MKSECOND
!     ===================

      subroutine mksecond(ps,poff)
!
      real ps,poff
      integer ivalues(8)
!
      call date_and_time(values=ivalues)
!
      zdays = ivalues(3) - 1
      zhour = 24.0 * zdays + ivalues(5)
      zmins = 60.0 * zhour + ivalues(6)
      zsecs = 60.0 * zmins + ivalues(7) + 0.001 * ivalues(8)
      ps = zsecs - poff
!
      return
      end

!     ==================
!     SUBROUTINE MKDHEAT
!     ==================

      subroutine mkdheat(zszt1,zszt2,zsdt1,zsdt2)
      use pumamod
#ifdef OMPSHARED
      use shtnsmod, only: sh_dv2uv, sh_sp2gp, sh_gp2sp, sh_slice
#endif
!
      real zszt1(NSPP,NLEV),zszt2(NSPP,NLEV)
      real zsdt1(NSPP,NLEV),zsdt2(NSPP,NLEV)
!
!     The six FULL spectral arrays this routine used to keep on the stack are
!     zhd, zhz, zhq, zhe and the three partial slots zhf1, zhf2, zhef, all in
!     pumamod. They were 14.1 MB a thread at T170 and one complete copy of the
!     global field per thread, which is the thing Stages A and B removed
!     everywhere else. What is left here is genuinely per-process: the grid
!     fields, which are this process's latitudes, and the NSPP partials.
      real zsdp(NSPP,NLEV),zszp(NSPP,NLEV),zsqp(NSPP,NLEV)
!
      real zsde(NSPP,NLEV)
      real zstt1(NSPP,NLEV)
      real zstt2(NSPP,NLEV)
!
      zsdp(:,:)=sdp(:,:)
      zszp(:,:)=szp(:,:)
      if (nqspec == 1) then
         zsqp(:,:)=sqp(:,:)
         call mpgallspp(zhq,zsqp,NLEV)
#ifdef OMPSHARED
      if (nshtns == 1) then
!        SHTns lands in GRID space, so sp2fl and its fc2gp both go.
!$omp barrier
         call sh_sp2gp(zhq, hdq_g, NLEV)
!$omp barrier
      else
         call sp2fl(zhq,hdq,NLEV)
         call fc2gp(hdq,NLON,NLPP*NLEV)
      endif
#else
         call sp2fl(zhq,hdq,NLEV)
         call fc2gp(hdq,NLON,NLPP*NLEV)
#endif
      else
         hdq(:,:) = 0.0
      endif
      call mpgallspp(zhd,zsdp,NLEV)
      call mpgallspp(zhz,zszp,NLEV)
#ifdef OMPSHARED
      if (nshtns == 1) then
!        One call for the pair, and it lands in grid space.
!$omp barrier
         call sh_dv2uv(zhd, zhz, hdu_g, hdv_g, NLEV)
!$omp barrier
      else
      call dv2uv(zhd,zhz,hdu,hdv)
      call fc2gp(hdu,NLON,NLPP*NLEV)
      call fc2gp(hdv,NLON,NLPP*NLEV)
      endif
#else
      call dv2uv(zhd,zhz,hdu,hdv)
      call fc2gp(hdu,NLON,NLPP*NLEV)
      call fc2gp(hdv,NLON,NLPP*NLEV)
#endif
!
      zsdp(:,:)=sdp(:,:)+zsdt1(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt1(:,:)*delt2
      call mpgallspp(zhd,zsdp,NLEV)
      call mpgallspp(zhz,zszp,NLEV)
#ifdef OMPSHARED
      if (nshtns == 1) then
!        One call for the pair, and it lands in grid space.
!$omp barrier
         call sh_dv2uv(zhd, zhz, hdun_g, hdvn_g, NLEV)
!$omp barrier
      else
      call dv2uv(zhd,zhz,hdun,hdvn)
      call fc2gp(hdun,NLON,NLPP*NLEV)
      call fc2gp(hdvn,NLON,NLPP*NLEV)
      endif
#else
      call dv2uv(zhd,zhz,hdun,hdvn)
      call fc2gp(hdun,NLON,NLPP*NLEV)
      call fc2gp(hdvn,NLON,NLPP*NLEV)
#endif
!
      do jlev = 1 , NLEV
       hdu(:,jlev)=cv*hdu(:,jlev)*SQRT(rcsq(:))
       hdv(:,jlev)=cv*hdv(:,jlev)*SQRT(rcsq(:))
       hdun(:,jlev)=cv*hdun(:,jlev)*SQRT(rcsq(:))
       hdvn(:,jlev)=cv*hdvn(:,jlev)*SQRT(rcsq(:))
       if (nqspec == 1) hdq(:,jlev)=hdq(:,jlev)*psurf/dp(:)

       hddt(:,jlev)=-(hdun(:,jlev)*hdun(:,jlev)                          &
     &                -hdu(:,jlev)*hdu(:,jlev)                            &
     &                +hdvn(:,jlev)*hdvn(:,jlev)                          &
     &                -hdv(:,jlev)*hdv(:,jlev))/deltsec2                  &
     &               *0.5/acpd/(1.+adv*dq(:,jlev))
      enddo
!
      hddt(:,:)=hddt(:,:)/ct/ww
#ifdef OMPSHARED
      if (nshtns == 1) then
!        The finished field, so the reduction has nothing to add and
!        each thread takes its slice instead.
!$omp barrier
         call sh_gp2sp(hddt_g, zhf1(1,1,0), NLEV)
!$omp barrier
         call sh_slice(zhf1(1,1,0), zstt1, NLEV)
!$omp barrier
      else
      call gp2fc(hddt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hddt(:,jlev),zhf1(1,jlev,mypart))
      enddo
      call mpsumscp(zhf1,zstt1,NLEV)
      endif
#else
      call gp2fc(hddt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hddt(:,jlev),zhf1(1,jlev,mypart))
      enddo
      call mpsumscp(zhf1,zstt1,NLEV)
#endif
!
      zsdp(:,:)=sdp(:,:)+zsdt2(:,:)*delt2
      zszp(:,:)=szp(:,:)+zszt2(:,:)*delt2
      call mpgallspp(zhd,zsdp,NLEV)
      call mpgallspp(zhz,zszp,NLEV)
#ifdef OMPSHARED
      if (nshtns == 1) then
!        One call for the pair, and it lands in grid space.
!$omp barrier
         call sh_dv2uv(zhd, zhz, hdun_g, hdvn_g, NLEV)
!$omp barrier
      else
      call dv2uv(zhd,zhz,hdun,hdvn)
      call fc2gp(hdun,NLON,NLPP*NLEV)
      call fc2gp(hdvn,NLON,NLPP*NLEV)
      endif
#else
      call dv2uv(zhd,zhz,hdun,hdvn)
      call fc2gp(hdun,NLON,NLPP*NLEV)
      call fc2gp(hdvn,NLON,NLPP*NLEV)
#endif
!
      do jlev = 1 , NLEV
       hdun(:,jlev)=cv*hdun(:,jlev)*SQRT(rcsq(:))
       hdvn(:,jlev)=cv*hdvn(:,jlev)*SQRT(rcsq(:))
       hdek(:,jlev)=(hdun(:,jlev)*hdun(:,jlev)                          &
     &                -hdu(:,jlev)*hdu(:,jlev)                            &
     &                +hdvn(:,jlev)*hdvn(:,jlev)                          &
     &                -hdv(:,jlev)*hdv(:,jlev))/deltsec2                  &
     &               *dp(:)/ga*dsigma(jlev)   
      enddo
#ifdef OMPSHARED
      if (nshtns == 1) then
!        The finished field, so the reduction has nothing to add and
!        each thread takes its slice instead.
!$omp barrier
         call sh_gp2sp(hdek_g, zhef(1,1,0), NLEV)
!$omp barrier
         call sh_slice(zhef(1,1,0), zsde, NLEV)
!$omp barrier
      else
      call gp2fc(hdek,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hdek(:,jlev),zhef(1,jlev,mypart))
      enddo
      call mpsumscp(zhef,zsde,NLEV)
      endif
#else
      call gp2fc(hdek,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hdek(:,jlev),zhef(1,jlev,mypart))
      enddo
      call mpsumscp(zhef,zsde,NLEV)
#endif
      call mpgallspp(zhe,zsde,NLEV)
!     Only the global mean survives, and zhe is one array the whole team reads,
!     so the clear goes slice by slice rather than whole-array from every
!     thread. mpzerosp is that, and on the MPI build it is the plain statement
!     it replaces.
      call mpzerosp(zhe,2,NLEV)
#ifdef OMPSHARED
      if (nshtns == 1) then
!$omp barrier
         call sh_sp2gp(zhe, hdek_g, NLEV)
!$omp barrier
      else
      call sp2fl(zhe,hdek,NLEV)
      call fc2gp(hdek,NLON,NLPP*NLEV)
      endif
#else
      call sp2fl(zhe,hdek,NLEV)
      call fc2gp(hdek,NLON,NLPP*NLEV)
#endif
      do jlev=1,NLEV
       hddt(:,jlev)=-hdek(:,jlev)                                    &
     &           *0.5/acpd/(1.+adv*dq(:,jlev))/dp(:)*ga/dsigma(jlev) 
      enddo
      hddt(:,:)=hddt(:,:)/ct/ww
#ifdef OMPSHARED
      if (nshtns == 1) then
!        The finished field, so the reduction has nothing to add and
!        each thread takes its slice instead.
!$omp barrier
         call sh_gp2sp(hddt_g, zhf2(1,1,0), NLEV)
!$omp barrier
         call sh_slice(zhf2(1,1,0), zstt2, NLEV)
!$omp barrier
      else
      call gp2fc(hddt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hddt(:,jlev),zhf2(1,jlev,mypart))
      enddo
      call mpsumscp(zhf2,zstt2,NLEV)
      endif
#else
      call gp2fc(hddt,NLON,NLPP*NLEV)
      do jlev=1,NLEV
       call fc2sp(hddt(:,jlev),zhf2(1,jlev,mypart))
      enddo
      call mpsumscp(zhf2,zstt2,NLEV)
#endif
!
!     energy diagnostics
!
      if(nenergy > 0) then
!     NOT OFF BY DEFAULT. config/planet.yaml declares energy_diagnostics and
!     every production namelist records NENERGY = 1, so this block runs on every
!     timestep of every run. It said "off by default" and its two transforms were
!     the last hot legmod callers in the model, outside every nshtns branch, so
!     the SHTns build did the Legendre transform twice: once through the wrappers
!     for the dynamics and once through legmod for this diagnostic. world-dr8.
!
!     zhd is finished with by here and is the right shape, so the diagnostic
!     borrows it rather than keeping a seventh full array alive all timestep.
!     The scaling goes on the PARTIAL before the gather, not on the gathered
!     array afterwards: zhd is shared, and every thread scaling the whole of it
!     would be a race on identical values. zsde is dead from the reduction above
!     and is the right shape.
       zsde(:,:)=zstt1(:,:)*ct*ww
       call mpgallspp(zhd,zsde,NLEV)
#ifdef OMPSHARED
       if (nshtns == 1) then
!        Same contract as the zhe transform above: the wrapper has no trailing
!        barrier, and hddt_g is what the whole team writes while hddt is this
!        thread's band of it. The LEADING barrier matters as much: hddt still
!        holds the heating rate the loop below this block's twin has been
!        reading.
!$omp barrier
          call sh_sp2gp(zhd, hddt_g, NLEV)
!$omp barrier
       else
       call sp2fl(zhd,hddt,NLEV)
       call fc2gp(hddt,NLON,NLPP*NLEV)
       endif
#else
       call sp2fl(zhd,hddt,NLEV)
       call fc2gp(hddt,NLON,NLPP*NLEV)
#endif
       denergy(:,23)=0.
       do jlev=1,NLEV
        denergy(:,23)=denergy(:,23)                                     &
     &               +hddt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,23)=hddt(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo
       zsde(:,:)=zstt2(:,:)*ct*ww
       call mpgallspp(zhd,zsde,NLEV)
#ifdef OMPSHARED
       if (nshtns == 1) then
!$omp barrier
          call sh_sp2gp(zhd, hddt_g, NLEV)
!$omp barrier
       else
       call sp2fl(zhd,hddt,NLEV)
       call fc2gp(hddt,NLON,NLPP*NLEV)
       endif
#else
       call sp2fl(zhd,hddt,NLEV)
       call fc2gp(hddt,NLON,NLPP*NLEV)
#endif
       denergy(:,25)=0.
       do jlev=1,NLEV
        denergy(:,25)=denergy(:,25)                                     &
     &               +hddt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,25)=hddt(:,jlev)                               &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        endif
       enddo

      endif
!
      stt(:,:)=stt(:,:)+zstt1(:,:)+zstt2(:,:)
!
      return
      end subroutine mkdheat
