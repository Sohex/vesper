#define kasting 0
#define realv 1
#define keeprecord 0
#define outputrain 0
#define precipweath 1

      module carbonmod
      use pumamod
!
!     version identifier (date)
!
      character(len=80) :: version = '07.14.2017 by Adiv'
!
!     Parameter   
!                          !Outgassing rates given in ubars/year
! #if realv == 1
!       parameter(VEARTH=5.0e-2) ! Rate based on volcanic outgassing estimates from Gerlach 2011
! #else
! #if kasting==1
!       parameter(VEARTH=6.3e-2)
! #else
!       parameter(VEARTH=7.0e-3) ! Earth volcanic CO2 outgassing rate (ubars), assumed to be equal to Earth
!                                ! weathering rate.
! #endif
! #endif
!     THIS WHOLE MODULE IS FITTED TO EARTH, and it is dormant at NCARBON = 0,
!     which is what every run of this project sets. `carbonstep` is still called
!     every timestep; everything below the ncarbon test is what a nonzero
!     NCARBON turns on. Recorded here because a dormant Earth constant reads
!     exactly like a live one to whoever throws the switch. world-9d1.
!
!       CO2EARTH, PEARTH and VEARTH are Earth reference states the weathering
!         law is expressed as a ratio against. They are the scheme's normalisation
!         and are not wrong on their own; what they mean for a different
!         atmosphere is undecided.
!       RAD_EARTH and RAD_EARTHSQ sit beside a correct `plarad` in pumamod and
!         are used for the surface area the outgassing is spread over. On this
!         world they are wrong by the radius ratio squared.
!       exp(kact*(tsurf - 288.0)) below takes 288 K as the reference surface
!         temperature. That is Earth's global mean, not this world's.
!       tune1 and tune2 are commented as a tuning adjustment to make the global
!         average match, so they are fitted to Earth by construction.
!
!     Turning NCARBON on is a re-derivation of all of these, not a namelist
!     change. dust-13's shape, for the carbon cycle.
      parameter(CO2EARTH=330.0)  ! Earth CO2 level in ubars
      parameter(RAD_EARTH=6371220.0) !Earth radius
      parameter(RAD_EARTHSQ=6371220.0*6371220.0) !Earth radius squared
!
!     namelist parameters
!

      character (256) :: carbon_namelist     = "carbonmod_namelist"
      
      integer :: ncarbon = 1 ! 1 = do weathering, 0 = constant pCO2
      integer :: nsupply = 0 ! Should it be supply-limited weathering?
      integer :: nco2evolve = 0 !Should CO2 and pressure evolve each year? (0/1)
      real :: volcanco2 = 0.02 ! Volcanic outgassing rate in units of VEARTH
      real :: kact = 0.09 ! activation energy factor
      real :: krun = 0.045 ! runoff efficiency factor
      real :: beta = 0.5 ! pCO2 dependence
      real :: frequency = 4.0 ! Number of times to compute weathering per day. Can be a float, i.e.
                              ! for once every 2 days, use frequency=0.5.
      real :: VEARTH = 5.0e-2 !ubars/year
      real :: PEARTH = 79.0 ! Annual precipitation on modern Earth that is relevant for weathering.
                           ! 79 cm/yr (Chen 2002 & Schneider 2014)
      real :: WMAX = 1.0 ! Maximum weathering rate for the supply-limited case, in ubar/yr
      
      real :: psurf0 = 101100.0 !Mean sea-level pressure
!
!     global arrays
!
      real :: dglobe(NHOR) = 0. ! Normalized cell area

      real :: localweathering(NHOR) = 0.0 ! Local instantaneous weathering rate relative to Earth standard
      real :: localavgweather(NHOR) = 0.0 ! Annual average weathering rate relative to Earth standard
      real :: localavgtemps(NHOR) = 0.0 ! Annual average surface temperature
      real :: localprecip(NHOR) = 0.0 !Local precipitation
      real :: localavgprecip(NHOR) = 0.0!Local net precipitation
      
      real :: aweathering(NHOR) = 0.0 !Monthly average weathering
      
#if outputrain==1      
      real :: cprecip(NHOR) = 0.0
#endif

!
!     global scalars
!
      integer :: interval = 8   ! Number of timesteps to go between weathering updates
      integer :: cstep = 0 ! Current timestep number in the cycle. When cstep = interval, do weathering
      integer :: istep = 0
!       integer :: nco2evolve = 0 ! Do we allow CO2 to evolve year-to-year?
      real :: timeweight = 0.0 ! Weight for computing annual averages
      real :: avgweathering = 0.0 ! Global annual average weathering rate
      real :: dpco2dt = 0.0 ! Change in pCO2 with respect to time. = (volcanco2 - avgweathering)*VEARTH
      real :: tune1 = 5.41 ! Tuning adjustment to make global average match for non-precip model.
      real :: tune2 = 2.20 ! Tuning adjustment to make precip model match non-precip model.
!

!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(avgweathering,aweathering,beta,carbon_namelist,cstep,dglobe,dpco2dt,frequency,&
!$omp&  interval,istep,kact,krun,localavgprecip,localavgtemps,localavgweather,localprecip,&
!$omp&  localweathering,ncarbon,nco2evolve,nsupply,pearth,psurf0,timeweight,tune1,tune2,vearth,&
!$omp&  version,volcanco2,wmax)

      end module carbonmod
      
!      
!     =============================
!

      subroutine carbonini
      use carbonmod
      
      namelist/carbonmod_nl/ncarbon,volcanco2,kact,krun,beta,frequency,VEARTH,PEARTH,nsupply,WMAX,nco2evolve
      
      if (mypid==NROOT) then
         open(23,file=carbon_namelist)
         read(23,carbonmod_nl)
         close(23)
         write(nud,'(/," *********************************************")')
         write(nud,'(" * CARBONMOD ",a34," *")') trim(version)
         write(nud,'(" *********************************************")')
         write(nud,'(" * Namelist CARBON_NL from <carbon_namelist> *")')
         write(nud,'(" *********************************************")')
         write(nud,carbonmod_nl)
      endif
      
      frequency = min(frequency,real(mtspd))
      interval = int(mtspd/frequency)
      timeweight = (1.0 / (frequency*m_days_per_year))
      
      psurf0 = psurf
      
      call makeareas
      
      call mpbci(interval)
      call mpbci(ncarbon)
      call mpbcr(volcanco2)
      call mpbcr(kact)
      call mpbcr(krun)
      call mpbcr(beta)
      call mpbcr(frequency)
      call mpbcr(timeweight)
      call mpbcr(PEARTH)
      call mpbcr(VEARTH)
      call mpbcr(WMAX)
      call mpbci(nsupply)
      call mpbci(nco2evolve)
      call mpbcr(psurf0)
      
      end subroutine carbonini
      
!      
!     =============================
!
      subroutine carbonstep
      use carbonmod
      use rainmod
      
      real tsurf(NHOR)
      real pco2
      real cweathering
      real landfraction
      real globalavgt

#if outputrain==1      
      character srainstep*5   
      real :: fullprecip(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(fullprecip)
      real :: fullweather(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(fullweather)
      real :: fullevap(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(fullevap)
      real :: fullqvi(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(fullqvi)
#endif   
      
      integer i
      
      if (ntime > 0) call mksecond(zsec,0.)
      
      if (mypid == NROOT) cstep = cstep + 1
      if (mypid == NROOT) istep = istep + 1
      
      call mpbci(cstep)
      
!     THE MODEL'S YEAR, not Earth's. 3.154e9 is Earth's seconds per year times
!     100 to reach cm, and it stood two lines from a timeweight that correctly
!     uses m_days_per_year, so the two disagreed by the ratio of the years --
!     about a factor of two on this world. dprc and dprl are in m/s, so the
!     conversion is the orbit in seconds times 100. world-9d1.
      localprecip(:) = (dprc(:) + dprl(:))*m_days_per_year*day_24hr*100.0 !cm/yr
      localavgprecip(:) = localavgprecip(:) + localprecip(:)*timeweight
      
      
      if (cstep .eq. interval) then
      
             
         if (mypid == NROOT) cstep = 0
         
         tsurf(:) = dt(:,NLEP)
         
         if (ncarbon .gt. 0.5) then !we're doing this
         
         localweathering(:) = 0.
         
         do i=1,NHOR
            if (dls(i) .gt. 0.5) then !land
               if (tsurf(i) .ge. 273.15) then !not frozen
                  pco2 = co2*dp(i)*1e-5 ! Convert ppmv to ubars (1e-6 for ppmv->ppv, and 10 for Pa->ubars)
#if precipweath == 1
                  localweathering(i) = ((pco2/CO2EARTH)**beta) * exp(kact*(tsurf(i)-288.0)) * &
     &                                 (localprecip(i)/PEARTH)**0.65 !* tune1*tune2 
#else
                  localweathering(i) = ((pco2/CO2EARTH)**beta) * exp(kact*(tsurf(i)-288.0)) * &
                                       (1+krun*(tsurf(i)-288.0))**0.65 * tune1
#endif
                  if (nsupply .eq. 1) then !Apply a weathering supply limit following Foley 2015
                     wmaxp = max((WMAX/VEARTH),2.0e-15)
                     localweathering(i) = wmaxp*(1-exp(-localweathering(i)/wmaxp))
                  endif
                  
               endif
            endif
#if outputrain==1
            cprecip(i) = (dprc(i) + dprl(i))
#endif
         enddo
         
         
         localweathering(:) = localweathering(:)*plarad**2 / (RAD_EARTHSQ*0.29) !Account for land fraction
         
!          endif
         
#if outputrain==1     
         call mpgagp(fullprecip,cprecip,1)
         call mpgagp(fullweather,localweathering*timeweight,1)
         call mpgagp(fullevap,devap,1)
         call mpgagp(fullqvi,dqvi,1)
         if (mypid == NROOT) then
            write(srainstep,'(I5.5)') istep
            call writegtextarray(fullprecip,NUGP,'precipmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullweather,NUGP,'weathrmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullevap,NUGP,'evaprtmp'//trim(adjustl(srainstep)))
            call writegtextarray(fullqvi,NUGP,'humdtymp'//trim(adjustl(srainstep)))
         endif
#endif        
        
         endif
         
         
#if keeprecord==1
         call write_short(localweathering,cweathering)
         call write_short(dls,landfraction)
         if (mypid==NROOT) then
            write(nud,'(" Land Fraction        = ",f10.2," [...]")') landfraction
            write(nud,'(" Average Weathering   = ",f10.2," [earth]")') cweathering
            cweathering = cweathering / (landfraction * PEARTH)
            open(unit=44,file="weathertrack.pso",position="append",status="unknown")
            write(44,'(1p8e13.5)') cweathering*ncarbon*VEARTH
            close(44)
         endif
#endif
         
          
         localavgtemps(:) = localavgtemps(:) + tsurf(:)*timeweight
         localavgweather(:) = localavgweather(:) + localweathering(:)*timeweight
         call mpbci(cstep)
         
         
      endif
      
      if (ntime > 0) then
         call mksecond(zsec,zsec)
         time4co2 = time4co2 + zsec
!        THE ACCUMULATION IS PER THREAD, THE REPORT IS ROOT'S. fluxstop,
!        rainstop, radstop and miscstop all print their time4* under
!        `mypid == NROOT .and. ntime == 1`; this block was the one that
!        printed from every thread, and nud is one Fortran unit for the
!        whole team. world-0ihs.
         if (mypid == NROOT) then
            write(nud,*)
            write(nud,*)'Carbon cycle routines (ROOT thread only): ',zsec,'s'
         endif
      endif
      
      return
      end subroutine carbonstep
      
!     
!     ============================
!

      subroutine carbonstop
      use carbonmod
      
      real landfraction
      real newco2
      real pco2
      real newpco2
      real globalavgt
      real :: globalweath(NUGP) = 0.0
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(globalweath)
      
      
      call mpgagp(globalweath,localavgweather,1)
      
      call write_short(localavgweather,avgweathering)
!       avgweathering = avgweathering / PEARTH
      call write_short(localavgtemps,globalavgt)
      call write_short(dls,landfraction)
      
      
      !
      ! We only computed weathering on land, but our global average included the ocean, so if we 
      ! used Earth conditions, we would so far have a weathering rate of W_earth * (land fraction). 
      ! Therefore, we have to divide by the land fraction. If we were to include seafloor weathering,
      ! then we would divide by the present-day (landfraction+seafraction), where the fraction that is
      ! sea ice is excluded. This way, the advance of sea ice would actually reduce the global weathering
      ! regardless of temperature--the fraction of the Earth undergoing weathering would go down from 
      ! present-day.
      !
      !
      ! Nevermind. We seem to get the right answer without dividing by the land fraction...
      ! Maybe the equation is already parameterized to account for land fraction...?
      !
      ! From Abbot, et al.... this isn't the global weathering rate; it's the continental weathering rate.
      ! This can be modulated by reduction in land due to glacial advance or sea level rise/fall, but the
      ! overall weathering rate is this PLUS the seafloor weathering rate, modulated by changes in sea surface
      ! area.
      !
      ! No, we actually should divide by the land fraction. Because we're summing area-weighted and time-weighted 
      ! estimates of the global weathering rate due to continental weathering--each cell's value before area-weighting
      ! is an estimate of what the global weathering rate would be if every cell on the planet had that weathering
      ! rate. So when we do the global average, if we don't divide by the land fraction, then we're actually giving
      ! the continental weathering rate multiplied by the land fraction. We have to divide to get the global overall
      ! weathering rate (factor of 3.33 increase, seems okay based on 1 AU tests) due to continental weathering. If 
      ! we were to add in seafloor weathering, then we'd add factors to account for changing land and sea area, but
      ! there would also be weights to determine the relative contributions of each, and those are parameterized (poorly).
!       if (mypid==NROOT) then 
!          if (landfraction > 0) then
!            avgweathering = avgweathering / landfraction
!          else
!            avgweathering = 0.0
!          endif
!       endif

      ! We should divide by Earth's land fraction as a constant pre-factor, and then upon integrating
      ! over the whole planet, we'll pick up a multiplicative factor of the actual land fraction.
      
      if (mypid==NROOT) then
         
         globalweath(:) = globalweath(:)*VEARTH*1000.0
         call writegtextarray(globalweath,NUGP,'annualweather')
         dpco2dt = ncarbon*VEARTH*(volcanco2 - avgweathering)
         pco2 = co2*1e-6*psurf0 ! Pa
         newco2 = (pco2 + dpco2dt*0.1)/(psurf0+dpco2dt*0.1) ! ppv [Pa/Pa]
         newpco2 = newco2*(psurf0+dpco2dt*0.1)*1.e-5 ! Convert to bars (from Pa)
         open(unit=43,file="weathering.pso",position="append",status="unknown")
         write(43,'(1p8e13.5)') pco2*1e-5,globalavgt,avgweathering,volcanco2,dpco2dt*1e-6,newpco2
         close(43)
         co2 = newco2*1e6 ! ppmv
         psurf = psurf0 + dpco2dt*0.1 ! Pa
      endif
      call mpbcr(avgweathering)
      call mpbcr(dpco2dt)
      call mpbcr(co2)
      call mpbcr(psurf)
      
      if (nco2evolve > 0.5) then
        call co2update !Change CO2 for following year
        call psurfupdate !Change surface pressure for following year
      endif
      
      return
      end subroutine carbonstop
      
!
!     ===============================
!
      
      subroutine co2update
      use radmod

!     Persist the evolved CO2 into <radmod_namelist>, editing only the key this
!     routine owns. Declaring a private copy of radmod_nl here and WRITING it
!     would emit exactly the keys the copy names and silently delete every other
!     key the file carries -- the stellar spectrum, the band albedo switch, the
!     ozone and water-vapour weights, the trace-gas bands. cons-15.

      character (len=32) :: yval

      if (mypid == NROOT) then
         write(yval,'(1pe24.16)') co2
         call nlsetkey(radmod_namelist,'CO2',adjustl(yval))
      endif
      return
      end subroutine co2update

!
!     ===============================
!

      subroutine psurfupdate
      use pumamod

!     Persist the evolved surface pressure into <plasim_namelist>, editing only
!     the key this routine owns. See co2update. cons-15.

      character (len=32) :: yval

      if (mypid == NROOT) then
         write(yval,'(1pe24.16)') psurf
         call nlsetkey(plasim_namelist,'PSURF',adjustl(yval))
      endif
      return
      end subroutine psurfupdate

!
!     ===============================
!

      subroutine nlsetkey(cfile,ckey,cvalue)
      use pumamod, only: nud

!     Assign <cvalue> to <ckey> in the namelist file <cfile>, leaving every other
!     line of the file byte for byte as it was. An existing assignment of <ckey>
!     is replaced in place; if the file does not assign it, the assignment is
!     inserted ahead of the group terminator. The file is rewritten only once the
!     edit has succeeded, so a missing file, or one with neither the key nor a
!     terminator, leaves the original standing and says so.
!
!     ROOT THREAD ONLY: <nud> is one Fortran unit for the whole thread team, and
!     the working file name is shared. Both callers guard on mypid == NROOT.

      implicit none
      character (len=*), intent(in) :: cfile,ckey,cvalue
      character (len=256) :: yline,yadj,ytmp
      character (len=80)  :: ykey,ywant
      integer :: ios,iunit,junit
      logical :: lset

      ytmp = trim(cfile)//'.new'
      ywant = ckey
      call nlupper(ywant)

      open(newunit=iunit,file=cfile,status='old',iostat=ios)
      if (ios /= 0) then
         write(nud,*) '*** nlsetkey: no file ',trim(cfile),             &
     &                ', ',trim(ckey),' not written'
         return
      endif
      open(newunit=junit,file=ytmp,status='replace')

      lset = .false.
      do
         read(iunit,'(a)',iostat=ios) yline
         if (ios /= 0) exit
         call nlkey(yline,ykey)
         if (ykey == ywant) then
            if (.not. lset) then
               write(junit,'(1x,a," = ",a)') trim(ckey),trim(cvalue)
               lset = .true.
            endif
            cycle                            ! and drop any later duplicate
         endif
         yadj = adjustl(yline)
         if (.not. lset .and. yadj(1:1) == '/') then
            write(junit,'(1x,a," = ",a)') trim(ckey),trim(cvalue)
            lset = .true.
         endif
         write(junit,'(a)') trim(yline)
      enddo
      close(iunit)
      close(junit)

      if (.not. lset) then
         open(newunit=iunit,file=ytmp,status='old',iostat=ios)
         if (ios == 0) close(iunit,status='delete')
         write(nud,*) '*** nlsetkey: ',trim(cfile),' carries neither ',  &
     &                trim(ckey),' nor a group terminator; left unchanged'
         return
      endif

!     the edit stands: put the working file in place

      open(newunit=iunit,file=ytmp,status='old',iostat=ios)
      if (ios /= 0) return
      open(newunit=junit,file=cfile,status='replace')
      do
         read(iunit,'(a)',iostat=ios) yline
         if (ios /= 0) exit
         write(junit,'(a)') trim(yline)
      enddo
      close(junit)
      close(iunit,status='delete')

      return
      end subroutine nlsetkey

!
!     ===============================
!

      subroutine nlkey(cline,ckey)

!     The namelist key assigned on one line of a namelist file, upper cased and
!     stripped of any array subscript. Blank if the line assigns nothing.

      implicit none
      character (len=*), intent(in)  :: cline
      character (len=*), intent(out) :: ckey
      integer :: ieq,ipar

      ckey = ' '
      ieq = index(cline,'=')
      if (ieq < 2) return
      ckey = adjustl(cline(1:ieq-1))
      ipar = index(ckey,'(')
      if (ipar > 1) ckey = ckey(1:ipar-1)
      if (ckey(1:1) == '&' .or. ckey(1:1) == '!' .or. ckey(1:1) == '/') then
         ckey = ' '
         return
      endif
      call nlupper(ckey)

      return
      end subroutine nlkey

!
!     ===============================
!

      subroutine nlupper(ckey)

!     Upper case <ckey> in place, ASCII only. Namelist keys are matched on it.

      implicit none
      character (len=*), intent(inout) :: ckey
      integer :: j,ic

      do j = 1 , len_trim(ckey)
         ic = ichar(ckey(j:j))
         if (ic >= 97 .and. ic <= 122) ckey(j:j) = char(ic-32)
      enddo

      return
      end subroutine nlupper
