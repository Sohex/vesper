!     =============
!     MODULE CALMOD
!     =============

      module calmod
      implicit none
      integer :: mondays(0:12) = (/0,31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: mona365(0:12) = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: monaccu(0:12)
      integer :: ny400d = 400 * 365 + 97
      integer :: ny100d = 100 * 365 + 24
      integer :: ny004d =   4 * 365 +  1
      integer :: ny001d =       365
      integer :: nud        =  6
      real    :: day_24hr         = 86400.0 ! [sec]

!     These values are copied from pumamod in subroutine calini

      integer :: n_days_per_month =  30
      integer :: n_days_per_year  = 360
      integer :: m_days_per_year  = 360
      integer :: m_days_per_month =  30
      integer :: mtspd            =   0
      integer :: n_start_step     =   0
      integer :: ntspd            =   0
      real    :: solar_day        = 86400.0 ! [sec]
      real    :: mpstep           =  45.0 ! Minutes per timestep
      real    :: tcalday          = 86400.0 ! seconds per 'calendar' (1/360th year) day


!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(day_24hr,m_days_per_month,m_days_per_year,mona365,monaccu,mondays,mpstep,mtspd,&
!$omp&  n_days_per_month,n_days_per_year,n_start_step,ntspd,nud,ny001d,ny004d,ny100d,ny400d,&
!$omp&  solar_day,tcalday)

      end module calmod


!     =================
!     SUBROUTINE CALINI
!     =================

      subroutine calini(k_days_per_month,k_days_per_year,k_start_step &
                       ,ktspd,psolday,kpid,kmpstep,kcal_days_per_year &
                       ,km_days_per_year,km_days_per_month,kmtspd)
      use calmod

!     calini has no implicit none of its own, and kmpstep would otherwise be
!     typed INTEGER by its initial letter while the caller passes REAL mpstep.
      real    kmpstep
      integer kcal_days_per_year
      integer km_days_per_year
      integer km_days_per_month
      integer kmtspd
      
      n_days_per_month = k_days_per_month
      n_days_per_year  = k_days_per_year
      n_start_step     = k_start_step
      ntspd            = ktspd
      solar_day        = psolday

!     ALL FIVE COPIED, and mtspd copied rather than derived. calmod declares its
!     own m_days_per_year, m_days_per_month and mtspd beside pumamod's, with the
!     same names, and this routine used to copy the first two of the five and
!     leave the calendar's own trio at Earth's 360 and 30 with an mtspd it
!     computed from them. The model then encoded a date with pumamod's numbers
!     and decoded it with Earth's, and cal2step and step2cal30 were not
!     inverses. world-x1k.
!
!     THE CALENDAR DAY IS THE 24-HOUR DAY. m_days_per_year is the count of
!     24-hour days in an orbit and mtspd is timesteps per 24-hour day, so
!     m_days_per_year * mtspd is the orbit in timesteps and tcalday below comes
!     out at day_24hr by construction rather than by coincidence.
      m_days_per_year  = km_days_per_year
      m_days_per_month = km_days_per_month
      mtspd            = kmtspd
      
      mpstep = kmpstep
      
      tcalday = mtspd * mpstep * 60.0 ! Timesteps/calendar day * minutes/timestep * 60 s/min

!     THE GREGORIAN TRAPDOOR, closed. n_days_per_year == 365 switches this
!     module to Earth's calendar entire -- the 400/100/4 leap rule, the twelve
!     named months of unequal length, weekdays -- at yday2mmdd, cal2step,
!     step2cal, ndayofyear and nweekday. It is derived from this world's flux
!     and rotation period, so nothing structurally keeps it off 365, and there
!     was no guard and no message. A world that lands there gets an abort and
!     the reason, not a different planet's calendar.
!     SERIALISED RATHER THAN ROOT-GUARDED. nud is unit 6, one Fortran unit for
!     the whole thread team, and calini is reached by every thread from prolog.
!     calmod is a leaf module with no pumamod dependency and so has no thread
!     identity to test; a critical section keeps the seven lines from
!     interleaving without giving the calendar a dependency on the parallel
!     layer. world-0ihs.
      if (n_days_per_year == 365) then
!$omp critical (nudwrite)
         write(nud,*) '*** calini: n_days_per_year is 365 ***'
         write(nud,*) 'That value switches this model to Earth''s Gregorian'
         write(nud,*) 'calendar: the 400/100/4 leap rule, twelve named months'
         write(nud,*) 'of unequal length, and weekdays. It is a rotation count'
         write(nud,*) 'for THIS world and the coincidence is meaningless.'
         write(nud,*) 'Change the rotation period or the orbit rather than'
         write(nud,*) 'running on a calendar that is not this planet''s.'
!$omp end critical (nudwrite)
         stop
      endif

      kcal_days_per_year = m_days_per_year

      if (kpid == 0) then
         write(nud,1050)
         write(nud,1060)
         write(nud,1050)
         write(nud,1000) n_days_per_year
         write(nud,1010) n_days_per_month
         write(nud,1030) n_start_step
         write(nud,1040) ntspd
         write(nud,1045) mtspd
         write(nud,1046) nint(tcalday)
         write(nud,1050)
         endif
      return
 1060 format(" * Calmod initialization             *")
 1000 format(" * Days per year:    ",i14," *")
 1010 format(" * Days per month:   ",i14," *")
 1030 format(" * Start step:",i21," *")
 1040 format(" * Timesteps per day:",i14," *")
 1045 format(" * Timesteps per calendar day:",i5," *")
 1046 format(" * Sec per calendar day:",i11," *")
 1050 format(" *************************************")

      end subroutine calini

!     ====================
!     SUBROUTINE YDAY2MMDD
!     ====================

      subroutine yday2mmdd(kyday,kmon,kday)
      use calmod
      if (n_days_per_year == 365) then
         kmon = 1
         do while (kyday > mona365(kmon))
            kmon = kmon + 1
         enddo
         kday = kyday - monaccu(kmon-1)
      else
         kmon = (kyday-1) / n_days_per_month
         kday = kyday - n_days_per_month * kmon
      endif
      return
      end

!     =================
!     FUNCTION NWEEKDAY
!     =================

      integer function nweekday(kday)
      nweekday = mod(kday+5,7)
      return
      end


!     ===================
!     SUBROUTINE STEP2CAL
!     ===================

      subroutine step2cal(kstep,ktspd,kdatim)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep     ! time step since simulation start
      integer, intent(IN ) :: ktspd     ! time steps per day
      integer, intent(OUT) :: kdatim(7) ! year,month,day,hour,min,weekday,leapyear

      integer :: iyea  ! current year   of simulation
      integer :: imon  ! current month  of simulation
      integer :: iday  ! current day    of simulation
      integer :: ihou  ! current hour   of simulation
      integer :: imin  ! current minute of simulation
      integer :: idall
      integer :: istp
      integer :: iy400,id400
      integer :: iy100,id100
      integer :: iy004,id004
      integer :: iy001,id001
      integer :: jmon

      logical :: leap

      idall = kstep   / ktspd

      iy400 = idall   / ny400d          ! segment [of 400 years]
      id400 = mod(idall,ny400d)

      if (id400 <= ny100d) then         ! century year is leap year
         iy100 =         0              ! century in segment [0]
         id100 =     id400
         iy004 =     id100 / ny004d     ! tetrade in century [0..24]
         id004 = mod(id100 , ny004d)
         leap  = (id004 <= ny001d)
         if (leap) then
            iy001 =     0               ! year in tetrade [0]
            id001 = id004
         else
            iy001 =    (id004-1)/ny001d ! year in tetrade [1,2,3]
            id001 = mod(id004-1, ny001d)
         endif
      else                              ! century year is not leap year
         iy100 =    (id400-1)/ny100d    ! century in segment [1,2,3]
         id100 = mod(id400-1, ny100d)
         if (id100 < ny004d-1) then
            iy004 = 0                   ! tetrade in century [0]
            id004 = id100
            leap  = .false.
            iy001 =     id004/ny001d    ! year in tetrade [1,2,3]
            id001 = mod(id004,ny001d)
         else
            iy004 = (id100+1)/ny004d    ! tetrade in century [1..24]
            id004 = mod(id100+1,ny004d)
            leap  = (id004 <= ny001d)
            if (leap) then
               iy001 =     0            ! year in tetrade [0]
               id001 = id004
            else
               iy001 =    (id004-1)/ny001d
               id001 = mod(id004-1, ny001d)
            endif
         endif
      endif

      iyea  = iy400 * 400 + iy100 * 100 + iy004 * 4 + iy001

      monaccu(0) = mondays(0)
      monaccu(1) = mondays(1)
      monaccu(2) = mondays(1) + mondays(2)
      if (leap) monaccu(2) = monaccu(2) + 1
      do jmon = 3 , 12
         monaccu(jmon) = monaccu(jmon-1) + mondays(jmon)
      enddo
      imon = 1
      id001 = id001 + 1
      do while (id001 > monaccu(imon))
         imon = imon + 1
      enddo
      iday = id001 - monaccu(imon-1)

      istp = mod(kstep,ktspd)
      imin = (istp * 1440) / ktspd
      ihou = imin / 60
      imin = mod(imin,60)

      kdatim(1) = iyea
      kdatim(2) = imon
      kdatim(3) = iday
      kdatim(4) = ihou
      kdatim(5) = imin
      kdatim(6) = mod(kstep/ktspd+5,7) ! day of week
      if (leap) then
         kdatim(7) = 1
      else
         kdatim(7) = 0
      endif

      return
      end subroutine step2cal

!     ===================
!     SUBROUTINE CAL2STEP
!     ===================

      subroutine cal2step(kstep,ktspd,kyea,kmon,kday,khou,kmin)
      use calmod
      implicit none
      integer, intent(OUT) :: kstep ! time step since simulation start
      integer, intent(IN ) :: ktspd ! time steps per day
      integer, intent(IN ) :: kyea  ! current year   of simulation
      integer, intent(IN ) :: kmon  ! current month  of simulation
      integer, intent(IN ) :: kday  ! current day    of simulation
      integer, intent(IN ) :: khou  ! current hour   of simulation
      integer, intent(IN ) :: kmin  ! current minute of simulation

      integer :: idall
      integer :: ilp
      integer :: iy400,id400
      integer :: iy100,id100
      integer :: iy004,id004
      integer :: jmon

      logical :: leap

!     m_days_*, NOT n_days_*: step2cal30 decodes with m_days_per_year and
!     m_days_per_month, and this is its inverse. n_days_per_year counts SIDEREAL
!     days and n_days_per_month is a twelfth of the 24-HOUR-day year, so the
!     pair does not even agree with itself -- twelve n_days_per_month is not
!     n_days_per_year. Callers pass mtspd as ktspd, which is now calmod's own.
!     world-x1k.
!     kyea-1, LIKE kmon AND kday. step2cal30 labels the first year 1, so year 1
!     month 1 day 1 is the epoch and encodes to step 0; with kyea taken as a
!     count of ELAPSED years it encoded to a whole year of steps and the two
!     were not inverses. That was not only a labelling defect: plasim's cold
!     start sets nstep = n_start_step from this call, and radmod takes the
!     orbital phase from mod(nstep,n_steps_per_year), so every cold start began
!     four fifths of an orbit past the meananomaly0 its own config declares.
      if (n_days_per_year /= 365) then ! simplified calendar
         kstep = ktspd * ((kyea-1) * m_days_per_year  &
                       +  (kmon-1) * m_days_per_month &
                       +   kday-1)
      return
      endif
      iy400 = kyea   /  400    ! segment [400]
      id400 = mod(kyea ,400)   ! year in segment [0..399]
      iy100 = id400   / 100    ! century [0,1,2,3]
      id100 = mod(id400,100)   ! year in century [0..99]
      iy004 = id100   /   4    ! tetrade [0..24]
      id004 = mod(id100,  4)   ! year in tetrade [0,1,2,3]

      leap  = (id004 == 0 .and. (id100 /= 0 .or. id400 == 0))

      ilp = -1
      if (id004 > 0) ilp = ilp + 1
      if (iy100 > 0 .and. id100 == 0 ) ilp = ilp + 1

      monaccu(0) = mondays(0)
      monaccu(1) = mondays(1)
      monaccu(2) = mondays(1) + mondays(2)
      if (leap) monaccu(2) = monaccu(2) + 1
      do jmon = 3 , 12
         monaccu(jmon) = monaccu(jmon-1) + mondays(jmon)
      enddo

      idall = iy400 * ny400d + iy100 * ny100d + iy004 * ny004d &
            + id004 * ny001d + monaccu(kmon-1)+ kday + ilp
      kstep = ktspd * idall + (ktspd * (khou * 60 + kmin)) / 1440

      return
      end subroutine cal2step

!     =====================
!     SUBROUTINE STEP2CAL30
!     =====================

      subroutine step2cal30(kstep,kdatim)
      use calmod
      implicit none
      integer, intent(IN ) :: kstep     ! time step since simulation start
      integer, intent(OUT) :: kdatim(7) ! year,month,day,hour,min,weekday,leapyear

      integer :: iyea  ! current year   of simulation
      integer :: imon  ! current month  of simulation
      integer :: iday  ! current day    of simulation
      integer :: ihou  ! current hour   of simulation
      integer :: imin  ! current minute of simulation
      integer :: idall
      integer :: istp

!       
!       idall = kstep / ktspd
!       iyea  = idall / n_days_per_year
!       idall = mod(idall,n_days_per_year)
!       imon  = idall / n_days_per_month + 1
! !       imon = mod(kstep/(n_days_per_year*ntspd / 12) + 1,12)
!       iday  = mod(idall,n_days_per_month) + 1
!       istp  = mod(kstep,ktspd)
!       imin  = (istp * solar_day) / (ktspd * 60)
!       ihou  = imin / 60
!       imin  = mod(imin,60)
      
      
      idall = kstep / mtspd                      !Which day is it (number)--steps/steps-per-1/360th
      iyea  = idall / m_days_per_year + 1        !Which year (number)--starting with 1: day/360+1
      idall = mod(idall,m_days_per_year)         !Which day of the year
      imon  = idall / m_days_per_month + 1       !Which month is it (day/30 + 1)
!       imon = mod(kstep/(n_days_per_year*ntspd / 12) + 1,12)
      iday  = mod(idall,m_days_per_month) + 1    !Which day of the month
      istp  = mod(kstep,mtspd)                   !Which step of the 1/360th day
      imin  = (istp * tcalday) / (mtspd * 60)   !Which minute of the day -- this might be wrong
              ! mtspd * 60 = seconds per 1/360th day
              ! istp * day_24hr = N steps * seconds per 24-hour day
              
              !  N dt   86400s    1/360th day   1 minute                  1/360th day
              ! ----- * ------ * ------------ * --------  =  M minutes * ------------
              !   1      24hr      mtspd dt       60 s                    24 hr day
              
              ! N dt   1/360th day      X s        1 minute
              ! ---- * ----------- * ----------- * -------- = M minutes  
              !  1      mtspd dt     1/360th day     60 s
              
      ihou  = imin / 60                          !Which hour of the day
      imin  = mod(imin,60)                       !Which minute of the hour

      kdatim(1) = iyea
      kdatim(2) = imon
      kdatim(3) = iday
      kdatim(4) = ihou
      kdatim(5) = imin
      kdatim(6) = 0 ! day of week
      kdatim(7) = 0 ! leap year

      return
      end subroutine step2cal30

!     =================
!     SUBROUTINE NTOMIN
!     =================

       subroutine ntomin(kstep,kmin,khou,kday,kmon,kyea)
       use calmod
       implicit none

       integer, intent(in ) :: kstep
       integer, intent(out) :: kmin,khou,kday,kmon,kyea
       integer :: idatim(7)

       if (n_days_per_year == 365) then
          call step2cal(kstep,ntspd,idatim)
       else
          call step2cal30(kstep,idatim)
       endif
       kyea = idatim(1)
       kmon = idatim(2)
       kday = idatim(3)
       khou = idatim(4)
       kmin = idatim(5)
       return
       end

!     =================
!     SUBROUTINE NTODAT
!     =================

      subroutine ntodat(istep,datch)
      character(len=18) datch
      character(len=3) mona(12)
      data mona /'Jan','Feb','Mar','Apr','May','Jun',                   &
     &           'Jul','Aug','Sep','Oct','Nov','Dec'/
      call ntomin(istep,imin,ihou,iday,imon,iyea)
      write (datch,20030) iday,mona(imon),iyea,ihou,imin
20030 format(i2,'-',a3,'-',i4.4,2x,i2,':',i2.2)
      end


!     =================
!     SUBROUTINE MOMINT
!     =================

!     Compute month indices and weights for time interpolation from
!     monthly mean data to current timestep

      subroutine momint(kperp,kstep,kmona,kmonb,pweight)
      use calmod
      implicit none
      integer, intent(in ) :: kperp   ! perpetual mode ?
      integer, intent(in ) :: kstep   ! target step
      integer, intent(out) :: kmona   ! current month (1-12)
      integer, intent(out) :: kmonb   ! next or previous month (0-13)
      real   , intent(out) :: pweight ! interpolation weight

      integer :: idatim(7) ! date time array
      integer :: iday
      integer :: ihour
      integer :: imin
      integer :: jmonb  ! next or previous month (1-12)

      real    :: zday   ! fractional day (including hour & minute)
      real    :: zdpma  ! days per month a
      real    :: zdpmb  ! days per month b
      real    :: zmeda  ! median day of the month a
      real    :: zmedb  ! median day of the month b

!     convert time step to date / time

      idatim(:) = 0
      if (kperp > 0) then                     ! perpetual date
         call yday2mmdd(kperp,idatim(2),idatim(3))
      else if (n_days_per_year == 365) then   ! real calendar
         call step2cal(kstep,ntspd,idatim)
      else                                    ! simple calendar
         call step2cal30(kstep,idatim)
      endif

      kmona = idatim(2)
      iday  = idatim(3)
      ihour = idatim(4)
      imin  = idatim(5)

!     set fractional day

      zday = iday + ((ihour * 60.0 + imin) * 60.0) / solar_day

!     compute median of month a

      zdpma = n_days_per_month
      if (n_days_per_year == 365) then
         zdpma = mondays(kmona)
         if (kmona == 2) zdpma = zdpma + idatim(7) ! leap year
      endif
      zmeda = 0.5 * (zdpma + 1.0) ! median day a

!     define neighbour month

      if (zday > zmeda) then
         kmonb = kmona + 1 !     next month (maybe 13)
      else
         kmonb = kmona - 1 ! previous month (maybe  0)
      endif

!     compute median of month b

      zdpmb = n_days_per_month
      if (n_days_per_year == 365) then
         jmonb = mod(kmonb+11,12) + 1 ! convert month (0-13) -> (1-12)
         zdpmb = mondays(jmonb)
         if (jmonb == 2) zdpmb = zdpmb + idatim(7) ! leap year
      endif
      zmedb = 0.5 * (zdpmb + 1.0) ! median day b

!     compute weight

      pweight = abs(zday - zmeda) / (zmeda + zmedb - 1.0)

      return
      end subroutine momint

