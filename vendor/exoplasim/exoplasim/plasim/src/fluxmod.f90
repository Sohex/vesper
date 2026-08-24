      module fluxmod
      use pumamod
!
!     version identifier
!

      character(len=80) :: version = '22.11.2002 by Larry'

!
!     global parameter
!

      parameter(vonkarman=0.4)    ! von karman const.

!
!     namelist parameter:
!

      integer :: nvdiff  = 1      ! switch for vertical diffusion
      integer :: nshfl   = 1      ! switch for sensible heat flux
      integer :: nevap   = 1      ! switch for evaporation/latend heat flux
      integer :: nstress = 1      ! switch for wind stress
      integer :: ntsa    = 2      ! flag for near surface temp calc.

!     THE BOUNDARY LAYER'S FIVE CONSTANTS, and every one of them is in
!     fluxmod_nl and has never been written by any run of this project: the
!     staged fluxmod_namelist is an empty group. They are named here rather than
!     left as bare literals so that a reader knows what a run would be changing.
!     world-e2k.
!
!     zumin is a gustiness floor standing in for the unresolved wind variance
!     that keeps the surface exchange from collapsing in a calm cell. It is in
!     m/s and it is used at two places against two DIFFERENT quantities: a
!     squared speed in mktcoe, where it is squared, and a wind DIFFERENCE
!     between adjacent levels in vdiff, where it is not. The comment used to say
!     m/s while the first site compared it against m2/s2; at the declared 1 m/s
!     the two forms coincide, which is why nothing noticed.
      real :: zumin      = 1.    ! gustiness floor on the exchange wind (m/s)
!
!     vdiff_lamm is the asymptotic mixing length, in metres: the single number
!     that sets free-tropospheric vertical diffusivity, through
!     zmixm = lamm*k*z/(lamm + k*z), and through zlamh below it the heat
!     version. 160 m is ECHAM's and is a fit to Earth's free troposphere. It is
!     not obviously transferable, and it is left at ECHAM's value because
!     nothing in this project has measured a replacement -- not because 160 m
!     has been shown to apply here.
      real :: vdiff_lamm = 160.  ! asymptotic mixing length (m)
!
!     The three Louis stability-function coefficients, cited in mktcoe to ECHAM
!     REPORT 218, which is a report and not a derivation. None of the three is
!     separately justified there or here; they are one fitted set and are used
!     as one.
      real :: vdiff_b    = 5.    ! Louis stability function, unstable slope
      real :: vdiff_c    = 5.    ! Louis stability function, unstable denominator
      real :: vdiff_d    = 5.    ! Louis stability function, stable branch
!
!     The reference temperature that converts a sigma level to a height, through
!     the hypsometric zzlev = -gascon*ztscal*ln(sigmah)/ga. gascon and ga are
!     this world's; 250 K is a representative mid-tropospheric temperature and
!     was a bare `parameter` inside vdiff with no comment and no namelist route.
      real :: ztscal     = 250.  ! reference temperature for sigma-to-height (K)
!
!     THE FREE-CONVECTION COEFFICIENT, derived in fluxini rather than declared.
!     Miller et al. (1992)'s unstable-ocean transfer collapses in the convective
!     limit to 0.0016*dth^(1/3)/(|U|*C_N), which is the convective velocity
!     scale w* = (g H z_i / T)^(1/3) with the gravity and the inversion height
!     folded into the coefficient. There is no ga anywhere in the expression, so
!     on a 12.81 m/s2 world the coefficient was low by (12.81/9.80665)^(1/3) =
!     1.093 and free-convection latent and sensible exchange over calm unstable
!     ocean was about 9 per cent weak -- over the warm ocean, where the
!     evaporation is. Scaled here, once, rather than at every gridpoint.
      real :: freeconv   = 0.0   ! set by fluxini from ga

!
!     arrays
!

      real :: dtransh(NHOR)       ! transfere coeficient for heat
      real :: dtransm(NHOR)       ! transfere coeficient for momentum

!
!     real scalars
!

      real :: time4fl = 0.        ! CPU time needed for fluxmod
      real :: time4sf = 0.        ! CPU time needed for surface fluxes
      real :: time4tr = 0.        ! CPU tiem needed for transfere coeff.
      real :: time4sh = 0.        ! CPU time needed for sensible heat flux
      real :: time4ev = 0.        ! CPU time needed for evaporation
      real :: time4st = 0.        ! CPU time needed for wind stress
      real :: time4vd = 0.        ! CPU time needed for vertical diffusion


!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!$omp threadprivate(dtransh,dtransm,nevap,nshfl,nstress,ntsa,nvdiff,time4ev,time4fl,time4sf,time4sh,&
!$omp&  time4st,time4tr,time4vd,vdiff_b,vdiff_c,vdiff_d,vdiff_lamm,version,zumin,ztscal,freeconv)

      end module fluxmod

!     ==================
!     SUBROUTINE FLUXINI
!     ==================

      subroutine fluxini
      use fluxmod
!
      namelist/fluxmod_nl/nvdiff,nshfl,nevap,nstress,ntsa                  &
     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d,ztscal
!
      if(mypid==NROOT) then
         open(11,file=fluxmod_namelist)
         read(11,fluxmod_nl)
         close(11)
         write(nud,'(/," ***********************************************")')
         write(nud,'(" * FLUXMOD ",a35," *")') trim(version)
         write(nud,'(" ***********************************************")')
         write(nud,'(" * Namelist FLUXMOD_NL from <fluxmod_namelist> *")')
         write(nud,'(" ***********************************************")')
         write(nud,fluxmod_nl)
      endif

      call mpbci(nvdiff)
      call mpbci(nshfl)
      call mpbci(nevap)
      call mpbci(nstress)
      call mpbci(ntsa)
      call mpbcr(zumin)
      call mpbcr(vdiff_lamm)
      call mpbcr(vdiff_b)
      call mpbcr(vdiff_c)
      call mpbcr(vdiff_d)
      call mpbcr(ztscal)

!     ONCE, FROM THIS WORLD'S GRAVITY. See the declaration: 0.0016 is Miller's
!     coefficient with Earth's g folded in, and w* goes as g^(1/3). world-e2k.
      freeconv = 0.0016 * (ga / 9.80665)**(1./3.)
      if (mypid == NROOT) then
         write(nud,'(" * free-convection coefficient ",f12.8," *")') freeconv
      endif

      return
      end subroutine fluxini

!     ===================
!     SUBROUTINE FLUXSTEP
!     ===================

      subroutine fluxstep
      use fluxmod

      if(ntime == 1) call mksecond(zsec,0.)
      call surflx
      if(ntime == 1) then
       call mksecond(zsec1,zsec)
       time4sf=time4sf+zsec1
       call mksecond(zsec1,0.)
      endif
      if(nvdiff > 0) call vdiff
      if(ntime == 1) then
       call mksecond(zsec1,zsec1)
       time4vd=time4vd+zsec1
       call mksecond(zsec,zsec)
       time4fl=time4fl+zsec
      endif

      return
      end subroutine fluxstep

!     ===================
!     SUBROUTINE FLUXSTOP
!     ===================

      subroutine fluxstop
      use fluxmod

      if(mypid == NROOT .and. ntime == 1) then
       write(nud,*)'*********************************************'
       write(nud,*)' CPU usage in FLUXSTEP (ROOT process only):  '
       write(nud,*)'    All routines      : ',time4fl,' s'
       write(nud,*)'    All surface fluxes: ',time4sf,' s'
       write(nud,*)'    Vertical diffusion: ',time4vd,' s'
       write(nud,*)'    Transfere coeff   : ',time4tr,' s'
       write(nud,*)'    Sensible heat flux: ',time4sh,' s'
       write(nud,*)'    Evaporation       : ',time4ev,' s'
       write(nud,*)'    Wind stress       : ',time4st,' s'
       write(nud,*)'*********************************************'
      endif

      return
      end subroutine fluxstop

!
!     =================
!     SUBROUTINE SURFLX
!     =================
!
      subroutine surflx
      use fluxmod
!
!     local parameters
!

!
!     local arrays
!

      real zabsu2(NHOR)    ! squared wind speed (m2/s2)
      real znl(NHOR)       ! z of lowermost level (m)
      real zbz0(NHOR)      ! z/z0
      real zri(NHOR)       ! bulk richardson number
      real zrifh(NHOR)     ! factor for heat flux transfere coef.
      real zrifm(NHOR)     ! factor for momentum flux transfere coef.
!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)

      if(ntime == 1) call mksecond(zsec,0.)

!
!*    set some const.:
!

      zexp=-akap
      zlnsig=ALOG(sigma(NLEV))

!
!*    surface air temperature
!

      if(ntsa==1) then
       dtsa(:)=dt(:,NLEV)*sigma(NLEV)**zexp
      elseif(ntsa==2) then
       dtsa(:)=dt(:,NLEV)*sigma(NLEV)**zexp                             &
     &        *(1.+(1./rdbrv-1.)*dq(:,NLEV))
      else
       write(nud,*) '!ERROR: wrong ntsa in fluxmod (only 1 or 2 valied)!'
       stop
      end if

!
!*    calculate transfere coef:
!
!     windspeed (squared):
!

!     SQUARED, because zabsu2 is a squared speed and zumin is a speed. At the
!     declared 1 m/s this is the same number; at any other value it was the
!     wrong one. world-e2k.
      zabsu2(:)=AMAX1(zumin*zumin,du(:,NLEV)*du(:,NLEV)+dv(:,NLEV)*dv(:,NLEV))

!
!     z off lowermost layer
!

      znl(:)=-gascon*0.5*(dt(:,NLEV)+dtsa(:))*zlnsig/ga
      
      
!
!     z/z0
!

      zbz0(:)=znl(:)/dz0(:)
      
      do jhor=1,NHOR
         if (zbz0(jhor) .le. -1.0) write(nud,*) dt(jhor,NLEV),dtsa(jhor),zlnsig,dz0(jhor)
      enddo
!
!     bulk richardson number
!

      if(ntsa==1) then
       zri(:)=ga*znl(:)*(dtsa(:)-dt(:,NLEP))/(zabsu2(:)*dtsa(:))
      elseif(ntsa==2) then
       zri(:)=ga*znl(:)/(zabsu2(:)*dtsa(:))                             &
     &       *(dtsa(:)-dt(:,NLEP)*(1.+(1./rdbrv-1.)*dq(:,NLEP)))
      endif

      do jhor=1,NHOR
       if (zbz0(jhor)+1. .le. -1.0) then
          write(nud,*) "negative z/z0"
          write(nud,*) zbz0(jhor)
       endif
       if (ALOG(zbz0(jhor)+1.) == 0.0) write(nud,*) zbz0(jhor)
       zkblnz2=(vonkarman/ALOG(zbz0(jhor)+1.))**2

!
!     richardson number dependent factor
!     see ECHAM REPORT 218
!     (unstable zrifh over ocean according to Miller et al. 1992)
!

       if(zri(jhor) <= 0.) then
        zdenom=1.+3.*vdiff_c*vdiff_b                                    &
     &        *sqrt(-1.*zri(jhor)*(zbz0(jhor)+1.))*zkblnz2
        zrifm(jhor)=1.-2.*vdiff_b*zri(jhor)/zdenom
        if(dls(jhor) < 1.) then
         zdth=-(dtsa(jhor)                                              &
     &         -dt(jhor,NLEP)*(1.+(1./rdbrv-1.)*dq(jhor,NLEP)))
         zrifh(jhor)=(1.+(freeconv*zdth**(1./3.)/SQRT(zabsu2(jhor))     &
     &                  /zkblnz2)**1.25)**0.8
        else
         zrifh(jhor)=1.-3.*vdiff_b*zri(jhor)/zdenom
        endif
       else
        zdenom=SQRT(1+vdiff_d*zri(jhor))
        zrifh(jhor)=1./(1.+3.*vdiff_b*zri(jhor)*zdenom)
        zrifm(jhor)=1./(1.+2.*vdiff_b*zri(jhor)/zdenom)
       endif

!
!     transfere coeff.
!

       ztrans=SQRT(zabsu2(jhor))*zkblnz2
       dtransh(jhor)=ztrans*zrifh(jhor)
       dtransm(jhor)=ztrans*zrifm(jhor)
      enddo

!
!*    dbug diagnostic output if nprint=2 (pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       call mpgagp(zprf1,zabsu2,1)
       call mpgagp(zprf2,znl,1)
       call mpgagp(zprf3,zri,1)
       call mpgagp(zprf4,zrifh,1)
       call mpgagp(zprf5,zrifm,1)
       call mpgagp(zprf6,dtransh,1)
       call mpgagp(zprf7,dtransm,1)
       if(mypid==NROOT) then
        write(nud,*)'In Surflx:'
        write(nud,*)'abs(u)= ',SQRT(zprf1(nprhor)),' znl= ', zprf2(nprhor)
        write(nud,*)'ri-number= ',zprf3(nprhor)
        write(nud,*)'rifh= ',zprf4(nprhor),' rifm= ', zprf5(nprhor)
        write(nud,*)'transh= ',zprf6(nprhor),' transm= ', zprf7(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
       deallocate(zprf6)
       deallocate(zprf7)
      endif

      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4tr=time4tr+zsec
       call mksecond(zsec,0.)
      endif

!
!*    calc. fluxes and tendencies
!

      if (nstress > 0 ) call mkstress   ! wind stress
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4st=time4st+zsec
       call mksecond(zsec,0.)
      endif
      if (nshfl   > 0 ) call mkshfl     ! sensible heat flux
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4sh=time4sh+zsec
       call mksecond(zsec,0.)
      endif
      if (nevap   > 0 ) call mkevap     ! evaporation/latent heat flux
      if(ntime == 1) then
       call mksecond(zsec,zsec)
       time4ev=time4ev+zsec
      endif

      return
      end subroutine surflx

!     ===================
!     SUBROUTINE MKSTRESS
!     ===================

      subroutine mkstress
      use fluxmod
!
!*    calculate u,v tendencies due to surface wind stress
!
!     implicit formulation is used:
!
!     (u(t+dt)-u(t))/dt=g/ps/dsigma* (cd*abs(u)*rho*u(t+dt))
!     => u(t+dt) => dudt => stress=cd*abs(u)*rho*u(t+dt)
!
!
!     local arrays
!

      real zun(NHOR),zvn(NHOR)         ! u,v at new time level
      real zkdiff(NHOR)                ! transfere coefficients
      real zdtdt(NHOR)                 ! heating due to dissipation

!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)
      real, allocatable :: zprf6(:)
      real, allocatable :: zprf7(:)

!
!*    set some const.
!

      zkonst1=ga*deltsec2/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/deltsec2/ga

!
!*    modified transfere coeficient (incl. some const.)
!

      zkdiff(:)=zkonst1*dtransm(:)/dt(:,NLEP)

!
!*    wind at t+dt
!

      zun(:)=du(:,NLEV)/(1.+zkdiff(:))
      zvn(:)=dv(:,NLEV)/(1.+zkdiff(:))

!
!*    add tendencies
!

      dudt(:,NLEV)=dudt(:,NLEV)+(zun(:)-du(:,NLEV))/deltsec2
      dvdt(:,NLEV)=dvdt(:,NLEV)+(zvn(:)-dv(:,NLEV))/deltsec2

!
!*    use loss of kinetic energy to warm the lower layer 
!     (for energy conservation)
!

      zdtdt(:)=-(zun(:)*zun(:)-du(:,NLEV)*du(:,NLEV)                    &
                +zvn(:)*zvn(:)-dv(:,NLEV)*dv(:,NLEV))/deltsec2          &
     &          *0.5/acpd/(1.+adv*dq(:,NLEV))
      if(ndheat > 0) then
       dtdt(:,NLEV)=dtdt(:,NLEV)+zdtdt(:)
      endif

!
!*    wind stress
!

      dtaux(:)=dp(:)*zkonst2*zkdiff(:)*zun(:)
      dtauy(:)=dp(:)*zkonst2*zkdiff(:)*zvn(:)

!
!*    ustar**3
!

      dust3(:)=SQRT(dtaux(:)**2+dtauy(:)**2)*gascon*dt(:,NLEP)/dp(:)
      dust3(:)=SQRT(dust3(:))
      dust3(:)=dust3(:)**3

!
!*    dbug diagnostic output if nprint=2 (pumamod)
!

      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       allocate(zprf6(NLON*NLAT))
       allocate(zprf7(NLON*NLAT))
       call mpgagp(zprf1,zun,1)
       call mpgagp(zprf2,zvn,1)
       call mpgagp(zprf3,du(1,NLEV),1)
       call mpgagp(zprf4,dv(1,NLEV),1)
       call mpgagp(zprf5,dtaux,1)
       call mpgagp(zprf6,dtauy,1)
       call mpgagp(zprf7,dust3,1)
       if(mypid==NROOT) then
        write(nud,*)'In mkstress:'
        write(nud,*)'u-new= ',zprf1(nprhor),' v-new= ', zprf2(nprhor)
        write(nud,*)'du/dt (m/s2)= ',(zprf1(nprhor)-zprf3(nprhor))/deltsec2  &
     &        ,' dv/dt (m/s2)= ',(zprf2(nprhor)-zprf4(nprhor))/deltsec2
        write(nud,*)'taux= ',zprf5(nprhor),' tauy= ', zprf6(nprhor)
        write(nud,*)'ust3= ',zprf7(nprhor)
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
!     energy diagnostics
!
      if(nenergy > 0) then
       denergy(:,21)=zdtdt(:)*acpd*(1.+adv*dq(:,NLEV))*dp(:)            &
     &              /ga*dsigma(NLEV) 
       if(nener3d > 0) then
        dener3d(:,1:NLEM,21)=0.
        dener3d(:,NLEV,21)=denergy(:,21)
       endif
      endif
!
      return
      end subroutine mkstress

!     =================
!     SUBROUTINE MKSHFL
!     =================

      subroutine mkshfl
      use fluxmod
!
!
!*    calculate t tendency due to surface sensible heat flux
!
!     implicit formulation is used:
!
!     (t(t+dt)-t(t))/dt=-g/cp/ps/dsigma*(cp*cd*abs(u)*rho*(t'(t+dt)-ts)
!     => t(t+dt) => dtdt => sensible heat flux=cp*cd*abs(u)*rho*(t'(t+dt)-ts)
!     =>d(sensible heat flux)/dts (needed for ts calculation)
!     t'= t*(1+(Rv/rd-1)*q)/sigma(NLEV)**akap
!
!
!     local arrays
!

      real ztn(NHOR)      ! new temperature
      real zkdiff(NHOR)   ! modified transfere coef.
      real zfac(NHOR)     ! factor to extrapolate temp. to surface
                          ! (needs to be independent of temp)
!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)

!
!*    set some const.
!

      zexp=-akap
      zkonst1=ga*deltsec2/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/deltsec2/ga
      zfac(:)=sigma(NLEV)**(zexp)

!
!*    modified transfere coeficien (incl. some const.)
!

      zkdiff(:)=zkonst1*dtransh(:)/dt(:,NLEP)

!
!*    temperatue at t+dt
!

      ztn(:)=(dt(:,NLEV)+zkdiff(:)*dt(:,NLEP))/(1.+zkdiff(:)*zfac(:))

!
!*    add tendency
!

      dtdt(:,NLEV)=dtdt(:,NLEV)+(ztn(:)-dt(:,NLEV))/deltsec2

!
!*    sensible heat flux and dshfl/dsurface_temperature
!

      dshfl(:)=zkonst2*zkdiff(:)*dp(:)*acpd*AMAX1(1.,1.+ADV*dq(:,NLEV)) &
     &        *(ztn(:)*zfac(:)-dt(:,NLEP))
      dshdt(:)=-1.*zkonst2*zkdiff(:)*dp(:)                              &
     &        *acpd*AMAX1(1.,1.+ADV*dq(:,NLEV))
!
      dtsa(:)=ztn(:)*zfac(:)
!
!
!*    dbug diagnostic output if nprint=2 (pumamod)
!
      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       call mpgagp(zprf1,ztn,1)
       call mpgagp(zprf2,dt(1,NLEV),1)
       call mpgagp(zprf3,zfac,1)
       call mpgagp(zprf4,dshfl,1)
       call mpgagp(zprf5,dshdt,1)
       if(mypid==NROOT) then
        zzts=zprf2(nprhor)*zprf3(nprhor)
        zztsn=zprf1(nprhor)*zprf3(nprhor)
        write(nud,*)'In mkshfl:'
        write(nud,*)'t-new= ',zprf1(nprhor),' tsa(old)= ',zzts,' tsa(new)= ' &
     &        ,zztsn
        write(nud,*)'dt/dt (K/s)= ',(zprf1(nprhor)-zprf2(nprhor))/deltsec2
        write(nud,*)'shfl= ',zprf4(nprhor),' dshdt= ', zprf5(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEM,3)=0.
       dgp3d(:,NLEV,3)=(ztn(:)-dt(:,NLEV))/deltsec2
      end if
!
!     energy diagnostics
!
      if(nenergy > 0) then
       denergy(:,7)=(ztn(:)-dt(:,NLEV))/deltsec2                        &
     &             *acpd*(1.+adv*dq(:,NLEV))*dp(:)/ga*dsigma(NLEV)
       if(nener3d > 0) then
        dener3d(:,1:NLEM,7)=0.
        dener3d(:,NLEV,7)=denergy(:,7)
       endif  
      endif
!
      return
      end subroutine mkshfl

!     =================
!     SUBROUTINE MKEVAP
!     =================

      subroutine mkevap
      use fluxmod
!
!*    calculate q tendency due to surface evaporation
!
!     implicit formulation is used:
!
!     (q(t+dt)-q(t))/dt=-g/ps/dsigma*(rhs*cd*abs(u)*rho*(q(t+dt)-qs)
!     => q(t+dt) => dqdt => evaporation=rhs*cd*abs(u)*rho*(q(t+dt)-qs)
!     =>latent heat flux=>d(latent heat flux)/dts (needed for ts calculation)
!     rhs= soil wetness factor
!
!     local arrays
!
!
      real zqn(NHOR)      ! new q
      real zkdiff(NHOR)   ! modified transfere coef.

!
!     allocatable arrays for dbug diagnostics
!

      real, allocatable :: zprf1(:)
      real, allocatable :: zprf2(:)
      real, allocatable :: zprf3(:)
      real, allocatable :: zprf4(:)
      real, allocatable :: zprf5(:)

!
!*    set some const.
!

      zkonst1=ga*deltsec2/(gascon*dsigma(NLEV))
      zkonst2=dsigma(NLEV)/deltsec2/ga

!
!*    modified transfere coeficient (incl. some const.)
!

      zkdiff(:)=drhs(:)*zkonst1*dtransh(:)/dt(:,NLEP)

!
!*    q at t+dt
!

      zqn(:)=AMAX1(dq(:,NLEV)                                           &
     &            ,(dq(:,NLEV)+zkdiff(:)*dq(:,NLEP))/(1.+zkdiff(:)))

!
!     limit evporation to available water
!
      where(dls(:) > 0.)
       zqn(:)=AMIN1(zqn(:)                                              &
     &             ,dwatc(:)/deltsec*1000./dp(:)/zkonst2+dq(:,NLEV))
      endwhere

!
!*    add tendency
!

      dqdt(:,NLEV)=dqdt(:,NLEV)+(zqn(:)-dq(:,NLEV))/deltsec2

!
!*    evaporation, latent heat flux and dlhfl/dtemperature
!

      devap(:)=-dp(:)*zkonst2/1000.*(zqn(:)-dq(:,NLEV))

!     THE SATURATION PHASE FOLLOWS THE ARM, not the temperature. The mask below
!     is the model's own statement of which phase is evaporating, and it is not
!     the temperature test alone: a cell with dls < 0.5 takes the liquid arm
!     however cold it is, because the water surface under a partial ice cover is
!     still water and icemod owns the ice. Selecting the coefficients on
!     dt < TMELT here would put an ice saturation under a liquid latent heat on
!     exactly those cells. world-ako.
      where(dt(:,NLEP) > TMELT .or. dls(:) < 0.5)
       dlhfl(:)=devap(:)*ALV*1000.
       dlhdt(:)=-1.*ALV*zkdiff(:)*zkonst2*dp(:)                         &
     &         *ra2*(TMELT-ra4)*dq(:,NLEP)/(dt(:,NLEP)-ra4)**2
      elsewhere
       dlhfl(:)=devap(:)*ALS*1000.
       dlhdt(:)=-1.*ALS*zkdiff(:)*zkonst2*dp(:)                         &
     &         *ra2i*(TMELT-ra4i)*dq(:,NLEP)/(dt(:,NLEP)-ra4i)**2
      endwhere
      where(dlhfl(:) == 0.) dlhdt(:)=0.
!
!*    dbug diagnostic output if nprint=2 (pumamod)
!
      if(nprint==2) then
       allocate(zprf1(NLON*NLAT))
       allocate(zprf2(NLON*NLAT))
       allocate(zprf3(NLON*NLAT))
       allocate(zprf4(NLON*NLAT))
       allocate(zprf5(NLON*NLAT))
       call mpgagp(zprf1,zqn,1)
       call mpgagp(zprf2,dq(1,NLEV),1)
       call mpgagp(zprf3,devap,1)
       call mpgagp(zprf4,dlhfl,1)
       call mpgagp(zprf5,dlhdt,1)
       if(mypid==NROOT) then
        write(nud,*)'In mkevap:'
        write(nud,*)'q-new= ',zprf1(nprhor)                                  &
     &        ,' dq/dt (1/s)= ',(zprf1(nprhor)-zprf2(nprhor))/deltsec2
        write(nud,*)'lhfl= ',zprf4(nprhor),' dlhdt= ', zprf5(nprhor)
        write(nud,*)'evap (mm/day)= ',1000.*ntspd*0.5*deltsec2*zprf3(nprhor)
       endif
       deallocate(zprf1)
       deallocate(zprf2)
       deallocate(zprf3)
       deallocate(zprf4)
       deallocate(zprf5)
      endif

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEM,13)=0.
       dgp3d(:,NLEV,13)=(zqn(:)-dq(:,NLEV))/deltsec2*dp(:)*dsigma(NLEV) &
     &                 /ga/1000.
      end if

!
!     entropy diagnostics
!

!
      return
      end subroutine mkevap

!     ================
!     SUBROUTINE VDIFF
!     ================

      subroutine vdiff
      use fluxmod
!
!     calculate t,q,u,v tendencies due to vertical diffusion
!     using the ECHAM semi-implicit scheme
!
!
      real zdtdt(NHOR,NLEV)
      real zdudt(NHOR,NLEV)
      real zdvdt(NHOR,NLEV)
      real zdqdt(NHOR,NLEV)
!
      real zkdiffm(NHOR,NLEV),zkdiffh(NHOR,NLEV)
      real zabsu(NHOR,NLEV)
      real ztn(NHOR,NLEV)
      real zun(NHOR,NLEV)
      real zvn(NHOR,NLEV)
      real zqn(NHOR,NLEV)
      real zt(NHOR,NLEV)
      real zu(NHOR,NLEV)
      real zv(NHOR,NLEV)
      real zq(NHOR,NLEV)
      real zebs(NHOR,NLEM)
      real zskap(NLEV),zskaph(NLEV)
      real zke(NHOR,NLEV),zken(NHOR,NLEV)
!
      zlamh=vdiff_lamm*SQRT(3.*vdiff_d*0.5)
      zkonst1=ga*deltsec2/gascon
      zkonst2=zkonst1*ga/gascon
      zkonst3=zkonst2*ga/gascon
!
      zskap(:)=sigma(:)**akap
      zskaph(:)=sigmah(:)**akap
!
      zt(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2
      zq(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
      zu(:,:)=du(:,1:NLEV)+dudt(:,1:NLEV)*deltsec2
      zv(:,:)=dv(:,1:NLEV)+dvdt(:,1:NLEV)*deltsec2

!
!     0) modified diffusion coefficents
!
      zabsu(:,:)=SQRT(zu(:,1:NLEV)**2+zv(:,1:NLEV)**2)
!
      do jlev=1,NLEM
       jlp=jlev+1
       zzlev=-gascon*ztscal*ALOG(sigmah(jlev))/ga
       zdz=gascon*ztscal*ALOG(sigma(jlp)/sigma(jlev))/ga
       zmixh=zlamh*vonkarman*zzlev/(zlamh+vonkarman*zzlev)
       zmixm=vdiff_lamm*vonkarman*zzlev/(vdiff_lamm+vonkarman*zzlev)
       zmixh2=zmixh*zmixh
       zmixm2=zmixm*zmixm
       zrdsig=1./(sigma(jlp)-sigma(jlev))
       zfac=zkonst3*sigmah(jlev)*sigmah(jlev)*sigmah(jlev)              &
     &     *zrdsig
       zzkfac=SQRT((((((zzlev+zdz)/zzlev)**(1./3.)-1.)/zdz)**3)/zzlev)
       do jhor=1,NHOR
        zrthl=(dsigma(jlp)+dsigma(jlev))                                &
     &       /(zt(jhor,jlev)*dsigma(jlp)+zt(jhor,jlp)*dsigma(jlev))
        zdvds=AMAX1(zumin,ABS(zabsu(jhor,jlev)-zabsu(jhor,jlp)))        &
     &       *zrdsig
        zqf1=1.+(1./rdbrv-1.)*zq(jhor,jlev)
        zqf2=1.+(1./rdbrv-1.)*zq(jhor,jlp)
        zdthds=(zt(jhor,jlev)*zqf1/zskap(jlev)                          &
     &         -zt(jhor,jlp)*zqf2/zskap(jlp))                           &
     &        *zrdsig
        zri=zskaph(jlev)*gascon*zdthds                                  &
     &     /(sigmah(jlev)*zdvds*zdvds)
        if(zri <= 0.) then
         zdenomh=1.+3.*vdiff_c*vdiff_b*sqrt(-1.*zri)*zzkfac*zmixh2
         zdenomm=1.+3.*vdiff_c*vdiff_b*sqrt(-1.*zri)*zzkfac*zmixm2
         zrifh=1.-3.*vdiff_b*zri/zdenomh
         zrifm=1.-2.*vdiff_b*zri/zdenomm
        else
         zdenom=SQRT(1+vdiff_d*zri)
         zrifh=1./(1.+3.*vdiff_b*zri*zdenom)
         zrifm=1./(1.+2.*vdiff_b*zri/zdenom)
        endif
        zkdiffh(jhor,jlev)=zfac*zmixh2*zrifh*zdvds*zrthl*zrthl*zrthl
        zkdiffm(jhor,jlev)=zfac*zmixm2*zrifm*zdvds*zrthl*zrthl*zrthl
       enddo
      enddo

!
!     2. semi implicit scheme
!
!     2a momentum
!
!     top layer elimination
!
      zebs(:,1)=zkdiffm(:,1)/(dsigma(1)+zkdiffm(:,1))
      zun(:,1)=dsigma(1)*zu(:,1)/(dsigma(1)+zkdiffm(:,1))
      zvn(:,1)=dsigma(1)*zv(:,1)/(dsigma(1)+zkdiffm(:,1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffm(:,jlev)                                     &
     &             /(dsigma(jlev)+zkdiffm(:,jlev)                       &
     &              +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
       zun(:,jlev)=(zu(:,jlev)*dsigma(jlev)+zkdiffm(:,jlem)*zun(:,jlem))&
     &            /(dsigma(jlev)+zkdiffm(:,jlev)                        &
     &             +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
       zvn(:,jlev)=(zv(:,jlev)*dsigma(jlev)+zkdiffm(:,jlem)*zvn(:,jlem))&
     &            /(dsigma(jlev)+zkdiffm(:,jlev)                        &
     &             +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
      enddo
!
!     bottom layer elimination
!
      zun(:,NLEV)=(zu(:,NLEV)*dsigma(NLEV)+zkdiffm(:,NLEM)*zun(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffm(:,NLEM)*(1.-zebs(:,NLEM)))
      zvn(:,NLEV)=(zv(:,NLEV)*dsigma(NLEV)+zkdiffm(:,NLEM)*zvn(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffm(:,NLEM)*(1.-zebs(:,NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zun(:,jlev)=zun(:,jlev)+zebs(:,jlev)*zun(:,jlep)
       zvn(:,jlev)=zvn(:,jlev)+zebs(:,jlev)*zvn(:,jlep)
      enddo
!
!     tendencies
!
      zdudt(:,1:NLEV)=(zun(:,1:NLEV)-zu(:,1:NLEV))/deltsec2
      zdvdt(:,1:NLEV)=(zvn(:,1:NLEV)-zv(:,1:NLEV))/deltsec2
      dudt(:,1:NLEV)=dudt(:,1:NLEV)+zdudt(:,1:NLEV)
      dvdt(:,1:NLEV)=dvdt(:,1:NLEV)+zdvdt(:,1:NLEV)
!
!     temp.-tendencies due to the loss of kin.energy (energy-conservation)
!
      if(ndheat > 0) then
!
!     compute diffusion of ekin (needs to be substracted)
!
       zke(:,:)=0.5*zabsu(:,:)*zabsu(:,:)
       zebs(:,1)=zkdiffm(:,1)/(dsigma(1)+zkdiffm(:,1))
       zken(:,1)=dsigma(1)*zke(:,1)/(dsigma(1)+zkdiffm(:,1))
       do jlev=2,NLEM
        jlem=jlev-1
        zebs(:,jlev)=zkdiffm(:,jlev)                                    &
     &              /(dsigma(jlev)+zkdiffm(:,jlev)                      &
     &               +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
        zken(:,jlev)=(zke(:,jlev)*dsigma(jlev)                          &
     &               +zkdiffm(:,jlem)*zken(:,jlem))                     &
     &              /(dsigma(jlev)+zkdiffm(:,jlev)                      &
     &               +zkdiffm(:,jlem)*(1.-zebs(:,jlem)))
       enddo
       zken(:,NLEV)=(zke(:,NLEV)*dsigma(NLEV)                           &
     &              +zkdiffm(:,NLEM)*zken(:,NLEM))                      &
     &             /(dsigma(NLEV)+zkdiffm(:,NLEM)*(1.-zebs(:,NLEM)))
       do jlev=NLEM,1,-1
        jlep=jlev+1
        zken(:,jlev)=zken(:,jlev)+zebs(:,jlev)*zken(:,jlep)
       enddo
!
!     temperature tendencies du to ekin dissipation
!
       zdtdt(:,1:NLEV)=-((zun(:,1:NLEV)*zun(:,1:NLEV)                   &
     &                   -zu(:,1:NLEV)*zu(:,1:NLEV)                     &
     &                   +zvn(:,1:NLEV)*zvn(:,1:NLEV)                   &
     &                   -zv(:,1:NLEV)*zv(:,1:NLEV))/deltsec2           &
     &                   -(zken(:,1:NLEV)-zke(:,1:NLEV))/deltsec2)      &
     &                 *0.5/acpd/(1.+adv*dq(:,1:NLEV))
       dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdt(:,1:NLEV)
       zdtdt(:,1:NLEV)=0.
      endif
!
!     2b moisture
!
!     top layer elimination
!
      zebs(:,1)=zkdiffh(:,1)/(dsigma(1)+zkdiffh(:,1))
      zqn(:,1)=dsigma(1)*zq(:,1)/(dsigma(1)+zkdiffh(:,1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffh(:,jlev)                                     &
     &             /(dsigma(jlev)+zkdiffh(:,jlev)                       &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)))
       zqn(:,jlev)=(zq(:,jlev)*dsigma(jlev)+zkdiffh(:,jlem)*zqn(:,jlem))&
     &            /(dsigma(jlev)+zkdiffh(:,jlev)                        &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)))
      enddo
!
!     bottom layer elimination
!
      zqn(:,NLEV)=(zq(:,NLEV)*dsigma(NLEV)+zkdiffh(:,NLEM)*zqn(:,NLEM)) &
     &           /(dsigma(NLEV)+zkdiffh(:,NLEM)*(1.-zebs(:,NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       zqn(:,jlev)=zqn(:,jlev)+zebs(:,jlev)*zqn(:,jlep)
      enddo
!
!     tendencies
!
      zdqdt(:,1:NLEV)=(zqn(:,1:NLEV)-zq(:,1:NLEV))/deltsec2
      do jlev=1,NLEV
       dqdt(:,jlev)=dqdt(:,jlev)+zdqdt(:,jlev)
      enddo
!
!     2c potential temperature
!
      do jlev=1,NLEM
       zkdiffh(:,jlev)=zkdiffh(:,jlev)*zskaph(jlev)
      enddo
!
!     semi implicit scheme
!
!     top layer elimination
!
      zebs(:,1)=zkdiffh(:,1)/(dsigma(1)+zkdiffh(:,1)/zskap(1))
      ztn(:,1)=dsigma(1)*zt(:,1)/(dsigma(1)+zkdiffh(:,1)/zskap(1))
!
!     middle layer elimination
!
      do jlev=2,NLEM
       jlem=jlev-1
       zebs(:,jlev)=zkdiffh(:,jlev)                                     &
     &             /(dsigma(jlev)+(zkdiffh(:,jlev)                      &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))     &
     &              /zskap(jlev))
       ztn(:,jlev)=(zt(:,jlev)*dsigma(jlev)                             &
     &             +zkdiffh(:,jlem)/zskap(jlem)*ztn(:,jlem))            &
     &            /(dsigma(jlev)+(zkdiffh(:,jlev)                       &
     &              +zkdiffh(:,jlem)*(1.-zebs(:,jlem)/zskap(jlem)))     &
     &              /zskap(jlev))
      enddo
!
!     bottom layer elimination
!
      ztn(:,NLEV)=(zt(:,NLEV)*dsigma(NLEV)                              &
     &            +zkdiffh(:,NLEM)*ztn(:,NLEM)/zskap(NLEM))             &
     &           /(dsigma(NLEV)+zkdiffh(:,NLEM)/zskap(NLEV)             &
     &                         *(1.-zebs(:,NLEM)/zskap(NLEM)))
!
!     back-substitution
!
      do jlev=NLEM,1,-1
       jlep=jlev+1
       ztn(:,jlev)=ztn(:,jlev)+zebs(:,jlev)*ztn(:,jlep)/zskap(jlep)
      enddo
!
!     tendencies
!
      zdtdt(:,1:NLEV)=(ztn(:,1:NLEV)-zt(:,1:NLEV))/deltsec2
      dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdt(:,1:NLEV)

!
!     franks dbug
!

      if(ndiaggp==1) then
       dgp3d(:,1:NLEV,4)=zdtdt(:,1:NLEV)
       do jlev=1,NLEV
        dgp3d(:,jlev,14)=zdqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if
!
!     energy diagnostics
!
      if(nenergy > 0) then
       denergy(:,8)=0.
       denergy(:,22)=0.
       do jlev=1,NLEV
        denergy(:,8)=denergy(:,8)                                       &
     &               +zdtdt(:,jlev)                                     &
     &               *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
        denergy(:,22)=denergy(:,22)                                     &
     &               -((zun(:,jlev)*zun(:,jlev)                         &
     &                 -zu(:,jlev)*zu(:,jlev)                           &
     &                 +zvn(:,jlev)*zvn(:,jlev)                         &
     &                 -zv(:,jlev)*zv(:,jlev))/deltsec2                 &
     &                 -(zken(:,jlev)-zke(:,jlev))/deltsec2)            &
     &               *0.5*dp(:)/ga*dsigma(jlev)
        if(nener3d > 0) then
         dener3d(:,jlev,8)=zdtdt(:,jlev)                                &
     &                   *acpd*(1.+adv*dq(:,jlev))*dp(:)/ga*dsigma(jlev)
         dener3d(:,jlev,22)=-((zun(:,jlev)*zun(:,jlev)                  &
     &                        -zu(:,jlev)*zu(:,jlev)                    &
     &                        +zvn(:,jlev)*zvn(:,jlev)                  &
     &                        -zv(:,jlev)*zv(:,jlev))/deltsec2          &
     &                       -(zken(:,jlev)-zke(:,jlev))/deltsec2)      &
     &                     *0.5*dp(:)/ga*dsigma(jlev)
        endif   
       enddo
      endif
!
      return
      end subroutine vdiff
