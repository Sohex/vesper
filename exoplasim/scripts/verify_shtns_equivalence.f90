      program shtns_equivalence
!     Does SHTns, as the model actually calls it, compute what legmod computes?
!
!     Worldbuilding frame: a correctness check on the Vesper climate model's
!     spectral transform. Nothing here is about the simulated planet.
!
!     THIS CALLS THE SHIPPED WRAPPERS, `sh_sp2gp` and `sh_dv2uv` in shtnsmod,
!     rather than SHTns directly. An earlier version held its own copy of the
!     conversion, which made it a check that two copies of a recipe agreed with
!     each other -- and shtnsmod exists precisely because a conversion written
!     twice is a conversion that will be wrong once. What ships is what is
!     tested.
!
!     Dense fields rather than one mode at a time. The per-mode probes measured
!     the conventions and MISSED the toroidal sign, because they drove
!     divergence alone; a dense mixed field is what caught it. Divergence and
!     vorticity are still driven separately BEFORE they are driven together, so
!     a failure says which of the two is wrong rather than only that one is.
!
!     THE CONTROL is applied by the shell to shtnsmod itself, not selected here:
!     it drops the Condon-Shortley phase, which negates every odd zonal
!     wavenumber and which no comparison of magnitudes can see. That is the
!     mistake this recipe was got wrong by once.
      use, intrinsic :: iso_fortran_env, only: real64
      use pumamod, only: NLAT, NLON, NLPP, NLEV, NTRU, NTP1, NCSP, NESP,        &
     &                   NUGP, sid, gwd, plavor, EZ, nfilter, ngptfilter,           &
     &                   nspvfilter, filterkappa, nfilterexp
      use shtnsmod, only: shtns_setup, sh_sp2gp, sh_dv2uv, sh_sp2grad,        &
     &                    sh_gp2sp, sh_uv2dv, sh_advtend, sh_dztend
      implicit none

      integer, parameter :: wp = real64
      integer :: jj, jcase, nbad
      real(wp) :: zsi(NLAT), zgw(NLAT), ztol, zerr
      real :: zspc(2,NCSP), zgpf(NLON,NLPP)
      real :: zsp1(NESP,1), zgp1(NUGP,1)
      real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
      real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
      real :: zsdf(NESP,1), zszf(NESP,1), zgu(NUGP,1), zgv(NUGP,1)
      real :: zgpj(NLON,NLPP), zgpm(NLON,NLPP), zgp(NLON,NLPP)
      real :: zdmu(NUGP,1), zdlam(NUGP,1)
!     the analysis direction
      real :: zga(NUGP,1), zgb(NUGP,1), zgc(NUGP,1)
      real :: zsl(2,NESP/2)
      real :: zfa(2,NLON/2,NLPP,NLEV), zfb(2,NLON/2,NLPP,NLEV)
      real :: zfc(2,NLON/2,NLPP,NLEV)
      real :: zrd(2,NESP/2,NLEV), zrz(2,NESP/2,NLEV), zrt(2,NESP/2,NLEV)
      real :: zhd(NESP,1), zhz(NESP,1), zht(NESP,1)
      integer :: jk
      integer :: jlon, jlat
      real :: zrm
      character(len=15) :: ycase(3) = ['divergence only','vorticity only ','both together  ']

      ztol = 1.0e-11_wp
      nbad = 0

      call inigau(NLAT,zsi,zgw)
      do jj = 1 , NLPP
         sid(jj) = zsi(jj)
         gwd(jj) = zgw(jj)
      enddo
!     THE CONFIGURATION UNDER TEST IS THE ONE THE MODEL RUNS, and both halves
!     of that sentence were learned the hard way. This driver used to set
!     plavor to zero and leave nfilter at its default of none, and it passed at
!     5e-14 while the model it was certifying disagreed by 100% at the first
!     step -- because a spectral filter and a rotating planet are exactly the
!     two things the wrappers were missing, and neither was switched on here.
!     A check may only simplify what it does not certify.
!
!     filterkappa and nfilterexp are `filter_kappa` and `filter_power` in
!     config/planet.yaml and must equal them; scripts/smoke_test.py compares the
!     two and fails on drift, because this pair went stale once and the gate
!     said nothing. The strength matters to the VERDICT and not only to the
!     realism of the setup: both arms carry the same skspgp(n), so the filter
!     divides out of a per-mode ratio but not out of this gate's field norm,
!     where it reweights the residual's spectrum against a denominator the low
!     modes own. exp(-kappa*x^16)/exp(-kappa*x^8) peaks at exp(kappa/4) -- 7.39
!     at kappa 8, at n/NTRU = (1/2)**(1/8) = 0.917 -- so a gate left at gamma 8
!     passes the mid-to-high band, where the SHTns residual is largest, at up to
!     a seventh of the amplitude the shipped configuration gives it.
      nfilter     = 2                  ! exponential, as the beds run it
      ngptfilter  = 1
      nspvfilter  = 1
      filterkappa = 8.0
      nfilterexp  = 16
      plavor      = EZ                 ! rotating, as the model runs it
      call legini
      call shtns_setup

      write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT, '  modes ', NCSP
      write(*,'(a,e9.2)') 'tolerance ', ztol
      write(*,*)

!     ---- scalar ------------------------------------------------------------
      do jj = 1 , NCSP
         zspc(1,jj) =  1.0 / real(jj)
         zspc(2,jj) = -0.5 / real(jj)
      enddo
      zspc(2,1:NTP1) = 0.0             ! m=0 is real; the transform discards an
                                       ! imaginary part there and is right to
      call sp2fc(zspc,zgpf)
      call fc2gp(zgpf,NLON,NLPP)
      zsp1(:,1) = 0.0
      do jj = 1 , NCSP
         zsp1(2*jj-1,1) = zspc(1,jj)
         zsp1(2*jj  ,1) = zspc(2,jj)
      enddo
      call sh_sp2gp(zsp1, zgp1, 1)
      zerr = maxval(abs(reshape(real(zgpf,wp),[NUGP]) - real(zgp1(:,1),wp)))    &
     &     / maxval(abs(real(zgpf,wp)))
      call verdict('scalar, dense spectrum', zerr, ztol, nbad)

!     ---- vector ------------------------------------------------------------
      do jcase = 1 , 3
         zsd(:,:,:) = 0.0 ; zsz(:,:,:) = 0.0
         do jj = 1 , NCSP
            if (jcase /= 2) then
               zsd(1,jj,1) =  1.0 / real(jj)
               zsd(2,jj,1) = -0.5 / real(jj)
            endif
            if (jcase /= 1) then
               zsz(1,jj,1) =  0.7 / real(jj)
               zsz(2,jj,1) =  0.3 / real(jj)
            endif
         enddo
         zsd(2,1:NTP1,1) = 0.0 ; zsz(2,1:NTP1,1) = 0.0

         call dv2uv(zsd,zsz,zfu,zfv)
         call fc2gp(zfu(1,1,1,1),NLON,NLPP)
         call fc2gp(zfv(1,1,1,1),NLON,NLPP)

         zsdf(:,1) = 0.0 ; zszf(:,1) = 0.0
         do jj = 1 , NCSP
            zsdf(2*jj-1,1) = zsd(1,jj,1) ; zsdf(2*jj,1) = zsd(2,jj,1)
            zszf(2*jj-1,1) = zsz(1,jj,1) ; zszf(2*jj,1) = zsz(2,jj,1)
         enddo
         call sh_dv2uv(zsdf, zszf, zgu, zgv, 1)

         zerr = maxval(abs(reshape(real(zfu(:,:,:,1),wp),[NUGP])                &
     &                   - real(zgu(:,1),wp)))                                  &
     &        / max(maxval(abs(real(zfu(:,:,:,1),wp))), 1.0e-30_wp)
         call verdict('u, '//ycase(jcase), zerr, ztol, nbad)
         zerr = maxval(abs(reshape(real(zfv(:,:,:,1),wp),[NUGP])                &
     &                   - real(zgv(:,1),wp)))                                  &
     &        / max(maxval(abs(real(zfv(:,:,:,1),wp))), 1.0e-30_wp)
         call verdict('v, '//ycase(jcase), zerr, ztol, nbad)
      enddo

!     ---- the two pressure derivatives ---------------------------------------
!     gpj as sp2fcdmu gives it, and gpmt as gridpointa builds it by hand from
!     gp's FOURIER coefficients by multiplying each by i*m -- the construction
!     that cannot survive a transform landing in grid space.
      call sp2fcdmu(zspc,zgpj)
      call fc2gp(zgpj,NLON,NLPP)
      call sp2fc(zspc,zgp)
      do jlat = 1 , NLPP
         do jlon = 1 , NLON-1 , 2
            zrm = real((jlon-1)/2)
            zgpm(jlon  ,jlat) = -zgp(jlon+1,jlat) * zrm
            zgpm(jlon+1,jlat) =  zgp(jlon  ,jlat) * zrm
         enddo
      enddo
      call fc2gp(zgpm,NLON,NLPP)

      call sh_sp2grad(zsp1, zdmu, zdlam, 1)
      zerr = maxval(abs(reshape(real(zgpj,wp),[NUGP]) - real(zdmu(:,1),wp)))    &
     &     / maxval(abs(real(zgpj,wp)))
      call verdict('gpj, the mu derivative ', zerr, ztol, nbad)
      zerr = maxval(abs(reshape(real(zgpm,wp),[NUGP]) - real(zdlam(:,1),wp)))   &
     &     / maxval(abs(real(zgpm,wp)))
      call verdict('gpmt, the zonal one    ', zerr, ztol, nbad)


!     ---- the analysis direction ---------------------------------------------
!     Three band-limited grid fields, made through the synthesis path that the
!     arms above have just certified, so the comparison measures the analysis
!     and not aliasing: a forward transform is a quadrature, and an arbitrary
!     grid field carries power above the truncation that the two sides fold
!     back differently.
      call mkgrid(1.0, -0.5, zga)
      call mkgrid(0.6,  0.2, zgb)
      call mkgrid(-0.3, 0.8, zgc)

!     scalar: gp2fc then fc2sp
      zfa(:,:,:,1) = reshape(zga(:,1),[2,NLON/2,NLPP])
      call gp2fc(zfa(1,1,1,1),NLON,NLPP)
      call fc2sp(zfa(1,1,1,1),zsl)
      call sh_gp2sp(zga, zhd, 1)
      call spcheck('gp2sp, the scalar analysis ', zsl, zhd, ztol, nbad)

!     vector: gp2fc on both then uv2dv
      zfa(:,:,:,1) = reshape(zga(:,1),[2,NLON/2,NLPP])
      zfb(:,:,:,1) = reshape(zgb(:,1),[2,NLON/2,NLPP])
      call gp2fc(zfa(1,1,1,1),NLON,NLPP)
      call gp2fc(zfb(1,1,1,1),NLON,NLPP)
      call uv2dv(zfa,zfb,zrd,zrz)
      call sh_uv2dv(zga, zgb, zhd, zhz, 1)
      call spcheck('uv2dv, divergence          ', zrd(1,1,1), zhd, ztol, nbad)
      call spcheck('uv2dv, vorticity           ', zrz(1,1,1), zhz, ztol, nbad)

!     qtend: minus the divergence of (uq,vq), plus the analysis of qn
      zfa(:,:,:,1) = reshape(zga(:,1),[2,NLON/2,NLPP])
      zfb(:,:,:,1) = reshape(zgb(:,1),[2,NLON/2,NLPP])
      zfc(:,:,:,1) = reshape(zgc(:,1),[2,NLON/2,NLPP])
      call gp2fc(zfa(1,1,1,1),NLON,NLPP)
      call gp2fc(zfb(1,1,1,1),NLON,NLPP)
      call gp2fc(zfc(1,1,1,1),NLON,NLPP)
      call qtend(zrt,zfc,zfa,zfb)
      call sh_advtend(zgc, zga, zgb, zht, 1)
      call spcheck('qtend, the advective one   ', zrt(1,1,1), zht, ztol, nbad)

!     mktend: d and z from (fu,fv) and the kinetic energy, t as qtend
      call mktend(zrd,zrt,zrz,zfc,zfa,zfb,zfc,zfa,zfb)
      call sh_dztend(zga, zgb, zgc, zhd, zhz, 1)
      call spcheck('mktend, divergence         ', zrd(1,1,1), zhd, ztol, nbad)
      call spcheck('mktend, vorticity          ', zrz(1,1,1), zhz, ztol, nbad)
      call sh_advtend(zgc, zga, zgb, zht, 1)
      call spcheck('mktend, temperature        ', zrt(1,1,1), zht, ztol, nbad)

      write(*,*)
      if (nbad == 0) then
         write(*,'(a)') '0 failed: the shipped wrappers compute this model'
      else
         write(*,'(i0,a)') nbad, ' failed'
         stop 1
      endif

      contains

      subroutine verdict(what, perr, ptol, kbad)
      character(len=*) :: what
      real (kind=8) :: perr, ptol
      integer :: kbad
      if (perr <= ptol) then
         write(*,'(a,a,a,e11.4)') '  [  ok  ] ', what, '  rel ', perr
      else
         write(*,'(a,a,a,e11.4)') '  [ FAIL ] ', what, '  rel ', perr
         kbad = kbad + 1
      endif
      end subroutine verdict


      subroutine mkgrid(pa, pb, pgrid)
!     A band-limited grid field from a dense spectral one, through the
!     synthesis path the arms above certify.
      real, intent(in)  :: pa, pb
      real, intent(out) :: pgrid(NUGP,1)
      real :: zs(NESP,1)
      integer :: jj
      zs(:,1) = 0.0
      do jj = 1 , NCSP
         zs(2*jj-1,1) = pa / real(jj)
         zs(2*jj  ,1) = pb / real(jj)
      enddo
      do jj = 1 , NTP1
         zs(2*jj,1) = 0.0
      enddo
!     Built with LEGMOD's synthesis, not SHTns's. The field is then exactly
!     band-limited in the representation legmod's own quadrature inverts, so the
!     analysis arms measure the analysis rather than whatever the other
!     library's synthesis left outside legmod's band.
      block
        real :: zf(NLON,NLPP)
        real :: zc(2,NCSP)
        integer :: k
        do k = 1 , NCSP
           zc(1,k) = zs(2*k-1,1) ; zc(2,k) = zs(2*k,1)
        enddo
        call sp2fc(zc,zf)
        call fc2gp(zf,NLON,NLPP)
        pgrid(:,1) = reshape(zf,[NUGP])
      end block
      end subroutine mkgrid

      subroutine spcheck(yname, pref, pgot, ptol, kbad)
!     Compare two packed spectral fields, EXCLUDING the m=0 imaginary slots.
!     A zonal mean has no imaginary part; legmod leaves whatever gp2fc put
!     there and nothing reads it, and the SHTns path sets it to zero. Comparing
!     them would measure which garbage each side happens to carry.
      character(len=*), intent(in) :: yname
      real, intent(in) :: pref(2*NCSP), pgot(2*NCSP)
      real(wp), intent(in) :: ptol
      integer, intent(inout) :: kbad
      real(wp) :: zerr, zden
      integer :: jj
      zerr = 0.0_wp
      zden = 0.0_wp
      do jj = 1 , NCSP
         zden = max(zden, abs(real(pref(2*jj-1),wp)))
         zerr = max(zerr, abs(real(pref(2*jj-1),wp) - real(pgot(2*jj-1),wp)))
         if (jj > NTP1) then
            zden = max(zden, abs(real(pref(2*jj),wp)))
            zerr = max(zerr, abs(real(pref(2*jj),wp) - real(pgot(2*jj),wp)))
         endif
      enddo
      if (zden > 0.0_wp) zerr = zerr / zden
      call verdict(yname, zerr, ptol, kbad)
      end subroutine spcheck

      end program shtns_equivalence
