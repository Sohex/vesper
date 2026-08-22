      program analysismargin
!     Why do the SHTns ANALYSIS arms miss their bar above T42 while synthesis is
!     at rounding scale everywhere? CLIM-61. Two hypotheses, one probe.
!
!     A. GAUSS QUADRATURE MARGIN. The grid satisfies NLAT >= (3T+1)/2, which is
!        the exactness condition for a SCALAR quadratic product. A vector
!        transform integrates products involving dP/dmu and m/(1-mu^2), whose
!        effective polynomial degree is higher, so at the scalar margin the
!        vector integrand is UNDER-INTEGRATED. legmod absorbs that with
!        derivative weight matrices and gwdc; SHTns uses analytic vector
!        harmonic relations. The residual aliasing differs between them, and the
!        difference would show only in the vector arms -- which is what happens.
!
!        THE TEST: band-limit the field well below the truncation. If
!        under-integration at high n is the mechanism, the disagreement collapses
!        as the cut comes down, and it collapses for the VECTOR arms while the
!        scalar arms barely move.
!
!     B. FFT RADIX. NLON is 2*NLAT, so T85 gives 256 and T170 gives 512, both
!        powers of two, while T127 gives 384 = 3*2^7 and pulls in radix-3
!        butterflies with different cancellation.
!
!        THE TEST: one SHTns configuration on the model's own nphi and a second
!        on the next power of two at the SAME nlat, round trip through both, and
!        compare. legmod is not involved, so nothing but the FFT differs.
!
!     The two are not exclusive and the measured pattern fits neither alone:
!     T85 is a power of two and still misses the bar, so radix cannot be the
!     whole story, and every resolution here sits at the scalar margin, so the
!     margin alone does not explain why T170 is six times better than T127.
      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: real64
      use pumamod, only: NLAT, NLON, NLPP, NLEV, NTRU, NTP1, NCSP, NESP,  &
     &                   NUGP, sid, gwd, plavor, EZ, nfilter, ngptfilter, &
     &                   nspvfilter, filterkappa, nfilterexp
      use shtnsmod, only: shtns_setup, sh_gp2sp, sh_uv2dv, shtcfg
      implicit none
      include 'shtns.f03'

      integer, parameter :: wp = real64
      integer :: jj, jm, jn, jlm, jc, kcut, kphi
      real(wp) :: zsi(NLAT), zgw(NLAT)
      real :: zspc(2,NCSP), zg(NLON,NLPP), zwork(NLON,NLPP)
      real :: zsl(2,NESP/2), zhs(NESP,1)
      real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
      real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
      real :: zrd(2,NESP/2,NLEV), zrz(2,NESP/2,NLEV)
      real :: zgu(NUGP,1), zgv(NUGP,1), zhd(NESP,1), zhz(NESP,1)
      real :: zflat(NUGP)
      integer :: kcuts(4)

      type(c_ptr) :: cfg2
      real(wp), allocatable :: za(:), zb(:)
      complex(wp), allocatable :: zl1(:), zl2(:)
      integer :: nphi2, knorm, klay
      real(wp) :: zeps, e1, e2

      call inigau(NLAT,zsi,zgw)
      do jj = 1 , NLPP
         sid(jj) = zsi(jj)
         gwd(jj) = zgw(jj)
      enddo
      nfilter     = 2
      ngptfilter  = 1
      nspvfilter  = 1
      filterkappa = 8.0
      nfilterexp  = 8
      plavor      = EZ
      call legini
      call shtns_setup

      write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT, '  NLON ', NLON
      write(*,'(a,f8.2,a,i0)') 'scalar exactness needs NLAT >= ',            &
     &      (3.0*NTRU+1.0)/2.0, ', grid has ', NLAT
      write(*,*)

!     ---- A: does the disagreement come from the top of the spectrum? --------
      kcuts = [NTRU, (NTRU*3)/4, NTRU/2, NTRU/4]
      write(*,'(a)') 'A. band-limiting the test field'
      write(*,'(a8,a18,a18)') 'n cut', 'scalar analysis', 'uv2dv divergence'
      do jc = 1 , 4
         kcut = kcuts(jc)
         zspc(:,:) = 0.0
         jlm = 0
         do jm = 0 , NTRU
            do jn = jm , NTRU
               jlm = jlm + 1
               if (jn > kcut) cycle
               zspc(1,jlm) =  1.0 / real(jlm)
               if (jm > 0) zspc(2,jlm) = -0.5 / real(jlm)
            enddo
         enddo

!        scalar: legmod synthesis, then both analyses
         call sp2fc(zspc,zg)
         call fc2gp(zg,NLON,NLPP)
         zwork = zg
         call gp2fc(zwork,NLON,NLPP)
         call fc2sp(zwork,zsl)
         zflat = reshape(zg,[NUGP])
         call sh_gp2sp(zflat, zhs, 1)

!        vector: a wind from this field as both divergence and vorticity
         zsd(:,:,:) = 0.0 ; zsz(:,:,:) = 0.0
         do jj = 1 , NCSP
            zsd(1,jj,1) = zspc(1,jj) ; zsd(2,jj,1) = zspc(2,jj)
            zsz(1,jj,1) = 0.7*zspc(1,jj) ; zsz(2,jj,1) = 0.3*zspc(2,jj)
         enddo
         plavor = 0.0
         call dv2uv(zsd,zsz,zfu,zfv)
         call fc2gp(zfu(1,1,1,1),NLON,NLPP)
         call fc2gp(zfv(1,1,1,1),NLON,NLPP)
         plavor = EZ
         zgu(:,1) = reshape(zfu(:,:,:,1),[NUGP])
         zgv(:,1) = reshape(zfv(:,:,:,1),[NUGP])
         zfu(:,:,:,1) = reshape(zgu(:,1),[2,NLON/2,NLPP])
         zfv(:,:,:,1) = reshape(zgv(:,1),[2,NLON/2,NLPP])
         call gp2fc(zfu(1,1,1,1),NLON,NLPP)
         call gp2fc(zfv(1,1,1,1),NLON,NLPP)
         call uv2dv(zfu,zfv,zrd,zrz)
         call sh_uv2dv(zgu, zgv, zhd, zhz, 1)

         write(*,'(i8,e18.4,e18.4)') kcut,                                   &
     &      relerr(zsl, zhs), relerr(zrd(1,1,1), zhd)
      enddo

!     ---- B: the FFT length, with legmod out of the picture ------------------
      knorm = SHT_ORTHONORMAL
      klay  = SHT_QUICK_INIT + SHT_PHI_CONTIGUOUS
      zeps  = 0.0_wp
      nphi2 = 1
      do while (nphi2 < NLON)
         nphi2 = nphi2 * 2
      enddo
      if (nphi2 == NLON) nphi2 = NLON * 2      ! already a power of two

      write(*,*)
      write(*,'(a)') 'B. SHTns round trip, same nlat, two FFT lengths'
      e1 = roundtrip(NLON)
      e2 = roundtrip(nphi2)
      write(*,'(a,i5,a,e12.4)') '   nphi ', NLON,  '  round trip ', e1
      write(*,'(a,i5,a,e12.4)') '   nphi ', nphi2, '  round trip ', e2
      write(*,'(a,f10.3)')      '   ratio (model/power-of-two) ', e1 / max(e2,1.0e-300_wp)

!     ---- C: the quadrature weights themselves -------------------------------
!     Weights enter the ANALYSIS and not the synthesis: legmod carries gwd
!     explicitly, SHTns applies its own. The synthesis arms agree to 3e-14, so
!     the NODES agree; nothing so far has checked the WEIGHTS. A field-blind,
!     band-limit-blind error that appears only on analysis is what a weight
!     mismatch would look like.
      block
        real(wp) :: zw(NLAT), zr, zrmin, zrmax
        integer :: k, nw
        call shtns_gauss_wts(shtcfg, zw)
        write(*,*)
        write(*,'(a)') 'C. quadrature weights, legmod gwd against SHTns'
!       shtns_gauss_wts fills the northern half only.
        nw = NLAT/2
        zrmin = 1.0e30_wp ; zrmax = -1.0e30_wp
        do k = 1 , nw
           if (zw(k) /= 0.0_wp) then
              zr = real(gwd(k),wp) / zw(k)
              zrmin = min(zrmin, zr) ; zrmax = max(zrmax, zr)
           endif
        enddo
        write(*,'(a,i5)')      '   weights compared      ', nw
        write(*,'(a,e20.12)')  '   ratio gwd/shtns min   ', zrmin
        write(*,'(a,e20.12)')  '   ratio gwd/shtns max   ', zrmax
        write(*,'(a,e12.4)')   '   spread of the ratio   ',                  &
     &     (zrmax - zrmin) / max(abs(zrmax), 1.0e-300_wp)
        write(*,'(a)') '   a CONSTANT ratio is a normalisation and harmless;'
        write(*,'(a)') '   a spread is the weights actually disagreeing.'
!       Printed at full precision so a third, independent implementation can
!       say which of the two is right.
        write(*,*)
        write(*,'(a)') '   lat        gwd (inigau)              shtns'
        do k = 1 , 5
           write(*,'(i6,2e26.17)') k, real(gwd(k),wp), zw(k)
        enddo
        do k = nw-2 , nw
           write(*,'(i6,2e26.17)') k, real(gwd(k),wp), zw(k)
        enddo
      end block

      contains

      function relerr(pref, pgot) result(z)
!     Max absolute difference over the modes, against the largest coefficient,
!     excluding the m=0 imaginary slots that neither side defines.
      real, intent(in) :: pref(2*NCSP), pgot(2*NCSP)
      real(wp) :: z, zden
      integer :: k
      z = 0.0_wp ; zden = 0.0_wp
      do k = 1 , NCSP
         zden = max(zden, abs(real(pref(2*k-1),wp)))
         z    = max(z, abs(real(pref(2*k-1),wp) - real(pgot(2*k-1),wp)))
         if (k > NTP1) then
            zden = max(zden, abs(real(pref(2*k),wp)))
            z    = max(z, abs(real(pref(2*k),wp) - real(pgot(2*k),wp)))
         endif
      enddo
      if (zden > 0.0_wp) z = z / zden
      end function relerr

      function roundtrip(kphi) result(z)
!     spat_to_SH(SH_to_spat(x)) against x, on a fresh configuration whose only
!     difference from the model's is the FFT length.
      integer, intent(in) :: kphi
      real(wp) :: z, zden
      type(c_ptr) :: cfg
      type(shtns_info), pointer :: inf
      complex(wp), allocatable :: zin(:), zout(:)
      real(wp), allocatable :: zsp(:)
      integer :: k

      call shtns_verbose(0)
      k = shtns_use_threads(1)
      cfg = shtns_create(NTRU, NTRU, 1, knorm)
      call shtns_set_grid(cfg, klay, zeps, NLAT, kphi)
      call c_f_pointer(cptr=cfg, fptr=inf)
      allocate(zin(inf%nlm), zout(inf%nlm), zsp(inf%nphi*inf%nlat))
      do k = 1 , inf%nlm
         zin(k) = cmplx(1.0_wp/real(k,wp), -0.5_wp/real(k,wp), kind=wp)
      enddo
!     m=0 has no imaginary part and the round trip does not preserve one, so
!     leaving it in measures the test's own construction and not the FFT.
      do k = 1 , NTP1
         zin(k) = cmplx(real(zin(k)), 0.0_wp, kind=wp)
      enddo
      call SH_to_spat(cfg, zin, zsp)
      call spat_to_SH(cfg, zsp, zout)
      z = 0.0_wp ; zden = 0.0_wp
      do k = 1 , inf%nlm
         zden = max(zden, abs(zin(k)))
         z    = max(z, abs(zin(k) - zout(k)))
      enddo
      if (zden > 0.0_wp) z = z / zden
      deallocate(zin, zout, zsp)
      call shtns_destroy(cfg)
      end function roundtrip

      end program analysismargin
