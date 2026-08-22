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
     &                   NUGP, sid, gwd, plavor
      use shtnsmod, only: shtns_setup, sh_sp2gp, sh_dv2uv
      implicit none

      integer, parameter :: wp = real64
      integer :: jj, jcase, nbad
      real(wp) :: zsi(NLAT), zgw(NLAT), ztol, zerr
      real :: zspc(2,NCSP), zgpf(NLON,NLPP)
      real :: zsp1(NESP,1), zgp1(NUGP,1)
      real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
      real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
      real :: zsdf(NESP,1), zszf(NESP,1), zgu(NUGP,1), zgv(NUGP,1)
      character(len=15) :: ycase(3) = ['divergence only','vorticity only ','both together  ']

      ztol = 1.0e-11_wp
      nbad = 0

      call inigau(NLAT,zsi,zgw)
      do jj = 1 , NLPP
         sid(jj) = zsi(jj)
         gwd(jj) = zgw(jj)
      enddo
      plavor = 0.0                     ! dv2uv would otherwise add it to mode 2
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

      end program shtns_equivalence
