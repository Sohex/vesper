      program shtns_equivalence
!     Does SHTns plus the conversion recipe compute what legmod computes?
!
!     Worldbuilding frame: a correctness check on the Vesper climate model's
!     spectral transform. Nothing here is about the simulated planet.
!
!     The conventions were measured one mode at a time by
!     `probe_shtns_conventions.f90` and `probe_shtns_vector_conventions.f90`.
!     This is the gate that comes after: the same DENSE field through both
!     paths, with the recipe applied, compared as fields. A per-mode ratio can
!     be right on every mode tried and still leave the recipe wrong -- a missed
!     m=0 special case, an index off by one somewhere in the middle, a factor
!     that only shows when modes interfere. This runs the whole spectrum at once.
!
!     THE RECIPE, measured and not read:
!
!       shtns_create(NTRU, NTRU, 1, SHT_ORTHONORMAL)   ! NOT SHT_NO_CS_PHASE
!       shtns_set_grid(SHT_GAUSS + SHT_PHI_CONTIGUOUS, ...)
!       shtns_robert_form(cfg, 1)
!       PlaSim lm  <->  SHTns lm-1
!       scalar :  model = +sqrt(2 pi) * SHTns
!       vector :  gu = -sqrt(2 pi)/(l(l+1)) * Vp ,  gv = +sqrt(2 pi)/(l(l+1)) * Vt
!
!     THE CONTROL is the mistake this recipe was got wrong by once: dropping the
!     Condon-Shortley phase. It negates every odd zonal wavenumber, which a
!     comparison of magnitudes cannot see and which leaves a model that runs.
!     Built here so it cannot recur quietly.
      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: real64
      use pumamod, only: NLAT, NLON, NLPP, NLEV, NTRU, NTP1, NCSP, NESP,        &
     &                   sid, gwd, plavor
      implicit none
      include 'shtns.f03'

      integer, parameter :: wp = real64
      type(shtns_info), pointer :: shc
      type(c_ptr) :: shp
      integer :: knorm, klay, jj, jm, jn, jlm, nbad, jcase
      character(len=15) :: ycase(3) = ['divergence only','vorticity only ','both together  ']
      real(wp) :: zeps, zsi(NLAT), zgw(NLAT), zroot, ztol, zerr, zscale
      real :: zspc(2,NCSP), zgpf(NLON,NLPP)
      real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
      real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
      real(wp), allocatable :: zg(:,:), zvt(:,:), zvp(:,:)
      complex(wp), allocatable :: zs(:), zt(:)
      integer, allocatable :: kdeg(:)
      logical :: lnocs

      ztol = 1.0e-11_wp
      nbad = 0
      zroot = sqrt(8.0_wp*atan(1.0_wp))          ! sqrt(2 pi)

!     the control is selected by an environment variable so one binary serves
!     both arms and they cannot drift apart
      call getflag(lnocs)

      call inigau(NLAT,zsi,zgw)
      do jj = 1 , NLPP
         sid(jj) = zsi(jj)
         gwd(jj) = zgw(jj)
      enddo
      plavor = 0.0                               ! it would land on mode 2
      call legini

      knorm = SHT_ORTHONORMAL
      if (lnocs) knorm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
      klay = SHT_GAUSS + SHT_PHI_CONTIGUOUS
      zeps = 0.0_wp
      call shtns_verbose(0)
      jj = shtns_use_threads(1)
      shp = shtns_create(NTRU, NTRU, 1, knorm)
      call shtns_set_grid(shp, klay, zeps, NLAT, NLON)
      call c_f_pointer(cptr=shp, fptr=shc)
      call shtns_robert_form(shp, 1)
      allocate( zg(shc%nphi,shc%nlat), zvt(shc%nphi,shc%nlat), zvp(shc%nphi,shc%nlat) )
      allocate( zs(shc%nlm), zt(shc%nlm), kdeg(NCSP) )

      jlm = 0                                    ! the degree l of each mode
      do jm = 0 , NTRU
         do jn = jm , NTRU
            jlm = jlm + 1
            kdeg(jlm) = jn
         enddo
      enddo

      write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT, '  modes ', NCSP
      if (lnocs) then
         write(*,'(a)') 'CONTROL ARM: Condon-Shortley phase dropped. Must FAIL.'
      else
         write(*,'(a,e9.2)') 'tolerance ', ztol
      endif
      write(*,*)

!     ---- scalar ----------------------------------------------------------
      do jj = 1 , NCSP
         zspc(1,jj) =  1.0 / real(jj)
         zspc(2,jj) = -0.5 / real(jj)
      enddo
      zspc(2,1:NTP1) = 0.0                       ! m=0 is real
      call sp2fc(zspc,zgpf)
      call fc2gp(zgpf,NLON,NLPP)
      do jj = 1 , NCSP
         zs(jj) = cmplx(real(zspc(1,jj),wp), real(zspc(2,jj),wp), kind=wp)
      enddo
      call SH_to_spat(shp, zs, zg)
      zerr = maxval(abs(real(zgpf,wp) - zroot*zg)) / maxval(abs(real(zgpf,wp)))
      call verdict('scalar, dense spectrum', zerr, ztol, nbad)

!     ---- vector -----------------------------------------------------------
!     Divergence and vorticity are driven SEPARATELY before they are driven
!     together. The per-mode probes measured the spheroidal convention with the
!     vorticity held at zero and took the toroidal on trust, which is exactly
!     the assumption a dense mixed field breaks. Separating them says WHICH of
!     the two is wrong instead of only that one is.
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

         do jj = 1 , NCSP
            zscale = 0.0_wp
            if (kdeg(jj) > 0) zscale = 1.0_wp / real(kdeg(jj)*(kdeg(jj)+1), wp)
            zs(jj) = cmplx(real(zsd(1,jj,1),wp),real(zsd(2,jj,1),wp),kind=wp)*zscale
!           THE TOROIDAL TAKES THE OPPOSITE SIGN, and it is not a fudge:
!           for V = grad(S) + curl(T r), divergence is the Laplacian of S while
!           vorticity is MINUS the Laplacian of T, so the two potentials sit on
!           opposite sides of their sources. Measured as a relative error of
!           exactly 2.000 on the vorticity-only arm, which is what a pure sign
!           flip looks like.
            zt(jj) = -cmplx(real(zsz(1,jj,1),wp),real(zsz(2,jj,1),wp),kind=wp)*zscale
         enddo
         call SHsphtor_to_spat(shp, zs, zt, zvt, zvp)

!        zfu is declared (2,NLON/2,NLPP,NLEV) and holds grid data after fc2gp,
!        while SHTns hands back (nphi,nlat). Same elements, different shape.
         zerr = maxval(abs(reshape(real(zfu(:,:,:,1),wp),[NLON*NLPP])           &
     &                   + zroot*reshape(zvp,[NLON*NLPP])))                     &
     &        / max(maxval(abs(real(zfu(:,:,:,1),wp))), 1.0e-30_wp)
         call verdict('u, '//ycase(jcase), zerr, ztol, nbad)
         zerr = maxval(abs(reshape(real(zfv(:,:,:,1),wp),[NLON*NLPP])           &
     &                   - zroot*reshape(zvt,[NLON*NLPP])))                     &
     &        / max(maxval(abs(real(zfv(:,:,:,1),wp))), 1.0e-30_wp)
         call verdict('v, '//ycase(jcase), zerr, ztol, nbad)
      enddo

      write(*,*)
      if (lnocs) then
         if (nbad == 0) then
            write(*,'(a)') 'CONTROL PASSED, which means this check cannot see a'
            write(*,'(a)') 'dropped Condon-Shortley phase. That is a failure.'
            stop 1
         else
            write(*,'(a)') 'control rejected, so the check has teeth'
         endif
      else
         if (nbad == 0) then
            write(*,'(a)') '0 failed: SHTns and the recipe compute this model'
         else
            write(*,'(i0,a)') nbad, ' failed'
            stop 1
         endif
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

      subroutine getflag(lflag)
      logical :: lflag
      character(len=8) :: zval
      integer :: kstat
      call get_environment_variable('SHTNS_DROP_CS', zval, status=kstat)
      lflag = (kstat == 0 .and. trim(zval) == '1')
      end subroutine getflag

      end program shtns_equivalence
