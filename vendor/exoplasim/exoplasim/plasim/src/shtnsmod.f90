!     ==================================================================
!     shtnsmod.f90
!     ------------
!     The SHTns spherical harmonic transform, wrapped for this model.
!
!     WHY A MODULE OF ITS OWN. It is the only C interoperation in the
!     model, and confining `iso_c_binding` and the shtns.f03 interface to
!     one file keeps that from spreading through legmod and plasim. It
!     also puts the CONVERSION in one place: the two libraries do not
!     agree on normalisation, phase or the sign of the toroidal
!     potential, and a conversion applied in twelve call sites is a
!     conversion that will be wrong in one of them.
!
!     THE RECIPE IS MEASURED, not read. Every constant here was obtained
!     by driving the same field through both transforms and dividing,
!     and `verify_shtns_equivalence.sh` re-derives all of it on dense
!     fields with a control. `exoplasim/notes/shtns-viability.md` has the
!     arithmetic. The three that are easy to get wrong:
!
!       the Condon-Shortley phase is INCLUDED -- SHT_NO_CS_PHASE would
!         negate every odd zonal wavenumber, which no comparison of
!         magnitudes can see;
!       Robert form is ON, because gu is the zonal wind times cos(phi);
!       the toroidal potential takes the OPPOSITE sign to the
!         spheroidal, because divergence is the Laplacian of S while
!         vorticity is minus the Laplacian of T.
!
!     WHAT IT REQUIRES. SHTns needs every latitude in one address space,
!     so it is available only to the threaded build, and it needs the
!     grid in LATITUDE order, so it is incompatible with LPAIRLAT --
!     which costs nothing, since LPAIRLAT exists to let legmod fold a
!     mirror pair together and SHTns replaces legmod. Both are checked
!     rather than assumed, in shtns_setup.
!     ==================================================================

      module shtnsmod
      use iso_c_binding
      use pumamod, only: NTRU, NLAT, NLON, NCSP, NLEV, NLPP, NPRO,      &
     &                   LPAIRLAT, nud, mypid, NROOT
      implicit none
      include 'shtns.f03'

!     One configuration for the whole team. Measured safe to call from
!     several threads at once on one config -- 1000 concurrent transforms
!     bit identical to the serial answer -- which is what lets the model
!     keep its own parallelism and call SHTns single-threaded.
!     `probe_shtns_concurrency.f90` is that measurement.
      type(c_ptr) :: shtcfg = c_null_ptr
      logical     :: lshtns = .false.     ! is the configuration usable

!     sqrt(2 pi). PlaSim's harmonics carry no 1/sqrt(2 pi) in phi, so its
!     grid is this times SHTns's for the same coefficient.
      real (kind=8), parameter :: SHTROOT = 2.5066282746310002_8

!     1/(l(l+1)) by MODE, zero at l=0. The factor between SHTns's
!     spheroidal and toroidal potentials and this model's divergence and
!     vorticity, and the one a naive substitution applies twice, since
!     legini already carries it inside fmu and fmv.
      real (kind=8), allocatable :: shtinv(:)

      contains

      subroutine shtns_setup
      integer :: knorm, klay, jm, jn, jlm, kok
      real (kind=8) :: zeps

      lshtns = .false.

!     Both of these are refusals rather than warnings: a model that
!     silently ran the wrong transform would be the expensive failure.
      if (NPRO > 1 .and. LPAIRLAT) then
         if (mypid == NROOT) then
            write(nud,*) '*** SHTns needs the grid in LATITUDE order and'
            write(nud,*) '*** LPAIRLAT permutes it. Rebuild without the'
            write(nud,*) '*** paired decomposition: SHTns replaces the'
            write(nud,*) '*** transform LPAIRLAT exists to accelerate.'
         endif
         call mpabort('nshtns with LPAIRLAT')
      endif

      knorm = SHT_ORTHONORMAL          ! Condon-Shortley INCLUDED
      klay  = SHT_GAUSS + SHT_PHI_CONTIGUOUS
      zeps  = 0.0_8                    ! no polar truncation: it is worth
                                       ! 7% and this model does not make
                                       ! that approximation elsewhere
      call shtns_verbose(0)
      kok = shtns_use_threads(1)       ! the model's team is the parallelism
      shtcfg = shtns_create(NTRU, NTRU, 1, knorm)
      call shtns_set_grid(shtcfg, klay, zeps, NLAT, NLON)
      call shtns_robert_form(shtcfg, 1)

      allocate(shtinv(NCSP))
      jlm = 0
      do jm = 0 , NTRU
         do jn = jm , NTRU
            jlm = jlm + 1
            if (jn > 0) then
               shtinv(jlm) = 1.0_8 / real(jn*(jn+1), 8)
            else
               shtinv(jlm) = 0.0_8
            endif
         enddo
      enddo

      lshtns = .true.
      if (mypid == NROOT) then
         write(nud,*) 'SHTns configured: T', NTRU, ' nlat', NLAT,        &
     &                ' nlon', NLON, ' modes', NCSP
      endif
      return
      end subroutine shtns_setup


      subroutine sh_sp2gp(psp, pgp, klev)
!     Spectral to GRID for klev levels of one scalar field, the whole globe.
!
!     Replaces sp2fc followed by fc2gp: SHTns does the Legendre transform and
!     the Fourier transform in one call, so there is no Fourier intermediate
!     and nothing downstream should look for one.
!
!     PARALLEL OVER LEVELS, with SHTns called single-threaded. The model's team
!     already exists and owns the parallelism; a thread takes whole levels and
!     writes disjoint columns of pgp, so the levels need no coordination beyond
!     the barrier that ends the worksharing. Calling one config from several
!     threads at once is measured safe -- probe_shtns_concurrency.f90.
      use pumamod, only: NESP, NUGP, NCSP
      integer, intent(in) :: klev
      real, intent(in)    :: psp(NESP,klev)      ! (2,NCSP) packed per level
      real, intent(out)   :: pgp(NUGP,klev)
      complex (kind=8) :: zlm(NCSP)
      real (kind=8) :: zg(NUGP)
      integer :: jlev, jm

!$omp do schedule(static)
      do jlev = 1 , klev
         do jm = 1 , NCSP
            zlm(jm) = cmplx(real(psp(2*jm-1,jlev),8),                   &
     &                      real(psp(2*jm  ,jlev),8), kind=8)
         enddo
         call SH_to_spat(shtcfg, zlm, zg)
         pgp(:,jlev) = real(SHTROOT * zg)
      enddo
!$omp end do
      return
      end subroutine sh_sp2gp


      subroutine sh_dv2uv(psd, psz, pgu, pgv, klev)
!     Divergence and vorticity to wind, the whole globe, for klev levels.
!
!     Replaces dv2uv followed by fc2gp on both components. The conversion is
!     the measured one and every part of it matters:
!
!       the potentials carry 1/(l(l+1)), which is the factor legini keeps
!         inside fmu and fmv and which a naive substitution applies twice;
!       the TOROIDAL takes the opposite sign to the spheroidal, because
!         divergence is the Laplacian of S and vorticity is MINUS that of T;
!       u takes a further minus, and both carry sqrt(2 pi).
!
!     Robert form is on, so SHTns returns the wind already multiplied by
!     cos(phi), which is what gu and gv are.
      use pumamod, only: NESP, NUGP, NCSP, NLEV
      integer, intent(in) :: klev
      real, intent(in)    :: psd(NESP,klev), psz(NESP,klev)
      real, intent(out)   :: pgu(NUGP,klev), pgv(NUGP,klev)
      complex (kind=8) :: zs(NCSP), zt(NCSP)
      real (kind=8) :: zvt(NUGP), zvp(NUGP)
      integer :: jlev, jm

!$omp do schedule(static)
      do jlev = 1 , klev
         do jm = 1 , NCSP
            zs(jm) =  cmplx(real(psd(2*jm-1,jlev),8),                   &
     &                      real(psd(2*jm  ,jlev),8), kind=8) * shtinv(jm)
            zt(jm) = -cmplx(real(psz(2*jm-1,jlev),8),                   &
     &                      real(psz(2*jm  ,jlev),8), kind=8) * shtinv(jm)
         enddo
         call SHsphtor_to_spat(shtcfg, zs, zt, zvt, zvp)
         pgu(:,jlev) = real(-SHTROOT * zvp)
         pgv(:,jlev) = real( SHTROOT * zvt)
      enddo
!$omp end do
      return
      end subroutine sh_dv2uv


      subroutine sh_sp2grad(psp, pgdmu, pgdlam, klev)
!     Both horizontal derivatives of a scalar, on the grid, in one call.
!
!     Replaces sp2fcdmu for the meridional derivative AND the model's own zonal
!     one, which gridpointa builds by hand from the FOURIER coefficients of gp
!     by multiplying each by i*m. That construction cannot survive a transform
!     that lands in grid space, and it does not have to: SHsph_to_spat returns
!     both components of the gradient together.
!
!     No 1/(l(l+1)) here, unlike sh_dv2uv. SHsph_to_spat takes the scalar
!     ITSELF and differentiates it, rather than a potential whose Laplacian is
!     the source, and the measurement confirms it -- the ratio is a pure
!     sqrt(2 pi) with no dependence on l.
      use pumamod, only: NESP, NUGP, NCSP
      integer, intent(in) :: klev
      real, intent(in)    :: psp(NESP,klev)
      real, intent(out)   :: pgdmu(NUGP,klev)     ! as sp2fcdmu gives it
      real, intent(out)   :: pgdlam(NUGP,klev)    ! the zonal derivative
      complex (kind=8) :: zlm(NCSP)
      real (kind=8) :: zvt(NUGP), zvp(NUGP)
      integer :: jlev, jm

!$omp do schedule(static)
      do jlev = 1 , klev
         do jm = 1 , NCSP
            zlm(jm) = cmplx(real(psp(2*jm-1,jlev),8),                   &
     &                      real(psp(2*jm  ,jlev),8), kind=8)
         enddo
         call SHsph_to_spat(shtcfg, zlm, zvt, zvp)
         pgdmu(:,jlev)  = real(-SHTROOT * zvt)
         pgdlam(:,jlev) = real( SHTROOT * zvp)
      enddo
!$omp end do
      return
      end subroutine sh_sp2grad


      subroutine shtns_teardown
      if (.not. lshtns) return
      call shtns_destroy(shtcfg)
      shtcfg = c_null_ptr
      if (allocated(shtinv)) deallocate(shtinv)
      lshtns = .false.
      return
      end subroutine shtns_teardown

      end module shtnsmod
