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


      subroutine shtns_teardown
      if (.not. lshtns) return
      call shtns_destroy(shtcfg)
      shtcfg = c_null_ptr
      if (allocated(shtinv)) deallocate(shtinv)
      lshtns = .false.
      return
      end subroutine shtns_teardown

      end module shtnsmod
