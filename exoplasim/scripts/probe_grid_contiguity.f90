      module gridprobe
!     Can a thread's band of a shared full-globe grid array be handed to the
!     physics without gfortran copying it?
!
!     This is the question that decides whether SHTns is a bounded change or an
!     eight-module rewrite. SHTns wants the grid full-globe and contiguous; the
!     physics wants each thread to address only its own latitudes. If a pointer
!     to a band can be passed to an explicit-shape dummy WITHOUT a copy, the 169
!     whole-array `where` statements never need touching, because the pointer
!     already addresses only this thread's band. If it copies, the copy is both
!     the traffic Stage A existed to remove and -- per Stage A's finding -- a
!     correctness hazard, because the copy-out lands after the callee's barrier.
!
!     Decided by ADDRESS, not by behaviour: if the callee sees a different
!     address than the caller's slice, a copy was made. Behaviour cannot tell,
!     because copy-out makes a copied argument look like a passed-through one.
      implicit none
      integer, parameter :: NLON = 512
      integer, parameter :: NLAT = 256
      integer, parameter :: NLEV = 10
      integer, parameter :: NTHR = 16
      integer, parameter :: NLPP = NLAT / NTHR
      integer, parameter :: NHORB = NLON * NLPP     ! a thread's band
      integer, parameter :: NHORG = NLON * NLAT     ! the whole globe

      real, target :: gfull(NHORG,NLEV) = 0.0       ! shared, laid out for SHTns
      contains

      subroutine takes_2d(a, addr)
      real :: a(NHORB,NLEV)                         ! explicit shape, two dims
      integer (kind=8) :: addr
      addr = loc(a)
      a(1,1) = a(1,1) + 1.0
      end subroutine takes_2d

      subroutine takes_1d(a, addr)
      real :: a(NHORB)                              ! explicit shape, one level
      integer (kind=8) :: addr
      addr = loc(a)
      a(1) = a(1) + 1.0
      end subroutine takes_1d

      subroutine takes_assumed(a, addr)
      real :: a(:,:)                                ! assumed shape
      integer (kind=8) :: addr
      addr = loc(a)
      a(1,1) = a(1,1) + 1.0
      end subroutine takes_assumed

      end module gridprobe


      program contig
      use gridprobe
      implicit none

      real, pointer :: gband(:,:)
      real, pointer, contiguous :: gcont(:)
      integer :: it, i0, i1
      integer (kind=8) :: a_slice, a_seen

      it = 3                                        ! any thread but the first
      i0 = it * NHORB + 1
      i1 = i0 + NHORB - 1

      write(*,'(a,i0,a,i0,a,i0)') 'NHORB=', NHORB, '  NHORG=', NHORG, '  NLEV=', NLEV
      write(*,*)

!     1. a 2-D band of the full array, to a two-dimensional dummy. The stride
!        between levels is NHORG and not NHORB, so this section is NOT
!        contiguous and the standard requires the compiler to make it so.
      gband => gfull(i0:i1,:)
      a_slice = loc(gfull(i0,1))
      call takes_2d(gband, a_seen)
      call verdict('2-D band -> explicit-shape (NHORB,NLEV)', a_slice, a_seen)

!     2. the same band ONE LEVEL at a time. Within a level the band is a
!        contiguous run, so there is nothing to copy.
      gcont => gfull(i0:i1,1)
      a_slice = loc(gfull(i0,1))
      call takes_1d(gcont, a_seen)
      call verdict('1-D band, one level  -> explicit-shape (NHORB)', a_slice, a_seen)

!     3. the 2-D band to an ASSUMED-SHAPE dummy, which can carry a stride and
!        so has nothing to copy -- at the cost of a descriptor in the callee.
      a_slice = loc(gfull(i0,1))
      call takes_assumed(gfull(i0:i1,:), a_seen)
      call verdict('2-D band -> assumed-shape (:,:)', a_slice, a_seen)

!     4. the base address itself, the Fortran 77 idiom the model already uses
!        everywhere: pass the first element and let the callee shape it.
      a_slice = loc(gfull(i0,1))
      call takes_1d(gfull(i0,1), a_seen)
      call verdict('base address gfull(i0,1) -> (NHORB)', a_slice, a_seen)

      contains

      subroutine verdict(what, want, got)
      character(len=*) :: what
      integer (kind=8) :: want, got
      if (want == got) then
         write(*,'(a,a,a)') '  [ no copy ] ', what, ''
      else
         write(*,'(a,a,a,i0,a)') '  [  COPY  ] ', what, '   (moved by ', got-want, ' bytes)'
      endif
      end subroutine verdict

      end program contig
