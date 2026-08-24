!     ==================================================================
!     Does mpimod_omp compute what mpimod computes?
!
!     Every collective it covers has an answer that can be written down in
!     advance, so this checks them against those answers rather than
!     against another run of the model. It is built and run by
!     verify_omp_collectives.sh. It is also the only reader of mpsumr, which
!     is why that routine survived world-ro6's sweep of the omp layer.
!
!     Worldbuilding frame: a correctness check on the Vesper climate
!     model's parallel layer. Nothing here is about the simulated planet.
!     ==================================================================
      program verify_omp_collectives
      use pumamod
!$    use omp_lib
      implicit none

      integer :: nfail
      nfail = 0

!$omp parallel num_threads(NPRO) default(shared) reduction(+:nfail)
      call runchecks(nfail)
!$omp end parallel

      write(*,*)
      if (nfail == 0) then
         write(*,*) "PASS: every collective returns what mpimod returns"
      else
         write(*,'(a,i4)') " FAIL: checks failed: ", nfail
         stop 1
      endif
      end program verify_omp_collectives


      subroutine runchecks(nfail)
      use pumamod
      implicit none
      integer :: nfail

      integer :: k(4), j, jl
      real    :: r(4)
      logical :: l(1)
      real    :: zgp(NUGP), zhor(NHOR), zback(NUGP)
      real    :: zsp(NESP), zspp(NSPP)
      real    :: zsum(NESP)
      real    :: zexp, zmax
      integer :: iwant

      call mpstart(-1)

!     ---- broadcast: the root's value reaches every thread -------------
      k(1) = -1
      if (mypid == NROOT) k(1) = 12345
      call mpbci(k)
      call expecti(k(1), 12345, 'mpbci', nfail)

      r(1:4) = -1.0
      if (mypid == NROOT) r(1:4) = (/1.5,2.5,3.5,4.5/)
      call mpbcrn(r,4)
      call expectr(r(4), 4.5, 'mpbcrn', nfail)

      l(1) = .false.
      if (mypid == NROOT) l(1) = .true.
      call mpbcl(l)
      if (.not. l(1)) call failed('mpbcl', nfail)

!     ---- grid round trip: scatter then gather must be the identity ----
!     Every cell carries its own global index, so a permutation that is
!     not its own inverse shows up as a wrong VALUE and not merely a
!     wrong total.
      if (mypid == NROOT) then
         do j = 1 , NUGP
            zgp(j) = real(j)
         enddo
      endif
      call mpscgp(zgp,zhor,1)
      zback(:) = 0.0
      call mpgagp(zback,zhor,1)
      if (mypid == NROOT) then
         do j = 1 , NUGP
            if (zback(j) /= real(j)) then
               call failed('mpscgp/mpgagp round trip', nfail)
               exit
            endif
         enddo
      endif

!     ---- a thread must receive a CONTIGUOUS block of latitudes ---------
!     Thread r holds global latitudes r*NLPP+1 .. (r+1)*NLPP, which is the
!     placement contract every routine that indexes a band by arithmetic
!     rather than by a lookup depends on. The round trip above cannot see
!     this: a scatter and a gather that share one permutation are still
!     mutually inverse, so they compose to the identity while every thread
!     holds latitudes other than the ones it is assumed to hold. Reading the
!     scattered array directly is the only thing that separates the two.
      do jl = 1 , NLPP
         iwant = mypid*NLPP + jl
         if (zhor(1+(jl-1)*NLON) /= real(1+(iwant-1)*NLON)) then
            call failed('mpscgp latitude placement', nfail)
            exit
         endif
      enddo

!     ---- spectral scatter and gather ----------------------------------
      if (mypid == NROOT) then
         do j = 1 , NESP
            zsp(j) = real(j)
         enddo
      endif
      call mpscsp(zsp,zspp,1)
      do j = 1 , NSPP
         if (zspp(j) /= real(mypid*NSPP + j)) then
            call expectr(zspp(j), real(mypid*NSPP + j), 'mpscsp', nfail)
            exit
         endif
      enddo

      zsp(:) = 0.0
      call mpgasp(zsp,zspp,1)
      if (mypid == NROOT) call expectr(zsp(NESP), real(NESP), 'mpgasp', nfail)

      zsp(:) = 0.0
      call mpgallsp(zsp,zspp,1)
      call expectr(zsp(NESP), real(NESP), 'mpgallsp', nfail)

!     ---- sum onto the root --------------------------------------------
!     Thread t contributes (t+1)*j at MODE j, so the sum is
!     j*NPRO(NPRO+1)/2 and it DEPENDS ON THE MODE INDEX. A uniform
!     contribution would let any index confusion cancel, and a check that
!     cannot fail on the mistake being made is not a check.
      zexp = real(NPRO*(NPRO+1)/2)
      do j = 1 , NESP
         zsum(j) = real(mypid + 1) * real(j)
      enddo
      call mpsum(zsum,1)
      if (mypid == NROOT) then
         do j = 1 , NESP
            if (abs(zsum(j) - zexp*real(j)) > 1.0e-6*zexp*real(j)) then
               call expectr(zsum(j), zexp*real(j), 'mpsum', nfail)
               exit
            endif
         enddo
      endif

!     ---- sum and scatter ----------------------------------------------
!     Thread t must come back with the sum for ITS OWN modes, which are
!     t*NSPP+1 .. t*NSPP+NSPP -- the one thing a reduce_scatter has to get
!     right and the one thing uniform data hides.
      do j = 1 , NESP
         zsum(j) = real(mypid + 1) * real(j)
      enddo
      zspp(:) = 0.0
      call mpsumsc(zsum,zspp,1)
      do j = 1 , NSPP
         if (abs(zspp(j) - zexp*real(mypid*NSPP+j))                      &
     &       > 1.0e-6*zexp*real(mypid*NSPP+j)) then
            call expectr(zspp(j), zexp*real(mypid*NSPP+j), 'mpsumsc', nfail)
            exit
         endif
      enddo

!     ---- the small reductions ------------------------------------------
      do j = 1 , 4
         r(j) = real(mypid + 1) * real(j)
      enddo
      call mpsumbcr(r,4)
      call expectr(r(1), zexp,       'mpsumbcr element 1', nfail)
      call expectr(r(4), zexp*4.0,   'mpsumbcr element 4', nfail)

      do j = 1 , 4
         r(j) = real(mypid + 1) * real(j)
      enddo
      call mpsumr(r,4)
      if (mypid == NROOT) then
         call expectr(r(1), zexp,     'mpsumr element 1', nfail)
         call expectr(r(4), zexp*4.0, 'mpsumr element 4', nfail)
      endif

      zhor(:) = real(mypid + 1)
      call mpmaxval(zhor,NHOR,1,zmax)
      call expectr(zmax, real(NPRO), 'mpmaxval', nfail)

      zhor(:) = 1.0
      call mpsumval(zhor,NHOR,1,zmax)
      call expectr(zmax, real(NUGP), 'mpsumval', nfail)

      return
      end subroutine runchecks


      subroutine expecti(igot,iwant,yn,nfail)
      use pumamod
      implicit none
      integer :: igot, iwant, nfail
      character (len=*) :: yn
      if (igot /= iwant) then
!$omp critical
         write(*,'(a,i3,a,a,a,i12,a,i12)') " thread",mypid," FAIL ",yn,      &
     &        " got",igot," want",iwant
!$omp end critical
         nfail = nfail + 1
      endif
      return
      end subroutine expecti


      subroutine expectr(zgot,zwant,yn,nfail)
      use pumamod
      implicit none
      real :: zgot, zwant
      integer :: nfail
      character (len=*) :: yn
      if (abs(zgot - zwant) > 1.0e-9 * max(abs(zwant),1.0)) then
!$omp critical
         write(*,'(a,i3,a,a,a,e16.8,a,e16.8)') " thread",mypid," FAIL ",yn,  &
     &        " got",zgot," want",zwant
!$omp end critical
         nfail = nfail + 1
      endif
      return
      end subroutine expectr


      subroutine failed(yn,nfail)
      use pumamod
      implicit none
      character (len=*) :: yn
      integer :: nfail
!$omp critical
      write(*,'(a,i3,a,a)') " thread",mypid," FAIL ",yn
!$omp end critical
      nfail = nfail + 1
      return
      end subroutine failed
