!     ==================================================================
!     mpimod_omp.f90
!     --------------
!     A third implementation of the MPI layer's contract, in which the
!     ranks are THREADS of one process. It replaces <mpimod.f90> the way
!     <mpimod_stub.f90> does, selected by ${MPIMOD} in make_plasim.
!
!     WHY. SHTns has no distributed mode -- its parallelism is OpenMP --
!     so a spectral transform library needs every latitude reachable in
!     one address space. Independently of that, in one address space
!     every gather and scatter is a memcpy rather than a wire transfer.
!
!     HOW IT IS THE SAME MODEL. Every module holding per-rank state is
!     declared !$omp threadprivate, so a thread owns an NHOR chunk
!     exactly as a rank owns one, and the physics is untouched: its
!     whole-array expressions address the thread's chunk as written.
!     `mypid` is the thread number, and the model's 345 tests of
!     `mypid == NROOT` are true on thread 0 and nowhere else, as before.
!
!     THE RULES THIS FILE OBEYS, and they are not style:
!
!     1. Every routine here that reaches a barrier must be reached by
!        EVERY thread. That is the same constraint an MPI collective
!        already imposes, and the model already satisfies it.
!     2. A staging buffer is written, read, and only then reused. That
!        is two barriers, not one: without the second, a thread racing
!        ahead overwrites a buffer another thread has not yet read.
!     3. Reductions are summed in a FIXED thread order, never in the
!        order threads arrive. Arrival order makes the model
!        non-deterministic run to run, which is archive CLIM-44 and has
!        been paid for here once already.
!     ==================================================================

      module mpiomp
      use pumamod, only: NUGP, NESP, NSPP, NLEV, NLAT, NLON, NLPP,      &
     &                   NPRO, NROOT
      implicit none

!     SHARED across the team. Nothing in this module may be threadprivate:
!     it IS the shared memory that stands in for the wire.

      integer, parameter :: NGEN = max(NUGP, NESP*NLEV)

      real         :: zbufgp(NUGP)                 ! one grid level
      real         :: zbufcs(NLAT)                 ! one cross section level
      real         :: zbufgen(NGEN)                ! broadcasts, small reductions
      real (kind=8):: zbufd(NGEN)                  ! the kind=8 scatters
      integer      :: kbufgen(NGEN)
      logical      :: lbufgen

!     Per-thread partial sums, for the reductions that must stay ordered.
      real         :: zbufred(NESP,NLEV,0:NPRO-1)

      end module mpiomp


!     ==================================================================
!     Broadcasts. Thread 0 publishes, everybody reads.
!     ==================================================================

      subroutine mpbci(k) ! broadcast 1 integer
      use pumamod
      use mpiomp
      integer :: k(*)
      if (mypid == NROOT) kbufgen(1) = k(1)
!$omp barrier
      k(1) = kbufgen(1)
!$omp barrier
      return
      end subroutine mpbci


      subroutine mpbcin(k,n) ! broadcast n integer
      use pumamod
      use mpiomp
      integer :: k(n)
      integer :: n
      if (n > NGEN) call mpabort('mpbcin: buffer too small')
      if (mypid == NROOT) kbufgen(1:n) = k(1:n)
!$omp barrier
      k(1:n) = kbufgen(1:n)
!$omp barrier
      return
      end subroutine mpbcin


      subroutine mpbcr(p) ! broadcast 1 real
      use pumamod
      use mpiomp
      real :: p(*)
      if (mypid == NROOT) zbufgen(1) = p(1)
!$omp barrier
      p(1) = zbufgen(1)
!$omp barrier
      return
      end subroutine mpbcr


      subroutine mpbcrn(p,n) ! broadcast n real
      use pumamod
      use mpiomp
      integer :: n
      real :: p(n)
      if (n > NGEN) call mpabort('mpbcrn: buffer too small')
      if (mypid == NROOT) zbufgen(1:n) = p(1:n)
!$omp barrier
      p(1:n) = zbufgen(1:n)
!$omp barrier
      return
      end subroutine mpbcrn


      subroutine mpbcl(l) ! broadcast 1 logical
      use pumamod
      use mpiomp
      logical :: l(*)
      if (mypid == NROOT) lbufgen = l(1)
!$omp barrier
      l(1) = lbufgen
!$omp barrier
      return
      end subroutine mpbcl


!     ==================================================================
!     Scatters of a latitude vector. Thread 0 holds the global array and
!     keeps it: mpi_scatter is in-place on the root's first chunk, so a
!     rank 0 that reads the vector globally afterwards still can -- and
!     oceanmod's hdiffo does exactly that. Preserve the behaviour.
!     ==================================================================

      subroutine mpscin(k,n) ! scatter n integer
      use pumamod
      use mpiomp
      integer :: n
      integer :: k(*)
      if (n*NPRO > NGEN) call mpabort('mpscin: buffer too small')
      if (mypid == NROOT) kbufgen(1:n*NPRO) = k(1:n*NPRO)
!$omp barrier
      if (mypid /= NROOT) k(1:n) = kbufgen(mypid*n+1:mypid*n+n)
!$omp barrier
      return
      end subroutine mpscin


      subroutine mpscrn(p,n) ! scatter n real
      use pumamod
      use mpiomp
      integer :: n
      real :: p(*)
      if (n*NPRO > NGEN) call mpabort('mpscrn: buffer too small')
      if (mypid == NROOT) zbufgen(1:n*NPRO) = p(1:n*NPRO)
!$omp barrier
      if (mypid /= NROOT) p(1:n) = zbufgen(mypid*n+1:mypid*n+n)
!$omp barrier
      return
      end subroutine mpscrn


      subroutine mpscdn(p,n) ! scatter n double precision
      use pumamod
      use mpiomp
      integer :: n
      real (kind=8) :: p(*)
      if (n*NPRO > NGEN) call mpabort('mpscdn: buffer too small')
      if (mypid == NROOT) zbufd(1:n*NPRO) = p(1:n*NPRO)
!$omp barrier
      if (mypid /= NROOT) p(1:n) = zbufd(mypid*n+1:mypid*n+n)
!$omp barrier
      return
      end subroutine mpscdn


!     ==================================================================
!     Grid point transfers. These carry the paired latitude permutation,
!     and here it is simpler than in the MPI version: a thread reads the
!     latitudes it owns straight out of the shared array at the global
!     offsets ilatperm gives, so the same code serves both layouts --
!     without LPAIRLAT the map is the identity and the slots are the
!     contiguous block the stock decomposition hands out.
!     ==================================================================

      subroutine mpscgp(pf,pp,klev) ! scatter gridpoint fields
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NUGP,klev)
      real :: pp(NHOR,klev)
      integer :: jlev, jlat, jg

      do jlev = 1 , klev
         if (mypid == NROOT) zbufgp(:) = pf(:,jlev)
!$omp barrier
         do jlat = 1 , NLPP
            jg = ilatperm(mypid*NLPP + jlat)
            pp(1+(jlat-1)*NLON:jlat*NLON,jlev) =                        &
     &         zbufgp(1+(jg-1)*NLON:jg*NLON)
         enddo
!$omp barrier
      enddo
      return
      end subroutine mpscgp


      subroutine mpgagp(pf,pp,klev) ! gather gridpoint fields
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NUGP,klev)
      real :: pp(NHOR,klev)
      integer :: jlev, jlat, jg

      do jlev = 1 , klev
         do jlat = 1 , NLPP
            jg = ilatperm(mypid*NLPP + jlat)
            zbufgp(1+(jg-1)*NLON:jg*NLON) =                             &
     &         pp(1+(jlat-1)*NLON:jlat*NLON,jlev)
         enddo
!$omp barrier
         if (mypid == NROOT) pf(:,jlev) = zbufgp(:)
!$omp barrier
      enddo
      return
      end subroutine mpgagp


      subroutine mpgallgp(pf,pp,klev) ! gather gridpoint to all
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NUGP,klev)
      real :: pp(NHOR,klev)
      integer :: jlev, jlat, jg

      do jlev = 1 , klev
         do jlat = 1 , NLPP
            jg = ilatperm(mypid*NLPP + jlat)
            zbufgp(1+(jg-1)*NLON:jg*NLON) =                             &
     &         pp(1+(jlat-1)*NLON:jlat*NLON,jlev)
         enddo
!$omp barrier
         pf(:,jlev) = zbufgp(:)
!$omp barrier
      enddo
      return
      end subroutine mpgallgp


      subroutine mpgacs(pcs) ! gather cross sections
      use pumamod
      use mpiomp
      real :: pcs(NLAT,NLEV)
      integer :: jlev, jlat

      do jlev = 1 , NLEV
         do jlat = 1 , NLPP
            zbufcs(ilatperm(mypid*NLPP + jlat)) = pcs(jlat,jlev)
         enddo
!$omp barrier
         if (mypid == NROOT) pcs(:,jlev) = zbufcs(:)
!$omp barrier
      enddo
      return
      end subroutine mpgacs


!     ==================================================================
!     Spectral transfers, on the NSPP mode decomposition.
!     ==================================================================

      subroutine mpscsp(pf,pp,klev) ! scatter spectral fields
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)
      integer :: jlev

      do jlev = 1 , klev
         if (mypid == NROOT) zbufgen(1:NESP) = pf(:,jlev)
!$omp barrier
         pp(:,jlev) = zbufgen(mypid*NSPP+1:mypid*NSPP+NSPP)
!$omp barrier
      enddo
      return
      end subroutine mpscsp


      subroutine mpgasp(pf,pp,klev) ! gather spectral fields
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)
      integer :: jlev

      do jlev = 1 , klev
         zbufgen(mypid*NSPP+1:mypid*NSPP+NSPP) = pp(:,jlev)
!$omp barrier
         if (mypid == NROOT) pf(:,jlev) = zbufgen(1:NESP)
!$omp barrier
      enddo
      return
      end subroutine mpgasp


      subroutine mpgallsp(pf,pp,klev) ! gather spectral to all
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)
      integer :: jlev

      do jlev = 1 , klev
         zbufgen(mypid*NSPP+1:mypid*NSPP+NSPP) = pp(:,jlev)
!$omp barrier
         pf(:,jlev) = zbufgen(1:NESP)
!$omp barrier
      enddo
      return
      end subroutine mpgallsp


!     ==================================================================
!     Reductions. Every one of these sums over threads 0,1,...,NPRO-1 in
!     that order, on every thread, so the answer does not depend on which
!     thread arrived first. Determinism is the requirement; speed is why
!     the partials go to a per-thread buffer rather than through a
!     critical section.
!     ==================================================================

      subroutine mpsumsc(psf,psp,klev) ! sum & scatter spectral
      use pumamod
      use mpiomp
      integer :: klev
      real :: psf(NESP,klev)
      real :: psp(NSPP,klev)
      integer :: jlev, j, w, it
      real :: z

      zbufred(:,1:klev,mypid) = psf(:,1:klev)
!$omp barrier
      do jlev = 1 , klev
         do j = 1 , NSPP
            w = mypid*NSPP + j
            z = zbufred(w,jlev,0)
            do it = 1 , NPRO-1
               z = z + zbufred(w,jlev,it)
            enddo
            psp(j,jlev) = z
         enddo
      enddo
!$omp barrier
      return
      end subroutine mpsumsc


      subroutine mpgallspp(pf,pp,klev) ! gather to all, into shared storage
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)

!     mpgallsp with the staging removed. pf is one array the whole team can
!     see, so a thread writing its own slice of it IS the gather, and the
!     buffer every partial used to be copied through is not needed.
!
!     Two barriers, and both are load bearing. The first keeps this write clear
!     of whatever the previous use of pf was still reading -- these arrays are
!     scratch and get reused. The second stops a reader running before every
!     slice is written.
!$omp barrier
      pf(mypid*NSPP+1:mypid*NSPP+NSPP,1:klev) = pp(:,1:klev)
!$omp barrier
      return
      end subroutine mpgallspp


      subroutine mpzerosp(pf,kfrom,klev) ! clear a replicated spectral array
      use pumamod
      use mpiomp
      integer :: kfrom, klev
      real :: pf(NESP,klev)
      integer :: lo, hi

!     pf is shared, so every thread clearing the whole of it would be a race --
!     on identical values, and still a race, and still reported. Each clears
!     only the slice it owns. Modes below kfrom survive wherever they live.
      lo = max(kfrom, mypid*NSPP + 1)
      hi = mypid*NSPP + NSPP
!$omp barrier
      if (lo <= hi) pf(lo:hi,1:klev) = 0.0
!$omp barrier
      return
      end subroutine mpzerosp


      subroutine mpsumscp(ppart,psp,klev) ! sum & scatter, partials in place
      use pumamod
      use mpiomp
      integer :: klev
      real :: ppart(NESP,klev,0:NPRO-1)
      real :: psp(NSPP,klev)
      integer :: jlev, j, w, it
      real :: z

!     The same reduction as mpsumsc and without its staging copy: the callers
!     wrote their partials into ppart already, so there is nothing to move
!     before summing. That copy was 68 GB over a 300-step T127 run and it is
!     the whole reason this routine exists beside the other one.
!
!     The barrier is still required. A thread reads every other thread's slot,
!     so it may not start until all of them are written.
!$omp barrier
      do jlev = 1 , klev
         do j = 1 , NSPP
            w = mypid*NSPP + j
            z = ppart(w,jlev,0)
            do it = 1 , NPRO-1
               z = z + ppart(w,jlev,it)
            enddo
            psp(j,jlev) = z
         enddo
      enddo
!$omp barrier
      return
      end subroutine mpsumscp


      subroutine mpsum(psp,klev) ! sum spectral fields onto the root
      use pumamod
      use mpiomp
      integer :: klev
      real :: psp(NESP,klev)
      integer :: jlev, j, it
      real :: z

      zbufred(:,1:klev,mypid) = psp(:,1:klev)
!$omp barrier
!     Split the sum across threads so it is not serial, then publish it.
      do jlev = 1 , klev
         do j = mypid*NSPP+1 , mypid*NSPP+NSPP
            z = zbufred(j,jlev,0)
            do it = 1 , NPRO-1
               z = z + zbufred(j,jlev,it)
            enddo
            zbufred(j,jlev,NROOT) = z
         enddo
      enddo
!$omp barrier
      if (mypid == NROOT) psp(:,1:klev) = zbufred(:,1:klev,NROOT)
!$omp barrier
      return
      end subroutine mpsum


      subroutine mpsumr(pr,kdim) ! sum kdim reals onto the root
      use pumamod
      use mpiomp
      integer :: kdim
      real :: pr(kdim)
      integer :: j, it

      if (kdim*NPRO > NGEN) call mpabort('mpsumr: buffer too small')
      zbufgen(mypid*kdim+1:mypid*kdim+kdim) = pr(1:kdim)
!$omp barrier
      if (mypid == NROOT) then
         do j = 1 , kdim
            pr(j) = zbufgen(j)
            do it = 1 , NPRO-1
               pr(j) = pr(j) + zbufgen(it*kdim+j)
            enddo
         enddo
      endif
!$omp barrier
      return
      end subroutine mpsumr


      subroutine mpsumbcr(pr,kdim) ! sum & broadcast kdim reals
      use pumamod
      use mpiomp
      integer :: kdim
      real :: pr(kdim)
      integer :: j, it

      if (kdim*NPRO > NGEN) call mpabort('mpsumbcr: buffer too small')
      zbufgen(mypid*kdim+1:mypid*kdim+kdim) = pr(1:kdim)
!$omp barrier
      do j = 1 , kdim
         pr(j) = zbufgen(j)
         do it = 1 , NPRO-1
            pr(j) = pr(j) + zbufgen(it*kdim+j)
         enddo
      enddo
!$omp barrier
      return
      end subroutine mpsumbcr


      subroutine mpmaxval(p,kdim,klev,pmax)
      use pumamod
      use mpiomp
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: pmax
      integer :: it

      zbufgen(mypid+1) = maxval(p(:,:))
!$omp barrier
      pmax = zbufgen(1)
      do it = 1 , NPRO-1
         pmax = max(pmax, zbufgen(it+1))
      enddo
!$omp barrier
      return
      end subroutine mpmaxval


      subroutine mpsumval(p,kdim,klev,psum)
      use pumamod
      use mpiomp
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: psum
      integer :: it

      zbufgen(mypid+1) = sum(p(:,:))
!$omp barrier
      psum = zbufgen(1)
      do it = 1 , NPRO-1
         psum = psum + zbufgen(it+1)
      enddo
!$omp barrier
      return
      end subroutine mpsumval


!     ==================================================================
!     Start, stop, identity.
!     ==================================================================

      subroutine mpstart(kworld) ! initialization
      use pumamod
      use mpiomp
!     The two runtime queries are declared EXTERNAL rather than taken from
!     omp_lib, because flang does not ship omp_lib.mod -- it carries seventeen
!     modules and not that one -- and this build has to compile under both
!     compilers so that an OpenMP-aware race detector can be pointed at it.
!     Both libgomp and LLVM's libomp export these under the Fortran name.
      integer :: omp_get_thread_num
      integer :: omp_get_num_threads
      external :: omp_get_thread_num, omp_get_num_threads
      integer :: kworld
      integer :: iteam

!     THE GUARD. nproc's declaration initialiser is NPRO, so reading NPRO
!     here means this thread received the module's initialisers through the
!     TLS initialisation image. That is what makes a copyin clause over
!     1194 threadprivate variables unnecessary -- and it is gfortran's
!     behaviour rather than anything OpenMP promises, so it is checked on
!     every thread instead of assumed. Without it, a compiler that zeroed
!     the copies instead would give every thread but the master a model
!     built from zeroed constants: a wrong climate, not a crash.
      if (nproc /= NPRO) then
         call mpabort('threadprivate initialisers did not reach this thread')
      endif

      iteam   = 1
      mypid   = 0
      myworld = 0
      mpinfo  = 0
!$    iteam = omp_get_num_threads()
!$    mypid = omp_get_thread_num()
      nproc = iteam
!     The threads share an address space, so each needs its own slot in the
!     tendency partials. Under MPI mypart stays 0 and there is one slot.
      mypart = mypid

      if (nproc /= NPRO) then
         if (mypid == 0) then
            write(nud,*)'Compiled for ',NPRO,' threads'
            write(nud,*)'Running on  ',nproc,' threads'
         endif
         call mpabort('thread count does not match NPRO')
      endif

      if (mypid == NROOT) then
         allocate(ympname(nproc)) ; ympname(:) = 'thread'
      endif
!$omp barrier

      call assoc_spectral
      return
      end subroutine mpstart


      subroutine mpstop
      return
      end subroutine mpstop


      subroutine ompi_info(nprocess,npid) ! get nproc and pid
      use pumamod
      integer :: nprocess, npid
      nprocess = nproc
      npid     = mypid
      return
      end subroutine ompi_info


!     ==================================================================
!     Serial file paths. Every one of these already runs on the root
!     alone in the MPI build, so they are the root thread's here too, and
!     the gather or scatter around them is the one above.
!     ==================================================================

      subroutine mpreadgp(ktape,p,kdim,klev)
      use pumamod
      integer :: ktape, kdim, klev
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      z = 0.0
      if (mypid == NROOT) read (ktape) z(:,:)
      if (kdim == NHOR) then
         call mpscgp(z,p,klev)
      else
         if (mypid == NROOT) p = z
      endif
      return
      end subroutine mpreadgp


      subroutine mpwritegp(ktape,p,kdim,klev)
      use pumamod
      integer :: ktape, kdim, klev
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      if (kdim == NHOR) then
         call mpgagp(z,p,klev)
         if (mypid == NROOT) write(ktape) z(1:NUGP,:)
      else
         if (mypid == NROOT) write(ktape) p(1:NUGP,:)
      endif
      return
      end subroutine mpwritegp


      subroutine mpwritegph(ktape,p,kdim,klev,ihead)
      use pumamod
      integer :: ktape, kdim, klev
      integer :: ihead(8)
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      real (kind=4) :: zz(NUGP,klev)
      real (kind=4) :: zp(kdim,klev)
      if (kdim == NHOR) then
         call mpgagp(z,p,klev)
         if (mypid == NROOT) then
            write(ktape) ihead
            zz(:,:) = z(1:NUGP,:)
            write(ktape) zz(1:NUGP,:)
         endif
      else
         if (mypid == NROOT) then
            write(ktape) ihead
            zp(:,:) = p(:,:)
            write(ktape) zp(1:NUGP,:)
         endif
      endif
      return
      end subroutine mpwritegph


      subroutine mpreadsp(ktape,p,kdim,klev)
      use pumamod
      integer :: ktape, kdim, klev
      real :: p(kdim,klev)
      real :: z(NESP,klev)
      integer :: i, j
      z = 0.0
      if (mypid == NROOT) read(ktape) ((z(i,j),i=1,NRSP),j=1,klev)
      if (kdim == NSPP) then
         call mpscsp(z,p,klev)
      else
         if (mypid == NROOT) p = z
      endif
      return
      end subroutine mpreadsp


      subroutine mpwritesp(ktape,p,kdim,klev)
      use pumamod
      integer :: ktape, kdim, klev
      real :: p(kdim,klev)
      real :: z(NESP,klev)
      integer :: i, j
      if (kdim == NSPP) then
         call mpgasp(z,p,klev)
         if (mypid == NROOT) write(ktape) ((z(i,j),i=1,NRSP),j=1,klev)
      else
         if (mypid == NROOT) write(ktape) ((p(i,j),i=1,NRSP),j=1,klev)
      endif
      return
      end subroutine mpwritesp


      subroutine mpgetsp(yn,p,kdim,klev)
      use pumamod
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: z(NESP,klev)
      z(:,:) = 0.0
      if (mypid == NROOT) call get_restart_array(yn,z,NRSP,NESP,klev)
      call mpscsp(z,p,klev)
      return
      end subroutine mpgetsp


      subroutine mpputsp(yn,p,kdim,klev)
      use pumamod
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: z(NESP,klev)
      call mpgasp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NRSP,NESP,klev)
      return
      end subroutine mpputsp


      subroutine mpgetgp(yn,p,kdim,klev)
      use pumamod
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      if (mypid == NROOT) call get_restart_array(yn,z,NUGP,NUGP,klev)
      call mpscgp(z,p,klev)
      return
      end subroutine mpgetgp


      subroutine mpputgp(yn,p,kdim,klev)
      use pumamod
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      call mpgagp(z,p,klev)
      if (mypid == NROOT) call put_restart_array(yn,z,NUGP,NUGP,klev)
      return
      end subroutine mpputgp


      subroutine mpsurfgp(yn,p,kdim,klev)
      use pumamod
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      real :: z(NUGP,klev)
      integer :: iread(1)
      iread(1) = 0
      if (mypid == NROOT) call get_surf_array(yn,z,NUGP,klev,iread)
      call mpbci(iread)
      if (iread(1) == 1) call mpscgp(z,p,klev)
      return
      end subroutine mpsurfgp


!     ==================================================================
!     Abort. Written by whichever thread reaches it; the file open is
!     serialised by libgfortran and the process stops either way.
!     ==================================================================

      subroutine mpabort(ym)
      use pumamod
      character (len=* ) :: ym
      character (len=64) :: ystar = ' '
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(ystar)
      character (len=64) :: ymess = ' '
!     Implicitly SAVE, so one copy shared by the whole team.
!$omp threadprivate(ymess)
      integer :: j, ilen

      ilen = 60
      do j = 1 , ilen+4
         ystar(j:j) = '*'
      enddo
      ymess(1:2) = '* '
      ymess(3:min(64,2+len_trim(ym))) = trim(ym)

      open (44,file='Abort_Message')
      write(44,'(A)') trim(ystar)
      write(44,'(A)') trim(ymess)
      write(44,'(A)') trim(ystar)
      close(44)
      write(nud,'(/,A)') trim(ystar)
      write(nud,'(A)')   trim(ymess)
      write(nud,'(A,/)') trim(ystar)

      stop
      end subroutine mpabort


!     ==================================================================
!     Multirun. Not supported by this variant, and it says so rather than
!     doing something plausible: two model instances in one process would
!     share every threadprivate copy by instance as well as by thread.
!     ==================================================================

      subroutine mrsum(k)
      integer :: k
      return
      end subroutine mrsum

      subroutine mrbci(k)
      integer :: k
      return
      end subroutine mrbci

      subroutine mrdiff(p,d,n)
      integer :: n
      real :: p(n)
      real :: d(n)
      return
      end subroutine mrdiff

      subroutine mrdimensions
      return
      end subroutine mrdimensions


!     ==================================================================
!     SUBROUTINE MPGATHERSP
!     ==================================================================

!     Make pf hold the full spectral field, given that each thread has
!     contributed its own slice in pp.
!
!     Here that is a BARRIER and nothing else: pp is a pointer into pf, so
!     the contribution was written in place and there is nothing to move.
!     This is where the 105 GB a run of staging traffic went.

      subroutine mpgathersp(pf,pp,klev)
      use pumamod
      use mpiomp
      integer :: klev
      real :: pf(NESP,klev)
      real :: pp(NSPP,klev)
!$omp barrier
      return
      end subroutine mpgathersp


!     ==================================================================
!     SUBROUTINE MPSYNCSP
!     ==================================================================

!     Publish the spectral state. Under the shared build each thread has
!     already written its own slice of the full array, so this is a barrier
!     and takes no arguments -- deliberately, because passing a partial
!     across a call boundary would make the compiler copy it in and out and
!     the copy-out would land AFTER this barrier.

      subroutine mpsyncsp
!$omp barrier
      return
      end subroutine mpsyncsp
