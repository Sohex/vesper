!     ==================================================================
!     What happens to a caller's array when the restart record it asked
!     for is not in the file? world-onw8.
!
!     Worldbuilding frame: a correctness check on the restart reader of
!     the Vesper climate model. Nothing here is about the simulated
!     planet.
!
!     `get_restart_array` returns without touching its argument when the
!     name is absent and `nexcheck` is 0, which reads like a caller can
!     pre-set a sentinel and test it afterwards. `mpgetgp` breaks that:
!     it reads into a buffer of its own and scatters that buffer into the
!     caller's array whichever way the lookup went, so the sentinel is
!     destroyed by the call that was supposed to leave it. `mpgetgp_found`
!     asks whether the record is there, skips the scatter when it is not,
!     and reports which happened.
!
!     WHAT IS REAL HERE AND WHAT IS NOT. `restartmod.f90` is compiled from
!     the model source, so `restart_ini`, `put_restart_array`,
!     `get_restart_array`, `reseek` and `has_restart_array` are the model's
!     own. The two reader shapes below are transcriptions of `mpgetgp` and
!     `mpgetgp_found` with `mpscgp` reduced to the copy it performs on a
!     one-thread team: the real ones live in `mpimod_omp.f90`, which needs
!     the whole of `pumamod` and a grid. Keep them in step with that file.
!     What the transcription cannot cover -- that every thread takes the
!     same branch -- is structural: the presence flag reaches the team
!     through `mpbci`, the same broadcast `mpsurfgp` uses.
!
!     THE CRITERIA, FIXED BEFORE THE RUN:
!       1. An absent name leaves `get_restart_array`'s argument alone.
!       2. `has_restart_array` answers 1 for a name in the file and 0 for
!          one that is not, and reads nothing while doing it.
!       3. THE CONTROL. The pre-fix `mpgetgp` shape must deliver a NaN
!          into the caller's array for an absent record, and the sentinel
!          test that shape invited -- `ALL(p < 0.0)` -- must come back
!          false. That is the defect, and a check nobody has seen fail is
!          not a check. It arms because the buffer is left unwritten and
!          `-finit-real=snan` fills it; a run at the shipped flag line
!          scatters ordinary stack instead, which is the same defect
!          without an instrument on it.
!       4. THE FIX. The `mpgetgp_found` shape must leave the caller's
!          array at the value set before the call and report 0 for an
!          absent record, and must deliver the record and report 1 for a
!          present one.
!
!     Built and run by verify_absent_restart_record.sh.
!     ==================================================================
      program verify_absent_restart_record
      use restartmod
      implicit none

      integer, parameter :: nugp = 8      ! stands in for the grid size
      real    :: zpresent(nugp,1)
      integer :: nfail

      nfail = 0
      zpresent(:,1) = 1.5

      call build_file(zpresent,nugp)
      call check_premise(nfail)
      call check_presence(nfail)
      call check_control(nfail,nugp)
      call check_fix(nfail,nugp,zpresent)
      call restart_stop

      write(*,*)
      if (nfail == 0) then
         write(*,'(a)') 'PASS: an absent record is reported, not scattered.'
      else
         write(*,'(a,i0,a)') 'FAIL: ',nfail,' check(s) failed.'
         stop 1
      endif
      end program verify_absent_restart_record


!     A restart file carrying one record, 'alpha'. 'beta' is the name that
!     is not in it.
      subroutine build_file(pa,kn)
      use restartmod
      implicit none
      integer :: kn
      real    :: pa(kn,1)
      logical :: lrestart

      call restart_prepare('verify_absent_restart_record.dat')
      call put_restart_array('alpha',pa,kn,kn,1)
      close(nwriunit)
      call restart_ini(lrestart,'verify_absent_restart_record.dat')
      if (.not. lrestart) then
         write(*,'(a)') '  [ FAIL ] the test restart was not written'
         stop 1
      endif
      return
      end subroutine build_file


!     1. An absent name leaves the argument alone. This is the premise the
!     sentinel idiom rests on, and it is true -- for a scalar or an array
!     the caller owns. It is only through mpgetgp that it stops being true.
      subroutine check_premise(kfail)
      use restartmod
      implicit none
      integer :: kfail
      integer, parameter :: kn = 8
      real :: zbuf(kn,1)

      zbuf(:,1) = 42.0
      nexcheck = 0
      call get_restart_array('beta',zbuf,kn,kn,1)
      nexcheck = 1
      if (all(zbuf(:,1) == 42.0)) then
         write(*,'(a)') '  [  ok  ] absent name: get_restart_array left the argument alone'
      else
         write(*,'(a)') '  [ FAIL ] absent name: get_restart_array wrote to the argument'
         kfail = kfail + 1
      endif
      return
      end subroutine check_premise


!     2. has_restart_array answers, and reads nothing while it answers.
      subroutine check_presence(kfail)
      use restartmod
      implicit none
      integer :: kfail
      integer :: ifound, ilast

      ilast = nlastrec
      call has_restart_array('beta',ifound)
      if (ifound == 0) then
         write(*,'(a)') '  [  ok  ] has_restart_array: absent name reports 0'
      else
         write(*,'(a,i0)') '  [ FAIL ] has_restart_array: absent name reported ',ifound
         kfail = kfail + 1
      endif
      call has_restart_array('alpha',ifound)
      if (ifound == 1) then
         write(*,'(a)') '  [  ok  ] has_restart_array: present name reports 1'
      else
         write(*,'(a,i0)') '  [ FAIL ] has_restart_array: present name reported ',ifound
         kfail = kfail + 1
      endif
      if (nlastrec == ilast) then
         write(*,'(a)') '  [  ok  ] has_restart_array read nothing'
      else
         write(*,'(a)') '  [ FAIL ] has_restart_array moved the read position'
         kfail = kfail + 1
      endif
      return
      end subroutine check_presence


!     THE PRE-FIX mpgetgp, transcribed: an unwritten buffer, scattered
!     whatever the lookup did. External and in its own subroutine so the
!     buffer is a real frame the poison flag can fill.
      subroutine getgp_prefix(yn,p,kdim,klev)
      use restartmod
      implicit none
      character (len=*) :: yn
      integer :: kdim, klev
      real :: p(kdim,klev)
      integer, parameter :: nugp = 8
      real :: z(nugp,klev)
      call get_restart_array(yn,z,nugp,nugp,klev)
      p(1:nugp,:) = z(1:nugp,:)     ! mpscgp on a one-thread team
      return
      end subroutine getgp_prefix


!     mpgetgp_found, transcribed. mypid == NROOT throughout here, so the
!     mpbci broadcast of iread has nothing to carry and is left out; what
!     it guarantees in the model is that the team agrees on the branch.
      subroutine getgp_found(yn,p,kdim,klev,kfound)
      use restartmod
      implicit none
      character (len=*) :: yn
      integer :: kdim, klev, kfound
      real :: p(kdim,klev)
      integer, parameter :: nugp = 8
      real :: z(nugp,klev)
      integer :: iread(1)
      iread(1) = 0
      z(:,:) = 0.0
      call has_restart_array(yn,iread(1))
      if (iread(1) == 1) call get_restart_array(yn,z,nugp,nugp,klev)
      if (iread(1) == 1) p(1:nugp,:) = z(1:nugp,:)
      kfound = iread(1)
      return
      end subroutine getgp_found


!     3. THE CONTROL. The pre-fix shape must destroy the sentinel.
      subroutine check_control(kfail,kn)
      use restartmod
      use ieee_arithmetic
      implicit none
      integer :: kfail, kn
      real :: zp(kn,1)

      zp(:,1) = -1.0                   ! landmod's dwatcl sentinel
      nexcheck = 0
      call getgp_prefix('beta',zp,kn,1)
      nexcheck = 1
      if (any(ieee_is_nan(zp(:,1)))) then
         write(*,'(a)') '  [  ok  ] control: the pre-fix shape scattered its unwritten buffer'
      else
         write(*,'(a)') '  [ FAIL ] control: the pre-fix shape did not scatter, so the'
         write(*,'(a)') '           poison never reached the buffer and this check is blind'
         kfail = kfail + 1
      endif
      if (.not. all(zp(:,1) < 0.0)) then
         write(*,'(a)') '  [  ok  ] control: the sentinel test came back false on an absent record'
      else
         write(*,'(a)') '  [ FAIL ] control: the sentinel survived, so the defect is not armed here'
         kfail = kfail + 1
      endif
      return
      end subroutine check_control


!     4. THE FIX.
      subroutine check_fix(kfail,kn,pexpect)
      use restartmod
      implicit none
      integer :: kfail, kn
      real :: pexpect(kn,1)
      real :: zp(kn,1)
      integer :: ifound

      zp(:,1) = -1.0
      call getgp_found('beta',zp,kn,1,ifound)
      if (ifound == 0 .and. all(zp(:,1) == -1.0)) then
         write(*,'(a)') '  [  ok  ] fix: absent record reports 0 and leaves the array alone'
      else
         write(*,'(a,i0)') '  [ FAIL ] fix: absent record reported ',ifound
         kfail = kfail + 1
      endif

      zp(:,1) = -1.0
      call getgp_found('alpha',zp,kn,1,ifound)
      if (ifound == 1 .and. all(zp(:,1) == pexpect(:,1))) then
         write(*,'(a)') '  [  ok  ] fix: present record reports 1 and delivers the values'
      else
         write(*,'(a,i0)') '  [ FAIL ] fix: present record reported ',ifound
         kfail = kfail + 1
      endif
      return
      end subroutine check_fix
