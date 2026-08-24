!     ====================================================================
!     cpl_stub.f90 -- the LSG coupling entry points for a build that has
!     no LSG in it.
!
!     CMakeLists.txt compiles this file and NOT cpl.f90, so plasim.x links
!     these three and nothing else answers to the names. They are non-module
!     externals with no explicit interface, which is why the signatures below
!     have to be maintained by hand against cpl.f90 and against the call
!     sites in oceanmod.f90: nothing diagnoses a mismatch and the link
!     succeeds either way. The version this replaced declared clsgini with
!     six arguments where oceanmod passes eight, and took the integers nlsg
!     and ngui into positions that were declared a real array, so NLSG=1 on
!     this build wrote through them.
!
!     They do not return. A build linked against this file cannot couple to
!     an ocean, so reaching one of them means NLSG > 0 was set on an
!     executable that has no LSG compiled in, and continuing would integrate
!     an uncoupled ocean under a namelist that says otherwise. The argument
!     lists are still declared to match cpl.f90 exactly, so that the failure
!     is the abort rather than the call.
!     ====================================================================
!
!     ==================
!     subroutine CLSGINI
!     ==================
!
      subroutine clsgini(kdatim,ktspd,kaomod,kcpl,kgui,pslm,kxa,kya)
!
      integer :: kdatim(7)  ! date and time
      integer :: ktspd      ! atm time steps per day
      integer :: kaomod     ! atm/ocean step ratio
      integer :: kcpl       ! coupling flag (nlsg)
      integer :: kgui       ! gui switch
      integer :: kxa,kya    ! atm dimensions
      real :: pslm(kxa,kya) ! atm land sea mask
!
      call mpabort('clsgini: NLSG > 0 needs cpl.f90, and this executable '   &
     &           //'was built with cpl_stub.f90')
!
      return
      end subroutine clsgini
!
!     ===================
!     subroutine CLSGSTEP
!     ===================
!
      subroutine clsgstep(kdatim,kstep,psst,ptaux,ptauy,ppme,proff,pice &
     &                   ,pheat,pfldo)
!
      real :: psst(*)      ! atm sst (input)
      real :: ptaux(*)     ! atm u-stress (input)
      real :: ptauy(*)     ! atm v-stress (input)
      real :: ppme(*)      ! atm p-e (input)
      real :: proff(*)     ! atm runoff (input)
      real :: pice(*)      ! atm ice thickness (incl. snow) (input)
      real :: pheat(*)     ! atm heat flux (input; not used yet)
      real :: pfldo(*)     ! atm deep ocean heat flux (output)
      integer :: kdatim(7) ! date and time
      integer :: kstep     ! current atm time step
!
      call mpabort('clsgstep: NLSG > 0 needs cpl.f90, and this executable '  &
     &           //'was built with cpl_stub.f90')
!
      return
      end subroutine clsgstep
!
!     ===================
!     subroutine CLSGSTOP
!     ===================
!
      subroutine clsgstop
!
      call mpabort('clsgstop: NLSG > 0 needs cpl.f90, and this executable '  &
     &           //'was built with cpl_stub.f90')
!
      return
      end subroutine clsgstop
