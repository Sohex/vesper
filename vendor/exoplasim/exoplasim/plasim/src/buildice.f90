

      program buildice
      use resmod

!     THE GRID COMES FROM resmod, like every other unit in this tree. It used to
!     be `NLAT = 32, NLON = 64` declared here, which is T21, and this program
!     reads its fields as unformatted records of NUGP reals with no dimension
!     check -- so at any other rung it either took a truncated field or hit an
!     I/O error, and could not tell which. resmod.f90 is GENERATED into the
!     build directory from the resolution the build was configured for, which is
!     precisely so a unit cannot be read as another rung's.
!
!     NOTHING BUILDS THIS YET. The offline glacier accelerator is not in
!     `config/pipeline.yaml` and is not in plasim/CMakeLists.txt's source list;
!     `notes/audits/dormant-exoplasim-modules.md` has the facility and CLIM-53
!     owns the verdict on whether it gets wired. Deleting it would settle that
!     verdict rather than answer this one, so the grid is fixed and the decision
!     is left where it lives. world-zsa.

      integer, parameter :: NLAT = NLAT_ATM
      integer, parameter :: NLON = NLAT + NLAT
      integer, parameter :: NUGP = NLAT*NLON

      real :: dsnowz(NUGP) = 0.
      real :: lsm(NUGP) = 0.
      
      open(40,file='newdsnow',form='unformatted')
      read(40) dsnowz(:)
      close(40)
      
      open(42,file='lsm',form='unformatted')
      read(42) lsm(:)
      close(42)
      
      where (lsm .gt. 0.5)
         dsnowz(:) = min(max(dsnowz(:) + 400.0,0.0),3.0e3)
      endwhere
      
      open(17,file='newdsnow',form='unformatted')  
      write(17) dsnowz(:)
      close(17)
      
      stop
      end
      