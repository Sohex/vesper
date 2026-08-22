      program roundtrip
!     Does the transform pair return what it was given?
!
!     Worldbuilding frame: a correctness check on the Vesper climate model's
!     spectral transform. Nothing here is about the simulated planet.
!
!     THE IDENTITY. Analysis composed with synthesis is the identity on a
!     Gaussian grid, exactly, because the quadrature is exact for the
!     polynomials the truncation admits. So
!
!         sp -> sp2fc -> fc2gp -> gp2fc -> fc2sp -> sp
!
!     must return the coefficients it started with, to rounding. That is a
!     RIGHT ANSWER and not a comparison: it needs no reference build, no second
!     binary, and no previous run.
!
!     WHY THAT MATTERS MORE THAN IT LOOKS. Every other check on this transform
!     compares one build against another -- threaded against MPI, factored
!     against eight matrices -- and all of them stop working the moment the
!     transform is REPLACED rather than rearranged. SHTns computes a different
!     algorithm, so nothing it produces will match legmod at rounding scale, and
!     the comparison that has caught every defect so far will simply go quiet.
!     This check does not care which implementation is underneath. It is the one
!     that survives the swap, which is why it is written before the swap.
!
!     NFILTER MUST BE 0, which is the default. The physics filter is folded into
!     the weights, so a filtered round trip returns the coefficient times
!     skspgp(n)*skgpsp(n) rather than the coefficient. That is not a defect and
!     this check would report it as one, so the filter is asserted rather than
!     assumed.
      use pumamod
      use legmod
      implicit none

      real (kind=8) :: zsi(NLAT), zgw(NLAT)
      real :: zspin(2,NCSP), zspout(2,NCSP)
      real :: zgrid(NLON,NLPP)
      integer :: j, jm, jn, lm, nbad
      real :: zerr, zworst, ztol
      character(len=16) :: zwhat

      ztol = 1.0e-11
      nbad = 0
      zworst = 0.0

      call inigau(NLAT,zsi,zgw)
      do j = 1 , NLPP
         sid(j) = zsi(j)
         gwd(j) = zgw(j)
      enddo
      if (nfilter /= 0) then
         write(*,*) 'refusing: nfilter is', nfilter, 'and must be 0 -- a'
         write(*,*) 'filtered round trip is not the identity, and this check'
         write(*,*) 'would report the filter as a defect.'
         stop 2
      endif
      call legini

      write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT, '  modes ', NCSP
      write(*,'(a,e9.2,a)') 'tolerance ', ztol, ', declared before the arms ran'
      write(*,*)

!     1. every mode on its own. A transform can be right on a dense spectrum
!        and wrong on one coefficient -- the planetary vorticity mode was
!        exactly that -- so each is driven separately and has to come back.
      do lm = 1 , NCSP
         zspin(:,:) = 0.0
         zspin(1,lm) = 1.0
!        The m=0 coefficients are REAL: the zonal mean has no imaginary part,
!        and the transform enforces that by discarding it. Driving one and
!        expecting it back is testing the transform for a property it is right
!        not to have. The first NTP1 modes are m=0 in this ordering, and the
!        check found this itself -- 43 of 946 failing at T42 is exactly NTP1.
         if (lm > NTP1) zspin(2,lm) = -0.5
         call sp2fc(zspin,zgrid)
         call fc2gp(zgrid,NLON,NLPP)
         call gp2fc(zgrid,NLON,NLPP)
         call fc2sp(zgrid,zspout)
         zerr = maxval(abs(zspout - zspin))
         if (zerr > zworst) then
            zworst = zerr
            write(zwhat,'(a,i0)') 'mode ', lm
         endif
         if (zerr > ztol) nbad = nbad + 1
      enddo
      if (nbad == 0) then
         write(*,'(a,i0,a)') '  [  ok  ] all ', NCSP,                           &
     &      ' modes return themselves through the round trip'
      else
         write(*,'(a,i0,a,i0,a)') '  [ FAIL ] ', nbad, ' of ', NCSP,            &
     &      ' modes do not come back'
      endif

!     2. a dense spectrum, which is what the model actually transforms, and
!        which exercises the accumulation rather than one term of it.
      do j = 1 , NCSP
         zspin(1,j) =  1.0 / real(j)
         zspin(2,j) = -0.5 / real(j)
      enddo
      zspin(2,1:NTP1) = 0.0          ! m=0 is real, as above
      call sp2fc(zspin,zgrid)
      call fc2gp(zgrid,NLON,NLPP)
      call gp2fc(zgrid,NLON,NLPP)
      call fc2sp(zgrid,zspout)
      zerr = maxval(abs(zspout - zspin))
      if (zerr > zworst) then
         zworst = zerr
         zwhat = 'dense spectrum'
      endif
      if (zerr <= ztol) then
         write(*,'(a,e11.4)') '  [  ok  ] a dense spectrum returns itself, worst ', zerr
      else
         write(*,'(a,e11.4)') '  [ FAIL ] a dense spectrum does not, worst ', zerr
         nbad = nbad + 1
      endif

      write(*,*)
      write(*,'(a,e11.4,a,a)') 'worst anywhere: ', zworst, '  at ', zwhat
      if (nbad == 0) then
         write(*,'(a)') '0 failed'
      else
         write(*,'(i0,a)') nbad, ' failed'
         stop 1
      endif
      end program roundtrip
