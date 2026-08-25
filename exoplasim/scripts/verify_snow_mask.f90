!     ==================================================================
!     Does the canopy/snow mask reduce to the scalar it replaced, and does it
!     behave at the limits where the answer can be written down?
!
!     snowmaskmod generalises upstream's one-scalar forest snow mask along
!     burial, plant area and canopy interception. The check that can fail is
!     the reduction: at the values the model ships, the operator has to return
!     the cover itself, bit for bit, because the land module then evaluates
!     exactly the expression it evaluated before.
!
!     Every other check here has an answer written down in advance rather than
!     taken from another run: a treeless cell hides nothing, a canopy buried to
!     its top hides nothing, a canopy loaded with snow gives the albedo of open
!     snow whatever its cover, and the hidden fraction never leaves 0 to 1 and
!     never rises with snow depth or falls with plant area.
!
!     It also reports the plant area index upstream's two masking fractions
!     imply under a gap-fraction reading, which is the comparison BIO-30 asks
!     for: the fraction and the law are two statements of the same quantity, so
!     one can be read as the other.
!
!     Built and run by verify_snow_mask.sh.
!
!     Worldbuilding frame: a correctness check on the surface albedo of a
!     simulated planet's snow-covered vegetated land.
!     ==================================================================
      program verify_snow_mask
      use snowmaskmod
      implicit none

      real, parameter :: covmx = 0.6153846   ! landmod's shipped values
      real, parameter :: covmn = 0.4
      real, parameter :: off   = -1.0        ! the inert canopy height and PAI
      real, parameter :: ext   = 1.0         ! landmod's shipped extinction

      integer :: nfail
      nfail = 0

      call reduction(nfail)
      call limits(nfail)
      call bounds(nfail)
      call monotone(nfail)
      call implied_pai()

      write(*,*)
      if (nfail == 0) then
         write(*,*) "PASS: the mask reduces to the cover and holds at every limit"
      else
         write(*,'(a,i4)') " FAIL: checks failed: ", nfail
         stop 1
      endif

      contains

      subroutine mask(pcover,psnowm,phcan,ppai,pint,                    &
     &                pfcovmx,pfcovmn,pfint,pkmx,pkmn)
      real, intent(in)  :: pcover,psnowm,phcan,ppai,pint
      real, intent(out) :: pfcovmx,pfcovmn,pfint,pkmx,pkmn
      call snowcanopymask(pcover,psnowm,phcan,ppai,ext,covmx,covmn,pint,&
     &                    pfcovmx,pfcovmn,pfint,pkmx,pkmn)
      end subroutine mask

!     THE REDUCTION. At the shipped values the operator must return the cover
!     itself and a zero interception weight, for every cover and at every snow
!     depth, so that the land module's blend is unchanged to the last bit.
      subroutine reduction(kfail)
      integer, intent(inout) :: kfail
      integer :: j,k
      real :: zcov,zsnow,zfmx,zfmn,zfint,zkmx,zkmn
      do j=0,20
         zcov = real(j)/20.0
         do k=0,10
            zsnow = real(k)*0.5
            call mask(zcov,zsnow,off,off,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
            if (zfmx /= zcov .or. zfmn /= zcov .or. zfint /= 0.0) then
               write(*,'(a,3f12.8)') " FAIL reduction: cover, fmx, fmn = ",  &
     &                               zcov,zfmx,zfmn
               kfail = kfail + 1
            endif
         enddo
      enddo
      write(*,*) "reduction: the cover is returned unchanged at every depth"
      end subroutine reduction

!     THE LIMITS BIO-30 NAMES. Each answer is written down here, not measured.
      subroutine limits(kfail)
      integer, intent(inout) :: kfail
      real :: zfmx,zfmn,zfint,zkmx,zkmn

!     Treeless: nothing masks anything, whatever else is set.
      call mask(0.0,0.2,20.0,4.0,1.0,zfmx,zfmn,zfint,zkmx,zkmn)
      call expect(zkmx,0.0,"treeless cell hides no ground",kfail)

!     Low vegetation buried to its top: the canopy is under the snow.
      call mask(1.0,0.5,0.4,4.0,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
      call expect(zkmx,0.0,"buried low vegetation hides no ground",kfail)

!     Half-buried, uniform plant area in the vertical: half of it is left.
      call mask(1.0,1.0,2.0,off,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
      call expect(zkmx,0.5*covmx,"a half-buried canopy hides half as much",kfail)

!     Open tall canopy on shallow snow: the plant-area law decides, and a
!     canopy this dense hides nearly all of the ground under it.
      call mask(1.0,0.1,20.0,8.0,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
      if (zkmn < 0.97 .or. zkmn > 1.0) then
         write(*,'(a,f12.8)') " FAIL: dense open canopy hides ",zkmn
         kfail = kfail + 1
      else
         write(*,*) "limit: a dense exposed canopy hides nearly all of it"
      endif

!     Complete canopy snow: the forested endmember becomes exposed snow, so the
!     blend gives open snow whatever the cover. Checked as the weight, which is
!     what the land module blends with.
      call mask(0.3,0.2,20.0,off,1.0,zfmx,zfmn,zfint,zkmx,zkmn)
      call expect(zfint,1.0,"a fully loaded canopy is the snow endmember",kfail)
      end subroutine limits

!     The hidden fraction is a fraction of the cell and cannot leave 0 to 1,
!     whatever it is handed. A sweep over the whole input space.
      subroutine bounds(kfail)
      integer, intent(inout) :: kfail
      integer :: j,k,l,m
      real :: zcov,zsnow,zhcan,zpai,zfmx,zfmn,zfint,zkmx,zkmn
      integer :: nbad
      nbad = 0
      do j=0,10
       zcov = real(j)/10.0
       do k=0,10
        zsnow = real(k)*0.4
        do l=0,6
         zhcan = real(l)*4.0 - 1.0
         do m=0,8
          zpai = real(m)*1.5 - 1.0
          call mask(zcov,zsnow,zhcan,zpai,0.5,zfmx,zfmn,zfint,zkmx,zkmn)
          if (zkmx < 0.0 .or. zkmx > 1.0 .or. zkmn < 0.0 .or. zkmn > 1.0)   &
     &       nbad = nbad + 1
         enddo
        enddo
       enddo
      enddo
      if (nbad > 0) then
         write(*,'(a,i6,a)') " FAIL: ",nbad," input combinations leave 0 to 1"
         kfail = kfail + 1
      else
         write(*,*) "bounds: the hidden fraction stays in 0 to 1 everywhere"
      endif
      end subroutine bounds

!     Deeper snow can only bury more canopy, and more plant area can only hide
!     more ground. Both directions are what makes this an operator rather than
!     a fit.
      subroutine monotone(kfail)
      integer, intent(inout) :: kfail
      integer :: k
      real :: zfmx,zfmn,zfint,zkmx,zkmn,zlast
      logical :: lbad

      lbad = .false.
      zlast = 2.0
      do k=0,40
         call mask(0.8,real(k)*0.25,6.0,off,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
         if (zkmx > zlast) lbad = .true.
         zlast = zkmx
      enddo
      if (lbad) then
         write(*,*) " FAIL: the hidden fraction rises with snow depth"
         kfail = kfail + 1
      else
         write(*,*) "monotone: deeper snow never hides more ground"
      endif

      lbad = .false.
      zlast = -1.0
      do k=1,40
         call mask(0.8,0.1,20.0,real(k)*0.25,0.0,zfmx,zfmn,zfint,zkmx,zkmn)
         if (zkmn < zlast) lbad = .true.
         zlast = zkmn
      enddo
      if (lbad) then
         write(*,*) " FAIL: the hidden fraction falls with plant area"
         kfail = kfail + 1
      else
         write(*,*) "monotone: more plant area never hides less ground"
      endif
      end subroutine monotone

!     THE GAP-FRACTION COMPARISON. A masking fraction and a gap-fraction law
!     are two statements of the same quantity, so upstream's two constants can
!     be read as plant area indices: L = -ln(1 - cov)/k. Reported rather than
!     asserted, because what it is worth is a modelling verdict and not a
!     property of this code.
      subroutine implied_pai()
      real :: zk
      integer :: j
      write(*,*)
      write(*,*) "the plant area index the shipped masking fractions imply:"
      do j=0,3
         zk = 0.5 + real(j)*0.25
         write(*,'(a,f4.2,a,f6.3,a,f6.3)')                                 &
     &      "   extinction ",zk,":  bright-snow end ",-log(1.0-covmx)/zk,  &
     &      "   aged-snow end ",-log(1.0-covmn)/zk
      enddo
      end subroutine implied_pai

      subroutine expect(pgot,pwant,ytext,kfail)
      real, intent(in) :: pgot,pwant
      character (len=*), intent(in) :: ytext
      integer, intent(inout) :: kfail
      if (abs(pgot-pwant) > 1.0e-6) then
         write(*,'(a,a,a,f12.8,a,f12.8)') " FAIL: ",ytext," got ",pgot,     &
     &                                    " want ",pwant
         kfail = kfail + 1
      else
         write(*,'(a,a)') " limit: ",ytext
      endif
      end subroutine expect

      end program verify_snow_mask
