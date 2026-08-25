!     ==================================================================
!     Does the threaded build's spectral analysis return the right answer
!     when the globe is cut into thread bands?
!
!     Worldbuilding frame: a correctness check on the Vesper climate
!     model's spectral transform. Nothing here is about the simulated
!     planet.
!
!     WHY THIS EXISTS. verify_threaded_numerics.sh was written as the
!     threaded build against an INDEPENDENT build, and world-38b left
!     one build, so its reference arm had nothing to compare against and
!     refused. world-d5l is the decision that replaced it: build a
!     standalone independent driver rather than keep only the
!     deterministic control arm or retire the gate. This is that driver.
!
!     WHAT IS THE SUBJECT AND WHAT IS THE REFERENCE, because everything
!     this links is the model's own source and that has to be said
!     plainly.
!
!     THE SUBJECT is the model's analysis chain, compiled from the model
!     tree and run at NPRO threads over a shared address space exactly as
!     the model runs it: `inigau` for the quadrature, `mpscdn` for the
!     scatter that gives a thread its latitudes, `legini` for the weight
!     matrices over those latitudes, `assoc_grid` for the band pointer
!     into the shared globe, `fc2sp` for the partial sum over the band,
!     and `mpsumsc` and `mpgallsp` for the reduction that adds the bands
!     up. That whole chain is under test.
!
!     THE REFERENCE is a file, written by
!     exoplasim/scripts/banded_transform_reference.py, and NOTHING in it
!     is computed here. Its associated Legendre functions come from
!     scipy, a different implementation of a different formulation; its
!     quadrature is the textbook construction in extended precision,
!     pinned by two identities that involve no method at all. That file's
!     header states what is independent, what is not, and what the choice
!     makes invisible. Read it before reading a result off this.
!
!     THE RIGHT ANSWER, and this is why it is a check rather than a
!     second opinion. The driver is handed Fourier coefficients that the
!     REFERENCE synthesised from a known spectral vector, and the model's
!     analysis of them must return that vector, because the reference
!     functions are orthonormal under the quadrature. Handing the model
!     coefficients IT synthesised would let an error in the weight
!     matrices cancel between the two halves and the round trip would
!     hold on a wrong table; that is the failure this arm is built to be
!     immune to, and it is why the synthesis is not done here.
!
!     WHAT THIS DOES NOT COVER, stated so that a pass is not read as more
!     than it is. It checks the ANALYSIS direction, the tables the
!     synthesis direction reads, and `dv2uv`. It does not integrate the
!     model, so it says nothing about the physics, about `mkdheat`, or
!     about any defect that needs a timestep to appear. `sp2fc` and
!     `sp2fcdmu` are not driven here: they are the same matrix-vector
!     product over the same rows that ARM A compares directly, they write
!     the thread's band the way `dv2uv` does and ARM D exercises that
!     placement, and their loops are checked one mode at a time by
!     exoplasim/scripts/verify_inverse_transform.py.
!
!     THE ARMS.
!
!     ARM A. The tables the thread ended up holding, against the
!     reference rows for the GLOBAL latitudes that thread owns. A thread
!     that received the wrong latitudes fails here even though its
!     arithmetic is perfect, which is the band question asked directly.
!
!     ARM B. The analysis itself, end to end, across the bands and
!     through the reduction, against the driving vector.
!
!     ARM C. ARM B repeated, and every repeat must agree BIT FOR BIT. A
!     reduction that summed in arrival order rather than thread order
!     gives a different answer run to run rather than a wrong one, so a
!     tolerance cannot see it and only an equality can. There is no
!     threshold here to choose.
!
!     ARM D. `dv2uv`, across the bands, against the reference's own wind
!     fields. This is the SYNTHESIS direction and it asks a different
!     question from ARM A: ARM A says the thread holds the right rows of
!     the weight matrices, ARM D says the routine puts the answer built
!     from those rows on the right rows of the shared wind globe. Nothing
!     in ARM A can see a loop that reads `qmat(2,l)` and writes latitude
!     l', and at one thread the two indices are the same integer, so
!     `verify_inverse_transform.py` cannot see it either. world-siq.
!
!     PLANETARY VORTICITY IS THE PART THAT NEEDED THIS. `dv2uv` takes
!     ABSOLUTE vorticity and removes the planetary part from its RESULT,
!     per latitude, at mode w=2 -- it cannot subtract it from the input,
!     which aliases the shared spectral state every thread reads. The
!     reference is not the model's formula for that removal but the
!     identity underneath it: `dv2uv(pd, pz)` with `plavor = P` is
!     `dv2uv(pd, pz - P at w=2)` with `P = 0`, so the reference computes
!     the winds from RELATIVE vorticity with no planetary term at all and
!     the driver hands the model that vorticity plus P. The mode, the
!     component, the sign and the factor are pinned by where P is put in,
!     and none of them is read off the model.
!
!     ONE HALF OF THAT REMOVAL IS DEAD ARITHMETIC AND NO ARM CAN GIVE IT
!     TEETH. Mode w=2 is m = 0, and `fmu` carries a factor of m, so
!     `fmu(2)` is zero and the model's `pv` planetary line adds exactly
!     nothing for any P. It is stated here so that a control on that line
!     is not written and read as passing.
!
!     THE BOUNDS. One arm has none, and the rest are derived from the
!     arithmetic or transplanted with their derivation named. None was
!     chosen after seeing what this produced.
!
!     ARM A1 HAS NO BOUND. Which global latitude a thread holds is an
!     INTEGER, so the check is which reference node each of the thread's
!     rows matches most closely, and the answer has to be the row the
!     decomposition says it owns. A threshold there would be answering a
!     different question than the one being asked.
!
!     THE NODE AND WEIGHT VALUES ARE NOT THIS GATE'S SUBJECT, and the
!     weights are compared here only because they reach the analysis
!     through the same scatter. `inigau` against an extended-precision
!     reference is exoplasim/scripts/verify_gauss_weights.sh, which
!     derives its bars and carries its own control, and ztolw is that
!     gate's weight bar rather than a second one invented here: 1e-12
!     relative, which its header derives as an order above the 2*n*eps
!     floor of evaluating P_n by the three-term recurrence. Two gates
!     with two bars on one quantity is how a quantity ends up with two
!     answers. The node VALUES are left to that gate entirely; ARM A1
!     uses the nodes to ask about placement and not about accuracy.
!
!     ztolc, gwd/cos^2, RELATIVE: ztolw plus eps/(1-mu^2) at the
!     pole-most node. `legini` forms cos^2 as 1 - sid*sid, and at the
!     pole-most latitude of the top rung sid is 0.99997, so that
!     subtraction loses four digits and the loss is what the second term
!     is. Both parts are the reference's own numbers, not the model's.
!
!     ztolp, the associated Legendre functions, ABSOLUTE for pmat and
!     scaled by the table's own largest magnitude for qmat:
!     4*NTRU^2*eps. THE FORM comes from the recurrence: `legini` walks an
!     NTRU-long chain in m and, inside it, an NTRU-long chain in n, and
!     each step is a subtract-and-scale contributing eps to values of
!     order one, so the compounded error is bounded by the product of the
!     two chain lengths. THE COEFFICIENT is an envelope rather than a
!     fit: legini's own recurrence differs from the reference table by
!     between 0.12*NTRU^2*eps and 0.35*NTRU^2*eps over the declared
!     ladder, so four clears every rung on it by between eleven and
!     thirty-three times. The ratio does not rise monotonically -- it is
!     0.16 at T21, dips to 0.12 at T106, and is 0.34 at T127 and 0.35 at
!     T170 -- so the envelope is a measured maximum and not a trend the
!     next rung can be read off. The rung-by-rung table is in
!     exoplasim/notes/banded-transform-reference.md, measured in Python
!     from legini's own recurrence and needing no build.
!
!     THE ENVELOPE IS OVER lib/rungs.py AND NOT BEYOND IT, and it is
!     still climbing at the top of the ladder. A rung above T170 does not
!     inherit this bound: re-measure legini's table error against the
!     reference at that rung, extend the table in the note, and either
!     confirm four still clears it or re-derive the coefficient, BEFORE
!     this gate is run there. Until that is done the failure at a new top
!     rung is a bound that was never argued for it, not a defective
!     build.
!
!     ztolf, the per-mode factors, ABSOLUTE: 4*eps. `fmu` and `fmv` are
!     m/(n(n+1)) and 1/(n(n+1)) at the unfiltered build this driver
!     refuses to run without, both sides form them by one division, and
!     neither exceeds one half. This is a rounding bar and not a
!     tolerance on a computation: it is here so that an ARM D failure
!     that is really a filtered build says so instead of reporting a
!     placement error.
!
!     ztold, dv2uv, ABSOLUTE: 4*NTP1*(eps*tscale + ztolp*(1+qscale)*
!     dvfmax). `tscale` and `dvfmax` are measured by the reference and
!     travel in the file: `tscale` is the largest sum of TERM MAGNITUDES
!     accumulated into any one output element, the planetary correction
!     included, and `dvfmax` is the largest |factor| times the largest
!     |coefficient| over the driving fields. THE FIRST TERM is the
!     accumulation's own rounding: at most NTP1 modes contribute and each
!     partial sum is bounded by `tscale`. THE SECOND is the weight
!     matrices the model built differing from the reference's by ARM A's
!     own bound -- ztolp absolute on `pmat` and ztolp*qscale on `qmat`,
!     which is what ARM A2 compares -- carried through NTP1 modes and
!     multiplied by the factor and coefficient each meets. The four is
!     the envelope: two terms per output element, doubled again, the same
!     construction ztolb uses.
!
!     ztold IS DOMINATED BY THE TABLE ENVELOPE AND IS LOOSE AGAINST WHAT
!     ARM D LOOKS FOR, and that is worth saying rather than leaving to be
!     noticed. ztolp is itself an envelope eleven to thirty-three times
!     above legini's measured table error, so ztold lands several orders
!     above what a correct build produces. The defect class ARM D exists
!     for is a term applied at the wrong latitude, which changes the
!     result by a whole term rather than by a rounding, and the margin is
!     printed beside the result so the comparison is visible rather than
!     asserted.
!
!     ztolb, the analysis, ABSOLUTE: 8*NLAT*eps*(1+pscale), where pscale
!     is the largest |pmat| in the reference table. The analysis is a sum
!     of NLAT terms each at most max|pmat| * gwd(l) * max|fc|; the
!     reference normalises max|fc| to one and the weights sum to two, so
!     the summands total at most 2*pscale and NLAT roundings of that is
!     2*NLAT*eps*pscale. The second term is the quadrature's own relative
!     error multiplying a result of order one. The whole is doubled
!     again, twice, to give the factor of eight.
!
!     WHAT DISCREPANCY THIS CAN SEE, against what the historical defects
!     were worth. ARM B's floor is ztolb, which is 2e-13 at T21 and 1e-12
!     at T170, on a result of order one. The weight pre-scaling defect
!     was a quadrature weight in the wrong place, which is an order-one
!     relative error: twelve orders clear. The `inigau` error that
!     `gaussmod.f90` records was 9.9e-11 relative, which is two orders
!     clear of ztolw. A band that overlaps by one latitude row changes
!     the sum by order one, and a band placed wrongly at all is caught by
!     ARM A1 with no margin to argue about. All of them are visible by a
!     wide margin, and the margin is printed beside the result rather
!     than asserted here.
!     ==================================================================

      module bandref
      implicit none

!     The reference, read once by the master thread before the parallel
!     region and read-only inside it. SHARED: nothing in this module may
!     become threadprivate, it is the second side of the comparison and
!     one copy is the point.
      integer :: nrlat, nrlon, nrtru, nrcsp, nrpro, nrlpp, nrlev
      real (kind=8) :: rscale, rversion, rplavor, tscale, dvfmax
      real (kind=8) :: pscale, qscale
      real (kind=8), allocatable :: rsid(:)      ! (nlat)
      real (kind=8), allocatable :: rgwd(:)      ! (nlat)
      real (kind=8), allocatable :: rpmat(:,:)   ! (ncsp,nlat)
      real (kind=8), allocatable :: rqmat(:,:)   ! (ncsp,nlat)
      real (kind=8), allocatable :: rspdrv(:)    ! (2*ncsp)
      real (kind=8), allocatable :: rfc(:)       ! (nlon*nlat)
!     ARM D. The synthesis half: the per-mode factors, the two spectral
!     fields dv2uv is driven with, and the winds it has to produce.
!     rspz is ABSOLUTE vorticity; rpu and rpv were computed from the
!     RELATIVE part with no planetary term, which is what makes them a
!     reference for that term rather than a copy of the model's formula.
      real (kind=8), allocatable :: rfmu(:)      ! (ncsp)
      real (kind=8), allocatable :: rfmv(:)      ! (ncsp)
      real (kind=8), allocatable :: rspd(:,:,:)  ! (2,ncsp,nlev)
      real (kind=8), allocatable :: rspz(:,:,:)  ! (2,ncsp,nlev)
      real (kind=8), allocatable :: rpu(:,:,:,:) ! (2,nlon/2,nlat,nlev)
      real (kind=8), allocatable :: rpv(:,:,:,:) ! (2,nlon/2,nlat,nlev)

!     Per-thread worst differences, so that one thread's verdict is not
!     written over another's and the master prints a table rather than
!     the threads interleaving lines.
      real (kind=8), allocatable :: wgwd(:), wgwc(:)
      real (kind=8), allocatable :: wpmt(:), wqmt(:)
      real (kind=8), allocatable :: wdvu(:), wdvv(:), wfac(:)
      integer, allocatable :: wrep(:), wlat(:)

!     The declared bounds. Set once, before any arm runs, and printed.
!     ztolw is verify_gauss_weights.sh's weight bar, not a second one.
      real (kind=8), parameter :: ZWEIGHTBAR = 1.0e-12_8
      real (kind=8) :: ztolw, ztolc, ztolp, ztolb, ztolf, ztold

!     The verdict. Raised under !$omp atomic from whichever thread sees a
!     failure; never lowered.
      integer :: nbad = 0

      integer, parameter :: NREP = 4   ! ARM C repeats
      end module bandref


      program bandcheck
      use pumamod
      use bandref
      implicit none

      character (len=256) :: yfile
      real (kind=8) :: zhead(12)
      real (kind=8) :: zeps
      integer :: io, j

      call get_command_argument(1,yfile)
      if (len_trim(yfile) == 0) then
         write(*,*) 'usage: verify_banded_transform.x <reference.bin>'
         write(*,*) 'the reference is written by'
         write(*,*) '  exoplasim/scripts/banded_transform_reference.py'
         stop 2
      endif

      open(71,file=trim(yfile),access='stream',form='unformatted',            &
     &     status='old',iostat=io)
      if (io /= 0) then
         write(*,*) 'cannot open the reference: ', trim(yfile)
         stop 2
      endif
      read(71) zhead
      nrlat = nint(zhead(1))
      nrlon = nint(zhead(2))
      nrtru = nint(zhead(3))
      nrcsp = nint(zhead(4))
      nrpro = nint(zhead(5))
      nrlpp = nint(zhead(6))
      rscale   = zhead(7)
      rversion = zhead(8)
      nrlev    = nint(zhead(9))
      rplavor  = zhead(10)
      tscale   = zhead(11)
      dvfmax   = zhead(12)

!     THE REFUSAL THAT STOPS A MISMATCHED FILE PASSING. A reference built
!     for another rung or another thread count would still read, still
!     fill arrays of some size, and still produce differences that mean
!     nothing. Every dimension is compared against the parameter this
!     executable was compiled with, and a disagreement is fatal here
!     rather than a number further down.
      if (nrlat /= NLAT .or. nrlon /= NLON .or. nrtru /= NTRU .or.            &
     &    nrcsp /= NCSP .or. nrpro /= NPRO .or. nrlpp /= NLPP .or.            &
     &    nrlev /= NLEV) then
         write(*,*) 'refusing: the reference does not describe this build.'
         write(*,'(a,7i7)') '  reference NLAT NLON NTRU NCSP NPRO NLPP NLEV:',&
     &                      nrlat, nrlon, nrtru, nrcsp, nrpro, nrlpp, nrlev
         write(*,'(a,7i7)') '  compiled                                    :',&
     &                      NLAT, NLON, NTRU, NCSP, NPRO, NLPP, NLEV
         stop 2
      endif
      if (nint(rversion) /= 2) then
         write(*,*) 'refusing: reference format version', rversion
         stop 2
      endif

      allocate(rsid(NLAT), rgwd(NLAT))
      allocate(rpmat(NCSP,NLAT), rqmat(NCSP,NLAT))
      allocate(rspdrv(2*NCSP), rfc(NLON*NLAT))
      allocate(rfmu(NCSP), rfmv(NCSP))
      allocate(rspd(2,NCSP,NLEV), rspz(2,NCSP,NLEV))
      allocate(rpu(2,NLON/2,NLAT,NLEV), rpv(2,NLON/2,NLAT,NLEV))
      read(71) rsid
      read(71) rgwd
      read(71) rpmat
      read(71) rqmat
      read(71) rspdrv
      read(71) rfc
      read(71) rfmu
      read(71) rfmv
      read(71) rspd
      read(71) rspz
      read(71) rpu
      read(71) rpv
      close(71)

      allocate(wgwd(0:NPRO-1), wgwc(0:NPRO-1))
      allocate(wpmt(0:NPRO-1), wqmt(0:NPRO-1))
      allocate(wdvu(0:NPRO-1), wdvv(0:NPRO-1), wfac(0:NPRO-1))
      allocate(wrep(0:NPRO-1), wlat(0:NPRO-1))
      wgwd = 0.0_8 ; wgwc = 0.0_8
      wpmt = 0.0_8 ; wqmt = 0.0_8 ; wrep = 0 ; wlat = 0
      wdvu = 0.0_8 ; wdvv = 0.0_8 ; wfac = 0.0_8

!     The two magnitudes the bounds are stated in terms of, taken from
!     the REFERENCE and not from anything the model produced.
      pscale = maxval(abs(rpmat))
      qscale = maxval(abs(rqmat))
      if (qscale <= 0.0_8) qscale = 1.0_8

      zeps  = epsilon(1.0_8)
      ztolw = ZWEIGHTBAR
      ztolc = ZWEIGHTBAR + zeps / (1.0_8 - maxval(rsid*rsid))
      ztolp = 4.0_8 * NTRU * NTRU * zeps
      ztolb = 8.0_8 * NLAT * zeps * (1.0_8 + pscale)
      ztolf = 4.0_8 * zeps
      ztold = 4.0_8 * NTP1 * (zeps * tscale                                   &
     &                      + ztolp * (1.0_8 + qscale) * dvfmax)

      write(*,'(a,i0,a,i0,a,i0,a,i0)')                                        &
     &   'T', NTRU, '  NLAT ', NLAT, '  modes ', NCSP, '  threads ', NPRO
      write(*,'(a)') 'declared before the arms ran:'
      write(*,'(a)')       '  band placement: an integer, no bound'
      write(*,'(a,es10.3)') '  weights, relative, verify_gauss_weights ', ztolw
      write(*,'(a,es10.3)') '  gw/cos2, relative, that + the cos2 loss ', ztolc
      write(*,'(a,es10.3)') '  P and Q, absolute, 4*NTRU^2*eps         ', ztolp
      write(*,'(a,es10.3)') '  analysis,absolute, 8*NLAT*eps*(1+Pmax)  ', ztolb
      write(*,'(a,es10.3)') '  fmu/fmv, absolute, 4*eps                ', ztolf
      write(*,'(a,es10.3)') '  dv2uv,   absolute, 4*NTP1*(...)         ', ztold
      write(*,'(a,es10.3)') '  largest |P| in the reference table      ', pscale
      write(*,'(a,es10.3)') '  largest |Q| in the reference table      ', qscale
      write(*,'(a,es10.3)') '  worst term sum into one dv2uv output    ', tscale
      write(*,'(a,es10.3)') '  largest |factor|*|coefficient|          ', dvfmax
      write(*,'(a,f8.3)')   '  planetary vorticity driven at           ', rplavor
      write(*,'(a,i0,a)')  '  ARM C repeats ', NREP, ', bit identical, no bound'
      write(*,*)

!$omp parallel
      call bandarm
!$omp end parallel

      write(*,*)
      if (nbad == 0) then
         write(*,'(a)') 'the banded analysis returns the reference answer.'
      else
         write(*,'(a,i0,a)') 'FAILED: ', nbad, ' arm(s) missed.'
         stop 1
      endif
      end program bandcheck


!     ==================================================================
!     One thread's share of every arm. Called from inside the parallel
!     region by every thread, which is what mpstart and the reductions
!     below require: each of them carries a barrier and a barrier every
!     thread does not reach is a hang.
!     ==================================================================

      subroutine bandarm
      use pumamod
      use legmod
      use bandref
      implicit none

      integer :: jl, jg, jw, jrep, j, it, jm, jv, ip
      real :: zpart(NESP)
      real :: zslice(NSPP)
      real :: zfull(NESP)
      real :: zkeep(NESP)
      real (kind=8) :: zew, zec, zep, zeq, zeb, zref, zbest
      real (kind=8) :: zdu, zdv, zdf
      integer :: idiff, iworld, imiss, inear

      iworld = 0
      call mpstart(iworld)

!     nfilter MUST be 0, and it is asserted rather than assumed. `fc2sp`
!     carries skgpsp in fgp, so a filtered analysis returns the
!     coefficient times the filter and the reference's right answer is
!     not the model's. The default is 0 and this driver reads no
!     namelist, so this fires only if the declaration changed.
      if (nfilter /= 0) then
         if (mypid == NROOT) then
            write(*,*) 'refusing: nfilter is', nfilter, 'and must be 0.'
            write(*,*) 'A filtered analysis returns the coefficient times'
            write(*,*) 'the filter, which is not what the reference states.'
         endif
!$omp barrier
         stop 2
      endif

!     The model's own grid setup, in the model's own order: the master
!     builds the global quadrature, the scatter hands each thread its
!     band of it, and legini builds the weight matrices from whatever
!     that thread received. Everything in these four lines is subject,
!     not reference.
      if (mypid == NROOT) then
         call inigau(NLAT,sid,gwd)
      endif
      call mpscdn(sid,NLPP)
      call mpscdn(gwd,NLPP)
      call legini

!     ------- ARM A1: which latitudes did this thread actually get? ----
!
!     Asked as an integer rather than as a difference. Each of the
!     thread's rows is matched against every reference node and the
!     nearest one wins; the winner has to be the global latitude the
!     decomposition says this thread owns. The nodes are separated by far
!     more than any rounding, so the nearest node is the identity of the
!     row and not a measurement of it, and a thread working the wrong
!     part of the globe with perfect arithmetic fails here.
      imiss = 0
      do jl = 1 , NLPP
         jg = mypid * NLPP + jl
         inear = 1
         zbest = abs(sid(jl) - rsid(1))
         do j = 2 , NLAT
            if (abs(sid(jl) - rsid(j)) < zbest) then
               zbest = abs(sid(jl) - rsid(j))
               inear = j
            endif
         enddo
         if (inear /= jg) imiss = imiss + 1
      enddo
      wlat(mypid) = imiss

!     ---------------- ARM A2: the tables this thread holds ------------
      zew = 0.0_8 ; zec = 0.0_8
      zep = 0.0_8 ; zeq = 0.0_8
      do jl = 1 , NLPP
         jg = mypid * NLPP + jl
         zew = max(zew, abs((gwd(jl) - rgwd(jg)) / rgwd(jg)))
!        gwd/cos^2, against the reference's own node and weight rather
!        than against the model's, so this is not the model's formula
!        checked against itself.
         zref = rgwd(jg) / (1.0_8 - rsid(jg) * rsid(jg))
         zec = max(zec, abs((gwdc(jl) - zref) / zref))
         do jw = 1 , NCSP
            zep = max(zep, abs(pmat(jw,jl) - rpmat(jw,jg)))
            zeq = max(zeq, abs(qmat(jw,jl) - rqmat(jw,jg)) / qscale)
         enddo
      enddo
      wgwd(mypid) = zew
      wgwc(mypid) = zec
      wpmt(mypid) = zep
      wqmt(mypid) = zeq
!$omp barrier

      if (mypid == NROOT) then
         write(*,'(a)') '==== ARM A1: the latitudes each thread got ===='
         imiss = sum(wlat)
         if (imiss /= 0) then
            do it = 0 , NPRO-1
               write(*,'(a,i0,a,i0,a,i0)') '    thread ', it, ': ', wlat(it), &
     &              ' of ', NLPP, ' rows are not the ones it owns'
            enddo
            write(*,'(a)') '  [ FAIL ] the bands are not the bands the'   //  &
     &                     ' decomposition describes.'
            call raisebad
         else
            write(*,'(a,i0,a)') '  [  ok  ] every thread holds its own '  //  &
     &           'NLPP rows, all ', NLAT, ' accounted for.'
         endif
         write(*,*)
         write(*,'(a)') '==== ARM A2: the tables each thread holds ===='
         write(*,'(a)') '  thread    weight     gw/cos2        P     '    //  &
     &                  '      Q'
         do it = 0 , NPRO-1
            write(*,'(i8,4es12.3)') it, wgwd(it), wgwc(it),                    &
     &                             wpmt(it), wqmt(it)
         enddo
         call judge('ARM A weights',   maxval(wgwd), ztolw)
         call judge('ARM A gw/cos2',   maxval(wgwc), ztolc)
         call judge('ARM A P      ',   maxval(wpmt), ztolp)
         call judge('ARM A Q      ',   maxval(wqmt), ztolp)
      endif
!$omp barrier

!     ------- ARM B and ARM C: the analysis, and its repeatability -----
!
!     The Fourier coefficients go into the SHARED globe and each thread
!     reads its band of them through the pointer assoc_grid gave it. That
!     is deliberate and it is what puts the band arithmetic itself under
!     test: a band that is slid, overlapped or short reads the wrong rows
!     here while holding the right weight matrices, which is exactly the
!     mistake this decomposition invites.
!
!     gp is NHOR long and fc2sp's dummy is fc(2,NLON/2,NLPP), which is
!     the same NHOR elements in the same order.
      idiff = 0
      do jrep = 1 , NREP
         if (mypid == NROOT) then
            do j = 1 , NUGP
               gp_g(j) = rfc(j)
            enddo
         endif
!$omp barrier
!        Cleared here and not left to fc2sp: it zeroes sp(2,NESP/2), which
!        is one element short of NESP when NPRO does not divide NRSP, and
!        mpsumsc reads the whole of it.
         do j = 1 , NESP
            zpart(j) = 0.0
         enddo
         call fc2sp(gp,zpart)
         call mpsumsc(zpart,zslice,1)
         call mpgallsp(zfull,zslice,1)
         if (jrep == 1) then
            do j = 1 , NESP
               zkeep(j) = zfull(j)
            enddo
         else
            do j = 1 , NESP
               if (zfull(j) /= zkeep(j)) idiff = idiff + 1
            enddo
         endif
      enddo

      zeb = 0.0_8
      do j = 1 , 2*NCSP
         zeb = max(zeb, abs(real(zkeep(j),8) - rspdrv(j)))
      enddo

      if (mypid == NROOT) then
         write(*,*)
         write(*,'(a)') '==== ARM B: the analysis, across the bands ===='
         write(*,'(a,es12.3)') '  worst difference from the driving vector ',  &
     &                        zeb
         call judge('ARM B        ', zeb, ztolb)
         write(*,*)
         write(*,'(a,i0,a)') '==== ARM C: ', NREP,                            &
     &                       ' repeats, bit identical ===='
      endif
!$omp barrier

!     Every thread reports its own count and the master reads the table:
!     a reduction that races may do so on one thread's slice and not
!     another's, so a single shared counter would lose which.
      wrep(mypid) = idiff
!$omp barrier
      if (mypid == NROOT) then
         idiff = sum(wrep)
         if (idiff /= 0) then
            do it = 0 , NPRO-1
               write(*,'(a,i0,a,i0,a)') '    thread ', it, ': ', wrep(it),    &
     &              ' coefficients changed between repeats'
            enddo
            write(*,'(a)') '  [ FAIL ] the answer depends on the run, so'//   &
     &                     ' the reduction is not summing in a fixed'
            write(*,'(a)') '           thread order.'
            call raisebad
         else
            write(*,'(a)') '  [  ok  ] every repeat agreed to the last bit.'
         endif
      endif

!     -------------- ARM D: dv2uv, across the bands -------------------
!
!     The model's own call, spelled the way `invlegd` spells it: the
!     shared spectral state in, the thread's POINTER into the shared wind
!     globe out. Every thread transforms the whole spectrum and writes
!     only its own band, so a routine that reads the right weight rows
!     and writes the wrong grid rows fails here and nowhere else.
!
!     THE ARRAY TEMPORARY THE CHECKED BUILD REPORTS HERE IS THE MODEL'S
!     OWN. `gu` points at a band of the shared globe, which is a strided
!     section, and `dv2uv` takes an explicit-shape dummy, so the compiler
!     packs it in and out. `plasim.f90` calls `dv2uv(sd,sz,gu,gv)` the
!     same way and gets the same temporary; it is what
!     probe_grid_contiguity.f90 exists to decide by address rather than
!     by behaviour. Two lines of runtime warning are the flag profile
!     saying so, not a defect in this driver.
!
!     plavor is threadprivate, so every thread sets it. It is the
!     reference's declared value rather than the planet's: what is under
!     test is where the term lands.
      plavor = real(rplavor)

!     The factors first, and separately, because an ARM D failure that is
!     really a filtered build has a different fix from a placement error
!     and must not be reported as one.
      zdf = 0.0_8
      do jw = 1 , NCSP
         zdf = max(zdf, abs(real(fmu(jw),8) - rfmu(jw)))
         zdf = max(zdf, abs(real(fmv(jw),8) - rfmv(jw)))
      enddo
      wfac(mypid) = zdf

      if (mypid == NROOT) then
         sd(:,:) = 0.0
         sz(:,:) = 0.0
         do jv = 1 , NLEV
            do jw = 1 , NCSP
!              sd(NESP,NLEV) is dv2uv's pd(2,NESP/2,NLEV) by sequence
!              association, so mode jw is elements 2*jw-1 and 2*jw.
               sd(2*jw-1,jv) = real(rspd(1,jw,jv))
               sd(2*jw  ,jv) = real(rspd(2,jw,jv))
               sz(2*jw-1,jv) = real(rspz(1,jw,jv))
               sz(2*jw  ,jv) = real(rspz(2,jw,jv))
            enddo
         enddo
         gu_g(:,:) = 0.0
         gv_g(:,:) = 0.0
      endif
!$omp barrier
      call dv2uv(sd,sz,gu,gv)
!$omp barrier

!     Read back out of the SHARED globe by global index rather than
!     through the pointer, so the comparison does not inherit whatever
!     the band pointer believes. Row mypid*NHOR + 2*(m-1) + NLON*(l-1)
!     is dv2uv's pu(:,m,l,v) and belongs to global latitude
!     mypid*NLPP + l.
      zdu = 0.0_8 ; zdv = 0.0_8
      do jv = 1 , NLEV
         do jl = 1 , NLPP
            jg = mypid * NLPP + jl
            do jm = 1 , NLON/2
               ip = mypid * NHOR + 2*(jm-1) + NLON*(jl-1)
               zdu = max(zdu, abs(real(gu_g(ip+1,jv),8) - rpu(1,jm,jg,jv)))
               zdu = max(zdu, abs(real(gu_g(ip+2,jv),8) - rpu(2,jm,jg,jv)))
               zdv = max(zdv, abs(real(gv_g(ip+1,jv),8) - rpv(1,jm,jg,jv)))
               zdv = max(zdv, abs(real(gv_g(ip+2,jv),8) - rpv(2,jm,jg,jv)))
            enddo
         enddo
      enddo
      wdvu(mypid) = zdu
      wdvv(mypid) = zdv
!$omp barrier

      if (mypid == NROOT) then
         write(*,*)
         write(*,'(a)') '==== ARM D: dv2uv, across the bands ===='
         write(*,'(a)') '  thread     fmu/fmv           u           v'
         do it = 0 , NPRO-1
            write(*,'(i8,3es12.3)') it, wfac(it), wdvu(it), wdvv(it)
         enddo
         call judge('ARM D fmu/fmv', maxval(wfac), ztolf)
         call judge('ARM D u      ', maxval(wdvu), ztold)
         call judge('ARM D v      ', maxval(wdvv), ztold)
      endif
!$omp barrier

      return
      end subroutine bandarm


!     A failure that has no bound to report: ARM A1 and ARM C.
      subroutine raisebad
      use bandref
      implicit none
!$omp atomic
      nbad = nbad + 1
      return
      end subroutine raisebad


!     One arm's verdict, so that the comparison against the declared
!     bound is written once and the same way for all of them.
      subroutine judge(ywhat,pgot,ptol)
      use bandref
      implicit none
      character (len=*) :: ywhat
      real (kind=8) :: pgot, ptol
      if (pgot > ptol) then
         write(*,'(a,a,es11.3,a,es11.3)') '  [ FAIL ] ', ywhat, pgot,           &
     &        '  over the declared bound ', ptol
!$omp atomic
         nbad = nbad + 1
      else
         write(*,'(a,a,es11.3,a,es11.3)') '  [  ok  ] ', ywhat, pgot,           &
     &        '  within                  ', ptol
      endif
      return
      end subroutine judge
