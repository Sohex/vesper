! ************************************************************************
! * module legmod - direct and indirect Legendre transformation routines *
! * E. Kirk 20-Feb-2009 tested for T21 - T341 resolutions 32 & 64 bit    *
! ************************************************************************

module legmod
! ************************
! * Legendre Polynomials *
! ************************
use pumamod, only:NTRU,NTP1,NCSP,NESP,NLON,NLPP,NLHP,NLAT,NHOR,NLEV,gwd,sid,plavor,nfilter
use pumamod, only:LPAIRLAT
use pumamod, only:ngptfilter, nspvfilter,landhoskn0,filterkappa,nfilterexp,nud,mypid,NROOT

! TWO MATRICES, NOT EIGHT.
!
! Every weight this module ever used is P or its mu-derivative Q times a factor
! that depends only on the MODE and a factor that depends only on the LATITUDE.
! Neither of those needs an NCSP x NLPP array. Storing the eight products cost
! 15.1 MB a thread at T170, 114.9 MB a die against CCD1's 32 MB of L3, and every
! one of them was streamed every timestep. Storing P and Q and applying the
! factors leaves 4.5 MB a thread and about 36 MB a die.
!
! It is the filter fold read backwards. That change moved a per-mode scalar INTO
! the matrices to save one multiply in three and returned 1.6%, because these
! loops are bound by streaming the matrices and not by the multiplier. If that
! is true then the multiply is nearly free and the streaming is what to cut.
!
! The per-mode vectors are NCSP long, are reused across every latitude, and stay
! in cache. The per-latitude factors are loop invariant in the mode loop and are
! hoisted into a scalar by the caller.
real :: pmat(NCSP,NLPP) ! P(m,n),      the associated Legendre function
real :: qmat(NCSP,NLPP) ! Q(m,n),      its derivative with respect to mu

! Per-mode factors. The comment on each names the matrix it replaces with the
! per-latitude factor listed beside it.
real :: fsp(NCSP) ! skspgp                     qi = pmat*fsp,   qj = qmat*fsp
real :: fgp(NCSP) ! skgpsp                     qc = pmat*fgp*gwd,  qe = qmat*fgp*gwdc
real :: fmu(NCSP) ! m/(n(n+1)) * skspgp        qu = pmat*fmu
real :: fmv(NCSP) ! 1/(n(n+1)) * skspgp        qv = qmat*fmv
real :: fmm(NCSP) ! m * skgpsp                 qm = pmat*fmm*gwdc
real :: fmq(NCSP) ! n(n+1)/2 * skgpsp          qq = pmat*fmq*gwdc

! Per-latitude factors. gwd is the Gaussian weight and lives in pumamod; this is
! the other one, and it is kept here because only these loops want it.
real :: gwdc(NLPP) ! gwd / cos2
! The physics filter, indexed by total wavenumber. legini carries it into the
! per-mode vectors above, which are indexed by MODE; these two are kept for the
! diagnostic print and are not read by the transform loops.
real :: skgpsp(NTP1) ! Physics filter for GP -> SP
real :: skspgp(NTP1) ! Physics filter for SP -> GP


!     Threads instead of ranks: a thread owns what a rank owned.
!     Inert without -fopenmp, so the MPI and serial builds are unchanged.
!     pmat and qmat stay threadprivate because their SECOND index is the
!     thread's own latitude, so the content differs per thread and sharing one
!     array would not save a byte. The per-mode vectors are the same on every
!     thread and could be shared; they are threadprivate too, because at 118 KB
!     each they are not worth a second storage rule, and because a thread
!     reading its own copy cannot false-share with a neighbour.
!$omp threadprivate(pmat,qmat,fsp,fgp,fmu,fmv,fmm,fmq,gwdc,skgpsp,skspgp)

end module legmod

! =================
! SUBROUTINE LEGINI
! =================

subroutine legini
use legmod
implicit none

integer :: jlat ! Latitude
integer :: lm
integer :: m
integer :: n

real (kind=8) :: amsq
real (kind=8) :: z1
real (kind=8) :: z2
real (kind=8) :: z3
real (kind=8) :: f1m
real (kind=8) :: f2m
real (kind=8) :: znn1
real (kind=8) :: zsin    ! sin
real (kind=8) :: zcsq    ! cos2
real (kind=8) :: zgwd    ! gw
real (kind=8) :: zgwdcsq ! gw / cos2

real (kind=8) :: zpli(NCSP)
real (kind=8) :: zpld(NCSP)

do jlat = 1 , NLPP

! set p(0,0) and p(0,1)

   zgwd    = gwd(jlat)            ! gaussian weight - from inigau
   zsin    = sid(jlat)            ! sin(phi) - from inigau
   zcsq    = 1.0_8 - zsin * zsin  ! cos(phi) squared
   zgwdcsq = zgwd / zcsq          ! weight / cos squared
   f1m     = sqrt(1.5_8)
   zpli(1) = sqrt(0.5_8)
   zpli(2) = f1m * zsin
   zpld(1) = 0.0
   lm      = 2

! loop over wavenumbers

   do m = 0 , NTRU
      if (m > 0) then
         lm  = lm + 1
         f2m = -f1m * sqrt(zcsq / (m+m))
         f1m =  f2m * sqrt(m+m + 3.0_8)
         zpli(lm) = f2m
         if (lm < NCSP) then
            lm = lm + 1
            zpli(lm  ) =       f1m * zsin
            zpld(lm-1) =  -m * f2m * zsin
         endif ! (lm < NCSP)
      endif ! (m > 0)

      amsq = m * m

      do n = m+2 , NTRU
         lm = lm + 1
         z1 = sqrt(((n-1)*(n-1) - amsq) / (4*(n-1)*(n-1)-1))
         z2 = zsin * zpli(lm-1) - z1 * zpli(lm-2)
         zpli(lm  ) = z2 * sqrt((4*n*n-1) / (n*n-amsq))
         zpld(lm-1) = (1-n) * z2 + n * z1 * zpli(lm-2)
      enddo ! n

      if (lm < NCSP) then ! mode (m,NTRU)
         z3 = sqrt((NTRU*NTRU-amsq) / (4*NTRU*NTRU-1))
         zpld(lm)=-NTRU*zsin*zpli(lm) + (NTRU+NTRU+1)*zpli(lm-1)*z3
      else                ! mode (NTRU,NTRU)
         zpld(lm)=-NTRU*zsin*zpli(lm)
      endif
   enddo ! m

   gwdc(jlat) = zgwdcsq

   lm = 0
   do m = 0 , NTRU
      do n = m , NTRU
           lm = lm + 1
           pmat(lm,jlat) = zpli(lm)
           qmat(lm,jlat) = zpld(lm)
      enddo ! n
   enddo ! m
   

   
   
enddo ! jlat

skgpsp(:) = 1.0
skspgp(:) = 1.0

   
do n=1,NTP1
   if (nfilter .eq. 0) then
     skgpsp(n) = 1.0    ! Filter for conversion from gp->sp
     skspgp(n) = 1.0    ! Filter for conversion from sp->gp
     
   else if (nfilter .eq. 1) then !Cesaro filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) +  ngptfilter*(1.0 - real(n)/(NTP1))
     skspgp(n) = (1-nspvfilter)*skspgp(n) +  nspvfilter*(1.0 - real(n)/(NTP1))
   
   else if (nfilter .eq. 2) then !Exponential filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + &
&                 ngptfilter*(exp(-filterkappa*(real(n)/NTRU)**nfilterexp))
     skspgp(n) = (1-nspvfilter)*skspgp(n) + &
&                 nspvfilter*(exp(-filterkappa*(real(n)/NTRU)**nfilterexp))
   
   else if (nfilter .eq. 3) then !Lander-Hoskins physics filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + ngptfilter*(exp(-(real(n)*real(n+1)/(landhoskn0*(landhoskn0+1)))**2))
     skspgp(n) = (1-nspvfilter)*skspgp(n) + nspvfilter*(exp(-(real(n)*real(n+1)/(landhoskn0*(landhoskn0+1)))**2))
   
   else if (nfilter .eq. 4) then !Riesz-2 filter
     skgpsp(n) = (1-ngptfilter)*skgpsp(n) + ngptfilter*((1.0 - real(n)/(NTP1))**2)
     skspgp(n) = (1-nspvfilter)*skspgp(n) + nspvfilter*((1.0 - real(n)/(NTP1))**2)
   
   else
     skgpsp(n) = 1.0
     skspgp(n) = 1.0
   endif
   
enddo

! ---------------------------------------------------------------------------
! Carry the physics filter into the per-mode factors.
!
! Every SP -> GP weight carries skspgp(n) and every GP -> SP weight carries
! skgpsp(n); no weight is ever used unfiltered. Both filters are functions of
! the total wavenumber alone, so they are constant for a spectral mode and
! belong here rather than in the innermost loop of every transform.
!
! It is now a per-MODE vector rather than a fold into the matrices, which is the
! same saving arrived at more cheaply: the filter is a function of n alone, so
! it belongs in something NCSP long that every latitude reads, not in something
! NCSP x NLPP that every latitude streams.
!
! This file counts m and n from 0 while the transform routines count from 1, so
! mode (m,n) here is skgpsp(n+1) there. The loop runs over lm because that is
! the fast index the transforms walk.
! ---------------------------------------------------------------------------

lm = 0
do m = 0 , NTRU
   do n = m , NTRU
      lm = lm + 1
      znn1 = 0.0
      if (n > 0) znn1 = 1.0_8 / (n*(n+1))
      fsp(lm) = skspgp(n+1)
      fgp(lm) = skgpsp(n+1)
      fmu(lm) = znn1 * m * skspgp(n+1)
      fmv(lm) = znn1 * skspgp(n+1)
      fmm(lm) = m * skgpsp(n+1)
      fmq(lm) = n * (n+1) * 0.5_8 * skgpsp(n+1)
   enddo ! n
enddo ! m

if (mypid==NROOT) then
   write(nud,*)"*************************************************"
   write(nud,*)"* GP->SP Physics Filter Coefficients:           *"
   do n=1,NTP1
      write(nud,*)"* n=",n,":  ",skgpsp(n)," *"
   enddo
   write(nud,*)"* SP->GP Physics Filter Coefficients:           *"
   do n=1,NTP1
      write(nud,*)"* n=",n,":  ",skspgp(n)," *"
   enddo
   write(nud,*)"*************************************************"
endif

return
end


! ================
! SUBROUTINE FC2SP
! ================

subroutine fc2sp(fc,sp)
use legmod
implicit none
real, intent(in ) :: fc(2,NLON/2,NLPP)
real, intent(out) :: sp(2,NESP/2)

integer :: l ! Index for latitude
integer :: m ! Index for zonal wavenumber
integer :: n ! Index for total wavenumber
integer :: w ! Index for spherical harmonic

sp(:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
  do l = 1 , NLPP
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        sp(1,w) = sp(1,w) + pmat(w,l)*fgp(w)*gwd(l) * fc(1,m,l)
        sp(2,w) = sp(2,w) + pmat(w,l)*fgp(w)*gwd(l) * fc(2,m,l)
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
  do l = 1 , NLHP
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        if (mod(m+n,2) == 0) then ! Symmetric modes
          sp(1,w) = sp(1,w) + pmat(w,l)*fgp(w)*gwd(l) * (fc(1,m,l) + fc(1,m,NLPP+1-l))
          sp(2,w) = sp(2,w) + pmat(w,l)*fgp(w)*gwd(l) * (fc(2,m,l) + fc(2,m,NLPP+1-l))
        else                      ! Antisymmetric modes
          sp(1,w) = sp(1,w) + pmat(w,l)*fgp(w)*gwd(l) * (fc(1,m,l) - fc(1,m,NLPP+1-l))
          sp(2,w) = sp(2,w) + pmat(w,l)*fgp(w)*gwd(l) * (fc(2,m,l) - fc(2,m,NLPP+1-l))
        endif
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
!----------------------------------------------------------------------
endif ! parallel ?
return
end


! ================
! SUBROUTINE SP2FC
! ================

subroutine sp2fc(sp,fc) ! Spectral to Fourier
use legmod
implicit none

real :: sp(2,NCSP)        ! Coefficients of spherical harmonics
real :: fc(2,NLON/2,NLPP) ! Fourier coefficients

integer :: j ! Loop index for spectral mode within one m
integer :: k ! Index for the mirror latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: w ! Index of the first spectral mode of one m

real :: ze1, ze2 ! partial sums over the symmetric modes
real :: zo1, zo2 ! partial sums over the antisymmetric modes

real :: zsp(2,NCSP) ! sp with the per-mode factor already in it

! The weight is pmat times a factor that depends only on the MODE, and the
! latitude loop is outside this, so the factor goes into the FIELD once rather
! than into the weight NLPP times. The inner loops below are then exactly the
! multiply they always were -- storing two matrices instead of eight costs
! nothing here at all.
do j = 1 , NCSP
   zsp(1,j) = fsp(j) * zsp(1,j)
   zsp(2,j) = fsp(j) * zsp(2,j)
enddo

fc(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do l = 1 , NLPP
   w = 1  
   do m = 1 , NTP1
      do n = m , NTP1
         fc(1,m,l) = fc(1,m,l) + pmat(w,l) * zsp(1,w)
         fc(2,m,l) = fc(2,m,l) + pmat(w,l) * zsp(2,w)
         w = w + 1
      enddo ! n
   enddo ! m
enddo ! l
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  qi carries P, so pmat(w,k) = s(w) * pmat(w,l) with s = (-1)**(m+n). Summing
!  the two parities apart therefore yields BOTH latitudes from one pass over
!  the modes: the mirror is the same two partial sums with the odd one
!  negated. Half the multiplies, and -- which matters more at this
!  resolution -- each element of qi is read once for the pair instead of
!  twice.
!
!  The parity is a property of the mode, so the n loop is split by STRIDE
!  rather than tested inside. A mod() in the innermost loop of the model's
!  hottest routine is a branch and a barrier to vectorisation both.
do l = 1 , NLHP
   k = NLPP + 1 - l
   w = 1
   do m = 1 , NTP1
      ze1 = 0.0
      ze2 = 0.0
      zo1 = 0.0
      zo2 = 0.0
      do j = w , w + NTP1 - m , 2       ! n = m, m+2, ...   symmetric
         ze1 = ze1 + pmat(j,l) * zsp(1,j)
         ze2 = ze2 + pmat(j,l) * zsp(2,j)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
         zo1 = zo1 + pmat(j,l) * zsp(1,j)
         zo2 = zo2 + pmat(j,l) * zsp(2,j)
      enddo ! j
      fc(1,m,l) = ze1 + zo1
      fc(2,m,l) = ze2 + zo2
      fc(1,m,k) = ze1 - zo1
      fc(2,m,k) = ze2 - zo2
      w = w + NTP1 - m + 1
   enddo ! m
enddo ! l
!----------------------------------------------------------------------
endif ! symmetric?
return
end


! ===================
! SUBROUTINE SP2FCDMU
! ===================

subroutine sp2fcdmu(sp,fc) ! Spectral to Fourier d/dmu
use legmod
implicit none

real :: sp(2,NCSP)        ! Coefficients of spherical harmonics
real :: fc(2,NLON/2,NLPP) ! Fourier coefficients

integer :: j ! Loop index for spectral mode within one m
integer :: k ! Index for the mirror latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: w ! Index of the first spectral mode of one m

real :: ze1, ze2 ! partial sums over the symmetric modes
real :: zo1, zo2 ! partial sums over the antisymmetric modes

real :: zsp(2,NCSP) ! sp with the per-mode factor already in it

! The weight is qmat times a factor that depends only on the MODE, and the
! latitude loop is outside this, so the factor goes into the FIELD once rather
! than into the weight NLPP times. The inner loops below are then exactly the
! multiply they always were -- storing two matrices instead of eight costs
! nothing here at all.
do j = 1 , NCSP
   zsp(1,j) = fsp(j) * zsp(1,j)
   zsp(2,j) = fsp(j) * zsp(2,j)
enddo

fc(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do l = 1 , NLPP
   w = 1  
   do m = 1 , NTP1
      do n = m , NTP1
         fc(1,m,l) = fc(1,m,l) + qmat(w,l) * zsp(1,w)
         fc(2,m,l) = fc(2,m,l) + qmat(w,l) * zsp(2,w)
         w = w + 1
      enddo ! n
   enddo ! m
enddo ! l
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  As sp2fc, with the OPPOSITE parity: qj carries Q = dP/dmu, and the
!  derivative of an even function is odd, so qmat(w,k) = -s(w) * qmat(w,l). The
!  mirror is therefore the odd partial sum MINUS the even one rather than
!  the other way round, and getting that sign the wrong way is the mistake
!  the single-mode check exists to catch.
do l = 1 , NLHP
   k = NLPP + 1 - l
   w = 1
   do m = 1 , NTP1
      ze1 = 0.0
      ze2 = 0.0
      zo1 = 0.0
      zo2 = 0.0
      do j = w , w + NTP1 - m , 2       ! n = m, m+2, ...   symmetric
         ze1 = ze1 + qmat(j,l) * zsp(1,j)
         ze2 = ze2 + qmat(j,l) * zsp(2,j)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
         zo1 = zo1 + qmat(j,l) * zsp(1,j)
         zo2 = zo2 + qmat(j,l) * zsp(2,j)
      enddo ! j
      fc(1,m,l) = ze1 + zo1
      fc(2,m,l) = ze2 + zo2
      fc(1,m,k) = zo1 - ze1
      fc(2,m,k) = zo2 - ze2
      w = w + NTP1 - m + 1
   enddo ! m
enddo ! l
!----------------------------------------------------------------------
endif ! symmetric?
return
end


! ================
! SUBROUTINE SP3FC
! ================

subroutine sp3fc
use pumamod, only:NLEV,sd,st,sz,gd,gt,gz
implicit none
integer :: v ! Loop index for level

do v = 1 , NLEV
   call sp2fc(sd(1,v),gd(1,v))
   call sp2fc(st(1,v),gt(1,v))
   call sp2fc(sz(1,v),gz(1,v))
enddo
return
end


! ================
! SUBROUTINE DV2UV        !SP->GP
! ================

subroutine dv2uv(pd,pz,pu,pv)
use legmod
implicit none

real :: pd(2,NESP/2,NLEV)
real :: pz(2,NESP/2,NLEV)
real :: pu(2,NLON/2,NLPP,NLEV)
real :: pv(2,NLON/2,NLPP,NLEV)

integer :: j ! Loop index for spectral mode within one m
integer :: k ! Index for the mirror latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: v ! Loop index for level
integer :: w ! Index of the first spectral mode of one m

! Eight products, each needing a symmetric and an antisymmetric partial sum.
! Named for what they carry: zv/zu the weight matrix, z/d the field, 1/2 the
! Fourier component, e/o the parity.
real :: zvz1e, zvz1o, zvz2e, zvz2o, zvd1e, zvd1o, zvd2e, zvd2o
real :: zuz1e, zuz1o, zuz2e, zuz2o, zud1e, zud1o, zud2e, zud2o

! The two weights are pmat and qmat times per-mode factors, and each weight
! multiplies BOTH fields, so the factors go into the fields once per level
! rather than into the weights NLPP times. Four combinations are needed and the
! inner loops are then the single multiply they always were.
real :: zzu(2,NCSP), zzv(2,NCSP) ! pz scaled by fmu and by fmv
real :: zdu(2,NCSP), zdv(2,NCSP) ! pd scaled by fmu and by fmv
integer :: jm                    ! mode index for the per-level scaling

pu(:,:,:,:) = 0.0
pv(:,:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do v = 1 , NLEV
  do jm = 1 , NCSP
     zzu(1,jm) = fmu(jm)*pz(1,jm,v) ; zzu(2,jm) = fmu(jm)*pz(2,jm,v)
     zzv(1,jm) = fmv(jm)*pz(1,jm,v) ; zzv(2,jm) = fmv(jm)*pz(2,jm,v)
     zdu(1,jm) = fmu(jm)*pd(1,jm,v) ; zdu(2,jm) = fmu(jm)*pd(2,jm,v)
     zdv(1,jm) = fmv(jm)*pd(1,jm,v) ; zdv(2,jm) = fmv(jm)*pd(2,jm,v)
  enddo
  do l = 1 , NLPP
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        pu(1,m,l,v)=pu(1,m,l,v)+qmat(w,l)*zzv(1,w)+pmat(w,l)*zdu(2,w)
        pu(2,m,l,v)=pu(2,m,l,v)+qmat(w,l)*zzv(2,w)-pmat(w,l)*zdu(1,w)
        pv(1,m,l,v)=pv(1,m,l,v)+pmat(w,l)*zzu(2,w)-qmat(w,l)*zdv(1,w)
        pv(2,m,l,v)=pv(2,m,l,v)-pmat(w,l)*zzu(1,w)-qmat(w,l)*zdv(2,w)
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
! Planetary vorticity rides on the single mode w=2, which is m=1,n=2, so its
! effect on the result is that one mode's contribution and nothing else. It is
! applied here rather than by subtracting it from pz in place: pz aliases the
! shared spectral state, every thread runs this routine over the whole of it,
! and a subtract-transform-restore around the loop is a write race on one
! element that costs run-to-run reproducibility. ThreadSanitizer reports it.
  do l = 1 , NLPP
    pu(1,1,l,v) = pu(1,1,l,v) - qmat(2,l) * fmv(2) * plavor
    pv(2,1,l,v) = pv(2,1,l,v) + pmat(2,l) * fmu(2) * plavor
  enddo ! l
enddo ! v
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  This one does NOT reduce to a single parity split, and pattern-matching
!  sp2fc onto it gives a wrong answer that looks right. dv2uv mixes pmat, carrying
!  P, with qmat, carrying dP/dmu, and the two have OPPOSITE parity:
!
!      pmat(w,k) =  s(w) * pmat(w,l)        s = (-1)**(m+n)
!      qmat(w,k) = -s(w) * qmat(w,l)
!
!  so u and v do not share a split and each of the eight products needs its
!  own pair of accumulators. Substituting the two relations into the four
!  outputs above gives the mirror as the same sixteen sums recombined, with
!  the qv terms and the odd terms each contributing a sign flip.
!
!  Sixteen live accumulators is the cost, and the risk: if they spill, the
!  saving goes to stack traffic instead. That is measured on the built code
!  rather than predicted.
do v = 1 , NLEV
  do jm = 1 , NCSP
     zzu(1,jm) = fmu(jm)*pz(1,jm,v) ; zzu(2,jm) = fmu(jm)*pz(2,jm,v)
     zzv(1,jm) = fmv(jm)*pz(1,jm,v) ; zzv(2,jm) = fmv(jm)*pz(2,jm,v)
     zdu(1,jm) = fmu(jm)*pd(1,jm,v) ; zdu(2,jm) = fmu(jm)*pd(2,jm,v)
     zdv(1,jm) = fmv(jm)*pd(1,jm,v) ; zdv(2,jm) = fmv(jm)*pd(2,jm,v)
  enddo
  do l = 1 , NLHP
    k = NLPP + 1 - l
    w = 1
    do m = 1 , NTP1
      zvz1e = 0.0
      zvz1o = 0.0
      zvz2e = 0.0
      zvz2o = 0.0
      zvd1e = 0.0
      zvd1o = 0.0
      zvd2e = 0.0
      zvd2o = 0.0
      zuz1e = 0.0
      zuz1o = 0.0
      zuz2e = 0.0
      zuz2o = 0.0
      zud1e = 0.0
      zud1o = 0.0
      zud2e = 0.0
      zud2o = 0.0
      do j = w , w + NTP1 - m , 2       ! n = m, m+2, ...   symmetric
        zvz1e = zvz1e + qmat(j,l)*zzv(1,j)
        zvz2e = zvz2e + qmat(j,l)*zzv(2,j)
        zvd1e = zvd1e + qmat(j,l)*zdv(1,j)
        zvd2e = zvd2e + qmat(j,l)*zdv(2,j)
        zuz1e = zuz1e + pmat(j,l)*zzu(1,j)
        zuz2e = zuz2e + pmat(j,l)*zzu(2,j)
        zud1e = zud1e + pmat(j,l)*zdu(1,j)
        zud2e = zud2e + pmat(j,l)*zdu(2,j)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
        zvz1o = zvz1o + qmat(j,l)*zzv(1,j)
        zvz2o = zvz2o + qmat(j,l)*zzv(2,j)
        zvd1o = zvd1o + qmat(j,l)*zdv(1,j)
        zvd2o = zvd2o + qmat(j,l)*zdv(2,j)
        zuz1o = zuz1o + pmat(j,l)*zzu(1,j)
        zuz2o = zuz2o + pmat(j,l)*zzu(2,j)
        zud1o = zud1o + pmat(j,l)*zdu(1,j)
        zud2o = zud2o + pmat(j,l)*zdu(2,j)
      enddo ! j
      pu(1,m,l,v) =  (zvz1e + zvz1o) + (zud2e + zud2o)
      pu(2,m,l,v) =  (zvz2e + zvz2o) - (zud1e + zud1o)
      pv(1,m,l,v) =  (zuz2e + zuz2o) - (zvd1e + zvd1o)
      pv(2,m,l,v) = -(zuz1e + zuz1o) - (zvd2e + zvd2o)
      pu(1,m,k,v) =  (zvz1o - zvz1e) + (zud2e - zud2o)
      pu(2,m,k,v) =  (zvz2o - zvz2e) + (zud1o - zud1e)
      pv(1,m,k,v) =  (zuz2e - zuz2o) + (zvd1e - zvd1o)
      pv(2,m,k,v) =  (zuz1o - zuz1e) + (zvd2e - zvd2o)
      w = w + NTP1 - m + 1
    enddo ! m
  enddo ! l
! The same planetary vorticity term as the branch above, and here it is not
! symmetric with itself: w=2 is n=2 against m=1, so it falls in the ODD loop
! and reaches zvz1o and zuz1o alone. Reading those two out of the four output
! recombinations gives the mirror latitude a sign the home latitude does not
! have -- taking the home form for both is the mistake this shape invites.
  do l = 1 , NLHP
    k = NLPP + 1 - l
    pu(1,1,l,v) = pu(1,1,l,v) - qmat(2,l) * fmv(2) * plavor
    pv(2,1,l,v) = pv(2,1,l,v) + pmat(2,l) * fmu(2) * plavor
    pu(1,1,k,v) = pu(1,1,k,v) - qmat(2,l) * fmv(2) * plavor
    pv(2,1,k,v) = pv(2,1,k,v) - pmat(2,l) * fmu(2) * plavor
  enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end


! ================
! SUBROUTINE UV2DV        !GP->SP
! ================

subroutine uv2dv(pu,pv,pd,pz)
use legmod
implicit none

real :: pd(2,NESP/2,NLEV)
real :: pz(2,NESP/2,NLEV)
real :: pu(2,NLON/2,NLPP,NLEV)
real :: pv(2,NLON/2,NLPP,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: v ! Loop index for level
integer :: w ! Loop index for spectral mode

pd(:,:,:) = 0.0
pz(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do v = 1 , NLEV
  do l = 1 , NLPP
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        pz(1,w,v) = pz(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*pu(1,m,l,v)-pmat(w,l)*fmm(w)*gwdc(l)*pv(2,m,l,v)
        pz(2,w,v) = pz(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*pu(2,m,l,v)+pmat(w,l)*fmm(w)*gwdc(l)*pv(1,m,l,v)
        pd(1,w,v) = pd(1,w,v)-qmat(w,l)*fgp(w)*gwdc(l)*pv(1,m,l,v)-pmat(w,l)*fmm(w)*gwdc(l)*pu(2,m,l,v)
        pd(2,w,v) = pd(2,w,v)-qmat(w,l)*fgp(w)*gwdc(l)*pv(2,m,l,v)+pmat(w,l)*fmm(w)*gwdc(l)*pu(1,m,l,v)
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
enddo ! v
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
do v = 1 , NLEV
  do l = 1 , NLHP
    k = NLPP+1-l
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        if (mod(m+n,2) == 0) then ! symmetric -----------------
          pz(1,w,v) = pz(1,w,v) + qmat(w,l)*fgp(w)*gwdc(l) * (pu(1,m,l,v)-pu(1,m,k,v)) &
                                - pmat(w,l)*fmm(w)*gwdc(l) * (pv(2,m,l,v)+pv(2,m,k,v))
          pz(2,w,v) = pz(2,w,v) + qmat(w,l)*fgp(w)*gwdc(l) * (pu(2,m,l,v)-pu(2,m,k,v)) &
                                + pmat(w,l)*fmm(w)*gwdc(l) * (pv(1,m,l,v)+pv(1,m,k,v))
          pd(1,w,v) = pd(1,w,v) - qmat(w,l)*fgp(w)*gwdc(l) * (pv(1,m,l,v)-pv(1,m,k,v)) &
                                - pmat(w,l)*fmm(w)*gwdc(l) * (pu(2,m,l,v)+pu(2,m,k,v))
          pd(2,w,v) = pd(2,w,v) - qmat(w,l)*fgp(w)*gwdc(l) * (pv(2,m,l,v)-pv(2,m,k,v)) &
                                + pmat(w,l)*fmm(w)*gwdc(l) * (pu(1,m,l,v)+pu(1,m,k,v))
        else ! ---------------- antisymmetric -----------------
          pz(1,w,v) = pz(1,w,v) + qmat(w,l)*fgp(w)*gwdc(l) * (pu(1,m,l,v)+pu(1,m,k,v)) &
                                - pmat(w,l)*fmm(w)*gwdc(l) * (pv(2,m,l,v)-pv(2,m,k,v))
          pz(2,w,v) = pz(2,w,v) + qmat(w,l)*fgp(w)*gwdc(l) * (pu(2,m,l,v)+pu(2,m,k,v)) &
                                + pmat(w,l)*fmm(w)*gwdc(l) * (pv(1,m,l,v)-pv(1,m,k,v))
          pd(1,w,v) = pd(1,w,v) - qmat(w,l)*fgp(w)*gwdc(l) * (pv(1,m,l,v)+pv(1,m,k,v)) &
                                - pmat(w,l)*fmm(w)*gwdc(l) * (pu(2,m,l,v)-pu(2,m,k,v))
          pd(2,w,v) = pd(2,w,v) - qmat(w,l)*fgp(w)*gwdc(l) * (pv(2,m,l,v)+pv(2,m,k,v)) &
                                + pmat(w,l)*fmm(w)*gwdc(l) * (pu(1,m,l,v)-pu(1,m,k,v))
        endif
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end 


! ================
! SUBROUTINE QTEND
! ================

subroutine qtend(q,qn,uq,vq)
use legmod
implicit none

real, intent(in) :: qn(2,NLON/2,NLPP,NLEV)
real, intent(in) :: uq(2,NLON/2,NLPP,NLEV)
real, intent(in) :: vq(2,NLON/2,NLPP,NLEV)

real, intent(out) :: q(2,NESP/2,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: v ! Loop index for level
integer :: w ! Loop index for spectral mode

q(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLPP
  w = 1
  do m = 1 , NTP1
   do n = m , NTP1
    q(1,w,v)=q(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*vq(1,m,l,v)+&
&            pmat(w,l)*fgp(w)*gwd(l)*qn(1,m,l,v)+pmat(w,l)*fmm(w)*gwdc(l)*uq(2,m,l,v)
    q(2,w,v)=q(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*vq(2,m,l,v)+&
&            pmat(w,l)*fgp(w)*gwd(l)*qn(2,m,l,v)-pmat(w,l)*fmm(w)*gwdc(l)*uq(1,m,l,v)
    w = w + 1
   enddo ! n
  enddo ! m
 enddo ! l
enddo ! v
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLHP
  k = NLPP+1-l
  w = 1
  do m = 1 , NTP1
   do n = m , NTP1
    if (mod(m+n,2) == 0) then ! symmetric -----------------
      q(1,w,v)=q(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vq(1,m,l,v)-vq(1,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(qn(1,m,l,v)+qn(1,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(uq(2,m,l,v)+uq(2,m,k,v))
      q(2,w,v)=q(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vq(2,m,l,v)-vq(2,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(qn(2,m,l,v)+qn(2,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(uq(1,m,l,v)+uq(1,m,k,v))
    else ! ---------------- antisymmetric -----------------
      q(1,w,v)=q(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vq(1,m,l,v)+vq(1,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(qn(1,m,l,v)-qn(1,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(uq(2,m,l,v)-uq(2,m,k,v))
      q(2,w,v)=q(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vq(2,m,l,v)+vq(2,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(qn(2,m,l,v)-qn(2,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(uq(1,m,l,v)-uq(1,m,k,v))
    endif
    w = w + 1
   enddo ! n
  enddo ! m
 enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end

! =================
! SUBROUTINE MKTEND
! =================

subroutine mktend(d,t,z,tn,fu,fv,ke,ut,vt)
use legmod
implicit none

real, intent(in) :: tn(2,NLON/2,NLPP,NLEV)
real, intent(in) :: fu(2,NLON/2,NLPP,NLEV)
real, intent(in) :: fv(2,NLON/2,NLPP,NLEV)
real, intent(in) :: ke(2,NLON/2,NLPP,NLEV)
real, intent(in) :: ut(2,NLON/2,NLPP,NLEV)
real, intent(in) :: vt(2,NLON/2,NLPP,NLEV)

real, intent(out) :: d(2,NESP/2,NLEV)
real, intent(out) :: t(2,NESP/2,NLEV)
real, intent(out) :: z(2,NESP/2,NLEV)

integer :: k ! Loop index for southern latitude
integer :: l ! Loop index for latitude
integer :: m ! Loop index for zonal wavenumber m
integer :: n ! Loop index for total wavenumber n
integer :: v ! Loop index for level
integer :: w ! Loop index for spectral mode

d(:,:,:) = 0.0
t(:,:,:) = 0.0
z(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLPP
  w = 1
  do m = 1 , NTP1
   do n = m , NTP1
    d(1,w,v)=d(1,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*ke(1,m,l,v) &
&                    -qmat(w,l)*fgp(w)*gwdc(l)*fv(1,m,l,v)-pmat(w,l)*fmm(w)*gwdc(l)*fu(2,m,l,v)
    d(2,w,v)=d(2,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*ke(2,m,l,v) &
&                    -qmat(w,l)*fgp(w)*gwdc(l)*fv(2,m,l,v)+pmat(w,l)*fmm(w)*gwdc(l)*fu(1,m,l,v)
    t(1,w,v)=t(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*vt(1,m,l,v) &
&                    +pmat(w,l)*fgp(w)*gwd(l)*tn(1,m,l,v)+pmat(w,l)*fmm(w)*gwdc(l)*ut(2,m,l,v)
    t(2,w,v)=t(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*vt(2,m,l,v) &
&                    +pmat(w,l)*fgp(w)*gwd(l)*tn(2,m,l,v)-pmat(w,l)*fmm(w)*gwdc(l)*ut(1,m,l,v)
    z(1,w,v)=z(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*fu(1,m,l,v) &
&                    -pmat(w,l)*fmm(w)*gwdc(l)*fv(2,m,l,v)
    z(2,w,v)=z(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*fu(2,m,l,v) &
&                    +pmat(w,l)*fmm(w)*gwdc(l)*fv(1,m,l,v)
    w = w + 1
   enddo ! n
  enddo ! m
 enddo ! l
enddo ! v
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
do v = 1 , NLEV
 do l = 1 , NLHP
  k = NLPP+1-l
  w = 1
  do m = 1 , NTP1
   do n = m , NTP1
    if (mod(m+n,2) == 0) then ! symmetric -----------------
      d(1,w,v)=d(1,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*(ke(1,m,l,v)+ke(1,m,k,v)) &
                       -qmat(w,l)*fgp(w)*gwdc(l)*(fv(1,m,l,v)-fv(1,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(fu(2,m,l,v)+fu(2,m,k,v))
      d(2,w,v)=d(2,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*(ke(2,m,l,v)+ke(2,m,k,v)) &
                       -qmat(w,l)*fgp(w)*gwdc(l)*(fv(2,m,l,v)-fv(2,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(fu(1,m,l,v)+fu(1,m,k,v))
      t(1,w,v)=t(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vt(1,m,l,v)-vt(1,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(tn(1,m,l,v)+tn(1,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(ut(2,m,l,v)+ut(2,m,k,v))
      t(2,w,v)=t(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vt(2,m,l,v)-vt(2,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(tn(2,m,l,v)+tn(2,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(ut(1,m,l,v)+ut(1,m,k,v))
      z(1,w,v)=z(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(fu(1,m,l,v)-fu(1,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(fv(2,m,l,v)+fv(2,m,k,v))
      z(2,w,v)=z(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(fu(2,m,l,v)-fu(2,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(fv(1,m,l,v)+fv(1,m,k,v))
    else ! ---------------- antisymmetric -----------------
      d(1,w,v)=d(1,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*(ke(1,m,l,v)-ke(1,m,k,v)) &
                       -qmat(w,l)*fgp(w)*gwdc(l)*(fv(1,m,l,v)+fv(1,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(fu(2,m,l,v)-fu(2,m,k,v))
      d(2,w,v)=d(2,w,v)+pmat(w,l)*fmq(w)*gwdc(l)*(ke(2,m,l,v)-ke(2,m,k,v)) &
                       -qmat(w,l)*fgp(w)*gwdc(l)*(fv(2,m,l,v)+fv(2,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(fu(1,m,l,v)-fu(1,m,k,v))
      t(1,w,v)=t(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vt(1,m,l,v)+vt(1,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(tn(1,m,l,v)-tn(1,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(ut(2,m,l,v)-ut(2,m,k,v))
      t(2,w,v)=t(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(vt(2,m,l,v)+vt(2,m,k,v)) &
                       +pmat(w,l)*fgp(w)*gwd(l)*(tn(2,m,l,v)-tn(2,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(ut(1,m,l,v)-ut(1,m,k,v))
      z(1,w,v)=z(1,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(fu(1,m,l,v)+fu(1,m,k,v)) &
                       -pmat(w,l)*fmm(w)*gwdc(l)*(fv(2,m,l,v)-fv(2,m,k,v))
      z(2,w,v)=z(2,w,v)+qmat(w,l)*fgp(w)*gwdc(l)*(fu(2,m,l,v)+fu(2,m,k,v)) &
                       +pmat(w,l)*fmm(w)*gwdc(l)*(fv(1,m,l,v)-fv(1,m,k,v))
    endif
    w = w + 1
   enddo ! n
  enddo ! m
 enddo ! l
enddo ! v
!----------------------------------------------------------------------
endif ! symmetric?
return
end


! ================
! SUBROUTINE SP2FL
! ================

subroutine sp2fl(psp,pfc,klev)
use legmod
implicit none

integer, intent(in ) :: klev
real,    intent(in ) :: psp(NESP,klev)
real,    intent(out) :: pfc(NHOR,klev)

integer :: jlev

do jlev = 1,klev
   call sp2fc(psp(1,jlev),pfc(1,jlev))
enddo

return
end 


! ==================
! SUBROUTINE INVLEGA
! ==================

subroutine invlega
use pumamod
implicit none

integer :: jlev

call dv2uv(sd,sz,gu,gv)

do jlev = 1,NLEV
  call sp2fc(sd(1,jlev),gd(1,jlev))
  call sp2fc(st(1,jlev),gt(1,jlev))
  call sp2fc(sz(1,jlev),gz(1,jlev))
  if (nqspec == 1) call sp2fc(sq(1,jlev),gq(1,jlev))
enddo

call sp2fc(sp,gp)
call sp2fcdmu(sp,gpj)

return
end


! ==================
! SUBROUTINE INVLEGD
! ==================

subroutine invlegd
use pumamod
implicit none

integer :: jlev

call dv2uv(sd,sz,gu,gv)

do jlev = 1,NLEV
  if (nqspec == 1) call sp2fc(sq(1,jlev),gq(1,jlev))
  call sp2fc(st(1,jlev),gt(1,jlev))
enddo

call sp2fc(sp,gp)

return
end

