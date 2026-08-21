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

real :: qi(NCSP,NLPP) ! P(m,n) * skspgp                     used in sp2fc
real :: qj(NCSP,NLPP) ! Q(m,n) * skspgp                     used in sp2fcdmu
real :: qc(NCSP,NLPP) ! P(m,n) * gwd * skgpsp               used in fc2sp, qtend, mktend
real :: qe(NCSP,NLPP) ! Q(m,n) * gwd / cos2 * skgpsp        used in mktend, qtend, uv2dv
real :: qm(NCSP,NLPP) ! P(m,n) * gwd / cos2 * m * skgpsp    used in mktend, qtend, uv2dv
real :: qq(NCSP,NLPP) ! P(m,n) * gwd / cos2 * n*(n+1)/2 * skgpsp   used in mktend
real :: qu(NCSP,NLPP) ! P(m,n) / (n*(n+1)) * m * skspgp     used in dv2uv
real :: qv(NCSP,NLPP) ! Q(m,n) / (n*(n+1)) * skspgp         used in dv2uv
! Filter folded in by legini: skgpsp into qc/qe/qm/qq, skspgp into
! qi/qj/qu/qv. Kept for reference and for the diagnostic print below; the
! transform loops no longer read them.
real :: skgpsp(NTP1) ! Physics filter for GP -> SP
real :: skspgp(NTP1) ! Physics filter for SP -> GP

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

   lm = 0
   do m = 0 , NTRU
      do n = m , NTRU
           lm = lm + 1
           znn1 = 0.0
           if (n > 0) znn1 = 1.0_8 / (n*(n+1))
           qi(lm,jlat) = zpli(lm)
           qj(lm,jlat) = zpld(lm)
           qc(lm,jlat) = zpli(lm) * zgwd
           qu(lm,jlat) = zpli(lm) * znn1 * m
           qv(lm,jlat) = zpld(lm) * znn1
           qe(lm,jlat) = zpld(lm) * zgwdcsq
           qq(lm,jlat) = zpli(lm) * zgwdcsq * n * (n+1) * 0.5_8
           qm(lm,jlat) = zpli(lm) * zgwdcsq * m
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
! Fold the physics filter into the weight matrices.
!
! Every use of qi, qj, qu and qv carries skspgp(n), and every use of qc, qe, qm
! and qq carries skgpsp(n); no weight matrix is ever used unfiltered. Both
! filters are functions of the total wavenumber alone, so they are constant for
! a spectral mode and belong here rather than in the innermost loop of every
! transform, where they cost one multiply in three.
!
! gfortran cannot do this itself. The filter is invariant in the LATITUDE loop,
! but hoisting it needs a temporary array, which the compiler will not invent.
!
! This file counts m and n from 0 while the transform routines count from 1, so
! mode (m,n) here is skgpsp(n+1) there. The inner loop runs over lm because that
! is the fast index of q(NCSP,NLPP).
! ---------------------------------------------------------------------------

do jlat = 1 , NLPP
   lm = 0
   do m = 0 , NTRU
      do n = m , NTRU
         lm = lm + 1
         qi(lm,jlat) = qi(lm,jlat) * skspgp(n+1)
         qj(lm,jlat) = qj(lm,jlat) * skspgp(n+1)
         qu(lm,jlat) = qu(lm,jlat) * skspgp(n+1)
         qv(lm,jlat) = qv(lm,jlat) * skspgp(n+1)
         qc(lm,jlat) = qc(lm,jlat) * skgpsp(n+1)
         qe(lm,jlat) = qe(lm,jlat) * skgpsp(n+1)
         qm(lm,jlat) = qm(lm,jlat) * skgpsp(n+1)
         qq(lm,jlat) = qq(lm,jlat) * skgpsp(n+1)
      enddo ! n
   enddo ! m
enddo ! jlat

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
        sp(1,w) = sp(1,w) + qc(w,l) * fc(1,m,l)
        sp(2,w) = sp(2,w) + qc(w,l) * fc(2,m,l)
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
          sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) + fc(1,m,NLPP+1-l))
          sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) + fc(2,m,NLPP+1-l))
        else                      ! Antisymmetric modes
          sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) - fc(1,m,NLPP+1-l))
          sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) - fc(2,m,NLPP+1-l))
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

fc(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do l = 1 , NLPP
   w = 1  
   do m = 1 , NTP1
      do n = m , NTP1
         fc(1,m,l) = fc(1,m,l) + qi(w,l) * sp(1,w)
         fc(2,m,l) = fc(2,m,l) + qi(w,l) * sp(2,w)
         w = w + 1
      enddo ! n
   enddo ! m
enddo ! l
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  qi carries P, so qi(w,k) = s(w) * qi(w,l) with s = (-1)**(m+n). Summing
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
         ze1 = ze1 + qi(j,l) * sp(1,j)
         ze2 = ze2 + qi(j,l) * sp(2,j)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
         zo1 = zo1 + qi(j,l) * sp(1,j)
         zo2 = zo2 + qi(j,l) * sp(2,j)
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

fc(:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do l = 1 , NLPP
   w = 1  
   do m = 1 , NTP1
      do n = m , NTP1
         fc(1,m,l) = fc(1,m,l) + qj(w,l) * sp(1,w)
         fc(2,m,l) = fc(2,m,l) + qj(w,l) * sp(2,w)
         w = w + 1
      enddo ! n
   enddo ! m
enddo ! l
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  As sp2fc, with the OPPOSITE parity: qj carries Q = dP/dmu, and the
!  derivative of an even function is odd, so qj(w,k) = -s(w) * qj(w,l). The
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
         ze1 = ze1 + qj(j,l) * sp(1,j)
         ze2 = ze2 + qj(j,l) * sp(2,j)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
         zo1 = zo1 + qj(j,l) * sp(1,j)
         zo2 = zo2 + qj(j,l) * sp(2,j)
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
real :: zsave

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

pu(:,:,:,:) = 0.0
pv(:,:,:,:) = 0.0

if (.not. LPAIRLAT) then ! Contiguous latitudes: no mirror is local
!----------------------------------------------------------------------
do v = 1 , NLEV
  zsave = pz(1,2,v)
  pz(1,2,v) = zsave - plavor
  do l = 1 , NLPP
    w = 1
    do m = 1 , NTP1
      do n = m , NTP1
        pu(1,m,l,v)=pu(1,m,l,v)+qv(w,l)*pz(1,w,v)+qu(w,l)*pd(2,w,v)
        pu(2,m,l,v)=pu(2,m,l,v)+qv(w,l)*pz(2,w,v)-qu(w,l)*pd(1,w,v)
        pv(1,m,l,v)=pv(1,m,l,v)+qu(w,l)*pz(2,w,v)-qv(w,l)*pd(1,w,v)
        pv(2,m,l,v)=pv(2,m,l,v)-qu(w,l)*pz(1,w,v)-qv(w,l)*pd(2,w,v)
        w = w + 1
      enddo ! n
    enddo ! m
  enddo ! l
  pz(1,2,v) = zsave
enddo ! v
else                     ! Paired latitudes: symmetry conserving
!----------------------------------------------------------------------
!  This one does NOT reduce to a single parity split, and pattern-matching
!  sp2fc onto it gives a wrong answer that looks right. dv2uv mixes qu, from
!  P, with qv, from dP/dmu, and the two have OPPOSITE parity:
!
!      qu(w,k) =  s(w) * qu(w,l)        s = (-1)**(m+n)
!      qv(w,k) = -s(w) * qv(w,l)
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
  zsave = pz(1,2,v)
  pz(1,2,v) = zsave - plavor
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
        zvz1e = zvz1e + qv(j,l)*pz(1,j,v)
        zvz2e = zvz2e + qv(j,l)*pz(2,j,v)
        zvd1e = zvd1e + qv(j,l)*pd(1,j,v)
        zvd2e = zvd2e + qv(j,l)*pd(2,j,v)
        zuz1e = zuz1e + qu(j,l)*pz(1,j,v)
        zuz2e = zuz2e + qu(j,l)*pz(2,j,v)
        zud1e = zud1e + qu(j,l)*pd(1,j,v)
        zud2e = zud2e + qu(j,l)*pd(2,j,v)
      enddo ! j
      do j = w + 1 , w + NTP1 - m , 2   ! n = m+1, m+3, ... antisymmetric
        zvz1o = zvz1o + qv(j,l)*pz(1,j,v)
        zvz2o = zvz2o + qv(j,l)*pz(2,j,v)
        zvd1o = zvd1o + qv(j,l)*pd(1,j,v)
        zvd2o = zvd2o + qv(j,l)*pd(2,j,v)
        zuz1o = zuz1o + qu(j,l)*pz(1,j,v)
        zuz2o = zuz2o + qu(j,l)*pz(2,j,v)
        zud1o = zud1o + qu(j,l)*pd(1,j,v)
        zud2o = zud2o + qu(j,l)*pd(2,j,v)
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
  pz(1,2,v) = zsave
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
        pz(1,w,v) = pz(1,w,v)+qe(w,l)*pu(1,m,l,v)-qm(w,l)*pv(2,m,l,v)
        pz(2,w,v) = pz(2,w,v)+qe(w,l)*pu(2,m,l,v)+qm(w,l)*pv(1,m,l,v)
        pd(1,w,v) = pd(1,w,v)-qe(w,l)*pv(1,m,l,v)-qm(w,l)*pu(2,m,l,v)
        pd(2,w,v) = pd(2,w,v)-qe(w,l)*pv(2,m,l,v)+qm(w,l)*pu(1,m,l,v)
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
          pz(1,w,v) = pz(1,w,v) + qe(w,l) * (pu(1,m,l,v)-pu(1,m,k,v)) &
                                - qm(w,l) * (pv(2,m,l,v)+pv(2,m,k,v))
          pz(2,w,v) = pz(2,w,v) + qe(w,l) * (pu(2,m,l,v)-pu(2,m,k,v)) &
                                + qm(w,l) * (pv(1,m,l,v)+pv(1,m,k,v))
          pd(1,w,v) = pd(1,w,v) - qe(w,l) * (pv(1,m,l,v)-pv(1,m,k,v)) &
                                - qm(w,l) * (pu(2,m,l,v)+pu(2,m,k,v))
          pd(2,w,v) = pd(2,w,v) - qe(w,l) * (pv(2,m,l,v)-pv(2,m,k,v)) &
                                + qm(w,l) * (pu(1,m,l,v)+pu(1,m,k,v))
        else ! ---------------- antisymmetric -----------------
          pz(1,w,v) = pz(1,w,v) + qe(w,l) * (pu(1,m,l,v)+pu(1,m,k,v)) &
                                - qm(w,l) * (pv(2,m,l,v)-pv(2,m,k,v))
          pz(2,w,v) = pz(2,w,v) + qe(w,l) * (pu(2,m,l,v)+pu(2,m,k,v)) &
                                + qm(w,l) * (pv(1,m,l,v)-pv(1,m,k,v))
          pd(1,w,v) = pd(1,w,v) - qe(w,l) * (pv(1,m,l,v)+pv(1,m,k,v)) &
                                - qm(w,l) * (pu(2,m,l,v)-pu(2,m,k,v))
          pd(2,w,v) = pd(2,w,v) - qe(w,l) * (pv(2,m,l,v)+pv(2,m,k,v)) &
                                + qm(w,l) * (pu(1,m,l,v)-pu(1,m,k,v))
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
    q(1,w,v)=q(1,w,v)+qe(w,l)*vq(1,m,l,v)+&
&            qc(w,l)*qn(1,m,l,v)+qm(w,l)*uq(2,m,l,v)
    q(2,w,v)=q(2,w,v)+qe(w,l)*vq(2,m,l,v)+&
&            qc(w,l)*qn(2,m,l,v)-qm(w,l)*uq(1,m,l,v)
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
      q(1,w,v)=q(1,w,v)+qe(w,l)*(vq(1,m,l,v)-vq(1,m,k,v)) &
                       +qc(w,l)*(qn(1,m,l,v)+qn(1,m,k,v)) &
                       +qm(w,l)*(uq(2,m,l,v)+uq(2,m,k,v))
      q(2,w,v)=q(2,w,v)+qe(w,l)*(vq(2,m,l,v)-vq(2,m,k,v)) &
                       +qc(w,l)*(qn(2,m,l,v)+qn(2,m,k,v)) &
                       -qm(w,l)*(uq(1,m,l,v)+uq(1,m,k,v))
    else ! ---------------- antisymmetric -----------------
      q(1,w,v)=q(1,w,v)+qe(w,l)*(vq(1,m,l,v)+vq(1,m,k,v)) &
                       +qc(w,l)*(qn(1,m,l,v)-qn(1,m,k,v)) &
                       +qm(w,l)*(uq(2,m,l,v)-uq(2,m,k,v))
      q(2,w,v)=q(2,w,v)+qe(w,l)*(vq(2,m,l,v)+vq(2,m,k,v)) &
                       +qc(w,l)*(qn(2,m,l,v)-qn(2,m,k,v)) &
                       -qm(w,l)*(uq(1,m,l,v)-uq(1,m,k,v))
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
    d(1,w,v)=d(1,w,v)+qq(w,l)*ke(1,m,l,v) &
&                    -qe(w,l)*fv(1,m,l,v)-qm(w,l)*fu(2,m,l,v)
    d(2,w,v)=d(2,w,v)+qq(w,l)*ke(2,m,l,v) &
&                    -qe(w,l)*fv(2,m,l,v)+qm(w,l)*fu(1,m,l,v)
    t(1,w,v)=t(1,w,v)+qe(w,l)*vt(1,m,l,v) &
&                    +qc(w,l)*tn(1,m,l,v)+qm(w,l)*ut(2,m,l,v)
    t(2,w,v)=t(2,w,v)+qe(w,l)*vt(2,m,l,v) &
&                    +qc(w,l)*tn(2,m,l,v)-qm(w,l)*ut(1,m,l,v)
    z(1,w,v)=z(1,w,v)+qe(w,l)*fu(1,m,l,v) &
&                    -qm(w,l)*fv(2,m,l,v)
    z(2,w,v)=z(2,w,v)+qe(w,l)*fu(2,m,l,v) &
&                    +qm(w,l)*fv(1,m,l,v)
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
      d(1,w,v)=d(1,w,v)+qq(w,l)*(ke(1,m,l,v)+ke(1,m,k,v)) &
                       -qe(w,l)*(fv(1,m,l,v)-fv(1,m,k,v)) &
                       -qm(w,l)*(fu(2,m,l,v)+fu(2,m,k,v))
      d(2,w,v)=d(2,w,v)+qq(w,l)*(ke(2,m,l,v)+ke(2,m,k,v)) &
                       -qe(w,l)*(fv(2,m,l,v)-fv(2,m,k,v)) &
                       +qm(w,l)*(fu(1,m,l,v)+fu(1,m,k,v))
      t(1,w,v)=t(1,w,v)+qe(w,l)*(vt(1,m,l,v)-vt(1,m,k,v)) &
                       +qc(w,l)*(tn(1,m,l,v)+tn(1,m,k,v)) &
                       +qm(w,l)*(ut(2,m,l,v)+ut(2,m,k,v))
      t(2,w,v)=t(2,w,v)+qe(w,l)*(vt(2,m,l,v)-vt(2,m,k,v)) &
                       +qc(w,l)*(tn(2,m,l,v)+tn(2,m,k,v)) &
                       -qm(w,l)*(ut(1,m,l,v)+ut(1,m,k,v))
      z(1,w,v)=z(1,w,v)+qe(w,l)*(fu(1,m,l,v)-fu(1,m,k,v)) &
                       -qm(w,l)*(fv(2,m,l,v)+fv(2,m,k,v))
      z(2,w,v)=z(2,w,v)+qe(w,l)*(fu(2,m,l,v)-fu(2,m,k,v)) &
                       +qm(w,l)*(fv(1,m,l,v)+fv(1,m,k,v))
    else ! ---------------- antisymmetric -----------------
      d(1,w,v)=d(1,w,v)+qq(w,l)*(ke(1,m,l,v)-ke(1,m,k,v)) &
                       -qe(w,l)*(fv(1,m,l,v)+fv(1,m,k,v)) &
                       -qm(w,l)*(fu(2,m,l,v)-fu(2,m,k,v))
      d(2,w,v)=d(2,w,v)+qq(w,l)*(ke(2,m,l,v)-ke(2,m,k,v)) &
                       -qe(w,l)*(fv(2,m,l,v)+fv(2,m,k,v)) &
                       +qm(w,l)*(fu(1,m,l,v)-fu(1,m,k,v))
      t(1,w,v)=t(1,w,v)+qe(w,l)*(vt(1,m,l,v)+vt(1,m,k,v)) &
                       +qc(w,l)*(tn(1,m,l,v)-tn(1,m,k,v)) &
                       +qm(w,l)*(ut(2,m,l,v)-ut(2,m,k,v))
      t(2,w,v)=t(2,w,v)+qe(w,l)*(vt(2,m,l,v)+vt(2,m,k,v)) &
                       +qc(w,l)*(tn(2,m,l,v)-tn(2,m,k,v)) &
                       -qm(w,l)*(ut(1,m,l,v)-ut(1,m,k,v))
      z(1,w,v)=z(1,w,v)+qe(w,l)*(fu(1,m,l,v)+fu(1,m,k,v)) &
                       -qm(w,l)*(fv(2,m,l,v)-fv(2,m,k,v))
      z(2,w,v)=z(2,w,v)+qe(w,l)*(fu(2,m,l,v)+fu(2,m,k,v)) &
                       +qm(w,l)*(fv(1,m,l,v)-fv(1,m,k,v))
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

