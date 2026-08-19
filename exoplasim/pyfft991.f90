
      
      ! =================
      ! SUBROUTINE LEGINI
      ! =================
      
      subroutine legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv,&
    &                   sfilt,nfs)
      implicit none
      
      integer :: jlat ! Latitude
      integer :: lm
      integer :: m
      integer :: n
      
      integer :: nfs
      
      integer :: NLAT, NLON, NTRU, NLEV
      
      real (kind=8) :: EZ
      real (kind=8) :: PI, TWOPI
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
      
      real (kind=8) :: zpli((NTRU+1)*(NTRU+2)/2)
      real (kind=8) :: zpld((NTRU+1)*(NTRU+2)/2)
      
      real (kind=8) :: sid(NLAT)     ! sin(phi)
      real (kind=8) :: gwd(NLAT)     ! Gaussian weights
      real (kind=8) :: csq(NLAT)     ! cos(phi)**2
      real (kind=8) :: cola(NLAT)    ! cos(phi)
      real (kind=8) :: rcs(NLAT)     ! 1 / cos(phi)
      real (kind=8) :: deglat(NLat)  ! latitude in degrees
            
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
          
      real (kind=8):: sfilt(NTRU+1) ! Physics filter
      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
          
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call inigau(NLAT,sid,gwd)
      
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
               zpld(lm)=-NTRU*zsin*zpli(lm) +  &
     &                  (NTRU+NTRU+1)*zpli(lm-1)*z3
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
         
      enddo! jlat
      
      sfilt(:) = 1.0
         
      do n=1,NTP1
         sfilt(n) = (1-nfs)*sfilt(n)+nfs*exp(-8*(real(n)/NTRU)**8)
      enddo
      
      return
      end
            
              
      ! ================
      ! SUBROUTINE FC2SP
      ! ================
      
      subroutine fc2sp(fc,sp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: nfs
      
      real (kind=8), intent(in ) :: fc(2,NLON/2,NLAT)
!f2py intent(in) :: fc
      real (kind=8), intent(out) :: sp(2,(NTRU+1)*(NTRU+2)/2)
!f2py intent(out) :: sp
      
      integer :: l ! Index for latitude
      integer :: m ! Index for zonal wavenumber
      integer :: n ! Index for total wavenumber
      integer :: w ! Index for spherical harmonic
       
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
           
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv,&
     &            sfilt,nfs)
      
           
      sp(:,:) = 0.0
      
      if (NLPP < NLAT) then  ! Universal (parallel executable) version
      !----------------------------------------------------------------------
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              sp(1,w) = sp(1,w) + qc(w,l) * fc(1,m,l)*sfilt(n)
              sp(2,w) = sp(2,w) + qc(w,l) * fc(2,m,l)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      else                   ! Single CPU version (symmetry conserving)
      !----------------------------------------------------------------------
        do l = 1 , NLAT/2
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              if (mod(m+n,2) == 0) then ! Symmetric modes
                sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) + &
     &                                    fc(1,m,NLAT+1-l))*sfilt(n)
                sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) + &
     &                                    fc(2,m,NLAT+1-l))*sfilt(n)
              else                      ! Antisymmetric modes
                sp(1,w) = sp(1,w) + qc(w,l) * (fc(1,m,l) - &
     &                                    fc(1,m,NLAT+1-l))*sfilt(n)
                sp(2,w) = sp(2,w) + qc(w,l) * (fc(2,m,l) - &
     &                                    fc(2,m,NLAT+1-l))*sfilt(n)
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
      
      subroutine sp2fc(sp,fc,NLAT,NLON,NTRU,NLEV,nfs) ! Spectral to Fourier
      implicit none
      
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: nfs
      
      real (kind=8), intent(in ) :: sp(2,(NTRU+1)*(NTRU+2)/2) ! Coefficients of spherical harmonics
!f2py intent(in) :: sp
      real (kind=8), intent(out) :: fc(2,NLON/2,NLAT) ! Fourier coefficients
!f2py intent(out) :: fc
      
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: w ! Loop index for spectral mode
      
       
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
      
      fc(:,:,:) = 0.0
      
      do l = 1 , NLPP
         w = 1  
         do m = 1 , NTP1
            do n = m , NTP1
               fc(1,m,l) = fc(1,m,l) + qi(w,l) * sp(1,w)*sfilt(n) !72
               fc(2,m,l) = fc(2,m,l) + qi(w,l) * sp(2,w)*sfilt(n)
               w = w + 1
            enddo ! n
         enddo ! m
      enddo ! l
      return
      end      
      
      
      ! ================
      ! SUBROUTINE SP3FC
      ! ================
      
      
      subroutine sp3fc(spp,fc,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real (kind=8), intent(in ) :: spp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(in ) :: spp
      real (kind=8):: sp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8), intent(out) :: fc(2,NLON/2,NLAT, NLEV) ! Fourier coefficients
!f2py intent(out) :: fc
      
      do v = 1 , NLEV
         do k = 1, (NTRU+1)*(NTRU+2)/2
           sp(1,k) = spp(2*k-1,v)
           sp(2,k) = spp(2*k  ,v)
         enddo
         call sp2fc(sp,fcc,NLAT,NLON,NTRU,NLEV,nfs)
         fc(:,:,:,v) = fcc(:,:,:)
      enddo
      return
      end      
            

      ! ===================
      ! SUBROUTINE SP2FCDMU
      ! ===================
      
      subroutine sp2fcdmu(sp,fc,NLAT,NLON,NTRU,NLEV,nfs) ! Spectral to Fourier d/dmu
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(in) :: sp(2,(NTRU+1)*(NTRU+2)/2)! Coefficients of spherical harmonics
!f2py intent(in ) :: sp
      real(kind=8), intent(out):: fc(2,NLON/2,NLAT) ! Fourier coefficients
!f2py intent(out) :: fc      
      
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: w ! Loop index for spectral mode
        
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      fc(:,:,:) = 0.0
      
      do l = 1 , NLPP
         w = 1  
         do m = 1 , NTP1
            do n = m , NTP1
               fc(1,m,l) = fc(1,m,l) + qj(w,l) * sp(1,w)*sfilt(n)
               fc(2,m,l) = fc(2,m,l) + qj(w,l) * sp(2,w)*sfilt(n)
               w = w + 1
            enddo ! n
         enddo ! m
      enddo ! l
      return
      end


      ! ================
      ! SUBROUTINE DV2UV        !SP->GP
      ! ================
      
      subroutine dv2uv(sd,sz,pu,pv,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(in)  :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in) :: sd
      real(kind=8), intent(in)  :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in) :: sz
      real(kind=8)  :: pd(2,(NTRU+1)*(NTRU+2)/2,NLEV)
      real(kind=8)  :: pz(2,(NTRU+1)*(NTRU+2)/2,NLEV)
      real(kind=8), intent(out) :: pu(2,NLON/2,NLAT,NLEV)
!f2py intent(out) :: pu
      real(kind=8), intent(out) :: pv(2,NLON/2,NLAT,NLEV)
!f2py intent(out) :: pv
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: v ! Loop index for level
      integer :: w ! Loop index for spectral mode
      integer :: k
         
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      do v=1,NLEV
        do k=1,(NTRU+1)*(NTRU+2)/2
          pd(1,k,v) = sd(2*k-1,v)
          pd(2,k,v) = sd(2*K  ,v)
          pz(1,k,v) = sz(2*k-1,v)
          pz(2,k,v) = sz(2*k  ,v)
        enddo
      enddo
     
      pu(:,:,:,:) = 0.0
      pv(:,:,:,:) = 0.0
      
      do v = 1 , NLEV
!         zsave = pz(1,2,v)
!         pz(1,2,v) = zsave - plavor
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              pu(1,m,l,v)=pu(1,m,l,v)+qv(w,l)*pz(1,w,v)*sfilt(n) &
     &                               +qu(w,l)*pd(2,w,v)*sfilt(n)
              pu(2,m,l,v)=pu(2,m,l,v)+qv(w,l)*pz(2,w,v)*sfilt(n) &
     &                               -qu(w,l)*pd(1,w,v)*sfilt(n)
              pv(1,m,l,v)=pv(1,m,l,v)+qu(w,l)*pz(2,w,v)*sfilt(n) &
     &                               -qv(w,l)*pd(1,w,v)*sfilt(n)
              pv(2,m,l,v)=pv(2,m,l,v)-qu(w,l)*pz(1,w,v)*sfilt(n) &
     &                               -qv(w,l)*pd(2,w,v)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
!         pz(1,2,v) = zsave
      enddo ! jv
      return
      end


      ! ================
      ! SUBROUTINE UV2DV        !GP->SP
      ! ================
      
      subroutine uv2dv(pu,pv,pd,pz,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs
      
      real(kind=8), intent(out) :: pd(2,(NTRU+1)*(NTRU+2)/2,NLEV)
!f2py intent(out) :: pd
      real(kind=8), intent(out) :: pz(2,(NTRU+1)*(NTRU+2)/2,NLEV)
!f2py intent(out) :: pz
      real(kind=8), intent(in) :: pu(2,NLON/2,NLAT,NLEV)
!f2py intent(in) :: pu
      real(kind=8), intent(in) :: pv(2,NLON/2,NLAT,NLEV)
!f2py intent(in) :: pv
      
      integer :: k ! Loop index for southern latitude
      integer :: l ! Loop index for latitude
      integer :: m ! Loop index for zonal wavenumber m
      integer :: n ! Loop index for total wavenumber n
      integer :: v ! Loop index for level
      integer :: w ! Loop index for spectral mode
         
      real (kind=8):: qi((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) = Associated Legendre Polynomials
      real (kind=8):: qj((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) = Used for d/d(mu)
      real (kind=8):: qc((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd              used in fc2sp
      real (kind=8):: qe((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(mn,) * gwd / cos2       used in mktend
      real (kind=8):: qm((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * m   used in mktend
      real (kind=8):: qq((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) * gwd / cos2 * n * (n+1) / 2  "
      real (kind=8):: qu((NTRU+1)*(NTRU+2)/2,NLAT) ! P(m,n) / (n*(n+1)) * m    used in dv2uv
      real (kind=8):: qv((NTRU+1)*(NTRU+2)/2,NLAT) ! Q(m,n) / (n*(n+1))        used in dv2uv
                      
      integer NLPP, NHOR, NUGP, NPGP, NLEM, NLEP, NLSQ, NTP1
      integer NRSP, NCSP, NSPP, NESP, NVCT
      
      real (kind=8) :: EZ
      real (kind=8) :: PI
      real (kind=8) :: TWOPI
      
      real (kind=8) :: sfilt(NTRU+1)
      
      EZ     = 1.63299310207D0
      PI     = 3.14159265359D0
      TWOPI  = PI * PI
      
      NLPP = NLAT         ! Latitudes per process
      NHOR = NLON * NLPP        ! Horizontal part
      NUGP = NLON * NLAT        ! Number of gridpoints
      NPGP = NLON * NLAT / 2    ! Dimension of packed fields
      NLEM = NLEV - 1           ! Levels - 1
      NLEP = NLEV + 1           ! Levels + 1
      NLSQ = NLEV * NLEV        ! Levels squared
      NTP1 = NTRU + 1           ! Truncation + 1
      NRSP =(NTRU+1)*(NTRU+2)   ! No of real global    modes
      NCSP = NRSP / 2           ! No of complex global modes
      NSPP = (NRSP+1-1)/1 ! Modes per process
      NESP = NSPP * 1        ! Dim of spectral fields
      NVCT = 2 * (NLEV+1)       ! Dim of Vert. Coord. Tab
      
      call legini(NLAT,NLON,NTRU,NLEV,qi,qj,qc,qe,qm,qq,qu,qv, &
     &            sfilt,nfs)
     
      pd(:,:,:) = 0.0
      pz(:,:,:) = 0.0
      
      if (NLPP < NLAT) then  ! Universal (parallel executable) version
      !----------------------------------------------------------------------
      do v = 1 , NLEV
        do l = 1 , NLPP
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              pz(1,w,v) = pz(1,w,v)+qe(w,l)*pu(1,m,l,v)*sfilt(n)&
     &                             -qm(w,l)*pv(2,m,l,v)*sfilt(n)
              pz(2,w,v) = pz(2,w,v)+qe(w,l)*pu(2,m,l,v)*sfilt(n)&
     &                             +qm(w,l)*pv(1,m,l,v)*sfilt(n)
              pd(1,w,v) = pd(1,w,v)-qe(w,l)*pv(1,m,l,v)*sfilt(n)&
     &                             -qm(w,l)*pu(2,m,l,v)*sfilt(n)
              pd(2,w,v) = pd(2,w,v)-qe(w,l)*pv(2,m,l,v)*sfilt(n)&
     &                             +qm(w,l)*pu(1,m,l,v)*sfilt(n)
              w = w + 1
            enddo ! n
          enddo ! m
        enddo ! l
      enddo ! v
      else                   ! Single CPU version (symmetry conserving)
      !----------------------------------------------------------------------
      do v = 1 , NLEV
        do l = 1 , NLAT/2
          k = NLAT+1-l
          w = 1
          do m = 1 , NTP1
            do n = m , NTP1
              if (mod(m+n,2) == 0) then ! symmetric -----------------
                pz(1,w,v) = pz(1,w,v) + qe(w,l) * &
     &                           (pu(1,m,l,v)-pu(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pv(2,m,l,v)+pv(2,m,k,v))*sfilt(n)
                pz(2,w,v) = pz(2,w,v) + qe(w,l) *  &
     &                           (pu(2,m,l,v)-pu(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pv(1,m,l,v)+pv(1,m,k,v))*sfilt(n)
                pd(1,w,v) = pd(1,w,v) - qe(w,l) *  &
     &                           (pv(1,m,l,v)-pv(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pu(2,m,l,v)+pu(2,m,k,v))*sfilt(n)
                pd(2,w,v) = pd(2,w,v) - qe(w,l) *  &
     &                           (pv(2,m,l,v)-pv(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pu(1,m,l,v)+pu(1,m,k,v))*sfilt(n)
              else ! ---------------- antisymmetric -----------------
                pz(1,w,v) = pz(1,w,v) + qe(w,l) *  &
     &                           (pu(1,m,l,v)+pu(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pv(2,m,l,v)-pv(2,m,k,v))*sfilt(n)
                pz(2,w,v) = pz(2,w,v) + qe(w,l) *  &
     &                           (pu(2,m,l,v)+pu(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pv(1,m,l,v)-pv(1,m,k,v))*sfilt(n)
                pd(1,w,v) = pd(1,w,v) - qe(w,l) *  &
     &                           (pv(1,m,l,v)+pv(1,m,k,v))*sfilt(n) &
     &               - qm(w,l) * (pu(2,m,l,v)-pu(2,m,k,v))*sfilt(n)
                pd(2,w,v) = pd(2,w,v) - qe(w,l) *  &
     &                           (pv(2,m,l,v)+pv(2,m,k,v))*sfilt(n) &
     &               + qm(w,l) * (pu(1,m,l,v)-pu(1,m,k,v))*sfilt(n)
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
      ! SUBROUTINE INIGAU
      ! =================
      
      subroutine inigau(klat,pz0,pzw)        ! pz0 & pzw are (kind=8) reals !!!
      implicit none
      integer ,intent(IN)        :: klat       ! Number of Gaussian latitudes
      real (kind=8), intent(out) :: pz0(klat)  ! Gaussian abscissas
      real (kind=8), intent(out) :: pzw(klat)  ! Gaussian weights
      integer                  :: jlat       ! Latitudinal loop index
      integer                  :: jiter      ! Iteration loop index
      integer      , parameter :: NITER = 50 ! Maximum # of iterations
      real (kind=8), parameter :: PI    =  3.14159265358979_8
      real (kind=8), parameter :: ZEPS  =  1.0e-16 ! Convergence criterion
      real (kind=8) :: z0,z1,z2,z3,z4,z5
      real (kind=8) :: ql,qld
      
      ! Compute Gaussian abscissas & weights
      
      z0 = PI / (2*klat+1)
      z1 = 1.0_8 / (klat*klat*8)
      z4 = 2.0_8 / (klat*klat)
      
      do jlat = 1 , klat/2
         z2 = z0 * (2*jlat - 0.5_8)
         z2 = cos(z2 + z1 / tan(z2))
         do jiter = 1 , NITER
            z3 = ql(klat,z2) * qld(klat,z2)
            z2 = z2 - z3
            if (abs(z3) < ZEPS) exit ! converged
         enddo ! jiter
         z5 = ql(klat-1,z2) / sqrt(klat - 0.5_8)
         pz0(jlat) = z2
         pzw(jlat) = z4 * (1.0_8 - z2 * z2) / (z5 * z5)
         pz0(klat-jlat+1) = -z2
         pzw(klat-jlat+1) = pzw(jlat)
      enddo ! jlat
      
      return
      end subroutine inigau
      
        
      ! ===========
      ! FUNCTION QL
      ! ===========
      
      real (kind=8) function ql(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8) :: z0,z1,z2,z3,z4
      integer :: j
      z0 = acos(p)
      z1 = 1.0
      z2 = 0.0
      do j = k , 0 , -2
         z3 = z1 * cos(z0 * j)
         z2 = z2 + z3
         z4 = (k-j+1) * (k+j) * 0.5_8
         z1 = z1 * z4 / (z4 + (j-1))
      enddo ! j
      if (mod(k,2) == 0) z2 = z2 - 0.5_8 * z3
      
      z0 = sqrt(2.0_8)
      do j = 1 ,k
         z0 = z0 * sqrt(1.0_8 - 0.25_8 / (j*j))
      enddo ! j
      ql = z0 * z2
      return
      end function ql
      
      ! ============
      ! FUNCTION QLD
      ! ============
      
      real (kind=8) function qld(k,p)
      implicit none
      integer      , intent(IN) :: k
      real (kind=8), intent(IN) :: p
      real (kind=8) :: z
      real (kind=8) :: ql
      
      z = p * ql(k,p) - sqrt((k + k + 1.0_8) / &
     &            (k + k - 1.0_8)) * ql(k-1,p)
      qld = (p * p - 1.0_8) / (k * z)
      
      return
      end function qld
      
!     =============
!     MODULE FFTMOD
!     =============

!     alternate module with factors 2, 3, 4, 5, 6, 8
!     this is a FORTRAN-90 version of the old FFT991 package
!     use this for resolutions not supported by the <fftmod>,
!     e.g. T63 (NLAT=96), T106 (NLAT=160)

      module fftmod
      parameter(NRES = 8)
      integer :: nallowed(NRES) = (/ 4, 64, 96, 128, 192, 256, 320, 384 /)
!     T1   - N4   : 4
!     T21  - N64  : 8-4-2
!     T31  - N96  : 8-4-3
!     T42  - N128 : 8-4-4
!     T63  - N192 : 8-6-4
!     T85  - N256 : 8-4-4-2
!     T106 - N320 : 8-5-4-2
!     T127 - N384 : 8-4-4-3
!     T170 - N512 : 8-4-4-4
      integer :: lastn = 0
      integer :: ifax(10)
      real (kind=8) ,allocatable :: trigs(:)
      end module fftmod
     

!     ================
!     SUBROUTINE GP2FC
!     ================

      subroutine gp2fc(a,c,n,lot)
      use fftmod
      real (kind=8), intent(in) :: a
      real (kind=8), intent(out):: c
      integer, intent(in) :: n
      integer, intent(in) :: lot
      
      real (kind=8) worka
      real (kind=8) workb

      dimension a(n,lot)
      dimension c(n,lot)
      dimension worka(n+2,lot)
      dimension workb(n+2,lot)
      
      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      worka(1:n,:) = a(:,:)
      worka(n+1:n+2,:) = 0.0

      call fft991(worka,workb,trigs,ifax,1,n+2,n,lot,-1)

      c(:,:) = worka(1:n,:)

      return
      end subroutine gp2fc
      
    
!     ================
!     SUBROUTINE FC2GP
!     ================

      subroutine fc2gp(a,c,n,lot)
      use fftmod
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out):: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      real (kind=8) worka
      real (kind=8) workb
      
      dimension a(n,lot)
      dimension c(n,lot)
      dimension worka(n+2,lot)
      dimension workb(n+2,lot)

      if (n /= lastn) then
         if (allocated(trigs)) deallocate(trigs)
         allocate(trigs(n))
         lastn = n
         call fftini(n)
      endif

      worka(1:n,:) = a(:,:)
      worka(n+1:n+2,:) = 0.0

      call fft991(worka,workb,trigs,ifax,1,n+2,n,lot,1)

      c(:,:) = worka(1:n,:)*1.4142135623730951 ! *sqrt(2)

      return
      end subroutine fc2gp
  
          
!     ================
!     SUBROUTINE FC3GP
!     ================
      
      
      subroutine fc3gp(a,c,n,lot,NLEV)
      implicit none
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out) :: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      integer, intent(in) :: NLEV
!f2py intent(in) :: NLEV

      integer jlev
      real (kind=8) :: dd
      
      dimension dd(n,lot)
      dimension a(n,lot,NLEV)
      dimension c(n,lot,NLEV)
      
      do jlev=1,NLEV
        call fc2gp(a(:,:,jlev),dd,n,lot)
        c(:,:,jlev) = dd(:,:)
      enddo
      
      return
      end subroutine
       
!     ================
!     SUBROUTINE GP3FC
!     ================
      
      
      subroutine gp3fc(a,c,n,lot,NLEV)
      implicit none
      real (kind=8), intent(in) :: a
!f2py intent(in) :: a
      real (kind=8), intent(out) :: c
!f2py intent(out) :: c
      integer, intent(in) :: n
!f2py intent(in) :: n
      integer, intent(in) :: lot
!f2py intent(in) :: lot
      integer, intent(in) :: NLEV
!f2py intent(in) :: NLEV

      integer jlev
      real (kind=8) :: dd
      
      dimension dd(n,lot)
      dimension a(n,lot,NLEV)
      dimension c(n,lot,NLEV)
      
      do jlev=1,NLEV
        call gp2fc(a(:,:,jlev),dd,n,lot)
        c(:,:,jlev) = dd(:,:)
      enddo
      
      return
      end subroutine
        
!     ================
!     SUBROUTINE SP2GP
!     ================
      
      subroutine sp2gp(sp,gp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs !1 or 0, whether use physics filter
      
      real (kind=8), intent(in ) :: sp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(in ) :: sp
      real (kind=8):: spp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8):: fcl(NLON,NLAT)
      real (kind=8):: gpp(NLON,NLAT)
      real (kind=8), intent(out) :: gp(NLON, NLAT, NLEV) ! Fourier coefficients
!f2py intent(out) :: gp
      
      do v = 1 , NLEV
         do k = 1, (NTRU+1)*(NTRU+2)/2
           spp(1,k) = sp(2*k-1,v)
           spp(2,k) = sp(2*k  ,v)
         enddo
         call sp2fc(spp,fcc,NLAT,NLON,NTRU,NLEV,nfs)
         do k = 1, NLON/2
           fcl(2*k-1,:) = fcc(1,k,:)
           fcl(2*k  ,:) = fcc(2,k,:)
         enddo
         call fc2gp(fcl,gpp,NLON,NLAT)
         gp(:,:,v) = gpp(:,:)
!          write(*,*) "Computed layer",v
      enddo
      
      return
      end subroutine    
      
   
!     ================
!     SUBROUTINE SPVGP
!     ================
      
      subroutine spvgp(sd,sz,rdcostheta,gu,gv,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k, l ! Loop index for level
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: nfs
      real (kind=8), intent(in ) :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in ) :: sd
      real (kind=8), intent(in ) :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(in ) :: sz
      real (kind=8), intent(in ) :: rdcostheta(NLAT) !radius / cos(latitude)
!f2py intent(in ) :: rdcostheta
      real (kind=8), intent(out) :: gu(NLON, NLAT, NLEV)
!f2py intent(out) :: gu
      real (kind=8), intent(out) :: gv(NLON, NLAT, NLEV)
!f2py intent(out) :: gv

!       real (kind=8) sdd(2, (NTRU+1)*(NTRU+2)/2, NLEV)
!       real (kind=8) szz(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) fuu(2, NLON/2, NLAT, NLEV)
      real (kind=8) fvv(2, NLON/2, NLAT, nLEV)
      real (kind=8) fup(NLON, NLAT)
      real (kind=8) fvp(NLON, NLAT)
      real (kind=8) gup(NLON, NLAT)
      real (kind=8) gvp(NLON, NLAT)
      
!       do v = 1, NLEV
!         do k = 1, (NTRU+1)*(NTRU+2)/2
!           sdd(1,k,v) = sd(2*k-1,v)
!           sdd(2,k,v) = sd(2*k  ,v)
!           szz(1,k,v) = sz(2*k-1,v)
!           szz(2,k,v) = sz(2*k  ,v)
!         enddo
!       enddo
      
      call dv2uv(sd,sz,fuu,fvv,NLAT,NLON,NTRU,NLEV,nfs)
      
      do v = 1, NLEV
        do l = 1, NLAT
          do k = 1, NLON/2
            fup(2*k-1,l) = fuu(1,k,l,v)
            fup(2*k  ,l) = fuu(2,k,l,v)
            fvp(2*k-1,l) = fvv(1,k,l,v)
            fvp(2*k  ,l) = fvv(2,k,l,v)
          enddo
        enddo
        call fc2gp(fup,gup,NLON,NLAT)
        call fc2gp(fvp,gvp,NLON,NLAT)
        do l = 1, NLAT
          gu(:,l,v) = gup(:,l)*rdcostheta(l)
          gv(:,l,v) = gvp(:,l)*rdcostheta(l)
        enddo
      enddo
      
      return
      end subroutine
      
      
!     ================
!     SUBROUTINE GPVSP
!     ================
      
      subroutine gpvsp(gu,gv,costhetadr,sd,sz,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k, l ! Loop index for level
      integer, intent(in ) :: NLEV
      integer, intent(in ) :: NLON
      integer, intent(in ) :: NTRU
      integer, intent(in ) :: NLAT
      integer, intent(in ) :: nfs  
      real (kind=8), intent(out) :: sd((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(out) :: sd
      real (kind=8), intent(out) :: sz((NTRU+1)*(NTRU+2),NLEV)
!f2py intent(out) :: sz
      real (kind=8), intent(in ) :: gu(NLON, NLAT, NLEV)
!f2py intent(in ) :: gu
      real (kind=8), intent(in ) :: gv(NLON, NLAT, NLEV)
!f2py intent(in ) :: gv
      real (kind=8), intent(in ) :: costhetadr(NLAT) !cos(latitude) / radius
!f2py intent(in ) :: costhetadr    

      real (kind=8) sdd(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) szz(2, (NTRU+1)*(NTRU+2)/2, NLEV)
      real (kind=8) guu(2, NLON/2, NLAT, NLEV)
      real (kind=8) gvv(2, NLON/2, NLAT, nLEV)
      real (kind=8) gur(NLON, NLAT)
      real (kind=8) gvr(NLON, NLAT)
      real (kind=8) fuu(NLON, NLAT)
      real (kind=8) fvv(NLON, NLAT)
      
      do v = 1, NLEV
        do l = 1, NLAT
          gur(:,l) = gu(:,l,v)*costhetadr(l)/1.4142135623730951
          gvr(:,l) = gv(:,l,v)*costhetadr(l)/1.4142135623730951
        enddo
        call gp2fc(gur,fuu,NLON,NLAT)
        call gp2fc(gvr,fvv,NLON,NLAT)
        do l = 1, NLAT
          do k = 1, NLON/2
            guu(1,k,l,v) = fuu(2*k-1,l)
            guu(2,k,l,v) = fuu(2*k  ,l)
            gvv(1,k,l,v) = fvv(2*k-1,l)
            gvv(2,k,l,v) = fvv(2*k  ,l)
          enddo
        enddo
      enddo
      
      call uv2dv(guu,gvv,sdd,szz,NLAT,NLON,NTRU,NLEV,nfs)
      
      do v = 1, NLEV
        do k = 1, (NTRU+1)*(NTRU+2)/2
          sd(2*k-1,v) = sdd(1,k,v)
          sd(2*k  ,v) = sdd(2,k,v)
          sz(2*k-1,v) = szz(1,k,v)
          sz(2*k  ,v) = szz(2,k,v)
        enddo
      enddo
      
      return
      end subroutine
            

      
!     ================
!     SUBROUTINE GP2SP
!     ================
      
      subroutine gp2sp(gp,sp,NLAT,NLON,NTRU,NLEV,nfs)
      implicit none
      integer :: v, k ! Loop index for level
      integer, intent(in) :: NLEV
      integer, intent(in) :: NLON
      integer, intent(in) :: NTRU
      integer, intent(in) :: NLAT
      integer, intent(in) :: nfs !1 or 0, whether use physics filter
      
      real (kind=8), intent(out ) :: sp((NTRU+1)*(NTRU+2), NLEV) 
!f2py intent(out) :: sp
      real (kind=8):: spp(2,(NTRU+1)*(NTRU+2)/2)
      real (kind=8):: fcc(2,NLON/2,NLAT)
      real (kind=8):: fcl(NLON,NLAT)
      real (kind=8):: gpp(NLON,NLAT)
      real (kind=8), intent(in) :: gp(NLON, NLAT, NLEV) ! Gridpoint variable
!f2py intent(in ) :: gp
      
      do v = 1 , NLEV
         call gp2fc(gp(:,:,v)/1.4142135623730951,fcl,NLON,NLAT)
         do k = 1, NLON/2
            fcc(1,k,:) = fcl(2*k-1,:) 
            fcc(2,k,:) = fcl(2*k  ,:) 
         enddo
         call fc2sp(fcc,spp,NLAT,NLON,NTRU,NLEV,nfs)
         do k = 1, (NTRU+1)*(NTRU+2)/2
            sp(2*k-1,v) = spp(1,k) 
            sp(2*k  ,v) = spp(2,k) 
         enddo
!          write(*,*) "Computed layer",v
      enddo
      
      return
      end subroutine 
      

!     =================
!     SUBROUTINE FFTINI
!     =================

      subroutine fftini(n)
      use fftmod
      integer, intent(in) :: n
      logical labort

!     check for allowed values of n

      labort = .true.
      do j = 1 , NRES
         if (n == nallowed(j)) labort = .false.
      enddo

      if (labort) then
         write (*,*) '*** FFT does not support n = ',n,' ***'
         write (*,*) 'Following resolutions may be used:'
         write (*,*) '----------------------------------'
         do j = 1 , NRES
            write (*,1000) nallowed(j), nallowed(j)/2, nallowed(j)/3
         enddo
         stop
      endif
 1000 format(' NLON=',I5,'  NLAT=',I5,'  NTRU=',I5)

      call set99(trigs,ifax,n) ! Factorization of n and sine/cosine values

      return
      end subroutine fftini
      
       

      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!     SUBROUTINE FFT991 - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     SUPERSEDES PREVIOUS ROUTINE 'FFT991'
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*MIN(LOT,64)
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
!     IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
      REAL (KIND=8), INTENT(INOUT) :: A
      REAL (KIND=8), INTENT(INOUT) :: WORK
      REAL (KIND=8), INTENT(IN) :: TRIGS
      INTEGER, INTENT(IN) :: IFAX
      INTEGER, INTENT(IN) :: JUMP
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: LOT
      INTEGER, INTENT(IN) :: ISIGN
      
      DIMENSION A(*),WORK(*),TRIGS(N),IFAX(10)
!
      NFAX=IFAX(1)
      NX=N+1
      IF (MOD(N,2).EQ.1) NX=N
      NBLOX=1+(LOT-1)/64
      NVEX=LOT-(NBLOX-1)*64
      IF (ISIGN.EQ.-1) GO TO 300
!
!     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
!     -----------------------------------------
100   CONTINUE
      ISTART=1
      DO 220 NB=1,NBLOX
        IA=ISTART
        I=ISTART
        DO 110 J=1,NVEX
          A(I+INC)=0.5*A(I)
          I=I+JUMP
110       CONTINUE
        IF (MOD(N,2).EQ.1) GO TO 130
        I=ISTART+N*INC
        DO 120 J=1,NVEX
          A(I)=0.5*A(I)
          I=I+JUMP
120       CONTINUE
130     CONTINUE
        IA=ISTART+INC
        LA=1
        IGO=+1
!
        DO 160 K=1,NFAX
          IFAC=IFAX(K+1)
          IERR=-1
          IF (IGO.EQ.-1) GO TO 140
          CALL RPASSM(A(IA),A(IA+LA*INC),WORK(1),WORK(IFAC*LA+1),TRIGS, &
                   INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
          GO TO 150
140       CONTINUE
          CALL RPASSM(WORK(1),WORK(LA+1),A(IA),A(IA+IFAC*LA*INC),TRIGS, &
                   1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
150       CONTINUE
          IF (IERR.NE.0) GO TO 500
          LA=IFAC*LA
          IGO=-IGO
          IA=ISTART
160       CONTINUE
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (MOD(NFAX,2).EQ.0) GO TO 190
        IBASE=1
        JBASE=IA
        DO 180 JJ=1,NVEX
          I=IBASE
          J=JBASE
          DO 170 II=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
170         CONTINUE
          IBASE=IBASE+NX
          JBASE=JBASE+JUMP
180       CONTINUE
190     CONTINUE
!
!     FILL IN ZEROS AT END
!     --------------------
        IX=ISTART+N*INC
        DO 210 J=1,NVEX
          A(IX)=0.0
          A(IX+INC)=0.0
          IX=IX+JUMP
210       CONTINUE
!
        ISTART=ISTART+NVEX*JUMP
        NVEX=64
220     CONTINUE
      RETURN
!
!     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
!     -----------------------------------------
300   CONTINUE
      ISTART=1
      DO 410 NB=1,NBLOX
        IA=ISTART
        LA=N
        IGO=+1
!
        DO 340 K=1,NFAX
          IFAC=IFAX(NFAX+2-K)
          LA=LA/IFAC
          IERR=-1
          IF (IGO.EQ.-1) GO TO 320
          CALL QPASSM(A(IA),A(IA+IFAC*LA*INC),WORK(1),WORK(LA+1),TRIGS, &
                   INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
          GO TO 330
320       CONTINUE
          CALL QPASSM(WORK(1),WORK(IFAC*LA+1),A(IA),A(IA+LA*INC),TRIGS, &
                   1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
330       CONTINUE
          IF (IERR.NE.0) GO TO 500
          IGO=-IGO
          IA=ISTART+INC
340       CONTINUE
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (MOD(NFAX,2).EQ.0) GO TO 370
        IBASE=1
        JBASE=IA
        DO 360 JJ=1,NVEX
          I=IBASE
          J=JBASE
          DO 350 II=1,N
            A(J)=WORK(I)
            I=I+1
            J=J+INC
350         CONTINUE
          IBASE=IBASE+NX
          JBASE=JBASE+JUMP
360       CONTINUE
370     CONTINUE
!
!     SHIFT A(0) & FILL IN ZERO IMAG PARTS
!     ------------------------------------
        IX=ISTART
        DO 380 J=1,NVEX
          A(IX)=A(IX+INC)
          A(IX+INC)=0.0
          IX=IX+JUMP
380       CONTINUE
        IF (MOD(N,2).EQ.1) GO TO 400
        IZ=ISTART+(N+1)*INC
        DO 390 J=1,NVEX
          A(IZ)=0.0
          IZ=IZ+JUMP
390       CONTINUE
400     CONTINUE
!
        ISTART=ISTART+NVEX*JUMP
        NVEX=64
410     CONTINUE
      RETURN
!
!     ERROR MESSAGES
!     --------------
500   CONTINUE
      GO TO (510,530,550) IERR
510   CONTINUE
      WRITE(6,520) NVEX
520   FORMAT(16H1VECTOR LENGTH =,I4,17H, GREATER THAN 64)
      GO TO 570
530   CONTINUE
      WRITE(6,540) IFAC
540   FORMAT( 9H1FACTOR =,I3,17H, NOT CATERED FOR)
      GO TO 570
550   CONTINUE
      WRITE(6,560) IFAC
560   FORMAT(9H1FACTOR =,I3,31H, ONLY CATERED FOR IF LA*IFAC=N)
570   CONTINUE
      RETURN
      END



      SUBROUTINE SET99(TRIGS,IFAX,N)
!     SUBROUTINE SET99 - COMPUTES FACTORS OF N & TRIGONOMETRIC
!     FUNCTINS REQUIRED BY FFT99 & FFT991
!
      REAL (KIND=8), INTENT(INOUT) :: TRIGS
      INTEGER, INTENT(INOUT) :: IFAX
      INTEGER, INTENT(IN) :: N
      DIMENSION TRIGS(N),IFAX(10),JFAX(10),LFAX(7)
      DATA LFAX/6,8,5,4,3,2,1/
      IXXX=1
!
      DEL=4.0*ASIN(1.0)/FLOAT(N)
      NIL=0
      NHL=(N/2)-1
      DO 10 K=NIL,NHL
        ANGLE=FLOAT(K)*DEL
        TRIGS(2*K+1)=COS(ANGLE)
        TRIGS(2*K+2)=SIN(ANGLE)
10      CONTINUE
!
!     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
!     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
      K=0
      L=1
20    CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
25    CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
30    CONTINUE
      L=L+1
      IFAC=LFAX(L)
      IF (IFAC.GT.1) GO TO 20
!
      WRITE(6,40) N
40    FORMAT(4H1N =,I4,27H - CONTAINS ILLEGAL FACTORS)
      RETURN
!
!     NOW REVERSE ORDER OF FACTORS
50    CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
        IFAX(NFAX+2-I)=JFAX(I)
60      CONTINUE
      IFAX(10)=N
      RETURN
      END




      SUBROUTINE QPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA,IERR)
!     SUBROUTINE QPASSM - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER ANALYSIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A(IFAC*LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA = N/(PRODUCT OF FACTORS USED SO FAR)
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!
!-----------------------------------------------------------------------
!
      REAL (KIND=8), INTENT(IN) :: A
      REAL (KIND=8), INTENT(IN) :: B
      REAL (KIND=8), INTENT(OUT) :: C
      REAL (KIND=8), INTENT(OUT) :: D
      REAL (KIND=8), INTENT(IN) :: TRIGS
      INTEGER, INTENT(IN) :: INC1
      INTEGER, INTENT(IN) :: INC2
      INTEGER, INTENT(IN) :: INC3
      INTEGER, INTENT(IN) :: INC4
      INTEGER, INTENT(IN) :: LOT
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: IFAC
      INTEGER, INTENT(IN) :: LA
      INTEGER, INTENT(OUT) :: IERR
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)
!
      REAL (KIND=8) SIN36,SIN72,QRT5,SIN60
      parameter(SIN36 = 0.587785252292473)
      parameter(SIN72 = 0.951056516295154)
      parameter(QRT5  = 0.559016994374947)
      parameter(SIN60 = 0.866025403784437)
!
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBAD=1
      IF (LOT.GT.64) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
!
!     CODING FOR FACTOR 2
!     -------------------
200   CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
!
      IF (LA.EQ.M) GO TO 290
!
      DO 220 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 210 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=A(IA+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
210       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
220     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JB) GO TO 260
      DO 250 K=LA,KSTOP,LA
        KB=K+K
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        JBASE=0
        DO 240 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 230 IJK=1,LOT
            C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
            C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
            D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
            D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
            I=I+INC3
            J=J+INC4
230         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
240       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB-JINK
250     CONTINUE
      IF (JA.GT.JB) GO TO 900
260   CONTINUE
      JBASE=0
      DO 280 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 270 IJK=1,LOT
          C(JA+J)=A(IA+I)
          D(JA+J)=-A(IB+I)
          I=I+INC3
          J=J+INC4
270       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
280     CONTINUE
      GO TO 900
!
290   CONTINUE
      Z=1.0/FLOAT(N)
      DO 294 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 292 IJK=1,LOT
          C(JA+J)=Z*(A(IA+I)+A(IB+I))
          C(JB+J)=Z*(A(IA+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
292       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
294     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 3
!     -------------------
300   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB
!
      IF (LA.EQ.M) GO TO 390
!
      DO 320 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 310 IJK=1,LOT
          C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
          C(JB+J)=A(IA+I)-0.5*(A(IB+I)+A(IC+I))
          D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
310       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
320     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JC) GO TO 360
      DO 350 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        JBASE=0
        DO 340 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 330 IJK=1,LOT
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
            A2=A(IA+I)-0.5*A1
            B2=B(IA+I)-0.5*B1
            A3=SIN60*((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
            B3=SIN60*((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
            C(JA+J)=A(IA+I)+A1
            D(JA+J)=B(IA+I)+B1
            C(JB+J)=A2+B3
            D(JB+J)=B2-A3
            C(JC+J)=A2-B3
            D(JC+J)=-(B2+A3)
            I=I+INC3
            J=J+INC4
330         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
340       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC-JINK
350     CONTINUE
      IF (JA.GT.JC) GO TO 900
360   CONTINUE
      JBASE=0
      DO 380 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 370 IJK=1,LOT
          C(JA+J)=A(IA+I)+0.5*(A(IB+I)-A(IC+I))
          D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
          C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
          I=I+INC3
          J=J+INC4
370       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
380     CONTINUE
      GO TO 900
!
390   CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 394 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 392 IJK=1,LOT
          C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
          C(JB+J)=Z*(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))
          D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
392       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
394     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 4
!     -------------------
400   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JB
!
      IF (LA.EQ.M) GO TO 490
!
      DO 420 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 410 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
          C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
          C(JB+J)=A(IA+I)-A(IC+I)
          D(JB+J)=A(ID+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
410       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
420     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JC) GO TO 460
      DO 450 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        JBASE=0
        DO 440 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 430 IJK=1,LOT
            A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
            A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
            A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
            B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
            B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
            B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
            C(JA+J)=A0+A1
            C(JC+J)=A0-A1
            D(JA+J)=B0+B1
            D(JC+J)=B1-B0
            C(JB+J)=A2+B3
            C(JD+J)=A2-B3
            D(JB+J)=B2-A3
            D(JD+J)=-(B2+A3)
            I=I+INC3
            J=J+INC4
430         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
440       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC-JINK
        JD=JD-JINK
450     CONTINUE
      IF (JB.GT.JC) GO TO 900
460   CONTINUE
      SIN45=SQRT(0.5)
      JBASE=0
      DO 480 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 470 IJK=1,LOT
          C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
          C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
          D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
          D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
          I=I+INC3
          J=J+INC4
470       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
480     CONTINUE
      GO TO 900
!
490   CONTINUE
      Z=1.0/FLOAT(N)
      DO 494 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 492 IJK=1,LOT
          C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JB+J)=Z*(A(IA+I)-A(IC+I))
          D(JB+J)=Z*(A(ID+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
492       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
494     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 5
!     -------------------
500   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB
!
      IF (LA.EQ.M) GO TO 590
!
      DO 520 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 510 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=A(IA+I)-0.25*(A1+A2)
          A6=QRT5*(A1-A2)
          C(JA+J)=A(IA+I)+(A1+A2)
          C(JB+J)=A5+A6
          C(JC+J)=A5-A6
          D(JB+J)=-SIN72*A3-SIN36*A4
          D(JC+J)=-SIN36*A3+SIN72*A4
          I=I+INC3
          J=J+INC4
510       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
520     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JD) GO TO 560
      DO 550 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        JBASE=0
        DO 540 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 530 IJK=1,LOT
            A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
            A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
            A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
            A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
            B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
            B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
            B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
            B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
            A5=A(IA+I)-0.25*(A1+A2)
            A6=QRT5*(A1-A2)
            B5=B(IA+I)-0.25*(B1+B2)
            B6=QRT5*(B1-B2)
            A10=A5+A6
            A20=A5-A6
            B10=B5+B6
            B20=B5-B6
            A11=SIN72*B3+SIN36*B4
            A21=SIN36*B3-SIN72*B4
            B11=SIN72*A3+SIN36*A4
            B21=SIN36*A3-SIN72*A4
            C(JA+J)=A(IA+I)+(A1+A2)
            C(JB+J)=A10+A11
            C(JE+J)=A10-A11
            C(JC+J)=A20+A21
            C(JD+J)=A20-A21
            D(JA+J)=B(IA+I)+(B1+B2)
            D(JB+J)=B10-B11
            D(JE+J)=-(B10+B11)
            D(JC+J)=B20-B21
            D(JD+J)=-(B20+B21)
            I=I+INC3
            J=J+INC4
530         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
540       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC+JINK
        JD=JD-JINK
        JE=JE-JINK
550     CONTINUE
      IF (JB.GT.JD) GO TO 900
560   CONTINUE
      JBASE=0
      DO 580 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 570 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=A(IA+I)+0.25*(A3-A4)
          A6=QRT5*(A3+A4)
          C(JA+J)=A5+A6
          C(JB+J)=A5-A6
          C(JC+J)=A(IA+I)-(A3-A4)
          D(JA+J)=-SIN36*A1-SIN72*A2
          D(JB+J)=-SIN72*A1+SIN36*A2
          I=I+INC3
          J=J+INC4
570       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
580     CONTINUE
      GO TO 900
!
590   CONTINUE
      Z=1.0/FLOAT(N)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      DO 594 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 592 IJK=1,LOT
          A1=A(IB+I)+A(IE+I)
          A3=A(IB+I)-A(IE+I)
          A2=A(IC+I)+A(ID+I)
          A4=A(IC+I)-A(ID+I)
          A5=Z*(A(IA+I)-0.25*(A1+A2))
          A6=ZQRT5*(A1-A2)
          C(JA+J)=Z*(A(IA+I)+(A1+A2))
          C(JB+J)=A5+A6
          C(JC+J)=A5-A6
          D(JB+J)=-ZSIN72*A3-ZSIN36*A4
          D(JC+J)=-ZSIN36*A3+ZSIN72*A4
          I=I+INC3
          J=J+INC4
592       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
594     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 6
!     -------------------
600   CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
!
      IF (LA.EQ.M) GO TO 690
!
      DO 620 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 610 IJK=1,LOT
          A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
          C(JA+J)=(A(IA+I)+A(ID+I))+A11
          C(JC+J)=(A(IA+I)+A(ID+I)-0.5*A11)
          D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
          A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
          C(JB+J)=(A(IA+I)-A(ID+I))-0.5*A11
          D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
          C(JD+J)=(A(IA+I)-A(ID+I))+A11
          I=I+INC3
          J=J+INC4
610       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
620     CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JC.EQ.JD) GO TO 660
      DO 650 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        KF=KE+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        C5=TRIGS(KF+1)
        S5=TRIGS(KF+2)
        JBASE=0
        DO 640 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 630 IJK=1,LOT
            A1=C1*A(IB+I)+S1*B(IB+I)
            B1=C1*B(IB+I)-S1*A(IB+I)
            A2=C2*A(IC+I)+S2*B(IC+I)
            B2=C2*B(IC+I)-S2*A(IC+I)
            A3=C3*A(ID+I)+S3*B(ID+I)
            B3=C3*B(ID+I)-S3*A(ID+I)
            A4=C4*A(IE+I)+S4*B(IE+I)
            B4=C4*B(IE+I)-S4*A(IE+I)
            A5=C5*A(IF+I)+S5*B(IF+I)
            B5=C5*B(IF+I)-S5*A(IF+I)
            A11=(A2+A5)+(A1+A4)
            A20=(A(IA+I)+A3)-0.5*A11
            A21=SIN60*((A2+A5)-(A1+A4))
            B11=(B2+B5)+(B1+B4)
            B20=(B(IA+I)+B3)-0.5*B11
            B21=SIN60*((B2+B5)-(B1+B4))
            C(JA+J)=(A(IA+I)+A3)+A11
            D(JA+J)=(B(IA+I)+B3)+B11
            C(JC+J)=A20-B21
            D(JC+J)=A21+B20
            C(JE+J)=A20+B21
            D(JE+J)=A21-B20
            A11=(A2-A5)+(A4-A1)
            A20=(A(IA+I)-A3)-0.5*A11
            A21=SIN60*((A4-A1)-(A2-A5))
            B11=(B5-B2)-(B4-B1)
            B20=(B3-B(IA+I))-0.5*B11
            B21=SIN60*((B5-B2)+(B4-B1))
            C(JB+J)=A20-B21
            D(JB+J)=A21-B20
            C(JD+J)=A11+(A(IA+I)-A3)
            D(JD+J)=B11+(B3-B(IA+I))
            C(JF+J)=A20+B21
            D(JF+J)=A21+B20
            I=I+INC3
            J=J+INC4
630         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
640       CONTINUE
        IBASE=IBASE+IJUMP
        JA=JA+JINK
        JB=JB+JINK
        JC=JC+JINK
        JD=JD-JINK
        JE=JE-JINK
        JF=JF-JINK
650     CONTINUE
      IF (JC.GT.JD) GO TO 900
660   CONTINUE
      JBASE=0
      DO 680 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 670 IJK=1,LOT
          C(JA+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))+ SIN60*(A(IB+I)-A(IF+I))
          D(JA+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))-SIN60*(A(IC+I)+A(IE+I))
          C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
          D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
          C(JC+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I)))-SIN60*(A(IB+I)-A(IF+I))
          D(JC+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I)))+SIN60*(A(IC+I)+A(IE+I))
          I=I+INC3
          J=J+INC4
670       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
680     CONTINUE
      GO TO 900
!
690   CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 694 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 692 IJK=1,LOT
          A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
          C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
          C(JC+J)=Z*((A(IA+I)+A(ID+I))-0.5*A11)
          D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
          A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
          C(JB+J)=Z*((A(IA+I)-A(ID+I))-0.5*A11)
          D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
          C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
          I=I+INC3
          J=J+INC4
692       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
694     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 8
!     -------------------
800   CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+LA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0/FLOAT(N)
      ZSIN45=Z*SQRT(0.5)
!
      DO 820 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 810 IJK=1,LOT
          C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+ &
                   ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
          C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))- &
                   ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
          C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
          D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
          C(JB+J)=Z*(A(IA+I)-A(IE+I)) &
                   +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
          C(JD+J)=Z*(A(IA+I)-A(IE+I)) &
                   -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
          D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I))) &
                   +Z*(A(IG+I)-A(IC+I))
          D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I))) &
                   -Z*(A(IG+I)-A(IC+I))
          I=I+INC3
          J=J+INC4
810       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
820     CONTINUE
!
!     RETURN
!     ------
900   CONTINUE
      IBAD=0
910   CONTINUE
      IERR=IBAD
      RETURN
      END



      SUBROUTINE RPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA,IERR)
!     SUBROUTINE 'RPASSM' - PERFORMS ONE PASS THROUGH DATA AS PART
!     OF MULTIPLE REAL FFT (FOURIER SYNTHESIS) ROUTINE
!
!     A IS FIRST REAL INPUT VECTOR
!         EQUIVALENCE B(1) WITH A (LA*INC1+1)
!     C IS FIRST REAL OUTPUT VECTOR
!         EQUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)
!     TRIGS IS A PRECALCULATED LIST OF SINES & COSINES
!     INC1 IS THE ADDRESSING INCREMENT FOR A
!     INC2 IS THE ADDRESSING INCREMENT FOR C
!     INC3 IS THE INCREMENT BETWEEN INPUT VECTORS A
!     INC4 IS THE INCREMENT BETWEEN OUTPUT VECTORS C
!     LOT IS THE NUMBER OF VECTORS
!     N IS THE LENGTH OF THE VECTORS
!     IFAC IS THE CURRENT FACTOR OF N
!     LA IS THE PRODUCT OF PREVIOUS FACTORS
!     IERR IS AN ERROR INDICATOR:
!              0 - PASS COMPLETED WITHOUT ERROR
!              1 - LOT GREATER THAN 64
!              2 - IFAC NOT CATERED FOR
!              3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
!
!-----------------------------------------------------------------------
!
      REAL (KIND=8), INTENT(IN) :: A
      REAL (KIND=8), INTENT(IN) :: B
      REAL (KIND=8), INTENT(OUT) :: C
      REAL (KIND=8), INTENT(OUT) :: D
      REAL (KIND=8), INTENT(IN) :: TRIGS
      INTEGER, INTENT(IN) :: INC1
      INTEGER, INTENT(IN) :: INC2
      INTEGER, INTENT(IN) :: INC3
      INTEGER, INTENT(IN) :: INC4
      INTEGER, INTENT(IN) :: LOT
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: IFAC
      INTEGER, INTENT(IN) :: LA
      INTEGER, INTENT(OUT) :: IERR
      
      REAL (KIND=8) A10,A11,A20,A21,B10,B11,B20,B21
      REAL (KIND=8) SIN36,SIN72,QRT5,SIN60
      
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(N)
!
      DIMENSION A10(64),A11(64),A20(64),A21(64),B10(64),B11(64),B20(64),B21(64)
      parameter(SIN36 = 0.587785252292473)
      parameter(SIN72 = 0.951056516295154)
      parameter(QRT5  = 0.559016994374947)
      parameter(SIN60 = 0.866025403784437)
!
      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      KSTOP=(N-IFAC)/(2*IFAC)
!
      IBAD=1
      IF (LOT.GT.64) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO
!
!     CODING FOR FACTOR 2
!     -------------------
200   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      JA=1
      JB=JA+JINK
!
      IF (LA.EQ.M) GO TO 290
!
      DO 220 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 210 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=A(IA+I)-A(IB+I)
          I=I+INC3
          J=J+INC4
210       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
220     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB-IINK
      IBASE=0
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IB) GO TO 260
      DO 250 K=LA,KSTOP,LA
        KB=K+K
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        IBASE=0
        DO 240 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 230 IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            D(JA+J)=B(IA+I)-B(IB+I)
            C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
            D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
            I=I+INC3
            J=J+INC4
230         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
240       CONTINUE
        IA=IA+IINK
        IB=IB-IINK
        JBASE=JBASE+JUMP
250     CONTINUE
      IF (IA.GT.IB) GO TO 900
260   CONTINUE
      IBASE=0
      DO 280 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 270 IJK=1,LOT
          C(JA+J)=A(IA+I)
          C(JB+J)=-B(IA+I)
          I=I+INC3
          J=J+INC4
270       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
280     CONTINUE
      GO TO 900
!
290   CONTINUE
      DO 294 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 292 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+A(IB+I))
          C(JB+J)=2.0*(A(IA+I)-A(IB+I))
          I=I+INC3
          J=J+INC4
292       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
294     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 3
!     -------------------
300   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
!
      IF (LA.EQ.M) GO TO 390
!
      DO 320 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 310 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=(A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))
          C(JC+J)=(A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))
          I=I+INC3
          J=J+INC4
310       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
320     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IC) GO TO 360
      DO 350 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        IBASE=0
        DO 340 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 330 IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
            D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
            C(JB+J)=C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+    &
             I)+B(IC+I))))   -S1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60 &
             *(A(IB+I)-A(IC+I))))
            D(JB+J)=S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+    &
             I)+B(IC+I))))   +C1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+(SIN60 &
             *(A(IB+I)-A(IC+I))))
            C(JC+J)=C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+    &
             I)+B(IC+I))))   -S2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60 &
             *(A(IB+I)-A(IC+I))))
            D(JC+J)=S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+    &
             I)+B(IC+I))))   +C2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-(SIN60 &
             *(A(IB+I)-A(IC+I))))
            I=I+INC3
            J=J+INC4
330         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
340       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC-IINK
        JBASE=JBASE+JUMP
350     CONTINUE
      IF (IA.GT.IC) GO TO 900
360   CONTINUE
      IBASE=0
      DO 380 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 370 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
          C(JC+J)=-(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
          I=I+INC3
          J=J+INC4
370       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
380     CONTINUE
      GO TO 900
!
390   CONTINUE
      SSIN60=2.0*SIN60
      DO 394 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 392 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+A(IB+I))
          C(JB+J)=(2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
          C(JC+J)=(2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
          I=I+INC3
          J=J+INC4
392       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
394     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 4
!     -------------------
400   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
!
      IF (LA.EQ.M) GO TO 490
!
      DO 420 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 410 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
          C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
          C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
          C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
          I=I+INC3
          J=J+INC4
410       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
420     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.IC) GO TO 460
      DO 450 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        IBASE=0
        DO 440 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 430 IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
            D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
            C(JC+J)=C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                   -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
            D(JC+J)=S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))) &
                   +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
            C(JB+J)=C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                   -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
            D(JB+J)=S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I))) &
                   +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
            C(JD+J)=C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                   -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
            D(JD+J)=S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I))) &
                   +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
            I=I+INC3
            J=J+INC4
430         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
440       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC-IINK
        ID=ID-IINK
        JBASE=JBASE+JUMP
450     CONTINUE
      IF (IB.GT.IC) GO TO 900
460   CONTINUE
      IBASE=0
      SIN45=SQRT(0.5)
      DO 480 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 470 IJK=1,LOT
          C(JA+J)=A(IA+I)+A(IB+I)
          C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
          C(JC+J)=B(IB+I)-B(IA+I)
          C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
          I=I+INC3
          J=J+INC4
470       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
480     CONTINUE
      GO TO 900
!
490   CONTINUE
      DO 494 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 492 IJK=1,LOT
          C(JA+J)=2.0*((A(IA+I)+A(IC+I))+A(IB+I))
          C(JB+J)=2.0*((A(IA+I)-A(IC+I))-B(IB+I))
          C(JC+J)=2.0*((A(IA+I)+A(IC+I))-A(IB+I))
          C(JD+J)=2.0*((A(IA+I)-A(IC+I))+B(IB+I))
          I=I+INC3
          J=J+INC4
492       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
494     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 5
!     -------------------
500   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC
      IE=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
!
      IF (LA.EQ.M) GO TO 590
!
      DO 520 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 510 IJK=1,LOT
          C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
          C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I &
           )))     -(SIN72*B(IB+I)+SIN36*B(IC+I))
          C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I &
           )))     -(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*(A(IB+I)-A(IC+I &
           )))     +(SIN36*B(IB+I)-SIN72*B(IC+I))
          C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*(A(IB+I)-A(IC+I &
           )))     +(SIN72*B(IB+I)+SIN36*B(IC+I))
          I=I+INC3
          J=J+INC4
510       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
520     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.ID) GO TO 560
      DO 550 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        IBASE=0
        DO 540 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 530 IJK=1,LOT
!
            A10(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))) &
             )         +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
            A20(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))) &
             )         -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
            B10(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))) &
             )         +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
            B20(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I))) &
             )         -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
            A11(IJK)=SIN72*(B(IB+I)+B(IE+I))+SIN36*(B(IC+I)+B(ID+I))
            A21(IJK)=SIN36*(B(IB+I)+B(IE+I))-SIN72*(B(IC+I)+B(ID+I))
            B11(IJK)=SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))
            B21(IJK)=SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))
!
            C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
            D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
            C(JB+J)=C1*(A10(IJK)-A11(IJK))-S1*(B10(IJK)+B11(IJK))
            D(JB+J)=S1*(A10(IJK)-A11(IJK))+C1*(B10(IJK)+B11(IJK))
            C(JE+J)=C4*(A10(IJK)+A11(IJK))-S4*(B10(IJK)-B11(IJK))
            D(JE+J)=S4*(A10(IJK)+A11(IJK))+C4*(B10(IJK)-B11(IJK))
            C(JC+J)=C2*(A20(IJK)-A21(IJK))-S2*(B20(IJK)+B21(IJK))
            D(JC+J)=S2*(A20(IJK)-A21(IJK))+C2*(B20(IJK)+B21(IJK))
            C(JD+J)=C3*(A20(IJK)+A21(IJK))-S3*(B20(IJK)-B21(IJK))
            D(JD+J)=S3*(A20(IJK)+A21(IJK))+C3*(B20(IJK)-B21(IJK))
!
            I=I+INC3
            J=J+INC4
530         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
540       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC+IINK
        ID=ID-IINK
        IE=IE-IINK
        JBASE=JBASE+JUMP
550     CONTINUE
      IF (IB.GT.ID) GO TO 900
560   CONTINUE
      IBASE=0
      DO 580 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 570 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
          C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+I &
           )))     -(SIN36*B(IA+I)+SIN72*B(IB+I))
          C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+(0.25*(A(IA+I)+A(IB+I))-A(IC+ &
           I)))    -(SIN36*B(IA+I)+SIN72*B(IB+I))
          C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+I &
           )))     -(SIN72*B(IA+I)-SIN36*B(IB+I))
          C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-(0.25*(A(IA+I)+A(IB+I))-A(IC+ &
           I)))    -(SIN72*B(IA+I)-SIN36*B(IB+I))
          I=I+INC3
          J=J+INC4
570       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
580     CONTINUE
      GO TO 900
!
590   CONTINUE
      QQRT5=2.0*QRT5
      SSIN36=2.0*SIN36
      SSIN72=2.0*SIN72
      DO 594 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 592 IJK=1,LOT
          C(JA+J)=2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))
          C(JB+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          C(JC+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JD+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-SSIN72*B(IC+I))
          C(JE+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I))) &
                 +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+SSIN36*B(IC+I))
          I=I+INC3
          J=J+INC4
592       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
594     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 6
!     -------------------
600   CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC+2*M*INC1
      IE=IC
      IF=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
!
      IF (LA.EQ.M) GO TO 690
!
      DO 620 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 610 IJK=1,LOT
          C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
          C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
          C(JB+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
                   -(SIN60*(B(IB+I)+B(IC+I)))
          C(JF+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I))) &
                   +(SIN60*(B(IB+I)+B(IC+I)))
          C(JC+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
                   -(SIN60*(B(IB+I)-B(IC+I)))
          C(JE+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I))) &
                   +(SIN60*(B(IB+I)-B(IC+I)))
          I=I+INC3
          J=J+INC4
610       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
620     CONTINUE
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IC.EQ.ID) GO TO 660
      DO 650 K=LA,KSTOP,LA
        KB=K+K
        KC=KB+KB
        KD=KC+KB
        KE=KD+KB
        KF=KE+KB
        C1=TRIGS(KB+1)
        S1=TRIGS(KB+2)
        C2=TRIGS(KC+1)
        S2=TRIGS(KC+2)
        C3=TRIGS(KD+1)
        S3=TRIGS(KD+2)
        C4=TRIGS(KE+1)
        S4=TRIGS(KE+2)
        C5=TRIGS(KF+1)
        S5=TRIGS(KF+2)
        IBASE=0
        DO 640 L=1,LA
          I=IBASE
          J=JBASE
!DIR$ IVDEP
          DO 630 IJK=1,LOT
!
            A11(IJK)= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
            A20(IJK)=(A(IA+I)+A(ID+I))-0.5*A11(IJK)
            A21(IJK)=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
            B11(IJK)= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
            B20(IJK)=(B(IA+I)-B(ID+I))-0.5*B11(IJK)
            B21(IJK)=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))
!
            C(JA+J)=(A(IA+I)+A(ID+I))+A11(IJK)
            D(JA+J)=(B(IA+I)-B(ID+I))+B11(IJK)
            C(JC+J)=C2*(A20(IJK)-B21(IJK))-S2*(B20(IJK)+A21(IJK))
            D(JC+J)=S2*(A20(IJK)-B21(IJK))+C2*(B20(IJK)+A21(IJK))
            C(JE+J)=C4*(A20(IJK)+B21(IJK))-S4*(B20(IJK)-A21(IJK))
            D(JE+J)=S4*(A20(IJK)+B21(IJK))+C4*(B20(IJK)-A21(IJK))
!
            A11(IJK)=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
            B11(IJK)=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
            A20(IJK)=(A(IA+I)-A(ID+I))-0.5*A11(IJK)
            A21(IJK)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            B20(IJK)=(B(IA+I)+B(ID+I))+0.5*B11(IJK)
            B21(IJK)=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))
!
            C(JD+J)=C3*((A(IA+I)-A(ID+I))+A11(IJK))-S3*((B(IA+I)+B(ID+I &
             ))-B11(IJK))
            D(JD+J)=S3*((A(IA+I)-A(ID+I))+A11(IJK))+C3*((B(IA+I)+B(ID+I &
             ))-B11(IJK))
            C(JB+J)=C1*(A20(IJK)-B21(IJK))-S1*(B20(IJK)-A21(IJK))
            D(JB+J)=S1*(A20(IJK)-B21(IJK))+C1*(B20(IJK)-A21(IJK))
            C(JF+J)=C5*(A20(IJK)+B21(IJK))-S5*(B20(IJK)+A21(IJK))
            D(JF+J)=S5*(A20(IJK)+B21(IJK))+C5*(B20(IJK)+A21(IJK))
!
            I=I+INC3
            J=J+INC4
630         CONTINUE
          IBASE=IBASE+INC1
          JBASE=JBASE+INC2
640       CONTINUE
        IA=IA+IINK
        IB=IB+IINK
        IC=IC+IINK
        ID=ID-IINK
        IE=IE-IINK
        IF=IF-IINK
        JBASE=JBASE+JUMP
650     CONTINUE
      IF (IC.GT.ID) GO TO 900
660   CONTINUE
      IBASE=0
      DO 680 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 670 IJK=1,LOT
          C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
          C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
          C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
          C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-(0.5*(B(IA+I)+B(IC+I))+B(IB+I))
          C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
          C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-(0.5*(A(IA+I)+A(IC+I))-A(IB+I))
          I=I+INC3
          J=J+INC4
670       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
680     CONTINUE
      GO TO 900
!
690   CONTINUE
      SSIN60=2.0*SIN60
      DO 694 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 692 IJK=1,LOT
          C(JA+J)=(2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))
          C(JD+J)=(2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))
          C(JB+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
                   -(SSIN60*(B(IB+I)+B(IC+I)))
          C(JF+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I))) &
                   +(SSIN60*(B(IB+I)+B(IC+I)))
          C(JC+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
                   -(SSIN60*(B(IB+I)-B(IC+I)))
          C(JE+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I))) &
                   +(SSIN60*(B(IB+I)-B(IC+I)))
          I=I+INC3
          J=J+INC4
692       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
694     CONTINUE
      GO TO 900
!
!     CODING FOR FACTOR 8
!     -------------------
800   CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+LA*INC1
      IC=IB+2*LA*INC1
      ID=IC+2*LA*INC1
      IE=ID+2*LA*INC1
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
      JG=JF+JINK
      JH=JG+JINK
      SSIN45=SQRT(2.0)
!
      DO 820 L=1,LA
        I=IBASE
        J=JBASE
!DIR$ IVDEP
        DO 810 IJK=1,LOT
          C(JA+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+(A(IB+I)+A(ID+I)))
          C(JE+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-(A(IB+I)+A(ID+I)))
          C(JC+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-(B(IB+I)-B(ID+I)))
          C(JG+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+(B(IB+I)-B(ID+I)))
          C(JB+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
                   +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JF+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I)) &
                   -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
          C(JD+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
                   -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          C(JH+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I)) &
                   +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
          I=I+INC3
          J=J+INC4
810       CONTINUE
        IBASE=IBASE+INC1
        JBASE=JBASE+INC2
820     CONTINUE
!
!     RETURN
!     ------
900   CONTINUE
      IBAD=0
910   CONTINUE
      IERR=IBAD
      RETURN
      END



       
    