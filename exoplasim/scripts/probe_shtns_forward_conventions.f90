      program fwdconventions
!     What relates SHTns's ANALYSIS to this model's forward transforms?
!
!     The inverse direction is settled and verified; this is its adjoint, and it
!     is not simply the reciprocal. Three things differ on the way back:
!
!       the FILTER is a different filter. legini folds skspgp into fsp, fmu and
!         fmv for spectral-to-grid, and skgpsp into fgp, fmm and fmq for
!         grid-to-spectral. They are set by separate namelist switches
!         (nspvfilter, ngptfilter) and are not the same array;
!       the QUADRATURE WEIGHT is explicit in legmod -- gwd in fc2sp, gwdc =
!         gwd/cos^2 in uv2dv and in mktend -- and implicit in SHTns, which
!         applies its own. So the constant here is not the inverse's;
!       uv2dv produces DIVERGENCE and VORTICITY, while spat_to_SHsphtor produces
!         spheroidal and toroidal POTENTIALS, which differ by l(l+1) -- the same
!         factor the inverse direction carries, and the same chance to apply it
!         twice.
!
!     THE TEST FIELD IS BAND LIMITED BY CONSTRUCTION. A forward transform is a
!     quadrature, so an arbitrary grid field carries power above the truncation
!     that the two sides alias differently, and the ratio would then measure the
!     aliasing rather than the convention. So the grid field is made by running
!     a dense SPECTRAL field through the inverse path, which is already verified.
!
!     SIGNED LEAST SQUARES, per mode, as for the inverse. Magnitudes hid the
!     Condon-Shortley phase once and a per-mode ratio that drove only one field
!     hid the toroidal sign once. Both mistakes are cheap to repeat here.
!
!     THE FILTER IS ON AND THE PLANET ROTATES, because a probe that switches
!     them off measures a configuration the model never runs. That is
!     docs/src/practice/failure-modes.md class 29, and it cost a session.
      use, intrinsic :: iso_c_binding
      use, intrinsic :: iso_fortran_env, only: real64
      use pumamod, only: NLAT, NLON, NLPP, NLEV, NTRU, NTP1, NCSP, NESP,  &
     &                   NUGP, sid, gwd, plavor, EZ, nfilter, ngptfilter, &
     &                   nspvfilter, filterkappa, nfilterexp
      use legmod, only: fgp, fsp
      use shtnsmod, only: shtns_setup, shtcfg, SHTROOT
      implicit none
      include 'shtns.f03'

      integer, parameter :: wp = real64
      integer :: jj, jm, jn, jlm, jlat, jlon, kshown
      real(wp) :: zsi(NLAT), zgw(NLAT)
      real :: zspc(2,NCSP), zgpf(NLON,NLPP), zwork(NLON,NLPP)
      real :: zsl(2,NESP/2)
      real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
      real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
      real :: zdl(2,NESP/2,NLEV), zzl(2,NESP/2,NLEV)
      real :: zgu(NLON,NLPP), zgv(NLON,NLPP)
      real(wp) :: zg(NUGP), zvt(NUGP), zvp(NUGP)
      real :: zflat(NUGP)
      complex(wp) :: zlm(NCSP), zs(NCSP), zt(NCSP)
      real(wp) :: rsc, rdv, rvo, znn1

      call inigau(NLAT,zsi,zgw)
      do jj = 1 , NLPP
         sid(jj) = zsi(jj)
         gwd(jj) = zgw(jj)
      enddo
      nfilter     = 2
      ngptfilter  = 1
      nspvfilter  = 1
      filterkappa = 8.0
      nfilterexp  = 8
      plavor      = EZ
      call legini
      call shtns_setup

      write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT,            &
     &                            '  modes ', NCSP
      write(*,'(a)') 'ratio = legmod coefficient / SHTns coefficient,'
      write(*,'(a)') 'then divided by the forward filter fgp, which legmod'
      write(*,'(a)') 'carries and SHTns does not. A CONSTANT column is the'
      write(*,'(a)') 'answer; a column that drifts with n means a factor of n'
      write(*,'(a)') 'is still unaccounted for. 1/sqrt(2 pi) is 0.39894228.'
      write(*,*)

!     ---- a band-limited grid field, built through the verified inverse -------
      do jj = 1 , NCSP
         zspc(1,jj) =  1.0 / real(jj)
         zspc(2,jj) = -0.5 / real(jj)
      enddo
      zspc(2,1:NTP1) = 0.0
      call sp2fc(zspc,zgpf)
      call fc2gp(zgpf,NLON,NLPP)

!     ---- scalar: gp2fc then fc2sp, against spat_to_SH ------------------------
      zwork = zgpf
      call gp2fc(zwork,NLON,NLPP)
      call fc2sp(zwork,zsl)

      zflat = reshape(zgpf,[NUGP])
      zg    = real(zflat, wp)
      call spat_to_SH(shtcfg, zg, zlm)

      write(*,'(a)') 'SCALAR, fc2sp against spat_to_SH'
      write(*,'(a4,a4,a8,a18,a18)') 'm','n','lm','ratio','/ fgp'
      kshown = 0
      jlm = 0
      do jm = 0 , NTRU
         do jn = jm , NTRU
            jlm = jlm + 1
            if (jm > 2 .or. jn > jm+2) cycle
            if (kshown >= 12) cycle
            if (abs(real(zlm(jlm))) < 1.0e-30_wp) cycle
            rsc = real(zsl(1,jlm),wp) / real(zlm(jlm))
            write(*,'(i4,i4,i8,f18.9,f18.9)') jm, jn, jlm, rsc,          &
     &                                        rsc / real(fgp(jlm),wp)
            kshown = kshown + 1
         enddo
      enddo

!     ---- vector: uv2dv against spat_to_SHsphtor ------------------------------
!     Winds made from a dense divergence AND a dense vorticity together, because
!     driving one alone leaves the other's sign untested -- which is exactly how
!     the toroidal sign survived the inverse probe.
      zsd(:,:,:) = 0.0 ; zsz(:,:,:) = 0.0
      do jj = 1 , NCSP
         zsd(1,jj,1) =  1.0 / real(jj)
         zsd(2,jj,1) = -0.5 / real(jj)
         zsz(1,jj,1) =  0.7 / real(jj)
         zsz(2,jj,1) =  0.3 / real(jj)
      enddo
      zsd(2,1:NTP1,1) = 0.0 ; zsz(2,1:NTP1,1) = 0.0
      plavor = 0.0            ! dv2uv would otherwise seed mode 2 of the wind
      call dv2uv(zsd,zsz,zfu,zfv)
      call fc2gp(zfu(1,1,1,1),NLON,NLPP)
      call fc2gp(zfv(1,1,1,1),NLON,NLPP)
      zgu = reshape(zfu(:,:,:,1),[NLON,NLPP])
      zgv = reshape(zfv(:,:,:,1),[NLON,NLPP])
      plavor = EZ

      zfu(:,:,:,1) = reshape(zgu,[2,NLON/2,NLPP])
      zfv(:,:,:,1) = reshape(zgv,[2,NLON/2,NLPP])
      call gp2fc(zfu(1,1,1,1),NLON,NLPP)
      call gp2fc(zfv(1,1,1,1),NLON,NLPP)
      call uv2dv(zfu,zfv,zdl,zzl)

!     SHTns takes the Robert-form components: theta grows southward, so the
!     meridional one is negated on the way in exactly as sh_dv2uv negates it on
!     the way out.
      zflat = reshape(zgv,[NUGP]) ; zvt = -real(zflat, wp)
      zflat = reshape(zgu,[NUGP]) ; zvp =  real(zflat, wp)
      call spat_to_SHsphtor(shtcfg, zvt, zvp, zs, zt)

      write(*,*)
      write(*,'(a)') 'VECTOR, uv2dv against spat_to_SHsphtor, both fields on'
      write(*,'(a4,a4,a8,a18,a18)') 'm','n','lm','div/(S*l(l+1))','vor/(T*l(l+1))'
      kshown = 0
      jlm = 0
      do jm = 0 , NTRU
         do jn = jm , NTRU
            jlm = jlm + 1
            if (jn == 0) cycle
            if (jm > 2 .or. jn > jm+2) cycle
            if (kshown >= 12) cycle
            znn1 = real(jn*(jn+1), wp)
            rdv = 0.0_wp ; rvo = 0.0_wp
            if (abs(real(zs(jlm))) > 1.0e-30_wp)                          &
     &         rdv = real(zdl(1,jlm,1),wp) / (real(zs(jlm)) * znn1        &
     &                                        * real(fgp(jlm),wp))
            if (abs(real(zt(jlm))) > 1.0e-30_wp)                          &
     &         rvo = real(zzl(1,jlm,1),wp) / (real(zt(jlm)) * znn1        &
     &                                        * real(fgp(jlm),wp))
            write(*,'(i4,i4,i8,f18.9,f18.9)') jm, jn, jlm, rdv, rvo
            kshown = kshown + 1
         enddo
      enddo

      end program fwdconventions
