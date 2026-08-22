program vecconventions
   ! What relates SHTns's vector transform to this model's dv2uv?
   !
   ! This is the one that can be wrong invisibly. SHTns's spheroidal and
   ! toroidal potentials differ from divergence and vorticity by l(l+1), and
   ! this model already carries 1/(n(n+1)) inside fmu and fmv, so a substitution
   ! that assumes the libraries agree applies the factor twice. The result runs,
   ! and looks like weather.
   !
   ! Three things are measured at once and none is assumed:
   !   the SCALE, and whether it depends on l the way l(l+1) would;
   !   the SIGN of the meridional component, since theta grows southward while
   !     latitude grows northward;
   !   whether Robert form matches, gu being "zonal wind (*cos(phi))".
   !
   ! plavor is forced to zero. dv2uv adds planetary vorticity to mode w=2, which
   ! would land on exactly the modes being compared and read as a convention.
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env, only: real64
   use pumamod, only: NLAT, NLON, NLPP, NLEV, NTRU, NCSP, NESP, sid, gwd, plavor
   implicit none
   include 'shtns.f03'

   integer, parameter :: wp = real64
   type(shtns_info), pointer :: shc
   type(c_ptr) :: shp
   integer :: knorm, klay, klmax, jj, jm, jn, jlm, kshown
   real(wp) :: zeps, zsi(NLAT), zgw(NLAT)
   real :: zsd(2,NESP/2,NLEV), zsz(2,NESP/2,NLEV)
   real :: zfu(2,NLON/2,NLPP,NLEV), zfv(2,NLON/2,NLPP,NLEV)
   real(wp), allocatable :: zvt(:,:), zvp(:,:)
   complex(wp), allocatable :: zs(:), zt(:)
   real(wp) :: ru, rv, expect

   call inigau(NLAT,zsi,zgw)
   do jj = 1 , NLPP
      sid(jj) = zsi(jj)
      gwd(jj) = zgw(jj)
   enddo
   plavor = 0.0                     ! see the header: it would land on mode 2
   call legini

   klmax = NTRU ; zeps = 0.0_wp
   knorm = SHT_ORTHONORMAL          ! Condon-Shortley INCLUDED, as PlaSim has it
   klay  = SHT_GAUSS + SHT_PHI_CONTIGUOUS
   call shtns_verbose(0)
   jj = shtns_use_threads(1)
   shp = shtns_create(klmax, klmax, 1, knorm)
   call shtns_set_grid(shp, klay, zeps, NLAT, NLON)
   call c_f_pointer(cptr=shp, fptr=shc)
   call shtns_robert_form(shp, 1)   ! gu is the wind times cos(phi)
   allocate( zvt(shc%nphi, shc%nlat), zvp(shc%nphi, shc%nlat) )
   allocate( zs(shc%nlm), zt(shc%nlm) )

   write(*,'(a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT
   write(*,'(a)') 'Robert form ON. ratios are model/SHTns on the matching component.'
   write(*,*)
   write(*,'(a4,a4,a8,a14,a14,a16)') 'm','n','lm','u ratio','v ratio','ratio*l(l+1)'

   kshown = 0
   jlm = 0
   do jm = 0 , NTRU
      do jn = jm , NTRU
         jlm = jlm + 1
         if (jn == 0) cycle                       ! l=0 has no vector part
         if (jm > 3) cycle
         if (jn > jm+1) cycle
         call pair(jlm, ru, rv)
         expect = ru * real(jn*(jn+1), wp)
         if (kshown < 12) then
            write(*,'(i4,i4,i8,f14.7,f14.7,f16.7)') jm, jn, jlm, ru, rv, expect
            kshown = kshown + 1
         endif
      enddo
   enddo

contains

   subroutine pair(klm, pu, pv)
   integer, intent(in) :: klm
   real(wp), intent(out) :: pu, pv
   real(wp) :: au, av, bu, bv
   ! divergence only, on level 1
   zsd(:,:,:) = 0.0 ; zsz(:,:,:) = 0.0
   zsd(1,klm,1) = 1.0
   call dv2uv(zsd,zsz,zfu,zfv)
   call fc2gp(zfu(1,1,1,1),NLON,NLPP)
   call fc2gp(zfv(1,1,1,1),NLON,NLPP)
   ! SIGNED, by least squares against the SHTns field rather than by comparing
   ! maxima. Magnitudes alone cannot see a flipped meridional component, and
   ! theta grows southward while latitude grows northward, so a sign difference
   ! is expected rather than hypothetical.

   zs(:) = (0.0_wp,0.0_wp) ; zt(:) = (0.0_wp,0.0_wp)
   zs(klm) = (1.0_wp,0.0_wp)
   call SHsphtor_to_spat(shp, zs, zt, zvt, zvp)
   pu = lstsq(reshape(real(zfu(:,:,:,1),wp),[NLON*NLPP]), reshape(zvp,[NLON*NLPP]))
   pv = lstsq(reshape(real(zfv(:,:,:,1),wp),[NLON*NLPP]), reshape(zvt,[NLON*NLPP]))
   au = 0.0_wp ; av = 0.0_wp ; bu = 0.0_wp ; bv = 0.0_wp
   end subroutine pair

   ! the scale s minimising |a - s*b|, which carries the sign a ratio of maxima
   ! throws away
   function lstsq(a, b) result(s)
   real(wp), intent(in) :: a(:), b(:)
   real(wp) :: s, den
   den = sum(b*b)
   if (den > 0.0_wp) then
      s = sum(a*b) / den
   else
      s = 0.0_wp
   endif
   end function lstsq

end program vecconventions
