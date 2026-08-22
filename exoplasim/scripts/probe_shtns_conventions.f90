program conventions
   ! What factor relates SHTns's transform to this model's, mode by mode?
   !
   ! Reading two normalisation conventions against each other and hoping is how
   ! a factor of sqrt(2) survives into a climate. This drives ONE spectral mode
   ! through both transforms onto the same grid and divides. If the conventions
   ! agree the ratio is 1 everywhere; if they differ by a constant it is that
   ! constant; if it depends on m or l, the dependence IS the convention and is
   ! visible rather than inferred.
   !
   ! The index map is already known and is the easy half: SHTns's LM(l,m) is
   ! 0-based with m outer and l inner, which is the order legini builds lm in,
   ! so PlaSim lm maps to SHTns lm-1. Only the scaling is in question.
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env, only: real64
   use pumamod, only: NLAT, NLON, NLPP, NTRU, NCSP, sid, gwd
   implicit none
   include 'shtns.f03'

   integer, parameter :: dp = real64
   type(shtns_info), pointer :: sh
   type(c_ptr) :: sh_c
   integer :: norm, layout, lmax, mmax, mres, j, m, n, lm, lmc, nshown
   real(dp) :: eps_polar
   real(dp) :: zsi(NLAT), zgw(NLAT)
   real :: zspc(2,NCSP), zgpf(NLON,NLPP)
   real(dp), allocatable :: Sh_g(:,:)
   complex(dp), allocatable :: Slm(:)
   real(dp) :: rmax, gmax, ratio, worstdev
   character(len=8) :: tag

   ! --- the model side, one process so the globe is local -------------------
   call inigau(NLAT,zsi,zgw)
   do j = 1 , NLPP
      sid(j) = zsi(j)
      gwd(j) = zgw(j)
   enddo
   call legini

   ! --- SHTns on the same grid and truncation -------------------------------
   lmax = NTRU ; mmax = NTRU ; mres = 1
   eps_polar = 0.0_dp
   norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
   layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS
   call shtns_verbose(0)
   j = shtns_use_threads(1)
   sh_c = shtns_create(lmax, mmax, mres, norm)
   call shtns_set_grid(sh_c, layout, eps_polar, NLAT, NLON)
   call c_f_pointer(cptr=sh_c, fptr=sh)

   if (sh%nlm /= NCSP) then
      write(*,*) 'mode counts differ:', sh%nlm, 'against NCSP', NCSP
      stop 2
   endif
   allocate( Sh_g(sh%nphi, sh%nlat), Slm(sh%nlm) )

   write(*,'(a,i0,a,i0,a,i0)') 'T', NTRU, '  NLAT ', NLAT, '  modes ', NCSP
   write(*,'(a)') 'ratio = (model grid) / (SHTns grid), same single mode'
   write(*,*)
   write(*,'(a4,a4,a8,a16,a16,a14)') 'm','n','lm','model max','shtns max','ratio'

   nshown = 0
   worstdev = 0.0_dp
   lm = 0
   do m = 0 , NTRU
      do n = m , NTRU
         lm = lm + 1
         ! only a readable sample, but the deviation below is over ALL of them
         if (.not. (m <= 2 .or. m == NTRU/2) .or. n > m+1) then
            call oneratio(lm, m, n, ratio, rmax, gmax)
            if (rmax > 1.0e-8_dp) worstdev = max(worstdev, abs(ratio-1.0_dp))
            cycle
         endif
         call oneratio(lm, m, n, ratio, rmax, gmax)
         if (rmax > 1.0e-8_dp) then
            worstdev = max(worstdev, abs(ratio-1.0_dp))
            if (nshown < 14) then
               write(*,'(i4,i4,i8,e16.6,e16.6,f14.8)') m, n, lm, rmax, gmax, ratio
               nshown = nshown + 1
            endif
         endif
      enddo
   enddo

   write(*,*)
   write(*,'(a,f14.8)') 'worst |ratio - 1| over every mode: ', worstdev
   if (worstdev < 1.0e-10_dp) then
      write(*,'(a)') '  the conventions AGREE as configured; no rescaling needed'
   else
      write(*,'(a)') '  the conventions DIFFER; the table above shows how'
   endif

contains

   subroutine oneratio(klm, km, kn, pratio, pmodel, pshtns)
   integer, intent(in) :: klm, km, kn
   real(dp), intent(out) :: pratio, pmodel, pshtns
   integer :: i
   zspc(:,:) = 0.0
   zspc(1,klm) = 1.0
   call sp2fc(zspc,zgpf)
   call fc2gp(zgpf,NLON,NLPP)
   Slm(:) = (0.0_dp, 0.0_dp)
   Slm(klm) = (1.0_dp, 0.0_dp)          ! PlaSim lm  <->  SHTns lm-1, 1-based here
   call SH_to_spat(sh_c, Slm, Sh_g)
   pmodel = maxval(abs(real(zgpf,dp)))
   pshtns = maxval(abs(Sh_g))
   if (pshtns > 0.0_dp) then
      pratio = pmodel / pshtns
   else
      pratio = 0.0_dp
   endif
   end subroutine oneratio

end program conventions
