program shtns_concurrency
   ! Can several threads call SHTns at once, on one config?
   !
   ! That is the single assumption the chosen parallel design rests on. The
   ! model transforms about five fields at ten levels a phase, so there are
   ! roughly fifty independent full-globe transforms to spread over sixteen
   ! threads -- IF a plan can be shared. SHTns says so only by implication:
   ! SHT_ALLOW_GPU documents that enabling the GPU means "the same plan cannot
   ! be used simultaneously by different threads ANYMORE", and the thread-safety
   ! warning in the changelog is scoped to MKL, which this build does not use.
   ! Implication is not measurement.
   !
   ! The check has a right answer rather than a comparison: each thread
   ! transforms a spectrum only it knows, and every result must equal what the
   ! same transform produces alone. A race would show as a wrong answer on some
   ! thread, and running many rounds is what makes an intermittent one likely to
   ! appear at least once.
   use iso_c_binding
   use iso_fortran_env, only: real64
   use omp_lib
   implicit none
   include 'shtns.f03'

   integer, parameter :: dp = real64
   integer, parameter :: NF = 50          ! independent transforms, as the model has
   integer, parameter :: NROUND = 20

   type(shtns_info), pointer :: shtns
   type(c_ptr) :: shtns_c
   integer :: lmax, mmax, mres, nlat, nphi, norm, layout, nthreads
   integer :: j, jf, jr, nbad
   real(dp) :: eps_polar, worst, e
   complex(dp), allocatable :: Slm(:,:)
   real(dp), allocatable :: Sh(:,:,:), Ref(:,:,:)
!  the analysis direction, which this probe did not originally cover
   complex(dp), allocatable :: Alm(:,:), Aref(:,:), Blm(:,:), Bref(:,:)
   real(dp), allocatable :: Vt(:,:,:), Vp(:,:,:)
   integer :: nbada, nbadv

   lmax = 170 ; mmax = 170 ; mres = 1
   nlat = 256 ; nphi = 512
   eps_polar = 0.0_dp
   norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE
   layout = SHT_GAUSS + SHT_PHI_CONTIGUOUS

   call shtns_verbose(0)
   nthreads = shtns_use_threads(1)        ! SHTns itself single threaded: the
                                          ! model's team provides the parallelism
   shtns_c = shtns_create(lmax, mmax, mres, norm)
   call shtns_set_grid(shtns_c, layout, eps_polar, nlat, nphi)
   call c_f_pointer(cptr=shtns_c, fptr=shtns)

   allocate( Slm(shtns%nlm, NF) )
   allocate( Sh(shtns%nphi, shtns%nlat, NF), Ref(shtns%nphi, shtns%nlat, NF) )

   ! a different spectrum per field, so a thread that writes another's output
   ! produces a visibly wrong answer rather than a plausible one
   do jf = 1, NF
      do j = 1, shtns%nlm
         Slm(j,jf) = cmplx(real(jf,dp)/real(j,dp), -0.5_dp/real(j*jf,dp), kind=dp)
      enddo
   enddo

   ! the reference, computed one at a time on one thread
   do jf = 1, NF
      call SH_to_spat(shtns_c, Slm(:,jf), Ref(:,:,jf))
   enddo

   nbad = 0 ; worst = 0.0_dp
   do jr = 1, NROUND
      Sh = 0.0_dp
      !$omp parallel do schedule(dynamic) private(jf) shared(Slm, Sh, shtns_c)
      do jf = 1, NF
         call SH_to_spat(shtns_c, Slm(:,jf), Sh(:,:,jf))
      enddo
      !$omp end parallel do
      do jf = 1, NF
         e = maxval(abs(Sh(:,:,jf) - Ref(:,:,jf)))
         if (e > worst) worst = e
         if (e /= 0.0_dp) nbad = nbad + 1
      enddo
   enddo

!  ---- ANALYSIS, which is a different question and was never asked ----------
!  Synthesis reads the configuration and writes only the caller's array. There
!  is no reason analysis has to be the same, and the model now calls
!  spat_to_SH and spat_to_SHsphtor from every thread at once. A shared scratch
!  inside the configuration would not make the answer wrong -- it would make it
!  differ in the last bits from run to run, which is what a threaded model must
!  never do.
   allocate( Alm(shtns%nlm, NF), Aref(shtns%nlm, NF) )
   allocate( Blm(shtns%nlm, NF), Bref(shtns%nlm, NF) )
   allocate( Vt(shtns%nphi, shtns%nlat, NF), Vp(shtns%nphi, shtns%nlat, NF) )
   do jf = 1, NF
      Vt(:,:,jf) = Ref(:,:,jf)
      Vp(:,:,jf) = Ref(:,:,jf) * 0.5_dp
   enddo
   do jf = 1, NF
      call spat_to_SH(shtns_c, Ref(:,:,jf), Aref(:,jf))
   enddo
   do jf = 1, NF
      call spat_to_SHsphtor(shtns_c, Vt(:,:,jf), Vp(:,:,jf), Bref(:,jf), Blm(:,jf))
   enddo

   nbada = 0 ; nbadv = 0
   do jr = 1, NROUND
      Alm = (0.0_dp, 0.0_dp)
      !$omp parallel do schedule(dynamic) private(jf) shared(Ref, Alm, shtns_c)
      do jf = 1, NF
         call spat_to_SH(shtns_c, Ref(:,:,jf), Alm(:,jf))
      enddo
      !$omp end parallel do
      do jf = 1, NF
         if (any(Alm(:,jf) /= Aref(:,jf))) nbada = nbada + 1
      enddo
   enddo

   write(*,'(a,i0,a,i0,a,i0)') 'lmax ', lmax, '  nlat ', nlat, '  nlm ', shtns%nlm
   write(*,'(a,i0,a,i0,a)') 'threads available: ', omp_get_max_threads(),        &
  &   ', transforms a round: ', NF, ''
   write(*,'(a,i0)')        'rounds:            ', NROUND
   write(*,*)
   if (nbad == 0) then
      write(*,'(a,i0,a)') '  [  ok  ] all ', NF*NROUND,                          &
  &      ' concurrent transforms are BIT IDENTICAL to the serial answer'
   else
      write(*,'(a,i0,a,i0,a)') '  [ FAIL ] ', nbad, ' of ', NF*NROUND,           &
  &      ' concurrent transforms differ from the serial answer'
      write(*,'(a,e12.5)')     '           worst difference ', worst
      stop 1
   endif
   if (nbada == 0) then
      write(*,'(a,i0,a)') '  [  ok  ] all ', NF*NROUND,                          &
  &      ' concurrent spat_to_SH are BIT IDENTICAL to the serial answer'
   else
      write(*,'(a,i0,a,i0,a)') '  [ FAIL ] ', nbada, ' of ', NF*NROUND,          &
  &      ' concurrent spat_to_SH differ from the serial answer'
      stop 1
   endif
end program shtns_concurrency
