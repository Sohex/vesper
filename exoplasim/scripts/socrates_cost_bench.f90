! Cost bench for the SOCRATES band-resolved scheme, for CLIM-61.
!
! Worldbuilding frame: a COMPUTE measurement of a candidate radiation code
! against the Vesper project's climate model on this desktop. Nothing here is
! about the simulated planet.
!
! Measures wall clock and retired work per COLUMN per RADIATION CALL for the
! runes interface at a configurable layer count and column count, so the number
! is directly comparable with radmod's swr/lwr cost per column per call.
!
! Environment:
!   SB_NLAYER    layers                       (default 10, the project's NLEV)
!   SB_NPROFILE  columns per runes call       (default 2048, T21's NHOR)
!   SB_NREP      calls, after a discarded warm-up
!   SB_SPECTRUM  spectral file stem directory (default .)
!   SB_SW        sw spectral file             (default sp_sw_ga7)
!   SB_LW        lw spectral file             (default sp_lw_ga7)
!   SB_CLOUD     0 clear-sky, 1 cloudy        (default 0)
!   SB_ZENITH    0 a full diurnal spread, so half the columns are dark and the
!                  scheme's own handling of a dark column is in the number;
!                1 every column lit at cos(zenith) = 0.5, which is the arm that
!                  cannot be flattered by a skipped column. Report both: the
!                  broadband scheme's `where` blocks evaluate both branches and
!                  select, so it does not save on a dark column either, and
!                  which arm is fair depends on whether the candidate skips.
!   SB_SETUPONLY 1 to do setup and exit       (the startup control)
program socrates_bench

  use socrates_set_spectrum, only: set_spectrum
  use socrates_runes, only: runes, StrDiag, &
                            ip_source_illuminate, ip_source_thermal, &
                            ip_cloud_representation_ice_water, &
                            ip_cloud_representation_off, &
                            ip_overlap_max_random, ip_inhom_homogeneous
  use realtype_rd, only: RealExt

  implicit none

  ! Vesper's own planet constants, config/planet.yaml. Declared here rather
  ! than defaulted because SOCRATES exposes them and the point of the candidate
  ! is that it does.
  real(RealExt), parameter :: grav_acc   = 11.18
  real(RealExt), parameter :: r_gas_dry  = 287.026
  real(RealExt), parameter :: cp_air_dry = 1.005e+03
  real(RealExt), parameter :: pi = 4.0*atan(1.0)

  integer :: n_profile, n_layer, n_rep, i_cloud, setup_only, i_zenith
  ! The parametrisation indices the ga7 spectral files themselves declare:
  ! block 10 gives droplet TYPE 5 and block 12 ice TYPE 8, which is what these
  ! index -- the type number in the file and not the parametrisation scheme's own
  ! index, which happens to be 5 for the droplets and 12 for the ice. A file
  ! carrying different types fails loudly rather than being evaluated wrong.
  integer :: i_param_water, i_param_ice
  character(len=256) :: sw_file, lw_file, buf

  type(StrDiag) :: sw_diag, lw_diag
  real(RealExt), allocatable, target :: sw_hr(:,:), lw_hr(:,:)
  real(RealExt), allocatable, target :: sw_up(:,:), sw_dn(:,:)
  real(RealExt), allocatable, target :: lw_up(:,:), lw_dn(:,:)

  real(RealExt), allocatable :: p_layer(:,:), t_layer(:,:)
  real(RealExt), allocatable :: p_level(:,:), t_level(:,:)
  real(RealExt), allocatable :: h2o(:,:), o3(:,:)
  real(RealExt), allocatable :: d_mass(:,:), density(:,:), lhc(:,:)
  real(RealExt), allocatable :: t_ground(:), cosz(:), sirr(:)
  real(RealExt), allocatable :: cld_frac(:,:), liq_frac(:,:), ice_frac(:,:)
  real(RealExt), allocatable :: liq_mmr(:,:), ice_mmr(:,:)
  real(RealExt), allocatable :: liq_dim(:,:), ice_dim(:,:)

  ! McClatchey mid-latitude summer, 32 layers, from SOCRATES's own runes_driver.
  ! Regridded below onto n_layer sigma levels by linear interpolation in log p.
  integer, parameter :: nm = 32
  real(RealExt) :: pm(nm) = (/ &
    0.337000E+01, 0.509050E+02, 0.135550E+03, 0.254500E+03, &
    0.492500E+03, 0.986000E+03, 0.204500E+04, 0.299500E+04, &
    0.349000E+04, 0.406500E+04, 0.473500E+04, 0.552500E+04, &
    0.645000E+04, 0.753500E+04, 0.881000E+04, 0.103000E+05, &
    0.120500E+05, 0.141500E+05, 0.166000E+05, 0.194000E+05, &
    0.226000E+05, 0.262000E+05, 0.302500E+05, 0.348000E+05, &
    0.399000E+05, 0.456500E+05, 0.520500E+05, 0.591000E+05, &
    0.669000E+05, 0.756000E+05, 0.852000E+05, 0.957500E+05 /)
  real(RealExt) :: tm(nm) = (/ &
    0.216982E+03, 0.262328E+03, 0.272545E+03, 0.263059E+03, &
    0.250428E+03, 0.238550E+03, 0.228094E+03, 0.223481E+03, &
    0.222481E+03, 0.220962E+03, 0.219481E+03, 0.218481E+03, &
    0.217481E+03, 0.216481E+03, 0.216000E+03, 0.216000E+03, &
    0.216000E+03, 0.216000E+03, 0.216000E+03, 0.219116E+03, &
    0.225632E+03, 0.232109E+03, 0.238624E+03, 0.245104E+03, &
    0.251619E+03, 0.258100E+03, 0.264097E+03, 0.270094E+03, &
    0.276092E+03, 0.282091E+03, 0.287573E+03, 0.292058E+03 /)
  real(RealExt) :: qm(nm) = (/ &
    0.399688E-05, 0.399530E-05, 0.399851E-05, 0.399700E-05, &
    0.399963E-05, 0.400241E-05, 0.400722E-05, 0.400994E-05, &
    0.400705E-05, 0.400353E-05, 0.399929E-05, 0.399791E-05, &
    0.399939E-05, 0.400000E-05, 0.400000E-05, 0.400058E-05, &
    0.400152E-05, 0.402072E-05, 0.485647E-05, 0.109264E-04, &
    0.349482E-04, 0.974304E-04, 0.199405E-03, 0.321272E-03, &
    0.509681E-03, 0.777969E-03, 0.114820E-02, 0.182544E-02, &
    0.305008E-02, 0.485372E-02, 0.722366E-02, 0.101064E-01 /)
  real(RealExt) :: om(nm) = (/ &
    0.606562E-06, 0.252165E-05, 0.469047E-05, 0.748127E-05, &
    0.957770E-05, 0.100812E-04, 0.814088E-05, 0.664711E-05, &
    0.603987E-05, 0.546986E-05, 0.480064E-05, 0.397211E-05, &
    0.319003E-05, 0.246208E-05, 0.181795E-05, 0.135296E-05, &
    0.102925E-05, 0.808670E-06, 0.612577E-06, 0.434212E-06, &
    0.328720E-06, 0.252055E-06, 0.198937E-06, 0.166297E-06, &
    0.139094E-06, 0.116418E-06, 0.981116E-07, 0.850660E-07, &
    0.743462E-07, 0.649675E-07, 0.577062E-07, 0.520021E-07 /)

  ! config/planet.yaml composition
  real(RealExt) :: co2_mmr = 6.8348e-4
  real(RealExt) :: ch4_mmr = 9.637e-7
  real(RealExt) :: n2o_mmr = 4.7255e-7
  real(RealExt) :: o2_mmr  = 0.23139

  real(RealExt) :: albedo_sw = 0.12, albedo_lw = 0.0
  real(RealExt) :: psurf = 101300.0, ptop = 30.0
  real(RealExt) :: zsk, s, sh_lo, sh_hi
  integer :: i, l, irep
  integer(kind=8) :: c0, c1, crate
  real(kind=8) :: t_setup, t_sw, t_lw

  call getenv_i('SB_NLAYER',    n_layer,    10)
  call getenv_i('SB_NPROFILE',  n_profile,  2048)
  call getenv_i('SB_NREP',      n_rep,      20)
  call getenv_i('SB_CLOUD',     i_cloud,    0)
  call getenv_i('SB_SETUPONLY', setup_only, 0)
  call getenv_i('SB_ZENITH',    i_zenith,   0)
  call getenv_i('SB_PARAM_WATER', i_param_water, 5)
  call getenv_i('SB_PARAM_ICE',   i_param_ice,   8)
  call get_environment_variable('SB_SW', buf)
  sw_file = trim(buf); if (len_trim(sw_file)==0) sw_file = 'sp_sw_ga7'
  call get_environment_variable('SB_LW', buf)
  lw_file = trim(buf); if (len_trim(lw_file)==0) lw_file = 'sp_lw_ga7'

  allocate(p_layer(n_profile,n_layer), t_layer(n_profile,n_layer))
  allocate(p_level(n_profile,0:n_layer), t_level(n_profile,0:n_layer))
  allocate(h2o(n_profile,n_layer), o3(n_profile,n_layer))
  allocate(d_mass(n_profile,n_layer), density(n_profile,n_layer))
  allocate(lhc(n_profile,n_layer))
  allocate(t_ground(n_profile), cosz(n_profile), sirr(n_profile))
  allocate(sw_hr(n_profile,n_layer), lw_hr(n_profile,n_layer))
  allocate(sw_up(n_profile,0:n_layer), sw_dn(n_profile,0:n_layer))
  allocate(lw_up(n_profile,0:n_layer), lw_dn(n_profile,0:n_layer))
  allocate(cld_frac(n_profile,n_layer), liq_frac(n_profile,n_layer))
  allocate(ice_frac(n_profile,n_layer), liq_mmr(n_profile,n_layer))
  allocate(ice_mmr(n_profile,n_layer), liq_dim(n_profile,n_layer))
  allocate(ice_dim(n_profile,n_layer))

  ! Half levels on the model's quartic sigma, plasim.f90:2150, with the top
  ! half-level pulled off the lid so the first layer has finite mass.
  do i = 0, n_layer
    zsk = real(i,RealExt)/real(n_layer,RealExt)
    s = 0.75*zsk + 1.75*zsk**3 - 1.5*zsk**4
    do l = 1, n_profile
      p_level(l,i) = max(ptop, s*psurf)
    end do
  end do
  do i = 1, n_layer
    do l = 1, n_profile
      p_layer(l,i) = 0.5*(p_level(l,i-1) + p_level(l,i))
    end do
    call interp_logp(pm, tm, nm, p_layer(1,i), t_layer(:,i), n_profile)
    call interp_logp(pm, qm, nm, p_layer(1,i), h2o(:,i),     n_profile)
    call interp_logp(pm, om, nm, p_layer(1,i), o3(:,i),      n_profile)
  end do
  do i = 0, n_layer
    call interp_logp(pm, tm, nm, p_level(1,i), t_level(:,i), n_profile)
  end do
  do i = 1, n_layer
    do l = 1, n_profile
      d_mass(l,i)  = (p_level(l,i)-p_level(l,i-1))/grav_acc
      density(l,i) = p_layer(l,i)/(r_gas_dry*t_layer(l,i))
      lhc(l,i)     = d_mass(l,i)*cp_air_dry
    end do
  end do
  t_ground(:) = t_level(:,n_layer)
  ! Spread the zenith cosine and the insolation over the columns so the sunlit
  ! and dark branches are both exercised, as they are over a real latitude band.
  do l = 1, n_profile
    if (i_zenith == 1) then
      cosz(l) = 0.5
    else
      cosz(l) = max(0.0_RealExt, cos( (real(l-1,RealExt)/real(n_profile,RealExt)) &
                                      * 2.0*pi ) )
    end if
    sirr(l) = 321.538
  end do

  cld_frac = 0.0; liq_frac = 0.0; ice_frac = 0.0
  liq_mmr = 0.0; ice_mmr = 0.0; liq_dim = 1.0e-5; ice_dim = 3.0e-5
  if (i_cloud == 1) then
    do i = 1, n_layer
      do l = 1, n_profile
        if (p_layer(l,i) > 0.4*psurf .and. p_layer(l,i) < 0.95*psurf) then
          cld_frac(l,i) = 0.5
          liq_frac(l,i) = 1.0
          liq_mmr(l,i)  = 5.0e-5
        else if (p_layer(l,i) <= 0.4*psurf .and. p_layer(l,i) > 0.15*psurf) then
          cld_frac(l,i) = 0.3
          ice_frac(l,i) = 1.0
          ice_mmr(l,i)  = 1.0e-5
        end if
      end do
    end do
  end if

  call system_clock(count_rate=crate)

  call system_clock(c0)
  call set_spectrum(spectrum_name='sw', spectral_file=trim(sw_file), &
                    l_all_gases=.true.)
  call set_spectrum(spectrum_name='lw', spectral_file=trim(lw_file), &
                    l_all_gases=.true.)
  call system_clock(c1)
  t_setup = real(c1-c0,8)/real(crate,8)

  sw_diag%heating_rate => sw_hr
  sw_diag%flux_up      => sw_up
  sw_diag%flux_down    => sw_dn
  lw_diag%heating_rate => lw_hr
  lw_diag%flux_up      => lw_up
  lw_diag%flux_down    => lw_dn

  write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'nprofile ', n_profile, '  nlayer ', &
    n_layer, '  nrep ', n_rep, '  cloud ', i_cloud, '  zenith ', i_zenith
  write(*,'(a,a,1x,a)') 'spectra  ', trim(sw_file), trim(lw_file)
  write(*,'(a,f10.4)') 'setup_s  ', t_setup
  if (setup_only == 1) stop

  ! Warm-up, discarded: first call touches the k-tables and the allocators.
  call do_sw(); call do_lw()

  call system_clock(c0)
  do irep = 1, n_rep
    call do_sw()
  end do
  call system_clock(c1)
  t_sw = real(c1-c0,8)/real(crate,8)

  call system_clock(c0)
  do irep = 1, n_rep
    call do_lw()
  end do
  call system_clock(c1)
  t_lw = real(c1-c0,8)/real(crate,8)

  write(*,'(a,es14.6)') 'sw_s_total    ', t_sw
  write(*,'(a,es14.6)') 'lw_s_total    ', t_lw
  write(*,'(a,es14.6)') 'sw_s_per_col_call ', t_sw/real(n_rep,8)/real(n_profile,8)
  write(*,'(a,es14.6)') 'lw_s_per_col_call ', t_lw/real(n_rep,8)/real(n_profile,8)
  write(*,'(a,es14.6)') 'both_s_per_col_call ', &
    (t_sw+t_lw)/real(n_rep,8)/real(n_profile,8)
  write(*,'(a,f12.6)') 'check_sw_dn_sfc ', sw_dn(1,n_layer)
  write(*,'(a,f12.6)') 'check_lw_up_toa ', lw_up(1,0)

contains

  subroutine do_sw()
    call runes(n_profile=n_profile, n_layer=n_layer, diag=sw_diag, &
      spectrum_name='sw', i_source=ip_source_illuminate, &
      p_layer=p_layer, t_layer=t_layer, p_level=p_level, &
      mass=d_mass, density=density, layer_heat_capacity=lhc, &
      h2o=h2o, o3=o3, co2_mix_ratio=co2_mmr, n2o_mix_ratio=n2o_mmr, &
      ch4_mix_ratio=ch4_mmr, o2_mix_ratio=o2_mmr, &
      cos_zenith_angle=cosz, solar_irrad=sirr, &
      l_grey_albedo=.true., grey_albedo=albedo_sw, &
      l_rayleigh=.true., l_invert=.false., &
      cloud_frac=cld_frac, liq_frac=liq_frac, ice_frac=ice_frac, &
      liq_mmr=liq_mmr, ice_mmr=ice_mmr, liq_dim=liq_dim, ice_dim=ice_dim, &
      i_cloud_representation=cloud_rep(), i_overlap=ip_overlap_max_random, &
      i_inhom=ip_inhom_homogeneous, &
      i_st_water=i_param_water, i_st_ice=i_param_ice)
  end subroutine do_sw

  subroutine do_lw()
    call runes(n_profile=n_profile, n_layer=n_layer, diag=lw_diag, &
      spectrum_name='lw', i_source=ip_source_thermal, &
      p_layer=p_layer, t_layer=t_layer, p_level=p_level, t_level=t_level, &
      t_ground=t_ground, &
      mass=d_mass, density=density, layer_heat_capacity=lhc, &
      h2o=h2o, o3=o3, co2_mix_ratio=co2_mmr, n2o_mix_ratio=n2o_mmr, &
      ch4_mix_ratio=ch4_mmr, &
      l_grey_albedo=.true., grey_albedo=albedo_lw, l_invert=.false., &
      cloud_frac=cld_frac, liq_frac=liq_frac, ice_frac=ice_frac, &
      liq_mmr=liq_mmr, ice_mmr=ice_mmr, liq_dim=liq_dim, ice_dim=ice_dim, &
      i_cloud_representation=cloud_rep(), i_overlap=ip_overlap_max_random, &
      i_inhom=ip_inhom_homogeneous, &
      i_st_water=i_param_water, i_st_ice=i_param_ice)
  end subroutine do_lw

  integer function cloud_rep()
    if (i_cloud == 1) then
      cloud_rep = ip_cloud_representation_ice_water
    else
      cloud_rep = ip_cloud_representation_off
    end if
  end function cloud_rep

  subroutine getenv_i(name, out, dflt)
    character(len=*), intent(in) :: name
    integer, intent(out) :: out
    integer, intent(in) :: dflt
    character(len=64) :: v
    integer :: ios
    call get_environment_variable(name, v)
    if (len_trim(v) == 0) then
      out = dflt
    else
      read(v,*,iostat=ios) out
      if (ios /= 0) out = dflt
    end if
  end subroutine getenv_i

  subroutine interp_logp(px, vx, n, ptarget, vout, m)
    integer, intent(in) :: n, m
    real(RealExt), intent(in) :: px(n), vx(n), ptarget
    real(RealExt), intent(out) :: vout(m)
    integer :: j
    real(RealExt) :: w, v
    if (ptarget <= px(1)) then
      v = vx(1)
    else if (ptarget >= px(n)) then
      v = vx(n)
    else
      do j = 1, n-1
        if (ptarget >= px(j) .and. ptarget <= px(j+1)) then
          w = (log(ptarget)-log(px(j)))/(log(px(j+1))-log(px(j)))
          v = vx(j)*(1.0-w) + vx(j+1)*w
          exit
        end if
      end do
    end if
    vout(:) = v
  end subroutine interp_logp

end program socrates_bench
