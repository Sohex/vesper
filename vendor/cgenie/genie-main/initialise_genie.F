      subroutine initialise_genie()
c      
      use genie_control
      use genie_global
      use genie_util, ONLY: check_unit, check_iostat, message, die
c
      implicit none
c
c     no args
c
c     locals
c
      integer ios
      integer :: unitNum = 8
c
c     Define members of the top-level control namelist
c
      namelist/genie_control_nml/expid_name,
     :      koverall_total,
     :      katm_loop,
     :      ksic_loop,
     :      kocn_loop,
     :      klnd_loop,
     :      conv_kocn_katchem,
     :      conv_kocn_kbiogem,
     :      conv_kocn_ksedgem,
     :      conv_kocn_krokgem,
     :      conv_kocn_kgoldlite,
     :      conv_kocn_kocnlite,
     :      conv_kocn_kecogem,
     :      conv_kocn_korb,
     :      kgemlite,
     :      flag_ebatmos,flag_plasimatmos,flag_fluxatmos,
     :      flag_goldsteinocean,
     :      flag_goldsteinseaice,
     :      flag_ents,
     :      flag_biogem,
     :      flag_atchem,
     :      flag_sedgem,
     :      flag_rokgem,
     :      flag_ecogem,
     :      flag_goldlite,
     :      flag_ocnlite,
     :      flag_gemlite,
     :      write_flag_atm,write_flag_ocn,write_flag_sic,
     :      outputdir_name,dt_write,
     :      fname_restart_main,
     :      fname_fluxrestart,
     :      flag_checkfluxes_sic,
     :      flag_checkfluxes_ocn,
     :      flag_checkfluxes_surf,
     :      flag_checkfluxes_atlantic,
     :      genie_timestep,
     :      lrestart_genie,
     :      genie_solar_constant,
     :      fname_topo,
     :      verbosity,
     :      debug_init,debug_end,debug_loop,
     :      gem_yr,gem_yr_min,gem_yr_max,
     :      gem_notyr,gem_notyr_min,gem_notyr_max,
     :      gem_dyr,
     :      gem_adapt_auto,gem_adapt_dpCO2dt,gem_adapt_DpCO2,
     :      gem_adapt_auto_unlimitedGEM,
     :      gem_adapt_diag_biogem_full
c
c     assign default values to members of the namelist
c
      koverall_total=0
      conv_kocn_katchem=-1
      conv_kocn_ksedgem=-1
      conv_kocn_kbiogem=-1
      conv_kocn_krokgem=-1
      conv_kocn_kgoldlite=-1
      conv_kocn_kocnlite=-1
      conv_kocn_kecogem=-1
      conv_kocn_korb=-1
      kgemlite=1
c
      flag_plasimatmos=.FALSE.
      flag_ebatmos=.FALSE.
      flag_fluxatmos=.FALSE.
c
      flag_goldsteinocean=.TRUE.
      flag_goldlite=.FALSE.
      flag_ocnlite=.FALSE.
c
      flag_goldsteinseaice=.TRUE.
c
      flag_ents=.FALSE.
c
      flag_biogem=.false.
      flag_atchem=.false.
      flag_sedgem=.false.
      flag_rokgem=.false.
      flag_ecogem=.false.
      flag_gemlite=.false.
c
      write_flag_atm=.true.
      write_flag_ocn=.true.
      write_flag_sic=.true.
      outputdir_name='~/genie_output/main'
      fname_restart_main='xxx'
      fname_fluxrestart='xxx'
      dt_write=1
c
      flag_checkfluxes_sic=.false.
      flag_checkfluxes_ocn=.false.
      flag_checkfluxes_surf=.false.
      flag_checkfluxes_atlantic=.false.
c
      genie_timestep = 3600.0
      genie_clock    = 0
c
      lrestart_genie=.false.
c
      genie_solar_constant=1368.0
c
      verbosity = 0
c
      debug_init = 1
c
      fname_topo='worbe2'
c
c open the file containing the namelist
c
      call check_unit(unitNum,__LINE__,__FILE__)
      open(unit=unitNum,file='data_genie',status='old',iostat=ios)
      if (ios /= 0) then
         call die("could not open GENIE control namelist file",
     &  __LINE__, __FILE__)
      end if
c
c read in namelist
c
      read(UNIT=unitNum,NML=genie_control_nml,IOSTAT=ios)
      if (ios /= 0) then
         call die('could not read GENIE control namelist',
     & __LINE__, __FILE__)
      else
         close(unitNum,iostat=ios)
         call check_iostat(ios,__LINE__,__FILE__)
      end if
c
      if (debug_init > 1)
     & call message("read values for genie_control_nml:",3)
      if(verbosity.ge.3) then
         write(UNIT=6,NML=genie_control_nml)
      end if
c
      if (debug_init > 1)
     & call message("read namelist: genie_control_nml",1)
c
c     The driven surface-flux path. flag_fluxatmos gates the third path in
c     genie.F: it advances istep_ocn and calls surflux_goldstein_seaice with
c     coupling fields supplied from outside, rather than with an atmosphere in
c     this executable computing them. It is an ATMOSPHERE flag and it names the
c     same slot the other two do, so more than one of the three is not a
c     configuration: istep_ocn would advance once per gate that is true.
c
      if (flag_fluxatmos.and.(flag_ebatmos.or.flag_plasimatmos)) then
         call die("flag_fluxatmos is an atmosphere and cannot be set "//
     &    "beside flag_ebatmos or flag_plasimatmos: istep_ocn would "//
     &    "advance once per gate", __LINE__, __FILE__)
      end if
c
c     The SUPPLY half is not in this tree. genie.F's gate fills nothing: the
c     coupling fields the path passes to surflux_goldstein_seaice, and the wind
c     stress GOLDSTEIN reads, are module arrays that no routine here writes when
c     no atmosphere runs. Running anyway would integrate an ocean forced by
c     whatever those arrays hold and produce output shaped exactly like a
c     result, so the flag is refused instead. Deleting this block is what
c     turning the path on consists of, once the forcing reader that fills the
c     fields listed in section 1d of
c     notes/audits/cgenie-embm-free-path-and-threading.md exists.
c
      if (flag_fluxatmos) then
         call die("flag_fluxatmos: the driven surface-flux path is "//
     &    "gated, flagged and registered, but the forcing supply that "//
     &    "fills its coupling fields is not in this tree", __LINE__,
     &    __FILE__)
      end if
c
      return
      end
