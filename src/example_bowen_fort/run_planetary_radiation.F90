PROGRAM RUN_PLANETARY_RADIATION
!driver file to generate pressure grid and state profiles as initial
!input to planetary_radiation, iterate the model, and save at specific
!timesteps
!N. Lombardo Jan 22 2020 - file created
!C. Keaveney Jul 30 2024 - file updated

  use planetary_radiation_driver_mod

  implicit none
 
  integer, parameter              :: nlay = 100      ! Modified by CK [Originally 100]
  real, parameter                 :: high = 145600.  ! Modified by CK [Originally 145600]
  real                            :: logstep, p_surf
  integer                         :: i, j, k, L_LEVELS, L_LAYERS
  real, dimension(1, 1, nlay+1)   :: phalf, zhalf
  real, dimension(1, 1, nlay)     :: pfull, zfull, t, t2
  real, dimension(2*nlay+3)       :: plev
  real, dimension(1, 1)           :: lat, lon, t_surf, cosz, H, cosz_t
  real, dimension(1, 1, nlay, 1)  :: tracers
  real                            :: low = 1.0
  integer                         :: num_i = 3000
  integer                         :: write_every = 1
  real                            :: time_step = 1000.
  integer                         :: is, js, ios
  logical                         :: have_surface_restart, have_top_tlev_restart
  real, dimension(1,1)            :: albedov, albedoi
  real, dimension(1,1)            :: sfc_flux
  real, dimension(1,1,nlay)       :: tdt, tdt_sw, tdt_lw
  real, dimension(1, 1, 41)       :: lw_spec
  real, dimension(1, 1, 43)       :: sw_spec
  real                            :: pi = 3.1415926535 ! Added by CK
  real                            :: testing_lat, init_t, total_time, arcsin, testing_ls, theta
  real                            :: hour_angle, declination ! Added by CK
  real                            :: H0_day, diurnal_arg      ! diurnal-average insolation
  logical                         :: diurnal = .false.        ! .false. = diurnally-averaged
  ! Planetary Parameters. Callable from namelist, defaults for Uranus (Added by CK on 6/26/24)
  real                            :: seconds_per_day = 62064.   ! [Uranus: 62064.  | Titan: 1375200.]
  real                            :: obliq = 97.77              ! [Uranus: 97.77   | Titan: 27]
  real                            :: grav = 8.94                ! [Uranus: 8.94    | Titan: 1.352]
  real                            :: Rdgas = 3149.2             ! [Uranus: 3149.2  | Titan: 287.05]
  real                            :: cp_air = 11022.2           ! [Uranus: 11022.2 | Titan: 1039.5]
  real                            :: albv = 0.0                 ! [Uranus: 0.0     | Titan: 0.2]
  real                            :: albi = 0.0                 ! [Uranus: 0.0     | Titan: 0.05]
  ! Extra net surface-energy term added to the radiative surface flux (W/m2).
  ! Positive values warm the surface directly.
  real                            :: surface_flux_offset = 0.042
  logical                         :: use_bulk_surface_flux = .false.
  real                            :: surface_wind_speed = 0.03
  real                            :: surface_drag_coeff = 0.001
  ! Huygens-like near-surface CH4: T = 93.5 K, vmr = 0.049 -> RH ~= 0.52
  ! under the simple saturation law used below.
  real                            :: surface_relative_humidity = 0.52
  real                            :: surface_methane_availability = 1.0
  logical                         :: use_analytic_ch4_profile = .false.
  real                            :: ch4_inventory_limit_pa = 0.0
  real                            :: ch4_tropopause_lapse_k_per_km = 0.276
  character(len=16)               :: analytic_ch4_saturation_scheme = 'antoine'
  real                            :: analytic_ch4_cc_prefactor = 1.0
  real                            :: surface_thermal_inertia = 601.0
  real                            :: surface_areal_heat_capacity = -1.0
  real                            :: surface_reference_period = -1.0
  real                            :: surface_initial_offset = 0.0
  real                            :: TITAN_HLV = 5.50e5
  real                            :: TITAN_RVGAS = 518.0
  real, parameter                 :: HAN2025_RHOCP = 3.0e6
  real                            :: t_air_low, p_air_low, rho_air_low
  real                            :: q_air_low, q_sat_air, q_sat_surface, q_surface_eq
  real                            :: ch4_effective_surface_rh, ch4_inventory_pa
  real                            :: ch4_lower_q, ch4_tropopause_q, ch4_tropopause_pressure
  real                            :: surface_beta
  real                            :: surface_sens_flux, surface_lat_flux, net_surface_flux
  real                            :: surface_dp
  real                            :: surface_omega, surface_active_layer_depth
  real                            :: surface_max_delta_t = 5.0
  real                            :: surface_delta_t
  real                            :: surface_rad_flux
  real                            :: surface_lw_emission
  real                            :: surface_emissivity
  real, parameter                 :: STEFAN = 5.670374419e-8
  ! Thermosphere conduction
  real                            :: Q_thermo = 1.36e-3
  real                            :: dTdz_top
  real                            :: dTdz_bot
  real                            :: k_top
  real                            :: k_bot
  ! Convective adjustment / diagnostics
  character(len=16)               :: convective_scheme = 'dry'
  real                            :: moist_conv_tau = 7200.0
  real                            :: moist_conv_rhbm = 0.8
  real                            :: moist_conv_min_cape = 0.0
  logical                         :: use_dry_convective_fallback = .true.
  logical                         :: write_convective_diagnostics = .true.
  real                            :: convective_cape
  real                            :: convective_cin
  real                            :: convective_p_lcl
  real                            :: convective_p_lfc
  real                            :: convective_p_lzb
  real                            :: convective_max_delta_t
  integer                         :: convective_flag
  real, dimension(3)              :: top_tlev_restart
  real, dimension(nlay)           :: ch4_profile_q, ch4_profile_vmr

  namelist/run_planetary_radiation_nml/     num_i, write_every, time_step, testing_lat, testing_ls, init_t, low, &
                                            diurnal, &
                                            seconds_per_day, obliq, grav, Rdgas, cp_air, albv, albi, &
                                            surface_flux_offset, use_bulk_surface_flux, surface_wind_speed, &
                                            surface_drag_coeff, surface_relative_humidity, &
                                            surface_methane_availability, use_analytic_ch4_profile, &
                                            ch4_inventory_limit_pa, ch4_tropopause_lapse_k_per_km, &
                                            analytic_ch4_saturation_scheme, analytic_ch4_cc_prefactor, &
                                            convective_scheme, moist_conv_tau, moist_conv_rhbm, &
                                            moist_conv_min_cape, use_dry_convective_fallback, &
                                            write_convective_diagnostics, &
                                            surface_thermal_inertia, &
                                            surface_areal_heat_capacity, surface_reference_period, &
                                            surface_initial_offset, surface_max_delta_t, Q_thermo
  
  ! Simpler namelist reading for testing
  open(5, file = 'namelist')
  read(5, run_planetary_radiation_nml)
  close(5)

  analytic_ch4_saturation_scheme = lowercase_string(trim(adjustl(analytic_ch4_saturation_scheme)))
  select case (trim(analytic_ch4_saturation_scheme))
  case ('antoine')
     continue
  case ('cc')
     analytic_ch4_cc_prefactor = max(0.0, analytic_ch4_cc_prefactor)
  case ('cc0.9')
     analytic_ch4_saturation_scheme = 'cc'
     analytic_ch4_cc_prefactor = 0.9
  case ('cc1.0')
     analytic_ch4_saturation_scheme = 'cc'
     analytic_ch4_cc_prefactor = 1.0
  case default
     stop 'Unsupported analytic_ch4_saturation_scheme; use antoine, cc, cc0.9, or cc1.0'
  end select

  convective_scheme = lowercase_string(trim(adjustl(convective_scheme)))
  select case (trim(convective_scheme))
  case ('dry', 'dry_adjustment')
     convective_scheme = 'dry'
  case ('moist_bm', 'moist', 'betts_miller')
     convective_scheme = 'moist_bm'
  case ('none')
     continue
  case default
     stop 'Unsupported convective_scheme; use dry, moist_bm, or none'
  end select

  if (trim(convective_scheme) .eq. 'moist_bm' .and. .not. use_analytic_ch4_profile) then
     print *, 'Moist BM convection requires analytic CH4; enabling use_analytic_ch4_profile.'
     use_analytic_ch4_profile = .true.
  end if

  if (surface_reference_period .le. 0.0) surface_reference_period = seconds_per_day
  surface_emissivity = max(0.0, min(1.0, 1.0 - albi))
  surface_omega = 2.0*pi / max(surface_reference_period, time_step)
  if (surface_areal_heat_capacity .le. 0.0) then
     ! Han et al. (2025): Huygens-region thermal inertia is O(600 TIU), which
     ! corresponds to an effective diurnally active layer of only a few 0.1 m.
     surface_areal_heat_capacity = surface_thermal_inertia / sqrt(surface_omega)
  endif
  surface_active_layer_depth = surface_areal_heat_capacity / HAN2025_RHOCP
  
  !equal log-space grid in Pa
  logstep = (log(high) - log(low)) / nlay
  phalf(1, 1, :) = (/(low * EXP(i*logstep) , i = 0, nlay)/)
  pfull(1, 1, :) = (/(low * EXP((i+0.5)*logstep) , i = 0, nlay - 1)/)

  p_surf = phalf(1, 1, nlay+1)
 
  ! Set init_t = 0. in namelist to read initial profile from file. Added by CK on 6/26/24
  if (init_t .eq. 0.) then
     inquire(file = 'initial_temperatures.txt', exist = have_surface_restart)
     have_top_tlev_restart = .false.
     if (have_surface_restart) then
        print*, 'Reading atmospheric and surface restart from initial_temperatures.txt...'
        open(10, file = 'initial_temperatures.txt', status = 'old', iostat = ios)
        if (ios .ne. 0) stop 'Unable to open initial_temperatures.txt'
        read(10, *, iostat = ios) (t(1, 1, i), i=1, nlay)
        if (ios .ne. 0) stop 'Unable to read atmospheric restart from initial_temperatures.txt'
        read(10, *, iostat = ios) t_surf(1, 1)
        if (ios .ne. 0) stop 'Unable to read surface restart from initial_temperatures.txt'
        read(10, *, iostat = ios) top_tlev_restart
        if (ios .eq. 0) then
           have_top_tlev_restart = .true.
           print*, 'Reading top ghost-layer restart from initial_temperatures.txt...'
        endif
        close(10)
     else
        print*, 'Reading in initial atmospheric temperature profile...'
        open(10, file = 'initial_temperature.txt', status = 'old', iostat = ios)
        if (ios .ne. 0) stop 'Unable to open initial_temperature.txt'
        read(10, *, iostat = ios) (t(1, 1, i), i=1, nlay)
        if (ios .ne. 0) stop 'Unable to read atmospheric restart from initial_temperature.txt'
        close(10)
        t_surf(1, 1) = t(1, 1, nlay) + surface_initial_offset
     endif
     call set_top_tlev_restart(top_tlev_restart, have_top_tlev_restart)
  else
     print*, 'Initializing with isothermal profile...'
     t(1, 1, :) = (/(init_t, i=1, nlay)/)
     t_surf(1, 1) = init_t + surface_initial_offset
     call set_top_tlev_restart(top_tlev_restart, .false.)
  endif

  print *, 'Initial t read'

  open (unit = 30, file = 'temperatures.txt')
  write(30, ' (  500(X, E12.5) ) ') pfull
  
  print *, 'Surface areal heat capacity (J/m2/K): ', surface_areal_heat_capacity
  print *, 'Equivalent active-layer depth (m): ', surface_active_layer_depth
  lat(1, 1) = testing_lat
  lon(1, 1) = 0.

  is = 1
  js = 1
  albedov(1, :) = (/(albv)/) 
  albedoi(1, :) = (/(albi)/) 
    
  L_LAYERS = nlay
  L_LEVELS = 2*nlay+3
    
  !Calculate plev
  do k = 1,L_LAYERS
    plev(2*k+1) = phalf(is,js,k)
    plev(2*k+2) = pfull(is,js,k)
  end do
  plev(L_LEVELS) = phalf(is, js, L_LAYERS+1)

  if (plev(3) .eq. 0.0) plev(3) = plev(4)/10.

  !plev extrapolation may be causing problem later, not actually following grid
  plev(2) = plev(3)*plev(3)/(plev(4))
  plev(1) = plev(2)*plev(2)/plev(3)

  tracers = 0.0
  ch4_profile_q = 0.0
  ch4_profile_vmr = 0.0
  ch4_effective_surface_rh = max(0.0, min(1.0, surface_relative_humidity))
  ch4_inventory_pa = 0.0
  ch4_lower_q = 0.0
  ch4_tropopause_q = 0.0
  ch4_tropopause_pressure = pfull(1,1,1)

  call set_analytic_ch4_mode(use_analytic_ch4_profile)
  
  call planetary_radiation_init(nlay, lat(:, 1), phalf(1,1,:), pfull(1,1,:))
  
  open (unit = 99, file = 'sw.txt')
  open (unit = 20, file = 'lw.txt')
  open (unit = 40, file = 'lwspec.txt')
  open (unit = 50, file = 'swspec.txt')
  open (unit = 60, file = 'tsurf.txt')
  open (unit = 61, file = 'surface_fluxes.txt')
  if (write_convective_diagnostics) then
     open (unit = 65, file = 'convective_diagnostics.txt')
     write(65, '(A)') '# time_days cape_jpkg cin_jpkg p_lcl_pa p_lfc_pa p_lzb_pa max_delta_t_k convflag'
  end if
  if (use_analytic_ch4_profile) then
     open (unit = 63, file = 'ch4_inventory.txt')
     open (unit = 64, file = 'ch4_vmr.txt')
     write(64, ' (  500(X, E12.5) ) ') pfull
  end if

  total_time = 0.
  surface_sens_flux = 0.0
  surface_lat_flux = 0.0
  net_surface_flux = 0.0
  convective_cape = 0.0
  convective_cin = 0.0
  convective_p_lcl = 0.0
  convective_p_lfc = 0.0
  convective_p_lzb = 0.0
  convective_max_delta_t = 0.0
  convective_flag = 0
  do  i=1, num_i
    if (modulo(i-1, write_every) .eq. 0) then
        print *, 'Calling planetary radiation, iteration ', i
    end if

    if (use_analytic_ch4_profile) then
       call build_simplified_ch4_profile(t(1,1,:), t_surf(1,1), pfull(1,1,:), phalf(1,1,:), p_surf, &
                                         surface_relative_humidity, ch4_inventory_limit_pa, ch4_profile_q, &
                                         ch4_profile_vmr, ch4_effective_surface_rh, ch4_inventory_pa, &
                                         ch4_lower_q, ch4_tropopause_q, ch4_tropopause_pressure)
       tracers(1,1,:,1) = ch4_profile_q
       q_air_low = ch4_profile_q(nlay)
    else
       tracers = 0.0
       q_air_low = 0.0
    end if

    t_air_low = t(1,1,nlay)
    p_air_low = pfull(1,1,L_LAYERS)
    surface_dp = max(pfull(1,1,L_LAYERS)-pfull(1,1,L_LAYERS-1), 1.0)
    
    !Calculate z grid from hydrostatic equilibrium, begin at surface
    zhalf(1, 1, nlay+1) = 0.
    zfull(1, 1, nlay) = zhalf(1, 1, nlay+1) + (Rdgas * (log(phalf(1, 1, nlay+1)) - log(pfull(1, 1, nlay)))*t(1, 1, nlay))/grav
    do j=nlay-1, 1, -1
        zhalf(1, 1, j+1) = zfull(1, 1, j+1) + (Rdgas * (log(pfull(1, 1, j+1))-log(phalf(1, 1, j+1)))*t(1, 1, j+1))/(grav)
        zfull(1, 1, j) = zhalf(1, 1, j+1) + (Rdgas * (log(phalf(1, 1, j+1)) - log(pfull(1, 1, j)))*t(1, 1, j))/(grav)
    end do
    zhalf(1, 1, 1) = zfull(1, 1, 1) + (Rdgas * (log(pfull(1, 1, 1))-log(phalf(1, 1, 1)))*t(1, 1, 1))/(grav)
    
    !---------------------
    ! Solar zenith angle
    !---------------------
    declination = obliq*pi/180.*real(testing_ls)
    if (diurnal) then
       ! instantaneous (resolves the diurnal cycle; needs sub-day timestep)
       hour_angle = 2*pi*total_time / seconds_per_day
       cosz = COS(testing_lat*pi/180.)*COS(declination)*COS(hour_angle) + SIN(declination)*SIN(testing_lat*pi/180.)
    else
       ! diurnally-AVERAGED insolation factor (constant) -- allows a full-day
       ! timestep and faster convergence to steady state
       diurnal_arg = -TAN(testing_lat*pi/180.)*TAN(declination)
       if (diurnal_arg >= 1.0) then
          H0_day = 0.0                       ! polar night
       else if (diurnal_arg <= -1.0) then
          H0_day = pi                        ! polar day
       else
          H0_day = ACOS(diurnal_arg)
       end if
       cosz = (1.0/pi)*( H0_day*SIN(testing_lat*pi/180.)*SIN(declination) &
                       + COS(testing_lat*pi/180.)*COS(declination)*SIN(H0_day) )
    end if
        
    call planetary_radiation(is, js, lat, lon, testing_ls, cosz, phalf, pfull, zhalf, zfull, albedov, albedoi, &
    t_surf, t, tracers, tdt, tdt_sw, tdt_lw, sfc_flux, lw_spec, sw_spec)

    if (use_bulk_surface_flux) then
       rho_air_low = p_air_low / (Rdgas * max(40.0, t_air_low))
       q_sat_air = methane_qsat(t_air_low, p_air_low, Rdgas, TITAN_RVGAS)
       if (.not. use_analytic_ch4_profile) then
          q_air_low = max(0.0, min(q_sat_air, surface_relative_humidity*q_sat_air))
       end if
       q_sat_surface = methane_qsat(t_surf(1,1), p_surf, Rdgas, TITAN_RVGAS)
       surface_beta = max(0.0, min(1.0, surface_methane_availability))
       ! Dry surface: beta = 0 -> q_surface_eq = q_air_low and latent flux vanishes.
       ! Fully wet surface: beta = 1 -> q_surface_eq = q_sat_surface.
       q_surface_eq = (1.0 - surface_beta) * q_air_low + surface_beta * q_sat_surface
       ! Positive sensible and latent fluxes remove energy from the surface.
       surface_sens_flux = rho_air_low * cp_air * surface_drag_coeff * max(0.0, surface_wind_speed) * &
                           (t_surf(1,1) - t_air_low)
       surface_lat_flux = rho_air_low * TITAN_HLV * surface_drag_coeff * max(0.0, surface_wind_speed) * &
                          (q_surface_eq - q_air_low)
       tdt(1,1,L_LAYERS) = tdt(1,1,L_LAYERS) + grav/cp_air * surface_sens_flux / surface_dp
       tdt(1,1,L_LAYERS) = tdt(1,1,L_LAYERS) + grav/cp_air * surface_lat_flux / surface_dp
    else
       surface_sens_flux = 0.0
       surface_lat_flux = 0.0
       q_sat_air = 0.0
       q_sat_surface = 0.0
       q_surface_eq = 0.0
    end if
    surface_lw_emission = surface_emissivity * STEFAN * t_surf(1,1)**4
    surface_rad_flux = sfc_flux(1,1) - surface_lw_emission
    net_surface_flux = surface_rad_flux + surface_flux_offset - surface_sens_flux - surface_lat_flux
    
    !------------------------
    ! Thermosphere conduction
    !------------------------
    tdt(1,1,1) = tdt(1,1,1) + grav/cp_air * Q_thermo/(pfull(1,1,2)-pfull(1,1,1))
    do k=2,L_LAYERS-1
        dTdz_top = (t(1,1,k-1)-t(1,1,k))/(zfull(1,1,k-1)-zfull(1,1,k))
        dTdz_bot = (t(1,1,k)-t(1,1,k+1))/(zfull(1,1,k)-zfull(1,1,k+1))
        k_top = 1.064e-3 * t(1,1,k)**0.906
        k_bot = 1.064e-3 * t(1,1,k+1)**0.906
        tdt(1,1,k) = tdt(1,1,k) + grav/cp_air * (k_top*dTdz_top - k_bot*dTdz_bot)/(phalf(1,1,k+1)-phalf(1,1,k-1))
    enddo

    if (ANY(ABS(time_step*tdt) .ge. 80.)) then
    !if (.false.) then
        t2 = t + 0.1*time_step * tdt
        print *, 'Decreasing timestep to prevent excessive heat/cooling'
    else
        t2 = t + time_step * tdt
    end if
    
    call apply_selected_convection(t2(1,1,:), pfull(1,1,:), phalf(1,1,:), zfull(1,1,:), ch4_profile_q, &
                                   convective_cape, convective_cin, convective_flag, convective_p_lcl, &
                                   convective_p_lfc, convective_p_lzb, convective_max_delta_t)

    if (use_analytic_ch4_profile) then
       call build_simplified_ch4_profile(t2(1,1,:), t_surf(1,1), pfull(1,1,:), phalf(1,1,:), p_surf, &
                                         surface_relative_humidity, ch4_inventory_limit_pa, ch4_profile_q, &
                                         ch4_profile_vmr, ch4_effective_surface_rh, ch4_inventory_pa, &
                                         ch4_lower_q, ch4_tropopause_q, ch4_tropopause_pressure)
       tracers(1,1,:,1) = ch4_profile_q
    end if

    if ((modulo(i-1, write_every) .eq. 0)) then
        write(99, ' (  500(X, E10.3) ) ') tdt_sw
        write(20, ' (  500(X, E10.3) ) ') tdt_lw
        write(30, ' (  500(X, F7.3) ) ') t2
        write(40, ' (  500(X, E10.3) ) ') lw_spec
        write(50, ' (  500(X, E10.3) ) ') sw_spec
        write(60, ' (  500(X, F7.3) ) ') t_surf
        write(61, ' (  6(X, E12.5) ) ') sfc_flux(1,1), surface_flux_offset, surface_sens_flux, &
                                         surface_lat_flux, net_surface_flux, q_air_low
        if (write_convective_diagnostics) then
            write(65, '(7(1X, E12.5), 1X, I4)') total_time/seconds_per_day, convective_cape, convective_cin, &
                                                 convective_p_lcl, convective_p_lfc, convective_p_lzb, &
                                                 convective_max_delta_t, convective_flag
        end if
        if (use_analytic_ch4_profile) then
            write(63, ' (  7(X, E12.5) ) ') total_time/seconds_per_day, ch4_effective_surface_rh, &
                                             ch4_inventory_pa, ch4_inventory_limit_pa, &
                                             ch4_lower_q, ch4_tropopause_q, ch4_tropopause_pressure
            write(64, ' (  500(X, E12.5) ) ') ch4_profile_vmr
        end if
    end if
    surface_delta_t = time_step * net_surface_flux / max(surface_areal_heat_capacity, 1.0)
    if (abs(surface_delta_t) .gt. surface_max_delta_t) then
        print *, 'Limiting surface temperature step from ', surface_delta_t, ' K'
        surface_delta_t = sign(surface_max_delta_t, surface_delta_t)
    end if
    t_surf(1,1) = t_surf(1,1) + surface_delta_t
    t = t2
    total_time = total_time + time_step
  end do

  call get_top_tlev_restart(top_tlev_restart, have_top_tlev_restart)

  open(unit = 62, file = 'initial_temperatures.txt', status = 'replace')
  write(62, ' (  500(X, E12.5) ) ') t(1, 1, :)
  write(62, ' (  X, E12.5 ) ') t_surf(1, 1)
  if (have_top_tlev_restart) write(62, ' (  3(X, E12.5) ) ') top_tlev_restart
  close(62)
  if (use_analytic_ch4_profile) then
     close(63)
     close(64)
  end if
  if (write_convective_diagnostics) close(65)

contains

  character(len=16) function lowercase_string(text_in)

    character(len=*), intent(in) :: text_in
    integer                      :: ii, ich

    lowercase_string = ' '
    do ii = 1, min(len(text_in), len(lowercase_string))
       ich = iachar(text_in(ii:ii))
       if (ich >= iachar('A') .and. ich <= iachar('Z')) then
          lowercase_string(ii:ii) = achar(ich + 32)
       else
          lowercase_string(ii:ii) = text_in(ii:ii)
       end if
    end do

  end function lowercase_string

  real function methane_sat_pressure(temp)

    real, intent(in) :: temp
    real, parameter  :: ref_pres = 10600.0
    real, parameter  :: ref_temp = 90.68
    real             :: temp_clip

    temp_clip = max(40.0, temp)
    methane_sat_pressure = 0.9 * ref_pres * exp(TITAN_HLV/TITAN_RVGAS * (1.0/ref_temp - 1.0/temp_clip))

  end function methane_sat_pressure

  real function analytic_ch4_sat_pressure(temp)

    real, intent(in) :: temp
    real, parameter  :: ref_pres = 10600.0
    real, parameter  :: ref_temp = 90.68
    real, parameter  :: antoine_a = 3.9895
    real, parameter  :: antoine_b = 443.028
    real, parameter  :: antoine_c = -0.49
    real             :: temp_clip

    temp_clip = max(40.0, temp)
    select case (trim(analytic_ch4_saturation_scheme))
    case ('antoine')
       analytic_ch4_sat_pressure = 1.0e5 * 10.0**(antoine_a - antoine_b / (temp_clip + antoine_c))
    case default
       analytic_ch4_sat_pressure = analytic_ch4_cc_prefactor * ref_pres * &
                                   exp(TITAN_HLV/TITAN_RVGAS * (1.0/ref_temp - 1.0/temp_clip))
    end select

  end function analytic_ch4_sat_pressure

  real function methane_qsat(temp, pressure, rdgas_local, rvgas_local)

    real, intent(in) :: temp, pressure, rdgas_local, rvgas_local
    real             :: es

    es = methane_sat_pressure(temp)
    es = min(es, 0.95*pressure)
    methane_qsat = rdgas_local/rvgas_local * es / (pressure - (1.0-rdgas_local/rvgas_local)*es)

  end function methane_qsat

  real function analytic_ch4_qsat(temp, pressure, rdgas_local, rvgas_local)

    real, intent(in) :: temp, pressure, rdgas_local, rvgas_local
    real             :: es

    es = analytic_ch4_sat_pressure(temp)
    es = min(es, 0.95*pressure)
    analytic_ch4_qsat = rdgas_local/rvgas_local * es / (pressure - (1.0-rdgas_local/rvgas_local)*es)

  end function analytic_ch4_qsat

  real function methane_lapse_rate_k_per_km(temp_profile, pfull_local, kk, grav_local, rdgas_local)

    real, intent(in), dimension(nlay) :: temp_profile, pfull_local
    integer, intent(in)               :: kk
    real, intent(in)                  :: grav_local, rdgas_local
    real                              :: dtdlnp

    if (kk <= 1) then
       dtdlnp = (temp_profile(2) - temp_profile(1)) / (log(pfull_local(2)) - log(pfull_local(1)))
    elseif (kk >= nlay) then
       dtdlnp = (temp_profile(nlay) - temp_profile(nlay-1)) / &
                (log(pfull_local(nlay)) - log(pfull_local(nlay-1)))
    else
       dtdlnp = (temp_profile(kk+1) - temp_profile(kk-1)) / &
                (log(pfull_local(kk+1)) - log(pfull_local(kk-1)))
    end if

    methane_lapse_rate_k_per_km = grav_local / max(rdgas_local * temp_profile(kk), 1.0) * dtdlnp * 1000.0

  end function methane_lapse_rate_k_per_km

  real function methane_vmr_from_q(qval, rdgas_local, rvgas_local)

    real, intent(in) :: qval, rdgas_local, rvgas_local
    real             :: eps

    eps = rdgas_local / rvgas_local
    methane_vmr_from_q = qval / (eps - (eps - 1.0)*qval)

  end function methane_vmr_from_q

  real function methane_inventory_from_profile(q_profile, phalf_local)

    real, intent(in), dimension(nlay)   :: q_profile
    real, intent(in), dimension(nlay+1) :: phalf_local
    integer                             :: kk

    methane_inventory_from_profile = 0.0
    do kk = 1, nlay
       methane_inventory_from_profile = methane_inventory_from_profile + &
                                        q_profile(kk) * max(phalf_local(kk+1) - phalf_local(kk), 0.0)
    end do

  end function methane_inventory_from_profile

  subroutine build_simplified_ch4_profile_for_rh(temp_profile, tsurf_local, pfull_local, phalf_local, psurf_local, &
                                                 surface_rh_use, q_profile, vmr_profile, inventory_pa, q0_out, &
                                                 q_trop_out, p_trop_out)

    real, intent(in), dimension(nlay)   :: temp_profile, pfull_local
    real, intent(in), dimension(nlay+1) :: phalf_local
    real, intent(in)                    :: tsurf_local, psurf_local, surface_rh_use
    real, intent(out), dimension(nlay)  :: q_profile, vmr_profile
    real, intent(out)                   :: inventory_pa, q0_out, q_trop_out, p_trop_out
    logical                             :: in_saturated_layer, in_upper_cap
    real                                :: qsat_layer, q0, rh_use, lapse_rate_local
    integer                             :: kk

    rh_use = max(0.0, min(1.0, surface_rh_use))
    ! Tie lower-layer methane to the near-surface air state, not the surface skin.
    q0 = rh_use * analytic_ch4_qsat(temp_profile(nlay), pfull_local(nlay), Rdgas, TITAN_RVGAS)
    q_profile = 0.0
    q_trop_out = q0
    p_trop_out = pfull_local(1)
    in_saturated_layer = .false.
    in_upper_cap = .false.

    do kk = nlay, 1, -1
       qsat_layer = analytic_ch4_qsat(temp_profile(kk), pfull_local(kk), Rdgas, TITAN_RVGAS)

       if (.not. in_saturated_layer) then
          if (q0 >= qsat_layer) then
             in_saturated_layer = .true.
             q_profile(kk) = qsat_layer
             q_trop_out = qsat_layer
             p_trop_out = pfull_local(kk)
             lapse_rate_local = methane_lapse_rate_k_per_km(temp_profile, pfull_local, kk, grav, Rdgas)
             if (lapse_rate_local <= ch4_tropopause_lapse_k_per_km) in_upper_cap = .true.
          else
             q_profile(kk) = q0
          end if
       elseif (.not. in_upper_cap) then
          q_profile(kk) = qsat_layer
          q_trop_out = qsat_layer
          p_trop_out = pfull_local(kk)
          lapse_rate_local = methane_lapse_rate_k_per_km(temp_profile, pfull_local, kk, grav, Rdgas)
          if (lapse_rate_local <= ch4_tropopause_lapse_k_per_km) in_upper_cap = .true.
       else
          q_profile(kk) = q_trop_out
       end if
    end do

    do kk = 1, nlay
       vmr_profile(kk) = methane_vmr_from_q(q_profile(kk), Rdgas, TITAN_RVGAS)
    end do

    inventory_pa = methane_inventory_from_profile(q_profile, phalf_local)
    q0_out = q0

  end subroutine build_simplified_ch4_profile_for_rh

  subroutine build_simplified_ch4_profile(temp_profile, tsurf_local, pfull_local, phalf_local, psurf_local, &
                                          target_surface_rh, inventory_limit_pa, q_profile, vmr_profile, &
                                          effective_surface_rh, inventory_pa, q0_out, q_trop_out, p_trop_out)

    real, intent(in), dimension(nlay)   :: temp_profile, pfull_local
    real, intent(in), dimension(nlay+1) :: phalf_local
    real, intent(in)                    :: tsurf_local, psurf_local, target_surface_rh, inventory_limit_pa
    real, intent(out), dimension(nlay)  :: q_profile, vmr_profile
    real, intent(out)                   :: effective_surface_rh, inventory_pa, q0_out, q_trop_out, p_trop_out
    real, dimension(nlay)               :: q_trial, vmr_trial
    real                                :: rh_low, rh_high, rh_mid
    real                                :: inventory_trial, q0_trial, q_trop_trial, p_trop_trial
    integer                             :: iter

    effective_surface_rh = max(0.0, min(1.0, target_surface_rh))

    call build_simplified_ch4_profile_for_rh(temp_profile, tsurf_local, pfull_local, phalf_local, psurf_local, &
                                             effective_surface_rh, q_profile, vmr_profile, inventory_pa, q0_out, &
                                             q_trop_out, p_trop_out)

    if (inventory_limit_pa <= 0.0 .or. inventory_pa <= inventory_limit_pa) return

    rh_low = 0.0
    rh_high = effective_surface_rh
    do iter = 1, 40
       rh_mid = 0.5 * (rh_low + rh_high)
       call build_simplified_ch4_profile_for_rh(temp_profile, tsurf_local, pfull_local, phalf_local, psurf_local, &
                                                rh_mid, q_trial, vmr_trial, inventory_trial, q0_trial, q_trop_trial, &
                                                p_trop_trial)
       if (inventory_trial > inventory_limit_pa) then
          rh_high = rh_mid
       else
          rh_low = rh_mid
       end if
    end do

    effective_surface_rh = rh_low
    call build_simplified_ch4_profile_for_rh(temp_profile, tsurf_local, pfull_local, phalf_local, psurf_local, &
                                             effective_surface_rh, q_profile, vmr_profile, inventory_pa, q0_out, &
                                             q_trop_out, p_trop_out)

  end subroutine build_simplified_ch4_profile

  subroutine apply_selected_convection(temp_profile, pfull_local, phalf_local, zfull_local, q_profile, &
                                       cape_out, cin_out, convflag_out, p_lcl_out, p_lfc_out, p_lzb_out, &
                                       max_delta_t_out)

    real, intent(inout), dimension(nlay) :: temp_profile
    real, intent(in), dimension(nlay)    :: pfull_local, zfull_local, q_profile
    real, intent(in), dimension(nlay+1)  :: phalf_local
    real, intent(out)                    :: cape_out, cin_out, p_lcl_out, p_lfc_out, p_lzb_out
    real, intent(out)                    :: max_delta_t_out
    integer, intent(out)                 :: convflag_out
    real, dimension(nlay)                :: temp_before, temp_after_moist

    temp_before = temp_profile
    cape_out = 0.0
    cin_out = 0.0
    p_lcl_out = 0.0
    p_lfc_out = 0.0
    p_lzb_out = 0.0
    convflag_out = 0

    select case (trim(convective_scheme))
    case ('none')
       continue
    case ('dry')
       call apply_dry_convective_adjustment(temp_profile, phalf_local, zfull_local)
       if (maxval(abs(temp_profile - temp_before)) > 1.0e-8) convflag_out = 1
    case ('moist_bm')
       call apply_moist_bm_adjustment(temp_profile, pfull_local, phalf_local, q_profile, cape_out, cin_out, &
                                      convflag_out, p_lcl_out, p_lfc_out, p_lzb_out)
       temp_after_moist = temp_profile
       if (use_dry_convective_fallback) then
          call apply_dry_convective_adjustment(temp_profile, phalf_local, zfull_local)
          if (maxval(abs(temp_profile - temp_after_moist)) > 1.0e-8) then
             if (convflag_out == 2) then
                convflag_out = 3
             else
                convflag_out = 1
             end if
          end if
       end if
    end select

    max_delta_t_out = maxval(abs(temp_profile - temp_before))

  end subroutine apply_selected_convection

  subroutine apply_dry_convective_adjustment(temp_profile, phalf_local, zfull_local)

    real, intent(inout), dimension(nlay) :: temp_profile
    real, intent(in), dimension(nlay+1)  :: phalf_local
    real, intent(in), dimension(nlay)    :: zfull_local
    real                                 :: lapse_limit, dtdz_local, dz_local
    real                                 :: dp_upper, dp_lower, temp_upper_new, temp_lower_new
    integer                              :: k, pass
    logical                              :: changed

    lapse_limit = -grav / cp_air

    do pass = 1, nlay
       changed = .false.
       do k = 1, nlay-1
          dz_local = zfull_local(k) - zfull_local(k+1)
          if (dz_local <= 0.0) cycle

          dtdz_local = (temp_profile(k+1) - temp_profile(k)) / (zfull_local(k+1) - zfull_local(k))
          if (dtdz_local <= lapse_limit) then
             dp_upper = max(phalf_local(k+1) - phalf_local(k), 1.0e-10)
             dp_lower = max(phalf_local(k+2) - phalf_local(k+1), 1.0e-10)
             temp_lower_new = (dp_upper*temp_profile(k) + dp_lower*temp_profile(k+1) - &
                               dp_upper*lapse_limit*dz_local) / (dp_upper + dp_lower)
             temp_upper_new = temp_lower_new + lapse_limit*dz_local
             if (abs(temp_upper_new - temp_profile(k)) > 1.0e-8 .or. &
                 abs(temp_lower_new - temp_profile(k+1)) > 1.0e-8) changed = .true.
             temp_profile(k) = temp_upper_new
             temp_profile(k+1) = temp_lower_new
          end if
       end do
       if (.not. changed) exit
    end do

  end subroutine apply_dry_convective_adjustment

  subroutine apply_moist_bm_adjustment(temp_profile, pfull_local, phalf_local, q_profile, cape_out, cin_out, &
                                       convflag_out, p_lcl_out, p_lfc_out, p_lzb_out)

    real, intent(inout), dimension(nlay) :: temp_profile
    real, intent(in), dimension(nlay)    :: pfull_local, q_profile
    real, intent(in), dimension(nlay+1)  :: phalf_local
    real, intent(out)                    :: cape_out, cin_out, p_lcl_out, p_lfc_out, p_lzb_out
    integer, intent(out)                 :: convflag_out
    real, dimension(nlay)                :: parcel_temp, tref, qref, q_work
    real                                 :: relax_factor, t_shift, dp_sum, energy_sum, qsat_local
    integer                              :: k, k_lcl, k_lfc, k_lzb

    convflag_out = 0
    cape_out = 0.0
    cin_out = 0.0
    p_lcl_out = 0.0
    p_lfc_out = 0.0
    p_lzb_out = 0.0
    parcel_temp = temp_profile
    tref = temp_profile
    qref = 0.0

    do k = 1, nlay
       qsat_local = analytic_ch4_qsat(temp_profile(k), pfull_local(k), Rdgas, TITAN_RVGAS)
       q_work(k) = max(0.0, min(q_profile(k), qsat_local))
    end do

    if (maxval(q_work) <= 1.0e-10) return

    call diagnose_methane_parcel(temp_profile, pfull_local, phalf_local, q_work, parcel_temp, cape_out, cin_out, &
                                 p_lcl_out, p_lfc_out, p_lzb_out, k_lcl, k_lfc, k_lzb)

    if (k_lzb <= 0 .or. cape_out <= moist_conv_min_cape) return

    do k = k_lzb, nlay
       tref(k) = parcel_temp(k)
       qsat_local = analytic_ch4_qsat(parcel_temp(k), pfull_local(k), Rdgas, TITAN_RVGAS)
       qref(k) = min(q_work(k), moist_conv_rhbm * qsat_local)
    end do

    energy_sum = 0.0
    dp_sum = 0.0
    do k = k_lzb, nlay
       energy_sum = energy_sum + ((tref(k) - temp_profile(k)) + TITAN_HLV/cp_air * (qref(k) - q_work(k))) * &
                                  max(phalf_local(k+1) - phalf_local(k), 0.0)
       dp_sum = dp_sum + max(phalf_local(k+1) - phalf_local(k), 0.0)
    end do

    if (dp_sum <= 0.0) return

    t_shift = -energy_sum / dp_sum
    tref(k_lzb:nlay) = tref(k_lzb:nlay) + t_shift

    relax_factor = min(1.0, time_step / max(moist_conv_tau, 1.0e-10))
    temp_profile(k_lzb:nlay) = temp_profile(k_lzb:nlay) + relax_factor * (tref(k_lzb:nlay) - temp_profile(k_lzb:nlay))
    convflag_out = 2

  end subroutine apply_moist_bm_adjustment

  subroutine diagnose_methane_parcel(temp_profile, pfull_local, phalf_local, q_profile, parcel_temp, cape_out, &
                                     cin_out, p_lcl_out, p_lfc_out, p_lzb_out, k_lcl, k_lfc, k_lzb)

    real, intent(in), dimension(nlay)    :: temp_profile, pfull_local, q_profile
    real, intent(in), dimension(nlay+1)  :: phalf_local
    real, intent(out), dimension(nlay)   :: parcel_temp
    real, intent(out)                    :: cape_out, cin_out, p_lcl_out, p_lfc_out, p_lzb_out
    integer, intent(out)                 :: k_lcl, k_lfc, k_lzb
    real, dimension(nlay)                :: parcel_r
    real                                 :: pref_local, r_surface, r_sat_surface, theta_parcel
    real                                 :: env_tv, parcel_tv, dlnp, t_lcl
    logical                              :: saturated, lcl_found, cape_found, lzb_found
    integer                              :: k

    parcel_temp = temp_profile
    parcel_r = 0.0
    cape_out = 0.0
    cin_out = 0.0
    p_lcl_out = 0.0
    p_lfc_out = 0.0
    p_lzb_out = 0.0
    k_lcl = 0
    k_lfc = 0
    k_lzb = 0

    pref_local = max(phalf_local(nlay+1), 1.0)
    r_surface = specific_to_mixing_ratio(max(q_profile(nlay), 0.0))
    if (r_surface <= 1.0e-12) return

    r_sat_surface = saturation_mixing_ratio(temp_profile(nlay), pfull_local(nlay))
    r_surface = min(r_surface, r_sat_surface)

    parcel_temp(nlay) = temp_profile(nlay)
    parcel_r(nlay) = r_surface
    saturated = r_surface >= (1.0 - 1.0e-6) * r_sat_surface
    lcl_found = saturated
    cape_found = .false.
    lzb_found = .false.

    if (lcl_found) then
       k_lcl = nlay
       p_lcl_out = pfull_local(nlay)
    end if

    do k = nlay, 1, -1
       env_tv = virtual_temp_from_q(temp_profile(k), q_profile(k))
       parcel_tv = virtual_temp_from_q(parcel_temp(k), mixing_to_specific_humidity(parcel_r(k)))
       dlnp = log(max(phalf_local(k+1), phalf_local(k) + 1.0e-10) / max(phalf_local(k), 1.0e-10))

       if (parcel_tv > env_tv) then
          cape_out = cape_out + Rdgas * (parcel_tv - env_tv) * dlnp
          if (.not. cape_found) then
             cape_found = .true.
             k_lfc = k
             p_lfc_out = pfull_local(k)
          end if
       else
          if (.not. cape_found) then
             cin_out = cin_out + Rdgas * (env_tv - parcel_tv) * dlnp
          else if (.not. lzb_found) then
             k_lzb = min(nlay, k + 1)
             p_lzb_out = pfull_local(k_lzb)
             lzb_found = .true.
             exit
          end if
       end if

       if (k == 1) exit

       if (.not. saturated) then
          theta_parcel = parcel_temp(k) * (pref_local / max(pfull_local(k), 1.0))**(Rdgas/cp_air)
          if (parcel_r(k) <= saturation_mixing_ratio(theta_parcel * (pfull_local(k-1)/pref_local)**(Rdgas/cp_air), &
                                                     pfull_local(k-1))) then
             parcel_temp(k-1) = theta_parcel * (pfull_local(k-1)/pref_local)**(Rdgas/cp_air)
             parcel_r(k-1) = parcel_r(k)
          else
             call find_lcl_between_levels(theta_parcel, parcel_r(k), pref_local, pfull_local(k), pfull_local(k-1), &
                                          p_lcl_out, t_lcl)
             saturated = .true.
             lcl_found = .true.
             k_lcl = k - 1
             call moist_ascent_step(t_lcl, p_lcl_out, pfull_local(k-1), parcel_temp(k-1), parcel_r(k-1))
          end if
       else
          call moist_ascent_step(parcel_temp(k), pfull_local(k), pfull_local(k-1), parcel_temp(k-1), parcel_r(k-1))
       end if
    end do

    if (.not. lcl_found) then
      cape_out = 0.0
      cin_out = 0.0
      k_lcl = 0
      k_lfc = 0
      k_lzb = 0
      p_lcl_out = 0.0
      p_lfc_out = 0.0
      p_lzb_out = 0.0
      parcel_temp = temp_profile
      return
    end if

    if (.not. cape_found) then
      cape_out = 0.0
      k_lfc = 0
      k_lzb = 0
      p_lfc_out = 0.0
      p_lzb_out = 0.0
      parcel_temp = temp_profile
      return
    end if

    if (.not. lzb_found) then
       k_lzb = 1
       p_lzb_out = pfull_local(1)
    end if

  end subroutine diagnose_methane_parcel

  subroutine find_lcl_between_levels(theta_parcel, parcel_r, pref_local, p_bottom, p_top, p_lcl, t_lcl)

    real, intent(in)  :: theta_parcel, parcel_r, pref_local, p_bottom, p_top
    real, intent(out) :: p_lcl, t_lcl
    real              :: lnp_low, lnp_high, lnp_mid, p_mid, t_mid, r_sat_mid
    integer           :: iter

    lnp_low = log(max(p_top, 1.0e-10))
    lnp_high = log(max(p_bottom, 1.0e-10))

    do iter = 1, 40
       lnp_mid = 0.5 * (lnp_low + lnp_high)
       p_mid = exp(lnp_mid)
       t_mid = theta_parcel * (p_mid/pref_local)**(Rdgas/cp_air)
       r_sat_mid = saturation_mixing_ratio(t_mid, p_mid)
       if (parcel_r > r_sat_mid) then
          lnp_low = lnp_mid
       else
          lnp_high = lnp_mid
       end if
    end do

    p_lcl = exp(0.5 * (lnp_low + lnp_high))
    t_lcl = theta_parcel * (p_lcl/pref_local)**(Rdgas/cp_air)

  end subroutine find_lcl_between_levels

  subroutine moist_ascent_step(t_start, p_start, p_end, t_end, r_end)

    real, intent(in)  :: t_start, p_start, p_end
    real, intent(out) :: t_end, r_end
    real              :: t_mid, p_mid, r_start, r_mid, dlnp
    real              :: a, b, dtdlnp

    dlnp = log(max(p_end, 1.0e-10) / max(p_start, 1.0e-10))
    p_mid = exp(0.5 * (log(max(p_start, 1.0e-10)) + log(max(p_end, 1.0e-10))))
    r_start = saturation_mixing_ratio(t_start, p_start)

    a = (Rdgas/cp_air) * t_start + (TITAN_HLV/cp_air) * r_start
    b = (TITAN_HLV**2) * r_start / (cp_air * TITAN_RVGAS * max(t_start, 40.0)**2)
    dtdlnp = a / (1.0 + b)

    t_mid = t_start + 0.5 * dlnp * dtdlnp
    r_mid = saturation_mixing_ratio(t_mid, p_mid)

    a = (Rdgas/cp_air) * t_mid + (TITAN_HLV/cp_air) * r_mid
    b = (TITAN_HLV**2) * r_mid / (cp_air * TITAN_RVGAS * max(t_mid, 40.0)**2)
    dtdlnp = a / (1.0 + b)

    t_end = t_start + dlnp * dtdlnp
    r_end = saturation_mixing_ratio(t_end, p_end)

  end subroutine moist_ascent_step

  real function saturation_mixing_ratio(temp, pressure)

    real, intent(in) :: temp, pressure
    real             :: qsat_local

    qsat_local = analytic_ch4_qsat(temp, pressure, Rdgas, TITAN_RVGAS)
    saturation_mixing_ratio = specific_to_mixing_ratio(qsat_local)

  end function saturation_mixing_ratio

  real function specific_to_mixing_ratio(qval)

    real, intent(in) :: qval

    specific_to_mixing_ratio = qval / max(1.0 - qval, 1.0e-10)

  end function specific_to_mixing_ratio

  real function mixing_to_specific_humidity(rval)

    real, intent(in) :: rval

    mixing_to_specific_humidity = rval / (1.0 + rval)

  end function mixing_to_specific_humidity

  real function virtual_temp_from_q(temp, qval)

    real, intent(in) :: temp, qval

    virtual_temp_from_q = temp * (1.0 + qval * (TITAN_RVGAS/Rdgas - 1.0))

  end function virtual_temp_from_q

END PROGRAM
