MODULE planetary_radiation_driver_mod
! ========================================================================================
! GENERALIZED RADIATIVE TRANSFER MODEL FOR PLANETARY GCMs
! Copyright (c) Juan M. Lora
! ========================================================================================
! Written by J.M. Lora
! ========================================================================================
! February 2021 - version 0.9
!   - Updates by Nick Lombardo
!   - Can now choose between using a climatology or prescribed abundance profiles
!   - Allows for spatially and temporally varying trace gases and aerosol properties
!       . trace gas climatology (TGC) in netCDF format in INPUT/DATA
!       . TGC must be on integer latitude grid, since latitudes are used as indices,
!           also must use integer ls value
!       . TGC includes saturate subroutine, which saturates molecules at low temperatures
!           constistent with calulated temperature profile
!       . haze climatology (HC) does not need to be on a special grid, interpolating
!           is done in haze.F90 module.
!       . haze_prec_file needs to be included even if HC is used for both vis and ir bands
!           as it still contains haze asymmetry and single-scattering albedos
!       . CIA is NOT functional with TGC yet, since N2, CH4, and H2 are not included in TGC
!   - Several bugs fixed throughout PR, PRD
!
!
! November 2017 - version 0.1 (untested 11/4/2017)
! 	- Allows for arbitrary radiatively active gases to be used, input through namelist
!		. k-coefficients must be provided as INPUT/DATA files
!		. ALL gas k-coefficients must be on the same temperature and pressure
!			grid, and share the same Gaussian points
!		. Number of temperatures, pressures, and Gaussian points used must be 
!			set in namelist (ktemps, kpress, kgauss)
!	- Allows for CIA from arbitrary pairs, input through namelist
!		. Transmission fits must be provided as INPUT/DATA files
!		. ALL fits must be on same temperature grid
!		. Number of fit points (tfits) and temperatures (ttemps) must be set in namelist
!	- Visible and IR bands can be set via namelist (k-coefficients and CIA 
!			parameters must match these bands)
!	- For any radiatively active gas (including CIA), it must either have a prescribed
!			profile, or be an advected tracer, and these must be set in namelist
!			with presc_gas and/or trace_gas
!		. For prescribed profiles, all must be on the SAME pressure grid
!		. Number of prescribed pressure levels for profiles must be set in 
!			namelist (ppress)
! ========================================================================================
 
        use planetary_radiation_mod,    only: setspv, setspi, get_taukcoeff, get_tauCIA, &
                                              optc, sfluxv, sfluxi, error_mesg 
        use haze_mod,                   only: haze_init, hazeclim_init, interp_haze_ls, clim_haze_profile, haze_end 
        use read_clim_mod,              only: readclim_MOD, saturate
        ! use netcdf  ! removed: no nf90_ calls reached in prescribed mode

        implicit none
        private
        save

        public :: planetary_radiation_init, planetary_radiation, planetary_radiation_end
        public :: set_top_tlev_restart, get_top_tlev_restart, set_analytic_ch4_mode

!-----------------------------------------------------------------------------------------
        ! Namelist (defaults defined here):
!-----------------------------------------------------------------------------------------

        real        :: pi           = 3.14159
        real        :: grav         = 8.94    ! [Titan: 1.352  | Uranus: 8.94]
        real        :: Rdgas        = 3149.2  ! [Titan: 287.05 | Uranus: 3149.2]
        real        :: Rvgas        = 461.50
        real        :: Cp_air       = 8600.   ! [Titan: 1039.5 | Uranus: 8600.]
        real        :: Cp_vapor     = 225.    ! Methane
        real        :: FATAL        = 0.
        real        :: NO_TRACER    = 0.
        integer     :: tmpsize

	real                    ::  sma = 9.58                 ! Orbital semi-major axis (AU)
        character(len=14)       ::  haze_presc_file            ! constant haze name
        character(len=14)       ::  haze_clim_file = 'haze.nc' ! clim haze name

        ! "Solar" band edges, defined in wavenumber (using wavelengths in nm)
        ! NOTE: original used the non-standard range syntax (/600:300:-50/);
        ! rewritten as explicit literals (identical values) so it compiles with
        ! stock gfortran.
  	real    :: bwnv(43) = 1.D6/[ 600.,550.,500.,450.,400.,350.,300., &
                                    250.,240.,230.,220.,210.,200.,190.,180.,170.,160., &
                                    150.,145.,140.,135.,130.,125.,120.,115.,110.,105.,100., &
                                    95.,90.,85.,80.,75.,70.,65.,60.,55.,50.,45.,40.,35.,30.,25. ]
        ! "Thermal infrared" band edges, defined in wavenumber
  	real    :: bwni(41) = [ 1.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500., &
                               550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.,1050., &
                               1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550., &
                               1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000. ]
  
        character(len=128)      :: solar_spec_file      = 'solar_spectrum_houghton.txt'
        character(len=128)      :: gas_profiles         = 'trace_gases.txt' ! get rid of this [CK]
        character(len=128)      :: rayleigh_file        = 'Rayleigh.txt' ! fix hard coding [CK]
 
        character(len=10), dimension(10) :: radgasv     = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10), dimension(10) :: radgasi     = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10), dimension(10) :: ciapairv    = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10), dimension(10) :: ciapairi    = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10), dimension(10) :: presc_gas   = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10), dimension(10) :: clim_gas    = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10)                :: haze_data   = 'presc'
        character(len=10), dimension(10) :: trace_gas   = (/'','','','','',&
                                                            '','','','',''/)
        character(len=10)                :: bkgnd_gas   = 'N2'
        character(len=10)                :: sphum_gas   = ''

        integer         :: ktemps       = 10    ! Number of temperatures for k-coeffs
        integer         :: kpress       = 20    ! Number of pressures for k-coeffs
        integer         :: kgauss       = 10    ! Number of Gaussian points for k-coeffs
        integer         :: ttemps       = 10    ! Number of temperatures for CIA transmission fits
        integer         :: tfits        = 3     ! Number of fit points for CIA transmission
        integer         :: ppress       = 100   ! Number of pressure levels in prescribed gas profiles
        logical         :: write_cia_pair_diagnostics = .false.

        namelist/planetary_radiation_nml/       sma, grav, Rdgas, Rvgas, Cp_air, Cp_vapor, &
                                                haze_presc_file, bwnv, bwni, &
                                                solar_spec_file, rayleigh_file, gas_profiles, &
                                                radgasv, radgasi, ciapairv, ciapairi, &
                                                presc_gas, trace_gas, clim_gas, bkgnd_gas, sphum_gas, &
                                                ktemps, kpress, kgauss, &
                                                ttemps, tfits, ppress, haze_clim_file, haze_data, &
                                                write_cia_pair_diagnostics

!-----------------------------------------------------------------------------------------
! "Visible" and "infrared" band variables:
!-----------------------------------------------------------------------------------------

	real, dimension(size(bwnv,1)-1)         ::      wnov, dwnv, tauray
	real, dimension(:), allocatable         ::      solarf
	real, dimension(size(bwni,1)-1)         ::      wnoi, dwni
	real, dimension(size(bwni,1)-1,8801)    ::      planckir ! CWK: originally 8501

!-----------------------------------------------------------------------------------------
! Correlated k-coefficient arrays:
!-----------------------------------------------------------------------------------------

	real, dimension(:,:,:,:,:), allocatable :: ckcv     ! "Visible" correlated k-coeffs
	real, dimension(:,:,:,:,:), allocatable :: ckci     ! "Infrared" correlated k-coeffs
	real, dimension(:), allocatable         :: ckc_temp ! Temps for k-coeffs
	real, dimension(:), allocatable         :: ckc_pres ! Pressures for k-coeffs
	real, dimension(:), allocatable         :: ckc_gwtv ! Gaussian weights for vis k-coeffs
	real, dimension(:), allocatable         :: ckc_gwti ! Gaussian weights for IR k-coeffs

!-----------------------------------------------------------------------------------------
! CIA pair transmission fit arrays:
!-----------------------------------------------------------------------------------------

	real, dimension(:,:,:,:), allocatable   :: ciav         ! "Visible" CIA fits
	real, dimension(:,:,:,:), allocatable   :: ciai         ! "Infrared" CIA fits
	real, dimension(:), allocatable         :: cia_temp     ! Temps for CIA fits

!-----------------------------------------------------------------------------------------
! Prescribed absorber profiles and radiatively active tracer indices
!-----------------------------------------------------------------------------------------

	real, dimension(:,:), allocatable       :: ref_gases
	real, dimension(:,:,:,:), allocatable   :: clim_data
	real, dimension(:, :, :), allocatable   :: cl_inst
	real, dimension(:), allocatable         :: clpress
	real, dimension(:), allocatable         :: clls            ! climatology ls (haze ls interp done in haze mod)
	real                                    :: cllsidx, cllswt ! climatology ls interpolation parameters
	real, dimension(:, :, :), allocatable   :: hazeclim
        integer, dimension(:), allocatable      :: tracer_indices
        integer                                 :: sphum_index
        integer, dimension(:), allocatable      :: clmolidx
	real, dimension(:, :), allocatable      :: tau_hazev, tau_hazei, ssa_hazev, ssa_hazei, g_hazev, g_hazei,&
                                                   tau_hazev_init, tau_hazei_init, dtau_hazev, dtau_hazei

!-----------------------------------------------------------------------------------------
! Define needed constants (not already in constants.F90)
!-----------------------------------------------------------------------------------------

	real, parameter         :: ubari        = 0.50
	real, parameter         :: Runiv        = 8.3144598
	real, parameter         :: mole_to_amg  = 0.02241397 ! amg m3 mol-1
	real, parameter         :: avo          = 6.02214086 ! Divided by 1e23!

!-----------------------------------------------------------------------------------------
! Other things:
!-----------------------------------------------------------------------------------------

        character(len=128), parameter :: version = '$ID: planetary_radiation.F90, v 0.1 $'
        character(len=128), parameter :: tagname = '$Name: testing $'

        integer         :: L_LAYERS, L_LEVELS, ww

	real            ::      missing_value           = -999.
	logical         ::      module_is_initialized   = .false.
	logical         ::      do_radgasv              = .false.
	logical         ::      do_radgasi              = .false.
	logical         ::      do_ciav                 = .false.
	logical         ::      do_ciai                 = .false.
	logical         ::      use_presc_gases         = .false.
        logical         ::      use_clim_gases          = .false.
        logical         ::      haze_or_clouds          = .false.
        logical         ::      use_haze_clim           = .false.
        logical         ::      use_top_tlev_restart    = .false.
        logical         ::      have_last_top_tlev      = .false.
        logical         ::      use_analytic_ch4_profile = .false.
        real, dimension(3) ::   top_tlev_restart        = 0.0
        real, dimension(3) ::   last_top_tlev           = 0.0

!-----------------------------------------------------------------------------------------
! Output ids:
!-----------------------------------------------------------------------------------------

        character(len=10), parameter    :: mod_name = 'radiation'

        integer                         :: id_sw_toa, id_lw_toa, id_sw_sfc, id_lwdn_sfc, id_lwup_sfc, &
                                           id_flux_sw, id_flux_lw, id_swdn, id_swup, id_lwdn, id_lwup, &
                                           id_tdt_rad, id_tdt_sw, id_tdt_lw, id_sw_wn, id_lw_wn, &
                                           id_sw_spec, id_lw_spec 

! ========================================================================================
        contains

subroutine set_top_tlev_restart(tlev_top, use_restart)

        real, intent(in), dimension(3) :: tlev_top
        logical, intent(in)            :: use_restart

        use_top_tlev_restart = use_restart
        if (use_restart) then
                top_tlev_restart = tlev_top
                last_top_tlev = tlev_top
                have_last_top_tlev = .true.
        else
                top_tlev_restart = 0.0
        end if

end subroutine set_top_tlev_restart

subroutine get_top_tlev_restart(tlev_top, have_top_tlev)

        real, intent(out), dimension(3) :: tlev_top
        logical, intent(out)            :: have_top_tlev

        tlev_top = last_top_tlev
        have_top_tlev = have_last_top_tlev

end subroutine get_top_tlev_restart

subroutine set_analytic_ch4_mode(use_analytic_ch4)

        logical, intent(in) :: use_analytic_ch4

        use_analytic_ch4_profile = use_analytic_ch4

end subroutine set_analytic_ch4_mode

subroutine apply_top_tlev_boundary(plev, tlev)

        real, intent(in), dimension(:)    :: plev
        real, intent(inout), dimension(:) :: tlev
        integer                            :: k

        if (use_top_tlev_restart) then
                tlev(1:3) = top_tlev_restart
        else
                do k = 3,1,-1
                        tlev(k) = tlev(k+1) + (tlev(k+1) - tlev(k+2)) * log(plev(k)/plev(k+1)) / log(plev(k+1)/plev(k+2))
                end do
        end if

        last_top_tlev = tlev(1:3)
        have_last_top_tlev = .true.

end subroutine apply_top_tlev_boundary
! ========================================================================================

SUBROUTINE planetary_radiation_init(npz, lat_grid, phalf, pfull)
        integer, intent(in)                     :: npz
        real, dimension(:), intent(in)          :: lat_grid
        real, intent(in), dimension(:)          :: phalf, pfull
        real, dimension(2*npz+3)                :: plev
        integer                                 :: ierr, io, unit, logunit
        integer                                 :: g, p, t, n, k
	real                                    :: check_var
	real, dimension(size(wnov,1))           :: check_wnov
	real, dimension(size(wnoi,1))           :: check_wnoi
	real, dimension(size(ckc_gwtv,1))       :: check_gwtv
	real, dimension(size(ckc_gwti,1))       :: check_gwti

!-----------------------------------------------------------------------------------------
! Read namelist and copy to logfile
!-----------------------------------------------------------------------------------------

        open(5, file = 'namelist')
        read(5, planetary_radiation_nml)
        close(5)

        if (.true.) then
           print *, 'Radiation: Using a semi-major axis of ',sma
           print *, 'Radiation: Stellar spectrum from: '//trim(solar_spec_file)
           print *, 'Radiation: Absorber profiles from: '//trim(gas_profiles)
           print *, 'Radiation: grav, Rdgas, Rvgas = ', grav, Rdgas, Rvgas
        end if

!-----------------------------------------------------------------------------------------
! Initialize layers & spectral intervals
!-----------------------------------------------------------------------------------------

        L_LAYERS = npz
        L_LEVELS = 2*npz + 3

        call setspv(wnov,dwnv,solarf,tauray,bwnv,solar_spec_file,rayleigh_file)
        call setspi(wnoi,dwni,planckir,bwni)

!-----------------------------------------------------------------------------------------
! Initialize k-coefficients
!-----------------------------------------------------------------------------------------
        do g = 1,size(radgasv,1)
                if (len(trim(radgasv(g))) .gt. 0) do_radgasv = .true.
        end do
        do g = 1,size(radgasi,1)
                if (len(trim(radgasi(g))) .gt. 0) do_radgasi = .true.
        end do

        ! CWK: Modified to handle the case where there are no radiative active gases in IR or vis
        if (do_radgasv .or. do_radgasi) then
                allocate (ckc_temp(ktemps))
                allocate (ckc_pres(kpress))
                allocate (ckc_gwtv(kgauss))
                allocate (ckc_gwti(kgauss))
        else
                allocate (ckc_temp(ktemps))
                allocate (ckc_pres(kpress))
                allocate (ckc_gwtv(kgauss))
                allocate (ckc_gwti(kgauss))
                ckc_temp(:) = 0.
                ckc_pres(:) = 0.
                ckc_gwtv(:) = 0.
                ckc_gwti    = (/ 1., 0., 0., 0., 0., 0., 0., 0., 0., 0. /)
        end if

        if (do_radgasv) then
                allocate ( ckcv( size(radgasv,1), kpress, ktemps, size(wnov,1), kgauss ), stat=io)
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for ckcv allocation', FATAL)
        end if

        if (do_radgasi) then
                allocate ( ckci( size(radgasi,1), kpress, ktemps, size(wnoi,1), kgauss ), stat=io)
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for ckci allocation', FATAL)
        end if
       
        if (do_radgasv) then
        do g = 1,size(radgasv,1)
                if (len(trim(radgasv(g))) .gt. 0) then
                        open(10,file='INPUT/DATA/ckc_'//trim(radgasv(g))//'vis.txt', status='old',iostat=io)
                        if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                                'planetary_radiation_init: Could not open vis ckc_'//trim(radgasv(g)), FATAL)

                        read(10,*) ! Read header text
                        read(10,*) check_wnov(:)
                        if ( .not. (all(abs(check_wnov - wnov) < 0.1))) then
                                call error_mesg('planetary_radiation_driver_mod', &
                                'planetary_radiation_init: CKCv bands do not match model bands', FATAL)
                        end if
                        if (g.eq.1) then
                                read(10,*) ckc_gwtv(:)
                        else
                                read(10,*) check_gwtv(:)
                                if ( .not. (all(check_gwtv .eq. ckc_gwtv))) then
                                        call error_mesg('planetary_radiation_driver_mod', &
                                                'planetary_radiation_init: CKCv Gaussian grids differ', FATAL)
                                end if
                        end if

                        do p = 1,kpress
                                if (g.eq.1) then
                                        read(10,*) ckc_pres(p)
                                else
                                        read(10,*) check_var
                                        if (check_var .ne. ckc_pres(p)) then
                                                call error_mesg('planetary_radiation_driver_mod', &
                                                'planetary_radiation_init: CKCv pressure grids differ', FATAL)
                                        end if
                                end if
                                do t = 1,ktemps
                                        if (g.eq.1) then
                                                read(10,*) ckc_temp(t)
                                        else
                                                read(10,*) check_var
                                                if (check_var .ne. ckc_temp(t)) then
                                                        call error_mesg('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: CKCv temp grids differ', FATAL)
                                                end if
                                        end if
                                        do n = 1,size(wnov,1)
                                                read(10,*) ckcv(g,p,t,n,:)
                                        end do
                                end do
                        end do
                        close(10)
            
                     if (.true.) then
                         print *, 'Radiation:  Using k-coeffs for visible '//trim(radgasv(g))
                     end if
                end if
        end do
        end if

        if (do_radgasi) then
                do g = 1,size(radgasi,1)
                        if (len(trim(radgasi(g))) .gt. 0) then
                                open(10,file='INPUT/DATA/ckc_'//trim(radgasi(g))//'ir.txt',status='old',iostat=io)
                                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                                'planetary_radiation_init: Could not open IR ckc_'//trim(radgasi(g)), FATAL)
                                
                                read(10,*) ! Read header text
                                read(10,*) check_wnoi(:) ! check_wnoi is ckc grid
                                if( .not. (all(abs(check_wnoi - wnoi) < 0.1))) then
                                    print *, check_wnoi
                                    print *, wnoi
                                    print *, 'ERROR - CKC AND MODEL BANDS DO NOT MATCH - PRD_L375'
                                    call error_mesg('planetary_radiation_driver_mod', &
                                            'planetary_radiation_init: CKCi bands do not match model bands',FATAL)
                                end if
                    
                                if (g.eq.1) then
                                        read(10,*) ckc_gwti(:)
                                else
                                        read(10,*) check_gwti(:)
                                        if ( .not. (all(check_gwti .eq. ckc_gwti))) then
                                                call error_mesg('planetary_radiation_driver_mod', &
                                                'planetary_radiation_init: CKCi Gaussian grids differ', FATAL)
                                        end if
                                end if

                                do p = 1,kpress
                                        if (g.eq.1) then
                                                read(10,*) ckc_pres(p)
                                        else
                                                read(10,*) check_var
                                                if (check_var .ne. ckc_pres(p)) then
                                                        call error_mesg('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: CKCi pressure grids differ', FATAL)
                                                end if
                                        end if
                                        do t = 1,ktemps
                                                if (g.eq.1) then
                                                        read(10,*) ckc_temp(t)
                                                else
                                                        read(10,*) check_var
                                                        if (check_var .ne. ckc_temp(t)) then
                                                                call error_mesg('planetary_radiation_driver_mod', &
                                                                'planetary_radiation_init: CKCi temp grids differ', FATAL)
                                                        end if
                                                end if
                                                do n = 1,size(wnoi,1)
                                                        read(10,*) ckci(g,p,t,n,:)
                                                end do
                                        end do
                                end do
                                close(10)
                                !NLEDIT - multiprocessor, uncomment for final
                                !if (mpp_pe() == mpp_root_pe()) then
                                if (.true.) then
                                        print *, 'Radiation:  Using k-coeffs for IR '//trim(radgasi(g))
                                end if
                        end if
                end do
        end if

!-----------------------------------------------------------------------------------------
! Initialize CIA transmission fit parameters
!-----------------------------------------------------------------------------------------
        do g = 1,size(ciapairv,1)
                if (len(trim(ciapairv(g))) .gt. 0) do_ciav = .true.
        end do
        do g = 1,size(ciapairi,1)
                if (len(trim(ciapairi(g))) .gt. 0) do_ciai = .true.
        end do
        if (do_ciav .or. do_ciai) then
                allocate (cia_temp(ttemps))
        else
                print *, 'CIA turned off - PRD_L416'
        end if
        if (do_ciav) then
                allocate ( ciav( size(ciapairv,1), ttemps, size(wnov,1), tfits * 2 ), stat = io) 
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for ciav allocation', FATAL)
        end if
        if (do_ciai) then
                allocate ( ciai( size(ciapairi,1), ttemps, size(wnoi,1), tfits * 2 ), stat = io)
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for ciai allocation', FATAL)
        end if
        if (do_ciav) then
                do g = 1,size(ciapairv,1)
                        if (len(trim(ciapairv(g))) .gt. 0) then
                                open(10,file='INPUT/DATA/trans_'//trim(ciapairv(g))//'.txt', status='old',iostat=io)
                                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: Could not open vis CIA '//trim(ciapairv(g)), FATAL)
                                read(10,*) ! Read header text
                                read(10,*) check_wnov(:)
                                if( .not. (all(abs(check_wnov - wnov) < 0.1))) then
                                        call error_mesg('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: CIAv bands do not match model bands',FATAL)
                                end if
                                do t = 1,ttemps
                                        if (g.eq.1) then
                                                read(10,*) cia_temp(t)
                                        else
                                                read(10,*) check_var
                                                if (check_var .ne. cia_temp(t)) then
                                                        call error_mesg('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: CIAv temp grids differ', FATAL)
                                                end if
                                        end if
                                        do n = 1,size(wnov,1)
                                                read(10,*) ciav(g,t,n,:)
                                        end do
                                end do
                                close(10)
                                !NLEDIT - multiprocessor, uncomment for final
                                !if (mpp_pe() == mpp_root_pe()) then
                                if (.true.) then
                                        print *, 'Using visible CIA for '//trim(ciapairv(g))
                                end if
                        end if
                end do
        end if
        if (do_ciai) then
                do g = 1,size(ciapairi,1)
                        if (len(trim(ciapairi(g))) .gt. 0) then
                                open(10,file='INPUT/DATA/trans_'//trim(ciapairi(g))//'.txt', status='old',iostat=io)
                                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: Could not open IR CIA '//trim(ciapairi(g)), FATAL)
                                read(10,*) ! Read header text
                                read(10,*) check_wnoi(:)
                                if( .not. (all(abs(check_wnoi - wnoi) < 0.1))) then
                                        call error_mesg('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: CIAi bands do not match model bands',FATAL)
                                end if
                                do t = 1,ttemps
                                        if (g.eq.1) then
                                                read(10,*) cia_temp(t)
                                        else
                                                read(10,*) check_var
                                                if (check_var .ne. cia_temp(t)) then
                                                        print *, check_var, cia_temp(t)
                                                        call error_mesg('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: CIAi temp grids differ', FATAL)
                                                end if
                                        end if
                                        do n = 1,size(wnoi,1)
                                                read(10,*) ciai(g,t,n,:)
                                        end do
                                end do
                                close(10)
                                !NLEDIT - multiprocessor, uncomment
                                !if (mpp_pe() == mpp_root_pe()) then
                                if (.true.) then
                                        print *, 'Radiation:  Using IR CIA for '//trim(ciapairi(g))
                                end if
                        end if
                end do
        end if
!-----------------------------------------------------------------------------------------
! Ensure there are no problematic gas absorber assignments
!-----------------------------------------------------------------------------------------
        if (do_radgasv .or. do_radgasi .or. do_ciav .or. do_ciai) then
                !NLEDIT - multiprocessor, uncomment
                !if (mpp_pe() .eq. mpp_root_pe()) then
                if (.true.) then
                        if (len(trim(bkgnd_gas)) .gt. 0) then
                                print *, 'Radiation:  Radiatively active background gas: '//trim(bkgnd_gas)
                        end if
                end if
                !if (mpp_pe() .eq. mpp_root_pe()) then - NLEDIT
                if (.true.) then
                        if (len(trim(sphum_gas)) .gt. 0) then
                                print *, 'Radiatively active gas from specific humidity tracer: '//trim(sphum_gas)
                        end if
                end if
                do g = 1,size(presc_gas)
                        if (len(trim(presc_gas(g))) .gt. 0) then
                                if (presc_gas(g) .eq. bkgnd_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(presc_gas(g))// &
                                        ' cannot be both prescribed and the background gas', FATAL)
                                elseif (presc_gas(g) .eq. sphum_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(presc_gas(g))// &
                                        ' cannot be both prescribed and the specific humidity', FATAL)
                                end if
                                do n = 1,size(clim_gas)
                                        if (len(trim(clim_gas(n))) .gt. 0) then
                                                if (presc_gas(g) .eq. clim_gas(n)) then
                                                        call error_mesg ('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: '//trim(presc_gas(g))// &
                                                        ' cannot be both prescribed and climatology', FATAL)
                                                end if
                                        end if
                                end do
                                do n = 1,size(trace_gas)
                                        if (len(trim(trace_gas(n))) .gt. 0) then
                                                if (presc_gas(g) .eq. trace_gas(n)) then
                                                        call error_mesg ('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: '//trim(presc_gas(g))// &
                                                        ' cannot be both prescribed and an advected tracer', FATAL)
                                                end if
                                        end if
                                end do
                        end if
                end do
                do g = 1,size(clim_gas)
                        if (len(trim(clim_gas(g))) .gt. 0) then
                                if (clim_gas(g) .eq. bkgnd_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(clim_gas(g))// &
                                        ' cannot be both climatology and the background gas', FATAL)
                                elseif (clim_gas(g) .eq. sphum_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(clim_gas(g))// &
                                        ' cannot be both climatology and the specific humidity', FATAL)
                                end if
                                do n = 1,size(trace_gas)
                                        if (len(trim(trace_gas(n))) .gt. 0) then
                                                if (clim_gas(g) .eq. trace_gas(n)) then
                                                        call error_mesg ('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: '//trim(clim_gas(g))// &
                                                        ' cannot be both climatology and an advected tracer', FATAL)
                                                end if
                                        end if
                                end do
                        end if
                end do
                do g = 1,size(trace_gas)
                        if (len(trim(trace_gas(g))) .gt. 0) then
                                if (trace_gas(n) .eq. bkgnd_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(trace_gas(n))// &
                                        ' cannot be both an advected tracer and the background gas', FATAL)
                                elseif (trace_gas(n) .eq. sphum_gas) then
                                        call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: '//trim(trace_gas(n))// &
                                        ' cannot be both an advected tracer and the spec humidity', FATAL)
                                end if
                        end if
                end do
        else
                if (.true.) print *, 'Radiation:  No radiatively active gases selected!'
        end if
!-----------------------------------------------------------------------------------------
! Initialize the vertical profiles for the prescribed absorbers
!-----------------------------------------------------------------------------------------
        if (do_radgasv .or. do_radgasi .or. do_ciav .or. do_ciai) then           
                allocate ( ref_gases( size(presc_gas,1)+1, ppress ), stat = io)
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for prescribed gas profiles', FATAL)
                do g = 1,size(presc_gas,1)
                        if (len(trim(presc_gas(g))) .gt. 0) then
                                if (use_analytic_ch4_profile .and. trim(presc_gas(g)) .eq. 'CH4') then
                                        ref_gases(g+1, :) = 0.0
                                        cycle
                                end if
                                open(10,file='INPUT/DATA/profile_'//trim(presc_gas(g))//'.txt', status='old',iostat=io)
                                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                                        'planetary_radiation_init: Could not open profile '//trim(presc_gas(g)), FATAL)
                                do p = 1,ppress
                                        if (.not. use_presc_gases) then
                                                read(10,*) ref_gases(1,p), ref_gases(g+1,p)
                                        else
                                                read(10,*) check_var, ref_gases(g+1,p)
                                                if (check_var .ne. ref_gases(1,p)) then
                                                        call error_mesg('planetary_radiation_driver_mod', &
                                                        'planetary_radiation_init: prescribed profile grids differ', FATAL)
                                                end if
                                        end if
                                end do
                                use_presc_gases = .true.
                        end if
                        close(10)
                end do

                ! Read climatology netCDF here - internal to init
                if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: Not enough memory for climatology gas data', FATAL)
                do g = 1,size(clim_gas,1)
                        if (len(trim(clim_gas(g))) .gt. 0) then
                                !call/read netcdf climatology, only need to do once
                                call readclim_MOD(clim_data, clpress, lat_grid, clls, clmolidx, clim_gas)
                                ALLOCATE( cl_inst(SIZE(clim_data(:, 1, 1, 1)), SIZE(clim_data(1, :, 1, 1)), SIZE(clim_data(1, 1, 1, :))) )
                                print *, 'Radiation:  Trace species climatology loaded'
                                use_clim_gases = .true.
                                exit
                        end if
                end do
        end if

!-----------------------------------------------------------------------------------------
! Get the indices of specific humidity and any other radiatively active advected tracers
!-----------------------------------------------------------------------------------------
        if (len(trim(sphum_gas)) .gt. 0) then
                !sphum_index = get_tracer_index(MODEL_ATMOS,'sphum')
                sphum_index = NO_TRACER
                if (sphum_index .eq. NO_TRACER) then
                        call error_mesg ('planetary_radiation_driver_mod', &
                        'planetary_radiation_init: specific humidity is not a tracer', FATAL)
                end if
        else
                sphum_index = NO_TRACER
        end if

        allocate ( tracer_indices( size(trace_gas,1) ), stat = io)
        if (io .ne. 0) call error_mesg ('planetary_radiation_driver_mod', &
                'planetary_radiation_init: Not enough memory for tracer index allocation', FATAL)
        do g = 1,size(trace_gas,1)
                if (len(trim(trace_gas(g))) .gt. 0) then
                        !tracer_indices(g) = get_tracer_index(MODEL_ATMOS,trim(trace_gas(g)))
                        tracer_indices(g) = NO_TRACER
                        if (tracer_indices(g) .eq. NO_TRACER) then
                                call error_mesg('planetary_radiation_driver_mod', &
                                'planetary_radiation_init: '//trim(trace_gas(g))//'is not a tracer',FATAL)
                        end if
                else
                        tracer_indices(g) = NO_TRACER
                end if
        end do

!-----------------------------------------------------------------------------------------
! Initialize Haze
!-----------------------------------------------------------------------------------------
        do k = 1,L_LAYERS
                plev(2*k+1) = phalf(k)
                plev(2*k+2) = pfull(k)
        end do
               
        if (plev(3) .eq. 0.0) plev(3) = plev(4)/10.

        plev(2) = plev(3)*plev(3)/(plev(4))
        plev(1) = plev(2)*plev(2)/plev(3)
        plev(L_LEVELS) = phalf( L_LAYERS+1)


        ALLOCATE(tau_hazev(size(plev), size(wnov)))
        ALLOCATE(tau_hazev_init(size(plev), size(wnov)))
        ALLOCATE(ssa_hazev(size(plev), size(wnov)))
        ALLOCATE(g_hazev(size(plev), size(wnov)))
        ALLOCATE(tau_hazei(size(plev), size(wnoi)))
        ALLOCATE(tau_hazei_init(size(plev), size(wnoi)))
        ALLOCATE(ssa_hazei(size(plev), size(wnoi)))
        ALLOCATE(g_hazei(size(plev), size(wnoi)))

        ALLOCATE(dtau_hazev(size(plev), size(wnov)))
        ALLOCATE(dtau_hazei(size(plev), size(wnoi)))

        if (TRIM(haze_data) .eq. 'clim') then
            call hazeclim_init(wnoi, wnov, lat_grid*180./pi, plev, haze_clim_file)
            call haze_init(tau_hazev_init, g_hazev, ssa_hazev, haze_presc_file, plev, wnov, 'v', doclimin = 1)
            call haze_init(tau_hazei_init, g_hazei, ssa_hazei, haze_presc_file, plev, wnoi, 'i', doclimin = 1)
            haze_or_clouds = .true.
        elseif (TRIM(haze_data) .eq. 'presc') then
            call haze_init(tau_hazev_init, g_hazev, ssa_hazev, haze_presc_file, plev, wnov, 'v')
            call haze_init(tau_hazei_init, g_hazei, ssa_hazei, haze_presc_file,  plev, wnoi, 'i')
            haze_or_clouds = .true.
        end if

        if (use_haze_clim) then
            print *, 'Radiation:  Haze climatology loaded'
        end if

        module_is_initialized = .true.
END SUBROUTINE planetary_radiation_init

! ========================================================================================

SUBROUTINE planetary_radiation(is, js, lat, lon, testing_ls, cosz, phalf, pfull, zhalf, zfull, &
                               albedov, albedoi, t_surf, t, tracers, tdt, tdt_sw, tdt_lw, sfc_flux, lw_spec, sw_spec)

        integer, intent(in)                             :: is, js
        !type(time_type), intent(in)			:: Time
	real, intent(in), dimension(:,:)                :: lat, lon, t_surf, albedov, albedoi, cosz
        !real, intent(in), dimension(:,:)               :: t_surf, albedov, albedoi    
        real, intent(in), dimension(:,:,:)              :: t, phalf, pfull, zhalf, zfull
	real, intent(in), dimension(:,:,:,:)            :: tracers
	real, intent(out), dimension(:,:)               :: sfc_flux
	real, intent(out), dimension(:,:,:)             :: tdt

	real, dimension(:, :), allocatable              :: tauz

        integer                                         :: i, j, k, g, n, X, ii, dd, l
	real                                            :: testing_ls
        integer, dimension(L_LEVELS)                    :: ref_index, clidx
	logical                                         :: used, use_column_ch4
	real, dimension(size(lat,1),size(lat,2))        :: fracday, toa_netv, toa_neti !cosz
	real                                            :: rrsun, orbital_time, rsdist, albv, albi
	real, dimension(size(wnov,1))                   :: sol
	real, dimension(L_LEVELS)                       :: plev, tlev, zlev, q, rho_moles, dist
	real, dimension(L_LEVELS)                       :: X_sphum, mixing_ratio, col_abund, CH4_vmr
	real, dimension(L_LEVELS)                       :: num_density1, num_density2
	real, dimension(L_LEVELS)                       :: ref_weight, clwt
	real, dimension(L_LEVELS,size(trace_gas))       :: rad_tracers
	real, dimension(L_LEVELS,size(wnov,1))          :: tau_rayleigh, tau_ciav
	real, dimension(L_LEVELS,size(wnov,1))          :: tau_hazev, tau_cloudv
	real, dimension(L_LEVELS,size(wnov,1))          :: ssa_cloudv
	real, dimension(L_LEVELS,size(wnov,1))          ::  g_cloudv
	real, dimension(L_LEVELS,size(wnov,1),kgauss)   :: tau_gasv
	real, dimension(L_LEVELS,size(wnoi,1))          :: tau_ciai, tau_hazei, tau_cloudi, tau_rayleighi
	real, dimension(L_LEVELS,size(wnoi,1),size(ciapairi,1)) :: tau_ciai_pair
	real, dimension(L_LEVELS,size(wnoi,1))          :: ssa_cloudi
	real, dimension(L_LEVELS,size(wnoi,1))          ::  g_cloudi
	real, dimension(L_LEVELS,size(wnoi,1),kgauss)   :: tau_gasi
        character(len=10)                               :: gas1, gas2
        character(len=64)                               :: cia_pair_file

	real, dimension(L_LAYERS+1,size(wnov,1),kgauss) :: dtauv, cosbv, wbarv
	real, dimension(L_LAYERS+2,size(wnov,1),kgauss) :: tauv
	real, dimension(L_LEVELS,size(wnov,1),kgauss)   :: taucumv
	real, dimension(L_LAYERS+1,size(wnoi,1),kgauss) :: dtaui, cosbi, wbari
	real, dimension(L_LAYERS+2,size(wnoi,1),kgauss) :: taui
	real, dimension(L_LEVELS,size(wnoi,1),kgauss)   :: taucumi
	real                                            :: nfluxtopv, nfluxtopi, diffv
	real, dimension(L_LAYERS+1)                     :: fmnetv, fluxupv, fluxdnv
	real, dimension(L_LAYERS+1)                     :: fmneti, fluxupi, fluxdni
	real, dimension(L_LAYERS)                       :: Cp_layer,heatingv,heatingi
	real, dimension(size(wnov,1))                   :: sw_spectrum
	real, dimension(size(wnoi,1))                   :: lw_spectrum

	real, dimension(size(lat,1),size(lat,2),L_LAYERS+1)                     :: netv,dnv,upv,neti,dni,upi
	real, dimension(size(lat,1),size(lat,2),L_LAYERS), intent(out)          :: tdt_sw,tdt_lw
	real, dimension(size(lat,1),size(lat,2),size(wnov,1)), intent(out)      :: sw_spec
	real, dimension(size(lat,1),size(lat,2),size(wnoi,1)), intent(out)      :: lw_spec

        integer                                         :: cls_ind
	real                                            :: cls_wt

!-----------------------------------------------------------------------------------------
! Verify that module is initialized
!-----------------------------------------------------------------------------------------

        if (.not. module_is_initialized) call error_mesg('planetary_radiation_driver_mod', &
                                                         'module has not been initialized', FATAL)

!-----------------------------------------------------------------------------------------
! Call astronomy to get cosz, fracday, rrsun, and orbital time; 
! calculate stellar flux at current planet-star distance
!-----------------------------------------------------------------------------------------

        !call diurnal_solar(lat,lon,Time,cosz,fracday,rrsun,orbital_time)
        rrsun = 1. !SA sun / 4 pi, meters
        rsdist = rrsun/(sma**2)
        sol(:) = solarf(:) * rsdist  !* 1.D13

!-----------------------------------------------------------------------------------------
! Begin loops over latitude (j) and longitude (i)
!-----------------------------------------------------------------------------------------
        !Interpolate TGC and HC onto current timestep ls
        if (use_clim_gases) then
                if (testing_ls .le. clls(1)) then
                        !interpolate over cycle
                        cls_ind = 1
                        cls_wt =  (testing_ls-clls(size(clls))-360.) / (clls(g) - clls(size(clls))-360.)
                        cl_inst(:, :, :) = clim_data(:, :, size(clls), :) + (clim_data(:, :, cls_ind, :) - clim_data(:, :, size(clls), :))*cls_wt
                elseif (testing_ls .ge. clls(size(clls))) then
                        !Do same thing as above
                        cls_ind = 1
                        cls_wt =  (testing_ls-clls(size(clls))-360.) / (clls(g) - clls(size(clls))-360.)
                        cl_inst(:, :, :) = clim_data(:, :, size(clls), :) + (clim_data(:, :, cls_ind, :) - clim_data(:, :, size(clls), :))*cls_wt
                else
                        do g = 1, size(clls)
                                if (clls(g) .ge. testing_ls) then
                                        cls_ind = g-1
                                        !Special case for first index
                                        if (cls_ind - 1 .ne. 0) then
                                                cls_wt =  (testing_ls-clls(g-1)) / (clls(g) - clls(g-1))
                                                cl_inst(:, :, :) = clim_data(:, :, cls_ind-1, :) + (clim_data(:, :, cls_ind, :) - clim_data(:, :, cls_ind-1, :))*cls_wt
                                                exit
                                        else
                                                cls_wt =  (testing_ls-clls(size(clls))) / (clls(g) - clls(size(clls)))
                                                cl_inst(:, :, :) = clim_data(:, :, size(clls), :) + (clim_data(:, :, cls_ind, :) - clim_data(:, :, size(clls), :))*cls_wt
                                                exit
                                        end if
                                end if
                        end do
                end if
        end if

        do j = 1,size(t,1)
                do i = 1,size(t,1)

!-----------------------------------------------------------------------------------------
! Set up the pressure (Pa) and temperature (K) grid needed by the radiation code;
! note that with a spectral dynamical core, p(1) might be 0 
!-----------------------------------------------------------------------------------------

                        do k = 1,L_LAYERS
                                plev(2*k+1) = phalf(i,j,k)
                                plev(2*k+2) = pfull(i,j,k)
                                tlev(2*k+2) = t(i,j,k)
                        end do
                       
                        if (plev(3) .eq. 0.0) plev(3) = plev(4)/10.
                        !plev extrapolation may be causing problem later, not actually following grid
                        plev(2) = plev(3)*plev(3)/plev(4)
                        plev(1) = plev(2)*plev(2)/plev(3) 
                        plev(L_LEVELS) = phalf(i, j, L_LAYERS+1)

                        do k = 5,L_LEVELS-2,2
                                tlev(k) = tlev(k+1) + (tlev(k-1) - tlev(k+1)) * log(plev(k)/plev(k+1)) / log(plev(k-1)/plev(k+1))
                        end do

                        tlev(L_LEVELS) = t_surf(i,j)

                        call apply_top_tlev_boundary(plev, tlev)

!-----------------------------------------------------------------------------------------
! Set up other layer variables (rho_moles in moles/m3; dist in km)
!-----------------------------------------------------------------------------------------

                        use_column_ch4 = use_analytic_ch4_profile .or. len(trim(sphum_gas)) .gt. 0

                        if (use_column_ch4) then
                                if (use_analytic_ch4_profile) then
                                        do k = 1,L_LAYERS
                                                q(2*k+1) = tracers(i,j,k,1)
                                        end do
                                else
                                        do k = 1,L_LAYERS
                                                q(2*k+1) = tracers(i,j,k,sphum_index)
                                        end do
                                end if
                                q(1) = q(3)
                                q(2) = q(3)
                                if (use_analytic_ch4_profile) then
                                        q(L_LEVELS) = tracers(i,j,L_LAYERS,1)
                                else
                                        q(L_LEVELS) = tracers(i,j,L_LAYERS,sphum_index)
                                end if
                                do k = 4,L_LEVELS-1,2
                                        q(k) = q(k+1) + (q(k-1) - q(k+1)) * log(plev(k)/plev(k+1)) / log(plev(k-1)/plev(k+1))
                                end do
                                X_sphum = q / ( (Rdgas/Rvgas) - ( (Rdgas/Rvgas) - 1.)*q )
                        end if

                        do g = 1,size(trace_gas,1)
                                if (tracer_indices(g) .ne. NO_TRACER) then
                                        do k = 1,L_LAYERS
                                                rad_tracers(2*k+1,g) = tracers(i,j,k,tracer_indices(g))
                                        end do
                                        rad_tracers(1,g) = rad_tracers(3,g)
                                        rad_tracers(2,g) = rad_tracers(3,g)
                                        rad_tracers(L_LEVELS,g) = tracers(i,j,L_LAYERS,tracer_indices(g))
                                        do k = 4,L_LEVELS-1,2
                                                rad_tracers(k,g) = rad_tracers(k+1,g) + (rad_tracers(k-1,g) - &
                                                                   rad_tracers(k+1,g)) * log(plev(k)/plev(k+1))/log(plev(k-1)/plev(k+1))
                                        end do
                                end if
                        end do

                        do k = 1,L_LEVELS
                                rho_moles(k)= plev(k)/(tlev(k)*Runiv) !THIS IS MOLAR VOLUME
                        end do   
                        
                        do k = 1,L_LAYERS
                                zlev(2*k+1) = zhalf(i,j,k)
                                zlev(2*k+2) = zfull(i,j,k)
                        end do
                        
                        zlev(L_LEVELS) = zhalf(i,j,L_LAYERS+1)
                        !zlev(2) = zlev(3)+(Rdgas * (plev(3)-plev(2))*tlev(2))/(grav)
                        !zlev(1) = zlev(2)+(Rdgas * (plev(2)-plev(1))*tlev(1))/(grav)
                        zlev(2) = zlev(3) - (log(plev(2))-log(plev(3)))*(Rdgas/grav)*(tlev(3))
                        zlev(1) = zlev(2) - (log(plev(1))-log(plev(2)))*(Rdgas/grav)*(tlev(2))
                        !zlev(3) + (rdgas*(log(plev(3))-log(plev(2))) * tlev(3)*(1.+(rvgas/rdgas-1.)*q(2)))/(grav*1000.)
                        
                        do k = 2,L_LEVELS
                                dist(k) = (zlev(k-1) - zlev(k))/1000.
                        end do
                        dist(1) = 0
                        
                        open (unit = 97, file = 'dist.txt')
                        do k=1, L_LEVELS
                                write (97, ' (  500(X, E12.5) ) ') plev(k), dist(k), zlev(k)
                        end do
                        close(97)
!-----------------------------------------------------------------------------------------
! Set surface albedos
!-----------------------------------------------------------------------------------------

                        albv = albedov(i,j)
                        albi = albedoi(i,j)

!-----------------------------------------------------------------------------------------
! Set up interpolation to prescribed profiles' pressure grid
!-----------------------------------------------------------------------------------------
                        if (use_presc_gases) then
                                do k = 1,L_LEVELS
                                        if (plev(k) .le. ref_gases(1,1)) then
                                                ref_index(k)  = 1
                                                ref_weight(k) = 0.
                                                cycle
                                        elseif (plev(k) .ge. ref_gases(1,ppress)) then
                                                ref_index(k)  = ppress
                                                ref_weight(k) = 0.
                                                cycle
                                        else
                                                do g = 1,size(ref_gases,2)
                                                        if (plev(k) .le. ref_gases(1, g)) then
                                                                ref_index(k)  = g
                                                                ref_weight(k) = log(ref_gases(1,g)/plev(k)) / &
                                                                                log(ref_gases(1,g)/ref_gases(1,g-1))
                                                                exit
                                                        end if
                                                end do
                                        end if
                                end do
                        end if

!-----------------------------------------------------------------------------------------
! Interpolate to climatology pressure
!-----------------------------------------------------------------------------------------

                        if (use_clim_gases) then
                                do k = 1,L_LEVELS
                                        if (plev(k) .le. clpress(1)) then
                                                clidx(k)  = 1
                                                clwt(k) = 0.
                                                cycle
                                        elseif (plev(k) .ge. clpress(SIZE(clpress))) then
                                                clidx(k)  = SIZE(clpress)
                                                clwt(k) = 0.
                                                cycle
                                        else
                                                do g = 1,size(clpress)
                                                        if (plev(k) .le. clpress(g)) then
                                                                clidx(k)  = g
                                                                clwt(k) = log(clpress(g)/plev(k))/log(clpress(g)/clpress(g-1))
                                                                exit
                                                        end if
                                                end do
                                        end if
                                end do
                        end if

                        !open(11, file = 'clwt.prof')
                        !do k = 1, L_LEVELS
                        !    write (11, ' (  500(X, E10.3) ) ') plev(k), clwt(k)
                        !end do
                        !close(11)
    
!-----------------------------------------------------------------------------------------
! If the sun is up, set up different opacities and begin calculation of solar fluxes
!-----------------------------------------------------------------------------------------

                        if (cosz(i,j) .ge. 1.0D-4) then
                                
                                tau_rayleigh = 0.D0
                                tau_gasv     = 0.D0
                                tau_ciav     = 0.D0
                                tau_hazev    = 0.D0
                                tau_cloudv   = 0.D0

!-----------------------------------------------------------------------------------------
! Get visible optical depths for Rayleigh scattering
!-----------------------------------------------------------------------------------------

                                do k = 2,L_LEVELS
                                        tau_rayleigh(k,:) = (plev(k) - plev(k-1))*tauray(:)
                                end do
                                ! Not defined at k = 0, set as next
                                tau_rayleigh(1,:) = tau_rayleigh(2,:)

!-----------------------------------------------------------------------------------------
! Get visible optical depths from k-coefficients
!-----------------------------------------------------------------------------------------

                                if (do_radgasv) then
                                        !print *, 'Do radgasv'
                                        col_abund = 0.0
                                        do g = 1,size(radgasv,1)
                                                if (len(trim(radgasv(g))) .gt. 0) then
                                                        if ((len(trim(sphum_gas)) .gt. 0 .and. radgasv(g) .eq. sphum_gas) .or. &
                                                            (use_analytic_ch4_profile .and. radgasv(g) .eq. 'CH4')) then
                                                                mixing_ratio = X_sphum
                                                                col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                                call get_taukcoeff( tau_gasv, col_abund, tlev, plev, &
                                                                                    ckcv(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                cycle
                                                        end if
                                                        if (radgasv(g) .eq. bkgnd_gas) then
                                                                if (use_column_ch4) then
                                                                        mixing_ratio = 1. - X_sphum
                                                                else
                                                                        mixing_ratio = 1.
                                                                end if
                                                                col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                                call get_taukcoeff( tau_gasv, col_abund, tlev, plev, &
                                                                                    ckcv(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                cycle
                                                        end if
                                                        do n = 1,size(presc_gas)
                                                                if (radgasv(g) .eq. presc_gas(n)) then
                                                                        do k = 1,L_LEVELS
                                                                                if (ref_weight(k) .eq. 0.) then
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k)) 
                                                                                else
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k))- &
                                                                                                          ref_weight(k)*(ref_gases(n+1,ref_index(k))- &
                                                                                                          ref_gases(n+1,ref_index(k)-1))
                                                                                end if
                                                                        end do
                                                                        col_abund = mixing_ratio*rho_moles*mole_to_amg*dist
                                                                        call get_taukcoeff( tau_gasv, col_abund, tlev, plev, &
                                                                                            ckcv(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                        exit
                                                                end if
                                                        end do
                                                        do n = 1,size(trace_gas)
                                                                if (radgasv(g) .eq. trace_gas(n)) then
                                                                        mixing_ratio = rad_tracers(:,n)
                                                                        col_abund = mixing_ratio*rho_moles*mole_to_amg*dist
                                                                        call get_taukcoeff( tau_gasv, col_abund, tlev, plev, &
                                                                                            ckcv(g,:,:,:,:),ckc_temp,ckc_pres)
                                                                        exit
                                                                end if
                                                        end do
                                                        ! Get latitude/time index for clim_data
                                                        ! indices are pres, time, lat, gas
                                                        do n = 1,size(clim_gas)
                                                                if (radgasv(g) .eq. clim_gas(n)) then
                                                                        print *, clim_gas(n)
                                                                        print *, cl_inst(n, :, :)
                                                                        do k = 1,L_LEVELS
                                                                                if (clwt(k) .eq. 0.) then
                                                                                        mixing_ratio(k) = cl_inst(n, j, clidx(k))
                                                                                else
                                                                                        mixing_ratio(k) = cl_inst(n, j, clidx(k))/&
                                                                                                          (cl_inst(n, j, clidx(k))/&
                                                                                                          cl_inst(n, j, clidx(k)-1))**clwt(k)
                                                                                end if
                                                                        end do
                                                                        col_abund = mixing_ratio*rho_moles*mole_to_amg*dist
                                                                        call get_taukcoeff( tau_gasv, col_abund, tlev, plev, &
                                                                                            ckcv(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                        exit
                                                                end if
                                                        end do
                                                end if
                                        end do

                                        ! save Gauss-weighted (mean) gas optical depth per band per level,
                                        ! same layout as tau_hazev.txt (row0 = plev, then one row per band)
                                        open (unit = 90, file = 'tau_gasv.txt', status = 'replace')
                                            write(90, ' (  500(X, E15.8) ) ') plev
                                            do dd = 1, size(tau_gasv, 2)
                                                 write(90, ' (  500(X, E15.8) ) ') matmul(tau_gasv(:, dd, :), ckc_gwtv)
                                            end do
                                        close(90)

                                end if

!-----------------------------------------------------------------------------------------
! Get visible optical depths from CIA pairs
!-----------------------------------------------------------------------------------------
                                ! CIA NOT INCLUDING CLIMATOLOGY YET
                                if (do_ciav) then
                                        do g = 1,size(ciapairv,1)
                                                num_density1 = 0.0
                                                num_density2 = 0.0
                                                if (len(trim(ciapairv(g))) .gt. 0) then
                                                        gas1 = ciapairv(g)(1:index(ciapairv(g),'-')-1)
                                                        gas2 = trim(ciapairv(g)(index(ciapairv(g),'-')+1:10))
                                                        if ((len(trim(sphum_gas)) .gt. 0 .and. gas1 .eq. sphum_gas) .or. &
                                                            (use_analytic_ch4_profile .and. gas1 .eq. 'CH4')) then
                                                                mixing_ratio = X_sphum
                                                                num_density1 = mixing_ratio * rho_moles * avo 
                                                                !NUMBER DENSITY IN MOLEC/M^3
                                                        elseif (gas1 .eq. bkgnd_gas) then
                                                                if (use_column_ch4) then
                                                                        mixing_ratio = 1. - X_sphum
                                                                else
                                                                        mixing_ratio = 1.
                                                                end if
                                                                num_density1 = mixing_ratio * rho_moles * avo
                                                        else
                                                                do n = 1,size(presc_gas)
                                                                        if (gas1 .eq. presc_gas(n)) then
                                                                                do k = 1,L_LEVELS
                                                                                        if (ref_weight(k) .eq. 0.) then
                                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k))
                                                                                        else
                                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k)) - &
                                                                                                                  ref_weight(k)*(ref_gases(n+1,ref_index(k))-&
                                                                                                                  ref_gases(n+1,ref_index(k)-1))
                                                                                        end if
                                                                                end do
                                                                                num_density1 = mixing_ratio * rho_moles * avo
                                                                                exit
                                                                        end if
                                                                end do
                                                                do n = 1,size(trace_gas)
                                                                        if (gas1 .eq. trace_gas(n)) then
                                                                                mixing_ratio = rad_tracers(:,n)
                                                                                num_density1 = mixing_ratio * rho_moles * avo
                                                                        end if
                                                                        exit
                                                                end do
                                                        end if
                                                        if ((len(trim(sphum_gas)) .gt. 0 .and. gas2 .eq. sphum_gas) .or. &
                                                            (use_analytic_ch4_profile .and. gas2 .eq. 'CH4')) then
                                                                mixing_ratio = X_sphum
                                                                num_density2 = mixing_ratio * rho_moles * avo
                                                        elseif (gas2 .eq. bkgnd_gas) then
                                                                if (use_column_ch4) then
                                                                        mixing_ratio = 1. - X_sphum
                                                                else
                                                                        mixing_ratio = 1.
                                                                end if
                                                                num_density2 = mixing_ratio * rho_moles * avo
                                                        else
                                                                do n = 1,size(presc_gas)
                                                                        if (gas2 .eq. presc_gas(n)) then
                                                                                do k = 1,L_LEVELS
                                                                                        if (ref_weight(k) .eq. 0.) then
                                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k))
                                                                                        else
                                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k)) - &
                                                                                                                  ref_weight(k)*(ref_gases(n+1,ref_index(k))-&
                                                                                                                  ref_gases(n+1,ref_index(k)-1))
                                                                                        end if
                                                                                end do
                                                                                num_density2 = mixing_ratio * rho_moles * avo
                                                                                exit
                                                                        end if
                                                                end do
                                                                do n = 1,size(trace_gas)
                                                                        if (gas2 .eq. trace_gas(n)) then
                                                                                mixing_ratio = rad_tracers(:,n)
                                                                                num_density2 = mixing_ratio * rho_moles * avo
                                                                        end if
                                                                        exit
                                                                end do
                                                        end if
                                                        call get_tauCIA(tau_ciav,num_density1*num_density2*dist*1.E-7, &
                                                                        tlev,ciav(g,:,:,:),tfits,cia_temp )
                                                        ! The above line uses number densities in molec/m^3, multiplied 
                                                        ! by distance in km, times 1e-7, to pass a column abundance
                                                        ! in molec^2/cm^5
                                                end if
                                        end do
                                end if 

!-----------------------------------------------------------------------------------------
! Get visible optical depths from hazes and clouds
!-----------------------------------------------------------------------------------------		 				 
                
                                if (TRIM(haze_data) .eq. 'clim') then
                                        call clim_haze_profile(tau_hazev, j, 1, wnov, plev)
                                else
                                        tau_hazev = tau_hazev_init
                                end if
                
                                do dd=2, size(tau_hazev,1)                                                                                           
                                        dtau_hazev(dd, :) = tau_hazev(dd,:) - tau_hazev(dd-1,:)
                                end do
                
                                dtau_hazev(1,:) = 2*dtau_hazev(2,:) - dtau_hazev(3,:)

!-----------------------------------------------------------------------------------------
! Calculate total optical constants, and compute visible fluxes
!-----------------------------------------------------------------------------------------		
                                open (unit = 90, file = 'tau_hazev.txt', status = 'replace')
                                        write(90, ' (  500(X, E15.8) ) ') plev
                                        do dd = 1, size(tau_hazev, 2)
                                                 write(90, ' (  500(X, E15.8) ) ') tau_hazev(:, dd)
                                        end do
                                close(90)

                                open (unit = 90, file = 'dtau_hazev.txt', status = 'replace')
                                        write(90, ' (  500(X, E15.8) ) ') plev
                                        do dd = 1, size(dtau_hazev, 2)
                                                 write(90, ' (  500(X, E15.8) ) ') dtau_hazev(:, dd)
                                        end do
                                close(90)
                
                                open (unit = 90, file = 'ssa_hazev.txt', status = 'replace')
                                        write(90, ' (  500(X, E15.8) ) ') plev
                                        do dd = 1, 103
                                                 write(90, ' (  500(X, E15.8) ) ') ssa_hazev(:, dd)
                                        end do
                                close(90)
       
                                if (haze_or_clouds) then
                                        call optc( dtauv,tauv,taucumv,wbarv,cosbv, &
                                                   tau_gasv,tau_ray_in=tau_rayleigh,tau_cia_in=tau_ciav,tau_haze_in=dtau_hazev,&
                                                   tau_cloud_in=tau_cloudv,ssa_haze_in=ssa_hazev,&
                                                   ssa_cloud_in=ssa_cloudv,g_haze_in=g_hazev,g_cloud_in=g_cloudv )
                                else
                                        call optc( dtauv,tauv,taucumv,wbarv,cosbv, &
                                                   tau_gasv,tau_ray_in = tau_rayleigh, tau_cia_in = tau_ciav)
                                end if
                                ! wbarv = 1 causes numerical crashing
                                WHERE (wbarv .ge. 0.999) wbarv = 0.999
                                call sfluxv( nfluxtopv,fmnetv,fluxupv,fluxdnv,diffv,sw_spectrum, &
                                             dtauv,tauv,taucumv,ckc_gwtv,albv,wbarv,cosbv,cosz(i,j),sol )
                        else
                                nfluxtopv = 0.0
                                fmnetv    = 0.0
                                fluxupv    = 0.0
                                fluxdnv   = 0.0
                        end if

!-----------------------------------------------------------------------------------------
! Begin calculation of IR fluxes
!-----------------------------------------------------------------------------------------				

                        tau_gasi     = 0.D0
                        tau_ciai     = 0.D0
                        tau_ciai_pair = 0.D0
                        tau_hazei    = 0.D0
                        tau_cloudi   = 0.D0

!-----------------------------------------------------------------------------------------
! Get IR optical depths from k-coefficients
!-----------------------------------------------------------------------------------------

                        if (do_radgasi) then
                                do g = 1,size(radgasi,1)
                                        col_abund = 0.0
                                        if (len(trim(radgasi(g))) .gt. 0) then
                                                if ((len(trim(sphum_gas)) .gt. 0 .and. radgasi(g) .eq. sphum_gas) .or. &
                                                    (use_analytic_ch4_profile .and. radgasi(g) .eq. 'CH4')) then
                                                        mixing_ratio = X_sphum
                                                        col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                        call get_taukcoeff( tau_gasi, col_abund, tlev, plev, &
                                                                            ckci(g,:,:,:,:),ckc_temp,ckc_pres )
                                                        cycle
                                                end if
                                                if (radgasi(g) .eq. bkgnd_gas) then
                                                        if (use_column_ch4) then
                                                                mixing_ratio = 1. - X_sphum
                                                        else
                                                                mixing_ratio = 1.
                                                        end if
                                                        col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                        call get_taukcoeff( tau_gasi, col_abund, tlev, plev, &
                                                                            ckci(g,:,:,:,:),ckc_temp,ckc_pres )
                                                        cycle
                                                end if
                                                do n = 1,size(presc_gas)
                                                        if (radgasi(g) .eq. presc_gas(n)) then
                                                                do k = 1,L_LEVELS
                                                                        if (ref_weight(k) .eq. 0.) then
                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k))
                                                                        else
                                                                                mixing_ratio(k) = ref_gases(n+1,ref_index(k)) - &
                                                                                                  ref_weight(k)*(ref_gases(n+1,ref_index(k)) - &
                                                                                                  ref_gases(n+1,ref_index(k)-1))
                                                                        end if
                                                                end do
                                                                col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                                
                                                                ! Temporary fix for N2 abundance
                                                                if (presc_gas(n) .eq. 'CH4') then
                                                                        CH4_vmr = mixing_ratio
                                                                end if
                                                                
                                                                call get_taukcoeff( tau_gasi, col_abund, tlev, plev, &
                                                                                    ckci(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                exit
                                                        end if
                                                end do

                                                ! indices are pres, time, lat, gas
                                                ! or define another clmolidx array that maps clim_gas molecule to clim index
                                                do n = 1,size(clim_gas)
                                                        if (radgasi(g) .eq. clim_gas(n)) then
                                                                do k = 1,L_LEVELS
                                                                        if (clwt(k) .eq. 0.) then
                                                                                mixing_ratio(k) = cl_inst(clmolidx(n), clidx(k), j)
                                                                        else
                                                                                mixing_ratio(k) = (cl_inst(clmolidx(n), clidx(k), j))/&
                                                                                                  ((cl_inst(clmolidx(n), clidx(k),j) /& 
                                                                                                  cl_inst(clmolidx(n), clidx(k)-1, j))**clwt(k))
                                                                        end if
                                                                end do
                                                                
                                                                ! Saturate molecules if they are oversaturated, RH defined in subroutine
                                                                call saturate(mixing_ratio, radgasi(g), tlev, plev)

                                                                col_abund = mixing_ratio*rho_moles*mole_to_amg*dist
                                                                call get_taukcoeff( tau_gasi, col_abund, tlev, plev, &
                                                                                    ckci(g,:,:,:,:),ckc_temp,ckc_pres )

                                                                ! Save molecule VMRs
                                                                if (.true.) then
                                                                        open (unit = 97, file = trim(clim_gas(n))//'vmr.txt')
                                                                        do k=1, L_LEVELS
                                                                                write (97, ' (  500(X, E12.5) ) ') plev(k), mixing_ratio(k), real(clidx(k))
                                                                        end do
                                                                        close(97)
                                                                end if
                                                                exit
                                                        end if
                                                end do
                                                do n = 1,size(trace_gas)
                                                        if (radgasi(g) .eq. trace_gas(n)) then
                                                                mixing_ratio = rad_tracers(:,n)
                                                                col_abund = mixing_ratio * rho_moles * mole_to_amg * dist
                                                                call get_taukcoeff( tau_gasi, col_abund, tlev, plev, &
                                                                                    ckci(g,:,:,:,:),ckc_temp,ckc_pres )
                                                                exit
                                                        end if
                                                end do
                                        end if
                                end do
                                ! save Gauss-weighted (mean) IR gas-line optical depth per band per level
                                open (unit = 90, file = 'tau_gasi.txt', status = 'replace')
                                    write(90, ' (  500(X, E15.8) ) ') plev
                                    do dd = 1, size(tau_gasi, 2)
                                             write(90, ' (  500(X, E15.8) ) ') matmul(tau_gasi(:, dd, :), ckc_gwti)
                                    end do
                                close(90)
                        end if
!-----------------------------------------------------------------------------------------
! Get IR optical depths from CIA pairs
!-----------------------------------------------------------------------------------------

                        if (do_ciai) then
                                do g = 1,size(ciapairi,1)
                                        num_density1 = 0.0
                                        num_density2 = 0.0
                                        if (len(trim(ciapairi(g))) .gt. 0) then
                                                !print *, 'Doing CIA pair for ', ciapairi(g)
                                                gas1 = trim(ciapairi(g)(1:index(ciapairi(g),'-')-1))
                                                gas2 = trim(ciapairi(g)(index(ciapairi(g),'-')+1:10))
                                                if ((len(trim(sphum_gas)) .gt. 0 .and. gas1 .eq. sphum_gas) .or. &
                                                    (use_analytic_ch4_profile .and. gas1 .eq. 'CH4')) then
                                                        mixing_ratio = X_sphum
                                                        num_density1 = mixing_ratio * rho_moles * avo
                                                elseif (gas1 .eq. bkgnd_gas) then
                                                        if (use_column_ch4) then
                                                                mixing_ratio = 1. - X_sphum
                                                        else
                                                                mixing_ratio = 1.
                                                        end if
                                                        num_density1 = mixing_ratio * rho_moles * avo
                                                else
                                                        do n = 1,size(presc_gas)
                                                                if (gas1 .eq. presc_gas(n)) then
                                                                        do k = 1,L_LEVELS
                                                                                if (ref_weight(k) .eq. 0.) then
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k))
                                                                                else
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k)) -&
                                                                                                          ref_weight(k)*(ref_gases(n+1,ref_index(k))- &
                                                                                                          ref_gases(n+1,ref_index(k)-1))
                                                                                end if
                                                                        end do
                                                                        num_density1 = mixing_ratio * rho_moles * avo
                                                                        !open (unit = 98, file = gas1//'vmr.txt')
                                                                                !write (98, ' (  500(X, E10.3) ) ') mixing_ratio
                                                                        !close(98)
                                                                        exit
                                                                end if
                                                        end do
                                                        do n = 1,size(trace_gas)
                                                                if (gas1 .eq. trace_gas(n)) then
                                                                        mixing_ratio = rad_tracers(:,n)
                                                                        num_density1 = mixing_ratio * rho_moles * avo
                                                                end if
                                                                exit
                                                        end do
                                                end if
                                                if ((len(trim(sphum_gas)) .gt. 0 .and. gas2 .eq. sphum_gas) .or. &
                                                    (use_analytic_ch4_profile .and. gas2 .eq. 'CH4')) then
                                                        mixing_ratio = X_sphum
                                                        num_density2 = mixing_ratio * rho_moles * avo
                                                elseif (gas2 .eq. bkgnd_gas) then
                                                        if (use_column_ch4) then
                                                                mixing_ratio = 1. - X_sphum
                                                        else
                                                                mixing_ratio = 1.
                                                        end if
                                                        num_density2 = mixing_ratio * rho_moles * avo
                                                else
                                                        do n = 1,size(presc_gas)
                                                                if (gas2 .eq. presc_gas(n)) then
                                                                        do k = 1,L_LEVELS
                                                                                if (ref_weight(k) .eq. 0.) then
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k))
                                                                                else
                                                                                        mixing_ratio(k) = ref_gases(n+1,ref_index(k)) -&
                                                                                                          ref_weight(k)*(ref_gases(n+1,ref_index(k))-&
                                                                                                          ref_gases(n+1,ref_index(k)-1))
                                                                                end if
                                                                        end do
                                                                        !open (unit = 97, file = trim(gas2)//'vmr.txt')
                                                                        !        write (97, ' (  500(X, E10.3) ) ') mixing_ratio
                                                                        !close(97)
                                                                        num_density2 = mixing_ratio * rho_moles * avo
                                                                        exit
                                                                end if
                                                        end do
                                                        do n = 1,size(trace_gas)
                                                                if (gas2 .eq. trace_gas(n)) then
                                                                        mixing_ratio = rad_tracers(:,n)
                                                                        num_density2 = mixing_ratio * rho_moles * avo
                                                                end if
                                                                exit
                                                        end do
                                                end if

                                                !write(90, ' (  500(X, E10.3) ) ') ((num_density1*num_density2)) *dist*1.E-7
                                                call get_tauCIA( tau_ciai_pair(:,:,g),((num_density1*num_density2)) * dist*1.E-7, &
                                                                 tlev,ciai(g,:,:,:),tfits,cia_temp )
                                                tau_ciai = tau_ciai + tau_ciai_pair(:,:,g)
                                        end if
                                end do
                        end if

                        open (unit = 90, file = 'tau_ciai.txt', status = 'replace')
                                write(90, ' (  500(X, E15.8) ) ') wnoi
                                write(90, ' (  500(X, E15.8) ) ') plev
                                do dd = 1, size(tau_ciai,2)
                                        write(90, ' (  500(X, E15.8) ) ') tau_ciai(:, dd)
                                end do
                                write(90, *) 'Closed'
                        close(90)
                        if (do_ciai .and. write_cia_pair_diagnostics) then
                                do g = 1,size(ciapairi,1)
                                        if (len(trim(ciapairi(g))) .gt. 0) then
                                                write(cia_pair_file,'(A,A,A)') 'tau_ciai_', trim(ciapairi(g)), '.txt'
                                                open (unit = 90, file = trim(cia_pair_file), status = 'replace')
                                                        write(90, ' (  500(X, E15.8) ) ') wnoi
                                                        write(90, ' (  500(X, E15.8) ) ') plev
                                                        do dd = 1, size(tau_ciai_pair,2)
                                                                write(90, ' (  500(X, E15.8) ) ') tau_ciai_pair(:, dd, g)
                                                        end do
                                                        write(90, *) 'Closed'
                                                close(90)
                                        end if
                                end do
                        end if
            
!-----------------------------------------------------------------------------------------
! Get IR optical depths from hazes and clouds
!-----------------------------------------------------------------------------------------		 				 

                        if (TRIM(haze_data) .eq. 'clim') then
                                call clim_haze_profile(tau_hazei, j, 0, wnoi, plev)
                        else
                                tau_hazei = tau_hazei_init
                        end if

                        do dd=2, size(tau_hazei, 1)
                                dtau_hazei(dd,:) = tau_hazei(dd, :) - tau_hazei(dd-1,:)
                        end do
                        dtau_hazei(1,:) = 2*dtau_hazei(2,:) - dtau_hazei(3,:)

                        open (unit = 90, file = 'tau_hazei.txt', status = 'replace')
                                write(90, ' (  500(X, E15.8) ) ') plev 
                                write(90, ' (  500(X, E15.8) ) ') wnoi
                                do dd = 1, size(dtau_hazei,2)
                                        write(90, ' (  500(X, E15.8) ) ') tau_hazei(:, dd)
                                end do
                        close(90)

!-----------------------------------------------------------------------------------------
! Calculate total optical constants, and compute IR fluxes
!-----------------------------------------------------------------------------------------

                        if (haze_or_clouds) then
                                call optc( dtaui,taui,taucumi,wbari,cosbi, &
                                           tau_gasi,tau_cia_in=tau_ciai,tau_haze_in=dtau_hazei,tau_cloud_in=tau_cloudi, &
                                           ssa_haze_in=ssa_hazei,ssa_cloud_in=ssa_cloudi,g_haze_in=g_hazei,g_cloud_in=g_cloudi )
                        else
                                call optc( dtaui,taui,taucumi,wbari,cosbi, tau_gasi, tau_cia_in = tau_ciai )
                        end if
                        
                        call sfluxi( nfluxtopi,fmneti,fluxupi,fluxdni,lw_spectrum, &
                                     plev,tlev,dtaui,taucumi,ckc_gwti,dwni,albi,   &
                                     wbari,cosbi,ubari,planckir )

!-----------------------------------------------------------------------------------------
! Calculate heating rates
!-----------------------------------------------------------------------------------------
                        do k = 1,L_LAYERS
                                if (use_column_ch4) then
                                        Cp_layer(k) = Cp_air * (1. - q(2*k+1)) + Cp_vapor * q(2*k+1)
                                else
                                        Cp_layer(k) = Cp_air
                                end if

                                heatingv(k) = ( fmnetv(k) - fmnetv(k+1) ) * (grav) / ( Cp_layer(k) * (plev(2*k+3) - plev(2*k+1)) )
                                heatingi(k) = ( fmneti(k+1) - fmneti(k) ) * (grav) / ( Cp_layer(k) * (plev(2*k+3) - plev(2*k+1)) )
                        end do

                        open (unit = 90, file = 'heating_calc.txt', status = 'replace')
                                do k = 1,L_LAYERS
                                        write(90, ' (  500(X, E12.5) ) ') plev(2*k+2), fmnetv(k), fmneti(k), fmnetv(k) - fmnetv(k+1), fmneti(k+1) - fmneti(k), Cp_layer(k)
                                end do
                        close(90)
            
                        open(11, status = 'replace', file = 'tdt_lw.txt')
                                write (11, ' (  1000(X, E10.3) ) ') pfull
                                write (11, ' (  1000(X, E10.3) ) ') heatingi
                        close(11)
            
                        ! Adjust top layer
                        heatingv(1) = ( nfluxtopv - fmnetv(2) ) * grav / ( Cp_layer(1) * (plev(5) - plev(3)) )
                        heatingi(1) = ( fmneti(2) - nfluxtopi ) * grav / ( Cp_layer(1) * (plev(5) - plev(3)) )
             
!-----------------------------------------------------------------------------------------
! Put everything on the GCM grid
!-----------------------------------------------------------------------------------------		

                        toa_netv(i,j)  = nfluxtopv
                        netv(i,j,:)    = fmnetv
                        upv(i,j,:)   = fluxupv
                        dnv(i,j,:)   = fluxdnv
                        tdt_sw(i,j,:)  = heatingv
                        sw_spec(i,j,:) = sw_spectrum/dwnv/pi

                        toa_neti(i,j)  = nfluxtopi
                        neti(i,j,:)   = fmneti
                        upi(i,j,:)   = fluxupi
                        dni(i,j,:)   = fluxdni
                        tdt_lw(i,j,:)  = heatingi
                        lw_spec(i,j,:) = lw_spectrum/dwni/pi
                        
                        tdt(i,j,:)   = tdt_sw(i,j,:) + tdt_lw(i,j,:)
                        sfc_flux(i,j)  = fmnetv(L_LAYERS+1) + fluxdni(L_LAYERS+1)

!-----------------------------------------------------------------------------------------
! End loops over latitude (j) and longitude (i)
!-----------------------------------------------------------------------------------------	

                end do
        end do

END SUBROUTINE planetary_radiation

! ========================================================================================

SUBROUTINE planetary_radiation_end
!NLEDIT - astronomy not called 
!	call astronomy_end

        deallocate (ckcv)
        deallocate (ckci)
        deallocate (ckc_temp)
        deallocate (ckc_pres)
        deallocate (ckc_gwtv)
        deallocate (ckc_gwti) 
        deallocate (ciav)
        deallocate (ciai)
        deallocate (cia_temp)
        deallocate (ref_gases)
        deallocate (tracer_indices)
        module_is_initialized = .false.
        
        call haze_end

END SUBROUTINE planetary_radiation_end

! ========================================================================================
END MODULE planetary_radiation_driver_mod
