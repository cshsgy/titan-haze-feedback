! ========================================================================================
! RECONSTRUCTED minimal modules: haze_mod + read_clim_mod
! ========================================================================================
! The uploaded TAM radiation code references haze_mod and read_clim_mod, whose
! source was NOT included in the upload or the data archive (only Mac-compiled
! .mod files were present).  These minimal reimplementations cover the
! PRESCRIBED-data path used by namelist_huygens_tropo_spinup
! (haze_data='presc', clim_gas=''), under which only haze_init / haze_end /
! saturate are actually called; the climatology routines are present only to
! satisfy the `use, only:` imports and stop if ever invoked.
!
! NOTE: saturate() is a documented pass-through (no condensation) and the
! climatology routines are stubs -- this reproduces the radiative configuration
! with prescribed haze + gas profiles, not the full TAM condensation/climatology
! physics.  Replace with the original haze.F90 / read_clim.F90 for that.
! ========================================================================================

MODULE haze_mod
        implicit none
        private
        public :: haze_init, hazeclim_init, interp_haze_ls, clim_haze_profile, haze_end
contains

! Read prescribed haze optical properties (tau, g, ssa) from the
! INPUT/DATA/<presc_file><band>_{tau,g,ssa}.txt files and interpolate them onto
! the model (plev, wno) grid.  File format (statclim):
!   line1: header; line2: nwns; line3: npls; line4: wn(nwns); line5: pl(npls);
!   then npls rows of nwns values.
SUBROUTINE haze_init(tau, g, ssa, presc_file, plev, wno, band, doclimin)
        real, dimension(:,:), intent(out)   :: tau, g, ssa
        character(len=*), intent(in)        :: presc_file
        real, dimension(:), intent(in)      :: plev, wno
        character(len=*), intent(in)        :: band
        integer, intent(in), optional       :: doclimin
        call read_presc(trim(presc_file)//trim(band)//'_tau.txt', plev, wno, tau)
        call read_presc(trim(presc_file)//trim(band)//'_g.txt',   plev, wno, g)
        call read_presc(trim(presc_file)//trim(band)//'_ssa.txt', plev, wno, ssa)
END SUBROUTINE haze_init

SUBROUTINE read_presc(fname, plev, wno, out)
        character(len=*), intent(in)        :: fname
        real, dimension(:), intent(in)      :: plev, wno
        real, dimension(:,:), intent(out)   :: out
        integer :: nwns, npls, i, j, ios
        real, allocatable :: wn(:), pl(:), dat(:,:)
        open(unit=77, file='INPUT/DATA/'//trim(fname), status='old', iostat=ios)
        if (ios /= 0) then
                print *, 'haze_init: cannot open INPUT/DATA/'//trim(fname)
                stop 1
        end if
        read(77,*)            ! header
        read(77,*) nwns
        read(77,*) npls
        allocate(wn(nwns), pl(npls), dat(npls,nwns))
        read(77,*) (wn(i), i=1,nwns)
        read(77,*) (pl(j), j=1,npls)
        do j=1,npls
                read(77,*) (dat(j,i), i=1,nwns)
        end do
        close(77)
        call bilin(pl, wn, dat, plev, wno, out)
        deallocate(wn, pl, dat)
END SUBROUTINE read_presc

! Bilinear interpolation in (log pressure, wavenumber); edge-clamped.
SUBROUTINE bilin(pl, wn, dat, pout, wout, out)
        real, dimension(:), intent(in)    :: pl, wn, pout, wout
        real, dimension(:,:), intent(in)  :: dat
        real, dimension(:,:), intent(out) :: out
        integer :: ip, iw, jp, jw
        real :: fp, fw, lp
        do ip = 1, size(pout)
                lp = log(max(pout(ip), 1.0e-30))
                jp = locate(log(pl), lp)
                fp = (lp - log(pl(jp))) / (log(pl(jp+1)) - log(pl(jp)))
                fp = max(0.0, min(1.0, fp))
                do iw = 1, size(wout)
                        jw = locate(wn, wout(iw))
                        fw = (wout(iw) - wn(jw)) / (wn(jw+1) - wn(jw))
                        fw = max(0.0, min(1.0, fw))
                        out(ip, iw) = (1-fp)*(1-fw)*dat(jp,jw)   + (1-fp)*fw*dat(jp,jw+1) &
                                    +     fp*(1-fw)*dat(jp+1,jw) +     fp*fw*dat(jp+1,jw+1)
                end do
        end do
END SUBROUTINE bilin

! Index j such that x(j) <= xi < x(j+1), clamped to [1, n-1] (x ascending).
INTEGER FUNCTION locate(x, xi)
        real, dimension(:), intent(in) :: x
        real, intent(in) :: xi
        integer :: j, n
        n = size(x)
        locate = 1
        do j = 1, n-1
                if (xi >= x(j)) locate = j
        end do
        locate = max(1, min(n-1, locate))
END FUNCTION locate

! --- climatology routines: not used in prescribed mode (stubs) ---
SUBROUTINE hazeclim_init(wnoi, wnov, lats, plev, fname)
        real, dimension(:), intent(in) :: wnoi, wnov, lats, plev
        character(len=*), intent(in)   :: fname
        print *, 'hazeclim_init: climatology haze not supported in this build'; stop 1
END SUBROUTINE hazeclim_init

SUBROUTINE interp_haze_ls()
        print *, 'interp_haze_ls: not supported in this build'; stop 1
END SUBROUTINE interp_haze_ls

SUBROUTINE clim_haze_profile(tau, j, band, wno, plev)
        real, dimension(:,:), intent(inout) :: tau
        integer, intent(in) :: j, band
        real, dimension(:), intent(in) :: wno, plev
        print *, 'clim_haze_profile: climatology haze not supported in this build'; stop 1
END SUBROUTINE clim_haze_profile

SUBROUTINE haze_end()
END SUBROUTINE haze_end

END MODULE haze_mod


MODULE read_clim_mod
        implicit none
        private
        public :: readclim_MOD, saturate
contains

! Trace-gas climatology read: not used with prescribed gas profiles (stub).
SUBROUTINE readclim_MOD(clim_data, clpress, lat_grid, clls, clmolidx, clim_gas)
        real, allocatable, intent(out)    :: clim_data(:,:,:,:)
        real, allocatable, intent(out)    :: clpress(:), clls(:)
        real, dimension(:), intent(in)    :: lat_grid
        integer, allocatable, intent(out) :: clmolidx(:)
        character(len=*), dimension(:), intent(in) :: clim_gas
        print *, 'readclim_MOD: gas climatology not supported in this build'; stop 1
END SUBROUTINE readclim_MOD

! Minimal stability cap standing in for the original condensation routine.
! The real read_clim saturate() limits each condensable to its saturation
! mixing ratio; without it, organic trace gases run away in the cold tropopause
! and the integration diverges.  Here we cap each gas at its (approximate)
! Clausius-Clapeyron saturation VMR -- enough for a stable run, but NOT the
! original quantitative condensation physics.
SUBROUTINE saturate(mixing_ratio, gasname, tlev, plev)
        real, dimension(:), intent(inout) :: mixing_ratio
        character(len=*), intent(in)      :: gasname
        real, dimension(:), intent(in)    :: tlev, plev
        integer :: k
        real :: L_R, T0, p0, psat, vmr_sat
        ! approximate (latent heat / R [K], ref T0 [K], ref psat p0 [Pa])
        select case (trim(gasname))
        case ('CH4');  L_R = 1010.;  T0 = 90.7;  p0 = 1.17e4
        case ('C2H6'); L_R = 1938.;  T0 = 90.0;  p0 = 1.1e0
        case ('C2H2'); L_R = 2613.;  T0 = 192.;  p0 = 1.28e5
        case ('C2H4'); L_R = 1880.;  T0 = 104.;  p0 = 1.3e3
        case ('HCN');  L_R = 4000.;  T0 = 260.;  p0 = 1.0e5
        case default;  mixing_ratio = min(mixing_ratio, 1.0e-3); return
        end select
        do k = 1, size(mixing_ratio)
                psat = p0 * exp(-L_R * (1.0/max(tlev(k),40.) - 1.0/T0))
                vmr_sat = psat / max(plev(k), 1.0e-30)
                mixing_ratio(k) = min(mixing_ratio(k), vmr_sat)
        end do
END SUBROUTINE saturate

END MODULE read_clim_mod
