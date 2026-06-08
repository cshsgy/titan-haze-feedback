MODULE planetary_radiation_mod
! ========================================================================================
! GENERALIZED RADIATIVE TRANSFER MODEL FOR PLANETARY GCMs
! Copyright (c) Juan M. Lora
! ========================================================================================
! Written by J.M. Lora, based on the TAM radiation code
! Some code is derived from the public Ames legacy Mars GCM radiation code
! ========================================================================================
! Contains subroutines needed for radiative transfer model, called from the 
! planetary_radiation_driver module.
! ========================================================================================

        !use fms_mod,             only: error_mesg, FATAL
        !use constants_mod,       only: Pi

        implicit none
        private
        save

        public :: setspv, setspi, get_taukcoeff, get_tauCIA, optc, sfluxv, sfluxi, error_mesg
    
        real :: pi = 3.14159
        real :: FATAL = 0.
    
! ========================================================================================
        contains
! ========================================================================================

SUBROUTINE error_mesg(inp1, inp2, inp3)

        character(len=*) :: inp1, inp2
        real, optional ::  inp3
        print *, "ERROR - "
        print *, inp1
        print *, ''
        print *, inp2
        print *, ''

END SUBROUTINE error_mesg

SUBROUTINE setspv(wnov,dwnv,model_solar,tauray,bwnv,solar_spec_file,rayleigh_file)
! ========================================================================================
!
!       INPUT:
!       bwnv    			- Wavenumber [cm^-1] at edges of intervals
!       solar_spec_file	                - Filename of solar (stellar) spectrum to be used
!
!       OUTPUT:
!       wnov    			- Wavenumbers at interval centers
!       dwnv    			- Widths of intervals
!       model_solar  			- Solar flux (W/m2) in each interval
!       tauray  			- Pressure-independent Rayleigh scattering in each interval
!
! ========================================================================================

	real, dimension(:), intent(in)               :: bwnv
	real, dimension(:), intent(out)              :: wnov, dwnv, tauray
        character(len=128), intent(in)               :: solar_spec_file, rayleigh_file
        integer                                      :: m, n, ios, solar_file_size, rayleigh_wn_size, k
	real, dimension(:), allocatable, intent(out) :: model_solar
	real, dimension(:), allocatable              :: check_wn, solarf, solarf_edge, taurayf, rayleigh_wn
	real                                         :: dummy, file_left

!-----------------------------------------------------------------------------------------
! Set up midpoint wavenumbers and wavenumber widths
!-----------------------------------------------------------------------------------------

        do m = 1,size(bwnv,1)-1
                wnov(m) = 0.5 * (bwnv(m+1)+bwnv(m))
                dwnv(m) = bwnv(m+1) - bwnv(m)
        end do

!-----------------------------------------------------------------------------------------
! Set up spectral interval and Rayleigh scattering
! Interpolate from file wavenumber grid to model wavenumber grid
!-----------------------------------------------------------------------------------------

        open(10,file='INPUT/DATA/'//trim(solar_spec_file),status='old',iostat=ios)
        if (ios .ne. 0) print *, 'planetary_radiation_mod setspv: Could not open'//trim(solar_spec_file)
        read(10,*) !Read header text
        read(10, *) solar_file_size
        allocate(solarf(solar_file_size))
        allocate(solarf_edge(solar_file_size-1))
        allocate(check_wn(solar_file_size))
        allocate(model_solar(size(wnov)))
        
        do m = 1, solar_file_size
                read(10,*) check_wn(m), solarf(m)
        end do
        close(10)

        ! Calculate band edges for solar file
        solarf_edge(1) = check_wn(1) - 0.5*(check_wn(2) - check_wn(1))
        solarf_edge(2) = (check_wn(1) - solarf_edge(1)) + check_wn(1)

        do m = 3, solar_file_size-1
                solarf_edge(m) = (check_wn(m-1) - solarf_edge(m-1)) + check_wn(m-1)
        end do

        if (solar_file_size .eq. size(wnov)) then
                if (ALL(check_wn .eq. wnov)) then
                        !print *, 'Solar spectrum wavenumber grid identical to model.'
                        model_solar = solarf
                else
                        !print *, 'Interpolating solar spectrum to model wavenumber grid'
                        do m = 1, size(bwnv)-1
                                file_left = bwnv(m)
                                model_solar(m) = 0.0
                                do n = 1, solar_file_size-1
                                        !print *, 'Doing solarf_edge = ', solarf_edge(n)
                                        if (solarf_edge(n) .ge. bwnv(m+1)) then
                                                !print *, 'option 1'
                                                model_solar(m) = model_solar(m) + (bwnv(m+1) - file_left)*solarf(n)
                                                exit
                                        elseif ( (solarf_edge(n) .le. bwnv(m+1)) .and. (solarf_edge(n) .ge. bwnv(m)) ) then
                                                !print *, 'option 2'
                                                model_solar(m) = model_solar(m) + (solarf_edge(n) - file_left) * solarf(n)
                                                file_left = solarf_edge(n)
                                        end if
                                end do
                        end do
                end if
        else
                !print *, 'Interpolating solar spectrum to model wavenumber grid'
                do m = 1, size(bwnv)-1
                        file_left = bwnv(m)
                        model_solar(m) = 0.0
                        do n = 1, solar_file_size-1
                                !print *, 'Doing solarf_edge = ', solarf_edge(n)
                                if (solarf_edge(n) .ge. bwnv(m+1)) then
                                        !print *, 'option 1'
                                        model_solar(m) = model_solar(m) + (bwnv(m+1) - file_left)*solarf(n)
                                        exit
                                elseif ( (solarf_edge(n) .le. bwnv(m+1)) .and. (solarf_edge(n) .ge. bwnv(m)) ) then
                                        !print *, 'option 2'
                                        model_solar(m) = model_solar(m) + (solarf_edge(n) - file_left) * solarf(n)
                                        file_left = solarf_edge(n)
                                end if
                        end do
                end do
        end if

        ! interpolate Rayleigh as well, at least check wn grid as in k table	
        ! print *, 'Starting rayleigh'
        open(10,file='INPUT/DATA/'//trim(rayleigh_file),status='old',iostat=ios)
        if (ios .ne. 0) print *, 'planetary_radiation_mod setspv: Could not open Rayleigh.txt'
        read(10,*) !Read header text
        read(10, *) rayleigh_wn_size
        ALLOCATE(rayleigh_wn(rayleigh_wn_size))
        ALLOCATE(taurayf(rayleigh_wn_size))
        do m = 1,rayleigh_wn_size
                read(10,*) rayleigh_wn(m), taurayf(m)
        end do
        close(10)

        do k = 1, SIZE(wnov)
                if (wnov(k) .lt. rayleigh_wn(1)) then
                        tauray(k) = taurayf(1)
                elseif (wnov(k) .gt. rayleigh_wn(SIZE(rayleigh_wn))) then
                        tauray(k) = taurayf(SIZE(taurayf)) + &
                                    (taurayf(SIZE(taurayf))-taurayf(SIZE(taurayf)-1))/(rayleigh_wn(SIZE(rayleigh_wn))-rayleigh_wn(SIZE(rayleigh_wn)-1)) &
                                    * (wnov(k)-rayleigh_wn(SIZE(rayleigh_wn)))
                        tauray(k) = tauray(k)*1.2
                else
                        do m = 1, rayleigh_wn_size
                                if (rayleigh_wn(m) .ge. wnov(k)) then
                                        tauray(k)=(wnov(k)-rayleigh_wn(m-1))*(taurayf(m)-taurayf(m-1))/(rayleigh_wn(m)-rayleigh_wn(m-1))+taurayf(m-1)
                                        exit
                                end if
                        end do
                end if
        end do

        open(10,file="model_interp_ray.txt")
        write(10, ' ( 500(X, E10.3) ) ') wnov
        write(10, ' ( 500(X, E10.3) ) ') tauray
        close(10)

        DEALLOCATE(solarf)
        DEALLOCATE(solarf_edge)
        DEALLOCATE(check_wn)
        DEALLOCATE(rayleigh_wn)
        DEALLOCATE(taurayf)

END SUBROUTINE setspv

! ========================================================================================

SUBROUTINE setspi(wnoi,dwni,planckir,bwni)
! ========================================================================================
!  Set up IR spectral intervals
!
!       INPUT:
!       bwni				- Wavenumber [cm^-1] at edges of intervals
!
!       OUTPUT:
!       wnoi				- Wavenumbers at interval centers
!       dwni    			- Widths of intervals
!       planckir			- Integral of the Planck function in each interval
!
! ========================================================================================
	real, dimension(:), intent(in)          :: bwni
	real, dimension(:), intent(out)         :: wnoi, dwni
	real, dimension(:,:), intent(inout)     :: planckir

	real    :: c1 = 3.741832D-16 !W m-2 [see Goody & Young (2nd Ed.)]
	real    :: c2 = 1.438786D-2  !m K   [see Goody & Young (2nd Ed.)]
	real    :: x(12) = [  -0.981560634246719D0,  -0.904117256370475D0, &
                          -0.769902674194305D0,  -0.587317954286617D0, &
                          -0.367831498998180D0,  -0.125233408511469D0, &
                           0.125233408511469D0,   0.367831498998180D0, &
                           0.587317954286617D0,   0.769902674194305D0, &
                           0.904117256370475D0,   0.981560634246719D0 ]
        real    :: w(12) = [   0.047175336386512D0,   0.106939325995318D0, &
                           0.160078328543346D0,   0.203167426723066D0, &
                           0.233492536538355D0,   0.249147045813403D0, &
                           0.249147045813403D0,   0.233492536538355D0, &
                           0.203167426723066D0,   0.160078328543346D0, &
                           0.106939325995318D0,   0.047175336386512D0 ]
        real    :: a, b, T, ans, y
        integer :: m, n, s

!-----------------------------------------------------------------------------------------
! Set up midpoint wavenumbers and wavenumber widths
!-----------------------------------------------------------------------------------------

        do m = 1,size(bwni,1)-1
                wnoi(m) = 0.5 * (bwni(m+1)+bwni(m))
                dwni(m) = bwni(m+1) - bwni(m)
        end do
        ! Manually set wnoi(1) to be band edge.  bwni(0) is 1, not 0, so the calculated wnoi is off by 0.5
        wnoi(1) = wnoi(1) - 0.5

!-----------------------------------------------------------------------------------------
! Compute integral of B(T) divided by interval wavenumber (W m-2 wavnum-1) for each band
!-----------------------------------------------------------------------------------------

        do m = 1,size(bwni)-1
                a = 1.0D-2/bwni(m+1)
                b = 1.0D-2/bwni(m)
                do n = 200,9000 ! CK: changed 500 to 200
                        T = dble(n)/1.0D+1
                        ans = 0.0D0
                        do s = 1,12
                                y = ((b-a)/2.0)*x(s) + ((b+a)/2.0)
                                ans = ans + w(s)*c1/(y**5*(exp(c2/(y*T))-1.0D0))
                        end do
                        planckir(m,n-199) = ans*((b-a)/2.0)/(pi*dwni(m)) ! CK: Originally 499
                end do
        end do

END SUBROUTINE setspi

! ========================================================================================

SUBROUTINE get_taukcoeff(tau,col_abund,temp,pres,k_coeff,tlim,plim)
! ========================================================================================
!  Calculate gas opacities with correlated k coefficients
!
!       INPUT
!       col_abund			- Column abundance at plev [km amg]
!       temp				- Temperature tlev
!       pres				- Pressure plev
!       k_coeff       	                - k coefficients table
!       tlim				- Temperatures for k coefficient table
!       plim				- Pressures for k coefficient table
!
!       OUTPUT
!       tau				- Optical depth due to absorber (added to input tau)
!
! ========================================================================================

	real, dimension(:), intent(in)          :: col_abund, temp, pres
	real, dimension(:,:,:,:), intent(in)    :: k_coeff
	real, dimension(:), intent(in)          :: tlim, plim
	real, dimension(:,:,:), intent(inout)   :: tau

	real, dimension(size(k_coeff,2),size(k_coeff,3),size(k_coeff,4))  :: p_interp
	real, dimension(size(k_coeff,3),size(k_coeff,4))                  :: k_abs
        integer                                                           :: l, i

!-----------------------------------------------------------------------------------------
! Interpolate k coefficient to appropriate pressure and temperature, calculate opacities
!-----------------------------------------------------------------------------------------

        do l = 1, size(col_abund,1)
                if (pres(l) .le. plim(1)) then
                        p_interp = k_coeff(1,:,:,:)
                else if (pres(l) .ge. plim(size(plim,1))) then
                        p_interp = k_coeff(size(plim,1),:,:,:)
                else
                        do i = 1, size(plim,1)
                                if (pres(l) .le. plim(i)) then
                                        p_interp = k_coeff(i,:,:,:) - log(plim(i)/pres(l))/log(plim(i)/plim(i-1)) * &
                                                   (k_coeff(i,:,:,:) - k_coeff(i-1,:,:,:))
                                        exit
                                end if
                        end do
                end if

                if (temp(l) .le. tlim(1)) then
                        k_abs = p_interp(1,:,:)
                else if (temp(l) .ge. tlim(size(tlim,1))) then
                        k_abs = p_interp(size(tlim,1),:,:)
                else
                        do i = 1, size(tlim,1)
                                if (temp(l) .le. tlim(i)) then
                                        k_abs = p_interp(i,:,:) - (tlim(i) - temp(l))/(tlim(i) - tlim(i-1)) * &
                                                (p_interp(i,:,:) - p_interp(i-1,:,:))
                                        exit
                                end if
                        end do
                end if

                tau(l,:,:) = tau(l,:,:) + k_abs(:,:) * col_abund(l)

                !Why is there an upper limit on tau?
                !where (tau(l,:,:) .gt. 1.D2)
                !	tau(l,:,:) = 1.D2
                !end where
        end do

END SUBROUTINE get_taukcoeff

! ========================================================================================

SUBROUTINE get_tauCIA(tau,col_abund,temp,cia_fits,nfits,tlim)
! ========================================================================================
!  Calculate opacities due to CIA
!
!       INPUT
!       col_abund			- Column abundance [molec^2/cm^5]
!       temp				- Temperature tlev
!       cia_fits 			- CIA fits table
!       nfits				- Number of fits points in table
!       tlim				- Temperatures for fits table
!
!       OUTPUT
!       tau				- Optical depth due to CIA (added to input tau)
!
! ========================================================================================

	real, dimension(:), intent(in)          :: col_abund, temp
	real, dimension(:,:,:), intent(in)      :: cia_fits
	real, dimension(:), intent(in)          :: tlim
        integer, intent(in)                     :: nfits
	real, dimension(:,:), intent(inout)     :: tau

	real, dimension(size(tau,2))            :: trans, trans_lo, trans_hi
        integer                                 :: l, i, t, n

!-----------------------------------------------------------------------------------------
! Interpolate to appropriate temperature, calculate transmission, calculate optical depth
!-----------------------------------------------------------------------------------------

        do l = 1, size(col_abund,1)
                if (temp(l) .le. tlim(1)) then
                        trans = 0.0
                        do i = 1, nfits
                                trans(:) = trans(:) + cia_fits(1,:,i) * &
                                           exp(-cia_fits(1,:,i+nfits) * col_abund(l))
                        end do
                        WHERE (trans(:) .gt. 1.) trans(:) = 1.
                        tau(l,:) = tau(l,:) - log(trans(:))
                        cycle
                else if (temp(l) .ge. tlim(size(tlim,1))) then
                        trans = 0.0
                        do i = 1, nfits
                                trans(:) = trans(:) + cia_fits(size(tlim,1),:,i) * &
                                           exp(-cia_fits(size(tlim,1),:,i+nfits) * col_abund(l))
                        end do
                        WHERE (trans(:) .gt. 1.) trans(:) = 1.
                        tau(l,:) = tau(l,:) - log(trans(:))
                        cycle
                else
                        do t = 1, size(tlim,1)
                                if (temp(l) .le. tlim(t)) then 
                                        trans_lo = 0.0
                                        trans_hi = 0.0
                                        ! Does linear interpolation on cia table temperature grid
                                        ! Calculates T at bounding CIA temperature points
                                        ! nfits number of exp
                                        do i = 1, nfits
                                                trans_lo(:) = trans_lo(:) + cia_fits(t-1,:,i) * &
                                                              exp(-cia_fits(t-1,:,i+nfits)*col_abund(l))
                                                trans_hi(:) = trans_hi(:) + cia_fits(t,:,i) * &
                                                              exp(-cia_fits(t,:,i+nfits)*col_abund(l))
                                        end do
                                        trans(:) = trans_hi - (tlim(t) - temp(l))/(tlim(t) - tlim(t-1)) * (trans_hi - trans_lo)
                                        ! CK 7/2/24: Model crashes when trans_hi is smaller than machine precision, which 
                                        !            occurs when the temperature at the model bottom gets too cold. To deal
                                        !            with this, I added a check for trans = 0, and if it finds one, it 
                                        !            replaces the zero with a super small number.
                                        do n = 1,size(trans,1)
                                                if (trans(n) .eq. 0.) then
                                                        trans(n) = 1.0D-307
                                                endif 
                                        end do

                                        WHERE (trans(:) .gt. 1.) trans(:) = 1.

                                        tau(l,:) = tau(l,:) - log(trans(:))
                                        exit
                                end if
                        end do
                end if
        end do

END SUBROUTINE get_tauCIA

! ========================================================================================

SUBROUTINE optc(dtau,tau,taucum,wbar,cosb,tau_gas,tau_ray_in,tau_cia_in,tau_haze_in, &
                tau_cloud_in,ssa_haze_in,ssa_cloud_in,g_haze_in,g_cloud_in)
! ========================================================================================
!  Calculate cumulative optical properties 
!
!       INPUT
!       tau_gas			- Gas optical depth
!       tau_ray_in		- Rayleigh optical depth
!       tau_cia_in		- CIA optical depth
!       tau_haze_in		- Haze optical depth
!       tau_cloud_in		- Cloud optical depth
!       ssa_haze_in		- Haze single scattering albedo
!       ssa_cloud_in		- Cloud single scattering albedo
!       g_haze_in		- Haze asymmetry parameter
!       g_cloud_in		- Cloud asymmetry parameter
!
!       OUTPUT
!       dtau			- Optical depth for each layer
!       tau			- Optical depth from top to layer bottom
!       taucum			- Optical depth from top to level bottom
!       wbar			- Single scattering albedo
!       cosb			- Asymmetry parameter
!
! ========================================================================================

	real, dimension(:,:,:), intent(in)                :: tau_gas
	real, dimension(:,:), intent(in), optional        :: tau_ray_in, tau_cia_in
	real, dimension(:,:), intent(in), optional        :: tau_haze_in,tau_cloud_in
	real, dimension(:,:), intent(in), optional        :: ssa_haze_in, ssa_cloud_in
	real, dimension(:,:), intent(in), optional        :: g_haze_in, g_cloud_in
	real, dimension(:,:,:), intent(out)               :: dtau, tau, taucum, wbar, cosb

	real, dimension(size(tau_gas,1),size(tau_gas,2))                  :: tau_ray, tau_cia
	real, dimension(size(tau_gas,1),size(tau_gas,2))                  :: tau_haze, tau_cloud
	real, dimension(size(tau_gas,1),size(tau_gas,2))                  :: ssa_haze, ssa_cloud
	real, dimension(size(tau_gas,1),size(tau_gas,2))                  :: g_haze, g_cloud
	real, dimension(size(tau_gas,1),size(tau_gas,2),size(tau_gas,3))  :: dtauk
        integer                                                           :: l, k, nw, ng, layers

!-----------------------------------------------------------------------------------------
! If no values passed in, set them to zero
!-----------------------------------------------------------------------------------------
        if (present(tau_ray_in)) then
                tau_ray = tau_ray_in
        else
                tau_ray = 0.0
        endif
        if (present(tau_cia_in)) then
                tau_cia = tau_cia_in
        else
                tau_cia = 0.0
        endif
        if (present(tau_haze_in)) then
                tau_haze = tau_haze_in
                ssa_haze = ssa_haze_in
                g_haze   = g_haze_in
        else
                tau_haze = 0.0
                ssa_haze = 0.0
                g_haze   = 0.0
        endif
        if (present(tau_cloud_in)) then
                tau_cloud = tau_cloud_in
                ssa_cloud = ssa_cloud_in
                g_cloud   = g_cloud_in
        else
                tau_cloud = 0.0
                ssa_cloud = 0.0
                g_cloud   = 0.0
        endif
!-----------------------------------------------------------------------------------------
! Determine the total level opacity at each spectral interval nw and gaussian point ng
!-----------------------------------------------------------------------------------------

        do ng = 1, size(tau_gas,3)
                dtauk(:,:,ng) = tau_gas(:,:,ng) + tau_ray + tau_cia + tau_haze + tau_cloud
        end do

!-----------------------------------------------------------------------------------------
! Calculate layer opacities, single scattering albedos, asymmetry factors
!-----------------------------------------------------------------------------------------

        layers = (size(tau_gas,1) - 3)/2

        do ng = 1, size(tau_gas,3)
                do nw = 1, size(tau_gas,2)
                        do l = 1, layers
                                k = 2*l + 1
                                dtau(l,nw,ng) = dtauk(k-1,nw,ng) + dtauk(k,nw,ng)
                                if (dtau(l,nw,ng) .gt. 1.D-9) then
                                        wbar(l,nw,ng) = (tau_ray(k-1,nw) + tau_ray(k,nw) + &
                                                                         ssa_haze(k-1,nw)*tau_haze(k-1,nw) + &
                                                                         ssa_haze(k,nw)*tau_haze(k,nw) + &
                                                                         ssa_cloud(k-1,nw)*tau_cloud(k-1,nw)+ &
                                                                         ssa_cloud(k,nw)*tau_cloud(k,nw)) / dtau(l,nw,ng)
                                        cosb(l,nw,ng) = (g_haze(k-1,nw)*ssa_haze(k-1,nw)*tau_haze(k-1,nw)   +&
                                                                         g_haze(k,nw)*ssa_haze(k,nw)*tau_haze(k,nw) + &
                                                                         g_cloud(k-1,nw)*ssa_cloud(k-1,nw)*tau_cloud(k-1,nw)+&
                                                                         g_cloud(k,nw)*ssa_cloud(k,nw)*tau_cloud(k,nw)) / &
                                                                         ( wbar(l,nw,ng) * dtau(l,nw,ng) )
                                        if (wbar(l,nw,ng) .eq. 0.0) cosb(l,nw,ng) = 0.0
                                        if (wbar(l,nw,ng) .eq. 1.0) wbar(l,nw,ng) = 0.9999999
                                else
                                        dtau(l,nw,ng) = 1.D-9    !1.D-9 NAL edit 11/17/21
                                        wbar(l,nw,ng) = 0.0
                                        cosb(l,nw,ng) = 0.0
                                end if
                        end do
!-----------------------------------------------------------------------------------------
! Special bottom layer
!-----------------------------------------------------------------------------------------
            
                        l = layers + 1
                        k = 2*l + 1
                        dtau(l,nw,ng) = dtauk(k,nw,ng)
                        if (dtau(l,nw,ng) .gt. 1.D-9) then
                                wbar(l,nw,ng) = (tau_ray(k,nw) + ssa_haze(k,nw)*tau_haze(k,nw) + &
                                                                 ssa_cloud(k,nw)*tau_cloud(k,nw)) / dtau(l,nw,ng)
                                cosb(l,nw,ng) = (g_haze(k,nw)*ssa_haze(k,nw)*tau_haze(k,nw) + &
                                                                 g_cloud(k,nw)*ssa_cloud(k,nw)*tau_cloud(k,nw)) / &
                                                                 ( wbar(l,nw,ng) * dtau(l,nw,ng) )
                                if (wbar(l,nw,ng) .eq. 0.0) cosb(l,nw,ng) = 0.0
                                if (wbar(l,nw,ng) .eq. 1.0) wbar(l,nw,ng) = 0.9999999
                        else
                                dtau(l,nw,ng) = 1.D-9     !1.D-9 NAL edit 11/17/21
                                wbar(l,nw,ng) = 0.0
                                cosb(l,nw,ng) = 0.0
                        end if

!-----------------------------------------------------------------------------------------
! Calculate total extinction optical depths
!-----------------------------------------------------------------------------------------

                        tau(1,nw,ng) = 0.0
                        do l = 1, layers+1
                                tau(l+1,nw,ng) = tau(l,nw,ng) + dtau(l,nw,ng)
                        end do
                        taucum(1,nw,ng) = 0.0
                        do k = 2, size(tau_gas,1)
                                taucum(k,nw,ng) = taucum(k-1,nw,ng) + dtauk(k,nw,ng)
                        end do
                end do
        end do


END SUBROUTINE optc

! ========================================================================================

SUBROUTINE sfluxv(nfluxtop,fmnet,fluxup,fluxdn,difft,sw_spectrum, &
                  dtau,tau,taucum,gweight,albedo,ssa,g,cosz,sol)
! ========================================================================================
!  Calculate radiative fluxes (solar wavelengths)
!
!       INPUT
!       dtau			- Opacity of layers
!       tau			- Opacity from top to bottom of layers
!       taucum			- Opacity from top to bottom of levels
!       gweight			- Gaussian weights for k coefficients
!       albedo			- Surface albedo
!       ssa			- Single scattering albedo
!       g			- Asymmetry parameter
!       cosz			- Cosine of incidence angle
!       sol		        - Solar flux in each band at TOA
!
!       OUTPUT
!       nfluxtop	        - Net flux at TOA
!       fmnet		        - Net flux at bottom of layers
!       fluxup			- Upward flux at layers
!       fluxdn        	        - Downward flux at layers
!       difft			- Diffuse downward flux
!       sw_spectrum		- Spectrum of reflected flux at TOA
!
! ========================================================================================

	real, dimension(:,:,:), intent(in)          :: dtau, tau, taucum, ssa, g
	real, dimension(:), intent(in)              :: gweight, sol
	real, intent(in)                            :: albedo, cosz
	real, dimension(:), intent(out)             :: fmnet, fluxup, fluxdn, sw_spectrum
	real, intent(out)                           :: nfluxtop, difft

	real, dimension(size(dtau,1))               :: fmup, fmdn
	real, dimension(size(dtau,2),size(dtau,3))  :: detau
	real                                        :: btop,bsurf,eterm,diff,flxup,flxdn, maxexp
        integer                                     :: nw, ng, l

!-----------------------------------------------------------------------------------------
! Initialize calculation
!-----------------------------------------------------------------------------------------

        nfluxtop    = 0.0
        difft       = 0.0
        fmnet       = 0.0
        fluxup      = 0.0
        fluxdn      = 0.0
        sw_spectrum = 0.0
        maxexp      = 35.0

!-----------------------------------------------------------------------------------------
! Begin loops over spectral intervals nw and Gaussian points ng
!-----------------------------------------------------------------------------------------

        do nw = 1, size(dtau,2)
                do ng = 1, size(dtau,3)

!-----------------------------------------------------------------------------------------
! Get scaled optical depth at the surface; set boundary conditions
!-----------------------------------------------------------------------------------------

                        call getdetau( dtau(:,nw,ng),tau(:,nw,ng),taucum(:,nw,ng), &
                                       ssa(:,nw,ng),g(:,nw,ng),detau(nw,ng) )
   
                        btop  = 0.0
                        eterm = min(detau(nw,ng)/cosz, maxexp)
                        bsurf = albedo * cosz * sol(nw) * exp(-eterm)

!-----------------------------------------------------------------------------------------
! Calculate fluxes for each interval, then cumulative fluxes
!-----------------------------------------------------------------------------------------

                        call gfluxv( dtau(:,nw,ng),tau(:,nw,ng),taucum(:,nw,ng), &
                                     ssa(:,nw,ng),g(:,nw,ng),cosz,sol(nw),albedo,btop,bsurf, &
                                     fmup,fmdn,diff,flxup,flxdn )

                        nfluxtop        = nfluxtop + (flxdn - flxup) * gweight(ng)
                        difft           = difft + diff * gweight(ng)
                        sw_spectrum(nw) = sw_spectrum(nw) + flxup * gweight(ng)

                        do l = 1, size(dtau,1)
                                fmnet(l)  = fmnet(l) + (fmdn(l) - fmup(l)) * gweight(ng)
                                fluxup(l) = fluxup(l) + fmup(l) * gweight(ng)
                                fluxdn(l) = fluxdn(l) + fmdn(l) * gweight(ng)
                        end do
                end do
        end do

END SUBROUTINE sfluxv

! ========================================================================================

SUBROUTINE sfluxi(nfluxtop,fmnet,fluxup,fluxdn,lw_spectrum, &
                  plev,tlev,dtau,taucum,gweight,dwni,albedo,ssa,g,cosz,planckir)
! ========================================================================================
!  Calculate radiative fluxes (thermal IR wavelengths)
!
!     INPUT
!     plev			- Pressure at levels
!     tlev			- Temperature at levels
!     dtau			- Opacity of layers
!     taucum			- Opacity from top to bottom of levels
!     gweight			- Gaussian weights for k coefficients
!     dwni			- Band widths in wavenumber
!     albedo			- Surface albedo
!     ssa			- Single scattering albedo
!     g				- Asymmetry parameter
!     cosz			- Cosine of incidence angle
!     planckir			- Integral of the Planck function in each interval
!
!     OUTPUT
!     nfluxtop			- Net flux at TOA
!     fmnet			- Net flux at bottom of layers
!     fluxup			- Upward flux at layers
!     fluxdn        	        - Downward flux at layers
!     lw_spectrum		- Spectrum of outgoing flux at TOA
!
! ========================================================================================

	real, dimension(:,:,:), intent(in)      :: dtau, taucum, ssa, g
	real, dimension(:,:), intent(in)        :: planckir
	real, dimension(:), intent(in)          :: plev, tlev, gweight, dwni
	real, intent(in)                        :: albedo, cosz
	real, dimension(:), intent(out)         :: fmnet, fluxup, fluxdn, lw_spectrum
	real, intent(out)                       :: nfluxtop

	real, dimension(size(dtau,1))           :: fmup, fmdn
	real                                    :: bsurf,btop,pltop,ftopup
        integer                                 :: ntt, nts, nw, ng, l, debug

!-----------------------------------------------------------------------------------------
! Initialize calculation
!-----------------------------------------------------------------------------------------

        nfluxtop    = 0.0
        fmnet       = 0.0
        fluxup      = 0.0
        fluxdn      = 0.0
        lw_spectrum = 0.0

!-----------------------------------------------------------------------------------------
! Set boundary temperature indices
!-----------------------------------------------------------------------------------------

        ntt = planck_index(tlev(2), size(planckir,2))
        nts = planck_index(tlev(size(tlev,1)), size(planckir,2))

!-----------------------------------------------------------------------------------------
! Begin loops over spectral intervals and Gaussian points; set boundary conditions
!-----------------------------------------------------------------------------------------

        do nw = 1, size(dtau,2)
                bsurf = (1. - albedo) * planckir(nw,nts)
                pltop = planckir(nw,ntt)
                do ng = 1, size(dtau,3)
                        btop = ( 1.0 - exp((-dtau(1,nw,ng)*plev(2)/(plev(4)-plev(2)))/cosz) )*pltop

!-----------------------------------------------------------------------------------------
! Calculate fluxes for each interval, then cumulative fluxes
!-----------------------------------------------------------------------------------------
                        call gfluxi( tlev,nw,planckir,dtau(:,nw,ng),taucum(:,nw,ng), &
                                     ssa(:,nw,ng),g(:,nw,ng),cosz,albedo,btop,bsurf, &
                                     ftopup,fmup,fmdn )

                        nfluxtop        = nfluxtop + ftopup * dwni(nw) * gweight(ng)
                        lw_spectrum(nw) = lw_spectrum(nw) + ftopup * dwni(nw) * gweight(ng)
                        do l = 1, size(dtau,1)
                                fmnet(l)  = fmnet(l) + (fmup(l) - fmdn(l)) * dwni(nw) * gweight(ng)
                                fluxup(l) = fluxup(l) + fmup(l) * dwni(nw) * gweight(ng)
                                fluxdn(l) = fluxdn(l) + fmdn(l) * dwni(nw) * gweight(ng)
                        end do
                end do
        end do

END SUBROUTINE sfluxi

! ========================================================================================

SUBROUTINE gfluxv(dtdel,tdel,taucumin,wdel,cdel,cosz,   &
                  f0pi,alb,btop,bsurf,fmidp,fmidm,diff, &
                  fluxup,fluxdn)
! ==================================================================================
!
!       INPUT
!       dtdel				- Optical depth of layers
!       tdel    			- Column optical depth of layers
!       taucumin 			- Column optical depth of levels
!       wdel     			- Single scattering albedo
!       cdel     			- Asymmetry factor
!       cosz     			- Cosine of the incidence angle
!       f0pi     			- Incident solar flux
!       alb      			- Surface albedo
!       btop     			- Upper boundary condition
!       bsurf    			- Lower boundary condition
!
!       OUTPUT
!       fmidp    			- Upward flux at layer midpoints
!       fmidm    			- Downward flux at layer midpoints
!       diff     			- Diffuse downward flux at surface
!       fluxup   			- Upward flux at toa
!       fluxdn   			- Downward flux at toa
!
! ==================================================================================

	real, dimension(:), intent(in)    :: dtdel,wdel,cdel,tdel,taucumin
  	real, intent(in)                  :: cosz,f0pi,alb,btop,bsurf
  	real, dimension(:), intent(out)   :: fmidp, fmidm
  	real, intent(out)                 :: fluxup,fluxdn,diff

  	real, dimension(size(dtdel,1))    :: w0, cosbar, dtau
  	real, dimension(size(tdel,1))     :: tau
  	real, dimension(size(dtdel,1))    :: lambda,alpha,xk1,xk2,g1,g2,g3,gama,cp,cm
  	real, dimension(size(dtdel,1))    :: cpm1,cmm1,e1,e2,e3,e4,exptrm
  	real, dimension(size(taucumin,1)) :: taucum
  	real                              :: g4,denom,am,ap,taumax,taumid,cpmid,cmmid,em,ep
  	real                              :: factor
        integer                           :: nayer, l, k

!-----------------------------------------------------------------------------------------
!  Delta-Eddington Scaling
!-----------------------------------------------------------------------------------------

        taumax    = 35.0
        nayer     = size(dtdel,1)
        factor    = 1.0 - wdel(1) * cdel(1)**2
        tau(1)    = tdel(1) * factor
        taucum(1) = 0.0
        taucum(2) = taucumin(2) * factor
        taucum(3) = taucum(2) + (taucumin(3)-taucumin(2))*factor
        do l = 1,nayer-1
                factor      = 1.0 - wdel(l) * cdel(l)**2
                w0(l)       = wdel(l)*(1.0 - cdel(l)**2) / factor
                cosbar(l)   = cdel(l) / (1.0 + cdel(l))
                dtau(l)     = dtdel(l) * factor
                tau(l+1)    = tau(l) + dtau(l)
                k           = 2*(l+1)
                taucum(k)   = tau(l+1)
                taucum(k+1) = taucum(k) + (taucumin(k+1)-taucumin(k))*factor
        end do
        l       = nayer
        !Adjusted versions, why is tau not adjusted as well -it's adjusted a few lines above
        factor        = 1.0 - wdel(l) * cdel(l)**2
        w0(l)         = wdel(l)*(1.0 - cdel(l)**2) / factor
        cosbar(l)     = cdel(l) / (1.0 + cdel(l))
        dtau(l)       = dtdel(l) * factor
        tau(l+1)      = tau(l) + dtau(l)
        taucum(2*l+1) = tau(l+1)

!-----------------------------------------------------------------------------------------
! Eddington method / Quadrature method
!-----------------------------------------------------------------------------------------

        do l = 1,nayer
                alpha(l) = sqrt( (1.0 - w0(l)) / (1.0 - w0(l)*cosbar(l)) )
                g1(l) = 0.25 * (7.0 - w0(l) * (4.0+3.0*cosbar(l)))
                g2(l) = -0.25 * (1.0 - w0(l) * (4.0-3.0*cosbar(l)))
                g3(l) = 0.25 * (2.0 - 3.0*cosbar(l)*cosz) 
                lambda(l) = sqrt( g1(l)**2 - g2(l)**2 )  !WARNING - g1, g2, are identical, lambda is 0
                gama(l)  = (g1(l) - lambda(l)) / g2(l)
        end do

        do l = 1,nayer
                g4 = 1.0 - g3(l)
                denom = lambda(l)**2 - 1./cosz**2
                if (denom .eq. 0.0) denom = 1.D-10
                am = f0pi*w0(l)*(g4 * (g1(l) + 1./cosz) + g2(l)*g3(l) ) / denom
                ap = f0pi*w0(l)*(g3(l) * (g1(l) - 1./cosz) + g2(l)*g4 ) / denom

!-----------------------------------------------------------------------------------------
! cpm and cmm are c+ and c- terms at top of layer; cp, cm are terms at bottom
!-----------------------------------------------------------------------------------------

                cpm1(l) = ap*exp(-min(tau(l)/cosz,taumax))
                cmm1(l) = am*exp(-min(tau(l)/cosz,taumax))
                cp(l) = ap*exp(-min(tau(l+1)/cosz,taumax))
                cm(l) = am*exp(-min(tau(l+1)/cosz,taumax))
        end do

!-----------------------------------------------------------------------------------------
! Calculate exponential terms for tridiagonal rotated layer method
!-----------------------------------------------------------------------------------------

        do l = 1,nayer
                exptrm(l) = min(taumax,lambda(l)*dtau(l))
                ep    = exp(exptrm(l))
                em    = 1.0/ep
                e1(l) = ep + gama(l)*em
                e2(l) = ep - gama(l)*em
                e3(l) = gama(l)*ep + em
                e4(l) = gama(l)*ep - em
        end do

!-----------------------------------------------------------------------------------------
! Call the tridiagonal dolver
!-----------------------------------------------------------------------------------------

        call dsolver(nayer,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4,btop,bsurf,alb,xk1,xk2)

!-----------------------------------------------------------------------------------------
! Calculate the fluxes at layer midpoints
!-----------------------------------------------------------------------------------------

        do l = 1,nayer-1
                if (l .eq. 1) then
                        exptrm(l) = min(taumax,lambda(l)*(taucum(2*l+1)-taucum(2*l)))
                else
                        exptrm(l) = min(taumax,lambda(l)*(taucum(2*l+1)-taucum(2*l)))
                end if
                ep = exp(exptrm(l))
                em = 1.0/ep
                g4 = 1.0 - g3(l)
                denom = lambda(l)**2 - 1./cosz**2
                
                if (denom .eq. 0.0) denom = 1.d-10
                am = f0pi*w0(l)*(g4 * (g1(l) + 1./cosz) + g2(l)*g3(l) ) / denom
                ap = f0pi*w0(l)*(g3(l) * (g1(l) - 1./cosz) + g2(l)*g4 ) / denom

                taumid = taucum(2*l+1)
                cpmid = ap*exp(-min(taumid/cosz,taumax))
                cmmid = am*exp(-min(taumid/cosz,taumax))

                fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
                fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

!-----------------------------------------------------------------------------------------
! Add the direct flux to the downwelling term
!-----------------------------------------------------------------------------------------

                fmidm(l) = fmidm(l) + cosz*f0pi*exp(-min(taumid/cosz,taumax))
        end do

!-----------------------------------------------------------------------------------------
! Repeat for the top layer
!-----------------------------------------------------------------------------------------

        ep = 1.0
        em = 1.0
        g4 = 1.0-g3(1)
        denom = lambda(1)**2 - 1./cosz**2
        if (denom .eq. 0.0) denom = 1.d-10

        am = f0pi*w0(1)*(g4 * (g1(1) + 1./cosz) + g2(1)*g3(1) ) / denom
        ap = f0pi*w0(1)*(g3(1) * (g1(1) - 1./cosz) + g2(1)*g4 ) / denom

        cpmid = ap
        cmmid = am

        fluxup = xk1(1)*ep + gama(1)*xk2(1)*em + cpmid
        fluxdn = xk1(1)*ep*gama(1) + xk2(1)*em + cmmid

        fluxdn = fluxdn + cosz*f0pi*exp(-min(taucum(1)/cosz,taumax))

!-----------------------------------------------------------------------------------------
! Repeat for the bottom layer
!-----------------------------------------------------------------------------------------

        l = nayer
        exptrm(l) = min(taumax,lambda(l)*(taucum(size(taucum,1))-taucum(size(taucum,1)-1)))
        ep = exp(exptrm(l))
        em = 1.0/ep
        g4 = 1.0-g3(l)
        denom = lambda(l)**2 - 1./cosz**2
        if (denom .eq. 0.0) denom = 1.d-10

        am = f0pi*w0(l)*(g4 * (g1(l) + 1./cosz) + g2(l)*g3(l) ) / denom
        ap = f0pi*w0(l)*(g3(l) * (g1(l) - 1./cosz) + g2(l)*g4 ) / denom
  
        taumid = min(taucum(size(taucum,1)),taumax)
        cpmid = ap*exp(-min(taumid/cosz,taumax))
        cmmid = am*exp(-min(taumid/cosz,taumax))

        fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
        fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

        diff = fmidm(l)

        fmidm(l) = fmidm(l) + cosz*f0pi*exp(-min(taumid/cosz,taumax))

END SUBROUTINE gfluxv

! ==================================================================================

SUBROUTINE gfluxi(tlev,nw,planckir,dtau,taucum,w0,cosbar,cosz,    &
                  alb,btop,bsurf,ftopup,fmidp,fmidm)
! ==================================================================================
!
!       INPUT
!       tlev     			- temperature at levels
!       dtau     			- optical depth of layers
!       taucum   			- column optical depth of levels
!       w0       			- single scattering albedo
!       cosbar   			- asymmetry factor
!       cosz     			- cosine of the incidence angle
!       alb      			- surface albedo
!       btop     			- upper boundary condition
!       bsurf    			- lower boundary condition
!
!       OUTPUT
!       fmidp    			- upward flux at layer midpoints
!       fmidm    			- downward flux at layer midpoints
!       ftopup   			- upward flux at toa
!
! ==================================================================================

	real, dimension(:), intent(in)          :: dtau,w0,cosbar,tlev,taucum
	real, dimension(:,:), intent(in)        :: planckir
	real, intent(in)                        :: cosz,alb,btop,bsurf
        integer, intent(in)                     :: nw
  	real, dimension(:), intent(out)         :: fmidp, fmidm
  	real, intent(out)                       :: ftopup

  	real, dimension(size(dtau,1))           :: b0,b1,alpha,lambda,xk1,xk2
  	real, dimension(size(dtau,1))           :: gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4
  	real                                    :: term,cpmid,cmmid,em,ep,dtauk,taumax,fluxup,fluxdn
        integer                                 :: nayer, l, nt,nt2

!-----------------------------------------------------------------------------------------
! Hemispheric constant method
!-----------------------------------------------------------------------------------------
        taumax = 35.0
        nayer  = size(dtau,1)
       
        do l = 1,nayer-1
        
                alpha(l) = sqrt( (1.0 - w0(l)) / (1.0 - w0(l)*cosbar(l)) )
                lambda(l) = alpha(l) * (1.0 - w0(l) * cosbar(l)) / cosz

                nt  = planck_index(tlev(2*l), size(planckir,2))
                nt2 = planck_index(tlev(2*l+2), size(planckir,2))

                b0(l) = planckir(nw,nt)
                b1(l) = (planckir(nw,nt2) - planckir(nw,nt)) / dtau(l)     

        end do
!-----------------------------------------------------------------------------------------
! Bottom layer
!-----------------------------------------------------------------------------------------

        l = nayer
        alpha(l) = sqrt( (1.0 - w0(l)) / (1.0 - w0(l)*cosbar(l)) )
        lambda(l) = alpha(l) * (1.0 - w0(l) * cosbar(l)) / cosz
        nt  = planck_index(tlev(2*l), size(planckir,2))
        nt2 = planck_index(tlev(2*l+1), size(planckir,2))
        b0(l) = planckir(nw,nt)
        b1(l) = (planckir(nw,nt2) - planckir(nw,nt)) / dtau(l)
  
!-----------------------------------------------------------------------------------------
! cpm and cmm are c+ and c- terms at top of layer; cp, cm are terms at bottom
!-----------------------------------------------------------------------------------------

        do l = 1,nayer
                gama(l) = (1.0 - alpha(l)) / (1.0 + alpha(l))
                term     = cosz / (1.0 - w0(l) * cosbar(l))

                cp(l) = b0(l) + b1(l)*dtau(l) + b1(l)*term
                cm(l) = b0(l) + b1(l)*dtau(l) - b1(l)*term
                cpm1(l) = b0(l) + b1(l)*term
                cmm1(l) = b0(l) - b1(l)*term 
        end do

!-----------------------------------------------------------------------------------------
! Calculate exponential terms for tridiagonal rotated layer method
!-----------------------------------------------------------------------------------------

        do l = 1,nayer
                ep    = exp(min((lambda(l)*dtau(l)),taumax))
                em    = 1.0/ep
                e1(l) = ep + gama(l)*em
                e2(l) = ep - gama(l)*em
                e3(l) = gama(l)*ep + em
                e4(l) = gama(l)*ep - em
        end do

!-----------------------------------------------------------------------------------------
! Call the tridiagonal dolver
!-----------------------------------------------------------------------------------------
        
        call dsolver(nayer,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4,btop,bsurf,alb,xk1,xk2)

!-----------------------------------------------------------------------------------------
! Calculate the fluxes at layer midpoints
!-----------------------------------------------------------------------------------------

        do l = 1,nayer-1
                dtauk = taucum(2*l+1) - taucum(2*l)
                ep = exp(min((lambda(l)*dtauk),taumax))
                em = 1.0/ep
                term = cosz / (1.0 - w0(l)*cosbar(l))

                cpmid = b0(l) + b1(l)*dtauk + b1(l)*term    
                cmmid = b0(l) + b1(l)*dtauk - b1(l)*term

                fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
                fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

!-----------------------------------------------------------------------------------------
! Integrate over the hemisphere
!-----------------------------------------------------------------------------------------
                fmidp(l) = fmidp(l)*pi
                fmidm(l) = fmidm(l)*pi
        end do

!-----------------------------------------------------------------------------------------
! Repeat for the top layer
!-----------------------------------------------------------------------------------------

        ep = 1.0
        em = 1.0
        term = cosz / (1.0 - w0(1)*cosbar(1))
        cpmid = b0(1) + b1(1)*term    
        cmmid = b0(1) - b1(1)*term

        fluxup = xk1(1)*ep + gama(1)*xk2(1)*em + cpmid
        fluxdn = xk1(1)*ep*gama(1) + xk2(1)*em + cmmid

        ftopup = (fluxup-fluxdn)*pi

!-----------------------------------------------------------------------------------------
! Repeat for the bottom layer
!-----------------------------------------------------------------------------------------

        l = nayer
        ep = exp(min((lambda(l)*dtau(l)),taumax))
        em = 1.0/ep
        term = cosz / (1.0 - w0(l)*cosbar(l))
        cpmid = b0(l) + b1(l)*dtau(l) + b1(l)*term    
        cmmid = b0(l) + b1(l)*dtau(l) - b1(l)*term

        fmidp(l) = xk1(l)*ep + gama(l)*xk2(l)*em + cpmid
        fmidm(l) = xk1(l)*ep*gama(l) + xk2(l)*em + cmmid

        fmidp(l) = fmidp(l)*pi
        fmidm(l) = fmidm(l)*pi

END SUBROUTINE gfluxi

! ==================================================================================

SUBROUTINE dsolver(nl,gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4, &
                   btop,bsurf,alb,xk1,xk2)
! ==================================================================================
!  Solve for the coefficients of the two-stream solution
!
!       INPUT
!       nl       			- number of layers
!       cp       			- c+ at TOA
!       cm       			- c- at TOA
!       cpm1     			- c+ at tau (bottom)
!       cmm1     			- c- at tau (bottom)
!       ep       			- exp(lambda*dtau)
!       em       			- 1/ep
!       e1       			- ep + gama * em
!       e2       			- ep - gama * em
!       e3       			- gama * ep + em
!       e4       			- gama * ep - em
!       btop     			- upper boundary condition
!       bsurf    			- lower boundary condition
!
!       OUTPUT
!       xk1      			- coefficient for positive exp term
!       xk2      			- coefficient for negative exp term
!
! ==================================================================================

	real, dimension(:), intent(in)  :: gama,cp,cm,cpm1,cmm1,e1,e2,e3,e4
  	real, intent(in)                :: btop, bsurf, alb
        integer, intent(in)             :: nl
  	real, dimension(:), intent(out) :: xk1, xk2

  	real, dimension(2*nl)           :: af, bf, cf, df, xk
        integer                         :: i, l, lm2, lm1, n

        l     = 2*nl
        af(1) = 0.0
        bf(1) = gama(1) + 1
        cf(1) = gama(1) - 1
        df(1) = btop - cmm1(1)
        n     = 0
        lm2   = l-2

!-----------------------------------------------------------------------------------------
! Even terms
!-----------------------------------------------------------------------------------------

        do i = 2,lm2,2
                n     = n+1
                af(i) = (e1(n)+e3(n)) * (gama(n+1) - 1.0)
                bf(i) = (e2(n)+e4(n)) * (gama(n+1) - 1.0)
                cf(i) = 2.0 * (1.0 - gama(n+1)**2)
                df(i) = (gama(n+1) - 1.0) * (cpm1(n+1) - cp(n)) + &
                        (1.0 - gama(n+1)) * (cm(n) - cmm1(n+1))
        end do

!-----------------------------------------------------------------------------------------
! Odd terms
!-----------------------------------------------------------------------------------------
  
        n   = 0
        lm1 = l-1

        do i = 3,lm1,2
                n     = n+1
                af(i) = 2.0 * (1.0 - gama(n)**2)
                bf(i) = (e1(n)-e3(n)) * (1.0 + gama(n+1))
                cf(i) = (e1(n)+e3(n)) * (gama(n+1) - 1.0)
                df(i) = e3(n) * (cpm1(n+1) - cp(n)) + e1(n) * (cm(n) - cmm1(n+1))
        end do

        af(l) = e1(nl) - alb*e3(nl)
        bf(l) = e2(nl) - alb*e4(nl)
        cf(l) = 0.0
        df(l) = bsurf - cp(nl) + alb*cm(nl)
    
        call dtridgl (l,af,bf,cf,df,xk)
  
!-----------------------------------------------------------------------------------------
! Unmix coefficients
!-----------------------------------------------------------------------------------------

        do n = 1,nl
                xk1(n) = xk(2*n-1) + xk(2*n)
                xk2(n) = xk(2*n-1) - xk(2*n)

                if (xk2(n) .eq. 0.0) cycle
                if (abs(xk2(n)/(xk(2*n-1) + 1.d-20)) .lt. 1.d-30) xk2(n) = 0.0
        end do

END SUBROUTINE dsolver

! ==================================================================================

SUBROUTINE dtridgl (l,af,bf,cf,df,xk)
! ==================================================================================
!       Solves a system of tridiagonal matrix equations:
!       A(i)*X(i-1) + B(i)*X(i) + C(i)*X(i+1) = D(i)
! ==================================================================================

	real, dimension(:), intent(in)  :: af,bf,cf,df
  	real, dimension(:), intent(out) :: xk
        integer, intent(in)             :: l

  	real, dimension(l)        :: as, ds
  	real, dimension(size(bf)) :: bf_adj
  	real                      :: x
        integer                   :: i

        bf_adj = bf

        as(l) = af(l)/bf_adj(l)
        ds(l) = df(l)/bf_adj(l)

        do i = 2,l
                x = 1.0/(bf_adj(l+1-i) - cf(l+1-i)*as(l+2-i))
                as(l+1-i) = af(l+1-i)*x
                ds(l+1-i) = (df(l+1-i) - cf(l+1-i)*ds(l+2-i))*x
        end do

        xk(1) = ds(1)
        do i = 2,l
                xk(i) = ds(i) - as(i)*xk(i-1)
        end do

END SUBROUTINE dtridgl

! ========================================================================================

SUBROUTINE getdetau(dtdel,tdel,taucumin,wdel,cdel,detau)
! ========================================================================================

  	real, dimension(:), intent(in)          :: dtdel,wdel,cdel,tdel,taucumin
  	real, intent(out)                       :: detau

  	real, dimension(size(dtdel,1))          :: w0, cosbar, dtau
  	real, dimension(size(tdel,1))           :: tau
  	real, dimension(size(taucumin,1))       :: taucum
  	real                                    :: factor
        integer                                 :: nayer, l, k

        nayer = size(dtdel,1)

        factor    = 1.0 - wdel(1) * cdel(1)**2
        tau(1)    = tdel(1) * factor
        taucum(1) = 0.0
        taucum(2) = taucum(1) * factor
        taucum(3) = taucum(2) + (taucumin(3)-taucumin(2))*factor

        do l = 1,nayer-1
                factor      = 1.0 - wdel(l) * cdel(l)**2
                w0(l)       = wdel(l)*(1.0 - cdel(l)**2) / factor
                cosbar(l)   = cdel(l) / (1.0 + cdel(l))
                dtau(l)     = dtdel(l) * factor
                tau(l+1)    = tau(l) + dtau(l)
                k           = 2*(l+1)
                taucum(k)   = tau(l+1)
                taucum(k+1) = taucum(k) + (taucumin(k+1)-taucumin(k))*factor
        end do

        l           = nayer
        factor      = 1.0 - wdel(l) * cdel(l)**2
        w0(l)       = wdel(l)*(1.0 - cdel(l)**2) / factor
        cosbar(l)   = cdel(l) / (1.0 + cdel(l))
        dtau(l)     = dtdel(l) * factor
        tau(l+1)    = tau(l) + dtau(l)
        taucum(2*l+1) = tau(l+1)
        detau         = taucum(2*l+1)

END SUBROUTINE getdetau

! ========================================================================================

INTEGER FUNCTION planck_index(temp, ntab)

        real, intent(in)    :: temp
        integer, intent(in) :: ntab
        real                :: temp_clip, temp_max

        temp_max = 20.0 + 0.1*(real(ntab) - 1.0)
        temp_clip = min(temp_max, max(20.0, temp))
        planck_index = int(temp_clip*10.0 - 199.0)
        planck_index = max(1, min(ntab, planck_index))

END FUNCTION planck_index

! ========================================================================================

END MODULE planetary_radiation_mod
