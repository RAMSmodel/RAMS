subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rte_rrtmgp_clouds stopping"
    error stop 1
  end if
end subroutine stop_on_err

! ----------------------------------------------------------------------------------
subroutine rte_rrtmgp_init()

  use rte_rrtmgp_rams
  use micro_prm, only:hucmfile
  use micphys, only:iaerorad
  use mo_load_cloud_coefficients, &
                             only: load_cld_lutcoeff
  use mo_load_aerosol_coefficients, &
                             only: load_aero_lutcoeff
  use mo_load_coefficients,  only: load_and_init
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_aerosol_optics,     only: ty_aerosol_optics

  implicit none

  CHARACTER(80) :: filename

  call stop_on_err(gas_concs%init(gas_names))
  !Initialize some gas concentrations that don't change in time
  !If you want them to change in time, they can easily be set along
  !with water vapor in the driver below.
  !ch4 and n2o are required by the code.
  call stop_on_err(gas_concs%set_vmr('ch4', 0.0))
  call stop_on_err(gas_concs%set_vmr('n2o', 0.0))
  call stop_on_err(gas_concs%set_vmr('o2 ', 0.209))
  call stop_on_err(gas_concs%set_vmr('co2', 420.e-6))

  filename = trim(hucmfile)//'/../RTE/rrtmgp-data-lw-g128-210809.nc'
  call load_and_init(k_dist_lw, filename, gas_concs)
  filename = trim(hucmfile)//'/../RTE/rrtmgp-data-sw-g112-210809.nc'
  call load_and_init(k_dist_sw, filename, gas_concs)

  filename = trim(hucmfile)//'/../RTE/mic2rrtmgp_lw.nc'
  call load_cld_lutcoeff (cloud_optics_lw, filename)
  filename = trim(hucmfile)//'/../RTE/mic2rrtmgp_sw.nc'
  call load_cld_lutcoeff (cloud_optics_sw, filename)

  if (iaerorad == 1) then
    filename = trim(hucmfile)//'/../RTE/aero2rrtmgp_lw.nc'
    call load_aero_lutcoeff (aerosol_optics_lw, filename)
    filename = trim(hucmfile)//'/../RTE/aero2rrtmgp_sw.nc'
    call load_aero_lutcoeff (aerosol_optics_sw, filename)
  endif

end subroutine rte_rrtmgp_init
! ----------------------------------------------------------------------------------
subroutine rte_rrtmgp_driver(nrad,ncat,aerocat,mu0,alb,p_lay,t_lay,rh_lay, &
                             z_lay,z_lev, &
                             rv_lay,t_sfc,wat_path,reff,rcat, &
                             acon,arad,atype,aodt, &
                             o3_lay,flux_up_sw,flux_dn_sw,flux_up_lw,flux_dn_lw,fthsw,fthlw)
  use micphys, only:iaerorad
  use rte_rrtmgp_rams !cloud_optics, k_dist
  use mo_rte_kind,           only: wp, i8, wl
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_aerosol_optics,     only: ty_aerosol_optics
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_source_functions,   only: ty_source_func_lw
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_heating_rates,      only: compute_heating_rate
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  use mo_rte_config,         only: rte_config_checks
  implicit none
  ! ----------------------------------------------------------------------------------
  ! Variables
  ! ----------------------------------------------------------------------------------
  ! Arrays: dimensions (col, lay)
  integer :: nrad, ncat, aerocat
  integer, parameter :: ncol=1
  real(wp), dimension(ncol,nrad) :: p_lev, t_lev
  real(wp), dimension(ncol,nrad-1) :: p_lay, t_lay, rv_lay, o3_lay, co2_lay, ch4_lay, n2o_lay, o2_lay
  real(wp), dimension(nrad) :: z_lev
  real(wp), dimension(nrad-1) :: z_lay
  real(wp), dimension(ncol,nrad-1,ncat) ::  wat_path, reff
  integer, dimension(ncol,nrad-1,ncat) :: rcat

  ! Aerosol only
  real(wp), dimension(ncol,nrad-1) :: rh_lay
  real(wp), dimension(ncol,nrad-1,aerocat) :: acon, arad
  integer, dimension(ncol,aerocat) :: atype
  real(wp), dimension(ncol) :: aodt

  !
  ! Longwave only
  !
  real(wp), dimension(ncol) :: t_sfc
  real(wp), dimension(:,:),   allocatable :: emis_sfc ! First dimension is band
  !
  ! Shortwave only
  !
  real(wp), dimension(ncol) :: mu0,alb
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Source functions
  !
  !   Longwave
  type(ty_source_func_lw), save         :: lw_sources
  !   Shortwave
  real(wp), dimension(:,:), allocatable :: toa_flux
  !
  ! Output variables
  !
  real(wp), dimension(ncol,nrad), target :: flux_up, flux_dn, flux_dir
  real, dimension(nrad) :: flux_up_sw, flux_up_lw, flux_dn_sw, flux_dn_lw
  real, dimension(nrad-1) ::  fthsw, fthlw
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  class(ty_optical_props_arry), &
                 allocatable :: atmos, clouds, aerosols
  type(ty_fluxes_broadband)  :: fluxes

  !
  ! Inputs to RRTMGP
  !
  integer :: k,top_at_1,iswlw
  logical :: is_sw, is_lw

  integer  :: nlay, nbnd, ngpt
  integer  :: icol, ilay, ibnd, iloop, igas

  logical :: use_luts = .true.

  ! Local variables
  integer :: nlm
  real(wp)                      :: H,dz,dzl
  real(wp), dimension(ncol,nrad-1) :: hrate,temp
  real :: start, finish
  ! NAR OpenMP CPU directives in compatible with OpenMP GPU directives
  !!$omp threadprivate( lw_sources, toa_flux, flux_up, flux_dn, flux_dir )
  ! ----------------------------------------------------------------------------------
  ! Code
  ! ----------------------------------------------------------------------------------
  !

  nlay = nrad - 1 

  ! Calculate level pressure and temperature
  do k = 2,nrad-1
     !Assuming exponential decrease between layers
     dzl = z_lay(k)-z_lay(k-1)
     H = dzl/log(p_lay(1,k-1)/p_lay(1,k))
     dz = z_lev(k)-z_lay(k-1)
     p_lev(1,k) = p_lay(1,k-1)*exp(-dz/H)

     !Assuming linear changes between layers
     t_lev(1,k) = t_lay(1,k-1) + (t_lay(1,k)-t_lay(1,k-1))/dzl*dz
  enddo
  ! Set surface values
  dzl = z_lay(2)-z_lay(1)
  dz = z_lev(1)-z_lay(1)
  H = dzl/log(p_lay(1,1)/p_lay(1,2))
  p_lev(1,1) = p_lay(1,1)*exp(-dz/H)
  t_lev(1,1) = t_sfc(1)
  ! Set top of atmosphere values
  dzl = z_lay(nlay) - z_lay(nlay-1)
  dz = z_lev(nrad)-z_lay(nlay)
  H = dzl/log(p_lay(1,nlay-1)/p_lay(1,nlay))
  p_lev(1,nrad) = p_lay(1,nlay)*exp(-dz/H)
  t_lev(1,nrad) = t_lay(1,nlay) + (t_lay(1,nlay)-t_lay(1,nlay-1))/dzl*dz

  !NEED VOLUME MIXING RATIO - Adele, yes, rv_lay is a volume mixing ratio 
  call stop_on_err(gas_concs%set_vmr('h2o', rv_lay))
  call stop_on_err(gas_concs%set_vmr('o3', o3_lay))

  ! ----------------------------------------------------------------------------
  ! load data into classes
  ! don't need gas_concentrations, just the list of available gases
 do iswlw = 1,2
  if (iswlw==1) then
     is_sw = .true.
     is_lw = .false.
  else
     is_sw = .false.
     is_lw = .true.
  endif
 
  !
  ! ----------------------------------------------------------------------------
  !
  ! Problem sizes
  !
  if (is_sw) then
    nbnd = k_dist_sw%get_nband()
    ngpt = k_dist_sw%get_ngpt()
  else
    nbnd = k_dist_lw%get_nband()
    ngpt = k_dist_lw%get_ngpt()
  endif

  top_at_1 = 0
  if(p_lay(1,1) < p_lay(1,nlay)) top_at_1 = 1

  ! ----------------------------------------------------------------------------
  ! LW calculations neglect scattering; SW calculations use the 2-stream approximation
  !   Here we choose the right variant of optical_props.
  !
  if(allocated(atmos)) deallocate(atmos)
  if(allocated(clouds)) deallocate(clouds)
  if(allocated(aerosols)) deallocate(aerosols)

  if(is_sw) then
    allocate(ty_optical_props_2str::atmos)
    allocate(ty_optical_props_2str::clouds)
    if(iaerorad==1) allocate(ty_optical_props_2str::aerosols)
  else
    allocate(ty_optical_props_1scl::atmos)
    allocate(ty_optical_props_1scl::clouds)
    if(iaerorad==1) allocate(ty_optical_props_1scl::aerosols)
  end if

  ! Clouds optical props are defined by band
  if(is_sw)  call stop_on_err(clouds%init(k_dist_sw%get_band_lims_wavenumber()))
  if(is_lw)  call stop_on_err(clouds%init(k_dist_lw%get_band_lims_wavenumber()))
  ! Aerosols optical props are defined by band
  if (iaerorad==1) then
    if(is_sw)  call stop_on_err(aerosols%init(k_dist_sw%get_band_lims_wavenumber()))
    if(is_lw)  call stop_on_err(aerosols%init(k_dist_lw%get_band_lims_wavenumber()))
  endif
  !
  ! Allocate arrays for the optical properties themselves.
  !
  select type(atmos)
  !I've hard coded here that LW is 1scl and SW is 2str
    class is (ty_optical_props_1scl)
      !$acc enter data copyin(atmos)
      call stop_on_err(atmos%alloc_1scl(ncol, nlay, k_dist_lw))
    class is (ty_optical_props_2str)
      call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist_sw))
    class default
      call stop_on_err("rte_rrtmgp_clouds: Don't recognize the kind of optical properties ")
  end select
  select type(clouds)
    class is (ty_optical_props_1scl)
      call stop_on_err(clouds%alloc_1scl(ncol, nlay))
    class is (ty_optical_props_2str)
      call stop_on_err(clouds%alloc_2str(ncol, nlay))
    class default
      call stop_on_err("rte_rrtmgp_clouds: Don't recognize the kind of optical properties ")
  end select
  if (iaerorad==1) then
    select type(aerosols)
      class is (ty_optical_props_1scl)
        call stop_on_err(aerosols%alloc_1scl(ncol, nlay))
      class is (ty_optical_props_2str)
        call stop_on_err(aerosols%alloc_2str(ncol, nlay))
      class default
        call stop_on_err("rte_rrtmgp_clouds: Don't recognize the kind of optical properties ")
    end select
  endif
  ! ----------------------------------------------------------------------------
  !  Boundary conditions depending on whether the k-distribution being supplied
  !   is LW or SW
  if(is_sw) then
    ! toa_flux is threadprivate
    allocate(toa_flux(ncol, ngpt))
    allocate(sfc_alb_dir(nbnd,ncol), sfc_alb_dif(nbnd,ncol))
    !Adele - we can do better
    sfc_alb_dir = alb(1)
    sfc_alb_dif = alb(1) 
  else
    ! lw_sources is threadprivate
    call stop_on_err(lw_sources%alloc(ncol, nlay, k_dist_lw))
    allocate(emis_sfc(nbnd, ncol))
    !Adele - couple with surface scheme here
    t_sfc = t_lev(1, merge(nlay+1, 1, top_at_1==1))
    emis_sfc = 1._wp!0.98_wp
  end if
  ! ----------------------------------------------------------------------------
  ! Multiple iterations for big problem sizes, and to help identify data movement
  !   For CPUs we can introduce OpenMP threading over loop iterations
  !
  !
!  call rte_config_checks(logical(.false., wl))

  if(is_sw)call stop_on_err(                                      &
   cloud_optics_sw%cloud_optics(wat_path, reff, rcat, clouds))
  if(is_lw)call stop_on_err(                                      &
   cloud_optics_lw%cloud_optics(wat_path, reff, rcat, clouds))

  aodt = 0.
  if(is_sw .and. iaerorad.eq.1)call stop_on_err(                                      &
   aerosol_optics_sw%aerosol_optics(acon, arad, atype, rh_lay, aerosols, aodt))
  if(is_lw .and. iaerorad.eq.1)call stop_on_err(                                      &
   aerosol_optics_lw%aerosol_optics(acon, arad, atype, rh_lay, aerosols))
    !
    ! Solvers
    !
  fluxes%flux_up => flux_up(:,:)
  fluxes%flux_dn => flux_dn(:,:)
  if(is_lw) then
    call stop_on_err(k_dist_lw%gas_optics(p_lay, p_lev, &
                                       t_lay, t_sfc, &
                                       gas_concs,    &
                                       atmos,        &
                                       lw_sources,   &
                                       tlev = t_lev))
    call stop_on_err(clouds%increment(atmos))
    if(iaerorad.eq.1)call stop_on_err(aerosols%increment(atmos))
    call stop_on_err(rte_lw(atmos, top_at_1, &
                            lw_sources,      &
                            emis_sfc,        &
                            fluxes))

    flux_up_lw = flux_up(1,:)
    flux_dn_lw = flux_dn(1,:)

    call stop_on_err(compute_heating_rate(fluxes%flux_up,fluxes%flux_dn,p_lev,hrate))
    fthlw = hrate(1,:)

  else
    !even though we're not using flux_dir, this statement also serves to allocate
    !fluxes%flux_dn_dir. Do not remove. 
    fluxes%flux_dn_dir => flux_dir(:,:)

    call stop_on_err(k_dist_sw%gas_optics(p_lay, p_lev, &
                                       t_lay,        &
                                       gas_concs,    &
                                       atmos,        &
                                       toa_flux))
    call stop_on_err(clouds%delta_scale())
    call stop_on_err(clouds%increment(atmos))
    if (iaerorad.eq.1) then
      call stop_on_err(aerosols%delta_scale())
      call stop_on_err(aerosols%increment(atmos))
    endif
    call stop_on_err(rte_sw(atmos, top_at_1, &
                            mu0,   toa_flux, &
                            sfc_alb_dir, sfc_alb_dif, &
                            fluxes))
    flux_up_sw = flux_up(1,:)
    flux_dn_sw = flux_dn(1,:)

    call stop_on_err(compute_heating_rate(fluxes%flux_up,fluxes%flux_dn,p_lev,hrate))    
    fthsw = hrate(1,:)
  end if
 enddo
 deallocate(toa_flux) 
end subroutine rte_rrtmgp_driver
