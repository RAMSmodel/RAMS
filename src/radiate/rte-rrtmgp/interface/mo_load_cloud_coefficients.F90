module mo_load_cloud_coefficients
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_cloud_optics,  only: ty_cloud_optics
  use mo_simple_netcdf, only: read_field, read_string, var_exists, get_dim_size, &
                              write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_cld_lutcoeff, load_cld_padecoeff
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read cloud optical property LUT coefficients from NetCDF file
  !
  subroutine load_cld_lutcoeff(cloud_spec, cld_coeff_file)
    class(ty_cloud_optics),     intent(inout) :: cloud_spec
    character(len=*),           intent(in   ) :: cld_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nsize_cld, nsize_drz, nsize_rain, nsize_scol, nsize_hcol, &
               nsize_plate, nsize_ros, nsize_agg, nsize_hail

    real(wp), dimension(:,:), allocatable                :: band_lims_wvn
    ! Lookup table interpolation constants
    real(wp) :: radcld_lwr          ! cloud droplet size lower bound for interpolation
    real(wp) :: radcld_upr          ! cloud droplet size upper bound for interpolation
    real(wp) :: radcld_fac          ! constant for calculating interpolation indices for droplets
    real(wp) :: raddrz_lwr          ! drizzle particle size lower bound for interpolation
    real(wp) :: raddrz_upr          ! drizzle particle size upper bound for interpolation
    real(wp) :: raddrz_fac          ! constant for calculating interpolation indices for drizzle
    real(wp) :: radrain_lwr          ! rain particle size lower bound for interpolation
    real(wp) :: radrain_upr          ! rain particle size upper bound for interpolation
    real(wp) :: radrain_fac          ! constant for calculating interpolation indices for rain
    real(wp) :: radscol_lwr          ! solid column size lower bound for interpolation
    real(wp) :: radscol_upr          ! solid column size upper bound for interpolation
    real(wp) :: radscol_fac          ! constant for calculating interpolation indices for solid columns
    real(wp) :: radhcol_lwr          ! hollow column size lower bound for interpolation
    real(wp) :: radhcol_upr          ! hollow column size upper bound for interpolation
    real(wp) :: radhcol_fac          ! constant for calculating interpolation indices for hollow columns
    real(wp) :: radplate_lwr          ! hexagonal plate size lower bound for interpolation
    real(wp) :: radplate_upr          ! hexagonal plate size upper bound for interpolation
    real(wp) :: radplate_fac          ! constant for calculating interpolation indices for hexagonal plates
    real(wp) :: radros_lwr          ! rosette size lower bound for interpolation
    real(wp) :: radros_upr          ! rosette size upper bound for interpolation
    real(wp) :: radros_fac          ! constant for calculating interpolation indices for rosettes
    real(wp) :: radagg_lwr          ! aggregate size lower bound for interpolation
    real(wp) :: radagg_upr          ! aggregate size upper bound for interpolation
    real(wp) :: radagg_fac          ! constant for calculating interpolation indices for aggregates
    real(wp) :: radhail_lwr          ! hail size lower bound for interpolation
    real(wp) :: radhail_upr          ! hail size upper bound for interpolation
    real(wp) :: radhail_fac          ! constant for calculating interpolation indices for hail
    ! LUT coefficients
    real(wp), dimension(:,:), allocatable :: lut_extcld   ! extinction
    real(wp), dimension(:,:), allocatable :: lut_ssacld   ! single scattering albedo
    real(wp), dimension(:,:), allocatable :: lut_asycld   ! asymmetry parameter
    real(wp), dimension(:,:), allocatable :: lut_extdrz   ! 
    real(wp), dimension(:,:), allocatable :: lut_ssadrz   ! 
    real(wp), dimension(:,:), allocatable :: lut_asydrz   ! 
    real(wp), dimension(:,:), allocatable :: lut_extrain   !
    real(wp), dimension(:,:), allocatable :: lut_ssarain   !
    real(wp), dimension(:,:), allocatable :: lut_asyrain   !
    real(wp), dimension(:,:), allocatable :: lut_extscol   !
    real(wp), dimension(:,:), allocatable :: lut_ssascol   !
    real(wp), dimension(:,:), allocatable :: lut_asyscol   !
    real(wp), dimension(:,:), allocatable :: lut_exthcol   !
    real(wp), dimension(:,:), allocatable :: lut_ssahcol   !
    real(wp), dimension(:,:), allocatable :: lut_asyhcol   !
    real(wp), dimension(:,:), allocatable :: lut_extplate   !
    real(wp), dimension(:,:), allocatable :: lut_ssaplate   !
    real(wp), dimension(:,:), allocatable :: lut_asyplate   !
    real(wp), dimension(:,:), allocatable :: lut_extros   ! 
    real(wp), dimension(:,:), allocatable :: lut_ssaros   ! 
    real(wp), dimension(:,:), allocatable :: lut_asyros   ! 
    real(wp), dimension(:,:), allocatable :: lut_extagg   ! 
    real(wp), dimension(:,:), allocatable :: lut_ssaagg   ! 
    real(wp), dimension(:,:), allocatable :: lut_asyagg   ! 
    real(wp), dimension(:,:), allocatable :: lut_exthail   ! 
    real(wp), dimension(:,:), allocatable :: lut_ssahail  ! 
    real(wp), dimension(:,:), allocatable :: lut_asyhail   ! 
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_lutcoeff(): can't open file " // trim(cld_coeff_file))

    ! Read LUT coefficient dimensions
    nband     = get_dim_size(ncid,'nbands')
    nsize_cld = get_dim_size(ncid,'ntabcld')
    nsize_drz = get_dim_size(ncid,'ntabdrz')
    nsize_rain = get_dim_size(ncid,'ntabrain')
    nsize_scol = get_dim_size(ncid,'ntabscol')
    nsize_hcol = get_dim_size(ncid,'ntabhcol')
    nsize_plate = get_dim_size(ncid,'ntabplate')
    nsize_ros = get_dim_size(ncid,'ntabros')
    nsize_agg = get_dim_size(ncid,'ntabagg')
    nsize_hail = get_dim_size(ncid,'ntabhail')

    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    ! Read LUT constants
    ! I believe the facs aren't used
    radcld_lwr = read_field(ncid, 'radcld_lwr')
    radcld_upr = read_field(ncid, 'radcld_upr')
    radcld_fac = 1.!read_field(ncid, 'radcld_fac')
    raddrz_lwr = read_field(ncid, 'raddrz_lwr')
    raddrz_upr = read_field(ncid, 'raddrz_upr')
    raddrz_fac = 1.!read_field(ncid, 'raddrz_fac')
    radrain_lwr = read_field(ncid, 'radrain_lwr')
    radrain_upr = read_field(ncid, 'radrain_upr')
    radrain_fac = 1.!read_field(ncid, 'radrain_fac')
    radscol_lwr = read_field(ncid, 'radscol_lwr')
    radscol_upr = read_field(ncid, 'radscol_upr')
    radscol_fac = 1.!read_field(ncid, 'radscol_fac')
    radhcol_lwr = read_field(ncid, 'radhcol_lwr')
    radhcol_upr = read_field(ncid, 'radhcol_upr')
    radhcol_fac = 1.!read_field(ncid, 'radhcol_fac')
    radplate_lwr = read_field(ncid, 'radplate_lwr')
    radplate_upr = read_field(ncid, 'radplate_upr')
    radplate_fac = 1.!read_field(ncid, 'radplate_fac')
    radros_lwr = read_field(ncid, 'radros_lwr')
    radros_upr = read_field(ncid, 'radros_upr')
    radros_fac = 1.!read_field(ncid, 'radros_fac')
    radagg_lwr = read_field(ncid, 'radagg_lwr')
    radagg_upr = read_field(ncid, 'radagg_upr')
    radagg_fac = 1.!read_field(ncid, 'radagg_fac')
    radhail_lwr = read_field(ncid, 'radhail_lwr')
    radhail_upr = read_field(ncid, 'radhail_upr')
    radhail_fac = 1.!read_field(ncid, 'radhail_fac')

    ! Allocate cloud property lookup table input arrays
    allocate(lut_extcld(nsize_cld, nband), &
             lut_ssacld(nsize_cld, nband), &
             lut_asycld(nsize_cld, nband), &
             lut_extdrz(nsize_drz, nband), &
             lut_ssadrz(nsize_drz, nband), &
             lut_asydrz(nsize_drz, nband), &
             lut_extrain(nsize_rain, nband), &
             lut_ssarain(nsize_rain, nband), &
             lut_asyrain(nsize_rain, nband), &
             lut_extscol(nsize_scol, nband), &
             lut_ssascol(nsize_scol, nband), &
             lut_asyscol(nsize_scol, nband), &
             lut_exthcol(nsize_hcol, nband), &
             lut_ssahcol(nsize_hcol, nband), &
             lut_asyhcol(nsize_hcol, nband), &
             lut_extplate(nsize_plate, nband), &
             lut_ssaplate(nsize_plate, nband), &
             lut_asyplate(nsize_plate, nband), &
             lut_extros(nsize_ros, nband), &
             lut_ssaros(nsize_ros, nband), &
             lut_asyros(nsize_ros, nband), &
             lut_extagg(nsize_agg, nband), &
             lut_ssaagg(nsize_agg, nband), &
             lut_asyagg(nsize_agg, nband), &
             lut_exthail(nsize_hail, nband), &
             lut_ssahail(nsize_hail, nband), &
             lut_asyhail(nsize_hail, nband))

    ! Read LUT coefficients
    lut_extcld = read_field(ncid, 'lut_extcld',  nsize_cld, nband)
    lut_ssacld = read_field(ncid, 'lut_ssacld',  nsize_cld, nband)
    lut_asycld = read_field(ncid, 'lut_asycld',  nsize_cld, nband)
    lut_extdrz = read_field(ncid, 'lut_extdrz',  nsize_drz, nband)
    lut_ssadrz = read_field(ncid, 'lut_ssadrz',  nsize_drz, nband)
    lut_asydrz = read_field(ncid, 'lut_asydrz',  nsize_drz, nband)
    lut_extrain = read_field(ncid, 'lut_extrain',  nsize_rain, nband)
    lut_ssarain = read_field(ncid, 'lut_ssarain',  nsize_rain, nband)
    lut_asyrain = read_field(ncid, 'lut_asyrain',  nsize_rain, nband)
    lut_extscol = read_field(ncid, 'lut_extscol',  nsize_scol, nband)
    lut_ssascol = read_field(ncid, 'lut_ssascol',  nsize_scol, nband)
    lut_asyscol = read_field(ncid, 'lut_asyscol',  nsize_scol, nband)
    lut_exthcol = read_field(ncid, 'lut_exthcol',  nsize_hcol, nband)
    lut_ssahcol = read_field(ncid, 'lut_ssahcol',  nsize_hcol, nband)
    lut_asyhcol = read_field(ncid, 'lut_asyhcol',  nsize_hcol, nband)
    lut_extplate = read_field(ncid, 'lut_extplate',  nsize_plate, nband)
    lut_ssaplate = read_field(ncid, 'lut_ssaplate',  nsize_plate, nband)
    lut_asyplate = read_field(ncid, 'lut_asyplate',  nsize_plate, nband)
    lut_extros = read_field(ncid, 'lut_extros',  nsize_ros, nband)
    lut_ssaros = read_field(ncid, 'lut_ssaros',  nsize_ros, nband)
    lut_asyros = read_field(ncid, 'lut_asyros',  nsize_ros, nband)
    lut_extagg = read_field(ncid, 'lut_extagg',  nsize_agg, nband)
    lut_ssaagg = read_field(ncid, 'lut_ssaagg',  nsize_agg, nband)
    lut_asyagg = read_field(ncid, 'lut_asyagg',  nsize_agg, nband)
    lut_exthail = read_field(ncid, 'lut_exthail',  nsize_hail, nband)
    lut_ssahail = read_field(ncid, 'lut_ssahail',  nsize_hail, nband)
    lut_asyhail = read_field(ncid, 'lut_asyhail',  nsize_hail, nband)

    ncid = nf90_close(ncid)
    call stop_on_err(cloud_spec%load(band_lims_wvn,                      &
                                     radcld_lwr, radcld_upr, radcld_fac, &
                                     raddrz_lwr, raddrz_upr, raddrz_fac, &
                                     radrain_lwr, radrain_upr, radrain_fac, &
                                     radscol_lwr, radscol_upr, radscol_fac, &
                                     radhcol_lwr, radhcol_upr, radhcol_fac, &
                                     radplate_lwr, radplate_upr, radplate_fac, &
                                     radros_lwr, radros_upr, radros_fac, &
                                     radagg_lwr, radagg_upr, radagg_fac, &
                                     radhail_lwr, radhail_upr, radhail_fac, &
                                     lut_extcld, lut_ssacld, lut_asycld, &
                                     lut_extdrz, lut_ssadrz, lut_asydrz, &
                                     lut_extrain, lut_ssarain, lut_asyrain, &
                                     lut_extscol, lut_ssascol, lut_asyscol, &
                                     lut_exthcol, lut_ssahcol, lut_asyhcol, &
                                     lut_extplate, lut_ssaplate, lut_asyplate, &
                                     lut_extros, lut_ssaros, lut_asyros, &
                                     lut_extagg, lut_ssaagg, lut_asyagg, &
                                     lut_exthail, lut_ssahail, lut_asyhail))
  end subroutine load_cld_lutcoeff
  !--------------------------------------------------------------------------------------------------------------------
  ! read cloud optical property Pade coefficients from NetCDF file
  !
  subroutine load_cld_padecoeff(cloud_spec, cld_coeff_file)
    class(ty_cloud_optics),       intent(inout) :: cloud_spec
    character(len=*),             intent(in   ) :: cld_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, nrghice, nsizereg, ncoeff_ext, ncoeff_ssa_g, nbound

    ! Spectral discretization
    real(wp), dimension(:,:), allocatable :: band_lims_wvn

    ! Pade coefficients
    real(wp), dimension(:,:,:),   allocatable :: pade_extliq   ! extinction: liquid
    real(wp), dimension(:,:,:),   allocatable :: pade_ssaliq   ! single scattering albedo: liquid
    real(wp), dimension(:,:,:),   allocatable :: pade_asyliq   ! asymmetry parameter: liquid
    real(wp), dimension(:,:,:,:), allocatable :: pade_extice   ! extinction: ice
    real(wp), dimension(:,:,:,:), allocatable :: pade_ssaice   ! single scattering albedo: ice
    real(wp), dimension(:,:,:,:), allocatable :: pade_asyice   ! asymmetry parameter: ice

    ! Pade particle size regime boundaries
    real(wp),  dimension(:),       allocatable :: pade_sizreg_extliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_asyliq
    real(wp),  dimension(:),       allocatable :: pade_sizreg_extice
    real(wp),  dimension(:),       allocatable :: pade_sizreg_ssaice
    real(wp),  dimension(:),       allocatable :: pade_sizreg_asyice
    ! -----------------
    ! Open cloud optical property coefficient file
    if(nf90_open(trim(cld_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_cld_padecoeff(): can't open file " // trim(cld_coeff_file))

    ! Read Pade coefficient dimensions
    nband        = get_dim_size(ncid,'nband')
    nrghice      = get_dim_size(ncid,'nrghice')
    nsizereg     = get_dim_size(ncid,'nsizereg')
    ncoeff_ext   = get_dim_size(ncid,'ncoeff_ext')
    ncoeff_ssa_g = get_dim_size(ncid,'ncoeff_ssa_g')
    nbound       = get_dim_size(ncid,'nbound')

    !
    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    ! Allocate cloud property Pade coefficient input arrays
    allocate(pade_extliq(nband, nsizereg, ncoeff_ext),   &
             pade_ssaliq(nband, nsizereg, ncoeff_ssa_g), &
             pade_asyliq(nband, nsizereg, ncoeff_ssa_g), &
             pade_extice(nband, nsizereg, ncoeff_ext,   nrghice), &
             pade_ssaice(nband, nsizereg, ncoeff_ssa_g, nrghice), &
             pade_asyice(nband, nsizereg, ncoeff_ssa_g, nrghice))

    pade_extliq  = read_field(ncid, 'pade_extliq', nband, nsizereg, ncoeff_ext)
    pade_ssaliq  = read_field(ncid, 'pade_ssaliq', nband, nsizereg, ncoeff_ssa_g)
    pade_asyliq  = read_field(ncid, 'pade_asyliq', nband, nsizereg, ncoeff_ssa_g)
    pade_extice  = read_field(ncid, 'pade_extice', nband, nsizereg, ncoeff_ext, nrghice)
    pade_ssaice  = read_field(ncid, 'pade_ssaice', nband, nsizereg, ncoeff_ssa_g, nrghice)
    pade_asyice  = read_field(ncid, 'pade_asyice', nband, nsizereg, ncoeff_ssa_g, nrghice)

    ! Allocate cloud property Pade coefficient particle size boundary input arrays
    allocate(pade_sizreg_extliq(nbound), &
             pade_sizreg_ssaliq(nbound), &
             pade_sizreg_asyliq(nbound), &
             pade_sizreg_extice(nbound), &
             pade_sizreg_ssaice(nbound), &
             pade_sizreg_asyice(nbound))

    pade_sizreg_extliq = read_field(ncid, 'pade_sizreg_extliq', nbound)
    pade_sizreg_ssaliq = read_field(ncid, 'pade_sizreg_ssaliq', nbound)
    pade_sizreg_asyliq = read_field(ncid, 'pade_sizreg_asyliq', nbound)
    pade_sizreg_extice = read_field(ncid, 'pade_sizreg_extice', nbound)
    pade_sizreg_ssaice = read_field(ncid, 'pade_sizreg_ssaice', nbound)
    pade_sizreg_asyice = read_field(ncid, 'pade_sizreg_asyice', nbound)

    ncid = nf90_close(ncid)

    call stop_on_err(cloud_spec%load(band_lims_wvn, &
                                     pade_extliq, pade_ssaliq, pade_asyliq, &
                                     pade_extice, pade_ssaice, pade_asyice, &
                                     pade_sizreg_extliq, pade_sizreg_ssaliq, pade_sizreg_asyliq, &
                                     pade_sizreg_extice, pade_sizreg_ssaice, pade_sizreg_asyice))
  end subroutine load_cld_padecoeff

  ! -----------------------------------------------------------------------------------
    subroutine stop_on_err(msg)
      !
      ! Print error message and stop
      !
      use iso_fortran_env, only : error_unit
      character(len=*), intent(in) :: msg
      if(len_trim(msg) > 0) then
        write (error_unit,*) trim(msg)
        error stop 1
      end if
    end subroutine

end module
