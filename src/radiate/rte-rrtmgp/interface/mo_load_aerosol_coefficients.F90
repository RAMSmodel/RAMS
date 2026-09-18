module mo_load_aerosol_coefficients
  use mo_rte_kind,      only: wp
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_aerosol_optics,  only: ty_aerosol_optics
  use mo_simple_netcdf, only: read_field, read_string, var_exists, get_dim_size, &
                              write_field, create_dim, create_var
  use netcdf

  implicit none
  private
  public :: load_aero_lutcoeff
  ! ----------------------------------------------------------------------------------

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! read cloud optical property LUT coefficients from NetCDF file
  !
  subroutine load_aero_lutcoeff(aero_spec, aero_coeff_file)
    class(ty_aerosol_optics),     intent(inout) :: aero_spec
    character(len=*),           intent(in   ) :: aero_coeff_file
    ! -----------------
    ! Local variables
    integer :: ncid, nband, ntype, ndust, ncarb, nbins, nradbins, nrh

    real(wp), dimension(:,:), allocatable                :: band_lims_wvn, radparam
    real(wp), dimension(:), allocatable :: rma, radlim
    ! LUT coefficients
    real(wp), dimension(:,:,:,:), allocatable :: lut_qext   ! 
    real(wp), dimension(:,:,:), allocatable :: lut_qext_carb   ! 
    real(wp), dimension(:,:,:), allocatable :: lut_qext_dust   ! 
    real(wp), dimension(:,:,:,:), allocatable :: lut_qscat   ! 
    real(wp), dimension(:,:,:), allocatable :: lut_qscat_carb   ! 
    real(wp), dimension(:,:,:), allocatable :: lut_qscat_dust   ! 
    real(wp), dimension(:,:,:,:), allocatable :: lut_gasym   !
    real(wp), dimension(:,:,:), allocatable :: lut_gasym_carb   !
    real(wp), dimension(:,:,:), allocatable :: lut_gasym_dust   !
    real(wp), dimension(:,:,:), allocatable :: lut_growth   !

    ! -----------------
    ! Open aerosol optical property coefficient file
    if(nf90_open(trim(aero_coeff_file), NF90_NOWRITE, ncid) /= NF90_NOERR) &
       call stop_on_err("load_aero_lutcoeff(): can't open file " // trim(aero_coeff_file))

    ! Read LUT coefficient dimensions
    nband     = get_dim_size(ncid,'nbands')
    nrh = get_dim_size(ncid,'nrh')
    nradbins = get_dim_size(ncid,'nradbins')
    nbins = get_dim_size(ncid,'nbins')
    ntype = get_dim_size(ncid,'ntype')
    ncarb = get_dim_size(ncid,'ncarb')
    ndust = get_dim_size(ncid,'ndust')

    allocate(band_lims_wvn(2, nband))
    band_lims_wvn = read_field(ncid, 'bnd_limits_wavenumber', 2, nband)

    allocate(rma(nradbins))
    rma = read_field(ncid,'rma',nradbins)

    allocate(radlim(nbins))
    radlim = read_field(ncid,'radlim',nbins)

    allocate(radparam(nbins,nradbins))
    radparam = read_field(ncid,'radparam',nbins, nradbins)

    ! Allocate aerosol property lookup table input arrays
    allocate(lut_qext(ntype, nband, nradbins,nrh), &
             lut_qext_carb(ncarb, nband, nradbins), &
             lut_qext_dust(ndust, nband, nradbins), &
             lut_qscat(ntype, nband, nradbins,nrh), &
             lut_qscat_carb(ncarb, nband, nradbins), &
             lut_qscat_dust(ndust, nband,nradbins), &
             lut_gasym(ntype, nband, nradbins, nrh), &
             lut_gasym_carb(ncarb, nband, nradbins), &
             lut_gasym_dust(ndust, nband,nradbins), &
             lut_growth(ntype, nradbins, nrh))

    ! Read LUT coefficients
    lut_qext = read_field(ncid, 'qext',  ntype, nband, nradbins, nrh)
    lut_qext_carb = read_field(ncid, 'qext_carb',  ncarb, nband, nradbins)
    lut_qext_dust = read_field(ncid, 'qext_dust',  ndust, nband, nradbins)
    lut_qscat = read_field(ncid, 'qscat',  ntype, nband, nradbins, nrh)
    lut_qscat_carb = read_field(ncid, 'qscat_carb',  ncarb, nband, nradbins)
    lut_qscat_dust = read_field(ncid, 'qscat_dust',  ndust, nband, nradbins)
    lut_gasym = read_field(ncid, 'gasym',  ntype, nband, nradbins, nrh)
    lut_gasym_carb = read_field(ncid, 'gasym_carb',  ncarb, nband, nradbins)
    lut_gasym_dust = read_field(ncid, 'gasym_dust',  ndust, nband, nradbins)
    lut_growth = read_field(ncid, 'growth',ntype,nradbins,nrh)

    ncid = nf90_close(ncid)
    call stop_on_err(aero_spec%load(band_lims_wvn, rma, radlim, radparam,  &
                                    lut_qext, lut_qext_carb, lut_qext_dust,    &
                                    lut_qscat, lut_qscat_carb, lut_qscat_dust, &
                                    lut_gasym, lut_gasym_carb, lut_gasym_dust, &
                                    lut_growth))
  end subroutine load_aero_lutcoeff

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
