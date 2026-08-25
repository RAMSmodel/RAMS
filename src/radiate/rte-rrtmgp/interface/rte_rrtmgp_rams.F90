module rte_rrtmgp_rams

 ! use mo_optical_props,      only: ty_optical_props, &
 !                                  ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_aerosol_optics,     only: ty_aerosol_optics
  use mo_gas_concentrations, only: ty_gas_concs
!  use mo_fluxes,             only: ty_fluxes_broadband

  integer, parameter :: ngas = 6
  !These 6 gases are "required". I've set ch4 and n2o to 0 and o2 to .209
  !You can turn off the strict requirement for these by modifying key_species in the rrtmgp data files in etc
  character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ','ch4','n2o','o2 ']
  type(ty_gas_optics_rrtmgp) :: k_dist_lw, k_dist_sw
  type(ty_aerosol_optics)    :: aerosol_optics_lw, aerosol_optics_sw
  type(ty_cloud_optics)      :: cloud_optics_lw, cloud_optics_sw
  type(ty_gas_concs)         :: gas_concs, gas_concs_garand, gas_concs_1col
!  class(ty_optical_props_arry), &
!                 allocatable :: atmos_sw, atmos_lw, clouds_sw, clouds_lw
end module
