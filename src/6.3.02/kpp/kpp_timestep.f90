!##############################################################################
Subroutine mckpp_ocean_mixing ()


use mem_leaf
use mem_turb
use mem_radiate
use mem_micro
use mem_basic
use micphys, only:idriz,irain,ipris,isnow,iaggr,igraup,ihail
use io_params, only:ssttime1,ssttime2,iupdsst
use mem_kpp, only:kpp_3d_fields,kpp_1d_fields,kpp_const_fields
use node_mod, only: ia,iz,ja,jz
use mem_grid, only:ngrid,runtype,initial,print_msg,nzg,time,dtlt
use kpp_parameters, only:frqkpp,hmixdepth,ikpp

implicit none

integer :: i,j
real :: flx1,flx2,flx3,flx4,flx5,flx6,kpp_total_pcp,timefac_sst

!Set parameter to zero to help determine the max domain HMIX
hmixdepth(ngrid) = 0.

!Check and make sure KPP timestep >= RAMS model grid timestep
if (kpp_const_fields%dto < dtlt) THEN
 if(print_msg)print*,'KPP:ocean model timestep must be >= DTLT from RAMS'
 if(print_msg)print*,'KPP-dto,frqkpp,dtlt',kpp_const_fields%dto,frqkpp,dtlt
 STOP
ENDIF

! Time interpolation factor for updating SST
if (iupdsst == 0) then
   timefac_sst = 0.
else
   timefac_sst = (time-ssttime1(ngrid)) / (ssttime2(ngrid)-ssttime1(ngrid))
endif

!Update SST from RAMS obs field for SST0 relaxation of prognosed SST.
!Note that seatp and seatf are in units of Kelvin, but SST0 needs celcius.
IF (kpp_const_fields%relax_sst .GT. 1.e-10)then
 do j = ja,jz
   do i = ia,iz
     kpp_3d_fields(ngrid)%SST0(i,j) = (leaf_g(ngrid)%seatp(i,j) &
        + (leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j)) &
        * timefac_sst - 273.15)
   enddo
 enddo
ENDIF

!Get parent model fluxes.
!Store surface fluxes in temporary variables to allow for
!coupling timesteps longer than the model timestep. That way we compute
!mean surface fluxes over the non-coupling timesteps, and use the mean
!fluxes to apply to KPP rather than instantaneous fluxes.
do j = ja,jz
 do i = ia,iz

  !Zonal wind stress (N/m2) (direction/sign dependent)
  flx1 = turb_g(ngrid)%sflux_u(i,j)

  !Meridional wind stress (N/m2) (direction/sign dependent)
  flx2 = turb_g(ngrid)%sflux_v(i,j)

  !net sfc shortwave radiation (swdn-swup)(W/m^2)
  flx3 = radiate_g(ngrid)%swdn(1,i,j) &         !shortwave down at sfc
       - radiate_g(ngrid)%swup(1,i,j)           !shortwave up at sfc

  !non-shortwave rad (lwdn-lwup-sens-latent) (W/m^2)
  flx4 = radiate_g(ngrid)%lwdn(1,i,j) &         !longwave down at sfc
       - radiate_g(ngrid)%lwup(1,i,j) &         !longwave up at sfc
       - (turb_g(ngrid)%sflux_t(i,j)*1004.) &   !sensible heat flx at sfc
       - (turb_g(ngrid)%sflux_r(i,j)*2.5e6)     !latent heat flx at sfc

  !melting of seaice(zero, not used)
  flx5 = 0.

  !net freshwater (precip - evaporation)(mm/sec)
  kpp_total_pcp = 0.
  if(idriz>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprd(i,j)
  if(irain>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprr(i,j)
  if(ipris>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprp(i,j)
  if(isnow>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprs(i,j)
  if(iaggr>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcpra(i,j)
  if(igraup>=1) kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprg(i,j)
  if(ihail>=1)  kpp_total_pcp=kpp_total_pcp+micro_g(ngrid)%pcprh(i,j)
  flx6 = kpp_total_pcp - turb_g(ngrid)%sflux_r(i,j)

  !Use instantaneous values for the first timestep
  IF (time .le. 0.001) THEN
   kpp_3d_fields(ngrid)%flx_ust(i,j) = flx1
   kpp_3d_fields(ngrid)%flx_vst(i,j) = flx2
   kpp_3d_fields(ngrid)%flx_nsw(i,j) = flx3
   kpp_3d_fields(ngrid)%flx_nlw(i,j) = flx4
   kpp_3d_fields(ngrid)%flx_ice(i,j) = flx5
   kpp_3d_fields(ngrid)%flx_pcp(i,j) = flx6
  ENDIF

  !Start storing mean SFLUX after initialization
  IF (time .gt. 0.001) THEN
   !Zonal wind stress
   kpp_3d_fields(ngrid)%flx_ust(i,j) = flx1 * dtlt/FRQKPP + &
                kpp_3d_fields(ngrid)%flx_ust(i,j)
   !Meridional wind stress
   kpp_3d_fields(ngrid)%flx_vst(i,j) = flx2 * dtlt/FRQKPP + &
                kpp_3d_fields(ngrid)%flx_vst(i,j)
   !Net downward shortwave flux
   kpp_3d_fields(ngrid)%flx_nsw(i,j) = flx3 * dtlt/FRQKPP + &
                kpp_3d_fields(ngrid)%flx_nsw(i,j)
   !Net downward longwave flux
   kpp_3d_fields(ngrid)%flx_nlw(i,j) = flx4 * dtlt/FRQKPP + &
                kpp_3d_fields(ngrid)%flx_nlw(i,j)
   !Melting of sea ice = 0.0
   kpp_3d_fields(ngrid)%flx_ice(i,j) = flx5
   !Precipitation minus evaporation
   kpp_3d_fields(ngrid)%flx_pcp(i,j) = flx6 * dtlt/FRQKPP + &
                kpp_3d_fields(ngrid)%flx_pcp(i,j)
  ENDIF

  !Run KPP if a coupling timestep and not initialization           
  IF (mod(time+0.001,frqkpp) .lt. dtlt) THEN
   !Proceed if RAMS LEAF3 has determine this to be an ocean point.
   !We allow coastal grid cells with partial ocean.
   IF (nint(leaf_g(ngrid)%leaf_class(i,j,1)) == 0. .and. &
            leaf_g(ngrid)%patch_area(i,j,1) > 0.50) THEN !Saleeby(2018)
    !Copy 3D fields stored arrays to temporary 1D column for physics
    CALL mckpp_fields_3dto1d (kpp_3d_fields(ngrid),kpp_1d_fields,i,j)
    !Run physics ocean timestep
    CALL mckpp_physics_ocnstep (kpp_1d_fields,kpp_const_fields,i,j)
    !Enforced overrides
    CALL mckpp_overrides_profile (kpp_1d_fields,kpp_const_fields,i,j)
    !Force bottom layer temperature to observations "bottomt"
    CALL mckpp_overrides_bottomtemp (kpp_1d_fields,kpp_const_fields,i,j)
    !Copy 1D updated fields from physics to 3D stored arrays
    CALL mckpp_fields_1dto3d (kpp_1d_fields,kpp_3d_fields(ngrid),i,j)

    !Output SST back to RAMS (X_t is in celcius already)
    !For now, setting top level ocean model temperature to SST
    leaf_g(ngrid)%soil_energy(nzg,i,j,1) = &
      334000. + 4186. * kpp_3d_fields(ngrid)%X_t(1,i,j)
   ENDIF
   !Reset surface fluxes to zero at all grid cells, not just ocean
   kpp_3d_fields(ngrid)%flx_ust(i,j)=0.
   kpp_3d_fields(ngrid)%flx_vst(i,j)=0.
   kpp_3d_fields(ngrid)%flx_nsw(i,j)=0.
   kpp_3d_fields(ngrid)%flx_nlw(i,j)=0.
   kpp_3d_fields(ngrid)%flx_ice(i,j)=0.
   kpp_3d_fields(ngrid)%flx_pcp(i,j)=0.
  ENDIF 

  !Set sub-domain maximum to help get total domain max HMIX
  !for printing to standard output to screen.
  if(kpp_3d_fields(ngrid)%hmix(i,j) .gt. hmixdepth(ngrid)) &
      hmixdepth(ngrid) = kpp_3d_fields(ngrid)%hmix(i,j)

 enddo
enddo

return
END SUBROUTINE mckpp_ocean_mixing

