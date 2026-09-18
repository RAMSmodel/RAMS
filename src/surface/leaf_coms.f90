!##############################################################################
Module leaf_coms

use grid_dims

implicit none

integer ::       &
    niter_leaf   & ! number of leaf timesteps in model long timestep
   ,niter_can    & ! number of canopy timesteps in leaf timestep
   ,ifreeslip    & ! flag for freeslip surface BC (set sflux to zero)
   ,icharnock      ! flag to use Charnock ocean roughness parameterization

real ::          &
    dtll         & ! leaf timestep
   ,dtll_factor  & ! leaf timestep factor (leaf timestep / model timestep)
   ,dtlc         & ! canopy timestep
   ,dtlc_factor  & ! canopy timestep factor (canopy timestep / leaf timestep)
   ,hcapcan      & ! canopy heat capacity
   ,wcapcan      & ! canopy vapor capacity
   ,hcapveg      & ! vegetation heat capacity
   ,dtllohcc     & ! leaf timestep / canopy heat capacity
   ,dtllowcc     & ! leaf timestep / canopy vapor capacity
   ,dtlcohcc     & ! canopy timestep / canopy heat capacity
   ,dtlcowcc     & ! canopy timestep / canopy vapor capacity
   ,dtlcohcv     & ! caonpy timestep / vegetation heat capacity
   
   ,ups          & ! U velocity at top of surface layer [up(2,i,j)]    
   ,vps          & ! V velocity at top of surface layer [vp(2,i,j)]
   ,ths          & ! potential temperature at top of surface layer [theta(2,i,j)]
   ,rvs          & ! vapor mixing ratio at top of surface layer [rv(2,i,j)]
   ,zts          & ! height at top of surface layer [zt(2)*rtgt(i,j)]
   ,pis          & ! Exner func at surface
   ,dens         & ! density at surface
   ,prss         & ! pressure at surface
   ,vels         & ! wind speed at top of surface layer
   ,vels_pat     & ! vels with patch-dependent ubmin imposed as minimum
   ,gzotheta     & ! (g*z/theta) at top of surface layer [for Richardson number]
   ,pcpgl        & ! precip mass from cuparm and/or micphys in leaf timestep
   ,qpcpgl       & ! precip energy from cuparm and/or micphys in leaf timestep
   ,dpcpgl       & ! precip depth from cuparm and/or micphys in leaf timestep
   ,pcpgc        & ! precip mass from cuparm and/or micphys in canopy timestep
   ,qpcpgc       & ! precip energy from cuparm and/or micphys in canopy timestep
   ,z0fac_water  & ! (.016 / g) factor of ustar^2 for z0 over water
   
   ,snowfac      & ! fraction of vegetation height covered by sfcwater
   ,vf           & ! product of veg_fracarea and (1-snowfac)
   ,thetacan     & ! canopy air potential temperature
   ,transp       & ! transpiration flux [kg/m2/s]
   ,snowrough    & ! snowcover roughness height
   ,timefac_sst  & ! time interpolation factor for SST
   
   ,rb           & ! vegetation aerodynamic resistance
   ,rd           & ! canopy to ground aerodynamic resistance
   ,rdi          & ! canopy to ground aerodynamic conductance
   
   ,rshort_g     & ! net SW radiation absorbed by grnd
   ,rshort_v     & ! net SW radiation absorbed by veg
   ,rshort_a     & ! net SW radiation reflected to atm by veg plus grnd
   
   ,rlonga_v     & ! net atm LW radiation absorbed by veg
   ,rlonga_gs    & ! net atm LW radiation absorbed by grnd OR snow
   ,rlongv_gs    & ! net veg LW radiation absorbed by grnd OR snow
   ,rlongv_a     & ! net veg LW radiation to atm
   ,rlonggs_v    & ! net grnd OR snow LW radiation absorbed by veg
   ,rlonggs_a    & ! net grnd OR snow LW radiation to atm
   ,rlonga_a     & ! net atm LW radiation reflected to atm by veg plus grnd OR snow

   ,hflxgc       & ! sensible heat from ground to canopy (J/m2)
   ,wflxgc       & ! water vapor from ground to canopy (kg/m2)
   ,hflxvc       & ! sensible heat from vegetation to canopy (J/m2)
   ,wflxvc       & ! water vapor from vegetation to canopy (kg/m2)
   
   ,wshed        & ! water shed from vegetation to ground (kg/m2)
   ,qwshed       & ! energy from shed water (J/m2)
   ,dewgnd         ! dew formation on ground (kg/m2)

real, dimension(nzgmax+nzsmax+1) ::  &
    dslz         & ! soil layer thickness at T point (nzg)
   ,dslzi        & ! (1. / soil layer thickness at T point) (nzg)
   ,dslzidt      & ! (dtll / soil layer thickness at T point) (nzg)
   ,slzt         & ! soil depth at T point (nzg)
   ,dslzt        & ! soil layer thickness at M point (nzg)
   ,dslzti       & ! (1. / soil layer thickness at M point) (nzg)
   ,dslztidt     & ! (dtll / soil layer thickness at M point) (nzg)

   ,rshort_s     & ! net SW radiation absorbed by snow layers (nzs)
   ,tempk        & ! diagnosed temp (K) of soil and sfcwater (nzg+nzs)
   ,fracliq      & ! diagnosed liquid fraction of soil_water and sfcwater_mass (nzg+nzs)
   
   ,hfluxgsc      & ! sensible heat flux between soil, sfcwater, canopy (nzg+nzs+1)
   ,psiplusz      & ! soil water potential plus geopotential [m] (nzg)
   ,half_soilair  & ! half of available airspace in soil [m] (nzg)
   ,rfactor       & ! soil, sfcwater thermal resistance (nzg+nzs) 
   ,wflux         & ! soil water flux [m] (nzg+1)
   ,soil_liq      & ! soil liquid water content [m] (nzg+1)
   ,qwflux          ! soil energy flux from water flux [J/m2] (nzg)

integer, parameter :: nstyp=12,nvtyp=20

real, dimension(nstyp)        :: slcpd,slbs,sfldcap  &
                                ,slmsts,slpots  &
                                ,soilcp,emisg
real, dimension(0:nvtyp)      :: albv_green,albv_brown,emisv,sr_max,tai_max  &
                                ,sai,veg_clump,veg_frac,veg_ht,glai_max  &
                                ,dead_frac,rcmin
integer, dimension(0:nvtyp)   :: kroot
real, dimension(nzgmax,nstyp) :: slcons1

END MODULE leaf_coms
