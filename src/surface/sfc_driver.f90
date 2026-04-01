!##############################################################################
Subroutine sfc_driver (mzp,mxp,myp,ia,iz,ja,jz)

use mem_all

implicit none

integer :: mzp,mxp,myp,ia,iz,ja,jz

integer :: ng

if (nstbot == 0) return

ng=ngrid

CALL leaf3_sib (mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz      &
   ,leaf_g (ng), basic_g (ng), turb_g (ng), radiate_g(ng)   &
   ,grid_g (ng), cuparm_g(ng), micro_g(ng), sib_g(ng))

return
END SUBROUTINE sfc_driver

!##############################################################################
Subroutine leaf3_sib (m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz   &
   ,leaf,basic,turb,radiate,grid,cuparm,micro,sib)

use mem_all
use leaf_coms
use rconstants

implicit none

integer :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz

type (leaf_vars)     leaf
type (basic_vars)    basic
type (turb_vars)     turb
type (radiate_vars)  radiate
type (grid_vars)     grid
type (cuparm_vars)   cuparm
type (micro_vars)    micro
type (sib_vars)      sib

real, allocatable, dimension(:,:) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2
integer :: i,j,ip,iter_leaf,runstars
integer :: jday,oyr,omn,ody,otm
integer, external :: julday
real,external :: rslif

! SIB2 - temporary diagnostic/output SIB2 variables
real :: cupr,lspr
integer :: ksn,nsoil,k

! To convert CO2 back and forth between ppm and mixing ratio
! (mass fraction) units.
real :: f
data f /1.51724e-6/  ! 1.51724e-6= 44/29/1.e6

allocate(ths2(m2,m3),rvs2(m2,m3),pis2(m2,m3),dens2(m2,m3) &
        ,ups2(m2,m3),vps2(m2,m3),zts2(m2,m3))

! Time interpolation factor for updating SST
if (iupdsst == 0) then
   timefac_sst = 0.
else
   timefac_sst = (time-ssttime1(ngrid)) / (ssttime2(ngrid)-ssttime1(ngrid))
endif

! Define LEAF3 and canopy time-split timesteps here.  This ensures that LEAF3
! will not use a timestep longer than about 40 seconds, and canopy will not
! use a timestep longer than about 15 seconds.  This allows values of 
! hcapcan = 2.e4, wcapcan = 2.e1, and hcapveg = 3.e4 as are now defined below.

!Saleeby(2010):Moved the computation of niter_can to after computation of dtll
niter_leaf  = max(1,nint(dtlt/40.+.4))
dtll_factor = 1. / float(niter_leaf)
dtll        = dtlt * dtll_factor

niter_can   = max(1,nint(dtll/15.+.4))
dtlc_factor = 1. / float(niter_can)
dtlc        = dtll * dtlc_factor

! Define SIB timestep as the same as dtlong (dtlt)
if(isfcl==2)then
 niter_leaf = 1
 niter_can  = 1
 dtll_factor = 1.
 dtlc_factor = 1.
 dtll=dtlt
 dtlc=dtlt
endif

hcapcan = 2.0e4
wcapcan = 2.0e1
hcapveg = 3.e4

dtllohcc = dtll / hcapcan
dtllowcc = dtll / wcapcan
dtlcohcc = dtlc / hcapcan
dtlcowcc = dtlc / wcapcan
dtlcohcv = dtlc / hcapveg

z0fac_water = .016 / g
snowrough = .001

! Copy surface atmospheric variables into 2d arrays for input to leaf
 CALL sfc_fields (m1,m2,m3,ia,iz,ja,jz,jdim                    &
      ,basic%theta(1,1,1) ,basic%rv (1,1,1) ,basic%up(1,1,1)  &
      ,basic%vp   (1,1,1) ,basic%dn0(1,1,1) ,basic%pp(1,1,1)  &
      ,basic%pi0  (1,1,1) ,grid%rtgt(1,1)   ,zt               &
      ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2                    )

do j = ja,jz
 do i = ia,iz

! Copy surface variables to single-column values
  ups = ups2(i,j)
  vps = vps2(i,j)
  ths = ths2(i,j)
  rvs = rvs2(i,j)
  zts = zts2(i,j)
  pis = pis2(i,j)
  dens = dens2(i,j)
  prss = pis ** cpor * p00
  vels = sqrt(ups ** 2 + vps ** 2)
  gzotheta = g * zts / ths

! Update water internal energy from time-dependent SST.
! Do this only if NOT running KPP ocean mixed layer model.
! If KPP is running, then KPP will update the SST / water-internal-energy.
! Note however that KPP may not be run every timestep and thus may only
! update SST/water-internal-energy at its timesteps.
  if(IKPP==0)then
   leaf%soil_energy(mzg,i,j,1) = 334000.  &
     + 4186. * (leaf%seatp(i,j) + (leaf%seatf(i,j) - leaf%seatp(i,j)) &
     * timefac_sst - 273.15)
  endif

! Fill surface precipitation arrays for input
  if(isfcl<=1) then     !for LEAF3
   CALL sfc_pcp (nnqparm(ngrid),level,i,j,cuparm,micro)
  elseif(isfcl==2) then !for SiB
   CALL sfc_pcp_sib (nnqparm(ngrid),level,i,j,cuparm,micro,cupr,lspr)
  endif

!Zero out albedo,upward surface longwave,and momentum,heat,and moisture
!flux arrays before summing over patches
  if (ilwrtyp > 0 .or. iswrtyp > 0) then
    radiate%albedt(i,j) = 0.
    radiate%rlongup(i,j) = 0.
  endif

  turb%sflux_u(i,j) = 0.
  turb%sflux_v(i,j) = 0.
  turb%sflux_w(i,j) = 0.
  turb%sflux_t(i,j) = 0.
  turb%sflux_r(i,j) = 0.

! For no soil model (patch 2) fill "canopy" temperature and moisture
  if (isfcl == 0) then
     leaf%can_temp  (i,j,2) = (ths - dthcon) * pis
     leaf%can_rvap  (i,j,2) = rvs - drtcon
     leaf%patch_area(i,j,1) = 1. - pctlcon
     leaf%patch_area(i,j,2) = pctlcon
  endif

! Begin patch loop
  do ip = 1,np

! Zero out SIB albedo,upward surface longwave,and momentum,heat,and moisture
! flux arrays before summing over patches.
! If running SIB land surface, use SiB surface albedo and surface
! upward longwave after first timestep.
   if(isfcl==2)then
    if (ilwrtyp > 0 .or. iswrtyp > 0) then
      sib%sfcswa(i,j,ip) = 0.
      sib%uplwrf(i,j,ip) = 0.
    endif
   endif

! Update time-dependent vegetation LAI and fractional coverage
   if (ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then
    if (isfcl>=1) CALL ndvi (ngrid                           &
        ,leaf%leaf_class  (i,j,ip)                           &
        ,leaf%veg_ndvip   (i,j,ip) ,leaf%veg_ndvic (i,j,ip)  &
        ,leaf%veg_ndvif   (i,j,ip))
    if (isfcl==1) CALL veg (ngrid                            &
        ,leaf%leaf_class  (i,j,ip)                           &
        ,leaf%veg_fracarea(i,j,ip) ,leaf%veg_lai   (i,j,ip)  &
        ,leaf%veg_tai     (i,j,ip) ,leaf%veg_rough (i,j,ip)  &
        ,leaf%veg_height  (i,j,ip) ,leaf%veg_albedo(i,j,ip)  &
        ,leaf%veg_ndvic   (i,j,ip))
   endif

! Begin leaf small timestep here.
   do iter_leaf = 1,niter_leaf

! Calculate radiative fluxes between atmosphere, vegetation, and ground/snow
! based on already-computed downward shortwave and longwave fluxes from
! the atmosphere.  Fill tempk array with soil and snow temperature (C) and
! fracliq array with liquid fraction of water content in soil and snow.
! Other snowcover properties are also computed here.

   if (iswrtyp > 0 .or. ilwrtyp > 0) then

      if (ip == 1 .or. (isfcl<=1 .and. leaf%patch_area(i,j,ip) >= .009)) then

         CALL sfcrad (mzg,mzs,ip                                           &
         ,leaf%soil_energy    (1,i,j,ip) ,leaf%soil_water      (1,i,j,ip)  &
         ,leaf%soil_text      (1,i,j,ip) ,leaf%sfcwater_energy (1,i,j,ip)  &
         ,leaf%sfcwater_depth (1,i,j,ip) ,leaf%patch_area      (i,j,ip)    &
         ,leaf%can_temp       (i,j,ip)   ,leaf%veg_temp        (i,j,ip)    &
         ,leaf%leaf_class     (i,j,ip)   ,leaf%veg_height      (i,j,ip)    &
         ,leaf%veg_fracarea   (i,j,ip)   ,leaf%veg_albedo      (i,j,ip)    &
         ,leaf%sfcwater_nlev  (i,j,ip)                                     &
         ,radiate%rshort      (i,j)      ,radiate%rlong        (i,j)       &
         ,radiate%albedt      (i,j)      ,radiate%rlongup      (i,j)       &
         ,radiate%cosz        (i,j)    )

      endif

   else

      if(ip==1) CALL qtk (leaf%soil_energy(mzg,i,j,ip),tempk(mzg),fracliq(mzg))

   endif

! For water surface (patch 1), compute surface saturation mixing ratio
! and roughness length based on previous ustar.
! For soil patches, compute roughness length based on vegetation and snow.

   if (ip == 1) then

      leaf%ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))   
      leaf%patch_rough(i,j,ip)  &
                  = max(z0fac_water * leaf%ustar(i,j,ip) ** 2,.0001)

   elseif (isfcl>=1) then

      if(isfcl==2) snowfac = &
       min(.99, leaf%sfcwater_depth(1,i,j,ip) / max(.001,leaf%veg_height(i,j,ip)))

      if (leaf%patch_area(i,j,ip) >= .009) then
        leaf%patch_rough(i,j,ip)   &
                     = max(grid%topzo(i,j),leaf%soil_rough(i,j,ip)  &
                        ,leaf%veg_rough(i,j,ip)) * (1. - snowfac)   &
                         + snowrough * snowfac
      endif

   endif

! Calculate turbulent fluxes between atmosphere and canopy (or "canopy").

   if (leaf%patch_area(i,j,ip) < .009 .and. ip >= 2 .and. &
      (isfcl == 1 .or. isfcl == 2)) then
        thetacan = ths
   else
        thetacan = leaf%can_temp(i,j,ip) / pis
   endif

! Set minimum wind speed for fluxes

   if(thetacan < ths) then   ! STABLE CASE
      ubmin = .1
   else                      ! UNSTABLE CASE
      ubmin = 1.
   endif

   !For RCE simulations, set a default ubmin
   if(irce == 1) ubmin = rce_ubmn

   !For KPP simulations, override ubmin if namelist value is larger
   if(ikpp > 0) ubmin = max(ubmin,ubmn_kpp)

   !Normal velocity used in computing stars
   vels_pat = max(vels,ubmin)

! Only run stars here if using LEAF3 or it is the first model timestep 
! or if its a water patch.

   runstars=0
   if    (isfcl<=1) then
      runstars=1
   elseif(isfcl==2) then
      if(ip==1 .or. (initial<=2 .and. time < 0.001) .or. &
         leaf%patch_area(i,j,ip) < .009) runstars=1
   endif
   if(runstars==1) &
    CALL stars (leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)  &
        ,leaf%rstar(i,j,ip),ths,rvs,thetacan,leaf%can_rvap(i,j,ip)  &
        ,zts,leaf%patch_rough(i,j,ip),vels_pat,dtllohcc,dens,dtlt)

! For water patches, update temperature and moisture of "canopy" from
! divergence of fluxes with water surface and atmosphere.  rdi = ustar/5
! is the viscous sublayer conductivity from Garratt (1992).

   if (ip == 1) then

      rdi = .2 * leaf%ustar(i,j,1)

      leaf%can_temp(i,j,1) = leaf%can_temp(i,j,1)        &
          + dtllohcc * dens * cp                          &
          * ((tempk(mzg) - leaf%can_temp(i,j,1)) * rdi    &
          + leaf%ustar(i,j,1) * leaf%tstar(i,j,1) * pis)

      leaf%can_rvap(i,j,1) = leaf%can_rvap(i,j,1) + dtllowcc * dens &
          * ((leaf%ground_rsat(i,j,1) - leaf%can_rvap(i,j,1)) * rdi  &
          + leaf%ustar(i,j,1) * leaf%rstar(i,j,1))

   endif

   !***************************************************************************
   ! SIB: For soil model patches,update temperature and moisture of soil,
   ! vegetation, and canopy using SiB model
   !***************************************************************************
   if (isfcl == 2 .and. ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then
     ! Diagnose soil temperature and liquid fraction
     ! as done in routine sfcrad
     do k = 1,mzg
        nsoil = nint(leaf%soil_text(k,i,j,ip))
        CALL qwtk (leaf%soil_energy(k,i,j,ip),  &
             leaf%soil_water(k,i,j,ip)*1.e3,   &
             slcpd(nsoil),tempk(k),fracliq(k))
     enddo
     ! Diagnose snow temperature as done in routine sfcrad
     ! SiB current allows only 1 surface water/snow layer
     ksn = nint(leaf%sfcwater_nlev(i,j,ip))
     do k = 1,ksn
        CALL qtk (leaf%sfcwater_energy(k,i,j,ip),  &
             tempk(k+mzg), fracliq(k+mzg))
     enddo
     !if there is not snow
     if(ksn == 0) then
        do k = 1,mzs
           tempk(k+mzg) = tempk(mzg)
        enddo
     endif
     !Compute the Julian calendar day
     CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time &
                                ,'s',oyr,omn,ody,otm)
     jday=julday(omn,ody,oyr)
     !Call to SIB main routine here.
     CALL sib_2pt5 (          &
     !THIS SECTION: INPUT GRID AND ATMOSPHERIC CONDITIONS FROM RAMS TO SIB
     i                        &!I grid point
     ,j                       &!J grid point
     ,mzg                     &!number soil layers
     ,mzs                     &!number snow layers
     ,dtlt                    &!timestep (sec)
     ,jday                    &!time in julian date (Day Of Year)
     ,ubmin                   &!minimum wind speed for surface fluxes
     ,grid%glat(i,j)          &!latitude
     ,ths                     &!lowest layer theta (potential temperature K)
     ,rvs                     &!lowest layer rv (water vapor mixing ratio g/kg)
     ,0.01*prss               &!lowest layer air pressure (hPa=mb) RAMS uses (Pa)
     ,dens                    &!lowest layer base state density (kg/m3)
     ,ths*pis                 &!lowest layer temperature (K)
     ,vels                    &!lowest layer wind speed (m/s)
     ,zts                     &!lowest model level height t-grid (m)
     ,cupr                    &!cuparm precip rate (conprr)(mm/sec)
     ,lspr/dtlt               &!microphysics precip rate (pcpg/dtlt)(mm/sec)
     ,radiate%rshort(i,j)     &!surface shortwave down (W/m2)
     ,radiate%rlong(i,j)      &!surface longwave down (W/m2)
     ,radiate%cosz(i,j)       &!cosine solar zenith angle
     !THIS SECTION: SURFACE FIELD INFO FROM RAMS TO SIB
     ,leaf%patch_area(i,j,ip)     &!BCs: -fractional patch area
     ,leaf%leaf_class(i,j,ip)     &!BCs: -leaf vegetation class
     ,leaf%veg_ndvic(i,j,ip)      &!BCs: -current month ndvi
     ,leaf%veg_ndvip(i,j,ip)      &!BCs: -previous/past month ndvi
     ,leaf%soil_text(mzg,i,j,ip)  &!BCs: -soil textural class
     !THIS SECTION: PROGNOSTIC VALUES FROM SIB TO FEEDBACK TO RAMS
     ,leaf%stom_resist(i,j,ip)  &!prog: -CANOPY/STOMATAL RESISTANCE (S M-1)
     ,leaf%can_temp(i,j,ip)     &!prog: -canopy temperature (K)
     ,leaf%can_rvap(i,j,ip)     &!prog: -canopy water vapor mixing ratio (kg/kg)
     ,leaf%veg_temp(i,j,ip)     &!prog: -vegetation temperature (K)
     ,tempk                     &!prog: -soil and snow temperature (K)
     ,leaf%soil_water(1,i,j,ip) &!prog: -volumetric soil moisture (m3/m3)
     ,sib%sfcswa(i,j,ip)        &!prog: -surface abledo (fraction)
     ,sib%uplwrf(i,j,ip)        &!prog: -upwelling longwave radiation (W/m2)
     ,sib%snow1(i,j,ip)         &!prog: -vegetation snow (kg/m2)
     ,sib%snow2(i,j,ip)         &!prog: -ground surface snow (kg/m2)
     ,sib%capac1(i,j,ip)        &!prog: -vegetation liquid store (kg/m^2)
     ,sib%capac2(i,j,ip)        &!prog: -grnd sfc liquid interception (kg/m^2)
     ,sib%pco2ap(i,j,ip)        &!prog: -canopy air space pCO2 (Pa)
     ,sib%co2flx(i,j,ip)        &!prog: -sfc CO2 flux (CAS to ref lev)(mol/m2/sec)
     ,leaf%ustar(i,j,ip)        &!prog: -u-star
     ,leaf%rstar(i,j,ip)        &!prog: -r-star
     ,leaf%tstar(i,j,ip)        &!prog: -t-star
     !THIS SECTION: INPUT CURRENT CO2 CONENTRATION TO SIB
     ,sib%rco2p(2,i,j)/f        &!lowest atm level CO2, in ppm (input to Sib)
     !THIS SECTION: DIAGNOSTIC QUANTITIES FROM SIB TO FILL STANDARD LEAF3 VARS
     ,leaf%sfcwater_nlev(i,j,ip)    &!diag: number of sfc water levels (0 or 1 from SiB)
     ,leaf%sfcwater_mass(1,i,j,ip)  &!diag: sfc water mass (snow+liquid) (kg/m2)
     ,leaf%sfcwater_depth(1,i,j,ip) &!diag: sfc water depth (snow+liquid) (m)
     ,leaf%veg_water(i,j,ip)        &!diag: vegetation water mass (snow+liquid) (kg/m2)
     !THIS SECTION: DIAGNOSTIC QUANTITIES UNIQUE TO SIB
     ,leaf%veg_albedo(i,j,ip)       &!diag: vegetation albedo (fraction)
     ,leaf%veg_fracarea(i,j,ip)     &!diag: vegetation fractional area
     ,leaf%veg_lai(i,j,ip)          &!diag: vegetation LAI
     ,leaf%veg_tai(i,j,ip)          &!diag: total LAI
     ,leaf%veg_height(i,j,ip)       &!diag: canopy top height (m)
     ,leaf%veg_rough(i,j,ip)        &!diag: vegetation roughness (m)
     ,sib%assimn(i,j,ip)  &!diag: net co2 assim by plants (umol/m^2/sec)
     ,sib%respg(i,j,ip)   &!diag: ground respiration flux (umol/m^2/sec)
     ,sib%rstfac1(i,j,ip) &!diag: CANOPY RESISTANCE STRESS1 :leaf sfc humidity
     ,sib%rstfac2(i,j,ip) &!diag: CANOPY RESISTANCE STRESS2 :soil moisture
     ,sib%rstfac3(i,j,ip) &!diag: CANOPY RESISTANCE STRESS3 :temperature
     ,sib%ect(i,j,ip)     &!diag: transpiration flux (W/m^2)
     ,sib%eci(i,j,ip)     &!diag: canopy interception flux (W/m^2)
     ,sib%egi(i,j,ip)     &!diag: ground interception flux (W/m^2)
     ,sib%egs(i,j,ip)     &!diag: ground surface layer evap (W/m^2)
     ,sib%hc(i,j,ip)      &!diag: canopy (veg) sensible heat flux (W/m^2)
     ,sib%hg(i,j,ip)      &!diag: ground surface sensible heat flux (W/m^2)
     ,sib%ra(i,j,ip)      &!diag: CAS-RAMS resistance (sec/m)
     ,sib%rb(i,j,ip)      &!diag: leaf sfc-CAS resistance (sec/m)
     ,sib%rc(i,j,ip)      &!diag: total canopy resistance (sec/m)
     ,sib%rd(i,j,ip)      &!diag: ground-CAS resistance (sec/m)
     ,sib%roff(i,j,ip)    &!diag: runoff (surface and subsurface) (mm)
     ,sib%green(i,j,ip)   &!diag: greenness fraction (-)
     ,sib%apar(i,j,ip)    &!diag: absorbed fraction of PAR
     ,sib%ventmf(i,j,ip)  &!diag: ventilation mass flux (kg/m^2/sec)
     ,sib%pco2c(i,j,ip)   &!diag: leaf chloroplast CO2 concentration (Pa)
     ,sib%pco2i(i,j,ip)   &!diag: leaf internal CO2 concentration (Pa)
     ,sib%pco2s(i,j,ip)   &!diag: leaf surface CO2 concentration (Pa)
     ,sib%pco2m(i,j,ip)   &!diag: lowest atm level CO2 concentration (Pa)
     ,sib%ea(i,j,ip)      &!diag: CAS water vapor pressure (hPa)
     ,sib%em(i,j,ip)      &!diag: ref level vapor pressure (hPa)
     ,sib%rha(i,j,ip)     &!diag: CAS relative humidity
     ,sib%radvbc(i,j,ip)  &!diag: radiation: visible beam (W/m^2)
     ,sib%radvdc(i,j,ip)  &!diag: radiation: visible diffuse (W/m^2)
     ,sib%radnbc(i,j,ip)  &!diag: radiation: nir beam (W/m^2)
     ,sib%radndc(i,j,ip)  &!diag: radiation: nir diffuse (W/m^2)
     ,sib%psy(i,j,ip)     )!diag: psychrometric constant (hPa deg^-1)
     !calculate back soil energy 
     do k=1,mzg
        nsoil = nint(leaf%soil_text(k,i,j,ip))
        leaf%soil_energy(k,i,j,ip) = (tempk(k)- 273.15)  &
             *(slcpd(nsoil)+leaf%soil_water(k,i,j,ip)*4.186e6) &
             + leaf%soil_water(k,i,j,ip)*3.34e8
     enddo
     !calculate back sfc water energy
     !(rough inverse of function qtk in therm_lib.f90)
     ksn = nint(leaf%sfcwater_nlev(i,j,ip))
     do k = 1,ksn
       if(sib%snow2(i,j,ip)>0.0) then
          leaf%sfcwater_energy(k,i,j,ip)=(tempk(k+mzg)-273.15)*2093.
       elseif(sib%capac2(i,j,ip)>0.0) then
          leaf%sfcwater_energy(k,i,j,ip)=(tempk(k+mzg)-193.36)*4186.
       else
          leaf%sfcwater_energy(k,i,j,ip)=0.0
       endif
     enddo
     !Compute ground saturation mixing ratios using LEAF3 method but based
     !on soil and canopy prediction from SiB. Diagnostic info only.
     CALL grndvap (mzs,leaf%soil_energy(mzg,i,j,ip),leaf%soil_water(mzg,i,j,ip) &
       ,leaf%soil_text(mzg,i,j,ip),leaf%sfcwater_energy(1,i,j,ip)         &
       ,leaf%sfcwater_nlev(i,j,ip),leaf%ground_rsat(i,j,ip)                 &
       ,leaf%ground_rvap(i,j,ip),leaf%can_rvap(i,j,ip),prss,i,j)

     !Update albedt and rlongup with SiB patch values since these will only
     !contain values for water patch thus far.
     if (ilwrtyp > 0 .or. iswrtyp > 0) then
       radiate%albedt (i,j) = radiate%albedt (i,j) &
         + leaf%patch_area(i,j,ip) * sib%sfcswa(i,j,ip)
       radiate%rlongup(i,j) = radiate%rlongup(i,j) &
         + leaf%patch_area(i,j,ip) * sib%uplwrf(i,j,ip)
     endif

   !***************************************************************************
   ! LEAF3: For soil model patches, update temperature and moisture of soil,
   ! vegetation, and canopy
   !***************************************************************************
   elseif (isfcl == 1 .and. ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then
     CALL leaftw (mzg,mzs  &
            ,leaf%soil_water     (1,i,j,ip) ,leaf%soil_energy    (1,i,j,ip)  &
            ,leaf%soil_text      (1,i,j,ip) ,leaf%sfcwater_mass  (1,i,j,ip)  &
            ,leaf%sfcwater_energy(1,i,j,ip) ,leaf%sfcwater_depth (1,i,j,ip)  &
            ,leaf%ustar            (i,j,ip) ,leaf%tstar            (i,j,ip)  &
            ,leaf%rstar            (i,j,ip) ,leaf%veg_fracarea     (i,j,ip)  &
            ,leaf%veg_lai          (i,j,ip) ,leaf%veg_tai          (i,j,ip)  &
            ,leaf%veg_rough        (i,j,ip) ,leaf%veg_height       (i,j,ip)  &
            ,leaf%leaf_class       (i,j,ip)  &
            ,leaf%soil_rough       (i,j,ip) ,leaf%sfcwater_nlev    (i,j,ip)  &
            ,leaf%stom_resist      (i,j,ip) ,leaf%ground_rsat      (i,j,ip)  &
            ,leaf%ground_rvap      (i,j,ip) ,leaf%veg_water        (i,j,ip)  &
            ,leaf%veg_temp         (i,j,ip) ,leaf%can_rvap         (i,j,ip)  &
            ,leaf%can_temp         (i,j,ip) ,radiate%rshort        (i,j)     &
            ,ip,i,j )
   endif !call LEAF3 or SIB

   !Compute turbulent fluxes for each patch after SiB has computed stars
   !Patch-1 stars comes from LEAF3
   !Patch-2 stars comes from SiB
   CALL sfclmcv (leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)     &
            ,leaf%rstar(i,j,ip),vels_pat,ups,vps,gzotheta  &
            ,leaf%patch_area(i,j,ip),turb%sflux_u(i,j)     &
            ,turb%sflux_v(i,j),turb%sflux_w(i,j)           &
            ,turb%sflux_t(i,j),turb%sflux_r(i,j))

   enddo !time iteration loop
  enddo !patch loop

 enddo !i loop
enddo !j loop
 
! Normalize accumulated fluxes and albedo seen by atmosphere over model
! timestep dtlt.

do j = ja,jz
   do i = ia,iz
     !If NOT using freeslip BC
     if(ifreeslip==0)then
      turb%sflux_u(i,j) = turb%sflux_u(i,j) * dens2(i,j) * dtll_factor
      turb%sflux_v(i,j) = turb%sflux_v(i,j) * dens2(i,j) * dtll_factor
      turb%sflux_w(i,j) = turb%sflux_w(i,j) * dens2(i,j) * dtll_factor
      turb%sflux_t(i,j) = turb%sflux_t(i,j) * dens2(i,j) * dtll_factor
      turb%sflux_r(i,j) = turb%sflux_r(i,j) * dens2(i,j) * dtll_factor
     endif
     !If using freeslip surface boundary conditions, set fluxes to zero
     if(ifreeslip==1)then
      turb%sflux_u(i,j) = 0.0
      turb%sflux_v(i,j) = 0.0
      turb%sflux_w(i,j) = 0.0
      turb%sflux_t(i,j) = 0.0
      turb%sflux_r(i,j) = 0.0
     endif
   enddo
enddo

if (ilwrtyp > 0 .or. iswrtyp > 0) then
   do j = ja,jz
      do i = ia,iz
         radiate%albedt (i,j) = radiate%albedt (i,j) * dtll_factor
         radiate%rlongup(i,j) = radiate%rlongup(i,j) * dtll_factor
      enddo
   enddo
endif

!Update CO2 tendencies with SiB surface CO2 flux
if(isfcl==2) &
 CALL co2_biosource (m1,m2,m3,ia,iz,ja,jz,ngrid,np  &
           ,tend%rco2t(1),basic%dn0,grid%rtgt)

deallocate(ths2,rvs2,pis2,dens2,ups2,vps2,zts2)

return
END SUBROUTINE leaf3_sib

!##############################################################################
Subroutine co2_biosource (m1,m2,m3,ia,iz,ja,jz,ng,np,rco2t,dn0,rtgt)

use mem_grid
use mem_sib
use mem_leaf

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ng,np
integer :: i,j,ip
real, dimension(m1,m2,m3) :: rco2t,dn0
real, dimension(m2,m3) :: rtgt
real f,dz
data f /.044/ !CO2 molar mass = .044kg/mol

!Update CO2 tendencies with SiB surface CO2 flux
DO j=ja,jz
 DO i=ia,iz
  DO ip=2,np !Only done for land patches
    dz = rtgt(i,j)/dzt(2) ! dzt=1/(z(k)-z(k-1))
    rco2t(2,i,j) = rco2t(2,i,j) &
     + f * sib_g(ng)%co2flx(i,j,ip) * leaf_g(ng)%patch_area(i,j,ip) &
     / (dz*dn0(2,i,j))
  ENDDO
 ENDDO
ENDDO

return
END SUBROUTINE co2_biosource

!##############################################################################
Subroutine sfc_pcp_sib (nqparm,level,i,j,cuparm,micro,cupr,lspr)

use mem_basic
use mem_micro
use mem_cuparm
use leaf_coms

implicit none

integer :: nqparm,level,i,j
real :: cupr,lspr
type (cuparm_vars)  cuparm
type (micro_vars)   micro

  cupr = 0.
  lspr = 0.
  if (nqparm > 0) cupr = cuparm%conprr(i,j) ! convective precip at mm/s
  if (level >= 3) lspr = micro%   pcpg(i,j) ! micro precip at kg/m2

return
END SUBROUTINE sfc_pcp_sib
