!##############################################################################
!SiB VERSION 2.5/3, MODIFIED FOR USE WITH RAMS

!Note: "prog" below indicates a SiB prognostic variable
!      "diag" below indicated a SiB diagnostic variable
!      "output" below indicates an output variable passed back to RAMS
!         for use in the main model carbon and surface fluxes.
!      "atmos" indicates fields given by LEAF and/or RAMS for use
!         in driving SiB.
!      "BCs" below indicates soil and vegetation conditions set by the
!         user or MAKESFC surface file conditions from surface characteristics
!         data files.

Subroutine sib_2pt5 ( &
  !THIS SECTION: INPUT GRID AND ATMOSPHERIC CONDITIONS FROM RAMS TO SIB
   igrp        & !I grid point
  ,jgrp        & !J grid point
  ,mzg         & !number of soil levels
  ,mzs         & !number of snow levels
  ,dt          & !timestep
  ,doy         & !time in julian date (Day Of Year)
  ,ubmin       & !minimum wind speed for surface fluxes (m/s)
  ,latitude    & !latitude
  ,thm         & !atmos: -lowest layer theta (potential temperature K)
  ,sh_x        & !atmos: -lowest layer rv (water vapor mixing ratio, g/kg)
  ,ps_x        & !atmos: -lowest layer air pressure (hPa=mb) RAMS uses (Pa)
  ,ros         & !atmos: -lowest layer base state density (kg/m3)
  ,ts_x        & !atmos: -lowest layer temperature (K)
  ,spdm        & !atmos: -lowest layer wind speed (m/s)
  ,zwind       & !atmos: -lowest model level height t-grid (m)
  ,cupr        & !atmos: -cuparm precip rate (conprr) (mm/sec)
  ,lspr        & !atmos: -microphysics precip rate (pcpg / dt) (mm/sec)
  ,dswbot      & !atmos: -surface incident shortwave radiation (W/m^2)
  ,dlwbot      & !atmos: -surface incident longwave radiation (W/m^2)
  ,cosz        & !atmos: -cosine solar zenith angle
  !THIS SECTION: SURFACE FIELD INFO FROM RAMS TO SIB
  ,patcharea   & !BCs: -fractional patch area
  ,biome_f     & !BCs: -leaf vegetation class
  ,cndvi       & !BCs: -current month ndvi
  ,pndvi       & !BCs: -previous/past month ndvi
  ,soiltype_f  & !BSs: -soil textural class
  !THIS SECTION: PROGNOSTIC VALUES FROM SIB TO FEEDBACK TO RAMS
  ,rst         & !prog: -CANOPY/STOMATAL RESISTANCE (S M-1)
  ,ta          & !prog: -canopy (CAS) temperature (K)
  ,sha         & !prog: -canopy (CAS) water vapor mixing ratio (kg/kg)
  ,tc          & !prog: -vegetation temperature (K)
  ,tempk       & !prog: -soil and snow temperature (K)
  ,soil_water  & !prog: -volumetric soil moisture (m3/m3)
  ,sfcswa      & !prog: -surface abledo (fraction)
  ,uplwrf      & !prog: -upwelling longwave radiation (W/m2)
  ,snow1       & !prog: -vegetation snow (kg/m2)
  ,snow2       & !prog: -ground surface snow (kg/m2)
  ,capac1      & !prog: -vegetation liquid store (kg/m^2)
  ,capac2      & !prog: -ground surface liquid interception store (kg/m^2)
  ,pco2ap      & !prog: -canopy air space pCO2 concentration (Pa)
  ,co2flx      & !prog: -sfc CO2 flux between CAS and ref lev(mol/m2/sec)
  ,ustar       & !prog: -u-star
  ,rstar       & !prog: -r-star
  ,tstar       & !prog: -t-star
  !THIS SECTION: INPUT CURRENT CO2 CONENTRATION TO SIB
  ,atm_co2_ppm & !lowest atm level CO2, in ppm
  !THIS SECTION: DIAGNOSTIC QUANTITIES FROM SIB TO FILL STANDARD LEAF3 VARS
  ,sfcwaterlev & !diag: number of sfc water levels (0 or 1 from SiB)
  ,sfcwatermas & !diag: sfc water mass (snow+liquid) (kg/m2)
  ,sfcwaterdep & !diag: sfc water depth (snow+liquid) (m)
  ,vegwatermas & !diag: vegetation water mass (snow+liquid) (kg/m2)
  !THIS SECTION: DIAGNOSTIC QUANTITIES UNIQUE TO SIB
  ,vegalb_out  & !diag: vegetation albedo (fraction)
  ,vcover_out  & !diag: vegetation fractional area
  ,vlt_out     & !diag: vegetation LAI
  ,zlt_out     & !diag: total LAI
  ,veght_out   & !diag: canopy top height (m)
  ,vrough_out  & !diag: vegetation roughness (m)
  ,assimn_out  & !diag: net co2 assim by plants (umol/m^2/sec)
  ,respg_out   & !diag: ground respiration flux (umol/m^2/sec)
  ,rstfac1_out & !diag: CANOPY RESISTANCE STRESS1 : leaf sfc humidity
  ,rstfac2_out & !diag: CANOPY RESISTANCE STRESS2 : soil moisture
  ,rstfac3_out & !diag: CANOPY RESISTANCE STRESS3 : temperature
  ,ect_out     & !diag: transpiration flux (W/m^2)
  ,eci_out     & !diag: canopy interception flux (W/m^2)
  ,egi_out     & !diag: ground interception flux (W/m^2)
  ,egs_out     & !diag: ground surface layer evap (W/m^2)
  ,hc_out      & !diag: canopy (veg) sensible heat flux (W/m^2)
  ,hg_out      & !diag: ground surface sensible heat flux (W/m^2)
  ,ra_out      & !diag: CAS-RAMS aerodynamic resistance (sec/m)
  ,rb_out      & !diag: leaf sfc-CAS resistance (sec/m)
  ,rc_out      & !diag: total canopy resistance (sec/m)
  ,rd_out      & !diag: ground-CAS resistance (sec/m)
  ,roff_out    & !diag: runoff (surface and subsurface) (mm)
  ,green_out   & !diag: greenness fraction (-)
  ,apar_out    & !diag: absorbed fraction of PAR
  ,ventmf_out  & !diag: ventilation mass flux (kg/m^2/sec)
  ,pco2c_out   & !diag: leaf chloroplast CO2 concentration (Pa)
  ,pco2i_out   & !diag: leaf internal CO2 concentration (Pa)
  ,pco2s_out   & !diag: leaf surface CO2 concentration (Pa)
  ,pco2m_out   & !diag: lowest atm level CO2 concentration (Pa)
  ,ea_out      & !diag: CAS water vapor pressure (hPa)
  ,em_out      & !diag: ref level vapor pressure (hPa)
  ,rha_out     & !diag: CAS relative humidity
  ,radvbc_out  & !diag: radiation: visible beam (W/m^2)
  ,radvdc_out  & !diag: radiation: visible diffuse (W/m^2)
  ,radnbc_out  & !diag: radiation: nir beam (W/m^2)
  ,radndc_out  & !diag: radiation: nir diffuse (W/m^2)
  ,psy_out     ) !diag: psychrometric constant (hPa deg^-1)

use mem_grid, only:ngrid
use mem_sib
use node_mod, only:mi0,mj0

implicit none

  !Variables with the '_x' subscript were modified to avoid
  !conflicts with RAMS common block variables.

  !-------------------------------------------------------------------
  ! REFERENCES: Sato, N., P. J. Sellers, D. A. Randall, E. K. Schneider,
  !     J. Shukla, J. L Kinter III, Y-T, Hou, and Albertazzi (1989)
  !     "Effects of implementing the simple biosphere model in a general
  !     circulation model. J. Atmos. Sci., 46, 2767-2782.
  !            Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
  !     C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996) A revise
  !     land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
  !     Model formulation. (accepted by JCL)
  !  MODIFICATIONS:
  !   - changed VQSAT call to VNQSAT.  kwitt 10/23
  !   - added in the prognostic stomatal conductance in addinc. changan
  !   - moved sib diagnostics accumulation from dcontrol to bldif
  !     dd 950202
  !  ROUTINES called:  VNQSAT, SNOW1, balan, VNTLAT
  !       DELHF, DELEF, NETRAD, SIBSLV, endtem, updat2, addinc
  !       inter2, balan, soilprop, soiltherm, begtem, rnload
  !  FUNCS called:
  !       none
  !-------------------------------------------------------------------

  !FIXING LEN and NSIB AT 1-WILL RUN AS A SINGLE POINT.
  INTEGER, PARAMETER :: len = 1  ! Run as single point with RAMS interface
  INTEGER, PARAMETER :: nsib = 1 ! Run as single point with RAMS interface
  INTEGER, PARAMETER :: nsoil = nzg_sib - 1  ! number of soil layers
  INTEGER, PARAMETER :: ioffset = 0 ! subdomain offset
  INTEGER :: sibprint,i,j,k,n,ksoil,l

  !PI
  REAL, PARAMETER :: num_pi = 3.1415926
  !gravity (m/s2)
  REAL, PARAMETER :: grav = 9.81
  !specific heat of air at const pres (J kg-1 deg-1)
  REAL, PARAMETER :: cp = 1004.
  !specific heat of air at const volume (J kg-1 deg-1)
  REAL, PARAMETER :: cv = 1952.
  !universal gas constant
  REAL, PARAMETER :: rgas = 287.
  REAL, PARAMETER :: hltm = 2.52E6
  REAL, PARAMETER :: delta = 0.608
  !conversion for kg water to snow depth (16.7)
  REAL, PARAMETER :: asnow = 16.7
  !R/cp
  REAL, PARAMETER :: kapa = 0.2861328125
  !latent heat of fusion for ice (J m^-3)
  REAL, PARAMETER :: snomel = 3.705185e8
  REAL, PARAMETER :: clai = 4.186*1000.0*0.2
  REAL, PARAMETER :: cww = 4.186*1000.0*1000.
  !von karman constant
  REAL, PARAMETER :: vkrmn = 0.35
  REAL, PARAMETER :: po2m = 20900.
  !stefan-boltzmann constant
  REAL, PARAMETER :: stefan = 5.67e-8
  REAL, PARAMETER :: grav2 = grav *0.01
  REAL, PARAMETER :: tice = 273.16
  REAL, PARAMETER :: snofac = hltm / ( hltm + SNOMEL * 1.E-3 )

  !VARIABLES INPUT/OUTPUT FROM/TO RAMS MODEL
  INTEGER :: mzs,mzg,igrp,jgrp,doy
  REAL, DIMENSION(mzg+mzs) :: tempk
  REAL, DIMENSION(mzg) :: soil_water
  REAL ::             &
        dt            &! time step (s)
       ,ubmin         &! minimum wind speed for surface fluxes (m/s)
       ,latitude      &! latitude of current grid cell
       ,thm(len)      &! mixed layer potential temperature (K)
       ,sh_x(len)     &! mixed layer water vapor mixing ratio (kg/kg)
       ,ps_x(len)     &! surface pressure (hPa=mb) RAMS uses (Pa)
       ,ros (len)     &! surface air density (kg/m^3)
       ,ts_x(len)     &! surface mixed layer air temperature (K)
       ,spdm(len)     &! boundary layer wind speed (m/s)
       ,zwind         &! lowest model level height t-grid (m)
       ,cupr(len)     &! convective-parm precipitation rate (mm/s)
       ,lspr(len)     &! microphysics precipitation rate (mm/s)
       ,dlwbot(len)   &! surface incident longwave radiation (W/m^2)
       ,dswbot(len)   &! surface incident shortwave radiation (W/m^2)
       ,cosz(len)     &! cosine of solar zenith angle
       ,patcharea     &! land patch fractional area
       ,biome_f       &! leaf vegetation class
       ,pndvi         &! past    value of ndvi for the gridcell
       ,cndvi         &! current value of ndvi for the gridcell
       ,soiltype_f    &! soil textural class
       ,rst(len)      &! canopy/stomatal resistance (S M-1)
       ,ta(len)       &! CAS temperature (K)
       ,sha(len)      &! CAS water vapor mixing ratio (kg/kg)
       ,tc(len)       &! canopy (vegetation) temperature (K)
       ,sfcswa        &! surface abledo (fraction)
       ,uplwrf        &! upwelling longwave radiation (W/m2)
       ,snow1(len)    &! vegetation snow (kg/m2)
       ,snow2(len)    &! ground surface snow (kg/m2)
       ,capac1(len)   &! vegetation liquid store (kg/m^2)
       ,capac2(len)   &! ground surface liquid interception store (kg/m^2)
       ,pco2ap(len)   &! canopy air space pCO2 (Pa)
       ,co2flx(len)   &! sfc CO2 flux between CAS and ref lev(mol/m^2/sec)
       ,ustar(len)    &! friction velocity (m/s)
       ,rstar(len)    &! moisture exchange parameter based on U* (kg/kg)/(m/s)
       ,tstar(len)    &! heat exchange parameter based on U* (K)/(m/s)
       ,atm_co2_ppm   &! lowest atm level CO2, in ppm
       ,sfcwaterlev   &! number of surface water layers (0 or 1)
       ,sfcwatermas   &! sfc water mass (snow+liquid) (kg/m2)
       ,sfcwaterdep   &! sfc water depth (snow+liquid) (m)
       ,vegwatermas   &! vegetation water mass (snow+liquid) (kg/m2)
       ,vegalb_out    &!diag: vegetation albedo (fraction)
       ,vcover_out    &!diag: vegetation fractional area
       ,vlt_out       &!diag: vegetation LAI
       ,zlt_out       &!diag: total LAI
       ,veght_out     &!diag: canopy top height (m)
       ,vrough_out    &!diag: vegetation roughness (m)
       ,assimn_out    &!diag: net co2 assim by plants (umol/m^2/sec)
       ,respg_out     &!diag: ground respiration flux (umol/m^2/sec)
       ,rstfac1_out   &!diag: CANOPY RESISTANCE STRESS1 :leaf sfc humidity
       ,rstfac2_out   &!diag: CANOPY RESISTANCE STRESS2 :soil moisture
       ,rstfac3_out   &!diag: CANOPY RESISTANCE STRESS3 :temperature
       ,ect_out       &!diag: transpiration flux (W/m^2)
       ,eci_out       &!diag: canopy interception flux (W/m^2)
       ,egi_out       &!diag: ground interception flux (W/m^2)
       ,egs_out       &!diag: ground surface layer evap (W/m^2)
       ,hc_out        &!diag: canopy (veg) sensible heat flux (W/m^2)
       ,hg_out        &!diag: ground surface sensible heat flux (W/m^2)
       ,ra_out        &!diag: CAS-RAMS resistance (sec/m)
       ,rb_out        &!diag: leaf sfc-CAS resistance (sec/m)
       ,rc_out        &!diag: total canopy resistance (sec/m)
       ,rd_out        &!diag: ground-CAS resistance (sec/m)
       ,roff_out      &!diag: runoff (surface and subsurface) (mm)
       ,green_out     &!diag: greenness fraction (-)
       ,apar_out      &!diag: absorbed fraction of PAR
       ,ventmf_out    &!diag: ventilation mass flux (kg/m^2/sec)
       ,pco2c_out     &!diag: leaf chloroplast CO2 concentration (Pa)
       ,pco2i_out     &!diag: leaf internal CO2 concentration (Pa)
       ,pco2s_out     &!diag: leaf surface CO2 concentration (Pa)
       ,pco2m_out     &!diag: lowest atm level CO2 concentration (Pa)
       ,ea_out        &!diag: CAS water vapor pressure (hPa)
       ,em_out        &!diag: ref level vapor pressure (hPa)
       ,rha_out       &!diag: CAS relative humidity
       ,radvbc_out    &!diag: radiation: visible beam (W/m^2)
       ,radvdc_out    &!diag: radiation: visible diffuse (W/m^2)
       ,radnbc_out    &!diag: radiation: nir beam (W/m^2)
       ,radndc_out    &!diag: radiation: nir diffuse (W/m^2)
       ,psy_out        !diag: psychrometric constant (hPa deg^-1)

  !SIB STATIC SURFACE PARAMETERS
  REAL ::             &
        z0d(len)      &! surface roughness length (m)
       ,z0(len)       &! sfc rough length corrected for canopy snow (m)
       ,zlt(len)      &! leaf area index
       ,z1(len)       &! canopy bottom height (meters)
       ,z2(len)       &! canopy top height (meters)
       ,cc1(len)      &! RB Coefficient (c1) = rbc
       ,cc2(len)      &! RC Coefficient (c2) = rdc
       ,dd(len)       &! Zero plane displacement
       ,poros(len)    &! soil porosity
       ,zdepth(len,3) &! porosity * soil hydrology model layer depths (m)
       ,phsat(len)    &! Soil tension at saturation (units)
       ,bee(len)      &! Clapp & Hornberge 'B' exponent
       ,respcp(len)   &! respiration fraction of Vmax
       ,vmax0(len)    &! rubisco velocity of sun leaf (mol/m2/s)
       ,green(len)    &! Canopy greeness fraction of LAI
       ,tran(len,2,2) &! Leaf transmittance
       ,ref(len,2,2)  &! Leaf reflectance
       ,gmudmu        &! Time-mean leaf projection (leaf orientation to par flux)
       ,trop(len)     &! temperature coefficient in GS-A model (K)
       ,phc(len)      &! one-half critical leaf-water potential limit(m)
       ,trda(len)     &! slope of high temp inhibition (leaf resp,1/K)
       ,trdm(len)     &! half point of high temp inhibition (leaf resp,K)
       ,slti(len)     &! slope of low temperature inhibition (1/K)
       ,shti(len)     &! slope of high temperature inhibition (1/K)
       ,hlti(len)     &! half piont of low temp inhibition (K)
       ,hhti(len)     &! half point of high temp inhibition (K)
       ,effcon(len)   &! quantum efficiency (mol/mol)
       ,binter(len)   &! conductance-photosynthesis intercept (mol/m2/s)
       ,gradm(len)    &! conductance-photosynthesis slope parm (mol/m2/s)
       ,atheta(len)   &! wc,we coupling parameter
       ,btheta(len)   &! wp,ws coupling parameter
       ,aparc(len)    &! Canopy absorbed fraction of PAR
       ,wopt(len)     &! Factor coeff for moisture effect on soil respiration
       ,zmx(len)      &! Power coeff for moisture effect on soil respiration
       ,wsat(len)     &! respiration at soil water saturation?
       ,vcover(len)   &! vegetation cover fraction
       ,sodep(len)    &! total soil depth (meters)
       ,rootd(len)    &! rooting depth (meters)
       ,soref(len,2)  &! Soil reflectance
       ,thermk(len)   &! canopy gap fraction for TIR radiation
       ,satco(len)    &! soil tension at 1/2 assimilation value (true?) units?
       ,slope(len)    &! slope
       ,chil(len)     &! leaf angle distribution factor
       ,ztdep(len,nsoil) ! soil thermal model layer depths (m)

  !    Intent: in/out and some variable local copies
  REAL ::                 & 
        cas_cap_heat(len) &! CAS heat capacity (J/K m^2)
       ,cas_cap_vap(len)  &! CAS vapor capacity
       ,cas_cap_co2(len)  &! CAS CO2 capacity (m/m^2) (moles air / m^2 in phosib)
       ,tg(len)           &! surface boundary temperature (K)
       ,td(len,nsoil)     &! deep soil temperature (K)
       ,www(len,3)        &! soil wetness
       ,wwwtem(len,3)     &! soil wetness copy
       ,snow(len,2)       &! snow cover (kg/m2)
       ,capac(len,2)      &! liquid interception store (kg/m^2)
       ,cuprt(len)        &! copy of cupr
       ,lsprt(len)        &! copy of lspr
       ,thmtem(len)       &! copy of thm
       ,shtem(len)        &! copy of sh
       ,zzwind(len)       &! intermediate surface wind for roughness length
       ,zztemp(len)       &! intermediate surface wind for roughness length
       ,ztemp              ! Used for ratio of reference height (zwind/ztemp)

  !respFactor is the annual total accumulation of carbon in the
  !previous year at each grid cell (annual total ASSIMN).
  !divided by the annual total of soilScale at the same grid pt.
  !respFactor*soilScale is the rate of release of CO2 by the soil.
  REAL :: respfactor(len,nzg_sib)
  !soilScale is a diagnostic of the instantaneous rate of
  !soil respiration (derived by Jim Collatz, similar to TEM)
  REAL :: soilscale(len,nzg_sib)
  REAL :: soilq10(len,nzg_sib)

  !MOSTLY VARIABLES FOR SURFACE FLUXES
  REAL ::                 & 
        fss(len)          &! surface sensible heat flux (W/m^2)
       ,fws(len)          &! surface evaporation (kg/m^2/s)
       ,cflux(len)        &! new formulation of CO2 flux (phosib)(mol/m^2/sec)
       ,cu(len)           &! momentum transfer coefficient (-)
       ,ct(len)           &! thermal transfer coefficient (-)
       ,ventmf(len)       &! ventilation mass flux (kg/m^2/sec)
       ,thvgm(len)        &! delta theta-v between atm and canopy
       ,bps(len)          &! (ps/1000)**kapa
       ,psy(len)          &! psychrometric 'constant'
       ,tha(len)          &! canopy airspace potential temperature (K)
       ,ea(len)           &! canopy airspace water vapor pressure (hPa)
       ,em(len)           &! mixed layer water vapor pressure (hPa)
       ,d(len)            &! dd corrected for snow covered canopy
       ,rbc(len)          &! cc1 corrected for snow covered canopy
       ,rdc(len)          &! cc2 corrected for snow covered canopy
       ,etmass(len)       &! evapotranspiration
       ,egmass(len)       &! ground evaporation (mm)
       ,ecmass(len)       &! canopy evaporation (mm)
       ,totwb(len)        &! total surface and soil water at begiN of timestep
       ,chf(len)          &! canopy heat flux (W/m^2)
       ,ahf(len)          &! CAS heat flux (W/m^2)
       ,shf(len)          &! soil heat flux (W/m^2)
       ,ect(len)          &! transpiration flux (J m^-2 for the timestep)
       ,eci(len)          &! canopy interception evap flux (veg-CAS) (J m^-2)
       ,egs(len)          &! soil/ground evaporation flux (J m^-2)
       ,egi(len)          &! ground interception evaporation flux (J m^-2)
       ,hc(len)           &! canopy (veg) sensible heat flux (W/m^2)
       ,hg(len)           &! ground surface sensible heat flux (W/m^2)
       ,hs(len)           &! snow surface sensible heat flux (W/m^2)
       ,heaten(len)       &! energy to heat snow to ground temp (J m^-2)
       ,hflux(len)        &! sensible heat flux (W/m2)
       ,tgs(len)          &! bare ground and snow surface mean temperature (K)
       ,tsnow(len)        &! Snow temp is lesser of ice or ground temp (K)
       ,czc(len)          &! canopy heat capacity
       ,etc(len)          &! vapor pressure (e*) of the canopy at (Tc) (Pa)
       ,etg(len)          &! vapor pressure (e*) of the ground sfc at (Tg) (Pa)
       ,etgs(len)         &! function result = E(TGS(i))
       ,btc(len)          &! derivatives of ETC
       ,btg(len)          &! derivatives of ETG
       ,rstfac(len,4)     &! canopy resistance stress factors
       !GROUND AND CANOPY STRESS AND PHOSIB VARIABLES
       ,rsoil(len)        &! SOIL SURFACE RESISTANCE (S M-1)
       ,hr(len)           &! SOIL SURFACE RELATIVE HUMIDITY
       ,wc(len)           &! CANOPY WETNESS FRACTION
       ,wg(len)           &! GROUND WETNESS FRACTION
       ,areas(len)        &! fractional snow coverage (0 to 1)
       ,csoil(len)        &! soil heat capacity (J m^-2 deg^-1)
       ,slamda(len,nsoil) &! soil thermal conductivities
       ,shcap(len,nsoil)  &! soil heat capacities
       ,czh(len)          &! surface layer heat capacity
       ,gect(len)         &! dry fraction of canopy/(Rst + 2Rb)
       ,geci(len)         &! wetted fraction of canopy/2Rb
       ,gegs(len)         &! dry fraction of ground/(fg*rsoil + Rd)
       ,gegi(len)         &! wet fraction of ground/Rd
       ,rb(len)           &! leaf sfc to CAS aerodynmaic resistance (sec/m)
       ,rd(len)           &! ground to CAS aerodynamic resistance (sec/m)
       ,rds(len)          &! rsoil + rd
       ,ra(len)           &! CAS to RAMS atmos aerodynamic resistance (sec/m)
       ,rib(len)          &! bulk richardson number
       ,rc(len)           &! total canopy resistance (sec/m)
       ,ggl(len)          &! overall leaf conductance
       ,hrr(len)          &! SOIL SURFACE LAYER RELATIVE HUMIDITY
       ,bintc(len)        &! (B*ZLT)  : EQUATION (35) , SE-92A
       ,aparkk(len)       &! (PI)     : EQUATION (31) , SE-92A
       ,wsfws(len)        &! Water stress
       ,wsfht(len)        &! High temperature stress
       ,wsflt(len)        &! Low temperature stress
       ,wci(len)          &! Intermediate assimilation weighted Ci
       ,whs(len)          &! Intermediate assimilation weighted RH stress factor
       ,wags(len)         &! Intermediate assimilation weighted stomatal conductance
       ,wegs(len)         &! Intermediate evaporation weighted stomatal conductance
       ,omepot(len)       &! Potential light limitation
       ,assim(len)        &! gross primary productivity (mol/m^2/s)
       ,assimpot(len)     &! Final Potential top leaf photosynthesis
       ,assimci(len)      &! Final stress limited top leaf photosynthesis
       ,assimnp(len)      &! Make assimn a top leaf, not the canopy
       ,antemp(len)       &! Bottom stopped assimilation
       ,ansqr(len)        &! Bottom stopped assimilation (squared)
       ,pfd(len)          &! 4.6E-6 * GMUDMU * (RADN(i,1,1)+RADN(i,1,2))
       ,zmstscale(len,2)  &! soil scaling parameter for shallow and root zone
       ,drst(len)         &! stomatal resistance increment
       !SENSIBLE HEAT FLUX DERIVATIVES AND SUCH
       ,dtg(len,2)        &! surface ground and snow temperature increments (K)
       ,dtc(len)          &! canopy temperature increment (K)
       ,dta(len)          &! CAS temperature increment (K)
       ,hgdtg(len)        &! dHG/dTG
       ,hgdta(len)        &! dHg/dTa
       ,hsdts(len)        &! dHS/dTS
       ,hsdta(len)        &! dHS/dTA
       ,hcdtc(len)        &! dHc/dTc
       ,hcdta(len)        &! dHc/dTa
       ,hadta(len)        &! dHA/dTA
       ,fc(len)           &! canopy range function? (not used here)
       ,fg(len)           &! ground humidity range function (0 or 1)
       !LONGWAVE RADIATIVE HEAT FLUX DERIVATIVES AND SUCH
       ,lcdtc(len)        &! dLC/dTC
       ,lcdtg(len)        &! dLC/dTG
       ,lcdts(len)        &! dLC/dTS
       ,lgdtg(len)        &! dLG/dTG
       ,lgdtc(len)        &! dLG/dTC
       ,lsdts(len)        &! dLS/dTS
       ,lsdtc(len)        &! dLS/dTC
       !LATENT HEAT FLUX GROUND AND CANOPY PARTIAL DERIVATIVE AND SUCH
       ,eg(len)           &! EGS + EGI
       ,ec(len)           &! ECT + ECI
       ,es(len)           &
       ,egdtg(len)        &! dEG/dTGS
       ,ecdtg(len)        &! dEC/dTGS
       ,ecdtc(len)        &! dEC/dTC
       ,egdtc(len)        &! dEG/dTC
       ,ecdea(len)        &! for the canopy leaves vapor pressure: W/ (m2* K)
       ,egdea(len)        &! for ground latent heat fluxes: W/ (m2* K)
       ,esdts(len)        &! for snow latent heat fluxes: W/ (m2* K)
       ,esdea(len)        &! for snow latent heat fluxes: W/ (m2 * Pa)
       ,eadea(len)        &! for CAS latent heat fluxes: W/ (m2* Pa)
       ,radt(len,3)       &! canopy, ground, and snow net radiation (W/m2)
       !INCREMENTS FOR INTEGRATION
       ,dea(len)          &! CAS moisture increment (Pa)
       ,dtd(len,nsoil)    &! deep soil temperature increments (K)
       ,q3l(len)          &! 'Liston' drainage from bottom of soillayer 3 (mm)
       ,q3o(len)          &! gravitational drainage out of soillayer 3 (mm)
       ,qqq(len,3)        &! soil layer drainage (mm m^-2 timestep)
       ,evt(len)          &
       ,eastar(len)       &! canopy saturation vapor pressure (hPa)
       ,rha(len)          &! canopy airspace relative humidity (%)
       !SURFACE WATER VARS
       ,zmelt(len)        &! total depth of melted water (m)
       ,zmelt1(len)       &! depth of melted water, main calculation (updat2)(m)
       ,zmelt2(len)       &! depth of melted water, from excess energy (inter2)(m)
       ,satcap(len,2)     &! saturation capacity of ground and vegetation (kg/m^2)
       ,exo(len)          &! total soil water excess of saturation (m)
       ,roffo(len)        &! runoff overland flow contribution (m)
       ,roff(len)         &! total runoff (surface and subsurface) (mm)
       !RADIATION VARIABLES
       ,salb(len,2,2)     &! surface albedos
       ,valb(len,2,2)     &! vegetation albedos
       ,nalb(len,2,2)     &! non-veg albedos
       ,canex(len)        &! [1.-( SNOWw(i,1)*5.-Z1(i))/(Z2(i)-Z1(i))]
       ,fac1(len)         &! effective ground cover for thermal radiation
       ,thgeff(len)       &! Tgeff(I) / BPS(I)
       ,tgeff(len)        &! effective (combined) skin temp from sfc thermal rad(K)
       ,shgeff(len)       &! saturation mixing ratio w.r.t tgeff
       ,tgeff4(len)       &! effective surface radiative temperature (K)
       ,radfac(len,2,2,2) &! radiation absorption factors
       ,radvbc(len)       &! surface incident visible direct beam (W/m^2)
       ,radnbc(len)       &! surface incident near IR direct beam (W/m^2)
       ,radvdc(len)       &! surface incident visible diffuse beam (W/m^2)
       ,radndc(len)       &! surface incident near IR diffuse beam (W/m^2)
       ,radc3(len,2)      &! SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
       ,radn(len,2,2)     &! INCIDENT RADIATION FLUXES (W M-2)
       ,closs(len)        &! vegetation IR loss
       ,gloss(len)        &! ground IR loss
       ,sloss(len)        &! snow IR loss
       ,dtc4(len)         &! 1st derivative of vegetation T^4
       ,dtg4(len)         &! 1st derivative of ground T^4
       ,dts4(len)         &! 1st derivative of snow T^4
       !RESPIRATION VARIABLES
       ,assimn(len)       &! net co2 assimilation by plants (mol/m^2/sec)
       ,respg(len)        &! ground respiration flux (mol/m^2/sec)
       ,pco2i(len)        &! leaf internal pCO2 (Pa)
       ,pco2c(len)        &! chloroplast pCO2 (Pa)
       ,pco2s(len)        &! leaf surface pCO2 (Pa)
       ,pco2m(len)        &! lowest atm level CO2 (Pa)
       ,co2cap(len)        ! moles of air in the canopy (moles/CAS)

  TYPE biome_morph_var
     REAL zc        ! Canopy inflection height (m)
     REAL lwidth    ! Leaf width
     REAL llength   ! Leaf length
     REAL laimax    ! Maximum LAI
     REAL stems     ! Stem area index
     REAL ndvimax   ! Maximum NDVI
     REAL ndvimin   ! Minimum NDVI
     REAL srmax     ! Maximum simple ratio
     REAL srmin     ! Minimum simple ratio
  END TYPE biome_morph_var
  TYPE(biome_morph_var) morphtab

  TYPE aero_var
     REAL zo       ! Canopy roughness coeff
     REAL zp_disp  ! Zero plane displacement
     REAL rbc      ! RB Coefficient
     REAL rdc      ! RC Coefficient
  END TYPE aero_var

  ! aerodynamic interpolation tables
  TYPE(aero_var),DIMENSION(50,50) :: aerovar

  TYPE time_dep_var
     REAL fpar    ! Canopy absorbed fraction of PAR
     REAL lai     ! Leaf-area index
     REAL green   ! Canopy greeness fraction of LAI
     REAL zo      ! Canopy roughness coeff
     REAL zp_disp ! Zero plane displacement
     REAL rbc     ! RB Coefficient (c1)
     REAL rdc     ! RC Coefficient (c2)
     REAL gmudmu  ! Time-mean leaf projection
  END TYPE time_dep_var
  TYPE(time_dep_var) timevar

  !Send in vegetation biosphere and soil class info and coordinate
  !with SiB equivalent classes of vegetation and soil. Note that
  !RAMS clases of vegetation go from 0:20 or 21 total classes, but
  !RAMS LEAF class=0 is ocean. Do not need this for SiB.
  !Zero is water. Should not be running SiB for water patches.
  !Alert model to report error.
  INTEGER :: biome   &  ! biome type, sent in as a single value
            ,soiltype   ! soil type, sent in as a single value
  INTEGER, DIMENSION(0:20) :: leaf_biome_map
  DATA leaf_biome_map /0,0,13,11,4,5,2,1,6,6,9,10,9,9,3,12,12,12,7,11,1/

  !LEAF soil classes do not match the order of classes in SiB so we need
  !to translate these correctly. LEAF uses the classes from Clapp & Hornberger
  !and SiB uses modified USDA classes. LEAF use of FAO data only assigns 
  !classes 2-8 (Clapp & Hornberger) via routine "datp_datsoil". Others
  !default to sandy-clay-loam.
  INTEGER, DIMENSION(1:12) :: soil_type_map
  DATA soil_type_map /1,2,3,4,6,7,10,9,8,11,12,7/

  !Minimum surface wind speed for generating surface fluxes (Louis)
  spdm=MAX(spdm,ubmin)

  !Conversions for sib from leaf_class to biome type
  biome = leaf_biome_map(INT(biome_f))
  if(biome==0) then
   print*,"SiB trying to run with BIOME = WATER. Cannot Do This!"
   stop
  endif
  !Conversion for sib from leaf soil type to sib soil type
  soiltype = soil_type_map(INT(soiltype_f))

  !Convert some RAMS soil model things to SiB soil model layers.
  !SiB has 6 soil levels + separate ground surface level to get nzg_sib=7.
  !RAMS has soil levels only.

  !Initialize permanent wetlands as assigned by LEAF3_INIT with
  !surface water mass of 100.0 kg/m2 (0.1 meters depth). Keep this
  !constant and do not let it dry out. SiB has either snow or water,
  !so place surface water in correct category.
  DO i=1,len
   if(INT(biome_f)==17 .or. INT(biome_f)==20) then
     if(snow2(i)>0.) then
      snow2(i) = max(snow2(i),100.)
     else
      capac2(i) = max(capac2(i),100.)
     endif
     !print*,'surfwater',igrp+mi0(ngrid),jgrp+mj0(ngrid),int(biome_f) &
     !  ,biome,capac2(i)
   endif
  ENDDO

  !SiB soil levels: nzg_sib=7 and nsoil=(nzg_sib-1)=6
  !Layer 1 is the deepest layer, layer nsoil is immediately below the surface.
  !Soil layer nsoil (top layer) use www(1)
  !Soil layers 3 through nsoil-1 (root zone) use www(2) and TD(3:nsoil-1)
  !Soil layers 1 and 2 (deepest layers) use www(3)
  DO i=1,len
     !Sib uses 3 soil water layer while RAMS is flexible
     !RAMS bottom was layer is layer-1. SiB bottom water layer
     !is layer-3. Exchange these accordingly.
     www(i,1) = soil_water(7) !Top soil moisture layer (m3/m3)
     www(i,2) = soil_water(5) !Middle soil moisture layer (m3/m3)
     www(i,3) = soil_water(1) !Bottom soil moisture layer (m3/m3)
     !Set Sib ground surface tempeature to RAMS top soil layer temp (K)
     tg(i) = tempk(mzg)
     !Set Sib soil layer temps to RAMS soil layer temps (K).
     !The k=1 level is the lowest soil level.
     DO k=1,nsoil
        td(i,k) = tempk(k)
     ENDDO
  ENDDO

  !Ppm in a gas is normally expressed on a mole fraction basis, so 385 
  !ppm of CO2 (in the atmosphere) is also 385 µmol/mol. By Dalton's law 
  !of partial pressure, the partial pressure of CO2 to total pressure 
  !will be in the same ratio. "Real" atmospheric pressure varies with 
  !altitude above sea level, and day-to-day with the weather; however,
  !"standard" atmospheric pressure is 101.325 kPa. For example, the 
  !partial pressure of CO2 is:
  !385 µmol/mol x 101.325 kPa x (unit conversion) = 39 Pa. (at least 
  !on a day when real and standard pressure are the same value). Here,
  !assign value from atm_co2_ppm (lowest atm level CO2 in ppm) to pco2m
  !convert PPM to PA assuming atm pressure of 1000mb (or 100.0 kPa)
  pco2m(1) = atm_co2_ppm * ps_x(1) * 1.e-4 !convert to Pa

  !Copy input arrays to local arrays
  DO i=1,len
     bps(i) = (ps_x(i)/1000.)**kapa  ! (sfcpressure/1000)**kapa
     snow(i,1) = snow1(i)    ! vegetation snow (kg/m2)
     snow(i,2) = snow2(i)    ! ground surface snow (kg/m2)
     capac(i,1) = capac1(i)  ! vegetation liquid store (kg/m^2)
     capac(i,2) = capac2(i)  ! ground surface liquid interception store (kg/m^2)
  ENDDO

  !Set up some key flags for runtime
  ztemp = zwind          ! Lowest level wind used for ratio of reference 
                         !   height (zwind/ztemp). Set equal for SiB-RAMS

  !First step: obtain TI BCs using biome type, etc
  !this won't be exactly like was done in Owen's
  !version, because I will already have a map
  !of biome,soil,ndvi, etc read in-single values
  !will be passed to here.

  DO i=1,len

     !Soil BC's assigned in sib2_init.f90
     bee(i)   = bee_sib(soiltype)    !Clapp & Hornberge 'B' exponent
     phsat(i) = phsat_sib(soiltype)  !Soil tension at saturation (units)?
     satco(i) = satco_sib(soiltype)  !soil tension at 1/2 assimilation value (units)?
     poros(i) = poros_sib(soiltype)  !porosity
     slope(i) = slope_sib(soiltype)  !slope
     wopt(i)  = wopt_sib(soiltype)   !Factor coeff for moisture effect on soil respiration
     zmx(i)   = skew_sib(soiltype)   !Power coeff for moisture effect on soil respiration
     wsat(i)  = respsat_sib(soiltype)!respiration at soil water saturation?

     !Vegetation BC'S assigned in sib2_init.f90
     z2(i)     = z2_sib(biome)      !canopy top height (meters)
     z1(i)     = z1_sib(biome)      !canopy bottom height (meters)
     vcover(i) = fvcover_sib(biome) !fractional vegetation cover
     chil(i)   = chil_sib(biome)    !leaf angle distribution factor
     sodep(i)  = sodep_sib(biome)   !total soil depth (meters)
     rootd(i)  = rootd_sib(biome)   !rooting depth (meters)
     phc(i)    = phc_sib(biome)     !one-half critical leaf-water potential limit(m)

     !Leaf transmittance and reflectance assigned in sib2_init.f90
     tran(i,1,1)  = tran_sib(biome,1,1)
     tran(i,2,1)  = tran_sib(biome,2,1)
     tran(i,1,2)  = tran_sib(biome,1,2)
     tran(i,2,2)  = tran_sib(biome,2,2)
     ref(i,1,1)   = ref_sib(biome,1,1)
     ref(i,2,1)   = ref_sib(biome,2,1)
     ref(i,1,2)   = ref_sib(biome,1,2)
     ref(i,2,2)   = ref_sib(biome,2,2)

     !Assigned in sib2_init.f90
     vmax0(i)     = vmax0_sib(biome)  !rubisco velocity of sun leaf (mol/m2/s)
     effcon(i)    = effcon_sib(biome) !quantum efficiency (mol/mol)
     gradm(i)     = gslope_sib(biome) !conductance-photosynthesis slope parm (mol/m2/s)
     binter(i)    = gsmin_sib(biome)  !conductance-photosynthesis intercept (mol/m2/s)
     atheta(i)    = atheta_sib(biome) !wc,we coupling parameter
     btheta(i)    = btheta_sib(biome) !wp,ws coupling parameter
     trda(i)      = trda_sib(biome)   !slope of high temp inhibition (leaf resp,1/K)
     trdm(i)      = trdm_sib(biome)   !half point of high temp inhibition (leaf resp,K)
     trop(i)      = trop_sib(biome)   !temperature coefficient in GS-A model (K)
     respcp(i)    = respcp_sib(biome) !respiration fraction of Vmax
     slti(i)      = slti_sib(biome)   !slope of low temperature inhibition (1/K)
     hlti(i)      = hlti_sib(biome)   !half piont of low temp inhibition (K)
     shti(i)      = shti_sib(biome)   !slope of high temperature inhibition (1/K)
     hhti(i)      = hhti_sib(biome)   !half point of high temp inhibition (K)

     !Soil reflectance assigned in sib2_init.f90
     soref(i,1)   = soref_sib(biome,1)
     soref(i,2)   = soref_sib(biome,2)

  ENDDO

  morphtab%zc = zc_sib(biome)        !Canopy inflection height(m) (nvtyp_sib)
  morphtab%lwidth = zlw_sib(biome)   !Leaf Width (nvtyp_sib)
  morphtab%llength = zlen_sib(biome) !Leaf Length (nvtyp_sib)
  morphtab%laimax = ltmax_sib(biome) !Max LAI (nvtyp_sib)
  morphtab%stems = stem_sib(biome)   !Stem area index (nvtyp_sib)
  morphtab%ndvimax = nd98_sib(biome) !Max NDVI (98%tile) (nvtyp_sib)
  morphtab%ndvimin = nd02_sib(biome) !Min NDVI (2%tile) (nvtyp_sib)
  morphtab%srmax = srmax_sib(biome)  !Maximum simple ratio (nvtyp_sib)
  morphtab%srmin = srmin_sib(biome)  !Minimum simple ratio (nvtyp_sib)

  DO j=1,50
     DO i=1,50
        aerovar(i,j)%zo      = a_zo_sib(biome,i,j)  ! Canopy roughness coeff
        aerovar(i,j)%zp_disp = a_zp_sib(biome,i,j)  ! Zero plane displacement
        aerovar(i,j)%rbc     = a_rbc_sib(biome,i,j) ! RB Coefficient (c1)
        aerovar(i,j)%rdc     = a_rdc_sib(biome,i,j) ! RC Coefficient (c2)
     ENDDO
  ENDDO

  !Calculates time dependant surface boundary condition variables
  i=1
  CALL mapper (latitude,doy,pndvi,cndvi,vcover(i) &
       ,chil(i),tran(i,1,1),ref(i,1,1),morphtab   &
       ,aerovar,laig_sib,fvcg_sib,timevar)

  DO i=1,len
     aparc(i) = timevar%fpar    ! Canopy absorbed fraction of PAR
     zlt(i)   = timevar%lai     ! Leaf-area index
     green(i) = timevar%green   ! Canopy greeness fraction of LAI
     z0d(i)   = timevar%zo      ! Canopy roughness coeff
     dd(i)    = timevar%zp_disp ! Zero plane displacement
     cc1(i)   = timevar%rbc     ! RB Coefficient (c1)
     cc2(i)   = timevar%rdc     ! RC Coefficient (c2)
  ENDDO

  !respFactor is the annual total accumulation of carbon in the
  !previous year at each grid cell (annual total ASSIMN)
  !divided by the annual total of soilScale at the same grid pt.
  !itb_cptec...now hardwiring for RJ
  DO i=1,len
     respfactor(i,1) = 0.0
     respfactor(i,2) = 0.0
     respfactor(i,3) = 3.65e-7/(8.-zlt(i))
     respfactor(i,4) = 7.62e-7/(8.-zlt(i))
     respfactor(i,5) = 8.03e-7/(8.-zlt(i))
     respfactor(i,6) = 4.33e-6/(8.-zlt(i))
     respfactor(i,7) = 3.40e-6/(8.-zlt(i))
  ENDDO

  !Get radiation components from single SW value
  CALL raddrv (nsib,dswbot,cosz,radvbc,radvdc,radnbc,radndc)

  !Set up soil information based on table in sib2_init.f90
  !where total soil depth (sodep) (m) and rooting depth (rootd) (m)
  !vary by vegetation type. Taller vegetation is given deeper soil
  !depth and root depth. 
  DO i=1,len
     rootd(i) = MIN(rootd(i),(sodep(i)*0.75))
     !Soil moisture layers
     ! zdepth = soil hydrology model layer depth * porosity (meters)
     zdepth(i,1) = 0.02 * poros(i)
     zdepth(i,2) = (rootd(i) - 0.02)*poros(i)
     zdepth(i,3) = poros(i)*sodep(i) - (zdepth(i,1)+zdepth(i,2))
     ! Quick patch to cover some underflow problems.
     IF(vcover(i) < zlt(i)/10.0)THEN
        vcover(i) = vcover(i) * 10.0
     ENDIF
     ! soil thermal model layer depths (m)
     ztdep(i,1) = 6.0 - sodep(i)
     ztdep(i,2) = sodep(i) - rootd(i)
     ztdep(i,3) = 8. * (rootd(i)-0.02) / 15.
     ztdep(i,4) = 4. * (rootd(i)-0.02) / 15.
     ztdep(i,5) = 2. * (rootd(i)-0.02) / 15.
     ztdep(i,6) = 1. * (rootd(i)-0.02) / 15.
  ENDDO

  !Now need to call rada2 to obtain albedos.
  CALL rada2 (snow,zlt,z1,z2                         &
             ,asnow,tg,cosz,tice,ref,tran,chil       &
             ,green,vcover,soref,radfac,salb,valb,nalb,thermk  &
             ,tgeff4,tc,len)

  !Some initialization, copy of soil wetness
  DO i = 1,len
     tsnow(i) = MIN(tg(i),tice) ! Snow temp is lesser of ice or ground temp (K)
     cuprt(i) = cupr(i) * 0.001 ! converting units to m/sec from mm/sec
     lsprt(i) = lspr(i) * 0.001 ! converting units to m/sec from mm/sec
     zmelt(i) = 0.0 ! depth of melted water (m)
     roff(i) = 0.0  ! total runoff (surface and subsurface) (mm)
     !RAMS sends in water like soil, deepest layer indexed = 1.
     !RAMS soil_water is m3/m3 and includes porosity so we need to take
     ! that factor out for Sib to just have soil wetness fraction.
     wwwtem(i,1) = www(i,1)/poros(i) !soil wetness fraction - top level
     wwwtem(i,2) = www(i,2)/poros(i) !soil wetness fraction - mid level
     wwwtem(i,3) = www(i,3)/poros(i) !soil wetness fraction - lowest level
     if(wwwtem(i,1)>1.001 .or. wwwtem(i,2)>1.001 .or. wwwtem(i,3)>1.001) then
       print*,'wwwoversat',igrp+mi0(ngrid),jgrp+mj0(ngrid),wwwtem(i,1) &
        ,wwwtem(i,2),wwwtem(i,3),biome_f,poros(i),int(soiltype_f),soiltype
       stop
     endif
     !For permanently irrigated crop land (LEAF3 biotype), keep soil saturated
     if(INT(biome_f)==16)then
       !print*,'wwwcropsat',igrp+mi0(ngrid),jgrp+mj0(ngrid),wwwtem(i,1) &
       ! ,int(biome_f),biome,poros(i),int(soiltype_f),soiltype
       wwwtem(i,1) = 1.0
       wwwtem(i,2) = 1.0
       wwwtem(i,3) = 1.0
     endif
  ENDDO

  !First guesses for ta (=canopy airspace potential temperature (K))
  ! and ea (=canopy airspace water vapor pressure (hPa)) and
  ! em (=mixed layer water vapor pressure (hPa)) (see temrec 120)
  DO I=1,len
     THA(I) = TA(I) / BPS(I)
  ENDDO
  DO I=1,len
     EA(I) = SHA(I) * PS_X(I) / (0.622 + SHA(I))
     EM(I) = SH_X(I) * PS_X(I) / (0.622 + SH_X(I))
  ENDDO

  !Distribute incident radiation between canopy and surface
  CALL rnload (len,nsib,radvbc,radvdc,radnbc,radndc,dlwbot,VCOVER &
              ,thermk,radfac,radn,radc3)

  DO i = 1,len
     CANEX(i)  = 1.-(SNOW(i,1)*5.-Z1(i))/(Z2(i)-Z1(i))
     CANEX(i)  = MAX( 0.1 , CANEX(i) )
     CANEX(i)  = MIN( 1.0 , CANEX(i) )
     !d is dd corrected for snow covered canopy
     D(i)  = Z2(i) - ( Z2(i)-DD(i) ) * CANEX(i)
     Z0(i) = Z0D(i)/( Z2(i)-DD(i) ) * ( Z2(i)-D(i) )
     RBC(i)    = CC1(i)/CANEX(i)
     RDC(i)    = CC2(i)*CANEX(i)
     AREAS(i)    = MIN(1. , ASNOW*SNOW(i,2))
     SATCAP(i,1) = ZLT(i)*0.0001 * CANEX(i)
     !Collatz-Bounoua change satcap(2) to 0.0002
     ! SATCAP(i,2) = 0.002
     satcap(i,2) = 0.0002  ! lahouari
  ENDDO

  !Initialize energy and water budgets
  CALL balan (1, 1.0, zdepth, wwwtem, capac, cupr,  &
       lspr, roff, etmass, totwb, radt, chf, shf,   &
       dt, ect, eci, egs, egi, hc, hg, heaten,      &
       hflux, snow, thm, tc, tg, tgs, td,           &
       ps_x, kapa, len, ioffset, nsoil )

  !Calculation of flux potentials and constants prior to heat
  !flux calculations.
  CALL begtem (tc, tg, cp, hltm, ps_x, snomel            &
       ,zlt, clai, cww, wwwtem, poros, num_pi, psy       &
       ,phsat, bee, czc, czh, phc                        &
       ,tgs, etc, etg, btc, btg, rstfac                  &
       ,rsoil, hr,  wc, wg                               &
       ,snow, capac, areas, satcap, csoil, tice, grav    &
       ,len, nsib )

  !Now that we have the new psy (psychrometric 'constant'),
  !calculate the new CAS capacities.
  !Used to make this max(4.,z2(:)), but we might have to boost the
  !min value upwards (itb)
  cas_cap_heat(:) = ros(:) * cp * MAX(4.,z2(:))
  !I think cas_cap_vap should use cv instead of cp (itb).
  cas_cap_vap(:)  = ros(:) * cv * MAX(4.,z2(:)) / psy(:)
  cas_cap_co2(:)  =               MAX(4.,z2(:))  ! this goes out to phosib

  !CALCULATE Net radiation (RADT) USING RADIATION FROM PHYSICS AND CURRENT
  !LOSSES FROM CANOPY AND GROUND
  CALL NETRAD (radc3, radt, stefan, fac1, vcover, thermk, tc,   &
       tg, tice, dtc4, dtg4, dts4, closs, gloss, sloss,         &
       tgeff, areas, len )

  DO I=1,len
     THgeff(I) = Tgeff(I) / BPS(I)
  ENDDO

  !Compute saturation mixing ratio
  CALL vnqsat (1,tgeff,PS_X,SHgeff,len)

  !GET RESISTANCES FOR SIB
  !Calls routines to get carbon fluxes
  CALL VNTLAT (igrp,jgrp,grav, tice,                           &
       vkrmn, delta, dt, tc, tg, ts_x, ps_x, zlt,              &
       wwwtem, tgs, etc, etg,snow,                             &
       rstfac, rsoil, hr, wc, wg, snofac,                      &
       sh_x, z0, spdm, sha, ros,cas_cap_co2,                   &
       cu, ra, thvgm, rib, ustar, rstar,tstar,                 &
       ventmf, thm, tha, z2, d,                                &
       fc, fg, rbc, rdc,gect,geci,gegs,gegi,                   &
       respcp, rb, rd, rds,bps, rst, rc, ecmass,               &
       ea, hrr, assimn, bintc, ta, pco2m, po2m, vmax0,         &
       green, tran, ref,TimeVar%gmudmu, trop, trda, trdm, slti,&
       shti, hlti, hhti, radn, effcon, binter, gradm,          &
       atheta, btheta, aparkk, wsfws, wsfht, wsflt, wci,       &
       whs, omepot, assimpot, assimci, antemp, assimnp,        &
       wags, wegs, aparc, pfd, assim, td, wopt, zmx, wsat,     &
       soilscale, zmstscale, drst,                             &
       soilq10, ansqr,                                         &
       nsib, len, nsoil,                                       &
       thgeff, shgeff, ct, zwind, ztemp,                       &
       respg, respfactor, pco2ap, pco2i, pco2c, pco2s,         &
       co2cap,cflux)

  !Calculate partial derivatives of the various heat fluxes
  !with respect to ground/canopy/snow temp, as well as
  !some other derivatives.

  CALL DELLWF (DT,dtc4,dtg4,dts4,fac1,areas                &
       ,       lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc, len )

  CALL DELHF ( DT,CP,bps,ts_x,tgs,tsnow,tc,ta,ros,ra,rb,rd           &
       ,                HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA,HADTA    &
       ,                hc, hg, hs, fss, len)

  CALL DELEF ( DT,CP,ps_x,em,ea,ros,HRr,fc,fg                         &
       ,                ra,rb,rd,rc,rsoil,snow,capac,wwwtem           &
       ,                ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                 &
       ,                ESDEA,EADEA                                   &
       ,                ec,eg,es,fws,hltm,cas_cap_vap                 &
       ,                etc,etg                                       &
       ,                btc,btg                                       &
       ,                areas, gect,geci,gegs,gegi, psy, snofac, hr   &
       ,                len )

  !Get soil thermal properties
  CALL soilprop ( td, tgs, slamda, shcap, wwwtem, poros, ztdep, &
    asnow, snow(1,2), areas, tice, snomel,                      &
    nsib, len, nsoil )

  !This routine sets up the coupled system of partial differential
  !equations described in Sato et al.,
  CALL SIBSLV ( DT,GRAV2,CP,HLTM,tgs,tsnow                            &
       ,                 td(1,nsoil),slamda(1,nsoil),num_pi           &
       ,                 areas,fac1                                   &
       ,                 VENTMF,BPS,ros,psy                           &
       ,                 czh,czc,cas_cap_heat,cas_cap_vap             &
       ,                 lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc    &
       ,                 HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA          &
       ,                 HADTA                                        &
       ,                 hc, hg, hs, fss                              &
       ,                 ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                &
       ,                 ESDEA,EADEA                                  &
       ,                 ec,eg,es,fws                                 &
       ,                 etc,etg                                      &
       ,                 btc,btg                                      &
       ,                 RADT                                         &
       ,                 dtc, dtg, dta, dea                           &
       ,                 len)

  DO i = 1,len
     radt(i,2) = (1.-areas(i))*radt(i,2)+areas(i)*radt(i,3)
  ENDDO

  !UPDATING OF ALL HYDROLOGICAL PROGNOSTIC VARIABLES
  CALL UPDAT2 (snow ,capac , snofac, ect, eci, egi,              &
       egs, hltm, wwwtem, num_pi, czh, dtg, dtc, ta, dta,   &
       dea, dt,                                                  &
       roff, tc, td, tg, bee, poros, satco,                      &
       slope, phsat, zdepth, ecmass, egmass,                     &
       shf, tice, snomel, asnow, czc, csoil, chf,                &
       hc, hg, areas, q3l, q3o,                                  &
       qqq, zmelt1, cww, len, nsib, nsoil,                       &
       etc,ea,btc,                                               &
       geci,ros,cp,psy,gect,etg,btg,                             &
       gegs,hr,fg,gegi,rd,rb,hcdtc,hcdta,                        &
       hgdta,slamda(1,nsoil))

  DO i = 1, len
     EVT(i) = (55.56 / dt) * (ECMASS(i) + EGMASS(i))
     wegs(i) = wegs(i) * evt(i)
  ENDDO

  !Get soil temperature increments
  CALL soiltherm ( td, dtd, tgs, dtg, slamda, shcap, ztdep, dt, &
       nsib, len, nsoil )

  !Update prognostic variables, get total latent and sensible fluxes
  DO i = 1,len
     thmtem(i) = thm(i)
     shtem(i) = sh_x(i)
  ENDDO

  CALL addinc ( grav2, cp, dt, hc, hg,                      &
       ps_x, bps, ecmass,psy, ros, hltm, cas_cap_heat,      &
       cas_cap_vap,                                         &
       egmass, fss, fws, hflux, etmass,                     &
       td, thmtem, ts_x, shtem,                             &
       tc, tg, ta, ea, ra, em, sha,                         &
       dtd, dtc, dtg, dta, dea,                             &
       drst, rst, bintc, len, nsib, nsoil)

  !Compute saturation vapor pressure
  CALL vnqsat (2, ta, ps_x, eastar, len)

  !inter2 replaces interc
  CALL INTER2 ( cuprt, lsprt, snow, capac, wwwtem , num_pi, &
       satcap, cww, tc, tg, clai, zlt, chil, roff,           &
       snomel, zdepth, ts_x, tice, asnow, csoil,             &
       satco, dt, vcover, roffo, zmelt2 ,len, nsib, exo )

  !     Calculate the surface flux of CO2 due to SiB (tracer T18)

  !     The net flux to the atmosphere is given by the release of CO2
  !     by soil respiration minus the uptake of CO2 by photosynthesis

  !     ASSIMN is the net assimilation of CO2 by the plants (from SiB2)
  !     respFactor*soilScale is the rate of release of CO2 by the soil

  !     soilScale is a diagnostic of the instantaneous rate of
  !        soil respiration (derived by Jim Collatz, similar to TEM)
  !     respFactor is the annual total accumulation of carbon in the
  !        previous year at each grid cell (annual total ASSIMN)
  !        divided by the annual total of soilScale at the same grid pt.

  !     Surface flux of CO2 used to be merely Assimn-Respg. With the
  !     prognostic CAS, the calculation becomes
  !
  !     co2flux =  (CO2A - CO2M)/ra
  !
  !     with a temperature correction thrown in. This calculation is
  !     performed in phosib.

  !sfc CO2 source flux between CAS and ref level, also known as
  !net ecosystem exchange (nee) (mol/m^2/sec) is same as co2flx here
  DO i = 1,len
     co2flx(i) = cflux(i)
  ENDDO

  !Some quantities for diagnostic output
  DO i = 1,len
     !Canopy relative humidity
     rha(i) = ea(i)/eastar(i)
     !Total snow melt from updat2 and inter2
     zmelt(i) = zmelt1(i)+zmelt2(i)
     !Calculate an overall leaf conductance, which is QP2(162)
     !and is used as a weighting func in QP2(162 and 163)
     ggl(i) = 1. / (rst(i) * (rb(i) + rc(i)))
  ENDDO

  !Copy soil moisture back, remember to invert for RAMS
  !Set surface water layers to 1 if there is ground snow or liquid
  DO i=1,len
     www(i,3) = wwwtem(i,3)*poros(i) !volumetric soil moisture (m3/m3) - lowest
     www(i,2) = wwwtem(i,2)*poros(i) !volumetric soil moisture (m3/m3) - middle
     www(i,1) = wwwtem(i,1)*poros(i) !volumetric soil moisture (m3/m3) - top layer
     snow1(i) = snow(i,1)   !vegetation snow (kg/m^2)
     snow2(i) = snow(i,2)   !ground surface snow (kg/m^2)
     capac1(i) = capac(i,1) !vegetation liquid store (kg/m^2)
     capac2(i) = capac(i,2) !ground surface liquid interception store (kg/m^2)
     !surface water layers (#)
     !Sib assumes snow density is 0.25 that of water
     if(snow2(i)>0.0 .or. capac2(i)>0.0) then
        sfcwaterlev=1.
        sfcwatermas=snow2(i)+capac2(i)
        sfcwaterdep=snow2(i)/250. + capac2(i)/1000.
     else 
        sfcwaterlev=0.
        sfcwatermas=0.
        sfcwaterdep=0.
     endif
     if(snow1(i)>0.0 .or. capac1(i)>0.0) then !vegetation water layers exists
        vegwatermas=snow1(i)+capac1(i)
     else 
        vegwatermas=0.
     endif
  ENDDO

  !Send soil water and temperature fields back to RAMS indexing
  DO i=1,len
     !Set RAMS top soil temp to Sib ground surface temperature(K)
     tempk(mzg) = tg(i)
     !Set RAMS soil temp to Sib soil layer termpature(K)
     DO k=1,nsoil
        tempk(k) = td(i,k)
     ENDDO
     !Remember that SiB level-1 = RAMS level-7 = top soil level.
     !Fill in extra RAMS soil moisture layers with similar values.
     soil_water(1) = www(i,3) !required
     soil_water(2) = www(i,3)
     soil_water(3) = www(i,2)
     soil_water(4) = www(i,2)
     soil_water(5) = www(i,2) !required
     soil_water(6) = www(i,1)
     soil_water(7) = www(i,1) !required
  ENDDO

  !Update surface albedo and outgoing longwave radiation
  i=1
  uplwrf = stefan * (tgeff(i) ** 4.0)
  sfcswa = salb(i,2,1) + salb(i,2,2)

  !Output diagnostics and change units as needed.
  DO i=1,len
  vlt_out     = green(i)*zlt(i) !vegetation LAI
  vegalb_out  = valb(i,2,1) + valb(i,2,2) !vegetation albedo (fraction)
  vcover_out  = vcover(i)  !vegetation fractional area
  zlt_out     = zlt(i)     !total LAI
  veght_out   = z2(i)      !canopy top height (m)
  vrough_out  = z0d(i)     !surface(vegetation) roughness length (m)
  assimn_out  = assimn(i)*1.e6  !uptake of CO2 by canopy plants (umol/m^2/s)
  respg_out   = respg(i)*1.e6   !ground respiration flux (umol/m^2/sec)
  rstfac1_out = rstfac(i,1) !CANOPY RESISTANCE STRESS1 : leaf sfc humidity
  rstfac2_out = rstfac(i,2) !CANOPY RESISTANCE STRESS2 : soil moisture
  rstfac3_out = rstfac(i,3) !CANOPY RESISTANCE STRESS3 : temperature
  ect_out     = ect(i) / dt !transpiration flux (W/m^2)
  eci_out     = eci(i) / dt !canopy interception flux (W/m^2)
  egi_out     = egi(i) / dt !ground interception flux (W/m^2)
  egs_out     = egs(i) / dt !ground surface layer evap (W/m^2)
  hc_out      = hc(i)  / dt !canopy (veg) sensible heat flux (W/m^2)
  hg_out      = hg(i)  / dt !ground surface sensible heat flux (W/m^2)
  ra_out      = ra(i)      !CAS to RAMS aerodynamic resistance (sec/m)
  rb_out      = rb(i)      !leaf sfc to CAS aerodynamic resistance (sec/m)
  rc_out      = rc(i)      !total canopy resistance (sec/m)
  rd_out      = rd(i)      !ground to CAS aerodynamic resistance (sec/m)
  roff_out    = roff(i)    !runoff (surface and subsurface) (kg/m^2)
  green_out   = green(i)   !greenness fraction (-)
  apar_out    = aparc(i)   !absorbed fraction of PAR
  ventmf_out  = ventmf(i)  !ventilation mass flux (kg/m^2/sec)
  pco2c_out   = pco2c(i)   !leaf chloroplast CO2 concentration (Pa)
  pco2i_out   = pco2i(i)   !leaf internal CO2 concentration (Pa)
  pco2s_out   = pco2s(i)   !leaf surface CO2 concentration (Pa)
  pco2m_out   = pco2m(i)   !lowest atm level CO2 concentration (Pa)
  ea_out      = ea(i)      !CAS water vapor pressure (hPa)
  em_out      = em(i)      !ref level vapor pressure (hPa)
  rha_out     = rha(i)     !CAS relative humidity
  radvbc_out  = radvbc(i)  !radiation: visible beam (W/m^2)
  radvdc_out  = radvdc(i)  !radiation: visible diffuse (W/m^2)
  radnbc_out  = radnbc(i)  !radiation: nir beam (W/m^2)
  radndc_out  = radndc(i)  !radiation: nir diffuse (W/m^2)
  psy_out     = psy(i)     !psychrometric constant (hPa deg^-1)
  ENDDO

  !Print initial info to check for correctness
  sibprint=1
  i=1
  if(sibprint==1.and.ngrid==1.and.igrp+mi0(ngrid)==15.and.jgrp+mj0(ngrid)==15) then
   print*,'GRID POINT:',igrp+mi0(ngrid),jgrp+mj0(ngrid)
   print*,'SIB-Lat/doy/dt:',latitude,doy,dt
   print*,'SIB-Patcharea,NDVI:',patcharea,pndvi,cndvi
   print*,'SIB-types:',soiltype,biome_f,biome
   print*,'SIB-soildep123',ztdep(1,1),ztdep(1,2),ztdep(1,3)
   print*,'SIB-soildep456',ztdep(1,4),ztdep(1,5),ztdep(1,6)
   print*,'SIB-soilwetdep123',zdepth(1,1),zdepth(1,2),zdepth(1,3)
   print*,'SIB-sodep,rootd,vcover',sodep,rootd,vcover(i)  
   print*,'SIB-soilwater765',soil_water(7),soil_water(6),soil_water(5)
   print*,'SIB-soilwater432',soil_water(4),soil_water(3),soil_water(2)
   print*,'SIB-soilwater1',soil_water(1)
   print*,'SIB-PartPres:',pco2m_out,pco2ap
   print*,'SIB-respiration:',assimn_out,respg_out
   print*,'SIB-co2flx,sens,lat flx:',co2flx,fss,fws*hltm
   print*,'SIB-ta,tc,tg:',ta,tc,tg
   print*,'SIB-resist:',rstfac1_out,rstfac2_out,rstfac3_out
   print*,'SIB-stars:',ustar,tstar,rstar
   print*,'SIB-deltatemp',ts_x-ta,thm-tha
   print*,'SIB-sfcrad',uplwrf,sfcswa
   print*,'SIB-Vlai,Tlai',vlt_out,zlt_out
   print*,'SIB-Talbedo',salb(i,1,1),salb(i,1,2),salb(i,2,1),salb(i,2,2)
   print*,'SIB-Valbedo',valb(i,1,1),valb(i,1,2),valb(i,2,1),valb(i,2,2)
   print*,'SIB-Nalbedo',nalb(i,1,1),nalb(i,1,2),nalb(i,2,1),nalb(i,2,2)
   !print*,DT,GRAV2,CP,HLTM,tgs,tsnow
   !print*,td(1,nsoil),slamda(1,nsoil),num_pi
   !print*, areas,fac1			
   !print*, VENTMF,BPS,ros,psy			
   !print*, czh,czc,cas_cap_heat,cas_cap_vap	
   !print*, lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc
   !print*, HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA,HADTA
   !print*, hc, hg, hs, fss			
   !print*, ECDTC,ECDEA,EGDTG,EGDEA,ESDTS,ESDEA,EADEA		
   !print*, ec,eg,es,fws				
   !print*, etc,etg,btc,btg,RADT					
   !print*, dtc, dtg, dta, dea  	
  endif

return
END SUBROUTINE sib_2pt5

!##############################################################################
Subroutine raddrv (nsib,swdown,sunang,radvbc,radvdc,radnbc,radndc)

implicit none

  !---------------------------------------------------------------------
  !  radiation radive code to use the downward sw at bottom
  !  and the formulation to estimate radvbc,radvdc, radndc, radndc
  !  radvbc(len) : surface incident visible direct beam (W/m^2)
  !  radnbc(len) : surface incident near IR direct beam (W/m^2)
  !  radvdc(len) : surface incident visible diffuse beam (W/m^2)
  !  radndc(len) : surface incident near IR diffuse beam (W/m^2)
  !---------------------------------------------------------------------

  INTEGER nsib
  REAL swdown(nsib)
  REAL sunang(nsib), stemp
  REAL radvbc(nsib),radvdc(nsib)     &
       ,      radnbc(nsib),radndc(nsib),c1,c2,c3,c4,c5,cloud,difrat   &
       ,      vnrat

  INTEGER i

  C1 = 580.
  C2 = 464.
  C3 = 499.
  C4 = 963.
  C5 = 1160.

  DO i=1,nsib
     sunang(i) = MAX( 0.001 , sunang(i) )
     stemp = swdown(i)
     stemp = MAX(stemp,0.01 )
     cloud = (c5 * sunang(i) - stemp) / (c4 * sunang(i))
     cloud = MAX(cloud,0.)
     cloud = MIN(cloud,1.)
     !         cloud = max(0.58,cloud)

     !z  use the real clouds here!
     !         cloud = cldtot(i)
     !         CLOUD = AMAX1(CLOUD,0.)
     !         CLOUD = AMIN1(CLOUD,1.)

     DIFRAT = 0.0604 / ( SUNANG(i)-0.0223 ) + 0.0683
     IF ( DIFRAT .LT. 0. ) DIFRAT = 0.
     IF ( DIFRAT .GT. 1. ) DIFRAT = 1.

     DIFRAT = DIFRAT + ( 1. - DIFRAT ) * CLOUD
     VNRAT = ( C1 - CLOUD*C2 ) / ( ( C1 - CLOUD*C3 )    &
          + ( C1 - CLOUD*C2 ) )

     radvbc(i) = (1.-DIFRAT)*VNRAT*stemp
     radvdc(i) = DIFRAT*VNRAT*stemp
     radnbc(i) = (1.-DIFRAT)*(1.-VNRAT)*stemp
     radndc(i) = DIFRAT*(1.-VNRAT)*stemp
  ENDDO

return
END SUBROUTINE raddrv

!##############################################################################
Subroutine mapper (      &
     lat,                &
     DOY,                &
     prevNDVI,           &
     curNDVI,            &
     fVCover,            &
     ChiL,               &
     LTran,              &
     LRef,               &
     MorphTab,           &
     AeroVar,            &
     LAIgrid,            &
     fVCovergrid,        &
     TimeVar             &
     )

! calculates time dependant boundary condition variables for SiB.

implicit none

  ! begin input variables
  REAL lat         ! center latitude of grid cell
  REAL curNDVI     ! FASIR NDVI values for a grid cell
  REAL prevNDVI    ! previous month's NDVI value
  REAL fVCover     !
  REAL ChiL        !
  REAL LTran(2,2)  !
  REAL LRef(2,2)   !
  INTEGER DOY      ! Day of Year (DOY) of ndvi input map
  !
  ! begin input biome dependant, physical morphology variables
  TYPE biome_morph_var
     REAL zc        ! Canopy inflection height (m)
     REAL LWidth    ! Leaf width
     REAL LLength   ! Leaf length	
     REAL LAImax    ! Maximum LAI
     REAL stems     ! Stem area index
     REAL NDVImax   ! Maximum NDVI
     REAL NDVImin   ! Minimum NDVI
     REAL SRmax     ! Maximum simple ratio
     REAL SRmin     ! Minimum simple ratio
  END TYPE biome_morph_var
  TYPE(biome_morph_var) MorphTab
  !
  ! begin input aerodynamic parameters
  TYPE aero_var
     REAL zo       ! Canopy roughness coeff
     REAL zp_disp  ! Zero plane displacement
     REAL RbC      ! RB Coefficient
     REAL RdC      ! RC Coefficient
  END TYPE aero_var

  TYPE(aero_var),DIMENSION(50,50) :: AeroVar ! aerodynamic
  !  interpolation tables

  REAL LAIgrid(50)   ! grid of LAI values for lookup table
  REAL fVCovergrid(50)! grid of fVCover values for interpolation table
  !
  ! begin time dependant, output variables
  TYPE time_dep_var
     REAL fPAR    ! Canopy absorbed fraction of PAR
     REAL LAI     ! Leaf-area index
     REAL Green   ! Canopy greeness fraction of LAI
     REAL zo      ! Canopy roughness coeff
     REAL zp_disp ! Zero plane displacement
     REAL RbC     ! RB Coefficient (c1)
     REAL RdC     ! RC Coefficient (c2)
     REAL gmudmu  ! Time-mean leaf projection
  END TYPE time_dep_var
  TYPE(time_dep_var) TimeVar
  !
  ! begin internal variables
  REAL prevfPAR    ! previous month's fPAR value
  REAL, PARAMETER :: fPARmax=0.95
  !                   ! Maximum possible FPAR corresponding to 98th percentile
  REAL, PARAMETER :: fPARmin=0.01
  !                   ! Minimum possible FPAR corresponding to 2nd percentile
  !     For more information on fPARmin and fPARmax, see
  !     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)

  !-----------------------------------------------------------------------
  ! Calculate time dependent variables
  !-----------------------------------------------------------------------

  ! Calculate first guess fPAR
  ! use average of Simple Ratio (SR) and NDVI methods.
  !
  ! print*,'call avgapar:',prevndvi,morphtab%ndvimin,morphtab%ndvimax

  CALL AverageAPAR (prevNDVI,                       &
       MorphTab%NDVImin, MorphTab%NDVImax,            &
       MorphTab%SRmin, MorphTab%SRmax,                &
       fPARmax, fParmin, prevfPAR)

  CALL AverageAPAR (curNDVI,                        &
       MorphTab%NDVImin, MorphTab%NDVImax,            &
       MorphTab%SRmin, MorphTab%SRmax,                &
       fPARmax, fParmin, TimeVar%fPAR)
  !
  !
  ! Calculate leaf area index (LAI) and greeness fraction (Green)
  !   See S. Los et al 1998 section 4.2.
  !
  !   Select previous month

  CALL laigrn (TimeVar%fPAR,                &
       prevfPAR,                            &
       fPARmax,                             &
       fVCover,                             &
       MorphTab%stems,                      &
       MorphTab%LAImax,                     &
       TimeVar%Green,                       &
       TimeVar%LAI)

  ! Interpolate to calculate aerodynamic, time varying variables

  !  PRINT*,'call aeroint:',timevar%lai,fvcover

  CALL AeroInterpolate (                           &
       TimeVar%LAI,                         &
       fVCover,                             &
       LAIgrid,                             &
       fVCovergrid,                         &
       AeroVar,                             &
       TimeVar%zo,                          &
       TimeVar%zp_disp,                     &
       TimeVar%RbC,                         &
       TimeVar%RdC)

  ! Calculate mean leaf orientation to par flux (gmudmu)
  CALL gmuder (lat,                        &
       DOY,                                &
       ChiL,                               &
       TimeVar%gmudmu)

  ! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
  ! and greeness fraction, and LAI

  CALL aparnew (TimeVar%LAI,                &
       TimeVar%Green,                       &
       LTran,                               &
       LRef,                                &
       TimeVar%gmudmu,                      &
       fVCover,                             &
       TimeVar%fPAR,                        &
       fPARmax,                             &
       fPARmin)

return
END SUBROUTINE mapper

!##############################################################################
Subroutine laigrn (fPAR,fPARm,fPARmax,fVCover,stems,   &
     LAImax,Green,LAI)

! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR.
! LAI is linear with vegetation fraction and exponential with fPAR.
! See Sellers et al (1994), Equations 7 through 13.

implicit none

  ! begin input variables
  REAL fPAR     ! fraction of PAR absorbed by plants at current time
  REAL fPARm    ! fraction of PAR absorbed by plants at previous time
  REAL fPARmax  ! maximum possible FPAR corresponding to 98th percentile
  REAL fVCover  ! vegetation cover fraction
  REAL stems    ! stem area index for the specific biome type
  REAL LAImax   ! maximum total leaf area index for specific biome type

  ! begin output variables
  REAL Green    ! greeness fraction of the total leaf area index
  REAL LAI      ! area average total leaf area index

  ! begin internal variables
  REAL LAIg     ! green leaf area index at current time
  REAL LAIgm    ! green leaf area index at previous time
  REAL LAId     ! dead leaf area index at current time

  ! Calculate current and previous green leaf area index (LAIg and LAIgm):
  ! LAIg is log-linear with fPAR.  Since measured fPAR is an area average,
  ! divide by fVCover to get portion due to vegetation.  Since fVCover can
  ! be specified, check to assure that calculated fPAR does not exceed fPARMax.

  !  print*,'lai_1',fpar,fvcover,fparmax

  IF(fPAR/fVCover.GE.fPARmax) THEN
     LAIg=LAImax
  ELSE
     LAIg=alog(1.-fPAR/fVCover)*LAImax/alog(1-fPARmax)
  ENDIF

  IF(fPARm/fVCover.GE.fPARmax) THEN
     LAIgm=LAImax
  ELSE
     LAIgm=alog(1.-fPARm/fVCover)*LAImax/alog(1-fPARmax)
  ENDIF

  ! Calculate dead leaf area index (LAId):
  ! If LAIg is increasing or unchanged, the vegetation is in growth mode.
  ! LAId is then very small (very little dead matter).
  ! If LAIg is decreasing, the peak in vegetation growth has passed and
  ! leaves have begun to die off.  LAId is then half the change in LAIg,
  ! assuming half the dead leaves fall off.

  ! Growth mode dead leaf area index:
  IF (LAIg.GE.LAIgm) LAId=0.0001

  ! die-off (post peak growth) dead leaf area index:
  IF (LAIg.LT.LAIgm) LAId=0.5*(LAIgm-LAIg)

  ! Calculate area average, total leaf area index (LAI):
  LAI=(LAIg+LAId+stems)*fVCover

  !print*,'laigrn1',laig,laid,stems,fvcover

  ! Calculate greeness fraction (Green):
  ! Greeness fraction=(green leaf area index)/(total leaf area index)
  Green=LAIg/(LAIg+LAId+stems)

  !PRINT*,'end laigrn',LAI,Green,laimax

return
END SUBROUTINE laigrn

!##############################################################################
Subroutine AeroInterpolate (LAI,fVCover,LAIgrid,fVCovergrid,AeroVar,zo &
                           ,zp_disp,RbC,RdC)

! This routine calculates the aerodynamic parameters by bi-linear
! interpolation from a lookup table of previously calculated values.
! The interpolation table is a numpts x numpts LAI/fVCover grid with
! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.

implicit none

  ! begin input variables
  REAL LAI            ! actual area averaged LAI for interpolation
  REAL fVCover        ! vegetation cover fraction for interpolation
  REAL LAIgrid(50)    ! grid of LAI values for lookup table
  REAL fVCovergrid(50)! grid of fVCover values for interpolation table
  TYPE aero_var
     REAL zo      ! Canopy roughness coeff
     REAL zp_disp ! Zero plane displacement
     REAL RbC     ! RB Coefficient
     REAL RdC     ! RC Coefficient
  END TYPE aero_var
  TYPE(aero_var) AeroVar(50,50) ! interpolation tables

  ! begin output variables
  REAL RbC            ! interpolated Rb coefficient
  REAL RdC            ! interpolated Rd coefficient
  REAL zo             ! interpolated roughness length
  REAL zp_disp        ! interpolated zero plane displacement

  ! begin internal variables
  INTEGER i           ! index for LAI grid location
  INTEGER j           ! index for fVCover grid location
  REAL LocLAI         ! local LAI var. to prevent changing main LAI value
  REAL LocfVCover     ! local fVCover var. to prevent changing fVCover value
  REAL DLAI           ! grid spacing between LAI values in tables
  REAL DfVCover       ! grid spacing between fVCover values in tables

  !print*,'aerointerp:lai,fvc=',lai,fvcover

  ! Calculate difference in LAI and veg-cover between array values.
  ! Both go from about 0 to 1 in 50 array bins.
  ! Used for interpolating between lookup table values below.
  DLAI=LAIgrid(2)-LAIgrid(1)
  DfVCover=fVCovergrid(2)-fVCovergrid(1)

  ! Assign input LAI and fVCover to local variables and make sure
  ! they lie within the limits of the interpolation tables, assuring
  ! the LAI and fVCover values returned from the routine are not modified.
  LocLAI=MAX(LAI,0.02)
  LocfVCover=MAX(fVCover,0.01)

  ! determine the nearest array location for the desired LAI and fVCover
  i=INT(LocLAI/DLAI+1)
  j=INT(LocfVCover/DfVCover+1)
  j=MIN(j,49)

  ! interpolate RbC variable
  CALL interpolate (                                      &
       LAIgrid(i),                                       &
       LocLAI,                                           &
       DLAI,                                             &
       fVCovergrid(j),                                   &
       LocfVCover,                                       &
       DfVCover,                                         &
       AeroVar(i,j)%RbC,                                 &
       AeroVar(i+1,j)%RbC,                               &
       AeroVar(i,j+1)%RbC,                               &
       AeroVar(i+1,j+1)%RbC,                             &
       RbC)

  ! interpolate RdC variable
  CALL interpolate (                                      &
       LAIgrid(i),                                       &
       LocLAI,                                           &
       DLAI,                                             &
       fVCovergrid(j),                                   &
       LocfVCover,                                       &
       DfVCover,                                         &
       AeroVar(i,j)%RdC,                                 &
       AeroVar(i+1,j)%RdC,                               &
       AeroVar(i,j+1)%RdC,                               &
       AeroVar(i+1,j+1)%RdC,                             &
       RdC)

  ! interpolate roughness length
  CALL interpolate (                                      &
       LAIgrid(i),                                       &
       LocLAI,                                           &
       DLAI,                                             &
       fVCovergrid(j),                                   &
       LocfVCover,                                       &
       DfVCover,                                         &
       AeroVar(i,j)%zo,                                  &
       AeroVar(i+1,j)%zo,                                &
       AeroVar(i,j+1)%zo,                                &
       AeroVar(i+1,j+1)%zo,                              &
       zo)

  ! interpolate zero plane displacement
  CALL interpolate (                                      &
       LAIgrid(i),                                       &
       LocLAI,                                           &
       DLAI,                                             &
       fVCovergrid(j),                                   &
       LocfVCover,                                       &
       DfVCover,                                         &
       AeroVar(i,j)%zp_disp,                             &
       AeroVar(i+1,j)%zp_disp,                           &
       AeroVar(i,j+1)%zp_disp,                           &
       AeroVar(i+1,j+1)%zp_disp,                         &
       zp_disp)

return
END SUBROUTINE AeroInterpolate

!##############################################################################
Subroutine interpolate (x1,x,Dx,y1,y,Dy,z11,z21,z12,z22,z)

! calculates the value of z=f(x,y) by linearly interpolating
! between the 4 closest data points on a uniform grid.  The routine
! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the
! 4 closest data points (z11, z21, z12, and z22).

implicit none

  ! begin input variables
  REAL x1  ! the x grid location of z11
  REAL x   ! x-value at which you will interpolate z=f(x,y)
  REAL Dx  ! grid spacing in the x direction
  REAL y1  ! the y grid location of z11
  REAL y   ! y-value at which you will interpolate z=f(x,y)
  REAL Dy  ! grid spacing in the y direction
  REAL z11 ! f(x1, y1)
  REAL z21 ! f(x1+Dx, y1)
  REAL z12 ! f(x1, y1+Dy)
  REAL z22 ! f(x1+Dx, y1+Dy)

  ! begin output variables
  REAL z   ! f(x,y), the desired interpolated value

  ! begin internal variables
  REAL zp  ! z'=first interpolated value at (x, y1)
  REAL zpp ! z''=second interpolated value at (x, Y1+Dy)

  ! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
  zp=z11+(x-x1)*(z21-z11)/Dx

  ! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
  zpp=z12+(x-x1)*(z22-z12)/Dx

  ! interpolate between zp and zpp to calculate z at (x,y)
  z=zp+(y-y1)*(zpp-zp)/Dy

return
END SUBROUTINE interpolate

!##############################################################################
Subroutine AverageAPAR (ndvi,NDVImin,NDVImax,SRmin,SRmax,fPARmax,fParmin,fPAR)

! calculates Canopy absorbed fraction of Photosynthetically
! Active Radiation (fPAR) using an average of the Simple Ratio (sr)
! and NDVI methods (Los et al. (1999), eqn 5-6).  The empirical
! SR method assumes a linear relationship between fPAR and SR.
! The NDVI method assumes a linear relationship between fPAR and NDVI.

implicit none

  ! begin input variables
  REAL ndvi     ! normalized difference vegetation index
  REAL NDVImin  ! minimum NDVI for vegetation type
  REAL NDVImax  ! maximum NDVI for vegetation type
  REAL SRmin    ! minimum NDVI for vegetation type
  REAL SRmax    ! maximum NDVI for vegetation type
  REAL fPARmax  ! Maximum possible FPAR corresponding to 98th percentile
  REAL fPARmin  ! Minimum possible FPAR corresponding to 2nd percentile

  ! begin output variables
  REAL fPAR     ! Canopy absorbed fraction of PAR

  ! begin internal variables
  REAL LocNDVI  ! local value of NDVI to prevent changes in input value
  REAL sr       ! simple ratio of near IR and visible radiances
  REAL NDVIfPAR ! fPAR from NDVI method
  REAL SRfPAR   ! fPAR from SR method

  ! switch to local value of ndvi to prevent any changes going back to main
  LocNDVI=NDVI

  ! Insure calculated NDVI value falls within physical limits for veg. type
  LocNDVI=MAX(LocNDVI,NDVImin)
  LocNDVI=MIN(LocNDVI,NDVImax)
  !print*,'fpar1:',locndvi,ndvimin,ndvimax

  ! Calculate simple ratio (SR)
  sr=(1.+LocNDVI)/(1.-LocNDVI)

  ! Calculate fPAR using SR method (Los et al. (1999), eqn 5)
  SRfPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
  !print*,'fpar2:',sr,srmin,fparmax,fparmin,srmin,srmax,fparmin

  ! Calculate fPAR using NDVI method (Los et al. (1999), eqn 6)
  NDVIfPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/(NDVImax-NDVImin)+fPARmin

  ! take average of two methods
  fPAR=0.5*(SRfPAR+NDVIfPAR)
  !print*,'fpar3:',fpar,srfpar,ndvifpar

return
END SUBROUTINE AverageAPAR

!##############################################################################
Subroutine aparnew (LAI,Green,LTran,LRef,gmudmu,fVCover,fPAR,fPARmax,fPARmin)

! recomputes the Canopy absorbed fraction of Photosynthetically
! Active Radiation (fPAR), adjusting for solar zenith angle and the
! vegetation cover fraction (fVCover) using a modified form of Beer's law..
! See Sellers et al. Part II (1996), eqns. 9-13.

implicit none

  ! begin input variables
  REAL LAI       ! Leaf Area Index
  REAL Green     ! Greeness fraction of Leaf Area Index
  REAL LTran(2,2)! Leaf transmittance for green/brown plants
  REAL LRef(2,2) ! Leaf reflectance for green/brown plants
  !                      For LTran and LRef:
  !                        (1,1)=shortwave, green plants
  !                        (2,1)=longwave, green plants
  !                        (1,2)=shortwave, brown plants
  !                        (2,2)=longwave, brown plants
  REAL gmudmu    ! daily Time-mean canopy optical depth
  REAL fVCover   ! Canopy cover fraction
  REAL fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
  REAL fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile

  ! begin output variables
  REAL fPAR      ! area average Canopy absorbed fraction of PAR

  ! begin internal variables
  REAL scatp     ! Canopy transmittance + reflectance coefficient wrt PAR
  REAL PARk      ! mean canopy absorption optical depth wrt PAR

  ! Calculate canopy transmittance + reflectance coefficient wrt PAR
  ! transmittance + reflectance coefficient=green plants + brown plants
  scatp=Green*(LTran(1,1)+LRef(1,1))+        &
       (1.-Green)*(LTran(1,2)+LRef(1,2))

  ! Calculate PAR absorption optical depth in canopy adjusting for
  ! variance in projected leaf area wrt solar zenith angle
  ! (Sellers et al. Part II (1996), eqn. 13b)
  ! PAR absorption coefficient=(1-scatp)
  PARk=SQRT(1.-scatp)*gmudmu

  ! Calculate the new fPAR (Sellers et al. Part II (1996), eqn. 9)
  fPAR=fVCover*(1.-EXP(-PARk*LAI/fVCover))

  ! Ensure calculated fPAR falls within physical limits
  fPAR=amax1(fPARmin,fPAR)
  fPAR=amin1(fPARmax,fPAR)

return
END SUBROUTINE aparnew

!##############################################################################
Subroutine gmuder (Lat, DOY, ChiL, gmudmu)

! calculates daily time mean optical depth of canopy relative to the Sun.
! Look to see if this matches the definition of G(mu)/mu described in
! Bonan (1996) and Sellers (1985)

implicit none

  ! begin input variables
  REAL Lat      ! latitude in degrees
  INTEGER DOY   ! day-of-year (typically middle day of the month)
  REAL ChiL     ! leaf angle distribution factor

  ! begin output variables
  REAL gmudmu   ! daily time mean canopy optical depth relative to Sun
  REAL test     ! test variable

  ! begin internal variables
  REAL mumax    ! max cosine of the Solar zenith angle (noon)
  REAL mumin    ! min cosine of the Solar zenith angle (rise/set)
  REAL dec      ! declination of the Sun (Solar Declination)
  REAL pi180    ! conversion factor from degrees to radians
  REAL aa       ! minimum possible LAI projection vs. cosine Sun angle
  REAL bb       ! slope leaf area projection vs. cosine Sun angle

  ! Calculate conversion factor from degrees to radians
  pi180=3.14159/180.

  ! Calculate solar declination in degrees
  dec=23.5*SIN(1.72e-2*(DOY-80))

  ! Calculate maximum cosine of zenith angle corresponding to noon
  mumax=COS((dec-lat)*pi180)
  mumax=MAX(0.02, mumax)

  ! Assign min cosine zenith angle corresponding to start disc set (cos(89.4))
  mumin=0.01

  ! The projected leaf area relative to the Sun is G(mu)=aa+bb*mu
  ! Calculate minimum projected leaf area
  aa=0.5-0.633*ChiL-0.33*ChiL*ChiL

  ! Calculate slope of projected leaf area wrt cosine sun angle
  bb=0.877*(1.-2.*aa)

  ! Calculate mean optical depth of canopy by integrating G(mu)/mu over
  ! all values of mu.  Since G(mu) has an analytical form, this comes to
  gmudmu=aa*alog(mumax/mumin)/(mumax-mumin)+bb

return
END SUBROUTINE gmuder

!##############################################################################
Subroutine rada2 (snoww,zlt,z1,z2                                  &
                 ,asnow,tg,sunang,tf,ref,tran,chil                 &
                 ,green,vcover,soref,radfac,salb,valb,nalb         &
                 ,thermk,tgeff4                                    &
                 ,tc,len)

! CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT
! AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY

implicit none

!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!       SALB(2,2)      SURFACE ALBEDOS
!       VALB(2,2)      VEGETATION ALBEDOS
!       NALB(2,2)      NON-VEGETATION ALBEDOS
!       TGEFF4         EFFECTIVE SURFACE RADIATIVE TEMPERATURE (K)
!       RADFAC(2,2,2)  RADIATION ABSORPTION FACTORS
!       THERMK         CANOPY GAP FRACTION FOR TIR RADIATION
!++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!       ALBEDO(2,2,2)  COMPONENT REFLECTANCES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len
  REAL snoww(len,2),                                   &
       ref(len,2,2),tran(len,2,2),soref(len,2),        &
       salb(len,2,2),valb(len,2,2),nalb(len,2,2),      &
       radfac(len,2,2,2),asnow,                        &
       zlt(len),z1(len),z2(len),tg(len),sunang(len),   &
       chil(len),green(len),vcover(len),               &
       tgeff4(len),tc(len), thermk(len)
  REAL tf

  !     local variables
  REAL TRANC1(2), TRANC2(2), TRANC3(2), satcap(len,2),             &
       f(len),fmelt(len),zmew(len),albedo(len,2,2,2),              &
       canex(len), areas(len), facs, scov, scat, chiv, aa, bb,     &
       fac2, fac1, zkat, tg4, tc4, tgs, hh6,                       &
       hh5, hh3, hh2, zmk, hh10, hh9, hh8, hh7, den, zp, f1, ge,   &
       ek, epsi, power2, power1, zat, psi, hh4, hh1, fe, de, bot,  &
       ce, be, betao, upscat, acss, extkb, proj, tran2, tran1,     &
       reff1, reff2
  INTEGER iwave, i, irad

  !----------------------------------------------------------------------
  !     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO
  !         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING
  !         SNOW TRANSMITTANCE = 0.20, 0.54
  !-----------------------------------------------------------------------

  DO i = 1,len                      ! loop over gridpoint
     !        this portion is snow1 inlined
     CANEX(i)  = 1.-( SNOWw(i,1)*5.-Z1(i))/(Z2(i)-Z1(i))
     CANEX(i)  = MAX( 0.1, CANEX(i) )
     CANEX(i)  = MIN( 1.0, CANEX(i) )
     AREAS(i)    = MIN(1.0 , ASNOW*SNOWw(i,2))
     SATCAP(i,1) = ZLT(i)*0.0001 * CANEX(i)

     !Collatz-Bounoua change satcap(2) to 0.0002
     !          SATCAP(i,2) = 0.002
     SATCAP(i,2) = 0.0002             ! lahouari
     !    end old snow1
     F(i) = MAX(0.01746,SUNANG(i))
     FACS  = ( TG(i)-TF ) * 0.04
     FACS  = MAX( 0.0 , FACS)
     FACS  = MIN( 0.4, FACS)
     FMELT(i) = 1. - FACS
  ENDDO

  DO IWAVE = 1,2

     !----------------------------------------------------------------------
     DO i = 1,len                      ! loop over gridpoint
        SCOV =  MIN( 0.5, SNOWw(i,1)/SATCAP(i,1) )
        REFF1 = ( 1. - SCOV ) * REF(i,IWAVE,1) + SCOV * ( 1.2 -           &
             IWAVE * 0.4 ) * FMELT(i)
        REFF2 = ( 1. - SCOV ) * REF(i,IWAVE,2) + SCOV * ( 1.2 -           &
             IWAVE * 0.4 ) * FMELT(i)
        TRAN1 = TRAN(i,IWAVE,1) * ( 1. - SCOV )                           &
             + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT(i) )         &
             * TRAN(i,IWAVE,1)
        TRAN2 = TRAN(i,IWAVE,2) * ( 1. - SCOV )                           &
             + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT(i) ) * 0.9   &
             * TRAN(i,IWAVE,2)

        !----------------------------------------------------------------------
        !     CALCULATE AVERAGE SCATTERING COEFFICIENT, LEAF PROJECTION AND
        !     OTHER COEFFICIENTS FOR TWO-STREAM MODEL.
        !
        !      SCAT  (OMEGA)        : EQUATION (1,2) , SE-85
        !      PROJ  (G(MU))        : EQUATION (13)  , SE-85
        !      EXTKB (K, G(MU)/MU)  : EQUATION (1,2) , SE-85
        !      ZMEW  (INT(MU/G(MU)) : EQUATION (1,2) , SE-85
        !      ACSS  (A-S(MU))      : EQUATION (5)   , SE-85
        !      EXTK  (K, VARIOUS)   : EQUATION (13)  , SE-85
        !      UPSCAT(OMEGA*BETA)   : EQUATION (3)   , SE-85
        !      BETAO (BETA-0)       : EQUATION (4)   , SE-85
        !----------------------------------------------------------------------

        SCAT = GREEN(i)*( TRAN1 + REFF1 ) +( 1.-GREEN(i) ) *     &
             ( TRAN2 + REFF2)
        CHIV = CHIL(i)
        !
        IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01
        AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV
        BB = 0.877 * ( 1. - 2. * AA )

        PROJ = AA + BB * F(i)
        EXTKB = ( AA + BB * F(i) ) / F(i)
        ZMEW(i) = 1. / BB * ( 1. - AA / BB                        &
             * LOG ( ( AA + BB ) / AA ) )
        ACSS = SCAT / 2. * PROJ / ( PROJ + F(i) * BB )
        ACSS = ACSS * ( 1. - F(i) * AA                               &
             / ( PROJ + F(i) * BB ) * LOG ( ( PROJ                &
             +   F(i) * BB + F(i) * AA ) / ( F(i) * AA ) ) )

        UPSCAT = GREEN(i) * TRAN1 + ( 1.-GREEN(i) ) * TRAN2
        UPSCAT = 0.5 * ( SCAT + ( SCAT - 2. * UPSCAT ) *             &
             (( 1. - CHIV ) / 2. ) ** 2 )
        BETAO = ( 1. + ZMEW(i) * EXTKB )                             &
             / ( SCAT * ZMEW(i) * EXTKB ) * ACSS

        !----------------------------------------------------------------------
        !     Intermediate variables identified in appendix of SE-85.
        !
        !      BE          (B)     : APPENDIX      , SE-85
        !      CE          (C)     : APPENDIX      , SE-85
        !      BOT         (SIGMA) : APPENDIX      , SE-85
        !      HH1         (H1)    : APPENDIX      , SE-85
        !      HH2         (H2)    : APPENDIX      , SE-85
        !      HH3         (H3)    : APPENDIX      , SE-85
        !      HH4         (H4)    : APPENDIX      , SE-85
        !      HH5         (H5)    : APPENDIX      , SE-85
        !      HH6         (H6)    : APPENDIX      , SE-85
        !      HH7         (H7)    : APPENDIX      , SE-85
        !      HH8         (H8)    : APPENDIX      , SE-85
        !      HH9         (H9)    : APPENDIX      , SE-85
        !      HH10        (H10)   : APPENDIX      , SE-85
        !      PSI         (H)     : APPENDIX      , SE-85
        !      ZAT         (L-T)   : APPENDIX      , SE-85
        !      EPSI        (S1)    : APPENDIX      , SE-85
        !      EK          (S2)    : APPENDIX      , SE-85
        !----------------------------------------------------------------------

        BE = 1. - SCAT + UPSCAT
        CE = UPSCAT
        BOT = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
        IF ( ABS(BOT) .LE. 1.E-10) THEN
           SCAT = SCAT* 0.98
           BE = 1. - SCAT + UPSCAT
           BOT = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
        END IF
        DE = SCAT * ZMEW(i) * EXTKB * BETAO
        FE = SCAT * ZMEW(i) * EXTKB * ( 1. - BETAO )
        HH1 = -DE * BE + ZMEW(i) * DE * EXTKB - CE * FE
        HH4 = -BE * FE - ZMEW(i) * FE * EXTKB - CE * DE

        PSI = SQRT(BE**2 - CE**2)/ZMEW(i)

        ZAT = ZLT(i)/VCOVER(i)*CANEX(i)

        POWER1 = MIN( PSI*ZAT, 50.E0 )
        POWER2 = MIN( EXTKB*ZAT, 50.E0 )
        EPSI = EXP( - POWER1 )
        EK = EXP ( - POWER2 )

        ALBEDO(i,2,IWAVE,1) = SOREF(i,IWAVE)*(1.-AREAS(i))            &
             + ( 1.2-IWAVE*0.4 )*FMELT(i) * AREAS(i)
        ALBEDO(i,2,IWAVE,2) = SOREF(i,IWAVE)*(1.-AREAS(i))            &
             + ( 1.2-IWAVE*0.4 )*FMELT(i) * AREAS(i)
        GE = ALBEDO(i,2,IWAVE,1)/ALBEDO(i,2,IWAVE,2)

        !----------------------------------------------------------------------
        !     CALCULATION OF DIFFUSE ALBEDOS
        !
        !     ALBEDO(1,IWAVE,2) ( I-UP ) : APPENDIX , SE-85
        !----------------------------------------------------------------------

        F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
        ZP = ZMEW(i) * PSI

        DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -                    &
             ( BE - ZP ) * ( F1 + ZP ) * EPSI
        HH7 = CE * ( F1 - ZP ) / EPSI / DEN
        HH8 = -CE * ( F1 + ZP ) * EPSI / DEN
        F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
        DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI

        HH9 = ( F1 + ZP ) / EPSI / DEN
        HH10 = - ( F1 - ZP ) * EPSI / DEN
        TRANC2(IWAVE) = HH9 * EPSI + HH10 / EPSI

        ALBEDO(i,1,IWAVE,2) =  HH7 + HH8

        !----------------------------------------------------------------------
        !     CALCULATION OF DIRECT ALBEDOS AND CANOPY TRANSMITTANCES.
        !
        !     ALBEDO(1,IWAVE,1) ( I-UP )   : EQUATION(11)   , SE-85
        !     TRANC(IWAVE)      ( I-DOWN ) : EQUATION(10)   , SE-85
        !----------------------------------------------------------------------

        F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
        ZMK = ZMEW(i) * EXTKB

        DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -           &
             ( BE - ZP ) * ( F1 + ZP ) * EPSI
        HH2 = ( DE - HH1 / BOT * ( BE + ZMK ) )              &
             * ( F1 - ZP ) / EPSI -                        &
             ( BE - ZP ) * ( DE - CE*GE - HH1 / BOT      &
             * ( F1 + ZMK ) ) * EK
        HH2 = HH2 / DEN
        HH3 = ( BE + ZP ) * (DE - CE*GE -                    &
             HH1 / BOT * ( F1 + ZMK ))* EK -               &
             ( DE - HH1 / BOT * ( BE + ZMK ) ) *           &
             ( F1 + ZP ) * EPSI
        HH3 = HH3 / DEN
        F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
        DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI
        HH5 = - HH4 / BOT * ( F1 + ZP ) / EPSI -            &
             ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) +            &
             HH4 / BOT*( ZMK-F1 ) ) * EK
        HH5 = HH5 / DEN
        HH6 =   HH4 / BOT * ( F1 - ZP ) * EPSI +            &
             ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) +            &
             HH4 / BOT*( ZMK-F1 ) ) * EK
        HH6 = HH6 / DEN
        TRANC1(IWAVE) = EK
        TRANC3(IWAVE) = HH4 / BOT * EK + HH5 * EPSI + HH6 / EPSI

        ALBEDO(i,1,IWAVE,1) = HH1 / BOT + HH2 + HH3

        !----------------------------------------------------------------------
        !     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES
        !     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND
        !
        !      RADFAC      (F(IL,IMU,IV)) : EQUATION (19,20) , SE-86
        !----------------------------------------------------------------------

        RADFAC(i,2,IWAVE,1) = ( 1.-VCOVER(i) )                   &
             * ( 1.-ALBEDO(i,2,IWAVE,1) ) + VCOVER(i)            &
             * ( TRANC1(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,1) )      &
             + TRANC3(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
        !
        RADFAC(i,2,IWAVE,2) = ( 1.-VCOVER(i) )                    &
             * ( 1.-ALBEDO(i,2,IWAVE,2) ) + VCOVER(i)             &
             *  TRANC2(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) )
        !
        RADFAC(i,1,IWAVE,1) = VCOVER(i)                           &
             * ( ( 1.-ALBEDO(i,1,IWAVE,1) )                       &
             - TRANC1(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,1) )         &
             - TRANC3(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
        !
        RADFAC(i,1,IWAVE,2) = VCOVER(i)                           &
             * ( ( 1.-ALBEDO(i,1,IWAVE,2) )                     &
             - TRANC2(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
     ENDDO

     !----------------------------------------------------------------------
     !     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB ) WITH WEIGHTING
     !     FOR COVER FRACTIONS.
     !----------------------------------------------------------------------

     DO IRAD = 1,2
        DO i = 1,len        !  loop over gridpoint
           SALB(i,IWAVE,IRAD) = ( 1.-VCOVER(i) )               &
                * ALBEDO(i,2,IWAVE,IRAD) +             &
                VCOVER(i) * ALBEDO(i,1,IWAVE,IRAD)
           VALB(i,IWAVE,IRAD) = ALBEDO(i,1,IWAVE,IRAD)
           NALB(i,IWAVE,IRAD) = ALBEDO(i,2,IWAVE,IRAD)
        ENDDO
     ENDDO
     !
     !----------------------------------------------------------------------
     !
  ENDDO
  !
  !----------------------------------------------------------------------
  !
  !     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND
  !
  !----------------------------------------------------------------------
  !
  DO i = 1,len                  !  loop over gridpoint
     TGS = MIN(TF,TG(i))*AREAS(i)            &
          + TG(i)*(1.-AREAS(i))
     TC4 = TC(i)**4
     TG4 = TGS**4
     !
     ZKAT = 1./ZMEW(i) * ZLT(i) / VCOVER(i)
     ZKAT = MIN( 50.E0 , ZKAT )
     ZKAT = MAX( 1.E-5, ZKAT )
     THERMK(i) = EXP(-ZKAT)
     !
     FAC1 =  VCOVER(i) * ( 1.-THERMK(i) )
     FAC2 =  1.
     TGEFF4(i) =  FAC1 * TC4                         &
          + (1. - FAC1 ) * FAC2 * TG4
  ENDDO

return
END SUBROUTINE rada2

!##############################################################################
Subroutine rnload (len,nsib,radvbc,radvdc,radnbc,radndc,dlwbot,VCOVER   &
                  ,thermk,radfac,radn,radc3)

! calculation of absorption of radiation by surface.  Note that
! output from this calculation (radc3) only accounts for the
! absorption of incident longwave and shortwave fluxes.  The
! total net radiation calculation is performed in routine
! netrad.

implicit none

  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       RADN(2,3)      INCIDENT RADIATION FLUXES (W M-2)
  !       RADC3(2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len, nsib
  REAL  radc3(len,2), radn(len,2,2), radfac(nsib,2,2,2),       &
       radvbc(len), radvdc(len), radnbc(len), radndc(len),    &
       dlwbot(len), VCOVER(len), thermk(len)

  INTEGER i, iveg, iwave, irad

  !-----------------------------------------------------------------------
  !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
  !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
  !     SOURCE FOR TRANSPIRATION.
  !
  !      RADN        (F(IW,IMU,O)) : EQUATION (19-22) , SE-86
  !      RADC3       (FC,FGS)      : EQUATION (21,22) , SE-86
  !-----------------------------------------------------------------------

  DO i = 1,len
     RADc3(i,1) = 0.
     RADc3(i,2) = 0.
     radn(i,1,1) = radvbc(i)
     radn(i,1,2) = radvdc(i)
     radn(i,2,1) = radnbc(i)
     radn(i,2,2) = radndc(i)
  ENDDO

  DO iveg=1,2
     DO iwave=1,2
        DO irad=1,2
           DO i = 1,len
              radc3(i,iveg) = radc3(i,iveg) +         &
                   radfac(i,iveg,iwave,irad) *          &
                   radn(i,iwave,irad)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !     absorbed downwelling radiation only

  DO i = 1,len
     RADc3(i,1) = RADc3(i,1) + dlwbot(i) *            &
          VCOVER(i) * (1.- THERMK(i))
     RADc3(i,2) = RADc3(i,2) + dlwbot(i) *            &
          (1.-VCOVER(i) * (1.-THERMK(i)))
  ENDDO

  RETURN
END SUBROUTINE rnload

!##############################################################################
real Function E (X)

!E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNC OF TEMPERATURE
! Func E - Defined by Alvaro

implicit none

REAL, INTENT(IN) :: X

E = EXP( 21.18123 - 5418. / X ) / .622

return
END FUNCTION E

!##############################################################################
real Function GE (X)

!GE(X) IS D E(X) / D ( TEMP )
! Func GE - Defined by Alvaro

implicit none

REAL, INTENT(IN) :: X

GE = EXP( 21.18123 - 5418. / X ) * 5418. / (X*X) / .622

return
END FUNCTION GE

!##############################################################################
Subroutine begtem (tc,tg,cpair,hlat,psur,snomel     &
                 ,zlt,clai,cw,www,poros,pie,psy    &
                 ,phsat,bee,ccx,cg,phc             &
                 ,tgs,etc,etgs,getc,getgs,rstfac   &
                 ,rsoil,hr,wc,wg,snoww,capac       &
                 ,areas,satcap,csoil,tf,g          &
                 ,len,nlen)

implicit none

  !========================================================================
  !     Calculation of flux potentials and constants prior to heat
  !         flux calculations.  Corresponds to first half of TEMREC
  !         in 1D model.
  !========================================================================
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !       RSTFAC(2)      SOIL MOISTURE STRESS FACTOR
  !       RSOIL          SOIL SURFACE RESISTANCE (S M-1)
  !       HR             SOIL SURFACE RELATIVE HUMIDITY
  !       WC             CANOPY WETNESS FRACTION
  !       WG             GROUND WETNESS FRACTION
  !       CCX            CANOPY HEAT CAPACITY (J M-2 K-1)
  !       CG             GROUND HEAT CAPACITY (J M-2 K-1)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

real, external :: e
real, external :: ge

integer len, nlen
real :: www(len,3),rstfac(len,4),satcap(len,2),            &
        snoww(nlen,2),capac(nlen,2),                       &
        tc(nlen),tg(nlen),tm(nlen),hlat,psur(nlen),        &
        zlt(nlen),poros(nlen),phsat(nlen),                 &
        bee(nlen),ccx(len),cg(nlen),tgs(len),etc(len),     &
        etgs(len),getgs(len),rsoil(len),                   &
        hr(len),wc(len),wg(len),areas(len),csoil(len),     &
        getc(len), phc(nlen), psy(len)
real :: cpair, snomel, cw, tf, g, clai, pie
!     local variables
integer :: i
real :: slamda,shcap,x,tsnow,rsnow,fac,psit,argg   &
       ,phroot(len),shcapf,shcapu

  !----------------------------------------------------------------------
  !     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNC OF TEMPERATURE
  !     GE(X) IS D E(X) / D ( TEMP )
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !     CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.
  !     N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN
  !     DEALING WITH VERY LARGE SNOWPACKS.
  !     HEAT CAPACITY OF THE SOIL, AS USED IN FORCE-RESTORE HEAT FLUX
  !     DESCRIPTION. DEPENDENCE OF CSOIL ON POROSITY AND WETNESS.
  !
  !      CG          (cg)    : EQUATION (?) , CS-81
  !----------------------------------------------------------------------
  DO i = 1,len
        ! new calculation for ground heat capacity cg - no longer force-restore
        ! now for the top 1cm of snow and soil, with phase change in the soil
        ! incorporated into the heat capacity from +0.5C to -0.5C
        SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,1)*              &
             POROS(i)) * 4.186E6
        shcapf =  0.5 * (1. + poros(i) * (www(i,1)-1.0)) * 4.186E6
        IF(tg(i).GE.tf+0.5) THEN
           csoil(i) = shcapu * 0.02
        ELSE IF(tg(i).LE.tf-0.5) THEN
           csoil(i) = shcapf * 0.02
        ELSE
           csoil(i) = (0.5*(shcapu+shcapf) +                     &
                snomel*poros(i)*www(i,1) ) * 0.02
        ENDIF

        CCX(i) = ZLT(i)*CLAI+                               &
             (0.5*SNOWw(i,1)+CAPAC(i,1))*CW
        CG(i)  = (1.-areas(i))*CSOIL(i) +                    &
             cw * (capac(i,2) + 0.01 * areas(i))
  ENDDO

  DO i = 1,len
     ! HLAT(i)  = ( 3150.19 - 2.378 * TM(i) ) * 1000. !use constant passed in
     PSY(i) = CPAIR / HLAT * PSUR(i) / .622
     !----------------------------------------------------------------------
     !      Calculation of ground surface temperature and wetness fractions
     !----------------------------------------------------------------------
     TSNOW = MIN ( TF-0.01, TG(i) )
     RSNOW = SNOWw(i,2) /                             &
          (SNOWw(i,2)+CAPAC(i,2)+1.E-10)

     TGS(i) = TSNOW*AREAS(i) + TG(i)*(1.-AREAS(i))
     IF(tgs(i).LT.0.0)tgs(i) = SQRT(tgs(i))

     ETC(i)   = E(TC(i))
     ETGS(i)  = E(TGS(i))
     GETC(i)  = GE(TC(i))
     GETGS(i) = GE(TGS(i))

     WC(i) = MIN(1.0,(CAPAC(i,1)+SNOWw(i,1))/SATCAP(i,1))
     WG(i) = MAX(0.0,CAPAC(i,2)/SATCAP(i,2))*0.25

     !-----------------------------------------------------------------------
     !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
     !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
     !     SOURCE FOR TRANSPIRATION.
     !
     !      PHROOT      (PSI-R) : EQUATION (48) , SE-86
     !      RSTFAC(2)  F(PSI-L) : MODIFICATION OF EQUATION (12), SE-89
     !-----------------------------------------------------------------------

     PHROOT(i) = PHSAT(i) * MAX(0.02,WWW(i,2))**(-BEE(i))
     PHROOT(i) = MAX ( PHROOT(i),-2.E3)
     RSTFAC(i,2) = 1./( 1. + EXP( 0.02*( PHC(i)-PHROOT(i)) ))
     RSTFAC(i,2) = MAX(0.0001,RSTFAC(i,2))
     RSTFAC(i,2) = MIN(1.0,   RSTFAC(i,2))

     !----------------------------------------------------------------------
     !      RSOIL FUNC FROM FIT TO FIFE-87 DATA.  Soil surface layer
     !         relative humidity.
     !
     !      RSOIL      (RSOIL) : HEISER 1992 (PERSONAL COMMUNICATION)
     !      HR         (Fh)    : EQUATION (66) , SE-86
     !----------------------------------------------------------------------

     FAC = MIN( WWW(i,1), 1.0 )
     FAC = MAX( FAC,0.02)

     ! Collatz-Bounoua change rsoil to FIFE rsoil formulation from eq(19) SE-92
     !cbl     RSOIL(i) =  MAX (0.1, 694. - FAC*1500.) + 23.6
     rsoil(i) =  EXP(8.206 - 4.255 * fac)      ! lahouari

     PSIT = PHSAT(i) * FAC ** (- BEE(i) )
     ARGG = MAX(-10.0,(PSIT*G/ (461.5*TGS(i)) ))
     HR(i) = EXP(ARGG)
  ENDDO

return
END SUBROUTINE begtem

!##############################################################################
Subroutine NETRAD (radc3, radt, stefan, fac1,                &
                   vcover, thermk, tc, tg, tf,               &
                   dtc4, dtg4, dts4, closs, gloss, sloss,    &
                   tgeff, areas, len)

! CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
! LOSSES FROM CANOPY AND GROUND

implicit none

!=======================================================================
!pl bands in sib: 0.2 to 0.7 micromets are VIS, then
!pl               0.7 to 4.0 is NIR, above 4.0 it is thermal
!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!       RADt (2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len
  REAL  radc3(len,2),radt(len,3),stefan, fac1(len),         &
       vcover(len),thermk(len),tc(len),tg(len), tf,        &
       areas(len), tgeff(len),                             &
       dtc4(len), dtg4(len), dts4(len),                    &
       closs(len), gloss(len), sloss(len)

  !     local variables
  INTEGER i
  REAL tc4, tg4, ts4, tsnow, radtbar, zlwup              &
       ,      feedfac  &!  feedback factor
       ,      stefanx4  !  stefan-boltzmann const X 4.0
  DATA feedfac/1.0/

  stefanx4 = 4.0 * stefan

  DO I=1,len
     tsnow = MIN(tf, tg(i))
     !itb...Joerg Kaduk assures me this is quicker than an exponent...
     TC4 = TC(i) * TC(i) * TC(i) * TC(i)
     TG4 = TG(i) * TG(i) * TG(i) * TG(i)
     ts4 = tsnow * tsnow * tsnow * tsnow

     !pl effective ground cover for thermal radiation
     FAC1(i) =  VCOVER(i) * ( 1.-THERMK(i) )

     !pl derivatives
     DTC4(i) = stefanX4 * TC(i)* TC(i) * TC(i)
     DTG4(i) = stefanX4 * TG(i)* TG(i) * TG(i)
     Dts4(i) = stefanX4 * tsnow * tsnow * tsnow

     CLOSS(i) =  2. * FAC1(i) * STEFAN * TC4
     !pl canopy leaves thermal radiation loss
     CLOSS(i) =  CLOSS(i) - FAC1(i) * STEFAN *             &
          ( (1.-areas(i))*TG4+areas(i)*ts4)
     !pl ground loss
     GLOSS(i) =  STEFAN * TG4 - FAC1(i) * STEFAN * TC4
     !pl snow loss
     SLOSS(i) =  STEFAN * Ts4 - FAC1(i) * STEFAN * TC4
     !pl canopy leaves net radiation
     RADT(i,1) = RADC3(i,1) - closs(i)
     !pl ground net radiation
     RADT(i,2) = RADC3(i,2) - gloss(i)
     !pl snow net radiation
     RADT(i,3) = RADC3(i,2) - sloss(i)
     !pl bulk, weighted net radiation from combined ground and snow
     radtbar = areas(i)*radt(i,3) + (1.-areas(i))*radt(i,2)
     !pl this is the exchange meant to help out exchanges between
     !pl ground and snow

     radt(i,2) = radtbar + (1.+feedfac)*(radt(i,2)-radtbar)
     radt(i,3) = radtbar + (1.+feedfac)*(radt(i,3)-radtbar)
     !pl total thermal radiation up from surface
     zlwup = fac1(i) * tc4 +                                    &
          (1.-fac1(i)) * (areas(i)*ts4+(1.-areas(i))*tg4)
     !pl effective (combined) skin temperature from surface thermal radiation
     TGEFF(i) =  ZLWUP ** 0.25
  ENDDO

return
END SUBROUTINE NETRAD

!##############################################################################
Subroutine vnqsat (iflag, TQS, PQS, QSS, IM)

! Computes saturation mixing ratio or saturation vapour pressure
! as a func of temperature and pressure for water vapour.

implicit none

! INPUT VARIABLES     TQS,PQS,iflag,im
! OUTPUT VARIABLES    QSS
! ROUTINES called  (AMAX1,AMIN1)
! iflag = 1 for saturation mixing ratio, otherwise for saturation
!         vapor pressure

!Modifications:
! - removed routine VHQSAT and made the call to it in c3vint.F
! a call to VNQSAT adding the iflag=1 argument.  changed
! routine name from VQSAT to VNQSAT and all calls to
! VQSAT (c1subs.F, c3subs.F, comp3.F, hstatc.F, ocean.F) to
! VNQSAT with the iflag=1 argument.  kwitt  10/23/91

! converted to fortran 90 syntax - dd 6/17/93

  !     argument declarations

  INTEGER im,ic(IM), iflag
  REAL TQS(IM), PQS(IM), QSS(IM)

  !     local declarations

  REAL EST(139), tq1(IM), es1(IM), es2(IM), epsinv

  DATA EST/ 0.0031195, 0.0036135, 0.0041800, 0.0048227, 0.0055571,     &
       0.0063934, 0.0073433, 0.0084286, 0.0096407, 0.011014,             &
       0.012582, 0.014353, 0.016341, 0.018574, 0.021095, 0.023926,       &
       0.027096, 0.030652, 0.034629, 0.039073, 0.044028, 0.049546,       &
       0.055691, 0.062508, 0.070077, 0.078700, 0.088128, 0.098477,       &
       0.10983, 0.12233, 0.13608, 0.15121, 0.16784, 0.18615, 0.20627,    &
       0.22837, 0.25263, 0.27923, 0.30838, 0.34030, 0.37520, 0.41334,    &
       0.45497, 0.50037, 0.54984, 0.60369, 0.66225, 0.72589, 0.79497,    &
       0.86991, 0.95113, 1.0391, 1.1343, 1.2372, 1.3484, 1.4684,         &
       1.5979, 1.7375, 1.8879, 2.0499, 2.2241, 2.4113, 2.6126, 2.8286,   &
       3.0604, 3.3091, 3.5755, 3.8608, 4.1663, 4.4930, 4.8423,           &
       5.2155, 5.6140, 6.0394, 6.4930, 6.9767, 7.4919, 8.0406, 8.6246,   &
       9.2457, 9.9061, 10.608, 11.353, 12.144, 12.983, 13.873, 14.816,   &
       15.815, 16.872, 17.992, 19.176, 20.428, 21.750, 23.148, 24.623,   &
       26.180, 27.822, 29.553, 31.378, 33.300, 35.324, 37.454, 39.696,   &
       42.053, 44.531, 47.134, 49.869, 52.741, 55.754, 58.916, 62.232,   &
       65.708, 69.351, 73.168, 77.164, 81.348, 85.725, 90.305, 95.094,   &
       100.10, 105.33, 110.80, 116.50, 122.46, 128.68, 135.17, 141.93,   &
       148.99, 156.34, 164.00, 171.99, 180.30, 188.95, 197.96, 207.33,   &
       217.08, 227.22, 237.76, 248.71/

  EPSINV = 1./1.622

  tq1(:) = MAX(1.00001E0,(tqs(:)-198.99999))
  tq1(:) = MIN(138.900001E0,tq1(:))
  ic(:) = tq1(:)

  es1(:) = est(ic(:))
  es2(:) = est(ic(:)+1)

  qss(:) = ic(:)
  qss(:) = es1(:) + (es2(:)-es1(:)) * (tq1(:)-qss(:))
  tq1(:) = pqs(:) * epsinv
  qss(:) = MIN(tq1(:),qss(:))

  IF (IFLAG .EQ. 1) qss(:) = 0.622 * qss(:) / (pqs(:)-qss(:))

return
END SUBROUTINE vnqsat

!##############################################################################
Subroutine VNTLAT                                         &
     (igrp,jgrp,grav, tice,                               &
     vkrmn, delta, dtt, tc, gt, ts, ps, zlt,              &
     www, tgs, etc, etg,snoww,                            &
     rstfac, rsoil, hr, wc, wg,snofac,                    &
     sh, z0, spdm, sha, ros,cas_cap_co2,                  &
     cu, ra, thvgm, rib, ustar, rstar,tstar,              &
     ventmf, thm, tha, z2, d,                             &
     fc, fg, rbc, rdc,gect,geci,gegs,gegi,                &
     respcp, rb, rd, rds, bps, rst, rc, ecmass,           &
     ea, hrr, assimn, bintc, ta, pco2m, po2m, vmax0,      &
     green, tran, ref, gmudmu, trop, trda, trdm, slti,    &
     shti, hlti, hhti, radn, effcon, binter, gradm,       &
     atheta, btheta, aparkk, wsfws, wsfht, wsflt, wci,    &
     whs, omepot, assimpot, assimci, antemp, assimnp,     &
     wags, wegs, aparc, pfd, assim, td, wopt, zm, wsat,   &
     soilscale, zmstscale, drst,                          &
     soilq10, ansqr,                                      &
     nsib, len, nsoil,                                    &
     thgeff, shgeff, ct, zwind, ztemp,                    &
     respg, respfactor, pco2ap, pco2i, pco2c, pco2s       &
     ,           co2cap,cflux)

!Optimized routine phosib
!dd 92.06.10

implicit none

  !    argument list declarations
  INTEGER nsib, len, nsoil,igrp,jgrp
  REAL grav, tice                                                       &
       ,    vkrmn, delta, dtt, zwind, ztemp                             &
       ,    snofac
  REAL tc(len), gt(len), ts(len), ps(len), zlt(len)                     &
       ,    www(len,3)                                                  &
       ,    tgs(len), etc(len), etg(len)                                &
       ,    rstfac(len,4), rsoil(len), hr(len)                          &
       ,    wc(len), wg(len)                                            &
       ,    sh(len),cflux(len)                                          &
       ,    z0(len), spdm(len), sha(len), ros(len)                      &
       ,    cu(len), ra(len), thvgm(len), rib(len), ustar(len)          &
       ,    rstar(len), tstar(len)                                      &
       ,    ventmf(len), thm(len), tha(len), z2(len), d(len)            &
       ,    fc(len), fg(len), rbc(len), rdc(len),rds(len)               &
       ,    respcp(len), cogr(len), cogs(len), rb(len), rd(len)         &
       ,    bps(len), rst(len), rc(len), ecmass(len), ea(len)           &
       ,    hrr(len), assimn(len), bintc(len), ta(len)                  &
       ,     pco2m(len), po2m, vmax0(len), green(len), tran(len,2,2)    &
       ,     ref(len,2,2), gmudmu, trop(len), trda(len)                 &
       ,     trdm(len), slti(len), shti(len), hlti(len)                 &
       ,     hhti(len), radn(len,2,2), effcon(len), binter(len)         &
       ,     gradm(len), atheta(len), btheta(len)                       &
       ,     aparkk(len), wsfws(len), wsfht(len), wsflt(len)            &
       ,     wci(len), whs(len), omepot(len), assimpot(len)             &
       ,     assimci(len), antemp(len), assimnp(len), wags(len)         &
       ,     wegs(len), aparc(len), pfd(len), assim(len), td(nsib,5)    &
       ,     wopt(len), zm(len), wsat(len), soilscale(len,nsoil+1)      &
       ,     zmstscale(len,2)                                           &
       ,     drst(len), soilq10(len,nsoil+1), ansqr(len)                &
       ,     respg(len)                                                 &
       ,     thgeff(len), shgeff(len), ct(len)                          &
       ,     pco2s(len),pco2i(len)    & !added for neil suits' programs
       ,     pco2ap(len),pco2c(len)   &! more added nsuits vars
       ,     co2cap(len)   & ! moles of air in canopy
       ,     snoww(len,2) & ! snow cover (veg and ground) (m)
       ,     cas_cap_co2(len)         &
       ,     gect(len)                &
       ,     geci(len)                &
       ,     gegs(len)                &
       ,     gegi(len)

  REAL respfactor(nsib, nsoil+1)

  !    local variables

  INTEGER i
  REAL  zln2, ghalf, dmin, pdamp, qdamp, dttin,eps
  REAL u2(len), cog1(len), cog2(len), vsib(len)                 &
       ,    thsib(len), coc(len), tprcor(len)                        &
       ,    cni(len), cuni(len), ctni(len), ctni3(len), cti(len)     &
       ,    cui(len) , temv(len), cun(len), ctn(len)                 &
       ,    z1z0Urt(len), z1z0Trt(len), zzwind(len), zztemp(len)     &
       ,    epsc(len),epsg(len)


  EPS    = 1. / SNOFAC

  !Calculate damping factors
  zln2 = 6.9314718e-1
  ghalf = 1.0257068e1
  dttin = 3.6e3
  dmin = 6.0e1
  pdamp = EXP(-1.0 * zln2*(dtt*dmin)/(dttin*ghalf))
  qdamp = 1.0 - pdamp

  !Calculate stars and ventilation mass flux
  DO i = 1,len
     zzwind(i) = z2(i)-d(i)+ zwind
     zztemp(i) = z2(i)-d(i)+ ztemp
  ENDDO
  CALL VMFCALZ (VKRMN,DELTA,GRAV                       &
       ,                SH,z0                                  &
       ,                SPDM,SHA,ROS,CU,ct,THVGM,RIB           &
       ,                USTAR,RSTAR,TSTAR,VENTMF,THM           &
       ,                tha, zzwind, zztemp                    &
       ,                cuni, cun, ctn, z1z0Urt, z1z0Trt, len )

  !Calculate AERODYNAMIC RESISTANCE
  DO I=1,len
     RA(I)    = ROS(i) / VENTMF(i)
     TEMV(I) = (Z2(i) - D(i)) / z0(i)
     U2(I)     = SPDM(i) / (CUNI(i) * VKRMN)
     !print*,'temv',u2(i),temv(I),Z2(i),D(i),z0(i)
     temv(i) = LOG(temv(i))
     U2(I) = U2(I) * TEMV(I)
  ENDDO

  FC(:) = 1.0
  FG(:) = 1.0

  CALL RBRD (tc,rbc,zlt,z2,u2,rd,rb,ta,grav,rdc,tgs,len)

  DO I=1,len

     !itb...here is inserted some PL prog CAS stuff...
     epsc(i) = 1.
     epsg(i) = 1.
     !itb...pl says this only makes sense for canopy leaves, since
     !itb...there can only be water OR snow, not both. switching epsc
     !itb...epsc to eps makes the hltm adapt to freezing/fusion.

     IF(snoww(i,1) .GT. 0.0) epsc(i) = eps
     IF(snoww(i,2) .GT. 0.0) epsg(i) = eps

     RC(i) = RST(i) + RB(i) + RB(i)

     RDS(i) = RSOIL(i) * FG(i) + RD(i)

     GECT(i) =  (1. - WC(i)) / RC(i)
     GECI(i) = epsc(i) * WC(i) / (RB(i) + RB(i))

     GEGS(i) =  (1. - WG(i)) / RDS(i)
     GEGI(i) = epsg(i) * WG(i) / RD(i)

     COC(i) = GECT(i) + GECI(i)


     VSIB(I)  = 1.0/RB(I) + 1.0/RD(I)
     THSIB(I) = (GT(i)/RD(I) + TC(i)/RB(I)) /    &
          (BPS(i)*VSIB(I))
     ENDDO

     !czzggrst calculate ecmass -- canopy evapotranspiration

  DO i=1,len
     ecmass(i) = (etc(i) - ea(i)) * coc(i) *              &
          ros (i) * 0.622e0 /ps(i) * dtt
  ENDDO

  ! Calculate the rate of CO2 efflux from soils, according to the "R-star"
  ! approach of Denning et al (1996), adapted for use with the Bonan
  ! 6-layer soil thermodynamics module.
  ! Changed soil Q10 value for soil respiration from 2.0 to 2.4
  ! following Raich and Schelsinger (1992, Tellus 44B, 81-89),
  ! Scott Denning, 9/14/95
  CALL respsib (len, nsib, nsoil, wopt, zm, www, wsat,   &
       tgs, td, respfactor, respg, soilscale,  &
       zmstscale, soilq10)

  ! Calculation of canopy conductance and photosynthesis
  CALL phosib (igrp,jgrp,pco2m,pco2ap,po2m,vmax0,tice,ps,green           &
       ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda        &
       ,          trdm,slti,shti,hlti,hhti                               &
       ,          radn,etc                                               &
       ,          ea,rb,ra,ts                                            &
       ,          effcon,rstfac,binter,gradm,assimn                      &
       ,          rst,atheta                                             &
       ,          btheta,tgs,respcp,aparkk, len, nsib,                   &
       omepot,assimpot,assimci,antemp,assimnp,                      &
       wsfws,wsfht,wsflt,wci,whs,                                   &
       wags,wegs,aparc,pfd,assim, td,                               &
       soilscale,                                                   &
       drst, pdamp,qdamp,ecmass,dtt,bintc,tprcor,ansqr,             &
       nsoil, respg, pco2c, pco2i,                                  &
       pco2s,co2cap,cflux)

  DO i = 1,len
     bintc(i) = bintc(i) * tc(i) / ( 44.6 * tprcor(i))
     IF(ea(i).GT.etc(i)) fc(i) = 0.0
     IF(ea(i).GT.etg(i)) fg(i) = 0.0
     HRR(I) = HR(I)
     IF (FG(I) .LT. 0.5) HRR(I) = 1.0
  ENDDO

return
END SUBROUTINE VNTLAT

!##############################################################################
Subroutine VMFCALZ (VKRMN,DELTA,GRAV                         &
     ,                SH,z0                                      &
     ,                SPDM,SHA,ROS,CU,ct,THVGM,RIB               &
     ,                USTAR,RSTAR,TSTAR,VENTMF,THM, tha          &
     ,                zzwind, zztemp                             &
     ,                cuni, cun, ctn, z1z0Urt, z1z0Trt, len )

implicit none

  !****************************************************************************
  !    VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum
  !    and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
  !    and Boville(1993), and by Beljaars and Holtslag (1991).
  !
  !     Rerences:
  !       Beljars and Holtslag (1991): Flux parameterization over land surfaces
  !              for atmospheric models. J. Appl. Meteo., 30, 327-341.
  !       Holtslag and Boville (1993): Local versus nonlocal boundary-layer
  !              diffusion in a global climate model. J. of Climate, 6, 1825-
  !              1842.
  !       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
  !              atmosphere. Boundary-Layer Meteo., 17, 187-202.
  !       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
  !              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
  !              Layer parameterization, ECMWF, 59-79.
  !
  !     General formulation:
  !        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.)
  !     Transfer coefficients for mommentum and heat fluxes are:
  !        CU = CUN*Fm, and
  !        CT = CTN*Fh
  !        where  CUN and CTN are neutral values of momentum and heat transfers,
  !           and Fm and Fh are stability funcs derived from surface
  !           similarity relationships.
  !****************************************************************************

  INTEGER len
  REAL  vkrmn, delta, grav                             &
       ,    z0(len),CU(len)                                 &
       ,    ROS(len),SH(len)                                &
       ,    THM(len),THVGM(len),RIB(len),THGM(len)          &
       ,    SHA(len),SPDM(len),USTAR(len)                   &
       ,    RSTAR(len),TSTAR(len)                           &
       ,    CUI(len),CTI(len),CT(len), tha(len)             &
       ,    VENTMF(len), cuni(len)

  !     local variables
  REAL TEMV(len), wgm(len)                                           &
       ,     bunstablM, bunstablT, cunstablM, cunstablT, bstabl, cstabl   &
       ,    zzwind(len), zztemp(len), zrib(len), cun(len), ctn(len)       &
       ,      z1z0U(len), z1z0Urt(len), z1z0T(len), z1z0Trt(len)          &
       ,      fmomn(len),fheat(len),ribtemp,dm,dh

  INTEGER i

  !Set U* to a minimum value
  real, parameter :: ustmin = .1

  !  condtants for surface flux funcs, according to Holtslag and
  !      Boville (1993, J. Climate)

  bunstablM = 10.       ! constants for unstable func
  bunstablT = 15.
  cunstablM = 75.
  cunstablT = 75.

  DO I=1,len
     zrib(i) = zzwind(i) **2 / zztemp(i)
     WGM(i)  = SHA(i) - SH(i)

     !        SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE
     !        WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.

     THGM(i)  = THA(i)  - THM(i)
     THVGM(i) = THGM(i) + THA(i) * DELTA * WGM(i)
  ENDDO

  !   Ratio of reference height (zwind/ztemp) and roughness length:
  DO i = 1, len                !for all grid points
     z1z0U(i) = zzwind(i)/ z0(i)
     z1z0Urt(i) = SQRT( z1z0U(i) )
     z1z0U(i) = LOG( z1z0U(i) )
     z1z0T(i) = zzwind(i)/ z0(i)
     z1z0Trt(i) = SQRT( z1z0T(i) )
     z1z0T(i) = LOG( z1z0T(i) )
  ENDDO

  !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

  DO i = 1, len
     !cun and ctn = drag coefficients in neutral conditions, here same for h/m
     !cun and ctn = a2 from LEAF3-model
     cun(i) = VKRMN*VKRMN / (z1z0U(i)*z1z0U(i) )   !neutral Cm & Ct
     ctn(i) = VKRMN*VKRMN / (z1z0T(i)*z1z0T(i) )
     cuni(i) = z1z0u(i) / vkrmn
  ENDDO

  !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
  !   RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.

  DO I=1,len
     TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
     temv(i) = MAX(0.000001E0,temv(i))
     RIB(I) = -THVGM(I) * GRAV * zrib(i) / TEMV(i)
  ENDDO

  !   The stability funcs for momentum and heat/moisture fluxes as
  !   derived from the surface-similarity theory by Louis (1979, 1982), and
  !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag
  !   (1991).

  DO I=1,len
     IF(rib(i).GE.0.0) THEN
        !
        !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
        !
        rib(i) = MIN( rib(i), 0.5E0)
        fmomn(i)=1.+10.*rib(i)/(sqrt(1.+5.*rib(i)))
        fmomn(i)=1./fmomn(i)
        fmomn(i)=MAX(0.0001,fmomn(i))
        fheat(i)=1.+15.*rib(i)*(sqrt(1.+5.*rib(i)))
        fheat(i)=1./fheat(i)
        fheat(i)=MAX(0.0001,fheat(i))
     ELSE
        !
        !        THE UNSTABLE CASE.
        !
        ribtemp = ABS(rib(i))
        ribtemp = SQRT( ribtemp )
        dm = 1. + cunstablM * cun(i) * z1z0Urt(i) * ribtemp
        dh = 1. + cunstablT * ctn(i) * z1z0Trt(i) * ribtemp
        fmomn(i) = 1. - (bunstablM * rib(i) ) / dm
        fheat(i) = 1. - (bunstablT * rib(i) ) / dh
     END IF
  ENDDO

  !   surface-air transfer coefficients for momentum CU, for heat and
  !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

  DO i = 1, len
     CU(i) = CUN(i) * fmomn(i)
     CT(i) = CTN(i) * fheat(i)
     CUI(i) = 1. / CU(i)
     CTI(i) = 1. / CT(i)
  ENDDO

  !   Ustar, Tstar, Rstar and ventlation mass flux via Louis (1979,1981)
  DO I=1,len
     USTAR(i) = SPDM(i)*SPDM(i)*CU(i)
     USTAR(i) = SQRT( USTAR(i) )
     USTAR(i) = max(USTAR(i),ustmin)
     RSTAR(i) =  -WGM(i)*CT(i)*SPDM(i)/USTAR(i)
     TSTAR(i) = -THGM(i)*CT(i)*SPDM(i)/USTAR(i)
     VENTMF(i)= ROS(i)*CT(i)* SPDM(i)
  ENDDO

  !   Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION
  !   VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
  !   free convection conditions.

return
END SUBROUTINE VMFCALZ

!##############################################################################
Subroutine respsib (len, nsib, nsoil, wopt, zm, www, wsat, &
             tg, td, respfactor, respg, soilscale,   &
             zmstscale, soilq10)

implicit none

  INTEGER len, nsib, nsoil,i,l
  REAL wopt(len), zm(len), www(len,3), wsat(len),     &
       tg(len), td(nsib,nsoil),woptzm
  REAL respfactor(nsib,nsoil+1)
  ! output
  REAL respg(len),soilscale(len,nsoil+1), zmstscale(len,2),    &
       soilq10(len,nsoil+1),                                    &
       ! local arrays
       b(len,2)

  !    Calculates the rate of CO2 efflux from soils, according to the "R-star"
  !    approach of Denning et al (1996), adapted for use with the Bonan
  !    6-layer soil thermodynamics module

  !    Changed soil Q10 value for soil respiration from 2.0 to 2.4
  !    following Raich and Schelsinger (1992, Tellus 44B, 81-89),
  !    Scott Denning, 9/14/95

     DO i = 1,len

        ! Moisture effect on soil respiration, shallow and root zone
        woptzm = wopt(i)**zm(i)
        b(i,1) = (((100.*www(i,1))**zm(i)-woptzm)/                &
             (woptzm - 100.**zm(i)))**2
        b(i,2) = (((100.*www(i,2))**zm(i)-woptzm)/                &
             (woptzm - 100.**zm(i)))**2
        b(i,1) = MIN(b(i,1),10.E0)
        b(i,2) = MIN(b(i,2),10.E0)
        zmstscale(i,1) = 0.8*wsat(i)**b(i,1) + 0.2
        zmstscale(i,2) = 0.8*wsat(i)**b(i,2) + 0.2

        ! Temperature effect is a simple Q10, with a reference T of 25 C

        ! Deepest soil layers do not respire (no carbon below root zone)
        DO L = 1, 2
           soilscale(i,L) = 0.0
           soilq10(i,L) = 0.0
        ENDDO

        ! Layers 3 through nsoil (root zone) use WWW(2) and TD(3:nSoil)
        DO l = 3, nsoil
           soilQ10(i,L) = EXP(0.087547 * (td(i,L) - 298.15))
           soilscale(i,L) = soilQ10(i,L) * zmstscale(i,2)
        ENDDO

        ! Surface soil uses TG and water layer 1
        soilQ10(i,nsoil+1) = EXP(0.087547 * (tg(i) - 298.15))
        soilscale(i,nsoil+1) = soilQ10(i,nsoil+1) * zmstscale(i,1)

        ! Dimensionalize soil resp flux to balance annual budget.
        ! respFactor*soilScale is the rate of release of CO2 by the soil.
        ! respg is the ground respiration flux (umol/m^2/sec).
        respg(i) = respfactor(i,1) * soilscale(i,1)

        DO L = 2, nsoil+1
           respg(i) = respg(i) + respfactor(i,L) * soilscale(i,L)
        ENDDO

     ENDDO

return
END SUBROUTINE respsib

!##############################################################################
Subroutine phosib (igrp,jgrp,pco2m,pco2ap,po2m,vmax0,tf,psur,green  &
     ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda     &
     ,          trdm,slti,shti,hlti,hhti,radn,etc                   &
     ,          ea,rb,ra,tm                                         &
     ,          effcon,rstfac,binter,gradm,assimn                   &
     ,          rst,atheta,btheta,tgs,respcp                        &
     ,          aparkk,len,nsib,                                    &
     omepot,assimpot,assimci,antemp,assimnp,                   &
     wsfws,wsfht,wsflt,wci,whs,                                &
     wags,wegs,aparc,pfd,assim,td,                             &
     soilscale,                                                &
     drst,pdamp,qdamp,ecmass,dtt,bintc,tprcor,ansqr,           &
     nsoil, respg, pco2c, pco2i,                               &
     pco2s,co2cap,cflux)

implicit none

  !=======================================================================
  !
  !     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED
  !     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
  !     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
  !     BASE REFERENCE IS SE-92A
  !
  !                          UNITS
  !                         -------
  !
  !      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
  !      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
  !      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
  !      EFFCON                                   : MOL CO2 MOL QUANTA-1
  !      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
  !      EVAPKG                                   : KG M-2 S-1
  !      Q                                        : KG KG-1
  !
  !                       CONVERSIONS
  !                      -------------
  !
  !      1 MOL H2O           = 0.018 KG
  !      1 MOL CO2           = 0.044 KG
  !      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
  !      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
  !pl the next line applies to the Ci to Cs pathway
  !      GS  (CO2)           = GS (H2O) * 1./1.6
  !pl 44.6 is the number of moles of air per cubic meter
  !      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
  !      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
  !      MM  (MOLAIR/MOLH2O) = 1.611
  !
  !                         OUTPUT
  !                      -------------
  !
  !      ASSIMN              = CANOPY NET ASSIMILATION RATE
  !      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
  !      1/RST               = CANOPY CONDUCTANCE
  !      PCO2I               = INTERNAL CO2 CONCENTRATION
  !      RESPC               = CANOPY RESPIRATION
  !      RESPG               = GROUND RESPIRATION
  !
  !----------------------------------------------------------------------
  !
  !         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
  !         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
  !         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
  !         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
  !
  !-----------------------------------------------------------------------

  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1)
  !       RST            CANOPY RESISTANCE (S M-1)
  !       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS
  !
  !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
  !
  !       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
  !       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
  !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
  !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
  !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !     Modifications:
  !       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
  !         dd 950221
  !
  !      Modified for multitasking - introduced gather/scatter indices
  !          - DD 951206
  !
  !      Added in pco2c (chloroplast partial co2) for neil's fractionation
  !     calculations - Sep99

  !     input arrays:
  INTEGER len, nsib, nsoil,igrp,jgrp

  REAL vmax0(len),psur(len),green(len),gmudmu,              &
       zlt(len),cas_cap_co2(len),tc(len),ta(len),           &
       trop(len),trda(len),                                 &
       slti(len),shti(len),hlti(len),hhti(len),             &
       ra(len),rb(len),                                     &
       cog1(len),cog2(len),tm(len),effcon(len),             &
       binter(len),gradm(len),atheta(len),btheta(len),      &
       tgs(len),respcp(len),tran(len,2,2),ref(len,2,2),     &
       radn(len,2,2),ecmass(len),trdm(len),etc(len),        &
       aparc(len),rst(len), cflux(len)
  REAL pdamp, qdamp, dtt, pco2m(len), pco2ap(len), tf, po2m

  !     output arrays:
  REAL assimn(len),ea(len),rstfac(len,4),            &
       pco2i(len),respc(len),respg(len),drst(len)
  
  ! output arrays
  REAL omepot(len),assimpot(len),assimci(len),        &
       assimnp(len),whs(len),antemp(len),             &
       wsfws(len),wsfht(len),wsflt(len),wci(len),     &
       wags(len),wegs(len),pfd(len),                  &
       td(nsib,nsoil),assim(len),                     &
       soilscale(len,nsoil+1),                        &
       tprcor(len),bintc(len),                        &
       ansqr(len),                                    &
       pco2c(len),                                    & !chloroplast pco2
       xgah2o(len),                                   &
       xgco2m(len)

  !     work arrays:
  REAL PCO2Y(len,6), EYY(len,6),assimny(len,6),      &
       assimy(len,6)
  REAL c3(len),c4(len),RANGE(len),gammas(len),       &
       aparkk(len),gah2o(len),                       &
       gbh2o(len),                                   &
       par(len),rrkk(len),                           &
       omss(len),vm(len),gsh2o(len),pco2s(len),      &
       templ(len),temph(len),                        &
       qt(len),co2s(len),scatp(len),scatg(len),      &
       park(len),respn(len),zkc(len),                &
       zko(len),spfy(len), co2a(len), co2m(len),co2cap(len)

  INTEGER icconv(len),igath(len)

  REAL soilfrz(len)
  REAL  cwsflt, cwsfht, cwsfws,                                  &
       ccoms, ccomc, ascitemp, dompdomc, omsci, ompci, omcci,    &
       omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2,   &
       gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti,         &
       pco2in, pco2a, soilfrztd, soilfrztg

  INTEGER i, ic1, ic, l

  !pl introduce a co2 capacity
  !pl this will basically be the mass of air under the top of the canopy (in
  !pl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

  !pl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
  !pl carbon are in m2 s / mol air

  !pl one mole of gas occupies 22.4 cubic dm
  !pl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
  !pl the units for the carbon capacity are mol air /m2.
  !pl (e.g. here 893 moles if thickness of the layer is 20m)
  !pl this means that the units for pc02ap should be mol co2 / mol air, but
  !pl it is also possible to keep just co2 pressure and convert

  DO i = 1,len                !  LOOP OVER GRIDPOINT

     TPRCOR(i) = TF*PSUR(i)*100./1.013E5
     co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air / m2

     !pl this needs to be modified as in sibslv3 to automatically use the
     !pl thickness of the canopy air space.
     !----------------------------------------------------------------------
     !pl        RESPG(i) = 0. E -6 ! fixed respiration at 5 micromoles
     !pl   no longe rused since we now have respsib
     !----------------------------------------------------------------------

     IF( EFFCON(i) .GT. 0.07 ) THEN
        C3(i) = 1.
     ELSE
        C3(i) = 0.
     ENDIF
     C4(i)     = 1. - C3(i)

     !-----------------------------------------------------------------------
     !
     !     CALCULATION OF CANOPY PAR USE PARAMETER.
     !
     !      APARKK      (PI)     : EQUATION (31) , SE-92A
     !-----------------------------------------------------------------------

     SCATP(I) =     GREEN(i)   *              &
          ( TRAN(i,1,1) + REF(i,1,1) )    &
          +( 1.-GREEN(i) ) *              &
          ( TRAN(i,1,2) + REF(i,1,2) )
     SCATG(i) = TRAN(i,1,1) + REF(i,1,1)
     PARK(i) = SQRT(1.-SCATP(i)) * GMUDMU

     ! Collatz-Bounoua commented the calculation of APARC
     ! replaced it with theone calculated in new_mapper.
     ! APARC(i) = 1. - EXP ( -PARK(i)*ZLT(i) )   ! lahouari

     APARKK(i)   = APARC(i) / PARK(i) * GREEN(i)
     !-----------------------------------------------------------------------
     !
     !     Q-10 AND STRESS TEMPERATURE EFFECTS
     !
     !      QT          (QT)    : TABLE (2)     , SE-92A
     !-----------------------------------------------------------------------
     !
     qt(i) = 0.1*( TC(i) - TROP(i) )
     RESPN(i) = RESPCP(i) * VMAX0(i) * RSTFAC(i,2)

     !itb...patch to prevent underflow if temp is too cool...
     !Saleeby-2016-0222 (check to see if we need the "else" statement since later
     ! versions do not use this conditional statement.)
     IF(TC(i) >= TRDM(i))THEN
        RESPC(i) = RESPN(i) * 2.0**qt(i)            &
             /( 1. + EXP( TRDA(i)*(TC(i)-TRDM(i))))
     ELSE
        RESPC(i) = RESPN(i) * 2.0**qt(i)
     ENDIF

     VM(i) = VMAX0(i) * 2.1**qt(i)
     TEMPL(i) = 1. + EXP(SLTI(i)*(HLTI(i)-TC(i)))
     TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
     RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))

     VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i)    &
          + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)

     !-----------------------------------------------------------------------
     !
     !     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
     !     COMPENSATION POINT
     !
     !      ZKC          (KC)     : TABLE (2)     , SE-92A
     !      ZKO          (KO)     : TABLE (2)     , SE-92A
     !      SPFY         (S)      : TABLE (2)     , SE-92A
     !      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
     !      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
     !      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
     !-----------------------------------------------------------------------

     ZKC(i) = 30. * 2.1**qt(i)
     ZKO(i) = 30000. * 1.2**qt(i)
     SPFY(i) = 2600. * 0.57**qt(i)
     GAMMAS(i) = 0.5 * PO2M/SPFY(i) * C3(i)
     PFD(i)    = 4.6E-6 * GMUDMU *               &
          ( RADN(i,1,1)+RADN(i,1,2) )

     !pl these here all go from being m/s to being mol/ (m2 sec)
     GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
     GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
     GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)

     xgah2o(i) = MAX(0.466E0, gah2o(i) )
     xgco2m(i) = 4000.0 * vmax0(i)

     RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)        &
          + VMAX0(i)/5.* ( 1.8**qt(i)) * C4(i)
     PAR(i)    = pfd(i)*EFFCON(i)*( 1.-SCATG(i) )
     soilfrztg = 1.+EXP(-1.5 * (MAX(270.0E0,tgs(i))-273.16) )
     soilfrztd = 1.+EXP(-1.5 * (MAX(270.0E0,td (i,nsoil))-273.16) )
     soilfrz(i) = MAX(1./soilfrztg, 1./soilfrztd)
     soilfrz(i) = MAX( soilfrz(i), 0.05E0)
     BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*                &
          RSTFAC(i,2) * soilfrz(i)
!     print'(a,5g16.6)','bintc:',binter(i),zlt(i),green(i),   &
!                rstfac(i,2),soilfrz(i)
     OMSS(i)   = ( VMAX0(i)/2.0 ) * ( 1.8**qt(i) )         &
          /TEMPL(i) * RSTFAC(i,2) * C3(i)       &
          + RRKK(i) * RSTFAC(i,2) * C4(i)

     !-----------------------------------------------------------------------
     !     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
     !     ASSIMILATION RATE.
     !-----------------------------------------------------------------------

     RANGE(i)    = PCO2M(i) * ( 1. - 1.6/GRADM(i) ) - GAMMAS(i)
     icconv(i) = 1

  ENDDO

  DO IC = 1,6
     DO i = 1,len        ! LOOP OVER GRIDPOINT
        PCO2Y(i,IC) = 0.
        EYY(i,IC) = 0.
     ENDDO
  ENDDO

  !pl beginning of PL's setup

  DO i=1,len
     gah2o(i) =  1. / MAX(0.446E0,GAH2O(i))
  ENDDO

  !Bio...HERE IS THE ITERATION LOOP.
  !Bio...
  !Bio...We iterate on PCO2C-sortin makes a 'first guess' at
  !Bio...then orders PCO2C/Error pairs on increasing error size,
  !Bio...then uses a combination of linear/quadratic fit to obtain
  !Bio...the 'next best guess' as iteration count increases.
  !Bio...CYCALC uses that value of PCO2C to get the next value
  !Bio...of ASSIMN. CO2A and CO2S follow.
  !Saleeby-2016-0222 (newer Lixin Sib3.2 iterates on pco2c instead of pco2i
  ! and updates equation for pco2s
  DO IC = 1,6

     CALL SORTIN ( EYY, PCO2Y, RANGE, GAMMAS, ic,len )

     CALL CYCALC ( APARKK, VM, ATHETA, BTHETA,par,           &
          GAMMAS, RESPC, RRKK, OMSS, C3, C4,        &
          PCO2Y(1,ic), assimny(1,ic), assimy(1,ic), &
          len  )

     DO i = 1,len
        !pl first diagnose the current CO2 flux in mol / (m2 s)
        !pl this is a modified ra that will get us the right units
        !pl in the conservation equation. its units are m2 s / (mol_air)

        !pl now prognose the new CAS CO2 according to flux divergence
        !pl we are going to do this in mol C / mol air (same as PaC/PaAir)

        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)

        CO2A(i)   = (  CO2A(i) + (dtt/co2cap(i)) *         &
             (respg(i) - assimny(i,ic)              &
             +co2m(i)*gah2o(i)        ) )           &
             / (1+dtt*gah2o(i)/ co2cap(i) )

        pco2a = co2a(i) * psur(i) * 100.

        PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic)      &
             * PSUR(i)*100.)

        !Saleeby-2016-0222 (Not sure which of the two options below to use.
        !The first is from earlier Sib2.5 and the latter from Lixin Sib3.2
        PCO2IN = PCO2S(i) - ASSIMNy(i,ic) * PSUR(i) * 100.0 *  &
                                (1.6/GSH2O(i))
        !PCO2IN = PCO2S(i) - ASSIMNy(i,ic) * PSUR(i) * 100.0 *  &
        !                        (1.6/gsh2o(i) + 1.0/xgco2m(i))
   
        EYY(i,IC) = PCO2Y(i,IC) - PCO2IN
     ENDDO

     IF(ic.GE.2) THEN
        ic1 = ic-1
        DO i = 1,len        ! LOOP OVER GRIDPOINT
           IF(ABS(eyy(i,ic1)).GE.0.1)THEN
              icconv(i) = ic
           ELSE
              eyy(i,ic) = eyy(i,ic1)
              pco2y(i,ic) = pco2y(i,ic1)
           ENDIF
        ENDDO
     ENDIF

  ENDDO !iteration loop

  !pl end of PL's setup

  DO i = 1,len        ! LOOP OVER GRIDPOINT
     icconv(i) = MIN(icconv(i),6)
     igath(i) = i+(icconv(i)-1)*len
  ENDDO

  DO i = 1,len         ! LOOP OVER GRIDPOINT

     pco2i(i) = pco2y(igath(i),1)
     assimn(i) = assimny(igath(i),1)
     assim(i)  = assimy(igath(i),1)

     pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0

     !pl now do the real C_A forecast with the iterated fluxes.
     CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
     co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)

     CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) *          &
          (respg(i) - assimn(i)               &
          +co2m(i)*gah2o(i) ) )               &
          / (1+dtt*gah2o(i)/co2cap(i))

     !pl go back from molC / mol air to Pascals
     pco2ap(i) = co2a(i) * psur(i) * 100.

     !itb...carbon flux between CAS and reference level (MOL M-2 S-1)
     cflux(i) = gah2o(i)*(co2a(i)-co2m(i))*0.012
  ENDDO

  dtti = 1./dtt
  DO i = 1,len        ! LOOP OVER GRIDPOINT

     H2OI   = ETC(i)/PSUR(i)
     H2OA   =  EA(i)/PSUR(i)
     ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
     H2OS = H2OA + ECMOLE / GBH2O(i)
     H2OS  = MIN( H2OS, H2OI )
     H2OS  = MAX( H2OS, 1.E-7)
     H2OSRH = H2OS/H2OI

     !pl I have relaxed this condition to 1/10 of previous. The old way made
     !pl the CO2 on top of the leaves always at least 1/2 of the value at the
     !pl reference level.

     CO2S(i) = MAX(PCO2S(I),PCO2M(i)*0.05) / (PSUR(i)*100.)

     !pl Ball-Berry equation right here !

     GSH2OINF = (GRADM(i) * MAX(1.0E-12,ASSIMN(i))            &
          * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)

     !pl this is the change in stomatal resistance

     DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/          &
          (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))

     !Saleeby-2016-0222 (check this statement)
     !This is the 'would be change' if we did not use the damping factors
     !qdamp and pdamp.
     !rstnew = (1./gsh2oinf) * 44.6 * tprcor(i)/tc(i)
     !DRST(i) = rstnew - RST(i)

     RSTFAC(i,1) = H2OS/H2OI
     RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
  ENDDO

  !Z CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
  !-----------------------------------------------------------------------
  !  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
  !    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
  !    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
  !    PFD(i)
  !

  sttp2 = 73.**0.2
  ohtp2 = 100.**0.2
  DO i = 1,len
     !-----------------------------------------------------------------------
     ! CALCULATION OF POTENTIAL ASSIMILATION
     !-----------------------------------------------------------------------

     ! Make assimn a top leaf, not the canopy.
     ASSIMNp(i) = ASSIMN(i) / APARKK(i)

     ! Bottom stopped assim.
     ANTEMP(i) = MAX(0.E0,ASSIMNp(i))

     ! Potential intercellular co2.
     PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i)/         &
          ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))

     ! Potential rubisco limitation.
     OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/        &
          (PCO2IPOT+RRKK(i))*C3(i) + C4(i))

     ! Potential light limitation.
     OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/                 &
          (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))

     ! Quad 1.
     SQRTIN = MAX(0.E0,((OMEPOT(i)+OMCPOT)**2-                &
          4.*ATHETA(i)*OMEPOT(i)*OMCPOT))

     ! Quad 1. Intermediate  top leaf photosynthesis.
     OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))

     ! Potential sink or pep limitation.
     OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i)                  &
          + RRKK(i)*PCO2IPOT*C4(i)

     ! Quad 2.
     SQRTIN=MAX(0.E0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)*          &
          OMPPOT*OMSPOT))

     ! Quad 2. Final Potential top leaf photosynthesis.
     ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))

     !-----------------------------------------------------------------------
     ! CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
     !-----------------------------------------------------------------------

     ! Stressed rubisco limitation.
     OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i)    &
          + C4(i))

     ! Quad 1.
     SQRTIN = MAX(0.E0,(OMEPOT(i)+OMCCI)**2-          &
          4.*ATHETA(i)*OMEPOT(i)*OMCCI)

     ! Quad 1. Intermediate stress limited top leaf photosynthesis.
     OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.*ATHETA(i))

     ! Stressed sink or pep limitation.
     OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))

     ! Quad 2.
     SQRTIN = MAX(0.E0,(OMPCI+OMSCI)**2-4.*BTHETA(i)*OMPCI*OMSCI)

     ! Quad 2. Final stress limited top leaf photosynthesis.
     ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))

     !-----------------------------------------------------------------------
     ! CALCULATION OF CONTROL COEFFICIENTS
     !-----------------------------------------------------------------------

     ! Intermediate.
     DOMPDOMC = (OMPCI-OMEPOT(i))/                       &
          (2.*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))

     ! Bottom stopped final stress limited top leaf photosynthesis.
     ASCITEMP = MAX(ASSIMCI(i),1.0E-12)

     ! Rubisco control coefficient.
     CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/                        &
          (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP

     ! Sink or pep control coefficient.
     CCOMS = ((ASSIMCI(i)-OMPCI)/            &
          (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP

     !-----------------------------------------------------------------------
     !  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
     !-----------------------------------------------------------------------
     ! Canopy values (overwrites top leaf).

     OMEPOT(i) = OMEPOT(i)*APARKK(i)
     ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
     ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
     ASSIM(i) = ASSIM(i)*APARKK(i)
     ANTEMP(i) = ANTEMP(i)*APARKK(i)
     ANSQR(i) = ANTEMP(i)*ANTEMP(i)
     ASSIMNp(i) = ASSIMNp(i)*APARKK(i)

     !-----------------------------------------------------------------------
     ! OUTPUT:WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
     !-----------------------------------------------------------------------

     ! Water stress.
     WSFWS(i) = ASSIMPOT(i)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)

     ! High temperature stress.
     WSFHT(i) = ASSIMPOT(i)*(1.-1./TEMPH(i))*CCOMC

     ! Low temperature stress.
     WSFLT(i) = ASSIMPOT(i)*(1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))

     !  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
     cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
     IF(cwsfws.GT.1. .OR. cwsfws.LT.0.) wsfws(i)=0.
     cwsfht = (1.-1./TEMPH(i))*CCOMC
     IF(cwsfht.GT.1. .OR. cwsfht.LT.0.) wsfht(i)=0.
     cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
     IF(cwsflt.GT.1. .OR. cwsflt.LT.0.) wsflt(i)=0.

     ! Intermediate assimilation weighted Ci.
     WCI(i) = ANTEMP(i)*PCO2I(i)

     ! Intermediate assimilation weighted relative humidty stress factor.
     WHS(i) = ANTEMP(i)*RSTFAC(i,1)

     ! Intermediate assimilation weighted stomatal conductance.
     WAGS(i) = GSH2O(i)*ANTEMP(i)

     ! Intermediate evaporation weighted stomatal conductance.(Step 1.
     !   Step 2 after routine updat2)
     WEGS(i) = GSH2O(i)

  ENDDO

return
END SUBROUTINE phosib

!##############################################################################
Subroutine SORTIN (EYY,PCO2Y,RANGE,GAMMAS,IC,len)

  !=======================================================================
  !
  !     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
  !       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
  !       QUADRATIC FITS.
  !
  !=======================================================================

implicit none

  INTEGER len, ic
  REAL EYY(len,6),PCO2Y(len,6),RANGE(len),gammas(len)

  !     work arrays

  REAL eyyi1(len),eyyi2(len),eyyi3(len),eyyis(len),    &
       eyyisp(len),pco2yis(len),pco2yisp(len),         &
       pco2b(len),                                     &
       pco2yi1(len),pco2yi2(len),pco2yi3(len)
  REAL aterm, bterm, cterm, pco2yq, cc1, cc2, bc1, bc2, ac1, ac2,  &
       pco2yl, a, b, pmin, emin
  INTEGER is(len)
  LOGICAL bitx(len)
  INTEGER i, ix, i1, i2, i3, isp, n, l, j

  IF( IC .LT. 4 ) THEN
     DO i = 1,len
        PCO2Y(i,1) = GAMMAS(i) + 0.5*RANGE(i)
        PCO2Y(i,2) = GAMMAS(i)                                    &
             + RANGE(i)*( 0.5 - 0.3*SIGN(1.0,EYY(i,1)) )
        PCO2Y(i,3) = PCO2Y(i,1)- (PCO2Y(i,1)-PCO2Y(i,2))          &
             /(EYY(i,1)-EYY(i,2)+1.E-10)*EYY(i,1)
        !
        PMIN = MIN( PCO2Y(i,1), PCO2Y(i,2) )
        EMIN = MIN(   EYY(i,1),   EYY(i,2) )
        IF ( EMIN .GT. 0. .AND. PCO2Y(i,3) .GT. PMIN )            &
             PCO2Y(i,3) = GAMMAS(i)
     ENDDO
  ELSE

     N = IC - 1
     DO l = 1,len
        bitx(l) = ABS(eyy(l,n)).GT.0.1
        IF(.NOT.bitx(l)) pco2y(l,ic) = pco2y(l,n)
     ENDDO
     DO l = 1,len
        IF(bitx(l)) THEN
           DO J = 2, N
              A = EYY(l,J)
              B = PCO2Y(l,J)
              DO I = J-1,1,-1
                 IF(EYY(l,I) .LE. A ) GO TO 100
                 EYY(l,I+1) = EYY(l,I)
                 PCO2Y(l,I+1) = PCO2Y(l,I)
              ENDDO
              i = 0
100           CONTINUE
              EYY(l,I+1) = A
              PCO2Y(l,I+1) = B
           ENDDO
        ENDIF
     ENDDO

     !-----------------------------------------------------------------------

     DO l = 1,len
        IF(bitx(l)) THEN
           PCO2B(l) = 0.
           IS(l)    = 1
        ENDIF
     ENDDO

     DO IX = 1, N
        DO l = 1,len
           IF(bitx(l)) THEN
              IF( EYY(l,IX) .LT. 0. )  THEN
                 PCO2B(l) = PCO2Y(l,IX)
                 IS(l) = IX
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     DO l = 1,len
        IF(bitx(l)) THEN
           I1 = IS(l)-1
           I1 = MAX(1, I1)
           I1 = MIN(N-2, I1)
           I2 = I1 + 1
           I3 = I1 + 2
           ISP   = IS(l) + 1
           ISP = MIN0( ISP, N )
           IS(l) = ISP - 1
           eyyisp(l) = eyy(l,isp)
           eyyis(l) = eyy(l,is(l))
           eyyi1(l) = eyy(l,i1)
           eyyi2(l) = eyy(l,i2)
           eyyi3(l) = eyy(l,i3)
           pco2yisp(l) = pco2y(l,isp)
           pco2yis(l) = pco2y(l,is(l))
           pco2yi1(l) = pco2y(l,i1)
           pco2yi2(l) = pco2y(l,i2)
           pco2yi3(l) = pco2y(l,i3)
        ENDIF
     ENDDO

     DO l = 1,len
        IF(bitx(l)) THEN

           !itb...Neil Suits' patch to check for zero in the denominator...
           IF(EYYis(l) /= EYYisp(l))THEN
              PCO2YL=PCO2Yis(l)                &
                   - (PCO2Yis(l)-PCO2Yisp(l))    &
                   /(EYYis(l)-EYYisp(l))*EYYis(l)
           ELSE
              PCO2YL = PCO2Yis(l) * 1.01
           ENDIF
           !
           !   METHOD USING A QUADRATIC FIT
           !
           AC1 = EYYi1(l)*EYYi1(l) - EYYi2(l)*EYYi2(l)
           AC2 = EYYi2(l)*EYYi2(l) - EYYi3(l)*EYYi3(l)
           BC1 = EYYi1(l) - EYYi2(l)
           BC2 = EYYi2(l) - EYYi3(l)
           CC1 = PCO2Yi1(l) - PCO2Yi2(l)
           CC2 = PCO2Yi2(l) - PCO2Yi3(l)

           !itb...Neil Suits' patch to prevent zero in denominator...
           IF(BC1*AC2-AC1*BC2 /= 0.0 .AND. AC1 /= 0.0)THEN
              BTERM = (CC1*AC2-CC2*AC1)/(BC1*AC2-AC1*BC2)
              ATERM = (CC1-BC1*BTERM)/AC1
              CTERM = PCO2Yi2(l)                              &
                   -ATERM*EYYi2(l)*EYYi2(l)-BTERM*EYYi2(l)
              PCO2YQ= CTERM
              PCO2YQ= MAX( PCO2YQ, PCO2B(l) )
              PCO2Y(l,IC) = ( PCO2YL+PCO2YQ)/2.
           ELSE
              PCO2Y(l,IC) = PCO2Y(l,IC) * 1.01
           ENDIF

        ENDIF
     ENDDO

  ENDIF
  DO i = 1,len
     pco2y(i,ic) = MAX(pco2y(i,ic),0.01E0)
  ENDDO

return
END SUBROUTINE SORTIN

!##############################################################################
Subroutine CYCALC ( APARKK, VM, ATHETA, BTHETA, par,     &
     GAMMAS, RESPC, RRKK, OMSS, C3, C4,   &
     PCO2I, ASSIMN, assim, len )

  !=======================================================================
  !
  !     CALCULATION EQUIVALENT TO STEPS IN FIGURE 4 OF SE-92A
  !     C4 CALCULATION BASED ON CO-92.
  !
  !=======================================================================

implicit none

  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
  !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
  !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
  !
  !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
  !
  !       OMC            RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1)
  !       OME            LIGHT LIMITED ASSIMILATION (MOL M-2 S-1)
  !       OMS            SINK LIMITED ASSIMILATION (MOL M-2 S-1)
  !       CO2S           CANOPY SURFACE CO2 CONCENTRATION (MOL MOL-1)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len
  REAL aparkk(len),vm(len),atheta(len),    &
       btheta(len),gammas(len),par(len),   &
       respc(len),rrkk(len),omss(len),     &
       c3(len),c4(len),pco2i(len),         &
       assimn(len),assim(len)

  !    local variables
  REAL ome, omc, omp, oms, sqrtin
  INTEGER i

  !-----------------------------------------------------------------------
  !     CALCULATE ASSIMILATION RATE
  !
  !      OMC         (OMEGA-C): EQUATION (11) , SE-92A
  !      OME         (OMEGA-E): EQUATION (12) , SE-92A
  !      OMS         (OMEGA-S): EQUATION (13) , SE-92A
  !      ASSIMN      (A-N)    : EQUATION (14,15), SE-92A
  !-----------------------------------------------------------------------

  DO i = 1,len
     OMC = VM(i) *(PCO2I(i)-GAMMAS(i))/(PCO2I(i) + RRKK(i))*C3(i)   &
          + VM(i) * C4(i)
     OME = PAR(i)*(PCO2I(i)-GAMMAS(i))/(PCO2I(i)+2.*GAMMAS(i))*C3(i) &
          + PAR(i) * C4(i)
     SQRTIN= MAX( 0.E0, ( (OME+OMC)**2 - 4.*ATHETA(i)*OME*OMC ) )
     OMP  = ( ( OME+OMC ) - SQRT( SQRTIN ) ) / ( 2.*ATHETA(i) )
     OMS  = OMSS(i) * C3(i) + OMSS(i)*PCO2I(i) * C4(i)
     SQRTIN= MAX( 0.E0, ( (OMP+OMS)**2 - 4.*BTHETA(i)*OMP*OMS ) )
     ASSIM(i) = ( ( OMS+OMP ) - SQRT( SQRTIN ) ) /    &
          ( 2.*BTHETA(i) )
     ASSIMN(i)= ( ASSIM(i) - RESPC(i)) * APARKK(i)

  ENDDO

return
END SUBROUTINE CYCALC

!##############################################################################
Subroutine DELLWF (DTA,dtc4,dtg4,dts4,fac1,areas    &
     ,                   lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc   &
     ,                   len )

  ! Calculation of partial derivatives of canopy and ground radiative
  ! heat fluxes with respect to Tc, Tgs.
  ! Here we are doing only the long wave radiative loss, which is the
  ! only radiative quantity we are trying to bring to the next time step.

  !------------------------------INPUT is coming from Netrad-------------
  !
  !       dtc4, dtg4, dts4, which are the derivatives of the LW loss
  !
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       LCDTC          dLC/dTC
  !       LCDTG          dLC/dTG
  !       LCDTS          dLC/dTS
  !       LGDTG          dLG/dTG
  !       LGDTC          dLG/dTC
  !       LSDTS          dLS/dTS
  !       LSDTC          dLS/dTC
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

  INTEGER len
  REAL                                                          &
       lcdtc(len),lcdtg(len),lcdts(len),lgdtg(len),lgdtc(len)   &
       ,    lsdts(len),lsdtc(len),areas(len),fac1(len)               &
       ,    dtc4(len),dtg4(len),dts4(len),dta
  !     local variables
  INTEGER i

  DO I=1,len

     !pl canopy leaves:
     LCDTC(I) =   2 * DTC4(i) * fac1(i)
     LCDTG(I) =     - DTG4(i) * fac1(i) * (1.-areas(i))
     LCDTS(I) =     - DTS4(i) * fac1(i) * (   areas(i))

     !pl ground:
     LGDTG(I) =   DTG4(i)
     LGDTC(I) = - DTC4(i) * fac1(i)

     !pl snow:
     LSDTS(I) =   DTS4(i)
     LSDTC(I) = - DTC4(i) * fac1(i)

  ENDDO

return
END SUBROUTINE DELLWF

!##############################################################################
Subroutine DELEF (DTA,CP,ps,em,ea,ros,HRr,fc,fg             &
     ,                ra,rb,rd,rc,rsoil,snow,capac,www      &
     ,                ECDTC,ECDEA,EGDTG,EGDEA,ESDTS         &
     ,                ESDEA,EADEA                           &
     ,                ec,eg,es,fws,hltm,cas_cap_vap         &
     ,                etc,etg                               &
     ,                btc,btg                               &
     ,                areas, gect,geci,gegs,gegi, psy, snofac, hr   &
     ,                len)

  ! Calculation of partial derivatives of canopy and ground latent
  ! heat fluxes with respect to Tc, Tgs, Theta-m, and Qm.
  ! Calculation of initial latent heat fluxes.

  ! The ETC, ETG and so on are the vapor pressure at temps TC, TG and so on
  ! The BTC, BTG are the derivatives of ETC, ETG with relation to TC, TG etc.
  
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       EC             ECT + ECI
  !       EG             EGS + EGI
  !       ECDTC          dEC/dTC
  !       ECDTG          dEC/dTGS
  !       EGDTC          dEG/dTC
  !       EGDTG          dEG/dTGS
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

  INTEGER len
  REAL                                                            &
       ECDTC(len),ECDEA(len),EGDTG(len),EGDEA(len)                &
       ,    ESDTS(len),ESDEA(len),EADEA(len)                      &
       ,    fc(len),fg(len),snow(len,2),capac(len,2),www(len,3)   &
       ,    hrr(len),rst(len),rb(len),rd(len),rc(len),rsoil(len)  &
       ,    ra(len),rds(len), cas_cap_vap(len)                    &
       ,    em(len),etc(len),etg(len),ea(len)                     &
       ,    ps(len),ec(len),eg(len),es(len)                       &
       ,    btc(len),btg(len)                                     &
       ,    sh(len),COG1(len),COC(len),D2(len)                    &
       ,    tc_0(len), tgs_0(len), intg(len), intc(len)           &
       ,    COG2(len), hr(len) ,fws(len),hltm, rcp(len), rcpg(len)&
       ,    ros(len), dta, cp, psy(len), snofac                   &
       ,    limci(len),limgi(len),limct(len),limgs(len)
  REAL cogr(len),cogs(len),areas(len), cpdpsy(len)
  REAL gect(len),geci(len),gegs(len),gegi(len),resrat

  INTEGER i

  !     MODIFICATION FOR SOIL DRYNESS : HR=REL. HUMIDITY IN TOP LAYER

  resrat = 0.5

  DO I=1,len
     hrr(I) = HR(I)
     IF(fg(i).LT.0.5) hrr(i) = 1.0
  ENDDO

  DO I=1,len
     rcp(i)  = ros(i) * cp
     rcpg(i) = rcp(i)/psy(i)              ! this is rho * cp / gamma
     cpdpsy(i) = cp / psy(i)
     rds(i) = rsoil(i) + rd(i)

     !-----------------------------------------------------------------------
     ! CALCULATION OF SURFACE RESISTANCE COMPONENTS, SEE EQUATIONS (64,66)
     !       OF SE-86
     ! the ge?? coefficients come all the way from VNTLAT and are common to
     ! all routines:
     !     gect(i)  =      (1. -wc(i)) /  rc(i)
     !	geci(i)  = epsc(i) * wc(i)  / (RB(I) + RB(I))
     !	gegs(i)  =        (1-wg(i)) / (rds(i))
     !	gegi(i)  = epsg(i) * wg(i)  /  rd(i)
     !-----------------------------------------------------------------------

     COC(I) =  gect(i) + geci(i)
     COG1(i) = (gegi(i) + gegs(i)*HRR(i))
     COG2(i) = (gegi(i) + gegs(i)       )

     !            D2(I)   = 1.0 / RA(I) + COC(I) + COG2(I)
     !-----------------------------------------------------------------------
     !
     !     FLUXES EXPRESSED IN JOULES M-2   CPL WHY ?????
     !
     !      ec         (EC)    : EQUATION (64) , SE-86
     !      eg         (EG)    : EQUATION (66) , SE-86
     !      es         (ES)    : EQUATION (66) , SE-86
     !      ea         (EA)    : EQUATION ????
     !-----------------------------------------------------------------------

     !pl these are the current time step fluxes in J/m2  WHY ?????

     !pl notice that the fluxes are already limited by the altered e*(T) values

     ec(I)  = (etc(I) - ea(i)) * COC(I) *       &
          ros(i)  * dta * cpdpsy(i)

     eg(I)  = ( etg(I) * COG1(I)                &
          - ea(i) * COG2(I)) *             &
          ros(i) * dta * cpdpsy(i)

     es(I)  = ((etg(I) - ea(i))/rd(i) )*                     &
          ros(i) * dta * cpdpsy(i)/snofac
     fws(I) = ((ea(I)  - em(i) ) / ra(i))                    &
          * ros(i) * dta * cpdpsy(i)

     !pl now we do the partial derivatives  these assume W/m2

     !pl for the canopy leaves vapor pressure: W/ (m2* K)
     ECDTC(I) =    btc(I) * COC(I)                       &
          * ros(i) * CPDPSY(i)
     ECDEA(I) = - COC(I) * ros(i) * CPDPSY(i)

     !pl for ground latent heat fluxes: W/ (m2* K)
     EGDTG(I) =   btg(I) * COG1(I)                    &
          * ros(i) * CPDPSY(i)
     EGDEA(I) = - (COG2(I)) * ros(i) * CPDPSY(i)

     !pl for snow latent heat fluxes: W/ (m2* K)
     ESDTS(I) =   btg(I) * ros(i) * CPDPSY(i)/RD(i)

     !pl for snow latent heat fluxes: W/ (m2 * Pa)
     ESDEA(I) = - ros(i) * CPDPSY(i)/RD(i)

     !pl for CAS latent heat fluxes: W/ (m2* Pa)
     EADEA(I) = ros(i) * CPDPSY(i) / ra(i)

  ENDDO

!print*,'delef: em,ea,etg,etc=',em,ea,etg,etc
!print*,'delef: coc,cog1,cog2,dta=',coc,cog1,cog2,dta
!print*,'delef: ros, cpdpsy=',ros,cpdpsy
!print*,'delef: ec,eg,es,fws=',ec,eg,es,fws

return
END SUBROUTINE DELEF

!##############################################################################
Subroutine DELHF (DTA,CP,bps,tm,tg,ts,tc,ta,ros,ra,rb,rd          &
     ,                HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA         &
     ,                HADTA                                       &
     ,                hc, hg, hs, fss, len )

  ! ATTENTION receiving tgs instead of tg !!!
  ! Calculation of partial derivatives of canopy and ground sensible
  ! heat fluxes with respect to Tc, Tgs, and Theta-m.
  ! Calculation of initial sensible heat fluxes.
  
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       HC             CANOPY SENSIBLE HEAT FLUX (J M-2)
  !       HG             GROUND SENSIBLE HEAT FLUX (J M-2)
  !       HS             SNOW   SENSIBLE HEAT FLUX (J M-2)
  !       HA             CAS    SENSIBLE HEAT FLUX (J M-2)
  !       HCDTC          dHC/dTC
  !       HCDTA          dHC/dTA
  !       HGDTG          dHG/dTG
  !       HGDTA          dHG/dTA
  !       HSDTS          dHS/dTS
  !       HSDTA          dHS/dTA
  !       HADTA          dHA/dTA
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

  INTEGER len
  REAL                                                         &
       bps(len),tm(len),tg(len),ts(len),HCDTC(len),HCDTA(len)  &
       ,    HGDTG(len),HGDTA(len),HSDTS(len),HSDTA(len)             &
       ,    HADTA(len),ra(len)                                      &
       ,    rb(len),rd(len),tc(len), ta(len), hc(len), hg(len)      &
       ,    ros(len), dta, cp, fss(len), hs(len)
  !     local variables
  REAL D1, d1i, rai, rbi,  rdi
  INTEGER i

  !-----------------------------------------------------------------------
  !   FLUXES EXPRESSED IN JOULES M-2, although in SIBSLV WE THEN WANT W/m2
  !pl WHY ????
  !
  !pl if we were to keep things simple, there is no need to separate
  !pl HG and HS, but it helps the derivatives keep clean.
  !
  !      HC          (HC)    : EQUATION (63) , SE-86
  !      HG          (HG)    : EQUATION (65) , SE-86
  !      HS          (HS)    : EQUATION (65) , SE-86
  !      HA          (HA)    : EQUATION ???
  !-----------------------------------------------------------------------

  DO I=1,len
     rai = 1.0 / ra(I)
     rbi = 1.0 / rb(I)
     rdi = 1.0 / rd(I)
     !        D1 = rai + rbi + rdi
     !        d1i = 1. / d1

     !pl these are the current time step fluxes in J/m2
     !pl can we change this to W/m2 ???

     HC(I)   = CP * ros(i) * (tc(i) - ta(i)) * rbi * DTA
     HG(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
     HS(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
     fss(I)  = CP * ros(i) * (ta(I) - tm(i)) * rai * DTA

     !pl now we do the partial derivatives

     !pl these are done assuming the fluxes in W/m2

     !pl for canopy leaves sensible heat flux: W/(m2 * K)
     HCDTC(I) =   CP * ros(i) * rbi
     HCDTA(I) = - HCDTC(I)
     !pl for ground and snow sensible heat fluxes: W/(m2 * K)
     HGDTG(I) =   CP * ros(i) * rdi
     HSDTS(I) =   HGDTG(I)
     HGDTA(I) = - HGDTG(I)
     HSDTA(I) = - HGDTG(I)
     !pl for the canopy air space (CAS) sensible heat flux: W/(m2 * K)
     HADTA(I) =   CP * ros(i) * rai

  ENDDO

  !itb...
  !      !print*,'delhf: rbi,rdi,rai,dta=',rbi,rdi,rai,dta
  !      !print*,'delhf: cp, ros=',cp,ros
  !      !print*,'delhf: tm,tc,tg,ta=',tm,tc,tg,ta
  !      !print*,'delhf: hc,hg,hs,fss=',hc,hg,hs,fss

return
END SUBROUTINE DELHF

!##############################################################################
Subroutine soilprop ( td, tg, slamda, shcap, www, poros, ztdep,  &
     asnow, snoww, areas, tf, snomel, nsib, len, nsoil )

implicit none

  !     this routine calculates the soil thermal properties heat capacity
  !         (shcap) and conductivity (slamda). Phase changes are incorporated
  !         into the heat capacity over a range of -0.5C to 0.5C.
  !     slamda(n) is the conductivity between layers n and n+1 / delta z
  !     layer 1 is the deepest layer, layer nsoil is immediately below the
  !         surface
  !     treatment based on Bonan

  !     argument list variables
  INTEGER len, nsib, nsoil
  REAL                                                       &
       td(nsib,nsoil), tg(len), slamda(len,nsoil),         &
       shcap(len,nsoil), www(len,3), poros(len), snomel,   &
       ztdep(nsib,nsoil), areas(len), snoww(len), asnow, tf

  !     local variables
  REAL                                                           &
       shcapu, shcapf, shsnow, klamda(len,nsoil),zsnow(len),   &
       kf, ku, ksnow, ksoil(len), kiw(len,3),kii(len,3),       &
       shcap1, klamda1, klamda2
  INTEGER i,n

  !  snow density assumed 0.25 that of water to convert from Bonan's constants
  !  for ksnow and max effective snow depth
  ksnow = 0.085
  DO i = 1,len
     ksoil(i) = 6.0**(1.-poros(i))
  ENDDO
  DO n = 1,3
     DO i = 1,len
        kiw(i,n) = 0.6**(www(i,n)*poros(i))
        kii(i,n) = 2.2**(www(i,n)*poros(i))
     ENDDO
  ENDDO
  !     heat capacity calculation and conductivity calculation

  !     soil layers 1 and 2 ( deepest layers, use www(3) )
  DO n = 1,2
     DO i = 1,len
        SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,3)*   &
             POROS(i)) * 4.186E6
        shcapf =  0.5 * (1. + poros(i) *            &
             (www(i,3)-1.0)) * 4.186E6
        ku = (ksoil(i)*kiw(i,3)-0.15)*www(i,3) + 0.15
        kf = (ksoil(i)*kii(i,3)-0.15)*www(i,3) + 0.15
        IF(td(i,n).GE.tf+0.5) THEN
           shcap(i,n) = shcapu * ztdep(i,n)
           klamda(i,n) = ku
        ELSE IF(td(i,n).LE.tf-0.5) THEN
           shcap(i,n) = shcapf * ztdep(i,n)
           klamda(i,n) = kf
        ELSE
           shcap(i,n) = (0.5*(shcapu+shcapf) +            &
                poros(i)*www(i,3)*snomel) * ztdep(i,n)
           klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5-tf)
        ENDIF
     ENDDO
  ENDDO

  !     soil layers 3,4,5 and nsoil ( intermediate layers, use www(2) )
  DO n = 3,nsoil
     DO i = 1,len
        SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,2)*       &
             POROS(i)) * 4.186E6
        shcapf =  0.5 * (1. + poros(i) *                &
             (www(i,2)-1.0)) * 4.186E6
        ku = (ksoil(i)*kiw(i,2)-0.15)*www(i,2) + 0.15
        kf = (ksoil(i)*kii(i,2)-0.15)*www(i,2) + 0.15
        IF(td(i,n).GE.tf+0.5) THEN
           shcap(i,n) = shcapu * ztdep(i,n)
           klamda(i,n) = ku
        ELSE IF(td(i,n).LE.tf-0.5) THEN
           shcap(i,n) = shcapf * ztdep(i,n)
           klamda(i,n) = kf
        ELSE
           shcap(i,n) = (0.5*(shcapu+shcapf) +           &
                poros(i)*www(i,2)*snomel) * ztdep(i,n)
           klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5-tf)
        ENDIF
     ENDDO
  ENDDO

  !     soil layer nsoil ( top layer, use www(1) )
  !     if snow covered add additional heat capacity due to snow
  DO i = 1,len
     SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,1)*       &
          POROS(i)) * 4.186E6
     shcapf =  0.5 * (1. + poros(i) *                &
          (www(i,1)-1.0)) * 4.186E6
     ku = (ksoil(i)*kiw(i,1)-0.15)*www(i,1) + 0.15
     kf = (ksoil(i)*kii(i,1)-0.15)*www(i,1) + 0.15
     IF(td(i,nsoil).GE.tf+0.5) THEN
        shcap1 = shcapu
        klamda1 = ku
     ELSE IF(td(i,nsoil).LE.tf-0.5) THEN
        shcap1 = shcapf
        klamda1 = kf
     ELSE
        shcap1 = 0.5*(shcapu+shcapf) +                &
             poros(i)*www(i,1)*snomel
        klamda1 = kf + (ku-kf)*(td(i,nsoil)+0.5-tf)
     ENDIF
     IF(tg(i).GE.tf+0.5) THEN
        klamda2 = ku
     ELSE IF(tg(i).LE.tf-0.5) THEN
        klamda2 = kf
     ELSE
        klamda2 = kf + (ku-kf)*(tg(i)+0.5-tf)
     ENDIF
     zsnow(i) = MIN( MAX(1./asnow,snoww(i)) - 0.02, 0.25E0 )
     shsnow = (0.5 * 4.186e6 * zsnow(i) + shcap1 * 0.02) /         &
          (zsnow(i)+0.02)
     shcap(i,nsoil) =  shcap(i,nsoil) * (1.-areas(i)) + areas(i) *   &
          shcap(i,nsoil)*shsnow*(ztdep(i,nsoil)+zsnow(i)+0.02) /       &
          (shsnow*ztdep(i,nsoil)+shcap(i,nsoil) *                      &
          (zsnow(i)+0.02)) * (ztdep(i,nsoil)+zsnow(i)+0.02)
     klamda1 = ksnow*klamda1 * (0.02+zsnow(i)) /                 &
          (ksnow*0.02 + klamda1*zsnow(i))

     !    soil conductivities / delta z

     slamda(i,nsoil) = (1.-areas(i)) * 2.*klamda(i,nsoil)*klamda2 /   &
          (klamda2*ztdep(i,nsoil)+0.02*klamda(i,nsoil)) +              &
          areas(i) * klamda(i,nsoil)*klamda1 /                &
          (klamda1*ztdep(i,nsoil)*0.05 +                               &
          klamda(i,nsoil)*(zsnow(i)+0.02))

  ENDDO

  DO n = 1,nsoil-1
     DO i = 1,len
        slamda(i,n) = 2.0 * klamda(i,n)*klamda(i,n+1) /          &
             (klamda(i,n)*ztdep(i,n+1)+                  &
             klamda(i,n+1)*ztdep(i,n))
     ENDDO
  ENDDO

return
END SUBROUTINE soilprop

!##############################################################################
Subroutine SIBSLV (DT,GRAV2,CP,HLTM,tg,tsnow,td,slamda,pi             &
     ,                 areas,fac1                                     &
     ,                 VENTMF,BPS,ros,psy                             &
     ,                 cg,ccx,cas_cap_heat,cas_cap_vap                &
     ,                 lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc      &
     ,                 HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA            &
     ,                 HADTA                                          &
     ,                 hc, hg, hs, fss                                &
     ,                 ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                  &
     ,                 ESDEA,EADEA                                    &
     ,                 ec,eg,es,fws                                   &
     ,                 etc,etg                                        &
     ,                 btc,btg                                        &
     ,                 RADT                                           &
     ,                 dtc, dtg, dta, dea                             &
     ,                 len)

! ATTENTION !!! I am calling this with the TGS temp instead of tg.

! Calculation of time increments in Tc, Tgs, Theta-m and Qm using an
! implicit backwards method with explicit coefficients.
! Similar to equations (10-15), SA-92B.
! Longwave feedbacks are now really included.

!pl ATTENTION !!! this is hardwired to work only with Bonan soil.
!pl if force-restore is wanted, we need to include approporiate if loops like
!pl in sibslv.F

  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       DTC            CANOPY TEMPERATURE INCREMENT (K)
  !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
  !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
  !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

  ! Declare parameters
  INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6)   ! Single
  INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=13)  ! Double

  INTEGER len

  REAL                                                           &
       lcdtc(len),lcdtg(len),lcdts(len),lgdtg(len),lgdtc(len)    &
       ,    lsdts(len),lsdtc(len),areas(len)                     &
       ,    dtc4(len),dtg4(len),dts4(len)
  REAL                                                           &
       VENTMF(len),BPS(len)                                      &
       ,    cg(len),ccx(len),z2(len),z1(len)                     &
       ,    HCDTC(len),HCDTA(len)                                &
       ,    HGDTG(len),HGDTA(len),HSDTS(len),HSDTA(len)          &
       ,    HADTA(len)                                           &
       ,    ECDTC(len),ECDEA(len),EGDTG(len),EGDEA(len)          &
       ,    ESDTS(len),ESDEA(len),EADEA(len)                     &
       ,    FSS(len),FWS(len),td(len),RADT(len,3)                &
       ,    tc(len), tg(len),tsnow(len), ta(len)                 &
       ,    ec(len),eg(len),es(len),ea(len)                      &
       ,    etc(len),etg(len)                                    &
       ,    btc(len),btg(len)                                    &
       ,    hc(len), hg(len), hs(len)                            &
       ,    dtc(len), dtg(len,2), dta(len), dea(len)             &
       ,    fac1(len)                                            &
       ,    dt, hltm, cp, grav2, slamda(len), pi                 &
       ,    cas_cap_heat(len),cas_cap_vap(len),ros(len)

  !     local variables
  REAL tmcn2,dti,ddthl,cpdt,TEMV(len),cpdpsy(len),psy(len)  &
       ,AVEC(len,5,6),BVEC(len,5)

  INTEGER i,j,k, error_flag

  !pl  this routine sets up the coupled system of partial
  !pl  differential equations described in Sato et al.,
  !pl  with the exception that now Ta and ea are prognostic
  !pl  variables, and so we have two more equations, reflected
  !pl  in two more lines and two more columns as compared
  !pl  to the old sibslv.F (used for no prognistic CAS calculations)

  !pl          J: /variables
  !pl J: equation/  1     2     3     4      5     6     7     8
  !pl              TC,   TG,   TS, TREF,  EREF,   TA,   EA,  FORCING past t..l.
  !pl 1: TC (vegetation temperature)
  !pl 2: TG (Ground or snow surface temperature)
  !pl 3: TS (Ground or snow surface temperature)
  !pl 4: TREF (lowest model layer/driver data temperature) (not used with RAMS)
  !pl 5: EREF (lowest model layer/driver data water vapor) (not used with RAMS)
  !pl 6: TA (CAS temperature) 
  !pl 7: EA (CAS water vapor)

  DTI  = 1.0 / DT

  DO I=1,len
     TEMV(I) = GRAV2 * VENTMF(i)
     cpdpsy(i) = cp / psy(i)
     !
     !     DTC EQUATION
     !
     AVEC(I,1,1) = ccx(I) * DTI + HCDTC(I) + ECDTC(I) + lcdtc(i)
     AVEC(I,1,2) = LCDTG(i)
     AVEC(I,1,3) = LCDTS(i)
     AVEC(I,1,4) = HCDTA(I)
     AVEC(I,1,5) = ECDEA(I)
     AVEC(I,1,6) = RADT(i,1) - HC(I) * DTI - ec(I) * DTI
     !
     !     DTG EQUATION
     !
     AVEC(i,2,1) = lgdtc(i)
     AVEC(i,2,2) = cg(i)   * DTI + HGDTG(I) + EGDTG(I)         &
          + lgdtg(i) + slamda(i)
     AVEC(i,2,3) = 0.
     AVEC(I,2,4) = hgdta(i)
     AVEC(i,2,5) = egdea(i)
     AVEC(I,2,6) = RADT(i,2) - HG(I) * DTI - eg(I) * DTI         &
          - slamda(i) * (tg(I) - td(i))
     !
     !     DTS EQUATION
     !
     AVEC(i,3,1) = lsdtc(i)
     AVEC(i,3,2) = 0.
     AVEC(i,3,3) = cg(i)   * DTI + HSDTS(I) + ESDTS(I)        &
          + lsdts(i) + slamda(i)
     AVEC(I,3,4) = hsdta(i)
     AVEC(i,3,5) = esdea(i)
     AVEC(I,3,6) = RADT(i,3) - HS(I) * DTI - es(I) * DTI       &
          - slamda(i) * (tg(I) - td(i))
     !
     !     DTA EQUATION
     !
     AVEC(i,4,1) = - HCDTC(i)
     AVEC(i,4,2) = - HGDTG(i) * (1.-areas(i))
     AVEC(i,4,3) = - HSDTS(i) * (   areas(i))
     AVEC(I,4,4) = cas_cap_heat(i)   * DTI                &
          + HADTA(I)  - HCDTA(i)                &
          - (1.-areas(i))*HGDTA(I) - areas(i)*HSDTA(I)
     AVEC(i,4,5) = 0.
     AVEC(I,4,6) = HC(I) * DTI - FSS(I) * DTI         &
          + (1.-areas(i))*HG(I) * DTI        &
          +     areas(i) *HS(I) * DTI
     !
     !     DEA EQUATION
     !
     AVEC(i,5,1) = - ECDTC(i)
     AVEC(i,5,2) = - EGDTG(i) * (1.-areas(i))
     AVEC(i,5,3) = - ESDTS(i) * (   areas(i))
     AVEC(I,5,4) = 0.
     AVEC(i,5,5) = cas_cap_vap(i)   * DTI            &
          + (EADEA(I)  - ECDEA(I)          &
          -  (1.-areas(i))*EGDEA(I)        &
          -      areas(i) *ESDEA(I))
     AVEC(I,5,6) = (EC(I) * DTI  -  FWS(I) * DTI      &
          +  (1.-areas(i))*EG(I) * DTI      &
          +      areas(i) *ES(I) * DTI)
  ENDDO

  !
  !     SOLVE 7 X 8 MATRIX EQUATION
  !

  CALL GAUSS (len,AVEC,5,6,BVEC)

  DO I=1,len
     DTC(I)   = BVEC(I,1)
     DTG(I,1) = BVEC(I,2)   ! this is DTG
     DTG(I,2) = BVEC(I,3)   ! this is DTS
     DTA(i)   = BVEC(I,4)
     DEA(i)   = BVEC(I,5)
  ENDDO

return
END SUBROUTINE SIBSLV

!##############################################################################
Subroutine GAUSS (len,WORK,N,NP1,X)

! SOLVE A LINEAR SYSTEM BY GAUSSIAN ELIMINATION.  DEVELOPED BY
! DR. CHIN-HOH MOENG.  A IS THE MATRIX OF COEFFICIENTS, WITH THE
! VECTOR OF CONSTANTS APPENDED AS AN EXTRA COLUMN.  X IS THE VECTOR
! CONTAINING THE RESULTS.  THE INPUT MATRIX IS NOT DESTROYED.

implicit none

  INTEGER len, n, np1
  REAL WORK(len,N,NP1),X(len,N)
  !     local variables
  REAL TEMV(len,2)
  INTEGER ii, j, i, k, kk, l

  DO II=2,N
     DO J=II,N
        DO I=1,len
           TEMV(I,1) = WORK(I,J,II-1) / WORK(I,II-1,II-1)
        ENDDO
        DO K=1,NP1
           DO I=1,len
              WORK(I,J,K) = WORK(I,J,K) - TEMV(I,1) * WORK(I,II-1,K)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO K=N,2,-1
     DO I=1,len
        TEMV(I,1) = WORK(I,K,NP1) / WORK(I,K,K)
     ENDDO
     KK = K-1
     DO L=KK,1,-1
        DO I=1,len
           TEMV(I,2) = TEMV(I,1) * WORK(I,L,K)
        ENDDO
        DO I=1,len
           WORK(I,L,NP1) = WORK(I,L,NP1) - TEMV(I,2)
        ENDDO
     ENDDO
  ENDDO

  DO II=1,N
     DO I=1,len
        X(I,II) = WORK(I,II,NP1) / WORK(I,II,II)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE GAUSS

!##############################################################################
Subroutine UPDAT2 (snoww,capac,snofac,ect,eci,egi                &
     ,             egs,hlat,www,pi,cg,dtg,dtc,ta,dta,dea         &
     ,             dtt                                           &
     ,             roff,tc,td,tg,bee, poros                      &
     ,             satco,slope,phsat,zdepth,ecmass               &
     ,             egmass,shf,tf,snomel,asnow                    &
     ,             ccx,csoil,chf,hc,hg,areas                     &
     ,             q3l,q3o,qqq,zmelt,cw, len, nsib               &
     ,             nsoil,etc,ea,btc                              &
     ,             geci,ros,cp,psy,gect,etg,btg                  &
     ,             gegs,hr,fg,gegi,rd,rb,hcdtc,hcdta             &
     ,             hgdta,slamda)

implicit none

  !====================================================================
  !
  !  MOVING STORAGE TERM UPDATES (SNOW, CAPAC) HERE FROM ENDTEM, WHICH
  !     NO LONGER EXISTS. FLUXES PREVIOUSLY CALCULATED IN ENDTEM ARE
  !     TAKEN CARE OF IN THE PROGNOSTIC C.A.S. CALCULATIONS, SO WE
  !     MERELY NEED TO TAKE CARE OF STORAGE TERMS NOW.
  !
  !     ITB November 2000
  !=======================================================================
  !
  !     UPDATING OF ALL HYDROLOGICAL PROGNOSTIC VARIABLES.  SNOW AND
  !        RUNOFF CALCULATIONS (SEE ALSO INTER2).  ROUTINES SNOW2 AND
  !        RUN2 OF 1D MODEL ARE INLINED IN THIS CODE.
  !
  !=======================================================================
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       DTC            CANOPY TEMPERATURE INCREMENT (K)
  !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
  !       WWW(3)         GROUND WETNESS
  !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
  !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
  !       ROFF           RUNOFF (MM)
  !
  !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
  !
  !       ECMASS         CANOPY EVAPOTRANSPIRATION (MM)
  !       EGMASS         GROUND EVAPOTRANSPIRATION (MM)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len, nsib, nsoil, i, l, iveg, ksoil,j

  ! INTENT = IN VARIABLES
  REAL, INTENT(in),DIMENSION(len) :: &
              cg             & ! surface layer heat capacity (J m^-2 deg^-1)
       ,      dta            & ! delta canopy air space (CAS) temperature (K)
       ,      dea            & ! delta CAS vapor pressure (Pa)
       ,      tc             & ! canopy temperature (K)
       ,      ta             & ! CAS temperature (K)
       ,      tg             & ! ground surface temperature (K)
       ,      bee            & ! Clapp&Hornberger 'b' exponent
       ,      poros          & ! soil porosity (fraction)
       ,      satco          & ! hydraulic conductivity at saturation (UNITS?)
       ,      slope          & ! cosine of mean slope
       ,      phsat          & ! soil tension at saturation  (UNITS?)
       ,      ccx            & ! canopy heat capacity (J m^-2 deg^-1)
       ,      csoil          & ! soil heat capacity (J m^-2 deg^-1)
       ,      etc            & ! 'E-star' of the canopy - vapor pressure (Pa)
       ,      ea             & ! CAS vapor pressure (Pa)
       ,      btc            & ! d(E(Tc))/d(Tc) - Clausius-Clapyron
       ,      gegs           & ! dry fraction of ground/(fg*rsoil + Rd)
       ,      hr             & ! soil surface relative humidity
       ,      fg             & ! flag indicating direction of vapor pressure
       ! deficit between CAS and ground: 0 => ea>e(Tg)
       !                                 1 => ea<e(Tg)
       ,      gegi           & ! wet fraction of ground/Rd
       ,      rd             & ! ground-CAS resistance
       ,      rb             & ! leaf-CAS resistance
       ,      hcdtc          & ! dHc/dTc
       ,      hcdta          & ! dHc/dTa
       ,      hgdta          & ! dHg/dTa
       ,      slamda         & !
       ,      etg            & ! 'E-star' of the ground surface  (Pa)
       ,      btg            & ! d(E(tg))/d(Tg) - Clausius-Clapyron
       ,      geci           & ! wetted fraction of canopy/2Rb
       ,      gect           & ! dry fraction of canopy/(Rst + 2Rb)
       ,      ros            & ! air density (kg m^-3)
       ,      psy              !

  REAL, INTENT(in)  ::&
       snofac         & !  ___(lat ht of vap)___     (unitless)
       ! (lat ht vap + lat ht ice)
       ,      hlat           & ! latent heat of vaporization of water (J kg^-1)
       ,      pi             & ! 3.1415...
       ,      dtt            & ! timestep (seconds)
       ,      tf             & ! freezing temperature (273.16 K)
       ,      snomel         & ! latent heat of fusion for ice (J m^-3)
       ,      asnow          & ! conversion for kg water to snow depth (16.7)
       ,      cw             & ! water heat capacity (J m^-3 deg^-1)
       ,      cp            ! specific heat of air at const pres (J kg-1 deg-1)

  REAL, INTENT(in),DIMENSION(nsib,3)  ::&
       zdepth          ! soil layer depth * porosity (meters)

  REAL, INTENT(in),DIMENSION(nsib,nsoil)         ::&
       td            ! deep soil temperature (K)

  ! INTENT = OUT VARIABLES
  REAL,INTENT(out),DIMENSION(len)                ::&
              ect        & ! transpiration flux (J m^-2 for the timestep)
       ,      eci        & ! interception flux (veg - CAS) (J m^-2)
       ,      egi        & ! ground interception flux (J m^-2)
       ,      egs        & ! ground evaporative flux (J m^-2)
       ,      ecmass     & ! canopy evaporation (mm)
       ,      egmass     & ! ground evaporation (mm)
       ,      shf        & ! soil heat flux (W m^-2)
       ,      chf        & ! canopy heat flux (W m^-2)
       ,      q3l        & ! 'Liston' drainage from bottom of soillayer 3 (mm)
       ,      q3o        & ! gravitational drainage out of soillayer 3 (mm)
       ,      zmelt        ! depth of melted water (m)

  REAL,INTENT(out),DIMENSION(len,3)              ::&
       QQQ             ! soil layer drainage (mm m^-2 timestep)

  ! INTENT = IN/OUT VARIABLES
  
  REAL,INTENT(inout),DIMENSION(len,3)            ::&
       www             ! soil moisture (% of saturation)

  REAL,INTENT(inout),DIMENSION(nsib,2)           ::&
       snoww           &! snow-interception (1-veg, 2-ground) (meters)
       ,      capac           ! liquid interception

  REAL,INTENT(inout),DIMENSION(len,2)            ::&
       dtg             ! delta ground surface temperature (K)

  REAL,INTENT(inout),DIMENSION(len)              ::&
       roff           & ! runoff (mm)
       ,      hc              &! canopy sensible heat flux (J m^-2)
       ,      hg              &! ground sensible heat flux (J m^-2)
       ,      areas           &! fractional snow coverage (0 to 1)
       ,      dtc             ! delta vegetation temperature (K)

  ! LOCAL VARIABLES
  REAL,DIMENSION(len)                            ::&
              dtgs           & ! snow/dry soil averaged temp increment (K)
       ,      zmelt_total    & ! total melting
       ,      cctt           &
       ,      cct            &
       ,      ts             &
       ,      dts            &
       ,      flux           &
       ,      fluxef         &
       ,      tsnow          &
       ,      dpdw           &
       ,      tgs            &
       ,      dts2           &
       ,      cpdpsy         &
       ,      heaten          ! energy to heat snow to ground temp (J m^-2)

  REAL,DIMENSION(len,3) :: &
            TEMW            &
       ,    TEMWP           &
       ,    TEMWPP

  REAL,DIMENSION(len,2) :: &
            AAA              &
       ,    BBB              &
       ,    CCC

  REAL :: &
            hlati            &
       ,    rsnow            &
       ,    facks            &
       ,    realc            &
       ,    realg            &
       ,    dpdwdz           &
       ,    qmax             &
       ,    qmin             &
       ,    rdenom           &
       ,    denom            &
       ,    props            &
       ,    avkmax           &
       ,    avkmin           &
       ,    div              &
       ,    rsame            &
       ,    pmin             &
       ,    wmin             &
       ,    pmax             &
       ,    wmax             &
       ,    egsdif           &
       ,    ectdif           &
       ,    extrak           &
       ,    facl             &
       ,    dtsg2            &
       ,    dtsg3            &
       ,    cool             &
       ,    zmelt2           &
       ,    dtsg             &
       ,    heat             &
       ,    exmelt           &
       ,    exheat           &
       ,    safe             &
       ,    avheat           &
       ,    avmelt           &
       ,    avk              &
       ,    pows             &
       ,    avex             &
       ,    freeze           &
       ,    ecpot

  REAL :: egpot        &
       ,    hrr        &
       ,    ecidif     & ! actual amount of canopy moisture put into CAS(J m^-2)
       ,    egidif     & ! actual amount of ground interception moisture
                         !  put into CAS (J m^-2))
       ,    egit       & ! temporary ground heat flux holder (J m^-2)
       ,    t1,t2      & ! snow depth measures (intermediate)
       ,    aven       & ! energy difference between actual snow and areas=1
       ,    darea      & ! adjustment to areas
       ,    arean      & ! adjustment to areas
       ,    ectmax     & ! upper bound for transpiratoin (J m^-2)
       ,    egsmax     & ! upper bound for soil evaporation (J m^-2)
       ,    dti        & ! 1/timestep
       ,    cogs1      & ! non snowcovered fraction * soil humidity
       ,    cogs2        ! non snowcovered fraction

  dti = 1./DTT

  DO i=1,len
     tsnow(i) = MIN(tf - 0.01, tg(i))
     cpdpsy(i) = cp/psy(i)
     tgs(i) = (1.0 - areas(i))*tg(i) + areas(i)*tsnow(i)
     dtgs(i) = (1.0-areas(i))*dtg(i,1) + areas(i)*dtg(i,2)
     rsnow = snoww(i,2) / max(snoww(i,2) + capac(i,2) , 1.0E-10)
     !pl this is the potential gradient in Pa
     !pl this WAS realized in sibslv

     ECPOT =     (etc(i) + btc(i)*DTC(i)) - (ea(I) + DEA(i))

     !pl and this is the  INTERCEPTION flux in J/m2
     ECI(i) = ECPOT * geci(i) * ros(i) * CPDPSY(i) * DTT

     !pl and this is the TRANSPIRATION flux in J/m2
     ECT(i) = ECPOT * gect(i) * ros(i) * CPDPSY(i) * DTT

     !pl this is the potential gradient in Pa
     !pl this WAS realized in sibslv25

     EGPOT =   (etg(i) + btg(i)*DTGS(i)) - ( ea(i) + DEA(i))

     !pl and this is the  INTERCEPTION flux in J/m2
     EGI(i) = EGPOT * (gegi(i) * (1.-areas(i)) + areas(i)/rd(i))     &
          * ros(i) * cpdpsy(i) * DTT

     HRR = HR(i)
     IF ( FG(i) .LT. .5 ) HRR = 1.
     COGS1    =  gegs(i) * HRR * (1.-AREAS(i))
     COGS2    =  gegs(i)       * (1.-AREAS(i))

     !pl and this is the EVAPORATIVE flux in J/m2
     EGS(i) =  (etg(i) + btg(i) * dtgs(i)) * COGS1     &
          -(ea(i) +           dea(i)) * COGS2
     EGS(i) = EGS(i) * ros(i) * CPDPSY(i) * DTT

     !itb...make sure you don't evap more than you have...
     EGSMAX = WWW(i,1) * ZDEPTH(i,1) * hlat * 1.e3 * 0.5
     EGS(i) = MIN ( EGS(i), EGSMAX )
     !itb...make sure you don't transpire more water than is in the soil
     ECTMAX = WWW(i,2) * ZDEPTH(i,2) * hlat * 1.e3 * 0.5
     ECT(i) = MIN ( ECT(i), ECTMAX )

     !itb...these fluxes were all realized in sibslv. If positive, they
     !itb...imply transfer of water vapor INTO the CAS. If negative,
     !itb...they imply transfer OUT OF the CAS. We need to adjust
     !itb...the various reserviors as well as the CAS vapor capacity,
     !itb...making sure that none go to negative values.

     !itb...the actual movement of vapor is taken care of in the
     !itb...equations in sibslv. All we do now is adjust the surface
     !itb...and vegetation storage reservoirs to reflect what we have
     !itb...already added or taken out.

     !pl this is the limitation to the ECI flux in J/m2

     ECIdif=MAX(0.0E0,(ECI(i)-(SNOWw(i,1)+CAPAC(i,1))     &
          *1.E3*hlat))

     ECI(i)   =MIN(ECI(i),                               &
          ( (SNOWw(i,1)+CAPAC(i,1))*1.E3*hlat))


     !pl this is the EGI flux in J/m2

     EGIdif=                                                    &
          MAX(0.0E0,EGI(i)-(SNOWw(i,2)+CAPAC(i,2))*1.E3*hlat)   &
          *(1.-RSNOW)


     EGIT  =                                                    &
          MIN(EGI(i), (SNOWw(i,2)+CAPAC(i,2))*1.E3*hlat )        &
          *(1.-RSNOW)


     !itb...print this stuff out, for grins
     !        print*,'updat2: eci,ect,egi,egs,ecidif,egidif'
     !        print*,eci(i),ect(i),egi(i),egs(i),ecidif,egidif

     !----------------------------------------------------------------------
     !     CALCULATION OF INTERCEPTION LOSS FROM GROUND-SNOW. IF SNOW PATCH
     !     SHRINKS, ENERGY IS TAKEN FROM EGI TO WARM EXPOSED SOIL TO TGS.
     !----------------------------------------------------------------------

     T1 = SNOWw(i,2) - 1./ASNOW
     T2 = MAX( 0.E0, T1 )
     AVEN = EGI(i) - T2*hlat*1.E3/SNOFAC
     IF ( (T1-T2)*EGI(i) .GT. 0. ) AVEN = EGI(i)
     DAREA = AVEN/( (TSNOW(i)-TG(i))*CSOIL(i)            &
          - 1./ASNOW*hlat*1.E3/SNOFAC)
     AREAN = AREAS(i) + DAREA
     EGIdif = EGIdif - MIN( 0.E0, AREAN )              &
          *hlat*1.E3/(asnow*SNOFAC)*RSNOW
     DAREA = MAX( DAREA, -AREAS(i) )
     DAREA = MIN( 1.-AREAS(i), DAREA )
     HEATEN(i) = (TSNOW(i)-TG(i))*CSOIL(i)*DAREA*RSNOW
     EGIT = EGIT + ( EGI(i) - HEATEN(i) - EGIdif )*RSNOW
     EGI(i) = EGIT

     !---------------------------------------------------------------------
     !     CALCULATION OF SENSIBLE HEAT FLUXES FOR THE END OF THE TIMESTEP.
     !        SEE FIGURE (2) OF SE-86.  NOTE THAT INTERCEPTION LOSS EXCESS
     !        ENERGIES (ECIDIF, EGIDIF) ARE ADDED.
     !
     !      HC          (HC)    : EQUATION (63) , SE-86
     !      HG          (HGS)   : EQUATION (65) , SE-86
     !----------------------------------------------------------------------
     !
     !        HC(i) = HC(i) + (HCDTC(i)*DTC(i)
     !     &                +  HCDTA(i)*dta(i))*DTT + ECIdif
     !        HG(i) = HG(i) + (HGDTC(i)*DTC(i)
     !     &                +  HGDTA(i)*dta(i))*DTT + EGIdif

     !itb...i've left the leaf one-sided, for now...
     HC(i) = ( (tc(i)+dtc(i)) - (ta(i)+dta(i)) ) /rb(i)         &
          * ros(i) * cp * DTT + ECIdif

     !itb...ground sensible heat flux includes soil and snow by using
     !itb...dtgs
     HG(i) = ( (tg(i)+dtgs(i)) - (ta(i)+dta(i)) ) /rd(i)          &
          * ros(i) * cp * DTT + EGIdif

     CHF(i) = CCX(i) * dti * DTC(i)

  ENDDO

  !----------------------------------------------------------------------
  !     CALCULATION OF STORAGE HEAT FLUXES
  !----------------------------------------------------------------------
  DO i = 1,len
        !  new soil thermodynamic model
        SHF(i) = dti * ( (1.-areas(i))*dtg(i,1) +                   &
             areas(i)*dtg(i,2) ) * cg(i)             &
             + slamda(i) *( TGS(i)+dtgs(i) - TD(i,nsoil) )
  ENDDO

  ksoil = 3
 
  !----------------------------------------------------------------------
  !    INTERCEPTION LOSSES APPLIED TO SURFACE WATER STORES.
  !    EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY
  !    ( HLAT*1000.) LOSS IS IN M M-2. MASS TERMS ARE IN KG M-2 DT-1
  !    INTERCEPTION AND DRAINAGE TREATED IN INTER2.
  !
  !      CAPAC/SNOWW(1) (M-C)   : EQUATION (3)  , SE-86
  !      CAPAC/SNOWW(2) (M-G)   : EQUATION (4)  , SE-86
  !----------------------------------------------------------------------

  hlati = 1. / hlat
  !PL HERE WE DO A CHECK FOR CONDENSATION AND MAKE SURE THAT IT ONLY
  !PL HAPPENS TRHOUGH ECI AND EGI

  DO i = 1,len
     RSNOW = SNOWW(i,1)/max(SNOWW(i,1)+CAPAC(i,1),1.E-10)
     FACKS = 1. + RSNOW * ( SNOFAC-1. )
     IF ( (ECT(i)+ECI(i)) .LE. 0.) THEN
        ECI(i) = ECT(i)+ECI(i)
        ECT(i) = 0.
        FACKS = 1. / FACKS
     ENDIF
     CAPAC(i,1) = CAPAC(i,1)-( 1.-RSNOW )*ECI(i)*FACKS*hlati*0.001
     SNOWW(i,1) = SNOWW(i,1)-     RSNOW  *ECI(i)*FACKS*hlati*0.001
     snoww(i,1) = MAX(snoww(i,1),0.0E0)
     capac(i,1) = MAX(capac(i,1),0.0E0)
     ECMASS(i) = ECI(i)*FACKS *hlati
     zmelt_total(i) = 0.0
  ENDDO
  !      do i = 1,len
  !         if(snoww(i,1).lt.0.0)
  !     *      print *,'snoww1 after updat2 100 ',i,snoww(i,1)
  !         if(capac(i,1).lt.0.0)
  !     *      print *,'capac1 after updat2 100 ',i,capac(i,1)
  !      enddo
  !
  DO i = 1,len
     RSNOW = SNOWW(i,2)/max(SNOWW(i,2)+CAPAC(i,2),1.e-10)
     FACKS = 1. + RSNOW * ( SNOFAC-1. )
     IF ( (EGS(i)+EGI(i)) .LE. 0. ) THEN
        EGI(i) = EGS(i)+EGI(i)
        EGS(i)= 0.
        FACKS = 1. / FACKS
     ENDIF
     CAPAC(i,2) = CAPAC(i,2)-( 1.-RSNOW )*EGI(i)*FACKS*hlati*0.001
     SNOWW(i,2) = SNOWW(i,2)-     RSNOW  *EGI(i)*FACKS*hlati*0.001
     snoww(i,2) = MAX(snoww(i,2),0.0E0)
     capac(i,2) = MAX(capac(i,2),0.0E0)
     EGMASS(i) = EGI(i)*FACKS *hlati
  ENDDO
  !      do i = 1,len
  !         if(snoww(i,2).lt.0.0)
  !     *      print *,'snoww2 after updat2 200 ',i,snoww(i,2)
  !         if(capac(i,2).lt.0.0)
  !     *      print *,'capac2 after updat2 200 ',i,capac(i,2)
  !      enddo
  !
  !----------------------------------------------------------------------
  !    DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE
  !----------------------------------------------------------------------

  DO IVEG = 1,2
     DO i = 1,len
        IF ( (SNOWW(i,iveg)+CAPAC(i,IVEG)) .LE. 0.00001 ) THEN
           WWW(i,1) = WWW(i,1) + (SNOWW(i,IVEG)+CAPAC(i,IVEG)) /     &
                ZDEPTH(i,1)
           CAPAC(i,IVEG) = 0.
           SNOWW(i,IVEG) = 0.
        ENDIF
     ENDDO
  ENDDO
  !      do i = 1,len
  !         if(www(i,1).le.0.0)print *,'www after updat2 1000 ',i,www(i,1)
  !        print*,rsnow , snoww(i,2) , capac(i,2),  snoww(i,1) , capac(i,1)
  !        stop
  !      enddo
  !
                                                               !
  !=======================================================================
  !------------------SNOW2 INLINED-------------------------------------
  !----------------------------------------------------------------------
  !    SNOWMELT / REFREEZE CALCULATION
  !----------------------------------------------------------------------
  !
  !     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES
  !
  !     MODIFICATION DEALS WITH SNOW PATCHES:
  !          TS < TF, TSNOW = TS
  !          TS > TF, TSNOW = TF
  !=======================================================================

  DO IVEG = 1,2

     REALC = (2 - IVEG)*1.
     REALG = (IVEG - 1)*1.

     DO i = 1,len
        CCTT(i) = REALC*CCX (i) +  REALG*CG(i)
        CCT(i)  = REALC*CCX(i)  +  REALG*CSOIL(i)
        TS(i)   = REALC*TC(i)   +  REALG*TG(i)
        DTS(i)  = REALC*DTC(i)  +  REALG*DTG(i,1)
        DTS2(i)  = REALC*DTC(i)  +  REALG*DTG(i,2)
        FLUX(i) = REALC*CHF(i)  +  REALG*              &
             ( (1.-areas(i))*DTG(i,1)+                   &
             areas(i)*dtg(i,2)  )*cg(i) /DTT
        !  fluxef moved up here to conserve energy
        fluxef(i) = ( shf(i) - flux(i)) * realg
        TSNOW(i) = MIN ( TF-0.01, TS(i) )
        ZMELT(i) = 0.
     ENDDO

     DO i = 1,len  ! this scalar loop needs vector optimization
        !itb   print*,'updat2:ts,dts,ts+dts,tf',ts(i),dts(i),ts(i)+dts(i),tf
        IF ( SNOWW(i,IVEG) .GT. 0. ) GO TO 102
        IF ( ( TS(i)+DTS(i)) .GT. TF ) GO TO 502
        !----------------------------------------------------------------------
        !
        !     NO SNOW  PRESENT, SIMPLE THERMAL BALANCE WITH POSSIBLE FREEZING.
        !
        !----------------------------------------------------------------------
        FREEZE = MIN ( 0.E0, (FLUX(i)*DTT-( TF-0.01 - TS(i))      &
             *CCTT(i)))
        SNOWW(i,IVEG) = MIN( CAPAC(i,IVEG), - FREEZE/SNOMEL )
        ZMELT(i) = CAPAC(i,IVEG) - SNOWW(i,IVEG)
        CAPAC(i,IVEG) = 0.
        DTS(i) = DTS(i) + SNOWW(i,IVEG)*SNOMEL/CCTT(i)
        GO TO 502
        !----------------------------------------------------------------------
        !
        !     SNOW PRESENT
        !
        !---------------------------------------------------------------------
102     CONTINUE
        !itb      IF ( TS(i) .LT. TF .AND. (TS(i)+DTS(i)) .LT. TF ) GO TO 502
        IF ( TS(i) .LT. TF .AND. (TS(i)+DTS2(i)) .LT. TF ) GO TO 502
        IF ( ts(i) .GT. TF ) GO TO 202
        !----------------------------------------------------------------
        !
        !     SNOW PRESENT : TS < TF,  TS+DTS > TF
        !
        !------------------------------------------------------------
        AVEX = FLUX(i) - ( TF-0.01 - TS(i) ) * CCTT(i)/DTT
        AVMELT = ( AVEX/SNOMEL * (AREAS(i)*REALG + REALC ) )*DTT
        ZMELT(i) = MIN( AVMELT, SNOWW(i,IVEG) )
        SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)
        AVHEAT = AVEX*( 1.-AREAS(i) )*REALG +                        &
             ( AVMELT-ZMELT(i) )*SNOMEL/DTT
        AREAS(i) = MIN( 0.999E0, ASNOW*SNOWW(i,2) )
        SAFE = MAX( ( 1.-AREAS(i)*REALG ), 1.E-8 )
        DTS(i) = TF-0.01 - TS(i) + AVHEAT / ( CCTT(i)*SAFE )*DTT
        GO TO 502
        !----------------------------------------------------------------------
        !
        !     SNOW PRESENT AND TS > TF : GROUND ONLY.
        !
        !------------------------------------------------------------------
202     CONTINUE

        EXHEAT = CCT(i)*( 1.001-MAX(0.1E0,AREAS(i))) * DTS(i)
        EXMELT = FLUX(i)*DTT - EXHEAT
        HEAT = EXHEAT
        DTSG = EXHEAT / ( CCT(i)*(1.001-AREAS(i) ))
        IF ( (TS(i)+DTSG) .GT. TF ) GO TO 302
        HEAT = ( TF-0.01 - TS(i) ) * ( CCT(i)*(1.-AREAS(i)) )
        DTSG = TF-0.01 - TS(i)

302     EXMELT = EXMELT + EXHEAT - HEAT

        IF( EXMELT .LT. 0. ) GO TO 402
        ZMELT(i) = EXMELT/SNOMEL
        IF( ASNOW*(SNOWW(i,IVEG)-ZMELT(i)) .LT. 1. )                     &
             ZMELT(i) = MAX( 0.E0, SNOWW(i,IVEG) - 1./ASNOW )
        SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)
        !print*,'XX' ,snoww(i,2) , snoww(i,2) , capac(i,2)
        EXMELT = EXMELT - ZMELT(i)*SNOMEL
        ZMELT2 = EXMELT/ ( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )
        ZMELT2 = MIN( ZMELT2, SNOWW(i,IVEG) )
        ZMELT(i) = ZMELT(i) + ZMELT2
        SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT2
        EXMELT = EXMELT - ZMELT2*( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )
        DTS(i)  = DTSG + EXMELT/CCT(i)
        GO TO 502

402     COOL = MIN( 0.E0, TF-0.01 -(TS(i)+DTSG)) * CCT(i)             &
             *(1.-AREAS(i))
        DTSG2 = MAX ( COOL, EXMELT ) / ( CCT(i)*( 1.001-AREAS(i) ) )
        EXMELT = EXMELT - DTSG2*CCT(i)*(1.-AREAS(i))
        DTSG3 =EXMELT/CCTT(i)
        DTS(i) = DTSG + DTSG2 + DTSG3

502     CONTINUE
     ENDDO

     DO i = 1,len
        !itb...patch
        IF(ZMELT(i) < 0.0 ) ZMELT(i) = 0.0
        !itb...patch
        WWW(i,1) = WWW(i,1) + ZMELT(i) /                     &
             ZDEPTH(i,1)

        IF(www(i,1).LT.0.0)zmelt_total(i) = SQRT(www(i,1))

        DTC(i) = DTC(i)*REALG + DTS(i)*REALC
        DTG(i,1) = DTG(i,1)*REALC + DTS(i)*REALG
        zmelt_total(i) = zmelt_total(i) + zmelt(i)
     ENDDO

  ENDDO

  !itb...put zmelt_total into zmelt
  zmelt(:) = zmelt_total(:)

  DO i = 1,len
     !------------------END SNOW2  -------------------------------------
     !----------------------------------------------------------------------
     !    EVAPOTRANSPIRATION LOSSES APPLIED TO SOIL MOISTURE STORE.
     !    EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE, SOIL EVAPORATION..
     !
     !      ECT         (E-DC)  : EQUATION (5,6), SE-86
     !      EGS         (E-S)   : EQUATION (5)  , SE-86
     !----------------------------------------------------------------------

     !PL STEP THREE part II
     !pl we have done the potential ECT and EGS inside ENDTEM25
     !pl now we limit these fluxes according to 1/2 of what is
     !pl in the soil reservoirs, WWW(i,1) for EGS and WWW(i,2) for ECT
     !pl we 'donate' the excess to HC and HG, if any.

     FACL   = hlati*0.001/ZDEPTH(i,2)
     EXTRAK = ECT(i)*FACL
     EXTRAK = MIN( EXTRAK, WWW(i,2)*0.5 )
     ECTDIF = ECT(i) - EXTRAK/FACL
     ECT(i)    = EXTRAK/FACL
     HC(i)     = HC(i) + ECTDIF
     ECMASS(i) = ECMASS(i) + ECT(i)*hlati
     WWW(i,2) = WWW(i,2) - ECT(i)*FACL

     FACL   = 0.001*hlati/ZDEPTH(i,1)
     EXTRAK = EGS(i)*FACL
     EXTRAK = MIN( EXTRAK, WWW(i,1) *0.5 )
     EGSDIF = EGS(i) - EXTRAK/FACL
     EGS(i)    = EXTRAK/FACL
     HG(i)     = HG(i) + EGSDIF
     EGMASS(i) = EGMASS(i) + EGS(i)*hlati
     WWW(i,1) = WWW(i,1) - EGS(i)*FACL

  ENDDO

  !========================================================================
  !------------------RUN2 INLINED-------------------------------------
  !========================================================================
  !----------------------------------------------------------------------
  !    CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO
  !    GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .
  !----------------------------------------------------------------------

  DO I = 1,3
     DO l = 1,len
        TEMW(l,I)   = MAX( 0.03E0, WWW(l,I) )
        TEMWP(l,I)  = TEMW(l,I) ** ( -BEE(l) )
        TEMWPP(l,I) = MIN( 1.E0, TEMW(l,I))**( 2.*BEE(l)+ 3. )
     ENDDO
  ENDDO

  !-----------------------------------------------------------------------
  !    CALCULATION OF GRAVITATIONALLY DRIVEN DRAINAGE FROM W(3) : TAKEN
  !    AS AN INTEGRAL OF TIME VARYING CONDUCTIVITY.
  !
  !     qqq(3) (Q3) : EQUATION (62) , SE-86
  !
  !    QQQ(3) IS AUGMENTED BY A LINEAR LOSS TERM RECOMMENDED BY LISTON (1992)
  !-----------------------------------------------------------------------

  DO i = 1,len
     POWS = 2.*BEE(i)+2.
     qqq(i,3) = TEMW(i,3)**(-POWS) + SATCO(i)/            &
          (ZDEPTH(i,3) )*                         &
          SLOPE(i)*POWS*DTT
     qqq(i,3) = qqq(i,3) ** ( 1. / POWS )
     qqq(i,3) = - ( 1. / qqq(i,3) - WWW(i,3) ) *           &
          ZDEPTH(i,3) / DTT
     qqq(i,3) = MAX( 0.E0, qqq(i,3) )
     q3o(i) = qqq(i,3) * dtt
     qqq(i,3) = MIN( qqq(i,3), WWW(i,3)*                    &
          ZDEPTH(i,3)/DTT )

     Q3l(i) = 0.002*ZDEPTH(i,3)*0.5 / 86400.*          &
          MAX(0.E0,(www(i,3)-0.01)/0.99 )
     qqq(i,3) = qqq(i,3) + q3l(i)
     q3l(i) = q3l(i) * dtt

     !----------------------------------------------------------------------
     !    CALCULATION OF INTER-LAYER EXCHANGES OF WATER DUE TO GRAVITATION
     !    AND HYDRAULIC GRADIENT. THE VALUES OF W(X) + DW(X) ARE USED TO
     !    CALCULATE THE POTENTIAL GRADIENTS BETWEEN LAYERS.
     !    MODIFIED CALCULATION OF MEAN CONDUCTIVITIES FOLLOWS MILLY AND
     !    EAGLESON (1982 ), REDUCES RECHARGE FLUX TO TOP LAYER.
     !
     !      DPDW           : ESTIMATED DERIVATIVE OF SOIL MOISTURE POTENTIAL
     !                       WITH RESPECT TO SOIL WETNESS. ASSUMPTION OF
     !                       GRAVITATIONAL DRAINAGE USED TO ESTIMATE LIKELY
     !                       MINIMUM WETNESS OVER THE TIME STEP.
     !
     !      QQQ  (Q     )  : EQUATION (61) , SE-86
     !             I,I+1
     !            -
     !      AVK  (K     )  : EQUATION (4.14) , ME-82
     !             I,I+1
     !----------------------------------------------------------------------

     WMAX = MAX( WWW(i,1), WWW(i,2), WWW(i,3), 0.05E0 )
     WMAX = MIN( WMAX, 1.E0 )
     PMAX = WMAX**(-BEE(i))
     WMIN = (PMAX-2.*poros(i)/( PHSAT(i)*                 &
          (ZDEPTH(i,1)+2.*ZDEPTH(i,2)+ZDEPTH(i,3))))          &
          **(-1./BEE(i))
     WMIN = MIN( WWW(i,1), WWW(i,2), WWW(i,3), WMIN )
     WMIN = MAX( WMIN, 0.02E0 )
     PMIN = WMIN**(-BEE(i))
     DPDW(i) = PHSAT(i)*( PMAX-PMIN )/( WMAX-WMIN )
  ENDDO

  DO I = 1,2

     DO l = 1,len
        RSAME = 0.
        AVK  = TEMWP(l,I)*TEMWPP(l,I) - TEMWP(l,I+1)*TEMWPP(l,I+1)
        DIV  = TEMWP(l,I+1) - TEMWP(l,I)
        IF ( ABS(DIV) .LT. 1.E-6 ) RSAME = 1.
        AVK = SATCO(l)*AVK /                                       &
             ( ( 1. + 3./BEE(l) ) * DIV + RSAME )
        AVKMIN = SATCO(l) * MIN( TEMWPP(l,I), TEMWPP(l,I+1) )
        AVKMAX = SATCO(l) * MAX( TEMWPP(l,I), TEMWPP(l,I+1) )*1.01
        AVK = MAX( AVK, AVKMIN )
        AVK = MIN( AVK, AVKMAX )

        ! Collatz-Bounoua change to effective hydraulic conductivity making
        ! it 10x harder for water to move up than down if the upper soil layer
        ! is wetter than lower soil layer.
        !----------------------------------------------------------------------
        IF (www(l,i) .LT. www(l,i+1)) avk = 0.1 * avk     ! lahouari

        !----------------------------------------------------------------------
        !     CONDUCTIVITIES AND BASE FLOW REDUCED WHEN TEMPERATURE DROPS BELOW
        !        FREEZING.
        !----------------------------------------------------------------------

        TSNOW(l) = MIN ( TF-0.01, TG(l) )
        TGS(l) = TSNOW(l)*AREAS(l) + TG(l)*(1.-AREAS(l))
        TS(l)    = TGS(l)*(2-I) + TD(l,ksoil)*(I-1)
        PROPS = ( TS(l)-(TF-10.) ) / 10.
        PROPS = MAX( 0.05E0, MIN( 1.0E0, PROPS ) )
        AVK  = AVK * PROPS
        qqq(l,3)  = qqq(l,3) * PROPS

        !----------------------------------------------------------------------
        !     BACKWARD IMPLICIT CALCULATION OF FLOWS BETWEEN SOIL LAYERS.
        !----------------------------------------------------------------------

        DPDWDZ = DPDW(l)*2.*poros(l)/                         &
             ( ZDEPTH(l,I) + ZDEPTH(l,I+1) )
        AAA(l,I) = 1. + AVK*DPDWDZ*                            &
             ( 1./ZDEPTH(l,I)+1./ZDEPTH(l,I+1) )         &
             *DTT
        BBB(l,I) =-AVK *   DPDWDZ * 1./ZDEPTH(l,2)*DTT
        CCC(l,I) = AVK * ( DPDWDZ * ( WWW(l,I)-WWW(l,I+1) ) + 1. +        &
             (I-1)*DPDWDZ*qqq(l,3)*1./ZDEPTH(l,3)*                  &
             DTT )
     ENDDO
  ENDDO

  DO i = 1,len
     DENOM  = ( AAA(i,1)*AAA(i,2) - BBB(i,1)*BBB(i,2) )
     RDENOM = 0.
     IF ( ABS(DENOM) .LT. 1.E-6 ) RDENOM = 1.
     RDENOM = ( 1.-RDENOM)/( DENOM + RDENOM )
     QQQ(i,1)  = ( AAA(i,2)*CCC(i,1) - BBB(i,1)*CCC(i,2) ) * RDENOM
     QQQ(i,2)  = ( AAA(i,1)*CCC(i,2) - BBB(i,2)*CCC(i,1) ) * RDENOM

     !-----------------------------------------------------------------------
     !     UPDATE WETNESS OF EACH SOIL MOISTURE LAYER DUE TO LAYER INTERFLOW
     !        AND BASE FLOW.
     !-----------------------------------------------------------------------

     WWW(i,3) = WWW(i,3) - qqq(i,3)*DTT/ZDEPTH(i,3)
     ROFF(i) = ROFF(i) + qqq(i,3) * DTT
  ENDDO

  DO I = 1,2

     DO l = 1,len
        QMAX   =  WWW(l,I)   *                       &
             (ZDEPTH(l,I)  /DTT) * 0.5
        QMIN   = -WWW(l,I+1) *                       &
             (ZDEPTH(l,I+1)/DTT) * 0.5
        QQQ(l,I) = MIN( QQQ(l,I),QMAX)
        QQQ(l,I) = MAX( QQQ(l,I),QMIN)
        WWW(l,I)   =   WWW(l,I)   -                  &
             QQQ(l,I)/ZDEPTH(l,I) *DTT
        WWW(l,I+1) =   WWW(l,I+1) +                  &
             QQQ(l,I)/ZDEPTH(l,I+1)*DTT
     ENDDO
  ENDDO

  !      do i = 1,len
  !         if(www(i,1).le.0.0)
  !	 print *,'www after updat2 ',i,www(i,1)
  !	print*,'XX' ,snoww(i,2) , snoww(i,2) , capac(i,2)
  !      enddo

return
END SUBROUTINE UPDAT2

!##############################################################################
Subroutine soiltherm (td,dtd,tgs,dtg,slamda,shcap,ztdep,dt,nsib,len,nsoil)

implicit none

! DTD            DEEP SOIL TEMPERATURE INCREMENT (K)
! This routine calculates the temperature increments dtd for
!  the soil, based on the soil thermodynamic model of Bonan.
! Layer 1 is the deepest layer, layer nsoil is immediately below the surface.
! The time step is crank-nicholson.
! A tridiagonal matrix system is solved.

  !     argument list variables
  INTEGER len, nsib, nsoil
  REAL                                                   &
       dtd(len,nsoil), tgs(len), slamda(len,nsoil),    &
       shcap(len,nsoil), dtg(len), ztdep(len,nsoil),   &
       td(nsib,nsoil), dt

  !     local variables    a(n)*dtd(n-1)+b(n)*dtd(n)+c(n)*dtd(n+1) = r(n)
  REAL &
             a(len,nsoil)    &! lower sub-diagonal
       ,     b(len,nsoil)    &! diagonal
       ,     c(len,nsoil)    &! upper sib-diagonal
       ,     r(len,nsoil)    &! right hand side
       ,     lamtem, rtem , dti, fac
  INTEGER i, n

  dti = 1. / dt   ! inverse time step

  !     construct matrix
  DO n = 1,nsoil
     DO i = 1,len
        b(i,n) = shcap(i,n) * dti
        r(i,n) = 0.0
     ENDDO
  ENDDO
  DO n = 1,nsoil-1
     DO i = 1,len
        lamtem = -0.5 * slamda(i,n)
        rtem = slamda(i,n) * (td(i,n+1) - td(i,n))
        a(i,n+1) = lamtem
        c(i,n) = lamtem
        b(i,n) = b(i,n) - c(i,n)
        b(i,n+1) = b(i,n+1) - a(i,n+1)
        r(i,n) = r(i,n) + rtem
        r(i,n+1) = r(i,n+1) - rtem
     ENDDO
  ENDDO
  DO i = 1,len
     r(i,nsoil) = r(i,nsoil) + slamda(i,nsoil) * (tgs(i)+dtg(i)   &
          - td(i,nsoil))
  ENDDO

  !     eliminate lower diagonal
  DO n = 1,nsoil - 1
     DO i = 1,len
        fac = a(i,n+1) / b(i,n)
        b(i,n+1) = b(i,n+1) - c(i,n) * fac
        r(i,n+1) = r(i,n+1) - r(i,n) * fac
     ENDDO
  ENDDO
  !     back-substitution
  DO i = 1,len
     dtd(i,nsoil) = r(i,nsoil) / b(i,nsoil)
  ENDDO
  DO n = nsoil-1,1,-1
     DO i = 1,len
        dtd(i,n) = (r(i,n) - c(i,n) * dtd(i,n+1)) / b(i,n)
     ENDDO
  ENDDO

return
END SUBROUTINE soiltherm

!##############################################################################
Subroutine addinc (grav2,cp,dtt,hc,hg,                &
     ps, bps,ecmass,psy,rho, hltm, cas_cap_heat,      &
     cas_cap_vap,                                     &
     egmass,fss,fws,hflux,etmass,                     &
     td,thm,ts,sh,tc,tg,ta,ea,ra,em,sha,              &
     dtd,dtc,dtg,dta,dea,drst,rst,bintc, len, nsib, nsoil)

implicit none

  !=======================================================================
  !        Add prognostic variable increments to prognostic variables
  !           and diagnose evapotranspiration and sensible heat flux.
  !
  !        Modified for multitasking - introduced gather/scatter indices
  !           - DD 951206
  !=======================================================================
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       TC             CANOPY TEMPERATURE (K)
  !       TG             GROUND SURFACE TEMPERATURE (K)
  !       TD             DEEP SOIL TEMPERATURE (K)
  !       THM            MIXED LAYER POTENTIAL TEMPERATURE (K)
  !       QM (gsh)       MIXED LAYER MIXING RATIO (KG KG-1)
  !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
  !pl now FWS mm/s, not mm !
  !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
  !       rst (FSS)      STOMATAL RESISTANCE (S M-1)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len, nsib, nsoil
  REAL                                                         &
       grav2,cp,dtt,hc(len),hg(len),bps(len),ecmass(len),     &
       egmass(len),fss(len),fws(len),ps(len),                 &
       ta(len),ea(len),ra(len),em(len),hflux(len),etmass(len),&
       td(nsib,nsoil),thm(len),sh(len),tc(len),tg(len),       &
       dtd(len,nsoil),dtc(len),dtg(len),DRST(len),            &
       rst(len),bintc(len),psy(len), rho(len), z1(len),z2(len), &
       sha(len),ts(len),dta(len),dea(len)
  !     local variables
  INTEGER i, n
  REAL ten, hltm, cas_cap_heat(len), cas_cap_vap(len)

  ten = 10.0

  DO I=1,len
     !	 print*,dtc(i),dtg(i),dta(i),len,tc(i)

     tc(i)   = tc(i) + dtc(i)
     tg(i)   = tg(i) + dtg(i)
     ta(i)   = ta(i) + dta(i)
     ea(i)   = ea(i) + dea(i)

     IF(tg(i).LT.0.0 .OR.   &
          ta(i).LT.0.0 .OR.   &
          ea(i).LT.0.0 .OR.   &
          tc(i).LT.0.0)THEN
        print*,'BAD Ta OR ea VALUE:'
        print*,'SiB point:',i
        print*,'ea:',ea(i)-dea(i),dea(i),ea(i)
        print*,'Ta:',ta(i)-dta(i),dta(i),ta(i)
        print*,'Tg:',tg(i)-dtg(i),dtg(i),tg(i)
        print*,'Tc:',tc(i)-dtc(i),dtc(i),tc(i)
        print*,' '
     ENDIF

     ! print*,'addinc: new tc,tg,ta,ea=',tc(1),tg(1),ta(1),ea(1)
     ! stop

     !Now do the flux from the CAS to the ref level
     !vidale et al. (1999) equation ??
     !Here we are using W/m2

     FSS(i) = rho(i)*cp * (ta(i)-ts(i)) / ra(i)
     hflux(i) = fss(i)

     SH(i)  = 0.622 / ( ps(i)/em(i) -1.)
     SHa(i) = 0.622 / ( ps(i)/ea(i) -1.)

     !This is the latent heat flux from the CAS
     !instead of using W/m2 we stick to kg/m^2/s
     !the conversion is then done at the output, for
     !instance in RAMS module control
     !vidale et al. (1999) equation ??

     !so, here we have want W/m2
     !in the next equation we need to multiply by cpdpsy

     fws(i) = (ea(i) - em(i)) / ra(i) * cp * rho(i) / psy(i)

     !But now let us go back to mm/s (or kg/(m2 s)) for the water flux,
     !in order to keep to the (confusing) system we had before

     fws(i) = fws(i) / hltm
     etmass(i) = fws(i)

     rst(i) = rst(i) + drst(i)
     ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
     ! ---- c.zhang, 2/3/93
     rst(i)=MIN( 1./bintc(i), rst(i) )
     rst(i)=MAX( ten, rst(i) )
  ENDDO
  DO n = 1,nsoil
     DO I=1,len
        TD(i,n)    = TD(i,n) + dtd(i,n)
     ENDDO
  ENDDO

return
END SUBROUTINE addinc

!##############################################################################
Subroutine INTER2 (ppc,ppl,snoww,capac,www                        &
     ,   pie,satcap,cw,tc,tg,clai,zlt,chil,roff                   &
     ,   snomel,zdepTH,tm,tf,asnow,csoil,satco,dtt,vcover,roffo   &
     ,   zmelt, len, nsib, exo)

implicit none

  INTEGER len, nsib,l
  REAL snoww(nsib,2),capac(nsib,2),satcap(len,2),WWW(len,3),     &
       ZDEPTH(nsib,3), ppc(len), ppl(len), pie, cw, tc(len),     &
       tg(len), clai, zlt(len), chil(len), roff(len), snomel,    &
       tm(len), tf, asnow, csoil(len), satco(len), dtt,          &
       vcover(len), roffo(len), zmelt(len), exo(len), excess

  !=======================================================================
  !
  !     CALCULATION OF  INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW
  !     INCORPORATING EFFECTS OF PATCHY SNOW COVER AND TEMPERATURE
  !     ADJUSTMENTS.
  !
  !----------------------------------------------------------------------
  !
  !     (1) NON-UNIFORM PRECIPITATION
  !         CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY
  !         RELATIONSHIP :-
  !
  !                   F(X) = A*EXP(-B*X)+C
  !
  !         THROUGHFALL, INTERCEPTION AND INFILTRATION
  !         EXCESS ARE PHUNCTIONAL ON THIS RELATIONSHIP
  !         AND PROPORTION OF LARGE-SCALE PPN.
  !         REFERENCE: SA-89B, APPENDIX.
  !
  !     (2) REORGANISATION OF SNOWMELT AND RAIN-FREEZE PROCEDURES.
  !               ROUTINE ADJUST
  !
  !     (3) ADDITIONAL CALCULATION FOR PARTIAL SNOW-COVER CASE.
  !               ROUTINE PATCHS
  !
  !     (4) REORGANISATION OF OVERLAND FLOW.
  !         REFERENCE: SA-89B, APPENDIX.
  !
  !     (5) MODIFIED CALCULATION OF SOIL HEAT CAPACITY AND
  !         CONDUCTIVITY.
  !
  !=======================================================================
  !      1D MODEL ROUTINE PATCHS INLINED.
  !----------------------------------------------------------------------

  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       ROFF           RUNOFF (MM)
  !       TC             CANOPY TEMPERATURE (K)
  !       TG             GROUND SURFACE TEMPERATURE (K)
  !       WWW(1)         GROUND WETNESS OF SURFACE LAYER
  !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
  !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !     local variables
  REAL                                                           &
       PCOEFS(2,2), bp, totalp, ap(len), cp(len), thru(len),    &
       fpi(len), realc, realg, xss, xsc, chiv, aa, bb,          &
       exrain, zload, xs, arg, tex(len), tti(len), snowwp(len),  &
       capacp(len), spechc(len), pinf(len), ts(len), equdep,     &
       dcap, tsd, ex, dareas, rhs, areas, snowhc, p0(len)
  INTEGER i, iveg

  DATA PCOEFS(1,1)/ 20. /, PCOEFS(1,2)/ .206E-8 /,                &
       PCOEFS(2,1)/ 0.0001 /, PCOEFS(2,2)/ 0.9999 /, BP /20. /

  !-----------------------------------------------------------------------
  !     PREC ( PI-X )   : EQUATION (C.3), SA-89B
  !-----------------------------------------------------------------------

  DO i = 1,len
     roffo(i) = 0.0
     zmelt(i) = 0.0
     AP(i) = PCOEFS(2,1)
     CP(i) = PCOEFS(2,2)
     TOTALP = (PPC(i) + PPL(i)) * dtt
     IF( SNOWW(i,1) .GT. 0. .OR. SNOWW(i,2) .GT. 0.      &
          .OR. TM(i) .LT. TF ) PPC(i) = 0.
     PPL(i) = TOTALP/dtt - PPC(i)
     IF(TOTALP.GE.1.E-8) THEN
        AP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,1) +           &
             PPL(i)*dtt/TOTALP * PCOEFS(2,1)
        CP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,2) +            &
             PPL(i)*dtt/TOTALP * PCOEFS(2,2)
     ENDIF

     THRU(i) = 0.
     FPI(i)  = 0.

     !----------------------------------------------------------------------
     !     PRECIP INPUT INTO INTER2 IN M/SEC; TOTALP IS IN METERS
     !----------------------------------------------------------------------

     P0(i) = TOTALP
  ENDDO

  DO IVEG = 1,2
     REALC = 2. - IVEG
     REALG = IVEG - 1.

     DO i = 1,len

        XSC = MAX(0.E0, CAPAC(i,IVEG) - SATCAP(i,IVEG) )
        CAPAC(i,IVEG) = CAPAC(i,IVEG) - XSC
        XSS = MAX(0.E0, SNOWW(i,IVEG) - SATCAP(i,IVEG) ) * REALC
        SNOWW(i,IVEG) = SNOWW(i,IVEG) - XSS
        ROFF(i) = ROFF(i) + XSC + XSS

        CAPACP(i) = CAPAC(i,IVEG)
        SNOWWP(i) = SNOWW(i,IVEG)

        SPECHC(i) =                                                    &
             MIN( 0.05E0, ( CAPAC(i,IVEG) + SNOWW(i,IVEG) ) ) * CW    &
             + REALC * ZLT(i) * CLAI + REALG * CSOIL(i)
        TS(i) = TC(i) * REALC + TG(i) * REALG

        !----------------------------------------------------------------------
        !    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)
        !
        !     TTI ( D-D )     : EQUATION (C.4), SA-89B
        !     XS  ( X-S )     : EQUATION (C.7), SA-89B
        !     TEX ( D-C )     : EQUATION (C.8), SA-89B
        !
        !----------------------------------------------------------------------

        CHIV = CHIL(i)
        IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01
        AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV
        BB = 0.877 * ( 1. - 2. * AA )
        EXRAIN = AA + BB

        ZLOAD = CAPAC(i,IVEG) + SNOWW(i,IVEG)
        FPI(i)=( 1.-EXP( - EXRAIN*ZLT(i)/VCOVER(i) ) )*            &
             VCOVER(i)*REALC + REALG
        TTI(i) = P0(i) * ( 1.-FPI(i) )
        XS = 1.
        IF ( P0(i) .GE. 1.E-9 ) THEN
           ARG =  ( SATCAP(i,IVEG)-ZLOAD )/             &
                ( P0(i)*FPI(i)*AP(i) ) -CP(i)/AP(i)
           IF ( ARG .GE. 1.E-9 ) THEN
              XS = -1./BP * LOG( ARG )
              XS = MIN( XS, 1.E0 )
              XS = MAX( XS, 0.E0 )
           ENDIF
        ENDIF
        TEX(i) = P0(i)*FPI(i) *                                &
             ( AP(i)/BP*( 1.- EXP( -BP*XS )) + CP(i)*XS ) -    &
             ( SATCAP(i,IVEG) - ZLOAD ) * XS
        TEX(i) = MAX( TEX(i), 0.E0 )
     ENDDO

     !-------------------------------------------------------------
     !    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION
     !-----------------------------------------------------------

     IF ( IVEG .EQ. 1 ) THEN

        !INLINE
        DO i = 1,len
           thru(i) = TTI(i) + TEX(i)
           PINF(i) = P0(i) - THRU(i)
           IF( TM(i).GT.TF ) THEN
              CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i)
           ELSE
              SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)
           ENDIF

           CALL ADJUST (Tc(i), SPECHC(i), CAPACP(i), SNOWWP(i),      &
                IVEG, capac(i,1), snoww(i,1), tm(i), tf,             &
                snomel, www(i,1), zdepth(i,1),                       &
                satcap(i,1), cw, nsib, len )

           P0(i) = THRU(i)
        ENDDO

     ELSE

        !INLINE
        DO i = 1,len
           IF ( TG(i) .GT. TF .AND. SNOWW(i,2) .GT. 0. ) THEN

              !=============================================================
              !------------------PATCHS INLINED-----------------------------
              !=============================================================
              !
              !CALCULATION OF EFFECT OF INTERCEPTED SNOW AND RAINFALL ON GROUND
              !PATCHY SNOWCOVER SITUATION INVOLVES COMPLEX TREATMENT TO KEEP
              !ENERGY CONSERVED.
              !
              !==============================================================
              !MARGINAL SITUATION: SNOW EXISTS IN PATCHES AT TEMPERATURE TF
              !WITH REMAINING AREA AT TEMPERATURE TG > TF.
              !--------------------------------------------------------------

              PINF(i) = P0(i)
              THRU(i) = 0.
              SNOWHC = MIN( 0.05E0, SNOWW(i,2) ) * CW
              areas = MIN( 1.E0,(ASNOW*SNOWW(i,2)) )
              IF( TM(i) .LE. TF ) THEN
                 !-----------------------------------------------------------
                 !     SNOW FALLING ONTO AREA
                 !---------------------------------------------------------
                 RHS = TM(i)*PINF(i)*CW + TF*(SNOWHC +          &
                      CSOIL(i)*areas)                         &
                      + TG(i)*CSOIL(i)*(1.-areas)
                 DAREAS = MIN( ASNOW*PINF(i), ( 1.-areas ) )
                 EX = RHS - TF*PINF(i)*CW -                      &
                      TF*(SNOWHC + CSOIL(i)*(areas + DAREAS))    &
                      - TG(i)*CSOIL(i)*(1.-areas-DAREAS)
                 IF( (areas+DAREAS) .GE. 0.999 )                 &
                      TG(i) = TF - 0.01
                 IF( EX .GE. 0. ) THEN
                    !----------------------------------------------------------
                    !EXCESS ENERGY IS POSITIVE, SOME SNOW MELTS AND INFILTRATES
                    !----------------------------------------------------------
                    ZMELT(i) = EX/SNOMEL
                    IF( ASNOW*(SNOWW(i,2) + PINF(i) - ZMELT(i))       &
                         .LE. 1. ) THEN
                       ZMELT(i) = 0.
                       IF( ASNOW*(SNOWW(i,2) + PINF(i)) .GE. 1. )      &
                            ZMELT(i) = ( ASNOW*(SNOWW(i,2) +          &
                            PINF(i)) - 1. ) / ASNOW
                       ZMELT(i) = ( EX - ZMELT(i)*SNOMEL )/            &
                            ( SNOMEL + ASNOW*CSOIL(i)*                 &
                            (TG(i)-TF) ) + ZMELT(i)
                    ENDIF
                    SNOWW(i,2) =  SNOWW(i,2) + PINF(i) - ZMELT(i)
                    WWW(i,1) = WWW(i,1) + ZMELT(i)/ZDEPTH(i,1)
                 ELSE
                    !----------------------------------------------------------
                    !EXCESS ENERGY IS NEGATIVE,
                    !BARE GROUND COOLS TO TF, THEN WHOLE
                    !AREA COOLS TOGETHER TO LOWER TEMPERATURE.
                    !----------------------------------------------------------
                    TSD = 0.
                    IF( (areas+DAREAS) .LE. 0.999 )               &
                         TSD = EX/(CSOIL(i)*( 1.-areas-DAREAS))   &
                         + TG(i)
                    IF( TSD .LE. TF )                             &
                         TSD = TF + ( EX - (TF-TG(i))*              &
                         CSOIL(i)*(1.-areas-DAREAS) )        &
                         /(SNOWHC+PINF(i)*CW+CSOIL(i))
                    TG(i) = TSD
                    SNOWW(i,2) = SNOWW(i,2) + PINF(i)
                 ENDIF
              ELSE
                 !-------------------------------------------------------------
                 !     RAIN FALLING ONTO AREA
                 !-------------------------------------------------------------
                 !-------------------------------------------------------------
                 !     RAIN FALLS ONTO SNOW-FREE SECTOR FIRST.
                 !-------------------------------------------------------------
                 TSD = TF - 0.01
                 IF ( areas .LT. 0.999 ) TSD =                &
                      ( TM(i)*PINF(i)*CW +               &
                      TG(i)*CSOIL(i) )                 &
                      /  ( PINF(i)*CW + CSOIL(i) )
                 TG(i) = TSD
                 WWW(i,1)= WWW(i,1)+PINF(i)*(1.-areas)/         &
                      ZDEPTH(i,1)
                 !-------------------------------------------------------------
                 !     RAIN FALLS ONTO SNOW-COVERED SECTOR NEXT.
                 !-------------------------------------------------------------
                 EX = ( TM(i) - TF )*PINF(i)*CW*areas
                 DCAP = -EX / ( SNOMEL + ( TG(i)-TF )*           &
                      CSOIL(i)*ASNOW )
                 IF( (SNOWW(i,2) + DCAP) .GE. 0. ) THEN
                    WWW(i,1) = WWW(i,1)+(PINF(i)*areas-DCAP)/     &
                         ZDEPTH(i,1)
                    SNOWW(i,2) = SNOWW(i,2) + DCAP
                 ELSE
                    TG(i) = ( EX - SNOMEL*SNOWW(i,2) -             &
                         ( TG(i)-TF )*CSOIL(i)*areas ) /          &
                         CSOIL(i) + TG(i)
                    WWW(i,1)=WWW(i,1)+(SNOWW(i,2)+PINF(i)*         &
                         areas)/zdepth(i,1)
                    CAPAC(i,2) = 0.
                    SNOWW(i,2) = 0.
                 ENDIF

              ENDIF
              !
              !----------------------------------------------------------------
              !---------------------END OF PATCHS -----------------------------
              !----------------------------------------------------------------

           ELSE

              THRU(i) = TTI(i) + TEX(i)
              IF ( TG(i) .LE. TF .OR. TM(i) .LE. TF )        &
                   THRU(i) = 0.
              PINF(i) = P0(i) - THRU(i)
              IF( TM(i) .GT. TF )THEN
                 CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i)

                 !-------------------------------------------------------------
                 !
                 !    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )
                 !
                 !     ROFF( R-I )     : EQUATION (C.13), SA-89B
                 !
                 !-------------------------------------------------------------

                 EQUDEP = SATCO(i) * DTT

                 XS = 1.
                 IF ( THRU(i) .GE. 1.E-9 ) THEN
                    ARG = EQUDEP / ( THRU(i) * AP(i) ) -CP(i)/AP(i)
                    IF ( ARG .GE. 1.E-9 ) THEN
                       XS = -1./BP * LOG( ARG )
                       XS = MIN( XS, 1.E0 )
                       XS = MAX( XS, 0.E0 )
                    ENDIF
                 ENDIF
                 ROFFO(i) = THRU(i) * ( AP(i)/BP *                 &
                      ( 1.-EXP( -BP*XS )) + CP(i)*XS ) -EQUDEP*XS
                 ROFFO(i) = MAX ( ROFFO(i), 0.E0 )
                 ROFF(i) = ROFF(i) + ROFFO(i)
                 WWW(i,1) = WWW(i,1) +                            &
                      (THRU(i) - ROFFO(i)) / ZDEPTH(i,1)
              ELSE
                 SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)
              ENDIF

              CALL ADJUST ( Tg(i), SPECHC(i), CAPACP(i), SNOWWP(i),   &
                   IVEG, capac(i,1), snoww(i,1), tm(i), tf,            &
                   snomel, www(i,1), zdepth(i,1),                      &
                   satcap(i,1), cw, nsib, len  )

           ENDIF
        ENDDO
     ENDIF   ! if(iveg.eq.1)

     !     make either all capac or all snow

     DO i = 1,len
        IF(capac(i,iveg).GT.snoww(i,iveg)) THEN
           capac(i,iveg) = capac(i,iveg) + snoww(i,iveg)
           snoww(i,iveg) = 0.0
        ELSE
           snoww(i,iveg) = snoww(i,iveg) + capac(i,iveg)
           capac(i,iveg) = 0.0
        ENDIF
     ENDDO

  ENDDO

  DO i = 1,len
     exo(i) = 0.0
  ENDDO
  DO I = 1,3
     DO l = 1,len
        EXCESS = MAX(0.E0,(WWW(l,I) - 1.))
        WWW(l,I) = WWW(l,I) - EXCESS
        exo(l) = exo(l) + EXCESS * ZDEPTH(l,I)

        ! Collatz-Bounoua put excess water into runoff according to
        ! original sib2 offline code .

        roff(l) = roff(l) + EXCESS * ZDEPTH(l,I)    ! lahouari
     ENDDO

  ENDDO

return
END SUBROUTINE INTER2

!##############################################################################
Subroutine ADJUST ( TS, SPECHC, CAPACP, SNOWWP, IVEG         &
     ,     capac,snoww,tm,tf,snomel,www,zdepth               &
     ,     satcap,cw,nsib, len)

implicit none

  INTEGER len, nsib, iveg
  REAL                                                           &
       capac(nsib,2),snoww(nsib,2),www,zdepth,satcap(len,2),     &
       ts, spechc, capacp, snowwp, tm, tf, snomel,               &
       cw

  !=======================================================================
  !
  !     TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION
  !
  !=======================================================================
  !
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       TC             CANOPY TEMPERATURE (K)
  !       WWW(1)         GROUND WETNESS OF SURFACE LAYER
  !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
  !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !     local variables
  REAL    &
       freeze, diff, ccp, cct, tsd, tta, ttb, cca, ccb, ccc, xs
  FREEZE = 0.
  DIFF = ( CAPAC(1,IVEG)+SNOWW(1,IVEG) - CAPACP-SNOWWP )*CW

  DIFF=MAX(DIFF,0.)

  CCP = SPECHC
  CCT = SPECHC + DIFF

  TSD = ( TS * CCP + TM * DIFF ) / CCT

  IF ( ( TS .GT. TF .AND. TM .LE. TF ) .OR.           &
       ( TS .LE. TF .AND. TM .GT. TF ) )THEN

     TTA = TS
     TTB = TM
     CCA = CCP
     CCB = DIFF
     IF ( TSD .LE. TF ) THEN

        !----------------------------------------------------------------------
        !    FREEZING OF WATER ON CANOPY OR GROUND
        !----------------------------------------------------------------------

        CCC = CAPACP * SNOMEL
        IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW
        TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT

        FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
        FREEZE = (MIN ( CCC, FREEZE )) / SNOMEL
        IF(TSD .GT. TF)TSD = TF - 0.01

     ELSE

        !----------------------------------------------------------------------
        !    MELTING OF SNOW ON CANOPY OR GROUND, WATER INFILTRATES.
        !----------------------------------------------------------------------

        CCC = - SNOWW(1,IVEG) * SNOMEL
        IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW
        !
        TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT
        !
        FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
        FREEZE = (MAX( CCC, FREEZE )) / SNOMEL
        IF(TSD .LE. TF)TSD = TF - 0.01
        !
     ENDIF
  ENDIF
  SNOWW(1,IVEG) = SNOWW(1,IVEG) + FREEZE
  CAPAC(1,IVEG) = CAPAC(1,IVEG) - FREEZE
  snoww(1,IVEG) = MAX(snoww(1,IVEG),0.0E0)
  capac(1,IVEG) = MAX(capac(1,IVEG),0.0E0)

  XS = MAX( 0.E0, ( CAPAC(1,IVEG) - SATCAP(1,IVEG) ) )
  IF( SNOWW(1,IVEG) .GE. 0.0000001 ) XS = CAPAC(1,IVEG)
  WWW = WWW + XS / ZDEPTH
  CAPAC(1,IVEG) = CAPAC(1,IVEG) - XS
  TS = TSD

return
END SUBROUTINE ADJUST

!##############################################################################
Subroutine balan ( IPLACE, tau, zdepth, www         &
     ,         capac, ppc                             &
     ,         ppl, roff, etmass, totwb, radt            &
     ,         chf, shf, dtt, ect, eci, egs, egi           &
     ,         hc, hg, heaten, hflux, snoww, thm, tc, tg, tgs, td, ps, kapa &
     ,         len, ioffset, nsoil )

implicit none

  INTEGER len,  iplace, nsoil
  INTEGER ioffset
  REAL                                                          &
       www(len,3),capac(len,2),zdepth(len,3),radt(len,2),       &
       snoww(len,2), kapa, tau, ppc(len), ppl(len)              &
       ,    roff(len), etmass(len), totwb(len), chf(len), shf(len)   &
       ,    dtt, ect(len), eci(len), egs(len), egi(len), hc(len)     &
       ,    hg(len), heaten(len), hflux(len), thm(len), tc(len)      &
       ,    tg(len), tgs(len), td(len, nsoil), ps(len)

  !     local variables
  INTEGER i, j, igp, jgp, nerror, indxerr(len), n
  REAL                                                        &
       endwb(len), errorw(len), pcmeter(len), plmeter(len)    &
       ,    emeter(len), cbal(len), gbal(len), errore(len)         &
       ,    zlhs(len), zrhs(len), tm
  !=======================================================================
  !
  !     ENERGY AND WATER BALANCE CHECK.
  !
  !-----------------------------------------------------------------------
  IF( IPLACE .EQ. 1 ) THEN

     DO i = 1,len
        ETMASS(i) = 0.
        ROFF(i)   = 0.

        TOTWB(i) = WWW(i,1) * ZDEPTH(i,1)                &
             + WWW(i,2) * ZDEPTH(i,2)                &
             + WWW(i,3) * ZDEPTH(i,3)                &
             + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2)
     ENDDO

  ELSE

     nerror = 0
     DO i = 1,len
        ENDWB(i) = WWW(i,1) * ZDEPTH(i,1)                         &
             + WWW(i,2) * ZDEPTH(i,2)                               &
             + WWW(i,3) * ZDEPTH(i,3)                               &
             + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2)    &
             - (PPL(i)+PPC(i))*0.001*dtt                            &
             + ETMASS(i)*0.001 + ROFF(i)
        ERRORW(i)= TOTWB(i) - ENDWB(i)
        pcmeter(i) = ppc(i) * 0.001*dtt
        plmeter(i) = ppl(i) * 0.001*dtt
        EMETER(i)= ETMASS(i) * 0.001
        !itb...trying a different error check, 1.e-6 is in the
        !itb...noise in 32-bit arithmetic, works fine in 64-bit
        !         if(abs(errorw(i)).gt.1.e-6) then
        IF(ABS(errorw(i)).GT.1.e-5*totwb(i))THEN
           nerror = nerror + 1
           indxerr(nerror) = i
        ENDIF
     ENDDO

     DO j = 1,nerror
        i = indxerr(j)
        igp = i+ioffset
        WRITE(6,900) tau,IGP, TOTWB(i), ENDWB(i), ERRORW(i),    &
             WWW(i,1), WWW(i,2), WWW(i,3),                 &
             CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2), &
             pcmeter(i),plmeter(i), EMETER(i), ROFF(i)
     ENDDO

     nerror = 0
     DO i = 1,len
        CBAL(i) = RADT(i,1) - CHF(i) -                           &
             (ECT(i)+HC(i)+ECI(i) )/DTT
        GBAL(i) = RADT(i,2) - SHF(i) - (EGS(i)+HG(i)+EGI(i)         &
             - HEATEN(i) )/DTT
        ZLHS(i) = RADT(i,1)+RADT(i,2) - CHF(i) - SHF(i)
        ZRHS(i) = HFLUX(i) + (ECT(i) + ECI(i) + EGI(i) + EGS(i)    &
             + HEATEN(i) ) /DTT

        ERRORE(i)= ZLHS(i) - ZRHS(i)
        IF(ABS(errore(i)).GT.1.0) THEN
           nerror = nerror + 1
           indxerr(nerror) = i
        ENDIF
     ENDDO

     DO j = 1,nerror
        i = indxerr(j)
        tm = thm(i) * (0.001*ps(i))**kapa
        igp = i+ioffset
        WRITE(6,910) tau,IGP, ZLHS(i), ZRHS(i),            &
             RADT(i,1), RADT(i,2),CHF(i), SHF(i),             &
             HFLUX(i), ECT(i), ECI(i), EGI(i), EGS(i),        &
             tm,tc(i),tg(i),tgs(i)
        WRITE(6,911)(td(i,n),n=nsoil,1,-1)
        WRITE(6,912) HC(i), HG(i), HEATEN(i), CBAL(i), GBAL(i)
        WRITE(6,901)  TOTWB(i), ENDWB(i), ERRORW(i),       &
             WWW(i,1), WWW(i,2), WWW(i,3),            &
             CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2),   &
             pcmeter(i),plmeter(i), EMETER(i), ROFF(i)
     ENDDO
  ENDIF

900 FORMAT(//,10X,'** WARNING: WATER BALANCE VIOLATION **  ',//,    &
         /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,                  &
         /,1X,'BEGIN, END, DIFF  ', 2(F10.7,1X),g13.5,                   &
         /,1X,'WWW1-3            ', 3(F10.7,1X),                         &
         /,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),                         &
         /,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )
910 FORMAT(//,10X,'** WARNING: ENERGY BALANCE VIOLATION **',//,    &
         /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,                 &
         /,1X,'RHS, LHS              ', 2G13.5,                        &
         /,1X,'RN1, RN2, CHF, SHF, H ', 5G13.5,                        &
         /,1X,'ECT, ECI, EGI, EGS    ', 4G13.5,                        &
         /,1X,'TM,  TC,  TG, TGS     ', 4G13.5,                        &
         /,1X,'HC        HG          ',  G12.5, 12X, G12.5,            &
         /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
911 FORMAT(1x,'TD ', 5G13.5)
912 FORMAT(1X,'HC        HG          ',  G12.5, 12X, G12.5,      &
         /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
901 FORMAT(10X,'WATER BALANCE'                            &
         /,1X,'BEGIN, END        ', 2(F10.7,1X),          &
         /,1X,'ERRORW            ', 1(g13.5,1X),          &
         /,1X,'WWW1-3            ', 3(F10.7,1X),          &
         /,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),          &
         /,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )

return
END SUBROUTINE balan

!##############################################################################
Subroutine RBRD (tc,rbc,zlt,z2,u2,rd,rb,ta,g,rdc,tgs,len)

implicit none

  !========================================================================
  !      CALCULATION OF RB AND RD AS FUNCS OF U2 AND TEMPERATURES
  !========================================================================
  !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
  !
  !       RB (GRB)       CANOPY TO CAS AERODYNAMIC RESISTANCE (S M-1)
  !       RD (GRD)       GROUND TO CAS AERODYNAMIC RESISTANCE (S M-1)
  !       TA (GTA)       CAS TEMPERATURE (K)
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  INTEGER len
  REAL tc(len),rbc(len),zlt(len),              &
      z2(len),u2(len),rd(len),rb(len),        &
      ta(len),rdc(len),tgs(len),              &
      temdif(len), g
  INTEGER i
  REAL fac, fih, d1(len)

  !-----------------------------------------------------------------------
  !      RB       (RB)       : EQUATION (A9), SE-86
  !-----------------------------------------------------------------------

  DO i = 1,len                      !  loop over gridpoint
     TEMDIF(i) = MAX( 0.1E0,  TC(i)-TA(i) )
     FAC = ZLT(i)/890.* (TEMDIF(i)*20.0)**0.25
     RB(i)  = 1.0/(SQRT(U2(i))/RBC(i)+FAC)

     !-----------------------------------------------------------------------
     !      RD       (RD)       : EQUATION (A15), SE-86
     !-----------------------------------------------------------------------

     TEMDIF(i) = MAX( 0.1E0, TGs(i)-TA(i) )
     FIH = SQRT( 1.+9.*G*TEMDIF(i)*Z2(i)/(TGS(i)*U2(i)*U2(i)) )
     RD(i)  = RDC(i) / (U2(i) * FIH)

  ENDDO

return
END SUBROUTINE RBRD
