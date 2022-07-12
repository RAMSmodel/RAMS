!##############################################################################
Subroutine sfcdata ()

use mem_grid
use mem_leaf
use mem_sib
use leaf_coms

implicit none

integer :: k,nnn

!Variables for LEAF3
real, dimension     (nstyp) :: xsand,xclay,xorgan,xrobulk
real, dimension   (8,nstyp) :: soilparms
real, dimension (12,0:nvtyp) :: bioparms

!Variables for SiB
real,dimension(nstyp_sib) :: sbee,sphsat,ssatco,sporos,sslope &
       ,swopt,sskew,srspsat
real,dimension(nvtyp_sib) :: sz1,sz2,sfvcover,schil,ssodep,srootd,sphihalf &
       ,stran11,stran12,stran21,stran22,sref11,sref12     &
       ,sref21,sref22,svmax0,seffcon,sgslope,sgsmin       &
       ,satheta,sbtheta,strda,strdm,strop,srespcp,sslti   & 
       ,sshti,shlti,shhti,ssoref1,ssoref2
character(len=strl1) :: cname

!LEAF soil classes do not match the order of classes in SiB so we need
!to translate these correctly. LEAF uses the classes from Clapp & Hornberger
!and SiB uses modified USDA classes. LEAF use of FAO data only assigns 
!classes 2-8 (Clapp & Hornberger) via routine "datp_datsoil". Others
!default to sandy-clay-loam.
INTEGER, DIMENSION(1:12) :: soil_type_map
DATA soil_type_map /1,2,3,4,6,7,10,9,8,11,12,7/

!******************************************************************************
!******** DATA TABLES FOR THE LEAF3 LAND SURFACE MODEL ************************
!******************************************************************************
!  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke,
!                        1981; Pielke, 1984; Tremback & Kessler, 1985)
!
!  slpots  - saturation moisture potential (m)
!  slmsts  - saturation volumetric moisture content (m3/m3)
!  slbs    - b exponent (dimensionless)
!  slcons  - saturation soil hydraulic conductivity (m/s)
!  slcons0 - surface value for slcons (m/s)
!  slcpd   - dry soil volumetric heat capacity (J/m3/K)
!  slden   - dry soil density (kg/m3) (also total soil porosity)

data soilparms/  &
!-----------------------------------------------------------------------------
!slpots        slbs          slcons0         slden       USDA SOIL CLASS
!      slmsts         slcons          slcpd              # AND NAME
!-----------------------------------------------------------------------------
 -.121, .395,  4.05, .18e-3, .50e-3, 1465.e3, 1600.,.135  & !  1 sand
,-.090, .410,  4.38, .16e-3, .60e-3, 1407.e3, 1600.,.150  & !  2 loamy sand
,-.218, .435,  4.9 , .34e-4, .77e-3, 1344.e3, 1600.,.195  & !  3 sandy loam
,-.786, .485,  5.3 , .72e-5, .11e-4, 1273.e3, 1600.,.255  & !  4 silt loam
,-.478, .451,  5.39, .69e-5, .22e-2, 1214.e3, 1600.,.240  & !  5 loam
,-.299, .420,  7.12, .63e-5, .15e-2, 1177.e3, 1600.,.255  & !  6 sandy clay loam
,-.356, .477,  7.75, .17e-5, .11e-3, 1319.e3, 1600.,.322  & !  7 silty clay loam
,-.630, .476,  8.52, .24e-5, .22e-2, 1227.e3, 1600.,.325  & !  8 clay loam
,-.153, .426, 10.4 , .22e-5, .22e-5, 1177.e3, 1600.,.310  & !  9 sandy clay
,-.490, .492, 10.4 , .10e-5, .10e-5, 1151.e3, 1600.,.370  & ! 10 silty clay
,-.405, .482, 11.4 , .13e-5, .13e-5, 1088.e3, 1600.,.367  & ! 11 clay 
,-.356, .863,  7.75, .80e-5, .80e-5,  874.e3,  300.,.535/   ! 12 peat

data  xsand  /.97,.92,.80,.57,.60,.65,.35,.48,.50,.30,.25,.20/
data  xclay  /.03,.07,.18,.40,.35,.31,.59,.45,.42,.61,.65,.20/
data  xorgan /.00,.01,.02,.03,.05,.04,.06,.07,.08,.09,.10,.60/

data  xrobulk/1200.,1250.,1300.,1400.,1350.,1350.  &
             ,1500.,1450.,1450.,1650.,1700., 500./

!         LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

data bioparms/  &
!-----------------------------------------------------------------------------
!albv_green     sr_max         veg_clump       rootdep             LEAF-3 CLASS #
!     albv_brown     tai_max        veg_frac        dead_frac      AND DESCRIPTION
!          emisv          sai            veg_ht         rcmin
!-----------------------------------------------------------------------------
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  0  Ocean
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  1  Lakes, rivers, streams
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  2  Ice cap/glacier
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  3  Desert, bare soil
 .14, .24, .97, 5.4, 8.0, 1.0, 1.0, .80, 20.0, 1.5, .0, 500., & !  4  Evergreen needleleaf tree
 .14, .24, .95, 5.4, 8.0, 1.0, 1.0, .80, 22.0, 1.5, .0, 500., & !  5  Deciduous needleleaf tree
 .20, .24, .95, 6.2, 7.0, 1.0,  .0, .80, 22.0, 1.5, .0, 500., & !  6  Deciduous broadleaf tree
 .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 1.5, .0, 500., & !  7  Evergreen broadleaf tree
 .21, .43, .96, 5.1, 4.0, 1.0,  .0, .75,   .3,  .7, .7, 100., & !  8  Short grass
 .24, .43, .96, 5.1, 5.0, 1.0,  .0, .80,  1.2, 1.0, .7, 100., & !  9  Tall grass
 .24, .24, .96, 5.1, 1.0,  .2, 1.0, .20,   .7, 1.0, .0, 500., & ! 10  Semi-desert
 .20, .24, .95, 5.1, 4.5,  .5, 1.0, .60,   .2, 1.0, .0,  50., & ! 11  Tundra
 .14, .24, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 12  Evergreen shrub
 .20, .28, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 13  Deciduous shrub
 .16, .24, .96, 6.2, 7.0, 1.0,  .5, .80, 22.0, 1.5, .0, 500., & ! 14  Mixed woodland
 .22, .40, .95, 5.1, 5.0,  .5,  .0, .85,  1.0, 1.0, .0, 100., & ! 15  Crop/mixed farming, C3 grassland
 .18, .40, .95, 5.1, 5.0,  .5,  .0, .80,  1.1, 1.0, .0, 500., & ! 16  Irrigated crop
 .12, .43, .98, 5.1, 7.0, 1.0,  .0, .80,  1.6, 1.0, .0, 500., & ! 17  Bog or marsh
 .20, .36, .96, 5.1, 6.0, 1.0,  .0, .80,  7.0, 1.0, .0, 100., & ! 18  Wooded grassland 
 .20, .36, .90, 5.1, 3.6, 1.0,  .0, .74,  6.0,  .8, .0, 500., & ! 19  Urban and built up
 .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 1.5, .0, 500./   ! 20  Wetland evergreen broadleaf tree

!******************************************************************************
!******** DATA TABLES FOR THE SIB2.5 LAND SURFACE MODEL ***********************
!******************************************************************************
  !...................SOIL CLASSES....................
  !
  !   The soil texture catagories are (based on the 12 USDA classes):
  !   
  !   Soil  Name            % clay  % sand
  !    1    sand               3       92
  !    2    loamy sand         5       82
  !    3    sandy loam         10      65
  !    4    silt loam          13      22
  !    5    silt               7       7
  !    6    loam               18      42
  !    7    sandy clay loam    28      58
  !    8    sandy clay         40      52
  !    9    clay loam          39      32
  !   10    silty clay loam    39      10
  !   11    silty clay         41      7
  !   12    clay               65      19
  !   
  !   Soil properties based on the approximate centroid of the soil texture 
  !   catagory within the USDA texture triangle.
  !   
  !   Modifications:
  !     Kevin Schaefer added respiration variables (Wopt, skew, RespSat) 
  !       using curve fits based on data from Raich et al., 1991 (6/19/00)
  !     Kevin Schaefer updated table using centroid % clay/sand (3/30/01)
  !   
  !   Variables:   
  !   
  !    Bee      : Clapp & Hornberge 'B' exponent
  !    phsat    : Soil tension at saturation (units)
  !    satco    : soil tension at 1/2 assimilation value (true?) units?
  !    poros    : porosity
  !    slope    : slope
  !    wopt     : optimum saturation percentage for respiration
  !    skew     : Power coeff for moisture effect on soil respiration
  !    respsat  : respiration at soil water saturation?

  data sbee/3.387,3.705,4.500,4.977,4.023,5.772,7.362,9.270,9.111,9.111 &
       ,9.429,13.245/

  data sphsat/-0.047,-0.064,-0.107,-0.391,-0.614,-0.214,-0.132,-0.158   &
       ,-0.289,-0.561,-0.614,-0.428/

  data ssatco/0.236E-4,0.166E-4,0.910E-5,0.200E-5,0.118E-5,0.405E-5     &
       ,0.711E-5,0.576E-5,0.285E-5,0.131E-5,0.118E-5,0.180E-5/

  data sporos/0.373,0.386,0.407,0.461,0.480,0.436,0.416,0.423,0.449     &
       ,0.476,0.480,0.465/

  data sslope/0.176,0.176,0.176,0.176,0.176,0.176,0.176,0.176,0.176     &
       ,0.176,0.176,0.176/

  data swopt/59.653,60.080,61.120,61.725,60.501,62.701,64.533,66.520    &
       ,66.363,66.363,66.675,69.920/

  data sskew/0.354,0.357,0.362,0.363,0.360,0.359,0.328,0.232,0.243      &
       ,0.243,0.221,-0.255/

  data srspsat/0.508,0.513,0.525,0.533,0.518,0.545,0.570,0.600,0.598    &
       ,0.598,0.603,0.663/

  !..................VEGETATION/BIOME TYPES............................
  !     The vegetation types are:
  !     # type  Name
  !     1  C3  Tall Broadleaf-Evergreen Trees
  !                  Ref: Stanford Group, Sellers et al. (1989)
  !     2  C3  Tall Broadleaf-Deciduous Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     3  C3  Tall Broadleaf and Needleleaf Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     4  C3  Tall Needleleaf Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Strebel et al. (1982)
  !     5  C3  Tall Needleleaf-DECIDUOUS Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Strebel et al. (1982)
  !     6  C4  Short Vegetation, Same as Types 6, 7, 
  !                       8, and 11 (Stanford-Carnegie)
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     7  C4  Short Vegetation: Ground Cover (Maize Optical Properties)
  !                  Ref: Klink and Willmott (1985), 
  !                       Miller (1972), Sellers (PC*)
  !     8  C4  Short Vegetation: Ground Cover (Maize Optical Properties)
  !                  Ref: Klink and Willmott (1985), 
  !                       Miller (1972), Sellers (PC*)
  !     9  C3  Short Broadleaf Shrubs with Bare Soil
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972),
  !                       Sellers (PC)
  !     10 C3  Short Ground Cover (Tundra)
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972),
  !                       Sellers (PC)
  !     11 C4  No Vegetation (Low Latitude Desert)
  !                  Ref: Sellers (PC) and DORMAN
  !     12 C3  Agriculture (Wheat) and C3 Grasslands
  !                  Ref: Sellers and Dorman (1987), 
  !                       Turner (1974), and Dorman
  !     13 C4  ice
  !                  Ref: Sellers (PC) and Dorman
  !     
  !     * personal communication
  !     
  !............VARIABLES........................
  !
  !     z2 --- canopy top height (meters)
  !     z1 --- canopy bottom height (meters)
  !     fvcover --- fractional vegetation cover (-)
  !     chil    --- leaf angle distribution factor
  !     sodep   --- total soil depth of 3 soil moist lyrs (meters)
  !     rootd   --- rooting depth (meters)
  !     phi_half --- one-half critical leaf-water potential limit(m)
  !     trans(2,2) --- leaf transmittance
  !                    (1,1) - SW green
  !                    (1,2) - LW green
  !                    (2,1) - SW brown
  !                    (2,2) - LW brown
  !     ref(2,2) --- leaf reflectance
  !                    (1,1) - SW green
  !                    (1,2) - LW green
  !                    (2,1) - SW brown
  !                    (2,2) - LW brown
  !     vmax0    - rubisco velocity of sun leaf (mol/m2/s)
  !     effcon   - quantum efficiency (mol/mol)
  !     gsslope  - conductance-photosynthesis slope parm (mol/m2/s)
  !     gsmin    - conductance-photosynthesis intercept (mol/m2/s)
  !     atheta   - wc,we coupling parameter
  !     btheta   - wp,ws coupling parameter
  !     trda     - slope of high temp inhibition function (leaf resp,1/K)
  !     trdm     - half point of high temp inhibition function (leaf resp,K)
  !     trop     - temperature coefficient in GS-A model (K)
  !     respcp   - respiration fraction of Vmax
  !     slti     - slope of low temperature inhibition (1/K)
  !     shti     - slope of high temperature inhibition (1/K)
  !     hlti     - half piont of low temp inhibition (K)
  !     hhti     - half point of high temp inhibition (K)
  !     soref(2) --- soil reflectance
  !                     (1) - visible
  !                     (2) - nir     

  data sz2/35.00,20.00,20.00,17.00,17.00,1.00,1.00,1.00,0.500       &
       ,0.600,1.00,1.00,1.00/

  data sz1/1.00,11.500,10.00,8.500,8.500,0.100,0.100,0.100          &
       ,0.100,0.100,0.100,0.100,0.100/

  data sfvcover/0.874,0.597,0.727,0.558,0.670,0.776,0.343,0.343     &
       ,0.136,0.402,0.055,0.553,0.055/

  data schil/0.100,0.250,0.125,0.010,0.010,-0.300,-0.300,-0.300     &
       ,0.010,0.200,-0.300,-0.300,-0.300/

  data ssodep/3.500,2.000,2.000,2.000,2.000,1.500,1.500,1.500       &
       ,1.500,1.500,1.500,1.500,1.500/

  data srootd/1.500,1.500,1.500,1.500,1.000,1.000,1.000,1.000       &
       ,1.000,1.000,1.000,1.000,1.000/

  data sphihalf/-200.0,-200.0,-200.0,-200.0,-200.0,-200.0           &
       ,-200.0,-200.0,-200.0,-200.0,-200.0,-200.0,-200.0/

  data stran11/0.050,0.050,0.050,0.050,0.050,0.070,0.070,0.070      &
       ,0.050,0.070,0.070,0.070,0.070/

  data stran12/0.250,0.250,0.150,0.100,0.100,0.248,0.248,0.248      &
       ,0.250,0.248,0.248,0.248,0.248/

  data stran21/0.001,0.001,0.001,0.001,0.001,0.220,0.220,0.220      &
       ,0.001,0.220,0.220,0.220,0.220/

  data stran22/0.001,0.001,0.001,0.001,0.001,0.375,0.375,0.375      &
       ,0.001,0.375,0.375,0.375,0.375/

  data sref11/0.060,0.070,0.070,0.080,0.080,0.060,0.100,0.120       &
       ,0.120,0.120,0.120,0.080,0.140/

  data sref12/0.390,0.390,0.380,0.370,0.370,0.400,0.400,0.400       &
       ,0.400,0.400,0.400,0.400,0.420/

  data sref21/0.160,0.160,0.160,0.160,0.160,0.160,0.220,0.220       &
       ,0.220,0.220,0.220,0.160,0.220/

  data sref22/0.430,0.430,0.420,0.410,0.410,0.440,0.480,0.480       &
       ,0.480,0.480,0.480,0.480,0.320/

  data svmax0/0.100E-3,0.100E-3,0.750E-4,0.600E-4,0.100E-3          &
       ,0.300E-4,0.300E-4,0.300E-4,0.600E-4,0.600E-4          &
       ,0.300E-4,0.100E-3,0.300E-4/

  data seffcon/0.080,0.080,0.080,0.080,0.080,0.050,0.050,0.050      &
       ,0.080,0.080,0.050,0.080,0.050/

  data sgslope/9.000,9.000,9.000,9.000,9.000,4.000,4.000,4.000      &
       ,9.000,9.000,4.000,9.000,4.000/

  data sgsmin/0.010,0.010,0.010,0.010,0.010,0.040,0.040,0.040       &
       ,0.010,0.010,0.040,0.010,0.040/

  data satheta/0.980,0.980,0.980,0.980,0.980,0.800,0.800,0.800      &
       ,0.980,0.980,0.800,0.980,0.800/

  data sbtheta/0.950,0.950,0.950,0.950,0.950,0.950,0.950,0.950      &
       ,0.950,0.950,0.950,0.950,0.950/

  data strda/1.300,1.300,1.300,1.300,1.300,1.300,1.300,1.300        &
       ,1.300,1.300,1.300,1.300,1.300/

  data strdm/328.16,328.16,328.16,328.16,328.16,328.16,328.16       &
       ,328.16,328.16,328.16,328.16,328.16,328.16/

  data strop/298.16,298.16,298.16,298.16,298.16,298.16,298.16       &
       ,298.16,298.16,298.16,298.16,298.16,298.16/

  data srespcp/0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015      &
       ,0.015,0.015,0.015,0.015,0.015/

  data sslti/0.200,0.200,0.200,0.200,0.200,0.300,0.300,0.300        &
       ,0.200,0.200,0.300,0.200,0.300/

  data shlti/288.16,283.16,280.16,278.16,278.16,288.16,288.16       &
       ,288.16,283.16,278.16,288.16,281.16,288.16/

  data sshti/0.300,0.300,0.300,0.300,0.300,0.300,0.300,0.300        &
       ,0.300,0.300,0.300,0.300,0.300/

  data shhti/313.16,311.16,307.16,303.16,303.16,313.16,313.16       &
       ,313.16,313.16,303.16,313.16,308.16,313.16/

  data ssoref1/0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.150      &
       ,0.300,0.110,0.300,0.100,0.300/

  data ssoref2/0.225,0.225,0.225,0.225,0.225,0.225,0.225,0.250      &
       ,0.350,0.225,0.350,0.150,0.350/
!******************************************************************************
!******** END OF SOIL AND VEGETATION CONSTANTS ******************************** 
!******************************************************************************

!Set up constants for LEAF3 and SiB. LEAF3 always runs water patch. SiB uses
! LEAF3 initialization, vegetation info, NDVI info, etc. SiB will run on
! vegetation patches only. If running SiB, use SiB soil porosity values to
! keep soil_water consistent between models, otherwise SiB can get soil moisture
! values > 1.0.
if(isfcl <= 2) then

  ! Soil constants
  slz(nzg+1) = 0.

 ! soil type ranges from 1 to nstyp
  do nnn = 1,nstyp
   do k = 1,nzg
      slcons1(k,nnn) = soilparms(4,nnn) ! ORIGINAL form - const with depth
   enddo
   slpots(nnn) = soilparms(1,nnn)
   if(isfcl<=1) slmsts(nnn) = soilparms(2,nnn)
   if(isfcl==2) slmsts(nnn) = sporos(soil_type_map(INT(nnn)))
   slbs(nnn)   = soilparms(3,nnn)
   slcpd(nnn)  = soilparms(6,nnn)
   sfldcap(nnn)= soilparms(8,nnn)
   emisg(nnn)  = .98
   soilcp(nnn) = 0.1 - 0.07 * xsand(nnn)
  enddo

  ! veg type ranges from 0 to nvtyp
  do nnn = 0,nvtyp
   albv_green(nnn) = bioparms(1,nnn)
   albv_brown(nnn) = bioparms(2,nnn)
   emisv(nnn)      = bioparms(3,nnn)
   sr_max(nnn)     = bioparms(4,nnn)
   tai_max(nnn)    = bioparms(5,nnn)
   sai(nnn)        = bioparms(6,nnn)
   veg_clump(nnn)  = bioparms(7,nnn)
   veg_frac(nnn)   = bioparms(8,nnn)
   veg_ht(nnn)     = bioparms(9,nnn)
   dead_frac(nnn)  = bioparms(11,nnn)
   rcmin(nnn)      = bioparms(12,nnn)
   glai_max(nnn)   = tai_max(nnn) - sai(nnn)

   kroot(nnn)   = nzg
   do k = nzg-1,1,-1
      if (slz(k+1) .gt. -bioparms(10,nnn)) kroot(nnn) = k
   enddo
  enddo
endif !If running LEAF3 or SiB

!Set up constants for SiB vegetation model
if(isfcl==2) then

  !Soil properties based on the approximate centroid of the soil texture 
  !catagory within the USDA texture triangle.
  do nnn=1,nstyp_sib
     bee_sib(nnn)    = sbee(nnn)    !Clapp & Hornberge 'B' exponent
     phsat_sib(nnn)  = sphsat(nnn)  !Soil tension at saturation (units)
     satco_sib(nnn)  = ssatco(nnn)  !soil tension at 1/2 assimilation value (units)?
     poros_sib(nnn)  = sporos(nnn)  !porosity
     slope_sib(nnn)  = sslope(nnn)  !slope
     wopt_sib(nnn)   = swopt(nnn)   !Factor coeff for moisture effect on soil respiration
     skew_sib(nnn)   = sskew(nnn)   !Power coeff for moisture effect on soil respiration
     respsat_sib(nnn)= srspsat(nnn) !respiration at soil water saturation?
  enddo

  do nnn=1,nvtyp_sib
     !Vegetation BC'S
     z2_sib(nnn)       = sz2(nnn)      !canopy top height (meters)
     z1_sib(nnn)       = sz1(nnn)      !canopy bottom height (meters)
     fvcover_sib(nnn)  = sfvcover(nnn) !fractional vegetation cover
     chil_sib(nnn)     = schil(nnn)    !leaf angle distribution factor
     sodep_sib(nnn)    = ssodep(nnn)   !total soil depth (meters)
     rootd_sib(nnn)    = srootd(nnn)   !rooting depth (meters)
     phc_sib(nnn)      = sphihalf(nnn) !one-half critical leaf-water potential limit(m)

     !Leaf transmittance and reflectance
     tran_sib(nnn,1,1) = stran11(nnn)
     tran_sib(nnn,1,2) = stran12(nnn)
     tran_sib(nnn,2,1) = stran21(nnn)
     tran_sib(nnn,2,2) = stran22(nnn)
     ref_sib(nnn,1,1)  = sref11(nnn)
     ref_sib(nnn,1,2)  = sref12(nnn)
     ref_sib(nnn,2,1)  = sref21(nnn)
     ref_sib(nnn,2,2)  = sref22(nnn)

     vmax0_sib(nnn)    = svmax0(nnn) !rubisco velocity of sun leaf (mol/m2/s)
     effcon_sib(nnn)   = seffcon(nnn)!quantum efficiency (mol/mol)
     gslope_sib(nnn)   = sgslope(nnn)!conductance-photosynthesis slope parm (mol/m2/s)
     gsmin_sib(nnn)    = sgsmin(nnn) !conductance-photosynthesis intercept (mol/m2/s)
     atheta_sib(nnn)   = satheta(nnn)!wc,we coupling parameter
     btheta_sib(nnn)   = sbtheta(nnn)!wp,ws coupling parameter
     trda_sib(nnn)     = strda(nnn)  !slope of high temp inhibition (leaf resp,1/K)
     trdm_sib(nnn)     = strdm(nnn)  !half point of high temp inhibition (leaf resp,K)
     trop_sib(nnn)     = strop(nnn)  !temperature coefficient in GS-A model (K)
     respcp_sib(nnn)   = srespcp(nnn)!respiration fraction of Vmax
     slti_sib(nnn)     = sslti(nnn)  !slope of low temperature inhibition (1/K)
     shti_sib(nnn)     = sshti(nnn)  !slope of high temperature inhibition (1/K)
     hlti_sib(nnn)     = shlti(nnn)  !half piont of low temp inhibition (K)
     hhti_sib(nnn)     = shhti(nnn)  !half point of high temp inhibition (K)

     !Soil reflectance
     soref_sib(nnn,1)  = ssoref1(nnn)
     soref_sib(nnn,2)  = ssoref2(nnn)
  enddo

  !Now need to read in data from the SiB-Veg-Morphology.txt file.
  !this data will change with different ndvi/soil/biome 
  !data, so it cannot be hardwired in. Contact Ian Baker
  !baker@atmos.colostate.edu with any questions.
  cname=sibfile(1:len_trim(sibfile))
  open(unit=34,file=cname,form='formatted',status='old')
  read(34,*)laig_sib   ! LAI values for lookup table (50)
  read(34,*)fvcg_sib   ! fVCover values for interpolation table (50)
  read(34,*)a_zo_sib   ! Canopy roughness coeff (nvtyp_sib,50,50)
  read(34,*)a_zp_sib   ! Zero plane displacement (nvtyp_sib,50,50)
  read(34,*)a_rbc_sib  ! RB coeff cc1 corrected for snow covered canopy (50x50)
  read(34,*)a_rdc_sib  ! RC coeff cc2 corrected for snow covered canopy (50x50)
  read(34,*)zc_sib     ! Canopy inflection height(m) (nvtyp_sib)
  read(34,*)zlw_sib    ! Leaf Width (nvtyp_sib)
  read(34,*)zlen_sib   ! Leaf Length (nvtyp_sib)
  read(34,*)ltmax_sib  ! Max LAI (nvtyp_sib)
  read(34,*)stem_sib   ! Stem area index (nvtyp_sib)
  read(34,*)nd98_sib   ! Max NDVI (98%tile) (nvtyp_sib)
  read(34,*)nd02_sib   ! Min NDVI (2%tile) (nvtyp_sib)
  read(34,*)srmax_sib  ! Maximum simple ratio (nvtyp_sib)
  read(34,*)srmin_sib  ! Minimum simple ratio (nvtyp_sib)
  close(34)

 if(iprntstmt>=1 .and. print_msg)then
  !Print morphology table values to verify correct loading
  print*,'SIB2 Morphology tables'
  print*,'laig:',laig_sib
  print*,'fvcg:',fvcg_sib
  do nnn=1,50
    print*,'z0:',nnn,a_zo_sib(6,25,nnn)
    print*,'zp:',nnn,a_zp_sib(6,25,nnn)
    print*,'RBC:',nnn,a_rbc_sib(6,25,nnn)
    print*,'RDC:',nnn,a_rdc_sib(6,25,nnn)
  enddo
  print*,'zc:',zc_sib
  print*,'zlw:',zlw_sib
  print*,'zlen:',zlen_sib
  print*,'ltmax:',ltmax_sib
  print*,'stem:',stem_sib
  print*,'nd98:',nd98_sib
  print*,'nd02:',nd02_sib
  print*,'srmax:',srmax_sib
  print*,'srmin:',srmin_sib

  !Print vegetation table values to verify correct loading
  print*,'SIB2 Vegetation tables'
  do nnn=1,nvtyp_sib
     print*,nnn,'z2_sib',z2_sib(nnn)       
     print*,nnn,'z1_sib',z1_sib(nnn)       
     print*,nnn,'fvcover_sib',fvcover_sib(nnn)  
     print*,nnn,'chil_sib',chil_sib(nnn)     
     print*,nnn,'sodep_sib',sodep_sib(nnn)    
     print*,nnn,'rootd_sib',rootd_sib(nnn)    
     print*,nnn,'phc_sib',phc_sib(nnn)      
     print*,nnn,'tran_sib-1-1',tran_sib(nnn,1,1) 
     print*,nnn,'tran_sib-1-2',tran_sib(nnn,1,2) 
     print*,nnn,'tran_sib-2-1',tran_sib(nnn,2,1) 
     print*,nnn,'tran_sib-2-2',tran_sib(nnn,2,2) 
     print*,nnn,'ref_sib-1-1',ref_sib(nnn,1,1)  
     print*,nnn,'ref_sib-1-2',ref_sib(nnn,1,2)  
     print*,nnn,'ref_sib-2-1',ref_sib(nnn,2,1)  
     print*,nnn,'ref_sib-2-2',ref_sib(nnn,2,2)  
     print*,nnn,'vmax0_sib',vmax0_sib(nnn)    
     print*,nnn,'effcon_sib',effcon_sib(nnn)   
     print*,nnn,'gslope_sib',gslope_sib(nnn)   
     print*,nnn,'gsmin_sib',gsmin_sib(nnn)    
     print*,nnn,'atheta_sib',atheta_sib(nnn)   
     print*,nnn,'btheta_sib',btheta_sib(nnn)   
     print*,nnn,'trda_sib',trda_sib(nnn)     
     print*,nnn,'trdm_sib',trdm_sib(nnn)     
     print*,nnn,'trop_sib',trop_sib(nnn)     
     print*,nnn,'respcp_sib',respcp_sib(nnn)   
     print*,nnn,'slti_sib',slti_sib(nnn)     
     print*,nnn,'shti_sib',shti_sib(nnn)     
     print*,nnn,'hlti_sib',hlti_sib(nnn)     
     print*,nnn,'hhti_sib',hhti_sib(nnn)     
     print*,nnn,'soref_sib-1',soref_sib(nnn,1)  
     print*,nnn,'soref_sib-2',soref_sib(nnn,2)  
  enddo
 endif

endif !If running SiB

return
END SUBROUTINE sfcdata

!##############################################################################
Subroutine sfcinit_file (n2,n3,mzg,npat,patch_area,leaf_class,soil_text)

use mem_leaf
use rconstants

implicit none

integer :: n2,n3,mzg,npat,i,j,k,ipat

real, dimension(mzg,n2,n3,npat) :: soil_text
real, dimension(n2,n3,npat) :: patch_area,leaf_class

! This routine fills the arrays PATCH_AREA, leaf_class, and SOIL_TEXT 
! horizontally homogeneously with default values that are defined in the 
! RAMSIN namelist file.  These fields comprise the land/sea surface data 
! types that are normally available on standard RAMS datasets.  The default 
! values assigned here may be overridden by (1) interpolation from coarser 
! grids, (2) specifying new values in routine sfcinit_user in the file
! ruser.f, or (3) reading data from the standard RAMS datasets.

do j = 1,n3
   do i = 1,n2

      patch_area(i,j,1) = 1. - pctlcon         ! patch 1
      leaf_class(i,j,1) = 0.                   ! patch 1

      patch_area(i,j,2) = pctlcon              ! patch 2
      leaf_class(i,j,2) = float(nvgcon)        ! patch 2

      do k = 1,mzg
         soil_text(k,i,j,1) = 0.               ! patch 1
         soil_text(k,i,j,2) = float(nslcon)    ! patch 2
      enddo

   enddo
enddo

do ipat = 3,npat
   do j = 1,n3
      do i = 1,n2

         patch_area(i,j,ipat) = 0.
         leaf_class(i,j,ipat) = leaf_class(i,j,2)

         do k = 1,mzg
            soil_text(k,i,j,ipat) = soil_text(k,i,j,2)
         enddo

      enddo
   enddo
enddo

return
END SUBROUTINE sfcinit_file

!##############################################################################
Subroutine sfcinit_nofile (n1,n2,n3,mzg,mzs,npat,ifm  &
   ,theta,pi0,pp,rv,seatp,seatf                       &
   ,soil_water     ,soil_energy      ,soil_text       &
   ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth  &
   ,veg_fracarea   ,veg_lai          ,veg_tai         &
   ,veg_rough      ,veg_height       ,veg_albedo      &
   ,patch_area     ,patch_rough      ,leaf_class      &
   ,soil_rough     ,sfcwater_nlev    ,stom_resist     &
   ,ground_rsat    ,ground_rvap      ,veg_water       &
   ,veg_temp       ,can_rvap         ,can_temp        & 
   ,veg_ndvip      ,veg_ndvic        ,veg_ndvif       &
   ,snow_mass      ,snow_depth       ,soil_moist_top  & 
   ,soil_moist_bot ,soil_temp_top    ,soil_temp_bot   &
   ,glat           ,glon             ,zot)

use mem_grid
use mem_leaf
use leaf_coms
use io_params
use rconstants
use mem_varinit

implicit none

logical :: there
character(len=7) :: cgrid
character(len=strl1) :: flnm
integer :: n1,n2,n3,mzg,mzs,npat,ifm,i,j,k,ipat &
 ,nveg,nsoil,initsurf,nc,runsoilingest,runsnowingest

real :: c1,airtemp,prsv,piv
real :: tsoil,fice !Saleeby: add these for frozen soil computation
real, dimension(n1,n2,n3) :: theta,pi0,pp,rv
real, dimension(n2,n3)    :: glat,glon,zot  &
                            ,seatp,seatf,snow_mass,snow_depth &
                            ,soil_moist_top,soil_moist_bot &
                            ,soil_temp_top,soil_temp_bot

real, dimension(mzg,n2,n3,npat) :: soil_water,soil_energy,soil_text
real, dimension(mzs,n2,n3,npat) :: sfcwater_mass,sfcwater_energy  &
                                  ,sfcwater_depth

real, dimension(n2,n3,npat) :: veg_fracarea ,veg_lai       ,veg_tai      &
                              ,veg_rough                                 &
                              ,veg_height   ,veg_albedo    ,patch_area   &
                              ,patch_rough  ,leaf_class   &
                              ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                              ,ground_rsat  ,ground_rvap   ,veg_water    &
                              ,veg_temp     ,can_rvap      ,can_temp     &
                              ,veg_ndvip    ,veg_ndvic     ,veg_ndvif

! This routine fills the primary LEAF3 arrays for which standard RAMS
! data files do not exist with initial default values.  Many of the 
! initial values are horizontally homogeneous, although some depend on
! atmospheric conditions.  The default values assigned here may be 
! overridden by (1) specification from coarser grids or (2) specifying new 
! values in routine sfcinit_nofile_user in the file ruser.f.

c1 = .5 * cpi

! Time interpolation factor for updating SST

if (iupdsst == 0) then
   timefac_sst = 0.
else
   timefac_sst = (time - ssttime1(ifm)) / (ssttime2(ifm) - ssttime1(ifm))
endif

!--------------------------------------------------------------------------------
!Saleeby (4-8-2013): Ingest soil moisture, soil temperature, snow depth
!and snow water from varfiles. Soil moisture ingest will only
!work at initialization and not for history restarts, even if adding a grid.
!--------------------------------------------------------------------------------
!Set soil/snow ingest flags to zero
initsurf=0
runsoilingest=0
runsnowingest=0
!Check for initial start for soil/snow ingest
if(trim(runtype) == 'INITIAL' .or. trim(runtype) == 'ERROR') initsurf=1
!Check to see which varfiles are present since these contain the soil data
if(initial==2)then
 write(cgrid,'(a2,i1,a3)') '-g',ifm,'.h5'
 nc=len_trim(fnames_varf(nvarffl))
 flnm=fnames_varf(nvarffl)(1:nc-4)//trim(cgrid)
 inquire(file=trim(flnm),exist=there)
endif
!Set flag for soil ingest
if(isoildat==1 .and. initsurf==1 .and. initial==2 .and. there) runsoilingest=1
!Set flag for snow ingest
if(isnowdat==1 .and. initsurf==1 .and. initial==2 .and. there) runsnowingest=1
!Print indication of using soil and snow from varfiles
if(runsoilingest==1 .or. runsnowingest==1)then
 if(iprntstmt>=1 .and. print_msg)then
  print*,'-------------------------------------------------------'
  print*,'SETTING UP SOIL AND/OR SNOW LEVEL VARIABLES ON GRID:',IFM
  print*,'-------------------------------------------------------'
 endif
endif

!Initialize surface properties
do j = 1,n3
   do i = 1,n2
      piv = c1 * (pi0(1,i,j) + pi0(2,i,j)   &
                      + pp(1,i,j) + pp(2,i,j))
      airtemp = theta(2,i,j) * piv
      prsv = piv ** cpor * p00

      patch_rough(i,j,1) = 0.001
      can_temp(i,j,1) = airtemp
      can_rvap(i,j,1) = rv(2,i,j)

      soil_energy(mzg,i,j,1) = 334000.  &
         + 4186. * (seatp(i,j) + (seatf(i,j) - seatp(i,j))  &
         * timefac_sst - 273.15)

      !Loop over patches
      do ipat = 2,npat

         nveg = nint(leaf_class(i,j,ipat))

         soil_rough(i,j,ipat) = zrough
         patch_rough(i,j,ipat) = max(zrough,zot(i,j))
         veg_height(i,j,ipat) = veg_ht(nveg)
         stom_resist(i,j,ipat) = 1.e6

         veg_temp(i,j,ipat) = airtemp
         can_temp(i,j,ipat) = airtemp

         veg_water(i,j,ipat) = 0.
         can_rvap(i,j,ipat) = rv(2,i,j)

         !Loop over soil layers
         do k = 1,mzg

            nsoil = nint(soil_text(k,i,j,ipat))
            soil_water(k,i,j,ipat) = max(soilcp(nsoil),slmstr(k) * slmsts(nsoil))

            !Saleeby(4-9-2013):If using variable soil initialization from ingest data
            if(runsoilingest==1)then
             if(k==mzg) then
              soil_water(k,i,j,ipat) = max(soilcp(nsoil),soil_moist_top(i,j))
             else
              soil_water(k,i,j,ipat) = max(soilcp(nsoil),soil_moist_bot(i,j))
             endif
            endif

            !For persistent wetlands (bogs, marshes, fens, swamps) and irrigated 
            !crops, initialize with saturated soil. Currently, this corresponds to 
            !leaf classes 16, 17, and 20.
            if (nint(leaf_class(i,j,ipat)) == 16 .or.  &
                nint(leaf_class(i,j,ipat)) == 17 .or.  &
                nint(leaf_class(i,j,ipat)) == 20) then
               soil_water(k,i,j,ipat) = slmsts(nsoil)
            endif

            !Limit soil water fraction maximum for soil type. Need to do this 
            !since input soil water from gridded data can differ from RAMS due 
            !to differences in soil classes and soil porosity.
            soil_water(k,i,j,ipat) = min(slmsts(nsoil),soil_water(k,i,j,ipat))

            !By default, initialize soil internal energy at a temperature equal 
            !to airtemp + stgoff(k) and with all water assumed to be liquid.  If 
            !the temperature is initially below 0C, this will immediately adjust 
            !to soil at 0C with part ice.  In order to begin with partially or 
            !totally frozen soil, reduce or remove the latent-heat-of-fusion 
            !term (the one with the factor of 3.34) from soil_energy below. If 
            !the soil is totally frozen and the temperature is below zero C, the 
            !factor of 4.186 should be changed to 2.093 to reflect the reduced 
            !heat capacity of ice compared to liquid. These changes may be 
            !alternatively be done in routine sfcinit_user in ruser.f

            !Saleeby: Old Method sometimes allows surface runaway cooling!
            !soil_energy(k,i,j,ipat) = (airtemp - 273.15 + stgoff(k))  &
            !   * (slcpd(nsoil) + soil_water(k,i,j,ipat) * 4.186e6)  &
            !   + soil_water(k,i,j,ipat) * 3.34e8

            !*************************************************************
            !Saleeby: New Method allows for tsoil < 0C and set initial soil 
            ! ice fraction (fice) to zero
            tsoil = (airtemp - 273.15 + stgoff(k))

            !Saleeby(4-9-2013):If we are ingesting soil temperature from varfiles
            if(runsoilingest==1)then
             if(k==mzg) then
              tsoil=max(tsoil,soil_temp_top(i,j)-273.15)
             else
              tsoil=max(tsoil,soil_temp_bot(i,j)-273.15)
             endif
            endif

            fice  = 0.0
            !For soil temperature > 0C
            if (tsoil.gt.0.) then
              soil_energy(k,i,j,ipat) = tsoil  &
                 * (slcpd(nsoil) + soil_water(k,i,j,ipat) * 4.186e6)  &
                 + soil_water(k,i,j,ipat) * 3.34e8
            !For soil temperature <= 0C
            elseif (tsoil.le.0.) then
              !This soil ice fraction should be improved
              !Just a linear computation with fice=0.0 at 273K
              ! up to fice=1.0 at 272K
              fice = min(1.0,max(0.0, -1.0*tsoil*0.1 ))
              soil_energy(k,i,j,ipat) = tsoil  &
                 * (slcpd(nsoil) + &
                    (1.-fice)*soil_water(k,i,j,ipat) * 4.186e6 +  &
                         fice*soil_water(k,i,j,ipat) * 2.093e6)  &
                 + (1.-fice)*soil_water(k,i,j,ipat) * 3.34e8
            endif
            !*************************************************************

         enddo !end looping over soil layers

         !Determine number of surface water or snow levels and set values of 
         !mass, depth, and energy for each level present.

         !Set number of inital surface water layers to zero
         sfcwater_nlev(i,j,ipat) = 0.

         !Loop over surface water / snow layers
         do k = 1,mzs

            sfcwater_mass(k,i,j,ipat)   = 0.
            sfcwater_energy(k,i,j,ipat) = 0.
            sfcwater_depth(k,i,j,ipat)  = 0.

            !For persistent wetlands (bogs, marshes, fens, swamps), initialize 
            !with 10 cm water depth. This corresponds to leaf classes 17 and 20.
            if (nint(leaf_class(i,j,ipat)) == 17 .or.  &
                nint(leaf_class(i,j,ipat)) == 20) then
               if (k .eq. 1) then
                  sfcwater_mass(k,i,j,ipat) = 100.
                  sfcwater_energy(k,i,j,ipat) = (airtemp - 193.36) * 4186.
                  sfcwater_depth(k,i,j,ipat) = .1
               endif
            endif

            !Initialize surface water layer quantities
            if (runsnowingest==1) then
             if (snow_mass(i,j) > 0.) then
               if (k .eq. 1) then
                  sfcwater_mass(k,i,j,ipat)   = sfcwater_mass(k,i,j,ipat)   &
                     + snow_mass(i,j)
                  sfcwater_energy(k,i,j,ipat) = sfcwater_energy(k,i,j,ipat) &
                     + min(0., (airtemp - 273.15) * 2093.)
                  if(runsnowingest==1) then 
                    sfcwater_depth(k,i,j,ipat) = sfcwater_depth(k,i,j,ipat) &
                        + snow_depth(i,j) ! From Varfiles
                  else
                    sfcwater_depth(k,i,j,ipat) = sfcwater_depth(k,i,j,ipat)  &
                        + snow_mass(i,j) * 5.e-3   ! 5x equivalent liquid depth
                  endif
               endif
             endif
            endif

            !Determine number of surface water layers
            if (sfcwater_mass(k,i,j,ipat) > 0.) sfcwater_nlev(i,j,ipat) = float(k)

         enddo !end loop over surface water layers

         !Call for all land surface options (isfcl)
         if (ipat >= 2) CALL ndvi (ifm  &
            ,leaf_g(ifm)%leaf_class (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvip  (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvic  (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvif  (i,j,ipat)   )

         !Only call this for LEAF3. SIB will use its own values.
         if (ipat >= 2 .and. isfcl<=1)           &
            CALL veg (ifm                        &
            ,leaf_g(ifm)%leaf_class   (i,j,ipat) &
            ,leaf_g(ifm)%veg_fracarea (i,j,ipat) &
            ,leaf_g(ifm)%veg_lai      (i,j,ipat) &
            ,leaf_g(ifm)%veg_tai      (i,j,ipat) &
            ,leaf_g(ifm)%veg_rough    (i,j,ipat) &
            ,leaf_g(ifm)%veg_height   (i,j,ipat) &
            ,leaf_g(ifm)%veg_albedo   (i,j,ipat) &
            ,leaf_g(ifm)%veg_ndvic    (i,j,ipat) )

         CALL grndvap (mzs  &
            ,soil_energy(mzg,i,j,ipat),soil_water(mzg,i,j,ipat)  &
            ,soil_text  (mzg,i,j,ipat),sfcwater_energy(1,i,j,ipat) &
            ,sfcwater_nlev(i,j,ipat),ground_rsat(i,j,ipat) &
            ,ground_rvap(i,j,ipat),can_rvap(i,j,ipat),prsv,i,j)

      enddo !end loop over patches
   enddo
enddo

return
END SUBROUTINE sfcinit_nofile

!##############################################################################
Subroutine datp_datq (datp,datq)

! This routine maps the input datp classes to a smaller set datq
! which represents the full set of LEAF-2 or LEAF-3 classes for which 
! LSP values are
! defined.

implicit none

integer datq,catb(0:94)
real datp

!  Olson Global Ecosystems dataset (94 classes) mapped to LEAF-3 classes
!  (see leaf3_document).

!-------------------------------------------!
data catb/ 0,                             & !
          19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  0
!!!       10, 2,17, 1, 0,12,13,14,18, 4,  & ! 10 (use for future Olson data)
          10, 2,20, 0, 0,12,13,14,18, 4,  & ! 10 (use for current Olson data)
           4, 4,14,14, 6, 6, 4, 7, 7,15,  & ! 20
          15, 6, 7, 7,15,16,16,16,16, 8,  & ! 30
           8, 8,18,17,17,12,12, 7,10, 3,  & ! 40
          10,10,11,14,18,18,18,18,13, 6,  & ! 50
           5, 4,11,12, 0, 0, 0, 0, 3, 2,  & ! 60
           3,20, 0,17,17,17, 4,14, 7, 3,  & ! 70
           3, 3, 3, 3, 3, 3, 8,12, 7, 6,  & ! 80
          18,15,15,15                  /    ! 90
!-------------------------------------------!
!          1  2  3  4  5  6  7  8  9 10

datq = catb(nint(datp))

return
END SUBROUTINE datp_datq

!##############################################################################
Subroutine datp_datsoil (datp,datsoil)

! This routine maps the input datp soil classes to a smaller set datsoil
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

implicit none

integer datsoil,catb(0:133)
real datp

! (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our 
! own missing value designated 133, to the USDA soil textural classes.  FAO 
! classes [0] (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not 
! used in original FAO dataset), [132] (water), and [133] (our missing value) 
! are all mapped to a default class of sandy clay loam in case they happen to
! correspond to a land surface area in the landuse dataset that RAMS uses to
! define land area.  We wrote missing value class 133 to the RAMS FAO files
! whenever a negative value (which is out of range of defined values) was
! found in the original FAO dataset, which occurred in about 2.6% of the
! pixels.  For the remaining FAO classes, a cross reference table to Zobler 
! soil texture classes that was provided, plus our own cross referencing table
! from Zobler to USDA classes listed below, provides the mapping from FAO to 
! USDA.  In this mapping, we use only organic USDA classes and omit nonorganic
! classes since most soils do contain organic matter, and organic content 
! information is not provided in the Zobler classes.

!  Zobler Class              USDA Class

!  1  Coarse                 2  Loamy sand
!  2  Medium                 4  Silt loam
!  3  Fine                   8  Clay loam
!  4  Coarse-medium          3  Sandy loam
!  5  Coarse-fine            6  Sandy clay loam
!  6  Medium-fine            7  Silty clay loam
!  7  Coarse-medium-fine     6  Sandy clay loam
!  8  Organic matter         5  Loam

!                            1  Sand (not used)
!                            9  Sandy clay (not used)
!                           10  Silty clay (not used)
!                           11  Clay (not used)
!                           12  Peat (not used)

data catb/ 6  &
         , 6, 4, 4, 7, 7, 8, 6, 4, 4, 4  &
         , 7, 4, 4, 4, 8, 4, 8, 4, 4, 8  &
         , 4, 2, 4, 4, 4, 4, 6, 8, 8, 8  &
         , 4, 8, 8, 2, 6, 4, 7, 4, 4, 3  &
         , 4, 6, 7, 4, 4, 4, 4, 4, 4, 4  &
         , 4, 4, 4, 4, 4, 4, 2, 4, 4, 2  &
         , 4, 3, 4, 2, 7, 6, 4, 4, 6, 8  &
         , 8, 7, 2, 5, 4, 5, 6, 6, 4, 2  &
         , 2, 2, 4, 6, 2, 2, 2, 2, 2, 4  &
         , 2, 2, 2, 4, 2, 4, 3, 6, 2, 7  &
         , 4, 4, 4, 8, 8, 8, 3, 7, 4, 4  &
         , 4, 3, 6, 4, 2, 4, 4, 4, 2, 2  &
         , 2, 4, 6, 4, 4, 7, 7, 6, 3, 2  &
         , 2, 6, 6 /

datsoil = catb(nint(datp))

return
END SUBROUTINE datp_datsoil
