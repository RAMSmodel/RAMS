!##############################################################################
Subroutine micro_master ()

use micphys
use rconstants

implicit none

integer :: lhcat,khcat,lcat
integer, dimension(8) :: lcat0

real, dimension(7,16) :: dstprms,dstprms1,dstprms2
real, dimension(16,16) :: jpairr,jpairc

data lcat0 /1,2,3,4,5,6,7,16/ ! lcat corressponding to lhcat

data dstprms1/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!---------------------------------------------------------------------- 
    .5,      524.,     3.,   1.26e7,   1.91,   2.e-6,   40.e-6,  & !cloud
    .5,      524.,     3.,     149.,     .5,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,  5.769e5,   1.88,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,  188.146,   .933,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,    3.084,     .2,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     93.3,     .5,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,     161.,     .5,   .8e-3,   10.e-3,  & !hail 
  .429,     .8854,    2.5,     316.,   1.01,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     316.,   1.01,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,  5.769e5,   1.88,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,   3.19e4,   1.66,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    4.836,    .25,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,    4.836,    .25,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,  188.146,   .933,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,  1348.38,  1.241,      00,       00,  & !snow ros
    .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !drizzle

data dstprms2/ &
!----------------------------------------------------------------------
! shape      cfmas   pwmas      cfvt    pwvt     dmb0      dmb1
!----------------------------------------------------------------------
    .5,      524.,     3.,   1.26e7,   1.91,   2.e-6,   40.e-6,  & !cloud
    .5,      524.,     3.,     144.,   .497,   .1e-3,    5.e-3,  & !rain
  .179,     110.8,   2.91,    1538.,   1.00,  15.e-6,  125.e-6,  & !pris col
  .179,  2.739e-3,   1.74,     27.7,   .484,   .1e-3,   10.e-3,  & !snow col
    .5,      .496,    2.4,     16.1,   .416,   .1e-3,   10.e-3,  & !aggreg
    .5,      157.,     3.,     332.,   .786,   .1e-3,    5.e-3,  & !graup
    .5,      471.,     3.,    152.1,   .497,   .8e-3,   10.e-3,  & !hail
  .429,     .8854,    2.5,   20801.,  1.377,      00,       00,  & !pris hex
 .3183,   .377e-2,     2.,     56.4,   .695,      00,       00,  & !pris den
 .1803,   1.23e-3,    1.8,   1617.9,   .983,      00,       00,  & !pris ndl
    .5,     .1001,  2.256,    6239.,   1.24,      00,       00,  & !pris ros
  .429,     .8854,    2.5,    30.08,   .563,      00,       00,  & !snow hex
 .3183,   .377e-2,     2.,     3.39,   .302,      00,       00,  & !snow den
 .1803,   1.23e-3,    1.8,     44.6,   .522,      00,       00,  & !snow ndl
    .5,     .1001,  2.256,    125.7,   .716,      00,       00,  & !snow ros
    .5,      524.,     3.,   1.26e7,   1.91,  65.e-6,  100.e-6/    !drizzle

data jpairr/  &
     0,  0,  0,  1,  2,  3,  4,  0,  0,  0,  0,  5,  6,  7,  8,  0,  &
     0,  0,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  0,  &
     0, 22, 23, 24,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
    25, 26, 27, 28,  0,  0,  0, 29, 30, 31, 32,  0,  0,  0,  0, 33,  &
    34, 35, 36, 37,  0,  0,  0, 38, 39, 40, 41, 42, 43, 44, 45, 46,  &
    47, 48, 49, 50, 51,  0,  0, 52, 53, 54, 55, 56, 57, 58, 59, 60,  &
    61, 62, 63, 64, 65, 66,  0, 67, 68, 69, 70, 71, 72, 73, 74, 75,  &
     0, 76,  0, 77,  0,  0,  0, 78,  0,  0,  0, 79, 80, 81, 82,  0,  &
     0, 83,  0, 84,  0,  0,  0,  0, 85,  0,  0, 86, 87, 88, 89,  0,  &
     0, 90,  0, 91,  0,  0,  0,  0,  0, 92,  0, 93, 94, 95, 96,  0,  &
     0, 97,  0, 98,  0,  0,  0,  0,  0,  0, 99,100,101,102,103,  0,  &
   104,105,106,  0,  0,  0,  0,107,108,109,110,111,  0,  0,  0,112,  &
   113,114,115,  0,  0,  0,  0,116,117,118,119,  0,120,  0,  0,121,  &
   122,123,124,  0,  0,  0,  0,125,126,127,128,  0,  0,129,  0,130,  &
   131,132,133,  0,  0,  0,  0,134,135,136,137,  0,  0,  0,138,139,  &
     0,  0,  0,140,141,142,143,  0,  0,  0,  0,144,145,146,147,  0/

data jpairc/  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
     4,  5,  6,  7,  0,  0,  0,  8,  9, 10, 11,  0,  0,  0,  0, 12,  &
    13, 14, 15, 16, 17,  0,  0, 18, 19, 20, 21, 22, 23, 24, 25, 26,  &
    27, 28, 29, 30, 31, 32,  0, 33, 34, 35, 36, 37, 38, 39, 40, 41,  &
    42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  &
     0, 58,  0,  0,  0,  0,  0, 59,  0,  0,  0,  0,  0,  0,  0,  0,  &
     0, 60,  0,  0,  0,  0,  0,  0, 61,  0,  0,  0,  0,  0,  0,  0,  &
     0, 62,  0,  0,  0,  0,  0,  0,  0, 63,  0,  0,  0,  0,  0,  0,  &
     0, 64,  0,  0,  0,  0,  0,  0,  0,  0, 65,  0,  0,  0,  0,  0,  &
    66, 67, 68,  0,  0,  0,  0, 69, 70, 71, 72, 73,  0,  0,  0, 74,  &
    75, 76, 77,  0,  0,  0,  0, 78, 79, 80, 81,  0, 82,  0,  0, 83,  &
    84, 85, 86,  0,  0,  0,  0, 87, 88, 89, 90,  0,  0, 91,  0, 92,  &
    93, 94, 95,  0,  0,  0,  0, 96, 97, 98, 99,  0,  0,  0,100,101,  &
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/

!  Define several parameters from above data list

do lhcat=1,nhcat
   !Using original RAMS 4.3 power laws
   if(iplaws==0) then
     dstprms(1,lhcat) = dstprms1(1,lhcat)
     dstprms(2,lhcat) = dstprms1(2,lhcat)
     dstprms(3,lhcat) = dstprms1(3,lhcat)
     dstprms(4,lhcat) = dstprms1(4,lhcat)
     dstprms(5,lhcat) = dstprms1(5,lhcat)
     dstprms(6,lhcat) = dstprms1(6,lhcat)
     dstprms(7,lhcat) = dstprms1(7,lhcat)
   !Using Carver/Mitchell 1996 power laws
   elseif(iplaws==1.or.iplaws==2) then
     dstprms(1,lhcat) = dstprms2(1,lhcat)
     dstprms(2,lhcat) = dstprms2(2,lhcat)
     dstprms(3,lhcat) = dstprms2(3,lhcat)
     dstprms(4,lhcat) = dstprms2(4,lhcat)
     dstprms(5,lhcat) = dstprms2(5,lhcat)
     dstprms(6,lhcat) = dstprms2(6,lhcat)
     dstprms(7,lhcat) = dstprms2(7,lhcat)
   endif

   shape(lhcat) = dstprms(1,lhcat)
   cfmas(lhcat) = dstprms(2,lhcat)
   pwmas(lhcat) = dstprms(3,lhcat)
   cfvt (lhcat) = dstprms(4,lhcat)
   pwvt (lhcat) = dstprms(5,lhcat)

   do khcat=1,nhcat
      ipairc(lhcat,khcat) = jpairc(lhcat,khcat)
      ipairr(lhcat,khcat) = jpairr(lhcat,khcat)
   enddo
enddo

do lcat=1,ncat
   lhcat = lcat0(lcat)
   emb0 (lcat) = cfmas(lhcat) * dstprms(6,lhcat) ** pwmas(lhcat)
   emb1 (lcat) = cfmas(lhcat) * dstprms(7,lhcat) ** pwmas(lhcat)
enddo

if (level .ne. 3) return

! Handle the creating of the main microphysics collection table
CALL mkcoltb ()

return
END SUBROUTINE micro_master

!##############################################################################
Subroutine init_ifn (n1,n2,n3,cifnp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real :: e_a,e_b,e_t,e_k,e_c
real, dimension(n1,n2,n3) :: cifnp,dn0
real :: cin_maxt

!CREATE YOUR OWN CUSTOM POTENTIAL IFN/INP PROFILE BELOW. OUR SORT OF
!DEFAULT IS THE EXPONENTIALLY DECREASING PROFILE, BUT THIS IS BY ALL
!MEANS NOT A DEFINITIVE PROFILE, JUST A NECESSARY PLACEHOLDER SO THAT
!SOME HETEROGENEOUS ICE NUCLEATION CAN OCCUR IF IIFN==1or2. IIFN==3
!USES ONE OF THE DEMOTT FORMULAS AND IS BASED ON LARGE AEROSOLS DETERMINED 
!FROM THE CLOUD DROPLET NUCLEATING AEROSOL DISTRIBUTIONS.

! Initialize IFN
if(iaeroprnt==1 .and. print_msg) then
 print*,'Start Initializing Ice Nuclei concentration'
endif

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Convert RAMSIN #/mg to #/kg
   cin_maxt = cin_max * 1.e6

   !Set up Vertical profile 

   !*************************************************************************** 
   !Exponential decrease that scales with pressure decrease.
   if (cin_max >= 0.0) then
    if(k<=2) cifnp(k,i,j)=cin_maxt
    if(k >2) cifnp(k,i,j)=cin_maxt*exp(-zt(k)/7000.)
    !Output initial sample profile
    if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' Ice Nuclei - init (k,zt,ifn/mg,ifn/L) on Grid:',ifm
     print'(a9,i5,f11.1,2f17.7)',' IFN-init' &
        ,k,zt(k),cifnp(k,i,j)/1.e6,cifnp(k,i,j)/1.e3*dn0(k,i,j)
    endif
   !*************************************************************************** 
   !Use SPICULE INP profile if cin_maxt ~ -1.0
   !Profile from aircraft obs during the SPICULE field campaign.
   !To mimic the SPICULE profile in magnitude, need to set the following parameters.
   !and use IIFN=2 for DeMott scheme, and IFN_FORMULA=2 for DeMott(2015) dust formula.
   !Perhaps using this profile is more realistic given that is has basis
   !in reality from observations of the central-western continental U.S. 
   elseif (cin_max > -1.01 .and. cin_max < -0.99) then
    e_t = 0.10 ! Controls the shape of the profile
    e_a = 2.50 ! approximately cin_maxt
    e_k = 3.80 
    e_b = 5.80 ! Level where you want profile to start decreasing
    e_c = 1.32 ! Necessary so you do not get negative numbers (changes with "a")
    cifnp(k,i,j)=(-1.0*(e_a/2.0)*(erf(((zt(k)/1000)-e_b)/sqrt(4.0*e_k*e_t)))+e_c)
    !NOTE THAT THIS CURVE FIT TO DATA IS BASED ON #/CM3 UNITS, SO CONVERT BELOW
    cifnp(k,i,j) = cifnp(k,i,j) * 1.e6 / dn0(k,i,j) ! convert #/cm3 to #/kg
    !Output initial sample profile
    if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' Ice Nuclei - init (k,zt,inp/cm3,inp/mg,inp/L) on Grid:',ifm
     print'(a9,i5,f11.1,3f12.3)',' IFN-init',k,zt(k) &
       ,cifnp(k,i,j)/1.e6*dn0(k,i,j),cifnp(k,i,j)/1.e6,cifnp(k,i,j)/1.e3*dn0(k,i,j)
    endif
   !*************************************************************************** 
   !Use MC3E INP profile if cin_maxt ~ 2.0
   !Profile from aircraft obs during the MC3E field project. (Marinescu et al. 2016)
   !To mimic the MC3E profile in magnitude, need to set the following parameters.
   !and use IIFN=2 for DeMott scheme, and IFN_FORMULA=1 for DeMott(2010) general formula.
   !Perhaps using this profile is more realistic given that is has basis
   !in reality from observations in the southern great plains of the continental U.S. 
   elseif (cin_max > -2.01 .and. cin_max < -1.99) then
    e_t = 0.35 ! Controls the shape of the profile
    e_a = 3.28 ! approximately cin_maxt
    e_k = 3.00
    e_b = 2.27 ! Level where you want profile to start decreasing
    e_c = 1.96 ! Necessary so you do not get negative numbers (changes with "a")
    cifnp(k,i,j)=(-1.0*(e_a/2.0)*(erf(((zt(k)/1000)-e_b)/sqrt(4.0*e_k*e_t)))+e_c)
    !NOTE THAT THIS CURVE FIT TO DATA IS BASED ON #/MG UNITS, SO CONVERT BELOW
    cifnp(k,i,j) = cifnp(k,i,j) * 1.e6 ! convert #/mg to #/kg
    !Output initial sample profile
    if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' Ice Nuclei - init (k,zt,inp/cm3,inp/mg,inp/L) on Grid:',ifm
     print'(a9,i5,f11.1,3f12.3)',' IFN-init',k,zt(k) &
       ,cifnp(k,i,j)/1.e6*dn0(k,i,j),cifnp(k,i,j)/1.e6,cifnp(k,i,j)/1.e3*dn0(k,i,j)
    endif
   !***************************************************************************
   endif

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_ifn

!##############################################################################
Subroutine init_ccn1 (n1,n2,n3,cn1np,cn1mp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: cn1np,cn1mp,dn0
real :: ccn1_maxt

! Initialize CCN mode 1
if(iaeroprnt==1 .and. print_msg) print*,'Start Initializing CCN mode 1 concen'

!Convert RAMSIN #/mg to #/kg
 ccn1_maxt = ccn1_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up Vertical profile
   if(k<=2) cn1np(k,i,j)=ccn1_maxt
   !Exponential decrease that scales with pressure decrease
   if(k>2)  cn1np(k,i,j)=ccn1_maxt*exp(-zt(k)/7000.)

   !Output initial sample profile
   if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' CCN-1-init (k,zt,ccn1/mg,ccn1/cc) on Grid:',ifm
     print'(a11,i5,f11.1,2f17.7)',' CCN-1-init' &
        ,k,zt(k),cn1np(k,i,j)/1.e6,cn1np(k,i,j)/1.e6*dn0(k,i,j)
   endif

   !Set up Field of CCN-mode-1 mass mixing ratio (kg/kg)
   cn1mp(k,i,j) = ((aero_medrad(1)*aero_rg2rm(1))**3.) &
                *cn1np(k,i,j)/(0.23873/aero_rhosol(1))

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_ccn1

!##############################################################################
Subroutine init_ccn2 (n1,n2,n3,cn2np,cn2mp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: cn2np,cn2mp,dn0
real :: ccn2_maxt

! Initialize CCN mode 2
if(iaeroprnt==1 .and. print_msg) print*,'Start Initializing CCN mode 2 concen'

!Convert RAMSIN #/mg to #/kg
 ccn2_maxt = ccn2_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Set up Vertical profile
   if(k<=2) cn2np(k,i,j)=ccn2_maxt
   ! Exponential decrease that scales with pressure decrease
   if(k>2)  cn2np(k,i,j)=ccn2_maxt*exp(-zt(k)/7000.)

   !Output initial sample profile
   if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' CCN-2-init (k,zt,ccn2/mg,ccn2/cc) on Grid:',ifm
     print'(a11,i5,f11.1,2f17.7)',' CCN-2-init' &
        ,k,zt(k),cn2np(k,i,j)/1.e6,cn2np(k,i,j)/1.e6*dn0(k,i,j)
   endif

   !Set up Field of CCN-mode-2 mass mixing ratio (kg/kg)
   cn2mp(k,i,j) = ((aero_medrad(2)*aero_rg2rm(2))**3.) &
                *cn2np(k,i,j)/(0.23873/aero_rhosol(2))

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_ccn2

!##############################################################################
Subroutine init_dust (n1,n2,n3,md1np,md2np,md1mp,md2mp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: md1np,md2np,md1mp,md2mp,dn0
real :: dust1_maxt,dust2_maxt

! Initialize Dust
if(iaeroprnt==1 .and. print_msg) then
 print*,'Start Initializing DUST concentration'
 print*,'idust,iaerorad',idust,iaerorad
endif

!Convert RAMSIN #/mg to #/kg
 dust1_maxt = dust1_max * 1.e6
 dust2_maxt = dust2_max * 1.e6  

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !If not using dust source model
   if(idust == 1) then
     !Set up concentration of SMALL MODE Mineral Dust (#/kg)
     if(k<=2) md1np(k,i,j)=dust1_maxt
     if(k>2)  md1np(k,i,j)=dust1_maxt*exp(-zt(k)/7000.)
     !Set up concentration of LARGE MODE Mineral Dust (#/kg)
     if(k<=2) md2np(k,i,j)=dust2_maxt
     if(k>2)  md2np(k,i,j)=dust2_maxt*exp(-zt(k)/7000.)

     !Set up Field of SMALL MODE DUST mass (kg/kg)
     md1mp(k,i,j) = ((aero_medrad(3)*aero_rg2rm(3))**3.) &
                    *md1np(k,i,j)/(0.23873/aero_rhosol(3))
     !Set up Field of LARGE MODE DUST mass (kg/kg)
     md2mp(k,i,j) = ((aero_medrad(4)*aero_rg2rm(4))**3.) &
                    *md2np(k,i,j)/(0.23873/aero_rhosol(4))

   !If using dust source model, do not initialize with background dust
   elseif(idust == 2) then
     md1np(k,i,j) = 0.
     md2np(k,i,j) = 0.
     md1mp(k,i,j) = 0.
     md2mp(k,i,j) = 0.
   endif

   !Output sample initial profile
   if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' Dust-init (k,zt,dust1/mg,dust2/mg) on Grid:',ifm
     print'(a10,i5,f11.1,2f17.7)',' DUST-init' &
        ,k,zt(k),md1np(k,i,j)/1.e6,md2np(k,i,j)/1.e6
   endif

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_dust

!##############################################################################
Subroutine init_absorbing_carbon (n1,n2,n3,abc1np,abc2np,abc1mp,abc2mp,dn0,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: abc1np,abc2np,abc1mp,abc2mp,dn0
real :: abc1_maxt,abc2_maxt

! Initialize Absorbing Carbon
if(iaeroprnt==1 .and. print_msg) then
 print*,'Start Initializing Absorbing Carbon concentration'
 print*,'iabcarb,iaerorad',iabcarb,iaerorad
endif

!Convert RAMSIN #/mg to #/kg
 abc1_maxt = abc1_max * 1.e6
 abc2_maxt = abc2_max * 1.e6

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !If not using absorbing carbon / smoke source model
   if(iabcarb == 1) then
     !Set up concentration of mode1 absorbing carbon (#/kg)
     if(k<=2) abc1np(k,i,j)=abc1_maxt
     if(k>2)  abc1np(k,i,j)=abc1_maxt*exp(-zt(k)/7000.)
     !Set up concentration of mode2 absorbing carbon (#/kg)
     if(k<=2) abc2np(k,i,j)=abc2_maxt
     if(k>2)  abc2np(k,i,j)=abc2_maxt*exp(-zt(k)/7000.)

     !Set up Field of mode1 absorbing carbon mass (kg/kg)
     abc1mp(k,i,j) = ((aero_medrad(8)*aero_rg2rm(8))**3.) &
                    *abc1np(k,i,j)/(0.23873/aero_rhosol(8))
     !Set up Field of mode2 absorbing carbon mass (kg/kg)
     abc2mp(k,i,j) = ((aero_medrad(9)*aero_rg2rm(9))**3.) &
                    *abc2np(k,i,j)/(0.23873/aero_rhosol(9))

   !Leaving this section here in case we add an absorbing carbon / smoke 
   !source in the future
   elseif(iabcarb == 2) then
     abc1np(k,i,j) = 0.
     abc2np(k,i,j) = 0.
     abc1mp(k,i,j) = 0.
     abc2mp(k,i,j) = 0.
   endif

   !Output sample initial profile
   if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' AbsorbingCarbon-init (k,zt,abc1/mg,abc2/mg) on Grid:',ifm
     print'(a21,i5,f11.1,2f17.7)',' AbsorbingCarbon-init' &
        ,k,zt(k),abc1np(k,i,j)/1.e6,abc2np(k,i,j)/1.e6
   endif

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_absorbing_carbon

!##############################################################################
Subroutine init_salt (n1,n2,n3,salt_film_np,salt_jet_np,salt_spum_np &
                       ,salt_film_mp,salt_jet_mp,salt_spum_mp,ifm)

use micphys
use rconstants
use mem_grid

implicit none

integer :: n1,n2,n3,i,j,k,ifm
real, dimension(n1,n2,n3) :: salt_film_np,salt_jet_np,salt_spum_np
real, dimension(n1,n2,n3) :: salt_film_mp,salt_jet_mp,salt_spum_mp
real :: saltf_maxt,saltj_maxt,salts_maxt

! Initialize Sea-salt
if(iaeroprnt==1 .and. print_msg) then
 print*,'Start Initializing SALT concentration'
 print*,'isalt,iaerorad',isalt,iaerorad
endif

!Convert RAMSIN #/mg to #/kg
 saltf_maxt = saltf_max * 1.e6
 saltj_maxt = saltj_max * 1.e6
 salts_maxt = salts_max * 1.e6

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !If not using dust source model
   if(isalt == 1) then
     !Set up concentration of FILM MODE SALT (#/kg)
     if(k<=2) salt_film_np(k,i,j)=saltf_maxt
     if(k>2)  salt_film_np(k,i,j)=saltf_maxt*exp(-zt(k)/7000.)
     !Set up concentration of JET MODE Mineral Dust (#/kg)
     if(k<=2) salt_jet_np(k,i,j)=saltj_maxt
     if(k>2)  salt_jet_np(k,i,j)=saltj_maxt*exp(-zt(k)/7000.)
     !Set up concentration of SPUME MODE Mineral Dust (#/kg)
     if(k<=2) salt_spum_np(k,i,j)=salts_maxt
     if(k>2)  salt_spum_np(k,i,j)=salts_maxt*exp(-zt(k)/7000.)

     !Set up 3D Field of FILM MODE SALT mass (kg/kg)
     salt_film_mp(k,i,j) = ((aero_medrad(5)*aero_rg2rm(5))**3.) &
                           *salt_film_np(k,i,j)/(0.23873/aero_rhosol(5))
     !Set up 3D Field of JET MODE SALT mass (kg/kg)
     salt_jet_mp(k,i,j)  = ((aero_medrad(6)*aero_rg2rm(6))**3.) &
                           *salt_jet_np(k,i,j) /(0.23873/aero_rhosol(6))
     !Set up 3D Field of SPUME MODE SALT mass (kg/kg)
     salt_spum_mp(k,i,j) = ((aero_medrad(7)*aero_rg2rm(7))**3.) &
                           *salt_spum_np(k,i,j)/(0.23873/aero_rhosol(7))

   !If using salt source model, do not initialize with background salt
   elseif(isalt == 2) then
     salt_film_np(k,i,j) = 0.
     salt_jet_np(k,i,j)  = 0.
     salt_spum_np(k,i,j) = 0.
     salt_film_mp(k,i,j) = 0.
     salt_jet_mp(k,i,j)  = 0.
     salt_spum_mp(k,i,j) = 0.
   endif

   !Output initial sample profile
   if(iaeroprnt==1 .and. i==1 .and. j==1 .and. print_msg) then
     if(k==1) print*,' Salt-init (k,zt,film/mg,jet/mg,spume/mg) on Grid:',ifm
     print'(a10,i5,f11.1,3f17.7)',' SALT-init',k,zt(k),salt_film_np(k,i,j)/1.e6 &
      ,salt_jet_np(k,i,j)/1.e6,salt_spum_np(k,i,j)/1.e6
   endif

  enddo
 enddo
enddo

if(iaeroprnt==1 .and. print_msg) print*,' '

return
END SUBROUTINE init_salt

!##############################################################################
Subroutine init_tracer (n1,n2,n3,tracerp,dn0,ifm,nsc)

use micphys
use rconstants
use mem_grid
use mem_tracer
use node_mod

implicit none

integer :: n1,n2,n3,i,j,k,ifm,nsc,ii,jj
real, dimension(n1,n2,n3) :: tracerp,dn0
real :: ccn1_maxt

! Initialize Tracers
if(print_msg) print*,'Start Initializing Tracers, Grid:',ifm,' Tracer:',nsc

!Convert RAMSIN #/mg to #/kg
 ccn1_maxt = ccn1_max * 1.e6 

do j = 1,n3
 do i = 1,n2
  do k = 1,n1

   !Get absolute grid points for parallel (& sequential) computation
   ii = i+mi0(ngrid)
   jj = j+mj0(ngrid)

   !Set up Vertical profile, Exponential decrease that scales with pressure
   if(nsc==1) then
    if(k<=2) tracerp(k,i,j)=ccn1_maxt
    if(k>2)  tracerp(k,i,j)=ccn1_maxt*exp(-zt(k)/7000.)
   endif
   !Set up Field of CCN mass mixing ratio (kg/kg)
   if(nsc==2) then
    tracerp(k,i,j) = ((aero_medrad(1)*aero_rg2rm(1))**3.) &
                *tracer_g(1,ifm)%tracerp(k,i,j)/(0.23873/aero_rhosol(1))
   endif

  enddo
 enddo
enddo

if(nsc.eq.itracer .and. print_msg) print*,' '

return
END SUBROUTINE init_tracer

!##############################################################################
Subroutine jnmbinit ()

use micphys
use mem_grid, only:print_msg

implicit none

if (level /= 3) then

   if (level <= 1) then
      jnmb(1) = 0
   else
      jnmb(1) = 4
   endif

   jnmb(2) = 0
   jnmb(3) = 0
   jnmb(4) = 0
   jnmb(5) = 0
   jnmb(6) = 0
   jnmb(7) = 0
   jnmb(8) = 0

else

   jnmb(1) = icloud
   jnmb(2) = irain
   jnmb(3) = ipris
   jnmb(4) = isnow
   jnmb(5) = iaggr
   jnmb(6) = igraup
   jnmb(7) = ihail
   jnmb(8) = idriz

   if (icloud .eq. 1) jnmb(1) = 4
   if (irain  .eq. 1) jnmb(2) = 2
   if (ipris  .ge. 1) jnmb(3) = 5
   if (isnow  .eq. 1) jnmb(4) = 2
   if (iaggr  .eq. 1) jnmb(5) = 2
   if (igraup .eq. 1) jnmb(6) = 2
   if (ihail  .eq. 1) jnmb(7) = 2
   if (idriz  .eq. 1) jnmb(8) = 4

   if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or.  &
      igraup == 5 .or. ihail == 5) then

      if (irain  >= 1) jnmb(2) = 5
      if (isnow  >= 1) jnmb(4) = 5
      if (iaggr  >= 1) jnmb(5) = 5
      if (igraup >= 1) jnmb(6) = 5
      if (ihail  >= 1) jnmb(7) = 5

   endif

endif

if(print_msg) then
 print*,''
 print*,'HYDROMETEOR SETTINGS'
 print*,'icloud ',icloud,'  ','jnmb-cloud   ',jnmb(1)
 print*,'idriz  ',idriz, '  ','jnmb-drizzle ',jnmb(8)
 print*,'irain  ',irain, '  ','jnmb-rain    ',jnmb(2)
 print*,'ipris  ',ipris, '  ','jnmb-pris    ',jnmb(3)
 print*,'isnow  ',isnow, '  ','jnmb-snow    ',jnmb(4)
 print*,'iaggr  ',iaggr, '  ','jnmb-aggr    ',jnmb(5)
 print*,'igraup ',igraup,'  ','jnmb-graup   ',jnmb(6)
 print*,'ihail  ',ihail, '  ','jnmb-hail    ',jnmb(7)
 print*,''
endif

return
END SUBROUTINE jnmbinit

!##############################################################################
Subroutine micinit ()

use micphys
use mem_grid, only:print_msg

implicit none

integer :: lhcat,lcat,ia
integer, dimension(16) :: lcat0
real :: cfmasi,c1,glg,glg1,glg2,glgm,glgc,flngi,dpsi,embsip,dnsip
real, external :: gammln
real, external :: gammp
real, external :: gammq

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/ ! lcat corressponding to lhcat

! Initialize arrays based on microphysics namelist parameters
! Note: pparm for pristine ice is obsolete since IPRIS = 0 or 5 only

parm(1) = cparm
parm(2) = rparm
parm(4) = sparm
parm(5) = aparm
parm(6) = gparm
parm(7) = hparm
parm(8) = dparm

if (icloud .eq. 1) parm(1) = .3e9
if (irain  .eq. 1) parm(2) = .1e-2
if (isnow  .eq. 1) parm(4) = .1e-2
if (iaggr  .eq. 1) parm(5) = .1e-2
if (igraup .eq. 1) parm(6) = .1e-2
if (ihail  .eq. 1) parm(7) = .3e-2
if (idriz  .eq. 1) parm(8) = .1e6  !# per kg ~ m^3 
                                   !(mid-range avg from Feingold(99)

if (icloud .eq. 0 .or. icloud .eq. 5) parm(1) = 0.
if (irain  .eq. 0 .or. irain  .eq. 5) parm(2) = 0.
if (isnow  .eq. 0 .or. isnow  .eq. 5) parm(4) = 0.
if (iaggr  .eq. 0 .or. iaggr  .eq. 5) parm(5) = 0.
if (igraup .eq. 0 .or. iaggr  .eq. 5) parm(6) = 0.
if (ihail  .eq. 0 .or. ihail  .eq. 5) parm(7) = 0.
if (idriz  .eq. 0 .or. idriz  .eq. 5) parm(8) = 0. 

if(print_msg) then
 print*,''
 print*,'More HYDROMETEOR SETTINGS'
 print*,'Not used for: (1)No micro, (2)two-moment micro, (3)Bin-micro'
 if(icloud==0 .or. icloud==5)then
  print*,'cparm - cloud   ','NOT USED'
 else
  print*,'cparm - cloud   ',parm(1)
 endif

 if(idriz==0 .or. idriz==5)then
  print*,'dparm - drizzle ','NOT USED'
 else
  print*,'dparm - drizzle ',parm(8)
 endif

 if(irain==0 .or. irain==5)then
  print*,'rparm - rain    ','NOT USED'
 else
  print*,'rparm - rain    ',parm(2)
 endif

 if(ipris==0 .or. ipris==5)then
  print*,'pparm - pris    ','NOT USED'
 else
  print*,'pparm - pris    ','NOT USED'
 endif

 if(isnow==0 .or. isnow==5)then
  print*,'sparm - snow    ','NOT USED'
 else
  print*,'sparm - snow    ',parm(4)
 endif

 if(iaggr==0 .or. iaggr==5)then
  print*,'aparm - aggr    ','NOT USED'
 else
  print*,'aparm - aggr    ',parm(5)
 endif

 if(igraup==0 .or. igraup==5)then
  print*,'gparm - graup   ','NOT USED'
 else
  print*,'gparm - graup   ',parm(6)
 endif

 if(ihail==0 .or. ihail==5)then
  print*,'hparm - hail    ','NOT USED'
 else
  print*,'hparm - hail    ',parm(7)
 endif

 print*,''
endif

! 125 micron max diameter ice crystal size for participating
! in Hallet-Mossop 2ndary ice splintering process.
dps = 125.e-6
dps2 = dps ** 2

rictmin = 1.0001
rictmax = 0.9999 * float(nembc)

do lhcat = 1,nhcat
   lcat=lcat0(lhcat)

   cfden(lhcat) = cfmas(lhcat) * 6.0 / 3.14159
   pwden(lhcat) = pwmas(lhcat) - 3.
   emb0log(lcat) = log(emb0(lcat))
   emb1log(lcat) = log(emb1(lcat))

! Define coefficients [frefac1, frefac2] used for terminal velocity
! and Reynolds number

   cfmasi = 1. / cfmas(lhcat)
   pwmasi(lhcat) = 1. / pwmas(lhcat)
   pwen0(lhcat) = 1. / (pwmas(lhcat) + 1.)
   pwemb0(lhcat) = pwmas(lhcat) / (pwmas(lhcat) + 1.)
   c1 = 1.5 + .5 * pwvt(lhcat)

   glg   = gammln(gnu(lcat))
   glg1  = gammln(gnu(lcat) + 1.)
   glg2  = gammln(gnu(lcat) + 2.)
   glgm  = gammln(gnu(lcat) + pwmas(lhcat))
   glgc  = gammln(gnu(lcat) + c1)

   if (jnmb(lcat) .eq. 3) then
      cfemb0(lhcat) = cfmas(lhcat) * exp(glgm - glg)  &
         ** pwen0(lhcat) * (1. / parm(lcat)) ** pwemb0(lhcat)
      cfen0(lhcat) = parm(lcat) * (exp(glg - glgm) / parm(lcat))  &
         ** pwen0(lhcat)
   endif

   dnfac(lhcat) = (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac1(lhcat) = shape(lhcat) * exp(glg1 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** pwmasi(lhcat)

   frefac2(lhcat) = shape(lhcat) * 0.229 * sqrt(cfvt(lcat))  &
      * (cfmasi * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)  &
      * exp(glgc - glg)

   sipfac(lhcat) = .785 * exp(glg2 - glg)  &
      * (cfmasi * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

   dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

   dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))
   if (lhcat .le. 4) gamm(lhcat) = exp(glg)
   if (lhcat .le. 4) gamn1(lhcat) = exp(glg1)

! gam1   :  the integral of the pristine distribution from dps to infty
! gam2   :  the integral of the snow dist. from 0 to dps
! gam3   :  values of the exponential exp(-dps/dn)

enddo

!***********************************************************************
!************** Secondary Ice Production Arrays ************************
!***********************************************************************
flngi = 1. / float(ngam)
do ia=1,ngam
   dpsi = dps * 1.e6 / float(ia)

   gam(ia,1) = gammq(gnu(3) + 1., dpsi)
   gam(ia,2) = gammp(gnu(4) + 1., dpsi)
   gam(ia,3) = exp(-dpsi)

   GAMINC(IA,1)=gammq(GNU(3),dpsi)
   GAMINC(IA,2)=gammp(GNU(4),dpsi)

   embsip = emb1(1) * float(ia) * flngi
   dnsip = dnfac(1) * embsip ** pwmasi(1)
   gamsip13(1,ia)=gammp(gnu(1),13.e-6/dnsip)
   gamsip24(1,ia)=gammq(gnu(1),24.e-6/dnsip)

   embsip = emb1(8) * float(ia) * flngi
   dnsip = dnfac(16) * embsip ** pwmasi(16)
   gamsip13(2,ia)=gammp(gnu(8),13.e-6/dnsip)
   gamsip24(2,ia)=gammq(gnu(8),24.e-6/dnsip)
enddo

return
END SUBROUTINE micinit

!##############################################################################
Subroutine setupF94 ()

!Routine to set up rain size bin partitions for F94 collection
!routine for rain-pris ice collection forming hail

use micphys

implicit none

integer :: nbins,i
real :: pi314,mind,maxd,m,ii3,rmin,mmin,xjo,massdobl,tmp1

!--local constants
nbins=94
pi314=4.*ATAN(1.)
ii3=(1./3.)
massdobl=2. !mass doubling bin factor for rain

!--set up rain diameter bin partitions
mind=3.4668048E-07            !min rain diam [m] from subr 'mkcoltb'
maxd=1.7334025E-02            !max rain diam [m] from subr 'mkcoltb'
rmin=mind/2.                  !minimum radius
mmin=(4./3.)*pi314*(rmin)**3  !minimum mass
tmp1=massdobl*3.              !mass is doubled every 'massdobl' # of bins
xjo=tmp1/LOG(2.)
do i=1,nbins+1 !*NOTE: nbins+1 needed to numerically integrate Fgam
 m=mmin*exp(3*float(i-1)/(xjo))
 rdbF94(i)=2.*((3.*m/4./pi314)**ii3)  !bin diameter bounds
enddo

!--set up avg bin diams, term vel's, and mass
do i=1,nbins
 radF94(i)=0.5*(rdbF94(i)+rdbF94(i+1))    !avg bin diam
 avtF94(i)=cfvt(2)*(radF94(i))**pwvt(2)   !term vel of avg diam
 ramF94(i)=cfmas(2)*(radF94(i))**pwmas(2) !mass of drop w/ avg diam
enddo

return
END SUBROUTINE setupF94
