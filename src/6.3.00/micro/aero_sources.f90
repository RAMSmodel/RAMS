!##############################################################################
Subroutine dust_sources (m1,m2,m3,ia,iz,ja,jz,rtgt,glat,glon         &
                        ,u,v,dn0,md1np,md2np,md1mp,md2mp             &
                        ,dustfrac,swater,parea,lclass,stext          &
                        ,vegrough                                    &
                        )

use micphys, only:idustloft,aero_rhosol
use mem_grid, only:dtlt,ngrid,nnxp,nnyp,nzg,npatch
use mem_varinit, only:nud_type,nudlat
use node_mod, only:mi0,mj0

implicit none

integer, parameter :: ncls=7
integer :: m1,m2,m3,ia,iz,ja,jz,i,j,z0limit,veglimit,dustbc,allowdustloft
real :: dust1_size,dust2_size,aux1,aux2,aux3
real, dimension(m2,m3) :: u10,v10,rtgt,num_dust_s,num_dust_l &
                         ,glat,glon,dustfrac
real, dimension(m1,m2,m3) :: u,v,md1np,md2np,md1mp,md2mp,dn0
real, dimension(nzg,m2,m3,npatch) :: swater
real, dimension(m2,m3,npatch) :: parea,lclass
real, dimension(nzg,m2,m3,npatch) :: stext
real, dimension(m2,m3,npatch) :: vegrough

 !******************************************************************
 !Dust lofting option default flags
 !******************************************************************
 !Limit dust flux due to vegetative roughness
 z0limit=1          !0 = off 1 = on (only for idustloft=0)
 !Limit dust lofting to specific bare or low grass surface types
 !3-Desert,bare soil, 8-Short-grass, 9-Tall-grass, 10-Semi-desert,
 !15-Crops, 18-Wooded-grassland
 veglimit=1         !0 = off 1 = on (only for idustloft=0)
 !******************************************************************

 !For Ginoux(2001) or NRL database dust lofting,
 !ignore surface and veg limits.
 if(idustloft==1 .or. idustloft==2 .or. idustloft==22) then 
   z0limit=0
   veglimit=0
 endif

 !Get the 10m wind components
 CALL get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,u,v,u10,v10,rtgt)

 !Compute dust sources one gridcell at a time
 CALL dustflux (u10,v10,swater,npatch,m1,m2,m3,ia,iz,ja,jz     &
       ,rtgt,dtlt,num_dust_s,num_dust_l,ncls,z0limit           &
       ,veglimit,glat,glon,dustfrac,parea,nzg,lclass,stext     &
       ,vegrough                                               &
       )

 !Set mean mass weighted radii of dust modes that corressponds to dust func.
 !Bin radii(um): r1=0.15,r2=0.265,r3=0.471,r4=0.838,r5=1.5,r6=2.65,r7=4.71
 !Apportion flux into size bins with sp, where sp are bin mass fractions.
 !For the 4 small clay size bins, each bin contains the following fractions:
 ! 1=0.9%, 2=8.1%, 3=23.4%, 4=67.6%
 !However, the 4 small clay bins represent only 1/10th the total dust silt,
 ! so the parameter "sp(1:4)" represents decimal-percentage/10
 !The three large classes are weight equally at 33.333% each.
 !sp(1:4)     = (/0.0009,0.0081,0.0234,0.0676/)
 !sp(5:ncls)  = 0.30
 !dust1_size = mass weighted mean fine dust radius (r1 - r4)
 !dust2_size = mass weighted mean course dust radius (r5 - r7)
 !fine: 0.15*0.009 + 0.265*0.081 + 0.471*0.234 + 0.838*0.676 = 0.699 um
 !coarse: 1.5*0.333 + 2.65*0.333 + 4.71*0.333 = 2.95038 microns
 dust1_size=0.699e-6 !small mode (meters)
 dust2_size=2.95e-6  !large mode (meters)

 !Update dust particles due to dust sources. Dust source model outputs in 
 !#/cm3 so we convert to #/kg and compute dust mass mixing ratio in kg/kg

 !First prevent dust lofting along boundaries due to potential constant
 !inflow issues causing dust to be continually lofted. Value here is the number
 !of grid cells in from lateral boundaries where lofting is allowed.
 dustbc=5
 
 do j = ja,jz
  do i = ia,iz

 !SaleebySCM:This condition modified in SCM version since it requires
 !knowledge of horizontal boundary conditions
   allowdustloft=0
   if(i+mi0(ngrid) > dustbc .and. i+mi0(ngrid) < nnxp(ngrid)-dustbc .and. &
      j+mj0(ngrid) > dustbc .and. j+mj0(ngrid) < nnyp(ngrid)-dustbc) then
    allowdustloft=1
   endif

   if(allowdustloft==1) then
    !Convert dust number (#/cm3) to (#/kg)
    num_dust_s(i,j) = num_dust_s(i,j) / dn0(2,i,j) * 1.e6
    num_dust_l(i,j) = num_dust_l(i,j) / dn0(2,i,j) * 1.e6

    !Update number of dust species
    md1np(2,i,j) = md1np(2,i,j) + num_dust_s(i,j)
    md1np(1,i,j) = md1np(2,i,j)

    md2np(2,i,j) = md2np(2,i,j) + num_dust_l(i,j)
    md2np(1,i,j) = md2np(2,i,j)

    !Update mass of dust species using mean mass radii rather than median radii
    md1mp(2,i,j) = md1mp(2,i,j) + &
        (dust1_size**3.)*num_dust_s(i,j)/(0.23873/aero_rhosol(3))
    md1mp(1,i,j) = md1mp(2,i,j)

    md2mp(2,i,j) = md2mp(2,i,j) + &
        (dust2_size**3.)*num_dust_l(i,j)/(0.23873/aero_rhosol(4))
    md2mp(1,i,j) = md2mp(2,i,j)
   endif

  enddo
 enddo

return
END SUBROUTINE dust_sources

!##############################################################################
Subroutine salt_sources (m1,m2,m3,ia,iz,ja,jz,rtgt,u,v,dn0       &
                       ,salt_film_np,salt_jet_np,salt_spum_np    &
                       ,salt_film_mp,salt_jet_mp,salt_spum_mp    &
                       ,parea,lclass)

use micphys, only:aero_rg2rm,aero_rhosol
use mem_grid, only:dtlt,npatch

implicit none

!Assumed salt masses based on radii below and density of salt 2.165 g/cm3
!Assumes all aerosols in each category are all the same size
!Need to alter these if we assume a certain distribution

!Salt mass and number arrays
real, dimension(m1,m2,m3) :: salt_film_np,salt_jet_np,salt_spum_np
real, dimension(m1,m2,m3) :: salt_film_mp,salt_jet_mp,salt_spum_mp

integer :: m1,m2,m3,ia,iz,ja,jz,i,j
real, dimension(m2,m3)    :: u10,v10,rtgt
real, dimension(m1,m2,m3) :: u,v,dn0
real, dimension(m2,m3)    :: film_source,jet_source,spm_source
real, dimension(m2,m3,npatch) :: parea
real, dimension(m2,m3,npatch) :: lclass
real :: sf_size,sj_size,ss_size

!Get the 10m wind components
 CALL get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,u,v,u10,v10,rtgt)

!Get Salt emission
 CALL salt_flux (u10,v10,m1,m2,m3,ia,iz,ja,jz,npatch,dtlt &
               ,salt_film_np,salt_jet_np,salt_spum_np     &
               ,film_source,jet_source,spm_source,dn0     &
               ,parea,lclass)

!Set median radii of salt modes that corresponds to salt func
sf_size=0.10e-6 !film mode radius (m3)
sj_size=1.00e-6 !jet mode radius (m3)
ss_size=6.00e-6 !spume mode radius (m3)

!Nudge diagnostic values over sea
do j = ja,jz
 do i = ia,iz

    !Convert dust number (#/m3) to (#/kg)
    film_source(i,j) = film_source(i,j) / dn0(2,i,j)
    jet_source(i,j)  = jet_source(i,j)  / dn0(2,i,j)
    spm_source(i,j)  = spm_source(i,j)  / dn0(2,i,j)

    !Update number of salt species
    salt_film_np(2,i,j) = salt_film_np(2,i,j) + film_source(i,j)
    salt_film_np(1,i,j) = salt_film_np(2,i,j)

    salt_jet_np(2,i,j)  = salt_jet_np(2,i,j)  + jet_source(i,j)
    salt_jet_np(1,i,j)  = salt_jet_np(2,i,j)

    salt_spum_np(2,i,j) = salt_spum_np(2,i,j) + spm_source(i,j)
    salt_spum_np(1,i,j) = salt_spum_np(2,i,j)
    
    !Update mass of salt species
    salt_film_mp(2,i,j) = salt_film_mp(2,i,j) + &
        ((sf_size*aero_rg2rm(5))**3.)*film_source(i,j)/(0.23873/aero_rhosol(5))
    salt_film_mp(1,i,j) = salt_film_mp(2,i,j)

    salt_jet_mp(2,i,j) = salt_jet_mp(2,i,j) + &
        ((sj_size*aero_rg2rm(6))**3.)*jet_source(i,j) /(0.23873/aero_rhosol(6))
    salt_jet_mp(1,i,j) = salt_jet_mp(2,i,j)

    salt_spum_mp(2,i,j) = salt_spum_mp(2,i,j) + &
        ((ss_size*aero_rg2rm(7))**3.)*spm_source(i,j) /(0.23873/aero_rhosol(7))
    salt_spum_mp(1,i,j) = salt_spum_mp(2,i,j)

 enddo
enddo

return
END SUBROUTINE salt_sources

!##############################################################################
Subroutine get_u10_v10 (m1,m2,m3,ia,iz,ja,jz,uu,vv,u10,v10,rtgt)

use mem_grid, only:jdim,zm
use micphys, only:iscm

implicit none

integer :: m1,m2,m3,i,j,ia,iz,ja,jz
real, parameter :: z0 = 0.05
real :: zref,aux1,aux2
real, dimension(m2,m3) :: rtgt
real, dimension(m2,m3) :: u10,v10
real, dimension(m1,m2,m3) :: uu,vv

aux2 = log(10./z0)

! Assumed roughness length only really valid for bare soil surfaces.
! Fine for now as the source function is only valid for bare surfacs,
! but if the source func is changed in the future, the roughness
! length from LEAF2 should be used instead of an assumed constant.

! Log law wind profile for the lowest 100m AGL. To get wind (V) at an 
! unknown height (Z) from a known wind (Vref) at reference height (Zref):
! V ~ Vref [ ln(Z/Z0) / ln(Zref/Z0) ]

do j=ja,jz
 do i=ia,iz

  zref = zm(2) * rtgt(i,j)
  aux1 = 1./log(zref/z0)
  u10(i,j) = aux1 * aux2 * ((uu(2,i-1,j) + uu(2,i,j)) * 0.5)
  v10(i,j) = aux1 * aux2 * ((vv(2,i,j-jdim) + vv(2,i,j)) * 0.5)

if(ISCM>=1)then
 !SaleebySCM:Override u10,v10 with approx wind at T point using single i,j
  u10(i,j) = aux1 * aux2 * uu(2,i,j)
  v10(i,j) = aux1 * aux2 * vv(2,i,j)
endif

 enddo
enddo

return
END SUBROUTINE get_u10_v10

!##############################################################################
Subroutine salt_flux (u10,v10,m1,m2,m3,ia,iz,ja,jz,npatch,dtlt  &
                     ,salt_film_np,salt_jet_np,salt_spum_np     &
                     ,film_source,jet_source,spm_source,dn0     &
                     ,xparea,xlclass)

!Code to calculate seasalt number flux on global scale based on the O'dowd
!1997, 1999 where u10 is the 10-meter wind speed m/s, salt_film is the number
!of particles /m3 in the submicron mode (ccn), and salt_jet is the number of
!particles /m3 in the supermicron mode (gccn), and salt_spume is the number
!particles /m3 in the ultra-giant mode (ultra-gccn), and film and jet medians
!are user defined in centimeters.
!The total number concentrations for each mode were observed to be:
! log Nfilm = 0.095 U10 + 6.2830, 0.1 microns mode radius ----> CCN
! log Njet = 0.0422 U10 + 5.7122, 1 micron mode radius -------> GCCN
! log Nspume= 0.069 U10 + 0.1900, 6 micron mode radius -------> Super giant
! O’Dowd C. D., Smith M. H. and Jennings S. G. (1993) Submicron aerosol, radon and
! soot carbon characteristics over the North East Atlantic. J. geophys. Res. 98,
! 1132-1136. Updated for O'Dowd 1997,1999.

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,npatch,itype
real :: dtlt,time_salt,patarea
real, dimension(m1,m2,m3) :: salt_film_np,salt_jet_np,salt_spum_np,dn0
real, dimension(m2,m3) :: salt_film_dia,salt_jet_dia,salt_spm_dia
real, dimension(m2,m3) :: wind,u10,v10,film_source,jet_source,spm_source
double precision, dimension(m2,m3)::salt_film_dia2,salt_jet_dia2 &
  ,salt_spm_dia2,wind_2
real, dimension(m2,m3,npatch) :: xparea
real, dimension(m2,m3,npatch) :: xlclass

time_salt = dtlt !relaxation time in seconds

do j = ja,jz
do i = ia,iz

 !Zero out the final source func prior to updated timestep computation
 film_source(i,j) = 0.0
 jet_source(i,j)  = 0.0
 spm_source(i,j)  = 0.0

 !Compute salt source term
 !Limit wind speed to 30m/s (Fan and Toon 2011)

 patarea=xparea(i,j,1)
 itype = nint(xlclass(i,j,1))

 if(itype .eq. 0 .and. patarea .gt. 0.009) then
  !Compute wind speed from 10m u,v winds in [m/s]
  wind(i,j) = min(30.,sqrt(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j)))
  wind_2(i,j) = dble(wind(i,j))
  !Compute predicted salt concentration based on wind speed (/m3)
  salt_film_dia2(i,j) = 10.**(0.0950*wind_2(i,j) + 6.2830)
  salt_film_dia(i,j)  = sngl(salt_film_dia2(i,j))
  salt_jet_dia2(i,j)  = 10.**(0.0422*wind_2(i,j) + 5.7122)
  salt_jet_dia(i,j)   = sngl(salt_jet_dia2(i,j))
  salt_spm_dia2(i,j)  = 10.**(0.0690*wind_2(i,j) + 0.1900)
  salt_spm_dia(i,j)   = sngl(salt_spm_dia2(i,j))
  !Note that sea salt predicted (/m3) while we carry (/kg) in the model.
  !So do correct density conversion here for adding particles.
  if(salt_film_dia(i,j) > salt_film_np(2,i,j)*dn0(2,i,j)) &
   film_source(i,j)=film_source(i,j) + (salt_film_dia(i,j) &
      - salt_film_np(2,i,j)*dn0(2,i,j)) * dtlt / time_salt * patarea
  if(salt_jet_dia(i,j)  > salt_jet_np(2,i,j)*dn0(2,i,j)) &
   jet_source(i,j) =jet_source(i,j)  + (salt_jet_dia(i,j)  &
      - salt_jet_np(2,i,j)*dn0(2,i,j))  * dtlt / time_salt * patarea
  if(salt_spm_dia(i,j)  > salt_spum_np(2,i,j)*dn0(2,i,j)) &
   spm_source(i,j) =spm_source(i,j)  + (salt_spm_dia(i,j)  &
      - salt_spum_np(2,i,j)*dn0(2,i,j)) * dtlt / time_salt * patarea
 endif

enddo
enddo

return
END SUBROUTINE salt_flux

!##############################################################################
Subroutine dustflux (u10,v10,swater,npatch,m1,m2,m3,ia,iz,ja,jz     &
         ,rtgt,dtlt,num_dust_s,num_dust_l,ncls,z0limit              &
         ,veglimit,glat,glon,dustfrac,xparea,nzg,xlclass,xstext     &
         ,xvegrough                                                 &
         )

!Code to calculate dust flux on regional-global scales based on dust code
!from Peter Colarco at Goddard. It calculates fluxes using a 1x1 degree
!horizontal resolution source func developed by Paul Ginoux. Source
!func is equivalent to the fraction of the grid cell emitting dust.
!Fluxes are calculated based on formulas given in Marticorena & Bermagetti,
!1995 using 10-m winds and soil moisture. Updated for Fecan et al.(1999)
!and Pierre et al.(2012).

use mem_grid, only:dzt
use micphys, only:iscm

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ncls,npatch,icls,i,j,np,itype,sclass
integer :: z0limit,veglimit,vegstop,isource,jsource
real :: patarea,netrough

!Set number of bins, air density (g/cm3), gravity (cm/s2)
!and "fact" as parameters. "fact" scales mass emissions to 1851 tg
real :: rhoa,grav,fact,pi,S
parameter(rhoa=1.25d-3, grav=981., fact=1.153d-16, pi=3.141593)

real :: dz,dtlt,rmrat,rbmin,rho_p,uth75,rmassmin,Ez0,gwet
real, dimension(ncls) :: sp,den,uth,diam,rmass,r,utwet &
                        ,mass_flux,mass_conc,num_conc
real, dimension(m2,m3) :: num_dust_s,num_dust_l,wind,u10,v10,rtgt &
                         ,glat,glon,dustfrac
real, dimension(12) :: clayfrac,slmsts,wprime
real, dimension(nzg,m2,m3,npatch) :: swater
real, dimension(m2,m3,npatch) :: xparea,xlclass
integer :: nzg
real, dimension(nzg,m2,m3,npatch) :: xstext
real, dimension(m2,m3,npatch) :: xvegrough

!Define clay percentage for soil textural classes
!  1 sand            5 loam                 9 sandy clay
!  2 loamy sand      6 sandy clay loam      10 silty clay
!  3 sandy loam      7 silty clay loam      11 clay
!  4 silt loam       8 clay loam            12 peat
data clayfrac/0.0,5.0,10.0,12.0,18.0,28.0,33.0,33.0,42.0,48.0,70.0,0.0/
data slmsts/.395,.410,.435,.485,.451,.420,.477,.476,.426,.492,.482,.863/
!Using Fecan et al.(1999) paramerization based on clay percentage
!w'=0.0014(clayfrac)^2 + 0.17(clayfrac)
data wprime/0.00000,0.885000,1.84000,2.24160,3.51360,5.85760,7.13460 &
           ,7.13460,9.60960,11.3856,18.7600,0.00000/

rmrat = (100**3)**(1./8.) !rmrat=volume ratio between bins
rbmin  = 1.0e-5*((1.+rmrat)/2.)**(1./3.) !radius of smallest bin [cm]
rho_p = 2.65 !particle density [g/cm3]

!Set up the bins/apportionment (Ginoux et al.2001)
!r1=0.15, r2=0.265, r3=0.471, r4=0.838, r5=1.5, r6=2.65, r7=4.71 microns
rmassmin = 4./3.*pi*rho_p*rbmin**3.
do icls = 1,ncls
  rmass(icls) = rmassmin*rmrat**(icls-1)
  r(icls) = (rmass(icls)/rho_p/(4./3.*pi))**(1./3.) !units of cm
enddo

!Apportion flux into size bins with sp, where sp are bin mass fractions
!For the 4 small clay size bins, each bin contains the following fractions:
! 1=0.9%, 2=8.1%, 3=23.4%, 4=67.6%
!However, the 4 small clay bins represent only 1/10th the total dust silt,
! so the parameter "sp(1:4)" represents decimal-percentage/10
!The three large classes are weight equally at 30%.
sp(1:4)     = (/0.0009,0.0081,0.0234,0.0676/)
sp(5:ncls)  = 0.30

diam        = 2*r
den(1:4)    = 2.50 !mass density of clay
den(5:ncls) = 2.65 !mass density of silt

!Threshold velocity calcs [cm/s] from
!Marticorena & Bermagetti(1995) Equations 3-6
do icls = 5,ncls
  uth(icls) = 0.13*sqrt(den(icls)*grav*diam(icls)/rhoa) &
              *sqrt(1. + 0.006/(den(icls)*grav*diam(icls)**2.5)) &
              /sqrt(1.928*(1331.*diam(icls)**1.56 + 0.38)**0.092 - 1.)
enddo

!For sub-micron particles assume uth is that of an reff=0.75 micron particle
!Gives higher uth for smaller particle sizes
uth75 = 0.13*sqrt(2.5*grav*(1.5e-4)/rhoa) &
        *sqrt(1. + 0.006/(2.5*grav*(1.5e-4)**2.5)) &
        /sqrt(1.928*(1331.*(1.5e-4)**1.56 + 0.38)**0.092 - 1.)
uth(1:4) = uth75

!Calculate the flux, mass and # concentration in the domain
!Compute fluxes and output as number distributions with median radii (cm)
do j = ja,jz
do i = ia,iz

 !Zero out the final source func prior to updated timestep computation
 num_dust_s(i,j) = 0.0
 num_dust_l(i,j) = 0.0

 !Set dust source erodible fraction variable S
 S = dustfrac(i,j)

 do np = 2,npatch

  patarea=xparea(i,j,np)

  if(patarea .gt. 0.009) then
   itype = nint(xlclass(i,j,np))
   sclass = floor(xstext(nzg,i,j,np))
   netrough = xvegrough(i,j,np)

if(ISCM>=1)then
 !SaleebySCM:Override land surface class info to be similar for SCM
 netrough = 3.00e-5
endif

   !If using Pierre et al.(2012) surface zoughness limit on dust lofting
   !if patch roughness > 3.10e-3 cm
   if(z0limit==1.and.netrough*100. > 3.10e-3) then
      Ez0 = 0.7304 - (0.0804 * log10(netrough*100.))
   else
      Ez0 = 1.0
   endif

   !If limiting dust lofting to only low or limited vegetation patches
   if(veglimit==1 .and. itype.ne.3  .and. itype.ne.8  .and. itype.ne.9 .and. & 
       itype.ne.10 .and. itype.ne.15 .and. itype.ne.18) then
     vegstop=1
   else
     vegstop=0
   endif

   !No dust lofting at all for water, ice, marsh, or urban
   if(itype==0.or.itype==1.or.itype==2.or.itype==17.or. &
      itype==19.or.itype==21) vegstop=1

   !Proceed to dust lofting if vegetation type allowed
   if(vegstop==0) then
    !Seigel (10-8-2010): Threshold velocity adjustment, accounting for soil 
    !moisture AND soil textural class. Note that swater is volumetric soil 
    !moisture m3/m3 and not total saturation fraction. Ginoux et al.(2001) uses 
    !total saturation fraction with values from 0.001 to 1.0. But Fecan et al.(1999)
    !uses volumetric soil moisture (m3/m3). So be careful switching between 
    !the two. For RAMS, the maximum volumetric soil moisture is slmsts(sclass).

    !For using Fecan et al.(1999) threshold velocity by soil type, then
    !use volumetric soil moisture (m3/m3) for velocity adjustment. This
    !adjustment is more constraing than Ginoux et al.(2001) and may help
    !prevent over-lofting in some circumstances.
    gwet = swater(nzg,i,j,np)
    do icls = 1,ncls
       if(gwet*100. .lt. wprime(sclass)) then
          utwet(icls) = uth(icls)
       else
          utwet(icls) = uth(icls)*sqrt(1 + 1.21* &
                           (gwet*100.-wprime(sclass))**0.68)
       endif
    enddo

    !Compute wind speed from 10m u,v winds in [cm/s]
    wind(i,j) = sqrt(u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j))*100.

    !Grid cell depth for lowest grid cell above ground
    dz = rtgt(i,j) / dzt(2)

    !Flux dust if soil is dry (soil moisture fraction (0->1) < 50%)
    !This lofts dust within the Ginoux et al.(2001) maximum soil moisture 
    !fraction limit for lofting just.
    gwet = max(0.0,min(1.0,swater(nzg,i,j,np)/slmsts(sclass)))
    if(gwet.lt.0.5) then
      do icls = 1,ncls
        !Calculate the mass flux in the domain [g cm-2 s-1]
        !Flux equation from Ginoux et al. 2001 equation 2.
        mass_flux(icls)=Ez0*S*fact*sp(icls)*(wind(i,j)**2.) &
                        *(wind(i,j)-utwet(icls))
        if(mass_flux(icls) <= 0.) then
          mass_flux(icls) = 0.
          mass_conc(icls) = 0.
          num_conc(icls)  = 0.
        endif
        if(mass_flux(icls) > 0.) then
          !Convert mass flux to mass concentration [g cm-3]
          !For mass concentration divide by grid cell depth (dz) in units of [cm]
          mass_conc(icls) = mass_flux(icls) * dtlt / (dz*100.)
          !Then convert to # concentration (particles/cm3) in each bin
          num_conc(icls)  = mass_conc(icls) / &
                               (den(icls)*(4./3.)*pi*r(icls)**3)
        endif
      enddo
      !Sum the sub and super micron bins for num_dust_s & num_dust_l
      do icls = 1,4
         num_dust_s(i,j) = num_dust_s(i,j) + num_conc(icls) * patarea
      enddo
      do icls = 5,ncls
         num_dust_l(i,j) = num_dust_l(i,j) + num_conc(icls) * patarea
      enddo
    endif !IF SOURCE EXISTS AND SOIL IS DRY 

   endif !IF patch is capable of releasing dust
  endif !IF patch area is big enough
 enddo !LOOP over surface patches

enddo !LOOP OVER J POINTS
enddo !LOOP OVER I POINTS

return
END SUBROUTINE dustflux
