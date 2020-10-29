!##############################################################################
Subroutine canopy (mzg,mzs,ksn,nveg  &
   ,soil_water,soil_text,sfcwater_mass  &
   ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height &
   ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class  &
   ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort)

use leaf_coms
use rconstants

implicit none

integer :: mzg,mzs,ksn,nveg

real, dimension(mzg) :: soil_water,soil_text
real, dimension(mzs) :: sfcwater_mass
real :: ustar,tstar,rstar,soil_rough,veg_rough,veg_height             &
       ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class    &
       ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort

integer :: k,nsoil,iter_can

real :: aux,brad,bswp,bthi,btlo,bvpd,c2,c3,c4,es,fac,factv         &
   ,fracliqv,fthi,ftlo,frad,fswp,fvpd,qwtot,rasgnd,rasveg,rleaf,rsatveg     &
   ,sigmaw ,slai,stai,slpotv,srad,sswp,sthi,stlo,svpd,swp,tveg,tvegc,tvegk  &
   ,vpd,wtemp,zognd,zoveg,zdisp,zveg,wflx,dewgndflx,ustar0         &
   ,transp_test,rc,rc_inf,wshed0,transp0
real, external :: eslf
!     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar)) 
!     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
!     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
!     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
!     total expression inside the radical in Eq. 3.37.
!bob      parameter(exar=3.5,covr=2.16,c1=98.8)
real, parameter :: exar=2.5,covr=2.16,c1=116.6

data   brad,  srad/    196.,   0.047/,  &
       btlo,  stlo/   281.5,    0.26/,  &
       bthi,  sthi/   310.1,  -0.124/,  &
       bvpd,  svpd/  4850.0, -0.0051/,  &
       bswp,  sswp/ -1.07e6, 7.42e-6/
!     save

! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
! the values (zognd,0) when veg covered.

zognd = soil_rough
zoveg = veg_rough * (1.-snowfac) + zognd * snowfac
zdisp = veg_height * (1.-snowfac)
zveg = zdisp / 0.63
!bob  rasgnd = log(zts/zognd)*log((zdisp+zoveg)/zognd)/(vonk*vonk*vels)

! The following value of ground-canopy resistance for the
! nonvegetated (bare soil or water) surface is from John Garratt.
! It is 5/ustar and replaces the one from old leaf.

rasgnd = 5. / ustar

if (veg_tai >= .1 .and. snowfac < .9) then

! If vegetation is sufficiently abundant and not covered by snow, compute
! heat and moisture fluxes from vegetation to canopy, and flux resistance
! from soil or snow to canopy.

   factv = log(zts / zoveg) / (vonk * vonk * vels)
   aux = exp(exar * (1. - (zdisp + zoveg) / zveg))
   rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
   c2 = max(0.,min(1., 1.1 * veg_tai / covr))
   rd = rasgnd * (1. - c2) + rasveg * c2
   wshed = 0.
   qwshed = 0.
   transp = 0.

else

! If the TAI is very small or if snow mostly covers the vegetation, bypass
! vegetation computations.  Set heat and moisture flux resistance rd between
! the "canopy" and snow or soil surface to its bare soil value.  Set shed
! precipitation heat and moisture to unintercepted values.

   wshed = pcpgl
   qwshed = qpcpgl
   transp = 0.
   rd = rasgnd

endif

! Compute sensible heat and moisture fluxes between top soil or snow surface
! and canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow
! evaporation and dewgnd is the mass of dew that forms on the snow/soil
! surface this timestep; both are defined as always positive or zero.

hflxgc = cp * dens * (tempk(mzg+ksn) - can_temp) / rd
wflx = dens * (ground_rsat - can_rvap) / rd
dewgndflx = max(0.,-wflx)
dewgnd = dewgndflx * dtll
if (ksn == 0) then
   wflxgc = max(0.,dens * (ground_rvap - can_rvap) / rd)
else
   wflxgc = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
endif

if (veg_tai >= .1 .and. snowfac < .9) then

! If vegetation is sufficiently abundant and not covered by snow, compute
! heat and moisture fluxes from vegetation to canopy, and flux resistance
! from soil or snow to canopy.

! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
! WITHOUT the LAI factor; the LAI factor is included later in computing
! fluxes that involve rb and/or rc.

   ustar0 = max(.1,ustar)
   c4 = .01 * sqrt(ustar0 * c1)

   rb  = (1. + .5 * veg_tai) / (.01 * sqrt(ustar0 * c1))

! Soil water potential factor for stomatal control

   swp = -200.0
   nveg = nint(leaf_class)

   do k = kroot(nveg),mzg
      nsoil = nint(soil_text(k))
      slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
      if (slpotv > swp) swp = slpotv
   enddo
   swp = swp * 9810.

! Begin canopy time-split iterations

   do iter_can = 1,niter_can

! Calculate the saturation vapor pressure at TVEG

      tveg = veg_temp
      tvegc = tveg - 273.15
      es = eslf(tvegc)
      rsatveg = .622 * es / (prss - es)

! Compute mixing ratio at leaf surface using previous rc

      rc = stom_resist
      rleaf = (rb * rsatveg + rc * can_rvap) / (rb + rc)
      vpd = max((es - rleaf * prss / (.622 + rleaf)),0.)

! Evaluate 5 environmental factors and new rc

      ftlo = 1. + exp(-stlo * (tveg - btlo))
      fthi = 1. + exp(-sthi * (tveg - bthi))
      frad = 1. + exp(-srad * (rshort - brad))
      fswp = 1. + exp(-sswp * (swp - bswp))
      fvpd = 1. + exp(-svpd * (vpd - bvpd))

! 15-minute response time for stomatal conductance (must have dtlc <= 900.) 

      rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(nveg)
      rc = 1. / (1. / rc + .0011 * dtlc * (1. / rc_inf - 1. / rc))

! Limit maximum transpiration to be <= 400 w m-2 by increasing rc if necessary.

      transp_test = alvl * dens * veg_lai * (rsatveg - can_rvap) / (rb + rc)
      if (transp_test > 400.) then
         rc = (rb + rc) * transp_test * .0025 - rb
      endif      
      
      stom_resist = rc

! Flux of heat and moisture from vegetation to canopy

      stai = veg_tai * (1. - snowfac)
      slai = veg_lai * (1. - snowfac)

      hflxvc = 2. * stai * cp * dens * (tveg - can_temp) / rb
      c3 = dens * (rsatveg - can_rvap)
      if (c3 >= 0.) then
         sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66667)
         wflxvc = min(c3 * 2. * stai * sigmaw / rb,veg_water/dtlc)
         transp0 = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
      else
         wflxvc = c3 * 2. * stai / rb
         transp0 = 0.
      endif

! Update vegetation moisture from precipitation, vapor flux with canopy,
! and shedding of excess moisture.

      wtemp = veg_water - wflxvc * dtlc + pcpgc * vf
      wshed0 = max(0.,wtemp - 0.2 * stai)
      wshed = wshed + wshed0
      transp = transp + transp0
      veg_water = max(0.,wtemp - wshed0)

! Update vegetation temperature from radiative fluxes and sensible and
! latent heat transfer with canopy air.  Heat capacity of the vegetation
! is defined as hcapveg; it requires further testing and needs to be
! increased for dense vegetation.  A value of hcapveg of 3.e4 corresponds
! to the heat capacity of a layer of water 7.5 mm thick.

      tvegc = veg_temp - 273.15 + dtlcohcv  &
         * (rshort_v + rlonga_v + rlonggs_v - rlongv_gs  &
         -  rlongv_a - hflxvc - (wflxvc + transp0) * alvl)

! Exchange heat between vegetation and precipitation in implicit scheme

      qwtot = qpcpgc * vf + hcapveg * tvegc
      CALL qwtk (qwtot,pcpgc * vf,hcapveg,tvegk,fracliqv)
      veg_temp = tvegk
      fac = 4186.
      if (fracliqv <= .0001) fac = 2093.
      qwshed = qwshed + (fac * (tvegk - 273.15) + fracliqv * 334000.) * wshed0

! Update temperature and moisture of canopy.  hcapcan [J/m2/K] and
! wcapcan [kg_air/m2] are
! the heat and moisture capacities of the canopy, while c2 and c3 are deltat
! times their inverses, respectively.  They require further testing.

      can_temp = can_temp + dtlcohcc * (hflxgc + hflxvc  &
               + dens * cp * ustar * tstar * pis)

      can_rvap = can_rvap + dtlcowcc * (wflxgc - dewgndflx + wflxvc  &
               + transp0 + dens * ustar * rstar)

   enddo

else

   can_temp = can_temp + dtllohcc * (hflxgc + dens * cp * ustar * tstar * pis)
   can_rvap = can_rvap + dtllowcc * (wflxgc - dewgndflx + dens * ustar * rstar)
   
endif

return
END SUBROUTINE canopy

