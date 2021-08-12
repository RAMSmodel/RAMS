!##############################################################################
Subroutine leaftw (mzg,mzs  &
   ,soil_water     ,soil_energy      ,soil_text       &
   ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth  &
   ,ustar          ,tstar            ,rstar           &
   ,veg_fracarea   ,veg_lai          ,veg_tai         &
   ,veg_rough      ,veg_height       ,leaf_class      &
   ,soil_rough     ,sfcwater_nlev    ,stom_resist     &
   ,ground_rsat    ,ground_rvap      ,veg_water       &
   ,veg_temp       ,can_rvap         ,can_temp        & 
   ,rshort,ip,i,j)

use leaf_coms
use mem_leaf
use rconstants
use mem_scratch

implicit none

integer :: mzg,mzs,ip,i,j

real, dimension(mzg) :: soil_water,soil_energy,soil_text
real, dimension(mzs) :: sfcwater_mass,sfcwater_energy,sfcwater_depth

real                 :: ustar        ,tstar         ,rstar        &
                       ,veg_fracarea ,veg_lai       ,veg_tai      &
                       ,veg_rough    ,veg_height    ,leaf_class   &
                       ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                       ,ground_rsat  ,ground_rvap   ,veg_water    &
                       ,veg_temp     ,can_rvap      ,can_temp     &
                       ,rshort

integer :: k,nveg,ksn,nsoil,ksnnew,newlayers,nlayers,kold,ktrans,nsl,kzs
integer, save :: ncall=0

!-----------------------------------------------------------------------------
! parameters for new soil heat conductivity (8/17/00):  Move to sfcdata later
real, dimension(12) :: soilcond0,soilcond1,soilcond2
data soilcond0/ 0.30, 0.30, 0.29, 0.27, 0.28, 0.28  &
               ,0.26, 0.27, 0.27, 0.25, 0.25, 0.06/
data soilcond1/ 4.80, 4.66, 4.27, 3.47, 3.63, 3.78  &
               ,2.73, 3.23, 3.32, 2.58, 2.40, 0.46/       
data soilcond2/-2.70,-2.60,-2.31,-1.74,-1.85,-1.96  &
              ,-1.20,-1.56,-1.63,-1.09,-0.96, 0.00/       
!-----------------------------------------------------------------------------

real :: stretch,thik,snden,vegfracc,qwfree,wfree          &
   ,depthgain,totsnow,qw,w,qwt,wt,soilhcap,fac,wfreeb,depthloss,soilcap    &
   ,sndenmax,sndenmin,snowmin,wtnew,wtold,wdiff,watermid,availwat,wg  &
   ,wloss,soilcond,waterfrac

real, save, dimension(20) :: thicknet
real, save, dimension(20,20) :: thick

do k = 1,mzg
   dslz   (k) = slz(k+1) - slz(k)
   dslzi  (k) = 1. / dslz(k)
   dslzidt(k) = dslzi(k) * dtll
   slzt   (k) = .5 * (slz(k) + slz(k+1))
enddo

do k = 2,mzg
   dslzt   (k) = slzt(k) - slzt(k-1)
   dslzti  (k) = 1. / dslzt(k)
   dslztidt(k) = dslzti(k) * dtll
enddo

! Initialize snow thickness scaling array

if (ncall /= 40) then
   ncall = 40
   stretch = 2.0
   do kzs = 1,mzs
      thik = 1.0
      thicknet(kzs) = 0.0
      do k = 1,(kzs+1)/2
         thick(k,kzs) = thik
         thick(kzs+1-k,kzs) = thik
         thicknet(kzs) = thicknet(kzs) + 2. * thik
         thik = thik * stretch
      enddo
      if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
      do k = 1,kzs
         thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
      enddo
   enddo
endif

nveg = nint(leaf_class)
ksn = nint(sfcwater_nlev)

! Evaluate any exchanges of heat and moisture to or from vegetation, apply
! moisture and heat changes to vegetation, and evaluate the resistance
! parameter rd between canopy air and the top soil or snow surface.

CALL canopy (mzg,mzs,ksn,nveg  &
   ,soil_water,soil_text,sfcwater_mass  &
   ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height  &
   ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class  &
   ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort)

! Compute soil and effective snow heat resistance times layer depth (rfactor).

do k = 1,mzg
   nsoil = nint(soil_text(k))
   waterfrac = soil_water(k) / slmsts(nsoil)
   soilcond = soilcond0(nsoil)  &
      + waterfrac * (soilcond1(nsoil) + waterfrac * soilcond2(nsoil))
   rfactor(k) = dslz(k) / soilcond
enddo

do k = 1,ksn
   snden = sfcwater_mass(k) / sfcwater_depth(k)
   rfactor(k+mzg) = sfcwater_depth(k)  &
      / (1.093e-3 * exp(.028 * tempk(k+mzg)) * (.030 + snden  &
      * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
enddo

! Find soil and snow internal sensible heat fluxes [W/m2]

hfluxgsc(1) = 0.

do k = 2,mzg+ksn
   hfluxgsc(k) = - (tempk(k) - tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)      
enddo

! Heat flux at soil or snow top from longwave, sensible, and
! upward latent heat fluxes [W/m^2]

hfluxgsc(ksn+1+mzg) = hflxgc + wflxgc * alvl - rlonga_gs  &
   - rlongv_gs + rlonggs_v + rlonggs_a

! Update soil Q values [J/m3] and snow Q values [J/kg] from sensible heat,
! upward water vapor (latent heat), longwave, and shortwave fluxes.
! This excludes effects of dew/frost formation, precipitation, shedding,
! and percolation.  Update top soil or snow moisture from evaporation only.

do k = 1,mzg
   soil_energy(k) = soil_energy(k) + dslzidt(k) * (hfluxgsc(k) - hfluxgsc(k+1))
enddo

soil_energy(mzg) = soil_energy(mzg) + dslzidt(mzg) * rshort_g

do k = 1,ksn
   sfcwater_energy(k) = sfcwater_energy(k) + dtll &
      * (hfluxgsc(k+mzg) - hfluxgsc(k+1+mzg) + rshort_s(k)) / sfcwater_mass(k)
enddo

if (ksn == 0) then
   soil_water(mzg) = soil_water(mzg) - 1.e-3 * wflxgc * dslzidt(mzg)
else
   sfcwater_mass(ksn) = max(0.,sfcwater_mass(ksn) - wflxgc * dtll)
endif

! New moisture, qw, and depth from dew/frost formation, precipitation,
! shedding, and percolation.  ksnnew is the layer that receives the new
! condensate that comes directly from the air above.  If there is no
! pre-existing snowcover, this is a temporary "snow" layer.

if (pcpgl + wshed + dewgnd > 1.e-9) then
   ksnnew = max(ksn,1)
   vegfracc = 1. - veg_fracarea
   qwfree = dewgnd * alvl + qpcpgl * vegfracc + qwshed * veg_fracarea
   wfree = dewgnd + pcpgl * vegfracc + wshed * veg_fracarea
   depthgain = dpcpgl * vegfracc + (dewgnd + wshed) * veg_fracarea * .001
else
   ksnnew = ksn
   qwfree = 0.
   wfree = 0.
   depthgain = 0.
endif

if (ksnnew > 0) then

! Transfer water downward through snow layers by percolation.
! Fracliq is the fraction of liquid in the snowcover or surface water.  wfree
! is the quantity of that liquid in kg/m2 which is free (not attached to
! snowcover) and therefore available to soak into the layer below).
! soilcap is the capacity of the top soil layer in kg/m2 to accept surface
! water.  wfree in the lowest snow layer is limited by this value.

   totsnow = 0.
   nsoil = nint(soil_text(mzg))

   do k = ksnnew,1,-1
      qw = sfcwater_energy(k) * sfcwater_mass(k) + qwfree
      w = sfcwater_mass(k) + wfree

! If (only) snow layer is too thin for computational stability, bring
! it to thermal equilibrium with the top soil layer by exchanging
! heat between them.

      if (ksnnew == 1 .and. sfcwater_mass(k) < 3.) then
         qwt = qw + soil_energy(mzg) * dslz(mzg)
         wt = w + soil_water(mzg) * 1.e3 * dslz(mzg)
         soilhcap = slcpd(nsoil) * dslz(mzg)
         CALL qwtk (qwt,wt,soilhcap,tempk(k+mzg),fracliq(k+mzg))
         fac = 4186.
         if (fracliq(k+mzg) <= .0001) fac = 2093.
         qw = (fac * (tempk(k+mzg) - 273.15) + fracliq(k+mzg) * 334000.) * w
         tempk(mzg) = tempk(k+mzg)
         fracliq(mzg) = fracliq(k+mzg)
         soil_energy(mzg) = (qwt - qw) * dslzi(mzg)
      else
         CALL qwtk (qw,w,100.,tempk(k+mzg),fracliq(k+mzg))
      endif

! Shed liquid in excess of a 1:9 liquid-to-ice ratio.  Limit this shed amount
! (wfreeb) in lowest snow layer to amount top soil layer can hold.

      wfreeb = max (0.,w * (fracliq(k+mzg) - .1) / 0.9)
      depthloss = wfreeb * 1.e-3
      if (k == 1) then
         soilcap = 1.e3 * max (0.,-slz(mzg) * (slmsts(nsoil) - soil_water(mzg)))
         wfreeb = min (wfreeb, soilcap)
         qwfree = wfreeb * 4186. * (tempk(k+mzg) - 193.36)
         soil_water(mzg) = soil_water(mzg) + 1.e-3 * wfreeb * dslzi(mzg)
         soil_energy(mzg) = soil_energy(mzg) + qwfree * dslzi(mzg)
      endif
      qwfree = wfreeb * 4186. * (tempk(k+mzg) - 193.36)
      sfcwater_mass(k) = w - wfreeb
      sfcwater_depth(k) = sfcwater_depth(k) + depthgain - depthloss

      totsnow = totsnow + sfcwater_mass(k)
      sfcwater_energy(k) = (qw - qwfree) / (max(1.e-4,sfcwater_mass(k)))
      sfcwater_energy(k) = max (-1.6e5, min (4.8e5, sfcwater_energy(k)))

! Temporary simple evolution of snow layer depth and density

      sfcwater_depth(k) = sfcwater_depth(k) * (1. - dtll / 1.e5)
      snden = sfcwater_mass(k) / max(1.e-6,sfcwater_depth(k))
      sndenmax = 1000.
      sndenmin = max(30.,200. * (wfree + wfreeb) / max(1.e-9,sfcwater_mass(k)))
      snden = min (sndenmax, max (sndenmin, snden))
      sfcwater_depth(k) = sfcwater_mass(k) / snden
      wfree = wfreeb
      depthgain = depthloss
   enddo

! Re-distribute snow layers to maintain prescribed distribution of mass

   if (totsnow < 1.e-9) then
      sfcwater_nlev = 0.
      sfcwater_mass(1) = 0.
      sfcwater_energy(1) = 0.
      sfcwater_depth(1) = 0.
   else
      nlayers = ksnnew
      snowmin = 3.0
      newlayers = 1
      do k = 2,mzs
         if (snowmin * thicknet(k) <= totsnow .and.  &
            sfcwater_energy(k) < 334000.) newlayers = newlayers + 1
      enddo
      newlayers = min (newlayers, mzs, nlayers+1)
      sfcwater_nlev = float(newlayers)
      kold = 1
      wtnew = 1.
      wtold = 1.
      do k = 1,newlayers
         vctr14(k) = totsnow * thick(k,newlayers)
         vctr16(k) = 0.
         vctr18(k) = 0.
10       continue
         wdiff = wtnew * vctr14(k) - wtold * sfcwater_mass(kold)
         if (wdiff > 0.) then
            vctr16(k) = vctr16(k) + wtold * sfcwater_mass(kold)  &
               * sfcwater_energy(kold)
            vctr18(k) = vctr18(k) + wtold * sfcwater_depth(kold)
            wtnew = wtnew - wtold * sfcwater_mass(kold) / vctr14(k)
            kold = kold + 1
            wtold = 1.
            if (kold <= ksn) go to 10
         else
            vctr16(k) = vctr16(k) + wtnew * vctr14(k) * sfcwater_energy(kold)
            vctr18(k) = vctr18(k) + wtnew * vctr14(k) * sfcwater_depth(kold)  &
               / max(1.e-9,sfcwater_mass(kold))
            wtold = wtold - wtnew * vctr14(k) / sfcwater_mass(kold)
            wtnew = 1.
         endif
      enddo

      do k = 1,newlayers
         sfcwater_mass(k) = vctr14(k)
         sfcwater_energy(k) = vctr16(k) / sfcwater_mass(k)
         sfcwater_depth(k) = vctr18(k)
      enddo

   endif
endif

! Compute gravitational potential plus moisture potential z + psi [m],
! liquid water content [m], and half the remaining water capacity [m].

do k = 1,mzg
   nsoil = nint(soil_text(k))
   psiplusz(k) = slzt(k) + slpots(nsoil)   &
               * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
   soil_liq(k) = dslz(k) * min(soil_water(k) - soilcp(nsoil)  &
      ,soil_water(k) * fracliq(k) )
   half_soilair(k) = (slmsts(nsoil) - soil_water(k)) * dslz(k) * .5
enddo

! Find amount of water transferred between soil layers (wflux) [m]
! modulated by the liquid water fraction

wflux(1) = 0.
wflux(mzg+1) = 0.
qwflux(1) = 0.
qwflux(mzg+1) = 0.

do k = 2,mzg
   nsoil = nint(soil_text(k))
   watermid = 0.5 * (soil_water(k) + soil_water(k-1))
   wflux(k) = dslztidt(k) * slcons1(k,nsoil)  &
      * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)  &
      * (psiplusz(k-1) - psiplusz(k)) * .5 * (fracliq(k) + fracliq(k-1))

! Limit water transfers to prevent over-saturation and over-depletion
! Compute q transfers between soil layers (qwflux) [J/m2]

   if (wflux(k) > 0.) then
      wflux(k) = min(wflux(k),soil_liq(k-1),half_soilair(k))
   else
      wflux(k) = - min(-wflux(k),soil_liq(k),half_soilair(k-1))
   endif

   qwflux(k) = wflux(k) * 4.186e6 * (tempk(k) - 193.36)

enddo

! Update soil moisture (impose minimum value of soilcp) and q value.

do k = 1,mzg
   nsoil = nint(soil_text(k))
   soil_water(k) = min(slmsts(nsoil), max(soilcp(nsoil)  &
                  ,soil_water(k)- dslzi(k) * (wflux(k+1) - wflux(k))))
   ! Limit soil moisture to the saturation level
   if (soil_water(k) > slmsts(nsoil)) then
      soil_water(k) = min(slmsts(nsoil),soil_water(k))
   endif
   soil_energy(k) = soil_energy(k) - dslzi(k) * (qwflux(k+1) - qwflux(k))
enddo

! Remove water only from moistest level in root zone for transpiration.
! Bottom k-index in root zone is kroot(nveg).  Limit soil moisture
! to values above soilcp.  More sophisticated use of root
! profile func and water extractibility func, and
! improved minimum value are being considered.
! Units of wloss are m3/m3, of transp are kg/m2/s.

wg = 0.
nsl = nint(soil_text(mzg))
ktrans = mzg
do k = kroot(nveg),mzg
   nsoil = nint(soil_text(k))
   availwat = soil_water(k) * fracliq(k)
   if (wg < availwat) then
      wg = availwat
      ktrans = k
      nsl = nsoil
   endif
enddo

wloss = min(transp * dslzidt(ktrans) * 1.e-3,wg   &
       ,soil_water(ktrans) - soilcp(nsl))

soil_water(ktrans) = soil_water(ktrans) - wloss
soil_energy(ktrans) = soil_energy(ktrans)  &
   - wloss * 4.186e6 * (tempk(ktrans) - 193.36)

! Compute ground vap mxrat for availability on next timestep; put into
! ground_rsat.
CALL grndvap (mzs,soil_energy(mzg),soil_water(mzg),soil_text(mzg)  &
   ,sfcwater_energy,sfcwater_nlev,ground_rsat            &
   ,ground_rvap,can_rvap,prss,i,j)

return
END SUBROUTINE leaftw

!##############################################################################
Subroutine grndvap (mzs,soil_energy,soil_water,soil_text,sfcwater_energy  &
   ,sfcwater_nlev,ground_rsat,ground_rvap,can_rvap,prsg,i,j)

use leaf_coms
use rconstants

implicit none

integer :: mzs
real :: soil_energy,soil_water,soil_text
real :: sfcwater_nlev,ground_rsat,ground_rvap,can_rvap,prsg
real, dimension(mzs) :: sfcwater_energy

integer :: ksn,nsoil,i,j

real :: gdrm,tempkk,fracliqq,slpotvn,alpha,beta
real, external :: rslif

gdrm = g / rm

ksn = nint(sfcwater_nlev)

! ground_rsat(i,j) is the saturation mixing ratio of the top soil/snow surface
! and is used for dew formation and snow evaporation.

if (ksn > 0) then
   CALL qtk (sfcwater_energy(ksn),tempkk,fracliqq)
   ground_rsat = rslif(prsg,tempkk)
else

! Without snowcover, ground_rvap is the effective saturation mixing
! ratio of soil and is used for soil evaporation.  First, compute the
! "alpha" term or soil "relative humidity" and the "beta" term.

   nsoil = nint(soil_text)
   
   CALL qwtk (soil_energy,soil_water*1.e3,slcpd(nsoil),tempkk,fracliqq)
   ground_rsat = rslif(prsg,tempkk)

   slpotvn = slpots(nsoil) * (slmsts(nsoil) / soil_water) ** slbs(nsoil)
   alpha = exp(gdrm * slpotvn / tempkk)
   beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nsoil)) * 3.14159)) ** 2
   ground_rvap = ground_rsat * alpha * beta + (1. - beta) * can_rvap

endif

return
END SUBROUTINE grndvap
