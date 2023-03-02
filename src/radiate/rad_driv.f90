!##############################################################################
Subroutine radiate (mzp,mxp,myp,ia,iz,ja,jz)

use mem_tend
use mem_grid
use mem_leaf
use mem_sib
use mem_radiate
use mem_basic
use mem_scratch
use mem_micro
use rconstants
use rrad3
use micphys
use ref_sounding

implicit none

integer :: mzp,mxp,myp,ia,iz,ja,jz

integer, save :: ncall=0

if (ilwrtyp + iswrtyp .eq. 0) return

CALL tend_copy (mzp,mxp,myp,ia,iz,ja,jz,radiate_g(ngrid)%fthrdp(1,1,1)  &
   ,radiate_g(ngrid)%fthrd(1,1,1))

!Proceed if this is a radiation timestep "radfrq"
if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. 0.001) then

   if(iprntstmt>=1 .and. print_msg) &
      print 90,time,time/3600.+(itime1/100+mod(itime1,100)/60.)
90      format(' Radiation Tendencies Updated Time =',F10.1,  &
        '  UTC TIME (HRS) =',F6.1)

   ! Compute solar zenith angle, multiplier for solar constant, sfc albedo,
   ! and surface upward longwave radiation. radprep provides zenith angle,
   ! surface downward shortwave and longwave for land-surface model in call
   ! to "sfcrad". Then "sfcrad" returns albedo, upward longwave, and upward
   ! shortwave from surface.

   CALL radprep (mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,jday   &

      ,leaf_g(ngrid)%soil_water      (1,1,1,1)  &
      ,leaf_g(ngrid)%soil_energy     (1,1,1,1)  &
      ,leaf_g(ngrid)%soil_text       (1,1,1,1)  &
      ,leaf_g(ngrid)%sfcwater_energy (1,1,1,1)  &
      ,leaf_g(ngrid)%sfcwater_depth  (1,1,1,1)  &
      ,leaf_g(ngrid)%leaf_class      (1,1,1)    &
      ,leaf_g(ngrid)%veg_fracarea    (1,1,1)    &
      ,leaf_g(ngrid)%veg_height      (1,1,1)    &
      ,leaf_g(ngrid)%veg_albedo      (1,1,1)    &
      ,leaf_g(ngrid)%patch_area      (1,1,1)    &
      ,leaf_g(ngrid)%sfcwater_nlev   (1,1,1)    &
      ,leaf_g(ngrid)%veg_temp        (1,1,1)    &
      ,leaf_g(ngrid)%can_temp        (1,1,1)    &
      ,solfac                                   &
      ,grid_g(ngrid)%glat            (1,1)      &
      ,grid_g(ngrid)%glon            (1,1)      &
      ,radiate_g(ngrid)%rshort       (1,1)      &
      ,radiate_g(ngrid)%rlong        (1,1)      &
      ,radiate_g(ngrid)%rlongup      (1,1)      &
      ,radiate_g(ngrid)%albedt       (1,1)      &
      ,radiate_g(ngrid)%cosz         (1,1)      &
      ,sib_g(ngrid)%sfcswa           (1,1,1)    &
      ,sib_g(ngrid)%uplwrf           (1,1,1)    )

   ! If using Mahrer-Pielke and/or Chen-Cotton radiation, call radcomp.

   if (ilwrtyp .le. 2 .or. iswrtyp .le. 2) then

      ! Zero out the radiative heating rate "fthrd" if this this a 
      ! radiation timestep and fthrd will be updated.

      CALL azero (mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))
   
      CALL radcomp (mzp,mxp,myp,ia,iz,ja,jz,solfac  &
         ,basic_g(ngrid)%theta     (1,1,1)  &
         ,basic_g(ngrid)%pi0       (1,1,1)  &
         ,basic_g(ngrid)%pp        (1,1,1)  &
         ,basic_g(ngrid)%rv        (1,1,1)  &
         ,basic_g(ngrid)%dn0       (1,1,1)  &
         ,basic_g(ngrid)%rtp       (1,1,1)  &
         ,radiate_g(ngrid)%fthrd   (1,1,1)  &
         ,grid_g(ngrid)%rtgt       (1,1)    &
         ,grid_g(ngrid)%f13t       (1,1)    &
         ,grid_g(ngrid)%f23t       (1,1)    &
         ,grid_g(ngrid)%glon       (1,1)    &
         ,radiate_g(ngrid)%rshort  (1,1)    &
         ,radiate_g(ngrid)%rlong   (1,1)    &
         ,radiate_g(ngrid)%albedt  (1,1)    &
         ,radiate_g(ngrid)%cosz    (1,1)    &
         ,radiate_g(ngrid)%rlongup (1,1))

   endif

   ! Using Harrington radiation

   if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then

      ! If first call for this node, initialize several quantities & Mclatchy
      ! sounding data.

      if (ncall .eq. 0) then
         CALL radinit (ng,nb,nsolb,npsb,nuum,prf,alpha,trf,beta  &
            ,xp,wght,wlenlo,wlenhi,solar0,ralcs,a0,a1,a2,a3  &
            ,exptabc,ulim,npartob,npartg,ncog,ncb  &
            ,ocoef,bcoef,gcoef,gnu)

         CALL mclatchy (1,mzp  &
            ,grid_g(ngrid)%glat       (1,1)  &
            ,grid_g(ngrid)%rtgt       (1,1)  &
            ,grid_g(ngrid)%topt       (1,1)  &
            ,radiate_g(ngrid)%rlongup (1,1)  &
            ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
            ,vctr8,vctr9,vctr10,vctr11,vctr12 &
            )

         ncall = ncall + 1
      endif

      ! For any call, interpolate the mclatchy sounding data by latitude and
      ! season.

      CALL mclatchy (2,mzp  &
         ,grid_g(ngrid)%glat       (1,1)  &
         ,grid_g(ngrid)%rtgt       (1,1)  &
         ,grid_g(ngrid)%topt       (1,1)  &
         ,radiate_g(ngrid)%rlongup (1,1)  &
         ,zm,zt,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7  &
         ,vctr8,vctr9,vctr10,vctr11,vctr12 &
         )

      ! If using Harrington radiation with moisture complexity LEVEL < 3,
      ! call radcomp3 which is a substitute driving structure to the bulk
      ! microphysics.

      if (level <= 2 .or. level==4) then
           ! Zero out the radiative heating rate "fthrd" if this this a 
           ! radiation timestep and fthrd will be updated.
           CALL azero (mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))
           ! Run the Harrington radiation for non-LEVEL=3 micro
           CALL radcomp3 (mzp,mxp,myp,ia,iz,ja,jz  &
            ,grid_g(ngrid)%glat       (1,1)    &
            ,grid_g(ngrid)%rtgt       (1,1)    &
            ,grid_g(ngrid)%topt       (1,1)    &
            ,radiate_g(ngrid)%albedt  (1,1)    &
            ,radiate_g(ngrid)%cosz    (1,1)    &
            ,radiate_g(ngrid)%rlongup (1,1)    &
            ,radiate_g(ngrid)%rshort  (1,1)    &
            ,radiate_g(ngrid)%rlong   (1,1)    &
            ,radiate_g(ngrid)%aodt    (1,1)    &
            ,basic_g(ngrid)%rv        (1,1,1)  &
            ,basic_g(ngrid)%dn0       (1,1,1)  &
            ,radiate_g(ngrid)%fthrd   (1,1,1)  &
            ,basic_g(ngrid)%pi0       (1,1,1)  &
            ,basic_g(ngrid)%pp        (1,1,1)  &
            ,basic_g(ngrid)%theta     (1,1,1)  &
            ,micro_g(ngrid)%rcp       (1,1,1)  &
            ,radiate_g(ngrid)%bext    (1,1,1)  &
            ,radiate_g(ngrid)%swup    (1,1,1)  &
            ,radiate_g(ngrid)%swdn    (1,1,1)  &
            ,radiate_g(ngrid)%lwup    (1,1,1)  &
            ,radiate_g(ngrid)%lwdn    (1,1,1)  &
            ,micro_g(ngrid)%cccnp     (1,1,1)  &
            ,micro_g(ngrid)%cccmp     (1,1,1)  &
            ,micro_g(ngrid)%gccnp     (1,1,1)  &
            ,micro_g(ngrid)%gccmp     (1,1,1)  &
            ,micro_g(ngrid)%md1np     (1,1,1)  &
            ,micro_g(ngrid)%md1mp     (1,1,1)  &
            ,micro_g(ngrid)%md2np     (1,1,1)  &
            ,micro_g(ngrid)%md2mp     (1,1,1)  &
            ,micro_g(ngrid)%salt_film_np (1,1,1)  &
            ,micro_g(ngrid)%salt_film_mp (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_np  (1,1,1)  &
            ,micro_g(ngrid)%salt_jet_mp  (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_np (1,1,1)  &
            ,micro_g(ngrid)%salt_spum_mp (1,1,1)  &
            ,micro_g(ngrid)%abc1np     (1,1,1)  &
            ,micro_g(ngrid)%abc1mp     (1,1,1)  &
            ,micro_g(ngrid)%abc2np     (1,1,1)  &
            ,micro_g(ngrid)%abc2mp     (1,1,1))
      endif

   endif

endif

return
END SUBROUTINE radiate

!##############################################################################
Subroutine tend_copy (m1,m2,m3,ia,iz,ja,jz,at,at2)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: at,at2

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         at(k,i,j) = at2(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE tend_copy

!##############################################################################
Subroutine radprep (m2,m3,mzg,mzs,np,ia,iz,ja,jz,jday   &
   ,soil_water       ,soil_energy      ,soil_text      &
   ,sfcwater_energy  ,sfcwater_depth   ,leaf_class     &
   ,veg_fracarea     ,veg_height       ,veg_albedo     &
   ,patch_area                                         &
   ,sfcwater_nlev    ,veg_temp         ,can_temp       & 
   ,solfac,glat,glon,rshort,rlong,rlongup,albedt,cosz  &
   ,sfcswa,uplwrf)

use rconstants
use mem_leaf, only: isfcl
use mem_grid, only: initial,time

implicit none

integer :: m2,m3,mzg,mzs,np,ia,iz,ja,jz,jday,runsfcrad
real :: solfac
real, dimension(m2,m3) :: glat,glon,rshort,rlong,rlongup,albedt,cosz
real, dimension(mzg,m2,m3,np) :: soil_water,soil_energy,soil_text
real, dimension(mzs,m2,m3,np) :: sfcwater_energy,sfcwater_depth
real, dimension(m2,m3,np) :: leaf_class,veg_fracarea,veg_height,veg_albedo  &
   ,patch_area,sfcwater_nlev,veg_temp,can_temp,sfcswa,uplwrf

integer :: ip,i,j

! Compute solar zenith angle [cosz(i,j)] & solar constant factr [solfac].

CALL zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

! Compute patch-averaged surface albedo [albedt(i,j)] and up longwave
! radiative flux [rlongup(i,j)].

CALL azero2 (m2*m3,albedt,rlongup)

do ip = 1,np

   !Only run sfcrad if using LEAF3 or it is the first model timestep 
   !or if its a water patch.
   !If running SIB land surface, use SiB surface albedo and surface
   !upward longwave after first timestep.
   runsfcrad=0
   if    (isfcl<=1) then
      runsfcrad=1
   elseif(isfcl==2) then
      if(ip==1 .or. (initial<=2 .and. time < 0.001)) runsfcrad=1
   endif

   do j = ja,jz
      do i = ia,iz

       if(runsfcrad==1) &
         CALL sfcrad (mzg,mzs,ip                                  &
          ,soil_energy    (1,i,j,ip) ,soil_water      (1,i,j,ip)  &
          ,soil_text      (1,i,j,ip) ,sfcwater_energy (1,i,j,ip)  &
          ,sfcwater_depth (1,i,j,ip) ,patch_area      (i,j,ip)    &
          ,can_temp       (i,j,ip)   ,veg_temp        (i,j,ip)    &
          ,leaf_class     (i,j,ip)   ,veg_height      (i,j,ip)    &
          ,veg_fracarea   (i,j,ip)   ,veg_albedo      (i,j,ip)    &
          ,sfcwater_nlev  (i,j,ip)   ,rshort          (i,j)       &
          ,rlong          (i,j)      ,albedt          (i,j)       &
          ,rlongup        (i,j)      ,cosz            (i,j)       )

       if(isfcl==2 .and. ip>=2) then
         albedt(i,j)  = albedt(i,j)  + patch_area(i,j,ip)*sfcswa(i,j,ip)
         rlongup(i,j) = rlongup(i,j) + patch_area(i,j,ip)*uplwrf(i,j,ip)
       endif

      enddo
   enddo
enddo

return
END SUBROUTINE radprep

!##############################################################################
Subroutine radcomp (m1,m2,m3,ia,iz,ja,jz,solfac  &
   ,theta,pi0,pp,rv,dn0,rtp,fthrd  &
   ,rtgt,f13t,f23t,glon,rshort,rlong,albedt,cosz,rlongup)

use mem_grid
use mem_scratch
use mem_radiate
use rconstants
use node_mod, only:mi0,mj0

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jday,oyr,omn,ody,otm
integer, external :: julday

real :: solfac,cdec,declin,dzsdx,dzsdy,dlon,a1,a2,dayhr,gglon,otmf  &
   ,dayhrr,hrangl,sinz,sazmut,slazim,slangl,cosi,pisolar,eqt_julian_nudge
real, dimension(nzpmax) :: rvr,rtr,dn0r,pird,prd,fthrl,dzmr,dztr,fthrs
real, dimension(nzpmax+1) :: temprd
real, dimension(m1,m2,m3) :: theta,pi0,pp,rv,dn0,rtp,fthrd
real, dimension(m2,m3) :: rtgt,f13t,f23t,glon,rshort,rlong,cosz  &
   ,albedt,rlongup
integer :: i,j,k

! Find the hour angle and julian day, then get cosine of zenith angle.
CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)
otmf=float(otm)
dayhr = (otmf-mod(otmf,10000.))/10000. + (mod(otmf,10000.)/6000.)
jday = julday(omn,ody,oyr)

! Find day of year equation of time adjustment (Saleeby2008: improve accuracy)
pisolar=3.1415927
eqt_julian_nudge=0.0
if(jday>=1 .and. jday<=106) then
 eqt_julian_nudge = -14.2 * sin(pisolar * (jday +   7) / 111.)
elseif(jday>=107 .and. jday<=166) then
 eqt_julian_nudge =   4.0 * sin(pisolar * (jday - 106) /  59.)
elseif(jday>=167 .and. jday<=246) then
 eqt_julian_nudge =  -6.5 * sin(pisolar * (jday - 166) /  80.)
elseif(jday>=247 .and. jday<=366) then
 eqt_julian_nudge =  16.4 * sin(pisolar * (jday - 247) / 113.)
endif
eqt_julian_nudge = eqt_julian_nudge / 60.0
! cdec - cosine of declination
declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
cdec = cos(declin)

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         pird(k) = (pp(k,i,j) + pi0(k,i,j)) / cp
         temprd(k) = theta(k,i,j) * pird(k)
         rvr(k) = max(0.,rv(k,i,j))
         rtr(k) = max(rvr(k),rtp(k,i,j))
         !Convert the next 7 variables to cgs for now.
         prd(k) = pird(k) ** cpor * p00 * 10.
         dn0r(k) = dn0(k,i,j) * 1.e-3
         dzmr(k) = dzm(k) / rtgt(i,j) * 1.e-2
         dztr(k) = dzt(k) / rtgt(i,j) * 1.e-2

         if(rvr(k)    <   0. .or.  &
            rtr(k)    <   0. .or.  &
            prd(k)    <   0. .or.  &
            dn0r(k)   <   0.) then
            print*, 'Negative value of density, vapor, or pressure'
            print*, 'when calling Chen-Cotton/Mahrer-Pielke radiation'
            print*, 'at (k,i,j),ngrid = ',k,i+mi0(ngrid),j,ngrid+mj0(ngrid)
            print*, 'rvr(k)  rtr(k)  prd(k)  dn0r(k)'
            print*, rvr(k), rtr(k), prd(k), dn0r(k)
            print*, 'stopping model'
            stop 'radiation call'
         endif

        if (temprd(k) < 160.) then
            print*, 'Temperature too low when calling'
            print*, 'Chen-Cotton/Mahrer-Pielke radiation'
            print*, 'at k,i,j = ',k,i+mi0(ngrid),j+mj0(ngrid)
            print*, 'ngrid,tempk=',ngrid,temprd(k)
            !stop
        endif

      enddo
      temprd(1) = (rlongup(i,j) / stefan) ** 0.25
      temprd(m1+1) = temprd(m1)

! Call the longwave parameterizations.

      if (ilwrtyp .eq. 2) then
         CALL lwradp (m1,temprd,rvr,dn0r,dztr,pird,vctr1,fthrl,rlong(i,j))
      elseif (ilwrtyp .eq. 1) then
         CALL lwradc (m1+1,rvr,rtr,dn0r,temprd,prd,pird,dztr,fthrl,rlong(i,j)  &
            ,vctr1,vctr2,vctr3,vctr4,vctr5,vctr6,vctr7,vctr8,vctr9,vctr10  &
            ,vctr11,vctr12,vctr13,vctr14,vctr15)
      endif

! The shortwave parameterizations are only valid if the cosine
!    of the zenith angle is greater than .03 .

      if (cosz(i,j) .gt. .03) then

         if (iswrtyp .eq. 2) then
            CALL shradp (m1,rvr,dn0r,dzmr,vctr1,pird,cosz(i,j)  &
               ,albedt(i,j),solar*1e3*solfac,fthrs,rshort(i,j))
         elseif (iswrtyp .eq. 1) then
            CALL shradc (m1+1,rvr,rtr,dn0r,dztr,prd,pird,vctr1  &
              ,albedt(i,j),solar*1.e3*solfac,cosz(i,j),fthrs,rshort(i,j))
         endif

! Modify the downward surface shortwave flux by considering
!    the slope of the topography.

         if (itopo .eq. 1) then
            dzsdx = f13t(i,j) * rtgt(i,j)
            dzsdy = f23t(i,j) * rtgt(i,j)

! The y- and x-directions must be true north and east for
! this correction. the following rotates the model y/x
! to the true north/east.

!Saleeby(2013):Check into this statement below
! The following rotation seems to be incorrect, so call this instead:
! routine uvtoueve (u,v,ue,ve,qlat,qlon,polelat,polelon)

            dlon = (polelon - glon(i,j)) * pi180
            a1 = dzsdx*cos(dlon) + dzsdy * sin(dlon)
            a2 = -dzsdx*sin(dlon) + dzsdy * cos(dlon)
            dzsdx = a1
            dzsdy = a2

            gglon = glon(i,j)
            if (lonrad .eq. 0) gglon = centlon(1)
!Saleeby(2009): Add eqt_julian_nudge to improve hour angle accuracy
            dayhrr = mod(dayhr+gglon/15.+24.,24.) + eqt_julian_nudge
            hrangl = 15. * (dayhrr - 12.) * pi180
            sinz = sqrt(1. - cosz(i,j) ** 2)
            sazmut = asin(max(-1.,min(1.,cdec*sin(hrangl)/sinz)))
            if (abs(dzsdx) .lt. 1e-15) dzsdx = 1.e-15
            if (abs(dzsdy) .lt. 1e-15) dzsdy = 1.e-15
            slazim = 1.571 - atan2(dzsdy,dzsdx)
            slangl = atan(sqrt(dzsdx*dzsdx+dzsdy*dzsdy))
            cosi = cos(slangl) * cosz(i,j) + sin(slangl) * sinz  &
               * cos(sazmut-slazim)
!Saleeby(08-10-2008): Check for cosi greater than zero
            if (cosi > 0.) then
               rshort(i,j) = rshort(i,j) * cosi / cosz(i,j)
            else
               rshort(i,j) =  0.
            endif
         endif

      else
         do k = 1,nzp
            fthrs(k) = 0.
         enddo
         rshort(i,j) = 0.
      endif
      
      ! Add fluxes
      do k = 2,m1-1
         fthrd(k,i,j) = fthrl(k) + fthrs(k)
      enddo

! Convert the downward flux at the ground to SI.

      rshort(i,j) = rshort(i,j) * 1.e-3 / (1. - albedt(i,j))
      rlong(i,j) = rlong(i,j) * 1.e-3
      fthrd(1,i,j) = fthrd(2,i,j)

   enddo
enddo

return
END SUBROUTINE radcomp

!##############################################################################
Subroutine radcomp3 (m1,m2,m3,ia,iz,ja,jz  &
   ,glat,rtgt,topt,albedt,cosz,rlongup,rshort,rlong,aodt  &
   ,rv,dn0,fthrd,pi0,pp,theta,rcp &
   ,bext,swup,swdn,lwup,lwdn &
   ,cccnp,cccmp,gccnp,gccmp,md1np,md1mp,md2np,md2mp &
   ,salt_film_np,salt_film_mp,salt_jet_np,salt_jet_mp &
   ,salt_spum_np,salt_spum_mp,abc1np,abc1mp,abc2np,abc2mp)

use mem_grid
use mem_micro
use mem_radiate
use rconstants
use rrad3
use micphys
use micro_prm, only:iceprocs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,mcat,i,j,k

real :: cfmasi,cparmi,glg,glgm,picpi
real, dimension(m2,m3) :: glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong,aodt
real, dimension(m1,m2,m3) :: dn0,rv,fthrd,pi0,pp,theta,rcp
real, dimension(m1,m2,m3) :: bext,swup,swdn,lwup,lwdn
real, dimension(m1,m2,m3) :: cccnp,cccmp,gccnp,gccmp,md1np,md1mp,md2np,md2mp &
  ,salt_film_np,salt_film_mp,salt_jet_np,salt_jet_mp,salt_spum_np,salt_spum_mp &
  ,abc1np,abc1mp,abc2np,abc2mp
real, external :: gammln

! Fill cloud parameters if not running microphysics
cparmi = 1. / cparm
if (level <= 1) then
   mcat = 0
elseif (level == 2) then
   mcat = 1
   pwmas(1) = 3.
   pwmasi(1) = 1. / pwmas(1)
   cfmasi = 1. / 524.
   emb0(1) = 5.24e-16
   emb1(1) = 3.35e-11
   glg  = gammln(gnu(1))
   glgm = gammln(gnu(1) + pwmas(1))
   dnfac(1) = (cfmasi * exp(glg - glgm)) ** pwmasi(1)
   do k = 2,m1-1
      jhcat(k,1) = 1
   enddo
endif

! Loop over columns

do j = ja,jz
   do i = ia,iz

      ! To be used in mclatchy call
      do k = 1,m1-1
         picpi = (pi0(k,i,j) + pp(k,i,j)) * cpi
         press(k) = p00 * picpi ** cpor
         tair(k) = theta(k,i,j) * picpi
      enddo

      ! Finish cloud parameter computation
      if (level == 2) then
         do k = 2,m1-1
            emb(k,1) = max(emb0(1),min(emb1(1),rcp(k,i,j) * cparmi))
            cx(k,1) = rcp(k,i,j) / emb(k,1)
         enddo
      elseif (level ==4) then
         CALL cloud_prep_lev4 (m1,i,j,ngrid)
         if (iceprocs == 1) then
            mcat = 11 
         else
            mcat = 2 
         endif
      endif

      ! Call the sub-driver
      if(iaerorad==1 .and. level .ne. 4)then
        CALL aero_copy (1,m1 &
           ,cccnp(1,i,j),cccmp(1,i,j) &
           ,gccnp(1,i,j),gccmp(1,i,j) &
           ,md1np(1,i,j),md1mp(1,i,j) &
           ,md2np(1,i,j),md2mp(1,i,j) &
           ,salt_film_np(1,i,j),salt_film_mp(1,i,j) &
           ,salt_jet_np(1,i,j) ,salt_jet_mp(1,i,j)  &
           ,salt_spum_np(1,i,j),salt_spum_mp(1,i,j) &
           ,abc1np(1,i,j),abc1mp(1,i,j) &
           ,abc2np(1,i,j),abc2mp(1,i,j))
      endif

      CALL radcalc3 (m1,i,j,ngrid,maxnzp,mcat,iswrtyp,ilwrtyp,zm,zt &
         ,glat(i,j),rtgt(i,j),topt(i,j),rv(1,i,j) &
         ,albedt(i,j)          &
         ,cosz(i,j)            &
         ,rlongup(i,j)         &
         ,rshort(i,j)          &
         ,rlong(i,j)           &
         ,aodt(i,j)            &
         ,fthrd(1,i,j)         &
         ,bext(1,i,j)          &
         ,swup(1,i,j)          &
         ,swdn(1,i,j)          &
         ,lwup(1,i,j)          &
         ,lwdn(1,i,j)          &
         ,dn0(1,i,j)           &
         )

   enddo
enddo

return
END SUBROUTINE radcomp3

!##############################################################################
Subroutine zen (m2,m3,ia,iz,ja,jz,jday,glat,glon,cosz,solfac)

use mem_grid
use mem_radiate
use rconstants

implicit none

integer :: m2,m3,ia,iz,ja,jz,jday,i,j,oyr,omn,ody,otm
integer, external :: julday

real :: solfac,sdec,cdec,declin,d0,d02,dayhr,radlat,cslcsd,snlsnd  &
   ,gglon,dayhrr,hrangl,pisolar,eqt_julian_nudge,otmf
real, dimension(m2,m3) :: glat,glon,cosz

! Find the hour angle and julian day, then get cosine of zenith angle.
CALL date_add_to (iyear1,imonth1,idate1,itime1*100,time,'s',oyr,omn,ody,otm)

! Fix the julian day of the year for RCE simulations for const radiation
if(irce == 1) then
 CALL date_add_to (iyear1,imonth1,idate1,itime1*100,0.,'s',oyr,omn,ody,otm)
endif

otmf=float(otm)
dayhr = (otmf-mod(otmf,10000.))/10000. + (mod(otmf,10000.)/6000.)
jday = julday(omn,ody,oyr)

!      sdec - sine of declination, cdec - cosine of declination
declin = -23.5 * cos(6.283 / 365. * (jday + 9)) * pi180
sdec = sin(declin)
cdec = cos(declin)

! Find the factor, solfac, to multiply the solar constant to correct
! for Earth's varying distance to the sun.

d0 = 6.2831853 * float(jday-1) / 365.
d02 = d0 * 2.
solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
   + 0.000719 * cos(d02) + 0.000077 * sin(d02)

! Find day of year equation of time adjustment (Saleeby2008: improve accuracy)
pisolar=3.1415927
eqt_julian_nudge=0.0
if(jday>=1 .and. jday<=106) then
 eqt_julian_nudge = -14.2 * sin(pisolar * (jday +   7) / 111.)
elseif(jday>=107 .and. jday<=166) then
 eqt_julian_nudge =   4.0 * sin(pisolar * (jday - 106) /  59.)
elseif(jday>=167 .and. jday<=246) then
 eqt_julian_nudge =  -6.5 * sin(pisolar * (jday - 166) /  80.)
elseif(jday>=247 .and. jday<=366) then
 eqt_julian_nudge =  16.4 * sin(pisolar * (jday - 247) / 113.)
endif
eqt_julian_nudge = eqt_julian_nudge / 60.0

do j = ja,jz
   do i = ia,iz
      radlat = glat(i,j) * pi180
      if (lonrad .eq. 0) radlat = centlat(1) * pi180
      if (radlat .eq. declin) radlat = radlat + 1.e-5
      cslcsd = cos(radlat) * cdec
      snlsnd = sin(radlat) * sdec
      gglon = glon(i,j)
      if (lonrad .eq. 0) gglon = centlon(1)

      !Saleeby(2009): Add eqt_julian_nudge to improve hour angle accuracy
      dayhrr = mod(dayhr+gglon/15.+24.,24.) + eqt_julian_nudge
      hrangl = 15. * (dayhrr - 12.) * pi180

      !Cosine of zenith angle (angle in radians)
      cosz(i,j) = snlsnd + cslcsd * cos(hrangl)

      !Hold zenith angle constant for RCE simulations
      if(irce == 1) cosz(i,j) = cos(rce_szen * pi180)

   enddo
enddo

return
END SUBROUTINE zen

!##############################################################################
Subroutine radcalc3 (m1,i,j,ngrid,maxnzp,mcat,iswrtyp,ilwrtyp,zm,zt &
   ,glat,rtgt,topt,rv,albedt,cosz,rlongup,rshort,rlong,aodt &
   ,fthrd,bext,swup,swdn,lwup,lwdn &
   ,dn0 &
   )

!-----------------------------------------------------------------------------
! radcalc3: column driver for two-stream radiation code
! variables used within routine radcalc3:
! ==================================================================
! Variables in rrad3 parameter statement
!  mb               : maximum allowed number of bands [=8]
!  mg               : maximum allowed number of gases [=3]
!  mk               : maximum number of pseudobands allowed for any gas [=7]
!  ncog             : number of fit coefficients (omega and asym) [=5]
!  ncb              : number of fit coefficients (extinction) [=2]
!  npartob          : number of hydrometeor categories (including different habits)
!  npartg           : number of hydrometeor categories used for gc coefficients [=7]
!  nrad             : total number of radiation levels used (m1 - 1 + narad)
!  narad            : number of radiation levels added above model
!  nsolb            : active number of solar bands
!  nb               : active number of bands
!  ng               : active number of gases
!  jday             : julian day
!  solfac           : solar constant multiplier for variable E-S distance
!  ralcs (mb)       : rayleigh scattering integration constants
!  solar1 (mb)      : solar fluxes at top of atmosphere - corrected for ES distance
!  solar0 (mb)      : solar fluxes at top of atmosphere - uncorrected for ES distance
!  nuum (mb)        :    continuum flags
!  a0,a1,a2,a3 (mb) : Planck func fit coefficients
!  npsb (mg,mb)     : number of pseudo bands
!  trf (mg,mb)      : reference temperature for xp and wght coefficients
!  prf (mg,mb)      : reference pressure for xp and wght coefficients
!  ulim (mg,mb)     : upper bound on pathlength for gases
!  xp (mg,mk,mb)    : coefficient used in computing gaseous optical depth
!  alpha (mg,mk,mb) : pressure correction factor exponent
!  beta (mg,mk,mb)  : temperature correction factor exponent
!  wght (mg,mk,mb)  : pseudo band weight
!  exptabc (150)    : table of exponential func values
!  ocoef(ncog,mb,npartob)  : fit coefficients for hyd. single scatter.
!  bcoef(ncb,mb ,npartob)  : fit coefficients for hyd. extinction coefficient.
!  gcoef(ncog,mb,npartg)   : fit coefficients for hyd. asymmetry parameter.
!
! Input variables from model
!
!  m1               : number of vertical levels in model grid
!  ncat             : max number of hydrometeor categories [=7]
!  mcat             : actual number of hydrometeor categories [= 0, 1, or 7]
!  nhcat            : number of hydrometeor categories including ice habits [=15]
!  iswrtyp          : shortwave radiation parameterization selection flag
!  ilwrtyp          : longwave radiation parameterization selection flag
!  glat             : latitude
!  rtgt             : terrain-following coordinate metric factor
!  topt             : topography height
!  albedt           : surface albedo
!  cosz             : solar zenith angle
!  rlongup          : upward longwave radiation at surface (W/m^2)
!  rshort           : downward shortwave radiation at surface (W/m^2)
!  rlong            : downward longwave radiation at surface (W/m^2)
!  aodt             : total aerosol optical depth (band=3)
!  jnmb (ncat)      : microphysics category flag
!  dnfac (nhcat)    : factor for computing dn from emb
!  pwmasi (nhcat)   : inverse of mass power law exponent for hydrometeors
!  zm (m1)          : model physical heights of W points (m)
!  zt (m1)          : model physical heights of T points (m)
!  press (nzpmax)   : model pressure (Pa)
!  tair (nzpmax)    : model temperature (K)
!  rv (m1)          : model vapor mixing ratio (kg/kg)
!  dn0 (m1)         : model air density (kg/m^3)
!  fthrd (m1)       : theta_il tendency from radiation
!  jhcat (nzpmax,ncat)  : microphysics category array
!  cx (nzpmax,ncat) : hydrometeor number concentration (#/kg)
!  emb (nzpmax,ncat): hydrometeor mean mass (kg)
!
! Variables input from model scratch space (redefined internally on each call)
!
!  zml (nrad)       : physical heights of W points of all radiation levels (m)
!  ztl (nrad)       : physical heights of T points of all radiation levels (m)
!  dzl (nrad)       : delta-z (m) of all radiation levels
!  pl (nrad)        : pressure (Pa)
!  tl (nrad)        : temperature (K)
!  dl (nrad)        : air density of all radiation levels (kg/m^3)
!  rl (nrad)        : vapor density of all radiation levels (kg/m^3)
!  vp (nrad)        : vapor pressure (Pa)
!  o3l (nrad)       : stores the calculated ozone profile (g/m^3)
!  flxu (nrad)      : Total upwelling flux (W/m^2)
!  flxd (nrad)      : Total downwelling flux (W/m^2)
!  t (nrad)         : layer transmission func
!  r (nrad)         : layer reflection func
!  tc (nrad)        : cumulative optical depth
!  sigu (nrad)      : upwelling layer source func
!  sigd (nrad)      : downwelling layer source func
!  re (nrad)        : cumulative reflection func
!  vd (nrad)        : multi-scat. diffuse downwelling contributions
!                         from source func
!  td (nrad)        : inverse of cumulative transmission fnct
!  vu (nrad)        : multi-scat. diffuse upwelling contributions
!                         from source func
!  tg (nrad)        : gaseous optical depth
!  tcr (nrad)       : continuum/Rayleigh optical depth
!  src (nrad)       : Planck func source for each band
!  fu (nrad*6)      : upwelling fluxes for pseudo-bands (W/m^2)
!  fd (nrad*6)      : downwelling fluxes for pseudo-bands (W/m^2)
!  u (nrad*mg)      : path-length for gases (H_2O, CO_2, O_3)  (Pa)
!  tp (nrad*mb)     : optical depth of hydrometeors (m^-1)
!  omgp (nrad*mb)   : Single scatter albedo of hydrometeors
!  gp (nrad*mb)     : Asymmetry factor of hydrometeors
!
! Locally-defined variables
!
!  ngass (mg) : flags indicating if H20,CO2,O3 are active for solar wavelengths
!  ngast (mg) : flags indicating if H20,CO2,O3 are active for long wavelengths
!
! Additional variables used only within routine mclatchy:
! ==================================================================
! namax : maximum allowed number of added rad levels above model top[=10]
!                       used for oc and bc coefficients [=13]
! mcdat (33,9,6)    : Mclatchy sounding data (33 levels, 9 soundings, 6 vars)
! mclat (33,9,6)    : mcdat interpolated by season to latitude bands
! mcol (33,6)       : mclat interpolated to lat-lon of grid column
!
! Additional variables used only within routine cloud_opt:
! ==================================================================
!  ib .......... band number
!  dn .......... characteristic diameter (m)
!  oc .......... scattering albedo fit coefficients
!  bc .......... extinction fit coefficients
!  gc .......... asymmetery fit coefficients
!  kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                   numbers as a func of 15 microphysics category numbers
!
! Particle Numbers describe the following particles:
!
!     Harrington radiation code             RAMS microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   rain                        2.  rain
!  3:   pristine ice columns        3.  pristine ice columns
!  4:   pristine ice rosettes       4.  snow columns
!  5:   pristine ice plates         5.  aggregates
!  6:   snow columns                6.  graupel
!  7:   snow rosettes               7.  hail
!  8:   snow plates                 8.  pristine ice hexagonal plates
!  9:   aggregates columns          9.  pristine ice dendrites
!  10:  aggregates rosettes        10.  pristine ice needles
!  11:  aggregates plates          11.  pristine ice rosettes
!  12:  graupel                    12.  snow hexagonal plates
!  13:  hail                       13.  snow dendrites
!                                  14.  snow needles
!                                  15.  snow rosettes
!
! for the asymmetery parameter, since we only have spherical
! particles, there are only 7 particle types...
!  1:   cloud drops
!  2:   rain
!  3:   pristine ice
!  4:   snow
!  5:   aggregates
!  6:   graupel
!  7:   hail
!---------------------------------------------------------------------------

use rconstants
use rrad3
use micphys
use node_mod

implicit none

integer m1,maxnzp,mcat,ngrid
integer :: iswrtyp,ilwrtyp
integer i,j,k
integer, save :: ncall = 0,nradmax
integer, save :: ngass(mg)=(/1, 1, 1/),ngast(mg)=(/1, 1, 1/)
!     one can choose the gases of importance here,
!       ngas = 1    gas active
!            = 0    gas not active
!
!       ngas(1) = H2O
!       ngas(2) = CO2
!       ngas(3) =  O3

real, save :: eps=1.e-15
real :: glat,rtgt,topt,cosz,albedt,rlongup,rshort,rlong,aodt
real :: zm(m1),zt(m1),dn0(m1),rv(m1),fthrd(m1)
real :: bext(m1),swup(m1),swdn(m1),lwup(m1),lwdn(m1)

real, allocatable, save, dimension(:)     :: zml,ztl,dzl,pl,tl,dl,rl,o3l  &
                                      ,vp,flxus,flxds,tg,tcr,src,t,r,tc  &
                                      ,flxul,flxdl  &
                                      ,sigu,sigd,re,vd,td,vu  &
                                      ,u,fu,fd,tp,omgp,gp

real :: exner(m1) !non-dimensional pressure

!Saleeby(2011):Variables for radiatively active aerosols
real :: relh(m1)
real, external :: rslf

if (ncall == 0) then
   ncall = 1
   nradmax = maxnzp + namax
   allocate(zml  (nradmax) ,ztl  (nradmax) ,dzl  (nradmax) ,pl (nradmax)  &
           ,tl   (nradmax) ,dl   (nradmax) ,rl   (nradmax) ,o3l(nradmax)  &
           ,vp   (nradmax) ,flxus(nradmax) ,flxds(nradmax) ,tg (nradmax)  &
           ,flxul(nradmax) ,flxdl(nradmax)                                &
           ,tcr  (nradmax) ,src  (nradmax) ,t    (nradmax) ,r  (nradmax)  &
           ,tc   (nradmax) ,sigu (nradmax) ,sigd (nradmax) ,re (nradmax)  &
           ,vd   (nradmax) ,td   (nradmax) ,vu   (nradmax))
   allocate(u(nradmax*mg),fu(nradmax*6),fd(nradmax*6))
   allocate(tp(nradmax*mb),omgp(nradmax*mb),gp(nradmax*mb))
   tg=0.
endif

nrad = m1 - 1 + narad

! rlongup used to set tl(1): stephan*tl^4=rlongup
 CALL mclatchy (3,m1  &
   ,glat,rtgt,topt,rlongup  &
   ,zm,zt,press,tair,dn0,rv,zml,ztl,pl,tl,dl,rl,o3l,dzl &
   )

! calculate non-dimensional pressure
do k=1,m1
  exner(k) = (press(k)*p00i)**rocp
enddo

! zero out scratch arrays
 CALL azero (nrad*mg,u)
 CALL azero (nrad*6,fu)
 CALL azero (nrad*6,fd)
 CALL azero (nrad*mb,tp)
 CALL azero (nrad*mb,omgp)
 CALL azero (nrad*mb,gp)
 CALL azero (nrad   ,tg)

! Saleeby(2011): Aerosol radiative impacts section
! This must be run before cloud_opt
! Only run this for level=3 microphysics
if(iaerorad==1 .and. level .ne. 4) then
 !Only do levels 1 through m1-1. If there are no additional radiation levels,
 !then pl and tl only get inititalized (mclatchy call with iaction == 3) up 
 !to level m1-1. Also, aerorad always computes from level 2 to m1-1.
 do k=1,m1-1
   relh(k) = rv(k)/rslf(pl(k),tl(k))
 enddo
 CALL aerorad (i,j,mb,nb,nrad,m1,dzl,relh,tp,omgp,gp,dn0,aodt)
endif

! Compute hydrometeor optical properties
 CALL cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp &
   ,ocoef,bcoef,gcoef,ncog,ncb,npartob,npartg)

! Sum up attenutation by aerosols and hydrometeors
 CALL sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

! Get the path lengths for the various gases...
 CALL path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

do k = 1,nrad
   if (rl(k) <   0. .or.  &
       dl(k) <   0. .or.  &
       pl(k) <   0. .or.  &
      o3l(k) <   0.) then
      print*, 'Negative value of density, vapor, pressure, or ozone'
      print*, 'when calling Harrington radiation'
      print*, 'at k,i,j = ',k,i+mi0(ngrid),j+mj0(ngrid)
      print*, 'ngrid=',ngrid
      print*, 'stopping model'
      print*, 'rad: rl(k), dl(k), pl(k), o3l(k)'
      print*, rv(k), dl(k), pl(k), o3l(k)
      stop
   endif
enddo
do k = 1,nrad
   if (tl(k) < 160.) then
      print*, 'Temperature too low when calling Harrington radiation' 
      print*, 'at k,i,j = ',k,i+mi0(ngrid),j+mj0(ngrid)
      print*, 'ngrid,tempk=',ngrid,tl(k)
      !stop
   endif
enddo   

! call shortwave and longwave schemes...

if (iswrtyp == 3 .and. cosz > 0.03) then
   CALL azero2 (nrad,flxus,flxds)

   CALL swrad (nrad,nb,nsolb,npsb,       & !counters
      u,pl,tl,dzl,                       & !model variables
      xp,alpha,beta,wght,prf,trf,ralcs,  & !band specifics
      solar1,ngass,                      & !coefficients
      albedt,cosz,                       & !boundaries
      tp,omgp,gp,fu,fd,                  & !optical properties  
      flxus,flxds,ulim)                    !sw fluxes

   rshort = flxds(1)

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxds(k) - flxds(k-1) + flxus(k-1) - flxus(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      swup(k) = flxus(k)
      swdn(k) = flxds(k)
    enddo
    !lower and upper boundary conditions on swup and swdn
    swup(1) = flxus(1)
    swup(m1) = flxus(nrad) ! use the top radiation value rather than m1 value
    swdn(1) = flxds(1)
    swdn(m1) = flxds(nrad) ! use the top radiation value rather than m1 value

else

   swup   = 0.
   swdn   = 0.
   rshort = 0.

endif

if (ilwrtyp == 3) then
   CALL azero2 (nrad,flxul,flxdl)

   CALL lwrad (i,j,nrad,nb,nsolb,npsb,nuum,   & !counters
      u,pl,tl,vp,                             & !model variables
      xp,alpha,beta,wght,prf,trf,ralcs,       & !band specifics
      a0,a1,a2,a3,                            & !coefficients
      exptabc,ngast,                          & !boundaries
      tp,omgp,gp,fu,fd,flxul,flxdl,ulim)        !fluxes, output

   !Set rlong to surface level downward longwave flux.
   rlong = flxdl(1)

   !Make lowest level upward longwave flux equal to rlongup
   !produced from land surface models (LEAF,SiB).
   flxul(1) = rlongup

   do k = 2,m1-1
      !divide by exner to get potential temp heating rate
      fthrd(k) = fthrd(k)  &
         + (flxdl(k) - flxdl(k-1) + flxul(k-1) - flxul(k)) &
            / (dl(k) * dzl(k) * cp * exner(k))
      lwup(k) = flxul(k)
      lwdn(k) = flxdl(k)
   enddo
   !lower and upper boundary conditions on lwup and lwdn
   lwup(1) = flxul(1)
   lwup(m1) = flxul(nrad) ! use the top radiation value rather than m1 value
   lwdn(1) = flxdl(1)
   lwdn(m1) = flxdl(nrad) ! use the top radiation value rather than m1 value

endif

return
END SUBROUTINE radcalc3

!##############################################################################
Subroutine cloud_opt (mb,nb,nrad,m1,mcat,dzl  &
   ,dn0,tp,omgp,gp,oc,bc,gc,ncog,ncb,npartob,npartg)

! computing properties of spherical liquid water and irregular ice
! using fits to adt theory
!
! ib .......... band number
! mb .......... maximum number of bands
! nb .......... total number of bands
! m1 .......... number of vertical levels
! dzl .......... delta z in each level (m)
! dn .......... characteristic diameter (m)
! emb ......... mean hydrometeor mass (kg)
! cx .......... hydrometeor concentration (#/kg)
! tp .......... optical depth
! omgp ........ scattering albedo
! gp .......... asymmetry parameter
! oc .......... scattering albedo fit coefficients
! bc .......... extinction fit coefficients
! gc .......... asymmetry fit coefficients
! ncog ........ number of fit coefficients (omega and asym)
! ncb ......... number of fit coefficients (extinction)
! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
!                 numbers as a func of 15 microphysics category numbers
! dn0 ......... model air density (kg/m^3)
! dnfac ....... factor for computing dn from emb
! pwmasi ...... inverse of power used in mass power law
! npartob ..... number of hydrometeor categories (including different habits)
!                 used for oc and bc coefficients
! npartg ...... number of hydrometeor categories used for gc coefficients
!
! Saleeby(2008): would need to modify dnmin,dnmax,kradcat for drizzle

use micphys
use micro_prm

implicit none

integer mb,nb,ib,nrad,m1,ncog,ncb,krc,npartob,npartg
integer icat,mcat,k,ihcat

integer kradcat(15)
real dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),dn0(m1)
real oc(ncog,mb,npartob),bc(ncb,mb,npartob),gc(ncog,mb,npartg)
real ext,om,gg,dn,gnu_tab

real dnmin(7),dnmax(7)
data dnmin /   1.,   10.,   1.,  125.,   10.,   10.,   10./
data dnmax /1000.,10000., 125.,10000.,10000.,10000.,10000./

data kradcat/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/

!For LEVEL=4
real dnb(nkr)
integer newcat, kr
integer,dimension(11) :: kradcat3,ncat3
real, dimension(11) :: cf,pw

data kradcat3 /1,2,3,6,5,8,5,8,10,12,13/
data ncat3 /1,2,3,4,3,4,3,4,5,6,7/
data pw /3.0,3.0,3.0,2.841,2.474,2.474,2.523,2.038,2.584,3.0,3.0/
data cf /523.4,523.4,22.1,2.924,0.884,0.884,1.124,0.023,4.843,209.4,471.2/

if (level <=3) then
   do icat = 1,mcat 
      if (jnmb(icat) .gt. 0) then
 
         do k = 2,m1-1

            if (cx(k,icat) .gt. 1.e-9) then
               ihcat = jhcat(k,icat)
               krc = kradcat(ihcat)
               dn = dnfac(ihcat) * emb(k,icat) ** pwmasi(ihcat) * 1.e6

               !Modification by Adele Igel (Oct 2020).
               !Use an effective Dn here constrainted by effective diameter.
               !Gnu from microphysics can be different from Gnu allowed in 
               !radiation due to lookup tables available only for Gnu=1,2.
               !The two lines below help mitigate this difference.
               !See Harrington's thesis for more detail.
               gnu_tab = real(max(1,min(2,nint(gnu(icat)))))
               dn = dn * (gnu(icat)+2.) / (gnu_tab + 2.)

               dn = max(dnmin(icat),min(dnmax(icat),dn))

               do ib = 1,nb

                  ext = cx(k,icat) * dn0(k) * dzl(k)  &
                     * bc(1,ib,krc) * dn ** bc(2,ib,krc)
                  om = oc(1,ib,krc)  &
                     + oc(2,ib,krc) * exp(oc(3,ib,krc) * dn)  &
                     + oc(4,ib,krc) * exp(oc(5,ib,krc) * dn)
                  gg = gc(1,ib,icat)  &
                     + gc(2,ib,icat) * exp(gc(3,ib,icat) * dn)  &
                     + gc(4,ib,icat) * exp(gc(5,ib,icat) * dn)

                  tp(k,ib)   = tp(k,ib)   + ext
                  omgp(k,ib) = omgp(k,ib) + om * ext
                  gp(k,ib)   = gp(k,ib)   + gg * om * ext
  
               enddo

            endif 
         enddo
      endif
   enddo
elseif (level == 4) then
!icat    limits      krc   newcat  species
!1      1,krdrop-1    1       1      cloud
!2      krdrop,nkr    2       2      rain
!3      1, 13         3       3      small columns
!4      14, nkr       6       4      large columns
!5      1, 15         5       3      small plates
!6      16, nkr       8       4      large plates
!7      1, 15         5       3      small dendrites
!8      16, nkr       8       4      large dendrites
!9      1, nkr        10      5      aggregates
!10     1, nkr        12      6      graupel
!11     1, nkr        13      7      hail

!Cutoff size for bulk micro between small and large for 
!ice crystals is 125 microns
!Newcat is corresponding bulk micro category

   do icat = 1,mcat
      do k = 2,m1-1
         if (cxb2(k,icat) .gt. 1.e-9) then
            krc = kradcat3(icat)
            newcat = ncat3(icat)
            dn = (rxb2(k,icat)/(cxb2(k,icat)*cf(icat)))**(1./pw(icat)) &
                 *1.e6/gnu(newcat)
            dn = max(dnmin(newcat), min(dnmax(newcat),dn))
            !This only really works for liquid - 
            !need to consider how to calculate 
            !characteristic diameter for ice species - 
            !can we use the bulk table values?
            do ib = 1,nb
               ext = cxb2(k,icat) * dn0(k) * dzl(k)  &
                  * bc(1,ib,krc) * dn ** bc(2,ib,krc)
               om = oc(1,ib,krc)  &
                  + oc(2,ib,krc) * exp(oc(3,ib,krc) * dn)  &
                  + oc(4,ib,krc) * exp(oc(5,ib,krc) * dn)
               gg = gc(1,ib,newcat)  &
                  + gc(2,ib,newcat) * exp(gc(3,ib,newcat) * dn)  &
                  + gc(4,ib,newcat) * exp(gc(5,ib,newcat) * dn)
               tp(k,ib)   = tp(k,ib)   + ext
               omgp(k,ib) = omgp(k,ib) + om * ext
               gp(k,ib)   = gp(k,ib)   + gg * om * ext
            enddo
         endif
     enddo
   enddo
endif

return
END SUBROUTINE cloud_opt

!##############################################################################
Subroutine sum_opt (mb,nrad,nb,m1,tp,omgp,gp,bext,dzl)

implicit none

integer nb,ib,m1,k,mb,nrad
real tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),bext(m1),dzl(m1)

! Combine the optical properties....

do ib = 1,nb
   do k = 2,m1-1

      if (tp(k,ib) .gt. 0.0) then
         gp(k,ib) = gp(k,ib) / omgp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'gp'
         ! from being unphysical. Needed for RCE simulations.
         gp(k,ib) = MAX(MIN(gp(k,ib),1.),0.)
         omgp(k,ib) = omgp(k,ib) / tp(k,ib)
         !Saleeby(2010): Need this min/max func to prevent 'omgp'
         ! from being unphysical. Needed for RCE simulations.
         omgp(k,ib) = MAX(MIN(omgp(k,ib),1.),0.)
      else
         omgp(k,ib) = 0.0
         gp(k,ib) = 0.0
      endif

      !Check for validity of opt values before calling radiation
      if (tp(k,ib) .lt. 0) then
         print*, 'tp(k,ib) less than zero for k,ib = ',k,ib
         print*, 'tp(k,ib) = ',tp(k,ib)
         stop 'opt1'
      endif
      if (omgp(k,ib) .lt. 0. .or. omgp(k,ib) .gt. 1.) then
         print*, 'omgp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'omgp(k,ib) = ',omgp(k,ib)
         stop 'opt2'
      endif
      if (gp(k,ib) .lt. 0. .or. gp(k,ib) .gt. 1.) then
         print*, 'gp(k,ib) out of range [0,1] for k,ib = ',k,ib
         print*, 'gp(k,ib) = ',gp(k,ib)
         stop 'opt3'
      endif

   enddo
enddo

! Calculating visual range (km) using the Koschmeider equation
! Units: tp(k,ib) : optical thickness in level k, band ib (dimensionless)
!        dzl(k) : [m]
!        final bext(k) : [km]
! Consider only over band #3 (245-700 nm)
 do k = 2,m1-1
  bext(k) = 0.
  bext(k) = tp(k,3) / dzl(k) !Compute extinction coefficient
  !Prevent infinite visibility: constrain max vis to 1000km
  if(bext(k) .lt. 3.912e-6) bext(k) = 3.912e-6
  bext(k) = 3.912/bext(k)/1000.
 enddo

return
END SUBROUTINE sum_opt

!##############################################################################
Subroutine path_lengths (nrad,u,rl,dzl,dl,o3l,vp,pl,eps)

! Get the path lengths for the various gases...

implicit none

integer :: nrad
real, dimension(nrad) :: rl,dzl,dl,o3l,vp,pl
real, dimension(nrad,3) :: u
real :: rvk0,rvk1,dzl9,rmix,eps
integer :: k

u(1,1) = .5 * (rl(2) + rl(1)) * 9.81 * dzl(1)
u(1,2) = .5 * (dl(2) + dl(1)) * .45575e-3 * 9.81 * dzl(1)
u(1,3) = o3l(1) * 9.81 * dzl(1)

rvk0 = rl(1)
do k = 2,nrad
   rvk1 = (rl(k) + 1.e-6)
   dzl9 = 9.81 * dzl(k)
   rmix = rvk1 / dl(k)
   vp(k) = pl(k) * rmix / (.622 + rmix)
   u(k,1) = (rvk1 - rvk0) / (log(rvk1 / rvk0) + eps) * dzl9
   u(k,2) = (dl(k) - dl(k-1)) / (log(dl(k) / dl(k-1)) + eps)  &
       * dzl9 * 0.45575e-3
   u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
   rvk0 = rvk1
enddo

return
END SUBROUTINE path_lengths

!##############################################################################
Subroutine cloud_prep_lev4 (m1,i,j,ng)

use mem_micro
use module_hujisbm, only: xl, xi, xs, xg, xh
use micro_prm, only: col, iceprocs, nkr, &
                     krdrop, krpris, rxb2, cxb2
use micphys, only:ipris, igraup, ihail

implicit none

integer :: m1,i,j,ng,kr
integer :: k,lcat

! Prepare cloud arrays for radiation. Requires LEVEL == 4

rxb2 = 0.
cxb2 = 0.

! fill scratch arrays for cloud water
do kr=1,krdrop-1
   do k = 2,m1-1
      rxb2(k,1) = rxb2(k,1) + micro_g(ng)%ffcd(k,i,j,kr) * col
      cxb2(k,1) = cxb2(k,1) + &
                  micro_g(ng)%ffcd(k,i,j,kr) * col / xl(kr) * 1000.
   enddo
enddo

! fill scratch arrays for rain water
do kr=krdrop,nkr
   do k = 2,m1-1
      rxb2(k,2) = rxb2(k,2) + micro_g(ng)%ffcd(k,i,j,kr) * col
      cxb2(k,2) = cxb2(k,2) + &
                  micro_g(ng)%ffcd(k,i,j,kr) * col / xl(kr) * 1000.
   enddo
enddo

if (iceprocs == 1) then
! fill scratch arrays for ice columns
if (ipris == 1 .or. ipris >= 4) then
do kr=1,krpris(1)-1
   do k = 2,m1-1
      rxb2(k,3) = rxb2(k,3) + micro_g(ng)%ffic(k,i,j,kr) * col
      cxb2(k,3) = cxb2(k,3) + &
                  micro_g(ng)%ffic(k,i,j,kr) * col / xi(kr,1) * 1000.
   enddo
enddo
do kr=krpris(1),nkr
   do k = 2,m1-1
      rxb2(k,4) = rxb2(k,4) + micro_g(ng)%ffic(k,i,j,kr) * col
      cxb2(k,4) = cxb2(k,4) + &
                  micro_g(ng)%ffic(k,i,j,kr) * col / xi(kr,1) * 1000.
   enddo
enddo
endif

! fill scratch arrays for ice plates

if (ipris == 2 .or. ipris >=4) then
do kr=1,krpris(2)-1
   do k = 2,m1-1
      rxb2(k,5) = rxb2(k,5) + micro_g(ng)%ffip(k,i,j,kr) * col
      cxb2(k,5) = cxb2(k,5) + &
                  micro_g(ng)%ffip(k,i,j,kr) * col / xi(kr,2) * 1000.
   enddo
enddo
do kr=krpris(2),nkr
   do k = 2,m1-1
      rxb2(k,6) = rxb2(k,6) + micro_g(ng)%ffip(k,i,j,kr) * col
      cxb2(k,6) = cxb2(k,6) + &
                  micro_g(ng)%ffip(k,i,j,kr) * col / xi(kr,2) * 1000.
   enddo
enddo
endif

! fill scratch arrays for ice dendrites

if (ipris >= 3) then
do kr=1,krpris(3)-1
   do k = 2,m1-1
      rxb2(k,7) = rxb2(k,7) + micro_g(ng)%ffid(k,i,j,kr) * col
      cxb2(k,7) = cxb2(k,7) + &
                  micro_g(ng)%ffid(k,i,j,kr) * col / xi(kr,3) * 1000.
   enddo
enddo
do kr=krpris(3),nkr
   do k = 2,m1-1
      rxb2(k,8) = rxb2(k,8) + micro_g(ng)%ffid(k,i,j,kr) * col
      cxb2(k,8) = cxb2(k,8) + &
                  micro_g(ng)%ffid(k,i,j,kr) * col / xi(kr,3) * 1000.
   enddo
enddo
endif

! fill scratch arrays for aggregates

do kr=1,nkr
   do k = 2,m1-1
      rxb2(k,9) = rxb2(k,9) + micro_g(ng)%ffsn(k,i,j,kr) * col
      cxb2(k,9) = cxb2(k,9) + &
                  micro_g(ng)%ffsn(k,i,j,kr) * col / xs(kr) * 1000.
   enddo
enddo

! fill scratch arrays for graupel

if (igraup > 0) then
do kr=1,nkr
   do k = 2,m1-1
      rxb2(k,10) = rxb2(k,10) + micro_g(ng)%ffgl(k,i,j,kr) * col
      cxb2(k,10) = cxb2(k,10) + &
                   micro_g(ng)%ffgl(k,i,j,kr) * col / xg(kr) * 1000.
   enddo
enddo
endif

! fill scratch arrays for hail

if (ihail > 0) then
do kr=1,nkr
   do k = 2,m1-1
      rxb2(k,11) = rxb2(k,11) + micro_g(ng)%ffhl(k,i,j,kr) * col
      cxb2(k,11) = cxb2(k,11) + &
                   micro_g(ng)%ffhl(k,i,j,kr) * col / xh(kr) * 1000.
   enddo
enddo
endif
endif !if Iceprocs

return
END SUBROUTINE cloud_prep_lev4
