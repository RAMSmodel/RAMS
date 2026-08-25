!##############################################################################
Subroutine eng_params ()

use mem_grid
use io_params

implicit none

!  Set some constants that were formerly defined in the namelist.
!  This set should be changed only for special tests or with code modification.

! Keep this list of variables in sync with broadcast_config() in mpass_init.f90

SSPCT   = 0.  ! Initial sound speed fraction (overwritten after t=0)
IMPL    = 1   ! Implicit flag for acoustic model  -  0=off, 1=on

return
END SUBROUTINE eng_params

!##############################################################################
Subroutine toptinit_user (n2,n3,ifm,topt,topzo)

use mem_grid

implicit none

integer :: i,j,n2,n3,ifm,flag_on
real :: hfwid,hfwid2,hgt
real, dimension(n2,n3) :: topt,topzo

!  This routine is the intended location for a user to customize TOPT,
!  the surface topography array.  It is called after all other types of
!  initialization of this field, so this routine has the last word.
!  By default the routine makes no change to the field.  The commented
!  lines below serve as a template for user-designed changes; the example
!  shown is the Witch of Agnesi mountain. Note that
!  this routine is called for each grid separately, so attention to the
!  current value of ngrid in this routine may be required by the user.

flag_on=0
if(flag_on==1)then

      if(ifm.eq.1)then
         hfwid=100000.
         hgt=10000.
         hfwid2=hfwid**2
         do j=1,n3
            do i=1,n2
               topt(i,j)=hgt*hfwid2/(hfwid2+xtn(i,1)**2)
               topt(i,j) = 0. + float(i)
               topzo(i,j) = 0.001 
            enddo
         enddo
         topt(5,1) = 10.
      elseif(ifm.eq.2)then
      endif

endif

return
END SUBROUTINE toptinit_user

!##############################################################################
Subroutine sstinit_user (n2,n3,ifm,seatf)

implicit none

integer :: n2,n3,ifm,i,j,flag_on
real, dimension(n2,n3) :: seatf

!  This routine is the intended location for a user to customize the
!  SEATP and SEATF arrays.  It is called after all other types of 
!  initialization of these fields, so this routine has the last word.  
!  By default the routine makes no change to the fields.  The commented 
!  lines below serve as a template for user-designed changes.   Note that 
!  this routine is called for each grid separately, so attention to the 
!  current value of ngrid in this routine may be required by the user.

flag_on=0
if(flag_on==1)then

 if (ifm .eq. 1) then
    do j = 1,n3
       do i = 1,n2
          seatf(i,j) = 0.
       enddo
    enddo
 endif

endif

return
END SUBROUTINE sstinit_user

!##############################################################################
Subroutine ndviinit_user (n2,n3,npat,ifm,veg_ndvif)

implicit none

integer :: n2,n3,ipat,npat,ifm,i,j,flag_on
real, dimension(n2,n3,npat) :: veg_ndvif

!  This routine is the intended location for a user to customize the
!  NDVIP and NDVIF arrays.  It is called after all other types of 
!  initialization of these fields, so this routine has the last word.  
!  By default the routine makes no change to the fields.  The commented 
!  lines below serve as a template for user-designed changes.   Note that 
!  this routine is called for each grid separately, so attention to the 
!  current value of ngrid in this routine may be required by the user.

flag_on=0
if(flag_on==1)then

 if (ifm .eq. 1) then
    do j = 1,n3
       do i = 1,n2
          veg_ndvif(i,j,1) = 0.
          veg_ndvif(i,j,2) = 0.
          do ipat = 3,npat
             veg_ndvif(i,j,ipat) = 0.
          enddo
       enddo
    enddo
 endif

endif

return
END SUBROUTINE ndviinit_user

!##############################################################################
Subroutine sfcinit_file_user (n2,n3,nzg,npat,ifm  &
   ,patch_area,leaf_class,soil_text)

use rconstants

implicit none

integer :: n2,n3,nzg,npat,ifm,i,j,k,ipat,flag_on
real, dimension(nzg,n2,n3,npat) :: soil_text
real, dimension(n2,n3,npat) :: patch_area,leaf_class

!  This routine is the intended location for a user to customize the
!  PATCH_AREA, leaf_class, and SOIL_TEXT arrays.  It is called after all 
!  other types of initialization of these fields, so this routine has 
!  the last word.  By default the routine makes no change to the
!  fields.  The commented lines below serve as a template for user-designed
!  changes.   Note that this routine is called for each grid separately, so
!  attention to the current value of ngrid in this routine may be required
!  by the user.

flag_on=0
if(flag_on==1)then

 if (ifm .eq. 1) then
  do j = 1,n3
    do i = 1,n2

       patch_area(i,j,1) = 0.        ! patch 1
       leaf_class(i,j,1) = 0.        ! patch 1

       patch_area(i,j,2) = 0.        ! patch 2
       leaf_class(i,j,2) = 0.        ! patch 2

       do k = 1,nzg
          soil_text(k,i,j,1) = 0.    ! patch 1
          soil_text(k,i,j,2) = 0.    ! patch 2
       enddo

    enddo
  enddo

  do ipat = 3,npat
    do j = 1,n3
       do i = 1,n2

          patch_area(i,j,ipat) = 0.
          leaf_class(i,j,ipat) = 0.

          do k = 1,nzg
             soil_text(k,i,j,ipat) = 0.
          enddo

       enddo
    enddo
  enddo

 endif

endif

return
END SUBROUTINE sfcinit_file_user

!##############################################################################
Subroutine sfcinit_nofile_user (n1,n2,n3,mzg,mzs,npat,ifm  &
   ,theta,pi0,pp,rv,seatp,seatf                       &
   ,soil_water     ,soil_energy      ,soil_text       &
   ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth  &
   ,veg_fracarea   ,veg_lai          ,veg_tai         &
   ,veg_rough      ,veg_height       ,veg_albedo      &
   ,patch_area     ,patch_rought     ,patch_roughm    &
   ,leaf_class      &
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

real :: c1,airtemp,prsv,piv,maxsoil
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
                              ,patch_rought ,patch_roughm  ,leaf_class   &
                              ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                              ,ground_rsat  ,ground_rvap   ,veg_water    &
                              ,veg_temp     ,can_rvap      ,can_temp     &
                              ,veg_ndvip    ,veg_ndvic     ,veg_ndvif

integer :: flag_on

! This routine fills the primary LEAF3 arrays for which standard RAMS
! data files do not exist with initial default values.  Many of the 
! initial values are horizontally homogeneous, although some depend on
! atmospheric conditions.  The default values assigned here may be 
! overridden by (1) specification from coarser grids or (2) specifying new 
! values in routine sfcinit_nofile_user in the file ruser.f.

!If this flag is turned on (==1) then proceed to override sfcinit_nofile
flag_on=0
if(flag_on==1)then 

c1 = .5 * cpi

! Time interpolation factor for updating SST

if (iupdsst == 0) then
   timefac_sst = 0.
else
   timefac_sst = (time - ssttime1(ifm)) / (ssttime2(ifm) - ssttime1(ifm))
endif

!--------------------------------------------------------------------------------
!Saleeby (4-8-2013): Ingest and check soil moisture, soil temperature, snow depth
!and snow water from varfiles. Also fill in water/land points when a disconnect
!exists between input dataset and RAMS grid. Soil moisture ingest will only
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
  print*,'SETTING UP SOIL AND SNOW LEVEL VARIABLES ON GRID:',IFM
  print*,'-------------------------------------------------------'
 endif
endif !end if soil moisture ingest from varfile

!Initialize surface properties
do j = 1,n3
   do i = 1,n2
      piv = c1 * (pi0(1,i,j) + pi0(2,i,j)   &
                      + pp(1,i,j) + pp(2,i,j))
      airtemp = theta(2,i,j) * piv
      prsv = piv ** cpor * p00

      patch_rought(i,j,1) = 0.001
      patch_roughm(i,j,1) = 0.001
      can_temp(i,j,1) = airtemp
      can_rvap(i,j,1) = rv(2,i,j)

      soil_energy(mzg,i,j,1) = 334000.  &
         + 4186. * (seatp(i,j) + (seatf(i,j) - seatp(i,j))  &
         * timefac_sst - 273.15)

      !Loop over patches
      do ipat = 2,npat

         nveg = nint(leaf_class(i,j,ipat))

         soil_rough(i,j,ipat) = ztrough
         patch_rought(i,j,ipat) = max(ztrough,zot(i,j))
         patch_roughm(i,j,ipat) = max(zmrough,zot(i,j))
         veg_rough(i,j,ipat) = .13 * veg_ht(nveg)

         veg_height(i,j,ipat) = veg_ht(nveg)
         veg_albedo(i,j,ipat) = albv_green(nveg)
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
            if(isfcl==2) &
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

            sfcwater_mass(k,i,j,ipat) = 0.
            sfcwater_energy(k,i,j,ipat) = 0.
            sfcwater_depth(k,i,j,ipat) = 0.

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
                  if(runsoilingest==1) then 
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
            ,soil_text  (mzg,i,j,ipat),sfcwater_energy(1,i,j,ipat)  &
            ,sfcwater_nlev(i,j,ipat),ground_rsat(i,j,ipat) &
            ,ground_rvap(i,j,ipat),can_rvap(i,j,ipat),prsv,i,j)

      enddo !end loop over patches
   enddo
enddo

endif

return
END SUBROUTINE sfcinit_nofile_user

!##############################################################################
Subroutine bubble (m1,m2,m3,i0,j0,thp,rtp)

use micphys
use mem_grid
use mem_radiate, only: irce,rce_bubl
use node_mod, only: my_rams_num, mainnum, nmachs

implicit none

integer :: m1,m2,m3,i0,j0
real, dimension(m1,m2,m3) :: thp,rtp
real :: R

integer bubtemp
integer :: i,j,k,id,jd
real acetmp1,acetmp2,acetmp3,acetmp4,acetmp5
real acetmp6,acetmp7,acetmp8,acetmp9
real bubctrx,bubctry,bubctrz
real bubradx,bubrady,bubradz
real, dimension(:,:), allocatable :: bub_rand_nums

if(ibubble==1) then
 if(print_msg) then
  print*,''
  print*,'INITIALIZING Square RAMSIN WARM BUBBLE HERE'
  print*,'On grid number=',IBUBGRD
  print*,'Bubble from I=',IBDXIA,'TO',IBDXIZ
  print*,'Bubble from J=',IBDYJA,'TO',IBDYJZ
  print*,'Bubble from K=',IBDZK1,'TO',IBDZK2
  print*,'THP perturbation=',BTHP
  print*,'RTP perturbation=',BRTP
  print*,''
 endif
 do jd = IBDYJA, IBDYJZ
   j = jd - j0
   if ((j .ge. 1) .and. (j .le. m3)) then
     do id = IBDXIA, IBDXIZ
       i = id - i0
       if ((i .ge. 1) .and. (i .le. m2)) then
         do k = IBDZK1,IBDZK2
              thp(k,i,j) = thp(k,i,j) + BTHP
              rtp(k,i,j) = rtp(k,i,j) * (1.0 + BRTP)
         enddo
       endif
     enddo
   endif
 enddo

elseif(ibubble==2) then
 if(print_msg) then
  print*,'INITIALIZING Cosine-Squared RAMSIN WARM BUBBLE HERE'
  print*,'On grid number=',IBUBGRD
  print*,'Use of I,J coordinates for Cosine-Squared bubble is'
  print*,'only accurate and valid if the domain POLE POINT'
  print*,'is the same at the CENTRAL LAT and LON'
 endif
 !Set up X location of bubble center relative to grid center
 bubctrx = deltax * ( (IBDXIA+IBDXIZ)/2.0 - NNXP(1)/2.0 )
 !Set up Y location of bubble center relative to grid center
 bubctry = deltax * ( (IBDYJA+IBDYJZ)/2.0 - NNYP(1)/2.0 )
 !Set up Z location of bubble center
 if(print_msg) then
  print*,'IBDXIA,IBDXIZ',IBDXIA,IBDXIZ
  print*,'IBDYJA,IBDYJZ',IBDYJA,IBDYJZ
  print*,'X-center,Y-center:',bubctrx,bubctry
  print*,''
 endif
 bubtemp=int(IBDZK1+IBDZK2)/2.0
 bubctrz=ZMN(bubtemp,1)
 !Set up length, width, and depth of bubble
 bubradx=(IBDXIZ-IBDXIA) * deltax * 0.5
 bubrady=(IBDYJZ-IBDYJA) * deltax * 0.5
 bubradz=(ZMN(IBDZK2,1)-ZMN(IBDZK1,1)) * 0.5
 !Set up gaussian bubble
 acetmp8=atan(1.0)*4.0/2.0 ! pi/2
 do j=1,m3
  do i=1,m2
   do k=1,m1
    acetmp1=(XMN(i+i0,1)+XMN(i+i0+1,1))*0.5
    acetmp2=(YMN(j+j0,1)+YMN(j+j0+1,1))*0.5
    acetmp3=(ZMN(k,1)+ZMN(k+1,1))*0.5
    acetmp4=(acetmp1-bubctrx)/bubradx
    acetmp4=acetmp4**2
    acetmp5=(acetmp2-bubctry)/bubrady
    acetmp5=acetmp5**2
    acetmp6=(acetmp3-bubctrz)/bubradz
    acetmp6=acetmp6**2
    if(jdim==0)acetmp5=acetmp4 !For 2D simulation let Y-dir = X-dir
    acetmp7=acetmp4+acetmp5+acetmp6
    acetmp7=SQRT(acetmp7)
    if(acetmp7.ge.1.0)then
     acetmp9=0.0
    else
     acetmp9=(COS(acetmp8*acetmp7))**2
    endif
    thp(k,i,j)=thp(k,i,j)+(BTHP*acetmp9) !pt pert
    rtp(k,i,j)=rtp(k,i,j)*(1.0+(acetmp9*BRTP)) !mix ratio pert
   enddo
  enddo
 enddo

elseif(ibubble==3) then
 if(print_msg) then
  print*,'Activating random temperature perturbation'
  print*,'at z=2, linearly-decreasing to 0K over 500m'
  print*,'On grid number=',IBUBGRD
  if(irce==1)then 
     print*,'With Max amplitude',rce_bubl
  else
     print*,'With Max amplitude in subroutine bubble.'
  endif
  print*,''
 endif

 ! SRH: parallel algorithm that will make the temperature perturbation
 ! initialization independent of the number of nodes. This is accomplished
 ! by having mainnum process do the assignment for the entire domain and
 ! broadcasting these results to the rest of the nodes.
 !
 ! In order to save memory space, do one level at a time.

 do k=2,m1
   ! select levels for temp. pert. based on altitude
   if( zt(k) <= (500.+zt(2)) ) then ! only over lowest 500 m
     ! allocate memory for entire horizontal domain
     allocate(bub_rand_nums(nnxp(ibubgrd), nnyp(ibubgrd)))

     ! if a parallel run and this is mainnum process, or if a sequential run,
     ! then calculate the temp pert for this level
     if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
       CALL calc_bub_rand_nums (nnxp(ibubgrd),nnyp(ibubgrd),bub_rand_nums)
     endif

     ! if parallel run then broadcast the random number array to everyone
     if (nmachs .gt. 1) then
       ! note the multiplication of nnxp*nnyp. this is done so that
       ! bub_rand_nums can be viewed as a 1D array in the broadcast routine
       CALL broadcast_bub_rand_nums (nnxp(ibubgrd)*nnyp(ibubgrd), bub_rand_nums)
     endif

     ! At this point bub_rand_nums covers the entire domain and is filled with
     ! random numbers (between -1 and 1) that were generated by only one
     ! process. Go through the subdomain and copy the corresponding section
     ! of bub_rand_nums when reading the random numbers.
     do j = 1, m3
       do i = 1, m2
         R = bub_rand_nums(i+i0,j+j0)

         !Default perturbations amplitude of 0.1K, largest at level 2,
         !decreased linearly over 500m.
         R = R*0.1*(500.+zt(2)-zt(k))/500.

         !RCE perturbations amplitude of "rce_bubl", largest at level 2,
         !decreased linearly over 500m. Override default.
         if(irce == 1) R = R*rce_bubl*(500.+zt(2)-zt(k))/500.

         thp(k,i,j)=thp(k,i,j) + R
         if(k==2) thp(1,i,j)=thp(k,i,j) ! set level 1 to level 2
       enddo ! i
     enddo ! j

     deallocate(bub_rand_nums)
   endif ! in lowest 500 m
 enddo ! k

endif

return
END SUBROUTINE bubble

!##################################################################################
Subroutine calc_bub_rand_nums (nx,ny,bub_rand_nums)

!Return array (entire horiz domain) filled with random numbers between -1 and 1

implicit none

integer :: nx, ny
real, dimension(nx, ny) :: bub_rand_nums
integer :: i, j
real :: rand_num

do j = 1, ny
  do i = 1, nx
    !Use lowercase "call" since this is a system call
    call random_number (rand_num)                ! random number between 0 and 1
    bub_rand_nums(i,j) = 1.0 - (2.0 * rand_num)  ! random number between -1 and 1
  enddo
enddo

END SUBROUTINE calc_bub_rand_nums

!##############################################################################
Subroutine conv_forcing (ut,vt)

use mem_grid
use mem_basic
use micphys
use node_mod
use mem_tend

implicit none
!   +------------------------------------------------------------------
!   ! This routine creates a convergenze zone for layer lifting
!   !   following Loftus et al. (2008, MWR), except convergence
!   !   zone is implemented as a tendency, similarly to
!   !   Schumacher (2009, JAS)
!   +------------------------------------------------------------------

integer :: i,j,k,ii,jj
real :: fz,camp,uprime,vprime,bot,top
real, dimension(mzp,mxp,myp) :: ut,vt

!check for termination type
if (ctmax.lt.0.0) then
 if(vertvel_max(ngrid) >= abs(ctmax)) iconv=0 !Turn off convergence
 if(iconv > 0 .and. print_msg) &
   print*,'Run Convergence Forcing, ICONV=',iconv,' Stop at W >=',abs(ctmax)
else
 if(time >= ctmax) iconv=0 !Turn off convergence
 if(iconv > 0 .and. print_msg) &
   print*,'Run Convergence Forcing, ICONV=',iconv,' Stop at Time=',ctmax
endif

if(iconv == 0) return

!calculate 1/2 wavelength boundaries if ICVERT=2
if (icvert.eq.2) then
  bot=zt(ickcent)-czrad ! m
  top=czrad+zt(ickcent) ! m
endif

! calculate camp from specified namelist settings and divergence
! (Loftus et al. 2008 eqn 7)  units: m^2/s
camp = -0.5*cdivmax / ( 1./(cxrad**2.) + 1./(cyrad**2.) )

do j=1,myp
  do i=1,mxp
    do k=2,mzp

      !create vertical decrease with height dependence - max at k=2
      if (icvert.eq.1 .and. k.le.ickmax) then
        fz=1.-(2.-k)/(2.-ickmax)
      endif
      !create 1/2 vertical wavelength elevated convergence (Schumacher 2009)
      if (icvert.eq.2 .and. (zt(k).ge.bot .and. zt(k).le.top)) then
        fz=cos(0.5*3.14159*abs((zt(k)-zt(ickcent))/czrad))**2.
      elseif (icvert.eq.2) then
        fz=0.
      endif

      !Get absolute grid points for parallel (& sequential) computation
      ii = i+mi0(ngrid)
      jj = j+mj0(ngrid)

      !create u' and v'. Units: m/s. From Loftus et al. eqns. 5 & 6
      if (iconv.eq.1) then !Gaussian in x and y, U and V
        uprime = (-2.*camp*float(ii-icicent)*deltax/cxrad**2.)*&
          exp( -1.*(float(ii-icicent)*deltax/cxrad)**2. )*&
          exp( -1.*(float(jj-icjcent)*deltax/cyrad)**2. )*fz
        vprime = (-2.*camp*float(jj-icjcent)*deltax/cyrad**2.)*&
          exp( -1.*(float(ii-icicent)*deltax/cxrad)**2. )*&
          exp( -1.*(float(jj-icjcent)*deltax/cyrad)**2. )*fz
      elseif (iconv.eq.2) then !Gaussian in x and y, U only
        uprime = (-2.*camp*float(ii-icicent)*deltax/cxrad**2.)*&
          exp( -1.*(float(ii-icicent)*deltax/cxrad)**2. )*&
          exp( -1.*(float(jj-icjcent)*deltax/cyrad)**2. )*fz
        vprime = 0.
      elseif (iconv.eq.3) then !Gaussian in x and y, V only
        uprime = 0.
        vprime = (-2.*camp*float(jj-icjcent)*deltax/cyrad**2.)*&
          exp( -1.*(float(ii-icicent)*deltax/cxrad)**2. )*&
          exp( -1.*(float(jj-icjcent)*deltax/cyrad)**2. )*fz
      elseif (iconv.eq.4) then !Gaussian in x, U only
        uprime = (-2.*camp*float(ii-icicent)*deltax/cxrad**2.)*&
          exp( -1.*(float(ii-icicent)*deltax/cxrad)**2. )*fz
        vprime = 0.
      elseif (iconv.eq.5) then !Gaussian in y, V only
        uprime = 0.
        vprime = (-2.*camp*float(jj-icjcent)*deltax/cyrad**2.)*&
          exp( -1.*(float(jj-icjcent)*deltax/cyrad)**2. )*fz
      endif ! iconv check

      ! add tendency over timescale specified in namelist
      ut(k,i,j) = ut(k,i,j) + uprime/ctau
      vt(k,i,j) = vt(k,i,j) + vprime/ctau

    enddo ! k = 2:mzp
  enddo ! i loop
enddo ! j loop

return
END SUBROUTINE conv_forcing
