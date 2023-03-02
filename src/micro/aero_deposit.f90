!##############################################################################
Subroutine deposition_driver (i,j,m1,xztop,rtgt       &
               ,xrv,xtheta,xup,xvp                    &
               ,dn0                                   &
               ,xpi0,xpp                              &
               ,xlclass,xparea                        &
               ,xustar,xprough,monthx                 &
               )        

!------------------Routine deposition_driver-----------------------------------
!This routine calculates dry deposition over water and vegetation and
!wet-precipitation scavenging via each hydrometeor species. Wet scavenging
!is allowed in-cloud and below-cloud by using 3D precipitation rates to
!apply to the calculation of scavenging coefficient. If anything we may
!still underestimate precipitation scavenging in-cloud by ice particles.
!Wet scavenging is parameterized according to Seinfeld and Pandis (2006).
!Dry deposition over water parameterized according to Slinn & Slinn (1980).
!Dry deposition over vegetation parameterized according to Slinn (1982).
!------------------------------------------------------------------------------

use mem_grid
use aero_include
use rconstants
use micphys

implicit none

! Some local constants
real, parameter::  AA = 2.53E11 , B =  5.42E3, smalle = 0.622, xubmin = 0.1

! Hydromet category distinctions
integer :: lhcat

! define extern variables (input)
integer :: k,rundep,lcat
integer :: i,j,m1,monthx
real ,dimension(m1) :: xztop,xrv,xpi0,xpp,xtheta,xup,xvp
real ,dimension(npatch)::xustar,xlclass,xparea,xprough
real :: netrough !net roughness length - patch roughness
real :: rtgt !vertical grid spacing correction over topography

! define extern variables (output)
real, dimension(8) :: acon_remove,amas_remove
real :: totalmasremove,totalconremove,removerat

!Storage term for particles settling into lower layers
real, dimension(m1) :: gaincon,gainmas

!Vertical variation of particle size and density
!rdry = dry particle radius (cm) as defined in mic_init
!rdryum = dry particle radius (microns)
!rwet = particle radius in wet conditions (microns)
real, dimension(m1) :: rdry,rdryum,rwet,dn0

real:: vd             ,&  ! dry deposition for each grid
       ttair          ,&  ! air temperature (in k)
       tmpvd          ,&  ! dry deposition for each patch
       es             ,&  ! saturated water vapor pressure
       eee            ,&  ! water pressure
       ppress         ,&  ! air pressure (pa)
       rh             ,&  ! relative humidity
       speed          ,&  ! wind speed for all levels above the first
       ustar          ,&  ! friction velocity of the surface
       picpi          ,&  ! pressure term to convert theta to temperature
       patarea        ,&  ! land patch area fraction 0.0 to 1.0
       xmu            ,&  ! air viscosity
       xlamda         ,&  ! air free path
       cc             ,&  ! slip correction factor
       diff           ,&  ! Brownian diffusion coefficient
       scarate        ,&  ! scavenging rate [s^-1]
       zthick         ,&  ! Thickness of each layer [m]
       xdeltat        ,&  ! timestep
       u_10           ,&  ! 10-meter U wind (m/s)
       v_10           ,&  ! 10-meter V wind (m/s)
       ref_10         ,&  ! wind speed (m/s) at 10m reference height
       totrr          ,&  ! Total rain rate = conprr + pcprr [mms^-1]
       Dpc            ,&  ! Hydrometer collector size (m)
       epsilonsol     ,&  ! particle solubility fraction 0.0 -> 1.0
       rhodry         ,&  ! dry particle density (kg/m3)
       rhowet             ! deliquesced particle density (kg/m3)
integer:: np,        &  ! surface patch loop indices
          abot,atop, &  ! vertical layer loop indices
          itype,     &  ! patch type
          iusetype,  &  ! sfc types (note itype 30, but iusetype only 15)
          seasoninx     ! seasonal index

!Loop deposition driver over all aerosol types
do acat=1,aerocat

 !Set default deposition flag to off until criteria are met
 rundep=0

 !Set up profile of aerosol properties if they exist
 if((acat==1 .and. iaerosol>0)  .or. &  ! CCN
    (acat==2 .and. iaerosol>0)  .or. &  ! GCCN
    (acat==3 .and. idust>0)    .or. &  ! Small dust mode
    (acat==4 .and. idust>0)    .or. &  ! Large dust mode
    (acat==5 .and. isalt>0)    .or. &  ! Salt film mode
    (acat==6 .and. isalt>0)    .or. &  ! Salt jet mode
    (acat==7 .and. isalt>0)    .or. &  ! Salt spume mode
    (acat==8 .and. iabcarb>0)  .or. &  ! Absorbing carbon 1 mode
    (acat==9 .and. iabcarb>0)  .or. &  ! Absorbing carbon 2 mode
    (acat==aerocat-1 .and. iccnlev>=2) .or. &  ! Small regenerated aerosol
    (acat==aerocat   .and. iccnlev>=2)) then   ! Large regenerated aerosol
    rundep=1
    epsilonsol = aero_epsilon(acat)

    !Initially compute dry radius (meters)
    do k = 2,m1-1
     if(aerocon(k,acat)>mincon .and. aeromas(k,acat)>=minmas) then
      rdry(k)=((0.23873/aero_rhosol(acat) &
        *aeromas(k,acat)/aerocon(k,acat))**(1./3.))/aero_rg2rm(acat)
     else
      rdry(k)=0.005e-6
     endif
    enddo
 endif

 if(rundep==1) then !Run aerosol deposition if deposition turned on

!#####Saleeby testing
!do k = 1,m1-1
!rdry(k) = 1.0e0 / 2.0 * 1.0e-6
!xup(k) = 15.0
!xvp(k) = 0.0
!enddo
!#####Saleeby testing

  !Must first convert aerosol from x/kg to x/m3 for falling
  !Similar to sedimentation of hydrometeors
  do k = 2,m1-1
    aerocon(k,acat) = aerocon(k,acat) * dn0(k)
    aeromas(k,acat) = aeromas(k,acat) * dn0(k)
    !Aerosol and solubility tracking
    if(iccnlev>=2.and.itrkepsilon==1.and.(acat==aerocat-1.or.acat==aerocat)) &
      regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) * dn0(k)
  enddo

  !Radii in mic_init are in meters: we need these to be in um
  !Also check for positive concentration and zero or negative mass
  do k = 2,m1-1
   rdryum(k) = rdry(k)*1.e6
   if((aerocon(k,acat)>mincon .and. aeromas(k,acat)<=0.) &
       .or.aeromas(k,acat)<0.0.or.aerocon(k,acat)<0.0) then
     print*,'Aerosol deposition check (ngrid,k,i,j,acat):',ngrid,k,i,j,acat
     print*,'aerocon,aeromas,r',aerocon(k,acat),aeromas(k,acat),rdry(k)
     stop
   endif
  enddo

  !Time step
  xdeltat = dtlt

  ! Check for bottom and top of aerosol layer; operations only on these layers
  abot = 2
  k = abot
  do while (k .le. m1-1 .and. aerocon(k,acat) .le. mincon)
    abot = k
    k = k + 1
  enddo
  atop = abot

  !Find top layer without assuming continuous dust in the vertical
  !This allows for discrete dust layers with clean air between
  if ((m1-1-abot) .ge. 1) then
    do while (k .le. m1-1)
      if (aerocon(k,acat) .gt. mincon) atop = k
      k = k + 1
    enddo
  endif

  !Loop wet and dry deposition over layers containing aerosol
  if ((abot .le. atop) .and. (abot .lt. m1-1)) then
   do k = abot,atop
    !calculate  layer thickness
    zthick = xztop(k) - xztop(k-1)
    !rh and temp
    picpi = (xpi0(k) + xpp(k)) * cpi
    ppress = p00 * picpi ** cpor
    ttair = xtheta(k) * picpi
    es = AA*exp(-1.*B/ttair)
    eee = ppress*xrv(k)/smalle
    rh = eee/es

    !Aerosol and solubility tracking
    if(iccnlev>=2 .and. itrkepsilon==1 .and. (acat==aerocat-1.or.acat==aerocat) &
      .and. aeromas(k,acat)>0.) &
        epsilonsol = min(1.0,regenmas(k,acat-(aerocat-2))/aeromas(k,acat))

    !Set default density (kg/m3) in k loop and update if particle deliquesced
    rhodry=aero_rhosol(acat)

!#####Saleeby testing
!rhodry = 1769. !(kg/m3)
!rh = 0.999
!epsilonsol = 1.0
!#####Saleeby testing

    !Using Fitzgerald (1975,JAM) to estimate wet particle size. His paper
    !requires 60% solubility. Paper predicts equilibrium radius above ~80%
    !RH, but we want instantaneous deliquesced size. We limit max radius
    !to 10um for use in gravitational settling equation. Return wet particle
    !size and volume weighted density of deliquesced particle.
    if (epsilonsol .lt. 0.60 .or. rh .lt. 0.75) then
      rwet(k) = rdryum(k)
      rhowet  = rhodry
    else
      CALL cal_dwet (ttair,rh,rdryum(k),rwet(k),aero_vanthoff(acat) &
                   ,epsilonsol,rhodry,rhowet)
    endif

!#####Saleeby testing
!rwet(k)=rdryum(k)
!#####Saleeby testing

    !calculate air viscosity, Schmidt number etc in ambient air
    !press is in pa, the input needs kpa. Input Rwet in microns.
    CALL calvar (ttair,ppress/1000.,rwet(k)*2.,xmu,xlamda,cc,diff)

    !Precip Scavenging for aerosols (Saleeby & van den Heever 2013)
    !Slinn (1983) method detailed in (X.Wang et al. 2010; Seinfeld & Pandis 2006)

    !Added in-cloud scavenging for aerosols that do not get nucleation scavenged.
    !Currently assuming spheres for snow/aggr scavenging
    !Literature suggests snow scavenging coefficients are higher than for
    !rain, so the spherical assumption likely sets a lower bound for snow
    !scavenging for now.

    if(level==3 .and. aerocon(k,acat) > mincon) then
     do lcat=2,8
      acon_remove(lcat)=0.0
      amas_remove(lcat)=0.0
      lhcat = jhcat(k,lcat)
      !Model level Hydromet Species precip rate (mm/s)
      totrr = pcpvx(k,lcat)
      ScaRate = 0.

!#####Saleeby testing
!totrr=2.777e-4 !1mm/hr pcpvx(k,lcat)
!rx(k,lcat)=4.192e-9 !(kg/kg) gives Dpc=.0002m,.2mmrain
!cx(k,lcat)=1. !(#/kg)
!#####Saleeby testing

      !Do scavenging if precip rate > 0.01 mm/hr (2.777e-6 mm/sec)
      if (totrr>2.777e-6 .and. rx(k,lcat)>rxmin .and. cx(k,lcat)>cxmin) then
        CALL wet_scavenge_slinn (k,rhowet,totrr,rwet(k),xmu,cc &
            ,diff,dn0(k),ScaRate,lcat,cfmas(lhcat),pwmas(lhcat) &
            ,cfvt(lhcat),pwvt(lhcat),Dpc)

!#####Saleeby testing
!if(acat==1.and.lcat==2.and.i==15.and.j==1.and.k==3)then
!print*,'Epsilonsol,RH,rhod,rhow',epsilonsol,rh,rhodry,rhowet
!print*,'Dry R,D(um),Epsilon,RH:   ',rdryum(k),rdryum(i)*2.0
!print*,'Wet R,D(um),Scavrate(/hr):',rwet(k),rwet(k)*2.0,ScaRate*3600.
!endif
!#####Saleeby testing

        acon_remove(lcat) = ScaRate*aerocon(k,acat)*xdeltat
        amas_remove(lcat) = ScaRate*aeromas(k,acat)*xdeltat
      endif
     enddo

     !Reduce scavenging if necessary so we do not over-remove aerosols
     totalconremove=0.0
     totalmasremove=0.0
     do lcat=2,8
       totalconremove = totalconremove + acon_remove(lcat)
       totalmasremove = totalmasremove + amas_remove(lcat)
     enddo
     if(totalconremove > aerocon(k,acat)) then
       removerat = 0.999 * aerocon(k,acat) / totalconremove
       do lcat=2,8
         acon_remove(lcat) = acon_remove(lcat) * removerat
       enddo
     endif
     if(totalmasremove > aeromas(k,acat)) then
       removerat = 0.999 * aeromas(k,acat) / totalmasremove
       do lcat=2,8
         amas_remove(lcat) = amas_remove(lcat) * removerat
       enddo
     endif

     !Check to see that we do not over-remove aerosols by scavenging
     totalconremove=0.0
     totalmasremove=0.0
     do lcat=2,8
       totalconremove = totalconremove + acon_remove(lcat)
       totalmasremove = totalmasremove + amas_remove(lcat)
     enddo

     !Finally remove the aerosols via wet scavenging
     do lcat=2,8
       aerocon(k,acat) = aerocon(k,acat) - acon_remove(lcat)
       aeromas(k,acat) = aeromas(k,acat) - amas_remove(lcat)
       !Aerosol and solubility tracking
       if(iccnlev>=2) then
         !Do density conversion to kg/kg for added aerosol mass to "lcat" dependent variables
         cnmhx(k,lcat) = cnmhx(k,lcat) + (amas_remove(lcat) / dn0(k))
         if(itrkepsilon==1) then
          snmhx(k,lcat) = snmhx(k,lcat) + (amas_remove(lcat) / dn0(k)) * epsilonsol
          if(acat==aerocat-1.or.acat==aerocat) &
            regenmas(k,acat-(aerocat-2))=regenmas(k,acat-(aerocat-2))-amas_remove(lcat)*epsilonsol
         endif
         if(itrkdust==1 .and. (acat==3.or.acat==4)) &
          dnmhx(k,lcat) = dnmhx(k,lcat) + (amas_remove(lcat) / dn0(k))
       endif
     enddo
    endif !wet scavenge if micro level=3 at least mincon

    !Dry deposition for lowest layer: includes vegetative effects
    Vd = 0.
    if (k .eq. 2)  then
     do np = 1,npatch
      patarea = xparea(np)

      if (patarea .gt. 0.009 ) then
        itype = nint(xlclass(np))
        ustar = max(xustar(np),xubmin)

!#####Saleeby testing
!itype = 0
!ustar = 0.43 !ocean
!ustar = 0.97 !vegetation
!#####Saleeby testing

        !Water: don't distinguish between ocean & fresh
        if (itype .eq. 0 .or. itype .eq. 1)   iusetype = 0

        !trees
        if (itype .eq. 4 .or. itype .eq. 20) iusetype = 1
        if (itype .eq. 5) iusetype = 3
        if (itype .eq. 6) iusetype = 4
        if (itype .eq. 7) iusetype = 2

        !short and tall grass are all grass
        if (itype .eq. 8  .or. itype .eq. 9  .or. &
            itype .eq. 18) iusetype = 6

        !semi desert and deserts are deserts
        if (itype .eq. 3 .or. itype .eq. 10) iusetype = 8

        !tundra
        if (itype .eq. 11) iusetype = 9

        !evergreen, deciduous and mixed shrub
        if (itype .eq. 12 .or. itype .eq. 13 .or. itype .eq. 14) iusetype = 10

        !crops/mixed farmer, irregular crop, and bog/marsh
        if (itype .eq. 15 .or. itype .eq. 16)  iusetype = 7

        !mixed cover
        !if ( itype .eq. 22) iusetype = 5

        !urban
        if ( itype .eq. 19) iusetype = 15

        !bog/marsh
        if ( itype .eq. 17) iusetype = 11

        !ice
        if ( itype .eq.  2) iusetype = 12

        if ( monthx .lt. 6 .and. monthx .ge. 3 ) seasoninx = 5 ! spring
        if ( monthx .lt. 3 .or.  monthx .gt. 11) seasoninx = 4 ! winter
        if ( monthx .ge. 6 .and. monthx .le. 8)  seasoninx = 1 ! summer
        if ( monthx .gt. 8 .and. monthx .le. 11) seasoninx = 2 ! summer
        !IF snow is present, must set seasoninx = 3

if(ISCM>=1)then
 !SaleebySCM:Override month to be similar for SCM
 seasoninx = 5
endif

        !Net roughness length for deposition velocity
        netrough = xprough(np)

if(ISCM>=1)then
 !SaleebySCM:Override net roughness length to be similar to SCM
 if(np == 1) netrough=0.005
 if(np  > 1) netrough=0.05
endif

        !Use vd_overwater for iusetype = 0 , otherwise use vd_overveg
        if (iusetype .eq. 0) then
         if (netrough .eq. 0.) then
           ref_10 = sqrt(xup(2)**2. + xvp(2)**2.)
         else
           u_10 = 1./log(zthick/netrough) * xup(1) * log(10./netrough)
           v_10 = 1./log(zthick/netrough) * xvp(1) * log(10./netrough)
           ref_10 = sqrt(u_10**2. + v_10**2.)
         endif

if(ISCM>=1)then
 !SaleebySCM:Override U* when not inputting U* from land-surface model
 ustar = max((ref_10*0.4)/(log(zthick/netrough)),xubmin)
endif

         CALL vd_overwater (k,rwet(k)*2.,rdryum(k)*2.,ustar,ref_10,xmu, &
                           cc,rhodry,rhowet,diff,dn0(k),tmpvd)
        else
         if (z0(iusetype,seasoninx) .eq. 0.) then
           ref_10 = sqrt(xup(2)**2. + xvp(2)**2.)
         else
           u_10 = 1/log(zthick/z0(iusetype,seasoninx)) * xup(1) &
                  * log(10./z0(iusetype,seasoninx))
           v_10 = 1/log(zthick/z0(iusetype,seasoninx)) * xvp(1) &
                  * log(10./z0(iusetype,seasoninx))
           ref_10 = sqrt(u_10**2. + v_10**2.)
         endif

if(ISCM>=1)then
 !SaleebySCM:Override U* when not inputting U* from land-surface model
 ustar = max((ref_10*0.4)/(log(zthick/z0(iusetype,seasoninx))),xubmin)
endif

         CALL vd_overveg (rwet(k)*2.,ref_10,ustar,xmu,cc &
                        ,diff,rhowet,iusetype,seasoninx,dn0(k),tmpvd)
        endif
        Vd = Vd + tmpVd * patarea
      endif ! if patch area is > 0
     enddo ! looping over all patches at lowest layer
    endif ! Dry deposition at lowest layer

    !Dry deposition for all layers above the lowest
    if (k .gt. 2) then
       speed = sqrt(xup(k)**2. + xvp(k)**2.)
       np = 1 !Tells code to use deposition over water for all upper layers
       ustar = max(xustar(np),xubmin)

if(ISCM>=1)then
 !SaleebySCM:Override U* when not inputting U* from land-surface model
 ustar = xubmin
endif

       CALL vd_overwater (k,rwet(k)*2.,rdryum(k)*2.,ustar,speed,xmu &
                        ,cc,rhodry,rhowet,diff,dn0(k),tmpvd)
       Vd = tmpvd
    endif !deposition for upper layers

!#####Saleeby testing
!if(acat==1.and.i==5.and.j==5.and.k==3) &
! print*,'Grav Settling V (cm/s):',Vd*100.
!if(acat==1.and.i==5.and.j==5.and.k==2) &
! print*,'Sfc  Settling V (cm/s):',Vd*100.,iusetype
!#####Saleeby testing

    !Cascade aerosols down through each vertical layer
    gaincon(k) = aerocon(k,acat)*Vd/zthick/rtgt*xdeltat
    gainmas(k) = aeromas(k,acat)*Vd/zthick/rtgt*xdeltat
    ! Limit gaincon to be <= aerocon. If gaincon > aerocon, it will
    ! cause aerocon to go negative, plus it will cause an artificial
    ! source on the k-1 layer in the loop below that adds the falling
    ! particles from the k+1 layer to the k layer.
    if (gaincon(k) .gt. aerocon(k,acat)) then
        gaincon(k) = aerocon(k,acat)
    endif
    aerocon(k,acat) = aerocon(k,acat) - gaincon(k)
    if (gainmas(k) .gt. aeromas(k,acat)) then
        gainmas(k) = aeromas(k,acat)
    endif
    aeromas(k,acat) = aeromas(k,acat) - gainmas(k)
    !Aerosol and solubility tracking
    if(iccnlev>=2.and.itrkepsilon==1.and.(acat==aerocat-1.or.acat==aerocat)) &
      regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2))-gainmas(k)*epsilonsol
   enddo !Wet and Dry deposition loop over all vertical layers
  endif !End conditional deposition over layers containing aerosols

  !Add falling particles to layer below if more than 1 layer is present
  if (abot .lt. atop) then
    if (abot .eq. 2) then
     do k = abot,atop-1
       aerocon(k,acat) = aerocon(k,acat) + gaincon(k+1)
       aeromas(k,acat) = aeromas(k,acat) + gainmas(k+1)
       !Aerosol and solubility tracking
       if(iccnlev>=2.and.itrkepsilon==1.and.(acat==aerocat-1.or.acat==aerocat) &
         .and.aeromas(k+1,acat)>0.) then 
         epsilonsol = min(1.0,regenmas(k+1,acat-(aerocat-2))/aeromas(k+1,acat))
         regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2))+gainmas(k+1)*epsilonsol
       endif
     enddo
    else
     do k = abot-1,atop-1
       aerocon(k,acat) = aerocon(k,acat) + gaincon(k+1)
       aeromas(k,acat) = aeromas(k,acat) + gainmas(k+1)
       !Aerosol and solubility tracking
       if(iccnlev>=2.and.itrkepsilon==1.and.(acat==aerocat-1.or.acat==aerocat) &
         .and.aeromas(k+1,acat)>0.) then 
         epsilonsol = min(1.0,regenmas(k+1,acat-(aerocat-2))/aeromas(k+1,acat))
         regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2))+gainmas(k+1)*epsilonsol
       endif
     enddo
    endif
  endif

  !Lastly, convert aerosol back from x/m3 to x/kg
  !Similar to sedimentation of hydrometeors
  do k = 2,m1-1
    aerocon(k,acat) = aerocon(k,acat) / dn0(k)
    aeromas(k,acat) = aeromas(k,acat) / dn0(k)
    if(aerocon(k,acat)<0. .or. aeromas(k,acat)<0.)then
      print*,'Aerosol Deposition Negative:',k,i,j,aerocon(k,acat),aeromas(k,acat)
      stop
    endif
    !Aerosol and solubility tracking
    if(iccnlev>=2.and.itrkepsilon==1.and.(acat==aerocat-1.or.acat==aerocat)) then
      regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) / dn0(k)
      if(regenmas(k,acat-(aerocat-2))<0.)then
        print*,'RegenNeg',k,i,j,regenmas(k,acat-(aerocat-2)),aeromas(k,acat)
        stop
      endif
    endif
  enddo

 endif !run deposition if turned on for aerosol type
enddo !end looping over each aerosol type

return
END SUBROUTINE deposition_driver

!##############################################################################
Subroutine cal_dwet (T,rh,ddry,dwet,vhoff,epsilonsol,rhodry,rhowet)

!This routine calculates the growth of size for ammonium sulfate
!follwoing formula Fitzgerald 1975, J. Appl. Meteo.,1044-1049.
!Equation r = A * rd ** B requires radius in microns
!r100 growth equation computes in meters and needs converstion
!to microns before applying to "ddry" which is in microns
!Then return radius in microns. ddry and dwet are actually radii
!and not diameter.

implicit none

!Extern variables
real :: t,rh,ddry,dwet,epsilonsol,rhodry,rhowet,volumeratio
integer :: vhoff

!Internal variables
real :: alpha,beta,phi,xk1,xk2,xe,estimat_dwet,rhmax

!Use rhmax for beta and alpha terms to prevent possible 
!division by zero or very small number
rhmax=min(rh,1.000)

!calculate beta
beta = exp(0.00077*rhmax/(1.009-rhmax))

! calculate phi
phi = 1.058
if (rh.gt.0.97) then
  phi = 1.058-0.0155*(rhmax-0.97)/(1.02-rhmax**1.4)
endif

!calculate alpha
alpha = 1.2*exp(0.066*rhmax/(phi-rhmax))
xk1 = 10.2 - 23.7*rhmax + 14.5*rhmax**2.
xk2 = -6.7 + 15.5*rhmax - 9.2*rhmax**2.
xe = epsilonsol
alpha = alpha*(1. - xk1*(1.-xe)-xk2*(1.-xe*xe))

if (rh .le. 0.001)  then
     dwet = ddry
endif

estimat_dwet = alpha * (ddry)**beta

! note when RH =1 the formla is (i*Mw*Rv*r0*T/(2*sigma*Ms))^0.5*rd^1.5
!   i: vant hoff constant, defined as number of inon versus total #
!      of solutions, NaCL(i=2), CaCl2(i=3), (NH4)2SO4(i=3)
!  sigma: surface tension, Newton/M, for water is 0.07
!      T: temperature (K)
!     Mw: molecular weight of water, 18g/mol
!     Ms: molecular weight of (NH4)2SO4, 132g/mol
! rhodry: density of particle, for ammonium sulfate is 1769kg/m3.
! rhowet: deliquesced solution droplet density
!     Rv: water vapor constant, 461 J/K-kg.
!    epsilonsol: epsilon solubility fraction 0.0->1.0
! Can use molecular weights in non-m.k.s. units since they are a ratio below

if (rh .ge. 1) then
 dwet=sqrt((vhoff*epsilonsol*18.*461.*rhodry*T)/(2.*0.07*132.*1.e6))*(ddry)**1.5
 if (dwet .lt. ddry) dwet = ddry
endif

! The following double checks if caculated wet radius greater than
! the one at rh = 1., then set it to 1. Because Fitzgerald 1975 only
! has the formula for r > 0.08um

if (rh.lt.1 .and. rh.gt.0.001) then
 dwet=sqrt((vhoff*epsilonsol*18.*461.*rhodry*T)/(2.*0.07*132.*1.e6))*(ddry)**1.5
 if (dwet .lt. ddry) dwet = ddry
 if (estimat_dwet .lt. dwet) dwet = max(ddry,estimat_dwet)
endif

!Limit max size to 20um diameter (10um radius) since this is used in
!gravitation settling for particles with Re<0.1 (~20micron diam particles)
!Plus this is the equilibrium size and may be overestimated for large
!particles since they take longer to reach equilibrium. What we really
!need is instantaneous deliquescent size.
dwet = min(10.,dwet)

!If particle swells, compute density of particle + water
volumeratio = min(1.00, ddry**3 / dwet**3)
rhowet = rhodry * volumeratio + 1000. * (1.-volumeratio)

if(rhowet<1000.0 .or. rhowet>2659.0 .or.rhodry<1000.0 .or. rhodry>2659.0) then
 print*,'Bad aerosol density in cal_dwet',rhodry,rhowet,ddry,dwet &
    ,volumeratio,estimat_dwet,rh,rhmax
 stop
endif

return
END SUBROUTINE cal_dwet

!##############################################################################
Subroutine calvar (T,P,dp,xmu,xlamda,cc,diff)

implicit none

 !For calculating viscosity, free path, slip correction factor
 !and Brownian diffusion coefficient

 !Extern variables
 real :: T,p,dp,xmu,xlamda,cc,diff

 !Internal Variables
 real :: tr, pr, s, xmur, xlamdar, pi, xk

 !Tr: reference temperature, S:Sutherland constant
 !xmur: air viscosity at the standard reference atm (kg/m/s)
 !both temperature and pressure are K
 Tr = 293.15 ! kelvin
 Pr = 101.3  ! kPa
 S = 110.4
 xmur=1.8203e-5 ! kg/m/s

 !see Baron and Willeke, equation 3-10,
 !in the book edits K. Willeke and P. Baron, p28. (kg/m/s)
 xmu = xmur*(Tr+S)/(T+S)*(T/Tr)**1.5

 !Seinfeld and Pandis (2006) eq.(9.7)
 !Mean free path length at reference atmosphere 0.0651um
 xlamdar = 0.0651
 !Baron and Willeke(2001) eq.(4-6)
 !Mean free path length at (t,p) in microns
 xlamda = xlamdar*(pr/p)*(T/Tr)*((1.+S/Tr)/(1.+S/T))

 !calcualte the slip correction factor CC
 !see Baron and Willeke(2001), equation 4-8
 !edit: using eq4-9 instead
 !cc = 1. + 1./p*dp*(15.60 + 7.00*exp(-0.059*p*dp))

 !slip correction factor from Seinfeld and Pandis (2006) eq.(9.34)
 cc = 1.0 + 2.0*xlamda/dp * (1.257 + 0.4*exp(-0.55*dp/xlamda))

 !Aerosol diffusivity (m2/s) from 
 !Seinfeld and Pandis (2006) eq. (9.73)
 !note dp is in unit of microns so we correct with a 1.e6 factor
 pi = 3.1415926
 !Boltzmann constant
 xk = 1.3807e-23
 !Aerosol diffusivity (m2/s)
 Diff =xk*T*Cc/(3.*pi*xmu*dp)*1.e6

return
END SUBROUTINE calvar

!##############################################################################
Subroutine vd_overwater (k,dpwet,dpdry,ustar,uh,xmu,cc,r0d,r0w,diff,rair,vd)

!***********************************************************************
! calculate the deposition velocity based on the Slinn and Slinn (1980)
!***********************************************************************

implicit none

! calcualte the dry deposition velocity based on slinn and slinn 1980.
! input:  dpwet : wet particle diameter
!         dpdry : dry particle diameter
!         uh    : wind speed at reference height
!         xmu   : air viscosity
!         xlamda: air free path
!         cc    : slip correction factor
!         diff  : diffusion coefficient
!         r0d   : dry particle density, kg/m3
!         r0w   : wet solution drop density, kg/m3
!         vgwet : wet gravitational settling
!         vgdry : dry gravitational settling
!         k   : layer, used to determine whether to include surface effects
! output:
!         vd    : dry deposition velocity

!Extern Variables
real, intent(IN) :: dpwet,dpdry,uh,xmu,cc,r0d,r0w,diff,rair,ustar
real, intent(OUT) :: vd
integer :: k

!Internal Variables
real :: cd,g,sc,vgwet,vgdry,st,xkcprime,xkdprime,xkc,xkd,kmpa

!  kmpa: 0.4 von Karmans constant
!  rair: density of air
!     g: gravity
!    sc: schmidt number
!    st: stokes number
!    cd: drag coefficient

!Constants
kmpa = 0.4
g = 9.8

!Note: Multiply gravitation settling velocity by 1.0e-12 since dpwet
!and dpdry are in microns and we need meters (squared) for m/s

if (k .gt. 2) then
  ! Wet gravitation setting velocity - Baron & Willeke(2001) eq.3-28
  vd = r0w*dpwet**2.*g*Cc/(18.*xmu)*1.0e-12
else
  ! Wet gravitation setting velocity - Baron & Willeke(2001) eq.3-28
  vgwet = r0w*dpwet**2*g*Cc/(18.*xmu)*1.0e-12
  ! Dry gravitation setting velocity - Baron & Willeke(2001) eq.3-28
  vgdry = r0d*dpdry**2*g*Cc/(18.*xmu)*1.0e-12
  ! Drag coefficient - Slinn & Slinn(1980) Text
  Cd = ustar**2 / Uh**2
  ! Stokes or impaction number - Slinn & Slinn(1980) Text
  st = vgwet/g*(ustar**2/(xmu/rair))
  ! Schmidt number - Slinn & Slinn(1980) Text
  Sc = xmu/(rair*diff)
  if (uh .le. 1.0e-10) then
    ! Gravitation setting velocity - Baron & Willeke(2001) eq.3-28
    vd = r0w*dpwet**2*g*Cc/(18.*xmu)*1.0e-12
  else
    ! kc' - Slinn & Slinn(1980) eq.6
    xkcprime = 1./(1.-kmpa)*Cd*uh
    ! kd' - Slinn & Slinn(1980) eq.5
    xkdprime = 1./kmpa*Cd*Uh*(1./sqrt(Sc)+10.**(-3./st))
    ! kc  - Slinn & Slinn(1980) eq.4                  
    xkc = xkcprime + vgdry
    ! kd  - Slinn & Slinn(1980) eq.4
    xkd = xkdprime + vgwet
    ! Deposition velocity -!Slinn & Slinn(1980) eq.4
    vd = 1./(1./xkc + 1./xkd - vgdry/(xkc*xkd))
  endif
endif

return
END SUBROUTINE vd_overwater

!##############################################################################
Subroutine vd_overveg (dpwet,uh,ustar,xmu,cc,diff,r0w,sfcinx,seasoninx,rair,vd)

!-----------------------------------------------------------------------------
! calculate the dry  deposition velocity on different surfaces
!      Vd = Vg + 1/(Ra+Rs), Vg: gratitational settling,
!                           Ra: aerodynamic  resistances
!                           Rs: surface resistance
!This routine serves deposition at the lowest layer of the model:
!Both surface and aerodynamic resistance are factored into the net
!Fallspeed
!Based on Slinn (1981)
!----------------------------------------------------------------------------

use aero_include

implicit none

! Input: dpwet:  wet particle diameter
!           uh:  wind speed
!        ustar:  friction velocity
!          xmu:  air visconsity
!           cc:  slip correction factor
!         diff:  diffusion coefficient
!          r0w:  wet solution drop density, kg/m3
!       sfcinx:  sfc type index
!    seasoninx:  seasonal index
!        vgwet:  gravitational settling of wet particle
!        vgdry:  dry particle gravitational settling
! Output:
!           vd:  dry deposition velocity

!Extern Variables
real, intent(IN) :: dpwet,uh,ustar,xmu,cc,diff,r0w,rair
real, intent(OUT) :: vd
integer, intent(IN) :: sfcinx,seasoninx

!Internal Variables
real :: g,e0,sc,vgwet,ra,eb,st,eim,ein,R1,rs,kmpa

!    kmpa: 0.4 von Karmans constant
!    Cd:  drag coefficient
!    e0: = 3 , an empirical constant
!    g: gravity
!    sc: schmidt number
!    Ra: aerodynamic resistance [ms^-1]
!    Rs: surface resistance
!    St: Stokes number
!    Eb: contribution to collection efficiency from Brownian diffusion
!    EIM: " " particle inertia/impaction
!    EIN: " " from interception
!    R1: Rebound fraction

!Constants
kmpa = 0.4
g = 9.8
e0 = 3.

!Note: Multiply gravitation settling velocity by 1.0e-12 since dpwet
!and dpdry are in microns and we need meters (squared) for m/s

   !Wet gravitation setting velocity - Baron & Willeke(2001) eq.3-28
   Vgwet = r0w*dpwet**2.*g*Cc/(18.*xmu)*1.0e-12

   Sc = xmu/(rair*diff)
   Ra = uh/(ustar**2.)
   Eb = Sc ** (-1.0 * ggamma(sfcinx,seasoninx))
   if (sfcinx.ne.8  .and. sfcinx.ne.9 .and. sfcinx.ne.12) then
     St = (Vgwet * Ustar) / (A(sfcinx,seasoninx)*g*1.0e-3)
   else
     St = (Vgwet * Ustar**2.) / (xmu/rair*g)
   endif
   Eim = (St / (alpha(sfcinx,seasoninx) + St))**2.
   if (sfcinx.ne.8 .and. sfcinx.ne.9 .and. sfcinx.ne.12) then
     Ein = 0.5 * (dpwet*1.0e-6 / (A(sfcinx, seasoninx)*1.0e-3))**2.
   else
     Ein = 0.0
   endif
   if (sfcinx.eq.0) then
     R1 = 1. !No rebound from wet surfaces
   else
     R1 = exp ( -1.*St**0.5 )
   endif
   Rs = 1./(e0*ustar*(Eb + Eim + Ein)*R1)

   Vd = Vgwet + 1./(Ra+ Rs)

return
END SUBROUTINE vd_overveg

!##############################################################################
Subroutine wet_scavenge_slinn (k,Rp,RR,Rwet,xmu,cc,diff &
          ,Rair,ScaRate,lcat,cfmasx,pwmasx,cfvtx,pwvtx,Dpc)

!--------------------------------------------------------------------
! Calculate the wet depostion below the cloud based on the scheme by
! Slinn (1983), also reference by Rosenfeld 1990s air pollutions.
!--------------------------------------------------------------------
! input parameters: particle size,temperature,pressure,rain rate
! output parameters: scavenging rate (s-1).
!--------------------------------------------------------------------

use micphys
use aero_include

implicit none

      ! externs (input)
      real Rwet    ! particle raidus after hydroscopic effect
      real Rp      ! particle density kgm-3
      real Rw      ! density of water kgm-3
      real RR      ! rain rate, mm/s
      real cc      ! correct factor
      real diff    ! diffusion coeffient due to Brownian effect
      real xmu     ! air dynamic viscosity
      real Rair    ! air density
      real cfmasx,pwmasx,cfvtx,pwvtx ! hydromet power law coefficients

      ! externs (output)
      real ScaRate  ! scavenging rate

      ! internals
      real Dpc           ! rain drop size (m)
      real Vt            ! vertical velocity, based on Dpc
      real xmuwater      ! water dynamic viscosity
      real Re            ! Reynolds #
      real Schm          ! Schmidt number
      real St            ! Stokes number
      real Stcorr        ! Stokes correction for collection efficiency
      real tau           ! relaxation time, tau*g = gravitational settling
      real Dratio        ! ratio between cloud and aerosol particle size
      real vcosityRatio  ! viscosity between water and air
      real Sstar         ! parameter used to calcualte collection efficient
      real Ecollect      ! collection coefficient.
      real vvt           ! particle gravitational velocity
      integer k          ! passing in the current vertical level
      integer lcat       ! passing in scavenging hydrometeor category 
      real collect1, collect2

      ! Get water dynamic viscosity (kg/m/s)
      ! Use a set value that does not vary with temperature
      ! Could add variability with temperature from include file
        xmuwater = dvcosity(3)*1.0e-6
      ! Density of water (kg/m3)
        Rw = 1.e3
      ! Calculate Schm = Schmidt number of collected particle
        Schm = xmu/Rair/diff
      ! Calculate the hydrometeor mean diameter (meters)
        Dpc = (rx(k,lcat)/cx(k,lcat)/cfmasx)**(1./pwmasx)

      ! Calculate Vt (m/s) from RAMS mic_init power laws
        vt = cfvtx*Dpc**pwvtx

      ! Calculate Re (Reynolds numbers), Dpc is in meters
        Re = (Dpc*vt*Rair)/(2.*xmu)

      ! Charaxteristic relaxation time of particle
      ! rwet is in unit of um, and dpc is in unit of meters
        tau = Rp * (2.0*rwet*1.0e-6)**2. * Cc / (18.*xmu)
        vvt = 1.19e-4*(rwet/2.)**2.
        St = 2. * tau * (Vt - vvt) / Dpc

      ! Ratio of particle diamater to hydrometeor diameter
        Dratio = (2.0*rwet*1.0e-6)/Dpc

      ! Viscosity ratio
        vcosityRatio = xmuwater/xmu

      ! Collection efficiency.
        Sstar = (1.2 + (1./12.)*alog(1.+Re)) / (1. + alog(1+Re))
        collect1 = (4./(Re*Schm)) &
              * ( 1.0 + 0.4*Re**0.5*Schm**0.333 + 0.16*Re**0.5*Schm**0.5 )
        collect2 = (4. * Dratio) * (1./vcosityRatio + (1. + 2.*Re**0.5)*Dratio)

        Ecollect = collect1 + collect2

        if (St.gt.Sstar) then
            Stcorr=((St-Sstar)/(St-Sstar+2./3.))**(3./2.)
            Ecollect = Ecollect + Stcorr * (Rp/Rw)**0.5
        endif
        if(Ecollect>1.0) Ecollect=1.0      

      ! Calulating wet deposition scavenging rate (/s)
      ! First convert RR in mm/s to m/s
        ScaRate = (3./2.)*Ecollect*(RR*1.0e-3)/Dpc

return
END SUBROUTINE wet_scavenge_slinn


