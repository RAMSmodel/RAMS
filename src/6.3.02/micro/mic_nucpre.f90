!##############################################################################
Subroutine prenuc_ccn (k,i,j)

use micphys
use mem_grid

implicit none

integer :: k,i,j

!*********** AEROSOL PROPERTIES CHECKING ******************************
do acat=1,aerocat

  !Set default values to override if aerosol type exists
  aero_rg(acat) = aero_medrad(acat) ! Default median radius 
      
  if((acat==1)                  .or. &  ! CCN
     (acat==2)                  .or. &  ! GCCN
     (acat==3 .and. idust>0)    .or. &  ! Small dust mode
     (acat==4 .and. idust>0)    .or. &  ! Large dust mode
     (acat==5 .and. isalt>0)    .or. &  ! Salt film mode
     (acat==6 .and. isalt>0)    .or. &  ! Salt jet mode
     (acat==7 .and. isalt>0)    .or. &  ! Salt spume mode
     (acat==8 .and. iabcarb>0)  .or. &  ! Absorbing carbon 1 mode
     (acat==9 .and. iabcarb>0)  .or. &  ! Absorbing carbon 2 mode
     (acat==aerocat-1 .and. iccnlev>=2) .or. &  ! Small regenerated aerosol
     (acat==aerocat   .and. iccnlev>=2)) then   ! Large regenerated aerosol

     !Assign aerosol specs to local arrays
     aeromass   = aeromas(k,acat)

     !Keep median radius and aerosol mass in bounds
     if(aerocon(k,acat) > mincon .and. aeromas(k,acat) > minmas) then

       rhosol=aero_rhosol(acat)
       aero_rg(acat)=((0.23873/rhosol*aeromas(k,acat)/aerocon(k,acat)) &
                    **(1./3.))/aero_rg2rm(acat)

       if(aero_rg(acat) < 0.01e-6) aero_rg(acat) = 0.01e-6
       if(aero_rg(acat) > 6.50e-6) aero_rg(acat) = 6.50e-6

       aeromas(k,acat) = ((aero_rg(acat)*aero_rg2rm(acat))**3.) &
                       *aerocon(k,acat)/(0.23873/rhosol)

       if(iccnlev>=2 .and. itrkepsilon==1 .and. (acat==aerocat-1.or.acat==aerocat)) &
         regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) * (aeromas(k,acat) / aeromass)

     endif

  endif

enddo

return
END SUBROUTINE prenuc_ccn

!##############################################################################
Subroutine prenuc_ifn (m1,k,dn0,rv)

use micphys
use mem_grid

implicit none

integer :: m1,k,in_thresh
real :: tot_in,ifnfrac
real :: vapnuc,vapnucr,availvap,immersed
real, dimension(m1) :: dn0,rv

!**********************************************************************
!*********** ICE NUCLEI DETERMIINATION ********************************
!**********************************************************************
 !For determining number concentration of aerosols with
 !diameter greater than 0.5 microns for use as IN
 tot_in = 0.0
 total_in(k) = 0.0

 !Loop over CCN, GCCN, SMALL DUST, LARGE DUST, REGEN 1&2 (1,2,3,4,8,9)
 !Not looping over salt species since these cannot act is ice nuclei
 do acat=1,aerocat

   in_thresh = 1
   totifnn(k,acat) = 0.0
   totifnm(k,acat) = 0.0

   if((acat==1)                  .or. &  ! CCN
      (acat==2)                  .or. &  ! GCCN
      (acat==3 .and. idust>0)    .or. &  ! Small dust mode
      (acat==4 .and. idust>0)    .or. &  ! Large dust mode
      (acat==8 .and. iabcarb>0)  .or. &  ! Absorbing carbon 1 mode
      (acat==9 .and. iabcarb>0)  .or. &  ! Absorbing carbon 2 mode
      (acat==aerocat-1 .and. iccnlev>=2) .or. &  ! Small regenerated aerosol
      (acat==aerocat   .and. iccnlev>=2)) then   ! Large regenerated aerosol

    concen_nuc = aerocon(k,acat)
    aeromass   = aeromas(k,acat)
    rg         = aero_rg(acat)
    rhosol     = aero_rhosol(acat)

    !If aerosols are present in the given category then proceed
    if(concen_nuc > mincon) then

     !Set up binned distribution mass and sizes
     if(rg<=0.015e-6) then
       rmsma = 1.0e-23
       rmlar = 2.0e-16
     elseif(rg>0.015e-6 .and. rg<=0.03e-6) then
       rmsma = 1.0e-22
       rmlar = 2.0e-15
     elseif(rg>0.03e-6 .and. rg<=0.04e-6) then
       rmsma = 1.0e-21
       rmlar = 2.0e-14
     elseif(rg>0.04e-6 .and. rg<=0.08e-6) then
       rmsma = 3.0e-21
       rmlar = 6.0e-14
     elseif(rg>0.08e-6 .and. rg<=0.16e-6) then
       rmsma = 5.0e-20
       rmlar = 1.0e-12
     elseif(rg>0.16e-6 .and. rg<=0.32e-6) then
       rmsma = 5.0e-19
       rmlar = 1.0e-11
     elseif(rg>0.32e-6 .and. rg<=0.64e-6) then
       rmsma = 1.0e-18
       rmlar = 2.0e-11
     elseif(rg>0.64e-6 .and. rg<=0.96e-6) then
       rmsma = 1.0e-17
       rmlar = 2.0e-10
     elseif(rg>0.96e-6 .and. rg<=2.00e-6) then
       rmsma = 3.0e-17
       rmlar = 6.0e-10
     elseif(rg>2.00e-6 .and. rg<=3.00e-6) then
       rmsma = 1.0e-16
       rmlar = 2.0e-09
     elseif(rg>3.00e-6 .and. rg<=4.00e-6) then
       rmsma = 7.0e-16
       rmlar = 1.5e-08
     elseif(rg>4.00e-6 .and. rg<=5.00e-6) then
       rmsma = 7.0e-16
       rmlar = 1.5e-08
     elseif(rg>5.00e-6 .and. rg<=6.00e-6) then
       rmsma = 1.0e-15
       rmlar = 2.0e-08
     elseif(rg>6.00e-6) then
       rmsma = 1.0e-15
       rmlar = 2.0e-08
     endif

     !Convert units for setting up lognormal distribution
     concen_nuc = concen_nuc * dn0(k) !Convert #/kg to #/m3
     aeromass   = aeromass   * dn0(k) !Convert kg/kg to kg/m3

     !Set up binned distribution concentration (#/m3)
     power = alog10(rmsma/rmlar) / float(itbin-1)
     do ic=1,itbin
      !bin radius equals [(3/4) x (1/pi) x (1/rho) x mass of bin] ^ (1/3)
      smass(ic) = rmlar * 10.**(float(ic-1) * power) !solute masses (kg)
      binrad(ic)=(0.23873/rhosol*smass(ic))**(.33333) !radius (meters)
     enddo
     do ic=1,itbin-1
      rcm=0.5*(binrad(ic)+binrad(ic+1))
      ccncon(ic) = concen_nuc/(1.47336267*rcm)*exp(-(alog(rcm/rg))**2/0.6909863)
      ccncon(ic) = 1.00032 * ccncon(ic)*(binrad(ic)-binrad(ic+1))
      ccnmas(ic) = 0.920 * smass(ic) * ccncon(ic)
      if(ic > 1) then
        ccncon(ic) = ccncon(ic) + ccncon(ic-1)
        ccnmas(ic) = ccnmas(ic) + ccnmas(ic-1)
      endif
      if(rcm<0.25e-6 .or. ccncon(ic)>=concen_nuc .or. ccnmas(ic)>=aeromass) then
        in_thresh=ic-1
        go to 101
      endif
     enddo
101  continue
     totifnn(k,acat) = ccncon(in_thresh) / dn0(k) !store in #/kg
     totifnm(k,acat) = ccnmas(in_thresh) / dn0(k) !store in kg/kg
     total_in(k) = total_in(k) + ccncon(in_thresh) / dn0(k) !store #/kg
    endif !if concen_nuc > 0

   endif !if aerosol species turned on

 enddo !acat do loop

 !DeMott IN activation based on total number of all aerosol
 !greater than 0.5 micron diameter.
 nifn(k) = 0.0
 if(tairc(k) .lt. 0.0 .and. total_in(k) .gt. mincon) then
   !For DeMott Eqn, the aerosol number to put into the eqn needs to be
   !the remaining of unprocessed aerosols > 0.5 microns + the remaining
   !aerosols > 0.5 microns contained in droplets + 
   !the tracked number that have already been ice nucleated
   immersed = immerhx(k,1) + immerhx(k,8) + immerhx(k,2)
   tot_in  = total_in(k) + immersed + ifnnucx(k)
   !Convert tot_in from #/kg to #/cm3
   tot_in = tot_in * dn0(k) / 1.e6

   !Convert to STP
   !Nc(STP) = Nc * (101300(Pa) * T(K)) / (Pressure * 273.2(K))
   !P=rho*Rd*T
   !Nc(STP) = Nc * (101300(Pa)/273.2(K)) * [1/(dn0*287(J/kg/k))] 
   tot_in  = tot_in*1.29/dn0(k)

   !Input aerosols in #/cm3 and outputs #/L activated
   if(iifn_formula==1) then
    !Original Demott(2010) formula
    nifn(k) = 0.0000594 * (-tairc(k))**3.33 &
            * (tot_in)**(0.0264*(-tairc(k))+0.0033)
   elseif(iifn_formula==2) then
    !Modified Demott(2010) for dust-dominated cases
    !Paul suggested an additional factor of 3 multiplier
    nifn(k) = 3.0 * 0.0008 * 10 ** (-0.2*(tairc(k)+9.7)) * tot_in ** 1.25
   endif

   !Adjust units and such
   nifn(k) = nifn(k)/1.29*dn0(k) !Convert FROM STP 
   tot_in  = tot_in /1.29*dn0(k) !Convert FROM STP
   tot_in  = tot_in  / dn0(k) * 1.e6 !Convert #/cm3 to #/kg
   nifn(k) = nifn(k) / dn0(k) * 1.e3 !Convert #/L to #/kg
   nifn(k) = min(tot_in,nifn(k)) !IFN in #/kg
 
   !Limit activation to only the number greater than that already
   !activated for the given grid cell parcel.
   nifn(k) = nifn(k) - ifnnucx(k)
   tot_in  = tot_in  - ifnnucx(k)
   if(nifn(k) < 1.0e-6) nifn(k) = 0.0
   if(tot_in  < 1.0e-6) tot_in  = 0.0
   ifnfrac = max(0.0,min(1.0,nifn(k)/tot_in))
   if(nifn(k)==0.0 .or. tot_in==0.0) ifnfrac = 0.0
 else
   nifn(k) = 0.0
   ifnfrac = 0.0
 endif

 !Pre-determine amount of IN that produce pristine ice and
 !reduce number of IN if excess vapor does not support nucleation
 !of all the pre-determined number of activate
 if(nifn(k) > 0.0 .and. tot_in > 0.0) then
   vapnuc = max(0.,nifn(k))
   vapnucr = vapnuc * emb0(3)
   if (vapnucr .gt. 0.) then
      availvap = .5 * (rv(k) - rvisair(k))
      if (vapnucr .gt. availvap) then
         vapnucr = min(vapnucr, max(0.,availvap))
      endif
   endif
   vapnuc = vapnucr / emb0(3)
   nifn(k) = vapnuc !IFN in #/kg
   ifnfrac = max(0.,min(1.,nifn(k)/tot_in))
 endif

 !Compute fraction of each aerosol species to reserve for icenuc
 !Then recompute the remainder of the aerosols for cldnuc
 if(nifn(k) > 0.0)then
 do acat=1,aerocat
   totifnn(k,acat) = totifnn(k,acat) * ifnfrac
   totifnm(k,acat) = totifnm(k,acat) * ifnfrac
   if(iccnlev>=1 .and. ifnfrac>0.0) then
    if((acat==1)                  .or. &  ! CCN
       (acat==2)                  .or. &  ! GCCN
       (acat==3 .and. idust>0)    .or. &  ! Small dust mode
       (acat==4 .and. idust>0)    .or. &  ! Large dust mode
       (acat==8 .and. iabcarb>0)  .or. &  ! Absorbing carbon 1 mode
       (acat==9 .and. iabcarb>0)  .or. &  ! Absorbing carbon 2 mode
       (acat==aerocat-1 .and. iccnlev>=2) .or. &  ! Small regenerated aerosol
       (acat==aerocat   .and. iccnlev>=2)) then   ! Large regenerated aerosol
      !Assign aerosol specs to local arrays
      epsil      = aero_epsilon(acat)

      !Aerosol and solubility tracking
      if(iccnlev>=2 .and. itrkepsilon==1 .and. (acat==aerocat-1.or.acat==aerocat) &
       .and. aeromas(k,acat)>0.) then
         epsil = min(1.0,regenmas(k,acat-(aerocat-2))/aeromas(k,acat))
      endif

      aerocon(k,acat) = aerocon(k,acat) - totifnn(k,acat)
      aeromas(k,acat) = aeromas(k,acat) - totifnm(k,acat)
      total_in(k)     = total_in(k)     - totifnn(k,acat)
      ifnnucx(k)      = ifnnucx(k)      + totifnn(k,acat)

      !Aerosol and solubility tracking
      !Store any aerosol mass in cnmhx arrays
      if(iccnlev>=2) then
       cnmhx(k,3) = cnmhx(k,3) + totifnm(k,acat)
       if(itrkepsilon==1) then
         snmhx(k,3) = snmhx(k,3) + totifnm(k,acat) * epsil
         if(acat==aerocat-1.or.acat==aerocat) &
          regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) - totifnm(k,acat) * epsil         
       endif
       if(itrkdust==1 .and. (acat==3 .or. acat==4)) &
         dnmhx(k,3) = dnmhx(k,3) + totifnm(k,acat)
       if(itrkdustifn==1 .and. (acat==3 .or. acat==4)) &
         dinhx(k,3) = dinhx(k,3) + totifnm(k,acat)
      endif

      aero_rg(acat) = aero_medrad(acat) ! Default median radius
      if(aerocon(k,acat) > mincon) then
       rhosol=aero_rhosol(acat)
       aero_rg(acat)=((0.23873/rhosol*aeromas(k,acat)/aerocon(k,acat)) &
             **(1./3.))/aero_rg2rm(acat)
       if(aero_rg(acat) < 0.01e-6) aero_rg(acat) = 0.01e-6
       if(aero_rg(acat) > 6.50e-6) aero_rg(acat) = 6.50e-6
       aeromas(k,acat) = ((aero_rg(acat)*aero_rg2rm(acat))**3.) &
             *aerocon(k,acat)/(0.23873/rhosol)
      endif
    endif
   endif
 enddo
 endif

return
END SUBROUTINE prenuc_ifn

