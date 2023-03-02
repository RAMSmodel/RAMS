!##############################################################################
Subroutine cldnuc (m1,k1cnuc,k2cnuc,k1dnuc,k2dnuc,rv,wp,i,j,dn0)

use micphys
use mem_grid

implicit none

integer :: m1,i,j,k,k1cnuc,k2cnuc,k1dnuc,k2dnuc &
          ,jtemp,jw,jconcen,rgccn1,ctc,epsnum,epstab,drop
real :: rnuc,excessrv,rcnew,vaprccn,num_ccn_ifn,epstemp
real, dimension(m1) :: rv,wp,dn0
real :: tairc_nuc,w_nuc,rg_nuc,tab,sfcareatotal
real :: rjw,wtw1,wtw2,rjconcen,wtcon1,wtcon2,jrg1,jrg2,eps1,eps2
real :: total_cld_nucc,total_drz_nucc,total_cld_nucr,total_drz_nucr
real, dimension(9) :: concen_tab

!Re-set cloud layer before nucleation
k1cnuc = 2
k2cnuc = 1
!Re-set drizzle layer before nucleation if number prognostic
if(jnmb(8)>=5)then
 k1dnuc = 2
 k2dnuc = 1
else
 k1dnuc = m1
 k2dnuc = 1
endif

!*********************************************************
!***** If NOT pronosing number concentration of cloud1****
!*********************************************************
if (jnmb(1) == 1 .or. jnmb(1) == 4) then
   rnuc = parm(1) * emb0(1)
   do k = 2,m1-1
      excessrv = rv(k) - 1.0001 * rvlsair(k)
      rcnew = 0.
      if (excessrv > 0.) then
         rcnew = min(rnuc,.5*excessrv)
         rx(k,1) = rx(k,1) + rcnew
         rv(k) = rv(k) - rcnew
         k2cnuc = k
         cx(k,1) = min(parm(1),rx(k,1) / emb0(1))
         if(imbudget >= 1) then
           xnuccldrt(k) = xnuccldrt(k) + rcnew * budget_scalet
         endif
      elseif (k2cnuc == 1) then
         k1cnuc = k + 1
      endif
   enddo

!*************************************************************************
!Saleeby(6/3/02) Prognosing number concentration of cloud and drizzle ****
!*************************************************************************
elseif (jnmb(1) >= 5) then
 do k = 2,m1-1

  !Do tolerance check for all aerosol types
  CALL prenuc_ccn (k,i,j)

  !Set temp variable to zero to carry potential IN concentration
  !as calculated in prenuc_ifn from aerosol distributions.
  total_in(k) = 0.0

  excessrv = rv(k) - 1.0001 * rvlsair(k)
  if (excessrv > 0.) then

   !Set aside ice nuclei for DeMott scheme and recomute aerosol stats 
   if(iifn==3) CALL prenuc_ifn (m1,k,dn0,rv)

!**********LOOP OVER CCN, GCCN, 2 DUST MODES, 3 SALT MODES ***************
   !Use the acat loop for turning on aerosol nucleation
   do acat=1,aerocat

      concen_tab(acat) = 0.0
      concen_nuc = 0.0

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
       concen_nuc = aerocon(k,acat)
       aeromass   = aeromas(k,acat)
       rg         = aero_rg(acat)
       epsil      = aero_epsilon(acat)

       !Aerosol and solubility tracking
       if(iccnlev>=2 .and. itrkepsilon==1 .and. (acat==aerocat-1.or.acat==aerocat) &
        .and. aeromas(k,acat)>0.) &
          epsil = min(1.0,regenmas(k,acat-(aerocat-2))/aeromas(k,acat))

       !Temporary check to make sure regensol is not > regen        
       if(acat==aerocat-1.or.acat==aerocat)then
        if(regenmas(k,acat-(aerocat-2)) > aeromas(k,acat))then
          print*,'Soluble regenerated mass > Regenerated mass' &
          ,acat,k,i,j,regenmas(k,acat-(aerocat-2))/aeromas(k,acat) &
          ,regenmas(k,acat-(aerocat-2)),aeromas(k,acat)
        endif
       endif

       !If not removing aerosols, subtract off Ice nuclei amount
       if(iifn==3 .and. iccnlev==0)then
        concen_nuc = concen_nuc - totifnn(k,acat)
       endif

       !Only nucleate if aerosols are populous
       if(concen_nuc > mincon .and. aeromass > minmas) then

        !**************TEMPERATURE DEPENDENCY*********************************
        tairc_nuc = tairc(k)
        if (tairc_nuc < -30.) then
           tairc_nuc = -30.
        elseif (tairc_nuc > 30.) then
           tairc_nuc = 30.
        endif
        jtemp = nint(.1 * (tairc_nuc + 30.)) + 1
        !*******VERTICAL VELOCITY DEPENDENCY**********************************
        w_nuc = wp(k)
        if (w_nuc < .010001) then 
           w_nuc = .010001
        elseif (w_nuc > 99.99) then
           w_nuc = 99.99
        endif
        rjw = 2. * log10(100. * w_nuc) + 1.
        jw = int(rjw)
        wtw2 = rjw - float(jw)
        wtw1 = 1. - wtw2
        !********** AEROSOL NUMBER, MASS, & MEDIAN RADIUS CONSTRAINTS ********
        rjconcen = max(1., min(7., 2. * log10(1.0e-7 * concen_nuc) + 1.))
        jconcen = int(rjconcen)
        wtcon2 = rjconcen - float(jconcen)
        wtcon1 = 1. - wtcon2
        !********** MEDIAN RADIUS DEPENDENCY FOR CCN *************************
        rg_nuc = rg
        if (rg_nuc < 0.01e-6) then
           rg_nuc = 0.011e-6
        elseif (rg_nuc > 0.96e-6) then
           rg_nuc = 0.959e-6
        endif
        do rgb=1,maxrg-1
         if((rg_nuc>=rg_ccn(rgb)) .and. (rg_nuc<=rg_ccn(rgb+1))) then
           rgccn1=rgb
           jrg2 = (rg_nuc-rg_ccn(rgb)) / (rg_ccn(rgb+1)-rg_ccn(rgb))
           jrg1 = 1. - jrg2
         endif
        enddo
        !********** EPSILON SOLUBILITY FRACTION FOR CCN **********************
        !Determine weights for interpolating between epsilon table values
        epstemp = epsil
        if (epsil < 0.05) then
           epsil = 0.05
        elseif (epsil > 1.00) then
           epsil = 1.00
        endif
        do epsnum=1,maxeps-1
         if((epsil>=epsfrac(epsnum)) .and. (epsil<=epsfrac(epsnum+1))) then
          epstab=epsnum
          eps2 = (epsil-epsfrac(epsnum)) / (epsfrac(epsnum+1)-epsfrac(epsnum))
          eps1 = 1. - eps2
         endif
        enddo
        !Reduce nucleation fraction for very small solubility
        epstemp = max(0.0,min(1.0,log10(epstemp*100. + 1.0)))
        !******************* DETERMINE LOOKUP TABLE VALUES *******************
        tab=0.0
        if(iaero_chem(acat)==1) then
          CALL aero_nuc_tab_nh42so4 (rv(k),rvlsair(k),eps1,eps2,wtw1,wtw2,wtcon1 &
            ,wtcon2,jrg1,jrg2,epstab,jw,jconcen,jtemp,rgccn1,tab)
        elseif(iaero_chem(acat)==2) then
          CALL aero_nuc_tab_nacl (rv(k),rvlsair(k),eps1,eps2,wtw1,wtw2,wtcon1 &
            ,wtcon2,jrg1,jrg2,epstab,jw,jconcen,jtemp,rgccn1,tab)
        endif
        concen_tab(acat) = concen_nuc * tab * epstemp

       endif !if concen_nuc > mincon
      endif !if acat aerosol type is valid
   enddo !looping over aerocat

!****Compute total surface areas and particle percentage of total*************
   sfcareatotal=0.0
   do acat=1,aerocat
    aero_ratio(acat) = 0.0  ! Aerosol fraction
    aero_vap(acat)   = 0.0  ! Total surface area of aerosol category
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
       if(aerocon(k,acat) > mincon) then
        aero_vap(acat) = (4.0 * 3.14159 * aero_rg(acat)**2) * concen_tab(acat)
        sfcareatotal = sfcareatotal + aero_vap(acat)
       endif
    endif
   enddo
   if(sfcareatotal > 0.0) then
    do acat=1,aerocat
      aero_ratio(acat) = aero_vap(acat) / sfcareatotal
    enddo
   endif

!*********FOR NUMBER CONCENTRATION PREDICTION OF CLOUD ***********************
   total_cld_nucc=0.0
   total_drz_nucc=0.0
   total_cld_nucr=0.0
   total_drz_nucr=0.0
   do acat=1,aerocat
     ctc = 0
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
      concen_nuc = aerocon(k,acat)
      aeromass   = aeromas(k,acat)
      rg         = aero_rg(acat)
      rhosol     = aero_rhosol(acat)
      cldrat     = aero_ratio(acat)
      epsil      = aero_epsilon(acat)

      !Aerosol and solubility tracking
      if(iccnlev>=2 .and. itrkepsilon==1 .and. (acat==aerocat-1.or.acat==aerocat) &
        .and. aeromas(k,acat)>0.) &
          epsil = min(1.0,regenmas(k,acat-(aerocat-2))/aeromas(k,acat))

      !If not removing aerosols, subtract off Ice nuclei amount
      if(iifn==3 .and. iccnlev==0)then
       concen_nuc = concen_nuc - totifnn(k,acat)
       aeromass   = aeromass   - totifnm(k,acat)
      endif

      if(concen_tab(acat) > mincon) ctc=1 !If #/kg > mincon, nucleate

      if(ctc==1) then

        !Vapor allocated to a given aerosol species
        vaprccn = 0.5*excessrv*cldrat !Sum of all nucleation <= 1/2 excessrv

        !Determine if nucleated droplets go to cloud or drizzle
        drop=1
        if(rg <= 0.96e-6) drop=1 !Nucleate to small droplet mode
        if(rg >  0.96e-6 .and. jnmb(8)>=5) drop=8 !Nucleate to drizzle mode

        !Keep nucleated droplets in size bounds 
        if(concen_tab(acat) > vaprccn / emb0(drop)) concen_tab(acat) = vaprccn / emb0(drop)
        if(concen_tab(acat) < vaprccn / emb1(drop)) vaprccn = concen_tab(acat) * emb1(drop)

!Nucleate to 2-micron diameter droplets whose mass is emb0
!if(concen_tab(acat) * emb0(drop) <= vaprccn) vaprccn=concen_tab(acat) * emb0(drop)

        !Accumulated nucleated particles if not removing them
        if(iccnlev==0) then
         if(drop==1)total_cld_nucc = total_cld_nucc + concen_tab(acat)
         if(drop==1)total_cld_nucr = total_cld_nucr + vaprccn
         if(drop==8)total_drz_nucc = total_drz_nucc + concen_tab(acat)
         if(drop==8)total_drz_nucr = total_drz_nucr + vaprccn
        endif

        !Add nucleated cloud water and number here if removing aerosol   
        if(iccnlev>=1)then
          cx(k,drop) = cx(k,drop) + concen_tab(acat)
          rx(k,drop) = rx(k,drop) + vaprccn

          !Nucleation budget diagnostics
          if(imbudget >= 1) then
            xnuccldrt(k) = xnuccldrt(k) + vaprccn * budget_scalet
          endif
          !Dust Budget diagnostics
          if(acat==3 .and. drop==1 .and. imbudget==3 .and. idust >= 1) &
             xdust1cldrt(k) = xdust1cldrt(k) + vaprccn * budget_scalet
          if(acat==3 .and. drop==8 .and. imbudget==3 .and. idust >= 1) &
             xdust1drzrt(k) = xdust1drzrt(k) + vaprccn * budget_scalet
          if(acat==4 .and. drop==1 .and. imbudget==3 .and. idust >= 1) &
             xdust2cldrt(k) = xdust2cldrt(k) + vaprccn * budget_scalet
          if(acat==4 .and. drop==8 .and. imbudget==3 .and. idust >= 1) &
             xdust2drzrt(k) = xdust2drzrt(k) + vaprccn * budget_scalet

          !Convert units for setting up lognormal distribution
          concen_tab(acat) = concen_tab(acat) * dn0(k) !Convert #/kg to #/m3
          concen_nuc       = concen_nuc       * dn0(k) !Convert #/kg to #/m3
          aeromass         = aeromass         * dn0(k) !Convert kg/kg to kg/m3

          !For determining mass of distirubtion to remove, use "concen_nuc"
          !Doing size specific removal. Preferentially remove larger particles.
          !Set up binned distribution masses(kg) and sizes(meters)
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

          !Set up binned distribution concentration
          power = alog10(rmsma/rmlar) / float(itbin-1)
          do ic=1,itbin
           !bin radius equals [(3/4) x (1/pi) x (1/rho) x mass of bin] ^ (1/3)
           smass(ic) = rmlar * 10.**(float(ic-1) * power) !solute masses (kg)
           binrad(ic)=(0.23873/rhosol*smass(ic))**(.33333) !radius (meters)
          enddo

          !Loop thru bins to determine amount of mass to remove based on number.
          !This removes from the large end of the distribution first and works toward
          !the smaller end with the assumption that larger particles activate first
          ccnmass = 0.0 !Variable ccnmass is kg/m3
          num_ccn_ifn = 0.0 !Number of particles > 0.25 microns radius
          do ic=1,itbin-1
           rcm=0.5*(binrad(ic)+binrad(ic+1))
           ccncon(ic) = concen_nuc/(1.47336267*rcm)*exp(-(alog(rcm/rg))**2/0.6909863)
           ccncon(ic) = 1.00032 * ccncon(ic)*(binrad(ic)-binrad(ic+1))
           ccnmas(ic) = 0.920 * smass(ic) * ccncon(ic)
           if(ic>1)then
            ccncon(ic) = ccncon(ic) + ccncon(ic-1)
            ccnmas(ic) = ccnmas(ic) + ccnmas(ic-1)
           endif
           !Track immersion freezing droplets that contain large CCN, GCCN, or DUST
           ! Do not track immersion freezing for salt species (acat=5,6,7)
           if(iifn==3.and.(acat==1.or.acat==2.or.acat==3.or.acat==4.or.acat==8.or.acat==9 &
               .or.acat==aerocat-1.or.acat==aerocat) &
               .and. rcm > 0.25e-6 .and. ic>1) num_ccn_ifn=ccncon(ic-1)
           !Track the amount of aerosol mass contained within new droplets 
           if(ccncon(ic)>=concen_tab(acat) .or. ccnmas(ic)>=aeromass .or. ic==itbin-1) then
             !Further immersion freezing tracking for (acat=1,2,3,4,8,9)
             if(iifn==3.and.(acat==1.or.acat==2.or.acat==3.or.acat==4.or.acat==8.or.acat==9 &
               .or.acat==aerocat-1.or.acat==aerocat) &
               .and. rcm > 0.25e-6 .and. ic>1) num_ccn_ifn=concen_tab(acat)
             ccnmass=ccnmas(ic-1)
             go to 111
           endif
          enddo
111       continue

          !If either is zero, set both to zero
          if(concen_tab(acat)==0.0 .or. ccnmass==0.0) then
             concen_tab(acat)=0.0
             ccnmass=0.0
          endif

          !Convert #/m3 back to #/kg and kg/m3 back to kg/kg
          concen_tab(acat) = concen_tab(acat) / dn0(k) !Convert #/m3 to #/kg
          ccnmass          = ccnmass          / dn0(k) !Convert kg/m3 to kg/kg
          num_ccn_ifn      = num_ccn_ifn      / dn0(k) !Convert #/m3 to #/kg

          !Subtract off aerosol mass and number and keep both + or zero
          aerocon(k,acat) = aerocon(k,acat) - concen_tab(acat)
          aeromas(k,acat) = aeromas(k,acat) - ccnmass
          if(aerocon(k,acat)<=0.0 .or. aeromas(k,acat)<=0.0) then
            aerocon(k,acat) = 0.0
            aeromas(k,acat) = 0.0
          endif

          !Aerosol and solubility tracking
          if(iccnlev>=2)then
            cnmhx(k,drop) = cnmhx(k,drop) + ccnmass
            if(itrkepsilon==1) then
             snmhx(k,drop) = snmhx(k,drop) + ccnmass * epsil
             if(acat==aerocat-1.or.acat==aerocat) &
              regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) - ccnmass * epsil
            endif
            if(itrkdust==1 .and. (acat==3 .or. acat==4)) &
              dnmhx(k,drop) = dnmhx(k,drop) + ccnmass
          endif

          !Store number of large particles for immersion freezing
          !if (1) 2-moment cloud, (2) iccnlev>=1, (3) iifn==3
          if(iifn==3) then
            immerhx(k,drop) = immerhx(k,drop) + num_ccn_ifn !#/kg
            total_in(k) = total_in(k) - num_ccn_ifn
          endif

        endif !if ICCNLEV>=1 then remove aerosol mass
      endif !if ctc==1
     endif !if acat aerosol type is valid
   enddo !loop over acat 1 to aerocat

   !If not removing aerosol, only add number in excess of droplet number
   if(iccnlev==0) then
     total_cld_nucc = max(0.,total_cld_nucc - cx(k,1))
     total_drz_nucc = max(0.,total_drz_nucc - cx(k,8))
     if(total_cld_nucc > total_cld_nucr / emb0(1)) &
        total_cld_nucc = total_cld_nucr / emb0(1)
     if(total_cld_nucc < total_cld_nucr / emb1(1)) &
        total_cld_nucr = total_cld_nucc * emb1(1)
     cx(k,1) = cx(k,1) + total_cld_nucc
     rx(k,1) = rx(k,1) + total_cld_nucr
     if(total_drz_nucc > total_drz_nucr / emb0(8)) &
        total_drz_nucc = total_drz_nucr / emb0(8)
     if(total_drz_nucc < total_drz_nucr / emb1(8)) &
        total_drz_nucr = total_drz_nucc * emb1(8)
     cx(k,8) = cx(k,8) + total_drz_nucc
     rx(k,8) = rx(k,8) + total_drz_nucr
     !Nucleation budget diagnostics
     if(imbudget >= 1) then
       xnuccldrt(k) = xnuccldrt(k) &
             + (total_cld_nucr + total_drz_nucr) * budget_scalet
     endif
   endif

  endif !if excess vapor

  !Update cloud layer following nucleation
  if (rx(k,1) .ge. rxmin) k2cnuc = k
  if (k2cnuc .eq. 1 .and. rx(k,1) .lt. rxmin) k1cnuc = k + 1
  !Update drizzle layer following nucleation if number prognostic
  if(jnmb(8)>=5) then
   if (rx(k,8) .ge. rxmin) k2dnuc = k
   if (k2dnuc .eq. 1 .and. rx(k,8) .lt. rxmin) k1dnuc = k + 1
  endif

 enddo !loop over all vertical levels

else
 print*, 'icloud not allowed to be 2 or 3'
 print*, 'stopping model '
 stop 'icloud'
endif !if number prediction of cloud droplets

return
END SUBROUTINE cldnuc

!##############################################################################
Subroutine icenuc (m1,kc1,kc2,kd1,kd2,k1pnuc,k2pnuc,ngr,rv,dn0,dtlt,i,j)

use rconstants
use micphys

implicit none

integer :: m1,kc1,kc2,kd1,kd2,k1pnuc,k2pnuc,ngr,k,idnc,itc,irhhz,ithz,lcat,i,j
real :: dn1,fraccld,ridnc,dtlt,ssi0,wdnc2,tc,ritc,wtc2  &
       ,pbvi,ptvi,pdvi,ptotvi,fracifn,cldnuc,rhhz,haznuc  &
       ,rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,ssi,diagni,heterofrac  &
       ,vapnuc,vapnucr,availvap,cont_nuc,homo_nuc,nucfrac,pcthaze,excessrv

real :: concen_tab,embtemp,frzc,frzr,tot_in,ifntemp,immerin,immersed

real, dimension(m1) :: rv,dn0

! Define ssi0 to be maximum supersaturation with respect to ice for
! determining total number of IFN that can nucleate in Meyers' formula
data ssi0/0.40/
save

! implement paul's immersion freezing of rain here.  This would
! replace mike's homogeneous freezing of rain which was in h03.

!************************************************************************
!************* CLOUD DROPLET HOMOGENEOUS ICE NUCLEATION******************
!************************************************************************
do k = kc1,kc2

 !If cloud water exists at a minimum quantity
 if (rx(k,1) .gt. rxmin) then

   !define dn locally from emb
   dn1 = dnfac(1) * emb(k,1) ** pwmasi(1)
   fraccld=0.
   if (tairc(k) .le. -30.01) then
      ridnc = max(1.,min(float(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)
      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)
      if(fraccld > 0.990) fraccld=1.0
   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
   CALL contnuc (rx(k,1),cx(k,1),tx(k,1),vap(k,1),press(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt)

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction
!::later   ptotvi = ptotvi * fracifn
! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.
   !Saleeby(2009): Need separate homogeneous freezing options for
   !1-moment and 2-moment cloud and drizzle droplet treatments
   if(jnmb(1) <  5) cldnuc = max(0.,fraccld * cx(k,1) - cx(k,3))
   if(jnmb(1) >= 5) cldnuc = max(0.,fraccld * cx(k,1))

   cont_nuc = ptotvi * emb(k,1)
   homo_nuc = fraccld * rx(k,1)

   if(cont_nuc + homo_nuc > rx(k,1)) then
     nucfrac=rx(k,1)/(cont_nuc + homo_nuc)
     cont_nuc=cont_nuc*nucfrac
     homo_nuc=homo_nuc*nucfrac
     cldnuc=cldnuc*nucfrac
     ptotvi=ptotvi*nucfrac
   endif

   !Aerosol and solubility tracking
   !Transfering aerosol mass from cloud to pristine ice
   if(iccnlev>=2) then
    rxferratio = min(1.0, (cont_nuc + homo_nuc) / rx(k,1))
    ccnmass  = cnmhx(k,1) * rxferratio
    cnmhx(k,1) = cnmhx(k,1) - ccnmass
    cnmhx(k,3) = cnmhx(k,3) + ccnmass
    if(itrkepsilon==1)then
     scnmass  = snmhx(k,1) * rxferratio
     snmhx(k,1) = snmhx(k,1) - scnmass
     snmhx(k,3) = snmhx(k,3) + scnmass
    endif
    if(itrkdust==1)then
     dcnmass  = dnmhx(k,1) * rxferratio
     dnmhx(k,1) = dnmhx(k,1) - dcnmass
     dnmhx(k,3) = dnmhx(k,3) + dcnmass
    endif
    if(itrkdustifn==1)then
     dinmass  = dinhx(k,1) * rxferratio
     dinhx(k,1) = dinhx(k,1) - dinmass
     dinhx(k,3) = dinhx(k,3) + dinmass
    endif
   endif

   !Remove immersion freezing nuclei from cloud 
   if(iifn==3 .and. iccnlev>=1 .and. cx(k,1)>0.0) then
    enxferratio = min(1.0,max(0.0,(cldnuc+ptotvi)/cx(k,1)))
    ccnnum = immerhx(k,1) * enxferratio
    immerhx(k,1) = immerhx(k,1) - ccnnum
   endif

   rx(k,3) = rx(k,3) + min(rx(k,1),cont_nuc + homo_nuc)
   rx(k,1) = rx(k,1) - min(rx(k,1),cont_nuc + homo_nuc)
   cx(k,3) = cx(k,3) + min(cx(k,1),cldnuc + ptotvi)
   cx(k,1) = cx(k,1) - min(cx(k,1),cldnuc + ptotvi)

   if(imbudget >= 1) then
     xnucicert(k) = xnucicert(k) + (cont_nuc + homo_nuc) * budget_scalet
   endif
   if(imbudget >= 2) then
     xinuchomrt(k)  = xinuchomrt(k)  + homo_nuc * budget_scalet
     xinuccontrt(k) = xinuccontrt(k) + cont_nuc * budget_scalet
   endif

 endif
enddo

!************************************************************************
!************* DRIZZLE DROPLET HOMOGENEOUS ICE NUCLEATION****************
!************************************************************************
do k = kd1,kd2

 !If drizzle water exists at a minimum quantity
 if (rx(k,8) .gt. rxmin) then

   !define dn locally from emb
   dn1 = dnfac(16) * emb(k,8) ** pwmasi(16)
   fraccld = 0.
   if (tairc(k) .le. -30.01) then
      ridnc = max(1.,min(float(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)
      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)
      if(fraccld > 0.990) fraccld=1.0
   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
   CALL contnuc (rx(k,8),cx(k,8),tx(k,8),vap(k,8),press(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt)

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction
!::later   ptotvi = ptotvi * fracifn
! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.
   !Saleeby(2009): Need separate homogeneous freezing options for
   !1-moment and 2-moment cloud and drizzle droplet treatments
   if(jnmb(8) <  5) cldnuc = max(0.,fraccld * cx(k,8) - cx(k,3))
   if(jnmb(8) >= 5) cldnuc = max(0.,fraccld * cx(k,8))

   cont_nuc = ptotvi * emb(k,8)
   homo_nuc = fraccld * rx(k,8)

   if(cont_nuc + homo_nuc > rx(k,8)) then
     nucfrac=rx(k,8)/(cont_nuc + homo_nuc)
     cont_nuc=cont_nuc*nucfrac
     homo_nuc=homo_nuc*nucfrac
     cldnuc=cldnuc*nucfrac
     ptotvi=ptotvi*nucfrac
   endif

   !Aerosol and solubility tracking
   !Transfering aerosol mass from drizzle to pristine ice
   if(iccnlev>=2) then
    rxferratio = min(1.0,(cont_nuc + homo_nuc) / rx(k,8))
    ccnmass  = cnmhx(k,8) * rxferratio
    cnmhx(k,8) = cnmhx(k,8) - ccnmass
    cnmhx(k,3) = cnmhx(k,3) + ccnmass
    if(itrkepsilon==1)then
     scnmass  = snmhx(k,8) * rxferratio
     snmhx(k,8) = snmhx(k,8) - scnmass
     snmhx(k,3) = snmhx(k,3) + scnmass
    endif
    if(itrkdust==1)then
     dcnmass  = dnmhx(k,8) * rxferratio
     dnmhx(k,8) = dnmhx(k,8) - dcnmass
     dnmhx(k,3) = dnmhx(k,3) + dcnmass
    endif
    if(itrkdustifn==1)then
     dinmass  = dinhx(k,8) * rxferratio
     dinhx(k,8) = dinhx(k,8) - dinmass
     dinhx(k,3) = dinhx(k,3) + dinmass
    endif
   endif

   !Remove immersion freezing nuclei from cloud 
   if(iifn==3 .and. iccnlev>=1 .and. cx(k,8)>0.0) then
    enxferratio = min(1.0,max(0.0,(cldnuc+ptotvi)/cx(k,8)))
    ccnnum = immerhx(k,8) * enxferratio
    immerhx(k,8) = immerhx(k,8) - ccnnum
   endif

   rx(k,3) = rx(k,3) + min(rx(k,8),cont_nuc + homo_nuc)
   rx(k,8) = rx(k,8) - min(rx(k,8),cont_nuc + homo_nuc)
   cx(k,3) = cx(k,3) + min(cx(k,8),cldnuc + ptotvi)
   cx(k,8) = cx(k,8) - min(cx(k,8),cldnuc + ptotvi)

   if(imbudget >= 1) then
     xnucicert(k) = xnucicert(k) + (cont_nuc + homo_nuc) * budget_scalet
   endif
   if(imbudget >= 2) then
     xinuchomrt(k)  = xinuchomrt(k)  + homo_nuc * budget_scalet
     xinuccontrt(k) = xinuccontrt(k) + cont_nuc * budget_scalet
   endif

 endif
enddo

!************************************************************************
!  Homogeneous nucleation of haze
!************************************************************************
k1pnuc = 2
k2pnuc = 1
do k = 2,m1-1
   rhhz = rv(k) / rvlsair(k)
   haznuc = 0.
   if (rhhz .gt. 0.82 .and. tairc(k) .le. -35.01) then
      rirhhz = min(0.1799,rhhz-0.82) / drhhz + 1.0
      irhhz = int(rirhhz)
      wrhhz2 = rirhhz - float(irhhz)
      thz = max(-59.99,tairc(k))
      rithz = (thz + 60.00) / dthz + 1.0
      ithz = int(rithz)
      wthz2 = rithz - float(ithz)
      frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  )  &
              +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  )  &
              + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1)  &
              +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
      frachaz = 1. - exp(-frachaz * dtlt)

      !Saleeby(2009): Haze nuclei can be too plentiful here compared
      ! to reality. For 2-moment cloud droplet prediction I scale the
      ! haze nuclei to the CCN concentration. Need better option here.
      ! Need haznuc in #/kg
      if(jnmb(1)>=5) haznuc = frachaz * aerocon(k,1)
      if(jnmb(1)< 5) haznuc = frachaz * 300.e6
   endif

!************************************************************************
! Heterogeneous nucleation by immersion deposition condensation freezing
!************************************************************************
   !use Meyers formula CIFNX(#/kg)
   diagni=0.0
   fracifn=0.0
   if(iifn==1) then
     ssi = min(ssi0,rv(k) / rvisair(k) - 1.)
     if (ssi .gt. 0. .and. tairc(k) .le. -5.) then
       !Meyers formula
       fracifn = exp(12.96 * (ssi - ssi0))
       !DeMott SPL modification
       !fracifn = (10 ** (-4.077421 * ssi + 0.097562)) * fracifn
       !Turn off traditional IN heterogeneous nucleation
       diagni = fracifn * cifnx(k)
     else
       diagni = 0.0
     endif
   !use DeMott(2010) IN nucleation NIFN(#/kg) from PPARM
   elseif(iifn==2) then
     excessrv = rv(k) - 1.0001 * rvlsair(k)
     if (excessrv > 0. .and. tairc(k) < 0.0) then
       !Convert tot_in from #/kg to #/cm3
       tot_in = cifnx(k) * dn0(k) / 1.e6
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
       tot_in  = tot_in * 1000. !Convert to #/liter
       nifn(k) = min(tot_in,nifn(k)) !IFN in #/liter
       diagni = nifn(k) * 1.e3 / dn0(k) !Convert to #/kg
     else
       diagni = 0.0
     endif
   !use DeMott(2010) IN nucleation NIFN(#/kg) from aerosol field
   elseif(iifn==3) then
     excessrv = rv(k) - 1.0001 * rvlsair(k)
     if (excessrv > 0.) then
       diagni = nifn(k)
     else
       diagni = 0.0
     endif
   endif

   !Orig Meyers formula: + diagni = exp(6.269 + 12.96 * ssi)
   !Combine nucleation types, and limit amounts
   !vapnuc is #/kg_air and vapnucr is kg/kg_air
   !BEGIN MIKE'S SECTION FOR LIMITING NUMBER OF CRYSTALS NUCLEATED
   !BY NUMBER OF ICE CRYSTALS PRESENT ALREADY FOR 1-MOMENT MICRO
   vapnuc=0.0
   if(jnmb(1)>=5)then
     if(iccnlev==0) then
       vapnuc = max(0.,haznuc + diagni - cx(k,3))
     elseif(iccnlev>=1 .and. iifn==3) then 
       vapnuc = max(0.,haznuc + diagni)
     elseif(iccnlev>=1 .and. iifn<=2) then 
       vapnuc = max(0.,haznuc)
       vapnuc = vapnuc + max(0.,diagni - cx(k,3))
     endif
   else
     vapnuc = max(0.,haznuc + diagni - cx(k,3))
   endif

   !Modify haze nucleation and IN nucleation if not enough vapor available
   vapnucr = vapnuc * emb0(3)
   if (vapnucr .gt. 0.) then
      availvap = .5 * (rv(k) - rvisair(k))
      if (vapnucr .gt. availvap) then
         vapnucr = min(vapnucr, max(0.,availvap))
      endif
   endif
   vapnuc = vapnucr / emb0(3)

   !(Saleeby10-17-2011) Determine haze nucleation and IN nucleation
   !DeMott(2010) formula already applied for ICLOUD >= 5
   !Allow large DeMott (diagni) particles to preferentially nucleate
   !over the haze particles (cccnp) for IIFN==3. Perhaps should include
   !all aerosols into potential haze nucleation?
   if((haznuc.gt.0.0 .or. diagni.gt.0.0) .and. jnmb(3).ge.5) then

    !Compute and test Haze fraction for DeMott(2010) aerosol micro
    if(iifn==3 .and. iccnlev>=1 .and. jnmb(1)>=5) then 
      heterofrac=(vapnuc-diagni)/max(1.0e-12,haznuc)
      heterofrac = min(1.0,max(0.0,heterofrac))
    !Compute and test Haze fraction for non-aerosol micro
    else
      heterofrac=(vapnuc)/(haznuc+diagni)
      heterofrac = min(1.0,max(0.0,heterofrac))
      diagni = diagni*heterofrac
    endif
    !Determine Haze fraction
    haznuc = haznuc*heterofrac

    !Subtract Haze particles from CCN field (acat==1) and add to tracking mass
    if(haznuc>0.0 .and. iccnlev>=1) then

     concen_nuc = aerocon(k,1)
     aeromass   = aeromas(k,1) 
     concen_tab = min(concen_nuc,haznuc)
     ccnmass = aeromass * (concen_tab/concen_nuc)
     concen_tab=max(0.0,min(concen_tab,aerocon(k,1)))
     ccnmass=max(0.0,min(ccnmass,aeromas(k,1)))
     aerocon(k,1) = aerocon(k,1) - concen_tab
     aeromas(k,1) = aeromas(k,1) - ccnmass

     !Aerosol and solubility tracking
     if(iccnlev>=2) then
      cnmhx(k,3) = cnmhx(k,3) + ccnmass
      if(itrkepsilon==1) snmhx(k,3) = snmhx(k,3) + ccnmass * aero_epsilon(1)
     endif

    endif

   endif

   !***********************************************************************
   !********** IMMERSION FREEZING OF DROPLETS *****************************
   !***********************************************************************
   !(Saleeby10-17-2011) Run DeMott formula again for immersed-in-cloud IN
   !DeMott IN activation based on total number of all aerosol
   !greater than 0.5 micron diameter. This is only done if:
   !(1) 2-moment cloud, (2) iccnlev>=1, (3) iifn==3
   immerin = 0.0
   if(iifn==3 .and. iccnlev>=1 .and. jnmb(1)>=5) then

    do lcat=1,ncat
     if(lcat==1.or.lcat==8.or.lcat==2) then

      frzc  = 0.0
      frzr  = 0.0
      ifntemp = 0.0

      !Check immersion nuclei following homogeneous freezing
      if(immerhx(k,lcat) > cx(k,lcat)) immerhx(k,lcat) = 0.9999 * cx(k,lcat)

      !If conditions allow for droplet immersion freezing
      if(tairc(k) < 0.0 .and. immerhx(k,lcat) > mincon) then
       !For DeMott Eqn, the aerosol number to put into the eqn needs to be
       !the remaining of unprocessed aerosols > 0.5 microns + the remaining
       !aerosols > 0.5 microns contained in droplets +
       !the tracked number that have already been ice nucleated
       immersed = immerhx(k,1) + immerhx(k,8) + immerhx(k,2)
       tot_in  = total_in(k) + immersed + ifnnucx(k)
       tot_in  = tot_in * dn0(k) / 1.e6 !Convert tot_in from #/kg to #/cm3
       tot_in  = tot_in*1.29/dn0(k) !Convert TO STP

       !Input aerosols in #/cm3 and outputs #/L activated
       if(iifn_formula==1) then
        !Original Demott(2010) formula
        ifntemp = 0.0000594 * (-tairc(k))**3.33 &
                * (tot_in)**(0.0264*(-tairc(k))+0.0033)
       elseif(iifn_formula==2) then
        !Modified Demott(2010) for dust-dominated cases
        !Paul suggested an additional factor of 3 multiplier
        ifntemp = 3.0 * 0.0008 * 10 ** (-0.2*(tairc(k)+9.7)) * tot_in ** 1.25
       endif

       !Adjust units and such
       ifntemp = ifntemp/1.29*dn0(k) !Convert FROM STP
       tot_in  = tot_in /1.29*dn0(k) !Convert FROM STP
       tot_in  = tot_in  / dn0(k) * 1.e6 !Convert #/cm3 to #/kg
       ifntemp = ifntemp / dn0(k) * 1.e3 !Convert #/L to #/kg
       ifntemp = min(tot_in,ifntemp)

       !Limit activation to only the number greater than that already
       !activated for the given grid cell parcel.
       ifntemp = ifntemp - ifnnucx(k)
       if(ifntemp < 0.0) ifntemp = 0.0

       !Determine freezing number and mixing ratio
       embtemp = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  & !mean mass
            / max(1.0e-6,cx(k,lcat))))
       frzc = min(cx(k,lcat),ifntemp)
       frzr = min(rx(k,lcat),frzc * embtemp)

       !Aerosol and solubility tracking
       !Transfering immersion freezing aerosol mass from cloud to pristine ice
       if(iccnlev>=2) then
        rxferratio = min(1.0,frzr/max(1.0e-12,rx(k,lcat)))
        ccnmass  = cnmhx(k,lcat) * rxferratio
        cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
        cnmhx(k,3)    = cnmhx(k,3)    + ccnmass
        if(itrkepsilon==1)then
         scnmass  = snmhx(k,lcat) * rxferratio
         snmhx(k,lcat) = snmhx(k,lcat) - scnmass
         snmhx(k,3)    = snmhx(k,3)    + scnmass
        endif
        if(itrkdust==1)then
         dcnmass  = dnmhx(k,lcat) * rxferratio
         dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
         dnmhx(k,3)    = dnmhx(k,3)    + dcnmass
        endif
        if(itrkdustifn==1)then
         !Add new immersion freezing dust mass to dust as IN tracking
         dinhx(k,3)    = dinhx(k,3)    + dcnmass
         !Xfer dust mass from liquid hydromet to ice if present
         dinmass  = dinhx(k,lcat) * rxferratio
         dinhx(k,lcat) = dinhx(k,lcat) - dinmass
         dinhx(k,3)    = dinhx(k,3)    + dinmass
        endif
       endif

       cx(k,lcat) = cx(k,lcat) - frzc
       rx(k,lcat) = rx(k,lcat) - frzr
       immerhx(k,lcat) = max(0.,immerhx(k,lcat)-frzc)
      endif

      !Update nucleation additive arrays for hazenuc and IN nuc
      vapnucr = vapnucr + frzr
      vapnuc  = vapnuc  + frzc
      immerin = immerin + frzc
      ifnnucx(k) = ifnnucx(k) + frzc

     endif
    enddo
   endif
   !***********************************************************************

   !Add new ice crystals to pristine ice category
   rx(k,3) = rx(k,3) + vapnucr
   cx(k,3) = cx(k,3) + vapnuc

   pcthaze = haznuc / max(1.0e-30,(haznuc + diagni + immerin))

   if(imbudget >= 1) then
     xnucicert(k) = xnucicert(k) + vapnucr * budget_scalet
   endif
   if(imbudget >= 2) then
     xinucifnrt(k) = xinucifnrt(k) + vapnucr * (1.0 - pcthaze) * budget_scalet
     xinuchazrt(k) = xinuchazrt(k) + vapnucr * pcthaze * budget_scalet
   endif

   if (rx(k,3) .ge. rxmin) k2pnuc = k
   if (k2pnuc .eq. 1 .and. rx(k,3) .lt. rxmin) k1pnuc = k + 1

enddo

! here mike has the habit diagnosis. option 1 is to use habit
! at cloud top, option 2 is to use new habit at each level.
! need to consider other options.  how about method of formation?
! my question about how much of habit is due to existing ice
! structure, and how much is due to current growth environment
! (temp and supsat). relevant supsat is wrt liquid?

return
END SUBROUTINE icenuc

!##############################################################################
Subroutine contnuc (rx,cx,tx,vap,press  &
   ,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt)

use micphys, only:rxmin

implicit none

real :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi  &
       ,dn1,dtlt,aka,raros,ana,akn,dfar,f1,f2,ft
data aka,raros/5.39e-3,3.e-7/

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
!
!  ana   = # IN per kg available for contact freezing (from Meyers et al. 1992
!          where ana was interpreted as # per m^3)
!  akn   = Knudsen number (Walko et al. 1995, Eq. 58)
!          [2.28e-5 = mfp * p00 / 293.15]
!  raros = aerosol radius = 3.e-7 m from Cotton et al. (1986)
!  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)
!          [7.32e-25 = Boltzmann constant / (6 pi)]
!  f1    = "func 1" (Walko et al. 1995 Eq. 55) multiplied by delta t
!           but now cld concen in #/kg_air so (pvbi, ptvi, pdvi) all per kg_air
!  f2    = "func 2" (Walko et al. 1995 Eq. 56)
!  ft    = "func ft" (Walko et al. 1995 Eq. 57)
!  pbvi  = Brownian motion nucleation amount this timestep [#/kg_air]
!  ptvi  = Thermophoretic nucleation amount this timestep [#/kg_air]
!  pdvi  = Diffusiophoretic nucleation amount this timestep [#/kg_air],
!          reformulated to use vapor diffusion directly.  Factor of 1.2
!          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102
!          divided by .622, the molecular weight ratio between water and air.

   ptotvi = 0.

   if (tx .le. -2. .and. rx .gt. rxmin) then

      ana = exp(4.11 - 0.262 * tx)
      akn = 2.28e-5 * tair / (press * raros)
      dfar = 7.32e-25 * tair * (1.+ akn) / (raros * dynvisc)
      f1 = 6.28318 * dn1 * cx * ana * dtlt
      f2 = thrmcon * (tairc - tx) / press
      ft = 0.4 * (1. + 1.45 * akn + 0.4 * akn * exp(-1. / akn))  &
         * (thrmcon + 2.5 * akn * aka)  &
         / ((1. + 3. * akn)  &
         * (2. * thrmcon + 5. * aka * akn + aka))
      pbvi = f1 * dfar
      ptvi = f1 * f2 * ft
      pdvi = 1.2 * ana * vap
      ptotvi = max(0.,pbvi + ptvi + pdvi)

   endif

return
END SUBROUTINE contnuc
