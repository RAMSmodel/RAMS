!##############################################################################
Subroutine thrmstr (m1,k1,k2,thp,theta,rtp,rv &
                   ,pp,pi0                    &
                   )

use rconstants
use micphys

implicit none

integer :: m1,k,lcat
real :: fracliq,tcoal,tairstr
integer, dimension(11) :: k1,k2
real, dimension(m1) :: thp,theta,rtp,rv
real, dimension(m1) :: pp,pi0

!OVER WHOLE COLUMN
do k = 1,m1
   pitot(k) = pi0(k) + pp(k)
   press(k) = p00 * (pitot(k) * cpi) ** cpor
   tair(k) = theta(k) * pitot(k) * cpi
enddo
!LOWER VAPOR LAYER
do k = 1,k1(11)-1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo
!UPPER VAPOR LAYER
do k = k2(11)+1,m1
   theta(k) = thp(k)
   rv(k) = rtp(k)
enddo
!ICE AND LIQUID LAYER
do k = k1(11),k2(11)
   til(k) = thp(k) * pitot(k) * cpi
   rliq(k) = 0.
   rice(k) = 0.
enddo
!CLOUD AND RAIN
do lcat = 1,2
   do k = k1(lcat),k2(lcat)
      rliq(k) = rliq(k) + rx(k,lcat)
   enddo
enddo
!DRIZZLE MODE
do lcat = 8,8
   do k = k1(lcat),k2(lcat)
      rliq(k) = rliq(k) + rx(k,lcat)
   enddo
enddo
!PRISTINE ICE, SNOW, AGGREGATES
do lcat = 3,5
   do k = k1(lcat),k2(lcat)
      rice(k) = rice(k) + rx(k,lcat)
   enddo
enddo
!GRAUPEL AND HAIL
!qtc gets temperature and fraction of liquid from q
do lcat = 6,7
   do k = k1(lcat),k2(lcat)
      CALL qtc (qx(k,lcat),tcoal,fracliq)
      rliq(k) = rliq(k) + rx(k,lcat) * fracliq
      rice(k) = rice(k) + rx(k,lcat) * (1. - fracliq)
   enddo
enddo
!FOR ICE/LIQUID LAYER get hydro heat and sat vap mix ratio
do k = k1(11),k2(11)
   qhydm(k) = alvl * rliq(k) + alvi * rice(k)
   rvstr(k) = rtp(k) - rliq(k) - rice(k)
   sa(k,1) = til(k) * qhydm(k) / max(1.0e-12,rliq(k) + rice(k))
enddo
!FOR ICE/LIQUID LAYER get temp of sat air in Celcius
do k = k1(11),k2(11)
   if (tair(k) .gt. 253.) then
      tairstr = 0.5 * (til(k)  &
         + sqrt(til(k) * (til(k) + cpi4 * qhydm(k))))
      sa(k,1) = sa(k,1) * cpi / (2. * tairstr - til(k))
   else
      tairstr = til(k) * (1. + qhydm(k) * cp253i)
      sa(k,1) = sa(k,1) * cp253i
  endif
  tairstrc(k) = tairstr - 273.15
enddo

return
END SUBROUTINE thrmstr

!##############################################################################
Subroutine diffprep (m1,lcat,k1,k2,rv,dn0)

use rconstants
use micphys

implicit none

integer :: m1,lcat,k1,k2,k,if1,if4,if6,if8,lhcat
real :: fre,scdei
real, dimension(m1) :: rv,dn0

!CODE BASED ON WALKO ET AL 2000
!EFFICIENT COMPUTATION OF VAPOR AND HEAT DIFFUSION BETWEEN HYDROMETEORS
!IN A NUMERICAL MODEL

!DETERMINES WHETHER TO CALCULATE USING LIQUID OR ICE "SA" ARRAYS
if ((lcat .le. 2) .or. (lcat .eq. 8)) then
   if1 = 1
   if4 = 4
   if6 = 6
   if8 = 8
else
   if1 = 2
   if4 = 5
   if6 = 7
   if8 = 9
endif

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) .lt. rxmin) go to 229

   fre = frefac1(lhcat) * emb(k,lcat) ** pwmasi(lhcat)  &
      + rdynvsci(k) * frefac2(lhcat) * emb(k,lcat) ** cdp1(lhcat)

   sb(k,lcat) = cx(k,lcat) * dn0(k) * fre * pi4dt
   su(k,lcat) = vapdif(k) * sb(k,lcat)

!zero for rain,graupel,hail
   sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
   se(k,lcat) = su(k,lcat) * sa(k,if6) + sb(k,lcat) * thrmcon(k)
   sf(k,lcat) = su(k,lcat) * sl(if1) - sb(k,lcat) * sa(k,2)
!q term used for rain,graupel,hail
   sg(k,lcat) = su(k,lcat) * sa(k,if8) + sb(k,lcat) * sa(k,3)  &
              + sj(lcat) * qr(k,lcat)
!reduces to 1/se if not rain,graupel,hail
   scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
   ss(k,lcat) = sf(k,lcat) * scdei
!reduces to sg*scdei if not all-liquid species
   sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
   ttest(k,lcat) = ss(k,lcat) * rv(k) + sw(k,lcat)

229    continue

enddo

!FOR ALL ICE HYDROS, "SM" IS 1
!IF PRISTINE,SNOW,AGG ABOVE OR= ZERO, "SH" IS 1
if (lcat .ge. 3 .and. lcat .le. 5) then
   do k = k1,k2
      if (rx(k,lcat) .lt. rxmin) go to 228
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
         sh(k,lcat) = 1.
         sd(k,lcat) = sh(k,lcat) * rx(k,lcat)
         scdei = 1. / (sc(if1) * sd(k,lcat) + se(k,lcat))
         ss(k,lcat) = sf(k,lcat) * scdei
         sw(k,lcat) = (sg(k,lcat) - sk(if1) * sd(k,lcat)) * scdei
      else
         sm(k,lcat) = 1.
      endif
228        continue
   enddo
endif

!FOR MIXED-PHASE HYDROMETEORS, "SM" IS 0
if (lcat .ge. 6 .and. lcat .le. 7) then
   do k = k1,k2
      if (rx(k,lcat) .lt. rxmin) go to 227
      if (ttest(k,lcat) .ge. 0.) then
         sm(k,lcat) = 0.
      else
         sm(k,lcat) = 1.
      endif
227        continue
   enddo
endif

do k = k1,k2
   if (rx(k,lcat) .lt. rxmin) go to 226
   sy(k,lcat) = rvsrefp(k,if1) * sm(k,lcat) * sw(k,lcat) - sa(k,if4)
   sz(k,lcat) = 1. - rvsrefp(k,if1) * ss(k,lcat) * sm(k,lcat)
   sumuy(k) = sumuy(k) + su(k,lcat) * sy(k,lcat)
   sumuz(k) = sumuz(k) + su(k,lcat) * sz(k,lcat)

226      continue
enddo

return
END SUBROUTINE diffprep

!##############################################################################
Subroutine vapdiff (m1,kf1,kf2,rv)

use micphys

implicit none

integer :: m1,kf1,kf2,k
real, dimension(m1) :: rv

!BASED ON FINAL EQUATION OF WALKO ET AL 2000
do k = kf1,kf2
  rv(k) = max(0.,(rvstr(k) + sumuy(k)) / (1.0 + sumuz(k)))
enddo

return
END SUBROUTINE vapdiff

!##############################################################################
Subroutine vapflux (m1,lcat,k1,k2,rv,i,j)

use micphys

implicit none

integer :: m1,lcat,k,k1,k2,if1,if4,i,j
real :: rxx
real, dimension(m1) :: rv

!UPDATES MIXING RATIO AND HEAT DUE TO FLUX OF VAPOR
!ALSO LINKED TO WALKO ET AL 2000

if ((lcat .le. 2) .or. (lcat .eq. 8)) then
   if1 = 1 
   if4 = 4
else
   if1 = 2
   if4 = 5
endif

do k = k1,k2

   if (rx(k,lcat) .lt. rxmin) go to 229
   rxx = 0.
   rxtemp = rx(k,lcat)
   cxtemp = cx(k,lcat)
   tx(k,lcat) = (ss(k,lcat) * rv(k) + sw(k,lcat)) * sm(k,lcat)
   vap(k,lcat) = su(k,lcat) * (rv(k) + sa(k,if4) - rvsrefp(k,if1) * tx(k,lcat))

   if (vap(k,lcat) .gt. -rx(k,lcat)) then

      rxx = rx(k,lcat) + vap(k,lcat)

      if (sm(k,lcat) .gt. .5) then
         qx(k,lcat) = sc(if1) * tx(k,lcat) + sk(if1)
         qr(k,lcat) = qx(k,lcat) * rxx
      else
         qx(k,lcat) = (rv(k) * sf(k,lcat) + sg(k,lcat)  &
                    - tx(k,lcat) * se(k,lcat)) / sd(k,lcat)
         qx(k,lcat) = min(350000.,max(-100000.,qx(k,lcat)))
         qr(k,lcat) = qx(k,lcat) * rxx
      endif

   endif

!Do the following section if pristine ice totally melts: evaporate it too.

!Saleeby(10/5/06): Aerosol restoration source.
!If total evap occurs, transfer aerosol mass within hydromet to regenerated 
!aerosol field and add aerosol number that is equal to the number of totally
!evaporated hydrometeors. We calculate the proportion of evaporated mass and
!return the proportion of aerosol mass contained with the hydrometeor species.
!We compute the median radius of restored particles with the assumption of
!log-normal distribution, and if the median radius < 0.96 microns this goes
!to small regenerated mode, otherwise large regenerated mode.

   rxferratio = 0.
   cxloss     = 0.

   !If we are doing FULL EVAPORATION
   if ((lcat .eq. 3 .and. qx(k,lcat) .gt. 330000.) .or.  &
      vap(k,lcat) .le. -rx(k,lcat)) then

      sumuy(k) = sumuy(k) - su(k,lcat) * sy(k,lcat)
      sumuz(k) = sumuz(k) - su(k,lcat) * sz(k,lcat)
      sumvr(k) = sumvr(k) + rx(k,lcat)
      rv(k) = max(0.,(rvstr(k) + sumuy(k) + sumvr(k)) / (1.0 + sumuz(k)))
      vap(k,lcat) = - rx(k,lcat)
      tx(k,lcat) = 0.
      rx(k,lcat) = 0.
      qx(k,lcat) = 0.
      qr(k,lcat) = 0.
      cxloss     = cx(k,lcat)
      cx(k,lcat) = 0.
      rxferratio = 1.

   !If we are doing VAPOR GROWTH or PARTIAL EVAPORATION
   else

      if(vap(k,lcat) .lt. 0.) then
        !Compute losses from partial evaporation
        fracmass = min(1.,-vap(k,lcat) / rx(k,lcat))
        cxloss = cx(k,lcat) * enmlttab( int(200.*fracmass)+1, jhcat(k,lcat) )
        cx(k,lcat) = cx(k,lcat) - cxloss
        rxferratio = rmlttab(int(200.*fracmass)+1,jhcat(k,lcat))
      endif

      rx(k,lcat) = rxx
   endif !if (full evaporation occurs) or (vapor growth or partial evap)

   !Aerosol and solubility tracking
   if(iccnlev>=2 .and. vap(k,lcat)<0.0 .and. rxferratio > .0001 .and. &
       cxloss > 0.0 .and. cnmhx(k,lcat) > 0.0) then
     ccnmass = cnmhx(k,lcat) * rxferratio
     cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass

     !For approximation, test rg based on regenerated aerosol
     acat=aerocat-1
     rg=((0.23873/aero_rhosol(acat)*ccnmass/cxloss)**(0.3333))/aero_rg2rm(acat)
     if(rg>0.96e-6) then
       acat=aerocat
     else
       acat=aerocat-1
     endif
     aeromas(k,acat) = aeromas(k,acat) + ccnmass
     aerocon(k,acat) = aerocon(k,acat) + cxloss

     if(itrkepsilon==1)then
      scnmass = snmhx(k,lcat) * rxferratio
      snmhx(k,lcat) = snmhx(k,lcat) - scnmass
      regenmas(k,acat-(aerocat-2)) = regenmas(k,acat-(aerocat-2)) + scnmass
     endif
     if(itrkdust==1)then
      dcnmass = dnmhx(k,lcat) * rxferratio
      dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
     endif
     if(itrkdustifn==1)then
      dinmass = dinhx(k,lcat) * rxferratio
      dinhx(k,lcat) = dinhx(k,lcat) - dinmass
     endif
   endif

   !Remove immersion freezing nuclei from cloud 
   if(iifn==3 .and. iccnlev>=1 .and. cxtemp>0.0) then
     if(lcat==1.or.lcat==8.or.lcat==2) then
       enxferratio = min(1.0,max(0.0,cxloss/cxtemp))
       ccnnum = immerhx(k,lcat) * enxferratio
       immerhx(k,lcat) = immerhx(k,lcat) - ccnnum
     endif
   endif

   !Vapor deposition and evaporation budgets for all species.
   !Deposition, condensation, evaporation, sublimation are all given positive
   ! values as process rates.
   if(imbudget >= 1) then
    if(lcat.eq.1 .or. lcat.eq.2 .or. lcat.eq.8) then
     if(vap(k,lcat) > 0.00) xvapliqt(k)  = xvapliqt(k)  + vap(k,lcat)*budget_scalet
     if(vap(k,lcat) < 0.00) xevapliqt(k) = xevapliqt(k) - vap(k,lcat)*budget_scalet
    endif
    if(lcat.ge.3 .and. lcat.le.7) then
     if(vap(k,lcat) > 0.00) xvapicet(k)  = xvapicet(k)  + vap(k,lcat)*budget_scalet
     if(vap(k,lcat) < 0.00) xevapicet(k) = xevapicet(k) - vap(k,lcat)*budget_scalet
    endif
   endif

   if(imbudget >= 2) then
    if(vap(k,lcat) > 0.00) then
     if(lcat==1) xvapcldt(k)  = xvapcldt(k)  + vap(k,lcat)*budget_scalet
     if(lcat==2) xvapraint(k) = xvapraint(k) + vap(k,lcat)*budget_scalet
     if(lcat==3) xvapprist(k) = xvapprist(k) + vap(k,lcat)*budget_scalet
     if(lcat==4) xvapsnowt(k) = xvapsnowt(k) + vap(k,lcat)*budget_scalet
     if(lcat==5) xvapaggrt(k) = xvapaggrt(k) + vap(k,lcat)*budget_scalet
     if(lcat==6) xvapgraut(k) = xvapgraut(k) + vap(k,lcat)*budget_scalet
     if(lcat==7) xvaphailt(k) = xvaphailt(k) + vap(k,lcat)*budget_scalet
     if(lcat==8) xvapdrizt(k) = xvapdrizt(k) + vap(k,lcat)*budget_scalet
    endif
    if(vap(k,lcat) < 0.00) then
     if(lcat==1) xevapcldt(k)  = xevapcldt(k)  - vap(k,lcat)*budget_scalet
     if(lcat==2) xevapraint(k) = xevapraint(k) - vap(k,lcat)*budget_scalet
     if(lcat==3) xevapprist(k) = xevapprist(k) - vap(k,lcat)*budget_scalet
     if(lcat==4) xevapsnowt(k) = xevapsnowt(k) - vap(k,lcat)*budget_scalet
     if(lcat==5) xevapaggrt(k) = xevapaggrt(k) - vap(k,lcat)*budget_scalet
     if(lcat==6) xevapgraut(k) = xevapgraut(k) - vap(k,lcat)*budget_scalet
     if(lcat==7) xevaphailt(k) = xevaphailt(k) - vap(k,lcat)*budget_scalet
     if(lcat==8) xevapdrizt(k) = xevapdrizt(k) - vap(k,lcat)*budget_scalet
    endif
   endif

229     continue

enddo

return
END SUBROUTINE vapflux

!##############################################################################
Subroutine psxfer (k13,k23,k14,k24,i,j)

use micphys

implicit none

integer :: k1,k2,k13,k23,k14,k24,i,j,k,lhcat,it
real :: embx,dn,xlim,dvap,dqr,dnum
real :: old_m,old_c,old_r,prelim_m,delta_c, delta_r

k1 = min(k13,k14)
k2 = max(k23,k24)

do k = k1,k2

   if (vap(k,3) .gt. 0. .or. vap(k,4) .lt. 0.) then

      if (vap(k,3) .gt. 0.) then
         lhcat = jhcat(k,3)
         embx = max(rxmin,rx(k,3)) / max(cxmin,cx(k,3))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = min(5000,max(1,nint(dn * 1.e6)))

         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.)  &
            / (gamn1(3) * pwmas(lhcat) * dn ** 2)

         dvap = min(rx(k,3),  &
                    vap(k,3) * (xlim + gam(it,1) / gamn1(3)))
         dqr = dvap * qx(k,3)
         dnum = min(cx(k,3),dvap * min(dpsmi(lhcat),1./embx))

         !Aerosol and solubility tracking
         !Transfer aerosol-in-hydrometeor masses between categories
         if(iccnlev>=2) then
          ccnmass=0.0
          scnmass=0.0
          dcnmass=0.0
          dinmass=0.0
          if(rx(k,3)>0.0 .and. dvap>0.0) then
           rxferratio = min(1.0 , dvap / rx(k,3))
           ccnmass = cnmhx(k,3) * rxferratio
           if(itrkepsilon==1) scnmass = snmhx(k,3) * rxferratio
           if(itrkdust==1)    dcnmass = dnmhx(k,3) * rxferratio
           if(itrkdustifn==1) dinmass = dinhx(k,3) * rxferratio
          endif
         endif

      else
         lhcat = jhcat(k,4)
         embx = max(rxmin,rx(k,4)) / max(cxmin,cx(k,4))
         dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
         it = min(5000,max(1,nint(dn * 1.e6)))
         
         xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.)  &
            / (gamn1(4) * pwmas(lhcat) * dn ** 2)

         dvap = max(-rx(k,4),vap(k,4) * xlim)
         dqr = dvap * qx(k,4)
         dnum = max(-cx(k,4),dvap * max(dpsmi(lhcat),1./embx))

         !Aerosol and solubility tracking
         !Transfer aerosol-in-hydrometeor masses between categories
         if(iccnlev>=2) then
          ccnmass=0.0
          scnmass=0.0
          dcnmass=0.0
          dinmass=0.0
          if(rx(k,4)>0.0 .and. dvap<0.0) then
           rxferratio = max(-1.0 , dvap / rx(k,4))
           ccnmass = cnmhx(k,4) * rxferratio
           if(itrkepsilon==1) scnmass = snmhx(k,4) * rxferratio
           if(itrkdust==1)    dcnmass = dnmhx(k,4) * rxferratio
           if(itrkdustifn==1) dinmass = dinhx(k,4) * rxferratio
          endif
         endif

      endif

      !Aerosol and solubility tracking
      !Transfer aerosol-in-hydrometeor masses between categories
      if(iccnlev>=2) then
        cnmhx(k,3) = cnmhx(k,3) - ccnmass
        cnmhx(k,4) = cnmhx(k,4) + ccnmass
        if(itrkepsilon==1) then
         snmhx(k,3) = snmhx(k,3) - scnmass
         snmhx(k,4) = snmhx(k,4) + scnmass
        endif
        if(itrkdust==1) then
         dnmhx(k,3) = dnmhx(k,3) - dcnmass
         dnmhx(k,4) = dnmhx(k,4) + dcnmass
        endif
        if(itrkdustifn==1) then
         dinhx(k,3) = dinhx(k,3) - dinmass
         dinhx(k,4) = dinhx(k,4) + dinmass
        endif
      endif

      rx(k,3) = rx(k,3) - dvap
      cx(k,3) = cx(k,3) - dnum
      qr(k,3) = qr(k,3) - dqr
      rx(k,4) = rx(k,4) + dvap
      cx(k,4) = cx(k,4) + dnum
      qr(k,4) = qr(k,4) + dqr

      !Make sure xfers do not create negative values.
      if(rx(k,3)<0.0.or.rx(k,4)<0.0.or.cx(k,3)<0.0.or.cx(k,4)<0.0) then
        print*,'Pristine ice or snow are negative at PSXFER',k,i,j
        print*,rx(k,3),rx(k,4),cx(k,3),cx(k,4),dvap,dnum,vap(k,3),vap(k,4)
        stop
      endif

      !Carrio 2003: Better xfer pristine to snow to prevent routine
      !enemb from artificially creating pristine number concentration
      !when a bounds problem appears with the pristine ice cutoff size.
      !Compare preliminary calculation pristine mean mass with maximum.
      delta_r = 0.0
      prelim_m = 0.0
      if(cx(k,3) > cxmin) prelim_m = rx(k,3) / cx(k,3)
      if (prelim_m .gt. emb1(3)) then
          old_m = (rx(k,3) + dvap) / (cx(k,3) + dnum)
          old_c = cx(k,3) + dnum
          old_r = rx(k,3) + dvap
          delta_r = rx(k,3) - old_c * emb1(3)
          delta_c = delta_r / emb1(3) ! delta_c only adds to snow

          !Aerosol and solubility tracking
          if(iccnlev>=2) then
           ccnmass=0.0
           scnmass=0.0
           dcnmass=0.0
           dinmass=0.0
          endif

          if(delta_r>0.0) then !xfer more Pristine ice to Snow
            delta_r = min(rx(k,3),delta_r)
            delta_c = delta_r / emb1(3)
            !Aerosol and solubility tracking
            if(iccnlev>=2) then 
             if(rx(k,3)>0.0) then
              rxferratio = min(1.0 , delta_r / rx(k,3))
              ccnmass = cnmhx(k,3) * rxferratio
              if(itrkepsilon==1) scnmass = snmhx(k,3) * rxferratio
              if(itrkdust==1)    dcnmass = dnmhx(k,3) * rxferratio
              if(itrkdustifn==1) dinmass = dinhx(k,3) * rxferratio
             endif
            endif
          elseif(delta_r<0.0) then !xfer more Snow to Pristine ice
            delta_r = max(-rx(k,4),delta_r)
            delta_c = delta_r / emb1(3)
            !Aerosol and solubility tracking
            if(iccnlev>=2) then 
             if(rx(k,4)>0.0) then
              rxferratio = max(-1.0 , delta_r / rx(k,4))
              ccnmass = cnmhx(k,4) * rxferratio
              if(itrkepsilon==1) scnmass = snmhx(k,4) * rxferratio
              if(itrkdust==1)    dcnmass = dnmhx(k,4) * rxferratio
              if(itrkdustifn==1) dinmass = dinhx(k,4) * rxferratio
             endif
            endif
          endif

          !Aerosol and solubility tracking
          if(iccnlev>=2 .and. delta_r /= 0) then
            cnmhx(k,3) = cnmhx(k,3) - ccnmass
            cnmhx(k,4) = cnmhx(k,4) + ccnmass
            if(itrkepsilon==1) then
             snmhx(k,3) = snmhx(k,3) - scnmass
             snmhx(k,4) = snmhx(k,4) + scnmass
            endif
            if(itrkdust==1) then
             dnmhx(k,3) = dnmhx(k,3) - dcnmass
             dnmhx(k,4) = dnmhx(k,4) + dcnmass
            endif
            if(itrkdustifn==1) then
             dinhx(k,3) = dinhx(k,3) - dinmass
             dinhx(k,4) = dinhx(k,4) + dinmass
            endif
          endif

          rx(k,3) = rx(k,3) - delta_r
          cx(k,3) = old_c
          rx(k,4) = rx(k,4) + delta_r
          cx(k,4) = cx(k,4) + delta_c
      endif

   endif
enddo

!Re-determine pris and snow levels since xfers occurred
k13=2
k14=2
k23=1
k24=1
do k = k1,k2
  if(rx(k,3) >= rxmin) then
    k23 = k
  else
    if(k23 == 1) k13 = k + 1
  endif
  if(rx(k,4) >= rxmin) then
    k24 = k
  else
    if(k24 == 1) k14 = k + 1
  endif
enddo

return
END SUBROUTINE psxfer

!##############################################################################
Subroutine newtemp (m1,kf1,kf2,rv,theta)

use rconstants
use micphys

implicit none

real, external :: rslf
real, external :: rsif

integer :: m1,kf1,kf2,k
real, dimension(m1) :: rv,theta

do k = kf1,kf2
   tairc(k) = tairstrc(k) + sa(k,1) * (rvstr(k) - rv(k))
   tair(k)  = tairc(k) + 273.15
   theta(k) = tair(k) * cp / pitot(k)

   rvlsair(k) = rslf(press(k),tair(k))
   rvisair(k) = rsif(press(k),tair(k))
enddo

return
END SUBROUTINE newtemp
