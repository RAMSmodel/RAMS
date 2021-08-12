!##############################################################################
Subroutine getict (k1,k2,lcat)

use micphys

implicit none

integer :: k1,k2,lcat,k
real :: rict,rictmm

do k = k1,k2
  if (rx(k,lcat) .ge. rxmin) then
   rict = dict(lcat) * (log(emb(k,lcat)) - emb0log(lcat)) + 1.
   rictmm = max(rictmin,min(rictmax,rict))
   ict1(k,lcat) = int(rictmm)
   ict2(k,lcat) = ict1(k,lcat) + 1
   wct2(k,lcat) = rictmm - float(ict1(k,lcat))
   wct1(k,lcat) = 1.0 - wct2(k,lcat)
  endif
enddo

return
END SUBROUTINE getict

!##############################################################################
Subroutine auto_accret (m1,k1,k2,k3,k4,dn0,dtlt,i,j)

use micphys

implicit none

integer :: m1,k,k1,k2,k3,k4,kbot,ktop &
          ,id1cc,id2dd,id1cd,id2cd,id1cr,id1crn,id2cr,id2crn &
          ,ir3cr,id3cr,i,j

real :: dtlt,dtlt3,dtlt6,dmb1cgs,dmb2cgs,dmb3cgs,r3cgs,en1cgs,en2cgs &
   ,ad1,ad2,ar3,ad3,bd1,bd2,br3,bd3,d3e  &    
   ,bd1cc,bd2dd,bd1cd,bd2cd,bd1cr,bd2cr,br3cr,bd3cr &    
   ,wd1cc,wd2dd,wd1cd,wd2cd,wd1cr,wd2cr,wr3cr &
   ,tm1cc,tn1cc,tn2cc,tm2dd,tn2dd,tn3dd,tm1cd,tn1cd,tm2cd,tn2cd &
   ,tm1cr,tn1cr,tm2cr,tn2cr &
   ,en1cgs_2,en2cgs_2,um1cc,un1cc,un2cc,um2dd,un2dd,un3dd,um1cd &
   ,un1cd,um2cd,un2cd,um1cr,un1cr,um2cr,un2cr,um12,um82 &
   ,um18,un1,un8,cfmasi1,cfmasi2,cfmasi3,pwmasi1,pwmasi2,pwmasi3

real, dimension(m1) :: dn0

dtlt3 = 1.e3 * dtlt
dtlt6 = 1.e6 * dtlt

cfmasi1 = 1. / cfmas(1)
cfmasi2 = 1. / cfmas(16)
cfmasi3 = 1. / cfmas(2)

pwmasi1 = 1. / pwmas(1)
pwmasi2 = 1. / pwmas(16)
pwmasi3 = 1. / pwmas(2)

if(jnmb(8) .ne. 0) then
  kbot=min(k1,k3)
  ktop=max(k2,k4)
else
  kbot=k1
  ktop=k2
endif

do k = kbot,ktop
if(rx(k,1) .ge. rxmin .or. rx(k,8) .ge. rxmin) then

! This routine works in cgs units, so convert inputs from mks
!mean diameter cloud,cloud2,rain - convert (meters) to (cm)
   dmb1cgs = 100. * (emb(k,1) * cfmasi1) ** pwmasi1
   dmb2cgs = 100. * (emb(k,8) * cfmasi2) ** pwmasi2
   dmb3cgs = 100. * (emb(k,2) * cfmasi3) ** pwmasi3

!mixing ratio rain - convert (kg/kg) to (g/cm3)
   r3cgs  = 1.0e-3 * rx(k,2) * dn0(k)
   if(rx(k,2) .lt. rxmin) r3cgs = 0.0

!number concentration cloud 1 & 2 - convert (#/kg) to (#/cm3)
   en1cgs = 1.0e-6 * cx(k,1) * dn0(k)
   en2cgs = 1.0e-6 * cx(k,8) * dn0(k)

!max diameter of cloud 1 & 2
   ad1 = max(d1min,min(d1max,dmb1cgs))
   ad2 = max(d2min,min(d2max,dmb2cgs))

!max mixing ratio rain
   ar3 = max(r3min,min(r3max,r3cgs))

!max diameter of rain
   d3minx = max(d3min,(r3cgs / (.1 * .5236)) ** pwmasi3)
   ad3 = max(d3minx,min(d3max,dmb3cgs))

!log of CLOUD 1 & 2 diameter range ratio
   bd1 = alog10(ad1/d1min)
   bd2 = alog10(ad2/d2min)
!log of RAIN mix ratio range ratio
   br3 = alog10(ar3/r3min)
!log of RAIN diamter range ratio
   bd3 = alog10(ad3/d3minx)
!log of RAIN diamter of ratio of given max to max from dmin and mix ratio
   d3e = alog10(d3max/d3minx)

!Calculate ratios of true diameter and mixratio to range in (mkautotab)
!to determine which bin on curve is closest
   bd1cc   = float(ndcc-1) * min(1.0,(ad1 - d1min) / (d1max - d1min)) + 1.
   bd2dd   = float(ndcc-1) * min(1.0,(ad2 - d2min) / (d2max - d2min)) + 1.
   bd1cd   = float(ndcd-1) * min(1.0,(ad1 - d1min) / (d1max - d1min)) + 1.
   bd2cd   = float(ndcd-1) * min(1.0,(ad2 - d2min) / (d2max - d2min)) + 1.
   bd1cr   = bd1 / d1ecr + 1.
   bd2cr   = bd2 / d2ecr + 1.
   br3cr   = br3 / r3ecr + 1.
   bd3cr   = bd3 / d3e * float(ndrcr-1) + 1.

!Find closest integer bin
   id1cc   =  int(bd1cc)
   id2dd   =  int(bd2dd)
   id1cd   =  int(bd1cd)
   id2cd   =  int(bd2cd)
   id1cr   =  int(bd1cr)
   id1crn  = nint(bd1cr)
   id2cr   =  int(bd2cr)
   id2crn  = nint(bd2cr)
   ir3cr   =  int(br3cr)
   id3cr   = nint(bd3cr)

!Find distance to closest integer bin
   wd1cc   = bd1cc   - float(id1cc)
   wd2dd   = bd2dd   - float(id2dd)
   wd1cd   = bd1cd   - float(id1cd)
   wd2cd   = bd2cd   - float(id2cd)
   wd1cr   = bd1cr   - float(id1cr)
   wd2cr   = bd2cr   - float(id2cr)
   wr3cr   = br3cr   - float(ir3cr)

!To fix array bounds exceptions
   if(id1cc.gt.ndcc-1) then
     id1cc=ndcc-1
     wd1cc=1.-wd1cc
   endif
   if(id2dd.gt.ndcc-1) then
     id2dd=ndcc-1
     wd2dd=1.-wd2dd
   endif
   if(id1cr.gt.ndccr-1) then
     id1cr=ndccr-1
     wd1cr=1.-wd1cr
   endif
   if(ir3cr.gt.nrrcr-1) then
     ir3cr=nrrcr-1
     wr3cr=1.-wr3cr
   endif
   if(id1cd.gt.ndcd-1)  then
     id1cd=ndcd-1
     wd1cd=1.-wd1cd
   endif
   if(id2cd.gt.ndcd-1)  then
     id2cd=ndcd-1
     wd2cd=1.-wd2cd
   endif
   if(id2cr.gt.ndccr-1) then
     id2cr=ndccr-1
     wd2cr=1.-wd2cr
   endif

!***************************************************************************
!cc-> effect m,n of cloud,drizzle and possibly n of rain
tm1cc =   (1.-wd1cc) * r1tabcc(id1cc) + wd1cc * r1tabcc(id1cc+1)
tn1cc =   (1.-wd1cc) * c1tabcc(id1cc) + wd1cc * c1tabcc(id1cc+1)
tn2cc =   (1.-wd1cc) * c2tabcc(id1cc) + wd1cc * c2tabcc(id1cc+1)
!***************************************************************************
!cr-> effect m,n of cloud on rain
tm1cr = (1.-wd1cr) * ((1.-wr3cr) * r1tabcr(id1cr  ,ir3cr  ,id3cr)  &
      +                   wr3cr  * r1tabcr(id1cr  ,ir3cr+1,id3cr)) &
      +     wd1cr  * ((1.-wr3cr) * r1tabcr(id1cr+1,ir3cr  ,id3cr)  &
      +                   wr3cr  * r1tabcr(id1cr+1,ir3cr+1,id3cr))

tn1cr =               (1.-wr3cr) * c1tabcr(id1crn,ir3cr  ,id3cr)   &
      +                   wr3cr  * c1tabcr(id1crn,ir3cr+1,id3cr)
!***************************************************************************

if(jnmb(8) .ne. 0) then
!***************************************************************************
!dd-> effect m,n of drizzle and n on rain
tm2dd =  (1.-wd2dd) * r2tabdd(id2dd) + wd2dd * r2tabdd(id2dd+1)
tn2dd =  (1.-wd2dd) * c2tabdd(id2dd) + wd2dd * c2tabdd(id2dd+1)
tn3dd =  (1.-wd2dd) * c3tabdd(id2dd) + wd2dd * c3tabdd(id2dd+1)
!***************************************************************************
!cd-> effect m,n of cloud,drizzle and possibly n on rain
!For cloud mixing ratio portion to rain
tm1cd = (1.-wd2cd) * ((1.-wd1cd) * r1tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * r1tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * r1tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * r1tabcd(id1cd+1,id2cd+1))
!For cloud loss of number
tn1cd = (1.-wd2cd) * ((1.-wd1cd) * c1tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * c1tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * c1tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * c1tabcd(id1cd+1,id2cd+1))
!For drizzle mixing ratio portion to rain
tm2cd = (1.-wd2cd) * ((1.-wd1cd) * r2tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * r2tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * r2tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * r2tabcd(id1cd+1,id2cd+1))
!For drizzle loss of number to rain
tn2cd = (1.-wd2cd) * ((1.-wd1cd) * c2tabcd(id1cd  ,id2cd  )  &
       +                  wd1cd  * c2tabcd(id1cd+1,id2cd  )) &
       +    wd2cd  * ((1.-wd1cd) * c2tabcd(id1cd  ,id2cd+1)  &
       +                  wd1cd  * c2tabcd(id1cd+1,id2cd+1))
!***************************************************************************
!dr-> effect m,n of drizzle
tm2cr =  (1.-wd2cr) * ((1.-wr3cr) * r2tabcr(id2cr  ,ir3cr  ,id3cr)  &
      +                    wr3cr  * r2tabcr(id2cr  ,ir3cr+1,id3cr)) &
      +      wd2cr  * ((1.-wr3cr) * r2tabcr(id2cr+1,ir3cr  ,id3cr)  &
      +                    wr3cr  * r2tabcr(id2cr+1,ir3cr+1,id3cr))

tn2cr =                (1.-wr3cr) * c2tabcr(id2crn,ir3cr  ,id3cr)  &
      +                    wr3cr  * c2tabcr(id2crn,ir3cr+1,id3cr)
!***************************************************************************
else
!Prevent uninitialized
tm2dd=0.; tn2dd=0.; tn3dd=0.; tm1cd=0.; tn1cd=0.
tm2cd=0.; tn2cd=0.; tm2cr=0.; tn2cr=0.
endif

!Put tables values in correct units
   en1cgs_2 = en1cgs ** 2
   en2cgs_2 = en2cgs ** 2

   um1cc = tm1cc * en1cgs_2 * dtlt3
   un1cc = tn1cc * en1cgs_2 * dtlt6
   un2cc = tn2cc * en1cgs_2 * dtlt6

   um1cr = 10. ** tm1cr * en1cgs * dtlt3
   un1cr = 10. ** tn1cr * en1cgs * dtlt6
   !un3rr = 10. ** tn3rr * dtlt6 !rain-rain done in "cols"
   !Saleeby(1-13-06) If rain mixing ratio is zero, do not collect cloud
   if(r3cgs == 0.) then
     um1cr=0.
     un1cr=0.
   endif

  if(jnmb(8) .ne. 0) then
   um2dd = tm2dd * en2cgs_2 * dtlt3
   un2dd = tn2dd * en2cgs_2 * dtlt6
   un3dd = tn3dd * en2cgs_2 * dtlt6
   um1cd = tm1cd * en1cgs * en2cgs * dtlt3
   un1cd = tn1cd * en1cgs * en2cgs * dtlt6
   um2cd = tm2cd * en1cgs * en2cgs * dtlt3
   un2cd = tn2cd * en1cgs * en2cgs * dtlt6
   um2cr = 10. ** tm2cr * en2cgs * dtlt3
   un2cr = 10. ** tn2cr * en2cgs * dtlt6
   !Saleeby(1-13-06) If rain mixing ratio is zero, do not collect drizzle
   if(r3cgs == 0.) then
     um2cr=0.
     un2cr=0.
   endif
  endif

! The above values are amounts in kg/m^3 or #/m^3 converted in the
! present timestep, but must still be corrected for the effect of
! density on fall velocity.  Thus, they must be multiplied by
! (dn0i ** .5) which fall velocity is proportional to.  Also, since
! rxfer and enxfer are in units of kg/kg and #/kg, respectively, the
! above transfer amounts must also be multiplied by dn0i.  Together,
! these factors make (dn0i ** 1.5).

!****************** FOR CLOUD1 - CLOUD2 - RAIN TRANSFERS *********************
if(jnmb(8) .eq. 0) then
   um12 = min(rx(k,1),(um1cc + um1cr) * dn0i(k))
   un1 = min(cx(k,1)*dn0(k),(un1cc + un1cr))
   rxfer(k,1,2)  =  rxfer(k,1,2) + um12
   qrxfer(k,1,2) = qrxfer(k,1,2) + um12 * qx(k,1)
   enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
   enxfer(k,1,2) = enxfer(k,1,2) + un2cc

   if(imbudget >= 1) xcld2raint(k) = xcld2raint(k) + um12 * budget_scalet
endif

if(jnmb(8) .ne. 0) then
   if(um2cd <= 0.0) then
     um12 = min(rx(k,1),(um1cr + um1cd + um2cd) * dn0i(k))
   else
     um12 = min(rx(k,1),(um1cr + um1cd) * dn0i(k))
   endif

   if(um2cd <= 0.0) then
     um82 = min(rx(k,8),(um2cr + um2dd) * dn0i(k))
   else
     um82 = min(rx(k,8),(um2cr + um2dd + um2cd) * dn0i(k))
   endif

   if(um2cd <= 0.0) then
     um18 = min(rx(k,1),(um1cc + abs(um2cd)) * dn0i(k))
   else
     um18 = min(rx(k,1),(um1cc) * dn0i(k))
   endif

    un1 = min(cx(k,1)*dn0(k),(un1cc + un1cd + un1cr))
    un8 = min(cx(k,8)*dn0(k),(un2dd + un2cd + un2cr))

    rxfer(k,1,2) =  rxfer(k,1,2) + um12
    rxfer(k,8,2) =  rxfer(k,8,2) + um82
    rxfer(k,1,8) =  rxfer(k,1,8) + um18
   qrxfer(k,1,2) = qrxfer(k,1,2) + um12 * qx(k,1)
   qrxfer(k,8,2) = qrxfer(k,8,2) + um82 * qx(k,8)
   qrxfer(k,1,8) = qrxfer(k,1,8) + um18 * qx(k,1)
   enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
   enxfer(k,1,8) = enxfer(k,1,8) + un2cc
   enxfer(k,8,8) = enxfer(k,8,8) + un8 - un3dd - un2cd
   enxfer(k,8,2) = enxfer(k,8,2) + un3dd + un2cd

   if(imbudget >= 1) xcld2raint(k) = xcld2raint(k) + (um12 + um82) * budget_scalet
endif

!USING CLOUD_2 TRANSFER TO TRANSFER CC2 COLLISION NUMBER TO
!RAIN SINCE YOU CAN'T JUST SUM THE LOSS OF NUMBER OF CLOUD1 AND CLOUD2
!CLOUD1 COLLISIONS WITH ITSELF DOESN'T AFFECT RAIN NUMBER.

endif !if cloud mixing ratio greater than min threshold
enddo !loop of vertical levels

return
END SUBROUTINE auto_accret

!##############################################################################
Subroutine auto_accret_ice (m1,jcat,lcat,k1,k2,dn0,dtlt)

use micphys

implicit none

integer :: m1,k,k1,k2,id1cr,id1crn,irici,idici &
          ,lcat,jhcaty,mx,it,ccat,jcat

real :: dtlt,dtlt3,dtlt6,cfmasi1,pwmasi1,dmb1cgs,dmbicgs &
   ,ricgs,en1cgs,ad1,ari,adi,bd1,bri,bdi,die,bd1cr,brici,bdici  &
   ,wd1cr,wrici,tm1ci,tn1ci,um1ci,un1ci,umcld,uncld,trime,urime,rimer &
   ,qcoal,tcoal,fracliq,area,cn13,cn24,sip,rsip,qrsip,coalliq   &
   ,rfinlz,xtoz,xtoy,ytoz,dcmin,dcmax,dcecr,rcoal

real, dimension(m1) :: dn0

integer, dimension(8) :: mcatc

data mcatc /0,0,0,6,6,7,7,0/

dtlt3 = 1.e3 * dtlt
dtlt6 = 1.e6 * dtlt

!Prevent uninitialized
cfmasi1=0.; pwmasi1=0.; dcmin=0.; dcmax=0.; dcecr=0.

if(jcat==1)then
  cfmasi1 = 1. / cfmas(1)
  pwmasi1 = 1. / pwmas(1)
  dcmin = d1min
  dcmax = d1max
  dcecr = d1ecr
elseif(jcat==8)then
  cfmasi1 = 1. / cfmas(16)
  pwmasi1 = 1. / pwmas(16)
  dcmin = d2min
  dcmax = d2max
  dcecr = d2ecr
endif

do k = k1,k2
if(rx(k,jcat).ge.rxmin .and. rx(k,lcat).ge.rxmin) then
if(  ((lcat==4 .or. lcat==5) .and. emb(k,jcat) .gt. 9.0e-13) .or. &
     ((lcat==6 .or. lcat==7) .and. emb(k,jcat) .gt. 3.4e-14) ) then

! This routine works in cgs units, so convert inputs from mks
!mean diameter cloud,cloud2,rain - convert (meters) to (cm)
   dmb1cgs = 100. * (emb(k,jcat) * cfmasi1) ** pwmasi1
   dmbicgs = 100. * (emb(k,lcat) / cfmas(lcat)) ** (1./pwmas(lcat))

!mixing ratio ice species - convert (kg/kg) to (g/cm3)
   ricgs  = 1.0e-3 * rx(k,lcat) * dn0(k)

!number concentration cloud 1 & 2 - convert (#/kg) to (#/cm3)
   en1cgs = 1.0e-6 * cx(k,jcat) * dn0(k)

!max diameter of cloud 1 & 2
   ad1 = max(dcmin,min(dcmax,dmb1cgs))

!max mixing ratio ice species
   ari = max(rimin,min(rimax,ricgs))

!max diameter of ice species
   diminx=max(dimin,((ricgs/1000./cfmas(lcat))**(1./pwmas(lcat)))*200.)
   adi = max(diminx,min(dimax,dmbicgs))

!log of CLOUD 1 & 2 diameter range ratio
   bd1 = alog10(ad1/dcmin)
!log of ICE mix ratio range ratio
   bri = alog10(ari/rimin)
!log of ICE diamter range ratio
   bdi = alog10(adi/diminx)
!log of ICE diamter of ratio of given max to max from dmin and mix ratio
   die = alog10(dimax/diminx)

!Calculate ratios of true diameter and mixratio to range in (mkautotab)
!to determine which bin on curve is closest
   bd1cr   = bd1 / dcecr + 1.
   brici   = bri / rieci + 1.
   bdici   = bdi / die * float(ndrcr-1)  + 1.

!Find closest integer bin
   id1cr   =  int(bd1cr)
   id1crn  = nint(bd1cr)
   irici   =  int(brici)
   idici   = nint(bdici)

!Find distance to closest integer bin
   wd1cr  = bd1cr  - float(id1cr)
   wrici  = brici  - float(irici)

!To fix array bounds exceptions
   if(id1cr.gt.ndccr-1) then
     id1cr=ndccr-1
     wd1cr=1.-wd1cr
   endif
   if(irici.gt.nrrcr-1) then
     irici=nrrcr-1
     wrici=1.-wrici
   endif

!Prevent uninitialized
tm1ci=0.; tn1ci=0.; trime=0.

!***************************************************************************
!ci-> effect m,n of cloud on ice
if(jcat==1)then
 tm1ci = (1.-wd1cr) * ((1.-wrici) * r1tabci(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r1tabci(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r1tabci(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r1tabci(id1cr+1,irici+1,idici,lcat-3))
 tn1ci =               (1.-wrici) * c1tabci(id1crn,irici  ,idici,lcat-3)   &
       +                   wrici  * c1tabci(id1crn,irici+1,idici,lcat-3)
 trime = (1.-wd1cr) * ((1.-wrici) * r1rimer(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r1rimer(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r1rimer(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r1rimer(id1cr+1,irici+1,idici,lcat-3))
endif
!***************************************************************************

!***************************************************************************
!c2i-> effect m,n of cloud2 on ice
if(jcat==8)then
 tm1ci = (1.-wd1cr) * ((1.-wrici) * r2tabci(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r2tabci(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r2tabci(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r2tabci(id1cr+1,irici+1,idici,lcat-3))
 tn1ci =               (1.-wrici) * c2tabci(id1crn,irici  ,idici,lcat-3)   &
       +                   wrici  * c2tabci(id1crn,irici+1,idici,lcat-3)
 trime = (1.-wd1cr) * ((1.-wrici) * r2rimer(id1cr  ,irici  ,idici,lcat-3)  &
       +                   wrici  * r2rimer(id1cr  ,irici+1,idici,lcat-3)) &
       +     wd1cr  * ((1.-wrici) * r2rimer(id1cr+1,irici  ,idici,lcat-3)  &
       +                   wrici  * r2rimer(id1cr+1,irici+1,idici,lcat-3))
endif
!***************************************************************************

 um1ci = 10. ** tm1ci * en1cgs * dtlt3
 un1ci = 10. ** tn1ci * en1cgs * dtlt6
 urime = 10. ** trime * en1cgs * dtlt3

!For Cloud-ice collisions 4=snow, 5=aggregates, 6=graupel, 7=hail
 umcld = min(rx(k,jcat),um1ci * dn0i(k))
 uncld = min(cx(k,jcat)*dn0(k),un1ci)
 rimer = min(rx(k,lcat),urime * dn0i(k))

!*****************************************************************
! Do secondary ice production and transfers due to collection
!****************************************************************
 mx=0
 if(jcat==1) mx=1 !uses gamsip(1,it) for cloud1 set in mic_init.f90
 if(jcat==8) mx=2 !uses gamsip(2,it) for cloud2 set in mic_init.f90

 !Compute fraction liquid of combined ice and cloud1 droplets
 rcoal = umcld + rimer
 qcoal = (qx(k,jcat) * umcld + qx(k,lcat) * rimer) / max(1.0e-20,rcoal)
 
 CALL qtc (qcoal,tcoal,fracliq)

 coalliq = rcoal * fracliq

 !Secondary Ice Production Based on Hydrometeor Internal Temperature
 if(tcoal.gt.-8.0.and.tcoal.lt.-3.0 .and. (jnmb(jcat)>=5.or.jnmb(lcat)>=5))then
   jhcaty=jhcat(k,lcat)
   area = cx(k,lcat) * dn0(k) * sipfac(jhcaty) * emb(k,lcat)  &
      ** (2.*pwmasi(jhcaty))
   it = nint(emb(k,jcat) / emb1(jcat) * 5000.)
   cn13 = uncld * gamsip13(mx,it) / (area * dtlt)
   cn24 = min(cx(k,jcat)*dn0(k),uncld) * gamsip24(mx,it)
   sip = 9.1e-10 * cn24 * cn13 ** .93
   if (tcoal .lt. -5.) then
      sip = 0.33333 * (tcoal + 8.) * sip
   else
      sip = -0.5 * (tcoal + 3.) * sip
   endif
   rsip = sip * emb0(3) * dn0i(k)
   qrsip = qcoal * rsip
   rcoal = rcoal - rsip
   enxfer(k,jcat,3) = enxfer(k,jcat,3) + sip
   rxfer(k,jcat,3)  = rxfer(k,jcat,3) + rsip
   qrxfer(k,jcat,3) = qrxfer(k,jcat,3) + qrsip
 endif

 !For Cloud-ice collisions 4=snow, 5=aggregates, 6=graupel, 7=hail

 !Saleeby(6-22-2015) Modified the amount to send to destination
 !category after implementing density weighting of rain-ice collisions
 !to determine destination species. "rfinlz" looks more like "col2" now.
 !This does not let snow/aggregates keep rimed water. Perhaps a
 !coalesced hydrometeor density calculation approach, similar to "col3"
 !would allow for this. Otherwise, snow/aggregates retaining water could
 !allow for too high of a density for snow/aggregates.
 if(lcat.eq.4 .or. lcat.eq.5) rfinlz = min(rcoal,coalliq+umcld)
 if(lcat.eq.6 .or. lcat.eq.7) rfinlz = min(rcoal,coalliq)
 if(lcat.eq.6 .and. tair(k).gt.273.15) rfinlz=0.

 xtoz = min(umcld,rfinlz)
 xtoy = umcld - xtoz
 ytoz = rfinlz - xtoz
 ccat=mcatc(lcat)

 if(imbudget >= 1) xrimecldt(k) = xrimecldt(k) + umcld * budget_scalet
 if(imbudget >= 2) then
    if(lcat.eq.4) xrimecldsnowt(k) = xrimecldsnowt(k) + umcld * budget_scalet
    if(lcat.eq.5) xrimecldaggrt(k) = xrimecldaggrt(k) + umcld * budget_scalet
    if(lcat.eq.6) xrimecldgraut(k) = xrimecldgraut(k) + umcld * budget_scalet
    if(lcat.eq.7) xrimecldhailt(k) = xrimecldhailt(k) + umcld * budget_scalet
 endif

 rxfer(k,jcat,ccat)  =  rxfer(k,jcat,ccat) + xtoz
 rxfer(k,jcat,lcat)  =  rxfer(k,jcat,lcat) + xtoy
 if(lcat.ne.ccat) rxfer(k,lcat,ccat) = rxfer(k,lcat,ccat) + ytoz

 qrxfer(k,jcat,ccat) = qrxfer(k,jcat,ccat) + qx(k,jcat) * xtoz
 qrxfer(k,jcat,lcat) = qrxfer(k,jcat,lcat) + qx(k,jcat) * xtoy
 if(lcat.ne.ccat) &
   qrxfer(k,lcat,ccat) = qrxfer(k,lcat,ccat) + qx(k,lcat) * ytoz

 enxfer(k,jcat,jcat) = enxfer(k,jcat,jcat) + min(uncld,cx(k,mx))
 if(lcat.ne.ccat) enxfer(k,lcat,ccat) = enxfer(k,lcat,ccat) &
   + ytoz * min(uncld,cx(k,lcat)) / max(1.0e-20,rx(k,lcat))

endif !if cloud mean mass is greater than min threshold
endif !if cloud mixing ratio greater than min threshold
enddo !loop of vertical levels

return
END SUBROUTINE auto_accret_ice

!##############################################################################
Subroutine effxy (m1,k1,k2)

use micphys

implicit none

integer :: m1,k,ncall7
integer, dimension(11) :: k1,k2
data ncall7/0/
save

! 1 = rp,rs,ra,rg,rh
if (ncall7 .eq. 0 .and. jnmb(2) .ge. 1 .and. jnmb(3) .ge. 1) then
   ncall7 = 7
   do k = 2,m1-1
      eff(k,1) = 1.0
   enddo
endif

! 2 = cs,ca
! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
! close to curve for 404 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(1),k2(1)
      if (emb(k,1) .gt. 9.e-13) then
         eff(k,2) = min(1.,30. * (emb(k,1) - 9.e-13) ** .15)
      else
         eff(k,2) = 0.
      endif
   enddo
endif

! 3 = rr
! Negative efficiencies used to handle rain drop breakup
! Mass(.113e-6kg) = Mean Diam(0.59mm)
! Mass(.158e-5kg) = Mean Diam(1.44mm)
! If Dm < 0.6mm, limit E <=  1.0
! If Dm > 0.6mm, limit E >= -5.0
if (jnmb(2) .ge. 1) then
   do k = k1(2),k2(2)
    if (rx(k,2) .ge. rxmin) then
      eff(k,3) = min( 1.0, max( -5.0 , 2.-exp(0.1326e7*(emb(k,2)-0.113e-6)) ) )
    endif
   enddo
endif

! 4 = pp,ps,pa
if (jnmb(5) .ge. 1) then
   do k = k1(3),k2(3)
      if (abs(tx(k,3)+14.) .le. 2.) then
         eff(k,4) = 1.4
      else
         eff(k,4) = min(0.2,10. ** (0.035 * tx(k,3) - 0.7))
      endif
   enddo

! 5 = ss,sa
   do k = k1(4),k2(4)
      if (abs(tx(k,4)+14.) .le. 2.) then
         eff(k,5) = 1.4
      else
         eff(k,5) = min(0.2,10. ** (0.035 * tx(k,4) - 0.7))
      endif
   enddo

! 6 = aa
   do k = k1(5),k2(5)
    if (rx(k,5) .ge. rxmin) then
      if (abs(tx(k,5)+14.) .le. 2.) then
         eff(k,6) = 1.4
      elseif (tx(k,5) .ge. -1.) then
         eff(k,6) = 1.
      else
         eff(k,6) = min(0.2,10. ** (0.035 * tx(k,5) - 0.7))
      endif
    endif
   enddo
endif

! 7 = pg,sg,ag,gg,gh
if (jnmb(6) .ge. 1) then
   do k = k1(6),k2(6)
      if (qr(k,6) .gt. 0.) then
         eff(k,7) = 1.0
      else
         eff(k,7) = min(0.2,10. ** (0.035 * tx(k,6) - 0.7))
      endif
   enddo
endif

! 8 = ph,sh,ah,gh
if (jnmb(7) .ge. 1) then
   do k = k1(7),k2(7)
    if (rx(k,7) .ge. rxmin) then
      if (qr(k,7) .gt. 0.) then
         eff(k,8) = 1.0
      else
         eff(k,8) = min(0.2,10. ** (0.035 * tx(k,7) - 0.7))
      endif
    endif
   enddo
endif

! 9 = cg,ch
! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
! close to curves for 142 and 305 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(1),k2(1)
      if (emb(k,1) .gt. 3.4e-14) then
         eff(k,9) = min(1.,1426. * (emb(k,1) - 3.4e-14) ** .28)
      else
         eff(k,9) = 0.
      endif
   enddo
endif

! 10 = hh (trial)
if (jnmb(7) .ge. 1) then
   do k = k1(7),k2(7)
      eff(k,10) = max(0.,.1 + .005 * tx(k,7))
   enddo
endif

! 11 = ds,da
! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
! close to curve for 404 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(8),k2(8)
      if (emb(k,8) .gt. 9.e-13) then
         eff(k,11) = min(1.,30. * (emb(k,8) - 9.e-13) ** .15)
      else
         eff(k,11) = 0.
      endif
   enddo
endif

! 12 = dg,dh
! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
! close to curves for 142 and 305 microns.  Replace with auto_accret eventually.
if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
   do k = k1(8),k2(8)
      if (emb(k,8) .gt. 3.4e-14) then
         eff(k,12) = min(1.,1426. * (emb(k,8) - 3.4e-14) ** .28)
      else
         eff(k,12) = 0.
      endif
   enddo
endif

return
END SUBROUTINE effxy

!##############################################################################
Subroutine cols (mx,mc1,k1,k2)

use micphys

implicit none

integer :: ipc,mx,mc1,k1,k2,k
real :: colnum,tabval

do k = k1,k2
if(rx(k,mx) .ge. rxmin) then
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

   tabval  &
   = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
   + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
   + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colnum = colfacc(k) * eff(k,mc1) * cx(k,mx) ** 2 * 10. ** (-tabval)
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colnum)
endif
enddo

return
END SUBROUTINE cols

!##############################################################################
Subroutine col3344 (mx,mz,mc1,k1,k2)

use micphys

implicit none

integer :: mx,mz,mc1,k1,k2,k,ip,ipc
real :: c1,tabvalx,colamt,tabvaln,colnum

do k = k1,k2
if(rx(k,mx) .ge. rxmin) then
   ip = ipairr(jhcat(k,mx),jhcat(k,mx))
   ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
   c1 = eff(k,mc1) * cx(k,mx) ** 2

   tabvalx  &
    = wct1(k,mx) ** 2               * coltabr(ict1(k,mx),ict1(k,mx),ip)  &
    + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ip)  &
    + wct2(k,mx) ** 2               * coltabr(ict2(k,mx),ict2(k,mx),ip)

   colamt = min(rx(k,mx),colfacr2(k) * c1 * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + colamt
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colamt * qx(k,mx)

   if(imbudget >= 1) xaggregatet(k) = xaggregatet(k) + colamt*budget_scalet
   if(imbudget >= 2) then
    if(mx.eq.3) xaggrselfprist(k) = xaggrselfprist(k) + colamt*budget_scalet
    if(mx.eq.4) xaggrselfsnowt(k) = xaggrselfsnowt(k) + colamt*budget_scalet
   endif

   if (jnmb(mz) >= 5) then

   tabvaln  &
   = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
   + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
   + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

   colnum = min(0.5 * cx(k,mx),colfacc2(k) * c1 * 10. ** (-tabvaln))
   enxfer(k,mx,mz) = enxfer(k,mx,mz) + colnum
   enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnum

   endif
endif
enddo

return
END SUBROUTINE col3344

!##############################################################################
Subroutine col3443 (mx,my,mz,k1,k2)

use micphys

implicit none

integer :: mx,my,mz,k1,k2,k,jhcatx,jhcaty,ipxy,ipyx,ipc
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum

do k = k1,k2
if(rx(k,mx) .ge. rxmin .and. rx(k,my) .ge. rxmin) then
   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)
   c1 = eff(k,4) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

   rxfer(k,my,mz) = rxfer(k,my,mz) + rcy
   qrxfer(k,my,mz) = qrxfer(k,my,mz) + rcy * qx(k,my)

   if(imbudget >= 1) xaggregatet(k) = xaggregatet(k) + (rcx + rcy)*budget_scalet
   if(imbudget >= 2) xaggrprissnowt(k) = xaggrprissnowt(k) + (rcx + rcy)*budget_scalet

   tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)
   colnum = c1 * colfacc(k) * 10. ** (-tabvaln)

   if (cx(k,mx) .gt. cx(k,my)) then
      enxfer(k,my,mz) = enxfer(k,my,mz) + min(cx(k,my),colnum)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(cx(k,mx),colnum)
   else
      enxfer(k,mx,mz) = enxfer(k,mx,mz) + min(cx(k,mx),colnum)
      enxfer(k,my,my) = enxfer(k,my,my) + min(cx(k,my),colnum)
   endif

! also loss for aerosol

endif
enddo

return
END SUBROUTINE col3443

!##############################################################################
Subroutine col1 (mx,my,mz,mc4,k1,k2)

use micphys

implicit none

integer :: mx,my,mz,mc4,k1,k2,k,ipxy,ipc
real :: c1,tabvalx,rcx,tabvaln,colnum

do k = k1,k2
if(rx(k,mx) .ge. rxmin .and. rx(k,my) .ge. rxmin) then
   ipxy = ipairr(jhcat(k,mx),jhcat(k,my))
   ipc  = ipairc(jhcat(k,mx),jhcat(k,my))
   c1 = eff(k,mc4) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))
   rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

   if (jnmb(mx) >= 5) then
      tabvaln  &
        = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
        + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
        + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
        + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))

! also loss for aerosol

   endif

endif
enddo

return
END SUBROUTINE col1

!##############################################################################
Subroutine col2 (m1,mx,my,mz,mc2,k1,k2,dn0,dtlt)

use rconstants
use micphys

implicit none

integer :: m1,mx,my,mz,mc2,k1,k2,k,jhcatx,jhcaty,ipxy,ipyx,ipc,it,mxx
real :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal  &
       ,qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area,cn13,cn24  &
       ,sip,rsip,qrsip,rfinlz,xtoz,dtlt

real, dimension(m1) :: dn0

real, dimension(16) ::  alpha,beta
!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
data alpha /00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1.,00./
data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3,00./

do k = k1,k2
if(rx(k,mx) .ge. rxmin .and. rx(k,my) .ge. rxmin) then
   jhcatx = jhcat(k,mx)
   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcatx,jhcaty)
   ipyx = ipairr(jhcaty,jhcatx)
   ipc  = ipairc(jhcatx,jhcaty)
   c2 = cx(k,mx) * cx(k,my)
   c1 = eff(k,mc2) * c2

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   if (jnmb(mx) >= 5 .or. jnmb(my) >= 5) then
      tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)
      colnum0 = c2 * colfacc(k) * 10. ** (-tabvaln)
      colnum = colnum0 * eff(k,mc2)
   else
      colnum0 = 0.
      colnum = 0.
   endif

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,mx)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy

   !Saleeby (3-15-06) Lowering the threshold to prevent over-freezing
   qcoal = qrcoal / max(1.0e-20,rcoal)

   CALL qtc (qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

! Secondary ice production: cn24 is the number fraction of collected cloud
! droplets larger than 24 microns and is obtained from an incomplete gamma
! func table. cn13 is the fraction of collected cloud droplets
! smaller than 13 microns. "area" is cross section area of collecting ice
! per m^3 of atmospheric volume.

! Saleeby(6/3/02): Hallett-Mossop is done for both cloud droplet modes, though
! contribution from the large droplet mode is minimal compared to the small
! droplet mode. Ice splintering is only done if number concentration is
! prognostic for at least one of the two hydromet species involved. This is
! specified above in the calculations for "colnum".

   if (tcoal.gt.-8.0 .and. tcoal.lt.-3.0 .and. (jnmb(mx)>=5.or.jnmb(my)>=5)) then
      mxx=0
      if(mx==1) mxx=1 !uses gamsip(1,it) for cloud1 set in mic_init.f90
      if(mx==8) mxx=2 !uses gamsip(2,it) for cloud2 set in mic_init.f90
      if(mxx<1.or.mxx>2)stop 'problem with 2ndary ice production'
      area = cx(k,my) * dn0(k) * sipfac(jhcaty) * emb(k,my)  &
         ** (2.*pwmasi(jhcaty))
      it = nint(emb(k,mx) / emb1(mx) * 5000.)
      cn13 = colnum * gamsip13(mxx,it) / (area * dtlt)
      cn24 = min(cx(k,mx)*dn0(k),colnum0) * gamsip24(mxx,it)
      sip = 9.1e-10 * cn24 * cn13 ** .93
      if (tcoal .lt. -5.) then
         sip = 0.33333 * (tcoal + 8.) * sip
      else
         sip = -0.5 * (tcoal + 3.) * sip
      endif
      rsip = sip * emb0(3) * dn0i(k)
      qrsip = qcoal * rsip
      rcoal = rcoal - rsip
      qrcoal = qrcoal - qrsip
      enxfer(k,mx,3) = enxfer(k,mx,3) + sip
      rxfer(k,mx,3) = rxfer(k,mx,3) + rsip
      qrxfer(k,mx,3) = qrxfer(k,mx,3) + qrsip
   endif

! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
! fracliq may be a little larger than fracx due to collected
! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
! or 1.2, or need ALPHA itself to be at least 1.0.

   rfinlz = min(rcoal,  &
      alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

   xtoz = min(rcx,rfinlz)

   if(imbudget >= 1) xrimecldt(k) = xrimecldt(k) + rcx * budget_scalet
   if(imbudget >= 2) then
    if(my.eq.4) xrimecldsnowt(k) = xrimecldsnowt(k) + rcx * budget_scalet
    if(my.eq.5) xrimecldaggrt(k) = xrimecldaggrt(k) + rcx * budget_scalet
    if(my.eq.6) xrimecldgraut(k) = xrimecldgraut(k) + rcx * budget_scalet
    if(my.eq.7) xrimecldhailt(k) = xrimecldhailt(k) + rcx * budget_scalet
   endif

   rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
   rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
   if (my .ne. mz) rxfer(k,my,mz) = rxfer(k,my,mz)  &
      + rfinlz - xtoz

   qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
   qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
   if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
      + qx(k,my) * (rfinlz - xtoz)

   if (jnmb(mx) >= 5 .or. jnmb(my) >= 5) then
     enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))
     if (my .ne. mz) enxfer(k,my,mz) = enxfer(k,my,mz)  &
        + (rfinlz - xtoz) * min(colnum,cx(k,my)) / max(1.0e-20,rcy)
   endif

! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y

endif
enddo

return
END SUBROUTINE col2

!##############################################################################
Subroutine col3 (m1,mx,my,k1,k2,dn0)

!Rain-ice collection routine that classifies destination category
!depending on density of resulting coalesced particles for r-s, r-a and r-g
!collection based on Milbradnt and Yau 2005 part II equations 39 and 40.
!Rain-pristine ice collisions go to hail using Ferrier 1994 formulation that
!allows for larger drops to collect more than one ice crystal.
! * r-s and r-a collisions go to either s,a,g or h category - compute coalesced
!   particle density to determine destination cat
! * r-g collisions go to either graupel or hail depending on resulting
!    density (denz) of coalesced particles [after Milbrandt and Yau 2005b]
! * r-p collisions go to hail using Ferrier 1994 formulation that allows
!   for larger drops to collect more than one ice crystal; special routine
!   'F94' called for r-p collisions

use micphys

implicit none

integer :: m1,mx,my,mz,k1,k2,k,ipxy,ipyx,ipc,jhcaty
real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum,colnumx,colnumy,coalnum  &
       ,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq,coalliq,coalice,xtoz  &
       ,rfinlz,tcoal,cfinlz
real, dimension(m1) :: dn0
real, dimension(16) :: alpha,beta
real :: cfinlzMY05,tmp1,tmp2,tmp3,denr,dens,deng,denz,denh,pwcoef,dmzz &
  ,dmeanr,dmeany,dengh,densg,ansn,ansr

!            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
data alpha /00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,00./
data beta  /00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2.,00./

!--Set local constants
denr=1000. !rain density [kg/m^3]
deng=300.  !graupel density [kg/m^3]
denh=900.  !hail density [kg/m^3]
dengh=600. !=0.5*(deng+denh) [kg/m^3]

!Note: mx=2 (rain), my=lcat, mz=7 (hail)
do k = k1,k2

 mz=7  !<--default destination category, needs to be reset for each level k
 if(rx(k,mx) .ge. rxmin .and. rx(k,my) .ge. rxmin) then

   jhcaty = jhcat(k,my)
   ipxy = ipairr(jhcat(k,mx),jhcaty)
   ipyx = ipairr(jhcaty,jhcat(k,mx))
   ipc  = ipairc(jhcat(k,mx),jhcaty)
   c1 = eff(k,1) * cx(k,mx) * cx(k,my)

   tabvalx  &
     = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
     + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
     + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
     + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

   rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

   tabvaly  &
     = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
     + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
     + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
     + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

   rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

   if (jnmb(mx) >= 5) then
      tabvaln  &
       = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
       + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
       + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
       + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
      colnumx = min(cx(k,mx),colnum)
      colnumy = min(cx(k,my),colnum)
      coalnum = min(colnumx,colnumy)
   else
      colnumx = 0.
      colnumy = 0.
      coalnum = 0.
   endif

   rcoal = rcx + rcy
   qrcx = rcx * qx(k,mx)
   qrcy = rcy * qx(k,my)
   qrcoal = qrcx + qrcy
   qcoal = qrcoal / max(1.0e-20,rcoal)

   CALL qtc (qcoal,tcoal,fracliq)

   coalliq = rcoal * fracliq
   coalice = rcoal - coalliq

   !Air temperature threshold for rain-ice coll (No new ice at T > 0 degC)
   if (fracliq.ge.0.99 .or. tair(k).gt.273.15) then

      rxfer(k,my,mx) = rxfer(k,my,mx) + rcy
      qrxfer(k,my,mx) = qrxfer(k,my,mx) + qrcy
      if (jnmb(mx) >= 5)  &
         enxfer(k,my,my) = enxfer(k,my,my) + colnumy

      if(imbudget >= 1) xice2raint(k) = xice2raint(k) + rcy * budget_scalet

   else

      rfinlz = min(rcoal,  &
         alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

      !For pristine ice, snow, aggregates, graupel: 
      !'xtoz' is ALWAYS equal to 'rcx' because 'rfinlz' is ALWAYS 
      !at least equal to 'rcx' or 'rcx+rcy'

      xtoz = min(rcx,rfinlz)

      !Partitioning for destination categories based
      !on Milbrandt and Yau 2005pII 3-component freezing (Eqns 39 & 40)
      if((my.ge.4).and.(my.le.6))then

       dmeanr=(emb(k,mx)/cfmas(mx))**(1./pwmas(mx)) !rain mean mass D [m]
       dmeany=(emb(k,my)/cfmas(my))**(1./pwmas(my)) !snow/aggr/grpl mean mass D
       dmzz=max(dmeanr,dmeany)        !mean mass diam of coalesc particles [m]
       dens=deng                      !default ice cat rho is that for grpl

       if((my.eq.4).or.(my.eq.5))then
        !coalesced particles will be hail based on pre-determined (offline)
        !calculations of resulting densities (lrg rain collects smaller ice)
        if(dmeanr.ge.dmeany)then
         denz=denh
         goto 33
        !compute densities of mean mass diams of s,a to determine dest cat
        else
         pwcoef=pwmas(my)-3.
         !snow/aggr density of mean mass diam [kg/m^3]
         dens=(6./3.14159)*cfmas(my)*(dmeany**pwcoef) 
         densg=0.5*(dens+deng)
        endif
       endif

       tmp1=dmzz*dmzz*dmzz
       tmp2=dmeanr*dmeanr*dmeanr
       tmp3=dmeany*dmeany*dmeany
       !rho_z from Eqn 39 (rho of coalesced particles)
       denz=(denr*tmp2+dens*tmp3)/tmp1

       !Determine destination category based on denz (rho_z)
       ! NOTE destination category is hail (mz=7) by default
       if((my.eq.4).or.(my.eq.5))then
        if(denz.le.densg)then
         mz=my                 !<--dest category is snow or aggr
        elseif((denz.gt.densg).and.(denz.le.dengh))then
         mz=6                  !<--dest category is graupel
        endif
       elseif((my.eq.6).and.(denz.le.dengh))then
        mz=6                   !<--dest category is graupel
       endif

       !Set dest category density
       denz=denh  !default dest. cat. density is that of hail
       if(mz.lt.6)then
        denz=dens
       elseif(mz.eq.6)then
        denz=deng
       endif

      elseif(my.eq.3)then !Use Ferrier 1994 formulation for r-p collisions
                          !to limit mass and # of new hailstones (frozen drops)
       dmeanr=(emb(k,mx)/cfmas(mx))**(1./3.) !mean mass diam of rain [m]
       dmeany=(emb(k,my)/cfmas(my))**(1./pwmas(my)) !pris ice mean mass diam [m]
       denz=denh    !coalesced particles have density of hail
       dmzz=dmeanr  !coalesced particles have mean mass diam of rain

       !diam threshholds for hail formation
       if((dmeanr.ge.8.e-4).and.(dmeany.ge.40.e-6))then 
        !ansn unitless, ansr [kg]
        CALL F94 (rx(k,2),cx(k,2),cx(k,3),dn0(k),ansn,ansr)
        ansn=ansn*cx(k,2) !#/conc of new hailstones [1/kg]
        ansr=ansr*cx(k,2) !mixing ratio of new hailstones [kg/kg]
        cfinlz=coalnum*rfinlz/max(rcoal,1.0e-20) !regular 2M value for # of new hail
        if(ansr.lt.rcx)then !recompute these variables
         rcx=ansr
         rcoal=rcx+rcy
         qrcx=rcx*qx(k,mx)
         xtoz=rcx
         rfinlz=rcoal
        endif
        if(ansn.lt.cfinlz)then
         colnumx=min(ansn,cx(k,2))
         coalnum=min(colnumx,colnumy)
        endif
       !rain and/or pris ice mean mass diams too small to form hail
       else
        xtoz=0.
        rfinlz=0.
        rcoal=0.
       endif !end of dmeanr >= 8.e-4 & dmeany >= 40e-6

      endif !end of my >= 4 and my < 7
      33 continue

      if(imbudget >= 1) xrain2icet(k) = xrain2icet(k) + rcx * budget_scalet
      if(imbudget >= 2) then
        if(my==3) xrain2prt(k) = xrain2prt(k) + rcx * budget_scalet
        if(my==4) xrain2snt(k) = xrain2snt(k) + rcx * budget_scalet
        if(my==5) xrain2agt(k) = xrain2agt(k) + rcx * budget_scalet
        if(my==6) xrain2grt(k) = xrain2grt(k) + rcx * budget_scalet
        if(my==7) xrain2hat(k) = xrain2hat(k) + rcx * budget_scalet
      endif

      rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz !from rain to cat z
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz

      rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
      qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)

      if (my .ne. mz) then !from cat y to cat z
       rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz
       qrxfer(k,my,mz) = qrxfer(k,my,mz) + qx(k,my) * (rfinlz - xtoz)
       tmp1=rcoal*(6./3.14159)
       tmp2=denz*dmzz*dmzz*dmzz
       cfinlzMY05=tmp1/tmp2  !<--Eqn 40 from Milbrandt and Yau 2005 part II
      endif                  !   # conc of newly formed cat z particles [g or h]

      if (jnmb(mx) >= 5) then
         cfinlz = coalnum * rfinlz / max(rcoal,1.0e-20)
         if (my .eq. mz) then
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
         elseif (colnumy .ge. colnumx) then
            enxfer(k,mx,mz) = enxfer(k,mx,mz) + min(cfinlzMY05,cfinlz)
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx - cfinlz &
                            + max(0.,cfinlz-cfinlzMY05)
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy
         else
            enxfer(k,my,mz) = enxfer(k,my,mz) + min(cfinlzMY05,cfinlz)
            enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy - cfinlz &
                            + max(0.,cfinlz-cfinlzMY05)
         endif
      endif !end of jnmb(mx) >= 5

   endif !end of 'fracliq >= 0.99' section
 endif   !end of rx(k,mx), rx(k,my) >= rxmin if-statement
enddo

return
END SUBROUTINE col3

!##############################################################################
Subroutine F94 (rxa,cxx,cxy,denk,ansn,ansr)

! call F94(rx(k,2),cx(k,2),cx(k,3),dn0(k),ansn,ansr)
! Compute mass and number of frozen drops for r-p collisions
! Based on Eqn 4.19 for mass and Eqn 4.23 for number conc in Ferrier 1994.
! This allows for larger raindrops to collect more than 1 ice crystal in
! one timestep dtlt when large numbers of crystals are present. Assumes
! that crystal fall speed is zero and rain-pris collection efficiencies are 1.
! This routine uses [1/m^3] and [kg/m^3] units

use mem_grid, only:dtlt
use micphys

implicit none

integer :: nbins,i,j
real :: rxa,cxx,cxy,dt,ipi4,pi314,mind,maxd,m
real :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,summ,sumn
real :: chdiam,meand,dmbodn,ansn,ansr,denk
real,dimension(95) :: fgam
real,dimension(94) :: Pri,Nri
real, external :: gammln

!--local constants
nbins=94
ipi4=3.1415927/4.
sumn=0.
summ=0.
ansn=0.
ansr=0.

!--compute characteristic diameter for rain (lcat=2)
dmbodn=(exp(gammln(gnu(2)+pwmas(2))-gammln(gnu(2))))**pwmasi(2)
tmp1=rxa/(max(cxx,1.0e-12))  !mean mass
meand=(tmp1/cfmas(2))**pwmasi(2)       !mean mass diam
chdiam=meand/dmbodn         !char. diam

tmp2=1./chdiam
tmp4=gnu(2)-1.
tmp5=exp(gammln(gnu(2)))
tmp5=1./tmp5

do i=1,nbins !numerically integrate Eqs 4.19 & 4.23 over rain diameters
!-- Compute Nri, Pri (see F94) Nri: # of crystals collected by 
!-- a drop with diam 'radF94' in one dtlt, assumed Eff=1
 Nri(i)=ipi4*cxy*denk*(radF94(i)*radF94(i))*avtF94(i)*dtlt !Eqn 4.21 in F94
 Nri(i)=Nri(i)/sqrt(denk)  !<--AML 2/4/11 account for densty effcts on Vt
 Pri(i)=min(1.,Nri(i))  !Probabilistic freezing of rain [Eqn 4.20 in F94]
 tmp1=radF94(i)*tmp2
 fgam(i)=tmp2*exp(-tmp1)*(tmp1**tmp4)*tmp5 !gamma pdf values
 if((fgam(i).lt.1.0e-35).or.(fgam(i).gt.1.0e+35))THEN
  fgam(i)=0.0
 endif
 tmp6=rdbF94(i+1)-rdbF94(i)  !size bin width
 sumn=sumn+Pri(i)*fgam(i)*tmp6  !eqn 4.23 in F94
 summ=summ+ramF94(i)*Pri(i)*fgam(i)*tmp6  !eqn 4.19 in F94
enddo

ansn=sumn !number of drops frozen
ansr=summ !change in rain mass owing to freezing of rain drops

return
END SUBROUTINE F94

!##############################################################################
Subroutine colxfers (m1,k1,k2,rloss,enloss)

use micphys

implicit none

integer :: m1,k,lcat,kd1,kd2,jcat
integer, dimension(11) :: k1,k2
real, dimension(m1) :: rloss,enloss
real :: cldnumratio

!  All rxfer values are nonnegative.
do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do k = kd1,kd2
         rloss(k) = 0.
         enloss(k) = 0.
      enddo

      do jcat = 1,8
! change this to include enxfer of the same categories
         if (jnmb(jcat) .ge. 1) then
            if (lcat .ne. jcat) then
               do k = kd1,kd2
                  rloss(k) = rloss(k) + rxfer(k,lcat,jcat)
               enddo
            endif
            do k = kd1,kd2
               enloss(k) = enloss(k) + enxfer(k,lcat,jcat)
            enddo
         endif
      enddo

      do k = kd1,kd2
         rloss(k) = min(1.,rx(k,lcat) / max(1.0e-20,rloss(k)))
         enloss(k) = min(1.,cx(k,lcat) / max(1.0e-10,enloss(k)))
      enddo

      do jcat = 1,8
         if (jnmb(jcat) .ge. 1) then
            if (lcat .ne. jcat) then
               do k = kd1,kd2
                  rxfer(k,lcat,jcat) = rxfer(k,lcat,jcat)*rloss(k)
                  qrxfer(k,lcat,jcat)=qrxfer(k,lcat,jcat)*rloss(k)
               enddo
            endif
            do k = kd1,kd2
               enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat)*enloss(k)
            enddo
         endif
      enddo
   endif
enddo

do lcat = 1,8

   if (jnmb(lcat) .ge. 1) then

      kd1 = k1(lcat)
      kd2 = k2(lcat)

      do jcat = 1,8 

         !TRANSFER MIXING RATIO AND NUMBER AND CCN MASS BETWEEN CATEGORIES
         if (jnmb(jcat) .ge. 1 .and. lcat .ne. jcat) then
            do k = kd1,kd2

               !Double check prevent transfer of more rx than exists
               if(rxfer(k,lcat,jcat)>rx(k,lcat)) then
                 rxfer(k,lcat,jcat)=rx(k,lcat)
               endif

               !Aerosol and solubility tracking
               !Transfer aerosol-in-hydrometeor masses between categories
               if(iccnlev>=2 .and. rx(k,lcat)>0.0) then
                  rxferratio = rxfer(k,lcat,jcat)/rx(k,lcat)
                  ccnmass  = cnmhx(k,lcat) * rxferratio
                  cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
                  cnmhx(k,jcat) = cnmhx(k,jcat) + ccnmass
                  if(itrkepsilon==1) then
                   scnmass  = snmhx(k,lcat) * rxferratio
                   snmhx(k,lcat) = snmhx(k,lcat) - scnmass
                   snmhx(k,jcat) = snmhx(k,jcat) + scnmass
                  endif
                  if(itrkdust==1)then
                   dcnmass  = dnmhx(k,lcat) * rxferratio
                   dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
                   dnmhx(k,jcat) = dnmhx(k,jcat) + dcnmass
                  endif
                  if(itrkdustifn==1)then
                   dinmass  = dinhx(k,lcat) * rxferratio
                   dinhx(k,lcat) = dinhx(k,lcat) - dinmass
                   dinhx(k,jcat) = dinhx(k,jcat) + dinmass
                  endif
               endif
               !Transfer immersion freezing nuclei between cloud, driz, rain
               if(iifn==3 .and. iccnlev>=1 .and. cx(k,lcat)>0.0) then
                 if(lcat==1.or.lcat==8.or.lcat==2) then
                   enxferratio = min(1.0,max(0.0,enxfer(k,lcat,jcat)/cx(k,lcat)))
                   ccnnum = immerhx(k,lcat) * enxferratio
                   immerhx(k,lcat) = immerhx(k,lcat) - ccnnum
                   if(jcat==8.or.jcat==2) &
                     immerhx(k,jcat) = immerhx(k,jcat) + ccnnum
                 endif
               endif
               !Update mixing ratio, number, and heat storage
               rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
               rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)
               qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
               qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)
               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
               cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
               if(rx(k,lcat)<=0.0 .or. cx(k,lcat)<=0.0) then
                 rx(k,lcat)=0.0
                 cx(k,lcat)=0.0
                 qx(k,lcat)=0.0
                 qr(k,lcat)=0.0
               endif

            enddo
         endif
      enddo

      !Section to handle self-collection number change
      if (jnmb(lcat) >= 5) then
         do k = kd1,kd2
            cxtemp = cx(k,lcat)
            cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,lcat)
            if(cx(k,lcat)<0.0) cx(k,lcat)=0.0
            !Adjust immersion freezing nuclei number
            if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==8.or.lcat==2)) then
              cldnumratio = min(1.0,max(0.0,cx(k,lcat)/max(1.0e-20,cxtemp)))
              immerhx(k,lcat) = max(0.,immerhx(k,lcat)*cldnumratio)
            endif
         enddo
      endif

   endif
enddo

return
END SUBROUTINE colxfers
