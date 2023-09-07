!##############################################################################
!  IMPORTANT ISSUE: k loop limits for the jnmb == 5 sections
!  need to consider collection efficiencies for different habits?
!  collection efficiency for hail too high.  big hail should not
!  coallesce.
!##############################################################################
Subroutine each_call (m1,dtlt)

use rconstants
use micphys

implicit none

integer :: m1,lcat,k,lhcat
real :: dtlt
integer, dimension(8) :: lcat0
data lcat0 /1,2,3,4,5,6,7,16/ ! lcat corressponding to lhcat

! Initialize constants for vapor diffusion and, for fixed diameter cases, emb.

colf = .785 * dtlt
pi4dt = pi4 * dtlt
sl(1) = alvl
sl(2) = alvi
sc(1) = 4186.
sc(2) = 2093.
sj(1) = 0
sj(2) = 1
sj(3) = 0
sj(4) = 0
sj(5) = 0
sj(6) = 1
sj(7) = 1
sj(8) = 0
sk(1) = alli
sk(2) = 0.

do lcat = 1,8
   lhcat = lcat0(lcat)
   if (jnmb(lcat) == 2) then
      do k = 2,m1-1
         emb(k,lcat) = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)
      enddo
   endif
   do k = 2,m1-1
      jhcat(k,lcat) = lhcat
   enddo
enddo

do k = 2,m1-1
   sh(k,1) = 0.
   sh(k,2) = 1.
   sh(k,6) = 1.
   sh(k,7) = 1.
   sh(k,8) = 0.

   sm(k,1) = 1.
   sm(k,2) = 1.
   sm(k,8) = 1.
enddo

return
END SUBROUTINE each_call

!##############################################################################
Subroutine range_check (m1,k1,k2,k3,i,j,frq,ngr,dtlt,time,micro)

use mem_micro
use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

type (micro_vars) :: micro

integer :: m1,i,j,k,lcat,l,jcat,ngr
integer, dimension(11) :: k1,k2,k3
real :: frq,time,dtlt

!Zero out microphysics scratch arrays for the present i,j column
do lcat = 1,ncat
   do k = 2,m1-1
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qr(k,lcat) = 0.
      qx(k,lcat) = 0.
      vap(k,lcat) = 0.
      tx(k,lcat) = 0.
   enddo

   if (jnmb(lcat) >= 3) then
      do k = 2,m1-1
         emb(k,lcat) = 0.
      enddo
   endif

   do jcat = 1,ncat
      do k = 2,m1-1
         rxfer(k,lcat,jcat) = 0.
         qrxfer(k,lcat,jcat) = 0.
         enxfer(k,lcat,jcat) = 0.
      enddo
   enddo
enddo

!Zero out aerosol scratch arrays
do acat = 1,aerocat
    do k = 1,m1
      aerocon(k,acat) = 0.0
      aeromas(k,acat) = 0.0
    enddo
enddo

!Zero out aerosol tracking mass scratch arrays
if(iccnlev>=2) then
 do lcat = 1,ncat
  pcpraerox(lcat) = 0.
  accpaerox(lcat) = 0.
  pcprdustx(lcat) = 0.
  accpdustx(lcat) = 0.
  do k = 1,m1
    cnmhx(k,lcat) = 0.
    snmhx(k,lcat) = 0.
    dnmhx(k,lcat) = 0.
    dinhx(k,lcat) = 0.
  enddo
 enddo
 do k = 1,m1
   regenmas(k,1)=0.0
   regenmas(k,2)=0.0
 enddo
endif

!Zero out IFN aerosol scratch arrays
if (iifn==3) then
 do k = 1,m1  
   nifn(k) = 0.
 enddo
 if(iccnlev>=1)then
   ifnnucx(k) = 0.
   do lcat = 1,ncat
     immerhx(k,lcat) = 0.
   enddo
 endif
endif

!Zero out top and bottom of hydrometeor layers
do l = 1,8
   k1(l) = 2
   k2(l) = 1
enddo

!************************************************************************
! FILL SCRATCH ARRAYS FOR TEMPORARY MICROPHYSICS VARIABLE IN COLUMNS
!************************************************************************

!Fill scratch arrays for aerosol modes
do k = 2,m1-1
  if (jnmb(3)>=5 .and. (iifn==1.or.iifn==2)) then
    cifnx(k) = micro%cifnp(k,i,j)
  endif
  if (iaerosol > 0) then
    aerocon(k,1)   = micro%cn1np(k,i,j)
    aeromas(k,1)   = micro%cn1mp(k,i,j)
    aerocon(k,2)   = micro%cn2np(k,i,j)
    aeromas(k,2)   = micro%cn2mp(k,i,j)
  endif
  if (idust > 0) then
    aerocon(k,3) = micro%md1np(k,i,j)
    aeromas(k,3) = micro%md1mp(k,i,j)
    aerocon(k,4) = micro%md2np(k,i,j)
    aeromas(k,4) = micro%md2mp(k,i,j)
  endif
  if (isalt > 0) then
    aerocon(k,5) = micro%salt_film_np(k,i,j)
    aeromas(k,5) = micro%salt_film_mp(k,i,j)
    aerocon(k,6) = micro%salt_jet_np(k,i,j)
    aeromas(k,6) = micro%salt_jet_mp(k,i,j)
    aerocon(k,7) = micro%salt_spum_np(k,i,j)
    aeromas(k,7) = micro%salt_spum_mp(k,i,j)
  endif
  if (iabcarb > 0) then
    aerocon(k,8) = micro%abc1np(k,i,j)
    aeromas(k,8) = micro%abc1mp(k,i,j)
    aerocon(k,9) = micro%abc2np(k,i,j)
    aeromas(k,9) = micro%abc2mp(k,i,j)
  endif
enddo

!Aerosol and solubility tracking variables scratch arrays
if (iccnlev>=2) then
 do k = 2,m1-1
   !Regenerated aerosol and solubility
   aerocon(k,aerocat-1) = micro%regen_aero1_np(k,i,j)
   aeromas(k,aerocat-1) = micro%regen_aero1_mp(k,i,j)
   aerocon(k,aerocat) = micro%regen_aero2_np(k,i,j)
   aeromas(k,aerocat) = micro%regen_aero2_mp(k,i,j)
   if(itrkepsilon==1) then
     regenmas(k,1) = micro%resol_aero1_mp(k,i,j)
     regenmas(k,2) = micro%resol_aero2_mp(k,i,j)
   endif
   !Aerosol masses in cloud
   if (jnmb(1) >= 1) then
     if (micro%rcp(k,i,j) >= rxmin) then
      cnmhx(k,1) = micro%cnmcp(k,i,j)
      if(itrkepsilon==1) snmhx(k,1) = micro%snmcp(k,i,j)
      if(itrkdust==1)    dnmhx(k,1) = micro%dnmcp(k,i,j)
      if(itrkdustifn==1) dinhx(k,1) = micro%dincp(k,i,j)
     endif
   endif
   !Aerosol masses in rain
   if (jnmb(2) >= 1) then
     if (micro%rrp(k,i,j) >= rxmin) then
      cnmhx(k,2) = micro%cnmrp(k,i,j)
      if(itrkepsilon==1) snmhx(k,2) = micro%snmrp(k,i,j)
      if(itrkdust==1)    dnmhx(k,2) = micro%dnmrp(k,i,j)
      if(itrkdustifn==1) dinhx(k,2) = micro%dinrp(k,i,j)
     endif
   endif
   !Aerosol masses in pristine ice 
   if (jnmb(3) >= 1) then
     if (micro%rpp(k,i,j) >= rxmin) then
      cnmhx(k,3) = micro%cnmpp(k,i,j)
      if(itrkepsilon==1) snmhx(k,3) = micro%snmpp(k,i,j)
      if(itrkdust==1)    dnmhx(k,3) = micro%dnmpp(k,i,j)
      if(itrkdustifn==1) dinhx(k,3) = micro%dinpp(k,i,j)
     endif
   endif
   !Aerosol masses in snow
   if (jnmb(4) >= 1) then
     if (micro%rsp(k,i,j) >= rxmin) then
      cnmhx(k,4) = micro%cnmsp(k,i,j)
      if(itrkepsilon==1) snmhx(k,4) = micro%snmsp(k,i,j)
      if(itrkdust==1)    dnmhx(k,4) = micro%dnmsp(k,i,j)
      if(itrkdustifn==1) dinhx(k,4) = micro%dinsp(k,i,j)
     endif
   endif
   !Aerosol masses in aggregates
   if (jnmb(5) >= 1) then
     if (micro%rap(k,i,j) >= rxmin) then
      cnmhx(k,5) = micro%cnmap(k,i,j)
      if(itrkepsilon==1) snmhx(k,5) = micro%snmap(k,i,j)
      if(itrkdust==1)    dnmhx(k,5) = micro%dnmap(k,i,j)
      if(itrkdustifn==1) dinhx(k,5) = micro%dinap(k,i,j)
     endif
   endif
   !Aerosol masses in graupel
   if (jnmb(6) >= 1) then
     if (micro%rgp(k,i,j) >= rxmin) then
      cnmhx(k,6) = micro%cnmgp(k,i,j)
      if(itrkepsilon==1) snmhx(k,6) = micro%snmgp(k,i,j)
      if(itrkdust==1)    dnmhx(k,6) = micro%dnmgp(k,i,j)
      if(itrkdustifn==1) dinhx(k,6) = micro%dingp(k,i,j)
     endif
   endif
   !Aerosol masses in hail
   if (jnmb(7) >= 1) then
     if (micro%rhp(k,i,j) >= rxmin) then
      cnmhx(k,7) = micro%cnmhp(k,i,j)
      if(itrkepsilon==1) snmhx(k,7) = micro%snmhp(k,i,j)
      if(itrkdust==1)    dnmhx(k,7) = micro%dnmhp(k,i,j)
      if(itrkdustifn==1) dinhx(k,7) = micro%dinhp(k,i,j)
     endif
   endif
   !Aerosol masses in drizzle
   if (jnmb(8) >= 1) then
     if (micro%rdp(k,i,j) >= rxmin) then
      cnmhx(k,8) = micro%cnmdp(k,i,j)
      if(itrkepsilon==1) snmhx(k,8) = micro%snmdp(k,i,j)
      if(itrkdust==1)    dnmhx(k,8) = micro%dnmdp(k,i,j)
      if(itrkdustifn==1) dinhx(k,8) = micro%dindp(k,i,j)
     endif
   endif
 enddo
endif

!For tracking immersion freezing nuclei
if(iifn==3 .and. iccnlev>=1) then
 do k = 2,m1-1
   if (jnmb(1) >= 5) ifnnucx(k)   = micro%ifnnucp(k,i,j)
   if (jnmb(1) >= 5) immerhx(k,1) = micro%immercp(k,i,j)
   if (jnmb(8) >= 5) immerhx(k,8) = micro%immerdp(k,i,j)
   if (jnmb(2) >= 5) immerhx(k,2) = micro%immerrp(k,i,j)
 enddo
endif

! fill scratch arrays for cloud water

if (jnmb(1) >= 1) then
   do k = 2,m1-1
      if (micro%rcp(k,i,j) >= rxmin) then
         k2(1) = k
         rx(k,1) = micro%rcp(k,i,j)
         if (jnmb(1) >= 5) cx(k,1) = micro%ccp(k,i,j)
      else
         if (k2(1) == 1) k1(1) = k + 1
      endif
   enddo
endif

! fill scratch arrays for rain

if (jnmb(2) >= 1) then
   do k = 2,m1-1
      if (micro%rrp(k,i,j) >= rxmin) then
         k2(2) = k
         rx(k,2) = micro%rrp(k,i,j)
         qx(k,2) = micro%q2(k,i,j)
         qr(k,2) = qx(k,2) * rx(k,2)
         if (jnmb(2) >= 5) cx(k,2) = micro%crp(k,i,j)
      else
         if (k2(2) == 1) k1(2) = k + 1
      endif
   enddo
endif

! fill scratch arrays for pristine ice

if (jnmb(3) >= 1) then
   do k = 2,m1-1
      if (micro%rpp(k,i,j) >= rxmin) then
         k2(3) = k
         rx(k,3) = micro%rpp(k,i,j)
         if (jnmb(3) >= 5) cx(k,3) = micro%cpp(k,i,j)
      else
         if (k2(3) == 1) k1(3) = k + 1
      endif
   enddo
endif

! fill scratch arrays for snow

if (jnmb(4) >= 1) then
   do k = 2,m1-1
      if (micro%rsp(k,i,j) >= rxmin) then
         k2(4) = k
         rx(k,4) = micro%rsp(k,i,j)
         if (jnmb(4) >= 5) cx(k,4) = micro%csp(k,i,j)
      else
         if (k2(4) == 1) k1(4) = k + 1
      endif
   enddo
endif

! fill scratch arrays for aggregates

if (jnmb(5) >= 1) then
   do k = 2,m1-1
      if (micro%rap(k,i,j) >= rxmin) then
         k2(5) = k
         rx(k,5) = micro%rap(k,i,j)
         if (jnmb(5) >= 5) cx(k,5) = micro%cap(k,i,j)
      else
         if (k2(5) == 1) k1(5) = k + 1
      endif
   enddo
endif

! fill scratch arrays for graupel

if (jnmb(6) >= 1) then
   do k = 2,m1-1
      if (micro%rgp(k,i,j) >= rxmin) then
         k2(6) = k
         rx(k,6) = micro%rgp(k,i,j)
         qx(k,6) = micro%q6(k,i,j)
         qr(k,6) = qx(k,6) * rx(k,6)
         if (jnmb(6) >= 5) cx(k,6) = micro%cgp(k,i,j)
      else
         if (k2(6) == 1) k1(6) = k + 1
      endif
   enddo
endif

! fill scratch arrays for hail

if (jnmb(7) >= 1) then
   do k = 2,m1-1
      if (micro%rhp(k,i,j) >= rxmin) then
         k2(7) = k
         rx(k,7) = micro%rhp(k,i,j)
         qx(k,7) = micro%q7(k,i,j)
         qr(k,7) = qx(k,7) * rx(k,7)
         if (jnmb(7) >= 5) cx(k,7) = micro%chp(k,i,j)
      else
         if (k2(7) == 1) k1(7) = k + 1
      endif
   enddo
endif

! fill scratch arrays for drizzle
if (jnmb(8) >= 1) then
   do k = 2,m1-1
      if (micro%rdp(k,i,j) >= rxmin) then
         k2(8) = k
         rx(k,8) = micro%rdp(k,i,j)
         if (jnmb(8) >= 5) cx(k,8) = micro%cdp(k,i,j)
      else
         if (k2(8) == 1) k1(8) = k + 1
      endif
   enddo
endif

k3(1) = k2(1)
k3(3) = k2(3)
k3(8) = k2(8)

k1(9) = min(k1(1),k1(2),k1(8))
k2(9) = max(k2(1),k2(2),k2(8))
k1(10) = min(k1(3),k1(4),k1(5),k1(6),k1(7))
k2(10) = max(k2(3),k2(4),k2(5),k2(6),k2(7))
k1(11) = min(k1(9),k1(10))
k2(11) = max(k2(9),k2(10))

!Microphysics budget arrays
if(imbudget>=1) then
 do k = 1,m1
   xlatheatvap(k)    = micro%latheatvap(k,i,j)
   xlatheatfrz(k)    = micro%latheatfrz(k,i,j)
   xnuccldrt(k)      = micro%nuccldrt(k,i,j)
   xcld2raint(k)     = micro%cld2raint(k,i,j)
   xice2raint(k)     = micro%ice2raint(k,i,j)
   xnucicert(k)      = micro%nucicert(k,i,j)
   xvapliqt(k)       = micro%vapliqt(k,i,j)
   xvapicet(k)       = micro%vapicet(k,i,j)
   xevapliqt(k)      = micro%evapliqt(k,i,j)
   xevapicet(k)      = micro%evapicet(k,i,j)
   xfreezingt(k)     = micro%freezingt(k,i,j)
   xmeltingt(k)      = micro%meltingt(k,i,j)
   xmelticet(k)      = micro%melticet(k,i,j)
   xrimecldt(k)      = micro%rimecldt(k,i,j)
   xaggregatet(k)    = micro%aggregatet(k,i,j)
   xrain2icet(k)     = micro%rain2icet(k,i,j)
   xlatheatvapt(k)   = micro%latheatvapt(k,i,j)
   xlatheatfrzt(k)   = micro%latheatfrzt(k,i,j)
 enddo
endif
if(imbudget>=2) then
 do k = 1,m1
   xinuchomrt(k)     = micro%inuchomrt(k,i,j)
   xinuccontrt(k)    = micro%inuccontrt(k,i,j)
   xinucifnrt(k)     = micro%inucifnrt(k,i,j)
   xinuchazrt(k)     = micro%inuchazrt(k,i,j)
   xvapcldt(k)       = micro%vapcldt(k,i,j)
   xvapraint(k)      = micro%vapraint(k,i,j)
   xvapprist(k)      = micro%vapprist(k,i,j)
   xvapsnowt(k)      = micro%vapsnowt(k,i,j)
   xvapaggrt(k)      = micro%vapaggrt(k,i,j)
   xvapgraut(k)      = micro%vapgraut(k,i,j)
   xvaphailt(k)      = micro%vaphailt(k,i,j)
   xvapdrizt(k)      = micro%vapdrizt(k,i,j)
   xevapcldt(k)      = micro%evapcldt(k,i,j)
   xevapraint(k)     = micro%evapraint(k,i,j)
   xevapprist(k)     = micro%evapprist(k,i,j)
   xevapsnowt(k)     = micro%evapsnowt(k,i,j)
   xevapaggrt(k)     = micro%evapaggrt(k,i,j)
   xevapgraut(k)     = micro%evapgraut(k,i,j)
   xevaphailt(k)     = micro%evaphailt(k,i,j)
   xevapdrizt(k)     = micro%evapdrizt(k,i,j)
   xmeltprist(k)     = micro%meltprist(k,i,j)
   xmeltsnowt(k)     = micro%meltsnowt(k,i,j)
   xmeltaggrt(k)     = micro%meltaggrt(k,i,j)
   xmeltgraut(k)     = micro%meltgraut(k,i,j)
   xmelthailt(k)     = micro%melthailt(k,i,j)
   xrimecldsnowt(k)  = micro%rimecldsnowt(k,i,j)
   xrimecldaggrt(k)  = micro%rimecldaggrt(k,i,j)
   xrimecldgraut(k)  = micro%rimecldgraut(k,i,j)
   xrimecldhailt(k)  = micro%rimecldhailt(k,i,j)
   xrain2prt(k)      = micro%rain2prt(k,i,j)
   xrain2snt(k)      = micro%rain2snt(k,i,j)
   xrain2agt(k)      = micro%rain2agt(k,i,j)
   xrain2grt(k)      = micro%rain2grt(k,i,j)
   xrain2hat(k)      = micro%rain2hat(k,i,j)
   xaggrselfprist(k) = micro%aggrselfprist(k,i,j)
   xaggrselfsnowt(k) = micro%aggrselfsnowt(k,i,j)
   xaggrprissnowt(k) = micro%aggrprissnowt(k,i,j)
 enddo
endif
if(imbudget==3 .and. idust>=1) then
 do k = 1,m1
   xdust1cldrt(k)    = micro%dust1cldrt(k,i,j)
   xdust2cldrt(k)    = micro%dust2cldrt(k,i,j)
   xdust1drzrt(k)    = micro%dust1drzrt(k,i,j)
   xdust2drzrt(k)    = micro%dust2drzrt(k,i,j)
 enddo
endif

!Microphysics budget arrays
!ZERO OUT INSTANTANEOUS dT FROM LATENT HEATING 
if(imbudget>=1) then
 do k = 1,m1
  xlatheatvap(k) = 0.
  xlatheatfrz(k) = 0.
 enddo
endif
!ZERO OUT MICRO BUDGET PROCESSES AFTER ANALYSIS WRITE
!THEN BEGIN ACCUMULATING AGAIN
if(mod(time+0.001,frq).lt.dtlt .or. time.lt.0.001)then
 if(iprntstmt>=1 .and. imbudget>=1.and.i==2.and.j==2 .and. print_msg) &
   print*,'Resetting micro budgets',time,ngr
 if(imbudget>=1) then
  do k = 1,m1
   xnuccldrt(k)      = 0.
   xcld2raint(k)     = 0.
   xice2raint(k)     = 0.
   xnucicert(k)      = 0.
   xvapliqt(k)       = 0.
   xvapicet(k)       = 0.
   xevapliqt(k)      = 0.
   xevapicet(k)      = 0.
   xfreezingt(k)     = 0.
   xmeltingt(k)      = 0.
   xmelticet(k)      = 0.
   xrimecldt(k)      = 0.
   xaggregatet(k)    = 0.
   xrain2icet(k)     = 0.
   xlatheatvapt(k)   = 0.
   xlatheatfrzt(k)   = 0.
  enddo
 endif
 if(imbudget>=2) then
  do k = 1,m1
   xinuchomrt(k)     = 0.
   xinuccontrt(k)    = 0.
   xinucifnrt(k)     = 0.
   xinuchazrt(k)     = 0.
   xvapcldt(k)       = 0.
   xvapraint(k)      = 0.
   xvapprist(k)      = 0.
   xvapsnowt(k)      = 0.
   xvapaggrt(k)      = 0.
   xvapgraut(k)      = 0.
   xvaphailt(k)      = 0.
   xvapdrizt(k)      = 0.
   xevapcldt(k)      = 0.
   xevapraint(k)     = 0.
   xevapprist(k)     = 0.
   xevapsnowt(k)     = 0.
   xevapaggrt(k)     = 0.
   xevapgraut(k)     = 0.
   xevaphailt(k)     = 0.
   xevapdrizt(k)     = 0.
   xmeltprist(k)     = 0.
   xmeltsnowt(k)     = 0.
   xmeltaggrt(k)     = 0.
   xmeltgraut(k)     = 0.
   xmelthailt(k)     = 0.
   xrimecldsnowt(k)  = 0.
   xrimecldaggrt(k)  = 0.
   xrimecldgraut(k)  = 0.
   xrimecldhailt(k)  = 0.
   xrain2prt(k)      = 0.
   xrain2snt(k)      = 0.
   xrain2agt(k)      = 0.
   xrain2grt(k)      = 0.
   xrain2hat(k)      = 0.
   xaggrselfprist(k) = 0.
   xaggrselfsnowt(k) = 0.
   xaggrprissnowt(k) = 0.
  enddo
 endif
 if(imbudget==3 .and. idust>=1) then
  do k = 1,m1
   xdust1cldrt(k)    = 0.
   xdust2cldrt(k)    = 0.
   xdust1drzrt(k)    = 0.
   xdust2drzrt(k)    = 0.
  enddo
 endif
endif !If time to reset at analysis write

return
END SUBROUTINE range_check

!##############################################################################
Subroutine each_column (m1,k1,k2,rv,dn0)

use rconstants
use micphys

implicit none

integer :: m1,k,nt,ns
integer, dimension(11) :: k1,k2
real :: ck1,ck2,ck3,elsref,elsrefp,dplinv,eisref,eisrefp,dpiinv,relhum
real, dimension(m1) :: rv,dn0
real, external :: rslf
real, external :: rsif
real, external :: eslf
real, external :: eslpf
real, external :: esif
real, external :: esipf

data ck1,ck2,ck3/-4.818544e-3,1.407892e-4,-1.249986e-7/

do k = 2,m1-1
   rvlsair(k) = rslf(press(k),tair(k))
   rvisair(k) = rsif(press(k),tair(k))
   dn0i(k) = 1. / dn0(k)
   tairc(k)   = tair(k) - 273.15
   tx(k,1) = tairc(k)
   thrmcon(k) = ck1 + (ck2 + ck3 * tair(k)) * tair(k)
   dynvisc(k) = .1718e-4 + .49e-7 * tairc(k)
   
   ! Diagnose habit of pristine ice and snow

   nt = max(1,min(31,-nint(tairc(k))))
   relhum = min(1.,rv(k) / rvlsair(k))
   ns = max(1,nint(100. * relhum))
   jhcat(k,3) = jhabtab(nt,ns,1)
   jhcat(k,4) = jhabtab(nt,ns,2)

enddo

do k = k1(11),k2(11)
   vapdif(k)     = 2.14 * (tair(k) / 273.15) ** 1.94 / press(k)
   rdynvsci(k) = sqrt(1. / dynvisc(k))
   denfac(k) = sqrt(dn0i(k))

   colfacr(k) = colf * denfac(k) * dn0(k)
   colfacr2(k) = 2. * colfacr(k)

!Saleeby(2010): Loftus: remove density from: colfacc(k) = colfacr(k) * dn0(k)
   colfacc(k) = colfacr(k)

   colfacc2(k) = 2. * colfacc(k)

   tref(k,1)   = tairc(k) - min(25.,700. * (rvlsair(k) - rv(k)))
   sa(k,2) = thrmcon(k) * sa(k,1)
   sa(k,3) = thrmcon(k) * (tairstrc(k) + sa(k,1) * rvstr(k))

   sumuy(k) = 0.
   sumuz(k) = 0.
   sumvr(k) = 0.
enddo

do k = k1(9),k2(9)
   elsref       = eslf(tref(k,1))
   elsrefp      = eslpf(tref(k,1))
   dplinv       = 1. / (press(k) - elsref)
   rvsref (k,1) = .622 * elsref * dplinv
   rvsrefp(k,1) = .622 * elsrefp * dplinv * (1. + elsref * dplinv)

   sa(k,4) = rvsrefp(k,1) * tref(k,1) - rvsref(k,1)
   sa(k,6) = alvl * rvsrefp(k,1)
   sa(k,8) = alvl * sa(k,4)
enddo

do k = k1(10),k2(10)
   tref(k,2)    = min(0.,tref(k,1))
   eisref       = esif(tref(k,2))
   eisrefp      = esipf(tref(k,2))
   dpiinv       = 1. / (press(k) - eisref)
   rvsref (k,2) = .622 * eisref * dpiinv
   rvsrefp(k,2) = .622 * eisrefp * dpiinv * (1. + eisref * dpiinv)
   rvs0(k)      = 379.4 / (press(k) - 610.)

   sa(k,5) = rvsrefp(k,2) * tref(k,2) - rvsref(k,2)
   sa(k,7) = alvi * rvsrefp(k,2)
   sa(k,9) = alvi * sa(k,5)
   sh(k,3) = 0.
   sh(k,4) = 0.
   sh(k,5) = 0.

enddo

return
END SUBROUTINE each_column

!##############################################################################
Subroutine enemb (m1,k1,k2,lcat,dn0)

use micphys

implicit none

integer, dimension(11) :: k1,k2
integer :: m1,lcat,k,lhcat
real :: embi,parmi,embtemp,temprx,rxdif
real, dimension(m1) :: dn0

if (jnmb(lcat) == 2) then
   embi = 1. / emb(2,lcat)
   do k = k1(lcat),k2(lcat)
      cx(k,lcat) = rx(k,lcat) * embi
   enddo
elseif (jnmb(lcat) == 3) then
   do k = k1(lcat),k2(lcat)
      lhcat = jhcat(k,lcat)
      emb(k,lcat) = cfemb0(lhcat) * (dn0(k) * rx(k,lcat)) ** pwemb0(lhcat)
      cx(k,lcat) = cfen0(lhcat) * dn0i(k)  &
         * (dn0(k) * rx(k,lcat)) ** pwen0(lhcat)
   enddo
elseif (jnmb(lcat) == 4) then
   parmi = 1. / parm(lcat)
   do k = k1(lcat),k2(lcat)
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat) * parmi))
      cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
   enddo
elseif (jnmb(lcat) == 5) then
   do k = k1(lcat),k2(lcat)
      embtemp=rx(k,lcat)/max(cxmin,cx(k,lcat))
      emb(k,lcat) = max(emb0(lcat),min(emb1(lcat),rx(k,lcat)  &
         / max(cxmin,cx(k,lcat))))
      !Saleeby(2011): Use of single precision here allows for enemb to produce
      !artificial number concentration even when mass is in bounds. This is
      !generally small but could accumulate over time. Added the IF statement
      !to stop adjustment if mean size is in bounds.
      if( (embtemp > 0.0 .and. (embtemp < 0.999*emb0(lcat)   .or.  &
                                embtemp > 1.001*emb1(lcat))) .or.  &
          (rx(k,lcat) == 0.0 .and. cx(k,lcat) > 0.0 )        .or.  &
          (embtemp >= emb0(lcat) .and. embtemp <= emb1(lcat) .and. &
               cx(k,lcat)<cxmin .and. rx(k,lcat)>0.0))      then
        cx(k,lcat) = rx(k,lcat) / emb(k,lcat)
      endif
   enddo
endif

return
END SUBROUTINE enemb

!##############################################################################
Subroutine x02 (m1,k1,k2,lcat,dn0)

use rconstants
use micphys

implicit none

integer :: m1,lcat,k,lhcat,inc,idns
integer, dimension(11) :: k1,k2
real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed  &
       ,qrmltshed,shedmass,fracmloss,dn
real, dimension(m1) :: dn0

k1(lcat) = k1(11)
k2(lcat) = 1
do k = k1(11),k2(11)
   if (rx(k,lcat) >= rxmin) k2(lcat) = k
   if (k2(lcat) == 1 .and. rx(k,lcat) < rxmin) k1(lcat) = k + 1
enddo

if ((lcat == 2 .or. lcat >= 4) .and. (lcat .ne. 8)) then
   CALL enemb (m1,k1,k2,lcat,dn0)
endif

if (lcat == 2) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
! limit rain to under 48C and over -80C
      qx(k,lcat) = max(0.,min(1.6*alli,qx(k,lcat)))

      endif

   enddo

elseif (lcat == 3) then
!Allow pristine ice to melt to cloud1 since we assume that smaller
! particles will melt first. Perhaps need a way in the future to treat
! melting of larger pristine ice into cloud2.

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin) then

      rinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rinv

      CALL qtc (qx(k,lcat),tcoal,fracliq)

      rmelt = rx(k,lcat) * fracliq
      cmelt = cx(k,lcat) * fracliq

      !Aerosol and solubility tracking
      if(iccnlev>=2) then
         rxferratio = min(1.0, rmelt / rx(k,lcat))
         ccnmass = cnmhx(k,lcat) * rxferratio
         cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
         cnmhx(k,1)    = cnmhx(k,1) + ccnmass
         if(itrkepsilon==1) then
          scnmass = snmhx(k,lcat) * rxferratio
          snmhx(k,lcat) = snmhx(k,lcat) - scnmass
          snmhx(k,1)    = snmhx(k,1) + scnmass
         endif
         if(itrkdust==1)then
          dcnmass = dnmhx(k,lcat) * rxferratio
          dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
          dnmhx(k,1)    = dnmhx(k,1) + dcnmass
         endif
         if(itrkdustifn==1)then
          dinmass = dinhx(k,lcat) * rxferratio
          dinhx(k,lcat) = dinhx(k,lcat) - dinmass
          dinhx(k,1)    = dinhx(k,1) + dinmass
         endif
      endif

      rx(k,lcat) = rx(k,lcat) - rmelt
      rx(k,1) = rx(k,1) + rmelt
      cx(k,lcat) = cx(k,lcat) - cmelt
      cx(k,1) = cx(k,1) + cmelt

      if(imbudget >= 1) xmelticet(k)  = xmelticet(k)  + rmelt * budget_scalet
      if(imbudget >= 2) xmeltprist(k) = xmeltprist(k) + rmelt * budget_scalet

      endif

   enddo
!
! meyers - source for cloud aerosol number here?
!
elseif (lcat == 4 .or. lcat == 5) then
!Allow snow and aggregates to melt to graupel.
!Perhaps snow should melt to cloud2 and aggregates to rain??
!Change this??? move to rain instead ??? look at melting decisions in col2

   do k = k1(lcat),k2(lcat)

     if (rx(k,lcat) >= rxmin) then

      rinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rinv
      CALL qtc (qx(k,lcat),tcoal,fracliq)

      if (fracliq > 1.0e-6) then
         !Compute liquid portion to xfer to graupel
         rmelt = rx(k,lcat) * fracliq
         !Compute ice portion of xfer to graupel
         ricetor6 = min(rx(k,lcat) - rmelt,rmelt)

         !Aerosol and solubility tracking
         if(iccnlev>=2) then
            rxferratio = min(1.0, rmelt+ricetor6 / rx(k,lcat))
            ccnmass = cnmhx(k,lcat) * rxferratio
            cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
            cnmhx(k,6)    = cnmhx(k,6) + ccnmass
            if(itrkepsilon==1) then
             scnmass = snmhx(k,lcat) * rxferratio
             snmhx(k,lcat) = snmhx(k,lcat) - scnmass
             snmhx(k,6)    = snmhx(k,6) + scnmass
            endif
            if(itrkdust==1)then
             dcnmass = dnmhx(k,lcat) * rxferratio
             dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
             dnmhx(k,6)    = dnmhx(k,6) + dcnmass
            endif
            if(itrkdustifn==1)then
             dinmass = dinhx(k,lcat) * rxferratio
             dinhx(k,lcat) = dinhx(k,lcat) - dinmass
             dinhx(k,6)    = dinhx(k,6) + dinmass
            endif
         endif

         rx(k,lcat) = rx(k,lcat) - rmelt - ricetor6
         rx(k,6) = rx(k,6) + rmelt + ricetor6
         qr(k,6) = qr(k,6) + rmelt * alli
         qx(k,lcat) = 0.

! keep the above the same with ricetor6
! meyers - use sa melt table here? yes
         fracmloss = (rmelt + ricetor6) * rinv
         closs = enmlttab(int(200. * fracmloss) + 1,jhcat(k,lcat)) * cx(k,lcat)
         cx(k,lcat) = cx(k,lcat) - closs
         cx(k,6) = cx(k,6) + closs

!        if(imbudget >= 1) xmelticet used to track ice melting to rain. Do not include
!                          melting of snow and aggregates since they are transferred
!                          to graupel, not rain
         if(imbudget >= 2) then
           if(lcat==4) xmeltsnowt(k) = xmeltsnowt(k) + (rmelt + ricetor6) * budget_scalet
           if(lcat==5) xmeltaggrt(k) = xmeltaggrt(k) + (rmelt + ricetor6) * budget_scalet
         endif

      endif

     endif
   enddo


elseif (lcat == 6) then

   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
      CALL qtc (qx(k,lcat),tcoal,fracliq)

      if (fracliq > 0.95) then
         rx(k,2) = rx(k,2) + rx(k,6)
         qr(k,2) = qr(k,2) + rx(k,6) * alli
         cx(k,2) = cx(k,2) + cx(k,6)

         !Aerosol and solubility tracking
         if(iccnlev>=2) then
            rxferratio = 1.0
            ccnmass = cnmhx(k,lcat) * rxferratio
            cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
            cnmhx(k,2)    = cnmhx(k,2) + ccnmass
            if(itrkepsilon==1) then
             scnmass = snmhx(k,lcat) * rxferratio
             snmhx(k,lcat) = snmhx(k,lcat) - scnmass
             snmhx(k,2)    = snmhx(k,2) + scnmass
            endif
            if(itrkdust==1)then
             dcnmass = dnmhx(k,lcat) * rxferratio
             dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
             dnmhx(k,2)    = dnmhx(k,2) + dcnmass
            endif
            if(itrkdustifn==1)then
             dinmass = dinhx(k,lcat) * rxferratio
             dinhx(k,lcat) = dinhx(k,lcat) - dinmass
             dinhx(k,2)    = dinhx(k,2) + dinmass
            endif
         endif

         if(imbudget >= 1) xmelticet(k)  = xmelticet(k) + rx(k,6) * budget_scalet
         if(imbudget >= 2) xmeltgraut(k) = xmeltgraut(k) + rx(k,6) * budget_scalet

         rx(k,6) = 0.
         qr(k,6) = 0.
         qx(k,6) = 0.
         cx(k,6) = 0.
      endif

      endif

   enddo

elseif (lcat == 7) then

   shedmass = 5.236e-7
   do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) >= rxmin) then

      rxinv = 1. / rx(k,lcat)
      qx(k,lcat) = qr(k,lcat) * rxinv
      CALL qtc (qx(k,lcat),tcoal,fracliq)

      if (fracliq > 0.95) then
         rx(k,2) = rx(k,2) + rx(k,7)
         qr(k,2) = qr(k,2) + rx(k,7) * alli
         cx(k,2) = cx(k,2) + cx(k,7)

         !Aerosol and solubility tracking
         if(iccnlev>=2) then
            rxferratio = 1.0
            ccnmass = cnmhx(k,lcat) * rxferratio
            cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
            cnmhx(k,2)    = cnmhx(k,2) + ccnmass
            if(itrkepsilon==1) then
             scnmass = snmhx(k,lcat) * rxferratio
             snmhx(k,lcat) = snmhx(k,lcat) - scnmass
             snmhx(k,2)    = snmhx(k,2) + scnmass
            endif
            if(itrkdust==1)then
             dcnmass = dnmhx(k,lcat) * rxferratio
             dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
             dnmhx(k,2)    = dnmhx(k,2) + dcnmass
            endif
            if(itrkdustifn==1)then
             dinmass = dinhx(k,lcat) * rxferratio
             dinhx(k,lcat) = dinhx(k,lcat) - dinmass
             dinhx(k,2)    = dinhx(k,2) + dinmass
            endif
         endif

         if(imbudget >= 1) xmelticet(k)  = xmelticet(k) + rx(k,7) * budget_scalet
         if(imbudget >= 2) xmelthailt(k) = xmelthailt(k) + rx(k,7) * budget_scalet

         rx(k,7) = 0.
         qr(k,7) = 0.
         qx(k,7) = 0.
         cx(k,7) = 0.

!  take out following IF statement?

      elseif (fracliq > 0.3) then

         lhcat = jhcat(k,lcat)
         inc = nint(200. * fracliq) + 1
         dn = dnfac(lhcat) * emb(k,lcat) ** pwmasi(lhcat)
         idns = max(1,nint(1.e3 * dn * gnu(lcat)))
         rshed = rx(k,lcat) * shedtab(inc,idns)
         rmltshed = rshed
         qrmltshed = rmltshed * alli

         !Aerosol and solubility tracking
         if(iccnlev>=2) then
            rxferratio = min(1.0, rmltshed / rx(k,lcat))
            ccnmass = cnmhx(k,lcat) * rxferratio
            cnmhx(k,lcat) = cnmhx(k,lcat) - ccnmass
            cnmhx(k,2)    = cnmhx(k,2) + ccnmass
            if(itrkepsilon==1) then
             scnmass = snmhx(k,lcat) * rxferratio
             snmhx(k,lcat) = snmhx(k,lcat) - scnmass
             snmhx(k,2)    = snmhx(k,2) + scnmass
            endif
            if(itrkdust==1)then
             dcnmass = dnmhx(k,lcat) * rxferratio
             dnmhx(k,lcat) = dnmhx(k,lcat) - dcnmass
             dnmhx(k,2)    = dnmhx(k,2) + dcnmass
            endif
            if(itrkdustifn==1)then
             dinmass = dinhx(k,lcat) * rxferratio
             dinhx(k,lcat) = dinhx(k,lcat) - dinmass
             dinhx(k,2)    = dinhx(k,2) + dinmass
            endif
         endif

         rx(k,2) = rx(k,2) + rmltshed
         qr(k,2) = qr(k,2) + qrmltshed
         cx(k,2) = cx(k,2) + rshed / shedmass

         rx(k,lcat) = rx(k,lcat) - rmltshed
         qr(k,lcat) = qr(k,lcat) - qrmltshed
         qx(k,lcat) = qr(k,lcat) * (1./rx(k,lcat))

         if(imbudget >= 1) xmelticet(k)  = xmelticet(k) + rmltshed * budget_scalet
         if(imbudget >= 2) xmelthailt(k) = xmelthailt(k) + rmltshed * budget_scalet

      endif

      endif

   enddo

endif

return
END SUBROUTINE x02

!##############################################################################
Subroutine sedim (m1,lcat,ngr,k1,k2  &
   ,rtp,thp,theta,dn0,alphasfc  &
   ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew  &
   ,pcpfillc,pcpfillr,sfcpcp,allpcp)

use rconstants
use micphys

implicit none

integer :: m1,lcat,ngr,k1,k2,k,lhcat,iemb,kkf,kk
real :: colddn0,rolddn0,qrolddn0,dispemb,riemb,psfc,qpcpg,pcpg &
   ,dpcpg,dtlti,qnew,alphasfc,psfcdust,cnmold,snmold,dnmold,dinold,immerold &
   ,pctemp,prtemp,psfctemp,psfcaero
real, dimension(m1) :: rtp,thp,theta,dn0,cnew,rnew,qrnew &
                      ,cnmnew,snmnew,immernew,dnmnew,dinnew
real, dimension(m1,maxkfall,nembfall,nhcat,ndensrtgt,nband) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat,ndensrtgt,nband) :: sfcpcp
real, dimension(m1,nembfall,nhcat,ndensrtgt,nband) :: allpcp

snmold=0.
cnmold=0.
dinold=0.
dnmold=0.
immerold=0.
psfc=0.
psfcaero=0.
psfcdust=0.
pcprx(lcat) = 0.
do k = 1,m1
   rnew(k) = 0.
   cnew(k) = 0.
   qrnew(k) = 0.
   pcpvx(k,lcat) = 0.
   cnmnew(k) = 0.
   snmnew(k) = 0.
   dnmnew(k) = 0.
   dinnew(k) = 0.
   immernew(k) = 0.
enddo

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) > rxmin) then
      colddn0 = cx(k,lcat) * dn0(k) !Convert #/kg to #/m3
      rolddn0 = rx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
      qrolddn0 = qx(k,lcat) * rolddn0

      !Aerosol and solubility tracking
      if(iccnlev>=2) then
       cnmold = cnmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
       if(itrkepsilon==1) snmold = snmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
       if(itrkdust==1)    dnmold = dnmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
       if(itrkdustifn==1) dinold = dinhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
      endif

      !For tracking immersion freezing nuclei
      if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
         immerold = immerhx(k,lcat) * dn0(k) !Convert #/kg to #/m3

      !Here determine which set of powerlaws to use: the original
      ! ones in RAMS or the ones from R.Carver from Mitchell 1996.
      !The Mitchell power laws are not based at sea level so we adjust the
      ! density factor based at 0.7 kg/m3 instead of 1.0 kg/m3.
      if(iplaws==0) then
        dispemb = ch1(lhcat)  &
          * (emb(k,lcat)/cfmas(lhcat)) ** ch3(lhcat) * sqrt(dn0i(k))
      else
        dispemb = ch1(lhcat)  &
          * (emb(k,lcat)/cfmas(lhcat)) ** ch3(lhcat) * (0.7*dn0i(k))**.362
      endif

      riemb = 1. + ch2(lhcat,ngr) * log10(dispemb / dispemb0(lhcat,ngr))

      !Limiting iemb to max of nembfall
      iemb = min(nint(riemb),nembfall)

      if (k <= maxkfall) then
         psfctemp = sfcpcp(k,iemb,lhcat,1,1)
         psfc = rolddn0 * psfctemp
         !Aerosol and solubility tracking
         if(iccnlev>=2)	psfcaero = cnmold * psfctemp
         if(iccnlev>=2.and.itrkdust==1) psfcdust = dnmold * psfctemp
      endif

      do kkf = 1,min(maxkfall,k-1)
         kk = k + 1 - kkf
         pctemp = pcpfillc(k,kkf,iemb,lhcat,1,1)
         prtemp = pcpfillr(k,kkf,iemb,lhcat,1,1)
         cnew(kk)  = cnew(kk)  +  colddn0 * dn0i(kk) * pctemp
         rnew(kk)  = rnew(kk)  +  rolddn0 * dn0i(kk) * prtemp
         qrnew(kk) = qrnew(kk) + qrolddn0 * dn0i(kk) * prtemp
         !Aerosol and solubility tracking
         if(iccnlev>=2) then
          cnmnew(kk) = cnmnew(kk) + cnmold * dn0i(kk) * prtemp
          if(itrkepsilon==1) snmnew(kk) = snmnew(kk) + snmold * dn0i(kk) * prtemp
          if(itrkdust==1)    dnmnew(kk) = dnmnew(kk) + dnmold * dn0i(kk) * prtemp
          if(itrkdustifn==1) dinnew(kk) = dinnew(kk) + dinold * dn0i(kk) * prtemp
         endif
         !For tracking immersion freezing nuclei
         if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
           immernew(kk) = immernew(kk) + immerold * dn0i(kk) * pctemp
      enddo

      !Surface Precip rate
      if (k <= maxkfall) then
         qpcpg = qpcpg + psfc * qx(k,lcat)
         pcprx(lcat) = pcprx(lcat) + psfc
         !Aerosol accumulation rate
         if(iccnlev>=2) &
           pcpraerox(lcat) = pcpraerox(lcat) + psfcaero
         if(iccnlev>=2 .and. itrkdust==1) &
           pcprdustx(lcat) = pcprdustx(lcat) + psfcdust
      endif

      !Precip rate at all levels (mm/s)
      pcpvx(k,lcat)= rolddn0 * allpcp(k,iemb,lhcat,1,1) * dtlti

   endif
enddo

pcpg = pcpg + pcprx(lcat)
accpx(lcat) = pcprx(lcat)
dpcpg = dpcpg + pcprx(lcat) * alphasfc
pcprx(lcat) = pcprx(lcat) * dtlti 

!Aerosol accumulation tracking  
if(iccnlev>=2) then
  !Accumulate total aerosol on surface (kg/m2)
  accpaerox(lcat) = pcpraerox(lcat)
  pcpraerox(lcat) = pcpraerox(lcat) * dtlti
  !Accumulate dust on surface (kg/m2)
  if(itrkdust==1) then
    accpdustx(lcat) = pcprdustx(lcat)
    pcprdustx(lcat) = pcprdustx(lcat) * dtlti
  endif
endif

do k = 2,k2
   rtp(k) = rtp(k) + rnew(k) - rx(k,lcat)
   qnew = qrnew(k) / max(1.0e-20, rnew(k))

   tairc(k) = tairc(k) - thp(k) * thp(k)  &
      * (2820. * (rnew(k) - rx(k,lcat))  &
      - cpi * (qrnew(k) - qx(k,lcat) * rx(k,lcat)))  &
      / (max(tair(k), 253.) * theta(k))

   rx(k,lcat) = rnew(k)
   cx(k,lcat) = cnew(k)
   qx(k,lcat) = qnew

   !Aerosol and solubility tracking
   if(iccnlev>=2) then
     cnmhx(k,lcat) = cnmnew(k)
     if(itrkepsilon==1) snmhx(k,lcat) = snmnew(k)
     if(itrkdust==1)    dnmhx(k,lcat) = dnmnew(k)
     if(itrkdustifn==1) dinhx(k,lcat) = dinnew(k)
   endif

   !For tracking immersion freezing nuclei
   if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
      immerhx(k,lcat) = immernew(k)

   if(rx(k,lcat) < rxmin) then
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
      pcpvx(k,lcat) = 0.

      !Aerosol and solubility tracking
      if(iccnlev>=2) then
        cnmhx(k,lcat) = 0.
        if(itrkepsilon==1) snmhx(k,lcat) = 0.
        if(itrkdust==1)    dnmhx(k,lcat) = 0.
        if(itrkdustifn==1) dinhx(k,lcat) = 0.
      endif
      !For tracking immersion freezing nuclei
      if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
        immerhx(k,lcat) = 0.
   endif
enddo

return
END SUBROUTINE sedim

!##############################################################################
Subroutine sedim_trubin (m1,lcat,ngr,k1,k2  &
   ,rtp,thp,theta,dn0,alphasfc  &
   ,pcpg,qpcpg,dpcpg,dtlti,cnew,rnew,qrnew  &
   ,pcpfillc,pcpfillr,sfcpcp,allpcp,rtgt)

use rconstants
use micphys

implicit none

integer :: m1,lcat,ngr,k1,k2,k,lhcat,iemb1,idensrtgt,iband,kkf,kk
real :: colddn0,rolddn0,qrolddn0,riemb,psfc,qpcpg,pcpg,pcpvxtemp1,pcpvxtemp2 &
   ,dpcpg,dtlti,qnew,alphasfc,psfcdust,cnmold,snmold,dnmold,dinold,immerold &
   ,rtgt,embwt1,embwt2,dmb,dmode,bwt1,bwt2,bdenpowfac &
   ,psfcaero,densrtgt,denwt1,denwt2,psfctemp1,psfctemp2,pctemp1,pctemp2 &
   ,prtemp1,prtemp2
real, dimension(m1) :: rtp,thp,theta,dn0,cnew,rnew,qrnew &
   ,cnmnew,snmnew,dnmnew,dinnew,immernew
real, dimension(m1,maxkfall,nembfall,nhcat,ndensrtgt,nband) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat,ndensrtgt,nband) :: sfcpcp
real, dimension(m1,nembfall,nhcat,ndensrtgt,nband) :: allpcp

snmold=0.
cnmold=0.
dinold=0.
dnmold=0.
immerold=0.
psfc=0.
psfcaero=0.
psfcdust=0.
idensrtgt=0
pcprx(lcat) = 0.
do k = 1,m1
   rnew(k) = 0.
   cnew(k) = 0.
   qrnew(k) = 0.
   pcpvx(k,lcat) = 0.
   cnmnew(k) = 0.
   snmnew(k) = 0.
   dnmnew(k) = 0.
   dinnew(k) = 0.
   immernew(k) = 0.
enddo

do k = k1,k2
   lhcat = jhcat(k,lcat)

   if (rx(k,lcat) > rxmin) then
      colddn0 = cx(k,lcat) * dn0(k) !Convert #/kg to #/m3
      rolddn0 = rx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
      qrolddn0 = qx(k,lcat) * rolddn0

      !Aerosol and solubility tracking
      if(iccnlev>=2) then
       cnmold = cnmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3 
       if(itrkepsilon==1) snmold = snmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
       if(itrkdust==1)    dnmold = dnmhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
       if(itrkdustifn==1) dinold = dinhx(k,lcat) * dn0(k) !Convert kg/kg to kg/m3
      endif

      !For tracking immersion freezing nuclei
      if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
         immerold = immerhx(k,lcat) * dn0(k) !Convert #/kg to #/m3

      !Determine which Vt power law band to use for IPLAWS=2
      dmb=(emb(k,lcat)/cfmas(lhcat))**pwmasi(lhcat)
      dmode=dmb*(pwmas(lhcat)+gnu(lcat)-1.)
      iband=0
      if(dmode < bdiam(1,lhcat)) then
         iband=1
         bwt2 = 0.0
         bwt1 = 1.0
      elseif(dmode >= bdiam(1,lhcat) .and. dmode < bdiam(2,lhcat)) then
         iband=1
         bwt2 = (dmode - bdiam(1,lhcat)) / (bdiam(2,lhcat)-bdiam(1,lhcat))
         bwt1 = 1.0 - bwt2
      elseif(dmode >= bdiam(2,lhcat) .and. dmode < bdiam(3,lhcat)) then
         iband=2
         bwt2 = (dmode - bdiam(2,lhcat)) / (bdiam(3,lhcat)-bdiam(2,lhcat))
         bwt1 = 1.0 - bwt2
      elseif(dmode >= bdiam(3,lhcat)) then
         iband=2
         bwt2 = 1.0
         bwt1 = 0.0
      else
         print*,'Sedimentation diameter weights out of range: stop'
         stop
      endif
      bdenpowfac=bwt1*bdenpow(iband,lhcat)+bwt2*bdenpow(iband+1,lhcat)

      if(iband==0) then
        print*,'IBAND not set in sedimentation'
        stop
      endif

      ! Determine which set of powerlaws to use: the original ones in RAMS or
      ! those from R.Carver from Mitchell 1996.
      ! The Mitchell power laws are not based at sea level so we adjust the
      ! density factor based at 0.7 kg/m3 instead of 1.0 kg/m3.
      ! Include density and rtgt topography factors. Density where hydrometeors
      ! are falling is typically between 0.3-1.2 kg/m3 and rtgt from 0.5-1.0.
      if(iplaws==0) densrtgt = 10.0 * sqrt(dn0i(k)) / rtgt
      if(iplaws==1) densrtgt = 10.0 * (0.7*dn0i(k))**0.362 / rtgt
      if(iplaws==2) densrtgt = 10.0 * (0.7*dn0i(k))**bdenpowfac / rtgt
      if(ISCM>=1)then
       !SaleebySCM: CCPP models can go to really high altitude and low density which
       !can cause "densrtgt" >> ndensrtgt, so we cap this value.
       densrtgt = min(ndensrtgt,int(densrtgt))
      endif
      if(int(densrtgt)<1 .or. int(densrtgt)>ndensrtgt)then
        print*,'Bad sedimentation densrtgt index.'
        stop
      endif
      
      !Get interpolation weights for density,rtgt indexing
      idensrtgt = int(densrtgt)
      if(idensrtgt==ndensrtgt) then
        idensrtgt = idensrtgt - 1
        denwt1 = 0.0
        denwt2 = 1.0
      else
        denwt2 = densrtgt - float(idensrtgt)
        denwt1 = float(idensrtgt+1) - densrtgt
      endif

      !Determine mean mass index
      riemb = 1. + ch2(lhcat,ngr) * log10(emb(k,lcat)/emb0(lcat))
      if(int(riemb)<1 .or. int(riemb)>nembfall)then
        print*,'Bad sedimentation mean mass index.'
        stop
      endif

      !Limiting iemb to max of nembfall + set interpolation weights
      iemb1 = int(riemb) !mean mass index
      if(iemb1==nembfall) then
        iemb1 = iemb1 - 1
        embwt1 = 0.0
        embwt2 = 1.0
      else
        embwt2 = riemb - float(iemb1)
        embwt1 = float(iemb1+1) - riemb
      endif

      if (k <= maxkfall) then
         psfctemp1 = denwt1 * &
               ( bwt1 * (embwt1 * sfcpcp(k,iemb1  ,lhcat,idensrtgt  ,iband)    &
                       + embwt2 * sfcpcp(k,iemb1+1,lhcat,idensrtgt  ,iband))   &
               + bwt2 * (embwt1 * sfcpcp(k,iemb1  ,lhcat,idensrtgt  ,iband+1)  &
                       + embwt2 * sfcpcp(k,iemb1+1,lhcat,idensrtgt  ,iband+1)) )
         psfctemp2 = denwt2 * &
               ( bwt1 * (embwt1 * sfcpcp(k,iemb1  ,lhcat,idensrtgt+1,iband)    &
                       + embwt2 * sfcpcp(k,iemb1+1,lhcat,idensrtgt+1,iband))   &
               + bwt2 * (embwt1 * sfcpcp(k,iemb1  ,lhcat,idensrtgt+1,iband+1)  &
                       + embwt2 * sfcpcp(k,iemb1+1,lhcat,idensrtgt+1,iband+1)) )
         psfc = rolddn0 * (psfctemp1 + psfctemp2)
         !Aerosol accumulation
         if(iccnlev>=2) psfcaero = cnmold * (psfctemp1 + psfctemp2)
         if(iccnlev>=2 .and. itrkdust==1) psfcdust = dnmold*(psfctemp1+psfctemp2)
      endif

      do kkf = 1,min(maxkfall,k-1)
         kk = k + 1 - kkf
         pctemp1 = denwt1 * &
          ( bwt1 * (embwt1 * pcpfillc(k,kkf,iemb1  ,lhcat,idensrtgt  ,iband)    &
                  + embwt2 * pcpfillc(k,kkf,iemb1+1,lhcat,idensrtgt  ,iband))   &
          + bwt2 * (embwt1 * pcpfillc(k,kkf,iemb1  ,lhcat,idensrtgt  ,iband+1)  &
                  + embwt2 * pcpfillc(k,kkf,iemb1+1,lhcat,idensrtgt  ,iband+1)) )
         pctemp2 = denwt2 * &
          ( bwt1 * (embwt1 * pcpfillc(k,kkf,iemb1  ,lhcat,idensrtgt+1,iband)    &
                  + embwt2 * pcpfillc(k,kkf,iemb1+1,lhcat,idensrtgt+1,iband))   &
          + bwt2 * (embwt1 * pcpfillc(k,kkf,iemb1  ,lhcat,idensrtgt+1,iband+1)  &
                  + embwt2 * pcpfillc(k,kkf,iemb1+1,lhcat,idensrtgt+1,iband+1)) )
         prtemp1 = denwt1 * &
          ( bwt1 * (embwt1 * pcpfillr(k,kkf,iemb1  ,lhcat,idensrtgt  ,iband)    &
                  + embwt2 * pcpfillr(k,kkf,iemb1+1,lhcat,idensrtgt  ,iband))   &
          + bwt2 * (embwt1 * pcpfillr(k,kkf,iemb1  ,lhcat,idensrtgt  ,iband+1)  &
                  + embwt2 * pcpfillr(k,kkf,iemb1+1,lhcat,idensrtgt  ,iband+1)) )
         prtemp2 = denwt2 * &
          ( bwt1 * (embwt1 * pcpfillr(k,kkf,iemb1  ,lhcat,idensrtgt+1,iband)    &
                  + embwt2 * pcpfillr(k,kkf,iemb1+1,lhcat,idensrtgt+1,iband))   &
          + bwt2 * (embwt1 * pcpfillr(k,kkf,iemb1  ,lhcat,idensrtgt+1,iband+1)  &
                  + embwt2 * pcpfillr(k,kkf,iemb1+1,lhcat,idensrtgt+1,iband+1)) )
         cnew(kk)  = cnew(kk)  +  colddn0  * dn0i(kk) * (pctemp1 + pctemp2)
         rnew(kk)  = rnew(kk)  +  rolddn0  * dn0i(kk) * (prtemp1 + prtemp2)
         qrnew(kk) = qrnew(kk) +  qrolddn0 * dn0i(kk) * (prtemp1 + prtemp2)
         !Aerosol and solubility tracking
         if(iccnlev>=2) then
           cnmnew(kk) = cnmnew(kk) + cnmold  * dn0i(kk) * (prtemp1 + prtemp2)
           if(itrkepsilon==1) snmnew(kk) = snmnew(kk) + snmold  * dn0i(kk) &
                              * (prtemp1 + prtemp2)
           if(itrkdust==1)    dnmnew(kk) = dnmnew(kk) + dnmold  * dn0i(kk) &
                              * (prtemp1 + prtemp2)
           if(itrkdustifn==1) dinnew(kk) = dinnew(kk) + dinold  * dn0i(kk) &
                              * (prtemp1 + prtemp2)
         endif
         !For tracking immersion freezing nuclei
         if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
           immernew(kk) = immernew(kk) + immerold * dn0i(kk) * (pctemp1 + pctemp2)
      enddo

      !Surface Precip rate
      if (k <= maxkfall) then
         qpcpg = qpcpg + psfc * qx(k,lcat)
         pcprx(lcat) = pcprx(lcat) + psfc
         !Aerosol accumulation
         if(iccnlev>=2)	&
           pcpraerox(lcat) = pcpraerox(lcat) + psfcaero
         if(iccnlev>=2 .and. itrkdust==1) &
           pcprdustx(lcat) = pcprdustx(lcat) + psfcdust
      endif

      !Precip rate at all levels (mm/s)
      pcpvxtemp1 = denwt1 * &
        ( bwt1 * (embwt1 * allpcp(k,iemb1  ,lhcat,idensrtgt  ,iband)    &
                + embwt2 * allpcp(k,iemb1+1,lhcat,idensrtgt  ,iband))   &
        + bwt2 * (embwt1 * allpcp(k,iemb1  ,lhcat,idensrtgt  ,iband+1)  &
                + embwt2 * allpcp(k,iemb1+1,lhcat,idensrtgt  ,iband+1)) )
      pcpvxtemp2 = denwt2 * &
        ( bwt1 * (embwt1 * allpcp(k,iemb1  ,lhcat,idensrtgt+1,iband)    &
                + embwt2 * allpcp(k,iemb1+1,lhcat,idensrtgt+1,iband))   &
        + bwt2 * (embwt1 * allpcp(k,iemb1  ,lhcat,idensrtgt+1,iband+1)  &
                + embwt2 * allpcp(k,iemb1+1,lhcat,idensrtgt+1,iband+1)) )
      pcpvx(k,lcat)= rolddn0 * dtlti * (pcpvxtemp1 + pcpvxtemp2)

   endif
enddo

pcpg = pcpg + pcprx(lcat)
accpx(lcat) = pcprx(lcat)
dpcpg = dpcpg + pcprx(lcat) * alphasfc
pcprx(lcat) = pcprx(lcat) * dtlti 

!Aerosol accumulation tracking
if(iccnlev>=2) then
  !Accumulate total aerosol on surface (kg/m2)
  accpaerox(lcat) = pcpraerox(lcat)
  pcpraerox(lcat) = pcpraerox(lcat) * dtlti
  !Accumulate dust on surface (kg/m2)
  if(itrkdust==1) then
    accpdustx(lcat) = pcprdustx(lcat)
    pcprdustx(lcat) = pcprdustx(lcat) * dtlti
  endif
endif

do k = 2,k2
   rtp(k) = rtp(k) + rnew(k) - rx(k,lcat)
   qnew = qrnew(k) / max(1.0e-20, rnew(k))

   tairc(k) = tairc(k) - thp(k) * thp(k)  &
      * (2820. * (rnew(k) - rx(k,lcat))  &
      - cpi * (qrnew(k) - qx(k,lcat) * rx(k,lcat)))  &
      / (max(tair(k), 253.) * theta(k))

   rx(k,lcat) = rnew(k)
   cx(k,lcat) = cnew(k)
   qx(k,lcat) = qnew

   !Aerosol and solubility tracking
   if(iccnlev>=2) then
     cnmhx(k,lcat) = cnmnew(k)
     if(itrkepsilon==1) snmhx(k,lcat) = snmnew(k)
     if(itrkdust==1)    dnmhx(k,lcat) = dnmnew(k)
     if(itrkdustifn==1) dinhx(k,lcat) = dinnew(k)
   endif

   !For tracking immersion freezing nuclei
   if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
      immerhx(k,lcat) = immernew(k)

   if(rx(k,lcat) < rxmin) then
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
      pcpvx(k,lcat) = 0.
      !Aerosol and solubility tracking
      if(iccnlev>=2) then
        cnmhx(k,lcat) = 0.
        if(itrkepsilon==1) snmhx(k,lcat) = 0.
        if(itrkdust==1)    dnmhx(k,lcat) = 0.
        if(itrkdustifn==1) dinhx(k,lcat) = 0.
      endif
      !For tracking immersion freezing nuclei
      if(iifn==3 .and. iccnlev>=1 .and. (lcat==1.or.lcat==2.or.lcat==8)) &
        immerhx(k,lcat) = 0.
   endif
enddo
  
return
END SUBROUTINE sedim_trubin
