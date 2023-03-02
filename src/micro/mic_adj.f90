!##############################################################################
Subroutine negadj1 ()

use mem_basic
use mem_micro
use mem_grid
use micphys, only:level
use node_mod, only:mzp,mxp,myp

implicit none

integer :: m1,m2,m3

if (level == 0) return

if (level <= 3) then
   CALL adj1 (mzp,mxp,myp,basic_g(ngrid)%rtp(1,1,1),micro_g(ngrid),ngrid)
elseif (level == 4) then
   CALL adj1_bin (mzp,mxp,myp,basic_g(ngrid)%rtp(1,1,1),micro_g(ngrid))
endif

return
END SUBROUTINE negadj1

!##############################################################################
Subroutine adj1 (m1,m2,m3,rtp,micro,ngr)

use mem_micro
use micphys
use node_mod, only:mi0,mj0

implicit none

integer :: m1,m2,m3,toomany,i,j,k,lcat,ngr
real :: frac,cnmhx_num,zerocheck,rxloss
real, dimension(m1,m2,m3) :: rtp
real, dimension(m1) :: rtemp
type (micro_vars) :: micro

if (level .eq. 0) return

do lcat = 1,ncat
   do k = 1,m1
      rx(k,lcat) = 0.
      cx(k,lcat) = 0.
      qx(k,lcat) = 0.
      pcpvx(k,lcat) = 0.
      cnmhx(k,lcat) = 0.
      snmhx(k,lcat) = 0.
      dnmhx(k,lcat) = 0.
      dinhx(k,lcat) = 0.
      immerhx(k,lcat) = 0.
   enddo
enddo

do j = 1,m3
   do i = 1,m2

   !Copy global variables mixing ratios to temporary variables
   if (jnmb(1) > 0) CALL ae1kmic (1,m1,rx(1,1),micro%rcp(1,i,j))
   if (jnmb(2) > 0) CALL ae1kmic (1,m1,rx(1,2),micro%rrp(1,i,j))
   if (jnmb(3) > 0) CALL ae1kmic (1,m1,rx(1,3),micro%rpp(1,i,j))
   if (jnmb(4) > 0) CALL ae1kmic (1,m1,rx(1,4),micro%rsp(1,i,j))
   if (jnmb(5) > 0) CALL ae1kmic (1,m1,rx(1,5),micro%rap(1,i,j))
   if (jnmb(6) > 0) CALL ae1kmic (1,m1,rx(1,6),micro%rgp(1,i,j))
   if (jnmb(7) > 0) CALL ae1kmic (1,m1,rx(1,7),micro%rhp(1,i,j))
   if (jnmb(8) > 0) CALL ae1kmic (1,m1,rx(1,8),micro%rdp(1,i,j))

   !Copy global variables internal hydrometeor energy to temp vars
   if (jnmb(2) > 0) CALL ae1kmic (1,m1,qx(1,2),micro%q2(1,i,j))
   if (jnmb(6) > 0) CALL ae1kmic (1,m1,qx(1,6),micro%q6(1,i,j))
   if (jnmb(7) > 0) CALL ae1kmic (1,m1,qx(1,7),micro%q7(1,i,j))

   !Copy global variables number concentrations to temporary variables
   if (jnmb(1) >= 5) CALL ae1kmic (1,m1,cx(1,1),micro%ccp(1,i,j))
   if (jnmb(2) >= 5) CALL ae1kmic (1,m1,cx(1,2),micro%crp(1,i,j))
   if (jnmb(3) >= 5) CALL ae1kmic (1,m1,cx(1,3),micro%cpp(1,i,j))
   if (jnmb(4) >= 5) CALL ae1kmic (1,m1,cx(1,4),micro%csp(1,i,j))
   if (jnmb(5) >= 5) CALL ae1kmic (1,m1,cx(1,5),micro%cap(1,i,j))
   if (jnmb(6) >= 5) CALL ae1kmic (1,m1,cx(1,6),micro%cgp(1,i,j))
   if (jnmb(7) >= 5) CALL ae1kmic (1,m1,cx(1,7),micro%chp(1,i,j))
   if (jnmb(8) >= 5) CALL ae1kmic (1,m1,cx(1,8),micro%cdp(1,i,j))

   !Copy global variable 3D precip rates to temporary variables
   if (jnmb(2) > 0) CALL ae1kmic (1,m1,pcpvx(1,2),micro%pcpvr(1,i,j))
   if (jnmb(3) > 0) CALL ae1kmic (1,m1,pcpvx(1,3),micro%pcpvp(1,i,j))
   if (jnmb(4) > 0) CALL ae1kmic (1,m1,pcpvx(1,4),micro%pcpvs(1,i,j))
   if (jnmb(5) > 0) CALL ae1kmic (1,m1,pcpvx(1,5),micro%pcpva(1,i,j))
   if (jnmb(6) > 0) CALL ae1kmic (1,m1,pcpvx(1,6),micro%pcpvg(1,i,j))
   if (jnmb(7) > 0) CALL ae1kmic (1,m1,pcpvx(1,7),micro%pcpvh(1,i,j))
   if (jnmb(8) > 0) CALL ae1kmic (1,m1,pcpvx(1,8),micro%pcpvd(1,i,j))

   !Aerosol and solubility tracking
   if(iccnlev>=2) then
     if(jnmb(1)>0) CALL ae1kmic (1,m1,cnmhx(1,1),micro%cnmcp(1,i,j))
     if(jnmb(2)>0) CALL ae1kmic (1,m1,cnmhx(1,2),micro%cnmrp(1,i,j))
     if(jnmb(3)>0) CALL ae1kmic (1,m1,cnmhx(1,3),micro%cnmpp(1,i,j))
     if(jnmb(4)>0) CALL ae1kmic (1,m1,cnmhx(1,4),micro%cnmsp(1,i,j))
     if(jnmb(5)>0) CALL ae1kmic (1,m1,cnmhx(1,5),micro%cnmap(1,i,j))
     if(jnmb(6)>0) CALL ae1kmic (1,m1,cnmhx(1,6),micro%cnmgp(1,i,j))
     if(jnmb(7)>0) CALL ae1kmic (1,m1,cnmhx(1,7),micro%cnmhp(1,i,j))
     if(jnmb(8)>0) CALL ae1kmic (1,m1,cnmhx(1,8),micro%cnmdp(1,i,j))
     if(itrkepsilon==1) then
       if(jnmb(1)>0) CALL ae1kmic (1,m1,snmhx(1,1),micro%snmcp(1,i,j))
       if(jnmb(2)>0) CALL ae1kmic (1,m1,snmhx(1,2),micro%snmrp(1,i,j))
       if(jnmb(3)>0) CALL ae1kmic (1,m1,snmhx(1,3),micro%snmpp(1,i,j))
       if(jnmb(4)>0) CALL ae1kmic (1,m1,snmhx(1,4),micro%snmsp(1,i,j))
       if(jnmb(5)>0) CALL ae1kmic (1,m1,snmhx(1,5),micro%snmap(1,i,j))
       if(jnmb(6)>0) CALL ae1kmic (1,m1,snmhx(1,6),micro%snmgp(1,i,j))
       if(jnmb(7)>0) CALL ae1kmic (1,m1,snmhx(1,7),micro%snmhp(1,i,j))
       if(jnmb(8)>0) CALL ae1kmic (1,m1,snmhx(1,8),micro%snmdp(1,i,j))
     endif
     if(itrkdust==1) then
       if(jnmb(1)>0) CALL ae1kmic (1,m1,dnmhx(1,1),micro%dnmcp(1,i,j))
       if(jnmb(2)>0) CALL ae1kmic (1,m1,dnmhx(1,2),micro%dnmrp(1,i,j))
       if(jnmb(3)>0) CALL ae1kmic (1,m1,dnmhx(1,3),micro%dnmpp(1,i,j))
       if(jnmb(4)>0) CALL ae1kmic (1,m1,dnmhx(1,4),micro%dnmsp(1,i,j))
       if(jnmb(5)>0) CALL ae1kmic (1,m1,dnmhx(1,5),micro%dnmap(1,i,j))
       if(jnmb(6)>0) CALL ae1kmic (1,m1,dnmhx(1,6),micro%dnmgp(1,i,j))
       if(jnmb(7)>0) CALL ae1kmic (1,m1,dnmhx(1,7),micro%dnmhp(1,i,j))
       if(jnmb(8)>0) CALL ae1kmic (1,m1,dnmhx(1,8),micro%dnmdp(1,i,j))
     endif
     if(itrkdustifn==1) then
       if(jnmb(1)>0) CALL ae1kmic (1,m1,dinhx(1,1),micro%dincp(1,i,j))
       if(jnmb(2)>0) CALL ae1kmic (1,m1,dinhx(1,2),micro%dinrp(1,i,j))
       if(jnmb(3)>0) CALL ae1kmic (1,m1,dinhx(1,3),micro%dinpp(1,i,j))
       if(jnmb(4)>0) CALL ae1kmic (1,m1,dinhx(1,4),micro%dinsp(1,i,j))
       if(jnmb(5)>0) CALL ae1kmic (1,m1,dinhx(1,5),micro%dinap(1,i,j))
       if(jnmb(6)>0) CALL ae1kmic (1,m1,dinhx(1,6),micro%dingp(1,i,j))
       if(jnmb(7)>0) CALL ae1kmic (1,m1,dinhx(1,7),micro%dinhp(1,i,j))
       if(jnmb(8)>0) CALL ae1kmic (1,m1,dinhx(1,8),micro%dindp(1,i,j))
     endif
   endif

   !Zero out very small mixing ratios
   do lcat = 1,ncat
     do k = 1,m1

       zerocheck=0
       cxloss=0
       rxloss=0

       !Zero out hydrometeor fields if they are below a min threshold for 2-moment
       if(jnmb(lcat)>=5 .and. (rx(k,lcat) < rxmin .or. cx(k,lcat) <= 0.0)) then
         zerocheck=1
         cxloss = cx(k,lcat)
         rxloss = rx(k,lcat)
         rx(k,lcat) = 0.
         cx(k,lcat) = 0.
         qx(k,lcat) = 0.
       endif
       !Finish zeroing fields or do this if 1-moment
       if(rx(k,lcat) < rxmin) then
         zerocheck=1
         rx(k,lcat) = 0.
         qx(k,lcat) = 0.
         pcpvx(k,lcat) = 0.
       endif

       !Aerosol and solubility tracking
       if(zerocheck==1 .and. iccnlev>=2 .and. cnmhx(k,lcat)>0.0) then

         !Determine how many number to restore
         if(cxloss > 0.0) then
           !If there are number to restore, compute median radius (rg)
           rg=((0.23873/aero_rhosol(aerocat)*cnmhx(k,lcat) / &
               max(1.e-10,cxloss))**(0.3333))/aero_rg2rm(aerocat)
           if(rg < 0.01e-6) rg = 0.01e-6
           if(rg > 6.50e-6) rg = 6.50e-6
         else
           !Compute a number to restore based on assumed 1 micron rg
           rg = 1.0e-6
         endif
         !Compute the number to restore
         cnmhx_num = cnmhx(k,lcat) * (0.23873/aero_rhosol(aerocat)) / &
              ((rg * aero_rg2rm(aerocat)) ** 3.)

         !Restore aerosols to regenerated category
         if(rg <= 0.96e-6) then
           micro%regen_aero1_mp(k,i,j) = micro%regen_aero1_mp(k,i,j) + cnmhx(k,lcat)
           micro%regen_aero1_np(k,i,j) = micro%regen_aero1_np(k,i,j) + cnmhx_num
         else
           micro%regen_aero2_mp(k,i,j) = micro%regen_aero2_mp(k,i,j) + cnmhx(k,lcat)
           micro%regen_aero2_np(k,i,j) = micro%regen_aero2_np(k,i,j) + cnmhx_num  
         endif

         !Statement to check the restoration of aerosol data
         !print*,'rmin',rxloss,cnmhx(k,lcat),cxloss,cnmhx_num,rg*1.e6

         !Zero out aerosol masses and such within hydrometeors
         cnmhx(k,lcat) = 0.
         if(itrkepsilon==1) snmhx(k,lcat) = 0.
         if(itrkdust==1)    dnmhx(k,lcat) = 0.
         if(itrkdustifn==1) dinhx(k,lcat) = 0.
       endif

       !Aerosol and solubility tracking
       if(iccnlev>=2) then
         if(cnmhx(k,lcat)<minmashydro) cnmhx(k,lcat) = 0.
         if(itrkepsilon==1 .and. snmhx(k,lcat)<minmashydro) snmhx(k,lcat) = 0.
         if(itrkdust==1    .and. dnmhx(k,lcat)<minmashydro) dnmhx(k,lcat) = 0.
         if(itrkdustifn==1 .and. dinhx(k,lcat)<minmashydro) dinhx(k,lcat) = 0.
         if(itrkepsilon==1 .and. snmhx(k,lcat)>cnmhx(k,lcat)) &
           snmhx(k,lcat)=0.99*cnmhx(k,lcat)
       endif

     enddo
   enddo

   !Prevent condensate from being > total water
   do k = 1,m1
     rtp(k,i,j) = max(0.,rtp(k,i,j))
     rtemp(k) = 1.001 * (rx(k,1)+ rx(k,2) + rx(k,3)  &
        + rx(k,4) + rx(k,5) + rx(k,6) + rx(k,7) + rx(k,8))
   enddo
   do k = 1,m1
     if (rtemp(k) > rtp(k,i,j)) then
       frac = rtp(k,i,j) / max(rxmin,rtemp(k))
       do lcat = 1,ncat
         rx(k,lcat) = rx(k,lcat) * frac
         !Aerosol and solubility tracking
         if(iccnlev>=2) then
           cnmhx(k,lcat) = cnmhx(k,lcat) * frac
           if(itrkepsilon==1) snmhx(k,lcat) = snmhx(k,lcat) * frac
           if(itrkdust==1)    dnmhx(k,lcat) = dnmhx(k,lcat) * frac
           if(itrkdustifn==1) dinhx(k,lcat) = dinhx(k,lcat) * frac
         endif
       enddo
     endif
   enddo

   !Reproportion number concentration if mixing ratio was reproportioned
   !so that mean diameters stay the same
   if(jnmb(1)>=5) CALL ae1mic (1,m1,cx(1,1),micro%rcp(1,i,j),rx(1,1))
   if(jnmb(2)>=5) CALL ae1mic (1,m1,cx(1,2),micro%rrp(1,i,j),rx(1,2))
   if(jnmb(3)>=5) CALL ae1mic (1,m1,cx(1,3),micro%rpp(1,i,j),rx(1,3))
   if(jnmb(4)>=5) CALL ae1mic (1,m1,cx(1,4),micro%rsp(1,i,j),rx(1,4))
   if(jnmb(5)>=5) CALL ae1mic (1,m1,cx(1,5),micro%rap(1,i,j),rx(1,5))
   if(jnmb(6)>=5) CALL ae1mic (1,m1,cx(1,6),micro%rgp(1,i,j),rx(1,6))
   if(jnmb(7)>=5) CALL ae1mic (1,m1,cx(1,7),micro%rhp(1,i,j),rx(1,7))
   if(jnmb(8)>=5) CALL ae1mic (1,m1,cx(1,8),micro%rdp(1,i,j),rx(1,8))

   !Copy mixing ratio back to global variables
   if(jnmb(1) > 0) CALL ae1kmic (1,m1,micro%rcp(1,i,j),rx(1,1))
   if(jnmb(2) > 0) CALL ae1kmic (1,m1,micro%rrp(1,i,j),rx(1,2))
   if(jnmb(3) > 0) CALL ae1kmic (1,m1,micro%rpp(1,i,j),rx(1,3))
   if(jnmb(4) > 0) CALL ae1kmic (1,m1,micro%rsp(1,i,j),rx(1,4))
   if(jnmb(5) > 0) CALL ae1kmic (1,m1,micro%rap(1,i,j),rx(1,5))
   if(jnmb(6) > 0) CALL ae1kmic (1,m1,micro%rgp(1,i,j),rx(1,6))
   if(jnmb(7) > 0) CALL ae1kmic (1,m1,micro%rhp(1,i,j),rx(1,7))
   if(jnmb(8) > 0) CALL ae1kmic (1,m1,micro%rdp(1,i,j),rx(1,8))

   !Copy internal hydrometeor energy back to global vars
   if(jnmb(2) > 0) CALL ae1kmic (1,m1,micro%q2(1,i,j),qx(1,2))
   if(jnmb(6) > 0) CALL ae1kmic (1,m1,micro%q6(1,i,j),qx(1,6))
   if(jnmb(7) > 0) CALL ae1kmic (1,m1,micro%q7(1,i,j),qx(1,7))

   !Copy number concentration back to global variables
   if(jnmb(1) >= 5) CALL ae1kmic (1,m1,micro%ccp(1,i,j),cx(1,1))
   if(jnmb(2) >= 5) CALL ae1kmic (1,m1,micro%crp(1,i,j),cx(1,2))
   if(jnmb(3) >= 5) CALL ae1kmic (1,m1,micro%cpp(1,i,j),cx(1,3))
   if(jnmb(4) >= 5) CALL ae1kmic (1,m1,micro%csp(1,i,j),cx(1,4))
   if(jnmb(5) >= 5) CALL ae1kmic (1,m1,micro%cap(1,i,j),cx(1,5))
   if(jnmb(6) >= 5) CALL ae1kmic (1,m1,micro%cgp(1,i,j),cx(1,6))
   if(jnmb(7) >= 5) CALL ae1kmic (1,m1,micro%chp(1,i,j),cx(1,7))
   if(jnmb(8) >= 5) CALL ae1kmic (1,m1,micro%cdp(1,i,j),cx(1,8))

   !Copy 3D precip rates back to global variables
   if(jnmb(2) > 0) CALL ae1kmic (1,m1,micro%pcpvr(1,i,j),pcpvx(1,2))
   if(jnmb(3) > 0) CALL ae1kmic (1,m1,micro%pcpvp(1,i,j),pcpvx(1,3))
   if(jnmb(4) > 0) CALL ae1kmic (1,m1,micro%pcpvs(1,i,j),pcpvx(1,4))
   if(jnmb(5) > 0) CALL ae1kmic (1,m1,micro%pcpva(1,i,j),pcpvx(1,5))
   if(jnmb(6) > 0) CALL ae1kmic (1,m1,micro%pcpvg(1,i,j),pcpvx(1,6))
   if(jnmb(7) > 0) CALL ae1kmic (1,m1,micro%pcpvh(1,i,j),pcpvx(1,7))
   if(jnmb(8) > 0) CALL ae1kmic (1,m1,micro%pcpvd(1,i,j),pcpvx(1,8))

   !Aerosol and solubility tracking; Copy back to global variables
   if(iccnlev>=2)then
     if(jnmb(1)>0) CALL ae1kmic (1,m1,micro%cnmcp(1,i,j),cnmhx(1,1))
     if(jnmb(2)>0) CALL ae1kmic (1,m1,micro%cnmrp(1,i,j),cnmhx(1,2))
     if(jnmb(3)>0) CALL ae1kmic (1,m1,micro%cnmpp(1,i,j),cnmhx(1,3))
     if(jnmb(4)>0) CALL ae1kmic (1,m1,micro%cnmsp(1,i,j),cnmhx(1,4))
     if(jnmb(5)>0) CALL ae1kmic (1,m1,micro%cnmap(1,i,j),cnmhx(1,5))
     if(jnmb(6)>0) CALL ae1kmic (1,m1,micro%cnmgp(1,i,j),cnmhx(1,6))
     if(jnmb(7)>0) CALL ae1kmic (1,m1,micro%cnmhp(1,i,j),cnmhx(1,7))
     if(jnmb(8)>0) CALL ae1kmic (1,m1,micro%cnmdp(1,i,j),cnmhx(1,8))
     if(itrkepsilon==1)then
      if(jnmb(1)>0) CALL ae1kmic (1,m1,micro%snmcp(1,i,j),snmhx(1,1))
      if(jnmb(2)>0) CALL ae1kmic (1,m1,micro%snmrp(1,i,j),snmhx(1,2))
      if(jnmb(3)>0) CALL ae1kmic (1,m1,micro%snmpp(1,i,j),snmhx(1,3))
      if(jnmb(4)>0) CALL ae1kmic (1,m1,micro%snmsp(1,i,j),snmhx(1,4))
      if(jnmb(5)>0) CALL ae1kmic (1,m1,micro%snmap(1,i,j),snmhx(1,5))
      if(jnmb(6)>0) CALL ae1kmic (1,m1,micro%snmgp(1,i,j),snmhx(1,6))
      if(jnmb(7)>0) CALL ae1kmic (1,m1,micro%snmhp(1,i,j),snmhx(1,7))
      if(jnmb(8)>0) CALL ae1kmic (1,m1,micro%snmdp(1,i,j),snmhx(1,8))
     endif
     if(itrkdust==1)then
      if(jnmb(1)>0) CALL ae1kmic (1,m1,micro%dnmcp(1,i,j),dnmhx(1,1))
      if(jnmb(2)>0) CALL ae1kmic (1,m1,micro%dnmrp(1,i,j),dnmhx(1,2))
      if(jnmb(3)>0) CALL ae1kmic (1,m1,micro%dnmpp(1,i,j),dnmhx(1,3))
      if(jnmb(4)>0) CALL ae1kmic (1,m1,micro%dnmsp(1,i,j),dnmhx(1,4))
      if(jnmb(5)>0) CALL ae1kmic (1,m1,micro%dnmap(1,i,j),dnmhx(1,5))
      if(jnmb(6)>0) CALL ae1kmic (1,m1,micro%dnmgp(1,i,j),dnmhx(1,6))
      if(jnmb(7)>0) CALL ae1kmic (1,m1,micro%dnmhp(1,i,j),dnmhx(1,7))
      if(jnmb(8)>0) CALL ae1kmic (1,m1,micro%dnmdp(1,i,j),dnmhx(1,8))
     endif
     if(itrkdustifn==1)then
      if(jnmb(1)>0) CALL ae1kmic (1,m1,micro%dincp(1,i,j),dinhx(1,1))
      if(jnmb(2)>0) CALL ae1kmic (1,m1,micro%dinrp(1,i,j),dinhx(1,2))
      if(jnmb(3)>0) CALL ae1kmic (1,m1,micro%dinpp(1,i,j),dinhx(1,3))
      if(jnmb(4)>0) CALL ae1kmic (1,m1,micro%dinsp(1,i,j),dinhx(1,4))
      if(jnmb(5)>0) CALL ae1kmic (1,m1,micro%dinap(1,i,j),dinhx(1,5))
      if(jnmb(6)>0) CALL ae1kmic (1,m1,micro%dingp(1,i,j),dinhx(1,6))
      if(jnmb(7)>0) CALL ae1kmic (1,m1,micro%dinhp(1,i,j),dinhx(1,7))
      if(jnmb(8)>0) CALL ae1kmic (1,m1,micro%dindp(1,i,j),dinhx(1,8))
     endif
   endif

   !Keep immersion-freezing IN positive-definite and keep
   !from exceeded total number of droplets
   if(iifn==3 .and. iccnlev>=1) then
    if (jnmb(1) >= 5) CALL ae1kmic (1,m1,ifnnucx(1),micro%ifnnucp(1,i,j))
    if (jnmb(1) >= 5) CALL ae1kmic (1,m1,immerhx(1,1),micro%immercp(1,i,j))
    if (jnmb(8) >= 5) CALL ae1kmic (1,m1,immerhx(1,8),micro%immerdp(1,i,j))
    if (jnmb(2) >= 5) CALL ae1kmic (1,m1,immerhx(1,2),micro%immerrp(1,i,j))
    do k = 1,m1
     if(ifnnucx(k) < minifn) ifnnucx(k) = 0.
     do lcat=1,ncat
     if(lcat==1.or.lcat==8.or.lcat==2)then !limit to liquid species
      cxtemp = 0.0
      if(lcat==1 .and. jnmb(1)>=5) cxtemp = micro%ccp(k,i,j)
      if(lcat==8 .and. jnmb(8)>=5) cxtemp = micro%cdp(k,i,j)
      if(lcat==2 .and. jnmb(2)>=5) cxtemp = micro%crp(k,i,j)
      if(immerhx(k,lcat) > cxtemp) immerhx(k,lcat) = 0.9999 * cxtemp
      if(immerhx(k,lcat) < minifn) immerhx(k,lcat) = 0.
     endif
     enddo
    enddo
    if (jnmb(1) >= 5) CALL ae1kmic (1,m1,micro%ifnnucp(1,i,j),ifnnucx(1))
    if (jnmb(1) >= 5) CALL ae1kmic (1,m1,micro%immercp(1,i,j),immerhx(1,1))
    if (jnmb(8) >= 5) CALL ae1kmic (1,m1,micro%immerdp(1,i,j),immerhx(1,8))
    if (jnmb(2) >= 5) CALL ae1kmic (1,m1,micro%immerrp(1,i,j),immerhx(1,2))
   endif

 !CHECK FOR NEGATIVE VALUES OF AEROSOL SPECIES OVER ALL LEVELS
 do k = 1,m1
  toomany=0

  !TOO LARGE OF VALUES OF AEROSOLS

  if(level==3 .and. ipris>=5 .and. (iifn==1.or.iifn==2)) then
   if(micro%cifnp(k,i,j) > maxaero) then
     print*,"Too many IFN:",micro%cifnp(k,i,j)
     toomany=1
   endif
  endif
  if(iaerosol>0)then
   if(micro%cccnp(k,i,j) > maxaero) then
     print*,"Too many CCN:",micro%cccnp(k,i,j)
     toomany=1  
   endif
   if(micro%gccnp(k,i,j) > maxaero) then
     print*,"Too many GCCN:",micro%gccnp(k,i,j)
     toomany=1
   endif
  endif
  if(idust>0) then
   if(micro%md1np(k,i,j) > maxaero) then
     print*,"Too many Dust1:",micro%md1np(k,i,j)
     toomany=1
   endif
   if(micro%md2np(k,i,j) > maxaero) then
     print*,"Too many Dust2:",micro%md2np(k,i,j)
     toomany=1
   endif
  endif
  if(iabcarb>0) then
   if(micro%abc1np(k,i,j) > maxaero) then
     print*,"Too many Absorbing Carbon 1:",micro%abc1np(k,i,j)
     toomany=1
   endif
   if(micro%abc2np(k,i,j) > maxaero) then
     print*,"Too many Absorbing Carbon 2:",micro%abc2np(k,i,j)
     toomany=1
   endif
  endif
  if(isalt>0) then
   if(micro%salt_film_np(k,i,j) > maxaero) then
     print*,"Too many Salt-Film:",micro%salt_film_np(k,i,j)
     toomany=1
   endif
   if(micro%salt_jet_np(k,i,j)  > maxaero) then
     print*,"Too many Salt-Jet:",micro%salt_jet_np(k,i,j)
     toomany=1
   endif
   if(micro%salt_spum_np(k,i,j) > maxaero) then
     print*,"Too many Salt-Spume:",micro%salt_spum_np(k,i,j)
     toomany=1
   endif
  endif
  if(level==3 .and. iccnlev>=2)then
   if(micro%regen_aero1_np(k,i,j) > maxaero) then
     print*,"Too many Small-Regen-aero:",micro%regen_aero1_np(k,i,j)
     toomany=1
   endif
   if(micro%regen_aero2_np(k,i,j) > maxaero) then
     print*,"Too many Large-Regen-aero:",micro%regen_aero2_np(k,i,j)
     toomany=1
   endif
  endif
  if(toomany==1) then
    print*,"Aerosols exceed max value at k,i,j: ",k,i+mi0(ngr),j+mj0(ngr)
    print*,"This error typically occurs due to advection/diffusion"
    print*,"around too steep topography in the sigma-z coordinate."
    print*,"Solutions are: (1)reduce timestep, (2)use IHORGRAD=2,"
    print*,"or (3)smooth your topography in RAMSIN."
    stop
  endif

  !SET SMALL VALUES TO ZERO
  if(iaerosol>0) then
   if(micro%cccnp(k,i,j)<mincon .or. micro%cccmp(k,i,j)<minmas) then
      micro%cccnp(k,i,j) = 0.0
      micro%cccmp(k,i,j) = 0.0
   endif
   if(micro%gccnp(k,i,j)<mincon .or. micro%gccmp(k,i,j)<minmas) then
      micro%gccnp(k,i,j) = 0.0
      micro%gccmp(k,i,j) = 0.0
   endif
  endif
  if(idust>0) then
   if(micro%md1np(k,i,j)<mincon .or. micro%md1mp(k,i,j)<minmas) then
      micro%md1np(k,i,j) = 0.0
      micro%md1mp(k,i,j) = 0.0
   endif
   if(micro%md2np(k,i,j)<mincon .or. micro%md2mp(k,i,j)<minmas) then
      micro%md2np(k,i,j) = 0.0
      micro%md2mp(k,i,j) = 0.0
   endif
  endif
  if(iabcarb>0) then
   if(micro%abc1np(k,i,j)<mincon .or. micro%abc1mp(k,i,j)<minmas) then
      micro%abc1np(k,i,j) = 0.0
      micro%abc1mp(k,i,j) = 0.0
   endif
   if(micro%abc2np(k,i,j)<mincon .or. micro%abc2mp(k,i,j)<minmas) then
      micro%abc2np(k,i,j) = 0.0
      micro%abc2mp(k,i,j) = 0.0
   endif
  endif
  if(isalt>0) then
   if(micro%salt_film_np(k,i,j)<mincon .or. micro%salt_film_mp(k,i,j)<minmas) then
      micro%salt_film_np(k,i,j) = 0.0
      micro%salt_film_mp(k,i,j) = 0.0
   endif
   if(micro%salt_jet_np(k,i,j)<mincon .or. micro%salt_jet_mp(k,i,j)<minmas) then
      micro%salt_jet_np(k,i,j) = 0.0
      micro%salt_jet_mp(k,i,j) = 0.0
   endif
   if(micro%salt_spum_np(k,i,j)<mincon .or. micro%salt_spum_mp(k,i,j)<minmas) then
      micro%salt_spum_np(k,i,j) = 0.0
      micro%salt_spum_mp(k,i,j) = 0.0
   endif
  endif
  if(level==3 .and. iccnlev>=2) then
   if(micro%regen_aero1_np(k,i,j)<mincon .or. micro%regen_aero1_mp(k,i,j)<minmas)then
      micro%regen_aero1_np(k,i,j) = 0.0
      micro%regen_aero1_mp(k,i,j) = 0.0
   endif
   if(micro%regen_aero2_np(k,i,j)<mincon .or. micro%regen_aero2_mp(k,i,j)<minmas)then
      micro%regen_aero2_np(k,i,j) = 0.0
      micro%regen_aero2_mp(k,i,j) = 0.0
   endif
   if(itrkepsilon==1) then
     if(micro%resol_aero1_mp(k,i,j)<minmas) micro%resol_aero1_mp(k,i,j) = 0.0
     if(micro%resol_aero2_mp(k,i,j)<minmas) micro%resol_aero2_mp(k,i,j) = 0.0
     if(micro%resol_aero1_mp(k,i,j) > micro%regen_aero1_mp(k,i,j)) &
        micro%resol_aero1_mp(k,i,j) = 0.99*micro%regen_aero1_mp(k,i,j)
     if(micro%resol_aero2_mp(k,i,j) > micro%regen_aero2_mp(k,i,j)) &
        micro%resol_aero2_mp(k,i,j) = 0.99*micro%regen_aero2_mp(k,i,j)
   endif
  endif
 enddo

!  Think about how thp should change here - should it be due to a change in
!     rtp or to a change in the condensate?
!
!               vctr10(k) = rrp(k,i,j +rpp(k,i,j)  + rsp(k,i,j) + rap(k,i,j)
!     +                   + rgp(k,i,j) + rhp(k,i,j)
!               thp(k,i,j) = thp(k,i,j)
!     +                    * (1. - aklv * (vctr8(k) - rtp(k,i,j))
!c or +                    * (1. - aklv * (vctr10(k) - vctr9(k,i,j))
!     +                    /(max(temp, 253.)))

   enddo
enddo

return
END SUBROUTINE adj1

!##############################################################################
Subroutine ae1mic (ka,kb,c1,r3,r1)

use micphys, only:rxmin

implicit none

integer :: ka,kb
real, dimension(kb) :: c1,r3,r1
integer :: k

do k = ka,kb
   c1(k) = c1(k) * r1(k) / max(rxmin,r3(k))
   if (c1(k) < 0.) c1(k) = 0.
enddo

return
END SUBROUTINE ae1mic

!##############################################################################
Subroutine ae1kmic (ka,kb,cr3,cr1)

implicit none

integer :: ka,kb
real, dimension(kb) :: cr3,cr1
integer :: k

do k = ka,kb
   cr3(k) = cr1(k)
enddo

return
END SUBROUTINE ae1kmic
