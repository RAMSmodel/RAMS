!##############################################################################
Subroutine cuparm ()

use mem_tend
use mem_cuparm
use mem_basic
use mem_grid
use node_mod

implicit none

   ! Zero out arrays at model start unless doing a history-init
   if(nint(time) == 0 .and. initial /= 3) then
      CALL azero (mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
      CALL azero (mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
      CALL azero (mxp*myp,cuparm_g(ngrid)%conprr(1,1))
   endif

   if (mod(time + .001,confrq) .lt. dtlt .or. time .lt. 0.001) then

      if(iprntstmt>=1 .and. print_msg)print 90,time+dtlt,(time+dtlt)/3600.  &
               +(itime1/100+mod(itime1,100)/60.)
      90   format(' Kuo Convective Tendencies Updated Time =',f10.1,  &
               '  UTC TIME (HRS) =',f6.1)

      CALL azero (mxp*myp*mzp,cuparm_g(ngrid)%thsrc(1,1,1))
      CALL azero (mxp*myp*mzp,cuparm_g(ngrid)%rtsrc(1,1,1))
      CALL azero (mxp*myp,cuparm_g(ngrid)%conprr(1,1))

      CALL conpar (mzp,mxp,myp,ia,iz,ja,jz  &
          ,basic_g(ngrid)%up      (1,1,1)  &
          ,basic_g(ngrid)%vp      (1,1,1)  &
          ,basic_g(ngrid)%wp      (1,1,1)  &
          ,basic_g(ngrid)%theta   (1,1,1)  &
          ,basic_g(ngrid)%pp      (1,1,1)  &
          ,basic_g(ngrid)%pi0     (1,1,1)  &
          ,basic_g(ngrid)%dn0     (1,1,1)  &
          ,basic_g(ngrid)%rv      (1,1,1)  &
          ,cuparm_g(ngrid)%thsrc  (1,1,1)  &
          ,cuparm_g(ngrid)%rtsrc  (1,1,1)  &
          ,grid_g(ngrid)%rtgt     (1,1)    &
          ,cuparm_g(ngrid)%conprr (1,1))

   endif
   
CALL accum (mxp*myp*mzp,tend%tht(1),cuparm_g(ngrid)%thsrc(1,1,1)    )
CALL accum (mxp*myp*mzp,tend%rtt(1),cuparm_g(ngrid)%rtsrc(1,1,1)    )

CALL update (mxp*myp,cuparm_g(ngrid)%aconpr(1,1)  &
                   ,cuparm_g(ngrid)%conprr(1,1)  ,dtlt)

return
END SUBROUTINE cuparm

!##############################################################################
Subroutine conpar (m1,m2,m3,ia,iz,ja,jz  &
     ,up,vp,wp,theta,pp,pi0,dn0,rv,thsrc,rtsrc,rtgt,conprr)

use conv_coms
use mem_grid
use mem_cuparm
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3),theta(m1,m2,m3)  &
         ,pp(m1,m2,m3),pi0(m1,m2,m3),dn0(m1,m2,m3),rv(m1,m2,m3)  &
         ,thsrc(m1,m2,m3),rtsrc(m1,m2,m3),conprr(m2,m3)  &
         ,rtgt(m2,m3)

integer :: icpcnt=0,i1,i2,j1,j2,i,j,k,iprtfrq,iqmax,jqmax,kqmax
real :: dthmax

!
!        FLAG TO CONTROL PRINTOUT
!          ICPRTFL=0 - NO PRINTOUT
!                  1 - BRIEF INFO ON UP/DOWN DRAFT AND ITERATIONS
!                  2 - 1 PLUS MODEL TENDENCIES
!                  3 - 2 PLUS FINAL CONVECTIVE STRUCTURE
!                  4 - 3 PLUS UP/DOWN DRAFT AND ENVIRONMENT

icprtfl=0
iprtfrq=8
icpltfl=0
icpcnt=icpcnt+1
if(mod(icpcnt-iprtfrq+1,iprtfrq).eq.0) then
  icprtfl=1
endif

i1 = ia
i2 = iz
j1 = ja
j2 = jz

!  If variable initialization, on a coarse grid,
!  do not run convective parameterization in the lateral boundary region.

dthmax=0.

do j = j1,j2
   do i = i1,i2

      do k = 1,m1
         ucon(k)=up(k,i,j)
         vcon(k)=vp(k,i,j)
         wcon(k)=wp(k,i,j)
         thtcon(k)=theta(k,i,j)
         picon(k)=(pp(k,i,j)+pi0(k,i,j))
         tmpcon(k)=thtcon(k)*picon(k)/cp
         dncon(k)=dn0(k,i,j)
         prcon(k)=(picon(k)/cp)**cpor*p00
         rvcon(k)=rv(k,i,j)
         thsrc(k,i,j)=0.
         rtsrc(k,i,j)=0.
      enddo
    
      ! Compute heights AGL
      do k = 1, m1
         zcon(k)=zt(k) *rtgt(i,j)
         zzcon(k)=zm(k) *rtgt(i,j)
      enddo
         
    conprr(i,j)=0.
    wconmin=wcldbs
    contim=confrq

    CALL cu_environ (2,m1-1)
    if(igo.ne.0) CALL kuocp ()

    if(igo.ne.0) then

      CALL cp2mod (2,m1)
      do k=2,m1-1
        thsrc(k,i,j)=ftcon(k)
        rtsrc(k,i,j)=frcon(k)

      enddo
      conprr(i,j)=cprecip

      do k=2,m1-1
        if(thsrc(k,i,j).gt.dthmax) then
          dthmax=thsrc(k,i,j)
          iqmax=i
          jqmax=j
          kqmax=k
        endif
      enddo

! This will not work in parallel unless we convert I,J to domain coords
!      if(icprtfl.gt.0) then
!        print 899,i,j,time
!899         format(' * CONVECTION AT',2I4,'  TIME = ',F10.2)
!      endif

    endif

  enddo
enddo

!This will not work in parallel unless we run a domain max func
!and have nodes pass information to a head node and then print.
!if(iprntstmt>=1 .and. print_msg) print 675,iqmax,jqmax,kqmax,dthmax*86400.
!675 format(' Max Convective Heating RAte at I,J,K,D/DAY -',3I5,F8.2)

return
END SUBROUTINE conpar

!##############################################################################
Subroutine cu_environ (k1,k2)

use conv_coms
use rconstants

implicit none

integer :: k1,k2
real :: hz(nkp),wcpmax,themax,tlll,plll,rlll,zlll,dzlll,dzdd,abe &
       ,thdu,tdu,rdsu,znz
integer :: k,nkmid,nk

!       Basic constants
dzlow=200.
dzhigh=500.
zmid=3000.

cdzmin=3000.

!         Compute moist static energy profile

do k=k1,k2
  hz(k)=cp*tmpcon(k)+g*zcon(k)+alvl*rvcon(k)
enddo

!         Check for conditional instability and any upward motion
!           greater than WCONMIN under ZMID

igo=0
do k=k1,k2
  if(hz(k).gt.hz(k+1))then
    igo=1
    go to 105
  endif
enddo
105 continue
if(igo.eq.0)return

igo=0
wcpmax=-1.e10
do k=k1,k2
  if(zcon(k).gt.zmid)go to 104
  wcpmax=max(wcpmax,wcon(k))
enddo
104 continue
if(wcpmax.gt.0.0.and.wcpmax.gt.wconmin)igo=1
if(igo.eq.0)return

!           INTERPOLATE MODEL SOUNDING (ENVIRONMENT) TO HIGHER
!             RESOLUTION GRID

nkmid=zmid/dzlow+1
zc(1)=0.
do k=2,nkmid
  zc(k)=zc(k-1)+dzlow
enddo
do k=nkmid+1,nkp
  zc(k)=zc(k-1)+dzhigh
enddo
ze(1)=0.
do k=2,nkp
  ze(k)=(zc(k)+zc(k-1))*.5
enddo
!                   FIND MODEL TOP ON CONVECTIVE GRID
znz=zcon(k2)
do k=nkp,1,-1
  if(ze(k).lt.znz)go to 13
enddo
stop ' envir stop 12'
13 continue
kmt=k
!                   DO ACTUAL INTERPOLATION

nk=k2-k1+1
CALL htint (nk,  ucon(k1-1), zcon(k1-1),kmt,upe,ze)
CALL htint (nk,  vcon(k1-1), zcon(k1-1),kmt,vpe,ze)
CALL htint (nk,  wcon(k1-1),zzcon(k1-1),kmt,wpe,ze)
CALL htint (nk,thtcon(k1-1), zcon(k1-1),kmt,the,ze)
CALL htint (nk, rvcon(k1-1), zcon(k1-1),kmt,rve,ze)
do k=1,kmt
  rve(k)=max(rve(k),1e-8)
enddo

!         COMPUTE THETA V, THETA E, AND GET PRESSURE PROFILE

pke(1)=picon(1)
do k=1,kmt
  thve(k)=the(k)*(1.+.61*rve(k))
enddo
do k=2,kmt
  pke(k)=pke(k-1)-g*2.*(ze(k)-ze(k-1))  &
        /(thve(k)+thve(k-1))
enddo
do k=1,kmt
  te(k)=the(k)*pke(k)/cp
  pe(k)=(pke(k)/cp)**cpor*p00
  rhoe(k)=pe(k)/(rgas*te(k)*(1.+.61*rve(k)))
enddo
do k=1,kmt
  CALL thetae (pe(k),te(k),rve(k),thee(k))
enddo

!         FIND THE MAIN SOURCE LEVEL OF THE UPDRAFT.

!           FIRST TEST - ANY INVERSION BELOW 1.2 KM

do k=3,nkmid
  if(te(k).gt.te(k-1).and.te(k).gt.te(k+1)  &
                     .and.ze(k).le.1200.)then
  kcon=k
  go to 77
  endif
enddo

!           IF THERE ISN'T AN INVERSION, USE THE LEVEL OF HIGHEST
!           THETA E .

themax=0.
do k=2,nkmid
  if(thee(k).gt.themax)then
    themax=thee(k)
    kcon=k
  endif
enddo

!         FIND THE LCL OF A LAYER AVERAGE AROUND THE SOURCE LEVEL

77 continue
tlll=(te(kcon)+te(kcon+1)+te(kcon-1))/3.
plll=pe(kcon)
rlll=(rve(kcon)+rve(kcon+1)+rve(kcon-1))/3.
zlll=ze(kcon)

CALL lcl (tlll,plll,rlll,tlcl,plcl,dzlcl)

!         FIND THE CLOSEST LEVEL ON THE CONVECTIVE GRID TO THE LCL

dzlll=1e20
do k=1,kmt
  dzdd=abs(ze(k)-(zlll+dzlcl))
  if(dzdd.lt.dzlll)then
    dzlll=dzdd
    klcl=k
  endif
enddo

!         IF THERE IS NOT UPWARD MOTION AT THE LCL, NO CONVECTION
!           (MUST BE GREATER THAN WCONMIN )

if(wpe(klcl).lt.0.0.or.wpe(klcl).lt.wconmin)then
  igo=0
  return
endif

!         LOCATE EQUILIBRIUM TEMPERATURE LEVEL OF AN UNENTRAINED PARCEL.
!         COMPUTE INITIAL ABE.  IF ABE IS LESS THAN 0, NO CONVECTION.

theu(klcl)=the(kcon)*exp(alvl*rve(kcon)/(cp*tlcl))

do k=klcl,kmt
  if(theu(klcl).le.thve(k))go to 66
enddo
print*,'convection above model top:',klcl,theu(klcl),thve(k)
ketl=kmt-2
!cccccc      stop 65
66 continue
ketl=k
if(ze(ketl)-ze(klcl).lt.cdzmin)then
  igo=0
  return
endif

abe=0.
do k=klcl,ketl
  CALL the2t (theu(klcl),pe(k),thdu,tdu,rdsu)
  abe=abe+(thdu*(1.+.61*rdsu)-thve(k))/thve(k)*(zc(k)-zc(k-1))
enddo
if(abe.le.0.)then
  igo=0
  return
endif

!     if(icprtfl.gt.0)then
!       print 899
! 899   format(///,' * convection is activated * ')
!     endif

return
END SUBROUTINE cu_environ

!##############################################################################
Subroutine kuocp ()

use conv_coms
use rconstants

implicit none

integer :: k,idownd,klfs,kdiv,kdet,kover,kcoolh,kheat
real :: supplyw,anegl,apos,anegh,dddt,dzdiv,wtlfs,wtlcl,wtdiv,wtgnd &
       ,bkuo,zdetr,dzdet,vhint,vmint,vdint,avgmin,avtdiff,overmax,factr &
       ,heatmx,coolhi,c1 

!         Downdraft flag - 0 - no downdrafts
!                          1 - simple downdraft model
idownd=1

do k=1,nkp
  ftcon(k)=0.
  frcon(k)=0.
enddo

!         Compute vertical moisture convergence into the cloud layer.
!           Vertical flux out cloud top is assumed small.

supplyw=rhoe(klcl)*rve(klcl)*(wpe(klcl)+wpe(klcl-1))*.5

supply=supplyw

if(supply.le.0.) then
  igo=0
  return
endif

!         This is the cloud model.  Updraft is constant THETA e and
!           saturated with respect to water.  There is no ice.
!           Cloud top is one level above ETL.
!
!         THETA e of the updraft

theu(klcl)=the(kcon)*exp(alvl*rve(kcon)/(cp*tlcl))

!         Equilibrium Temperature Level of the source level air.

igo=0
do k=klcl,kmt
  CALL the2t (theu(klcl),pe(k),thu(k),tu(k),rsu(k))
  if(thu(k).gt.the(k).and.igo.eq.0) then
    igo=1
    klfc=k
  endif
  if(thu(k).le.the(k).and.igo.eq.1)go to 66
enddo
if(igo.eq.0) return
PRINT*,' Convection beyond model top - THup, THenv ',THU(KMT)  &
      ,THE(KMT)
k=kmt-1
66 continue
ketl=min(k,kmt)
kct=min(ketl+1,kmt)
CALL the2t (theu(klcl),pe(kct),thu(kct),tu(kct),rsu(kct))
do k=1,klcl-1
  thu(k)=the(k)
enddo

!         If the cloud is not at least CDZMIN deep or cloud top is
!           under 500 mb, no convection.

if(ze(ketl)-ze(klfc).lt.cdzmin.or.pe(kct).gt.50000.)then
  igo=0
  return
endif

!         Require the positive area be 50% greater than the negative
!           area below the LFC and  5% greater in total.

anegl=0.
do k=klcl,klfc-1
  anegl=anegl+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
apos=0.
do k=klfc,ketl-1
  apos=apos+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
anegh=0.
do k=ketl,kct
  anegh=anegh+(thu(k)-the(k))*(zc(k)-zc(k-1))
enddo
if(apos.lt.abs(anegl)*1.5.or.apos.lt.abs(anegl+anegh)*1.05) then
  igo=0
  return
endif

if(idownd.eq.1) then

!         The downdraft model - starts at THETA e minimum (LFS).
!             Downdraft is 2 degrees colder than
!             environment at cloud base increasing to 5 degrees
!             colder at the ground.


!         Find LFS as THETA e minimum

do k=kct,2,-1
  if(thee(k).lt.thee(k+1).and.thee(k).lt.thee(k-1))go to 11
enddo
k=2
11 continue
klfs=k
if(klfs.le.klcl)klfs=klcl+1
thd(klfs)=the(klfs)

!        Limit dd deficit at the ground to the maximum of positive
!          temperature difference of updraft if less than 2.5 degrees.

dddt=0.
do k=klcl,kct
  dddt=max(dddt,thu(k)-the(k))
enddo
if(dddt.gt.2.5) dddt=5.

thd(2)=the(2)-dddt
thd(klcl)=the(klcl)-dddt*.2
do k=klcl,klfs
  thd(k)=thd(klcl)+(thd(klfs)-thd(klcl))/(ze(klfs)-ze(klcl))  &
    *(ze(k)-ze(klcl))
enddo
do k=3,klcl-1
  thd(k)=thd(2)+(thd(klcl)-thd(2))/(ze(klcl)-ze(2))  &
    *(ze(k)-ze(2))
enddo

!         Now we need to weight the downdraft relative to the updraft.
!           Assume that the dd weight is zero at the LFS, 1/2 of
!           updraft at cloud base, and equal to the updraft at cloud
!           base at the ground.

kdiv=0
dzdiv=1e20
do k=1,kmt
  if(abs(ze(k)-800.).lt.dzdiv)then
    kdiv=k
    dzdiv=abs(ze(k)-800.)
  endif
enddo
kdiv=max(min(klcl,kdiv),2)
if(kdiv.eq.klcl) kdiv=klcl-1

do k=1,nkp
  wtd(k)=0.
enddo
wtlfs=0.
wtlcl=.1
wtdiv=.2
wtgnd=1.
do k=klcl+1,klfs
  wtd(k)=wtlcl+(wtlfs-wtlcl)/(ze(klfs)-ze(klcl))  &
    *(ze(k)-ze(klcl))
enddo
do k=kdiv,klcl
  wtd(k)=wtdiv+(wtlcl-wtdiv)/(ze(klcl)-ze(kdiv))  &
    *(ze(k)-ze(kdiv))
enddo
do k=2,kdiv-1
  wtd(k)=wtgnd+(wtdiv-wtgnd)/(ze(kdiv)-ze(2))  &
    *(ze(k)-ze(2))
enddo

else

  do k=1,nkp
    wtd(k)=0.
  enddo
  do k=2,klcl-1
    thu(k)=the(k)
  enddo

endif

!         Compute infamous b parameter.  Use Fritsch/Chappell's
!           precipitation efficiency.

envshr=sqrt((upe(kct)-upe(klfc))**2  &
           +(vpe(kct)-vpe(klfc))**2)  &
           /(ze(kct)-ze(klfc))*1e3
if(envshr.gt.1.35) then
  preff=1.591-.639*envshr+.0953*envshr**2-.00496*envshr**3
else
  preff=.9
endif
bkuo=1.-preff

!         Vertical profiles of convective heating and moistening

do k=2,kmt
  vheat(k)=0.
  vmois(k)=0.
  vmdry(k)=0.
enddo

!         Find the weighted THETA to use for the convection.

do k=2,kct
  thcon(k)=wtd(k)*thd(k)+(1.-wtd(k))*thu(k)
enddo

!         Heating profile is difference between convective THETAs and
!           environment.

do k=2,kct
  vheat(k)=thcon(k)-the(k)
enddo

!         Moisture profile is difference between vapor's of updraft and
!           environment in the cloud layer.  Below cloud base, air is
!           dried by SUPPLY.  Downdrafts are assumed to have no effect
!           on this.

kdet=0
zdetr=.66667*ze(kct)
dzdet=1000000.
do k=klcl,kct
  if(abs(ze(k)-zdetr).lt.dzdet)then
    dzdet=abs(ze(k)-zdetr)
    kdet=k
  endif
enddo

do k=kdet,kct
  vmois(k)=1.
enddo
!      do k=klcl,kct
!        vmois(k)=rsu(k)-rve(k)
!      enddo

do k=2,klcl-1
  vmdry(k)=rve(k)
enddo

vhint=0.
vmint=0.
vdint=0.
do k=2,kmt
  vhint=vhint+vheat(k)*(zc(k)-zc(k-1))
  vmint=vmint+vmois(k)*(zc(k)-zc(k-1))
  vdint=vdint+vmdry(k)*(zc(k)-zc(k-1))
enddo

!         If VHINT is less than 0, there is more negative area than
!           positive area.  No convection allowed.

if(vhint.le.0.) then
  igo=0
  return
endif

!         Also require that there is a minimum average
!           temperature difference between the updraft and environment
!           from the LFC to the ETL.  This eliminates the cases where
!           VHINT is very small and the heating and cooling rates get
!           astronomically large.

avgmin=.10
avtdiff=0.
do k=klfc,ketl-1
  avtdiff=avtdiff+(thcon(k)-the(k))
enddo
avtdiff=avtdiff/max(1,ketl-klfc)
if(avtdiff.lt.avgmin) then
  igo=0
  return
endif

!         Heating and moistening rates

3100 continue
do k=2,kmt
  ftcon(k)=alvl*preff*supply*vheat(k)  &
     /(pke(k)*rhoe(k)*vhint)
enddo
do k=klcl,kct
  frcon(k)=bkuo*supply*vmois(k)/(rhoe(k)*vmint)
enddo
do k=2,klcl-1
  frcon(k)=-supply*vmdry(k)/(rhoe(k)*vdint)
enddo

do k=klfc,ketl-1
  qvct1(k)=the(k)+contim*ftcon(k)
enddo
overmax=0.
do k=klfc,ketl-1
  if(qvct1(k)-thu(k).gt.overmax)then
    overmax=(qvct1(k)-thu(k))/(ftcon(k)*contim)
    kover=k
  endif
enddo

if(overmax.gt.0.) then
  factr=1.-overmax
  supply=factr*supply
!        if(icprtfl.ge.1)print*,' reducing supply ',kover,factr
!     +   ,qvct1(kover),thu(kover)
  go to 3100
endif

cprecip=preff*supply

if(icprtfl.gt.0) then
!        PRINT*,' ----------------------------------------------------'
!        PRINT 898,ZE(KCON),ZE(KLCL),ZE(KLFC),ZE(KETL),ZE(KCT)
! 898   FORMAT(' CLOUD LEVELS - SOURCE,LCL,LFC,ETL,TOP(KM) ',-3P,5F5.1)
!        PRINT 896,SUPPLY/SUPPLYW
! 896   FORMAT(' SUPPLIES ' ,F8.4)
!        PRINT 897,THEU(KLCL),RSU(KLCL)*1E3
!897   FORMAT(' CLOUD PROPERTIES -  THETA E, RS AT LCL',F6.1,F8.2)
  coolhi=100000.
  heatmx=-10000.
  kcoolh=0
  kheat=0
  do k=ketl,kct
    if(ftcon(k).lt.coolhi) then
      coolhi=ftcon(k)
      kcoolh=k
    endif
  enddo
  do k=klcl,kct
    if(ftcon(k).gt.heatmx) then
      heatmx=ftcon(k)
      kheat=k
    endif
  enddo
  C1=86400.
!        PRINT 905,PE(KHEAT),HEATMX*C1,PE(KCOOLH),COOLHI*C1
! 905   FORMAT(' MAX-MIN HEATING- P,(K/DAY)', 2(-2PF8.1,0PF7.1))
!        PRINT 906,PREFF,PREFF*SUPPLY*3600.*10.
!     +           ,(WPE(KLCL)+WPE(KLCL-1))*.5
!     +           ,RVE(KLCL)*1E3
!906   FORMAT(' PRECIPITATION - EFFICIENCY,RATE(MM/HR)',F5.2,F6.2,  &
!  '  LCL W,RV',F8.4,F7.2)
ENDIF


IF(ICPRTFL.GE.2) THEN
!      PRINT 95,(K,ZE(K),PE(K),TE(K),THE(K),THEE(K),RVE(K),UCON(K)
!     +   ,VCON(K),WPE(K),K=1,KMT)
!95 FORMAT(//' ENVIRONMENT-K,Z,P,TE,THE,THEE,RVE,UP,VP,WP'/,  &
!   (I3, -2P,F8.1,-3P,F8.2,0P,3F7.2,3P,F6.2, -2P,3F8.2))
!      PRINT96,(K,THE(K),THU(K)   ,FTCON(K)*86400.,RVE(K),RSU(K),
!     +       FRCON(K)*86400.,VHEAT(K)*(ZC(K)-ZC(K-1))/VHINT,
!     +                       VMOIS(K)*(ZC(K)-ZC(K-1))/VMINT,
!     +                       VMDRY(K)*(ZC(K)-ZC(K-1))/VDINT,
!     +  K=2,KMT)
!96 FORMAT(//' HEATING/MOISTENING'  &
! ,'-K,THE,THU,THSRC,RVE,RSU,RTSRC,HEAT%,MOIST%,DRY%'/,  &
!   (I3,0P,3F7.1,3P,3F7.1,0P,3F6.2))
ENDIF

return
END SUBROUTINE kuocp

!##############################################################################
Subroutine cp2mod (k1,k2)

use conv_coms
use mem_scratch
use rconstants

implicit none

integer :: k1,k2
real, external :: ssumvect
integer :: k
real :: tftc,tftm,tfrc,tfrm,ftres,frres

!        Compute integrated heating and moistening tendencies


do k=2,kmt
  qvct1(k)=rhoe(k)*ftcon(k)*pke(k)
  qvct2(k)=rhoe(k)*alvl*frcon(k)
  qvct3(k)=(zc(k)-zc(k-1))*qvct1(k)
  qvct4(k)=(zc(k)-zc(k-1))*qvct2(k)
enddo
tftc = ssumvect(kmt-1,qvct3(2),1)
tfrc = ssumvect(kmt-1,qvct4(2),1)

!         Transfer tendencies to model grid

CALL vertmap2 (qvct1,zc,kmt,vctr5,zzcon,k2)
CALL vertmap2 (qvct2,zc,kmt,vctr6,zzcon,k2)

do k=k1,k2
  vctr5(k)=vctr5(k)*(zzcon(k)-zzcon(k-1))
  vctr6(k)=vctr6(k)*(zzcon(k)-zzcon(k-1))
enddo


!         Make sure the transfer from the convective grid to the model
!           grid happened correctly.

tftm = ssumvect(k2-k1+1,vctr5(k1),1)
tfrm = ssumvect(k2-k1+1,vctr6(k1),1)
!
ftres=tftm-tftc
frres=tfrm-tfrc
if(abs(ftres) > .01*abs(tftc)) then
  print*,' energy error in grid tranfser in convective param.'
  print*,' tftm,tftc ',tftm,tftc
endif

!         Change energy tendencies to temperature and mixing ratio
!           tendencies.

do k=k1,k2
  ftcon(k)=vctr5(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*picon(k))
  frcon(k)=vctr6(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*alvl)

enddo

return
END SUBROUTINE cp2mod

!##############################################################################
Subroutine vertmap2 (datin,zin,n3in,datout,zout,n3out)

implicit none

integer :: n3in,n3out
real, dimension(n3in) :: datin,zin
real, dimension(n3out) :: datout,zout(n3out)

real, allocatable :: qvct(:),vctr(:)
integer :: k,l
real :: dzlft

!  This routine assumes that output vertical structure will be lower than
!   input vertical structure!!!!

!         Transfer quantity from input grid levels to output grid

allocate(qvct(n3in),vctr(n3out))

do k=1,n3out
   vctr(k)=0.
enddo

dzlft=0.
l=2
do k=2,n3out
   !     print *,'******************************* working on output layer ',k
   if(dzlft.ne.0.) then
      if(zin(l) .gt. zout(k)) then
         vctr(k)=vctr(k)+datin(l)*(zout(k)-zout(k-1))
         dzlft=zin(l)-zout(k)
          !  print*,'dzlft2 layer:',k,l,datin(l),dzlft
         go to 61
      else
         vctr(k)=vctr(k)+datin(l)*dzlft
         !   print*,'dzlft layer:',k,l,datin(l),dzlft
         l=l+1
         if(l > n3in) exit
         dzlft=0.
      endif
   endif
60   continue
   if(zin(l) <= zout(k)) then
      vctr(k)=vctr(k)+datin(l)*(zin(l)-zin(l-1))
      !   print*,'full layer:',k,l,datin(l),zin(L),zin(L-1)
      l=l+1
      if(l > n3in) exit
      dzlft=0.
      if(zin(l-1) == zout(k)) go to 61
      go to 60
   else
      vctr(k)=vctr(k)+datin(l)*(zout(k)-zin(l-1))
      !    print*,'part layer:',k,l,vctr(k),datin(l),zout(k),zin(l-1)
      dzlft=zin(l)-zout(k)
   endif
61    continue
      !  print *,'****************************** done with output layer ',k
enddo


!         Change energy tendencies to temperature and mixing ratio
!           tendencies.

do k=2,n3out
  datout(k)=vctr(k)/(zout(k)-zout(k-1))
enddo

deallocate(qvct,vctr)

return
END SUBROUTINE vertmap2

!##############################################################################
Subroutine lcl (t0,pp0,r0,tlcl,plcl,dzlcl)

use rconstants

implicit none

real :: t0,pp0,r0,tlcl,plcl,dzlcl
real, parameter :: cpg=102.45
integer :: nitt,ip
real :: p0k,pi0i,ttth0,ttd,dz,pki,pppi,ti,rvs
real, external :: tdewpt
real, external :: rsatmix

ip=0
11 continue

plcl=pp0
tlcl=t0
p0k=pp0**rocp
pi0i=p0k/p00k*cp
ttth0=t0*p00k/p0k
ttd=tdewpt(pp0,r0)
dz=cpg*(t0-ttd)
if(dz.le.0.)then
dzlcl=0.
return
endif
do 100 nitt=1,50
pki=pi0i-g*dz/(ttth0*(1.+.61*r0))
pppi=(pki/cp)**cpor*p00
ti=ttth0*pki/cp
rvs=rsatmix(pppi,ti)
if(abs(rvs-r0).lt..00003)go to 110
ttd=tdewpt(pppi,r0)
dz=dz+cp/g*(ti-ttd)
100 continue
print*, 'no converge in LCL:',t0,pp0,r0
ip=ip+1
if(ip==1)go to 11
stop 'LCL no convergence'

110 continue
plcl=pppi
tlcl=ti
dzlcl=dz

return
END SUBROUTINE lcl
