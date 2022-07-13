!##############################################################################
Subroutine obs_isen (m1,m2,tsnd,psnd,zsnd,rsnd,usndz,vsndz  &
                    ,zsndz,lpst,lzst  &
                    ,obsu,obsv,obsp,obss,obsr,no,chstid)
use isan_coms
use rconstants

implicit none

integer :: m1,m2,no
real :: tsnd(m1,m2),psnd(m1,m2),zsnd(m1,m2)  &
         ,rsnd(m1,m2),usndz(m1,m2)  &
         ,vsndz(m1,m2),zsndz(m1,m2)  &
         ,obsu(no,*),obsv(no,*),obsp(no,*),obss(no,*),obsr(no,*)
integer, dimension(m1) :: lpst,lzst
character(len=8), dimension(*) :: chstid

real,dimension(maxlev) :: pk,pth,zi
integer ns,lp,lz,k,lmbot,nglev,lmtop,l,lbchyd,lbc,lbcp
real :: th,wt,pkn,syo,po,tho

do ns=1,nsta
   lp=lpst(ns)
   lz=lzst(ns)

   do k=1,nisn
      obsp(ns,k)=1E30
      obsr(ns,k)=1E30
      obsu(ns,k)=1E30
      obsv(ns,k)=1E30
      obss(ns,k)=1E30
   enddo

   lmbot=0
   nglev=0
   do k=1,lp
      pk(k)=1.e30
      if(psnd(ns,k).lt.1.e29) pk(k)=psnd(ns,k)**rocp
      pth(k)=1.e30
      if(tsnd(ns,k).lt.1.e29.and.psnd(ns,k).lt.1.e29)  &
           pth(k)=tsnd(ns,k)*(p00/psnd(ns,k))**rocp
      if(pth(k).lt.1.e30.and.lmbot.eq.0) lmbot=k
      if(pth(k).lt.1.e29) nglev=nglev+1
   enddo
   lmtop=0
   do k=lp,1,-1
      if(pth(k).lt.1.e30.and.lmtop.eq.0) lmtop=k
   enddo
   if(nglev.le.1) goto 600

   l=2
   do k=1,nisn
      th=levth(k)

      if(th.lt.pth(lmbot).or.th.gt.pth(lmtop)) then
         obsp(ns,k)=1e30
         obsr(ns,k)=1e30
         goto 510
      endif

      511 if(th.ge.pth(l-1).and.th.lt.pth(l)) goto 515
      l=l+1
      if(l.gt.lp) goto 500
      goto 511

      515 continue
      wt=(th-pth(l-1))/(pth(l)-pth(l-1))
      pkn=(pk(l)-pk(l-1))*wt+pk(l-1)
      obsp(ns,k)=pkn**cpor
      if(rsnd(ns,l).lt.1e19.and.rsnd(ns,l-1).lt.1e19) then
         obsr(ns,k)=rsnd(ns,l-1)+(rsnd(ns,l)-rsnd(ns,l-1))*wt
      else
         obsr(ns,k)=1e30
      endif
      510 continue
   enddo
   500 continue

   ! Find a boundary condition as first level at or below 360K

   lbchyd=360
   do k=nisn,1,-1
      !print*,'mjb p',k,ns,levth(k),lbchyd,obsp(ns,k)
      if(levth(k).le.lbchyd.and.obsp(ns,k).lt.1e19) then
         lbc=k
         goto 52
      endif
   enddo
   print*,'---Could not find 1st obs isentropic boundary level ',ns,' ',chstid(ns)
   !stop 'obs_isen'
   !MJW Seems like just set everything to missing and continue rather than stop
   do k=1,nisn
      obss(ns,k)=1.e30
      obsu(ns,k)=1.e30
      obsv(ns,k)=1.e30
   enddo
   cycle

   52 continue

   do k=lp,1,-1
      if(pth(k).le.float(levth(lbc))) then
         lbcp=k
         goto 51
      endif
   enddo
   print*,'---Could not find 2nd obs isentropic boundary level ',ns,' ',chstid(ns)
   !stop 'obs_isen'
   !MJW Seems like just set everything to missing and continue rather than stop
   do k=1,nisn
      obss(ns,k)=1.e30
      obsu(ns,k)=1.e30
      obsv(ns,k)=1.e30
   enddo
   cycle 

   51 continue
   syo=cp*pth(lbcp)*pk(lbcp)/p00**rocp+g*zsnd(ns,lbcp)
   obss(ns,lbc)=syo+cp*(pk(lbcp)+obsp(ns,lbc)**rocp)  &
        *.5/p00**rocp *(levth(lbc)-pth(lbcp))
   po=obsp(ns,lbc)
   syo=obss(ns,lbc)
   tho=levth(lbc)
   do k=lbc+1,nisn
      obss(ns,k)=1e30
      if(obsp(ns,k).lt.1e19) then
         obss(ns,k)=syo+cp*(po**rocp+obsp(ns,k)**rocp)  &
              *.5/p00**rocp *(levth(k)-tho)
         syo=obss(ns,k)
         po=obsp(ns,k)
         tho=levth(k)
      endif
   enddo

   syo=obss(ns,lbc)
   po=obsp(ns,lbc)
   tho=levth(lbc)
   do k=lbc-1,1,-1
      obss(ns,k)=1e30
      if(obsp(ns,k).lt.1e19) then
         obss(ns,k)=syo+cp*(po**rocp+obsp(ns,k)**rocp)*.5  &
              /p00**rocp*(levth(k)-tho)
         syo=obss(ns,k)
         po=obsp(ns,k)
         tho=levth(k)
      endif
   enddo

   do k=1,nisn
      if(obss(ns,k).lt.1e19) then
         zi(k)=(obss(ns,k)-cp*levth(k)*(obsp(ns,k)*p00i)**rocp)/g
      else
         zi(k)=1e30
      endif
   enddo

   l=2
   do k=1,nisn

      if(lz==0 .or. zi(k) < zsndz(ns,1) .or. zi(k) > zsndz(ns,lz)) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle
      endif

      611 if(zi(k).ge.zsndz(ns,l-1).and.zi(k).lt.zsndz(ns,l)) goto 615
      l=l+1
      if(l.gt.lz) goto 600
      goto 611
      615 continue
      wt=(zi(k)-zsndz(ns,l-1))/(zsndz(ns,l)-zsndz(ns,l-1))
      if(usndz(ns,l).lt.1e19.and.usndz(ns,l-1).lt.1e19  &
           .and.vsndz(ns,l).lt.1e19.and.vsndz(ns,l-1).lt.1e19) then
         obsu(ns,k)=usndz(ns,l-1)+(usndz(ns,l)-usndz(ns,l-1))*wt
         obsv(ns,k)=vsndz(ns,l-1)+(vsndz(ns,l)-vsndz(ns,l-1))*wt
      else
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
      endif
      
   enddo

   600 continue

enddo

return
END SUBROUTINE obs_isen

!##############################################################################
Subroutine obs_sigz (m1,m2,tsnd,psnd,zsnd,rsnd,usndz,vsndz  &
                    ,zsndz,lpst,lzst,sndtopg  &
                    ,obsu,obsv,obsp,obst,obsr,no,ztop)
use isan_coms
use rconstants

implicit none

integer :: m1,m2,no,lpst(m1),lzst(m1)
real ::   tsnd(m1,m2),psnd(m1,m2),zsnd(m1,m2)  &
         ,rsnd(m1,m2),usndz(m1,m2),vsndz(m1,m2),zsndz(m1,m2),sndtopg(m1)  &
         ,obsu(no,*),obsv(no,*),obsp(no,*),obst(no,*),obsr(no,*)
         
real :: pk(maxlev),pth(maxlev),sigzr(maxsigz)  &
         ,tsz(maxsigz),rsz(maxsigz),psz(maxsigz)

integer :: ns,k,lmbot,lp,lz,nglev,lmtop,l,lbc,lbcp,kbcu
real :: ztop,wt,bchyd,pio,zso,tho
real, external :: rsatmix


do ns=1,nsta
   lp=lpst(ns)
   lz=lzst(ns)
   
   !DO K=1,nsigz
   !   print*,'mjb a5',ns,nsigz,k,zsnd(ns,k)
   !enddo

   ! Fill sigzr from interpolated grid topography height

   DO K=1,nsigz
      sigzr(K)=sndtopg(ns)+sigz(k)*(1.-sndtopg(ns)/ztop)
   ENDDO

   do k=1,nsigz
      obsp(ns,k)=1E30
      obsr(ns,k)=1E30
      obst(ns,k)=1E30
   enddo

   lmbot=0
   nglev=0
   do k=1,lp
      pk(k)=1.e30
      if(psnd(ns,k).lt.1.e29) pk(k)=psnd(ns,k)**rocp
      pth(k)=1.e30
      if(tsnd(ns,k).lt.1.e29.and.psnd(ns,k).lt.1.e29)  &
         pth(k)=tsnd(ns,k)*(p00/psnd(ns,k))**rocp
      if(pth(k).lt.1.e30.and.lmbot.eq.0) lmbot=k
      if(pth(k).lt.1.e29) nglev=nglev+1
   enddo
   lmtop=0
   do k=lp,1,-1
      if(pth(k).lt.1.e30.and.lmtop.eq.0) lmtop=k
   enddo
   if(nglev.le.1) goto 1000

   l=2

   !print*,'mjb a',ns,zsnd(ns,1)
   do k=1,nsigz
      if(sigzr(k).lt.zsnd(ns,1)) then
         obst(ns,k)=pth(1)
         obsr(ns,k)=rsnd(ns,1)
         goto 510
      endif
      if(sigzr(k).gt.zsnd(ns,lp)) then
         obst(ns,k)=1e30
         obsr(ns,k)=1e30
         goto 510
      endif

      511 if(sigzr(k).ge.zsnd(ns,l-1).and.sigzr(k).lt.zsnd(ns,l)) goto 515
      l=l+1
      if(l.gt.lp) goto 500
      goto 511

      515 continue
      wt=(sigzr(k)-zsnd(ns,l-1))/(zsnd(ns,l)-zsnd(ns,l-1))
      obst(ns,k)=(pth(l)-pth(l-1))*wt+pth(l-1)

      if(rsnd(ns,l).lt.1e19.and.rsnd(ns,l-1).lt.1e19) then
         obsr(ns,k)=rsnd(ns,l-1)+(rsnd(ns,l)-rsnd(ns,l-1))*wt
      else
         obsr(ns,k)=1e30
      endif

      510 continue
      !print*,'mjb b',k,ns,lp,obst(ns,k),obsr(ns,k)
   enddo
   500 continue

   ! Find a boundary condition as first level at or below 10000m

   bchyd=10000.
   do k=nsigz,1,-1
      !print*,'mjb c',ns,k,sigzr(k),bchyd,obst(ns,k)
      if(sigzr(k).le.bchyd.and.obst(ns,k).lt.1e19) then
         lbc=k
         goto 52
      endif
   enddo
   print*,'---Could not find 1st obs sigma-z boundary level',ns
   !stop 'obs_sigz'
   !MJW Seems like just set everything to missing and continue rather than stop
   do k=1,nsigz
      obsp(ns,k)=1.e30
      obst(ns,k)=1.e30
      obsu(ns,k)=1.e30
      obsv(ns,k)=1.e30
   enddo
   cycle
   52 continue

   do k=lp,1,-1
      if(zsnd(ns,k).le.sigzr(lbc)) then
         lbcp=k
         goto 51
      endif
   enddo
   print*,'---Could not find 2nd obs sigma-z boundary level',ns
   !stop 'obs_sigz'
   !MJW Seems like just set everything to missing and continue rather than stop
   do k=1,nsigz
      obsp(ns,k)=1.e30
      obst(ns,k)=1.e30
      obsu(ns,k)=1.e30
      obsv(ns,k)=1.e30
   enddo
   cycle
   51 continue

   pio=cp*(psnd(ns,lbcp)/p00)**rocp
   obsp(ns,lbc)=pio-(sigzr(lbc)-zsnd(ns,lbcp))*g/((pth(lbcp)+obst(ns,lbc))*.5)
   pio=obsp(ns,lbc)
   zso=sigzr(lbc)
   tho=obst(ns,lbc)
   do k=lbc+1,nsigz
      obsp(ns,k)=1e30
      if(obst(ns,k).lt.1e19) then
         obsp(ns,k)=pio-(sigzr(k)-zso)*g/((obst(ns,k)+tho)*.5)
         zso=sigzr(k)
         pio=obsp(ns,k)
         tho=obst(ns,k)
      endif
   enddo

   zso=sigzr(lbc)
   pio=obsp(ns,lbc)
   tho=obst(ns,lbc)
   do k=lbc-1,1,-1
      obsp(ns,k)=1e30
      if(obst(ns,k).lt.1e19) then
         obsp(ns,k)=pio-(sigzr(k)-zso)*g/((obst(ns,k)+tho)*.5)
         zso=sigzr(k)
         pio=obsp(ns,k)
         tho=obst(ns,k)
      endif
   enddo

   ! Compute virtual temp

   DO K=1,nsigz
      IF(obsp(ns,K)+obst(ns,k).LT.1E19) THEN
         tsz(k)=obst(ns,k)*obsp(ns,k)/cp
         psz(K)=(obsp(ns,k)/cp)**cpor*p00
         IF(obsr(ns,k).LT.1E19) THEN
            rsz(k)=rsatmix(psz(k),tsz(k)) * obsr(ns,k)
            tsz(k)=obst(ns,k)*(1.+.61*rsz(k))
         else
            tsz(k)=obst(ns,k)
         endif
      endif
   enddo

   ! Recompute pressure profile with virtual temp

   do k=nsigz,1,-1
      if(obsp(ns,k).lt.1e19) then
         kbcu=k
         goto 551
      endif
   enddo
   print*, '---Could not find good obs sigma-z boundary level 3',ns
   !stop 'obs_sigz'
   !MJW Seems like just set everything to missing and continue rather than stop
   do k=1,nsigz
      obsp(ns,k)=1.e30
      obst(ns,k)=1.e30
      obsu(ns,k)=1.e30
      obsv(ns,k)=1.e30
   enddo
   cycle
   551 continue

   zso=sigzr(kbcu)
   pio=obsp(ns,kbcu)
   tho=tsz(kbcu)
   do k=kbcu-1,1,-1
      obsp(ns,k)=1e30
      if(obst(ns,k).lt.1e19) then
         obsp(ns,k)=pio-(sigzr(k)-zso)*g/((tsz(k)+tho)*.5)
         zso=sigzr(k)
         pio=obsp(ns,k)
         tho=tsz(k)
      endif
   enddo
   do k=1,nsigz
      if(obsp(ns,k).lt.1e19) obsp(ns,k)=(obsp(ns,k)/cp)**cpor*p00
   enddo

   ! Vertically interpolate winds in height

   1000 continue

   l=2
   do k=1,nsigz

      if(lz==0 .or. sigzr(k) < zsndz(ns,1) .or. sigzr(k) > zsndz(ns,lz)) then
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
         cycle
      endif

      611 if(sigzr(k).ge.zsndz(ns,l-1).and.sigzr(k).lt.zsndz(ns,l)) goto 615
      l=l+1
      if(l.gt.lz) goto 600
      goto 611
      615 continue
      wt=(sigzr(k)-zsndz(ns,l-1))/(zsndz(ns,l)-zsndz(ns,l-1))
      if(usndz(ns,l).lt.1e19.and.usndz(ns,l-1).lt.1e19  &
           .and.vsndz(ns,l).lt.1e19.and.vsndz(ns,l-1).lt.1e19) then
         obsu(ns,k)=usndz(ns,l-1)+(usndz(ns,l)-usndz(ns,l-1))*wt
         obsv(ns,k)=vsndz(ns,l-1)+(vsndz(ns,l)-vsndz(ns,l-1))*wt
      else
         obsu(ns,k)=1e30
         obsv(ns,k)=1e30
      endif
      
   enddo

   600 continue

enddo

return
END SUBROUTINE obs_sigz

!##############################################################################
Subroutine vterpp_i (np1,np2,np3,npi3,un,vn,tn,zn,rn,ui2,vi2,pi2,si2,ri2)
     
use isan_coms
use rconstants

implicit none

integer :: np1,np2,np3,npi3
real :: un(np1,np2,np3),vn(np1,np2,np3),tn(np1,np2,np3)  &
         ,zn(np1,np2,np3),rn(np1,np2,np3)  &
         ,ui2(np1,np2,npi3),vi2(np1,np2,npi3),pi2(np1,np2,npi3)  &
         ,si2(np1,np2,npi3),ri2(np1,np2,npi3)

integer, parameter :: npr=maxpr+2
real :: ppd(npr),thd(npr),pkd(npr),ud(npr),vd(npr),rd(npr)

integer :: i,j,k,mcnt,npd,lp,kpbc,kibc
real :: thl,wt,pkn,pbc,sy,dtemp
real, external :: rsatmix

do j=1,np2
   do i=1,np1

      DO K=1,NISN
         UI2(I,J,K)=1.E30
         VI2(I,J,K)=1.E30
         RI2(I,J,K)=1.E30
         PI2(I,J,K)=1.E30
         SI2(I,J,K)=1.E30
      ENDDO

      ! determine if this column is all missing. if so, just leave
      !   isentropic data as missing

      mcnt=0
      DO K=1,NPRZ
         if(tn(i,j,k).gt.1000.) mcnt=mcnt+1
      ENDDO
      if(mcnt.eq.nprz) goto 4500

      DO K=1,NPRZ
         pnpr(k)=levpr(k)*100.                  !pressure(Pa)
         UD(K+2)=UN(I,J,K)                      !u-wind(m/s)
         VD(K+2)=VN(I,J,K)                      !v-wind(m/s)
         RD(K+2)=RN(I,J,K)                      !RH-wind(m/s)
         PPD(K+2)=PNPR(K)                       !pressure(Pa)
         THD(K+2)=TN(I,J,K)*(P00/PNPR(K))**ROCP !Theta(K)
      ENDDO

      ! Define two phony levels for isentropes underground

      UD(1)=UD(3)
      UD(2)=UD(3)
      VD(1)=VD(3)
      VD(2)=VD(3)
      RD(1)=RD(3)
      RD(2)=RD(3)
      PPD(1)=120000.                             
      PPD(2)=110000. 
      THD(2)=(TN(I,J,1)-2.25)*(P00/PNPR(1))**ROCP
      THD(1)=180.

      !Saleeby:(2020) Experimenting Here:
      !Try to make better extrapolation below 1000mb level
      !Extrapolate RH(2) rather than assume same moisture (keep between 10%-95% RH)
      !RD(2) = ( (RD(3)-RD(4))/(PPD(3)-PPD(4)) * (PPD(2)-PPD(3)) ) + RD(3)
      !RD(2) = min(0.95,max(0.10,RD(2)))
      !RD(1) = RD(2)
      !Saleeby:(2020) Extrapolate TN,THD(2) rather than assume an arbitrary lapse rate
      !dtemp = (TN(I,J,1)-TN(I,J,2))/(pnpr(1)-pnpr(2)) * (PPD(2)-PPD(3))
      !THD(2)=(TN(I,J,1)-dtemp)*(P00/PPD(2))**ROCP

      npd=nprz+2
      DO K=1,NPD !Up to pressure levels +2
         PKD(K)=PPD(K)**ROCP !Equals p00(T/Theta)
      ENDDO

      !Check some output here if needed
      !if(i==74.and.j==135)then
      ! DO K=1,NPRZ
      !   write(*,'(a,2(I8,1X),4(F10.3,1X))') 'plevs' &
      !    ,K,levpr(K),ZN(I,J,K),RN(I,J,K),TN(I,J,K) &
      !    ,1000.*RN(I,J,K)*rsatmix(PNPR(K),TN(I,J,K))
      ! ENDDO
      !endif

      lp = npd

      do k = nisn,1,-1 !Compute from top to bottom
         35 continue
         thl = levth(k)

         if (thl .gt. thd(npd)) then

            pi2(i,j,k) = (thd(npd) / thl) ** cpor * ppd(npd)
            ui2(i,j,k) = ud(npd)
            vi2(i,j,k) = vd(npd)
            ri2(i,j,k) = 0.

         elseif (thl .le. thd(lp) .and. thl .gt. thd(lp-1)) then

            wt = (thl - thd(lp-1)) / (thd(lp) - thd(lp-1))    !interpolation weights
            pkn = pkd(lp-1) + (pkd(lp) - pkd(lp-1)) * wt      !p00*(T/Theta)
            pi2(i,j,k) = pkn ** cpor                          !pressure on isen surface
            ui2(i,j,k) = ud(lp-1) + (ud(lp) - ud(lp-1)) * wt  !u-wind
            vi2(i,j,k) = vd(lp-1) + (vd(lp) - vd(lp-1)) * wt  !v-wind
            ri2(i,j,k) = rd(lp-1) + (rd(lp) - rd(lp-1)) * wt  !RH on isen surface

         else
            lp = lp - 1
            if (lp .le. 1) then
               print*, 'vterpp_i interpolation tried to go below'
               print*, 'lowest pressure level.'
               stop 'vterpp_i'
            endif
            goto 35

         endif

      enddo

      ! *****************************************************
      ! Find a pressure b.c. as second isentrope above 400 mb
      ! *****************************************************

      pbc=40000.
      kpbc=1
      do k=1,npd
         if(ppd(k).lt.pbc) then
            kpbc=k
            exit
         endif
      enddo

      kibc=2
      do k=2,nisn
         if(pi2(i,j,k).lt.ppd(kpbc)) then
            kibc=k
            exit
         endif
      enddo

      !Montgomery streamfunction (PSI or "sy") value at level "kpbc"
      sy=cp*thd(kpbc)*pkd(kpbc)/p00k+g*zn(i,j,kpbc-2)

      !Compute streamfunction is kibc-1 isentropic level
      si2(i,j,kibc-1)=sy-cp*(pi2(i,j,kibc-1)**rocp  &
           +ppd(kpbc)**rocp)/(2.*p00k)*(thd(kpbc)-levth(kibc-1))
      !Integrate streamfunction downward based on level above it using hydrostatics
      do k=kibc-2,1,-1
         si2(i,j,k)=si2(i,j,k+1)+cp*(pi2(i,j,k+1)**rocp  &
              +pi2(i,j,k)**rocp)/(2.*p00k)*(levth(k)-levth(k+1))
      enddo

      !Compute streamfunction is kibc isentropic level
      si2(i,j,kibc)=sy+cp*(pi2(i,j,kibc)**rocp  &
           +ppd(kpbc)**rocp)/(2.*p00k)*(levth(kibc)-thd(kpbc))
      !Integrate streamfunction upward based on level below it using hydrostatics
      do k=kibc+1,nisn
         si2(i,j,k)=si2(i,j,k-1)+cp*(pi2(i,j,k-1)**rocp  &
              +pi2(i,j,k)**rocp)/(2.*p00k)*(levth(k)-levth(k-1))
      enddo

      4500 continue

      !Check some output here if needed
      !if(i==74.and.j==135)then
      ! do k = nisn,1,-1 !Compute from top to bottom
      !  write(*,'(a,2(I4,1X),4(F8.3,1X),2(F11.3,1X))') 'vterpi' &
      !    ,k,lp,thd(lp-1),thl,thd(lp),ri2(i,j,k),pi2(i,j,k),si2(i,j,k)
      ! enddo
      !endif

   enddo
enddo

! If any level is above 100 mb, compute geostrophic winds

!GDKM=2.*SPCON
!DO J=2,NP2-1
!   GDLAT=(XSWLAT+(J-1)*GDATDY)*PI180
!   FCORI=1./(2.*7.292E-5*SIN(GDLAT))
!   DO I=2,NP1-1
!      DO K=NISN,1,-1
!         IF(PI2(I,J,K).LT.10000.) THEN
!            UI2(I,J,K)=-FCORI*(SI2(I,J+1,K)-SI2(I,J-1,K))/(GDKM*GDATDY)
!            VI2(I,J,K)= FCORI*(SI2(I+1,J,K)-SI2(I-1,J,K))  &
!                 /(GDKM*GDATDX*COS(GDLAT))
!         ENDIF
!      ENDDO
!   ENDDO
!ENDDO

return
END SUBROUTINE vterpp_i

!##############################################################################
Subroutine vterpp_s (np1,np2,np3,npi3,un,vn,tn,zn,rn  &
                    ,ui2,vi2,pi2,ti2,ri2,topt,rtgt)
     
use isan_coms
use rconstants

implicit none

integer :: np1,np2,np3,npi3
real :: un(np1,np2,np3),vn(np1,np2,np3),tn(np1,np2,np3)  &
         ,zn(np1,np2,np3),rn(np1,np2,np3)  &
         ,ui2(np1,np2,npi3),vi2(np1,np2,npi3),pi2(np1,np2,npi3)  &
         ,ti2(np1,np2,npi3),ri2(np1,np2,npi3)  &
         ,topt(np1,np2),rtgt(np1,np2)

integer, parameter :: npr=maxpr+2
real :: ppd(npr),thetd(npr),pkd(npr),ud(npr),vd(npr),zd(npr)  &
         ,rd(npr),pid(npr),tempd(npr),rtd(npr),thvd(npr)  &
         ,sigzr(maxsigz),vvv(maxsigz)

integer :: i,j,k,mcnt,npd,kl,kpbc,kibc
real :: pbc,thvp,pii,dtemp
real, external :: rsatmix

do j=1,np2
   do i=1,np1

      DO K=1,npi3
         UI2(I,J,K)=1.E30
         VI2(I,J,K)=1.E30
         RI2(I,J,K)=1.E30
         PI2(I,J,K)=1.E30
         TI2(I,J,K)=1.E30
      ENDDO

      ! determine if this column is all missing.
      ! if so, just leave data as missing

      mcnt=0
      DO K=1,NPRZ
         if(tn(i,j,k).gt.1000.) mcnt=mcnt+1
      ENDDO
      if(mcnt.eq.nprz) goto 4500

      DO K=1,NPRZ
         PNPR(K)=levpr(k)*100.
         UD(K+2)=UN(I,J,K)
         VD(K+2)=VN(I,J,K)
         RD(K+2)=RN(I,J,K)
         PPD(K+2)=PNPR(K)
         THETD(K+2)=TN(I,J,K)*(P00/PNPR(K))**ROCP
         ZD(K+2)=ZN(I,J,K)
      ENDDO

      ! Define two phony levels for isentropes underground

      UD(1)=UD(3)
      UD(2)=UD(3)
      VD(1)=VD(3)
      VD(2)=VD(3)
      RD(1)=RD(3)
      RD(2)=RD(3)
      PPD(1)=120000.
      PPD(2)=110000.
      THETD(2)=(TN(I,J,1)-2.25)*(P00/PNPR(1))**ROCP
      THETD(1)=220.

      !Saleeby:(2020) Experimenting here:
      !Try to make better extrapolation below 1000mb level
      !Extrapolate RH(2) rather than assume same moisture (keep between 10%-95% RH)
      !RD(2) = ( (RD(3)-RD(4))/(PPD(3)-PPD(4)) * (PPD(2)-PPD(3)) ) + RD(3)
      !RD(2) = min(0.95,max(0.10,RD(2)))
      !RD(1) = RD(2)
      !Saleeby:(2020) Extrapolate TN,THD(2) rather than assume an arbitrary lapse rate
      !dtemp = (TN(I,J,1)-TN(I,J,2))/(pnpr(1)-pnpr(2)) * (PPD(2)-PPD(3))
      !THETD(2)=(TN(I,J,1)-dtemp)*(P00/PPD(2))**ROCP

      npd=nprz+2
      DO K=1,NPD
         PKD(K)=PPD(K)**ROCP                    !p00*(T/Theta)
         pid(k)=cp*(ppd(k)/p00)**rocp           !PI
         tempd(k)=thetd(k)*pid(k)/cp            !Temperature
         rtd(k)=rd(k)*rsatmix(ppd(k),tempd(k))  !vapor mixing ratio
         thvd(k)=thetd(k)*(1.+.61*rtd(k))       !Theta-V
      ENDDO

      !Compute height of fake levels using Theta-V and PI
      !Recall: PI = cp * (p/p0)**(R/cp).
      !This equation comes from natural-log laws and 
      !geo height / hypsometric equations.
      zd(2)=zd(3)+(thvd(3)+thvd(2))*.5*(pid(3)-pid(2))/g
      zd(1)=zd(2)+(thvd(2)+thvd(1))*.5*(pid(2)-pid(1))/g
      do k=1,npi3
         sigzr(k)=topt(i,j)+sigz(k)*rtgt(i,j)
      enddo 

      !Vert Interp vars from pressure geoheight levs to RAMS sigma levs
      CALL htint (npd,ud,zd,npi3,vvv,sigzr)    !for u-wind
      CALL psfill (npi3,vvv,ui2,np1,np2,i,j)   !u-wind to ui2 (U-wind)
      CALL htint (npd,vd,zd,npi3,vvv,sigzr)    !for v-wind
      CALL psfill (npi3,vvv,vi2,np1,np2,i,j)   !v-wind to vi2 (V-wind)
      CALL htint (npd,thetd,zd,npi3,vvv,sigzr) !for thetd
      CALL psfill (npi3,vvv,ti2,np1,np2,i,j)   !thetd to ti2  (potential temp)
      CALL htint (npd,rd,zd,npi3,vvv,sigzr)    !for rd
      CALL psfill (npi3,vvv,ri2,np1,np2,i,j)   !rd to ri2     (RH)
      CALL htint (npd,pid,zd,npi3,vvv,sigzr)   !for pid 
      CALL psfill (npi3,vvv,pi2,np1,np2,i,j)   !pid to pi2    (PI)
      do k=1,npi3
         pi2(i,j,k)=(pi2(i,j,k)/cp)**cpor*p00  !convert PI back to Pressure
         !Check some output here if needed
         !if(i==74.and.j==135)then
         ! write(*,'(a,1(I4,1X),4(F11.3,1X))') &
         !   'vterps1',k,sigzr(k),ri2(i,j,k),ti2(i,j,k),pi2(i,j,k)
         !endif
      enddo

      DO KL=Npi3,1,-1
         IF(TI2(I,J,KL).LT.1E19) GOTO 402
      ENDDO
      STOP 'ST2-MISS'
      402 CONTINUE
      KL=KL+1
      DO K=KL,Npi3
         TI2(I,J,K)=tempd(npd)*(p00/pi2(i,j,k))**rocp
         UI2(I,J,K)=UD(NPD)
         VI2(I,J,K)=VD(NPD)
         RI2(I,J,K)=0.
      ENDDO

      ! *************************************************
      ! Find a pressure b.c. as second level above 700 mb
      ! *************************************************
      pbc=70000.
      kpbc=1
      DO K=1,npd
         if(ppd(k).lt.pbc) then
            kpbc=k
            goto 320
         endif
      ENDDO
      320 continue

      kibc=2
      DO K=2,npi3
         if(ti2(i,j,k).gt.thetd(kpbc)) then
            kibc=k
            goto 321
         endif
      ENDDO
      print*,'ISAN error: domain top not high enough'
      stop 'vterpp_s'
      321 continue

      do k=1,npi3
         !Convert theta to temperature
         vvv(k)=ti2(i,j,k)*(pi2(i,j,k)/p00)**rocp
         !Compute theta-v
         vvv(k)=ti2(i,j,k)*(1.+.61*ri2(i,j,k)*rsatmix(pi2(i,j,k),vvv(k)))
      enddo
      !Theta-v at the pressure B.C. defined above
      thvp=thetd(kpbc)*(1.+.61*rd(kpbc)*rsatmix(ppd(kpbc),tempd(kpbc)))

      !PI at pressure B.C. level
      pii=cp*pkd(kpbc)/p00**rocp
      !PI at sigma-z level interpolated from B.C. pressure level but 
      ! based on theta-v
      pi2(i,j,kibc-1)=pii+(zd(kpbc)-sigzr(kibc-1))*g/(.5*(thvp+vvv(kibc-1)))
      !Integrate PI downward hydrostatically
      do k=kibc-2,1,-1
         pi2(i,j,k)=pi2(i,j,k+1)+(sigzr(k+1)-sigzr(k))*g/(.5*(vvv(k+1)+vvv(k)))
      enddo
      !Integrate PI upward hydrostatically
      do k=kibc,npi3
         pi2(i,j,k)=pi2(i,j,k-1)-(sigzr(k)-sigzr(k-1))*g/(.5*(vvv(k-1)+vvv(k)))
      enddo
      !Convert PI back to Pressure but now based on Theta-V and extrapolation
      !Pressure was first interpolated to sigma levels and then put in hydrostatic
      !balance in the columns.
      do k=1,npi3
         pi2(i,j,k)=(pi2(i,j,k)/cp)**cpor*p00
         !if(i==74.and.j==135)then
         ! write(*,'(a,1(I4,1X),4(F11.3,1X))') &
         !   'vterps2',k,sigzr(k),ri2(i,j,k),ti2(i,j,k),pi2(i,j,k)
         !endif
      enddo

      4500 continue
   ENDDO
ENDDO

! If any level is above 100 mb, compute geostrophic winds

!GDKM=2.*SPCON
!DO J=2,NPRY-1
!   GDLAT=(XSWLAT+(J-1)*GDATDY)*PI180
!   FCORI=1./(2.*7.292E-5*SIN(GDLAT))
!   DO I=2,NPRX-1
!      DO K=NISN,1,-1
!         IF(PI2(I,J,K).LT.10000.) THEN
!            UI2(I,J,K)=-FCORI*(SI2(I,J+1,K)-SI2(I,J-1,K))/(GDKM*GDATDY)
!            VI2(I,J,K)= FCORI*(SI2(I+1,J,K)-SI2(I-1,J,K))  &
!                        /(GDKM*GDATDX*COS(GDLAT))
!         ENDIF
!      ENDDO
!   ENDDO
!ENDDO


return
END SUBROUTINE vterpp_s

!##############################################################################
Subroutine psfill (nz,vvv,aaa,np1,np2,i,j)

implicit none

integer :: nz,np1,np2,i,j,k
real :: vvv(nz),aaa(np1,np2,nz)

do k=1,nz
   aaa(i,j,k)=vvv(k)
enddo

return
END SUBROUTINE psfill
