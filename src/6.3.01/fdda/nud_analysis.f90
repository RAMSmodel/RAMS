!##############################################################################
Subroutine datassim ()

use mem_tend
use mem_basic
use mem_grid
use mem_varinit
use node_mod
use mem_leaf

implicit none

integer :: il,ir,jl,jr

!     Set bounds for nudging this sub-domain

il=ia  ;  ir=iz  ;  jl=ja  ;  jr=jz

!     West and east boundaries.
if(iand(ibcon,1) /= 0)  il=1
if(iand(ibcon,2) /= 0)  ir=mxp


!     South and north boundaries.
if(jdim == 1) then
   if(iand(ibcon,4) /= 0) jl=1
   if(iand(ibcon,8) /= 0) jr=myp
endif


! Basic boundary and analysis nudging scheme

CALL nudge (mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1) &
    ,varinit_g(ngrid)%varup(1,1,1) ,varinit_g(ngrid)%varvp(1,1,1)  &
    ,varinit_g(ngrid)%varpp(1,1,1) ,varinit_g(ngrid)%vartp(1,1,1)  &
    ,varinit_g(ngrid)%varrp(1,1,1)                                 & 
    ,varinit_g(ngrid)%varuf(1,1,1) ,varinit_g(ngrid)%varvf(1,1,1)  &
    ,varinit_g(ngrid)%varpf(1,1,1) ,varinit_g(ngrid)%vartf(1,1,1)  &
    ,varinit_g(ngrid)%varrf(1,1,1)                                 &
    ,basic_g(ngrid)%up(1,1,1)      ,basic_g(ngrid)%vp(1,1,1)       &
    ,basic_g(ngrid)%theta(1,1,1)   ,basic_g(ngrid)%rtp(1,1,1)      &
    ,basic_g(ngrid)%pp(1,1,1)                                      &
    ,basic_g(ngrid)%th00(1,1,1)    ,basic_g(ngrid)%rvt00(1,1,1)    &
    ,tend%ut(1),tend%vt(1),tend%tht(1),tend%rtt(1),tend%pt(1))


! Condensate nudging scheme

if (nud_cond == 1 .and. time >= tcond_beg .and. time <= tcond_end) &
 CALL nudge_cond (mzp,mxp,myp,il,ir,jl,jr,varinit_g(ngrid)%varwts(1,1,1) &
        ,varinit_g(ngrid)%varrph(1,1,1),varinit_g(ngrid)%varcph(1,1,1)   &
        ,varinit_g(ngrid)%varrfh(1,1,1),varinit_g(ngrid)%varcfh(1,1,1)   &
        ,basic_g(ngrid)%rtp(1,1,1),tend%rtt(1))

! Soil moisture nudging scheme

if (snudcent > 0.) &
 CALL nudge_soil (mxp,myp,il,ir,jl,jr &
        ,varinit_g(ngrid)%varsoilm1p(1,1)   &
        ,varinit_g(ngrid)%varsoilm2p(1,1)   &
        ,varinit_g(ngrid)%varsoilm1f(1,1)   &
        ,varinit_g(ngrid)%varsoilm2f(1,1)   &
        ,varinit_g(ngrid)%varsoilt1p(1,1)   &
        ,varinit_g(ngrid)%varsoilt2p(1,1)   &
        ,varinit_g(ngrid)%varsoilt1f(1,1)   &
        ,varinit_g(ngrid)%varsoilt2f(1,1)   &
        ,leaf_g(ngrid)%soil_text(1,1,1,1)   &
        ,leaf_g(ngrid)%soil_energy(1,1,1,1) &
        ,leaf_g(ngrid)%soil_water(1,1,1,1)  &
        ,leaf_g(ngrid)%leaf_class(1,1,1)    )

return
END SUBROUTINE datassim

!##############################################################################
Subroutine nudge (m1,m2,m3,ia,iz,ja,jz,varwts  &
     ,varup,varvp,varpp,vartp,varrp  &
     ,varuf,varvf,varpf,vartf,varrf  &
     ,up,vp,theta,rtp,pp,th00,rvt00,ut,vt,tht,rtt,pt)

use mem_grid
use mem_varinit
use micphys, only:idiffperts

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: varup,varvp,vartp,varrp,varpp  &
                            ,varuf,varvf,vartf,varrf,varpf  &
                            ,varwts,up,vp,theta,rtp,pp,th00,rvt00  &
                            ,ut,vt,tht,rtt,pt             

integer :: i,j,k
real :: tfact, wt_uv, wt_th, wt_pi, wt_rt
real :: var_tend, tot_weight

!         Linearly interpolate values in time, then nudge.

if (nud_type == 1) then
   tfact=(time-vtime1)/(vtime2-vtime1)
else
   stop 'should not be nudging if nud_type=0'
endif

wt_uv=wt_nudge_g(ngrid)* wt_nudge_uv
wt_th=wt_nudge_g(ngrid)* wt_nudge_th
wt_pi=wt_nudge_g(ngrid)* wt_nudge_pi
wt_rt=wt_nudge_g(ngrid)* wt_nudge_rt

do j=ja,jz
   do i=ia,iz
     do k=1,m1
         var_tend   = varup(k,i,j)+(varuf(k,i,j)-varup(k,i,j))*tfact
         tot_weight = (varwts(k,i,j)+varwts(k,min(m2,i+1),j))*.5* wt_uv
         ut(k,i,j)  = ut(k,i,j) + tot_weight*(var_tend-up(k,i,j))

         var_tend   = varvp(k,i,j)+(varvf(k,i,j)-varvp(k,i,j))*tfact
         tot_weight = (varwts(k,i,j)+varwts(k,i,min(m3,j+jdim)))*.5* wt_uv
         vt(k,i,j)  = vt(k,i,j) + tot_weight*(var_tend-vp(k,i,j))

         var_tend   = vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact
         tot_weight = varwts(k,i,j)* wt_th
         tht(k,i,j) = tht(k,i,j) + tot_weight*(var_tend-theta(k,i,j))

         var_tend   = varpp(k,i,j)+(varpf(k,i,j)-varpp(k,i,j))*tfact
         tot_weight = varwts(k,i,j)* wt_pi
         pt(k,i,j)  = pt(k,i,j) + tot_weight*(var_tend-pp(k,i,j))

         var_tend   = varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact
         tot_weight = varwts(k,i,j)* wt_rt
         rtt(k,i,j) = rtt(k,i,j) + tot_weight*(var_tend-rtp(k,i,j))
     enddo
   enddo
enddo

!Do this if diffusing THP,RTP from varfile base state perturbations
if(idiffperts==3) then
 do j=1,m3
   do i=1,m2
     do k=1,m1
         var_tend     = vartp(k,i,j)+(vartf(k,i,j)-vartp(k,i,j))*tfact
         th00(k,i,j)  = var_tend
         var_tend     = varrp(k,i,j)+(varrf(k,i,j)-varrp(k,i,j))*tfact
         rvt00(k,i,j) = var_tend
     enddo
   enddo
 enddo
endif

return
END SUBROUTINE nudge

!##############################################################################
Subroutine nudge_cond (m1,m2,m3,ia,iz,ja,jz,varwts  &
                     ,varrph,varcph,varrfh,varcfh,rtp,rtt)

use mem_grid
use mem_varinit

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: varrph,varcph,varrfh,varcfh  &
                            ,varwts,rtp,rtt             

integer :: i,j,k
real :: tfact, wt_rc
real :: varr_tend, varc_tend, nud_wt

tfact=(time-vtime1)/(vtime2-vtime1)

wt_rc=wt_nudgec(ngrid)/t_nudge_rc

do j=ja,jz
   do i=ia,iz

      do k=1,m1
         varr_tend = varrph(k,i,j)+(varrfh(k,i,j)-varrph(k,i,j))*tfact
         varc_tend = varcph(k,i,j)+(varcfh(k,i,j)-varcph(k,i,j))*tfact

         !Choose a nudging weight method below. By default we are only
         !applying the condensate nudging weight, but to the whole domain

         !This method sets nudging weights in conjunction with standard
         !analysis nudging weights in space and time         
         nud_wt = varwts(k,i,j)*wt_rc

         !This method applies only condensate nudging weights only for
         !all the domain. Currently use this to override option just above.
         nud_wt = wt_rc

        ! Only nudging total water RTP where condensate mixing ratio exists.
        ! Condensate will have a value of -9999 for traditional varfiles using
        ! "MAKEVFILE" since we do not ingest condensate from "dp" files. We do
        ! ingest condensate from History-Varfiles (MAKEHFILE), so values should
        ! be positive-definite.
         if (varc_tend > 0.)  &
            rtt(k,i,j)=rtt(k,i,j) + nud_wt*(varr_tend-rtp(k,i,j))
      enddo

   enddo
enddo

return
END SUBROUTINE nudge_cond

!##############################################################################
Subroutine nudge_soil (m2,m3,ia,iz,ja,jz                           &
                      ,varsoilm1p,varsoilm2p,varsoilm1f,varsoilm2f &
                      ,varsoilt1p,varsoilt2p,varsoilt1f,varsoilt2f &
                      ,soil_text,soil_energy,soil_water,leaf_class)

use mem_grid
use mem_varinit
use leaf_coms, only:slmsts,soilcp,tempk

implicit none

integer :: i,j,k,m2,m3,ip,ia,iz,ja,jz,nsoil
real, dimension(m2,m3) :: varsoilm1p,varsoilm2p,varsoilm1f,varsoilm2f &
                         ,varsoilt1p,varsoilt2p,varsoilt1f,varsoilt2f
real, dimension(nzg,m2,m3,npatch) :: soil_water,soil_text,soil_energy
real, dimension(m2,m3,npatch) :: leaf_class
real :: tfact,nud_wt,soilmtend,varsoilm2_tend,varsoilm1_tend &
       ,soilttend,varsoilt2_tend,varsoilt1_tend

tfact=(time-vtime1)/(vtime2-vtime1)

nud_wt=1.0/snudcent

do j=ja,jz
 do i=ia,iz
  do ip=2,npatch
   do k=1,nzg
    !Limit the soil moisture to the allowed maximum for soil textural class.
    !A similar limitation is placed on soil moisture at initialization time.
    nsoil = nint(soil_text(k,i,j,ip))
    varsoilm2p(i,j) = max(soilcp(nsoil),min(slmsts(nsoil),varsoilm2p(i,j)))
    varsoilm2f(i,j) = max(soilcp(nsoil),min(slmsts(nsoil),varsoilm2f(i,j)))
    varsoilm1p(i,j) = max(soilcp(nsoil),min(slmsts(nsoil),varsoilm1p(i,j)))
    varsoilm1f(i,j) = max(soilcp(nsoil),min(slmsts(nsoil),varsoilm1f(i,j)))
    !For persistent wetlands (bogs, marshes, fens, swamps) and irrigated
    !crops, retain saturated soil. Currently, this corresponds to
    !leaf classes 16, 17, and 20.
    if (nint(leaf_class(i,j,ip)) == 16 .or.  &
        nint(leaf_class(i,j,ip)) == 17 .or.  &
        nint(leaf_class(i,j,ip)) == 20) then
          varsoilm2p(i,j) = slmsts(nsoil)
          varsoilm2f(i,j) = slmsts(nsoil)
          varsoilm1p(i,j) = slmsts(nsoil)
          varsoilm1f(i,j) = slmsts(nsoil)
    endif
    !Compute soil moisture nudging tendencies
    if(k==nzg)then !top soil level (k==nzg)
     varsoilm2_tend = varsoilm2p(i,j)+(varsoilm2f(i,j)-varsoilm2p(i,j))*tfact
     soilmtend = nud_wt*(varsoilm2_tend-soil_water(k,i,j,ip))
    else
     varsoilm1_tend = varsoilm1p(i,j)+(varsoilm1f(i,j)-varsoilm1p(i,j))*tfact
     soilmtend = nud_wt*(varsoilm1_tend-soil_water(k,i,j,ip))
    endif
    !Apply nudging to soil moisture based on tendency and allowed maximum.
    !Update soil energy accordingly as well.
    soil_water(k,i,j,ip) = soil_water(k,i,j,ip) + (soilmtend * dtlt)
    soil_energy(k,i,j,ip) = soil_energy(k,i,j,ip) &
       + (soilmtend * dtlt) * 4.186e6 * (tempk(k) - 193.36)
   enddo
  enddo
 enddo
enddo

return
END SUBROUTINE nudge_soil

!##############################################################################
Subroutine varweight (n1,n2,n3,varwts,topt,rtgt)

use mem_grid
use mem_varinit
use node_mod

implicit none

integer :: n1,n2,n3
real :: varwts(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

integer :: i,j,k
real :: tnudcenti,tnudtopi,tnudlati,rown,rows,rowe,roww,zloc,wttop &
       ,wtlat,delzi

!         Get weights for large scale and model tendencies

if (nudlat .le. 0) return

tnudcenti=0.
if(tnudcent.gt. .01) tnudcenti=1./tnudcent
tnudtopi=0.
if(tnudtop.gt. .01) tnudtopi=1./tnudtop-tnudcenti
tnudlati=0.
if(tnudlat.gt. .01) tnudlati=1./tnudlat-tnudcenti

delzi=0.0 !Variable initialized
if(ztop.gt.znudtop) then
   delzi=1./(ztop-znudtop)
elseif(tnudtop.gt. .01) then
   print*,'Incorrect specification of znudtop ! ! !'
   print*,' znudtop = ',znudtop
   print*,'    ztop = ',ztop
   stop 'varwt-znud'
endif


do j=1,n3
   do i=1,n2

!                       quadratic weight func for lateral boundaries

      ! Calculate row[nsew] considering the full domain. n1,n2,n3 will be
      ! set to mzp,mxp,myp (sub-domain sizes). In the case of just one node,
      ! m[xyz]p will be equal to n[xyz]p.
      if (nyp > 1) then
         rown=max(0.,float((j+j0)+nudlat-nyp))
         rows=max(0.,float(nudlat+1-(j+j0)))
      else
         rown=0.
         rows=0.
      endif

      rowe=max(0.,float((i+i0)+nudlat-nxp))
      roww=max(0.,float(nudlat+1-(i+i0)))
      wtlat=max(rown*rown,rows*rows,rowe*rowe,roww*roww)  &
           /float(nudlat*nudlat)

!                       linear weight func for top boundary

      do k=1,n1
         zloc=ztn(k,1)*rtgt(i,j)+topt(i,j)
         wttop=max(0.,(zloc-znudtop)*delzi)

!                       full 3-D weight func

         varwts(k,i,j)=tnudcenti  &
              +max(tnudlati*wtlat,tnudtopi*wttop)

      enddo

   enddo
enddo

return
END SUBROUTINE varweight
