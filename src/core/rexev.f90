
!##############################################################################
Subroutine exevolve (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,key,thvlast)

use mem_tend
use mem_basic
use mem_grid


implicit none

character*(*) :: key
integer :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,i,j,k
real,dimension(mzp,mxp,myp) :: thvlast

if(key.eq.'ADV')then

  do i=ia,iz
    do j=ja,jz
      do k=1,mzp
        thvlast(k,i,j)=0.0
      enddo
    enddo
  enddo

  CALL exadvlf (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv                  &
    ,grid_g(ngrid)%rtgu   (1,1)   ,grid_g(ngrid)%fmapui (1,1)    &
    ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%fmapvi (1,1)    &
    ,grid_g(ngrid)%f13t   (1,1)   ,grid_g(ngrid)%f23t   (1,1)    &
    ,grid_g(ngrid)%rtgt   (1,1)   ,grid_g(ngrid)%fmapt  (1,1)    &
    ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%dyt    (1,1)    &
    ,basic_g(ngrid)%uc    (1,1,1) ,basic_g(ngrid)%dn0u  (1,1,1)  &
    ,basic_g(ngrid)%vc    (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
    ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%wc    (1,1,1)  &
    ,tend%pt              (1)     ,basic_g(ngrid)%pc    (1,1,1)  )

  CALL excondiv (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv                 &
    ,basic_g(ngrid)%uc    (1,1,1) ,basic_g(ngrid)%vc    (1,1,1)  &
    ,basic_g(ngrid)%wc    (1,1,1) ,basic_g(ngrid)%pc    (1,1,1)  &
    ,tend%pt              (1)                                    &
    ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%dyt    (1,1)    &
    ,grid_g(ngrid)%rtgt   (1,1)   ,grid_g(ngrid)%rtgu   (1,1)    &
    ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%f13t   (1,1)    &
    ,grid_g(ngrid)%f23t   (1,1)   ,grid_g(ngrid)%fmapt  (1,1)    &
    ,grid_g(ngrid)%fmapui (1,1)   ,grid_g(ngrid)%fmapvi (1,1)    )

  CALL fthvlast (mzp,mxp,myp,ia,iz,ja,jz,thvlast                 &
    ,basic_g(ngrid)%theta(1,1,1)  ,basic_g(ngrid)%rtp(1,1,1)     &
    ,basic_g(ngrid)%rv(1,1,1))

elseif(key.eq.'THV')then

  CALL exheat (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,thvlast)

endif

return
END SUBROUTINE exevolve

!##############################################################################
Subroutine exadvlf (m1,m2,m3,ia,iz,ja,jz,izu,jzv,rtgu,fmapui,rtgv,fmapvi &
         ,f13t,f23t,rtgt,fmapt,dxt,dyt,uc,dn0u,vc,dn0v,dn0,wc,pt,pc)

use mem_basic
use mem_grid

implicit none

integer i,j,k,m1,m2,m3,jm,im,ja,jz,ia,iz,jzv,izu
real :: c1z,c1x,c1y
real, dimension(m2,m3) :: rtgu,fmapui,rtgv,fmapvi,f13t,f23t,rtgt,fmapt,dxt,dyt
real, dimension(m1,m2,m3) :: flxu,flxv,uc,dn0u,vc,dn0v,dn0,wc,flxw,pt,pc

! Compute momentum fluxes flxu, flxv, flxw

do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         flxu(k,i,j) = uc(k,i,j) * dn0u(k,i,j) * rtgu(i,j)  &
            * fmapui(i,j)
         flxv(k,i,j) = vc(k,i,j) * dn0v(k,i,j) * rtgv(i,j)  &
            * fmapvi(i,j)
      enddo
   enddo
enddo

if(itopo.eq.0) then
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1-1
            flxw(k,i,j) = wc(k,i,j)  &
               * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
         enddo
      enddo
   enddo
else
   do j = 1,m3
      jm = max(j-1,1)
      do i = 1,m2
         im = max(i-1,1)
         do k = 1,m1-1
            flxw(k,i,j) = wc(k,i,j)  &
               * .5 * (dn0(k,i,j) + dn0(k+1,i,j))  &
               + hw4(k) * ((flxu(k,i,j) + flxu(k+1,i,j)  &
               + flxu(k,im,j) + flxu(k+1,im,j)) * f13t(i,j)  &
               + (flxv(k,i,j) + flxv(k+1,i,j)  &
               + flxv(k,i,jm) + flxv(k+1,i,jm)) * f23t(i,j))
         enddo
      enddo
   enddo
endif

! Compute advection contribution to PI tendency

do j = ja,jz
   do i = ia,izu
      c1z = 0.5 / rtgt(i,j)
      c1x = c1z * fmapt(i,j) * dxt(i,j)
      c1y = c1z * fmapt(i,j) * dyt(i,j)

      do k = 2,m1-1
         pt(k,i,j) = pt(k,i,j) - c1x / dn0(k,i,j) * (  &
              flxu(k,i,j)  &
               * (pc(k,i,j) + pc(k,i+1,j))  &
            - flxu(k,i-1,j)  &
               * (pc(k,i,j) + pc(k,i-1,j))  &
            - (flxu(k,i,j) - flxu(k,i-1,j)) * 2.* pc(k,i,j) )
      enddo
   enddo
enddo

do j=ja,jzv
  do i=ia,iz
      do k=2,m1-1
        pt(k,i,j)=pt(k,i,j) - c1y /dn0(k,i,j) * ( &
          flxv(k,i,j)  &
           * (pc(k,i,j)+pc(k,i,j+jdim))  &
         -flxv(k,i,j-jdim)  &
           * (pc(k,i,j)+pc(k,i,j-jdim))  &
           -  (flxv(k,i,j)-flxv(k,i,j-jdim))*2.*pc(k,i,j) )
      enddo
   enddo
enddo

do j=ja,jz
  do i=ia,iz
      do k=2,m1-1
        pt(k,i,j)=pt(k,i,j) - c1z * dzt(k) /dn0(k,i,j) * ( &
          flxw(k,i,j)  &
           * (pc(k,i,j)+pc(k+1,i,j))  &
         -flxw(k-1,i,j)  &
           * (pc(k,i,j)+pc(k-1,i,j))  &
           -  (flxw(k,i,j)-flxw(k-1,i,j))*2.*pc(k,i,j) )
      enddo
   enddo
enddo

return
END SUBROUTINE exadvlf

!##############################################################################
Subroutine excondiv (m1,m2,m3,ia,iz,ja,jz,izu,jzv,uc,vc,wc,pc,pt  &
    ,dxt,dyt,rtgt,rtgu,rtgv,f13t,f23t,fmapt,fmapui,fmapvi)

use mem_basic
use rconstants
use mem_grid

implicit none

integer :: m1,m2,m3,i,j,k,ia,iz,ja,jz,izu,jzv,im,jm
real :: c1z,c1x,c1y
real, dimension(m1,m2,m3) :: flxu,flxv,flxw,uc,vc,wc,pt,pc
real, dimension(m2,m3) :: rtgu,fmapui,rtgv,fmapvi,rtgt,fmapt,dxt,dyt,f13t,f23t
                             
! Compute divergence
!-----------
! Prep Fluxes
!-----------

!  These are:  (transformed velocities) times (a) times (mapfactor)
do j=1,m3
  do i=1,m2
    do k=1,m1
      flxu(k,i,j)=uc(k,i,j)*rtgu(i,j)*fmapui(i,j)
      flxv(k,i,j)=vc(k,i,j)*rtgv(i,j)*fmapvi(i,j)
    enddo
  enddo
enddo

if(itopo.eq.0)then
  do j=1,m3
    do i=1,m2
      do k=1,m1-1
        flxw(k,i,j)=wc(k,i,j)
      enddo
    enddo
  enddo
else
  do j=1,m3
    jm=max(j-1,1)
    do i=1,m2
      im=max(i-1,1)
      do k=1,m1-1
        flxw(k,i,j)=wc(k,i,j) &
           + hw4(k) * ( (flxu(k,i,j)+flxu(k+1,i,j) &
           +flxu(k,im,j) + flxu(k+1,im,j)) * f13t(i,j) &
           + (flxv(k,i,j)+flxv(k+1,i,j) &
           +flxv(k,i,jm) + flxu(k+1,i,jm)) * f23t(i,j)) 
      enddo
    enddo
  enddo
endif

do j=ja,jz
  do i=ia,izu
    c1x=fmapt(i,j)*dxt(i,j)/rtgt(i,j)
    do k=2,m1-1
      pt(k,i,j)=pt(k,i,j)     & 
        - c1x * ( flxu(k,i,j)-flxu(k,i-1,j) )   &
        * pc(k,i,j) * rocv
    enddo
  enddo
enddo

do j=ja,jzv
  do i=ia,iz
    c1y=fmapt(i,j)*dyt(i,j)/rtgt(i,j)
    do k=2,m1-1
      pt(k,i,j)=pt(k,i,j)  &
        - c1y * (flxv(k,i,j)-flxv(k,i,j-jdim) )  &
        * pc(k,i,j) *rocv
    enddo
  enddo
enddo

do j=ja,jz
  do i=ia,iz
    c1z=1.0/rtgt(i,j)
    do k=2,m1-1
      pt(k,i,j)=pt(k,i,j)  &
        - c1z * dzm(k) * (flxw(k,i,j)-flxw(k-1,i,j) )   &
        * pc(k,i,j) * rocv
    enddo
  enddo
enddo

return
END SUBROUTINE excondiv

!##############################################################################
Subroutine fthvlast (m1,m2,m3,ia,iz,ja,jz,thvlast,theta,rtp,rv)

use mem_basic
use mem_grid
use micphys, only:level

implicit none

integer :: m1,m2,m3,i,j,k,ia,iz,ja,jz
real,dimension(m1,m2,m3) :: thvlast,theta,rtp,rv

if(level.gt.0)then
  do j=ja,jz
    do i=ia,iz
      do k=1,m1
        thvlast(k,i,j)=theta(k,i,j)*(1.0+1.61*rv(k,i,j))/(1.0+rtp(k,i,j))
      enddo
    enddo
  enddo
else
  do j=ja,jz
    do i=ia,iz
      do k=1,m1
        thvlast(k,i,j)=theta(k,i,j)
      enddo
    enddo
  enddo
endif

return
END SUBROUTINE fthvlast

!##############################################################################
Subroutine exheat (m1,m2,m3,ia,iz,ja,jz,izu,jzv,thvlast)

use mem_tend
use mem_basic
use mem_grid, only:dtlt,ngrid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,mxyzp,i,j,k
real, dimension(m1,m2,m3) :: thvtend,thvadv,thvlast

CALL thetvtend (m1,m2,m3,ia,iz,ja,jz,izu,jzv,basic_g(ngrid)%theta(1,1,1) &
    ,basic_g(ngrid)%rtp(1,1,1),basic_g(ngrid)%rv(1,1,1),thvtend,thvlast)

CALL exthvadv (m1,m2,m3,ia,iz,ja,jz,izu,jzv,thvadv)

CALL exhtend (m1,m2,m3,ia,iz,ja,jz,basic_g(ngrid)%pi0(1,1,1) &
    ,basic_g(ngrid)%pc    (1,1,1) ,basic_g(ngrid)%rtp(1,1,1) &
    ,basic_g(ngrid)%theta (1,1,1) ,basic_g(ngrid)%rv (1,1,1) &
    ,tend%pt,thvtend,thvadv)

return
END SUBROUTINE exheat

!##############################################################################
Subroutine thetvtend (m1,m2,m3,ia,iz,ja,jz,izu,jzv,theta,rtp,rv &
                     ,thvtend,thvlast)

use mem_basic
use mem_grid, only:dtlt,time
use micphys, only:level

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,mxyzp,i,j,k
real :: dtlti
real, dimension(m1,m2,m3) :: theta,rtp,rv,thvtend,thvlast

if(time.le.0.1)then
    do j=ja,jz
      do i=ia,iz
        do k=1,m1
          thvtend(k,i,j)=0.0
        enddo
      enddo
    enddo
  return
endif

dtlti=1.0/dtlt

if(level.gt.0)then
  do j=ja,jz
    do i=ia,iz
      do k=2,m1-1
        thvtend(k,i,j) = (  theta(k,i,j) * (1.0+1.61*rv(k,i,j))  &
          / (1.0 + rtp(k,i,j))   -  thvlast(k,i,j)  )  * dtlti
      enddo
    enddo
  enddo
else
  do j=ja,jz
    do i=ia,iz
      do k=1,m1
        thvtend(k,i,j) = ( theta(k,i,j) - thvlast(k,i,j) ) * dtlti
      enddo
    enddo
  enddo
endif  

return
END SUBROUTINE thetvtend

!##############################################################################
Subroutine exthvadv (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,thvadv)

use mem_basic
use mem_grid
use mem_scratch

implicit none

integer :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,i,j,k,ind,imono
real :: dtlto2
real, dimension(mzp,mxp,myp) :: thvadv,thetav,srthtv

dtlto2 = .5 * dtlt
ind = 0
do j = 1,myp
   do i = 1,mxp
      do k = 1,mzp
         ind = ind + 1
         scratch%vt3da(ind) = (basic_g(ngrid)%uc(k,i,j)  &
            + basic_g(ngrid)%uc(k,i,j)) * dtlto2
         scratch%vt3db(ind) = (basic_g(ngrid)%vc(k,i,j)  &
            + basic_g(ngrid)%vc(k,i,j)) * dtlto2
         scratch%vt3dc(ind) = (basic_g(ngrid)%wc(k,i,j)  &
            + basic_g(ngrid)%wc(k,i,j)) * dtlto2
      enddo
   enddo
enddo

CALL prep_thetv (mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%theta(1,1,1) &
         ,basic_g(ngrid)%rtp(1,1,1),basic_g(ngrid)%rv(1,1,1),thetav)

CALL fa_preptc (mzp,mxp,myp        &
      ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
      ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
      ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
      ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
      ,scratch%vt3dj        (1)     ,scratch%vt3dk        (1)      &
      ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%dn0u  (1,1,1)  &
      ,basic_g(ngrid)%dn0v  (1,1,1) ,grid_g(ngrid)%rtgt   (1,1)    &
      ,grid_g(ngrid)%rtgu   (1,1)   ,grid_g(ngrid)%rtgv   (1,1)    &
      ,grid_g(ngrid)%fmapt  (1,1)   ,grid_g(ngrid)%fmapui (1,1)    &
      ,grid_g(ngrid)%fmapvi (1,1)   ,grid_g(ngrid)%f13t   (1,1)    &
      ,grid_g(ngrid)%f23t   (1,1)   ,grid_g(ngrid)%dxu    (1,1)    &
      ,grid_g(ngrid)%dyv    (1,1)   ,grid_g(ngrid)%dxt    (1,1)    &
      ,grid_g(ngrid)%dyt    (1,1))

do i=1,mxp
  do j=1,myp
    do k=1,mzp
      srthtv(k,i,j)=thetav(k,i,j)
    enddo
  enddo
enddo

imono = 1 !flag for monotonic flux limiter in Z only

CALL fa_xc (mzp,mxp,myp,ia,iz                     &
           ,thetav,srthtv                         &
           ,scratch%vt3da (1) ,scratch%vt3dd (1)  &
           ,scratch%vt3dg (1) ,scratch%vt3dh (1)  &
           ,scratch%vt3di (1))


if (jdim .eq. 1) &
 CALL fa_yc (mzp,mxp,myp,ia,iz,ja,jz              &
           ,thetav,srthtv                         &
           ,scratch%vt3db (1) ,scratch%vt3de (1)  &
           ,scratch%vt3dg (1) ,scratch%vt3dj (1)  &
           ,scratch%vt3di (1) ,jdim)


if (imono.eq.1) then
 CALL fluxlimits (mzp,mxp,myp,ia,iz,ja,jz         &
           ,srthtv,thetav                         &
           ,scratch%vt3dc(1), scratch%vt3df(1)    &
           ,scratch%vt3dg(1), scratch%vt3dk(1)    &
           ,vctr1,vctr2,basic_g(ngrid)%dn0)
else         
 CALL fa_zc (mzp,mxp,myp,ia,iz,ja,jz              &
           ,thetav,srthtv                         &
           ,scratch%vt3dc (1) ,scratch%vt3df (1)  &
           ,scratch%vt3dg (1) ,scratch%vt3dk (1)  &
           ,vctr1,vctr2)
endif

do j=ja,jz
  do i=ia,iz
    do k=2,mzp-1
      thvadv(k,i,j)=0.0
    enddo
  enddo
enddo

CALL advtndc (mzp,mxp,myp,ia,iz,ja,jz,thetav,srthtv,thvadv,dtlt)

do j=ja,jz
  do i=ia,iz
    do k=2,mzp-1
      thvadv(k,i,j)=-1.0*thvadv(k,i,j)
    enddo
  enddo
enddo

return
END SUBROUTINE exthvadv

!##############################################################################
Subroutine prep_thetv (mzp,mxp,myp,ia,iz,ja,jz,theta,rtp,rv,thetav)

use micphys, only:level

implicit none

integer :: mzp,mxp,myp,k,i,j,ia,iz,ja,jz
real, dimension(mzp,mxp,myp) :: theta,rtp,rv,thetav

if(level.gt.0)then
  do i=1,mxp
    do j=1,myp
      do k=1,mzp
        thetav(k,i,j)=theta(k,i,j)*(1.0+1.61*rv(k,i,j))/(1.0+rtp(k,i,j))
      enddo
    enddo
  enddo
else
  do i=1,mxp
    do j=1,myp
      do k=1,mzp
        thetav(k,i,j)=theta(k,i,j)
      enddo
    enddo
  enddo
endif

return
END SUBROUTINE prep_thetv

!##############################################################################
Subroutine exhtend (mzp,mxp,myp,ia,iz,ja,jz,pi0,pc,rtp,theta,rv,pt &
                   ,thvtend,thvadv)

use mem_tend
use mem_basic
use rconstants
use micphys, only:level

implicit none

integer mzp,mxp,myp,ia,iz,ja,jz,i,j,k
real, dimension(mzp,mxp,myp) :: pi0,pc,rtp,theta,pt,thvtend,thvadv,rv

if(level.gt.0)then
  do j=ja,jz
    do i=ia,iz
      do k=2,mzp-1
        pt(k,i,j) = pt(k,i,j) + rocv * (pi0(k,i,j) + pc(k,i,j))  &
          / ( theta(k,i,j) * (1.0+1.61*rv(k,i,j)) / (1.0+rtp(k,i,j)) )   &
          * (  thvtend(k,i,j) + thvadv(k,i,j)  )
      enddo
    enddo
  enddo
else
  do j=ja,jz
    do i=ia,iz
      do k=2,mzp-1
        pt(k,i,j) = pt(k,i,j) + rocv * (pi0(k,i,j) + pc(k,i,j))  &
          / theta(k,i,j) * ( thvtend(k,i,j) + thvadv(k,i,j) )
      enddo
    enddo
  enddo
endif

return
END SUBROUTINE exhtend

