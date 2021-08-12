!##############################################################################
Subroutine advectc (varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv)

use mem_tend
use var_tables
use mem_scratch
use mem_grid
use mem_basic

implicit none

integer :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,n,mxyzp
integer :: imono
character(len=*) :: varn
real :: dtlto2
real, dimension(maxgrds), save :: save_dtlt
integer :: i,j,k,ind
real, pointer :: scalarp, scalart
integer, dimension(maxgrds), save :: ncall
data ncall/maxgrds*0/

if (ncall(ngrid) == 0 .or. dtlt .ne. save_dtlt(ngrid) ) then
  ncall(ngrid) = 1
  save_dtlt(ngrid) = dtlt
endif

mxyzp = mxp * myp * mzp
 !  basic_g(ngrid)%up    (1:nzp,1:nxp,1:nyp) = 10.
 !  basic_g(ngrid)%uc    (1:nzp,1:nxp,1:nyp) = 10.
 !  basic_g(ngrid)%vp    (1:nzp,1:nxp,1:nyp) = 0.
 !  basic_g(ngrid)%vc    (1:nzp,1:nxp,1:nyp) = 0.
 !  basic_g(ngrid)%wp    (1:nzp,1:nxp,1:nyp) = 0.
 !  basic_g(ngrid)%wc    (1:nzp,1:nxp,1:nyp) = 0.

if (varn .eq. 'V' .or. varn .eq. 'ALL') then

 ! Advect  U, V, and W
 CALL vel_advectc (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv  &
         ,basic_g(ngrid)%uc    (1,1,1) ,basic_g(ngrid)%vc    (1,1,1)  &
         ,basic_g(ngrid)%wc    (1,1,1) ,tend%ut              (1)      &
         ,tend%vt              (1)     ,tend%wt              (1)      &
         ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%dn0u  (1,1,1)  &
         ,basic_g(ngrid)%dn0v  (1,1,1) ,scratch%vt3da        (1)      &
         ,scratch%vt3db        (1)     ,scratch%vt3dc        (1)      &
         ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%dxu    (1,1)    &
         ,grid_g(ngrid)%dxv    (1,1)   ,grid_g(ngrid)%dyt    (1,1)    &
         ,grid_g(ngrid)%dyu    (1,1)   ,grid_g(ngrid)%dyv    (1,1)    &
         ,grid_g(ngrid)%rtgt   (1,1)   ,grid_g(ngrid)%rtgu   (1,1)    &
         ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%f13t   (1,1)    &
         ,grid_g(ngrid)%f23t   (1,1)   ,grid_g(ngrid)%fmapt  (1,1)    &
         ,grid_g(ngrid)%fmapu  (1,1)   ,grid_g(ngrid)%fmapv  (1,1)    &
         ,grid_g(ngrid)%fmapui (1,1)   ,grid_g(ngrid)%fmapvi (1,1)    )

endif

if (varn .eq. 'T' .or. varn .eq. 'ALL') THEN

   ! Advect  scalars

   dtlto2 = .5 * dtlt
   ind = 0
   do j = 1,myp
      do i = 1,mxp
         do k = 1,mzp
            ind = ind + 1
            scratch%vt3da(ind) = (basic_g(ngrid)%up(k,i,j)  &
               + basic_g(ngrid)%uc(k,i,j)) * dtlto2
            scratch%vt3db(ind) = (basic_g(ngrid)%vp(k,i,j)  &
               + basic_g(ngrid)%vc(k,i,j)) * dtlto2
            scratch%vt3dc(ind) = (basic_g(ngrid)%wp(k,i,j)  &
               + basic_g(ngrid)%wc(k,i,j)) * dtlto2
         enddo
      enddo
   enddo

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

   !flag for monotonic flux limiter in Z only.
   !imono=1 option should be considered experimental.
   !May need bug fixing (Saleeby Oct 5, 2020).
   !pp and pm terms might not have density multiplied correctly.
   imono = 0

   do n=1,num_scalar(ngrid)
      
      scalarp => scalar_tab(n,ngrid)%var_p
      scalart => scalar_tab(n,ngrid)%var_t
    
      CALL atob (mxyzp,scalarp,scratch%scr1)

      CALL fa_xc (mzp,mxp,myp,ia,iz                &
            ,scalarp,scratch%scr1  (1)             &
            ,scratch%vt3da (1) ,scratch%vt3dd (1)  &
            ,scratch%vt3dg (1) ,scratch%vt3dh (1)  &
            ,scratch%vt3di (1))

      if (jdim .eq. 1)  &
      CALL fa_yc (mzp,mxp,myp,ia,iz,ja,jz           &
            ,scalarp,scratch%scr1  (1)             &
            ,scratch%vt3db (1) ,scratch%vt3de (1)  &
            ,scratch%vt3dg (1) ,scratch%vt3dj (1)  &
            ,scratch%vt3di (1) ,jdim)

      if (imono.eq.1) then
        CALL fluxlimits (mzp,mxp,myp,ia,iz,ja,jz       &
               ,scratch%scr1(1),scalarp               &
               ,scratch%vt3dc(1), scratch%vt3df(1)    &
               ,scratch%vt3dg(1), scratch%vt3dk(1)    &
               ,vctr1,vctr2,basic_g(ngrid)%dn0)
      else
        CALL fa_zc (mzp,mxp,myp,ia,iz,ja,jz            &
               ,scalarp,scratch%scr1  (1)             &
               ,scratch%vt3dc (1) ,scratch%vt3df (1)  &
               ,scratch%vt3dg (1) ,scratch%vt3dk (1)  &
               ,vctr1,vctr2)
      endif

      CALL advtndc (mzp,mxp,myp,ia,iz,ja,jz    &
            ,scalarp,scratch%scr1 (1)  &
            ,scalart,dtlt)

   enddo
 
endif

return
END SUBROUTINE advectc

!##############################################################################
Subroutine vel_advectc (m1,m2,m3,ia,iz,ja,jz,izu,jzv  &
   ,uc,vc,wc,ut,vt,wt,dn0,dn0u,dn0v,flxu,flxv,flxw  &
   ,dxt,dxu,dxv,dyt,dyu,dyv,rtgt,rtgu,rtgv,f13t,f23t  &
   ,fmapt,fmapu,fmapv,fmapui,fmapvi)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,izu,jzv
real, dimension(m1,m2,m3) :: uc,vc,wc,ut,vt,wt,dn0,dn0u,dn0v,flxu,flxv,flxw
real, dimension(m2,m3) :: dxt,dxu,dxv,dyt,dyu,dyv,rtgt,rtgu,rtgv,f13t,f23t  &
   ,fmapt,fmapu,fmapv,fmapui,fmapvi
integer :: j,i,k,jm,im
real :: c1z,c1x,c1y

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

! Compute advection contribution to U tendency
do j = ja,jz
   do i = ia,izu
      c1z = .25 / rtgu(i,j)
      c1x = c1z * fmapu(i,j) * dxu(i,j)
      c1y = c1z * fmapu(i,j) * dyu(i,j)

      do k = 2,m1-1
         ut(k,i,j) = ut(k,i,j) + c1x / dn0u(k,i,j) * (  &
              (flxu(k,i,j) + flxu(k,i-1,j))  &
               * (uc(k,i,j) + uc(k,i-1,j))  &
            - (flxu(k,i,j) + flxu(k,i+1,j))  &
               * (uc(k,i,j) + uc(k,i+1,j))  &
            + (flxu(k,i+1,j) - flxu(k,i-1,j)) * 2.* uc(k,i,j) )
      enddo

      do k = 2,m1-1
         ut(k,i,j) = ut(k,i,j) + c1y / dn0u(k,i,j) * (  &
              (flxv(k,i,j-jdim) + flxv(k,i+1,j-jdim))  &
               * (uc(k,i,j) + uc(k,i,j-jdim))  &
            - (flxv(k,i,j) + flxv(k,i+1,j))  &
               * (uc(k,i,j) + uc(k,i,j+jdim))&
            + (flxv(k,i,j) + flxv(k,i+1,j) - flxv(k,i,j-jdim)  &
            - flxv(k,i+1,j-jdim)) * 2.* uc(k,i,j) )
      enddo

      do k = 2,m1-1
         ut(k,i,j) = ut(k,i,j) + c1z * dzt(k) / dn0u(k,i,j) * (  &
              (flxw(k-1,i,j) + flxw(k-1,i+1,j))  &
               * (uc(k,i,j) + uc(k-1,i,j))  &
            - (flxw(k,i,j) + flxw(k,i+1,j))  &
               * (uc(k,i,j) + uc(k+1,i,j))   &
            + (flxw(k,i,j) + flxw(k,i+1,j) - flxw(k-1,i,j)  &
            - flxw(k-1,i+1,j)) * 2.* uc(k,i,j) )
      enddo
   enddo
enddo

! Compute advection contribution to V tendency

do j = ja,jzv
   do i = ia,iz
      c1z = .25 / rtgv(i,j)
      c1x = c1z * fmapv(i,j) * dxv(i,j)
      c1y = c1z * fmapv(i,j) * dyv(i,j)

      do k = 2,m1-1
         vt(k,i,j) = vt(k,i,j) + c1x / dn0v(k,i,j) * (  &
              (flxu(k,i-1,j) + flxu(k,i-1,j+jdim))  &
               * (vc(k,i,j) + vc(k,i-1,j))  &
            - (flxu(k,i,j) + flxu(k,i,j+jdim))  &
               * (vc(k,i,j) + vc(k,i+1,j))  &
            + (flxu(k,i,j) + flxu(k,i,j+jdim) - flxu(k,i-1,j)  &
            - flxu(k,i-1,j+jdim)) * 2.* vc(k,i,j) )
      enddo

      do k = 2,m1-1
         vt(k,i,j) = vt(k,i,j) + c1y / dn0v(k,i,j) * (  &
              (flxv(k,i,j) + flxv(k,i,j-jdim))  &
               * (vc(k,i,j) + vc(k,i,j-jdim))  &
            - (flxv(k,i,j) + flxv(k,i,j+jdim))  &
               * (vc(k,i,j) + vc(k,i,j+jdim))  &
            + (flxv(k,i,j+jdim) - flxv(k,i,j-jdim))  &
            * 2.* vc(k,i,j) )
      enddo

      do k = 2,m1-1
         vt(k,i,j) = vt(k,i,j) + c1z * dzt(k) / dn0v(k,i,j) * (  &
              (flxw(k-1,i,j) + flxw(k-1,i,j+jdim))  &
               * (vc(k,i,j) + vc(k-1,i,j))  &
            - (flxw(k,i,j) + flxw(k,i,j+jdim))  &
               * (vc(k,i,j) + vc(k+1,i,j))  &
            + (flxw(k,i,j) + flxw(k,i,j+jdim) - flxw(k-1,i,j)  &
            - flxw(k-1,i,j+jdim)) * 2.* vc(k,i,j) )
      enddo
   enddo
enddo
! Compute advection contribution to W tendency

do j = ja,jz
   do i = ia,iz
      c1z = .5 / rtgt(i,j)
      c1x = c1z * fmapt(i,j) * dxt(i,j)
      c1y = c1z * fmapt(i,j) * dyt(i,j)

      do k = 2,m1-2
         wt(k,i,j) = wt(k,i,j)  &
            + c1x / (dn0(k,i,j) + dn0(k+1,i,j)) * (  &
              (flxu(k,i-1,j) + flxu(k+1,i-1,j))  &
               * (wc(k,i,j) + wc(k,i-1,j))  &
            - (flxu(k,i,j) + flxu(k+1,i,j))  &
               * (wc(k,i,j) + wc(k,i+1,j))  &
            + (flxu(k,i,j) + flxu(k+1,i,j) - flxu(k,i-1,j)  &
            - flxu(k+1,i-1,j)) * 2.* wc(k,i,j) )
      enddo

      do k = 2,m1-2
         wt(k,i,j) = wt(k,i,j)  &
            + c1y / (dn0(k,i,j) + dn0(k+1,i,j)) * (  &
              (flxv(k,i,j-jdim) + flxv(k+1,i,j-jdim))  &
               * (wc(k,i,j) + wc(k,i,j-jdim))  &
            - (flxv(k,i,j) + flxv(k+1,i,j))  &
               * (wc(k,i,j) + wc(k,i,j+jdim))  &
            + (flxv(k,i,j) + flxv(k+1,i,j) - flxv(k,i,j-jdim)  &
            - flxv(k+1,i,j-jdim)) * 2.* wc(k,i,j) )
      enddo

      do k = 2,m1-2
         wt(k,i,j) = wt(k,i,j)  &
            + c1z * dzm(k) / (dn0(k,i,j) + dn0(k+1,i,j)) * (  &
              (flxw(k,i,j) + flxw(k-1,i,j))  &
               * (wc(k,i,j) + wc(k-1,i,j))  &
            - (flxw(k,i,j) + flxw(k+1,i,j))  &
               * (wc(k,i,j) + wc(k+1,i,j))   &
            + (flxw(k+1,i,j) - flxw(k-1,i,j)) * 2.* wc(k,i,j) )
      enddo
   enddo
enddo

return
END SUBROUTINE vel_advectc

!##############################################################################
Subroutine fa_preptc (m1,m2,m3,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
   ,vt3dh,vt3di,vt3dj,vt3dk,dn0,dn0u,dn0v  &
   ,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
   ,dxu,dyv,dxt,dyt)

use mem_grid
use mem_scratch

implicit none

integer :: m1,m2,m3,j,i,k,im,ip,jm,jp
real :: c1,c2,c3,c4,rtgti
real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
   ,vt3dh,vt3di,vt3dj,vt3dk,dn0,dn0u,dn0v
real, dimension(m2,m3) :: rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
   ,dxu,dyv,dxt,dyt

! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
! between past and current time levels) times dtlt.

! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
!    and include 1/rtgt factor in VT3DC
! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
! Compute weight at scalar point: VT3DH
! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2
do j = 1,m3
   jm = max(1,j-1)
   jp = min(m3,j+1)
   do i = 1,m2
      im = max(1,i-1)
      ip = min(m2,i+1)
      rtgti = 1. / rtgt(i,j)
      c1 = .5 * dxu(i,j)
      c2 = .5 * dyv(i,j)
      c3 = dxt(i,j) * fmapt(i,j) * rtgti
      c4 = dyt(i,j) * fmapt(i,j) * rtgti

      do k = 1,m1-1
         vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
            + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
            + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
            + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
            + vt3dc(k,i,j) * rtgti
         vt3dd(k,i,j) = c1 * vt3da(k,i,j)
         vt3de(k,i,j) = c2 * vt3db(k,i,j)
         vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
         vctr3(k) = 1. / dn0(k,i,j)
         vt3dh(k,i,j) = c3 * vctr3(k)
         vt3dj(k,i,j) = c4 * vctr3(k)
         vt3dk(k,i,j) = dzt(k) * vctr3(k)
      enddo

!            vt3di(1,i,j) = dxu(i,j) / (dxu(i,j) + dxt(ip,j))
!            vt3di(2,i,j) = dxu(i,j) / (dxu(i,j) + dxt(i,j))
!            vt3di(3,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,jp))
!            vt3di(4,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,j))
! temporary override
      vt3di(1,i,j) = .5
      vt3di(2,i,j) = .5
      vt3di(3,i,j) = .5
      vt3di(4,i,j) = .5
   enddo
enddo

do k = 1,m1-1
   vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
   vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
enddo

! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
! into mass fluxes times dtlt.

do j = 1,m3
   do i = 1,m2
      c1 = fmapui(i,j) * rtgu(i,j)
      c2 = fmapvi(i,j) * rtgv(i,j)
      do k = 1,m1-1
         vt3da(k,i,j) = vt3da(k,i,j) * c1 * dn0u(k,i,j)
         vt3db(k,i,j) = vt3db(k,i,j) * c2 * dn0v(k,i,j)
         vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
            * (dn0(k,i,j) + dn0(k+1,i,j))
      enddo
   enddo
enddo

return
END SUBROUTINE fa_preptc

!##############################################################################
Subroutine fluxlimits (m1,m2,m3,ia,iz,ja,jz,scr1,scp,vt3dc,vt3df,vt3dg &
            ,vt3dk,vctr1,vctr2,dn)

implicit none

integer:: m1,m2,m3
integer:: k,i,j,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: scp,scr1
real, dimension(m1,m2,m3) :: vt3dc,vt3df,vt3dg,vt3dk
real, dimension(m1,m2,m3) :: dn
real, dimension(*) :: vctr1, vctr2
real, dimension(:,:,:), allocatable :: F1,Fcorz,s_max,s_min,rp,rm &
   ,scrf1,scr1rho,scrf1rho
real :: pp,qp,pm,qm,cfac 

! Scratch arrays
allocate(F1(m1,m2,m3))
allocate(Fcorz(m1,m2,m3))
allocate(s_max(m1,m2,m3))
allocate(s_min(m1,m2,m3))
allocate(rp(m1,m2,m3))
allocate(rm(m1,m2,m3))
allocate(scrf1(m1,m2,m3))
allocate(scr1rho(m1,m2,m3))
allocate(scrf1rho(m1,m2,m3))

!Flux correction algorithms based on Zalesk (1979) and Durran (1999)

!z advection ---------------------------------------------------

rp = 1.
rm = 1.
do j = ja,jz
   do i = ia,iz

! Compute scalar flux VT3DG

      do k = 1,m1-1
         vt3dg(k,i,j) = vt3dc(k,i,j)  &
            * (vctr1(k) * scr1(k,i,j)  &
            +  vctr2(k) * scr1(k+1,i,j)  &
            +  vt3df(k,i,j) * (scr1(k,i,j) - scr1(k+1,i,j)))
      
! Compute 1st order fluxes
         F1(k,i,j) = (vt3dc(k,i,j)*(scr1(k,i,j)+scr1(k+1,i,j)) &
               - abs(vt3dc(k,i,j))*(scr1(k+1,i,j)-scr1(k,i,j)))*0.5

! Compute anti-diffusive fluxes
         Fcorz(k,i,j) = vt3dg(k,i,j)-F1(k,i,j)
      enddo
   enddo
enddo
!Do monotonic correction
do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
            scrf1(k,i,j) = scr1(k,i,j)  &
            + vt3dk(k,i,j) * (F1(k-1,i,j) - F1(k,i,j)) 

          !Additional optional correction - Zalesk (1979) calls it "cosmetic." Durran (1999)
          !argues it could be more important.
            !if (Fcorz(k,i,j)*(scrf1(k+1,i,j)-scrf1(k,i,j))<0.) then
            !   if (Fcorz(k,i,j)*(scrf1(k+2,i,j)-scrf1(k+1,i,j))<0. .or. &
            !      Fcorz(k,i,j)*(scrf1(k,i,j)-scrf1(k,i-1,j))<0.) &
            !      Fcorz(k,i,j) = 0.
            !endif
     enddo
   scrf1(1,i,j) = scrf1(2,i,j)
   scrf1(m1,i,j) = scrf1(m1-1,i,j)
   enddo
enddo
scr1rho = scr1*dn
scrf1rho = scrf1*dn
do j=ja,jz
   do i=ia,iz
      do k=2,m1-1
         s_max(k,i,j) = max(0.,maxval(scr1rho(k-1:k+1,i,j)), &
                            maxval(scrf1rho(k-1:k+1,i,j)))/dn(k,i,j)
         s_min(k,i,j) = max(0.,min(minval(scr1rho(k-1:k+1,i,j)), &
                            minval(scrf1rho(k-1:k+1,i,j))))/dn(k,i,j)
         !s_max(k,i,j) = max(0.,maxval(scr1rho(k-1:k+1,i,j)), &
         !                   scrf1rho(k,i,j))/dn(k,i,j)
         !s_min(k,i,j) = max(0.,min(minval(scr1rho(k-1:k+1,i,j)), &
         !                   scrf1rho(k,i,j)))/dn(k,i,j)
         pp = max(0.,Fcorz(k-1,i,j)) - min(0.,Fcorz(k,i,j)) * vt3dk(k,i,j) 
         qp = s_max(k,i,j) - scrf1(k,i,j)
         if (pp > 0. .and. qp > 0.) then
            rp(k,i,j) = min(1.,qp/pp)
         else
            rp(k,i,j) = 0.
         endif

         pm = max(0.,Fcorz(k,i,j)) - min(0.,Fcorz(k-1,i,j)) * vt3dk(k,i,j) 
         qm = scrf1(k,i,j) - s_min(k,i,j)
         if (pm > 0. .and. qm > 0.) then
            rm(k,i,j) = min(1., qm/pm)
         else
            rm(k,i,j) = 0.
         endif
      enddo
   rp(1,i,j) = rp(2,i,j)
   rm(1,i,j) = rm(2,i,j)
   enddo
enddo
do j=ja,jz
   do i=ia,iz
      do k=1,m1-1
         if (Fcorz(k,i,j) >= 0.) then
            cfac = min(rp(k+1,i,j),rm(k,i,j))
         else
            cfac = min(rp(k,i,j),rm(k+1,i,j))
         endif
         Fcorz(k,i,j) = Fcorz(k,i,j)*cfac
      enddo
   enddo
enddo
!Final update
do j=ja,jz
   do i=ia,iz
      do k=2,m1-1
         scr1(k,i,j) = scrf1(k,i,j) + vt3dk(k,i,j) * (Fcorz(k-1,i,j) - Fcorz(k,i,j) &
            + scp(k,i,j) * (vt3dc(k,i,j) - vt3dc(k-1,i,j)))
      enddo
   enddo
enddo

! Clean up
deallocate(F1)
deallocate(Fcorz)
deallocate(s_max)
deallocate(s_min)
deallocate(rp)
deallocate(rm)
deallocate(scrf1)
deallocate(scr1rho)
deallocate(scrf1rho)

return
END SUBROUTINE fluxlimits

!##############################################################################
Subroutine fa_xc (m1,m2,m3,ia,iz  &
   ,scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di)

implicit none

integer :: m1,m2,m3,ia,iz,i,j,k
real :: dfact
real, dimension(m1,m2,m3) :: scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di

dfact = .5
do j = 1,m3
   do i = 1,iz

! Compute scalar flux times dtlt [VT3DG]

      do k = 2,m1-1
         vt3dg(k,i,j) = vt3da(k,i,j)  &
            * (vt3di(1,i,j) * scr1(k,i,j)  &
            +  vt3di(2,i,j) * scr1(k,i+1,j)  &
            +  vt3dd(k,i,j) * (scr1(k,i,j) - scr1(k,i+1,j)))
      enddo

! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

      do k = 2,m1-1
         if (vt3da(k,i,j) .gt. 0.) then
            if (vt3dg(k,i,j) * vt3dh(k,i,j) .gt.  &
               dfact * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3da(k,i,j) .lt. 0.) then
            if (-vt3dg(k,i,j) * vt3dh(k,i+1,j) .gt.  &
               dfact * scr1(k,i+1,j)) then
               vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i+1,j)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = 1,m3
   do i = ia,iz
      do k = 2,m1-1
         scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dh(k,i,j) * (vt3dg(k,i-1,j) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3da(k,i,j) - vt3da(k,i-1,j)))
      enddo
   enddo
enddo

return
END SUBROUTINE fa_xc

!##############################################################################
Subroutine fa_yc (m1,m2,m3,ia,iz,ja,jz  &
   ,scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jdim,i,j,k
real :: dfact
real, dimension(m1,m2,m3) :: scp,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di

dfact = .5
do j = 1,jz
   do i = ia,iz

! Compute scalar flux VT3DG

      do k = 2,m1-1
         vt3dg(k,i,j) = vt3db(k,i,j)  &
            * (vt3di(3,i,j) * scr1(k,i,j)  &
            +  vt3di(4,i,j) * scr1(k,i,j+jdim)  &
            +  vt3de(k,i,j) * (scr1(k,i,j) - scr1(k,i,j+jdim)))
      enddo

!      Modify fluxes to retain positive-definiteness on scalar quantities.
!         If a flux will remove 1/2 quantity during a timestep,
!         reduce to first order flux. This will remain positive-definite
!         under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!         both fluxes are evacuating the box.

      do k = 2,m1-1
         if (vt3db(k,i,j) .gt. 0.) then
            if (vt3dg(k,i,j) * vt3dj(k,i,j) .gt.  &
               dfact * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3db(k,i,j) .lt. 0.) then
            if (-vt3dg(k,i,j) * vt3dj(k,i,j+jdim) .gt.  &
               dfact * scr1(k,i,j+jdim)) then
               vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j+jdim)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dj(k,i,j) * (vt3dg(k,i,j-jdim) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3db(k,i,j) - vt3db(k,i,j-jdim)))
      enddo
   enddo
enddo

return
END SUBROUTINE fa_yc

!##############################################################################
Subroutine fa_zc (m1,m2,m3,ia,iz,ja,jz  &
   ,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: dfact
real, dimension(m1,m2,m3) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dk
real, dimension(*) :: vctr1,vctr2

dfact = .5
do j = ja,jz
   do i = ia,iz

! Compute scalar flux VT3DG

      do k = 1,m1-1
         vt3dg(k,i,j) = vt3dc(k,i,j)  &
            * (vctr1(k) * scr1(k,i,j)  &
            +  vctr2(k) * scr1(k+1,i,j)  &
            +  vt3df(k,i,j) * (scr1(k,i,j) - scr1(k+1,i,j)))
      enddo

! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

      do k = 1,m1-1
         if (vt3dc(k,i,j) .gt. 0.) then
            if (vt3dg(k,i,j) * vt3dk(k,i,j) .gt.  &
               dfact * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3dc(k,i,j) .lt. 0.) then
            if (-vt3dg(k,i,j) * vt3dk(k+1,i,j) .gt.  &
               dfact * scr1(k+1,i,j)) then
               vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k+1,i,j)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dk(k,i,j) * (vt3dg(k-1,i,j) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3dc(k,i,j) - vt3dc(k-1,i,j)))
      enddo
   enddo
enddo

return
END SUBROUTINE fa_zc

!##############################################################################
Subroutine advtndc (m1,m2,m3,ia,iz,ja,jz,scp,sca,sct,dtl)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: dtl,dtli
real, dimension(m1,m2,m3) :: scp,sca,sct

dtli = 1. / dtl
do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         sct(k,i,j) = sct(k,i,j) + (sca(k,i,j)-scp(k,i,j)) * dtli
      enddo
   enddo
enddo

return
END SUBROUTINE advtndc
