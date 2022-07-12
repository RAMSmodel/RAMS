!##############################################################################
Subroutine acoustic ()
!--------------------------------------------------------------------
!                  Acoustic terms small time-step driver
!
!     This routine calls all the necessary routines to march the model
!     through the small timesteps.
!-----------------------------------------------------------------------
use mem_tend
use mem_grid
use mem_basic
use node_mod

implicit none

CALL acoust (mzp,mxp,myp &
           ,basic_g(ngrid)%dn0,basic_g(ngrid)%pi0  &
           ,basic_g(ngrid)%th0,basic_g(ngrid)%up  &
           ,basic_g(ngrid)%vp,basic_g(ngrid)%wp  &
           ,basic_g(ngrid)%pp  &
           ,tend%ut,tend%vt,tend%wt,tend%pt  &
           ,grid_g(ngrid)%topt,grid_g(ngrid)%topu  &
           ,grid_g(ngrid)%topv,grid_g(ngrid)%rtgt  &
           ,grid_g(ngrid)%rtgu,grid_g(ngrid)%f13u  &
           ,grid_g(ngrid)%dxu,grid_g(ngrid)%rtgv  &
           ,grid_g(ngrid)%dyv  &
           ,grid_g(ngrid)%f23v,grid_g(ngrid)%f13t  &
           ,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapui  &
           ,grid_g(ngrid)%fmapvi,grid_g(ngrid)%dxt  &
           ,grid_g(ngrid)%dyt,grid_g(ngrid)%fmapt)

return
END SUBROUTINE acoustic

!##############################################################################
Subroutine acoust (m1,m2,m3  &
                 ,dn0,pi0,th0,up,vp,wp,pp,ut,vt,wt,pt  &
                 ,topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv,dyv  &
                 ,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt)
!--------------------------------------------------------------------
!                  Acoustic terms small time-step driver
!
!     This routine calls all the necessary routines to march the model
!     through the small timesteps.
!-----------------------------------------------------------------------
use mem_grid
use node_mod

implicit none

integer :: m1,m2,m3
real, dimension(m1,m2,m3) ::dn0,pi0,th0,up,vp,wp,pp
real, dimension(m2,m3) ::   topt,topu,topv,rtgt,rtgu,f13u,dxu,rtgv,dyv  &
                           ,f23v,f13t,f23t,fmapui,fmapvi,dxt,dyt,fmapt
real, dimension(*) ::      ut,vt,wt,pt

real, dimension(:,:,:), allocatable :: acoc, acof, acog, amoe, &
                                       amof, acoaa, amog
real, dimension(:,:), allocatable :: heatfx1

real :: a1da2

integer :: iter

! Temporary coefficient values
!  coefz() sets these
!  various prdct routines read these
allocate(acoc(m1,m2,m3))
allocate(acof(m1,m2,m3))
allocate(acog(m1,m2,m3))
allocate(amoe(m1,m2,m3))
allocate(amof(m1,m2,m3))
allocate(acoaa(m1,m2,m3))

! Temporary coefficient values
!   prdctw2() sets amog
!   prdctw3() reads amog
allocate(amog(m1,m2,m3))

! Temporary heat flux (at k == 1 level)
!   prdctp1() sets heatfx1
!   prdctw2() reads heatfx1
allocate(heatfx1(m2,m3))

do iter=1,nnacoust(ngrid)

!     Get coefficients for computations

dts = 2. * dtlt / nnacoust(ngrid)

if (iter .eq. 1)  &
   CALL coefz (mzp,mxp,myp,ia,iz,ja,jz  &
      ,acoc,acof,acog,dn0,pi0,th0  &
      ,rtgt,a1da2,amoe,amof,acoaa)

CALL prdctu (mzp,mxp,myp,ia,izu,ja,jz  &
   ,up,ut,pp,th0,f13u,rtgu,rtgt,dxu,topu)

! Update u cyclic boundaries
if (ngrid .eq. 1) CALL update_cyclic (LBC_UP)

CALL prdctv (mzp,mxp,myp,ia,iz,ja,jzv  &
   ,vp,vt,pp,th0,f23v,rtgv,rtgt,dyv,topv)

! Update v cyclic boundaries
if (ngrid .eq. 1) CALL update_cyclic (LBC_VP)

CALL prdctw1 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp,wt,pp,acoc,a1da2,rtgt,topt)

CALL prdctp1 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,pp,up,vp,pi0,dn0,th0,pt,f13t,f23t,rtgt  &
   ,rtgu,rtgv,heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)

CALL prdctw2 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)

CALL prdctw3 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,wp,amog,amoe,impl)

! Update w cyclic boundaries
if (ngrid .eq. 1) CALL update_cyclic (LBC_WP)

CALL prdctp2 (mzp,mxp,myp,ia,iz,ja,jz  &
   ,pp,wp,acof,acog)

! Update p cyclic boundaries
if (ngrid .eq. 1) CALL update_cyclic (LBC_PP)

enddo

deallocate(acoc)
deallocate(acof)
deallocate(acog)
deallocate(amoe)
deallocate(amof)
deallocate(acoaa)

deallocate(amog)

deallocate(heatfx1)

return
END SUBROUTINE acoust

!##############################################################################
Subroutine prdctu (m1,m2,m3,ia,iz,ja,jz  &
   ,up,ut,pp,th0,f13u,rtgu,rtgt,dxu,topu)

use mem_grid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: dxl
real, dimension(m1,m2,m3) :: up,ut,th0,pp
real, dimension(m2,m3) ::  f13u,rtgu,rtgt,dxu,topu
real :: dummy
real, dimension(:,:,:), allocatable :: dpdx, dpdx_top

!     U prediction

allocate(dpdx(m1,m2,m3))
allocate(dpdx_top(m1,m2,m3))

CALL azero (m1*m2*m3,dpdx)

!     Calculate acoustic tendency (horizontal pressure gradient)

do j = ja,jz
   do i = ia,iz
      do k = 1,m1-1
         dpdx_top(k,i,j) = (pp(k,i,j) + pp(k+1,i,j)  &
            + pp(k,i+1,j) + pp(k+1,i+1,j)) * hw4(k)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      dxl = dxu(i,j) / rtgu(i,j)
      do k = 2,m1-1
         dpdx(k,i,j) = -(th0(k,i,j) + th0(k,i+1,j)) * .5  &
            * ((pp(k,i+1,j) * rtgt(i+1,j) - pp(k,i,j) * rtgt(i,j)) * dxl  &
            + (dpdx_top(k,i,j) - dpdx_top(k-1,i,j)) * dzt(k) * f13u(i,j))
      enddo
   enddo
enddo

if (distim .ne. 0.) then
   ! 3rd to last argument is theta tendency which is only used
   ! during rayleigh friction calculation for theta (1st arg == 4)
   CALL rayf (1,m1,m2,m3,ia,iz,ja,jz,up,th0,dummy,rtgu,topu)
endif

do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         up(k,i,j) = up(k,i,j) + dts * (dpdx(k,i,j) + ut(k,i,j))
      enddo
   enddo
enddo

if (nstbot .eq. 1 .and. itopo .eq. 1)  &
     CALL botset (m1,m2,m3,up,'U')

deallocate(dpdx)
deallocate(dpdx_top)

! Update u overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,up,3)

return
END SUBROUTINE prdctu

!##############################################################################
Subroutine prdctv (m1,m2,m3,ia,iz,ja,jz  &
   ,vp,vt,pp,th0,f23v,rtgv,rtgt,dyv,topv)
   
use mem_grid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k

real :: dyl
real, dimension (m1,m2,m3) :: vp,vt,th0,pp
real, dimension(m2,m3) :: f23v,rtgv,rtgt,dyv,topv
real :: dummy
real, dimension(:,:,:), allocatable :: dpdy, dpdy_top

!     V prediction

if (jdim .eq. 1) then

   allocate(dpdy(m1,m2,m3))
   allocate(dpdy_top(m1,m2,m3))

   CALL azero (m1*m2*m3,dpdy)

!       calculate acoustic tendency (horizontal pressure gradient)

   do j = ja,jz
      do i = ia,iz
         do k = 1,m1-1
            dpdy_top(k,i,j) =(pp(k,i,j) + pp(k+1,i,j)  &
               + pp(k,i,j+1) + pp(k+1,i,j+1)) * hw4(k)
         enddo
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         dyl = dyv(i,j) / rtgv(i,j)
         do k = 2,m1-1
            dpdy(k,i,j) = -(th0(k,i,j) + th0(k,i,j+1)) * .5  &
               * ((pp(k,i,j+1) * rtgt(i,j+1) - pp(k,i,j) * rtgt(i,j)) * dyl  &
               + (dpdy_top(k,i,j) - dpdy_top(k-1,i,j)) * dzt(k) * f23v(i,j))
         enddo
      enddo
   enddo

   if (distim .ne. 0.) then
      ! 3rd to last argument is theta tendency which is only used
      ! during rayleigh friction calculation for theta (1st arg == 4)
      CALL rayf (2,m1,m2,m3,ia,iz,ja,jz,vp,th0,dummy,rtgv,topv)
   endif

   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            vp(k,i,j) = vp(k,i,j) + dts * (dpdy(k,i,j) + vt(k,i,j))
         enddo
      enddo
   enddo

   if (nstbot .eq. 1 .and. itopo .eq. 1)  &
        CALL botset (m1,m2,m3,vp,'V')

   deallocate(dpdy)
   deallocate(dpdy_top)

   ! Update v overlap regions
   if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,vp,3)

endif

return
END SUBROUTINE prdctv

!##############################################################################
Subroutine prdctw1 (m1,m2,m3,ia,iz,ja,jz  &
   ,wp,wt,pp,acoc,a1da2,rtgt,topt)

use mem_grid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: a1da2
real, dimension(m1,m2,m3) :: wp,wt,pp,acoc
real, dimension(m2,m3) :: rtgt,topt
real :: dummy

!     First part of prediction at I,J point

!     Compute forward part of Crank-Nickelson scheme. This will be total
!     W prediction for explicit case.

if (distim .ne. 0.) then
   ! 3rd to last argument is theta tendency which is only used
   ! during rayleigh friction calculation for theta (1st arg == 4)
   CALL rayf (3,m1,m2,m3,ia,iz,ja,jz,wp,wp,dummy,rtgt,topt)
endif

!      do j=ja,jz
!         do i=ia,iz

do j = 1,m3
   do i = 1,m2
      do k = 1,m1-2
         !This is Wk* in RAMS tech manual
         wp(k,i,j) = wp(k,i,j) + dts * wt(k,i,j)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      do k = 1,m1-2
         !This is dk(PIk') in RAMS tech manual (dk = dt * theta0 / dz)
         wp(k,i,j) = wp(k,i,j) + a1da2 * acoc(k,i,j) * (pp(k,i,j)-pp(k+1,i,j))
      enddo
   enddo
enddo

! Update w overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,wp,3)

return
END SUBROUTINE prdctw1

!##############################################################################
Subroutine prdctw2 (m1,m2,m3,ia,iz,ja,jz  &
   ,wp,pp,acoc,amof,amog,acoaa,rtgt,heatfx1)

use mem_grid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: wp,pp,acoc,amof,amog,acoaa
real, dimension(m2,m3) :: heatfx1,rtgt(m2,m3)

if (nstbot .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         wp(1,i,j) = -heatfx1(i,j) * rtgt(i,j)
      enddo
   enddo
endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    new trial fix 1/11/96

if (nsttop .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         wp(nz,i,j) = 0.
      enddo
   enddo
endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

if (impl .eq. 1) then

!         First implicit part of the w prediction

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-2
            wp(k,i,j) = wp(k,i,j) - (pp(k+1,i,j) - pp(k,i,j)) * acoc(k,i,j)
         enddo
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         amog(1,i,j) = -wp(1,i,j) / amof(1,i,j)
      enddo
      do k = 2,m1-2
         do i = ia,iz
            amog(k,i,j) = (-wp(k,i,j) - acoaa(k,i,j) * amog(k-1,i,j))  &
               / amof(k,i,j)
         enddo
      enddo
   enddo

endif

! Update w overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,wp,3)

return
END SUBROUTINE prdctw2

!##############################################################################
Subroutine prdctw3 (m1,m2,m3,ia,iz,ja,jz,wp,amog,amoe,impl)

use mem_grid, only:ngrid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,impl,i,j,k,isflag,vnum

real, dimension(m1,m2,m3) :: wp,amog,amoe

!     Conclusion of implicit w prediction

if (impl .eq. 1) then
   do k = m1-2,2,-1
      do j = ja,jz
         do i = ia,iz
            wp(k,i,j) = amog(k,i,j) - amoe(k,i,j) * wp(k+1,i,j)
         enddo
      enddo
   enddo
endif

! Update w overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,wp,3)

return
END SUBROUTINE prdctw3

!##############################################################################
Subroutine prdctp1 (m1,m2,m3,ia,iz,ja,jz  &
   ,pp,up,vp,pi0,dn0,th0,pt,f13t,f23t,rtgt,rtgu,rtgv  &
   ,heatfx1,fmapui,fmapvi,dxt,dyt,fmapt)

use mem_grid
use node_mod, only:nmachs
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real :: rocvpct
real, dimension(m1,m2,m3) :: pp,up,vp,pi0,pt,dn0,th0
real, dimension(m2,m3) :: f13t,f23t,rtgt,heatfx1,rtgu,rtgv,fmapui  &
   ,fmapvi,dxt,dyt,fmapt
real, dimension(:,:,:), allocatable :: heatdv, heatfx

allocate(heatdv(m1,m2,m3))
allocate(heatfx(m1,m2,m3))

CALL azero (m1*m2*m3,heatdv)
rocvpct =rocv *sspct ** 2

!     Divergence calculations for topographical transformation

!     First calculate vertically transformed heat flux

do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         heatfx(k,i,j) = ((up(k,i,j) + up(k,i-1,j)) * f13t(i,j)  &
            + (vp(k,i,j) + vp(k,i,j-jdim)) * f23t(i,j)  &
            ) * dn0(k,i,j) * th0(k,i,j)
      enddo
   enddo
enddo
do j = ja,jz
   do i = ia,iz
      do k = 1,m1-1
         heatfx(k,i,j) = (heatfx(k,i,j) + heatfx(k+1,i,j)) * hw4(k)
      enddo
      heatfx1(i,j) = heatfx(1,i,j) / (.5 * (dn0(1,i,j) * th0(1,i,j)  &
         + dn0(2,i,j) * th0(2,i,j)))
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         heatdv(k,i,j) = (heatfx(k,i,j) - heatfx(k-1,i,j)) * dzt(k)
      enddo
   enddo
enddo
do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         heatfx(k,i,j) = dn0(k,i,j) * th0(k,i,j)
      enddo
   enddo
enddo
do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         heatdv(k,i,j) = -rocvpct * pi0(k,i,j) / heatfx(k,i,j)  &
            * (heatdv(k,i,j) + fmapt(i,j)  &
            * ((up(k,i,j) * rtgu(i,j) * fmapui(i,j)  &
            * (heatfx(k,i,j) + heatfx(k,i+1,j))  &
            - up(k,i-1,j) * rtgu(i-1,j) * fmapui(i-1,j)  &
            * (heatfx(k,i,j) + heatfx(k,i-1,j))) * dxt(i,j) * .5  &
            + (vp(k,i,j) * rtgv(i,j) * fmapvi(i,j)  &
            * (heatfx(k,i,j) + heatfx(k,i,j+jdim))  &
            - vp(k,i,j-jdim) * rtgv(i,j-jdim)  &
            * fmapvi(i,j-jdim)  &
            * (heatfx(k,i,j) + heatfx(k,i,j-jdim)))  &
            * dyt(i,j) * .5) / rtgt(i,j))
      enddo
   enddo
enddo


do j = ja,jz
   do i = ia,iz
      do k = 1,m1
         pp(k,i,j) = pp(k,i,j) + (pt(k,i,j) + heatdv(k,i,j)) * dts
      enddo
   enddo
enddo

deallocate(heatdv)
deallocate(heatfx)

! Update p overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,pp,3)

return
END SUBROUTINE prdctp1

!##############################################################################
Subroutine prdctp2 (m1,m2,m3,ia,iz,ja,jz,pp,wp,acof,acog)

use mem_grid
use node_mod, only:nmachs

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: pp,wp,acof,acog

!           Finish pressure prediction

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         pp(k,i,j) = pp(k,i,j)  &
            + (wp(k,i,j) * acof(k,i,j) + wp(k-1,i,j) * acog(k,i,j))
      enddo
   enddo
enddo

if (nstbot .eq. 1) CALL botset (m1,m2,m3,pp,'P')

! Update p overlap regions
if (nmachs .gt. 1) CALL update_lbc_var (m1,m2,m3,ngrid,pp,3)

return
END SUBROUTINE prdctp2

!##############################################################################
Subroutine coefz (m1,m2,m3,ia,iz,ja,jz  &
   ,acoc,acof,acog,dn0,pi0,th0,rtgt,a1da2,amoe,amof,acoaa)

use mem_grid
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k

real :: dt2al2,a1da2,rdto2cv,dt2al2r,rdtr
real, dimension(m1,m2,m3) :: th0,pi0,dn0,acoc,acof,acog,amoe,amof,acoaa
real, dimension(m2,m3) :: rtgt
real, dimension(:), allocatable :: acobb,acocc,mconn_coeff,rho_theta_base

allocate(acobb(m1))
allocate(acocc(m1))
allocate(mconn_coeff(m1))
allocate(rho_theta_base(m1))

! +--------------------------------------------------------------------+
! \   Calculate coefficients for the vertical pressure gradient        \
! \     and divergence terms.  These will be combined later for the    \
! \     implicit computations.                                         \
! +--------------------------------------------------------------------+

if (impl .eq. 1) then
   dt2al2 = dts * .75
   a1da2 = 1. / 3.
else
   dt2al2 = dts
   a1da2 = 1.
endif
rdto2cv = sspct ** 2 * rgas * dts / (2.0 * cv)

do j = ja,jz
   do i = ia,iz

!         Coefficient for the vertical pressure gradient term

      dt2al2r = .5 * dt2al2 / rtgt(i,j)
      do k = 1,m1-1
         acoc(k,i,j) = dt2al2r * dzm(k) * (th0(k,i,j) + th0(k+1,i,j))
      enddo

!         Coefficients for the vertical divergence term

      rdtr = rdto2cv / rtgt(i,j)
      do k = 2,m1
         rho_theta_base(k) = dn0(k,i,j) * th0(k,i,j)
         mconn_coeff(k) = rdtr * pi0(k,i,j) * dzt(k) / rho_theta_base(k)
      enddo
      rho_theta_base(1) = dn0(1,i,j) * th0(1,i,j)
      do k = 2,m1-1
         acof(k,i,j) = -mconn_coeff(k) * (rho_theta_base(k) + rho_theta_base(k+1))
         acog(k,i,j) = mconn_coeff(k) * (rho_theta_base(k) + rho_theta_base(k-1))
      enddo
      acog(m1,i,j) = mconn_coeff(m1) * (rho_theta_base(m1) + rho_theta_base(nz))

      do k = 2,m1-1
         acoaa(k,i,j) = acoc(k,i,j) * acog(k,i,j)
         acobb(k) = acoc(k,i,j) * (acof(k,i,j) - acog(k+1,i,j)) - 1.
         acocc(k) = -acoc(k,i,j) * acof(k+1,i,j)
      enddo
      acobb(1) = -1.
      acocc(1) = 0.
      acoaa(m1,i,j) = 0.
      acobb(m1) = -1.

      amof(1,i,j) = acobb(1)
      amoe(1,i,j) = acocc(1) / amof(1,i,j)
      do k = 2,m1
         amof(k,i,j) = acobb(k) - acoaa(k,i,j) * amoe(k-1,i,j)
         amoe(k,i,j) = acocc(k) / amof(k,i,j)
      enddo

   enddo
enddo

deallocate(acobb)
deallocate(acocc)
deallocate(mconn_coeff)
deallocate(rho_theta_base)

return
END SUBROUTINE coefz

!##############################################################################
Subroutine buoyancy ()

use mem_tend
use mem_basic
use mem_grid
use node_mod
use micphys

implicit none

CALL boyanc (mzp,mxp,myp,ia,iz,ja,jz,tend%wt &
   ,basic_g(ngrid)%theta,basic_g(ngrid)%rtp  &
   ,basic_g(ngrid)%rv,basic_g(ngrid)%th0,dtlt)

return
END SUBROUTINE buoyancy

!##############################################################################
Subroutine boyanc (m1,m2,m3,ia,iz,ja,jz,wt,theta,rtc,rv,th0,dtlt)

use rconstants
use micphys
use mem_basic
use mem_grid, only:ngrid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m1,m2,m3) :: wt,theta,rtc,rv,th0
real, dimension(:,:,:), allocatable :: vtemp,wpbuoytheta,wpbuoycond,wpadvdif
real :: dtlt

allocate(vtemp(m1,m2,m3))

if(imbudget>=1) then
 allocate(wpbuoytheta(m1,m2,m3))
 allocate(wpbuoycond(m1,m2,m3))
 allocate(wpadvdif(m1,m2,m3))
 CALL azero (m1*m2*m3,wpbuoytheta)
 CALL azero (m1*m2*m3,wpbuoycond)
 CALL azero (m1*m2*m3,wpadvdif)
endif

if (level .ge. 1) then
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vtemp(k,i,j) = gg * ((theta(k,i,j) * (1. + .61 * rv(k,i,j))  &
               - th0(k,i,j)) / th0(k,i,j) - (rtc(k,i,j) - rv(k,i,j)) )
            !calculate partial W buoyancy budgets
            if(imbudget>=1) then
              wpbuoytheta(k,i,j) = gg * ((theta(k,i,j)*(1.+.61*rv(k,i,j)) &
                - th0(k,i,j)) / th0(k,i,j))
              wpbuoycond(k,i,j)  = gg * (-1.0*(rtc(k,i,j) - rv(k,i,j)))
            endif
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vtemp(k,i,j) = gg * (theta(k,i,j) / th0(k,i,j) - 1.)
            !calculate partial W buoyancy budgets
            if(imbudget>=1) wpbuoytheta(k,i,j) = vtemp(k,i,j)
         enddo
      enddo
   enddo
endif

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-2
         !Calculate W buoyancy budgets (m/s)
         !Calculate W due to advection and diffusion (current wt)
         !Multiply by 2*dt for leapfrog timestep t-dt to t+dt
         if(imbudget>=1) then         
           wpadvdif(k,i,j)    = 2.0 * dtlt * wt(k,i,j)
           wpbuoytheta(k,i,j) = 2.0 * dtlt * (wpbuoytheta(k,i,j) &
                                            + wpbuoytheta(k+1,i,j))
           wpbuoycond(k,i,j)  = 2.0 * dtlt * (wpbuoycond(k,i,j) &
                                            + wpbuoycond(k+1,i,j))
         endif
         wt(k,i,j) = wt(k,i,j) + vtemp(k,i,j) + vtemp(k+1,i,j)
      enddo
   enddo
enddo

deallocate(vtemp)

!Copy local W buoyancy budgets (m/s) to global budget variables
if(imbudget>=1) then
 basic_g(ngrid)%wp_buoy_theta = wpbuoytheta
 basic_g(ngrid)%wp_buoy_cond  = wpbuoycond
 basic_g(ngrid)%wp_advdif     = wpadvdif
 deallocate(wpbuoytheta)
 deallocate(wpbuoycond)
 deallocate(wpadvdif)
endif

return
END SUBROUTINE boyanc
