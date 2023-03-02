!##############################################################################
Subroutine diffuse ()

! +-----------------------------------------------------------------+
! \     this routine is the subdriver to compute tendencies due to  \
! \       subgrid-scale turbulence.                                 \
! +-----------------------------------------------------------------+

use mem_tend
use mem_basic
use var_tables
use mem_turb
use mem_grid
use mem_leaf
use mem_micro
use mem_scratch
use node_mod
use micphys
use micro_prm, only:krdrop

implicit none

integer :: ind,n
real :: s1,s2,s3
real, pointer :: scalarp, scalart
integer :: i,j,k,ksf,izero

! ldgrant(2016): added lots of comments and notes
! From strain, if idiffk = 1 or 2:
!  vt3da=du/dx  vt3db=dv/dx  vt3dc=dv/dy  vt3dd=du/dz
!  vt3de=dv/dz  vt3df=dw/dx  vt3dg=dw/dy  
!  vt3dh = 2*[ (du/dx)^2 + (dv/dy)^2 ] + (dv/dx + du/dy)^2
!   = horizontal strain rate tensor magnitude squared 
!  vt3di = (du/dz)^2 + (dv/dz)^2 = vertical wind shear magnitude squared
!  vt3dn=du/dy  
!  scr2 is not set to anything
! if idiffk >= 3:
!  vt3dh = vt3di = 3D strain rate tensor magnitude squared 
!  scr2 = dw/dz
CALL strain (mzp,mxp,myp,ia,iz,ja,jz                          &
      ,ia_1,ja_1,iz1,jz1,jdim                                &
      ,basic_g(ngrid)%up (1,1,1) ,basic_g(ngrid)%vp (1,1,1)  &
      ,basic_g(ngrid)%wp (1,1,1) ,scratch%vt3da     (1)      &
      ,scratch%vt3db     (1)     ,scratch%vt3dc     (1)      &
      ,scratch%vt3dd     (1)     ,scratch%vt3de     (1)      &
      ,scratch%vt3df     (1)     ,scratch%vt3dg     (1)      &
      ,scratch%vt3dh     (1)     ,scratch%vt3di     (1)      &
      ,scratch%vt3dn     (1)     ,scratch%scr2      (1)      &
      ,idiffk(ngrid))

if (level <= 1) CALL azero (mxyzp,scratch%vt3dp(1))
! vt3dp set equal to rcp
if (level == 2 .or. level == 3) &
  CALL ae1 (mxyzp,scratch%vt3dp(1),micro_g(ngrid)%rcp(1,1,1))
izero=0
if (level == 4) CALL sum_bins (micro_g(ngrid)%ffcd,scratch%vt3dp, &
                              mzp,mxp,myp,1,krdrop-1,izero)

! vt3dj = N^2 = moist brunt vaisala frequency
CALL bruvais (mzp,mxp,myp,ia,iz,ja,jz                          &
   ,basic_g(ngrid)%theta (1,1,1) ,basic_g(ngrid)%rtp (1,1,1)  &
   ,basic_g(ngrid)%rv    (1,1,1) ,scratch%vt3dp(1)            &
   ,basic_g(ngrid)%pp    (1,1,1) ,basic_g(ngrid)%pi0 (1,1,1)  &
   ,scratch%vt3dj        (1)     ,grid_g(ngrid)%rtgt (1,1))

! vt3dh is passed in as the (horizontal or 3D) sqd magnitude of 
!   strain rate tensor
! vt3dh = Khz = vertical diffusion coefficient for scalars on output
!   if idiffk = 2 or 3 (Smagorinsky)
! vt3di = vertical wind shear or 3D strain tensor magnitude sqd from routine strain
! vt3dj = N^2 = moist brunt vaisala frequency from subroutine bruvais
! vt3dk = Richardson #, calculated here
! scr1 = Kmz = vertical diffusion coef for momentum vars for idiffk = 2, 3
! scr2 = Kmx = horiz diffusion coef for momentum vars, set for all idiffk = 1,2,3
if (idiffk(ngrid) <= 3) then
   CALL mxdefm (mzp,mxp,myp,ia,iz,ja,jz                      &
      ,scratch%vt3dh      (1)     ,scratch%vt3di      (1)    &
      ,scratch%vt3dj      (1)     ,scratch%vt3dk      (1)    &
      ,scratch%scr1       (1)     ,scratch%scr2       (1)    &
      ,basic_g(ngrid)%dn0 (1,1,1) ,grid_g(ngrid)%rtgt (1,1)  &
      ,grid_g(ngrid)%dxt  (1,1))
endif

if (idiffk(ngrid) == 1) then
   CALL tkemy (mzp,mxp,myp,ia,iz,ja,jz,jdim  &
      ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
      ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
      ,scratch%vt3dj        (1)     ,scratch%scr1         (1)      &
      ,grid_g(ngrid)%rtgt   (1,1)   ,basic_g(ngrid)%theta (1,1,1)  &
      ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%up    (1,1,1)  &
      ,basic_g(ngrid)%vp    (1,1,1) ,basic_g(ngrid)%wp    (1,1,1)  &
      ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
      ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1),vctr34)
endif

if (idiffk(ngrid) == 4) then
   CALL mxtked (mzp,mxp,myp,ia,iz,ja,jz  &
      ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
      ,scratch%vt3da        (1)     ,scratch%vt3dc        (1)      &
      ,scratch%vt3dh        (1)     ,scratch%vt3dj        (1)      &
      ,scratch%scr1         (1)     ,scratch%scr2         (1)      &
      ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%rtgt   (1,1))
endif

 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%scr1 (1),basic_g(ngrid)%dn0(1,1,1))
 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%scr2 (1),basic_g(ngrid)%dn0(1,1,1))
 CALL klbnd (mzp,mxp,myp,ibcon,jdim  &
    ,scratch%vt3dh(1),basic_g(ngrid)%dn0(1,1,1))

!bob  swap new hkm, vkm, and vkh with past time level:  lagged K's have
!bob  internal lateral boundary values from neighboring nodes
! this is for numerical stability

! vt3dh = Khz = vkh [kg/m3 * m2/s]
! scr1 = Kmz = vkm
! scr2 = Kmx = hkm
ind = 0
do j = 1,mmyp(ngrid)
   do i = 1,mmxp(ngrid)
      do k = 1,mmzp(ngrid)
         ind = ind + 1
         s1 = scratch%scr2(ind) ! hkm
         s2 = scratch%scr1(ind) ! vkm
         s3 = scratch%vt3dh(ind) ! vkh
         scratch%scr2(ind) = turb_g(ngrid)%hkm(k,i,j)
         scratch%scr1(ind) = turb_g(ngrid)%vkm(k,i,j)
         scratch%vt3dh(ind) = turb_g(ngrid)%vkh(k,i,j)
!! also for vt3di = K(tke) ?????    22 March 02
!!         scratch%vt3di(ind) = turb_g(ngrid)%vke(k,i,j)
         turb_g(ngrid)%hkm(k,i,j) = s1 ! now these are current time level K's
         turb_g(ngrid)%vkm(k,i,j) = s2
         turb_g(ngrid)%vkh(k,i,j) = s3
      enddo
   enddo
enddo

CALL diffvel (mzp,mxp,myp,ia,iz,ja,jz,jdim                         &
      ,iz1,jz1,izu,jzv,idiffk(ngrid)                               &
      ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
      ,basic_g(ngrid)%wp    (1,1,1) ,tend%ut              (1)      &
      ,tend%vt              (1)     ,tend%wt              (1)      &
      ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
      ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
      ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
      ,scratch%vt3dg        (1)     ,scratch%vt3dj        (1)      &
      ,scratch%vt3dk        (1)     ,scratch%vt3dl        (1)      &
      ,scratch%vt3dm        (1)     ,scratch%vt3dn        (1)      &
      ,scratch%vt3do        (1)     ,grid_g(ngrid)%rtgu   (1,1)    &
      ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
      ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
      ,turb_g(ngrid)%sflux_w(1,1)   ,basic_g(ngrid)%dn0   (1,1,1)  &
      ,basic_g(ngrid)%dn0u  (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
      ,scratch%scr1         (1)     ,scratch%scr2         (1)      )

! Convert momentum K's to scalar K's, if necessary

if (idiffk(ngrid) <= 3) then
   ! now scr2 = Khx = hkh
   do ind = 1,mxyzp
      scratch%scr2(ind) = scratch%scr2(ind) * xkhkm(ngrid)
   enddo
elseif (idiffk(ngrid) == 4) then
   do ind = 1,mxyzp
      scratch%vt3di(ind) = 2. * scratch%scr1(ind)
   enddo
endif

do n = 1,num_scalar(ngrid)

   scalarp => scalar_tab(n,ngrid)%var_p ! scalar
   scalart => scalar_tab(n,ngrid)%var_t ! scalar tendency

   CALL azero (mxp*myp,scratch%vt2da(1)) ! vt2da: zeroed out
   if (nstbot == 1) then
      if (scalar_tab(n,ngrid)%name == 'THP') then
         ! vt2da: include sflux_t
         CALL atob (mxp*myp,turb_g(ngrid)%sflux_t(1,1),scratch%vt2da(1))
      elseif (scalar_tab(n,ngrid)%name == 'RTP') then
         ! vt2da: include sflux_r
         CALL atob (mxp*myp,turb_g(ngrid)%sflux_r(1,1),scratch%vt2da(1))
      endif
   endif

! 3/10/01 - Define ksf below, the "K scalar flag", to let routine diffsclr
! know which vertical K is being passed to it.  If diffsclr sees that it's
! a different K from the previous one, diffsclr will re-compute the tridiff
! matrix coefficients.  In order to use vertical scalar K's other than
! vt3dh and vt3di, use ksf = 3, ksf = 4, etc. for each different K.


   ! Need to adjust the vkkh and hkkh arguments of diffsclr() according
   ! to which scheme we are using. These are the last two arguments:
   !
   !    diffsclr(..., vkkh, hkkh)
   !
   if (scalar_tab(n,ngrid)%name == 'TKEP') then
      ksf = 1   
      if (idiffk(ngrid) == 4) then
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%vt3di(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf           &
         ,scalarp,scalart             ,scratch%vt3da(1)             &
         ,scratch%vt3db      (1)      ,scratch%vt3df      (1)       &
         ,scratch%vt3dg      (1)      ,scratch%vt3dj      (1)       &
         ,scratch%vt3dk      (1)      ,scratch%vt3do      (1)       &
         ,scratch%vt3dc      (1)      ,scratch%vt3dd      (1)       &
         ,scratch%vt3dl      (1)      ,scratch%vt3dm      (1)       &
         ,scratch%vt2db      (1)      ,grid_g(ngrid)%rtgt (1,1)     &
         ,scratch%vt2da      (1)      ,basic_g(ngrid)%dn0 (1,1,1)   &
         ,basic_g(ngrid)%th0 (1,1,1)  ,basic_g(ngrid)%rvt0 (1,1,1)  &
         ,basic_g(ngrid)%th00 (1,1,1) ,basic_g(ngrid)%rvt00 (1,1,1) &
         ,scratch%vt3di      (1)      ,scratch%vt3di      (1)      )
      else
        ! vkkh <- scratch%vt3di(1)
        ! hkkh <- scratch%scr2(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf           &
         ,scalarp,scalart             ,scratch%vt3da(1)             &
         ,scratch%vt3db      (1)      ,scratch%vt3df      (1)       &
         ,scratch%vt3dg      (1)      ,scratch%vt3dj      (1)       &
         ,scratch%vt3dk      (1)      ,scratch%vt3do      (1)       &
         ,scratch%vt3dc      (1)      ,scratch%vt3dd      (1)       &
         ,scratch%vt3dl      (1)      ,scratch%vt3dm      (1)       &
         ,scratch%vt2db      (1)      ,grid_g(ngrid)%rtgt (1,1)     &
         ,scratch%vt2da      (1)      ,basic_g(ngrid)%dn0 (1,1,1)   &
         ,basic_g(ngrid)%th0 (1,1,1)  ,basic_g(ngrid)%rvt0 (1,1,1)  &
         ,basic_g(ngrid)%th00 (1,1,1) ,basic_g(ngrid)%rvt00 (1,1,1) &
         ,scratch%vt3di      (1)      ,scratch%scr2       (1)      )
      endif
   else ! scalars other than TKEP
      ksf = 2
      if (idiffk(ngrid) == 4) then
        ! vkkh <-  scratch%vt3dh(1)
        ! hkkh <-  scratch%vt3dh(1)
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf           &
         ,scalarp,scalart             ,scratch%vt3da(1)             &
         ,scratch%vt3db      (1)      ,scratch%vt3df      (1)       &
         ,scratch%vt3dg      (1)      ,scratch%vt3dj      (1)       &
         ,scratch%vt3dk      (1)      ,scratch%vt3do      (1)       &
         ,scratch%vt3dc      (1)      ,scratch%vt3dd      (1)       &
         ,scratch%vt3dl      (1)      ,scratch%vt3dm      (1)       &
         ,scratch%vt2db      (1)      ,grid_g(ngrid)%rtgt (1,1)     &
         ,scratch%vt2da      (1)      ,basic_g(ngrid)%dn0 (1,1,1)   &
         ,basic_g(ngrid)%th0 (1,1,1)  ,basic_g(ngrid)%rvt0 (1,1,1)  &
         ,basic_g(ngrid)%th00 (1,1,1) ,basic_g(ngrid)%rvt00 (1,1,1) &
         ,scratch%vt3dh      (1)      ,scratch%vt3dh      (1)      )
      else ! This is the call for idiffk = 1, 2, or 3
        ! ksf = a flag to tell whether to compute vertical diffusion matrix coeffs
        ! scalarp = scalar variable 
        ! scalart = scalar tendency
        ! vt3da = [something from diffvel call], flux divergence x-dir
        ! vt3db = [something from diffvel call], flux divergence y-dir
        ! vt3df = 
        ! vt3dg = 
        ! for the next 3 vars: if rearrange the vertical part of the scalar 
        !  tendency eqn so LHS has scp(k) at time t+dt, then the following are 
        !  coefficients: 
        ! vt3dj = coef on scp(k-1) [kg/m3]
        ! vt3dk = coef on scp(k+1) [kg/m3]
        ! vt3do = coef on scp(k) [kg/m3]
        ! vt3dc = rho0*scalarp [kg/m3 * scalar units]
        ! vt3dd = 1/rho0 [m3/kg]
        ! vt3dl = 
        ! vt3dm =
        ! vt2db = 1/[coef on scp(k) at k=2] [m3/kg]
        ! vt2da = sflux_r if scalar = RTP or sflux_t if scalar = THP or 0 otherwise
        ! vt3dh = vkh at past time level [kg/m3 * m2/s]
        ! scr2 = hkh at past time level [kg/m3 * m2/s]
        ! scalart is updated here
        CALL diffsclr (mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf           &
         ,scalarp,scalart             ,scratch%vt3da(1)             &
         ,scratch%vt3db      (1)      ,scratch%vt3df      (1)       &
         ,scratch%vt3dg      (1)      ,scratch%vt3dj      (1)       &
         ,scratch%vt3dk      (1)      ,scratch%vt3do      (1)       &
         ,scratch%vt3dc      (1)      ,scratch%vt3dd      (1)       &
         ,scratch%vt3dl      (1)      ,scratch%vt3dm      (1)       &
         ,scratch%vt2db      (1)      ,grid_g(ngrid)%rtgt (1,1)     &
         ,scratch%vt2da      (1)      ,basic_g(ngrid)%dn0 (1,1,1)   &
         ,basic_g(ngrid)%th0 (1,1,1)  ,basic_g(ngrid)%rvt0 (1,1,1)  &
         ,basic_g(ngrid)%th00 (1,1,1) ,basic_g(ngrid)%rvt00 (1,1,1) &
         ,scratch%vt3dh      (1)      ,scratch%scr2       (1)      )
      endif
   endif

enddo

return
END SUBROUTINE diffuse

!##############################################################################
Subroutine strain (m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1  &
   ,jd,up,vp,wp,dudx,dvdx,dvdy,dudz,dvdz  &
   ,dwdx,dwdy,straintens1,straintens2,dudy,dwdz,idiffk)

! ldgrant(2016): changed scratch variable names and added lots of comments

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jd,idiffk,i,j,k
real, dimension(m1,m2,m3) :: up,vp,wp,dudx,dvdx,dvdy,dudz,dvdz,dwdx &
   ,dwdy,straintens1,straintens2,dudy,dwdz

CALL grad (m1,m2,m3,ia,iz1,ja,jz,up,dudx,'XDIR','UPNT') ! du/dx
CALL grad (m1,m2,m3,ia_1,iz,ja_1,jz,vp,dvdx,'XDIR','VPNT') ! dv/dx
CALL grad (m1,m2,m3,ia_1,iz,ja,jz,wp,dwdx,'XDIR','WPNT') ! dw/dx

CALL grad (m1,m2,m3,ia_1,iz,ja_1,jz,up,dudy,'YDIR','UPNT') ! du/dy
CALL grad (m1,m2,m3,ia,iz,ja,jz1,vp,dvdy,'YDIR','VPNT') ! dv/dy
CALL grad (m1,m2,m3,ia,iz,ja_1,jz,wp,dwdy,'YDIR','WPNT') ! dw/dy

CALL grad (m1,m2,m3,ia_1,iz,ja,jz,up,dudz,'ZDIR','UPNT') ! du/dz
CALL grad (m1,m2,m3,ia,iz,ja_1,jz,vp,dvdz,'ZDIR','VPNT') ! dv/dz
if(idiffk.ge.3)then 
   ! dwdz used here only for dx~dz turb schemes
   CALL grad (m1,m2,m3,ia,iz,ja,jz,wp,dwdz,'ZDIR','WPNT') ! dw/dz
endif

if (idiffk .le. 2) then
   ! code for idiffk=1,2 when dx>dz
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            ! horiz strain rate tensor magnitude sqd =
            ! 2*[ (du/dx)^2 + (dv/dy)^2 ] + (dv/dx + du/dy)^2
            straintens1(k,i,j) =2. * (dudx(k,i,j) * dudx(k,i,j)  & ! 2*(du/dx)^2
                  + dvdy(k,i,j) * dvdy(k,i,j))  & ! 2*(dv/dy)^2
               + .0625 * (dvdx(k,i,j) + dvdx(k,i-1,j)  &
                  + dvdx(k,i,j-jd) + dvdx(k,i-1,j-jd)  & ! (avg dv/dx)^2
                  + dudy(k,i,j) + dudy(k,i-1,j)  &
                  + dudy(k,i,j-jd) + dudy(k,i-1,j-jd)) ** 2 ! (avg du/dy)^2
            ! vertical shear of horizontal wind magnitude squared =
            ! (du/dz)^2 + (dv/dz)^2
            ! .0625=1/16
            straintens2(k,i,j) = .0625 * ((dudz(k,i,j) + dudz(k-1,i,j)  &
                  + dudz(k,i-1,j) + dudz(k-1,i-1,j)) ** 2  & ! (avg du/dz)^2
               + (dvdz(k,i,j) + dvdz(k-1,i,j)  &
                  + dvdz(k,i,j-jd) + dvdz(k-1,i,j-jd)) ** 2) ! (avg dv/dz)^2
         enddo
      enddo
   enddo
else
   ! code for idiffk=3,4 when dx~dz
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            dudx(k,i,j) = 2. * dudx(k,i,j) ! 2*du/dx
            dvdy(k,i,j) = 2. * dvdy(k,i,j) ! 2*dv/dy
            dwdz(k,i,j) = 2. * dwdz(k,i,j) ! 2*dw/dz
            dvdx(k,i,j) = dvdx(k,i,j) + dudy(k,i,j) ! dv/dx + du/dy
            dudy(k,i,j) = dvdx(k,i,j) ! dv/dx + du/dy
            dudz(k,i,j) = dudz(k,i,j) + dwdx(k,i,j) ! du/dz + dw/dx
            dvdz(k,i,j) = dvdz(k,i,j) + dwdy(k,i,j) ! dv/dz + dw/dy
            ! 2/3*( du/dx + dv/dy + dw/dz )
            straintens2(k,i,j) = 0.333333  &
               * (dudx(k,i,j) + dvdy(k,i,j) + dwdz(k,i,j)) 
         enddo
      enddo

      ! east and west boundaries
      do k = 2,m1-1
         dudx(k,iz1,j) = 2. * dudx(k,iz1,j) ! 2*du/dx along east edge
         dvdx(k,ia_1,j) = dvdx(k,ia_1,j) + dudy(k,ia_1,j) ! dv/dx + du/dy at west edge
         dudy(k,ia_1,j) = dvdx(k,ia_1,j) ! dv/dx + du/dy at west edge
         dudz(k,ia_1,j) = dudz(k,ia_1,j) + dwdx(k,ia_1,j) ! du/dz + dw/dx at west edge
      enddo
   enddo ! j loop

   ! south and north boundaries
   do i = ia_1,iz
      do k = 2,m1-1
         dvdy(k,i,jz1) = 2. * dvdy(k,i,jz1) ! 2*dv/dy at north edge
         dvdx(k,i,ja_1) = dvdx(k,i,ja_1) + dudy(k,i,ja_1) ! dv/dx + du/dy at south edge
         dudy(k,i,ja_1) = dvdx(k,i,ja_1) ! dv/dx + du/dy at south edge
         dvdz(k,i,ja_1) = dvdz(k,i,ja_1) + dwdy(k,i,ja_1) ! dv/dz + dw/dy at south edge
      enddo
   enddo

   ! 3D strain rate tensor magnitude squared calculation - 
   ! is first term magnitude right? Needs checking
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            ! 1/2*[ (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2 ] + [(dv/dx)+(du/dy)]^2 
            !   + [(du/dz)+(dw/dx)]^2 + [(dv/dz)+(dw/dy)]^2
            straintens1(k,i,j) = .5 * (  & ! 1/2 * [ (du/dx)^2 + (dv/dy)^2 + (dw/dz)^2 ]
                  (dudx(k,i,j) - straintens2(k,i,j)) ** 2  &
                  + (dvdy(k,i,j) - straintens2(k,i,j)) ** 2  &
                  + ( dwdz(k,i,j) - straintens2(k,i,j)) ** 2)  &
               + .0625 * ((dvdx(k,i,j) + dvdx(k,i-1,j)  & ! [avg(dv/dx) + avg(du/dy)]^2
                     + dvdx(k,i,j-jd) + dvdx(k,i-1,j-jd)) ** 2  &
                  + (dudz(k,i,j) + dudz(k,i-1,j)  & ! [avg(du/dz) + avg(dw/dx)]^2
                     + dudz(k-1,i,j) + dudz(k-1,i-1,j)) ** 2  &
                  + (dvdz(k,i,j) + dvdz(k-1,i,j)  & ! [avg(dv/dz) + avg(dw/dy)]^2
                     + dvdz(k,i,j-jd) + dvdz(k-1,i,j-jd)) ** 2)
            straintens2(k,i,j) = straintens1(k,i,j)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE strain

!##############################################################################
Subroutine bruvais (m1,m2,m3,ia,iz,ja,jz,theta,rtp,rv,rcp,pp,pi0,en2,rtgt)

use mem_scratch
use micphys
use mem_grid
use rconstants
use mem_turb, only:fracsat

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: theta,rtp,rcp,rv,pp,pi0,en2
real, dimension(m2,m3) :: rtgt
integer :: i,j,k,iweten,ki
real :: c1,c2,c3,ci1,ci2,ci3,rvlsi,rvii
real, dimension(nzpmax) :: pi,temp,prt,rvls,rc

! Calculate brunt-vaisala frequency squared (en2).
! From Hill(1974) / Durran and Klemp(1982) (DK82) moist BV modification.
! ldgrant(2016): added comments

!IWETEN:
!Flag for using saturated-BV (SBV) when cloud water exists (RC>0) and
!water vapor mixing ratio (RV) is greater than a certain fraction of
!saturated water vapor mixing ratio (RVLS). Otherwise, use moist-BV (MBV).
!We use the RV, RVLS condition to prevent applying SBV in cloudy but
!very sub-saturated conditions where evaporation is strong.
iweten = 1

!Some constants used in the BV calculations
c1 = alvl / rgas ! Llv / Rd [K]
c2 = ep * alvl ** 2 / (cp * rgas) ! Rd/Rv * Llv^2 / (cp*Rd) [K^2]
c3 = alvl / cp ! Llv / cp [K]
ci1 = alvi / rgas ! Liv / Rd [K]
ci2 = ep * alvi ** 2 / (cp * rgas) ! Rd/Rv * Liv^2 / (cp*Rd) [K^2]
ci3 = alvi / cp ! Liv / cp [K]
ki=m1+1

!     calculate potential temperature profile

do j = ja,jz
   do i = ia,iz

      do k = 1,m1
         vctr11(k) = theta(k,i,j)
         vctr12(k) = theta(k,i,j)
         vctr32(k) = 0.
      enddo
      if (level .ge. 1) then
         do k = 1,m1
            vctr12(k) = vctr11(k) * (1. + .61 * rv(k,i,j)) ! theta-v
            vctr32(k) = (rtp(k,i,j) - rv(k,i,j)) ! condensate mixing ratio
         enddo
      endif

!     check for saturation if level is 2 or greater.

      if (level .ge. 2 .and. iweten .eq. 1) then
         do k = 1,m1
            pi(k) = (pp(k,i,j) + pi0(k,i,j)) / cp ! non-dim. exner
            temp(k) = theta(k,i,j) * pi(k) ! temp [K]
            prt(k) = p00 * pi(k) ** cpor ! pressure [Pa]
         enddo
         CALL mrsl (m1,prt(1),temp(1),rvls(1)) ! saturation m.r. wrt liq at all levels
         do k = 1,m1
            vctr2(k) = c1 ! Llv/Rd [K]
            vctr3(k) = c2 ! Rd/Rv * Llv^2 / (cp*Rd) [K^2]
            vctr4(k) = c3 ! Llv/cp [K]
         enddo
         ki = m1 + 1

!     if any ice phase microphysics are activated ....

         if ((ipris .ge. 1 .or. isnow .ge. 1 .or.  &
             igraup .ge. 1 .or. iaggr .ge. 1 .or.  &
             ihail .ge. 1) .and. level .eq. 3) then

!            find level of -20 c.  assume ice saturation above this
!              level.

            do k = 1,m1
               if (temp(k) .le. 253.16) then
                  ki = k ! first level where T<=-20C
                  go to 10
               endif
            enddo
            ki = m1 + 1 ! if no T's below -20C, ki=model top level
10               continue

            ! rvls is sat m.r. wrt ice at all levels above ki
            CALL mrsi (m1-ki+1,prt(ki),temp(ki),rvls(ki))  
            ! rvlsi is sat m.r. wrt ice at level ki-1
            if(ki > 1) CALL mrsi (1,prt(ki-1),temp(ki-1),rvlsi)
            do k = ki,m1
               vctr2(k) = ci1 ! Liv/Rd [K]
               vctr3(k) = ci2 ! Rd/Rv * Liv^2 / (cp*Rd) [K^2]
               vctr4(k) = ci3 ! Liv/cp [K]
            enddo
         endif

         if (level .eq. 3) then
            do k = 1,m1
               rc(k) = rcp(k,i,j) ! cloud mixing ratio
            enddo
         else
            do k = 1,m1 ! "saturation adjustment"
               rc(k) = max(rv(k,i,j) / rvls(k) - .999,0.)
            enddo
         endif

      endif

      do k = 2,m1-1
         ! g/(2*dz)
         vctr1(k) = g / ((zt(k+1) - zt(k-1)) * rtgt(i,j))
      enddo
      if (level .ge. 2 .and. iweten .eq. 1) then
         do k = 2,m1-1
            !Do saturated BV equation if cloud water is present and water vapor
            !mixing ratio is greater than user chosen fraction of saturation
            !water vapor mixing ratio. 
            if ( rc(k) .gt. 0.0 .and. rv(k,i,j) >= fracsat*rvls(k) ) then
               rvii = rvls(k-1) ! sat m.r. at level below
               if (k .eq. ki) rvii = rvlsi ! at level where T<=-20C
               ! g/(2*dz) * [ (1+Llv/Rd*rvs/T) / (1+eps*Llv^2/cp/Rd*rvs/T^2)
               !             * (dtheta/theta+Llv/cp/T*drvs) - drt ]
               ! need 2*dz in denominator since dtheta and drvs are calculated 
               ! with a centered difference.
               ! Durran and Klemp (1982) Eq.36.
               en2(k,i,j) = vctr1(k) * (  &
                  (1. + vctr2(k) * rvls(k) / temp(k))  &
                  / (1. + vctr3(k) * rvls(k) / temp(k) ** 2)  &
                  * ((vctr11(k+1) - vctr11(k-1)) / vctr11(k)  &
                     + vctr4(k) / temp(k) * (rvls(k+1) - rvii))  &
                  - (rtp(k+1,i,j) - rtp(k-1,i,j)))
            else ! no cloud water
               ! g/(2*dz) * [ dthetav/thetav - drt ]
               en2(k,i,j) = vctr1(k)*((vctr12(k+1)-vctr12(k-1))  &
                  / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
            endif
         enddo
      else
         do k = 2,m1-1
            en2(k,i,j) = vctr1(k) * ((vctr12(k+1)-vctr12(k-1))  &
               / vctr12(k) - (vctr32(k+1) - vctr32(k-1)))
         enddo
      endif
      ! lower and upper boundary conditions
      en2(1,i,j) = en2(2,i,j)
      en2(nzp,i,j)=en2(nz,i,j)

   enddo
enddo

return
END SUBROUTINE bruvais

!##############################################################################
Subroutine mxdefm (m1,m2,m3,ia,iz,ja,jz  &
   ,strain_Khz,straintens,en2,richnum,Kmz,Kmx,dn0,rtgt,dxt)

!     +-------------------------------------------------------------+
!     \   this routine calculates the mixing coefficients with a    \
!     \     smagorinsky-type deformational based k with an optional \
!     \     unstable brunt-vaisala enhancement and an optional      \
!     \     richardson number modification.                         \
!     +-------------------------------------------------------------+

! strain_Khz (input) = (horiz. or 3D) strain rate tensor magnitude sqd from subroutine strain
! strain_Khz (output) = Khz = vertical diffsn coef for scalars for idiffk = 2, 3
! straintens = (du/dz)^2 + (dv/dz)^2 if idiffk = 1 or 2
!   or = 3D strain rate tensor squared if idiffk = 3, from subroutine strain
! en2 = N^2 = moist brunt vaisala frequency from subroutine bruvais
! richnum = Richardson #, calculated here
! Kmz = vertical diffusion coef for momentum vars 
!   set here for idiffk = 2, 3
! Kmx = horiz diffusion coef for momentum vars, 
!   set for all idiffk = 1,2,3

use mem_scratch
use mem_grid
use mem_turb
use rconstants
use micphys, only:idiffperts

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: strain_Khz,straintens,en2,richnum,Kmz,Kmx,dn0
real, dimension(m2,m3) :: rtgt,dxt
integer :: i,j,k,irich,ienfl
real :: cxsqd,enfl,rchmax,rtgtsqd,deltaxsqd,cxsqd_dxsqd,akm,rmin,rmax
real :: Lilly_RiN_mod

irich = 1 ! flag for maximum Richardson number modification
ienfl = 1 ! flag for Hill(1974)/Durran and Klemp(1982) moist BV modification

if(idiffperts > 0) irich = 0 
!Note: setting irich to 0 does not turn off the Ri# mod by itself, 
!need to also set rmax to 0 below.

! csx is a namelist parameter, dimensionless constant (Smagorinsky constant)
cxsqd = csx(ngrid) * csx(ngrid) ! csx sqd
if (idiffk(ngrid) .eq. 2 .or. idiffk(ngrid) .eq. 3) then
   ! Smagorinsky-vertical scheme necessary calculations
   rmin = -100.
   rmax = 1. / zkhkm(ngrid)
   if(irich==0) rmax=0.0
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            ! calculate Richardson #
            richnum(k,i,j) = max(min(en2(k,i,j)  &
               / max(straintens(k,i,j),1.e-15),rmax),rmin)
         enddo
      enddo
   enddo
   enfl = float(ienfl) ! =0 or 1
   rchmax = 1.0 + 9.0 * float(irich) ! =1 if irich=0, =10 if irich=1
   do k = 2, m1
      vctr1(k) = csz(ngrid) * (zm(k) - zm(k-1)) ! csz*dz
      vctr2(k) = vctr1(k) * vctr1(k) ! (csz*dz)^2
   enddo
endif

if (idiffk(ngrid) .eq. 1) then 
   ! Smagorinsky-type horizontal, but MY vertical (calculated elsewhere)
   do j = ja,jz
      do i = ia,iz
         deltaxsqd = 1.0 / (dxt(i,j) * dxt(i,j)) ! dx^2
         cxsqd_dxsqd = cxsqd * deltaxsqd ! (csx*dx)^2
         !0.075 * akmin * dx^4/3
         akm = akmin(ngrid) * 0.075 * deltaxsqd ** (0.666667)
         do k = 2,m1-1
            ! horizontal diffusion coef for momentum vars [kg/m3 * m2/s]
            Kmx(k,i,j) = dn0(k,i,j)  &
             * max(akm,cxsqd_dxsqd*sqrt(strain_Khz(k,i,j))) !strain rate tensor mag
         enddo
      enddo
   enddo
elseif (idiffk(ngrid) .eq. 2) then 
   ! Smagorinsky anisotropic scheme
   do j = ja,jz
      do i = ia,iz
         rtgtsqd = rtgt(i,j) * rtgt(i,j)
         deltaxsqd = 1.0 / (dxt(i,j) * dxt(i,j)) ! dx^2
         cxsqd_dxsqd = cxsqd * deltaxsqd ! (csx*dx)^2
         !0.075 * akmin * dx^4/3
         akm = akmin(ngrid) * 0.075 * deltaxsqd ** (0.666667)
         do k = 2,m1-1

         !Kmz: vertical diffusion coef. for momentum vars [kg/m3 * m2/s]
         !Kmz = rho0 * (csz*dz)^2 * {sqrt[(du/dz)^2+(dv/dz)^2] 
         !      + FlagN2*sqrt(max(0,-N^2))}
         Kmz(k,i,j) = dn0(k,i,j) * rtgtsqd * vctr2(k) &
           * ( sqrt(straintens(k,i,j))  &  !strain tensor
               + enfl * sqrt(max(0.,-en2(k,i,j))) ) !Hill-Durran stability mod (N^2) 
         !Lilly Richardson # modification: sqrt(max(0,1-[Kh/Km]*Ri))
         Lilly_RiN_mod = min(rchmax,sqrt(max(0.,(1.-zkhkm(ngrid)*richnum(k,i,j))))) 
         Kmz(k,i,j) = Kmz(k,i,j) * Lilly_RiN_mod
         !Kmx: horizontal diffusion coef. for momentum vars
         Kmx(k,i,j) = dn0(k,i,j)  &
          * max(akm,cxsqd_dxsqd*sqrt(strain_Khz(k,i,j))) !strain rate tensor mag
         !strain_Khz: here this variable becomes Khz:
         !vertical diffusion coef for scalars [kg/m3 * m2/s]
         strain_Khz(k,i,j) = Kmz(k,i,j) * zkhkm(ngrid)
         enddo
      enddo
   enddo

elseif (idiffk(ngrid) .eq. 3) then 
   ! isotropic Smagorinksy scheme
   do j = ja,jz
      do i = ia,iz
         rtgtsqd = rtgt(i,j) * rtgt(i,j)
         do k = 2,m1-1
            ! Kmz = vertical diffusion coef for momentum vars [kg/m3 * m2/s]
            ! ldgrant(2016): bug fix: vctr2=(csz*dz)^2 but the coefficient
            !  here should be (csx*dx)*(csz*dz). Changed vctr2 to 
            !  vctr1*csx*/dxt where vctr1=csz*dz and 1/dxt=deltax
            Kmz(k,i,j) = dn0(k,i,j) * rtgtsqd * vctr1(k)*csx(ngrid)/dxt(i,j)  &

               * ( sqrt(strain_Khz(k,i,j))  & ! strain rate tensor magnitude
                  + enfl * sqrt(max(0.,-en2(k,i,j))) ) ! Hill stability mod
            ! ldgrant(2016): separate out Lilly Ri# modification - comment out
            !  next two lines here to turn it off for idiffk=3
            Lilly_RiN_mod = min(rchmax,sqrt(max(0.,(1.-zkhkm(ngrid)*richnum(k,i,j))))) 
            Kmz(k,i,j) = Kmz(k,i,j) * Lilly_RiN_mod
            ! Kmx = horizontal diffusion coef. for momentum vars, set equal to Kmz
            Kmx(k,i,j) = Kmz(k,i,j)
            ! Khz = vertical diffusion coef for scalars
            strain_Khz(k,i,j) = Kmz(k,i,j) * zkhkm(ngrid) ! Khz
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE mxdefm

!##############################################################################
Subroutine klbnd (m1,m2,m3,ibcon,jd,akay,dn0)

implicit none

integer :: m1,m2,m3,ibcon,jd
real, dimension(m1,m2,m3) :: akay,dn0
integer ::i,j,k

!     boundary conditions on a mixing coefficient

do j = 1,m3
   do i = 1,m2
      akay(1,i,j) = akay(2,i,j) * dn0(1,i,j) / dn0(2,i,j)
      akay(m1,i,j) = akay(m1-1,i,j) * dn0(m1,i,j)  &
         / dn0(m1-1,i,j)
   enddo
enddo

if (iand(ibcon,1) .ne. 0) then
   do j = 1,m3
      do k = 1,m1
         akay(k,1,j) = akay(k,2,j)
      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = 1,m3
      do k = 1,m1
         akay(k,m2,j) = akay(k,m2-1,j)
      enddo
   enddo
endif

if (jd .eq. 1) then
   if (iand(ibcon,4) .ne. 0) then
      do i = 1,m2
         do k = 1,m1
            akay(k,i,1) = akay(k,i,2)
         enddo
      enddo
   endif

   if (iand(ibcon,8) .ne. 0) then
      do i = 1,m2
         do k = 1,m1
            akay(k,i,m3) = akay(k,i,m3-1)
         enddo
      enddo
   endif
endif

return
END SUBROUTINE klbnd
