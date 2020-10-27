!##############################################################################
Subroutine diffvel (m1,m2,m3,ia,iz,ja,jz,jd  &
   ,iz1,jz1,izu,jzv,idiffkk  &
   ,up,vp,wp,ut,vt,wt,vt3da,vt3db,vt3dc  &
   ,vt3dd,vt3de,vt3df,vt3dg,vt3dj,vt3dk  &
   ,vt3dl,vt3dm,vt3dn,vt3do,rtgu,rtgv,rtgt  &
   ,sflux_u,sflux_v,sflux_w,dn0,dn0u,dn0v,scr1,scr2)

use mem_grid
use mem_scratch
use mem_turb

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd,iz1,jz1,izu,jzv,idiffkk,i,j,k

real :: akn,ako,akp,cross,c1,c2,dtlvi
real, dimension(m1,m2,m3) :: up,vp,wp,ut,vt,wt,vt3da,vt3db,vt3dc,vt3dd  &
   ,vt3de,vt3df,vt3dg,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn,vt3do,scr1,scr2  &
   ,dn0,dn0u,dn0v
real, dimension(m2,m3) :: sflux_u,sflux_v,sflux_w,rtgt,rtgu,rtgv

!          compute fluxes with k's and previous strains

if (ihorgrad .eq. 1) then

   do j = 1,jz1
      do i = 1,iz1
         do k = 1,m1-1
            vt3da(k,i,j) = -vt3da(k,i,j) * scr2(k,i,j)
            vt3dc(k,i,j) = -vt3dc(k,i,j) * scr2(k,i,j)
         enddo
      enddo
   enddo

   do j = 1,jz
      do i = 1,iz
         do k = 1,m1-1
            akn = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )  &
                       + scr2(k,i  ,j+jd) + scr2(k+1,i  ,j+jd))
            ako = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )  &
                       + scr2(k,i+1,j   ) + scr2(k+1,i+1,j   ))
            akp = .25 * (scr2(k,i  ,j   ) + scr2(k  ,i+1,j   )  &
                       + scr2(k,i  ,j+jd) + scr2(k  ,i+1,j+jd))

            vt3db(k,i,j) = -vt3db(k,i,j) * akp
            vt3dn(k,i,j) = -vt3dn(k,i,j) * akp
            vt3dd(k,i,j) = -vt3dd(k,i,j) * ako
            vt3df(k,i,j) = -vt3df(k,i,j) * ako
            vt3de(k,i,j) = -vt3de(k,i,j) * akn
            vt3dg(k,i,j) = -vt3dg(k,i,j) * akn
         enddo
      enddo
   enddo

elseif(ihorgrad.eq.2)then

   do j = 1,jz
      do i = 1,iz
         do k = 1,m1-1
            akn = .25 * (scr2(k,i  ,j ) + scr2(k+1,i  ,j )  &
                       + scr2(k,i  ,j+jd) + scr2(k+1,i  ,j+jd))
            ako = .25 * (scr2(k,i  ,j ) + scr2(k+1,i  ,j )  &
                       + scr2(k,i+1,j ) + scr2(k+1,i+1,j ))

            vt3df(k,i,j) = -vt3df(k,i,j) * ako
            vt3dg(k,i,j) = -vt3dg(k,i,j) * akn
         enddo
      enddo
   enddo

endif

cross = 0.
if (idiffkk .ge. 3) cross = 1.

!------------------------ vertical u-component diffusion ----------------

do j = ja,jz
   do i = ia,izu
      do k = 1,m1-1
         vctr1(k) = cross * vt3df(k,i,j)
         vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
                  + scr1(k,i+1,j) + scr1(k+1,i+1,j))
      enddo

      if (nstbot .eq. 1) then
         vctr1(1) = .5 * (sflux_u(i,j) + sflux_u(i+1,j))
      else
         vctr1(1) = vctr1(1) + up(1,i,j)  &
            * .25 * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i+1,j)  &
            + scr1(2,i+1,j)) * dzm(1) / rtgu(i,j)
      endif

      if (nsttop .eq. 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = vctr1(m1-1) - up(m1-1,i,j)  &
            * .25 * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
            + scr1(m1,i+1,j) + scr1(m1-1,i+1,j)) * dzm(m1-1)  &
            / rtgu(i,j)
      endif

      c1 = dtlv / rtgu(i,j)
      c2 = .25 * dtlv / (rtgu(i,j) * rtgu(i,j))
      do k = 2,m1-1
         vt3dl(k,i,j) = up(k,i,j) * dn0u(k,i,j)  &
                      + c1 * dzt(k) * (vctr1(k-1) - vctr1(k))
         vt3dj(k,i,j) = - c2 * dzt(k) * vctr2(k-1)
         vt3dk(k,i,j) = - c2 * dzt(k) * vctr2(k)
         vt3do(k,i,j) = dn0u(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
      enddo
      vt3dj(2,i,j) = 0.
      if (nstbot .eq. 1) vt3do(2,i,j) = dn0u(2,i,j) - vt3dk(2,i,j)
      vt3dk(m1-1,i,j) = 0.
      if (nsttop .eq. 1) vt3do(m1-1,i,j) = dn0u(m1-1,i,j)  &
         - vt3dj(m1-1,i,j)

   enddo
enddo

CALL tridiff1 (m1,m2,m3,ia,izu,ja,jz,m1-1  &
   ,vt3dj,vt3do,vt3dk,vt3dl,vt3do,vt3dk)

!---------------------horizontal u-component diffusion ----------------

if(ihorgrad.eq.1)then
   CALL divcart (m1,m2,m3,ia,izu,ja,jz,vt3da,vt3dj,'XDIR','TPNT')
   CALL divcart (m1,m2,m3,ia,izu,ja,jz,vt3dn,vt3dk,'YDIR','PPNT')

elseif(ihorgrad.eq.2)then
!        average the dn0*hkh to the velocity points
   CALL avgvel (m1,m2,m3,ia,izu,ja,jz,'xdir',jd,vt3da,scr2)
   CALL truhor (m1,m2,m3,ia,izu,ja,jz  &
              ,up,vt3dj,'xdir','dxt',grid_g(ngrid)%dxt(1,1)  &
              ,grid_g(ngrid)%topu(1,1),grid_g(ngrid)%rtgu(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,vt3da,'NONE')


   CALL avgvel (m1,m2,m3,ia,izu,ja,jz,'ydir',jd,vt3da,scr2)
   CALL truhor (m1,m2,m3,ia,izu,ja,jz  &
              ,up,vt3dk,'ydir','dym',grid_g(ngrid)%dym(1,1)  &
              ,grid_g(ngrid)%topu(1,1),grid_g(ngrid)%rtgu(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,vt3da,'NONE')
endif

dtlvi = 1.0 / dtlv

if (jd .eq. 1) then
   do j = ja,jz
      do i = ia,izu
         do k = 2,m1-1
            ut(k,i,j) = ut(k,i,j)  &
               + dtlvi * (vt3do(k,i,j) - up(k,i,j))  &
               - (vt3dj(k,i,j) + vt3dk(k,i,j)) / dn0u(k,i,j)
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,izu
         do k = 2,m1-1
            ut(k,i,j) = ut(k,i,j)  &
               + dtlvi * (vt3do(k,i,j)-up(k,i,j))  &
               - vt3dj(k,i,j) / dn0u(k,i,j)
         enddo
      enddo
   enddo
endif

!------------------------ vertical v-component diffusion ----------------

if (jd .eq. 0) go to 99

do j = ja,jzv
   do i = ia,iz
      do k = 1,m1-1
         vctr1(k) = cross * vt3dg(k,i,j)
         vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
            + scr1(k,i,j+jd) + scr1(k+1,i,j+jd))
      enddo

      if (nstbot .eq. 1) then
         vctr1(1) = .5 * (sflux_v(i,j) + sflux_v(i,j+jd))
      else
         vctr1(1) = vctr1(1) + vp(1,i,j)  &
            * .25 * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i,j+jd)  &
            + scr1(2,i,j+jd)) * dzm(1) / rtgv(i,j)
      endif
      if (nsttop .eq. 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = vctr1(m1-1) - vp(m1-1,i,j)  &
            * .25 * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
            + scr1(m1,i,j+jd) + scr1(m1-1,i,j+jd)) * dzm(m1-1)  &
            / rtgv(i,j)
      endif

      c1 = dtlv / rtgv(i,j)
      c2 = .25 * dtlv / (rtgv(i,j) * rtgv(i,j))
      do k = 2,m1-1
         vt3dm(k,i,j) = vp(k,i,j) * dn0v(k,i,j)  &
                      + c1 * dzt(k) * (vctr1(k-1) - vctr1(k))
         vt3dj(k,i,j) = -c2 * dzt(k) * vctr2(k-1)
         vt3dk(k,i,j) = -c2 * dzt(k) * vctr2(k)
         vt3do(k,i,j) = dn0v(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
      enddo
      vt3dj(2,i,j) = 0.
      if (nstbot .eq. 1) vt3do(2,i,j) = dn0v(2,i,j) - vt3dk(2,i,j)

      vt3dk(m1-1,i,j) = 0.
      if (nsttop .eq. 1)  &
           vt3do(m1-1,i,j) = dn0v(m1-1,i,j) - vt3dj(m1-1,i,j)
   enddo
enddo

CALL tridiff1 (m1,m2,m3,ia,iz,ja,jzv,m1-1  &
   ,vt3dj,vt3do,vt3dk,vt3dm,vt3do,vt3dk)

!--------------------- horizontal v-component diffusion ----------------

if (ihorgrad .eq. 1) then
   CALL divcart (m1,m2,m3,ia,iz,ja,jzv,vt3dc,vt3dj,'YDIR','TPNT')
   CALL divcart (m1,m2,m3,ia,iz,ja,jzv,vt3db,vt3dk,'XDIR','PPNT')

elseif(ihorgrad.eq.2)then
!        average the dn0*hkh to the velocity points
   CALL avgvel (m1,m2,m3,ia,iz,ja,jzv,'ydir',jd,vt3da,scr2)
   CALL truhor (m1,m2,m3,ia,iz,ja,jzv  &
              ,vp,vt3dj,'ydir','dyt',grid_g(ngrid)%dyt(1,1)  &
              ,grid_g(ngrid)%topv(1,1),grid_g(ngrid)%rtgv(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,vt3da,'NONE')

   CALL avgvel (m1,m2,m3,ia,iz,ja,jzv,'xdir',jd,vt3da,scr2)
   CALL truhor (m1,m2,m3,ia,iz,ja,jzv  &
              ,vp,vt3dk,'xdir','dxm',grid_g(ngrid)%dxm(1,1)  &
              ,grid_g(ngrid)%topv(1,1),grid_g(ngrid)%rtgv(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,vt3da,'NONE')
endif

do j = ja,jzv
   do i = ia,iz
      do k = 2,m1-1
         vt(k,i,j) = vt(k,i,j) + dtlvi * (vt3do(k,i,j)-vp(k,i,j))  &
             - (vt3dj(k,i,j) + vt3dk(k,i,j)) / dn0v(k,i,j)
      enddo
   enddo
enddo

99   continue

!------------------------ vertical w-component diffusion ----------------

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         vctr1(k) = 0.
         vctr2(k) = dzt(k) * scr1(k,i,j)
      enddo

      if (nstbot .eq. 1) then
         vctr1(2) = sflux_w(i,j)
      else
         vctr1(2) = (1.0 + cross) * wp(1,i,j)  &
            * scr1(2,i,j) * dzt(2) / rtgt(i,j)
      endif
      if (nsttop .eq. 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = -(1.0 + cross) * wp(m1-1,i,j)  &
            * scr1(m1-1,i,j) * dzt(m1-1) / rtgt(i,j)
      endif

      c1 = dtlv / rtgt(i,j)
      c2 = 2.0 * dtlv / (rtgt(i,j) * rtgt(i,j))
      do k = 2,m1-2
         vt3dn(k,i,j) = wp(k,i,j)  &
            * (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
            + c1 * dzm(k) * (vctr1(k) - vctr1(k+1))
         vt3dj(k,i,j) = - c2 * dzm(k) * vctr2(k)
         vt3dk(k,i,j) = - c2 * dzm(k) * vctr2(k+1)
         vt3do(k,i,j) = (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
            - vt3dj(k,i,j) - vt3dk(k,i,j)
      enddo
      vt3dj(2,i,j) = 0.
      if (nstbot .eq. 1) vt3do(2,i,j) =  &
         (dn0(2,i,j) + dn0(3,i,j)) * .5 - vt3dk(2,i,j)
      vt3dk(m1-2,i,j) = 0.
      if (nsttop .eq. 1) vt3do(m1-2,i,j) =  &
           (dn0(m1-2,i,j) + dn0(m1-1,i,j)) * .5  &
           -vt3dj(m1-2,i,j)
   enddo
enddo

CALL tridiff1 (m1,m2,m3,ia,iz,ja,jz,m1-2  &
   ,vt3dj,vt3do,vt3dk,vt3dn,vt3do,vt3dk)

!--------------------- horizontal w-component diffusion ----------------

if (idiffkk .ge. 3) then
   if (ihorgrad .eq. 1) then
      CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3dd,vt3dj,'XDIR','OPNT')
      CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3de,vt3dk,'YDIR','NPNT')

   elseif(ihorgrad.eq.2)then
      CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
                 ,wp,vt3dj,'xdir','dxu',grid_g(ngrid)%dxu(1,1)  &
                 ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
                 ,zm,vctr1,vctr2,vctr3,vctr4,vctr5  &
                 ,vctr6,vctr7,jd,scr2,'NONE')

      CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
                 ,wp,vt3dk,'ydir','dyv',grid_g(ngrid)%dyv(1,1)  &
                 ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
                 ,zm,vctr1,vctr2,vctr3,vctr4,vctr5  &
                 ,vctr6,vctr7,jd,scr2,'NONE')
   endif

else
   if(ihorgrad.eq.1)then
      CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3df,vt3dj,'XDIR','OPNT')
      CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3dg,vt3dk,'YDIR','NPNT')

   elseif(ihorgrad.eq.2)then
      CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
                 ,wp,vt3dj,'xdir','dxu',grid_g(ngrid)%dxu(1,1)  &
                 ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
                 ,zm,vctr1,vctr2,vctr3,vctr4,vctr5  &
                 ,vctr6,vctr7,jd,scr2,'NONE')

      CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
                 ,wp,vt3dk,'ydir','dyv',grid_g(ngrid)%dyv(1,1)  &
                 ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
                 ,zm,vctr1,vctr2,vctr3,vctr4,vctr5  &
                 ,vctr6,vctr7,jd,scr2,'NONE')
   endif
endif

if (jd .eq. 1) then
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-2
            wt(k,i,j) = wt(k,i,j)  &
               + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
               - (vt3dj(k,i,j) + vt3dk(k,i,j))  &
               / ((dn0(k,i,j) + dn0(k+1,i,j)) * .5)
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-2
             wt(k,i,j) = wt(k,i,j)  &
               + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
               - vt3dj(k,i,j)  &
               / ((dn0(k,i,j) + dn0(k+1,i,j)) * .5)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE diffvel

!##############################################################################
Subroutine diffsclr (m1,m2,m3,ia,iz,ja,jz,jd  &
   ,n,ksf  & 
   ,scp,sct,fluxdivx,fluxdivy,vt3df,vt3dg  &
   ,coefz_km1,coefz_kp1,coefz_k,rhoscp,dn03i,vt3dl,vt3dm,icoefz_k_2,rtgt &
   ,sfcflx,dn0,th0,rvt0,th00,rvt00,vkkh,hkkh) 

! ldgrant(2016): added lots of comments and changed scratch variable names

! ksf = a flag to tell whether to compute vertical diffusion matrix coeffs
! scp = scalar variable (at past time level)
! sct = scalar tendency
! fluxdivx = flux divergence in x-direction [kg/m3 * [scp]/s]
! fluxdivy = flux divergence in y-direction [kg/m3 * [scp]/s]
! vt3df = 
! vt3dg = 
! for the next 3 vars: if rearrange the vertical part of scalar tendency eqn 
!  so LHS has scp(k) at time t+dt, then the following are coefficients: 
! coefz_km1 = coef on scp(k-1) [kg/m3]
! coefz_kp1 = coef on scp(k+1) [kg/m3]
! coefz_k = coef on scp(k) [kg/m3] 
! rhoscp = rho0*scp [kg/m3 * scalar units]
! dn03i = 1/rho0 [m3/kg]
! vt3dl = 
! vt3dm =
! icoefz_k_2 = 1/coefz_k(2) [m3/kg]
! sfcflx = sflux_r if scalar = RTP or sflux_t if scalar = THP or 0 otherwise
! vkkh = vkh at past time level [kg/m3 * m2/s]
! hkkh = hkh at past time level [kg/m3 * m2/s]

use mem_grid
use mem_scratch
use mem_turb
use var_tables
use mem_basic
use micphys, only:idiffperts
use node_mod, only:mi0,mj0

implicit none

integer :: n,ksf

integer :: m1,m2,m3,ia,iz,ja,jz,jd,i,j,k
integer, save :: ksf_save = 0
real, dimension(:,:,:), allocatable :: scp_diffuse

real :: dto2,dtlti
real, dimension(m1,m2,m3) :: scp,sct,fluxdivx,fluxdivy,rhoscp,vt3df,vt3dg &
   ,coefz_km1,coefz_kp1,coefz_k,vkkh,hkkh,dn0,dn03i,vt3dl,vt3dm,th0,rvt0  &
   ,th00,rvt00
real, dimension(m2,m3) :: sfcflx,rtgt,icoefz_k_2

allocate(scp_diffuse(m1,m2,m3))

! ldgrant(2016): subtract off base state from thp and rtp scalars to 
!  diffuse perturbation save it in variable scp_diffuse and use 
!  scp_diffuse instead of scp for the gradient, flux divergence, 
!  tridiagonal matrix solver etc. calculations.
!  note: scp_diffuse=scp for all scalars except THP and RTP
!
!  th0 is set to the base state theta-v. At time zero (when the base state 
!  is established), theta and theta-il (THP) will be equal since there is 
!  no condensate. Theta-v at time zero, however, will hold different values 
!  than theta (or theta-il) since there will be vapor at time zero. This is
!  why th0 is converted from theta-v to theta before subtracting from THP.
!
! For idiffperts=2 we compute the domain average state for the "new" base state.
! This means th00 = theta/thp and rvt00 = rtp
! smsaleeb(2017): th00 is set to the the domain average THP (theta-il), simply
!  subtract it from scp (when scp is set to THP).
!
! For idiffperts=3 we use varfiles to get interpolated 3D state from
! variables THETA and RV, so th00=theta and rvt00=rv from varfiles.

!IDIFFPERTS
!0=normal scalar diffusion
!1=diffuse perturbations THP,RTP from 3D base state (THETA-V,RV)
!2=diffuse perturbations THP,RTP from 3D current domain-mean state
!3=diffuse perturbations THP,RTP from 3D current varfile state (THETA,RV)

if( scalar_tab(n,ngrid)%name == 'THP' .and. idiffperts >= 1) then
 do j = 1,m3
  do i = 1,m2
   do k = 1,m1
    if(idiffperts==1)scp_diffuse(k,i,j) = scp(k,i,j)-(th0(k,i,j) &
                                         /(1.+0.61*rvt0(k,i,j)))
    if(idiffperts==2)scp_diffuse(k,i,j) = scp(k,i,j)-th00(k,i,j)
    if(idiffperts==3)scp_diffuse(k,i,j) = scp(k,i,j)-th00(k,i,j)
   enddo
  enddo
 enddo
elseif( scalar_tab(n,ngrid)%name == 'RTP' .and. idiffperts >= 1) then
 do j = 1,m3
  do i = 1,m2
   do k = 1,m1
    if(idiffperts==1)scp_diffuse(k,i,j) = scp(k,i,j)-rvt0(k,i,j)
    if(idiffperts==2)scp_diffuse(k,i,j) = scp(k,i,j)-rvt00(k,i,j)
    if(idiffperts==3)scp_diffuse(k,i,j) = scp(k,i,j)-rvt00(k,i,j)
   enddo
  enddo
 enddo
else ! all other scalars: diffuse full scalar field
 do j = 1,m3
  do i = 1,m2
   do k = 1,m1
    scp_diffuse(k,i,j) = scp(k,i,j)
   enddo
  enddo
 enddo
endif

!          compute vertical diffusion matrix coefficients for scalars
if (n == 1 .or. ksf /= ksf_save) then
   ksf_save = ksf
   do j = ja,jz
      do i = ia,iz
         do k = 1,m1-1
            ! 2*avg(vkkh)/dzm (m levs) [kg/m3*m/s]
            vctr1(k) = dzm(k) * (vkkh(k,i,j) + vkkh(k+1,i,j))
            ! Note - this doesn't account for stretched grid!
         enddo
         dto2 = .5 * dtlt / (rtgt(i,j) * rtgt(i,j)) ! dt/2
         do k = 2,m1-1
            ! -[avg(vkkh)/dzm from m level k-1] * [1/dz from t level k] * dt [kg/m3]
            coefz_km1(k,i,j) = -dto2 * dzt(k) * vctr1(k-1) 
            ! -[avg(vkkh)/dzm from m level k] * [1/dz from t level k] * dt [kg/m3]
            coefz_kp1(k,i,j) = -dto2 * dzt(k) * vctr1(k) 
            ! rho0 + avg(vkkh(k-1))/dzm(k-1)*dt/dzt(k) + avg(vkkh(k))/dzm(k)*dt/dzt(k)
            coefz_k(k,i,j) = dn0(k,i,j) - coefz_km1(k,i,j) - coefz_kp1(k,i,j) 
            dn03i(k,i,j) = 1. / dn0(k,i,j) ! [m3/kg] (t levs)
         enddo
         ! z=2 conditions
         coefz_km1(2,i,j) = 0.
         if (nstbot .eq. 1)  & ! indicates if this grid goes to the sfc
            ! exclude coefz_km1 from just below zt=2
            coefz_k(2,i,j) = dn0(2,i,j) - coefz_kp1(2,i,j)
         ! z=m1-1 conditions
         coefz_kp1(m1-1,i,j) = 0.
         if (nsttop .eq. 1)  & ! indicates if this grid goes to the model top
            ! exclude coefz_kp1 from just above zt=m1-1
            coefz_k(m1-1,i,j) = dn0(m1-1,i,j) - coefz_km1(m1-1,i,j)

         icoefz_k_2(i,j) = 1. / coefz_k(2,i,j) ! 1/coefz_k at zt=2
         ! coefz_kp1/coefz_k at zt=2 [-]
         vt3dl(2,i,j) = coefz_kp1(2,i,j) * icoefz_k_2(i,j)
         ! tri diagonal matrix algorithm coefficients
         do k = 3,m1-1
            ! [m3/kg ~ 1/rho]
            vt3dm(k,i,j) = 1. / (coefz_k(k,i,j) - coefz_km1(k,i,j) * vt3dl(k-1,i,j))
            vt3dl(k,i,j) = coefz_kp1(k,i,j) * vt3dm(k,i,j) ! [-]
         enddo
      enddo
   enddo

endif

!     compute 2 horizontal scalar gradients needed for dscp/dt

if (ihorgrad .eq. 1 .and. scalar_tab(n,ngrid)%name /= 'THP') then
   ! only for ihorgrad = 1 for all scalars except thp
   ! vt3df=d(scp)/dx(?)
   CALL grad (m1,m2,m3,1,iz,ja,jz,scp_diffuse,vt3df,'XDIR','TPNT')
   ! vt3dg=d(scp)/dy(?)
   CALL grad (m1,m2,m3,ia,iz,1,jz,scp_diffuse,vt3dg,'YDIR','TPNT')

   do j = ja,jz
      do i = 1,iz
         do k = 1,m1-1
            ! -d(scp)/dx*avg(hkkh) on u points [scp*kg/m2/s] = scalar mass flux
            vt3df(k,i,j) = -vt3df(k,i,j)  &
               * .5 * (hkkh(k,i,j) + hkkh(k,i+1,j))
         enddo
      enddo
   enddo

   do j = 1,jz
      do i = ia,iz
         do k = 1,m1-1
            ! -d(scp)/dy*avg(hkkh) on v points
            vt3dg(k,i,j) = -vt3dg(k,i,j)  &
               * .5 * (hkkh(k,i,j) + hkkh(k,i,j+jd))
         enddo
      enddo
   enddo
endif

!         horizontal flux divergence for scalars

if (ihorgrad .eq. 1 .and. scalar_tab(n,ngrid)%name /= 'THP') then
   ! This is called for ihorgrad=1 for all scalars except thp
   ! fluxdivx = scalar flux divergence in x direction [scp*kg/m3/s]
   ! fluxdivy = scalar flux divergence in y direction
   CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3df,fluxdivx,'XDIR','UPNT')

   CALL divcart (m1,m2,m3,ia,iz,ja,jz,vt3dg,fluxdivy,'YDIR','VPNT')

elseif (ihorgrad .eq. 2 .or. scalar_tab(n,ngrid)%name == 'THP') then
   ! This is called for ihorgrad=2 for all scalars, or just for 
   !  thp if ihorgrad=1 - haven't gone through this subroutine
   CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
              ,scp_diffuse,fluxdivx,'xdir','dxu',grid_g(ngrid)%dxu(1,1)  &
              ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,hkkh,scalar_tab(n,ngrid)%name)

   CALL truhor (m1,m2,m3,ia,iz,ja,jz  &
              ,scp_diffuse,fluxdivy,'ydir','dyv',grid_g(ngrid)%dyv(1,1)  &
              ,grid_g(ngrid)%topt(1,1),grid_g(ngrid)%rtgt(1,1)  &
              ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
              ,vctr6,vctr7,jd,hkkh,scalar_tab(n,ngrid)%name)
endif

!         finish matrix coefficients

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         ! set rhoscp [kg/m3 * scalar units]
         rhoscp(k,i,j) = scp_diffuse(k,i,j) * dn0(k,i,j)
      enddo

      ! lower boundary condition on scalar*rho0 - include surface fluxes
      if (nstbot .eq. 1) then ! indicates if this grid goes to the sfc
         ! units: K*rho (for THP); kg/kg*rho (for RTP)
         rhoscp(2,i,j) = scp_diffuse(2,i,j) * dn0(2,i,j)  &
            + sfcflx(i,j) * dtlt * dzt(2) / rtgt(i,j)
      else ! ?
         rhoscp(2,i,j) = scp_diffuse(2,i,j) * dn0(2,i,j)  &
            + .5 * dtlt * (vkkh(1,i,j) + vkkh(2,i,j))  &
            * scp_diffuse(1,i,j) * dzm(2) * dzt(2) / rtgt(i,j) ** 2
      endif

      ! top boundary condition (m1-1) on scalar*rho0 - not sure what this is doing
      if (nsttop .eq. 0) then ! is this check right?
         rhoscp(m1-1,i,j) = scp_diffuse(m1-1,i,j) * dn0(m1-1,i,j)  &
            - .5 * dtlt * (vkkh(m1-1,i,j) + vkkh(m1,i,j))  &
            * scp_diffuse(m1,i,j) * dzm(m1-1) * dzt(m1) / rtgt(i,j) ** 2
      endif
   enddo
enddo

! vertical diffusion loop - tridiagonal matrix algorithm! but what is on the LHS - 
!  scp(t+dt)? It must be implicit in time?
! trapezoidal scheme (i.e. Crank-Nicolson method) (e.g. Durran text p. 137): 
! [scp(k,t+dt)-scp(k,t)]/dt = M/2*{d^2[scp(k,t+dt)]/dx^2+d^2[scp(k,t)]/dx^2}
! Durran p. 142
do j = ja,jz
   do i = ia,iz

      ! vt3df is reset here. First: initialize k=2
      ! sclr*rho * (1/coef with rho units)
      vt3df(2,i,j) = rhoscp(2,i,j) * icoefz_k_2(i,j)

      ! loop up from k=3 to m1-1 - tridiagonal solver with coefficients
      do k = 3,m1-1
         vt3df(k,i,j) = (rhoscp(k,i,j) - coefz_km1(k,i,j) * vt3df(k-1,i,j))  &
            * vt3dm(k,i,j)
      enddo

      ! loop down from k=m1-2 to 2 and back out solution which is scp at t+dt(?)
      ! vt3df(m1-1) is the solution at that level and the others follow from there
      do k = m1-2,2,-1
         ! This solves for the new variable, scp(k) at t+dt ?
         vt3df(k,i,j) = vt3df(k,i,j) - vt3dl(k,i,j) * vt3df(k+1,i,j)
      enddo

   enddo
enddo

! Now vt3df has scp at t+dt due to vertical diffusion, so to get the tendency,
! we must back it out with d(scp)/dt = [scp(t+dt)-scp(t)]/dt
! based on below, since sct calculation has scp/dt subtracted out in the horizontal 
! flux divergence part, then vt3df must be scp(t+dt) due to vertical diffusion

dtlti = 1.0 / dtlt
 
if (jd .eq. 1) then ! 3D
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            ! d(scp)/dt = d(scp)/dt - scp/dt 
            !   - 1/rho0*[d(rho0*Khx*dscp/dx)/dx + d(rho0*Khx*dscp/dy)/dy]
            ! fluxdivx: horiz flux div in x dir; fluxdivy: y dir
            sct(k,i,j) = sct(k,i,j) - scp_diffuse(k,i,j) * dtlti  & 
               - (fluxdivx(k,i,j) + fluxdivy(k,i,j)) * dn03i(k,i,j)
         enddo
      enddo
   enddo
else ! 2D
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            sct(k,i,j) = sct(k,i,j) - scp_diffuse(k,i,j) * dtlti  &
               - fluxdivx(k,i,j) * dn03i(k,i,j)
         enddo
      enddo
   enddo
endif

do j = ja,jz
   do i = ia,iz
      do k = 2,m1-1
         ! vt3df: new scalar val: back out tend
         sct(k,i,j) = sct(k,i,j) + vt3df(k,i,j) * dtlti
      enddo
   enddo
enddo

deallocate(scp_diffuse)

return
END SUBROUTINE diffsclr

!##############################################################################
Subroutine tridiff1 (m1,m2,m3,ia,iz,ja,jz,kz,cim1,ci,cip1,rhs,cj,cjp1)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,kz,i,j,k
real :: cji
real, dimension(m1,m2,m3) :: cim1,ci,cip1,rhs,cj,cjp1

do j = ja,jz
   do i = ia,iz

      cjp1(2,i,j) = cip1(2,i,j) / ci(2,i,j)
      rhs(2,i,j) = rhs(2,i,j) / ci(2,i,j)

      do k = 3,kz
         cj(k,i,j) = ci(k,i,j) - cim1(k,i,j) * cjp1(k-1,i,j)
         cji = 1. / cj(k,i,j)
         cjp1(k,i,j) = cip1(k,i,j) * cji
         rhs(k,i,j) = (rhs(k,i,j) - cim1(k,i,j) * rhs(k-1,i,j))  &
                    * cji
      enddo

      cj(kz,i,j) = rhs(kz,i,j)

      do k = kz-1,2,-1
         cj(k,i,j) = rhs(k,i,j) - cjp1(k,i,j) * cj(k+1,i,j)
      enddo

   enddo
enddo

return
END SUBROUTINE tridiff1

!##############################################################################
Subroutine truhor (m1,m2,m3,ia,iz,ja,jz  &
                 ,vc3da,vc3db,dir,gpnt,dxy,topo,rtg  &
                 ,z,vctr1,vctr2,vctr3,zintl,zintr  &
                 ,dxynu,dxytem,jd,dn0hkh,sclrname)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd,nz,jaa,jzz,je,jf,i,j,k
real :: delz1,delz3,distnu,distold,scale,aknu,akdiff,dthdz,lapse_rate_max
real, dimension(m1,m2,m3) :: vc3da,vc3db,dn0hkh
real, dimension(m2,m3) :: dxy,topo,rtg
real, dimension(m1) :: z,vctr1,vctr2,vctr3,zintl,zintr,dxynu,dxytem
character(len=*) :: dir,gpnt,sclrname

nz=m1

jaa=ja
jzz=jz
if(jd.eq.0) then
   jaa=1
   jzz=1
endif

  if(dir.eq.'xdir')then
    if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
      je=0
      jf=1
    else
      je=-1
      jf=0
    endif
    do j=jaa,jzz
      do i=ia,iz
        do k=1,m1
          vctr1(k)=topo(i-1,j)+z(k)*rtg(i-1,j)
          vctr2(k)=topo(i  ,j)+z(k)*rtg(i ,j)
          vctr3(k)=topo(i+1,j)+z(k)*rtg(i+1,j)
          dxytem(k)=vc3da(k,i-1,j)
          dxynu(k)=vc3da(k,i+1,j)
        enddo

        !extrapolate the gradient from 3--->2 to 2--->1.
        delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
        delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
        dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
        dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

        CALL htint (nz,dxytem,vctr1,nz,zintl,vctr2)
        CALL htint (nz,dxynu, vctr3,nz,zintr,vctr2)

        !check intersection of cartesian sfc. on left side of mtn.
        CALL topobnd (m1,m2,m3,i,j, 1,0,'l',vc3da,dxy(i+je,j)  &
                    ,vctr1,vctr2,vctr3,zintr,dxytem)
        !check intersection of cartesian sfc. on right side of mtn.
        CALL topobnd (m1,m2,m3,i,j,-1,0,'r',vc3da,dxy(i+jf,j)  &
                    ,vctr1,vctr2,vctr3,zintl,dxynu)

        !Saleeby(2009):Limit lapse rate to prevent runaway cooling
        if (sclrname == 'THP') then
          lapse_rate_max = 0.010
          if (zintr(2) < vc3da(2,i+je,j)) then
             dthdz = (vc3da(2,i+je,j)-zintr(2))/(vctr3(2)-vctr3(1))
             if (dthdz > lapse_rate_max) dthdz = lapse_rate_max
             zintr(2) = vc3da(2,i+je,j) - dthdz*(vctr3(2)-vctr3(1))
          endif
          if (zintl(2) < vc3da(2,i+jf,j)) then
             dthdz = (vc3da(2,i+jf,j)-zintl(2))/(vctr1(2)-vctr1(1))
             if (dthdz > lapse_rate_max) dthdz = lapse_rate_max
             zintl(2) = vc3da(2,i+jf,j) - dthdz*(vctr1(2)-vctr1(1))
          endif
          !if(ngrid==3.and.i==50.and.j==50) &
          !print*,'Lapse1:',sclrname,ngrid,zintr(2),zintl(2)
        endif 

        do k=2,m1
          distnu=1./dxytem(k)+1./dxynu(k)
          distold=1./dxy(i+je,j)+1./dxy(i+jf,j)
          scale=distnu**2/distold**2
          if(scale.gt.1)stop ' scale incorrect-xdir--'
          !Need to multiply hkh by 2 since going over 2dx
          aknu=dn0hkh(k,i,j)*scale*2.
          akdiff=-aknu/distnu
          vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                              -(vc3da(k,i,j)-zintl(k))*dxynu(k) )
        enddo

        vc3db(1,i,j)=vc3db(2,i,j)
        if(gpnt.eq.'wpnt')then
          vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
        else
          vc3db(m1,i,j)=vc3db(m1-1,i,j)
        endif
      enddo
    enddo

  elseif(dir.eq.'ydir' .and. jd==1)then
    if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
      je=0
      jf=jd
    else
      je=-jd
      jf=0
    endif
    do j=jaa,jzz
      do i=ia,iz
        do k=1,m1
          vctr1(k)=topo(i,j-jd)+z(k)*rtg(i,j-jd)
          vctr2(k)=topo(i,j)   +z(k)*rtg(i,j)
          vctr3(k)=topo(i,j+jd)+z(k)*rtg(i,j+jd)
          dxytem(k)=vc3da(k,i,j-jd)
          dxynu(k)=vc3da(k,i,j+jd)
        enddo

        !extrapolate the gradient from 3--->2 to 2--->1.
        delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
        delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
        dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
        dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

        CALL htint (nz,dxytem,vctr1,nz,zintl,vctr2)
        CALL htint (nz,dxynu, vctr3,nz,zintr,vctr2)

        !check intersection of cartesian sfc. on left side of mtn.
        CALL topobnd (m1,m2,m3,i,j,0, 1,'l',vc3da,dxy(i,j+je)  &
                    ,vctr1,vctr2,vctr3,zintr,dxytem)
        !check intersection of cartesian sfc. on right side of mtn.
        CALL topobnd (m1,m2,m3,i,j,0,-1,'r',vc3da,dxy(i,j+jf)  &
                    ,vctr1,vctr2,vctr3,zintl,dxynu)

        !Saleeby(2009):Limit lapse rate to prevent runaway cooling
        if (sclrname == 'THP') then
          lapse_rate_max = 0.010
          if (zintr(2) < vc3da(2,i,j+je)) then
             dthdz = (vc3da(2,i,j+je)-zintr(2))/(vctr3(2)-vctr3(1))
             if (dthdz > lapse_rate_max) dthdz = lapse_rate_max
             zintr(2) = vc3da(2,i,j+je) - dthdz*(vctr3(2)-vctr3(1))
          endif
          if (zintl(2) < vc3da(2,i,j+jf)) then
             dthdz = (vc3da(2,i,j+jf)-zintl(2))/(vctr1(2)-vctr1(1))
             if (dthdz > lapse_rate_max) dthdz = lapse_rate_max
             zintl(2) = vc3da(2,i,j+jf) - dthdz*(vctr1(2)-vctr1(1))
          endif
          !if(ngrid==3.and.i==50.and.j==50) &
          !print*,'Lapse1:',sclrname,ngrid,zintr(2),zintl(2)
        endif

        do k=2,m1
          distnu=1./dxytem(k)+1./dxynu(k)
          distold=1./dxy(i,j+je)+1./dxy(i,j+jf)
          scale=distnu**2/distold**2
          if(scale.gt.1)stop ' scale incorrect-ydir--'
          !Need to multiply hkh by 2 since going over 2dx
          aknu=dn0hkh(k,i,j)*scale*2.
          akdiff=-aknu/distnu
          vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                              -(vc3da(k,i,j)-zintl(k))*dxynu(k) )
        enddo

        vc3db(1,i,j)=vc3db(2,i,j)
        if(gpnt.eq.'npnt')then
          vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
        else
          vc3db(m1,i,j)=vc3db(m1-1,i,j)
        endif
      enddo
    enddo

  endif

return
END SUBROUTINE truhor

!##############################################################################
Subroutine topobnd (m1,m2,m3,i,j,ipm1,jpm1,sid,vc3da,dxy  &
                  ,vctr1,vctr2,vctr3,zint,dxynu)

implicit none

integer :: m1,m2,m3,i,j,ipm1,jpm1,k

real :: dxy,tanth,pct,delz1,delz2,delz3,vc1,vc3,vc2
real, dimension(m1,m2,m3) :: vc3da
real, dimension(m1) :: vctr1,vctr2,zint,dxynu,vctr3
character(len=*) :: sid

!     keep track of the distance between scalar points.
CALL ae0 (m1,dxynu,dxy)

!                        |
!                    ----|----
!                   /    |    \
!                  /     |     \
!                 /      |      \
!                /       |       \
!               /        |        \
!              /         |         \
!             /          |          \
!            /           |           \
!      ------            |            -------
!     |      |      |    |   |     |     |
!     vctr1  |   vctr3   | vctr1   |     vctr3
!           vctr2        |       vctr2
!     t(i-1) t(i) t(i+1) | t(i-1) t(i)   t(i+1)
!
!     if the cartesian slice intersects topo.,
!       find the slope of the topography (tanth), then calculate the
!       horizontal cartesian distance to sigma=2 level (1/dxynu).
!       then interpolate along the sigma surface to find the value of
!       the variable at this intersection.

!     check left side of hill for rhs intersection.
if(sid.eq.'l')then

  do k=3,m1
    if(vctr3(2).le.vctr2(k))goto 10
    tanth=(vctr3(2)-vctr2(2))*dxy
    dxynu(k)=tanth/(vctr2(k)-vctr2(2))
    pct=dxy/dxynu(k)
    zint(k)=vc3da(2,i+ipm1,j+jpm1)*pct+vc3da(2,i,j)*(1.-pct)
  enddo
10     continue
  if(vctr3(1).le.vctr2(2))return

  tanth=(vctr3(1)-vctr2(1))*dxy
  dxynu(2)=tanth/(vctr2(2)-vctr2(1))
  pct=dxy/dxynu(2)

!       extrapolate the gradient from 3--->2 to 2--->1.
  delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
  delz2=(vctr2(2)-vctr2(1))/(vctr2(3)-vctr2(2))
  delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
  vc1=vc3da(2,i-ipm1,j-jpm1)-(vc3da(3,i-ipm1,j-jpm1)  &
                             -vc3da(2,i-ipm1,j-jpm1))*delz1
  vc3=vc3da(2,i+ipm1,j+jpm1)-(vc3da(3,i+ipm1,j+jpm1)  &
                             -vc3da(2,i+ipm1,j+jpm1))*delz3
!-----------------------------------------------------------------------
!       use the first vc2 formulation for strictly hh runs, as this is
!         the 'correct' value to use as an endpoint when interpolating
!         along a sigma surface. the second formulation for vc2 is a
!         'fix' around the excessive gridpoint cooling caused by
!         extrapolating a locally cold gridpoint to k=1, and then using
!         this excessively cold value as an endpoint. this then forces
!         the gridpoint cooler at k=2 since the interpolation along the
!         k=1 sigma surface will results in a value colder than the k=2
!         gridpoint value. while a hh run with radiate and sfclyr
!         commented out, lsflg=3, nudlat=0, initial=1 and us,vs=0 will
!         not stay hh, it appears that simulations with the standard
!         rturb.f and this version produce very similar results. this
!         gives confidence to the solution when using the substantially
!         higher resolution allowed by this formulation.
!       note that this scheme does not conserve energy in its present
!         formulation since it diffuses according to a simple 1:2:1
!         scheme.
!       it also appears that while vertical resolution can be as high
!         as 20m with no winds, experiments with the real time
!         forecasting system indicate a resolution nearer 70m is needed
!         since problems with horizontal gradients in a sigma coordinate
!         system near steep topography exist in other routines as well.
!-----------------------------------------------------------------------
!        vc2=vc3da(2,i,j)-(vc3da(3,i,j)-vc3da(2,i,j))*delz2
  vc2=0.5*(vc3+vc1)
  zint(2)=vc3*pct+vc2*(1.-pct)

!     check right side of hill for lhs intersection.
elseif(sid.eq.'r')then

  do k=3,m1
    if(vctr1(2).le.vctr2(k))goto 20
    tanth=(vctr1(2)-vctr2(2))*dxy
    dxynu(k)=tanth/(vctr2(k)-vctr2(2))
    pct=dxy/dxynu(k)
    zint(k)=vc3da(2,i+ipm1,j+jpm1)*pct+vc3da(2,i,j)*(1.-pct)
  enddo
20     continue
  if(vctr1(1).le.vctr2(2))return

  tanth=(vctr1(1)-vctr2(1))*dxy
  dxynu(2)=tanth/(vctr2(2)-vctr2(1))
  pct=dxy/dxynu(2)

!       extrapolate the gradient from 3--->2 to 2--->1.
  delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
  delz2=(vctr2(2)-vctr2(1))/(vctr2(3)-vctr2(2))
  delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
  vc1=vc3da(2,i+ipm1,j+jpm1)-(vc3da(3,i+ipm1,j+jpm1)  &
                             -vc3da(2,i+ipm1,j+jpm1))*delz1
  vc3=vc3da(2,i-ipm1,j-jpm1)-(vc3da(3,i-ipm1,j-jpm1)  &
                             -vc3da(2,i-ipm1,j-jpm1))*delz3
!        vc2=vc3da(2,i,j)-(vc3da(3,i,j)-vc3da(2,i,j))*delz2
  vc2=0.5*(vc3+vc1)
  zint(2)=vc1*pct+vc2*(1.-pct)

endif

return
END SUBROUTINE topobnd

!##############################################################################
Subroutine avgvel (m1,m2,m3,ia,iz,ja,jz,dir,jd,a,b)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd,jaa,jzz,i,j,k
real, dimension(m1,m2,m3) :: a,b
character(len=*) :: dir

jaa=ja
jzz=jz
if(jd.eq.0) then
   jaa=1
   jzz=1
endif

if(dir.eq.'xdir')then
  do j=jaa,jzz
    do i=ia,iz
      do k=1,m1
        a(k,i,j)=0.5*(b(k,i,j)+b(k,i+1,j))
      enddo
    enddo
  enddo
elseif(dir.eq.'ydir')then
  do j=jaa,jzz
    do i=ia,iz
      do k=1,m1
        a(k,i,j)=0.5*(b(k,i,j)+b(k,i,j+jd))
      enddo
    enddo
  enddo
endif

return
END SUBROUTINE avgvel
