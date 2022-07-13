!##############################################################################
Subroutine latbnd ()

use mem_tend
use mem_basic
use mem_grid
use node_mod

implicit none

!     This routine drives the computation of the radiative lateral
!     boundary condition for normal velocity on the coarsest grid
!     and the recomputation of all boundary tendencies on nested grids
!     after the first nested grid timestep.

if (nxtnest(ngrid) .eq. 0) then

!         Radiative and/or mesoscale compensation region lateral
!            boundary conditions.

   if (ibnd .eq. 1 .or. jbnd .eq. 1) then

      CALL latnormv (mzp,mxp,myp,ibcon                  &
         ,basic_g(ngrid)%up  (1,1,1)  ,tend%ut            (1)   &      
         ,basic_g(ngrid)%vp  (1,1,1)  ,tend%vt            (1)   &
         ,grid_g(ngrid)%dxt  (1,1)    ,grid_g(ngrid)%dyt  (1,1) )

   endif

endif

return
END SUBROUTINE latbnd

!##############################################################################
Subroutine latnormv (m1,m2,m3,ibcon,up,ut,vp,vt,dxt,dyt)
   
use mem_grid

implicit none

integer :: m1,m2,m3,ibcon,i,j,k

real :: dxl,dxr,cphx,cphy
real, dimension(m1,m2,m3) :: up,ut,vp,vt
real, dimension(m2,m3) :: dxt,dyt

!     This routine ultimately updates tendencies at lateral boundaries
!     after first diagnosing appropriate phase speeds.
!
!     IBND and JBND are flags for the radiative type in the X and Y
!     direction. Their meaning is:
!
!        IBND=1......Klemp-Wilhelmson (1978) type; phase speed given
!                    by CPHAS
!
!     If this is the first call to a routine, initialize the phase
!       speed arrays if necessary.

if (ibcon.eq.0) return

!     first compute "X" boundaries.

if (iand(ibcon,1) .ne. 0) then
   do j = 1,m3
      dxl = 1. / (dtlv * dxt(2,j))
      do k = 2,m1

         cphx = min(0.,max(-dxl,(up(k,1,j)-cphas)))
         ut(k,1,j) = ut(k,1,j) - cphx * dxt(2,j)  &
            * (up(k,2,j) + ut(k,2,j) * dtlv - up(k,1,j))

      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = 1,m3
      dxr = 1. / (dtlv * dxt(m2-1,j))
      do k = 2,m1

         cphx = max(0.,min(dxr,(up(k,m2-1,j)+cphas)))
         ut(k,m2-1,j) = ut(k,m2-1,j) - cphx * dxt(m2-1,j)  &
            * (up(k,m2-1,j) - (up(k,m2-2,j) + ut(k,m2-2,j) * dtlv))

      enddo
   enddo
endif

!     South and north boundaries.

if (jdim .eq. 1) then

   if (iand(ibcon,4) .ne. 0) then
      do i = 1,m2
         dxl = 1. / (dtlv * dyt(i,2))
         do k = 2,m1
            cphy = min(0.,max(-dxl,(vp(k,i,1)-cphas)))
            vt(k,i,1) = vt(k,i,1) - cphy * dyt(i,2)  &
               * (vp(k,i,2) + vt(k,i,2) * dtlv - vp(k,i,1))
         enddo
      enddo
   endif

   if (iand(ibcon,8) .ne. 0) then
      do i = 1,m2
         dxr = 1. / (dtlv * dyt(i,m3-1))
         do k = 2,m1
            cphy = max(0.,min(dxr,(vp(k,i,m3-1)+cphas)))
            vt(k,i,m3-1) = vt(k,i,m3-1) - cphy * dyt(i,m3-1)  &
               * (vp(k,i,m3-1) - (vp(k,i,m3-2) + vt(k,i,m3-2) * dtlv))
         enddo
      enddo
   endif

endif

return
END SUBROUTINE latnormv

!##############################################################################
Subroutine update_cyclic (isflag)

! This routine will update cyclic boundary values for either sequential or
! parallel modes. The argument isflag denotes which variables to update.
!
! update_cyclic is called from trsets() for updating scalar vars, and from 
! acoust()  for updating u, v, p, and w.

use var_tables
use mem_grid
use mem_basic
use mem_turb
use node_mod

implicit none

  integer :: isflag

  integer :: ivar
  real, pointer :: scalarp

  if (nmachs .eq. 1) then
    ! sequential run
    select case (isflag)
      case(LBC_ALL_SCALARS)
        do ivar = 1,num_scalar(1)      
          scalarp => scalar_tab(ivar,1)%var_p
          CALL cyclic (nnzp(1),nnxp(1),nnyp(1),scalarp,'T')
        enddo
      case(LBC_UP)
        CALL cyclic (nnzp(1),nnxp(1),nnyp(1),basic_g(1)%up,'U')
      case(LBC_VP)
        CALL cyclic (nnzp(1),nnxp(1),nnyp(1),basic_g(1)%vp,'V')
      case(LBC_PP)
        CALL cyclic (nnzp(1),nnxp(1),nnyp(1),basic_g(1)%pp,'P')
      case(LBC_WP)
        CALL cyclic (nnzp(1),nnxp(1),nnyp(1),basic_g(1)%wp,'W')
      case default
        print*, 'ERROR: update_cyclic: isflag value out of allowed'
        print*, 'range (1-5): ', isflag
        stop
    endselect
  else
    ! parallel run
    
    ! Want to mimic the algorithm in cyclic() which first processes in 
    ! the x-direction followed by the y-direction. The second arg to 
    ! sendcyclic/getcyclic specifies the direction (1 -> x, 2 -> y).
    if (ibnd .eq. 2) then
      CALL node_sendcyclic (isflag,1)
      CALL node_getcyclic (isflag,1)

      ! update overlap regions
      select case (isflag)
        case(LBC_ALL_SCALARS)
          CALL update_lbc_vgroup (ngrid,LBC_ALL_SCALARS)
        case(LBC_UP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%up,3)
        case(LBC_VP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%vp,3)
        case(LBC_PP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%pp,3)
        case(LBC_WP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%wp,3)
      endselect
    endif

    if ((jbnd .eq. 2) .and. (jdim .eq. 1)) then
      CALL node_sendcyclic (isflag,2)
      CALL node_getcyclic (isflag,2)

      ! update overlap regions
      select case (isflag)
        case(LBC_ALL_SCALARS)
          CALL update_lbc_vgroup (ngrid,LBC_ALL_SCALARS)
        case(LBC_UP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%up,3)
        case(LBC_VP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%vp,3)
        case(LBC_PP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%pp,3)
        case(LBC_WP)
          CALL update_lbc_var (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid),ngrid &
                             ,basic_g(ngrid)%wp,3)
      endselect
    endif
  endif

  return
END SUBROUTINE update_cyclic

!##############################################################################
Subroutine cyclic (n1,n2,n3,var,vpnt)

use mem_grid

implicit none

integer :: n1,n2,n3,ncycle,nshift,i,j,k
real, dimension(n1,n2,n3) :: var
character(len=*) :: vpnt

!     This routine inputs a 3d variable VAR dimensioned by N1,N2,N3
!     and sets the boundaries to cyclic symmetry.  This version is
!     set up for second order advection only.

if (ibnd .eq. 2) then
   ncycle = n2 - 3
   nshift = 0
   if (vpnt .eq. 'U') nshift = 1
   do j = 1,n3
      do k = 1,n1
         var(k,1,j) = var(k,ncycle+1,j)
         var(k,n2-nshift,j) = var(k,n2-nshift-ncycle,j)
      enddo
   enddo

   if (vpnt .ne. 'U') then
      do j = 1,n3
         do k = 1,n1
            var(k,2,j) = 0.5 * (var(k,2,j) + var(k,n2-1,j))
            var(k,n2-1,j) = var(k,2,j)
         enddo
      enddo
   endif
endif

if (jbnd .eq. 2 .and. jdim .eq. 1) then
   ncycle = n3 - 3
   nshift = 0
   if (vpnt .eq. 'V') nshift = 1
   do i = 1,n2
      do k = 1,n1
         var(k,i,1) = var(k,i,ncycle+1)
         var(k,i,n3-nshift) = var(k,i,n3-nshift-ncycle)
      enddo
   enddo

   if (vpnt .ne. 'V') then
      do i = 1,n2
         do k = 1,n1
            var(k,i,2) = 0.5 * (var(k,i,2) + var(k,i,n3-1))
            var(k,i,n3-1) = var(k,i,2)
         enddo
      enddo
   endif

endif

return
END SUBROUTINE cyclic

!##############################################################################
Subroutine vpsets ()

use mem_basic
use mem_grid
use node_mod

implicit none

if (nxtnest(ngrid) .eq. 0) then
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'              &
      ,basic_g(ngrid)%up (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'              &
      ,basic_g(ngrid)%vp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'W'              &
      ,basic_g(ngrid)%wp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'              &
      ,basic_g(ngrid)%pp (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
endif

if (nsttop .eq. 1) then
   CALL topset (mzp,mxp,myp,basic_g(ngrid)%up(1,1,1),'U')
   CALL topset (mzp,mxp,myp,basic_g(ngrid)%vp(1,1,1),'V')
endif

if (nstbot .eq. 1) then
      CALL botset (mzp,mxp,myp,basic_g(ngrid)%up(1,1,1),'U')
      CALL botset (mzp,mxp,myp,basic_g(ngrid)%vp(1,1,1),'V')
endif

CALL topset (mzp,mxp,myp,basic_g(ngrid)%pp(1,1,1),'P')

CALL botset (mzp,mxp,myp,basic_g(ngrid)%pp(1,1,1),'P')

CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%wp(1,1,1),'W')
CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%up(1,1,1),'U')
CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%vp(1,1,1),'V')

if (nxtnest(ngrid) .eq. 0) then
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'  &
      ,basic_g(ngrid)%uc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'  &
      ,basic_g(ngrid)%vc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'W'  &
      ,basic_g(ngrid)%wc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
   CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'  &
      ,basic_g(ngrid)%pc (1,1,1)  ,basic_g(ngrid)%up (1,1,1)  &
      ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
      ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
      ,grid_g(ngrid)%dym (1,1)    ,'NULL')
endif

if (nsttop .eq. 1) then
   CALL topset (mzp,mxp,myp,basic_g(ngrid)%uc(1,1,1),'U')
   CALL topset (mzp,mxp,myp,basic_g(ngrid)%vc(1,1,1),'V')
endif

if (nstbot .eq. 1) then
      CALL botset (mzp,mxp,myp,basic_g(ngrid)%uc(1,1,1),'U')
      CALL botset (mzp,mxp,myp,basic_g(ngrid)%vc(1,1,1),'V')
endif

CALL topset (mzp,mxp,myp,basic_g(ngrid)%pc(1,1,1),'P')

CALL botset (mzp,mxp,myp,basic_g(ngrid)%pc(1,1,1),'P')

CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%wc(1,1,1),'W')
CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%uc(1,1,1),'U')
CALL dumset (mzp,mxp,myp,ibcon,basic_g(ngrid)%vc(1,1,1),'V')

return
END SUBROUTINE vpsets

!##############################################################################
Subroutine trsets ()

use var_tables
use mem_basic
use mem_grid
use mem_turb
use node_mod

implicit none

integer :: n
real, pointer :: scalarp, scalart

!     Apply lateral, top, and bottom boundary conditions.

! Update scalar overlap regions
if (nmachs .gt. 1) CALL update_lbc_vgroup (ngrid,LBC_ALL_SCALARS)

! Update scalar cyclic boundaries
if (ngrid .eq. 1) CALL update_cyclic (LBC_ALL_SCALARS)

do n = 1,num_scalar(ngrid)
   scalarp => scalar_tab(n,ngrid)%var_p
   scalart => scalar_tab(n,ngrid)%var_t

   if (nxtnest(ngrid) .eq. 0) then
      CALL latset (mzp,mxp,myp,ia,iz,ja,jz,ibcon,'TR'           &
        ,scalarp                    ,basic_g(ngrid)%up (1,1,1)  &
        ,basic_g(ngrid)%vp (1,1,1)  ,grid_g(ngrid)%dxu (1,1)    &
        ,grid_g(ngrid)%dxm (1,1)    ,grid_g(ngrid)%dyv (1,1)    &
        ,grid_g(ngrid)%dym (1,1)    ,scalar_tab(n,ngrid)%name)
   endif

   if (nsttop .eq. 1)  &
      CALL topset (mzp,mxp,myp,scalarp,'T')
   if (nstbot .eq. 1)  then
      CALL botset (mzp,mxp,myp,scalarp,'T')
   endif
enddo

!       Make sure all positive definite quantities remain such.

CALL tkeinit (mzp,mxp,myp)

return
END SUBROUTINE trsets

!##############################################################################
Subroutine latset (m1,m2,m3,ia,iz,ja,jz,ibcon,vnam,ap,uc,vc,dxu,dxm,dyv,dym &
                   ,sname)

use mem_grid
use micphys
use mem_leaf, only:isfcl

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,lbw,lbe,lbs,lbn,setlbcinit
real :: thresh,dtlx,c1,dxr,dyr,atau
real, dimension(m1,m2,m3) :: ap,uc,vc
real, dimension(m2,m3) :: dxu,dxm,dyv,dym
character(len=*) :: vnam,sname
real, dimension(:), allocatable :: temp_v1, temp_v2 

lbw=0 !Variable initialized
lbe=0 !Variable initialized
lbs=0 !Variable initialized
lbn=0 !Variable initialized

if (iand(ibcon,1) .gt. 0) lbw = ia - 1
if (iand(ibcon,2) .gt. 0) lbe = iz + 1
if (iand(ibcon,4) .gt. 0) lbs = ja - 1
if (iand(ibcon,8) .gt. 0) lbn = jz + 1

thresh = 0.
if (vnam.eq.'U' .or. vnam.eq.'V' .or. vnam.eq.'W' .or. vnam.eq.'P') then
   dtlx = dtlv
else
   dtlx = dtlt
endif

!This is in place to be able to set 1 interior boundary points
!toward the initial conditions for certain variables so we do not
!lose certain fields over time. Only used for zero gradient or
!radiative BC (LSFLG = 0 or 1).
setlbcinit=0
atau=max(bctau(ngrid),dtlt) !prevent from being less than dtlt

!LBC reset for aerosols
if(iaerolbc(ngrid)==1)then
 if( ((sname == 'CIFNP') .and. &
       (jnmb(3)>=5 .and. (iifn==1.or.iifn==2))) .or. &
     ((sname == 'CCCNP' .or. sname == 'CCCMP' .or. &
       sname == 'GCCNP' .or. sname == 'GCCMP') .and. &
       iaerosol==1) .or. &
     ((sname == 'MD1NP' .or. sname == 'MD1MP' .or. &
       sname == 'MD2NP' .or. sname == 'MD2MP') .and. &
       idust==1) .or. &
     ((sname == 'ABC1NP' .or. sname == 'ABC1MP' .or. &
       sname == 'ABC2NP' .or. sname == 'ABC2MP') .and. &
       iabcarb==1) .or. &
     ((sname == 'SALT_FILM_NP' .or. sname == 'SALT_FILM_MP' .or. &
       sname == 'SALT_JET_NP'  .or. sname == 'SALT_JET_MP'  .or. &
       sname == 'SALT_SPUM_NP' .or. sname == 'SALT_SPUM_MP') .and. &
       isalt==1) ) then
    setlbcinit=1
 endif
endif

!LBC reset for CO2 if SiB land surface model is used
if(ico2lbc(ngrid)==1)then
 if(sname == 'RCO2P' .and. isfcl==2) then
    setlbcinit=1
 endif
endif

allocate(temp_v1(m1))
allocate(temp_v2(m1))

if (ibnd .ne. 2 .and. vnam .ne. 'U' .and. lsflg .ne. 3) then

!     Western and Eastern boundaries for zero gradient option

   if (lsflg .eq. 0) then
      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               if(setlbcinit==1) then
                ap(k,ia,j) = ap(k,ia,j)+(ap(k,lbw,j)-ap(k,ia,j))*(dtlt/atau)
               endif
               ap(k,lbw,j) = ap(k,ia,j)
            enddo
         enddo
      endif
      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               if(setlbcinit==1) then
                ap(k,iz,j) = ap(k,iz,j)+(ap(k,lbe,j)-ap(k,iz,j))*(dtlt/atau)
               endif
               ap(k,lbe,j) = ap(k,iz,j)
            enddo
         enddo
      endif
   else

!     Western boundary for lsflg = 1 or 2

      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3
            if (vnam .eq. 'V') then
               dxr = dxm(ia,j) / dxm(lbw,j)
               c1 = .5 * dtlx * dxm(lbw,j)
               do k = 1,m1
                 if(j==m3)then
                  temp_v1(k) = -c1 * uc(k,lbw,j)
                 else
                  temp_v1(k) = -c1 * (uc(k,lbw,j) + uc(k,lbw,j+jdim))
                 endif
!if(k==2.and.j==m3)print*,'westbc-V',lbw,temp_v1(k),uc(k,lbw,j),uc(k,lbw,j+jdim)
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = .5 * dtlx * dxu(lbw,j)
               do k = 1,m1
                 if(k==m1)then
                  temp_v1(k) = -c1 * uc(k,lbw,j)
                 else
                  temp_v1(k) = -c1 * (uc(k,lbw,j) + uc(k+1,lbw,j))
                 endif
!if(k==m1.and.j==m3)print*,'westbc-W',lbw,temp_v1(k),uc(k,lbw,j),uc(k+1,lbw,j)
               enddo
            else
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = dtlx * dxu(lbw,j)
               do k = 1,m1
                  temp_v1(k) = -c1 * uc(k,lbw,j)
               enddo
            endif
            do k = 1,m1
               temp_v2(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
            enddo
            do k = 1,m1
               if (temp_v1(k) .ge. thresh) then
                  if(setlbcinit==1) then
                   ap(k,ia,j) = ap(k,ia,j)+(ap(k,lbw,j)-ap(k,ia,j))*(dtlt/atau)
                   ap(k,lbw,j) = ap(k,ia,j)
                  else !radiative BC
                   ap(k,lbw,j) = temp_v2(k)
                  endif
               elseif (lsflg .eq. 1) then
                  if(setlbcinit==1) then
                   ap(k,ia,j) = ap(k,ia,j)+(ap(k,lbw,j)-ap(k,ia,j))*(dtlt/atau)
                  endif
                  ap(k,lbw,j) = ap(k,ia,j)
               endif
            enddo
         enddo
      endif

!     Eastern Boundary for LSFLG = 1 or 2

      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3
            if (vnam .eq. 'V') then
               dxr = dxm(iz-1,j) / dxm(iz,j)
               c1 = .5 * dtlx * dxm(iz,j)
               do k = 1,m1
                 if(j==m3)then
                  temp_v1(k) = c1 * uc(k,iz,j)
                 else
                  temp_v1(k) = c1 * (uc(k,iz,j) + uc(k,iz,j+jdim))
                 endif
!if(k==2.and.j==m3)print*,'eastbc-V',iz,temp_v1(k),uc(k,iz,j),uc(k,iz,j+jdim)
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = .5 * dtlx * dxu(iz,j)
               do k = 1,m1
                 if(k==m1)then
                  temp_v1(k) = c1 * uc(k,iz,j)
                 else
                  temp_v1(k) = c1 * (uc(k,iz,j) + uc(k+1,iz,j))
                 endif
!if(k==m1.and.j==m3)print*,'eastbc-W',iz,temp_v1(k),uc(k,iz,j),uc(k+1,iz,j)
               enddo
            else
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = dtlx * dxu(iz,j)
               do k = 1,m1
                  temp_v1(k) = c1 * uc(k,iz,j)
               enddo
            endif
            do k = 1,m1
               temp_v2(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
            enddo
            do k = 1,m1
               if (temp_v1(k) .ge. thresh) then
                  if(setlbcinit==1) then
                   ap(k,iz,j) = ap(k,iz,j)+(ap(k,lbe,j)-ap(k,iz,j))*(dtlt/atau)
                   ap(k,lbe,j) = ap(k,iz,j)
                  else !radiative BC
                   ap(k,lbe,j) = temp_v2(k)
                  endif
               elseif (lsflg .eq. 1) then
                  if(setlbcinit==1) then
                   ap(k,iz,j) = ap(k,iz,j)+(ap(k,lbe,j)-ap(k,iz,j))*(dtlt/atau)
                  endif
                  ap(k,lbe,j) = ap(k,iz,j)
               endif
            enddo
         enddo
      endif
   endif
endif

if(jdim.eq.1.and.jbnd.ne.2.and.vnam.ne.'V'.and.lsflg.ne.3)then

!     Southern and Northern boundaries for zero gradient option

  if (lsflg .eq. 0) then
     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              if(setlbcinit==1) then
               ap(k,i,ja) = ap(k,i,ja)+(ap(k,i,lbs)-ap(k,i,ja))*(dtlt/atau)
              endif
              ap(k,i,lbs) = ap(k,i,ja)
           enddo
        enddo
     endif
     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              if(setlbcinit==1) then
               ap(k,i,jz) = ap(k,i,jz)+(ap(k,i,lbn)-ap(k,i,jz))*(dtlt/atau)
              endif
              ap(k,i,lbn) = ap(k,i,jz)
           enddo
        enddo
     endif
  else

!     Southern boundary for LSFLG = 1 or 2

     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2
           if (vnam .eq. 'U') then
              dyr = dym(i,ja) / dym(i,lbs)
              c1 = .5 * dtlx * dym(i,lbs)
              do k = 1,m1
                if(i==m2)then
                 temp_v1(k) = -c1 * vc(k,i,lbs)
                else
                 temp_v1(k) = -c1 * (vc(k,i,lbs) + vc(k,i+1,lbs))
                endif
!if(k==2.and.i==m2)print*,'southbc-U',lbs,temp_v1(k),vc(k,i,lbs),vc(k,i+1,lbs)
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = .5 * dtlx * dyv(i,lbs)
              do k = 1,m1
                if(k==m1)then
                 temp_v1(k) = -c1 * vc(k,i,lbs)
                else
                 temp_v1(k) = -c1 * (vc(k,i,lbs) + vc(k+1,i,lbs))
                endif
!if(k==m1.and.i==m2)print*,'southbc-W',lbs,temp_v1(k),vc(k,i,lbs),vc(k+1,i,lbs)
              enddo
           else
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = dtlx * dyv(i,lbs)
              do k = 1,m1
                 temp_v1(k) = -c1 * vc(k,i,lbs)
              enddo
           endif
           do k = 1,m1
              temp_v2(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
           enddo
           do k = 1,m1
              if (temp_v1(k) .ge. thresh) then
                 if(setlbcinit==1) then
                  ap(k,i,ja) = ap(k,i,ja)+(ap(k,i,lbs)-ap(k,i,ja))*(dtlt/atau)
                  ap(k,i,lbs) = ap(k,i,ja)                
                 else !radiative BC
                  ap(k,i,lbs) = temp_v2(k)
                 endif
              elseif (lsflg .eq. 1) then
                 if(setlbcinit==1) then
                  ap(k,i,ja) = ap(k,i,ja)+(ap(k,i,lbs)-ap(k,i,ja))*(dtlt/atau)
                 endif
                 ap(k,i,lbs) = ap(k,i,ja)
              endif
           enddo
        enddo
     endif

!     Northern Boundary for LSFLG = 1 or 2

     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2
           if (vnam .eq. 'U') then
              dyr = dym(i,jz-1) / dym(i,jz)
              c1 = .5 * dtlx * dym(i,jz)
              do k = 1,m1
                if(i==m2)then
                 temp_v1(k) = c1 * vc(k,i,jz)
                else
                 temp_v1(k) = c1 * (vc(k,i,jz) + vc(k,i+1,jz))
                endif
!if(k==2.and.i==m2)print*,'northbc-U',jz,temp_v1(k),vc(k,i,jz),vc(k,i+1,jz)
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = .5 * dtlx * dyv(i,jz)
              do k = 1,m1
                if(k==m1)then
                 temp_v1(k) = c1 * vc(k,i,jz)
                else
                 temp_v1(k) = c1 * (vc(k,i,jz) + vc(k+1,i,jz))
                endif
!if(k==m1.and.i==m2)print*,'northbc-W',jz,temp_v1(k),vc(k,i,jz),vc(k+1,i,jz)
              enddo
           else
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = dtlx * dyv(i,jz)
              do k = 1,m1
                 temp_v1(k) = c1 * vc(k,i,jz)
              enddo
           endif
           do k = 1,m1
              temp_v2(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
           enddo
           do k = 1,m1
              if (temp_v1(k) .ge. thresh) then
                 if(setlbcinit==1) then
                  ap(k,i,jz) = ap(k,i,jz)+(ap(k,i,lbn)-ap(k,i,jz))*(dtlt/atau)
                  ap(k,i,lbn) = ap(k,i,jz)
                 else !radiative BC
                  ap(k,i,lbn) = temp_v2(k)
                 endif
              elseif (lsflg .eq. 1) then
                 if(setlbcinit==1) then
                  ap(k,i,jz) = ap(k,i,jz)+(ap(k,i,lbn)-ap(k,i,jz))*(dtlt/atau)
                 endif
                 ap(k,i,lbn) = ap(k,i,jz)
              endif
           enddo
        enddo
     endif
  endif
endif

deallocate(temp_v1)
deallocate(temp_v2)

return
END SUBROUTINE latset

!##############################################################################
Subroutine topset (m1,m2,m3,ap,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,i,j
real :: dzmr,dztr
real, dimension(m1,m2,m3) :: ap
character(len=*) :: vnam

dzmr = dzm(m1-2) / dzm(m1-1)
dztr = dzt(m1-2) / dzt(m1-1)

!     Computation of all prognostic variables (other than W) at
!       level NZP by extrapolation from below

if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
      enddo
   enddo
endif
if (vnam .eq. 'T') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = max(0.,ap(m1-1,i,j)+dzmr*(ap(m1-1,i,j)-ap(m1-2,i,j)))
      enddo
   enddo
endif

return
END SUBROUTINE topset

!##############################################################################
Subroutine botset (m1,m2,m3,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,i,j
real :: dzmr
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'P') then
   dzmr = dzm(2) / dzm(1)
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j) + (aa(2,i,j) - aa(3,i,j)) * dzmr
      enddo
   enddo
else
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j)
      enddo
   enddo
endif

return
END SUBROUTINE botset

!##############################################################################
Subroutine dumset (m1,m2,m3,ibcon,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ibcon,i,j,k
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'U' .and. iand(ibcon,2) .gt. 0) then
   do j = 1,m3
      do k = 1,m1
         aa(k,m2,j) = aa(k,m2-1,j)
      enddo
   enddo
elseif (vnam .eq. 'V' .and. iand(ibcon,8) .gt. 0) then
   do i = 1,m2
      do k = 1,m1
         aa(k,i,m3) = aa(k,i,m3-jdim)
      enddo
   enddo
elseif (vnam .eq. 'W') then
   do j = 1,m3
      do i = 1,m2
         aa(m1,i,j) = aa(m1-1,i,j)
      enddo
   enddo
endif

return
END SUBROUTINE dumset

!##############################################################################
Subroutine rayft ()

use mem_tend
use mem_basic
use mem_grid
use node_mod
use micphys

implicit none

integer :: i,j,k
real, dimension(:,:,:), allocatable :: theta_v

!     This routine is the rayleigh friction driver for the
!     theta friction and is called from the long timestep.

if (nfpt .eq. 0 .or. distim .eq. 0.) return

!     First load past virtual theta into temporary.
allocate(theta_v(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

if (level .ge. 1) then
   do j = 1,mmyp(ngrid)
      do i = 1,mmxp(ngrid)
         do k = 1,mmzp(ngrid)
            theta_v(k,i,j) = basic_g(ngrid)%theta(k,i,j)  &
               * (1. + .61 * basic_g(ngrid)%rv(k,i,j))
         enddo
      enddo
   enddo
else
  CALL atob (mxyzp,basic_g(ngrid)%theta(1,1,1),theta_v(1,1,1))
endif

!     Now get rayleigh friction tendency

 CALL rayf (4,mzp,mxp,myp,ia,iz,ja,jz                         &
      ,theta_v            (1,1,1) ,basic_g(ngrid)%th0 (1,1,1)  &
      ,tend%tht           (1)     ,grid_g(ngrid)%rtgt (1,1)    &
      ,grid_g(ngrid)%topt (1,1)                                )
      
deallocate(theta_v)

return
END SUBROUTINE rayft

!##############################################################################
Subroutine rayf (ifrom,m1,m2,m3,ia,iz,ja,jz,var,th0,tht,rtgx,topx)

use mem_grid
use ref_sounding

implicit none

integer :: ifrom,m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: var,th0,tht
real, dimension(m2,m3) :: rtgx,topx
real, dimension(:), allocatable :: temp_v1, temp_v2

real :: zmkf,c1,c2
integer :: kf,i,j,k

!     This routine calculates rayleigh friction terms velocity and theta_il

if (nfpt .eq. 0 .or. distim .le. 0) return
kf = nnz(1) - nfpt
zmkf = zmn(kf,1)
c1 = 1. / (distim * (ztop - zmkf))
c2 = dts * c1

allocate(temp_v1(nzp))
allocate(temp_v2(nzp))

goto(100,200,300,400) ifrom
100   continue

!     u friction

do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         temp_v1(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      CALL htint (nzp,u01dn(1,ngrid),zt,nzp,temp_v2,temp_v1)
      do k = nz,2,-1
         if (temp_v1(k) .le. zmkf) go to 10

         var(k,i,j) = var(k,i,j) + c2 * (temp_v1(k) - zmkf)  &
            * (temp_v2(k) - var(k,i,j))

      enddo
10         continue
   enddo
enddo

deallocate(temp_v1)
deallocate(temp_v2)

return
200   continue

!     V friction

if (jdim .eq. 0 .and. icorflg .eq. 0) then
  deallocate(temp_v1)
  deallocate(temp_v2)
  return
endif
do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         temp_v1(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      CALL htint (nzp,v01dn(1,ngrid),zt,nzp,temp_v2,temp_v1)
      do k = nz,2,-1
         if (temp_v1(k) .le. zmkf) go to 20
         var(k,i,j) = var(k,i,j) + c2 * (temp_v1(k) - zmkf)  &
            * (temp_v2(k) - var(k,i,j))
      enddo
20         continue
   enddo
enddo

deallocate(temp_v1)
deallocate(temp_v2)

return
300   continue

!     W friction

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         temp_v1(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (temp_v1(k) .le. zmkf) go to 30
         var(k,i,j) = var(k,i,j) - c2 * (temp_v1(k) - zmkf) * var(k,i,j)
      enddo
30         continue
   enddo
enddo

deallocate(temp_v1)
deallocate(temp_v2)

return
400   continue

!     THETA FRICTION

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         temp_v1(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (temp_v1(k) .le. zmkf) go to 40
         tht(k,i,j) = tht(k,i,j) + c1 * (temp_v1(k) - zmkf)  &
              * (th0(k,i,j) - var(k,i,j))
      enddo
40         continue
   enddo
enddo

deallocate(temp_v1)
deallocate(temp_v2)

return
END SUBROUTINE rayf
