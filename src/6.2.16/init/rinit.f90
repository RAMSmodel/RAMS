!##############################################################################
Subroutine fldinit ()

use mem_grid
use mem_basic
use mem_turb
use node_mod

implicit none

integer :: initflg

! finish initializing past time level variables

CALL atob (mxyzp,basic_g(ngrid)%uc(1,1,1),basic_g(ngrid)%up(1,1,1))
CALL atob (mxyzp,basic_g(ngrid)%vc(1,1,1),basic_g(ngrid)%vp(1,1,1))
CALL atob (mxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1))
CALL atob (mxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1))

! initialize TKE to a minimum default value

CALL tkeinit (mzp,mxp,myp)

! set initial lateral boundary conditions for momentum

CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%uc(1,1,1),'U')
CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%vc(1,1,1),'V')
CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%wc(1,1,1),'W')
CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%up(1,1,1),'U')
CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%vp(1,1,1),'V')
CALL dumset (mzp,mxp,myp,15,basic_g(ngrid)%wp(1,1,1),'W')

return
END SUBROUTINE fldinit

!##############################################################################
Subroutine gridloc_prt ()

use mem_grid
use node_mod
use rconstants

implicit none

real :: gscr(maxgrds,19),centx,centy
integer :: ngr,midx,midy,i,j,ng
real, dimension(:), allocatable :: all_gscr

! gscr hold lat/lon extent of the current grid
!   index     value
!     1         sw lat
!     3         se lat
!     5         ne lat
!     7         nw lat
!
!     2         sw lon
!     4         se lon
!     6         ne lon
!     8         nw lon
!     

do ngr=1,ngrids
   CALL newgrid (ngr)
   ! Find extent of sub domain
   gscr(ngr,1)=grid_g(ngr)%glat(1,1)
   gscr(ngr,2)=grid_g(ngr)%glon(1,1)
   gscr(ngr,3)=grid_g(ngr)%glat(mxp,1)
   gscr(ngr,4)=grid_g(ngr)%glon(mxp,1)
   gscr(ngr,5)=grid_g(ngr)%glat(mxp,myp)
   gscr(ngr,6)=grid_g(ngr)%glon(mxp,myp)
   gscr(ngr,7)=grid_g(ngr)%glat(1,myp)
   gscr(ngr,8)=grid_g(ngr)%glon(1,myp)

   ! get the full domain corner values
   !  args: igrid, sw_val, nw_val, se_val, ne_val
   CALL get_fd_corner_vals (ngr,gscr(ngr,1),gscr(ngr,7),gscr(ngr,3) &
                           ,gscr(ngr,5)) ! lat
   CALL get_fd_corner_vals (ngr,gscr(ngr,2),gscr(ngr,8),gscr(ngr,4) &
                           ,gscr(ngr,6)) ! lon

   midx = (nxp+1)/2
   midy = (nyp+1)/2
   if (nxp/2 .eq. midx) then
      centx = xm(nxp/2)
   else
      centx = xt(midx)
   endif
   if (nyp/2 .eq. midy) then
      centy = ym(nyp/2)
   else
      centy = yt(midy)
   endif
   CALL xy_ll (gscr(ngr,9),gscr(ngr,10),polelat,polelon  &
      ,centx,centy)

   ! Find the inner and outer extents of each grid side
   !
   ! Do this by running the calculations on the sub domain. If there
   ! was only one node, then the sub domain matches the full domain
   ! and we are done.
   !
   ! If there were more than one node, then gather all of the calculations
   ! onto mainnum. Then mainnum will continue the calculations by taking 
   ! the min/max of the results coming from all of the other nodes.
   !
   !
   ! Compare the sub domain indices against the full domain values in order
   ! to only look at lat/lon values when on the full domain boundaries.
   ! If a sub domain does not contain part of the full domain boundary, 
   ! then the initial +/- 1000 value will be left in the gscr entry and will
   ! not come into play when mainnum does the min/max of the subdomain results.
   ! Otherwise the sections from all of the sub domains that do contain
   ! a particluar full domain boundary (E, W, S, or N) will be strung together
   ! in a min/max comparison which will result in the correct min/max value
   ! being calculated.
   !
   gscr(ngr,12)=-1000.
   gscr(ngr,13)=1000.
   gscr(ngr,14)=1000.
   gscr(ngr,15)=-1000.
   gscr(ngr,16)=-1000.
   gscr(ngr,17)=1000.
   gscr(ngr,18)=-1000.
   gscr(ngr,19)=1000.
   do j=1,myp
      do i=1,mxp
         if(j+j0.eq.nyp)  &
            gscr(ngr,12)=max(gscr(ngr,12),grid_g(ngr)%glat(i,j))  ! outer s
         if(j+j0.eq.1)  &
            gscr(ngr,13)=min(gscr(ngr,13),grid_g(ngr)%glat(i,j))  ! outer n
         if(i+i0.eq.1)  &
            gscr(ngr,14)=min(gscr(ngr,14),grid_g(ngr)%glon(i,j))  ! outer w
         if(i+i0.eq.nxp)  & 
            gscr(ngr,15)=max(gscr(ngr,15),grid_g(ngr)%glon(i,j))  ! outer e
         ! provide a conservative estimate so interpolation and
         ! round off error don't result in post process nans
         if(j+j0.eq.3)  &
            gscr(ngr,16)=max(gscr(ngr,16),grid_g(ngr)%glat(i,j))  ! inner s
         if(j+j0.eq.nyp-2)  &
            gscr(ngr,17)=min(gscr(ngr,17),grid_g(ngr)%glat(i,j))  ! inner n
         if(i+i0.eq.3)  &
            gscr(ngr,18)=max(gscr(ngr,18),grid_g(ngr)%glon(i,j))  ! inner w
         if(i+i0.eq.nxp-2)  &
            gscr(ngr,19)=min(gscr(ngr,19),grid_g(ngr)%glon(i,j))  ! inner e
       enddo
   enddo

   ! At this point gscr values have been set from a sub domain perspective.
   ! If the are more than 1 node, then gather the results on mainnum and
   ! finish off the min/max calculation. Only main will be printing the
   ! report, so there is no need to broadcast the results back to all of
   ! the other nodes.
   !
   if(nmachs .gt. 1) then
     allocate(all_gscr(nmachs))

     if(my_rams_num .eq. mainnum) then
       ! receive from other nodes
       CALL par_gather_floats (gscr(ngr,12), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer s
       gscr(ngr,12) = maxval(all_gscr)
       CALL par_gather_floats (gscr(ngr,13), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer n
       gscr(ngr,13) = minval(all_gscr)
       CALL par_gather_floats (gscr(ngr,14), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer w
       gscr(ngr,14) = minval(all_gscr)
       CALL par_gather_floats (gscr(ngr,15), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer e
       gscr(ngr,15) = maxval(all_gscr)

       CALL par_gather_floats (gscr(ngr,16), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner s
       gscr(ngr,16) = maxval(all_gscr)
       CALL par_gather_floats (gscr(ngr,17), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner n
       gscr(ngr,17) = minval(all_gscr)
       CALL par_gather_floats (gscr(ngr,18), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner w
       gscr(ngr,18) = maxval(all_gscr)
       CALL par_gather_floats (gscr(ngr,19), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner e
       gscr(ngr,19) = minval(all_gscr)
     else
       ! send to mainnum
       CALL par_gather_floats (gscr(ngr,12), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer s
       CALL par_gather_floats (gscr(ngr,13), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer n
       CALL par_gather_floats (gscr(ngr,14), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer w
       CALL par_gather_floats (gscr(ngr,15), all_gscr, 1 &
                              ,machnum(mainnum)) ! outer e

       CALL par_gather_floats (gscr(ngr,16), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner s
       CALL par_gather_floats (gscr(ngr,17), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner n
       CALL par_gather_floats (gscr(ngr,18), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner w
       CALL par_gather_floats (gscr(ngr,19), all_gscr, 1 &
                              ,machnum(mainnum)) ! inner e
     endif

     deallocate(all_gscr)
   endif

ENDDO

do ng=1,ngrids
  if(print_msg) then
    write(6,100)ng
    write(6,101)gscr(ng,7),gscr(ng,8),gscr(ng,5),gscr(ng,6)
    write(6,102)gscr(ng,9),gscr(ng,10)
    write(6,105)gscr(ng,1),gscr(ng,2),gscr(ng,3),gscr(ng,4)
    print*, ' '
    write(6,106)xtn(1,ng),xtn(nnxp(ng),ng)
    write(6,107)ytn(1,ng),ytn(nnyp(ng),ng)
    write(6,109)deltaxn(ng)
    print*, ' '
    write(6,108)zmn(1,ng),zmn(nnzp(ng)-1,ng)
    write(6,110)zmn(2,ng)-zmn(1,ng)
    print*, ' '
    write(6,113) gscr(ng,12)
    write(6,114) gscr(ng,17)
    write(6,115) gscr(ng,14),gscr(ng,15)
    write(6,116) gscr(ng,18),gscr(ng,19)
    write(6,117) gscr(ng,16)
    write(6,118) gscr(ng,13)
  endif
enddo
!
100  FORMAT(/,'---------------Location and Dimensions of GRID',I4  &
        ,'---------------------------',/)
101  FORMAT('  NW lat/lon      NE lat/lon (deg) = '  &
   F8.4,',',F9.4,4X,F8.4,',',F9.4)
102  FORMAT('       Center lat/lon (deg)        =  '  &
   ,11X,F8.4,',',F9.4)
105  FORMAT('  SW lat/lon      SE lat/lon (deg) = '  &
   F8.4,',',F9.4,4X,F8.4,',',F9.4)
106  FORMAT('     West PS coord (km) =',-3PF10.3  &
     ,'               East PS coord (km) =',-3PF10.3)
107  FORMAT('    South PS coord (km) =',-3PF10.3  &
     ,'              North PS coord (km) =',-3PF10.3)
109  FORMAT('         Delta-X (m) =',F10.1)
108  FORMAT('    Bottom coord (m) =',F10.1  &
     ,'            Top coordinate (m) =',F10.1)
110  FORMAT('  Bottom Delta-Z (m) =',F10.1)
!111  format(/,' max N lat=',f9.3,'  min S lat=',f9.3  &
!      ,/,' min W lon=',f9.3,' max E lon=',f9.3 )
113  format('      Outer N lat (deg) =       ',F8.4)
114  format('      Inner N lat (deg) =       ',F8.4)
115  format('      Outer W lon (deg) = ',F9.4,2X,F9.4,'  = Outer E lon (deg)')
116  format('      Inner W lon (deg) = ',F9.4,2X,F9.4,'  = Inner E lon (deg)')
117  format('      Inner S lat (deg) =       ',F8.4)
118  format('      Outer S lat (deg) =       ',F8.4,/)
112  format(80('-'),//)

if(print_msg) write(6,112)

return
END SUBROUTINE gridloc_prt

!##############################################################################
Subroutine refs3d (n1,n2,n3,pi0,dn0,dn0u,dn0v,th0,topt,rtgt)

use mem_grid
use ref_sounding
use rconstants
use node_mod

implicit none

integer :: n1,n2,n3
real :: pi0(n1,n2,n3),dn0(n1,n2,n3),dn0u(n1,n2,n3)  &
   ,dn0v(n1,n2,n3),th0(n1,n2,n3),topt(n2,n3),rtgt(n2,n3)

integer :: i,j,k
real :: c1,c2,c3
real, dimension(:), allocatable :: temp_vec

! +---------------------------------------------------------------------
! _    This routine initializes the 3-D reference state arrays from the
!        1-D reference state.
! +---------------------------------------------------------------------

allocate(temp_vec(n1))

do j=1,n3
   do i=1,n2
 
      do k = 1,n1
         temp_vec(k) = zt(k) * rtgt(i,j) + topt(i,j)
      enddo
      CALL htint (nzp,pi01dn(1,ngrid),zt,nzp,pi0(1,i,j),temp_vec)
      CALL htint (nzp,th01dn(1,ngrid),zt,nzp,th0(1,i,j),temp_vec)
      c1 = g * 2. * (1. - topt(i,j) / ztop)

      c2 = 1. - cpor
      c3 = cp ** c2
      do k = n1-1,1,-1
         pi0(k,i,j) = pi0(k+1,i,j)  &
                    + c1 / ((th0(k,i,j) + th0(k+1,i,j)) * dzm(k))
       enddo

      do k = 1,n1
         dn0(k,i,j) = (c3 * p00) / (rgas * th0(k,i,j) * pi0(k,i,j) ** c2)
      enddo

   enddo
enddo

CALL fill_dn0uv (n1,n2,n3,dn0,dn0u,dn0v)

deallocate(temp_vec)

return
END SUBROUTINE refs3d

!##############################################################################
Subroutine fill_dn0uv (n1,n2,n3,dn0,dn0u,dn0v)

implicit none

integer :: n1,n2,n3
real :: dn0(n1,n2,n3),dn0u(n1,n2,n3),dn0v(n1,n2,n3)
integer :: i,j,k,i1,j1

do j = 1,n3
   j1 = min(j+1,n3)
   do i = 1,n2
      i1 = min(i+1,n2)
      do k = 1,n1
         dn0u(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i1,j))
         dn0v(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i,j1))
      enddo
   enddo
enddo

return
END SUBROUTINE fill_dn0uv
