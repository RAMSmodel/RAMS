!############################################################################
! update_lbc_vgroup()
!
! node_sendlbc_vgroup() and node_getlbc_vgroup() always have to be called in pairs.
! This routine provides a convenient utility to do that.
!
Subroutine update_lbc_vgroup (ngr,isflag)
  implicit none

  integer :: ngr, isflag

  CALL node_sendlbc_vgroup (ngr,isflag)
  CALL node_getlbc_vgroup (ngr,isflag)

  return
END SUBROUTINE update_lbc_vgroup

!############################################################################
! node_sendlbc_vgroup()
!
! This routine along with node_getlbc_vgroup() handles the updating of internal
! (internodal) overlap regions. These routines will run on a group of variables
! as opposed to a single variable. For single variable updating, use the
! update_lbc_var(), node_sendlbc_var() and node_getlbc_var() routines.
!
! The variables involved in these transfers are of different sizes, but can
! all be cast into a general 4D array. The idim_type element of the variable
! table structures identifies the size of that variable. So for a 4D variable
! that is organized as var(n1,n2,n3,n4), use the following sizes for
! each dimension:
!
!   idim_type   n1    n2    n3     n4                 note
!
!       2        1   mxp   myp        1           2D atmos variable
!       3      mzp   mxp   myp        1           3D atmos variable
!       4      nzg   mxp   myp   npatch           4D soil variable
!       5      nzs   mxp   myp   npatch           4D surface water variable
!       6        1   mxp   myp   npatch           3D leaf variable
!       7      mzp   mxp   myp      nkr           3D atmos variable, bin micro
!      10    nkppz   mxp   myp        1           3D KPP ocean model variable     
!
! Use the following sizes for calculating the size of the buffer that
! holds the overlap piece (tile). This buffer will be a 1D buffer
! with size = n1 * n2 * n3 * n4
!
!   idim_type   n1    n2    n3     n4                 note
!
!       2        1   nxt   nyt        1           2D atmos variable
!       3      mzp   nxt   nyt        1           3D atmos variable
!       4      nzg   nxt   nyt   npatch           4D soil variable
!       5      nzs   nxt   nyt   npatch           4D surface water variable
!       6        1   nxt   nyt   npatch           3D leaf variable
!       7      mzp   nxt   nyt      nkr           3D atmos variable, bin micro
!      10    nkppz   nxt   nyt        1           3D KPP ocean model variable
!
! where nxt and nyt are the horizontal sizes of the tile.
!
Subroutine node_sendlbc_vgroup (ngr,isflag)

use mem_grid
use node_mod
use var_tables
use mem_basic
use mem_varinit
use micro_prm, only:nkr
use kpp_parameters, only:nkppz

implicit none

integer :: isflag,ngr
integer :: itype,nm,i1,i2,j1,j2,nv,mtp2,mtp3,mtp4,mtp5,mtp6,mtp7,mtp10
real, pointer :: scalarp
real, dimension(:), allocatable :: lbc_buff

itype=1
!______________________
!
!   First, before we send anything, let's post the receives. Also, make sure
!     any pending sends are complete.

do nm=1,nmachs
   if (iget_paths(itype,ngr,nm).ne.not_a_node) then
      CALL par_get_noblock (node_buffs(nm)%lbc_recv_buff(1)  &
          ,node_buffs(nm)%nrecv ,LBOUND_TAG_BASE+ngr,machnum(nm),irecv_req(nm) )
   endif
enddo


!______________________
!
!   Now we can actually go on to sending the stuff

do nm=1,nmachs

   if(ipaths(5,itype,ngr,nm).ne.not_a_node) then

      i1=ipaths(1,itype,ngr,nm)
      i2=ipaths(2,itype,ngr,nm)
      j1=ipaths(3,itype,ngr,nm)
      j2=ipaths(4,itype,ngr,nm)

      CALL par_init_put (node_buffs(nm)%lbc_send_buff(1)  &
                       ,node_buffs(nm)%nsend )
      CALL par_put_int (i1,1)
      CALL par_put_int (i2,1)
      CALL par_put_int (j1,1)
      CALL par_put_int (j2,1)
      CALL par_put_int (my_rams_num,1)

      mtp2  = (i2-i1+1)*(j2-j1+1)
      mtp3  = nnzp(ngr) * mtp2
      mtp4  = nzg * mtp2 * npatch
      mtp5  = nzs * mtp2 * npatch
      mtp6  = mtp2 * npatch
      mtp7  = nkr * mtp3
      mtp10 = nkppz * mtp2

      select case (isflag)
        case (LBC_ALL_VARS)
          do nv = 1,num_var(ngr)
             !This applies for impt1=1 which is for advected quantities
             if ( vtab_r(nv,ngr)%impt1 == 1) then
                if ( vtab_r(nv,ngr)%idim_type == 2) then
                   allocate(lbc_buff(mtp2))
                   CALL mklbcbuff (1,mxp,myp,1,mtp2,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp2)
                elseif ( vtab_r(nv,ngr)%idim_type == 3) then
                   allocate(lbc_buff(mtp3))
                   CALL mklbcbuff (mzp,mxp,myp,1,mtp3,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp3)
                elseif ( vtab_r(nv,ngr)%idim_type == 7) then
                   allocate(lbc_buff(mtp7))
                   CALL mklbcbuff (mzp,mxp,myp,nkr,mtp7,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp7)
                endif
                deallocate(lbc_buff)
             endif
          enddo
        case (LBC_ALL_INIT_VARS)
          do nv = 1,num_var(ngr)
             !This applies to all variables that need to be initialized
             if ( vtab_r(nv,ngr)%impti == 1) then

                if ( vtab_r(nv,ngr)%idim_type == 2) then
                   allocate(lbc_buff(mtp2))
                   CALL mklbcbuff (1,mxp,myp,1,mtp2,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp2)
                elseif ( vtab_r(nv,ngr)%idim_type == 3) then
                   allocate(lbc_buff(mtp3))
                   CALL mklbcbuff (mzp,mxp,myp,1,mtp3,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp3)
                elseif ( vtab_r(nv,ngr)%idim_type == 4) then
                   allocate(lbc_buff(mtp4))
                   CALL mklbcbuff (nzg,mxp,myp,npatch,mtp4,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp4)
                elseif ( vtab_r(nv,ngr)%idim_type == 5) then
                   allocate(lbc_buff(mtp5))
                   CALL mklbcbuff (nzs,mxp,myp,npatch,mtp5,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp5)
                elseif ( vtab_r(nv,ngr)%idim_type == 6) then
                   allocate(lbc_buff(mtp6))
                   CALL mklbcbuff (1,mxp,myp,npatch,mtp6,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp6)
                elseif ( vtab_r(nv,ngr)%idim_type == 7) then
                   allocate(lbc_buff(mtp7))
                   CALL mklbcbuff (mzp,mxp,myp,nkr,mtp7,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp7)
                elseif ( vtab_r(nv,ngr)%idim_type == 10) then
                   allocate(lbc_buff(mtp10))
                   CALL mklbcbuff (nkppz,mxp,myp,1,mtp10,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                   CALL par_put_float (lbc_buff,mtp10)

                endif
                deallocate(lbc_buff)

             endif
          enddo
        case (LBC_ALL_SCALARS)
            ! scalar vars (all scalars are 3D)
            do nv = 1, num_scalar(ngr)
              scalarp => scalar_tab(nv,ngr)%var_p
              allocate(lbc_buff(mtp3))
              CALL mklbcbuff (mzp,mxp,myp,1,mtp3,scalarp  &
                  ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
              CALL par_put_float (lbc_buff,mtp3)
              deallocate(lbc_buff)
            enddo
      endselect

      CALL par_send_noblock (machnum(ipaths(5,itype,ngr,nm))  &
           ,LBOUND_TAG_BASE+ngr,isend_req(nm))

   endif

enddo

return
END SUBROUTINE node_sendlbc_vgroup

!############################################################################
! node_getlbc_vgroup()
!
! This routine along with node_getlbc_vgroup() handles the updating of internal
! (internodal) overlap regions.
!
Subroutine node_getlbc_vgroup (ngr,isflag)

use mem_grid
use node_mod
use var_tables
use mem_basic
use mem_varinit
use micro_prm, only:nkr
use kpp_parameters, only:nkppz

implicit none

integer :: isflag,ngr
integer :: itype,nm,ibytes,msgid,ihostnum,i1,i2,j1,j2  &
          ,nv,node_src,mtp2,mtp3,mtp4,mtp5,mtp6,mtp7,mtp10
real, pointer :: scalarp
real, dimension(:), allocatable :: lbc_buff

itype=1

!_____________________________________________________________________
!
!  First, let's make sure our sends are all finished and de-allocated

do nm=1,nmachs
   if(ipaths(5,itype,ngr,nm).ne.not_a_node) then
      CALL par_wait (isend_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  Now, let's wait on our receives

do nm=1,nmachs
   if (iget_paths(itype,ngr,nm).ne.not_a_node) then
      CALL par_wait (irecv_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  We got all our stuff. Now unpack it into appropriate space.

do nm=1,nmachs

   if (iget_paths(itype,ngr,nm).ne.not_a_node) then


      CALL par_assoc_buff (node_buffs(nm)%lbc_recv_buff(1)  &
                         ,node_buffs(nm)%nrecv) 


      CALL par_get_int (i1,1)
      CALL par_get_int (i2,1)
      CALL par_get_int (j1,1)
      CALL par_get_int (j2,1)
      CALL par_get_int (node_src,1)

      mtp2  = (i2-i1+1)*(j2-j1+1)
      mtp3  = nnzp(ngr) * mtp2
      mtp4  = nzg * mtp2 * npatch
      mtp5  = nzs * mtp2 * npatch
      mtp6  = mtp2 * npatch
      mtp7  = nkr * mtp3
      mtp10 = nkppz * mtp2

      select case (isflag)
        case (LBC_ALL_VARS)
          do nv = 1,num_var(ngr)
             if ( vtab_r(nv,ngr)%impt1 == 1) then
                if ( vtab_r(nv,ngr)%idim_type == 2) then
                   allocate(lbc_buff(mtp2))
                   CALL par_get_float (lbc_buff,mtp2)
                   CALL exlbcbuff (1,mxp,myp,1,mtp2,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 3) then
                   allocate(lbc_buff(mtp3))
                   CALL par_get_float (lbc_buff,mtp3)
                   CALL exlbcbuff (mzp,mxp,myp,1,mtp3,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 7) then
                   allocate(lbc_buff(mtp7))
                   CALL par_get_float (lbc_buff,mtp7)
                   CALL exlbcbuff (mzp,mxp,myp,nkr,mtp7,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                endif
                deallocate(lbc_buff)
             endif
          enddo
        case (LBC_ALL_INIT_VARS)
          do nv = 1,num_var(ngr)
             if ( vtab_r(nv,ngr)%impti == 1) then
    
                if ( vtab_r(nv,ngr)%idim_type == 2) then
                   allocate(lbc_buff(mtp2))
                   CALL par_get_float (lbc_buff,mtp2)
                   CALL exlbcbuff (1,mxp,myp,1,mtp2,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 3) then
                   allocate(lbc_buff(mtp3))
                   CALL par_get_float (lbc_buff,mtp3)
                   CALL exlbcbuff (mzp,mxp,myp,1,mtp3,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 4) then
                   allocate(lbc_buff(mtp4))
                   CALL par_get_float (lbc_buff,mtp4)
                   CALL exlbcbuff (nzg,mxp,myp,npatch,mtp4,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 5) then
                   allocate(lbc_buff(mtp5))
                   CALL par_get_float (lbc_buff,mtp5)
                   CALL exlbcbuff (nzs,mxp,myp,npatch,mtp5,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 6) then
                   allocate(lbc_buff(mtp6))
                   CALL par_get_float (lbc_buff,mtp6)
                   CALL exlbcbuff (1,mxp,myp,npatch,mtp6,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 7) then
                   allocate(lbc_buff(mtp7))
                   CALL par_get_float (lbc_buff,mtp7)
                   CALL exlbcbuff (mzp,mxp,myp,nkr,mtp7,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                elseif ( vtab_r(nv,ngr)%idim_type == 10) then
                   allocate(lbc_buff(mtp10))
                   CALL par_get_float (lbc_buff,mtp10)
                   CALL exlbcbuff (nkppz,mxp,myp,1,mtp10,vtab_r(nv,ngr)%var_p  &
                       ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
                endif
                deallocate(lbc_buff)

             endif
          enddo
        case (LBC_ALL_SCALARS)
            ! scalar vars (all scalars are 3D)
            do nv = 1, num_scalar(ngr)
              allocate(lbc_buff(mtp3))
              scalarp => scalar_tab(nv,ngr)%var_p
              CALL par_get_float (lbc_buff,mtp3)
              CALL exlbcbuff (mzp,mxp,myp,1,mtp3,scalarp  &
                  ,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
              deallocate(lbc_buff)
            enddo
      endselect
   endif
enddo

return
END SUBROUTINE node_getlbc_vgroup

!############################################################################
! update_lbc_var()
!
! node_sendlbc_var() and node_getlbc_var() always have to be called in pairs.
! This routine provides a convenient utility to do that.
!
Subroutine update_lbc_var (nz,nx,ny,ngr,var,dimtype)
  implicit none

  integer :: nz, nx, ny, ngr, dimtype
  real, dimension(nz,nx,ny) :: var

  CALL node_sendlbc_var (nz,nx,ny,ngr,var,dimtype)
  CALL node_getlbc_var (nz,nx,ny,ngr,var,dimtype)

  return
END SUBROUTINE update_lbc_var

!############################################################################
! node_sendlbc_var()
!
! This routine along with node_getlbc_var() handles the updating of internal
! (internodal) overlap regions. These routines will handle a single variable.
!
Subroutine node_sendlbc_var (nz,nx,ny,ngr,var,dimtype)

use node_mod
use var_tables
use mem_basic
use mem_varinit

implicit none

integer :: nz, nx, ny, ngr, dimtype
real, dimension(nz,nx,ny) :: var

integer :: itype,nm,i1,i2,j1,j2,npts
real, dimension(:), allocatable :: lbc_buff

itype=1
!______________________
!
!   First, before we send anything, let's post the receives. Also, make sure
!     any pending sends are complete.

do nm=1,nmachs
   if (iget_paths(itype,ngr,nm).ne.not_a_node) then
      CALL par_get_noblock (node_buffs(nm)%lbc_recv_buff(1)  &
          ,node_buffs(nm)%nrecv ,LBOUND_TAG_BASE+ngr,machnum(nm),irecv_req(nm) )
   endif
enddo


!______________________
!
!   Now we can actually go on to sending the stuff

do nm=1,nmachs

   if(ipaths(5,itype,ngr,nm).ne.not_a_node) then

      i1=ipaths(1,itype,ngr,nm)
      i2=ipaths(2,itype,ngr,nm)
      j1=ipaths(3,itype,ngr,nm)
      j2=ipaths(4,itype,ngr,nm)

      CALL par_init_put (node_buffs(nm)%lbc_send_buff(1)  &
                       ,node_buffs(nm)%nsend )
      CALL par_put_int (i1,1)
      CALL par_put_int (i2,1)
      CALL par_put_int (j1,1)
      CALL par_put_int (j2,1)
      CALL par_put_int (my_rams_num,1)

      select case (dimtype)
        case (2)
          npts = ((i2-i1)+1)*((j2-j1)+1)
        case (3)
          npts = nz*((i2-i1)+1)*((j2-j1)+1)
      endselect

      allocate(lbc_buff(npts))
      CALL mklbcbuff (nz,nx,ny,1,npts,var,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
      CALL par_put_float (lbc_buff,npts)
      deallocate(lbc_buff)

      CALL par_send_noblock (machnum(ipaths(5,itype,ngr,nm))  &
           ,LBOUND_TAG_BASE+ngr,isend_req(nm))

   endif

enddo

return
END SUBROUTINE node_sendlbc_var

!############################################################################
! node_getlbc_var()
!
! This routine along with node_getlbc() handles the updating of internal
! (internodal) overlap regions. These routines handle a single variable.
!
Subroutine node_getlbc_var (nz,nx,ny,ngr,var,dimtype)

use node_mod
use var_tables
use mem_basic
use mem_varinit

implicit none

integer :: nz, nx, ny, ngr, dimtype
real, dimension(nz,nx,ny) :: var

integer :: itype,nm,ibytes,msgid,ihostnum,i1,i2,j1,j2,node_src,npts
real, dimension(:), allocatable :: lbc_buff

itype=1

!_____________________________________________________________________
!
!  First, let's make sure our sends are all finished and de-allocated

do nm=1,nmachs
   if(ipaths(5,itype,ngr,nm).ne.not_a_node) then
      CALL par_wait (isend_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  Now, let's wait on our receives

do nm=1,nmachs
   if (iget_paths(itype,ngr,nm).ne.not_a_node) then
      CALL par_wait (irecv_req(nm),ibytes,msgid,ihostnum)
   endif
enddo
!_____________________________________________________________________
!
!  We got all our stuff. Now unpack it into appropriate space.

do nm=1,nmachs

   if (iget_paths(itype,ngr,nm).ne.not_a_node) then


      CALL par_assoc_buff (node_buffs(nm)%lbc_recv_buff(1)  &
                         ,node_buffs(nm)%nrecv) 


      CALL par_get_int (i1,1)
      CALL par_get_int (i2,1)
      CALL par_get_int (j1,1)
      CALL par_get_int (j2,1)
      CALL par_get_int (node_src,1)

      select case (dimtype)
        case (2)
          npts = ((i2-i1)+1)*((j2-j1)+1)
        case (3)
          npts = nz*((i2-i1)+1)*((j2-j1)+1)
      endselect

      allocate(lbc_buff(npts))
      CALL par_get_float (lbc_buff,npts)
      CALL exlbcbuff (nz,nx,ny,1,npts,var,lbc_buff,i1-i0,i2-i0,j1-j0,j2-j0)
      deallocate(lbc_buff)
   endif
enddo

return
END SUBROUTINE node_getlbc_var

!##############################################################################
Subroutine mklbcbuff (n1,n2,n3,n4,nb,a,b,il,ir,jb,jt)

implicit none

integer :: n1,n2,n3,n4,nb,il,ir,jb,jt
real :: a(n1,n2,n3,n4),b(nb)
integer :: i,j,k,l,ind

ind=0
do l=1,n4
  do j=jb,jt
    do i=il,ir
      do k=1,n1
        ind=ind+1
        b(ind)=a(k,i,j,l)
      enddo
    enddo
  enddo
enddo

return
END SUBROUTINE mklbcbuff

!##############################################################################
Subroutine exlbcbuff (n1,n2,n3,n4,nb,a,b,il,ir,jb,jt)

implicit none

integer :: n1,n2,n3,n4,nb,il,ir,jb,jt
real :: a(n1,n2,n3,n4),b(nb)
integer :: i,j,k,l,ind

ind=0
do l=1,n4
  do j=jb,jt
    do i=il,ir
      do k=1,n1
        ind=ind+1
        a(k,i,j,l)=b(ind)
      enddo
    enddo
  enddo
enddo

return
END SUBROUTINE exlbcbuff
