!##############################################################################
! update_ctile()
!
! This routine will make a node_send_ctile(), node_get_ctile() pair call.
!
Subroutine update_ctile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile &
                        ,ic0,jc0,ifm,icm,dimtype)

  implicit none

  integer :: nzc, nxc, nyc, nzct, nxct, nyct, ifm, icm, ic0, jc0, dimtype
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzct,nxct,nyct) :: varc_tile

  CALL node_send_ctile (nzc,nxc,nyc,varc,ic0,jc0,ifm,icm,dimtype)
  CALL node_get_ctile (nzct,nxct,nyct,varc_tile,ifm,dimtype)

  return
END SUBROUTINE update_ctile

!##############################################################################
! node_send_ctile()
!
! This routine, along with node_get_ctile, will form a coarse grid tile, which will
! be used for interpolation to the fine grid, and exchange the tile pieces among
! nodes so that each fine grid subdomain will end up with their associated coarse
! grid tile.
!
! The argument 'dimtype' describes what the actual variable organization looks like,
! while varc is always treated as a 3D variable (nzc,nxc,nyc). The caller of node_send_ctile()
! and node_get_ctile() needs to send in nzc,nxc,nyc,ic0,jc0 according to below:
!
!                 actual         argument settings
!  dimtype     organization    nzc   nxc   nyc   ic0  jc0
!
!    1              var(nzc)    nzc     1     1      0     0
!    2          var(nxc,nyc)      1   nxc   nyc    ic0   jc0
!    3      var(nzc,nxc,nyc)    nzc   nxc   nyc    ic0   jc0
!
Subroutine node_send_ctile (nzc,nxc,nyc,var,ic0,jc0,ifm,icm,dimtype)

  use node_mod

  implicit none

  integer :: nzc,nxc,nyc,ic0,jc0,ifm,icm,dimtype
  real, dimension(nzc,nxc,nyc) :: var

  integer :: nm,i1,i2,j1,j2,k1,k2,ng,itype,npts,ctile_tag
  real, allocatable :: buffnest(:)

  itype=5
  ctile_tag = CTILE_TAG_BASE + (icm * 1000) + (ifm * 100)

  !______________________
  !
  !   First, before we send anything, let's post the receives.

  do nm=1,nmachs
    irecv_req(nm)=0
    if (iget_paths(itype,ifm,nm).ne.not_a_node) then
      CALL par_get_noblock (node_buffs(nm)%lbc_recv_buff(1)  &
          ,node_buffs(nm)%nrecv ,ctile_tag,machnum(nm),irecv_req(nm) )
    endif
  enddo

  ! Send coarse grid points necessary for fine grid boundary interpolation
  !   to fine grid nodes. Note that even though coarse grid points are sent,
  !   ipaths is referenced by the fine grid, since all nests only have one
  !   parent, not vice versa.

  do nm=1,nmachs

    isend_req(nm)=0

    if(ipaths(5,itype,ifm,nm).ne. not_a_node) then

      i1=ipaths(1,itype,ifm,nm)
      i2=ipaths(2,itype,ifm,nm)
      j1=ipaths(3,itype,ifm,nm)
      j2=ipaths(4,itype,ifm,nm)
      k1=1
      k2=nzc

      select case (dimtype)
         case (1)
           ! var(nzc,1,1)
           ! force i1,i2,j1,j2 to 1 while creating the buffer
           npts = (k2-k1+1)
           allocate(buffnest(npts))
           CALL mknest_buff (var,buffnest,nzc,nxc,nyc,ic0,jc0,1,1,1,1,k1,k2)
         case (2)
           ! var(1,nxc,nyc)
           ! force k1,k2 to 1 while creating the buffer
           npts = (i2-i1+1)*(j2-j1+1)
           allocate(buffnest(npts))
           CALL mknest_buff (var,buffnest,nzc,nxc,nyc,ic0,jc0,i1,i2,j1,j2,1,1)
         case (3)
           ! var(nzc,nxc,nyc)
           npts = (k2-k1+1)*(i2-i1+1)*(j2-j1+1)
           allocate(buffnest(npts))
           CALL mknest_buff (var,buffnest,nzc,nxc,nyc,ic0,jc0,i1,i2,j1,j2,k1,k2)
      endselect

      ! Fill up MPI buffer and send
      CALL par_init_put (node_buffs(nm)%lbc_send_buff(1)  &
                       ,node_buffs(nm)%nsend )

      CALL par_put_int (i1,1)
      CALL par_put_int (i2,1)
      CALL par_put_int (j1,1)
      CALL par_put_int (j2,1)
      CALL par_put_int (k1,1)
      CALL par_put_int (k2,1)
      CALL par_put_int (npts,1)

      CALL par_put_float (buffnest,npts)

      CALL par_send_noblock (machnum(ipaths(5,itype,ifm,nm)),ctile_tag,isend_req(nm))

      deallocate(buffnest)

    endif

  enddo

  return

END SUBROUTINE node_send_ctile

!*************************************************************************************
! node_get_ctile()
!
! This routine, along with node_send_ctile, will form the coarse grid tiles to be
! used for interpolation to the fine grid.
!
Subroutine node_get_ctile (nzct, nxct, nyct, varc_tile, ifm, dimtype)
  
  use node_mod
  
  implicit none
  
  integer :: nzct, nxct, nyct, ifm, dimtype
  real, dimension(nzct,nxct,nyct) :: varc_tile
  
  integer, dimension(maxmach) :: i1c,i2c,j1c,j2c,k1c,k2c,nodim,iptc
  integer :: itype,nm,ibytes,msgid,ihostnum,iptr,nwords,nv,ict0,jct0
  real, allocatable :: buffnest(:)
  integer :: nbuff_size

  itype=5
  
  ! i[12]c, j[12]c, k[12]c contain dimensions for 3D vars, which will not
  ! work for 2D and 1D vars. nodim is an array that holds '1's which can
  ! be substituted for i[12]c, or j[12]c, or k[12] for non 3D vars.
  do nm = 1, nmachs
    nodim(nm) = 1
  enddo
  
  !_____________________________________________________________________
  !
  !  First, let's make sure our sends are all finished and de-allocated
  
  do nm=1,nmachs
     if(ipaths(5,itype,ifm,nm).ne.not_a_node)then
        CALL par_wait (isend_req(nm),ibytes,msgid,ihostnum)
     endif
  enddo
  !_____________________________________________________________________
  !
  !  Now, let's wait on our receives
  
  do nm=1,nmachs
     if(iget_paths(itype,ifm,nm).ne.not_a_node)then
        CALL par_wait (irecv_req(nm),ibytes,msgid,ihostnum)
     endif
  enddo
  !_____________________________________________________________________
  !
  
  ! Compute size of buffer needed and allocate if necessary
  ! Data from all incoming nodes will be placed into buffnest so the
  ! size of buffnest must be the sum of all the buffer sizes from the
  ! incoming nodes.
  
  nbuff_size = 0
  do nm = 1, nmachs
   if(iget_paths(itype,ifm,nm).ne.not_a_node)then
     nbuff_size = nbuff_size + node_buffs(nm)%nrecv
   endif
  enddo
  
  allocate (buffnest(nbuff_size))
  
  !     From the fine grid nodes, get the coarse grid buffers,
  !      interpolate the boundaries, and put them in the "b" array.
  
  iptr = 0
  ict0 = -1
  jct0 = -1
  do nm=1,nmachs
  
     if(iget_paths(itype,ifm,nm).ne.not_a_node) then
  
        CALL par_assoc_buff (node_buffs(nm)%lbc_recv_buff(1)  &
                           ,node_buffs(nm)%nrecv)
  
        CALL par_get_int (i1c(nm),1)
        CALL par_get_int (i2c(nm),1)
        CALL par_get_int (j1c(nm),1)
        CALL par_get_int (j2c(nm),1)
        CALL par_get_int (k1c(nm),1)
        CALL par_get_int (k2c(nm),1)
        CALL par_get_int (nwords,1)
  
        CALL par_get_float (buffnest(1+iptr),nwords)
        iptc(nm)=1+iptr
  
        iptr=iptr+nwords
  
        if ((ict0 .eq. -1) .or. (i1c(nm) .lt. ict0)) then
          ict0 = i1c(nm)
        endif
        if ((jct0 .eq. -1) .or. (j1c(nm) .lt. jct0)) then
          jct0 = j1c(nm)
        endif
     endif
  enddo
  ict0 = ict0 - 1
  jct0 = jct0 - 1
  
  ! At this point, buffnest contains the data for all of the tile pieces for
  ! the variable (vnum) from all of the nodes affecting this node. The data are
  ! organized sequentially in buffnest as follows:
  !
  !     <buffnest> := <node_block><node_block>...<node_block>
  !
  ! where <node_block> consiststs of coarse tile pieces for the variable from
  ! one of the nodes that affects this node.
  !
  ! Note that different node blocks can be different lengths.
  !
  ! iptc(nm) points to the beginning of the data from node nm in buffnest
  !

  select case (dimtype)
    case(1)
      CALL extract_var_from_buff (nzct,nxct,nyct,varc_tile,nbuff_size,buffnest,iptc, &
                                  nmachs,nodim,nodim,nodim,nodim,k1c,k2c,0,0,ifm,itype)
    case(2)
      CALL extract_var_from_buff (nzct,nxct,nyct,varc_tile,nbuff_size,buffnest,iptc, &
                                  nmachs,i1c,i2c,j1c,j2c,nodim,nodim,ict0,jct0,ifm,itype)
    case(3)
      CALL extract_var_from_buff (nzct,nxct,nyct,varc_tile,nbuff_size,buffnest,iptc, &
                                  nmachs,i1c,i2c,j1c,j2c,k1c,k2c,ict0,jct0,ifm,itype)
  endselect

  deallocate(buffnest)
  
  return

END SUBROUTINE node_get_ctile

!*************************************************************************************
Subroutine mknest_buff (ac,acs,m1,m2,m3,i0,j0,i1,i2,j1,j2,k1,k2)

implicit none

integer :: m1,m2,m3,i0,j0,i1,i2,j1,j2,k1,k2
real :: ac(m1,m2,m3),acs(0:k2-k1,0:i2-i1,0:j2-j1)
integer :: i,j,k

do j=j1,j2
   do i=i1,i2
      do k=k1,k2
         acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)
      enddo
   enddo
enddo

return
END SUBROUTINE mknest_buff

!*************************************************************************************
! extract_var_from_buffnest()
!
! This routine will piece together the tiles for a single variable by collecting
! and assimilating the corresponding var blocks from each node block.
!
Subroutine extract_var_from_buff (nzc,nxc,nyc,ctile,nbuff,buffnest,iptc, &
                                  nmachs,i1c,i2c,j1c,j2c,k1c,k2c,ict0,jct0,ifm,itype)

  use node_mod, only:IGET_PATHS,NOT_A_NODE

  implicit none

  integer :: nzc, nxc, nyc, nbuff, nmachs, ict0, jct0, ifm, itype
  real, dimension(nzc,nxc,nyc) :: ctile
  real, dimension(nbuff) :: buffnest
  integer, dimension(nmachs) :: iptc, i1c, i2c, j1c, j2c, k1c, k2c

  integer :: nm, numx, numy, numz, mtp

  do nm=1,nmachs
    if(iget_paths(itype,ifm,nm).ne.not_a_node) then
      numz=k2c(nm)-k1c(nm)+1
      numx=i2c(nm)-i1c(nm)+1
      numy=j2c(nm)-j1c(nm)+1
      mtp=numz*numx*numy
  
      ! Copy the var block data into the place where they fit within the ctile
      CALL unmkbuff (ctile,buffnest(iptc(nm)),nzc,nxc,nyc,numz,numx,numy,  &
             i1c(nm)-ict0,i2c(nm)-ict0,j1c(nm)-jct0,j2c(nm)-jct0,k1c(nm),k2c(nm))
  
    endif
  enddo
  
  return
END SUBROUTINE extract_var_from_buff

!##############################################################################
Subroutine unmkbuff (ac,buff,max1,max2,max3,m1,m2,m3,i1,i2,j1,j2,k1,k2)

implicit none

integer :: max1,max2,max3,m1,m2,m3,i1,i2,j1,j2,k1,k2
real :: ac(max1,max2,max3),buff(0:m1-1,0:m2-1,0:m3-1)
integer :: i,j,k

do j=j1,j2
   do i=i1,i2
      do k=k1,k2
         ac(k,i,j)=buff(k-k1,i-i1,j-j1)
      enddo
   enddo
enddo

return
END SUBROUTINE unmkbuff

