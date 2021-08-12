!##############################################################################
! This routine will collect elapsed cpu and wall clock times from all of the
! nodes and report these to stdout. mainnum node will handle the collecting
! and reporting so that only one copy of the message is output. MPI_Gather()
! is used to do the collecting. MPI_Gather() will place the data into the 
! receive buffers in the order of the node numbering.
!##############################################################################
Subroutine report_time_allnodes (msg,cpu_etime,wall_etime)

use mem_grid
use node_mod

implicit none

character (len=*) :: msg
real :: cpu_etime,wall_etime,sum_cpu_time,sum_wall_time,avg_cpu_time &
        ,avg_wall_time
real, dimension(:), allocatable :: cpu_times,wall_times
integer :: i

allocate(cpu_times(nmachs),wall_times(nmachs))

  if (my_rams_num .eq. mainnum) then
    ! mainnum collects the data and prints out the messages.
    if (nmachs .gt. 1) then
      CALL par_gather_floats (cpu_etime,cpu_times,1,machnum(mainnum))
      CALL par_gather_floats (wall_etime,wall_times,1,machnum(mainnum))

      sum_cpu_time = 0.0
      sum_wall_time = 0.0
      do i = 1, nmachs
        sum_cpu_time = sum_cpu_time + cpu_times(i)
        sum_wall_time = sum_wall_time + wall_times(i)
      enddo
      avg_cpu_time = sum_cpu_time / float(nmachs)
      avg_wall_time = sum_wall_time / float(nmachs)
    else
      cpu_times(1) = cpu_etime
      wall_times(1) = wall_etime
      sum_cpu_time = cpu_etime
      sum_wall_time = wall_etime
      avg_cpu_time = cpu_etime
      avg_wall_time = wall_etime
    endif

!    print 200,msg,sum_cpu_time,avg_cpu_time,sum_wall_time,avg_wall_time
!200 format(a,'  CPU time- Sum:Avg(sec)=(',f10.2,':',f10.2,') &
!       , Wall time- Sum:Avg(sec)=(',f10.2,':',f10.2')')
    print 200,msg,avg_wall_time
200 format(a,'  Wall time(sec)=',f6.2)

    ! If printing is enabled, add in the individual node times
    if(iprntstmt>=2 .and. print_msg) then
      print 201,(cpu_times(i),i=1,nmachs)
      print 202,(wall_times(i),i=1,nmachs)
201   format(' Node---CPU(sec)=',2000f8.3)
202   format(' Node--wall(sec)=',2000f8.3)
    endif
  else
    if (nmachs .gt. 1) then
      ! Nodes other than mainnum send their data. Note that the machine
      ! number of mainnum is used instead of the node of this process.
      CALL par_gather_floats (cpu_etime,cpu_times,1,machnum(mainnum))
      CALL par_gather_floats (wall_etime,wall_times,1,machnum(mainnum))
    endif
  endif

  deallocate(cpu_times,wall_times)

return
END SUBROUTINE report_time_allnodes

!##############################################################################
! This routine will collect cflxy and cflz numbers from all of the nodes
! nodes to mainnum, determine the maximum of these and broadcast the result
! back to all of the nodes. MPI_Gather() is used to do the collecting.
! MPI_Gather() will place the data into the receive buffers in the order of
! the node numbering.
!
! Note that the cflxy_all and cflz_all arrays must be organized as (ngrids,nmachs)
! since the gather func will place the cflxy(1:ngrids) from node 1 in the first
! portion of the buffer, followed by clfxy(1:ngrids) from node 2, etc.

Subroutine set_cflmax_allnodes ()

use mem_grid
use node_mod

implicit none

real, dimension(:,:), allocatable :: cflxy_all, cflz_all
integer :: igrid,inode,nwords
real, dimension(:), allocatable :: buff

allocate(cflxy_all(ngrids,nmachs),cflz_all(ngrids,nmachs))

  ! Prep for broadcast
  if (nmachs .gt. 1) then
    nwords = 10 + (2 * ngrids)
    allocate (buff(nwords))
  endif

  ! Collect the cfl numbers from all the nodes
  if (my_rams_num .eq. mainnum) then
    ! mainnum collects the data and prints out the messages.
    if (nmachs .gt. 1) then
      CALL par_gather_floats (cflxy,cflxy_all,ngrids,machnum(mainnum))
      CALL par_gather_floats (cflz,cflz_all,ngrids,machnum(mainnum))
    else
      cflxy_all(1:ngrids,1) = cflxy(1:ngrids)
      cflz_all(1:ngrids,1) = cflz(1:ngrids)
    endif
  else
    if (nmachs .gt. 1) then
      ! Nodes other than mainnum send their data. Note that the machine
      ! number of mainnum is used instead of the node of this process.
      CALL par_gather_floats (cflxy,cflxy_all,ngrids,machnum(mainnum))
      CALL par_gather_floats (cflz,cflz_all,ngrids,machnum(mainnum))
    endif
  endif

  ! At this point, cflxy_all and cflz_all are filled up with the
  ! cfl numbers from all of the nodes (including the case where there
  ! is only one node.
  !
  ! If this is the mainnum node, calculate the cflxy_max and cflz_max values
  ! and broadcast to the other nodes.
  if (my_rams_num .eq. mainnum) then
    do igrid = 1, ngrids
      cflxy_max(igrid) = cflxy_all(igrid,1)
      cflz_max(igrid) = cflz_all(igrid,1)

      do inode = 2, nmachs
        if (cflxy_all(igrid,inode) .gt. cflxy_max(igrid)) then
          cflxy_max(igrid) = cflxy_all(igrid,inode)
        endif
        if (cflz_all(igrid,inode) .gt. cflz_max(igrid)) then
          cflz_max(igrid) = cflz_all(igrid,inode)
        endif
      enddo
    enddo

    if (nmachs .gt. 1) then
      ! broadcast the cflmax numbers to all the nodes
      CALL par_init_put (buff,nwords)
      CALL par_put_float (cflxy_max,ngrids)
      CALL par_put_float (cflz_max,ngrids)
      CALL par_broadcast (machnum(mainnum))
    endif
  else
    if (nmachs .gt. 1) then
      ! Receiving nodes for broadcast
      CALL par_init_recv_bcast (buff,nwords)
      CALL par_broadcast (machnum(mainnum))
      CALL par_get_float (cflxy_max,ngrids)
      CALL par_get_float (cflz_max,ngrids)
    endif
  endif

  deallocate(cflxy_all,cflz_all)

  if (nmachs .gt. 1) then
    deallocate(buff)
  endif

return
END SUBROUTINE set_cflmax_allnodes

!##############################################################################
! This routine will collect vertical velocity from all of the nodes
! nodes to mainnum, determine the maximum of these and broadcast the result
! back to all of the nodes. MPI_Gather() is used to do the collecting.
! MPI_Gather() will place the data into the receive buffers in the order of
! the node numbering.
!
! Note that the vertvel_all array must be organized as (ngrids,nmachs) since
! the gather func will place the vertvel(1:ngrids) from node 1 in the first
! portion of the buffer, followed by vertvel(1:ngrids) from node 2, etc.

Subroutine set_vertvelmax_allnodes ()

use mem_grid
use node_mod

implicit none

real, dimension(:,:), allocatable :: vertvel_all
integer :: igrid,inode,nwords
real, dimension(:), allocatable :: buff

  allocate(vertvel_all(ngrids,nmachs))
  CALL azero (ngrids*nmachs,vertvel_all)

  ! Prep for broadcast
  if (nmachs .gt. 1) then
    nwords = 10 + (1 * ngrids)
    allocate (buff(nwords))
  endif

  ! Collect the vertical velocity from all the nodes
  if (my_rams_num .eq. mainnum) then
    ! mainnum collects the data and prints out the messages.
    if (nmachs .gt. 1) then
      CALL par_gather_floats (vertvel,vertvel_all,ngrids,machnum(mainnum))
    else
      vertvel_all(1:ngrids,1) = vertvel(1:ngrids)
    endif
  else
    if (nmachs .gt. 1) then
      ! Nodes other than mainnum send their data. Note that the machine
      ! number of mainnum is used instead of the node of this process.
      CALL par_gather_floats (vertvel,vertvel_all,ngrids,machnum(mainnum))
    endif
  endif

  ! At this point, vertvel_all is filled up with the
  ! vertvel from all of the nodes (including the case where there
  ! is only one node.
  !
  ! If this is the mainnum node, calculate the vertvel_max values
  ! and broadcast to the other nodes.
  if (my_rams_num .eq. mainnum) then
    do igrid = 1,ngrids
      vertvel_max(igrid) = vertvel_all(igrid,1)

      do inode = 2,nmachs
        if (vertvel_all(igrid,inode) .gt. vertvel_max(igrid)) then
          vertvel_max(igrid) = vertvel_all(igrid,inode)
        endif
      enddo
    enddo

    if (nmachs .gt. 1) then
      ! broadcast the max vertical velocity to all the nodes
      CALL par_init_put (buff,nwords)
      CALL par_put_float (vertvel_max,ngrids)
      CALL par_broadcast (machnum(mainnum))
    endif
  else
    if (nmachs .gt. 1) then
      ! Receiving nodes for broadcast
      CALL par_init_recv_bcast (buff,nwords)
      CALL par_broadcast (machnum(mainnum))
      CALL par_get_float (vertvel_max,ngrids)
    endif
  endif

  deallocate(vertvel_all)

  if (nmachs .gt. 1) then
    deallocate(buff)
  endif

return
END SUBROUTINE set_vertvelmax_allnodes

!##############################################################################
! This routine will collect ocean mixed-layer depth (HMIX) from all the nodes
! nodes to mainnum, determine the maximum of these and broadcast the result
! back to all of the nodes. MPI_Gather() is used to do the collecting.
! MPI_Gather() will place the data into the receive buffers in the order of
! the node numbering.
!
! Note that the hmixdepth_all array must be organized as (ngrids,nmachs) since
! the gather func will place the hmixdepth(1:ngrids) from node 1 in the first
! portion of the buffer, followed by hmixdepth(1:ngrids) from node 2, etc.

Subroutine set_hmixmax_allnodes ()

use mem_grid
use node_mod
use kpp_parameters, only:hmixdepth,hmixdepth_max,ikpp

implicit none

real, dimension(:,:), allocatable :: hmixdepth_all
integer :: igrid,inode,nwords
real, dimension(:), allocatable :: buff

  allocate(hmixdepth_all(ngrids,nmachs))
  CALL azero (ngrids*nmachs,hmixdepth_all)

  ! Prep for broadcast
  if (nmachs .gt. 1) then
    nwords = 10 + (1 * ngrids)
    allocate (buff(nwords))
  endif

  ! Collect the vertical velocity from all the nodes
  if (my_rams_num .eq. mainnum) then
    ! mainnum collects the data and prints out the messages.
    if (nmachs .gt. 1) then
      CALL par_gather_floats (hmixdepth,hmixdepth_all,ngrids,machnum(mainnum))
    else
      hmixdepth_all(1:ngrids,1) = hmixdepth(1:ngrids)
    endif
  else
    if (nmachs .gt. 1) then
      ! Nodes other than mainnum send their data. Note that the machine
      ! number of mainnum is used instead of the node of this process.
      CALL par_gather_floats (hmixdepth,hmixdepth_all,ngrids,machnum(mainnum))
    endif
  endif

  ! At this point, hmixdepth_all is filled up with the
  ! hmixdepth from all of the nodes (including the case where there
  ! is only one node.
  !
  ! If this is the mainnum node, calculate the hmixdepth_max values
  ! and broadcast to the other nodes.
  if (my_rams_num .eq. mainnum) then
    do igrid = 1,ngrids
      hmixdepth_max(igrid) = hmixdepth_all(igrid,1)

      do inode = 2,nmachs
        if (hmixdepth_all(igrid,inode) .gt. hmixdepth_max(igrid)) then
          hmixdepth_max(igrid) = hmixdepth_all(igrid,inode)
        endif
      enddo
    enddo

    if (nmachs .gt. 1) then
      ! broadcast the max vertical velocity to all the nodes
      CALL par_init_put (buff,nwords)
      CALL par_put_float (hmixdepth_max,ngrids)
      CALL par_broadcast (machnum(mainnum))
    endif
  else
    if (nmachs .gt. 1) then
      ! Receiving nodes for broadcast
      CALL par_init_recv_bcast (buff,nwords)
      CALL par_broadcast (machnum(mainnum))
      CALL par_get_float (hmixdepth_max,ngrids)
    endif
  endif

  deallocate(hmixdepth_all)

  if (nmachs .gt. 1) then
    deallocate(buff)
  endif

  if(iprntstmt>=1 .and. print_msg) &
    print*,'Ngrid, Time, Domain Max HMIX: ',ngrid,time,hmixdepth_max(ngrid)

return
END SUBROUTINE set_hmixmax_allnodes

