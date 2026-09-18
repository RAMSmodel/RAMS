!##############################################################################
Subroutine rams_node ()

  use mem_grid
  use node_mod
  use micphys, only:idiffperts,iconv
  use kpp_parameters, only:IKPP
  use ref_sounding, only:iugforce

  implicit none

  real :: wtime_ref,wtime_beg,wtime_end,ctime_beg,ctime_end,begtime,timeh
  real, external :: walltime
  character (len=200) :: timestep_msg
  integer :: npass,igrid,icm,ifm,nfeed

  ! Record a reference time for the walltime routine.
  wtime_ref = walltime(0.0)

  ! Initialize counters
  istp = 0

  ! Run the time integration.
  do while (time < timmax)

    ! Record start cpu and wall times.
    CALL timing (1,ctime_beg)
    wtime_beg = walltime(wtime_ref)

    ! Record beginning time for this integration step. This is used to
    ! offset time to the approriate value when working on nested grid
    ! integration steps.
    begtime = time

    istp = istp + 1

    ! Examine Courant numbers in case model needs to be stopped or
    CALL dtset ()

!----------------------------------------------------------------------------
!   Update base state winds if doing geostrophic wind forcing
!----------------------------------------------------------------------------
    if(iugforce==2)then
      CALL update_base_winds()
    endif
!----------------------------------------------------------------------------
!   Compute new "temporary base state" for THP and RTP for diffusing
!   perturbations from the mean state at current time.
!----------------------------------------------------------------------------
    if(idiffperts==1 .or. idiffperts==2)then
     do ifm = 1,ngrids
      CALL newgrid (ifm)
      CALL new_base_state (mzp,mxp,myp)
     enddo
    endif

!----------------------------------------------------------------------------
!                  Loop through all grids and advance a 'DTLONG' timestep.
!----------------------------------------------------------------------------
!                  Start the timestep schedule

    do npass = 1, nsubs
      !print*,'-->starting schedule :',mynum,npass,nt,nsubs

      isstp=isched(npass,3)
      igrid=isched(npass,1)
      CALL newgrid (igrid)

      CALL node_index ()

!----------------------------------------------------------------------------
!                  Advance this grid forward by the appropriate timestep.

      time=begtime + (isched(npass,5)-1) * dtlt

      !Call main timestep driver
      CALL timestep ()

      ngbegun(igrid)=1

!----------------------------------------------------------------------------
!                  Is it time to send the coarse domain points to the nested
!                  nodes to interpolate a nested grid's boundaries?

      if(isched(npass,2).ne.0) then
        ifm=isched(npass,2)
        CALL newgrid (ifm)
        isstp=isched(npass,3)
        icm=nxtnest(ifm)
        CALL interp_nest_bounds (icm,ifm)
      endif

    enddo ! do npass = 1, nsubs

! All nesting operations (interpolation) have been carried out.

!----------------------------------------------------------------------------
! Compute Courant numbers cflxy and cflz and do averaging.

    do ifm = 1,ngrids
      CALL newgrid (ifm)
      CALL cfl (mzp,mxp,myp,i0,j0,my_rams_num)
      !THETAM and RVM have not been updated after nesting feedback
      !This means that these variables are really a timestep
      !behind the instantaneous variables.
      !Calculate the mean/average each of the analysis variables over time
      CALL anlavg ()
    enddo

    ! get cfl max numbers for checking during the next iteration (by dtset())
    CALL set_cflmax_allnodes ()

    ! Make the following call to get the domain max vertical velocity for the
    ! convergence forcing code used to initiate idealized convection.
    if(iprntstmt>=1 .or. ICONV > 0) &
     CALL set_vertvelmax_allnodes ()

    ! Make the following call to get the domain max mixed layer depth if
    ! KPP mixed layer ocean model is running
    if(iprntstmt>=1 .and. IKPP > 0) &
     CALL set_hmixmax_allnodes ()

    ! Move to next integration step. Make sure you are stepping by the
    ! primary (grid 1) dtlong value.
    time = begtime + dtlongn(1)

    ! Report time for this integration step.
    CALL timing (2,ctime_end)
    wtime_end = walltime(wtime_ref)

    timeh = time / 3600.0
    write(timestep_msg, 900) istp, time, timeh
900 format(' Timestep-', i8, '   Sim time(sec)=', f10.1, '  (hr)=', f9.2)
    CALL report_time_allnodes (trim(timestep_msg),ctime_end-ctime_beg &
                              ,wtime_end-wtime_beg)
    if(print_msg) then
      print*,'----------------------------------------------'
    endif

    ! Write analysis files. rams_output() will check to see if it's the
    ! proper time before issuing the file write.
    CALL rams_output ()

  enddo ! do while (time < timmax)

return
END SUBROUTINE rams_node

!##############################################################################
Subroutine init_fields ()

use mem_grid
use node_mod
use var_tables

implicit none

integer :: ng,nm,itype,i1,j1,i2,j2,memf,npvar,nv


!     Can we use existing memory for the nesting communication buffers?
!       If not, allocate new buffers or compute buffer sizes.

!       Check feedback buffer.

itype=6
nbuff_feed=0
do ng=1,ngrids
   do nm=1,nmachs
      i1=ipaths(1,itype,ng,nm)
      i2=ipaths(2,itype,ng,nm)
      j1=ipaths(3,itype,ng,nm)
      j2=ipaths(4,itype,ng,nm)
      memf=(i2-i1+1)*(j2-j1+1)*(nnzp(ng))  &
           *(NUM_NEST_VARS+num_scalar(ng))
      nbuff_feed=max(nbuff_feed,memf)
   enddo
enddo

!____________________________________________
!
!    Allocate long time step send and receive buffers
do nm=1,nmachs
   ! send buffers
   if (allocated(node_buffs(nm)%lbc_send_buff) ) &
       deallocate(node_buffs(nm)%lbc_send_buff)
   if (node_buffs(nm)%nsend > 0) &
       allocate(node_buffs(nm)%lbc_send_buff(node_buffs(nm)%nsend))

   ! receive buffers
   if (allocated(node_buffs(nm)%lbc_recv_buff) ) &
       deallocate(node_buffs(nm)%lbc_recv_buff)
   if (node_buffs(nm)%nrecv > 0) &
       allocate(node_buffs(nm)%lbc_recv_buff(node_buffs(nm)%nrecv))
   if(iprntstmt>=2)print'(a,6i10)','LBC node alloc:',my_rams_num,nm &
       ,node_buffs(nm)%nsend, node_buffs(nm)%nrecv
enddo

! If using cyclic boundary conditions, initialize parallel communication
! for them
!  Find number of lbc variables to be communicated.
   npvar=0
   do nv = 1,num_var(1)
      if(vtab_r(nv,1)%impt1 == 1 ) then
         npvar=npvar+1 
      endif
   enddo

if (ibnd .eq. 2 .or. jbnd .eq. 2) then
   CALL node_cycinit (nnzp(1),npvar,nmachs,my_rams_num)
endif

return
END SUBROUTINE init_fields

!##############################################################################
Subroutine node_index ()

use node_mod

implicit none

ia_1=max(ia-1,1)
ia_2=max(ia-2,1)
ia_3=max(ia-3,1)
ia1=ia+1
ia2=ia+2
ia3=ia+3
iz_1=iz-1
iz_2=iz-2
iz_3=iz-3
iz1=min(iz+1,mxp)
iz2=min(iz+2,mxp)
iz3=min(iz+3,mxp)

izu=iz
if(iand(ibcon,2).ne.0) izu=iz-1

if(myp.gt.1) then
   ja_1=max(ja-1,1)
   ja_2=max(ja-2,1)
   ja_3=max(ja-3,1)
   ja1=ja+1
   ja2=ja+2
   ja3=ja+3
   jz_1=jz-1
   jz_2=jz-2
   jz_3=jz-3
   jz1=min(jz+1,myp)
   jz2=min(jz+2,myp)
   jz3=min(jz+3,myp)

   jzv=jz
   if(iand(ibcon,8).ne.0) jzv=jz-1

else
   if (nmachs .gt. 1) then
     print*,'Trying to do 2-dimensional run?'
     stop 'no parallel 2d'
   endif
   ja=1
   jz=1
   ja_1=1
   ja_2=1
   ja_3=1
   ja1=1
   ja2=1
   ja3=1
   jz_1=1
   jz_2=1
   jz_3=1
   jz1=1
   jz2=1
   jz3=1
   jzv=1
endif

return
END SUBROUTINE node_index
