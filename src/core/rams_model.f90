!##############################################################################
Subroutine rams_model (nl_fname)

use node_mod
use mem_all
use isan_coms

implicit none

  integer :: stopparallel

  character (len=*) :: nl_fname
  character(len=strl1) :: anamelh
  integer :: i,ifm,force1node
  real :: wtime_ref,wtime_beg,wtime_end,ctime_beg,ctime_end
  real, external :: walltime

  ! runtype name RAMS or REVU distinction
  ramsorrevu='RAMS'
  
  ! This variable is accessed by a number of routines and is intended
  ! to help de-clutter messages being printed to stdout. Be careful
  ! however - only use this in cases where all nodes would be printing
  ! the same thing. Don't use this when different nodes need to
  ! print out messages specific to each node.
  print_msg = (my_rams_num .eq. mainnum)

  ! record the start time for subsequent reporting
  wtime_ref = walltime(0.0)
  CALL timing (1,ctime_beg)
  wtime_beg = walltime(wtime_ref)

  ! High level flow:
  !
  !  1. Establish configuration
  !  2. Set up domain
  !  3. Decompose domain into sub-domains
  !  4. Set up sub-domains
  !  5. Initialize sub-domains
  !  6. Run the model
  !
  ! Steps 1, 2 and 3 are done by RAMS node 1 (mainnum), then node 1
  ! broadcasts the results to all other nodes. This is
  ! done to help efficiency (ie, read namelist file once instead
  ! of potentially thousands of times - once per node) and to
  ! guarantee that the variables describing the domain and
  ! decomposition of the domain match on all nodes.
  !
  ! Then steps 4, 5 and 6 are done on all nodes.

  !---------------- Configuration ---------------
  CALL establish_config (nl_fname)

  ! Set flags regarding history restart or initialization and set some
  ! initial condition flags with not a history restart.
  hrestart=0
  if(trim(runtype)=='HISTORY') then  
     CALL read_distribute_hheader (hfilin)
     initial=3
     if(ngrids <= ngridsh) hrestart=1
     if(ngrids  > ngridsh) hrestart=2
  else
     hrestart=0
     if(initial==3) hrestart=2
     initorig=0
     time=0.
     ngbegun(1:ngrids) = 0
  endif

  ! If this is a history initialization, history restart with added grids,
  ! or a recycle simulation, then we need to run this sequentially due to 
  ! history vs current grid domain mismatches that currently cannot be 
  ! handled in parallel.
  stopparallel=0
  if(trim(runtype)=='MAKESFC' .or. trim(runtype)=='MAKEVFILE' .or. &
     trim(runtype)=='MAKEHFILE') then
    stopparallel=1
  elseif(trim(runtype)=='INITIAL' .and. &
     (initial==3 .or. ((initial==1.or.initial==2).and.ipast_sfc==1))) then
    stopparallel=2
  elseif(trim(runtype)=='HISTORY' .and. hrestart==2) then
    stopparallel=3
  endif

  !Here we reset nmachs to 1 for sequential and return the other nodes.
  force1node=0
  if (stopparallel > 0)then
   if (nmachs .gt. 1 .and. my_rams_num .ne. mainnum) then
    return
   endif
   if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
    nmachs=1
    force1node=1
    print*,''
    print*,'**************************************************************'
    print*,'FORCING NUMBER OF PROCESSORS TO 1 FOR:'
    if(stopparallel==1) then
      print*,'MAKESFC, MAKEVFILE, or MAKEHFILE RUN'
    elseif(stopparallel==2) then
      print*,'HISTORY INITIALIZATION OR RECYCLING OF VARIABLES'
      print*,'ON POTENTIALLY DIFFERENT GRIDS'
    elseif(stopparallel==3) then
      print*,'HISTORY RESTART WITH ADDED GRIDS'
    endif
    print*,'**************************************************************'
    print*,''
   endif
  endif

  !---------------- Domain-wide Grid Setup -------------------------------
  CALL set_full_domain_grid ()

  !---------------- Domain Decomposition into Sub-Domains ---------------
  CALL decompose_full_domain ()

  ! At this point all nodes have enough information to proceed on their own.

  !---------------- Sub-domain Setup ---------------
  ! Fill in the description of the sub domain (mmxp, mmyp, mmzp, etc. plus node paths)
  CALL node_init ()

  ! Assign hydrometeor modes
  ! This needs to be called before rams_mem_alloc() since alloc_micro()
  ! uses the hydrometeor mode numbers.
  CALL jnmbinit ()

  ! Allocate variable memory for INITIAL and HISTORY starts. Make surface and
  ! varfiles have their own smaller memory allocation.
  ! This needs to be called before finishing node decomposition and grid setup.
  if (stopparallel/=1) CALL rams_mem_alloc ()
  if (stopparallel==1) CALL sfc_mem_alloc ()

  ! Set the buffer sizes for internodal communication
  if(force1node==0) CALL par_set_comm_buff_sizes ()

  !-------Finish off MAKESFC or MAKEVFILE or MAKEHFILE--------------------
  ! No more initialization is required to complete either a MAKESFC run
  ! or a MAKEVFILE run, so these can be completed at this time.
  !
  ! Any run (INITIAL, HISTORY, MAKESFC, MAKEVFILE, MAKEHFILE, ERROR) needs 
  ! a set of surface files in place so always call update_sfc_var_files() at 
  ! this point.

  ! Make sure that the correct surface files and varfiles are in place
  CALL update_sfc_var_files ()

  if(trim(runtype) .eq. 'MAKESFC') then
    if(print_msg) then
      print*, 'MAKESFC run complete'
    endif
    return
  endif

  if(trim(runtype) .eq. 'MAKEVFILE') then
    if(print_msg) then
       print*, ' ISAN complete'
    endif
    return
  endif

  if(trim(runtype) .eq. 'MAKEHFILE') then
    if(print_msg) then
       print*, ' History Nudging Varfiles Complete'
    endif
    return
  endif

  !------------------ Initialization ------------------------
  ! Allocate memory for communicating boundary data between nodes
  if (nmachs .gt. 1) then
    CALL init_fields ()
  endif

  ! main initialization driver
  CALL initlz ()
  if(force1node==1)then
    print*,'SINGLE PROCESSOR INITIALIZATION WAS FORCED.'
    print*,'INITIAL ANALYSIS FILE WAS OUTPUT FOR:'
    if(stopparallel==2) then
      print*,'HISTORY INITIALIZATION OR RECYCLING OF VARIABLES'
      print*,'ON POTENTIALLY DIFFERENT GRIDS'
    elseif(stopparallel==3) then
      print*,'HISTORY RESTART WITH ADDED GRIDS'
    endif
    print*,'PLEASE PROCEED BY RUNNING HISTORY RESTART ON THE NEWLY'
    print*,'OUTPUT ANALYSIS FILE. THIS CAN BE DONE IN PARALLEL MODE.'  
    return
  endif

  ! Compute Courant numbers cflxy and cflz.
  do ifm = 1, ngrids
    CALL newgrid (ifm)
    CALL cfl (mzp,mxp,myp,i0,j0,my_rams_num)
  enddo

  ! Initialize dtlongn, nndtrat, and nnacoust, and compute the timestep
  ! schedule for all grid operations.
  CALL dtset ()
  CALL modsched (isched,maxsched,ngrids,nxtnest,nndtrat,nsubs)

  ! report time for initialization
  CALL timing (2,ctime_end)
  wtime_end = walltime(wtime_ref)
  CALL report_time_allnodes ('++++++++ Model Initialization:' &
      ,ctime_end-ctime_beg, wtime_end-wtime_beg)
  if(print_msg) print*,''

  ! Exit if doing a zero time run or if doing and ERROR run
  if(time >= timmax) then
    return
  elseif (trim(runtype) == 'ERROR') then
    if(print_msg) then
      print*,'---------------------------------------------------------'
      print*,'|  ERROR run completed successfully. No fatal errors.'
      print*,'---------------------------------------------------------'
      print*,''
    endif
    return
  endif

  ! Run the time integration driver
  CALL rams_node ()

  ! Report total run times
  if(print_msg) then
    print*,''
    print*,''
  endif
  CALL timing (2,ctime_end)
  wtime_end = walltime(wtime_ref)
  CALL report_time_allnodes ('----Total run time: ' &
                  ,ctime_end-ctime_beg,wtime_end-wtime_beg)

return
END SUBROUTINE rams_model

!##################################################################################
Subroutine establish_config (nl_fname)

! This routine will set the rams configuration on all nodes. It does this by
! having the mainnum node read the namelist file, set eng_params and then
! broadcast the results to all the other nodes.

use node_mod
use mem_all

implicit none

  character (len=*) :: nl_fname

  integer :: i

  if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
    ! read the configuration specs and set parameters
    CALL read_nl (nl_fname)
    CALL eng_params ()

    ! print initial banner
    write(6,'(a1,78a1)') ' ',('*',i=1,78)
    write(6,'(2a1,a42)') ' ','*','    RAMS - Version 6.0'
    write(6,'(2a1)') ' ','*'
    write(6,'(3a1,A)') ' ','*',' ',trim(expnme)
    write(6,'(2a1)') ' ','*'
    write(6,'(2a1,A,A)') ' ','*',' RUNTYPE = ',trim(runtype)
    write(6,'(a1,78a1)') ' ',('*',i=1,78)
    if(iprntstmt>=1) then
      do i = 1, nmachs
        print*,'CPU:',machnum(i),mainnum
      enddo
      print*,'FILES:',nl_fname(1:len_trim(nl_fname))
      print*, ''
    endif

    ! Run the checks on the configuration variables
    CALL opspec1 ()
    CALL opspec3 ()
    if (trim(runtype) .eq. 'MAKEVFILE' .or. trim(runtype) .eq. 'MAKEHFILE') then
      CALL opspec4 ()
    endif

    ! Print report of var settings
    CALL nameout ()
  endif

  ! if more than the mainnum node, then send out results
  if (nmachs .gt. 1) then
    CALL broadcast_config ()
  endif

return
END SUBROUTINE establish_config

!##################################################################################
Subroutine set_full_domain_grid ()

! This routine will set the full domain grid specs for all nodes. It does this by
! having the mainnum node set up the grid and then broadcast the results 
! to all the other nodes.

use node_mod
use mem_grid

implicit none

  if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
    ! Fill in coordinate arrays, eg.
    CALL grid_setup (1)

    ! Run checks on the coarse/fine grid relationships
    CALL opspec2 ()
  endif

  ! if more than the mainnum node, then send out results
  if (nmachs .gt. 1) then
    CALL broadcast_grid ()
  endif

return
END SUBROUTINE set_full_domain_grid

!##################################################################################
Subroutine decompose_full_domain ()

! This routine will set the full domain grid specs for all nodes. It does this by
! having the mainnum node set up the grid and then broadcast the results 
! to all the other nodes.

use node_mod
use mem_grid

implicit none

  if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
    CALL assign_node_subdomain ()
  endif

  ! if more than the mainnum node, then send out results
  if (nmachs .gt. 1) then
    CALL broadcast_decomp ()
  endif

return
END SUBROUTINE decompose_full_domain

!##################################################################################
Subroutine update_sfc_var_files ()

! This routine will check any existing surface files to make sure that they
! are compatible with the RAMSIN grid specs. If not, or if files are missing,
! then new files will be created.
!
! Due to complications in getting individual subdomains to properly transform
! from the O-grid (in the lat/lon block files) to the R-grid (RAMS) with the intent
! to match the result from one node doing the entire R-grid, it was decided
! to continue having just one grid create the files. So if there are more than
! one node, then have the mainnum node call make_sfcfiles (which is coded
! to do the full domain) and force the other nodes to wait (par_pause) until
! the file creation is finished.
!
! In order to accompilish having one node do the surface file creation when there
! is actually more than one node, need to temporarily make mainnum look like it's
! the only node. Do this by saving the values of:
!    nmachs
!    mm[xyz]p
!    mem_read
!    mem_write
!    file_read
!    file_write
!    file_xchunk
!    file_ychunk
!
!    ct[ij]0
!    ctn[xyz]p
!    ct[ijk]pm
!
! and setting these to values according to the full domain. Once back from
! the make_sfcfiles call, restore these variables to their original values
! so that a subsequent simulation can run with multiple nodes.

use grid_dims
use node_mod
use mem_grid

implicit none

  integer :: orig_nmachs
  integer, dimension(maxgrds) :: orig_mmxp, orig_mmyp
  integer, dimension(maxgrds) :: orig_mmzp, orig_file_xchunk, orig_file_ychunk
  integer, dimension(maxgrds) :: orig_mmxyzp, orig_mmxysp, orig_mmxyp
  integer, dimension(maxgrds) :: orig_mi0, orig_mj0, orig_mibcon
  integer, dimension(maxgrds) :: orig_ctnxp, orig_ctnyp, orig_ctnzp
  integer, dimension(maxgrds) :: orig_cti0, orig_ctj0
  integer, dimension(nxpmax,maxgrds) :: orig_ctipm
  integer, dimension(nypmax,maxgrds) :: orig_ctjpm
  integer, dimension(nzpmax,maxgrds) :: orig_ctkpm
  type (io_descrip), dimension(maxgrds) :: orig_mem_read, orig_mem_write
  type (io_descrip), dimension(maxgrds) :: orig_file_read, orig_file_write
  integer :: igrid
  integer :: ctbeg, ctend
  integer :: i, j, k

  if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
    ! save the variables
    orig_nmachs      = nmachs
    orig_mmxp        = mmxp
    orig_mmyp        = mmyp
    orig_mmzp        = mmzp
    orig_mmxyzp      = mmxyzp
    orig_mmxysp      = mmxysp
    orig_mmxyp       = mmxyp
    orig_mi0         = mi0
    orig_mj0         = mj0
    orig_mibcon      = mibcon
    orig_mem_read    = mem_read
    orig_mem_write   = mem_write
    orig_file_read   = file_read
    orig_file_write  = file_write
    orig_file_xchunk = file_xchunk
    orig_file_ychunk = file_ychunk

    orig_ctnxp = ctnxp
    orig_ctnyp = ctnyp
    orig_ctnzp = ctnzp
    orig_cti0  = cti0
    orig_ctj0  = ctj0
    orig_ctipm = ctipm
    orig_ctjpm = ctjpm
    orig_ctkpm = ctkpm

    ! set to full domain values
    nmachs     = 1
    mmxp       = nnxp
    mmyp       = nnyp
    mmzp       = nnzp
    do igrid = 1, ngrids
      mmxyzp(igrid) = nnxp(igrid) * nnyp(igrid) * nnzp(igrid)
      mmxysp(igrid) = nnxp(igrid) * nnyp(igrid) * (nzg+nzs+3) * npatch
      mmxyp(igrid)  = nnxp(igrid) * nnyp(igrid)

      mi0(igrid)    = 0
      mj0(igrid)    = 0
      mibcon(igrid) = 15

      mem_read(igrid)%xblock = nnxp(igrid)
      mem_read(igrid)%yblock = nnyp(igrid)
      mem_read(igrid)%xoff   = 0
      mem_read(igrid)%yoff   = 0

      mem_write(igrid)%xblock = nnxp(igrid)
      mem_write(igrid)%yblock = nnyp(igrid)
      mem_write(igrid)%xoff   = 0
      mem_write(igrid)%yoff   = 0

      file_read(igrid)%xblock = nnxp(igrid)
      file_read(igrid)%yblock = nnyp(igrid)
      file_read(igrid)%xoff   = 0
      file_read(igrid)%yoff   = 0

      file_write(igrid)%xblock = nnxp(igrid)
      file_write(igrid)%yblock = nnyp(igrid)
      file_write(igrid)%xoff   = 0
      file_write(igrid)%yoff   = 0
      
      ! the minimum flag here is to make sure that we are writing
      ! a varfile and/or surface file that doesn't have chunks that are too
      ! big to be read into memory. 50x50 chunking is arbitrary, but in
      ! tests done by Sean Freeman (August 2020), it didn't seem to make
      ! a substantial difference in write or read time as long as the 
      ! chunk size is not too big. 
      file_xchunk(igrid) = min(50, nnxp(igrid))
      file_ychunk(igrid) = min(50, nnyp(igrid))

      if (igrid .eq. 1) then
         ! dummy values since grid 1 does not have a parent grid
         ctnxp(igrid) = nnxp(igrid)
         cti0(igrid)  = 0
         ctipm(1:nnxp(igrid),igrid) = 0

         ctnyp(igrid) = nnzp(igrid)
         ctj0(igrid)  = 0
         ctjpm(1:nnyp(igrid),igrid) = 0

         ctnzp(igrid) = nnzp(igrid)
         ctkpm(1:nnzp(igrid),igrid) = 0
      else
        ctbeg = ipm(1,igrid) - 2
        ctend = ipm(nnxp(igrid),igrid) + 1

        ctnxp(igrid) = (ctend - ctbeg) + 1
        cti0(igrid)  = ctbeg - 1
        do i = 1, nnxp(igrid)
          ctipm(i,igrid) = ipm(i,igrid) - cti0(igrid)
        enddo

        ctbeg = jpm(1,igrid) - 2
        ctend = jpm(nnyp(igrid),igrid) + 1

        ctnyp(igrid) = (ctend - ctbeg) + 1
        ctj0(igrid) = ctbeg - 1
        do j = 1, nnyp(igrid)
          ctjpm(j,igrid) = jpm(j,igrid) - ctj0(igrid)
        enddo

        ctnzp(igrid) = nnzp(igrid)
      endif
    enddo

    CALL make_sfcfiles ()
    if(trim(runtype) == 'MAKEVFILE') CALL isan_driver ()
    if(trim(runtype) == 'MAKEHFILE') CALL nudh_driver ()

    ! restore the variables
    nmachs      = orig_nmachs
    mmxp        = orig_mmxp
    mmyp        = orig_mmyp
    mmzp        = orig_mmzp
    mmxyzp      = orig_mmxyzp
    mmxysp      = orig_mmxysp
    mmxyp       = orig_mmxyp
    mi0         = orig_mi0
    mj0         = orig_mj0
    mibcon      = orig_mibcon
    mem_read    = orig_mem_read
    mem_write   = orig_mem_write
    file_read   = orig_file_read
    file_write  = orig_file_write
    file_xchunk = orig_file_xchunk
    file_ychunk = orig_file_ychunk

    ctnxp = orig_ctnxp
    ctnyp = orig_ctnyp
    ctnzp = orig_ctnzp
    cti0  = orig_cti0
    ctj0  = orig_ctj0
    ctipm = orig_ctipm
    ctjpm = orig_ctjpm
    ctkpm = orig_ctkpm
  endif

  ! If more than the mainnum node, then make everyone wait for the files
  ! to be built before proceding.
  if (nmachs .gt. 1) then
    CALL par_pause (my_rams_num, 201)
  endif

return
END SUBROUTINE update_sfc_var_files

!##############################################################################
Subroutine rams_output ()

use mem_leaf
use mem_varinit
use mem_cuparm
use io_params
use mem_grid
use node_mod

implicit none

integer :: ierr,ifm
logical :: analwrite
real :: timmaxr


do ifm = 1,ngrids
   CALL newgrid (ifm)
   timmaxr = timmax - .01*dtlongn(ifm)

! Implement call to NON_SCALAR_BC if calling ANLWRT.
! This sets all BC's for non-scalar (non-advected) quantities such
! as THETA, RV, and variables in LEAF, RADIATION, TURBULENCE.

!        For testing the mean fields, skip this rediagnosis
!        GOTO 10

   if (ioutput  ==  0) go to 10

   analwrite=.false.

   !FOR STATE FILE OUTPUT
   if(mod(time,frqstate(ifm)) <  dtlongn(ifm)               .or.  &
                        time  >= timmaxr                    .or.  &
                        iflag == 1                                &
      ) analwrite=.true.
   !FOR LITE FILE OUTPUT
   if(frqlite > 0.) then
      if( mod(time,frqlite) < dtlongn(1).or.  &
                         time  >=  timmaxr ) then
        analwrite=.true.
      endif
   endif
   !FOR MEAN FILE OUTPUT
   if (frqmean > 0.) then
      if(avgtim > 0.0.and.mod(time-avgtim/2.,frqmean) < dtlongn(1)  &
         .and.time >= avgtim) then
        analwrite=.true.
      endif
      if(avgtim < 0.0.and.mod(time,frqmean) < dtlongn(1)) then
        analwrite=.true.
      endif
   endif
   !FOR BOTH LITE AND MEAN FILE OUTPUT
   if (frqboth>0.) then
      if(avgtim > 0.0.and.mod(time-avgtim/2.,frqboth) < dtlongn(1)  &
         .and.time >= avgtim) analwrite=.true.
      if(avgtim < 0.0.and.mod(time,frqboth) < dtlongn(1))  &
         analwrite=.true.
   endif

   !If we are writing any output files, then do the non-scalar BC's.
   !We have to do this call within the ifm grid loop after call to 
   !"newgrid" so we have the grid appropriate nxp,nyp,nzp.
   if (analwrite) then
     CALL non_scalar_bc (mzp,mxp,myp,nzg,nzs,ibcon)
   endif

10      continue

enddo

timmaxr = timmax - .01*dtlongn(1)

analwrite=.false.
do ifm = 1,ngrids
   if(mod(time,frqstate(ifm)) < dtlongn(1).or.  &
      time  >=  timmaxr .or.  &
      iflag == 1) analwrite=.true.
enddo
if (analwrite) CALL anal_write ('INST')

!     call the analysis writing routine again for the other var types
if(frqlite > 0.) then
   if( mod(time,frqlite) < dtlongn(1).or.  &
                      time  >=  timmaxr ) then
    CALL anal_write ('LITE')
   endif
endif

if (frqmean > 0.) then
   if(iprntstmt>=1.and.print_msg)print*,'check avg',time,frqmean,avgtim,dtlongn(1)
   if(avgtim > 0.0.and.mod(time-avgtim/2.,frqmean) < dtlongn(1)  &
      .and.time >= avgtim) then
      CALL anal_write ('MEAN')
   endif
      
   if(avgtim < 0.0.and.mod(time,frqmean) < dtlongn(1)) then
      CALL anal_write ('MEAN')
   endif
endif

if (frqboth>0.) then
   if(avgtim > 0.0.and.mod(time-avgtim/2.,frqboth) < dtlongn(1)  &
      .and.time >= avgtim) CALL anal_write ('BOTH')
   if(avgtim < 0.0.and.mod(time,frqboth) < dtlongn(1))  &
      CALL anal_write ('BOTH')
endif

timmaxr = time+.00001

if (iupdsst  ==  1 .and. timmaxr < timmax) then
   do ifm = 1,ngrids
      CALL sst_read (3,ifm,ierr)
   enddo
endif

if (iupdndvi  ==  1 .and. timmaxr < timmax) then
   do ifm = 1,ngrids
      CALL ndvi_read (3,ifm,ierr)
   enddo   
endif

if(nud_type == 1 .and. time >= vtime2 .and. timmaxr < timmax) then
   CALL varf_read (2)
endif

if (iflag  ==  1) stop 'IFLAG'

return
END SUBROUTINE rams_output

!##################################################################################
Subroutine update_for_hist_grid (histgridflg,ngrids1,nnzp1,nnxp1,nnyp1 &
                                ,nzg1,nzs1,npatch1)

use grid_dims
use node_mod
use mem_grid

implicit none

integer :: histgridflg,ngrids1,nzg1,nzs1,npatch1
integer, save :: orig_nmachs
integer, dimension(maxgrds) :: nnzp1,nnxp1,nnyp1
integer, dimension(maxgrds),save :: orig_mmxp, orig_mmyp
integer, dimension(maxgrds),save :: orig_mmzp, orig_file_xchunk, orig_file_ychunk
integer, dimension(maxgrds),save :: orig_mmxyzp, orig_mmxysp, orig_mmxyp
type (io_descrip), dimension(maxgrds),save :: orig_mem_read,  orig_mem_write
type (io_descrip), dimension(maxgrds),save :: orig_file_read, orig_file_write
integer :: igrid

if(histgridflg==1) then

    ! save the variables
    orig_nmachs      = nmachs
    orig_mmxp        = mmxp
    orig_mmyp        = mmyp
    orig_mmzp        = mmzp
    orig_mmxyzp      = mmxyzp
    orig_mmxysp      = mmxysp
    orig_mmxyp       = mmxyp
    orig_mem_read    = mem_read
    orig_mem_write   = mem_write
    orig_file_read   = file_read
    orig_file_write  = file_write
    orig_file_xchunk = file_xchunk
    orig_file_ychunk = file_ychunk

    ! set to full domain values
    nmachs     = 1
    mmxp       = nnxp1
    mmyp       = nnyp1
    mmzp       = nnzp1
    do igrid = 1, ngrids1
      mmxyzp(igrid) = nnxp1(igrid) * nnyp1(igrid) * nnzp1(igrid)
      mmxysp(igrid) = nnxp1(igrid) * nnyp1(igrid) * (nzg1+nzs1+3) * npatch1
      mmxyp(igrid)  = nnxp1(igrid) * nnyp1(igrid)

      mem_read(igrid)%xblock = nnxp1(igrid)
      mem_read(igrid)%yblock = nnyp1(igrid)
      mem_read(igrid)%xoff   = 0
      mem_read(igrid)%yoff   = 0

      mem_write(igrid)%xblock = nnxp1(igrid)
      mem_write(igrid)%yblock = nnyp1(igrid)
      mem_write(igrid)%xoff   = 0
      mem_write(igrid)%yoff   = 0

      file_read(igrid)%xblock = nnxp1(igrid)
      file_read(igrid)%yblock = nnyp1(igrid)
      file_read(igrid)%xoff   = 0
      file_read(igrid)%yoff   = 0

      file_write(igrid)%xblock = nnxp1(igrid)
      file_write(igrid)%yblock = nnyp1(igrid)
      file_write(igrid)%xoff   = 0
      file_write(igrid)%yoff   = 0

      file_xchunk(igrid) = nnxp1(igrid)
      file_ychunk(igrid) = nnyp1(igrid)
    enddo

elseif(histgridflg==2) then

    ! restore the variables
    nmachs      = orig_nmachs
    mmxp        = orig_mmxp
    mmyp        = orig_mmyp
    mmzp        = orig_mmzp
    mmxyzp      = orig_mmxyzp
    mmxysp      = orig_mmxysp
    mmxyp       = orig_mmxyp
    mem_read    = orig_mem_read
    mem_write   = orig_mem_write
    file_read   = orig_file_read
    file_write  = orig_file_write
    file_xchunk = orig_file_xchunk
    file_ychunk = orig_file_ychunk

endif

return
END SUBROUTINE update_for_hist_grid
