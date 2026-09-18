!##############################################################################
Subroutine inithis (inithisflg)

use var_tables
use an_header
use mem_basic
use mem_grid
use ref_sounding
use io_params, only:hfilin
use micphys, only:level
use grid_struct
use mem_leaf, only:nvegpat
use hdf5_utils
use mem_varinit, only:igrid_match
use node_mod
use micro_prm, only:nkr

implicit none

integer :: ngrids1,nzg1,nzs1,npatch1,nvegpat1,ierr,ng,nc &
           ,ie,maxarr1,maxarr2,ngr,maxx1,maxy1,maxz1,npts,nv,nvh,nzpg1 &
           ,iyearh,imonthh,idateh,itimeh,ihtran1,checkhist,inithisflg,goahead
integer, external :: cio_i,cio_f
integer,save :: iunhd=11
integer, allocatable, dimension(:) :: nnxp1,nnyp1,nnzp1
real :: time1,ztop1,polelat1,polelon1
real, allocatable, dimension(:,:) :: xmn1,xtn1,ymn1,ytn1,zmn1,ztn1,topt1
real, allocatable, dimension(:) :: u01dn1,v01dn1,rt01dn1,th01dn1  &
                         ,pi01dn1,dn01dn1,scr,scr1
character(len=strl1) :: hnameinh
character(len=2) :: cng
character(len=1) :: cgrid
logical :: exists

integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

type (head_table), allocatable :: hr_table(:)

type(grid_def), allocatable :: grdefh(:)
type(grid_def), allocatable :: grdefn(:)

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

if(print_msg) print*,'Doing HISTORY FILE INITIALIZATION option.'

33 format(a30,2i5,3x,a18,i8)
34 format(a,i3,3x,a,i3,a)

!Get length of history file name
nc=len_trim(hfilin)

!Make sure new grid start time is the same as history file time
!If doing recycle run or history restart, times below do not need to match.
if(inithisflg==1)then
 read(hfilin(nc-14:nc-9),'(i6)') itimeh
 read(hfilin(nc-17:nc-16),'(i2)') idateh
 read(hfilin(nc-20:nc-19),'(i2)') imonthh
 read(hfilin(nc-25:nc-22),'(i4)') iyearh
 if(trim(runtype)/='HISTORY')then
  if(iyear1 /= iyearh .or. imonth1 /= imonthh .or. &
     idate1 /= idateh .or. itime1*100 /= itimeh) then
   if(print_msg) then
    print*,''
    print*,'Initialization date is not the same as history file date.'
    print*,'These need to be equal for consistency.'
    print*,'Current init time: ',iyear1,imonth1,idate1,itime1*100
    print*,'History file time: ',iyearh,imonthh,idateh,itimeh
   endif
   stop
  endif
 endif
endif

!Open RAMS history header file
CALL rams_f_open (iunhd,hfilin,'FORMATTED','OLD','READ',0)

!Get the sounding state (needed just for header file output)
ie=cio_i(iunhd,1,'iref',iref,1)
ie=cio_i(iunhd,1,'jref',jref,1)
ie=cio_f(iunhd,1,'topref',topref,1)
ie=cio_i(iunhd,1,'nsndg',nsndg,1)
ie=cio_f(iunhd,1,'us',us,nsndg)
ie=cio_f(iunhd,1,'vs',vs,nsndg)
ie=cio_f(iunhd,1,'ts',ts,nsndg)
ie=cio_f(iunhd,1,'thds',thds,nsndg)
ie=cio_f(iunhd,1,'rts',rts,nsndg)
ie=cio_f(iunhd,1,'ps',ps,nsndg)
ie=cio_f(iunhd,1,'hs',hs,nsndg)
if(io3flg==1)ie=cio_f(iunhd,1,'o3s',o3s,nsndg)

!Get original simulation type
ie=cio_i(iunhd,1,'initorig',initorig,1)

!Get grid specs for grid comparison
ie=cio_i(iunhd,1,'ngrids',ngrids1,1)

allocate (nnxp1(ngrids1),nnyp1(ngrids1),nnzp1(ngrids1))

!Get history grid structure info so we can allocate space
ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
ie=cio_i(iunhd,1,'npatch',npatch1,1)
ie=cio_i(iunhd,1,'nvegpat',nvegpat1,1)
ie=cio_i(iunhd,1,'nzg',nzg1,1)
ie=cio_i(iunhd,1,'nzs',nzs1,1)
ie=cio_f(iunhd,1,'time',time1,1)
ie=cio_f(iunhd,1,'ztop',ztop1,1)
ie=cio_f(iunhd,1,'polelat',polelat1,1)
ie=cio_f(iunhd,1,'polelon',polelon1,1)
ie=cio_i(iunhd,1,'ihtran',ihtran1,1)

!Make sure new grid polelat and polelon match original grid
if(polelat /= polelat1 .or. polelon /= polelon1) then
  if(print_msg) then
   print*,'The POLELAT and POLELON of the new grid much match original grid.'
   print*,'Past POLELAT:',polelat1,' New POLELAT:',polelat
   print*,'Past POLELON:',polelon1,' New POLELON:',polelon
  endif
  stop
endif

!Make sure new grid top does not exceed old grid top
if(ztop > ztop1) then
  if(print_msg) then
   print*,'New interpolated grid top higher than past grid top!'
   print*,'Cannot interpolate. Extrapolation not recommended.'
   print*,'Past top:',ztop1,' New top:',ztop
  endif
  stop
endif

!If adding a grid on history restart or history initializing from history
!file that had progressed beyond the first timestep, then set ngbegun to 1.
!If you history initialize from a previous history initialization rathan than
!an original simulations, then the condition below might not be met if
!performing a new history initilization from the time=0 of the previous
!history initialization. This would be a rather odd configuration but would
!force the use of the default first-timestep "filter parameter" "epsu" used
!in the Asselin Filter for acoustic timestep fields.
if(time1 > 0.0) ngbegun(1:ngrids) = 1

!Check to make sure surface parameters are the same
if (nzg /= nzg1 .or. nzs /= nzs1 .or. npatch /= npatch1 .or. &
    nvegpat /= nvegpat1) then
  if(print_msg) then
   print*,'LEAF parameters must be same for initial-history start'
   print*,'npatch: ',npatch,npatch1
   print*,'nvegpat: ',nvegpat,nvegpat1
   print*,'nzg: ',nzg,nzg1
   print*,'nzs: ',nzs,nzs1
  endif
  stop
endif
!Check to make sure grid type is the same
if (ihtran /= ihtran1) then
  if(print_msg) then
   print*,'GRID TYPE (Cartesian or Polar Stereographic) must be the same'
   print*,'between history grid and current grid. May need to consider'
   print*,'wind rotation and such if we allow a change.'
  endif
  stop
endif
!Make sure original simulations type is the same
if (initorig /= 1 .and. initorig /= 2) then
  if(print_msg) then
   print*,'We do not have info on initial simulation type.'
   print*,'Need to know if history was HH or Variable-Init.'
  endif
  stop
endif

!Find maximum size of any array on history file. Allocate scratch arrays.
maxarr1=0
maxarr2=0
maxx1=maxval(nnxp1(1:ngrids1))
maxy1=maxval(nnyp1(1:ngrids1))
maxz1=maxval(nnzp1(1:ngrids1))
do ngr=1,ngrids1
   maxarr1=max(maxarr1,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr)*nkr  &
         ,nnxp1(ngr)*nnyp1(ngr)*nzg1*npatch1 &
         ,nnxp1(ngr)*nnyp1(ngr)*nzs1*npatch1)
   maxarr2=max(maxarr2,nnxp1(ngr)*nnyp1(ngr))
enddo
allocate (scr(maxarr1),scr1(maxarr1))

!Allocate and read in grid specifications
allocate(xmn1(maxx1,ngrids1),xtn1(maxx1,ngrids1))
allocate(ymn1(maxy1,ngrids1),ytn1(maxy1,ngrids1))
allocate(zmn1(maxz1,ngrids1),ztn1(maxz1,ngrids1))
do ngr=1,ngrids1
   write(cng,'(i2.2)') ngr
   ie=cio_f(iunhd,1,'xmn'//cng,xmn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'xtn'//cng,xtn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'ymn'//cng,ymn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'ytn'//cng,ytn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'zmn'//cng,zmn1(1,ngr),nnzp1(ngr))
   ie=cio_f(iunhd,1,'ztn'//cng,ztn1(1,ngr),nnzp1(ngr))
enddo

!Read variable header info
rewind(iunhd)

!Top of header file. Total number of variables on all grids.
read(iunhd,*) nvbtab
!Allocate space for header file variable information
allocate (hr_table(nvbtab))
!Reader info on each variables in header file.
!(ex. Name  IfAPointer  Dimension   Grid   #GridPoints)  
!(ex.  UP       0           3        1        49000)
do nv=1,nvbtab
  read(iunhd,*)              &
     hr_table(nv)%string     &
    ,hr_table(nv)%npointer   &
    ,hr_table(nv)%idim_type  &
    ,hr_table(nv)%ngrid      &
    ,hr_table(nv)%nvalues
enddo

!Check to see that all history grids are present at this time
checkhist=1
do ngr=1,ngrids1
  write(cgrid,'(i1)') ngr
  hnameinh=hfilin(1:nc-9)//'-g'//cgrid//'.h5'
  inquire(file=hnameinh,exist=exists)
  if(.not. exists)then
   checkhist=0
  endif
enddo
if(checkhist==0)then
 if(print_msg) then
  print*,'Not all original grids are present at this time: ',hfilin(1:nc-9)
  print*,'Choose a history initialization time in which all original'
  print*,'simulation grids are available.'
  print*,''
 endif
 stop
endif

!Allocate history topography and get TOPT
allocate (topt1(maxarr2,ngrids1))
do ngr=1,ngrids1
   write(cgrid,'(i1)') ngr
   hnameinh=hfilin(1:nc-9)//'-g'//cgrid//'.h5'
   CALL update_for_hist_grid (1,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
   if (nmachs .gt. 1) then
     iphdf5 = 1
   else
     iphdf5 = 0
   endif
   CALL shdf5_open (hnameinh,'R',iphdf5,h5_fid)
   CALL shdf5_set_hs_select (2,'R',ngr &
          ,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'TOPT' &
          ,mem_select,file_select,rvara=topt1(1,ngr))
   CALL shdf5_close (h5_fid)   
   CALL update_for_hist_grid (2,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
enddo

!******************************************************************************
!********* START CHECK GRID MATCHING ******************************************
!******************************************************************************
! Set a flag array (for each grid on the history file) to determine:
!  = 1 = This grid is identical to a current grid
!  = 0 = This grid is different.

igrid_match(1:maxgrds,1:maxgrds)=0

if(print_msg) then
 print*,'###################################################################'
 print*,'# History Grids:',ngrids1
 print*,'# Current Grids:',ngrids
endif

!Allocate grid structures
allocate (grdefn(ngrids))
allocate (grdefh(ngrids1))

!Define the History Grid
do ngr=1,ngrids1
   CALL alloc_grid_def ( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) )
   CALL fill_grid_def  ( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr)   &
                       ,nzg1,nzs1,npatch1,nvegpat1,polelat1,polelon1  &
                       ,xtn1(1,ngr),xmn1(1,ngr),ytn1(1,ngr),ymn1(1,ngr) &
                       ,ztn1(1,ngr),zmn1(1,ngr),topt1(1,ngr) )
enddo
!Define the Current Initial Grid
do ngr=1,ngrids
   CALL alloc_grid_def ( grdefn(ngr),mmxp(ngr),mmyp(ngr),mmzp(ngr) )
   CALL fill_grid_def  ( grdefn(ngr),mmxp(ngr),mmyp(ngr),mmzp(ngr)  &
                       ,nzg,nzs,npatch,nvegpat,polelat,polelon    &
                       ,xtn(1,ngr),xmn(1,ngr),ytn(1,ngr),ymn(1,ngr) &
                       ,ztn(1,ngr),zmn(1,ngr),grid_g(ngr)%topt )
enddo

! See if the history grids match any of the new grids...assuming 1:1 grid number
!   correspondence for now
do ngr=1,ngrids1
 do ng=1,ngrids
   CALL compare_grid_def (grdefh(ngr),grdefn(ng),'nud_update',ierr,ngr,ng)
   if (ierr /= 0) then
      ! No match...
      if(print_msg) print*,'These grids do NOT match (Hist,New):',ngr,ng
   else
      ! We have a match...
      igrid_match(ngr,ng)=1
      if(print_msg) print*,'These grids DO match (Hist,New):',ngr,ng
   endif
 enddo
enddo

if(print_msg) then
 print*,'###################################################################'
endif

!******************************************************************************
!******************* END CHECK GRID MATCHING **********************************
!******************************************************************************

!******************* START READING IN THE DATA ***************************
! Finally, process the fields...

!*************************************************************************
!The code below (not including reference state section) interpolates all 
!history grids to all current grids in the model if a current grid fits 
!fully inside a history grid. Variables specific to the current grid that
!come from surface files (via MAKESFC) are not over-written or interpolated
!from history files. These include TOPT, TOPZO, LEAF_CLASS, PATCH_AREA,
!SOIL_TEXT, and VEG_NDVIF. Further, GLAT and GLON from the current grids
!are not over-written.
!*************************************************************************

!*************************************************************************
!Loop over history grids "ngr" (ngrids1 = number of grids in history file)
!*************************************************************************
do ngr=1,ngrids1

  !Get history file name and open the file
  write(cgrid,'(i1)') ngr
  hnameinh=hfilin(1:nc-9)//'-g'//cgrid//'.h5'

  !*************************************************************************
  !Loop over current run grids "ng" (ngrids = number of grids in current run)
  !*************************************************************************
  do ng=1,ngrids

   if (xtn(1,ng) < xtn1(1,ngr) .or. xtn(mmxp(ng),ng) > xtn1(nnxp1(ngr),ngr) .or. &
       ytn(1,ng) < ytn1(1,ngr) .or. ytn(mmyp(ng),ng) > ytn1(nnyp1(ngr),ngr) ) then

    if(print_msg) then
     print*,'###################################################################'
     print 34,' Current Grid:',ng,'outside bounds of History Grid:',ngr,': NO Interp'
     print*,'Will not interpolate to a new grid bigger than hist grid.'
     print*,'Current grid larger than all history grids.'
     print*,'###################################################################'
    endif

    if(ngr==1)stop

   else 

    if(print_msg) then
     print*,'###################################################################'
     print 34,' Current Grid:',ng,' within bounds of History Grid:',ngr,': Interp'
     print*,'###################################################################'
    endif

    !Loop through all variables in the history files (for all grids)
    varloop: do nvh=1,nvbtab

      !Cycle if the history variable is not on the grid we are looping over
      if(ngr /= hr_table(nvh)%ngrid) cycle varloop

      !num_var = Global variable of number of variables in the current run.
      !vtab_r = Current run variables info. List includes extra variables not
      !         contained in the history file that are computed within.
      !Find which is the corresponding variable in the current run
      do nv = 1,num_var(ng)

       !Continue if doing full history initialization or recycling surface vars.
       if( (inithisflg == 0 .and. vtab_r(nv,ng)%irecycle_sfc == 1) .or. &
           (inithisflg == 1) ) then

         !Number of grid points for the current run variable
         npts=vtab_r(nv,ng)%npts

         !See if this variable is active in the current run and interpolate
         !to new grid stucture. Some variables are excluded. For surface
         !variables we only copy/interpolate the ones that are "recycleable" as
         !specified in memory routines; the idea is not to override values from
         !surface data files.

         if(hr_table(nvh)%string == vtab_r(nv,ng)%name) then

          !Copy current grid structure to temporary arrays and tell
          !HDF5 to look for history file grid structure to read.
          CALL update_for_hist_grid (1,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
          if (nmachs .gt. 1) then
            iphdf5 = 1
          else
            iphdf5 = 0
          endif
          !Open file, Read data from the history files, Close file
          CALL shdf5_open (hnameinh,'R',iphdf5,h5_fid)
          CALL shdf5_set_hs_select (vtab_r(nv,ng)%idim_type,'R',ngr &
                 ,mem_select,file_select,file_chunks)
          CALL shdf5_irec (h5_fid,iphdf5,trim(hr_table(nvh)%string) &
                 ,mem_select,file_select,rvara=scr1)
          CALL shdf5_close (h5_fid)
          !Copy current grid structure back to main grid arrays.
          CALL update_for_hist_grid (2,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)

          !Arrange data for RAMS to work with to copy or interpolate based
          !on the dimensions of the variable ingested.

          !If this is a history INITIALIZATION we cannot run in parallel
          !since history grids can be different from current grids.
          !NMACHS will have already been forced to = 1, so we can proceed
          !below using the number of grid points on history grid (nnxp1,nnyp1).
          select case(hr_table(nvh)%idim_type)
            case(2,6)
               CALL atob (hr_table(nvh)%nvalues,scr1,scr)
            case(3)
               CALL unarrange (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),scr1,scr)
            case(4)
               CALL unarrange_p (nnxp1(ngr),nnyp1(ngr),nzg1,npatch1,scr1,scr)
            case(5)
               CALL unarrange_p (nnxp1(ngr),nnyp1(ngr),nzs1,npatch1,scr1,scr)
            case(7)
               CALL unarrange_p (nnxp1(ngr),nnyp1(ngr),nnzp1(ngr),nkr,scr1,scr)
          end select

          !*******************************************************************
          if (vtab_r(nv,ng)%idim_type == 2) then
               goahead=0
               !For history RESTART:
               if(trim(runtype)=='HISTORY' .and. ng <= ngrids1) goahead=1
               !For history INITIALIZATION or history restart added grid:
               !Do not override Topography from Topo surface files.
               !Do not override current grid lat and lon.
               !Do not bring over history accumulated precipitation. Reset.
               if(trim(runtype)/='HISTORY' .or. ng > ngrids1) then
                 if(hr_table(nvh)%string /= 'TOPT'   .and. &
                    hr_table(nvh)%string /= 'TOPZO'  .and. &
                    hr_table(nvh)%string /= 'GLAT'   .and. &
                    hr_table(nvh)%string /= 'GLON'   .and. &
                    hr_table(nvh)%string /= 'ACCPD'  .and. &
                    hr_table(nvh)%string /= 'ACCPR'  .and. &
                    hr_table(nvh)%string /= 'ACCPP'  .and. &
                    hr_table(nvh)%string /= 'ACCPS'  .and. &
                    hr_table(nvh)%string /= 'ACCPA'  .and. &
                    hr_table(nvh)%string /= 'ACCPG'  .and. &
                    hr_table(nvh)%string /= 'ACCPH'  .and. &
                    hr_table(nvh)%string /= 'ACONPR') goahead=1
               endif
               if(goahead==1) then
                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 2',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 2',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (1,nnxp1(ngr),nnyp1(ngr),1,1,scr(1)       &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)    &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1         &
                         ,1,mmxp(ng),mmyp(ng),1,1,vtab_r(nv,ng)%var_p        &
                         ,ng,ngr,vtab_r(nv,ng)%name,2)
                  endif
               endif
          !*******************************************************************
          elseif (vtab_r(nv,ng)%idim_type == 3) then

                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 3',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 3',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,1,scr(1) &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)       &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1            &
                         ,mmzp(ng),mmxp(ng),mmyp(ng),1,1,vtab_r(nv,ng)%var_p    &
                         ,ng,ngr,vtab_r(nv,ng)%name,3)
                  endif

          !**********************************************************************
          elseif (vtab_r(nv,ng)%idim_type == 4 .and. &
                  vtab_r(nv,ng)%irecycle_sfc  == 1) then

                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 4-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 4-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (nzg1,nnxp1(ngr),nnyp1(ngr),npatch1,1,scr(1) &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)       &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1            &
                         ,nzg,mmxp(ng),mmyp(ng),npatch,1,vtab_r(nv,ng)%var_p    &
                         ,ng,ngr,vtab_r(nv,ng)%name,4)
                  endif

          !*********************************************************************
          elseif (vtab_r(nv,ng)%idim_type == 5 .and. &
                  vtab_r(nv,ng)%irecycle_sfc  == 1) then

                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 5-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 5-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (nzs1,nnxp1(ngr),nnyp1(ngr),npatch1,1,scr(1) &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)       &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1            &
                         ,nzs,mmxp(ng),mmyp(ng),npatch,1,vtab_r(nv,ng)%var_p    &
                         ,ng,ngr,vtab_r(nv,ng)%name,5)
                  endif

          !**********************************************************************
          elseif (vtab_r(nv,ng)%idim_type == 6 .and. &
                  vtab_r(nv,ng)%irecycle_sfc  == 1) then

                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 6-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 6-1',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (1,nnxp1(ngr),nnyp1(ngr),npatch1,1,scr(1)  &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)     &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1          &
                         ,1,mmxp(ng),mmyp(ng),npatch,1,vtab_r(nv,ng)%var_p    &
                         ,ng,ngr,vtab_r(nv,ng)%name,6)
                  endif

          !*******************************************************************
          elseif (vtab_r(nv,ng)%idim_type == 7) then

                  if(igrid_match(ngr,ng)==1) then
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: copy: 7',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL atob (npts,scr(1),vtab_r(nv,ng)%var_p)
                  else
                    if(iprntstmt>=1 .and. print_msg) &
                      print 33,'init_update: interp: 7',ngr,ng,vtab_r(nv,ng)%name,nv
                    CALL hi_interp (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,nkr,scr(1) &
                         ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)         &
                         ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1              &
                         ,mmzp(ng),mmxp(ng),mmyp(ng),1,nkr,vtab_r(nv,ng)%var_p    &
                         ,ng,ngr,vtab_r(nv,ng)%name,7)
                  endif

          endif !end if dimension type
          exit
         endif !end if variable is active on current grid

       endif !if a variable to recycle or inithis
      enddo !do looping over variables on current grid

    enddo varloop !do looping over variables on history file
   endif !if current grid fully contained in history grid
  enddo !do looping over current run grids

enddo !do looping over history grids

!*************************************************************************
!The following code for reference state variables is set up for
!grid-1. If there is more than one grid, the model will interpolate
!from grid-1 to additional grids. This is done in rdint.f90 routine "initlz".
!*************************************************************************

!Prepare 1D reference sounding for history grid-1 and apply to
! new grid-1. Will let interp_fine_grid interpolate to new nested grids.
if(inithisflg == 1) then

 nzpg1=nnzp1(1)
 allocate(u01dn1(nzpg1), v01dn1(nzpg1),rt01dn1(nzpg1)  &
        ,th01dn1(nzpg1),pi01dn1(nzpg1),dn01dn1(nzpg1) )

 cng='01'
 ie=cio_f(iunhd,1,'u01dn'//cng,  u01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'v01dn'//cng,  v01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'pi01dn'//cng,pi01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'th01dn'//cng,th01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'dn01dn'//cng,dn01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'rt01dn'//cng,rt01dn1(1),nzpg1)

 CALL htint (nzpg1,th01dn1,ztn1(1,1),mmzp(1),th01dn(1,1),ztn(1,1))
 CALL htint (nzpg1,u01dn1 ,ztn1(1,1),mmzp(1),u01dn(1,1) ,ztn(1,1))
 CALL htint (nzpg1,v01dn1 ,ztn1(1,1),mmzp(1),v01dn(1,1) ,ztn(1,1))
 CALL htint (nzpg1,pi01dn1,ztn1(1,1),mmzp(1),pi01dn(1,1),ztn(1,1))
 CALL htint (nzpg1,dn01dn1,ztn1(1,1),mmzp(1),dn01dn(1,1),ztn(1,1))

 if (level .ge. 1) then
    CALL htint (nzpg1,rt01dn1,ztn1(1,1),mmzp(1),rt01dn(1,1),ztn(1,1))
 else
    rt01dn(1:mmzp(1),1) = 0.
 endif

 u01dn(1,1) = u01dn(2,1)
 v01dn(1,1) = v01dn(2,1)
 rt01dn(1,1) = rt01dn(2,1)
 th01dn(1,1) = th01dn(2,1)

 !Close the input history header file
 close(iunhd)

 !Compute 3d reference state for grid 1
 CALL newgrid (1)
 CALL refs3d (mmzp(1),mmxp(1),mmyp(1)  &
   ,basic_g(1)%pi0  (1,1,1)  ,basic_g(1)%dn0  (1,1,1)  &
   ,basic_g(1)%dn0u (1,1,1)  ,basic_g(1)%dn0v (1,1,1)  &
   ,basic_g(1)%th0  (1,1,1)  ,grid_g(1)%topt  (1,1)    &
   ,grid_g(1)%rtgt  (1,1)                                )

 !Deallocate temporary arrays
 deallocate(u01dn1,v01dn1,rt01dn1,th01dn1,pi01dn1,dn01dn1)

 !Calculate number of levels for radiation profiles using a combination 
 !of the model prognostic grid and the input sounding. If ozone is present,
 !also fill the ozone reference profile above the prognostic model top
 CALL calc_refs()

 !Calculate large-scale subsidence
 CALL calc_wsub()
endif !If doing full history initialization

!Deallocate temporary arrays
deallocate(topt1,hr_table,xmn1,ymn1,zmn1,xtn1,ytn1,ztn1,scr,scr1 &
          ,nnxp1,nnyp1,nnzp1)

do ngr=1,ngrids1
   CALL dealloc_grid_def (grdefh(ngr))
enddo
do ngr=1,ngrids
   CALL dealloc_grid_def (grdefn(ngr))
enddo

return
END SUBROUTINE inithis

!##############################################################################
Subroutine sfcinit_hstart ()
                     
use mem_leaf
use mem_basic
use mem_grid
use leaf_coms
use rconstants
use mem_varinit, only: igrid_match
use node_mod

implicit none

integer :: i,j,ifm,ipat,k,nveg,nsoil,hifm,maxmatch
real :: c1,hpis,hprss

if (isfcl == 3) return
! This routine fills the LEAF arrays for a history-initial start.

! Refill many of the LEAF variables, as the interpolated values may not be relevant.

c1 = .5 * cpi

ifmloop: do ifm=1,ngrids

! Only re-init if no history grids match current grids.
! If there is a match then cycle loop and go to next current grid.
maxmatch=0
do hifm=1,maxgrds
 maxmatch=max(maxmatch,igrid_match(hifm,ifm))
enddo

if(maxmatch /= 0) cycle ifmloop

if(iprntstmt>=1 .and. print_msg) then
 print*,'##############################################################'
 print*,'Reinitializing some surface variables for grid:',ifm
 print*,'##############################################################'
endif

do j = 1,mmyp(ifm)
   do i = 1,mmxp(ifm)
      hpis = c1 * (basic_g(ifm)%pi0(1,i,j) + basic_g(ifm)%pi0(2,i,j)   &
                      + basic_g(ifm)%pp(1,i,j) + basic_g(ifm)%pp(2,i,j))
 
      hprss = hpis ** cpor * p00

      leaf_g(ifm)%patch_rought(i,j,1) = 0.001
      leaf_g(ifm)%patch_roughm(i,j,1) = 0.001

      do ipat = 2,npatch

         !Set vegetation height since leaf_init is not run for hist init.
         nveg = nint(leaf_g(ifm)%leaf_class(i,j,ipat))
         leaf_g(ifm)%veg_height(i,j,ipat) = veg_ht(nveg)

         !Reset soil roughness, patch roughness, and stomatal resistance
         !since leaf class may change due to interpolation. Land patch may
         !now be a small fraction (<.009) or could be significantly different
         !such as water patch. These values will be correctly set 1st timestep.
         leaf_g(ifm)%soil_rough(i,j,ipat) = ztrough
         leaf_g(ifm)%patch_rought(i,j,ipat) = max(ztrough,grid_g(ifm)%topzo(i,j))
         leaf_g(ifm)%patch_roughm(i,j,ipat) = max(zmrough,grid_g(ifm)%topzo(i,j))
         leaf_g(ifm)%stom_resist(i,j,ipat) = 1.e6

         !Recompute soil moisture since we changed soil classes. Must make
         !sure soil moisture does not exceed 100% due to difference in
         !soil class on different history initialization grid.
         do k = 1,nzg
           nsoil = nint(leaf_g(ifm)%soil_text(k,i,j,ipat))
           leaf_g(ifm)%soil_water(k,i,j,ipat) = &
             min(slmsts(nsoil),leaf_g(ifm)%soil_water(k,i,j,ipat))
         enddo

         !Recompute surface water levels since this variable is a float
         !but needs and integer value. Intepolation from former grid
         !would create floats rather than integers.
         do k = 1,nzs
            leaf_g(ifm)%sfcwater_nlev(i,j,ipat) = 0.
            if (leaf_g(ifm)%sfcwater_mass(k,i,j,ipat) > 0.)  &
               leaf_g(ifm)%sfcwater_nlev(i,j,ipat) = float(k)
         enddo

         !Initialize for most land surface options (isfcl=0,1,2,not 3)
         !Update this since leaf classes and ndvi have changes due to
         !history initialization on different grid.
         if (ipat >= 2) CALL ndvi (ifm  &
            ,leaf_g(ifm)%leaf_class (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvip  (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvic  (i,j,ipat)   &
            ,leaf_g(ifm)%veg_ndvif  (i,j,ipat)   )

         !Only call this for LEAF3. SIB will use its own values.
         !Update this since leaf classes are likely different following
         !history initialization on different grid.
         if (ipat >= 2 .and. isfcl<=1)           &
            CALL veg (ifm                        &
            ,leaf_g(ifm)%leaf_class   (i,j,ipat) &
            ,leaf_g(ifm)%veg_fracarea (i,j,ipat) &
            ,leaf_g(ifm)%veg_lai      (i,j,ipat) &
            ,leaf_g(ifm)%veg_tai      (i,j,ipat) &
            ,leaf_g(ifm)%veg_rough    (i,j,ipat) &
            ,leaf_g(ifm)%veg_height   (i,j,ipat) &
            ,leaf_g(ifm)%veg_albedo   (i,j,ipat) &
            ,leaf_g(ifm)%veg_ndvic    (i,j,ipat) )

         !Must recompute here to prevent over saturation and such
         !since soil_water and sfcwater_nlev could have been updated
         !due to change in soil class and snow integer level interpolation
         CALL grndvap (nzs  &
            ,leaf_g(ifm)%soil_energy    (nzg,i,j,ipat)  &
            ,leaf_g(ifm)%soil_water     (nzg,i,j,ipat)  &
            ,leaf_g(ifm)%soil_text      (nzg,i,j,ipat)  &
            ,leaf_g(ifm)%sfcwater_energy(1,i,j,ipat)    &
            ,leaf_g(ifm)%sfcwater_nlev  (i,j,ipat)      &
            ,leaf_g(ifm)%ground_rsat    (i,j,ipat)      &
            ,leaf_g(ifm)%ground_rvap    (i,j,ipat)      &
            ,leaf_g(ifm)%can_rvap       (i,j,ipat)      &
            ,hprss,i,j)

      enddo
   enddo
enddo

enddo ifmloop

return
END SUBROUTINE sfcinit_hstart

!##############################################################################
Subroutine hi_interp (n1,n2,n3,n4,n5,vn,xm1,xt1,ym1,yt1,zm1,zt1  &
                    ,topt1,ztop1,m1,m2,m3,m4,m5,vm,ngm,ngr1,vname,idim)
                      
use mem_grid
use mem_scratch

implicit none

integer :: n1,n2,n3,n4,n5,ngr1,idim
real, dimension(n1,n2,n3,n4,n5) :: vn
real, dimension(n2,n3) :: topt1
real, dimension(n1) :: zm1,zt1
real, dimension(n2) :: xm1,xt1
real, dimension(n3) :: ym1,yt1
real :: ztop1

integer :: m1,m2,m3,m4,m5,ngm
real, dimension(m1,m2,m3,m4,m5) :: vm

character(len=*) :: vname

integer :: i,j,k,np,nb
real :: fixxm,fiyym,topoh,rtgth
real, allocatable :: scr(:,:,:,:,:)

!This routine will interpolate the new simulations grids
!Note here: n1,n2,n3,n4 are the original grid specs and
! m1,m2,m3,m4 are the new grid specs we are interpolating to.

if(allocated(scr)) deallocate(scr) ; allocate(scr(n2,n3,n1,n4,n5))

!We are going to swap indices due to excessive memory copies of
!non-contiguous memory
do nb=1,n5
 do np=1,n4
   do k=1,n1
      scr(1:n2,1:n3,k,np,nb)=vn(k,1:n2,1:n3,np,nb)
   enddo
 enddo
enddo

do j=1,m3
   do i=1,m2

      !See if point is on this grid.
      if (xtn(1,ngm) < xt1(1) .or. xtn(m2,ngm) > xt1(n2) .or. &
          ytn(1,ngm) < yt1(1) .or. ytn(m3,ngm) > yt1(n3) ) then
         !This is the input coarse grid and point is not on this grid. Stop.
         if (ngr1 == 1) then
            print*,''
            print*,'Hist grid: ',ngr1,' : ','Current grid: ',ngm
            print*,'His_init: grid point not on history file grids'
            print*,i,xtn(1,ngm),xt1(1),xtn(m2,ngm),xt1(n2)
            print*,j,ytn(1,ngm),yt1(1),ytn(m3,ngm),yt1(n3)
            stop 'His_init: hi_interp'
         else
            !Otherwise, go to next point
            cycle
         endif
      endif
          
      !We are okay horizontally, now interpolate vertical column from field 
      
      !Find x,y grid point locations on input field.
      !Assuming constant spacing and deal with stagger for winds
      
      if(vname == 'UP' .or. vname == 'UC') then
         !Handle U winds on M-grid in horizontal in X direction
         fixxm=1.+(xmn(i,ngm)-xm1(1))/(xm1(2)-xm1(1))
      else
         !Handle scalars on T-grid in horizontal in X direction
         fixxm=1.+(xtn(i,ngm)-xt1(1))/(xt1(2)-xt1(1))
      endif
      if(vname == 'VP' .or. vname == 'VC') then
         !Handle V winds on M-grid in horizontal in Y direction
         fiyym=1.+(ymn(j,ngm)-ym1(1))/(ym1(2)-ym1(1))
      else
         !Handle scalars on T-grid in horizontal in X direction
         fiyym=1.+(ytn(j,ngm)-yt1(1))/(yt1(2)-yt1(1))
      endif

      !Set this to 1.0 for 2D simulations
      if(m3==1)fiyym=1.

      do nb=1,n5
       do np=1,n4
        do k=1,n1
         CALL gdtost2 (scr(1,1,k,np,nb),n2,n3,fixxm,fiyym,vctr1(k))
        enddo

        if (idim == 3 .or. idim == 7) then
         !Interpolate this column vertically to actual grid if 3d variable
         CALL gdtost2 (topt1(1,1),n2,n3,fixxm,fiyym,topoh)    
         rtgth=1.-topoh/ztop1
         if(vname == 'WP' .or. vname == 'WC') then
           !Adjust for topography on M-grid in the vertical
           do k=1,m1
            !Actual grid level heights
            vctr2(k)= grid_g(ngm)%topt(i,j)   &
                      + zmn(k,ngm) *grid_g(ngm)%rtgt(i,j)
           enddo
           do k=1,n1
            !History grid level heights
            vctr3(k)= topoh + zm1(k) *rtgth
           enddo
         else
           !Adjust for topography on T-grid in the vertical
           do k=1,m1
            !Actual grid level heights
            vctr2(k)= grid_g(ngm)%topt(i,j)   &
                      + ztn(k,ngm) *grid_g(ngm)%rtgt(i,j)
           enddo
           do k=1,n1
            !History grid level heights
            vctr3(k)= topoh + zt1(k) *rtgth
           enddo
         endif

         !Interpolate vertically
         CALL htint (n1,vctr1(1),vctr3(1),m1,vctr10(1),vctr2(1))

         vm(1:m1,i,j,np,nb)=vctr10(1:m1)

        elseif (idim == 2) then
         vm(1,i,j,np,nb)=vctr1(1)
        elseif (idim == 4 .or. idim == 5) then
         vm(1:m1,i,j,np,nb)=vctr1(1:m1)
        elseif (idim == 6) then
         vm(1,i,j,np,nb)=vctr1(1)
        endif

       enddo
      enddo

   enddo
enddo

deallocate(scr)

return
END SUBROUTINE hi_interp
