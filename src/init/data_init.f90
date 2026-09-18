!##############################################################################
Subroutine handle_dustsource (n2,n3,glat,glon,dustfrac,ifm)

use micphys, only:dustfile,idustloft
use node_mod
use mem_grid

implicit none

! This routine will manage reading of the dust source data.
! Have only mainnum do the file I/O. After the file is read by mainnum, 
! then broadcast the data to the other nodes.

integer :: i,j,isource,jsource,n2,n3,ifm,ii,jj,nx_source,ny_source,setnxny
character(len=1) :: dgrid
real, dimension(n2,n3) :: glat,glon,dustfrac
character(len=strl1) :: cname,dname
real :: glat1,glon1,glat2,glon2,sumsource,sumnum
real, allocatable, dimension(:) :: lat_source,lon_source
real, allocatable, dimension(:,:) :: source

!Standard DUSTFILE name from RAMSIN. This is used when not needing to
!create grid-specific dust file or reading in pre-created grid-specific
!dust file.
cname=dustfile(1:len_trim(dustfile))

!For reading in grid specific pre-created dust files with grid number
!added to end of file name for potential nested grids.
if(idustloft==22) then
 write(dgrid,'(i1)') ifm
 dname='NRL-dustdata-g'//dgrid
 cname=dustfile(1:len_trim(dustfile))//dgrid
endif

setnxny=0
!Ginoux 1deg dust erodible fraction.
if(idustloft==1) then
 nx_source=360
 ny_source=180
 setnxny=1
!NRL middle east raw data 1's or 0's for dust lofting.
!Reads data, interpolates to RAMS grid, writes out data, and stops.
!This can take a while if RAMS grid is large.
elseif(idustloft==2) then
 nx_source=5600
 ny_source=3400
 setnxny=1
!NRL middle east data precreated, using IDUSTLOFT=2, for grid-specific 
!read-in. This option exists since interpolating the raw data to a specific 
!grid can take a while, and we do not want to wait every time.
elseif(idustloft==22) then
 nx_source=nnxp(ngrid)
 ny_source=nnyp(ngrid)
 setnxny=1
endif

!Allocate if not using IDUSTLOFT=0.
if(setnxny==1) then
 allocate(lat_source(ny_source))
 allocate(lon_source(nx_source))
 allocate(source(nx_source,ny_source))
endif

! Process mainnum handles the file I/O and reading data file
if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
 ! If reading Ginoux(2001) 1-deg dataset
 if(idustloft==1) then
  open(91,file=cname,form='formatted',status='old')
  print*,' '
  print*,'GINOUX 1-DEGREE DUST SOURCE ERODIBLE FRACTION DATA READ IN'
  print*,' '
  do i = 1,nx_source
   do j = 1,ny_source
     read(91,*)lat_source(j),lon_source(i),source(i,j)
   enddo
  enddo
  close(91)
 endif
 ! If reading NRL Middle-East RAW DATA dust database
 if(idustloft==2 .or. idustloft==22) then
  open(91,file=cname,form='formatted',status='old')
  print*,' '
  print*,'NRL 1-KM MIDDLE EAST DUST SOURCE LOCATION DATA READ IN'
  print*,' '
  do i = 1,nx_source
   do j = 1,ny_source
     read(91,*)lat_source(j),lon_source(i),source(i,j)
   enddo
  enddo
  close(91)
  print*,'FINISHED READING NRL 1-KM MIDDLE EAST DUST SOURCE'
 endif
endif

! Process mainnum now has the desired dust data. If there
! are model nodes, broadcast that information to them.
if(setnxny==1 .and. nmachs .gt. 1) then
  CALL broadcast_dustsource (source,nx_source,ny_source)
endif

! Initialize dust fraction on grid
if(print_msg) then
  print*,' '
  print*,'ASSIGNING DUST DATABASE TO MODEL GRID'
  print*,' '
endif

!FOR NRL middle east raw dust data, open file for writing out 
!grid-specfic datafile.
if(idustloft==2) then
 open(92,file=dname,form='formatted',status='unknown')
endif

do i = 1,n2
 do j = 1,n3

  !Default idealized lofting. Every grid cell can potentially
  !loft dust if soil conditions and vegetation permit this.
  if(idustloft==0) then
    dustfrac(i,j) = 1.0
  endif

  !Ginoux(2001) datasource.
  !Look for the source to be used in gridcell (i,j)
  !Source data is in 1 degree increment
  !I-points are incremented westward from -179.5 longitude
  !J-points are incremented northward from -89.5 latitude
  if(idustloft==1) then
    isource = nint(glon(i,j)-(-179.50))+1
    jsource = nint(glat(i,j)-(-089.50))+1
    dustfrac(i,j) = source(isource,jsource)
  endif

  !Use NRL 1km Middle East yes/no data source database to create
  !grid-specific output to read in for later run. This takes a
  !long time since database is very large. Prints grid point
  !assignment to standard output. Trim standard output to create a
  !grid-specific file to read in. Not set up for nested grids yet.
  if(idustloft==2) then
   dustfrac(i,j)=0.0
   if(j>1.and.j<n3.and.i>1.and.i<n2)then
    glat1=(glat(i,j)+glat(i,j-1))/2.0
    glat2=(glat(i,j)+glat(i,j+1))/2.0
    glon1=(glon(i,j)+glon(i-1,j))/2.0
    glon2=(glon(i,j)+glon(i+1,j))/2.0
    sumsource=0.0
    sumnum=0.0
    do ii = 1,nx_source
     do jj = 1,ny_source
       if(lat_source(jj)>glat1 .and. lat_source(jj)<=glat2 .and. &
          lon_source(ii)>glon1 .and. lon_source(ii)<=glon2) then
         sumsource=sumsource+source(ii,jj)
         sumnum=sumnum+1.0
       endif
     enddo
    enddo
    if(sumnum>0)dustfrac(i,j)=sumsource/sumnum
   endif
   write(92,*) j,i,dustfrac(i,j)
  endif

  !If assigning predetermined NRL database grid-specific setup.
  !Assign dustfrac from pre-created NRL 1km Middle East yes/no data 
  !source database grid-specific file.
  if(idustloft==22) then
    dustfrac(i,j) = source(i+mi0(ngrid),j+mj0(ngrid))
  endif

 enddo
enddo

if(idustloft==2) close(92)

if(setnxny==1) then
  deallocate(lat_source,lon_source,source)
endif

if(idustloft==2 .and. ifm==ngrids)then
 if(print_msg) then
  print*,"NRL 1km Middle East dust database written to grid-specific"
  print*," files with prefix 'NRL-dustdata-g'."
  print*," Use new files for new model start."
  print*," Set IDUSTLOFT=22 and DUSTFILE=<path>/NRL-dustdata-g"
  stop
 endif
endif

return
END SUBROUTINE handle_dustsource

!##########################################################################
Subroutine readforcing ()

use ref_sounding
use node_mod
use mem_grid
use netcdf
use io_params, only:iupdsst

implicit none

integer :: ivar, dimid
integer :: ncid, stat, varid, nlev, nt
integer :: ispresent(2)=0

! Forcing file should conform to the DEPHY format for variable names

!Initialize nlev in case it isn't read from the file
nlev = 1
if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
   if (nf90_open(trim(forcingfile), NF90_NOWRITE, ncid) .eq. NF90_NOERR) then

      !Get dimensions and read them
     if(iugforce==2) then
      stat = nf90_inq_dimid(ncid, "lev", dimid)
      if (stat .eq. NF90_NOERR) then
         stat = nf90_inquire_dimension(ncid, dimid, len=nlev)
      else
         stop "Can't find dimension 'lev' in forcing file"
      endif
      allocate(forc_lev(nlev))
      stat = nf90_inq_varid(ncid,"lev",varid)
      stat = nf90_get_var(ncid, varid, forc_lev)
    endif
 
    if(iugforce==2 .or. iupdsst==2) then
      stat = nf90_inq_dimid(ncid, "time", dimid)
      if (stat .eq. NF90_NOERR) then
         stat = nf90_inquire_dimension(ncid, dimid, len=nt)
      else
         stop "Can't find dimension 'time' in forcing file"
      endif
      allocate(forc_time(nt))
      stat = nf90_inq_varid(ncid,"time",varid)
      stat = nf90_get_var(ncid, varid, forc_time)
    endif

      !Now start checking for possible forcing variables and reading them
      stat = nf90_inq_varid(ncid,"ug",varid)
      if (stat == NF90_NOERR) then
         allocate(ug(nlev,nt), vg(nlev,nt))
         stat = nf90_get_var(ncid, varid, ug)
         
         stat = nf90_inq_varid(ncid,"vg",varid)
         stat = nf90_get_var(ncid,varid, vg)
         if (stat .ne. NF90_NOERR) stop "can't read geostrophic wind in forcing file"
      elseif (stat .ne. NF90_NOERR .and. iugforce==2) then
         stop "geostrophic wind not present in forcing file"
      endif

      stat = nf90_inq_varid(ncid,"ts",varid)
      if (stat == NF90_NOERR) then
         allocate(forc_ts(nt))
         stat = nf90_get_var(ncid,varid,forc_ts)
         if (stat .ne. NF90_NOERR) stop "can't read surface temperature forcing"
      elseif (stat .ne. NF90_NOERR .and. iupdsst==2) then
         stop "SSTs not present in forcing file"
      endif
         
   else
      print*,'Error opening forcing file!!!!'
      stop
   endif
endif

! Process mainnum now has the desired data. If there
! are model nodes, broadcast that information to them.
if(nmachs .gt. 1) then
  CALL broadcast_forcing_dims(nlev,nt)

  if (my_rams_num .ne. mainnum) then 
     if (iugforce==2) then
        allocate(forc_lev(nlev),forc_time(nt))
        allocate(ug(nlev,nt),vg(nlev,nt))
     endif
     if (iupdsst==2) then
        if (.not. allocated(forc_time)) allocate(forc_time(nt))
        allocate(forc_ts(nt))
     endif
  endif

  if (iugforce == 2) then     
     CALL broadcast_forcing (ug,nlev,nt)
     CALL broadcast_forcing (vg,nlev,nt)
     CALL broadcast_forcing (forc_lev,nlev,1)
  endif
  if (iupdsst == 2) then
     CALL broadcast_forcing (forc_ts,1,nt)
  endif
  if (iugforce == 2 .or. iupdsst == 2) then
     CALL broadcast_forcing (forc_time,1,nt)
  endif
endif

END SUBROUTINE readforcing
