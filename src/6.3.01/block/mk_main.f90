Program mk_main

use hdf5_utils

implicit none

integer :: i,j,nx,ny,iwlon,islat,iblsize,ibldim,ndims,idims(4) &
  ,nsqx,nsqy,nsx,nsy,lat,lon,npt,jj,ii,nc,numarg,iargc
real :: tres,offlat,offlon
real, allocatable :: datain(:,:),latd(:),lond(:),scr(:) 
character(len=80) :: icent,iyear,imonth,idate,itime
character(len=256) :: source,hpath,fpref,subdir,field_name &
  ,title1,title2,title3,fname,datatype

numarg=iargc()

if (numarg < 1) then
  print*, 'usage: '
  print*, 'executable  infile  datatype  output_directory  prefix'
  print*, '(ie. mk-blkfiles-6.1.0 131231_ 1 output  S)'  
  print*, '(see SST example in directory example_modis_sst)'
  print*, '(file 131231_ is ASCII lat,lon,SST(C) MODIS data)'
  print*, '(files S131231_90S000E.h5 S131231_90S180W.h5 are the output)'
  print*, '(file SHEADER is the header file)'
  print*, ' '
  print*, 'where:'
  print*, '  input_file  - name of input data file'
  print*, '  header_path - directory where header file will be written'
  print*, '  prefix      - filename prefix for header and data files'
  stop ' usage'
endif

! Get the command line arguments
CALL ugetarg (1,source)   ! input file name
CALL ugetarg (2,datatype) ! datatype number indicator
CALL ugetarg (3,hpath)    ! output file path (and location of header file)
CALL ugetarg (4,fpref)    ! output file name prefix for header and blocks

! Set the attributes of the input data. These could be set in an input namelist
!    but let's hardwire these in here for now:

!      Input parameters:
!
!         nx      - Number of x (longitude) gripoints in input array  
!         ny      - Number of y (latitude) gripoints in input array 
!         iwlon   - West longitude of data (-180 to 180 degrees)
!         islat   - South latitude of data (-90 to 90 degrees)
!         tres    - data resolution in degrees
!         iblsize - Size in degrees of output files, they will have the
!                   same number of points in x and y.
!         field_name - name given to field. Only relevant for hdf5 files
!
!     Note that iwlon, islat, and iblsize are all integers.

!Set grid and/or parameter sizes
if(trim(datatype)=='1')then
 nx = 4321
 ny = 2161
 iwlon = -180
 islat = -90
 tres = .08333333
 iblsize = 180
 offlat=0.
 offlon=0.                                             
 field_name='MODIS_9km_sst'
elseif(trim(datatype)=='2')then
 nx = 361
 ny = 181
 iwlon = -180
 islat = -90
 tres = 1.0
 iblsize = 180
 offlat=0.
 offlon=0.                                             
 field_name='ReynoldsV2_1deg_sst'
elseif(trim(datatype)=='3')then
 nx = 10801
 ny = 5401
 iwlon = -180
 islat = -90
 tres = .03333333
 iblsize = 180
 offlat=0.
 offlon=0.
 field_name='Etopo2-Bathymetry-3.72km'
else
 stop 'Back Datatype Command line ID number'
endif

! Allocate array to hold input data and scratch
ibldim=int(float(iblsize)/tres+.001)+1
allocate (scr((ibldim+1)**2))
allocate (datain(nx+1,ny+1))
allocate (latd(ny),lond(nx))

! Read the input data. 
! This will have to be modified depending on actual files....
open(30,status='old',file=source,form='formatted')
4 FORMAT(f8.2,f8.2,E11.3)
DO I=1,nx
 DO J=1,ny
    read(30,4) LATD(J),LOND(I),datain(I,J)
    if(trim(datatype)=='1' .or. trim(datatype)=='2')then
     datain(I,J) = datain(I,J) + 273.16
    endif
 enddo
enddo
close(30)

! Any modifications to input? This example takes sst, overlaps a row in EW-NS,
!    then converts to K
!datain(nx,1:ny) = datain(1,1:ny)
!datain(1:nx,ny) = datain(1:nx,2160)
!datain = datain + 273.16

! Compute number of grid points per block
nsqx=nx/(ibldim-1)
nsqy=ny/(ibldim-1)

! Make sure hpath and subdir have slashes on the end...
nc=len_trim(hpath)  ; if(hpath(nc:nc) /= '/') hpath(nc+1:nc+1)='/'
nc=len_trim(subdir) ; if(subdir(nc:nc) /= '/') subdir(nc+1:nc+1)='/'

! Write the header file for this. If you do more than one time, the header 
! file will have to be manually combined for now.
title3=trim(hpath)//trim(fpref)//'HEADER'
open(29,status='replace',file=title3,form='formatted')
print*, 'making file- ',trim(title3),' ',trim(field_name)
write(29,'(4i5,2f10.6,1X,A)')iblsize,ibldim,islat,iwlon,offlat &
     ,offlon,trim(field_name)
write(29,'(i4)')1
icent='20'
iyear=source(1:2)
imonth=source(3:4)
idate=source(5:6)
itime='00'
write(29,'(1X,a7,1X,a2,a2,3x,a2,3x,a2,3x,a2)') &
  source,icent,iyear,imonth,idate,itime
close(29)

title3=trim(fpref)//source
do nsx=1,nsqx
 do nsy=1,nsqy
   i=(nsx-1)*(ibldim-1)+1
   j=(nsy-1)*(ibldim-1)+1
   lat=islat+(nsy-1)*iblsize
   lon=iwlon+(nsx-1)*iblsize
   if(lat >= 0)then
       write(title1,'(i2.2,a1)')lat,'N'
   else
       write(title1,'(i2.2,a1)')abs(lat),'S'
   endif
   if(lon >= 0)then
       write(title2,'(i3.3,a1)')lon,'E'
   else
       write(title2,'(i3.3,a1)')abs(lon),'W'
   endif
   fname=trim(hpath)//trim(title3)//title1(1:3)//title2(1:4)//'.h5'
   print*, 'making file- ',trim(fname)
   npt=1
   do jj=j,(ibldim-1)+j
    do ii=i,(ibldim-1)+i
       scr(npt)=datain(ii,jj)
       npt=npt+1
    enddo
   enddo
   ! Write the data file block
   ! Write the HDF5 file
   call shdf5_open(fname,'W',1)
   ndims=2 ; idims(1)=ibldim ; idims(2)=ibldim
   call shdf5_orec(ndims,idims,field_name,rvara=scr)
   call shdf5_close()
 enddo
enddo

deallocate (scr,datain,latd,lond)

END PROGRAM mk_main
