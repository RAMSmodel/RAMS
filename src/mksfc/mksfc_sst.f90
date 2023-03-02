!##############################################################################
Subroutine sst_read_dataheader (ifm)

use mem_mksfc
use mem_grid
use io_params
use node_mod

implicit none

integer :: ifm
integer :: itime,nsst,issty,isstm,isstd,issth
character(len=strl1) :: flnm,line,line2
character(len=1) :: dummy
logical :: there
character(len=14) :: totdate_init,totdatem,totdatesst

! Read header file for all sstdata files (all times and locations).  The header
! file contains:
! first line:       geographic block size (degrees), file size, geographic
!                   starting point for the dataset (south pole), and offsets
! second line:      number of data times (NDTIM)
! next NDTIM lines: file prefix, year, month, day, and hour for one data time
! last 3 lines:     comments describing the first, second, and following lines
!                   above

! Construct header file name

flnm=trim(isstfn(ifm))//'HEADER'

print*, 'isstfn(ifm):',trim(flnm)

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' sst data. '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)
if (.not.there) then
   print*,'SSTDATA header file for grid ',ifm,' not there.'
   stop 'sst_read_fileheader-1'
endif

! Read this header file

CALL rams_f_open (25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! read number of data times in dataset

read(25,*) dummy
read(25,*) nsst

if (nsst <= 0) then
   print*, 'No SST input files found with specified prefix or incorrect header'
   close(25)
   stop 'sst_read_fileheader-2'
endif

! read prefix list and times 

CALL date_make_big (iyear1,imonth1,idate1,itime1*100,totdate_init)
CALL date_add_to_big (totdate_init,timmax,'s',totdatem)
   
nvsstf(ifm)=0
do itime = 1,nsst
   read(25,'(A256)') line
   CALL char_strip_var (line,flnm,line2)
   read (line2,*) issty,isstm,isstd,issth
   
   CALL date_make_big (issty,isstm,isstd,issth*100,totdatesst)
   
   ! assumes data in header file is chronologically ordered oldest to newest
   ! if issty=0, assume climo data and do all times
   
   nvsstf(ifm)=nvsstf(ifm)+1
   if(issty /= 0 .and. totdatesst < totdate_init) nvsstf(ifm)=1
   
   vsstfil(nvsstf(ifm),ifm)=trim(isstfn(ifm))//trim(flnm)
   iyearvs(nvsstf(ifm),ifm)=issty
   imonthvs(nvsstf(ifm),ifm)=isstm
   idatevs(nvsstf(ifm),ifm)=isstd
   ihourvs(nvsstf(ifm),ifm)=issth
   
   
   if(issty /= 0 .and. totdatesst > totdatem) exit

! For testing, print out nsstf, ssttime, and vsstfil.
!      print*, 'ifm,itime,nvsstf(ifm)',ifm,itime,nvsstf(ifm)
!      print*, 'vsstfil(itime,ifm)',vsstfil(nvsstf(ifm),ifm)  &
!              ,iyearvs(nvsstf(ifm),ifm),imonthvs(nvsstf(ifm),ifm)  &
!              ,idatevs(nvsstf(ifm),ifm),ihourvs(nvsstf(ifm),ifm)
!      stop

enddo

close(25)

return
END SUBROUTINE sst_read_dataheader

!##############################################################################
Subroutine sstnest (ifm,ivtime)

use mem_mksfc
use mem_grid
use io_params
use node_mod

implicit none

integer :: ifm,icm,ivtime

icm = nxtnest(ifm)

! Initialize SEATP and SEATF in routine sstinit

CALL sstinit (mmxp(ifm),mmyp(ifm),sfcfile_p(ifm)%seatf(1,1))

if (icm >= 1 .and. isstflg(ifm) == 0) then

! Interpolate SEATF from coarser grid
   CALL newgrid (ifm)
   CALL interp_fine_grid_sfcvar (ifm,icm,IFG_SST_INIT)

   nvsstf(ifm) = nvsstf(icm)
   iyearvs (1:nvsstf(ifm),ifm) = iyearvs (1:nvsstf(ifm),icm)
   imonthvs(1:nvsstf(ifm),ifm) = imonthvs(1:nvsstf(ifm),icm)
   idatevs (1:nvsstf(ifm),ifm) = idatevs (1:nvsstf(ifm),icm)
   ihourvs (1:nvsstf(ifm),ifm) = ihourvs (1:nvsstf(ifm),icm)

elseif (isstflg(ifm) == 1) then

! Interpolate SEATF from standard dataset

   CALL geodat (mmxp(ifm),mmyp(ifm),sfcfile_p(ifm)%seatf(1,1)  &
      ,isstfn(ifm),vsstfil(ivtime,ifm),vt2da,vt2db,ifm,'SST')

else

   iyearvs (1,ifm) = iyear1 ; imonthvs(1,ifm) = imonth1
   idatevs (1,ifm) = idate1 ; ihourvs (1,ifm) = ihour1        

endif

! If desired, override current values of SEATF with user-defined
! changes to routine sstinit_user.

CALL sstinit_user (mmxp(ifm),mmyp(ifm),ifm ,sfcfile_p(ifm)%seatf(1,1))

return
END SUBROUTINE sstnest

!##############################################################################
Subroutine sstinit (n2,n3,seatf)

use mem_leaf

implicit none

integer :: n2,n3,i,j
real, dimension(n2,n3) :: seatf

! Fill the SEATF array with a default value of seatmp.  This 
! default is used only when a standard RAMS sst dataset is not used and when 
! no overrides to sea temperature are defined in routine sstinit_user 
! in the file ruser.f90.

do j = 1,n3
   do i = 1,n2
      seatf(i,j) = seatmp
   enddo
enddo 

return
END SUBROUTINE sstinit

!##############################################################################
Subroutine sst_write (ifm,ivt)

use mem_mksfc
use mem_grid
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: ifm,ivt
real :: glatr,glonr
character(len=strl1) :: flnm
character(len=2) :: cgrid
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! Write sst data to sst file for one grid and one time

write(cgrid,'(a1,i1)') 'g',ifm
CALL makefnam (flnm,sstfpfx,0.,iyearvs(ivt,ifm),imonthvs(ivt,ifm) &
      ,idatevs(ivt,ifm),ihourvs (ivt,ifm)*10000,'W',cgrid,'h5')

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'W',iphdf5,h5_fid,iclobber)

! Scalar vars
CALL shdf5_set_hs_select (1,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'year',mem_select,file_select    &
                ,file_chunks,ivars=iyearvn(ivt,ifm))
CALL shdf5_orec (h5_fid,iphdf5,'month',mem_select,file_select   &
                ,file_chunks,ivars=imonthvn(ivt,ifm))
CALL shdf5_orec (h5_fid,iphdf5,'day',mem_select,file_select     &
                ,file_chunks,ivars=idatevn(ivt,ifm))
CALL shdf5_orec (h5_fid,iphdf5,'hour',mem_select,file_select    &
                ,file_chunks,ivars=ihourvn(ivt,ifm))
CALL shdf5_orec (h5_fid,iphdf5,'nx',mem_select,file_select      &
                ,file_chunks,ivars=nnxp(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'ny',mem_select,file_select      &
                ,file_chunks,ivars=nnyp(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'dx',mem_select,file_select      &
                ,file_chunks,rvars=deltaxn(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'polelat',mem_select,file_select &
                ,file_chunks,rvars=polelat)
CALL shdf5_orec (h5_fid,iphdf5,'polelon',mem_select,file_select &
                ,file_chunks,rvars=polelon)
CALL shdf5_orec (h5_fid,iphdf5,'sw_lat',mem_select,file_select  &
                ,file_chunks,rvars=glatr)
CALL shdf5_orec (h5_fid,iphdf5,'sw_lon',mem_select,file_select  &
                ,file_chunks,rvars=glonr)

! Atmos 2D vars
CALL shdf5_set_hs_select (2,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'SEATF',mem_select,file_select &
                ,file_chunks,rvara=sfcfile_p(ifm)%seatf)

CALL shdf5_close (h5_fid)

return
END SUBROUTINE sst_write
