!##############################################################################
Subroutine ndvi_read_dataheader (ifm)

use mem_mksfc
use io_params
use node_mod

implicit none

integer :: ifm
integer :: itime,nc
character(len=strl1) :: flnm,line,line2
character(len=1) :: dummy
logical :: there
integer, external :: lastslash

! Read header file for all ndvidata files (all times and locations).  The header
! file contains:
! first line:       geographic block size (degrees), file size, geographic
!                   starting point for the dataset (south pole), and offsets
! second line:      number of data times (NDTIM)
! next NDTIM lines: file prefix, year, month, day, and hour for one data time
! last 3 lines:     comments describing the first, second, and following lines
!                   above

! Construct header file name

flnm=trim(ndvifn(ifm))//'HEADER'

print*, 'ndvifn(ifm):',trim(flnm)

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' ndvi data. '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)
if (.not.there) then
   print*,'ndvifn data header file for grid ',ifm,' not there.'
   stop 'ndvi_read_fileheader-1'
endif

! Read this header file

CALL rams_f_open (25,flnm,'FORMATTED','OLD','READ',0)
rewind 25

! read number of data times in dataset

read(25,*) dummy
read(25,*) nvndvif(ifm)

if (nvndvif(ifm) <= 0) then
   print*, 'No ndvi input files found with specified prefix or incorrect header'
   close(25)
   stop 'ndvi_read_fileheader-2'
endif

! read prefix list and times 

do itime = 1,nvndvif(ifm)
   read(25,'(A80)') line
   CALL char_strip_var (line,flnm,line2)
   read (line2,*) iyearvn(itime,ifm),imonthvn(itime,ifm)  &
      ,idatevn(itime,ifm),ihourvn(itime,ifm)

   nc = lastslash(ndvifn(ifm))
   vndvifil(itime,ifm)=ndvifn(ifm)(1:nc)//trim(flnm)

! For testing, print out nndvif, ndvitime, and vndvifil.
!   print*, 'ifm,itime,nvndvif(ifm)',ifm,itime,nvndvif(ifm)
!      print*, 'vndvifil(itime,ifm)'  &
!              ,vndvifil(itime,ifm),iyearvs(itime,ifm),imonthvs(itime,ifm)  &
!      ,idatevs(itime,ifm),ihourvs(itime,ifm)
!      stop

enddo

close(25)

return
END SUBROUTINE ndvi_read_dataheader

!##############################################################################
Subroutine ndvinest (ifm,ivtime)

use mem_mksfc
use mem_grid
use mem_leaf
use io_params
use node_mod

implicit none

integer :: ifm,icm,ivtime,i,j,k,ic,jc,ipat

icm = nxtnest(ifm)

! Initialize SEATP and SEATF in routine ndviinit

CALL ndviinit (mmxp(ifm),mmyp(ifm),npatch,sfcfile_p(ifm)%veg_ndvif(1,1,1))

if (icm >= 1 .and. ndviflg(ifm) == 0) then

   ! Assign NDVIF from coarser grid

   do ipat = 2,npatch
      do k = 1,nzg
         do j = 1,mmyp(ifm)
            do i = 1,mmxp(ifm)
               ic = ipm(i,ifm)
               jc = jpm(j,ifm)

               sfcfile_p(ifm)%veg_ndvif(i,j,ipat)  =  &
                  sfcfile_p(icm)%veg_ndvif(ic,jc,ipat)

            enddo
         enddo
      enddo
   enddo
   
   nvndvif(ifm) = nvndvif(icm)
   iyearvn (1:nvndvif(ifm),ifm) = iyearvn (1:nvndvif(ifm),icm)
   imonthvn(1:nvndvif(ifm),ifm) = imonthvn(1:nvndvif(ifm),icm)
   idatevn (1:nvndvif(ifm),ifm) = idatevn (1:nvndvif(ifm),icm)
   ihourvn (1:nvndvif(ifm),ifm) = ihourvn (1:nvndvif(ifm),icm)

elseif (ndviflg(ifm) == 1) then

   ! Assign NDVIF from standard dataset:

   CALL landuse_opqr (mmxp(ifm),mmyp(ifm),nzg,npatch,nvegpat  &
      ,ivegtfn(ifm),isoilfn(ifm) &
      ,ndvifn(ifm),vndvifil(ivtime,ifm)  &
      ,'ndvi' &
      ,sfcfile_p(ifm)%soil_text(1,1,1,1)  &
      ,sfcfile_p(ifm)%patch_area(1,1,1)   &
      ,sfcfile_p(ifm)%leaf_class(1,1,1)   &
      ,sfcfile_p(ifm)%veg_ndvif(1,1,1))
else

   iyearvn (1,ifm) = iyear1 ; imonthvn(1,ifm) = imonth1
   idatevn (1,ifm) = idate1 ; ihourvn (1,ifm) = ihour1        

endif

! If desired, override current values of NDVIP and NDVIF in user-specified
! changes to routine ndviinit_user in the file ruser.f.

CALL ndviinit_user (mmxp(ifm),mmyp(ifm),npatch,ifm  &
   ,sfcfile_p(ifm)%veg_ndvif(1,1,1))

return
END SUBROUTINE ndvinest

!##############################################################################
Subroutine ndviinit (n2,n3,npat,veg_ndvif)

implicit none

integer :: n2,n3,npat,i,j,ipat
real, dimension(n2,n3,npat) :: veg_ndvif

! Fill the veg_ndvif array with a default value of .5.  This default is 
! used only when a standard RAMS ndvi dataset is not used and when no 
! overrides to ndvi values are defined in routine ndviinit_user in the
! file ruser.f.

do j = 1,n3
   do i = 1,n2
      veg_ndvif(i,j,1) = 0.
      veg_ndvif(i,j,2) = .5

      do ipat = 3,npat
         veg_ndvif(i,j,ipat) = veg_ndvif(i,j,2)
      enddo

   enddo
enddo

return
END SUBROUTINE ndviinit

!##############################################################################
Subroutine ndvi_write (ifm,ivt)

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

! Write ndvi data to ndvi file for one grid and one time

write(cgrid,'(a1,i1)') 'g',ifm
CALL makefnam (flnm,ndvifpfx,0.,iyearvn(ivt,ifm),imonthvn(ivt,ifm) &
      ,idatevn(ivt,ifm),ihourvn (ivt,ifm)*10000,'N',cgrid,'h5')

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
CALL shdf5_orec (h5_fid,iphdf5,'npatch',mem_select,file_select  &
                ,file_chunks,ivars=npatch)
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

! Leaf3 2D vars
CALL shdf5_set_hs_select (6,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'VEG_NDVIF',mem_select,file_select &
                ,file_chunks,rvara=sfcfile_p(ifm)%veg_ndvif)
CALL shdf5_close (h5_fid)

return
END SUBROUTINE ndvi_write
