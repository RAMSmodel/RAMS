!##############################################################################
Subroutine sfc_read (ifm)

use mem_grid
use mem_leaf
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: ifm,ip,i,j,k
logical :: there

character(len=strl1) :: flnm
character(len=2) :: cgrid
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

real, dimension(:), allocatable :: r_scratch

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! read the "sfc" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'SFC_read: file for grid ',ifm,' not there.'
   print*,'SFC_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'sfc_read: no file'
endif

CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Leaf3 2D vars
CALL shdf5_set_hs_select (6,'R',ifm,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'PATCH_AREA',mem_select     &
                ,file_select,rvara=leaf_g(ifm)%patch_area)
CALL shdf5_irec (h5_fid,iphdf5,'LEAF_CLASS',mem_select     &
                ,file_select,rvara=leaf_g(ifm)%leaf_class)

! Leaf3 Soil 3D vars
allocate(r_scratch(mmxysp(ifm)))
CALL shdf5_set_hs_select (4,'R',ifm,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'SOIL_TEXT',mem_select &
                ,file_select,rvara=r_scratch)
CALL unarrange_p (mmxp(ifm),mmyp(ifm),nzg,npatch,r_scratch &
                 ,leaf_g(ifm)%soil_text)
deallocate(r_scratch)

CALL shdf5_close (h5_fid)

!Saleeby(2014): Set grid boundaries as in leaf3 for surface characteristics
!Do this here so that boundary assignments are made before the first
!analysis write and first timestep.
do ip = 1,npatch
  do j = 1,mmyp(ifm)
     leaf_g(ifm)%leaf_class(1,j,ip)=leaf_g(ifm)%leaf_class(2,j,ip)
     leaf_g(ifm)%patch_area(1,j,ip)=leaf_g(ifm)%patch_area(2,j,ip)
     do k = 1,nzg
      leaf_g(ifm)%soil_text(k,1,j,ip)=leaf_g(ifm)%soil_text(k,2,j,ip)
     enddo
     leaf_g(ifm)%leaf_class(mmxp(ifm),j,ip)=leaf_g(ifm)%leaf_class(mmxp(ifm)-1,j,ip)
     leaf_g(ifm)%patch_area(mmxp(ifm),j,ip)=leaf_g(ifm)%patch_area(mmxp(ifm)-1,j,ip)
     do k = 1,nzg
      leaf_g(ifm)%soil_text(k,mmxp(ifm),j,ip)=leaf_g(ifm)%soil_text(k,mmxp(ifm)-1,j,ip)
     enddo
  enddo
  if (jdim == 1) then
   do i = 1,mmxp(ifm)
      leaf_g(ifm)%leaf_class(i,1,ip)=leaf_g(ifm)%leaf_class(i,2,ip)
      leaf_g(ifm)%patch_area(i,1,ip)=leaf_g(ifm)%patch_area(i,2,ip)
      do k = 1,nzg
       leaf_g(ifm)%soil_text(k,i,1,ip)=leaf_g(ifm)%soil_text(k,i,2,ip)
      enddo
      leaf_g(ifm)%leaf_class(i,mmyp(ifm),ip)=leaf_g(ifm)%leaf_class(i,mmyp(ifm)-1,ip)
      leaf_g(ifm)%patch_area(i,mmyp(ifm),ip)=leaf_g(ifm)%patch_area(i,mmyp(ifm)-1,ip)
      do k = 1,nzg
       leaf_g(ifm)%soil_text(k,i,mmyp(ifm),ip)=leaf_g(ifm)%soil_text(k,i,mmyp(ifm)-1,ip)
      enddo
   enddo
  endif
enddo

return
END SUBROUTINE sfc_read

!##############################################################################
Subroutine sfc_check (ifm,ierr)

! This routine checks for the existence of a surface file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ifileok is returned with a value of 0.  If the file
! exists and is ok, ifileok is returned with a value of 1.

use mem_grid
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: ifm,ierr
integer :: lc,nsfx,nsfy,nsfzg ,nsivegtflg,nsisoilflg,nspatch
real ::  sfdx,sfplat,sfplon,sflat,sflon,glatr,glonr

character(len=strl1) :: flnm
character(len=2) :: cgrid
logical there
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

lc=len_trim(sfcfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' sfc file... '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'SFCfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Scalar variables (single values, no array)
CALL shdf5_set_hs_select (1,'R',ifm,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'nx',mem_select,file_select,ivars=nsfx)
CALL shdf5_irec (h5_fid,iphdf5,'ny',mem_select,file_select,ivars=nsfy)
CALL shdf5_irec (h5_fid,iphdf5,'nzg',mem_select,file_select,ivars=nsfzg)
CALL shdf5_irec (h5_fid,iphdf5,'npatch',mem_select,file_select,ivars=nspatch)
CALL shdf5_irec (h5_fid,iphdf5,'dx',mem_select,file_select,rvars=sfdx)
CALL shdf5_irec (h5_fid,iphdf5,'polelat',mem_select,file_select,rvars=sfplat)
CALL shdf5_irec (h5_fid,iphdf5,'polelon',mem_select,file_select,rvars=sfplon)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lat',mem_select,file_select,rvars=sflat)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lon',mem_select,file_select,rvars=sflon)
CALL shdf5_irec (h5_fid,iphdf5,'ivegtflg',mem_select,file_select,ivars=nsivegtflg)
CALL shdf5_irec (h5_fid,iphdf5,'isoilflg',mem_select,file_select,ivars=nsisoilflg)

CALL shdf5_close (h5_fid)


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    nsfzg                      .ne. nzg           .or.  &
    nspatch                    .ne. npatch        .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfplat-polelat)        .gt. .001          .or.  &
    abs(sfplon-polelon)        .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsivegtflg                 .ne. ivegtflg(ifm) .or.  &
    nsisoilflg                 .ne. isoilflg(ifm) ) then

   ierr = 1

   print*,'SFCfile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'deltax:',deltaxn(ifm),sfdx
   print*,'polelat:',polelat,sfplat
   print*,'polelon:',polelon,sfplon
   print*,'SW lat:',glatr,sflat
   print*,'SW lon:',glonr,sflon
   print*,'ivegtflg:',ivegtflg(ifm),nsivegtflg
   print*,'isoilflg:',isoilflg(ifm),nsisoilflg
   print*,'-------------------'

else

   ierr = 0
   print*,'---> Grid:',ifm,' surface file data okay. '
   print*,'------------------------------------------------'

endif

return
END SUBROUTINE sfc_check

!##############################################################################
Subroutine sfc_write (ifm)

use mem_mksfc
use mem_grid
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: ifm
real :: glatr,glonr
character(len=strl1) :: flnm
character(len=2) :: cgrid
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

real, dimension(:), allocatable :: r_scratch

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

!     write surface characteristics, one file for each grid


write(cgrid,'(a1,i1)') 'g',ifm

flnm=trim(sfcfiles)//'-S-'//cgrid//'.h5'

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))


CALL shdf5_open (flnm,'W',iphdf5,h5_fid,iclobber)

! Scalar variables
CALL shdf5_set_hs_select (1,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'nx',mem_select,file_select       &
                ,file_chunks,ivars=nnxp(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'ny',mem_select,file_select       &
                ,file_chunks,ivars=nnyp(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'nzg',mem_select,file_select      &
                ,file_chunks,ivars=nzg)
CALL shdf5_orec (h5_fid,iphdf5,'npatch',mem_select,file_select   &
                ,file_chunks,ivars=npatch)
CALL shdf5_orec (h5_fid,iphdf5,'dx',mem_select,file_select       &
                ,file_chunks,rvars=deltaxn(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'polelat',mem_select,file_select  &
                ,file_chunks,rvars=polelat)
CALL shdf5_orec (h5_fid,iphdf5,'polelon',mem_select,file_select  &
                ,file_chunks,rvars=polelon)
CALL shdf5_orec (h5_fid,iphdf5,'sw_lat',mem_select,file_select   &
                ,file_chunks,rvars=glatr)
CALL shdf5_orec (h5_fid,iphdf5,'sw_lon',mem_select,file_select   &
                ,file_chunks,rvars=glonr)
CALL shdf5_orec (h5_fid,iphdf5,'ivegtflg',mem_select,file_select &
                ,file_chunks,ivars=ivegtflg(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'isoilflg',mem_select,file_select &
                ,file_chunks,ivars=isoilflg(ifm))

! Leaf3 2D vars
CALL shdf5_set_hs_select (6,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'PATCH_AREA',mem_select,file_select &
                ,file_chunks,rvara=sfcfile_p(ifm)%patch_area)
CALL shdf5_orec (h5_fid,iphdf5,'LEAF_CLASS',mem_select,file_select &
                ,file_chunks,rvara=sfcfile_p(ifm)%leaf_class)

! Leaf3 Soil 3D vars
allocate(r_scratch(mmxysp(ifm)))
CALL rearrange_p (mmxp(ifm),mmyp(ifm),nzg,npatch,sfcfile_p(ifm)%soil_text &
                 ,r_scratch)
CALL shdf5_set_hs_select (4,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'SOIL_TEXT',mem_select,file_select &
                ,file_chunks,rvara=r_scratch)
deallocate(r_scratch)

CALL shdf5_close (h5_fid)

return
END SUBROUTINE sfc_write
