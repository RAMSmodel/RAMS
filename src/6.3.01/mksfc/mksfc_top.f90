!##############################################################################
Subroutine top_read (ifm)

use mem_grid
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: ifm
character(len=strl1) :: flnm
character(len=2) :: cgrid
logical :: there
integer*8 :: h5_fid 
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! read the "top" file

write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

inquire(file=flnm,exist=there)

if(.not.there) then
   print*,'------------------------------------------------'
   print*,'TOP_read: file for grid ',ifm,' not there.'
   print*,'TOP_read: file:',trim(flnm)
   print*,'------------------------------------------------'
   stop 'top_read: no file'
endif

CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Atmos 2D vars
CALL shdf5_set_hs_select (2,'R',ifm,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'TOPT',mem_select     &
                ,file_select,rvara=grid_g(ifm)%topt)
CALL shdf5_irec (h5_fid,iphdf5,'TOPZO',mem_select    &
                ,file_select,rvara=grid_g(ifm)%topzo)

CALL shdf5_close (h5_fid)

return
END SUBROUTINE top_read

!##############################################################################
Subroutine top_check (ifm,ierr)

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

integer :: lc,nsfx,nsfy,nsitoptflg,nsitopsflg,nsiz0flg
real ::  sfdx,sfplat,sfplon,sflat,sflon,stoptenh,stoptwvl  &
   ,sz0max,sz0fact,glatr,glonr

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

lc=len_trim(topfiles)
write(cgrid,'(a1,i1)') 'g',ifm
flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

print*,'------------------------------------------------'
print*,'---> Check grid:',ifm,' top file... '
print*,'--->   Filename:',trim(flnm)

inquire(file=flnm,exist=there)

if(.not.there) then
   ierr = 1
   print*,'TOPfile for grid ',ifm,' not there.'
   print*,'------------------------------------------------'
   return
endif

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Scalar vars
CALL shdf5_set_hs_select (1,'R',ifm,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'nx',mem_select,file_select,ivars=nsfx)
CALL shdf5_irec (h5_fid,iphdf5,'ny',mem_select,file_select,ivars=nsfy)
CALL shdf5_irec (h5_fid,iphdf5,'dx',mem_select,file_select,rvars=sfdx)
CALL shdf5_irec (h5_fid,iphdf5,'polelat',mem_select,file_select,rvars=sfplat)
CALL shdf5_irec (h5_fid,iphdf5,'polelon',mem_select,file_select,rvars=sfplon)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lat',mem_select,file_select,rvars=sflat)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lon',mem_select,file_select,rvars=sflon)
CALL shdf5_irec (h5_fid,iphdf5,'itoptflg',mem_select,file_select,ivars=nsitoptflg)
CALL shdf5_irec (h5_fid,iphdf5,'itopsflg',mem_select,file_select,ivars=nsitopsflg)
CALL shdf5_irec (h5_fid,iphdf5,'toptenh',mem_select,file_select,rvars=stoptenh)
CALL shdf5_irec (h5_fid,iphdf5,'toptwvl',mem_select,file_select,rvars=stoptwvl)
CALL shdf5_irec (h5_fid,iphdf5,'iz0flg',mem_select,file_select,ivars=nsiz0flg)
CALL shdf5_irec (h5_fid,iphdf5,'z0max',mem_select,file_select,rvars=sz0max)
CALL shdf5_irec (h5_fid,iphdf5,'z0fact',mem_select,file_select,rvars=sz0fact)

CALL shdf5_close (h5_fid)


if (nsfx                       .ne. nnxp(ifm)     .or.  &
    nsfy                       .ne. nnyp(ifm)     .or.  &
    abs(sfdx-deltaxn(ifm))     .gt. .001          .or.  &
    abs(sfplat-polelat)        .gt. .001          .or.  &
    abs(sfplon-polelon)        .gt. .001          .or.  &
    abs(sflat-glatr)           .gt. .001          .or.  &
    abs(sflon-glonr)           .gt. .001          .or.  &
    nsitoptflg                 .ne. itoptflg(ifm) .or.  &
    nsitopsflg                 .ne. itopsflg(ifm) .or.  &
    abs(stoptenh-toptenh(ifm)) .gt. .001          .or.  &
    abs(stoptwvl-toptwvl(ifm)) .gt. .001          .or.  &
    nsiz0flg                   .ne. iz0flg(ifm)   .or.  &
    abs(sz0max-z0max(ifm))     .gt. .001          .or.  &
    abs(sz0fact-z0fact)        .gt. .00001) then

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
   print*,'itoptflg:',itoptflg(ifm),nsitoptflg
   print*,'itopsflg:',itopsflg(ifm),nsitopsflg
   print*,'toptenh:',toptenh(ifm),stoptenh
   print*,'toptwvl:',toptwvl(ifm),stoptwvl
   print*,'iz0flg:',iz0flg(ifm),nsiz0flg
   print*,'z0max:',z0max(ifm),sz0max
   print*,'z0fact:',z0fact,sz0fact
   print*,'-------------------'

else

   ierr = 0
   print*,'---> Grid:',ifm,' topography file data okay. '
   print*,'------------------------------------------------'

endif

return
END SUBROUTINE top_check

!##############################################################################
Subroutine toptinit (n2,n3,topt,topzo)

implicit none

integer :: n2,n3,i,j
real, dimension(n2,n3) :: topt,topzo

! Fill the TOPT array with a default value of 0.  This default is used only
! when a standard RAMS topography dataset is not used and when no overrides
! to topography heights are defined in routine toptinit_user in the
! file ruser.f.

do j = 1,n3
   do i = 1,n2
      topt(i,j) = 0.
      topzo(i,j) = .0001
   enddo
enddo

return
END SUBROUTINE toptinit

!##############################################################################
Subroutine top_write (ifm)

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

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

!     write surface characteristics, one file for each grid

write(cgrid,'(a1,i1)') 'g',ifm

flnm=trim(topfiles)//'-S-'//cgrid//'.h5'

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))

CALL shdf5_open (flnm,'W',iphdf5,h5_fid,iclobber)

! Scalar vars
CALL shdf5_set_hs_select (1,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'nx',mem_select,file_select       &
                ,file_chunks,ivars=nnxp(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'ny',mem_select,file_select       &
                ,file_chunks,ivars=nnyp(ifm))
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
CALL shdf5_orec (h5_fid,iphdf5,'itoptflg',mem_select,file_select &
                ,file_chunks,ivars=itoptflg(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'itopsflg',mem_select,file_select &
                ,file_chunks,ivars=itopsflg(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'toptenh',mem_select,file_select  &
                ,file_chunks,rvars=toptenh(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'toptwvl',mem_select,file_select  &
                ,file_chunks,rvars=toptwvl(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'iz0flg',mem_select,file_select   &
                ,file_chunks,ivars=iz0flg(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'z0max',mem_select,file_select    &
                ,file_chunks,rvars=z0max(ifm))
CALL shdf5_orec (h5_fid,iphdf5,'z0fact',mem_select,file_select   &
                ,file_chunks,rvars=z0fact)

! Atmos 2D vars
CALL shdf5_set_hs_select (2,'W',ifm,mem_select,file_select,file_chunks)
CALL shdf5_orec (h5_fid,iphdf5,'TOPT',mem_select,file_select  &
                ,file_chunks,rvara=sfcfile_p(ifm)%topt)
CALL shdf5_orec (h5_fid,iphdf5,'TOPZO',mem_select,file_select &
                ,file_chunks,rvara=sfcfile_p(ifm)%topzo)

CALL shdf5_close (h5_fid)

return
END SUBROUTINE top_write
