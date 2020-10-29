!##############################################################################
integer Function rams_getvar (string,itype,ngrd,a,flnm)

use mem_grid
use node_mod
use an_header
use hdf5_utils
use rcommons

implicit none

real :: a(*)
integer :: itype,ngrd,ni,npts,iword
character(len=*) :: flnm,string
character(len=1) :: cgrid
character(len=strl1) :: flng,errmsg
logical :: there

!File and parallel info (should only do this sequentially, not parallel)
integer :: iphdf5
integer*8 :: h5_fid
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

!Not running in parallel so use default info for sequential
if (nmachs .gt. 1) then
  stop 'Should not be attempting to run REVU in parallel'
else
  iphdf5 = 0
endif

! First see if data file for this grid/time exists...

write(cgrid,'(i1)') ngrd
flng=trim(flnm)//'-g'//cgrid//'.h5'
inquire(file=flng,exist=there)

if(.not.there) then
   errmsg='File not found - '//flng
   CALL error_mess (trim(errmsg))
   rams_getvar=2
   return
endif

! Now search table for actual variable
do ni=1,nvbtab
   if(string == anal_table(ni)%string.and.ngrd == anal_table(ni)%ngrid) then

      npts=anal_table(ni)%nvalues
      itype=anal_table(ni)%idim_type
      iword=anal_table(ni)%npointer

         CALL shdf5_open (trim(flng),'R',iphdf5,h5_fid)
         !call shdf5_info (h5_fid,string,ndims,idims)
         !npts_chk=product(idims(1:ndims)) 
         !if (npts /= npts_chk) then
         !   print*,'No. of points in anal table and in hdf5 file do not match.'
         !   print*,'   anal field:',string
         !   print*,'   anal table:',npts
         !   print*,'   hdf5 file :',idims(1:ndims)
         !   print*,'   hdf5 file :',npts_chk
         !  stop
         !endif
         CALL shdf5_set_hs_select (itype,'R',ngrd,mem_select,file_select,file_chunks)
         CALL shdf5_irec (h5_fid,iphdf5,trim(string),mem_select,file_select,rvara=a)
         CALL shdf5_close (h5_fid)

      rams_getvar=0
      ifound=ifound+1
      return

   endif
enddo

errmsg='Variable not available in this run - '//trim(string)
CALL error_mess (trim(errmsg))
rams_getvar=1
ierr_getvar=1

return
END FUNCTION rams_getvar
