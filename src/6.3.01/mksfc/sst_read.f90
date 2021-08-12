!##############################################################################
Subroutine sst_read (runflag,ifm,ierr)

use mem_grid
use mem_leaf
use io_params
use node_mod

implicit none

integer :: runflag,ifm,ierr,nf,ng
character(len=14), save  :: totdate_start,totdate_init
character(len=14) :: totdate,totdatem
real(kind=8) :: secs_init,secs1,secs2

ierr = 0

if (runflag == 1 .or. runflag == 2) then   ! Initialization(1) or file check(2)

   ! Inventory all sst surface files. 
   CALL sst_file_inv (sstfpfx,ierr)
   if(ierr == 1) then
      if(runflag == 2) return
      if(runflag == 1) stop 'sst_read: error on init'
   endif

   ! Find init and start date
   CALL date_make_big (iyear1,imonth1,idate1,itime1*100  &
                ,totdate_init)
   if (trim(runtype) == 'HISTORY') then
      CALL date_add_to_big (totdate_init,time,'s',totdate_start)
   elseif (trim(runtype) == 'INITIAL'   .or. trim(runtype) == 'ERROR' .or. &
           trim(runtype) == 'MAKEVFILE' .or. trim(runtype) == 'MAKEHFILE') then
      totdate_start=totdate_init
   endif

   ! Do some checks on times
   do ng = 1, ngrids
      if(itotdate_sst(1,ng) > totdate_start .and. isstcycdata == 0) then
         print*, 'sst_read: initial sst file time for grid ',ng  &
               ,' later than beginning of run'
         ierr = 1
      endif
         
      if(iupdsst == 1 .and. nsstfiles(ng) == 1) then
         print*, 'sst_read: updating SST values but only one SST file'  &
                  ,' for grid',ng
         ierr = 1
      endif
      
      CALL date_add_to_big (totdate_init,timmax,'s',totdatem)
      if(iupdsst == 1 .and. isstcycdata == 0 .and. &
               itotdate_sst(nsstfiles(ng),ng) < totdatem) then
         print*, 'sst_read: final sst file time for grid ',ng  &
                  ,'earlier than end of run - making new sst files'
         ierr = 1
      endif
   enddo
   
   ! Return if errors
   if(ierr == 1 .and. runflag == 1) return
   
   ! If we are only checking, we're done.
   if(runflag == 2) return 

   do ng = 1, ngrids
   
      ! Change the sst file dates to the current year. We will increment
      !   when we need to read a future file.
      if(isstcycdata == 1) then
         do nf=1,nsstfiles(ng)
            itotdate_sst(nf,ng)(1:4)=totdate_start(1:4)
         enddo
      endif
         
      ! Find past time file. The files are ordered in time, so we only 
      !    need to find when start time is greater than a file.
      isstflp(ng)=0
      do nf=nsstfiles(ng),1,-1
         if(totdate_start >= itotdate_sst(nf,ng) ) then
            isstflp(ng)=nf
            exit
         endif
      enddo
      
      isstflf(ng)=isstflp(ng)+1
   
      ! If we are cyclic, we possibly didn't find a start time later than
      !   a file time. I think it is safe to assume that we will use the 
      !   last file for the start file. Also see if future time will be on 
      !   next cycle.
      if(isstcycdata == 1 ) then
         if(isstflp(ng) == 0) isstflp(ng) = nsstfiles(ng)
         if(isstflf(ng) > nsstfiles(ng)) isstflf(ng)=1
      endif
   
      if(print_msg) print*,'sst starting at file:',ng,isstflp(ng),isstflf(ng)

      ! Read past/future time sst field

      CALL sst_update (0,isstflp(ng))
      
      if(iupdsst == 1 ) then
         ! Read future time sst field if updating
         CALL sst_update (1,isstflf(ng))
     
         ! Compute times as number of seconds past 1 Jan 1900
         CALL date_abs_secs (totdate_init,secs_init)
         
         ! Get model time of past file
         totdatem=itotdate_sst(isstflp(ng),ng)
         if(isstcyclic == 1) then
            ! If month of past file > current month, subtract a year
            if(totdatem(5:6) > totdate_start(5:6))   &
               CALL date_add_to_big (totdatem,-365.,'d',totdatem)
         endif
         CALL date_abs_secs (totdatem,secs1)
   
         totdatem=itotdate_sst(isstflf(ng),ng)
         if(isstcyclic == 1) then
            if(totdatem < totdate_start) then
               ! Future file is in next year. Update all file names for a new year
               CALL date_add_to_big (totdatem, 365.,'d',totdatem)
               do nf=1,nsstfiles(ng)
                  itotdate_sst(nf,ng)(1:4)=totdatem(1:4)
               enddo
            endif
         endif
         CALL date_abs_secs (totdatem,secs2)
      
         ssttime1(ng) = secs1 - secs_init
         ssttime2(ng) = secs2 - secs_init
      else 
         leaf_g(ng)%seatp(1:mmxp(ng),1:mmyp(ng))=  &
         leaf_g(ng)%seatf(1:mmxp(ng),1:mmyp(ng))
         ssttime1(ng) = 0.
         ssttime2(ng) = 0.
      endif
   enddo
   
   return

elseif (runflag == 3) then   ! Runtime file increment
   
   if ( time >= ssttime2(ifm) ) then
   
      ! Update sst fields
      isstflp(ifm) = isstflf(ifm)
      isstflf(ifm) = isstflp(ifm) + 1
      
      ! Compute times as number of seconds past 1 Jan 1900
      !   If cyclic, modify file date 
      CALL date_abs_secs (totdate_init,secs_init)
      CALL date_add_to_big (totdate_init,time,'s',totdate)

      totdatem=itotdate_sst(isstflp(ifm),ifm)
      CALL date_abs_secs (totdatem,secs1)
      
      ! Need to deal with cyclic. 
      if(isstcyclic == 1 .and. isstflf(ifm) > nsstfiles(ifm)) then
         isstflf(ifm)=1
         ! Update all file names for a new year
         totdatem=itotdate_sst(isstflf(ifm),ifm)
         CALL date_add_to_big (totdatem, 365.,'d',totdatem)
         do nf=1,nsstfiles(ifm)
            itotdate_sst(nf,ifm)(1:4)=totdatem(1:4)
         enddo
      endif
         
      totdatem=itotdate_sst(isstflf(ifm),ifm)
      CALL date_abs_secs (totdatem,secs2)
      
      ssttime1(ifm) = secs1 - secs_init
      ssttime2(ifm) = secs2 - secs_init
      
      ! Finally read the actual field           
      CALL sst_update (1,isstflf(ifm))

      print*,'Switched sst files:',ssttime1(ifm),ssttime2(ifm)  &
               ,secs1,secs2,secs_init
   
   else
   
      return
   
   endif
   
endif

return
END SUBROUTINE sst_read

!##############################################################################
Subroutine sst_file_inv (sfilin,ierr)

use mem_grid
use io_params

implicit none

character(len=*) :: sfilin
integer :: ierr,nf,lnf,nftot,ng,inyear,inmonth,indate,inhour
character(len=strl1) :: fnames(maxsstfiles)
character(len=14) :: itotdate
character(len=1) :: cgrid
real(kind=8) :: secs_init

ierr=0

! Get abs seconds of run start

CALL date_abs_secs2 (iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through sst files and make inventory. We unfortunately have to do this
!   for all grids.

isstcyclic=0
isstcycdata=0

do ng = 1, ngrids
   
   nftot = -1
   write(cgrid,'(i1)') ng
   CALL rams_filelist (fnames,trim(sfilin)//'-W-*-g'//cgrid//'.h5',nftot)
   
   if(nftot <= 0) then
      print*,'No sst files for grid '//cgrid
      ierr=1
      return
   endif   
   
   if(nftot > maxsstfiles) then
      print*,'too many sst files'
      stop 'sst_file_inv: lots_of_sst'
   endif
   

   ! We need to see if the data is to be considered "cyclic". Count how many files
   !  have the year 0000. If all of them do not, we probably have an error.
   !  isstcycdata = 1 when data is cyclic
   !  isstcyclic =1 when we are have cyclic data and will be updating in time

   if(print_msg)print*,'-----------------------------------------------------'
   nsstfiles(ng)=0
   do nf=1,nftot
      lnf=len_trim(fnames(nf))
      read(fnames(nf)(lnf-22:lnf-6),20) inyear,inmonth,indate,inhour
      20 format(i4,1x,i2,1x,i2,1x,i6)

      ! Check the file headers
      CALL sst_check_header (ng,fnames(nf),ierr)
      if(ierr == 1) return
   
      if(inyear == 0) isstcycdata=isstcycdata+1

      CALL date_make_big (inyear,inmonth,indate,inhour,itotdate)

      nsstfiles(ng)=nsstfiles(ng)+1
      fnames_sst(nsstfiles(ng),ng)=fnames(nf)
      itotdate_sst(nsstfiles(ng),ng)=itotdate

   enddo

   CALL rams_dintsort (nsstfiles(ng),itotdate_sst(1,ng),fnames_sst(1,ng))


   !  start printing section
   !--------------------------------------------------------------

   if(print_msg) then
    print*,'-------------------------------------------------------------'
    print*,'-----------  SST Input File Inventory: Grid '//cgrid
    print*,'-------------------------------------------------------------'
    do nf=1,nsstfiles(ng)
       print*,  itotdate_sst(nf,ng),'   ',trim(fnames_sst(nf,ng))
    enddo
    print*,'------------------------------------------------------'
   endif

enddo

! Check the cyclic data condition. 
!   WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS
if(isstcycdata > 0) then
   if(isstcycdata /= sum(nsstfiles(1:ngrids))) then
      print*, 'All sst surface files do not have year 0000'
      print*, 'This confuses the gods and can not occur.'
      stop 'sst_inv'
   endif
   isstcycdata=1
else
   isstcycdata=0
endif


! Set the main cyclic flag. Only relevant if we are updating in time.
if(iupdsst == 1 .and. isstcycdata == 1) isstcyclic=1

return
END SUBROUTINE sst_file_inv

!##############################################################################
Subroutine sst_update (iswap,nfile)

use mem_grid
use mem_leaf
use node_mod
use io_params
use hdf5_utils
use node_mod

implicit none

integer :: iswap,nfile

integer :: ng,nc
character(len=1) :: cgrid
character(len=strl1) :: flnm
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! Put new fields into future arrays. If iswap == 1, 
!     swap future into past first

if (iswap == 1) then
   do ng=1,ngrids
      leaf_g(ng)%seatp(1:mmxp(ng),1:mmyp(ng))=  &
         leaf_g(ng)%seatf(1:mmxp(ng),1:mmyp(ng))
   enddo
endif


! Open the input file for each grid and read field.
do ng=1,ngrids
   ! Contruct file name for this grid
   write(cgrid, '(i1)') ng
   flnm=fnames_sst(nfile,ng)
   nc=len_trim(flnm)-3
   flnm(nc:nc)=cgrid

   CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

   ! Atmos 2D vars
   CALL shdf5_set_hs_select (2,'R',ng,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'SEATF',mem_select &
                   ,file_select,rvara=leaf_g(ng)%seatf)

   CALL shdf5_close (h5_fid)

enddo

return
END SUBROUTINE sst_update

!##############################################################################
Subroutine sst_check_header (ifm,flnm,ierr)

! This routine checks for the existence of an sst file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ierr is returned with a value of 1.  If the file
! exists and is ok, ierr is returned with a value of 0.

use mem_grid
use node_mod
use hdf5_utils

implicit none

integer :: ifm,ierr,nsfx,nsfy
character(len=*) :: flnm
real :: sfdx,sfplat,sfplon,sflat,sflon,glatr,glonr
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

ierr = 0

if(print_msg) &
 print '(a,i3,1x,a,a,$)',' ---> Check grid:',ifm,' sst filename:',trim(flnm)


CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))


CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Scalar vars
CALL shdf5_set_hs_select (1,'R',1,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'nx',mem_select,file_select,ivars=nsfx)
CALL shdf5_irec (h5_fid,iphdf5,'ny',mem_select,file_select,ivars=nsfy)
CALL shdf5_irec (h5_fid,iphdf5,'dx',mem_select,file_select,rvars=sfdx)
CALL shdf5_irec (h5_fid,iphdf5,'polelat',mem_select,file_select,rvars=sfplat)
CALL shdf5_irec (h5_fid,iphdf5,'polelon',mem_select,file_select,rvars=sfplon)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lat',mem_select,file_select,rvars=sflat)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lon',mem_select,file_select,rvars=sflon)

CALL shdf5_close (h5_fid)

if (nsfx                   .ne. nnxp(ifm) .or.  &
    nsfy                   .ne. nnyp(ifm) .or.  &
    abs(sfdx-deltaxn(ifm)) .gt. .001      .or.  &
    abs(sfplat-polelat)    .gt. .001      .or.  &
    abs(sfplon-polelon)    .gt. .001      .or.  &
    abs(sflat-glatr)       .gt. .001      .or.  &
    abs(sflon-glonr)       .gt. .001) then

   ierr = 1

   print*
   print*,'SSTfile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'deltax:',deltaxn(ifm),sfdx
   print*,'polelat:',polelat,sfplat
   print*,'polelon:',polelon,sfplon
   print*,'SW lat:',glatr,sflat
   print*,'SW lon:',glonr,sflon
   print*,'-------------------'

else

   ierr = 0
   if(print_msg) print*,' ---> good '

endif

return
END SUBROUTINE sst_check_header

