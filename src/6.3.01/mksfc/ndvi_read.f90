!##############################################################################
Subroutine ndvi_read (runflag,ifm,ierr)

use mem_grid
use mem_leaf
use io_params
use node_mod

implicit none

integer :: runflag,ifm,ierr

character(len=14), save  :: totdate_start,totdate_init
character(len=14) :: totdate,totdatem
integer :: nf,ng,i,j,ip
real :: timefac_ndvi
real(kind=8) :: secs_init,secs1,secs2

ierr = 0

if (runflag == 1 .or. runflag == 2) then   ! Initialization(1) or file check(2)

   ! Inventory all ndvi surface files. 
   CALL ndvi_file_inv (ndvifpfx,ierr)
   if(ierr == 1) then
      if(runflag == 2) return
      if(runflag == 1) stop 'ndvi_read: error on init'
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
      if(itotdate_ndvi(1,ng) > totdate_start .and. indvicycdata == 0) then
         print*, 'initial ndvi file time for grid ',ng  &
               ,' later than beginning of run'
         ierr = 1
      endif
         
      if(iupdndvi == 1 .and. nndvifiles(ng) == 1) then
         print*, 'updating ndvi values but only one ndvi file'  &
                  ,' for grid',ng
         ierr = 1
      endif
      
      CALL date_add_to_big (totdate_init,timmax,'s',totdatem)
      if(iupdndvi == 1 .and. indvicycdata == 0 .and. &
               itotdate_ndvi(nndvifiles(ng),ng) < totdatem) then
         print*, 'ndvi_read: final ndvi file time for grid ',ng  &
                  ,'earlier than end of run - making new ndvi files'
         ierr = 1
      endif
   enddo
   
   ! Return if errors
   if(ierr == 1 .and. runflag == 1) return
   
   ! If we are only checking, we're done.
   if(runflag == 2) return 

   do ng = 1, ngrids
   
      ! Change the ndvi file dates to the current year. We will increment
      !   when we need to read a future file.
      if(indvicyclic == 1) then
         do nf=1,nndvifiles(ng)
            itotdate_ndvi(nf,ng)(1:4)=totdate_start(1:4)
         enddo
      endif
         
      ! Find past time file. The files are ordered in time, so we only 
      !    need to find when start time is greater than a file.
      indviflp(ng)=0
      do nf=nndvifiles(ng),1,-1
         if(totdate_start >= itotdate_ndvi(nf,ng) ) then
            indviflp(ng)=nf
            exit
         endif
      enddo
      
      indviflf(ng)=indviflp(ng)+1
   
      ! If we are cyclic, we possibly didn't find a start time later than
      !   a file time. I think it is safe to assume that we will use the 
      !   last file for the start file. Also see if future time will be on 
      !   next cycle.
      if(indvicycdata == 1 ) then
         if(indviflp(ng) == 0) indviflp(ng) = nndvifiles(ng)
         if(indviflf(ng) > nndvifiles(ng)) indviflf(ng)=1
      endif
   
      if(print_msg) print*,'ndvi starting at file:',ng,indviflp(ng),indviflf(ng)

      ! Read past/future time ndvi field

      CALL ndvi_update (0,indviflp(ng))
      
      if(iupdndvi == 1) then
         ! Read future time ndvi field if updating
         CALL ndvi_update (1,indviflf(ng))

         ! Compute times as number of seconds past 1 Jan 1900
         CALL date_abs_secs (totdate_init,secs_init)
         
         ! Get model time of past file
         totdatem=itotdate_ndvi(indviflp(ng),ng)
         if(indvicyclic == 1) then
            ! If month of past file > current month, subtract a year
            if(totdatem(5:6) > totdate_start(5:6))   &
               CALL date_add_to_big (totdatem,-365.,'d',totdatem)
         endif
         CALL date_abs_secs (totdatem,secs1)
   
         totdatem=itotdate_ndvi(indviflf(ng),ng)
         if(indvicyclic == 1) then
            if(totdatem < totdate_start) then
               ! Future file is in next year. Update all file names for a new year
               CALL date_add_to_big (totdatem, 365.,'d',totdatem)
               do nf=1,nndvifiles(ng)
                  itotdate_ndvi(nf,ng)(1:4)=totdatem(1:4)
               enddo
            endif
         endif
         CALL date_abs_secs (totdatem,secs2)
      
         ndvitime1(ng) = secs1 - secs_init
         ndvitime2(ng) = secs2 - secs_init
      else
         do ip=1,npatch
            leaf_g(ng)%veg_ndvip(1:mmxp(ng),1:mmyp(ng),ip)=  &
            leaf_g(ng)%veg_ndvif(1:mmxp(ng),1:mmyp(ng),ip)
         enddo
         ndvitime1(ng) = 0.
         ndvitime2(ng) = 0.
      endif

      ! Fill "current" time ndvi values for INITIAL start
      ! History restart will already have the current value stored
      if (trim(runtype) /= 'HISTORY')then
       if (iupdndvi == 0) then
         timefac_ndvi = 0.
       else
         timefac_ndvi = (time - ndvitime1(ng))   &
                      / (ndvitime2(ng) - ndvitime1(ng))
       endif

       do ip=1,npatch
         do j=1,mmyp(ng)
            do i=1,mmxp(ng)
               leaf_g(ng)%veg_ndvic(i,j,ip) = leaf_g(ng)%veg_ndvip(i,j,ip)  &
                  + (leaf_g(ng)%veg_ndvif(i,j,ip)  &
                  -  leaf_g(ng)%veg_ndvip(i,j,ip)) * timefac_ndvi
                 ! print*,'ndvi:',i,j,ip,leaf_g(ng)%veg_ndvic(i,j,ip)
            enddo
         enddo
       enddo     
      endif

   enddo
   
   return

elseif (runflag == 3) then   ! Runtime file increment
   
   CALL date_add_to_big (totdate_init,time,'s',totdate)

   if ( time >= ndvitime2(ifm) ) then
   
      ! Update ndvi fields
      indviflp(ifm) = indviflf(ifm)
      indviflf(ifm) = indviflp(ifm) + 1
      
      ! Compute times as number of seconds past 1 Jan 1900
      !   If cyclic, modify file date 
      CALL date_abs_secs (totdate_init,secs_init)
      CALL date_add_to_big (totdate_init,time,'s',totdate)

      totdatem=itotdate_ndvi(indviflp(ifm),ifm)
      CALL date_abs_secs (totdatem,secs1)
      
      ! Need to deal with cyclic. 
      if(indvicyclic == 1 .and. indviflf(ifm) > nndvifiles(ifm)) then
         indviflf(ifm)=1
         ! Update all file names for a new year
         totdatem=itotdate_ndvi(indviflf(ifm),ifm)
         CALL date_add_to_big (totdatem, 365.,'d',totdatem)
         do nf=1,nndvifiles(ifm)
            itotdate_ndvi(nf,ifm)(1:4)=totdatem(1:4)
         enddo
      endif
         
      totdatem=itotdate_ndvi(indviflf(ifm),ifm)
      CALL date_abs_secs (totdatem,secs2)
      
      ndvitime1(ifm) = secs1 - secs_init
      ndvitime2(ifm) = secs2 - secs_init
   
      ! Finally read the actual field           
      CALL ndvi_update (1,indviflf(ifm))

      print*,'Switched ndvi files:',ndvitime1(ifm),ndvitime2(ifm)  &
               ,secs1,secs2,secs_init

   else
   
      return
   
   endif
   
endif

return
END SUBROUTINE ndvi_read

!##############################################################################
Subroutine ndvi_file_inv (sfilin,ierr)

use mem_grid
use io_params

implicit none

character(len=*) :: sfilin
integer :: ierr,nf,lnf,nftot,ng,inyear,inmonth,indate,inhour
character(len=strl1) :: fnames(maxndvifiles)
character(len=14) :: itotdate
character(len=1) :: cgrid
real(kind=8) :: secs_init

ierr=0

! Get abs seconds of run start

CALL date_abs_secs2 (iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through ndvi files and make inventory. We unfortunately have to do this
!   for all grids.

indvicyclic=0
indvicycdata=0

do ng = 1, ngrids
   
   nftot = -1
   write(cgrid,'(i1)') ng
   CALL rams_filelist (fnames,trim(sfilin)//'-N-*-g'//cgrid//'.h5',nftot)
   
   if(nftot <= 0) then
      print*,'No ndvi files for grid '//cgrid
      ierr=1
      return
   endif   
   
   if(nftot > maxndvifiles) then
      print*,'too many ndvi files'
      stop 'ndvi_file_inv: lots_of_ndvi'
   endif
   

   ! We need to see if the data is to be considered "cyclic". Count how many files
   !  have the year 0000. If all of them do not, we probably have an error.
   !  indvicycdata = 1 when data is cyclic
   !  indvicyclic =1 when we are have cyclic data and will be updating in time

   nndvifiles(ng)=0
   do nf=1,nftot
      lnf=len_trim(fnames(nf))
      read(fnames(nf)(lnf-22:lnf-6),20) inyear,inmonth,indate,inhour
      20 format(i4,1x,i2,1x,i2,1x,i6)

      ! Check the file headers
      CALL ndvi_check_header (ng,fnames(nf),ierr)
      if(ierr == 1) return
   
      if(inyear == 0) indvicycdata=indvicycdata+1

      CALL date_make_big (inyear,inmonth,indate,inhour,itotdate)

      nndvifiles(ng)=nndvifiles(ng)+1
      fnames_ndvi(nndvifiles(ng),ng)=fnames(nf)
      itotdate_ndvi(nndvifiles(ng),ng)=itotdate

   enddo

   CALL rams_dintsort (nndvifiles(ng),itotdate_ndvi(1,ng),fnames_ndvi(1,ng))


   !  start printing section
   !--------------------------------------------------------------

   if(print_msg) then
    print*,'-------------------------------------------------------------'
    print*,'-----------  NDVI Input File Inventory: Grid '//cgrid
    print*,'-------------------------------------------------------------'
    do nf=1,nndvifiles(ng)
       print*,  itotdate_ndvi(nf,ng),'   ',trim(fnames_ndvi(nf,ng))
    enddo
    print*,'------------------------------------------------------'
   endif

enddo

! Check the cyclic condition. Only relevant if we are updating in time.
!   WE ARE ONLY ALLOWING CYCLIC ON ALL GRIDS
if(indvicycdata > 0) then
   if(indvicycdata /= sum(nndvifiles(1:ngrids))) then
      print*, 'All ndvi surface files do not have year 0000'
      print*, 'This confuses the gods and can not occur.'
      stop 'ndvi_inv'
   endif
   indvicycdata=1
else
   indvicycdata=0
endif


! Set the main cyclic flag. Only relevant if we are updating in time.
if(iupdndvi == 1 .and. indvicycdata == 1) indvicyclic=1

return
END SUBROUTINE ndvi_file_inv

!##############################################################################
Subroutine ndvi_update (iswap,nfile)

use mem_leaf
use mem_grid
use node_mod
use io_params
use hdf5_utils

implicit none

integer :: iswap,nfile,i,j
integer :: ng,nc,ip
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
      do ip = 1,npatch
         leaf_g(ng)%veg_ndvip(1:mmxp(ng),1:mmyp(ng),ip)=  &
            leaf_g(ng)%veg_ndvif(1:mmxp(ng),1:mmyp(ng),ip)
      enddo
   enddo
endif


! Open the input file for each grid and read field.
do ng=1,ngrids
   ! Contruct file name for this grid
   write(cgrid, '(i1)') ng
   flnm=fnames_ndvi(nfile,ng)
   nc=len_trim(flnm)-3
   flnm(nc:nc)=cgrid

   CALL shdf5_open (flnm,'R',iphdf5,h5_fid)
   ! Leaf 3 2D vars
   CALL shdf5_set_hs_select (6,'R',ng,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'VEG_NDVIF',mem_select &
                   ,file_select,rvara=leaf_g(ng)%veg_ndvif)
   CALL shdf5_close (h5_fid)
enddo

!Saleeby(2010): Set grid boundarys as in leaf3 for NDVI
!Do this here so that boundary assignments are made before the first
!analysis write and first timestep.
do ng=1,ngrids
 do ip = 1,npatch
   do j = 1,mmyp(ng)
     leaf_g(ng)%veg_ndvif(1,j,ip)=leaf_g(ng)%veg_ndvif(2,j,ip)
     leaf_g(ng)%veg_ndvif(mmxp(ng),j,ip)=leaf_g(ng)%veg_ndvif(mmxp(ng)-1,j,ip)
   enddo
   if (jdim == 1) then
    do i = 1,mmxp(ng)
      leaf_g(ng)%veg_ndvif(i,1,ip)=leaf_g(ng)%veg_ndvif(i,2,ip)
      leaf_g(ng)%veg_ndvif(i,mmyp(ng),ip)=leaf_g(ng)%veg_ndvif(i,mmyp(ng)-1,ip)
    enddo
   endif
 enddo
enddo

return
END SUBROUTINE ndvi_update

!##############################################################################
Subroutine ndvi_check_header (ifm,flnm,ierr)

! This routine checks for the existence of an ndvi file for
! grid number ifm, and if it exists, also checks for agreement of
! grid configuration between the file and the current model run.
! If the file does not exist or does not match grid configuration,
! the flag ierr is returned with a value of 1.  If the file
! exists and is ok, ierr is returned with a value of 0.

use mem_grid
use node_mod
use hdf5_utils

implicit none

integer :: ifm,ierr
character(len=*) :: flnm

integer :: nsfx,nsfy,nsfpat
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
 print '(a,i3,1x,a,a,$)',' ---> Check grid:',ifm,' ndvi filename:',trim(flnm)

CALL xy_ll (glatr,glonr,polelat,polelon,xtn(1,ifm),ytn(1,ifm))


CALL shdf5_open (flnm,'R',iphdf5,h5_fid)

! Scalar vars (set 3rd argument, igrid, to 1 since always have at least 1 grid)
CALL shdf5_set_hs_select (1,'R',1,mem_select,file_select,file_chunks)
CALL shdf5_irec (h5_fid,iphdf5,'nx',mem_select,file_select,ivars=nsfx)
CALL shdf5_irec (h5_fid,iphdf5,'ny',mem_select,file_select,ivars=nsfy)
CALL shdf5_irec (h5_fid,iphdf5,'npatch',mem_select,file_select,ivars=nsfpat)
CALL shdf5_irec (h5_fid,iphdf5,'dx',mem_select,file_select,rvars=sfdx)
CALL shdf5_irec (h5_fid,iphdf5,'polelat',mem_select,file_select,rvars=sfplat)
CALL shdf5_irec (h5_fid,iphdf5,'polelon',mem_select,file_select,rvars=sfplon)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lat',mem_select,file_select,rvars=sflat)
CALL shdf5_irec (h5_fid,iphdf5,'sw_lon',mem_select,file_select,rvars=sflon)

CALL shdf5_close (h5_fid)

if (nsfx                   .ne. nnxp(ifm) .or.  &
    nsfy                   .ne. nnyp(ifm) .or.  &
    nsfpat                 .ne. npatch    .or.  &
    abs(sfdx-deltaxn(ifm)) .gt. .001      .or.  &
    abs(sfplat-polelat)    .gt. .001      .or.  &
    abs(sfplon-polelon)    .gt. .001      .or.  &
    abs(sflat-glatr)       .gt. .001      .or.  &
    abs(sflon-glonr)       .gt. .001) then

   ierr = 1

   print*
   print*,'ndvifile mismatch on grid:',ifm
   print*,'Values: model, file'
   print*,'-------------------'
   print*,'nnxp:',nnxp(ifm),nsfx
   print*,'nnyp:',nnyp(ifm),nsfy
   print*,'npatch:',npatch,nsfpat
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
END SUBROUTINE ndvi_check_header

