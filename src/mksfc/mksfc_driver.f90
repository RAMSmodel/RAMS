!##############################################################################
Subroutine make_sfcfiles ()

!TAKE INTO ACCOUNT ALL OPTIONS, HISTORY STARTS, ADDED GRIDS, ETC.

use mem_mksfc
use mem_grid
use io_params
use node_mod
use mem_leaf, only:isfcl

implicit none

integer :: ifm,nvtime,ivtime,ng1,ng2,ng1t,ng2t,ng1s,ng2s
integer :: isfcerr,itoperr,issterr,indvierr

! This routine makes sure that all surface, topo, sst, and ndvi
! files required for the present run exist and are correct. 
 
! The logical choices made are as follows:

! For runtype = 'MAKESFC':
!    Make all surface, topo, sst, and ndvi files for grids 1:NGRIDS.
! For runtype = 'INITIAL', 'ERROR', or 'MAKEVFILE', 'MAKEHFILE, 'HISTORY':
!    Check for existence and correctness of all surface, topo, sst,
!       and ndvi files for grids 1:NGRIDS.  Remake all grids if the files
!       are incorrect for the set of files: topo, sfc/ndvi, sst

isfcerr = 0
itoperr = 0
issterr = 0
indvierr = 0

! Allocate memory needed for initializing sfcfiles

!print*,'before alloc_sfcfile',ngrids
allocate( sfcfile_p(ngrids) )
do ifm = 1,ngrids
   CALL alloc_sfcfile (sfcfile_p(ifm),mmxp(ifm),mmyp(ifm),nzg,npatch)
enddo

allocate (scr1 (nxpmax,nypmax))
allocate (scr2 (nxpmax,nypmax))
allocate (vt2da(nxpmax,nypmax))
allocate (vt2db(nxpmax,nypmax))
   
!print*,'after allocations',runtype,nxpmax,nypmax,maxval(nnxyp(1:ifm))*nzg*npatch

ng1=1  ; ng2 =ngrids !Variable initialized
ng1t=1 ; ng2t=ngrids !Variable initialized
ng1s=1 ; ng2s=ngrids !Variable initialized

if (trim(runtype) == 'MAKESFC') then

   if (isfcl .ne. 3) then
     print*, 'MAKESFC run: Making all surface, topo, sst, and ndvi files.'
     itoperr = 1
     isfcerr = 1
     issterr = 1
     indvierr = 1
   else
     print*,'MAKESFC run: Making topo file.'
     itoperr = 1
   endif
   
   ng1=1  ; ng2 =ngrids    ! sst,ndvi grid bounds
   ng1t=1 ; ng2t=ngrids    ! topo grid bounds
   ng1s=1 ; ng2s=ngrids    ! sfc grid bounds

elseif (trim(runtype) == 'INITIAL' .or. trim(runtype) == 'MAKEVFILE'  &
                                   .or. trim(runtype) == 'MAKEHFILE'  &
                                   .or. trim(runtype) == 'HISTORY'    &
                                   .or. trim(runtype) == 'ERROR') then
   
   ! Check topo files 
   do ifm = 1,ngrids
      CALL top_check (ifm,itoperr)
      if(itoperr == 1) exit
   enddo
   
   if (isfcl .ne. 3) then
     ! Check sfc files 
     do ifm = 1,ngrids
        CALL sfc_check (ifm,isfcerr)
        if(isfcerr == 1) exit
     enddo
   
     ! Check sst files
     CALL sst_read (2,ifm,issterr)

     ! Check ndvi files
     CALL ndvi_read (2,ifm,indvierr)
   
     ! If we are making ndvi files, we must also make the sfc files (and vice versa)
     if(indvierr == 1) isfcerr = 1
     if(isfcerr == 1) indvierr = 1
   
   endif
   if (isfcerr==0 .and. issterr==0 .and.  &
       itoperr==0 .and.indvierr==0) then
      print*, 'Surface, topo, sst, and ndvi files all ok for'
      print*, '   RUNTYPE = ',trim(runtype)

      ! Deallocate memory needed for initializing sfcfiles
      do ifm = 1,ngrids
         CALL dealloc_sfcfile (sfcfile_p(ifm))
      enddo
      deallocate (sfcfile_p)

      deallocate (scr1,scr2,vt2da,vt2db)

      return
   else
      print*, 'Nonexistent or incorrect surface files for'
      print*, '   RUNTYPE = ',trim(runtype),'...(re)making:'
      if(isfcerr == 1 ) print*, '   sfc files'
      if(itoperr == 1 ) print*, '   top files'
      if(issterr == 1 ) print*, '   sst files'
      if(indvierr == 1) print*, '   ndvi files'
   endif
   
   ng1=1 ; ng2=ngrids
   ng1t=1 ; ng2t=ngrids
   ng1s=1 ; ng2s=ngrids

endif

! If we got here, at least one set of files are bad. Re-make the bad ones.

!------------------------------------------
!  TOP (topo and roughness) file creation
if(itoperr == 1) then
   ! do topography, topo roughness on all grids
   CALL toptnest (ng1t,ng2t)
   do ifm = 1,ngrids
      CALL top_write (ifm)
   enddo
endif

!------------------------------------------
!  SFC (veg class, patch area, soil type) and NDVI file creation
!      If we are making ndvi files, we must also make the sfc files (and vice versa)
if(isfcerr == 1 .or. indvierr == 1) then
   
   ! If iupdndvi = 1, require that:
   !    (1) ndviflg = 1 for a grid that has nxtnest = 0,
   !    (2) ndviflg /= 2 for all other grids.

   if (iupdndvi == 1) then
      do ifm = 1,ngrids
         if (ndviflg(ifm) /= 1 .and. nxtnest(ifm) == 0) then
            print*, 'iupdndvi = 1 and ndviflg /= 1 for grid ', ifm
            stop 'iupdndvi'
         endif

         if (ndviflg(ifm) == 2) then
            print*, 'iupdndvi = 1 and ndviflg = 2 for grid ', ifm
            stop 'iupdndvi'
         endif
      enddo
   endif
   
   
   do ifm = ng1s,ng2s

      if (ivegtflg(ifm) == 1 .or. isoilflg(ifm) == 1 .or. ndviflg(ifm) == 1) then   

         ! Find size of patch arrays
         CALL patch_array_size (npq,(xtn(2,ifm)-xtn(1,ifm))  &
               ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm)  &
               ,ndviflg(ifm),ndvifn(ifm) )

         ! Allocate arrays that need npq dimension
         allocate (glatp(npq,npq,mmxp(ifm),mmyp(ifm))  &
                  ,glonp(npq,npq,mmxp(ifm),mmyp(ifm))  &
                  ,datq_patch(npq,npq,mmxp(ifm),mmyp(ifm))  &
                  ,datp(npq,npq,mmxp(ifm),mmyp(ifm)) )

         ! Fill lat-lon
         ! Need to line up XTN and YTN with this sub-domain. Do this by setting the
         ! addresses for the xt and yt arguments to 1 plus the x and y offsets associated
         ! with the sub-domain.
         CALL patch_latlon (mmxp(ifm),mmyp(ifm)  &
               ,xtn(mi0(ifm)+1,ifm),ytn(mj0(ifm)+1,ifm),deltaxn(ifm)  &
               ,polelat,polelon ) 
      endif
   

      ! do sfcfile
      CALL geonest_file (ifm)
      CALL sfc_write (ifm)

      
      ! do ndvifile
      if (ndviflg(ifm) == 1) then
         CALL ndvi_read_dataheader (ifm)
         nvtime = nvndvif(ifm)
      elseif (ndviflg(ifm) == 0) then
         nvtime = nvndvif(nxtnest(ifm))
      else
         nvtime = 1
      endif
   
      do ivtime = 1,nvtime
         CALL ndvinest (ifm,ivtime)
         CALL ndvi_write (ifm,ivtime)
      enddo
      
      ! Deallocate arrays that needed npq dimension
      if(allocated(glatp)) deallocate (glatp,glonp,datq_patch,datp)

   enddo
   
endif
   
!------------------------------------------
!  SST file creation
if(issterr == 1) then
   
   ! If iupdsst = 1, require that:
   !    (1) isstflg = 1 for a grid that has nxtnest = 0,
   !    (2) isstflg /= 2 for all other grids.

   if (iupdsst == 1) then
      do ifm = 1,ngrids
         if (isstflg(ifm) /= 1 .and. nxtnest(ifm) == 0) then
            print*, 'iupdsst = 1 and isstflg /= 1 for grid ', ifm
            stop 'iupdsst'
         endif

         if (isstflg(ifm) == 2) then
            print*, 'iupdsst = 1 and isstflg = 2 for grid ', ifm
            stop 'iupdsst'
         endif
      enddo
   endif
   
   do ifm = ng1,ng2
      ! do sstfile
      if (isstflg(ifm) == 1) then
         CALL sst_read_dataheader (ifm)
         nvtime = nvsstf(ifm)
      elseif (isstflg(ifm) == 0) then
         nvtime = nvsstf(nxtnest(ifm))
      else
         nvtime = 1
      endif
   
      do ivtime = 1,nvtime
         CALL sstnest (ifm,ivtime)
         CALL sst_write (ifm,ivtime)
      enddo
   enddo

endif


! Deallocate memory needed for initializing sfcfiles
do ifm = 1,ngrids
   CALL dealloc_sfcfile (sfcfile_p(ifm))
enddo
deallocate (sfcfile_p)

deallocate (scr1,scr2,vt2da,vt2db)

return
END SUBROUTINE make_sfcfiles

