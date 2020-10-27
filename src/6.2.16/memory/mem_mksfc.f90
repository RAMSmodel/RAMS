!##############################################################################
Module mem_mksfc

use grid_dims

implicit none

   Type sfcfile_vars

!(nxp,nyp,nzg,npatch)
   real, allocatable, dimension(:,:,:,:) :: soil_text

!(nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:) :: patch_area,leaf_class,veg_ndvif

!(nxp,nyp)
   real, allocatable, dimension(:,:) :: topt,seatf,topzo

   End Type
   


   type (sfcfile_vars), allocatable :: sfcfile_p(:)
  
!(nxpmax,nypmax)
   real, allocatable, dimension(:,:) :: scr1,scr2,vt2da,vt2db

!(np,np,nxp,nyp)
   real, allocatable, dimension(:,:,:,:) :: glatp,glonp,datp

!(np,np,nxp,nyp)
   integer, allocatable, dimension(:,:,:,:) :: datq_patch

!(np*np*nxp*nyp)
   integer, allocatable, dimension(:) :: ptable

!(iblksizo,iblksizo)
   real, allocatable, dimension(:,:) :: dato

!(iblksizo,iblksizo)
   character(len=1), allocatable, dimension(:,:) :: cdato
   integer, allocatable, dimension(:,:) :: idato

!(ifile_max,jfile_max)
   integer, allocatable, dimension(:,:) :: nump,numpind,numpind1,numpind2

   integer :: npq
   
   ! SST file creation variables
   integer, dimension(maxsstdata,maxgrds) :: iyearvs,imonthvs,idatevs,ihourvs
   integer,dimension(maxgrds) :: nvsstf
   character(len=strl1), dimension(maxsstdata,maxgrds) :: vsstfil
   
   ! NDVI file creation variables
   integer, dimension(maxndvidata,maxgrds) :: iyearvn,imonthvn,idatevn,ihourvn
   integer,dimension(maxgrds) :: nvndvif
   character(len=strl1), dimension(maxndvidata,maxgrds) :: vndvifil

Contains

!##############################################################################
Subroutine alloc_sfcfile (sfcfile,nx,ny,nzg,npat)

implicit none

   type (sfcfile_vars) :: sfcfile
   integer, intent(in) :: nx,ny,nzg,npat

   !print*,'alloc_sfc:',nx,ny,nzg,npat
   allocate (sfcfile%soil_text    (nzg,nx,ny,npat))

   allocate (sfcfile%patch_area   (nx,ny,npat))
   allocate (sfcfile%leaf_class  (nx,ny,npat))
   allocate (sfcfile%veg_ndvif    (nx,ny,npat))

   allocate (sfcfile%topt         (nx,ny))
   allocate (sfcfile%seatf        (nx,ny))
   allocate (sfcfile%topzo        (nx,ny))
      
return
END SUBROUTINE alloc_sfcfile

!##############################################################################
Subroutine dealloc_sfcfile (sfcfile)

implicit none

   type (sfcfile_vars) :: sfcfile

   deallocate (sfcfile%soil_text)

   deallocate (sfcfile%patch_area)
   deallocate (sfcfile%leaf_class)
   deallocate (sfcfile%veg_ndvif)
   
   deallocate (sfcfile%topt)
   deallocate (sfcfile%seatf)
   deallocate (sfcfile%topzo)
   
return
END SUBROUTINE dealloc_sfcfile

END MODULE mem_mksfc
