!##############################################################################
Module mem_leaf

use grid_dims

implicit none

   Type leaf_vars
   
      ! Variables to be dimensioned by (nzg,nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:,:) :: &
                  soil_water, soil_energy, soil_text

      ! Variables to be dimensioned by (nzs,nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:,:) :: &
                  sfcwater_mass, sfcwater_energy, sfcwater_depth

      ! Variables to be dimensioned by (nxp,nyp,npatch)
   real, allocatable, dimension(:,:,:) :: &
                  ustar,tstar,rstar  &  
                 ,veg_fracarea,veg_lai,veg_rough,veg_height  &
                 ,veg_albedo,veg_tai  &
                 ,patch_area,patch_rought,patch_roughm,leaf_class  &
                 ,soil_rough,sfcwater_nlev,stom_resist  &
                 ,ground_rsat,ground_rvap  &
                 ,veg_water,veg_temp,can_rvap,can_temp &
                 ,veg_ndvip,veg_ndvic,veg_ndvif

      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                  snow_mass,snow_depth,seatp,seatf &
                 ,soil_moist_top,soil_moist_bot,soil_temp_top,soil_temp_bot
   End Type
   
   type (leaf_vars), allocatable :: leaf_g(:), leafm_g(:)

   !-------------------------------------------------------------------------------
   integer                   :: nslcon,nvgcon,nvegpat,isfcl,isoildat,isnowdat
   real                      :: ztrough,zmrough=0.005,pctlcon,ubmin,albedo,drtcon,dthcon,seatmp
   real, dimension(nzgmax)   :: stgoff,slmstr
   real, dimension(nzgmax+1) :: slz
   character(len=strl1)      :: sibfile

Contains

!##############################################################################
Subroutine alloc_leaf (leaf,nx,ny,nzg,nzs,np)

implicit none

   type (leaf_vars) :: leaf
   integer, intent(in) :: nx,ny,nzg,nzs,np


! Allocate arrays based on options (if necessary)

   allocate (leaf%soil_water     (nzg,nx,ny,np))
   allocate (leaf%soil_energy    (nzg,nx,ny,np))
   allocate (leaf%soil_text      (nzg,nx,ny,np))

   allocate (leaf%sfcwater_mass  (nzs,nx,ny,np))
   allocate (leaf%sfcwater_energy(nzs,nx,ny,np))
   allocate (leaf%sfcwater_depth (nzs,nx,ny,np))

   allocate (leaf%ustar        (nx,ny,np))
   allocate (leaf%tstar        (nx,ny,np))
   allocate (leaf%rstar        (nx,ny,np))
   
   allocate (leaf%veg_fracarea (nx,ny,np))
   allocate (leaf%veg_lai      (nx,ny,np))
   allocate (leaf%veg_rough    (nx,ny,np))
   allocate (leaf%veg_height   (nx,ny,np))
   allocate (leaf%veg_albedo   (nx,ny,np))
   allocate (leaf%veg_tai      (nx,ny,np))

   allocate (leaf%patch_area   (nx,ny,np))
   allocate (leaf%patch_rought (nx,ny,np))
   allocate (leaf%patch_roughm (nx,ny,np))
   allocate (leaf%leaf_class   (nx,ny,np))
   
   allocate (leaf%soil_rough   (nx,ny,np))
   allocate (leaf%sfcwater_nlev(nx,ny,np))
   allocate (leaf%stom_resist  (nx,ny,np))

   allocate (leaf%ground_rsat  (nx,ny,np))
   allocate (leaf%ground_rvap  (nx,ny,np))
   
   allocate (leaf%veg_water    (nx,ny,np))
   allocate (leaf%veg_temp     (nx,ny,np))

   allocate (leaf%can_rvap     (nx,ny,np))
   allocate (leaf%can_temp     (nx,ny,np))

   allocate (leaf%veg_ndvip    (nx,ny,np))
   allocate (leaf%veg_ndvic    (nx,ny,np))
   allocate (leaf%veg_ndvif    (nx,ny,np))
   
   allocate (leaf%snow_mass    (nx,ny))
   allocate (leaf%snow_depth   (nx,ny))
   allocate (leaf%seatp        (nx,ny))
   allocate (leaf%seatf        (nx,ny))
   allocate (leaf%soil_moist_top(nx,ny))
   allocate (leaf%soil_moist_bot(nx,ny))
   allocate (leaf%soil_temp_top(nx,ny))
   allocate (leaf%soil_temp_bot(nx,ny))

return
END SUBROUTINE alloc_leaf

!##############################################################################
Subroutine dealloc_leaf (leaf)

implicit none

  type (leaf_vars) :: leaf

  if(allocated(leaf%soil_water))      deallocate (leaf%soil_water)
  if(allocated(leaf%soil_energy))     deallocate (leaf%soil_energy)
  if(allocated(leaf%soil_text))       deallocate (leaf%soil_text)

  if(allocated(leaf%sfcwater_mass))   deallocate (leaf%sfcwater_mass)
  if(allocated(leaf%sfcwater_energy)) deallocate (leaf%sfcwater_energy)
  if(allocated(leaf%sfcwater_depth))  deallocate (leaf%sfcwater_depth)

  if(allocated(leaf%ustar))           deallocate (leaf%ustar)
  if(allocated(leaf%tstar))           deallocate (leaf%tstar)
  if(allocated(leaf%rstar))           deallocate (leaf%rstar)
 
  if(allocated(leaf%veg_fracarea))    deallocate (leaf%veg_fracarea)
  if(allocated(leaf%veg_lai))         deallocate (leaf%veg_lai)
  if(allocated(leaf%veg_rough))       deallocate (leaf%veg_rough)
  if(allocated(leaf%veg_height))      deallocate (leaf%veg_height)
  if(allocated(leaf%veg_albedo))      deallocate (leaf%veg_albedo)
  if(allocated(leaf%veg_tai))         deallocate (leaf%veg_tai)

  if(allocated(leaf%patch_area))      deallocate (leaf%patch_area)
  if(allocated(leaf%patch_rought))    deallocate (leaf%patch_rought)
  if(allocated(leaf%patch_roughm))    deallocate (leaf%patch_roughm)
  if(allocated(leaf%leaf_class))      deallocate (leaf%leaf_class)
 
  if(allocated(leaf%soil_rough))      deallocate (leaf%soil_rough)
  if(allocated(leaf%sfcwater_nlev))   deallocate (leaf%sfcwater_nlev)
  if(allocated(leaf%stom_resist))     deallocate (leaf%stom_resist)

  if(allocated(leaf%ground_rsat))     deallocate (leaf%ground_rsat)
  if(allocated(leaf%ground_rvap))     deallocate (leaf%ground_rvap)
 
  if(allocated(leaf%veg_water))       deallocate (leaf%veg_water)
  if(allocated(leaf%veg_temp))        deallocate (leaf%veg_temp)

  if(allocated(leaf%can_rvap))        deallocate (leaf%can_rvap)
  if(allocated(leaf%can_temp))        deallocate (leaf%can_temp)

  if(allocated(leaf%veg_ndvip))       deallocate (leaf%veg_ndvip)
  if(allocated(leaf%veg_ndvic))       deallocate (leaf%veg_ndvic)
  if(allocated(leaf%veg_ndvif))       deallocate (leaf%veg_ndvif)
   
  if(allocated(leaf%snow_mass))       deallocate (leaf%snow_mass)
  if(allocated(leaf%snow_depth))      deallocate (leaf%snow_depth)
  if(allocated(leaf%seatp))           deallocate (leaf%seatp)
  if(allocated(leaf%seatf))           deallocate (leaf%seatf)
  if(allocated(leaf%soil_moist_top))  deallocate (leaf%soil_moist_top)
  if(allocated(leaf%soil_moist_bot))  deallocate (leaf%soil_moist_bot)
  if(allocated(leaf%soil_temp_top))   deallocate (leaf%soil_temp_top)
  if(allocated(leaf%soil_temp_bot))   deallocate (leaf%soil_temp_bot)

return
END SUBROUTINE dealloc_leaf

!##############################################################################
Subroutine filltab_leaf (leaf,leafm,imean,nx,ny,nzg,nzs,np,ng)

use var_tables

implicit none

   type (leaf_vars) :: leaf,leafm
   integer, intent(in) :: imean,nx,ny,nzg,nzs,np,ng
   integer :: npts

! Fill arrays into variable tables

   npts=nzg*nx*ny*np
   CALL vtables2 (leaf%soil_water(1,1,1,1),leafm%soil_water(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_WATER :4:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%soil_energy(1,1,1,1),leafm%soil_energy(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_ENERGY :4:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%soil_text(1,1,1,1),leafm%soil_text(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_TEXT :4:anal:mpti')

   npts=nzs*nx*ny*np
   CALL vtables2 (leaf%sfcwater_mass(1,1,1,1),leafm%sfcwater_mass(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'SFCWATER_MASS :5:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%sfcwater_energy(1,1,1,1), leafm%sfcwater_energy(1,1,1,1) &
                 ,ng, npts, imean,  &
                 'SFCWATER_ENERGY :5:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%sfcwater_depth(1,1,1,1),leafm%sfcwater_depth(1,1,1,1)  &
                 ,ng, npts, imean,  &
                 'SFCWATER_DEPTH :5:anal:mpti:recycle_sfc')

   npts=nx*ny*np
   CALL vtables2 (leaf%ustar(1,1,1),leafm%ustar(1,1,1)  &
                 ,ng, npts, imean,  &
                 'USTAR :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%tstar(1,1,1),leafm%tstar(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TSTAR :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%rstar(1,1,1),leafm%rstar(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RSTAR :6:anal:mpti:recycle_sfc')
            
   CALL vtables2 (leaf%veg_fracarea(1,1,1),leafm%veg_fracarea(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_FRACAREA :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_lai(1,1,1),leafm%veg_lai(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_LAI :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_rough(1,1,1),leafm%veg_rough(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_ROUGH :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_height(1,1,1),leafm%veg_height(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_HEIGHT :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_albedo(1,1,1),leafm%veg_albedo(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_ALBEDO :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_tai(1,1,1),leafm%veg_tai(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_TAI :6:anal:mpti:recycle_sfc')

   CALL vtables2 (leaf%patch_area(1,1,1),leafm%patch_area(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PATCH_AREA :6:anal:mpti')
   CALL vtables2 (leaf%patch_rought(1,1,1),leafm%patch_rought(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PATCH_ROUGHT :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%patch_roughm(1,1,1),leafm%patch_roughm(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PATCH_ROUGHM :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%leaf_class(1,1,1),leafm%leaf_class(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LEAF_CLASS :6:anal:mpti')

   CALL vtables2 (leaf%soil_rough(1,1,1),leafm%soil_rough(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_ROUGH :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%sfcwater_nlev(1,1,1),leafm%sfcwater_nlev(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SFCWATER_NLEV :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%stom_resist(1,1,1),leafm%stom_resist(1,1,1)  &
                 ,ng, npts, imean,  &
                 'STOM_RESIST :6:anal:mpti:recycle_sfc')

   CALL vtables2 (leaf%ground_rsat(1,1,1),leafm%ground_rsat(1,1,1)  &
                 ,ng, npts, imean,  &
                 'GROUND_RSAT :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%ground_rvap(1,1,1),leafm%ground_rvap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'GROUND_RVAP :6:anal:mpti:recycle_sfc')

   CALL vtables2 (leaf%veg_water(1,1,1),leafm%veg_water(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_WATER :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_temp(1,1,1),leafm%veg_temp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_TEMP :6:anal:mpti:recycle_sfc')

   CALL vtables2 (leaf%can_rvap(1,1,1),leafm%can_rvap(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CAN_RVAP :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%can_temp(1,1,1),leafm%can_temp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'CAN_TEMP :6:anal:mpti:recycle_sfc')

   CALL vtables2 (leaf%veg_ndvip(1,1,1),leafm%veg_ndvip(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_NDVIP :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_ndvic(1,1,1),leafm%veg_ndvic(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_NDVIC :6:anal:mpti:recycle_sfc')
   CALL vtables2 (leaf%veg_ndvif(1,1,1),leafm%veg_ndvif(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VEG_NDVIF :6:anal:mpti')

   npts=nx*ny
   CALL vtables2 (leaf%snow_mass(1,1),leafm%snow_mass(1,1)  &
                 ,ng, npts, imean,  &
                 'SNOW_MASS :2:mpti')
   CALL vtables2 (leaf%snow_depth(1,1),leafm%snow_depth(1,1)  &
                 ,ng, npts, imean,  &
                 'SNOW_DEPTH :2:mpti')
   CALL vtables2 (leaf%seatp(1,1),leafm%seatp(1,1)  &
                 ,ng, npts, imean,  &
                 'SEATP :2:mpti')
   CALL vtables2 (leaf%seatf(1,1),leafm%seatf(1,1)  &
                 ,ng, npts, imean,  &
                 'SEATF :2:mpti')
   CALL vtables2 (leaf%soil_moist_top(1,1),leafm%soil_moist_top(1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_MOIST_TOP :2:mpti')
   CALL vtables2 (leaf%soil_moist_bot(1,1),leafm%soil_moist_bot(1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_MOIST_BOT :2:mpti')
   CALL vtables2 (leaf%soil_temp_top(1,1),leafm%soil_temp_top(1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_TEMP_TOP :2:mpti')
   CALL vtables2 (leaf%soil_temp_bot(1,1),leafm%soil_temp_bot(1,1)  &
                 ,ng, npts, imean,  &
                 'SOIL_TEMP_BOT :2:mpti')

return
END SUBROUTINE filltab_leaf

END MODULE mem_leaf
