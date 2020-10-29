!##############################################################################
Module mem_oda

use grid_dims

implicit none

   Type oda_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                            uk,vk,tk,rk,ukv,vkv,tkv,rkv
                           
   End Type
   
   type (oda_vars), allocatable :: oda_g(:), odam_g(:)

character(len=strl1), dimension(maxodafiles) :: fnames_upa, fnames_sfc
character(len=14), dimension(maxodafiles) :: itotdate_upa,itotdate_sfc
character(len=16), dimension(maxodasta) :: staid_sfc, staid_upa
integer, dimension(maxodasta) :: ntimes_sfc, ntimes_upa
integer :: maxtimes_sfc, maxtimes_upa

! Namelist inputs

character(len=strl1) :: oda_upapref,oda_sfcpref
integer :: if_oda
real :: frqoda,todabeg,todaend,tnudoda,wt_oda_grid(maxodagrids)  &
       ,oda_sfc_til,oda_sfc_tel,oda_upa_til,oda_upa_tel  &
       ,wt_oda_uv,wt_oda_th,wt_oda_pi,wt_oda_rt
real, dimension(maxodagrids) :: roda_sfce, roda_sfc0  &
                              , roda_upae, roda_upa0  &
                              , roda_zfac, roda_hgt


integer ::  nupafiles,nsfcfiles,num_oda_sfc,num_oda_upa

! Surface data 

type oda_sfc_info_type
   character(len=16) :: id
   integer :: intid
   integer :: ntimes
   integer, dimension(maxodagrids) :: iactive
   real, dimension(maxodagrids) :: xista, xjsta
   real :: xlat,xlon,xsta,ysta,stopo
end type   
   
   
type oda_sfc_type
   real, allocatable, dimension(:) :: temp, dewpt, us, vs, ps,u,v 
   real, allocatable, dimension(:) :: time 
   real                        :: psref
end type   
   
type(oda_sfc_info_type), allocatable :: oda_sfc_info(:)
type(oda_sfc_type)     , allocatable :: oda_sfc_obs(:)
   

! Upper air info

type oda_upa_info_type
   character(len=16) :: id
   integer :: intid
   integer :: ntimes
   integer, dimension(maxodagrids) :: iactive
   real, dimension(maxodagrids) :: xista, xjsta
   real :: xlat,xlon,xsta,ysta,stopo
end type   
   
   ! Upper air data
type oda_upa_type
   character(len=14) :: ctotdate
   real, allocatable, dimension(:,:) :: theta, rv, us, vs, u, v, zz, pi, zgeo
   real, allocatable, dimension(:) :: time 
   integer, allocatable, dimension(:) :: lp,lz
end type   
   
type(oda_upa_info_type), allocatable :: oda_upa_info(:)
type(oda_upa_type)     , allocatable :: oda_upa_obs(:)

! Krigging routine info

real :: rmaxkrg(maxodanzp,maxodagrids)  &
       ,ckrg(3,maxodagrids),akrg(maxodanzp,maxodagrids)  &
       ,caxkrg(9,maxodagrids),caykrg(9,maxodagrids),cazkrg(9,maxodagrids)
integer :: nstkrg(maxodagrids)


! Filled obs arrays for an analysis time

integer, dimension(maxkobs) :: idobs
real, dimension(maxkobs) :: ukobs,vkobs,tkobs,rkobs,pkobs  &
                           ,xkobs,ykobs,zkobs,ekobs  &
                           ,ikobs,jkobs

contains

!##############################################################################
Subroutine alloc_oda (oda,n1,n2,n3)

implicit none

   type (oda_vars) :: oda
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)

   if( if_oda == 1 ) then
                     allocate (oda%uk(n1,n2,n3))
                     allocate (oda%vk(n1,n2,n3))
                     allocate (oda%tk(n1,n2,n3))
                     allocate (oda%rk(n1,n2,n3))
                     allocate (oda%ukv(n1,n2,n3))
                     allocate (oda%vkv(n1,n2,n3))
                     allocate (oda%tkv(n1,n2,n3))
                     allocate (oda%rkv(n1,n2,n3))
   endif
  
return
END SUBROUTINE alloc_oda

!##############################################################################
Subroutine dealloc_oda (oda)

implicit none

   type (oda_vars) :: oda

   if (allocated(oda%uk))      deallocate (oda%uk)
   if (allocated(oda%vk))      deallocate (oda%vk)
   if (allocated(oda%tk))     deallocate (oda%tk)
   if (allocated(oda%rk))     deallocate (oda%rk)
   if (allocated(oda%ukv))     deallocate (oda%ukv)
   if (allocated(oda%vkv))     deallocate (oda%vkv)
   if (allocated(oda%tkv))    deallocate (oda%tkv)
   if (allocated(oda%rkv))    deallocate (oda%rkv)

return
END SUBROUTINE dealloc_oda

!##############################################################################
Subroutine filltab_oda (oda,odam,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (oda_vars) :: oda,odam
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(oda%uk))  &
      CALL vtables2 (oda%uk(1,1,1),odam%uk(1,1,1)  &
                 ,ng, npts, imean,  &
                 'UKODA :3:')
   if (allocated(oda%vk))  &
      CALL vtables2 (oda%vk(1,1,1),odam%vk(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VKODA :3:')
   if (allocated(oda%tk))  &
      CALL vtables2 (oda%tk(1,1,1),odam%tk(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TKODA :3:')
   if (allocated(oda%rk))  &
      CALL vtables2 (oda%rk(1,1,1),odam%rk(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RKODA :3:')
   if (allocated(oda%ukv))  &
      CALL vtables2 (oda%ukv(1,1,1),odam%ukv(1,1,1)  &
                 ,ng, npts, imean,  &
                 'UVODA :3:')
   if (allocated(oda%vkv))  &
      CALL vtables2 (oda%vkv(1,1,1),odam%vkv(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VVODA :3:')
   if (allocated(oda%tkv))  &
      CALL vtables2 (oda%tkv(1,1,1),odam%tkv(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TVODA :3:')
   if (allocated(oda%rkv))  &
      CALL vtables2 (oda%rkv(1,1,1),odam%rkv(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RVODA :3:')
                 
return
END SUBROUTINE filltab_oda

END MODULE mem_oda   
