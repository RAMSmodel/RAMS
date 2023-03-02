!##############################################################################
Module mem_nestb

use grid_dims

implicit none

type nest_bounds

   real, allocatable, dimension(:,:,:)   :: bux,buy,buz,bvx,bvy,bvz  &
                                       ,bwx,bwy,bwz,bpx,bpy,bpz
   real, allocatable, dimension(:,:,:,:) :: bsx,bsy,bsz
   
end type


type (nest_bounds) :: nbounds(maxgrds)


Contains

!##############################################################################
Subroutine alloc_nestb (ng,nx,ny,nz)

use var_tables

implicit none

integer :: ng,nx,ny,nz

!  Allocate "b" array components. All grids will be allocated,
!     only to 1's if nesting isn't done.

allocate( nbounds(ng)%bux(nz,ny,2) )
allocate( nbounds(ng)%buy(nz,nx,2) )
allocate( nbounds(ng)%buz(nx,ny,2) )

allocate( nbounds(ng)%bvx(nz,ny,2) )
allocate( nbounds(ng)%bvy(nz,nx,2) )
allocate( nbounds(ng)%bvz(nx,ny,2) )

allocate( nbounds(ng)%bwx(nz,ny,2) )
allocate( nbounds(ng)%bwy(nz,nx,2) )
allocate( nbounds(ng)%bwz(nx,ny,2) )

allocate( nbounds(ng)%bpx(nz,ny,2) )
allocate( nbounds(ng)%bpy(nz,nx,2) )
allocate( nbounds(ng)%bpz(nx,ny,2) )

allocate( nbounds(ng)%bsx(nz,ny,2,num_scalar(ng)) )
allocate( nbounds(ng)%bsy(nz,nx,2,num_scalar(ng)) )
allocate( nbounds(ng)%bsz(nx,ny,2,num_scalar(ng)) )

return
END SUBROUTINE alloc_nestb

END MODULE mem_nestb
