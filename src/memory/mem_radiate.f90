!##############################################################################
Module mem_radiate

implicit none

   Type radiate_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          fthrd,fthrdp,bext,swup,swdn,lwup,lwdn,fthrdsw,fthrdlw
                          !GRL 2024-03-22 added variables for lw and sw heating rates

                          
      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          rshort,rlong,rlongup,rlontop,albedt,cosz,aodt 

   End Type
   
   type (radiate_vars), allocatable :: radiate_g(:), radiatem_g(:)
   
   integer :: lonrad,ilwrtyp,iswrtyp,irce
   real    :: radfrq,rce_ubmn,rce_bubl,rce_solc,rce_szen
  
Contains

!##############################################################################
Subroutine alloc_radiate (radiate,n1,n2,n3)

implicit none

   type (radiate_vars) :: radiate
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)
      
      if(ilwrtyp+iswrtyp > 0)  then
                         allocate (radiate%fthrd(n1,n2,n3))
                         allocate (radiate%fthrdp(n1,n2,n3))
                         allocate (radiate%rshort(n2,n3))
                         allocate (radiate%rlong(n2,n3))
                         allocate (radiate%rlongup(n2,n3))
                         allocate (radiate%rlontop(n2,n3))
                         allocate (radiate%albedt(n2,n3))
                         allocate (radiate%cosz(n2,n3))
                         allocate (radiate%aodt(n2,n3))
                         !GRL 2024-03-22 added variables for lw and sw heating rates
                         allocate (radiate%fthrdsw(n1,n2,n3))
                         allocate (radiate%fthrdlw(n1,n2,n3))
      endif
      if(ilwrtyp >= 3 .or. iswrtyp >= 3) then
         allocate (radiate%bext(n1,n2,n3))
      endif
      if(ilwrtyp >= 3)  then
         allocate (radiate%lwup(n1,n2,n3))
         allocate (radiate%lwdn(n1,n2,n3))
      endif
      if(iswrtyp >= 3)  then
         allocate (radiate%swup(n1,n2,n3))
         allocate (radiate%swdn(n1,n2,n3))
      endif
   
return
END SUBROUTINE alloc_radiate

!##############################################################################
Subroutine dealloc_radiate (radiate)

implicit none

   type (radiate_vars) :: radiate

   if (allocated(radiate%fthrd))    deallocate (radiate%fthrd)
   if (allocated(radiate%fthrdp))   deallocate (radiate%fthrdp)
   if (allocated(radiate%rshort))   deallocate (radiate%rshort)
   if (allocated(radiate%rlong))    deallocate (radiate%rlong)
   if (allocated(radiate%rlongup))  deallocate (radiate%rlongup)
   if (allocated(radiate%rlontop))  deallocate (radiate%rlontop)
   if (allocated(radiate%albedt))   deallocate (radiate%albedt)
   if (allocated(radiate%cosz))     deallocate (radiate%cosz)
   if (allocated(radiate%aodt))     deallocate (radiate%aodt)
   if (allocated(radiate%bext))     deallocate (radiate%bext)
   if (allocated(radiate%swup))     deallocate (radiate%swup)
   if (allocated(radiate%swdn))     deallocate (radiate%swdn)
   if (allocated(radiate%lwup))     deallocate (radiate%lwup)
   if (allocated(radiate%lwdn))     deallocate (radiate%lwdn)
   !GRL 2024-03-22 added variables for lw and sw heating rates
   if (allocated(radiate%fthrdsw))  deallocate (radiate%fthrdsw)
   if (allocated(radiate%fthrdlw))  deallocate (radiate%fthrdlw)

return
END SUBROUTINE dealloc_radiate

!##############################################################################
Subroutine filltab_radiate (radiate,radiatem,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (radiate_vars) :: radiate,radiatem
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(radiate%fthrd))  &
      CALL vtables2 (radiate%fthrd(1,1,1),radiatem%fthrd(1,1,1)  &
                 ,ng, npts, imean,  &
                 'FTHRD :3:anal:mpti')
   if (allocated(radiate%bext))  &
      CALL vtables2 (radiate%bext(1,1,1),radiatem%bext(1,1,1)  &
                 ,ng, npts, imean,  &
                 'BEXT :3:anal:mpti')
   if (allocated(radiate%swup))  &
      CALL vtables2 (radiate%swup(1,1,1),radiatem%swup(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SWUP :3:anal:mpti')
   if (allocated(radiate%swdn))  &
      CALL vtables2 (radiate%swdn(1,1,1),radiatem%swdn(1,1,1)  &
                 ,ng, npts, imean,  &
                 'SWDN :3:anal:mpti')
   if (allocated(radiate%lwup))  &
      CALL vtables2 (radiate%lwup(1,1,1),radiatem%lwup(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LWUP :3:anal:mpti')
   if (allocated(radiate%lwdn))  &
      CALL vtables2 (radiate%lwdn(1,1,1),radiatem%lwdn(1,1,1)  &
                 ,ng, npts, imean,  &
                 'LWDN :3:anal:mpti')
   !GRL 2024-03-22 added variables for lw and sw heating rates
   if (allocated(radiate%fthrdlw))  &
      CALL vtables2 (radiate%fthrdlw(1,1,1),radiatem%fthrdlw(1,1,1)  &
                 ,ng, npts, imean,  &
                 'FTHRDLW :3:anal:mpti')
   if (allocated(radiate%fthrdsw))  &
      CALL vtables2 (radiate%fthrdsw(1,1,1),radiatem%fthrdsw(1,1,1)  &
                 ,ng, npts, imean,  &
                 'FTHRDSW :3:anal:mpti')

   npts=n2*n3
   if (allocated(radiate%rshort))  &
      CALL vtables2 (radiate%rshort(1,1),radiatem%rshort(1,1)  &
                 ,ng, npts, imean,  &
                 'RSHORT :2:anal:mpti')
   if (allocated(radiate%rlong))  &
      CALL vtables2 (radiate%rlong(1,1),radiatem%rlong(1,1)  &
                 ,ng, npts, imean,  &
                 'RLONG :2:anal:mpti')
   if (allocated(radiate%rlongup))  &
      CALL vtables2 (radiate%rlongup(1,1),radiatem%rlongup(1,1)  &
                 ,ng, npts, imean,  &
                 'RLONGUP :2:anal:mpti')
   if (allocated(radiate%rlontop))  &
      CALL vtables2 (radiate%rlontop(1,1),radiatem%rlontop(1,1)  &
                 ,ng, npts, imean,  &
                 'RLONTOP :2:anal:mpti')
   if (allocated(radiate%albedt))  &
      CALL vtables2 (radiate%albedt(1,1),radiatem%albedt(1,1)  &
                 ,ng, npts, imean,  &
                 'ALBEDT :2:anal:mpti')
   if (allocated(radiate%cosz))  &
      CALL vtables2 (radiate%cosz(1,1),radiatem%cosz(1,1)  &
                 ,ng, npts, imean,  &
                 'COSZ :2:anal:mpti')
   if (allocated(radiate%aodt))  &
      CALL vtables2 (radiate%aodt(1,1),radiatem%aodt(1,1)  &
                 ,ng, npts, imean,  &
                 'AODT :2:anal:mpti')

return
END SUBROUTINE filltab_radiate

END MODULE mem_radiate
