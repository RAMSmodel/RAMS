!##############################################################################
Module mem_cuparm

use grid_dims

implicit none

   Type cuparm_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          thsrc,rtsrc,thsrcf,rtsrcf,thsrcp,rtsrcp &
                         ,rcsrc,rrsrc,rssrc,rpsrc,w0avg,w0avglt

      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          aconpr,conprr,conprrp,conprrf &
                         ,nca,convgo
   

   End Type

   type (cuparm_vars), allocatable :: cuparm_g(:), cuparmm_g(:)

   integer, dimension(maxgrds) :: nnqparm
   real :: wcldbs,confrq
  
Contains

!##############################################################################
Subroutine alloc_cuparm (cuparm,n1,n2,n3,ng)

implicit none

   type (cuparm_vars) :: cuparm
   integer, intent(in) :: n1,n2,n3,ng

! Allocate arrays based on options (if necessary)
      
      if( nnqparm(ng)>= 1)  then
                         allocate (cuparm%thsrc(n1,n2,n3))
                         allocate (cuparm%rtsrc(n1,n2,n3))
                         allocate (cuparm%aconpr(n2,n3))
                         allocate (cuparm%conprr(n2,n3))
                         if (nnqparm(ng) == 2) then
                            allocate (cuparm%rcsrc(n1,n2,n3))
                            allocate (cuparm%rrsrc(n1,n2,n3))
                            allocate (cuparm%rssrc(n1,n2,n3))
                            allocate (cuparm%rpsrc(n1,n2,n3))
                            allocate (cuparm%w0avg(n1,n2,n3))
                            allocate (cuparm%w0avglt(n1,n2,n3))
                            allocate (cuparm%nca(n2,n3))
                            allocate (cuparm%convgo(n2,n3))
                         endif                           
      endif
                         
return
END SUBROUTINE alloc_cuparm

!##############################################################################
Subroutine dealloc_cuparm (cuparm)

implicit none

   type (cuparm_vars) :: cuparm

   if (allocated(cuparm%thsrc))    deallocate (cuparm%thsrc)
   if (allocated(cuparm%rtsrc))    deallocate (cuparm%rtsrc)
   if (allocated(cuparm%thsrcp))    deallocate (cuparm%thsrcp)
   if (allocated(cuparm%rtsrcp))    deallocate (cuparm%rtsrcp)
   if (allocated(cuparm%thsrcf))    deallocate (cuparm%thsrcf)
   if (allocated(cuparm%rtsrcf))    deallocate (cuparm%rtsrcf)
   if (allocated(cuparm%aconpr))   deallocate (cuparm%aconpr)
   if (allocated(cuparm%conprr))   deallocate (cuparm%conprr)
   if (allocated(cuparm%conprrp))   deallocate (cuparm%conprrp)
   if (allocated(cuparm%conprrf))   deallocate (cuparm%conprrf)

   if (allocated(cuparm%rcsrc)  )   deallocate (cuparm%rcsrc)
   if (allocated(cuparm%rrsrc)  )   deallocate (cuparm%rrsrc)
   if (allocated(cuparm%rssrc)  )   deallocate (cuparm%rssrc)
   if (allocated(cuparm%rpsrc)  )   deallocate (cuparm%rpsrc)
   if (allocated(cuparm%w0avg)  )   deallocate (cuparm%w0avg)
   if (allocated(cuparm%w0avglt))   deallocate (cuparm%w0avglt)
   if (allocated(cuparm%nca)    )   deallocate (cuparm%nca)
   if (allocated(cuparm%convgo) )   deallocate (cuparm%convgo)
   
return
END SUBROUTINE dealloc_cuparm

!##############################################################################
Subroutine filltab_cuparm (cuparm,cuparmm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (cuparm_vars) :: cuparm,cuparmm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(cuparm%thsrc))  &
      CALL vtables2 (cuparm%thsrc(1,1,1),cuparmm%thsrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'THSRC :3:anal:mpti')
   if (allocated(cuparm%rtsrc))  &
      CALL vtables2 (cuparm%rtsrc(1,1,1),cuparmm%rtsrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RTSRC :3:anal:mpti')
   if (allocated(cuparm%thsrcp))  &
      CALL vtables2 (cuparm%thsrcp(1,1,1),cuparmm%thsrcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'THSRCP :3:mpti:')
   if (allocated(cuparm%rtsrcp))  &
      CALL vtables2 (cuparm%rtsrcp(1,1,1),cuparmm%rtsrcp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RTSRCP :3:mpti:')
   if (allocated(cuparm%thsrcf))  &
      CALL vtables2 (cuparm%thsrcf(1,1,1),cuparmm%thsrcf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'THSRCF :3:mpti:')
   if (allocated(cuparm%rtsrcf))  &
      CALL vtables2 (cuparm%rtsrcf(1,1,1),cuparmm%rtsrcf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RTSRCF :3:mpti:')

   npts=n2*n3
   if (allocated(cuparm%aconpr))  &
      CALL vtables2 (cuparm%aconpr(1,1),cuparmm%aconpr(1,1)  &
                 ,ng, npts, imean,  &
                 'ACONPR :2:anal:mpti')
   if (allocated(cuparm%conprr))  &
      CALL vtables2 (cuparm%conprr(1,1),cuparmm%conprr(1,1)  &
                 ,ng, npts, imean,  &
                 'CONPRR :2:anal:mpti')
   if (allocated(cuparm%conprrp))  &
      CALL vtables2 (cuparm%conprrp(1,1),cuparmm%conprrp(1,1)  &
                 ,ng, npts, imean,  &
                 'CONPRRP :2:mpti')
   if (allocated(cuparm%conprrf))  &
      CALL vtables2 (cuparm%conprrf(1,1),cuparmm%conprrf(1,1)  &
                 ,ng, npts, imean,  &
                 'CONPRRF :2:mpti')

   npts=n1*n2*n3
   if (allocated(cuparm%rcsrc)  )  &   
      CALL vtables2 (cuparm%rcsrc(1,1,1),cuparmm%rcsrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RCSRC :3:anal:mpti')
   if (allocated(cuparm%rrsrc)  )  &   
      CALL vtables2 (cuparm%rrsrc(1,1,1),cuparmm%rrsrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RRSRC :3:anal:mpti')
   if (allocated(cuparm%rssrc)  )  &   
      CALL vtables2 (cuparm%rssrc(1,1,1),cuparmm%rssrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RSSRC :3:anal:mpti')
   if (allocated(cuparm%rpsrc)  )  &   
      CALL vtables2 (cuparm%rpsrc(1,1,1),cuparmm%rpsrc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RPSRC :3:anal:mpti')
   if (allocated(cuparm%w0avg)  )  &   
      CALL vtables2 (cuparm%w0avg(1,1,1),cuparmm%w0avg(1,1,1)  &
                 ,ng, npts, imean,  &
                 'W0AVG :3:anal:mpti')
   if (allocated(cuparm%w0avglt))  &   
      CALL vtables2 (cuparm%w0avglt(1,1,1),cuparmm%w0avglt(1,1,1)  &
                 ,ng, npts, imean,  &
                 'W0AVGLT :3:anal:mpti')
   npts=n2*n3
   if (allocated(cuparm%nca)    )   &  
      CALL vtables2 (cuparm%nca(1,1),cuparmm%nca(1,1)  &
                 ,ng, npts, imean,  &
                 'NCA :2:anal:mpti')
   if (allocated(cuparm%convgo) )  &   
      CALL vtables2 (cuparm%convgo(1,1),cuparmm%convgo(1,1)  &
                 ,ng, npts, imean,  &
                 'CONVGO :2:anal:mpti')
                 
return
END SUBROUTINE filltab_cuparm

END MODULE mem_cuparm
