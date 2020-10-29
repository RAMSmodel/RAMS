!##############################################################################
Module mem_tracer

use grid_dims

implicit none

   Type tracer_vars

      ! Variables to be dimensioned by (nzp,nxp,nyp)
      real, allocatable, dimension(:,:,:) :: tracerp
      real, allocatable, dimension(:) :: tracert

   End Type
   
   ! tracerp allocated by (maxsclr,ngrids)
   type (tracer_vars), allocatable :: tracer_g(:,:), tracerm_g(:,:)

   integer :: itracer,itrachist

Contains

!##############################################################################
Subroutine alloc_tracer (tracer,n1,n2,n3,ng)

implicit none

   type (tracer_vars) :: tracer(*)
   integer, intent(in) :: n1,n2,n3,ng
   integer :: nsc

! Allocate arrays based on options (if necessary)

   do nsc=1,itracer
      allocate (tracer(nsc)%tracerp(n1,n2,n3))
   enddo

return
END SUBROUTINE alloc_tracer

!##############################################################################
Subroutine dealloc_tracer (tracer,ng)

implicit none

   type (tracer_vars) :: tracer(*)

   integer, intent(in) :: ng
   integer :: nsc

   do nsc=1,itracer
     if (allocated(tracer(nsc)%tracerp))  deallocate (tracer(nsc)%tracerp)
   enddo

return
END SUBROUTINE dealloc_tracer

!##############################################################################
Subroutine filltab_tracer (tracer,tracerm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (tracer_vars) :: tracer(*),tracerm(*)
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: nsc,npts
   character (len=10) :: sname

! Fill arrays into variable tables

   npts=n1*n2*n3
   do nsc=1,itracer
     if (allocated(tracer(nsc)%tracerp)) then
      write(sname,'(a7,i3.3)') 'TRACERP',nsc
      CALL vtables2 (tracer(nsc)%tracerp(1,1,1),tracerm(nsc)%tracerp(1,1,1) &
         ,ng, npts, imean, sname//' :3:anal:mpti:mpt1')
     endif
   enddo

return
END SUBROUTINE filltab_tracer

END MODULE mem_tracer
