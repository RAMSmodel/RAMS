!##############################################################################
Module mem_basic

implicit none

   Type basic_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                          up,uc,vp,vc,wp,wc,pp,pc  &
                         ,rv,theta,thp,rtp &
                         ,pi0,th0,rvt0,dn0,dn0u,dn0v &
                         ,wp_buoy_theta,wp_buoy_cond,wp_advdif

      ! These were used for testing perturbations of updated
      ! domain-mean base state quantities. Could be useful in testing.
   real, allocatable, dimension(:,:,:) :: th00,rvt00

      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                          fcoru,fcorv

   End Type
   
   type (basic_vars), allocatable :: basic_g(:), basicm_g(:)
  
Contains

!##############################################################################
Subroutine alloc_basic (basic,n1,n2,n3)

use micphys

implicit none

   type (basic_vars) :: basic
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)

      allocate (basic%up(n1,n2,n3))
      allocate (basic%uc(n1,n2,n3))
      allocate (basic%vp(n1,n2,n3))
      allocate (basic%vc(n1,n2,n3))
      allocate (basic%wp(n1,n2,n3))
      allocate (basic%wc(n1,n2,n3))
      allocate (basic%pp(n1,n2,n3))
      allocate (basic%pc(n1,n2,n3))
      allocate (basic%thp(n1,n2,n3))
      allocate (basic%rtp(n1,n2,n3))
      allocate (basic%rv(n1,n2,n3))
      allocate (basic%theta(n1,n2,n3))
      allocate (basic%pi0(n1,n2,n3))
      allocate (basic%th0(n1,n2,n3))
      allocate (basic%rvt0(n1,n2,n3))
      allocate (basic%dn0(n1,n2,n3))
      allocate (basic%dn0u(n1,n2,n3))
      allocate (basic%dn0v(n1,n2,n3))
      
      allocate (basic%fcoru(n2,n3))
      allocate (basic%fcorv(n2,n3))

      allocate (basic%th00(n1,n2,n3))
      allocate (basic%rvt00(n1,n2,n3))

      if(imbudget>=1) then
        allocate (basic%wp_buoy_theta(n1,n2,n3))
        allocate (basic%wp_buoy_cond(n1,n2,n3))
        allocate (basic%wp_advdif(n1,n2,n3))
      endif

return
END SUBROUTINE alloc_basic

!##############################################################################
Subroutine dealloc_basic (basic)

implicit none

   type (basic_vars) :: basic
   
   if (allocated(basic%up   ))    deallocate (basic%up   )
   if (allocated(basic%uc   ))    deallocate (basic%uc   )
   if (allocated(basic%vp   ))    deallocate (basic%vp   )
   if (allocated(basic%vc   ))    deallocate (basic%vc   )
   if (allocated(basic%wp   ))    deallocate (basic%wp   )
   if (allocated(basic%wc   ))    deallocate (basic%wc   )
   if (allocated(basic%pp   ))    deallocate (basic%pp   )
   if (allocated(basic%pc   ))    deallocate (basic%pc   )
   if (allocated(basic%thp  ))    deallocate (basic%thp  )
   if (allocated(basic%rtp  ))    deallocate (basic%rtp  )
   if (allocated(basic%rv   ))    deallocate (basic%rv   )
   if (allocated(basic%theta))    deallocate (basic%theta)
   if (allocated(basic%pi0  ))    deallocate (basic%pi0  )
   if (allocated(basic%th0  ))    deallocate (basic%th0  )
   if (allocated(basic%rvt0  ))   deallocate (basic%rvt0 )
   if (allocated(basic%dn0  ))    deallocate (basic%dn0  )
   if (allocated(basic%dn0u ))    deallocate (basic%dn0u )
   if (allocated(basic%dn0v ))    deallocate (basic%dn0v )
   if (allocated(basic%fcoru ))   deallocate (basic%fcoru )
   if (allocated(basic%fcorv ))   deallocate (basic%fcorv )

   if (allocated(basic%th00 ))    deallocate (basic%th00  )
   if (allocated(basic%rvt00 ))   deallocate (basic%rvt00 )

   if (allocated(basic%wp_buoy_theta)) deallocate (basic%wp_buoy_theta)
   if (allocated(basic%wp_buoy_cond))  deallocate (basic%wp_buoy_cond)
   if (allocated(basic%wp_advdif))     deallocate (basic%wp_advdif)

return
END SUBROUTINE dealloc_basic

!##############################################################################
Subroutine filltab_basic (basic,basicm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (basic_vars) :: basic,basicm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3
   
   if (allocated(basic%up))  &
      CALL vtables2 (basic%up(1,1,1),basicm%up(1,1,1)  &
                 ,ng, npts, imean,  &
                 'UP :3:anal:mpti')
   if (allocated(basic%vp))  &
      CALL vtables2 (basic%vp(1,1,1),basicm%vp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VP :3:anal:mpti')
   if (allocated(basic%wp))  &
      CALL vtables2 (basic%wp(1,1,1),basicm%wp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'WP :3:anal:mpti')
   if (allocated(basic%pp))  &
      CALL vtables2 (basic%pp(1,1,1),basicm%pp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PP :3:anal:mpti')
   if (allocated(basic%uc))  &
      CALL vtables2 (basic%uc(1,1,1),basicm%uc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'UC :3:anal:mpti')
   if (allocated(basic%vc))  &
      CALL vtables2 (basic%vc(1,1,1),basicm%vc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VC :3:anal:mpti')
   if (allocated(basic%wc))  &
      CALL vtables2 (basic%wc(1,1,1),basicm%wc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'WC :3:anal:mpti')
   if (allocated(basic%pc))  &
      CALL vtables2 (basic%pc(1,1,1),basicm%pc(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PC :3:anal:mpti')


   if (allocated(basic%thp)) &
      CALL vtables2 (basic%thp(1,1,1),basicm%thp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'THP :3:anal:mpti:mpt1')
   if (allocated(basic%rtp)) &
      CALL vtables2 (basic%rtp(1,1,1),basicm%rtp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RTP :3:anal:mpti:mpt1')

   if (allocated(basic%theta)) &
      CALL vtables2 (basic%theta(1,1,1),basicm%theta(1,1,1)  &
                 ,ng, npts, imean,  &
                 'THETA :3:anal:mpti')
   if (allocated(basic%rv)) &
      CALL vtables2 (basic%rv(1,1,1),basicm%rv(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RV :3:anal:mpti')
                 
   if (allocated(basic%pi0)) &
      CALL vtables2 (basic%pi0(1,1,1),basicm%pi0(1,1,1)  &
                 ,ng, npts, imean,  &
                 'PI0 :3:mpti')
   if (allocated(basic%th0)) &
      CALL vtables2 (basic%th0(1,1,1),basicm%th0(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TH0 :3:mpti')
   if (allocated(basic%rvt0)) &
      CALL vtables2 (basic%rvt0(1,1,1),basicm%rvt0(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RVT0 :3:mpti')
   if (allocated(basic%dn0)) &
      CALL vtables2 (basic%dn0(1,1,1),basicm%dn0(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DN0 :3:anal:mpti')
   if (allocated(basic%dn0u)) &
      CALL vtables2 (basic%dn0u(1,1,1),basicm%dn0u(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DN0U :3:mpti')
   if (allocated(basic%dn0v)) &
      CALL vtables2 (basic%dn0v(1,1,1),basicm%dn0v(1,1,1)  &
                 ,ng, npts, imean,  &
                 'DN0V :3:mpti')

   if (allocated(basic%th00)) &
      CALL vtables2 (basic%th00(1,1,1),basicm%th00(1,1,1)  &
                 ,ng, npts, imean,  &
                 'TH00 :3:mpti')
   if (allocated(basic%rvt00)) &
      CALL vtables2 (basic%rvt00(1,1,1),basicm%rvt00(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RVT00 :3:mpti')

   if (allocated(basic%wp_buoy_theta))  &
      CALL vtables2 (basic%wp_buoy_theta(1,1,1),basicm%wp_buoy_theta(1,1,1)  &
                 ,ng, npts, imean,  &
                 'WP_BUOY_THETA :3:anal:mpti')
   if (allocated(basic%wp_buoy_cond))  &
      CALL vtables2 (basic%wp_buoy_cond(1,1,1),basicm%wp_buoy_cond(1,1,1)  &
                 ,ng, npts, imean,  &
                 'WP_BUOY_COND :3:anal:mpti')
   if (allocated(basic%wp_advdif))  &
      CALL vtables2 (basic%wp_advdif(1,1,1),basicm%wp_advdif(1,1,1)  &
                 ,ng, npts, imean,  &
                 'WP_ADVDIF :3:anal:mpti')

   ! 2D CORIOLIS INFO FOR VTABLES
   npts=n2*n3
   if (allocated(basic%fcoru)) &
      CALL vtables2 (basic%fcoru(1,1),basicm%fcoru(1,1)  &
                 ,ng, npts, imean,  &
                 'FCORU :2:mpti')      
   if (allocated(basic%fcorv)) &
      CALL vtables2 (basic%fcorv(1,1),basicm%fcorv(1,1)  &
                 ,ng, npts, imean,  &
                 'FCORV :2:mpti')
 
return
END SUBROUTINE filltab_basic

END MODULE mem_basic
