!##############################################################################
Module mem_varinit

use grid_dims

implicit none

   ! Memory for varfile, history, and condensate nudging

   Type varinit_vars
   
      ! Variables to be dimensioned by (nzp,nxp,nyp)
   real, allocatable, dimension(:,:,:) :: &
                            varup,varvp,varpp,vartp,varrp  &
                           ,varuf,varvf,varpf,vartf,varrf  &
                           ,varwts &
                           ,varrph,varrfh,varcph,varcfh

   real, allocatable, dimension(:,:) :: &
            varsoilm1p,varsoilm2p,varsoilt1p,varsoilt2p,varsnowmp,varsnowdp &
           ,varsoilm1f,varsoilm2f,varsoilt1f,varsoilt2f,varsnowmf,varsnowdf

   End Type
   
   type (varinit_vars), allocatable :: varinit_g(:), varinitm_g(:)

   integer :: nud_type, nnudfiles, nnudfl, nudlat

   character(len=strl1), dimension(maxnudfiles) :: fnames_nud
   character(len=14) , dimension(maxnudfiles) :: itotdate_nud
   real, dimension(maxnudfiles) :: nud_times
   real :: tnudlat,tnudcent,tnudtop,znudtop,snudcent
   real :: wt_nudge_uv,wt_nudge_th,wt_nudge_pi,wt_nudge_rt
   real :: wt_nudge_g(maxgrds)
   real :: htime1, htime2
   
   integer :: igrid_match(maxgrds,maxgrds)

   !-------------------------------------------------------------------------------
   integer :: nud_cond
   real :: tcond_beg, tcond_end, wt_nudgec(maxgrds),t_nudge_rc

   !-------------------------------------------------------------------------------
   character(len=strl1), dimension(maxnudfiles) :: fnames_varf
   character(len=14) , dimension(maxnudfiles) :: itotdate_varf
   real, dimension(maxnudfiles) :: varf_times
   character(len=strl1) :: varfpfx
   real :: vtime1,vtime2,vwait1,vwaittot
   integer :: nvarffiles,nvarffl
   
   !-------------------------------------------------------------------------------
  
Contains

!##############################################################################
Subroutine alloc_varinit (varinit,n1,n2,n3)
   
use mem_grid
use mem_leaf, only:isoildat

implicit none

   type (varinit_vars) :: varinit
   integer, intent(in) :: n1,n2,n3

! Allocate arrays based on options (if necessary)
      
      if(nud_type == 1 .or. initial == 2) then
                         allocate (varinit%varup(n1,n2,n3))
                         allocate (varinit%varvp(n1,n2,n3))
                         allocate (varinit%varpp(n1,n2,n3))
                         allocate (varinit%vartp(n1,n2,n3))
                         allocate (varinit%varrp(n1,n2,n3))
                         allocate (varinit%varuf(n1,n2,n3))
                         allocate (varinit%varvf(n1,n2,n3))
                         allocate (varinit%varpf(n1,n2,n3))
                         allocate (varinit%vartf(n1,n2,n3))
                         allocate (varinit%varrf(n1,n2,n3))                      
                         allocate (varinit%varwts(n1,n2,n3))
      
                         if (nud_cond == 1) then
                          allocate (varinit%varcph(n1,n2,n3))
                          allocate (varinit%varcfh(n1,n2,n3))                      
                          allocate (varinit%varrph(n1,n2,n3))
                          allocate (varinit%varrfh(n1,n2,n3))                      
                         endif

                         if (isoildat == 1 .or. snudcent > 0.) then
                          allocate (varinit%varsoilm1p(n2,n3))
                          allocate (varinit%varsoilm2p(n2,n3))
                          allocate (varinit%varsoilt1p(n2,n3))
                          allocate (varinit%varsoilt2p(n2,n3))
                          allocate (varinit%varsnowmp(n2,n3))
                          allocate (varinit%varsnowdp(n2,n3))
                          allocate (varinit%varsoilm1f(n2,n3))
                          allocate (varinit%varsoilm2f(n2,n3))
                          allocate (varinit%varsoilt1f(n2,n3))
                          allocate (varinit%varsoilt2f(n2,n3))
                          allocate (varinit%varsnowmf(n2,n3))
                          allocate (varinit%varsnowdf(n2,n3))
                         endif
      endif

return
END SUBROUTINE alloc_varinit

!##############################################################################
Subroutine dealloc_varinit (varinit)

implicit none

   type (varinit_vars) :: varinit

   if (allocated(varinit%varup))     deallocate (varinit%varup)
   if (allocated(varinit%varvp))     deallocate (varinit%varvp)
   if (allocated(varinit%varpp))     deallocate (varinit%varpp)
   if (allocated(varinit%vartp))     deallocate (varinit%vartp)
   if (allocated(varinit%varrp))     deallocate (varinit%varrp)
   if (allocated(varinit%varuf))     deallocate (varinit%varuf)
   if (allocated(varinit%varvf))     deallocate (varinit%varvf)
   if (allocated(varinit%varpf))     deallocate (varinit%varpf)
   if (allocated(varinit%vartf))     deallocate (varinit%vartf)
   if (allocated(varinit%varrf))     deallocate (varinit%varrf)
   if (allocated(varinit%varwts))    deallocate (varinit%varwts)

   if (allocated(varinit%varcph))     deallocate (varinit%varcph)
   if (allocated(varinit%varcfh))     deallocate (varinit%varcfh)
   if (allocated(varinit%varrph))     deallocate (varinit%varrph)
   if (allocated(varinit%varrfh))     deallocate (varinit%varrfh)

   if (allocated(varinit%varsoilm1p)) deallocate (varinit%varsoilm1p)
   if (allocated(varinit%varsoilm2p)) deallocate (varinit%varsoilm2p)
   if (allocated(varinit%varsoilt1p)) deallocate (varinit%varsoilt1p)
   if (allocated(varinit%varsoilt2p)) deallocate (varinit%varsoilt2p)
   if (allocated(varinit%varsnowmp))  deallocate (varinit%varsnowmp)
   if (allocated(varinit%varsnowdp))  deallocate (varinit%varsnowdp)
   if (allocated(varinit%varsoilm1f)) deallocate (varinit%varsoilm1f)
   if (allocated(varinit%varsoilm2f)) deallocate (varinit%varsoilm2f)
   if (allocated(varinit%varsoilt1f)) deallocate (varinit%varsoilt1f)
   if (allocated(varinit%varsoilt2f)) deallocate (varinit%varsoilt2f)
   if (allocated(varinit%varsnowmf))  deallocate (varinit%varsnowmf)
   if (allocated(varinit%varsnowdf))  deallocate (varinit%varsnowdf)

return
END SUBROUTINE dealloc_varinit

!##############################################################################
Subroutine filltab_varinit (varinit,varinitm,imean,n1,n2,n3,ng)

use var_tables

implicit none

   type (varinit_vars) :: varinit,varinitm
   integer, intent(in) :: imean,n1,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n1*n2*n3

   if (allocated(varinit%varup))  &
      CALL vtables2 (varinit%varup(1,1,1),varinitm%varup(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARUP :3:mpti')
   if (allocated(varinit%varvp))  &
      CALL vtables2 (varinit%varvp(1,1,1),varinitm%varvp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARVP :3:mpti')
   if (allocated(varinit%varpp))  &
      CALL vtables2 (varinit%varpp(1,1,1),varinitm%varpp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARPP :3:mpti')
   if (allocated(varinit%vartp))  &
      CALL vtables2 (varinit%vartp(1,1,1),varinitm%vartp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARTP :3:mpti')
   if (allocated(varinit%varrp))  &
      CALL vtables2 (varinit%varrp(1,1,1),varinitm%varrp(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARRP :3:mpti')
   if (allocated(varinit%varuf))  &
      CALL vtables2 (varinit%varuf(1,1,1),varinitm%varuf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARUF :3:mpti')
   if (allocated(varinit%varvf))  &
      CALL vtables2 (varinit%varvf(1,1,1),varinitm%varvf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARVF :3:mpti')
   if (allocated(varinit%varpf))  &
      CALL vtables2 (varinit%varpf(1,1,1),varinitm%varpf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARPF :3:mpti')
   if (allocated(varinit%vartf))  &
      CALL vtables2 (varinit%vartf(1,1,1),varinitm%vartf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARTF :3:mpti')
   if (allocated(varinit%varrf))  &
      CALL vtables2 (varinit%varrf(1,1,1),varinitm%varrf(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARRF :3:mpti')
   if (allocated(varinit%varwts))  &
      CALL vtables2 (varinit%varwts(1,1,1),varinitm%varwts(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARWTS :3:mpti')

   if (allocated(varinit%varcph))  &
      CALL vtables2 (varinit%varcph(1,1,1),varinitm%varcph(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARCPH :3:mpti')
   if (allocated(varinit%varcfh))  &
      CALL vtables2 (varinit%varcfh(1,1,1),varinitm%varcfh(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARCFH :3:mpti')
   if (allocated(varinit%varrph))  &
      CALL vtables2 (varinit%varrph(1,1,1),varinitm%varrph(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARRPH :3:mpti')
   if (allocated(varinit%varrfh))  &
      CALL vtables2 (varinit%varrfh(1,1,1),varinitm%varrfh(1,1,1)  &
                 ,ng, npts, imean,  &
                 'VARRFH :3:mpti')

   npts=n2*n3

   if (allocated(varinit%varsoilm1p))  &
      CALL vtables2 (varinit%varsoilm1p(1,1),varinitm%varsoilm1p(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILM1P :2:mpti')
   if (allocated(varinit%varsoilm2p))  &
      CALL vtables2 (varinit%varsoilm2p(1,1),varinitm%varsoilm2p(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILM2P :2:mpti')
   if (allocated(varinit%varsoilt1p))  &
      CALL vtables2 (varinit%varsoilt1p(1,1),varinitm%varsoilt1p(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILT1P :2:mpti')
   if (allocated(varinit%varsoilt2p))  &
      CALL vtables2 (varinit%varsoilt2p(1,1),varinitm%varsoilt2p(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILT2P :2:mpti')
   if (allocated(varinit%varsnowmp))  &
      CALL vtables2 (varinit%varsnowmp(1,1),varinitm%varsnowmp(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSNOWMP :2:mpti')
   if (allocated(varinit%varsnowdp))  &
      CALL vtables2 (varinit%varsnowdp(1,1),varinitm%varsnowdp(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSNOWDP :2:mpti')
   if (allocated(varinit%varsoilm1f))  &
      CALL vtables2 (varinit%varsoilm1f(1,1),varinitm%varsoilm1f(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILM1F :2:mpti')
   if (allocated(varinit%varsoilm2f))  &
      CALL vtables2 (varinit%varsoilm2f(1,1),varinitm%varsoilm2f(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILM2F :2:mpti')
   if (allocated(varinit%varsoilt1f))  &
      CALL vtables2 (varinit%varsoilt1f(1,1),varinitm%varsoilt1f(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILT1F :2:mpti')
   if (allocated(varinit%varsoilt2f))  &
      CALL vtables2 (varinit%varsoilt2f(1,1),varinitm%varsoilt2f(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSOILT2F :2:mpti')
   if (allocated(varinit%varsnowmf))  &
      CALL vtables2 (varinit%varsnowmf(1,1),varinitm%varsnowmf(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSNOWMF :2:mpti')
   if (allocated(varinit%varsnowdf))  &
      CALL vtables2 (varinit%varsnowdf(1,1),varinitm%varsnowdf(1,1)  &
                 ,ng, npts, imean,  &
                 'VARSNOWDF :2:mpti')

return
END SUBROUTINE filltab_varinit

END MODULE mem_varinit
