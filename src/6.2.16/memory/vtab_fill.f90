!##############################################################################
Subroutine vtables2 (var,varm,ng,npts,imean,tabstr)

use var_tables
      
implicit none

   real, target :: var,varm
   integer, intent(in) :: ng,npts,imean
   character (len=*), intent(in) :: tabstr
   character (len=1) :: toksep=':'
   character (len=32) :: tokens(10),ctab
   integer :: ntok,nt,nv
     
   CALL tokenize1 (tabstr,tokens,ntok,toksep)
   
   num_var(ng)=num_var(ng)+1
   nv=num_var(ng)

   vtab_r(nv,ng)%var_p => var
   vtab_r(nv,ng)%var_m => varm

  
   vtab_r(nv,ng)%name=tokens(1)
   vtab_r(nv,ng)%npts=npts
   read(tokens(2),*) vtab_r(nv,ng)%idim_type
   !print*,'tab:',nv,ng,vtab_r(nv,ng)%name ,vtab_r(nv,ng)%npts
     
   vtab_r(nv,ng)%ianal=0
   vtab_r(nv,ng)%imean=imean
   vtab_r(nv,ng)%ilite=0
   vtab_r(nv,ng)%impti=0
   vtab_r(nv,ng)%impt1=0
   vtab_r(nv,ng)%irecycle_sfc=0

   do nt=3,ntok
      ctab=tokens(nt)         
      
      if(ctab == 'anal' ) then
         vtab_r(nv,ng)%ianal=1
      elseif(ctab == 'lite' ) then
         vtab_r(nv,ng)%ilite=1
      elseif(ctab == 'mpti' ) then
         vtab_r(nv,ng)%impti=1
      elseif(ctab == 'mpt1' ) then
         vtab_r(nv,ng)%impt1=1
      elseif(ctab == 'recycle_sfc' ) then
         vtab_r(nv,ng)%irecycle_sfc=1
      else
         print*, 'Illegal table specification for var:', tokens(1),ctab
         stop 'bad var table'
      endif

   enddo
  
return
END SUBROUTINE vtables2

!##############################################################################
Subroutine lite_varset ()

use var_tables
use io_params
use mem_grid, only:print_msg

implicit none

integer :: nv,ng,nvl,ifound

! Loop over each variable input in namelist "LITE_VARS" and set
!   lite flag in var_tables

do ng = 1,nvgrids   
   vtab_r(1:num_var(ng),ng)%ilite = 0
enddo

do nvl=1,nlite_vars
   ifound=0
      
   do ng=1,nvgrids
   
      do nv=1,num_var(ng)
   
         if (vtab_r(nv,ng)%name == lite_vars(nvl) ) then
            vtab_r(nv,ng)%ilite = 1
            ifound=1
         endif
         
      enddo
      
   enddo

   if(print_msg)then   
    if(ifound == 0) then
      print*,'!---------------------------------------------------------'
      print*,'! LITE_VARS variable does not exist in main variable table'
      print*,'!    variable name-->',lite_vars(nvl),'<--'
      print*,'!---------------------------------------------------------'
    else
      print*,'!---------------------------------------------------------'
      print*,'! LITE_VARS variable added--->',trim(lite_vars(nvl))
      print*,'!---------------------------------------------------------'
    endif
   endif
   
enddo

return
END SUBROUTINE lite_varset

!##############################################################################
Subroutine vtables_scalar (varp,vart,ng,tabstr)

use var_tables
use mem_grid, only:iprntstmt,print_msg
      
implicit none

   real, target :: varp,vart
   integer, intent(in) :: ng
   character (len=*), intent(in) :: tabstr
   character (len=1) :: toksep=':'
   character (len=32) :: tokens(10),cname
   integer :: ntok,nv

   CALL tokenize1 (tabstr,tokens,ntok,toksep)
   cname=tokens(1)
   
!  Fill in existing table slot

   num_scalar(ng)=num_scalar(ng)+1
   nv=num_scalar(ng)
   scalar_tab(nv,ng)%name = cname

   scalar_tab(nv,ng)%var_p => varp
   scalar_tab(nv,ng)%var_t => vart
   if(iprntstmt>=1 .and. print_msg)print*,'Scalars: ',ng,nv,scalar_tab(nv,ng)%name
  
return
END SUBROUTINE vtables_scalar
