!##############################################################################
Subroutine sfc_mem_alloc ()
 
 use mem_leaf
 use mem_grid
 use node_mod
 
 implicit none 
 
 integer :: ng
 
 ! Allocate Leaf type
 if(iprntstmt>=1 .and. print_msg)print*,'start leaf alloc'
 allocate(leaf_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_leaf (leaf_g(ng))
    CALL alloc_leaf (leaf_g(ng),mmxp(ng),mmyp(ng),nzg,nzs,npatch) 
 enddo
 
 ! Allocate grid variables data type. 
 if(iprntstmt>=1 .and. print_msg)print*,'start grid alloc'
 allocate(grid_g(ngrids))
 do ng=1,ngrids
    CALL dealloc_grid (grid_g(ng))
    CALL alloc_grid (grid_g(ng),mmxp(ng),mmyp(ng)) 
 enddo
 
 return
 END SUBROUTINE sfc_mem_alloc 

!##############################################################################
Subroutine rams_mem_alloc ()

use mem_all
use node_mod

implicit none 

integer :: ng,nv,imean

! First, depending on type of process, define grid point pointers correctly...

if(iprntstmt>=1 .and. print_msg)print*,' %% mem_alloc %%:',my_rams_num

!  If we are doing time-averaging for output, set flag ...
imean=0
if (avgtim /= 0.) imean=1

! Zero variable counts
num_var(1:maxgrds)=0
nvgrids=ngrids
num_scalar(1:maxgrds)=0

! Allocate Basic variables data type
if(iprntstmt>=1 .and. print_msg)print*,'start basic alloc'
allocate(basic_g(ngrids),basicm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_basic (basic_g(ng))
   CALL dealloc_basic (basicm_g(ng))
   CALL alloc_basic (basic_g(ng),mmzp(ng),mmxp(ng),mmyp(ng)) 
   if (imean == 1) then  
      CALL alloc_basic (basicm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng))
   elseif (imean == 0) then
      CALL alloc_basic (basicm_g(ng),1,1,1)
   endif
   
   CALL filltab_basic (basic_g(ng),basicm_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate Cuparm variables data type
if(iprntstmt>=1 .and. print_msg)print*,'start cuparm alloc'
allocate(cuparm_g(ngrids),cuparmm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_cuparm (cuparm_g(ng))
   CALL dealloc_cuparm (cuparmm_g(ng))
   CALL alloc_cuparm (cuparm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),ng) 
   if (imean == 1) then  
      CALL alloc_cuparm (cuparmm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),ng)
   elseif (imean == 0) then
      CALL alloc_cuparm (cuparmm_g(ng),1,1,1,ng)
   endif
   
   CALL filltab_cuparm (cuparm_g(ng),cuparmm_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate Leaf type

if(iprntstmt>=1 .and. print_msg)print*,'start leaf alloc'
allocate(leaf_g(ngrids),leafm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_leaf (leaf_g(ng))
   CALL dealloc_leaf (leafm_g(ng))
   CALL alloc_leaf (leaf_g(ng),mmxp(ng),mmyp(ng)  &
       ,nzg,nzs,npatch) 
   if (imean == 1) then  
      CALL alloc_leaf (leafm_g(ng),mmxp(ng),mmyp(ng)  &
          ,nzg,nzs,npatch)
   elseif (imean == 0) then
      CALL alloc_leaf (leafm_g(ng),1,1,1,1,1)
   endif
   
   CALL filltab_leaf (leaf_g(ng),leafm_g(ng),imean  &
          ,mmxp(ng),mmyp(ng),nzg,nzs,npatch,ng) 
enddo

! Allocate KPP Ocean mixed layer model type
if(iprntstmt>=1 .and. print_msg)print*,'start KPP-Ocean alloc'
allocate(kpp_const_fields%wmt(0:891,0:49))
allocate(kpp_const_fields%wst(0:891,0:49))
allocate(kpp_const_fields%tri(0:nkppz,0:1))
allocate(kpp_3d_fields(ngrids),kpp_3d_fieldsm(ngrids))
do ng=1,ngrids
   CALL dealloc_kpp (kpp_3d_fields(ng))
   CALL dealloc_kpp (kpp_3d_fieldsm(ng))
   CALL alloc_kpp (kpp_3d_fields(ng),mmxp(ng),mmyp(ng),nkppz,ng) 
   if (imean == 1) then  
      CALL alloc_kpp (kpp_3d_fieldsm(ng),mmxp(ng),mmyp(ng),nkppz,ng)
   elseif (imean == 0) then
      CALL alloc_kpp (kpp_3d_fieldsm(ng),1,1,1,ng)
   endif
   
   CALL filltab_kpp (kpp_3d_fields(ng),kpp_3d_fieldsm(ng),imean  &
     ,mmxp(ng),mmyp(ng),nkppz,ng) 
enddo

! Allocate SiB-RAMS
if(iprntstmt>=1 .and. print_msg)print*,'start Sib alloc'
allocate(sib_g(ngrids),sibm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_sib (sib_g(ng))
   CALL dealloc_sib (sibm_g(ng))
   CALL alloc_sib (sib_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),npatch)
   if (imean==1) then
      CALL alloc_sib (sibm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),npatch)
   elseif (imean==0) then
      CALL alloc_sib (sibm_g(ng),1,1,1,1)
   endif

   CALL filltab_sib (sib_g(ng),sibm_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),npatch,ng)
enddo

! Allocate Micro variables data type
if(iprntstmt>=1 .and. print_msg)print*,'start micro alloc'
allocate(micro_g(ngrids),microm_g(ngrids),pcp_tab(ngrids))
do ng=1,ngrids
   CALL dealloc_micro (micro_g(ng))
   CALL dealloc_micro (microm_g(ng))
   CALL dealloc_sedim (pcp_tab(ng))
   CALL alloc_micro (micro_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),nkr)
   CALL alloc_sedim (pcp_tab(ng),mmzp(ng))
   if (imean == 1) then  
      CALL alloc_micro (microm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),nkr)
   elseif (imean == 0) then
      CALL alloc_micro (microm_g(ng),1,1,1,1)
   endif
   
   CALL filltab_micro (micro_g(ng),microm_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate radiate variables data type
if(iprntstmt>=1 .and. print_msg)print*,'start radiate alloc'
allocate(radiate_g(ngrids),radiatem_g(ngrids))
do ng=1,ngrids
   CALL dealloc_radiate (radiate_g(ng))
   CALL dealloc_radiate (radiatem_g(ng))
   CALL alloc_radiate (radiate_g(ng),mmzp(ng),mmxp(ng),mmyp(ng)) 
   if (imean == 1) then  
      CALL alloc_radiate (radiatem_g(ng),mmzp(ng),mmxp(ng),mmyp(ng))
   elseif (imean == 0) then
      CALL alloc_radiate (radiatem_g(ng),1,1,1)
   endif
   
   CALL filltab_radiate (radiate_g(ng),radiatem_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate turb variables data type
if(iprntstmt>=1 .and. print_msg)print*,'start turb alloc'
allocate(turb_g(ngrids),turbm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_turb (turb_g(ng))
   CALL dealloc_turb (turbm_g(ng))
   CALL alloc_turb (turb_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),ngrids) 
   if (imean == 1) then  
      CALL alloc_turb (turbm_g(ng),mmzp(ng),mmxp(ng),mmyp(ng),ngrids)
   elseif (imean == 0) then
      CALL alloc_turb (turbm_g(ng),1,1,1,ngrids)
   endif
   
   CALL filltab_turb (turb_g(ng),turbm_g(ng),imean  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate varinit variables data type. 
!    These do not need "mean" type ever.
if(iprntstmt>=1 .and. print_msg)print*,'start varinit alloc'
allocate(varinit_g(ngrids),varinitm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_varinit (varinit_g(ng))
   CALL dealloc_varinit (varinitm_g(ng))
   CALL alloc_varinit (varinit_g(ng),mmzp(ng),mmxp(ng),mmyp(ng)) 
   CALL alloc_varinit (varinitm_g(ng),1,1,1)
   
   CALL filltab_varinit (varinit_g(ng),varinitm_g(ng),0  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate oda variables data type. 
!    These do not need "mean" type ever.
if(iprntstmt>=1 .and. print_msg)print*,'start oda alloc'
allocate(oda_g(ngrids),odam_g(ngrids))
do ng=1,ngrids
   CALL dealloc_oda (oda_g(ng))
   CALL dealloc_oda (odam_g(ng))
   CALL alloc_oda (oda_g(ng),mmzp(ng),mmxp(ng),mmyp(ng)) 
   CALL alloc_oda (odam_g(ng),1,1,1)
   
   CALL filltab_oda (oda_g(ng),odam_g(ng),0  &
          ,mmzp(ng),mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate grid variables data type. 

if(iprntstmt>=1 .and. print_msg)print*,'start grid alloc'
allocate(grid_g(ngrids),gridm_g(ngrids))
do ng=1,ngrids
   CALL dealloc_grid (grid_g(ng))
   CALL dealloc_grid (gridm_g(ng))
   CALL alloc_grid (grid_g(ng),mmxp(ng),mmyp(ng)) 
   if (imean == 1) then
      CALL alloc_grid (gridm_g(ng),mmxp(ng),mmyp(ng)) 
   elseif (imean == 0) then
      CALL alloc_grid (gridm_g(ng),1,1)
   endif

   CALL filltab_grid (grid_g(ng),gridm_g(ng),imean,mmxp(ng),mmyp(ng),ng) 
enddo

! Allocate any added Scalar types

if(iprntstmt>=1 .and. print_msg)print*,'start scalar alloc'
allocate(tracer_g(itracer,ngrids),tracerm_g(itracer,ngrids))
do ng=1,ngrids
   CALL dealloc_tracer (tracer_g(:,ng),ng)
   CALL dealloc_tracer (tracerm_g(:,ng),ng)
   CALL alloc_tracer (tracer_g(:,ng),mmzp(ng),mmxp(ng),mmyp(ng),ng)
   if (imean == 1) then
    CALL alloc_tracer (tracerm_g(:,ng),mmzp(ng),mmxp(ng),mmyp(ng),ng)
   elseif (imean == 0) then
    CALL alloc_tracer (tracerm_g(:,ng),1,1,1,ng)
   endif

   CALL filltab_tracer (tracer_g(:,ng),tracerm_g(:,ng),imean  &
                       ,mmzp(ng),mmxp(ng),mmyp(ng),ng)
enddo

! Allocate Tendency data type,  filltab_tendency is responsible 
!   for filling the main atmospheric model variables in the scalar table,
!   so make sure to call any routines that define scalar variables first.

! Assuming same scalars on all grids!!!!!

if(iprntstmt>=1 .and. print_msg)print*,'start tendency alloc'
CALL dealloc_tend (ngrids)
CALL alloc_tend (mmzp,mmxp,mmyp,ngrids)
do ng=1,ngrids
   CALL filltab_tend (basic_g(ng),micro_g(ng),turb_g(ng),sib_g(ng) &
      ,tracer_g(:,ng),ng) 
enddo

! Allocate Scratch data type, This also fills the max's that are needed
!    by nesting stuff.

if(iprntstmt>=1 .and. print_msg)print*,'start scratch alloc'
CALL dealloc_scratch ()
CALL alloc_scratch (mmzp,mmxp,mmyp,nnzp,nnxp,nnyp,ngrids  &
                  ,nzg,nzs,npatch,maxnxp,maxnyp,maxnzp)

! Allocate nested boundary interpolation arrays. All grids will be allocated.

if(iprntstmt>=1 .and. print_msg)print*,'start nestb alloc'
do ng=1,ngrids
   if(nxtnest(ng) == 0 ) then
      CALL alloc_nestb (ng,1,1,1)
   else
      CALL alloc_nestb (ng,nnxp(ng),nnyp(ng),nnzp(ng))
   endif
enddo

! Set "Lite" variable flags according to namelist input LITE_VARS.

CALL lite_varset ()


! Set ALL variables in the vtab_r variable table to zero by default. These
!  are variables processed in the filltab_* routines with a call to vtables2.
!  This does NOT include scratch arrays, tendencies, or mean arrays.

do ng = 1, ngrids
 do nv = 1,num_var(ng)
  if(iprntstmt>=1 .and. print_msg) &
     print*,'Zeroing out array:',ng,nv,vtab_r(nv,ng)%name
  CALL azero (vtab_r(nv,ng)%npts,vtab_r(nv,ng)%var_p)
 enddo
enddo

if(iprntstmt>=1 .and. print_msg)print*,'end alloc'

return
END SUBROUTINE rams_mem_alloc

