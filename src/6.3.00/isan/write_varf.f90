!##############################################################################
Subroutine write_varf ()

use node_mod
use isan_coms
use mem_grid
use io_params
use hdf5_utils

implicit none

integer :: ng
character(len=strl1) :: locfn
character(len=3) :: csuff
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

do ng=1,nigrids
   nxyzp=nnxp(ng)*nnyp(ng)*nnzp(ng)
   nxyp =nnxp(ng)*nnyp(ng)
   write(csuff,'(a1,i1)') 'g',ng

   CALL makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                 ,ihour*100,'V',csuff,'h5 ')
 
   CALL shdf5_open (locfn,'W',iphdf5,h5_fid,iclobber)

   ! scalar values
   CALL shdf5_set_hs_select (1,'W',ng,mem_select,file_select,file_chunks)
   CALL shdf5_orec (h5_fid,iphdf5,'year'   ,mem_select,file_select &
       ,file_chunks,ivars=iyear)     
   CALL shdf5_orec (h5_fid,iphdf5,'month'  ,mem_select,file_select &
       ,file_chunks,ivars=imonth)    
   CALL shdf5_orec (h5_fid,iphdf5,'day'    ,mem_select,file_select &
       ,file_chunks,ivars=idate)     
   CALL shdf5_orec (h5_fid,iphdf5,'hour'   ,mem_select,file_select &
       ,file_chunks,ivars=ihour)     
   CALL shdf5_orec (h5_fid,iphdf5,'nx'     ,mem_select,file_select &
       ,file_chunks,ivars=nnxp(ng))  
   CALL shdf5_orec (h5_fid,iphdf5,'ny'     ,mem_select,file_select &
       ,file_chunks,ivars=nnyp(ng))  
   CALL shdf5_orec (h5_fid,iphdf5,'nz'     ,mem_select,file_select &
       ,file_chunks,ivars=nnzp(ng))  
   CALL shdf5_orec (h5_fid,iphdf5,'polelat',mem_select,file_select &
       ,file_chunks,rvars=polelat)        
   CALL shdf5_orec (h5_fid,iphdf5,'polelon',mem_select,file_select &
       ,file_chunks,rvars=polelon)       
   CALL shdf5_orec (h5_fid,iphdf5,'dx'     ,mem_select,file_select &
       ,file_chunks,rvars=deltaxn(ng))         
   CALL shdf5_orec (h5_fid,iphdf5,'dz'     ,mem_select,file_select &
       ,file_chunks,rvars=deltazn(ng))     
   CALL shdf5_orec (h5_fid,iphdf5,'dzrat'  ,mem_select,file_select &
       ,file_chunks,rvars=dzrat)     
   CALL shdf5_orec (h5_fid,iphdf5,'dzmax'  ,mem_select,file_select &
       ,file_chunks,rvars=dzmax)

   ! Atmospheric 3D vars
   CALL shdf5_set_hs_select (3,'W',ng,mem_select,file_select,file_chunks)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'UP',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_v,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'VP',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_p,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'PI',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_t,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'THETA',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)
   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_r,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'RV',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)

   CALL rearrange (nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_cond,rr_scr3)
   CALL shdf5_orec (h5_fid,iphdf5,'COND',mem_select,file_select &
       ,file_chunks,rvara=rr_scr3)
   
   ! Atmospheric 2D vars
   CALL shdf5_set_hs_select (2,'W',ng,mem_select,file_select,file_chunks)
   CALL vmissw (is_grids(ng)%rr_soilmoist1(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SOILMOIST1',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soilmoist2(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SOILMOIST2',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soiltemp1(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SOILTEMP1',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_soiltemp2(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SOILTEMP2',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_snowmass(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SNOWMASS',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)
   CALL vmissw (is_grids(ng)%rr_snowdepth(1,1),nxyp,rr_vt2da(1),1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'SNOWDEPTH',mem_select,file_select &
       ,file_chunks,rvara=rr_vt2da)

   CALL shdf5_close (h5_fid)

enddo

! Write the header file, only allow rammain to do this so the other
! processes don't collide during the write.
if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then
  CALL makefnam (locfn,varpfx,0,iyear,imonth,idate  &
                ,ihour*100,'V','$','tag')
  CALL rams_f_open (2,locfn,'FORMATTED','REPLACE','WRITE',iclobber)
  write(2,*) nigrids
  close(2)
endif

return
END SUBROUTINE write_varf
