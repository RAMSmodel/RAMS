!##############################################################################
Subroutine masterput_oda ()

use grid_dims
use mem_oda
use node_mod

implicit none

!   +------------------------------------------------------------------
!   ! This routine gives obs data to the nodes for ODA. All obs are sent 
!   !    to all nodes.
!   +------------------------------------------------------------------

real, allocatable :: buff(:)
integer :: nwords, nm, ns

nwords = 12 + 7 * maxodagrids + 6
allocate (buff(nwords))

CALL par_init_put (buff,nwords)

! Namelist info
CALL par_put_float (wt_oda_grid,maxodagrids)
CALL par_put_float (frqoda,1)
CALL par_put_float (todabeg,1)
CALL par_put_float (todaend,1)
CALL par_put_float (tnudoda,1)
CALL par_put_float (wt_oda_uv,1)
CALL par_put_float (wt_oda_th,1)
CALL par_put_float (wt_oda_pi,1)
CALL par_put_float (wt_oda_rt,1)
CALL par_put_float (oda_sfc_til,1)
CALL par_put_float (oda_sfc_tel,1)
CALL par_put_float (oda_upa_til,1)
CALL par_put_float (oda_upa_tel,1)
CALL par_put_float (roda_sfce,maxodagrids)
CALL par_put_float (roda_sfc0,maxodagrids)
CALL par_put_float (roda_upae,maxodagrids)
CALL par_put_float (roda_upa0,maxodagrids)
CALL par_put_float (roda_zfac,maxodagrids)
CALL par_put_float (roda_hgt,maxodagrids)

! Sizes for arrays
CALL par_put_int (num_oda_sfc,1)
CALL par_put_int (nsfcfiles,1)
CALL par_put_int (maxtimes_sfc,1)
CALL par_put_int (num_oda_upa,1)
CALL par_put_int (nupafiles,1)
CALL par_put_int (maxtimes_upa,1)


do nm=1,nmachs
   CALL par_send (machnum(nm),ODA_TAG_BASE+1)
enddo

deallocate (buff)

! Surface obs

nwords= num_oda_sfc*( 10 + 3*maxodagrids + 8*maxtimes_sfc)
allocate (buff(nwords))

CALL par_init_put (buff,nwords)


do ns=1,num_oda_sfc
   CALL par_put_char  (oda_sfc_info(ns)%id,8)
   CALL par_put_int   (oda_sfc_info(ns)%intid,1)
   CALL par_put_int   (oda_sfc_info(ns)%ntimes,1)
   CALL par_put_int   (oda_sfc_info(ns)%iactive(1),maxodagrids)
   CALL par_put_float (oda_sfc_info(ns)%xista(1),maxodagrids)
   CALL par_put_float (oda_sfc_info(ns)%xjsta(1),maxodagrids)
   CALL par_put_float (oda_sfc_info(ns)%xlat,1)
   CALL par_put_float (oda_sfc_info(ns)%xlon,1)
   CALL par_put_float (oda_sfc_info(ns)%xsta,1)
   CALL par_put_float (oda_sfc_info(ns)%ysta,1)
   CALL par_put_float (oda_sfc_info(ns)%stopo,1)
   
   CALL par_put_float (oda_sfc_obs(ns)%time(1) ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%temp(1) ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%dewpt(1),maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%us(1)   ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%vs(1)   ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%u (1)   ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%v (1)   ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%ps(1)   ,maxtimes_sfc)
   CALL par_put_float (oda_sfc_obs(ns)%psref   ,1)
enddo

do nm=1,nmachs
   CALL par_send (machnum(nm),ODA_TAG_BASE+2)
enddo

deallocate (buff)

! Upper air obs

nwords= num_oda_upa*( 9 + 3*maxodagrids + 3*maxtimes_upa  &
                                + 9*maxupalevs*maxtimes_upa)
allocate (buff(nwords))

CALL par_init_put (buff,nwords)


do ns=1,num_oda_upa
   CALL par_put_char  (oda_upa_info(ns)%id,8)
   CALL par_put_int   (oda_upa_info(ns)%intid,1)
   CALL par_put_int   (oda_upa_info(ns)%ntimes,1)
   CALL par_put_int   (oda_upa_info(ns)%iactive(1),maxodagrids)
   CALL par_put_float (oda_upa_info(ns)%xista(1),maxodagrids)
   CALL par_put_float (oda_upa_info(ns)%xjsta(1),maxodagrids)
   CALL par_put_float (oda_upa_info(ns)%xlat,1)
   CALL par_put_float (oda_upa_info(ns)%xlon,1)
   CALL par_put_float (oda_upa_info(ns)%xsta,1)
   CALL par_put_float (oda_upa_info(ns)%ysta,1)
   CALL par_put_float (oda_upa_info(ns)%stopo,1)

   CALL par_put_float (oda_upa_obs(ns)%time(1),maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%lp  (1),maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%lz  (1),maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%theta(1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%rv   (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%us   (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%vs   (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%zz   (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%u    (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%v    (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%pi   (1,1),maxupalevs*maxtimes_upa)
   CALL par_put_float (oda_upa_obs(ns)%zgeo (1,1),maxupalevs*maxtimes_upa)
enddo

do nm=1,nmachs
   CALL par_send (machnum(nm),ODA_TAG_BASE+3)
enddo

deallocate (buff)

return
END SUBROUTINE masterput_oda

!##############################################################################
Subroutine nodeget_oda ()

use node_mod
use mem_oda

implicit none

real, allocatable :: buff(:)
integer :: nwords,ibytes,msgtype,ihostnum, ns

nwords = 12 + 7 * maxodagrids + 6
allocate (buff(nwords))

CALL par_get_new (buff,nwords,ODA_TAG_BASE+1,ibytes,msgtype,ihostnum)

! Namelist info
CALL par_get_float (wt_oda_grid,maxodagrids)
CALL par_get_float (frqoda,1)
CALL par_get_float (todabeg,1)
CALL par_get_float (todaend,1)
CALL par_get_float (tnudoda,1)
CALL par_get_float (wt_oda_uv,1)
CALL par_get_float (wt_oda_th,1)
CALL par_get_float (wt_oda_pi,1)
CALL par_get_float (wt_oda_rt,1)
CALL par_get_float (oda_sfc_til,1)
CALL par_get_float (oda_sfc_tel,1)
CALL par_get_float (oda_upa_til,1)
CALL par_get_float (oda_upa_tel,1)
CALL par_get_float (roda_sfce,maxodagrids)
CALL par_get_float (roda_sfc0,maxodagrids)
CALL par_get_float (roda_upae,maxodagrids)
CALL par_get_float (roda_upa0,maxodagrids)
CALL par_get_float (roda_zfac,maxodagrids)
CALL par_get_float (roda_hgt,maxodagrids)

! Sizes for arrays
CALL par_get_int (num_oda_sfc,1)
CALL par_get_int (nsfcfiles,1)
CALL par_get_int (maxtimes_sfc,1)
CALL par_get_int (num_oda_upa,1)
CALL par_get_int (nupafiles,1)
CALL par_get_int (maxtimes_upa,1)

deallocate (buff)

! Allocate main obs arrays

CALL oda_obs_alloc ()

! Surface obs

nwords= num_oda_sfc*( 10 + 3*maxodagrids + 8*maxtimes_sfc)
allocate (buff(nwords))

CALL par_get_new (buff,nwords,ODA_TAG_BASE+2,ibytes,msgtype,ihostnum)

do ns=1,num_oda_sfc
   CALL par_get_char  (oda_sfc_info(ns)%id,8)
   CALL par_get_int   (oda_sfc_info(ns)%intid,1)
   CALL par_get_int   (oda_sfc_info(ns)%ntimes,1)
   CALL par_get_int   (oda_sfc_info(ns)%iactive(1),maxodagrids)
   CALL par_get_float (oda_sfc_info(ns)%xista(1),maxodagrids)
   CALL par_get_float (oda_sfc_info(ns)%xjsta(1),maxodagrids)
   CALL par_get_float (oda_sfc_info(ns)%xlat,1)
   CALL par_get_float (oda_sfc_info(ns)%xlon,1)
   CALL par_get_float (oda_sfc_info(ns)%xsta,1)
   CALL par_get_float (oda_sfc_info(ns)%ysta,1)
   CALL par_get_float (oda_sfc_info(ns)%stopo,1)
   
   CALL par_get_float (oda_sfc_obs(ns)%time(1) ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%temp(1) ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%dewpt(1),maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%us(1)   ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%vs(1)   ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%u (1)   ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%v (1)   ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%ps(1)   ,maxtimes_sfc)
   CALL par_get_float (oda_sfc_obs(ns)%psref   ,1)
enddo

deallocate (buff)

! Upper air obs

nwords= num_oda_upa*( 9 + 3*maxodagrids + 3*maxtimes_upa  &
                                + 9*maxupalevs*maxtimes_upa)
allocate (buff(nwords))

CALL par_get_new (buff,nwords,ODA_TAG_BASE+3,ibytes,msgtype,ihostnum)

do ns=1,num_oda_upa
   CALL par_get_char  (oda_upa_info(ns)%id,8)
   CALL par_get_int   (oda_upa_info(ns)%intid,1)
   CALL par_get_int   (oda_upa_info(ns)%ntimes,1)
   CALL par_get_int   (oda_upa_info(ns)%iactive(1),maxodagrids)
   CALL par_get_float (oda_upa_info(ns)%xista(1),maxodagrids)
   CALL par_get_float (oda_upa_info(ns)%xjsta(1),maxodagrids)
   CALL par_get_float (oda_upa_info(ns)%xlat,1)
   CALL par_get_float (oda_upa_info(ns)%xlon,1)
   CALL par_get_float (oda_upa_info(ns)%xsta,1)
   CALL par_get_float (oda_upa_info(ns)%ysta,1)
   CALL par_get_float (oda_upa_info(ns)%stopo,1)

   CALL par_get_float (oda_upa_obs(ns)%time(1),maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%lp  (1),maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%lz  (1),maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%theta(1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%rv   (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%us   (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%vs   (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%zz   (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%u    (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%v    (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%pi   (1,1),maxupalevs*maxtimes_upa)
   CALL par_get_float (oda_upa_obs(ns)%zgeo (1,1),maxupalevs*maxtimes_upa)
enddo

deallocate (buff)

return
END SUBROUTINE nodeget_oda
