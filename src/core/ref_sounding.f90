!##############################################################################
Module ref_sounding

use grid_dims

implicit none

   integer :: iref,jref,nzref
   real :: topref,divls

   !Base state variables in 1D (constant in time)
   real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn,o3ref,wsub
   
   !Base state forcing variables
   character(strl1) :: sound_file,forcingfile
   integer :: iugforce
   real, allocatable, dimension(:) :: forc_lev, forc_time, forc_ts
   real, allocatable, dimension(:,:) :: ug, vg
   
   integer                    :: ipsflg,itsflg,irtsflg,iusflg,io3flg,nsndg
   real, dimension(maxsndg)   :: us,vs,ts,thds,ps,hs,rts,o3s

END MODULE ref_sounding
