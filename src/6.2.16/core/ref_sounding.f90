!##############################################################################
Module ref_sounding

use grid_dims

implicit none

   integer :: iref,jref
   real :: topref

   !Base state variables in 1D (constant in time)
   real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn
   
   integer                    :: ipsflg,itsflg,irtsflg,iusflg,nsndg
   real, dimension(maxsndg)   :: us,vs,ts,thds,ps,hs,rts

END MODULE ref_sounding
