!##############################################################################
Subroutine new_base_state (n1,n2,n3)

use mem_grid
use mem_basic
use node_mod
use ref_sounding, only:rt01dn
use micphys, only:idiffperts

implicit none

integer :: i,j,k,n1,n2,n3,nxm2,nym2
double precision :: grwt
real, allocatable :: zttopt(:),basethp(:),basertp(:), th_avg_1d(:), rt_avg_1d(:)
double precision, allocatable :: avgth(:),avgrt(:), dom_avgth(:), dom_avgrt(:)

!Note that this new base state is based on theta-il and NOT
!theta-v as is done for th01dn and th0. This is done so we do not
!have to recompute theta-il (THP) before diffusion.

allocate(zttopt(n1))

if(idiffperts == 1) then

 !Make 3D base state for RTP (rvt0) from 1D base state. TH0 exists for THP.
 !This new base state will be used for numerical diffusion of THP and RTP only.
 do j = 1,n3
  do i = 1,n2

   !Adjusted Grid point height + topography
   do k = 1,n1
    zttopt(k)=zt(k)*grid_g(ngrid)%rtgt(i,j)+grid_g(ngrid)%topt(i,j)
   enddo

   CALL htint (n1,rt01dn(1,ngrid),zt,n1,basic_g(ngrid)%rvt0(1,i,j),zttopt)

  enddo
 enddo

elseif(idiffperts == 2) then

 !Set the weight for each point in domain interior
 nxm2 = nnxp(ngrid) - 2
 nym2 = nnyp(ngrid) - 2
 if(jdim==1)grwt=1.0/float(nxm2*nym2)
 if(jdim==0)grwt=1.0/float(nxm2)

 allocate(avgth(n1))
 allocate(avgrt(n1))
 allocate(dom_avgth(n1))
 allocate(dom_avgrt(n1))
 allocate(th_avg_1d(n1))
 allocate(rt_avg_1d(n1))
 allocate(basethp(n1))
 allocate(basertp(n1))

 !Zero out some quantities
 do k = 1,n1
   avgth(k) = 0.0d0
   avgrt(k) = 0.0d0
   dom_avgth(k) = 0.0d0
   dom_avgrt(k) = 0.0d0
 enddo

 !Loop over subdomain to compute this subdomain's piece of the domain
 !mean thp and rtp while adjusting for topography to get a new "base state".
 do j = ja,jz
  do i = ia,iz

   !Adjusted Grid point height + topography
   do k = 1,n1    
    zttopt(k)=zt(k)*grid_g(ngrid)%rtgt(i,j)+grid_g(ngrid)%topt(i,j) 
   enddo

   !Back interpolate theta-il (THP) and total water (RTP) to a base 
   !state like profile for averaging across a domain that might have 
   !topography
   CALL htint2 (n1,basic_g(ngrid)%thp(1,i,j),zttopt,n1,basethp,zt)
   CALL htint2 (n1,basic_g(ngrid)%rtp(1,i,j),zttopt,n1,basertp,zt)

   !Sum all basestate profiles to get domain mean profile
   !Will use this to get a "domain mean" basestate, rather than using
   !a basestate based on a single domain i,j point. Use double precision
   !for grid averaging for best precision.
   do k = 1,n1
    avgth(k) = avgth(k) + (dble(basethp(k)) * grwt)
    avgrt(k) = avgrt(k) + (dble(basertp(k)) * grwt)
   enddo

  enddo
 enddo

 ! Avgth and avgrt contain the sums that represent this subdomain's piece
 ! of the domaon average. If nmachs == 1, then we already have the
 ! domain average in avgth and avgrt. Otherwise, we need to run an
 ! MPI reduce in order to complete the sum for the whole domain. Running
 ! an MPI allreduce will do the sum and broadcast to all nodes in
 ! one call.
 if (nmachs .eq. 1) then
   do k = 1,n1
     dom_avgth(k) = avgth(k)
     dom_avgrt(k) = avgrt(k)
   enddo
 else
   CALL par_allreduce_sum_doubles (avgth, dom_avgth, n1)
   CALL par_allreduce_sum_doubles (avgrt, dom_avgrt, n1)
 endif

 ! Now dom_avgth and dom_avgrt contain the whole domain average. Copy
 ! these to the TH00 and RT00 variables.
 !
 !Convert averages to single precision values
 do k = 1,n1
   th_avg_1d(k) = sngl(dom_avgth(k))
   rt_avg_1d(k) = sngl(dom_avgrt(k))
 enddo

 !Take the new "domain mean" new 1D base state and make 3D base state.
 !This new base state will be used for numerical diffusion of THP and RTP only.
 do j = 1,n3
  do i = 1,n2

   !Adjusted Grid point height + topography
   do k = 1,n1    
    zttopt(k)=zt(k)*grid_g(ngrid)%rtgt(i,j)+grid_g(ngrid)%topt(i,j) 
   enddo

   CALL htint (n1,th_avg_1d,zt,n1,basic_g(ngrid)%th00(1,i,j),zttopt)
   CALL htint (n1,rt_avg_1d,zt,n1,basic_g(ngrid)%rvt00(1,i,j),zttopt)
 
  enddo
 enddo
 deallocate(basethp,basertp,avgth,avgrt,dom_avgth,dom_avgrt)

endif !If idiffperts==2

deallocate(zttopt)

return
END SUBROUTINE new_base_state

!#######################################################################
Subroutine update_base_winds()

use ref_sounding
use node_mod, only:mmzp
use mem_grid, only:ztn,ngrid,time

implicit none

integer :: it1, it2,nlev
real :: timefac
real, dimension(size(ug,1)) :: newug, newvg

nlev = size(ug,1)
it1 = findloc(time-forc_time >= 0, .TRUE., DIM=1, BACK=.TRUE.) 
if (time - forc_time(it1) == 0.) then
   CALL htint(nlev,ug(:,it1),forc_lev,mmzp(ngrid),u01dn(:,ngrid),ztn(:,ngrid) )
   CALL htint(nlev,vg(:,it1),forc_lev,mmzp(ngrid),v01dn(:,ngrid),ztn(:,ngrid) )
else
   it2 = it1 + 1
   if (it2>size(forc_time)) stop "No future UG/VG for GEOSTROPHIC WIND FORCING"
   timefac = (time - forc_time(it1))/(forc_time(it2)-forc_time(it1))
   newug = ug(:,it1) + (ug(:,it2) - ug(:,it1)) * timefac 
   newvg = vg(:,it1) + (vg(:,it2) - vg(:,it1)) * timefac 

   CALL htint(nlev,newug,forc_lev,mmzp(ngrid),u01dn(:,ngrid),ztn(:,ngrid) )
   CALL htint(nlev,newvg,forc_lev,mmzp(ngrid),v01dn(:,ngrid),ztn(:,ngrid) )
endif

END SUBROUTINE update_base_winds
