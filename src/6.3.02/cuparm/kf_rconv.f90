!##############################################################################
! KAIN-FRITSCH CUMULUS PARAMETERIZATION SCHEME 
! Updated July 2002 
!
! RAMS-KF interface originally written, modified, and tested by:
!
! Christopher L. Castro*
! William Y.Y. Cheng
! Adriana B. Beltran
! 
! Department of Atmospheric Science 
! Colorado State University
! Fort Collins, CO  80523
!
! *E-mail contact: chris@atmos.colostate.edu
!
! References:
!
! Castro, C.L., W.Y.Y. Cheng, A.B. Beltran, R.A. Pielke, Sr., and
! W.R. Cotton, 2002. The Incorporation of the Kain-Fritsch Cumulus
! Parameterization Scheme in RAMS.  In preparation.
!
! Kain, J.S., and J.M. Fritsch, 1993. Convective Parameterization
! for Mesoscale Models, The Kain-Fritsch Scheme.  The Representation
! of Convection in Numerical Models. Meteor. Monogr., No. 24,
! Amer. Meteor. Soc., 165-170.
!
!  
! Note: For a thorough description of variables used in the KF interface
! code, please see kfdriver.f90

!##############################################################################
Subroutine kf_main ()

use mem_tend
use mem_cuparm
use mem_micro
use mem_basic
use mem_grid
use mem_tend
use mem_turb
use node_mod
use micphys

implicit none

integer :: j_lu,ng
integer, dimension(mxp,myp) :: nca_int,convgo_int

ng=ngrid

 !Zero out arrays at model start unless doing a history-init
 if(nint(time) == 0 .and. initial /= 3) then
   if(iprntstmt>=1 .and. print_msg) &
      print*,'Zeroing out arrays for K-F cuparm'

   CALL azero (mxp*myp*mzp,cuparm_g(ng)%thsrc(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%rtsrc(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%rcsrc(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%rrsrc(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%rssrc(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%rpsrc(1,1,1))
 
   CALL azero (mxp*myp,cuparm_g(ng)%conprr(1,1))
 
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%w0avg(1,1,1))
   CALL azero (mxp*myp*mzp,cuparm_g(ng)%w0avglt(1,1,1))
   CALL azero (mxp*myp,cuparm_g(ng)%nca(1,1))
   CALL azero (mxp*myp,cuparm_g(ng)%convgo(1,1))
 
endif

CALL var_change (mxp*myp,cuparm_g(ng)%nca(1,1),nca_int(1,1),1)
CALL var_change (mxp*myp,cuparm_g(ng)%convgo(1,1),convgo_int(1,1),1)

! Check for convective activity at the specified frequency in
! RAMSIN.  A value of 600 seconds for CONFRQ is recommended

 if (mod(time + .001,confrq) .lt. dtlt .or. time .lt. 0.001) then

   if(iprntstmt>=1 .and. print_msg) &
      print*,'Check For Convection in K-F at ', time, ' sec'

! Call the KF driver (see kfdriver.f90)

   j_lu=1
   CALL kfdrive (mzp,mxp,myp,ia,iz,ja,jz,level,icloud,irain  &
       ,ipris,isnow,idiffk(ng)        &             
       ,basic_g(ng)%up(1,1,1),   basic_g(ng)%vp(1,1,1)  &
       ,basic_g(ng)%theta(1,1,1),basic_g(ng)%thp(1,1,1)  &
       ,basic_g(ng)%pp(1,1,1),   basic_g(ng)%pi0(1,1,1)  &
       ,tend%pt(1)  &
       ,basic_g(ng)%rv(1,1,1)    &
       ,grid_g(ng)%rtgt(1,1),  grid_g(ng)%topt(1,1),dzt,zt   &
       ,grid_g(ng)%dxt(1,1)   &
       ,nca_int(1,1),          convgo_int(1,1)  &
       ,cuparm_g(ng)%w0avg(1,1,1),cuparm_g(ng)%w0avglt(1,1,1) &
       ,turb_g(ng)%tkep(1,1,1),dtlt         &               
       ,cuparm_g(ng)%thsrc(1,1,1),cuparm_g(ng)%rtsrc(1,1,1)  &
       ,cuparm_g(ng)%rcsrc(1,1,1),cuparm_g(ng)%rrsrc(1,1,1) &
       ,cuparm_g(ng)%rssrc(1,1,1),cuparm_g(ng)%rpsrc(1,1,1)  &
       ,cuparm_g(ng)%conprr(1,1),j_lu)

endif

! Update the vertical velocity parameters needed by the KF scheme

CALL update_w0avgkf (mzp,mxp,myp,ia,iz,ja,jz   &
       ,grid_g(ng)%f13t(1,1),grid_g(ng)%f23t(1,1)  &
       ,grid_g(ng)%rtgt(1,1),grid_g(ng)%topt(1,1),ztop,zt  &
       ,basic_g(ng)%up(1,1,1),basic_g(ng)%vp(1,1,1)  &
       ,basic_g(ng)%wp(1,1,1)  &
       ,cuparm_g(ng)%w0avg(1,1,1),cuparm_g(ng)%w0avglt(1,1,1) &
       ,dtlt)
! Update convective tendencies

! Theta_il tendency

CALL accum_kf (mzp,mxp,myp,tend%tht(1)              &
            ,cuparm_g(ng)%thsrc(1,1,1),nca_int(1,1))  

! Total water mixing ratio tendency

CALL accum_kf (mzp,mxp,myp,tend%rtt(1)              &
            ,cuparm_g(ng)%rtsrc(1,1,1),nca_int(1,1))  

! Cloud water mixing ratio tendency

if(level == 3) then

! Cloud water mixing ratio tendency
  if(jnmb(1) >= 1) &
   CALL accum_kf (mzp,mxp,myp,tend%rct(1)              &
            ,cuparm_g(ng)%rcsrc(1,1,1),nca_int(1,1))  

! Rain water mixing ratio tendency
  if(jnmb(2) >= 1) &
   CALL accum_kf (mzp,mxp,myp,tend%rrt(1)              &
            ,cuparm_g(ng)%rrsrc(1,1,1),nca_int(1,1))  

! Snow mixing ratio tendency
  if(jnmb(4) >= 1) &
   CALL accum_kf (mzp,mxp,myp,tend%rst(1)              &
            ,cuparm_g(ng)%rssrc(1,1,1),nca_int(1,1))  

! Pristine ice ratio tendency
  if(jnmb(3) >= 1) &
   CALL accum_kf (mzp,mxp,myp,tend%rpt(1)              &
            ,cuparm_g(ng)%rpsrc(1,1,1),nca_int(1,1))  

endif

! Convective precipitation (mm or kg m^-2)

CALL update_kf (mxp*myp,cuparm_g(ng)%aconpr(1,1),      &
             cuparm_g(ng)%conprr(1,1),nca_int(1,1),dtlt)           
                                           
! Updates the NCA array

CALL nca_update (mxp*myp,cuparm_g(ng)%nca(1,1),nca_int(1,1))                
CALL var_change (mxp*myp,cuparm_g(ng)%convgo(1,1),convgo_int(1,1),2)
                                               
return
END SUBROUTINE kf_main

!##############################################################################
Subroutine var_change (npoints,var,var_int,vartype)

implicit none

integer :: npoints,j
real :: var(npoints)
integer :: var_int(npoints),vartype

      do j=1,npoints
       if(vartype==1) then
        var_int(j)=int(var(j))
       else if(vartype==2) then
        var(j)=float(var_int(j))
       endif
      enddo

return
END SUBROUTINE var_change

!##############################################################################
!  Routine to compute the contravariant component of vertical
!  velocity. Two new vertical velocity components created here:
!
!     w0avg: running mean vertical velocity (for t=600s)
!     w0avglt: running mean horizontal components of contravariant vertical
!              velocity, scaled by elevation
!
!  The scheme was originally tuned for a July 2001 simulation of the 
!  North American Monsoon at 33 km grid spacing.  Users may wish to modify
!  the variable topocoef in if precipitation amounts 
!  in areas of significant terrain are physically unreasonable.
!
!  RAMS USERS SHOULD INCLUDE THIS OR A SIMILAR ROUTINE WHEN INCORPORATING
!  ANY ADDITIONAL CONVECTION SCHEMES IN RAMS GIVEN THE PRESENT VERTICAL 
!  COORDINATE!
!------------------------------------------------------------------------------
Subroutine update_w0avgkf (mzp,mxp,myp,ia,iz,ja,jz,   &
  ter_gradx,ter_grady,rtgt,topot,ztop,zt,up,vp,wp,w0avg,w0avglt,dt)

implicit none

integer ntst,k,i,j,ja,jz,ia,iz,mzp,mxp,myp
real dt,ztop
real zt(mzp)  
real up(mzp,mxp,myp), vp(mzp,mxp,myp)
real wp(mzp,mxp,myp), w0avglt(mzp,mxp,myp), w0avg(mzp,mxp,myp)
real topot(mxp,myp)
real ter_gradx(mxp,myp), ter_grady(mxp,myp), rtgt(mxp,myp)
real ucond, vcond
real w0_parm, w0_parmlt

      ntst=nint(600./dt)         ! number of time steps in 600 sec


! The horizontal components of the contravariant vertical velocity.  
! Note that here (in contrast to same routine in Kuo) the terrain related 
! components must be cubed because the trigger func in KF is proportional 
! to w^(1/3)


      do k=2,mzp-1
       do j=ja,jz
        do i=ia,iz

         ucond=up(k,i,j)*((((1.-(zt(k)-topot(i,j))/(ztop))*ter_gradx(i,j))  &
              /(rtgt(i,j))))**3
         vcond=vp(k,i,j)*((((1.-(zt(k)-topot(i,j))/(ztop))*ter_grady(i,j))  &
              /(rtgt(i,j))))**3  
       
         w0_parmlt=(ucond+vcond)
	 
         w0_parm=wp(k,i,j)
         

! Running mean vertical velocity

         w0avg(k,i,j)=(w0avg(k,i,j)*(real(ntst)-1.)+w0_parm)/real(ntst)       
         if(w0avg(k,i,j)>0.0 .and. w0avg(k,i,j)<1.e-6)w0avg(k,i,j)=0.0
         if(w0avg(k,i,j)<0.0 .and. w0avg(k,i,j)>-1.e-6)w0avg(k,i,j)=0.0

! Running mean horizontal components of contravariant velocity

         w0avglt(k,i,j)=(w0avglt(k,i,j)*(real(ntst)-1.)+w0_parmlt)/real(ntst)   
         if(w0avglt(k,i,j)>0.0 .and. w0avglt(k,i,j)<1.e-6)w0avglt(k,i,j)=0.0
         if(w0avglt(k,i,j)<0.0 .and. w0avglt(k,i,j)>-1.e-6)w0avglt(k,i,j)=0.0

        enddo
       enddo
      enddo

return
END SUBROUTINE update_w0avgkf

!##############################################################################
Subroutine nca_update (npoints,nca,nca_int)

implicit none

integer :: j,npoints
integer :: nca_int(npoints)
real :: nca(npoints)

      do j=1,npoints
       nca(j)=float(nca_int(j)-1)
      enddo

return
END SUBROUTINE nca_update

!##############################################################################
Subroutine accum_kf (mzp,mxp,myp,tht,thsrc,nca)

!routine for adding source term to tendency

implicit none

integer :: k,i,j,mzp,mxp,myp
real :: tht(mzp,mxp,myp),thsrc(mzp,mxp,myp)
integer :: nca(mxp,myp)

      do j=1,myp
       do i=1,mxp
       if (nca(i,j).gt.0) then
        do k=2,mzp
         tht(k,i,j)=tht(k,i,j)+thsrc(k,i,j)
        enddo
       endif
       enddo
      enddo

return
END SUBROUTINE accum_kf

!##############################################################################
Subroutine update_kf (npoints,acccon,conpcp,nca,dt)

!routine for updating accumulated precip

implicit none

integer :: j,npoints
real :: dt
real :: acccon(npoints),conpcp(npoints)
integer :: nca(npoints)

      do j=1,npoints
       if (nca(j).gt.0) then
        if(conpcp(j).gt.0.0) then
          acccon(j)=acccon(j)+conpcp(j)*dt
        endif
        if(acccon(j).lt.0.0) then
          acccon(j)=0.0
        endif
       endif
      enddo

return
END SUBROUTINE update_kf

