!##############################################################################
Subroutine oda_nudge ()

use mem_oda
use mem_grid
use mem_scratch
use mem_tend
use mem_basic
use io_params
use node_mod

implicit none

integer :: ncall=0,nobs,ng

if(iprntstmt>=1)print*,'start nud'

! Namelist variables:

! - Flag to turn on oda                     (IF_ODA)
! - Model start and end time for oda        (TODABEG,TODAEND)
! - Frequency of tendency updates           (FRQODA) (0. for every coarse grid dt)
! - For each grid:
!     - Nudging timescale                   (TNUDODA)
!     - Number of vertical levels to apply  (no)
!     - Flag to turn on for this grid       (IG_ODA)


! Covariance stuff:

! For each grid:
! - Surface radii - 1) e-2  2) zero         (RODA_SFCE, RODA_SFC0) (SFC0 <= UPA0)
! - "mid-level" radii - same                (RODA_UPAE, RODA_UPA0)
! - "mid-level" radii height                (RODA_HGT)
!         (linear change from surface to mid, constant above)
!
! - vertical "100" factor                   (RODA_ZFAC)

! Obs processing:
! - Time interpolate limit (TIL)- if the future-past obs time 
!    is > this limit, do not use to interpolate (ODA_SFC_TIL,ODA_UPA_TIL)
!
! - Time extrapolate limit (TEL)- if past/future obs is greater than TIL,
!    but less than TEL, use obs                 (ODA_SFC_TEL,ODA_UPA_TEL)



if (ncall == 0) then
   CALL oda_nudge_init (ngrids,mmzp(1),mmxp(1),mmyp(1))
   ncall = 1
endif

if(iprntstmt>=1)print*,'time,frqoda:',my_rams_num,time,frqoda,mod(time,frqoda),dtlongn(1)

if (frqoda > 0. .and. mod(time,frqoda) >= dtlongn(1) ) go to 20


if ( (ngrid == 1 .and. time >= todabeg .and. time <= todaend) ) then
      
   ! Compute new krigged fields and variances

   do ng=1,ngrids
   
      if (wt_oda_grid(ng) > 0.0) then
   
      ! Process observations and call krigging 

      CALL oda_proc_obs (mmzp(ng),mmxp(ng),mmyp(ng),mi0(ng),mj0(ng),ng,nobs)
            
      if(iprntstmt>=1)print*,'oda-sfc nobs:',ng,nobs

      !do i=1,nobs
      !print*,'uobs:',i,ukobs(i),oda_sfc_obs(i)%u(1:oda_sfc_info(i)%ntimes)
      !enddo
      oda_g(ng)%uk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%ukv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      CALL krig (mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,ukobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt(1,1)  &
               ,oda_g(ng)%uk(1,1,1),oda_g(ng)%ukv(1,1,1),1.)

      oda_g(ng)%vk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%vkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      CALL krig (mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,vkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt(1,1)  &
               ,oda_g(ng)%vk(1,1,1),oda_g(ng)%vkv(1,1,1),1.)

      !do i=1,nobs
      !print*,my_rams_num,'---tobs:',ng,i,tkobs(i)
      !enddo
      oda_g(ng)%tk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%tkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      CALL krig (mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,tkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt(1,1)  &
               ,oda_g(ng)%tk(1,1,1),oda_g(ng)%tkv(1,1,1),1.)
      
      oda_g(ng)%rk(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      oda_g(ng)%rkv(1:mmzp(ng),1:mmxp(ng),1:mmyp(ng))=0.
      CALL krig (mmzp(ng),mmxp(ng),mmyp(ng)  &
               ,xtn(1+mi0(ng),ng),ytn(1+mj0(ng),ng),ztn(1,ng)  &
               ,nobs,xkobs,ykobs,zkobs,ekobs,rkobs  &
               ,ng,nnzp(ng),grid_g(ng)%topt(1,1)  &
               ,oda_g(ng)%rk(1,1,1),oda_g(ng)%rkv(1,1,1),1.)
      endif
      
   enddo

endif

20 continue

! Compute and apply tendencies

ng=ngrid
if (wt_oda_grid(ng) > 0.0 .and. time >= todabeg .and. time <= todaend) then

CALL oda_tendency (mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%up(1,1,1)  &
                 ,tend%ut(1),oda_g(ng)%uk(1,1,1),oda_g(ng)%ukv(1,1,1)  &
                 ,wt_oda_uv * wt_oda_grid(ng)/tnudoda)
CALL oda_tendency (mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%vp(1,1,1)  &
                 ,tend%vt(1),oda_g(ng)%vk(1,1,1),oda_g(ng)%vkv(1,1,1)  &
                 ,wt_oda_uv * wt_oda_grid(ng)/tnudoda)
CALL oda_tendency (mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%theta(1,1,1)  &
                 ,tend%tht(1),oda_g(ng)%tk(1,1,1),oda_g(ng)%tkv(1,1,1)  &
                 ,wt_oda_th * wt_oda_grid(ng)/tnudoda)
CALL oda_tendency (mmzp(ng),mmxp(ng),mmyp(ng),mia(ng),miz(ng),mja(ng),mjz(ng)  &
                 ,basic_g(ng)%rtp(1,1,1)  &
                 ,tend%rtt(1),oda_g(ng)%rk(1,1,1),oda_g(ng)%rkv(1,1,1)  &
                 ,wt_oda_rt * wt_oda_grid(ng)/tnudoda)

endif

return
END SUBROUTINE oda_nudge

!##############################################################################
Subroutine oda_nudge_init (ngg,m1m,m2m,m3m)

use mem_oda
use mem_grid

implicit none

integer :: ngg
integer, dimension(ngg) :: m1m,m2m,m3m

integer ::ng,k

! Turn namelist parameters into krigging routine parameters
if(iprntstmt>=1)print*,'nud_init:'

do ng=1,ngrids
   if (wt_oda_grid(ng) > 0.0) then
      nstkrg(ng)=1
      ckrg(1:nstkrg(ng),ng) = 1
      do k=1,nnzp(ng)
         if(ztn(k,ng) < roda_hgt(ng)) then
            rmaxkrg(k,ng)=   &
               roda_sfc0(ng)+(roda_upa0(ng)-roda_sfc0(ng))  &
                             * ztn(k,ng)/roda_hgt(ng)
            ! Kriging routine needs the following to be negative...
            akrg(k,ng)   = -  &
               (roda_sfce(ng)+(roda_upae(ng)-roda_sfce(ng))  &
                             * ztn(k,ng)/roda_hgt(ng)) 
         else
            rmaxkrg(k,ng)=roda_upa0(ng)
            akrg(k,ng) = -roda_upae(ng)
         endif
         !print*,k,ztn(k,ng), roda_hgt(ng),rmaxkrg(k,ng),akrg(k,ng)
      enddo
      
      caxkrg(1,ng)= 1. ; caykrg(1,ng)= 0. ; cazkrg(1,ng)= 0.
      caxkrg(2,ng)= 0. ; caykrg(2,ng)= 1. ; cazkrg(2,ng)= 0.
      caxkrg(3,ng)= 0. ; caykrg(3,ng)= 0. ; cazkrg(3,ng)= roda_zfac(ng)
      
      ! Initialize analysis and variance fields so zero tendency will be computed

      oda_g(ng)%uk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%vk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%tk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%rk(1:m1m(ng),1:m2m(ng),1:m3m(ng))=0.
      oda_g(ng)%ukv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%vkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%tkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      oda_g(ng)%rkv(1:m1m(ng),1:m2m(ng),1:m3m(ng))=2.
      
   endif
enddo

return
END SUBROUTINE oda_nudge_init

!##############################################################################
Subroutine oda_tendency (n1,n2,n3,ia,iz,ja,jz,ap,at,ak,akv,tscalei)

use mem_oda
use mem_grid, only:iprntstmt

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz
real, dimension(n1,n2,n3) :: ap,at,ak,akv
real :: tscalei
real :: ttmin,ttmax

integer :: i,j,k
real :: fnna

! Compute and apply tendencies based on krigged field/variance 

if(iprntstmt>=1)print*,'oda_tend:',ia,iz,ja,jz,n1,n2,n3

ttmax=-1e20; ttmin=1e20
do j=ja,jz
   do i=ia,iz
      do k=1,n1
         fnna=(1.0-akv(k,i,j)/2.0)*(ak(k,i,j)-ap(k,i,j))*tscalei
         at(k,i,j)=at(k,i,j)+fnna
      enddo
   enddo
enddo

return
END SUBROUTINE oda_tendency
