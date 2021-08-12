!##############################################################################
Subroutine mxtked (m1,m2,m3,ia,iz,ja,jz  &
   ,tkep,tket,vt3da,vt3dc,vt3dh,vt3dj,scr1,scr2,dxt,rtgt)
!  +-------------------------------------------------------------------+
!  \  this routine calculates the mixing coefficients based on the     \
!  \    prognostic tke and diagnostic length scale, according to the   \
!  \    deardorff scheme.                                              \
!  +-------------------------------------------------------------------+

use mem_grid
use mem_scratch
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz
real, dimension(m1,m2,m3) :: tkep,tket,vt3da,vt3dc,vt3dh,vt3dj,scr1,scr2
real, dimension(m2,m3) :: dxt,rtgt
integer :: i,j,k
real :: sqrttkep,tket2,c1,sclu,scl

do j = ja,jz
   do i = ia,iz
      tket2 = tket(2,i,j)
      c1 = 3.375 * rtgt(i,j) / dxt(i,j) ** 2  ! 3.375 = 1.5 ** 3
      do k = 2,m1-1
         sqrttkep = sqrt(tkep(k,i,j))
         sclu = (c1 * (zm(k) - zm(k-1))) ** .333333
         if (vt3dj(k,i,j) .lt. 1.e-10) then
            scl = sclu                                   ! unstable case
         else
            scl = .76 * sqrt(tkep(k,i,j) / vt3dj(k,i,j))   !stable case
            if (scl .gt. sclu) scl = sclu
         endif

! now, vt3dh is strain rate and scr2 is dw/dz

         scr1(k,i,j) = .1 * scl * sqrttkep
         vctr25(k) = scr1(k,i,j) * vt3dh(k,i,j)  &
            - .333333 * tkep(k,i,j)  &
            * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
         vt3dh(k,i,j) = scr1(k,i,j) * (1. + 2. * scl / sclu)
         scr2(k,i,j) = scr1(k,i,j)

! now, vt3dh is vkh and scr2 is hkm

         vctr26(k) = .19 + .51 * scl / sclu
         vctr27(k) = - vt3dh(k,i,j) * vt3dj(k,i,j)
         vctr28(k) = vctr26(k) * tkep(k,i,j) * sqrttkep / scl
         tket(k,i,j) = tket(k,i,j) + vctr25(k) + vctr27(k)  &
            - vctr28(k)
      enddo
      tket(2,i,j) = tket2 + vctr25(3) + vctr27(2) - vctr28(2)

   enddo
enddo

return
END SUBROUTINE mxtked

!##############################################################################
Subroutine tkemy (m1,m2,m3,ia,iz,ja,jz,jd  &
   ,tkep,tket,vt3dh,vt3di,vt3dj,scr1,rtgt  &
   ,theta,dn0,up,vp,wp,sflux_u,sflux_v,sflux_w,sflux_t,tkep2)

use mem_grid
use mem_scratch
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd
real, dimension(m1,m2,m3) :: tkep,tket,vt3dh,vt3di,vt3dj,scr1,theta,dn0  &
   ,up,vp,wp
real, dimension(m2,m3) :: rtgt,sflux_u,sflux_v,sflux_w,sflux_t
real, dimension(m1) :: tkep2
integer :: i,j,k
real :: a1,a2,b1,b2,c1,aux1,aux2,rf1,rf2,rf3,rf4,wght1,wght3,sumtkz,sumtk  &
   ,al0,tket2,ri,rf,shr,smr,tker,qq,ssmf,shf,sh0,ssm,aux,gm,gh,sm1,sm2  &
   ,sh1,sh2,dzloc
real, external :: ssumvect

data a1,a2,b1,b2,c1/0.92,0.74,16.6,10.1,0.08/
data aux1,aux2/0.758964199,2.58286747/
data rf1,rf2,rf3,rf4/1.,0.191232309,0.223117196,0.234067819/

!        7 - mellor and yamada (after andre et al, 1978)

!Saleeby(2014):Check older code for possible experimental values

wght3=1.
wght1=1.0-wght3

do j=ja,jz
   do i=ia,iz
     do k=2,m1-1
         tkep2(k) = 2.0 * tkep(k,i,j)
         vctr30(k) = sqrt(tkep2(k))
!Upper limit for length scale "L" in stable conditions (Andre et al 1978)
         vctr31(k)=.75*vctr30(k)/sqrt(max(1.e-20,vt3dj(k,i,j)))
         vctr1(k)=(zt(k)-zm(1))*rtgt(i,j)
         dzloc=(zm(k)-zm(k-1))*rtgt(i,j)
         vctr33(k)=vctr30(k)*dzloc
         vctr32(k)=vctr33(k)*vctr1(k)
      enddo

      sumtkz = ssumvect(m1-2,vctr32(2),1)
      sumtk  = ssumvect(m1-2,vctr33(2),1)
      al0=.1*sumtkz/sumtk
      tket2=tket(2,i,j)

      do k=2,m1-1
!Turbulent length scale "L" after Mellor-Yamada (1982)
         vctr9(k)=min(vonk*vctr1(k)/(1.+vonk*vctr1(k)/al0)  &
              ,vctr31(k))

! --- for growing turbulence use helfand and labraga's modified sm and sh (sh0)
!Brunt-Vaisala contribution over vertical strain
         ri=min(vt3dj(k,i,j)/max(vt3di(k,i,j),1.e-11),0.190)
         rf=min(0.6588*(ri+0.1776-sqrt(ri*(ri-0.3221)  &
              +0.03156)),0.16)
         shr=aux2*(rf-rf2)/(rf-rf1)
         smr=aux1*(rf-rf4)/(rf-rf3)*shr
         tker=max(16.6*vctr9(k)*vctr9(k)*(smr*vt3di(k,i,j)  &
            -shr*vt3dj(k,i,j)),2.*tkmin)
         if(tker.gt.tkep2(k) )then
            qq=sqrt(tkep2(k)/tker)
            ssmf=qq*smr
            shf=qq*shr
            sh0=shf
            ssm=ssmf
         else

! --- for decaying turbulence use mellor and yamada's sm and sh (sh0)
            aux=vctr9(k)*vctr9(k)/tkep2(k) ! L^2 / 2e
            gm=aux*vt3di(k,i,j) ! L^2 / 2e  * vertical strain
            gh=-aux*vt3dj(k,i,j) ! L^2 / 2e  * Brunt-vaisala contribution

            sm1=0.6992-9.33948672*gh
            sm2=1.-(36.7188-187.4408515*gh+88.83949824*gm)*gh  &
                 +5.0784*gm
            ssm=sm1/max(sm2,1e-10) !Sm nondimensional eddy diffusivity

            sh1=0.74-4.0848*ssm*gm
            sh2=1.-30.5916*gh
            sh0=sh1/max(sh2,1e-10) !Sh nondimensional eddy diffusivity
         endif

!Vertical eddy diffusivity for momentum VKM
         scr1(k,i,j)=vctr9(k)*vctr30(k)*ssm
!Vertical eddy diffusivity for heat (scalars) VKH
         vt3dh(k,i,j)=vctr9(k)*vctr30(k)*sh0
!Ps + Pb or (shear production term + buoyancy production term in
!prognostic TKE equation
!(VKM*[(du/dz)^2 + (dv/dz)^2])  - (VKH*(-g/theta * dtheta/dz))
         vctr5(k)=scr1(k,i,j)*vt3di(k,i,j)  &
            -vt3dh(k,i,j)*vt3dj(k,i,j)
         tket(k,i,j)=tket(k,i,j) + 0.5 * vctr5(k)  &
            -tkep2(k) * sqrt(max(1.e-20,tkep2(k)))/(vctr9(k)*16.6)
!Density weighting for VKM and VKH
         scr1(k,i,j)=scr1(k,i,j)*dn0(k,i,j)
         vt3dh(k,i,j)=vt3dh(k,i,j)*dn0(k,i,j)
!Vertical eddy diffusivity for TKE
         vt3di(k,i,j)=0.2*vctr9(k)*vctr30(k)*dn0(k,i,j)
!print*,'VHM= ',scr1(k,i,j),'  VKH= ',vt3dh(k,i,j)
      enddo
      if(nstbot.eq.1)then

         tket(2,i,j) =tket2  &
            + 0.5 * wght3 * vctr5(2+1)   &
            + wght1 *  &
             (-(sflux_u(i,j)*(up(2,i,j)+up(2,i-1,j))  &
               +sflux_v(i,j)*(vp(2,i,j)+vp(2,i,j-jd))  &
               +sflux_w(i,j)*wp(2,i,j))/vctr1(2)  &
               +sflux_t(i,j)*g/theta(2,i,j))/dn0(2,i,j) &
            
            - tkep2(2)*sqrt(max(1.e-20,tkep2(2)))/(vctr9(2)*16.6)

      endif
      
   enddo
enddo

return
END SUBROUTINE tkemy

!##############################################################################
Subroutine tkeinit (n1,n2,n3)

use mem_grid
use mem_turb
use rconstants

implicit none

integer :: n1,n2,n3
integer :: i,j,k

!        limit the values to the minimum

if( allocated(turb_g(ngrid)%tkep) ) then
   do j = 1,n3
      do i = 1,n2
         do k = 1,n1
            turb_g(ngrid)%tkep(k,i,j) = max(tkmin,turb_g(ngrid)%tkep(k,i,j))
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE tkeinit
