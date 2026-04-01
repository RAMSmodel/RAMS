!##############################################################################
Subroutine trsets_ns ()

use mem_grid
use mem_micro
use mem_radiate
use micphys
use node_mod, only:mzp,mxp,myp

implicit none

!Here we set top and bottom boundaries for non-advected
!3D vars or "non-scalars" since our definition of scalar refers
!to 3D vars that are advected and diffused.

if(ilwrtyp+iswrtyp > 0) then
 CALL trsetns (mzp,mxp,myp,radiate_g(ngrid)%fthrd,1)
endif

if(ilwrtyp == 3 .or. iswrtyp == 3) then
 CALL trsetns (mzp,mxp,myp,radiate_g(ngrid)%bext,2)
endif

if(jnmb(2) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvr,2)
endif

if(jnmb(3) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvp,2)
endif

if(jnmb(4) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvs,2)
endif

if(jnmb(5) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpva,2)
endif

if(jnmb(6) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvg,2)
endif

if(jnmb(7) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvh,2)
endif

if(jnmb(8) >= 1) then
 CALL trsetns (mzp,mxp,myp,micro_g(ngrid)%pcpvd,2)
endif

return
END SUBROUTINE trsets_ns

!##############################################################################
Subroutine trsetns (m1,m2,m3,ap,flg)

use mem_grid

implicit none

real, dimension(m1,m2,m3) :: ap
integer :: m1,m2,m3,i,j,flg
real :: dzmr

dzmr = dzm(m1-2) / dzm(m1-1)

!Top and bottom boundary setting
do j = 1,m3
 do i = 1,m2

  if(flg==1) &
   ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))

  if(flg==2) &
   ap(m1,i,j) = max(0., ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j)))

  ap(1,i,j)  = ap(2,i,j)

 enddo
enddo

return
END SUBROUTINE trsetns
