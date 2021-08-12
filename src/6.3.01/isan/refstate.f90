!##############################################################################
Subroutine varfile_refstate (n1,n2,n3,thp,pc,pi0,th0,rtp,dn0  &
                 ,dn0u,dn0v,topt,rtgt,zt,ztop,piref,thref,dnref,rtref)

use rconstants
                 
implicit none

integer :: n1,n2,n3,ir,jr,i,j,k,i1,j1
real, dimension(n1,n2,n3) :: thp,pc,pi0,rtp,dn0,th0,dn0u,dn0v
real, dimension(n2,n3)    :: topt,rtgt
real, dimension(n1)       :: zt,piref,thref,rtref,dnref
real :: ztop,topref,c1,c2,c3
real, allocatable :: vctr1(:),vctr2(:)

!                Reference sounding is point with lowest topography

allocate(vctr1(n1),vctr2(n1))
ir=1
jr=1
topref=1.e10
do j=1,n3
   do i=1,n2
      if(topt(i,j).lt.topref) then
         ir=i
         jr=j
         topref=topt(i,j)
      endif
   enddo
enddo

do k=1,n1
   vctr2(k)=zt(k)*(1.-topref/ztop)+topref
enddo

CALL htint2 (n1,thp(1,ir,jr),vctr2,n1,vctr1,zt)
CALL htint2 (n1,rtp(1,ir,jr),vctr2,n1,rtref(1),zt)

do k = 1,n1
   thref(k) = vctr1(k) * (1. + .61 * rtref(k))
enddo
rtref(1) = rtref(2)
thref(1) = thref(2)

piref(1) = pc(1,ir,jr) + g * (vctr2(1) - zt(1))  &
                / (.5 * (thref(1)  &
                + thp(1,ir,jr) * (1. + .61 * rtp(1,ir,jr))))
do k = 2,n1
   piref(k) = piref(k-1) - g * (zt(k)-zt(k-1))  &
             / ( .5 * (thref(k) + thref(k-1)))
enddo

do k = 1,n1
  vctr1(k) = (piref(k) / cp) ** cpor * p00
  dnref(k) = cp * vctr1(k)  &
     / (rgas * thref(k) * piref(k))
enddo

!        Compute 3-D reference state from 1-D reference state

do j=1,n3
  do i=1,n2

    do k=1,n1
      vctr2(k)=zt(k)*rtgt(i,j)+topt(i,j)
    enddo
    CALL htint (n1,piref(1),zt,n1,pi0(1,i,j),vctr2)
    CALL htint (n1,thref(1),zt,n1,th0(1,i,j),vctr2)

    c1=g*2.*(1.-topt(i,j)/ztop)
    c2=(1-cpor)
    c3=cp**c2
    do k=n1-1,1,-1
      pi0(k,i,j)=pi0(k+1,i,j)  &
                +c1*(zt(k+1)-zt(k))/(th0(k,i,j)+th0(k+1,i,j))
    enddo

    do k=1,n1
      dn0(k,i,j)=(c3*p00)/(rgas*th0(k,i,j)*pi0(k,i,j)**c2)
    enddo

  enddo
enddo

do j = 1,n3
   j1 = min(j+1,n3)
   do i = 1,n2
      i1 = min(i+1,n2)
      do k = 1,n1
         dn0u(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i1,j))
         dn0v(k,i,j) = .5 * (dn0(k,i,j) + dn0(k,i,j1))
      enddo
   enddo
enddo

deallocate(vctr1,vctr2)

return
END SUBROUTINE varfile_refstate
