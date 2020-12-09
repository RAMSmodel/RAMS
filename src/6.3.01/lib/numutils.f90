!##############################################################################
Subroutine azero (n1,a1)

implicit none

integer :: n,n1
real :: a1(n1)

do n=1,n1
   a1(n)=0.
enddo

return
END SUBROUTINE azero

!##############################################################################
Subroutine azero2 (n1,a1,a2)

implicit none

integer :: n,n1
real :: a1(n1),a2(n1)

do n=1,n1
   a1(n)=0.
   a2(n)=0.
enddo

return
END SUBROUTINE azero2
!##############################################################################
Subroutine azero_int (n1,a1)

implicit none

integer :: n,n1
integer :: a1(n1)

do n=1,n1
   a1(n)=0
enddo

return
END SUBROUTINE azero_int

!##############################################################################
Subroutine ae1 (npts,a,b)

implicit none

integer :: npts
real :: a(npts),b(npts)
integer :: i

do i=1,npts
  a(i)=b(i)
enddo

return
END SUBROUTINE ae1

!##############################################################################
Subroutine ae2 (n2,n3,i1,i2,j1,j2,a,b)

implicit none

integer :: n2,n3,i1,i2,j1,j2
real :: a(n2,n3),b(n2,n3)
integer :: i,j

do j=j1,j2
  do i=i1,i2
    a(i,j)=b(i,j)
  enddo
enddo

return
END SUBROUTINE ae2

!##############################################################################
Subroutine ae0 (npts,a,b)

implicit none

integer :: npts
real :: a(npts),b
integer :: i

do i=1,npts
  a(i)=b
enddo

return
END SUBROUTINE ae0

!##############################################################################
Subroutine update (n,a,fa,dt)

implicit none

integer :: n,nn
real :: a(n),fa(n),dt

do 10 nn=1,n
  a(nn)=a(nn)+fa(nn)*dt
10 continue

return
END SUBROUTINE update

!##############################################################################
Subroutine accum (nxyz,arr1,arr2)

implicit none

integer :: nxyz,n
real :: arr1(nxyz),arr2(nxyz)

do n=1,nxyz
  arr1(n)=arr1(n)+arr2(n)
enddo

return
END SUBROUTINE accum

!##############################################################################
Subroutine atob (n,a,b)

implicit none

integer :: n,i
real :: a(n),b(n)

do 100 i=1,n
b(i)=a(i)
100 continue

return
END SUBROUTINE atob

!##############################################################################
real Function valugp (n1,n2,n3,k,i,j,a)

implicit none

integer :: n1,n2,n3,k,i,j
real :: a(n1,n2,n3)

valugp=a(k,i,j)

return
END FUNCTION valugp

!##############################################################################
Subroutine sort3 (a,b,c,n)

implicit none

integer :: n
real :: a(n),b(n),c(n)
integer :: np1,k,i
integer, external :: ismin
real :: at,bt,ct

np1=n+1
do 10 k=1,n
i=ismin(np1-k,a(k),1)+k-1
at=a(i)
bt=b(i)
ct=c(i)
a(i)=a(k)
b(i)=b(k)
c(i)=c(k)
a(k)=at
b(k)=bt
c(k)=ct
10 continue

return
END SUBROUTINE sort3

!##############################################################################
real Function ssumvect (nn,vctr,inc)

!Return sum of vector.

implicit none

integer :: nn,inc
real :: vctr(*)
integer :: n,nnn
real :: sum

sum=0.
nnn=nn*inc
do 10 n=1,nnn,inc
  sum=sum+vctr(n)
10 continue
ssumvect = sum

return
END FUNCTION ssumvect

!##############################################################################
integer Function ismin (nn,vctr,inc)

!Return index of minimum of vector

implicit none

integer :: nn,inc
real :: vctr(*)
integer :: ism,nnn
real :: smin

ism=0
smin=1e10
do 10 nnn=1,nn,inc
  if(vctr(nnn).lt.smin)then
    ism=nnn
    smin=vctr(nnn)
  endif
10 continue
ismin = ism

return
END FUNCTION ismin

!##############################################################################
real Function cvmgp (vct1,vct2,vct3)

!Return vct1 if vct3 => 0., else vct2.

implicit none

real :: vct1,vct2,vct3

if(vct3.ge.0)then
cvmgp = vct1
else
cvmgp = vct2
endif

return
END FUNCTION cvmgp

!##############################################################################
real Function cvmgm (vct1,vct2,vct3)

!Return vct1 if vct3 <= 0., else vct2.

implicit none

real :: vct1,vct2,vct3

if(vct3.lt.0)then
cvmgm = vct1
else
cvmgm = vct2
endif

return
END FUNCTION cvmgm

!##############################################################################
integer Function check_real (xx,x,nxx,nx)

! Check two corresponding real arrays and see values are close enough

implicit none

integer :: nx,nxx,i
real :: xx(nxx), x(nx), tol

tol = min( (maxval(xx)-minval(xx)),(maxval(x)-minval(x)) )*.0001
do i = 1, min(nxx,nx)
   if (abs(xx(i)-x(i)) > tol) then
      check_real = i
      return
   endif
enddo

check_real = 0

return
END FUNCTION check_real
