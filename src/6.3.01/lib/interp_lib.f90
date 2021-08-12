!##############################################################################
Subroutine binom (x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)

implicit none

real :: x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy
real :: wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
integer :: istend

yyy=1e30;yz22=0.;yz23=0. !Variable initialized

if(x2.gt.1.e19.or.x3.gt.1.e19.or.  &
   y2.gt.1.e19.or.y3.gt.1.e19)return
wt1=(xxx-x3)/(x2-x3)
wt2=1.0-wt1
istend=0
if(y4.lt.1.e19.and.x4.lt.1.e19) go to 410
yz22=wt1
yz23=wt2
yz24=0.0
istend= 1
410   if(y1.lt.1.e19.and.x1.lt.1.e19) go to 430
yz11=0.0
yz12=wt1
yz13=wt2
if(istend.eq.1)go to 480
go to 450
430   yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
if(istend.eq.  1    ) go to 470
450   yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
470   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
 go to 490
480      yyy=wt1*y2+wt2*y3
490   yoo=yyy

return
END SUBROUTINE binom

!##############################################################################
Subroutine gdtost (a,ix,iy,stax,stay,staval)

implicit none

integer :: ix,iy
real :: a(ix,iy),r(4),scr(4),stax,stay,staval

!     ROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
!     FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
!     GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
!     IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
!     LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
!     AND STATION COLUMN.
!     VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.

integer :: iy1,iy2,ix1,ix2,ii,i,jj,j
real :: fiym2,fixm2,yy,xx

iy1=int(stay)-1
iy2=iy1+3
ix1=int(stax)-1
ix2=ix1+3
staval=1e30
fiym2=float(iy1)-1
fixm2=float(ix1)-1
ii=0
do 100 i=ix1,ix2
ii=ii+1
if(i.ge.1.and.i.le.ix) go to 101
scr(ii)=1e30
go to 100
101   jj=0
do 111 j=iy1,iy2
jj=jj+1
if(j.ge.1.and.j.le.iy) go to 112
r(jj)=1e30
go to 111
112   r(jj)=a(i,j)
111   continue
yy=stay-fiym2
CALL binom (1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
100   continue
xx=stax-fixm2
CALL binom (1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)

return
END SUBROUTINE gdtost

!##############################################################################
Subroutine gdtost2 (a,ix,iy,stax,stay,staval)

!     Routine to return stations back-interpolated values (staval)
!     from uniform grid points using bi-linear interpolation.
!     Gridded values of input array a dimensioned a(ix,iy), where
!     ix=grid points in x, iy = grid points in y.  Station
!     location given in terms of grid relative station x,y (stax,stay).

implicit none

! passed variables

integer :: ix, iy
real :: stax,stay,staval
real, dimension(ix,iy) :: a

! internal variables

integer :: i,j
real :: wtx1,wtx2,wty1,wty2
staval = 1.e30
i = int(stax)
j = int(stay)

if (iy > 1) then !2D case
   if(i < 1) i=1
   if(i > ix-1) i=ix-1
   if(j < 1) j=1
   if(j > iy-1) j=iy-1

   if(a(i,j)   > 1.e19 .or. a(i,j+1)   > 1.e19 .or. &
      a(i+1,j) > 1.e19 .or. a(i+1,j+1) > 1.e19 ) return

   wtx2 = stax - float(i)
   wty2 = stay - float(j)
   wtx1 = 1. - wtx2
   wty1 = 1. - wty2

   staval = wtx1 * (wty1 * a(i  ,j  )+ wty2 * a(i  ,j+1)) &
          + wtx2 * (wty1 * a(i+1,j  )+ wty2 * a(i+1,j+1))
else !1D case
   if(i < 1) i=1
   if(i > ix-1) i=ix-1

   if(a(i,j)   > 1.e19 .or. a(i+1,j) > 1.e19) return

   wtx2 = stax - float(i)
   wtx1 = 1. - wtx2

   staval = wtx1 * a(i,j) + wtx2 * a(i+1,j)
endif

return
END SUBROUTINE gdtost2

!##############################################################################
Subroutine htint (nzz1,vctra,eleva,nzz2,vctrb,elevb)

implicit none

integer :: nzz1,nzz2
real :: vctra(nzz1),vctrb(nzz2),eleva(nzz1),elevb(nzz2)
integer :: l,k,kk
real :: wt

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 35
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
if(elevb(k).gt.eleva(nzz1))go to 36
l=l+1
if(l.eq.nzz1) then
  print *,'htint:nzz1',nzz1
  do kk=1,l
    print*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
  enddo
  stop 'htint'
endif
go to 30
35 continue
if(elevb(k).eq.eleva(l+1) .and. (l+1).lt.nzz1) l=l+1
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
!Override above arithmetic for numerical precision when elevations match
if(elevb(k).eq.eleva(l+1)) vctrb(k)=vctra(l+1)
go to 20
36 continue
wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
20 continue

return
END SUBROUTINE htint

!##############################################################################
Subroutine htint2 (nzz1,vctra,eleva,nzz2,vctrb,elevb)

implicit none

integer :: nzz1,nzz2
real :: vctra(nzz1),vctrb(nzz2),eleva(nzz1),elevb(nzz2)
integer :: l,k
real :: wt

!      htint for holding values of vctrb constant under eleva(1)

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 34
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
if(elevb(k).gt.eleva(nzz1))go to 36
l=l+1
if(l.eq.nzz1)stop 'htint2'
go to 30
34   continue
vctrb(k)=vctra(1)
go to 20
35 continue
if(elevb(k).eq.eleva(l+1) .and. (l+1).lt.nzz1) l=l+1
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
!Override above arithmetic for numerical precision when elevations match
if(elevb(k).eq.eleva(l+1)) vctrb(k)=vctra(l+1)
go to 20
36 continue
wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
20 continue

return
END SUBROUTINE htint2
