!##############################################################################
Subroutine rams_comp_empty (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=0.
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_empty

!##############################################################################
Subroutine rams_comp_dir (n1,n2,n3,a,b,ngrd)

use mem_grid, only:polelat,polelon,xtn,ytn

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: a,b
real :: qlat,qlon,u,v,ff

do k=1,n3
   do j=1,n2
      do i=1,n1
         CALL xy_ll (qlat,qlon,polelat,polelon  &
            ,xtn(i,ngrd),ytn(j,ngrd))
         u=a(i,j,k)
         v=b(i,j,k)
         CALL uvtoueve (u,v,a(i,j,k),b(i,j,k)  &
                      ,qlat,qlon,polelat,polelon)
         CALL winddf (a(i,j,k),ff,a(i,j,k),b(i,j,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_dir

!##############################################################################
Subroutine rams_comp_dewK (n1,n2,n3,a,b,c)

use rconstants, only:cp,cpor,p00

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c
real :: xpress,xtemp,xwatsat
real, external :: tdewpt,rsatmix

do k=1,n3
   do j=1,n2
      do i=1,n1
         xpress=(b(i,j,k)/cp)**cpor*p00
         xtemp=c(i,j,k)*b(i,j,k)/cp
         xwatsat=rsatmix(xpress,xtemp)
         a(i,j,k)=tdewpt(xpress,min(a(i,j,k),xwatsat) )
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_dewK

!##############################################################################
Subroutine rams_comp_thete (n1,n2,n3,a,b,c)

use rconstants, only:cp,cpor,p00,alvl

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c
real :: xpress,xtemp,xwatsat
real, external :: tdewpt,rsatmix

do k=1,n3
   do j=1,n2
      do i=1,n1
         xpress=(b(i,j,k)/cp)**cpor*p00
         xtemp=c(i,j,k)*b(i,j,k)/cp
         xwatsat=rsatmix(xpress,xtemp)
         a(i,j,k)=c(i,j,k)*exp( alvl*xwatsat  &
              /(cp*tdewpt(xpress,min(a(i,j,k),xwatsat) )) )
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_thete

!##############################################################################
Subroutine rams_comp_thetv (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*(1. + .61 * b(i,j,k))
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_thetv

!##############################################################################
Subroutine rams_comp_thetrho (n1,n2,n3,a,b,c,d)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c,d

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*(1. + .608 * b(i,j,k) - c(i,j,k) - d(i,j,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_thetrho

!##############################################################################
Subroutine rams_comp_thetrho_buoy (n1,n2,n3,a,b,c,d)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c,d
real, dimension(n3) :: emean

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*(1. + .608 * b(i,j,k) - c(i,j,k) - d(i,j,k))
      enddo
   enddo
enddo

do k=1,n3
   emean(k)=0.0
   do j=1,n2
      do i=1,n1
         emean(k)=emean(k)+a(i,j,k)
      enddo
   enddo
   emean(k)=emean(k)/(n1*n2)
enddo

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=9.81*(a(i,j,k)-emean(k))/a(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_thetrho_buoy

!##############################################################################
Subroutine rams_comp_bowen (n1,n2,n3,a,b)

use rconstants, only:cp,alvl

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)/max(1.0e-12,b(i,j,k))*cp/alvl
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_bowen

!##############################################################################
Subroutine rams_comp_rh (n1,n2,n3,a,b,c)

use rconstants, only:cp,cpor,p00

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c
real :: xtemp,xpress
real, external :: rsatmix

do k=1,n3
   do j=1,n2
      do i=1,n1
         xtemp=c(i,j,k)*b(i,j,k)/cp
         xpress=(b(i,j,k)/cp)**cpor*p00
         a(i,j,k)=100.*a(i,j,k)/rsatmix(xpress,xtemp)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_rh

!##############################################################################
Subroutine rams_comp_relvortx (n1,n2,n3,a,b,c,topt,ngrd)

use mem_grid
use mem_scratch, only:vctr1,vctr2

implicit none

integer :: n1,n2,n3,i,j,k,ngrd,j1,j2,k1,k2
real, dimension(n1,n2,n3) :: a,b,c
real, dimension(n1,n2) :: topt
real :: factor

factor = ztn(1,ngrd) / ztn(2,ngrd)
do j=1,n2
   do i=1,n1
      a(i,j,1) = a(i,j,2) * factor
   enddo
enddo

CALL gradr (n1,n2,n3,2,n1-1,1,n2-1,b,c,'ydir','wpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)
CALL gradr (n1,n2,n3,2,n1-1,1,n2-1,a,b,'zdir','vpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)

do k=1,n3
   do j=1,n2
      do i=1,n1
         b(i,j,k) = c(i,j,k) - b(i,j,k)
      enddo
   enddo
enddo

do j = 1,n2
   j1 = max(j-1,1)
   j2 = min(j,n2-1)
   do i = 1,n1
      do k = 1,n3
         k1 = max(k-1,1)
         k2 = min(k,n3-1)
         a(i,j,k) =0.25 * (b(i,j1,k1) + b(i,j1,k2)  &
                         + b(i,j2,k1) + b(i,j2,k2))
      enddo
   enddo
enddo
  
return
END SUBROUTINE rams_comp_relvortx

!##############################################################################
Subroutine rams_comp_relvorty (n1,n2,n3,a,b,c,topt,ngrd)

use mem_grid
use mem_scratch, only:vctr1,vctr2

implicit none

integer :: n1,n2,n3,i,j,k,ngrd,i1,i2,k1,k2
real, dimension(n1,n2,n3) :: a,b,c
real, dimension(n1,n2) :: topt
real :: factor

factor = ztn(1,ngrd) / ztn(2,ngrd)
do j=1,n2
   do i=1,n1
      a(i,j,1) = a(i,j,2) * factor
   enddo
enddo

CALL gradr (n1,n2,n3,1,n1-1,2,n2-1,b,c,'xdir','wpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)
CALL gradr (n1,n2,n3,1,n1-1,2,n2-1,a,b,'zdir','upnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)

do k=1,n3
   do j=1,n2
      do i=1,n1
         b(i,j,k) = b(i,j,k) - c(i,j,k)
      enddo
   enddo
enddo

do j = 1,n2
   do i = 1,n1
      i1 = max(i-1,1)
      i2 = min(i,n1-1)
      do k = 1,n3
         k1 = max(k-1,1)
         k2 = min(k,n3-1)
         a(i,j,k) = 0.25 * (b(i1,j,k1) + b(i1,j,k2)  &
                          + b(i2,j,k1) + b(i2,j,k2))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_relvorty

!##############################################################################
Subroutine rams_comp_relvortz (n1,n2,n3,a,b,c,topt,ngrd)

use mem_grid
use mem_scratch, only:vctr1,vctr2

implicit none

integer :: n1,n2,n3,i,j,k,ngrd,j1,j2,i1,i2
real, dimension(n1,n2,n3) :: a,b,c
real, dimension(n1,n2) :: topt
real :: factor

factor = ztn(1,ngrd) / ztn(2,ngrd)
do j=1,n2
   do i=1,n1
      a(i,j,1) = a(i,j,2) * factor
      b(i,j,1) = b(i,j,2) * factor
   enddo
enddo

CALL gradr (n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)
CALL gradr (n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)

do k=1,n3
   do j=1,n2
      do i=1,n1
         b(i,j,k) = c(i,j,k) - b(i,j,k)
      enddo
   enddo
enddo

do j = 1,n2
   j1 = max(j-1,1)
   j2 = min(j,n2-1)
   do i = 1,n1
      i1 = max(i-1,1)
      i2 = min(i,n1-1)
      do k = 1,n3
         a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                          + b(i2,j1,k) + b(i2,j2,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_relvortz

!##############################################################################
Subroutine rams_comp_totvortz (n1,n2,n3,a,b,c,topt,ngrd)

use mem_grid
use mem_scratch, only:vctr1,vctr2
use rconstants, only:pi180

implicit none

integer :: n1,n2,n3,i,j,k,ngrd,j1,j2,i1,i2
real, dimension(n1,n2,n3) :: a,b,c
real, dimension(n1,n2) :: topt
real :: factor,fcor,omega2,xlat,xlon

factor = ztn(1,ngrd) / ztn(2,ngrd)
do j=1,n2
   do i=1,n1
      a(i,j,1) = a(i,j,2) * factor
      b(i,j,1) = b(i,j,2) * factor
   enddo
enddo

CALL gradr (n1,n2,n3,1,n1-1,1,n2-1,b,c,'xdir','vpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)
CALL gradr (n1,n2,n3,1,n1-1,1,n2-1,a,b,'ydir','upnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)

do k=1,n3
   do j=1,n2
      do i=1,n1
         b(i,j,k) = c(i,j,k) - b(i,j,k)
      enddo
   enddo
enddo

do j = 1,n2
   j1 = max(j-1,1)
   j2 = min(j,n2-1)
   do i = 1,n1
      i1 = max(i-1,1)
      i2 = min(i,n1-1)
      do k = 1,n3
         a(i,j,k) = 0.25 * (b(i1,j1,k) + b(i1,j2,k)  &
                          + b(i2,j1,k) + b(i2,j2,k))
      enddo
   enddo
enddo

omega2 = 2. * 7.292e-5
do j = 1,n2
   do i = 1,n1
      CALL xy_ll (xlat,xlon,polelat,polelon  &
           ,xtn(i,ngrd),ytn(j,ngrd))
      fcor = omega2 * sin(xlat * pi180)
      do k = 1,n3
         a(i,j,k) = a(i,j,k) + fcor
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_totvortz

!##############################################################################
Subroutine rams_comp_potvortz (n1,n2,n3,a,b,c,e,topt,ngrd)

use mem_grid
use mem_scratch, only:vctr1,vctr2

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: a,b,c,e
real, dimension(n1,n2) :: topt

CALL gradr (n1,n2,n3,1,n1-1,1,n2-1,b,e,'zdir','tpnt',topt  &
   ,xmn(1,ngrd),xtn(1,ngrd),ymn(1,ngrd),ytn(1,ngrd)  &
   ,zmn(1,ngrd),ztn(1,ngrd),deltaxn(ngrd)  &
   ,dzmn(1,ngrd),dztn(1,ngrd),vctr1,vctr2,zmn(nnzp(1)-1,1)  &
   ,jdim,ihtran,polelat,polelon)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k) = a(i,j,k) * e(i,j,k) / (9.8 * c(i,j,k))
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_potvortz

!##############################################################################
Subroutine rams_comp_vegclass (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i
real, dimension(n1,n2,n3) :: a

do i=1,n1
   a(i,1,1) = a(i,1,1) + .1
enddo
   
return
END SUBROUTINE rams_comp_vegclass

!##############################################################################
Subroutine rams_comp_horizdiv (n1,n2,n3,a,topt,ngrd)

use mem_grid, only:ztop,zmn,nnzp,dztn

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: a
real, dimension(n1,n2) :: topt
real :: rtgt

ztop = zmn(nnzp(1)-1,1)
do k=n3,2,-1
   do j=1,n2
      do i=1,n1
         rtgt = 1. - topt(i,j) / ztop
         a(i,j,k)=-(a(i,j,k)-a(i,j,k-1))*dztn(k,ngrd)/rtgt
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_horizdiv

!##############################################################################
Subroutine rams_comp_vertmin (n1,n2,n3,a,c)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c

do j = 1,n2
   do i = 1,n1
      a(i,j,1) = 0.
      do k = 2,n3-5
         a(i,j,1) = min(a(i,j,1),c(i,j,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_vertmin

!##############################################################################
Subroutine rams_comp_vertmax (n1,n2,n3,a,c)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c

do j = 1,n2
   do i = 1,n1
      a(i,j,1) = 0.
      do k = 2,n3-5
         a(i,j,1) = max(a(i,j,1),c(i,j,k))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_vertmax

!##############################################################################
Subroutine rams_comp_vertavg (n1,n2,n3,a,c)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c

do j = 1,n2
   do i = 1,n1
      a(i,j,1) = 0.
      do k = 2,n3-5
         a(i,j,1) = a(i,j,1) + c(i,j,k)
      enddo
      a(i,j,1) = a(i,j,1) / (n3-5-1)
   enddo
enddo

return
END SUBROUTINE rams_comp_vertavg

!##############################################################################
Subroutine rams_comp_vertint (n1,n2,n3,a,topt,ngrd)

use mem_grid, only:ztop,zmn,nnzp

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: a
real, dimension(n1,n2) :: topt
real :: rtgt

ztop = zmn(nnzp(1)-1,1)
do j = 1,n2
   do i = 1,n1
      rtgt = 1. - topt(i,j) / ztop
      a(i,j,1) = 0.
      do k = 2,n3-1
         a(i,j,1) = a(i,j,1) + a(i,j,k) * (zmn(k,ngrd)-zmn(k-1,ngrd)) * rtgt
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_vertint

!##############################################################################
Subroutine rams_comp_cloudtop_sstbase (n1,n2,n3,a,b,c)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b,c

do j = 1,n2
   do i = 1,n1
      a(i,j,1)= c(i,j,1)
      do k = n3,2,-1
        if(a(i,j,k) > 5.e-4) then
          a(i,j,1) = b(i,j,k)
          goto 11
        endif
      enddo
      11 continue
   enddo
enddo

return
END SUBROUTINE rams_comp_cloudtop_sstbase

!##############################################################################
Subroutine rams_comp_cloudtop_nobase (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,b

do j = 1,n2
   do i = 1,n1
      a(i,j,1)= -999.
      do k = n3,2,-1
        if(a(i,j,k) > 5.e-4) then
          a(i,j,1) = b(i,j,k)
          goto 12
        endif
      enddo
      12 continue
   enddo
enddo

return
END SUBROUTINE rams_comp_cloudtop_nobase

!##############################################################################
Subroutine rams_comp_dn0 (n1,n2,n3,a,b,c,topt,ngrd) 

use mem_grid, only:zmn,ztn,ztop,dzmn,nnzp
use mem_scratch, only:vctr2,vctr11,vctr12
use rconstants, only:cp,cpor,p00,g,rgas
use ref_sounding, only:pi01dn,th01dn

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: a,b,c
real, dimension(n1,n2) :: topt
real :: c1,c2,c3

ztop = zmn(nnzp(1)-1,1)
do j=1,n2
   do i=1,n1
      do k=1,n3
         vctr2(k)=ztn(k,ngrd)*(1.-topt(i,j)/ztop)+topt(i,j)
      enddo
      CALL htint (n3,pi01dn(1,ngrd),ztn(1,ngrd),n3,vctr11,vctr2)
      CALL htint (n3,th01dn(1,ngrd),ztn(1,ngrd),n3,vctr12,vctr2)

      do k=1,n3
         b(i,j,k)=vctr12(k)
      enddo
      a(i,j,n3) = vctr11(n3)

      c1=g*2.*(1.-topt(i,j)/ztop)
      c2=(1-cpor)
      c3=cp**c2
      do k=n3-1,1,-1
         a(i,j,k)=a(i,j,k+1)  &
             +c1/((b(i,j,k)+b(i,j,k+1))*dzmn(k,ngrd))
      enddo

      do k=1,n3
         c(i,j,k)=(c3*p00)/(rgas*b(i,j,k)*a(i,j,k)**c2)
      enddo

   enddo
enddo
   
return
END SUBROUTINE rams_comp_dn0

!##############################################################################
Subroutine rams_comp_ppress (n1,n2,n3,a,c)

use rconstants, only:cp,cpor

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k) = 1000. * (a(i,j,k)/cp) ** cpor  &
                  - 1000. * (c(i,j,k)/cp) ** cpor
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_ppress

!##############################################################################
Subroutine rams_comp_raintemp (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k) = (a(i,j,k) - 334000.) / 4186.
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_raintemp

!##############################################################################
Subroutine rams_comp_qwtc (n1,n2,n3,a,b,c)

use leaf_coms, only:slcpd

implicit none

integer :: n1,n2,n3,i,j,k,nsoil
real, dimension(n1,n2,n3) :: a,b,c
real :: qwliq0,dryhcap

do k=1,n3
   do j=1,n2
      do i=1,n1

         !a=soilenergy(J/m3), b=soilwater(m3/m3), c=soilclass
         !volumetric specific heat of liquid (4.186e6) (J/K/m3)
         !volumetric specific heat of ice (2.093e6) (J/K/m3)

         !soil_water(m3/m3) * volumetric latent heat of freezing(J/m3)
         qwliq0 = b(i,j,k) * 3.34e8 

         !soil class
         nsoil = nint(c(i,j,k))

         !volumetric dry heat capacity of soil (J/K/m3)
         dryhcap = slcpd(nsoil)

         if (a(i,j,k) .le. 0.) then
            a(i,j,k) = a(i,j,k) / (2.093e6 * b(i,j,k) + dryhcap)
            b(i,j,k) = 0.

         elseif (a(i,j,k) .ge. qwliq0) then
            a(i,j,k) = (a(i,j,k) - qwliq0) / (4.186e6 * b(i,j,k) + dryhcap)
            b(i,j,k) = 1.

         else
            b(i,j,k) = a(i,j,k) / qwliq0
            a(i,j,k) = 0.
         endif

      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_qwtc

!##############################################################################
Subroutine rams_comp_copysst (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k) = (a(i,j,n3) - 334000.) / 4186.
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_copysst

!##############################################################################
Subroutine rams_comp_qtcpcp (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         if (a(i,j,k) .le. 0.) then
            a(i,j,k) = a(i,j,k) / 2093.
         elseif (a(i,j,k) .le. 334000.) then
            a(i,j,k) = 0.
         else
            a(i,j,k) = (a(i,j,k) - 334000.) / 4186.
         endif
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_qtcpcp

!##############################################################################
Subroutine rams_comp_qtcpcp_ps (n1,n2,n3,n4,a)

implicit none

integer :: n1,n2,n3,n4,k,ip,i,j
real, dimension(n1,n2,n3,n4) :: a

do j=1,n2
  do i=1,n1
    do ip=2,n4
      do k=1,n3
         if (a(i,j,k,ip) .le. 0.) then
            a(i,j,k,ip) = a(i,j,k,ip) / 2093.
         elseif (a(i,j,k,ip) .le. 334000.) then
            a(i,j,k,ip) = 0.
         else
            a(i,j,k,ip) = (a(i,j,k,ip) - 334000.) / 4186.
         endif
      enddo
    enddo

  enddo
enddo
   
return
END SUBROUTINE rams_comp_qtcpcp_ps

!##############################################################################
Subroutine rams_comp_fracliq (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         if (a(i,j,k) .le. 0.) then
            a(i,j,k) = 0.
         elseif (a(i,j,k) .ge. 334000.) then
            a(i,j,k) = 1.
         else
            a(i,j,k) = a(i,j,k) / 334000.
         endif
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_fracliq

!##############################################################################
Subroutine rams_comp_fracice (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         if (a(i,j,k) .le. 0.) then
            a(i,j,k) = 1.
         elseif (a(i,j,k) .ge. 334000.) then
            a(i,j,k) = 0.
         else
            a(i,j,k) = 1. - a(i,j,k) / 334000.
         endif
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_fracice

!##############################################################################
Subroutine rams_comp_hydrodiam (n1,n2,n3,a,c,ccfmas,ppwmas)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c
real :: ccfmas,ppwmas,rpwmas

rpwmas = 1. / ppwmas
do k=1,n3
   do j=1,n2
      do i=1,n1
         if(a(i,j,k) .gt. 1.e-10 .and. c(i,j,k).gt.1.e-10)then
            a(i,j,k) = (a(i,j,k) / (c(i,j,k) * ccfmas))**rpwmas
         else
            a(i,j,k) = 0.
         endif
      enddo
   enddo
enddo
 
return
END SUBROUTINE rams_comp_hydrodiam

!##############################################################################
Subroutine rams_comp_hydrogamma (n1,n2,n3,a,c,ccfmas,ppwmas,gnu,z)

implicit none

integer :: n1,n2,n3,i,j,k,ibin,z
integer, parameter :: nbins=1000
real, dimension(n1,n2,n3) :: a,c
real, dimension(nbins) :: db,fmg
real, external :: gammln
real :: ccfmas,ppwmas,rpwmas,SumMd,SumMdD,dD,dmean,nc,gammaa,dn,totfmg &
       ,gnu,totcx,totrx,Md,D,Dm,mixr,numc,sigma_m,Sum2nd,D0,Nw,liwc &
       ,pwcoef,dens

!Output radar parameters based on Williams et al. (2014)
!"Describing the shape of raindrop size distributions using
!uncorrelated raindrop mass spectrum parameters."
!JAMC, 2014, Vol. 53, 1282-1296.

!Saleeby:some older simulations required fixed gnu because gnu was
!not in the header files
!gnu=2.

do k=1,n3
   do j=1,n2
      do i=1,n1

 !zero out gamma variables
 Md      = 0.
 D       = 0.
 SumMd   = 0.
 SumMdD  = 0.
 Dm      = 0.
 Sum2nd  = 0.
 Sigma_m = 0.
 D0      = 0.
 Nw      = 0.
 totfmg  = 0.
 totcx   = 0.
 totrx   = 0.

 !"a" is hydrometeor mass concentration in kg/m3
 !"c" is hydrometeor number concentration in #/m3

 !If mixing ratio is large enough (>= .001 g/kg), proceed
 if(a(i,j,k) .ge. 0.000001)then
 
  !Mean mass diameter (meters)
  dmean = (a(i,j,k) / (c(i,j,k) * ccfmas))**(1./ppwmas)
  !convert to Williams et al. units
  dmean = dmean * 1000. !(mm)

  !Hydrometeor density (kg/m3)
  dens = (6./3.14159)*ccfmas*(dmean**(ppwmas-3))
  !convert to Williams et al. units
  dens = dens / 1000. !(g/cm3)

  !Liquid/Ice water content temporary variable (g/m3)
  liwc  = a(i,j,k) * 1000.
  !Number mixing ratio temporary variable (#/m3)
  numc  = c(i,j,k)

  !****************************************************************************
  !Compute gamma distribution info
  !compute characteristic diameter (mm)
  dn = dmean / gnu
  !compute gamma value
  gammaa = exp(gammln(gnu))
  !compute distribution bin increment (mm)
  dD = 0.02 * dn
  !determine fraction per bin based on mean diameters and shape parameter
  do ibin = 1,nbins
    db(ibin) = dD * (float(ibin) - 0.5)
    fmg(ibin) = dD * (db(ibin) / dn) ** (gnu - 1.)  &
           / (dn * gammaa) * exp(-db(ibin) / dn)
  enddo
  !**************************************************************************
  !Work on bin values to get Dm
  do ibin = 1,nbins
    Md      = liwc * fmg(ibin)      !Bin(mm) liquid/ice (g/m3)
    D       = db(ibin)              !Bin(mm) diameter (mm)
    SumMd   = SumMd + (Md*dD)
    SumMdD  = SumMdD + (Md*D*dD)
  enddo

  !Mass weighted mean diam (mm)
  Dm = SumMdD/SumMd

  !Work on bin values of mass mixing ratio and compute total from distribution
  do ibin = 1,nbins
    Md      = liwc * fmg(ibin)      !Bin(mm) liquid/ice (g/m3)
    D       = db(ibin)              !Bin(mm) diameter (mm)
    Sum2nd  = Sum2nd + (((D-Dm)**2)*Md*dD)    
  enddo

  !Standard deviation of mass spectrum (mm)
  Sigma_m = sqrt(Sum2nd/SumMd)

  !Compute volumetric mean diameter (mm)
  !Use (gnu-1) since RAMS' gnu is different than Williams et al.
  D0 = Dm * (3.67 + (gnu-1.)) / (4.0 + (gnu-1.))

  !Compute normalized intercept parameter (1/mm * 1/m3) (Dolan et al. 2018)
  Nw = log10( 3.67**4. * 1000. / (3.1415 * dens) * (liwc/D0**4.) )

  !***************************************************************************
  !Assign values to output variable "a"
  if(z==1) a(i,j,k) = Dm !mass weighted mean diameter (mm)
  if(z==2) a(i,j,k) = D0 !volumetric mean diameter (mm)
  if(z==3) a(i,j,k) = Nw !normalized intercept parameter (1/mm * 1/m3)
  if(z==4) a(i,j,k) = Sigma_m !standard deviation of mass spectrum (mm)

 ! if(z==1 .and. Dm > 2.2 .and. liwc>0.0005) print*,'Dm',Dm,D0,liwc,Nw
 ! if(z==3 .and. Nw < 2.)print*,'log10(Nw)',k,Nw,liwc
 ! if(z==4 .and. Sigma_m > 1.5)print*,'Sigma',k,i,j,Sigma_m

  !***************************************************************************
  !Compute total values of mass mixing ratio and number from distribution
  !do ibin = 1,nbins
  !  totcx = totcx + numc * fmg(ibin)
  !  totrx = totrx + liwc * fmg(ibin)
  !  totfmg = totfmg + fmg(ibin)
  !enddo

 else
  a(i,j,k) = 0.0
 endif

      enddo
   enddo
enddo
 
return
END SUBROUTINE rams_comp_hydrogamma

!##############################################################################
Subroutine rams_comp_aeromedrad (n1,n2,n3,a,c,aerodens)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c
real :: aerodens,aero_rg2rm,sigma

do k=1,n3
   do j=1,n2
      do i=1,n1
         if(a(i,j,k) .gt. 1.e-30 .and. c(i,j,k).gt.1.e-7)then

            !Compute numerator of Eq.5 in Saleeby and vdh (2013)
            a(i,j,k) = a(i,j,k)/c(i,j,k)
            a(i,j,k) = a(i,j,k)*(3.0/(4.0*3.14159*aerodens))
            a(i,j,k) = a(i,j,k)**0.333333

            !Compute demonenator of Eq.5 in Saleeby and vdh (2013)
            sigma = 1.8 !lognormal distribution width for RAMS (2000-2018)
            aero_rg2rm = exp(1.5 * alog(sigma)**2)

            !Compute final values of Eq.5
            a(i,j,k) = a(i,j,k) / aero_rg2rm

         else

            a(i,j,k) = 0.

         endif
      enddo
   enddo
enddo
 
return
END SUBROUTINE rams_comp_aeromedrad

!##############################################################################
Subroutine rams_comp_ctopdiam (n1,n2,n3,a,c,ccfmas,ppwmas)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: a,c
real :: rpwmas,ppwmas,ccfmas

rpwmas = 1. / ppwmas
do j=1,n2
   do i=1,n1
      a(i,j,1) = -9999.
      do k = n3,2,-1
         if(a(i,j,k) .gt. 1.e-5 .and. c(i,j,k).gt.1.e-2)then
            a(i,j,1) = 1.e6 * (a(i,j,k) / (c(i,j,k) * ccfmas))**rpwmas
            goto 14
         endif
      enddo
      14 continue
   enddo
enddo

return
END SUBROUTINE rams_comp_ctopdiam

!##############################################################################
Subroutine rams_comp_slmstf (n1,n2,n3,a,c)

implicit none

integer :: n1,n2,n3,i
real, dimension(n1,n2,n3) :: a,c
real :: slmsts0(12)
data slmsts0/0.395, 0.410, 0.435, 0.485, 0.451, 0.420  &
            ,0.477, 0.476, 0.426, 0.492, 0.482, 0.863/

do i=1,n1
   a(i,1,1) = a(i,1,1) / max(1.e-6,slmsts0(nint(c(i,1,1))))
enddo
   
return
END SUBROUTINE rams_comp_slmstf

!##############################################################################
Subroutine rams_sum_snowlayers_ps (n1,n2,n3,n4,a,f)

implicit none

integer :: n1,n2,n3,n4,k,ip,i,j
real, dimension(n1,n2,n3,n4) :: a
real, dimension(n1,n2,n4) :: ksum,f
real, dimension(n1,n2) :: psum

!n3 = number of snow layers (nzs), n4 = number of patches (npatch)
!a = sfcwater variable, f = patch_area

do j=1,n2
  do i=1,n1
    psum(i,j) = 0.
    do ip=2,n4
      ksum(i,j,ip) = 0.
      do k=1,n3
        ksum(i,j,ip) = ksum(i,j,ip) + a(i,j,k,ip)
      enddo
      if (f(i,j,1) .lt. .991) then
        psum(i,j) = psum(i,j) + f(i,j,ip) * ksum(i,j,ip) / (1.0-f(i,j,1))
      else
        psum(i,j) = a(i,j,1,2)
      endif
    enddo
  enddo
enddo

! Copy psum into f, which was passed in as a(1).  n3 may exceed n4 but this
! should be ok.
do j = 1,n2
   do i = 1,n1
      f(i,j,1) = psum(i,j)
   enddo
enddo
   
return
END SUBROUTINE rams_sum_snowlayers_ps

!##############################################################################
Subroutine rams_fill_sst (n1,n2,n3,kp,a,c)

implicit none

integer :: n1,n2,n3,i,j,kp
real, dimension(n1,n2,n3) :: a,c

do j=1,n2
   do i = 1,n1
      a(i,j,1) = (c(i,j,kp)- 334000.) / 4186.
   enddo
enddo
   
return
END SUBROUTINE rams_fill_sst

!##############################################################################
Subroutine rams_comp_pbl (n1,n2,n3,a,b,c,ngrd)

use mem_grid, only:ztn,zmn,nnzp

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real, dimension(n1,n2,n3) :: b,c
real, dimension(n1,n2) :: a
real :: tkethrsh,pblht

tkethrsh=0.001   ! tke threshold for PBL height in m2/s2
do j=1,n2
   do i=1,n1
      pblht=0.
      do k=2,n3
         pblht=ztn(k,ngrd)*(1.-c(i,j,1)/zmn(nnzp(1)-1,1))
         !if(i.ge.10.and.i.le.25.and.j.ge.13.and.j.le.25)  &
         !   print*,'i,j,k,z,pbl=',i,j,k,ztn(k,ngrd),pblht
         if(b(i,j,k).le.tkethrsh) goto 10
      enddo
      10 continue
      a(i,j)=pblht
   enddo
enddo

return
END SUBROUTINE rams_comp_pbl

!##############################################################################
Subroutine rams_comp_slpress (n1,n2,n3,theta,pp,z,slp)

! This routine calculates the pressure at level zlev. it
!   is hardwired here to calculate mean sea level pressure,
!   but can be easily changed to calculate pressure at any level
!   by just changing zlev.
! A standard atmosphere lapse rate of 6.5 c/km is used to
!   interpolate temperature down from level 2 in the model.

implicit none

integer :: n1,n2,n3,i,j,k,kbot,ktop
real, dimension(n1,n2,n3) :: theta,pp,z
real, dimension(n1,n2) :: slp
real :: sl_p00,sl_g,sl_cp,sl_r,sl_cpor,rlap,zlev,thbar,ddz

sl_p00=1000.
sl_g=9.8
sl_cp=1004.
sl_r=287.
sl_cpor=sl_cp/sl_r
!rlap=-.0065   ! standard temp lapse rate
rlap=.0025     ! approx standard theta lapse rate
zlev=0.
ktop=0
kbot=0

do j=1,n2
   do i=1,n1
      do k=2,n3
         if(z(i,j,k).ge.zlev) then
            ktop=k
            kbot=k-1
            goto 31
         endif
      enddo
      31 continue

      !if(i.eq.1.and.j.eq.1)  &
      !   print*,'kbot:',kbot,ktop,z(i,j,kbot),z(i,j,ktop)  &
      !         ,pp(i,j,kbot),theta(i,j,kbot)
      ddz=zlev-z(i,j,kbot)
      if(zlev.lt.z(i,j,kbot))then
         thbar=(theta(i,j,kbot)-.5*ddz*rlap)
      else
         thbar=.5*(theta(i,j,kbot)+theta(i,j,ktop))
      endif
      slp(i,j)=pp(i,j,kbot)-ddz*sl_g/thbar
      slp(i,j)=(slp(i,j)/sl_cp)**sl_cpor*sl_p00
   enddo
enddo
   
return
END SUBROUTINE rams_comp_slpress

!##############################################################################
Subroutine rams_comp_pprime (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k,band,iters,jminus,jplus,iminus,iplus,j1,i1           
real, dimension(n1,n2,n3) :: a,b,dif,bp1,bp2,iter1,iter2
real :: x1,c3,g3,sum1,sum2,ewt,wt,r2

!Set to default initialization
c3=0.;g3=0.;ewt=0.

x1=119*(b(1,2,1)-b(1,1,1))
do band=1,2
if(band.eq.1) then
 c3=-4.0*5000
 g3=0.3
endif
if(band.eq.2) then
 c3=-4.0*40000
 g3=0.4
endif
do iters=1,2

do k=1,1
do j=1,n2
 jminus=j-10
 jplus=j+10
 if(jminus.lt.1) jminus=1
 if(jplus.gt.n2) jplus=n2
 do i=1,n1
  sum1=0.0
  sum2=0.0
  iminus=i-10
  iplus=i+10
  if(iminus.lt.1) iminus=1
  if(iplus.gt.n1) iplus=n1
  do j1=jminus,jplus
   do i1=iminus,iplus
     if(j1.ne.j .and. i1.ne.i) then
      r2 = ((x1*abs(i-i1))**2.0 + (x1*abs(j-j1))**2.0)
      if(iters.eq.1) ewt=(r2/c3)
      if(iters.eq.2) ewt=(r2/(c3*g3))
      if(ewt>-300.0) then
        wt=exp(ewt)
        if(iters.eq.1) sum1=sum1 + wt*a(i1,j1,k)
        if(iters.eq.2) sum1=sum1 + wt*dif(i1,j1,k)
        sum2=sum2 + wt
      endif
     endif
   enddo
  enddo
  if(iters.eq.1) then
    iter1(i,j,k) = sum1/sum2
    dif(i,j,k) = a(i,j,k) - iter1(i,j,k)
  endif
  if(iters.eq.2) then
    iter2(i,j,k) = sum1/sum2
    if(band.eq.1) bp1(i,j,k) = iter1(i,j,k) + iter2(i,j,k)
    if(band.eq.2) then
      bp2(i,j,k) = iter1(i,j,k) + iter2(i,j,k)
      a(i,j,k) = 1.25*(bp1(i,j,k) - bp2(i,j,k))
    endif
  endif
 enddo
enddo
enddo

enddo
enddo

return
END SUBROUTINE rams_comp_pprime

!##############################################################################
Subroutine rams_comp_zero (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=0.
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_zero

!##############################################################################
Subroutine rams_comp_1plus (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_1plus

!##############################################################################
Subroutine rams_comp_1minus (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=1.-a(i,j,k)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_1minus

!##############################################################################
Subroutine rams_comp_mults (n1,n2,n3,a,s)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),s

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k) * s
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_mults

!##############################################################################
Subroutine rams_comp_accum (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)+b(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_accum

!##############################################################################
Subroutine rams_comp_noneg (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=max(a(i,j,k),0.)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_noneg

!##############################################################################
Subroutine rams_comp_nonegm (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
       !If value is very small make it zero to save output space.
       !For native hydrometeor mixing ratio this is 1.e-6 kg/kg or 1.e-3 g/kg.
       !Note that 3 rain drops of 1mm diameter ~ 0.0015 grams, so 3 rain
       !drops per kg of air would meet the sampling criteria for mixing ratio.
       if(a(i,j,k) .lt. 0.000001)then
         a(i,j,k)=0.
       endif
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_nonegm

!##############################################################################
Subroutine rams_comp_nonegp (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
       !Like "nonegm" above, we set a minimum for the microphysical process
       !rate to output from revu so we can keep files smaller. No need to 
       !analyze super small values that are typically thresholded out anyway.
       !Applying minimum of 0.0000001 kg/kg/time = 0.0001 g/kg/time

       !Apply threshold to both negative and positive values to account for
       !evaporation (negative) and deposition (positive) in the budget 
       !variables like "vapliqt".
       if( a(i,j,k) .gt. -0.0000001  .and.  a(i,j,k) .lt. 0.0000001 ) then
         a(i,j,k)=0.
       endif
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_nonegp

!##############################################################################
Subroutine rams_comp_nonegn (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3) !number concentration
real :: b(n1,n2,n3) !mixing ratio

do k=1,n3
   do j=1,n2
      do i=1,n1
       !If value is very small make it zero to save output space.
       !For native hydrometeor mixing ratio this is 1.e-6 kg/kg or 1.e-3 g/kg.
       !Note that 3 rain drops of 1mm diameter ~ 0.0015 grams, so 3 rain
       !drops per kg of air would meet the sampling criteria for mixing ratio.
       if(b(i,j,k) .lt. 0.000001)then  !hydrometeor mixing ratio (kg/kg)
         a(i,j,k)=0. !hydrometeor number concentration (#/kg)
       endif
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_nonegn

!##############################################################################
Subroutine rams_comp_subt (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)-b(i,j,k)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_subt

!##############################################################################
Subroutine rams_comp_aeroepsilon (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         !Only computes solubility fraction if aerosol mass
         !in hydrometeors is at least ~1ug/m3. Otherwise zero.
         if(b(i,j,k)>=1.0e-12) then
            a(i,j,k)=a(i,j,k)/b(i,j,k)
         else
            a(i,j,k)=0.0
         endif
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_aeroepsilon

!##############################################################################
Subroutine rams_comp_mult (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*b(i,j,k)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_mult

!##############################################################################
Subroutine rams_comp_z (n1,n2,n3,a,b,ngrd)

use mem_grid, only:ztn,zmn,nnzp

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=b(i,j,1)  &
              +ztn(k,ngrd)*(1.-b(i,j,1)/zmn(nnzp(1)-1,1))
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_z

!##############################################################################
Subroutine rams_comp_rotate (n1,n2,n3,a,b,ngrd)

use mem_grid, only:polelat,polelon,xtn,ytn

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real :: a(n1,n2,n3),b(n1,n2,n3)
real :: qlat,qlon,u,v

do k=1,n3
   do j=1,n2
      do i=1,n1
         CALL xy_ll (qlat,qlon,polelat,polelon  &
            ,xtn(i,ngrd),ytn(j,ngrd))
         u=a(i,j,k)
         v=b(i,j,k)
         CALL uvtoueve (u,v,a(i,j,k),b(i,j,k)  &
                      ,qlat,qlon,polelat,polelon)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_rotate

!##############################################################################
Subroutine rams_comp_tempK (n1,n2,n3,a,b)

use rconstants, only:cp

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*b(i,j,k)/cp
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_tempK

!##############################################################################
Subroutine rams_comp_tempC (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)-273.16
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_tempC

!##############################################################################
Subroutine rams_comp_tempF (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=(a(i,j,k)-273.16)*1.8+32.
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_tempF

!##############################################################################
Subroutine rams_comp_rslf (n1,n2,n3,a,b,c)

use rconstants, only:cp,p00,cpor

implicit none

real, external :: rslf
integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),tempk,p

do k=1,n3
   do j=1,n2
      do i=1,n1
        tempk=a(i,j,k)*b(i,j,k)/cp
        tempk=max(193.15,tempk)
        p=(b(i,j,k)/cp)**cpor*p00
        a(i,j,k) = (c(i,j,k) / rslf(p,tempk)) - 1.
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_rslf

!##############################################################################
Subroutine rams_comp_rsif (n1,n2,n3,a,b,c)

use rconstants, only:cp,p00,cpor

implicit none

real, external :: rsif
integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3),c(n1,n2,n3),tempk,p

do k=1,n3
   do j=1,n2
      do i=1,n1

       tempk=a(i,j,k)*b(i,j,k)/cp
       tempk=max(193.15,tempk)
       p=(b(i,j,k)/cp)**cpor*p00

       a(i,j,k) = (c(i,j,k) / rsif(p,tempk)) - 1.

      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_rsif

!##############################################################################
Subroutine rams_comp_vaporpress (n1,n2,n3,a,b)

use rconstants, only:cp,p00,cpor

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=(a(i,j,k)/cp)**cpor*p00*.01
         a(i,j,k)=b(i,j,k) * a(i,j,k) / (b(i,j,k) + 0.622)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_vaporpress

!##############################################################################
Subroutine rams_comp_press (n1,n2,n3,a)

use rconstants, only:cp,p00,cpor

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=(a(i,j,k)/cp)**cpor*p00*.01
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_press

!##############################################################################
Subroutine rams_comp_wcms (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=a(i,j,k)*100.
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_wcms

!##############################################################################
Subroutine rams_comp_avgw (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=n3,2,-1
   do j=1,n2
      do i=1,n1
         a(i,j,k)=0.5*(a(i,j,k)+a(i,j,k-1))
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_avgw

!##############################################################################
Subroutine rams_comp_avgu (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=n1,2,-1
         a(i,j,k)=0.5*(a(i,j,k)+a(i-1,j,k))
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_avgu

!##############################################################################
Subroutine rams_comp_avgv (n1,n2,n3,a)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3)

do k=1,n3
   do j=n2,2,-1
      do i=1,n1
         a(i,j,k)=0.5*(a(i,j,k)+a(i,j-1,k))
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_avgv

!##############################################################################
Subroutine rams_comp_sfcdiv (n1,n2,n3,a,ngrd)

use mem_grid, only:dztn

implicit none

integer :: n1,n2,n3,i,j,ngrd
real :: a(n1,n2,n3)

do j=1,n2
   do i=1,n1
      a(i,j,1)=-(a(i,j,2)-a(i,j,1))*dztn(2,ngrd)
   enddo
enddo
   
return
END SUBROUTINE rams_comp_sfcdiv

!##############################################################################
Subroutine rams_comp_speed (n1,n2,n3,a,b)

implicit none

integer :: n1,n2,n3,i,j,k
real :: a(n1,n2,n3),b(n1,n2,n3)

do k=1,n3
   do j=1,n2
      do i=1,n1
         a(i,j,k)=sqrt(a(i,j,k)**2+b(i,j,k)**2)
      enddo
   enddo
enddo
   
return
END SUBROUTINE rams_comp_speed

!##############################################################################
Subroutine rams_comp_patchsum (n1,n2,n3,n4,a,f,psum)

implicit none

integer :: n1,n2,n3,n4,i,j,k,ip
real :: a(n1,n2,n3,n4),f(n1,n2,n4),psum(n1,n2,n3)

! This routine is for quantities such as net roughness that are defined
! for all patches

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         psum(i,j,k) = 0.
         do ip = 1,n4
            psum(i,j,k) = psum(i,j,k) + f(i,j,ip) * a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo

! Copy psum into f, which was passed in as a(1).  n3 may exceed n4 but this
! should be ok.

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         f(i,j,k) = psum(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_patchsum

!##############################################################################
Subroutine rams_comp_patchsum_l (n1,n2,n3,n4,a,f,psum)

implicit none

integer :: i,j,k,n1,n2,n3,n4,ip
real :: a(n1,n2,n3,n4),f(n1,n2,n4),psum(n1,n2,n3)

! This routine is for quantities such as veg roughness that are not
! defined for water patches

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         if (f(i,j,1) .lt. .991) then
            psum(i,j,k) = 0.
            do ip = 2,n4
               psum(i,j,k) = psum(i,j,k) + f(i,j,ip) * a(i,j,k,ip)  &
                           / (1. - f(i,j,1))
            enddo
         else
            psum(i,j,k) = a(i,j,k,2)
         endif
      enddo
   enddo
enddo

! Copy psum into f, which was passed in as a(1).  n3 may exceed n4 but this
! should be ok.

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         f(i,j,k) = psum(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_patchsum_l

!##############################################################################
Subroutine rams_comp_bigpatch (n1,n2,n3,n4,a,f,b)

implicit none

integer :: n1,n2,n3,n4,i,j,k
real :: a(n1,n2,n3,n4),f(n1,n2,n4),b(n1,n2,n3)

! Extract LSP value from largest patch

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         if (f(i,j,2) .ge. f(i,j,1)) then
            b(i,j,k) = a(i,j,k,2)
         else
            b(i,j,k) = a(i,j,k,1)
         endif
      enddo
   enddo
enddo

! Copy b into f, which was passed in as a(1).  n3 may exceed n4 but this
! should be ok.

do k = 1,n3
   do j = 1,n2
      do i = 1,n1
         f(i,j,k) = b(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE rams_comp_bigpatch

!##############################################################################
Subroutine rams_comp_5050 (n1,n2,n3,a,d)

implicit none 

integer :: n1,n2,n3,i,j
real :: a(n1,n2),d(n1,n2,n3)

do j = 1,n2
   do i = 1,n1
      a(i,j) = .5 * (a(i,j) + d(i,j,2))
   enddo
enddo

return
END SUBROUTINE rams_comp_5050

!##############################################################################
Subroutine rams_reduced_temp (n1,n2,n3,n4,tempnew,speed,ustar  &
                             ,tstar,znew,zold,zrough,patfrac  &
                             ,cantemp,theta,pi,topo,ztop)

use rconstants, only:cp,g,vonk
        
implicit none

integer :: n1,n2,n3,n4,i,j,np
real :: tempnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
       ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
       ,theta(n1,n2,n3),pi(n1,n2,n3),topo(n1,n2),ztop,tstar(n1,n2,n4)
real :: richno,rtgt,zagl,rtemp,rtempw,z0,a2,spd,cantheta,sfcpi


do j=1,n2
   do i=1,n1
      
      rtgt=1.-topo(i,j)/ztop
      zagl=zold*rtgt
      sfcpi=.5*(pi(i,j,1)+pi(i,j,2))
      rtempw=0.
      
      do np=1,n4
      
         z0=zrough(i,j,np)
         if(np==1) z0=.001
         spd=max(speed(i,j,2),.25)
         cantheta=cantemp(i,j,np)*cp/sfcpi

         richno=g*zagl*(theta(i,j,2)-cantheta)  &
                     /(theta(i,j,2)*spd**2)
         a2 = (vonk / log(znew / z0)) ** 2

         if(richno.gt.0.) then
            rtemp=cantheta                            &
             +(ustar(i,j,np)*tstar(i,j,np)*0.74)/(a2*spd)  &
                    *(1.+15.*richno*sqrt(1+5*richno))  
            rtemp=min(max(rtemp, cantheta),theta(i,j,2))
         else
            rtemp=cantheta                              &
             +((ustar(i,j,np)*tstar(i,j,np)*0.74)/(a2*spd))  &
               / (1.- 15.*richno/(1.+75.*a2   &
                             * sqrt(-znew*richno/z0)))
            rtemp=max(min(rtemp, cantheta),theta(i,j,2))
         endif
         
         !if((i==50.and.j==25)) then
         !   print*,'====tempf2m:',i,j
         !   print*,np,patfrac(i,j,np),cantheta
         !   print*,np,ustar(i,j,np),zrough(i,j,np),tstar(i,j,np)
         !   print*,np,theta(i,j,2),speed(i,j,2),rtemp
         !endif
         
         rtempw=rtempw+rtemp*patfrac(i,j,np)
      
      enddo
      
      tempnew(i,j)=rtempw
      
   enddo
enddo

return
END SUBROUTINE rams_reduced_temp

!##############################################################################
Subroutine rams_reduced_wind (n1,n2,n3,n4,velnew,speed,ustar  &
                             ,znew,zold,zrough,patfrac,cantemp  &
                             ,theta,pi,topo,ztop)

use rconstants, only:cp,g,vonk

implicit none

integer :: n1,n2,n3,n4,i,j,np
real :: velnew(n1,n2),speed(n1,n2,n3),ustar(n1,n2,n4),znew,zold  &
          ,zrough(n1,n2,n4),patfrac(n1,n2,n4),cantemp(n1,n2,n4)  &
          ,theta(n1,n2,n3),pi(n1,n2,n3),topo(n1,n2),ztop
real:: richno,rtgt,zagl,rwind,rwindw,z0,a2,spd,cantheta,sfcpi

do j=1,n2
   do i=1,n1
      
      rtgt=1.-topo(i,j)/ztop
      zagl=zold*rtgt
      sfcpi=.5*(pi(i,j,1)+pi(i,j,2))
      
      rwindw=0.
      
      do np=1,n4
      
         z0=zrough(i,j,np)
         if(np==1) z0=.001
         spd=max(speed(i,j,2),.25)
         cantheta=cantemp(i,j,np)*cp/sfcpi

         richno=g*zagl*(theta(i,j,2)-cantheta)  &
                      /(theta(i,j,2)*spd**2)
         a2 = (vonk / log(znew / z0)) ** 2

         if(richno.gt.0.) then
            rwind=sqrt(ustar(i,j,np)**2/a2   &
                     *(1.+10.*richno/sqrt(1+5*richno)) )
         else
            rwind=sqrt( ustar(i,j,np)**2/a2  &
                / (1.- 10.*richno/(1.+75.*a2  &
                              * sqrt(-znew*richno/z0))))
         endif
         
         rwind=max(min(rwind,speed(i,j,2)),0.)
         
         !if(i==50.and.j==25) then
         !   print*,'====speed10m'
         !   print*,np,patfrac(i,j,np),cantemp(i,j,np)
         !   print*,np,ustar(i,j,np),zrough(i,j,np)
         !   print*,np,theta(i,j,2),speed(i,j,2),rwind
         !endif
         
         rwindw=rwindw+rwind*patfrac(i,j,np)
      
      enddo
      
      velnew(i,j)=rwindw
      
   enddo
enddo

return
END SUBROUTINE rams_reduced_wind

!##############################################################################
Subroutine rams_net_rad_flx (n1,n2,n3,netflx,swup,lwup,swdn,rshort,rlong  &
                            ,rlongup,albedt)

implicit none

integer :: n1,n2,n3,i,j
real, dimension(n1,n2) :: netflx,rshort,rlong,rlongup,albedt
real, dimension(n1,n2,n3) :: swup,lwup,swdn
real :: netuptoa,netupsfc,rshortup

!Note: + values indicate cooling
!Note: - values indicate warming

do j=1,n2
   do i=1,n1
     netuptoa = swup(i,j,n3) + lwup(i,j,n3) - swdn(i,j,n3)
     rshortup = rshort(i,j) * albedt(i,j)
     netupsfc = (rshortup + lwup(i,j,1)) - (rshort(i,j) + rlong(i,j))
     netflx(i,j)= netuptoa - netupsfc
   enddo
enddo

return
END SUBROUTINE rams_net_rad_flx

!##############################################################################
Subroutine rams_heatrate (n1,n2,n3,heatrt,lswup,lswdn,dn0,pi,ngrd)

use mem_grid, only:dztn
use rconstants, only:cp

implicit none

integer :: n1,n2,n3,i,j,k,ngrd
real :: divisor
real, dimension(n1,n2,n3) :: heatrt,lswup,lswdn,dn0,pi

!Note: summing as (UP - DOWN)

do j=1,n2
 do i=1,n1
  do k=2,n3-2
     divisor = dn0(i,j,k) * (1./dztn(k,ngrd)) * cp * pi(i,j,k)/cp
     heatrt(i,j,k) = &
        (lswdn(i,j,k)-lswdn(i,j,k-1)+lswup(i,j,k-1)-lswup(i,j,k))/divisor
  enddo
  heatrt(i,j,1)=heatrt(i,j,2)
  heatrt(i,j,n3)=heatrt(i,j,n3-2)!issue with not having McLatchy levels post-run
  heatrt(i,j,n3-1)=heatrt(i,j,n3-2)!issue with not having McLatchy levels post-run
 enddo
enddo

return
END SUBROUTINE rams_heatrate

!##############################################################################
Subroutine rams_sum_rad_flx (n1,n2,n3,netflx,swup,lwup,swdn,lwdn)

implicit none

integer :: n1,n2,n3,i,j,k
real, dimension(n1,n2,n3) :: netflx,swup,lwup,swdn,lwdn

!Note: summing as (UP - DOWN)

do j=1,n2
 do i=1,n1
  do k=1,n3
     netflx(i,j,k)=swup(i,j,k)+lwup(i,j,k)-swdn(i,j,k)-lwdn(i,j,k)
  enddo
 enddo
enddo

return
END SUBROUTINE rams_sum_rad_flx

!##############################################################################
Subroutine cldfraction (n1,n2,n3,frac,pi,rh)

! calculate cloud fraction - actually returns clear sky fraction

implicit none

integer :: i,j,k,kmax,n1,n2,n3
real :: frac(n1,n2),pi(n1,n2,n3),rh(n1,n2,n3)
real, allocatable :: rhc(:), cs(:)
real :: kappai,c_1,c_2,c_junk,pop2,csmax

c_1    = 2.
c_junk = 3.
c_2    = c_junk**0.5
kappai = (1./.286)

allocate (rhc(n3),cs(n3))

do j=1,n2
   do i=1,n1
      frac(i,j) = 0.
      csmax = 0.
      kmax  = 0
      do k = 1, n3
         rhc(k)= 0.
         cs(k)= 0.
      enddo

      do k = 1, n3
         pop2 = (pi(i,j,k)/pi(i,j,2))**kappai

         rhc(k) = 100. - (100.*c_1*pop2)*(1.-pop2)*(1.+c_2*(pop2-0.5))

         if(rh(i,j,k) .ge. rhc(k))then
            if(rhc(k).eq.100.)rhc(k)=rhc(k)+0.0000001
            cs(k) = ( (rh(i,j,k)-rhc(k))/(100.-rhc(k)) ) **2. 
         else
            cs(k) = 0.
         endif
         if(cs(k).gt.csmax)then
            csmax=cs(k)
            kmax = k
         endif
         frac(i,j) = frac(i,j) + cs(k)*(1./float(k))
         !if(i==20.and.j==20) print*,k,pi(i,j,k),rh(i,j,k),frac(i,j)
      enddo
      
      csmax=max(csmax,0.)

      !frac(i,j) = 1.-min(1.,max(0.,csmax))
      frac(i,j) = 1.-min(1.,max(0.,frac(i,j)))
   enddo
enddo

deallocate (rhc,cs)

return
END SUBROUTINE cldfraction

!##############################################################################
Subroutine reflectivity_zero (kace,iace,jace &
  ,temprmix,temprNt,tempgmix,tempgNt,temphmix,temphNt &
  ,temppmix,temppNt,tempsmix,tempsNt,tempamix,tempaNt)

implicit none

integer :: i,j,k,iace,jace,kace

real,dimension(kace,iace,jace) :: temprmix,temprNt
real,dimension(kace,iace,jace) :: tempgmix,tempgNt
real,dimension(kace,iace,jace) :: temphmix,temphNt
real,dimension(kace,iace,jace) :: temppmix,temppNt
real,dimension(kace,iace,jace) :: tempsmix,tempsNt
real,dimension(kace,iace,jace) :: tempamix,tempaNt

do j=1,jace
 do i=1,iace
  do k=1,kace
      temprmix(k,i,j)=0.0
      tempgmix(k,i,j)=0.0
      temphmix(k,i,j)=0.0
      temppmix(k,i,j)=0.0
      tempsmix(k,i,j)=0.0
      tempamix(k,i,j)=0.0
      temprNt(k,i,j)=0.0
      tempgNt(k,i,j)=0.0
      temphNt(k,i,j)=0.0
      temppNt(k,i,j)=0.0
      tempsNt(k,i,j)=0.0
      tempaNt(k,i,j)=0.0
  enddo
 enddo
enddo

return
END SUBROUTINE reflectivity_zero

!##############################################################################
Subroutine reflectivity_all (kace,iace,jace &
  ,temprmix,temprNt,tempgmix,tempgNt,temphmix,temphNt &
  ,temppmix,temppNt,tempsmix,tempsNt,tempamix,tempaNt &
  ,tempdens,reflcp)

implicit none

! DESCRIPTION:
! This routine diagnoses radar reflectivity (dBZ) for rain, graupel
! and hail distributions with constant gamma distribution shape
! parameters (gnu) only and assumes hydrometeors are spherical!!
! Reflectivity is related to 6th moment of hydrometeor distribution
! by the following formula:
!
!   Z = A*F(v)*((r/alpha_m)^2)*(rho/Nt)   [mm^6/m^3]
! where
!  A = 1*10^18 to obtain correct units for Z
!  v = gnu (gamma dist. shape parameter)    [non-dimensional]
!  F(v) = ((5+v)(4+v)(3+v))/((2+v)(1+v)v)
!  r = hydrometeor mixing ratio       [kg/kg]
!  rho = air density (computed from pi & theta) [kg/m^3]
!  alpha_m = mass coeff for hydromet (cfmas) [hydromet dependnt]
!  Nt = hydromet # concentration            [#/kg]
!
! Radar reflectivity is simply:
!   reflcp = 10*LOG10(Z)   [dBZ]
!
! Author: Adrian Loftus 8/22/07
!
! For g & h, need to use equivalent radar reflectivity?
!  --> Ze = [Ki]^2/[Kw]^2 * (alpha_x/alpha_r)^2 * Z
!      [Ki]^2/[Kw]^2 = 0.224, ratio of dielectric constants
!      alpha_x is mass coeff for graupel or hail
! *NOTE* Ze IS CURRENTLY NOT USED! 

integer :: i,j,k,iace,jace,kace
real :: gnur,gnug,gnuh      !gamma distr shape params for rain,graupel,hail
real :: gnup,gnus,gnua      !gamma distr shape params for pris,snow,aggregates
real :: alpha_mr            !cfmas for rain as in mic_init.f90 file
real :: alpha_mg            !cfmas for graupel as in mic_init.f90 file
real :: alpha_mh            !cfmas for hail as in mic_init.f90 file
real :: alpha_mp            !cfmas for rain as in mic_init.f90 file
real :: alpha_ms            !cfmas for graupel as in mic_init.f90 file
real :: alpha_ma            !cfmas for hail as in mic_init.f90 file
real :: F_gnu1,F_gnu2       !used to compute F(v) as above
real :: F_gnur,F_gnug,F_gnuh !F(v) for rain, graupel, hail
real :: F_gnup,F_gnus,F_gnua !F(v) for pris, snow, aggregates

real :: tmp91,tmp92,tmp93,tmp94,tmp95,tmp96,tmp98,tmp99,q,qn,D,M !temporary variables
real,dimension(kace,iace,jace) :: temprmix,temprNt
real,dimension(kace,iace,jace) :: tempgmix,tempgNt
real,dimension(kace,iace,jace) :: temphmix,temphNt
real,dimension(kace,iace,jace) :: temppmix,temppNt
real,dimension(kace,iace,jace) :: tempsmix,tempsNt
real,dimension(kace,iace,jace) :: tempamix,tempaNt
real,dimension(kace,iace,jace) :: reflcp    !reflectivity values [dBZ]
real,dimension(kace,iace,jace) :: tempdens  !grid point density [kg/m^3]

!-- Set up some local parameters
!***NOTE: The shape params should match those in RAMSIN
gnur=2.0  !gamma shape parameter for rain
gnug=2.0  !gamma shape parameter for graupel
gnuh=2.0  !gamma shape parameter for hail
gnup=2.0  !gamma shape parameter for pris
gnus=2.0  !gamma shape parameter for snow
gnua=2.0  !gamma shape parameter for aggregates

q=1.E-10 !kg/kg of mixing ratio
qn=1.E-3 !#/kg of hydromets

!-- Setup distr params to use (category depndnt)
F_gnu1=(5.+gnur)*(4.+gnur)*(3.+gnur)
F_gnu2=(2.+gnur)*(1.+gnur)*gnur
F_gnur=F_gnu1/F_gnu2
F_gnu1=(5.+gnug)*(4.+gnug)*(3.+gnug)
F_gnu2=(2.+gnug)*(1.+gnug)*gnug
F_gnug=F_gnu1/F_gnu2
F_gnu1=(5.+gnuh)*(4.+gnuh)*(3.+gnuh)
F_gnu2=(2.+gnuh)*(1.+gnuh)*gnuh
F_gnuh=F_gnu1/F_gnu2
F_gnu1=(5.+gnup)*(4.+gnup)*(3.+gnup)
F_gnu2=(2.+gnup)*(1.+gnup)*gnup
F_gnup=F_gnu1/F_gnu2
F_gnu1=(5.+gnus)*(4.+gnus)*(3.+gnus)
F_gnu2=(2.+gnus)*(1.+gnus)*gnus
F_gnus=F_gnu1/F_gnu2
F_gnu1=(5.+gnua)*(4.+gnua)*(3.+gnua)
F_gnu2=(2.+gnua)*(1.+gnua)*gnua
F_gnua=F_gnu1/F_gnu2

do j=1,jace
 do i=1,iace
  do k=1,kace

   alpha_mr=524.0    !rain mass coeff
   alpha_mg=157.0    !graupel mass coeff
   alpha_mh=471.0    !hail mass coeff
   alpha_mp=110.8    !pris mass coeff
   alpha_ms=2.739e-3 !snow mass coeff
   alpha_ma=0.496    !aggregates mass coeff

   !The alpha mass-diameter coefficients above are the inital values. For 
   !rain, hail and graupel they assume spheres (ie. beta=3.0). Not so for 
   !pristine ice, snow, and aggregates. As such, below we compute an equivalent
   !alpha using an assumption that beta=3.0. This appears to eliminate the
   !problem of ice particles using too small of an alpha which provide a 
   !proxy for particle density. If particle density is too small then the
   !equations make Diameter too big. This compute equivalent alpha helps
   !solve this problem and leaves us with realistic looking dBZ.
   !Consider this a first-look dBZ. Best to use a radar simulator (Quickbeam).

   reflcp(k,i,j)=-1000. !Set radar reflctvty [dBZ] to default value initially.
                        !This is for grid points where mixing ratio
                        !or # concen of hydromet's are zero
   tmp99=0.0 !Set this to zero to start since it is accumulated

   if((temprmix(k,i,j).le.q).OR.(temprNt(k,i,j).le.qn))then
    tmp99=0.
   else
    tmp98=temprmix(k,i,j)/alpha_mr
    tmp98=F_gnur*(tmp98*tmp98)*tempdens(k,i,j)
    tmp99=1.0E18*(tmp98/temprNt(k,i,j)) ! 'Z' due to rain
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of rain only
   endif

   if((tempgmix(k,i,j).le.q).OR.(tempgNt(k,i,j).le.qn))then
    tmp92=0.0
   else
    tmp91=tempgmix(k,i,j)/alpha_mg
    tmp91=F_gnug*(tmp91*tmp91)*tempdens(k,i,j)
    tmp92=1.0E18*(tmp91/tempgNt(k,i,j)) ! 'Z' due to graupel
    tmp99=tmp99+tmp92
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of 
   endif                              !graupel and rain (if present)

   if((temphmix(k,i,j).le.q).OR.(temphNt(k,i,j).le.qn))then
    tmp93=0.0
   else
    tmp91=temphmix(k,i,j)/alpha_mh
    tmp91=F_gnuh*(tmp91*tmp91)*tempdens(k,i,j)
    tmp93=1.0E18*(tmp91/temphNt(k,i,j)) ! 'Z' due to hail
    tmp99=tmp99+tmp93
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of hail,
   endif                       !rain (if present), and graupel (if present)

   if((temppmix(k,i,j).le.q).OR.(temppNt(k,i,j).le.qn))then
    tmp94=0.0
   else
    M=temppmix(k,i,j)/temppNt(k,i,j)
    D=(M/alpha_mp)**(1.0/2.91)
    alpha_mp=M/(D**3.0)
    tmp91=temppmix(k,i,j)/alpha_mp
    tmp91=F_gnup*(tmp91*tmp91)*tempdens(k,i,j)
    tmp94=1.0E18*(tmp91/temppNt(k,i,j)) ! 'Z' due to pris
    tmp99=tmp99+tmp94
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of pris,
   endif   !rain + graupel + hail

   if((tempsmix(k,i,j).le.q).OR.(tempsNt(k,i,j).le.qn))then
    tmp95=0.0
   else
    M=tempsmix(k,i,j)/tempsNt(k,i,j)
    D=(M/alpha_ms)**(1.0/1.74)
    alpha_ms=M/(D**3.0)
    tmp91=tempsmix(k,i,j)/alpha_ms
    tmp91=F_gnus*(tmp91*tmp91)*tempdens(k,i,j)
    tmp95=1.0E18*(tmp91/tempsNt(k,i,j)) ! 'Z' due to snow
    tmp99=tmp99+tmp95
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of snow,
   endif   !rain + graupel + hail + pris

   if((tempamix(k,i,j).le.q).OR.(tempaNt(k,i,j).le.qn))then
    tmp96=0.0
   else
    M=tempamix(k,i,j)/tempaNt(k,i,j)
    D=(M/alpha_ma)**(1.0/2.40)
    alpha_ma=M/(D**3.0)
    tmp91=tempamix(k,i,j)/alpha_ma
    tmp91=F_gnua*(tmp91*tmp91)*tempdens(k,i,j)
    tmp96=1.0E18*(tmp91/tempaNt(k,i,j)) ! 'Z' due to aggr
    tmp99=tmp99+tmp96
    reflcp(k,i,j)=10.0*(LOG10(tmp99)) !reflectivity [dBZ] of aggregates,
   endif   !rain + graupel + hail + pris + snow

  enddo
 enddo
enddo

return
END SUBROUTINE reflectivity_all

!##############################################################################
Subroutine density4reflc (kace,iace,jace,den,pres,tempk)

use rconstants, only:rgas

implicit none

! DESCRIPTION:
! This routine computes gridpoint density values (not reference
! density) for use in calculating reflectivity values
!
! INPUT: 'pres' which is pressure (in mb)
!        'tempk' which is air temperature [K]
!
! OUTPUT: 'den' which is actual grid point density [kg/m^3]

integer :: i,j,k,iace,jace,kace
real,dimension(iace,jace,kace) :: tempk,pres
real,dimension(kace,iace,jace) :: den

do j=1,jace
 do i=1,iace
  do k=1,kace
   den(k,i,j)=(pres(i,j,k)*100.)/(rgas*tempk(i,j,k))
  enddo
 enddo
enddo

return
END SUBROUTINE density4reflc

!##############################################################################
Subroutine arrayswap (nn3,nn1,nn2,a_in,a_out,swaptype)

! DESCRIPTION: routine to swap arrays of dimensions (i,j,k)
!   to arrays of dimensions (k,i,j) or vice versa
!   [arrays for revu have dimensions (i,j,k) whereas
!   reflectivity routine requires dimensions to be (k,i,j)]

implicit none

integer :: i,j,k,nn1,nn2,nn3,swaptype
real,dimension(nn1,nn2,nn3) :: a_in
real,dimension(nn3,nn1,nn2) :: a_out

if(swaptype.eq.1)then
 do j=1,nn2
  do i=1,nn1
   do k=1,nn3
    a_out(k,i,j)=a_in(i,j,k)
   enddo
  enddo
 enddo
elseif(swaptype.eq.2)then
 do k=1,nn3
  do j=1,nn2
   do i=1,nn1
    a_in(i,j,k)=a_out(k,i,j)
   enddo
  enddo
 enddo
endif

return
END SUBROUTINE arrayswap

!##############################################################################
Subroutine set2zero (nn3,nn1,nn2,ar1,ar2,ar3,ar4)

implicit none 

integer :: nn1,nn2,nn3,i,j,k
real,dimension(nn3,nn1,nn2) :: ar1,ar2,ar3,ar4

do j=1,nn2
 do i=1,nn1
  do k=1,nn3
   ar1(k,i,j)=0.0
   ar2(k,i,j)=0.0
   ar3(k,i,j)=0.0
   ar4(k,i,j)=0.0
  enddo
 enddo
enddo

return
END SUBROUTINE set2zero

