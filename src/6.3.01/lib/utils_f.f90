!##############################################################################
Subroutine rams_mm (indata,ni1,omin,omax)

implicit none

integer :: ni1
real :: indata(ni1),omin,omax
integer :: i

omax=indata(1)
omin=indata(1)
   do i=2,ni1
      omax=max(indata(i),omax)
      omin=min(indata(i),omin)
   enddo

return
END SUBROUTINE rams_mm

!##############################################################################
real Function walltime (wstart)

implicit none

real :: wstart
integer :: ii,ir

!Use lowercase "call" since this is a system call
call system_clock (count=ii,count_rate=ir)
walltime = float(ii)/float(ir) - wstart

return
END FUNCTION walltime

!##############################################################################
Subroutine rearrange (nzp,nxp,nyp,a,b)

implicit none

integer :: nzp,nxp,nyp
real :: a(nzp,nxp,nyp),b(nxp,nyp,nzp)
integer :: k,i,j

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(i,j,k)=a(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE rearrange

!##############################################################################
Subroutine unarrange (nzp,nxp,nyp,a,b)

implicit none

integer :: nzp,nxp,nyp
real :: a(nxp,nyp,nzp),b(nzp,nxp,nyp)
integer :: k,i,j

do i=1,nxp
   do j=1,nyp
      do k=1,nzp
         b(k,i,j)=a(i,j,k)
      enddo
   enddo
enddo

return
END SUBROUTINE unarrange

!##############################################################################
Subroutine rearrange_p (n2,n3,n4,n5,a,b)

implicit none

integer :: n2,n3,n4,n5
real :: a(n4,n2,n3,n5),b(n2,n3,n4,n5)
integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(i,j,k,ip) = a(k,i,j,ip)
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE rearrange_p

!##############################################################################
Subroutine unarrange_p (n2,n3,n4,n5,a,b)

implicit none

integer :: n2,n3,n4,n5
real :: a(n2,n3,n4,n5),b(n4,n2,n3,n5)
integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(k,i,j,ip) = a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo

return
END SUBROUTINE unarrange_p

!##############################################################################
Subroutine makefnam (fname,prefix,tinc,iyr,imn,idy,itm,type,post,fmt)

! creates standard timestamped filename

implicit none

integer :: iyr,imn,idy,itm,oyr,omn,ody,otm,ib1,ib2
character(len=*) :: fname,prefix,post,fmt
character(len=1) :: type
real :: tinc
character(len=40) :: dstring

!   print*,iyr,imn,idy,itm,tinc
if(tinc == 0.) then
   oyr=iyr ; omn=imn ; ody=idy ; otm=itm
else
   CALL date_add_to (iyr,imn,idy,itm,tinc,'s',oyr,omn,ody,otm)
!   print*,oyr,omn,ody,otm
endif

write(dstring,100) '-',type,'-',oyr,'-',omn,'-',ody,'-',otm
100 format(3a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

ib1=len_trim(prefix)
fname=prefix(1:ib1)//dstring(1:20)
if (post(1:1) /= '$') then
   ib1=len_trim(fname)
   ib2=len_trim(post)
   fname=fname(1:ib1)//'-'//post(1:ib2)
endif
ib1=len_trim(fname)
fname=fname(1:ib1)//'.'//trim(fmt)

return
END SUBROUTINE makefnam

!##############################################################################
Subroutine rams_f_open (iunit,filenm,formt,stat,act,iclob)

! replaces old jclopen and jclget
! files are overwritten unless iclob (ICLOBBER) set to 1

use mem_grid, only:print_msg

implicit none

integer :: iunit,iclob
character(len=*) :: filenm,formt,stat,act
logical :: exans

!print*,'filenm,formt,stat1=',filenm,formt,stat

inquire(FILE=filenm,EXIST=exans)

if(exans.and.iclob.eq.0.and.  &
     (act(1:4).eq.'WRIT'.or.act(1:4).eq.'writ') .and. print_msg) then
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'!!!   trying to open file name :'
   print*,'!!!       ',filenm
   print*,'!!!   but it already exists. run is ended.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop 'rams_f_open - exists'
endif

!print*,'filenm,formt,stat2=',filenm(1:len_trim(filenm)),formt,stat
open(iunit,STATUS=stat,FILE=trim(filenm),FORM=formt)
if(print_msg) print*,'F_open - ',trim(filenm)

return
END SUBROUTINE rams_f_open
