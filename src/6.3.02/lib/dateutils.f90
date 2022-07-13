!##############################################################################
Subroutine date_abs_secs (indate1,seconds)

implicit none

character(len=14) :: indate1
real(kind=8) :: seconds

! compute number of seconds past 1 January 1900 12:00 am

real(kind=8) :: s1,s2,s3,s4
integer :: year1,month1,date1,hour1,iy,ndays
integer, external :: julday

CALL date_unmake_big (year1,month1,date1,hour1,indate1)

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1= dble(ndays) *86400.
s2= dble(hour1/10000)*3600.
s3= dble(mod(hour1,10000)/100)*60.
s4= dble(mod(hour1,100))
seconds= s1+s2+s3+s4

return
END SUBROUTINE date_abs_secs

!##############################################################################
Subroutine date_abs_secs2 (year1,month1,date1,hour1,seconds)

implicit none

real(kind=8) :: seconds

! compute number of seconds past 1 January 1900 12:00 am

real(kind=8) :: s1,s2,s3,s4
integer :: year1,month1,date1,hour1,iy,ndays
integer, external :: julday

iy = year1 - 1900
ndays = iy * 365 + (iy-1)/4 + julday(month1,date1,iy)
s1= dble(ndays) *86400.
s2= dble(hour1/10000)*3600.
s3= dble(mod(hour1,10000)/100)*60.
s4= dble(mod(hour1,100))
seconds= s1+s2+s3+s4

return
END SUBROUTINE date_abs_secs2

!##############################################################################
Subroutine date_secs_ymdt (seconds,iyear1,imonth1,idate1,ihour1)

implicit none

real(kind=8) :: seconds,s1
integer :: iyear1,imonth1,idate1,ihour1

! compute real time given number of seconds past 1 January 1900 12:00 am  

integer :: ny,nyr,ileap,nm,ihr,imn,isc

integer :: mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/

! Get what year it is
nyr=0
s1=seconds
do ny=0,10000
   ileap=0
   if(mod(1900+ny,4) == 0) ileap=1
   s1=s1-(365.+ileap)*86400.
   if(s1 < 0.) then
      nyr=ny
      s1=s1+(365.+ileap)*86400.
      exit
   endif
enddo
iyear1=1900+nyr

! s1 is now number of secs into the year
!   Get month
do nm=1,12
   ileap=0
   if(mod(1900+ny,4) == 0 .and. nm == 2) ileap=1
   s1=s1-(mondays(nm)+ileap)*86400.
   if(s1 < 0.) then
      s1=s1+(mondays(nm)+ileap)*86400.
      exit
   endif
enddo
imonth1=nm

! s1 is now number of secs into the month
!   Get date and time

idate1=int(s1/86400.)
s1=s1-idate1*86400.
idate1=idate1+1 ! Since date starts at 1

ihr=int(s1/3600.)
s1=s1-ihr*3600.
imn=int(s1/60.)
s1=s1-imn*60.
isc=s1
ihour1=ihr*10000+imn*100+isc

return
END SUBROUTINE date_secs_ymdt

!##############################################################################
Subroutine date_add_to_big (cindate,tinc,tunits,coutdate)

implicit none

character(len=14) :: cindate,coutdate
real :: tinc
character(len=1) :: tunits

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year

real(kind=8) :: ttinc,secs
integer :: inyear,inmonth,indate,inhour  &
          ,outyear,outmonth,outdate,outhour

! convert input time to seconds

ttinc=tinc
if(tunits.eq.'m') ttinc=tinc*60.
if(tunits.eq.'h') ttinc=tinc*3600.
if(tunits.eq.'d') ttinc=tinc*86400.
!print*,'inc:',tinc,tunits,ttinc
!print*,'big:',cindate
CALL date_unmake_big (inyear,inmonth,indate,inhour,cindate)
!print*,'big:',inyear,inmonth,indate,inhour
CALL date_abs_secs2 (inyear,inmonth,indate,inhour,secs)

!print*,'big:',secs,ttinc
secs=secs+ttinc

CALL date_secs_ymdt (secs,outyear,outmonth,outdate,outhour)
!print*,'big:',outyear,outmonth,outdate,outhour
CALL date_make_big (outyear,outmonth,outdate,outhour,coutdate)

!print*,'out stuff:',coutdate
!stop

return
END SUBROUTINE date_add_to_big

!##############################################################################
Subroutine date_add_to (inyear,inmonth,indate,inhour  &
                        ,tinc,tunits,outyear,outmonth,outdate,outhour)

implicit none

integer inyear,inmonth,indate,inhour  &
       ,outyear,outmonth,outdate,outhour
real tinc
character(len=1) :: tunits

! adds/subtracts a time increment to a date and output new date
! -> uses hhmmss for hours, 4 digit year


real(kind=8) :: ttinc,secs

! convert input time to seconds

ttinc=tinc
if(tunits.eq.'m') ttinc=tinc*60.
if(tunits.eq.'h') ttinc=tinc*3600.
if(tunits.eq.'d') ttinc=tinc*86400.
!print*,'inc:',tinc,tunits,ttinc


CALL date_abs_secs2 (inyear,inmonth,indate,inhour,secs)

secs=secs+ttinc

CALL date_secs_ymdt (secs,outyear,outmonth,outdate,outhour)

!print*,'out stuff:',outyear,outmonth,outdate,outhour

return
END SUBROUTINE date_add_to

!##############################################################################
Subroutine date_make_big (inyear,inmonth,indate,inhour,outdate)

implicit none

integer :: inyear,inmonth,indate,inhour
character(len=14) :: outdate

write(outdate(1:4),10) inyear
write(outdate(5:6),11) inmonth
write(outdate(7:8),11) indate
write(outdate(9:14),12) inhour
10 format (i4.4)
11 format (i2.2)
12 format (i6.6)

return
END SUBROUTINE date_make_big

!##############################################################################
Subroutine date_unmake_big (inyear,inmonth,indate,inhour,outdate)

implicit none

integer :: inyear,inmonth,indate,inhour
character(len=14) :: outdate

read(outdate(1:4),10) inyear
read(outdate(5:6),11) inmonth
read(outdate(7:8),11) indate
read(outdate(9:14),12) inhour
10 format (i4)
11 format (i2)
12 format (i6)

return
END SUBROUTINE date_unmake_big

!##############################################################################
Subroutine rams_dintsort (ni,chnums,cstr)

use grid_dims

implicit none

integer :: ni
character(len=14) :: chnums(*)
character(len=*) :: cstr(*)

! sort an array of character strings by an associated character field

character(len=strl1) :: cscr
character(len=14) :: mini,nscr

integer :: n,nm,nmm

nmm=0
do n=1,ni
   mini='99999999999999'
   do nm=n,ni
      if(chnums(nm).lt.mini) then
         nmm=nm
         mini=chnums(nm)
      endif
   enddo
   nscr=chnums(n)
   chnums(n)=chnums(nmm)
   chnums(nmm)=nscr
   cscr=cstr(n)
   cstr(n)=cstr(nmm)
   cstr(nmm)=cscr
enddo

return
END SUBROUTINE rams_dintsort

!##############################################################################
Subroutine rams_sort_dint3 (n1,ia1,n2,ia2,n3,ia3,nt,iall)

implicit none

integer :: n1,n2,n3,nt
character(len=14) :: ia1(*),ia2(*),ia3(*),iall(*)

!     sort 3 arrays of char's, put back in 1 array
!     copy all to output array

character(len=14) :: mini,nscr
integer :: n,nm,nmm

nmm=0
nt=0
do n=1,n1
   nt=nt+1
   iall(nt)=ia1(n)
enddo
do n=1,n2
   nt=nt+1
   iall(nt)=ia2(n)
enddo
do n=1,n3
   nt=nt+1
   iall(nt)=ia3(n)
enddo

do n=1,nt
   mini='99999999999999'
   do nm=n,nt
      if(iall(nm).lt.mini) then
         nmm=nm
         mini=iall(nm)
      endif
   enddo
   nscr=iall(n)
   iall(n)=iall(nmm)
   iall(nmm)=nscr
enddo

return
END SUBROUTINE rams_sort_dint3

!##############################################################################
Subroutine rams_unique_dint (n1,ia1)

implicit none

integer :: n1
character(len=14) :: ia1(*)
integer :: n,nt,nn

! reduce an array to get rid of duplicate entries

nt=n1
10 continue
do n=2,nt
   if(ia1(n).eq.ia1(n-1)) then
      do nn=n,nt
         ia1(nn-1)=ia1(nn)
      enddo
      nt=nt-1
      goto 10
   endif
enddo
n1=nt

return
END SUBROUTINE rams_unique_dint

!##############################################################################
integer Function julday (imonth,iday,iyear)

implicit none

integer :: imonth,iday,iyear

! compute the julian day from a normal date

julday = iday  &
      + min(1,max(0,imonth-1))*31  &
      + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4))))  &
      + min(1,max(0,imonth-3))*31  &
      + min(1,max(0,imonth-4))*30  &
      + min(1,max(0,imonth-5))*31  &
      + min(1,max(0,imonth-6))*30  &
      + min(1,max(0,imonth-7))*31  &
      + min(1,max(0,imonth-8))*31  &
      + min(1,max(0,imonth-9))*30  &
      + min(1,max(0,imonth-10))*31  &
      + min(1,max(0,imonth-11))*30  &
      + min(1,max(0,imonth-12))*31

return
END FUNCTION julday
