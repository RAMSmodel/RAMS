!##############################################################################
Subroutine rams_anal_init (nfile,file_prefix)

use an_header
use mem_grid
use rcommons

implicit none

integer :: nfile,nv,nfn,ln,ihr1,imin1
integer :: iitime,iiday,iimonth,iiyear,iihour,iimin,iisecs,nc
character(len=*) :: file_prefix 
character(len=4) :: ctime1
character(len=strl1) :: fpref
character(len=strl1) :: cgrid

! get the files for the chosen grid to output
print*
fpref=file_prefix
write(cgrid,'(a3,i1,a3)') '*-g',igrid,'.h5'
fpref(len_trim(fpref)+1:)=trim(cgrid)
CALL rams_filelist (fnames,fpref,nfile)

! construct arrays
do nfn=1,nfile
   nc=len_trim(fnames(nfn))-6
   fnames(nfn)=fnames(nfn)(1:nc)//'-head.txt'
   open(10,file=fnames(nfn),form='formatted')
   read(10,*) nvbtab
   allocate (anal_table(nvbtab))
   do nv=1,nvbtab
      read(10,*) anal_table(nv)%string   &
                ,anal_table(nv)%npointer  &
                ,anal_table(nv)%idim_type  &
                ,anal_table(nv)%ngrid  &
                ,anal_table(nv)%nvalues
   enddo
   CALL commio ('READ',10)
   close(10)
   
   if(nfn==1) then
      write(ctime1,'(i4.4)') itime1
      read(ctime1(1:2),*) ihr1
      read(ctime1(3:4),*) imin1
      startutc=float(ihr1)+float(imin1)/60.
   endif
   
   ln=len_trim(fnames(nfn))
   read(fnames(nfn)(ln-10:ln-9),*) iisecs
   read(fnames(nfn)(ln-12:ln-11),*) iimin
   read(fnames(nfn)(ln-14:ln-13),*) iihour
   read(fnames(nfn)(ln-14:ln-11),*) iitime
   read(fnames(nfn)(ln-17:ln-16),*) iiday
   read(fnames(nfn)(ln-20:ln-19),*) iimonth
   read(fnames(nfn)(ln-25:ln-22),*) iiyear
   ftimes(nfn)=time
   ifdates(nfn)=iiyear*10000+iimonth*100+iiday
   iftimes(nfn)=iitime

   write(*,'(a,i4,f12.2,3x,i8,a,i4.4,a,i4.4,5(a,i2.2))') &
    ' files- ',nfn,ftimes(nfn),ifdates(nfn),'/',iitime,' ',iiyear,' ',iimonth &
    ,' ',iiday,' ',iihour,' ',iimin,' ',iisecs
   nfgrids(nfn)=ngrids
   
   close(10)
   deallocate (anal_table)
enddo

return
END SUBROUTINE rams_anal_init

!##############################################################################
Subroutine rams_get_fdata (nopt,nfl,fdata)

use rcommons

implicit none

integer :: nopt,nfl
real :: fdata(*)

if(nopt.eq.1) then
   fdata(1)=ftimes(nfl)
elseif(nopt.eq.2) then
   fdata(1)=startutc
endif

return
END SUBROUTINE rams_get_fdata

!##############################################################################
Subroutine rams_get_idata (nopt,nfl,idata)

use rcommons

implicit none

integer :: nopt,nfl,idata(*)

if(nopt.eq.1) then
   idata(1)=nfgrids(nfl)
elseif(nopt.eq.2) then
   idata(1)=ifdates(nfl)
elseif(nopt.eq.3) then
   idata(1)=iftimes(nfl)
endif

return
END SUBROUTINE rams_get_idata

!##############################################################################
Subroutine rams_get_cdata (nopt,nfl,cdata)

use rcommons

implicit none

integer :: nopt,nfl
character(len=*) :: cdata(*)

if(nopt.eq.0) then
   cdata(1)=fnames(nfl)
endif

return
END SUBROUTINE rams_get_cdata
