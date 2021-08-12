!##############################################################################
Subroutine commio (io,iun)

use micphys
use ref_sounding
use mem_grid
use mem_leaf
use leaf_coms
use kpp_parameters, only:nkppz
use mem_kpp, only:kpp_const_fields

implicit none

integer :: iun
character(len=*) :: io

!     This routine reads or writes the analysis file header.

integer, external :: cio_i
integer, external :: cio_i1
integer, external :: cio_f
integer, external :: cio_f1
integer, external :: cio_c
character(len=2) :: cng
integer :: irw,ie,ng

IF(IO.EQ.'READ') irw=1
IF(IO.EQ.'WRITE') irw=2

!print*,'in commio ',io,' ',iun

!Below are variables sent to analysis header file. 
!Some are needed for history restart. You can see which ones are 
!need for history restart in routine "history_start".
!Additional ones are needed for REVU post-processing.

!Experiment name and type
ie=cio_c(iun,irw,'expnme',expnme,1)
ie=cio_i1(iun,irw,'initial',initial,1)
ie=cio_i1(iun,irw,'initorig',initorig,1)
ie=cio_i1(iun,irw,'jdim',jdim,1)  !0 - NOT 2D, 1 - 2D simulations

!Time sets
ie=cio_i1(iun,irw,'iyear1',iyear1,1)
ie=cio_i1(iun,irw,'imonth1',imonth1,1)
ie=cio_i1(iun,irw,'idate1',idate1,1)
ie=cio_i1(iun,irw,'itime1',itime1,1)
ie=cio_f1(iun,irw,'time',time,1)

!Sounding State
ie=cio_i1(iun,irw,'nsndg',nsndg,1) !Number of input Sounding levels
ie=cio_f(iun,irw,'us',us,nsndg)
ie=cio_f(iun,irw,'vs',vs,nsndg)
ie=cio_f(iun,irw,'ts',ts,nsndg)
ie=cio_f(iun,irw,'thds',thds,nsndg)
ie=cio_f(iun,irw,'ps',ps,nsndg)
ie=cio_f(iun,irw,'hs',hs,nsndg)

!Grid navigation and structures
ie=cio_i1(iun,irw,'ihtran',ihtran,1)
ie=cio_i1(iun,irw,'ngrids',ngrids,1)
ie=cio_i1(iun,irw,'nzg',nzg,1)
ie=cio_i1(iun,irw,'nzs',nzs,1)
ie=cio_i1(iun,irw,'npatch',npatch,1)
ie=cio_i1(iun,irw,'nvegpat',nvegpat,1)
ie=cio_i1(iun,irw,'nkppz',nkppz,1)
ie=cio_f1(iun,irw,'ztop',ztop,1)
ie=cio_f1(iun,irw,'polelat',polelat,1)
ie=cio_f1(iun,irw,'polelon',polelon,1)
ie=cio_i(iun,irw,'nnxp',nnxp,ngrids)
ie=cio_i(iun,irw,'nnyp',nnyp,ngrids)
ie=cio_i(iun,irw,'nnzp',nnzp,ngrids)
ie=cio_i(iun,irw,'nstratx',nstratx,ngrids)
ie=cio_i(iun,irw,'nstraty',nstraty,ngrids)
ie=cio_i(iun,irw,'nxtnest',nxtnest,ngrids)
ie=cio_i(iun,irw,'ninest',ninest,ngrids)
ie=cio_i(iun,irw,'njnest',njnest,ngrids)
ie=cio_i(iun,irw,'nknest',nknest,ngrids)
ie=cio_f(iun,irw,'deltaxn',deltaxn,ngrids)
ie=cio_f(iun,irw,'deltazn',deltazn,ngrids)
ie=cio_i1(iun,irw,'nestz',nestz,1)
ie=cio_i(iun,irw,'nstratz',nstratz,nnzp(1))

!Needed for acoustic routines
ie=cio_i(iun,irw,'ngbegun',ngbegun,ngrids) !0 at initial start, 1 thereafter

!Reference state and grid stagger points
ie=cio_i1(iun,irw,'iref',iref,1)     !I-grid point used for reference state
ie=cio_i1(iun,irw,'jref',jref,1)     !J-grid point used for reference state
ie=cio_f1(iun,irw,'topref',topref,1) !Topography of reference state point
do ng=1,ngrids
   write(cng,1) ng
1       format(i2.2)
   ie=cio_f(iun,irw,'xmn'//cng,xmn(1,ng),nnxp(ng))
   ie=cio_f(iun,irw,'xtn'//cng,xtn(1,ng),nnxp(ng))
   ie=cio_f(iun,irw,'ymn'//cng,ymn(1,ng),nnyp(ng))
   ie=cio_f(iun,irw,'ytn'//cng,ytn(1,ng),nnyp(ng))
   ie=cio_f(iun,irw,'zmn'//cng,zmn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'ztn'//cng,ztn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'dzmn'//cng,dzmn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'dztn'//cng,dztn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'u01dn'//cng,u01dn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'v01dn'//cng,v01dn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'pi01dn'//cng,pi01dn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'th01dn'//cng,th01dn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'dn01dn'//cng,dn01dn(1,ng),nnzp(ng))
   ie=cio_f(iun,irw,'rt01dn'//cng,rt01dn(1,ng),nnzp(ng))
enddo

!KPP ocean model grid levels
ie=cio_f(iun,irw,'kppzm',kpp_const_fields%zm,nkppz)

!Microphysics constants for hydromets and aerosols
ie=cio_f(iun,irw,'gnu',gnu,ncat)
ie=cio_i(iun,irw,'iaero_chem',iaero_chem,aerocat)
ie=cio_f(iun,irw,'cfmas',cfmas,nhcat)
ie=cio_f(iun,irw,'pwmas',pwmas,nhcat)

!Leaf-3 surface characteristics
ie=cio_f(iun,irw,'slcpd',slcpd,nstyp)
ie=cio_f(iun,irw,'slmsts',slmsts,nstyp)
ie=cio_f(iun,irw,'slz',slz,nzg)

return
END SUBROUTINE commio

!##############################################################################
Subroutine cio_pos_file (iun,cstr,ierr)

use mem_grid

implicit none

integer :: iun,ierr
character(len=*) :: cstr
character(len=strl2) :: line,csearch

integer :: nl,nc,iend,ilen

!print*,'cio_pos:',iun,cstr

iend=0
1    continue
do nl=1,1000000
   read(iun,10,end=100) line
10      format(a)
   ilen=len(cstr)
   csearch='__'//cstr(1:ilen)
   nc=index(line,csearch(1:ilen+2) )
   !print*,'cio_pos:',nl,nc,line
   if(nc.eq.1) then
      ierr=0
      if(iprntstmt>=1 .and. print_msg) &
        print*,'---- Name found on header file:',cstr
      return
   endif
enddo

100  continue
if(iend.eq.1) then
   ierr=1
   if(iprntstmt>=1 .and. print_msg) &
     print*,'---- Name NOT found on header file:',cstr
   rewind(iun)
   return
endif
rewind(iun)
iend=1
goto 1

END SUBROUTINE cio_pos_file

!##############################################################################
integer Function cio_i (iun,irw,cstr,ia,n)

implicit none

integer :: iun,irw,n
integer :: ia(n)
character(len=*) :: cstr
integer :: nn,i

if (irw.eq.1) then
   CALL cio_pos_file (iun,cstr,cio_i)
   if(cio_i.eq.1) return
   read(iun,*) nn
   read(iun,*) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) (ia(i),i=1,n)
11      format(i6)
   cio_i = 0
endif

return
END FUNCTION cio_i

!##############################################################################
integer Function cio_i1 (iun,irw,cstr,ia,n)

implicit none

integer :: iun,irw,n
integer :: ia
character(len=*) :: cstr
integer :: nn,i

if (irw.eq.1) then
   CALL cio_pos_file (iun,cstr,cio_i1)
   if(cio_i1.eq.1) return
   read(iun,*) nn
   read(iun,*) ia
   if(nn .ne. n) stop 'Variable array length not what expected in header'
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,11) ia
11      format(i6)
   cio_i1 = 0
endif

return
END FUNCTION cio_i1

!##############################################################################
integer Function cio_f (iun,irw,cstr,ia,n)

implicit none

integer :: iun,irw,n
real ia(*)
character(len=*) :: cstr
integer :: nn,i

if (irw.eq.1) then
   CALL cio_pos_file (iun,cstr,cio_f)
   if(cio_f.eq.1) return
   read(iun,*) nn
   read(iun,*) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
!Saleeby(2012-02-07):Changed format from e16.8 to e26.18
!Need more precision in header file for history restart
   write(iun,11) (ia(i),i=1,n)
11      format(e26.18)
   cio_f = 0
endif

return
END FUNCTION cio_f

!##############################################################################
integer Function cio_f1 (iun,irw,cstr,ia,n)

implicit none

integer :: iun,irw,n
real :: ia
character(len=*) :: cstr
integer :: nn,i

if (irw.eq.1) then
   CALL cio_pos_file (iun,cstr,cio_f1)
   if(cio_f1.eq.1) return
   read(iun,*) nn
   read(iun,*) ia
   if(nn .ne. n) stop 'Variable array length not what expected in header'
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
!Saleeby(2012-02-07):Changed format from e16.8 to e26.18
!Need more precision in header file for history restart
   write(iun,11) ia
11      format(e26.18)
   cio_f1 = 0
endif

return
END FUNCTION cio_f1

!##############################################################################
integer Function cio_c (iun,irw,cstr,ia,n)

implicit none

integer :: iun,irw,n
character(len=*) :: ia(*)
character(len=*) :: cstr
integer :: nn,i

if (irw.eq.1) then
   CALL cio_pos_file (iun,cstr,cio_c)
   if(cio_c.eq.1) return
   read(iun,*) nn
   read(iun,10) (ia(i),i=1,nn)
elseif(irw.eq.2) then
   write(iun,20) cstr
20      format('__',a)
   write(iun,*) n
   write(iun,10) (trim(ia(i)),i=1,n)
10      format(a)
   cio_c = 0
endif

return
END FUNCTION cio_c
