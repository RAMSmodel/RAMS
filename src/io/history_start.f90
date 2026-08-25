!##############################################################################
Subroutine history_start ()

! This routine initializes the model from the history file

use io_params
use mem_grid
use mem_basic
use node_mod
use micro_prm, only:nkr
use kpp_parameters, only:nkppz
use ref_sounding, only:iugforce


implicit none

integer :: maxarr,ngr
integer :: ifm,icm

! Find maximum size of any array on history file. Allocate scratch array of
! this size.

maxarr=0
do ngr=1,ngridsh
   maxarr=max(maxarr,mmxp(ngr)*mmyp(ngr)*mmzp(ngr)  &
         ,mmxp(ngr)*mmyp(ngr)*nzg*npatch &
         ,mmxp(ngr)*mmyp(ngr)*nzs*npatch &
         ,mmxp(ngr)*mmyp(ngr)*mmzp(ngr)*nkr &
         ,mmxp(ngr)*mmyp(ngr)*nkppz)
enddo

! read stuff here

CALL hist_read (maxarr,trim(hfilin))

if(print_msg) print*,'back from read'

if (iugforce==2 .or. iupdsst==2) CALL readforcing ()

CALL calc_refs()

CALL calc_wsub()

do ifm = 1,ngrids
   icm = nxtnest(ifm)
   if (icm  ==  0) then
      CALL newgrid (ifm)
      CALL refs3d (mzp,mxp,myp  &
      ,basic_g(ifm)%pi0  (1,1,1),basic_g(ifm)%dn0  (1,1,1)  &
      ,basic_g(ifm)%dn0u (1,1,1),basic_g(ifm)%dn0v (1,1,1)  &
      ,basic_g(ifm)%th0  (1,1,1),grid_g(ifm)%topt  (1,1)    &
      ,grid_g(ifm)%rtgt  (1,1)  )
   endif
enddo

return
END SUBROUTINE history_start

!##############################################################################
Subroutine read_distribute_hheader (hnamein)

use an_header
use mem_grid
use ref_sounding
use node_mod

implicit none

  character (len=*) :: hnamein

  integer :: ngrids1,nnxp1(maxgrds),nnyp1(maxgrds),nnzp1(maxgrds) &
            ,nzg1,nzs1,npatch1,nkppz1

  integer :: ie,ngr,nv
  character(len=2) :: cng
  integer, external :: cio_i
  integer, external :: cio_f
  integer :: iunhd=11

! Read the input history header file and collect:
!   Environmental sounding
!   Simulation time
!   Number of grids and grid specs
!     In case RAMSIN is specifying more grids
!   List of variables contained in the history file
!
! Note that read_distribute_hheader() will set the global variable
! ngridsh which is used in the following loop and in hist_read().
! read_distribute_hheader() also sets the global variables nvbtab
! and anal_table of which read_hist uses.
!
! Note that read_distribute_hheader() will allocate anal_table, and
! hist_read() will deallocate anal_table.

  ! Have node mainnum read the file then broadcast the data to the other nodes

  if ((my_rams_num .eq. mainnum) .or. (nmachs .eq. 1)) then

    CALL rams_f_open (iunhd,trim(hnamein),'FORMATTED','OLD','READ',0)

    !Get history grid structure info so we can allocate space
    ie=cio_i(iunhd,1,'ngrids',ngrids1,1)
    ngridsh=ngrids1
    ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
    ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
    ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
    ie=cio_i(iunhd,1,'npatch',npatch1,1)
    ie=cio_i(iunhd,1,'nzg',nzg1,1)
    ie=cio_i(iunhd,1,'nzs',nzs1,1)
    ie=cio_i(iunhd,1,'nkppz',nkppz1,1)
    ie=cio_f(iunhd,1,'time',time,1)
    
    !Flag to determine if we are past first timesetp
    ie=cio_i(iunhd,1,'ngbegun',ngbegun,ngrids)
    
    !Get the sounding state (needed just for standard output)
    ie=cio_i(iunhd,1,'iref',iref,1)
    ie=cio_i(iunhd,1,'jref',jref,1)
    ie=cio_f(iunhd,1,'topref',topref,1)
    ie=cio_i(iunhd,1,'nsndg',nsndg,1)
    ie=cio_f(iunhd,1,'us',us,nsndg)
    ie=cio_f(iunhd,1,'vs',vs,nsndg)
    ie=cio_f(iunhd,1,'ts',ts,nsndg)
    ie=cio_f(iunhd,1,'thds',thds,nsndg)
    ie=cio_f(iunhd,1,'rts',rts,nsndg)
    ie=cio_f(iunhd,1,'ps',ps,nsndg)
    ie=cio_f(iunhd,1,'hs',hs,nsndg)
    if(io3flg==1)ie=cio_f(iunhd,1,'o3s',o3s,nsndg)

    !Get original simulation type
    ie=cio_i(iunhd,1,'initorig',initorig,1)
    
    ! Get the 1-d reference state
    do ngr=1,ngridsh
       write(cng,1) ngr
1           format(i2.2)
       ie=cio_f(iunhd,1,'u01dn'//cng,u01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'v01dn'//cng,v01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'pi01dn'//cng,pi01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'th01dn'//cng,th01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'dn01dn'//cng,dn01dn(1,ngr),mmzp(ngr))
       ie=cio_f(iunhd,1,'rt01dn'//cng,rt01dn(1,ngr),mmzp(ngr))
    enddo

    !  Read variable header info

    rewind(iunhd)

    read(iunhd,*) nvbtab
    allocate (anal_table(nvbtab))
    do nv=1,nvbtab
       read(iunhd,*)  anal_table(nv)%string   &
                     ,anal_table(nv)%npointer  &
                     ,anal_table(nv)%idim_type  &
                     ,anal_table(nv)%ngrid  &
                     ,anal_table(nv)%nvalues
    enddo

    close(iunhd)
  endif

  ! If more than one node, broadcast the data to the other nodes
  if (nmachs .gt. 1) then
    CALL broadcast_hist_header (ngrids1,nnxp1,nnyp1,nnzp1,nzg1 &
                               ,nzs1,nkppz1,npatch1)
  endif

return
END SUBROUTINE read_distribute_hheader

!##############################################################################
Subroutine hist_read (maxarr,hnamein)

use an_header
use var_tables
use mem_grid
use hdf5_utils
use node_mod
use micro_prm, only:nkr
use kpp_parameters, only:nkppz

implicit none

integer :: maxarr,checkhist
character(len=*) :: hnamein

integer :: ngr,npts,nc,nv,nvh,ndims,idims(4)
character(len=1) :: cgrid
character(len=strl1) :: hname
real, allocatable :: scr(:)
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks
logical :: exists

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

allocate (scr(maxarr))

!Check to see that all history grids are present at this time
checkhist=1
do ngr=1,ngridsh
  write(cgrid,'(i1)') ngr
  nc=len_trim(hnamein)
  hname=hnamein(1:nc-9)//'-g'//cgrid//'.h5'
  inquire(file=hname,exist=exists)
  if(.not. exists)then
   checkhist=0
  endif
enddo
if(checkhist==0)then
 if(print_msg)then
 print*,'Not all original grids are present at this time: ',hnamein(1:nc-9)
 print*,'Choose a history restart time in which all original'
 print*,'simulation grids are available.'
 print*,''
 endif
 stop
endif

do ngr=1,ngridsh

   ! Open file
   write(cgrid,'(i1)') ngr
   nc=len_trim(hnamein)
   hname=hnamein(1:nc-9)//'-g'//cgrid//'.h5'

   CALL shdf5_open (hname,'R',iphdf5,h5_fid)

   ! Loop through all variables
   varloop: do nvh=1,nvbtab
      if(ngr /= anal_table(nvh)%ngrid) cycle varloop
      
      ! See if variable should be read and stored
      do nv = 1,num_var(ngr)
         if(anal_table(nvh)%string == vtab_r(nv,ngr)%name) then
            ! there is a match on this grid. see if hist flag is set...
            if(iprntstmt>=1 .and. print_msg) &
               print*,'found: ',trim(anal_table(nvh)%string)
            if(vtab_r(nv,ngr)%ianal /= 1) cycle varloop
            if(iprntstmt>=1 .and. print_msg) &
               print*,'read : ',trim(anal_table(nvh)%string)
            
            ! We want it...read, maybe rearrange, and store it

            ! call and output variable array names and size in hist file
            if(iprntstmt>=1 .and. print_msg) then
              CALL shdf5_info (h5_fid,trim(anal_table(nvh)%string),ndims,idims)
              print*,'name,ndims,idims: ', trim(anal_table(nvh)%string)
              print*,ndims,idims(1:ndims)
            endif

            CALL shdf5_set_hs_select (vtab_r(nv,ngr)%idim_type,'R',ngr &
                ,mem_select,file_select,file_chunks)
            CALL shdf5_irec (h5_fid,iphdf5,trim(anal_table(nvh)%string) &
                ,mem_select,file_select,rvara=scr)
            
            npts = vtab_r(nv,ngr)%npts
            
            select case(vtab_r(nv,ngr)%idim_type)
               case(2,6) ; CALL atob (npts,scr,vtab_r(nv,ngr)%var_p)
               case(3)
                  CALL unarrange (mmzp(ngr),mmxp(ngr),mmyp(ngr) &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(4)
                  CALL unarrange_p (mmxp(ngr),mmyp(ngr),nzg,npatch &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(5)
                  CALL unarrange_p (mmxp(ngr),mmyp(ngr),nzs,npatch &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(7)
                  CALL unarrange_p (mmxp(ngr),mmyp(ngr),mmzp(ngr),nkr &
                                 ,scr,vtab_r(nv,ngr)%var_p)
               case(10)
                  CALL unarrange (nkppz,mmxp(ngr),mmyp(ngr) &
                                 ,scr,vtab_r(nv,ngr)%var_p)
            end select
            cycle varloop
         endif
      enddo
   enddo varloop
   
   CALL shdf5_close (h5_fid)

enddo

deallocate(scr,anal_table)

return
END SUBROUTINE hist_read
