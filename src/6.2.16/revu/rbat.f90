!##############################################################################
Subroutine plotspc (nfl,action)

use an_header
use rhdf5_utils
use rcommons
use mem_grid
use node_mod, only:nmachs
use hdf5_utils

implicit none

character(len=*) :: action
character(len=1) :: toksepfr
character(len=2) :: cgrid
character(len=5) :: ctime
character(len=64) :: frtokens(50),cvar
character(len=strl1) :: flnm
character(len=RHDF5_MAX_STRING) :: rh5_file_name,rh5_file_acc
integer :: lastslash,ngds,ngd,ivtime,ibegs,iends,nplot,ntokfr &
  ,nfl,nfile,cvar_ok,var_ok(maxrevu)
integer*8:: rh5_file

!File and parallel info (should only do this sequentially, not parallel)
integer :: h5_fid, iphdf5

data toksepfr/'/'/  

!Not running in parallel so use default info for sequential
if (nmachs .gt. 1) then
  stop 'Should not be attempting to run REVU in parallel'
else
  iphdf5 = 0
endif

CALL rams_get_idata (1,1,ngds)

!Determine if we are outputting all grids for a single grid
if(igrid.gt.0 .and. igrid.le.ngds) then
   ngd=igrid
else
   print*,'The model grids you are trying to process does not exist'
   stop
endif

! set type of z coordinate for file name prefix
!  according to value in iztran
!    1 - sigma      --> 'S'
!    2 - cartesian  --> 'C'
!    3 - pressure   --> 'P'
!    4 - ground     --> 'G'
if(iztran.eq.1) then
   ftran='S'               !for sigma
elseif(iztran.eq.2) then
   ftran='C'               !for cartesian
elseif(iztran.eq.3) then
   ftran='P'               !for pressure level
elseif(iztran.eq.4) then
   ftran='G'               !for ground
endif

if(action(1:4)=='TEXT' .or. action(1:4)=='HDF5') then
   !Output some grid type and time information for HDF5 output
   if(action(1:4)=='HDF5') then
     ! open the output file
     CALL rams_get_cdata (0,1,flnm)
     write(rh5_file_name,'(2a,2a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6,a2,i0,a3)') &
       revpref(1:len_trim(revpref))  &
      ,flnm(lastslash(flnm)+1:len_trim(flnm)-27),ftran,'-'  &
      ,iyear1,'-',imonth1,'-',idate1,'-',itime1*100,'-g',ngd,'.h5'
     print*
     print*,'===='
     print*,'HDF5 file: ',rh5_file_name(1:len_trim(rh5_file_name))
     print*
     ! open in write mode    --> 2nd arg = 'W'
     ! delete file if exists --> 3rd arg = 1
     rh5_file_acc = 'W'
     CALL shdf5_open (rh5_file_name,rh5_file_acc,iphdf5,rh5_file,1)
   endif
endif
   
!LOOP OVER ALL REQUESTED TIMES UP TO MAX NUMBER OF TIMES AVAILABLE
ivtime=0
ibegs=max(itbeg,1)   !beginning time to process >=1
if(itend.eq.0)iends=nfl !last time to process
if(itend.gt.0)iends=min(itend,nfl) !last time to process
iends=max(1,iends)   !last time to process >=1
itstep=max(itstep,1) !output timestepping

do nfile=ibegs,iends,itstep
     write(*,'(a,4i5)') 'Doing file,grid - ',nfile,ngd
     ivtime=ivtime+1 !increment number of times processed

     write(cgrid,'(i2)') ngd
     write(ctime,'(i5)') nfile

     !loop through max number of plots
     do nplot=1,maxrevu

        var_ok(nplot)=0

        !only do this if have readable primary var
        if(revuvar(nplot)(1:1).ne.'/') goto 210
  
        frtokens(1)=' '
        CALL tokenize1 (revuvar(nplot),frtokens,ntokfr,toksepfr)
        if(frtokens(1)=='') frtokens(1)='none'
        if(frtokens(1)=='none') goto 210

        cvar=frtokens(1)

        ! print action summary
        print*
        if(cvar(1:4).ne.'none')  &
           print*,'Doing var- ',cvar(1:len(cvar))

        ! prohibit use of interpollated wind directions
        if(cvar(1:9)=='direction'.and.iztran.ne.1) then
           print*,'Cannot output wind direction at interpollated grid points'
           print*,'Skipping variable'
           goto 210
        endif

        CALL read_rams (action,ivtime,cgrid,ctime,iztran,cvar,rh5_file,cvar_ok)
        var_ok(nplot)=cvar_ok

        210 continue

     enddo !loop thru variables
enddo !loop thru times to process

!SEND GRID NAVIGATION AND COORDINATE INFORMATION TO HDF5 OUTPUT FILE
if(action(1:4)=='HDF5') then
  print*
  print*, 'Attaching dimension (coordinate) specs to variables in hdf5 file: '
  print*
  print*, '  HDF5 file: ', trim(rh5_file_name)
  print*

  ! declare the datasets that hold coordinate values for dimensions
  CALL rhdf5_set_var_to_dim (rh5_file,'x_coords','x')
  CALL rhdf5_set_var_to_dim (rh5_file,'y_coords','y')
  CALL rhdf5_set_var_to_dim (rh5_file,'z_coords','z')
  CALL rhdf5_set_var_to_dim (rh5_file,'t_coords','t')

  do nplot = 1,maxrevu
    ! only the variable specs start with a '/'
    ! if the first char is a '/' than get the variable name and attach the dimension
    if (revuvar(nplot)(1:1) .eq. '/') then
     frtokens(1)=' '
     CALL tokenize1 (revuvar(nplot),frtokens,ntokfr,toksepfr)
     if(frtokens(1) .ne. '' .and. var_ok(nplot)==1) then
        cvar=frtokens(1)
        print*, '  Doing var- ', trim(cvar)
        CALL rhdf5_attach_dims_to_var (rh5_file,cvar)
        print*
     endif
    endif
  enddo

  print*, 'Closing hdf5 file: ', trim(rh5_file_name)
  print*
  CALL shdf5_close (rh5_file)
  print*
endif

return
END SUBROUTINE plotspc

!##############################################################################
Subroutine backset ()

use rcommons

implicit none

CALL var_parse (xvar,ixbeg,ixend,ixstep)
CALL var_parse (yvar,iybeg,iyend,iystep)
CALL var_parse (zvar,izbeg,izend,izstep)
CALL var_parse (tvar,itbeg,itend,itstep)

ixstep=max(1,ixstep)
iystep=max(1,iystep)
izstep=max(1,izstep)

return
END SUBROUTINE backset

!##############################################################################
Subroutine var_parse (string,isbeg,isend,istep)

implicit none

! This routine takes the strings in XVAR,YVAR,ZVAR and TVAR and
! parses them out into more usable components which are used throughout.

integer, external :: numberchk
integer :: ntok,isbeg,isend,istep
character(len=*) :: string
character(len=8) :: tokens(50)
character(len=1) :: toksep(3)
data toksep/'/',':',','/

CALL tokenize (string,tokens,ntok,toksep,3)

! Start parsing the range information, which is delimeted
! by colons and must start at the 2nd token

CALL ch2int (tokens(2),isbeg)
if(ntok.gt.4 .and. numberchk(tokens(4))==1) then
   CALL ch2int (tokens(4),isend)
else
   isend=isbeg
endif

if(ntok.gt.6 .and. numberchk(tokens(6))==1) then
   CALL ch2int (tokens(6),istep)
else
   istep=1
endif

return
END SUBROUTINE var_parse
