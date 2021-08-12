!##############################################################################
Module hdf5_utils

implicit none

integer, parameter :: HDF5_MAX_DIMS = 10

type hdf5_select_type
  integer :: ndims
  ! description of memory space
  integer, dimension(HDF5_MAX_DIMS) :: dims, block, count, stride, offset
endtype hdf5_select_type

Contains

!##############################################################################
Subroutine shdf5_open (locfn,access,iphdf5,fileid,idelete)

! This routine will open an HDF5 file specified by the fname, facc, fdelete
! arguments, and pass back the file id number in the fid argument.
!
!  facc: 'R'  --> read
!        'RW' --> read/write
!        'W'  --> write
!
!  fdelete: 0 --> do not delete if file exists
!           1 --> delete if opening in write mode
!
! If read mode
!    file exists --> open (fh5f_open)
!    file does not exist --> error
!
! If read/write mode
!    file exists --> open (fh5f_open)
!    file does not exist --> error
!
! If write mode
!    file exists:
!      fdelete flag is 1 --> truncate file (fh5f_create)
!      fdelete fla0 is 0 --> error
!    file does not exists --> create (fh5f_create)
!
! For the 'iaccess' argument to fh5f_open and fh5f_create routines
! (iaccess is 2nd argument to both routines)
!   If calling fh5f_open:
!      read only  --> set iaccess to 1
!      read/write --> set iaccess to 2
!   If calling fh5f_create
!      truncate the file if file exists --> set iaccess to 1
!      fail if file exists              --> set iaccess to 2

use node_mod
        
implicit none

character(len=*) :: locfn     ! file name
character(len=*) :: access    ! File access ('R','W','RW')
integer :: iphdf5             ! 1 -> use PHDF5, 0 -> use HDF5
integer*8 :: fileid             ! file id number from HDF5
integer, optional :: idelete  ! If W, delete/overwrite file if exists?
                              ! 1=yes, 0=no
                              ! Only needed when access='W'

integer :: hdferr ! Error flag
integer :: iaccess ! int access flag
character(len=2) :: caccess ! File access ('R ','W ','RW')

logical :: exists ! File existence

caccess = access

! Check for existence of RAMS file.

inquire(file=trim(locfn),exist=exists)

! Create a new file or open an existing RAMS file.
if (access(1:1) == 'R') then
   if (.not.exists) then
      print*,'shdf5_open:'
      print*,'   Attempt to open a file for reading that does not exist.'
      print*,'   Filename: ',trim(locfn)
      stop 'shdf5_open: no file'
   else
      if (caccess == 'R ') iaccess = 1
      if (caccess == 'RW') iaccess = 2

      ! this odd construction of //char(0) is because C requires strings
      ! to be null-terminated- this appends a null char at the end.
      CALL fh5f_open (trim(locfn)//char(0), iaccess, iphdf5, fileid, hdferr)
      
      if (hdferr < 0) then
         print*,'shdf5_open:'
         print*,'   Error opening hdf5 file - error -',hdferr
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open: open error'      
      endif
   endif
elseif (access(1:1) == 'W') then
   if (.not.exists) then
      iaccess=2
      CALL fh5f_create (trim(locfn)//char(0), iaccess, iphdf5, fileid, hdferr)
   else
      ! File exists and we want to create a new file
      ! With multiple processes, we have to be careful not to remove the
      ! file before every node has had a chance to see it existing. We also
      ! have to be careful to not allow anyone to call fh5f_create() until
      ! the file is removed so that all nodes have opened the same file
      ! (in a Linux inode sense).
      !
      ! To accomplish this, only allow mainnum to delete the file while calling
      ! par_pause() before and after the delete command. This defers the delete
      ! until everyone has had a chance to see the file and note its existence.
      ! It also forces everyone to wait until the file is gone before trying
      ! to open the new one.

      ! Check for proper idelete setting
      if(.not.present(idelete) ) then
         print*,'shdf5_open: idelete not specified when access=W'
         stop 'shdf5_open: no idelete'
      endif
      
      if(idelete == 0) then
         print*,'In shdf5_open:'
         print*,'   Attempt to open an existing file for writing, '
         print*,'      but overwrite is disabled. idelete=',idelete
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open'
      endif

      if (nmachs .gt. 1) then
        CALL par_pause (my_rams_num, 101)
        ! Only mainnum removes the file
        if (my_rams_num .eq. mainnum) then
           CALL usystem ('rm -f '//trim(locfn)//char(0))
        endif
        CALL par_pause (my_rams_num, 102)
      else
        CALL usystem ('rm -f '//trim(locfn)//char(0))
      endif

      iaccess=1
      CALL fh5f_create (trim(locfn)//char(0), iaccess, iphdf5, fileid, hdferr)
   endif
   if(hdferr < 0) then
      print*,'HDF5 file create failed:',hdferr
      print*,'file name:',trim(locfn),' ',trim(access), idelete
      stop 'shdf5_open: bad create'
   endif
endif

return
END SUBROUTINE shdf5_open

!##############################################################################
Subroutine shdf5_info (fileid,dsetname,ndims,dims)

implicit none

integer*8 :: fileid
character(len=*) :: dsetname ! Dataset name
integer :: dims(*)
integer :: ndims ! Dataset rank (in file)
integer, dimension(6) :: dimsh
integer :: hdferr ! Error flag

! Reading the dimension sizes from the C interface will return them in the
! reverse order that FORTRAN needs them (see notes in shdf5_orec).
CALL fh5d_info (fileid, trim(dsetname)//char(0), ndims, dimsh, hdferr)
if (hdferr < 0) then
   print*,'shdf5_info: call fh5d_info:',trim(dsetname),'hdf5 error =',hdferr
   stop   'shdf5_info'
endif

CALL shdf5_reverse_array (ndims, dimsh, dims)

return
END SUBROUTINE shdf5_info

!##############################################################################
Subroutine shdf5_orec (fileid,iphdf5,dsetname,mem_select,file_select &
                      ,file_chunks &
                      ,ivara,rvara,cvara,dvara,lvara  &
                      ,ivars,rvars,cvars,dvars,lvars)

implicit none

character(len=*) :: dsetname ! Variable label
type (hdf5_select_type) :: mem_select, file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks
integer*8 :: fileid
integer :: iphdf5

! Array and scalar arguments for diff types. Only specify one in each call.
integer, optional :: ivara(*),ivars
real,    optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
double precision,    optional :: dvara(*),dvars
logical,   optional :: lvara(*),lvars

integer:: h5_type   ! Local type designator

integer, dimension(:), allocatable :: m_sel, f_sel
integer :: m_ndims, f_ndims
integer, dimension(HDF5_MAX_DIMS) :: f_chunks

character(len=2) :: ctype    ! Variable type: int, real, char
integer :: hdferr ! Error flag

! Find which data type is input
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_orec'
   stop 'shdf5_orec: bad data field'
endif

if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4
if (ctype(1:1) == 'l') h5_type=5

! Set up hyperslab selection specs for memory and file spaces
!
! Since we are working in FORTRAN and switching to C for the HDF5 I/O,
! we need to reverse the specification of dimension sizes. This is 
! necessary to keep the data in the HDF5 file consistent. This need 
! comes about due to the column-major array storage order in FORTRAN
! versus the row-major storage in C. Eg., 2D array in FORTRAN that 
! has 3 rows and 2 columns:
!    1 4
!    2 5
!    3 6
!
! It's typical to associate i with the rows and j with the columns,
! ie the array is indexed using (i,j).
!
! FORTRAN stores this in linear memory as:
!    1 2 3 4 5 6
!
! If we tell C that this linear storage represents an array that is
! 3 rows and 2 columns, C will see this as:
!    1 2
!    3 4
!    5 6
!
! Note that the association with (i,j) is now broken. Eg, A(i,j) where
! i = 1,2,3 and j is held at 1 gives 1 3 5, which in FORTRAN would have 
! given 1 2 3.
!
! However, if we reverse the specs of the dimension sizes (tell C that 
! the array is 2 rows by 3 columns), then C will see this as:
!    1 2 3
!    4 5 6
! Along with reversing the dimension sizes we also reverse the association 
! with i and j: now the rows are associated with j and the columns are 
! associated with i. With this scheme, when C is asked to give A(j,i) 
! where i - 1,2,3 and j is held at 1, this results in 1 2 3 which matches 
! what you would get in FORTRAN.
!
! Note: This could be also be accomplished by transposing the data in 
! FORTRAN right before handing it to C, but it is way more efficient
! to reverse the dimension specs.

CALL shdf5_build_select (mem_select, m_ndims, m_sel)
CALL shdf5_build_select (file_select, f_ndims, f_sel)
CALL shdf5_reverse_array (f_ndims, file_chunks, f_chunks)

! Write the dataset.
if (ctype == 'is') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,ivars,hdferr)
elseif (ctype == 'rs') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,rvars,hdferr)
elseif (ctype == 'cs') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,cvars,hdferr)
elseif (ctype == 'ds') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,dvars,hdferr)
elseif (ctype == 'ls') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,lvars,hdferr)
elseif (ctype == 'ia') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,ivara,hdferr)
elseif (ctype == 'ra') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,rvara,hdferr)
elseif (ctype == 'ca') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,cvara,hdferr)
elseif (ctype == 'da') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,dvara,hdferr)
elseif (ctype == 'la') then
   CALL fh5d_write (fileid,trim(dsetname)//char(0),h5_type,iphdf5 &
                   ,m_ndims,m_sel,f_ndims,f_sel,f_chunks,lvara,hdferr)
endif

if (hdferr /= 0) then
   print*,'In shdf5_orec: hdf5 write error =',hdferr
   stop 'shdf5_orec: hdf5 write error'
endif

! Deallocate arrays allocated by shdf5_build_select()
deallocate(m_sel,f_sel)

return
END SUBROUTINE shdf5_orec

!##############################################################################
Subroutine shdf5_irec (fileid,iphdf5,dsetname,mem_select,file_select &
                      ,ivara,rvara,cvara,dvara,lvara  &
                      ,ivars,rvars,cvars,dvars,lvars)
        
implicit none

integer*8 :: fileid
integer :: iphdf5
character(len=*) :: dsetname ! Dataset name
type (hdf5_select_type) :: mem_select, file_select

! Array and scalar arguments for diff types. Only specify one in each call.
integer, optional :: ivara(*),ivars
real,    optional :: rvara(*),rvars
character(len=*), optional :: cvara(*),cvars
double precision,    optional :: dvara(*),dvars
logical,   optional :: lvara(*),lvars

integer, dimension(:), allocatable :: m_sel, f_sel
integer :: m_ndims, f_ndims

integer:: h5_type   ! Local type designator

integer :: hdferr ! Error flag

character(len=2) :: ctype

! Find which data type will be read
    if(present(ivars)) then ; ctype='is'
elseif(present(rvars)) then ; ctype='rs'
elseif(present(cvars)) then ; ctype='cs'
elseif(present(dvars)) then ; ctype='ds'
elseif(present(lvars)) then ; ctype='ls'
elseif(present(ivara)) then ; ctype='ia'
elseif(present(rvara)) then ; ctype='ra'
elseif(present(cvara)) then ; ctype='ca'
elseif(present(dvara)) then ; ctype='da'
elseif(present(lvara)) then ; ctype='la'
else
   print*,'Incorrect or missing data field argument in shdf5_irec'
   stop 'shdf5_irec: bad data field'
endif

! Figure out the file hyperslab selection
    
! Read data from hyperslab in the file into the memory buffer
if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4
if (ctype(1:1) == 'l') h5_type=5

! Set up hyperslab selection specs for memory and file spaces
CALL shdf5_build_select (mem_select,m_ndims,m_sel)
CALL shdf5_build_select (file_select,f_ndims,f_sel)

if (ctype == 'is') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,ivars,hdferr)
elseif (ctype == 'rs') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,rvars,hdferr)
elseif (ctype == 'cs') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,cvars,hdferr)
elseif (ctype == 'ds') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,dvars,hdferr)
elseif (ctype == 'ls') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,lvars,hdferr)
elseif (ctype == 'ia') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,ivara,hdferr)
elseif (ctype == 'ra') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,rvara,hdferr)
elseif (ctype == 'ca') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,cvara,hdferr)
elseif (ctype == 'da') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,dvara,hdferr)
elseif (ctype == 'la') then
   CALL fh5d_read (fileid,trim(dsetname)//char(0),h5_type &
                  ,m_ndims,m_sel,f_ndims,f_sel,lvara,hdferr)
endif

if (hdferr /= 0) then
   print*,'shdf5_irec: call fh5d_read: hdf5 error =',hdferr
   stop
endif

! Deallocate arrays allocated by shdf5_build_select()
deallocate(m_sel,f_sel)

return
END SUBROUTINE shdf5_irec

!##############################################################################
Subroutine shdf5_close (fileid)
        
implicit none

integer*8 :: fileid

integer :: hdferr  ! Error flags

! Close RAMS hdf file.

CALL fh5f_close (fileid, hdferr)

return
END SUBROUTINE shdf5_close

!#############################################################################
Subroutine shdf5_reverse_array (n,array,rev_array)

implicit none

integer :: n
integer, dimension(n) :: array, rev_array
integer :: i, irev

!Array ordering section
do i = 1,n
  irev = (n-i) + 1
  rev_array(irev) = array(i)
enddo

!We used array swap starting in RAMS 6.1.12. To process old data with REVU
!just comment out the do-loop above and uncomment the line below for rev_array
!Saleeby(2018).
!rev_array(:) = array(:)

return
END SUBROUTINE shdf5_reverse_array

!#############################################################################
! This routine will build a linear array with the hdf5 hyperslab selection
! specs loaded into it. To facilitate passing arguments to fh5d_read and
! fh5d_write, a linear array containing the dimension sizes, and hyperslab
! block, count, offset, and stride values is created. The order of loading
! these data into the array is:
!
!     dims
!     block
!     count
!     offset
!     stride
!
! Note that this order needs to be kept in sync with the extraction of these
! data in fh5d_read and fh5d_write.
!
! Also, since we are going from FORTRAN to C we need to reverse the order
! of the entries we find in these arrays (dims, block, count, etc.)
!
! It is assumed that the caller will deallocate select_array

Subroutine shdf5_build_select (hs_select,ndims,select_array)

implicit none

  type (hdf5_select_type) :: hs_select
  integer :: ndims
  integer, dimension(:), allocatable :: select_array

  integer :: i, isa
  integer, dimension(HDF5_MAX_DIMS) :: temp_array

  ! Need to allocate 5 * ndims, where 5 is the number of quantities
  ! we are loading from the hs_select struct.
  allocate(select_array(5 * hs_select%ndims))

  ! Go in the proper order: dims, block, count, offset, stride
  isa = 1
  ndims = hs_select%ndims

  ! 1: dims
  CALL shdf5_reverse_array (ndims,hs_select%dims,temp_array)
  do i = 1, ndims
    select_array(isa) = temp_array(i)
    isa = isa + 1
  enddo
  
  ! 2: block
  CALL shdf5_reverse_array (ndims,hs_select%block,temp_array)
  do i = 1, ndims
    select_array(isa) = temp_array(i)
    isa = isa + 1
  enddo

  ! 3: count
  CALL shdf5_reverse_array (ndims,hs_select%count,temp_array)
  do i = 1, ndims
    select_array(isa) = temp_array(i)
    isa = isa + 1
  enddo

  ! 4: offset
  CALL shdf5_reverse_array (ndims,hs_select%offset,temp_array)
  do i = 1, ndims
    select_array(isa) = temp_array(i)
    isa = isa + 1
  enddo

  ! 5: stride
  CALL shdf5_reverse_array (ndims,hs_select%stride,temp_array)
  do i = 1, ndims
    select_array(isa) = temp_array(i)
    isa = isa + 1
  enddo

return
END SUBROUTINE shdf5_build_select

!######################################################################
! This routine determines the hyperslab selection specs for the given io
! transaction.
!
! io_type:
!   'R' -> read, select entire sub-domain (ie include all overlap)
!   'W' -> write, select entire sub-domain minus the internal overlap
!
! ivar_type:
!    negative -> Vector variable (i), i is abs(ivar_type)
!       Note a scalar variable (single value) can be obtained by
!       setting ivar_type to either +1 or -1. Allow the +1 value
!       for backward compatibility.
!
! 2 -> Atmospheric, ISAN 2D variable (i,j)
! 3 -> Atmospheric 3D variable (i,j,k), k is nnzp (atmos. levels)
! 4 -> Leaf3 4D soil variable (i,j,k,p), k is nzg (soil levels), p is patch
! 5 -> Leaf3 4D snow variable (i,j,k,p), k is nzs (snow levels), p is patch
! 6 -> Leaf3 3D variable (i,j,p), p is patch
! 7 -> BIN 4D atmospheric variable (i,j,k,kr), kr is nkr (number of bins)
! 8 -> ISAN 3D isentropic variable (i,j,k), k is nisn (isentropic levels)
! 9 -> ISAN 3D sigma-z variable (i,j,k), k is nsigz (sigma-z levels)
!
! This routine requires that it be called in between any 
! arrange/rearrage/unarrange calls and the shdf5_irec/shdf5_orec calls.
! It only needs to be called once before shdf5_irec/shdf5_orec calls
! as long as the variables being read/written are the same dimensions.
! If you change hard-coded dimensions from say a 2D to 3D variables,
! you need to re-call this routine.

Subroutine shdf5_set_hs_select (ivar_type,io_type,igrid,m_select &
                               ,f_select,f_chunks)

use node_mod
use mem_grid
use isan_coms
use micro_prm, only:nkr
use kpp_parameters, only:nkppz

implicit none

  integer :: ivar_type, igrid
  character (len=1) :: io_type
  type (hdf5_select_type) :: m_select, f_select
  integer, dimension(HDF5_MAX_DIMS) :: f_chunks

  integer :: m_xblk, m_yblk, m_xoff, m_yoff
  integer :: f_xblk, f_yblk, f_xoff, f_yoff
  integer :: iv_type, vlen

  ! If ivar_type is negative, it is signaling that a vector (1D) variable
  ! is being requested. Eg. the coordinate values of one of the dimensions
  ! of a field. The length is the absolute value of ivar_type, and then
  ! iv_type will be set to 1 to signal to subsequent code that a
  ! vector variable is being accessed.
  vlen = abs(ivar_type)
  if (ivar_type .lt. 0) then
    iv_type = 1
  else
    iv_type = ivar_type
  endif

  if (trim(ramsorrevu) == 'REVU') then
   mmxp(igrid)       = nnxp(igrid)
   mmyp(igrid)       = nnyp(igrid)
   mmzp(igrid)       = nnzp(igrid)
   mem_read(igrid)%xblock = nnxp(igrid)
   mem_read(igrid)%yblock = nnyp(igrid)
   mem_read(igrid)%xoff   = 0
   mem_read(igrid)%yoff   = 0
   mem_write(igrid)%xblock = nnxp(igrid)
   mem_write(igrid)%yblock = nnyp(igrid)
   mem_write(igrid)%xoff   = 0
   mem_write(igrid)%yoff   = 0
   file_read(igrid)%xblock = nnxp(igrid)
   file_read(igrid)%yblock = nnyp(igrid)
   file_read(igrid)%xoff   = 0
   file_read(igrid)%yoff   = 0
   file_write(igrid)%xblock = nnxp(igrid)
   file_write(igrid)%yblock = nnyp(igrid)
   file_write(igrid)%xoff   = 0
   file_write(igrid)%yoff   = 0
   file_xchunk(igrid) = nnxp(igrid)
   file_ychunk(igrid) = nnyp(igrid)
  endif

  ! First determine impact of read vs write
  if (io_type .eq. 'R') then
    ! memory select
    m_xblk = mem_read(igrid)%xblock
    m_yblk = mem_read(igrid)%yblock
    m_xoff = mem_read(igrid)%xoff
    m_yoff = mem_read(igrid)%yoff
    ! file select
    f_xblk = file_read(igrid)%xblock
    f_yblk = file_read(igrid)%yblock
    f_xoff = file_read(igrid)%xoff
    f_yoff = file_read(igrid)%yoff
  elseif (io_type .eq. 'W') then
    m_xblk = mem_write(igrid)%xblock
    m_yblk = mem_write(igrid)%yblock
    m_xoff = mem_write(igrid)%xoff
    m_yoff = mem_write(igrid)%yoff
    f_xblk = file_write(igrid)%xblock
    f_yblk = file_write(igrid)%yblock
    f_xoff = file_write(igrid)%xoff
    f_yoff = file_write(igrid)%yoff
  else
    print*, 'ERROR: shdf5_set_hs_select: Unknown io type: ', trim(io_type)
    stop
  endif

  ! Set selection according to variable type
  select case (iv_type)
    case (1)
      ! Vector variable
      ! memory
      m_select%ndims = 1
      m_select%dims(1) = vlen
      m_select%block(1) = vlen
      m_select%count(1) = 1
      m_select%offset(1) = 0
      m_select%stride(1) = 1
  
      ! file
      f_select%ndims = 1
      f_select%dims(1) = vlen
      f_select%block(1) = vlen
      f_select%count(1) = 1
      f_select%offset(1) = 0
      f_select%stride(1) = 1
  
      f_chunks(1) = vlen
    case (2)
      ! Atmos, ISAN 2D var
      ! memory
      m_select%ndims = 2
      m_select%dims(1:2) = (/ mmxp(igrid), mmyp(igrid) /)
      m_select%block(1:2) = (/ m_xblk, m_yblk /)
      m_select%count(1:2) = (/ 1, 1 /)
      m_select%offset(1:2) = (/ m_xoff, m_yoff /)
      m_select%stride(1:2) = (/ 1, 1 /)
  
      ! file
      f_select%ndims = 2
      f_select%dims(1:2) = (/ nnxp(igrid), nnyp(igrid) /)
      f_select%block(1:2) = (/ f_xblk, f_yblk /)
      f_select%count(1:2) = (/ 1, 1 /)
      f_select%offset(1:2) = (/ f_xoff, f_yoff /)
      f_select%stride(1:2) = (/ 1, 1 /)
  
      f_chunks(1:2) = (/ file_xchunk(igrid), file_ychunk(igrid) /)
    case (3)
      ! Atmos 3D var
      ! memory
      m_select%ndims = 3
      m_select%dims(1:3) = (/ mmxp(igrid), mmyp(igrid), mmzp(igrid) /)
      m_select%block(1:3) = (/ m_xblk, m_yblk, mmzp(igrid) /)
      m_select%count(1:3) = (/ 1, 1, 1 /)
      m_select%offset(1:3) = (/ m_xoff, m_yoff, 0 /)
      m_select%stride(1:3) = (/ 1, 1, 1 /)
  
      ! file
      f_select%ndims = 3
      f_select%dims(1:3) = (/ nnxp(igrid), nnyp(igrid), nnzp(igrid) /)
      f_select%block(1:3) = (/ f_xblk, f_yblk, mmzp(igrid) /)
      f_select%count(1:3) = (/ 1, 1, 1 /)
      f_select%offset(1:3) = (/ f_xoff, f_yoff, 0 /)
      f_select%stride(1:3) = (/ 1, 1, 1 /)
  
      f_chunks(1:3) = (/ file_xchunk(igrid), file_ychunk(igrid), mmzp(igrid) /)
    case (4)
      ! Leaf3 4D soil var
      ! memory
      m_select%ndims = 4
      m_select%dims(1:4) = (/ mmxp(igrid), mmyp(igrid), nzg, npatch /)
      m_select%block(1:4) = (/ m_xblk, m_yblk, nzg, npatch /)
      m_select%count(1:4) = (/ 1, 1, 1, 1 /)
      m_select%offset(1:4) = (/ m_xoff, m_yoff, 0, 0 /)
      m_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      ! file
      f_select%ndims = 4
      f_select%dims(1:4) = (/ nnxp(igrid), nnyp(igrid), nzg, npatch /)
      f_select%block(1:4) = (/ f_xblk, f_yblk, nzg, npatch /)
      f_select%count(1:4) = (/ 1, 1, 1, 1 /)
      f_select%offset(1:4) = (/ f_xoff, f_yoff, 0, 0 /)
      f_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      f_chunks(1:4) = (/ file_xchunk(igrid), file_ychunk(igrid), nzg, npatch /)
    case (5)
      ! Leaf3 4D surface water var
      ! memory
      m_select%ndims = 4
      m_select%dims(1:4) = (/ mmxp(igrid), mmyp(igrid), nzs, npatch /)
      m_select%block(1:4) = (/ m_xblk, m_yblk, nzs, npatch /)
      m_select%count(1:4) = (/ 1, 1, 1, 1 /)
      m_select%offset(1:4) = (/ m_xoff, m_yoff, 0, 0 /)
      m_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      ! file
      f_select%ndims = 4
      f_select%dims(1:4) = (/ nnxp(igrid), nnyp(igrid), nzs, npatch /)
      f_select%block(1:4) = (/ f_xblk, f_yblk, nzs, npatch /)
      f_select%count(1:4) = (/ 1, 1, 1, 1 /)
      f_select%offset(1:4) = (/ f_xoff, f_yoff, 0, 0 /)
      f_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      f_chunks(1:4) = (/ file_xchunk(igrid), file_ychunk(igrid), nzs, npatch /)
    case (6)
      ! Leaf3 3D var
      ! memory
      m_select%ndims = 3
      m_select%dims(1:3) = (/ mmxp(igrid), mmyp(igrid), npatch /)
      m_select%block(1:3) = (/ m_xblk, m_yblk, npatch /)
      m_select%count(1:3) = (/ 1, 1, 1 /)
      m_select%offset(1:3) = (/ m_xoff, m_yoff, 0 /)
      m_select%stride(1:3) = (/ 1, 1, 1 /)
  
      ! file
      f_select%ndims = 3
      f_select%dims(1:3) = (/ nnxp(igrid), nnyp(igrid), npatch /)
      f_select%block(1:3) = (/ f_xblk, f_yblk, npatch /)
      f_select%count(1:3) = (/ 1, 1, 1 /)
      f_select%offset(1:3) = (/ f_xoff, f_yoff, 0 /)
      f_select%stride(1:3) = (/ 1, 1, 1 /)
  
      f_chunks(1:3) = (/ file_xchunk(igrid), file_ychunk(igrid), npatch /)
    case (7)
      ! Bin micro 4D var
      ! memory
      m_select%ndims = 4
      m_select%dims(1:4) = (/ mmxp(igrid), mmyp(igrid), mmzp(igrid), nkr /)
      m_select%block(1:4) = (/ m_xblk, m_yblk, mmzp(igrid), nkr /)
      m_select%count(1:4) = (/ 1, 1, 1, 1 /)
      m_select%offset(1:4) = (/ m_xoff, m_yoff, 0, 0 /)
      m_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      ! file
      f_select%ndims = 4
      f_select%dims(1:4) = (/ nnxp(igrid), nnyp(igrid), nnzp(igrid), nkr /)
      f_select%block(1:4) = (/ f_xblk, f_yblk, mmzp(igrid), nkr /)
      f_select%count(1:4) = (/ 1, 1, 1, 1 /)
      f_select%offset(1:4) = (/ f_xoff, f_yoff, 0, 0 /)
      f_select%stride(1:4) = (/ 1, 1, 1, 1 /)
  
      f_chunks(1:4) = (/ file_xchunk(igrid), file_ychunk(igrid), mmzp(igrid), nkr /)
    case (8)
      ! ISAN 3D isentropic var
      ! memory
      m_select%ndims = 3
      m_select%dims(1:3) = (/ mmxp(igrid), mmyp(igrid), nisn /)
      m_select%block(1:3) = (/ m_xblk, m_yblk, nisn /)
      m_select%count(1:3) = (/ 1, 1, 1 /)
      m_select%offset(1:3) = (/ m_xoff, m_yoff, 0 /)
      m_select%stride(1:3) = (/ 1, 1, 1 /)
  
      ! file
      f_select%ndims = 3
      f_select%dims(1:3) = (/ nnxp(igrid), nnyp(igrid), nisn /)
      f_select%block(1:3) = (/ f_xblk, f_yblk, nisn /)
      f_select%count(1:3) = (/ 1, 1, 1 /)
      f_select%offset(1:3) = (/ f_xoff, f_yoff, 0 /)
      f_select%stride(1:3) = (/ 1, 1, 1 /)
  
      f_chunks(1:3) = (/ file_xchunk(igrid), file_ychunk(igrid), nisn /)
    case(9)
      ! ISAN 3D sigma-z var
      ! memory
      m_select%ndims = 3
      m_select%dims(1:3) = (/ mmxp(igrid), mmyp(igrid), nsigz /)
      m_select%block(1:3) = (/ m_xblk, m_yblk, nsigz /)
      m_select%count(1:3) = (/ 1, 1, 1 /)
      m_select%offset(1:3) = (/ m_xoff, m_yoff, 0 /)
      m_select%stride(1:3) = (/ 1, 1, 1 /)
  
      ! file
      f_select%ndims = 3
      f_select%dims(1:3) = (/ nnxp(igrid), nnyp(igrid), nsigz /)
      f_select%block(1:3) = (/ f_xblk, f_yblk, nsigz /)
      f_select%count(1:3) = (/ 1, 1, 1 /)
      f_select%offset(1:3) = (/ f_xoff, f_yoff, 0 /)
      f_select%stride(1:3) = (/ 1, 1, 1 /)
  
      f_chunks(1:3) = (/ file_xchunk(igrid), file_ychunk(igrid), nsigz /)
    case (10)
      ! Ocean 3D var
      ! memory
      m_select%ndims = 3
      m_select%dims(1:3) = (/ mmxp(igrid), mmyp(igrid), nkppz /)
      m_select%block(1:3) = (/ m_xblk, m_yblk, nkppz /)
      m_select%count(1:3) = (/ 1, 1, 1 /)
      m_select%offset(1:3) = (/ m_xoff, m_yoff, 0 /)
      m_select%stride(1:3) = (/ 1, 1, 1 /)

      ! file
      f_select%ndims = 3
      f_select%dims(1:3) = (/ nnxp(igrid), nnyp(igrid), nkppz /)
      f_select%block(1:3) = (/ f_xblk, f_yblk, nkppz /)
      f_select%count(1:3) = (/ 1, 1, 1 /)
      f_select%offset(1:3) = (/ f_xoff, f_yoff, 0 /)
      f_select%stride(1:3) = (/ 1, 1, 1 /)

      f_chunks(1:3) = (/ file_xchunk(igrid), file_ychunk(igrid), nkppz /)
    case default
      print*, 'ERROR: shdf5_set_hs_select: Unknown variable type &
        (iv_type): ', iv_type
      stop
  endselect

return
END SUBROUTINE shdf5_set_hs_select

END MODULE hdf5_utils
