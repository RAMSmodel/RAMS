!##############################################################################
Module rhdf5_utils

implicit none

! max limits for arrays, strings, etc, keep these in sync with like named
! defines in rhdf5_f2c.c
integer, parameter :: RHDF5_MAX_STRING = 128
integer, parameter :: RHDF5_MAX_DIMS   =  10

! integer coding for HDF5 types, keep these in sync with like named defines
! in rhdf5_f2c.c
integer, parameter :: RHDF5_TYPE_STRING  = 0
integer, parameter :: RHDF5_TYPE_INTEGER = 1
integer, parameter :: RHDF5_TYPE_FLOAT   = 2
integer, parameter :: RHDF5_TYPE_CHAR    = 3

! structure to hold variable information
type Rhdf5Var
  character(len=RHDF5_MAX_STRING) :: vname
  integer :: ndims
  integer, dimension(RHDF5_MAX_DIMS) :: dims
  character(len=RHDF5_MAX_STRING) :: units
  character(len=RHDF5_MAX_STRING) :: descrip
  character(len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames
  real, dimension(:), allocatable :: vdata
end type Rhdf5Var

Contains

!******************************************************************************
! Routines for HDF5 REVU IO
!******************************************************************************

!******************************************************************************
! Low level routines
!******************************************************************************

!******************************************************************************
!rhdf5_write_variable type routines

! This routine will create and populate the dataset and attributes
! for an hdf5 file variable.
!
! We are going through the C interface since the HDF5 Fortran interface does
! not work with the pgf90 compiler. This interface will store the array data
! in row major fashion (first dimension changes the slowest when deteriming the 
! linear storage order for the array data) whereas Fortran stores array data
! in column major fashion (first dimension changes the fastest).
!
! The data coming into this routine is ordered (x,y,z) and is in column 
! major fashion.
!
! This routine will handle two styles of variable storage.
!   1. Store the variable once as is. This is intended for grid information 
!      such as lat/lon/levels values. It is not necessary to store the same 
!      thing over for each time step, just once when setting up the grid for 
!      the first time. Setting the itstep argument to zero denotes this type 
!      of variable. This variable will be stored with fixed dimensions in the 
!      HDF5 file.
!
!   2. Store the variable each time step. This is intended for RAMS variable 
!      (2D or 3D field) storage where for each time step a version of the 
!      field is stored. The storage will be optimized toward reading out the 
!      variable in a time step by time step fashion. This routine will allow 
!      for an unlimited number of time steps as it will create the HDF5 
!      dataset with the time dimension having an unlimited max dimension. 
!      Along with this HDF5 chunking will be used so that we don't have
!      to preallocate contiguous file space for all time steps before we start
!      writing in data. Making the chunk size match the size of the data coming
!      into this routine (x,y,z) will allow for writing each time step's data
!      into one chunk (which is a contiguous file space). Because of this 
!      scheme, we want time to be the slowest changing dimension which makes 
!      it so that for each time step only one chunk needs to be read out of 
!      the file to get the (x,y,z) data. Note that if time was not the slowest
!      changing dimension, then multiple chunks would have to read out of the 
!      file in order to fill in one time step's (x,y,z) data.
!
! Since we want to make t the slowest changing dimension (first dimension due 
! to the C interface with the row major array storage) just tagging on the 
! (x,y,z) data that came into this routine would make for an awkward storage 
! in the HDF5 file, where you are taking memory containing (x,y,z) in column 
! major order and writing that out to the HDF5 file as if it were row major 
! order.
!
! MATLAB and IDL use a clever trick to deal with this issue. Leave the linear 
! storage alone, and simply tell the C HDF5 write routine that the dimensions 
! for (x,y,z) go in the reverse order. This automatically changes the 
! perspective on the linear storage to row major and the data in the HDF5 file
! will be consistent. That is, the value at a given x,y,z location will remain
! at that location in the file making it possible for a C routine, using row 
! major order, to read in the data and have it match what was written out by 
! this routine.
!
! For example if the data coming in has (x,y,z) with dimension sizes (4,3,2) 
! then send the data as is to the C hdf5 write routine, but tell that routine 
! that the order of the data is now (z,y,x) with dimension sizes (2,3,4).
!
! If attaching time as a variable always place it first in the variable order.
! It is the slowest changing variable so placing it first is consistent with 
! row major ordering.
!
! Note that sdata is an array of strings where each string gets cast into a 
! C style string terminated with a null byte. cdata on the other hand gets 
! cast into a character array with a fixed length (RHDF5_MAX_STRING).
!******************************************************************************

!##############################################################################
Subroutine rhdf5_write_variable_rdata (id,vname,ndims,itstep,dims &
           ,units,descrip,dimnames,deflvl,rdata)

implicit none

  integer :: i,ndims,itstep,sizeof,hdferr,dtype,deflvl,dset_ndims &
            ,dsize
  integer*8 :: id, dsetid
  integer, dimension(ndims+1) :: dset_dims,ext_dims,chunk_dims
  integer :: idummy   ! to get the size of an integer
  real :: rdummy      ! to get the size of a real
  character :: cdummy ! to get the size of a character
  character(len=*) :: vname,units,descrip
  integer, dimension(ndims) :: dims
  character(len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
  character(len=RHDF5_MAX_STRING) :: dnstring,arrayorg

  real, dimension(*) :: rdata

  ! Check for valid argument values
  if (itstep .lt. 0) then
    print*,'ERROR: hdf5_write_variable: itstep argument must be >= zero'
    stop 'hdf5_write_variable: bad variable write'
  endif

  !For input type rdata
  dtype = RHDF5_TYPE_FLOAT
  dsize = sizeof(rdummy)

  ! If ndims is zero then the data coming in is time data. This is okay as
  ! long as itstep is not zero as well (which would mean that there is no 
  ! data to write). Print a warning and return if ndims and itstep are both
  ! zero.
  if ((ndims .eq. 0) .and. (itstep .eq. 0)) then
    print*,'WARNING: hdf5_write_variable:'
    print*,'Both ndims and itstep are zero,'
    print*,'skipping write for variable: ',trim(vname)
    return
  endif

  ! Array organization is always row major for now
  arrayorg = 'row major'

  ! Figure out what the "real" dimensions for the data will be. Takes into 
  ! account that the x,y,z dimensions need to be reversed for row major 
  ! ordering
  CALL rhdf5_build_dims_for_write (itstep,ndims,dims,dimnames &
                      ,dset_ndims,dset_dims,ext_dims,chunk_dims,dnstring)
  CALL rhdf5_adjust_chunk_sizes (dset_ndims,chunk_dims,dsize)

  CALL rh5d_setup_and_write (id,trim(vname)//char(0),dtype,0,dset_ndims &
       ,dset_dims,ext_dims,chunk_dims,deflvl,dsetid,rdata,hdferr)

  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_write_variable:'
    print*,'cannot write data for variable: ',trim(vname)
    stop 'hdf5_write_variable: bad variable write'
  endif

  ! write out the attributes
  CALL rh5a_write_anyscalar (dsetid,'units'//char(0) &
                       ,trim(units)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'long_name'//char(0) &
                       ,trim(descrip)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'ArrayOrg'//char(0) &
                       ,trim(arrayorg)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'DimNames'//char(0) &
                       ,trim(dnstring)//char(0),RHDF5_TYPE_STRING,hdferr)

  if (trim(vname) .eq. 'time') then
    ! helper attribute for ncdump
    CALL rh5a_write_anyscalar (dsetid,'C_format'//char(0),'%c'//char(0) &
                              ,RHDF5_TYPE_STRING,hdferr)
  endif

  ! close the dataset
  CALL rh5d_close (dsetid, hdferr)

return
END SUBROUTINE rhdf5_write_variable_rdata

!##############################################################################
Subroutine rhdf5_write_variable_fsngl (id,vname,ndims,itstep,dims &
           ,units,descrip,dimnames,deflvl,fsngl)

implicit none

  integer :: i,ndims,itstep,sizeof,hdferr,dtype,deflvl,dset_ndims &
            ,dsize
  integer*8 :: id, dsetid
  integer, dimension(ndims+1) :: dset_dims,ext_dims,chunk_dims
  integer :: idummy   ! to get the size of an integer
  real :: rdummy      ! to get the size of a real
  character :: cdummy ! to get the size of a character
  character(len=*) :: vname,units,descrip
  integer, dimension(ndims) :: dims
  character(len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
  character(len=RHDF5_MAX_STRING) :: dnstring,arrayorg

  real :: fsngl

  ! Check for valid argument values
  if (itstep .lt. 0) then
    print*,'ERROR: hdf5_write_variable: itstep argument must be >= zero'
    stop 'hdf5_write_variable: bad variable write'
  endif

  !For input type rdata
  dtype = RHDF5_TYPE_FLOAT
  dsize = sizeof(rdummy)

  ! If ndims is zero then the data coming in is time data. This is okay as
  ! long as itstep is not zero as well (which would mean that there is no 
  ! data to write). Print a warning and return if ndims and itstep are both
  ! zero.
  if ((ndims .eq. 0) .and. (itstep .eq. 0)) then
    print*,'WARNING: hdf5_write_variable:'
    print*,'Both ndims and itstep are zero,'
    print*,'skipping write for variable: ',trim(vname)
    return
  endif

  ! Array organization is always row major for now
  arrayorg = 'row major'

  ! Figure out what the "real" dimensions for the data will be. Takes into 
  ! account that the x,y,z dimensions need to be reversed for row major 
  ! ordering
  CALL rhdf5_build_dims_for_write (itstep,ndims,dims,dimnames &
                      ,dset_ndims,dset_dims,ext_dims,chunk_dims,dnstring)
  CALL rhdf5_adjust_chunk_sizes (dset_ndims,chunk_dims,dsize)

  CALL rh5d_setup_and_write (id,trim(vname)//char(0),dtype,0,dset_ndims &
       ,dset_dims,ext_dims,chunk_dims,deflvl,dsetid,fsngl,hdferr)

  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_write_variable:'
    print*,'cannot write data for variable: ',trim(vname)
    stop 'hdf5_write_variable: bad variable write'
  endif

  ! write out the attributes
  CALL rh5a_write_anyscalar (dsetid,'units'//char(0) &
                       ,trim(units)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'long_name'//char(0) &
                       ,trim(descrip)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'ArrayOrg'//char(0) &
                       ,trim(arrayorg)//char(0),RHDF5_TYPE_STRING,hdferr)
  CALL rh5a_write_anyscalar (dsetid,'DimNames'//char(0) &
                       ,trim(dnstring)//char(0),RHDF5_TYPE_STRING,hdferr)

  if (trim(vname) .eq. 'time') then
    ! helper attribute for ncdump
    CALL rh5a_write_anyscalar (dsetid,'C_format'//char(0),'%c'//char(0) &
                              ,RHDF5_TYPE_STRING,hdferr)
  endif

  ! close the dataset
  CALL rh5d_close (dsetid, hdferr)

return
END SUBROUTINE rhdf5_write_variable_fsngl

!##############################################################################
Subroutine rhdf5_build_dims_for_write (itstep,ndims,dims,dimnames &
                 ,dset_ndims,dset_dims,ext_dims,chunk_dims,dnstring)

! This routine will figure out the dimensions we want for the upcoming write.
!
! If itstep is > 0, then the time dimension will be tagged onto the front of
! dims creating one extra dimension to the data.
!
! ext_dims gets set to show which dimensions, if any, are extendable (only
! extending time dimension for now).
!
! chunk_dims shows how to organize chunking in the HDF5 file.

implicit none

integer :: itstep, ndims
integer, dimension(ndims) :: dims
character(len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
integer :: dset_ndims
integer, dimension(ndims+1) :: dset_dims, ext_dims, chunk_dims
character(len=RHDF5_MAX_STRING) :: dnstring

integer :: i, irev
integer :: dummy_int
real :: dummy_real
character :: dummy_char

! The array ext_dims holds a 1 (extendable) or 0 (not extendable) to 
! mark which, if any, dimensions of the variable are to be extendable. For 
! now we just want the time dimension, if it exists, to be extendable.
!
! The array chunk_dims holds a description of how big to make the chunks 
! in the HDF file. We want one chunk per time step. Always make the chunk 
! size match the size of the field (set the individual chunk sizes to the 
! size of the dimensions of the field). If the variable has a time dimension
! set the corresponding chunk_size for the time dimension (first dimension)
! to one which keep the overall chunk size matching the size of the field.
!
! Always place time as the first dimension. Then reverse the order of the 
! field dimensions which will automatically make the data appear to be row 
! major ordered.

  if (itstep .eq. 0) then
    dset_ndims = ndims
    do i = 1, ndims
      ! reverse the ordering
      irev = (ndims - i) + 1
      dset_dims(irev) = dims(i)
      ext_dims(irev) = 0
      chunk_dims(irev) = dims(i)
      if (i .eq. 1) then
        dnstring = trim(dimnames(irev))
      else
        dnstring = trim(dnstring) // ' ' // trim(dimnames(irev))
      endif
    enddo
  elseif (itstep .gt. 0) then
    dset_ndims = ndims + 1
    dset_dims(1) = itstep
    ext_dims(1) = 1
    chunk_dims(1) = 1
    dnstring = 't'
    do i = 1, ndims
      ! reverse the ordering, time is in the first slot
      ! so shift this up by one
      irev = (ndims - i) + 2
      dset_dims(irev) = dims(i)
      ext_dims(irev) = 0
      chunk_dims(irev) = dims(i)
      ! note that the dimnames array is not to be shifted
      dnstring = trim(dnstring) // ' ' // trim(dimnames(irev-1))
    enddo
  endif

return
END SUBROUTINE rhdf5_build_dims_for_write

!##############################################################################
Subroutine rhdf5_adjust_chunk_sizes (ndims,chunk_dims,dsize)

! This routine will attempt to keep the chunk size at optimal
! values for the dataset.
!
! First of all, the HDF5 interface uses a default chunk cache (hash
! table) with 521 entries where each entry is 1MB in size. For small
! datasets (data buffer for each write under 512KB) there is no
! need to change the chunk sizes since they default to the data
! dimension sizes making the entire data buffer one chunk.
!
! For large datasets (data buffer for each write > 512KB) diffferent
! sources have different recommendations for the optimal chunk size.
! The HDF5 documentation recommendation is: < 512KB. They 
! all agree that the thing that really matters is the access order.
!
! Since we are just streaming out the data time step by time step
! with no subset selection, we need to make the chunks fit in a
! contiguous line through the memory. We also would like to make
! all the chunks add up to the total data buffer size without
! any extra padding so that we keep the file size at a minimum.
! Should be able to accomplish this by dividing up each dimension
! starting with the first (slowest changing) dimension into smaller
! pieces (that don't leave any remainder) until the chunk size is
! at the right size. For exmaple if you have 4-byte data organized as
! (20, 50, 1000), 4 million bytes, then allow the chunk size
! (2, 50, 1000), 400 thousand bytes, but don't allow (3, 50, 1000),
! 600 thousand bytes (closer to 512KB) since 3 does not divide into
! 20 evenly resulting in unused space in the dataset.

implicit none

integer, parameter :: KB_512 = 524288
integer :: ndims, dsize
integer, dimension(ndims) :: chunk_dims
integer :: i, chunk_size, new_dim, shrink
integer :: rhdf5_find_even_divisor

  chunk_size = dsize
  do i = 1, ndims
    chunk_size = chunk_size * chunk_dims(i)
  enddo

  ! only adjust if the default chunk size is greater than 512KB
  if (chunk_size .gt. KB_512) then
    shrink = ceiling(float(chunk_size) / float(KB_512))

    ! shrink is the factor that we want to reduce the chunk size by,
    ! but we want the actual reduction to be an even divisor of the
    ! given chunk dimensions so that the chunks fit over the data
    ! exactly with no wasted space.
    !
    ! Walk through the chunk dimensions. Keep setting the current
    ! chunk dimension to 1 (and dividing the shrink factor by that
    ! dimension) while the shrink factor remains larger than that
    ! dimension. As soon as the dimension is greater than the shrink
    ! factor, divide up the dimension with the nearest even divisor
    ! less than the current dimension divided by the shrink factor.
    !
    ! The intent of this algorithm is to yield exact fitting tiles
    ! (chunks) over the data that have their size as close to 512KB
    ! as possible without going over 512KB.
    do i = 1, ndims
      if (shrink .gt. chunk_dims(i)) then
        shrink = ceiling(float(shrink) / float(chunk_dims(i)))
        chunk_dims(i) = 1
      elseif (shrink .gt. 1) then
        new_dim = floor(float(chunk_dims(i))/float(shrink))
        chunk_dims(i) = rhdf5_find_even_divisor(chunk_dims(i),new_dim)
        shrink = 1
      endif
    enddo
  endif

return
END SUBROUTINE rhdf5_adjust_chunk_sizes

!##############################################################################
Subroutine rhdf5_c2f_string (cstring,fstring,ssize)

! This routine will convert the C style string (null terminated) in
! cstring and convert it to a FORTRAN style string in fstring.

implicit none

integer :: ssize
character(len=ssize) :: cstring
character(len=ssize) :: fstring
integer :: i
integer :: null_pos

  do i = 1, ssize
    if (ichar(cstring(i:i)) .eq. 0) then
      null_pos = i
      exit
    endif
  enddo

  if (null_pos .gt. 1) then
    fstring = cstring(1:null_pos-1)
  else
    fstring = ''
  endif

return
END SUBROUTINE rhdf5_c2f_string

!##############################################################################
Subroutine rhdf5_read_dim_name_string (ndims,dimnames,dnstring)

! This routine will read the dimension name string (dnstring) and load
! up the entries in the dimnames array. The different names in dnstring
! are space separated.

implicit none

integer :: ndims,i,j
character(len=RHDF5_MAX_STRING), dimension(ndims) :: dimnames
character(len=RHDF5_MAX_STRING) :: dnstring

  ! initialize dimnames to null strings
  do j = 1, ndims
    dimnames(j) = ''
  enddo

  j = 1 ! index into dimnames
  do i = 1, len_trim(dnstring)
    if (dnstring(i:i) .ne. ' ') then
       ! keep appending dnstring(i) to dimnames(j)
       dimnames(j) = trim(dimnames(j)) // dnstring(i:i)
    else
       ! hit a space, bump j to next entry
       j = j + 1
    endif
  enddo

return
END SUBROUTINE rhdf5_read_dim_name_string

!##############################################################################
Subroutine rhdf5_set_var_to_dim (fileid,vname,dname)

! This routine will set a "dimension scale" on the given variable so that
! variable can be used as dimensions coordinate values in the hdf5 file.

implicit none

integer :: hdferr
integer*8 :: fileid, dsetid
character(len=*) :: vname,dname

  CALL rh5d_open (fileid, trim(vname)//char(0), dsetid, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim:'
    print*,'cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_set_var_to_dim: bad variable read'
  endif

  ! set the scale on the variable
  ! add a new attribute "axis" to the variable for GRADS
  CALL rh5ds_set_scale (dsetid, trim(dname)//char(0), hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim:'
    print*,'cannot set coordinate variable to a dimension'
    stop 'hdf5_set_var_to_dim: bad variable set dimension'
  endif

  CALL rh5a_write_anyscalar (dsetid,'axis'//char(0),trim(dname)//char(0) &
                            ,RHDF5_TYPE_STRING, hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_set_var_to_dim:'
    print*,'cannot add "axis" attribute to variable: ', trim(vname)
    stop 'hdf5_set_var_to_dim: bad variable set dimension'
  endif

  ! close the dataset
  CALL rh5d_close (dsetid,hdferr)

return
END SUBROUTINE rhdf5_set_var_to_dim

!##############################################################################
Subroutine rhdf5_attach_dims_to_var (fileid,vname)

! This routine will attach dimensions specs to the given variable. This
! is done via the HDF5 dimension scaling feature, and the purpose of
! doing the attachment is so that GRADS can read in the resulting
! HDF5 file directly without the use of a descriptor file.

implicit none

integer :: i,ndims,hdferr
integer*8 :: fileid, dsetid, dsclid
character(len=*) :: vname
integer, dimension(RHDF5_MAX_DIMS) :: dims
character(len=RHDF5_MAX_STRING) :: units,descrip,dnstring,stemp,coord_name
character(len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames

  ! open dataset
  CALL rh5d_open (fileid,trim(vname)//char(0),dsetid,hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_attach_dims_to_var:'
    print*,'cannot open dataset for variable: ',trim(vname)
    stop 'hdf5_attach_dims_to_var: bad variable read'
  endif

  ! dimensions
  CALL rh5d_read_get_dims (dsetid,ndims,dims,hdferr)
  if (hdferr .ne. 0) then
    print*,'ERROR: hdf5_attach_dims_to_var:'
    print*,'cannot obtain dimension information for variable: ',trim(vname)
    stop 'hdf5_attach_dims_to_var: bad variable read'
  endif

  ! read the DimNames attribute
  CALL rh5a_read_anyscalar (dsetid,'DimNames'//char(0),stemp &
                           ,RHDF5_TYPE_STRING,hdferr)
  CALL rhdf5_c2f_string (stemp,dnstring,RHDF5_MAX_STRING)
  CALL rhdf5_read_dim_name_string (ndims,dimnames,dnstring)

  ! There is no need to reverse the dimensions since we are going immediately
  ! back into the C interface to do the dimension scale attach.
  !
  ! Assume that the coordinate information is stored in datasets with names
  ! that reflect the dimension names as shown below:
  !     dim name            dataset name with coordinates
  !       x                     x_coords
  !       y                     y_coords
  !       z                     z_coords
  !       t                     t_coords
  !      <n>                    <n>_coords  (in general)

  do i = 1, ndims
    coord_name = trim(dimnames(i)) // '_coords'
    CALL rh5d_open (fileid,trim(coord_name)//char(0),dsclid,hdferr)
    if (hdferr .ne. 0) then
      print*,'ERROR: hdf5_attach_dims_to_var:'
      print*,'cannot open dataset for coordinates: ',trim(coord_name)
      stop 'hdf5_attach_dims_to_var: bad variable read'
    endif

    ! do the attach, C indices start with zero so send (i-1) to this routine
    CALL rh5ds_attach_scale (dsetid,dsclid,i-1,hdferr)
    if (hdferr .ne. 0) then
      print*,'ERROR: hdf5_attach_dims_to_var:'
      print*,'cannot attach coordinate data to variable:'
      print*,'ERROR:  Variable dataset: ',trim(vname)
      print*,'ERROR:  Coordinate dataset: ',trim(coord_name)
      stop 'hdf5_attach_dims_to_var: bad variable attach'
    endif

    CALL rh5d_close (dsclid,hdferr)
  enddo


  ! close the dataset
  CALL rh5d_close (dsetid,hdferr)

return
END SUBROUTINE rhdf5_attach_dims_to_var

END MODULE rhdf5_utils

!##############################################################################
integer Function rhdf5_find_even_divisor (num, div)

implicit none

integer :: num,div,i

rhdf5_find_even_divisor = 1
do i = div, 1, -1
  if (mod(num, i) .eq. 0) then
    rhdf5_find_even_divisor = i
    exit
  endif
enddo

return
END FUNCTION rhdf5_find_even_divisor
