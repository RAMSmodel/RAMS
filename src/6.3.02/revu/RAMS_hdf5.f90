!##############################################################################
Subroutine rams_hdf5 (rh5_file,a,iztrans,ivtype,nngd,n1,n2,n3  &
                  ,zlevs,ns,idate,itime,ivtime,fcstsec  &
                  ,cvar,cdname,cdunits,sfclat,sfclon)


!-----------------------------------------------------------------
! Routine to write RAMS fields in HDF5 format as chosen through
! REVUIN. This is intended to be modified at will to produce the
! desired output. Basically, one 3 or 2-dimensional variable at a
! time will be sent here.

! Arguments:
! ----------
! cdname - name of variable
! cdunits - units of variable
! a - data
! n1,n2,n3 - actual dimensions of field
! iztrans - type of vertical transformation 1-sigma-z, 2-Cartesian, 3-pressure
! ivtype - type of variable 2-2d surface, 3-3d atmospheric
! nib,nie,njb,nje - horizontal or vertical "window" as chosen in namelist.
! nii,njj - number of windowed points
! n3 - number of atmospheric height levels (if sigma_z or Cartesian)
! zlevs(n3) - atmospheric coordinate levels (if sigma_z or Cartesian)
! ns - number of soil coordinate levels
! slz(ns) - soil coordinate levels
! iyear1 - year of model start
! imonth1 - month of model start
! idate1 - date of model start
! itime1 - time of model start (hhmm)
! fcstsec - seconds into run
!-----------------------------------------------------------------

! DESIGN NOTE:
!
! The resulting file needs to end up with:
! one dataset per variable that was requested in the revuvar() control
! one dataset per dimension that contains 1D arrays with the coordinate values
!
! The dimensions are named:
!     t    (time)
!     x    (longitude)
!     y    (latitude)
!     z    (levels)
!
! and the corresponding datasets (that contain the coordinate values)
! are named:
!     t_coords
!     x_coords
!     y_coords
!     z_coords
!
!     Note that the code in the dimension attaching routine will
!     take the name of the dimension ('t' for example) and append
!     '_coords' (forming 't_coords') to create the name of the
!     dataset that contains the coordinates.
!
! For GRADS:
!   the dimension datasets need to be turned into HDF5 dimension scales
!   and attached to their corresponding dimensions within the variable
!   datasets.
!
!   all datasets must be flat (ie, live under the root group)
!
!   all datasets must be made of atomic (ie, no compound) types
!
!   follow COARDS conventions for attribute names and values on the
!   coordinate datasets
! 
! Example: have zonal velocity field (3D) in the variable 'u' for multiple
! time steps --> 4D variable. This requires:
!
!     u be organized as: (t,z,y,x)
!
!     create dimension datasets:
!
!       t_coords - number of seconds into the simulation for each time step
!       z_coords - level values for grid
!       y_coords - latitude values for grid
!       x_coords - longitude values for grid
!
!     then t_coords gets attached to the first dimension of u
!          x_coords gets attached to the second dimension of u, etc.
!
! This allows the resulting hdf5 file to be directly readable by GRADS,
! as well as MATLAB and IDL.
!
! keep RHDF5_TSTRING_LEN in sync (matching length) with write format
! in rams_hdf5_create_time_string, which is currently in unused_routines

use an_header
use rhdf5_utils
use rcommons
use mem_grid

implicit none

integer, parameter :: RHDF5_TSTRING_LEN = 20
integer :: iztrans,ivtype,n1,n2,n3,ns,idate,itime,ivtime,nngd,i,j,k &
          ,deflvl
integer*8 :: rh5_file
real :: fcstsec,a(n1,n2,n3),zlevs(n3),sfclat(n1,n2),sfclon(n1,n2)
character(len=*) :: cvar,cdname,cdunits
character(len=strl1) :: flnm1
logical,dimension(maxgrds),save :: first_call_for_grid=(/(.true.,i=1,maxgrds)/)
integer,dimension(maxgrds),save :: prior_ivtime=(/(0,i=1,maxgrds)/)
integer, save :: prior_grid = 0
real, dimension(nplevs) :: plevs
real, dimension(:), allocatable :: vdata,lat_data,lon_data,xdata,ydata,zdata
integer, dimension(RHDF5_MAX_DIMS) :: sdims,lat_dims,lon_dims,zdims,cdims,tdims
real, dimension(max(nzg,nzs,npatch)) :: glevs
character(len=RHDF5_MAX_STRING) :: zunits, zdescrip, tunits
character(len=RHDF5_TSTRING_LEN) :: tstring
character(len=RHDF5_MAX_STRING), dimension(RHDF5_MAX_DIMS) :: dimnames

!DEBUG
!print*,'===> RAMS_hdf5 => fcstsec: ',fcstsec
!print*,'DEBUG: ivtime: ', ivtime
!print*,'DEBUG: cdname,n1,n2,n3: ',cdname(1:len_trim(cdname)-1),n1,n2,n3
!print*,'DEBUG: iztrans,ivtype: ',iztrans,ivtype
!print*,'DEBUG: nii,nib,nie,niinc: ',nii,nib,nie,niinc
!print*,'DEBUG: njj,njb,nje,njinc: ',njj,njb,nje,njinc
!print*,'DEBUG: nnb,nne,nninc: ',nnb,nne,nninc
!print*,'DEBUG: zlevs: ', zlevs
!print*,'DEBUG: iplevs: ', iplevs
!print*,'DEBUG: cvar, cdname, cdunits: ', cvar, cdname, cdunits
!print*,'DEBUG: nngd, first_call_for_grid: ', nngd, first_call_for_grid
!flush(6)

! Data set for this call is stored in: /<grid>/<variable>
! which gets organized in the file as:
!    Group: /<grid>
!    Dataset: <variable>
!
! It is assumed that the grid for this variable has been opened, 
! ie. it is up to the caller to open the grid

! Set the variable name, description (cdname) and units. Note that cdname 
! and cdunits come with a trailing semicolon which needs to be stripped off;
! Then select the data (XVAR, YVAR, TVAR specs) and write it out attached to
! the <grid> group.The 5th argument to rhdf5_write_variable is the time step.
! rhdf5_write_variable will tag on the time step as the size of an additional
! dimension to what is specified with ndims (4th arg set to ivtype). It is 
! expected that the time step will increase by one for each successive call.
dimnames(1) = 'x'
dimnames(2) = 'y'
dimnames(3) = 'z'

! write out the data. deflvl is the compression level (1-9)
! note that compression levels higher than 1 from REVU output
! can produce odd results in some fields such as reflectivity.
! Please use caution with this variable.
deflvl=6

CALL rams_hdf5_select_var_data (a,n1,n2,n3,nib,nie,niinc &
         ,njb,nje,njinc,nnb,nne,nninc,ivtype,sdims,vdata)
CALL rhdf5_write_variable_rdata (rh5_file,cvar,ivtype,ivtime,sdims &
         ,cdunits(1:len_trim(cdunits)-1),cdname(1:len_trim(cdname)-1) &
         ,dimnames,deflvl,vdata)
deallocate(vdata)

! If it is the first time in this routine for this grid, then store lat, 
! lon and zlevs, in the grid group. Last arg on rams_hdf5_select_var_data is 
! number of dimensions.
if (first_call_for_grid(nngd)) then
  ! Select the lat, lon values that correspond to the data selected 
  ! from the variables
  CALL rams_hdf5_select_var_data (sfclat,n1,n2,n3,nib,nie,niinc &
                ,njb,nje,njinc,nnb,nne,nninc,2,lat_dims,lat_data)
  CALL rams_hdf5_select_var_data (sfclon,n1,n2,n3,nib,nie,niinc &
                ,njb,nje,njinc,nnb,nne,nninc,2,lon_dims,lon_data)
  ! lat and lon have the same dimensions at this point
  CALL rams_hdf5_build_xycoords (lon_dims(1),lon_dims(2),lon_data &
                ,lat_data,xdata,ydata)

  ! Select the z values
  !
  ! rams_hdf5_select_var_data assumes that 1D data is being stored in the
  ! "i" position, that is: the first n1 entries of the input array hold the 
  ! data and ni* describe how to step through the data.
  !
  ! For sigma and cartesian vertical coordinates, zlevs(n3) holds the data for
  ! all levels; for pressure vertical coordinates, iplevs(nplevs) holds the
  ! data for all levels. For all types of vertical coordinates, nn* describe 
  ! how to step through the data. Therefore place the n3 (or nplevs) and the
  ! nn* variables in the positions where n1 and ni* would normally go. The 
  ! n2,n3,nj*,nn* arguments are ignored for the 1D case so it doesn't really
  ! matter what they are set to.
  !
  ! iztrans describes what the vertical levels are
  !    1 - sigma-z
  !    2 - Cartesian
  !    3 - pressure
  if (iztrans .eq. 1) then
    zunits = 'meter'
    zdescrip = 'sigma-z'
    CALL rams_hdf5_select_var_data (zlevs,n3,n2,n3,nnb,nne,nninc &
                     ,njb,nje,njinc,nnb,nne,nninc,1,zdims,zdata)
  elseif (iztrans .eq. 2) then
    zunits = 'meter'
    zdescrip = 'standard height values'
    CALL rams_hdf5_select_var_data (zlevs,n3,n2,n3,nnb,nne,nninc &
                     ,njb,nje,njinc,nnb,nne,nninc,1,zdims,zdata)
  elseif (iztrans .eq. 3) then
    zunits = 'millibar'
    zdescrip = 'standard pressure values'
    ! record iplevs as real values
    plevs = float(iplevs)
    CALL rams_hdf5_select_var_data (plevs,nplevs,n2,n3,nnb,nne,nninc &
                     ,njb,nje,njinc,nnb,nne,nninc,1,zdims,zdata)
  elseif (iztrans .eq. 4) then
    zunits = 'none'
    zdescrip = 'integer levels or patches'
    glevs=(/(i,i=1,max(nzg,nzs,npatch))/)
    CALL rams_hdf5_select_var_data (glevs,max(nzg,nzs,npatch),n2,n3 &
     ,nnb,max(nzg,nzs,npatch),nninc,njb,nje,njinc,nnb,nne,nninc,1,zdims,zdata)
  else
    stop 'iztrans should only be 1,2,3,4'
  endif

  ! Write out the lat, lon and z data, 5th argument is the time step. Setting
  ! to 0 tells rhdf5_write_variable that these variables do not have a time
  ! dimension, just (x,y,z).
  !
  ! lat_data and lon_data are 2D arrays that contain the latitude and
  ! longitude values at every (selected) grid point. These are useful for
  ! applications that need to preserve locations from the Earth's surface.
  !
  ! Most applications that will be used to produce plots or generate 
  ! diagnostics want 1D versions of the latitude and longitude (ie, values 
  ! for the x and y coordinates of the grid). Use the middle column of 
  ! lon_data for the x coordinate values, and the middle row of lat_data for
  ! the y coordinate values. This should produce good enough values for a lot 
  ! of circumstances. When this isn't good enough the 2D versions of lat and
  ! lon can be used.
  !
  ! Place the 1D versions of lon and lat into x_coords and y_coords
  ! respectively; and place the 2D versions of lon and lat in sfclon and
  ! sfclat respectively.
  dimnames(1) = 'x'
  CALL rhdf5_write_variable_rdata (rh5_file, 'x_coords', 1, 0, lon_dims(1) &
              , 'degrees_east', 'longitude', dimnames, deflvl, xdata)

  dimnames(1) = 'y'
  CALL rhdf5_write_variable_rdata (rh5_file, 'y_coords', 1, 0, lat_dims(2) &
              , 'degrees_north', 'latitude', dimnames, deflvl, ydata)

  dimnames(1) = 'z'
  CALL rhdf5_write_variable_rdata (rh5_file, 'z_coords', 1, 0, zdims, zunits &
              , zdescrip, dimnames, deflvl, zdata)

  deallocate(zdata)
  deallocate(lat_data)
  deallocate(lon_data)
  deallocate(xdata)
  deallocate(ydata)

  first_call_for_grid(nngd) = .false.
endif

! Write out the time values. Each grid is in a separate file so keep track
! of which time step you are on per grid. The calls to this routine will go
! in order from time step (ivtime) 1...Nt for each grid and each grid has 
! it's own file. This works out nice because we never shrink the time
! dimension (which is extendable). Note that the hdf5 interface does not 
! support shrinking extendable dimensions at this time so if this sequence 
! changes, there may be more work to be done to support shrinking the time 
! dimension.
!
! Present the time as real numbers respresenting the number of seconds into
! the simulation (fcstsec). Make the dataset 1D (t dimension only) and use
! the name 't_coords' which allows the code that attaches this dimension data
! to the variable datasets to a

if (ivtime .ne. prior_ivtime(nngd)) then
  ! ndims (third arg) set to zero indicates to rhdf5_write_variable that
  ! these are time values

  CALL rams_hdf5_create_time_units (iyear1, imonth1, idate1, itime1, tunits)
  dimnames(1) = ''
  tdims(1) = 0
  CALL rhdf5_write_variable_fsngl (rh5_file, 't_coords', 0, ivtime, tdims, tunits &
       , 'simulation time', dimnames, deflvl, fcstsec)

  prior_ivtime(nngd) = ivtime
endif

return

!**************************************************************************
! Routines called by rams_hdf5 go here
! use the 'contains' keyword so that the compiler can do more checking
! (argument types, etc.) 'contains' also allows an allocatable array to be
! passed to a routine where that routine allocates the array
!**************************************************************************
contains

!***********************************************************************
! This routine will select data out of the whole field (input "a") according
! to the revu specs XVAR, YVAR, and ZVAR. The information of how to step
! through the "a" array is contained in the nib,nie,niinc,njb,nje,njinc,nnb
! ,nne,nninc inputs.
!
! "a" is considered to be organized as a(i,j,k) regardless of 1D, 2D or 3D
! data meaning that the caller must place the input data into:
!    the first n1 entries for 1D data
!    the first n1,n2 entries for 2D data
!    the n1,n2,n3 (all) entries for 3D data
! This will work out since FORTRAN accesses arrays in a column major fashion
! (that is, the first subscript varies the fastest in the memory image)
!
! The stepping information for selecting the data out of the "a" array are
! given:
!    "i" direction: nib, nie, niinc (beginning, end, increment)
!    "j" direction: njb, nje, njinc
!    "k" direction: nnb, nne, nninc
!
! ndim discerns between 1D, 2D and 3D data which will determine how the
! dataset in rh5_file is allocated

Subroutine rams_hdf5_select_var_data (a,n1,n2,n3,nib,nie,niinc,njb,nje,njinc &
                                    ,nnb,nne,nninc,ndims,dims,sdata)

implicit none

  integer :: n1,n2,n3,nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc,ndims
  real, dimension(n1,n2,n3) :: a
  integer, dimension(*) :: dims
  real, dimension(:), allocatable :: sdata

  ! Want the data store in the rh5_file dataset to be contiguous so we can
  ! simply write it out into the hdf5 file in one call. Figure out the new
  ! dimensions of the data store after the selection has taken place and 
  ! allocate accordingly.
  !
  ! The formula for calculating the number of points in each dimension, using
  ! integer arithmetic, is:
  !
  !      num_points = ((end - beginning)/increment) + 1
  !
  dims(1) = ((nie - nib)/niinc) + 1
  dims(2) = ((nje - njb)/njinc) + 1
  dims(3) = ((nne - nnb)/nninc) + 1

  ! Allocate the space for the rh5_file dataset, and then copy the data from
  ! the "a" array. The dataset in rh5_file is a 1D array so it can be made 
  ! to look like any of a 1D, 2D or 3D structure. This depends on calling a 
  ! routine where you can cast the 1D, 2D or 3D structure on top of the 
  ! 1D allocated storage.
  if (allocated(sdata)) then
    deallocate(sdata)
  endif

  if (ndims .eq. 1) then
    ! 1D data
    allocate(sdata(dims(1)))
    CALL rams_hdf5_select_1D_data (sdata,a,n1,dims(1),nib,nie,niinc)
  elseif (ndims .eq. 2) then
    ! 2D data (idim_type=2)
    allocate(sdata(dims(1)*dims(2)))
    CALL rams_hdf5_select_2D_data (sdata,a,n1,n2,dims(1),dims(2),nib,nie &
                                 ,niinc,njb,nje,njinc)
  elseif (ndims .eq. 3) then
    ! 3D data (idim_type=3)
    allocate(sdata(dims(1)*dims(2)*dims(3)))
    CALL rams_hdf5_select_3D_data (sdata,a,n1,n2,n3,dims(1),dims(2),dims(3), &
      nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc)
  elseif (ndims .eq. 4) then
    ! Snow layer data (patch sum) (idim_type=4)
    ndims=3
    allocate(sdata(dims(1)*dims(2)*dims(3)))
    CALL rams_hdf5_select_3D_data (sdata,a,n1,n2,n3,dims(1),dims(2),dims(3), &
      nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc)
  elseif (ndims .eq. 5) then
    ! Soil layer data (patch sum) (idim_type=5)
    ndims=3
    allocate(sdata(dims(1)*dims(2)*dims(3)))
    CALL rams_hdf5_select_3D_data (sdata,a,n1,n2,n3,dims(1),dims(2),dims(3), &
      nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc)
  elseif (ndims .eq. 6) then
    ! Patch specific variables (idim_type=6)
    ndims=3
    allocate(sdata(dims(1)*dims(2)*dims(3)))
    CALL rams_hdf5_select_3D_data (sdata,a,n1,n2,n3,dims(1),dims(2),dims(3), &
      nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc)
  else
    print*,'ERROR: rams_hdf5_select_var_data: ndims = ',ndims,dims(1),dims(2),dims(3)
    print*,'Only know how to handle 1D, 2D and 3D arrays'
    stop 'rams_hdf5_select_var_data: bad dimensions'
  endif

return
END SUBROUTINE rams_hdf5_select_var_data

!##############################################################################
Subroutine rams_hdf5_select_1D_data (dat,a,n1,num_i,nib,nie,niinc)

! This routine will select out 1D data from the "a" array into the
! rh5_file dataset

implicit none

integer :: n1,num_i,nib,nie,niinc,ii
real, dimension(n1) :: a
real, dimension(num_i) :: dat

ii = 1
do i = nib, nie, niinc
  dat(ii) = a(i)
  ii = ii + 1
enddo

return
END SUBROUTINE rams_hdf5_select_1D_data

!##############################################################################
Subroutine rams_hdf5_select_2D_data (dat,a,n1,n2,num_i,num_j,nib,nie &
                                   ,niinc,njb,nje,njinc)

! This routine will select out 2D data from the "a" array into the
! rh5_file dataset

implicit none

integer :: n1,n2,num_i,num_j,nib,nie,niinc,njb,nje,njinc,ii,jj
real, dimension(n1,n2) :: a
real, dimension(num_i,num_j) :: dat

jj = 1
do j = njb, nje, njinc
  ii = 1
  do i = nib, nie, niinc
    dat(ii,jj) = a(i,j)
    ii = ii + 1
  enddo
  jj = jj + 1
enddo

return
END SUBROUTINE rams_hdf5_select_2D_data

!##############################################################################
Subroutine rams_hdf5_select_3D_data (dat,a,n1,n2,n3,num_i,num_j,num_k &
                         ,nib,nie,niinc,njb,nje,njinc,nnb,nne,nninc)

! This routine will select out 3D data from the "a" array into the
! rh5_file dataset

implicit none

integer :: n1,n2,n3,num_i,num_j,num_k,nib,nie,niinc &
          ,njb,nje,njinc,nnb,nne,nninc,ii,jj,kk
real, dimension(n1,n2,n3) :: a
real, dimension(num_i,num_j,num_k) :: dat

kk = 1
do k = nnb, nne, nninc
  jj = 1
  do j = njb, nje, njinc
    ii = 1
    do i = nib, nie, niinc
      dat(ii,jj,kk) = a(i,j,k)
      ii = ii + 1
    enddo
    jj = jj + 1
  enddo
  kk = kk + 1
enddo

return
END SUBROUTINE rams_hdf5_select_3D_data

!##############################################################################
Subroutine rams_hdf5_create_time_units (iyear1,imonth1,idate1,itime1,tunits)

implicit none

integer :: iyear1,imonth1,idate1,itime1,ihour1,imin1
character(len=*) :: tunits

! itime1 is in the format: hhmm
ihour1 = itime1 / 100
imin1 = itime1 - (ihour1 * 100)

write(tunits, '(a14,i0.4,a1,i0.2,a1,i0.2,a1,i0.2,a1,i0.2,a9)') &
 'seconds since ', iyear1, '-', imonth1, '-', idate1, ' ' &
 , ihour1, ':', imin1, ':00 00:00'

return
END SUBROUTINE rams_hdf5_create_time_units

!##############################################################################
Subroutine rams_hdf5_build_xycoords (n1,n2,lon_data,lat_data,xdata,ydata)

! This routine will select a representative column from the 2D longitude
! data (lon_data) to place in xdata, and a representative row from
! the 2D latitude data (lat_data) to place in ydata.
!
! Use the middle column (row) for xdata (ydata) for now.
!   formula for "middle"
!     if the "other" dimension is < 3, then set the middle to 1
!     else
!       integer divide the "other" dimension by 2
!       add 1 to result
!       (this gets the middle if the other dim is odd, next entry
!       past middle if the other dim is even)

implicit none

integer :: n1,n2,i,mid
real, dimension(n1,n2) :: lon_data, lat_data
real, dimension(:), allocatable :: xdata,ydata

mid = 1
if (n2 .gt. 2) then
  mid = (n2 / 2) + 1
endif
allocate(xdata(n1))
! xdata gets first column of lon_data
do i = 1, n1
  xdata(i) = lon_data(i,mid)
enddo

mid = 1
if (n1 .gt. 2) then
  mid = (n1 / 2) + 1
endif
allocate(ydata(n2))
! ydata gets first row of lat_data
do i = 1, n2
  ydata(i) = lat_data(mid,i)
enddo

return
END SUBROUTINE rams_hdf5_build_xycoords

END SUBROUTINE rams_hdf5

