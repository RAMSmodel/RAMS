!##############################################################################
Module grid_dims

implicit none

! This contains very basic specification of grid dimensions and other 
!    parameters that will be used to dimension arrays and allocate memory.

! Set maximum values of parameters:

integer, parameter ::  &
  maxargs      = 10    & ! Max # of command line arguments (using only 3 now)
 ,strl1        = 256   & ! Long character string length
 ,strl2        = 256   & ! Longer line read character string length
 ,maxgrds      = 8     & ! Max # of grids
 ,nxpmax       = 3030  & ! Max # of points in x-direction
 ,nypmax       = 3030  & ! Max # of points in y-direction
 ,nzpmax       = 332   & ! Max # of points in z-direction
 ,nzgmax       = 20    & ! Max # of soil levels
 ,nzsmax       = 20    & ! Max # of surface water / snow levels
 ,maxkppz      = 200   & ! Max # of KPP ocean model levels
 ,maxdim       = 3030  & ! Largest of NXPMAX,NYPMAX,NZPMAX+10,NZGMAX
 ,maxsclr      = 350   & ! Max # of scalars
 ,maxrevu      = 1000  & ! Max # of revu variables in REVU_IN list
 ,maxfiles     = 10000 & ! Max # of analysis files for REVU to read
 ,maxvars      = 1000  & ! Max # of variables (3d + 2d + leaf)
 ,maxrec       = 1000  & ! Max record length (lines) of namelists
 ,maxvalues    = 300   & ! Max # of tokens to be read in from namelists
 ,maxmach      = 2048  & ! Max # of parallel processors
 ,maxlite      = 99    & ! Max # of lite variables
 ,maxsstfiles  = 2000  & ! Max # of SST total files
 ,maxndvifiles = 2000  & ! Max # of NDVI total files
 ,maxsstdata   = 100   & ! Max # of SST file times
 ,maxndvidata  = 100   & ! Max # of NDVI file times
 ,maxmiss      = 1000  & ! Max # of missing surface data block files
 !
 !Configuration for isentropic data analysis package maximum sizes
 ,maxpr        = 100   & ! Max # of vertical levels allowed in the pressure data
 ,maxisn       = 200   & ! Max # of vertical levels allowed in the isentropic analysis
 ,maxx         = 1000  & ! Max # of X (west-east) grid points in RAMS or pressure grids
 ,maxy         = 1000  & ! Max # of Y (north-south) grid points in RAMS or pressure grids
 ,maxtimes     = 200   & ! Max # of data analysis times that can be processed in a single run.
 ,maxagrds     = 10    & ! Max # of RAMS grids that can have varfiles generated
 ,maxsigz      = 200   & ! Max # of vertical levels allowed in sigma-z analysis.
 ,maxlev       = 500   & ! Max # of levels in an input rawinsonde
 ,maxsname     = 10000 & ! Max # of input observations
 ,maxisfiles   = 1000  & ! Max # of input data times
 !
 !Configuration for observational data assimilation maximum sizes
 ,maxnudfiles  = 500   & ! Max # of nudging variable initialization files
 ,maxodafiles  = 1000  & ! Max # of ODA files
 ,maxodasta    = 5000  & ! Max # of station files
 ,maxodagrids  = 10    & ! Max # of ODA grids
 ,maxodanzp    = 200   & ! Max # of ODA NZP vertical levels
 ,maxupalevs   = 100   & ! Max # of upper air data levels
 ,maxkobs      = 10000 & ! Max # of vertical observations
 ,maxodatimes  = 3000  & ! Max # of ODA times (3*maxodafiles)
 !
 !Configuration for sounding input
 ,maxsndg      = 20000   ! Max # of sounding levels for input

END MODULE grid_dims

