!##############################################################################
Module kpp_parameters

use grid_dims, only:maxgrds !from RAMS

implicit none

integer :: IKPP !flag to activate KPP from RAMSIN

real :: & !KPP flags from RAMSIN
         FRQKPP & !KPP timestep (seconds)
     ,RELAX_SST & !timescale (days) for ocean column salinity relax/nudging
    ,RELAX_OCNT & !timescale (days) for ocean column temp relax/nudging
     ,RELAX_SAL & !timescale (days) for SST relax/nudging
       ,DMAXKPP & !maximum ocean depth (meters)
     ,DSCALEKPP & !exponential for stretched grid
     ,UBMN_KPP    !ubmin minimum wind speed used for wind stress

integer :: & !KPP flags from RAMSIN
           NKPPZ & !number of kpp grid point levels
     ,KPPITERMAX & !max number of iterations allowed
         ,KPPRNT   !KPP data print flag

integer :: & !parameters dynamically set at initialization
           kppNZ & !number of kpp layers used in KPP (=nkppz-1)
        ,kppNZP1   !number of kpp grid points used in KPP (=nkppz)

integer, parameter :: &
         kppNVEL = 2  & !number of velocity components, i.e. 2 for U and V
       ,kppNSCLR = 2    !number of scalars, i.e. T, S

integer, parameter :: &
             iprnt=50 &
            ,jprnt=85


real, dimension(maxgrds) :: hmixdepth,hmixdepth_max

END MODULE kpp_parameters
