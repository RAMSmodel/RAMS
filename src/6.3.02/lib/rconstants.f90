!##############################################################################
Module rconstants

implicit none

real, parameter ::                    &
        rgas     = 287.               & !Dry gas constant (J/kg/K)
    ,   cp       = 1004.              & !Specific heat at const pressure (J/kg/K)
    ,   cv       = 717.               & !Specific heat at const volumne (J/kg/K)
    ,   rm       = 461.               &
    ,   p00      = 1.e5               & !Standard pressure (Pa)
    ,   t00      = 273.16             & !Freezing point (K)
    ,   g        = 9.80               & !Gravity (m/s^2)
    ,   pi180    = 3.1415927 / 180.   &
    ,   pi4      = 3.1415927 * 4.     &
    ,   spcon    = 111120.            &
    ,   erad     = 6367000.           & !Earth radius (m)
    ,   vonk     = 0.40               & !von Karman constant
    ,   tkmin    = 5.e-4              &
    ,   alvl     = 2.50e6             & !Latent heat of vaporization at 0C (J/kg)
    ,   alvi     = 2.834e6            & !Latent heat of sublimation at 0C (J/kg)
    ,   alli     = 0.334e6            & !Latent heat of fusion at 0C (J/kg)
    ,   alvl2    = 6.25e12            &
    ,   alvi2    = 8.032e12           &
    ,   solar    = 1.3533e3           &
    ,   stefan   = 5.6696e-8          &
    ,   cww      = 4218.              &
    ,   c0       = 752.55 * 4.18684e4 &
    ,   viscos   = .15e-4             &
    ,   rowt     = 1.e3               &
    ,   dlat     = 111120.            &
    ,   omega    = 7.292e-5           &
    ,   rocp     = rgas / cp          &
    ,   p00i     = 1. / p00           &
    ,   cpor     = cp / rgas          &
    ,   rocv     = rgas / cv          &
    ,   cpi      = 1. / cp            &
    ,   cpi4     = 4. * cpi           &
    ,   cp253i   = cpi / 253.         & 
    ,   allii    = 1. / alli          &
    ,   aklv     = alvl / cp          &
    ,   akiv     = alvi / cp          &
    ,   gama     = cp / cv            &
    ,   gg       = .5 * g             &
    ,   ep       = rgas / rm          & 
    ,   p00k     = 26.870941          &  !  = p00 ** rocp  
    ,   p00ki    = 1. / p00k

END MODULE rconstants

