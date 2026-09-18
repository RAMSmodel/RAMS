

! CVS:  $Id: bugs_lwr.F,v 1.7 2006/11/16 18:45:09 norm Exp $
! CVS:  $Name:  $
! Modified for new version of planck function

!-----------------------------------------------------------------------

      subroutine bugs_lwr &
                     (    ncol ,    nlm ,    pp ,    ppl, &
                            dp ,     tt ,  rmix ,  cwrho, &
                           cwn ,  &
                         cirho ,  o3mix ,    ts , cldamt, &
                        cldmax ,     b1 ,    b2 ,     b3, &
                            b4 ,  umco2 , umch4 ,  umn2o, &
                          fdlw , fdlwcl ,  fulw , fulwcl, &
                     sel_rules &
                     )
     
      use bugs_kinds
      use bugsrad_planck, only:  planck
      use gases_ckd, only:  gases, stanpir,pscale
      use continuum
      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! bugs_lwr replaces crclwr written by G. Stephens. bugs_lwr computes the
! downward and upward longwave radiative fluxes, and longwave heating
! rates.
! Laura D. Fowler (slikrock/08-23-96).

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu.

! MODIFICATIONS:
! * moved the computation of the all-sky and clear-sky radiative heating
!   rates to bugs_rad.
!   Laura D. Fowler and Phil Partain/slikrock (01-27-99).

! * added effective radii of cloud droplets and ice crystals that are
!   dependent on the cloud water and cloud ice contents.
!   Laura D. Fowler/slikrock (06-08-00).

! * cleaned up the argument list to remove variables related to short
!   wave radiative transfer.
!   Laura D. Fowler/slikrock (02-01-00).

! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:

!     pscale         : pressure scaling.
!     cloudg          : computes cloud optical properties of water/ice
!                       clouds.
!     gascon_ckd_parm : water vapor continuum absorption.
!     plank           : computes planck function.
!     comscp1         : combines optical properties for gray absorption
!                       (clouds and water vapor continuum).
!     comscp2         : combines optical properties for non-gray gaseous
!                       absorption.
!     gases           : computes gases absorption.
!     two_rt_lw       : two-stream parameterization.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.

!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
       sel_rules

      integer (kind=int_kind), intent(in):: &
       ncol,  &!Length of sub-domain.
       nlm   !Number of layers.

      real (kind=dbl_kind), intent(in):: &
       umco2, & !Concentration of CO2                              (ppm).
       umch4, & !Concentration of CH4                              (???).
       umn2o !Concentration of N2o                              (???).

      real (kind=dbl_kind), intent(in),  dimension(ncol):: &
       ts,  &  !Surface temperature                                 (K).
       cldmax !Maximum cloud fraction                              (-). 

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
       ppl, & !Layer pressure                                 (hPa).
       dp,    & !Layer thickness                               (hPa).
       tt,    & !Temperature                                   (K).
       rmix,  & !Water vapor mixing ratio                      (kg/kg).
       cwrho, & !Cloud water mixing ratio                      (g/m^3).
       cwn,   & !Cloud droplet concentration                   (#/kg) .
       cirho, & !Cloud ice mixing ratio                        (g/m^3).
       o3mix, & !Ozone mixing ratio                            (kg/kg).
       cldamt, & !Cloud fraction                               (-).
       b1,  &  !Cloud overlap parameter                        (-).
       b2,  &  !Cloud overlap parameter                        (-).
       b3,  &  !Cloud overlap parameter                        (-).
       b4    !Cloud overlap parameter                         (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
       pp    !Level pressure                                    (hPa).


!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
       fdlw,   & !Downward LW flux                           (W/m^2).
       fdlwcl, & !Downward clear-ksy LW flux                 (W/m^2).
       fulw,   & !Upward LW flux                             (W/m^2).
       fulwcl  !Upward clear-sky LW flux                     (W/m^2).

! LOCAL VARIABLES:

      integer (kind=int_kind), parameter :: &
       mb = 18, &  !Total number of spectral intervals.
       mbs = 6,  &  !Number of shortwave (SW) spectral intervals.
       mbir = 12   !Number of longwave (LW) spectral intervals.
     

      integer (kind=int_kind)::    &
       i,    &  !Horizontal index.
       l,    &  !Vertical index.
       ib,   &  !Index of spectral interval.
       ig,   &  !Index of k-distribution.
       ibmbs   !Index of LW spectral interval.

      integer (kind=int_kind), dimension(ncol,nlm):: &
       ip1,  & !Used in conjunction with pressure weigthing.
       ip2   !Used in conjunction with pressure weigthing.

      real (kind=dbl_kind) :: &
       hk,    &   !Weighted spectral solar constant             (W/m^2).
       eps,   &   !Threshold for cloud optical properties              .       
       tmax,  &   !Temperature threshold                            (K).
       pdist
      data eps,tmax,pdist /1.e-05,340.,2./
          
      real (kind=dbl_kind), dimension(mbir):: &
       kg       !Nb of k-distributions per spectral intervals.  
      data kg /2,3,4,4,3,5,2,10,12,7,7,8/

      real (kind=dbl_kind), dimension(mbir):: &
        asym_wat,  &!Spectral asymmetry factor of water clouds.
        asym_ice !Spectral asymmetry factor of ice clouds.

      real (kind=dbl_kind), dimension(mb):: &
       cnrw,  &  !Real part of refractive index (Water clouds).
       cniw,  &  !Imaginary part of refractive index (Water clouds).
       cnri,  &  !Real part of refractive index (Ice clouds).
       cnii,  &  !Imaginary part of refractive indec (Ice clouds).
       xlam    !Center of spectral band.

      real (kind=dbl_kind), dimension(ncol,mbir):: &
      es       !Spectral surface emissivity                       (-).

      real (kind=dbl_kind), dimension(ncol,nlm):: &
       rew,     & !Effective radius for cloud water                 (mu).
       rei,     & !Effective radius for cloud ice                   (mu).
       ttem,    & !Local temperature                                 (K).
       pkd,     & !
       tau1,    & !All-sky optical depth                             (-).
       tauclr1, & !Clear-sky optical depth                           (-).
       tau,     & !All-sky optical depth                             (-).
       tauclr,  & !Clear-sky optical depth                           (-).
       taer,    & !Aerosol optical depth                             (-).
       tray,    & !Rayley optical depth                              (-).
       tg,      & !Gases optical depth                               (-).
       tgm,     & !WV continuum optical depth                        (-).
       tcldi,   & !Ice cloud optical depth                           (-).
       tcldw,   & !Water cloud optical depth                         (-).
       wc,      & !All-sky single scattering albedo                  (-).
       wcclr,   & !Clear-sky single scattering albedo                (-).
       waer,    & !Aerosol single scattering albedo                  (-).
       wray,    & !Rayley single scattering albedo                   (-).
       wcldi,   & !Ice cloud single scattering albedo                (-).
       wcldw,   & !Water cloud single scattering albedo              (-).
       asym,    & !All-sky asymmetry factor                          (-).
       asyclr,  & !Clear-sky asymmetry factor                        (-).
       asyaer,  & !Aerosol asymmetry factor                          (-).
       asycldi, & !Ice cloud asymmetry factor                        (-).
       asycldw, & !Water cloud asymmetry factor                      (-).
       fwclr,   & !
       fwcld   !

      real (kind=dbl_kind), dimension(ncol,nlm+1):: &
       bf,    &   !Planck function for layers                    (W/m^2).
       fdg,   &   !Spectral downward flux                        (W/m^2).
       fdgcl, &   !Spectral clear-sky downward flux              (W/m^2).
       fug,   &   !Spectral upward flux                          (W/m^2).
       fugcl   !Spectral clear-sky upward flux                (W/m^2).

!     longwave asymmetry parameters:
!     (assumes: re=10 for water; re=30 for ice)
      data asym_wat /0.8200, 0.8547, 0.8619, 0.8683, 0.8723, 0.8703, &
                     0.8566, 0.8040, 0.7463, 0.6579, 0.5103, 0.1279 /
      data asym_ice /0.8524, 0.8791, 0.9022, 0.8797, 0.8637, 0.8722, &
                   0.8609, 0.8168, 0.7663, 0.6584, 0.6172, 0.3585 /

!---  cnrw and cniw (water clouds):
      data cnrw/1.3422,1.3281,1.3174,1.2901,1.3348,1.3700,1.3191,1.2821, &
                1.3160,1.3030,1.2739,1.2319,1.1526,1.1981,1.3542,1.4917, &
                1.5463,1.8718/
      data cniw/6.4790e-9,1.3417e-06,1.2521e-4,7.1533e-4,4.2669e-2, &
                4.3785e-3,1.3239e-2 ,1.5536e-2,5.3894e-2,3.4346e-2, &
                3.7490e-2,4.7442e-2 ,1.2059e-1,3.3546e-1,4.1698e-1, &
                4.0674e-1,3.6362e-1 ,5.2930e-1/

!--- cnri and cnii (ice clouds):
      data cnri/1.3266,1.2986,1.2826,1.2556,1.2963,1.3956, & 
                1.3324,1.2960,1.3121,1.3126,1.2903,1.2295, &
                1.1803,1.5224,1.5572,1.5198,1.4993,1.7026/
      data cnii/7.0696e-9,9.1220e-7,1.2189e-4,5.7648e-4,4.3144e-2, &
                8.2935e-3,1.5540e-2,2.5594e-2,5.9424e-2,5.1511e-2, &
                4.0325e-2,4.7994e-2,2.3834e-1,3.0697e-1,1.1852e-1, &
                4.3048e-2,6.3218e-2,1.5843e-1/

!---- spectral band center:
      data xlam/0.45  ,1.0   ,1.6  ,2.2  ,3.0   ,3.75  ,4.878 ,5.556, &
                6.452 ,7.547 ,8.511,9.615,11.236,13.605,16.529,21.277, & 
                29.412,71.403/

!-----------------------------------------------------------------------

!---- 0. initialize output arrays:

      fdlw(:,:)   = 0.
      fdlwcl(:,:) = 0.
      fulw(:,:)   = 0.
      fulwcl(:,:) = 0.

      !rew(:,:)    = 3.
      !assuming gamma shape parameter of 4
      !rew(:,:) = 10.
      rew(:,:) = 6.e6*(cwrho/cwn*0.0019894e-6)**(1./3.)
      rei(:,:)    = 30.

      do l = 1, nlm
         do i = 1, ncol
            ttem(i,l) = min(tmax,tt(i,l))
          enddo
      enddo

!---- note: this will be changed to accomodate the spectral dependence
!     the surface emissivity:

      do ib = 1, mbir
         do i = 1, ncol
            es(i,ib) = 1.
         enddo
      enddo

!--   pressure scaling:

       call pscale(ncol,nlm,ppl,stanpir,pkd,ip1,ip2)

!---- 1. loop over the mbir spectral intervals starts here:

      do ib = mbs+1, mb                 
         ibmbs = ib - mbs
         
         tray(:,:)   = 0.
         wray(:,:)   = 0.
         taer(:,:)   = 0.
         waer(:,:)   = 0.
         asyaer(:,:) = 1.

!---- 1.1 optical properties of water and ice clouds (as in crclwr for
!        now):

         call cloudg &
                 (   ncol ,    nlm  ,    mb ,    ib, &
                       pp ,     tt  , cwrho ,   rew, &
                    pdist ,   cnrw  ,  cniw ,  cnri, &
                     cnii ,   xlam  , tcldw , wcldw, &
                  asycldw , .false. &
                 )

         call cloudg &
                 (   ncol ,   nlm   ,    mb ,    ib, &
                       pp ,    tt   , cirho ,   rei, &
                    pdist ,  cnrw   ,  cniw ,  cnri, &
                     cnii ,  xlam   , tcldi , wcldi, &
                  asycldi , .true. &
                 )

!     the asymmetry factor for water and ice clouds are fixed as for now
!     functions of the spectral intervals:

         do l = 1, nlm
            do i = 1, ncol
               if(cwrho(i,l) .ge. eps) asycldw(i,l) = asym_wat(ibmbs)
               if(cirho(i,l) .ge. eps) asycldi(i,l) = asym_ice(ibmbs)
            enddo
         enddo    

!---- 1.2 water vapor continuum:             

         call gascon &
                 (ncol , nlm, ib ,   pp, &
                   ppl ,  dp, tt , rmix, &
                   tgm &
                 )

!---- 1.3 planck function:
         call planck(ncol,nlm,ibmbs,ts,tt,bf)

!---- 1.4 combines single-scattering properties for gray absorption:

         call comscp1 &
                 (   ncol ,     nlm ,  taer ,   tcldi, &
                    tcldw ,     tgm ,  tray ,    waer, &
                    wcldi ,   wcldw ,  wray ,  asyaer, &
                  asycldi , asycldw ,  tau1 , tauclr1, &
                     asym ,  asyclr , fwcld ,   fwclr &
                 )
           
!---- loop over the k-probability distributions starts here:

         do ig = 1, kg(ibmbs)

!---- 1.5 gaseous absorption:

            call gases &
                    ( ncol ,   nlm ,    ib ,    ig, &
                        pp ,    dp ,    tt ,  rmix, &
                     o3mix , umco2 , umch4 , umn2o, &
                        hk ,    tg ,   pkd ,   ip1, &
                       ip2 &
                    )

!---- 1.6 combines all single-scattering properties:

            call comscp2 &
                    (  ncol ,  nlm ,      tg , fwcld, &
                      fwclr , tau1 , tauclr1 ,   tau, &
                     tauclr ,   wc ,   wcclr &
                    )


      !NBW - Minor (?) bug fix
      !With near-zero CO2 and very low water vapor amounts, the
      !correlated-K parameterization can generate negative optical
      !depths in the CO2-H2O overlap bands.  Here's a quick fix:
            where (tau .lt. 0)
               tau = 0.
            endwhere
            where (tauclr .lt. 0)
               tauclr = 0.
            endwhere

!---- 1.7 two-stream approximation:
! No overlap
            call two_rt_lw &
                    (     ncol , nlm ,  mbs , mbir, &
                            ib ,  wc , asym ,  tau, &
                            es ,  bf ,  fug ,  fdg, &
                     sel_rules &
                    )

            call two_rt_lw &
                    (     ncol ,   nlm ,   mbs  ,   mbir, &
                            ib , wcclr , asyclr , tauclr, &
                            es ,    bf ,  fugcl ,  fdgcl, &
                     sel_rules &
                    )

            fdlw(:,:)   = fdlw(:,:)   + fdg(:,:)*hk
            fulw(:,:)   = fulw(:,:)   + fug(:,:)*hk
            fdlwcl(:,:) = fdlwcl(:,:) + fdgcl(:,:)*hk
            fulwcl(:,:) = fulwcl(:,:) + fugcl(:,:)*hk
      !print'(5e20.10)',xlam(ib),fdg(1,2),fdg(1,nlm),fug(1,1),fug(1,nlm)
         enddo ! end k-distribution

      enddo ! end spectral interval
      

      return
      end subroutine bugs_lwr
      
!-----------------------------------------------------------------------


! CVS:  $Id: bugs_swr.F,v 1.4 2005/11/22 21:55:48 norm Exp $
! CVS:  $Name:  $

!-----------------------------------------------------------------------

      subroutine bugs_swr &
                      (  ncol ,     nlm ,       pp ,      ppl, &
                           dp ,      tt ,     rmix ,    cwrho, &
                          cwn , &
                        cirho ,   o3mix ,       ts ,     amu0,  &
                          slr ,   alvdf ,    alndf ,    alvdr,  &
                        alndr ,  cldamt ,   cldmax ,    umco2 , &
                        umch4 ,   umn2o ,       b1 ,       b2 , &
                           b3 ,      b4 ,     fdsw ,   fdswcl , &
                         fusw ,  fuswcl ,   radvbc , radvbccl , &
                       radvdc ,radvdccl ,   radnbc , radnbccl , &
                       radndc ,radndccl ,sel_rules &
                      )


      use bugs_kinds
      use gases_ckd, only: gases, stanps, pscale
      use rayleigh, only: rayle
      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! bugs_swr replaces crcswr written by G. Stephens. BUGSswr computes the
! downward and upward SW radiative fluxes, and SW heating rates.
! Laura D. Fowler (slikrock/08-20-97).

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu.

! MODIFICATIONS:
! * moved the computation of the all-sky and clear-sky radiative heating
!   rates to bugs_rad.
!   Laura D. Fowler and Phil Partain(slikrock/01-27-99).

! * added effective radii of cloud droplets and ice crystals that are
!   dependent on the cloud water and cloud ice contents.
!   Laura D. Fowler/slikrock (06-08-00).

! * cleaned up the argument list to remove variables related to short
!   wave radiative transfer.
!   Laura D. Fowler/slikrock (02-01-00).

! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:

!     pscale         : pressure scaling.
!     cloudg          : computes cloud optical properties of water/ice
!                       clouds.
!     rayle           : Computes Rayleigh scattering properties.
!     comscp1         : combines optical properties for gray absorption
!                       (clouds and water vapor continuum).
!     comscp2         : combines optical properties for non-gray gaseous
!                       absorption.
!     gases           : computes gases absorption.
!     two_rt_sw       : two-stream parameterization.

! FUNCTIONS CALLED:
!     none. 

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.

!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        sel_rules

      integer (kind=int_kind), intent(in):: &
        ncol,  & !Length of sub-domain.
        nlm    !Number of layers.

      real (kind=dbl_kind), intent(in):: &
        umco2, & !Concentration of CO2                              (ppm).
        umch4, & !Concentration of CH4                              (???).
        umn2o  !Concentration of N2o                              (???).

      real (kind=dbl_kind), intent(in), dimension(ncol):: &
        ts, &!Surface temperature                                  (K).
        amu0, &   !Cosine of solar zenith angle                    (-).
        slr, &    !Fraction of daylight                            (-).
        alvdr, &  !Visible direct surface albedo                   (-).
        alndr, &  !Near-IR direct surface albedo                   (-).
        alvdf, &  !Visible diffuse surface albedo                  (-).
        alndf, & !Near-IR diffuse surface albedo                   (-).
        cldmax !Maximum cloud fraction                             (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        ppl, &     !Layer pressure                                   (hPa).
        dp, &      !Layer thickness                                  (hPa).
        tt, &      !Temperature                                        (K).
        rmix, &    !Water vapor mixing ratio                       (kg/kg).
        cwrho, &   !Cloud water water content                     (g/m^-3).
        cwn, &     !Cloud droplet concentration                     (#/kg).
        cirho, &   !Cloud ice content                             (g/m^-3).
        o3mix, &   !Ozone mixing ratio                             (kg/kg).
        cldamt, &  !Cloud fraction                                     (-).
        b1, &      !Cloud overlap parameter                            (-).
        b2, &      !Cloud overlap parameter                            (-).
        b3, &      !Cloud overlap parameter                            (-).
        b4         !Cloud overlap parameter                            (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
        pp      !Level pressure                                   (hPa).



!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol):: &
        radvbc, &   !SFC all-sky visible direct net SW radiation   (W/m^-2).
        radvbccl, & !SFC clear-sky visible direct net SW radiation (W/m^-2).
        radvdc, &   !SFC all-sky visible direct net SW radiation   (W/m^-2).
        radvdccl, & !SFC clear-sky visible direct net SW radiation (W/m^-2).
        radnbc, &   !SFC all-sky near-ir direct net SW radiation   (W/m^-2).
        radnbccl, & !SFC clear-sky near-ir direct net SW radiation (W/m^-2).
        radndc, &   !SFC all-sky near-ir direct net SW radiation   (W/m^-2).
        radndccl    !SFC clear-sky near-ir direct net SW radiation (W/m^-2).

      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fdsw , &   !Downward SW flux                              (W/m^-2).
        fdswcl, &  !Downward clear-ksy SW flux                    (W/m^-2).
        fusw, &    !Upward SW flux                                (W/m^-2).
        fuswcl     !Upward clear-sky SW flux                      (W/m^-2).
     
! LOCAL VARIABLES:  
      integer (kind=int_kind), parameter :: &
       mb = 18, &  !Total number of spectral intervals.
       mbs = 6,  &  !Number of shortwave (SW) spectral intervals.
       mbir = 12   !Number of longwave (LW) spectral intervals.    

      integer (kind=int_kind) ::  &  
        i, &       !Horizontal index.
        l, &       !Vertical index.
        ib, &      !Index of spectral interval.
        ig      !Index of k-distribution.

      integer (kind=int_kind), dimension(ncol,nlm):: &
        ip1, &     !Used in conjunction with pressure weigthing.
        ip2     !Used in conjunction with pressure weigthing.

      real (kind=dbl_kind) :: &
        hk, &      !Weighted spectral solar constant              (W/m^-2).       .
        tmax, &    !Temperature threshold                              (K).
        eps, &     !Threshold for cloud optical properties                .
        pdist
      data eps,tmax,pdist /1.e-05,340.,2./     
      
      real (kind=dbl_kind), dimension(mbs):: &
        kg      !Nb of k-distributions per spectral intervals.  
      data kg /10,8,12,7,12,5/  

      real (kind=dbl_kind), dimension(mbs):: &
        asym_wat, & !Spectral asymmetry factor of water clouds.
        asym_ice, & !Spectral asymmetry factor of ice clouds.
        ri          !Coefficients related to Rayleigh absorption.     
      data ri / 0.9022e-5, 0.5282e-6, 0.5722e-7, &
                0.1433e-7, 0.4526e-8, 0.1529e-8 /

      real (kind=dbl_kind), dimension(mb):: &
        cnrw, &    !Real part of refractive index (Water clouds).
        cniw, &    !Imaginary part of refractive index (Water clouds).
        cnri, &    !Real part of refractive index (Ice clouds).
        cnii, &    !Imaginary part of refractive indec (Ice clouds).
        xlam       !Center of spectral band.

      real (kind=dbl_kind), dimension(ncol,mbs):: &
        asdir, &   !Spectral direct surface albedo                     (-).
        asdif      !Spectral diffuse surface albedo                    (-).

      real (kind=dbl_kind), dimension(ncol,nlm):: &
        rew, &     !Effective radius for cloud water                  (mu).
        rei, &     !Effective radius for cloud ice                    (mu).
        ttem, &    !Local temperature                                  (K).
        pkd, &     !
        tau1, &    !All-sky optical depth                              (-).
        tauclr1, & !Clear-sky optical depth                            (-).
        tau, &     !All-sky optical depth                              (-).
        tauclr, &  !Clear-sky optical depth                            (-).
        taer, &    !Aerosol optical depth                              (-).
        tray, &    !Rayley optical depth                               (-).
        tg, &      !Gases optical depth                                (-).
        tgm, &     !WV continuum optical depth                         (-).
        tcldi, &   !Ice cloud optical depth                            (-).
        tcldw, &   !Water cloud optical depth                          (-).
        wc, &      !All-sky single scattering albedo                   (-).
        wcclr, &   !Clear-sky single scattering albedo                 (-).
        waer, &    !Aerosol single scattering albedo                   (-).
        wray, &    !Rayley single scattering albedo                    (-).
        wcldi, &   !Ice cloud single scattering albedo                 (-).
        wcldw, &   !Water cloud single scattering albedo               (-).
        asym, &    !All-sky asymmetry factor                           (-).
        asyclr, &  !Clear-sky asymmetry factor                         (-).
        asyaer, &  !Aerosol asymmetry factor                           (-).
        asycldi, & !Ice cloud asymmetry factor                         (-).
        asycldw, & !Water cloud asymmetry factor                       (-).
        fwclr, &   !
        fwcld   !

      real (kind=dbl_kind), dimension(ncol,nlm+1):: &
        fdgdir, & !Spectral direct downward flux                    (W/m^2).
        fdgcldir, & !Spectral direct clear-sky downward flux        (W/m^2).
        fdgdif, & !Spectral diffuse downward flux                   (W/m^2).
        fdgcldif, & !Spectral diffuse clear-sky downward flux       (W/m^2).
        fugdif, & !Spectral diffuse upward flux                     (W/m^2).
        fugcldif   !Spectral diffuse clear-sky upward flux         (W/m^2).

!     shortwave asymmetry parameters:
!     (assumes: re=10 for water; re=30 for ice)
      data asym_wat / 0.8625, 0.8469, 0.8287, 0.8182, 0.9472, 0.7630 /
      data asym_ice / 0.8678, 0.8640, 0.8653, 0.8615, 0.9526, 0.8293 /

!---  cnrw and cniw (water clouds):
      data cnrw/1.3422,1.3281,1.3174,1.2901,1.3348,1.3700,1.3191,1.2821, &
                1.3160,1.3030,1.2739,1.2319,1.1526,1.1981,1.3542,1.4917, &
                1.5463,1.8718/
      data cniw/6.4790e-9,1.3417e-06,1.2521e-4,7.1533e-4,4.2669e-2, &
                4.3785e-3,1.3239e-2 ,1.5536e-2,5.3894e-2,3.4346e-2, &
                3.7490e-2,4.7442e-2 ,1.2059e-1,3.3546e-1,4.1698e-1, &
                4.0674e-1,3.6362e-1 ,5.2930e-1/

!--- cnri and cnii (ice clouds):
      data cnri/1.3266,1.2986,1.2826,1.2556,1.2963,1.3956, &
                1.3324,1.2960,1.3121,1.3126,1.2903,1.2295, &
                1.1803,1.5224,1.5572,1.5198,1.4993,1.7026/
      data cnii/7.0696e-9,9.1220e-7,1.2189e-4,5.7648e-4,4.3144e-2, &
                8.2935e-3,1.5540e-2,2.5594e-2,5.9424e-2,5.1511e-2, &
                4.0325e-2,4.7994e-2,2.3834e-1,3.0697e-1,1.1852e-1, &
                4.3048e-2,6.3218e-2,1.5843e-1/

!---- spectral band center:
      data xlam/0.45  ,1.0   ,1.6  ,2.2  ,3.0   ,3.75  ,4.878 ,5.556, &
                6.452 ,7.547 ,8.511,9.615,11.236,13.605,16.529,21.277, &
                29.412,71.403/

!-----------------------------------------------------------------------

!---- 0. initialize local and output arrays:

      radvbc(:)   = 0.
      radvbccl(:) = 0.
      radvdc(:)   = 0.
      radvdccl(:) = 0.
      radnbc(:)   = 0.
      radnbccl(:) = 0.
      radndc(:)   = 0.
      radndccl(:) = 0.

      fdsw(:,:)   = 0.
      fdswcl(:,:) = 0.
      fusw(:,:)   = 0.
      fuswcl(:,:) = 0.

      !rew(:,:)  = 10.
      rew(:,:) = 6.e6*(cwrho/cwn*0.0019894e-6)**(1./3.)
      rei(:,:)  = 30.

      fdgdir(:,:) = 0.0
      fdgcldir(:,:) = 0.0
      fdgdif(:,:) = 0.0
      fdgcldif(:,:) = 0.0
      fugdif(:,:) = 0.0
      fugcldif(:,:) = 0.0

      do l = 1, nlm
         do i = 1, ncol
            ttem(i,l) = min(tmax,tt(i,l))
          enddo
      enddo

      do i = 1, ncol
        asdir(i,1)   = alvdr(i)
        asdir(i,2:6) = alndr(i)
        asdif(i,1)   = alvdf(i)
        asdif(i,2:6) = alndf(i)
      enddo

!--   pressure scaling:

       call pscale(ncol,nlm,ppl,stanps,pkd,ip1,ip2)

!---- 1. loop over the mbs spectral intervals starts here:

      do ib = 1, mbs         
      
         tgm(:,:)      = 0.
         taer(:,:)     = 0.
         waer(:,:)     = 0.
         asyaer(:,:)   = 1.
!        fdswband(:,:) = 0.
!        fuswband(:,:) = 0.

!---- 1.1 rayleigh absorption:

          call rayle ( &
                    nlm, &
                    ib, &
                    pp, &
                    tray, &
                    wray)

!---- 1.2 optical properties of water and ice clouds (as in crcswr for
!        now):

         call cloudg &
                 (   ncol ,     nlm ,    mb ,    ib, &
                       pp ,      tt , cwrho ,   rew, &
                    pdist ,    cnrw ,  cniw ,  cnri, &
                     cnii ,    xlam , tcldw , wcldw, &
                  asycldw , .false. &
                 )

         call cloudg &
                 (   ncol ,   nlm   ,    mb ,    ib, &
                       pp ,    tt   , cirho ,   rei, &
                    pdist ,  cnrw   ,  cniw ,  cnri, &
                     cnii ,  xlam   , tcldi , wcldi, &
                  asycldi , .true. &
                 )

!     the asymmetry factor for water and ice clouds are fixed as
!     functions of the spectral intervals.

         do l = 1, nlm
            do i = 1, ncol
               if(cwrho(i,l).ge.eps) asycldw(i,l) = asym_wat(ib)
               if(cirho(i,l).ge.eps) asycldi(i,l) = asym_ice(ib)
            enddo
         enddo

!---- 1.3 combines single-scattering properties for gray absorption:

        call comscp1 &
                (   ncol ,     nlm ,  taer ,   tcldi, &
                   tcldw ,     tgm ,  tray ,    waer, &
                   wcldi ,   wcldw ,  wray ,  asyaer, &
                 asycldi , asycldw ,  tau1 , tauclr1, &
                    asym ,  asyclr , fwcld ,   fwclr &
                )

!---- loop over the k-probability distributions starts here:

         do ig = 1, kg(ib) 

!---- 1.4 non-gray gaseous absorption:         

            call gases &
                    ( ncol ,   nlm ,    ib ,    ig, &
                        pp ,    dp ,  ttem ,  rmix, &
                     o3mix , umco2 , umch4 , umn2o, &
                        hk ,    tg ,   pkd ,   ip1, &
                       ip2 &
                    )

!---- 1.5 combines single-scattering properties:

            call comscp2 &
                    (  ncol ,  nlm ,      tg , fwcld, &
                      fwclr , tau1 , tauclr1 ,   tau, &
                     tauclr ,   wc ,   wcclr &
                    )

!---- 1.6 two-stream approximation:
! No overlap
            call two_rt_sw &
                    (  ncol ,    nlm ,       mbs ,     ib, &
                        slr ,   amu0 ,        wc ,   asym, &
                        tau ,  asdir ,     asdif , fugdif, &
                     fdgdir , fdgdif , sel_rules &
                    )

            call two_rt_sw &
                    (    ncol ,     nlm ,       mbs ,       ib, &
                          slr ,    amu0 ,     wcclr ,   asyclr, &
                       tauclr ,   asdir ,     asdif , fugcldif, &
                     fdgcldir ,fdgcldif , sel_rules &
                    )

            fdsw(:,:)     = fdsw(:,:)   &
                            + (fdgdir(:,:)+fdgdif(:,:))*hk
            fusw(:,:)     = fusw(:,:)   + fugdif(:,:)*hk
            fdswcl(:,:)   = fdswcl(:,:) &
                            + (fdgcldir(:,:)+fdgcldif(:,:)) * hk
            fuswcl(:,:)   = fuswcl(:,:) + fugcldif(:,:)*hk

!---- 1.7 computes the surface visible and near infrared net radiation.
            select case (ib)

               case(1)
                  radvbc(:)   = radvbc(:) + fdgdir(:,nlm+1)*hk
                  radvbccl(:) = radvbccl(:) + fdgcldir(:,nlm+1)*hk
                  radvdc(:)   = radvdc(:) + fdgdif(:,nlm+1)*hk
                  radvdccl(:) = radvdccl(:) + fdgcldif(:,nlm+1)*hk

               case(2:6)
                  radnbc(:) = radnbc(:) + fdgdir(:,nlm+1)*hk
                  radnbccl(:) = radnbccl(:) + fdgcldir(:,nlm+1)*hk
                  radndc(:) = radndc(:) + fdgdif(:,nlm+1)*hk
                  radndccl(:) = radndccl(:) + fdgcldif(:,nlm+1)*hk

            end select 

         enddo ! end k-distribution
      enddo ! end spectral interval

      return
      end subroutine bugs_swr

!-----------------------------------------------------------------------





! CVS: $Id: cloudg.F,v 1.8 2006/11/16 19:54:12 norm Exp $
! CVS: $Name:  $

!-----------------------------------------------------------------------

      subroutine cloudg &
                    (  ncol ,  nlm ,    mb ,   ib, &
                         pp ,   tt , wcont ,   re, &
                      pdist , cnrw ,  cniw , cnri, &
                       cnii , xlam ,  tcld , wcld, &
                     asycld , flag &
                    )

      use bugs_kinds

      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! Cleaned up version of cloud.f from G. Stephens. Computes the cloud
! optical properties.

! This routine has been modified to include the the corrections to
! ADT for spherical particles based upon the work of David Mitchell
! DRI.  All the derivations have been carried out for the modified
! gamma distribution assuming that m=0.5 (em in program), a 
! parameter in eqn (5) of Mitchell (1994).
 
! tcld, wcld, asycld are the optical depth, single scattering albedo,
! and asymmetry parameter of cloud particles based on the use of
! ADT theory as used by Stephens et al (1990). Effective radius re  
! is input (in microns) and the water content is in g/m3.  The logical
! variable flag is .false. for water and .true. for ice.

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu.

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        flag   !If true, computes optical properties of ice clouds, of
!              of water clouds otherwise.

      integer (kind=int_kind), intent(in):: &
        ncol, &   !Length of sub-domain.
        nlm, &    !Number of layers.
        mb, &     !Total number of spectral intervals.
        ib       !Index of spectral interval.

      real (kind=dbl_kind), intent(in), dimension(mb):: &
        cnrw, &   !Real part of refractive index (Water clouds).
        cniw, &   !Imaginary part of refractive index (Water clouds).
        cnri, &   !Real part of refractive index (Ice clouds).
        cnii, &   !Imaginary part of refractive index (Ice clouds).
        xlam      !Center of spectral band.

      real (kind=dbl_kind), intent(in):: &
        pdist  !

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        tt    , & !Temperature                                        (K).
        wcont , & !Cloud water/ice content                       (g/m^-3).
        re       !Cloud effective radius                            (mu).      

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
        pp     !Level pressure                                   (hPa).
     
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm):: &
        tcld, &   !Cloud optical depth                                (-).
        wcld, &   !Cloud single scattering albedo                     (-).
        asycld !Cloud asymmetry factor                             (-).

! LOCAL VARIABLES:
      complex (kind=dbl_kind):: &
        cm,um
          
      integer (kind=int_kind):: &
        i, l

      real (kind=dbl_kind):: &
          abs ,  area ,    c0 ,  c1, &
          cnr ,   cni ,    dz , eps, &
          ext ,    f2 ,    f3 ,  no, &
           p0 ,    p1 ,    p2 ,  pi, &
           rm ,    xm ,    vm ,  rho_water
     
!-----------------------------------------------------------------------

!---- initialize local and output arrays:

      tcld(:,:)   = 0.
      wcld(:,:)   = 0.
      asycld(:,:) = 1.

!--   initialize miscellaneous constants and indices of refraction:

      pi  = acos(-1.)
      eps = 1.e-14!1.e-5
      
      if(flag) then
         !Ice case
         cnr = cnri(ib)
         cni = cnii(ib)
         rho_water = 1.0e6   !gm^-3
      else
         cnr = cnrw(ib)
         cni = cniw(ib)
         rho_water = 1.0e6   !gm^-3
      endif

!--   constants depending upon the characteristic width of the distribu
!     tion.(these may be made to vary with hydrometeor species and thus 
!     pdist could be made to depend upon level and column numbers).

!     p0    = 0.
      p0    = pdist
      p1    = p0 + 1.
      p2    = p0 + 2.
      f2    = p1 * p0
      f3    = p2 * f2

!---- calculate cloud optical properties:

      do l = 1, nlm
         do i = 1, ncol 
            if(wcont(i,l) .gt. eps) then
               dz=29.286*log(pp(i,l+1)/pp(i,l)) * tt(i,l)
               rm = re(i,l)/p2
               no = wcont(i,l) / ( (4.*pi/3.)*f3*rho_water*rm**3 )  !Particles per cubic micrometer
               area = pi*f2*no*rm**2*1.0e6                          !The factor converts inverse micrometers,
                                                                    !i.e., micrometer^2/micrometer^3, to inverse meters
               c0 = 2.*area
               c1 = c0/f2
               xm = 2.*pi*rm/xlam(ib)
               cm = cmplx(cnr,-cni)

               if (ib .eq. 1 .or. ib .eq. 2) then
                  !For band 1 (0.5 um) only compute the extinction.
                  um = 2.*xm*(cnr-1.)*cmplx(0.d0,1.d0)
                  ext = c0 + 2.*c1*real(p0/(um*(um+1.)**p1) &
                          + 1./(um**2*(um+1.)**p0)-1./um**2)
                  tcld(i,l) = ext*dz
                  wcld(i,l) = 0.999999
                  asycld(i,l) = 0.85
               else
                  !Compute both extinction and absorption coefficients for all other bands.
                  um = 2.*xm*(cm-1.)*cmplx(0.d0,1.d0)
                  ext = c0 + 2.*c1*real( p0/(um*(um+1.)**p1) &
                        + 1./(um**2*(um+1.)**p0)-1./um**2)
                  vm = 4.*xm*cni
                  abs = area + c1*sngl( p0/(vm*(vm+1.)**p1) &
                          + 1./(vm**2*(vm+1.)**p0) - 1./vm**2 )
!                 abs = area + c1*sngl( &
!                    p0/(dble(vm)*(dble(vm)+1.)**dble(p1)) &
!                  + 1./(dble(vm)**2*(dble(vm)+1.)**dble(p0)) &
!                  - 1./dble(vm)**2)
                 tcld(i,l) = ext*dz

                 if (ext.lt.abs) ext = abs
                 wcld(i,l) = (ext-abs)/ext

                 if(wcld(i,l) .lt. 0.) then
                      print *,wcld(i,l),ext,abs,wcont(i,l)
                      print *,pp(i,l),pp(i,l+1)
                      print *,tt(i,l)
                      print *,re(i,l)
                      stop
                 endif
               endif
!               if(flag) then
!        if(ib==2) print*,ib,tcld(i,l),wcld(i,l),rm,xlam(ib),re(i,l),p2
!        if(ib==2) print*,um,xm,cm,ext,c0,c1,p0,p1,area,vm
!                else
!               !if(ib<=8) print*,ib,tcld(i,l),wcld(i,l)
!               endif
               asycld(i,l)=0.85  !default, overridden in bugs_swr(), bugs_lwr()
             endif
         enddo
      enddo

      return
      end subroutine cloudg
!-----------------------------------------------------------------------


! CVS:  $Id: comscp1.F,v 1.3 2001/04/30 08:48:56 norm Exp $
! CVS:  $Name:  $

!-----------------------------------------------------------------------

      subroutine comscp1 &
                    (   ncol ,     nlm ,   taer ,  tcldi, &
                       tcldw ,     tgm ,   tray ,   waer, &
                       wcldi ,   wcldw ,   wray , asyaer, &
                     asycldi , asycldw , tccld1 , tcclr1, &
                      asycld ,  asyclr ,  fwcld ,  fwclr &
                    )

      use bugs_kinds

      implicit none
             
!-----------------------------------------------------------------------
! REFERENCES:
! comscp1 combines single scattering properties due to Rayleigh absorp
! tion, aerosols, water continuum, gray gaseous absorption, ice crystals
! and water droplets. 
! Laura D. Fowler/slikrock (08-12-97).

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. In this subroutine, all the arrays are defined as 
! local arrays in BUGSswr.

!     INPUT ARGUMENTS:
!     ----------------
      integer (kind=int_kind), intent(in):: &
        ncol, &    !Length of sub-domain..
        nlm        !Number of layers.
     
      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        asyaer, &  !Asymmetry factor of aerosols                      (-).
        asycldi, & !Asymmetry factor of ice clouds                    (-).
        asycldw, & !Asymmetry factor of water clouds                  (-).
        taer, &    !Optical depth of aerosols                         (-).
        tcldi, &   !Optical depth of ice clouds                       (-).
        tcldw, &   !Optical depth of water clouds                     (-).
        tgm, &     !Optical depth of water vapor continuum            (-).
        tray, &    !Optical depth due to Rayleigh absorption          (-).
        waer, &    !Single scattering albedo of aerosols              (-).
        wcldi, &   !Single scattering albedo of ice clouds            (-).
        wcldw, &   !Single scattering albedo of water clouds          (-).
        wray       !Single scattering albedo due to Rayleigh
!                   absorption                                        (-).

!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm):: &
        asyclr, &  !Clear-sky asymmetry factor                        (-).
        asycld, &  !All-sky asymmetry factor                          (-).
        tcclr1, &  !Clear-sky optical depth                           (-).
        tccld1, &  !All-sky optical depth                             (-).
        fwclr, &   !Total clear-sky single-scattering albedo          (-).
        fwcld   !Total cloudy single-scattering albedo             (-).

! LOCAL LIST VARIABLES:    

      integer (kind=int_kind):: &
        i, &       !Horizontal index.
        l          !Vertical index.

      real (kind=dbl_kind):: &
        wwray,wwaer,wwcldi,wwcldw
     
!-----------------------------------------------------------------------

      do l = 1, nlm
         do i = 1, ncol

            tcclr1(i,l) = tgm(i,l) + tray(i,l) + taer(i,l)
            tccld1(i,l) = tcclr1(i,l)+ tcldi(i,l) + tcldw(i,l)

            wwray  = wray(i,l)*tray(i,l)
            wwaer  = waer(i,l)*taer(i,l)
            wwcldi = wcldi(i,l)*tcldi(i,l)
            wwcldw = wcldw(i,l)*tcldw(i,l)

            fwclr(i,l)  = wwray+wwaer
            fwcld(i,l)  = fwclr(i,l)+wwcldi+wwcldw
 
            if(fwclr(i,l).gt.1.e-10) then
               asyclr(i,l) = (asyaer(i,l)*wwaer)/fwclr(i,l)
            else
               asyclr(i,l) = 1.
            endif
 
            if(fwcld(i,l).gt.1.e-10) then
               asycld(i,l) = (asyaer(i,l)*wwaer+asycldi(i,l)*wwcldi &
                           +  asycldw(i,l)*wwcldw) &
                           / fwcld(i,l)
            else
               asycld(i,l) = 1.
            endif

     	 enddo
      enddo

      return
      end subroutine comscp1

!-----------------------------------------------------------------------


! CVS:  $Id: comscp2.F,v 1.3 2001/04/30 08:47:14 norm Exp $
! CVS:  $Name:  $


!-----------------------------------------------------------------------

      subroutine comscp2 &
                    ( ncol ,    nlm ,     tg , fwcld, &
                     fwclr , tccld1 , tcclr1 , tccld, &
                     tcclr ,  wccld ,  wcclr &
                    )

      use bugs_kinds
      
      implicit none
             
!-----------------------------------------------------------------------
! REFERENCES:
! comscp2 combines the single scattering properties computed in comscp1
! to the single scattering properties due to non-gray absorption.
! Laura D. Fowler (slikrock. 08-12-97).

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. In this subroutine, all the arrays are defined as 
! local arrays in BUGSswr.

!     INPUT ARGUMENTS:
!     ----------------
      integer (kind=int_kind), intent(in):: &
        ncol, &   !Length of sub-domain..
        nlm       !Number of layers.
     
      real (kind=dbl_kind), dimension(ncol,nlm):: &
        tg, &     !Optical depth of non-gray gases                    (-).
        fwclr, &  !Clear-sky single scattering albedo from comscp1    (-).
        fwcld     !

!     INPUT/OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), dimension(ncol,nlm):: &
        tcclr1, & !Clear-sky optical depth                            (-).
        tccld1, & !All-sky optical depth                              (-).
        tcclr, &  !Clear-sky optical depth                            (-).
        tccld     !All-sky optical depth                              (-).

!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm):: &
        wcclr, &  !Clear-sky single scattering albedo                 (-).
        wccld     !All-sky single scattering albedo                   (-).

! LOCAL LIST VARIABLES:    

      integer (kind=int_kind):: &
        i, &      !Horizontal index.
        l      !Vertical index.
     
!-----------------------------------------------------------------------

      do l = 1, nlm
         do i = 1, ncol

            tcclr(i,l) = tcclr1(i,l) + tg(i,l)
            tccld(i,l) = tccld1(i,l) + tg(i,l)
 
            if(tcclr(i,l).gt.0.) then
              wcclr(i,l) = fwclr(i,l)/tcclr(i,l)
            else
              wcclr(i,l) = 0.
            endif
            wcclr(i,l)  = min(.999999_dbl_kind,wcclr(i,l))
 
            if(tccld(i,l).gt.0.) then
               wccld(i,l) = fwcld(i,l)/tccld(i,l)
            else
               wccld(i,l) = 0.
            endif
            wccld(i,l) = min(.999999_dbl_kind,wccld(i,l))

     	   enddo
	enddo

  return
	end subroutine comscp2

!-----------------------------------------------------------------------


! CVS:  $Id: two_rt_lw.F,v 1.7 2003/11/11 21:55:13 norm Exp $
! CVS:  $Name:  $

!-----------------------------------------------------------------------
 
      subroutine two_rt_lw &
                    (     ncol , nlm,  mbs , mbir, &
                            ib ,  wc, asym ,  tau, &
                            es ,  bf,   fu ,   fd, &
                     sel_rules &
                    )
 
      use bugs_kinds



      implicit none
!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_lw replaces two_rt and add written by G. Stephens. two_rt_lw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        sel_rules

      integer (kind=int_kind), intent(in):: &
        ncol, &  !Length of sub-domain.
        nlm, &   !Number of layers.
        mbs, &   !Number of SW spectral intervals.
        mbir, &  !Number of IR spectral intervals.
        ib       !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol,mbir):: &
        es    !Spectral surface emissivity                         (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        wc, &    !Single scattering albedo                            (-).
        asym, &  !Asymmetry factor                                    (-).
        tau   !Optical depth                                       (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
        bf    !Planck function                                 (W/m^2).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fd, &    !Spectral downward flux                          (W/m^2).
        fu       !Spectral upward flux                            (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind):: &
        i, &     !Horizontal index.
        l, &     !Vertical index.
        ibms  !Index of spectral interval.
 
      real (kind=dbl_kind), dimension(nlm):: &
        rr, &    !
        tr, &    !
        sigu, &  !
        sigd  !

      real (kind=dbl_kind):: &
            aa ,    bb , beta0 ,     cc, &
        diffac , denom ,  fact , eggtau, &
         ggtau, &
         kappa ,   oms ,  prop ,r, rinf, &
             t ,  taus
      data diffac /2./

      real (kind=dbl_kind), dimension(nlm):: &
        td , vu , exptau

      real (kind=dbl_kind), dimension(nlm+1):: &
        re, vd
 
! SELECTION RULE VARIABLES

      logical (kind=log_kind):: &
        fail

      real (kind=dbl_kind):: &
        tausthresh, &
        wcthresh, &
        tauscat

      data tausthresh / 0.001 /
      data wcthresh   / 0.975 /
 
!----------------------------------------------------------------------

      fd(:,1) = 0.
 
      ibms = ib - mbs
       
      do 1000 i = 1, ncol
      
         if(sel_rules) then

            fail = .false.
            tauscat = 0.0
            do l = nlm, 1, -1
               if (wc(i,l).gt.wcthresh) fail = .true.
               tauscat = tauscat + wc(i,l)*tau(i,l)
            enddo
            if (fail.and.tauscat.ge.tausthresh) goto 2000
 
!>> BEGIN SELECTION RULES <<
!           print *,'selection rules'
            do l = 1, nlm
               exptau(l) = exp(-2.0*tau(i,l))
               if(tau(i,l) .lt. 0.8e-2) then
                  sigu(l) = (bf(i,l)+bf(i,l+1))*tau(i,l)
                  sigd(l) = sigu(l)
               else
                  prop = (1.-exptau(l))/tau(i,l)
                  aa = 2.-prop
                  bb = -2.*exptau(l)+prop
                  cc = 0.5
                  sigu(l) = (aa*bf(i,l)+bb*bf(i,l+1))*cc
                  sigd(l) = (bb*bf(i,l)+aa*bf(i,l+1))*cc
               endif
               fd(i,l+1) = sigd(l) + exptau(l) * fd(i,l)
            enddo
 
            fu(i,nlm+1) = bf(i,nlm+1)*es(i,ibms)
!    &                  + fd(i,nlm+1)*(1.0-es(i,ibms))
 
            do l = nlm , 1, -1
               fu(i,l) = sigu(l) + exptau(l) * fu(i,l+1)
            enddo
 
            cycle

         endif

!>> END SELECTION RULES <<
 
!>> BEGIN FULL CALCULATION <<
2000  re(1) = 0.
      vd(1) = 0.
!     print *,'full up calculation'
 
      do l = 1, nlm
         fact = asym(i,l)*asym(i,l)
         oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
         taus   = (1.-fact*wc(i,l))*tau(i,l)
 
         beta0 = (4.+asym(i,l))/(8.*(1.+asym(i,l)))
         t = diffac*(1.-oms*(1.-beta0))     !-0.25
         r = diffac*oms*beta0               !-0.25
         kappa = sqrt(t**2-r**2)
         rinf   = r/(kappa+t)
         ggtau  = kappa*taus
         eggtau = exp(-ggtau)
         denom  = (1.-rinf**2*eggtau**2)
         tr(l) = (1.-rinf**2)*eggtau/denom
         rr(l) = rinf*(1.-eggtau**2)/denom

         if(taus .lt. 0.8e-2) then
            sigu(l) = 0.5*diffac*(bf(i,l)+bf(i,l+1))*taus
            sigd(l) = sigu(l)
         else
            aa =  (t+r)*(1.-rr(l))-(1.+rr(l)-tr(l))/taus
            bb = -(t+r)*tr(l)+(1.+rr(l)-tr(l))/taus
            cc = diffac*(1.-oms)/kappa**2
            sigu(l) = cc*(aa*bf(i,l)+bb*bf(i,l+1))
            sigd(l) = cc*(bb*bf(i,l)+aa*bf(i,l+1))
         endif
      enddo
 
!---- 1. do adding, going from top down:

        do l = 1, nlm
           prop = 1. / (1. - re(l)*rr(l))
           re(l+1) = rr(l) + tr(l)**2*re(l)*prop
           vd(l+1) = sigd(l) + (tr(l)*vd(l) &
                   + tr(l)*re(l)*sigu(l))*prop
           vu(l)   = (rr(l)*vd(l) + sigu(l))*prop
           td(l)   = prop
        enddo
 
!---- 2. calculate fluxes going up through the layers:

        fu(i,nlm+1) = es(i,ibms)*bf(i,nlm+1)
 
        do l = nlm+1, 2, -1
           fd(i,l)   = re(l)*fu(i,l) + vd(l)
           fu(i,l-1) = tr(l-1)*fu(i,l)*td(l-1) + vu(l-1)
        enddo

!>> END FULL CALCULATION <<
 
 1000 continue
 
      return
      end subroutine two_rt_lw
 
!------------------------------------------------------------------------


! CVS:  $Id: two_rt_lw_iter.F,v 1.2 2003/11/11 21:55:13 norm Exp $
! CVS:  $Name:  $

!-----------------------------------------------------------------------
 
      subroutine two_rt_lw_iter &
                    ( &
                          ncol ,    nlm ,  mbs ,   mbir, &
                            ib , cldamt ,   wc ,  wcclr, &
                          asym , asyclr ,  tau , tauclr, &
                            es ,     bf ,   fu ,     fd, &
                     sel_rules ,     b1 ,   b2 ,     b3, &
                            b4 &
                    )

      use bugs_kinds



      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_lw replaces two_rt and add written by G. Stephens. two_rt_lw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        sel_rules

      integer (kind=int_kind), intent(in):: &
        ncol, &   !Length of sub-domain.
        nlm, &    !Number of layers.
        mbs, &    !Number of SW spectral intervals.
        mbir, &   !Number of IR spectral intervals.
        ib     !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol,mbir):: &
        es    !Spectral surface emissivity                         (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        cldamt, & !Cloud fraction                                     (-).
        wc, &     !All sky single scattering albedo                   (-).
        wcclr, &  !Clear sky single scattering albedo                 (-).
        asym, &   !All sky asymmetry factor                           (-).
        asyclr, & !Clear sky asymmetry factor                         (-).
        tau, &    !All sky optical depth                              (-).
        tauclr, & !Clear sky optical depth                            (-).
        b1, &     !Cloud overlap parameter                            (-).
        b2, &     !Cloud overlap parameter                            (-).
        b3, &     !Cloud overlap parameter                            (-).
        b4        !Cloud overlap parameter                            (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
        bf    !Planck function                                 (W/m^2).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fd, &     !Spectral downward flux                         (W/m^2).
        fu     !Spectral upward flux                           (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind):: &
        i, &      !Horizontal index.
        l, &      !Vertical index.
        ibms, &   !Index of spectral interval.
        j, &
        nsr, &
        nsmx, &
        n, &
        ii, &
        jj, &
        kk, &
        ir, &
        iter
      integer (kind=int_kind), dimension(16*nlm-6):: &
        idc

      integer (kind=int_kind), dimension(4*nlm+2):: &
        nir

      real (kind=dbl_kind), dimension(nlm):: &
        rrcld, &    !All sky global reflection                     (-).
        rrclr, &    !Clear sky global reflection                   (-).
        trcld, &    !All sky global transmission                   (-).
        trclr, &    !Clear sky global transmission                 (-).
        sigucld, &  !All sky upwelling source                      (-).
        siguclr, &  !Clear sky upwelling source                    (-).
        sigdcld, &  !All sky downwelling source                    (-).
        sigdclr, &  !Clear sky downwelling source                  (-).
        exptau   !

      real (kind=dbl_kind), dimension(4*nlm+2):: &
        b, &
        fvc, &
        error !, &
!        old_fvc

      real (kind=dbl_kind), dimension(16*nlm-6):: &
        smx

      real (kind=dbl_kind):: &
            aa ,    bb , beta0 ,     cc, &
        diffac , denom ,  fact , eggtau, &
         ggtau, &
         kappa ,   oms ,  prop ,r, rinf, &
             t ,  taus , omega
      data diffac /2./
 
! SELECTION RULE VARIABLES

      logical (kind=log_kind):: &
        fail

      real (kind=dbl_kind):: &
        tausthresh, &
        wcthresh, &
        tauscat
      data tausthresh / 0.001 /
      data wcthresh   / 0.975 /
 
!----------------------------------------------------------------------
      nsr = 4*nlm + 2
      nsmx = 16*nlm - 6

      fd(:,1) = 0.
 
      ibms = ib - mbs
 
      do i = 1, ncol
 
!---- 2. DO LONGWAVE:
         do l = 1, nlm
            !ALL SKY
            fact = asym(i,l)*asym(i,l)
            oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
            taus   = (1.-fact*wc(i,l))*tau(i,l)
 
            beta0 = (4.+asym(i,l))/(8.*(1.+asym(i,l)))
            t = diffac*(1.-oms*(1.-beta0))     !-0.25
            r = diffac*oms*beta0               !-0.25
            kappa = sqrt(t**2-r**2)
            rinf   = r/(kappa+t)
            ggtau  = kappa*taus

            eggtau = exp(-ggtau)

            denom  = (1.-rinf**2*eggtau**2)
            trcld(l) = (1.-rinf**2)*eggtau/denom
            rrcld(l) = rinf*(1.-eggtau**2)/denom

            if(taus .lt. 0.8e-2) then
               sigucld(l) = cldamt(i,l)*0.5*diffac*(bf(i,l)+ &
                            bf(i,l+1))*taus
               sigdcld(l) = cldamt(i,l)*sigucld(l)
            else
               aa =  (t+r)*(1.-rrcld(l))-(1.+rrcld(l)-trcld(l))/taus
               bb = -(t+r)*trcld(l)+(1.+rrcld(l)-trcld(l))/taus
               cc = diffac*(1.-oms)/kappa**2
               sigucld(l) = cldamt(i,l)*cc*(aa*bf(i,l)+bb*bf(i,l+1))
               sigdcld(l) = cldamt(i,l)*cc*(bb*bf(i,l)+aa*bf(i,l+1))
            endif

            !CLEAR SKY
            fact = asyclr(i,l)*asyclr(i,l)
            oms  = ((1.-fact)*wcclr(i,l))/(1.-fact*wcclr(i,l))
            taus   = (1.-fact*wcclr(i,l))*tauclr(i,l)
 
            beta0 = (4.+asyclr(i,l))/(8.*(1.+asyclr(i,l)))
            t = diffac*(1.-oms*(1.-beta0))     !-0.25
            r = diffac*oms*beta0               !-0.25
            kappa = sqrt(t**2-r**2)
            rinf   = r/(kappa+t)
            ggtau  = kappa*taus

            eggtau = exp(-ggtau)

            denom  = (1.-rinf**2*eggtau**2)
            trclr(l) = (1.-rinf**2)*eggtau/denom
            rrclr(l) = rinf*(1.-eggtau**2)/denom

            if(taus .lt. 0.8e-2) then
               siguclr(l) = (1.0-cldamt(i,l))*0.5*diffac*(bf(i,l)+ &
                            bf(i,l+1))*taus
               sigdclr(l) = (1.0-cldamt(i,l))*siguclr(l)
            else
               aa =  (t+r)*(1.-rrclr(l))-(1.+rrclr(l)-trclr(l))/taus
               bb = -(t+r)*trclr(l)+(1.+rrclr(l)-trclr(l))/taus
               cc = diffac*(1.-oms)/kappa**2
               siguclr(l) = (1.0-cldamt(i,l))*cc*(aa*bf(i,l)+ &
                            bb*bf(i,l+1))
               sigdclr(l) = (1.0-cldamt(i,l))*cc*(bb*bf(i,l)+ &
                            aa*bf(i,l+1))
            endif

         enddo

 
!---- 1. LOAD SMX VECTOR
        nir(:) = 4

        idc(1) = 5
        idc(2) = 6
        smx(1) = -trcld(1) * b4(i,1)
        smx(2) = -trcld(1) * (1.-b2(i,1))
        nir(1) = 2

        idc(3) = 5
        idc(4) = 6
        smx(3) = -trclr(1) * (1.-b4(i,1))
        smx(4) = -trclr(1) * b2(i,1)
        nir(2) = 2

        idc(5) = 5
        idc(6) = 6
        smx(5) = -rrcld(1) * b4(i,1)
        smx(6) = -rrcld(1) * (1.-b2(i,1))
        nir(3) = 2

        idc(7) = 5
        idc(8) = 6
        smx(7) = -rrclr(1) * (1.-b4(i,1))
        smx(8) = -rrclr(1) * b2(i,1)
        nir(4) = 2

        do l = 1,nlm-1
          n = (l-1)*16 + 9
          ir = 4*l

          idc(n)   = ir-1
          idc(n+1) = ir
          idc(n+2) = ir+5
          idc(n+3) = ir+6
          smx(n)   = -rrcld(l+1) * b3(i,l+1)
          smx(n+1) = -rrcld(l+1) * (1.-b1(i,l+1))
          smx(n+2) = -trcld(l+1) * b4(i,l+1)
          smx(n+3) = -trcld(l+1) * (1.-b2(i,l+1))

          idc(n+4) = ir-1
          idc(n+5) = ir
          idc(n+6) = ir+5
          idc(n+7) = ir+6
          smx(n+4) = -rrclr(l+1) * (1.-b3(i,l+1))
          smx(n+5) = -rrclr(l+1) * b1(i,l+1)
          smx(n+6) = -trclr(l+1) * (1.-b4(i,l+1))
          smx(n+7) = -trclr(l+1) * b2(i,l+1)

          idc(n+8) = ir-1
          idc(n+9) = ir
          idc(n+10) = ir+5
          idc(n+11) = ir+6
          smx(n+8) = -trcld(l+1) * b3(i,l+1)
          smx(n+9) = -trcld(l+1) * (1.-b1(i,l+1))
          smx(n+10) = -rrcld(l+1) * b4(i,l+1)
          smx(n+11) = -rrcld(l+1) * (1.-b2(i,l+1))

          idc(n+12) = ir-1
          idc(n+13) = ir
          idc(n+14) = ir+5
          idc(n+15) = ir+6
          smx(n+12) = -trclr(l+1) * (1.-b3(i,l+1))
          smx(n+13) = -trclr(l+1) * b1(i,l+1)
          smx(n+14) = -rrclr(l+1) * (1.-b4(i,l+1))
          smx(n+15) = -rrclr(l+1) * b2(i,l+1)
        enddo

        ir = 4*nlm

        idc(16*nlm-7) = 4*nlm-1
        idc(16*nlm-6) = 4*nlm
        smx(16*nlm-7) = 0.0 !1.-es(i,ibms)
        smx(16*nlm-6) = 0.0 !1.-es(i,ibms)
        nir(ir+1) = 1
        nir(ir+2) = 1

        b(:) = 0.0
        do l = 1,nlm
          b(l*4-3) = sigucld(l)
          b(l*4-2) = siguclr(l)
          b(l*4-1) = sigdcld(l)
          b(l*4)   = sigdclr(l)
        enddo
        b(nlm*4+1) = cldamt(i,nlm)*es(i,ibms)*bf(i,nlm+1)
        b(nlm*4+2) = (1.-cldamt(i,nlm))*es(i,ibms)*bf(i,nlm+1)
        



!-------------- GAUSS SEIDEL W/ OVERRELAXATION ------------------
        omega = 1.0

        fvc(:) = b(:)

        do iter=1,200
           kk = 1
           do ii=1,nsr
             t = 0.0
             do j=1,nir(ii)
               jj = idc(kk)
               t = t + smx(kk) * fvc(jj)
               kk = kk + 1
             enddo
             t = t + fvc(ii)
             fvc(ii) = fvc(ii) + omega * (b(ii)-t)
             error(ii) = b(ii) - t
           enddo

           if (maxval(abs(error)) .le. 0.05) then
             !print *,omega,iter,' iterations'
             exit
           endif
        enddo


!---- 3. SUM CLEAR AND CLOUDY FLUXES
        do l = 1,nlm+1
          fu(i,l) = fvc(l*4-3)+fvc(l*4-2)
        enddo
        do l = 1,nlm
          fd(i,l+1) = fvc(l*4-1)+fvc(l*4)
        enddo
 
      end do  !i=1,ncol
 
      return
      end subroutine two_rt_lw_iter


! CVS:  $Id: two_rt_lw_sel.F,v 1.3 2003/11/11 21:55:13 norm Exp $
! CVS:  $Name:  $


!-----------------------------------------------------------------------
 
      subroutine two_rt_lw_sel &
                       ( &
                          ncol , nlm ,  mbs , mbir , ib, &
                        tauclr ,  es ,   bf ,   fu , fd &
                       )
 
      use bugs_kinds




      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_lw replaces two_rt and add written by G. Stephens. two_rt_lw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      integer (kind=int_kind), intent(in):: &
        ncol, &   !Length of sub-domain.
        nlm, &    !Number of layers.
        mbs, &    !Number of SW spectral intervals.
        mbir, &   !Number of IR spectral intervals.
        ib        !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol,mbir):: &
        es     ! Spectral surface emissivity                       (-).
      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        tauclr !Optical depth                                      (-).
      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1):: &
        bf     !Planck function                                    (-).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fd, &    !Spectral downward flux                          (W/m^2).
        fu    !Spectral upward flux                            (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind):: &
        i, &     !Horizontal index.
        l, &     !Vertical index.
        ibms     !Index of spectral interval.
 
      real (kind=dbl_kind), dimension(nlm):: &
        sigu, &
        sigd, &
        exptau

      real (kind=dbl_kind):: &
        aa, &
        bb, &
        cc, &
        prop

!----------------------------------------------------------------------

 
      ibms = ib - mbs
 
      do 1000 i = 1, ncol
 
          !TOA/BOA initializations
          fu(i,nlm+1) = bf(i,nlm+1)*es(i,ibms)
          fd(i,1) = 0.

          do l=1,nlm
            exptau(l) = exp(-2.*tauclr(i,l))
            if(tauclr(i,l) .lt. 0.8e-2) then
              sigu(l) = (bf(i,l)+bf(i,l+1))*tauclr(i,l)
              sigd(l) = sigu(l)
            else
              prop = (1.-exptau(l))/tauclr(i,l)
              aa = 2.-prop
              bb = -2.*exptau(l)+prop
              cc = 0.5
              sigu(l) = (aa*bf(i,l)+bb*bf(i,l+1))*cc
              sigd(l) = (bb*bf(i,l)+aa*bf(i,l+1))*cc
            endif
            fd(i,l+1) = sigd(l) + exptau(l) * fd(i,l)
          enddo
 
          do l=nlm,1,-1
            fu(i,l) = sigu(l) + exptau(l) * fu(i,l+1)
          enddo

1000  continue

      return
      end subroutine two_rt_lw_sel


! CVS:  $Id: two_rt_sw_bs.F,v 1.3 2003/11/11 21:55:13 norm Exp $
! CVS:  $Name:  $

!-----------------------------------------------------------------------
 
      subroutine two_rt_sw_bs &
                     ( ncol ,      nlm ,   mbs ,     ib, &
                        slr ,     amu0 ,    wc ,  wcclr, &
                       asym ,   asyclr ,   tau , tauclr, &
                      asdir ,    asdif , fudif ,  fddir, &
                      fddif ,sel_rules ,    b1 ,     b2, &
                         b3 ,       b4 &
                     )
 
 
      use bugs_kinds
      use bandsolve 



      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_sw replaces two_rt and add written by G. Stephens. two_rt_sw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        sel_rules

      integer (kind=int_kind), intent(in):: &
        ncol, &  !Length of sub-domain.
        nlm, &   !Number of layers.
        mbs, &   !Number of SW spectral intervals.
        ib       !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol):: &
        slr, &   !Fraction of daylight                             (-).
        amu0     !Cosine of the solar zenith angle                 (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,mbs):: &
        asdir, & !Spectral direct surface albedo                   (-).
        asdif    !Spectral diffuse surface albedo                  (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        wc, &     !Single scattering albedo                           (-).
        wcclr, &  !Single scattering albedo                           (-).
        asym, &   !Asymmetry factor                                   (-).
        asyclr, & !Asymmetry factor                                   (-).
        tau, &    !Optical depth                                      (-).
        tauclr, & !Optical depth                                      (-).
        b1, &     !Cloud overlap parameter                            (-).
        b2, &     !Cloud overlap parameter                            (-).
        b3, &     !Cloud overlap parameter                            (-).
        b4        !Cloud overlap parameter                            (-).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fddir, & !Spectral direct downward flux                   (W/m^2).
        fddif, & !Spectral diffuse downward flux                  (W/m^2).
        fudif    !Spectral diffuse upward flux                    (W/m^2).
 
! LOCAL VARIABLES:
! ----------------
      integer (kind=int_kind):: &
        i, &     ! Horizontal index.
        l     ! Vertical index.

      integer (kind=int_kind), dimension(nlm*4+2):: &
        indx    !Vector used by bandec and banbks
 
      real (kind=dbl_kind), dimension(nlm):: &
        rrcld, &   !All sky global reflection
        rrclr, &   !Clear sky global reflection
        trcld, &   !All sky global transmission
        trclr, &   !Clear sky global transmission
        sigucld, & !All sky upwelling source
        siguclr, & !Clear sky upwelling source
        sigdcld, & !All sky downwelling source
        sigdclr    !Clear sky downwelling source

      real (kind=dbl_kind), dimension(nlm*4+2,11):: &
        a, &       !Diagonal matrix for bandec and banbks
        al      !Matrix used by bandec and banbks

      real (kind=dbl_kind), dimension(nlm*4+2):: &
        b       !Vector of sources (for A*x = b)

      real (kind=dbl_kind):: &
        exptaucld, &
        exptauclr, &
        directcld(nlm+1), &
        directclr(nlm+1), &
        d

      real (kind=dbl_kind):: &
            aa ,    bb ,   cc , denom, &
        eggtau ,   eps ,   g3 ,    g4, &
         ggtau , kappa ,    r ,  rinf, &
             t ,   oms , taus ,  fact, &
           asy
      data eps /1.e-02/

      real (kind=dbl_kind), dimension(nlm+1):: &    
        re , vd
 
! SELECTION RULE VARIABLES:
! -------------------------
      logical (kind=log_kind):: &
        fail

      real (kind=dbl_kind):: &
        tausthresh, &
        wcthresh, &
        tauscat

!#ifdef usenewexp
!      real(kind=dbl_kind), external :: exp
!#endif

!     data tausthresh / 0.001 /
!     data wcthresh   / 0.975 /
      data tausthresh / 0.01 /
      data wcthresh   / 0.98 /
 
!----------------------------------------------------------------------
 
      do 1000 i = 1, ncol
 
!         if (sel_rules) then
!            fail = .false.
!            tauscat = 0.0
!            do l = nlm, 1, -1
!               if (wc(i,l).gt.wcthresh) fail = .true.
!               tauscat = tauscat + wc(i,l)*tau(i,l)
!            enddo
!            if (fail.and.tauscat.ge.tausthresh) goto 2000
! 
!!>> BEGIN SELECTION RULES <<
!!            print *,'selection rules'
!            fddir(i,1) = amu0(i)*slr(i)
!            fddif(i,:) = 0.0
!            do l=1,nlm
!               fddir(i,l+1) = fddir(i,l) * exp(-1.*tau(i,l)/amu0(i))
!            enddo
!
!            fudif(i,nlm+1) = fddir(i,nlm+1) * asdir(i,ib)
! 
!            do l=nlm,1,-1
!               fudif(i,l) = fudif(i,l+1) * exp(-2.*tau(i,l))
!            enddo
! 
!            cycle
!         endif
!!>> END SELECTION RULES <<

!>> BEGIN FULL CALCULATION <<
2000     directcld(1) = 0.
         directclr(1) = 1.
         re(1) = 0.
         vd(1) = 0.
!         print *,'full up calculation'
 
!---- 1. DO SHORTWAVE:
         !ALL SKY
         do l = 1, nlm
            fact = asym(i,l)*asym(i,l)
            oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
            taus   = (1.-fact*wc(i,l))*tau(i,l)
            asy = asym(i,l)/(1.+asym(i,l))
 
            exptaucld = exp(-taus/amu0(i))
 
!--- local coefficients:  delta-eddington
            t     = 0.25 * (7. - oms*(4.+3.*asy))
            r     = -0.25 * (1. - oms*(4.-3.*asy))
            kappa  = sqrt(t**2-r**2)
            rinf   = r/(kappa+t)
            ggtau  = kappa*taus

            eggtau = exp(-ggtau)
            denom  = (1.-rinf**2*eggtau**2)
            trcld(l) = (1.-rinf**2)*eggtau/denom
            rrcld(l) = rinf*(1.-eggtau**2)/denom
 
            if(abs(kappa**2-1./amu0(i)**2) .lt. eps) then
               fact = 1./eps
            else
               fact = 1./(kappa**2-1./amu0(i)**2)
            endif
            !cc = oms*slr(i)*fact
            cc = oms*fact
            g3 = 0.5-0.75*asy*amu0(i)
            g4 = 1.-g3
            aa = g3*(t-1./amu0(i))+g4*r
            bb = g4*(t+1./amu0(i))+g3*r
 
            sigucld(l)=cc*((aa-rrcld(l)*bb)-aa*trcld(l)*exptaucld)* &
                    (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))
            sigdcld(l)=cc*(-bb*trcld(l)+(bb-rrcld(l)*aa)*exptaucld)* &
                    (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))

            !CLEAR SKY
            fact = asyclr(i,l)*asyclr(i,l)
            oms  = ((1.-fact)*wcclr(i,l))/(1.-fact*wcclr(i,l))
            taus   = (1.-fact*wcclr(i,l))*tauclr(i,l)
            asy = asyclr(i,l)/(1.+asyclr(i,l))
 
            exptauclr = exp(-taus/amu0(i))

 
!--- local coefficients:  delta-eddington
            t      = 0.25 * (7. - oms*(4.+3.*asy))
            r      = -0.25 * (1. - oms*(4.-3.*asy))
            kappa  = sqrt(t**2-r**2)
            rinf   = r/(kappa+t)
            ggtau  = kappa*taus

            eggtau = exp(-ggtau)

            denom  = (1.-rinf**2*eggtau**2)
            trclr(l) = (1.-rinf**2)*eggtau/denom
            rrclr(l) = rinf*(1.-eggtau**2)/denom
 
            if(abs(kappa**2-1./amu0(i)**2) .lt. eps) then
               fact = 1./eps
            else
               fact = 1./(kappa**2-1./amu0(i)**2)
            endif
            !cc = oms*slr(i)*fact
            cc = oms*fact
            g3 = 0.5-0.75*asy*amu0(i)
            g4 = 1.-g3
            aa = g3*(t-1./amu0(i))+g4*r
            bb = g4*(t+1./amu0(i))+g3*r

            siguclr(l)=cc*((aa-rrclr(l)*bb)-aa*trclr(l)*exptauclr)* &
                    (b1(i,l)*directclr(l)+(1.-b3(i,l))*directcld(l))
            sigdclr(l)=cc*(-bb*trclr(l)+(bb-rrclr(l)*aa)*exptauclr)* &
                    (b1(i,l)*directclr(l)+(1.-b3(i,l))*directcld(l))
 
            directclr(l+1) = exptauclr * &
                ((1.-b3(i,l))*directcld(l) + b1(i,l)*directclr(l))
            directcld(l+1) = exptaucld * &
                (b3(i,l)*directcld(l) + (1.-b1(i,l))*directclr(l))
         enddo
 
!---- 2. LOAD A MATRIX, B MATRIX
         a(:,:) = 0.0
 
         a(1,6)  = -1.0
         a(1,10) = trcld(1) * b4(i,1)
         a(1,11) = trcld(1) * (1.-b2(i,1))
 
         a(2,6)  = -1.0
         a(2,9)  = trclr(1) * (1.-b4(i,1))
         a(2,10) = trclr(1) * b2(i,1)
 
         a(3,6)  = -1.0
         a(3,8)  = rrcld(1) * b4(i,1)
         a(3,9)  = rrcld(1) * (1.-b2(i,1))
 
         a(4,6)  = -1.0
         a(4,7)  = rrclr(1) * (1.-b4(i,1))
         a(4,8)  = rrclr(1) * b2(i,1)
 
         do l = 2,nlm
            a(l*4-3,4)  = rrcld(l) * b3(i,l)
            a(l*4-3,5)  = rrcld(l) * (1.-b1(i,l))
            a(l*4-3,6)  = -1.0
            a(l*4-3,10) = trcld(l) * b4(i,l)
            a(l*4-3,11) = trcld(l) * (1.-b2(i,l))
 
            a(l*4-2,3)  = rrclr(l) * (1.-b3(i,l))
            a(l*4-2,4)  = rrclr(l) * b1(i,l)
            a(l*4-2,6)  = -1.0
            a(l*4-2,9)  = trclr(l) * (1.-b4(i,l))
            a(l*4-2,10) = trclr(l) * b2(i,l)
 
            a(l*4-1,2)  = trcld(l) * b3(i,l)
            a(l*4-1,3)  = trcld(l) * (1.-b1(i,l))
            a(l*4-1,6)  = -1.0
            a(l*4-1,8)  = rrcld(l) * b4(i,l)
            a(l*4-1,9)  = rrcld(l) * (1.-b2(i,l))
 
            a(l*4,1)    = trclr(l) * (1.-b3(i,l))
            a(l*4,2)    = trclr(l) * b1(i,l)
            a(l*4,6)    = -1.0
            a(l*4,7)    = rrclr(l) * (1.-b4(i,l))
            a(l*4,8)    = rrclr(l) * b2(i,l)
         enddo
 
         a(nlm*4+1,4) = asdif(i,ib)
         a(nlm*4+1,6) = -1.0
         a(nlm*4+2,4) = asdif(i,ib)
         a(nlm*4+2,6) = -1.0

         b(:) = 0.0
         do l = 1,nlm
!            b(l*4-3) = -1.*(b3(i,l)*sigucld(l) + (1.-b1(i,l)) *siguclr(l))
!            b(l*4-2) = -1.*((1.-b3(i,l))*sigucld(l) + b1(i,l) *siguclr(l))
!            b(l*4-1) = -1.*(b3(i,l)*sigdcld(l) + (1.-b1(i,l)) *sigdclr(l))
!            b(l*4)   = -1.*((1.-b3(i,l))*sigdcld(l) + b1(i,l) *sigdclr(l))
            b(l*4-3) = -sigucld(l)
            b(l*4-2) = -siguclr(l)
            b(l*4-1) = -sigdcld(l)
            b(l*4)   = -sigdclr(l)
         enddo
         b(nlm*4+1) = -asdir(i,ib)*directcld(nlm+1)*amu0(i)
         b(nlm*4+2) = -asdir(i,ib)*directclr(nlm+1)*amu0(i)
 
!         do l = 1,nlm*4+2
!            print '(15E15.7)',(a(l,k),k=1,11)
!         enddo
 
!         do l = 1,nlm*4+2
!            print *,l,b(l)
!         enddo
!
         call bandec(a,nlm*4+2,5,5,nlm*4+2,11,al,11,indx,d)
         call banbks(a,nlm*4+2,5,5,nlm*4+2,11,al,11,indx,b)
 
!---- 3. SUM CLEAR AND CLOUDY FLUXES
         do l = 1,nlm+1
            fudif(i,l) = slr(i)*(b(l*4-3)+b(l*4-2))
         enddo
         do l = 1,nlm
            fddif(i,l+1) = slr(i)*(b(l*4-1)+b(l*4))
         enddo
 
         fddir(i,1) = amu0(i)*slr(i)
         do l = 1, nlm
            fddir(i,l+1) = amu0(i)*slr(i)* &
                           (directcld(l+1)+directclr(l+1))
         enddo
!         do l = 1, nlm+1
!            print *,ib,l,fddir(l),fddif(l),fudif(l)
!         enddo
!>> END FULL CALCULATION <<

 1000 continue
 
      return
      end subroutine two_rt_sw_bs


! CVS:  $Id: two_rt_sw.F,v 1.7 2006/11/16 21:19:38 norm Exp $
! CVS:  $Name:  $

!NBW This has been modified to apply delta-M descaling
!NBW to the direct and diffuse downwelling fluxes


!-----------------------------------------------------------------------
 
      subroutine two_rt_sw &
                     ( ncol ,   nlm ,      mbs ,    ib, &
                        slr ,  amu0 ,       wc ,  asym, &
                        tau , asdir ,    asdif , fudif, &
                      fddir , fddif ,sel_rules &
                     )
 
      use bugs_kinds



      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_sw replaces two_rt and add written by G. Stephens. two_rt_sw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in):: &
        sel_rules

      integer (kind=int_kind), intent(in):: &
        ncol, &  !Length of sub-domain.
        nlm, &   !Number of layers.
        mbs, &   !Number of SW spectral intervals.
        ib    !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol):: &
        slr, &   !Fraction of daylight                                (-).
        amu0  !Cosine of the solar zenith angle                    (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,mbs):: &
        asdir, & !Spectral direct surface albedo                      (-).
        asdif !Spectral diffuse surface albedo                     (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm):: &
        wc, &    !Single scattering albedo                            (-).
        asym, &  !Asymmetry factor                                    (-).
        tau   !Optical depth                                       (-).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1):: &
        fddir, & !Spectral direct downward flux                (W/m^2).
        fddif, & !Spectral diffuse downward flux               (W/m^2).
        fudif    !Spectral diffuse upward flux                 (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind):: &
        i, &     ! Horizontal index.
        l     ! Vertical index.
 
      real (kind=dbl_kind), dimension(nlm):: &
        rr, &      !
        tr, &      !
        sigu, &    !
        sigd       !

      real (kind=dbl_kind):: &
        exptau_s(nlm),direct_s(nlm+1),exptau_us(nlm),direct_us(nlm+1)

      real (kind=dbl_kind):: &
            aa ,    bb ,   cc , denom, &
        eggtau ,   eps ,   g3 ,    g4, &
         ggtau , kappa ,    r ,  rinf, &
             t ,   oms , taus ,  fact, &
           asy ,  prop, fd_total
      data eps /1.e-02/

      real (kind=dbl_kind), dimension(nlm):: &
        td , vu

      real (kind=dbl_kind), dimension(nlm+1):: &    
        re , vd
 
! SELECTION RULE VARIABLES

      logical (kind=log_kind):: &
        fail

      real (kind=dbl_kind):: &
        tausthresh, &
        wcthresh, &
        tauscat


!     data tausthresh / 0.001 /
!     data wcthresh   / 0.975 /
      data tausthresh / 0.01 /
      data wcthresh   / 0.98 /
 
!----------------------------------------------------------------------
 
      do 1000 i = 1, ncol
 
         if(sel_rules) then

            fail = .false.
            tauscat = 0.0
            do l = nlm, 1, -1
               if (wc(i,l).gt.wcthresh) fail = .true.
               tauscat = tauscat + wc(i,l)*tau(i,l)
            enddo
            
            if (fail.and.tauscat.ge.tausthresh) goto 2000
 
!>> BEGIN SELECTION RULES <<
           print *,'selection rules'
            fddir(i,1) = amu0(i)*slr(i)
            fddif(i,:) = 0.0
            do l = 1, nlm
               fddir(i,l+1) = fddir(i,l) * exp(-1.*tau(i,l)/amu0(i))
            enddo
            fudif(i,nlm+1) = fddir(i,nlm+1) * asdir(i,ib)
 
            do l = nlm, 1, -1
               fudif(i,l) = fudif(i,l+1) * exp(-2.*tau(i,l))
            enddo

            cycle

         endif
!>> END SELECTION RULES <<

!>> BEGIN FULL CALCULATION <<
2000     direct_s(1) = 1.
         direct_us(1) = 1.
         re(1) = 0.
         vd(1) = 0.
 
         do l = 1, nlm
            !delta-M on
            fact = asym(i,l)*asym(i,l)
            asy = asym(i,l)/(1.+asym(i,l))    !Orig w/ delta M
            !This is the delta-M scaling approx
            !delta-M off
            !fact = 0.
            !asy = asym(i,l)
            !End delta-M on/off
            oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
            taus   = (1.-fact*wc(i,l))*tau(i,l)
            
 
            exptau_s(l) = exp(-taus/amu0(i))         !delta-M scaled
            direct_s(l+1) = exptau_s(l)*direct_s(l)  !delta-M scaled
            exptau_us(l) = exp(-tau(i,l)/amu0(i))            !unscaled
            direct_us(l+1) = exptau_us(l)*direct_us(l)       !unscaled
 
!---- local coefficients:  delta-eddington
            t      = 0.25 * (7. - oms*(4.+3.*asy))
            r      = -0.25 * (1. - oms*(4.-3.*asy))
            kappa  = sqrt(t**2-r**2)
            rinf   = r/(kappa+t)
            ggtau  = kappa*taus
            eggtau = exp(-ggtau)
            denom  = (1.-rinf**2*eggtau**2)
            tr(l)  = (1.-rinf**2)*eggtau/denom
            rr(l)  = rinf*(1.-eggtau**2)/denom
 
            if(abs(kappa**2-1./amu0(i)**2) .lt. eps) then
               fact = 1./eps
            else
               fact = 1./(kappa**2-1./amu0(i)**2)
            endif
            cc = oms*slr(i)*fact
            g3 = 0.5-0.75*asy*amu0(i)
            g4 = 1.-g3
            aa = g3*(t-1./amu0(i))+g4*r
            bb = g4*(t+1./amu0(i))+g3*r
            sigu(l) = cc*((aa-rr(l)*bb)-aa*tr(l)*exptau_s(l)) &
                    * direct_s(l)
            sigd(l) = cc*(-bb*tr(l)+(bb-rr(l)*aa)*exptau_s(l)) &
                    * direct_s(l)
         enddo
 
!---- 1. do adding, going from top down:

         do l = 1, nlm
            prop = 1. / (1. - re(l)*rr(l))
            re(l+1) = rr(l) + tr(l)**2*re(l)*prop
            vd(l+1) = sigd(l) + (tr(l)*vd(l) &
                    + tr(l)*re(l)*sigu(l))*prop
            vu(l)   = (rr(l)*vd(l) + sigu(l))*prop
            td(l)   = prop
         enddo
 
!---- 2. calculate diffuse fluxes:

         fddif(i,1) = 0.
         fudif(i,nlm+1) = (asdif(i,ib)*vd(nlm+1) &
                        + asdir(i,ib)*slr(i)*amu0(i)*direct_s(nlm+1)) &
                        / (1.-asdif(i,ib)*re(nlm+1))
 
         do l = nlm+1, 2, -1
            fddif(i,l)   = re(l)*fudif(i,l) + vd(l)
            fudif(i,l-1) = tr(l-1)*fudif(i,l)*td(l-1) + vu(l-1)
         enddo
 
!---- 3. Compute direct flux and descale the direct and diffuse fluxes
         fddir(i,1) = amu0(i)*slr(i)
         do l = 1, nlm
            fd_total = amu0(i)*slr(i)*direct_s(l+1) + fddif(i,l+1)
            fddir(i,l+1) = amu0(i)*slr(i)*direct_us(l+1)
            fddif(i,l+1) = fd_total - fddir(i,l+1)
         enddo
         
!>> END FULL CALCULATION <<

 1000 continue
      return
      end subroutine two_rt_sw
 
!------------------------------------------------------------------------
