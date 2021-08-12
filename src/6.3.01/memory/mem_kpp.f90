!##############################################################################
Module mem_kpp
  
use kpp_parameters, only:kppNVEL,kppNSCLR
use grid_dims, only:maxkppz

implicit none

  Type kpp_3d_type

    !The necessary z,x,y dimension categories are given here, but for RAMS
    !simplicity, we will allocate over NKPPZ for all 3D vars. This prevents
    !having to have too many IDIM_TYPE categories, which would be a pain.
    !Note: RAMSIN-> NKPPZ = number of KPP layers
    ! IN KPP code, kppNZ = NKPPZ-1 and kppNZP1 = NKPPZ
    ! So we cover the 0th array slot for some variables that loop to kppNZ

    !nxp,nyp
    real, allocatable, dimension(:,:) ::        &
               old & !index to ID past value of Us,Xs
              ,new & !index to ID new  value of Us,Xs
           ,jerlov & !ocean optical clarity category for radiation effects
          ,ocdepth & !ocean depth (meters)
             ,hmix & !mixed-layer / boundary layer depth (meters)
                ,f & !coriolis parameter
          ,bottomt & !ocean bottom temperature (celcius in KPP)
             ,SST0 & !temp to which to relax mixed-layer temp (celcius)
             ,Sref & !reference salinity
       ,freez_flag & !fraction of levels prevented from freezing, set to -1.8C
       ,reset_flag & !flag to indicate isothermal column, reset T/S to climo
          ,flx_ust & !sflux(1) zonal surface wind stress (N/m2)
          ,flx_vst & !sflux(2) meridional surface wind stress (N/m2)
          ,flx_nsw & !sflux(3) net sfc shortwave radiation (swdn-swup)(W/m^2)
          ,flx_nlw & !sflux(4) non-shortwave rad (lwdn-lwup-sens-latent) (W/m^2)
          ,flx_ice & !sflux(5) melting of seaice(zero, not used)
          ,flx_pcp   !sflux(6) net freshwater (precip - evaporation)(mm/sec)

    !0:nz,nxp,nyp
    !Variables on grid cell interaces with 0 as the ocean surface
    real, allocatable, dimension(:,:,:) ::     &
          swdk_opt & !solar shortwave flux fraction on DM depths
               ,wB & !w'B' total kinematic buoyancy flux (m**2/s**3)
            ,wXNTt & !w'T'(NT) (non-turbulent temperature flux, degC m/s)
               ,wU & !w'u'(turbulent zonal velocity flux, m2/s2) 
               ,wV & !and w'v' (turbulent meridional velocity flux, m2/s2)
              ,wXt & !w'T' (turbulent temperature flux, degC m/s) 
              ,wXs   !and w'S' (turbulent salinity flux, o/oo 1/s)

    !nzp1,nxp,nyp
    !Variables on grid cell center levels with NZP1 as the model bottom.
    real, allocatable, dimension(:,:,:) ::     &
            swfrac & !solar shortwave flux fraction on ZM depths
       ,tinc_fcorr & !temperature increment from flux corrections with depth,(K)
       ,sinc_fcorr & !salinity increment from flux corrections with depth,(o/oo)
         ,sal_clim & !3D salinity climatology (can be time updated)
        ,ocnT_clim & !3D temperature climatology (can be time updated)
             ,buoy & !buoyancy (m/s2)
              ,rho & !density (kg/m3)
               ,cp & !specific heat capacity (J/kg/K)
                ,U & !latest value of ocean currents velocity U (m/s)
                ,V & !latest value of ocean currents velocity V (m/s)
           ,U_init & !initial value of ocean currents velocity U
           ,V_init & !initial value of ocean currents velocity V
              ,Us0 & !contain old/new U ocean currents (m/s)
              ,Us1 & !contain old/new U ocean currents (m/s)
              ,Vs0 & !contain old/new V ocean currents (m/s)
              ,Vs1 & !contain old/new V ocean currents (m/s)
              ,X_t & !latest ocean column temperature (Celcius)
              ,X_s & !latest ocean column salinity (o/oo) (+Sref)
            ,Xs_t0 & !contain old/new Temperature
            ,Xs_t1 & !contain old/new Temperature
            ,Xs_s0 & !contain old/new Salinity
            ,Xs_s1   !contain old/new Salinity

  End Type kpp_3d_type

  !******************************************************************
  Type kpp_1d_type

    real    ::                           &
               old                       &
              ,new                       &
           ,jerlov                       &
          ,ocdepth                       &
             ,hmix                       &
             ,kmix                       &
                ,f                       &
          ,bottomt                       &
             ,SST0                       &
             ,Tref                       &
             ,uref                       &
             ,vref                       &
             ,Sref                       &
       ,freez_flag                       &
       ,reset_flag                       &
       ,dampu_flag                       &
       ,dampv_flag                       &
            ,dbloc(maxkppz)              &
             ,ghat(maxkppz)              &
              ,Rig(maxkppz)              &
             ,Shsq(maxkppz)              &
           ,swfrac(maxkppz)              &
       ,tinc_fcorr(maxkppz)              &
       ,sinc_fcorr(maxkppz)              &
         ,sal_clim(maxkppz)              &
        ,ocnT_clim(maxkppz)              &
             ,buoy(maxkppz)              &
              ,rho(maxkppz)              &
               ,cp(maxkppz)              &
           ,talpha(maxkppz)              & !coefficient for buoyancy calculation
            ,sbeta(maxkppz)              & !coefficient for buoyancy calculation
           ,U_init(maxkppz,kppNVEL)      &
                ,U(maxkppz,kppNVEL)      &
               ,Us(maxkppz,kppNVEL,0:1)  &
                ,X(maxkppz,kppNSCLR)     &
               ,Xs(maxkppz,kppNSCLR,0:1) &
         ,swdk_opt(0:maxkppz)            &
             ,difm(0:maxkppz)            &
             ,difs(0:maxkppz)            &
             ,dift(0:maxkppz)            &
               ,wB(0:maxkppz)            &
             ,wXNT(0:maxkppz)            &
               ,wU(0:maxkppz,kppNVEL)    &
               ,wX(0:maxkppz,kppNSCLR)   &
            ,sflux(6)

   logical :: comp_flag !flag to decide if more convergence iterations needed

  End Type kpp_1D_type

  !******************************************************************
  Type kpp_const_type

   real ::   dto &
      ,relax_sal &
     ,relax_ocnT &
      ,relax_sst &
      ,dt_uvdamp &
           ,DMAX &
         ,dscale &
    ,hmixtolfrac &
           ,sice &
     ,iso_thresh

   integer :: iso_bot

   logical :: LKPP               &
             ,LRI                &
             ,LDD                &
             ,L_JERLOV           &
             ,L_STRETCHGRID      &
             ,L_NO_FREEZE        &
             ,L_NO_ISOTHERM      &
             ,L_BOTTOM_TEMP      &
             ,L_DAMP_CURR

   real :: &
      zm(maxkppz) &  !ocean layer grid center depth (m) (1,kppnzp1)
     ,hm(maxkppz) &  !ocean layer thickness (m)         (1,kppnzp1)
     ,dm(0:maxkppz)  !ocean layer interface depths (m) (0,kppnz)

   real, allocatable :: wmt(:,:),wst(:,:),tri(:,:)

  End Type kpp_const_type

  !******************************************************************
  type (kpp_3d_type), allocatable :: kpp_3d_fields(:), kpp_3d_fieldsm(:)
  type (kpp_1d_type) :: kpp_1d_fields
  type (kpp_const_type) :: kpp_const_fields

Contains

!##############################################################################
Subroutine alloc_kpp (kpp_3d_fields,nx,ny,nkppz,ng)

use kpp_parameters, only:IKPP

implicit none

type (kpp_3d_type) :: kpp_3d_fields
integer, intent(in) :: nx,ny,nkppz,ng

!The necessary z,x,y dimension categories are given in "module mem_kpp",
!but for RAMS simplicity, we will allocate over 0:NZP1 for all 3D vars. 
!This prevents having to have too many IDIM_TYPE categories.
!KPP still runs only over necessary levels/layers.

if(IKPP > 0) then
  allocate (kpp_3d_fields%Sref        (nx,ny))
  allocate (kpp_3d_fields%jerlov      (nx,ny))
  allocate (kpp_3d_fields%ocdepth     (nx,ny))
  allocate (kpp_3d_fields%f           (nx,ny))
  allocate (kpp_3d_fields%bottomt     (nx,ny))
  allocate (kpp_3d_fields%SST0        (nx,ny))
  allocate (kpp_3d_fields%old         (nx,ny))
  allocate (kpp_3d_fields%new         (nx,ny))
  allocate (kpp_3d_fields%flx_ust     (nx,ny))
  allocate (kpp_3d_fields%flx_vst     (nx,ny))
  allocate (kpp_3d_fields%flx_nsw     (nx,ny))
  allocate (kpp_3d_fields%flx_nlw     (nx,ny))
  allocate (kpp_3d_fields%flx_ice     (nx,ny))
  allocate (kpp_3d_fields%flx_pcp     (nx,ny))
  allocate (kpp_3d_fields%hmix        (nx,ny)) !Diagnostic
  allocate (kpp_3d_fields%freez_flag  (nx,ny)) !Diagnostic
  allocate (kpp_3d_fields%reset_flag  (nx,ny)) !Diagnostic
  allocate (kpp_3d_fields%swdk_opt    (nkppz,nx,ny))
  allocate (kpp_3d_fields%U_init      (nkppz,nx,ny))
  allocate (kpp_3d_fields%V_init      (nkppz,nx,ny))
  allocate (kpp_3d_fields%swfrac      (nkppz,nx,ny))
  allocate (kpp_3d_fields%sal_clim    (nkppz,nx,ny))
  allocate (kpp_3d_fields%ocnT_clim   (nkppz,nx,ny))
  allocate (kpp_3d_fields%U           (nkppz,nx,ny))
  allocate (kpp_3d_fields%V           (nkppz,nx,ny))
  allocate (kpp_3d_fields%Us0         (nkppz,nx,ny))
  allocate (kpp_3d_fields%Us1         (nkppz,nx,ny))
  allocate (kpp_3d_fields%Vs0         (nkppz,nx,ny))
  allocate (kpp_3d_fields%Vs1         (nkppz,nx,ny))
  allocate (kpp_3d_fields%X_t         (nkppz,nx,ny))
  allocate (kpp_3d_fields%X_s         (nkppz,nx,ny))
  allocate (kpp_3d_fields%Xs_t0       (nkppz,nx,ny))
  allocate (kpp_3d_fields%Xs_t1       (nkppz,nx,ny))
  allocate (kpp_3d_fields%Xs_s0       (nkppz,nx,ny))
  allocate (kpp_3d_fields%Xs_s1       (nkppz,nx,ny))
  allocate (kpp_3d_fields%buoy        (nkppz,nx,ny)) !Diagnostic
endif
if(IKPP==2) then
  allocate (kpp_3d_fields%wB          (nkppz,nx,ny))
  allocate (kpp_3d_fields%wXNTt       (nkppz,nx,ny))
  allocate (kpp_3d_fields%wU          (nkppz,nx,ny))
  allocate (kpp_3d_fields%wV          (nkppz,nx,ny))
  allocate (kpp_3d_fields%wXt         (nkppz,nx,ny))
  allocate (kpp_3d_fields%wXs         (nkppz,nx,ny))
  allocate (kpp_3d_fields%tinc_fcorr  (nkppz,nx,ny))
  allocate (kpp_3d_fields%sinc_fcorr  (nkppz,nx,ny))
  allocate (kpp_3d_fields%rho         (nkppz,nx,ny))
  allocate (kpp_3d_fields%cp          (nkppz,nx,ny))
endif

return
END SUBROUTINE alloc_kpp

!##############################################################################
Subroutine dealloc_kpp (kpp_3d_fields)

implicit none

type (kpp_3d_type) :: kpp_3d_fields

if(allocated(kpp_3d_fields%old))        deallocate (kpp_3d_fields%old)
if(allocated(kpp_3d_fields%new))        deallocate (kpp_3d_fields%new)
if(allocated(kpp_3d_fields%jerlov))     deallocate (kpp_3d_fields%jerlov)
if(allocated(kpp_3d_fields%ocdepth))    deallocate (kpp_3d_fields%ocdepth)
if(allocated(kpp_3d_fields%hmix))       deallocate (kpp_3d_fields%hmix)
if(allocated(kpp_3d_fields%f))          deallocate (kpp_3d_fields%f)
if(allocated(kpp_3d_fields%bottomt))    deallocate (kpp_3d_fields%bottomt)
if(allocated(kpp_3d_fields%SST0))       deallocate (kpp_3d_fields%SST0)
if(allocated(kpp_3d_fields%Sref))       deallocate (kpp_3d_fields%Sref)
if(allocated(kpp_3d_fields%freez_flag)) deallocate (kpp_3d_fields%freez_flag)
if(allocated(kpp_3d_fields%reset_flag)) deallocate (kpp_3d_fields%reset_flag)
if(allocated(kpp_3d_fields%flx_ust))    deallocate (kpp_3d_fields%flx_ust)
if(allocated(kpp_3d_fields%flx_vst))    deallocate (kpp_3d_fields%flx_vst)
if(allocated(kpp_3d_fields%flx_nsw))    deallocate (kpp_3d_fields%flx_nsw)
if(allocated(kpp_3d_fields%flx_nlw))    deallocate (kpp_3d_fields%flx_nlw)
if(allocated(kpp_3d_fields%flx_ice))    deallocate (kpp_3d_fields%flx_ice)
if(allocated(kpp_3d_fields%flx_pcp))    deallocate (kpp_3d_fields%flx_pcp)
if(allocated(kpp_3d_fields%swdk_opt))   deallocate (kpp_3d_fields%swdk_opt)
if(allocated(kpp_3d_fields%wB))         deallocate (kpp_3d_fields%wB)
if(allocated(kpp_3d_fields%wXNTt))      deallocate (kpp_3d_fields%wXNTt)
if(allocated(kpp_3d_fields%wU))         deallocate (kpp_3d_fields%wU)
if(allocated(kpp_3d_fields%wV))         deallocate (kpp_3d_fields%wV)
if(allocated(kpp_3d_fields%wXt))        deallocate (kpp_3d_fields%wXt)
if(allocated(kpp_3d_fields%wXs))        deallocate (kpp_3d_fields%wXs)
if(allocated(kpp_3d_fields%swfrac))     deallocate (kpp_3d_fields%swfrac)
if(allocated(kpp_3d_fields%tinc_fcorr)) deallocate (kpp_3d_fields%tinc_fcorr)
if(allocated(kpp_3d_fields%sinc_fcorr)) deallocate (kpp_3d_fields%sinc_fcorr)
if(allocated(kpp_3d_fields%sal_clim))   deallocate (kpp_3d_fields%sal_clim)
if(allocated(kpp_3d_fields%ocnT_clim))  deallocate (kpp_3d_fields%ocnT_clim)
if(allocated(kpp_3d_fields%buoy))       deallocate (kpp_3d_fields%buoy)
if(allocated(kpp_3d_fields%rho))        deallocate (kpp_3d_fields%rho)
if(allocated(kpp_3d_fields%cp))         deallocate (kpp_3d_fields%cp)
if(allocated(kpp_3d_fields%U))          deallocate (kpp_3d_fields%U)
if(allocated(kpp_3d_fields%V))          deallocate (kpp_3d_fields%V)
if(allocated(kpp_3d_fields%U_init))     deallocate (kpp_3d_fields%U_init)
if(allocated(kpp_3d_fields%V_init))     deallocate (kpp_3d_fields%V_init)
if(allocated(kpp_3d_fields%Us0))        deallocate (kpp_3d_fields%Us0)
if(allocated(kpp_3d_fields%Us1))        deallocate (kpp_3d_fields%Us1)
if(allocated(kpp_3d_fields%Vs0))        deallocate (kpp_3d_fields%Vs0)
if(allocated(kpp_3d_fields%Vs1))        deallocate (kpp_3d_fields%Vs1)
if(allocated(kpp_3d_fields%X_t))        deallocate (kpp_3d_fields%X_t)
if(allocated(kpp_3d_fields%X_s))        deallocate (kpp_3d_fields%X_s)
if(allocated(kpp_3d_fields%Xs_t0))      deallocate (kpp_3d_fields%Xs_t0)
if(allocated(kpp_3d_fields%Xs_t1))      deallocate (kpp_3d_fields%Xs_t1)
if(allocated(kpp_3d_fields%Xs_s0))      deallocate (kpp_3d_fields%Xs_s0)
if(allocated(kpp_3d_fields%Xs_s1))      deallocate (kpp_3d_fields%Xs_s1)

return
END SUBROUTINE dealloc_kpp

!##############################################################################
Subroutine filltab_kpp (kpp_3d_fields,kpp_3d_fieldsm,imean,nx,ny,nkppz,ng)

use var_tables

implicit none

type (kpp_3d_type) :: kpp_3d_fields,kpp_3d_fieldsm
integer, intent(in) :: imean,nx,ny,ng,nkppz
integer :: npts

!Fill arrays into variable tables

!The necessary z,x,y dimension categories are given in "module mem_kpp",
!but for RAMS simplicity, we will allocate over 0:NZP1 for all 3D vars. 
!This prevents having to have too many IDIM_TYPE categories.
!KPP still runs only over necessary levels/layers.

 npts=nx*ny !2D fields horizontal fields

 if (allocated(kpp_3d_fields%old)) &
  CALL vtables2 (kpp_3d_fields%old,kpp_3d_fieldsm%old &
               ,ng, npts, imean,  &
               'KPP_OLD :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%new)) &
  CALL vtables2 (kpp_3d_fields%new,kpp_3d_fieldsm%new &
               ,ng, npts, imean,  &
               'KPP_NEW :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%jerlov)) &
  CALL vtables2 (kpp_3d_fields%jerlov,kpp_3d_fieldsm%jerlov &
               ,ng, npts, imean,  &
               'KPP_JERLOV :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%ocdepth)) &
  CALL vtables2 (kpp_3d_fields%ocdepth,kpp_3d_fieldsm%ocdepth &
               ,ng, npts, imean,  &
               'KPP_OCDEPTH :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%hmix)) &
  CALL vtables2 (kpp_3d_fields%hmix,kpp_3d_fieldsm%hmix &
               ,ng, npts, imean,  &
               'KPP_HMIX :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%f)) &
  CALL vtables2 (kpp_3d_fields%f,kpp_3d_fieldsm%f &
               ,ng, npts, imean,  &
               'KPP_F :2:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%bottomt)) &
  CALL vtables2 (kpp_3d_fields%bottomt,kpp_3d_fieldsm%bottomt &
               ,ng, npts, imean,  &
               'KPP_BOTTOMT :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%SST0)) &
  CALL vtables2 (kpp_3d_fields%SST0,kpp_3d_fieldsm%SST0 &
               ,ng, npts, imean,  &
               'KPP_SST0 :2:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Sref)) &
  CALL vtables2 (kpp_3d_fields%Sref,kpp_3d_fieldsm%Sref &
               ,ng, npts, imean,  &
               'KPP_SREF :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%freez_flag)) &
  CALL vtables2 (kpp_3d_fields%freez_flag,kpp_3d_fieldsm%freez_flag &
               ,ng, npts, imean,  &
               'KPP_FREEZ_FLAG :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%reset_flag)) &
  CALL vtables2 (kpp_3d_fields%reset_flag,kpp_3d_fieldsm%reset_flag &
               ,ng, npts, imean,  &
               'KPP_RESET_FLAG :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_ust)) &
  CALL vtables2 (kpp_3d_fields%flx_ust,kpp_3d_fieldsm%flx_ust &
               ,ng, npts, imean,  &
               'KPP_FLX_UST :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_vst)) &
  CALL vtables2 (kpp_3d_fields%flx_vst,kpp_3d_fieldsm%flx_vst &
               ,ng, npts, imean,  &
               'KPP_FLX_VST :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_nsw)) &
  CALL vtables2 (kpp_3d_fields%flx_nsw,kpp_3d_fieldsm%flx_nsw &
               ,ng, npts, imean,  &
               'KPP_FLX_NSW :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_nlw)) &
  CALL vtables2 (kpp_3d_fields%flx_nlw,kpp_3d_fieldsm%flx_nlw &
               ,ng, npts, imean,  &
               'KPP_FLX_NLW :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_ice)) &
  CALL vtables2 (kpp_3d_fields%flx_ice,kpp_3d_fieldsm%flx_ice &
               ,ng, npts, imean,  &
               'KPP_FLX_ICE :2:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%flx_pcp)) &
  CALL vtables2 (kpp_3d_fields%flx_pcp,kpp_3d_fieldsm%flx_pcp &
               ,ng, npts, imean,  &
               'KPP_FLX_PCP :2:anal:mpti:recycle_sfc')

 !******************************************************************************
 npts=nkppz*nx*ny !fields that are allocated 0:nkppz+1,nx,ny

 if (allocated(kpp_3d_fields%swdk_opt)) &
  CALL vtables2 (kpp_3d_fields%swdk_opt,kpp_3d_fieldsm%swdk_opt &
               ,ng, npts, imean,  &
               'KPP_SWDK_OPT :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wB)) &
  CALL vtables2 (kpp_3d_fields%wB,kpp_3d_fieldsm%wB &
               ,ng, npts, imean,  &
               'KPP_wB :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wU)) &
  CALL vtables2 (kpp_3d_fields%wU,kpp_3d_fieldsm%wU &
               ,ng, npts, imean,  &
               'KPP_wU :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wV)) &
  CALL vtables2 (kpp_3d_fields%wV,kpp_3d_fieldsm%wV &
               ,ng, npts, imean,  &
               'KPP_wV :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wXt)) &
  CALL vtables2 (kpp_3d_fields%wXt,kpp_3d_fieldsm%wXt &
               ,ng, npts, imean,  &
               'KPP_wXt :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wXs)) &
  CALL vtables2 (kpp_3d_fields%wXs,kpp_3d_fieldsm%wXs &
               ,ng, npts, imean,  &
               'KPP_wXs :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%wXNTt)) &
  CALL vtables2 (kpp_3d_fields%wXNTt,kpp_3d_fieldsm%wXNTt &
               ,ng, npts, imean,  &
               'KPP_wXNTt :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%swfrac)) &
  CALL vtables2 (kpp_3d_fields%swfrac,kpp_3d_fieldsm%swfrac &
               ,ng, npts, imean,  &
               'KPP_SWFRAC :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%tinc_fcorr)) &
  CALL vtables2 (kpp_3d_fields%tinc_fcorr,kpp_3d_fieldsm%tinc_fcorr &
               ,ng, npts, imean,  &
               'KPP_TINC_FCORR :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%sinc_fcorr)) &
  CALL vtables2 (kpp_3d_fields%sinc_fcorr,kpp_3d_fieldsm%sinc_fcorr &
               ,ng, npts, imean,  &
               'KPP_SINC_FCORR :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%sal_clim)) &
  CALL vtables2 (kpp_3d_fields%sal_clim,kpp_3d_fieldsm%sal_clim &
               ,ng, npts, imean,  &
               'KPP_SAL_CLIM :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%ocnT_clim)) &
  CALL vtables2 (kpp_3d_fields%ocnT_clim,kpp_3d_fieldsm%ocnT_clim &
               ,ng, npts, imean,  &
               'KPP_OCNT_CLIM :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%buoy)) &
  CALL vtables2 (kpp_3d_fields%buoy,kpp_3d_fieldsm%buoy &
               ,ng, npts, imean,  &
               'KPP_BUOY :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%rho)) &
  CALL vtables2 (kpp_3d_fields%rho,kpp_3d_fieldsm%rho &
               ,ng, npts, imean,  &
               'KPP_RHO :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%cp)) &
  CALL vtables2 (kpp_3d_fields%cp,kpp_3d_fieldsm%cp &
               ,ng, npts, imean,  &
               'KPP_CP :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%U)) &
  CALL vtables2 (kpp_3d_fields%U,kpp_3d_fieldsm%U &
               ,ng, npts, imean,  &
               'KPP_U :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%V)) &
  CALL vtables2 (kpp_3d_fields%V,kpp_3d_fieldsm%V &
               ,ng, npts, imean,  &
               'KPP_V :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%U_init)) &
  CALL vtables2 (kpp_3d_fields%U_init,kpp_3d_fieldsm%U_init &
               ,ng, npts, imean,  &
               'KPP_U_init :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%V_init)) &
  CALL vtables2 (kpp_3d_fields%V_init,kpp_3d_fieldsm%V_init &
               ,ng, npts, imean,  &
               'KPP_V_init :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Us0)) &
  CALL vtables2 (kpp_3d_fields%Us0,kpp_3d_fieldsm%Us0 &
               ,ng, npts, imean,  &
               'KPP_US0 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Us1)) &
  CALL vtables2 (kpp_3d_fields%Us1,kpp_3d_fieldsm%Us1 &
               ,ng, npts, imean,  &
               'KPP_US1 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Vs0)) &
  CALL vtables2 (kpp_3d_fields%Vs0,kpp_3d_fieldsm%Vs0 &
               ,ng, npts, imean,  &
               'KPP_VS0 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Vs1)) &
  CALL vtables2 (kpp_3d_fields%Vs1,kpp_3d_fieldsm%Vs1 &
               ,ng, npts, imean,  &
               'KPP_VS1 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%X_t)) &
  CALL vtables2 (kpp_3d_fields%X_t,kpp_3d_fieldsm%X_t &
               ,ng, npts, imean,  &
               'KPP_X_T :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%X_s)) &
  CALL vtables2 (kpp_3d_fields%X_s,kpp_3d_fieldsm%X_s &
               ,ng, npts, imean,  &
               'KPP_X_S :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Xs_t0)) &
  CALL vtables2 (kpp_3d_fields%Xs_t0,kpp_3d_fieldsm%Xs_t0 &
               ,ng, npts, imean,  &
               'KPP_XS_T0 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Xs_t1)) &
  CALL vtables2 (kpp_3d_fields%Xs_t1,kpp_3d_fieldsm%Xs_t1 &
               ,ng, npts, imean,  &
               'KPP_XS_T1 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Xs_s0)) &
  CALL vtables2 (kpp_3d_fields%Xs_s0,kpp_3d_fieldsm%Xs_s0 &
               ,ng, npts, imean,  &
               'KPP_XS_S0 :10:anal:mpti:recycle_sfc')
 if (allocated(kpp_3d_fields%Xs_s1)) &
  CALL vtables2 (kpp_3d_fields%Xs_s1,kpp_3d_fieldsm%Xs_s1 &
               ,ng, npts, imean,  &
               'KPP_XS_S1 :10:anal:mpti:recycle_sfc')

return
END SUBROUTINE filltab_kpp

END MODULE mem_kpp

