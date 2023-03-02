!##############################################################################
Subroutine non_scalar_bc (mzp,mxp,myp,mzg,mzs,ibcon)

use mem_grid
use mem_basic
use mem_radiate
use mem_leaf
use mem_sib
use mem_turb
use mem_micro
use mem_kpp
use kpp_parameters, only:ikpp,nkppz
use micphys, only:level

implicit none

integer :: mzp,mxp,myp,mzg,mzs,ibcon

if (ilwrtyp+iswrtyp > 0) CALL rad_bcond (mzp,mxp,myp,radiate_g(ngrid),ibcon)

if (ilwrtyp+iswrtyp > 0) CALL sfcrad_bcond (mxp,myp,radiate_g(ngrid),ibcon)

CALL leaf_bcond (mxp,myp,mzg,mzs,npatch,leaf_g(ngrid),ibcon)

CALL turb_bcond (mxp,myp,turb_g(ngrid),ibcon)

if (level==3) CALL micro_bcond (mzp,mxp,myp,micro_g(ngrid),ibcon)

if (isfcl==2) CALL sib_bcond (mxp,myp,npatch,sib_g(ngrid),ibcon)

if (ikpp > 0) CALL kpp_bcond (nkppz,mxp,myp,kpp_3d_fields(ngrid),ibcon)

return
END SUBROUTINE non_scalar_bc

!##############################################################################
Subroutine kpp_bcond (mkpp,m2,m3,kpp_fields,ibcon)

use mem_grid
use mem_kpp
use kpp_parameters, only:ikpp

implicit none

type (kpp_3d_type) :: kpp_fields

integer :: mkpp,m2,m3,ibcon,i,j,k

!Lateral boundary set
! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do j = 1,m3
  do k = 1,mkpp
   kpp_fields%swdk_opt(k,1,j)  = kpp_fields%swdk_opt(k,2,j)
   kpp_fields%U_init(k,1,j)    = kpp_fields%U_init(k,2,j)
   kpp_fields%V_init(k,1,j)    = kpp_fields%V_init(k,2,j)
   kpp_fields%swfrac(k,1,j)    = kpp_fields%swfrac(k,2,j)
   kpp_fields%sal_clim(k,1,j)  = kpp_fields%sal_clim(k,2,j)
   kpp_fields%ocnT_clim(k,1,j) = kpp_fields%ocnT_clim(k,2,j)
   kpp_fields%U(k,1,j)         = kpp_fields%U(k,2,j)
   kpp_fields%V(k,1,j)         = kpp_fields%V(k,2,j)
   kpp_fields%Us0(k,1,j)       = kpp_fields%Us0(k,2,j)
   kpp_fields%Us1(k,1,j)       = kpp_fields%Us1(k,2,j)
   kpp_fields%Vs0(k,1,j)       = kpp_fields%Vs0(k,2,j)
   kpp_fields%Vs1(k,1,j)       = kpp_fields%Vs1(k,2,j)
   kpp_fields%X_t(k,1,j)       = kpp_fields%X_t(k,2,j)
   kpp_fields%X_s(k,1,j)       = kpp_fields%X_s(k,2,j)
   kpp_fields%Xs_t0(k,1,j)     = kpp_fields%Xs_t0(k,2,j)
   kpp_fields%Xs_t1(k,1,j)     = kpp_fields%Xs_t1(k,2,j)
   kpp_fields%Xs_s0(k,1,j)     = kpp_fields%Xs_s0(k,2,j)
   kpp_fields%Xs_s1(k,1,j)     = kpp_fields%Xs_s1(k,2,j)
   kpp_fields%buoy(k,1,j)      = kpp_fields%buoy(k,2,j)
   if(IKPP==2) then
    kpp_fields%wB(k,1,j)         = kpp_fields%wB(k,2,j)
    kpp_fields%wXNTt(k,1,j)      = kpp_fields%wXNTt(k,2,j)
    kpp_fields%wU(k,1,j)         = kpp_fields%wU(k,2,j)
    kpp_fields%wV(k,1,j)         = kpp_fields%wV(k,2,j)
    kpp_fields%wXt(k,1,j)        = kpp_fields%wXt(k,2,j)
    kpp_fields%wXs(k,1,j)        = kpp_fields%wXs(k,2,j)
    kpp_fields%tinc_fcorr(k,1,j) = kpp_fields%tinc_fcorr(k,2,j)
    kpp_fields%sinc_fcorr(k,1,j) = kpp_fields%sinc_fcorr(k,2,j)
    kpp_fields%rho(k,1,j)        = kpp_fields%rho(k,2,j)
    kpp_fields%cp(k,1,j)         = kpp_fields%cp(k,2,j)
   endif
  enddo
  kpp_fields%Sref(1,j)       = kpp_fields%Sref(2,j)
  kpp_fields%jerlov(1,j)     = kpp_fields%jerlov(2,j)
  kpp_fields%ocdepth(1,j)    = kpp_fields%ocdepth(2,j)
  kpp_fields%f(1,j)          = kpp_fields%f(2,j)
  kpp_fields%bottomt(1,j)    = kpp_fields%bottomt(2,j)
  kpp_fields%SST0(1,j)       = kpp_fields%SST0(2,j)
  kpp_fields%old(1,j)        = kpp_fields%old(2,j)
  kpp_fields%new(1,j)        = kpp_fields%new(2,j)
  kpp_fields%flx_ust(1,j)    = kpp_fields%flx_ust(2,j)
  kpp_fields%flx_vst(1,j)    = kpp_fields%flx_vst(2,j)
  kpp_fields%flx_nsw(1,j)    = kpp_fields%flx_nsw(2,j)
  kpp_fields%flx_nlw(1,j)    = kpp_fields%flx_nlw(2,j)
  kpp_fields%flx_ice(1,j)    = kpp_fields%flx_ice(2,j)
  kpp_fields%flx_pcp(1,j)    = kpp_fields%flx_pcp(2,j)
  kpp_fields%hmix(1,j)       = kpp_fields%hmix(2,j)
  kpp_fields%freez_flag(1,j) = kpp_fields%freez_flag(2,j)
  kpp_fields%reset_flag(1,j) = kpp_fields%reset_flag(2,j)
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do j = 1,m3
  do k = 1,mkpp
   kpp_fields%swdk_opt(k,m2,j)  = kpp_fields%swdk_opt(k,m2-1,j)
   kpp_fields%U_init(k,m2,j)    = kpp_fields%U_init(k,m2-1,j)
   kpp_fields%V_init(k,m2,j)    = kpp_fields%V_init(k,m2-1,j)
   kpp_fields%swfrac(k,m2,j)    = kpp_fields%swfrac(k,m2-1,j)
   kpp_fields%sal_clim(k,m2,j)  = kpp_fields%sal_clim(k,m2-1,j)
   kpp_fields%ocnT_clim(k,m2,j) = kpp_fields%ocnT_clim(k,m2-1,j)
   kpp_fields%U(k,m2,j)         = kpp_fields%U(k,m2-1,j)
   kpp_fields%V(k,m2,j)         = kpp_fields%V(k,m2-1,j)
   kpp_fields%Us0(k,m2,j)       = kpp_fields%Us0(k,m2-1,j)
   kpp_fields%Us1(k,m2,j)       = kpp_fields%Us1(k,m2-1,j)
   kpp_fields%Vs0(k,m2,j)       = kpp_fields%Vs0(k,m2-1,j)
   kpp_fields%Vs1(k,m2,j)       = kpp_fields%Vs1(k,m2-1,j)
   kpp_fields%X_t(k,m2,j)       = kpp_fields%X_t(k,m2-1,j)
   kpp_fields%X_s(k,m2,j)       = kpp_fields%X_s(k,m2-1,j)
   kpp_fields%Xs_t0(k,m2,j)     = kpp_fields%Xs_t0(k,m2-1,j)
   kpp_fields%Xs_t1(k,m2,j)     = kpp_fields%Xs_t1(k,m2-1,j)
   kpp_fields%Xs_s0(k,m2,j)     = kpp_fields%Xs_s0(k,m2-1,j)
   kpp_fields%Xs_s1(k,m2,j)     = kpp_fields%Xs_s1(k,m2-1,j)
   kpp_fields%buoy(k,m2,j)      = kpp_fields%buoy(k,m2-1,j)
   if(IKPP==2) then
    kpp_fields%wB(k,m2,j)         = kpp_fields%wB(k,m2-1,j)
    kpp_fields%wXNTt(k,m2,j)      = kpp_fields%wXNTt(k,m2-1,j)
    kpp_fields%wU(k,m2,j)         = kpp_fields%wU(k,m2-1,j)
    kpp_fields%wV(k,m2,j)         = kpp_fields%wV(k,m2-1,j)
    kpp_fields%wXt(k,m2,j)        = kpp_fields%wXt(k,m2-1,j)
    kpp_fields%wXs(k,m2,j)        = kpp_fields%wXs(k,m2-1,j)
    kpp_fields%tinc_fcorr(k,m2,j) = kpp_fields%tinc_fcorr(k,m2-1,j)
    kpp_fields%sinc_fcorr(k,m2,j) = kpp_fields%sinc_fcorr(k,m2-1,j)
    kpp_fields%rho(k,m2,j)        = kpp_fields%rho(k,m2-1,j)
    kpp_fields%cp(k,m2,j)         = kpp_fields%cp(k,m2-1,j)
   endif
  enddo
  kpp_fields%Sref(m2,j)       = kpp_fields%Sref(m2-1,j)
  kpp_fields%jerlov(m2,j)     = kpp_fields%jerlov(m2-1,j)
  kpp_fields%ocdepth(m2,j)    = kpp_fields%ocdepth(m2-1,j)
  kpp_fields%f(m2,j)          = kpp_fields%f(m2-1,j)
  kpp_fields%bottomt(m2,j)    = kpp_fields%bottomt(m2-1,j)
  kpp_fields%SST0(m2,j)       = kpp_fields%SST0(m2-1,j)
  kpp_fields%old(m2,j)        = kpp_fields%old(m2-1,j)
  kpp_fields%new(m2,j)        = kpp_fields%new(m2-1,j)
  kpp_fields%flx_ust(m2,j)    = kpp_fields%flx_ust(m2-1,j)
  kpp_fields%flx_vst(m2,j)    = kpp_fields%flx_vst(m2-1,j)
  kpp_fields%flx_nsw(m2,j)    = kpp_fields%flx_nsw(m2-1,j)
  kpp_fields%flx_nlw(m2,j)    = kpp_fields%flx_nlw(m2-1,j)
  kpp_fields%flx_ice(m2,j)    = kpp_fields%flx_ice(m2-1,j)
  kpp_fields%flx_pcp(m2,j)    = kpp_fields%flx_pcp(m2-1,j)
  kpp_fields%hmix(m2,j)       = kpp_fields%hmix(m2-1,j)
  kpp_fields%freez_flag(m2,j) = kpp_fields%freez_flag(m2-1,j)
  kpp_fields%reset_flag(m2,j) = kpp_fields%reset_flag(m2-1,j)
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
 do i = 1,m2
  do k = 1,mkpp
   kpp_fields%swdk_opt(k,i,1)  = kpp_fields%swdk_opt(k,i,2)
   kpp_fields%U_init(k,i,1)    = kpp_fields%U_init(k,i,2)
   kpp_fields%V_init(k,i,1)    = kpp_fields%V_init(k,i,2)
   kpp_fields%swfrac(k,i,1)    = kpp_fields%swfrac(k,i,2)
   kpp_fields%sal_clim(k,i,1)  = kpp_fields%sal_clim(k,i,2)
   kpp_fields%ocnT_clim(k,i,1) = kpp_fields%ocnT_clim(k,i,2)
   kpp_fields%U(k,i,1)         = kpp_fields%U(k,i,2)
   kpp_fields%V(k,i,1)         = kpp_fields%V(k,i,2)
   kpp_fields%Us0(k,i,1)       = kpp_fields%Us0(k,i,2)
   kpp_fields%Us1(k,i,1)       = kpp_fields%Us1(k,i,2)
   kpp_fields%Vs0(k,i,1)       = kpp_fields%Vs0(k,i,2)
   kpp_fields%Vs1(k,i,1)       = kpp_fields%Vs1(k,i,2)
   kpp_fields%X_t(k,i,1)       = kpp_fields%X_t(k,i,2)
   kpp_fields%X_s(k,i,1)       = kpp_fields%X_s(k,i,2)
   kpp_fields%Xs_t0(k,i,1)     = kpp_fields%Xs_t0(k,i,2)
   kpp_fields%Xs_t1(k,i,1)     = kpp_fields%Xs_t1(k,i,2)
   kpp_fields%Xs_s0(k,i,1)     = kpp_fields%Xs_s0(k,i,2)
   kpp_fields%Xs_s1(k,i,1)     = kpp_fields%Xs_s1(k,i,2)
   kpp_fields%buoy(k,i,1)      = kpp_fields%buoy(k,i,2)
   if(IKPP==2) then
    kpp_fields%wB(k,i,1)         = kpp_fields%wB(k,i,2)
    kpp_fields%wXNTt(k,i,1)      = kpp_fields%wXNTt(k,i,2)
    kpp_fields%wU(k,i,1)         = kpp_fields%wU(k,i,2)
    kpp_fields%wV(k,i,1)         = kpp_fields%wV(k,i,2)
    kpp_fields%wXt(k,i,1)        = kpp_fields%wXt(k,i,2)
    kpp_fields%wXs(k,i,1)        = kpp_fields%wXs(k,i,2)
    kpp_fields%tinc_fcorr(k,i,1) = kpp_fields%tinc_fcorr(k,i,2)
    kpp_fields%sinc_fcorr(k,i,1) = kpp_fields%sinc_fcorr(k,i,2)
    kpp_fields%rho(k,i,1)        = kpp_fields%rho(k,i,2)
    kpp_fields%cp(k,i,1)         = kpp_fields%cp(k,i,2)
   endif
  enddo
  kpp_fields%Sref(i,1)       = kpp_fields%Sref(i,2)
  kpp_fields%jerlov(i,1)     = kpp_fields%jerlov(i,2)
  kpp_fields%ocdepth(i,1)    = kpp_fields%ocdepth(i,2)
  kpp_fields%f(i,1)          = kpp_fields%f(i,2)
  kpp_fields%bottomt(i,1)    = kpp_fields%bottomt(i,2)
  kpp_fields%SST0(i,1)       = kpp_fields%SST0(i,2)
  kpp_fields%old(i,1)        = kpp_fields%old(i,2)
  kpp_fields%new(i,1)        = kpp_fields%new(i,2)
  kpp_fields%flx_ust(i,1)    = kpp_fields%flx_ust(i,2)
  kpp_fields%flx_vst(i,1)    = kpp_fields%flx_vst(i,2)
  kpp_fields%flx_nsw(i,1)    = kpp_fields%flx_nsw(i,2)
  kpp_fields%flx_nlw(i,1)    = kpp_fields%flx_nlw(i,2)
  kpp_fields%flx_ice(i,1)    = kpp_fields%flx_ice(i,2)
  kpp_fields%flx_pcp(i,1)    = kpp_fields%flx_pcp(i,2)
  kpp_fields%hmix(i,1)       = kpp_fields%hmix(i,2)
  kpp_fields%freez_flag(i,1) = kpp_fields%freez_flag(i,2)
  kpp_fields%reset_flag(i,1) = kpp_fields%reset_flag(i,2)
 enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
 do i = 1,m2
  do k = 1,mkpp
   kpp_fields%swdk_opt(k,i,m3)  = kpp_fields%swdk_opt(k,i,m3-1)
   kpp_fields%U_init(k,i,m3)    = kpp_fields%U_init(k,i,m3-1)
   kpp_fields%V_init(k,i,m3)    = kpp_fields%V_init(k,i,m3-1)
   kpp_fields%swfrac(k,i,m3)    = kpp_fields%swfrac(k,i,m3-1)
   kpp_fields%sal_clim(k,i,m3)  = kpp_fields%sal_clim(k,i,m3-1)
   kpp_fields%ocnT_clim(k,i,m3) = kpp_fields%ocnT_clim(k,i,m3-1)
   kpp_fields%U(k,i,m3)         = kpp_fields%U(k,i,m3-1)
   kpp_fields%V(k,i,m3)         = kpp_fields%V(k,i,m3-1)
   kpp_fields%Us0(k,i,m3)       = kpp_fields%Us0(k,i,m3-1)
   kpp_fields%Us1(k,i,m3)       = kpp_fields%Us1(k,i,m3-1)
   kpp_fields%Vs0(k,i,m3)       = kpp_fields%Vs0(k,i,m3-1)
   kpp_fields%Vs1(k,i,m3)       = kpp_fields%Vs1(k,i,m3-1)
   kpp_fields%X_t(k,i,m3)       = kpp_fields%X_t(k,i,m3-1)
   kpp_fields%X_s(k,i,m3)       = kpp_fields%X_s(k,i,m3-1)
   kpp_fields%Xs_t0(k,i,m3)     = kpp_fields%Xs_t0(k,i,m3-1)
   kpp_fields%Xs_t1(k,i,m3)     = kpp_fields%Xs_t1(k,i,m3-1)
   kpp_fields%Xs_s0(k,i,m3)     = kpp_fields%Xs_s0(k,i,m3-1)
   kpp_fields%Xs_s1(k,i,m3)     = kpp_fields%Xs_s1(k,i,m3-1)
   kpp_fields%buoy(k,i,m3)      = kpp_fields%buoy(k,i,m3-1)
   if(IKPP==2) then
    kpp_fields%wB(k,i,m3)         = kpp_fields%wB(k,i,m3-1)
    kpp_fields%wXNTt(k,i,m3)      = kpp_fields%wXNTt(k,i,m3-1)
    kpp_fields%wU(k,i,m3)         = kpp_fields%wU(k,i,m3-1)
    kpp_fields%wV(k,i,m3)         = kpp_fields%wV(k,i,m3-1)
    kpp_fields%wXt(k,i,m3)        = kpp_fields%wXt(k,i,m3-1)
    kpp_fields%wXs(k,i,m3)        = kpp_fields%wXs(k,i,m3-1)
    kpp_fields%tinc_fcorr(k,i,m3) = kpp_fields%tinc_fcorr(k,i,m3-1)
    kpp_fields%sinc_fcorr(k,i,m3) = kpp_fields%sinc_fcorr(k,i,m3-1)
    kpp_fields%rho(k,i,m3)        = kpp_fields%rho(k,i,m3-1)
    kpp_fields%cp(k,i,m3)         = kpp_fields%cp(k,i,m3-1)
   endif
  enddo
  kpp_fields%Sref(i,m3)       = kpp_fields%Sref(i,m3-1)
  kpp_fields%jerlov(i,m3)     = kpp_fields%jerlov(i,m3-1)
  kpp_fields%ocdepth(i,m3)    = kpp_fields%ocdepth(i,m3-1)
  kpp_fields%f(i,m3)          = kpp_fields%f(i,m3-1)
  kpp_fields%bottomt(i,m3)    = kpp_fields%bottomt(i,m3-1)
  kpp_fields%SST0(i,m3)       = kpp_fields%SST0(i,m3-1)
  kpp_fields%old(i,m3)        = kpp_fields%old(i,m3-1)
  kpp_fields%new(i,m3)        = kpp_fields%new(i,m3-1)
  kpp_fields%flx_ust(i,m3)    = kpp_fields%flx_ust(i,m3-1)
  kpp_fields%flx_vst(i,m3)    = kpp_fields%flx_vst(i,m3-1)
  kpp_fields%flx_nsw(i,m3)    = kpp_fields%flx_nsw(i,m3-1)
  kpp_fields%flx_nlw(i,m3)    = kpp_fields%flx_nlw(i,m3-1)
  kpp_fields%flx_ice(i,m3)    = kpp_fields%flx_ice(i,m3-1)
  kpp_fields%flx_pcp(i,m3)    = kpp_fields%flx_pcp(i,m3-1)
  kpp_fields%hmix(i,m3)       = kpp_fields%hmix(i,m3-1)
  kpp_fields%freez_flag(i,m3) = kpp_fields%freez_flag(i,m3-1)
  kpp_fields%reset_flag(i,m3) = kpp_fields%reset_flag(i,m3-1)
 enddo
endif

return
END SUBROUTINE kpp_bcond

!##############################################################################
Subroutine micro_bcond (m1,m2,m3,micro,ibcon)

use mem_grid
use mem_micro
use micphys

implicit none

type (micro_vars) :: micro

integer :: m1,m2,m3,ibcon,i,j,k

!Lateral boundary set
! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do j = 1,m3
  do k = 1,m1
   if(idriz  >= 1) micro%pcpvd(k,1,j) = 0.
   if(irain  >= 1) micro%pcpvr(k,1,j) = 0.
   if(ipris  >= 1) micro%pcpvp(k,1,j) = 0.
   if(isnow  >= 1) micro%pcpvs(k,1,j) = 0.
   if(iaggr  >= 1) micro%pcpva(k,1,j) = 0.
   if(igraup >= 1) micro%pcpvg(k,1,j) = 0.
   if(ihail  >= 1) micro%pcpvh(k,1,j) = 0.
   if(level==3) then
    if(imbudget>=1) then
     micro%latheatvap(k,1,j) = 0.
     micro%latheatfrz(k,1,j) = 0.
     micro%nuccldrt(k,1,j) = 0.
     micro%cld2raint(k,1,j) = 0.
     micro%ice2raint(k,1,j) = 0.
     micro%nucicert(k,1,j) = 0.
     micro%vapliqt(k,1,j) = 0.
     micro%vapicet(k,1,j) = 0.
     micro%evapliqt(k,1,j) = 0.
     micro%evapicet(k,1,j) = 0.
     micro%freezingt(k,1,j) = 0.
     micro%meltingt(k,1,j) = 0.
     micro%melticet(k,1,j) = 0.
     micro%rimecldt(k,1,j) = 0.
     micro%rain2icet(k,1,j) = 0.
     micro%aggregatet(k,1,j) = 0.
     micro%latheatvapt(k,1,j) = 0.
     micro%latheatfrzt(k,1,j) = 0.
    endif
    if(imbudget>=2) then
     micro%inuchomrt(k,1,j) = 0.
     micro%inuccontrt(k,1,j) = 0.
     micro%inucifnrt(k,1,j) = 0.
     micro%inuchazrt(k,1,j) = 0.
     micro%vapcldt(k,1,j) = 0.
     micro%vapraint(k,1,j) = 0.
     micro%vapprist(k,1,j) = 0.
     micro%vapsnowt(k,1,j) = 0.
     micro%vapaggrt(k,1,j) = 0.
     micro%vapgraut(k,1,j) = 0.
     micro%vaphailt(k,1,j) = 0.
     micro%vapdrizt(k,1,j) = 0.
     micro%evapcldt(k,1,j) = 0.
     micro%evapraint(k,1,j) = 0.
     micro%evapprist(k,1,j) = 0.
     micro%evapsnowt(k,1,j) = 0.
     micro%evapaggrt(k,1,j) = 0.
     micro%evapgraut(k,1,j) = 0.
     micro%evaphailt(k,1,j) = 0.
     micro%evapdrizt(k,1,j) = 0.
     micro%meltprist(k,1,j) = 0.
     micro%meltsnowt(k,1,j) = 0.
     micro%meltaggrt(k,1,j) = 0.
     micro%meltgraut(k,1,j) = 0.
     micro%melthailt(k,1,j) = 0.
     micro%rimecldsnowt(k,1,j) = 0.
     micro%rimecldaggrt(k,1,j) = 0.
     micro%rimecldgraut(k,1,j) = 0.
     micro%rimecldhailt(k,1,j) = 0.
     micro%rain2prt(k,1,j) = 0.
     micro%rain2snt(k,1,j) = 0.
     micro%rain2agt(k,1,j) = 0.
     micro%rain2grt(k,1,j) = 0.
     micro%rain2hat(k,1,j) = 0.
     micro%aggrselfprist(k,1,j) = 0.
     micro%aggrselfsnowt(k,1,j) = 0.
     micro%aggrprissnowt(k,1,j) = 0.
    endif
    if(imbudget==3 .and. idust>=1) then
     micro%dust1cldrt(k,1,j) = 0.
     micro%dust2cldrt(k,1,j) = 0.
     micro%dust1drzrt(k,1,j) = 0.
     micro%dust2drzrt(k,1,j) = 0.
    endif
   endif
  enddo
  if(level==3 .and. iccnlev>=2)then
   micro%accpaero(1,j) = 0.
   micro%pcpraero(1,j) = 0.
  endif
  if(idriz  >= 1) micro%pcprd(1,j) = 0.
  if(irain  >= 1) micro%pcprr(1,j) = 0.
  if(ipris  >= 1) micro%pcprp(1,j) = 0.
  if(isnow  >= 1) micro%pcprs(1,j) = 0.
  if(iaggr  >= 1) micro%pcpra(1,j) = 0.
  if(igraup >= 1) micro%pcprg(1,j) = 0.
  if(ihail  >= 1) micro%pcprh(1,j) = 0.
  if(idriz  >= 1) micro%accpd(1,j) = 0.
  if(irain  >= 1) micro%accpr(1,j) = 0.
  if(ipris  >= 1) micro%accpp(1,j) = 0.
  if(isnow  >= 1) micro%accps(1,j) = 0.
  if(iaggr  >= 1) micro%accpa(1,j) = 0.
  if(igraup >= 1) micro%accpg(1,j) = 0.
  if(ihail  >= 1) micro%accph(1,j) = 0.
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do j = 1,m3
  do k = 1,m1
   if(idriz  >= 1) micro%pcpvd(k,m2,j) = 0.
   if(irain  >= 1) micro%pcpvr(k,m2,j) = 0.
   if(ipris  >= 1) micro%pcpvp(k,m2,j) = 0.
   if(isnow  >= 1) micro%pcpvs(k,m2,j) = 0.
   if(iaggr  >= 1) micro%pcpva(k,m2,j) = 0.
   if(igraup >= 1) micro%pcpvg(k,m2,j) = 0.
   if(ihail  >= 1) micro%pcpvh(k,m2,j) = 0.
   if(level==3) then
    if(imbudget>=1) then
     micro%latheatvap(k,m2,j) = 0.
     micro%latheatfrz(k,m2,j) = 0.
     micro%nuccldrt(k,m2,j) = 0.
     micro%cld2raint(k,m2,j) = 0.
     micro%ice2raint(k,m2,j) = 0.
     micro%nucicert(k,m2,j) = 0.
     micro%vapliqt(k,m2,j) = 0.
     micro%vapicet(k,m2,j) = 0.
     micro%evapliqt(k,m2,j) = 0.
     micro%evapicet(k,m2,j) = 0.
     micro%freezingt(k,m2,j) = 0.
     micro%meltingt(k,m2,j) = 0.
     micro%melticet(k,m2,j) = 0.
     micro%rimecldt(k,m2,j) = 0.
     micro%rain2icet(k,m2,j) = 0.
     micro%aggregatet(k,m2,j) = 0.
     micro%latheatvapt(k,m2,j) = 0.
     micro%latheatfrzt(k,m2,j) = 0.
    endif
    if(imbudget>=2) then
     micro%inuchomrt(k,m2,j) = 0.
     micro%inuccontrt(k,m2,j) = 0.
     micro%inucifnrt(k,m2,j) = 0.
     micro%inuchazrt(k,m2,j) = 0.
     micro%vapcldt(k,m2,j) = 0.
     micro%vapraint(k,m2,j) = 0.
     micro%vapprist(k,m2,j) = 0.
     micro%vapsnowt(k,m2,j) = 0.
     micro%vapaggrt(k,m2,j) = 0.
     micro%vapgraut(k,m2,j) = 0.
     micro%vaphailt(k,m2,j) = 0.
     micro%vapdrizt(k,m2,j) = 0.
     micro%evapcldt(k,m2,j) = 0.
     micro%evapraint(k,m2,j) = 0.
     micro%evapprist(k,m2,j) = 0.
     micro%evapsnowt(k,m2,j) = 0.
     micro%evapaggrt(k,m2,j) = 0.
     micro%evapgraut(k,m2,j) = 0.
     micro%evaphailt(k,m2,j) = 0.
     micro%evapdrizt(k,m2,j) = 0.
     micro%meltprist(k,m2,j) = 0.
     micro%meltsnowt(k,m2,j) = 0.
     micro%meltaggrt(k,m2,j) = 0.
     micro%meltgraut(k,m2,j) = 0.
     micro%melthailt(k,m2,j) = 0.
     micro%rimecldsnowt(k,m2,j) = 0.
     micro%rimecldaggrt(k,m2,j) = 0.
     micro%rimecldgraut(k,m2,j) = 0.
     micro%rimecldhailt(k,m2,j) = 0.
     micro%rain2prt(k,m2,j) = 0.
     micro%rain2snt(k,m2,j) = 0.
     micro%rain2agt(k,m2,j) = 0.
     micro%rain2grt(k,m2,j) = 0.
     micro%rain2hat(k,m2,j) = 0.
     micro%aggrselfprist(k,m2,j) = 0.
     micro%aggrselfsnowt(k,m2,j) = 0.
     micro%aggrprissnowt(k,m2,j) = 0.
    endif
    if(imbudget==3 .and. idust>=1) then
     micro%dust1cldrt(k,m2,j) = 0.
     micro%dust2cldrt(k,m2,j) = 0.
     micro%dust1drzrt(k,m2,j) = 0.
     micro%dust2drzrt(k,m2,j) = 0.
    endif
   endif
  enddo
  if(level==3 .and. iccnlev>=2)then
   micro%accpaero(m2,j) = 0.
   micro%pcpraero(m2,j) = 0.
  endif
  if(idriz  >= 1) micro%pcprd(m2,j) = 0.
  if(irain  >= 1) micro%pcprr(m2,j) = 0.
  if(ipris  >= 1) micro%pcprp(m2,j) = 0.
  if(isnow  >= 1) micro%pcprs(m2,j) = 0.
  if(iaggr  >= 1) micro%pcpra(m2,j) = 0.
  if(igraup >= 1) micro%pcprg(m2,j) = 0.
  if(ihail  >= 1) micro%pcprh(m2,j) = 0.
  if(idriz  >= 1) micro%accpd(m2,j) = 0.
  if(irain  >= 1) micro%accpr(m2,j) = 0.
  if(ipris  >= 1) micro%accpp(m2,j) = 0.
  if(isnow  >= 1) micro%accps(m2,j) = 0.
  if(iaggr  >= 1) micro%accpa(m2,j) = 0.
  if(igraup >= 1) micro%accpg(m2,j) = 0.
  if(ihail  >= 1) micro%accph(m2,j) = 0.
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
    do k = 1,m1
     if(idriz  >= 1) micro%pcpvd(k,i,1) = 0.
     if(irain  >= 1) micro%pcpvr(k,i,1) = 0.
     if(ipris  >= 1) micro%pcpvp(k,i,1) = 0.
     if(isnow  >= 1) micro%pcpvs(k,i,1) = 0.
     if(iaggr  >= 1) micro%pcpva(k,i,1) = 0.
     if(igraup >= 1) micro%pcpvg(k,i,1) = 0.
     if(ihail  >= 1) micro%pcpvh(k,i,1) = 0.
     if(level==3) then
      if(imbudget>=1) then
       micro%latheatvap(k,i,1) = 0.
       micro%latheatfrz(k,i,1) = 0.
       micro%nuccldrt(k,i,1) = 0.
       micro%cld2raint(k,i,1) = 0.
       micro%ice2raint(k,i,1) = 0.
       micro%nucicert(k,i,1) = 0.
       micro%vapliqt(k,i,1) = 0.
       micro%vapicet(k,i,1) = 0.
       micro%evapliqt(k,i,1) = 0.
       micro%evapicet(k,i,1) = 0.
       micro%freezingt(k,i,1) = 0.
       micro%meltingt(k,i,1) = 0.
       micro%melticet(k,i,1) = 0.
       micro%rimecldt(k,i,1) = 0.
       micro%rain2icet(k,i,1) = 0.
       micro%aggregatet(k,i,1) = 0.
       micro%latheatvapt(k,i,1) = 0.
       micro%latheatfrzt(k,i,1) = 0.
      endif
      if(imbudget>=2) then
       micro%inuchomrt(k,i,1) = 0.
       micro%inuccontrt(k,i,1) = 0.
       micro%inucifnrt(k,i,1) = 0.
       micro%inuchazrt(k,i,1) = 0.
       micro%vapcldt(k,i,1) = 0.
       micro%vapraint(k,i,1) = 0.
       micro%vapprist(k,i,1) = 0.
       micro%vapsnowt(k,i,1) = 0.
       micro%vapaggrt(k,i,1) = 0.
       micro%vapgraut(k,i,1) = 0.
       micro%vaphailt(k,i,1) = 0.
       micro%vapdrizt(k,i,1) = 0.
       micro%evapcldt(k,i,1) = 0.
       micro%evapraint(k,i,1) = 0.
       micro%evapprist(k,i,1) = 0.
       micro%evapsnowt(k,i,1) = 0.
       micro%evapaggrt(k,i,1) = 0.
       micro%evapgraut(k,i,1) = 0.
       micro%evaphailt(k,i,1) = 0.
       micro%evapdrizt(k,i,1) = 0.
       micro%meltprist(k,i,1) = 0.
       micro%meltsnowt(k,i,1) = 0.
       micro%meltaggrt(k,i,1) = 0.
       micro%meltgraut(k,i,1) = 0.
       micro%melthailt(k,i,1) = 0.
       micro%rimecldsnowt(k,i,1) = 0.
       micro%rimecldaggrt(k,i,1) = 0.
       micro%rimecldgraut(k,i,1) = 0.
       micro%rimecldhailt(k,i,1) = 0.
       micro%rain2prt(k,i,1) = 0.
       micro%rain2snt(k,i,1) = 0.
       micro%rain2agt(k,i,1) = 0.
       micro%rain2grt(k,i,1) = 0.
       micro%rain2hat(k,i,1) = 0.
       micro%aggrselfprist(k,i,1) = 0.
       micro%aggrselfsnowt(k,i,1) = 0.
       micro%aggrprissnowt(k,i,1) = 0.
      endif
      if(imbudget==3 .and. idust>=1) then
       micro%dust1cldrt(k,i,1) = 0.
       micro%dust2cldrt(k,i,1) = 0.
       micro%dust1drzrt(k,i,1) = 0.
       micro%dust2drzrt(k,i,1) = 0.
      endif
     endif
    enddo
    if(level==3 .and. iccnlev>=2)then
     micro%accpaero(i,1) = 0.
     micro%pcpraero(i,1) = 0.
    endif
    if(idriz  >= 1) micro%pcprd(i,1) = 0.
    if(irain  >= 1) micro%pcprr(i,1) = 0.
    if(ipris  >= 1) micro%pcprp(i,1) = 0.
    if(isnow  >= 1) micro%pcprs(i,1) = 0.
    if(iaggr  >= 1) micro%pcpra(i,1) = 0.
    if(igraup >= 1) micro%pcprg(i,1) = 0.
    if(ihail  >= 1) micro%pcprh(i,1) = 0.
    if(idriz  >= 1) micro%accpd(i,1) = 0.
    if(irain  >= 1) micro%accpr(i,1) = 0.
    if(ipris  >= 1) micro%accpp(i,1) = 0.
    if(isnow  >= 1) micro%accps(i,1) = 0.
    if(iaggr  >= 1) micro%accpa(i,1) = 0.
    if(igraup >= 1) micro%accpg(i,1) = 0.
    if(ihail  >= 1) micro%accph(i,1) = 0.
  enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
    do k = 1,m1
     if(idriz  >= 1) micro%pcpvd(k,i,m3) = 0.
     if(irain  >= 1) micro%pcpvr(k,i,m3) = 0.
     if(ipris  >= 1) micro%pcpvp(k,i,m3) = 0.
     if(isnow  >= 1) micro%pcpvs(k,i,m3) = 0.
     if(iaggr  >= 1) micro%pcpva(k,i,m3) = 0.
     if(igraup >= 1) micro%pcpvg(k,i,m3) = 0.
     if(ihail  >= 1) micro%pcpvh(k,i,m3) = 0.
     if(level==3) then
      if(imbudget>=1) then
       micro%latheatvap(k,i,m3) = 0.
       micro%latheatfrz(k,i,m3) = 0.
       micro%nuccldrt(k,i,m3) = 0.
       micro%cld2raint(k,i,m3) = 0.
       micro%ice2raint(k,i,m3) = 0.
       micro%nucicert(k,i,m3) = 0.
       micro%vapliqt(k,i,m3) = 0.
       micro%vapicet(k,i,m3) = 0.
       micro%evapliqt(k,i,m3) = 0.
       micro%evapicet(k,i,m3) = 0.
       micro%freezingt(k,i,m3) = 0.
       micro%meltingt(k,i,m3) = 0.
       micro%melticet(k,i,m3) = 0.
       micro%rimecldt(k,i,m3) = 0.
       micro%rain2icet(k,i,m3) = 0.
       micro%aggregatet(k,i,m3) = 0.
       micro%latheatvapt(k,i,m3) = 0.
       micro%latheatfrzt(k,i,m3) = 0.
      endif
      if(imbudget>=2) then
       micro%inuchomrt(k,i,m3) = 0.
       micro%inuccontrt(k,i,m3) = 0.
       micro%inucifnrt(k,i,m3) = 0.
       micro%inuchazrt(k,i,m3) = 0.
       micro%vapcldt(k,i,m3) = 0.
       micro%vapraint(k,i,m3) = 0.
       micro%vapprist(k,i,m3) = 0.
       micro%vapsnowt(k,i,m3) = 0.
       micro%vapaggrt(k,i,m3) = 0.
       micro%vapgraut(k,i,m3) = 0.
       micro%vaphailt(k,i,m3) = 0.
       micro%vapdrizt(k,i,m3) = 0.
       micro%evapcldt(k,i,m3) = 0.
       micro%evapraint(k,i,m3) = 0.
       micro%evapprist(k,i,m3) = 0.
       micro%evapsnowt(k,i,m3) = 0.
       micro%evapaggrt(k,i,m3) = 0.
       micro%evapgraut(k,i,m3) = 0.
       micro%evaphailt(k,i,m3) = 0.
       micro%evapdrizt(k,i,m3) = 0.
       micro%meltprist(k,i,m3) = 0.
       micro%meltsnowt(k,i,m3) = 0.
       micro%meltaggrt(k,i,m3) = 0.
       micro%meltgraut(k,i,m3) = 0.
       micro%melthailt(k,i,m3) = 0.
       micro%rimecldsnowt(k,i,m3) = 0.
       micro%rimecldaggrt(k,i,m3) = 0.
       micro%rimecldgraut(k,i,m3) = 0.
       micro%rimecldhailt(k,i,m3) = 0.
       micro%rain2prt(k,i,m3) = 0.
       micro%rain2snt(k,i,m3) = 0.
       micro%rain2agt(k,i,m3) = 0.
       micro%rain2grt(k,i,m3) = 0.
       micro%rain2hat(k,i,m3) = 0.
       micro%aggrselfprist(k,i,m3) = 0.
       micro%aggrselfsnowt(k,i,m3) = 0.
       micro%aggrprissnowt(k,i,m3) = 0.
      endif
      if(imbudget==3 .and. idust>=1) then
       micro%dust1cldrt(k,i,m3) = 0.
       micro%dust2cldrt(k,i,m3) = 0.
       micro%dust1drzrt(k,i,m3) = 0.
       micro%dust2drzrt(k,i,m3) = 0.
      endif
     endif
    enddo
    if(level==3 .and. iccnlev>=2)then
     micro%accpaero(i,m3) = 0.
     micro%pcpraero(i,m3) = 0.
    endif
    if(idriz  >= 1) micro%pcprd(i,m3) = 0.
    if(irain  >= 1) micro%pcprr(i,m3) = 0.
    if(ipris  >= 1) micro%pcprp(i,m3) = 0.
    if(isnow  >= 1) micro%pcprs(i,m3) = 0.
    if(iaggr  >= 1) micro%pcpra(i,m3) = 0.
    if(igraup >= 1) micro%pcprg(i,m3) = 0.
    if(ihail  >= 1) micro%pcprh(i,m3) = 0.
    if(idriz  >= 1) micro%accpd(i,m3) = 0.
    if(irain  >= 1) micro%accpr(i,m3) = 0.
    if(ipris  >= 1) micro%accpp(i,m3) = 0.
    if(isnow  >= 1) micro%accps(i,m3) = 0.
    if(iaggr  >= 1) micro%accpa(i,m3) = 0.
    if(igraup >= 1) micro%accpg(i,m3) = 0.
    if(ihail  >= 1) micro%accph(i,m3) = 0.
  enddo
endif

return
END SUBROUTINE micro_bcond

!##############################################################################
Subroutine leaf_bcond (m2,m3,mzg,mzs,npat,leaf,ibcon)

use mem_grid
use mem_leaf

implicit none

type (leaf_vars) :: leaf

integer :: m2,m3,mzg,mzs,ibcon,npat,ipat,i,j,k

!Set BCs for Leaf variables. Do not include veg_ndvif, leaf_class,
! patch_area, and soil_text since these are read from the surface
! files and are boundary set when files are read in by the model.

! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do ipat = 1,npat
   do j = 1,m3
      leaf%ustar          (1,j,ipat) = leaf%ustar            (2,j,ipat)
      leaf%tstar          (1,j,ipat) = leaf%tstar            (2,j,ipat)
      leaf%rstar          (1,j,ipat) = leaf%rstar            (2,j,ipat)
      leaf%veg_albedo     (1,j,ipat) = leaf%veg_albedo       (2,j,ipat)
      leaf%veg_fracarea   (1,j,ipat) = leaf%veg_fracarea     (2,j,ipat)
      leaf%veg_lai        (1,j,ipat) = leaf%veg_lai          (2,j,ipat)
      leaf%veg_tai        (1,j,ipat) = leaf%veg_tai          (2,j,ipat)
      leaf%veg_rough      (1,j,ipat) = leaf%veg_rough        (2,j,ipat)
      leaf%veg_height     (1,j,ipat) = leaf%veg_height       (2,j,ipat)
      leaf%patch_rough    (1,j,ipat) = leaf%patch_rough      (2,j,ipat)
      leaf%soil_rough     (1,j,ipat) = leaf%soil_rough       (2,j,ipat)
      leaf%sfcwater_nlev  (1,j,ipat) = leaf%sfcwater_nlev    (2,j,ipat)
      leaf%stom_resist    (1,j,ipat) = leaf%stom_resist      (2,j,ipat)
      leaf%ground_rsat    (1,j,ipat) = leaf%ground_rsat      (2,j,ipat)
      leaf%ground_rvap    (1,j,ipat) = leaf%ground_rvap      (2,j,ipat)
      leaf%veg_water      (1,j,ipat) = leaf%veg_water        (2,j,ipat)
      leaf%veg_temp       (1,j,ipat) = leaf%veg_temp         (2,j,ipat)
      leaf%can_rvap       (1,j,ipat) = leaf%can_rvap         (2,j,ipat)
      leaf%can_temp       (1,j,ipat) = leaf%can_temp         (2,j,ipat)
      leaf%veg_ndvip      (1,j,ipat) = leaf%veg_ndvip        (2,j,ipat)
      leaf%veg_ndvic      (1,j,ipat) = leaf%veg_ndvic        (2,j,ipat)
   
      do k = 1,mzg
       leaf%soil_water       (k,1,j,ipat) = leaf%soil_water      (k,2,j,ipat)
       leaf%soil_energy      (k,1,j,ipat) = leaf%soil_energy     (k,2,j,ipat)
      enddo

      do k = 1,mzs
       leaf%sfcwater_mass    (k,1,j,ipat) = leaf%sfcwater_mass   (k,2,j,ipat)
       leaf%sfcwater_energy  (k,1,j,ipat) = leaf%sfcwater_energy (k,2,j,ipat)
       leaf%sfcwater_depth   (k,1,j,ipat) = leaf%sfcwater_depth  (k,2,j,ipat)
      enddo
   enddo   
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do ipat = 1,npat
   do j = 1,m3
      leaf%ustar         (m2,j,ipat) = leaf%ustar         (m2-1,j,ipat)
      leaf%tstar         (m2,j,ipat) = leaf%tstar         (m2-1,j,ipat)
      leaf%rstar         (m2,j,ipat) = leaf%rstar         (m2-1,j,ipat)
      leaf%veg_albedo    (m2,j,ipat) = leaf%veg_albedo    (m2-1,j,ipat)
      leaf%veg_fracarea  (m2,j,ipat) = leaf%veg_fracarea  (m2-1,j,ipat)
      leaf%veg_lai       (m2,j,ipat) = leaf%veg_lai       (m2-1,j,ipat)
      leaf%veg_tai       (m2,j,ipat) = leaf%veg_tai       (m2-1,j,ipat)
      leaf%veg_rough     (m2,j,ipat) = leaf%veg_rough     (m2-1,j,ipat)
      leaf%veg_height    (m2,j,ipat) = leaf%veg_height    (m2-1,j,ipat)
      leaf%patch_rough   (m2,j,ipat) = leaf%patch_rough   (m2-1,j,ipat)
      leaf%soil_rough    (m2,j,ipat) = leaf%soil_rough    (m2-1,j,ipat)
      leaf%sfcwater_nlev (m2,j,ipat) = leaf%sfcwater_nlev (m2-1,j,ipat)
      leaf%stom_resist   (m2,j,ipat) = leaf%stom_resist   (m2-1,j,ipat)
      leaf%ground_rsat   (m2,j,ipat) = leaf%ground_rsat   (m2-1,j,ipat)
      leaf%ground_rvap   (m2,j,ipat) = leaf%ground_rvap   (m2-1,j,ipat)
      leaf%veg_water     (m2,j,ipat) = leaf%veg_water     (m2-1,j,ipat)
      leaf%veg_temp      (m2,j,ipat) = leaf%veg_temp      (m2-1,j,ipat)
      leaf%can_rvap      (m2,j,ipat) = leaf%can_rvap      (m2-1,j,ipat)
      leaf%can_temp      (m2,j,ipat) = leaf%can_temp      (m2-1,j,ipat)
      leaf%veg_ndvip     (m2,j,ipat) = leaf%veg_ndvip     (m2-1,j,ipat)
      leaf%veg_ndvic     (m2,j,ipat) = leaf%veg_ndvic     (m2-1,j,ipat)
   
      do k = 1,mzg
       leaf%soil_water    (k,m2,j,ipat) = leaf%soil_water   (k,m2-1,j,ipat)
       leaf%soil_energy   (k,m2,j,ipat) = leaf%soil_energy  (k,m2-1,j,ipat)
      enddo

      do k = 1,mzs
       leaf%sfcwater_mass  (k,m2,j,ipat) = leaf%sfcwater_mass  (k,m2-1,j,ipat)
       leaf%sfcwater_energy(k,m2,j,ipat) = leaf%sfcwater_energy(k,m2-1,j,ipat)
       leaf%sfcwater_depth (k,m2,j,ipat) = leaf%sfcwater_depth (k,m2-1,j,ipat)
      enddo
   enddo   
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
  do ipat = 1,npat
      do i = 1,m2
         leaf%ustar          (i,1,ipat) = leaf%ustar            (i,2,ipat)
         leaf%tstar          (i,1,ipat) = leaf%tstar            (i,2,ipat)
         leaf%rstar          (i,1,ipat) = leaf%rstar            (i,2,ipat)
         leaf%veg_albedo     (i,1,ipat) = leaf%veg_albedo       (i,2,ipat)
         leaf%veg_fracarea   (i,1,ipat) = leaf%veg_fracarea     (i,2,ipat)
         leaf%veg_lai        (i,1,ipat) = leaf%veg_lai          (i,2,ipat)
         leaf%veg_tai        (i,1,ipat) = leaf%veg_tai          (i,2,ipat)
         leaf%veg_rough      (i,1,ipat) = leaf%veg_rough        (i,2,ipat)
         leaf%veg_height     (i,1,ipat) = leaf%veg_height       (i,2,ipat)
         leaf%patch_rough    (i,1,ipat) = leaf%patch_rough      (i,2,ipat)
         leaf%soil_rough     (i,1,ipat) = leaf%soil_rough       (i,2,ipat)
         leaf%sfcwater_nlev  (i,1,ipat) = leaf%sfcwater_nlev    (i,2,ipat)
         leaf%stom_resist    (i,1,ipat) = leaf%stom_resist      (i,2,ipat)
         leaf%ground_rsat    (i,1,ipat) = leaf%ground_rsat      (i,2,ipat)
         leaf%ground_rvap    (i,1,ipat) = leaf%ground_rvap      (i,2,ipat)
         leaf%veg_water      (i,1,ipat) = leaf%veg_water        (i,2,ipat)
         leaf%veg_temp       (i,1,ipat) = leaf%veg_temp         (i,2,ipat)
         leaf%can_rvap       (i,1,ipat) = leaf%can_rvap         (i,2,ipat)
         leaf%can_temp       (i,1,ipat) = leaf%can_temp         (i,2,ipat)
         leaf%veg_ndvip      (i,1,ipat) = leaf%veg_ndvip        (i,2,ipat)
         leaf%veg_ndvic      (i,1,ipat) = leaf%veg_ndvic        (i,2,ipat)
   
         do k = 1,mzg
          leaf%soil_water     (k,i,1,ipat) = leaf%soil_water      (k,i,2,ipat)
          leaf%soil_energy    (k,i,1,ipat) = leaf%soil_energy     (k,i,2,ipat)
         enddo

         do k = 1,mzs
          leaf%sfcwater_mass  (k,i,1,ipat) = leaf%sfcwater_mass   (k,i,2,ipat)
          leaf%sfcwater_energy(k,i,1,ipat) = leaf%sfcwater_energy (k,i,2,ipat)
          leaf%sfcwater_depth (k,i,1,ipat) = leaf%sfcwater_depth  (k,i,2,ipat)
         enddo
      enddo   
  enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
   do ipat = 1,npat
      do i = 1,m2
         leaf%ustar         (i,m3,ipat) = leaf%ustar         (i,m3-1,ipat)
         leaf%tstar         (i,m3,ipat) = leaf%tstar         (i,m3-1,ipat)
         leaf%rstar         (i,m3,ipat) = leaf%rstar         (i,m3-1,ipat)
         leaf%veg_albedo    (i,m3,ipat) = leaf%veg_albedo    (i,m3-1,ipat)
         leaf%veg_fracarea  (i,m3,ipat) = leaf%veg_fracarea  (i,m3-1,ipat)
         leaf%veg_lai       (i,m3,ipat) = leaf%veg_lai       (i,m3-1,ipat)
         leaf%veg_tai       (i,m3,ipat) = leaf%veg_tai       (i,m3-1,ipat)
         leaf%veg_rough     (i,m3,ipat) = leaf%veg_rough     (i,m3-1,ipat)
         leaf%veg_height    (i,m3,ipat) = leaf%veg_height    (i,m3-1,ipat)
         leaf%patch_rough   (i,m3,ipat) = leaf%patch_rough   (i,m3-1,ipat)
         leaf%soil_rough    (i,m3,ipat) = leaf%soil_rough    (i,m3-1,ipat)
         leaf%sfcwater_nlev (i,m3,ipat) = leaf%sfcwater_nlev (i,m3-1,ipat)
         leaf%stom_resist   (i,m3,ipat) = leaf%stom_resist   (i,m3-1,ipat)
         leaf%ground_rsat   (i,m3,ipat) = leaf%ground_rsat   (i,m3-1,ipat)
         leaf%ground_rvap   (i,m3,ipat) = leaf%ground_rvap   (i,m3-1,ipat)
         leaf%veg_water     (i,m3,ipat) = leaf%veg_water     (i,m3-1,ipat)
         leaf%veg_temp      (i,m3,ipat) = leaf%veg_temp      (i,m3-1,ipat)
         leaf%can_rvap      (i,m3,ipat) = leaf%can_rvap      (i,m3-1,ipat)
         leaf%can_temp      (i,m3,ipat) = leaf%can_temp      (i,m3-1,ipat)
         leaf%veg_ndvip     (i,m3,ipat) = leaf%veg_ndvip     (i,m3-1,ipat)
         leaf%veg_ndvic     (i,m3,ipat) = leaf%veg_ndvic     (i,m3-1,ipat)
   
         do k = 1,mzg
          leaf%soil_water     (k,i,m3,ipat)=leaf%soil_water   (k,i,m3-1,ipat)
          leaf%soil_energy    (k,i,m3,ipat)=leaf%soil_energy  (k,i,m3-1,ipat)
         enddo

         do k = 1,mzs
          leaf%sfcwater_mass  (k,i,m3,ipat)=leaf%sfcwater_mass  (k,i,m3-1,ipat)
          leaf%sfcwater_energy(k,i,m3,ipat)=leaf%sfcwater_energy(k,i,m3-1,ipat)
          leaf%sfcwater_depth (k,i,m3,ipat)=leaf%sfcwater_depth (k,i,m3-1,ipat)
         enddo
      enddo   
   enddo
endif

return
END SUBROUTINE leaf_bcond

!##############################################################################
Subroutine rad_bcond (m1,m2,m3,radiate,ibcon)

use mem_grid
use mem_radiate

implicit none

type (radiate_vars) :: radiate

integer :: m1,m2,m3,ibcon,i,j,k

!Lateral boundary set
! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do j = 1,m3
  do k = 1,m1
   radiate%fthrd(k,1,j)  = radiate%fthrd(k,2,j)
   if(ilwrtyp == 3 .or. iswrtyp == 3) then
     radiate%bext(k,1,j)   = radiate%bext(k,2,j)
   endif
   if(iswrtyp == 3) then
     radiate%swup(k,1,j) = radiate%swup(k,2,j)
     radiate%swdn(k,1,j) = radiate%swdn(k,2,j)
   endif
   if(ilwrtyp == 3) then
     radiate%lwup(k,1,j) = radiate%lwup(k,2,j)
     radiate%lwdn(k,1,j) = radiate%lwdn(k,2,j)
   endif
  enddo
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do j = 1,m3
  do k = 1,m1
   radiate%fthrd(k,m2,j) = radiate%fthrd(k,m2-1,j)
   if(ilwrtyp == 3 .or. iswrtyp == 3) then
     radiate%bext(k,m2,j)  = radiate%bext(k,m2-1,j)
   endif
   if(iswrtyp == 3) then
     radiate%swup(k,m2,j) = radiate%swup(k,m2-1,j)
     radiate%swdn(k,m2,j) = radiate%swdn(k,m2-1,j)
   endif
   if(ilwrtyp == 3) then
     radiate%lwup(k,m2,j) = radiate%lwup(k,m2-1,j)
     radiate%lwdn(k,m2,j) = radiate%lwdn(k,m2-1,j)
   endif
  enddo
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
   do k = 1,m1
     radiate%fthrd(k,i,1)   = radiate%fthrd(k,i,2)
     if(ilwrtyp == 3 .or. iswrtyp == 3) then
       radiate%bext(k,i,1)    = radiate%bext(k,i,2)
     endif
     if(iswrtyp == 3) then
       radiate%swup(k,i,1) = radiate%swup(k,i,2)
       radiate%swdn(k,i,1) = radiate%swdn(k,i,2)
     endif
     if(ilwrtyp == 3) then
       radiate%lwup(k,i,1) = radiate%lwup(k,i,2)
       radiate%lwdn(k,i,1) = radiate%lwdn(k,i,2)
     endif
   enddo
  enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
   do k = 1,m1
     radiate%fthrd(k,i,m3)  = radiate%fthrd(k,i,m3-1)
     if(ilwrtyp == 3 .or. iswrtyp == 3) then
       radiate%bext(k,i,m3)   = radiate%bext(k,i,m3-1)
     endif
     if(iswrtyp == 3) then
       radiate%swup(k,i,m3) = radiate%swup(k,i,m3-1)
       radiate%swdn(k,i,m3) = radiate%swdn(k,i,m3-1)
     endif
     if(ilwrtyp == 3) then
       radiate%lwup(k,i,m3) = radiate%lwup(k,i,m3-1)
       radiate%lwdn(k,i,m3) = radiate%lwdn(k,i,m3-1)
     endif
   enddo
  enddo
endif

return
END SUBROUTINE rad_bcond

!##############################################################################
Subroutine turb_bcond (m2,m3,turb,ibcon)

use mem_grid
use mem_turb

implicit none

type (turb_vars) :: turb

integer :: m2,m3,ibcon,i,j

! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do j = 1,m3
  turb%sflux_u(1,j) = turb%sflux_u(2,j)
  turb%sflux_v(1,j) = turb%sflux_v(2,j)
  turb%sflux_w(1,j) = turb%sflux_w(2,j)
  turb%sflux_t(1,j) = turb%sflux_t(2,j)
  turb%sflux_r(1,j) = turb%sflux_r(2,j)
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do j = 1,m3
  turb%sflux_u(m2,j) = turb%sflux_u(m2-1,j)
  turb%sflux_v(m2,j) = turb%sflux_v(m2-1,j)
  turb%sflux_w(m2,j) = turb%sflux_w(m2-1,j)
  turb%sflux_t(m2,j) = turb%sflux_t(m2-1,j)
  turb%sflux_r(m2,j) = turb%sflux_r(m2-1,j)
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
     turb%sflux_u(i,1) = turb%sflux_u(i,2)
     turb%sflux_v(i,1) = turb%sflux_v(i,2)
     turb%sflux_w(i,1) = turb%sflux_w(i,2)
     turb%sflux_t(i,1) = turb%sflux_t(i,2)
     turb%sflux_r(i,1) = turb%sflux_r(i,2)
  enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
     turb%sflux_u(i,m3) = turb%sflux_u(i,m3-1)
     turb%sflux_v(i,m3) = turb%sflux_v(i,m3-1)
     turb%sflux_w(i,m3) = turb%sflux_w(i,m3-1)
     turb%sflux_t(i,m3) = turb%sflux_t(i,m3-1)
     turb%sflux_r(i,m3) = turb%sflux_r(i,m3-1)
  enddo
endif

return
END SUBROUTINE turb_bcond

!##############################################################################
Subroutine sfcrad_bcond (m2,m3,radiate,ibcon)

use mem_grid
use mem_radiate

implicit none

type (radiate_vars) :: radiate

integer :: m2,m3,ibcon,i,j

! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do j = 1,m3
  radiate%rshort(1,j)  = radiate%rshort(2,j)
  radiate%rlong(1,j)   = radiate%rlong(2,j)
  radiate%rlongup(1,j) = radiate%rlongup(2,j)
  radiate%albedt(1,j)  = radiate%albedt(2,j)
  radiate%cosz(1,j)    = radiate%cosz(2,j)
  radiate%aodt(1,j)    = radiate%aodt(2,j)
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do j = 1,m3
  radiate%rshort(m2,j)  = radiate%rshort(m2-1,j)
  radiate%rlong(m2,j)   = radiate%rlong(m2-1,j)
  radiate%rlongup(m2,j) = radiate%rlongup(m2-1,j)
  radiate%albedt(m2,j)  = radiate%albedt(m2-1,j)
  radiate%cosz(m2,j)    = radiate%cosz(m2-1,j)
  radiate%aodt(m2,j)    = radiate%aodt(m2-1,j)
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
     radiate%rshort(i,1)  = radiate%rshort(i,2)
     radiate%rlong(i,1)   = radiate%rlong(i,2)
     radiate%rlongup(i,1) = radiate%rlongup(i,2)
     radiate%albedt(i,1)  = radiate%albedt(i,2)
     radiate%cosz(i,1)    = radiate%cosz(i,2)
     radiate%aodt(i,1)    = radiate%aodt(i,2)
   enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
  do i = 1,m2
     radiate%rshort(i,m3)  = radiate%rshort(i,m3-1)
     radiate%rlong(i,m3)   = radiate%rlong(i,m3-1)
     radiate%rlongup(i,m3) = radiate%rlongup(i,m3-1)
     radiate%albedt(i,m3)  = radiate%albedt(i,m3-1)
     radiate%cosz(i,m3)    = radiate%cosz(i,m3-1)
     radiate%aodt(i,m3)    = radiate%aodt(i,m3-1)
   enddo
endif

return
END SUBROUTINE sfcrad_bcond

!##############################################################################
Subroutine sib_bcond (m2,m3,npat,sib,ibcon)

use mem_grid
use mem_sib

implicit none

type (sib_vars) :: sib

integer :: m2,m3,ibcon,npat,ipat,i,j

!Set BCs for SiB variables.

! West side, ibcon 1's place bit set
if (iand(ibcon,1) .ne. 0) then
 do ipat = 1,npat
   do j = 1,m3
      sib%snow1(1,j,ipat)   = sib%snow1(2,j,ipat)
      sib%snow2(1,j,ipat)   = sib%snow2(2,j,ipat)
      sib%capac1(1,j,ipat)  = sib%capac1(2,j,ipat)
      sib%capac2(1,j,ipat)  = sib%capac2(2,j,ipat)
      sib%pco2ap(1,j,ipat)  = sib%pco2ap(2,j,ipat)
      sib%co2flx(1,j,ipat)  = sib%co2flx(2,j,ipat)
      sib%sfcswa(1,j,ipat)  = sib%sfcswa(2,j,ipat)
      sib%uplwrf(1,j,ipat)  = sib%uplwrf(2,j,ipat)
      sib%assimn(1,j,ipat)  = sib%assimn(2,j,ipat)
      sib%respg(1,j,ipat)   = sib%respg(2,j,ipat)
      sib%rstfac1(1,j,ipat) = sib%rstfac1(2,j,ipat)
      sib%rstfac2(1,j,ipat) = sib%rstfac2(2,j,ipat)
      sib%rstfac3(1,j,ipat) = sib%rstfac3(2,j,ipat)
      sib%ect(1,j,ipat)     = sib%ect(2,j,ipat)
      sib%eci(1,j,ipat)     = sib%eci(2,j,ipat)
      sib%egi(1,j,ipat)     = sib%egi(2,j,ipat)
      sib%egs(1,j,ipat)     = sib%egs(2,j,ipat)
      sib%hc(1,j,ipat)      = sib%hc(2,j,ipat)
      sib%hg(1,j,ipat)      = sib%hg(2,j,ipat)
      sib%ra(1,j,ipat)      = sib%ra(2,j,ipat)
      sib%rb(1,j,ipat)      = sib%rb(2,j,ipat)
      sib%rc(1,j,ipat)      = sib%rc(2,j,ipat)
      sib%rd(1,j,ipat)      = sib%rd(2,j,ipat)
      sib%roff(1,j,ipat)    = sib%roff(2,j,ipat)
      sib%green(1,j,ipat)   = sib%green(2,j,ipat)
      sib%apar(1,j,ipat)    = sib%apar(2,j,ipat)
      sib%ventmf(1,j,ipat)  = sib%ventmf(2,j,ipat)
      sib%pco2c(1,j,ipat)   = sib%pco2c(2,j,ipat)
      sib%pco2i(1,j,ipat)   = sib%pco2i(2,j,ipat)
      sib%pco2s(1,j,ipat)   = sib%pco2s(2,j,ipat)
      sib%pco2m(1,j,ipat)   = sib%pco2m(2,j,ipat)
      sib%ea(1,j,ipat)      = sib%ea(2,j,ipat)
      sib%em(1,j,ipat)      = sib%em(2,j,ipat)
      sib%rha(1,j,ipat)     = sib%rha(2,j,ipat)
      sib%radvbc(1,j,ipat)  = sib%radvbc(2,j,ipat)
      sib%radvdc(1,j,ipat)  = sib%radvdc(2,j,ipat)
      sib%radnbc(1,j,ipat)  = sib%radnbc(2,j,ipat)
      sib%radndc(1,j,ipat)  = sib%radndc(2,j,ipat)
      sib%psy(1,j,ipat)     = sib%psy(2,j,ipat)
   enddo   
 enddo
endif

! East side, ibcon 2's place bit set
if (iand(ibcon,2) .ne. 0) then
 do ipat = 1,npat
   do j = 1,m3
      sib%snow1(m2,j,ipat)   = sib%snow1(m2-1,j,ipat)
      sib%snow2(m2,j,ipat)   = sib%snow2(m2-1,j,ipat)
      sib%capac1(m2,j,ipat)  = sib%capac1(m2-1,j,ipat)
      sib%capac2(m2,j,ipat)  = sib%capac2(m2-1,j,ipat)
      sib%pco2ap(m2,j,ipat)  = sib%pco2ap(m2-1,j,ipat)
      sib%co2flx(m2,j,ipat)  = sib%co2flx(m2-1,j,ipat)
      sib%sfcswa(m2,j,ipat)  = sib%sfcswa(m2-1,j,ipat)
      sib%uplwrf(m2,j,ipat)  = sib%uplwrf(m2-1,j,ipat)
      sib%assimn(m2,j,ipat)  = sib%assimn(m2-1,j,ipat)
      sib%respg(m2,j,ipat)   = sib%respg(m2-1,j,ipat)
      sib%rstfac1(m2,j,ipat) = sib%rstfac1(m2-1,j,ipat)
      sib%rstfac2(m2,j,ipat) = sib%rstfac2(m2-1,j,ipat)
      sib%rstfac3(m2,j,ipat) = sib%rstfac3(m2-1,j,ipat)
      sib%ect(m2,j,ipat)     = sib%ect(m2-1,j,ipat)
      sib%eci(m2,j,ipat)     = sib%eci(m2-1,j,ipat)
      sib%egi(m2,j,ipat)     = sib%egi(m2-1,j,ipat)
      sib%egs(m2,j,ipat)     = sib%egs(m2-1,j,ipat)
      sib%hc(m2,j,ipat)      = sib%hc(m2-1,j,ipat)
      sib%hg(m2,j,ipat)      = sib%hg(m2-1,j,ipat)
      sib%ra(m2,j,ipat)      = sib%ra(m2-1,j,ipat)
      sib%rb(m2,j,ipat)      = sib%rb(m2-1,j,ipat)
      sib%rc(m2,j,ipat)      = sib%rc(m2-1,j,ipat)
      sib%rd(m2,j,ipat)      = sib%rd(m2-1,j,ipat)
      sib%roff(m2,j,ipat)    = sib%roff(m2-1,j,ipat)
      sib%green(m2,j,ipat)   = sib%green(m2-1,j,ipat)
      sib%apar(m2,j,ipat)    = sib%apar(m2-1,j,ipat)
      sib%ventmf(m2,j,ipat)  = sib%ventmf(m2-1,j,ipat)
      sib%pco2c(m2,j,ipat)   = sib%pco2c(m2-1,j,ipat)
      sib%pco2i(m2,j,ipat)   = sib%pco2i(m2-1,j,ipat)
      sib%pco2s(m2,j,ipat)   = sib%pco2s(m2-1,j,ipat)
      sib%pco2m(m2,j,ipat)   = sib%pco2m(m2-1,j,ipat)
      sib%ea(m2,j,ipat)      = sib%ea(m2-1,j,ipat)
      sib%em(m2,j,ipat)      = sib%em(m2-1,j,ipat)
      sib%rha(m2,j,ipat)     = sib%rha(m2-1,j,ipat)
      sib%radvbc(m2,j,ipat)  = sib%radvbc(m2-1,j,ipat)
      sib%radvdc(m2,j,ipat)  = sib%radvdc(m2-1,j,ipat)
      sib%radnbc(m2,j,ipat)  = sib%radnbc(m2-1,j,ipat)
      sib%radndc(m2,j,ipat)  = sib%radndc(m2-1,j,ipat)
      sib%psy(m2,j,ipat)     = sib%psy(m2-1,j,ipat)
   enddo   
 enddo
endif

! South side, ibcon 4's place bit set
if ((iand(ibcon,4) .ne. 0) .and. (jdim .eq. 1)) then
 do ipat = 1,npat
   do i = 1,m2
      sib%snow1(i,1,ipat)   = sib%snow1(i,2,ipat)
      sib%snow2(i,1,ipat)   = sib%snow2(i,2,ipat)
      sib%capac1(i,1,ipat)  = sib%capac1(i,2,ipat)
      sib%capac2(i,1,ipat)  = sib%capac2(i,2,ipat)
      sib%pco2ap(i,1,ipat)  = sib%pco2ap(i,2,ipat)
      sib%co2flx(i,1,ipat)  = sib%co2flx(i,2,ipat)
      sib%sfcswa(i,1,ipat)  = sib%sfcswa(i,2,ipat)
      sib%uplwrf(i,1,ipat)  = sib%uplwrf(i,2,ipat)
      sib%assimn(i,1,ipat)  = sib%assimn(i,2,ipat)
      sib%respg(i,1,ipat)   = sib%respg(i,2,ipat)
      sib%rstfac1(i,1,ipat) = sib%rstfac1(i,2,ipat)
      sib%rstfac2(i,1,ipat) = sib%rstfac2(i,2,ipat)
      sib%rstfac3(i,1,ipat) = sib%rstfac3(i,2,ipat)
      sib%ect(i,1,ipat)     = sib%ect(i,2,ipat)
      sib%eci(i,1,ipat)     = sib%eci(i,2,ipat)
      sib%egi(i,1,ipat)     = sib%egi(i,2,ipat)
      sib%egs(i,1,ipat)     = sib%egs(i,2,ipat)
      sib%hc(i,1,ipat)      = sib%hc(i,2,ipat)
      sib%hg(i,1,ipat)      = sib%hg(i,2,ipat)
      sib%ra(i,1,ipat)      = sib%ra(i,2,ipat)
      sib%rb(i,1,ipat)      = sib%rb(i,2,ipat)
      sib%rc(i,1,ipat)      = sib%rc(i,2,ipat)
      sib%rd(i,1,ipat)      = sib%rd(i,2,ipat)
      sib%roff(i,1,ipat)    = sib%roff(i,2,ipat)
      sib%green(i,1,ipat)   = sib%green(i,2,ipat)
      sib%apar(i,1,ipat)    = sib%apar(i,2,ipat)
      sib%ventmf(i,1,ipat)  = sib%ventmf(i,2,ipat)
      sib%pco2c(i,1,ipat)   = sib%pco2c(i,2,ipat)
      sib%pco2i(i,1,ipat)   = sib%pco2i(i,2,ipat)
      sib%pco2s(i,1,ipat)   = sib%pco2s(i,2,ipat)
      sib%pco2m(i,1,ipat)   = sib%pco2m(i,2,ipat)
      sib%ea(i,1,ipat)      = sib%ea(i,2,ipat)
      sib%em(i,1,ipat)      = sib%em(i,2,ipat)
      sib%rha(i,1,ipat)     = sib%rha(i,2,ipat)
      sib%radvbc(i,1,ipat)  = sib%radvbc(i,2,ipat)
      sib%radvdc(i,1,ipat)  = sib%radvdc(i,2,ipat)
      sib%radnbc(i,1,ipat)  = sib%radnbc(i,2,ipat)
      sib%radndc(i,1,ipat)  = sib%radndc(i,2,ipat)
      sib%psy(i,1,ipat)     = sib%psy(i,2,ipat)
   enddo   
 enddo
endif

! North side, ibcon 8's place bit set
if ((iand(ibcon,8) .ne. 0) .and. (jdim .eq. 1)) then
 do ipat = 1,npat
   do i = 1,m2
      sib%snow1(i,m3,ipat)   = sib%snow1(i,m3-1,ipat)
      sib%snow2(i,m3,ipat)   = sib%snow2(i,m3-1,ipat)
      sib%capac1(i,m3,ipat)  = sib%capac1(i,m3-1,ipat)
      sib%capac2(i,m3,ipat)  = sib%capac2(i,m3-1,ipat)
      sib%pco2ap(i,m3,ipat)  = sib%pco2ap(i,m3-1,ipat)
      sib%co2flx(i,m3,ipat)  = sib%co2flx(i,m3-1,ipat)
      sib%sfcswa(i,m3,ipat)  = sib%sfcswa(i,m3-1,ipat)
      sib%uplwrf(i,m3,ipat)  = sib%uplwrf(i,m3-1,ipat)
      sib%assimn(i,m3,ipat)  = sib%assimn(i,m3-1,ipat)
      sib%respg(i,m3,ipat)   = sib%respg(i,m3-1,ipat)
      sib%rstfac1(i,m3,ipat) = sib%rstfac1(i,m3-1,ipat)
      sib%rstfac2(i,m3,ipat) = sib%rstfac2(i,m3-1,ipat)
      sib%rstfac3(i,m3,ipat) = sib%rstfac3(i,m3-1,ipat)
      sib%ect(i,m3,ipat)     = sib%ect(i,m3-1,ipat)
      sib%eci(i,m3,ipat)     = sib%eci(i,m3-1,ipat)
      sib%egi(i,m3,ipat)     = sib%egi(i,m3-1,ipat)
      sib%egs(i,m3,ipat)     = sib%egs(i,m3-1,ipat)
      sib%hc(i,m3,ipat)      = sib%hc(i,m3-1,ipat)
      sib%hg(i,m3,ipat)      = sib%hg(i,m3-1,ipat)
      sib%ra(i,m3,ipat)      = sib%ra(i,m3-1,ipat)
      sib%rb(i,m3,ipat)      = sib%rb(i,m3-1,ipat)
      sib%rc(i,m3,ipat)      = sib%rc(i,m3-1,ipat)
      sib%rd(i,m3,ipat)      = sib%rd(i,m3-1,ipat)
      sib%roff(i,m3,ipat)    = sib%roff(i,m3-1,ipat)
      sib%green(i,m3,ipat)   = sib%green(i,m3-1,ipat)
      sib%apar(i,m3,ipat)    = sib%apar(i,m3-1,ipat)
      sib%ventmf(i,m3,ipat)  = sib%ventmf(i,m3-1,ipat)
      sib%pco2c(i,m3,ipat)   = sib%pco2c(i,m3-1,ipat)
      sib%pco2i(i,m3,ipat)   = sib%pco2i(i,m3-1,ipat)
      sib%pco2s(i,m3,ipat)   = sib%pco2s(i,m3-1,ipat)
      sib%pco2m(i,m3,ipat)   = sib%pco2m(i,m3-1,ipat)
      sib%ea(i,m3,ipat)      = sib%ea(i,m3-1,ipat)
      sib%em(i,m3,ipat)      = sib%em(i,m3-1,ipat)
      sib%rha(i,m3,ipat)     = sib%rha(i,m3-1,ipat)
      sib%radvbc(i,m3,ipat)  = sib%radvbc(i,m3-1,ipat)
      sib%radvdc(i,m3,ipat)  = sib%radvdc(i,m3-1,ipat)
      sib%radnbc(i,m3,ipat)  = sib%radnbc(i,m3-1,ipat)
      sib%radndc(i,m3,ipat)  = sib%radndc(i,m3-1,ipat)
      sib%psy(i,m3,ipat)     = sib%psy(i,m3-1,ipat)
   enddo   
 enddo
endif

return
END SUBROUTINE sib_bcond

