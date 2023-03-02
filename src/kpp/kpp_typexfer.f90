!##############################################################################
Subroutine mckpp_fields_3dto1d (kpp_fields_3d,kpp_fields_1d,i,j)

use mem_kpp, only:kpp_3d_type,kpp_1d_type
use kpp_parameters, only:kppNZ,kppNZP1,kpprnt,iprnt,jprnt
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid

implicit none

!Accepts a 3D variable of the KPP derived type.
!Returns a 1D variable of the KPP derived type, extracted from the 3D variable
!at a specified i,j.
!Variables are labeled as:
!(H) = needs to be in analysis file for restart (never re-initialize)
!(init) = static initializations from RAMS/KPP
!(update) = fields that are updated over time, such as SST

 TYPE(kpp_3d_type),intent(in) :: kpp_fields_3d
 TYPE(kpp_1d_type),intent(out) :: kpp_fields_1d
 INTEGER :: f,k,i,j,m,p

 !2D fields set by RAMS
 kpp_fields_1d%Sref        = kpp_fields_3d%Sref(i,j)         !Set At Initial (H)
 kpp_fields_1d%jerlov      = nint(kpp_fields_3d%jerlov(i,j)) !Set by RAMS (init)
 kpp_fields_1d%ocdepth     = kpp_fields_3d%ocdepth(i,j)      !Set by RAMS (init)
 kpp_fields_1d%f           = kpp_fields_3d%f(i,j)            !Set by RAMS (init)
 kpp_fields_1d%bottomt     = kpp_fields_3d%bottomt(i,j)      !Set by RAMS (update)
 kpp_fields_1d%SST0        = kpp_fields_3d%SST0(i,j)         !Set by RAMS (update)
 kpp_fields_1d%old         = nint(kpp_fields_3d%old(i,j))    !prognostic (H)
 kpp_fields_1d%new         = nint(kpp_fields_3d%new(i,j))    !prognostic (H)
 !Insert surface level fluxes to ocean model variables
 kpp_fields_1d%sflux(1)    = kpp_fields_3d%flx_ust(i,j)   !Set by RAMS (update)
 kpp_fields_1d%sflux(2)    = kpp_fields_3d%flx_vst(i,j)   !Set by RAMS (update)
 kpp_fields_1d%sflux(3)    = kpp_fields_3d%flx_nsw(i,j)   !Set by RAMS (update)
 kpp_fields_1d%sflux(4)    = kpp_fields_3d%flx_nlw(i,j)   !Set by RAMS (update)
 kpp_fields_1d%sflux(5)    = kpp_fields_3d%flx_ice(i,j)   !Set by RAMS (update)
 kpp_fields_1d%sflux(6)    = kpp_fields_3d%flx_pcp(i,j)   !Set by RAMS (update)

 !Fluxes across DM interfaces
 DO k=0,kppNZ
  kpp_fields_1d%swdk_opt(k) = kpp_fields_3d%swdk_opt(k+1,i,j)  !Set by RAMS (init)
 ENDDO

 !Fields set at ZM grid center levels
 DO k=1,kppNZP1
  kpp_fields_1d%U_init(k,1)  = kpp_fields_3d%U_init(k,i,j)    !Set At Initial (H)
  kpp_fields_1d%U_init(k,2)  = kpp_fields_3d%V_init(k,i,j)    !Set At Initial (H)
  kpp_fields_1d%swfrac(k)    = kpp_fields_3d%swfrac(k,i,j)    !Set by RAMS (init)
  kpp_fields_1d%sal_clim(k)  = kpp_fields_3d%sal_clim(k,i,j)  !Set by RAMS (update)
  kpp_fields_1d%ocnT_clim(k) = kpp_fields_3d%ocnT_clim(k,i,j) !Set by RAMS (update)
  kpp_fields_1d%U(k,1)       = kpp_fields_3d%U(k,i,j)         !prognostic (H)
  kpp_fields_1d%U(k,2)       = kpp_fields_3d%V(k,i,j)         !prognostic (H)
  kpp_fields_1d%Us(k,1,0)    = kpp_fields_3d%Us0(k,i,j)       !prognostic (H)
  kpp_fields_1d%Us(k,1,1)    = kpp_fields_3d%Us1(k,i,j)       !prognostic (H)
  kpp_fields_1d%Us(k,2,0)    = kpp_fields_3d%Vs0(k,i,j)       !prognostic (H)
  kpp_fields_1d%Us(k,2,1)    = kpp_fields_3d%Vs1(k,i,j)       !prognostic (H)
  kpp_fields_1d%X(k,1)       = kpp_fields_3d%X_t(k,i,j)       !prognostic (H)
  kpp_fields_1d%X(k,2)       = kpp_fields_3d%X_s(k,i,j)       !prognostic (H)
  kpp_fields_1d%Xs(k,1,0)    = kpp_fields_3d%Xs_t0(k,i,j)     !prognostic (H)
  kpp_fields_1d%Xs(k,1,1)    = kpp_fields_3d%Xs_t1(k,i,j)     !prognostic (H)
  kpp_fields_1d%Xs(k,2,0)    = kpp_fields_3d%Xs_s0(k,i,j)     !prognostic (H)
  kpp_fields_1d%Xs(k,2,1)    = kpp_fields_3d%Xs_s1(k,i,j)     !prognostic (H)
 ENDDO

 !Print some things for debugging (Saleeby2018)
 if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
  print*,'SFLUX-1-3',kpp_fields_1d%sflux(1),kpp_fields_1d%sflux(2) &
    ,kpp_fields_1d%sflux(3)
  print*,'SFLUX-4-6',kpp_fields_1d%sflux(4),kpp_fields_1d%sflux(5) &
    ,kpp_fields_1d%sflux(6)
 endif

return
END SUBROUTINE mckpp_fields_3dto1d

!##############################################################################
Subroutine mckpp_fields_1dto3d (kpp_fields_1d,kpp_fields_3d,i,j)

use mem_kpp, only: kpp_1d_type,kpp_3d_type
use kpp_parameters, only:kppNZ,kppNZP1,kpprnt,iprnt,jprnt,ikpp
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid,time

implicit none

 TYPE(kpp_3d_type),intent(inout)  :: kpp_fields_3d
 TYPE(kpp_1d_type),intent(in) :: kpp_fields_1d
 INTEGER :: f,k,i,j,m,p

 kpp_fields_3d%old(i,j)        = nint(kpp_fields_1d%old)       !prognostic (H)
 kpp_fields_3d%new(i,j)        = nint(kpp_fields_1d%new)       !prognostic (H)
 kpp_fields_3d%hmix(i,j)       = kpp_fields_1d%hmix            !DIAGNOSTIC
 kpp_fields_3d%freez_flag(i,j) = kpp_fields_1d%freez_flag      !DIAGNOSTIC
 kpp_fields_3d%reset_flag(i,j) = kpp_fields_1d%reset_flag      !DIAGNOSTIC

 !Print some things for debugging (Saleeby2018)
 if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
  print*,'HMIX',time,kpp_fields_1d%hmix,kpp_fields_1d%kmix,kpp_fields_1d%Tref
 endif

 DO k=0,kppNZ
  if(IKPP==2)then !Fluxes
   kpp_fields_3d%wB(k+1,i,j)       = kpp_fields_1d%wB(k)        !DIAGNOSTIC
   kpp_fields_3d%wXNTt(k+1,i,j)    = kpp_fields_1d%wXNT(k)      !DIAGNOSTIC
   kpp_fields_3d%wU(k+1,i,j)       = kpp_fields_1d%wU(k,1)      !DIAGNOSTIC
   kpp_fields_3d%wV(k+1,i,j)       = kpp_fields_1d%wU(k,2)      !DIAGNOSTIC
   kpp_fields_3d%wXt(k+1,i,j)      = kpp_fields_1d%wX(k,1)      !DIAGNOSTIC
   kpp_fields_3d%wXs(k+1,i,j)      = kpp_fields_1d%wX(k,2)      !DIAGNOSTIC
  endif
  !Print some things for debugging (Saleeby2018)
  if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1.and.k<=5) then
   print*,'WX',kpp_fields_1d%wX(k,1),kpp_fields_1d%wX(k,2),kpp_fields_1d%wB(k)
  endif
 ENDDO

 DO k=1,kppNZP1
  if(IKPP==2)then
   kpp_fields_3d%tinc_fcorr(k,i,j) = kpp_fields_1d%tinc_fcorr(k) !DIAGNOSTIC
   kpp_fields_3d%sinc_fcorr(k,i,j) = kpp_fields_1d%sinc_fcorr(k) !DIAGNOSTIC
   kpp_fields_3d%rho(k,i,j)        = kpp_fields_1d%rho(k)        !DIAGNOSTIC
   kpp_fields_3d%cp(k,i,j)         = kpp_fields_1d%cp(k)         !DIAGNOSTIC
  endif
  kpp_fields_3d%buoy(k,i,j)       = kpp_fields_1d%buoy(k)       !DIAGNOSTIC
  kpp_fields_3d%U(k,i,j)          = kpp_fields_1d%U(k,1)        !prognostic (H)
  kpp_fields_3d%V(k,i,j)          = kpp_fields_1d%U(k,2)        !prognostic (H)
  kpp_fields_3d%Us0(k,i,j)        = kpp_fields_1d%Us(k,1,0)     !prognostic (H)
  kpp_fields_3d%Us1(k,i,j)        = kpp_fields_1d%Us(k,1,1)     !prognostic (H)
  kpp_fields_3d%Vs0(k,i,j)        = kpp_fields_1d%Us(k,2,0)     !prognostic (H)
  kpp_fields_3d%Vs1(k,i,j)        = kpp_fields_1d%Us(k,2,1)     !prognostic (H)
  kpp_fields_3d%X_t(k,i,j)        = kpp_fields_1d%X(k,1)        !prognostic (H)
  kpp_fields_3d%X_s(k,i,j)        = kpp_fields_1d%X(k,2)        !prognostic (H)
  kpp_fields_3d%Xs_t0(k,i,j)      = kpp_fields_1d%Xs(k,1,0)     !prognostic (H)
  kpp_fields_3d%Xs_t1(k,i,j)      = kpp_fields_1d%Xs(k,1,1)     !prognostic (H)
  kpp_fields_3d%Xs_s0(k,i,j)      = kpp_fields_1d%Xs(k,2,0)     !prognostic (H)
  kpp_fields_3d%Xs_s1(k,i,j)      = kpp_fields_1d%Xs(k,2,1)     !prognostic (H)
  !Print some things for debugging (Saleeby2018)
  if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1.and.k<=5) then
   print*,'oceanT',kpp_fields_1d%X(k,1),kpp_fields_1d%buoy(k),kpp_fields_1d%Rig(k)
  endif
 ENDDO

 DO k=1,kppNZP1,2
  if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==2.and.k<=20) then
   print*,'oceanT',k,kpp_fields_1d%X(k,1),kpp_fields_1d%X(k,2)+kpp_fields_1d%Sref
  endif
 ENDDO
return
END SUBROUTINE mckpp_fields_1dto3d

