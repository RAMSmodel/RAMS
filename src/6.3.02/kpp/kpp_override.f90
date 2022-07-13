!##############################################################################
Subroutine mckpp_overrides_bottomtemp (kpp_1d_fields,kpp_const_fields,i,j)

use mem_kpp, only: kpp_1d_type,kpp_const_type
use kpp_parameters, only: kppNZP1

implicit none

TYPE(kpp_1d_type) :: kpp_1d_fields
TYPE(kpp_const_type) :: kpp_const_fields

integer :: i,j

IF (kpp_const_fields%L_BOTTOM_TEMP) THEN
   kpp_1d_fields%tinc_fcorr(kppNZP1) = &
         kpp_1d_fields%bottomt - kpp_1d_fields%X(kppNZP1,1)
   kpp_1d_fields%X(kppNZP1,1)=kpp_1d_fields%bottomt

   !Compute ocntcorr if doing flux correction method, which RAMS does not do
   !at this time. This is more for a multi-year climate mode. See comments
   !in kpp_ocint.f90 for for info on ocntcorr.
   !kpp_1d_fields%ocnTcorr(kppNZP1)  = &
   !      kpp_1d_fields%tinc_fcorr(kppNZP1) &
   !           * kpp_1d_fields%rho(kppNZP1) &
   !           * kpp_1d_fields%cp(kppNZP1) &
   !           / kpp_const_fields%dto
ENDIF

return
END SUBROUTINE mckpp_overrides_bottomtemp

!##############################################################################
Subroutine mckpp_overrides_profile (kpp_1d_fields,kpp_const_fields,i,j)

use mem_kpp, only: kpp_1d_type,kpp_const_type
use kpp_parameters, only: kppNZP1
use mem_grid, only:ngrid
use node_mod, only:mi0,mj0

implicit none

TYPE(kpp_1d_type) :: kpp_1d_fields
TYPE(kpp_const_type) :: kpp_const_fields

INTEGER :: z,k,j,i
REAL :: dz_total,dtdz_total,dz

  ! If the integration has failed because of unrealistic values in T, S, U or V
  ! or very high RMS difference between the old and new profiles, then reset
  ! T and S to climatology (if available) and U and V to the initial profiles.
  ! NPK 17/5/13.
  IF (kpp_1d_fields%comp_flag) THEN !Saleeby(2019) Make sure we have ocnt/sal profs
   kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
   kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
   kpp_1d_fields%U=kpp_1d_fields%U_init(:,:)
   print*,'Resetting point to climatology(i,j): ',i+mi0(ngrid),j+mj0(ngrid)
   !WRITE(*,*) 'T = ',kpp_1d_fields%ocnT_clim(:)
   !WRITE(*,*) 'S = ',kpp_1d_fields%sal_clim(:)
   !WRITE(*,*) 'U = ',kpp_1d_fields%U_init(:,1)
   !WRITE(*,*) 'V = ',kpp_1d_fields%U_init(:,2)
   kpp_1d_fields%reset_flag=999
  ENDIF
  
  ! Check whether the temperature at any (x,z) point is less than the
  ! threshold for sea ice (-1.8C).  If it is, reset it to -1.8C and
  ! set a flag.
  ! Note that the value of the flag is equal to the *fraction* of levels
  ! at that point that were < -1.8C.
  IF (kpp_const_fields%L_NO_FREEZE) THEN
   DO z=1,kppNZP1
    IF (kpp_1d_fields%X(z,1) .lt. -1.8) THEN
     kpp_1d_fields%tinc_fcorr(z)=kpp_1d_fields%tinc_fcorr(z)+&
          (-1.8-kpp_1d_fields%X(z,1))
     kpp_1d_fields%X(z,1)=-1.8
     kpp_1d_fields%freez_flag=kpp_1d_fields%freez_flag+1.0/REAL(kppNZP1)
     print*,'Setting freezing pt(k,i,j): ',z,i+mi0(ngrid),j+mj0(ngrid)
    ENDIF
   ENDDO
  ENDIF
        
  ! Check whether the temperature difference between the surface
  ! and a user-specified level (presumably deep) is less than a user-specified 
  ! threshold (presumably small).  If so, reset the temperature and salinity
  ! profiles to climatological values.  Added to prevent spurious very
  ! deep mixing that creates unrealistic isothermal (and isohaline) layers.
  ! NPK 15/5/2013 for R4.
  IF (kpp_const_fields%L_NO_ISOTHERM) THEN
     dtdz_total=0.
     dz_total=0.
     DO k=2,kpp_const_fields%iso_bot
        dz=kpp_const_fields%zm(k)-kpp_const_fields%zm(k-1)
        dtdz_total=dtdz_total+ABS((kpp_1d_fields%X(k,1)-&
             kpp_1d_fields%X(k-1,1)))*dz
        dz_total=dz_total+dz
     ENDDO
     dtdz_total=dtdz_total/dz_total
     
     ! If resetting to climo because of isothermal layer (rather than because of 
     ! computational instability trap), then set reset_flag to a negative
     ! value (-1*number of interations in of semi-implicit integration).
     IF (ABS(dtdz_total).lt.kpp_const_fields%iso_thresh) THEN
        kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
        kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
        kpp_1d_fields%reset_flag=(-1.)*kpp_1d_fields%reset_flag
        print*,'Isothermal layer detected(i,j): ',i+mi0(ngrid),j+mj0(ngrid)
     ENDIF
  ELSE
     kpp_1d_fields%reset_flag=0
  ENDIF
  
return
END SUBROUTINE mckpp_overrides_profile
