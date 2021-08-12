!##############################################################################
Subroutine toptnest (ngra,ngrb)

use mem_mksfc
use mem_grid
use io_params
use node_mod

implicit none

integer :: ngra,ngrb 
integer :: ifm,icm,nc1

do ifm = ngra,ngrb
   icm = nxtnest(ifm)
! Initialize TOPOGRAPHY in toptinit.

   CALL toptinit (mmxp(ifm),mmyp(ifm)  &
      ,sfcfile_p(ifm)%topt(1,1),sfcfile_p(ifm)%topzo(1,1))

   if (icm .ge. 1 .and. itoptflg(ifm) .eq. 0) then
      CALL newgrid (ifm)
      CALL interp_fine_grid_sfcvar (ifm,icm,IFG_TOPT1_INIT)

   elseif (itoptflg(ifm) .eq. 1) then

      ! Interpolate TOPO from standard dataset:
      CALL geodat (mmxp(ifm),mmyp(ifm),sfcfile_p(ifm)%topt(1,1)  &
                 ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'TOP')
      ! Interpolate TOPO ZO from standard dataset:
      CALL geodat (mmxp(ifm),mmyp(ifm),sfcfile_p(ifm)%topzo(1,1)  &
               ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'ZOT')

   endif
   
! If desired, override current values of TOPOGRAPHY in ruser.f routine.

   CALL toptinit_user (mmxp(ifm),mmyp(ifm),ifm  &
         ,sfcfile_p(ifm)%topt(1,1) ,sfcfile_p(ifm)%topzo(1,1))

enddo

if (ngra .eq. ngrb) return

! In case topography data have been independently reassigned on any grid,
! average fine mesh topography sequentially to the coarser grids.

do ifm = ngrb,ngra,-1
   if (nxtnest(ifm) .gt. ngridsh .and. ifm .ge. 2) then
      icm = nxtnest(ifm)

      CALL fdback (sfcfile_p(icm)%topt(1,1),sfcfile_p(ifm)%topt(1,1)  &
         ,vt2da,scr2,1,mmxp(icm),mmyp(icm)  &
         ,1,mmxp(ifm),mmyp(ifm),ifm,'terr',vt2db)

   endif
enddo

! In case terrain heights have been independently reassigned on
! any grid, we want to smooth out (re-interpolate) the boundary
! region (roughly one nstratx worth of fine grid points). Do this
! by:
!   1) copy fine topt to scr1
!   2) interpolate topt
!   3) copy interpolated fine topt to scr2
!   4) copy the interior (call to ae2 subroutine) of
!      scr1 (orignial topt) to scr2 (interpolated topt)
!   5) copy scr2 to fine topt

do ifm = ngra,ngrb
   icm = nxtnest(ifm)
   if (icm .ge. 1) then
      ! Copy fine topt to scr1
      CALL ae1 (mmxp(ifm)*mmyp(ifm),scr1,sfcfile_p(ifm)%topt(1,1))

      ! Interpolate topt
      CALL newgrid (ifm)
      CALL interp_fine_grid_sfcvar (ifm,icm,IFG_TOPT2_INIT)

      ! Copy interpolated topt to scr2
      CALL ae1 (mmxp(ifm)*mmyp(ifm),scr2,sfcfile_p(ifm)%topt(1,1))

      ! Copy interior of scr1 to scr2
      nc1 = jdim * (nstraty(ifm) + 1)
      CALL ae2 (mmxp(ifm),mmyp(ifm),2+nstratx(ifm)  &
         ,mmxp(ifm)-1-nstratx(ifm),1+nc1,mmyp(ifm)-nc1  &
         ,scr2,scr1)

      ! Copy scr2 to topt
      CALL ae1 (mmxp(ifm)*mmyp(ifm),sfcfile_p(ifm)%topt(1,1),scr2)

   endif
enddo

return
END SUBROUTINE toptnest

!##############################################################################
Subroutine geonest_file (ifm)

use mem_mksfc
use mem_grid
use io_params
use mem_leaf
use node_mod

implicit none

integer :: ifm 
integer :: icm,ipat,i,j,k

icm = nxtnest(ifm)

! Initialize PATCH AREA, LANDUSE CLASS, and SOIL TEXTURAL CLASS
! in routine sfcinit.

CALL sfcinit_file (mmxp(ifm),mmyp(ifm),nzg,npatch  &
   ,sfcfile_p(ifm)%patch_area(1,1,1)   &
   ,sfcfile_p(ifm)%leaf_class(1,1,1)   &
   ,sfcfile_p(ifm)%soil_text(1,1,1,1))
     

print*, ' '
print*,'====================================================='
print*,'Starting landuse data input on grid ',ifm
print*,'====================================================='

if (icm .ge. 1 .and. ivegtflg(ifm) .eq. 0) then

! Assign PATCH AREAS and PATCH CLASSES from coarser grid:

   do ipat = 1,npatch
      do j = 1,mmyp(ifm)
         do i = 1,mmxp(ifm)

            sfcfile_p(ifm)%patch_area(i,j,ipat) =  &
               sfcfile_p(icm)%patch_area(ipm(i,ifm),jpm(j,ifm),ipat)
            sfcfile_p(ifm)%leaf_class(i,j,ipat) =  &
               sfcfile_p(icm)%leaf_class(ipm(i,ifm),jpm(j,ifm),ipat)
                        
         enddo
      enddo
   enddo
   

elseif (ivegtflg(ifm) .eq. 1) then

! Assign PATCH AREAS and PATCH CLASSES from standard dataset:

   CALL landuse_opqr (mmxp(ifm),mmyp(ifm),nzg,npatch,nvegpat  &
      ,ivegtfn(ifm),isoilfn(ifm) &
      ,ndvifn(ifm),vndvifil(1,ifm)  &
      ,'veg' &
      ,sfcfile_p(ifm)%soil_text(1,1,1,1)  &
      ,sfcfile_p(ifm)%patch_area(1,1,1)   &
      ,sfcfile_p(ifm)%leaf_class(1,1,1)   &
      ,sfcfile_p(ifm)%veg_ndvif(1,1,1))
   

endif

if (icm .ge. 1 .and. isoilflg(ifm) .eq. 0) then

! Assign SOIL TEXTURE CLASS from coarser grid

   do ipat = 2,npatch
      do k = 1,nzg
         do j = 1,mmyp(ifm)
            do i = 1,mmxp(ifm)
               sfcfile_p(ifm)%soil_text(k,i,j,ipat) =  &
               sfcfile_p(icm)%soil_text(k,ipm(i,ifm),jpm(j,ifm),ipat)
            enddo
         enddo
      enddo
   enddo

elseif (isoilflg(ifm) .eq. 1) then

! Assign SOIL TEXTURE CLASS from standard dataset:

   CALL landuse_opqr (mmxp(ifm),mmyp(ifm),nzg,npatch,nvegpat  &
      ,ivegtfn(ifm),isoilfn(ifm) &
      ,ndvifn(ifm),vndvifil(1,ifm)  &
      ,'soil' &
      ,sfcfile_p(ifm)%soil_text(1,1,1,1)  &
      ,sfcfile_p(ifm)%patch_area(1,1,1)   &
      ,sfcfile_p(ifm)%leaf_class(1,1,1)   &
      ,sfcfile_p(ifm)%veg_ndvif(1,1,1))

endif

! If desired, override current values of PATCH AREA, PATCH CLASS, 
! LEAF-2 VEGETATION CLASS, SOIL TEXTURAL CLASS, and/or
! NDVI in ruser.f routines.

CALL sfcinit_file_user (mmxp(ifm),mmyp(ifm),nzg,npatch,ifm &
   ,sfcfile_p(ifm)%patch_area  (1,1,1)    &
   ,sfcfile_p(ifm)%leaf_class (1,1,1)     &
   ,sfcfile_p(ifm)%soil_text   (1,1,1,1) )

! As a final initialization step, eliminate any land patch area that is less 
! than 1% of the total grid cell area.  Set its area to zero, and compensate
! by enlarging areas of remaining patches.

CALL patch_minsize (mmxp(ifm),mmyp(ifm),npatch  &
   ,sfcfile_p(ifm)%patch_area(1,1,1))

return
END SUBROUTINE geonest_file

!##############################################################################
Subroutine geonest_nofile (ngra,ngrb)

use mem_leaf
use mem_sib
use mem_basic
use mem_scratch
use mem_grid
use io_params
use node_mod

implicit none

integer :: ngra,ngrb
integer :: isiz,ifm,icm,ipat,i,j,k,ic,jc

! Initialization/interpolation of leaf-2 variables for which standard RAMS
! datasets never exist.

isiz = maxnxp * maxnyp

do ifm = ngra,ngrb
   icm = nxtnest(ifm)

! First, fill NOFILE LEAF-2 variables with default values in SFCINIT.

   CALL sfcinit_nofile (mmzp(ifm),mmxp(ifm),mmyp(ifm),nzg,nzs    &
      ,npatch,ifm                                               &
      ,basic_g(ifm)%theta          (1,1,1)    &
      ,basic_g(ifm)%pi0            (1,1,1)    &
      ,basic_g(ifm)%pp             (1,1,1)    &
      ,basic_g(ifm)%rv             (1,1,1)    &

      ,leaf_g(ifm)%seatp           (1,1)      &
      ,leaf_g(ifm)%seatf           (1,1)      &

      ,leaf_g(ifm)%soil_water      (1,1,1,1)  & 
      ,leaf_g(ifm)%soil_energy     (1,1,1,1)  &
      ,leaf_g(ifm)%soil_text       (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_mass   (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_energy (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_depth  (1,1,1,1)  &
      ,leaf_g(ifm)%veg_fracarea    (1,1,1)    &
      ,leaf_g(ifm)%veg_lai         (1,1,1)    &
      ,leaf_g(ifm)%veg_tai         (1,1,1)    &
      ,leaf_g(ifm)%veg_rough       (1,1,1)    &
      ,leaf_g(ifm)%veg_height      (1,1,1)    &
      ,leaf_g(ifm)%veg_albedo      (1,1,1)    &
      ,leaf_g(ifm)%patch_area      (1,1,1)    &
      ,leaf_g(ifm)%patch_rough     (1,1,1)    &
      ,leaf_g(ifm)%leaf_class      (1,1,1)    &
      ,leaf_g(ifm)%soil_rough      (1,1,1)    &
      ,leaf_g(ifm)%sfcwater_nlev   (1,1,1)    &
      ,leaf_g(ifm)%stom_resist     (1,1,1)    &
      ,leaf_g(ifm)%ground_rsat     (1,1,1)    &
      ,leaf_g(ifm)%ground_rvap     (1,1,1)    &
      ,leaf_g(ifm)%veg_water       (1,1,1)    &
      ,leaf_g(ifm)%veg_temp        (1,1,1)    &
      ,leaf_g(ifm)%can_rvap        (1,1,1)    &
      ,leaf_g(ifm)%can_temp        (1,1,1)    & 
      ,leaf_g(ifm)%veg_ndvip       (1,1,1)    &
      ,leaf_g(ifm)%veg_ndvic       (1,1,1)    & 
      ,leaf_g(ifm)%veg_ndvif       (1,1,1)    &
      ,leaf_g(ifm)%snow_mass       (1,1)      &
      ,leaf_g(ifm)%snow_depth      (1,1)      &
      ,leaf_g(ifm)%soil_moist_top  (1,1)      &
      ,leaf_g(ifm)%soil_moist_bot  (1,1)      &
      ,leaf_g(ifm)%soil_temp_top   (1,1)      &
      ,leaf_g(ifm)%soil_temp_bot   (1,1)      &
      ,grid_g(ifm)%glat            (1,1)      &
      ,grid_g(ifm)%glon            (1,1)      &
      ,grid_g(ifm)%topzo           (1,1))

! Assignment section for NOFILE surface variables (LEAF3 and/or SiB).
! If needed, this section of code could be used to assign parent
! grid variable to nested grid variables in the following manner.
!   if (icm > 0) then
!    do ipat = 1,npatch
!     do j = 1,nnyp(ifm)
!      do i = 1,nnxp(ifm)
!       ic = ipm(i,ifm)
!       jc = jpm(j,ifm)
!       do k = 1,nzg
!        leaf_g(ifm)%soil_water           (k,i,j,ipat) = &
!             leaf_g(icm)%soil_water      (k,ic,jc,ipat)
!        leaf_g(ifm)%soil_energy          (k,i,j,ipat) = &
!             leaf_g(icm)%soil_energy     (k,ic,jc,ipat)
!       enddo
!       do k = 1,nzs
!        leaf_g(ifm)%sfcwater_mass        (k,i,j,ipat) = &
!             leaf_g(icm)%sfcwater_mass   (k,ic,jc,ipat)
!        leaf_g(ifm)%sfcwater_energy      (k,i,j,ipat) = &
!             leaf_g(icm)%sfcwater_energy (k,ic,jc,ipat)
!        leaf_g(ifm)%sfcwater_depth       (k,i,j,ipat) = &
!             leaf_g(icm)%sfcwater_depth  (k,ic,jc,ipat)
!       enddo
!       leaf_g(ifm)%veg_fracarea         (i,j,ipat) = &
!            leaf_g(icm)%veg_fracarea    (ic,jc,ipat)
!       leaf_g(ifm)%veg_lai              (i,j,ipat) = &
!            leaf_g(icm)%veg_lai         (ic,jc,ipat)
!       leaf_g(ifm)%veg_tai              (i,j,ipat) = &
!            leaf_g(icm)%veg_tai         (ic,jc,ipat)
!       leaf_g(ifm)%veg_rough            (i,j,ipat) = &
!            leaf_g(icm)%veg_rough       (ic,jc,ipat)
!       leaf_g(ifm)%veg_height           (i,j,ipat) = &
!            leaf_g(icm)%veg_height      (ic,jc,ipat)
!       leaf_g(ifm)%veg_albedo           (i,j,ipat) = &
!            leaf_g(icm)%veg_albedo      (ic,jc,ipat)
!       leaf_g(ifm)%patch_rough          (i,j,ipat) = &
!            leaf_g(icm)%patch_rough     (ic,jc,ipat)
!       leaf_g(ifm)%soil_rough           (i,j,ipat) = &
!            leaf_g(icm)%soil_rough      (ic,jc,ipat)
!       leaf_g(ifm)%sfcwater_nlev        (i,j,ipat) = &
!            leaf_g(icm)%sfcwater_nlev   (ic,jc,ipat)
!       leaf_g(ifm)%stom_resist          (i,j,ipat) = &
!            leaf_g(icm)%stom_resist     (ic,jc,ipat)
!       leaf_g(ifm)%ground_rsat          (i,j,ipat) = &
!            leaf_g(icm)%ground_rsat     (ic,jc,ipat)
!       leaf_g(ifm)%ground_rvap          (i,j,ipat) = &
!            leaf_g(icm)%ground_rvap     (ic,jc,ipat)
!       leaf_g(ifm)%veg_water            (i,j,ipat) = &
!            leaf_g(icm)%veg_water       (ic,jc,ipat)
!       leaf_g(ifm)%veg_temp             (i,j,ipat) = &
!            leaf_g(icm)%veg_temp        (ic,jc,ipat)
!       leaf_g(ifm)%can_rvap             (i,j,ipat) = &
!            leaf_g(icm)%can_rvap        (ic,jc,ipat)
!       leaf_g(ifm)%can_temp             (i,j,ipat) = &
!            leaf_g(icm)%can_temp        (ic,jc,ipat)
!       leaf_g(ifm)%veg_ndvic            (i,j,ipat) = &
!            leaf_g(icm)%veg_ndvic       (ic,jc,ipat)
!       !If running SIB land surface scheme
!       if(isfcl==2)then
!        sib_g(ifm)%snow1   (i,j,ipat) = sib_g(icm)%snow1   (ic,jc,ipat)
!        sib_g(ifm)%snow2   (i,j,ipat) = sib_g(icm)%snow2   (ic,jc,ipat)
!        sib_g(ifm)%capac1  (i,j,ipat) = sib_g(icm)%capac1  (ic,jc,ipat)
!        sib_g(ifm)%capac2  (i,j,ipat) = sib_g(icm)%capac2  (ic,jc,ipat)
!        sib_g(ifm)%pco2ap  (i,j,ipat) = sib_g(icm)%pco2ap  (ic,jc,ipat)
!        sib_g(ifm)%co2flx  (i,j,ipat) = sib_g(icm)%co2flx  (ic,jc,ipat)
!        sib_g(ifm)%sfcswa  (i,j,ipat) = sib_g(icm)%sfcswa  (ic,jc,ipat)
!        sib_g(ifm)%uplwrf  (i,j,ipat) = sib_g(icm)%uplwrf  (ic,jc,ipat)
!        sib_g(ifm)%assimn  (i,j,ipat) = sib_g(icm)%assimn  (ic,jc,ipat)
!        sib_g(ifm)%respg   (i,j,ipat) = sib_g(icm)%respg   (ic,jc,ipat)
!        sib_g(ifm)%rstfac1 (i,j,ipat) = sib_g(icm)%rstfac1 (ic,jc,ipat)
!        sib_g(ifm)%rstfac2 (i,j,ipat) = sib_g(icm)%rstfac2 (ic,jc,ipat)
!        sib_g(ifm)%rstfac3 (i,j,ipat) = sib_g(icm)%rstfac3 (ic,jc,ipat)
!        sib_g(ifm)%ect     (i,j,ipat) = sib_g(icm)%ect     (ic,jc,ipat)
!        sib_g(ifm)%eci     (i,j,ipat) = sib_g(icm)%eci     (ic,jc,ipat)
!        sib_g(ifm)%egi     (i,j,ipat) = sib_g(icm)%egi     (ic,jc,ipat)
!        sib_g(ifm)%egs     (i,j,ipat) = sib_g(icm)%egs     (ic,jc,ipat)
!        sib_g(ifm)%hc      (i,j,ipat) = sib_g(icm)%hc      (ic,jc,ipat)
!        sib_g(ifm)%hg      (i,j,ipat) = sib_g(icm)%hg      (ic,jc,ipat)
!        sib_g(ifm)%ra      (i,j,ipat) = sib_g(icm)%ra      (ic,jc,ipat)
!        sib_g(ifm)%rb      (i,j,ipat) = sib_g(icm)%rb      (ic,jc,ipat)
!        sib_g(ifm)%rc      (i,j,ipat) = sib_g(icm)%rc      (ic,jc,ipat)
!        sib_g(ifm)%rd      (i,j,ipat) = sib_g(icm)%rd      (ic,jc,ipat)
!        sib_g(ifm)%roff    (i,j,ipat) = sib_g(icm)%roff    (ic,jc,ipat)
!        sib_g(ifm)%green   (i,j,ipat) = sib_g(icm)%green   (ic,jc,ipat)
!        sib_g(ifm)%apar    (i,j,ipat) = sib_g(icm)%apar    (ic,jc,ipat)
!        sib_g(ifm)%ventmf  (i,j,ipat) = sib_g(icm)%ventmf  (ic,jc,ipat)
!        sib_g(ifm)%pco2c   (i,j,ipat) = sib_g(icm)%pco2c   (ic,jc,ipat)
!        sib_g(ifm)%pco2i   (i,j,ipat) = sib_g(icm)%pco2i   (ic,jc,ipat)
!        sib_g(ifm)%pco2s   (i,j,ipat) = sib_g(icm)%pco2s   (ic,jc,ipat)
!        sib_g(ifm)%pco2m   (i,j,ipat) = sib_g(icm)%pco2m   (ic,jc,ipat)
!        sib_g(ifm)%ea      (i,j,ipat) = sib_g(icm)%ea      (ic,jc,ipat)
!        sib_g(ifm)%em      (i,j,ipat) = sib_g(icm)%em      (ic,jc,ipat)
!        sib_g(ifm)%rha     (i,j,ipat) = sib_g(icm)%rha     (ic,jc,ipat)
!        sib_g(ifm)%radvbc  (i,j,ipat) = sib_g(icm)%radvbc  (ic,jc,ipat)
!        sib_g(ifm)%radvdc  (i,j,ipat) = sib_g(icm)%radvdc  (ic,jc,ipat)
!        sib_g(ifm)%radnbc  (i,j,ipat) = sib_g(icm)%radnbc  (ic,jc,ipat)
!        sib_g(ifm)%radndc  (i,j,ipat) = sib_g(icm)%radndc  (ic,jc,ipat)
!        sib_g(ifm)%psy     (i,j,ipat) = sib_g(icm)%psy     (ic,jc,ipat)
!       endif
!      enddo
!     enddo
!    enddo
!   endif
      
! Override any of the above variable assignments by user-specified changes
! to routine sfcinit_nofile_user.

   CALL sfcinit_nofile_user (mmzp(ifm),mmxp(ifm),mmyp(ifm),nzg,nzs    &
      ,npatch,ifm                                               &
      ,basic_g(ifm)%theta          (1,1,1)    &
      ,basic_g(ifm)%pi0            (1,1,1)    &
      ,basic_g(ifm)%pp             (1,1,1)    &
      ,basic_g(ifm)%rv             (1,1,1)    &

      ,leaf_g(ifm)%seatp           (1,1)      &
      ,leaf_g(ifm)%seatf           (1,1)      &

      ,leaf_g(ifm)%soil_water      (1,1,1,1)  & 
      ,leaf_g(ifm)%soil_energy     (1,1,1,1)  &
      ,leaf_g(ifm)%soil_text       (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_mass   (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_energy (1,1,1,1)  &
      ,leaf_g(ifm)%sfcwater_depth  (1,1,1,1)  &
      ,leaf_g(ifm)%veg_fracarea    (1,1,1)    &
      ,leaf_g(ifm)%veg_lai         (1,1,1)    &
      ,leaf_g(ifm)%veg_tai         (1,1,1)    &
      ,leaf_g(ifm)%veg_rough       (1,1,1)    &
      ,leaf_g(ifm)%veg_height      (1,1,1)    &
      ,leaf_g(ifm)%veg_albedo      (1,1,1)    &
      ,leaf_g(ifm)%patch_area      (1,1,1)    &
      ,leaf_g(ifm)%patch_rough     (1,1,1)    &
      ,leaf_g(ifm)%leaf_class      (1,1,1)    &
      ,leaf_g(ifm)%soil_rough      (1,1,1)    &
      ,leaf_g(ifm)%sfcwater_nlev   (1,1,1)    &
      ,leaf_g(ifm)%stom_resist     (1,1,1)    &
      ,leaf_g(ifm)%ground_rsat     (1,1,1)    &
      ,leaf_g(ifm)%ground_rvap     (1,1,1)    &
      ,leaf_g(ifm)%veg_water       (1,1,1)    &
      ,leaf_g(ifm)%veg_temp        (1,1,1)    &
      ,leaf_g(ifm)%can_rvap        (1,1,1)    &
      ,leaf_g(ifm)%can_temp        (1,1,1)    & 
      ,leaf_g(ifm)%veg_ndvip       (1,1,1)    &
      ,leaf_g(ifm)%veg_ndvic       (1,1,1)    & 
      ,leaf_g(ifm)%veg_ndvif       (1,1,1)    &
      ,leaf_g(ifm)%snow_mass       (1,1)      &
      ,leaf_g(ifm)%snow_depth      (1,1)      &
      ,leaf_g(ifm)%soil_moist_top  (1,1)      &
      ,leaf_g(ifm)%soil_moist_bot  (1,1)      &
      ,leaf_g(ifm)%soil_temp_top   (1,1)      &
      ,leaf_g(ifm)%soil_temp_bot   (1,1)      &
      ,grid_g(ifm)%glat            (1,1)      &
      ,grid_g(ifm)%glon            (1,1)      &
      ,grid_g(ifm)%topzo           (1,1))

enddo

return
END SUBROUTINE geonest_nofile

!##############################################################################
Subroutine patch_minsize (n2,n3,npat,patch_area)

implicit none

integer :: n2,n3,npat,i,j,ipat,jpat
real :: orig_size
real, dimension(n2,n3,npat) :: patch_area

do j = 1,n3
   do i = 1,n2
      do ipat = 2,npat
         if (patch_area(i,j,ipat) .gt. 0. .and.  &
             patch_area(i,j,ipat) .lt. .01) then

            orig_size = patch_area(i,j,ipat)
            patch_area(i,j,ipat) = 0.

            do jpat = 1,npat
               if (jpat .ne. ipat) then
                  patch_area(i,j,jpat) = patch_area(i,j,jpat)  &
                     / (1. - orig_size)
               endif
            enddo
         endif
      enddo
   enddo
enddo

return
END SUBROUTINE patch_minsize

