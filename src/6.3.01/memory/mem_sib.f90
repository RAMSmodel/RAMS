!##############################################################################
Module mem_sib

use grid_dims, only:nzpmax

implicit none

! Here are the arrays for using SiB-2 in RAMS

  real, dimension(nzpmax) :: co2_init ! from RAMSIN

  integer, parameter:: nstyp_sib = 12, nvtyp_sib = 13, nzg_sib = 7

  real, dimension(nstyp_sib) :: bee_sib,phsat_sib,satco_sib       &
       ,poros_sib,slope_sib,wopt_sib,skew_sib,respsat_sib

  real, dimension(nvtyp_sib,2,2) :: tran_sib, ref_sib
  real, dimension(nvtyp_sib,2)   :: soref_sib

  real, dimension(nvtyp_sib) :: z2_sib,z1_sib,fvcover_sib,chil_sib  &
       ,sodep_sib,rootd_sib,phc_sib,vmax0_sib,effcon_sib    &
       ,gslope_sib,gsmin_sib,atheta_sib,btheta_sib,trda_sib &
       ,trdm_sib,trop_sib,respcp_sib,slti_sib,shti_sib      &
       ,hlti_sib,hhti_sib

  real,dimension(50) :: laig_sib,fvcg_sib

  real,dimension(nvtyp_sib,50,50) :: a_zo_sib,a_zp_sib,a_rbc_sib    &
       ,a_rdc_sib

  real,dimension(nvtyp_sib) :: zc_sib,zlw_sib,zlen_sib,ltmax_sib    &
       ,stem_sib,nd98_sib,nd02_sib,srmax_sib,srmin_sib

  Type sib_vars

     ! Variables to be dimensioned by (nxp,nyp,npatch)
     real, allocatable, dimension(:,:,:) ::              &
      !Prognostic: required for history restart
       snow1,snow2,capac1,capac2,pco2ap,co2flx           &
      ,sfcswa,uplwrf                                     &
      !Diagnostic: not required
      ,assimn,respg,rstfac1,rstfac2,rstfac3              &
      ,ect,eci,egi,egs,hc,hg                             &
      ,ra,rb,rc,rd,roff,green,apar                       &
      ,ventmf,pco2c,pco2i,pco2s,pco2m,ea,em,rha,radvbc   &
      ,radvdc,radnbc,radndc,psy

     ! Variables to be dimensioned by (nzp,nyp,nxp)
     real, allocatable, dimension(:,:,:) :: rco2p  !mixing ratio of CO2

  End Type

  type (sib_vars), allocatable :: sib_g(:), sibm_g(:)

CONTAINS

!##############################################################################
Subroutine alloc_sib (sib,nz,nx,ny,np)

use mem_leaf, only:isfcl

implicit none

   type (sib_vars) :: sib
   integer, intent(in) :: nz,nx,ny,np

! Allocate arrays based on options (if necessary)
 if (isfcl == 2) then
   !3D atmospheric
   allocate (sib%rco2p     (nz,nx,ny))

   !2D + patches prognostic
   allocate (sib%snow1     (nx,ny,np))
   allocate (sib%snow2     (nx,ny,np))
   allocate (sib%capac1    (nx,ny,np))
   allocate (sib%capac2    (nx,ny,np))
   allocate (sib%pco2ap    (nx,ny,np))
   allocate (sib%co2flx    (nx,ny,np))
   allocate (sib%sfcswa    (nx,ny,np))
   allocate (sib%uplwrf    (nx,ny,np))

   !2D + patches diagnostic
   allocate (sib%assimn    (nx,ny,np))
   allocate (sib%respg     (nx,ny,np))
   allocate (sib%rstfac1   (nx,ny,np))
   allocate (sib%rstfac2   (nx,ny,np))
   allocate (sib%rstfac3   (nx,ny,np))
   allocate (sib%ect       (nx,ny,np))
   allocate (sib%eci       (nx,ny,np))
   allocate (sib%egi       (nx,ny,np))
   allocate (sib%egs       (nx,ny,np))
   allocate (sib%hc        (nx,ny,np))
   allocate (sib%hg        (nx,ny,np))
   allocate (sib%ra        (nx,ny,np))
   allocate (sib%rb        (nx,ny,np))
   allocate (sib%rc        (nx,ny,np))
   allocate (sib%rd        (nx,ny,np))
   allocate (sib%roff      (nx,ny,np))
   allocate (sib%green     (nx,ny,np))
   allocate (sib%apar      (nx,ny,np))
   allocate (sib%ventmf    (nx,ny,np))
   allocate (sib%pco2c     (nx,ny,np))
   allocate (sib%pco2i     (nx,ny,np))
   allocate (sib%pco2s     (nx,ny,np))
   allocate (sib%pco2m     (nx,ny,np))
   allocate (sib%ea        (nx,ny,np))
   allocate (sib%em        (nx,ny,np))
   allocate (sib%rha       (nx,ny,np))
   allocate (sib%radvbc    (nx,ny,np))
   allocate (sib%radvdc    (nx,ny,np))
   allocate (sib%radnbc    (nx,ny,np))
   allocate (sib%radndc    (nx,ny,np))
   allocate (sib%psy       (nx,ny,np))
 endif

return
END SUBROUTINE alloc_sib

!##############################################################################
Subroutine dealloc_sib (sib)

implicit none

  type (sib_vars) :: sib
  !3D atmospheric
  if(allocated(sib%rco2p))     deallocate (sib%rco2p)
  !2D + patches prognostic
  if(allocated(sib%snow1))     deallocate (sib%snow1)
  if(allocated(sib%snow2))     deallocate (sib%snow2)
  if(allocated(sib%capac1))    deallocate (sib%capac1)
  if(allocated(sib%capac2))    deallocate (sib%capac2)
  if(allocated(sib%pco2ap))    deallocate (sib%pco2ap)
  if(allocated(sib%co2flx))    deallocate (sib%co2flx)
  if(allocated(sib%sfcswa))    deallocate (sib%sfcswa)
  if(allocated(sib%uplwrf))    deallocate (sib%uplwrf)
  !2D + patches diagnostic
  if(allocated(sib%assimn))    deallocate (sib%assimn)
  if(allocated(sib%respg))     deallocate (sib%respg)
  if(allocated(sib%rstfac1))   deallocate (sib%rstfac1)
  if(allocated(sib%rstfac2))   deallocate (sib%rstfac2)
  if(allocated(sib%rstfac3))   deallocate (sib%rstfac3)
  if(allocated(sib%ect))       deallocate (sib%ect)
  if(allocated(sib%eci))       deallocate (sib%eci)
  if(allocated(sib%egi))       deallocate (sib%egi)
  if(allocated(sib%egs))       deallocate (sib%egs)
  if(allocated(sib%hc))        deallocate (sib%hc)
  if(allocated(sib%hg))        deallocate (sib%hg)
  if(allocated(sib%ra))        deallocate (sib%ra)
  if(allocated(sib%rb))        deallocate (sib%rb)
  if(allocated(sib%rc))        deallocate (sib%rc)
  if(allocated(sib%rd))        deallocate (sib%rd)
  if(allocated(sib%roff))      deallocate (sib%roff)
  if(allocated(sib%green))     deallocate (sib%green)
  if(allocated(sib%apar))      deallocate (sib%apar)
  if(allocated(sib%ventmf))    deallocate (sib%ventmf)
  if(allocated(sib%pco2c))     deallocate (sib%pco2c)
  if(allocated(sib%pco2i))     deallocate (sib%pco2i)
  if(allocated(sib%pco2s))     deallocate (sib%pco2s)
  if(allocated(sib%pco2m))     deallocate (sib%pco2m)
  if(allocated(sib%ea))        deallocate (sib%ea)
  if(allocated(sib%em))        deallocate (sib%em)
  if(allocated(sib%rha))       deallocate (sib%rha)
  if(allocated(sib%radvbc))    deallocate (sib%radvbc)
  if(allocated(sib%radvdc))    deallocate (sib%radvdc)
  if(allocated(sib%radnbc))    deallocate (sib%radnbc)
  if(allocated(sib%radndc))    deallocate (sib%radndc)
  if(allocated(sib%psy))       deallocate (sib%psy)

return
END SUBROUTINE dealloc_sib

!##############################################################################
Subroutine filltab_sib (sib,sibm,imean,nz,nx,ny,np,ng)

use var_tables

implicit none

   type (sib_vars) :: sib,sibm
   integer, intent(in) :: imean,nz,nx,ny,np,ng
   integer :: npts

! Fill arrays into variable tables

   npts=nz*nx*ny

!SIB SPECIFIC PROGNOSTIC 3D VARIABLES
   if (allocated(sib%rco2p))   &
      CALL vtables2 (sib%rco2p(1,1,1),sibm%rco2p(1,1,1)  &
                 ,ng, npts, imean,  &
                 'RCO2P :3:anal:mpti:mpt1')

   npts=nx*ny*np

!SIB SPECIFIC PROGNOSTIC 2D PATCH VARIABLES
   if (allocated(sib%snow1))   &
      CALL vtables2 (sib%snow1(1,1,1),sibm%snow1(1,1,1) &
                 ,ng, npts, imean,  &
                 'SNOW1 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%snow2))   &
      CALL vtables2 (sib%snow2(1,1,1),sibm%snow2(1,1,1) &
                 ,ng, npts, imean,  &
                 'SNOW2 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%capac1))   &
      CALL vtables2 (sib%capac1(1,1,1),sibm%capac1(1,1,1) &
                 ,ng, npts, imean,  &
                 'CAPAC1 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%capac2))   &
      CALL vtables2 (sib%capac2(1,1,1),sibm%capac2(1,1,1) &
                 ,ng, npts, imean,  &
                 'CAPAC2 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%pco2ap))   &
      CALL vtables2 (sib%pco2ap(1,1,1),sibm%pco2ap(1,1,1) &
                 ,ng, npts, imean,  &
                 'PCO2AP :6:anal:mpti:recycle_sfc')
   if (allocated(sib%co2flx))   &
      CALL vtables2 (sib%co2flx(1,1,1),sibm%co2flx(1,1,1) &
                 ,ng, npts, imean,  &
                 'CO2FLX :6:anal:mpti:recycle_sfc')
   if (allocated(sib%sfcswa))   &
      CALL vtables2 (sib%sfcswa(1,1,1),sibm%sfcswa(1,1,1) &
                 ,ng, npts, imean,  &
                 'SFCSWA :6:anal:mpti:recycle_sfc')
   if (allocated(sib%uplwrf))   &
      CALL vtables2 (sib%uplwrf(1,1,1),sibm%uplwrf(1,1,1) &
                 ,ng, npts, imean,  &
                 'UPLWRF :6:anal:mpti:recycle_sfc')

!SIB SPECIFIC DIAGNOSTIC 2D PATCH VARIABLES
   if (allocated(sib%assimn))   &
      CALL vtables2 (sib%assimn(1,1,1),sibm%assimn(1,1,1) &
                 ,ng, npts, imean,  &
                 'ASSIMN :6:anal:mpti:recycle_sfc')
   if (allocated(sib%respg))   &
      CALL vtables2 (sib%respg(1,1,1),sibm%respg(1,1,1) &
                 ,ng, npts, imean,  &
                 'RESPG :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rstfac1))   &
      CALL vtables2 (sib%rstfac1(1,1,1),sibm%rstfac1(1,1,1) &
                 ,ng, npts, imean,  &
                 'RSTFAC1 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rstfac2))   &
      CALL vtables2 (sib%rstfac2(1,1,1),sibm%rstfac2(1,1,1) &
                 ,ng, npts, imean,  &
                 'RSTFAC2 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rstfac3))   &
      CALL vtables2 (sib%rstfac3(1,1,1),sibm%rstfac3(1,1,1) &
                 ,ng, npts, imean,  &
                 'RSTFAC3 :6:anal:mpti:recycle_sfc')
   if (allocated(sib%ect))   &
      CALL vtables2 (sib%ect(1,1,1),sibm%ect(1,1,1) &
                 ,ng, npts, imean,  &
                 'ECT :6:anal:mpti:recycle_sfc')
   if (allocated(sib%eci))   &
      CALL vtables2 (sib%eci(1,1,1),sibm%eci(1,1,1) &
                 ,ng, npts, imean,  &
                 'ECI :6:anal:mpti:recycle_sfc')
   if (allocated(sib%egi))   &
      CALL vtables2 (sib%egi(1,1,1),sibm%egi(1,1,1) &
                 ,ng, npts, imean,  &
                 'EGI :6:anal:mpti:recycle_sfc')
   if (allocated(sib%egs))   &
      CALL vtables2 (sib%egs(1,1,1),sibm%egs(1,1,1) &
                 ,ng, npts, imean,  &
                 'EGS :6:anal:mpti:recycle_sfc')
   if (allocated(sib%hc))   &
      CALL vtables2 (sib%hc(1,1,1),sibm%hc(1,1,1) &
                 ,ng, npts, imean,  &
                 'HC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%hg))   &
      CALL vtables2 (sib%hg(1,1,1),sibm%hg(1,1,1) &
                 ,ng, npts, imean,  &
                 'HG :6:anal:mpti:recycle_sfc')
   if (allocated(sib%ra))   &
      CALL vtables2 (sib%ra(1,1,1),sibm%ra(1,1,1) &
                 ,ng, npts, imean,  &
                 'RA :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rb))   &
      CALL vtables2 (sib%rb(1,1,1),sibm%rb(1,1,1) &
                 ,ng, npts, imean,  &
                 'RB :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rc))   &
      CALL vtables2 (sib%rc(1,1,1),sibm%rc(1,1,1) &
                 ,ng, npts, imean,  &
                 'RC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rd))   &
      CALL vtables2 (sib%rd(1,1,1),sibm%rd(1,1,1) &
                 ,ng, npts, imean,  &
                 'RD :6:anal:mpti:recycle_sfc')
   if (allocated(sib%roff))   &
      CALL vtables2 (sib%roff(1,1,1),sibm%roff(1,1,1) &
                 ,ng, npts, imean,  &
                 'ROFF :6:anal:mpti:recycle_sfc')
   if (allocated(sib%green))   &
      CALL vtables2 (sib%green(1,1,1),sibm%green(1,1,1) &
                 ,ng, npts, imean,  &
                 'GREEN :6:anal:mpti:recycle_sfc')
   if (allocated(sib%apar))   &
      CALL vtables2 (sib%apar(1,1,1),sibm%apar(1,1,1) &
                 ,ng, npts, imean,  &
                 'APAR :6:anal:mpti:recycle_sfc')
   if (allocated(sib%ventmf))   &
      CALL vtables2 (sib%ventmf(1,1,1),sibm%ventmf(1,1,1) &
                 ,ng, npts, imean,  &
                 'VENTMF :6:anal:mpti:recycle_sfc')
   if (allocated(sib%pco2c))   &
      CALL vtables2 (sib%pco2c(1,1,1),sibm%pco2c(1,1,1) &
                 ,ng, npts, imean,  &
                 'PCO2C :6:anal:mpti:recycle_sfc')
   if (allocated(sib%pco2i))   &
      CALL vtables2 (sib%pco2i(1,1,1),sibm%pco2i(1,1,1) &
                 ,ng, npts, imean,  &
                 'PCO2I :6:anal:mpti:recycle_sfc')
   if (allocated(sib%pco2s))   &
      CALL vtables2 (sib%pco2s(1,1,1),sibm%pco2s(1,1,1) &
                 ,ng, npts, imean,  &
                 'PCO2S :6:anal:mpti:recycle_sfc')
   if (allocated(sib%pco2m))   &
      CALL vtables2 (sib%pco2m(1,1,1),sibm%pco2m(1,1,1) &
                 ,ng, npts, imean,  &
                 'PCO2M :6:anal:mpti:recycle_sfc')
   if (allocated(sib%ea))   &
      CALL vtables2 (sib%ea(1,1,1),sibm%ea(1,1,1) &
                 ,ng, npts, imean,  &
                 'EA :6:anal:mpti:recycle_sfc')
   if (allocated(sib%em))   &
      CALL vtables2 (sib%em(1,1,1),sibm%em(1,1,1) &
                 ,ng, npts, imean,  &
                 'EM :6:anal:mpti:recycle_sfc')
   if (allocated(sib%rha))   &
      CALL vtables2 (sib%rha(1,1,1),sibm%rha(1,1,1) &
                 ,ng, npts, imean,  &
                 'RHA :6:anal:mpti:recycle_sfc')
   if (allocated(sib%radvbc))   &
      CALL vtables2 (sib%radvbc(1,1,1),sibm%radvbc(1,1,1) &
                 ,ng, npts, imean,  &
                 'RADVBC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%radvdc))   &
      CALL vtables2 (sib%radvdc(1,1,1),sibm%radvdc(1,1,1) &
                 ,ng, npts, imean,  &
                 'RADVDC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%radnbc))   &
      CALL vtables2 (sib%radnbc(1,1,1),sibm%radnbc(1,1,1) &
                 ,ng, npts, imean,  &
                 'RADNBC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%radndc))   &
      CALL vtables2 (sib%radndc(1,1,1),sibm%radndc(1,1,1) &
                 ,ng, npts, imean,  &
                 'RADNDC :6:anal:mpti:recycle_sfc')
   if (allocated(sib%psy))   &
      CALL vtables2 (sib%psy(1,1,1),sibm%psy(1,1,1) &
                 ,ng, npts, imean,  &
                 'PSY :6:anal:mpti:recycle_sfc')

return
END SUBROUTINE filltab_sib

END MODULE mem_sib
