!##############################################################################
Module isan_coms

use grid_dims

implicit none

!---------------------------------------------------------------------------
!    Configuration RAMS isentropic data analysis package.
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
integer :: ioflgisz,ioflgvar,natime,iszstage,ivrstage,iyear,imonth,idate  &
          ,ihour,isan_inc,i1st_flg,iupa_flg,isfc_flg,idatain
!---------------------------------------------------------------------------
character(len=strl1) :: innpr,inrawi,varpfx,insrfce,iapr,iarawi,iasrfce,var_hfile
character(len=8) :: pdata

!---------------------------------------------------------------------------
!     Input pressure file header
!---------------------------------------------------------------------------
integer :: iyy,imm,idd,ihh,itinc,inproj
real    :: xnelat,xnelon,cntlat,cntlon

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Input pressure data memory

real, allocatable, dimension(:,:,:) :: p_u,p_v,p_t,p_z,p_r,p_ur,p_vr
real, allocatable, dimension(:,:)   :: p_lat,p_lon
real, allocatable, dimension(:,:)   :: p_soilmoist1,p_soilmoist2,p_soiltemp1 &
                                      ,p_soiltemp2,p_snowmass,p_snowdepth
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Polar-stereo/pressure grid memory

real, allocatable, dimension(:,:,:) :: pp_u,pp_v,pp_t,pp_z,pp_r
real, allocatable, dimension(:,:)   :: pp_sglob
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Polar-stereo/isentropic grid memory

real, allocatable, dimension(:,:,:) :: pi_u,pi_v,pi_p,pi_s,pi_r  &
                                      ,pi_scra,pi_scrb
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Polar-stereo/sigma-z grid memory
!                         :: 
real, allocatable, dimension(:,:,:) :: ps_u,ps_v,ps_p,ps_t,ps_r  &
                                      ,ps_scra,ps_scrb
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Polar-stereo/surface grid memory

real, allocatable, dimension(:,:) :: rs_u,rs_v,rs_p,rs_t,rs_r,rs_s  &
                                    ,rs_top,rs_qual  &
                                    ,rs_soilmoist1,rs_soilmoist2,rs_soiltemp1 &
                                    ,rs_soiltemp2,rs_snowmass,rs_snowdepth
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Data type to replace A array memory use in ISAN. 

type isan_grids
   real, allocatable, dimension(:,:,:) :: rr_u,rr_v,rr_t,rr_p,rr_r,rr_cond  
   real, allocatable, dimension(:,:,:) :: rr_ug,rr_vg,rr_tg,rr_pg,rr_rg  
   real, allocatable, dimension(:,:,:) :: rr_pi0,rr_th0,rr_dn0,rr_dn0u,rr_dn0v  
   real, allocatable, dimension(:,:)   :: rr_soilmoist1,rr_soilmoist2,rr_soiltemp1 &
                                     ,rr_soiltemp2,rr_snowmass,rr_snowdepth  
end type
real, allocatable, dimension(:)    :: rr_scr1,rr_scr2,rr_scr3,rr_vt2da

type (isan_grids)                  :: is_grids(maxagrds)

real, dimension(maxsigz,maxagrds)  :: piref,thref,dnref,rtref

integer                            :: maxix,maxiy,maxiz

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Input observation data memory
!
real, allocatable, dimension(:,:)  :: up_uz,up_vz,up_ur,up_vr,up_zz  &
                                     ,up_p,up_t,up_z,up_r
real, allocatable, dimension(:)    :: up_lat,up_lon,up_top
real, allocatable, dimension(:,:)  :: up_topg
integer, allocatable, dimension(:) :: up_lp, up_lz
character(len=8), allocatable, dimension(:) :: up_chstid

real, allocatable, dimension(:)    :: sf_u,sf_v,sf_p,sf_t,sf_s,sf_r
real, allocatable, dimension(:)    :: sf_ur,sf_vr
real, allocatable, dimension(:)    :: sf_lat,sf_lon,sf_top,sf_scra
character(len=8), allocatable, dimension(:) :: sf_chstid
character(len=14), allocatable, dimension(:) :: sf_date
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!  Upper air-isentropic/sigma-z memory
!
real, allocatable, dimension(:,:)  :: upi_u,upi_v,upi_p,upi_s,upi_r
real, allocatable, dimension(:,:)  :: ups_u,ups_v,ups_p,ups_t,ups_r
!---------------------------------------------------------------------------

integer                          :: npdates
integer, dimension(maxisfiles,4) :: iproc_flag
!---------------------------------------------------------------------------
character(len=strl1), dimension(maxisfiles,4) :: iproc_names
!---------------------------------------------------------------------------
character(len=14), dimension(maxisfiles) :: iproc_dates
!---------------------------------------------------------------------------
integer                   :: nprx,npry,nprz,idatelin,iglobew,iglobs,iglobn
integer, dimension(maxpr) :: levpr
real                      :: xswlon,xswlat,gdatdx,gdatdy
real, dimension(maxpr)    :: pnpr
!---------------------------------------------------------------------------
integer                   :: nsta,nssfc,notsta  &
                            ,maxsta,maxsfc,iobswin
real                      :: stasep
real, dimension(maxagrds) :: wvlnth,respon,swvlnth
!---------------------------------------------------------------------------
character(len=8), dimension(50) :: notid
character(len=strl1) :: used_file
!---------------------------------------------------------------------------
integer                    :: nisx,nisy,nisn,interp,igridfl,nigrids,nsigz  &
                             ,nfeedvar
integer, dimension(maxisn) :: levth
real                       :: gobsep,gobrad,topsigz,hybbot,hybtop,sfcinf  &
                             ,sigzwt
real, dimension(maxsigz)   :: sigz
real, dimension(maxagrds)  :: gridwt
!---------------------------------------------------------------------------

END MODULE isan_coms
