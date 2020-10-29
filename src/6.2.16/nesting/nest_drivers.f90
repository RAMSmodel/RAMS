!##############################################################################
Subroutine nstbdriv ()

use mem_tend
use var_tables
use mem_basic
use mem_nestb
use node_mod
use mem_grid

implicit none

integer :: n
real :: tymeinvv,tymeinvs

tymeinvv = 1.0 / (dtlv * 0.5 * float(nndtrat(ngrid)+2-isstp))
tymeinvs = 1.0 / (dtlt * float(nndtrat(ngrid)+1-isstp))

!         print*, 'calling nstbtnd for u'
!         print*, 'jdim,ibcon,ia,iz,ja,jz'
!         print*,  jdim,ibcon,ia,iz,ja,jz


CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%up(1,1,1),tend%ut(1)  &
   ,nbounds(ngrid)%bux(1,1,1),nbounds(ngrid)%buy(1,1,1)  &
   ,nbounds(ngrid)%buz(1,1,1)  &
   ,'u',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%vp(1,1,1),tend%vt(1)  &
   ,nbounds(ngrid)%bvx(1,1,1),nbounds(ngrid)%bvy(1,1,1)  &
   ,nbounds(ngrid)%bvz(1,1,1)  &
   ,'v',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%wp(1,1,1),tend%wt(1)  &
   ,nbounds(ngrid)%bwx(1,1,1),nbounds(ngrid)%bwy(1,1,1)  &
   ,nbounds(ngrid)%bwz(1,1,1)  &
   ,'w',tymeinvv,nstbot,nsttop,jdim)

CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
   ,basic_g(ngrid)%pp(1,1,1),tend%pt(1)  &
   ,nbounds(ngrid)%bpx(1,1,1),nbounds(ngrid)%bpy(1,1,1)  &
   ,nbounds(ngrid)%bpz(1,1,1)  &
   ,'p',tymeinvv,nstbot,nsttop,jdim)

do n = 1,num_scalar(ngrid)
   CALL nstbtnd (mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t  &
      ,nbounds(ngrid)%bsx(1,1,1,n),nbounds(ngrid)%bsy(1,1,1,n)  &
      ,nbounds(ngrid)%bsz(1,1,1,n)  &
      ,'t',tymeinvs,nstbot,nsttop,jdim)
enddo

return
END SUBROUTINE nstbdriv

!##############################################################################
Subroutine interp_fine_grid (ifm,icm,ifg_mode)

! This routine will perform initialization of a fine grid by interpolating
! from the parent (coarse) grid.
!
! The argument ifg_mode denotes which subsystem (consequently which variable
! set) to initialize. The following modes are currently supported, and the
! corresponding parameters are found in node_mode.f90.
!
!    parameter           mode
!
!  IFG_HH_INIT       horizonally homogeneous initialization
!  IFG_VF_INIT       var file initialization
!  IFG_HIST_INIT     history file initialization
!  IFG_MKVF_INIT     initialization for MAKEVFILE
!  IFG_MKHF_INIT     initialization for MAKEHFILE
!  IFG_TOPT1_INIT    initialization for MAKESFC, topography
!  IFG_TOPT2_INIT    initialization for MAKESFC, topography
!  IFG_SST_INIT      initialization for MAKESFC, sea surface temperature
!
! It is the caller's responsibility to check that ifm and icm are a valid
! coarse-fine grid pair. The caller also needs to do a call to newgrid(ifm).

  use node_mod
  use mem_grid
  use mem_basic
  use mem_mksfc
  use ref_sounding
  use var_tables
  use isan_coms

  implicit none

  integer :: ifm, icm
  integer :: ifg_mode

  integer :: nzc, nxc, nyc
  integer :: nzct, nxct, nyct
  integer :: nzf, nxf, nyf
  integer :: ic0, jc0
  integer :: ict0, jct0
  integer :: if0, jf0
  integer :: ibconf
  integer :: ifileok
  integer :: i, j, k, ivar
  real, dimension(:,:,:), allocatable :: dn0ct, dn0uct, dn0vct
  real, dimension(:), allocatable :: dn01dnct
  real, dimension(:,:), allocatable :: toptf_int

  ! These are for convenience
  nxf = mmxp(ifm)
  nyf = mmyp(ifm)
  nzf = mmzp(ifm)

  nxc = mmxp(icm)
  nyc = mmyp(icm)
  nzc = mmzp(icm)

  nxct = ctnxp(ifm)
  nyct = ctnyp(ifm)
  nzct = ctnzp(ifm)

  ic0 = mi0(icm)
  jc0 = mj0(icm)

  ict0 = cti0(ifm)
  jct0 = ctj0(ifm)

  if0 = mi0(ifm)
  jf0 = mj0(ifm)

  ibconf = mibcon(ifm)

  allocate(dn01dnct(nzct))
  allocate(toptf_int(nxf,nyf))
  allocate(dn0ct(nzct,nxct,nyct))
  allocate(dn0uct(nzct,nxct,nyct))
  allocate(dn0vct(nzct,nxct,nyct))

  ! Fill the "tiles" for the density vars
  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT) &
                                  .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL fill_density_tiles (nzc,nxc,nyc,dn01dn(1,icm), &
        basic_g(icm)%dn0,basic_g(icm)%dn0u,basic_g(icm)%dn0v, &
        nzct,nxct,nyct,dn01dnct,dn0ct,dn0uct,dn0vct,ic0,jc0,ict0,jct0,ifm,icm)
  elseif ((ifg_mode .eq. IFG_MKVF_INIT) .or. (ifg_mode .eq. IFG_MKHF_INIT)) then
    CALL fill_density_tiles (nzc,nxc,nyc,dnref(1,icm), &
        is_grids(icm)%rr_dn0,is_grids(icm)%rr_dn0u,is_grids(icm)%rr_dn0v, &
        nzct,nxct,nyct,dn01dnct,dn0ct,dn0uct,dn0vct,ic0,jc0,ict0,jct0,ifm,icm)
  endif


  ! Interpolate topography for subsequent interp_var calls, or for MAKESFC
  if ((ifg_mode .eq. IFG_TOPT1_INIT) .or. (ifg_mode .eq. IFG_TOPT2_INIT)) then
    CALL interp_var (1,nxc,nyc,sfcfile_p(icm)%topt,1,nxct,nyct,dn0ct, &
                    1,nxf,nyf,sfcfile_p(ifm)%topt,basic_g(ifm)%dn0,  &
                    grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm), &
                    ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                    ibconf,'t',0,0,2)
    if (ifg_mode .eq. IFG_TOPT1_INIT) then
      CALL interp_var (1,nxc,nyc,sfcfile_p(icm)%topzo,1,nxct,nyct,dn0ct, &
                      1,nxf,nyf,sfcfile_p(ifm)%topzo,basic_g(ifm)%dn0,  &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),  &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                      ibconf,'t',0,0,2)
    endif
  else
    CALL interp_var (1,nxc,nyc,grid_g(icm)%topt,1,nxct,nyct,dn0ct, &
                    1,nxf,nyf,toptf_int,basic_g(ifm)%dn0,grid_g(ifm)%topt, &
                    toptf_int,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm), &
                    ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',0,0,2)
  endif

  ! Interpolate SST for MAKESFC
  if (ifg_mode .eq. IFG_SST_INIT) then
    CALL interp_var (1,nxc,nyc,sfcfile_p(icm)%seatf,1,nxct,nyct,dn0ct, &
                    1,nxf,nyf,sfcfile_p(ifm)%seatf,basic_g(ifm)%dn0, &
                    grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm), &
                    ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                    ibconf,'t',0,0,2)
  endif


  ! Interpolate the 1D reference state
  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT) &
                                  .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL interp_refs_1d (nzc,dn01dn(1,icm),th01dn(1,icm),u01dn(1,icm),          &
           v01dn(1,icm),rt01dn(1,icm),pi01dn(1,icm),nzct,dn01dnct,              &
           nzf,nxf,nyf,dn01dn(1,ifm),th01dn(1,ifm),u01dn(1,ifm),v01dn(1,ifm),   &
           rt01dn(1,ifm),pi01dn(1,ifm),grid_g(ifm)%topt,toptf_int,              &
           ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ibconf,ifg_mode)
  elseif (ifg_mode .eq. IFG_MKVF_INIT) then
    CALL interp_refs_1d (nzc,dnref(1,icm),thref(1,icm),u01dn(1,icm),            &
           v01dn(1,icm),rtref(1,icm),piref(1,icm),nzct,dn01dnct,                &
           nzf,nxf,nyf,dnref(1,ifm),thref(1,ifm),u01dn(1,ifm),v01dn(1,ifm),     &
           rtref(1,ifm),piref(1,ifm),grid_g(ifm)%topt,toptf_int,                &
           ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ibconf,ifg_mode)
  endif


  ! Interpolate the 3D reference state
  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT) &
                                  .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL interp_refs_3d (nzc,nxc,nyc,basic_g(icm)%dn0,basic_g(icm)%th0,nzct,     &
          nxct,nyct,dn0ct,nzf,nxf,nyf,basic_g(ifm)%dn0,basic_g(ifm)%dn0u,       &
          basic_g(ifm)%dn0v,basic_g(ifm)%th0,basic_g(ifm)%pi0,basic_g(ifm)%pp,  &
          grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),         &
          ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,ifg_mode)
  elseif ((ifg_mode .eq. IFG_MKVF_INIT) .or. (ifg_mode .eq. IFG_MKHF_INIT)) then
    CALL interp_refs_3d (nzc,nxc,nyc,is_grids(icm)%rr_dn0,is_grids(icm)%rr_th0,  &
          nzct,nxct,nyct,dn0ct,nzf,nxf,nyf,is_grids(ifm)%rr_dn0,                &
          is_grids(ifm)%rr_dn0u,is_grids(ifm)%rr_dn0v,is_grids(ifm)%rr_th0,     &
          is_grids(ifm)%rr_pi0,is_grids(ifm)%rr_p,grid_g(ifm)%topt,toptf_int,   &
          ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0, &
          jf0,ibconf,ifg_mode)
  endif


  ! Interpolate the prognostic variables (u, v, w, p, and scalars), then do the
  ! terrain adjustment.
  !   Skip this step if doing history file initialization, as well as any file
  ! creation processes (MAKESFC, etc)
  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL interp_var (nzc,nxc,nyc,basic_g(icm)%uc,nzct,nxct,nyct,dn0uct, &
      nzf,nxf,nyf,basic_g(ifm)%uc,basic_g(ifm)%dn0u,grid_g(ifm)%topt,toptf_int, &
      ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
      ibconf,'u',1,1,3)

    CALL interp_var (nzc,nxc,nyc,basic_g(icm)%vc,nzct,nxct,nyct,dn0vct, &
      nzf,nxf,nyf,basic_g(ifm)%vc,basic_g(ifm)%dn0v,grid_g(ifm)%topt,toptf_int, &
      ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
      ibconf,'v',1,1,3)

    CALL interp_var (nzc,nxc,nyc,basic_g(icm)%wc,nzct,nxct,nyct,dn0ct, &
      nzf,nxf,nyf,basic_g(ifm)%wc,basic_g(ifm)%dn0,grid_g(ifm)%topt,toptf_int, &
      ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
      ibconf,'w',1,1,3)

    CALL interp_var (nzc,nxc,nyc,basic_g(icm)%pc,nzct,nxct,nyct,dn0ct, &
      nzf,nxf,nyf,basic_g(ifm)%pc,basic_g(ifm)%dn0,grid_g(ifm)%topt,toptf_int, &
      ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
      ibconf,'t',0,1,3)

    do ivar = 1,num_scalar(ifm)
      CALL interp_var (nzc,nxc,nyc,scalar_tab(ivar,icm)%var_p,nzct,nxct,nyct, &
        dn0ct,nzf,nxf,nyf,scalar_tab(ivar,ifm)%var_p,basic_g(ifm)%dn0, &
        grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm), &
        ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',1,1,3)
    enddo
  endif

  ! If initializing from var files, then make the attempt to read the var files. 
  ! If var files are missing, the code above has horizontally homogeneously 
  ! initialized the prognostic vars and we will continue with that initialization.
  if (ifg_mode .eq. IFG_VF_INIT) then
    CALL varf_update (0,ifileok,1)

    if (ifileok  ==  1) then
      ! Everything's cool...
      if (print_msg) print*,'VarRead Initial varfile read of grid-',ifm
    else
      ! Using interpolated nudging arrays from parent grid.
      if (print_msg) print*,'VarRead Initial interpolation of grid-',ifm
    endif
  endif

  ! Interpolate 3D density, considering terrain this time
  ! For parallel runs, calc_dn0uv() uses the data in the overlap region 
  ! for dn0 so need to update
  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT) &
                                  .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL interp_topo (nzf,nxf,nyf,basic_g(ifm)%dn0,toptf_int,grid_g(ifm)%topt, &
                     ifm,'t',ibconf,3)
    CALL calc_dn0uv (nzf,nxf,nyf,basic_g(ifm)%dn0,basic_g(ifm)%dn0u, &
                    basic_g(ifm)%dn0v,ibconf,ifm)
  elseif (ifg_mode .eq. IFG_MKVF_INIT) then
    CALL interp_topo (nzf,nxf,nyf,is_grids(ifm)%rr_dn0,toptf_int, &
                     grid_g(ifm)%topt,ifm,'t',ibconf,3)
    CALL calc_dn0uv (nzf,nxf,nyf,is_grids(ifm)%rr_dn0,is_grids(ifm)%rr_dn0u, &
                    is_grids(ifm)%rr_dn0v,ibconf,ifm)
  endif
 

  if (((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT)) &
     .and. (print_msg)) then
    print'(a)',    '----------------------------------------------------'            
    print'(a,i5)', 'Model Start Initial interpolation of grid-',ifm
    print'(a)',    '----------------------------------------------------'
    print'(a)',    ''
  endif

  deallocate(dn01dnct)
  deallocate(toptf_int)
  deallocate(dn0ct)
  deallocate(dn0uct)
  deallocate(dn0vct)

  return
END SUBROUTINE interp_fine_grid

!##############################################################################
Subroutine interp_fine_grid_sfcvar (ifm,icm,ifg_mode)

! This routine will perform initialization of a fine grid by interpolating
! from the parent (coarse) grid.
!
! The argument ifg_mode denotes which subsystem (consequently which variable
! set) to initialize. The following modes are currently supported, and the
! corresponding parameters are found in node_mode.f90.
!
!    parameter           mode
!
!  IFG_MKVF_INIT     initialization for MAKEVFILE
!  IFG_MKHF_INIT     initialization for MAKEHFILE
!  IFG_TOPT1_INIT    initialization for MAKESFC, topography
!  IFG_TOPT2_INIT    initialization for MAKESFC, topography
!  IFG_SST_INIT      initialization for MAKESFC, sea surface temperature
!
! It is the caller's responsibility to check that ifm and icm are a valid
! coarse-fine grid pair. The caller also needs to do a call to newgrid(ifm).

  use node_mod
  use mem_grid
  use mem_basic
  use mem_mksfc
  use ref_sounding
  use var_tables
  use isan_coms

  implicit none

  integer :: ifm, icm
  integer :: ifg_mode

  integer :: nzc, nxc, nyc
  integer :: nzct, nxct, nyct
  integer :: nzf, nxf, nyf
  integer :: ic0, jc0
  integer :: ict0, jct0
  integer :: if0, jf0
  integer :: ibconf
  integer :: ifileok
  integer :: i, j, k, ivar
  real, dimension(:,:,:), allocatable :: dn0ct, dn0uct, dn0vct
  real, dimension(:), allocatable :: dn01dnct
  real, dimension(:,:), allocatable :: toptf_int

  ! These are for convenience
  nxf = mmxp(ifm)
  nyf = mmyp(ifm)
  nzf = mmzp(ifm)

  nxc = mmxp(icm)
  nyc = mmyp(icm)
  nzc = mmzp(icm)

  nxct = ctnxp(ifm)
  nyct = ctnyp(ifm)
  nzct = ctnzp(ifm)

  ic0 = mi0(icm)
  jc0 = mj0(icm)

  ict0 = cti0(ifm)
  jct0 = ctj0(ifm)

  if0 = mi0(ifm)
  jf0 = mj0(ifm)

  ibconf = mibcon(ifm)

  allocate(dn01dnct(nzct))
  allocate(toptf_int(nxf,nyf))
  allocate(dn0ct(nzct,nxct,nyct))
  allocate(dn0uct(nzct,nxct,nyct))
  allocate(dn0vct(nzct,nxct,nyct))

  ! Fill the "tiles" for the density vars
  if ((ifg_mode .eq. IFG_MKVF_INIT) .or. (ifg_mode .eq. IFG_MKHF_INIT)) then
    CALL fill_density_tiles (nzc,nxc,nyc,dnref(1,icm), &
        is_grids(icm)%rr_dn0,is_grids(icm)%rr_dn0u,is_grids(icm)%rr_dn0v, &
        nzct,nxct,nyct,dn01dnct,dn0ct,dn0uct,dn0vct,ic0,jc0,ict0,jct0,ifm,icm)
  endif


  ! Interpolate topography for subsequent interp_var calls, or for MAKESFC
  if ((ifg_mode .eq. IFG_TOPT1_INIT) .or. (ifg_mode .eq. IFG_TOPT2_INIT)) then
    CALL interp_varsfc (1,nxc,nyc,sfcfile_p(icm)%topt,1,nxct,nyct, &
                    1,nxf,nyf,sfcfile_p(ifm)%topt,  &
                    ifm,icm,ctipm(1,ifm), &
                    ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                    ibconf,'t',2)
    if (ifg_mode .eq. IFG_TOPT1_INIT) then
      CALL interp_varsfc (1,nxc,nyc,sfcfile_p(icm)%topzo,1,nxct,nyct, &
                      1,nxf,nyf,sfcfile_p(ifm)%topzo,  &
                      ifm,icm,ctipm(1,ifm),  &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                      ibconf,'t',2)
    endif
  else
    CALL interp_varsfc (1,nxc,nyc,grid_g(icm)%topt,1,nxct,nyct, &
                    1,nxf,nyf,toptf_int, &
                    ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm), &
                    ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',2)
  endif

  ! Interpolate SST for MAKESFC
  if (ifg_mode .eq. IFG_SST_INIT) then
    CALL interp_varsfc (1,nxc,nyc,sfcfile_p(icm)%seatf,1,nxct,nyct, &
                    1,nxf,nyf,sfcfile_p(ifm)%seatf, &
                    ifm,icm,ctipm(1,ifm), &
                    ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                    ibconf,'t',2)
  endif


  ! Interpolate the 1D reference state
  if (ifg_mode .eq. IFG_MKVF_INIT) then
    CALL interp_refs_1d (nzc,dnref(1,icm),thref(1,icm),u01dn(1,icm),            &
           v01dn(1,icm),rtref(1,icm),piref(1,icm),nzct,dn01dnct,                &
           nzf,nxf,nyf,dnref(1,ifm),thref(1,ifm),u01dn(1,ifm),v01dn(1,ifm),     &
           rtref(1,ifm),piref(1,ifm),grid_g(ifm)%topt,toptf_int,                &
           ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ibconf,ifg_mode)
  endif


  ! Interpolate the 3D reference state
  if ((ifg_mode .eq. IFG_MKVF_INIT) .or. (ifg_mode .eq. IFG_MKHF_INIT)) then
    CALL interp_refs_3d (nzc,nxc,nyc,is_grids(icm)%rr_dn0,is_grids(icm)%rr_th0,  &
          nzct,nxct,nyct,dn0ct,nzf,nxf,nyf,is_grids(ifm)%rr_dn0,                &
          is_grids(ifm)%rr_dn0u,is_grids(ifm)%rr_dn0v,is_grids(ifm)%rr_th0,     &
          is_grids(ifm)%rr_pi0,is_grids(ifm)%rr_p,grid_g(ifm)%topt,toptf_int,   &
          ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0, &
          jf0,ibconf,ifg_mode)
  endif

  ! Interpolate 3D density, considering terrain this time
  ! For parallel runs, calc_dn0uv() uses the data in the overlap region 
  ! for dn0 so need to update
  if (ifg_mode .eq. IFG_MKVF_INIT) then
    CALL interp_topo (nzf,nxf,nyf,is_grids(ifm)%rr_dn0,toptf_int, &
                     grid_g(ifm)%topt,ifm,'t',ibconf,3)
    CALL calc_dn0uv (nzf,nxf,nyf,is_grids(ifm)%rr_dn0,is_grids(ifm)%rr_dn0u, &
                    is_grids(ifm)%rr_dn0v,ibconf,ifm)
  endif
 

  deallocate(dn01dnct)
  deallocate(toptf_int)
  deallocate(dn0ct)
  deallocate(dn0uct)
  deallocate(dn0vct)

  return
END SUBROUTINE interp_fine_grid_sfcvar

!##############################################################################
Subroutine interp_varf (ifm, ifflag)

! This routine will perform initialization of fine grids by interpolating
! from the parent (coarse) grid for variables in the var file system
! (varu[pf], varv[pf], varp[pf], vart[pf], varr[pf]).

  use node_mod
  use mem_grid
  use mem_basic
  use mem_varinit

  implicit none

  integer :: ifm, ifflag
  integer :: icm
  integer :: nzc, nxc, nyc
  integer :: nzct, nxct, nyct
  integer :: nzf, nxf, nyf
  integer :: ic0, jc0
  integer :: ict0, jct0
  integer :: if0, jf0
  integer :: ibconf
  integer :: ifileok
  real :: c1, c2
  integer :: i, j, k, ivar
  real, dimension(:,:,:), allocatable :: dn0ct, dn0uct, dn0vct
  real, dimension(:,:), allocatable :: toptf_int

  icm = nxtnest(ifm)
  if (icm == 0) return

  CALL newgrid (ifm)

  ! These are for convenience
  nxf = mmxp(ifm)
  nyf = mmyp(ifm)
  nzf = mmzp(ifm)

  nxc = mmxp(icm)
  nyc = mmyp(icm)
  nzc = mmzp(icm)

  nxct = ctnxp(ifm)
  nyct = ctnyp(ifm)
  nzct = ctnzp(ifm)

  ic0 = mi0(icm)
  jc0 = mj0(icm)

  ict0 = cti0(ifm)
  jct0 = ctj0(ifm)

  if0 = mi0(ifm)
  jf0 = mj0(ifm)

  ibconf = mibcon(ifm)

  allocate(toptf_int(nxf,nyf))
  allocate(dn0ct(nzct,nxct,nyct))
  allocate(dn0uct(nzct,nxct,nyct))
  allocate(dn0vct(nzct,nxct,nyct))

  ! Grab the coarse tile density vars
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0,nzct,nxct,nyct,dn0ct, &
                         ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0u,nzct,nxct,nyct,dn0uct, &
                         ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0v,nzct,nxct,nyct,dn0vct, &
                         ic0,jc0,ict0,jct0,ifm,icm,3)

  ! Interpolate topography for subsequent interp_var calls
  CALL interp_var (1,nxc,nyc,grid_g(icm)%topt,1,nxct,nyct,dn0ct, &
                  1,nxf,nyf,toptf_int,basic_g(ifm)%dn0,grid_g(ifm)%topt, &
                  toptf_int,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm),ctkpm(1,ifm), &
                  ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',0,0,2)

  ! ifflag == 1 --> interpolate the var weights
  ! ifflag == 2 --> interpolate the future level atmospheric variables

  select case (ifflag)
    case (1)
      ! var weights
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%varwts,nzct,nxct,nyct,dn0ct, &
                      nzf,nxf,nyf,varinit_g(ifm)%varwts,basic_g(ifm)%dn0, &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm), &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0, &
                      ibconf,'t',0,1,3)

    case (2)
      ! future level vars

      ! U
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%varuf,nzct,nxct,nyct,dn0uct, &
                      nzf,nxf,nyf,varinit_g(ifm)%varuf,basic_g(ifm)%dn0u,      &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),         &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,     &
                      ibconf,'u',1,1,3)

      ! V
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%varvf,nzct,nxct,nyct,dn0vct, &
                      nzf,nxf,nyf,varinit_g(ifm)%varvf,basic_g(ifm)%dn0v,      &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),         &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,     &
                      ibconf,'v',1,1,3)

      ! P
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%varpf,nzct,nxct,nyct,dn0ct, &
                      nzf,nxf,nyf,varinit_g(ifm)%varpf,basic_g(ifm)%dn0,      &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),        &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,    &
                      ibconf,'t',0,1,3)

      ! T
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%vartf,nzct,nxct,nyct,dn0ct, &
                      nzf,nxf,nyf,varinit_g(ifm)%vartf,basic_g(ifm)%dn0,      &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),        &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,    &
                      ibconf,'t',1,1,3)

      ! R
      CALL interp_var (nzc,nxc,nyc,varinit_g(icm)%varrf,nzct,nxct,nyct,dn0ct, &
                      nzf,nxf,nyf,varinit_g(ifm)%varrf,basic_g(ifm)%dn0,      &
                      grid_g(ifm)%topt,toptf_int,ifm,icm,ctipm(1,ifm),        &
                      ctjpm(1,ifm),ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,    &
                      ibconf,'t',1,1,3)

  endselect

  deallocate(toptf_int)
  deallocate(dn0ct)
  deallocate(dn0uct)
  deallocate(dn0vct)

  return
END SUBROUTINE interp_varf

!##############################################################################
Subroutine interp_nest_bounds (icm,ifm)

! This routine will perform initialization of fine grids by interpolating
! from the parent (coarse) grid.

  use node_mod
  use mem_grid
  use mem_basic
  use mem_nestb
  use var_tables

  implicit none

  integer :: icm, ifm

  integer :: nxf, nyf, nzf, nxc, nyc, nzc, nxct, nyct, nzct, ic0, jc0, ict0, jct0, if0, jf0, ivar, ibconf
  real, dimension(:,:,:), allocatable :: dn0ct, dn0uct, dn0vct, varf

  ! These are for convenience
  nxf = mmxp(ifm)
  nyf = mmyp(ifm)
  nzf = mmzp(ifm)

  nxc = mmxp(icm)
  nyc = mmyp(icm)
  nzc = mmzp(icm)

  nxct = ctnxp(ifm)
  nyct = ctnyp(ifm)
  nzct = ctnzp(ifm)

  ic0 = mi0(icm)
  jc0 = mj0(icm)

  ict0 = cti0(ifm)
  jct0 = ctj0(ifm)

  if0 = mi0(ifm)
  jf0 = mj0(ifm)

  ibconf = mibcon(ifm)

  allocate(dn0ct(nzct,nxct,nyct))
  allocate(dn0uct(nzct,nxct,nyct))
  allocate(dn0vct(nzct,nxct,nyct))

  allocate(varf(nzf,nxf,nyf))

  ! Interpolate the prognostic variables (u, v, w, p, and scalars), then copy
  ! tile edges to domain boundaries.
  ! Since the topology interpolation is not being used (irtgflg == 0) in the 
  ! calls to interp_var, use grid_g(ifm)%topt for both topt arguments.
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0,nzct,nxct,nyct,dn0ct   &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0u,nzct,nxct,nyct,dn0uct &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,basic_g(icm)%dn0v,nzct,nxct,nyct,dn0vct &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)

  ! U
  CALL interp_var (nzc,nxc,nyc,basic_g(icm)%uc,nzct,nxct,nyct,dn0uct        &
                 ,nzf,nxf,nyf,varf,basic_g(ifm)%dn0u,grid_g(ifm)%topt       &
                 ,grid_g(ifm)%topt,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm)        &
                 ,ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'u',1,0,3)
  CALL copy_nest_bounds (nzf,nxf,nyf,varf,nbounds(ifm)%bux,nbounds(ifm)%buy &
                        ,nbounds(ifm)%buz,'u',nstbot,nsttop,ibcon)

  ! V
  CALL interp_var (nzc,nxc,nyc,basic_g(icm)%vc,nzct,nxct,nyct,dn0vct        &
                 ,nzf,nxf,nyf,varf,basic_g(ifm)%dn0v,grid_g(ifm)%topt       &
                 ,grid_g(ifm)%topt,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm)        &
                 ,ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'v',1,0,3)
  CALL copy_nest_bounds (nzf,nxf,nyf,varf,nbounds(ifm)%bvx,nbounds(ifm)%bvy &
                        ,nbounds(ifm)%bvz,'v',nstbot,nsttop,ibcon)

  ! W
  CALL interp_var (nzc,nxc,nyc,basic_g(icm)%wc,nzct,nxct,nyct,dn0ct         &
                 ,nzf,nxf,nyf,varf,basic_g(ifm)%dn0,grid_g(ifm)%topt        &
                 ,grid_g(ifm)%topt,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm)        &
                 ,ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'w',1,0,3)
  CALL copy_nest_bounds (nzf,nxf,nyf,varf,nbounds(ifm)%bwx,nbounds(ifm)%bwy &
                        ,nbounds(ifm)%bwz,'w',nstbot,nsttop,ibcon)

  ! P
  CALL interp_var (nzc,nxc,nyc,basic_g(icm)%pc,nzct,nxct,nyct,dn0ct         &
                 ,nzf,nxf,nyf,varf,basic_g(ifm)%dn0,grid_g(ifm)%topt        &
                 ,grid_g(ifm)%topt,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm)        &
                 ,ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',0,0,3)
  CALL copy_nest_bounds (nzf,nxf,nyf,varf,nbounds(ifm)%bpx,nbounds(ifm)%bpy &
                        ,nbounds(ifm)%bpz,'t',nstbot,nsttop,ibcon)

  do ivar = 1,num_scalar(ifm)
    CALL interp_var (nzc,nxc,nyc,scalar_tab(ivar,icm)%var_p,nzct,nxct,nyct   &
                   ,dn0ct,nzf,nxf,nyf,varf,basic_g(ifm)%dn0,grid_g(ifm)%topt &
                   ,grid_g(ifm)%topt,ifm,icm,ctipm(1,ifm),ctjpm(1,ifm)       &
                   ,ctkpm(1,ifm),ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',1,0,3)
    CALL copy_nest_bounds (nzf,nxf,nyf,varf,nbounds(ifm)%bsx(1,1,1,ivar)     &
            ,nbounds(ifm)%bsy(1,1,1,ivar),nbounds(ifm)%bsz(1,1,1,ivar),'t'   &
            ,nstbot,nsttop,ibcon)
  enddo

  deallocate(dn0ct)
  deallocate(dn0uct)
  deallocate(dn0vct)

  deallocate(varf)

  return
END SUBROUTINE interp_nest_bounds

