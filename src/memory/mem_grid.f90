!##############################################################################
Module mem_grid

use grid_dims

implicit none

   Type grid_vars
                            
      ! Variables to be dimensioned by (nxp,nyp)
   real, allocatable, dimension(:,:) :: &
                 topt,  topu,  topv,  topm  &
              ,  rtgt,  rtgu,  rtgv,  rtgm  &
              ,  f13t,  f13u,  f13v,  f13m, f23t , f23u , f23v , f23m  &
              ,   dxt,   dxu,   dxv,   dxm,  dyt ,  dyu ,  dyv ,  dym  & 
              , fmapt, fmapu, fmapv, fmapm,fmapti,fmapui,fmapvi,fmapmi &
              ,  glat,  glon, topzo
              
   End Type
   
   type (grid_vars), allocatable :: grid_g(:), gridm_g(:)

   !-------------------------------------------------------------------------------
   character(len=strl1) :: expnme,ramsorrevu
   !-------------------------------------------------------------------------------
   integer :: ngrids,ngridsh
   !-------------------------------------------------------------------------------
   integer, dimension(maxgrds) :: nnxp,nnyp,nnzp                              

   !-------------------------------------------------------------------------------
   integer, dimension(maxgrds) :: nnx,nnx1,nnx2,nny,nny1,nny2,nnz,nnz1  &
                              ,nnxyzp,nnxysp,nnxyp
   !-------------------------------------------------------------------------------
   real, dimension(nzpmax,maxgrds) :: htn,ht2n,ht4n,hwn,hw2n,hw4n
   !-------------------------------------------------------------------------------
   real, dimension(nzpmax,maxgrds) :: dztn,dzmn,dzt2n,dzm2n,ztn,zmn
   real, dimension(nxpmax,maxgrds) :: xtn,xmn
   real, dimension(nypmax,maxgrds) :: ytn,ymn
   !-------------------------------------------------------------------------------
   integer :: nxp,nx,nx1,nx2,nyp,ny,ny1,ny2,nzp,nzpp,nz,nz1  &
             ,nxyzp,nxyp,nxysp,nsttop,nstbot,ndtrat,jdim
   !-------------------------------------------------------------------------------
   real                    :: deltax,deltaz
   real, dimension(nzpmax) :: ht,ht2,ht4,hw,hw2,hw4,zt,zm,dzt,dzm,dzt2,dzm2
   real, dimension(nxpmax) :: xt,xm
   real, dimension(nypmax) :: yt,ym
   
   integer :: nzg,nzs,npatch 

   integer :: itopo,ihtran

   !-------------------------------------------------------------------------------
   integer :: ngrid,ngridc,ngrido,iscr1,iscr2
   !-------------------------------------------------------------------------------
   integer                    :: memsize,iounit,maxpro,memscr,memind  &
                                ,iogrid,maxpts,maxnzp,maxnxp,maxnyp,i2dvar
   integer,dimension(maxgrds) :: memgrd

   !-------------------------------------------------------------------------------
   real                     :: time,ztop,dzrat,dzmax,zdelay
   integer                  :: impl,iyear1,imonth1,idate1,ihour1  &
                              ,itime1,nacoust,initial,initorig,iflag,hrestart
   integer, dimension(maxgrds) :: nnacoust                              
   real, dimension(maxgrds) :: deltazn,deltaxn,dtlongn
   real, dimension(nzpmax)  :: zz

   real                     :: dtlong,sspct,polelat,polelon
   real, dimension(maxgrds) :: centlat,centlon  &
                              ,cflxy,cflz,cflxy_max,cflz_max &
                              ,vertvel,vertvel_max
   !-------------------------------------------------------------------------------
   character (len=16) :: runtype
   character(len=1)   :: timeunit
   !-------------------------------------------------------------------------------
   integer :: isstp,istp
   real    :: timmax,dts,dtlt,dtlv

   integer                   :: nestz,iprntstmt
   logical                   :: print_msg
   integer, dimension(maxgrds) :: nndtrat,nstratx,nstraty  &
                                 ,ngbegun,nxtnest,nnsttop,nnstbot  &
                                 ,ninest,njnest,nknest
   integer, dimension(nzpmax)  :: nstratz
   
   !-------------------------------------------------------------------------------
   integer, parameter                     :: maxsched=2000,maxschent=5
   integer                                :: nsubs
   integer, dimension(maxsched,maxschent) :: isched

   !---------------------------------------------------------------------------
   integer, dimension(nxpmax,maxgrds) :: ipm
   integer, dimension(nypmax,maxgrds) :: jpm 
   integer, dimension(nzpmax,maxgrds) :: nrz,kpm 
   real, dimension(nxpmax,maxgrds)    :: ei1,ei2,ei3,ei4,ei5,ei6,ei7 
   real, dimension(nypmax,maxgrds)    :: ej1,ej2,ej3,ej4,ej5,ej6,ej7
   real, dimension(nzpmax,maxgrds)    :: ek1,ek2,ek3,ek4,ek5,ek6,ek7 
   real, dimension(nzpmax,maxgrds,4)  :: fbcf
   !---------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   integer                  :: lsflg,ibnd,jbnd,icorflg,nfpt
   integer, dimension(maxgrds) :: isponge_pts
   !-------------------------------------------------------------------------------
   real                     :: distim ,cphas
   real, dimension(maxgrds) :: sponge_tau
   !-------------------------------------------------------------------------------

Contains

!##############################################################################
Subroutine alloc_grid (grid,n2,n3)

implicit none

   type (grid_vars) :: grid
   integer, intent(in) :: n2,n3

! Allocate arrays based on options (if necessary)
      
     allocate (grid%topt(n2,n3))
     allocate (grid%topu(n2,n3))
     allocate (grid%topv(n2,n3))
     allocate (grid%topm(n2,n3))
     allocate (grid%rtgt(n2,n3))
     allocate (grid%rtgu(n2,n3))
     allocate (grid%rtgv(n2,n3))
     allocate (grid%rtgm(n2,n3))
     allocate (grid%f13t(n2,n3))
     allocate (grid%f13u(n2,n3))
     allocate (grid%f13v(n2,n3))
     allocate (grid%f13m(n2,n3))
     allocate (grid%f23t(n2,n3))
     allocate (grid%f23u(n2,n3))
     allocate (grid%f23v(n2,n3))
     allocate (grid%f23m(n2,n3))
     allocate (grid%dxt(n2,n3))
     allocate (grid%dxu(n2,n3))
     allocate (grid%dxv(n2,n3))
     allocate (grid%dxm(n2,n3))
     allocate (grid%dyt(n2,n3))
     allocate (grid%dyu(n2,n3))
     allocate (grid%dyv(n2,n3))
     allocate (grid%dym(n2,n3))
     allocate (grid%fmapt(n2,n3))
     allocate (grid%fmapu(n2,n3))
     allocate (grid%fmapv(n2,n3))
     allocate (grid%fmapm(n2,n3))
     allocate (grid%fmapti(n2,n3))
     allocate (grid%fmapui(n2,n3))
     allocate (grid%fmapvi(n2,n3))
     allocate (grid%fmapmi(n2,n3))
     allocate (grid%glat(n2,n3))
     allocate (grid%glon(n2,n3))
     allocate (grid%topzo(n2,n3))
    
return
END SUBROUTINE alloc_grid

!##############################################################################
Subroutine dealloc_grid (grid)

implicit none

   type (grid_vars) :: grid

! Deallocate arrays
      
     if(allocated(grid%topt   ))    deallocate (grid%topt)
     if(allocated(grid%topu   ))    deallocate (grid%topu)
     if(allocated(grid%topv   ))    deallocate (grid%topv)
     if(allocated(grid%topm   ))    deallocate (grid%topm)
     if(allocated(grid%rtgt   ))    deallocate (grid%rtgt)
     if(allocated(grid%rtgu   ))    deallocate (grid%rtgu)
     if(allocated(grid%rtgv   ))    deallocate (grid%rtgv)
     if(allocated(grid%rtgm   ))    deallocate (grid%rtgm)
     if(allocated(grid%f13t   ))    deallocate (grid%f13t)
     if(allocated(grid%f13u   ))    deallocate (grid%f13u)
     if(allocated(grid%f13v   ))    deallocate (grid%f13v)
     if(allocated(grid%f13m   ))    deallocate (grid%f13m)
     if(allocated(grid%f23t   ))    deallocate (grid%f23t)
     if(allocated(grid%f23u   ))    deallocate (grid%f23u)
     if(allocated(grid%f23v   ))    deallocate (grid%f23v)
     if(allocated(grid%f23m   ))    deallocate (grid%f23m)
     if(allocated(grid%dxt    ))    deallocate (grid%dxt)
     if(allocated(grid%dxu    ))    deallocate (grid%dxu)
     if(allocated(grid%dxv    ))    deallocate (grid%dxv)
     if(allocated(grid%dxm    ))    deallocate (grid%dxm)
     if(allocated(grid%dyt    ))    deallocate (grid%dyt)
     if(allocated(grid%dyu    ))    deallocate (grid%dyu)
     if(allocated(grid%dyv    ))    deallocate (grid%dyv)
     if(allocated(grid%dym    ))    deallocate (grid%dym)
     if(allocated(grid%fmapt  ))    deallocate (grid%fmapt)
     if(allocated(grid%fmapu  ))    deallocate (grid%fmapu)
     if(allocated(grid%fmapv  ))    deallocate (grid%fmapv)
     if(allocated(grid%fmapm  ))    deallocate (grid%fmapm)
     if(allocated(grid%fmapti ))    deallocate (grid%fmapti)
     if(allocated(grid%fmapui ))    deallocate (grid%fmapui)
     if(allocated(grid%fmapvi ))    deallocate (grid%fmapvi)
     if(allocated(grid%fmapmi ))    deallocate (grid%fmapmi)
     if(allocated(grid%glat   ))    deallocate (grid%glat)
     if(allocated(grid%glon   ))    deallocate (grid%glon)
     if(allocated(grid%topzo  ))    deallocate (grid%topzo)
    
return
END SUBROUTINE dealloc_grid

!##############################################################################
Subroutine filltab_grid (grid,gridm,imean,n2,n3,ng)

use var_tables

implicit none

   type (grid_vars) :: grid,gridm
   integer, intent(in) :: imean,n2,n3,ng
   integer :: npts

! Fill arrays into variable tables

   npts=n2*n3
   if (allocated(grid%topt)) &
      CALL vtables2 (grid%topt(1,1),gridm%topt(1,1),ng,npts,imean,  &
                 'TOPT :2:anal:mpti')      
   if (allocated(grid%topu)) &
      CALL vtables2 (grid%topu(1,1),gridm%topu(1,1),ng, npts, imean,  &
                 'TOPU :2:mpti')      
   if (allocated(grid%topv)) &
      CALL vtables2 (grid%topv(1,1),gridm%topv(1,1),ng, npts, imean,  &
                 'TOPV :2:mpti')      
   if (allocated(grid%topm)) &
      CALL vtables2 (grid%topm(1,1),gridm%topm(1,1),ng, npts, imean,  &
                 'TOPM :2:mpti')      
   if (allocated(grid%rtgt)) &
      CALL vtables2 (grid%rtgt(1,1),gridm%rtgt(1,1),ng, npts, imean,  &
                 'RTGT :2:mpti')      
   if (allocated(grid%rtgu)) &
      CALL vtables2 (grid%rtgu(1,1),gridm%rtgu(1,1),ng, npts, imean,  &
                 'RTGU :2:mpti')      
   if (allocated(grid%rtgv)) &
      CALL vtables2 (grid%rtgv(1,1),gridm%rtgv(1,1),ng, npts, imean,  &
                 'RTGV :2:mpti')      
   if (allocated(grid%rtgm)) &
      CALL vtables2 (grid%rtgm(1,1),gridm%rtgm(1,1),ng, npts, imean,  &
                 'RTGM :2:mpti')      
   if (allocated(grid%f13t)) &
      CALL vtables2 (grid%f13t(1,1),gridm%f13t(1,1),ng, npts, imean,  &
                 'F13T :2:mpti')      
   if (allocated(grid%f13u)) &
      CALL vtables2 (grid%f13u(1,1),gridm%f13u(1,1),ng, npts, imean,  &
                 'F13U :2:mpti')      
   if (allocated(grid%f13v)) &
      CALL vtables2 (grid%f13v(1,1),gridm%f13v(1,1),ng, npts, imean,  &
                 'F13V :2:mpti')      
   if (allocated(grid%f13m)) &
      CALL vtables2 (grid%f13m(1,1),gridm%f13m(1,1),ng, npts, imean,  &
                 'F13M :2:mpti')      
   if (allocated(grid%f23t)) &
      CALL vtables2 (grid%f23t(1,1),gridm%f23t(1,1),ng, npts, imean,  &
                 'F23T :2:mpti')      
   if (allocated(grid%f23u)) &
      CALL vtables2 (grid%f23u(1,1),gridm%f23u(1,1),ng, npts, imean,  &
                 'F23U :2:mpti')      
   if (allocated(grid%f23v)) &
      CALL vtables2 (grid%f23v(1,1),gridm%f23v(1,1),ng, npts, imean,  &
                 'F23V :2:mpti')      
   if (allocated(grid%f23m)) &
      CALL vtables2 (grid%f23m(1,1),gridm%f23m(1,1),ng, npts, imean,  &
                 'F23M :2:mpti')      
   if (allocated(grid%dxt)) &
      CALL vtables2 (grid%dxt(1,1),gridm%dxt(1,1),ng, npts, imean,  &
                 'DXT :2:mpti')      
   if (allocated(grid%dxu)) &
      CALL vtables2 (grid%dxu(1,1),gridm%dxu(1,1),ng, npts, imean,  &
                 'DXU :2:mpti')      
   if (allocated(grid%dxv)) &
      CALL vtables2 (grid%dxv(1,1),gridm%dxv(1,1),ng, npts, imean,  &
                 'DXV :2:mpti')      
   if (allocated(grid%dxm)) &
      CALL vtables2 (grid%dxm(1,1),gridm%dxm(1,1),ng, npts, imean,  &
                 'DXM :2:mpti')      
   if (allocated(grid%dyt)) &
      CALL vtables2 (grid%dyt(1,1),gridm%dyt(1,1),ng, npts, imean,  &
                 'DYT :2:mpti')      
   if (allocated(grid%dyu)) &
      CALL vtables2 (grid%dyu(1,1),gridm%dyu(1,1),ng, npts, imean,  &
                 'DYU :2:mpti')      
   if (allocated(grid%dyv)) &
      CALL vtables2 (grid%dyv(1,1),gridm%dyv(1,1),ng, npts, imean,  &
                 'DYV :2:mpti')      
   if (allocated(grid%dym)) &
      CALL vtables2 (grid%dym(1,1),gridm%dym(1,1),ng, npts, imean,  &
                 'DYM :2:mpti')      
   if (allocated(grid%fmapt)) &
      CALL vtables2 (grid%fmapt(1,1),gridm%fmapt(1,1),ng, npts, imean,  &
                 'FMAPT :2:mpti')      
   if (allocated(grid%fmapu)) &
      CALL vtables2 (grid%fmapu(1,1),gridm%fmapu(1,1),ng, npts, imean,  &
                 'FMAPU :2:mpti')      
   if (allocated(grid%fmapv)) &
      CALL vtables2 (grid%fmapv(1,1),gridm%fmapv(1,1),ng, npts, imean,  &
                 'FMAPV :2:mpti')      
   if (allocated(grid%fmapm)) &
      CALL vtables2 (grid%fmapm(1,1),gridm%fmapm(1,1),ng, npts, imean,  &
                 'FMAPM :2:mpti')      
   if (allocated(grid%fmapti)) &
      CALL vtables2 (grid%fmapti(1,1),gridm%fmapti(1,1),ng, npts, imean,  &
                 'FMAPTI :2:mpti')      
   if (allocated(grid%fmapui)) &
      CALL vtables2 (grid%fmapui(1,1),gridm%fmapui(1,1),ng, npts, imean,  &
                 'FMAPUI :2:mpti')      
   if (allocated(grid%fmapvi)) &
      CALL vtables2 (grid%fmapvi(1,1),gridm%fmapvi(1,1),ng, npts, imean,  &
                 'FMAPVI :2:mpti')      
   if (allocated(grid%fmapmi)) &
      CALL vtables2 (grid%fmapmi(1,1),gridm%fmapmi(1,1),ng, npts, imean,  &
                 'FMAPMI :2:mpti')      
   if (allocated(grid%glat)) &
      CALL vtables2 (grid%glat(1,1),gridm%glat(1,1),ng, npts, imean,  &
                 'GLAT :2:mpti:anal')      
   if (allocated(grid%glon)) &
      CALL vtables2 (grid%glon(1,1),gridm%glon(1,1),ng, npts, imean,  &
                 'GLON :2:mpti:anal')      
   if (allocated(grid%topzo)) &
      CALL vtables2 (grid%topzo(1,1),gridm%topzo(1,1),ng, npts, imean,  &
                 'TOPZO :2:mpti:anal')

return
END SUBROUTINE filltab_grid

END MODULE mem_grid
