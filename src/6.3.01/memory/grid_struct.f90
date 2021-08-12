!##############################################################################
Module grid_struct

implicit none

Type grid_def

integer :: nxp,nyp,nzp,nzg,nzs,npatch,nvegpat
real :: plat, plon
real, allocatable, dimension(:) :: xtn,ytn,xmn,ymn,ztn,zmn
real, allocatable, dimension(:,:) :: topo

End Type

Contains

!##############################################################################
Subroutine alloc_grid_def (grd,nx,ny,nz)

implicit none

type(grid_def) :: grd

integer :: nx,ny,nz

allocate (grd%xtn(nx),grd%ytn(ny),grd%xmn(nx),grd%ymn(ny),grd%ztn(nz),grd%zmn(nz))
allocate (grd%topo(nx,ny))

return
END SUBROUTINE alloc_grid_def

!##############################################################################
Subroutine dealloc_grid_def (grd)

implicit none

type(grid_def) :: grd

deallocate (grd%xtn,grd%ytn,grd%xmn,grd%ymn,grd%ztn,grd%zmn)
deallocate (grd%topo)

return
END SUBROUTINE dealloc_grid_def

!##############################################################################
Subroutine fill_grid_def (grd,nx,ny,nz,ng,ns,np,npv,plt,pln  &
                           ,xt,xm,yt,ym,zt,zm,topo)

! Fill grid_def structure - must call alloc_grid_def first

implicit none

type(grid_def) :: grd

integer :: nx,ny,nz,ng,ns,np,npv
real    :: plt,pln, xt(nx),yt(ny),zt(nz)  &
                   ,xm(nx),ym(ny),zm(nz),topo(nx,ny)

grd%nxp=nx   ; grd%nyp=ny   ; grd%nzp=nz 
grd%nzg=ng   ; grd%nzs=ns   ; grd%npatch=np ; grd%nvegpat=npv
grd%plat=plt ; grd%plon=pln

grd%xtn(1:nx)=xt(1:nx)
grd%xmn(1:nx)=xm(1:nx)
grd%ytn(1:ny)=yt(1:ny)
grd%ymn(1:ny)=ym(1:ny)
grd%ztn(1:nz)=zt(1:nz)
grd%zmn(1:nz)=zm(1:nz)

grd%topo(1:nx,1:ny)=topo(1:nx,1:ny)

return
END SUBROUTINE fill_grid_def

!##############################################################################
Subroutine compare_grid_def (g1,g2,string,ierr,hngr,cngr)

! Compare two grids to see if they are the same structure 
!   Input: g1,g2 - grid_def structures filled with call to fill_grid_def
!          string - identifier string to be printed on errors

implicit none

type(grid_def) :: g1, g2
integer :: ierr,ict,nd,hngr,cngr
character(len=*) :: string
integer, external :: check_real

ict=0
ierr=0
print*,'###################################################################'
print*,'### COMPARING HISTORY AND CURRENT GRIDS FOR MATCHING'
print*,'### HISTORY GRID:',hngr,' CURRENT GRID:',cngr
print*,'###################################################################'
if (g1%nxp /= g2%nxp ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nxp different',g1%nxp,g2%nxp
endif

if (g1%nyp /= g2%nyp ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nyp different',g1%nyp,g2%nyp
endif

if (g1%nzp /= g2%nzp ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nzp different',g1%nzp,g2%nzp
endif

if (g1%nzg /= g2%nzg ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nzg different',g1%nzg,g2%nzg
endif

if (g1%nzs /= g2%nzs ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nzs different',g1%nzs,g2%nzs
endif

if (g1%npatch /= g2%npatch ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': npatch different',g1%npatch,g2%npatch
endif

if (g1%nvegpat /= g2%nvegpat ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': nvegpat different',g1%nvegpat,g2%nvegpat
endif

nd=check_real(g1%plat,g2%plat,1,1)
if ( nd > 0 ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': plat different',g1%plat,g2%plat
endif

nd=check_real(g1%plon,g2%plon,1,1)
if ( nd > 0 ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': plon different',g1%plon,g2%plon
endif

nd=check_real(g1%xtn,g2%xtn,g1%nxp,g2%nxp)
if ( nd > 0 ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': xtn different',g1%xtn(nd),g2%xtn(nd)
endif

nd=check_real(g1%ytn,g2%ytn,g1%nyp,g2%nyp)
if ( nd > 0 ) then
   ict=ict+1
   print*,'Compare_grid: ',trim(string),': ytn different',g1%ytn(nd),g2%ytn(nd)
endif

if (g1%nzp == g2%nzp) then
   nd=check_real(g1%ztn(1),g2%ztn(1),g1%nzp,g2%nzp)
   if ( nd > 0 ) then
      ict=ict+1
      print*,'Compare_grid: ',trim(string),': ztn different',g1%ztn(nd),g2%ztn(nd)
   endif
else
   print*,'Compare_grid: ',trim(string),': not checking: ztn, nzp different'
endif

if(g1%nxp == g2%nxp .and. g1%nyp == g2%nyp) then
   nd=check_real(g1%topo(1,1),g2%topo(1,1),g1%nxp*g1%nyp,g2%nxp*g2%nyp)
   if ( nd > 0 ) then
      ict=ict+1
      print*,'Compare_grid: ',trim(string),': topo different'
   endif
else
   print*,'Compare_grid: ',trim(string),': not checking: topo, nxp or nyp different'
endif

if (ict > 0) then
   print*,'Compare_grid: ',trim(string),': ',ict,' mismatches'
   ierr=1
endif

return
END SUBROUTINE compare_grid_def

END MODULE grid_struct
