!##############################################################################
Module rcommons

use grid_dims

implicit none

!**************************************************
! For rainit.f90 routines
!**************************************************
character(len=strl1), save :: fnames(maxfiles)
integer, save :: nfgrids(maxfiles),ifdates(maxfiles),iftimes(maxfiles)
real, save :: ftimes(maxfiles),startutc
!**************************************************
! FROM REVUIN
!**************************************************
character(len=strl1) :: anpref,revpref
character(len=8) :: anatype
character(len=20) :: xvar,yvar,zvar,tvar
character(len=strl1) :: revuvar(maxrevu)
integer :: igrid,iztran
!**************************************************

character(len=1) :: ftran
integer,save :: maxmem
integer :: ierr_getvar,ifound,ivar_type,iany,iptrans_lim

integer :: nii,nib,nie,niinc,njj,njb,nje,njinc,nnb,nne,nninc &
          ,itbeg,itstep,itend,ixbeg,ixstep,ixend  &
          ,iybeg,iystep,iyend,izbeg,izstep,izend

integer, dimension(nxpmax,maxgrds) :: ipm,jpm,nrz,kpm
integer, parameter :: nplevs=41
integer, save :: iplevs(nplevs)

data iplevs/1050, 1025, 1000, &
       975,  950,  925,  900, &
       875,  850,  825,  800, &
       775,  750,  725,  700, &
       675,  650,  625,  600, &
       575,  550,  525,  500, &
       475,  450,  425,  400, &
       375,  350,  325,  300, &
       275,  250,  225,  200, &
       175,  150,  125,  100, &
        75,   50/

END MODULE rcommons
