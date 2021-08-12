!##############################################################################
Module rrad3

implicit none

!---------------------------------------------------------------------------
integer                         :: nrad,narad,nsolb,nb,ng,jday
integer, parameter              :: mb=8,mg=3,mk=7,ncog=5,ncb=2,npartob=13  &
                                  ,npartg=7,namax=10
integer, dimension(mg,mb)       :: npsb
real                            :: solfac
real, dimension(mb)             :: nuum,ralcs,solar1,solar0,a0,a1,a2,a3  &
                                  ,wlenlo,wlenhi
real, dimension(150)            :: exptabc
real, dimension(mg,mb)          :: prf,trf,ulim
real, dimension(mg,mk,mb)       :: wght,xp,alpha,beta
real, dimension(ncog,mb,npartob) :: ocoef
real, dimension(ncb,mb,npartob) :: bcoef
real, dimension(ncog,mb,npartg) :: gcoef
!---------------------------------------------------------------------------

END MODULE rrad3
