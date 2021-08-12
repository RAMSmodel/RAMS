!##############################################################################
Module gdf_input

! Surface variables

use grid_dims

implicit none

integer, parameter :: max_sfc_vars=50

type gdf_sfc_obs
   integer :: iqflags(max_sfc_vars,3),ihgtflg
   character(len=16) :: id
   real :: lat,lon,elev,hgt
   integer :: jyear,jmonth,jdate,jtime
   real :: sdata(max_sfc_vars)
end type

type(gdf_sfc_obs) :: rsfc_obs

! Upper air variables

integer, parameter :: max_up_levs=100,max_upa_vars=12

type gdf_upa_obs
   real, dimension(max_up_levs) :: p,t,z,r,zz,dz,fz
   character(len=16) :: id
   real :: lat,lon,elev
   integer :: lp,lz
   integer :: iqflagsp(max_up_levs,4,3),iqflagsz(max_up_levs,3,3)
   integer :: jyear,jmonth,jdate,jtime
end type

type(gdf_upa_obs) :: rupa_obs


! Header info

integer, parameter :: max_head_vars=10

type obs_header
   character(len=strl1) :: head_string(max_head_vars),head_units(max_head_vars) &
                       ,sfc_string(max_sfc_vars),sfc_units(max_sfc_vars)  &
                       ,therm_string(max_upa_vars),wind_string(max_upa_vars)
                        
   integer :: iun,iver,nhead,nvsfc,nvtherm,nvwind
end type

type(obs_header) :: header(1)

integer :: curfile


! Good station variables
integer, parameter :: maxgood_gdf=1000
integer :: ngood_gdf
character(len=64) :: good_id_gdf(maxgood_gdf)

! Bad station variables
integer, parameter :: maxbad_gdf=1000
integer :: nbad_gdf
character(len=64) :: bad_id_gdf(maxbad_gdf)


Contains

!##############################################################################
Subroutine gdf_good_obs ()

! Read a file of station names. If a station in the list is read,
!   use it. If a station is read that is not in the list, skip it.

! File name is "GDF_Good_Obs" and must exist in current directory.
! File format is one station ID (character string) per line.
! '#' character in column 1 designates a comment, blank lines are ignored
! Reading will stop upon file end or error. 

! If both good and bad files are present, the bad will override 
!    the good specification.

! A further assumption: if the good file exists, but is empty,
!   assume all stations are 

implicit none

character(len=64), save :: bname='GDF_Good_Obs'
character(len=64) :: station
logical :: there

inquire(file=bname, exist=there)

if (there) then

   open(188,file=bname, status='old')
   
   ngood_gdf=0
   do while(.true.)
   
      read(188,*,end=100,err=100) station
      if(station(1:1) == '#') cycle
      if(len_trim(station) == 0) cycle
      
      ngood_gdf=ngood_gdf+1
      good_id_gdf(ngood_gdf)=station
      
      print*,'Good obs setup:',ngood_gdf,trim(good_id_gdf(ngood_gdf))
      
      if(ngood_gdf >= maxgood_gdf) then
         print*,'Number of good stations in:',trim(bname)
         print*,'  exceeds maximum in gdf_good_obs. max=',maxgood_gdf
         stop 'gdf_good_obs: Too many good stations'
      endif
      
   enddo

else

   ngood_gdf=0
   
endif


100 continue
close(188)

return
END SUBROUTINE gdf_good_obs

!##############################################################################
Subroutine gdf_bad_obs ()

! Read a file of station names. If a station in the list is read,
!   skip it.
! File name is "GDF_Bad_Obs" and must exist in current directory.
! File format is one station ID (character string) per line.
! '#' character in column 1 designates a comment, blank lines are ignored
! Reading will stop upon file end or error. 

! If both good and bad files are present, the bad will override 
!    the good specification.

implicit none

character(len=64), save :: bname='GDF_Bad_Obs'
character(len=64) :: station
logical :: there

inquire(file=bname, exist=there)

if (there) then

   open(188,file=bname, status='old')
   
   nbad_gdf=0
   do while(.true.)
   
      read(188,*,end=100,err=100) station
      if(station(1:1) == '#') cycle
      if(len_trim(station) == 0) cycle

      nbad_gdf=nbad_gdf+1
      bad_id_gdf(nbad_gdf)=station
      
      print*,'Bad obs setup:',nbad_gdf,trim(bad_id_gdf(nbad_gdf))
      
      if(nbad_gdf >= maxbad_gdf) then
         print*,'Number of bad stations in:',trim(bname)
         print*,'  exceeds maximum in gdf_bad_obs. max=',maxbad_gdf
         stop 'gdf_bad_obs: Too many bad stations'
      endif
      
   enddo

else

   nbad_gdf=0
   
endif


100 continue
close(188)

return
END SUBROUTINE gdf_bad_obs

END MODULE gdf_input
