!##############################################################################
Subroutine oda_read (n1,n2,n3,pi0)

use mem_oda
use mem_grid
use rconstants

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: pi0

integer :: iyears,imonths,idates,ihours,iyearf,imonthf,idatef,ihourf
integer :: ns,i,j,ng
real, dimension(:,:), allocatable :: prs

! Inventory all observation files. Due to possible mismatches with 
!   start/end time and obs file times, we are going to start 12 hours
!   before the actual start and end 12 hours after timmax

CALL date_add_to (iyear1,imonth1,idate1,itime1*100  &
                ,-43200.,'s',iyears,imonths,idates,ihours)
CALL date_add_to (iyear1,imonth1,idate1,itime1*100  &
                ,timmax+43200.,'s',iyearf,imonthf,idatef,ihourf)

CALL oda_file_inv (iyears,imonths,idates,ihours  &
                  ,iyearf,imonthf,idatef,ihourf)
   !print*,'++++++++++++++++: ',nupafiles

! First pass through the files: find number of unique station ids

CALL oda_sta_count (polelat,polelon)
   !print*,'++++++++++++++++: ',nupafiles

! Allocate obs structures

CALL oda_obs_alloc ()

! Fill data structures

CALL oda_sta_input (polelat,polelon,ngrids)

! We are assuming the files are consecutive in time,
! but if not, sort obs according to time

! Need a pressure to convert dewpoint (or whatever) to mixing ratio, as almost no
!   obs report mixing ratio. Problems though: many stations do not report pressure.
!   Tried to use model pressure at the time, but because of grids and
!   domain decomposition, don't have that for stations just outside grid/subdomain
!   boundaries. 
!
!  So only for this purpose, use a reference state pressure for the stations that
!  do not report, since the conversion is not overly sensitive to pressure.

allocate(prs(n2,n3))

ng=1
do j=1,n3
   do i=1,n2
      prs(i,j)=( (pi0(1,i,j)+pi0(2,i,j))*.5  &
                 *cpi) ** cpor * p00
  enddo
enddo

! Interpolate to station locations

do ns=1,num_oda_sfc
   CALL gdtost (prs, nnxp(ng), nnyp(ng)  &
               ,oda_sfc_info(ns)%xista(ng) &
               ,oda_sfc_info(ns)%xjsta(ng)  &
               ,oda_sfc_obs(ns)%psref)
   if (oda_sfc_obs(ns)%psref > 1e10) oda_sfc_obs(ns)%psref=-999.
enddo            

deallocate(prs)

return
END SUBROUTINE oda_read

!##############################################################################
Subroutine oda_obs_alloc ()

use mem_oda

implicit none

integer :: ns

allocate (oda_sfc_info(num_oda_sfc), oda_sfc_obs(num_oda_sfc))
allocate (oda_upa_info(num_oda_upa), oda_upa_obs(num_oda_upa))

!  Assuming one data time per ralph file for the allocations
!    This wasn't enough!!!!! Actual count of times have been done

do ns=1,num_oda_sfc
   allocate(oda_sfc_obs(ns)%time (maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%temp (maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%dewpt(maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%us(maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%vs(maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%u (maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%v (maxtimes_sfc))
   allocate(oda_sfc_obs(ns)%ps(maxtimes_sfc))
enddo

! Upper air
do ns=1,num_oda_upa
   allocate(oda_upa_obs(ns)%time(maxtimes_upa))
   allocate(oda_upa_obs(ns)%lp  (maxtimes_upa))
   allocate(oda_upa_obs(ns)%lz  (maxtimes_upa))
   allocate(oda_upa_obs(ns)%theta(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%rv(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%us(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%vs(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%zz(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%u (maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%v (maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%pi(maxupalevs,maxtimes_upa))
   allocate(oda_upa_obs(ns)%zgeo (maxupalevs,maxtimes_upa))
enddo

return
END SUBROUTINE oda_obs_alloc

!##############################################################################
Subroutine oda_file_inv (iyear1,imonth1,idate1,itime1  &
                        ,iyear2,imonth2,idate2,itime2)

use mem_oda

implicit none

integer :: iyear1,imonth1,idate1,itime1,iyear2,imonth2,idate2,itime2
integer :: nc,nf,lnf,nsfctot,nupatot
integer :: inyear,inmonth,indate,inhour
character(len=strl1) :: fnames(maxodafiles)
character(len=14) :: itotdate,itotdate_start,itotdate_end

!          Go through upper air and surface input files
!            and make inventory

!print*,'st:',iyear1,imonth1,idate1,itime1
!print*,'en:',iyear2,imonth2,idate2,itime2
CALL date_make_big (iyear1,imonth1,idate1,itime1,itotdate_start)
CALL date_make_big (iyear2,imonth2,idate2,itime2,itotdate_end)

if(oda_upapref(1:1) /= ' ' .and. oda_upapref(1:1) /= char(0) ) then

   nc=len_trim(oda_upapref)
   nupatot=-1
   CALL rams_filelist (fnames,oda_upapref(1:nc)//'????-??-??-????',nupatot)

   if(nupatot > maxodafiles) then
      print*,'too many oda upper air files'
      stop 'lots_of_oda_upper_air'
   endif

   nupafiles=0
   do nf=1,nupatot
      lnf=len_trim(fnames(nf))
   !   print*,lnf ,nupatot,fnames(nf)
   !   print*,lnf ,fnames(nf)(lnf-14:lnf)
      read(fnames(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour
      20 format(i4,1x,i2,1x,i2,1x,i4)

      CALL date_make_big (inyear,inmonth,indate,inhour*100,itotdate)
    ! print*, inyear,inmonth,indate,inhour
!print*,itotdate,itotdate_start,itotdate_end
      if(itotdate >= itotdate_start .and. itotdate <= itotdate_end) then
         nupafiles=nupafiles+1
         fnames_upa(nupafiles)=fnames(nf)
         itotdate_upa(nupafiles)=itotdate
      endif

   enddo

   CALL rams_dintsort (nupafiles,itotdate_upa,fnames_upa)

  ! do nf=1,nupafiles
  !    print*,'up files:',nf,itotdate_upa(nf),fnames_upa(nf)
  ! enddo

endif


if(oda_sfcpref(1:1) /= ' '.and. oda_sfcpref(1:1) /= char(0) ) then

   nc=len_trim(oda_sfcpref)
   nsfctot=-1
   CALL rams_filelist (fnames,oda_sfcpref(1:nc)//'????-??-??-????',nsfctot)

   if(nsfctot > maxodafiles) then
      print*,'too many oda surface air files'
      stop 'lots_of_oda_surface'
   endif

   nsfcfiles=0
   do nf=1,nsfctot
      lnf=len_trim(fnames(nf))
    !  print*,lnf ,nsfctot,fnames(nf)
    !  print*,lnf ,fnames(nf)(lnf-14:lnf)
      read(fnames(nf)(lnf-14:lnf),20) inyear,inmonth,indate,inhour

      CALL date_make_big (inyear,inmonth,indate,inhour*100,itotdate)

      if(itotdate >= itotdate_start .and. itotdate <= itotdate_end) then
         nsfcfiles=nsfcfiles+1
         fnames_sfc(nsfcfiles)=fnames(nf)
         itotdate_sfc(nsfcfiles)=itotdate
      endif

   enddo

   CALL rams_dintsort (nsfcfiles,itotdate_sfc,fnames_sfc)

  ! do nf=1,nsfcfiles
  !    print*,'sf files:',nf,itotdate_sfc(nf),fnames_sfc(nf)
  ! enddo

endif


!  start printing section
!--------------------------------------------------------------

print*,' '
print*,' '
print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  Obs 4DDA Input File Date Inventory -------------'
print*,'-------------------------------------------------------------'
print*,'---- Upper air   files:'
do nf=1,nupafiles
   print*,  itotdate_upa(nf),'   ',trim(fnames_upa(nf))
enddo
print*,'---- Surface obs files:'
do nf=1,nsfcfiles
   print*,  itotdate_sfc(nf),'   ',trim(fnames_sfc(nf))
enddo
print*,'------------------------------------------------------'

return
END SUBROUTINE oda_file_inv
