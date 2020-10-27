!##############################################################################
Subroutine input_rawi ()

use isan_coms
use gdf_input
use node_mod

implicit none

character(len=strl1) :: fnameu
real, dimension(maxlev) :: p,t,z,h,f,pp,zp,tp,hp,vz,uz,zz,rp,dz,fz
integer :: iqflagsp(maxlev,4,3),iqflagsz(maxlev,3,3)
character(len=16) :: idsta
character(len=14) :: obsdate
real(kind=8) :: nasecs,obssecs
integer :: k,ierr,ifile,nnn
integer :: jt,lp,lz,jyr,jmo,jdy,kl
real :: xlat,xlon,elev

CALL date_abs_secs (iproc_dates(natime),nasecs)
CALL date_unmake_big (iyear,imonth,idate,ihour,iproc_dates(natime))
ihour=ihour/100

print 559,inrawi
559 format(/'   Acquiring rawindsonde input file    ',A)
CALL rams_f_open (31,inrawi,'FORMATTED','OLD','READ',0)

! optionally output a file containing the lat lon locations
!    of each profile used for this analysis time

if (used_file(1:4) /= 'none') then 
   fnameu=trim(used_file)//'-upa-'  &
               //inrawi(len_trim(inrawi)-18:len_trim(inrawi))
   open(91,file=fnameu,status='unknown')
   rewind(91)
   write(91,'(a)') '999999  3 '
   write(91,'(a)') '3 HEADER  '
   write(91,'(a)') 'STATION_ID'
   write(91,'(a)') 'LATITUDE  '
   write(91,'(a)') 'LONGITUDE '
   write(91,'(a)') '0 DATA    '
endif

ifile=1
header(ifile)%iun=31
CALL gdf_read_upa_ver (ifile)


do while(.true.)
   
   CALL gdf_read_upa_obs (ifile,'yes',ierr)
   if(ierr == 1) exit

   do k=1,maxlev
      p(k)=0.  ;      t(k)=0.
      z(k)=0.  ;      h(k)=0.
      f(k)=0.  ;      pp(k)=0.
      zp(k)=0. ;      tp(k)=0.
      hp(k)=0. ;      vz(k)=0.
      uz(k)=0. ;      zz(k)=0.
      rp(k)=0. ;      dz(k)=0.
      fz(k)=0.
   enddo
   
   jyr   =rupa_obs%jyear
   jmo   =rupa_obs%jmonth
   jdy   =rupa_obs%jdate
   jt    =rupa_obs%jtime
   idsta =rupa_obs%id
   lp    =rupa_obs%lp
   lz    =rupa_obs%lz
   xlat  =rupa_obs%lat
   xlon  =rupa_obs%lon
   elev  =rupa_obs%elev

   if(lp > maxlev .or. lz > maxlev) then
      print*,'Rawindsonde read error!'
      print*,'  Number of input levels greater than MAXLEV'
      print*,'  lp,lz,maxlev:',lp,lz,maxlev
      stop 'rawin-maxlev'
   endif

   CALL gdf_upa_get_profile (p(1),nnn,'p_pas',  'p')
   CALL gdf_upa_get_profile (z(1),nnn,'geo',    'p')
   CALL gdf_upa_get_profile (t(1),nnn,'tempk',  'p')
   CALL gdf_upa_get_profile (h(1),nnn,'relfrac','p')
   
   iqflagsp(1:lp,1:4,1:3)=rupa_obs%iqflagsp(1:lp,1:4,1:3)

   CALL gdf_upa_get_profile (zz(1),nnn,'zz',       'z')
   CALL gdf_upa_get_profile (fz(1),nnn,'speed',    'z')
   CALL gdf_upa_get_profile (dz(1),nnn,'direction','z')

   iqflagsz(1:lz,1:3,1:3)=rupa_obs%iqflagsz(1:lz,1:3,1:3)
   
   ! Check date
   CALL date_make_big (jyr,jmo,jdy,jt*100,obsdate)
   CALL date_abs_secs (obsdate,obssecs)
   !print*,'===> ',natime,iobswin
   !print*,'===> ',iproc_dates(natime),' ',nasecs ,iyear,imonth,idate,ihour
   !print*,'===> ',obsdate            ,' ',obssecs,jyr  ,jmo   ,jdy  ,jt*100
   if((iobswin < 0 .and.  &
      (obssecs < nasecs+iobswin .or. obssecs > nasecs)).or.  &
      (iobswin > 0 .and.  &
      (obssecs < nasecs-iobswin .or. obssecs > nasecs+iobswin)).or.  &
      (iobswin .eq. 0 .and. nasecs .ne. obssecs)) then
      print*,'station ',idsta(1:len_trim(idsta))  &
            ,' excluded - out of obs window (s)',iobswin
      cycle
   !else
   !   print*,'station ',idsta(1:len_trim(idsta))  &
   !         ,' within obs window (s)',iobswin
   endif

   ! ---------------------------------------------------------
   ! Convert input variables to SI units and appropriate form
   !    Set missing values to 1E30.
   !      PP - Pressure (Pascals)
   !      TP - Temperature (Kelvin)
   !      RP - Relative humidity (fraction)
   !      ZP - Height of pressure levels (meters)
   !      ZZ - Height of wind levels (meters)
   !      DZ - Wind direction (degrees)
   !      FZ - Wind speed (m/s)
   !      XLAT - Station latitude
   !      XLON - Station longitude
   !      ELEV - Station elevation (meters)
   !      LP   - Number of pressure levels
   !      LZ   - Number of wind levels in height
   !      IDSTA- Station ID
   ! ---------------------------------------------------------

   do kl=1,lp
      if (p(kl) > -998.) then
         pp(kl)=p(kl)
      else
         pp(kl)=1.e30
      endif
     
      if (z(kl) > -998.) then
         zp(kl)=z(kl)
      else
         zp(kl)=1.e30
      endif
     
      if (t(kl) > -998.) then
         tp(kl)=t(kl)
      else
         tp(kl)=1.e30
      endif
     
      if (h(kl) > -998.) then
         rp(kl)=h(kl)
      else
         rp(kl)=1.e30
      endif
      
      !MJW Make sure that either pressure or height is non-missing
      if(pp(kl) > 1.e20 .and. zp(kl) > 1.e20) then
         tp(kl)=1.e30 ; rp(kl)=1.e30
      endif

   enddo

   do kl=1,lz
      if (zz(kl) > -998.) then
         zz(kl)=zz(kl)
      else
         zz(kl)=1.e30
      endif

      if (fz(kl) > -998.) then
         fz(kl)=fz(kl)
      else
         fz(kl)=1.e30
      endif

      if (dz(kl) > -998.) then
         dz(kl)=dz(kl)
      else
         dz(kl)=1.e30
      endif

      !MJW Make sure that height is non-missing
      if(zz(kl) > 1.e20) then
         fz(kl)=1.e30  ;  dz(kl)=1.e30
      endif
      
   enddo
   
   !----------------------------------

   CALL sndproc (lp,lz,xlat,xlon,elev,idsta  &
                ,tp,zp,pp,rp,dz,fz,zz)

enddo

print*,'End of rawindsonde file:',inrawi
close(31)

if (used_file(1:4) /= 'none') close(91) 

return
END SUBROUTINE input_rawi

!##############################################################################
Subroutine sndproc (lp,lz,xlat,xlon,elev,idsta,tp,zp,pp,rp,dz,fz,zz)

use isan_coms
use rconstants

implicit none
                   
real, dimension(*) :: tp,zp,pp,rp,dz,fz,zz
integer :: lp,lz
real :: xlat,xlon,elev
character(len=*) :: idsta

real, dimension(maxlev) :: uz,vz,thz

integer :: not,llp,k,llz,kk,lbc
real :: bchyd,rss,pio,tho,zso
real, external :: rsatmix

do not=1,notsta
   if('r'//idsta(1:5).eq.notid(not)(1:6)) then
      print 31,idsta,xlat,xlon,elev,0,0
      return
   endif
enddo
31 format(' ----- Station omitted   ',1X,A8,' Lat:',F8.2,'  Lon:',F8.2  &
         ,'  Elev:',F10.2,'  P Levels:',I5,'  Wind levels:',I5)

nsta=nsta+1

!--------------------------------------------------------------------
!mjb added the following from 4a as the following fix does not handle
!    pressure levels flagged as suspect by qc (501) 

! Throw out any levels that have suspect pressure
! ALSO, don't read in any levels above the gridded data, since that leads
!  to bulls-eyes around all soundings near model top
      
llp=0
do k=1,lp
   if(pp(k).lt.1.e20) then
      if(pp(k) < levpr(nprz)*100.)then
         print 35,pp(k)/100.,levpr(nprz)
35       format('   discarding ',f6.1, &
                'mb level since it is above gridded data top at ',i4,'mb')
         cycle
      endif
      llp=llp+1
      tp(llp)=tp(k)
      zp(llp)=zp(k)
      pp(llp)=pp(k)
      rp(llp)=rp(k)
   endif
enddo
lp=llp

! Throw out any levels that have suspect winds
llz=0
do k=1,lz
   if(zz(k).lt.1.e20.and.dz(k).lt.1.e20.and.fz(k).lt.1.e20) then
      llz=llz+1
      zz(llz)=zz(K)
      fz(llz)=fz(K)
      dz(llz)=dz(K)
      !print*,'mjb w1',llz,zz(llz),fz(llz),dz(llz)
   endif
enddo
!print*,'mjb - wind levels reduced: ',lz,llz
lz=llz
!--------------------------------------------------------------------

! arrange winds on height levels
do k=1,lz
   if(dz(k).lt.1e19.and.fz(k).lt.1e19) CALL winduv (dz(k),fz(k),uz(k),vz(k))
enddo
CALL sort3 (zz,uz,vz,lz)

up_lat(nsta)=xlat
up_lon(nsta)=xlon
up_top(nsta)=elev
up_lp(nsta)=lp
up_lz(nsta)=lz
up_chstid(nsta)=idsta(1:8)

! fix missing values if possible
do k=1,lp
   ! missing temperatures are:
   ! set to valid upper temperature if none below,
   ! interpolated linearly in logP if between valid values,
   ! left missing if no valid temperatures are above
   if(tp(k).gt.1e29.and.k.ne.lp) then
      !if(k > 1)then
      !   print*,'missing temp',k,lp,tp(k-1)
      !else
      !   print*,'missing temp',k,lp,tp(k)
      !endif
      kk=k+1
      do while (kk.lt.lp.and.tp(kk).gt.1e29)
         kk=kk+1
      enddo
      if (k.eq.1.and.tp(kk).lt.1e30) then
         tp(k)=tp(kk)
      elseif(tp(kk).lt.1e30) then
         tp(k)=tp(k-1)+(tp(kk)-tp(k-1))*log(pp(k)/pp(k-1))/log(pp(kk)/pp(k-1))
      endif
   endif
   ! same for humidities
   if(rp(k).gt.1e29.and.k.ne.lp) then
      !if(k > 1)then
      !   print*,'missing rh',k,lp,rp(k-1)
      !else
      !   print*,'missing rh',k,lp,rp(k)
      !endif
      kk=k+1
      do while (kk.lt.lp.and.rp(kk).gt.1e29)
         kk=kk+1
      enddo
      if (k.eq.1.and.rp(kk).lt.1e30) then
         rp(k)=rp(kk)
      elseif(rp(kk).lt.1e30) then
         rp(k)=rp(k-1)+(rp(kk)-rp(k-1))*log(pp(k)/pp(k-1))/log(pp(kk)/pp(k-1))
      endif
   endif

enddo

! Now that we filled in as many temps and rh's we could, we will recompute heights.
!        Find a boundary condition as first level at or below 500 mb. If we do 
!        not find a good level, we don't use the sounding.
!
!  If there is no thermodynamic info (maybe only winds), skip this section.

if (lp > 0) then

   bchyd=50000.
   lbc=0
   do k=lp,1,-1
      if(pp(k) < 1.e19 .and. pp(k) >= bchyd .and.  &
         tp(k) < 1.e19 .and. zp(k) < 1.e19) then
         lbc=k
         go to 52
      endif
   enddo

   ! Didn't find a good level, so let's look up from 500mb...
   do k=1,lp
      if(pp(k) < 1.e19 .and. pp(k) < bchyd .and.  &
         tp(k) < 1.e19 .and. zp(k) < 1.e19) then
         lbc=k
         go to 52
      endif
   enddo
   nsta=nsta-1
   print*,' sndproc: Could not find good rawindsonde z boundary level.'
   print*,'    Discarding sounding:',idsta(1:8)
   return

   52 continue

   ! Compute theta - if rh non-missing, compute vitual theta.

   do k=1,lp
      thz(k)=1.e30
      if(tp(k) < 1.e19 .and. pp(k) < 1.e19) then
         thz(k)=tp(k)*(p00/pp(k))**rocp
         if(rp(k) < 1.e19) then
            rss=rsatmix(pp(k),tp(k))*rp(k)
            thz(k)=thz(k)*(1.+.61*rss)
         endif
      endif
   enddo

   ! Recompute heights

   pio=cp*(pp(lbc)/p00)**rocp
   zso=zp(lbc)
   tho=thz(lbc)
   do k=lbc+1,lp
      zp(k)=1.e30
      if(pp(k) < 1e19 .and. thz(k) < 1e19)then
         zp(k)=zso+.5*(thz(k)+tho)*(pio-cp*(pp(k)/p00)**rocp)/g
         zso=zp(k)
         pio=cp*(pp(k)/p00)**rocp
         tho=thz(k)
      endif
   enddo

   zso=zp(lbc)
   pio=cp*(pp(lbc)/p00)**rocp
   tho=thz(lbc)
   do k=lbc-1,1,-1
      zp(k)=1.e30
      if(pp(k) < 1e19 .and. thz(k) < 1e19)then
         zp(k)=zso+.5*(thz(k)+tho)*(pio-cp*(pp(k)/p00)**rocp)/g
         zso=zp(k)
         pio=cp*(pp(k)/p00)**rocp
         tho=thz(k)
      endif
   enddo

endif

! Do one last filter and fill final arrays. 
!    If any heights are missing, discard the level.

llp=0
do k=1,lp
   if(zp(k) < 1.e19) then
      llp=llp+1
      up_t(nsta,llp)=tp(k)
      up_z(nsta,llp)=zp(k)
      up_p(nsta,llp)=pp(k)
      up_r(nsta,llp)=rp(k)
   endif
enddo

up_lp(nsta)=llp

do k=1,lz
   up_uz(nsta,k)=uz(k)
   up_vz(nsta,k)=vz(k)
   up_zz(nsta,k)=zz(k)
enddo

print 1001,nsta,idsta,up_lat(nsta),up_lon(nsta),up_top(nsta),lp,lz
1001 format(' Rawindsonde station found ',I5,' ID:',1X,A8,' Lat:'  &
     ,F8.2,'  Lon:', F8.2,'  Elev:',0PF10.2,'  P levels:',I5  &
     ,'  Wind levels:',I5)


if (used_file(1:4) /= 'none') then
   ! write to file containing the lat lon locations of each profile
   write(91,'(a8,2f12.4)') idsta,up_lat(nsta),up_lon(nsta)
endif
   

if(nsta.ge.maxsta) then
   print 92,maxsta
   92 format(' Number of stations found is greater than MAXSTA',I6)
   stop 'st3-maxsta'
endif

return
END SUBROUTINE sndproc

!##############################################################################
Subroutine input_sfc ()

use isan_coms
use gdf_input
use node_mod

implicit none

character(len=16) :: idsta
character(len=strl1) :: fnameu
integer :: iqflags(max_sfc_vars,3)
real*8 nasecs

integer, parameter :: num_convars=5
character(len=16) :: varn(num_convars)

real :: vars(num_convars)
data varn/'direction','speed','tempc','dewptc','press_pa'/

integer :: ifile,ierr
integer :: nv,jyr,jmo,jdy,jt

real :: xlat,xlon,zx,ffx,ddx,tx,tdx,px
   
CALL date_abs_secs (iproc_dates(natime),nasecs)
CALL date_unmake_big (iyear,imonth,idate,ihour,iproc_dates(natime))
ihour=ihour/100

if(used_file(1:4) /= 'none') then
   ! output a file containing the lat lon locations of each sfc obs
   fnameu=trim(used_file)//'-sfc-'  &
      //insrfce(len_trim(insrfce)-18:len_trim(insrfce))
   open(92,file=fnameu,status='unknown')
   write(92,'(a)') '999999  3 '
   write(92,'(a)') '3 HEADER  '
   write(92,'(a)') 'STATION_ID'
   write(92,'(a)') 'LATITUDE  '
   write(92,'(a)') 'LONGITUDE '
   write(92,'(a)') '0 DATA    '
endif

! Read surface observations from a file

print 559,insrfce
559 format(/'   Acquiring surface input file    ',A40)
CALL rams_f_open (31,insrfce,'FORMATTED','OLD','READ',0)

ifile=1
header(ifile)%iun=31
CALL gdf_read_sfc_ver (ifile)


do while (.true.)
   
   ! Read a surface obs

   CALL gdf_read_sfc_obs (ifile,'yes',ierr)
   if(ierr==1) exit
   
   jyr=rsfc_obs%jyear
   jmo=rsfc_obs%jmonth
   jdy=rsfc_obs%jdate
   jt=rsfc_obs%jtime
   idsta=rsfc_obs%id
   xlat=rsfc_obs%lat
   xlon=rsfc_obs%lon
   zx=rsfc_obs%elev

   CALL gdf_sfc_data_convert (vars,varn,num_convars)
            
   ffx=  vars(2)
   ddx=  vars(1)
   tx=   vars(3)
   tdx=  vars(4)
   px=   vars(5)

   do nv=1,max_sfc_vars
      iqflags(nv,1:3)=rsfc_obs%iqflags(nv,1:3)
   enddo
 
   !-------------------------------------------------------------------
   !         Section for conversion to SI  units
   !-------------------------------------------------------------------
   !   At the end of this section, the following variables are
   !   required to be in SI  units.
   !
   !     XLAT   - STATION LATITUDE
   !     XLON   - STATION LONGITUDE
   !       ZX   - STATION ELEVATION
   !      DDX   - WIND DIRECTION
   !      FFX   - WIND SPEED
   !       TX   - TEMPERATURE IN KELVIN
   !      TDX   - DEWPOINT TEMPERATURE IN KELVIN
   !       PX   - SURFACE PRESSURE
   !    IDSTA   - Station ID (char*5)
   !
   !              Any missing values must be set to 1.E30
   !--------------------------------------------------------------------

   if(zx < -998.) cycle

   if (ffx > -998.) then
      ffx=ffx
   else
      ffx=1.e30
   endif

   if (ddx > -998.) then
      ddx=ddx
   else
      ddx=1.e30
   endif

   if (tx > -998.) then
      tx=tx+273.15
   else
      tx=1.e30
   endif

   if (tdx > -998.) then
      tdx=tdx+273.15
   else
      tdx=1.e30
   endif

   if (px > -998.) then
      px=px
   else
      px=1.e30
   endif

   CALL sfcproc (nasecs,jyr,jmo,jdy,jt,idsta,xlat,xlon  &
                ,zx,ddx,ffx,tx,tdx,px)

enddo

print*,'End of surface file:',insrfce
close(31)

if (used_file(1:4) /= 'none') close(92) 

return
END SUBROUTINE input_sfc

!##############################################################################
Subroutine sfcproc (nasecs,jyr,jmo,jdy,jt,idsta,xlat,xlon  &
                   ,zx,ddx,ffx,tx,tdx,px)
                   
use isan_coms
use rconstants

implicit none

integer :: jyr,jmo,jdy,jt
real :: xlat,xlon,zx,ddx,ffx,tx,tdx,px
real*8 nasecs
character(len=*) :: idsta

real :: t1000,t900,t850,t700,t600,t500,p1000,p900,p850,p700,p600,p500  &
       ,z1000,z900,z850,z700,z600,z500
real*8 obssecs,obssecsns
character(len=14) :: obsdate
character(len=8) :: csfc

integer :: nsu,ns,not,mis1,mis2,irepl,isfc
real, external :: rsatmix

CALL date_make_big (jyr,jmo,jdy,jt*100,obsdate)
CALL date_abs_secs (obsdate,obssecs)

isfc=0
nssfc=nssfc+1
nsu=nssfc
sf_chstid(nsu)=idsta
sf_top(nsu)=zx
sf_lat(nsu)=xlat
sf_lon(nsu)=xlon
if(ddx.lt.1.e19.and.ffx.lt.1.e19) then
   CALL winduv (ddx,ffx,sf_u(nsu),sf_v(nsu))
else
   sf_u(nsu)=1.e30
   sf_v(nsu)=1.e30
endif
sf_date(nsu)=obsdate

! check whether to omit station

do not=1,notsta
   if('s'//idsta(1:5).eq.notid(not)(1:6)) then
      isfc=5
      nssfc=nssfc-1
      goto 110
   endif
enddo

! Check date against iobswin

if((iobswin < 0 .and.  &
   (obssecs .lt. nasecs+iobswin .or. obssecs .gt. nasecs)).or.  &
   (iobswin > 0 .and.  &
   (obssecs .lt. nasecs-iobswin .or. obssecs .gt. nasecs+iobswin)).or.  &
   (iobswin .eq. 0 .and. nasecs .ne. obssecs)) then
   isfc=3
   nssfc=nssfc-1
   goto 110
endif

t1000=287.6
t900=281.9
t850=278.6
t700=268.6
t600=260.8
t500=252.0
p1000=100000.
p900=90000.
p850=85000.
p700=70000.
p600=60000.
p500=50000.
z1000=111.
z900=988.
z850=1457.
z700=3012.
z600=4206.
z500=5574.
if(px.gt.1.e20.and.tx.lt.1.e20) then
   if(2.*zx.lt.z1000+z900) then
      px=p1000*exp(-g*(zx-z1000)/(rgas*.5*(tx+t1000)))
   elseif(2.*zx.lt.z900+z850) then
      px=p900*exp(-g*(zx-z900)/(rgas*.5*(tx+t900)))
   elseif(2.*zx.lt.z850+z700) then
      px=p850*exp(-g*(zx-z850)/(rgas*.5*(tx+t850)))
   elseif(2.*zx.lt.z700+z600) then
      px=p700*exp(-g*(zx-z700)/(rgas*.5*(tx+t700)))
   elseif(2.*zx.lt.z600+z500) then
      px=p600*exp(-g*(zx-z600)/(rgas*.5*(tx+t600)))
   endif
endif

if(tx.lt.1.e19.and.px.lt.1.e19) then
   if(tdx.lt.1.e19) then
      sf_r(nsu)=rsatmix(px,tdx)/rsatmix(px,tx)
   else
      sf_r(nsu)=1.e30
   endif
   if(zx.lt.1.e19) then
      sf_s(nsu)=cp*tx+g*zx
   else
      sf_s(nsu)=1.e30
   endif
   sf_t(nsu)=tx*(p00/px)**rocp
   sf_p(nsu)=px
else
   sf_t(nsu)=1.e30
   sf_p(nsu)=1.e30
   sf_r(nsu)=1.e30
   sf_s(nsu)=1.e30
endif

do ns=1,nsu-1

   ! determine if current date is closer for idsta
   
   if(sf_chstid(ns)==idsta) then
      CALL date_abs_secs (sf_date(ns),obssecsns)
      if(abs(nasecs-obssecsns) > abs(nasecs-obssecs)) then
         isfc=7
         nssfc=nssfc-1
         mis1=0
         if(sf_u(ns).gt.1.e19) mis1=mis1+1
         if(sf_v(ns).gt.1.e19) mis1=mis1+1
         if(sf_t(ns).gt.1.e19) mis1=mis1+1
         if(sf_p(ns).gt.1.e19) mis1=mis1+1
         if(sf_r(ns).gt.1.e19) mis1=mis1+1
         if(sf_top(ns).gt.1.e19) mis1=mis1+1
         mis2=0
         if(sf_u(nsu).gt.1.e19) mis2=mis2+1
         if(sf_v(nsu).gt.1.e19) mis2=mis2+1
         if(sf_t(nsu).gt.1.e19) mis2=mis2+1
         if(sf_p(nsu).gt.1.e19) mis2=mis2+1
         if(sf_r(nsu).gt.1.e19) mis2=mis2+1
         if(sf_top(nsu).gt.1.e19) mis2=mis2+1
         if(mis2.ge.mis1) then
            sf_u(ns)=sf_u(nsu)
            sf_v(ns)=sf_v(nsu)
            sf_t(ns)=sf_t(nsu)
            sf_s(ns)=sf_s(nsu)
            sf_p(ns)=sf_p(nsu)
            sf_top(ns)=sf_top(nsu)
            sf_r(ns)=sf_r(nsu)
            sf_lat(ns)=sf_lat(nsu)
            sf_lon(ns)=sf_lon(nsu)
            csfc=sf_chstid(ns)
            sf_chstid(ns)=sf_chstid(nsu)
            sf_date(ns)=sf_date(nsu)
            isfc=4
            irepl=ns
            goto 110
         endif
      endif
   endif

enddo

do ns=1,nsu-1

   ! check for redundant surface obs

   if(abs(xlon-sf_lon(ns))<stasep.and.abs(xlat-sf_lat(ns))<stasep) then
      isfc=2
      nssfc=nssfc-1
      csfc=sf_chstid(ns)
      irepl=ns
      
      mis1=0
      if(sf_u(ns).gt.1.e19) mis1=mis1+1
      if(sf_v(ns).gt.1.e19) mis1=mis1+1
      if(sf_t(ns).gt.1.e19) mis1=mis1+1
      if(sf_p(ns).gt.1.e19) mis1=mis1+1
      if(sf_r(ns).gt.1.e19) mis1=mis1+1
      if(sf_top(ns).gt.1.e19) mis1=mis1+1
      mis2=0
      if(sf_u(nsu).gt.1.e19) mis2=mis2+1
      if(sf_v(nsu).gt.1.e19) mis2=mis2+1
      if(sf_t(nsu).gt.1.e19) mis2=mis2+1
      if(sf_p(nsu).gt.1.e19) mis2=mis2+1
      if(sf_r(nsu).gt.1.e19) mis2=mis2+1
      if(sf_top(nsu).gt.1.e19) mis2=mis2+1
      if(mis2.lt.mis1) then
         sf_u(ns)=sf_u(nsu)
         sf_v(ns)=sf_v(nsu)
         sf_t(ns)=sf_t(nsu)
         sf_s(ns)=sf_s(nsu)
         sf_p(ns)=sf_p(nsu)
         sf_top(ns)=sf_top(nsu)
         sf_r(ns)=sf_r(nsu)
         sf_lat(ns)=sf_lat(nsu)
         sf_lon(ns)=sf_lon(nsu)
         sf_chstid(ns)=sf_chstid(nsu)
         sf_date(ns)=sf_date(nsu)
         isfc=1
      endif
   endif

enddo

if(nssfc.ge.maxsfc.or.nssfc.gt.maxsname) then
   print 92,maxsfc
   92 format(' Number of sfc obs found greater','than MAXSFC or maxsname ',I6)
   stop 'st3-maxsfc'
endif

110 continue

if(isfc==0) print 100,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu)
100 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2)
   
if(isfc==1) print 101,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl        
101 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' miss - repl ',A8,1X,I5)
     
if(isfc==2) print 102,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl    
102 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' redundant w/',A8,1X,I5)
 
if(isfc==3) print 103,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),iobswin          
103 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' out of obswin',I5)
    
if(isfc==4) print 104,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl       
104 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' newer- repl ',A8,1X,I5)
     
if(isfc==5) print 105,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu)     
105 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' ommitted in namelist')
           
if(isfc==6) print 106,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl   
106 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' redundant w/',A8,1X,I5)
    
if(isfc==7) print 107,nssfc,idsta,sf_lat(nsu),sf_lon(nsu),sf_top(nsu),csfc,irepl       
107 format('  Sfc obs ',I5,' - ',A8,2F8.2,F10.2  &
           ,' newer - miss ',A8,1X,I5)

if(isfc == 0 .and. used_file(1:4) /= 'none') then
   ! write to file containing the lat lon locations of each used sfc obs
   write(92,'(a8,2f12.4)') idsta,sf_lat(nsu),sf_lon(nsu)
endif

return
END SUBROUTINE sfcproc
