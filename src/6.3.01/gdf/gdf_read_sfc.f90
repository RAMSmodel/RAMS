!##############################################################################
Subroutine gdf_read_sfc_ver (ifile)

! Reads Ralph sfc file version and header

use gdf_input

implicit none

integer :: ifile

integer :: nh
character(len=strl1) :: line,tokens(100)
integer :: ntok

header(ifile)%head_string(1:max_head_vars)=''
header(ifile)%sfc_string(1:max_sfc_vars)=''

read(header(ifile)%iun,'(a)') line

CALL parse (line,tokens,ntok)

if(trim(tokens(1)) == '999999') then
   read(tokens(2),*) header(ifile)%iver
else
   header(ifile)%iver=1
endif

! Partial implementation of new v3

if(header(ifile)%iver==3) then

   ! Read header info
   read(header(ifile)%iun,'(a)') line
   CALL parse (line,tokens,ntok)
   read(tokens(1),*) header(ifile)%nhead
   do nh=1,header(ifile)%nhead
      read(header(ifile)%iun,'(a)') line
      tokens=' ' ; CALL parse (line,tokens,ntok)
      header(ifile)%head_string(nh) = tokens(1)
      header(ifile)%head_units(nh) = ' '
      if(ntok > 1)  header(ifile)%head_units(nh) = tokens(2)
   enddo
endif

! Read variable labels and units
read(header(ifile)%iun,'(a)') line
CALL parse (line,tokens,ntok)
read(tokens(1),*) header(ifile)%nvsfc
do nh=1,header(ifile)%nvsfc
   read(header(ifile)%iun,'(a)') line
   tokens=' ' ; CALL parse (line,tokens,ntok)
   header(ifile)%sfc_string(nh) = tokens(1)
   header(ifile)%sfc_units(nh) = ' '
   if(ntok > 1)  header(ifile)%sfc_units(nh) = tokens(2)
enddo

if(header(ifile)%iver == 1) read(header(ifile)%iun,*)

! See if there is a "good" station list present
CALL gdf_good_obs ()

! See if there is a "bad" station list present
CALL gdf_bad_obs ()

!print*,'done read header'

return
END SUBROUTINE gdf_read_sfc_ver

!##############################################################################
Subroutine gdf_read_sfc_obs (ifile,qcheck,ierr)

! Reads one surface obs from gdf sfc file

use gdf_input

implicit none

integer :: ifile,ierr
character(len=*) :: qcheck
integer :: ntok,iqfl(max_sfc_vars)
character(len=strl1) :: var_string
character(len=16) :: cflags,tokens(100)
character(len=strl2) :: line
integer :: jd,nvar,nv,ic,k,iq,iqf,nh,nt
integer :: nbad,ngood

100 continue

ierr=0
nvar=5
tokens=''

! Set file id to be used for some subsequent calls
curfile = ifile

read(header(ifile)%iun,'(a)',end=20,err=20) line
if(len_trim(line)==0) goto 10

! Check to see if this is a header line. If so, (re)read header.

!!!!!! Not working for v3 yet.

if(line(1:6) == '999999') then
   read(line(7:),*) header(ifile)%iver
   read(header(ifile)%iun,*) nvar
   do nv=1,nvar
      read(header(ifile)%iun,*) var_string
   enddo
   read(header(ifile)%iun,'(a)',end=20,err=20) line
endif

CALL parse (line,tokens,ntok)
   
if(header(ifile)%iver==1) then

   read(tokens(1),*,err=1,end=1) jd
   read(tokens(2),*,err=1,end=1) rsfc_obs%jtime
   read(tokens(3),'(a)',err=1,end=1) rsfc_obs%id
   read(tokens(4),*,err=1,end=1) rsfc_obs%lat
   read(tokens(5),*,err=1,end=1) rsfc_obs%lon
   read(tokens(6),*,err=1,end=1) rsfc_obs%elev
   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   read(tokens(7),*,err=1,end=1) rsfc_obs%sdata(1)
   read(tokens(8),*,err=1,end=1) rsfc_obs%sdata(2)
   read(tokens(9),*,err=1,end=1) rsfc_obs%sdata(3)
   read(tokens(10),*,err=1,end=1) rsfc_obs%sdata(4)
   read(tokens(11),*,err=1,end=1) rsfc_obs%sdata(5)
   read(tokens(12),'(a)',err=1,end=1) cflags
   1 continue

   ic=1
   do nv=1,header(ifile)%nvsfc
      read(cflags(ic:),'(1x,3i1)') (rsfc_obs%iqflags(nv,k),k=1,3)
      ic=ic+4
   enddo
   
   rsfc_obs%jyear=jd/10000
   rsfc_obs%jmonth=mod(jd,10000)/100
   rsfc_obs%jdate=mod(jd,100)

elseif(header(ifile)%iver==2) then
   read(tokens(1),*,err=2,end=2) rsfc_obs%jyear
   read(tokens(2),*,err=2,end=2) rsfc_obs%jmonth
   read(tokens(3),*,err=2,end=2) rsfc_obs%jdate
   read(tokens(4),*,err=2,end=2) rsfc_obs%jtime
   read(tokens(5),'(a)',err=2,end=2) rsfc_obs%id
   read(tokens(6),*,err=2,end=2) rsfc_obs%lat
   read(tokens(7),*,err=2,end=2) rsfc_obs%lon
   read(tokens(8),*,err=2,end=2) rsfc_obs%elev
   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   
   nt=9
   do nv=1,header(ifile)%nvsfc
      read(tokens(nt  ),*,err=2,end=2) rsfc_obs%sdata(nv)
      read(tokens(nt+1),*,err=2,end=2) iqfl(nv)
      nt = nt + 2
   enddo
   
   2 continue
   
   do nv=1,header(ifile)%nvsfc
      rsfc_obs%iqflags(nv,1)=iqfl(nv)/100
      rsfc_obs%iqflags(nv,2)=mod(iqfl(nv)/10,10)
      rsfc_obs%iqflags(nv,3)=mod(iqfl(nv),10)
   enddo

elseif(header(ifile)%iver==3) then

   rsfc_obs%hgt=0.
   rsfc_obs%ihgtflg=1
   nt=1
   do nh=1,header(ifile)%nhead
      if(header(ifile)%head_string(nh)=='YEAR') then
         read(tokens(nt),*,err=3) rsfc_obs%jyear
      elseif(header(ifile)%head_string(nh)=='MONTH') then
         read(tokens(nt),*,err=3) rsfc_obs%jmonth
      elseif(header(ifile)%head_string(nh)=='DAY') then
         read(tokens(nt),*,err=3) rsfc_obs%jdate
      elseif(header(ifile)%head_string(nh)=='HOUR') then
         read(tokens(nt),*,err=3) rsfc_obs%jtime
      elseif(header(ifile)%head_string(nh)=='STATION_ID') then
         read(tokens(nt),'(a)',err=3) rsfc_obs%id
      elseif(header(ifile)%head_string(nh)=='LATITUDE') then
         read(tokens(nt),*,err=3) rsfc_obs%lat
      elseif(header(ifile)%head_string(nh)=='LONGITUDE') then
         read(tokens(nt),*,err=3) rsfc_obs%lon
      elseif(header(ifile)%head_string(nh)=='ELEVATION') then
         read(tokens(nt),*,err=3) rsfc_obs%elev
      elseif(header(ifile)%head_string(nh)=='HEIGHT') then
         read(tokens(nt),*,err=3) rsfc_obs%hgt
      elseif(header(ifile)%head_string(nh)=='HEIGHT_FLAG') then
         read(tokens(nt),*,err=3) rsfc_obs%ihgtflg
      else
         print*,'unknown header record:', trim(header(ifile)%head_string(nh))
      endif
      nt=nt+1
   enddo
   
   do nv=1,header(ifile)%nvsfc
      read(tokens(nt  ),*,err=3,end=3) rsfc_obs%sdata(nv)
      read(tokens(nt+1),*,err=3,end=3) iqfl(nv)
      !print*,'(((((((read:',nt,nv,rsfc_obs%sdata(nv)
      nt = nt + 2
   enddo

   3 continue
   
   do nv=1,header(ifile)%nvsfc
      rsfc_obs%iqflags(nv,1)=iqfl(nv)/100
      rsfc_obs%iqflags(nv,2)=mod(iqfl(nv)/10,10)
      rsfc_obs%iqflags(nv,3)=mod(iqfl(nv),10)
   enddo

endif

! See if this station is on the "good" list. If we have the list 
!   and the station isn't on it, go get next obs.

if(ngood_gdf > 0) then
   do ngood=1,ngood_gdf
      if(trim(rsfc_obs%id) == trim(good_id_gdf(ngood)) ) then
         print*,'Sfc good obs found:',trim(good_id_gdf(ngood))
         go to 101
      endif
   enddo
   go to 100
endif

101 continue

! See if this station is on the "bad" list. If so, go get next obs.

if(nbad_gdf > 0) then
   do nbad=1,nbad_gdf
      if(trim(rsfc_obs%id) == trim(bad_id_gdf(nbad)) ) then
         print*,'Sfc bad obs found:',trim(bad_id_gdf(nbad))
         go to 100
      endif
   enddo
endif

! If desired, check QC flags and set appropriate values to missing.

if(qcheck == 'yes') then
   do nv=1,header(ifile)%nvsfc
      iqf=1
      do iq=1,3
         if((rsfc_obs%iqflags(nv,iq) /= 5 .and.  &
             rsfc_obs%iqflags(nv,iq) /= 0 .and.  &
             rsfc_obs%iqflags(nv,iq) /= 8) .or. &
             rsfc_obs%sdata(nv) < -998.) iqf=0
      enddo
      if(iqf == 0) rsfc_obs%sdata(nv)=-999.
   enddo

endif

10 continue

return

20 continue
ierr=1

return
END SUBROUTINE gdf_read_sfc_obs

!##############################################################################
Subroutine gdf_sfc_data_convert (varn,cvars,nvars)

use gdf_input

implicit none

integer :: nvars
real :: varn(nvars)
character(len=*) :: cvars(nvars)
character(len=16) :: cvar
real :: v1,v2,v3,vv
integer :: nv,ierr
real, external :: gdf_get_sfc_val
real, external :: rsatmix


! Convert units and type of sfc input data

do nv=1,nvars
   cvar=cvars(nv)
   varn(nv)=-999.
   ierr = 0
   if(trim(cvar)=='ue') then
      ! earth-relative u in m/s
      ! Check first for actual component
      varn(nv) = gdf_get_sfc_val('WIND_U_COMPONENT','m/s',ierr)
      if (ierr == 1) then
         ierr = 0
         ! Try to compute from dir and speed
         v1 = gdf_get_sfc_val('WIND_DIRECTION','m/s',ierr)
         v2 = gdf_get_sfc_val('WINDSPEED','m/s',ierr)
         if(v1 > -998. .and. v2 > -998.)  &
            CALL winduv (v1,v2,varn(nv),vv)
      endif
   
   elseif(trim(cvar)=='ve') then  
      ! earth-relative v in m/s
      ! Check first for actual component
      varn(nv) = gdf_get_sfc_val('WIND_V_COMPONENT','m/s',ierr)
      if (ierr == 1) then
         ierr = 0
         ! Try to compute from dir and speed
         v1 = gdf_get_sfc_val('WIND_DIRECTION','m/s',ierr)
         v2 = gdf_get_sfc_val('WINDSPEED','m/s',ierr)
         if(v1 > -998. .and. v2 > -998.)  &
            CALL winduv (v1,v2,vv,varn(nv))
      endif
   
   elseif(trim(cvar)=='speed') then  
      ! wind speed in m/s
      varn(nv) = gdf_get_sfc_val('WINDSPEED','m/s',ierr)
      if (ierr == 1) then
         ierr = 0
         v1 = gdf_get_sfc_val('WIND_U_COMPONENT','m/s',ierr)
         v2 = gdf_get_sfc_val('WIND_V_COMPONENT','m/s',ierr)
         if(v1 > -998. .and. v2 > -998.)  &
            CALL winddf (vv,varn(nv),v1,v2)
      endif
   
   elseif(trim(cvar)=='direction') then  
      ! wind direction
      varn(nv) = gdf_get_sfc_val('WIND_DIRECTION','deg',ierr)
      if (ierr == 1) then
         ierr = 0
         v1 = gdf_get_sfc_val('WIND_U_COMPONENT','m/s',ierr)
         v2 = gdf_get_sfc_val('WIND_V_COMPONENT','m/s',ierr)
         if(v1 > -998. .and. v2 > -998.)  &
            CALL winddf (varn(nv),vv,v1,v2)
      endif
   
   elseif(trim(cvar)=='tempc') then  
      ! temperature in C
      varn(nv) = gdf_get_sfc_val('TEMPERATURE','C',ierr)
   
   elseif(trim(cvar)=='tempf') then  
      ! temperature in F
      v1 = gdf_get_sfc_val('TEMPERATURE','C',ierr)
      if(v1 > -998.)  varn(nv)=v1*1.8+32.
   
   elseif(trim(cvar)=='dewptc') then  
      ! dewpoint in C
      varn(nv) = gdf_get_sfc_val('DEWPOINT','C',ierr)
   
   elseif(trim(cvar)=='dewptf') then  
      ! dewpoint in F
      v1 = gdf_get_sfc_val('DEWPOINT','C',ierr)
      if(v1 > -998.) varn(nv)=v1*1.8+32.
   
   elseif(trim(cvar)=='press') then  
      ! pressure in mb
      v1 = gdf_get_sfc_val('STN_PRES','Pa',ierr)
      if(v1 > -998.) varn(nv)=v1*.01
   
   elseif(trim(cvar)=='press_pa') then  
      ! pressure in mb
      varn(nv) = gdf_get_sfc_val('STN_PRES','Pa',ierr)
   
   elseif(trim(cvar)=='relhum') then  
      ! rh in percent
      v1 = gdf_get_sfc_val('TEMPERATURE','C',ierr)
      v2 = gdf_get_sfc_val('DEWPOINT','C',ierr)
      v3 = gdf_get_sfc_val('STN_PRES','Pa',ierr)
      if(v1 > -998. .and. v2 > -998. .and. v3 > -998.) &
      varn(nv)=100.*min(1.,max(0.,rsatmix(v3,v2+273.16)/rsatmix(v3,v1+273.16)))
               
   elseif(trim(cvar)=='ch-so2-ug') then  
      ! SO2
      varn(nv) = gdf_get_sfc_val('SO2','ug/m3',ierr)
   
   elseif(trim(cvar)=='ch-no2-ug') then  
      ! NO2
      varn(nv) = gdf_get_sfc_val('NO2','ug/m3',ierr)
   
   elseif(trim(cvar)=='ch-no-ug') then  
      ! NO
      varn(nv) = gdf_get_sfc_val('NO','ug/m3',ierr)
   
   elseif(trim(cvar)=='ch-pm10-ug') then  
      ! PM10
      varn(nv) = gdf_get_sfc_val('PM10','ug/m3',ierr)
   
   elseif(trim(cvar)=='ch-o3-ug') then  
      ! O3
      varn(nv) = gdf_get_sfc_val('O3','ug/m3',ierr)
      !print*,'======here:',ierr,nv,varn(nv)
   elseif(trim(cvar)=='ch-voc-ug') then  
      ! O3
      varn(nv) = gdf_get_sfc_val('VOC','ug/m3',ierr)
      !print*,'======here:',ierr,nv,varn(nv)
   elseif(trim(cvar)=='em-no') then  
      varn(nv) = gdf_get_sfc_val('EM-NO','moles/hr',ierr)
   elseif(trim(cvar)=='em-no2') then  
      varn(nv) = gdf_get_sfc_val('EM-NO2','moles/hr',ierr)
   elseif(trim(cvar)=='em-so2') then  
      varn(nv) = gdf_get_sfc_val('EM-SO2','moles/hr',ierr)
   
   
   else
      print*,'UNKNOWN CONVERT VARIABLE in sfc_data_convert !!!!:',trim(cvar)
      stop 'sfc_data_convert'
   endif

   if (ierr > 0) then
      print*,'gdf_sfc_data_convert: unable to find necessary info on GDF file'
      print*,'gdf_sfc_data_convert: for variable: ',trim(cvar)
   endif

enddo

return
END SUBROUTINE gdf_sfc_data_convert

!##############################################################################
real Function gdf_get_sfc_val (varn,units,ierr)

use gdf_input

implicit none

character(len=*) :: varn,units
integer :: ierr

integer :: nv

gdf_get_sfc_val = -999.

! Find which slot var is in

do nv = 1, header(curfile)%nvsfc
!print*,'||',trim(varn),'||',trim(header(curfile)%sfc_string(nv)),'||'
   if (trim(varn) == trim(header(curfile)%sfc_string(nv))) then
      gdf_get_sfc_val = rsfc_obs%sdata(nv)
      print*,'|| Found:',nv,' ',trim(varn),trim(units),gdf_get_sfc_val
      return
   endif
enddo

ierr = 1

return
END FUNCTION gdf_get_sfc_val
