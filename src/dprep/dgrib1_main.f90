!##############################################################################
Program griber

implicit none

character*10 c_date,c_hr,c_datatype
character*128 filein,fileout
character*2 cflag(10)
character*20 projection
character(len=20), allocatable :: fields(:),var3d(:), var2d(:), varsd(:)
integer, allocatable :: levs(:),irecs(:),idates(:),ifhrs(:)
real, allocatable :: a(:,:)

integer :: nx,ny,nrec,longdate,mode &
          ,ia,iargc,istr,iend,inc,i,j,lcg,ndate,nhr &
          ,nlev,slev,inproj,miss,startlev,endlev &
          ,iyyyy,imm,idd,ihh,fyyyy,fmm,fdd,fhh,leapyear &
          ,nvar3d,nvar2d,nvarsd,iplevs(1000),islevs(1000),datatype &
          ,writesoillevs,writesnowlevs,mx,plev1(1000),plev2(1000)
real :: alat1,alon1,alat2,alon2,alov,aorient,dx,dy,reflat1,reflat2 &
       ,tinc,amiss,alatscr
logical :: llisglobal,matched
data amiss /-999./
llisglobal = .true.
!******************************************************************************
!Check for command line arguments
!******************************************************************************
IA = IARGC ( )
IF (IA .lt. 8) THEN
   WRITE(*,'("Usage:dgrib [-t #] [-d YYYYMMDDHH] [-h fcsthr] [-f gribfile]")')
   WRITE(*,'("-t option:")')
   WRITE(*,'("  # = 1: NcepReanalysis 1 or 2 Global Lat/lon at 2.5-deg")')
   WRITE(*,'("         NO soil/snow; Using R / RELH")')
   WRITE(*,'("         Upper air names: UGRD,VGRD,TMP,HGT,RH")')
   WRITE(*,'("         SWE,SnowDepth names: NOTHING,NOTHING")')
   WRITE(*,'("         Soil Moisture,Temp names: NOTHING,NOTHING")')
   WRITE(*,'("  # = 2: GDAS-FNL Global Lat/Lon at 1.0deg UNTIL Jan 14, 2015 00Z")')
   WRITE(*,'("         YES soil/snow; Using R / RELH")')
   WRITE(*,'("         Upper air names: UGRD,VGRD,TMP,HGT,RH")')
   WRITE(*,'("         SWE,SnowDepth names: WEASD,WEASD")')
   WRITE(*,'("         Soil Moisture,Temp names: SOILW,TMP")')
   WRITE(*,'("  # = 3: GDAS-FNL Global Lat/Lon at 1.0deg AFTER Jan 14, 2015 00Z")')
   WRITE(*,'("  # = 3: GDAS-FNL Global Lat/Lon at 0.25deg after Jul 8, 2015")')
   WRITE(*,'("  # = 3: Forecast grids for GFS and HRRR Lat/Lon Grid")')
   WRITE(*,'("  # = 3: RAP Rapid Refresh Analysis 32km Awips Grid 221")')
   WRITE(*,'("         YES soil/snow; Using R / RELH")')
   WRITE(*,'("         Upper air names: UGRD,VGRD,TMP,HGT,RH")')
   WRITE(*,'("         SWE,SnowDepth names: WEASD,SNOD")')
   WRITE(*,'("         Soil Moisture,Temp names: SOILW,TSOIL")')
   WRITE(*,'("  # = 4: NARR 32km Lambert Conformal Grid")')
   WRITE(*,'("         YES soil/snow; Using Q / SPFH")')
   WRITE(*,'("         Upper air names: UGRD,VGRD,TMP,HGT,SPFH")')
   WRITE(*,'("         SWE,SnowDepth names: WEASD,SNOD")')
   WRITE(*,'("         Soil Moisture,Temp names: SOILW,TSOIL")')
   WRITE(*,'("  # = 5: ERA-Interim Global Lat/Lon at various resolutions")')
   WRITE(*,'("  # = 5: ERA5 Global Lat/Lon or area subset at 0.25deg")')
   WRITE(*,'("         YES soil/snow; Using R / RELH")')
   WRITE(*,'("         Upper air names: U,V,T,Z,R")')
   WRITE(*,'("         SWE,SnowDepth names: SD,SD")')
   WRITE(*,'("         Soil Moisture,Temp names: SWVL2,SWVL1,STL2,STL1")')
   WRITE(*,'("  # = 6: ERA5 Lat/Lon Global or area subset at 0.25deg")')
   WRITE(*,'("         YES soil/snow; Using Q / SPFH")')
   WRITE(*,'("         Upper air names: U,V,T,Z,Q")')
   WRITE(*,'("         SWE,SnowDepth names: SD,SD")')
   WRITE(*,'("         Soil Moisture,Temp names: SWVL2,SWVL1,STL2,STL1")')
   WRITE(*,'("-d option: YYYYMMDDHH is the time")')
   WRITE(*,'("-h option: forecast hour (should be 0 for reanalysis)")')
   WRITE(*,'("-f option: gribfile is the name of the file to extract data")')
   WRITE(*,'("See Dprep Documentation with RAMS for info on data in Grib files.")')
   STOP
ENDIF

!Set default date and forecast hour
c_date='99999999'
c_hr  ='99999999'

!Read command line arguments
do i=1,ia,2
   CALL ugetarg (i,cflag(i))
   if(cflag(i).eq.'-t')then
      CALL ugetarg (i+1,c_datatype)
   elseif(cflag(i).eq.'-d')then
      CALL ugetarg (i+1,c_date)
   elseif(cflag(i).eq.'-h')then
      CALL ugetarg (i+1,c_hr)
   elseif(cflag(i).eq.'-f')then
      CALL ugetarg (i+1,filein)
   endif
enddo

!wgrib only prints out the 2 digit year in short inventory
read(c_date(3:),'(i10)')ndate
read(c_hr,'(i10)')nhr
read(c_datatype,'(i10)')datatype
!Print out command line inputs
lcg=index(filein,' ')-1
print*,'gribfile=',filein(1:lcg)
print*,'c_date=',c_date,ndate
print*,'c_hr=',c_hr,nhr
print*,'c_datatype=',c_datatype,datatype

!******************************************************************************
! Datsets numbering described above in "WRITE" statements.
!
! Note the following:
!
! 1. ERA5 relative humidity is with respect to water or ice depending on 
!    temperature. But RAMS does initialization assuming RH with respect to 
!    water, so estimates can be high at colder temperatures. Better to use 
!    ERA5 specific humidity. Yes soil/snow.
!
! 2. Note that soil/snow data MUST be on the same grid structure (e.g. lat-lon,
!    lambert-conformal) and resolution (e.g. 1.0-deg, 0.25-deg) as the atmospheric
!    pressure level data.
!
! 3. ERA5 data found at the following on the Copernicus data server:
!
! Pressure level reanalysis:
! "ERA5 hourly data on pressure levels from 1979 to present":
! https://cds.climate.copernicus.eu/cdsapp#!/dataset
!  /reanalysis-era5-pressure-levels?tab=form
!
! Snow, soil moisture, soil temperatuere reanalysis:
! "ERA5 hourly data on single levels from 1979 to present":
! https://cds.climate.copernicus.eu/cdsapp#!/dataset
!  /reanalysis-era5-single-levels?tab=form
!******************************************************************************

!******************************************************************************
!NCEPreanalysis global 2.5deg lat-lon
if(datatype==1)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='UGRD'
 var3d(2)='VGRD'
 var3d(3)='TMP'
 var3d(4)='HGT'
 var3d(5)='RH'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='NOTHING'        !Accum. snow [kg/m^2] / SWE(mm)
 var2d(2)='NOTHING'        !Snow depth [m]
 !SOIL moisture variables here if available
 nvarsd=2
 allocate (varsd(nvarsd))
 varsd(1)='NOTHING'        !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='NOTHING'        !Soil temperature (K)

!******************************************************************************
!GDAS-FNL global 1.0-deg grids
! (Valid for GDAS before and including 2015 Jan 14 00Z)
! (1 soil level thru 2005053106. Added 3 more 2005053112.)
! Layers 0-10cm, 10-40cm, 40-100cm, 100-200cm soil (after 2005 Jun 01 00Z)
! Layer 0-10cm only (before 2005 Jun 01 00Z)
! depths so write out top 2 levels expected by RAMS
! No SNOD(snow depth) variable up to this time. Will estimate SNOD(m depth).
elseif(datatype==2)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='UGRD'
 var3d(2)='VGRD'
 var3d(3)='TMP'
 var3d(4)='HGT'
 var3d(5)='RH'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='WEASD'          !Accum. snow [kg/m^2] / SWE(mm)
 var2d(2)='WEASD'          !Snow depth [m]
 !SOIL moisture variables here if available
 nvarsd=2
 allocate (varsd(nvarsd))
 varsd(1)='SOILW'          !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='TMP'            !Soil temperature (K)

!******************************************************************************
!GDAS-FNL global 1.0-deg,0.25-deg; HRRR,GFS forecast grids (lat-lon)
! (Valid for GDAS after and including 2015 Jan 14 06Z)
! These data have 0-10cm, 10-40cm, 40-100cm, 100-200cm soil
! depths so write out top 2 levels expected by RAMS
!RAP Rapid Refresh Analysis 32km (awips grid 221)
! This data has 0cm,1cm,4cm,10cm,30cm,60cm,100cm,160cm,300cm soil
! depths so write out top 2 levels expected by RAMS
elseif(datatype==3)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='UGRD'
 var3d(2)='VGRD'
 var3d(3)='TMP'
 var3d(4)='HGT'
 var3d(5)='RH'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='WEASD'          !Accum. snow [kg/m^2] / SWE(mm)
 var2d(2)='SNOD'           !Snow depth [m]
 !SOIL moisture variables here if available
 nvarsd=2
 allocate (varsd(nvarsd))
 varsd(1)='SOILW'          !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='TSOIL'          !Soil temperature (K)

!******************************************************************************
!NARR reanalysis 32km (lambert-conformal grid)
! Make sure RAMSIN namelist knows you are using Q/SPFH (IDATAIN=1).
! The NARR 32km data has 0-10cm, 10-40cm, 40-100cm, 100-200cm soil
! depths so write out top 2 levels expected by RAMS
elseif(datatype==4)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='UGRD'
 var3d(2)='VGRD'
 var3d(3)='TMP'
 var3d(4)='HGT'
 var3d(5)='SPFH'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='WEASD'          !Accum. snow [kg/m^2] / SWE(mm)
 var2d(2)='SNOD'           !Snow depth [m]
 !SOIL moisture variables here if available
 nvarsd=2
 allocate (varsd(nvarsd))
 varsd(1)='SOILW'          !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='TSOIL'          !Soil temperature (K)

!******************************************************************************
!ECMWF-ERA-Interim global 1.5deg (lat-lon)
!ERA5 0.25-deg (global lat-lon or subset with pressure levels and soil/snow on
! same grid resolution. Using R/RELH rather than Q/SPFH for moisture.
! Top 2 ERA-Interim soil levels are 7-28cm and 0-7cm down
! Note that SD/SNOD here is meters of water equivalent (be careful)
elseif(datatype==5)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='U'
 var3d(2)='V'
 var3d(3)='T'
 var3d(4)='Z'    !Geopotential, not Geopotential Height(m)
 var3d(5)='R'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='SD'           !WEASD estimated from SNOD
 var2d(2)='SD'           !Snow depth [m water equivalent](be careful)
 !SOIL moisture variables here if available
 nvarsd=4
 allocate (varsd(nvarsd))
 varsd(1)='SWVL2'         !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='SWVL1'
 varsd(3)='STL2'          !Soil temperature (K)
 varsd(4)='STL1'

!******************************************************************************
!ECMWF-ERA5 (lat-lon) (Global or subset)
! Using specific humidity (Q or SPFH) instead of relative humidity (R or RELH)
! Make sure RAMSIN namelist knows you are using Q/SPFH (IDATAIN=1).
! Top 2 ERA-Interim soil levels are 7-28cm and 0-7cm down
! Note that SD/SNOD here is meters of water equivalent (be careful)
elseif(datatype==6)then
 !3D atmospheric variables
 nvar3d=5
 allocate (var3d(nvar3d))
 var3d(1)='U'
 var3d(2)='V'
 var3d(3)='T'
 var3d(4)='Z'    !Geopotential, not Geopotential Height(m)
 var3d(5)='Q'
 !2D snow water and depth variables here if available
 nvar2d=2
 allocate (var2d(nvar2d))
 var2d(1)='SD'           !WEASD estimated from SNOD
 var2d(2)='SD'           !Snow depth [m water equivalent](be careful)
 !SOIL moisture variables here if available
 nvarsd=4
 allocate (varsd(nvarsd))
 varsd(1)='SWVL2'         !Volumetric Soil Moisture (m3/m3)(fraction)
 varsd(2)='SWVL1'
 varsd(3)='STL2'          !Soil temperature (K)
 varsd(4)='STL1'

endif
!******************************************************************************

!Read some of the grib metadata
CALL grib_query  (filein(1:lcg))            !Open the grib file
CALL grib_queryc ('projection',projection)  !get Grid projection
CALL grib_queryf ('lat1',alat1)             !get SW lat
CALL grib_queryf ('lon1',alon1)             !get SW lon
CALL grib_queryf ('lat2',alat2)             !get NE lat
CALL grib_queryf ('lon2',alon2)             !get NE lon
CALL grib_queryf ('dx',dx)                  !get Delta-X (km or deg)
CALL grib_queryf ('dy',dy)                  !get Delta-Y (km or deg)
CALL grib_queryi ('nx',nx)                  !get # of X grid points
CALL grib_queryi ('ny',ny)                  !get # of Y grid points
allocate (a(nx,ny))
CALL grib_queryi ('nrec',nrec)              !get # of records
CALL grib_queryi ('longdate',longdate)      !get the date/time
CALL grib_queryi ('mode',mode)              !get grid mode info (see grid tables)
CALL grib_queryf ('orient',aorient)         
CALL grib_queryf ('lov',alov)               !get grid orientation for LC projection
CALL grib_queryf ('latin1',reflat1)         !get reference lat 1 for LC projection
CALL grib_queryf ('latin2',reflat2)         !get reference lat 2 for LC projection
allocate (irecs(nrec))
allocate (levs(nrec))
allocate (idates(nrec))
allocate (ifhrs(nrec))
allocate (fields(nrec))
CALL grib_queryi ('irecs',irecs)
CALL grib_queryi ('levs',levs)
CALL grib_queryi ('idates',idates)
CALL grib_queryi ('ifhrs',ifhrs)
CALL grib_queryc ('fields',fields)

print*,'DX,DY: ',dx,dy
print*,'NX,NY: ',nx,ny
print*,'LAT1,LON1: ',alat1,alon1
print*,'LAT2,LON2: ',alat2,alon2
print*,'PROJECTION: ',projection
print*,'LOV: ',alov                   !??? grid rotation
print*,'Orient: ',aorient             !??? grid rotation
print*,'Mode: ',mode
print*,'Nrec: ',nrec
print*,'Irecs: ',irecs(nrec)







!Do the dprep stuff

!First some grid navigation bookeeping
if (projection(1:6)=="latlon") then
   inproj=1
   !Will have to flip the latlon projections,
   !Will modify start lat/lon to be -90;-180
   !RAMS wants coordinates -180 to 180 longitude (not 0 to 360)
   if(mode==128)then
    if(alat1==90.0)then !For truly global lat-lon data
     llisglobal = .true.
     alat1=-alat1
     alon1=alon1-180.
     alat2=-alat2
     if(alon2 <  180.) alon2=alon2+180.
     if(alon2 >= 180.) alon2=alon2-180.
    else
     !If we have a subset of the data, we still need to switch lat1 and lat2. 
     !We should also save that this is not global.
     alatscr = alat1
     alat1 = alat2
     alat2 = alatscr
     llisglobal = .false.
    endif
   endif
elseif (projection(1:7)=="Lambert") then
   inproj=2
   alat2=0.0 !NE lat
   alon2=0.0 !NE lon
   dx=dx*1000.
   dy=dy*1000.
   aorient=alov
elseif (projection(1:5)=="polar") then
   inproj=3
   !Assume the polar stereo is true polar stereo
   reflat1=90.
   reflat2=amiss
endif
print*,'PERHAPS MODIFIED LAT/LON REPORTED HERE'
print*,'LAT1,LON1: ',alat1,alon1
print*,'LAT2,LON2: ',alat2,alon2

!Quick check to see if there are any fields matching time wanted
j=0
do i=1,nrec
   if(('HGT'==fields(i).or.'GP'==fields(i).or.'Z'==fields(i)) &
       .and.ndate==idates(i) &
       .and.nhr==ifhrs(i))then
      j=1
      exit
   endif
enddo
if(j==0)stop "NOTHING FOUND"

!**************************************************************************
!Determine the number of Pressure levels we are going to use
nlev=0
do i=1,nrec
   if(levs(i)==0)cycle  !don't count the surface fields here
   if((fields(i)=='HGT'.or.fields(i)=='GP'.or.fields(i)=='Z') &
       .and.ndate==idates(i) &
       .and.nhr==ifhrs(i))then
      nlev=nlev+1
      iplevs(nlev)=levs(i)
   endif
enddo

if(iplevs(1) < iplevs(nlev))then
   nlev=0
   do i=nrec,1,-1
      if(levs(i)==0)cycle  !don't count the surface fields here
      if((fields(i)=='HGT'.or.fields(i)=='GP'.or.fields(i)=='Z') &
          .and.ndate==idates(i) &
          .and.nhr==ifhrs(i))then
         nlev=nlev+1
         iplevs(nlev)=levs(i)
      endif
   enddo
endif
print*,'Number of Pressure Levels= ',nlev

!**************************************************************************
!Make sure the pressure levels are in descending order for correct output
!if not ordered properly in grib file - plev1, plev23, and mx added to
!integer variables declarations (Stephen Noble, SRNL)
plev1=iplevs
do i=1,nlev
   mx=maxval(plev1)
   plev2(i)=mx
   do j=1,nlev
      if (plev1(j)==mx) then
         plev1(j)=-999
      endif
   enddo
enddo
iplevs=plev2

!**************************************************************************
!Determine the number of SOIL levels we are going to use.
!Only add soil (moisture and/or temperature) if the soil moisture is there.
slev=0
do i=1,nrec
   if((fields(i)=='SOILW'.or.fields(i)=='VSOILM').and.ndate==idates(i) &
       .and.nhr==ifhrs(i))then
      slev=slev+1
      islevs(slev)=levs(i)
   endif
enddo
!Invert soil levels so level=1 is the bottom soil level (similar to RAMS)
if(abs(islevs(1)) < abs(islevs(slev)))then
   slev=0
   do i=nrec,1,-1
      if((fields(i)=='SOILW'.or.fields(i)=='VSOILM').and.ndate==idates(i) &
          .and.nhr==ifhrs(i))then
         slev=slev+1
         islevs(slev)=levs(i)
      endif
   enddo
endif
!Check for independently named soil levels such as for ERA-Interim
!Order to that level=1 is the bottom soil level (similar to RAMS)
do i=1,nrec
   if((fields(i)=='SWVL1').and.ndate==idates(i).and.nhr==ifhrs(i))then
      slev=slev+1
      islevs(slev)=levs(i)
      print*,'ERA',slev,islevs(slev)
   endif
   if((fields(i)=='SWVL2').and.ndate==idates(i).and.nhr==ifhrs(i))then
      slev=slev+1
      islevs(slev)=levs(i)
      print*,'ERA',slev,islevs(slev)
   endif
enddo

print*,'Number of Soil Levels= ',slev

!**************************************************************************
!Set times for writing to file and dp header info
read(c_date(1:4) ,'(i4)')iyyyy
read(c_date(5:6) ,'(i2)')imm
read(c_date(7:8) ,'(i2)')idd
read(c_date(9:10),'(i2)')ihh

leapyear=0
!leap year check. Every 4th year is a leap year except years divisible by 100
!unless they are also divisible by 400.
if(  mod(iyyyy, 4)==0) then
  !we could be a leap year. Keep checking.
  leapyear=1
  if(mod(iyyyy, 100)==0) then
    !maybe not
    leapyear=0
    if(mod(iyyyy,400)==0) then
      !yes leap year
      leapyear=1
    endif
  endif
endif
print*,'Leapyear= ',leapyear
!adjust date if using forecast grids
fyyyy = iyyyy
fmm = imm
fhh = ihh + mod(nhr,24)
fdd = idd + (nhr/24)
if(fhh>=24) then
  fdd = fdd+fhh/24
  fhh = mod(fhh,24)
endif
if(imm==1 .and. fdd>31)then
 fmm=2
 fdd=fdd-31
elseif(imm==2 .and. fdd>28 .and. leapyear==0)then
 fmm=3
 fdd=fdd-28
elseif(imm==2 .and. fdd>29 .and. leapyear==1)then
 fmm=3
 fdd=fdd-29
elseif(imm==3 .and. fdd>31)then
 fmm=4
 fdd=fdd-31
elseif(imm==4 .and. fdd>30)then
 fmm=5
 fdd=fdd-30
elseif(imm==5 .and. fdd>31)then
 fmm=6
 fdd=fdd-31
elseif(imm==6 .and. fdd>30)then
 fmm=7
 fdd=fdd-30
elseif(imm==7 .and. fdd>31)then
 fmm=8
 fdd=fdd-31
elseif(imm==8 .and. fdd>31)then
 fmm=9
 fdd=fdd-31
elseif(imm==9 .and. fdd>30)then
 fmm=10
 fdd=fdd-30
elseif(imm==10 .and. fdd>31)then
 fmm=11
 fdd=fdd-31
elseif(imm==11 .and. fdd>30)then
 fmm=12
 fdd=fdd-30
elseif(imm==12 .and. fdd>31)then
 fmm=1
 fdd=fdd-31
 fyyyy=iyyyy+1
endif

write(fileout,45 ) 'dp-p',fyyyy,'-',fmm,'-',fdd,'-',fhh,'00'

print*,'Writing file:',fileout
45   format(a4,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a2)

open(10,file=fileout,status='unknown',form='formatted')

!Do the following 2 lines to make values between -180 to +180 for RAMS
if(alon1>180.) alon1=alon1-360.
if(aorient>180.) aorient=aorient-360.

WRITE(10,900)iyyyy,imm,idd,ihh*100,nhr*100,nlev,nx,ny &
     ,inproj,dx,dx,alat1,alon1,alat2,alon2,reflat1,aorient &
     ,(iplevs(i),i=1,nlev)
900  FORMAT(8I6,/,i1,2x,8(F10.4,2x),/,(12I6))

!*********************************************************************
!*********************************************************************
!     Go through each variable at each level for writing.

!Uncomment if you are pulling 3D variables from grib (Saleeby 2-28-19)
!Must output all 5 3D vars even if null data from missing RH.
do j=1,nlev
 do i=1,nvar3d
   a=amiss
   matched = .false.
   CALL getfield (a,var3d(i),nrec,irecs,fields,levs,idates,ifhrs &
                ,iplevs(j),ndate,nhr,matched)
   CALL prepfield (a,nx,ny,projection,llisglobal,var3d(i),mode &
                 ,amiss,'3d',datatype,i)
!write(10,*)'3d NEXT ',var3d(i)
   write(10,'(8f12.5)')a
 enddo
enddo

!Uncomment if you are pulling SOIL variables from grib (Saleeby 2-28-19)
!Datasets vary in the number of soil levels they have. Write out the top
!2 levels if we have 2+. If we only have 1 soil level, write it 2 times.
!RAMS wants to ingest 2 soil levels. If data fields are empty RAMS uses
!a default.
writesoillevs=0
do i=1,nvarsd
  startlev=max(1,slev-1)
  endlev=max(1,slev)
  do j=startlev,endlev
     a=amiss
     matched = .false.
     CALL getfield (a,varsd(i),nrec,irecs,fields,levs,idates,ifhrs &
                  ,islevs(j),ndate,nhr,matched)
     if(matched) then
      CALL prepfield (a,nx,ny,projection,llisglobal,varsd(i),mode &
                    ,amiss,'sd',datatype,i)
!write(10,*)'SOIL NEXT ',varsd(i),islevs(j)
      write(10,'(8f12.5)')a
      writesoillevs=writesoillevs+1
     endif
  enddo
  !Write out at least 2 levels. If only 1 soil level exists, repeat it.
  if(startlev==endlev)then
    do j=startlev,endlev
      a=amiss
      matched = .false.
      CALL getfield (a,varsd(i),nrec,irecs,fields,levs,idates,ifhrs &
                   ,islevs(j),ndate,nhr,matched)
      if(matched) then
       CALL prepfield (a,nx,ny,projection,llisglobal,varsd(i),mode &
                     ,amiss,'sd',datatype,i)
!write(10,*)'SOIL NEXT EXTRA ',varsd(i)
       write(10,'(8f12.5)')a
       writesoillevs=writesoillevs+1
      endif
    enddo
  endif
enddo
!If no soil data found, then input place filler
!If partial soil data found, then stop cause it is not written correctly.
if(writesoillevs < 4) then
 if(writesoillevs == 0)then
  do i=1,4 !2 soil moisture layers and then 2 soil temperature layers
    a=amiss
    print*,'SOIL Place Holder ',i
!write(10,*)'SOIL NEXT ',i
    write(10,'(8f12.5)')a
  enddo
 else
   print*,'INCONSISTENCY IN SOIL LAYERS - CHECK DATA'
   stop
 endif
endif

!Uncomment if you are pulling 2D SNOW variables from grib (Saleeby 2-28-19)
!These are SWE(kg/m2,mm) and snow depth(m)
writesnowlevs=0
do i=1,nvar2d
   a=amiss
   j=0 !specify j=0 for the surface
   CALL getfield (a,var2d(i),nrec,irecs,fields,levs,idates,ifhrs &
                ,j,ndate,nhr,matched)
   if(matched) then
    CALL prepfield (a,nx,ny,projection,llisglobal,var2d(i),mode &
                 ,amiss,'2d',datatype,i)
!write(10,*)'SNOW NEXT ',var2d(i)
    write(10,'(8f12.5)')a
    writesnowlevs=writesnowlevs+1
   endif
enddo
!If no snow data found, then input place filler
!If partial snow data found, then stop cause it is not written correctly.
if(writesnowlevs < 2) then
 if(writesnowlevs == 0)then
  do i=1,2 !snow water equivalent(kg/m2,mm) and snow depth(m)
    a=amiss
    print*,'SNOW Place Holder ',i
!write(10,*)'SNOW NEXT ',i
    write(10,'(8f12.5)')a
  enddo
 else
   print*,'INCONSISTENCY IN SNOW DATA - CHECK DATA'
   stop
 endif
endif

!*********************************************************************
!*********************************************************************

close(10)

!Write the tag file to indicate file is complete
open(10,file=fileout(1:len_trim(fileout))//'.tag',status='unknown' &
    ,form='formatted')
write(10,'(a)')'done'
close(10)

deallocate (a)
deallocate (fields)
deallocate (ifhrs)
deallocate (idates)
deallocate (levs)
deallocate (irecs)
deallocate (var3d)
if(allocated(var2d)) deallocate (var2d)
if(allocated(varsd)) deallocate (varsd)

END

!************************************************************************
Subroutine getfield (a,type,nrec,irecs,fields,levs,idates,ifhrs &
                     ,iplev,ndate,nhr,matched)

implicit none

real, dimension(*) :: a
character*(*) type, fields(*)
integer, dimension(*) :: irecs,levs,idates,ifhrs
integer :: iplev,ndate,nhr,nrec,i
logical :: matched

matched = .false.

do i=1,nrec,1

   if(type==fields(i).and.iplev==levs(i).and.ndate==idates(i).and. &
      nhr==ifhrs(i))then
      CALL grib_get (a,irecs(i))
      print*,'matched ',fields(i),levs(i),idates(i),ifhrs(i)
      matched = .true.
      return
   endif

enddo

return
END SUBROUTINE getfield

!************************************************************************
Subroutine prepfield (a,nx,ny,projection,llisglobal,type,mode,amiss,dimtype &
                     ,datatype,vars)

implicit none

integer nx,ny,datatype,vars
real, dimension(nx*ny) :: a
character*(*) type,projection,dimtype

real, allocatable :: b(:)
real amiss
integer mode
logical :: llisglobal
integer i

if(projection=='latlon'.and.mode==128)then
   allocate (b(nx*ny))
   if(llisglobal) then
    print*,'   flipping the projection'
    CALL flip (nx,ny,a,b)
   else
    print*,'   flipping the projection (latsonly)'
    CALL flip_lats (nx,ny,a,b)
   endif
   a=b
   deallocate (b)
endif

do i=1,nx*ny
   if(a(i) >= 1.e10)a(i)=amiss
   if(a(i) <= -9999.)then
    print*,'Some data values less than -9999. Stopping.'
    stop
   endif
enddo

!Turn from relative humidity in percent to a fraction
if(type=='RH' .or. type=='R')then
   print*,'   factoring RH'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/100.
   enddo

!ECMWF uses GEOPOTENTIAL as the geopotential height multiplied by the WMO
!defined gravity constant of 9.80665 m/s**2 which is constant for all 
!latitudes and all heights. We must adjust for geopotential height.
elseif(type=='GP' .or. type=='Z')then
   print*,'   convert geopotential to geopotential-height'
   do i=1,nx*ny
      a(i)=a(i)/9.80665
      if(a(i) < amiss) then
       print*,'Bad negative geopotential height. Check data',i,a(i)
       stop
      endif
   enddo

!Turn from Pa to hPa or millibars
elseif(type=='PRMSL'.or.type=='PRES')then
   print*,'   factoring to mb'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/100.
   enddo

!For GDAS-FNL snowdepth, Convert the 2nd WEASD to a snow-depth(m) 
!assuming density of 10:1 (Convert mm-liquid-equiv to m-snowdepth)
elseif(datatype==2.and.dimtype=='2d'.and.vars==2.and.type=='WEASD')then
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)/100.
   enddo

!For ERA-Int/ERA5, Convert the SD/SNOD snow-depth (meters of water equivalent)
!(Convert m-liquid-equiv to mm-liquid-equiv (kg/m2))(vars=1=WEASD)
elseif((datatype==5.or.datatype==6).and.dimtype=='2d'.and.vars==1.and.type=='SD')then
   print*,'   convert SD/SNOD(SWE in meters) to kg/m2(SWE in mm)'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)*1000.
   enddo

!For ERA-Int/ERA5, Convert the SD/SNOD snow depth (meters of water equivalent)
!assuming density of 10:1 (Convert m-liquid-equiv to m-snowdepth)(vars=2=SNOD)
elseif((datatype==5.or.datatype==6).and.dimtype=='2d'.and.vars==2.and.type=='SD')then
   print*,'   convert SD/SNOD(SWE in meters) to snow depth(m) assuming 10:1'
   do i=1,nx*ny
      if(a(i) > amiss)a(i)=a(i)*10.
   enddo
endif

return
END SUBROUTINE prepfield

!************************************************************************
Subroutine flip (nx,ny,vin,vout)

implicit none

integer nx,ny,i,j,ii,jj
real, dimension(nx,ny) :: vin, vout

!Adjust starting longitude. Grib starts at 0.0deg lon
!while RAMS starts at -180.0deg
do i=1,nx
  ii=nx/2+i
  if(ii.gt.nx)ii=i-nx/2
  !Flip N/S for grib1 format
  !wgrib says: WE:NS winds(N/S) 
  do j=1,ny
    jj=ny-j+1
    vout(i,j)=vin(ii,jj)
  enddo
enddo

return
END SUBROUTINE flip

!************************************************************************
Subroutine flip_lats (nx,ny,vin,vout)
!Flip GRIB1 data (which is in N/S format) to RAMS
!which expects this data to go from S to N. This only flips N/S
!and is used for GLOBAL DATA SUBSETS.
!Testing has shown we need to do this for lat-lon global data subset of
!ERA5 data. Limited area subset of the global data can be retrieved 
!by the API rather than by the Copernicus web server.

implicit none

integer nx,ny,i,j,jj
real, dimension(nx,ny) :: vin, vout

!Flip the latitudes only, do not flip the longitudes. 
do i=1,nx
  do j=1,ny
    jj=ny-j+1
    vout(i,j)=vin(i,jj)
  enddo
enddo

return
END SUBROUTINE flip_lats

