!##############################################################################
Subroutine rams_text (a,iztrans,ivtype,nngd,n1,n2,n3,fcstsec &
                   ,swlat,swlon,nelat,nelon,polelatn,polelonn,cvar)

! Saleeby: GEMPAK text format
!-----------------------------------------------------------------
! Routine to write RAMS fields in text format as chosen through
! REVUIN. This is intended to be modified at will to produce the
! desired output. Basically, one 3 or 2-dimensional variable at a
! time will be sent here.

! Arguments:
! ----------
! cdname - name of variable
! cdunits - units of variable
! a - data
! n1,n2,n3 - actual dimensions of field
! iztrans - type of vertical transformation 1-sigma-z, 2-Cartesian, 3-pressure
! ivtype - type of variable 2-2d surface, 3-3d atmospheric
! nib,nie,njb,nje - horizontal or vertical "window" as chosen in namelist.
! nplevs - number of atmospheric coordinate levels (if pressure transformation)
! iplevs(nplevs)- atmospheric coordinate levels (if pressure transformation)
! zlev - atmospheric coordinate levels (if sigma_z or Cartesian)
! iyear1 - year of model start
! imonth1 - month of model start
! idate1 - date of model start
! itime1 - time of model start (hhmm)
! fcstsec - seconds into run
!-----------------------------------------------------------------

use an_header
use rcommons
use mem_grid

implicit none

integer :: iztrans,ivtype,n1,n2,n3,nngd,ihour2,imin2,idate2,imonth2,iyear2 &
          ,jcount,i,j,k,lastslash,iyr2,typlev,jj,kk,zl,z
real :: fcstsec,swlat,swlon,nelat,nelon,polelatn,polelonn,fcstmin,fcsthrs &
          ,fcstday,fcstsec1,zlev
real :: a(n1,n2,n3)
character(len=4) :: coordn
character(len=6) :: ngem
character(len=40) :: cvar
character(len=strl1) :: flnm,flnm1,out1,out2
integer, save :: ncall(maxgrds),iun,iuntag,numgrids &
                ,xbeg,xend,ybeg,yend,xinc,yinc,xpts,ypts,iyr1
real, save :: swlat1,swlon1,nelat1,nelon1

data ncall/maxgrds*0/

!print*,'===> text out =>',fcstsec
!print*,iyear1,imonth1,idate1,itime1,nngd
!print*,cdname(1:len_trim(cdname)-1),n1,n2,n3
!print*,iztrans,ivtype
!print*,nib,nie,njb,nje,nnb,nne

! If it is the first time into this routine, make and open a file.
if(ncall(nngd).eq.0) then

  if(niinc > 2 .or. njinc > 2) then
    print*,'Can only use increment of 1 or 2 right now'   
    stop
  endif

  CALL rams_get_cdata (0,1,flnm1)
  write(flnm,'(2a,2a1,i4.4,a1,i2.2,a1,i2.2,a1,i6.6,a2,i1)' )  &
      revpref(1:len_trim(revpref))  &
     ,flnm1(lastslash(flnm1)+1:len_trim(flnm1)-27),ftran,'-'  &
     ,iyear1,'-',imonth1,'-',idate1,'-',itime1*100,'-g',nngd

  write(out1,'(a,a4)') trim(flnm),'.txt'
  iun=79
  open(iun,file=out1,status='unknown')  
  rewind iun

  write(out2,'(a,a4)') trim(flnm),'.tag'
  iuntag=78
  open(iuntag,file=out2,status='unknown')
  rewind iuntag

  swlat1 = swlat
  swlon1 = swlon
  nelat1 = nelat
  nelon1 = nelon
  xbeg = nib
  ybeg = njb
  xend = nie
  yend = nje
  xinc = niinc
  yinc = njinc

  xpts=0
  do z=xbeg,xend,xinc
   xpts=xpts+1
  enddo
  ypts=0
  do z=ybeg,yend,yinc
   ypts=ypts+1
  enddo

  if(nje==1) then
    ybeg = 1
    yend = 5
    yinc = 1
    ypts = 5
    nelat1 = swlat1 + 1.0
  endif

  !Saleeby(2020)
  !If using a tiny LES type grid, expand for GEMPAK.
  !Note that this could be problem for grid draped over the poles.
  if(nje>1 .and. abs(nelat1-swlat1)<0.1 .and. abs(nelon1-swlon1)<0.1)then
   nelat1 = swlat1 + 0.5
   nelon1 = swlon1 + 0.5 / cos(polelatn*3.14159/180.)
  endif

  if(iyear1 >= 2000) iyr1=iyear1-2000
  if(iyear1 <  2000) iyr1=iyear1-1900

  !write navigation file for creating gempak grid file
  write(iuntag,'(a7,a80)') &
   ,'GDEFIL=',out1
  write(iuntag,'(a7,3i2.2,i4.4,a1,i1,a9)') &
   ,'GDOUTF=',iyr1,imonth1,idate1,itime1,'g',nngd,'_rams.gem'
  write(iuntag,'(a7,3i2.2,i4.4,a1,i1,a9)') &
   ,'GDFILE=',iyr1,imonth1,idate1,itime1,'g',nngd,'_rams.gem'
  write(iuntag,'(a10,f8.2,a1,f8.2,a5)')'PROJ="str/',polelatn,';',polelonn,';0.0"'
  write(iuntag,'(a9,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,a1)') &
   ,'GRDAREA="',swlat1,';',swlon1,';',nelat1,';',nelon1,'"'
  write(iuntag,'(a7,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,a1)') &
   ,'GAREA="',swlat1,';',swlon1,';',nelat1,';',nelon1,'"'
  write(iuntag,'(a6,i5,a1,i5,a1)') 'KXKY="',xpts,';',ypts,'"'

  ncall(nngd)=1
  numgrids=0

endif

!************************************************************
!***** sets correct forecast time for gempak grid file *****
!************************************************************
ihour2=int(itime1/100.0)
imin2=itime1-ihour2*100.0
iyear2=iyear1
imonth2=imonth1
idate2=idate1
fcstsec1=fcstsec
fcstday=fcstsec/86400.0
if(fcstday >=1) then
 idate2=idate2+int(fcstday)
 fcstsec1=fcstsec1 - int(fcstday)*86400.0
endif

fcsthrs=fcstsec1/3600.0  
if(fcsthrs >=1) then
 ihour2=ihour2+int(fcsthrs)
 fcstsec1=fcstsec1 - int(fcsthrs)*3600.0
endif

fcstmin=fcstsec1/60.0
if(fcstmin >=1)then
 imin2=imin2+int(fcstmin)
 fcstsec1=fcstsec1 - int(fcstmin)*60.0
endif

if(imin2>=60) then
   imin2=imin2-60
   ihour2=ihour2+1
endif
if(ihour2>=24) then
   ihour2=ihour2-24
   idate2=idate2+1
endif  
if(imonth2 == 1.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 2.and.idate2 > 28.and.mod(iyear2,4)>0) then
   idate2=idate2-28
   imonth2=imonth2+1
endif
if(imonth2 == 2.and.idate2 > 29.and.mod(iyear2,4)==0) then
   idate2=idate2-29
   imonth2=imonth2+1
endif
if(imonth2 == 3.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 4.and.idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 5.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 6.and.idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 7.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 8.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 9.and.idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 10.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
endif
if(imonth2 == 11.and.idate2 > 30) then
   idate2=idate2-30
   imonth2=imonth2+1
endif
if(imonth2 == 12.and.idate2 > 31) then
   idate2=idate2-31
   imonth2=imonth2+1
   if(imonth2>12) then
     iyear2=iyear2+1
     imonth2=1
   endif
endif

if(fcstsec1 > 0.0) print*,'BAD TIME PROGRESSION: TEST CODE'

if(iyear2 >= 2000) iyr2=iyear2-2000
if(iyear2 <  2000) iyr2=iyear2-1900

!Specify the grid coordinate that the variable is native to.
zl=1 !Default use of ztn levels on T-grid

!EMPTY VARIABLES - 1 variables
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'empty3d')                 ngem='EMT3'

!AEROSOL AOD - 18 variables
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'ccn_dry_AOD_550')         ngem='ACND'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'ccn_wet_AOD_550')         ngem='ACNW'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust1_dry_AOD_550')       ngem='AD1D'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust1_wet_AOD_550')       ngem='AD1W'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust2_dry_AOD_550')       ngem='AD2D'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust2_wet_AOD_550')       ngem='AD2W'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero1_dry_AOD_550') ngem='AR1D'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero1_wet_AOD_550') ngem='AR1W'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero2_dry_AOD_550') ngem='AR2D'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero2_wet_AOD_550') ngem='AR2W'
if(len_trim(cvar).eq.21.and.cvar(1:21).eq.'salt_film_dry_AOD_550')   ngem='ASFD'
if(len_trim(cvar).eq.21.and.cvar(1:21).eq.'salt_film_wet_AOD_550')   ngem='ASFW'
if(len_trim(cvar).eq.20.and.cvar(1:20).eq.'salt_jet_dry_AOD_550')    ngem='ASJD'
if(len_trim(cvar).eq.20.and.cvar(1:20).eq.'salt_jet_wet_AOD_550')    ngem='ASJW'
if(len_trim(cvar).eq.22.and.cvar(1:22).eq.'salt_spume_dry_AOD_550')  ngem='ASSD'
if(len_trim(cvar).eq.22.and.cvar(1:22).eq.'salt_spume_wet_AOD_550')  ngem='ASSW'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'Total_dry_AOD_550')       ngem='ATOD'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'Total_wet_AOD_550')       ngem='ATOW'

!AEROSOL EXTINCTION COEFFICIENT - 18 variables
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'ccn_dry_ext_550')         ngem='ECND'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'ccn_wet_ext_550')         ngem='ECNW'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust1_dry_ext_550')       ngem='ED1D'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust1_wet_ext_550')       ngem='ED1W'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust2_dry_ext_550')       ngem='ED2D'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dust2_wet_ext_550')       ngem='ED2W'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero1_dry_ext_550') ngem='ER1D'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero1_wet_ext_550') ngem='ER1W'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero2_dry_ext_550') ngem='ER2D'
if(len_trim(cvar).eq.23.and.cvar(1:23).eq.'regen_aero2_wet_ext_550') ngem='ER2W'
if(len_trim(cvar).eq.21.and.cvar(1:21).eq.'salt_film_dry_ext_550')   ngem='ESFD'
if(len_trim(cvar).eq.21.and.cvar(1:21).eq.'salt_film_wet_ext_550')   ngem='ESFW'
if(len_trim(cvar).eq.20.and.cvar(1:20).eq.'salt_jet_dry_ext_550')    ngem='ESJD'
if(len_trim(cvar).eq.20.and.cvar(1:20).eq.'salt_jet_wet_ext_550')    ngem='ESJW'
if(len_trim(cvar).eq.22.and.cvar(1:22).eq.'salt_spume_dry_ext_550')  ngem='ESSD'
if(len_trim(cvar).eq.22.and.cvar(1:22).eq.'salt_spume_wet_ext_550')  ngem='ESSW'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'Total_dry_ext_550')       ngem='ETOD'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'Total_wet_ext_550')       ngem='ETOW'

!3D VELOCITY AND VORTICITY VARIABLES - 21 variables
if(len_trim(cvar).eq.1 .and.cvar(1:1) .eq.'u')                       ngem='UWND'
if(len_trim(cvar).eq.1 .and.cvar(1:1) .eq.'v')                       ngem='VWND'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'u_avg')                   ngem='UWDA' 
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'v_avg')                   ngem='VWDA'
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'ue')                      ngem='UEWD'
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'ve')                      ngem='VEWD'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'ue_avg')                  ngem='UEWA'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'ve_avg')                  ngem='VEWA'
if(len_trim(cvar).eq.1 .and.cvar(1:1) .eq.'w')                       ngem='WWND'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'wcms')                    ngem='WCMS'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'w_avg')                   ngem='WAVG'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'speed')                   ngem='SPED'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'speed_mph')               ngem='SMPH'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'speed10m')                ngem='SP10'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'direction')               ngem='DRCT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'relvortx')                ngem='XVOR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'relvorty')                ngem='YVOR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'relvortz')                ngem='ZVOR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'absvortz')                ngem='AVOR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'potvortz')                ngem='PVOR'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'horiz_div')               ngem='HDIV'

!3D THERMODYNAMIC PROPERTIES OF AIR - 18 variables
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'pi')                      ngem='XNER'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'press')                   ngem='PRES'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pprime')                  ngem='PPRM'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'theta_il')                ngem='THIL'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'theta')                   ngem='THTA'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'dn0')                     ngem='DEN0'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'pi0')                     ngem='XNR0'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'th0')                     ngem='THV0'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'pert_pressure')           ngem='PERT'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'tempk')                   ngem='TMPK'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'tempc')                   ngem='TMPC'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'tempf')                   ngem='TMPF'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'theta_e')                 ngem='THTE'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'theta_v')                 ngem='THTV'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'theta_rho')               ngem='THTR'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'buoyancy_liquid')         ngem='BOYL'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'tempf2m')                 ngem='TMPF2'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'tempc2m')                 ngem='TMPC2'

!3D HYDROMETEOR GAMMA DISTRIBUTION INFO - 26 variables
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'cloud_gam_dm')            ngem='CGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'cloud_gam_d0')            ngem='CGD0'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'rain_gam_dm')             ngem='RGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'rain_gam_d0')             ngem='RGD0'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'rain_gam_lognw')          ngem='RGNW'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'rain_gam_sigma')          ngem='RGSG'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'pris_gam_dm')             ngem='PGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'pris_gam_d0')             ngem='PGD0'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'snow_gam_dm')             ngem='SGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'snow_gam_d0')             ngem='SGD0'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'snow_gam_lognw')          ngem='SGNW'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'snow_gam_sigma')          ngem='SGSG'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'aggr_gam_dm')             ngem='AGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'aggr_gam_d0')             ngem='AGD0'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'aggr_gam_lognw')          ngem='AGNW'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'aggr_gam_sigma')          ngem='AGSG'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'grau_gam_dm')             ngem='GGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'grau_gam_d0')             ngem='GGD0'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'grau_gam_lognw')          ngem='GGNW'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'grau_gam_sigma')          ngem='GGSG'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'hail_gam_dm')             ngem='HGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'hail_gam_d0')             ngem='HGD0'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'hail_gam_lognw')          ngem='HGNW'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'hail_gam_sigma')          ngem='HGSG'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'driz_gam_dm')             ngem='DGDM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'driz_gam_d0')             ngem='DGD0'

!3D MOISTURE MASS MIXING RATIOS AND HUMIDITY - 38 variables
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'vapr_press')              ngem='VPRS'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'rslf')                    ngem='RSLF'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'rsif')                    ngem='RSIF'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'vapor')                   ngem='VMIX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapor_m3')                ngem='VMXV'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'cloud')                   ngem='CMIX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'cloud_m3')                ngem='CMXV'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'rain')                    ngem='RMIX'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'rain_m3')                 ngem='RMXV'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'pristine')                ngem='PMIX'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'pristine_m3')             ngem='PMXV'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'snow')                    ngem='SMIX'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'snow_m3')                 ngem='SMXV'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'aggregates')              ngem='AMIX'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'aggregates_m3')           ngem='AMXV'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'graupel')                 ngem='GMIX'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'graupel_m3')              ngem='GMXV'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'hail')                    ngem='HMIX'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'hail_m3')                 ngem='HMXV'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'drizzle')                 ngem='DMIX'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'drizzle_m3')              ngem='DMXV'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'prissnowagg')             ngem='PSAM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'grauphail')               ngem='GHMX'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'liquid')                  ngem='LMIX'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'ice')                     ngem='IMIX'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'total_cond')              ngem='TMIX'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'total_cond_m3')           ngem='TMXV'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'total_mixr')              ngem='MIXS'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'total_mixr_m3')           ngem='MIXV'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'ctop_tempc_sstbase')      ngem='CTST'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'ctop_tempc_nobase')       ngem='CTOP'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'dewptk')                  ngem='DWPK'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'dewptf')                  ngem='DWPF'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'dewptc')                  ngem='DWPC'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'relhum')                  ngem='RELH'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'relhum_frac')             ngem='RHFR'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'clear_frac')              ngem='CLRF'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'cloud_frac')              ngem='CLDF'

!3D HYDROMETEOR NUMBER CONCENTRATIONS - 22 variables
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'cloud_concen_mg')         ngem='CNMG'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'cloud_concen_kg')         ngem='CNKG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'rain_concen_kg')          ngem='RNKG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'pris_concen_mg')          ngem='PNMG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'pris_concen_kg')          ngem='PNKG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'snow_concen_kg')          ngem='SNKG'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'agg_concen_kg')           ngem='ANKG'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'graup_concen_kg')         ngem='GNKG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'hail_concen_kg')          ngem='HNKG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'drizzle_concen_mg')       ngem='DNMG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'drizzle_concen_kg')       ngem='DNKG'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'cloud_concen_cm3')        ngem='CNC3'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'rain_concen_m3')          ngem='RNM3'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'rain_concen_dm3')         ngem='RND3'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'pris_concen_m3')          ngem='PNM3'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'pris_concen_cm3')         ngem='PNC3'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'snow_concen_m3')          ngem='SNM3'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'snow_concen_cm3')         ngem='SNC3'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'agg_concen_m3')           ngem='ANM3'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'graup_concen_m3')         ngem='GNM3'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'hail_concen_m3')          ngem='HNM3'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'drizzle_concen_cm3')      ngem='DNC3'

!HUCM-SBM SPECIFIC MICROPHYSICS – 18 variables
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'ice_plates')              ngem='IPMX'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'ice_columns')             ngem='ICMX'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'ice_dendrites')           ngem='IDMX'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'plates_concen_mg')        ngem='PCMG'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'plates_concen_kg')        ngem='PCKG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'columns_concen_mg')       ngem='CCMG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'columns_concen_kg')       ngem='CCKG'
if(len_trim(cvar).eq.19.and.cvar(1:19).eq.'dendrites_concen_mg')     ngem='DCMG'
if(len_trim(cvar).eq.19.and.cvar(1:19).eq.'dendrites_concen_kg')     ngem='DCKG'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcpvip')                  ngem='PVIP'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcpvic')                  ngem='PVIC'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcpvid')                  ngem='PVID'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcprip')                  ngem='PRIP'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcpric')                  ngem='PRIC'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pcprid')                  ngem='PRID'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'accpip')                  ngem='ACIP'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'accpic')                  ngem='ACIC'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'accpid')                  ngem='ACID'

!3D AEROSOLS NUMBER, MASS, SIZE, SOLUBILITY - 37 variables
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'ifn_concen_mg')           ngem='IFNM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'ifn_concen_cm3')          ngem='IFNC'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'ccn_concen_mg')           ngem='CCNM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'ccn_concen_cm3')          ngem='CCNC'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'gccn_concen_mg')          ngem='GCNM'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'gccn_concen_cm3')         ngem='GCNC'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'dust1_concen')            ngem='D1CN'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'dust2_concen')            ngem='D2CN'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'salt_film_concen')        ngem='SFCN'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'salt_jet_concen')         ngem='SJCN'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'salt_spume_concen')       ngem='SSCN'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'regen_aero1_concen')      ngem='R1CN'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'regen_aero2_concen')      ngem='R2CN'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'ccn_mass')                ngem='CCCM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'gccn_mass')               ngem='GCCM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust1_mass')              ngem='D1CM'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'dust1_massd10')           ngem='D1CM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust2_mass')              ngem='D2CM'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'dust2_massd10')           ngem='D2CM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'salt_film_mass')          ngem='SFCM'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'salt_jet_mass')           ngem='SJCM'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'salt_spume_mass')         ngem='SSCM'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'regen_aero1_mass')        ngem='R1CM'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'regen_aero2_mass')        ngem='R2CM'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'resol_aero1_mass')        ngem='R1SM'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'resol_aero2_mass')        ngem='R2SM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'regen1_epsilon')          ngem='R1EP'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'regen2_epsilon')          ngem='R2EP'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'ccn_medrad')              ngem='CCCR'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'gccn_medrad')             ngem='GCCR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'dust1_medrad')            ngem='D1CR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'dust2_medrad')            ngem='D2CR'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'salt_film_medrad')        ngem='SFCR'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'salt_jet_medrad')         ngem='SJCR'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'salt_spume_medrad')       ngem='SSCR'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'regen_aero1_medrad')      ngem='R1CR'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'regen_aero2_medrad')      ngem='R2CR'

!3D AEROSOL TRACKING VARIABLES - 41 variables
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'aerosol_cloud_mass')      ngem='ARMC'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_rain_mass')       ngem='ARMR'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_pris_mass')       ngem='ARMP'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_snow_mass')       ngem='ARMS'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_aggr_mass')       ngem='ARMA'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_grau_mass')       ngem='ARMG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_hail_mass')       ngem='ARMH'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'aerosol_driz_mass')       ngem='ARMD'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'aerosol_hydro_mass')      ngem='ARHY'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'soluble_cloud_mass')      ngem='SLMC'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_rain_mass')       ngem='SLMR'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_pris_mass')       ngem='SLMP'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_snow_mass')       ngem='SLMS'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_aggr_mass')       ngem='SLMA'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_grau_mass')       ngem='SLMG'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_hail_mass')       ngem='SLMH'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soluble_driz_mass')       ngem='SLMD'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'soluble_hydro_mass')      ngem='SLHY'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'aero_epsilon')            ngem='EPSI'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'dust_cloud_mass')         ngem='DUMC'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_rain_mass')          ngem='DUMR'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_pris_mass')          ngem='DUMP'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_snow_mass')          ngem='DUMS'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_aggr_mass')          ngem='DUMA'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_grau_mass')          ngem='DUMG'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_hail_mass')          ngem='DUMH'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'dust_driz_mass')          ngem='DUMD'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'dust_hydro_mass')         ngem='DUHY'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'dustifn_cloud_mass')      ngem='DINC'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_rain_mass')       ngem='DINR'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_pris_mass')       ngem='DINP'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_snow_mass')       ngem='DINS'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_aggr_mass')       ngem='DINA'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_grau_mass')       ngem='DING'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_hail_mass')       ngem='DINH'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'dustifn_driz_mass')       ngem='DIND'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'dustifn_hydro_mass')      ngem='DIHY'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'ifn_nuc_numtrack')        ngem='INTR'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'ifn_incloud')             ngem='CICN'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'ifn_indriz')              ngem='DICN'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'ifn_inrain')              ngem='RICN'

!3D VERTICAL VELOSITY AND MICROPHYSICAL INSTANTANEOUS BUDGETS - 15 variables
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'wp_advdif')               ngem='WPAD'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'wp_buoy_theta')           ngem='WPTH'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'wp_buoy_cond')            ngem='WPCD'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'latheatvap')              ngem='LHVP'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'latheatfrz')              ngem='LHFZ'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'nuccldr')                 ngem='NUCR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'cld2rain')                ngem='CL2R'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'ice2rain')                ngem='IC2R'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'nucicer')                 ngem='NUIR'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'vapliq')                  ngem='VAPL'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'vapice')                  ngem='VAPI'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'meltice')                 ngem='MELT'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'rimecld')                 ngem='RIMC'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2ice')                ngem='R2IC'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'aggregate')               ngem='AGGR'

!3D MICROPHYSICAL TOTAL BUDGETS - 88 variables
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'nuccldrt')                ngem='NUCRT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'cld2raint')               ngem='CL2RT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'ice2raint')               ngem='IC2RT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'nucicert')                ngem='NUIRT'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'vapliqt')                 ngem='VAPLT'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'vapicet')                 ngem='VAPIT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'melticet')                ngem='MELTT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rimecldt')                ngem='RIMCT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'rain2icet')               ngem='R2ICT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'aggregatet')              ngem='AGGRT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'latheatvapt')             ngem='LHVPT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'latheatfrzt')             ngem='LHFZT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'inuchomrt')               ngem='IHMRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'inuccontrt')              ngem='ICORT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'inucifnrt')               ngem='IINRT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'inuchazrt')               ngem='IHZRT'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'vapcldt')                 ngem='VAPCT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapraint')                ngem='VAPRT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapprist')                ngem='VAPPT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapsnowt')                ngem='VAPST'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapaggrt')                ngem='VAPAT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapgraut')                ngem='VAPGT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vaphailt')                ngem='VAPHT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vapdrizt')                ngem='VAPDT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'meltprist')               ngem='MELPT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'meltsnowt')               ngem='MELST'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'meltaggrt')               ngem='MELAT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'meltgraut')               ngem='MELGT'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'melthailt')               ngem='MELHT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'rimecldsnowt')            ngem='RIMST'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'rimecldaggrt')            ngem='RIMAT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'rimecldgraut')            ngem='RIMGT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'rimecldhailt')            ngem='RIMHT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2prt')                ngem='R2PRT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2snt')                ngem='R2SNT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2agt')                ngem='R2AGT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2grt')                ngem='R2GRT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rain2hat')                ngem='R2HAT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'aggrselfprist')           ngem='AGPPT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'aggrselfsnowt')           ngem='AGSST'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'aggrprissnowt')           ngem='AGPST'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust1cldrt')              ngem='D1CRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust2cldrt')              ngem='D2CRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust1drzrt')              ngem='D1DRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'dust2drzrt')              ngem='D2DRT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_nuccldrt')             ngem='VNUCRT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_cld2raint')            ngem='VCL2RT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_ice2raint')            ngem='VIC2RT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_nucicert')             ngem='VNUIRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'vt_vapliqt')              ngem='VVAPLT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'vt_vapicet')              ngem='VVAPIT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_melticet')             ngem='VMELTT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rimecldt')             ngem='VRIMCT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_rain2icet')            ngem='VR2ICT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_aggregatet')           ngem='VAGGRT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_inuchomrt')            ngem='VIHMRT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_inuccontrt')           ngem='VICORT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_inucifnrt')            ngem='VIINRT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_inuchazrt')            ngem='VIHZRT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'vt_vapcldt')              ngem='VVAPCT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapraint')             ngem='VVAPRT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapprist')             ngem='VVAPPT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapsnowt')             ngem='VVAPST'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapaggrt')             ngem='VVAPAT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapgraut')             ngem='VVAPGT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vaphailt')             ngem='VVAPHT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_vapdrizt')             ngem='VVAPDT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_meltprist')            ngem='VMELPT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_meltsnowt')            ngem='VMELST'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_meltaggrt')            ngem='VMELAT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_meltgraut')            ngem='VMELGT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vt_melthailt')            ngem='VMELHT'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'vt_rimecldsnowt')         ngem='VRIMST'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'vt_rimecldaggrt')         ngem='VRIMAT'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'vt_rimecldgraut')         ngem='VRIMGT'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'vt_rimecldhailt')         ngem='VRIMHT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rain2prt')             ngem='VR2PRT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rain2snt')             ngem='VR2SNT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rain2agt')             ngem='VR2AGT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rain2grt')             ngem='VR2GRT'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vt_rain2hat')             ngem='VR2HAT'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'vt_aggrselfprist')        ngem='VAGPPT'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'vt_aggrselfsnowt')        ngem='VAGSST'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'vt_aggrprissnowt')        ngem='VAGPST'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_dust1cldrt')           ngem='VD1CRT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_dust2cldrt')           ngem='VD2CRT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_dust1drzrt')           ngem='VD1DRT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vt_dust2drzrt')           ngem='VD2DRT'

!3D HYDROMETEOR DIAMETERS - 9 variables
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'cloudtop_diam')           ngem='TDIM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'cloud_diam')              ngem='CDIM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'rain_diam')               ngem='RDIM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'pris_diam')               ngem='PDIM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'snow_diam')               ngem='SDIM'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'agg_diam')                ngem='ADIM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'graup_diam')              ngem='GDIM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'hail_diam')               ngem='HDIM'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'drizzle_diam')            ngem='DDIM'

!3D HYDROMETEOR TEMPERATURE, ENERGY, LIQUID FRACTION - 11 variables
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'q2')                      ngem='Q2RA'
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'q6')                      ngem='Q6GR'
if(len_trim(cvar).eq.2 .and.cvar(1:2) .eq.'q7')                      ngem='Q7HA'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'rain_temp')               ngem='RTMP'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'graup_temp')              ngem='GTMP'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'hail_temp')               ngem='HTMP'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'rain_air_tempdif')        ngem='RATD'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'graup_air_tempdif')       ngem='GATD'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'hail_air_tempdif')        ngem='HATD'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'graup_fracliq')           ngem='GLIQ'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'hail_fracliq')            ngem='HLIQ'

!3D MISCELLANEOUS FIELDS - 4 variables
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'geo')                     ngem='HGHT'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'tke')                     ngem='TKET'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'pbl_ht')                  ngem='PBLH'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'reflect_all')             ngem='DBZZ'

!CUMULUS PARM - RADIATION - TURBULENCE PARAMETERS - 15 variables
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'cuparm_thetasrc')         ngem='CVHR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'cuparm_rtsrc')            ngem='CVMR'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'khh')                     ngem='KHHC'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'khv')                     ngem='KHVC'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'visibility')              ngem='VISB'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'swup')                    ngem='SWUP'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'swdn')                    ngem='SWDN'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'lwup')                    ngem='LWUP'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'lwdn')                    ngem='LWDN'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'aodt')                    ngem='AODT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'rad_thetasrc')            ngem='RAHR'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'column_net_rad_flx')      ngem='NETR'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'sum_rad_flx')             ngem='NETF'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'sw_heat_rate')            ngem='SWHT'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'lw_heat_rate')            ngem='LWHT'

!2D SURFACE PRECIPITATION and VERTICALLY INTEGRATED FIELDS - 55 variables
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accpr')                   ngem='ACCR'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accpp')                   ngem='ACCP'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accps')                   ngem='ACCS'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accpa')                   ngem='ACCA'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accpg')                   ngem='ACCG'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accph')                   ngem='ACCH'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'accpd')                   ngem='ACCD'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'accpaero')                ngem='ACTA'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'accpdust')                ngem='ACDU'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'dustfrac')                ngem='DFRC'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'totpcp')                  ngem='TRPM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'totpcp_in')               ngem='TRPI'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'precip')                  ngem='TAPM'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'precip_in')               ngem='TAPI'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprr')                   ngem='PCRR'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvr')                   ngem='PCVR'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprp')                   ngem='PCRP'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvp')                   ngem='PCVP'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprs')                   ngem='PCRS'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvs')                   ngem='PCVS'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpra')                   ngem='PCRA'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpva')                   ngem='PCVA'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprg')                   ngem='PCRG'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvg')                   ngem='PCVG'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprh')                   ngem='PCRH'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvh')                   ngem='PCVH'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcprd')                   ngem='PCRD'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'pcpvd')                   ngem='PCVD'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'pcpg')                    ngem='PCPG'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'qpcpg')                   ngem='PCPQ'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'dpcpg')                   ngem='PCPD'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'pcprate')                 ngem='PRRM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'pcprate_in')              ngem='PRRI'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'precipr')                 ngem='PRTM'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'precipr_in')              ngem='PRTI'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'conpcp')                  ngem='CNPR'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'acccon')                  ngem='ACON'

if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'vertavg_w')               ngem='VAVW'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'vertmax_w')               ngem='VMXW'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'vertmin_w')               ngem='VMNW'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vertmax_cloud')           ngem='VMXC'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_rain')            ngem='VMXR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_pris')            ngem='VMXP'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_snow')            ngem='VMXS'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_aggr')            ngem='VMXA'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_grau')            ngem='VMXG'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_hail')            ngem='VMXH'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertmax_driz')            ngem='VMXD'

if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_cond')            ngem='COND'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'vertint_rt')              ngem='WATR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_orig')            ngem='VERT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vertint_vapor')           ngem='VRTV'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vertint_liq')             ngem='VRTL'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'vertint_ice')             ngem='VRTI'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'vertint_cloud')           ngem='VRTC'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_driz')            ngem='VRTD'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_rain')            ngem='VRTR'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_pris')            ngem='VRTP'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_snow')            ngem='VRTS'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_aggr')            ngem='VRTA'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'vertint_graupel')         ngem='VRTG'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_hail')            ngem='VRTH'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'vertint_dust')            ngem='VTDU'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'vertint_dust_hydro')      ngem='VTDH'

!2D SEA ICE COVERAGE, DEPTH, ROUGHNESS, TEMP, SNOW COVER - 5 variables
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'snowdepthonice')          ngem='DEPS'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'cicedepth')               ngem='DEPI'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'cicefract')               ngem='ICEF'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'cicetemp')                ngem='ICET'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'cicerough')               ngem='ICER'

!2D SURFACE HEAT, MOISTURE, MOMENTUM AND RADIATIVE FLUXES - 12 variables
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'sens_flux')               ngem='SFLX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'lat_flux')                ngem='LFLX'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'etrans')                  ngem='EVAP'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'etrans_in')               ngem='ETRI'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'umom_flx')                ngem='UFLX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'vmom_flx')                ngem='VFLX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'wmom_flx')                ngem='WFLX'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'bowen')                   ngem='BOWN'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'rshort')                  ngem='RSHT'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'rlong')                   ngem='RLON'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'rlongup')                 ngem='RLNU'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'albedt')                  ngem='ALBE'

!2D TOPOGRAPHY AND GEOGRAPHIC VALUES - 3 variables
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'topt')                    ngem='TOPT'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'lat')                     ngem='LATI'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'lon')                     ngem='LONG'

!2D MISCELLANEOUS FIELDS - 3 variables
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'sea_press')               ngem='MSLP'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'sfc_div')                 ngem='SDIV'
if(len_trim(cvar).eq.3 .and.cvar(1:3) .eq.'sst')                     ngem='SSTC'

!LEAF/SIB VARIABLES SECTION - 34 variables
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'patch_area')              ngem='PFRA'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'water')                   ngem='OCEN'
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'land')                    ngem='LAND'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'snow_levels')             ngem='SNOL'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'snow_depth_ps')           ngem='SNOD'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'snow_mass_ps')            ngem='SNOM'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'snow_temp_ps')            ngem='SNOT'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'topo_z0_ps')              ngem='TRUF'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'net_z0_ps')               ngem='NRUF'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'soil_z0_ps')              ngem='SRUF'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'veg_z0_ps')               ngem='VRUF'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'veg_ndvi_ps')             ngem='NDVI'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'veg_class_bp')            ngem='VEGC'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'veg_albedo_ps')           ngem='VEGA'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'veg_fracarea_ps')         ngem='VEGF'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'veg_lai_ps')              ngem='LAIF'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'veg_disp_ps')             ngem='VDIS'
if(len_trim(cvar).eq.16.and.cvar(1:16).eq.'canopy_mixrat_ps')        ngem='CANM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'grnd_mixrat_ps')          ngem='GRDM'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'soil_mixrat_ps')          ngem='SOIM'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'veg_moist_ps')            ngem='VEGM'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'veg_temp_ps')             ngem='VEGT'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'canopy_tempc_ps')         ngem='CANC'
if(len_trim(cvar).eq.15.and.cvar(1:15).eq.'canopy_tempf_ps')         ngem='CANF'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'ustar_ps')                ngem='USTR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'tstar_ps')                ngem='TSTR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'rstar_ps')                ngem='RSTR'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'sltex_bp')                ngem='SLTX'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'soilq_ps')                ngem='SOIQ'
if(len_trim(cvar).eq.12.and.cvar(1:12).eq.'soil_temp_ps')            ngem='SOIT'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'soil_moist_ps')           ngem='SLMS'
if(len_trim(cvar).eq.17.and.cvar(1:17).eq.'soil_moistfrac_ps')       ngem='SLMF'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'5050_tempc_ps')           ngem='50TC'
if(len_trim(cvar).eq.13.and.cvar(1:13).eq.'5050_tempf_ps')           ngem='50TF'

!SIB VARIABLES SECTION - 40 variables
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'co2_concen')              ngem='CO2C'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'snow1_ps')                ngem='SNO1'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'snow2_ps')                ngem='SNO2'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'capac1_ps')               ngem='CAP1'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'capac2_ps')               ngem='CAP2'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'pco2ap_ps')               ngem='PCOA'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'co2flx_ps')               ngem='CO2F'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'sfcswa_ps')               ngem='SFAL'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'uplwrf_ps')               ngem='SFUP'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'assimn_ps')               ngem='ASSM'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'respg_ps')                ngem='RESP'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'rstfac1_ps')              ngem='RST1'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'rstfac2_ps')              ngem='RST2'
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'rstfac3_ps')              ngem='RST3'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'ect_ps')                  ngem='ECTF'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'eci_ps')                  ngem='ECIF'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'egi_ps')                  ngem='EGIF'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'egs_ps')                  ngem='EGSF'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'hc_ps')                   ngem='HCFX'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'hg_ps')                   ngem='HGFX'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'ra_ps')                   ngem='RAST'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'rb_ps')                   ngem='RBST'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'rc_ps')                   ngem='RCST'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'rd_ps')                   ngem='RDST'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'roff_ps')                 ngem='ROFF'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'green_ps')                ngem='GREN'
if(len_trim(cvar).eq.7 .and.cvar(1:7) .eq.'apar_ps')                 ngem='APAR'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'ventmf_ps')               ngem='VENT'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'pco2c_ps')                ngem='PCOC'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'pco2i_ps')                ngem='PCOI'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'pco2s_ps')                ngem='PCOS'
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'pco2m_ps')                ngem='PCOM'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'ea_ps')                   ngem='EAPR'
if(len_trim(cvar).eq.5 .and.cvar(1:5) .eq.'em_ps')                   ngem='EMPR'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'rha_ps')                  ngem='RHAC'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'radvbc_ps')               ngem='RVDR'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'radvdc_ps')               ngem='RVDF'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'radnbc_ps')               ngem='RNDR'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'radndc_ps')               ngem='RNDF'
if(len_trim(cvar).eq.6 .and.cvar(1:6) .eq.'psy_ps')                  ngem='PSYC'

!KPP OCEAN MIXED LAYER MODEL VARIABLES - 10 variables
if(len_trim(cvar).eq.8 .and.cvar(1:8) .eq.'kpp_hmix')                ngem='KHMX'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_ocdepth')             ngem='KOCD'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_ust')             ngem='KFUS'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_vst')             ngem='KFVS'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_nsw')             ngem='KNSW'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_nlw')             ngem='KNLW'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_ice')             ngem='KICE'
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'kpp_flx_pcp')             ngem='KPCP'
if(len_trim(cvar).eq.14.and.cvar(1:14).eq.'kpp_depth_temp')          ngem='KDTP'
if(len_trim(cvar).eq.18.and.cvar(1:18).eq.'kpp_depth_salinity')      ngem='KDSL'

!TRACER FIELDS SECTION - 6 variables
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer001')               ngem='T001'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer002')               ngem='T002'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer003')               ngem='T003'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer004')               ngem='T004'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer005')               ngem='T005'
if(len_trim(cvar).eq.9 .and.cvar(1:9) .eq.'tracer006')               ngem='T006'


!SET SPECIAL FLAG FOR VARIABLES OUTPUT ON ZM LEVELS
if(len_trim(cvar).eq.10.and.cvar(1:10).eq.'visibility')              zl=2 
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'swup')                    zl=2 
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'swdn')                    zl=2 
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'lwup')                    zl=2 
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'lwdn')                    zl=2 
if(len_trim(cvar).eq.11.and.cvar(1:11).eq.'sum_rad_flx')             zl=2 
if(len_trim(cvar).eq.1 .and.cvar(1:1) .eq.'w')                       zl=2     
if(len_trim(cvar).eq.4 .and.cvar(1:4) .eq.'wcms')                    zl=2 

!*****************************************************************************
!**** Begin extracting the data **********************************************
!*****************************************************************************
if(ivtype.ge.2) then

 if(iztrans.eq.1.or.iztrans.eq.2) coordn='HGHT' !Sigma or cartesian coordinate
 if(iztrans.eq.3) coordn='PRES' !Pressure coordinate

 ! If this is a 3-dimensional atmospheric variable      
 do k=nnb,nne,nninc

  !2D variables
  if(ivtype==2) then
    typlev=0
    coordn='NONE'
  endif
  !LEAF surface model variables
  if(ivtype==4 .or. ivtype==5) then
    typlev=k
    coordn='NONE' !outputs over number of soil or snow layers
  endif
  !3D variables
  if(ivtype==3) then
    if(iztrans.eq.1.or.iztrans.eq.2) then
      if(zl==1) typlev=int(ztn(k,nngd)) ! T-vertical stagger (scalars)
      if(zl==1 .and. typlev==0) typlev=1 !do not allow zt level = 0m
      if(zl==2) typlev=int(zmn(k,nngd)) ! M-vertical stagger (W,radfluxes)
    endif
    if(iztrans.eq.3) typlev=int(iplevs(k))
  endif
  !Patch variables
  if(ivtype==6) then
    typlev=k
    coordn='NONE' !outputs over number of patches
  endif

  !Advance the number of grid output thus far
  numgrids=numgrids+1

  !Write gempak grid header information
  write(iun,*)
  write(iun,*)
  write(iun,'(a12,3i2.2,i4.4,a9)') &
    ,' Grid file: ',iyr1,imonth1,idate1,itime1,'_rams.gem'
  write(iun,'(a17)'),' GRID IDENTIFIER:'
  write(iun,'(a60)') &
    ,'    TIME1             TIME2         LEVL1 LEVL2   VCORD PARM'
  write(iun,'(3i2.2,a1,2i2.2,24X,i6,10X,a4,1X,a5)') &
    ,iyr2,imonth2,idate2,'/',ihour2,imin2,typlev,coordn,ngem
  write(iun,'(a6,f7.2,a1,f8.2,a1,f7.2,a1,f8.2,18X,a13,2i5)') &
    ,' AREA:',swlat1,';',swlon1,';',nelat1,';',nelon1 &
    ,'GRID SIZE: ',xpts,ypts

!##################################################################
!DIFFERENT WAYS TO WRITE THE COLUMNS
  write(iun,'(a17,i3,a18,i3)') &
    ,' COLUMNS:     1  ',xpts,'     ROWS:     1  ',ypts
!  write(iun,'(a12,i3,a2,i3,a13,i3,a2,i3)') &
!    ,' COLUMNS:   ',xbeg,'  ',xend,'     ROWS:   ',ybeg,'  ',yend
!##################################################################

  write(iun,*)
  write(iun,'(a21)') ' Scale factor: 10** 0'
  write(iun,*)
  write(iun,*)
  write(iun,'(a8,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7)') &
    ' COLUMN:',(i,i=1,8)
  write(iun,'(8X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7,2X,i7)') &
    (i,i=9,xpts)

  !Determine how to write out arrays depending on vtype/dimension
  !and if simulation is 2D or 3D
  if(ivtype.eq.2) kk=1
  if(ivtype.ge.3) kk=k
  jj=1

  jcount=ypts
  do j=yend,ybeg,-yinc
    if(nje>1)jj=j
    write(iun,'(a4,i3,a4,8e14.6)') &
        ' ROW',jcount,'    ',(a(i,jj,kk),i=xbeg,(7*xinc+xbeg),xinc)
    write(iun,'(11X,8e14.6)') (a(i,jj,kk),i=(7*xinc+xbeg+xinc),xend,xinc)
    jcount=jcount-1
  enddo

 enddo !loop k-levels

 write(iuntag,'(a9,i5,a1)'),'NUMGRDS="',numgrids,'"'
 backspace iuntag

else
   print*,'unknown ivtype',ivtype
   stop 'RAMS_text'
endif

return
END SUBROUTINE rams_text
