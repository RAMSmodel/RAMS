!##############################################################################
Subroutine latlon_ps (n1,n2,np1,np2,np3,nvar,ipos,aps,psglob,all  &
     ,glat,glon,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)

implicit none

integer :: n1,n2,np1,np2,np3,nvar,ipos,idatelin,iglobew,iglobs,iglobn
real :: xswlat,xswlon,gdatdx,gdatdy
real :: aps(n1,n2,np3,*),psglob(np1+3,np2+2),all(np1,np2,np3,*)  &
     ,glon(n1,n2),glat(n1,n2)
real :: fioff,fjoff,grx,gry,gglon
integer :: np1h,ioff,joff,nv,i,j,k

np1h = np1 / 2

! Determine if psglob scratch array needs to be offset to west and/or south.

ioff = 0
joff = 0

if (iglobew .eq. 1) then
   ioff = 1
   if (iglobs .eq. 1) then
      joff = 1
   endif
endif
fioff = float(ioff)
fjoff = float(joff)

! Loop through all variable types and vertical levels

do nv = 1,nvar
   do k = 1,np3

! Fill main part of scratch array psglob from data

      do j = 1,np2
         do i = 1,np1
            psglob(i+ioff,j+joff) = all(i,j,k,nv)
         enddo
      enddo

! If all longitudes exist in input data, add 1 row on west boundary and 2
! rows on east boundary of input data in scratch array.

      if (iglobew .eq. 1) then
         do j = 1,np2
            psglob(    1,j+joff) = psglob(np1+1,j+joff)
            psglob(np1+2,j+joff) = psglob(    2,j+joff)
            psglob(np1+3,j+joff) = psglob(    3,j+joff)
          enddo
      endif

! If input data includes the south pole, add 1 row to south in scratch array.

      if (iglobs .eq. 1) then
         do i = 1,np1h
            psglob(i,1) = psglob(i+np1h,3)
         enddo
         do i = np1h+1,np1+3
            psglob(i,1) = psglob(i-np1h,3)
         enddo
      endif

! If input data includes the north pole, add 1 row to north in scratch array.

      if (iglobn .eq. 1) then
         do i = 1,np1h
            psglob(i,np2+2) = psglob(i+np1h,np2)
         enddo
         do i = np1h+1,np1+3
            psglob(i,np2+2) = psglob(i-np1h,np2)
         enddo
      endif

! Loop through all i,j points on model grid

      do j = 1,n2
         do i = 1,n1
            gry = (glat(i,j) - xswlat) / gdatdy + 1. + fjoff

            gglon = glon(i,j)
            if (idatelin .eq. 1 .and. gglon .lt. xswlon)  &
               gglon = glon(i,j) + 360.
            grx = (gglon - xswlon) / gdatdx + 1. + fioff

!     horizontally interpolate variable

            CALL gdtost (psglob,np1+3,np2+2,grx,gry,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
END SUBROUTINE latlon_ps

!##############################################################################
Subroutine lambcon_ps (n1,n2,np1,np2,n3,nvar,ipos,aps,all,glat,glon  &
     ,xswlat,xswlon,gdatdx,cntlat,cntlon)

!lambert-conformal to polar stereographic

implicit none

integer :: n1,n2,np1,np2,n3,nvar,ipos
real :: xswlat,xswlon,gdatdx,cntlat,cntlon
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)
real :: spacing,ri,rj
integer :: i,j,k,nv

!c     real spacing
!c     spacing=  ! (meters, for AWIPS grid 212 this should be 40635.25)
!c
  spacing = gdatdx 

do i = 1, n1
   do j = 1, n2


      CALL ll2lc (cntlat, cntlon, glat(i,j), glon(i,j), ri, rj  &
         ,xswlat,xswlon, spacing)


      if(ri.lt.1.or.ri.gt.np1.or.rj.lt.1.or.rj.gt.np2)  &
           then
         print*,' gridpoint outside of lambert-con data coverage'
         print*,'   polar i,j,lat,lon-',i,j,glat(i,j),glon(i,j)
         print*,'   computed datafile-relative point:',ri,rj
         stop 'vi-ob'
      endif
!
!     horizontally interpolate variables
!
      do nv=1,nvar
         do k=1,n3
            CALL gdtost (all(1,1,k,nv),np1,np2  &
                 ,ri,rj,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
END SUBROUTINE lambcon_ps

!##############################################################################
Subroutine ll2lc (dcentlat,dcentlon,dinlat,dinlon,ri,rj & 
                 ,doriglat,doriglon,spacing)

!routine to convert lat/lon pair to real i,j pair on Lambert
!conformal Conic grid

implicit none

!     ***** Passed variables ********

!     dcentlat:tangent latitude (deg)
!     dcentlon:longitude of y axis on lcc (deg)
!     dinlat:latitude of input point (deg)
!     dinlon:longitude of input point (deg)
!     ri:place for output i-index (altered)
!     rj:place for output j-index (altered)
!     doriglat: latitude of (1,1) point on lcc grid (deg)
!     doriglon: longitude of (1,1) point on lcc grid (deg)
!     spacing: grid spacing on lcc grid (meters)

real dcentlat, dcentlon, dinlat, dinlon, ri, rj, doriglat,  &
     doriglon, spacing

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     r*: radian equivalents of d* above
!     others are scratch variables (see code)

real cone, rcentlat, rcentlon, rinlat, rinlon, roriglat,  &
     roriglon, xorig, yorig, rho, theta, Ffactor

!     ***** Parameters **************

!     radearth: mean radius of earth (m)

real radearth, pi

parameter (radearth=6370997) ! (USGS Paper 1395)
parameter (pi=3.14159265359)

!     ***** Code ********************

!     Convert degrees to radians

rcentlat=dcentlat*pi/180.0
rcentlon=dcentlon*pi/180.0
rinlat=dinlat*pi/180.0
rinlon=dinlon*pi/180.0
roriglat=doriglat*pi/180.0
roriglon=doriglon*pi/180.0

!     Compute projection constants

cone=sin(rcentlat)
Ffactor=cos(rcentlat)*tan(pi/4.0 + rcentlat/2.0)**cone/cone


!     Compute location of origin of lcc grid (1,1)

rho= Ffactor/tan(pi/4.0 + roriglat/2.0)**cone
theta= cone*(roriglon-rcentlon)
xorig= rho*sin(theta)
yorig= -1.0*rho*cos(theta)


!     Compute location of input lat and lon

rho= Ffactor/tan(pi/4.0 + rinlat/2.0)**cone
theta= cone*(rinlon-rcentlon)
ri= rho*sin(theta)
rj= -1.0*rho*cos(theta)

!     Compute location relative to origin

ri= ri - xorig
rj= rj - yorig

!     Convert to index

ri= ri*radearth/spacing + 1.0
rj= rj*radearth/spacing + 1.0

!     Done

return
END SUBROUTINE ll2lc

!##############################################################################
Subroutine lc2ll (dcentlat,dcentlon,qlat,qlon,ri,rj,doriglat,doriglon,spacing)

!routine to convert real i,j pair on Lambert
!conformal Conic grid to lat/lon

implicit none

!     ***** Passed variables ********

!     dcentlat:tangent latitude (deg)
!     dcentlon:longitude of y axis on lcc (deg)
!     qlat:output latitude (deg)
!     qlon:output longitude (deg)
!     ri:input i-index
!     rj:input j-index
!     doriglat: latitude of (1,1) point on lcc grid (deg)
!     doriglon: longitude of (1,1) point on lcc grid (deg)
!     spacing: grid spacing on lcc grid (meters)

real dcentlat, dcentlon, qlat, qlon, ri, rj, doriglat,  &
     doriglon, spacing

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     r*: radian equivalents of d* above
!     others are scratch variables (see code)

real cone, rcentlat, rcentlon, roriglat,  &
     roriglon, xorig, yorig, rho, theta, Ffactor, rho0

!     ***** Parameters **************

!     radearth: mean radius of earth (m)

real radearth, pi

parameter (radearth=6370997) ! (USGS Paper 1395)
parameter (pi=3.14159265359)

!     ***** Code ********************

!     Convert degrees to radians

rcentlat=dcentlat*pi/180.0
rcentlon=dcentlon*pi/180.0
roriglat=doriglat*pi/180.0
roriglon=doriglon*pi/180.0

!     Compute projection constants

cone=sin(rcentlat)
Ffactor=cos(rcentlat)*(tan(pi/4.0 + rcentlat/2.0))**cone/cone

!print*, 'cone,ffactor',cone,ffactor


!     Compute location of origin of lcc grid (1,1)

rho0= Ffactor/(tan(pi/4.0 + rcentlat/2.0))**cone
rho = Ffactor/(tan(pi/4.0 + roriglat/2.0))**cone
theta= cone*(roriglon-rcentlon)
xorig= rho*sin(theta)
yorig= rho0 - rho*cos(theta)

!print*, 'rho,theta,xorig,yorig',rho,theta,xorig,yorig


! convert ri,rj to radians

ri = (ri - 1.) * spacing / radearth
rj = (rj - 1.) * spacing / radearth

! get ri,rj relative to rcentlat,rcentlon

ri = ri + xorig
rj = rj + yorig

!     Compute location of input lat and lon

rho = sqrt(ri ** 2 + (rho0 - rj) ** 2)
theta = atan2(ri,(rho0 - rj))
qlon = theta / cone + rcentlon
qlat = 2. * atan((Ffactor/rho) ** (1./cone)) - pi / 2.

! convert to degrees
 
qlat = qlat * 180. / pi
qlon = qlon * 180. / pi

return
END SUBROUTINE lc2ll

!##############################################################################
Subroutine rotate_winds (rot_type,n1,n2,n3,u1,v1,u2,v2  &
                       ,qlat,qlon,reflat1,reflon1,reflat2,reflon2,idatain)

!    Rotate 3-D wind field depending on value of rot_type
!       u1,v1           :  input winds
!       reflat1,reflon1 :  reference lat-lon of input wind projection
!       u2,v2           :  output winds
!       reflat2,reflon2 :  reference lat-lon of output wind projection
!       qlat, qlon      :  lat-lon arrays at grid points
                       
implicit none

character(len=*) :: rot_type
integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: u1,v1,u2,v2
real, dimension(n1,n2) :: qlat,qlon
real :: reflat1,reflon1,reflat2,reflon2

real :: uu,vv
integer :: i,j,k,idatain

print*,'in rot:',rot_type
print*,'in rot:',n1,n2,n3,reflat1,reflon1,reflat2,reflon2
if(n1 > 1 .and. n2 > 1)print*,'in rot:',qlat(n1/2,n2/2),qlon(n1/2,n2/2)

print*,'DATA INGEST(idatain): ',idatain

if (rot_type == 'll_rps') then
   !   lat-lon (earth-relative) winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            CALL uevetouv (u2(i,j,k),v2(i,j,k),u1(i,j,k),v1(i,j,k)  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
         enddo
      enddo
   enddo
elseif (rot_type == 'rps_ll') then
   !   RAMS polar-stereo  to lat-lon (earth-relative) winds
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            CALL uvtoueve (u1(i,j,k),v1(i,j,k),u2(i,j,k),v2(i,j,k)  &
                         ,qlat(i,j),qlon(i,j),reflat1,reflon1)
         enddo
      enddo
   enddo
elseif (rot_type == 'lc_rps') then
   !   lambert-conformal-relative winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            !Saleeby(2010): NARR is on lambert-conformal grid, but the winds
            !are already rotated to earth-relative. So we want to go straight
            !from earth-relative to RAMS rps projection.
            if(idatain == 1) then
              CALL uevetouv (u2(i,j,k),v2(i,j,k),u1(i,j,k),v1(i,j,k)  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
            else
              CALL uvlc_uvll (u1(i,j,k),v1(i,j,k),uu,vv  &
                          ,qlat(i,j),qlon(i,j),reflat1,reflon1)
              CALL uevetouv (u2(i,j,k),v2(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
            endif
            if(i==n1/2.and.j==n2/2.and.k==8) then
               print*,'input:',u1(i,j,k),v1(i,j,k),qlat(i,j),qlon(i,j)
               print*,'earth:',uu,vv,reflat1,reflon1
               print*,'ps:',u2(i,j,k),v2(i,j,k),reflat2,reflon2
            endif
         enddo
      enddo
   enddo
elseif (rot_type == 'ps_rps') then
   !   polar-stereo-relative winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            CALL uvtoueve (u1(i,j,k),v1(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat1,reflon1)
            CALL uevetouv (u2(i,j,k),v2(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
         enddo
      enddo
   enddo
endif

return
END SUBROUTINE rotate_winds

!##############################################################################
Subroutine uvlc_uvll (ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon)

!    Rotate lambert-conformal-relative winds to earth-relative components

implicit none

real :: ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon,x0,y0,x1,y1,angle

CALL ll_rotate_lc (qlat,qlon   ,cntlat,cntlon,x0,y0)
CALL ll_rotate_lc (qlat,qlon+.1,cntlat,cntlon,x1,y1)

angle = -atan2(y1-y0,x1-x0)
ull = ulc * cos(angle) - vlc * sin(angle)
vll = ulc * sin(angle) + vlc * cos(angle)

return
END SUBROUTINE uvlc_uvll

!##############################################################################
Subroutine ll_rotate_lc (qlat,qlon,centlat,centlon,x,y)

!routine to convert qlat,qlon to x,y (relative to centlat/centlon)
!on Lambert conformal Conic grid

!     ***** Passed variables ********

!     qlat    : latitude of input point (deg)
!     qlon    : longitude of input point (deg)
!     centlat : tangent latitude (deg)
!     centlon : longitude of y axis on lcc (deg)
!     x       : x-coordinate of point on Lambert Conformal grid
!     y       : y-coordinate of point on Lambert Conformal grid

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     f
!     others are scratch variables (see code)

implicit none

real :: qlat,qlon,centlat,centlon,x,y

real, parameter :: radearth = 6370997           ! (USGS Paper 1395)
real, parameter :: pi180 = 3.14159265359 / 180.

real :: rcentlat,rcentlon,rqlat,rqlon,cone,Ffactor,rho,rho8999,theta

rcentlat = centlat * pi180
rcentlon = centlon * pi180
rqlat = qlat * pi180
rqlon = qlon * pi180

cone = sin(rcentlat)
Ffactor = cos(rcentlat) * tan(45. * pi180 + rcentlat/2.) ** cone / cone

if (qlat .le. 89.99) then
   rho = Ffactor/tan(45. * pi180 + rqlat / 2.) ** cone
else
   rho8999 = Ffactor/tan(45. * pi180 + 89.99 * pi180 / 2.) ** cone
   rho = 100. * (90. - qlat) * rho8999
endif
   
theta = cone * (rqlon - rcentlon)
x = rho * sin(theta)
y = - rho * cos(theta)

x = x * radearth
y = y * radearth

return
END SUBROUTINE ll_rotate_lc

!##############################################################################
Subroutine trueps60_ps (n1,n2,np1,np2,n3,nvar,ipos,aps,all,glat,glon  &
                       ,xswlat,xswlon,gdatdx,cntlon)

implicit none

integer :: n1,n2,np1,np2,n3,nvar,ipos
real :: xswlat,xswlon,gdatdx,cntlon
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)
real :: spacing,ri,rj
integer :: i,j,k,nv
     
! True polar stereographic to RAMS rotated polar stereographic

!-----------------------------------------------------------------------
! This code could be used if we made new ll_xy and xy_ll routines
!   that used NCEP constants, mainly
!   DATA  RERTH /6.3712E+6/,PI/3.1416/
!   DATA  SS60 /1.86603/

!  Call ll_xy (xswlat,xswlon,cntlat,cntlon,xw,ys)
!  Convert distance true at 60N to distance true at polepoint
!  fact=0.5*(1.+sqrt(3.)/2.)
!  xe=xw+float(np1-1)*splamx/fact
!  yn=ys+float(np2-1)*splamy/fact
!  print*,'xw,xe,ys,yn=',xw,xe,ys,yn
!  Call xy_ll (xnelat,xnelon,cntlat,cntlon,xe,yn)
!  print*,'SW,NE corner=',xswlat,xswlon,xnelat,xnelon
!
!  Call xy_ll (xnwlat,xnwlon,cntlat,cntlon,xw,yn)
!  Call xy_ll (xselat,xselon,cntlat,cntlon,xe,ys)
!  print*,'NW,SE corner=',xnwlat,xnwlon,xselat,xselon
!-----------------------------------------------------------------------

!(meters, for AWIPS grid 212 this should be 40635.25)
spacing = gdatdx

do i = 1, n1
   do j = 1, n2

      ! RAMS lat/lon coordinates on the other polar stereo coord
      CALL w3fb06 (glat(i,j),glon(i,j),xswlat,xswlon,spacing,cntlon,ri,rj)

      if(ri.lt.1.or.ri.gt.np1.or.rj.lt.1.or.rj.gt.np2) then
         print*,' gridpoint outside true Polr Ster data coverage'
         print*,'   polar i,j,lat,lon-',i,j,glat(i,j),glon(i,j)
         print*,'   computed datafile-relative point:',ri,rj
         stop 'vi-ob in sub trueps60_ps in asti.f90'
      endif

     ! horizontally interpolate variables
      do nv=1,nvar
         do k=1,n3
            CALL gdtost (all(1,1,k,nv),np1,np2,ri,rj,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
END SUBROUTINE trueps60_ps

!##############################################################################
Subroutine w3fb07 (XI,XJ,ALAT1,ALON1,DX,ALONV,ALAT,ALON)

implicit none
 
! ROUTINE:  W3FB07        GRID COORDS TO LAT/LON FOR GRIB
!   PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-04-05
!
! ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION ON EARTH GIVEN IN A
!   GRID COORDINATE OVERLAID ON A POLAR STEREOGRAPHIC MAP PRO-
!   JECTION TRUE AT 60 DEGREES N OR S LATITUDE TO THE
!   NATURAL COORDINATE OF LATITUDE/LONGITUDE
!   W3FB07 IS THE REVERSE OF W3FB06.
!   USES GRIB SPECIFICATION OF THE LOCATION OF THE GRID
!
! HISTORY LOG:
!   88-01-01  ORIGINAL AUTHOR:  STACKPOLE, W/NMC42
!
! USAGE:  W3FB07(XI,XJ,ALAT1,ALON1,DX,ALONV,ALAT,ALON)
!   INPUT ARGUMENT LIST:
!     XI       - I COORDINATE OF THE POINT  REAL*4
!     XJ       - J COORDINATE OF THE POINT  REAL*4
!     ALAT1    - LATITUDE  OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                LATITUDE <0 FOR SOUTHERN HEMISPHERE; REAL*4
!     ALON1    - LONGITUDE OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                  EAST LONGITUDE USED THROUGHOUT; REAL*4
!     DX       - MESH LENGTH OF GRID IN METERS AT 60 DEG LAT
!                 MUST BE SET NEGATIVE IF USING
!                 SOUTHERN HEMISPHERE PROJECTION; REAL*4
!                   190500.0 LFM GRID,
!                   381000.0 NH PE GRID, -381000.0 SH PE GRID, ETC.
!     ALONV    - THE ORIENTATION OF THE GRID.  I.E.,
!                THE EAST LONGITUDE VALUE OF THE VERTICAL MERIDIAN
!                WHICH IS PARALLEL TO THE Y-AXIS (OR COLUMNS OF
!                THE GRID) ALONG WHICH LATITUDE INCREASES AS
!                THE Y-COORDINATE INCREASES.  REAL*4
!                   FOR EXAMPLE:
!                   255.0 FOR LFM GRID,
!                   280.0 NH PE GRID, 100.0 SH PE GRID, ETC.
!
!   OUTPUT ARGUMENT LIST:
!     ALAT     - LATITUDE IN DEGREES (NEGATIVE IN SOUTHERN HEMI.)
!     ALON     - EAST LONGITUDE IN DEGREES, REAL*4
!
!   REMARKS: FORMULAE AND NOTATION LOOSELY BASED ON HOKE, HAYES,
!     AND RENNINGER'S "MAP PROJECTIONS AND GRID ...", MARCH 1981
!     AFGWC/TN-79/003
!
! ATTRIBUTES:
!   LANGUAGE: IBM VS FORTRAN
!   MACHINE:  NAS

real, parameter :: RERTH=6.3712E+6,PI=3.1416, SS60=1.86603
real ::xi,xj,alat1,alon1,dx,alonv,alat,alon,h,dxl,reflon,radpd,degprd &
      ,rebydx,ala1,rmll,alo1,polei,polej,xx,yy,r2,gi2,arccos

! PRELIMINARY VARIABLES AND REDIFINITIONS

! H = 1 FOR NORTHERN HEMISPHERE; = -1 FOR SOUTHERN

! REFLON IS LONGITUDE UPON WHICH THE POSITIVE X-COORDINATE
! DRAWN THROUGH THE POLE AND TO THE RIGHT LIES
! ROTATED AROUND FROM ORIENTATION (Y-COORDINATE) LONGITUDE
! DIFFERENTLY IN EACH HEMISPHERE

IF(DX.LT.0) THEN
   H = -1.
   DXL = -DX
   REFLON = ALONV - 90.
ELSE
   H = 1.
   DXL = DX
   REFLON = ALONV - 270.
ENDIF

RADPD = PI/180.0
DEGPRD = 180./PI
REBYDX = RERTH/DXL

! RADIUS TO LOWER LEFT HAND (LL) CORNER

ALA1 =  ALAT1 * RADPD
RMLL = REBYDX * COS(ALA1) * SS60/(1. + H * SIN(ALA1))

! USE LL POINT INFO TO LOCATE POLE POINT

ALO1 = (ALON1 - REFLON) * RADPD
POLEI = 1. - RMLL * COS(ALO1)
POLEJ = 1. - H * RMLL * SIN(ALO1)

! RADIUS TO THE I,J POINT (IN GRID UNITS)

XX =  XI - POLEI
YY = (XJ - POLEJ) * H
R2 =  XX**2 + YY**2

! NOW THE MAGIC FORMULAE

IF(R2.EQ.0) THEN
   ALAT = H * 90.
   ALON = REFLON
ELSE
   GI2 = (REBYDX * SS60)**2
   ALAT = DEGPRD * H * ASIN((GI2 - R2)/(GI2 + R2))
   ARCCOS = ACOS(XX/SQRT(R2))
   IF(YY.GT.0) THEN
      ALON = REFLON + DEGPRD * ARCCOS
   ELSE
      ALON = REFLON - DEGPRD * ARCCOS
   ENDIF
ENDIF
!IF(ALON.LT.0) ALON = ALON + 360.

return
END SUBROUTINE w3fb07

!##############################################################################
Subroutine w3fb06 (ALAT,ALON,ALAT1,ALON1,DX,ALONV,XI,XJ)

implicit none

real :: ALAT,ALON,ALAT1,ALON1,DX,ALONV,XI,XJ
real :: h,reflon,radpd,rebydx,ala1,alo1,polei,polej,ala,rm,alo,dxl,rmll

! ROUTINE:  W3FB06        LAT/LON TO POLA (I,J) FOR GRIB
!   PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-04-05

real, parameter :: RERTH=6.3712E+6,PI=3.1416, SS60=1.86603

IF(DX.LT.0) THEN
   H = -1.
   DXL = -DX
   REFLON = ALONV - 90.
ELSE
   H = 1.
   DXL = DX
   REFLON = ALONV - 270.
ENDIF

RADPD = PI/180.0
REBYDX = RERTH/DXL

ALA1 =  ALAT1 * RADPD
RMLL = REBYDX * COS(ALA1) * SS60/(1. + H * SIN(ALA1))

ALO1 = (ALON1 - REFLON) * RADPD
POLEI = 1. - RMLL * COS(ALO1)
POLEJ = 1. - H * RMLL * SIN(ALO1)

ALA =  ALAT * RADPD
RM = REBYDX * COS(ALA) * SS60/(1. + H * SIN(ALA))

ALO = (ALON - REFLON) * RADPD
XI = POLEI + RM * COS(ALO)
XJ = POLEJ + H * RM * SIN(ALO)

return
END SUBROUTINE w3fb06
