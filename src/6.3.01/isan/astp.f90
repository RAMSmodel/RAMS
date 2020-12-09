!##############################################################################
Subroutine pressure_stage (n1,n2,glat,glon)

use isan_coms
use rconstants

implicit none

integer :: n1,n2
real :: glat(n1,n2),glon(n1,n2)
real :: fnprx,fnpry,grx,gry,gglat,gglon,thmax,thmin  &
       ,xswlon_east,cntlon_east,rr,rhmax,tinc
integer :: i,j,k,lv,n,iunit,iyyo,immo,iddo,ihho
real, external :: rsatmix

! Read the header of input pressure file.

iunit=11

print 91,trim(innpr)
91 format(//,' Reading pressure gridded data',/,a,//)

open(iunit,file=innpr)
print*,'Doing RALPH 2 format'
read(iunit,*) iyy,imm,idd,ihh,itinc,nprz,nprx,npry
read(iunit,*) inproj,gdatdx,gdatdy,xswlat,xswlon  &
             ,xnelat,xnelon,cntlat,cntlon
read(iunit,*) (levpr(lv),lv=1,nprz)

!Saleeby(2010): Set itinc=0 since this variable seems un-necessary
! and interferes with making of varfiles
!itinc=0

! Check for consistency between file parameters and processing time
tinc = (itinc/100 * 60) + mod(itinc,100)
CALL date_add_to (iyy,imm,idd,ihh*100,tinc, 'm',iyyo,immo,iddo,ihho)
ihho = ihho/100

if(iyyo /= iyear .or. immo /= imonth  &
     .or. iddo /= idate .or. ihho /= ihour) then
   print*,'Pressure file dates not the same as namelist!'
   print*,'Year :',iyyo,iyear
   print*,'Month:',immo,imonth
   print*,'Day  :',iddo,idate
   print*,'Hour :',ihho,ihour
   stop 'pr_dates'
endif

!Sean Freeman 2/20/19
!if we have a subset of lat lon data, we do not need to do any of the extra.
!check if we do not span the whole globe.
if(inproj.eq.1 .and. abs(xnelat-xswlat)<179 .and. abs(xswlon-xnelon)<357) then
  ! If necessary, convert longitude specification to values in the range
  ! [-180.,180.].

   xswlon=mod(xswlon+900.,360.)-180.
   xnelon=mod(xnelon-900.,-360.)+180.
   iglobew=0
   iglobn=0
   iglobs=0
   idatelin =0


!if we are latlon and span the globe
elseif (inproj.eq.1) then

  ! If necessary, convert longitude specification to values in the range
  ! [-180.,180.].

   xswlon=mod(xswlon+900.,360.)-180.
   xnelon=mod(xnelon-900.,-360.)+180.

   if(xswlon.lt.-180..or.xswlon.ge.180.  .or.  &
      xnelon.lt.-180..or.xnelon.gt.180.01.or.  &
      xswlat.lt.-90. .or.xnelat.gt.90.01) then
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!! BAD DOMAIN SPECIFICATION   !!'
      print*,'!!! xswlat,xswlon - ',xswlat,xswlon
      print*,'!!! xnelat,xnelon-xswlon - ',xnelat,xnelon-xswlon
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop 'bad_domain'
   endif

   fnprx=float(nprx)
   fnpry=float(npry)

   ! Set "global domain" flags to determine whether 2 extra rows will be
   ! added to borders of the input gridded data.  iglobew = 1 if doing
   ! full 360 degrees of longitude.  iglobs = 1 if including the south
   ! pole and if iglobew = 1, and iglobn = 1 if including
   ! the north pole and if iglobew = 1.

   iglobew=0
   iglobn=0
   iglobs=0

   if(fnprx*gdatdx.gt.359.9) then
      if(abs(xswlon+180.).gt.1.e-3) then
         print*,'When extracting all longitudes of pressure data,'
         print*,'xswlon must be set to -180.'
         stop 'xswlon'
      endif
      iglobew=1
   endif
   if(iglobew.eq.1.and.xswlat.lt.-89.9) iglobs=1
   if(iglobew.eq.1.and.xnelat.gt. 89.9) iglobn=1
   idatelin=0
   if(xnelon.lt.xswlon+.01) idatelin=1

   ! Search coarse grid points for any that are outside the bounds of
   ! the pressure data.  If any are found, stop.
   ! Saleeby(2014): May need to taylor this section to your lat/lon dateset

   do j = 1,n2
     do i = 1,n1
            gglat=glat(i,j)
            gglon=glon(i,j)
            gry = (gglat - xswlat) / gdatdy + 1.
            if(gry.lt.2-iglobs.or.gry.gt.fnpry+iglobn-1) then
               print*,'Model grid point latitude must be'
               print*,'at least 1 gridpoint inside'
               print*,'pressure data area'
               print*,'i,j,glat,xswlat,xnelat'  &
                     ,i,j,gglat,xswlat,xnelat
               stop 'isan: outside p-data 1'
            endif
            if(idatelin.eq.1.and.gglon.lt.xswlon) gglon=gglon+360.
            grx=(gglon-xswlon)/gdatdx+1.
            if(grx.lt.2-iglobew.or.grx.gt.fnprx+iglobew-1) then
             !Saleeby(2014):May have to fix this for your global lat/lon data
             if( abs(xswlon)-xnelon <= gdatdx .and. &
                 abs(xswlon) >= gglon .and. gglon >= xnelon ) then
               print*,'Model grid point longitude outside of pressure'
               print*,'data area, but this is an artificial artifact.'
               print*,'i,j,gglon,xswlon,xnelon,gdatdx,iglobew,iglobs,iglobn'  &
                      ,i,j,gglon,xswlon,xnelon,gdatdx,iglobew,iglobs,iglobn
             else
               print*,'Model grid point longitude must be'
               print*,'at least 1 gridpoint inside'
               print*,'pressure data area'
               print*,'i,j,gglon,xswlon,xnelon,gdatdx,iglobew,iglobs,iglobn'  &
                      ,i,j,gglon,xswlon,xnelon,gdatdx,iglobew,iglobs,iglobn
               stop 'isan: outside p-data 2'
             endif
            endif
     enddo
   enddo

endif !inproj==1

! Deallocate memory for the pressure data 

if(allocated(p_u))   deallocate(p_u)
if(allocated(p_v))   deallocate(p_v)
if(allocated(p_t))   deallocate(p_t)
if(allocated(p_z))   deallocate(p_z)
if(allocated(p_r))   deallocate(p_r)
if(allocated(p_ur))  deallocate(p_ur)
if(allocated(p_vr))  deallocate(p_vr)
if(allocated(p_lat)) deallocate(p_lat)
if(allocated(p_lon)) deallocate(p_lon)
if(allocated(p_soilmoist1)) deallocate(p_soilmoist1)
if(allocated(p_soilmoist2)) deallocate(p_soilmoist2)
if(allocated(p_soiltemp1)) deallocate(p_soiltemp1)
if(allocated(p_soiltemp2)) deallocate(p_soiltemp2)
if(allocated(p_snowmass))  deallocate(p_snowmass)
if(allocated(p_snowdepth)) deallocate(p_snowdepth)

! Allocate memory for the pressure data

allocate(p_u(nprx,npry,nprz))
allocate(p_v(nprx,npry,nprz))
allocate(p_t(nprx,npry,nprz))
allocate(p_z(nprx,npry,nprz))
allocate(p_r(nprx,npry,nprz))

! p_ur,p_vr arrays for rotated winds used later

allocate(p_ur(nprx,npry,nprz))
allocate(p_vr(nprx,npry,nprz))

allocate(p_lat(nprx,npry))
allocate(p_lon(nprx,npry))

allocate(p_soilmoist1(nprx,npry))
allocate(p_soilmoist2(nprx,npry))
allocate(p_soiltemp1(nprx,npry))
allocate(p_soiltemp2(nprx,npry))
allocate(p_snowmass(nprx,npry))
allocate(p_snowdepth(nprx,npry))


! Fill with missing in case they are not present
p_soilmoist1 (1:nprx,1:npry)=1.e30
p_soilmoist2 (1:nprx,1:npry)=1.e30
p_soiltemp1 (1:nprx,1:npry)=1.e30
p_soiltemp2 (1:nprx,1:npry)=1.e30
p_snowmass  (1:nprx,1:npry)=1.e30
p_snowdepth (1:nprx,1:npry)=1.e30

print 92
92 format(//,1x,70('*')/  &
       ,'  Access coarse resolution pressure data',/,1X,70('*')//)

do k=1,nprz
   pnpr(k)=levpr(k)*100.
   print*,'Pressure levels (mb): ',pnpr(k)/100.
   if(k==nprz)print*,''
enddo
!Really need lowest level input to be 1000mb. Otherwise this can create
!some extrapolation issues to the lowest RAMS sigma levels causing temperature
!and moisture to be over- or under-estimated compared to reanalysis.
if(levpr(1) .ne. 1000) then
 print*,'Lowest reanalysis input levels needs to mb 1000mb.'
 print*,'Your lowest level input is ',levpr(1)
 stop 'pressure level missing'
endif

! Call routine to fill pressure arrays from the chosen dataset.

CALL get_press (iunit)


!Checking RH-MAX here to see if values are unreasonable.
!Also can see if we are working with RH or specific humidity (SPFH)
!Note that some data such as ERA5 (circa 2019) does not use 
!RH w.r.t. liquid for all levels, but RAMS deals with this w.r.t liquid.
!This is a problem, and the best solution is to extract and use SPFH instead.
rhmax=maxval(p_r, mask=p_r.lt.1.e29)
print*,"RHMAX or SpecificHumidity before= ",rhmax

!Saleeby(2014): Do this for NARR which gives specific humidity (SPFH) in kg/kg.
!Convert to RH from SPFH. If you know you are using SPFH instead of RH, then
!set the IDATAIN flag in RAMSIN to reflect this difference. RH is default.
if (IDATAIN == 1)then
   print*,'------------------------------------------------'
   print*,' Converting specific humidity to rh',nprx,npry,nprz,rhmax
   print*,'------------------------------------------------'
   do k=1,nprz
      do j=1,npry
         do i=1,nprx
            !Convert to mixing ratio assuming input specific humidity in kg/kg
            rr=p_r(i,j,k)/(1.-p_r(i,j,k))
            !Convert from mixing ratio (kg/kg) to RH (fraction, not %)
            p_r(i,j,k)=rr/rsatmix(pnpr(k),p_t(i,j,k))
            !Limit RH to > 10%)
            if(p_r(i,j,k) <= 0.1) p_r(i,j,k)=.1
         enddo
      enddo
   enddo
endif

!Re-check RH-MAX and see if conversion helped.
rhmax=maxval(p_r, mask=p_r.lt.1.e29)
print*,"RHMAX-after= ",rhmax

! Find max-min theta at bottom and top levels
thmax=1.
thmin=1000.
do j=1,npry
   do i=1,nprx
      if(p_t(i,j,nprz).lt.1.e20)  &
           thmax=max(thmax,p_t(i,j,nprz)*(p00/pnpr(nprz))**rocp)
      if(p_t(i,j,1).lt.1.e20)  &
           thmin=min(thmin,p_t(i,j,1)*(p00/pnpr(1))**rocp)
   enddo
enddo

print 300,levpr(1),thmin,levpr(nprz),thmax
300 format(//,' Minimum THETA at ',I4,' mb - ',F8.2/  &
          ,' Maximum THETA at ',I4,' mb - ',F8.2)


!  Compute lat-lon at input pressure data points

if (inproj==1) then
   do j=1,npry
      do i=1,nprx
         p_lat(i,j)=xswlat+(j-1)*gdatdy
         p_lon(i,j)=xswlon+(i-1)*gdatdx
      enddo
   enddo
elseif (inproj==2) then
   print*,'lamb-con:',xswlat,xswlon,cntlat,cntlon,gdatdx,gdatdy
   do j=1,npry
      do i=1,nprx
         CALL lc2ll (cntlat,cntlon,p_lat(i,j),p_lon(i,j)  &
                   ,float(i),float(j),xswlat,xswlon,gdatdx)
!         print*,i,j,p_lat(i,j),p_lon(i,j)
      enddo
   enddo
elseif (inproj==3) then
   print*,'polar:',cntlat,cntlon,gdatdx,gdatdy,xswlat,xswlon
   xswlon_east=xswlon
   if(xswlon<0.) xswlon_east=xswlon+360.
   cntlon_east=cntlon
   if(cntlon<0.) cntlon_east=cntlon+360.
   do j=1,npry
      do i=1,nprx
         CALL w3fb07 (float(i),float(j),xswlat,xswlon_east,gdatdx,cntlon_east  &
                    ,p_lat(i,j),p_lon(i,j))
!         print*,i,j,p_lat(i,j),p_lon(i,j)
      enddo
   enddo
endif

close(iunit)

return
END SUBROUTINE pressure_stage

!##############################################################################
Subroutine get_press (iunit)

use isan_coms

implicit none

integer :: iunit
real,allocatable::as(:,:)
integer :: ithere(maxpr,5),isfthere(6)
character(len=4) :: idat(5)
integer :: i,j,k,nv,nvar,misstot,lv,n

data idat/'T','R','U','V','H'/

! Allocate array to hold one level data for one variable

allocate(as(nprx,npry))

!  Read upper air fields

do lv=1,nprz
   do nvar=1,5

      read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)

      if(nvar.eq.1) then
         CALL prfill (nprx,npry,as,p_u(1,1,lv))
      elseif(nvar.eq.2) then
         CALL prfill (nprx,npry,as,p_v(1,1,lv))
      elseif(nvar.eq.3) then
         CALL prfill (nprx,npry,as,p_t(1,1,lv))
      elseif(nvar.eq.4) then
         CALL prfill (nprx,npry,as,p_z(1,1,lv))
      elseif(nvar.eq.5) then
         CALL prfill (nprx,npry,as,p_r(1,1,lv))
      endif

      print 555,levpr(lv),nvar,imonth,idate,iyear,ihour
      555 format(' ==  Read pressure field  ',2I4,' at ',I2,'/',I2  &
                ,'/',I4,I6.4,' UTC')
   enddo
enddo

!  Read surface fields

do nvar=1,6
   read(iunit,*,end=70,err=70) ((as(i,j),i=1,nprx),j=1,npry)

   if(nvar.eq.1) then
      CALL prfill (nprx,npry,as,p_soilmoist1(1,1))
   elseif(nvar.eq.2) then
      CALL prfill (nprx,npry,as,p_soilmoist2(1,1))
   elseif(nvar.eq.3) then
      CALL prfill (nprx,npry,as,p_soiltemp1(1,1))
   elseif(nvar.eq.4) then
      CALL prfill (nprx,npry,as,p_soiltemp2(1,1))
   elseif(nvar.eq.5) then
      CALL prfill (nprx,npry,as,p_snowmass(1,1))
   elseif(nvar.eq.6) then
      CALL prfill (nprx,npry,as,p_snowdepth(1,1))
   endif
enddo

goto 71

70 CONTINUE
print*,'Premature end of file or error in pressure input file!'
print*,'  (or optional surface fields are not present)'
print*,'We''ll close our eyes and pretend it didn''t happen!'
71 continue

deallocate(as)

! Check for levels that may be all missing

ithere(1:nprz,1:5)=0
isfthere(1:6)=0

do k=1,nprz
   do j=1,npry
      do i=1,nprx
         if(p_t(i,j,k) > 1e20) ithere(k,1)=ithere(k,1)+1
         if(p_r(i,j,k) > 1e20) ithere(k,2)=ithere(k,2)+1
         if(p_u(i,j,k) > 1e20) ithere(k,3)=ithere(k,3)+1
         if(p_v(i,j,k) > 1e20) ithere(k,4)=ithere(k,4)+1
         if(p_z(i,j,k) > 1e20) ithere(k,5)=ithere(k,5)+1
      enddo
   enddo
enddo
do nv=1,5
   do k=1,nprz
      if(ithere(k,nv) < nprx*npry) then
         ithere(k,nv)=1
      else
         ithere(k,nv)=0
      endif
   enddo
enddo

misstot=0
do nv=1,5
   do k=1,nprz
      if(ithere(k,nv) == 0) misstot=misstot+1
   enddo
enddo

do j=1,npry
   do i=1,nprx
      if(p_soilmoist1 (i,j) > 1e20) isfthere(1)=isfthere(1)+1
      if(p_soilmoist2 (i,j) > 1e20) isfthere(2)=isfthere(2)+1
      if(p_soiltemp1 (i,j) > 1e20) isfthere(3)=isfthere(3)+1
      if(p_soiltemp2 (i,j) > 1e20) isfthere(4)=isfthere(4)+1
      if(p_snowmass  (i,j) > 1e20) isfthere(5)=isfthere(5)+1
      if(p_snowdepth (i,j) > 1e20) isfthere(6)=isfthere(6)+1
   enddo
enddo
do nv=1,6
   if(isfthere(nv) < nprx*npry) then
      isfthere(nv)=1
   else
      isfthere(nv)=0
   endif
enddo




print*,'------------------------------------------------'
print*,' Missing parameters: 0 = all missing'
print*,'------------------------------------------------'
print '(t20,5(a1,6x))',(idat(n),n=1,5)

do k=1,nprz
   print '(f10.1,t14,5(i7))',pnpr(k),(ithere(k,n),n=1,5)
enddo

print*,'--------------------------------------------------------------'
print*,' Missing surface parameters: 0 = all missing'
print*,'--------------------------------------------------------------'
print*,'    SOILM1   SOILM2   SOILT1   SOILT2  SNOWMASS   SNOWDEPTH'

print '(6(i9))',(isfthere(n),n=1,6)


if(misstot.gt.0) then
   ! Let's see if we can get creative and make up data for the missing fields

   CALL press_miss (nprx,npry,nprz,p_u,p_v,p_t,p_z,p_r  &
                  ,ithere,maxpr,levpr)

  ! Check again for missing fields

   misstot=0
   do nv=1,5
      do k=1,nprz
         if(ithere(k,nv).eq.0) misstot=misstot+1
      enddo
   enddo

   print*,'------------------------------------------------'
   print*,' After missing parameters check: 0 = all missing'
   print*,'------------------------------------------------'
   print '(t20,5(a1,6x))',(idat(n),n=1,5)

   do k=1,nprz
      print '(f10.1,t14,5(i7))',pnpr(k),(ithere(k,n),n=1,5)
   enddo

endif

return
END SUBROUTINE get_press

!##############################################################################
Subroutine prfill (nprx,npry,xx,dn)

implicit none

integer :: nprx,npry
real :: xx(nprx,npry),dn(nprx,npry)
integer :: i,j

do j=1,npry
   do i=1,nprx
      dn(i,j)=xx(i,j)
      !Saleeby(2019): May need to change to -9998.0 if value between dprep
      ! of -999 for missing or bad data has to be lowered.
      if(dn(i,j).le.-998.) then
        dn(i,j)=1.e30
      endif
   enddo
enddo

return
END SUBROUTINE prfill

!##############################################################################
Subroutine press_miss (n1,n2,n3,un,vn,tn,zn,rn,ithere,maxpr,levpr)

implicit none

integer :: n1,n2,n3,maxpr
real, dimension(n1,n2,n3) ::  un,vn,tn,zn,rn
integer :: ithere(maxpr,5),levpr(*)
real :: prs(n3),prsln(n3)
integer :: it=1,ir=2,iu=3,iv=4,iz=5
integer :: ierr,k

do k=1,n3
   prs(k)=float(levpr(k))
   prsln(k)=log(prs(k))
enddo

! first do moisture. since there are no physical relationships we can
!   use to help us, we will simply interpolate to a missing level. If
!   the top or bottom level is missing, we will fill it with the relative
!   humidity above or below.

CALL pr_miss_fill (n1,n2,n3,rn,ithere(1,ir),ithere(1,ir)  &
     ,prsln,'rel hum',ierr)
if(ierr.eq.1) return

! do temperature in a similar manner, but only if there are levels
!   where height is missing also

CALL pr_miss_fill (n1,n2,n3,tn,ithere(1,it),ithere(1,iz)  &
     ,prsln,'temp',ierr)
if(ierr.eq.1) return


! check height at each level. It can be computed hydrostatically
!   from temp as long as temp is not missing at this level and temp and
!   height is available on a level above or below.
!   Make a downward sweep first, then an upward.


do k=n3-1,1,-1
   if(ithere(k,iz).eq.0) then
      ! check for temp, z, above
      if(ithere(k+1,iz).eq.1.and.ithere(k+1,it).eq.1.and.  &
         ithere(k,it).eq.1) then
         !ok, we can get this
         CALL pr_hystatic_z (n1*n2,zn(1,1,k),zn(1,1,k+1)  &
                            ,tn(1,1,k),tn(1,1,k+1),rn(1,1,k),rn(1,1,k+1)  &
                            ,prs(k),prs(k+1))
         ithere(k,iz)=1
         print*,'-->Computing hydrostatic pressure at level:',k
      endif
   endif
enddo

do k=2,n3
   if(ithere(k,iz).eq.0) then
      ! check for temp, z, below
      if(ithere(k-1,iz).eq.1.and.ithere(k-1,it).eq.1.and.  &
           ithere(k,it).eq.1) then
         ! ok, we can get this
         CALL pr_hystatic_z (n1*n2,zn(1,1,k),zn(1,1,k-1)  &
                            ,tn(1,1,k),tn(1,1,k-1),rn(1,1,k),rn(1,1,k-1)  &
                            ,prs(k),prs(k-1))
         print*,'-->Computing hydrostatic pressure at level:',k
         ithere(k,iz)=1
      endif
   endif
enddo

! try temperature in a similar manner. It can also be computed
!   from height as long as height is not missing at this level and temp and
!   height is available on a level above or below.
!   Note that vapor mixing ratio (used to compute virtual temperature)
!   is somewhat difficult to compute from relative humidity if temperature
!   is missing. Therefore, the vitual temperature factor at the
!   non-missing level is assumed at the missing level.

do k=n3-1,1,-1
   if(ithere(k,it).eq.0) then
      ! check for temp, z, above
      if(ithere(k+1,iz).eq.1.and.ithere(k+1,it).eq.1.and.  &
         ithere(k,iz).eq.1) then
         ! ok, we can get this
         CALL pr_hystatic_t (n1*n2,zn(1,1,k),zn(1,1,k+1)  &
                           ,tn(1,1,k),tn(1,1,k+1),rn(1,1,k),rn(1,1,k+1)  &
                           ,prs(k),prs(k+1))
         ithere(k,it)=1
         print*,'-->Computing hydrostatic temperature at level:',k
      endif
   endif
enddo

do k=2,n3
   if(ithere(k,it).eq.0) then
      ! check for temp, z, below
      if(ithere(k-1,iz).eq.1.and.ithere(k-1,it).eq.1.and.  &
         ithere(k,iz).eq.1) then
         ! ok, we can get this
         CALL pr_hystatic_t (n1*n2,zn(1,1,k),zn(1,1,k-1)  &
                            ,tn(1,1,k),tn(1,1,k-1),rn(1,1,k),rn(1,1,k-1)  &
                            ,prs(k),prs(k-1))
         ithere(k,it)=1
         print*,'-->Computing hydrostatic temperature at level:',k
      endif
   endif
enddo

! For the u and v components, do a straight interpolation again like we
!   did for rel humidity.

CALL pr_miss_fill (n1,n2,n3,un,ithere(1,iu),ithere(1,iu),prsln,'u-comp',ierr)
if(ierr.eq.1) return

CALL pr_miss_fill (n1,n2,n3,vn,ithere(1,iv),ithere(1,iv),prsln,'v-comp',ierr)
if(ierr.eq.1) return

return
END SUBROUTINE press_miss

!##############################################################################
Subroutine pr_hystatic_z (np,z1,z2,t1,t2,r1,r2,p1,p2)

implicit none

integer :: np
real :: z1(np),z2(np),t1(np),t2(np),r1(np),r2(np)
integer :: n
real :: p1,p2,tv1,tv2
real, external :: rslf

do n=1,np
   if(z2(n).lt.1.e20.and.t1(n).lt.1.e20.and.  &
      t2(n).lt.1.e20.and.r1(n).lt.1.e20.and.  &
      r2(n).lt.1.e20 ) then
      tv1=t1(n)*(1.+.61*rslf(p1*100.,t1(n))*r1(n))
      tv2=t2(n)*(1.+.61*rslf(p2*100.,t2(n))*r2(n))
      z1(n)=z2(n)- 287.*.5*(tv1+tv2)*(log(p1*100.)-log(p2*100.))/9.8
      !print*,z1(n),z2(n),t1(n),t2(n),tv1,tv2,p1,p2
   else
      z1(n)=1.e30
   endif
enddo

return
END SUBROUTINE pr_hystatic_z

!##############################################################################
Subroutine pr_hystatic_t (np,z1,z2,t1,t2,r1,r2,p1,p2)

implicit none

integer :: np
real :: z1(np),z2(np),t1(np),t2(np),r1(np),r2(np)
integer :: n
real :: p1,p2,vtfact
real, external :: rslf

do n=1,np
   if(t2(n).lt.1.e20.and.z1(n).lt.1.e20  &
        .and.z2(n).lt.1.e20.and.r1(n).lt.1.e20  &
        .and.r2(n).lt.1.e20 ) then
      vtfact=(1.+.61*rslf(p2*100.,t2(n))*r2(n))
      t1(n)=-t2(n)- (2.*9.8*(z1(n)-z2(n))  &
           /(287. *(log(p1*100.)-log(p2*100.))) )  &
           /vtfact
      !if(n.eq.36*(13-1)+13) print*,t1(n),t2(n),z1(n),z2(n),vtfact,p1,p2 
   else
      t1(n)=1.e30
   endif
enddo

return
END SUBROUTINE pr_hystatic_t

!##############################################################################
Subroutine pr_miss_fill (n1,n2,n3,a,ithere,ithere2,prsln,varname,ierr)

! simple filling for missing pressure fields. Assuming there are no
!   physical relationships we can
!   use to help us, we will simply interpolate to a missing level. If
!   the top or bottom level is missing, we will fill it with the next
!   non-missing field above or below.

implicit none

integer :: n1,n2,n3,ithere(n3),ithere2(n3),ierr
real :: a(n1,n2,n3),prsln(n3)
character(len=*) :: varname
integer :: k,kk,ilev

ierr=0
k=1
if(ithere(k).eq.0.and.ithere2(k).eq.0) then
   ! find first non-missing level above
   ilev=0
   do kk=2,n3
      if (ithere(kk).eq.1) then
         ilev=kk
         goto 10
      endif
   enddo
   print*,'All data missing; fixing stopped for:',varname
   ierr=1
   return
   10 continue
   CALL atob (n1*n2,a(1,1,ilev),a(1,1,k))
   ithere(k)=1
   print*,'-->Filling:',varname,' at level:',k,' from level:',ilev
endif

k=n3
if(ithere(k).eq.0.and.ithere2(k).eq.0) then
   ! find first non-missing level below
   ilev=0
   do kk=n3-1,1,-1
      if (ithere(kk).eq.1) then
         ilev=kk
         goto 11
      endif
   enddo
   print*,'All data missing; fixing stopped for:',varname
   ierr=1
   return
   11 continue
   CALL atob (n1*n2,a(1,1,ilev),a(1,1,k))
   ithere(k)=1
   print*,'-->Filling:',varname,' at level:',k,' from level:',ilev
endif

! Now interpolate to internal missing levels

do k=2,n3-1
   if (ithere(k).eq.0.and.ithere2(k).eq.0) then
      ! find first non-missing level above
      ilev=0
      do kk=k+1,n3
         if (ithere(kk).eq.1) then
            ilev=kk
            goto 12
         endif
      enddo
      stop 'pr_miss_fill'
      12 continue
      CALL pr_interp (n1*n2,a(1,1,k-1),a(1,1,k),a(1,1,ilev)  &
                     ,prsln(k-1),prsln(k),prsln(ilev))
      ithere(k)=1
      print*,'-->Interpolating:',varname,' at level:',k  &
            ,' from levels:',k-1,' and:',ilev

   endif
enddo

return
END SUBROUTINE pr_miss_fill

!##############################################################################
Subroutine pr_interp (np,a1,a2,a3,pln1,pln2,pln3)

implicit none

integer :: np
real :: a1(np),a2(np),a3(np),pln1,pln2,pln3
integer :: n

do n=1,np
   if(a3(n).lt.1.e20.and.a1(n).lt.1.e20) then
      a2(n)=a1(n)+  (pln2-pln1)* (a3(n)-a1(n))/(pln3-pln1)
   else
      a2(n)=1.e30
   endif
enddo

return
END SUBROUTINE pr_interp

