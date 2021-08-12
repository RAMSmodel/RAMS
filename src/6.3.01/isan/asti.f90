!##############################################################################
Subroutine isnstage ()

use isan_coms
use mem_grid

implicit none

character(len=8) :: rot_type
real :: reflat1,reflon1,olat1,olat2,olon1,olon2  &
       ,olatmin,olatmax,olonmin,olonmax
integer :: ng

! ---------------------------------------------start  ngrid =1 section
if(ngrid.eq.1) then

!-----------------------------------------------------------------
!      Read rawinsondes and surface observations
!-----------------------------------------------------------------
!      Default to RAMS grid 1 for obs access and plotting
!-----------------------------------------------------------------

   CALL rams_mm (grid_g(ngrid)%glat(1,1),nnxp(1)*nnyp(1),olatmin,olatmax)
   CALL rams_mm (grid_g(ngrid)%glon(1,1),nnxp(1)*nnyp(1),olonmin,olonmax)

   olat1=olatmin-5.
   olat2=olatmax+5.
   olon1=olonmin-5.
   olon2=olonmax+5.

   nsta=0

   if(iproc_flag(natime,3).eq.1) then

      if(.not.allocated(up_p))   allocate(up_p(maxsta,maxlev))
      if(.not.allocated(up_t))   allocate(up_t(maxsta,maxlev))
      if(.not.allocated(up_z))   allocate(up_z(maxsta,maxlev))
      if(.not.allocated(up_r))   allocate(up_r(maxsta,maxlev))
      if(.not.allocated(up_uz))  allocate(up_uz(maxsta,maxlev))
      if(.not.allocated(up_vz))  allocate(up_vz(maxsta,maxlev))
      if(.not.allocated(up_ur))  allocate(up_ur(maxsta,maxlev))
      if(.not.allocated(up_vr))  allocate(up_vr(maxsta,maxlev))
      if(.not.allocated(up_zz))  allocate(up_zz(maxsta,maxlev))

      if(.not.allocated(up_lat)) allocate(up_lat(maxsta))
      if(.not.allocated(up_lon)) allocate(up_lon(maxsta))
      if(.not.allocated(up_top)) allocate(up_top(maxsta))
      if(.not.allocated(up_lp))  allocate(up_lp(maxsta))
      if(.not.allocated(up_lz))  allocate(up_lz(maxsta))
      if(.not.allocated(up_topg))allocate(up_topg(maxsta,nigrids))

      if(.not.allocated(up_chstid)) allocate (up_chstid(maxsta))

      inrawi=iproc_names(natime,2)

      CALL input_rawi ()
           
      do ng=1,nigrids

         if(ng.gt.1) up_topg(1:maxsta,ng)=up_topg(1:maxsta,nxtnest(ng)) 

         CALL soundtopo (ng,nsta,maxsta,up_topg(1,ng)  &
                        ,up_lat,up_lon,up_top  &
                        ,grid_g(ng)%topt(1,1) &
                        ,nnxp(ng),nnyp(ng),polelat,polelon  &
                        ,xtn(1,ng),ytn(1,ng),deltaxn(ng))
      enddo


   endif

!-----------------------------------------------------------------
!     Read surface observations
!-----------------------------------------------------------------

   nssfc=0

   if(iproc_flag(natime,4).eq.1) then

      
      if(.not.allocated(sf_u)) allocate(sf_u(maxsfc))
      if(.not.allocated(sf_v)) allocate(sf_v(maxsfc))
      if(.not.allocated(sf_ur)) allocate(sf_ur(maxsfc))
      if(.not.allocated(sf_vr)) allocate(sf_vr(maxsfc))
      if(.not.allocated(sf_p)) allocate(sf_p(maxsfc))
      if(.not.allocated(sf_t)) allocate(sf_t(maxsfc))
      if(.not.allocated(sf_s)) allocate(sf_s(maxsfc))
      if(.not.allocated(sf_r)) allocate(sf_r(maxsfc))
      if(.not.allocated(sf_lat)) allocate(sf_lat(maxsfc))
      if(.not.allocated(sf_lon)) allocate(sf_lon(maxsfc))
      if(.not.allocated(sf_top)) allocate(sf_top(maxsfc))
      if(.not.allocated(sf_scra)) allocate(sf_scra(maxsfc))
      
      if(.not.allocated(sf_chstid)) allocate (sf_chstid(maxsfc))
      if(.not.allocated(sf_date)) allocate (sf_date(maxsfc))

      insrfce=iproc_names(natime,3)

      CALL input_sfc ()
  
   endif

endif


! end of ngrid =1 section

!-----------------------------------------------------------------
!    Interpolate pressure data to RAMS grid
!-----------------------------------------------------------------

if(igridfl.gt.0) then

!         First, interpolate first guess pressure data horizontally to the 
!                RAMS polar-stereo grid
!         --------------------------------------------------------

      print*,'Allocating polar/pressure grid-'  &
           ,nnxp(ngrid),nnyp(ngrid),nprz,ngrid

      if(.not.allocated(pp_u)) allocate(pp_u(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_v)) allocate(pp_v(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_t)) allocate(pp_t(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_z)) allocate(pp_z(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_r)) allocate(pp_r(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_sglob)) allocate(pp_sglob(nprx+3,npry+2))
      
!         Rotate winds to the RAMS polar-stereo grid
!         --------------------------------------------------------
      if(inproj.eq.1) then
         rot_type='ll_rps'
         reflat1=0.
         reflon1=0.
      elseif(inproj.eq.2) then
         rot_type='lc_rps'
         reflat1=cntlat
         reflon1=cntlon
      elseif(inproj.eq.3) then
         rot_type='ps_rps'
         reflat1=cntlat
         reflon1=cntlon
      endif
         
      CALL rotate_winds (rot_type,nprx,npry,nprz  &
                       ,p_u(1,1,1),p_v(1,1,1),p_ur(1,1,1),p_vr(1,1,1)  &
                       ,p_lat(1,1),p_lon(1,1)  &
                       ,reflat1,reflon1,polelat,polelon,idatain)

      if(inproj.eq.1) then
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,pp_sglob,p_ur,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,pp_sglob,p_vr,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,pp_sglob,p_t,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,pp_sglob,p_z,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,pp_sglob,p_r,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
              ! Surface fields 
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist1,pp_sglob,p_soilmoist1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist2,pp_sglob,p_soilmoist2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp1,pp_sglob,p_soiltemp1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp2,pp_sglob,p_soiltemp2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowmass,pp_sglob,p_snowmass,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         CALL latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowdepth,pp_sglob,p_snowdepth,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)

      elseif(inproj.eq.2) then
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,p_ur,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,p_vr,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,p_t,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,p_z,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,p_r,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
              ! Surface fields 
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist1,p_soilmoist1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist2,p_soilmoist2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp1,p_soiltemp1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp2,p_soiltemp2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowmass,p_snowmass,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)
         CALL lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowdepth,p_snowdepth,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlat,cntlon)

      elseif(inproj.eq.3) then
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,p_ur,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,p_vr,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,p_t,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,p_z,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,p_r,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
              ! Surface fields 
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist1,p_soilmoist1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soilmoist2,p_soilmoist2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp1,p_soiltemp1,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_soiltemp2,p_soiltemp2,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowmass,p_snowmass,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
         CALL trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snowdepth,p_snowdepth,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
             ,xswlat,xswlon,gdatdx,cntlon)
      endif
      
!         Next, interpolate vertically to the specified isentropes
!         --------------------------------------------------------
      CALL vterpp_i (nnxp(ngrid),nnyp(ngrid),nprz,nisn  &
                    ,pp_u,pp_v,pp_t,pp_z,pp_r  &
                    ,pi_u,pi_v,pi_p,pi_s,pi_r)

!         Last, interpolate vertically to the sigma-z levels
!         --------------------------------------------------------

      CALL vterpp_s (nnxp(ngrid),nnyp(ngrid),nprz,nsigz  &
                    ,pp_u,pp_v,pp_t,pp_z,pp_r  &
                    ,ps_u,ps_v,ps_p,ps_t,ps_r  &
                    ,grid_g(ngrid)%topt(1,1)  &
                    ,grid_g(ngrid)%rtgt(1,1))

      if(allocated(pp_u)) deallocate(pp_u)
      if(allocated(pp_v)) deallocate(pp_v)
      if(allocated(pp_t)) deallocate(pp_t)
      if(allocated(pp_z)) deallocate(pp_z)
      if(allocated(pp_r)) deallocate(pp_r)
      if(allocated(pp_sglob)) deallocate(pp_sglob)

endif

!-----------------------------------------------------------------
!     Do objective analysis or interpolation to RAMS grids
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!     If only using gridded p data, fill rs_qual, then we're done...
!-----------------------------------------------------------------

if(igridfl.eq.3) then

   CALL strmfun (nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%topt(1,1)  &
                    ,grid_g(ngrid)%rtgt(1,1))

   rs_qual(1:nnxp(ngrid),1:nnyp(ngrid)) = 0.

   goto 1200
endif

!-----------------------------------------------------------------
!           Perform upper air objective analysis/interpolation
!-----------------------------------------------------------------

!    Rotate obs winds to RAMS polar-stereo grid  
!        (reflat1,reflon1) doesn't matter since obs winds 
!         are earth relative
   
if(iproc_flag(natime,3).eq.1) then
   CALL rotate_winds ('ll_rps',maxsta,1,maxlev  &
                    ,up_uz,up_vz,up_ur,up_vr  &
                    ,up_lat,up_lon  &
                    ,reflat1,reflon1,polelat,polelon,idatain)
endif

if(iproc_flag(natime,4).eq.1) then
   CALL rotate_winds ('ll_rps',nssfc,1,1 &
                    ,sf_u,sf_v,sf_ur,sf_vr  &
                    ,sf_lat,sf_lon  &
                    ,reflat1,reflon1,polelat,polelon,idatain)
endif

!                 Isentropic coordinates
!-----------------------------------------------------------------

!         Interpolate rawindsondes vertically to isentropes

   print*,'Allocating obs/isen array-',nsta,nisn
   
   allocate(upi_u(nsta,nisn))
   allocate(upi_v(nsta,nisn))
   allocate(upi_p(nsta,nisn))
   allocate(upi_s(nsta,nisn))
   allocate(upi_r(nsta,nisn))

   if (nsta > 0) CALL obs_isen (maxsta,maxlev  &
                 ,up_t,up_p,up_z,up_r,up_ur,up_vr,up_zz,up_lp,up_lz  &
                 ,upi_u,upi_v,upi_p,upi_s,upi_r,nsta,up_chstid)

   CALL obj_anal ('isen',ngrid,nnxp(ngrid),nnyp(ngrid)  &
                 ,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
                 ,polelat,polelon,xtn(1,ngrid),ytn(1,ngrid)  &
                 ,deltaxn(ngrid))
              
   deallocate (upi_u,upi_v,upi_p,upi_s,upi_r)

!                 Sigma-z coordinates
!-----------------------------------------------------------------

!         Interpolate rawindsondes vertically to sigma-z

print*,'Allocating obs/sigma-z array-',nsta,nsigz

allocate(ups_u(nsta,nsigz))
allocate(ups_v(nsta,nsigz))
allocate(ups_p(nsta,nsigz))
allocate(ups_t(nsta,nsigz))
allocate(ups_r(nsta,nsigz))
   
if (nsta > 0)  CALL obs_sigz (maxsta,maxlev,up_t,up_p,up_z &
              ,up_r,up_ur,up_vr,up_zz,up_lp,up_lz,up_topg(1,ngrid) &
              ,ups_u,ups_v,ups_p,ups_t,ups_r,nsta,zmn(nnzp(1)-1,1))

CALL obj_anal ('sigz',ngrid,nnxp(ngrid),nnyp(ngrid)  &
              ,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
              ,polelat,polelon,xtn(1,ngrid),ytn(1,ngrid)  &
              ,deltaxn(ngrid))

deallocate (ups_u,ups_v,ups_p,ups_t,ups_r)

!-----------------------------------------------------------------
!          Find Montgomery streamfunc on isentropic surfaces
!-----------------------------------------------------------------

CALL strmfun (nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%topt(1,1)  &
            ,grid_g(ngrid)%rtgt(1,1))

!-----------------------------------------------------------------
!          Perform ground surface objective analysis.
!-----------------------------------------------------------------

CALL obj_anal ('surf',ngrid,nnxp(ngrid),nnyp(ngrid) &
           ,grid_g(ngrid)%glat(1,1),grid_g(ngrid)%glon(1,1)  &
           ,polelat,polelon,xtn(1,ngrid),ytn(1,ngrid)  &
           ,deltaxn(ngrid))


! Check if "quality" of surface grid point is good

CALL sfcqual (nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%glat(1,1)  &
             ,grid_g(ngrid)%glon(1,1)  &
             ,rs_qual,sf_lat,sf_lon,sf_p,nssfc,sf_scra)

1200 continue

return
END SUBROUTINE isnstage

!##############################################################################
Subroutine sfcqual (nxp,nyp,glat,glon,quals3  &
                   ,sflt,sfln,pss,nssfc,scra)

use rconstants

implicit none

integer :: nxp,nyp,nssfc
real ::  glat(nxp,nyp),glon(nxp,nyp),sflt(*),sfln(*),scra(*)  &
        ,quals3(nxp,nyp),pss(*)

integer :: iqf(4),i,j,ns
real :: gobkm,gobrd

gobkm=2.*111.120
gobrd=4.*111.120

do i=1,nxp
   do j=1,nyp
      quals3(i,j)=0.
      iqf(1)=0
      iqf(2)=0
      iqf(3)=0
      iqf(4)=0
      do ns=1,nssfc
         scra(ns)=sqrt(((glat(i,j)-sflt(ns))*111.12)**2  &
                 +((glon(i,j)-sfln(ns))*111.12  &
                 *cos((glat(i,j)+sflt(ns))*.5*.01745))**2)
      enddo
      do ns=1,nssfc
         if(pss(ns).lt.1e19.and.scra(ns).lt.gobkm) then
            quals3(i,j)=1.
         endif
      enddo
      do ns=1,nssfc
         if(scra(ns).le.gobrd) then
            if(sfln(ns).le.glon(i,j).and.sflt(ns).ge.glat(i,j)) iqf(1)=1
            if(sfln(ns).ge.glon(i,j).and.sflt(ns).ge.glat(i,j)) iqf(2)=1
            if(sfln(ns).ge.glon(i,j).and.sflt(ns).le.glat(i,j)) iqf(3)=1
            if(sfln(ns).le.glon(i,j).and.sflt(ns).le.glat(i,j)) iqf(4)=1
         endif
      enddo
      if(iqf(1)+iqf(2)+iqf(3)+iqf(4).ge.3) quals3(i,j)=quals3(i,j)+2.
   enddo
enddo

return
END SUBROUTINE sfcqual

!##############################################################################
Subroutine strmfun (nxp,nyp,topt,rtgt)

use isan_coms
use rconstants

implicit none

integer :: nxp,nyp
real, dimension(nxp,nyp) :: topt,rtgt

real, dimension(maxsigz) :: sigzr,temp,thv

integer :: k,i,j,lbchyd,lbc,lbcp
real, external :: rsatmix
real :: syo,po,tho,thvo,sigo,bcpr

   ! Find a boundary condition as first level at or below 360K
 
   lbchyd=360
   do k=nisn,1,-1
      if(levth(k).le.lbchyd) then
         do i=1,nxp
            do j=1,nyp
               if(pi_p(i,j,k).gt.1e19) goto 50
            enddo
         enddo
         lbc=k
         print 54,lbc,levth(lbc)
         54 format(///,' isentropic hydrostatic boundary set at level',2i6,' k')
         goto 52
         50 continue
      endif
   enddo
   print 53,i,j,k
   53 format(' could not find good isentropic boundary level ',3i6)
   stop 'st3-under'
   52 continue
 
   do i=1,nxp
      do j=1,nyp
 
        syo=pi_s(i,j,lbc)
        po=pi_p(i,j,lbc)
        tho=levth(lbc)
        !Integrate from isen hydrostatic B.C. to top
        do k=lbc+1,nisn
           pi_s(i,j,k)=1e30
           if(pi_p(i,j,k).lt.1e19) then
               pi_s(i,j,k)=syo+cp*(po**rocp+pi_p(i,j,k)**rocp)  &
                    *.5/p00**rocp *(levth(k)-tho)
               syo=pi_s(i,j,k)
               po=pi_p(i,j,k)
               tho=levth(k)
            endif
         enddo
 
         syo=pi_s(i,j,lbc)
         po=pi_p(i,j,lbc)
         tho=levth(lbc)
         !Integrate from isen hydrostatic B.C. to bottom
         do k=lbc-1,1,-1
            pi_s(i,j,k)=1e30
            if(pi_p(i,j,k).lt.1e19) then
               pi_s(i,j,k)=syo+cp*(po**rocp+pi_p(i,j,k)**rocp)  &
                    *.5/p00**rocp*(levth(k)-tho)
               syo=pi_s(i,j,k)
               po=pi_p(i,j,k)
               tho=levth(k)
            endif
         enddo
 
      enddo
   enddo

bcpr=10000.
do k=nsigz,1,-1
   if(sigz(k).le.bcpr) then
      do i=1,nxp
         do j=1,nyp
            if(ps_t(i,j,k).gt.1e19) goto 60
         enddo
      enddo
      lbcp=k
      print 64,lbcp,sigz(lbcp)
      64 format(///,' Sigma-z hydrostatic boundary set at level',i6,f10.1,' m')
      goto 62
      60 continue
   endif
enddo
print 63
63 format(' Could not find good sigma-z boundary level ')
stop 'st3-unders'
62 continue

do i=1,nxp
   do j=1,nyp

      !Compute sigma-z heights and theta-v
      do k=1,nsigz
         temp(k)=ps_t(i,j,k)*(ps_p(i,j,k)/p00)**rocp
         sigzr(k)=topt(i,j)+sigz(k)*rtgt(i,j)
         thv(k)=ps_t(i,j,k)*(1.+.61*ps_r(i,j,k)  &
              *rsatmix(ps_p(i,j,k),temp(k)))
      enddo

      !Compute PI from pressure at sigma-z hydrostatic boundary level
      ps_p(i,j,lbcp)=cp*(ps_p(i,j,lbcp)/p00)**rocp

      thvo=thv(lbcp)
      po=ps_p(i,j,lbcp)
      sigo=sigzr(lbcp)
      !Integrate PI from upper hydrostatic boundary to top sigma lev
      do k=lbcp+1,nsigz
         ps_p(i,j,k)=1e30
         if(ps_t(i,j,k).lt.1e19.and.ps_r(i,j,k).lt.1.e19) then
            ps_p(i,j,k)=po-g*(sigzr(k)-sigo)/((thvo+thv(k))*.5)
            thvo=thv(k)
            po=ps_p(i,j,k)
            sigo=sigzr(k)
         endif
      enddo

      thvo=thv(lbcp)
      po=ps_p(i,j,lbcp)
      sigo=sigzr(lbcp)
      !Integrate PI from upper hydrostatic boundary to bottom sigma lev
      do k=lbcp-1,1,-1
         ps_p(i,j,k)=1e30
         if(ps_t(i,j,k).lt.1e19.and.ps_r(i,j,k).lt.1.e19) then
            ps_p(i,j,k)=po+g*(sigo-sigzr(k))/((thvo+thv(k))*.5)
            thvo=thv(k)
            po=ps_p(i,j,k)
            sigo=sigzr(k)
         endif
      enddo
      !Convert PI back to pressure
      do k=1,nsigz
         ps_p(i,j,k)=(ps_p(i,j,k)/cp)**cpor*p00
      enddo

      !Check some output here if needed
      !if(i==74.and.j==135)then
      ! do k=nsigz,1,-1
      !  write(*,'(a,1(I4,1X),4(F7.2,1X),1(F9.2,1X),1(F6.2,1X))') 'strmfun1' &
      !    ,k,ps_r(i,j,k),temp(k),ps_t(i,j,k),thv(k) &
      !    ,ps_p(i,j,k)/100.,(1.e3)*rsatmix(ps_p(i,j,k),temp(k))
      ! enddo
      ! do k=nisn,1,-1
      !  write(*,'(a,2(I4,1X),2(F11.3,1X))') 'strmfun2' &
      !    ,k,levth(k),pi_s(i,j,k),pi_p(i,j,k)/100.
      ! enddo
      !endif

   enddo
enddo

return
END SUBROUTINE strmfun

!##############################################################################
Subroutine soundtopo (ngr,nst,m1,sndtopg,stlt,stln,topsnd,topt  &
                     ,n1,n2,polat,polon,swx,swy,delx)

implicit none                     

integer ::     ngr,nst,m1 ,n1,n2
real :: stlt(m1),stln(m1),sndtopg(m1),topsnd(m1)  &
     ,topt(n1,n2),sdat(2)
real :: polat,polon,swx,swy,delx

integer :: ns
real :: topo,stx,sty

do ns=1,nst
   sdat(1)=0.
   CALL ll_xy (stlt(ns),stln(ns),polat,polon,stx,sty)
   CALL stainterp (topt,n1,n2,1,1,1,sdat,stx,sty  &
                 ,1,swx,swy,delx,-1.e30)
   topo=-sdat(2)
   if(ngr.eq.1) then
      if(topo.lt.1e20) then
         sndtopg(ns)=topo
      else
         sndtopg(ns)=topsnd(ns)
      endif
   else
      if(topo.lt.1e20) sndtopg(ns)=topo
   endif
enddo

return
END SUBROUTINE soundtopo

