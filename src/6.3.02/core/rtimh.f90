!##############################################################################
Subroutine timestep ()

use mem_basic
use node_mod
use mem_radiate
use mem_cuparm
use mem_varinit
use mem_turb
use mem_oda,   only:if_oda
use micphys,   only:level,icheckmic,iscm,scmtime,iscmx,iscmy
use mem_grid
use kpp_parameters, only:IKPP

implicit none

integer :: callmassflux,massfluxfreq
real, dimension(mzp,mxp,myp) :: thvlast

 CALL acctimes ('INIT')

!        +-------------------------------------------------------------+
!        |   Timestep driver for the hybrid non-hydrostatic time-split |
!        |      model.                                                 |
!        +-------------------------------------------------------------+

!------------------------------------------------------------------------------
!  Zero out all tendency arrays.   
!------------------------------------------------------------------------------
 CALL tend0 ()          
 CALL acctimes ('TEND0')

!------------------------------------------------------------------------------
!  Thermodynamic diagnosis
!------------------------------------------------------------------------------
 CALL thermo () 
 CALL acctimes ('THERMO')

!------------------------------------------------------------------------------
!For computing FULL EXNER FUNCTION (after thermo update)
!This is experimental. Do not use yet. (Nov 19, 2017)
!------------------------------------------------------------------------------
! CALL exevolve (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,'ADV',thvlast)
! CALL acctimes ('EXEVOLVE1')

!------------------------------------------------------------------------------
!  Radiation parameterization
!  If running Harrington radiation, this radiation call only updates
!  theta-il tendency, zeros out fthrd, runs radprep and then mclatchy.
!  Actual Harrinton call done from microphysics driver for Level=3 micro.
!  For Chen-Cotton or Mahrer-Pielke radiation options, the radiation
!  physics is done at this point. It updates the tendency but does not
!  apply it until the following timestep. Updates fluxes needed for 
!  land-surface model which is called next. Radiation is typically
!  called on infrequent timesteps given by namelist variable "radfrq".
!------------------------------------------------------------------------------
 CALL radiate (mzp,mxp,myp,ia,iz,ja,jz) 
 CALL acctimes ('RADIATE')

!------------------------------------------------------------------------------
!  Surface layer, soil, veggie, urban models
!------------------------------------------------------------------------------
 CALL sfc_driver (mzp,mxp,myp,ia,iz,ja,jz) 
 CALL acctimes ('SFC_DRIVER')

!------------------------------------------------------------------------------
!  Coriolis terms
!------------------------------------------------------------------------------
 CALL corlos (mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv) 
 CALL acctimes ('CORLOS')

!------------------------------------------------------------------------------
!  Velocity advection
!------------------------------------------------------------------------------
 CALL advectc ('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv)
 CALL acctimes ('ADVECTv')

!------------------------------------------------------------------------------
!  Cumulus parameterization
!------------------------------------------------------------------------------
 IF(NNQPARM(ngrid) == 1 ) CALL cuparm ()      
 IF(NNQPARM(ngrid) == 2 ) CALL kf_main ()      
 CALL acctimes ('CUPARM')

!------------------------------------------------------------------------------
!  Analysis nudging and boundary condition
!------------------------------------------------------------------------------
 IF(NUD_TYPE > 0) CALL datassim ()  
 CALL acctimes ('DATASSIM')

!------------------------------------------------------------------------------
!  Observation data assimilation
!------------------------------------------------------------------------------
 IF(IF_ODA == 1) CALL oda_nudge ()  
 CALL acctimes ('DATASSIM')

!------------------------------------------------------------------------------
!  Nested grid boundaries
!------------------------------------------------------------------------------
 if(nxtnest(ngrid) >= 1) CALL nstbdriv ()  
 CALL acctimes ('NSTBDRIV')

!------------------------------------------------------------------------------
!  Rayleigh friction for theta
!------------------------------------------------------------------------------
 CALL rayft ()           
 CALL acctimes ('RAYFT')

!------------------------------------------------------------------------------
!  Update the overlap region between parallel nodes
!------------------------------------------------------------------------------
 if(nmachs .gt. 1) CALL update_lbc_vgroup (ngrid,LBC_ALL_VARS)  
 CALL acctimes ('UPDATElbc')

!------------------------------------------------------------------------------
!  Sub-grid diffusion terms
!------------------------------------------------------------------------------
 CALL diffuse ()
 CALL acctimes ('DIFFUSE')

!------------------------------------------------------------------------------
!  Scalar advection
!------------------------------------------------------------------------------
 CALL advectc ('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv)
 CALL acctimes ('ADVECTs')

!------------------------------------------------------------------------------
!  Update scalars
!  Skip the update call for SCM initial condition output.
!------------------------------------------------------------------------------
 if(iscm == 0 .or. (iscm >= 1 .and. scmtime /= 0.0))then
  CALL predtr ()          
  CALL acctimes ('PREDTR')
 endif

!##############################################################################
!SINGLE COLUMN MODEL (SCM) SECTION NEEDED FOR COMPARISON
!##############################################################################
!------------------------------------------------------------------------------
!  Output SCM column BEFORE call to physics packages
!------------------------------------------------------------------------------
 if(iscm==ngrid .and. time==scmtime .and. my_rams_num == my_scm_num) then
  CALL simdata ("parent_sim_info.txt","scm.in/")
  CALL readwrite_scm (2,"scm.in/")
 endif

!------------------------------------------------------------------------------
! Update radiation tendency to THP (theta-il) separately for including this
! tendency update to the SCM.
!------------------------------------------------------------------------------
 CALL predthp ()
 CALL acctimes ('PREDTHP')

!------------------------------------------------------------------------------
!  Moisture variables positive definite (negative adjustment rtp, etc,)
!  to be done following tendency updates and before call to microphysics.
!------------------------------------------------------------------------------
 CALL negadj1 () 
 CALL acctimes ('NEGADJ-1')

!------------------------------------------------------------------------------
!  Thermodynamic diagnosis (update theta,rv after negative adjustment and
!  after updates to THP, RTP, and hydrometeors from tendencies)
!------------------------------------------------------------------------------
 CALL thermo () 
 CALL acctimes ('THERMO')

!------------------------------------------------------------------------------
!  Aerosol Initial & Sources of dust and salt before microphysics call.
!------------------------------------------------------------------------------
 if (level <= 3) then
   CALL aerosols ()
   CALL acctimes ('AEROSOLS')
 endif

!------------------------------------------------------------------------------
!  Full 1 or 2 Moment Microphysics call.
!------------------------------------------------------------------------------
 if (level == 3) then
   CALL micro ()
 elseif (level == 4) then
   CALL micro_bin ()
 endif
 CALL acctimes ('MICRO')

!------------------------------------------------------------------------------
!  Apply scalar b.c.'s
!------------------------------------------------------------------------------
 CALL trsets ()          
 CALL acctimes ('TRSETS')

!------------------------------------------------------------------------------
!  Apply non-advected scalar top,bottom b.c.'s
!------------------------------------------------------------------------------
 CALL trsets_ns ()
 CALL acctimes ('TRSETS_NS')

!------------------------------------------------------------------------------
!  Moisture variables positive definite following microphysics update.
!------------------------------------------------------------------------------
 CALL negadj1 ()
 CALL acctimes ('NEGADJ-2')

!------------------------------------------------------------------------------
!  Thermodynamic diagnosis (update theta,rv after negative adjustment
!  and potentially before writing output after timestep is finished)
!------------------------------------------------------------------------------
 CALL thermo ()
 CALL acctimes ('THERMO')

!------------------------------------------------------------------------------
!  Check for negative micro and Nans
!------------------------------------------------------------------------------
 if(icheckmic == 1) then
  CALL checkmicro ('MICRO')
 endif

!------------------------------------------------------------------------------
!  Output SCM column AFTER call to physics packages
!------------------------------------------------------------------------------
 if(iscm==ngrid .and. time==scmtime .and. my_rams_num == my_scm_num) then
  CALL readwrite_scm (3,"scm.out/")
 endif
!##############################################################################
!##############################################################################

!------------------------------------------------------------------------------
!  KPP Ocean mixed-layer model
!------------------------------------------------------------------------------
 IF(IKPP > 0 ) CALL mckpp_ocean_mixing ()
 CALL acctimes ('KPP-Ocean')

!------------------------------------------------------------------------------
!  For computing FULL EXNER FUNCTION (after thermo update)
!  This is experimental. Do not use yet. (Nov 19, 2017)
!------------------------------------------------------------------------------
! CALL exevolve (mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,'THV',thvlast)
! CALL acctimes ('EXEVOLVE2')

!------------------------------------------------------------------------------
!  Lateral velocity boundaries - radiative
!------------------------------------------------------------------------------
 CALL latbnd ()         
 CALL acctimes ('LATBND')

!------------------------------------------------------------------------------
!  First stage Asselin filter
!------------------------------------------------------------------------------
 CALL hadvance (1)     
 CALL acctimes ('HADVANCE')

!------------------------------------------------------------------------------
!  Buoyancy term for w equation
!------------------------------------------------------------------------------
 CALL buoyancy ()
 CALL acctimes ('BUOYANCY')

!------------------------------------------------------------------------------
!  Acoustic small timesteps
!------------------------------------------------------------------------------
 CALL acoustic ()
 CALL acctimes ('ACOUSTIC')

!------------------------------------------------------------------------------
!  Last stage of Asselin filter
!------------------------------------------------------------------------------
 CALL hadvance (2)      
 CALL acctimes ('HADVANCE')

!------------------------------------------------------------------------------
!  Velocity/pressure boundary conditions
!------------------------------------------------------------------------------
 CALL vpsets ()          
 CALL acctimes ('VPSETS')

callmassflux=0    !flag for output BC mass flux: (0=off, 1==on)
massfluxfreq=300. !frequency of BC mass flux (seconds)
if(callmassflux==1) &
 CALL mass_flux_bc (mzp,mxp,myp &
      ,basic_g(ngrid)%up(1,1,1),basic_g(ngrid)%vp(1,1,1)  &
      ,basic_g(ngrid)%dn0(1,1,1) &
      ,basic_g(ngrid)%pp(1,1,1),basic_g(ngrid)%pi0(1,1,1) &
      ,grid_g(ngrid)%rtgu(1,1) ,grid_g(ngrid)%rtgv(1,1)    &
      ,grid_g(ngrid)%dyu(1,1)  ,grid_g(ngrid)%dyv(1,1))

return
END SUBROUTINE timestep

!##############################################################################
Subroutine acctimes (string)

use mem_all
use node_mod

implicit none

real, external :: valugp
integer :: ip,jp,kp,bp,patch
integer :: i, j, k
character(len=*) :: string

!Could call checkmicro here to narrow nans to a certain routine
if(icheckmic == 2) then
  CALL checkmicro (string)
endif

if(ngrid==1000 .and. string.eq.'VPSETS') then
do kp=29,13,-1
do jp=101,101
do ip=52,52
  i = ip - mi0(ngrid)
  j = jp - mj0(ngrid)
  k = kp
  if ((i.ge.1) .and. (i.le.mmxp(ngrid)) .and. &
      (j.ge.1) .and. (j.le.mmyp(ngrid))) then
   !if(micro_g(ngrid)%cccnp(kp,ip,jp)>631.0) &
   print'(a,3i5,2e17.9)','Steve',k,j,i,micro_g(ngrid)%rcp(k,i,j)
  endif
enddo
enddo
enddo
endif

!  only here for debugging purposes
kp=19
ip=25
jp=17
bp=1

101    format ('DEBUG: NODE',i0,': LOC(',i0,',',i0,',',i0,') --> (' &
   ,i0,',',i0,',',i0,') ',a10,': ',a10,f9.1,100e20.10)

if(ngrid==2000)then
do i=ia,iz
do j=ja,jz
 k=kp
 patch=2
  if(string .eq. 'MICRO')then
   print 101, my_rams_num, kp, ip, jp, k, i, j, 'NUCCLD', string, time &
     ,micro_g(ngrid)%nuccldrt(k,i,j)
   !print 101, my_rams_num, kp, ip, jp, k, i, j, 'FFCD', string, time &
   !  ,micro_g(ngrid)%ffcd(k,i,j,bp)
 endif
enddo
enddo
endif

if(ngrid==1000) then
  k = kp
  i = ip - mi0(ngrid)
  j = jp - mj0(ngrid)
  patch = 2
  if ((i.ge.1) .and. (i.le.mmxp(ngrid)) .and. &
      (j.ge.1) .and. (j.le.mmyp(ngrid))) then
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'THETA', string, time &
!     ,basic_g(ngrid)%theta(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'THP', string, time &
!     ,basic_g(ngrid)%thp(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'NUCCLD', string, time &
!     ,micro_g(ngrid)%nuccldrt(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'RTP', string, time &
!     ,basic_g(ngrid)%rtp(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'RV', string, time &
!     ,basic_g(ngrid)%rv(k,i,j)

   print 101, my_rams_num, kp, ip, jp, k, i, j, 'RRP', string, time &
     ,micro_g(ngrid)%rrp(k,i,j)
   print 101, my_rams_num, kp, ip, jp, k, i, j, 'rrt', string, time &
     ,valugp(mzp,mxp,myp,kp,ip,jp,tend%rrt(1))
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'RV', string, time &
!     ,basic_g(ngrid)%rv(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'THT', string, time &
!     ,valugp(mzp,mxp,myp,kp,ip,jp,tend%tht(1))
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'THETA', string, time &
!     ,basic_g(ngrid)%theta(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'flx_nsw', string, time &
!     ,kpp_3d_fields(ngrid)%flx_nsw(i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'khmx', string, time &
!     ,kpp_3d_fields(ngrid)%hmix(i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'cccnt', string, time &
!     ,valugp(mzp,mxp,myp,kp,ip,jp,tend%cccnt(1))
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'frtrd', string, time &
!     ,radiate_g(ngrid)%fthrd(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'frtrdp', string, time &
!     ,radiate_g(ngrid)%fthrdp(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'swdn', string, time &
!     ,radiate_g(ngrid)%swdn(k,i,j)
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'cifnt', string, time &
!     ,valugp(mzp,mxp,myp,kp,ip,jp,tend%cifnt(1))
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'WT', string, time &
!     ,valugp(mzp,mxp,myp,kp,ip,jp,tend%wt(1))
!   print 101, my_rams_num, kp, ip, jp, k, i, j, 'PATCH_AREA', string, time &
!     ,leaf_g(ngrid)%patch_area(i,j,patch)
  endif
endif

          !    , basic_g(ngrid)%up(1:mzp,ip,jp)
          !    , basic_g(ngrid)%pp(1:mzp,ip,jp)
          !    , micro_g(ngrid)%inuchazrt(kp,ip,jp)         &
          !    , micro_g(ngrid)%inucifnrt(kp,ip,jp)         &
          !    , radiate_g(ngrid)%rshort(ip,jp)             &
          !    , radiate_g(ngrid)%rlongup(ip,jp)            &
          !    , basic_g(ngrid)%wp_buoy_theta(kp,ip,jp)     &
          !    , basic_g(ngrid)%wp_buoy_cond(kp,ip,jp)      &
          !    , basic_g(ngrid)%wp(kp,ip,jp)                &
          !    , micro_g(ngrid)%rcp(kp,ip,jp)               &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%rct(1))   &
          !    , micro_g(ngrid)%ccp(kp,ip,jp)               &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%cct(1))   &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%ut(1))    &
          !    , basic_g(ngrid)%vp(kp,ip,jp)                &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%vt(1))    &
          !    , basic_g(ngrid)%thp(kp,ip,jp)               &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%vt(1))    &
          !    , basic_g(ngrid)%wp(kp,ip,jp)                &
          !    , valugp(mzp,mxp,myp,kp,ip,jp,tend%wt(1))    &
          !    , basic_g(ngrid)%pi0(kp,ip,jp)                 &
          !    , basicm_g(ngrid)%pi0(kp,ip,jp)         

return
END SUBROUTINE acctimes

!##############################################################################
Subroutine mass_flux_bc (m1,m2,m3,up,vp,dn0,pp,pi0,rtgu,rtgv,dyu,dxv)

use mem_grid
use mem_basic
use rconstants
use node_mod

implicit none

integer :: m1,m2,m3,i,j,k
real,dimension(m1,m2,m3) :: up,vp,dn0,pp,pi0
real,dimension(m2,m3) :: rtgu,rtgv,dyu,dxv
real :: wmass,emass,smass,nmass,prtot,tmass,ppp,area
real, save :: aintmass=0.

if(mod(time,300.).gt.0.1) return

if(nmachs .gt. 1)then
 print*,'Cannot compute boundary condition mass flux in parallel.'
 stop 'See TIMESTEP and MASS_FLUX_BC in rtimh.f90'
endif

wmass=0.
emass=0.
smass=0.
nmass=0.
prtot=0.
area=0.
!************************************************************
if(jdim==1) then  !for 3D simulations
!  west/east bound
do j=2,nyp-1
   do k=2,nzp-1
      i=1
      wmass=wmass +  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
      i=nxp-1
      emass=emass -  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
   enddo
enddo
!  north/south bound
do i=2,nxp-1
   do k=2,nzp-1
      j=1
      smass=smass +  &
           vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i,j+1))*.5
      j=nyp-1
      nmass=nmass -  &
           vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i,j+1))*.5
   enddo
enddo
k=2
do j=2,nyp-1
   do i=2,nxp-1
      ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
      prtot=prtot+ppp/(dyu(i,j)*dxv(i,j))
   enddo
enddo
area=(nxp-2)*deltax*(nyp-2)*deltax
endif

!************************************************************
if(jdim==0) then  !for 2D simulations
!  west/east bound
do j=1,1
   do k=2,nzp-1
      i=1
      wmass=wmass +  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
      i=nxp-1
      emass=emass -  &
           up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
           *(dn0(k,i,j)+dn0(k,i+1,j))*.5
   enddo
enddo
k=2
do j=1,1
   do i=2,nxp-1
      ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
      prtot=prtot+ppp/dyu(i,j)
   enddo
enddo
area=(nxp-2)*deltax
endif

tmass=wmass+emass+smass+nmass
aintmass=aintmass+tmass*dtlong
print*,'==============================='
print*,' Mass flux - W, E, S, N'
print*,  wmass,emass,smass,nmass
print*, 'total (kg/(m2 s)):   ',tmass/area
print*, 'total (kg/m2):       ',aintmass/area
print*, 'total pr change (pa):',aintmass/area*9.8
print*, 'computed mean press: ',prtot/area
print*,'==============================='

return
END SUBROUTINE mass_flux_bc
