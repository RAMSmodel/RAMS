!##############################################################################
Subroutine micro ()

use mem_basic
use mem_micro
use mem_grid
use mem_radiate
use mem_leaf
use node_mod
use rrad3
use micphys
use io_params

implicit none

integer :: ngr,lhcat,i,j,k
integer, dimension(11)  :: k1,k2,k3
integer, save :: ncall = 0
integer, save, dimension(maxgrds) :: ncallg
integer, save, dimension(16)  :: lcat0
real :: dtlti,sumtotal,sumcheck
data ncallg/maxgrds*0/
data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/

if (level .ne. 3) return

!******************************************************************************
! INITIALIZE MICROPHYSICS LOOKUP TABLES AND VARIABLES AT START ONLY
!******************************************************************************
!Run non-grid-dependent initialization at model startup for first timestep
if(ncall == 0) then
   ncall = 1
   !Initialize many microphysics parameters and lookup tables
   CALL micinit ()
   CALL setupF94 ()
   CALL make_autotab ()
   CALL haznuc ()
   CALL tabmelt ()
   CALL tabhab ()
   !Initialize several hydrometeor specific mass and fall speed coefficients
   do lhcat = 1,nhcat
      ch3(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
      cdp1(lhcat) = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
      pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
   enddo
endif

!Run grid-dependent initialization at model start for parent and nested grids
if(ncallg(ngrid) == 0) then
   ncallg(ngrid) = 1
   !Let user know when sedimentation tables are being created
   if(iprntstmt>=1 .and. print_msg) &
       print*,'Start Making Sedimentation Tables On (ngrid,dt):',ngrid,dtlt
   !Make grid-dependent sedimentation lookup table for quasi-bin implementation
   if(isedim==0) &
    CALL mksedim_tab (mzp,ngrid,zm,dzt  &
        ,pcp_tab(ngrid)%pcpfillc(1,1,1,1,1,1)  &
        ,pcp_tab(ngrid)%pcpfillr(1,1,1,1,1,1)  &
        ,pcp_tab(ngrid)%sfcpcp(1,1,1,1,1)      &
        ,pcp_tab(ngrid)%allpcp(1,1,1,1,1)      &
        ,dtlt)
   !Make grid-dependent sedimentation lookup table for true bin implementation
   if(isedim==1) &
    CALL mksedim_tab_trubin (mzp,zm,dzt  &
        ,pcp_tab(ngrid)%pcpfillc(1,1,1,1,1,1)  &
        ,pcp_tab(ngrid)%pcpfillr(1,1,1,1,1,1)  &
        ,pcp_tab(ngrid)%sfcpcp(1,1,1,1,1)      &
        ,pcp_tab(ngrid)%allpcp(1,1,1,1,1)      &
        ,dtlt)
   !Let user know when sedimentation tables finished being created
   if(iprntstmt>=1 .and. print_msg) &
       print*,'Done  Making Sedimentation Tables On (ngrid,dt):',ngrid,dtlt
   !Make some sedimentation variable for each hydrometeor category
   do lhcat = 1,nhcat
     if(isedim==0) &
       ch2(lhcat,ngrid) = float(nembfall-1) &
                      / log10(dispemb1(lhcat,ngrid) / dispemb0(lhcat,ngrid))
     if(isedim==1) &
       ch2(lhcat,ngrid) = float(nembfall-1) &
                      / log10(1.*emb1(lcat0(lhcat))/emb0(lcat0(lhcat)))
   enddo
   !Make lookup table for homogeneous freezing that is grid and time dependent
   CALL homfrzcl (dtlt,ngrid)
endif

!Get time-step dependent coefficients for current grid-nest and get initial
! mean mass if jnmb==2 (fixed diameter), and get jhcat arrays, and set some 
! initial arrays for vapor deposition.
CALL each_call (mzp,dtlt)
!Set inverse timestep to use as multiplier and avoid division by timestep
dtlti = 1. / dtlt

!******************************************************************************
! RUN THE MICROPHYSICS WITH IMBEDDED call TO HARRINGTON RADIATION
!******************************************************************************
!Going forth, ngr is used as abbreviation for ngrid
ngr = ngrid

!Zero out the radiative heating rate "fthrd" if this this a radiation timestep
!and if running Harrington Radiation as called below with microphysics loop.
if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then
  if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. .001) then
    CALL azero (mzp*mxp*myp,radiate_g(ngrid)%fthrd(1,1,1))
  endif
endif

!Loop through domain or parallel sub-domain
do j = ja,jz
   do i = ia,iz
      !Do micro flag check and copy global variabes to local variables
      CALL range_check (mzp,k1,k2,k3,i,j,frqstate(ngr),ngr,dtlt,time,micro_g(ngr))
      !Run the microphysics and Harrington radiation code 
      CALL mcphys (mzp,k1,k2,k3,i,j,ngr,maxnzp   &
         ,dtlt,dtlti,time,zm                     &
         ,zt,radiate_g(ngr)                      &
         ,basic_g(ngr)%up      (1,i,j)           &
         ,basic_g(ngr)%vp      (1,i,j)           &
         ,basic_g(ngr)%thp     (1,i,j)           &
         ,basic_g(ngr)%theta   (1,i,j)           &
         ,basic_g(ngr)%rtp     (1,i,j)           &
         ,basic_g(ngr)%rv      (1,i,j)           &
         ,basic_g(ngr)%wp      (1,i,j)           &
         ,grid_g(ngr)%rtgt     (i,j)             &
         ,micro_g(ngr)%pcpg    (i,j)             &
         ,micro_g(ngr)%qpcpg   (i,j)             &
         ,micro_g(ngr)%dpcpg   (i,j)             &
         ,pcp_tab(ngr)%pcpfillc(1,1,1,1,1,1)     &
         ,pcp_tab(ngr)%pcpfillr(1,1,1,1,1,1)     &
         ,pcp_tab(ngr)%sfcpcp  (1,1,1,1,1)       &
         ,pcp_tab(ngr)%allpcp  (1,1,1,1,1)       &
         ,grid_g(ngr)%glat     (i,j)             &
         ,grid_g(ngr)%topt     (i,j)             &
         ,basic_g(ngr)%dn0     (1,i,j)           &
         ,basic_g(ngr)%pi0     (1,i,j)           &
         ,basic_g(ngr)%pp      (1,i,j)           &
         ,leaf_g(ngr)%leaf_class(i,j,1:npatch)   &
         ,leaf_g(ngr)%patch_area(i,j,1:npatch)   &
         ,npatch                                 &
         !LEAF Variables needed for aerosol deposition
         ,leaf_g(ngr)%ustar(i,j,1:npatch)        &
         ,leaf_g(ngr)%patch_rough(i,j,1:npatch)  &
         ,imonth1                                &
         )
      !Copy local variables back to global variables
      CALL copyback (mzp,k2,k3,i,j,micro_g(ngr))
      !Quick way to view total accum precip
      !Good for testing in SEQUENTIAL mode - not for testing in PARALLEL
      sumcheck=0
      if(sumcheck==1)then
       if(i==ia.and.j==ja)sumtotal=0.0
       if(idriz>0)  sumtotal=sumtotal + micro_g(ngr)%accpd(i,j)
       if(irain>0)  sumtotal=sumtotal + micro_g(ngr)%accpr(i,j)
       if(ipris>0)  sumtotal=sumtotal + micro_g(ngr)%accpp(i,j)
       if(isnow>0)  sumtotal=sumtotal + micro_g(ngr)%accps(i,j)
       if(iaggr>0)  sumtotal=sumtotal + micro_g(ngr)%accpa(i,j)
       if(igraup>0) sumtotal=sumtotal + micro_g(ngr)%accpg(i,j)
       if(ihail>0)  sumtotal=sumtotal + micro_g(ngr)%accph(i,j)
       if(i==iz.and.j==jz)print*,"Grid total precip:",sumtotal
      endif
      !Write some output as need for testing and troubleshooting
      !if(i+mi0(ngr)==12 .and. j+mj0(ngr)==14) then
      ! do k=14,mzp,1
      !  print'(5i4,1f10.1,2X1f20.10)',k,i,ia,j,ja,time &
      !   ,micro_g(ngr)%q7(k,i,j)
      ! enddo
      !endif
   enddo
enddo

return
END SUBROUTINE micro

!##############################################################################
Subroutine mcphys (m1,k1,k2,k3,i,j,ngr,maxnzp          &
   ,dtlt,dtlti,time,zm,zt                              &
   ,radiate                                            &
   ,uup,vvp,thp,theta,rtp,rv,wp                        &
   ,rtgt,pcpg,qpcpg,dpcpg                              &
   ,pcpfillc,pcpfillr,sfcpcp,allpcp                    &
   ,glat,topt                                          &
   ,dn0                                                &
   ,pi0,pp                                             &
   !LEAF Variables needed for aerosol deposition onto surface
   ,lclass,parea,npatch                                &
   ,ustar,prough,imonthx                               &
   )

use mem_radiate
use rconstants
use rrad3
use micphys

implicit none

type (radiate_vars) :: radiate

integer :: i,j,k,lcat,jcat,icv,icx,mc1,mc2,mc3,mc4,m1  &
          ,ngr  &
          ,maxnzp,mcat  &
          ,k1cnuc,k2cnuc,k1dnuc,k2dnuc,k1pnuc,k2pnuc,lhcat

real,    dimension(8)   :: dpcp0
integer, dimension(8)   :: mcats,mivap,mix02
integer, dimension(9,4) :: mcat1
integer, dimension(8,2) :: mcat2,mcat22
integer, dimension(4)   :: mcat33
integer, dimension(11)  :: k1,k2,k3
real                    :: dtlt,dtlti,time
real                    :: rtgt,pcpg,qpcpg,dpcpg
real, dimension(m1)     :: zm,thp,theta,rtp,rv,wp,uup,vvp
real, dimension(m1)     :: dn0
real, dimension(m1)     :: pi0,pp

! Variables needed for Harrington radiation scheme
real                    :: glat,topt
real, dimension(m1)     :: zt

! Variables needed for hydrometeor sedimentation
real, dimension(m1,maxkfall,nembfall,nhcat,ndensrtgt,nband)::pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat,ndensrtgt,nband)::sfcpcp
real, dimension(m1,nembfall,nhcat,ndensrtgt,nband)::allpcp

! LEAF Variables needed for aerosol surface deposition
integer :: npatch
integer :: imonthx
real, dimension(npatch) :: lclass,parea
real, dimension(npatch) :: ustar,prough

! (mcats) is for rain, agg, graupel, hail self-collection (cols)
data mcats /0,3,0,0,6,7,10,0/     ! (effxy) number for coll efficiency
! (mcat1) is for ice-ice interactions (col1)
data mcat1 /3,3,3,4,4,4,5,5,6  &  ! 1st variable to interact
           ,5,6,7,5,6,7,6,7,7  &  ! 2nd variable to interact
           ,5,6,7,5,6,7,6,7,7  &  ! transfer variable for (r) & (q)
           ,4,7,8,5,7,8,7,8,8/    ! (effxy) number for coll efficiency
! (mcat2) is for cloud_1 - ice interactions (col2)
data mcat2 /0,0,0,6,6,7,7,0  &
           ,0,0,0,2,2,9,9,0/
! (mcat22) is for cloud_2 - ice interactions (col2)
data mcat22 /0,0,0,6,6,7,7,0  &
            ,0,0,0,11,11,12,12,0/
! (mcat33) is for pristine, snow self-collection (col3344)
data mcat33 /0,0,4,5/  ! (effxy) number for coll efficiency
! Order of vapor flux calculation
data mivap /1,8,3,4,5,2,6,7/
! Order of melting calculation and transfer
! though Cloud and Drizzle not considered
data mix02 /3,1,8,4,5,6,7,2/
! Multiplier in (sedim)
data dpcp0 /.001,.001,.010,.010,.010,.003,.001,.001/
save

!Output some data for diagnostic purposes
!if(i==4)then
! do k=2,38
! print*,'ice_start',k,(theta(k) * (pi0(k)+pp(k)) / 1004.) - 273.16 &
!   ,rx(k,3)+rx(k,4),(cx(k,3)+cx(k,4))*dn0(k)/1000.
! enddo
!endif

! Compute pressure, temperature, and moisture for vapor diffusion
 CALL thrmstr (m1,k1,k2,thp(1),theta(1),rtp(1),rv(1) &
              ,pp(1),pi0(1)                          &
              )

! Compute vapor diffusion terms
 CALL each_column (m1,k1,k2,rv(1),dn0(1))

! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.
do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      CALL enemb (m1,k1,k2,lcat,dn0(1))
   endif
enddo

! Evaluate radiative heating rates if using Harrington radiation scheme
if (iswrtyp .eq. 3 .or. ilwrtyp .eq. 3) then
   if (mod(time + .001,radfrq) .lt. dtlt .or. time .lt. .001) then
      !Saleeby(2008): Change passing of 7 to 8 if adding drizzle mode
      ! and modify locations in radcalc3 and radcomp3 to match
      CALL radcalc3 (m1,i,j,ngr,maxnzp,7,iswrtyp,ilwrtyp,zm,zt &
         ,glat,rtgt,topt,rv(1)    &
         ,radiate%albedt (i,j)    &
         ,radiate%cosz   (i,j)    &
         ,radiate%rlongup(i,j)    &
         ,radiate%rshort (i,j)    &
         ,radiate%rlong  (i,j)    &
         ,radiate%aodt   (i,j)    &
         ,radiate%fthrd(1,i,j)    &
         ,radiate%bext (1,i,j)    &
         ,radiate%swup (1,i,j)    &
         ,radiate%swdn (1,i,j)    &
         ,radiate%lwup (1,i,j)    &
         ,radiate%lwdn (1,i,j)    &
         ,dn0(1)                  &
         )
   endif
endif

! Save rx and qx before vapor diffusion.
rx_lhr = rx
qx_lhr = qx

! Pre-compute variables for vapor diffusion
do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      CALL diffprep (m1,lcat,k1(lcat),k2(lcat),rv(1),dn0(1))
   endif
enddo

! Pre-compute variables for vapor diffusion
 CALL vapdiff (m1,k1(11),k2(11),rv(1))

! Calculate vapor flux in order of species given by (mivap)
do icv = 1,8
   lcat = mivap(icv)
   if (jnmb(lcat) .ge. 1) then
      CALL vapflux (m1,lcat,k1(lcat),k2(lcat),rv(1),i,j)
   endif
enddo

! Pristine ice to snow transfer
if (jnmb(4) .ge. 1) then
   CALL psxfer (k1(3),k2(3),k1(4),k2(4),i,j)
endif

! Update temperature and moisture following vapor growth
 CALL newtemp (m1,k1(11),k2(11),rv(1),theta(1))

! Computed updated mean mass and weights for collision lookup tables
do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
      CALL enemb (m1,k1,k2,lcat,dn0(1))
      CALL getict (k1(lcat),k2(lcat),lcat)
   endif
enddo

! Add call to update LHR theta after vapor diffusion
 CALL calc_lhr_vap (k1,k2)

! Auto-accretion considered only for rain
if (jnmb(2) .ge. 1) then
   CALL auto_accret (m1,k1(1),k2(1),k1(8),k2(8),dn0(1),dtlt,i,j)
endif

! Binned riming of cloud species by the ice species
if(irime==1)then
 do jcat = 1,8,7
 do lcat = 4,7
  if (jnmb(jcat) .ge. 1 .and. jnmb(lcat) .ge. 1) then
   CALL auto_accret_ice (m1,jcat,lcat,k1(jcat),k2(jcat),dn0(1),dtlt)
  endif
 enddo
 enddo
endif

! Calculate collision/coalescence efficiencies before collection routines
 CALL effxy (m1,k1,k2)

! Self collection of rain, aggregates, graupel, hail:  number change only
do lcat = 2,7
   if (lcat .eq. 3 .or. lcat .eq. 4) go to 29
   mc1 = mcats(lcat)
   if (jnmb(lcat) >= 5) then
      CALL cols (lcat,mc1,k1(lcat),k2(lcat))
   endif
29 continue
enddo

! Self collection of pristine ice, snow (transfers to aggregates)
do lcat = 3,4
   mc1 = mcat33(lcat)
   if (jnmb(lcat) .ge. 1 .and. jnmb(5) .ge. 1) then
      CALL col3344 (lcat,5,mc1,k1(lcat),k2(lcat))
   endif
enddo

! Collection between pristine ice and snow
if (jnmb(5) .ge. 1) then
    CALL col3443 (3,4,5,max(k1(3),k1(4)),min(k2(3),k2(4)))
endif

! Ice-ice collisions
do icx = 1,9
   mc1 = mcat1(icx,1)
   mc2 = mcat1(icx,2)
   mc3 = mcat1(icx,3)
   mc4 = mcat1(icx,4)
   if (jnmb(mc1) .ge. 1 .and. jnmb(mc3) .ge. 1) then
      CALL col1 (mc1,mc2,mc3,mc4,max(k1(mc1),k1(mc2))  &
                ,min(k2(mc1),k2(mc2)))
   endif
enddo

! Ice-cloud collisions
if(irime==0) then
do jcat = 1,8,7
 do lcat = 4,7
   if(jcat==1)then
     mc1=mcat2(lcat,1)
     mc2=mcat2(lcat,2)
   endif
   if(jcat==8)then
     mc1=mcat22(lcat,1)
     mc2=mcat22(lcat,2)
   endif
   if (jnmb(jcat) .ge. 1 .and. jnmb(lcat).ge.1 .and. jnmb(mc1).ge.1) then
      CALL col2 (m1,jcat,lcat,mc1,mc2 &
    ,max(k1(jcat),k1(lcat)),min(k2(jcat),k2(lcat)),dn0(1),dtlt)
   endif
 enddo
enddo
endif

! Ice-rain collisions
do lcat = 3,7
   if (jnmb(2) .ge. 1 .and. jnmb(lcat) .ge. 1 .and. jnmb(7) .ge. 1) then
      CALL col3 (m1,2,lcat,max(k1(2),k1(lcat)),min(k2(2),k2(lcat)),dn0(1))
   endif
enddo

! Save rx and qx before collisions and melting...
! (colxfers is where mixing ratio is changed after collisions)
rx_lhr = rx
qx_lhr = qx

! Make hydrometeor transfers due to collision-coalescence
 CALL colxfers (m1,k1,k2,scrmic1,scrmic2)

! Pristine ice to snow transfer done after collision-coalescence to
! avoid any mass/number adjustments that impact cloud-ice number
if (jnmb(4) .ge. 1) then
   CALL psxfer (k1(3),k2(3),k1(4),k2(4),i,j)
endif

! Calcs r,q,c for each category considering melting processes
! in the order of pristine,cloud,drizzle,snow,agg,graupel,hail,rain
! though nothing done for cloud or drizzle in (x02) routine
do mcat = 1,8
   lcat = mix02(mcat)
   if (jnmb(lcat) .ge. 1) then
      CALL x02 (m1,k1,k2,lcat,dn0(1))
   endif
enddo

! Add call to update LHR theta after collisions and melting
 CALL calc_lhr_collmelt (m1)

! Save rx and qx before cloud nucleation...
rx_lhr = rx
qx_lhr = qx

! Calcs cloud mix ratio and concentration for cloud nucleation from CCN
if (jnmb(1) .ge. 1) then
   CALL cldnuc (m1,k1cnuc,k2cnuc,k1dnuc,k2dnuc,rv(1),wp(1),i,j,dn0(1))
endif

! Finds bottom and top layer of cloud water
if(jnmb(1) .ge. 1) then
 k1(1) = min(k1(1),k1cnuc)
 k2(1) = max(k2(1),k2cnuc)
 k3(1) = max(k2(1),k3(1))
endif

! Finds bottom and top layer of drizzle water
if(jnmb(8) .ge. 1) then
 k1(8) = min(k1(8),k1dnuc)
 k2(8) = max(k2(8),k2dnuc)
 k3(8) = max(k2(8),k3(8))
endif

! Keeps cloud and drizzle mass and number in bounds
if (jnmb(1) .ge. 3) CALL enemb (m1,k1,k2,1,dn0(1))
if (jnmb(8) .ge. 3) CALL enemb (m1,k1,k2,8,dn0(1))

! Add call to update LHR theta after cloud nucleation...
 CALL calc_lhr_cldnuc (k1,k2)

! Save rx and qx before ice nucleation.
rx_lhr = rx
qx_lhr = qx

! Ice nucleation
if (jnmb(3) .ge. 1) then
   CALL icenuc (m1,k1(1),k2(1),k1(8),k2(8),k1pnuc,k2pnuc,ngr,rv(1)  &
   ,dn0(1),dtlt,i,j)
endif

!Output some data for diagnostic purposes
!if(i==4)then
! do k=2,38
! print*,'ice_nuc',k,(theta(k) * (pi0(k)+pp(k)) / 1004.) - 273.16 &
!   ,rx(k,3)+rx(k,4),(cx(k,3)+cx(k,4))*dn0(k)/1000.
! enddo
!endif

! Finds bottom and top later of pristine ice
if (jnmb(3) .ge. 1) then
 k1(3) = min(k1(3),k1pnuc)
 k2(3) = max(k2(3),k2pnuc)
 k3(3) = max(k2(3),k3(3))
endif

! Keeps pristine, cloud, and drizzle mass and number in bounds
if (jnmb(3) .ge. 3) CALL enemb (m1,k1,k2,3,dn0(1))
if (jnmb(1) .ge. 3) CALL enemb (m1,k1,k2,1,dn0(1))
if (jnmb(8) .ge. 3) CALL enemb (m1,k1,k2,8,dn0(1))

!if(i==4)then
! do k=2,38
! print*,'ice_presed',k,(theta(k) * (pi0(k)+pp(k)) / 1004.) - 273.16 &
!   ,rx(k,3)+rx(k,4),(cx(k,3)+cx(k,4))*dn0(k)/1000.
! enddo
!endif

! Update latent heating budgets after ice nucleation
 CALL calc_lhr_icenuc (k1,k2)

! Zero out precip arrays that feed into land surface and soil
pcpg  = 0.
qpcpg = 0.
dpcpg = 0.

! tairc is used here to accumulate changes to thp from sedim
do k = 2,m1
   tairc(k) = 0.
enddo

! Compute hydrometeor specific coefficient that goes into sedimentation
do lhcat = 1,nhcat
   ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt
enddo

! Hydrometeor Sedimentation
do lcat = 1,8
   if (jnmb(lcat) .ge. 1) then
     if(isedim==0) &
      CALL sedim (m1,lcat,ngr  &
         ,k1(lcat),k2(lcat)  &
         ,rtp(1),thp(1),theta(1),dn0(1),dpcp0(lcat)  &
         ,pcpg,qpcpg,dpcpg,dtlti,scrmic1,scrmic2,scrmic3  &
         ,pcpfillc,pcpfillr,sfcpcp,allpcp)
     if(isedim==1) &
      CALL sedim_trubin (m1,lcat,ngr  &
         ,k1(lcat),k2(lcat)  &
         ,rtp(1),thp(1),theta(1),dn0(1),dpcp0(lcat)  &
         ,pcpg,qpcpg,dpcpg,dtlti,scrmic1,scrmic2,scrmic3  &
         ,pcpfillc,pcpfillr,sfcpcp,allpcp,rtgt)
   endif
enddo

! Update theta-il after microphysics is run
do k = 2,m1
   thp(k) = thp(k) + tairc(k)
enddo

! Dust and sea-salt dry and wet deposition driver
! applied after sedimentation for computation of precip rate
if(iaerodep==1) &
 CALL deposition_driver (i,j,m1,zm,rtgt               &
    ,rv(1),theta(1),uup(1),vvp(1)                     &
    ,dn0(1)                                           &
    ,pi0(1),pp(1)                                     &
    ,lclass(1),parea(1)                               &
    ,ustar(1),prough(1),imonthx                       &
    )

!Output some data for diagnostic purposes
!if(i==4)then
! do k=2,38
! print*,'ice_last',k,(theta(k) * (pi0(k)+pp(k)) / 1004.) - 273.16 &
!   ,rx(k,3)+rx(k,4),(cx(k,3)+cx(k,4))*dn0(k)/1000.
! enddo
!endif

return
END SUBROUTINE mcphys

!##############################################################################
Subroutine copyback (m1,k2,k3,i,j,micro)

use mem_micro
use micphys

implicit none

type (micro_vars) :: micro

integer, dimension(11)  :: k2,k3
integer :: m1,i,j,dtop,lcat

!Set drizzle top for 1 or 2 moment
if (jnmb(8) >= 1) then
   if(jnmb(8) >= 5) then
      dtop=k3(8)
   else
      dtop=k2(11)
   endif
endif

!Note that "cifnp" does not need a copyback since the IN concentration
!is not changed by microphysics when used.

!Copyback AEROSOLS
if (iaerosol > 0) then
   CALL ae1kmic (2,m1-1,micro%cn1np(1,i,j),aerocon(1,1))
   CALL ae1kmic (2,m1-1,micro%cn1mp(1,i,j),aeromas(1,1))
   CALL ae1kmic (2,m1-1,micro%cn2np(1,i,j),aerocon(1,2))
   CALL ae1kmic (2,m1-1,micro%cn2mp(1,i,j),aeromas(1,2))
endif
if (idust > 0) then
   CALL ae1kmic (2,m1-1,micro%md1np(1,i,j),aerocon(1,3))
   CALL ae1kmic (2,m1-1,micro%md1mp(1,i,j),aeromas(1,3))
   CALL ae1kmic (2,m1-1,micro%md2np(1,i,j),aerocon(1,4))
   CALL ae1kmic (2,m1-1,micro%md2mp(1,i,j),aeromas(1,4))
endif
if (isalt > 0) then
   CALL ae1kmic (2,m1-1,micro%salt_film_np(1,i,j),aerocon(1,5))
   CALL ae1kmic (2,m1-1,micro%salt_film_mp(1,i,j),aeromas(1,5))
   CALL ae1kmic (2,m1-1,micro%salt_jet_np(1,i,j) ,aerocon(1,6))
   CALL ae1kmic (2,m1-1,micro%salt_jet_mp(1,i,j) ,aeromas(1,6))
   CALL ae1kmic (2,m1-1,micro%salt_spum_np(1,i,j),aerocon(1,7))
   CALL ae1kmic (2,m1-1,micro%salt_spum_mp(1,i,j),aeromas(1,7))
endif
if (iabcarb > 0) then
   CALL ae1kmic (2,m1-1,micro%abc1np(1,i,j),aerocon(1,8))
   CALL ae1kmic (2,m1-1,micro%abc1mp(1,i,j),aeromas(1,8))
   CALL ae1kmic (2,m1-1,micro%abc2np(1,i,j),aerocon(1,9))
   CALL ae1kmic (2,m1-1,micro%abc2mp(1,i,j),aeromas(1,9))
endif

!Copyback AEROSOL TRACKING VARIABLES
if (iccnlev>=2) then
   if (jnmb(1) >= 1) then
    CALL ae1kmic (2,k3(1),micro%cnmcp(1,i,j),cnmhx(1,1))
    if(itrkepsilon==1) CALL ae1kmic (2,k3(1),micro%snmcp(1,i,j),snmhx(1,1))
    if(itrkdust==1)    CALL ae1kmic (2,k3(1),micro%dnmcp(1,i,j),dnmhx(1,1))
    if(itrkdustifn==1) CALL ae1kmic (2,k3(1),micro%dincp(1,i,j),dinhx(1,1))
   endif
   if (jnmb(2) >= 1) then
    CALL ae1kmic (2,k2(11),micro%cnmrp(1,i,j),cnmhx(1,2))
    if(itrkepsilon==1) CALL ae1kmic (2,k2(11),micro%snmrp(1,i,j),snmhx(1,2))
    if(itrkdust==1)    CALL ae1kmic (2,k2(11),micro%dnmrp(1,i,j),dnmhx(1,2))
    if(itrkdustifn==1) CALL ae1kmic (2,k2(11),micro%dinrp(1,i,j),dinhx(1,2))
   endif
   if (jnmb(3) >= 1) then
    CALL ae1kmic (2,k3(3),micro%cnmpp(1,i,j),cnmhx(1,3))
    if(itrkepsilon==1) CALL ae1kmic (2,k3(3),micro%snmpp(1,i,j),snmhx(1,3))
    if(itrkdust==1)    CALL ae1kmic (2,k3(3),micro%dnmpp(1,i,j),dnmhx(1,3))
    if(itrkdustifn==1) CALL ae1kmic (2,k3(3),micro%dinpp(1,i,j),dinhx(1,3))
   endif
   if (jnmb(4) >= 1) then
    CALL ae1kmic (2,k2(11),micro%cnmsp(1,i,j),cnmhx(1,4))
    if(itrkepsilon==1) CALL ae1kmic (2,k2(11),micro%snmsp(1,i,j),snmhx(1,4))
    if(itrkdust==1)    CALL ae1kmic (2,k2(11),micro%dnmsp(1,i,j),dnmhx(1,4))
    if(itrkdustifn==1) CALL ae1kmic (2,k2(11),micro%dinsp(1,i,j),dinhx(1,4))
   endif
   if (jnmb(5) >= 1) then
    CALL ae1kmic (2,k2(11),micro%cnmap(1,i,j),cnmhx(1,5))
    if(itrkepsilon==1) CALL ae1kmic (2,k2(11),micro%snmap(1,i,j),snmhx(1,5))
    if(itrkdust==1)    CALL ae1kmic (2,k2(11),micro%dnmap(1,i,j),dnmhx(1,5))
    if(itrkdustifn==1) CALL ae1kmic (2,k2(11),micro%dinap(1,i,j),dinhx(1,5))
   endif
   if (jnmb(6) >= 1) then
    CALL ae1kmic (2,k2(11),micro%cnmgp(1,i,j),cnmhx(1,6))
    if(itrkepsilon==1) CALL ae1kmic (2,k2(11),micro%snmgp(1,i,j),snmhx(1,6))
    if(itrkdust==1)    CALL ae1kmic (2,k2(11),micro%dnmgp(1,i,j),dnmhx(1,6))
    if(itrkdustifn==1) CALL ae1kmic (2,k2(11),micro%dingp(1,i,j),dinhx(1,6))
   endif
   if (jnmb(7) >= 1) then
    CALL ae1kmic (2,k2(11),micro%cnmhp(1,i,j),cnmhx(1,7))
    if(itrkepsilon==1) CALL ae1kmic (2,k2(11),micro%snmhp(1,i,j),snmhx(1,7))
    if(itrkdust==1)    CALL ae1kmic (2,k2(11),micro%dnmhp(1,i,j),dnmhx(1,7))
    if(itrkdustifn==1) CALL ae1kmic (2,k2(11),micro%dinhp(1,i,j),dinhx(1,7))
   endif
   if (jnmb(8) >= 1) then
    CALL ae1kmic (2,dtop,micro%cnmdp(1,i,j),cnmhx(1,8))
    if(itrkepsilon==1) CALL ae1kmic (2,dtop,micro%snmdp(1,i,j),snmhx(1,8))
    if(itrkdust==1)    CALL ae1kmic (2,dtop,micro%dnmdp(1,i,j),dnmhx(1,8))
    if(itrkdustifn==1) CALL ae1kmic (2,dtop,micro%dindp(1,i,j),dinhx(1,8))
   endif
   !Regenerated aerosol variables
   CALL ae1kmic (2,m1-1,micro%regen_aero1_np(1,i,j),aerocon(1,aerocat-1))
   CALL ae1kmic (2,m1-1,micro%regen_aero1_mp(1,i,j),aeromas(1,aerocat-1))
   CALL ae1kmic (2,m1-1,micro%regen_aero2_np(1,i,j),aerocon(1,aerocat))
   CALL ae1kmic (2,m1-1,micro%regen_aero2_mp(1,i,j),aeromas(1,aerocat))
   if(itrkepsilon==1) then
    CALL ae1kmic (2,m1-1,micro%resol_aero1_mp(1,i,j),regenmas(1,1))
    CALL ae1kmic (2,m1-1,micro%resol_aero2_mp(1,i,j),regenmas(1,2))
   endif
   !Total aerosol accumulation variables
   micro%pcpraero(i,j) = 0.0
   do lcat=1,ncat
    micro%accpaero(i,j) = micro%accpaero(i,j) + accpaerox(lcat)
    micro%pcpraero(i,j) = micro%pcpraero(i,j) + pcpraerox(lcat)
   enddo
   !Dust accumulation variables
   if(itrkdust==1) then
    micro%pcprdust(i,j) = 0.0
    do lcat=1,ncat
     micro%accpdust(i,j) = micro%accpdust(i,j) + accpdustx(lcat)
     micro%pcprdust(i,j) = micro%pcprdust(i,j) + pcprdustx(lcat)
    enddo
   endif
endif

!Copyback Immersion freezing nuclei variables
if (iifn==3 .and. iccnlev>=1) then
   if (jnmb(1) >= 5) CALL ae1kmic (2,m1-1,micro%ifnnucp(1,i,j),ifnnucx(1))
   if (jnmb(1) >= 5) CALL ae1kmic (2,k3(1),micro%immercp(1,i,j),immerhx(1,1))
   if (jnmb(2) >= 5) CALL ae1kmic (2,k2(11),micro%immerrp(1,i,j),immerhx(1,2))
   if (jnmb(8) >= 5) CALL ae1kmic (2,dtop,micro%immerdp(1,i,j),immerhx(1,8))
endif

if (jnmb(1) >= 1) then
   CALL ae1kmic (2,k3(1),micro%rcp(1,i,j),rx(1,1))
   if (jnmb(1) >= 5) CALL ae1kmic (2,k3(1),micro%ccp(1,i,j),cx(1,1))
endif

if (jnmb(2) >= 1) then
   CALL ae1kmic (2,k2(11),micro%rrp(1,i,j),rx(1,2))
   CALL ae1kmic (2,k2(11),micro%q2(1,i,j),qx(1,2))
   micro%accpr(i,j) = micro%accpr(i,j) + accpx(2)
   micro%pcprr(i,j) = pcprx(2)
   CALL ae1kmic (2,k2(11),micro%pcpvr(1,i,j),pcpvx(1,2))
   if (jnmb(2) >= 5) CALL ae1kmic (2,k2(11),micro%crp(1,i,j),cx(1,2))
endif

if (jnmb(3) >= 1) then
   CALL ae1kmic (2,k3(3),micro%rpp(1,i,j),rx(1,3))
   micro%accpp(i,j) = micro%accpp(i,j) + accpx(3)
   micro%pcprp(i,j) = pcprx(3)
   CALL ae1kmic (2,k3(3),micro%pcpvp(1,i,j),pcpvx(1,3))
   if (jnmb(3) >= 5) CALL ae1kmic (2,k3(3),micro%cpp(1,i,j),cx(1,3))
endif

if (jnmb(4) >= 1) then
   CALL ae1kmic (2,k2(11),micro%rsp(1,i,j),rx(1,4))
   micro%accps(i,j) = micro%accps(i,j) + accpx(4)
   micro%pcprs(i,j) = pcprx(4)
   CALL ae1kmic (2,k2(11),micro%pcpvs(1,i,j),pcpvx(1,4))
   if (jnmb(4) >= 5) CALL ae1kmic (2,k2(11),micro%csp(1,i,j),cx(1,4))
endif

if (jnmb(5) >= 1) then
   CALL ae1kmic (2,k2(11),micro%rap(1,i,j),rx(1,5))
   micro%accpa(i,j) = micro%accpa(i,j) + accpx(5)
   micro%pcpra(i,j) = pcprx(5)
   CALL ae1kmic (2,k2(11),micro%pcpva(1,i,j),pcpvx(1,5))
   if (jnmb(5) >= 5) CALL ae1kmic (2,k2(11),micro%cap(1,i,j),cx(1,5))
endif

if (jnmb(6) >= 1) then
   CALL ae1kmic (2,k2(11),micro%rgp(1,i,j),rx(1,6))
   CALL ae1kmic (2,k2(11),micro%q6(1,i,j),qx(1,6))
   micro%accpg(i,j) = micro%accpg(i,j) + accpx(6)
   micro%pcprg(i,j) = pcprx(6)
   CALL ae1kmic (2,k2(11),micro%pcpvg(1,i,j),pcpvx(1,6))
   if (jnmb(6) >= 5) CALL ae1kmic (2,k2(11),micro%cgp(1,i,j),cx(1,6))
endif

if (jnmb(7) >= 1) then
   CALL ae1kmic (2,k2(11),micro%rhp(1,i,j),rx(1,7))
   CALL ae1kmic (2,k2(11),micro%q7(1,i,j),qx(1,7))
   micro%accph(i,j) = micro%accph(i,j) + accpx(7)
   micro%pcprh(i,j) = pcprx(7)
   CALL ae1kmic (2,k2(11),micro%pcpvh(1,i,j),pcpvx(1,7))
   if (jnmb(7) >= 5) CALL ae1kmic (2,k2(11),micro%chp(1,i,j),cx(1,7)) 
endif

if (jnmb(8) >= 1) then
   CALL ae1kmic (2,dtop,micro%rdp(1,i,j),rx(1,8))
   micro%accpd(i,j) = micro%accpd(i,j) + accpx(8)
   micro%pcprd(i,j) = pcprx(8)
   CALL ae1kmic (2,dtop,micro%pcpvd(1,i,j),pcpvx(1,8))
   if (jnmb(8) >= 5) CALL ae1kmic (2,dtop,micro%cdp(1,i,j),cx(1,8))  
endif

!Copyback Microphysical Budget Variables
if(imbudget>=1) then
  CALL ae1kmic (1,m1,micro%latheatvap(1,i,j),xlatheatvap(1))
  CALL ae1kmic (1,m1,micro%latheatfrz(1,i,j),xlatheatfrz(1))
  CALL ae1kmic (1,m1,micro%nuccldrt(1,i,j),xnuccldrt(1))
  CALL ae1kmic (1,m1,micro%cld2raint(1,i,j),xcld2raint(1))
  CALL ae1kmic (1,m1,micro%ice2raint(1,i,j),xice2raint(1))
  CALL ae1kmic (1,m1,micro%nucicert(1,i,j),xnucicert(1))
  CALL ae1kmic (1,m1,micro%vapliqt(1,i,j),xvapliqt(1))
  CALL ae1kmic (1,m1,micro%vapicet(1,i,j),xvapicet(1))
  CALL ae1kmic (1,m1,micro%evapliqt(1,i,j),xevapliqt(1))
  CALL ae1kmic (1,m1,micro%evapicet(1,i,j),xevapicet(1))
  CALL ae1kmic (1,m1,micro%freezingt(1,i,j),xfreezingt(1))
  CALL ae1kmic (1,m1,micro%meltingt(1,i,j),xmeltingt(1))
  CALL ae1kmic (1,m1,micro%melticet(1,i,j),xmelticet(1))
  CALL ae1kmic (1,m1,micro%rimecldt(1,i,j),xrimecldt(1))
  CALL ae1kmic (1,m1,micro%rain2icet(1,i,j),xrain2icet(1))
  CALL ae1kmic (1,m1,micro%aggregatet(1,i,j),xaggregatet(1))
  CALL ae1kmic (1,m1,micro%latheatvapt(1,i,j),xlatheatvapt(1))
  CALL ae1kmic (1,m1,micro%latheatfrzt(1,i,j),xlatheatfrzt(1))
endif
if(imbudget>=2) then
  CALL ae1kmic (1,m1,micro%inuchomrt(1,i,j),xinuchomrt(1))
  CALL ae1kmic (1,m1,micro%inuccontrt(1,i,j),xinuccontrt(1))
  CALL ae1kmic (1,m1,micro%inucifnrt(1,i,j),xinucifnrt(1))
  CALL ae1kmic (1,m1,micro%inuchazrt(1,i,j),xinuchazrt(1))
  CALL ae1kmic (1,m1,micro%vapcldt(1,i,j),xvapcldt(1))
  CALL ae1kmic (1,m1,micro%vapraint(1,i,j),xvapraint(1))
  CALL ae1kmic (1,m1,micro%vapprist(1,i,j),xvapprist(1))
  CALL ae1kmic (1,m1,micro%vapsnowt(1,i,j),xvapsnowt(1))
  CALL ae1kmic (1,m1,micro%vapaggrt(1,i,j),xvapaggrt(1))
  CALL ae1kmic (1,m1,micro%vapgraut(1,i,j),xvapgraut(1))
  CALL ae1kmic (1,m1,micro%vaphailt(1,i,j),xvaphailt(1))
  CALL ae1kmic (1,m1,micro%vapdrizt(1,i,j),xvapdrizt(1))
  CALL ae1kmic (1,m1,micro%evapcldt(1,i,j),xevapcldt(1))
  CALL ae1kmic (1,m1,micro%evapraint(1,i,j),xevapraint(1))
  CALL ae1kmic (1,m1,micro%evapprist(1,i,j),xevapprist(1))
  CALL ae1kmic (1,m1,micro%evapsnowt(1,i,j),xevapsnowt(1))
  CALL ae1kmic (1,m1,micro%evapaggrt(1,i,j),xevapaggrt(1))
  CALL ae1kmic (1,m1,micro%evapgraut(1,i,j),xevapgraut(1))
  CALL ae1kmic (1,m1,micro%evaphailt(1,i,j),xevaphailt(1))
  CALL ae1kmic (1,m1,micro%evapdrizt(1,i,j),xevapdrizt(1))
  CALL ae1kmic (1,m1,micro%meltprist(1,i,j),xmeltprist(1))
  CALL ae1kmic (1,m1,micro%meltsnowt(1,i,j),xmeltsnowt(1))
  CALL ae1kmic (1,m1,micro%meltaggrt(1,i,j),xmeltaggrt(1))
  CALL ae1kmic (1,m1,micro%meltgraut(1,i,j),xmeltgraut(1))
  CALL ae1kmic (1,m1,micro%melthailt(1,i,j),xmelthailt(1))
  CALL ae1kmic (1,m1,micro%rimecldsnowt(1,i,j),xrimecldsnowt(1))
  CALL ae1kmic (1,m1,micro%rimecldaggrt(1,i,j),xrimecldaggrt(1))
  CALL ae1kmic (1,m1,micro%rimecldgraut(1,i,j),xrimecldgraut(1))
  CALL ae1kmic (1,m1,micro%rimecldhailt(1,i,j),xrimecldhailt(1))
  CALL ae1kmic (1,m1,micro%rain2prt(1,i,j),xrain2prt(1))
  CALL ae1kmic (1,m1,micro%rain2snt(1,i,j),xrain2snt(1))
  CALL ae1kmic (1,m1,micro%rain2agt(1,i,j),xrain2agt(1))
  CALL ae1kmic (1,m1,micro%rain2grt(1,i,j),xrain2grt(1))
  CALL ae1kmic (1,m1,micro%rain2hat(1,i,j),xrain2hat(1))
  CALL ae1kmic (1,m1,micro%aggrselfprist(1,i,j),xaggrselfprist(1))
  CALL ae1kmic (1,m1,micro%aggrselfsnowt(1,i,j),xaggrselfsnowt(1))
  CALL ae1kmic (1,m1,micro%aggrprissnowt(1,i,j),xaggrprissnowt(1))
endif
if(imbudget==3 .and. idust>=1) then
  CALL ae1kmic (1,m1,micro%dust1cldrt(1,i,j),xdust1cldrt(1))
  CALL ae1kmic (1,m1,micro%dust2cldrt(1,i,j),xdust2cldrt(1))
  CALL ae1kmic (1,m1,micro%dust1drzrt(1,i,j),xdust1drzrt(1))
  CALL ae1kmic (1,m1,micro%dust2drzrt(1,i,j),xdust2drzrt(1))
endif

!Set Micro Budgets bottom level with first level above ground
if(imbudget>=1) then
  micro%latheatvap(1,i,j)  = micro%latheatvap(2,i,j)
  micro%latheatfrz(1,i,j)  = micro%latheatfrz(2,i,j)
  micro%nuccldrt(1,i,j)    = micro%nuccldrt(2,i,j)
  micro%cld2raint(1,i,j)   = micro%cld2raint(2,i,j)
  micro%ice2raint(1,i,j)   = micro%ice2raint(2,i,j)
  micro%nucicert(1,i,j)    = micro%nucicert(2,i,j)
  micro%vapliqt(1,i,j)     = micro%vapliqt(2,i,j)
  micro%vapicet(1,i,j)     = micro%vapicet(2,i,j)
  micro%evapliqt(1,i,j)    = micro%evapliqt(2,i,j)
  micro%evapicet(1,i,j)    = micro%evapicet(2,i,j)
  micro%freezingt(1,i,j)   = micro%freezingt(2,i,j)
  micro%meltingt(1,i,j)    = micro%meltingt(2,i,j)
  micro%melticet(1,i,j)    = micro%melticet(2,i,j)
  micro%rimecldt(1,i,j)    = micro%rimecldt(2,i,j)
  micro%rain2icet(1,i,j)   = micro%rain2icet(2,i,j)
  micro%aggregatet(1,i,j)  = micro%aggregatet(2,i,j)
  micro%latheatvapt(1,i,j) = micro%latheatvapt(2,i,j)
  micro%latheatfrzt(1,i,j) = micro%latheatfrzt(2,i,j)
endif
if(imbudget>=2) then
  micro%inuchomrt(1,i,j)     = micro%inuchomrt(2,i,j)
  micro%inuccontrt(1,i,j)    = micro%inuccontrt(2,i,j)
  micro%inucifnrt(1,i,j)     = micro%inucifnrt(2,i,j)
  micro%inuchazrt(1,i,j)     = micro%inuchazrt(2,i,j)
  micro%vapcldt(1,i,j)       = micro%vapcldt(2,i,j)
  micro%vapraint(1,i,j)      = micro%vapraint(2,i,j)
  micro%vapprist(1,i,j)      = micro%vapprist(2,i,j)
  micro%vapsnowt(1,i,j)      = micro%vapsnowt(2,i,j)
  micro%vapaggrt(1,i,j)      = micro%vapaggrt(2,i,j)
  micro%vapgraut(1,i,j)      = micro%vapgraut(2,i,j)
  micro%vaphailt(1,i,j)      = micro%vaphailt(2,i,j)
  micro%vapdrizt(1,i,j)      = micro%vapdrizt(2,i,j)
  micro%evapcldt(1,i,j)      = micro%evapcldt(2,i,j)
  micro%evapraint(1,i,j)     = micro%evapraint(2,i,j)
  micro%evapprist(1,i,j)     = micro%evapprist(2,i,j)
  micro%evapsnowt(1,i,j)     = micro%evapsnowt(2,i,j)
  micro%evapaggrt(1,i,j)     = micro%evapaggrt(2,i,j)
  micro%evapgraut(1,i,j)     = micro%evapgraut(2,i,j)
  micro%evaphailt(1,i,j)     = micro%evaphailt(2,i,j)
  micro%evapdrizt(1,i,j)     = micro%evapdrizt(2,i,j)
  micro%meltprist(1,i,j)     = micro%meltprist(2,i,j)
  micro%meltsnowt(1,i,j)     = micro%meltsnowt(2,i,j)
  micro%meltaggrt(1,i,j)     = micro%meltaggrt(2,i,j)
  micro%meltgraut(1,i,j)     = micro%meltgraut(2,i,j)
  micro%melthailt(1,i,j)     = micro%melthailt(2,i,j)
  micro%rimecldsnowt(1,i,j)  = micro%rimecldsnowt(2,i,j)
  micro%rimecldaggrt(1,i,j)  = micro%rimecldaggrt(2,i,j)
  micro%rimecldgraut(1,i,j)  = micro%rimecldgraut(2,i,j)
  micro%rimecldhailt(1,i,j)  = micro%rimecldhailt(2,i,j)
  micro%rain2prt(1,i,j)      = micro%rain2prt(2,i,j)
  micro%rain2snt(1,i,j)      = micro%rain2snt(2,i,j)
  micro%rain2agt(1,i,j)      = micro%rain2agt(2,i,j)
  micro%rain2grt(1,i,j)      = micro%rain2grt(2,i,j)
  micro%rain2hat(1,i,j)      = micro%rain2hat(2,i,j)
  micro%aggrselfprist(1,i,j) = micro%aggrselfprist(2,i,j)
  micro%aggrselfsnowt(1,i,j) = micro%aggrselfsnowt(2,i,j)
  micro%aggrprissnowt(1,i,j) = micro%aggrprissnowt(2,i,j)
endif
if(imbudget==3 .and. idust>=1) then
  micro%dust1cldrt(1,i,j)     = micro%dust1cldrt(2,i,j)
  micro%dust2cldrt(1,i,j)     = micro%dust2cldrt(2,i,j)
  micro%dust1drzrt(1,i,j)     = micro%dust1drzrt(2,i,j)
  micro%dust2drzrt(1,i,j)     = micro%dust2drzrt(2,i,j)
endif

return
END SUBROUTINE copyback

!##############################################################################
Subroutine calc_lhr_vap (k1,k2)

! This routine computes the change to potential temperature associated
! with all of the RAMS microphysics vapor diffusion processes. It assumes
! that the heat storage and mixing ratios were stored before calls to vapor
! diffusion and accounts for fractional liquid/ice on hail and graupel. It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between liquid and
! ice diffusion, nor does it keep track of condensation vs. evaporation.

use rconstants
use micphys

implicit none

integer :: lcat,k
integer, dimension(11) :: k1,k2
real :: temp,fracliq1,fracliq2,latheat,rxchange

if(imbudget>=1)then

 ! For each category, compute the change due to vapor transfer and apply to theta
 do lcat=1,ncat
  if (jnmb(lcat) .ge. 1) then
    do k=k1(lcat),k2(lcat)

      ! Vapor to liquid
      if (lcat .eq. 1 .or. lcat .eq. 2 .or. lcat .eq. 8) then
        if (lhrtheta) then
            latheat =(alvl/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat)) 
          else
            latheat = alvl * cpi * (rx(k,lcat)-rx_lhr(k,lcat))
        endif
        xlatheatvap(k) = xlatheatvap(k) + latheat
        xlatheatvapt(k) = xlatheatvapt(k) + latheat
      endif

      ! Vapor to ice
      if (lcat .eq. 3 .or. lcat .eq. 4 .or. lcat .eq. 5) then
        if (lhrtheta) then
            latheat =(alvi/pitot(k))*(rx(k,lcat)-rx_lhr(k,lcat)) 
          else
            latheat = alvi * cpi *(rx(k,lcat)-rx_lhr(k,lcat)) 
        endif
        xlatheatvap(k) = xlatheatvap(k) + latheat
        xlatheatvapt(k) = xlatheatvapt(k) + latheat
      endif

      ! Account for fractional liquid in graupel and hail categories
      if (lcat .eq. 6 .or. lcat .eq. 7) then
        ! Get fraction of liquid from old values
        CALL qtc (qx_lhr(k,lcat),temp,fracliq1)
        ! Get fraction of liquid from new values
        CALL qtc (qx(k,lcat),    temp,fracliq2)

        ! First calculate vapor depsition 
        if (lhrtheta) then
            latheat = (alvi/pitot(k)) * &
              (rx(k,lcat)*(1.-fracliq1)-rx_lhr(k,lcat)*(1.-fracliq1))
        else
            latheat = alvi * cpi * &
              (rx(k,lcat)*(1.-fracliq1)-rx_lhr(k,lcat)*(1.-fracliq1))
        endif
        xlatheatvap(k) = xlatheatvap(k) + latheat 
        xlatheatvapt(k) = xlatheatvapt(k) + latheat

        ! Now calculate melting
        if (lhrtheta) then
            latheat = (alli/pitot(k)) * &
              (rx(k,lcat)*(1.-fracliq2)-rx(k,lcat)*(1.-fracliq1))
        else
            latheat = alli * cpi * &
              (rx(k,lcat)*(1.-fracliq2)-rx(k,lcat)*(1.-fracliq1))
        endif
        xlatheatfrz(k) = xlatheatfrz(k) + latheat
        xlatheatfrzt(k) = xlatheatfrzt(k) + latheat

        ! Compute total freezing or melting and enter as positive values 
        rxchange = (rx(k,lcat)*(1.-fracliq2)-rx(k,lcat)*(1.-fracliq1))
        if (rxchange >= 0.0) then
          xfreezingt(k) = xfreezingt(k) + rxchange
        else
          xmeltingt(k)  = xmeltingt(k)  - rxchange
        endif
      endif

    enddo
  endif
 enddo

endif

return
END SUBROUTINE calc_lhr_vap

!##############################################################################
Subroutine calc_lhr_collmelt (m1)

! This routine computes the change to potential temperature associated
! with all of the RAMS microphysics collision freezing/melting processes. It 
! assumes that the heat storage and mixing ratios were stored before calls to 
! colxfers and accounts for fractional liquid/ice on hail and graupel. It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between freezing and
! melting--it only computes the net effect of all freezing/melting processes
! on the potential temperature.

use rconstants
use micphys

implicit none

integer :: m1,lcat,k
real :: fracliq1, fracliq2, temp
real, dimension(m1) :: rliq1,rice1,rliq2,rice2

if(imbudget>=1)then

 ! Compute total liquid and ice amounts before and after collisions/melting
 rliq1 = 0.
 rice1 = 0.
 rliq2 = 0.
 rice2 = 0.
 do lcat=1,ncat
  if (jnmb(lcat) .ge. 1) then
    do k=2,m1
      if (lcat.eq.1 .or. lcat.eq.2 .or. lcat.eq.8) then ! liquid only
        rliq1(k) = rliq1(k) + rx_lhr(k,lcat)
        rliq2(k) = rliq2(k) + rx(k,lcat)
      elseif (lcat.eq.3 .or. lcat.eq.4 .or. lcat.eq.5) then ! ice only
        rice1(k) = rice1(k) + rx_lhr(k,lcat)
        rice2(k) = rice2(k) + rx(k,lcat)
      elseif (lcat.eq.6 .or. lcat.eq.7) then !mixed phase
        !Get liquid fractions
        CALL qtc (qx_lhr(k,lcat),temp,fracliq1)
        CALL qtc (qx(k,lcat),    temp,fracliq2)
        !Accumulate liquid and ice
        rliq1(k) = rliq1(k) + rx_lhr(k,lcat) * fracliq1
        rliq2(k) = rliq2(k) + rx(k,lcat)     * fracliq2
        rice1(k) = rice1(k) + rx_lhr(k,lcat) * (1.-fracliq1)
        rice2(k) = rice2(k) + rx(k,lcat)     * (1.-fracliq2)      
      endif
    enddo
  endif
 enddo

 ! Heating/cooling is due to change in net liquid vs. ice mass only
 do k=2,m1
   if (lhrtheta) then
    xlatheatfrz(k) = xlatheatfrz(k) + (alli/pitot(k))*(rice2(k)-rice1(k))
    xlatheatfrzt(k) = xlatheatfrzt(k) + (alli/pitot(k))*(rice2(k)-rice1(k))
   else
    xlatheatfrz(k) = xlatheatfrz(k) +  alli * cpi * (rice2(k)-rice1(k))
    xlatheatfrzt(k) = xlatheatfrzt(k) +  alli * cpi * (rice2(k)-rice1(k))
   endif
   ! Compute total freezing or melting and enter as positive values
   if (rice2(k)-rice1(k) >= 0.0) then
     xfreezingt(k) = xfreezingt(k) + (rice2(k)-rice1(k))
   else
     xmeltingt(k)  = xmeltingt(k)  - (rice2(k)-rice1(k))
   endif
 enddo

endif

return
END SUBROUTINE calc_lhr_collmelt

!##############################################################################
Subroutine calc_lhr_cldnuc (k1,k2)

! This routine computes the change to potential temperature associated
! with cloud nucleation in RAMS microphysics. It assumes that the mass mixing
! ratios were stored before calls to cldnuc (hail/graupel are not considered
! here, and so there is no need to store the heat storage vectors). It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between nucleation and
! evaporation--it only computes the net effect of the nucleation processes
! on the potential temperature.

use rconstants
use micphys

implicit none

integer :: k
integer, dimension(11) :: k1,k2

if(imbudget>=1)then

 ! Cloud nucleation only changes the mixing ratio 
 ! (and number concentration) of cloud1 and cloud2
 if (jnmb(1) .ge. 1) then
  do k=k1(1),k2(1)
    if (lhrtheta) then
      xlatheatvap(k) = xlatheatvap(k) + (alvl/pitot(k))*(rx(k,1)-rx_lhr(k,1))
      xlatheatvapt(k) = xlatheatvapt(k) + (alvl/pitot(k))*(rx(k,1)-rx_lhr(k,1))
    else
      xlatheatvap(k) = xlatheatvap(k) +  alvl * cpi * (rx(k,1)-rx_lhr(k,1))
      xlatheatvapt(k) = xlatheatvapt(k) +  alvl * cpi * (rx(k,1)-rx_lhr(k,1))
    endif
  enddo
 endif
 if (jnmb(8) .ge. 1) then
  do k=k1(8),k2(8)
    if (lhrtheta) then
      xlatheatvap(k) = xlatheatvap(k) + (alvl/pitot(k))*(rx(k,8)-rx_lhr(k,8))
      xlatheatvapt(k) = xlatheatvapt(k) + (alvl/pitot(k))*(rx(k,8)-rx_lhr(k,8))
    else
      xlatheatvap(k) = xlatheatvap(k) +  alvl * cpi * (rx(k,8)-rx_lhr(k,8))
      xlatheatvapt(k) = xlatheatvapt(k) +  alvl * cpi * (rx(k,8)-rx_lhr(k,8))
    endif
  enddo
 endif

endif

return
END SUBROUTINE calc_lhr_cldnuc

!##############################################################################
Subroutine calc_lhr_icenuc (k1,k2)

! This routine computes the change to potential temperature associated
! with ice nucleation in RAMS microphysics. It assumes that the mass mixing
! ratios were stored before calls to icenuc (hail/graupel are not considered
! here, and so there is no need to store the heat storage vectors). It uses
! the "jnmb" vector to determine whether changes to a given species are
! predicted. At this time, it does not differentiate between nucleation and
! sublimation--it only computes the net effect of the nucleation processes
! on the potential temperature.

use rconstants
use micphys

implicit none

integer :: k
integer, dimension(11) :: k1,k2
real, dimension(nzpmax) :: vapdep

if(imbudget>=1)then

 ! Ice nucleation only changes the mixing ratio 
 ! (and number concentration) of pristine ice
 if (jnmb(3) .ge. 1) then
  ! However, some comes from contact nucleation of cloud1 and cloud2, 
  ! and some from vapor deposition...
  ! These need to be treated seperately due to different latent heats...
  
  ! First, changes in cloud1 and cloud2 (contact nucleation)
  if (jnmb(1) .ge. 1) then
    do k=k1(1),k2(1) ! Loss of cloud1
      if (lhrtheta) then
        xlatheatvap(k) = xlatheatvap(k) - (alli/pitot(k))*(rx(k,1)-rx_lhr(k,1))
        xlatheatvapt(k) = xlatheatvapt(k) - (alli/pitot(k))*(rx(k,1)-rx_lhr(k,1))
      else
        xlatheatvap(k) = xlatheatvap(k) - alli * cpi * (rx(k,1)-rx_lhr(k,1))
        xlatheatvapt(k) = xlatheatvapt(k) - alli * cpi * (rx(k,1)-rx_lhr(k,1))
      endif
    enddo
  endif
  if (jnmb(8) .ge. 1) then
    do k=k1(8),k2(8) ! Loss of cloud2
      if (lhrtheta) then
        xlatheatvap(k) = xlatheatvap(k) - (alli/pitot(k))*(rx(k,8)-rx_lhr(k,8))
        xlatheatvapt(k) = xlatheatvapt(k) - (alli/pitot(k))*(rx(k,8)-rx_lhr(k,8))
      else
        xlatheatvap(k) = xlatheatvap(k) - alli * cpi * (rx(k,8)-rx_lhr(k,8))
        xlatheatvapt(k) = xlatheatvapt(k) - alli * cpi * (rx(k,8)-rx_lhr(k,8))
      endif
    enddo
  endif

  ! Then, changes in pristine ice minus changes in 
  ! cloud1 and cloud2 (vapor deposition)
  do k=k1(3),k2(3)
    ! First, compute total change in pristine ice (should be positive)
    vapdep(k) = (rx(k,3)-rx_lhr(k,3))
    ! Then, subtract off contributions from cloud1 and cloud2
    if (jnmb(1) .ge. 1) vapdep(k) = vapdep(k) - (rx_lhr(k,1)-rx(k,1))
    if (jnmb(8) .ge. 1) vapdep(k) = vapdep(k) - (rx_lhr(k,8)-rx(k,8))
    ! Finally, compute change to theta due to vapor deposition...
    if (lhrtheta) then
      xlatheatvap(k) = xlatheatvap(k) + (alvi/pitot(k)) * vapdep(k)
      xlatheatvapt(k) = xlatheatvapt(k) + (alvi/pitot(k)) * vapdep(k)
    else
      xlatheatvap(k) = xlatheatvap(k) +  alvi * cpi * vapdep(k)
      xlatheatvapt(k) = xlatheatvapt(k) +  alvi * cpi * vapdep(k)
    endif
  enddo

 endif ! End test for whether pristine ice is predicted

endif

return
END SUBROUTINE calc_lhr_icenuc
