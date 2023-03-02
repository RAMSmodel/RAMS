!##############################################################################
Module micphys

use grid_dims

implicit none

!--------------------------------------------------------------------------
!     The product [(nthz-1)  * dthz ] must equal 25.0.
!     The product [(nrhhz-1) * drhhz] must equal 0.18.
!     The product [(ntc-1)   * dtc  ] must equal 20.0.
!     The product [(ndnc-1)  * ddnc ] must equal 20.e-6.

integer, parameter :: nthz=26,nrhhz=10,ngam=5000,ninc=201   &
                     ,ndns=15,ntc=21,ndnc=11                &
                     ,ndcc=60,ndcd=20,ndccr=15,nrrcr=30     &
                     ,ndrcr=30                              &
                     ,ncat=8,nhcat=16,npairc=101,npairr=147 &
                     ,nembc=20
real, parameter    :: dtc=1.,ddnc=2.e-6 ,dthz=1.,drhhz=.02
real, parameter    :: budget_scalet=1.
real, parameter    :: rxmin=1.e-9,cxmin=1.e-3

!IDIFFPERTS
!0=normal scalar diffusion
!1=diffuse perturbations THP,RTP from 3D base state
!2=diffuse perturbations THP,RTP from 3D current domain-mean state
!3=diffuse perturbations THP,RTP from 3D varfile state
integer :: idiffperts

!--------------------------------------------------------------------------
integer :: level,icloud,idriz,irain,ipris,isnow,iaggr,igraup,ihail      &
  ,irime,iplaws,iaerosol,idust,idustloft,iabcarb,isalt,iaerorad,iifn    &
  ,imbudget,isedim,itrkepsilon,itrkdust,itrkdustifn,iaerodep,icheckmic  &
  ,iaeroprnt,iaerohist,iifn_formula,iscm,iscmx,iscmy

integer, dimension(maxgrds) :: iaerolbc,ico2lbc
real, dimension(maxgrds) :: bctau

integer, dimension(ncat)        :: jnmb
integer, dimension(nhcat,nhcat) :: ipairc,ipairr
integer, dimension(31,100,2)    :: jhabtab
integer, dimension(nzpmax,ncat) :: jhcat,ict1,ict2

logical, parameter :: lhrtheta = .true.

real :: cparm,dparm,rparm,pparm,sparm,aparm,gparm,hparm         &
       ,rictmin,rictmax,dps,dps2                                &
       ,d1min,d1max,d2min,d2max,d3min,d3minx,d3max,r3min,r3max  &
       ,d1ecr,d2ecr,d3ecr,r3ecr                                 &
       ,colf,pi4dt,sedtime0,sedtime1                            &
       ,dimin,diminx,dimax,rimin,rimax,dieci,rieci,scmtime

real, dimension(ncat)  :: emb0,emb1,gnu,parm,emb0log,emb1log,dict
real, dimension(nhcat) :: shape,cfmas,pwmas,cfvt,pwvt,dpsmi,cfden,pwden  &
                         ,cfemb0,cfen0,pwemb0,pwen0,frefac1,frefac2  &
                         ,dnfac,sipfac,pwmasi,ch1,ch3,cdp1,pwvtmasi
real, dimension(nzpmax) :: tair,tairc,tairstrc,til,rvstr,press,pitot  &
                          ,rliq,rice,qhydm,rvlsair,rvisair,rvs0,thrmcon  &
                          ,vapdif,dynvisc,rdynvsci,denfac,dn0i,colfacr  &
                          ,colfacr2,colfacc,colfacc2,sumuy,sumuz,sumvr  &
                          ,scrmic1,scrmic2,scrmic3
real, dimension(nzpmax,ncat) :: rx,cx,qr,qx,tx,pcpvx,emb,vap,ttest,wct1  &
                               ,wct2,sb,sd,se,sf,sg,sh,sm,ss,su,sw,sy,sz &
                               ,rx_lhr,qx_lhr
real, dimension(nzpmax,2)  :: tref,rvsref,rvsrefp
real, dimension(nzpmax,9)  :: sa
real, dimension(nzpmax,12) :: eff

!Arrays for pristine ice - rain collisions
real,dimension(95) :: rdbF94
real,dimension(94) :: radF94,ramF94,avtF94

real, dimension(nzpmax,ncat,ncat) :: rxfer,qrxfer,enxfer
real, dimension(nhcat,maxgrds)    :: dispemb0,dispemb1,ch2

real, dimension(nembc,nembc,npairc) :: coltabc
real, dimension(nembc,nembc,npairr) :: coltabr

real, dimension(nrhhz,nthz)       :: frachz
real, dimension(ndnc,ntc,maxgrds) :: fracc
real, dimension(4)                :: gamm,gamn1
real, dimension(ngam,3)           :: gam
real, dimension(ngam,2)           :: gaminc
real, dimension(2,ngam)           :: gamsip13,gamsip24
real, dimension(ninc,nhcat)       :: rmlttab
real, dimension(ninc,nhcat)       :: enmlttab
real, dimension(ninc,ndns)        :: shedtab
real, dimension(2)                :: sc,sk,sl
real, dimension(8)                :: sj,pcprx,accpx

!Lookup tables for droplet autoconversion and liquid collection
real, dimension(ndcc)        :: r1tabcc,c1tabcc,c2tabcc
real, dimension(ndcc)        :: r2tabdd,c2tabdd,c3tabdd
real, dimension(ndcd,ndcd)   :: r1tabcd,c1tabcd,r2tabcd,c2tabcd
real, dimension(ndccr,nrrcr,ndrcr)  :: r1tabcr,c1tabcr,r2tabcr,c2tabcr

character(len=strl1) :: dustfile

!Lookup table arrays for Binned Riming Scheme
real, dimension(ndccr,nrrcr,ndrcr,4) :: r1tabci,c1tabci,r2tabci,c2tabci &
                                       ,r1rimer,r2rimer

!******Variables Needed for COMPUTING BUDGETS ******************************
!For imbudget>=1
real, dimension(nzpmax) :: xlatheatvap,xlatheatfrz,xnuccldrt,xcld2raint &
,xice2raint,xnucicert,xvapliqt,xvapicet,xevapliqt,xevapicet &
,xmelticet,xrimecldt,xaggregatet,xfreezingt,xmeltingt &
,xrain2icet,xlatheatvapt,xlatheatfrzt

!For imbudget>=2
real, dimension(nzpmax) :: xinuchomrt,xinuccontrt,xinucifnrt,xinuchazrt   &
,xvapcldt,xvapraint,xvapprist,xvapsnowt,xvapaggrt,xvapgraut,xvaphailt     &
,xvapdrizt,xevapcldt,xevapraint,xevapprist,xevapsnowt,xevapaggrt          &
,xevapgraut,xevaphailt,xevapdrizt                                         &
,xmeltprist,xmeltsnowt,xmeltaggrt,xmeltgraut,xmelthailt                   &
,xrimecldsnowt,xrimecldaggrt,xrimecldgraut,xrimecldhailt,xrain2prt        &
,xrain2snt,xrain2agt,xrain2grt,xrain2hat,xaggrselfprist                   &
,xaggrselfsnowt,xaggrprissnowt

!For imbudget>=3
real, dimension(nzpmax) :: xdust1cldrt,xdust2cldrt,xdust1drzrt,xdust2drzrt

!******Variables Needed for BUBBLE SIMULATION ******************************
integer :: ibubble,ibubgrd,ibdxia,ibdxiz,ibdyja,ibdyjz,ibdzk1,ibdzk2
real :: bthp,brtp

!******Variables Needed for CONVERGENCE FORCING ****************************
integer :: iconv,icongr,icicent,icjcent,icvert,ickmax,ickcent
real :: cxrad,cyrad,czrad,cdivmax,ctau,ctmax

!******Variables Needed for CCN nucleation and restore *********************
integer :: iccnlev,ic,rgb
real :: cin_max,ccn_max,gccn_max,dust1_max,dust2_max,saltf_max,saltj_max &
 ,salts_max,enxferratio,rxferratio,ccnmass,ccnnum,rxtemp,cxtemp,fracmass &
 ,cxloss,concen_nuc,aeromass,rg,rhosol,cldrat,epsil,ant,rcm,rmlar,rmsma &
 ,power,scnmass,dcnmass,dinmass,abc1_max,abc2_max
real, dimension(nzpmax) :: nifn

!Tracking total aerosol mass, immersion freezing number in hydrometeors cats
real, dimension(nzpmax,ncat) :: cnmhx,immerhx,snmhx,dnmhx,dinhx
real, dimension(ncat) :: pcpraerox,accpaerox,pcprdustx,accpdustx
real, dimension(nzpmax) :: ifnnucx,total_in

!Number of bins in lognormal aerosol distribution
!If you change itbin, you need to recheck its use in micro
!Current rg thresholds are set for itbin of 100
integer, parameter :: itbin=100
real, dimension(itbin) :: binrad,smass,ccncon,ccnmas

!Solubility fraction (epsilon)
integer,parameter :: maxeps=7
real, dimension(maxeps) :: epsfrac
data epsfrac / 0.05,0.1,0.2,0.4,0.6,0.8,1.0 /

!Median radii (meters) for CCN
integer, parameter :: maxrg=20
real, dimension(maxrg) :: rg_ccn
data rg_ccn / 0.01e-6,0.02e-6,0.04e-6,0.08e-6 &
             ,0.16e-6,0.32e-6,0.48e-6,0.64e-6 &
             ,0.96e-6,1.50e-6,2.00e-6,2.50e-6 &
             ,3.00e-6,3.50e-6,4.00e-6,4.50e-6 &
             ,5.00e-6,5.50e-6,6.00e-6,6.50e-6 /

!Number of aerosol species being used & Ice nuclei arrays
!Make sure you change both if you alter number of species
!Each category can have the soluble component be either
!Ammonium sulfate (NH4-2SO4) or Sodium chloride (NaCl)
! 1 = Sub-micron CCN
! 2 = Super-micron GCCN
! 3 = Small mode mineral dust (soluble coating)
! 4 = Large mode mineral dust (soluble coating)
! 5 = Film mode sea salt
! 6 = Jet mode sea salt
! 7 = Spume mode sea salt
! 8 = Sub-micron radius regenerated mixed aerosols
! 9 = Super-micron radius regenerated mixed aerosols
integer :: acat
integer, parameter :: aerocat=11
real, dimension(nzpmax) :: cifnx
real, dimension(nzpmax,2) :: regenmas
real, dimension(nzpmax,aerocat) :: totifnn,totifnm,aerocon,aeromas
real, dimension(aerocat) :: aero_rhosol,aero_epsilon,aero_medrad &
      ,aero_sigma,aero_rg2rm,aero_rg,aero_vap,aero_ratio
integer, dimension(aerocat) :: iaero_chem,aero_vanthoff

!Minimum aerosol concentration (#/kg) and mass (kg/kg) 
!values for condition statements involving aerosols
real, parameter :: mincon=1.0e-1         &  
                  ,minmas=1.0e-21        &
                  ,maxaero=20000.e6      &
                  ,minmashydro=1.0e-27   &
                  ,minifn=1.0e-14

!Aerosol distribution spectral width (sigma)
!Note: do not change this unless you update cloud nucleation lookup
!tables with parcel model runs using a new sigma. The median radii bins
!in nucleation and pre_nucleation routines (rg, rmsma, rmlar) are also
!specifically set for sigma=1.8. These would need to be updated as well.
data aero_sigma  / 1.80 &       !CCN 
                  ,1.80 &       !GCCN 
                  ,1.80 &       !small mineral dust
                  ,1.80 &       !large mineral dust
                  ,1.80 &       !salt film mode 
                  ,1.80 &       !salt jet mode 
                  ,1.80 &       !salt spume mode
                  ,1.80 &       !absorbing carbon mode-1
                  ,1.80 &       !absorbing carbon mode-2
                  ,1.80 &       !sub-micro regenerated aerosol (mixed)
                  ,1.80 /       !super-micro regenerated aerosol (mixed)
!Set the relationship between median radius and mean mass radius 
!based on aerosol distribution spectral width
!exp(1.5 * (alog(sigma))**2)
data aero_rg2rm  / 1.6791 &     !CCN 
                  ,1.6791 &     !GCCN 
                  ,1.6791 &     !small mineral dust
                  ,1.6791 &     !large mineral dust
                  ,1.6791 &     !salt film mode 
                  ,1.6791 &     !salt jet mode 
                  ,1.6791 &     !salt spume mode
                  ,1.6791 &     !absorbing carbon mode-1
                  ,1.6791 &     !absorbing carbon mode-2
                  ,1.6791 &     !sub-micro regenerated aerosol (mixed)
                  ,1.6791 /     !super-micro regenerated aerosol (mixed)

!*********************************************************************
! R.W. Carver's banded sedimentation scheme, defines band number, and 
! banded variables from data derived from Mitchell (1996)
! Last Updated on 03-30-2005
!---------------------------------------------------------------------
! Notes on banded data
! The idea is simple, looking at Mitchell's power-law data, on a 
! log-log scale, we see that vt(diam) is a curve that can be 
! well approximated by a set of power-laws.  The starting point
! for the valid range of each law is given by rc_diam.  For diams <
! the first band cutoff, we just use the data for the first band.  Not
! much of a problem since the first cutoffs are near dmb0
!---------------------------------------------------------------------
integer, parameter :: nband=3,nembfall=40,ndensrtgt=40,maxkfall=4
real, dimension(nband,nhcat) :: bcfvt,bpwvt,bdenpow,bdiam

!Band Data by Habit
               !Small     !Medium   !Large
data bcfvt /   1.26E7,    1.26E7,   234101.5, & ! 1 cloud 
               2032.,     143.9,    143.9,    & ! 2 rain
               3207543.,  72523.,   1538.,    & ! 3 pristine columns
               871.,      27.7,     5.00,     & ! 4 snow columns
               427.3,     16.1,     3.26,     & ! 5 aggregates
               22133.2,   332.4,    34.9,     & ! 6 graupel
               256914.,   2183.,    152.1,    & ! 7 hail
               106861.,   106861.,  20801.,   & ! 8 pristine hex
               74256.,    754,      56.4,     & ! 9 pristine dendrites (fix)
               37880.,    37880.,   1617.9,   & ! 10 pristine needles
               183098.,   183098.,  6239.,    & ! 11 pristine rosette
               968.1,     30.08,    5.33,     & ! 12 snow hex
               753.9,     56.4,     3.39,     & ! 13 snow dendrites
               1617.,     44.6,     7.26,     & ! 14 snow needles
               6239.,     125.7,    16.3,     & ! 15 snow rosettes
               1.26E7,    1.26E7,   234101.5 /  ! 16 drizzle droplets

data bpwvt /   1.91,      1.91,     1.493,    &
               0.914,     0.497,    0.497,    &
               1.82,      1.42,     1.00,     &
               0.933,     0.484,    0.161,    &
               0.844,     0.416,    0.108,    &
               1.33,      0.786,    0.397,    &
               1.493,     0.914,    0.497,    &
               1.522,     1.522,    1.377,    &
               1.492,     0.978,    0.695,    &
               1.31,      1.31,     0.983,    &
               1.61,      1.61,     1.24,     &
               1.04,      0.563,    0.222,    &
               0.979,     0.695,    0.302,    &
               0.982,     0.522,    0.191,    &
               1.24,      0.716,    0.342,    &
               1.91,      1.91,     1.493 /

data bdenpow / 0.0,       0.0,      0.0,      &
               0.169,     0.362,    0.501,    &
               0.0,       0.0,      0.0,      &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.0,       0.0,      0.0,      &
               0.0,       0.0,      0.0,      &
               0.0,       0.0,      0.0,      &
               0.0,       0.0,      0.0,      &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.169,     0.362,    0.501,    &
               0.0,       0.0,      0.0 /

data bdiam /   15E-6,     25E-6,    71E-6,    &
               275E-6,    1768E-6,  11768E-6, &
               15E-6,     86E-6,    101E-6,   &
               301E-6,    462E-6,   5091E-6,  &
               76E-6,     471E-6,   5820E-6,  &
               99E-6,     422E-6,   3099E-6,  &
               69E-6,     265E-6,   1704E-6,  &
               15E-6,     25E-6,    101E-6,   &
               15E-6,     101E-6,   108E-6,   &
               15E-6,     25E-6,    75E-6,    &
               15E-6,     25E-6,    120E-6,   &
               123E-6,    646E-6,   6314E-6,  &
               101E-6,    108E-6,   789E-6,   &
               75E-6,     410E-6,   4254E-6,  &
               120E-6,    541E-6,   4311E-6,  &
               15E-6,     25E-6,    71E-6 /
!---------------------------------------------------------------
!lcat Species Category
!  1  cloud
!  2  rain
!  3  pristine
!  4  snow
!  5  aggregate
!  6  graupel
!  7  hail
!  8  large cloud droplets
!---------------------------------------------------------------
!---------------------------------------------------------------
!lhcat Habit Category
! 1 cloud
! 2 rain
! 3 pristine columns
! 4 snow columns
! 5 aggregates
! 6 graupel
! 7 hail
! 8 pristine hex
! 9 pristine dendrites
! 10 pristine needles
! 11 pristine rosette
! 12 snow hex
! 13 snow dendrites
! 14 snow needles
! 15 snow rosettes 
! 16 large cloud droplets
!---------------------------------------------------------------

END MODULE micphys
