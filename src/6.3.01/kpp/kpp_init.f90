!##############################################################################
Subroutine mckpp_initialize_fields ()

use mem_kpp, only:kpp_3d_fields,kpp_const_fields
use kpp_parameters, only:kppNZ,kppNZP1,nkppz
use mem_grid, only:runtype,initial,print_msg

implicit none

  kppNZ   = nkppz-1  !number of layers from KPP used in KPP 
  kppNZP1 = kppNZ+1  !number of grid point levels used in KPP

  !Call routine to copy constants and logicals needed for ocean
  !physics into the kpp_const_fields derived type.
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_CONSTANTS'
  CALL mckpp_initialize_constants (kpp_const_fields) 

  !Call routine to set up lat/lon kpp grid and ocean depth
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_BATHYMETRY'
  CALL mckpp_initialize_bathymetry ()

  !Initialize the vertical grid
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_GEOGRAPHY'
  CALL mckpp_initialize_geography ()

  !Initialize water type for optical properties of seawater
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_OPTICS'
  CALL mckpp_initialize_optics ()

  !Initialize ocean profiles (most of this not done for HISTORY runs)
  !Here setup routine to set: (must be interpolated to KPP levels)
  !1. "ocnt_clim" (climo temperature profiles)
  !2. "sal_clim" (climo salinity profiles)
  !3. "bottomt" (bottom level ocean temperature)
  !4. U and V ocean currents
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_OCEAN_PROFILES'
  CALL mckpp_initialize_ocean_profiles ()

  !Generate lookup tables for wmt,wst used to get turbulent velocity
  !scales at normalized depth (sigma)(d/hbl)
  if(print_msg) print*,'Calling MCKPP_PHYSICS_LOOKUP'
  CALL mckpp_physics_lookup ()

  !Initialize ocean model (most of this not done for HISTORY runs)
  !and get initial HMIX and fluxes with zero initial Atmos forcing.
  if(print_msg) print*,'Calling MCKPP_INITIALIZE_OCEAN_MODEL'
  CALL mckpp_initialize_ocean_model ()

return
END SUBROUTINE mckpp_initialize_fields

!##############################################################################
Subroutine mckpp_initialize_constants (kpp_const_fields)

use mem_kpp, only:kpp_const_type
use mem_grid, only:print_msg,dtlt
use kpp_parameters
  
implicit none

TYPE(kpp_const_type),intent(inout) :: kpp_const_fields

  !Initialize constants
  kpp_const_fields%dto=FRQKPP

  !Salinity relaxation (1/sec)
  if(relax_sal < 1.0e-12) then
   kpp_const_fields%relax_sal  = 1.0e-12
  else
   kpp_const_fields%relax_sal  = 1. / (relax_sal*86400.)
  endif

  !Ocean temperature profiles relaxation (1/sec)
  if(relax_ocnT < 1.0e-12) then
   kpp_const_fields%relax_ocnT  = 1.0e-12
  else
   kpp_const_fields%relax_ocnT  = 1. / (relax_ocnT*86400.)
  endif

  !SST relaxation (1/sec)
  if(relax_sst < 1.0e-12) then
   kpp_const_fields%relax_sst  = 1.0e-12
  else
   kpp_const_fields%relax_sst  = 1. / (relax_sst*86400.)
  endif

  !Relaxation timescale for damping currents (1/sec)
  kpp_const_fields%dt_uvdamp  = 1. / (360. * 86400.)
  !Maximum ocean depth (meters)
  kpp_const_fields%DMAX       = DMAXKPP
  !Exponential for stretched grid
  kpp_const_fields%dscale     = DSCALEKPP

  !convergence tolerance for hmix(new)-hmix(old) iteration in ocnstep
  !frac of layer thickness hm(kmix)
  kpp_const_fields%hmixtolfrac=0.1
  !salinity of ice
  kpp_const_fields%sice=4.0
  !deepest model layer to check for isothermal (90% of levels)
  kpp_const_fields%iso_bot=nint(nkppz - nkppz*0.1)
  !threshold for resetting isothermal layer
  kpp_const_fields%iso_thresh=0.002
  
  !Flag for running boundary layer mixing scheme
  kpp_const_fields%LKPP=.TRUE.
  !Flag for Ri mixing below the diagnosed mixing depth
  kpp_const_fields%LRI=.TRUE.
  !Flag for double diffusion calculations (this is slow)
  kpp_const_fields%LDD=.FALSE.
  !Flag for using 2D ocean optical data
  kpp_const_fields%L_JERLOV=.TRUE.
  !Vertical grid flags
  kpp_const_fields%L_STRETCHGRID=.TRUE.
  !Flag to prevent unrealistic ocean freezing
  kpp_const_fields%L_NO_FREEZE=.TRUE.
  !Flag for Prevent isothermal layers
  kpp_const_fields%L_NO_ISOTHERM=.TRUE.
  !Flags to initialize Bottom-temp
  kpp_const_fields%L_BOTTOM_TEMP=.TRUE.
  !Damp KPP currents at all levels to zero using non-linear relaxation.
  !Recommended for climate-length simulations
  kpp_const_fields%L_DAMP_CURR=.TRUE.

  !Make sure an ocean depth is given
  IF (kpp_const_fields%DMAX .LE. 0.0) THEN 
   if(print_msg)print*,'KPP:You must specify a depth for the domain'
   STOP
  ENDIF

  !Make sure stretched grid is given a scale
  IF ((kpp_const_fields%L_STRETCHGRID) .AND. &
      (kpp_const_fields%dscale .EQ. 0.0)) THEN
   if(print_msg)print*,'KPP:You cannot have dscale=0 for stretched grids'
   STOP
  ENDIF

  !Do not do SST and OCNT corrections together
  IF (kpp_const_fields%relax_sst  > 1.e-10 .and. &
      kpp_const_fields%relax_ocnT > 1.e-10)THEN         
   if(print_msg)print*,'KPP:Set RELAX_SST>0 or RELAX_OCNT>0; Not both'
   STOP
  ENDIF

  !Isothermal detection routine requires 3D ocean temperature and salinity
  IF (kpp_const_fields%L_NO_ISOTHERM .AND.       &
     (kpp_const_fields%relax_sal  <= 1.e-10 .OR. &
      kpp_const_fields%relax_ocnT <= 1.e-10))THEN         
   if(print_msg)print*,'KPP:L_NO_ISOTHERM requires RELAX_SAL and RELAX_OCNT'
   STOP
  ENDIF

return
END SUBROUTINE mckpp_initialize_constants

!##############################################################################
Subroutine mckpp_initialize_bathymetry ()

use mem_kpp, only:kpp_3d_fields
use mem_grid, only:grid_g,ngrid
use node_mod, only:ia,iz,ja,jz

implicit none

integer :: i,j

  ! Initialize ocean depth / bathymetry.
  ! Ocean depth values should be NEGATIVE!
  do j=ja,jz
   do i=ia,iz
     kpp_3d_fields(ngrid)%ocdepth(i,j)=-1000. !from bathymetry !Saleeby(2018)
   enddo
  enddo

return
END SUBROUTINE mckpp_initialize_bathymetry

!##############################################################################
Subroutine mckpp_initialize_geography ()

use mem_kpp, only:kpp_const_fields,kpp_3d_fields
use mem_grid, only:grid_g,ngrid,print_msg
use node_mod, only:ia,iz,ja,jz
use kpp_parameters, only:kppNZ,kppNZP1

implicit none

  ! Local Variables
  REAL sumh,hsum,dfac,sk,twopi
  INTEGER i,j

  twopi = 8.*atan(1.)

  IF (kpp_const_fields%L_STRETCHGRID) THEN
     sumh = 0.0
     dfac = 1.0 - exp(-kpp_const_fields%dscale)
     DO i = 1,kppNZ
        sk = - (float(i)-0.5)/float(kppNZ)
        kpp_const_fields%hm(i) = kpp_const_fields%DMAX*dfac/&
             float(kppNZ)/kpp_const_fields%dscale / ( 1.0 + sk*dfac )
        sumh = sumh + kpp_const_fields%hm(i)
     ENDDO
  ENDIF

  ! layer thickness h, layer grids zgrid, interface depths d
  hsum = 0.0
  kpp_const_fields%dm(0) = 0.0
  DO i=1,kppNZ
     if(kpp_const_fields%L_STRETCHGRID) then
        kpp_const_fields%hm(i) = kpp_const_fields%hm(i) &
           * kpp_const_fields%DMAX / sumh 
     else   
        kpp_const_fields%hm(i) = kpp_const_fields%DMAX / real(kppNZ) 
     endif
     kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * kpp_const_fields%hm(i) )
     hsum = hsum + kpp_const_fields%hm(i)
     kpp_const_fields%dm(i) = hsum
  ENDDO

  kpp_const_fields%dm(kppNZP1) = -9999.
  kpp_const_fields%hm(kppNZP1) = 1.e-10 
  kpp_const_fields%zm(kppNZP1) = -kpp_const_fields%DMAX

  ! Calculate Coriolis parameter
  ! Enforce minimum value of Coriolis parameter equal to 2.5 degrees latitude
  do j=ja,jz
   do i=ia,iz
    IF (ABS(grid_g(ngrid)%glat(i,j)) .lt. 2.5) THEN
      kpp_3d_fields(ngrid)%f(i,j) = 2. * (twopi/86164.) * &
           sin(2.5*twopi/360.) &
         * SIGN(1.,grid_g(ngrid)%glat(i,j))
    ELSE
      kpp_3d_fields(ngrid)%f(i,j) = 2. * (twopi/86164.) * &
           sin(grid_g(ngrid)%glat(i,j)*twopi/360.)
    ENDIF
   ENDDO
  ENDDO

  !Print some header info to the RAMS standard output
  if(print_msg) then
   do i=1,kppNZP1
    if(i==1) then
     print*,''
     print*,'KPP Ocean Mixed Layer Model Layer Info: 0th level = sea-surface'
     print*,'  KPP-level     LayerInterface     GridCenter     LayerThickness'
    endif
    print*,i,'   ',kpp_const_fields%dm(i),'  ',kpp_const_fields%zm(i) &
            ,' ',kpp_const_fields%hm(i)
   enddo
  endif
  
return
END SUBROUTINE mckpp_initialize_geography

!##############################################################################
Subroutine mckpp_initialize_optics ()

use mem_kpp, only:kpp_3d_fields,kpp_const_fields
use node_mod, only:ia,iz,ja,jz
use mem_grid, only:ngrid,print_msg
use kpp_parameters, only:kppNZ,kppNZP1

implicit none

  integer nwtype,i,j
  parameter(nwtype=5) ! max number of different water types 

  ! Input
  real hbf                 ! scale factor to apply to depth array
  
  real  rfac(nwtype),a1(nwtype),a2(nwtype)
  real rmin,r1,r2
  integer l
  
  ! jerlov water type :  I       IA      IB      II      III
  !            jwtype    1       2       3       4       5
  !
  data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
  data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
  data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
  data rmin         / -80. /
  data hbf          /  1.0    /

  !Saleeby(2018): Could add input of Jerlov ocean optical regime here.
  !Might have to assign data over ia,iz and ja,jz parallel loops.
  !If data not available use default values of 3.
  !Jerlov: 1=I, 2=IA, 3=IB, 4=II, 5=III
  IF (kpp_const_fields%L_JERLOV) THEN
    do j=ja,jz
     do i=ia,iz
        kpp_3d_fields(ngrid)%jerlov(i,j)=1 !read jerlov; Saleeby(2018)
     enddo
    enddo
  ELSE
    do j=ja,jz
     do i=ia,iz
        kpp_3d_fields(ngrid)%jerlov(i,j)=1
     enddo
    enddo
  ENDIF

  !compute fraction of solar short-wave flux penetrating to specified
  !depth (times hbf) due to exponential decay in  Jerlov water type
  !reference : two band solar absorption model of simpson and 
  !paulson (1977)
  do j=ja,jz
   do i=ia,iz
    !set shortwave fraction on ZM depths (zm values are negative)
    do l=1,kppNZP1
     r1=MAX(kpp_const_fields%zm(l)*hbf/a1(nint(kpp_3d_fields(ngrid)%jerlov(i,j))),rmin)
     r2=MAX(kpp_const_fields%zm(l)*hbf/a2(nint(kpp_3d_fields(ngrid)%jerlov(i,j))),rmin)
     kpp_3d_fields(ngrid)%swfrac(l,i,j)=rfac(nint(kpp_3d_fields(ngrid)%jerlov(i,j))) &
          * exp(r1) + (1.-rfac(nint(kpp_3d_fields(ngrid)%jerlov(i,j)))) * exp(r2)
    enddo
    !set shortwave fraction on DM depths (dm values are positive)
    do l=0,kppNZ
     r1=MAX(-kpp_const_fields%dm(l)/a1(nint(kpp_3d_fields(ngrid)%jerlov(i,j))),rmin)
     r2=MAX(-kpp_const_fields%dm(l)/a2(nint(kpp_3d_fields(ngrid)%jerlov(i,j))),rmin)
     kpp_3d_fields(ngrid)%swdk_opt(l+1,i,j)=rfac(nint(kpp_3d_fields(ngrid)%jerlov(i,j))) &
          * exp(r1) + (1.-rfac(nint(kpp_3d_fields(ngrid)%jerlov(i,j)))) * exp(r2) 
    enddo
   enddo
  enddo

  !Print some header info to the RAMS standard output
  if(print_msg) then
   do l=1,kppNZP1
    if(l==1) then
     print*,''
     print*,'KPP Ocean Mixed Layer Model Layer Info: 0th level = sea-surface'
     print*,'Profile of penetrating Shortwave fraction at a sampled point'
     print*,'  KPP-level      swfrac(zm)    swdk_opt(dm)'
    endif
    print*,l,'   ',kpp_3d_fields(ngrid)%swfrac(l,ia,ja) &
                  ,kpp_3d_fields(ngrid)%swdk_opt(l,ia,ja)
   enddo
  endif

return
END SUBROUTINE mckpp_initialize_optics

!##############################################################################
Subroutine mckpp_initialize_ocean_profiles ()

use mem_kpp, only:kpp_3d_fields,kpp_const_fields
use node_mod, only:ia,iz,ja,jz
use mem_grid, only:ngrid,runtype,initial,print_msg,time
use kpp_parameters, only:kppNZ,kppNZP1,nkppz,ikpp
use mem_leaf
use io_params, only:ssttime1,ssttime2,iupdsst

implicit none

integer :: i,j,k
real :: sstinit,timefac_sst
real, dimension(36)::otz,ot,os,ou,ov,otratio
real, dimension(kppNZP1)::oceant,oceans,oceanu,oceanv

  !Only call this for RUNTYPE=INITIAL and INITIAL/=3 simulation starting
  !at time=0; not for RAMS history restart (ie. RUNTYPE=HISTORY) or 
  !RAMS history initialization (ie. RUNTYPE=INITIAL, INITIAL=3), both of
  !which interpolate or use previous data at model start.

  !History restart will read in and set
  ! U, U_init, X, Sref, Tref which are initialized below

 !ocean depth (meters) array for idealized ocean temp profile
 !Emily Riley Changed 4/14/19 from Steve's profile to ERD's GLORY 
 !monthly reanalysis profile for PISTON Luzon simulation. Area averaged
 !over ocean area similar to Riley-Delaripa et al.(2019) domain.
 !The u/v ocean currents are in m/s. I downloaded the 2016 July and Aug monthly 
 !ocean reanalysis and took the average of July and August. For example, if the 
 !July top level temperature was 30C and the August top level temperature was 
 !31C, the top level temperature value I would use is 30.5C. I averaged in 
 !x and y for the same area as the model domain, where I made sure the model 
 !points were included (i.e., the reanalysis points had to be equal to or 
 !inside the model domain; not outside the model domain). 
 !Ocean reanalysis x/y bounds: 116.333°E - 125.667°E & 10.5°N - 19.5°N
 !Domain x/y bounds: 116.354°E - 125.646°E & 10.512°N - 19.488°N
 otz(1)=0.49402499
 otz(2)=1.5413750
 otz(3)=2.6456690
 otz(4)=3.8194950
 otz(5)=5.0782242
 otz(6)=6.4406142
 otz(7)=7.9295602
 otz(8)=9.5729971
 otz(9)=11.405000
 otz(10)=13.467140
 otz(11)=15.810070
 otz(12)=18.495560
 otz(13)=21.598820
 otz(14)=25.211411
 otz(15)=29.444731
 otz(16)=34.434151
 otz(17)=40.344051
 otz(18)=47.373692
 otz(19)=55.764290
 otz(20)=65.807266
 otz(21)=77.853851
 otz(22)=92.326073
 otz(23)=109.72930
 otz(24)=130.66600
 otz(25)=155.85069
 otz(26)=186.12560
 otz(27)=222.47520
 otz(28)=266.04031
 otz(29)=318.12741
 otz(30)=380.21301
 otz(31)=453.93771
 otz(32)=541.08893
 otz(33)=643.56677
 otz(34)=763.33313
 otz(35)=902.33929
 otz(36)=1062.4399
 ot(1)=30.066090
 ot(2)=30.024014
 ot(3)=29.997210
 ot(4)=29.985439
 ot(5)=29.979153
 ot(6)=29.975258
 ot(7)=29.970398
 ot(8)=29.966806
 ot(9)=29.962917
 ot(10)=29.958921
 ot(11)=29.953201
 ot(12)=29.944405
 ot(13)=29.930378
 ot(14)=29.902710
 ot(15)=29.838472
 ot(16)=29.683643
 ot(17)=29.347584
 ot(18)=28.751217
 ot(19)=27.928467
 ot(20)=26.898888
 ot(21)=25.651855
 ot(22)=24.047022
 ot(23)=22.163759
 ot(24)=20.229223
 ot(25)=18.464861
 ot(26)=16.844282
 ot(27)=15.242334
 ot(28)=13.742815
 ot(29)=12.298230
 ot(30)=10.732548
 ot(31)=9.2441196
 ot(32)=7.8945026
 ot(33)=6.8418884
 ot(34)=5.8521953
 ot(35)=4.8802743
 ot(36)=4.0321999
 os(1)=34.152065
 os(2)=34.154835
 os(3)=34.157471
 os(4)=34.159668
 os(5)=34.161785
 os(6)=34.164009
 os(7)=34.167671
 os(8)=34.171135
 os(9)=34.175179
 os(10)=34.179420
 os(11)=34.184959
 os(12)=34.192741
 os(13)=34.203957
 os(14)=34.221706
 os(15)=34.251175
 os(16)=34.296738
 os(17)=34.360058
 os(18)=34.431709
 os(19)=34.499954
 os(20)=34.565399
 os(21)=34.634483
 os(22)=34.702820
 os(23)=34.748833
 os(24)=34.754417
 os(25)=34.736214
 os(26)=34.696487
 os(27)=34.625751
 os(28)=34.543259
 os(29)=34.471127
 os(30)=34.412819
 os(31)=34.381081
 os(32)=34.382336
 os(33)=34.414345
 os(34)=34.454685
 os(35)=34.494591
 os(36)=34.533295
 ou(1)=0.019100454
 ou(2)=0.016321890
 ou(3)=0.014060168
 ou(4)=0.011857715
 ou(5)=0.0096536130
 ou(6)=0.0074455901
 ou(7)=0.0050409529
 ou(8)=0.0027234163
 ou(9)=0.00023679016
 ou(10)=-0.0023469469
 ou(11)=-0.0052551287
 ou(12)=-0.0085307565
 ou(13)=-0.012290464
 ou(14)=-0.016672498
 ou(15)=-0.021679830
 ou(16)=-0.027013715
 ou(17)=-0.030701458
 ou(18)=-0.031918801
 ou(19)=-0.033096995
 ou(20)=-0.034624822
 ou(21)=-0.035509244
 ou(22)=-0.034193810
 ou(23)=-0.030427990
 ou(24)=-0.026155602
 ou(25)=-0.024032710
 ou(26)=-0.021403963
 ou(27)=-0.018873321
 ou(28)=-0.016758438
 ou(29)=-0.013506316
 ou(30)=-0.0093148304
 ou(31)=-0.0050458922
 ou(32)=-0.0015908780
 ou(33)=-0.00048178932
 ou(34)=0.00017648222
 ou(35)=0.00033577811
 ou(36)=0.00055340666
 ov(1)=0.033420175
 ov(2)=0.029788986
 ov(3)=0.026523601
 ov(4)=0.023757031
 ov(5)=0.021405082
 ov(6)=0.019356415
 ov(7)=0.017362187
 ov(8)=0.015635069
 ov(9)=0.014054944
 ov(10)=0.012625669
 ov(11)=0.011338353
 ov(12)=0.010397403
 ov(13)=0.0099184299
 ov(14)=0.010435903
 ov(15)=0.012561996
 ov(16)=0.017077940
 ov(17)=0.022756575
 ov(18)=0.026075317
 ov(19)=0.025187567
 ov(20)=0.021348305
 ov(21)=0.015851133
 ov(22)=0.0082361754
 ov(23)=0.0013649256
 ov(24)=-0.0036659446
 ov(25)=-0.0073529691
 ov(26)=-0.0092373220
 ov(27)=-0.0094265882
 ov(28)=-0.0094937831
 ov(29)=-0.011583607
 ov(30)=-0.014184179
 ov(31)=-0.016872544
 ov(32)=-0.014716862
 ov(33)=-0.010820009
 ov(34)=-0.0066377516
 ov(35)=-0.0038285339
 ov(36)=-0.0023947731

 !Compute ratio of ocean depth T to SST (top level T)
 !The idea is to use the above default ocean T profile (ot) for any given
 !SST. Keep the shape of the profile but adjust according to SST with
 !the caveat that we will not initialize with T < -1.0C.
 otratio=ot/ot(1)
 ot=otratio
 !do k=1,36
 ! write(*,'(10X,a8,i2,a2,f12.8)') 'otratio(',k,')=',otratio(k)
 !enddo

 !Interpolate ocean profiles from data to KPP depth levels
 CALL htint (36,ot,otz,kppNZP1,oceant,(-1.*kpp_const_fields%zm))
 CALL htint (36,os,otz,kppNZP1,oceans,(-1.*kpp_const_fields%zm))
 CALL htint (36,ou,otz,kppNZP1,oceanu,(-1.*kpp_const_fields%zm))
 CALL htint (36,ov,otz,kppNZP1,oceanv,(-1.*kpp_const_fields%zm))

 !Only run this for true INITIAL "cold start" of simulations.
 !Do not even do this for HISTORY INITIALIZATION since interpolation is done.
 !Do not do this for HISTORY restart.
 IF(trim(runtype) == 'INITIAL' .and. INITIAL /= 3) THEN     

  !Time interpolation factor for updating SST
  if (iupdsst == 0) then
    timefac_sst = 0.
  else
    timefac_sst = (time-ssttime1(ngrid)) / (ssttime2(ngrid)-ssttime1(ngrid))
  endif

  do j=ja,jz
   do i=ia,iz
    !Determine initial SST to apply to ocean column
    sstinit = (leaf_g(ngrid)%seatp(i,j) &
        + (leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j)) &
        * timefac_sst - 273.15)
    !Loop over ocean depth levels
    do k=1,kppNZP1
     !###### Initialized zonal ocean current velocity
     !input from RAMS (Saleeby)
     kpp_3d_fields(ngrid)%U(k,i,j) = oceanu(k)
     kpp_3d_fields(ngrid)%U_init(k,i,j)=kpp_3d_fields(ngrid)%U(k,i,j)

     !###### Initialized meridional ocean current velocity  
     !input from RAMS (Saleeby)
     kpp_3d_fields(ngrid)%V(k,i,j) = oceanv(k)
     kpp_3d_fields(ngrid)%V_init(k,i,j)=kpp_3d_fields(ngrid)%V(k,i,j)

     !###### Initialized ocean temperature profile
     !KPP requires temperatures in CELSIUS.  If initial conditions
     !are in Kelvin, subtract 273.15
     !input from RAMS (Saleeby)
     kpp_3d_fields(ngrid)%ocnT_clim(k,i,j) = max(-1.0,oceant(k)*sstinit)
     kpp_3d_fields(ngrid)%X_t(k,i,j) = max(-1.0,oceant(k)*sstinit)

     !###### Initialize salinity profile
     !input from RAMS (Saleeby)
     kpp_3d_fields(ngrid)%sal_clim(k,i,j)   = oceans(k)
     kpp_3d_fields(ngrid)%X_s(k,i,j) = oceans(k)
    enddo

    !Calculate reference salinity      
    kpp_3d_fields(ngrid)%Sref(i,j)= &
     (kpp_3d_fields(ngrid)%X_s(1,i,j)+kpp_3d_fields(ngrid)%X_s(kppNZP1,i,j))/2.

    !Remove reference salinity from X and sal_clim
    do k=1,kppNZP1
     kpp_3d_fields(ngrid)%X_s(k,i,j)= &
      kpp_3d_fields(ngrid)%X_s(k,i,j)-kpp_3d_fields(ngrid)%Sref(i,j)
     kpp_3d_fields(ngrid)%sal_clim(k,i,j)= &
      kpp_3d_fields(ngrid)%sal_clim(k,i,j)-kpp_3d_fields(ngrid)%Sref(i,j)
    enddo

    !Initial bottom temperature
    kpp_3d_fields(ngrid)%bottomt(i,j) = kpp_3d_fields(ngrid)%X_t(kppNZP1,i,j)
   enddo
  enddo

 ENDIF !end U,X,init initialization if NOT a HISTORY run

 !Print some header info to the RAMS standard output
 if(print_msg) then
  do k=1,kppNZP1
   if(k==1) then
    print*,''
    print*,'KPP Ocean Mixed Layer Model Layer Info: 0th level = sea-surface'
    print*,'Profile of ocean temperature (C) and salinity at a sampled point'
    print*,'  KPP-level      ocean-Depth      ocean-T        ocean-S'
   endif
   print*,k,'   ',kpp_const_fields%zm(k),kpp_3d_fields(ngrid)%X_t(k,ia,ja) &
      ,kpp_3d_fields(ngrid)%X_s(k,ia,ja)+kpp_3d_fields(ngrid)%Sref(ia,ja)
  enddo
  do k=1,kppNZP1
   if(k==1) then
    print*,''
    print*,'KPP Ocean Mixed Layer Model Layer Info: 0th level = sea-surface'
    print*,'Profile of ocean currents at a sampled point'
    print*,'  KPP-level      ocean-Depth    U-current      V-current'
   endif
   print*,k,'   ',kpp_const_fields%zm(k),kpp_3d_fields(ngrid)%U(k,ia,ja) &
      ,kpp_3d_fields(ngrid)%V(k,ia,ja)
  enddo
 endif

return
END SUBROUTINE mckpp_initialize_ocean_profiles

!##############################################################################
Subroutine mckpp_physics_lookup ()

use mem_kpp, only:kpp_const_fields

implicit none

  !von karman constant
  real, parameter :: vonk=0.4

  real zmin,zmax,umin,umax,usta,zeta,zehat,epsln,&
       am,cm,c1,c2,zetam,as,cs,c3,zetas,cstar,deltau,deltaz
  integer i,j,ni,nj
  parameter ( ni = 890,&     ! number of values for zehat
       nj = 48)             ! number of values for ustar      
  data epsln             /   1.e-20 /
  data c1                /   5.0    /
  data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
  data umin,umax  /  0.   , .04   / ! m/s
  data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
  data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
  data cstar             /    5.    /
  
  deltaz = (zmax-zmin)/(ni+1) 
  deltau = (umax-umin)/(nj+1)
  
  DO i=0,ni+1
     zehat = deltaz*(i) + zmin
     DO j=0,nj+1
      usta = deltau*(j) + umin
      zeta = zehat/(usta**3+epsln)

      if(zehat.ge.0.) then
         kpp_const_fields%wmt(i,j) = vonk*usta/(1.+c1*zeta)
         kpp_const_fields%wst(i,j) = kpp_const_fields%wmt(i,j)
      else
        if(zeta.gt.zetam) then
         kpp_const_fields%wmt(i,j) = vonk * usta * (1.-c2*zeta)**(1./4.)
        else
         kpp_const_fields%wmt(i,j) = vonk * (am*usta**3 - cm*zehat)**(1./3.)
        endif
        if(zeta.gt.zetas) then
         kpp_const_fields%wst(i,j) = vonk * usta * (1.-c3*zeta)**(1./2.)
        else
         kpp_const_fields%wst(i,j) = vonk * (as*usta**3 - cs*zehat)**(1./3.)
        endif
      endif
     ENDDO
  ENDDO
      
return
END SUBROUTINE mckpp_physics_lookup

!##############################################################################
Subroutine mckpp_initialize_ocean_model ()

use mem_kpp, only:kpp_3d_fields,kpp_const_fields,kpp_1d_fields
use mem_leaf, only:leaf_g
use node_mod, only:ia,iz,ja,jz
use mem_grid, only:ngrid,runtype,initial,print_msg
use kpp_parameters, only:kppNZ,kppNZP1,kppNVEL,kppNSCLR

implicit none

  ! Initialize ocean model:
  ! Set coefficients for tridiagonal matrix solver.
  ! Compute hmix and diffusivity profiles for initial profile.
  ! Prepare for first time step.
  
  ! Local
  real dzb(kppNZ)              ! diff. between grid-levels below z(j)
  integer k,kmix0,n,l,i,j
  real hmix0,deltaz

  !gravity constant
  real, parameter :: grav=9.816

  ! Compute factors for coefficients of tridiagonal matrix elements.
  ! tri(0     ,1) : dt/h(1) factor for rhs flux
  ! tri(k=1:kppNZ,0) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
  ! tri(k=1:kppNZ,1) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}

 !Set layer differences in meters
 DO k=1,kppNZ
    dzb(k) = kpp_const_fields%zm(k) - kpp_const_fields%zm(k+1)
 ENDDO

 !Set matrix values
 kpp_const_fields%tri(0,1)=kpp_const_fields%dto/kpp_const_fields%hm(1)
 kpp_const_fields%tri(1,1)=kpp_const_fields%dto/kpp_const_fields%hm(1)/dzb(1)
 DO k=2,kppNZ
  kpp_const_fields%tri(k,1)=kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k)
  kpp_const_fields%tri(k,0)=kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k-1)
 ENDDO

 !Initialize if NOT doing history restart or history initialization
 ! with RAMS (Saleeby 2018)
 IF(trim(runtype) == 'INITIAL' .and. INITIAL /= 3) THEN

  !Determine hmix for initial profile:
  do j=ja,jz
   do i=ia,iz
     !Proceed if RAMS LEAF3 has determine this to be an ocean point.
     !We allow coastal grid cells with partial ocean.
     IF (nint(leaf_g(ngrid)%leaf_class(i,j,1)) == 0. .and. &
              leaf_g(ngrid)%patch_area(i,j,1) > 0.50) THEN !Saleeby(2018)

        CALL mckpp_fields_3dto1d  (kpp_3d_fields(ngrid),kpp_1d_fields,i,j)

        CALL MCKPP_PHYSICS_VERTICALMIXING (kpp_1d_fields,kpp_const_fields &
                                          ,hmix0,kmix0,i,j)

        kpp_1d_fields%hmix = hmix0
        kpp_1d_fields%kmix = kmix0
        kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)

        !Evaluate initial fluxes (to write to output data file)
        DO k=1,kppNZ
         deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
         DO n=1,kppNSCLR
           kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*&
                ((kpp_1d_fields%X(k,n)-kpp_1d_fields%X(k+1,n))/deltaz-&
                kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
         ENDDO
         IF(kpp_const_fields%LDD) then 
           kpp_1d_fields%wX(k,1)=-kpp_1d_fields%dift(k)*&
              ((kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1)) &
              /deltaz-kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,1))
         ENDIF
         kpp_1d_fields%wB(k)= grav &
           * (kpp_1d_fields%talpha(k)*kpp_1d_fields%wX(k,1) - &
              kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
         DO  n=1,kppNVEL
            kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*&
                 (kpp_1d_fields%U(k,n)-kpp_1d_fields%U(k+1,n))/deltaz
         ENDDO
        ENDDO

        ! Prepare for first time step
           
        ! indices for extrapolation
        kpp_1d_fields%old = 0
        kpp_1d_fields%new = 1               
        ! initialize array for extrapolating Us,Xs
        DO k=1,kppNZP1
           DO l=1,kppNVEL
              kpp_1d_fields%Us(k,l,0)=kpp_1d_fields%U(k,l)
              kpp_1d_fields%Us(k,l,1)=kpp_1d_fields%U(k,l)
           ENDDO
           DO l=1,kppNSCLR
              kpp_1d_fields%Xs(k,l,0)=kpp_1d_fields%X(k,l)
              kpp_1d_fields%Xs(k,l,1)=kpp_1d_fields%X(k,l)
           ENDDO
        ENDDO

        CALL mckpp_fields_1dto3d (kpp_1d_fields,kpp_3d_fields(ngrid),i,j)

     ENDIF
   ENDDO
  ENDDO

 ENDIF

return
END SUBROUTINE mckpp_initialize_ocean_model

