!##############################################################################
! KAIN-FRITSCH CUMULUS PARAMETERIZATION SCHEME
! RAMS VERSION 4.3.0
! KFDRIVER.F90
! Updated July 2002
!
! INPUT VARIABLES:
!
! mzp: k-dim of array - 0-D
! mxp: i-dim of array - 0-D
! myp: j-dim of array - 0-D
! ia: beginning i-index - 0-D
! iz: ending i-index - 0-D
! ja: beginning j-index - 0-D
! jz: ending j-index - 0-D
! imphys: microphysics flag from RAMSIN
! icloud: cloud water flag from RAMSIN
! irain: rain water flag from RAMSIN
! ipris: pristine ice flag from RAMSIN
! isnow: snow flag from RAMSIN
! idiffk: turbulence parameterization flag
! up: u-wind (m/s) - 3-D
! vp: v-wind (m/s) - 3-D
! wp: w-wind (m/s) - 3-D
! theta: pot'l temp (K) - 3-D
! theta_il: ice-liquid pot'l temp (K) - 3-D
! pi_p: Exner func perturbation - 3-D
! dpidt: Exner func tendency
! rv: water vapor mixing ratio (kg/kg) - 3-D
! r_cloud: cloud water mixing ratio (kg/kg) - 3-D
! r_rain: rain water mixing ratio (kg/kg) - 3-D
! r_snow: snow mixing ratio (kg/kg) - 3-D
! r_hail: hail mixing ratio (kg/kg) - 3-D
! r_p_ice: pristine ice mixing ratio (kg/kg) - 3-D
! r_agg: aggregate mixing ratio (kg/kg) - 3-D
! r_graup: graupel mixing ratio (kg/kg) - 3-D
! rtgt: 1-z_topo/z_top - 2-D
! topo: z-topo - 2-D
! dzt_i: inverse of sigma_z spacing (m^(-1)) - 1-D
! zt: model height (m)
! dxt_i: inverse of delta x spacing (m^(-1)) - 2-D
! dyt_i: inverse of delta y spacing (m^(-1)) - 2-D
! nca_rams: integer counter keeping track of number of time
!      steps that convective tendencies maintained     - 2-D
!      if nca > 0: maintain given convective tendencies
!        from previous timestep and do not call KF scheme
!      if nca < 0: KF scheme is called, if necessary, and
!        convective tendencies updated
! convgo: integer which keeps track if pre-convection checks 
!         satisfied.  If they are, convgo=1 and dumpbucket
!         scheme is not activated at the same grid point. - 2-D              
! w0avg_rams: a running mean average of vertical velocity - 3-D
! w0avglt_rams: a running mean average of the horizontal
!      components of the contravariant vertical velocity - 3-D
! tke: turbulent kinetic energy (only used for Mellor-Yamada option)
! dt:    RAMS model time step in sec - 0-D
!
! OUTPUT:
!
! dthildt: theta_il tendency from K-F CPS (K/s)
! dqtdt: total water mixing ratio tendency from K-F CPS (kg/(kg s))
! dqcdt**: cloud water mixing ratio tendency from K-F CPS (kg/(kg s))
! dqrdt**: rain water mixing ratio tendency from K-F CPS (kg/(kg s))
! dqsdt**: snow mixing ratio tendency from K-F CPS (kg/(kg s))
! dqidt**: ice mixing ratio tendency from K-F CPS (kg/(kg s))
! conprr: convective precip rate (kg/(m^2 s))
!
! ** Used only for full microphysics

!##############################################################################
Subroutine kfdrive (mzp,mxp,myp,ia,iz,ja,jz,imphys,icloud, &
           irain,ipris,isnow,idiffk,       &
           up,vp,theta,theta_il,pi_p,                   &
           pi_0,dpidt,rv,rtgt,topo,dzt_i,zt,       &
           dxt_i,nca_rams,convgo_rams,w0avg_rams,w0avglt_rams,   &
           tke,dt,             &
           dthildt,dqtdt,dqcdt,dqrdt,dqsdt,dqidt,conprr,j_lu)

!The following code to drive the version of the KF scheme used is
!similar to that in the WRF model in 1-D (stand-alone) mode..jsk 8/11/00

use cu_kfeta
use mem_grid, only:iprntstmt,ngrid,print_msg
use mem_micro
use micphys, only:jnmb

implicit none

      integer :: j_lu,imphys,icloud,irain,ipris,isnow,idiffk
      integer :: mzp,mxp,myp,ia,iz,ja,jz
      real, dimension(mzp,mxp,myp) :: up,vp,theta,theta_il,   &
                      pi_p,pi_0,dpidt,rv,w0avg_rams,w0avglt_rams, &
                      tke,dthildt,dqtdt,dqcdt,dqrdt,dqsdt,dqidt

      real, dimension(:,:,:), allocatable :: r_cloud,r_rain,r_p_ice,r_snow &
                                            ,r_agg,r_graup,r_hail

      real, dimension(mxp,myp) :: rtgt, topo, dxt_i, conprr

      integer, dimension(mxp,myp) :: nca_rams,convgo_rams

      real, dimension(mzp) :: dzt_i,zt

      integer :: k200, flag2, k400, flag4, isink, kfgo

      integer :: ka
      real :: topo_t

! intermediate variables

      real :: r_liq, r_sol, press, dpdt, temp_k
! ---------------------------------------------------------------
!
! The following variable declarations for variable which are directly
! ingested by the KF scheme follow the format used for the WRF
! interface
!  
      INTEGER, PARAMETER :: kts=1,kte=500 !max number of k levels
      INTEGER, PARAMETER :: ids=1,ide=1
      INTEGER, PARAMETER :: jds=1,jde=1
      INTEGER, PARAMETER :: kds=1,kde=1
      INTEGER, PARAMETER :: ims=1,ime=1
      INTEGER, PARAMETER :: jms=1,jme=1
      INTEGER, PARAMETER :: kms=1,kme=500 !max number of k levels
      INTEGER, PARAMETER :: its=1,ite=1
      INTEGER, PARAMETER :: jts=1,jte=1
      INTEGER, PARAMETER :: P_QI=1                    ! mods for 1-d version
      INTEGER, PARAMETER :: P_QS=0                    ! mods for 1-d version
      REAL, PARAMETER :: XLV0         = 3.15E6
      REAL, PARAMETER :: XLV1         = 2370.
      REAL, PARAMETER :: XLS0         = 2.905E6
      REAL, PARAMETER :: XLS1         = 259.532
      REAL, PARAMETER :: CP=1005.7,R=287.,G=9.81

      INTEGER :: NTST        ! number of time steps in 600 sec
                             ! formula taken from MM5

      REAL :: DT,DX,DXSQ

      REAL :: TOPO1D

      REAL, DIMENSION( kts:kte ) ::                        &
                                                   DQDT1D, &
                                                  DQIDT1D, &
                                                  DQCDT1D, &
                                                  DQRDT1D, &
                                                  DQSDT1D, &
                                                   DTDT1D

      INTEGER, DIMENSION( ims:ime , jms:jme ) ::      NCA

      REAL, DIMENSION( ims:ime , jms:jme ) ::      RAINCV 

      REAL, DIMENSION( kts:kte ) ::                          &
                                                        U1D, &
                                                        V1D, &
                                                        T1D, &
                                                       DZ1D, &
                                                       QV1D, &
                                                        P1D, &
                                                     DPDT1D, &
                                                    W0AVG1D, &
                                                  W0AVG1DLT, &
                                                      TKE1D, &
                                                       HZ1D, &
                                                       SZ1D
                                                       
      REAL, DIMENSION( kts:kte ) :: DP

      INTEGER :: K
      INTEGER :: II, JJ     ! add integers for outer-loop: 2001-05-16

      REAL :: T00
      REAL    , PARAMETER ::  SVP1=0.6112
      REAL    , PARAMETER ::  SVP2=17.67
      REAL    , PARAMETER ::  SVP3=29.65
      REAL    , PARAMETER ::  SVPT0=273.15

      ! Scratch arrays
      allocate(r_cloud(mzp,mxp,myp))
      allocate(r_rain(mzp,mxp,myp))
      allocate(r_p_ice(mzp,mxp,myp))
      allocate(r_snow(mzp,mxp,myp))
      allocate(r_agg(mzp,mxp,myp))
      allocate(r_graup(mzp,mxp,myp))
      allocate(r_hail(mzp,mxp,myp))

      if(iprntstmt>=1 .and. print_msg) print*,'Start Kain-Fritsch Driver'

      NTST=NINT(600./DT)         ! number of time steps in 600 sec

      T00=273.15


! Loop through all gridpoints to
!  1) Check for convection
!  2) Update convective tendencies

! ----- see if nca_rams is negative, if so call KFPARA ---------

do jj=ja,jz
   do ii=ia,iz

      do k = 1, mzp
        r_cloud(k,ii,jj) = 0.0
        r_rain(k,ii,jj)  = 0.0
        r_p_ice(k,ii,jj) = 0.0
        r_snow(k,ii,jj)  = 0.0
        r_agg(k,ii,jj)   = 0.0
        r_graup(k,ii,jj) = 0.0
        r_hail(k,ii,jj)  = 0.0
        if (jnmb(1) >= 1) r_cloud(k,ii,jj) = micro_g(ngrid)%rcp(k,ii,jj)
        if (jnmb(2) >= 1) r_rain(k,ii,jj)  = micro_g(ngrid)%rrp(k,ii,jj)
        if (jnmb(3) >= 1) r_p_ice(k,ii,jj) = micro_g(ngrid)%rpp(k,ii,jj)
        if (jnmb(4) >= 1) r_snow(k,ii,jj)  = micro_g(ngrid)%rsp(k,ii,jj)
        if (jnmb(5) >= 1) r_agg(k,ii,jj)   = micro_g(ngrid)%rap(k,ii,jj)
        if (jnmb(6) >= 1) r_graup(k,ii,jj) = micro_g(ngrid)%rgp(k,ii,jj)
        if (jnmb(7) >= 1) r_hail(k,ii,jj)  = micro_g(ngrid)%rhp(k,ii,jj)
      enddo
     
      if (nca_rams(ii,jj).gt.0) then
         convgo_rams(ii,jj)=1
      endif

!  If the NCA array is less than 0, then convection is NOT active.
!  If pre-convection checks satified, KF scheme is activated.

      if (nca_rams(ii,jj).le.0) then

         convgo_rams(ii,jj)=0
         dx=1./dxt_i(ii,jj)
         DXSQ = DX*DX

         if (j_lu.eq.1) then
            CALL kf_lutab (SVP1,SVP2,SVP3,SVPT0)
            j_lu=0
         endif

!--------------------------------------------------------------------

! Define arrays to be passed to KF_eta_PARA (1D arrays)
! All arrays interpolated to T-point according to the
! staggered Arakawa C-grid in RAMS. Pressure terms created
! using Exner func

         topo1d=topo(ii,jj)
         topo_t = 0.

         do k = 1, mzp-3
            ka = k+1

            p1d (k)=1.0e5*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)**(cp/r)

            !Saleeby(2006): Included absolute value func and
            ! following if statement to keep positive
            dpdt1d (k)=1.0e5*(abs(dpidt(ka,ii,jj))/cp)**(cp/r)
            if(dpidt(ka,ii,jj).lt.0.) dpdt1d(k)=-dpdt1d(k)

            t1d (k)=theta(ka,ii,jj)*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)
            qv1d(k)=rv(ka,ii,jj)
            if (qv1d(k).lt.0.) qv1d(k)=0.
            u1d (k)=0.5*(up(ka,ii,jj)+up(ka,ii-1,jj))
            v1d (k)=0.5*(vp(ka,ii,jj)+vp(ka,ii,jj-1))
            w0avg1d(k)=0.5*(w0avg_rams(ka-1,ii,jj)+   &
                              w0avg_rams(ka,ii,jj))
            w0avg1dlt(k)=0.5*(w0avglt_rams(ka-1,ii,jj)+  &
                              w0avglt_rams(ka,ii,jj))

!Saleeby(2009):Only try to initialize tke1d if TKE is predicted by Mellor-Yamada
! otherwise segmentation faults can occur since array allocation is not executed
            if(idiffk==1) then
               tke1d(k)=tke(ka,ii,jj)
            else
               tke1d(k)=0.0
            endif

            dz1d(k) =rtgt(ii,jj)/dzt_i(ka)

! Dry and Moist static energy 
            sz1d(k)=cp*t1d(k)+G*(zt(ka)*rtgt(ii,jj)-topo_t)
            hz1d(k)=cp*t1d(k)+G*(zt(ka)*rtgt(ii,jj)-topo_t)+(2.5e6)*qv1d(k)

! Formulation for dp
            if(ka.eq.mzp) then
               dp(k)=(1.0e5*  &
                  ((pi_0(ka-1,ii,jj)+pi_p(ka-1,ii,jj))/cp)**(cp/r)- &
                   1.0e5*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)**(cp/r))
            else
               dp(k)=(1.0e5*  &
                  ((pi_0(ka-1,ii,jj)+pi_p(ka-1,ii,jj))/cp)**(cp/r)+ &
                   1.0e5*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)**(cp/r))/2 -  &
                 (1.0e5*((pi_0(ka+1,ii,jj)+pi_p(ka+1,ii,jj))/cp)**(cp/r) + &
                 1.0e5*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)**(cp/r))/2
            endif

         enddo   ! end k-loop and 1D variable assignments

!---------------------------------------------------------------

! Convection initition checks: if satisfied kfgo=1

         kfgo=0
         flag2=0
         flag4=0
         k200=0
         k400=0

         do k=2,mzp-3
     
!   Determine height of lowest 200-mb

            if(p1d(k).lt.p1d(1)-2.0e4.and.flag2.ne.1) then
               k200=k
               flag2=1
            endif

!   Determine height of the lowest 400-mb

            if(p1d(k).lt.p1d(1)-4.0e4.and.flag4.ne.1) then
               k400=k
               flag4=1
            endif

         enddo

! Determine if air is sinking everywhere in the lowest 400-mb

         isink=1

         do k=1,k400
            if(w0avg1d(k).gt.0) then
               isink=0
            endif
         enddo

! If air is sinking everywhere, require a superadiabatic layer in the
! lowest 200-mb

         if(isink.eq.1) then
            do k=2,k200
               if(sz1d(k).gt.sz1d(k+1)) then
                  kfgo=1
               endif
            enddo

         elseif(isink.eq.0) then
            ! If upward motion, require conditional instability in the lowest 400-mb
            do k=2,k400
               if(hz1d(k).gt.hz1d(k+1)) then
                  kfgo=1
               endif
            enddo
         endif   

         convgo_rams(ii,jj)=kfgo
         NCA(1,1)=nca_rams(ii,jj)

! If pre-convection checks satified AND NCA array < 0, then 
! activate KF scheme

         if(kfgo.eq.1) then
            CALL kf_eta_para (1, 1,                  &
               U1D,V1D,T1D,QV1D,P1D,DZ1D,   &
               W0AVG1D,W0AVG1DLT,TKE1D,TOPO1D,DT,DX,DXSQ,        &
               IMPHYS,ICLOUD,IRAIN,IPRIS,ISNOW,IDIFFK,   &
               XLV0,XLV1,XLS0,XLS1,CP,R,G, &
               DQDT1D,DQIDT1D,DQCDT1D,DQRDT1D,DQSDT1D,DTDT1D, &
               RAINCV,NCA,NTST,      &
               ims,ime, jms,jme,         &
               kts,mzp-3,dp)                

            nca_rams(ii,jj)=NCA(1,1)
         endif

! Update the tendencies if NCA > 0

         if (nca_rams(ii,jj).gt.0) then
 
            do k=1,mzp-2
               ka = k+1

               press=1.0e5*((pi_0(ka,ii,jj)+pi_p(ka,ii,jj))/cp)**(cp/r)

               !Saleeby(2017): Included absolute value func and
               ! following if statement to keep positive in the event
               ! of small negative pressure changes
               dpdt=(abs(dpidt(ka,ii,jj))/cp)**(cp/r)
               if(dpidt(ka,ii,jj).lt.0.) dpdt=-dpdt

               temp_k=theta(ka,ii,jj)*(press/1.0e5)**(r/cp)

               r_liq=r_cloud(ka,ii,jj)+r_rain(ka,ii,jj)

               if (r_liq.lt.0.) r_liq=0.

               r_sol=r_snow(ka,ii,jj)+r_hail(ka,ii,jj)+     &
                 r_p_ice(ka,ii,jj)+r_agg(ka,ii,jj)+     &
                 r_graup(ka,ii,jj)

               if (r_sol.lt.0.) r_sol=0.


! Heating and moisture tendency terms: Full microphysics option
!    Equation for theta_il tendency derived from equation 2.52 
!    in Storm and Cloud Dynamics, by Cotton and Anthes

               if(imphys.ge.2) then

                  if(temp_k.le.253) then

                     dthildt(ka,ii,jj)=((((1.0e5/press)**(r/cp))*         &
                        (temp_k*(-r/cp)*(1/press)*dpdt+dtdt1d(k)))-       &
                        theta_il(ka,ii,jj)*(2.5e6*(dqcdt1d(k)+dqrdt1d(k)) &
                        +2.83e6*(dqidt1d(k)+dqsdt1d(k)))/                 &
                              (cp*amax1(temp_k,253.)))/                   &
                        (1+(2.5e6*r_liq+2.83e6*r_sol)/                    &
                              (cp*amax1(temp_k,253.)))

                  else

                     dthildt(ka,ii,jj)=((((1.0e5/press)**(r/cp))*         &
                        (temp_k*(-r/cp)*(1/press)*dpdt+dtdt1d(k)))-       &
                        theta_il(ka,ii,jj)*((2.5e6*(dqcdt1d(k)+dqrdt1d(k))&
                        +2.83e6*(dqidt1d(k)+dqsdt1d(k)))/(cp*temp_k))     &
                        +theta_il(ka,ii,jj)*((2.5e6*r_liq*dtdt1d(k)+      &
                        2.83e6*r_sol*dtdt1d(k))/                          &
                        (cp*temp_k*temp_k)))/                             &
                        (1+(2.5e6*r_liq+2.83e6*r_sol)/                    &
                        (cp*amax1(temp_k,253.)))

                  endif

                  dqtdt(ka,ii,jj)=dqdt1d(k)+dqidt1d(k)+dqcdt1d(k)  &
                          +dqrdt1d(k)+dqsdt1d(k)

                  dqidt(ka,ii,jj)=dqidt1d(k)
                  dqcdt(ka,ii,jj)=dqcdt1d(k)
                  dqrdt(ka,ii,jj)=dqrdt1d(k)
                  dqsdt(ka,ii,jj)=dqsdt1d(k)
            
               else

! Heating and moisture tendency terms: No microphysics

                  dthildt(ka,ii,jj)=((((1.0e5/press)**(r/cp))*(temp_k*    &
                     (-r/cp)*(1/press)*dpdt+dtdt1d(k))))

                  dqtdt(ka,ii,jj)=dqdt1d(k)
                  dqidt(ka,ii,jj)=0.
                  dqcdt(ka,ii,jj)=0.
                  dqrdt(ka,ii,jj)=0.
                  dqsdt(ka,ii,jj)=0.

               endif
 
            enddo  ! end k-loop

! Convective precipitation rate from KF
            if (raincv(1,1).ge.0.) conprr(ii,jj)=raincv(1,1)

            dthildt(1,ii,jj)=0.
            dthildt(mzp,ii,jj)=0.
            dqtdt(1,ii,jj)=0.
            dqtdt(mzp,ii,jj)=0.

         endif   ! endif for NCA > 0 condition
                 ! when K-F have been activated
                 ! therefore should update tendencies

      endif   ! endif for first NCA < 0 condition

! If NCA < 0 after checking for convection in KF, then tendency
! terms are = 0.

      if (nca_rams(ii,jj).le.0) then
         do k = 1,mzp
            dthildt(k,ii,jj)=0.
            dqtdt(k,ii,jj)=0.
            dqidt(k,ii,jj)=0.
            dqcdt(k,ii,jj)=0.
            dqrdt(k,ii,jj)=0.
            dqsdt(k,ii,jj)=0.
         enddo
         conprr(ii,jj)=0.

      endif

   enddo  ! end ii-loop
enddo   ! end jj-loop

deallocate(r_cloud)
deallocate(r_rain)
deallocate(r_p_ice)
deallocate(r_snow)
deallocate(r_agg)
deallocate(r_graup)
deallocate(r_hail)

return
END SUBROUTINE kfdrive
