!##############################################################################
! KAIN-FRITSCH CUMULUS PARAMETERIZATION SCHEME
! Updated July 2002
!
! This is the actual 1-D version of the Kain-Fritsch CPS code
! used in the ETA model at the National Severe Storms Laboratory.
! 
! NOTES/MODIFICATIONS for RAMS 4.3.0:
!
! The modifications to the convective trigger func formulation
! by Christopher Castro are noted by "CLC RAMS trigger"
! 
! Modifications by Adriana Beltran (abb) are necessary 
! to execute this code in parallel (tested on a linux cluster).  Do not
! remove them!!

Module cu_kfeta

implicit none

save

!--------------------------------------------------------------------
! Lookup table variables:
      INTEGER, PARAMETER :: KFNT=250,KFNP=220
      REAL, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
      REAL, DIMENSION(1:KFNP) :: THE0K
      REAL, DIMENSION(1:200) :: ALU
      REAL :: RDPR,RDTHK,PLUTOP


! Note:  KF Lookup table is used by routines kf_eta_para, TPMIX2,
!        TPMIX2DD, ENVIRTHT
! End of Lookup table variables:

CONTAINS

!##############################################################################
Subroutine kf_eta_para (I,J,U0,V0,T0,QV0,P0,DZQ,W0AVG1D,W0AVG1DLT &
          ,TKE,TOPO,DT,DX,DXSQ,IMPHYS,ICLOUD,IRAIN,IPRIS,ISNOW,IDIFFK  &
          ,XLV0,XLV1,XLS0,XLS1,CP,R,G,DQDT,DQIDT,DQCDT,DQRDT,DQSDT     &
          ,DTDT,RAINCV,NCA,NTST                &
          ,ims,ime,jms,jme,kts,kte,dp)

implicit none

      INTEGER, INTENT(IN   ) :: ims,ime, jms,jme, &
                                kts,kte, &
                                I,J,IMPHYS,ICLOUD,IRAIN,IPRIS, &
                                ISNOW,IDIFFK,NTST

      REAL, DIMENSION( kts:kte ),                          &
            INTENT(IN   ) ::                           U0, &
                                                       V0, &
                                                       T0, &
                                                      QV0, &
                                                       P0, &
                                                      DZQ, &
                                                  W0AVG1D, &
                                                W0AVG1DLT
      REAL,  INTENT(IN   ) :: DT,DX,DXSQ,TOPO

      REAL,  INTENT(IN   ) :: XLV0,XLV1,XLS0,XLS1,CP,R,G
      REAL, DIMENSION( kts:kte ), INTENT(INOUT) ::         &
                                                     DQDT, &
                                                    DQIDT, &
                                                    DQCDT, &
                                                    DQRDT, &
                                                    DQSDT, &
                                                     DTDT

      INTEGER, DIMENSION( ims:ime , jms:jme ),             &
            INTENT(INOUT) ::                          NCA

      REAL, DIMENSION( ims:ime , jms:jme ),                &
            INTENT(INOUT) ::                       RAINCV

!...DEFINE LOCAL VARIABLES...

      REAL, DIMENSION( kts:kte ) ::                        &
            Q0,Z0,TV0,TU,TVU,QU,TZ,TVD,                    &
            QD,QES,TG,TVG,QG,WU,WD,EMS,EMSD,      &
            UMF,UER,UDR,DMF,DER,DDR,UMF2,UER2,             &
            UDR2,DMF2,DER2,DDR2,DZA,THTA0,THETEE,          &
            THTAU,THETEU,THTAD,THETED,QLIQ,QICE,           &
            QLQOUT,QICOUT,PPTLIQ,PPTICE,DETLQ,DETIC,       &
            DETLQ2,DETIC2,RATIO2


      REAL, DIMENSION( kts:kte ) ::                        &
            DOMGDP,EXN,RHOE,TVQU,DP,RH,EQFRC,WSPD,         &
            QDT,FXM,THTAG,THPA,THFXOUT,                    &
            THFXIN,QPA,QFXOUT,QFXIN,QLPA,QLFXIN,           &
            QLFXOUT,QIPA,QIFXIN,QIFXOUT,QRPA,              &
            QRFXIN,QRFXOUT,QSPA,QSFXIN,QSFXOUT,            &
            QL0,QLG,QI0,QIG,QR0,QRG,QS0,QSG


      REAL, DIMENSION( kts:kte+1 ) :: OMG
      REAL, DIMENSION( kts:kte ) :: RAINFB,SNOWFB
      REAL, DIMENSION( kts:kte ) ::                        &
            CLDHGT,QSD,DILFRC,DDILFRC,TKE,TGU,QGU,THTEEG

! LOCAL VARS

      REAL    :: P00,T00,RLF,RHIC,RHBC,PIE,                &
                 TTFRZ,TBFRZ,C5,RATE
      REAL    :: GDRY,ROCP,ALIQ,BLIQ,CLIQ,DLIQ
      REAL    :: FBFRC,P300,DPTHMX,QMIX,ZMIX,PMIX,   &
                 TMIX,EMIX,TLOG,TDPT,TLCL,TVLCL,     &
                 PLCL,ES,DLP,TENV,QENV,TVEN,   &
                 ZLCL,WKL,WKLORIG,WKLTOPO,WKLELEV, &
                 TOPOCOEF,ELEVCOEF,     &
                 WABS,TRPPT,DTLCL,GDT,WLCL,WLCLMIN, WLCLELEV, &
                 DWLCLELEV,PRESTHRS,PRESTHRS2,DPRESTHRS, SLOPETHRS, &
                 TVAVG,WTW,RHOLCL,AU0,VMFLCL,UPOLD,   &
                 UPNEW,ABE,WKLCL,TTEMP,FRC1,   &
                 QNEWIC,RL,BE,BOTERM,ENTERM,&
                 DZZ,REI,EE2,UD2,TTMP,F1,F2,         &
                 THTTMP,QTMP,TMPLIQ,TMPICE,TU95,TU10,EE1,  &
                 UD1,DPTT,QNEWLQ,DUMFDP,           &
                 VCONV,TIMEC,SHSIGN,VWS,PEF, &
                 CBH,RCBH,PEFCBH,PEFF,PEFF2,TDER,   &
                 TADVEC,DPDD,DPT,RDD,A1,     &
                 DSSDT,DTMP,T1RH,QSRH,PPTFLX,CPR,          &
                 AINCM2,            &
                 DDINC,AINCMX,AINCM1, &
                 AINC,TDER2,PPTFL2,FABE,STAB,DTT,DTT1,     &
                 DTIME,TMA,TMB,TMM,BCOEFF,ACOEFF,QVDIFF,   &
                 TOPOMG,CPM,DQ,ABEG,DABE,DFDA,FRC2,DR,     &
                 UDFRC,TUC,QGS,RH0,RHG,QINIT,QINIT1,QINIT2,QFNL,QFNL2,QFNL1,ERR2,    &
                 RLC,RLS,RNC,FABEOLD,AINCOLD,UEFRC, &
                 DDFRC,TDC,DEFRC,RHBAR,DMFFRC,DPMIN,DILBE
   REAL    ::    ASTRT,TP,VALUE,AINTRP,TKEMAX,QFRZ,&
                 QSS,PPTMLT,DTMELT,RHH,EVAC,BINC

   INTEGER :: INDLU,NU,NUCHM,NNN,KLFS
   REAL    :: CHMIN,PM15,CHMAX,DTRH,RAD,DPPP
   REAL    :: U00,QSLCL,RHLCL,DQSSDT

      INTEGER :: KX,K,KL,nkj,ij

      integer abb2(50)
      real abb

      INTEGER :: NCHECK
      INTEGER, DIMENSION (kts:kte) :: KCHECK

      INTEGER :: ISTOP,ML,L5,KMIX,LOW,                     &
                 LC,LLFC,NLAYRS,NK,                 &
                 KPBL,KLCL,LCL,LET,IFLAG,                  &
                 NK1,LTOP,NJ,LTOP1,                        &
                 LTOPM1,KSTART,LFS,               &
                 ND,NIC,LDB,LDT,ND1,NDK,                   &
                 LMAX,NCOUNT,NOITR,                     &
                 NSTEP,NTC,NCHM,ISHALL
      LOGICAL :: IPRNT=.FALSE.
   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPT0=273.15
!
      DATA P00,T00/1.E5,273.15/
      DATA RLF/3.339E5/
      DATA RHIC,RHBC/1.,0.90/
      DATA PIE,TTFRZ,TBFRZ,C5/3.141592654,268.16,248.16,1.0723E-3/
      DATA RATE/0.03/
!-----------------------------------------------------------
      GDRY=-G/CP
      ROCP=R/CP
      KL=kte
      KX=kte
!
!     ALIQ = 613.3
!     BLIQ = 17.502
!     CLIQ = 4780.8
!     DLIQ = 32.19

      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3

! Prevent uninitialized variables by setting to null value
      LDB=0
      L5=0
      KLFS=0
      wklorig=0.
      wkltopo=0.
      wkl=0.
      dtlcl=0.
      chmin=0.
      aincm1=0.
      fabeold=0.
      aincold=0.

!****************************************************************************
!                                                      ! PPT FB MODS
!...OPTION TO FEED CONVECTIVELY GENERATED RAINWATER    ! PPT FB MODS
!...INTO GRID-RESOLVED RAINWATER (OR SNOW/GRAUPEL)     ! PPT FB MODS
!...FIELD.  'FBFRC' IS THE FRACTION OF AVAILABLE       ! PPT FB MODS
!...PRECIPITATION TO BE FED BACK (0.0 - 1.0)...        ! PPT FB MODS
!     print *,'Inside KFPARA...'

! abb modification:
     abb=0

     do ij=1,50
      abb2(ij)=0.
     end do

      FBFRC=0.0                                        ! PPT FB MODS
!...mods to allow shallow convection...
      NCHM = 0
      ISHALL = 0
      DPMIN = 5.E3
!...
      P300=P0(1)-30000.

!...INPUT A VERTICAL SOUNDING ... NOTE THAT MODEL LAYERS ARE NUMBERED
!...FROM BOTTOM-UP IN THE KF SCHEME...
!
      ML=0
      LLFC=0
!SUE  tmprpsb=1./PSB(I,J)
!SUE  CELL=PTOP*tmprpsb

      DO K=1,KX
!
!...IF Q0 IS ABOVE SATURATION VALUE, REDUCE IT TO SATURATION LEVEL...
!
         ES=ALIQ*EXP((BLIQ*T0(K)-CLIQ)/(T0(K)-DLIQ))
         QES(K)=0.622*ES/(P0(K)-ES)
         Q0(K)=AMIN1(QES(K),QV0(K))
         Q0(K)=AMAX1(0.000001,Q0(K))
         QL0(K)=0.
         QI0(K)=0.
         QR0(K)=0.
         QS0(K)=0.
         RH(K) = Q0(K)/QES(K)
         DILFRC(K) = 1.
         TV0(K)=T0(K)*(1.+0.608*Q0(K))
         RHOE(K)=P0(K)/(R*TV0(K))

! IF Turbulent Kinetic Energy (TKE) is available from turbulent mixing scheme
! use it for shallow convection...For now, assume it is not available....
!         TKE(K) = Q2(I,J,NK)
         IF(IDIFFK.ne.1) THEN
           TKE(K) = 0.
         ENDIF
         CLDHGT(K) = 0.
         IF(P0(K).GE.500E2)L5=K
         IF(P0(K).GE.P300)LLFC=K
         IF(T0(K).GT.T00)ML=K
      ENDDO
!
!...DZQ IS DZ BETWEEN SIGMA SURFACES, DZA IS DZ BETWEEN MODEL HALF LEVEL
        Z0(1)=.5*DZQ(1)

        DO K=2,KL
          Z0(K)=Z0(K-1)+.5*(DZQ(K)+DZQ(K-1))
          DZA(K-1)=Z0(K)-Z0(K-1)
        ENDDO   

        DZA(KL)=0.
!
!
!  To save time, specify a pressure interval to move up in sequential
!  check of different ~50 mb deep groups of adjacent model layers in
!  the process of identifying updraft source layer (USL).  Note that 
!  this search is terminated as soon as a buoyant parcel is found and 
!  this parcel can produce a cloud greater than specifed minimum depth
!  (CHMIN)...For now, set interval at 15 mb...
!
       NCHECK = 1
       KCHECK(NCHECK)=1
       PM15 = P0(1)-15.E2
       DO K=2,LLFC
         IF(P0(K).LT.PM15)THEN
           NCHECK = NCHECK+1
           KCHECK(NCHECK) = K
           PM15 = PM15-15.E2
         ENDIF
       ENDDO
!
       NU=0
       NUCHM=0
usl:   DO
           NU = NU+1
           IF(NU.GT.NCHECK)THEN 
             IF(ISHALL.EQ.1)THEN
               CHMAX = 0.
               NCHM = 0
               DO NK = 1,NCHECK
                 NNN=KCHECK(NK)
                 IF(CLDHGT(NNN).GT.CHMAX)THEN
                   NCHM = NNN
                   NUCHM = NK
                   CHMAX = CLDHGT(NNN)
                 ENDIF
               ENDDO
               NU = NUCHM-1
               FBFRC=1.
               CYCLE usl
             ELSE
               RETURN
             ENDIF
           ENDIF      
           KMIX = KCHECK(NU)
           LOW=KMIX
!...
           LC = LOW
!
!...ASSUME THAT IN ORDER TO SUPPORT A DEEP UPDRAFT YOU NEED A LAYER OF
!...UNSTABLE AIR AT LEAST 50 mb DEEP...TO APPROXIMATE THIS, ISOLATE A
!...GROUP OF ADJACENT INDIVIDUAL MODEL LAYERS, WITH THE BASE AT LEVEL
!...LC, SUCH THAT THE COMBINED DEPTH OF THESE LAYERS IS AT LEAST 50 mb..
!   
           NLAYRS=0
           DPTHMX=0.
           NK=LC-1
           DO 
             NK=NK+1   
             DPTHMX=DPTHMX+DP(NK)
             NLAYRS=NLAYRS+1
             IF(DPTHMX.GT.DPMIN)THEN
               EXIT 
             ENDIF
           END DO    
           IF(DPTHMX.LT.DPMIN)THEN 
             RETURN
           ENDIF
           KPBL=LC+NLAYRS-1   
!
!...********************************************************
!...for computational simplicity without much loss in accuracy,
!...mix temperature instead of theta for evaluating convective
!...initiation (triggering) potential...
!          THMIX=0.
           TMIX=0.
           QMIX=0.
           ZMIX=0.
           PMIX=0.

!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
           DO NK=LC,KPBL
             TMIX=TMIX+DP(NK)*T0(NK)
             QMIX=QMIX+DP(NK)*Q0(NK)
             ZMIX=ZMIX+DP(NK)*Z0(NK)
             PMIX=PMIX+DP(NK)*P0(NK)
           ENDDO   
!         THMIX=THMIX/DPTHMX
          TMIX=TMIX/DPTHMX
          QMIX=QMIX/DPTHMX
          ZMIX=ZMIX/DPTHMX
          PMIX=PMIX/DPTHMX
          EMIX=QMIX*PMIX/(0.622+QMIX)

!
!...FIND THE TEMPERATURE OF THE MIXTURE AT ITS LCL...
!
!        TLOG=ALOG(EMIX/ALIQ)
! ...calculate dewpoint using lookup table...
!
          astrt=1.e-3
          ainc=0.075
          a1=emix/aliq
          tp=(a1-astrt)/ainc
          indlu=int(tp)+1
          value=(indlu-1)*ainc+astrt
          aintrp=(a1-value)/ainc
          tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
          TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
          TLCL=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(TMIX-T00))*(TMIX-TDPT)
          TLCL=AMIN1(TLCL,TMIX)
          TVLCL=TLCL*(1.+0.608*QMIX)
!          TVQU(K) = TVLCL
          ZLCL = ZMIX+(TLCL-TMIX)/GDRY
          NK = LC-1
          DO 
            NK = NK+1
            KLCL=NK
            IF(ZLCL.LE.Z0(NK) .or. NK.GT.KL)THEN
              EXIT
            ENDIF 
          ENDDO   
          IF(NK.GT.KL)THEN
            RETURN  
          ENDIF
          K=KLCL-1
          DLP=(ZLCL-Z0(K))/(Z0(KLCL)-Z0(K))
!     
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!     
          TENV=T0(K)+(T0(KLCL)-T0(K))*DLP
          QENV=Q0(K)+(Q0(KLCL)-Q0(K))*DLP
          TVEN=TENV*(1.+0.608*QENV)
!   
!   CONVECTIVE TRIGGER: ORIGINAL DESCRIPTION
!  
!...CHECK TO SEE IF CLOUD IS BUOYANT USING FRITSCH-CHAPPELL TRIGGER
!...FUNC DESCRIBED IN KAIN AND FRITSCH (1992)...W0 IS AN
!...APROXIMATE VALUE FOR THE RUNNING-MEAN GRID-SCALE VERTICAL
!...VELOCITY, WHICH GIVES SMOOTHER FIELDS OF CONVECTIVE INITIATION
!...THAN THE INSTANTANEOUS VALUE...FORMULA RELATING TEMPERATURE
!...PERTURBATION TO VERTICAL VELOCITY HAS BEEN USED WITH THE MOST
!...SUCCESS AT GRID LENGTHS NEAR 25 km.  FOR DIFFERENT GRID-LENGTHS,
!...ADJUST VERTICAL VELOCITY TO EQUIVALENT VALUE FOR 25 KM GRID
!...LENGTH, ASSUMING LINEAR DEPENDENCE OF W ON GRID LENGTH...
! 
!   CLC RAMS TRIGGER MODIFICATIONS
!
!   THE ORIGINAL FORMULATION FOR PARAMETER WKL IS NOTED AS WKLORIG.
!   WKLTOPO TAKES INTO ACCOUNT THE HORIZONTAL COMPONENTS OF THE 
!   CONTRAVARIANT VERTICAL VELOCITY ONLY IF THEY (mathematically)
!   DECREASE THE VALUE OF WKLORIG.  ANY TERRAIN CONTRIBUTION WHICH
!   WOULD INCREASE THE VALUE OF WKL IS NEGLECTED.  WKL CAN BE NEGATIVE
!   IF THE ABSOLUTE VALUE OF WKLTOPO EXCEEDS WKLORIG.  
!
!   RAMS USERS CAN ADJUST THE PARAMETER TOPOCOEF, AS IT WAS ORIGINALLY
!   TUNED FOR A JULY 2001 NORTH AMERICAN MONSOON SIMULATION AT 33 KM
!   GRID SPACING.  TOPOCOEF IMPLICITY ACCOUNTS FOR GRID SPACING SCALING
!   (e.g. ratio DX/25.0 instead of DX/25.E3)

          TOPOCOEF=98000
          ELEVCOEF=0.0
   
          IF(ZLCL.LT.2.E3)THEN
            WKLCL=(0.02*ZLCL/2.E3)
          ELSE
            WKLCL=0.02
          ENDIF



   if(w0avg1dlt(k).le.0.0.and.w0avg1dlt(KLCL).le.0.0) then
     WKLORIG=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
     WKLTOPO=TOPOCOEF*W0AVG1DLT(K)*DLP*(DX/25.E3)*(DX/33.E3)
   endif    


   if(w0avg1dlt(k).le.0.0.and.w0avg1dlt(KLCL).gt.0.0) then
     WKLORIG=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
     WKLTOPO=(-TOPOCOEF*W0AVG1DLT(KLCL)+TOPOCOEF*W0AVG1DLT(K))*DLP &
              *(DX/25.E3)*(DX/33.E3)
   endif

   if(w0avg1dlt(k).gt.0.0.and.w0avg1dlt(KLCL).le.0.0) then
     WKLORIG=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
     WKLTOPO=(-TOPOCOEF*W0AVG1DLT(K))*(DX/25.E3)*(DX/33.E3)          
   endif

   if(w0avg1dlt(k).gt.0.0.and.w0avg1dlt(KLCL).gt.0.0) then
     WKLORIG=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
     WKLTOPO=(-TOPOCOEF*W0AVG1DLT(K)+(-TOPOCOEF*W0AVG1DLT(KLCL))*DLP) &
              *(DX/25.E3)*(DX/33.E3)
   endif


          PRESTHRS=850*1.E2
          DPRESTHRS=25*1.E2
          SLOPETHRS=ABS(WKLTOPO/WKLORIG)

          WKLELEV=ELEVCOEF*(0.02*(1.E5/P0(1))-0.02)*exp((1.E5/P0(1)-1))

          if(P0(1).LE.PRESTHRS-DPRESTHRS.OR.SLOPETHRS.GE.0.5) then
             WKLELEV=WKLELEV
          elseif(P0(1).GE.PRESTHRS+DPRESTHRS) then
             WKLELEV=0.0
          else
             WKLELEV=(-(1/(2*DPRESTHRS))*P0(1)+0.5+(1/(2*DPRESTHRS)) &
                  *PRESTHRS)*WKLELEV
          endif

          PRESTHRS2=925*1.E2
        
          if(P0(1).GE.PRESTHRS2+DPRESTHRS.AND.SLOPETHRS.LE.0.01) then
             WKLTOPO=0.
          elseif(P0(1).LE.PRESTHRS2-DPRESTHRS) then
             WKLTOPO=WKLTOPO
          else
             if(SLOPETHRS.LE.0.01) then
               WKLTOPO=(-(1/(2*DPRESTHRS))*P0(1)+0.5+(1/(2*DPRESTHRS)) &
                   *PRESTHRS2)*WKLTOPO
             else
               WKLTOPO=WKLTOPO
             endif
          endif


          if(WKLTOPO.gt.0) then
             WKLTOPO=0
          endif
          if(WKLORIG.LE.0) then
             WKL=WKLTOPO*exp(3*(1.E5/P0(1)-1))
          endif
          if(WKLORIG.GT.0) then
             WKL=(WKLORIG+WKLTOPO*exp(3*(1.E5/P0(1)-1)))
          endif

          WKL=WKL-WKLELEV
          if(WKL.LE.0) then
            DTLCL=-4.64*(ABS(WKL))**0.33
          endif
          if(WKL.GT.0) then
            DTLCL=4.64*WKL**0.33
          endif

! Original code for convective trigger func
            
!             WKL=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
!             WABS=ABS(WKL)
!             IF(WABS.EQ.0.)THEN
!               WSIGNE=1.
!               DTLCL=0.
!             ELSE 
!               WSIGNE=WKL/WABS
!               DTLCL=4.64*WKL**0.33
!               DTLCL=4.64*WSIGNE*WABS**0.33
!               DTLCL = AMAX1(DTLCL,0.)
!             ENDIF



!...for ETA model, give parcel an extra temperature perturbation based
!...the threshold RH for condensation (U00)...
!
!...for now, just assume U00=0.75...
!...!!!!!! for MM5, SET DTRH = 0. !!!!!!!!
         U00 = 0.75
         IF(U00.lt.1.)THEN
           QSLCL=QES(K)+(QES(KLCL)-QES(K))*DLP
           RHLCL = QENV/QSLCL
           DQSSDT = QMIX*(CLIQ-BLIQ*DLIQ)/((TLCL-DLIQ)*(TLCL-DLIQ))
           IF(RHLCL.ge.0.75 .and. RHLCL.le.0.95)then
             DTRH = 0.25*(RHLCL-0.75)*QMIX/DQSSDT
           ELSEIF(RHLCL.GT.0.95)THEN
             DTRH = (1./RHLCL-1.)*QMIX/DQSSDT
           ELSE
               DTRH = 0.
           ENDIF
         ENDIF   
!         IF(ISHALL.EQ.1)IPRNT=.TRUE.


          IPRNT=.TRUE.

       IF(IPRNT)THEN
        TVAVG=0.5*(TV0(KLCL)+TENV*(1.+0.608*QENV))
        PLCL=P0(KLCL)*EXP(G/(R*TVAVG)*(Z0(KLCL)-ZLCL))
       ENDIF


!         IF(TLCL+DTLCL.GT.TENV)GOTO 45
!

trigger:  IF(TLCL+DTLCL+DTRH.LT.TENV)THEN   
!trigger:  IF(TLCL+DTLCL+DTRH.LT.TENV)THEN   
!
! Parcel not buoyant, CYCLE back to start of trigger and evaluate next potential USL...
!
! Check to see if in the lowest 300-mb            
           IF(KMIX.LE.LLFC) then
             CYCLE usl
           ENDIF
!
          ELSE                            ! Parcel is buoyant, determine updraft
!     
!...CONVECTIVE TRIGGERING CRITERIA HAS BEEN SATISFIED...COMPUTE
!...EQUIVALENT POTENTIAL TEMPERATURE
!...(THETEU) AND VERTICAL VELOCITY OF THE RISING PARCEL AT THE LCL...
!     
            CALL envirtht (PMIX,TMIX,QMIX,THETEU(K),ALIQ,BLIQ,CLIQ,DLIQ)
!
!...modify calculation of initial parcel vertical velocity...jsk 11/26/97

! MM4 formulation for GDT: Used for RAMS

              GDT=G*(DTLCL+DTRH)*(ZLCL-Z0(LC))/(TV0(LC)+TVEN)

! CLC RAMS TRIGGER: Formulation for WLCL
!
! New formulation for WLCL considers the minumum value of WLCL to vary
! approximately as a step func with an elevation threshold of 1500 m
! (WLCLELEV), so it is effectively zero above this level.
! Necessary because the previous minimum value of WLCL caused too much
! precipitation at high elevations.  Users may wish to modify the
! threshold if high elevation precipitation is over or underestimated.


              WLCLELEV=1500
              DWLCLELEV=150

              if(TOPO.GE.WLCLELEV+DWLCLELEV) then
                WLCLMIN=1.E-3
              elseif(TOPO.LE.WLCLELEV-DWLCLELEV) then
                WLCLMIN=1.0
              else
                WLCLMIN=-(1/(2*DWLCLELEV))*TOPO+0.5+(1/(2*DWLCLELEV))*WLCLELEV
              endif

              if(GDT.LE.0) then
                WLCL=WLCLMIN
              else
                WLCL=(WLCLMIN+.5*SQRT(ABS(GDT)+1.E-10))
              endif

! Original formulation for WLCL
!
!                 WLCL=(1.+.5*SQRT(ABS(GDT)+1.E-10))


! Eta formulation for GDT and WLCL
!            GDT=2.*G*(DTLCL+DTRH)*500./TVEN
!            WLCL=1.+0.5*WSIGNE*SQRT(ABS(GDT)+1.E-10)
!            WLCL=1.+0.5*SQRT(ABS(GDT)+1.E-10)

! MM4 does not use these
!           WLCL = AMIN1(WLCL,5.)
!           WLCL = AMIN1(WLCL,4.)
!           WLCL = AMIN1(WLCL,3.)
            PLCL=P0(K)+(P0(KLCL)-P0(K))*DLP

            WTW=WLCL*WLCL

            TVLCL=TLCL*(1.+0.608*QMIX)
            RHOLCL=PLCL/(R*TVLCL)
!        
            LCL=KLCL
            LET=LCL
! make RAD a func of background vertical velocity...

! CLC RAMS TRIGGER
! Also, if WLCL is less than one, set the updraft radius near zero
! so no convection occurs.

          IF(WLCL.ge.1.0) THEN
            IF(WKL.LT.0)THEN
               RAD = 1000
            ELSEIF(WKL.GT.0.1)THEN
               RAD = 2000
            ELSE
              RAD = (1000.+1000*WKL/0.1)
            ENDIF
          ELSE
            RAD=1.E-4
          ENDIF
!     
!*******************************************************************
!                                                                  *
!                 COMPUTE UPDRAFT PROPERTIES                       *
!                                                                  *
!*******************************************************************
!     
!     
!...
!...ESTIMATE INITIAL UPDRAFT MASS FLUX (UMF(K))...
!     
            WU(K)=WLCL
            AU0=PIE*RAD*RAD
            UMF(K)=RHOLCL*AU0
            VMFLCL=UMF(K)
            UPOLD=VMFLCL
            UPNEW=UPOLD
!     
!...RATIO2 IS THE DEGREE OF GLACIATION IN THE CLOUD (0 TO 1),
!...UER IS THE ENVIR ENTRAINMENT RATE, ABE IS AVAILABLE
!...BUOYANT ENERGY, TRPPT IS THE TOTAL RATE OF PRECIPITATION
!...PRODUCTION...
!     
            RATIO2(K)=0.
            UER(K)=0.
            ABE=0.
            TRPPT=0.
            TU(K)=TLCL
            TVU(K)=TVLCL
            QU(K)=QMIX
            EQFRC(K)=1.
            QLIQ(K)=0.
            QICE(K)=0.
            QLQOUT(K)=0.
            QICOUT(K)=0.
            DETLQ(K)=0.
            DETIC(K)=0.
            PPTLIQ(K)=0.
            PPTICE(K)=0.

            IFLAG=0
!     
!...TTEMP IS USED DURING CALCULATION OF THE LINEAR GLACIATION
!...PROCESS; IT IS INITIALLY SET TO THE TEMPERATURE AT WHICH
!...FREEZING IS SPECIFIED TO BEGIN.  WITHIN THE GLACIATION
!...INTERVAL, IT IS SET EQUAL TO THE UPDRAFT TEMP AT THE
!...PREVIOUS MODEL LEVEL...
!     
            TTEMP=TTFRZ
!     
!...ENTER THE LOOP FOR UPDRAFT CALCULATIONS...CALCULATE UPDRAFT TEMP,
!...MIXING RATIO, VERTICAL MASS FLUX, LATERAL DETRAINMENT OF MASS AND
!...MOISTURE, PRECIPITATION RATES AT EACH MODEL LEVEL...
!     
!     
            EE1=1.
            UD1=0.
            REI = 0.
            DILBE = 0.
updraft:    DO NK=K,KL-1
              NK1=NK+1
              RATIO2(NK1)=RATIO2(NK)
              FRC1=0.
              TU(NK1)=T0(NK1)
              THETEU(NK1)=THETEU(NK)
              QU(NK1)=QU(NK)
              QLIQ(NK1)=QLIQ(NK)
              QICE(NK1)=QICE(NK)

 CALL tpmix2 (p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1),   &
                     qice(nk1),qnewlq,qnewic,XLV1,XLV0)

!...CHECK TO SEE IF UPDRAFT TEMP IS ABOVE THE TEMPERATURE AT WHICH
!...GLACIATION IS ASSUMED TO INITIATE; IF IT IS, CALCULATE THE
!...FRACTION OF REMAINING LIQUID WATER TO FREEZE...TTFRZ IS THE
!...TEMP AT WHICH FREEZING BEGINS, TBFRZ THE TEMP BELOW WHICH ALL
!...LIQUID WATER IS FROZEN AT EACH LEVEL...
!
              IF(TU(NK1).LE.TTFRZ)THEN
                IF(TU(NK1).GT.TBFRZ)THEN
                  IF(TTEMP.GT.TTFRZ)TTEMP=TTFRZ
                  FRC1=(TTEMP-TU(NK1))/(TTEMP-TBFRZ)
                ELSE
                  FRC1=1.
                  IFLAG=1
                ENDIF
                TTEMP=TU(NK1)
!
!  DETERMINE THE EFFECTS OF LIQUID WATER FREEZING WHEN TEMPERATURE
!...IS BELOW TTFRZ...
!
                QFRZ = (QLIQ(NK1)+QNEWLQ)*FRC1
                QNEWIC=QNEWIC+QNEWLQ*FRC1
                QNEWLQ=QNEWLQ-QNEWLQ*FRC1
                QICE(NK1) = QICE(NK1)+QLIQ(NK1)*FRC1
                QLIQ(NK1) = QLIQ(NK1)-QLIQ(NK1)*FRC1

                CALL dtfrznew (TU(NK1),P0(NK1),THETEU(NK1),QU(NK1),QFRZ,         &
                          QICE(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
              ENDIF

              TVU(NK1)=TU(NK1)*(1.+0.608*QU(NK1))
!
!  CALCULATE UPDRAFT VERTICAL VELOCITY AND PRECIPITATION FALLOUT...
!
              IF(NK.EQ.K)THEN
                BE=(TVLCL+TVU(NK1))/(TVEN+TV0(NK1))-1.
                BOTERM=2.*(Z0(NK1)-ZLCL)*G*BE/1.5
! mm4 mods
                ENTERM=0.
                DZZ=Z0(NK1)-ZLCL
              ELSE
                BE=(TVU(NK)+TVU(NK1))/(TV0(NK)+TV0(NK1))-1.
                BOTERM=2.*DZA(NK)*G*BE/1.5
! mm4 mods
                ENTERM=2.*UER(NK)*WTW/UPOLD
                DZZ=DZA(NK)
              ENDIF
! Eta mods
!              ENTERM=2.*REI*WTW/UPOLD


              CALL condload (QLIQ(NK1),QICE(NK1),WTW,DZZ,BOTERM,ENTERM,      &
                        RATE,QNEWLQ,QNEWIC,QLQOUT(NK1),QICOUT(NK1),G)

              WABS=SQRT(ABS(WTW))
              WU(NK1)=WTW/WABS

!
!...IF VERT VELOCITY IS LESS THAN ZERO, EXIT THE UPDRAFT LOOP AND,
!...IF CLOUD IS TALL ENOUGH, FINALIZE UPDRAFT CALCULATIONS...
!
              IF(WU(NK1).LT.0.)THEN
                EXIT
              ENDIF
!
!...Calculate value of THETA-E in environment to entrain into updraft...
!
              CALL envirtht (P0(NK1),T0(NK1),Q0(NK1),THETEE(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
!
!...REI IS THE RATE OF ENVIRONMENTAL INFLOW...
!
! abb Original REI value
              REI=VMFLCL*DP(NK1)*0.03/RAD

! abb Doubled REI value, Mapes suggestion....

!               REI=(VMFLCL*DP(NK1)*0.03/RAD)*4.0

! abb End 
              TVQU(NK1)=TU(NK1)*(1.+0.608*QU(NK1)-QLIQ(NK1)-QICE(NK1))

              IF(NK.EQ.K)THEN
                DILBE=((TVLCL+TVQU(NK1))/(TVEN+TV0(NK1))-1.)*DZZ
              ELSE
                DILBE=((TVQU(NK)+TVQU(NK1))/(TV0(NK)+TV0(NK1))-1.)*DZZ
              ENDIF
              IF(DILBE.GT.0.)ABE=ABE+DILBE*G
!
!...IF CLOUD PARCELS ARE VIRTUALLY COLDER THAN THE ENVIRONMENT, MINIMAL 
!...ENTRAINMENT (0.5*REI) IS IMPOSED...
!
              IF(TVQU(NK1).LE.TV0(NK1))THEN    ! Entrain/Detrain IF BLOCK
                EE2=0.5
                UD2=1.
                EQFRC(NK1)=0.
              ELSE
                LET=NK1
                TTMP=TVQU(NK1)
!
!...DETERMINE THE CRITICAL MIXED FRACTION OF UPDRAFT AND ENVIRONMENTAL AIR...
!
                F1=0.95
                F2=1.-F1
                THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
                QTMP=F1*Q0(NK1)+F2*QU(NK1)
                TMPLIQ=F2*QLIQ(NK1)
                TMPICE=F2*QICE(NK1)
                CALL tpmix2 (p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,        &
                           qnewlq,qnewic,XLV1,XLV0)
                TU95=TTMP*(1.+0.608*QTMP-TMPLIQ-TMPICE)
                IF(TU95.GT.TV0(NK1))THEN
                  EE2=1.
                  UD2=0.
                  EQFRC(NK1)=1.0
                ELSE
                  F1=0.10
                  F2=1.-F1
                  THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
                  QTMP=F1*Q0(NK1)+F2*QU(NK1)
                  TMPLIQ=F2*QLIQ(NK1)
                  TMPICE=F2*QICE(NK1)
                  CALL tpmix2 (p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,        &
                               qnewlq,qnewic,XLV1,XLV0)
                  TU10=TTMP*(1.+0.608*QTMP-TMPLIQ-TMPICE)
                  IF(TU10.EQ.TVQU(NK1))THEN
                    EE2=1.
                    UD2=0.
                    EQFRC(NK1)=1.0
                  ELSE
                    EQFRC(NK1)=(TV0(NK1)-TVQU(NK1))*F1/(TU10-TVQU(NK1))
                    EQFRC(NK1)=AMAX1(0.0,EQFRC(NK1))
                    EQFRC(NK1)=AMIN1(1.0,EQFRC(NK1))
                    IF(EQFRC(NK1).EQ.1)THEN
                      EE2=1.
                      UD2=0.
                    ELSEIF(EQFRC(NK1).EQ.0.)THEN
                      EE2=0.
                      UD2=1.
                    ELSE
!
!...ROUTINE PROF5 INTEGRATES OVER THE GAUSSIAN DIST TO DETERMINE THE
!   FRACTIONAL ENTRAINMENT AND DETRAINMENT RATES...
!
                      CALL prof5 (EQFRC(NK1),EE2,UD2)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF                            ! End of Entrain/Detrain IF BLOCK
!
!
!...NET ENTRAINMENT AND DETRAINMENT RATES ARE GIVEN BY THE AVERAGE FRACTIONAL
!   VALUES IN THE LAYER...
!
              EE2 = AMAX1(EE2,0.5)
              UD2 = 1.5*UD2
              UER(NK1)=0.5*REI*(EE1+EE2)
              UDR(NK1)=0.5*REI*(UD1+UD2)
!
!...IF THE CALCULATED UPDRAFT DETRAINMENT RATE IS GREATER THAN THE TOTAL
!   UPDRAFT MASS FLUX, ALL CLOUD MASS DETRAINS, EXIT UPDRAFT CALCULATIONS...
!
              IF(UMF(NK)-UDR(NK1).LT.10.)THEN
!
!...IF THE CALCULATED DETRAINED MASS FLUX IS GREATER THAN THE TOTAL UPD MASS
!   FLUX, IMPOSE TOTAL DETRAINMENT OF UPDRAFT MASS AT THE PREVIOUS MODEL LVL..
!   First, correct ABE calculation if needed...
!
                IF(DILBE.GT.0.)THEN
                  ABE=ABE-DILBE*G
                ENDIF
                LET=NK
                EXIT 
              ELSE
                EE1=EE2
                UD1=UD2
                UPOLD=UMF(NK)-UDR(NK1)
                UPNEW=UPOLD+UER(NK1)
                UMF(NK1)=UPNEW
                DILFRC(NK1) = UPNEW/UPOLD
!
!...DETLQ AND DETIC ARE THE RATES OF DETRAINMENT OF LIQUID AND
!...ICE IN THE DETRAINING UPDRAFT MASS...
!
                DETLQ(NK1)=QLIQ(NK1)*UDR(NK1)
                DETIC(NK1)=QICE(NK1)*UDR(NK1)
                QDT(NK1)=QU(NK1)
                QU(NK1)=(UPOLD*QU(NK1)+UER(NK1)*Q0(NK1))/UPNEW
                THETEU(NK1)=(THETEU(NK1)*UPOLD+THETEE(NK1)*UER(NK1))/UPNEW
                QLIQ(NK1)=QLIQ(NK1)*UPOLD/UPNEW
                QICE(NK1)=QICE(NK1)*UPOLD/UPNEW
!
!...PPTLIQ IS THE RATE OF GENERATION (FALLOUT) OF
!...LIQUID PRECIP AT A GIVEN MODEL LVL, PPTICE THE SAME FOR ICE,
!...TRPPT IS THE TOTAL RATE OF PRODUCTION OF PRECIP UP TO THE
!...CURRENT MODEL LEVEL...
!
                PPTLIQ(NK1)=QLQOUT(NK1)*UMF(NK)
                PPTICE(NK1)=QICOUT(NK1)*UMF(NK)
!
                TRPPT=TRPPT+PPTLIQ(NK1)+PPTICE(NK1)
       
                IF(NK1.LE.KPBL)UER(NK1)=UER(NK1)+VMFLCL*DP(NK1)/DPTHMX
              ENDIF
!
            END DO updraft
!
!...CHECK CLOUD DEPTH...IF CLOUD IS TALL ENOUGH, ESTIMATE THE EQUILIBRIU
!   TEMPERATURE LEVEL (LET) AND ADJUST MASS FLUX PROFILE AT CLOUD TOP SO
!   THAT MASS FLUX DECREASES TO ZERO AS A LINEAR FUNC OF PRESSURE BE
!   THE LET AND CLOUD TOP...
!     
!...LTOP IS THE MODEL LEVEL JUST BELOW THE LEVEL AT WHICH VERTICAL VELOC
!   FIRST BECOMES NEGATIVE...
!     
            LTOP=NK


            CLDHGT(LC)=Z0(LTOP)-ZLCL 
!...Instead of using the same minimum cloud height (for deep convection)
!...everywhere, try specifying minimum cloud depth as a func of TLCL...
!
            IF(TLCL.GT.293.)THEN
              CHMIN = 4.E3
            ELSEIF(TLCL.LE.293. .and. TLCL.GE.273)THEN
              CHMIN = 2.E3 + 100.*(TLCL-273.)
            ELSEIF(TLCL.LT.273.)THEN
              CHMIN = 2.E3
            ENDIF

!     
!...If cloud top height is less than the specified minimum for deep 
!...convection, save value to consider this level as source for 
!...shallow convection, go back up to check next level...
!     
!...Try specifying minimum cloud depth as a func of TLCL...
!
!
!...DO NOT ALLOW ANY CLOUD FROM THIS LAYER IF:
!
!...		1.) if there is no CAPE, or 
!...		2.) cloud top is at model level just above LCL, or
!...            3.) cloud top is within updraft source layer, or
!...            4.) cloud-top detrainment layer begins within 
!...                updraft source layer.
!
! No Convection Allowed
  IF(LTOP.LE.KLCL .or. LTOP.LE.KPBL .or. LET+1.LE.KPBL)THEN 
              CLDHGT(LC)=0.
              DO NK=K,LTOP
                UMF(NK)=0.
                UDR(NK)=0.
                UER(NK)=0.
                DETLQ(NK)=0.
                DETIC(NK)=0.
                PPTLIQ(NK)=0.
                PPTICE(NK)=0.
              ENDDO
!        
   ELSEIF(CLDHGT(LC).GT.CHMIN .and. ABE.GT.1)THEN      ! Deep Convection allowed
              ISHALL=0
              EXIT usl
            ELSE
!
!...TO DISALLOW SHALLOW CONVECTION, COMMENT OUT NEXT LINE !!!!!!!!
              ISHALL = 1
              IF(NU.EQ.NUCHM)THEN
                EXIT usl               ! Shallow Convection from this layer
              ELSE
! Remember this layer (by virtue of non-zero CLDHGT) 
! as potential shallow-cloud layer
                DO NK=K,LTOP
                  UMF(NK)=0.
                  UDR(NK)=0.
                  UER(NK)=0.
                  DETLQ(NK)=0.
                  DETIC(NK)=0.
                  PPTLIQ(NK)=0.
                  PPTICE(NK)=0.
                ENDDO
              ENDIF
            ENDIF
          ENDIF trigger
        END DO usl
    IF(ISHALL.EQ.1)THEN
      KSTART=MAX0(KPBL,KLCL)
      LET=KSTART
    endif
!     
!...IF THE LET AND LTOP ARE THE SAME, DETRAIN ALL OF THE UPDRAFT MASS FL
!   THIS LEVEL...
!     
    IF(LET.EQ.LTOP)THEN
      UDR(LTOP)=UMF(LTOP)+UDR(LTOP)-UER(LTOP)
      DETLQ(LTOP)=QLIQ(LTOP)*UDR(LTOP)*UPNEW/UPOLD
      DETIC(LTOP)=QICE(LTOP)*UDR(LTOP)*UPNEW/UPOLD

      UER(LTOP)=0.
      UMF(LTOP)=0.
    ELSE 
!     
!   BEGIN TOTAL DETRAINMENT AT THE LEVEL ABOVE THE LET...
!     
      DPTT=0.
      DO NJ=LET+1,LTOP
        DPTT=DPTT+DP(NJ)
      ENDDO
      DUMFDP=UMF(LET)/DPTT
!     
!...ADJUST MASS FLUX PROFILES, DETRAINMENT RATES, AND PRECIPITATION FALL
!   RATES TO REFLECT THE LINEAR DECREASE IN MASS FLX BETWEEN THE LET AND
!     
      DO NK=LET+1,LTOP
!
!...entrainment is allowed at every level except for LTOP, so disallow
!...entrainment at LTOP and adjust entrainment rates between LET and LTOP
!...so the the dilution factor due to entyrianment is not changed but 
!...the actual entrainment rate will change due due forced total 
!...detrainment in this layer...
!
        IF(NK.EQ.LTOP)THEN
          UDR(NK) = UMF(NK-1)
          UER(NK) = 0.
          DETLQ(NK) = UDR(NK)*QLIQ(NK)*DILFRC(NK)
          DETIC(NK) = UDR(NK)*QICE(NK)*DILFRC(NK)
        ELSE
          UMF(NK)=UMF(NK-1)-DP(NK)*DUMFDP
          UER(NK)=UMF(NK)*(1.-1./DILFRC(NK))
          UDR(NK)=UMF(NK-1)-UMF(NK)+UER(NK)
          DETLQ(NK)=UDR(NK)*QLIQ(NK)*DILFRC(NK)
          DETIC(NK)=UDR(NK)*QICE(NK)*DILFRC(NK)
        ENDIF
        IF(NK.GE.LET+2)THEN
          TRPPT=TRPPT-PPTLIQ(NK)-PPTICE(NK)
          PPTLIQ(NK)=UMF(NK-1)*QLQOUT(NK)
          PPTICE(NK)=UMF(NK-1)*QICOUT(NK)
          TRPPT=TRPPT+PPTLIQ(NK)+PPTICE(NK)
        ENDIF
      ENDDO
    ENDIF
!     
! Initialize some arrays below cloud base and above cloud top...
!
    DO NK=1,K
      IF(NK.GE.LC)THEN
        IF(NK.EQ.LC)THEN
          UMF(NK)=VMFLCL*DP(NK)/DPTHMX
          UER(NK)=VMFLCL*DP(NK)/DPTHMX
        ELSEIF(NK.LE.KPBL)THEN
          UER(NK)=VMFLCL*DP(NK)/DPTHMX
          UMF(NK)=UMF(NK-1)+UER(NK)
        ELSE
          UMF(NK)=VMFLCL
          UER(NK)=0.
        ENDIF
        TU(NK)=TMIX+(Z0(NK)-ZMIX)*GDRY
        QU(NK)=QMIX
        WU(NK)=WLCL
      ELSE
        TU(NK)=0.
        QU(NK)=0.
        UMF(NK)=0.
        WU(NK)=0.
        UER(NK)=0.
      ENDIF
      UDR(NK)=0.
      QDT(NK)=0.
      QLIQ(NK)=0.
      QICE(NK)=0.
      QLQOUT(NK)=0.
      QICOUT(NK)=0.
      PPTLIQ(NK)=0.
      PPTICE(NK)=0.
      DETLQ(NK)=0.
      DETIC(NK)=0.
      RATIO2(NK)=0.
      CALL envirtht (P0(NK),T0(NK),Q0(NK),THETEE(NK),ALIQ,BLIQ,CLIQ,DLIQ)
      EQFRC(NK)=1.0
    ENDDO
!     
      LTOP1=LTOP+1
      LTOPM1=LTOP-1
!     
!...DEFINE VARIABLES ABOVE CLOUD TOP...
!     
      DO NK=LTOP1,KX
        UMF(NK)=0.
        UDR(NK)=0.
        UER(NK)=0.
        QDT(NK)=0.
        QLIQ(NK)=0.
        QICE(NK)=0.
        QLQOUT(NK)=0.
        QICOUT(NK)=0.
        DETLQ(NK)=0.
        DETIC(NK)=0.
        PPTLIQ(NK)=0.
        PPTICE(NK)=0.
        IF(NK.GT.LTOP1)THEN
          TU(NK)=0.
          QU(NK)=0.
          WU(NK)=0.
        ENDIF
        THTA0(NK)=0.
        THTAU(NK)=0.
        EMS(NK)=0.
        EMSD(NK)=0.
        TG(NK)=T0(NK)
        QG(NK)=Q0(NK)
        QLG(NK)=0.
        QIG(NK)=0.
        QRG(NK)=0.
        QSG(NK)=0.
        OMG(NK)=0.

      ENDDO
        OMG(KX+1)=0.

        DO NK=1,LTOP
          EMS(NK)=DP(NK)*DXSQ/G
          EMSD(NK)=1./EMS(NK)
!     
!...INITIALIZE SOME VARIABLES TO BE USED LATER IN THE VERT ADVECTION SCH
!     
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*QDT(NK)))
          THTAU(NK)=TU(NK)*EXN(NK)
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*Q0(NK)))
          THTA0(NK)=T0(NK)*EXN(NK)
          DDILFRC(NK) = 1./DILFRC(NK)
          OMG(NK)=0.
        ENDDO
!     
!...COMPUTE CONVECTIVE TIME SCALE(TIMEC). THE MEAN WIND AT THE LCL
!...AND MIDTROPOSPHERE IS USED.
!     
        WSPD(KLCL)=SQRT(U0(KLCL)*U0(KLCL)+V0(KLCL)*V0(KLCL))
        WSPD(L5)=SQRT(U0(L5)*U0(L5)+V0(L5)*V0(L5))
        WSPD(LTOP)=SQRT(U0(LTOP)*U0(LTOP)+V0(LTOP)*V0(LTOP))
        VCONV=.5*(WSPD(KLCL)+WSPD(L5))
!...for ETA model, DX is a func of location...
!       TIMEC=DX(I,J)/VCONV
        TIMEC=DX/VCONV
        TADVEC=TIMEC
        TIMEC=AMAX1(1800.,TIMEC)
        TIMEC=AMIN1(3600.,TIMEC)
        IF(ISHALL.EQ.1)TIMEC=2400.
        NIC=NINT(TIMEC/DT)
        TIMEC=FLOAT(NIC)*DT
!     
!...COMPUTE WIND SHEAR AND PRECIPITATION EFFICIENCY.
!     
        IF(WSPD(LTOP).GT.WSPD(KLCL))THEN
          SHSIGN=1.
        ELSE
          SHSIGN=-1.
        ENDIF
        VWS=(U0(LTOP)-U0(KLCL))*(U0(LTOP)-U0(KLCL))+(V0(LTOP)-V0(KLCL))*   &
            (V0(LTOP)-V0(KLCL))
        VWS=1.E3*SHSIGN*SQRT(VWS)/(Z0(LTOP)-Z0(LCL))
        PEF=1.591+VWS*(-.639+VWS*(9.53E-2-VWS*4.96E-3))
        PEF=AMAX1(PEF,.2)
        PEF=AMIN1(PEF,.9)
!     
!...PRECIPITATION EFFICIENCY IS A FUNC OF THE HEIGHT OF CLOUD BASE.
!     
        CBH=(ZLCL-Z0(1))*3.281E-3
        IF(CBH.LT.3.)THEN
          RCBH=.02
        ELSE
          RCBH=.96729352+CBH*(-.70034167+CBH*(.162179896+CBH*(-            &
               1.2569798E-2+CBH*(4.2772E-4-CBH*5.44E-6))))
        ENDIF
        IF(CBH.GT.25)RCBH=2.4
        PEFCBH=1./(1.+RCBH)
        PEFCBH=AMIN1(PEFCBH,.9)
!     
!... MEAN PEF. IS USED TO COMPUTE RAINFALL.
!     
        PEFF=.5*(PEF+PEFCBH)
        PEFF2 = PEFF                                ! JSK MODS
       IF(IPRNT)THEN  
       endif     
!        WRITE(98,1035)PEF,PEFCBH,LC,LET,WKL,VWS
!*****************************************************************
!                                                                *
!                  COMPUTE DOWNDRAFT PROPERTIES                  *
!                                                                *
!*****************************************************************
!     
!     
       TDER=0.
 devap:IF(ISHALL.EQ.1)THEN
         LFS = 1
       ELSE
!
!...start downdraft about 150 mb above cloud base...
!
!        KSTART=MAX0(KPBL,KLCL)
!        KSTART=KPBL                                  ! Changed 7/23/99
         KSTART=KPBL+1                                ! Changed 7/23/99
         DO NK = KSTART+1,KL
           DPPP = P0(KSTART)-P0(NK)
!          IF(DPPP.GT.200.E2)THEN
           IF(DPPP.GT.150.E2)THEN
             KLFS = NK
             EXIT 
           ENDIF
         ENDDO
         KLFS = MIN0(KLFS,LET-1)
         LFS = KLFS

!
!...if LFS is not at least 50 mb above cloud base (implying that the 
!...level of equil temp, LET, is just above cloud base) do not allow a
!...downdraft...
!
        IF((P0(KSTART)-P0(LFS)).GT.50.E2)THEN
          THETED(LFS) = THETEE(LFS)
          QD(LFS) = Q0(LFS)
!
!...call tpmix2dd to find wet-bulb temp, qv...
!
          CALL tpmix2dd (p0(lfs),theted(lfs),tz(lfs),qss)
          THTAD(LFS)=TZ(LFS)*(P00/P0(LFS))**(0.2854*(1.-0.28*QSS))
!     
!...TAKE A FIRST GUESS AT THE INITIAL DOWNDRAFT MASS FLUX...
!     
          TVD(LFS)=TZ(LFS)*(1.+0.608*QSS)
          RDD=P0(LFS)/(R*TVD(LFS))
          A1=(1.-PEFF)*AU0
          DMF(LFS)=-A1*RDD
          DER(LFS)=DMF(LFS)
          DDR(LFS)=0.
          RHBAR = RH(LFS)*DP(LFS)
          DPTT = DP(LFS)
          DO ND = LFS-1,KSTART,-1
            ND1 = ND+1
            DER(ND)=DER(LFS)*EMS(ND)/EMS(LFS)
            DDR(ND)=0.
            DMF(ND)=DMF(ND1)+DER(ND)
            THETED(ND)=(THETED(ND1)*DMF(ND1)+THETEE(ND)*DER(ND))/DMF(ND)
            QD(ND)=(QD(ND1)*DMF(ND1)+Q0(ND)*DER(ND))/DMF(ND)    
            DPTT = DPTT+DP(ND)
            RHBAR = RHBAR+RH(ND)*DP(ND)
          ENDDO
          RHBAR = RHBAR/DPTT
          DMFFRC = 2.*(1.-RHBAR)
          DPDD = 0.
!...Calculate melting effect
!... first, compute total frozen precipitation generated...
!
          pptmlt = 0.
          DO NK = KLCL,LTOP
            PPTMLT = PPTMLT+PPTICE(NK)
          ENDDO
          if(lc.lt.ml)then
!            DTMELT = RLF*PPTMLT/(CP*DMFFRC*UMF(KLCL))
!...
!...For now, calculate melting effect as if DMF = -UMF at KLCL, i.e., as
!...if DMFFRC=1.  Otherwise, for small DMFFRC, DTMELT gets too large!
!...12/14/98 jsk...
            DTMELT = RLF*PPTMLT/(CP*UMF(KLCL))
          else
            DTMELT = 0.
          endif
          LDT = MIN0(LFS-1,KSTART-1)
          CALL tpmix2dd (p0(kstart),theted(kstart),tz(kstart),qss)
          tz(kstart) = tz(kstart)-dtmelt
          ES=ALIQ*EXP((BLIQ*TZ(KSTART)-CLIQ)/(TZ(KSTART)-DLIQ))
          QSS=0.622*ES/(P0(KSTART)-ES)
          THETED(KSTART)=TZ(KSTART)*(1.E5/P0(KSTART))**(0.2854*(1.-0.28*QSS))*    &
                EXP((3374.6525/TZ(KSTART)-2.5403)*QSS*(1.+0.81*QSS))
!....  
          LDT = MIN0(LFS-1,KSTART-1)
          DO ND = LDT,1,-1
            DPDD = DPDD+DP(ND)
            THETED(ND) = THETED(KSTART)
            QD(ND)     = QD(KSTART)       
!
!...call tpmix2dd to find wet bulb temp, saturation mixing ratio...
!
            CALL tpmix2dd (p0(nd),theted(nd),tz(nd),qss)
            qsd(nd) = qss
!
!...specify RH decrease of 20%/km in downdraft...
!
            RHH = 1.-0.2/1000.*(Z0(KSTART)-Z0(ND))
!
!...adjust downdraft TEMP, Q to specified RH:
!
            IF(RHH.LT.1.)THEN
              DSSDT=(CLIQ-BLIQ*DLIQ)/((TZ(ND)-DLIQ)*(TZ(ND)-DLIQ))
              RL=XLV0-XLV1*TZ(ND)
              DTMP=RL*QSS*(1.-RHH)/(CP+RL*RHH*QSS*DSSDT)
              T1RH=TZ(ND)+DTMP
              ES=RHH*ALIQ*EXP((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
              QSRH=0.622*ES/(P0(ND)-ES)
!
!...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
!...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...
!
              IF(QSRH.LT.QD(ND))THEN
                QSRH=QD(ND)
                T1RH=TZ(ND)+(QSS-QSRH)*RL/CP
              ENDIF
              TZ(ND)=T1RH
              QSS=QSRH
              QSD(ND) = QSS
            ENDIF         
            TVD(nd) = tz(nd)*(1.+0.608*qsd(nd))
            IF(TVD(ND).GT.TV0(ND).OR.ND.EQ.1)THEN
              LDB=ND
              EXIT
            ENDIF
          ENDDO
          IF((P0(LDB)-P0(LFS)) .gt. 50.E2)THEN   ! minimum Downdraft depth! 
            DO ND=LDT,LDB,-1
              ND1 = ND+1
              DDR(ND) = -DMF(KSTART)*DP(ND)/DPDD
              DER(ND) = 0.
              DMF(ND) = DMF(ND1)+DDR(ND)
              TDER=TDER+(QSD(nd)-QD(ND))*DDR(ND)
              QD(ND)=QSD(nd)
              THTAD(ND)=TZ(ND)*(P00/P0(ND))**(0.2854*(1.-0.28*QD(ND)))
            ENDDO
          ENDIF
        ENDIF
      ENDIF devap
!     
!...IF DOWNDRAFT DOES NOT EVAPORATE ANY WATER FOR SPECIFIED RELATIVE
!...HUMIDITY, NO DOWNDRAFT IS ALLOWED...
!     
d_mf:   IF(TDER.LT.1.)THEN
          PPTFLX=TRPPT
          CPR=TRPPT
          TDER=0.
          LDB=LFS
          DO NDK=1,LTOP
            DMF(NDK)=0.
            DER(NDK)=0.
            DDR(NDK)=0.
            THTAD(NDK)=0.
            WD(NDK)=0.
            TZ(NDK)=0.
            QD(NDK)=0.
          ENDDO
          AINCM2=100.
        ELSE 
          DDINC = -DMFFRC*UMF(KLCL)/DMF(KSTART)
!         DDINC = 0.
          IF(TDER*DDINC.GT.TRPPT)THEN
            DDINC = TRPPT/TDER
          ENDIF
          TDER = TDER*DDINC
    
          DO NK=LDB,LFS
            DMF(NK)=DMF(NK)*DDINC
            DER(NK)=DER(NK)*DDINC
            DDR(NK)=DDR(NK)*DDINC
          ENDDO
         CPR=TRPPT
         PPTFLX = TRPPT-TDER
         PEFF=PPTFLX/TRPPT
         IF(IPRNT)THEN
         ENDIF
!
!...ZERO OUT THE ARRAYS FOR DOWNDRAFT DATA AT LEVELS ABOVE AND BELOW THE
!...DOWNDRAFT...
!     
         IF(LDB.GT.1)THEN
           DO NK=1,LDB-1
             DMF(NK)=0.
             DER(NK)=0.
             DDR(NK)=0.
             WD(NK)=0.
             TZ(NK)=0.
             QD(NK)=0.
             THTAD(NK)=0.
           ENDDO
         ENDIF
         DO NK=LFS+1,KX
           DMF(NK)=0.
           DER(NK)=0.
           DDR(NK)=0.
           WD(NK)=0.
           TZ(NK)=0.
           QD(NK)=0.
           THTAD(NK)=0.
         ENDDO
         DO NK=LDT+1,LFS-1
           TZ(NK)=0.
           QD(NK)=0.
           THTAD(NK)=0.
         ENDDO
       ENDIF d_mf
!
!...SET LIMITS ON THE UPDRAFT AND DOWNDRAFT MASS FLUXES SO THAT THE INFL
!   INTO CONVECTIVE DRAFTS FROM A GIVEN LAYER IS NO MORE THAN IS AVAILAB
!   IN THAT LAYER INITIALLY...
!     
       AINCMX=1000.
       LMAX=MAX0(KLCL,LFS)
       DO NK=LC,LMAX
         IF((UER(NK)-DER(NK)).GT.0.)AINCM1=EMS(NK)/((UER(NK)-DER(NK))*TIMEC)
          AINCMX=AMIN1(AINCMX,AINCM1)
       ENDDO
       AINC=1.
       IF(AINCMX.LT.AINC)AINC=AINCMX
!     
!...SAVE THE RELEVENT VARIABLES FOR A UNIT UPDRAFT AND DOWNDRAFT...THEY WILL 
!...BE ITERATIVELY ADJUSTED BY THE FACTOR AINC TO SATISFY THE STABILIZATION
!...CLOSURE...
!     
       TDER2=TDER
       PPTFL2=PPTFLX
       DO NK=1,LTOP
         DETLQ2(NK)=DETLQ(NK)
         DETIC2(NK)=DETIC(NK)
         UDR2(NK)=UDR(NK)
         UER2(NK)=UER(NK)
         DDR2(NK)=DDR(NK)
         DER2(NK)=DER(NK)
         UMF2(NK)=UMF(NK)
         DMF2(NK)=DMF(NK)
       ENDDO


       FABE=1.
       STAB=0.95
       NOITR=0
       ISTOP=0

!
        IF(ISHALL.EQ.1)THEN                              ! First for shallow convection
!
! No iteration for shallow convection; if turbulent kinetic energy (TKE) is available
! from a turbulence parameterization, scale cloud-base updraft mass flux as a func
! of TKE, but for now, just specify shallow-cloud mass flux using TKEMAX = 5...
!
!...find the maximum TKE value between LC and KLCL...
!         TKEMAX = 0.
          TKEMAX = 5.
          IF(IDIFFK.eq.1) then
            DO 173 K = LC,KLCL
              NK = KX-K+1
              TKEMAX = AMAX1(TKEMAX,TKE(NK))
 173        CONTINUE
            TKEMAX = AMIN1(TKEMAX,10.)
            TKEMAX = AMAX1(TKEMAX,5.)
          ENDIF
!c         TKEMAX = 10.
!c...3_24_99...DPMIN was changed for shallow convection so that it is the
!c...          the same as for deep convection (5.E3).  Since this doubles
!c...          (roughly) the value of DPTHMX, add a factor of 0.5 to calcu-
!c...          lation of EVAC...
!c         EVAC  = TKEMAX*0.1
           EVAC  = 0.5*TKEMAX*0.1
! CLC test
!          EVAC=0
!         AINC = 0.1*DPTHMX*DXIJ*DXIJ/(VMFLCL*G*TIMEC)
!          AINC = EVAC*DPTHMX*DX(I,J)*DX(I,J)/(VMFLCL*G*TIMEC)
          AINC = EVAC*DPTHMX*DXSQ/(VMFLCL*G*TIMEC)
          TDER=TDER2*AINC
          PPTFLX=PPTFL2*AINC

          DO NK=1,LTOP
            UMF(NK)=UMF2(NK)*AINC
            DMF(NK)=DMF2(NK)*AINC
            DETLQ(NK)=DETLQ2(NK)*AINC
            DETIC(NK)=DETIC2(NK)*AINC
            UDR(NK)=UDR2(NK)*AINC
            UER(NK)=UER2(NK)*AINC
            DER(NK)=DER2(NK)*AINC
            DDR(NK)=DDR2(NK)*AINC
          ENDDO
        ENDIF                                           ! Otherwise for deep convection

! use iterative procedure to find mass fluxes...
iter:     DO NCOUNT=1,10
!     
!*****************************************************************
!                                                                *
!           COMPUTE PROPERTIES FOR COMPENSATIONAL SUBSIDENCE     *
!                                                                *
!*****************************************************************
!     
!...DETERMINE OMEGA VALUE NECESSARY AT TOP AND BOTTOM OF EACH LAYER TO
!...SATISFY MASS CONTINUITY...
!     
            DTT=TIMEC
    
            DO NK=1,LTOP

              DOMGDP(NK)=-(UER(NK)-DER(NK)-UDR(NK)-DDR(NK))*EMSD(NK)

              IF(NK.GT.1)THEN
                OMG(NK)=OMG(NK-1)-DP(NK-1)*DOMGDP(NK-1)
                DTT1=0.75*DP(NK-1)/(ABS(OMG(NK))+1.E-10)
                DTT=AMIN1(DTT,DTT1)
              ENDIF

            ENDDO

            DO NK=1,LTOP
              THPA(NK)=THTA0(NK)
              QPA(NK)=Q0(NK)
!              NSTEP=NINT(TIMEC/DTT+1)
!              DTIME=TIMEC/FLOAT(NSTEP)
              FXM(NK)=OMG(NK)*DXSQ/G
            ENDDO

              NSTEP=NINT(TIMEC/DTT+1)
              DTIME=TIMEC/FLOAT(NSTEP)
!     
!...DO AN UPSTREAM/FORWARD-IN-TIME ADVECTION OF THETA, QV...
!     
        DO NTC=1,NSTEP
!     
!...ASSIGN THETA AND Q VALUES AT THE TOP AND BOTTOM OF EACH LAYER BASED
!...SIGN OF OMEGA...
!     
            DO  NK=1,LTOP
              THFXIN(NK)=0.
              THFXOUT(NK)=0.
              QFXIN(NK)=0.
              QFXOUT(NK)=0.
            ENDDO
            DO NK=2,LTOP
              IF(OMG(NK).LE.0.)THEN
                THFXIN(NK)=-FXM(NK)*THPA(NK-1)
                QFXIN(NK)=-FXM(NK)*QPA(NK-1)
                THFXOUT(NK-1)=THFXOUT(NK-1)+THFXIN(NK)
                QFXOUT(NK-1)=QFXOUT(NK-1)+QFXIN(NK)
              ELSE
                THFXOUT(NK)=FXM(NK)*THPA(NK)
                QFXOUT(NK)=FXM(NK)*QPA(NK)
                THFXIN(NK-1)=THFXIN(NK-1)+THFXOUT(NK)
                QFXIN(NK-1)=QFXIN(NK-1)+QFXOUT(NK)
              ENDIF
            ENDDO
!     
!...UPDATE THE THETA AND QV VALUES AT EACH LEVEL...
!     
            DO NK=1,LTOP
              THPA(NK)=THPA(NK)+(THFXIN(NK)+UDR(NK)*THTAU(NK)+DDR(NK)*      &
                       THTAD(NK)-THFXOUT(NK)-(UER(NK)-DER(NK))*THTA0(NK))*  &
                       DTIME*EMSD(NK)
              QPA(NK)=QPA(NK)+(QFXIN(NK)+UDR(NK)*QDT(NK)+DDR(NK)*QD(NK)-    &
                      QFXOUT(NK)-(UER(NK)-DER(NK))*Q0(NK))*DTIME*EMSD(NK)
            ENDDO   
          ENDDO   
          DO NK=1,LTOP
            THTAG(NK)=THPA(NK)
            QG(NK)=QPA(NK)
          ENDDO
!     
!...CHECK TO SEE IF MIXING RATIO DIPS BELOW ZERO ANYWHERE;  IF SO, BORRO
!...MOISTURE FROM ADJACENT LAYERS TO BRING IT BACK UP ABOVE ZERO...
!     
        DO NK=1,LTOP
          IF(QG(NK).LT.0.)THEN
            IF(NK.EQ.1)THEN                             ! JSK MODS
              PRINT *,'!!!!! PROBLEM WITH KF SCHEME:  ' ! JSK MODS
              PRINT *,'QG = 0 AT THE SURFACE!!!!!!!'    ! JSK MODS
              STOP 'QG'                                 ! JSK MODS
            ENDIF                                       ! JSK MODS
            NK1=NK+1
            IF(NK.EQ.LTOP)THEN
              NK1=KLCL
            ENDIF
            TMA=QG(NK1)*EMS(NK1)
            TMB=QG(NK-1)*EMS(NK-1)
            TMM=(QG(NK)-1.E-9)*EMS(NK  )
            BCOEFF=-TMM/((TMA*TMA)/TMB+TMB)
            ACOEFF=BCOEFF*TMA/TMB
            TMB=TMB*(1.-BCOEFF)
            TMA=TMA*(1.-ACOEFF)
            IF(NK.EQ.LTOP)THEN
              QVDIFF=(QG(NK1)-TMA*EMSD(NK1))*100./QG(NK1)
              IF(ABS(QVDIFF).GT.1.)THEN
             PRINT *,'!!!WARNING!!! CLOUD BASE WATER VAPOR CHANGES BY ',     &
                      QVDIFF,                                                &
                     '% WHEN MOISTURE IS BORROWED TO PREVENT NEGATIVE ',     &
                     'VALUES IN KAIN-FRITSCH'
              ENDIF
            ENDIF
            QG(NK)=1.E-9
            QG(NK1)=TMA*EMSD(NK1)
            QG(NK-1)=TMB*EMSD(NK-1)
          ENDIF
        ENDDO
        TOPOMG=(UDR(LTOP)-UER(LTOP))*DP(LTOP)*EMSD(LTOP)

        IF(ABS(TOPOMG-OMG(LTOP)).GT.1.E-3)THEN
! abb Checking ISTOP =1
       write(6,*)'I, J, =',i,j
  WRITE(6,*)'ERROR:  MASS DOES NOT BALANCE IN KF SCHEME; TOPOMG, OMG =',TOPOMG,OMG(LTOP)
          ISTOP=1
          IPRNT=.TRUE.
          EXIT iter
        ENDIF

!     
!...CONVERT THETA TO T...
!     
        DO NK=1,LTOP
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*QG(NK)))
          TG(NK)=THTAG(NK)/EXN(NK)
          TVG(NK)=TG(NK)*(1.+0.608*QG(NK))
        ENDDO
        IF(ISHALL.EQ.1)THEN
          EXIT iter
        ENDIF
!     
!*******************************************************************
!                                                                  *
!     COMPUTE NEW CLOUD AND CHANGE IN AVAILABLE BUOYANT ENERGY.    *
!                                                                  *
!*******************************************************************
!     
!...THE FOLLOWING COMPUTATIONS ARE SIMILAR TO THAT FOR UPDRAFT
!     
!        THMIX=0.
          TMIX=0.
          QMIX=0.
!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
          DO NK=LC,KPBL
            TMIX=TMIX+DP(NK)*TG(NK)
            QMIX=QMIX+DP(NK)*QG(NK)  
          ENDDO
          TMIX=TMIX/DPTHMX
          QMIX=QMIX/DPTHMX
          ES=ALIQ*EXP((TMIX*BLIQ-CLIQ)/(TMIX-DLIQ))
          QSS=0.622*ES/(PMIX-ES)
!     
!...REMOVE SUPERSATURATION FOR DIAGNOSTIC PURPOSES, IF NECESSARY...
!     
          IF(QMIX.GT.QSS)THEN
            RL=XLV0-XLV1*TMIX
            CPM=CP*(1.+0.887*QMIX)
            DSSDT=QSS*(CLIQ-BLIQ*DLIQ)/((TMIX-DLIQ)*(TMIX-DLIQ))
            DQ=(QMIX-QSS)/(1.+RL*DSSDT/CPM)
            TMIX=TMIX+RL/CP*DQ
            QMIX=QMIX-DQ
            TLCL=TMIX
          ELSE
            QMIX=AMAX1(QMIX,0.)
            EMIX=QMIX*PMIX/(0.622+QMIX)
            astrt=1.e-3
            binc=0.075
            a1=emix/aliq
            tp=(a1-astrt)/binc
            indlu=int(tp)+1
            value=(indlu-1)*binc+astrt
            aintrp=(a1-value)/binc
            tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
            TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
            TLCL=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(TMIX-T00))*(TMIX-TDPT)
            TLCL=AMIN1(TLCL,TMIX)
          ENDIF
          TVLCL=TLCL*(1.+0.608*QMIX)
          ZLCL = ZMIX+(TLCL-TMIX)/GDRY
          DO NK = LC,KL
            KLCL=NK
            IF(ZLCL.LE.Z0(NK))THEN
              EXIT 
            ENDIF
          ENDDO
          K=KLCL-1
          DLP=(ZLCL-Z0(K))/(Z0(KLCL)-Z0(K))
!     
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!     
          TENV=TG(K)+(TG(KLCL)-TG(K))*DLP
          QENV=QG(K)+(QG(KLCL)-QG(K))*DLP
          TVEN=TENV*(1.+0.608*QENV)
          PLCL=P0(K)+(P0(KLCL)-P0(K))*DLP
          THETEU(K)=TMIX*(1.E5/PMIX)**(0.2854*(1.-0.28*QMIX))*             &
                  EXP((3374.6525/TLCL-2.5403)*QMIX*(1.+0.81*QMIX))
!     
!...COMPUTE ADJUSTED ABE(ABEG).
!     
          ABEG=0.
          DO NK=K,LTOPM1
            NK1=NK+1
            THETEU(NK1) = THETEU(NK)
            CALL tpmix2dd (p0(nk1),theteu(nk1),tgu(nk1),qgu(nk1))
            TVQU(NK1)=TGU(NK1)*(1.+0.608*QGU(NK1)-QLIQ(NK1)-QICE(NK1))
            IF(NK.EQ.K)THEN
              DZZ=Z0(KLCL)-ZLCL
              DILBE=((TVLCL+TVQU(NK1))/(TVEN+TVG(NK1))-1.)*DZZ
            ELSE
              DZZ=DZA(NK)
              DILBE=((TVQU(NK)+TVQU(NK1))/(TVG(NK)+TVG(NK1))-1.)*DZZ
            ENDIF
            IF(DILBE.GT.0.)ABEG=ABEG+DILBE*G
!
!...DILUTE BY ENTRAINMENT BY THE RATE AS ORIGINAL UPDRAFT...
!
    CALL envirtht (P0(NK1),TG(NK1),QG(NK1),THTEEG(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
    THETEU(NK1)=THETEU(NK1)*DDILFRC(NK1)+THTEEG(NK1)*(1.-DDILFRC(NK1))
          ENDDO
!     
!...ASSUME AT LEAST 90% OF CAPE (ABE) IS REMOVED BY CONVECTION DURING
!...THE PERIOD TIMEC...
!     
          IF(NOITR.EQ.1)THEN
            EXIT iter
          ENDIF
          DABE=AMAX1(ABE-ABEG,0.1*ABE)
          FABE=ABEG/(ABE+1.E-8)
          IF(FABE.GT.1. .and. ISHALL.EQ.0)THEN
!          WRITE(98,*)'UPDRAFT/DOWNDRAFT COUPLET INCREASES CAPE AT THIS
!     *GRID POINT; NO CONVECTION ALLOWED!'
            RETURN  
          ENDIF
          IF(NCOUNT.NE.1)THEN
            DFDA=(FABE-FABEOLD)/(AINC-AINCOLD)
            IF(DFDA.GT.0.)THEN
              NOITR=1
              AINC=AINCOLD
              CYCLE iter
            ENDIF
          ENDIF
          AINCOLD=AINC
          FABEOLD=FABE
          IF(AINC/AINCMX.GT.0.999.AND.FABE.GT.1.05-STAB)THEN
            EXIT
          ENDIF
          IF((FABE.LE.1.05-STAB.AND.FABE.GE.0.95-STAB) .or. NCOUNT.EQ.10)THEN
            EXIT iter
          ELSE
            IF(NCOUNT.GT.10)THEN
              EXIT
            ENDIF
!     
!...IF MORE THAN 10% OF THE ORIGINAL CAPE REMAINS, INCREASE THE CONVECTI
!...MASS FLUX BY THE FACTOR AINC:
!     
            IF(FABE.EQ.0.)THEN
              AINC=AINC*0.5
            ELSE
              AINC=AINC*STAB*ABE/(DABE+1.E-8)
            ENDIF
!           AINC=AMIN1(AINCMX,AINC)
            AINC=AMIN1(AINCMX,AINC)
!...IF AINC BECOMES VERY SMALL, EFFECTS OF CONVECTION ! JSK MODS
!...WILL BE MINIMAL SO JUST IGNORE IT...              ! JSK MODS
            IF(AINC.LT.0.05)then
              RETURN                          ! JSK MODS
            ENDIF
!            AINC=AMAX1(AINC,0.05)                        ! JSK MODS
            TDER=TDER2*AINC
            PPTFLX=PPTFL2*AINC
!           IF (XTIME.LT.10.)THEN
!           WRITE(98,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,
!          *              FABEOLD,AINCOLD 
!           ENDIF
            DO NK=1,LTOP
              UMF(NK)=UMF2(NK)*AINC
              DMF(NK)=DMF2(NK)*AINC
              DETLQ(NK)=DETLQ2(NK)*AINC
              DETIC(NK)=DETIC2(NK)*AINC
              UDR(NK)=UDR2(NK)*AINC
              UER(NK)=UER2(NK)*AINC
              DER(NK)=DER2(NK)*AINC
              DDR(NK)=DDR2(NK)*AINC
            ENDDO
!     
!...GO BACK UP FOR ANOTHER ITERATION...
!     
          ENDIF
        ENDDO iter
!     
!...COMPUTE HYDROMETEOR TENDENCIES AS IS DONE FOR T, QV...
!     
!...FRC2 IS THE FRACTION OF TOTAL CONDENSATE      !  PPT FB MODS
!...GENERATED THAT GOES INTO PRECIPITIATION       !  PPT FB MODS
!
!  Redistribute hydormeteors according to the final mass-flux values:
!
        IF(CPR.GT.0.)THEN 
          FRC2=PPTFLX/(CPR*AINC)                    !  PPT FB MODS
        ELSE
           FRC2=0.
        ENDIF
        DO NK=1,LTOP
          QLPA(NK)=QL0(NK)
          QIPA(NK)=QI0(NK)
          QRPA(NK)=QR0(NK)
          QSPA(NK)=QS0(NK)
          RAINFB(NK)=PPTLIQ(NK)*AINC*FBFRC*FRC2   !  PPT FB MODS
          SNOWFB(NK)=PPTICE(NK)*AINC*FBFRC*FRC2   !  PPT FB MODS
        ENDDO
        DO NTC=1,NSTEP
!     
!...ASSIGN HYDROMETEORS CONCENTRATIONS AT THE TOP AND BOTTOM OF EACH LAY
!...BASED ON THE SIGN OF OMEGA...
!     
          DO NK=1,LTOP
            QLFXIN(NK)=0.
            QLFXOUT(NK)=0.
            QIFXIN(NK)=0.
            QIFXOUT(NK)=0.
            QRFXIN(NK)=0.
            QRFXOUT(NK)=0.
            QSFXIN(NK)=0.
            QSFXOUT(NK)=0.
          ENDDO   
          DO NK=2,LTOP
            IF(OMG(NK).LE.0.)THEN
              QLFXIN(NK)=-FXM(NK)*QLPA(NK-1)
              QIFXIN(NK)=-FXM(NK)*QIPA(NK-1)
              QRFXIN(NK)=-FXM(NK)*QRPA(NK-1)
              QSFXIN(NK)=-FXM(NK)*QSPA(NK-1)
              QLFXOUT(NK-1)=QLFXOUT(NK-1)+QLFXIN(NK)
              QIFXOUT(NK-1)=QIFXOUT(NK-1)+QIFXIN(NK)
              QRFXOUT(NK-1)=QRFXOUT(NK-1)+QRFXIN(NK)
              QSFXOUT(NK-1)=QSFXOUT(NK-1)+QSFXIN(NK)
            ELSE
              QLFXOUT(NK)=FXM(NK)*QLPA(NK)
              QIFXOUT(NK)=FXM(NK)*QIPA(NK)
              QRFXOUT(NK)=FXM(NK)*QRPA(NK)
              QSFXOUT(NK)=FXM(NK)*QSPA(NK)
              QLFXIN(NK-1)=QLFXIN(NK-1)+QLFXOUT(NK)
              QIFXIN(NK-1)=QIFXIN(NK-1)+QIFXOUT(NK)
              QRFXIN(NK-1)=QRFXIN(NK-1)+QRFXOUT(NK)
              QSFXIN(NK-1)=QSFXIN(NK-1)+QSFXOUT(NK)
            ENDIF
          ENDDO   
!     
!...UPDATE THE HYDROMETEOR CONCENTRATION VALUES AT EACH LEVEL...
!     
          DO NK=1,LTOP
            QLPA(NK)=QLPA(NK)+(QLFXIN(NK)+DETLQ(NK)-QLFXOUT(NK))*DTIME*EMSD(NK)
            QIPA(NK)=QIPA(NK)+(QIFXIN(NK)+DETIC(NK)-QIFXOUT(NK))*DTIME*EMSD(NK)
            QRPA(NK)=QRPA(NK)+(QRFXIN(NK)-QRFXOUT(NK)+RAINFB(NK))*DTIME*EMSD(NK)         !  PPT FB MODS
            QSPA(NK)=QSPA(NK)+(QSFXIN(NK)-QSFXOUT(NK)+SNOWFB(NK))*DTIME*EMSD(NK)         !  PPT FB MODS
          ENDDO     
        ENDDO
        DO NK=1,LTOP
          QLG(NK)=QLPA(NK)
          QIG(NK)=QIPA(NK)
          QRG(NK)=QRPA(NK)
          QSG(NK)=QSPA(NK)
        ENDDO   
!
!...CLEAN THINGS UP, CALCULATE CONVECTIVE FEEDBACK TENDENCIES FOR THIS
!...GRID POINT...
!     
       IF(IPRNT)THEN  
       endif  
!     
!...SEND FINAL PARAMETERIZED VALUES TO OUTPUT FILES...
!     
           abb=ltop
	   
           DO NKj=1,LTOP
!            abb=ltop
             K=abb-NKj+1
             DTT=(TG(K)-T0(K))*86400./TIMEC
             RL=XLV0-XLV1*TG(K)
             DR=-(QG(K)-Q0(K))*RL*86400./(TIMEC*CP)
             UDFRC=UDR(K)*TIMEC*EMSD(K)
             UEFRC=UER(K)*TIMEC*EMSD(K)
             DDFRC=DDR(K)*TIMEC*EMSD(K)
             DEFRC=-DER(K)*TIMEC*EMSD(K)

! abb Reassignation of ltop (or abb) to abb2
             abb2(k)=abb
           end do
            
           DO NK=1,KL
             K=KX-NK+1
             DTT=TG(K)-T0(K)
             TUC=TU(K)-T00
             IF(K.LT.LC.OR.K.GT.LTOP)TUC=0.
             TDC=TZ(K)-T00
             IF((K.LT.LDB.OR.K.GT.LDT).AND.K.NE.LFS)TDC=0.
             IF(T0(K).LT.T00)THEN
               ES=ALIQ*EXP((BLIQ*TG(K)-CLIQ)/(TG(K)-DLIQ))
             ELSE
               ES=ALIQ*EXP((BLIQ*TG(K)-CLIQ)/(TG(K)-DLIQ))
             ENDIF  
             QGS=ES*0.622/(P0(K)-ES)
             RH0=Q0(K)/QES(K)
             RHG=QG(K)/QGS
           ENDDO
!     
!...IF CALCULATIONS ABOVE SHOW AN ERROR IN THE MASS BUDGET, PRINT OUT A
!...TO BE USED LATER FOR DIAGNOSTIC PURPOSES, THEN ABORT RUN...
!    

         IF(ISTOP.EQ.1.or.ISHALL.EQ.1)THEN
         IF(ISHALL.NE.1)THEN
            IF(ISTOP.EQ.1)THEN
              print 90, istop
90            format(' stop in KAIN-FRITSCH because an error in the ' &
                     ,'mass budget value of istop is 1 ', I2)
              STOP 'KAIN-FRITSCH'
            ENDIF
         ENDIF
       ENDIF

!     
!  EVALUATE MOISTURE BUDGET...
!     

        QINIT=0.
        QFNL=0.
        DPT=0.

        do nk=1,ltop
          dpt=dpt+dp(nk)
          qinit=qinit+q0(nk)*ems(nk)
          qinit1=qinit
          qinit=qinit+(ql0(nk)+qi0(nk)+qr0(nk)+qs0(nk))*ems(nk)
          qinit2=qinit
          qfnl=qfnl+qg(nk)*ems(nk)
          qfnl1=qfnl
          qfnl=qfnl+(qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
          qfnl2=qfnl
        end do

        QFNL=QFNL+PPTFLX*TIMEC*(1.-FBFRC)       !  PPT FB MODS
!        QFNL=QFNL+PPTFLX*TIMEC                 !  PPT FB MODS
        ERR2=(QFNL-QINIT)*100./QINIT
      
      IF(ABS(ERR2).GT.0.05 .AND. ISTOP.EQ.0)THEN 
        IPRNT=.TRUE.
        ISTOP=1
         print*, 'QFNL=', QFNL
         print*, 'QINIT=', QINIT
         print*, 'QINIT1=', QINIT1
         print*, 'QINIT2 msg=', QINIT2
         print*, 'QFNL1=', QFNL1
         print*, 'QFNL2=', QFNL2
         STOP 'QVERR'
      ENDIF
!        RELERR=ERR2*QINIT/(PPTFLX*TIMEC+1.E-10)
!     
!...FEEDBACK TO RESOLVABLE SCALE TENDENCIES.
!     
!...IF THE ADVECTIVE TIME PERIOD (TADVEC) IS LESS THAN SPECIFIED MINIMUM
!...TIMEC, ALLOW FEEDBACK TO OCCUR ONLY DURING TADVEC...
!     
        IF(TADVEC.LT.TIMEC)NIC=NINT(TADVEC/DT)
        NCA(I,J)=NIC
        IF(ISHALL.EQ.1)THEN
          TIMEC = 2400.
          NCA(I,J) = NTST
        ENDIF 
        DO K=1,KX
!         IF(IMOIST(INEST).NE.2)THEN
!
!...IF HYDROMETEORS ARE NOT ALLOWED, THEY MUST BE EVAPORATED OR SUBLIMAT
!...AND FED BACK AS VAPOR, ALONG WITH ASSOCIATED CHANGES IN TEMPERATURE.
!...NOTE:  THIS WILL INTRODUCE CHANGES IN THE CONVECTIVE TEMPERATURE AND
!...WATER VAPOR FEEDBACK TENDENCIES AND MAY LEAD TO SUPERSATURATED VALUE
!...OF QG...
!
!   FOR RAMS ACTIVATE THIS OPTION FOR NO MICROPHYSICS OR DUMPBUCKET

          IF(IMPHYS.LT.2.OR.ICLOUD.EQ.0.OR.IRAIN.EQ.0)THEN
           RLC=XLV0-XLV1*TG(K)
           RLS=XLS0-XLS1*TG(K)
           CPM=CP*(1.+0.887*QG(K))
           TG(K)=TG(K)-(RLC*(QLG(K)+QRG(K))+RLS*(QIG(K)+QSG(K)))/CPM
           QG(K)=QG(K)+(QLG(K)+QRG(K)+QIG(K)+QSG(K))
           DQIDT(K)=0.
           DQRDT(K)=0.
           DQSDT(K)=0.
          END IF

! RAMS MICROPHYSICS INTERACTIVE-KF OPTIONS

         IF(IMPHYS.GT.1) THEN

!   LIQUID WATER PHASES ONLY ALLOWED (IPRIS and ISNOW = 0)

          IF(IPRIS.EQ.0.OR.ISNOW.EQ.0)THEN
!
!...IF ICE PHASE IS NOT ALLOWED, MELT ALL FROZEN HYDROMETEORS...
!
            CPM=CP*(1.+0.887*QG(K))
            TG(K)=TG(K)-(QIG(K)+QSG(K))*RLF/CPM
            DQCDT(K)=(QLG(K)+QIG(K)-QL0(K)-QI0(K))/TIMEC
            DQIDT(K)=0.
            DQRDT(K)=(QRG(K)+QSG(K)-QR0(K)-QS0(K))/TIMEC
            DQSDT(K)=0.

! FOR RAMS: CONSIDER ALL TENDENCIES OF HYDROMETEORS IF MICROPHYSICS
! PARAMETERS > 1

          ELSE
!
!...IF MIXED PHASE HYDROMETEORS ARE ALLOWED, FEED BACK CONVECTIVE TENDEN
!...OF HYDROMETEORS DIRECTLY...
!
            DQCDT(K)=(QLG(K)-QL0(K))/TIMEC
            DQIDT(K)=(QIG(K)-QI0(K))/TIMEC
            DQRDT(K)=(QRG(K)-QR0(K))/TIMEC
            DQSDT(K)=(QSG(K)-QS0(K))/TIMEC

          END IF          
       END IF


          DTDT(K)=(TG(K)-T0(K))/TIMEC
          DQDT(K)=(QG(K)-Q0(K))/TIMEC
        ENDDO

! comment out line below to convert to a precip tendency
!        RAINCV(I,J)=DT*PPTFLX*(1.-FBFRC)/DXSQ     !  PPT FB MODS

! Final formulation for convective rainfall (rate) 

         RAINCV(I,J)=PPTFLX*(1.-FBFRC)/DXSQ
 
!        RAINCV(I,J)=.1*.5*DT*PPTFLX/DXSQ               !  PPT FB MODS
!         RNC=0.1*TIMEC*PPTFLX/DXSQ
        RNC=RAINCV(I,J)*NIC
!      WRITE(98,*)'at NTSD =',NTSD,',No. of KF points activated =',
!     *            NCCNT

END SUBROUTINE kf_eta_para

!##############################################################################
Subroutine tpmix2 (ppa,tthes,ttu,qui,qqliq,qqice,qqnewlq,qqnewic,xxlv1,xxlv0)

implicit none

   REAL, INTENT(IN)   :: ppa,tthes,xxlv1,xxlv0
   REAL, INTENT(OUT)   :: ttu,qqnewlq,qqnewic
   REAL, INTENT(INOUT)   :: qui,qqliq,qqice
   REAL    ::    tp,qq,bth,tth,pp,t00,t10,t01,t11,q00,q10,q01,q11,          &
                 TEMP,qnew,dq,qqtot,RLL,CPP
   INTEGER ::    IPTB,ITHTB
   real :: var1, var2,a1,a2,a3

!***********************************************************************
!     scaling pressure & tt table index                         
!***********************************************************************

      var2=qui

      tp=(ppa-plutop)*rdpr
      qq=tp-aint(tp)
      iptb=int(tp)+1

!***********************************************************************
!              base and scaling factor for the                           
!***********************************************************************
!
!  scaling the & tt table index                                        
      bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
      tth=(tthes-bth)*rdthk
      pp   =tth-aint(tth)
      ithtb=int(tth)+1

      t00=ttab(ithtb,iptb)
      t10=ttab(ithtb+1,iptb)
!            var2=var2+qtOt
      t01=ttab(ithtb,iptb+1)
      t11=ttab(ithtb+1,iptb+1)

      q00=qstab(ithtb,iptb)
      q10=qstab(ithtb+1,iptb)
      q01=qstab(ithtb,iptb+1)
      q11=qstab(ithtb+1,iptb+1)

      a1=(q10-q00)*pp
      a2=(q01-q00)*qq
      a3=(q00-q10-q01+q11)*pp*qq
  
!
!***********************************************************************
!              parcel temperature                                        
!***********************************************************************
!
      temp=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
!
! abb original formulation of qs 
!      qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)

      var1=a3+a2
      var1=var1+a1
      var1=var1+q00

      if(var1.le.var2) then
        qnew=var2-var1
         var2=var1
        else

!   IF THE PARCEL IS SUBSATURATED, TEMPERATURE AND MIXING RATIO MUST BE
!   ADJUSTED...IF LIQUID WATER IS PRESENT, IT IS ALLOWED TO EVAPORATE

        qnew=0.
        dq=var1-var2
        qqtOt=qqliq+qqice

!   IF THERE IS ENOUGH LIQUID OR ICE TO SATURATE THE PARCEL, TEMP STAYS AT ITS
!   WET BULB VALUE, VAPOR MIXING RATIO IS AT SATURATED LEVEL, AND THE MIXING
!   RATIOS OF LIQUID AND ICE ARE ADJUSTED TO MAKE UP THE ORIGINAL SATURATION
!   DEFICIT... OTHERWISE, ANY AVAILABLE LIQ OR ICE VAPORIZES AND APPROPRIATE
!   ADJUSTMENTS TO PARCEL TEMP; VAPOR, LIQUID, AND ICE MIXING RATIOS ARE MADE.
!
!...subsaturated values only occur in calculations involving various mixtures of
!...updraft and environmental air for estimation of entrainment and detrainment.
!...For these purposes, assume that reasonable estimates can be given using 
!...liquid water saturation calculations only - i.e., ignore the effect of the
!...ice phase in this process only...will not affect conservative properties...

! abb In the original code the if condition was:
!        IF(qqtot.ge.DQ)THEN
! abb When changed to "gt" the code ran OK in cluster

        IF(qqtot.gt.DQ)THEN
          qqliq=qqliq-dq*qqliq/qqtot
          qqice=qqice-dq*qqice/qqtot
          var2=var1
        ELSE
          RLL=xxlv0-xxlv1*TEMP
          CPP=1005.7*(1.+0.89*var2)

          IF(qqtot.LT.1.E-10)THEN

!...IF NO LIQUID WATER OR ICE IS AVAILABLE, TEMPERATURE IS GIVEN BY:
            TEMP=TEMP+RLL*(DQ/(1.+DQ))/CPP
            ELSE

!...IF SOME LIQ WATER/ICE IS AVAILABLE, BUT NOT ENOUGH TO ACHIEVE SATURATION,
!   THE TEMPERATURE IS GIVEN BY:

            TEMP=TEMP+RLL*((DQ-qqtot)/(1+DQ-qqtot))/CPP
            var2=var2+qqtot
            qqtOt=0.
            qqliq=0.
            qqice=0.
          ENDIF
        ENDIF
      ENDIF
      ttu=TEMP
      qqnewlq=qnew
      qqnewic=0.

! abb Reassign var2
    qui=var2

END SUBROUTINE tpmix2

!##############################################################################
Subroutine dtfrznew (TU,P,THTEU,QU,QFRZ,QICE,ALIQ,BLIQ,CLIQ,DLIQ)

implicit none

   REAL,         INTENT(IN   )   :: P,QFRZ,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(INOUT)   :: TU,THTEU,QU,QICE
   REAL    ::    RLC,RLS,RLF,CPP,A,DTFRZ,ES,QS,DQEVAP,PII

!...ALLOW THE FREEZING OF LIQUID WATER IN THE UPDRAFT TO PROCEED AS AN 
!...APPROXIMATELY LINEAR FUNC OF TEMPERATURE IN THE TEMPERATURE RANGE 
!...TTFRZ TO TBFRZ...
!...FOR COLDER TERMPERATURES, FREEZE ALL LIQUID WATER...
!...THERMODYNAMIC PROPERTIES ARE STILL CALCULATED WITH RESPECT TO LIQUID WATER
!...TO ALLOW THE USE OF LOOKUP TABLE TO EXTRACT TMP FROM THETAE...

      RLC=2.5E6-2369.276*(TU-273.15)
      RLS=2833922.-259.532*(TU-273.15)
      RLF=RLS-RLC
      CPP=1005.7*(1.+0.89*QU)

!  A = D(es)/DT IS THAT CALCULATED FROM BUCK'S (1981) EMPERICAL FORMULAS
!  FOR SATURATION VAPOR PRESSURE...

      A=(CLIQ-BLIQ*DLIQ)/((TU-DLIQ)*(TU-DLIQ))
      DTFRZ = RLF*QFRZ/(CPP+RLS*QU*A)
      TU = TU+DTFRZ
      
      ES = ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
      QS = ES*0.622/(P-ES)

!...FREEZING WARMS THE AIR AND IT BECOMES UNSATURATED...ASSUME THAT SOME OF THE 
!...LIQUID WATER THAT IS AVAILABLE FOR FREEZING EVAPORATES TO MAINTAIN SATURA-
!...TION...SINCE THIS WATER HAS ALREADY BEEN TRANSFERRED TO THE ICE CATEGORY,
!...SUBTRACT IT FROM ICE CONCENTRATION, THEN SET UPDRAFT MIXING RATIO AT THE NEW
!...TEMPERATURE TO THE SATURATION VALUE...

      DQEVAP = QS-QU
      QICE = QICE-DQEVAP
      QU = QU+DQEVAP
      PII=(1.E5/P)**(0.2854*(1.-0.28*QU))
      THTEU=TU*PII*EXP((3374.6525/TU-2.5403)*QU*(1.+0.81*QU))

END SUBROUTINE dtfrznew

!##############################################################################
Subroutine condload (QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ,           &
                          QNEWIC,QLQOUT,QICOUT,G)

implicit none

!-----------------------------------------------------------------------
!  9/18/88...THIS PRECIPITATION FALLOUT SCHEME IS BASED ON THE SCHEME US
!  BY OGURA AND CHO (1973).  LIQUID WATER FALLOUT FROM A PARCEL IS CAL-
!  CULATED USING THE EQUATION DQ=-RATE*Q*DT, BUT TO SIMULATE A QUASI-
!  CONTINUOUS PROCESS, AND TO ELIMINATE A DEPENDENCY ON VERTICAL
!  RESOLUTION THIS IS EXPRESSED AS Q=Q*EXP(-RATE*DZ).

      REAL, INTENT(IN   )   :: G
      REAL, INTENT(IN   )   :: DZ,BOTERM,ENTERM,RATE
      REAL, INTENT(INOUT)   :: QLQOUT,QICOUT,WTW,QLIQ,QICE,QNEWLQ,QNEWIC
      REAL :: QTOT,QNEW,QEST,G1,WAVG,CONV,RATIO3,OLDQ,RATIO4,DQ,PPTDRG


!  9/18/88...THIS PRECIPITATION FALLOUT SCHEME IS BASED ON THE SCHEME US
!  BY OGURA AND CHO (1973).  LIQUID WATER FALLOUT FROM A PARCEL IS CAL- 
!  CULATED USING THE EQUATION DQ=-RATE*Q*DT, BUT TO SIMULATE A QUASI-   
!  CONTINUOUS PROCESS, AND TO ELIMINATE A DEPENDENCY ON VERTICAL        
!  RESOLUTION THIS IS EXPRESSED AS Q=Q*EXP(-RATE*DZ).                   

      QTOT=QLIQ+QICE                                                    
      QNEW=QNEWLQ+QNEWIC                                                

!                                                                       
!  ESTIMATE THE VERTICAL VELOCITY SO THAT AN AVERAGE VERTICAL VELOCITY 
!  BE CALCULATED TO ESTIMATE THE TIME REQUIRED FOR ASCENT BETWEEN MODEL 
!  LEVELS...                                                            
!                                                                       

      QEST=0.5*(QTOT+QNEW)                                              
      G1=WTW+BOTERM-ENTERM-2.*G*DZ*QEST/1.5                             
      IF(G1.LT.0.0)G1=0.                                                
      WAVG=0.5*(SQRT(WTW)+SQRT(G1))                                      
      CONV=RATE*DZ/WAVG                                                 
!                                                                       
!  RATIO3 IS THE FRACTION OF LIQUID WATER IN FRESH CONDENSATE, RATIO4 IS
!  THE FRACTION OF LIQUID WATER IN THE TOTAL AMOUNT OF CONDENSATE INVOLV
!  IN THE PRECIPITATION PROCESS - NOTE THAT ONLY 60% OF THE FRESH CONDEN
!  SATE IS IS ALLOWED TO PARTICIPATE IN THE CONVERSION PROCESS...       
!                                                                       
      RATIO3=QNEWLQ/(QNEW+1.E-10)                                       
!     OLDQ=QTOT                                                         
      QTOT=QTOT+0.6*QNEW                                                
      OLDQ=QTOT                                                         
      RATIO4=(0.6*QNEWLQ+QLIQ)/(QTOT+1.E-10)                            
      QTOT=QTOT*EXP(-CONV)                                              
!                                                                       
!  DETERMINE THE AMOUNT OF PRECIPITATION THAT FALLS OUT OF THE UPDRAFT  
!  PARCEL AT THIS LEVEL...                                              
!                                                                       
      DQ=OLDQ-QTOT                                                      
      QLQOUT=RATIO4*DQ                                                  
      QICOUT=(1.-RATIO4)*DQ                                             
!                                                                       
!  ESTIMATE THE MEAN LOAD OF CONDENSATE ON THE UPDRAFT IN THE LAYER, CAL
!  LATE VERTICAL VELOCITY                                               
!                                                                       
      PPTDRG=0.5*(OLDQ+QTOT-0.2*QNEW)                                   
      WTW=WTW+BOTERM-ENTERM-2.*G*DZ*PPTDRG/1.5                          

!                                                                       
!  DETERMINE THE NEW LIQUID WATER AND ICE CONCENTRATIONS INCLUDING LOSSE
!  DUE TO PRECIPITATION AND GAINS FROM CONDENSATION...                  
!                                                                       
      QLIQ=RATIO4*QTOT+RATIO3*0.4*QNEW                                  
      QICE=(1.-RATIO4)*QTOT+(1.-RATIO3)*0.4*QNEW                        
      QNEWLQ=0.                                                         
      QNEWIC=0.                                                         

END SUBROUTINE condload

!##############################################################################
Subroutine prof5 (EQ,EE,UD)                                        

!***********************************************************************
!*****    GAUSSIAN TYPE MIXING PROFILE....******************************
!  THIS BROUTINE INTEGRATES THE AREA UNDER THE CURVE IN THE GAUSSIAN  
!  DISTRIBUTION...THE NUMERICAL APPROXIMATION TO THE INTEGRAL IS TAKEN F
!  "HANDBOOK OF MATHEMATICAL FUNC WITH FORMULAS, GRAPHS AND MATHEMA
!  TABLES" ED. BY ABRAMOWITZ AND STEGUN, NAT'L BUREAU OF STANDARDS APPLI
!  MATHEMATICS SERIES.  JUNE, 1964., MAY, 1968.                         
!                                     JACK KAIN                         
!                                     7/6/89                            
!-----------------------------------------------------------------------

implicit none

   REAL,         INTENT(IN   )   :: EQ
   REAL,         INTENT(INOUT)   :: EE,UD
   REAL ::       SQRT2P,A1,A2,A3,P,SIGMA,FE,X,Y,EY,E45,T1,T2,C1,C2

      DATA SQRT2P,A1,A2,A3,P,SIGMA,FE/2.506628,0.4361836,-0.1201676,       &
           0.9372980,0.33267,0.166666667,0.202765151/                        
      X=(EQ-0.5)/SIGMA                                                  
      Y=6.*EQ-3.                                                        
      EY=EXP(Y*Y/(-2))                                                  
      E45=EXP(-4.5)                                                     
      T2=1./(1.+P*ABS(Y))                                               
      T1=0.500498                                                       
      C1=A1*T1+A2*T1*T1+A3*T1*T1*T1                                     
      C2=A1*T2+A2*T2*T2+A3*T2*T2*T2                                     
      IF(Y.GE.0.)THEN                                                   
        EE=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*EQ*EQ/2.
        UD=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*(0.5+EQ*EQ/2.-    &
           EQ)                                                          
      ELSE                                                              
        EE=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*EQ*EQ/2.       
        UD=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*(0.5+EQ*   &
           EQ/2.-EQ)                                                    
      ENDIF                                                             
      EE=EE/FE                                                          
      UD=UD/FE                                                          

END SUBROUTINE prof5

!##############################################################################
Subroutine tpmix2dd (p,thes,ts,varqs)

implicit none

!-----------------------------------------------------------------------
   REAL,         INTENT(IN)   :: P,THES
   REAL,         INTENT(OUT)   :: TS,varqs
   REAL    ::    TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11
   INTEGER ::    IPTB,ITHTB
!-----------------------------------------------------------------------

!***********************************************************************
!     scaling pressure & tt table index                         
!***********************************************************************

      tp=(p-plutop)*rdpr
      qq=tp-aint(tp)
      iptb=int(tp)+1

!***********************************************************************
!              base and scaling factor for the                           
!***********************************************************************

!  scaling the & tt table index                                        
      bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
      tth=(thes-bth)*rdthk
      pp   =tth-aint(tth)
      ithtb=int(tth)+1

      t00=ttab(ithtb  ,iptb  )
      t10=ttab(ithtb+1,iptb  )
      t01=ttab(ithtb  ,iptb+1)
      t11=ttab(ithtb+1,iptb+1)

      q00=qstab(ithtb  ,iptb  )
      q10=qstab(ithtb+1,iptb  )
      q01=qstab(ithtb  ,iptb+1)
      q11=qstab(ithtb+1,iptb+1)

!     parcel temperature and saturation mixing ratio                                        
      ts=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
      varqs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)

END SUBROUTINE tpmix2dd

!##############################################################################
Subroutine envirtht (P1,T1,Q1,THT1,ALIQ,BLIQ,CLIQ,DLIQ)                       

implicit none

   REAL,         INTENT(IN   )   :: P1,T1,Q1,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(INOUT)   :: THT1
   REAL    ::    EE,TLOG,ASTRT,AINC,A1,TP,VALUE,AINTRP,TDPT,TSAT,THT,      &
                 T00,P00,C1,C2,C3,C4,C5
   INTEGER ::    INDLU

   DATA T00,P00,C1,C2,C3,C4,C5/273.15,1.E5,3374.6525,2.5403,3114.834,   &
           0.278296,1.0723E-3/                                          
                                                                     
!  CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMPERATURE...
!  NOTE: Calculations for mixed/ice phase no longer used...jsk 8/00

      EE=Q1*P1/(0.622+Q1)                                             
!     TLOG=ALOG(EE/ALIQ)                                              
! ...calculate LOG term using lookup table...

      astrt=1.e-3
      ainc=0.075
      a1=ee/aliq
      tp=(a1-astrt)/ainc
      indlu=int(tp)+1
      value=(indlu-1)*ainc+astrt
      aintrp=(a1-value)/ainc
      tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)

      TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)                               
      TSAT=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(T1-T00))*(T1-TDPT) 
      THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))                          
      THT1=THT*EXP((C1/TSAT-C2)*Q1*(1.+0.81*Q1))                      

END SUBROUTINE envirtht

!##############################################################################
Subroutine kf_lutab (SVP1,SVP2,SVP3,SVPT0)

!This routine is a lookup table.
!Given a series of series of saturation equivalent potential 
!temperatures, the temperature is calculated.

implicit none

!Lookup table variables:
!     INTEGER, SAVE, PARAMETER :: KFNT=250,KFNP=220
!     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
!     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
!     REAL, SAVE, DIMENSION(1:200) :: ALU
!     REAL, SAVE :: RDPR,RDTHK,PLUTOP
!End of Lookup table variables:

     INTEGER :: KP,IT,ITCNT,I
     REAL :: DTH,TMIN,TOLER,PBOT,DPR,                               &
             TEMP,P,ES,QS,PI,THES,TGUES,THGUES,F0,T1,T0,F1,DT, &
             ASTRT,AINC,A1,THTGS
     REAL    :: ALIQ,BLIQ,CLIQ,DLIQ,SVP1,SVP2,SVP3,SVPT0

! equivalent potential temperature increment
      data dth/1./
! minimum starting temp 
      data tmin/150./
! tolerance for accuracy of temperature 
      data toler/0.001/
! top pressure (pascals)
      plutop=100.0
! bottom pressure (pascals)
      pbot=110000.0

      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3

! compute parameters

! 1./(sat. equiv. theta increment)
      rdthk=1./dth
! pressure increment
      dpr=(pbot-plutop)/float(kfnp-1)
! 1./(pressure increment)
      rdpr=1./dpr

! compute the spread of thes
!     thespd=dth*(kfnt-1)

! calculate the starting sat. equiv. theta
      temp=tmin 
      p=plutop-dpr
      es=aliq*exp((bliq*temp-cliq)/(temp-dliq))

      do kp=1,kfnp
        p=p+dpr
!        es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
        qs=0.622*es/(p-es)
        pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
        the0k(kp)=temp*pi*exp((3374.6525/temp-2.5403)*qs*(1.+0.81*qs))
      enddo   

! compute temperatures for each sat. equiv. potential temp.

      p=plutop-dpr
      do kp=1,kfnp
        thes=the0k(kp)-dth
        p=p+dpr
        do it=1,kfnt
! define sat. equiv. pot. temp.
          thes=thes+dth
! iterate to find temperature
! find initial guess
          if(it.eq.1) then
            tgues=tmin
          else
            tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          qs=0.622*es/(p-es)
          pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
          thgues=tgues*pi*exp((3374.6525/tgues-2.5403)*qs*      &
               (1.+0.81*qs))
          f0=thgues-thes
          t1=tgues-0.5*f0
          t0=tgues
          itcnt=0
! iteration loop

          do itcnt=1,11
            es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
            qs=0.622*es/(p-es)
            pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
            thtgs=t1*pi*exp((3374.6525/t1-2.5403)*qs*(1.+0.81*qs))
            f1=thtgs-thes
            if(abs(f1).lt.toler)then
              exit
            endif
!           itcnt=itcnt+1
            dt=f1*(t1-t0)/(f1-f0)
            t0=t1
            f0=f1
            t1=t1-dt
!            if(itcnt.eq.11) then
!            endif
          enddo 

          ttab(it,kp)=t1 
          qstab(it,kp)=qs
        enddo
      enddo   

! lookup table for tlog(emix/aliq)

! set up intial values for lookup tables

       astrt=1.e-3
       ainc=0.075

       a1=astrt-ainc
       do i=1,200
         a1=a1+ainc
         alu(i)=alog(a1)
       enddo   

END SUBROUTINE kf_lutab

END MODULE cu_kfeta

