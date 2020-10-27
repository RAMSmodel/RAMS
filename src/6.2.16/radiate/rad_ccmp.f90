!##############################################################################
Subroutine lwradc (NZPP,RVR,RTR,DN0R,TEMPRD,PRD,PIRD,DZZR,FTHR,RLONG  &
                 ,SCR1,SCR2,SCR3,VPR,DMR,CO2,CLD,UF1,UF2,DF1  &
                 ,DF2,EM1,EM2,BB1,BB2)

!--------------------------------------------------------------------
!  Longwave radiation parameterization based on Rodgers and
!  Stephens and discussed in Chen and Cotton (1983).  First written
!  by Chen, later modified by Bjorn.  All variables are in cgs so
!  as to confuse people.
!
!  Upward and downward fluxes will be calculated at w points while
!  the heating rates will be at thermo points. The routine takes as
!  input the first five arrays passed in:
!  RVR ......... vapor mixing ratio
!  RTR ......... total water mixing ratio
!  DN0R ........ density in (cgs units)
!  TEMPRD ...... temperature in K (TEMPRD(1) is surface temperature)
!  PRD ......... pressure at thermo points in (cgs units)
!  PIRD   - Exner func (p/p00)**(R/Cp) at each level
!  DZTR ........ inverse distance between w points (cgs units)
!
!  and uses the following 17 entries as output or scratch space.
!
!  FTHRL ....... longwave heating rate at w point (K/s)
!  RLONG ....... net flux into ground (or lowest level)
!  BB1 ......... source func for water vapor
!  BB2 ......... source func for CO2
!  VPR ......... water vapor path
!  DMR ......... water vapor path * vapor pressure for dimer correction
!  CO2 ......... CO2 path
!  CLD ......... liquid water path
!  EM1 & 2 ..... scratch arrays for emissivities
!  SCR1 2 & 3 .. general scratch arrays
!  FU1,FU2 ..... upwelling fluxes (1-vapor) (2-CO2)
!  FD1,FD2 ..... downwelling fluxes (1-vapor) (2-CO2)
!--------------------------------------------------------------------

implicit none

INTEGER NZPP
REAL RVR(NZPP),RTR(NZPP),DN0R(NZPP),PRD(NZPP),TEMPRD(NZPP)
REAL FTHR(NZPP),DZZR(NZPP),RLONG
REAL SCR1(NZPP),SCR2(NZPP),SCR3(NZPP),PIRD(NZPP)
REAL VPR(NZPP),DMR(NZPP),CO2(NZPP),CLD(NZPP)
REAL UF1(NZPP),UF2(NZPP),DF1(NZPP),DF2(NZPP)
REAL EM1(NZPP),EM2(NZPP),BB1(NZPP),BB2(NZPP)

REAL PATH_FACT,PRES_WGHT,SIGMAT_FACT,TRANS
REAL STEFAN,X,G,SX2,SX1,FHL,RCO2,CP
INTEGER K,KK,KL,NZ,NZP,LCLDBS,LCLDTP
SAVE
!
!     Water vapor band. The vibration rotation and continuum effects
!     of the water vapor are considered.
!
REAL AD(4),BD(5),AU(4),BU(5),EU(5),ED(5)
REAL C11,C21,C31,C02,C12,C22,C32,B1,B2,BNDWID,BNDI
DATA AD/8.857,-332.8,14607.,-261900./
DATA BD/.6558,.12175,1.4976E-2,1.4981E-3,.49E-4/
DATA AU/9.329,-446.4,824,259700./
DATA BU/.5983,.15068,3.4041E-2,6.5535E-3,4.887E-4/
DATA EU/.21699,-9.185E-2,-7.971E-2,-1.502E-2,-8.754E-4/
DATA ED/.2837,-.1231,-.1057,-.0199,-1.16E-3/
DATA B1,B2/7.345,142.47/
DATA C11,C21,C31/160.87,-326.5,-158.22/
DATA C02,C12,C22,C32/74.103,19.632,0.821,-0.11834/
DATA BNDWID/200./
!
! the CO2 concentration is assumed to be 330 ppm from surface to 40km:
!          CO2(K) = 330./1000.*(44.011/28.966)*(22415./44.011)
! also downwelling flux at top of model is given by FHL
!
DATA RCO2/.25537/,FHL/0./
DATA CP/1.004E7/,STEFAN/5.6696E-5/,G/980./
!
NZP=NZPP-1
NZ=NZP-1
BNDI=1./BNDWID
!
! calculation of optical path lengths
!
DO K=2,NZ
   PRES_WGHT=(PRD(K)/1.01325E6)**.86
   PATH_FACT=PRES_WGHT*DN0R(K)/DZZR(K)
   VPR(K)= RVR(K)*PATH_FACT
   CLD(K)=(RTR(K)-RVR(K))*PATH_FACT
   DMR(K)=VPR(K)*PRD(K)*RVR(K)/(.622*1.01325E6)
   CO2(K)=RCO2*PATH_FACT
ENDDO
VPR(NZP)=VPR(NZ)
VPR(NZPP)=VPR(NZ)
CLD(NZP)=0.
CLD(NZPP)=0.
DMR(NZP)=DMR(NZ)
DMR(NZPP)=DMR(NZ)
CO2(NZP)=CO2(NZ)
CO2(NZPP)=RCO2*PRD(NZP)/G
!
! computation of black body source func with weightings given
! by sigma t factor of 0.87*(1/567-1/767)*TEMP
!
DO K=1,NZP
   SIGMAT_FACT=.0004001*TEMPRD(K)
   BB1(K)=STEFAN*TEMPRD(K)**4*(1.-SIGMAT_FACT)
   BB2(K)=BB1(K)*SIGMAT_FACT/(1.-SIGMAT_FACT)
ENDDO
BB1(NZPP)=BB1(NZP)
BB2(NZPP)=BB2(NZP)
UF1(1)=BB1(1)
UF2(1)=BB2(1)
!
! here the level of the lowest and highest cloud levels are computed
! so as to control the region for which mixed emmissivities need be
! computed
!
LCLDBS=NZPP+1
LCLDTP=0
DO K=1,NZP
   IF(RTR(K)-RVR(K).GT.0.)LCLDBS=MIN(K,LCLDBS)
   IF(RTR(K)-RVR(K).GT.0.)LCLDTP=MAX(K,LCLDTP)
ENDDO
!
! ----------------------------- Upward computations
!
EM1(1)=0.
EM2(1)=0.
DO K=2,NZP
!
! Sum the vapor, dimer corrected vapor and co2 path lengths for use
! in the emissivity polynomial fits
!
   SCR1(1)=0.
   SCR2(1)=0.
   SCR3(1)=0.
   DO KK=2,K
      SCR1(KK)=SCR1(KK-1)+VPR(K-KK+2)
      SCR2(KK)=SCR2(KK-1)+DMR(K-KK+2)
      SCR3(KK)=SCR3(KK-1)+CO2(K-KK+2)
   ENDDO
!
! Find level of path length of 1E-3 gm/cm^2 and compute upward
! emissivities for water vapor.  store in array indexed IUE1
!
   DO KK=2,K
      IF(SCR1(KK).GT.1.E-3)GO TO 23
   ENDDO
   KL=K+1
   GO TO 24
23      KL=KK
24      CONTINUE

   DO KK=2,KL-1
      X=SQRT(SCR1(KK))
      EM1(KK)=X*(AU(1)+X*(AU(2)+X*(AU(3)+X*AU(4))))
   ENDDO
   DO KK=KL,K
      X=LOG(SCR1(KK))
      EM1(KK)=BU(1)+X*(BU(2)+X*(BU(3)+X*(BU(4)+X*BU(5))))
   ENDDO
!
! Correct vapor emissivities for dimer path length
!
   DO KK=2,K
      IF(SCR2(KK).GT.1E-3)GO TO 33
   ENDDO
   GO TO 34
33      KL=KK

   DO KK=KL,K
      X=LOG(MIN(1.,SCR2(KK)))
      EM1(KK)=EM1(KK)  &
             +EU(1)+X*(EU(2)+X*(EU(3)+X*(EU(4)+X*EU(5))))
      EM1(KK)=MIN(1.,EM1(KK))
   ENDDO
34      CONTINUE
!
! Compute upward emissivities for CO2, storing in IUE2, again finding
! the level of the critical path length
!
   DO KK=2,K
      IF(SCR3(KK).GT.1.E-2)GO TO 43
   ENDDO
   KL=K+1
   GO TO 44
43      KL=KK
44      CONTINUE

   DO KK=2,KL-1
      X=SQRT(SCR3(KK))
      EM2(KK)=1.-(X*(C11+X*(C21+X*C31)))*BNDI
   ENDDO
   DO KK=KL,K
      X=LOG(SCR3(KK))
      EM2(KK)=1.-(C02+X*(C12+X*(C22+X*C32)))*BNDI
   ENDDO
!
! Calculate the CO2-H2O overlap emissivity using the transmittance of
! water vapor given by the exponential form
!
   DO KK=2,K
      TRANS=EXP(-B1*SCR1(KK)/SQRT(1.+B2*SCR1(KK)))
      EM2(KK)=1.-EM2(KK)*TRANS
   ENDDO
!
! if at a level greater than that of the lowest cloud compute upward
! emissivity for clouds and mixed emissivities as per Goody.
!
   IF(K.GE.LCLDBS)THEN
      SCR1(1)=0.
      DO KK=2,K
         SCR1(KK)=SCR1(KK-1)+CLD(K-KK+2)
         SCR2(KK)=1.-EXP(-.13E4*SCR1(KK))
      ENDDO
      DO KK=2,K
         EM1(KK)=1.-(1.-EM1(KK))*(1.-SCR2(KK))
         EM2(KK)=1.-(1.-EM2(KK))*(1.-SCR2(KK))
      ENDDO
   ENDIF
!
! compute terms in (RTE), yielding net upward fluxes
!
   SX1=0.
   SX2=0.
   DO KK=2,K
      SX1 = SX1 + BB1(K-KK+2)*(EM1(KK)-EM1(KK-1))
      SX2 = SX2 + BB2(K-KK+2)*(EM2(KK)-EM2(KK-1))
   ENDDO
!
   UF1(K)=UF1(1)*(1.-EM1(K))+SX1
   UF2(K)=UF2(1)*(1.-EM2(K))+SX2
ENDDO
!
! ----------------------------- Downward computations
!
EM1(1)=0.
EM2(1)=0.
DO K=1,NZP
!
! Sum the vapor, dimer corrected vapor and co2 path lengths for use
! in the emissivity polynomial fits
!
   SCR1(1)=0.
   SCR2(1)=0.
   SCR3(1)=0.
   DO KK=2,K+1
      SCR1(KK)=SCR1(KK-1)+VPR(KK-K+NZP)
      SCR2(KK)=SCR2(KK-1)+DMR(KK-K+NZP)
      SCR3(KK)=SCR3(KK-1)+CO2(KK-K+NZP)
   ENDDO
!
! Find level of path length of 1E-3 gm/cm^2 and compute upward
! emissivities for water vapor.  store in array indexed IDE1
!
   DO KK=2,K+1
      IF(SCR1(KK).GT.1.E-3)GO TO 123
   ENDDO
   KL=K+2
   GO TO 124
123     CONTINUE
   KL=KK
124     CONTINUE


   DO KK=2,KL-1
      X=SQRT(SCR1(KK))
      EM1(KK)=X*(AD(1)+X*(AD(2)+X*(AD(3)+X*AD(4))))
   ENDDO
   DO KK=KL,K+1
      X=LOG(SCR1(KK))
      EM1(KK)=BD(1)+X*(BD(2)+X*(BD(3)+X*(BD(4)+X*BD(5))))
   ENDDO
!
! Correct vapor emissivities for dimer path length
!
   DO KK=2,K+1
      IF(SCR2(KK).GT.1E-3)GO TO 133
   ENDDO
   GO TO 134
133     KL=KK
!

   DO KK=KL,K+1
      X=LOG(MIN(1.,SCR2(KK)))
      EM1(KK)=EM1(KK)  &
             +ED(1)+X*(ED(2)+X*(ED(3)+X*(ED(4)+X*ED(5))))
      EM1(KK)=MIN(1.,EM1(KK))
   ENDDO
134     CONTINUE
!
! Compute upward emissivities for CO2, storing in IUE2, again finding
! the level of the critical path length
!
   DO KK=2,K+1
      IF(SCR3(KK).GT.1.E-2)GO TO 143
   ENDDO
   KL=K+2
   GO TO 144
143     CONTINUE
   KL=KK
144     CONTINUE

   DO KK=2,KL-1
      X=SQRT(SCR3(KK))
      EM2(KK)=1.-(X*(C11+X*(C21+X*C31)))*BNDI
   ENDDO
   DO KK=KL,K+1
      X=LOG(SCR3(KK))
      EM2(KK)=1.-(C02+X*(C12+X*(C22+X*C32)))*BNDI
   ENDDO
!
! Calculate the CO2-H2O overlap emissivity using the transmittance of
! water vapor given by the exponential form
!
   DO KK=2,K+1
      TRANS=EXP(-B1*SCR1(KK)/SQRT(1.+B2*SCR1(KK)))
      EM2(KK)=1.-EM2(KK)*TRANS
   ENDDO
!
! if at a level less than that of the highest cloud compute downward
! emissivity for clouds and mixed emissivity as per Goody.
!
   IF(NZP+2-K.LE.LCLDTP)THEN
      SCR1(1)=0.
      DO KK=2,K+1
         SCR1(KK)=SCR1(KK-1)+CLD(KK-K+NZP)
         SCR2(KK)=1.-EXP(-.158E4*SCR1(KK))
      ENDDO
      DO KK=2,K+1
         EM1(KK)=1.-(1.-EM1(KK))*(1.-SCR2(KK))
         EM2(KK)=1.-(1.-EM2(KK))*(1.-SCR2(KK))
      ENDDO
   ENDIF
!
! compute terms in radiative transfer equation
!
   SX1=0.
   SX2=0.
   DO KK=2,K+1
      SX1 = SX1 + BB1(KK-K+NZP)*(EM1(KK)-EM1(KK-1))
      SX2 = SX2 + BB2(KK-K+NZP)*(EM2(KK)-EM2(KK-1))
   ENDDO
!
! Compute downward fluxes for water vapor (iuf1) and co2 (iuf2)
!
   DF1(NZPP-K)=FHL*(1.-EM1(K+1))+SX1
   DF2(NZPP-K)=FHL*(1.-EM2(K+1))+SX2
ENDDO
!
! ----------------------------- Net flux & heating rates
!
DO K=1,NZP
   SCR1(K)=UF1(K)-DF1(K)+UF2(K)-DF2(K)
ENDDO
DO K=2,NZP
   FTHR(K)=-(SCR1(K)-SCR1(K-1))*DZZR(K)/(CP*DN0R(K)*PIRD(K))
ENDDO
RLONG=DF1(1)+DF2(1)

return
END SUBROUTINE lwradc

!##############################################################################
Subroutine shradp (NZP,RVR,DN0R,DZR,SC,PIRD,COSZ,ALBEDO  &
                 ,SOLAR,FTHR,RSHORT)

!----------------------------------------------------------------------
!     Shortwave radiation parameterization described in Mahrer and
!     Pielke(1977).
!
!       Arguments:
!       ----------
!       Input:  NZP    - number of vertical levels
!               RVR    - water vapor at each level
!               DN0R   - air density at each level
!               DZR    - inverse of delta z = 1./(Z(K)-Z(K-1))
!                        where Z is at the vapor levels
!               PIRD   - Exner func (p/p00)**(R/Cp) at each level
!               SC     - scratch array at least 2*NZP long
!               COSZ   - cosine of the zenith angle
!               ALBEDO - albedo of the ground surface
!               SOLAR  - the solar constant
!
!       Output: FTHR   - radiation tendency on potential temperture.
!               RSHORT - downward shortwave flux on a flat surface at
!                        the ground
!----------------------------------------------------------------------

implicit none

integer :: nzp
real :: COSZ,ALBEDO,SOLAR,RSHORT
real :: RVR(NZP),DN0R(NZP),PIRD(NZP),FTHR(NZP),DZR(NZP),SC(NZP,2)
SAVE
integer, parameter :: IV1=1,IV2=2
real, parameter :: CP=1.004E7

integer :: nz,k
real :: raysct,rdcon1,vabs
real, external :: ssumvect

NZ=NZP-1

!     Rayleigh scattering (numerator in SQRT should be
!        SQRT((.000949*P+.051)/COSZ), but is ignored. (P in mb)

RAYSCT=1.021-.0824*SQRT(1./COSZ)

!     Vapor path length

DO K=1,NZ
  SC(K,IV1)=(RVR(K)*DN0R(K)+RVR(K+1)*DN0R(K+1))*.5/DZR(K)
  SC(K,IV1)=MAX(SC(K,IV1),1E-10)
ENDDO
DO K=1,NZ
  SC(K,IV2)=ssumvect(NZP-K,SC(K,IV1),1)+1E-10
ENDDO
SC(NZP,IV2)=0.

!     Shortwave heating by vapor absorbtion

RDCON1=.0231*SOLAR/CP
DO K=2,NZ
  FTHR(K)=(RDCON1*(SC(K,IV2)/COSZ)**(-.7)*COSZ*RVR(K))/PIRD(K)
ENDDO
VABS=.077*(SC(1,IV2)/COSZ)**.3

!     Shortwave on a flat surface

RSHORT=MAX(SOLAR*COSZ*(1.-ALBEDO)*(RAYSCT-VABS),0.)

return
END SUBROUTINE shradp

!##############################################################################
Subroutine lwradp (NZP,TEMPRD,RVR,DN0R,DZZR,PIRD,SC,FTHR,RLONG)

!-----------------------------------------------------------------------
!     Longwave radiation parameterization described in Mahrer and
!     Pielke (1977).  Does not include any cloud effects.
!
!       Arguments:
!       ----------
!
!       Input:  NZP    - number of vertical levels
!               TEMPRD - temperature in Kelvin at each level
!               RVR    - water vapor at each level
!               DN0R   - air density at each level
!               DZZR   - inverse of delta z = 1./(ZZ(K)-ZZ(K-1))
!                        where ZZ is staggered midway between T levels
!               PIRD   - Exner func (p/p00)**(R/Cp) at each level
!               SC     - scratch array at least 20*NZP long
!
!       Output: FTHR   - radiation tendency on potential temperture.
!               RLONG  - downward longwave flux at the ground
!
!-----------------------------------------------------------------------

implicit none

integer :: nzp
real :: RLONG
real :: RVR(NZP),DN0R(NZP),TEMPRD(NZP),FTHR(NZP),DZZR(NZP),PIRD(NZP)
real :: SC(NZP,18)
SAVE
real, parameter :: G=980.,CP=1.004E7,STEFAN=5.6696E-5,R=.287E7,P00=1E6
integer, parameter :: IV1=1,IV2=2,IV3=3,IV4=4,IV5=5,IV6=6,IV7=7,IV8=8  &
                     ,IV9=9,IV10=10,IV11=11,IV12=12,IV13=13,IV14=14  &
                     ,IV15=15,IV16=16,IV17=17,IV18=18
integer :: nz,nz1,k
real :: c1,c2                     

NZ=NZP-1
NZ1=NZ-1

!                   COMPUTE UPWARD AND DOWNWARD VAPOR PATH
DO K=2,NZ
  SC(K,IV1)=RVR(K)*DN0R(K)/DZZR(K)
  SC(K,IV1)=MAX(SC(K,IV1),1E-10)
ENDDO
SC(NZP,IV1)=SC(NZ,IV1)
!
SC(1,IV2)=0.
DO K=2,NZ
  SC(K,IV2)=SC(K-1,IV2)+SC(K,IV1)
ENDDO
SC(NZ,IV3)=SC(NZP,IV1)
DO K=NZ1,1,-1
  SC(K,IV3)=SC(K+1,IV3)+SC(K+1,IV1)
ENDDO

!                          WATER VAPOR EMISSIVITY CALCULATION
DO K=1,NZ
  SC(K,IV4)=LOG10(SC(K,IV2)+1E-30)
  SC(K,IV5)=LOG10(SC(K,IV3)+1E-30)
ENDDO

DO K=1,NZ
  IF(SC(K,IV4).LE.-4.)  &
    SC(K,IV6)=.1129*LOG10(1.+12.63*SC(K,IV2))
  IF(SC(K,IV4).LE.-3.0.AND.SC(K,IV4).GT.-4.0)  &
    SC(K,IV6)=.104*SC(K,IV4)+.440
  IF(SC(K,IV4).LE.-1.5.AND.SC(K,IV4).GT.-3.0)  &
    SC(K,IV6)=.121*SC(K,IV4)+.491
  IF(SC(K,IV4).LE.-1.0.AND.SC(K,IV4).GT.-1.5)  &
    SC(K,IV6)=.146*SC(K,IV4)+.527
  IF(SC(K,IV4).LE. 0.0.AND.SC(K,IV4).GT.-1.0)  &
    SC(K,IV6)=.161*SC(K,IV4)+.542
  IF(SC(K,IV4).GT. 0.)  &
    SC(K,IV6)=.136*SC(K,IV4)+.542

  IF(SC(K,IV5).LE.-4.)  &
    SC(K,IV7)=.1129*LOG10(1.+12.63*SC(K,IV3))
  IF(SC(K,IV5).LE.-3.0.AND.SC(K,IV5).GT.-4.0)  &
    SC(K,IV7)=.104*SC(K,IV5)+.440
  IF(SC(K,IV5).LE.-1.5.AND.SC(K,IV5).GT.-3.0)  &
    SC(K,IV7)=.121*SC(K,IV5)+.491
  IF(SC(K,IV5).LE.-1.0.AND.SC(K,IV5).GT.-1.5)  &
    SC(K,IV7)=.146*SC(K,IV5)+.527
  IF(SC(K,IV5).LE. 0.0.AND.SC(K,IV5).GT.-1.0)  &
    SC(K,IV7)=.161*SC(K,IV5)+.542
  IF(SC(K,IV5).GT. 0.)  &
    SC(K,IV7)=.136*SC(K,IV5)+.542
ENDDO

!                           CO2 path lengths and emissivities
C1=.0004148239
C2=C1*G
DO K=2,NZ
  SC(K,IV11)=C2*DN0R(K)/DZZR(K)
ENDDO
SC(NZP,IV11)=C1*PIRD(NZP)**(CP/R)*P00

SC(1,IV12)=0.
DO K=2,NZ
  SC(K,IV12)=SC(K-1,IV12)+SC(K,IV11)
ENDDO
SC(NZ,IV13)=SC(NZP,IV11)
DO K=NZ1,1,-1
  SC(K,IV13)=SC(K+1,IV13)+SC(K+1,IV11)
ENDDO

DO K=1,NZ
  SC(K,IV8)=.185*(1.-EXP(-.3919*SC(K,IV12)**.4))
  SC(K,IV9)=.185*(1.-EXP(-.3919*SC(K,IV13)**.4))
ENDDO

!                        Add CO2 and H2O emissivities, find SIG(T**4)
DO K=1,NZP
  SC(K,IV14)=SC(K,IV8)+SC(K,IV6)
  SC(K,IV15)=SC(K,IV9)+SC(K,IV7)
  SC(K,IV16)=STEFAN*TEMPRD(K)**4
ENDDO
SC(NZP,IV16)=SC(NZP,IV16)*SC(NZ,IV15)

!                       Calculate upward and downward divergences
DO K=2,NZ
  SC(K,IV17)=(SC(K,IV16)-SC(1,IV16))*(SC(K,IV14)-SC(K-1,IV14))
  SC(K,IV18)=(SC(NZP,IV16)-SC(K,IV16))*(SC(K,IV15)-SC(K-1,IV15))
ENDDO
SC(1,IV18)=SC(NZP,IV16)*(1.-SC(1,IV15))+SC(2,IV16)*SC(2,IV15)

DO K=2,NZ
  FTHR(K)=-(SC(K,IV17)+SC(K,IV18))*DZZR(K)/(CP*DN0R(K)*PIRD(K))
ENDDO
RLONG=SC(1,IV18)

!------------------------------------------------------------------
!      PRINT 6667,(K,TEMPRD(K),SC(K,IV6),SC(K,IV7),SC(K,IV8),SC(K,IV9)
!     +  ,SC(K,IV13),FTHR(K)*24.*3600.,K=NZP,1,-1)
!6667 FORMAT(' LONGWAVE-T,EMISSUR,EMISSDR,EMISSUC,EMISSDC,PATHD,D/DAY',  &
!      /,(I3,7E10.3))
!      PRINT 6668,(K,TEMPRD(K),SC(K,IV2),SC(K,IV3),SC(K,IV12),SC(K,IV13)
!     +  ,K=NZP,1,-1)
!6668 FORMAT(' LONGWAVE-T,UP VAP,DN VAP,UP CO2, DN CO2',  &
!      /,(I3,5E10.3))
!      PRINT*,'  LONGWAVE DOWN ',-SC(1,IV18)*1E-3

return
END SUBROUTINE lwradp

!##############################################################################
Subroutine shradc (NZPP,RVR,RTR,DN0R,DZZR,PRD,PIRD,SC,ALBEDO,  &
                   SOLAR,COSZ,FTHR,RSHORT)

implicit none
                  
integer :: nzpp
real :: ALBEDO,SOLAR,COSZ,RSHORT
real :: RVR(NZPP),RTR(NZPP),DN0R(NZPP),PRD(NZPP)  &
               ,FTHR(NZPP),DZZR(NZPP),SC(NZPP,40),PIRD(NZPP)

! PIRD   - Exner func (p/p00)**(R/Cp) at each level

!  Tak changed model from 3 to 8 bands
!     DIMENSION SFCT1(3),SFCT2(3)
integer, PARAMETER :: NZMAX=200
real :: SFCT1(8),SFCT2(8),O3(NZMAX)

REAL ALBRAY   ! effective albedo including rayleigh scatter

SAVE
integer, parameter :: IV1=1,IV2=2,IV3=3,IV4=4,IV5=5,IV6=6,IV7=7,IV8=8  &
                     ,IV9=9,IV10=10,IV11=11,IV12=12,IV13=13,IV14=14  &
                     ,IV15=15,IV16=16,IV17=17,IV18=18,IV19=19,IV20=20  &
    ,IV21=21,IV22=22,IV23=23,IV24=24,IV25=25,IV26=26,IV27=27,IV28=28 &
    ,IAOZ=29,IV30=30,IV31=31,IV32=32,IREA=33,ITR1=34,ITR2=35,IV36=36 &
    ,IAB2=37,IREFS=38,IREB=39,IV40=40
integer, parameter :: IIV1=1, IIV2=2, IIV3=3

!  Tak changed to 8 bands
!     DATA SFCT1/.19649,.00132,7.8179/, SFCT2 /.12096,.80556,.07348/

DATA SFCT1/0.00004,0.002,0.035,0.377,1.95,9.40,44.6,190.0/
DATA SFCT2/0.647,0.0698,0.1443,0.0584,0.0335,0.0225,0.0158,0.0087/

real, parameter :: CP=1.004E7

integer :: nzp,nz,k,k1,nbnd,icldfl
real :: radc1,pfct,rabar,rabarbar,trsmt
real, external :: ssumvect
real, external :: cvmgp
real, external :: cvmgm

NZP=NZPP-1
NZ=NZP-1

!     O3: Ozone mixing ratio

SC(1,IV1)=0.
DO K=2,NZP
  SC(K,IV1)=SC(K-1,IV1)+1./DZZR(K)
ENDDO

SC(1,IV11)=0.
DO K=2,NZPP
  SC(K,IV11)=.4+.4*EXP(-4.)/(1.+EXP((SC(NZP+2-K,IV1)  &
  *1.E-5-20.)/5.))
  O3(NZP+3-K)=SC(K,IV11)-SC(K-1,IV11)
ENDDO

!     Compute ozone absorptance
!     REF.....LACIS,HANSEN,1974,J.A.S. P118

RADC1=35./SQRT(1224.*COSZ*COSZ+1.)
SC(1,IV2)=0.
SC(1,IV3)=0.
SC(1,IV4)=0.
SC(1,IV5)=0.
DO K=2,NZPP
  SC(K,IV2)=RADC1*O3(NZP-K+3)
ENDDO
DO K=2,NZPP
  SC(K,IV3)=ssumvect(K,SC(1,IV2),1)
  SC(K,IV4)=.02118*SC(K,IV3)/(1.+.042*SC(K,IV3)  &
    +.000323*SC(K,IV3)*SC(K,IV3))
  SC(K,IV5)=1.082*SC(K,IV3)/((1.+138.6*SC(K,IV3))**.805)  &
           +.0658*SC(K,IV3)/(1.+(103.6*SC(K,IV3))**3)
ENDDO
DO K=2,NZPP
  SC(NZP-K+3,IAOZ)=SC(K,IV4)-SC(K-1,IV4)+SC(K,IV5)-SC(K-1,IV5)
ENDDO
!
!     Precomputation of reflectance,transmittance,absorptance in cloudy
!     REF.....STEPHENS,1978,J.A.S.P2123
!
!     --- Cloud fractional coverage
!  Old way all or nothing cloud
DO K=2,NZP
  SC(K,IV36)=cvmgp(1.,0.,RTR(K)-RVR(K)-1.E-5)
ENDDO
SC(NZPP,IV36)=SC(NZP,IV36)
ICLDFL=MIN(1.,ssumvect(NZ,SC(2,IV36),1))

IF(ICLDFL.EQ.0)GO TO 600
!
!     .75UM is a line of demarcation
!
!==================================================
!  If passing liquid water mixing ratio from cloud1d change this
DO K=2,NZP
  SC(K,IV1)=1.E4*(RTR(K)-RVR(K))*DN0R(K)/DZZR(K)
ENDDO
!==================================================
SC(NZPP,IV1)=0.

!     IV2: TN1
!     IV3: TN2
DO K=2,NZPP

!  if W < 10 g/m^2 then use top two equations so tau linearly goes to zero
!  old way
!       SC(K,IV2)=.2633*SC(K,IV1)
!       SC(K,IV3)=.3492*SC(K,IV1)
!  new way from tripoli
  SC(K,IV2)=.1833*SC(K,IV1)
  SC(K,IV3)=.2234*SC(K,IV1)
  SC(K,IV4)=MAX(10.,SC(K,IV1))
  SC(K,IV5)=10.**(.2633+1.7095*LOG(LOG10(SC(K,IV4))))
  SC(K,IV6)=10.**(.3492+1.6518*LOG(LOG10(SC(K,IV4))))
ENDDO
DO K=2,NZPP
  SC(K,IV2)=cvmgp(SC(K,IV5),SC(K,IV2),SC(K,IV1)-10.)
  SC(K,IV3)=cvmgp(SC(K,IV6),SC(K,IV3),SC(K,IV1)-10.)
ENDDO
!
DO K=2,NZPP
  CALL stable (1,1,COSZ,SC(K,IV2),SC(K,IV4),SC(K,IV6))
  CALL stable (2,1,COSZ,SC(K,IV3),SC(K,IV5),SC(K,IV6))
  CALL stable (3,1,COSZ,SC(K,IV3),SC(K,IV6),SC(K,IV6))
  SC(K,IV6)=MIN(SC(K,IV6),.99999)
ENDDO
!
!     LANDA.LT..75UM
!
DO K=2,NZPP
  SC(K,IV30)=SC(K,IV4)*SC(K,IV2)
  SC(K,IREA)=SC(K,IV30)/(COSZ+SC(K,IV30))
  SC(K,IREA)=MIN(1.,MAX(0.,SC(K,IREA)))
  SC(K,ITR1)=1.-SC(K,IREA)
  SC(K,IV11)=(1.-SC(K,IV6)+2.*SC(K,IV5)*SC(K,IV6))  &
    /(1.-SC(K,IV6))
  SC(K,IV12)=SQRT(SC(K,IV11))
  SC(K,IV13)=(1.-SC(K,IV6))*(1.-SC(K,IV6)+2.*SC(K,IV5)  &
    *SC(K,IV6))
  SC(K,IV13)=SQRT(SC(K,IV13))*SC(K,IV3)/COSZ
  SC(K,IV30)=EXP(SC(K,IV13))
  SC(K,IV31)=1./SC(K,IV30)
  SC(K,IV14)=(SC(K,IV12)+1.)*(SC(K,IV12)+1.)*SC(K,IV30)  &
           -(SC(K,IV12)-1.)*(SC(K,IV12)-1.)*SC(K,IV31)
!
!       LANDA.GT..75UM
!
  SC(K,IREB)=(SC(K,IV11)-1.)*(SC(K,IV30)-SC(K,IV31))  &
    /SC(K,IV14)
  SC(K,ITR2)=4.*SC(K,IV12)/SC(K,IV14)
  SC(K,IAB2)=1.-SC(K,IREB)-SC(K,ITR2)
ENDDO
!
600 CONTINUE
!              Limit quantities and multiply by appropriate factors
!                    for the flux computations
DO K=2,NZPP
  SC(K,IV1)  = SC(K,IREB)-1.
  SC(K,IREB) = cvmgp(1.,SC(K,IREB),SC(K,IV1))
  SC(K,ITR2) = cvmgp(0.,SC(K,ITR2),SC(K,IV1))
  SC(K,IAB2) = cvmgp(0.,SC(K,IAB2),SC(K,IV1))
  SC(K,IV30) = .483*SC(K,IV36)
  SC(K,IV31) = .517*SC(K,IV36)
  SC(K,IREB) = cvmgm(0.,SC(K,IREB),SC(K,IREB))*SC(K,IV30)
  SC(K,ITR2) = cvmgm(1.,SC(K,ITR2),SC(K,IREB))*SC(K,IV30)
  SC(K,IAB2) = cvmgm(0.,SC(K,IAB2),SC(K,IREB))*SC(K,IV30)
  SC(K,ITR1) = SC(K,ITR1)*SC(K,IV31)
  SC(K,IREA) = SC(K,IREA)*SC(K,IV31)
ENDDO
!
!     Compute short wave fluxes
!
!-----------------------------------------------------
!  old way for rayleigh scatter
!     .07 reflectance due to molecular scattering
!
!          Rayleigh scattering
!             VCT1=.219/(1.+.816*COSZ)*.517/3039.E3

!     VCT1=3.7257E-8/(1.+.816*COSZ)
!     DO K=1,NZP
!       SC(K,IREFS)=PRD(K)*VCT1*(1.-SC(K,IV36))
!     ENDDO
!     SC(NZPP,IREFS)=SC(NZP,IREFS)
!
!-----------------------------------------------------
!  tripoli new way include rayleigh scatter in an effective albedo

!     Lacis and Hansen Least square fit
!     for Rayleigh reflectance gives:

!  k1 is the pressure at the first grid level above the surface
!  his model was eta type coordinate where 1 could be below ground
k1 = 2             ! this is the first level above the surface
pfct=PRD(k1)*1.e-6
Rabar=pfct*0.219/(1+0.816*cosz)
Rabarbar=pfct*0.144
albray=rabar+(1.-rabar)*(1.-rabarbar)*  &
       albedo/(1.-rabarbar*albedo)

!     The above represents the effective albedo of the lower atmosphere
!     due to Rayleigh scattering

!-----------------------------------------------------
!     Compute water vapor path
!
!  Vapor water path for clear atmosphere?
DO K=2,NZP
  SC(K,IV40)=RVR(K)*(PRD(K)/1.01325E6)**.86*RADC1  &
    *DN0R(K)/DZZR(K)
ENDDO
SC(NZPP,IV40)=SC(NZP,IV40)
!
!     REF.....STEPHENS,1977
!
DO K=1,NZPP
  SC(K,IV3)=0.
  SC(K,IV4)=0.
!  old way with rayleigh reflection added in
!       SC(K,IV12)=SC(K,IREA)+SC(K,IREB)+SC(K,IREFS)
!  new way with out rayleigh reflection (albedo is changed)
  SC(K,IV12)=SC(K,IREA)+SC(K,IREB)
  SC(K,IV30)=1.-SC(K,IV12)-SC(K,IAB2)-SC(K,IAOZ)
  SC(K,ITR1)=SC(K,ITR1)+SC(K,ITR2)
ENDDO

!  Tak switched to 8 bands
!                 Loop through 3 "pseudo-bands"
!     DO NBND=1,3

!                 Loop through 8 "pseudo-bands"
DO NBND=1,8

!
!     Compute AB,REF,TRP,TRN
!
  TRSMT=1.
  DO K=2,NZPP

!  new tripoli way with Trans of clear * cloud
!  because cloud tr does not include water vapor
    SC(K,IV1)=(SC(K,ITR1)+(1.-SC(K,IV36)))  &
       *EXP(-SFCT1(NBND)*SC(K,IV40))

!  old chen way with just liquid water effects
!         SC(K,IV1)=SC(K,ITR1)+(1.-SC(K,IV36))
!    +       *EXP(-SFCT1(NBND)*SC(K,IV40))

    SC(K,IV16)=(1.-SC(K,IV36))*(1.-SC(K,IV1))
  ENDDO

  DO K=NZPP,2,-1
    SC(K,IV13)=SC(K,IV30)-SC(K,IV16)*TRSMT
    TRSMT=TRSMT*MAX(0.,SC(K,IV13))
  ENDDO
  DO K=2,NZPP
    SC(K,IV14)=SC(K,IV30)-SC(K,IV16)*TRSMT
    TRSMT=TRSMT*MAX(0.,SC(K,IV14))
  ENDDO
!
!     REF.....STEPHENS,1979,J.A.S. P1542
  DO K=1,NZP
    SC(K,IV21)=SC(NZPP-K+1,IV12)
    SC(K,IV22)=SC(NZPP-K+1,IV13)
    SC(K,IV23)=SC(NZPP-K+1,IV14)
  ENDDO
  SC(1,IV26)=0.
  SC(1,IV28)=SOLAR*COSZ
  DO K=1,NZP
    SC(K,IV18)=1./(1.-SC(K,IV26)*SC(K,IV21))
    SC(K+1,IV26)=SC(K,IV21)+SC(K,IV22)*SC(K,IV26)  &
      *SC(K,IV23)*SC(K,IV18)
    SC(K+1,IV27)=SC(K,IV21)*SC(K,IV28)*SC(K,IV18)
    SC(K+1,IV28)=SC(K,IV22)*SC(K,IV28)*SC(K,IV18)
  ENDDO

!  old way with old rayleigh scatter
!       SC(NZPP,IV24)=SC(NZPP,IV28)/(1.-SC(NZPP,IV26)*ALBEDO)
!       SC(NZPP,IV25)=ALBEDO*SC(NZPP,IV24)

!  new way with effective albedo from tripoli
  SC(NZPP,IV24)=SC(NZPP,IV28)/(1.-SC(NZPP,IV26)*ALBRAY)
  SC(NZPP,IV25)=ALBRAY*SC(NZPP,IV24)
  DO K=NZP,1,-1
    SC(K,IV25)=SC(K,IV23)*SC(K+1,IV25)  &
             /(1.-SC(K,IV26)*SC(K,IV21))+SC(K+1,IV27)
  ENDDO
  DO K=2,NZPP
    SC(K,IV24)=SC(K,IV26)*SC(K,IV25)+SC(K,IV28)
  ENDDO
!       IV3=FLXU   IV4=FLXD
  DO K=2,NZPP
    SC(NZPP+1-K,IV3)=SC(NZPP+1-K,IV3)+SFCT2(NBND)*SC(K,IV25)
    SC(NZPP+1-K,IV4)=SC(NZPP+1-K,IV4)+SFCT2(NBND)*SC(K,IV24)
  ENDDO
ENDDO
!
!     Compute shortwave radiative tendency
!
DO K=1,NZP
  SC(K,IV11)=SC(K,IV3)-SC(K,IV4)
ENDDO
DO K=2,NZP
  FTHR(K)=-(SC(K,IV11)-SC(K-1,IV11))*DZZR(K)/(CP*DN0R(K)*PIRD(K))
!     write(9,*) ' short k fthr fc ',k,fthr(k),fcr(k)
ENDDO
RSHORT=-SC(1,IV11)
!
!      DO K=1,NZP
!        PRINT*,' D/DAY ',K,FTHR(K)*86400.,SC(K,IV3),SC(K,IV4)
!      ENDDO
return
END SUBROUTINE shradc
