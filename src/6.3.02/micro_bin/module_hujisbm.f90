!##############################################################################
Module module_hujisbm

use micro_prm
use micphys, only: imbudget,ipris,igraup,ihail

! YWLL_1000MB(nkr,nkr) - input array of kernels for pressure 1000mb
! YWLL_750MB(nkr,nkr) - input array of kernels for pressure 750mb
! YWLL_500MB(nkr,nkr) - input array of kernels for pressure 500mb
REAL, SAVE :: &
! CRYSTALS 
   YWLI(NKR,NKR,ICEMAX) &
! MIXTURES
   ,YWIL(NKR,NKR,ICEMAX),YWII(NKR,NKR,ICEMAX,ICEMAX) &
   ,YWIS(NKR,NKR,ICEMAX),YWIG(NKR,NKR,ICEMAX) &
   ,YWIH(NKR,NKR,ICEMAX),YWSI(NKR,NKR,ICEMAX) &
   ,YWGI(NKR,NKR,ICEMAX),YWHI(NKR,NKR,ICEMAX)
!
REAL, DIMENSION(NKR,NKR),SAVE :: &
    YWLL_1000MB,YWLL_750MB,YWLL_500MB,YWLL,YWLS,YWLG,YWLH &
! SNOW :
    ,YWSL,YWSS,YWSG,YWSH &
! GRAUPELS :
    ,YWGL,YWGS,YWGG,YWGH &
! HAIL :
    ,YWHL,YWHS,YWHG,YWHH
REAL,SAVE :: &
       XI(NKR,ICEMAX) &
      ,RADXX(NKR,NHYDR-1),MASSXX(NKR,NHYDR-1),DENXX(NKR,NHYDR-1) &
      ,RADXXO(NKR,NHYDRO),MASSXXO(NKR,NHYDRO),DENXXO(NKR,NHYDRO) &
      ,RIEC(NKR,ICEMAX),COEFIN(NKR),SLIC(NKR,6),TLIC(NKR,2) &
      ,RO2BL(NKR,ICEMAX)
REAL, SAVE :: VR1(NKR),VR2(NKR,ICEMAX),VR3(NKR) &
      ,VR4(NKR),VR5(NKR),VRX(NKR)
REAL,DIMENSION(NKR),SAVE ::  &
       XL,RLEC,XX,XCCN,XS,RSEC &
      ,XG,RGEC,XH,RHEC,RO1BL,RO3BL,RO4BL,RO5BL &
      ,ROCCN,RCCN,DROPRADII
REAL, SAVE ::  ICEN(NKR)
REAL, SAVE ::  FCCNR0(NKR), FCCNR_mp(NKR)

REAL :: C2,C3,C4
double precision,save ::  cwll(nkr,nkr)
double precision,save::  &
      xl_mg(0:nkr),xs_mg(0:nkr),xg_mg(0:nkr),xh_mg(0:nkr) &
     ,xi1_mg(0:nkr),xi2_mg(0:nkr),xi3_mg(0:nkr) &
     ,chucm(nkr,nkr),ima(nkr,nkr) &
     ,cwll_1000mb(nkr,nkr),cwll_750mb(nkr,nkr),cwll_500mb(nkr,nkr) &
     ,cwli_1(nkr,nkr),cwli_2(nkr,nkr),cwli_3(nkr,nkr) &
     ,cwls(nkr,nkr),cwlg(nkr,nkr),cwlh(nkr,nkr) &

     ,cwil_1(nkr,nkr),cwil_2(nkr,nkr),cwil_3(nkr,nkr) &

     ,cwii_1_1(nkr,nkr),cwii_1_2(nkr,nkr),cwii_1_3(nkr,nkr) &
     ,cwii_2_1(nkr,nkr),cwii_2_2(nkr,nkr),cwii_2_3(nkr,nkr) &
     ,cwii_3_1(nkr,nkr),cwii_3_2(nkr,nkr),cwii_3_3(nkr,nkr) &

     ,cwis_1(nkr,nkr),cwis_2(nkr,nkr),cwis_3(nkr,nkr) &
     ,cwig_1(nkr,nkr),cwig_2(nkr,nkr),cwig_3(nkr,nkr) &
     ,cwih_1(nkr,nkr),cwih_2(nkr,nkr),cwih_3(nkr,nkr) &

     ,cwsl(nkr,nkr) &
     ,cwsi_1(nkr,nkr),cwsi_2(nkr,nkr),cwsi_3(nkr,nkr)&
     ,cwss(nkr,nkr),cwsg(nkr,nkr),cwsh(nkr,nkr) &
     ,cwgl(nkr,nkr)&
     ,cwgi_1(nkr,nkr),cwgi_2(nkr,nkr),cwgi_3(nkr,nkr)&
     ,cwgs(nkr,nkr),cwgg(nkr,nkr),cwgh(nkr,nkr) &

     ,cwhl(nkr,nkr) &
     ,cwhi_1(nkr,nkr),cwhi_2(nkr,nkr),cwhi_3(nkr,nkr) &
     ,cwhs(nkr,nkr),cwhg(nkr,nkr),cwhh(nkr,nkr) &
     ,dlnr &
     ,CTURBLL(KRMAX_LL,KRMAX_LL)&
     ,CTURB_LL(K0_LL,K0_LL)&
     ,CTURBGL(KRMAXG_GL,KRMAXL_GL)&
     ,CTURB_GL(K0G_GL,K0L_GL)

DOUBLE PRECISION, save :: &
     BRKWEIGHT(JBREAK),PKIJ(JBREAK,JBREAK,JBREAK), &
     QKJ(JBREAK,JBREAK),ECOALMASSM(NKR,NKR)

contains

!##############################################################################
real Function sum_mass (speciesbin,xbin,dens,krs,kre)

implicit none

integer :: krs,kre
real, dimension(nkr) :: speciesbin,xbin
real :: dens

sum_mass = col3/dens*sum(speciesbin(krs:kre)*xbin(krs:kre)*xbin(krs:kre))

return
END FUNCTION sum_mass

!##############################################################################
Subroutine BREAKINIT ()

!     GT    : MASS DISTRIBUTION function
!     XT_MG : MASS OF BIN IN MG
!     JMAX  : NUMBER OF BINS

IMPLICIT NONE

character(len=strl1) :: hname
INTEGER,PARAMETER :: hujisbm_unit1=22
LOGICAL, PARAMETER :: PRINT_diag=.FALSE.
INTEGER AP,IE,JE,KE
PARAMETER (AP = 1)
INTEGER I,J,K,JDIFF
REAL PI,D0
DOUBLE PRECISION M(0:JBREAK),ALM
INTEGER IP,KP,JP,KQ,JQ
CHARACTER*20 FILENAME_P,FILENAME_Q

FILENAME_P = 'coeff_p.asc'
FILENAME_Q = 'coeff_q.asc'

IE = JBREAK
JE = JBREAK
KE = JBREAK
PI    = 3.1415927
D0    = 0.0101593
M(1)  = PI/6.0 * D0**3

!.....IN CGS
!.....SHIFT BETWEEN COAGULATION AND BREAKUP GRID

JDIFF = JMAX - JBREAK

!.....INITIALIZATION
!........CALCULATING THE BREAKUP GRID
ALM  = 2.d0
M(0)  = M(1)/ALM
DO K=1,KE-1
   M(K+1) = M(K)*ALM
ENDDO
DO K=1,KE
   BRKWEIGHT(K) = 2./(M(K)**2 - M(K-1)**2)
ENDDO

!........OUTPUT
!WRITE (*,*) 'COLL_BREAKUP_INI: COAGULATION AND BREAKUP GRID'
!WRITE (*,'(2A5,5A15)') 'ICOAG','IBREAK','XCOAG','DCOAG', &
!    'XBREAK','DBREAK','MWEIGHT'

!........READ DER BREAKUP COEFFICIENTS FROM INPUT FILE

!WRITE (*,*) 'COLL_BREAKUP: READ THE BREAKUP COEFFS'
!WRITE (*,*) '              FILE PKIJ: ', FILENAME_P
hname=trim(hucmfile)//'/coeff_p.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
DO K=1,KE
   DO I=1,IE
      DO J=1,I
         READ(hujisbm_unit1,'(3I6,1E16.8)') KP,IP,JP,PKIJ(KP,IP,JP)
      ENDDO
   ENDDO
ENDDO
CLOSE(hujisbm_unit1)

hname=trim(hucmfile)//'/coeff_q.asc'
OPEN(UNIT=hujisbm_unit1,FILE=hname,FORM="FORMATTED",STATUS="OLD")
DO K=1,KE
   DO J=1,JE
      READ(hujisbm_unit1,'(2I6,1E16.8)') KQ,JQ,QKJ(KQ,JQ)
   ENDDO
ENDDO
CLOSE(hujisbm_unit1)

DO I=1,JMAX
   DO J=1,JMAX
      ECOALMASSM(I,J)=1.0D0
   ENDDO
ENDDO

DO I=1,JMAX
   DO J=1,JMAX
      ECOALMASSM(I,J)=ECOALMASS(XL(I),XL(J))
   ENDDO
ENDDO

RETURN
END SUBROUTINE BREAKINIT

!##############################################################################
real Function ECOALMASS (ETA,KSI)

IMPLICIT NONE
REAL PI
PARAMETER (PI = 3.1415927)

REAL ETA,KSI
REAL KPI,RHO
REAL DETA,DKSI

PARAMETER (RHO  = 1.0)

!     REAL ECOALDIAM
!     EXTERNAL ECOALDIAM

KPI = 6./PI

DETA = (KPI*ETA/RHO)**(1./3.)
DKSI = (KPI*KSI/RHO)**(1./3.)

ECOALMASS = ECOALDIAM(DETA,DKSI)

RETURN
END FUNCTION ECOALMASS

!#############################################################################
real Function ECOALDIAM (DETA,DKSI)

!------------------------------------------------
!     COALESCENCE EFFICIENCY AS FUNC OF DIAMETERS
!------------------------------------------------

REAL DETA,DKSI
REAL DGR,DKL,RGR,RKL,P,Q,E,X,QMIN,QMAX
REAL ZERO,ONE,EPS,PI

PARAMETER (ZERO = 0.0)
PARAMETER (ONE  = 1.0)
PARAMETER (EPS  = 1.0E-30)
PARAMETER (PI   = 3.1415927)

DGR = MAX(DETA,DKSI)
DKL = MIN(DETA,DKSI)

RGR = 0.5*DGR
RKL = 0.5*DKL

P = (RKL / RGR)
Q = (RKL * RGR)**0.5
Q = 0.5 * (RKL + RGR)

qmin = 250e-4
qmax = 400e-4  
if (q.lt.qmin) then
   e = max(ECOALOCHS(Dgr,Dkl),ecoalBeard(Dgr,Dkl)) 
elseif (q.ge.qmin.and.q.lt.qmax) then
   x = (q - qmin) / (qmax - qmin)
   e = sin(pi/2.0*x)**2 * ECOALLOWLIST(Dgr,Dkl) &
          + sin(pi/2.0*(1 - x))**2 * ECOALOCHS(Dgr,Dkl)
elseif (q.ge.qmax) then
   e = ECOALLOWLIST(Dgr,Dkl)
else
   e  = 1.0
endif

ECOALDIAM  = MAX(MIN(ONE,E),EPS)

RETURN
END FUNCTION ECOALDIAM

!############################################################################
real Function ECOALLOWLIST (DGR,DKL)

!--------------------------------------------------
!     COALESCENCE EFFICIENCY (LOW&LIST)
!--------------------------------------------------
IMPLICIT NONE
REAL PI,SIGMA,KA,KB,EPSI
REAL DGR,DKL,RGR,RKL
REAL ST,SC,ET,DSTSC,CKE,W1,W2,DC,ECL
REAL QQ0,QQ1,QQ2

PARAMETER (EPSI=1.E-20)

PI = 3.1415927
SIGMA = 72.8
KA = 0.778
KB = 2.61E-4

RGR = 0.5*DGR
RKL = 0.5*DKL

CALL COLLENERGY (DGR,DKL,CKE,ST,SC,W1,W2,DC)

DSTSC = ST-SC
ET = CKE+DSTSC
IF (ET .LT. 50.0) THEN
   QQ0=1.0+(DKL/DGR)
   QQ1=KA/QQ0**2
   QQ2=KB*SIGMA*(ET**2)/(SC+EPSI)
   ECL=QQ1*EXP(-QQ2)
ELSE
   ECL=0.0
ENDIF

ECOALLOWLIST = ECL

RETURN
END FUNCTION ECOALLOWLIST

!###########################################################################
real Function ECOALOCHS (D_L,D_S)

!--------------------------------------------------
!     COALESCENCE EFFICIENCY (BEARD AND OCHS)
!--------------------------------------------------
IMPLICIT NONE
REAL D_L,D_S
REAL PI,SIGMA,N_W,R_S,R_L,DV,P,G,X,E
REAL EPSF,FPMIN
PARAMETER (EPSF  = 1.E-30)
PARAMETER (FPMIN = 1.E-30)

PI = 3.1415927
SIGMA = 72.8

R_S = 0.5 * D_S
R_L = 0.5 * D_L
P   = R_S / R_L

DV  = ABS(VTBEARD(D_L) - VTBEARD(D_S))
IF (DV.LT.FPMIN) DV = FPMIN
N_W = R_S * DV**2 / SIGMA
G   = 2**(3./2.)/(6.*PI) * P**4 * (1.+ P) / ((1.+P**2)*(1.+P**3))
X   = N_W**(0.5) * G
E   = 0.767 - 10.14 * X

ECOALOCHS = E

RETURN
END FUNCTION ECOALOCHS

!##############################################################################
Subroutine COLLENERGY (DGR,DKL,CKE,ST,SC,W1,W2,DC)

!-----------------------------------------
!     CALCULATING THE COLLISION ENERGY
!-----------------------------------------

IMPLICIT NONE

REAL DGR,DKL,DC
REAL K10,PI,SIGMA,RHO
REAL CKE,W1,W2,ST,SC
REAL DGKA3,DGKB3,DGKA2
REAL V1,V2,DV
REAL EPSF,FPMIN
PARAMETER (EPSF  = 1.E-30)
PARAMETER (FPMIN = 1.E-30)

PI    = 3.1415927
RHO   = 1.0
SIGMA = 72.8

K10=RHO*PI/12.0D0

DGR = MAX(DGR,EPSF)
DKL = MAX(DKL,EPSF)

DGKA2=(DGR**2)+(DKL**2)

DGKA3=(DGR**3)+(DKL**3)

IF (DGR.NE.DKL) THEN
   V1 = VTBEARD(DGR)
   V2 = VTBEARD(DKL)
   DV = (V1-V2)
   IF (DV.LT.FPMIN) DV = FPMIN
   DV = DV**2
   IF (DV.LT.FPMIN) DV = FPMIN
   DGKB3=(DGR**3)*(DKL**3)
   CKE = K10 * DV * DGKB3/DGKA3
ELSE
   CKE = 0.0D0
ENDIF
ST = PI*SIGMA*DGKA2
SC = PI*SIGMA*DGKA3**(2./3.)

W1=CKE/(SC+EPSF)
W2=CKE/(ST+EPSF)

DC=DGKA3**(1./3.)

RETURN
END SUBROUTINE COLLENERGY

!##############################################################################
real Function VTBEARD (DIAM)

!--------------------------------------------------
!     CALCULATING TERMINAL VELOCITY (BEARD-FORMULA)
!--------------------------------------------------

IMPLICIT NONE

REAL DIAM,AA
REAL ROP,RU,AMT,PP,RL,TT,ETA,DENS,CD,D,A
REAL ALA,GR,SI,BOND,PART,XX,YY,RE,VT
REAL B00,B11,B22,B33,B44,B55,B0,B1,B2,B3,B4,B5,B6
INTEGER ID

DATA B00,B11,B22,B33,B44,B55,B0,B1,B2,B3,B4,B5,B6/-5.00015, &
     5.23778,-2.04914,.475294,-.0542819,.00238449,-3.18657,.992696, &
     -.153193E-2,-.987059E-3,-.578878E-3,.855176E-4,-.327815E-5/

AA   = DIAM/2.0
ROP  = 1.0
RU   = 8.3144E+7
AMT  = 28.9644
ID   = 10000
PP   = FLOAT(ID)*100.
RL   = RU/AMT
TT   = 283.15
ETA  = (1.718+.0049*(TT-273.15))*1.E-4
DENS = PP/TT/RL
ALA  = 6.6E-6*1.01325E+6/PP*TT/293.15
GR   = 979.69
SI   = 76.1-.155*(TT-273.15)

IF (AA.GT.500.E-4) THEN
   BOND = GR*(ROP-DENS)*AA*AA/SI
   PART = (SI**3*DENS*DENS/(ETA**4*GR*(ROP-DENS)))**(1./6.)
   XX = LOG(16./3.*BOND*PART)
   YY = B00+B11*XX+B22*XX*XX+B33*XX**3+B44*XX**4+B55*XX**5
   RE = PART*EXP(YY)
   VT = ETA*RE/2./DENS/AA
ELSEIF (AA.GT.1.E-3) THEN
   CD = 32.*AA*AA*AA*(ROP-DENS)*DENS*GR/3./ETA/ETA
   XX = LOG(CD)
   RE = EXP(B0+B1*XX+B2*XX*XX+B3*XX**3+B4*XX**4+B5*XX**5+B6*XX**6)
   D  = CD/RE/24.-1.
   VT = ETA*RE/2./DENS/AA
ELSE
   A  = 1.+1.26*ALA/AA
   A  = A*2.*AA*AA*GR*(ROP-DENS)/9./ETA
   CD = 12*ETA/A/AA/DENS
   VT = A
ENDIF

VTBEARD = VT

RETURN
END FUNCTION VTBEARD

!########################################################################
real Function ecoalBeard (D_l,D_s)

!-------------------------------------------------- 
!     function f. Coalescence-Efficiency 
!     Eq. (7) of Beard and Ochs (1995)
!--------------------------------------------------
 
IMPLICIT NONE 
REAL D_l,D_s
REAL R_s,R_l
REAL rcoeff
REAL epsf
PARAMETER (epsf  = 1.e-30) 
INTEGER its
COMPLEX acoeff(4),x

R_s = 0.5 * D_s
R_l = 0.5 * D_l

rcoeff = 5.07 - log(R_s*1e4) - log(R_l*1e4/200.0)

acoeff(1) = CMPLX(rcoeff)
acoeff(2) = CMPLX(-5.94)
acoeff(3) = CMPLX(+7.27)
acoeff(4) = CMPLX(-5.29)

x = (0.50,0)

CALL laguer (acoeff,3,x,its)

ecoalBeard = REAL(x)

RETURN 
END FUNCTION ecoalBeard 

!##############################################################################
Subroutine laguer (a,m,x,its)

INTEGER m,its,MAXIT,MR,MT
REAL EPSS
COMPLEX a(m+1),x
PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
INTEGER iter,j
REAL abx,abp,abm,err,frac(MR)
COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
SAVE frac
DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
do iter=1,MAXIT
   its=iter
   b=a(m+1)
   err=abs(b)
   d=cmplx(0.,0.)
   f=cmplx(0.,0.)
   abx=abs(x)
   do j=m,1,-1
      f=x*f+d
      d=x*d+b
      b=x*b+a(j)
      err=abs(b)+abx*err
   enddo
   err=EPSS*err
   if(abs(b).le.err) then
      return
   else
      g=d/b
      g2=g*g
      h=g2-2.*f/b
      sq=sqrt((m-1)*(m*h-g2))
      gp=g+sq
      gm=g-sq
      abp=abs(gp)
      abm=abs(gm)
      if(abp.lt.abm) gp=gm
      if (max(abp,abm).gt.0.) then
         dx=m/gp
      else
         dx=exp(cmplx(log(1.+abx),float(iter)))
      endif
   endif
   x1=x-dx
   if(x.eq.x1)return
   if (mod(iter,MT).ne.0) then
      x=x1
   else
      x=x-dx*frac(iter/MT)
   endif
enddo
return
END SUBROUTINE laguer

!##############################################################################
Subroutine courant_bott ()

implicit none
integer k,kk,j,i
double precision x0
! ima(i,j) - k-category number,
! chucm(i,j)   - courant number :
! logarithmic grid distance(dlnr) :
xl_mg(0)=xl_mg(1)/2
do i=1,nkr
   do j=i,nkr
      x0=xl_mg(i)+xl_mg(j)
      do k=j,nkr
         kk=k
         if(xl_mg(k).ge.x0.and.xl_mg(k-1).lt.x0) then
           chucm(i,j)=dlog(x0/xl_mg(k-1))/(3.d0*dlnr)
           if(chucm(i,j).gt.1.-1.d-08) then
             chucm(i,j)=0.
             kk=kk+1
           endif
           ima(i,j)=min(nkr-1,kk-1)

           goto 2000
         endif
      enddo
 2000       continue
      chucm(j,i)=chucm(i,j)
      ima(j,i)=ima(i,j)
   enddo
enddo
return
END SUBROUTINE courant_bott

!##############################################################################
Subroutine kernals (DTIME)

! KHAIN30/07/99
IMPLICIT NONE
INTEGER I,J
REAL PI
!******************************************************************
data pi/3.141592654/
! dtime - timestep of integration (calculated in main program) :
! dlnr - logarithmic grid distance
! ima(i,j) - k-category number, c(i,j) - courant number 
! cw*(i,j) (in cm**3) - multiply help kernel with constant 
! timestep(dt) and logarithmic grid distance(dlnr) :
REAL DTIME
! logarithmic grid distance(dlnr) :
! dlnr=dlog(2.d0)/(3.d0*scal)
! scal is micro.prm file parameter(scal=1.d0 for x(k+1)=x(k)*2)
! calculation of cw*(i,j) (in cm**3) - multiply help kernel 
! with constant timestep(dtime) and logarithmic grid distance(dlnr) :
!     print*,'dlnr in kernal = ',dlnr,dtime
  DO I=1,NKR
     DO J=1,NKR
        CWLL_1000MB(I,J)=DTIME*DLNR*YWLL_1000MB(I,J)
        CWLL_750MB(I,J)=DTIME*DLNR*YWLL_750MB(I,J)
        CWLL_500MB(I,J)=DTIME*DLNR*YWLL_500MB(I,J)

        CWLL(I,J)=DTIME*DLNR*YWLL(I,J)
        CWLS(I,J)=DTIME*DLNR*YWLS(I,J)
        CWLG(I,J)=DTIME*DLNR*YWLG(I,J)
        CWLH(I,J)=DTIME*DLNR*YWLH(I,J)

        CWSL(I,J)=DTIME*DLNR*YWSL(I,J)
        CWSS(I,J)=DTIME*DLNR*YWSS(I,J)
        CWSG(I,J)=DTIME*DLNR*YWSG(I,J)
        CWSH(I,J)=DTIME*DLNR*YWSH(I,J)

        CWGL(I,J)=DTIME*DLNR*YWGL(I,J)
        IF(RADXXO(I,6).LT.2.0D-2) THEN
          IF(RADXXO(J,1).LT.1.0D-3) THEN
            IF(RADXXO(J,1).GE.7.0D-4) THEN
              CWGL(I,J)=DTIME*DLNR*YWGL(I,J)/1.5D0
            ELSE
              CWGL(I,J)=DTIME*DLNR*YWGL(I,J)/3.0D0
            ENDIF
          ENDIF
        ENDIF
        IF(I.LE.14.AND.J.LE.7) CWGL(I,J)=0.0D0
        CWGS(I,J)=DTIME*DLNR*YWGS(I,J)
        CWGG(I,J)=DTIME*DLNR*YWGG(I,J)
        CWGH(I,J)=DTIME*DLNR*YWGH(I,J)

        CWHL(I,J)=DTIME*DLNR*YWHL(I,J)
        CWHS(I,J)=DTIME*DLNR*YWHS(I,J)
        CWHG(I,J)=DTIME*DLNR*YWHG(I,J)
        CWHH(I,J)=DTIME*DLNR*YWHH(I,J)

        CWLI_1(I,J)=DTIME*DLNR*YWLI(I,J,1)
        CWLI_2(I,J)=DTIME*DLNR*YWLI(I,J,2)
        CWLI_3(I,J)=DTIME*DLNR*YWLI(I,J,3)
        
        CWIL_1(I,J)=DTIME*DLNR*YWIL(I,J,1)
        CWIL_2(I,J)=DTIME*DLNR*YWIL(I,J,2)
        CWIL_3(I,J)=DTIME*DLNR*YWIL(I,J,3)

        CWIS_1(I,J)=DTIME*DLNR*YWIS(I,J,1)
        CWIS_2(I,J)=DTIME*DLNR*YWIS(I,J,2)
        CWIS_3(I,J)=DTIME*DLNR*YWIS(I,J,3)

        CWSI_1(I,J)=DTIME*DLNR*YWSI(I,J,1)
        CWSI_2(I,J)=DTIME*DLNR*YWSI(I,J,2)
        CWSI_3(I,J)=DTIME*DLNR*YWSI(I,J,3)

        CWIG_1(I,J)=DTIME*DLNR*YWIG(I,J,1)
        CWIG_2(I,J)=DTIME*DLNR*YWIG(I,J,2)
        CWIG_3(I,J)=DTIME*DLNR*YWIG(I,J,3)

        CWGI_1(I,J)=DTIME*DLNR*YWGI(I,J,1)
        CWGI_2(I,J)=DTIME*DLNR*YWGI(I,J,2)
        CWGI_3(I,J)=DTIME*DLNR*YWGI(I,J,3)

        CWIH_1(I,J)=DTIME*DLNR*YWIH(I,J,1)
        CWIH_2(I,J)=DTIME*DLNR*YWIH(I,J,2)
        CWIH_3(I,J)=DTIME*DLNR*YWIH(I,J,3)

        CWHI_1(I,J)=DTIME*DLNR*YWHI(I,J,1)
        CWHI_2(I,J)=DTIME*DLNR*YWHI(I,J,2)
        CWHI_3(I,J)=DTIME*DLNR*YWHI(I,J,3)

        CWII_1_1(I,J)=DTIME*DLNR*YWII(I,J,1,1)
        CWII_1_2(I,J)=DTIME*DLNR*YWII(I,J,1,2)
        CWII_1_3(I,J)=DTIME*DLNR*YWII(I,J,1,3)

        CWII_2_1(I,J)=DTIME*DLNR*YWII(I,J,2,1)
        CWII_2_2(I,J)=DTIME*DLNR*YWII(I,J,2,2)
        CWII_2_3(I,J)=DTIME*DLNR*YWII(I,J,2,3)

        CWII_3_1(I,J)=DTIME*DLNR*YWII(I,J,3,1)
        CWII_3_2(I,J)=DTIME*DLNR*YWII(I,J,3,2)
        CWII_3_3(I,J)=DTIME*DLNR*YWII(I,J,3,3)
     ENDDO
  ENDDO
  CALL TURBCOEF ()
  DO J=1,7
     DO I=15,24-J
        CWGL(I,J)=0.0D0
     ENDDO
  ENDDO
  DO I=1,NKR
     DO J=1,NKR
        CWLG(J,I)=CWGL(I,J)
     ENDDO
  ENDDO
  DO I=KRMING_GL,KRMAXG_GL
     DO J=KRMINL_GL,KRMAXL_GL
       IF (ICETURB.EQ.1)THEN
         CWGL(I,J)=CTURBGL(I,J)*CWGL(I,J)
       ELSE
         CWGL(I,J)=CWGL(I,J)
       END IF
     ENDDO
  ENDDO
  DO I=KRMING_GL,KRMAXG_GL
     DO J=KRMINL_GL,KRMAXL_GL
        CWLG(J,I)=CWGL(I,J)
     ENDDO
  ENDDO

RETURN
END SUBROUTINE kernals

!##############################################################################
Subroutine TURBCOEF ()

IMPLICIT NONE
INTEGER I,J
DOUBLE PRECISION X_KERN,Y_KERN
DOUBLE PRECISION RL_LL(K0_LL),RL_GL(K0L_GL),RG_GL(K0G_GL)
RL_LL(1)=RADXXO(KRMIN_LL,1)*1.E4
RL_LL(2)=10.0D0
RL_LL(3)=20.0D0
RL_LL(4)=30.0D0
RL_LL(5)=40.0D0
RL_LL(6)=50.0D0
RL_LL(7)=60.0D0
RL_LL(8)=RADXXO(KRMAX_LL,1)*1.E4
DO J=1,K0_LL
   DO I=1,K0_LL
      CTURB_LL(I,J)=1.0D0
   ENDDO
ENDDO 
  CTURB_LL(1,1)=4.50D0
  CTURB_LL(1,2)=4.50D0
  CTURB_LL(1,3)=3.00D0
  CTURB_LL(1,4)=2.25D0
  CTURB_LL(1,5)=1.95D0
  CTURB_LL(1,6)=1.40D0
  CTURB_LL(1,7)=1.40D0
  CTURB_LL(1,8)=1.40D0

  CTURB_LL(2,1)=4.50D0
  CTURB_LL(2,2)=4.50D0
  CTURB_LL(2,3)=3.00D0
  CTURB_LL(2,4)=2.25D0
  CTURB_LL(2,5)=1.95D0
  CTURB_LL(2,6)=1.40D0
  CTURB_LL(2,7)=1.40D0
  CTURB_LL(2,8)=1.40D0

  CTURB_LL(3,1)=3.00D0
  CTURB_LL(3,2)=3.00D0
  CTURB_LL(3,3)=2.70D0
  CTURB_LL(3,4)=2.25D0
  CTURB_LL(3,5)=1.65D0
  CTURB_LL(3,6)=1.40D0
  CTURB_LL(3,7)=1.40D0
  CTURB_LL(3,8)=1.40D0

  CTURB_LL(4,1)=2.25D0
  CTURB_LL(4,2)=2.25D0
  CTURB_LL(4,3)=2.25D0
  CTURB_LL(4,4)=1.95D0
  CTURB_LL(4,5)=1.65D0
  CTURB_LL(4,6)=1.40D0
  CTURB_LL(4,7)=1.40D0
  CTURB_LL(4,8)=1.40D0

  CTURB_LL(5,1)=1.95D0
  CTURB_LL(5,2)=1.95D0
  CTURB_LL(5,3)=1.65D0
  CTURB_LL(5,4)=1.65D0
  CTURB_LL(5,5)=1.65D0
  CTURB_LL(5,6)=1.40D0
  CTURB_LL(5,7)=1.40D0
  CTURB_LL(5,8)=1.40D0

  CTURB_LL(6,1)=1.40D0
  CTURB_LL(6,2)=1.40D0
  CTURB_LL(6,3)=1.40D0
  CTURB_LL(6,4)=1.40D0
  CTURB_LL(6,5)=1.40D0
  CTURB_LL(6,6)=1.40D0
  CTURB_LL(6,7)=1.40D0
  CTURB_LL(6,8)=1.40D0

  CTURB_LL(7,1)=1.40D0
  CTURB_LL(7,2)=1.40D0
  CTURB_LL(7,3)=1.40D0
  CTURB_LL(7,4)=1.40D0
  CTURB_LL(7,5)=1.40D0
  CTURB_LL(7,6)=1.40D0
  CTURB_LL(7,7)=1.40D0
  CTURB_LL(7,8)=1.40D0

  CTURB_LL(8,1)=1.40D0
  CTURB_LL(8,2)=1.40D0
  CTURB_LL(8,3)=1.40D0
  CTURB_LL(8,4)=1.40D0
  CTURB_LL(8,5)=1.40D0
  CTURB_LL(8,6)=1.40D0
  CTURB_LL(8,7)=1.40D0
  CTURB_LL(8,8)=1.40D0
DO J=1,K0_LL
   DO I=1,K0_LL
      CTURB_LL(I,J)=(CTURB_LL(I,J)-1.0D0)/1.5D0+1.0D0
   ENDDO
ENDDO
DO I=KRMIN_LL,KRMAX_LL
   DO J=KRMIN_LL,KRMAX_LL
      CTURBLL(I,J)=1.0D0
   ENDDO
ENDDO
DO I=KRMIN_LL,KRMAX_LL
   X_KERN=RADXXO(I,1)*1.0D4
   IF(X_KERN.LT.RL_LL(1)) X_KERN=RL_LL(1)
   IF(X_KERN.GT.RL_LL(K0_LL)) X_KERN=RL_LL(K0_LL) 
   DO J=KRMIN_LL,KRMAX_LL
      Y_KERN=RADXXO(J,1)*1.0D4
      IF(Y_KERN.LT.RL_LL(1)) Y_KERN=RL_LL(1)
      IF(Y_KERN.GT.RL_LL(K0_LL)) Y_KERN=RL_LL(K0_LL)
      CTURBLL(I,J)=F(X_KERN,Y_KERN,RL_LL,RL_LL,CTURB_LL &
                     ,K0_LL,K0_LL) 
   ENDDO
ENDDO
  RL_GL(1) = RADXXO(1,1)*1.E4 
  RL_GL(2) = 8.0D0
  RL_GL(3) = 10.0D0
  RL_GL(4) = 16.0D0
  RL_GL(5) = 20.0D0
  RL_GL(6) = 30.0D0
  RL_GL(7) = 40.0D0
  RL_GL(8) = 50.0D0
  RL_GL(9) = 60.0D0
  RL_GL(10)= 70.0D0
  RL_GL(11)= 80.0D0
  RL_GL(12)= 90.0D0
  RL_GL(13)=100.0D0
  RL_GL(14)=200.0D0
  RL_GL(15)=300.0D0
  RL_GL(16)=RADXXO(24,1)*1.0D4
! TURBULENCE GRAUPEL BULK RADII IN MKM
  RG_GL(1) = RADXXO(1,6)*1.0D4 
  RG_GL(2) = 30.0D0  
  RG_GL(3) = 60.0D0 
  RG_GL(4) = 100.0D0 
  RG_GL(5) = 200.0D0 
  RG_GL(6) = 300.0D0
  RG_GL(7) = 400.0D0
  RG_GL(8) = 500.0D0
  RG_GL(9) = 600.0D0
  RG_GL(10)= 700.0D0
  RG_GL(11)= 800.0D0
  RG_GL(12)= 900.0D0
  RG_GL(13)=1000.0D0
  RG_GL(14)=2000.0D0
  RG_GL(15)=3000.0D0
  RG_GL(16)=RADXXO(33,6)*1.0D4
DO I=KRMING_GL,KRMAXG_GL
   DO J=KRMINL_GL,KRMAXL_GL
      CTURBGL(I,J)=1.0D0
   ENDDO
ENDDO
DO I=1,K0G_GL
   DO J=1,K0L_GL
      CTURB_GL(I,J)=1.0D0
   ENDDO
ENDDO 
IF(IEPS_400.EQ.1) THEN
    CTURB_GL(1,1)=0.0D0
    CTURB_GL(1,2)=0.0D0
    CTURB_GL(1,3)=1.2D0
    CTURB_GL(1,4)=1.3D0
    CTURB_GL(1,5)=1.4D0
    CTURB_GL(1,6)=1.5D0
    CTURB_GL(1,7)=1.5D0
    CTURB_GL(1,8)=1.5D0
    CTURB_GL(1,9)=1.5D0
    CTURB_GL(1,10)=1.5D0
    CTURB_GL(1,11)=1.5D0
    CTURB_GL(1,12)=1.0D0
    CTURB_GL(1,13)=1.0D0
    CTURB_GL(1,14)=1.0D0
    CTURB_GL(1,15)=1.0D0

    CTURB_GL(2,1)=1.0D0
    CTURB_GL(2,2)=1.4D0
    CTURB_GL(2,3)=1.8D0
    CTURB_GL(2,4)=2.2D0
    CTURB_GL(2,5)=2.6D0
    CTURB_GL(2,6)=3.0D0
    CTURB_GL(2,7)=2.85D0
    CTURB_GL(2,8)=2.7D0
    CTURB_GL(2,9)=2.55D0
    CTURB_GL(2,10)=2.4D0
    CTURB_GL(2,11)=2.25D0
    CTURB_GL(2,12)=1.0D0
    CTURB_GL(2,13)=1.0D0
    CTURB_GL(2,14)=1.0D0

    CTURB_GL(3,1)=7.5D0
    CTURB_GL(3,2)=7.5D0
    CTURB_GL(3,3)=4.5D0
    CTURB_GL(3,4)=4.5D0
    CTURB_GL(3,5)=4.65D0
    CTURB_GL(3,6)=4.65D0
    CTURB_GL(3,7)=4.5D0
    CTURB_GL(3,8)=4.5D0
    CTURB_GL(3,9)=4.0D0
    CTURB_GL(3,10)=3.0D0
    CTURB_GL(3,11)=2.0D0
    CTURB_GL(3,12)=1.5D0
    CTURB_GL(3,13)=1.3D0
    CTURB_GL(3,14)=1.0D0
    
    CTURB_GL(4,1)=5.5D0
    CTURB_GL(4,2)=5.5D0
    CTURB_GL(4,3)=4.5D0
    CTURB_GL(4,4)=4.5D0
    CTURB_GL(4,5)=4.65D0
    CTURB_GL(4,6)=4.65D0
    CTURB_GL(4,7)=4.5D0
    CTURB_GL(4,8)=4.5D0
    CTURB_GL(4,9)=4.0D0
    CTURB_GL(4,10)=3.0D0
    CTURB_GL(4,11)=2.0D0
    CTURB_GL(4,12)=1.5D0
    CTURB_GL(4,13)=1.35D0
    CTURB_GL(4,14)=1.0D0
 
    CTURB_GL(5,1)=4.5D0
    CTURB_GL(5,2)=4.5D0
    CTURB_GL(5,3)=3.3D0
    CTURB_GL(5,4)=3.3D0
    CTURB_GL(5,5)=3.3D0
    CTURB_GL(5,6)=3.4D0
    CTURB_GL(5,7)=3.8D0
    CTURB_GL(5,8)=3.8D0
    CTURB_GL(5,9)=3.8D0
    CTURB_GL(5,10)=3.6D0
    CTURB_GL(5,11)=2.5D0
    CTURB_GL(5,12)=2.0D0
    CTURB_GL(5,13)=1.4D0
    CTURB_GL(5,14)=1.0D0
 
    CTURB_GL(6,1)=4.0D0
    CTURB_GL(6,2)=4.0D0
    CTURB_GL(6,3)=2.8D0
    CTURB_GL(6,4)=2.8D0
    CTURB_GL(6,5)=2.85D0
    CTURB_GL(6,6)=2.9D0
    CTURB_GL(6,7)=3.0D0
    CTURB_GL(6,8)=3.1D0
    CTURB_GL(6,9)=2.9D0
    CTURB_GL(6,10)=2.6D0
    CTURB_GL(6,11)=2.5D0
    CTURB_GL(6,12)=2.0D0
    CTURB_GL(6,13)=1.3D0
    CTURB_GL(6,14)=1.1D0

    CTURB_GL(7,1)=3.5D0
    CTURB_GL(7,2)=3.5D0
    CTURB_GL(7,3)=2.5D0
    CTURB_GL(7,4)=2.5D0
    CTURB_GL(7,5)=2.6D0
    CTURB_GL(7,6)=2.7D0
    CTURB_GL(7,7)=2.8D0
    CTURB_GL(7,8)=2.8D0
    CTURB_GL(7,9)=2.8D0
    CTURB_GL(7,10)=2.6D0
    CTURB_GL(7,11)=2.3D0
    CTURB_GL(7,12)=2.0D0
    CTURB_GL(7,13)=1.3D0
    CTURB_GL(7,14)=1.1D0

    CTURB_GL(8,1)=3.25D0
    CTURB_GL(8,2)=3.25D0
    CTURB_GL(8,3)=2.3D0
    CTURB_GL(8,4)=2.3D0
    CTURB_GL(8,5)=2.35D0
    CTURB_GL(8,6)=2.37D0
    CTURB_GL(8,7)=2.55D0
    CTURB_GL(8,8)=2.55D0
    CTURB_GL(8,9)=2.55D0
    CTURB_GL(8,10)=2.3D0
    CTURB_GL(8,11)=2.1D0
    CTURB_GL(8,12)=1.9D0
    CTURB_GL(8,13)=1.3D0
    CTURB_GL(8,14)=1.1D0

    CTURB_GL(9,1)=3.0D0
    CTURB_GL(9,2)=3.0D0
    CTURB_GL(9,3)=3.1D0
    CTURB_GL(9,4)=2.2D0
    CTURB_GL(9,5)=2.2D0
    CTURB_GL(9,6)=2.2D0
    CTURB_GL(9,7)=2.3D0
    CTURB_GL(9,8)=2.3D0
    CTURB_GL(9,9)=2.5D0
    CTURB_GL(9,10)=2.5D0
    CTURB_GL(9,11)=2.2D0
    CTURB_GL(9,12)=1.8D0
    CTURB_GL(9,13)=1.25D0
    CTURB_GL(9,14)=1.1D0

    CTURB_GL(10,1)=2.75D0
    CTURB_GL(10,2)=2.75D0
    CTURB_GL(10,3)=2.0D0
    CTURB_GL(10,4)=2.0D0
    CTURB_GL(10,5)=2.0D0
    CTURB_GL(10,6)=2.1D0
    CTURB_GL(10,7)=2.2D0
    CTURB_GL(10,8)=2.2D0
    CTURB_GL(10,9)=2.3D0
    CTURB_GL(10,10)=2.3D0
    CTURB_GL(10,11)=2.3D0
    CTURB_GL(10,12)=1.8D0
    CTURB_GL(10,13)=1.2D0
    CTURB_GL(10,14)=1.1D0

    CTURB_GL(11,1)=2.6D0
    CTURB_GL(11,2)=2.6D0
    CTURB_GL(11,3)=1.95D0
    CTURB_GL(11,4)=1.95D0
    CTURB_GL(11,5)=1.95D0
    CTURB_GL(11,6)=2.05D0
    CTURB_GL(11,7)=2.15D0
    CTURB_GL(11,8)=2.15D0
    CTURB_GL(11,9)=2.25D0
    CTURB_GL(11,10)=2.25D0
    CTURB_GL(11,11)=1.9D0
    CTURB_GL(11,12)=1.8D0
    CTURB_GL(11,13)=1.2D0
    CTURB_GL(11,14)=1.1D0

    CTURB_GL(12,1)=2.4D0
    CTURB_GL(12,2)=2.4D0
    CTURB_GL(12,3)=1.85D0
    CTURB_GL(12,4)=1.85D0
    CTURB_GL(12,5)=1.85D0
    CTURB_GL(12,6)=1.75D0
    CTURB_GL(12,7)=1.85D0
    CTURB_GL(12,8)=1.85D0
    CTURB_GL(12,9)=2.1D0
    CTURB_GL(12,10)=2.1D0
    CTURB_GL(12,11)=1.9D0
    CTURB_GL(12,12)=1.8D0 
    CTURB_GL(12,13)=1.3D0
    CTURB_GL(12,14)=1.1D0

    CTURB_GL(13,1)=1.67D0
    CTURB_GL(13,2)=1.67D0
    CTURB_GL(13,3)=1.75D0
    CTURB_GL(13,4)=1.83D0
    CTURB_GL(13,5)=1.87D0
    CTURB_GL(13,6)=2.0D0
    CTURB_GL(13,7)=2.1D0
    CTURB_GL(13,8)=2.12D0
    CTURB_GL(13,9)=2.15D0
    CTURB_GL(13,10)=2.18D0
    CTURB_GL(13,11)=2.19D0
    CTURB_GL(13,12)=1.67D0
    CTURB_GL(13,13)=1.28D0
    CTURB_GL(13,14)=1.0D0

    CTURB_GL(14,1)=1.3D0
    CTURB_GL(14,2)=1.3D0
    CTURB_GL(14,3)=1.35D0
    CTURB_GL(14,4)=1.4D0
    CTURB_GL(14,5)=1.6D0
    CTURB_GL(14,6)=1.7D0
    CTURB_GL(14,7)=1.7D0
    CTURB_GL(14,8)=1.7D0
    CTURB_GL(14,9)=1.7D0
    CTURB_GL(14,10)=1.7D0
    CTURB_GL(14,11)=1.7D0
    CTURB_GL(14,12)=1.4D0
    CTURB_GL(14,13)=1.25D0
    CTURB_GL(14,14)=1.0D0

    CTURB_GL(15,1)=1.17D0
    CTURB_GL(15,2)=1.17D0
    CTURB_GL(15,3)=1.17D0
    CTURB_GL(15,4)=1.25D0
    CTURB_GL(15,5)=1.3D0
    CTURB_GL(15,6)=1.35D0
    CTURB_GL(15,7)=1.4D0
    CTURB_GL(15,8)=1.4D0
    CTURB_GL(15,9)=1.45D0
    CTURB_GL(15,10)=1.47D0
    CTURB_GL(15,11)=1.44D0
    CTURB_GL(15,12)=1.3D0
    CTURB_GL(15,13)=1.12D0
    CTURB_GL(15,14)=1.0D0

    CTURB_GL(16,1)=1.17D0
    CTURB_GL(16,2)=1.17D0
    CTURB_GL(16,3)=1.17D0
    CTURB_GL(16,4)=1.25D0
    CTURB_GL(16,5)=1.3D0
    CTURB_GL(16,6)=1.35D0
    CTURB_GL(16,7)=1.4D0
    CTURB_GL(16,8)=1.45D0
    CTURB_GL(16,9)=1.45D0
    CTURB_GL(16,10)=1.47D0
    CTURB_GL(16,11)=1.44D0
    CTURB_GL(16,12)=1.3D0
    CTURB_GL(16,13)=1.12D0
    CTURB_GL(16,14)=1.0D0
ENDIF
IF(IEPS_800.EQ.1) THEN
    CTURB_GL(1,1) =0.00D0
    CTURB_GL(1,2) =0.00D0
    CTURB_GL(1,3) =1.00D0
    CTURB_GL(1,4) =1.50D0
    CTURB_GL(1,5) =1.40D0
    CTURB_GL(1,6) =1.30D0
    CTURB_GL(1,7) =1.20D0
    CTURB_GL(1,8) =1.10D0
    CTURB_GL(1,9) =1.00D0
    CTURB_GL(1,10)=1.00D0
    CTURB_GL(1,11)=1.00D0
    CTURB_GL(1,12)=1.00D0
    CTURB_GL(1,13)=1.00D0
    CTURB_GL(1,14)=1.00D0
    CTURB_GL(1,15)=1.00D0
    CTURB_GL(1,16)=1.00D0

    CTURB_GL(2,1) =0.00D0
    CTURB_GL(2,2) =0.00D0
    CTURB_GL(2,3) =1.00D0
    CTURB_GL(2,4) =2.00D0
    CTURB_GL(2,5) =1.80D0
    CTURB_GL(2,6) =1.70D0
    CTURB_GL(2,7) =1.60D0
    CTURB_GL(2,8) =1.50D0
    CTURB_GL(2,9) =1.50D0
    CTURB_GL(2,10)=1.50D0
    CTURB_GL(2,11)=1.50D0
    CTURB_GL(2,12)=1.50D0
    CTURB_GL(2,13)=1.50D0
    CTURB_GL(2,14)=1.00D0
    CTURB_GL(2,15)=1.00D0
    CTURB_GL(2,16)=1.00D0

    CTURB_GL(3,1) =0.00D0
    CTURB_GL(3,2) =0.00D0
    CTURB_GL(3,3) =4.00D0
    CTURB_GL(3,4) =7.65D0
    CTURB_GL(3,5) =7.65D0
    CTURB_GL(3,6) =8.00D0
    CTURB_GL(3,7) =8.00D0
    CTURB_GL(3,8) =7.50D0
    CTURB_GL(3,9) =6.50D0
    CTURB_GL(3,10)=6.00D0
    CTURB_GL(3,11)=5.00D0
    CTURB_GL(3,12)=4.50D0
    CTURB_GL(3,13)=4.00D0
    CTURB_GL(3,14)=2.00D0
    CTURB_GL(3,15)=1.30D0
    CTURB_GL(3,16)=1.00D0

    CTURB_GL(4,1) =7.50D0
    CTURB_GL(4,2) =7.50D0
    CTURB_GL(4,3) =7.50D0
    CTURB_GL(4,4) =7.65D0
    CTURB_GL(4,5) =7.65D0
    CTURB_GL(4,6) =8.00D0
    CTURB_GL(4,7) =8.00D0
    CTURB_GL(4,8) =7.50D0
    CTURB_GL(4,9) =6.50D0
    CTURB_GL(4,10)=6.00D0
    CTURB_GL(4,11)=5.00D0
    CTURB_GL(4,12)=4.50D0
    CTURB_GL(4,13)=4.00D0
    CTURB_GL(4,14)=2.00D0
    CTURB_GL(4,15)=1.30D0
    CTURB_GL(4,16)=1.00D0
    
    CTURB_GL(5,1) =5.50D0
    CTURB_GL(5,2) =5.50D0
    CTURB_GL(5,3) =5.50D0
    CTURB_GL(5,4) =5.75D0
    CTURB_GL(5,5) =5.75D0
    CTURB_GL(5,6) =6.00D0
    CTURB_GL(5,7) =6.25D0
    CTURB_GL(5,8) =6.17D0
    CTURB_GL(5,9) =5.75D0
    CTURB_GL(5,10)=5.25D0
    CTURB_GL(5,11)=4.75D0
    CTURB_GL(5,12)=4.25D0
    CTURB_GL(5,13)=4.00D0
    CTURB_GL(5,14)=2.00D0
    CTURB_GL(5,15)=1.35D0
    CTURB_GL(5,16)=1.00D0
 
    CTURB_GL(6,1) =4.50D0
    CTURB_GL(6,2) =4.50D0
    CTURB_GL(6,3) =4.50D0
    CTURB_GL(6,4) =4.75D0
    CTURB_GL(6,5) =4.75D0
    CTURB_GL(6,6) =5.00D0
    CTURB_GL(6,7) =5.25D0
    CTURB_GL(6,8) =5.25D0
    CTURB_GL(6,9) =5.00D0
    CTURB_GL(6,10)=4.75D0
    CTURB_GL(6,11)=4.50D0
    CTURB_GL(6,12)=4.00D0
    CTURB_GL(6,13)=3.75D0
    CTURB_GL(6,14)=2.00D0
    CTURB_GL(6,15)=1.40D0
    CTURB_GL(6,16)=1.00D0
 
    CTURB_GL(7,1) =4.00D0
    CTURB_GL(7,2) =4.00D0
    CTURB_GL(7,3) =4.00D0
    CTURB_GL(7,4) =4.00D0
    CTURB_GL(7,5) =4.00D0
    CTURB_GL(7,6) =4.25D0
    CTURB_GL(7,7) =4.50D0
    CTURB_GL(7,8) =4.67D0
    CTURB_GL(7,9) =4.50D0
    CTURB_GL(7,10)=4.30D0
    CTURB_GL(7,11)=4.10D0
    CTURB_GL(7,12)=3.80D0
    CTURB_GL(7,13)=3.50D0
    CTURB_GL(7,14)=2.00D0
    CTURB_GL(7,15)=1.30D0
    CTURB_GL(7,16)=1.10D0

    CTURB_GL(8,1) =3.50D0
    CTURB_GL(8,2) =3.50D0
    CTURB_GL(8,3) =3.50D0
    CTURB_GL(8,4) =3.65D0
    CTURB_GL(8,5) =3.65D0
    CTURB_GL(8,6) =3.80D0
    CTURB_GL(8,7) =4.1D02
    CTURB_GL(8,8) =4.17D0
    CTURB_GL(8,9) =4.17D0
    CTURB_GL(8,10)=4.00D0
    CTURB_GL(8,11)=3.80D0
    CTURB_GL(8,12)=3.67D0
    CTURB_GL(8,13)=3.40D0
    CTURB_GL(8,14)=2.00D0
    CTURB_GL(8,15)=1.30D0
    CTURB_GL(8,16)=1.10D0

    CTURB_GL(9,1) =3.25D0
    CTURB_GL(9,2) =3.25D0
    CTURB_GL(9,3) =3.25D0
    CTURB_GL(9,4) =3.25D0
    CTURB_GL(9,5) =3.25D0
    CTURB_GL(9,6) =3.50D0
    CTURB_GL(9,7) =3.75D0
    CTURB_GL(9,8) =3.75D0
    CTURB_GL(9,9) =3.75D0
    CTURB_GL(9,10)=3.75D0
    CTURB_GL(9,11)=3.60D0
    CTURB_GL(9,12)=3.40D0
    CTURB_GL(9,13)=3.25D0
    CTURB_GL(9,14)=2.00D0
    CTURB_GL(9,15)=1.30D0
    CTURB_GL(9,16)=1.10D0
    
    CTURB_GL(10,1) =3.00D0
    CTURB_GL(10,2) =3.00D0
    CTURB_GL(10,3) =3.00D0
    CTURB_GL(10,4) =3.10D0
    CTURB_GL(10,5) =3.10D0
    CTURB_GL(10,6) =3.25D0
    CTURB_GL(10,7) =3.40D0
    CTURB_GL(10,8) =3.50D0
    CTURB_GL(10,9) =3.50D0
    CTURB_GL(10,10)=3.50D0
    CTURB_GL(10,11)=3.40D0
    CTURB_GL(10,12)=3.25D0
    CTURB_GL(10,13)=3.15D0
    CTURB_GL(10,14)=1.90D0
    CTURB_GL(10,15)=1.30D0
    CTURB_GL(10,16)=1.10D0

    CTURB_GL(11,1) =2.75D0
    CTURB_GL(11,2) =2.75D0
    CTURB_GL(11,3) =2.75D0
    CTURB_GL(11,4) =2.75D0
    CTURB_GL(11,5) =2.75D0
    CTURB_GL(11,6) =3.00D0
    CTURB_GL(11,7) =3.25D0
    CTURB_GL(11,8) =3.25D0
    CTURB_GL(11,9) =3.25D0
    CTURB_GL(11,10)=3.25D0
    CTURB_GL(11,11)=3.25D0
    CTURB_GL(11,12)=3.15D0
    CTURB_GL(11,13)=3.00D0
    CTURB_GL(11,14)=1.80D0
    CTURB_GL(11,15)=1.30D0
    CTURB_GL(11,16)=1.10D0

    CTURB_GL(12,1) =2.60D0
    CTURB_GL(12,2) =2.60D0
    CTURB_GL(12,3) =2.60D0
    CTURB_GL(12,4) =2.67D0
    CTURB_GL(12,5) =2.67D0
    CTURB_GL(12,6) =2.75D0
    CTURB_GL(12,7) =3.00D0
    CTURB_GL(12,8) =3.17D0
    CTURB_GL(12,9) =3.17D0
    CTURB_GL(12,10)=3.17D0
    CTURB_GL(12,11)=3.10D0
    CTURB_GL(12,12)=2.90D0
    CTURB_GL(12,13)=2.80D0
    CTURB_GL(12,14)=1.87D0
    CTURB_GL(12,15)=1.37D0
    CTURB_GL(12,16)=1.10D0

    CTURB_GL(13,1) =2.40D0
    CTURB_GL(13,2) =2.40D0
    CTURB_GL(13,3) =2.40D0
    CTURB_GL(13,4) =2.50D0
    CTURB_GL(13,5) =2.50D0
    CTURB_GL(13,6) =2.67D0
    CTURB_GL(13,7) =2.83D0
    CTURB_GL(13,8) =2.90D0
    CTURB_GL(13,9) =3.00D0
    CTURB_GL(13,10)=2.90D0
    CTURB_GL(13,11)=2.85D0
    CTURB_GL(13,12)=2.80D0
    CTURB_GL(13,13)=2.75D0
    CTURB_GL(13,14)=1.83D0
    CTURB_GL(13,15)=1.30D0
    CTURB_GL(13,16)=1.10D0

    CTURB_GL(14,1) =1.67D0
    CTURB_GL(14,2) =1.67D0
    CTURB_GL(14,3) =1.67D0
    CTURB_GL(14,4) =1.75D0
    CTURB_GL(14,5) =1.75D0
    CTURB_GL(14,6) =1.83D0
    CTURB_GL(14,7) =1.87D0
    CTURB_GL(14,8) =2.00D0
    CTURB_GL(14,9) =2.10D0
    CTURB_GL(14,10)=2.12D0
    CTURB_GL(14,11)=2.15D0
    CTURB_GL(14,12)=2.18D0
    CTURB_GL(14,13)=2.19D0
    CTURB_GL(14,14)=1.67D0
    CTURB_GL(14,15)=1.28D0
    CTURB_GL(14,16)=1.00D0

    CTURB_GL(15,1) =1.30D0
    CTURB_GL(15,2) =1.30D0
    CTURB_GL(15,3) =1.30D0
    CTURB_GL(15,4) =1.35D0
    CTURB_GL(15,5) =1.35D0
    CTURB_GL(15,6) =1.40D0
    CTURB_GL(15,7) =1.60D0
    CTURB_GL(15,8) =1.70D0
    CTURB_GL(15,9) =1.70D0
    CTURB_GL(15,10)=1.70D0
    CTURB_GL(15,11)=1.70D0
    CTURB_GL(15,12)=1.70D0
    CTURB_GL(15,13)=1.70D0
    CTURB_GL(15,14)=1.40D0
    CTURB_GL(15,15)=1.25D0
    CTURB_GL(15,16)=1.00D0

    CTURB_GL(16,1) =1.17D0
    CTURB_GL(16,2) =1.17D0
    CTURB_GL(16,3) =1.17D0
    CTURB_GL(16,4) =1.17D0
    CTURB_GL(16,5) =1.17D0
    CTURB_GL(16,6) =1.25D0
    CTURB_GL(16,7) =1.30D0
    CTURB_GL(16,8) =1.35D0
    CTURB_GL(16,9) =1.40D0
    CTURB_GL(16,10)=1.45D0
    CTURB_GL(16,11)=1.45D0
    CTURB_GL(16,12)=1.47D0
    CTURB_GL(16,13)=1.44D0
    CTURB_GL(16,14)=1.30D0
    CTURB_GL(16,15)=1.12D0
    CTURB_GL(16,16)=1.00D0
ENDIF
IF(IEPS_800.EQ.1.AND.IEPS_1600.EQ.1) THEN
   DO I=1,K0G_GL
      DO J=1,K0L_GL
         CTURB_GL(I,J)=CTURB_GL(I,J)*1.7D0
      ENDDO
    ENDDO 
ENDIF
DO J=1,K0L_GL
   DO I=1,K0G_GL
      CTURB_GL(I,J)=(CTURB_GL(I,J)-1.0D0)/1.5D0+1.0D0
   ENDDO
ENDDO
DO I=KRMING_GL,KRMAXG_GL
   DO J=KRMINL_GL,KRMAXL_GL
      CTURBGL(I,J)=1.
   ENDDO
ENDDO
DO I=KRMING_GL,KRMAXG_GL 
   X_KERN=RADXXO(I,6)*1.0D4
   IF(X_KERN.LT.RG_GL(1)) X_KERN=RG_GL(1)
   IF(X_KERN.GT.RG_GL(K0G_GL)) X_KERN=RG_GL(K0G_GL) 
   DO J=KRMINL_GL,KRMAXL_GL
      Y_KERN=RADXXO(J,1)*1.0D4
      IF(Y_KERN.LT.RL_GL(1)) Y_KERN=RL_GL(1)
      IF(Y_KERN.GT.RL_GL(K0L_GL)) Y_KERN=RL_GL(K0L_GL)
      CTURBGL(I,J)=F(X_KERN,Y_KERN,RG_GL,RL_GL,CTURB_GL &
                     ,K0G_GL,K0L_GL)
   ENDDO
ENDDO
IF(IEPS_800.EQ.1) THEN
   DO I=KRMING_GL,15
      DO J=KRMINL_GL,13
         IF(CTURBGL(I,J).LT.3.0D0) CTURBGL(I,J)=3.0D0
      ENDDO
   ENDDO
ENDIF
IF(IEPS_1600.EQ.1) THEN
   DO I=KRMING_GL,15
      DO J=KRMINL_GL,13
         IF(CTURBGL(I,J).LT.5.1D0) CTURBGL(I,J)=5.1D0
      ENDDO
   ENDDO
ENDIF
DO I=1,33
   DO J=1,24
      IF(I.LE.14.AND.J.EQ.8) CTURBGL(I,J)=1.0D0
      IF(I.GT.14.AND.J.LE.8) CTURBGL(I,J)=1.2D0
   ENDDO
ENDDO     
RETURN
END SUBROUTINE TURBCOEF

!##############################################################################
real Function f (x,y,x0,y0,table,k0,kk0)

! two-dimensional linear interpolation of the collision efficiency
! with help table(k0,kk0)

 implicit none
 integer k0,kk0,k,ir,kk,iq
 double precision x,y,p,q,ec,ek
 double precision x0(k0),y0(kk0),table(k0,kk0)

do k=2,k0
   if(x.le.x0(k).and.x.ge.x0(k-1)) then
      ir=k     
   elseif(x.gt.x0(k0)) then
      ir=k0+1
   elseif(x.lt.x0(1)) then
      ir=1
   endif
enddo
do kk=2,kk0
   if(y.le.y0(kk).and.y.ge.y0(kk-1)) iq=kk
enddo
if(ir.lt.k0+1) then
   if(ir.ge.2) then
      p =(x-x0(ir-1))/(x0(ir)-x0(ir-1))
      q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
      ec=(1.d0-p)*(1.d0-q)*table(ir-1,iq-1)+ &
             p*(1.d0-q)*table(ir,iq-1)+ &
             q*(1.d0-p)*table(ir-1,iq)+ &
                  p*q*table(ir,iq)    
    else
      q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
      ec=(1.d0-q)*table(1,iq-1)+q*table(1,iq)    
    endif
  else
    q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
    ek=(1.d0-q)*table(k0,iq-1)+q*table(k0,iq)
    ec=min(ek,1.d0) 
  endif
  f=ec

return
END FUNCTION f

!##############################################################################
Subroutine FREEZ (FF1,FF2,FF5,TIN,DT,RO,K,I,J) 

use mem_micro
use mem_grid, only:ngrid

IMPLICIT NONE 
INTEGER KR,ICE_TYPE
REAL DT,RO,YKK,PF,PF_1,DEL_T,ARG_1,DF1,BF,ARG_M, & 
     TT_DROP_AFTER_FREEZ,CFREEZ,TIN,TTIN,AF
REAL FF1(NKR),FF2(NKR,ICEMAX),FF5(NKR)
REAL SUM_NUM,SUM_ICE
INTEGER K,I,J

TTIN=TIN
DEL_T=TTIN-273.15
!Frozen cloud droplets sent to plates
ICE_TYPE=2

IF(DEL_T.LT.0.AND.any(FF1>0.)) THEN
   SUM_ICE=0.0
   SUM_NUM = 0.0 
   AF=AFREEZMY
   CFREEZ=(BFREEZMAX-BFREEZMY)/XL(NKR)
   DO KR=1,NKR
      ARG_M=XL(KR)
      BF=BFREEZMY+CFREEZ*ARG_M
      PF_1=AF*EXP(-BF*DEL_T)
      PF=ARG_M*PF_1
      YKK=EXP(-PF*DT)
      DF1=FF1(KR)*(1.-YKK)
      FF1(KR)=FF1(KR)*YKK
      IF(KR.LE.KRFREEZ)  THEN
         FF2(KR,ICE_TYPE)=FF2(KR,ICE_TYPE)+DF1
      ELSE
         FF5(KR)=FF5(KR)+DF1
      ENDIF
      SUM_ICE=SUM_ICE+DF1*XI(KR,1)*XI(KR,1)*COL3/RO
      SUM_NUM = SUM_NUM + DF1*XI(KR,1)*COL3/RO*1000.
      IF (imbudget >= 1) THEN
          micro_g(ngrid)%nucicert(k,i,j) = micro_g(ngrid)%nucicert(k,i,j) + SUM_ICE
          micro_g(ngrid)%nucicect(k,i,j) = micro_g(ngrid)%nucicect(k,i,j) + SUM_NUM
      END IF
      IF (imbudget >= 2) THEN
         micro_g(ngrid)%inuchomrt(k,i,j) = micro_g(ngrid)%inuchomrt(k,i,j) + SUM_ICE
         micro_g(ngrid)%inuchomct(k,i,j) = micro_g(ngrid)%inuchomct(k,i,j) + SUM_NUM
      END IF    
   ENDDO

   !NEW TEMPERATURE
   ARG_1=333.*SUM_ICE ! 333=334000/1004 [K]
   TT_DROP_AFTER_FREEZ=TTIN+ARG_1
   TIN=TT_DROP_AFTER_FREEZ
ENDIF

RETURN     
END SUBROUTINE FREEZ 

!##############################################################################
Subroutine MELT (FF1,FF2,FF3,FF4,FF5,TIN,RO,K,I,J)

use mem_micro
use mem_grid, only:ngrid

IMPLICIT NONE

INTEGER K,I,J
INTEGER KR,ICE
REAL ARG_M,ARG_1,RO,SUMICE,TIN
REAL FF1(NKR),FF2(NKR,ICEMAX),FF3(NKR),FF4(NKR),FF5(NKR)

IF (IMBUDGET >= 2) THEN !Everything melts - just add it all up
   IF (IPRIS >=2) THEN
      micro_g(ngrid)%meltprist(k,i,j) = micro_g(ngrid)%meltprist(k,i,j) &
                        + SUM_MASS(FF2(:,1),XI(:,1),RO,1,KRPRIS(1)-1) &
                        + SUM_MASS(FF2(:,2),XI(:,2),RO,1,KRPRIS(2)-1) &
                        + SUM_MASS(FF2(:,3),XI(:,3),RO,1,KRPRIS(3)-1)

      micro_g(ngrid)%meltsnowt(k,i,j) = micro_g(ngrid)%meltsnowt(k,i,j) &
                        + SUM_MASS(FF2(:,1),XI(:,1),RO,KRPRIS(1),NKR) &
                        + SUM_MASS(FF2(:,2),XI(:,2),RO,KRPRIS(2),NKR) &
                        + SUM_MASS(FF2(:,3),XI(:,3),RO,KRPRIS(3),NKR)

   ELSEIF (IPRIS==1 .or. IPRIS==2 .or. IPRIS==3) THEN
      micro_g(ngrid)%meltprist(k,i,j) = micro_g(ngrid)%meltprist(k,i,j) &
              + SUM_MASS(FF2(:,IPRIS),XI(:,IPRIS),RO,1,KRPRIS(IPRIS)-1)
      micro_g(ngrid)%meltsnowt(k,i,j) = micro_g(ngrid)%meltsnowt(k,i,j) &
              + SUM_MASS(FF2(:,IPRIS),XI(:,IPRIS),RO,KRPRIS(IPRIS),NKR)
   ENDIF
   micro_g(ngrid)%meltaggrt(k,i,j) = micro_g(ngrid)%meltaggrt(k,i,j) &
           + SUM_MASS(FF3,XS,RO,1,NKR)
   IF(IGRAUP>0) micro_g(ngrid)%meltgraut(k,i,j) = micro_g(ngrid)%meltgraut(k,i,j) &
                        + SUM_MASS(FF4,XG,RO,1,NKR)
   IF(IHAIL>0)  micro_g(ngrid)%melthailt(k,i,j) = micro_g(ngrid)%melthailt(k,i,j) &
                        + SUM_MASS(FF5,XH,RO,1,NKR)
ENDIF 

SUMICE=0.
DO KR=1,NKR
   ARG_M=FF3(KR)+FF4(KR)+FF5(KR)
   DO ICE=1,ICEMAX
      ARG_M=ARG_M+FF2(KR,ICE)
      FF2(KR,ICE)=0.
   ENDDO
   FF1(KR)=FF1(KR)+ARG_M
   FF3(KR)=0.
   FF4(KR)=0.
   FF5(KR)=0.
   SUMICE=SUMICE+ARG_M*3.*XL(KR)*XL(KR)*COL
ENDDO
! NEW TEMPERATURE :
ARG_1=333.*SUMICE/RO
TIN=TIN-ARG_1
IF (IMBUDGET >= 1) micro_g(ngrid)%melticet(k,i,j) = &
    micro_g(ngrid)%melticet(k,i,j) + SUMICE/RO

RETURN     
END SUBROUTINE MELT 

!##############################################################################
Subroutine JERNUCL01 (PSI1,PSI2,FCCNR,DTT,ROR,DSUP1,DSUP2,SUP2_OLD,DT,finr,K,I,J)

IMPLICIT NONE
 
INTEGER K,I,J
REAL SUP2_OLD, &
     FCCNR(NKR),FINR(NKR)
DOUBLE PRECISION DTT,DSUP1,DSUP2
REAL TT,ROR, &
     TPC,SUP1,SUP2,DEL1N,DEL2N,AL1,AL2, &
     TEMP1,TEMP2,TEMP3,A1,B1,A2,B2,DT
REAL PSI1(NKR),PSI2(NKR,ICEMAX)
REAL alwc
DATA A1,B1,A2,B2/-0.639,0.1296,-2.8,0.262/
DATA TEMP1,TEMP2,TEMP3/-5.,-2.,-20./
DATA AL1/2500./,AL2/2834./
SUP1=DSUP1
SUP2=DSUP2

TT=DTT

DEL1N=100.*SUP1
TPC=TT-273.15

IF(DEL1N.GT.0.AND.TPC.GT.-30.) THEN
   CALL WATER_NUCL (PSI1,FCCNR,XL,TT,SUP1  &
       ,RCCN,DROPRADII,NKR,ror,K,I,J)
ENDIF

IF (ICEPROCS.EQ.1)THEN
   DEL2N=100.*SUP2
   IF(TPC.LT.0..AND.TPC.GE.-35..AND.DEL2N.GT.0.) THEN
         CALL ICE_NUCL (PSI2,XI,TT,ROR,SUP2,SUP2_OLD,SUP1,DT,finr,K,I,J)
   ENDIF
ENDIF

RETURN
END SUBROUTINE JERNUCL01

!##############################################################################
Subroutine WATER_NUCL (PSI1,FCCNR,X1,TT,SUP1 &
      ,RCCN,DROPRADII,NKR,dens,K,I,J)

use mem_micro
use mem_grid, only:ngrid

IMPLICIT NONE

INTEGER K,I,J
INTEGER NDROPMAX,KR,NKR
REAL PSI1(NKR),FCCNR(NKR),X1(NKR)
REAL DROPCONCN(NKR)
REAL RCCN(NKR),DROPRADII(NKR)
REAL TT,SUP1,DX
real dens

CALL NUCLEATION (SUP1,TT,FCCNR,DROPCONCN  &
     ,NDROPMAX,RCCN,DROPRADII,NKR)

! NEW WATER SIZE DISTRIBUTION function (BEGIN)
  DO KR=1,NDROPMAX
     DX=3.*COL*X1(KR)
     PSI1(KR)=PSI1(KR)+DROPCONCN(KR)/DX

     IF (imbudget >= 1) THEN
        micro_g(ngrid)%nuccldrt(k,i,j) = micro_g(ngrid)%nuccldrt(k,i,j) &
            + dropconcn(kr)*x1(kr)/dens
        micro_g(ngrid)%nuccldct(k,i,j) = micro_g(ngrid)%nuccldct(k,i,j) &
            + dropconcn(kr)/dens*1000.
     ENDIF
  ENDDO

RETURN
END SUBROUTINE WATER_NUCL

!##############################################################################
Subroutine ICE_NUCL (PSI2,X2,TT,ROR,SUP2,SUP2_OLD &
                     ,SUP1,DT,finr,K,I,J)

!     ICE_NUCL Modifcation History
!     (April 2007 J. Comstock)
!     modified to include classical theory heterogeneous nucleation
!     via the condensation/immersion freezing mode
!     Added ICEFLAG=0 Use Meyers param, ICEFLAG=1 Use Classical Theory

use mem_micro
use mem_grid, only:ngrid

IMPLICIT NONE

  INTEGER K,I,J
  INTEGER ITYPE,KR,ICE,NRGI,K1,ki
  REAL DEL2N,SUP1,SUP2,C1,C2,TPC,TT,ROR
  REAL DX
  REAL HELEK1,HELEK2,TPCC,DEL2NN,FF1BN
  REAL FACT,DSUP2N,SUP2_OLD,DELTACD,DELTAF,ADDF
  REAL X2(NKR,ICEMAX),PSI2(NKR,ICEMAX)
  REAL DT,finr(NKR)
  REAL QHET(NKR),NUMHET(NKR)

  real rcn(nkr)
  real inic(nkr)

  REAL A1,B1,A2,B2
  DATA A1,B1,A2,B2/-0.639,0.1296,-2.8,0.262/
  REAL TEMP1,TEMP2,TEMP3
  DATA TEMP1,TEMP2,TEMP3/-5.,-2.,-20./

  C1=C1_MEY
  C2=C2_MEY
! TYPE OF ICE WITH NUCLEATION (BEGIN)

  TPC=TT-273.15
  ITYPE=0

  IF((TPC.GT.-4.0).OR.(TPC.LE.-8.1.AND.TPC.GT.-12.7).OR.&
     (TPC.LE.-17.8.AND.TPC.GT.-22.4)) THEN
    ITYPE=2
  ELSE
    IF((TPC.LE.-4.0.AND.TPC.GT.-8.1).OR.(TPC.LE.-22.4)) THEN
      ITYPE=1
    ELSE
      ITYPE=3
    ENDIF
  ENDIF

! NEW CRYSTAL SIZE DISTRIBUTION function (BEGIN)

  ICE=ITYPE

  IF (ICEFLAG .EQ. 0) THEN  !USE MEYERS ICE NUCLEATION SCHEME
     NRGI=2
     IF(TPC.LT.TEMP1) THEN
        DEL2N=100.*SUP2
        DEL2NN=DEL2N
        IF(DEL2N.GT.50.0) DEL2NN=50.
           HELEK1=C1*EXP(A1+B1*DEL2NN)
      ELSE
           HELEK1=0.
      ENDIF

      IF(TPC.LT.TEMP2) THEN
         TPCC=TPC
         IF(TPCC.LT.TEMP3) TPCC=TEMP3
            HELEK2=C2*EXP(A2-B2*TPCC)
      ELSE
         HELEK2=0.
      ENDIF

      FF1BN=HELEK1+HELEK2

      FACT=1.
      DSUP2N=(SUP2-SUP2_OLD)*100.

      SUP2_OLD=SUP2

      IF(DSUP2N.GT.50.) DSUP2N=50.

      DELTACD=FF1BN*B1*DSUP2N

      IF(DELTACD.GE.FF1BN) DELTACD=FF1BN
      IF(DELTACD.GT.0.) THEN
         DELTAF=DELTACD*FACT
         DO KR=1,NRGI-1
            DX=X2(KR,ICE)*COL3
            ADDF=DELTAF/DX
            PSI2(KR,ICE)=PSI2(KR,ICE)+ADDF
            IF (imbudget >=1) THEN
               micro_g(ngrid)%nucicert(k,i,j) = micro_g(ngrid)%nucicert(k,i,j) &
                   + DELTAF*X2(KR,ICE)/ROR
               micro_g(ngrid)%nucicect(k,i,j) = micro_g(ngrid)%nucicect(k,i,j) &
                   + DELTAF/ROR*1000.
            ENDIF
            IF (imbudget >= 2) THEN
               micro_g(ngrid)%inucifnrt(k,i,j) = micro_g(ngrid)%inucifnrt(k,i,j) &
                   + DELTAF*X2(KR,ICE)/ROR
               micro_g(ngrid)%inucifnct(k,i,j) = micro_g(ngrid)%inucifnct(k,i,j) &
                   + DELTAF/ROR*1000.
            ENDIF
         ENDDO
      ENDIF
   ELSE   !USE CLASSICAL THEORY NUCLEATION SCHEME

!===================CCN-version ice nucleation (added by J. Fan)=======================

      DO KR=1,NKR
         rcn(kr) = rccn(kr)*1.e-2              ! transform to m from cm
      ENDDO
          
      DO KR=1,NKR
         inic(kr) = finr(kr)*col*1.e6  ! m-3
      ENDDO

      CALL heteronuc_ccn (TT,SUP1,DT,NKR,inic,rcn,QHET,NUMHET)
      ! Distribute nucleated aerosols to appropriate ice size bin (QHET:kg; NUMHET: m-3)

      ! Update IN
      DO KR=1, NKR
         finr(kr) = inic(kr)*1.e-6/col
      ENDDO

      if (sum(NUMHET(:)) > 0.) then
         do ki=1,NKR
            do kr=1,NKR
               if (QHET (ki)*1.e3 .LE. X2(1,ice) )  k1=1
               if (QHET (ki)*1.e3 .GE. X2(nkr,ice) )  k1=nkr
               if ((kr.GT.1) .and. (QHET (ki)*1.e3 .LE. X2(kr,ice)).and.  &
                   (QHET (ki)*1.e3 .GT. X2(kr-1,ice)))  k1=kr
            enddo
            DX=3.*X2(k1,ICE)*COL
            PSI2(K1,ICE) = PSI2(K1,ICE) + NUMHET(ki)*1.e-6/DX
            IF (imbudget >=1) THEN
               micro_g(ngrid)%nucicert(k,i,j) = micro_g(ngrid)%nucicert(k,i,j) &
                   + numhet(ki)/DX*COL3*X2(k1,ICE)*X2(k1,ICE)/ROR*1.e-6
               micro_g(ngrid)%nucicect(k,i,j) = micro_g(ngrid)%nucicect(k,i,j) &
                   + numhet(ki)/DX*COL3*X2(k1,ICE)/ROR*1.e-3
            ENDIF
            IF (imbudget >=2) THEN
               micro_g(ngrid)%inucifnrt(k,i,j) = micro_g(ngrid)%inucifnrt(k,i,j) &
                   + numhet(ki)/DX*COL3*X2(k1,ICE)*X2(k1,ICE)/ROR*1.e-6
               micro_g(ngrid)%inucifnct(k,i,j) = micro_g(ngrid)%inucifnct(k,i,j) &
                   + numhet(ki)/DX*COL3*X2(k1,ICE)/ROR*1.e-3
            ENDIF
              
         enddo
      endif
   ENDIF

RETURN
END SUBROUTINE ICE_NUCL

!##############################################################################
Subroutine heteronuc_ccn (temp,Sw,dt,INBINMAX,icenuclei, &
           draeros,qhet,numhet)

! Computes heterogeneous freezing rate, crystal number and mass assuming
! immersion freezing from ccn distribution (assuming ccn in aqueous soln)
! Assumes nucleation is immersion/condensation freezing from solution drops+ccn

!     This is main subroutine for heterogeneous nucleation
!     described in Morrison et al. 2005 JAS.

!     modified to interface with HUCM Spectral Bin Model (18 Apr 2007 JMC)

!     Germ formation rates due to condensation freezing are described in
!     Khvorostyanov and Curry 2005 JAS.
!     This formulation passes the aerosol (IN) size distribution into subroutine
!     and tracks IN loss due to nucleation scavenging.

!     modified to pass the drop size distribution rather than the IN dry aerosol
!     distribution and removed aerosol growth. (Oct 2007 JMC)

implicit none

!     input:
integer INBINMAX    ! number of aerosol (IN) bins
real temp,        &  ! temperature (K)
    Sw  ,        &  ! water supersaturation ratio (dimensionless)
    dt  ,        &  ! time step (s)
    icenuclei(INBINMAX),& ! ice nuclei (from ccn distribution) (m-3)
    draeros(INBINMAX)    ! dry aerosol radius (m)

!     output:
!     note: numhet and qhet are the total number aerosol nucleated during time step,
!     not nucleation rate!!!
real numhet(INBINMAX)     ! number of aerosol nucleated for aerosol bin (m-3)
real qhet(INBINMAX)       ! mass of single ice crystal nucleated per aerosol bin (kg)

!     internal:
real raeros,          &    ! radius of insoluble substrate (m)
    wraeros,         &    ! wet aerosol radius (m)
    ajlsfrz,         &    ! het frz rate for particle
    probhet,         &    ! probability of het frz for particle
    Tc,              &    ! temperature (C)
    hrel                 ! saturation ratio (unitless)
integer kr                ! counter for aerosol size bin

!     other parameter
real rhow,             &   !density of water
     rhoi                  !density of ice
real pi

Tc = temp-273.15          !! K to C
hrel = Sw + 1.0           ! saturation ratio

!     define constants
pi = 4.*atan(1.)

!      qvaero = 0.85   !aerosol solubility in terms of volume fraction
!                      this quantity is assumed, may want to upgrade formulation
!     wettability
!      wetcoef = 0.9
!     relative area of active sites
!      alf = 0.5e-5
!     misfit strain parameter
!      epsil = 2.5e-2            ! e=2.5% Turnbull vonegut, 1952
!      epsil = 0.1e-2            !! e~0; Turnbull vonegut, 1952
!      epsil=0.01  !e=1%
!     density of water (kg/m3)
rhow = 1000.0
!     density of ice (kg/m3)
rhoi = 900.0
!     aerosol parameter describing composition (Khvorostyanov and Curry 1999)
!      betafr = 0.5

!-----------------------------------------------------------------------
!     numerical integration of size bins

!     initialize nucleation number and mass

do kr = 1,INBINMAX
   numhet(kr) = 0.
   qhet(kr) = 0.
end do

!     main loop around aerosol size bins, starting with largest bin

do kr=INBINMAX,1,-1

!     determine radius of insoluble portion of aerosol raeros
!     size of insoluble substrate determines germ formation rate
!     See Khvorostyanov and Curry 2004

   raeros = draeros(kr) * (1.-qvaero)**(1./3.)  !m

!   print*,'raeros=',raeros
!     call to get germ formation rate

   if (raeros.gt.0.) then
      CALL jhetfrz2 (Tc,hrel,raeros,ajlsfrz)
   else
      ajlsfrz = 0.
   end if

!     calculate probability of aerosol particle freezing within dt

   probhet = 1.-exp(-ajlsfrz*dt)
!   if (probhet > 1.e-3) print *,'probhet',probhet,ajlsfrz,raeros

!     print *,probhet,ajlsfrz,raeros

!     if there is no freezing associated with largest size bin, then
!     exit loop, there will be no nucleation
!     use probability of 1.e-10 as cutoff

   if (kr.eq.INBINMAX) then
      if (probhet.lt.1.e-10) return
   end if

!     nucleate if probability is significant

  if (probhet.ge.1.e-10) then

!     number of ice nucleated for each drop size bin (1/m3)

      numhet(kr) = probhet * AMAX1(icenuclei(kr),0.0)
!    print*,'numhet',kr,numhet(kr),ajlsfrz,probhet,icenuclei(kr)

!----------------------------------------------------------------------------------
!     note, use wetted aerosol for radius to calculate mass
!     of nucleated crystals for size bin (kg m-3)
!     this is done by multiplying nucleation rate (m-3) by
!     4/3*pi*rhoi*r^3
!     rhoi = 900 kg m-3 is assumed bulk density of ice

      CALL req_haze (temp,Sw,draeros(kr),wraeros)
!      print*, 'wet rad', kr, draeros(kr), wraeros

      qhet(kr) = 4./3.*pi*rhoi*(wraeros)**3.0 !mass of a single particle in (kg)

      icenuclei(kr) = icenuclei(kr) - numhet(kr)

   end if                 !! significant freezing probability
end do                    !! end loop around size bins

return
END SUBROUTINE heteronuc_ccn

!##############################################################################
Subroutine jhetfrz2 (ttc,relh,raeros,ajlshet)

!     this subroutine calculates germ formation rates on insoluble subsrate for
!     given T, S, etc.
!     Khvorostyanov and Curry, 2005 JAS
!     All units are mks (JMC) 18 Apr 2007

implicit none

!     input data

real ttc,     &      ! temperature (C!!!!)
     relh,    &      ! saturation ratio over water
     sig_iw0, &      !surface tension ice/water interface
     raeros          !radius of the insolouble substrate (m)

!     output data

real  ajlshet       !nucleation rate (s-1)
real  sw_th    !threshold saturation

!     internal variables
!!!!!!! Modified by Yi Wang !!!!!!!
!!!! Using c1s (c,one,s) is same as cls (c,l,s)
!!!! Change c1s to cleans

REAL cturn,rhoimean,akbol,rhoi,akplank,cleans,almmeanc  &
          ,a_lm0,a_lm1,a_lm2,ttze,a_lmtn(4),rlog,cor1,cor2   &
          ,almmean,ttk,tt0,humf,expind,sw_thresh,fgerms,Rv
REAL factivh,xxx,fisize,fisizeSS,xmfi,ccc1,ccc2,ccc3, &
     ffshape2,fshape,fgerms0,dfgerms,coefnucl,rgerms

!     Turnb-Vonneg strain const, [NT/m^-2], 24 June 2000, PK97, p.343
CTurn=1.7e10
cleans=1.e19     ! [m^-2], conc. molec. per 1 m^2 of surface, PK97, p.342
              ! PK97, p. 206, was 10^18!!!  !m^-2, number of molec. per 1
              ! [m^2], contact. ice germ
RHOIMEAN=900.0      !kg/m3, mean ice density
RHOI=900.0          !kg/m3, ice density


akbol    =1.380622e-23    ![J/K], Boltzmann constant
akplank  =6.626176e-34    ![J*s], Planck constant
Rv       =4.615e2         !specific gas constant for water vapor (J/kg/K)

sig_iw0=(28.0+0.25*TTC)*1.e-3 !surface tension, [NT/m] ice-water interface


TTZE=273.15

A_LM0    =70.0            ! Lm=70 cal/g !Prup'95
A_LM1    =0.708
A_LM2    =-2.5e-3         !L''m, Prup.'78, p.89
ALMMEANC =79.7 + 0.693*TTC-2.5E-3*TTC**2 !old PK78, corr. 2 June 1998
A_LMTN(1)=79.7 + A_LM1*TTC+A_LM2*(TTC**2) ! Lm, Prup'78, p.89;
A_LMTN(2)=A_LM0 + A_LM1*TTC+A_LM2*TTC**2 ! Pr95, fig.8, PK97, fig. 3-12

!     below:  Temperature correction  fi(T) from
!     Integr. melt. heat  -temp. corr.

RLOG    =alog(TTZE/(TTC+TTZE))
cor1    =-A_LM1/A_LM0*(TTZE+TTC/RLOG)
cor2    = A_LM2/A_LM0*(TTZE**2.0-TTC/RLOG*(TTC-2.0*TTZE)/2.)

A_LMTN(4)=1.+cor1+cor2
A_LMTN(3)=A_LM0 * A_LMTN(4)

!     effective melting heat

ALMMEANC=A_LMTN(1)        ![kcal/kg] or [cal/g]
ALMMEAN=ALMMEANC*4.18e3   ![J/kg]

TTK=TTC+273.15            !Temperature Kelvin
TT0=273.15                !0 C (Kelvin=273.15)
!     param. G=Rv*T/Lm (divided by corr. Lm(TTC), 1 June 1998
humf=0.11*TTK/ALMMEANC


expind=CTurn*epsil**2.0/(RHOI*ALMMEAN)

!     thershold saturation

Sw_thresh=(TTK/TT0*exp(expind))**(1./humf)
Sw_th=Sw_thresh           !for output

!     if sw less than theshold then exit out, no nucleation
!      print*,' Sw_th=',Sw_th
if (Sw_thresh.GE.relh) then !no nucleation, go out
   AJLSHet=1.e-20
   FGERMS=1.e-14
   rgerms=1.e-08
   go to 98765
endif

!     activation energy

!******Activation energy at T=-30 C is 10 kkal/mole (PK97, Fig. 3-11, p.95);
!***   this is =0.694*10^-12 erg; starting point for this fit at T=-30 C.
!***   It decreases with T and coincides with Jensen94 at T=-90 C

FACTIVH= &           ! Activation energy, erg
          0.694E-12*(1.+0.027*(TTC+30.0)) !Linear fit to Prupp.'95  4 Nov 1997
if (TTC.gt.-30.0) FACTIVH=0.385e-12*  &
          exp(-8.423e-3*TTC+6.384e-4*TTC**2.0+7.891e-6*TTC**3.0) !p.96

if (TTC.le.-40.) FACTIVH=  & !for`low T<-40 C
          0.694E-12*(1.+0.027*(TTC+30.0)*exp(0.010*(TTC+30.0))) !MY CHOICE
FACTIVH = FACTIVH * 1.E-7  !convert from erg to Joules

!     rgerms in m
rgerms=2.0*sig_iw0/(ALMMEAN*RHOIMEAN*(ALOG(273.15/TTK)   &
          +humf*ALOG(relh))-CTurn*epsil**2.0) !corr 24 June 2000 for misfit strain

!     this check is to make sure that radius of germ is greater than 0

if (rgerms.gt.0.) go to 12345 !continue if rgerms>0

!*****************go out with AJLSHet=0, if rgerms<0  **********************

if (rgerms.le.0.) then
   AJLSHet=1.e-20         !no nucleation, go out
   FGERMS=1.e-24
   rgerms=1.e-18
!   print*,'rgerms.le.0.'
   go to 98765      !go out with AJLSHet=0, if rgerms<0
endif

12345 xxx=raeros/rgerms

!     shape factor

fisizeSS=(1.0-2.0*wetcoef*xxx+xxx**2.0)**0.5 !12 Nov 1998, PK97, p. 302
fisize=fisizeSS           !18 July 2000, derevn.
xmfi=(xxx-wetcoef)/fisizeSS

ccc1=((1-wetcoef*xxx)/fisizeSS)**3.0
ccc2=(xxx**3.0)*(2.0-3.0*xmfi+xmfi**3.0)
ccc3=3*wetcoef*(xxx**2)*(xmfi-1)
ffshape2=0.5*(1+ccc1+ccc2+ccc3)
fshape=ffshape2
!print*,'fshape=',fshape
FGERMS0=4./3.*3.1416*sig_iw0*(rgerms**2.0)*fshape !Germ energy with shape
!print*,'fgerms0=',fgerms0,sig_iw0,rgerms

!*****!Fletcher's correction for active site, PK97, p.345, !11 April 1999
DFGERMS=alf*(raeros**2.0)*sig_iw0*(1.-wetcoef) !11 April 1999, PK97, p.345
!print*,'dfgerms=',dfgerms

!     correction for active site, PK97, p.345; 24 June 2000

FGERMS=FGERMS0-DFGERMS
if (FGERMS.le.0.) FGERMS=0. ! no zero values ???

!     output to screen
!     if (FGERMS.le.0.) write (*,*) 'FGERMS LT 0.',' TTK=',TTK,
!     # ' relh=',relh,' FGERMS0=',FGERMS0,' DFGERMS=',DFGERMS,'alf=',alf,
!     # ' raeros=',raeros

!     preexpon. fact. PK97, p.342, heterogen. 12 Nov 1998

coefnucl=(akbol*TTK/akplank)*    & !12 Nov 1998, PK97, p. 342
          (4.0*3.1416*raeros**2.0)*cleans
!print*,'coefnucl=',coefnucl,akbol,ttk,akplank,raeros,cleans
!     calculate freezing rate

AJLSHet=           &     !s^-1, heterogen. salt nucleation rate, PK97, p. 342
          coefnucl*exp(-(FACTIVH+FGERMS)/(akbol*TTK)) !KS.98

if (AJLSHet.lt.1.e-20) AJLSHet=1.e-20

98765 return  !go here if rgerms<0
END SUBROUTINE jhetfrz2

!##############################################################################
Subroutine req_haze (temp,sw,draeros,wraeros)

!     this subroutine calculates germ formation rates on insoluble subsrate for
!     given T, S, etc.
!     Khvorostyanov and Curry, 2005 JAS
!     All units are mks (JMC) 18 Apr 2007

!       purpose
!       input   temp            [K]
!               sw              saturation ratio
!               draeros         [m] dry aerosol radius
!               betafr          aerosol parameter describing composition
!               qvaeros         aerosol solubility in terms of volume fraction (unitless)
!       output  wraeros         [m]
!
!       comment
!               dilute solution approx.
!               works for RHW <=100.0
!       Note
!               output wraeros >= draeros
!       source:
!               Hugh Morrison, Jennifer Comstock
!       reference:
!               Khvorostyanov and Curry (1999) JGR
!               http://mathworld.wolfram.com
!===============================================================

IMPLICIT NONE

!      INCLUDE "aerosol_prop_mks.inc"

real      temp            ! temperature [K}
real      sw              ! saturation ratio (water)
real      draeros         ! dry aerosol radius [m]
real      Tc,T0
real      H
real      sigvl,rhow
real      b, A, BB
real      wr100
real      a0,a1,a2,s1,s2
real      Q,Q3,R,R2,D,D12
real      theta
real      wraeros
real      Mw,Rv
!real      rhos2

!print*,'Req_haze...'

!rhos2 = 1770.0  !density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
!Ms2 = 0.132     !molecular weight of dry aerosol (g/mol) (NH4)2SO4 KC1999
Mw = 1.8e-2      !molecular weight of water (kg/mol)
Rv = 4.615e2     !specific gas constant for water vapor (J/kg/K)
T0 = 273.15
Tc= temp-T0      ! K to C
H = sw+1.0       ! saturation ratio
rhow=  1000.0    ! density of water (kg/m3)

sigvl= 0.0761 - 1.55e-4 * Tc

!aerosol parameter describing composition,  Eq (11), KC1999
b =    2.0 * qvaero*(rhos2/rhow)*(Mw/Ms2)
!Kelvin parameter, KC1999
BB =   2.0 * sigvl/(Rv*temp*rhow)
A =    b * (draeros)**(2.0+2.0*betafr)
wr100= (A/BB)**0.5   !wet radius at RH=100%

!print*,'BB kel=',BB
!if (h > 0.0) print*,'check para', qvaero, rhos2, Ms2, draeros, betafr,h

if (H .lt. 1.0) then

!     cubic formula, http://mathworld.wolfram.com

   a0     = A/(H-1.)
   a1     = 0.
   a2     = -BB/(H-1.)
!   print*,a0,a1,a2
!     cubic root solution

   Q= (3.*a1   -(a2*a2))/9.0
   R= (9.*a1*a2-27.*a0-2.*a2*a2*a2)/54.0
   Q3     = Q*Q*Q
   R2     = R*R
   D= Q3+R2
!   print*,'Q&R',Q,R
! print*,'Dvalue', D, Q, R, a0, a2, A, BB

   if (D .gt. 0.0) then
      D12 = D**0.5

      if ( (R + D12).gt.0.0) then
         s1 = (R + D12)**(1./3.)
      else
         s1 = -(abs(R + D12)**(1./3.))
      endif

      if ( (R - D12).gt.0.0) then
         s2 = (R - D12)**(1./3.)
      else
         s2 = -(abs(R - D12))**(1./3.)
      endif
!      print*,'s1&s2=',s1,s2

      wraeros = (s1+s2)-a2/3.
!      print*,'D gt 0:',wraeros
   else if (D .lt. 0.0) then
!     3 real solutions
!     choose the first solution (most realistic)
      theta = acos(R/sqrt(-Q3))

      wraeros =   2 * ((-Q)**0.5) * cos(  theta / 3. ) - a2/3.0
!      print*,'D lt 0:',wraeros
   endif
else
!     for water saturated or greater

   wraeros = wr100
!     make sure aerosol size does not exceed value when RH = 100%

   if (wraeros .gt. wr100) then
      wraeros = wr100
   endif

!     make sure wraeros is not smaller than draeros
   wraeros = AMAX1(wraeros,draeros)
endif

RETURN
END SUBROUTINE req_haze

!##############################################################################
Subroutine NUCLEATION (SUP1,TT,FCCNR,DROPCONCN,NDROPMAX,RCCN,DROPRADII,NKR)

! DROPCONCN(KR), 1/cm^3 - drop bin concentrations, KR=1,...,NKR

! determination of new size spectra due to drop nucleation

IMPLICIT NONE
INTEGER NDROPMAX,IDROP,INEXT,ISMALL,KR,NCRITI
INTEGER IMIN,IMAX,NKR,I
REAL &
 SUP1,TT,RACTMAX,RCRITI,BKOE, &
 AKOE,DEG01
REAL CCNCONC(NKR)
REAL CCNCONC_BFNUCL,CCNCONC_AFNUCL, DEl_CCNCONC
REAL RCCN(NKR),DROPRADII(NKR),FCCNR(NKR)
REAL RACT(NKR),DROPCONCN(NKR)

  DEG01=1./3.


! calculation initial value of NDROPMAX - maximal number of drop bin
! which is activated

! initial value of NDROPMAX

  NDROPMAX=0

  DO KR=1,NKR
! initialization of bin radii of activated drops
     RACT(KR)=0.
! initialization of aerosol(CCN) bin concentrations
     CCNCONC(KR)=0.
! initialization of drop bin concentrations
     DROPCONCN(KR)=0.
  ENDDO

! CCNCONC_BFNUCL - concentration of aerosol particles before
!nucleation

  CCNCONC_BFNUCL=0.
  DO I=1,NKR
     CCNCONC_BFNUCL=CCNCONC_BFNUCL+FCCNR(I)
  ENDDO

  CCNCONC_BFNUCL=CCNCONC_BFNUCL*COL

  IF(CCNCONC_BFNUCL.EQ.0.) THEN
     RETURN    
  ELSE
     CALL BOUNDARY (IMIN,IMAX,FCCNR,NKR)
     CALL CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)
!     print*, 'rcriti',RCRITI,imax,RCCN(IMAX)
 
     IF(RCRITI.GE.RCCN(IMAX))  RETURN
  END IF

! calculation of CCNCONC(I) - aerosol(CCN) bin concentrations;
!     I=IMIN,...,IMAX
! determination of NCRITI - number bin in which is located RCRITI
  IF (IMIN.EQ.1)THEN
   CALL CCNIMIN (IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
   CALL CCNLOOP (IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
  ELSE
   CALL CCNLOOP (IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC, &
      FCCNR,NKR)
  END IF

! calculation CCNCONC_AFNUCL - ccn concentration after nucleation

 CCNCONC_AFNUCL=0.

 DO I=IMIN,IMAX
    CCNCONC_AFNUCL=CCNCONC_AFNUCL+FCCNR(I)
 ENDDO

 CCNCONC_AFNUCL=CCNCONC_AFNUCL*COL

 DEL_CCNCONC=CCNCONC_BFNUCL-CCNCONC_AFNUCL

 CALL ACTIVATE (IMIN,IMAX,AKOE,BKOE,RCCN,RACTMAX,NKR)
 CALL DROPMAX (DROPRADII,RACTMAX,NDROPMAX,NKR)
! put nucleated droplets into the drop bin according to radius
! change in drop concentration due to activation DROPCONCN(IDROP)
 ISMALL=NCRITI

 INEXT=ISMALL

 DO IDROP=1,NDROPMAX
    DROPCONCN(IDROP)=0.
    DO I=ISMALL,IMAX
       IF(RACT(I).LE.DROPRADII(IDROP)) THEN
          DROPCONCN(IDROP)=DROPCONCN(IDROP)+CCNCONC(I)
          INEXT=I+1
       ENDIF
    ENDDO
    ISMALL=INEXT
 ENDDO

RETURN
END SUBROUTINE NUCLEATION

!##############################################################################
Subroutine BOUNDARY (IMIN,IMAX,FCCNR,NKR)

! IMIN - left CCN spectrum boundary
IMPLICIT NONE
INTEGER I,IMIN,IMAX,NKR
REAL FCCNR(NKR)

IMIN=0

DO I=1,NKR
   IF(FCCNR(I).NE.0.) THEN
     IMIN=I
     EXIT
   ENDIF
ENDDO

! IMAX - right CCN spectrum boundary
IMAX=0

DO I=NKR,1,-1
   IF(FCCNR(I).NE.0.) THEN
     IMAX=I
     EXIT
   ENDIF
ENDDO

RETURN
END SUBROUTINE BOUNDARY

!##############################################################################
Subroutine CRITICAL (AKOE,BKOE,TT,RCRITI,SUP1,DEG01)

! AKOE & BKOE - constants in Koehler equation
IMPLICIT NONE
REAL AKOE,BKOE,TT,RCRITI,SUP1,DEG01

  AKOE=3.3E-05/TT
  BKOE=ions*4.3/mwaero
  BKOE=BKOE*(4./3.)*3.141593*RO_SOLUTE                  

! RCRITI, cm - critical radius of "dry" aerosol

  RCRITI=(AKOE/3.)*(4./BKOE/SUP1/SUP1)**DEG01

RETURN
END SUBROUTINE CRITICAL

!##############################################################################
Subroutine CCNIMIN (IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,FCCNR,NKR)

IMPLICIT NONE
INTEGER IMIN,II,IMAX,NCRITI,NKR
REAL RCRITI
REAL RCCN(NKR),FCCNR(NKR),CCNCONC(NKR)
REAL RCCN_MIN
REAL DLN1,DLN2
! rccn_min - minimum aerosol(ccn) radius
  RCCN_MIN=RCCN(1)/10000.
! calculation of ccnconc(ii)=fccnr(ii)*col - aerosol(ccn) bin
!                                            concentrations,
!                                            ii=imin,...,imax
! determination of ncriti   - number bin in which is located rcriti
! calculation of ccnconc(ncriti)=fccnr(ncriti)*dln1/(dln1+dln2),
! where,    
! dln1=Ln(rcriti)-Ln(rccn_min)
! dln2=Ln(rccn(1)-Ln(rcriti)

IF(RCRITI.LE.RCCN_MIN) THEN
   NCRITI=1
   DO II=NCRITI+1,IMAX
      CCNCONC(II)=COL*FCCNR(II)     
      FCCNR(II)=0.
   ENDDO
ELSEIF(RCRITI.GT.RCCN_MIN.AND.RCRITI.LT.RCCN(IMIN)) THEN
   NCRITI=1
   DO II=NCRITI+1,IMAX
      CCNCONC(II)=COL*FCCNR(II)
      FCCNR(II)=0.
   ENDDO
   DLN1=ALOG(RCRITI)-ALOG(RCCN_MIN)
   DLN2=ALOG(RCCN(1))-ALOG(RCRITI)
   CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
   FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/(DLN1+DLN2)
ENDIF

RETURN
END SUBROUTINE CCNIMIN

!##############################################################################
Subroutine CCNLOOP (IMIN,IMAX,RCRITI,NCRITI,RCCN,CCNCONC,FCCNR,NKR)

IMPLICIT NONE
INTEGER I,IMIN,IMAX,NKR,II,NCRITI
REAL RCRITI,RCCN(NKR),CCNCONC(NKR),FCCNR(NKR)
REAL DLN1,DLN2

  IF(IMIN.GT.1) THEN
    IF(RCRITI.LE.RCCN(IMIN-1)) THEN
      NCRITI=IMIN
      DO II=NCRITI,IMAX
         CCNCONC(II)=COL*FCCNR(II)
         FCCNR(II)=0.
      ENDDO
      RETURN
    ENDIF
    IF(RCRITI.LT.RCCN(IMIN).AND.RCRITI.GT.RCCN(IMIN-1)) THEN
       NCRITI=IMIN
      
       DO II=NCRITI+1,IMAX
          CCNCONC(II)=COL*FCCNR(II)
          FCCNR(II)=0.
       ENDDO
       DLN1=ALOG(RCRITI)-ALOG(RCCN(IMIN-1))
       DLN2=COL-DLN1
       CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
       FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
       RETURN
    ENDIF
  ENDIF

  DO I=IMIN,IMAX-1
     IF(RCRITI.EQ.RCCN(I)) THEN
       NCRITI=I+1
       DO II=I+1,IMAX
          CCNCONC(II)=COL*FCCNR(II)
          FCCNR(II)=0.
       ENDDO
       RETURN
     ENDIF
     IF(RCRITI.GT.RCCN(I).AND.RCRITI.LT.RCCN(I+1)) THEN
       NCRITI=I+1
       IF(I.NE.IMAX-1) THEN
         DO II=NCRITI+1,IMAX
            CCNCONC(II)=COL*FCCNR(II)
            FCCNR(II)=0.
         ENDDO
       ENDIF
       DLN1=ALOG(RCRITI)-ALOG(RCCN(I))
       DLN2=COL-DLN1
       CCNCONC(NCRITI)=DLN2*FCCNR(NCRITI)
       FCCNR(NCRITI)=FCCNR(NCRITI)*DLN1/COL
       RETURN
     END IF
  ENDDO

RETURN
END SUBROUTINE CCNLOOP

!##############################################################################
Subroutine ACTIVATE (IMIN,IMAX,AKOE,BKOE,RCCN,RACTMAX,NKR)

IMPLICIT NONE

INTEGER IMIN,IMAX,NKR
INTEGER I,I0,I1
REAL RCCN(NKR)
REAL  R03,SUPCRITI,RACT(NKR),XKOE
REAL AKOE,BKOE,AKOE23,RACTMAX

! Spectrum of activated drops (begin) 
  DO I=IMIN,IMAX
  ! critical water supersaturations appropriating CCN radii
     XKOE=(4./27.)*(AKOE**3./BKOE)
     AKOE23=AKOE*2./3.
     R03=RCCN(I)**3.
     SUPCRITI=SQRT(XKOE/R03)

! RACT(I) - radii of activated drops, I=IMIN,...,IMAX

     IF(RCCN(I).LE.(0.3E-5)) RACT(I)=AKOE23/SUPCRITI
     IF(RCCN(I).GT.(0.3E-5)) RACT(I)=5.*RCCN(I)
  ENDDO

  I0=IMIN

   DO I=IMIN,IMAX-1
      IF(RACT(I+1).LT.RACT(I)) THEN
         I0=I+1
         EXIT
      ENDIF
   ENDDO

   I1=I0-1

   IF(I0.EQ.IMIN) GOTO 47

   IF(I0.EQ.IMAX) THEN
      RACT(IMAX)=RACT(IMAX-1)
      GOTO 47
   ENDIF

   IF(RACT(IMAX).LE.RACT(I0-1)) THEN
      DO I=I0,IMAX
         RACT(I)=RACT(I0-1)
      ENDDO
      GOTO 47
   ENDIF

   DO I=I0+1,IMAX
      IF(RACT(I).GE.RACT(I0-1)) THEN
         I1=I
         EXIT
      ENDIF
   ENDDO

! line interpolation RACT(I) for I=I0,...,I1

   DO I=I0,I1
      RACT(I)=RACT(I0-1)+(I-I0+1)*(RACT(I1)-RACT(I0-1)) &
                           /(I1-I0+1)
   ENDDO

 47    CONTINUE

 RACTMAX=0.

 DO I=IMIN,IMAX
    RACTMAX=AMAX1(RACTMAX,RACT(I))
 ENDDO

RETURN
END SUBROUTINE ACTIVATE

!##############################################################################
Subroutine DROPMAX (DROPRADII,RACTMAX,NDROPMAX,NKR)

IMPLICIT NONE
INTEGER IDROP,NKR,NDROPMAX
REAL RACTMAX,DROPRADII(NKR)
! calculation of NDROPMAX - maximal number of drop bin which
! is activated

  NDROPMAX=1

  DO IDROP=1,NKR
     IF(RACTMAX.LE.DROPRADII(IDROP)) THEN
       NDROPMAX=IDROP
       RETURN
     ENDIF
  ENDDO

RETURN
END SUBROUTINE DROPMAX

!##############################################################################
Subroutine ONECOND1 (TT,QQ,PP,ROR &
      ,VR1,PSINGLE &
      ,DEL1N,DEL2N,DIV1,DIV2 &
      ,FF1,PSI1,R1,RLEC,RO1BL &
      ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
      ,COL,DTCOND,ICEMAX,NKR &
      ,vapcldt,vapraint)

IMPLICIT NONE

INTEGER NKR,ICEMAX
REAL    COL,VR1(NKR),PSINGLE &
        ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
        ,DTCOND

INTEGER KR,ITIME,KCOND,NR,NRM,KLIMIT,KM  
REAL AL1,AL2,D,GAM,POD, &
     RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY,DT0LREF,DTLREF, &
     DT,DTT, &
     TPC1, TPC2, TPC3, TPC4, TPC5, &
     EPSDEL, EPSDEL2,DT0L, DT0I,&
     ROR, CWHUCM,B6,B8L, &
     DEL1,DEL2,DEL1S,DEL2S, &
     TIMENEW,TIMEREV,SFN11,SFN12, &
     SFNL,SFNI,B5L,B5I,B7L,B7I,DOPL,DOPI,RW,RI,QW,PW, &
     PI,QI,DTNEWL,DTNEWL1,D1N,D2N, &
     DT0L0,DT0I0

INTEGER K

REAL  FF1_OLD(NKR),SUPINTW(NKR)
DOUBLE PRECISION DSUPINTW(NKR),DD1N,DB11_MY
DOUBLE PRECISION COL3,DAL1,RORI,TPN,TPS,QPN,QPS,TOLD,QOLD &
                 ,FI1_K,R1_K,FI1R1 &
                 ,RMASSLAA,RMASSLBB &
                 ,ES1N,ES2N,EW1N,ARGEXP &
                 ,TT,QQ,PP &
                 ,DEL1N,DEL2N,DIV1,DIV2 &
                 ,OPER2,OPER3,AR1,AR2

DOUBLE PRECISION DELMASSL1

! DROPLETS 
                                                                       
REAL R1(NKR) &
     ,RLEC(NKR),RO1BL(NKR) &
     ,FI1(NKR),FF1(NKR),PSI1(NKR) &
     ,B11_MY(NKR),B12_MY(NKR)

REAL DTIMEO(NKR),DTIMEL(NKR) &
     ,TIMESTEPD(NKR)

! CCN regeneration
real tot_before, tot_after

!condensation budgets
real vapcldt, vapraint,vapc,vapr

! NEW ALGORITHM (NO TYPE OF ICE)

OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)

DATA AL1 /2500./, AL2 /2834./, D /0.211/ &
     ,GAM /1.E-4/, POD /10./ 
           
DATA RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY &
      /461.5,0.24E-1,0.211E-4,2.5E6,2.834E6/

DATA TPC1, TPC2, TPC3, TPC4, TPC5 &
     /-4.0,-8.1,-12.7,-17.8,-22.4/ 

DATA EPSDEL, EPSDEL2 /0.1E-03,0.1E-03/  
    
DATA DT0L, DT0I /1.E20,1.E20/

vapc=0.
vapr=0.

! CONTROL OF DROP SPECTRUM IN subroutine ONECOND

COL3=3.0*COL
KCOND=0
DT0LREF=0.2
DTLREF=0.4

NR=NKR
DT=DTCOND
DTT=DTCOND

CWHUCM=0.
B6=CWHUCM*GAM
B8L=1./ROR
RORI=1./ROR

!!! add ccn regeneration (calculate the total droplets before cond-evap)
tot_before = 0.0
do k = 1, nkr
   tot_before = tot_before+psi1(k)*3.0*r1(k)*col
enddo
!!!

DO KR=1,NKR
   FF1_OLD(KR)=FF1(KR)
   SUPINTW(KR)=0.
   DSUPINTW(KR)=0.
ENDDO
TPN=TT
QPN=QQ
DO KR=1,NKR
   FI1(KR)=FF1(KR)
ENDDO
TIMENEW=0.
ITIME=0
TOLD=TPN
QOLD=QPN
!----------------Begin small timestep loop----------------------------
56   ITIME=ITIME+1
  TIMEREV=DT-TIMENEW
  DEL1=DEL1N
  DEL2=DEL2N
  DEL1S=DEL1N
  DEL2S=DEL2N
  TPS=TPN
  QPS=QPN
  CALL JERRATE (R1,TPS,PP,ROR,VR1,PSINGLE &
             ,RLEC,RO1BL,B11_MY,B12_MY,1,1,ICEMAX,NKR)

  CALL JERTIMESC (FI1,R1,SFN11,SFN12 &
               ,B11_MY,B12_MY,B8L,1,COL,NKR)        

  SFNL=SFN11+SFN12
  SFNI=0.       

  B5L=BB1_MY/TPS/TPS
  B5I=BB2_MY/TPS/TPS
  B7L=B5L*B6                                                     
  B7I=B5I*B6
  DOPL=1.+DEL1S                                                     
  DOPI=1.+DEL2S                                                     
  RW=(OPER2(QPS)+B5L*AL1)*DOPL*SFNL                                                 
  RI=(OPER2(QPS)+B5L*AL2)*DOPL*SFNI
  QW=B7L*DOPL
  PW=(OPER2(QPS)+B5I*AL1)*DOPI*SFNL
  PI=(OPER2(QPS)+B5I*AL2)*DOPI*SFNI
  QI=B7I*DOPI

  KCOND=10

  IF(DEL1.GT.0) KCOND=11

  IF(KCOND.EQ.11) THEN !Condensation
     DTNEWL=ABS(R1(ITIME)/(B11_MY(ITIME)*DEL1N &
                      -B12_MY(ITIME)))
     IF(DTNEWL.GT.DT) DTNEWL=DT
     TIMESTEPD(ITIME)=DTNEWL

     IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1))  & 
       DTNEWL=DT-TIMENEW
     IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW

     TIMESTEPD(ITIME)=DTNEWL
     TIMENEW=TIMENEW+DTNEWL
     DTT=DTNEWL

! SOLVING FOR SUPERSATURATION 
     CALL JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                   ,RW,PW,RI,PI,QW,QI &
                   ,DTT,D1N,D2N,DT0L,DT0I)

! DROPLET DISTRIBUTION function 

     CALL JERDFUN (R1,B11_MY,B12_MY &
                 ,FI1,PSI1,D1N &
                 ,1,6,COL,NKR,vapc,vapr,rori)

  ELSE ! EVAPORATION - ONLY WATER 

     IF (DEL1N.EQ.0)THEN
        DTIMEO(1)=DT
        DO KR=2,NKR
           DTIMEO(KR)=DT
        ENDDO
     ELSE
        DTIMEO(1)=-R1(1)/(B11_MY(1)*DEL1N-B12_MY(1))
        DO KR=2,NKR
           KM=KR-1
           DTIMEO(KR)=(R1(KM)-R1(KR))/(B11_MY(KR)*DEL1N &
                                      -B12_MY(KR))
        ENDDO
     END IF

     KLIMIT=1

     DO KR=1,NKR
        IF(DTIMEO(KR).GT.TIMEREV) GOTO 55
        KLIMIT=KR
     ENDDO

55   KLIMIT=KLIMIT-1

     IF(KLIMIT.LT.1) KLIMIT=1

     DTNEWL1=AMIN1(DTIMEO(KLIMIT),TIMEREV)
     IF(DTNEWL1.LT.DTLREF) DTNEWL1=AMIN1(DTLREF,TIMEREV)
     DTNEWL=DTNEWL1
     TIMESTEPD(ITIME)=DTNEWL

     IF(DTNEWL.GT.DT) DTNEWL=DT
     IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1))  &
         DTNEWL=DT-TIMENEW
     IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW

     TIMESTEPD(ITIME)=DTNEWL
     TIMENEW=TIMENEW+DTNEWL
     DTT=DTNEWL

! SOLVING FOR SUPERSATURATION 
     CALL JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                    ,RW,PW,RI,PI,QW,QI &
                    ,DTT,D1N,D2N,DT0L0,DT0I0)
       
     CALL JERDFUN (R1,B11_MY,B12_MY &
                  ,FI1,PSI1,D1N &
                  ,1,6,COL,NKR,vapc,vapr,rori)

  ENDIF

  RMASSLBB=0.
  RMASSLAA=0.

! BEFORE JERNEWF (ONLY WATER) 
   DO K=1,NKR
      FI1_K=FI1(K)
      R1_K=R1(K)
      FI1R1=FI1_K*R1_K*R1_K
      RMASSLBB=RMASSLBB+FI1R1
   ENDDO
   RMASSLBB=RMASSLBB*COL3*RORI
   IF(RMASSLBB.LE.0.) RMASSLBB=0.
   DO K=1,NKR
      FI1_K=PSI1(K)
      R1_K=R1(K)
      FI1R1=FI1_K*R1_K*R1_K
      RMASSLAA=RMASSLAA+FI1R1
   ENDDO
   RMASSLAA=RMASSLAA*COL3*RORI
   IF(RMASSLAA.LE.0.) RMASSLAA=0.
   DELMASSL1=RMASSLAA-RMASSLBB
   QPN=QPS-DELMASSL1
   DAL1=AL1
   TPN=TPS+DAL1*DELMASSL1
! SUPERSATURATION (ONLY WATER)
   ARGEXP=-BB1_MY/TPN
   ES1N=AA1_MY*DEXP(ARGEXP)
   ARGEXP=-BB2_MY/TPN
   ES2N=AA2_MY*DEXP(ARGEXP)
   EW1N=OPER3(QPN,PP)
   IF(ES1N.EQ.0)THEN
    DEL1N=0.5
    DIV1=1.5
   ELSE
    DIV1=EW1N/ES1N
    DEL1N=EW1N/ES1N-1.
   END IF
   IF(ES2N.EQ.0)THEN
    DEL2N=0.5
    DIV2=1.5
   ELSE
    DEL2N=EW1N/ES2N-1.
    DIV2=EW1N/ES2N
   END IF
   DO KR=1,NKR
     SUPINTW(KR)=SUPINTW(KR)+B11_MY(KR)*D1N
     DD1N=D1N
     DB11_MY=B11_MY(KR)
     DSUPINTW(KR)=DSUPINTW(KR)+DB11_MY*DD1N
   ENDDO
! REPEAT TIME STEP (ONLY WATER: CONDENSATION OR EVAPORATION) 
   IF(TIMENEW.LT.DT) GOTO 56
!------------------------------------------------------------
CALL JERDFUN_NEW (R1,DSUPINTW &
                 ,FF1_OLD,PSI1,D1N &
                 ,1,1,COL,NKR,vapcldt,vapraint,rori)

RMASSLAA=0.0
RMASSLBB=0.0
! BEFORE JERNEWF
DO K=1,NKR
   FI1_K=FF1_OLD(K)
   R1_K=R1(K)
   FI1R1=FI1_K*R1_K*R1_K
   RMASSLBB=RMASSLBB+FI1R1
ENDDO
RMASSLBB=RMASSLBB*COL3*RORI

IF(RMASSLBB.LT.0.0) RMASSLBB=0.0
! AFTER  JERNEWF
!add CCN regeneration. Calculate the total droplets after condensation -J. Fan
tot_after = 0.
DO K=1,NKR
   FI1_K=PSI1(K)
   R1_K=R1(K)
   FI1R1=FI1_K*R1_K*R1_K
   RMASSLAA=RMASSLAA+FI1R1
   tot_after = tot_after+psi1(k)*3.0*r1(k)*col
ENDDO

! CCN regeneration from evaporation
ccnreg = max((tot_before-tot_after),0.0)

RMASSLAA=RMASSLAA*COL3*RORI
IF(RMASSLAA.LT.0.0) RMASSLAA=0.0
DELMASSL1=RMASSLAA-RMASSLBB

QPN=QOLD-DELMASSL1
DAL1 = AL1
TPN=TOLD+DAL1*DELMASSL1
! SUPERSATURATION
ARGEXP=-BB1_MY/TPN
ES1N=AA1_MY*DEXP(ARGEXP)
ARGEXP=-BB2_MY/TPN
ES2N=AA2_MY*DEXP(ARGEXP)
EW1N=OPER3(QPN,PP)
DIV1=EW1N/ES1N
DEL1N=EW1N/ES1N-1.
DEL2N=EW1N/ES2N-1.
DIV2=EW1N/ES2N
TT=TPN
QQ=QPN
DO KR=1,NKR
   FF1(KR)=PSI1(KR)
ENDDO

RETURN
END SUBROUTINE ONECOND1

!##############################################################################
Subroutine JERDFUN (R2,B21_MY,B22_MY,FI2,PSI2,DEL2N,IND,ITYPE,COL,NKR, &
                    vapcldt,vapraint,rori)

IMPLICIT NONE
REAL COL,DEL2N
INTEGER IND,ITYPE,KR,ICE,ITYP,NRM,NR,NKR
REAL R2(NKR,IND),R2N(NKR,IND) &
     ,FI2(NKR,IND),PSI2(NKR,IND) &
     ,B21_MY(NKR,IND),B22_MY(NKR,IND) &
     ,DEL_R2M(NKR,IND)
DOUBLE PRECISION R2R(NKR),R2NR(NKR),FI2R(NKR),PSI2R(NKR)
DOUBLE PRECISION DR2(NKR,IND),DR2N(NKR,IND),DDEL2N, &
          DB21_MY(NKR,IND)
DOUBLE PRECISION CHECK
REAL vapcldt, vapraint
DOUBLE PRECISION rori

CHECK=0.D0
DO KR=1,NKR
   CHECK=B21_MY(1,1)*B21_MY(KR,1)
   IF (CHECK.LT.0)print*,'CHECK < 0',KR,B21_MY(1,1),B21_MY(KR,1)
   IF (CHECK.LT.0)STOP
END DO

IF(IND.NE.1) THEN
  ITYP=ITYPE
  ELSE
  ITYP=1
ENDIF

DDEL2N=DEL2N
DO KR=1,NKR
   PSI2R(KR)=FI2(KR,ITYP)
   FI2R(KR)=FI2(KR,ITYP)
   DR2(KR,ITYP)=R2(KR,ITYP)
   DB21_MY(KR,ITYP)=B21_MY(KR,ITYP)
ENDDO

NR=NKR
NRM=NKR-1

! NEW DISTRIBUTION function
DO ICE=1,IND
  IF(ITYP.EQ.ICE) THEN
     DO KR=1,NKR
        DR2N(KR,ICE)=DR2(KR,ICE)+DDEL2N*DB21_MY(KR,ICE)
        R2N(KR,ICE)=DR2N(KR,ICE)
     ENDDO
  ENDIF
ENDDO
    
IF (ITYPE==6) THEN
  DO KR=1,KRDROP-1
   vapcldt=vapcldt+3.*FI2R(KR)*R2(KR,1)*(max(0.,R2N(KR,1))-R2(KR,1))*COL*RORI
  ENDDO
  DO KR=KRDROP,NKR
   vapraint=vapraint+3.*FI2R(KR)*R2(KR,1)*(max(0.,R2N(KR,1))-R2(KR,1))*COL*RORI
  ENDDO
ENDIF

DO ICE=1,IND
   IF(ITYP.EQ.ICE) THEN
      DO KR=1,NKR
         R2R(KR)=DR2(KR,ICE)
         R2NR(KR)=DR2N(KR,ICE)   
      ENDDO

      CALL JERNEWF (NR,NRM,R2R,FI2R,PSI2R,R2NR,COL,NKR)
      DO KR=1,NKR
         PSI2(KR,ICE)=PSI2R(KR)
      ENDDO
   ENDIF
ENDDO

RETURN
END SUBROUTINE JERDFUN

!##############################################################################
Subroutine JERDFUN_NEW (R2,B21_MY &
                         ,FI2,PSI2,DEL2N &
                         ,IND,ITYPE,COL,NKR,vapcldt,vapraint,rori)

IMPLICIT NONE
INTEGER IND,ITYPE,KR,ICE,NRM,NR,NKR
REAL R2(NKR,IND),FI2(NKR,IND),PSI2(NKR,IND),COL,DEL2N
DOUBLE PRECISION  B21_MY(NKR,IND)
DOUBLE PRECISION R2R(NKR),R2NR(NKR),FI2R(NKR),PSI2R(NKR)
REAL vapcldt,vapraint
DOUBLE PRECISION rori

DO KR=1,NKR
   FI2R(KR)=FI2(KR,ITYPE)
ENDDO

NR=NKR
NRM=NKR-1
ICE = ITYPE

DO KR=1,NKR
   R2R(KR)=R2(KR,ICE)
   R2NR(KR)=R2(KR,ICE)+B21_MY(KR,ICE)
ENDDO

  DO KR=1,KRDROP-1
   vapcldt=vapcldt+3.*FI2R(KR)*R2R(KR)*(max(0.,R2NR(KR))-R2R(KR))*COL*RORI
  ENDDO
  DO KR=KRDROP,NKR
   vapraint=vapraint+3.*FI2R(KR)*R2R(KR)*(max(0.,R2NR(KR))-R2R(KR))*COL*RORI
  ENDDO

CALL JERNEWF (NR,NRM,R2R,FI2R,PSI2R,R2NR,COL,NKR)
DO KR=1,NKR
   PSI2(KR,ICE)=PSI2R(KR)
ENDDO

RETURN
END SUBROUTINE JERDFUN_NEW

!##############################################################################
Subroutine JERNEWF (NRX,NRM,RR,FI_OLD,PSI,RN,COL,NKR)
 
IMPLICIT NONE

INTEGER  I,K,KM,NRXP,IM,IP,IFIN,IIN,NKR
REAL COL 
DOUBLE PRECISION RNTMP,RRTMP,RRP,RRM,RNTMP2,RRTMP2,RRP2,RRM2, &
        GN1,GN1P,GN2,GN3,GMAT2

DOUBLE PRECISION DRP,FNEW,FIK,PSINEW,DRM,GMAT,R1,R2,R3,RRI,RNK

INTEGER NRX,NRM

DOUBLE PRECISION RR(NRX),FI(NRX),PSI(NRX),RN(NRX) &
       ,RRS(NKR+1),RNS(NKR+1),PSIN(NKR+1),FIN(NKR+1)

DOUBLE PRECISION FI_OLD(NRX)
DOUBLE PRECISION PSI_IM,PSI_I,PSI_IP
 
NRXP=NRX+1

DO K=1,NRX
   FI(K)=FI_OLD(K)
ENDDO
 
DO K=1,NRX
   PSI(K)=0.0D0
ENDDO

IF(RN(NRX).NE.RR(NRX)) THEN
! Kovetz-Olund method     (start)
   IF(RN(1).GT.RR(1)) THEN !CONDENSATION
      RNS(NRXP)=1024.0D0*RR(NRX)
      RRS(NRXP)=1024.0D0*RR(NRX)
      PSIN(NRXP)=0.0D0
      FIN(NRXP)=0.0D0

      DO K=1,NRX
         RNS(K)=RN(K)
         RRS(K)=RR(K)
         PSIN(K)=0.0D0
! FIN(K) - initial(before condensation) concentration of hydrometeors
         FIN(K)=3.0D0*FI(K)*RR(K)*COL
      ENDDO

! NUMBER OF NEW RADII POSITION IN REGULAR GRID 
! RNK - new first bin mass(after condensation)
      RNK=RNS(1)

      DO I=1,NRX
         RRI=RRS(I)
         IF(RRI.GT.RNK) EXIT
      ENDDO
      IIN=I-1

      IFIN=NRX

      DO I=IIN,IFIN
         IP=I+1
         IM=MAX(1,I-1)

         R1=RRS(IM)
         R2=RRS(I)
         R3=RRS(IP)

         DRM=R2-R1
         DRP=R3-R2

         FNEW=0.0D0

         DO K=1,I
            FIK=FIN(K)

            IF(FIK.NE.0.0D0) THEN
               KM=K-1

! RNK - new bin mass(after condensation)
               RNK=RNS(K)

               IF(RNK.NE.R2) THEN
                  GMAT=0.0D0
                  IF(RNK.GT.R1.AND.RNK.LT.R3) THEN
                     IF(RNK.LT.R2) THEN
                        GMAT=(RNK-R1)/DRM
                     ELSE
                        GMAT=(R3-RNK)/DRP
                      ENDIF
                  ENDIF
                ELSE
                    GMAT=1.0D0
                ENDIF

                FNEW=FNEW+FIK*GMAT
            ENDIF
         ENDDO
! PSIN(I)) - new concentration of hydrometeors after condensation
         PSIN(I)=FNEW
      ENDDO
! PSI(K) - new size distribution function of hydrometeors after 
!          condensation, K=1,...,NRX=NKR
      DO K=1,NRX
         PSI(K)=PSIN(K)/3./RR(K)/COL
      ENDDO

   ELSE !EVAPORATION

      RNS(1)=0.0D0
      RRS(1)=0.0D0
      FIN(1)=0.0D0
      PSIN(1)=0.0D0

! FIN(K) - initial(before evaporation) concentration of hydrometeors
      DO K=2,NRXP
         KM=K-1
         RNS(K)=RN(KM)
         RRS(K)=RR(KM)
         PSIN(K)=0.0D0
         FIN(K)=3.0D0*FI(KM)*RR(KM)*COL
      ENDDO

      DO I=1,NRXP
         IM=MAX(I-1,1)
         IP=MIN(I+1,NRXP)

         R1=RRS(IP)
         R2=RRS(I)
         R3=RRS(IM)

         DRM=R1-R2
         DRP=R2-R3

         FNEW=0.0D0

         DO K=I,NRXP
            RNK=RNS(K)
            IF(RNK.GE.R1) EXIT
            IF(RNK.GT.R3)THEN
               IF(RNK.GT.R2) THEN
                   FNEW=FNEW+FIN(K)*(R1-RNK)/DRM
               ELSE
                   FNEW=FNEW+FIN(K)*(RNK-R3)/DRP
               ENDIF
            ENDIF
         ENDDO

! PSIN(I) - new concentration of hydrometeors after evaporation

          PSIN(I)=FNEW

      ENDDO
! NEW VALUES OF DISTRIBUTION function (start)
! PSI(K), 1/g/cm^3 - new size distribution function of hydrometeors 
!                    after evaporation, K=1,...,NRX
      DO K=2,NRXP
         KM=K-1
         R1=PSIN(K)*RR(KM)
         PSINEW=PSIN(K)/3.0D0/RR(KM)/COL
         IF(R1.LT.1.0D-20) PSINEW=0.0D0
         PSI(KM)=PSINEW
      ENDDO
   ENDIF !End if condensation or evaporation

   IF(I3POINT.NE.0) THEN

      DO K=1,NKR
         RRS(K)=RR(K)
      ENDDO

      RRS(NKR+1)=RRS(NKR)*1024.0D0

      DO I=1,NKR 
         PSI(I)=PSI(I)*RR(I)
! PSI(I) - concenration hydrometeors after KO divided on COL*3.0D0
! RN(I), g - new masses after condensation or evaporation
         IF(RN(I).LT.0.0D0) THEN 
            RN(I)=1.0D-50
            FI(I)=0.0D0
         ENDIF

      ENDDO
 
      DO K=1,NKR
         IF(FI(K).NE.0.0D0) THEN
            IF(RRS(2).LT.RN(K)) THEN
               I=2

               DO WHILE (.NOT.(RRS(I).LT.RN(K).AND.RRS(I+1).GT.RN(K)) &
                          .AND.I.LT.NKR)
                  I=I+1
               ENDDO
               IF(I.LT.NKR-2) THEN
                  RNTMP=RN(K)
                  RRTMP=RRS(I)
                  RRP=RRS(I+1)
                  RRM=RRS(I-1)
 
                  RNTMP2=RN(K+1)

                  RRTMP2=RRS(I+1)
                  RRP2=RRS(I+2)
                  RRM2=RRS(I)
 
                  GN1=(RRP-RNTMP)*(RRTMP-RNTMP)/(RRP-RRM)/ &
                      (RRTMP-RRM)

                  GN1P=(RRP2-RNTMP2)*(RRTMP2-RNTMP2)/ &
                       (RRP2-RRM2)/(RRTMP2-RRM2)

                  GN2=(RRP-RNTMP)*(RNTMP-RRM)/(RRP-RRTMP)/ &
                      (RRTMP-RRM)
 
                  GMAT=(RRP-RNTMP)/(RRP-RRTMP)
                  GN3=(RRTMP-RNTMP)*(RRM-RNTMP)/(RRP-RRM)/ &
                                                (RRP-RRTMP)
                  GMAT2=(RNTMP-RRTMP)/(RRP-RRTMP)

                  PSI_IM=PSI(I-1)+GN1*FI(K)*RR(K)
                  PSI_I=PSI(I)+(GN1P+GN2-GMAT)*FI(K+1)*RR(K+1)
                  PSI_IP=PSI(I+1)+(GN3-GMAT2)*FI(K)*RR(K)
                    
                  IF(PSI_IM.GT.0.0D0) THEN
                     IF(PSI_IP.GT.0.0D0) THEN
                        IF(I.GT.2) THEN
! smoothing criteria
                           IF(PSI_IM.GT.PSI(I-2).AND.PSI_IM.LT.PSI_I &
                              .AND.PSI(I-2).LT.PSI(I).OR.PSI(I-2) &
                                 .GE.PSI(I)) THEN

                              PSI(I-1)=PSI_IM

                              PSI(I)=PSI(I)+FI(K)*RR(K)*(GN2-GMAT)

                              PSI(I+1)=PSI_IP
                           ENDIF 
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
!!! CCN generation from evaporation J.Fan Oct.2007
!               else          ! RRS(2).GE.RN(K)
!               ccnreg = ccnreg+fi(k)*3*col*rr(k)
!               ccnreg = ccnreg+fi(k)*3*col*xccn(k)
!!!
            ENDIF
         ENDIF
      ENDDO
! PSI(K) - new hydrometeor size distribution function
      DO K=1,NKR
         PSI(K)=PSI(K)/RR(K)
      ENDDO
   ENDIF !If using 3 point redistribution
ELSE !If droplet sizes haven't changed - don't need to do anything
   DO K=1,NKR
      PSI(K)=FI(K)
   ENDDO
ENDIF

RETURN 
END SUBROUTINE JERNEWF

!##############################################################################
Subroutine JERRATE (R1S,TP,PP,ROR,VR1,PSINGLE,RIEC,RO1BL &
                   ,B11_MY,B12_MY,ID,IN,ICEMAX,NKR)

IMPLICIT NONE
INTEGER ID,IN,KR,ICE,NRM,ICEMAX,NKR
DOUBLE PRECISION TP,PP
REAL DETL,FACTPL,VENTPL,CONSTL,RVT,D_MY,CONST
REAL VR1(NKR,ID),PSINGLE &
     ,R1S(NKR,ID),B11_MY(NKR,ID),B12_MY(NKR,ID) &
     ,RO1BL(NKR,ID),RIEC(NKR,ID) &
     ,FD1(NKR,ICEMAX),FK1(NKR,ICEMAX),FACTRL(NKR,ICEMAX) &
     ,R11_MY(NKR,ICEMAX) &
     ,R1_MY1(NKR,ICEMAX) &
     ,AL1(2),AL1_MY(2),A1_MY(2),BB1_MY(2),ESAT1(2),CONSTLI(ICEMAX)
DOUBLE PRECISION TZERO
REAL PZERO,CF_MY,D_MYIN,RV_MY,DEG01,DEG03
REAL COEFF_VISCOUS,SHMIDT_NUMBER,A,B
REAL REINOLDS_NUMBER,RESHM,ROR
PARAMETER (TZERO=273.150,PZERO=1.013E6)
DATA AL1/2500.,2833./

DEG01=1./3.     
DEG03=1./3.     
CONST=12.566372 !4*pi
AL1_MY(1)=2.5E10
AL1_MY(2)=2.834E10
A1_MY(1)=2.53E12
A1_MY(2)=3.41E13
BB1_MY(1)=5.42E3
BB1_MY(2)=6.13E3
CF_MY=2.4E3
D_MYIN=0.221
RV_MY=461.5E4
NRM=NKR-1
! rhs for "maxwell" equation
! coefficient of diffusion
D_MY=D_MYIN*(PZERO/PP)*(TP/TZERO)**1.94
! coefficient of viscousity
COEFF_VISCOUS=1.72E-2*SQRT(TP/273.)*393./(TP-120.)/ROR
! Shmidt number
SHMIDT_NUMBER=COEFF_VISCOUS/D_MY
! Constants used for calculation of Reynolds number
A=2.*(3./4./3.141593)**DEG01
B=A/COEFF_VISCOUS
  
RVT=RV_MY*TP
ESAT1(IN)=A1_MY(IN)*EXP(-BB1_MY(IN)/TP)
DO ICE=1,ID
   DO KR=1,NKR
! Reynolds numbers
      REINOLDS_NUMBER=B*VR1(KR,ICE)*SQRT(1.E6/PSINGLE)* &
       (R1S(KR,ICE)/RO1BL(KR,ICE))**DEG03
      RESHM=SQRT(REINOLDS_NUMBER)*SHMIDT_NUMBER**DEG03
      !Adele - see P&K Ch. 13
      IF(REINOLDS_NUMBER.LT.2.5) THEN
         VENTPL=1.+0.108*RESHM*RESHM
      ELSE
         VENTPL=0.78+0.308*RESHM
      ENDIF
      CONSTL=CONST*RIEC(KR,ICE) 
      CONSTLI(ICE)=CONSTL
      FACTPL=1.     
      FACTRL(KR,ICE)=FACTPL     
      FD1(KR,ICE)=RVT/D_MY/ESAT1(IN)/FACTPL 
      FK1(KR,ICE)=(AL1_MY(IN)/RVT-1.)*AL1_MY(IN)/CF_MY/TP
      R1_MY1(KR,ICE)=VENTPL*CONSTL
      R11_MY(KR,ICE)=R1_MY1(KR,ICE)
! growth rate 
      DETL=FK1(KR,ICE)+FD1(KR,ICE)
      B11_MY(KR,ICE)=R11_MY(KR,ICE)/DETL
      B12_MY(KR,ICE)=0.
   ENDDO
ENDDO

RETURN
END SUBROUTINE JERRATE

!##############################################################################
Subroutine JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                     ,RW,PW,RI,PI,QW,QI &
                     ,DT,DEL1INT,DEL2INT,DT0L,DT0I)

IMPLICIT NONE
REAL DEL1,DEL2,RW,PW,RI,PI,QW,QI, &
     DT,DEL1INT,DEL2INT,DT0L,DT0I,DTLIN,DTIIN
REAL DETER,ALFAMX,ALFAPX, &
     R1RES,R2RES,R1,R2,R3,R4,R21,R11,R10,R41,R31,R30
DOUBLE PRECISION DEL1N,DEL2N
DOUBLE PRECISION DEL1N_2P,DEL1INT_2P,DEL2N_2P,DEL2INT_2P 
DOUBLE PRECISION EXPP_2P,EXPM_2P,ARGEXP     
DOUBLE PRECISION RW_DP,PW_DP,PI_DP,RI_DP,X_DP,ALFA_DP
DTLIN=1000.E17
DTIIN=1000.E17
DETER=RW*PI-PW*RI

IF(DETER.EQ.0)  THEN
  IF(RW.EQ.0.AND.RI.EQ.0) THEN ! WITHOUT WATER & ICE
    DEL1N_2P=DEL1
    DEL2N_2P=DEL2
    DEL1INT_2P=DEL1*DT
    DEL2INT_2P=DEL2*DT
  ELSE ! (WATER OR ICE) 
    IF(RW.NE.0) THEN ! ONLY WATER
      ARGEXP=-RW*DT
      DEL1N_2P=DEL1*DEXP(ARGEXP)+QW*(1.-DEXP(ARGEXP))
      DEL1INT_2P=(DEL1-DEL1N_2P)/RW
      DEL2N_2P=DEL2-PW*DEL1INT_2P
      DEL2INT_2P= &
      (DEL2N_2P-PW*DEL1N_2P/RW)*DT+PW*DEL1INT_2P/RW
    ELSE ! ONLY ICE 
      ARGEXP=-PI*DT
      DEL2N_2P=DEL2*DEXP(ARGEXP)+QI*(1.-DEXP(ARGEXP))
      DEL2INT_2P=(DEL2-DEL2N_2P)/PI
      DEL1N_2P=DEL1-RI*DEL2INT_2P
      DEL1INT_2P= (DEL1N_2P-RI*DEL2N_2P/PI)*DT+RI*DEL2INT_2P/PI
     ENDIF
  ENDIF
ELSE ! IN CASE: DETER.NE.0 ! COMPLETE SOLUTION 
  RW_DP=RW
  RI_DP=RI
  PI_DP=PI
  PW_DP=PW
  IF (RW.LE.0)PRINT*,'RW = ',RW
  IF (PW.LE.0)PRINT*,'PW = ',PW
  IF (RI.LE.0)PRINT*,'RI = ',RI
  IF (PI.LE.0)PRINT*,'PI = ',PI,DEL1,DEL2,QW,QI
  IF (RW.LE.0 .OR. PW.LE.0 .OR. RI.LE.0 .OR. PI.LE.0)STOP
  ALFA_DP=SQRT((RW_DP-PI_DP)*(RW_DP-PI_DP)+4.*PW_DP*RI_DP) 
  X_DP=RW_DP+PI_DP
  ALFAPX=.5*(ALFA_DP+X_DP)
  IF (ALFAPX.LE.0)PRINT*,'ALFAPX=',ALFAPX
  IF (ALFAPX.LE.0)STOP
  ALFAMX=.5*(ALFA_DP-X_DP)
 
  ARGEXP=-ALFAPX*DT
  EXPP_2P=DEXP(ARGEXP)
  ARGEXP=ALFAMX*DT
  EXPM_2P=DEXP(ARGEXP)
! DROPLETS 
  R10=RW*DEL1+RI*DEL2
  R11=R10-ALFAPX*DEL1
  R21=R10+ALFAMX*DEL1
  DEL1N_2P=(R21*EXPP_2P-R11*EXPM_2P)/ALFA_DP
  IF(ALFAMX.NE.0) THEN
    R1=-R11/ALFAMX
    R2=R21/ALFAPX
    DEL1INT_2P=(R1*(EXPM_2P-1.)-R2*(EXPP_2P-1.))/ALFA_DP
  ELSE
    DEL1INT_2P = 0.
  ENDIF

  R1RES=0.
  IF(R11.NE.0) R1RES=R21/R11
  IF(R1RES.GT.0) DTLIN=ALOG(R1RES)/ALFA_DP
! ICE 
  R30=PW*DEL1+PI*DEL2
  R31=R30-ALFAPX*DEL2
  R41=R30+ALFAMX*DEL2

  DEL2N_2P=(R41*EXPP_2P-R31*EXPM_2P)/ALFA_DP
  IF(ALFAMX.NE.0.AND.ALFAPX.NE.0) THEN
    R3=-R31/ALFAMX
    R4=R41/ALFAPX
    DEL2INT_2P=(R3*(EXPM_2P-1.)-R4*(EXPP_2P-1.))/ALFA_DP
  ELSE
    DEL2INT_2P=0.
  ENDIF
  R2RES=0.
  IF(R31.NE.0) R2RES=R41/R31
  IF(R2RES.GT.0) DTIIN=ALOG(R2RES)/ALFA_DP
! END OF COMPLETE SOLUTION
ENDIF
DEL1N=DEL1N_2P
DEL2N=DEL2N_2P
 
DEL1INT=DEL1INT_2P
DEL2INT=DEL2INT_2P
DT0L=DTLIN
IF(DT0L.LT.0) DT0L=1.E20
DT0I=DTIIN
IF(DT0I.LT.0) DT0I=1.E20

RETURN
END SUBROUTINE JERSUPSAT

!##############################################################################
Subroutine JERTIMESC (FI1,X1,SFN11,SFN12 &
                     ,B11_MY,B12_MY,CF,ID,COL,NKR)

IMPLICIT NONE
INTEGER NRM,KR,ICE,ID,NKR
REAL B12,B11,FUN,DELM,FK,CF,SFN12S,SFN11S
REAL  COL,X1(NKR,ID),FI1(NKR,ID),B11_MY(NKR,ID),B12_MY(NKR,ID) &
     ,SFN11,SFN12

NRM=NKR-1
DO ICE=1,ID  
   SFN11S=0.
   SFN12S=0.
   SFN11=CF*SFN11S      
   SFN12=CF*SFN12S
   DO KR=1,NRM
! VALUE OF DISTRIBUTION function
      FK=FI1(KR,ICE)
      DELM=X1(KR,ICE)*3.*COL
! INTEGRAL'S EXPRESSION 
      FUN=FK*DELM
! VALUES OF INTEGRALS
      B11=B11_MY(KR,ICE)
      B12=B12_MY(KR,ICE)
      SFN11S=SFN11S+FUN*B11 
      SFN12S=SFN12S+FUN*B12
   ENDDO
     SFN11=CF*SFN11S
     SFN12=CF*SFN12S
ENDDO

RETURN
END SUBROUTINE JERTIMESC

!##############################################################################
Subroutine JERTIMESC_ICE (FI1,X1,SFN11,SFN12 &
                         ,B11_MY,B12_MY,CF,ID,COL,NKR)

IMPLICIT NONE
INTEGER NRM,KR,ICE,ID,NKR
REAL B12,B11,FUN,DELM,FK,CF,SFN12S,SFN11S
REAL  COL,X1(NKR,ID),FI1(NKR,ID),B11_MY(NKR,ID),B12_MY(NKR,ID) &
     ,SFN11(ID),SFN12(ID)

NRM=NKR-1
DO ICE=1,ID  
   SFN11S=0.
   SFN12S=0.
   SFN11(ICE)=CF*SFN11S
   SFN12(ICE)=CF*SFN12S
   DO KR=1,NRM
! VALUE OF DISTRIBUTION function
      FK=FI1(KR,ICE)
      DELM=X1(KR,ICE)*3.*COL
! INTEGRAL'S EXPRESSION 
      FUN=FK*DELM
! VALUES OF INTEGRALS
      B11=B11_MY(KR,ICE)
      B12=B12_MY(KR,ICE)
      SFN11S=SFN11S+FUN*B11 
      SFN12S=SFN12S+FUN*B12
   ENDDO
   SFN11(ICE)=CF*SFN11S
   SFN12(ICE)=CF*SFN12S
ENDDO

RETURN
END SUBROUTINE JERTIMESC_ICE

!##############################################################################
Subroutine ONECOND2 (TT,QQ,PP,ROR  &
,VR2,VR3,VR4,VR5,PSINGLE &
,DEL1N,DEL2N,DIV1,DIV2 &
,FF2,PSI2,R2,RIEC,RO2BL &
,FF3,PSI3,R3,RSEC,RO3BL &
,FF4,PSI4,R4,RGEC,RO4BL &
,FF5,PSI5,R5,RHEC,RO5BL &
,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
,COL,DTCOND,ICEMAX,NKR &
,ISYM2,ISYM3,ISYM4,ISYM5)

IMPLICIT NONE
INTEGER NKR,ICEMAX
REAL VR2(NKR,ICEMAX),VR3(NKR),VR4(NKR) &
     ,VR5(NKR),PSINGLE &
     ,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
     ,FI1_K,FI2_K,FI3_K,FI4_K,FI5_K &
     ,R1_K,R2_K,R3_K,R4_K,R5_K &
     ,DTCOND,COL
INTEGER KR,ICE,ITIME,ICM,KCOND,NR,NRM, &
     ISYM2,ISYM3,ISYM4,ISYM5
REAL AL1,AL2,D,GAM,POD, &
     RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY,DT0LREF,DTLREF, &
     A1_MYN, BB1_MYN, A2_MYN, BB2_MYN,DT,DTT,XRAD, &
     TPC1, TPC2, TPC3, TPC4, TPC5, &
     EPSDEL, DT0L, DT0I,ROR,CWHUCM,B6,B8L,B8I, &
     DEL1,DEL2,DEL1S,DEL2S,TIMENEW,TIMEREV, &
     SFNL,SFNI,B5L,B5I,B7L,B7I,DOPL,DOPI,OPERQ,RW,RI,QW,PW, &
     PI,QI,DTNEWL,D1N,D2N,DT0L0,DT0I0,SFN31,SFN32,SFN52, &
     SFNII1,SFN21,SFN22,DTNEWI3,DTNEWI4,DTNEWI5,DTNEWI2_1, &
     DTNEWI2_2,DTNEWI1,SFN41,SFN42,DTNEWI2_3,DTNEWI2, &
     SFN51,SFNII2
REAL DT_WATER_COND,DT_WATER_EVAP,DT_ICE_COND,DT_ICE_EVAP, &
     DT_MIX_COND,DT_MIX_EVAP,DT_MIX_BERGERON,DT_MIX_ANTIBERGERON
INTEGER K
DOUBLE PRECISION DAL2
DOUBLE PRECISION RORI,TPN,TPS,QPN,QPS &
                 ,FI2R2,FI3R3,FI4R4,FI5R5 &
                 ,RMASSIAA,RMASSIBB &
                 ,ES1N,ES2N,EW1N,ARGEXP &
                 ,TT,QQ,PP &
                 ,DEL1N,DEL2N,DIV1,DIV2 &
                 ,OPER2,OPER3,AR1,AR2  
DOUBLE PRECISION DELMASSI1
! CRYSTALS
                                                                       
REAL R2(NKR,ICEMAX),RIEC(NKR,ICEMAX),RO2BL(NKR,ICEMAX) &
     ,FI2(NKR,ICEMAX),PSI2(NKR,ICEMAX),FF2(NKR,ICEMAX) &
     ,B21_MY(NKR,ICEMAX),B22_MY(NKR,ICEMAX)
! SNOW  
REAL R3(NKR),RSEC(NKR),RO3BL(NKR) &
     ,FI3(NKR),FF3(NKR),PSI3(NKR) &
     ,B31_MY(NKR),B32_MY(NKR)

! GRAUPELS
REAL R4(NKR),RGEC(NKR),RO4BL(NKR) &
     ,FI4(NKR),FF4(NKR),PSI4(NKR) &
     ,B41_MY(NKR),B42_MY(NKR)  

! HAIL  
REAL R5(NKR),RHEC(NKR),RO5BL(NKR) &
     ,FI5(NKR),FF5(NKR),PSI5(NKR) &
     ,B51_MY(NKR),B52_MY(NKR)  

REAL DEL2D(ICEMAX),SFNI1(ICEMAX),SFNI2(ICEMAX) &
    ,TIMESTEPD(NKR) 

! For IN regeneration
real totin_before, totin_after

!dummy budgets
real vapc,vapr
 
OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)

DATA AL1 /2500./, AL2 /2834./, D /0.211/ &
     ,GAM /1.E-4/, POD /10./ 
     
DATA RV_MY,CF_MY,D_MYIN,AL1_MY,AL2_MY &
     /461.5,0.24E-1,0.211E-4,2.5E6,2.834E6/

DATA A1_MYN, BB1_MYN, A2_MYN, BB2_MYN &
     /2.53,5.42,3.41E1,6.13/

DATA TPC1, TPC2, TPC3, TPC4, TPC5 &
     /-4.0,-8.1,-12.7,-17.8,-22.4/ 

DATA EPSDEL/0.1E-03/
    
DATA DT0L, DT0I /1.E20,1.E20/

! SOME CONSTANTS 
COL3=3.0*COL
ICM=ICEMAX
ITIME=0
KCOND=0
DT_WATER_COND=0.4
DT_WATER_EVAP=0.4
DT_ICE_COND=0.4
DT_ICE_EVAP=0.4
DT_MIX_COND=0.4
DT_MIX_EVAP=0.4
DT_MIX_BERGERON=0.4
DT_MIX_ANTIBERGERON=0.4
DT0LREF=0.2
DTLREF=0.4

NR=NKR
NRM=NKR-1
DT=DTCOND
DTT=DTCOND
XRAD=0.

CWHUCM=0.
XRAD=0.
B6=CWHUCM*GAM-XRAD
B8L=1./ROR
B8I=1./ROR
RORI=1./ROR

TPN=TT
QPN=QQ

DO ICE=1,ICEMAX
   SFNI1(ICE)=0.
   SFNI2(ICE)=0.
   DEL2D(ICE)=0.
ENDDO

TIMENEW=0.
ITIME=0
! IN regeneration: calculate the total ice particles before cond - J. Fan
totin_before = 0.0
do kr = 1, nkr
   totin_before = totin_before+psi2(kr,1)*r2(kr,1)*col3  &
                              +psi2(kr,2)*r2(kr,2)*col3  &
                              +psi2(kr,3)*r2(kr,3)*col3  &
                              +psi3(kr)*r3(kr)*col3      &
                              +psi4(kr)*r4(kr)*col3      &
                              +psi5(kr)*r5(kr)*col3
enddo

!-----------Start small time step loop----------------------------
46 ITIME=ITIME+1

   TIMEREV=DT-TIMENEW

   DEL1=DEL1N
   DEL2=DEL2N
   DEL1S=DEL1N
   DEL2S=DEL2N
   DEL2D(1)=DEL2N
   DEL2D(2)=DEL2N
   DEL2D(3)=DEL2N
   TPS=TPN
   QPS=QPN
   DO KR=1,NKR
      FI3(KR)=PSI3(KR)
      FI4(KR)=PSI4(KR)
      FI5(KR)=PSI5(KR)
      DO ICE=1,ICEMAX
         FI2(KR,ICE)=PSI2(KR,ICE)
      ENDDO
   ENDDO
! TIME-STEP GROWTH RATE: 
   CALL JERRATE (R2,TPS,PP,ROR,VR2,PSINGLE &
                ,RIEC,RO2BL,B21_MY,B22_MY,3,2,ICEMAX,NKR)   
   CALL JERRATE (R3,TPS,PP,ROR,VR3,PSINGLE &
                ,RSEC,RO3BL,B31_MY,B32_MY,1,2,ICEMAX,NKR)
   CALL JERRATE (R4,TPS,PP,ROR,VR4,PSINGLE &
                ,RGEC,RO4BL,B41_MY,B42_MY,1,2,ICEMAX,NKR)
   CALL JERRATE (R5,TPS,PP,ROR,VR5,PSINGLE &
                ,RHEC,RO5BL,B51_MY,B52_MY,1,2,ICEMAX,NKR)

! INTEGRALS IN DELTA EQUATION

   CALL JERTIMESC_ICE (FI2,R2,SFNI1,SFNI2,B21_MY,B22_MY,B8I,ICM,COL,NKR) 
   CALL JERTIMESC (FI3,R3,SFN31,SFN32,B31_MY,B32_MY,B8I,1,COL,NKR)  
   CALL JERTIMESC (FI4,R4,SFN41,SFN42,B41_MY,B42_MY,B8I,1,COL,NKR) 
   CALL JERTIMESC (FI5,R5,SFN51,SFN52,B51_MY,B52_MY,B8I,1,COL,NKR)
   SFNII1=SFNI1(1)+SFNI1(2)+SFNI1(3)
   SFNII2=SFNI2(1)+SFNI2(2)+SFNI2(3)
   SFN21=SFNII1+SFN31+SFN41+SFN51  
   SFN22=SFNII2+SFN32+SFN42+SFN52 
   SFNL=0.
   SFNI=SFN21+SFN22 
! SOME CONSTANTS 
   B5L=BB1_MY/TPS/TPS
   B5I=BB2_MY/TPS/TPS
   B7L=B5L*B6     
   B7I=B5I*B6
   DOPL=1.+DEL1S     
   DOPI=1.+DEL2S     
   OPERQ=OPER2(QPS)  
   RW=(OPERQ+B5L*AL1)*DOPL*SFNL  
   QW=B7L*DOPL
   PW=(OPERQ+B5I*AL1)*DOPI*SFNL
   RI=(OPERQ+B5L*AL2)*DOPL*SFNI
   PI=(OPERQ+B5I*AL2)*DOPI*SFNI
   QI=B7I*DOPI
   KCOND=20

   IF(DEL2.GT.0) KCOND=21

   IF(KCOND.EQ.21)  THEN ! ONLY_ICE: CONDENSATION
      DT0I=1.E20
      DTNEWI1=DTCOND
      DTNEWL=DTNEWI1
      IF(ITIME.GE.NKR) THEN
         PRINT *, 'ONLY_ICE: CONDENSATION'
         STOP
      ENDIF
      TIMESTEPD(ITIME)=DTNEWL
! NEW TIME STEP (ONLY_ICE: CONDENSATION)
      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1))  &
         DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMESTEPD(ITIME)=DTNEWL
      TIMENEW=TIMENEW+DTNEWL
      DTT=DTNEWL

! SOLVING FOR SUPERSATURATION (ONLY ICE: CONDENSATION) 
      CALL JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                     ,RW,PW,RI,PI,QW,QI &
                     ,DTT,D1N,D2N,DT0L0,DT0I0)

      IF(ISYM2.NE.0) THEN
        CALL JERDFUN (R2,B21_MY,B22_MY &
                     ,FI2,PSI2,D2N &
                     ,ICM,1,COL,NKR,vapc,vapr,rori)

        CALL JERDFUN (R2,B21_MY,B22_MY &
                     ,FI2,PSI2,D2N &
                     ,ICM,2,COL,NKR,vapc,vapr,rori)

        CALL JERDFUN (R2,B21_MY,B22_MY &
                     ,FI2,PSI2,D2N &
                     ,ICM,3,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM3.NE.0) THEN !Snow
        CALL JERDFUN (R3,B31_MY,B32_MY &
                     ,FI3,PSI3,D2N &
                     ,1,4,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM4.NE.0) THEN  !Graupel
        CALL JERDFUN (R4,B41_MY,B42_MY &
                     ,FI4,PSI4,D2N &
                     ,1,4,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM5.NE.0) THEN !Hail
        CALL JERDFUN (R5,B51_MY,B52_MY &
                     ,FI5,PSI5,D2N &
                     ,1,5,COL,NKR,vapc,vapr,rori)
      ENDIF
   ELSE
      DT0I=1.E20
      IF (DEL2N.EQ.0)THEN
         DTNEWL=DT
      ELSE
         DTNEWI3=-R3(3)/(B31_MY(3)*DEL2N-B32_MY(3))
         DTNEWI4=-R4(3)/(B41_MY(3)*DEL2N-B42_MY(3))
         DTNEWI5=-R5(3)/(B51_MY(3)*DEL2N-B52_MY(3))
         DTNEWI2_1=-R2(3,1)/(B21_MY(1,1)*DEL2N-B22_MY(1,1))
         DTNEWI2_2=-R2(3,2)/(B21_MY(1,2)*DEL2N-B22_MY(1,2))
         DTNEWI2_3=-R2(3,3)/(B21_MY(1,3)*DEL2N-B22_MY(1,3))
         DTNEWI2=AMIN1(DTNEWI2_1,DTNEWI2_2,DTNEWI2_3)
         DTNEWI1=AMIN1(DTNEWI2,DTNEWI3,DTNEWI4 &
                      ,DTNEWI5,DT0I,TIMEREV)
         DTNEWI1=AMIN1(DTNEWI2,DTNEWI4,DTNEWI5,DT0I,TIMEREV)
         DTNEWL=DTNEWI1
         IF(DTNEWL.LT.DTLREF) DTNEWL=AMIN1(DTLREF,TIMEREV)
      END IF
      IF(ITIME.GE.NKR) THEN
         PRINT *, 'ONLY_ICE: EVAPORATION'
         STOP
      ENDIF
      TIMESTEPD(ITIME)=DTNEWL

      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1))  &
         DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMENEW=TIMENEW+DTNEWL
      TIMESTEPD(ITIME)=DTNEWL
      DTT=DTNEWL
! SOLVING FOR SUPERSATURATION (ONLY_ICE: EVAPORATION) 
      CALL JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                     ,RW,PW,RI,PI,QW,QI &
                     ,DTT,D1N,D2N,DT0L0,DT0I0)
      IF(ISYM2.NE.0) THEN !CRYSTALS
          CALL JERDFUN (R2,B21_MY,B22_MY &
                       ,FI2,PSI2,D2N &
                       ,ICM,1,COL,NKR,vapc,vapr,rori)

          CALL JERDFUN (R2,B21_MY,B22_MY &
                       ,FI2,PSI2,D2N &
                       ,ICM,2,COL,NKR,vapc,vapr,rori)

          CALL JERDFUN (R2,B21_MY,B22_MY &
                       ,FI2,PSI2,D2N &
                       ,ICM,3,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM3.NE.0) THEN !SNOW
         CALL JERDFUN (R3,B31_MY,B32_MY &
                      ,FI3,PSI3,D2N &
                      ,1,3,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM4.NE.0) THEN !GRAUPEL
          CALL JERDFUN (R4,B41_MY,B42_MY &
                       ,FI4,PSI4,D2N &
                       ,1,4,COL,NKR,vapc,vapr,rori)
      ENDIF
      IF(ISYM5.NE.0) THEN !HAIL
          CALL JERDFUN (R5,B51_MY,B52_MY &
                       ,FI5,PSI5,D2N &
                       ,1,5,COL,NKR,vapc,vapr,rori)
      ENDIF
   ENDIF

! MASSES
   RMASSIBB=0.0
   RMASSIAA=0.0
! BEFORE JERNEWF
   DO K=1,NKR
      DO ICE =1,ICEMAX
         FI2_K=FI2(K,ICE)
         R2_K=R2(K,ICE)
         FI2R2=FI2_K*R2_K*R2_K
         RMASSIBB=RMASSIBB+FI2R2
      ENDDO
      FI3_K=FI3(K)
      FI4_K=FI4(K)
      FI5_K=FI5(K)
      R3_K=R3(K)
      R4_K=R4(K)
      R5_K=R5(K)
      FI3R3=FI3_K*R3_K*R3_K
      FI4R4=FI4_K*R4_K*R4_K
      FI5R5=FI5_K*R5_K*R5_K
      RMASSIBB=RMASSIBB+FI3R3
      RMASSIBB=RMASSIBB+FI4R4
      RMASSIBB=RMASSIBB+FI5R5
   ENDDO
   RMASSIBB=RMASSIBB*COL3*RORI
   IF(RMASSIBB.LT.0.0) RMASSIBB=0.0
   DO K=1,NKR
      DO ICE =1,ICEMAX
         FI2_K=PSI2(K,ICE)
         R2_K=R2(K,ICE)
         FI2R2=FI2_K*R2_K*R2_K
         RMASSIAA=RMASSIAA+FI2R2
      ENDDO
      FI3_K=PSI3(K)
      FI4_K=PSI4(K)
      FI5_K=PSI5(K)
      R3_K=R3(K)
      R4_K=R4(K)
      R5_K=R5(K)
      FI3R3=FI3_K*R3_K*R3_K
      FI4R4=FI4_K*R4_K*R4_K
      FI5R5=FI5_K*R5_K*R5_K
      RMASSIAA=RMASSIAA+FI3R3
      RMASSIAA=RMASSIAA+FI4R4
      RMASSIAA=RMASSIAA+FI5R5
   ENDDO
   RMASSIAA=RMASSIAA*COL3*RORI
   IF(RMASSIAA.LT.0.0) RMASSIAA=0.0
   DELMASSI1=RMASSIAA-RMASSIBB
   QPN=QPS-DELMASSI1
   DAL2=AL2
   TPN=TPS+DAL2*DELMASSI1
! SUPERSATURATION
   ARGEXP=-BB1_MY/TPN
   ES1N=AA1_MY*DEXP(ARGEXP)
   ARGEXP=-BB2_MY/TPN
   ES2N=AA2_MY*DEXP(ARGEXP)
   EW1N=OPER3(QPN,PP)
   IF(ES1N.EQ.0)THEN
    DEL1N=0.5
    DIV1=1.5
    print*,'es1n onecond2 = 0'
    stop
   ELSE
    DIV1=EW1N/ES1N
    DEL1N=EW1N/ES1N-1.
   END IF
   IF(ES2N.EQ.0)THEN
    DEL2N=0.5
    DIV2=1.5
    print*,'es2n onecond2 = 0'
    stop
   ELSE
    DEL2N=EW1N/ES2N-1.
    DIV2=EW1N/ES2N
   END IF

!  END OF TIME SPLITTING 
IF(TIMENEW.LT.DT) GOTO 46

! IN regeneration: calculate the total ice particles after cond - J. Fan
totin_after = 0.0
do kr = 1, nkr
totin_after = totin_after+psi2(kr,1)*r2(kr,1)*3.*col  &
                          +psi2(kr,2)*r2(kr,2)*3.*col  &
                          +psi2(kr,3)*r2(kr,3)*3.*col  &
                          +psi3(kr)*r3(kr)*3.*col      &
                          +psi4(kr)*r4(kr)*3.*col      &
                          +psi5(kr)*r5(kr)*3.*col
enddo
inreg = max((totin_before - totin_after),0.0)

TT=TPN
QQ=QPN
DO KR=1,NKR
   DO ICE=1,ICEMAX
      FF2(KR,ICE)=PSI2(KR,ICE)
   ENDDO
   FF3(KR)=PSI3(KR)
   FF4(KR)=PSI4(KR)
   FF5(KR)=PSI5(KR)
ENDDO

RETURN                                          
END SUBROUTINE ONECOND2

!##############################################################################
Subroutine ONECOND3 (TT,QQ,PP,ROR &
,VR1,VR2,VR3,VR4,VR5,PSINGLE &
,DEL1N,DEL2N,DIV1,DIV2 &
,FF1,PSI1,R1,RLEC,RO1BL &
,FF2,PSI2,R2,RIEC,RO2BL &
,FF3,PSI3,R3,RSEC,RO3BL &
,FF4,PSI4,R4,RGEC,RO4BL &
,FF5,PSI5,R5,RHEC,RO5BL &
,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
,COL,DTCOND,ICEMAX,NKR &
,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5 &
,vapcldt,vapraint)

IMPLICIT NONE
INTEGER NKR,KR,ITIME,ICE,KCOND,K,ICEMAX &
        ,ISYM1,ISYM2,ISYM3,ISYM4,ISYM5
INTEGER KLIMITL,KLIMITG,KLIMITH,KLIMITI_1,KLIMITI_2,KLIMITI_3
REAL ROR,VR1(NKR),VR2(NKR,ICEMAX),VR3(NKR),VR4(NKR) &
     ,VR5(NKR),PSINGLE,AA1_MY,BB1_MY,AA2_MY,BB2_MY &
     ,DTCOND,COL

! DROPLETS    
REAL R1(NKR),RLEC(NKR),RO1BL(NKR) &
    ,FI1(NKR),FF1(NKR),PSI1(NKR) &
    ,B11_MY(NKR),B12_MY(NKR)

! CRYSTALS
REAL R2(NKR,ICEMAX) &
    ,RIEC(NKR,ICEMAX) &
    ,RO2BL(NKR,ICEMAX) &
    ,FI2(NKR,ICEMAX),PSI2(NKR,ICEMAX) &
    ,FF2(NKR,ICEMAX) &
    ,B21_MY(NKR,ICEMAX),B22_MY(NKR,ICEMAX) 

! SNOW  
REAL R3(NKR) &
    ,RSEC(NKR),RO3BL(NKR) &
    ,FI3(NKR),FF3(NKR),PSI3(NKR) &
    ,B31_MY(NKR),B32_MY(NKR) 

! GRAUPELS 
REAL R4(NKR) &
    ,RGEC(NKR),RO4BL(NKR) &
    ,FI4(NKR),FF4(NKR),PSI4(NKR) &
    ,B41_MY(NKR),B42_MY(NKR) 

! HAIL  
REAL R5(NKR) &
    ,RHEC(NKR),RO5BL(NKR) &
    ,FI5(NKR),FF5(NKR),PSI5(NKR) &
    ,B51_MY(NKR),B52_MY(NKR) 

DOUBLE PRECISION DAL1,DAL2
DOUBLE PRECISION RORI,TPN,TPS,QPN,QPS &
                 ,FI1R1,FI2R2,FI3R3,FI4R4,FI5R5 &
                 ,FI1_K,FI2_K,FI3_K,FI4_K,FI5_K &
                 ,R1_K,R2_K,R3_K,R4_K,R5_K &
                 ,RMASSLAA,RMASSLBB,RMASSIAA,RMASSIBB &
                 ,ES1N,ES2N,EW1N,ARGEXP &
                 ,TT,QQ,PP,DEL1N0,DEL2N0 &
                 ,DEL1N,DEL2N,DIV1,DIV2 &
                 ,OPER2,OPER3,AR1,AR2

DOUBLE PRECISION DELTAQ1,DELMASSI1,DELMASSL1

REAL A1_MYN, BB1_MYN, A2_MYN, BB2_MYN
DATA A1_MYN, BB1_MYN, A2_MYN, BB2_MYN &
     /2.53,5.42,3.41E1,6.13/
REAL B8L,B8I,SFN11,SFN12,SFNL,SFNI
REAL B5L,B5I,B7L,B7I,B6,DOPL,DEL1S,DEL2S,DOPI,RW,QW,PW, &
 RI,PI,QI,SFNI1(ICEMAX),SFNI2(ICEMAX),AL1,AL2
REAL D1N,D2N,DT0L, DT0I,D1N0,D2N0
REAL SFN21,SFN22,SFNII1,SFNII2,SFN31,SFN32,SFN41,SFN42,SFN51,SFN52
REAL DEL1,DEL2
REAL  TIMEREV,DT,DTT,TIMENEW
REAL DTIMEG(NKR),DTIMEH(NKR)

REAL DEL2D(ICEMAX),DTIMEL(NKR) &
    ,DTIMEI_1(NKR),DTIMEI_2(NKR),DTIMEI_3(NKR)
REAL DT_WATER_COND,DT_WATER_EVAP,DT_ICE_COND,DT_ICE_EVAP, &
 DT_MIX_COND,DT_MIX_EVAP,DT_MIX_BERGERON,DT_MIX_ANTIBERGERON
REAL DTNEWL0,DTNEWL1,DTNEWI1,DTNEWI2_1,DTNEWI2_2,DTNEWI2_3, &
DTNEWI2,DTNEWI3,DTNEWI4,DTNEWI5, &
DTNEWL,DTNEWL2,DTNEWG1,DTNEWH1
REAL TIMESTEPD(NKR)

! CCN/IN regeneration
real totccn_before, totccn_after, totin_before, totin_after

!condensation budgets
real vapcldt, vapraint,vapc,vapr

DATA AL1 /2500./, AL2 /2834./
REAL EPSDEL,EPSDEL2
DATA EPSDEL, EPSDEL2 /0.1E-03,0.1E-03/

OPER2(AR1)=0.622/(0.622+0.378*AR1)/AR1
OPER3(AR1,AR2)=AR1*AR2/(0.622+0.378*AR1)

vapc=0.
vapr=0.

DT_WATER_COND=0.4
DT_WATER_EVAP=0.4
DT_ICE_COND=0.4
DT_ICE_EVAP=0.4
DT_MIX_COND=0.4
DT_MIX_EVAP=0.4
DT_MIX_BERGERON=0.4
DT_MIX_ANTIBERGERON=0.4

ITIME = 0
TIMENEW=0.
DT=DTCOND
DTT=DTCOND

B6=0.
B8L=1./ROR
B8I=1./ROR
RORI=1.D0/ROR
TPN=TT
QPN=QQ
! Calculate total CCN and IN before cond - J. Fan
totccn_before = 0.0
totin_before = 0.0
do kr = 1, nkr
   totccn_before = totccn_before+psi1(kr)*r1(kr)*col3
   totin_before = totin_before+psi2(kr,1)*r2(kr,1)*col3  &
                              +psi2(kr,2)*r2(kr,2)*col3  &
                              +psi2(kr,3)*r2(kr,3)*col3  &
                              +psi3(kr)*r3(kr)*col3      &
                              +psi4(kr)*r4(kr)*col3      &
                              +psi5(kr)*r5(kr)*col3
enddo

!-------------Time-splitting Loop begins here-----------------------
16 ITIME=ITIME+1
   IF((TPN-273.15).GE.-0.187) GO TO 17
   TIMEREV=DT-TIMENEW
   DEL1=DEL1N
   DEL2=DEL2N
   DEL1S=DEL1N
   DEL2S=DEL2N
! NEW ALGORITHM (NO TYPE ICE)
   DEL2D(1)=DEL2N
   DEL2D(2)=DEL2N
   DEL2D(3)=DEL2N
   TPS=TPN
   QPS=QPN
   DO KR=1,NKR
      FI1(KR)=PSI1(KR)
      FI3(KR)=PSI3(KR)
      FI4(KR)=PSI4(KR)
      FI5(KR)=PSI5(KR)
      DO ICE=1,ICEMAX
         FI2(KR,ICE)=PSI2(KR,ICE)
      ENDDO
   ENDDO
! TIME-STEP GROWTH RATE
   CALL JERRATE (R1,TPS,PP,ROR,VR1,PSINGLE &
               ,RLEC,RO1BL,B11_MY,B12_MY,1,1,ICEMAX,NKR)
   CALL JERRATE (R2,TPS,PP,ROR,VR2,PSINGLE &
               ,RIEC,RO2BL,B21_MY,B22_MY,3,2,ICEMAX,NKR)
   CALL JERRATE (R3,TPS,PP,ROR,VR3,PSINGLE &
               ,RSEC,RO3BL,B31_MY,B32_MY,1,2,ICEMAX,NKR)
   CALL JERRATE (R4,TPS,PP,ROR,VR4,PSINGLE &
               ,RGEC,RO4BL,B41_MY,B42_MY,1,2,ICEMAX,NKR)
   CALL JERRATE (R5,TPS,PP,ROR,VR5,PSINGLE &
               ,RHEC,RO5BL,B51_MY,B52_MY,1,2,ICEMAX,NKR)

   CALL JERTIMESC (FI1,R1,SFN11,SFN12 &
                  ,B11_MY,B12_MY,B8L,1,COL,NKR)
   CALL JERTIMESC_ICE (FI2,R2,SFNI1,SFNI2 &
                      ,B21_MY,B22_MY,B8I,ICEMAX,COL,NKR)
   CALL JERTIMESC (FI3,R3,SFN31,SFN32 &
                  ,B31_MY,B32_MY,B8I,1,COL,NKR)
   CALL JERTIMESC (FI4,R4,SFN41,SFN42 &
                  ,B41_MY,B42_MY,B8I,1,COL,NKR)
   CALL JERTIMESC (FI5,R5,SFN51,SFN52 &
                  ,B51_MY,B52_MY,B8I,1,COL,NKR)

   SFNII1=SFNI1(1)+SFNI1(2)+SFNI1(3)
   SFNII2=SFNI2(1)+SFNI2(2)+SFNI2(3)
   SFN21=SFNII1+SFN31+SFN41+SFN51
   SFN22=SFNII2+SFN32+SFN42+SFN52
   SFNL=SFN11+SFN12
   SFNI=SFN21+SFN22
! SOME CONSTANTS (QW,QI=0,since B6=0.)
   B5L=BB1_MY/TPS/TPS
   B5I=BB2_MY/TPS/TPS
   B7L=B5L*B6
   B7I=B5I*B6
   DOPL=1.+DEL1S
   DOPI=1.+DEL2S
   RW=(OPER2(QPS)+B5L*AL1)*DOPL*SFNL
   QW=B7L*DOPL
   PW=(OPER2(QPS)+B5I*AL1)*DOPI*SFNL
   RI=(OPER2(QPS)+B5L*AL2)*DOPL*SFNI
   PI=(OPER2(QPS)+B5I*AL2)*DOPI*SFNI
   QI=B7I*DOPI
! SOLVING FOR TIMEZERO
   CALL JERSUPSAT (DEL1,DEL2,DEL1N0,DEL2N0 &
                  ,RW,PW,RI,PI,QW,QI &
                  ,DTT,D1N0,D2N0,DT0L,DT0I)
! DEL1 > 0, DEL2 < 0    (ANTIBERGERON MIXED PHASE - KCOND=50)
! DEL1 < 0 AND DEL2 < 0 (EVAPORATION MIXED_PHASE - KCOND=30)
! DEL1 > 0 AND DEL2 > 0 (CONDENSATION MIXED PHASE - KCOND=31)
! DEL1 < 0, DEL2 > 0    (BERGERON MIXED PHASE - KCOND=32)
   KCOND=50
   IF(DEL1.LT.0.AND.DEL2.LT.0) KCOND=30
   IF(DEL1.GT.0.AND.DEL2.GT.0) KCOND=31
   IF(DEL1.LT.0.AND.DEL2.GT.0) KCOND=32

!----------SSw>0, SSi<0----------------------------------
!Adele-Physically impossible, why does this section exist?
   IF(KCOND.EQ.50) THEN 
      IF(DT0L.EQ.0) THEN
         DTNEWL=DT
      ELSE
         DTNEWL=AMIN1(DT,DT0L)
      ENDIF
      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1)) &
          DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMENEW=TIMENEW+DTNEWL
      DTT=DTNEWL
      IF(ITIME.GE.NKR) THEN
         PRINT *, 'ANTIBERGERON MIXED PHASE'
         STOP
      ENDIF
      TIMESTEPD(ITIME)=DTNEWL
   ENDIF
!------------SSw>0, SSi>0  Everything grows -----------------------------
   IF(KCOND.EQ.31) THEN
      IF (DEL1N.EQ.0)THEN
         DTNEWL0=DT
      ELSE
         DTNEWL0=ABS(R1(ITIME)/(B11_MY(ITIME)*DEL1N- &
                                B12_MY(ITIME)))
      END IF

      IF (DEL2N.EQ.0)THEN
         DTNEWI2_1=DT
         DTNEWI2_2=DT
         DTNEWI2_3=DT
         DTNEWI2=DT
         DTNEWI3=DT
         DTNEWI4=DT
         DTNEWI5=DT
      ELSE
         DTNEWI2_1=ABS(R2(ITIME,1)/ &
          (B21_MY(ITIME,1)*DEL2N-B22_MY(ITIME,1)))
         DTNEWI2_2=ABS(R2(ITIME,2)/ &
          (B21_MY(ITIME,2)*DEL2N-B22_MY(ITIME,2))) 
         DTNEWI2_3=ABS(R2(ITIME,3)/ &
          (B21_MY(ITIME,3)*DEL2N-B22_MY(ITIME,3)))  
         DTNEWI2=AMIN1(DTNEWI2_1,DTNEWI2_2,DTNEWI2_3)

         DTNEWI3=ABS(R3(ITIME)/(B31_MY(ITIME)*DEL2N- &
                                B32_MY(ITIME)))
         DTNEWI4=ABS(R4(ITIME)/(B41_MY(ITIME)*DEL2N- &
                                B42_MY(ITIME)))
         DTNEWI5=ABS(R5(ITIME)/(B51_MY(ITIME)*DEL2N- &
                                B52_MY(ITIME)))
      END IF
      DTNEWI1=AMIN1(DTNEWI2,DTNEWI4,DTNEWI5,DT0I)
      IF(DT0L.NE.0) THEN
         IF(ABS(DT0L).LT.DT_MIX_COND) THEN
            DTNEWL1=AMIN1(DT_MIX_COND,DTNEWL0)
         ELSE
            DTNEWL1=AMIN1(DT0L,DTNEWL0)
         ENDIF
      ELSE
         DTNEWL1=DTNEWL0
      ENDIF
      DTNEWL=AMIN1(DTNEWL1,DTNEWI1)
      IF(ITIME.GE.NKR) THEN
         PRINT *, 'CONDENSATION MIXED PHASE'
         STOP
      ENDIF
      TIMESTEPD(ITIME)=DTNEWL
      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1)) &
         DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMENEW=TIMENEW+DTNEWL
      TIMESTEPD(ITIME)=DTNEWL
      DTT=DTNEWL
   ENDIF
!-------------SSw<0, SSi<0  All evaporation -------------------------------
   IF(KCOND.EQ.30) THEN
      DO KR=1,NKR
         DTIMEL(KR)=0.
         DTIMEG(KR)=0.
         DTIMEH(KR)=0.
         DTIMEI_1(KR)=0.
         DTIMEI_2(KR)=0.
         DTIMEI_3(KR)=0.
      ENDDO
      DO KR=1,NKR
         IF (DEL1N.EQ.0) THEN
            DTIMEL(KR)=DT
            DTIMEG(KR)=DT
            DTIMEH(KR)=DT
         ELSE
            DTIMEL(KR)=-R1(KR)/(B11_MY(KR)*DEL1N- &
                                B12_MY(KR))
            DTIMEG(KR)=-R4(KR)/(B41_MY(KR)*DEL1N- &
                                B42_MY(KR))
            DTIMEH(KR)=-R5(KR)/(B51_MY(KR)*DEL1N- &
                            B52_MY(KR))
! NEW ALGORITHM (NO TYPE OF ICE)
         END IF
         IF (DEL2N.EQ.0) THEN
            DTIMEI_1(KR)=DT
            DTIMEI_2(KR)=DT
            DTIMEI_3(KR)=DT
         ELSE
            DTIMEI_1(KR)=-R2(KR,1)/ &
                 (B21_MY(KR,1)*DEL2N-B22_MY(KR,1))
            DTIMEI_2(KR)=-R2(KR,2)/ &
                 (B21_MY(KR,2)*DEL2N-B22_MY(KR,2))
            DTIMEI_3(KR)=-R2(KR,3)/ &
                 (B21_MY(KR,3)*DEL2N-B22_MY(KR,3))
         END IF
      ENDDO
! WATER
      KLIMITL=1
      DO KR=1,NKR
         IF(DTIMEL(KR).GT.TIMEREV) EXIT
         KLIMITL=KR
      ENDDO
      KLIMITL=KLIMITL-1
      IF(KLIMITL.LT.1) KLIMITL=1
      DTNEWL1=AMIN1(DTIMEL(KLIMITL),DT0L,TIMEREV)

! GRAUPELS
      KLIMITG=1
      DO KR=1,NKR
         IF(DTIMEG(KR).GT.TIMEREV) EXIT
         KLIMITG=KR
      ENDDO
      KLIMITG=KLIMITG-1
      IF(KLIMITG.LT.1) KLIMITG=1
      DTNEWG1=AMIN1(DTIMEG(KLIMITG),TIMEREV)

! HAIL
      KLIMITH=1
      DO KR=1,NKR
         IF(DTIMEH(KR).GT.TIMEREV) EXIT
         KLIMITH=KR
      ENDDO
      KLIMITH=KLIMITH-1
      IF(KLIMITH.LT.1) KLIMITH=1
      DTNEWH1=AMIN1(DTIMEH(KLIMITH),TIMEREV)

! ICE CRYSTALS
      KLIMITI_1=1
      KLIMITI_2=1
      KLIMITI_3=1
      DO KR=1,NKR
         IF(DTIMEI_1(KR).GT.TIMEREV) EXIT
         KLIMITI_1=KR
      ENDDO
 
      DO KR=1,NKR
         IF(DTIMEI_2(KR).GT.TIMEREV) EXIT
         KLIMITI_2=KR
      ENDDO

      DO KR=1,NKR
         IF(DTIMEI_3(KR).GT.TIMEREV) EXIT
         KLIMITI_3=KR
      ENDDO

      KLIMITI_1=KLIMITI_1-1
      IF(KLIMITI_1.LT.1) KLIMITI_1=1
      DTNEWI2_1=AMIN1(DTIMEI_1(KLIMITI_1),TIMEREV)

      KLIMITI_2=KLIMITI_2-1
      IF(KLIMITI_2.LT.1) KLIMITI_2=1
      DTNEWI2_2=AMIN1(DTIMEI_2(KLIMITI_2),TIMEREV)

      KLIMITI_3=KLIMITI_3-1
      IF(KLIMITI_3.LT.1) KLIMITI_3=1
      DTNEWI2_3=AMIN1(DTIMEI_3(KLIMITI_3),TIMEREV)

      DTNEWI2=AMIN1(DTNEWI2_1,DTNEWI2_2,DTNEWI2_3)
      DTNEWI1=AMIN1(DTNEWI2,DTNEWG1,DTNEWH1,DT0I)

      IF(ABS(DEL2N).LT.EPSDEL2) &
         DTNEWI1=AMIN1(DTNEWI2,DTNEWG1,DTNEWH1)
      DTNEWL2=AMIN1(DTNEWL1,DTNEWI1)
      DTNEWL=DTNEWL2
      IF(DTNEWL.LT.DT_MIX_EVAP) &
         DTNEWL=AMIN1(DT_MIX_EVAP,TIMEREV)  
      IF(ITIME.GE.NKR) THEN
         PRINT *, 'EVAPORATION MIXED PHASE'
         STOP
      ENDIF

      TIMESTEPD(ITIME)=DTNEWL
      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1)) &
          DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMESTEPD(ITIME)=DTNEWL
      TIMENEW=TIMENEW+DTNEWL
      DTT=DTNEWL
   ENDIF

!-------------SSw<0, SSi>0    Bergeron case  ------------------------------------
   IF(KCOND.EQ.32) THEN
      IF (DEL1N.EQ.0)THEN
         DTNEWL0=DT
      ELSE
         DTNEWL0=-R1(1)/(B11_MY(1)*DEL1N-B12_MY(1))
      END IF

      IF (DEL2N.EQ.0)THEN
         DTNEWI2_1=DT
         DTNEWI2_2=DT
         DTNEWI2_3=DT
      ELSE
         DTNEWI2_1=R2(1,1)/(B21_MY(1,1)*DEL2N-B22_MY(1,1))
         DTNEWI2_2=R2(1,2)/(B21_MY(1,2)*DEL2N-B22_MY(1,2))
         DTNEWI2_3=R2(1,3)/(B21_MY(1,3)*DEL2N-B22_MY(1,3))
      END IF
         DTNEWI2=AMIN1(DTNEWI2_1,DTNEWI2_2,DTNEWI2_3)
      IF (DEL2N.EQ.0)THEN
         DTNEWI3=DT
         DTNEWI4=DT
         DTNEWI5=DT
      ELSE
         DTNEWI3=R3(1)/(B31_MY(1)*DEL2N-B32_MY(1))
         DTNEWI4=R4(1)/(B41_MY(1)*DEL2N-B42_MY(1))
         DTNEWI5=R5(1)/(B51_MY(1)*DEL2N-B52_MY(1))
      END IF
      DTNEWL1=AMIN1(DTNEWL0,DT0L,TIMEREV)
      DTNEWI1=AMIN1(DTNEWI2,DTNEWI3,DTNEWI4 &
                    ,DTNEWI5,DT0I,TIMEREV)
      DTNEWI1=AMIN1(DTNEWI2,DTNEWI4,DTNEWI5,DT0I,TIMEREV)
      DTNEWL=AMIN1(DTNEWL1,DTNEWI1)

      IF(DTNEWL.LT.DT_MIX_BERGERON) &
         DTNEWL=AMIN1(DT_MIX_BERGERON,TIMEREV)
      TIMESTEPD(ITIME)=DTNEWL

      IF(DTNEWL.GT.DT) DTNEWL=DT
      IF((TIMENEW+DTNEWL).GT.DT.AND.ITIME.LT.(NKR-1)) &
         DTNEWL=DT-TIMENEW
      IF(ITIME.EQ.(NKR-1)) DTNEWL=DT-TIMENEW
      TIMESTEPD(ITIME)=DTNEWL
      TIMENEW=TIMENEW+DTNEWL
      DTT=DTNEWL
   ENDIF
!------------end kcond cases----------------------------------------------

! SOLVING FOR SUPERSATURATION 
         
   CALL JERSUPSAT (DEL1,DEL2,DEL1N,DEL2N &
                  ,RW,PW,RI,PI,QW,QI &
                  ,DTT,D1N,D2N,DT0L,DT0I)

!------------- DROPLETS----------------------------------------------- 
   IF(ISYM1.NE.0) THEN
      CALL JERDFUN (R1,B11_MY,B12_MY &
                   ,FI1,PSI1,D1N &
                   ,1,6,COL,NKR,vapcldt,vapraint,rori)
   ENDIF                     
!------------ CRYSTALS ---------------------------------------------
   IF(ISYM2.NE.0) THEN
      CALL JERDFUN (R2,B21_MY,B22_MY &
                   ,FI2,PSI2,D2N &
                   ,ICEMAX,1,COL,NKR,vapc,vapr,rori)

      CALL JERDFUN (R2,B21_MY,B22_MY &
                   ,FI2,PSI2,D2N &
                   ,ICEMAX,2,COL,NKR,vapc,vapr,rori)

      CALL JERDFUN (R2,B21_MY,B22_MY &
                   ,FI2,PSI2,D2N &
                   ,ICEMAX,3,COL,NKR,vapc,vapr,rori)
   ENDIF
! ---------------SNOW -----------------------------------------
   IF(ISYM3.NE.0) THEN
      CALL JERDFUN (R3,B31_MY,B32_MY &
                   ,FI3,PSI3,D2N &
                   ,1,3,COL,NKR,vapc,vapr,rori)
   ENDIF
!------------- GRAUPELS -------------------------------------
   IF(ISYM4.NE.0) THEN
      CALL JERDFUN (R4,B41_MY,B42_MY &
                   ,FI4,PSI4,D2N &
                   ,1,4,COL,NKR,vapc,vapr,rori)
   ENDIF
! ---------------HAIL --------------------------------------
   IF(ISYM5.NE.0) THEN
      CALL JERDFUN (R5,B51_MY,B52_MY &
                   ,FI5,PSI5,D2N &
                   ,1,5,COL,NKR,vapc,vapr,rori)
   ENDIF

!------------Update other variables-----------------------------
! MASSES
   RMASSLBB=0.D0
   RMASSIBB=0.D0
   RMASSLAA=0.D0
   RMASSIAA=0.D0
! BEFORE JERNEWF
   DO K=1,NKR
      FI1_K=FI1(K)
      R1_K=R1(K)
      FI1R1=FI1_K*R1_K*R1_K
      RMASSLBB=RMASSLBB+FI1R1
      DO ICE =1,ICEMAX
         FI2_K=FI2(K,ICE)
         R2_K=R2(K,ICE)
         FI2R2=FI2_K*R2_K*R2_K
         RMASSIBB=RMASSIBB+FI2R2
      ENDDO
      FI3_K=FI3(K)
      FI4_K=FI4(K)
      FI5_K=FI5(K)
      R3_K=R3(K)
      R4_K=R4(K)
      R5_K=R5(K)
      FI3R3=FI3_K*R3_K*R3_K
      FI4R4=FI4_K*R4_K*R4_K
      FI5R5=FI5_K*R5_K*R5_K
      RMASSIBB=RMASSIBB+FI3R3
      RMASSIBB=RMASSIBB+FI4R4
      RMASSIBB=RMASSIBB+FI5R5
   ENDDO
   RMASSIBB=RMASSIBB*COL3*RORI
   IF(RMASSIBB.LT.0.0) RMASSIBB=0.0
   RMASSLBB=RMASSLBB*COL3*RORI
   IF(RMASSLBB.LT.0.0) RMASSLBB=0.0
! AFTER  JERNEWF
   DO K=1,NKR
      FI1_K=PSI1(K)
      R1_K=R1(K)
      FI1R1=FI1_K*R1_K*R1_K
      RMASSLAA=RMASSLAA+FI1R1
      DO ICE =1,ICEMAX
         FI2(K,ICE)=PSI2(K,ICE)
         FI2_K=FI2(K,ICE)
         R2_K=R2(K,ICE)
         FI2R2=FI2_K*R2_K*R2_K
         RMASSIAA=RMASSIAA+FI2R2
      ENDDO
      FI3_K=PSI3(K)
      FI4_K=PSI4(K)
      FI5_K=PSI5(K)
      R3_K=R3(K)
      R4_K=R4(K)
      R5_K=R5(K)
      FI3R3=FI3_K*R3_K*R3_K
      FI4R4=FI4_K*R4_K*R4_K
      FI5R5=FI5_K*R5_K*R5_K
      RMASSIAA=RMASSIAA+FI3R3
      RMASSIAA=RMASSIAA+FI4R4
      RMASSIAA=RMASSIAA+FI5R5
   ENDDO
   RMASSIAA=RMASSIAA*COL3*RORI
   IF(RMASSIAA.LE.0.0) RMASSIAA=0.0
   RMASSLAA=RMASSLAA*COL3*RORI
   IF(RMASSLAA.LT.0.0) RMASSLAA=0.0
   DELMASSL1=RMASSLAA-RMASSLBB
   DELMASSI1=RMASSIAA-RMASSIBB
   DELTAQ1=DELMASSL1+DELMASSI1
   QPN=QPS-DELTAQ1
   DAL1=AL1
   DAL2=AL2
   TPN=TPS+DAL1*DELMASSL1+DAL2*DELMASSI1
! SUPERSATURATION
   ARGEXP=-BB1_MY/TPN
   ES1N=AA1_MY*DEXP(ARGEXP)
   ARGEXP=-BB2_MY/TPN
   ES2N=AA2_MY*DEXP(ARGEXP)
   EW1N=OPER3(QPN,PP)
   IF(ES1N.EQ.0)THEN
    DEL1N=0.5
    DIV1=1.5
   ELSE
    DIV1=EW1N/ES1N
    DEL1N=EW1N/ES1N-1.
   END IF
   IF(ES2N.EQ.0)THEN
    DEL2N=0.5
    DIV2=1.5
   ELSE
    DEL2N=EW1N/ES2N-1.
    DIV2=EW1N/ES2N
   END IF
IF(TIMENEW.LT.DT) GOTO 16
! -----------------END OF TIME SPLITTING-------------------------------
17 CONTINUE

! Calculate total CCN and IN after cond (CCN/IN regeneration)- J. Fan
totccn_after = 0.0
totin_after = 0.0
do kr = 1, nkr
   totccn_after = totccn_after+psi1(kr)*r1(kr)*col3
   totin_after = totin_after+psi2(kr,1)*r2(kr,1)*col3  &
                            +psi2(kr,2)*r2(kr,2)*col3  &
                            +psi2(kr,3)*r2(kr,3)*col3  &
                            +psi3(kr)*r3(kr)*col3      &
                            +psi4(kr)*r4(kr)*col3      &
                            +psi5(kr)*r5(kr)*col3
enddo
ccnreg = max((totccn_before-totccn_after),0.0)
inreg = max((totin_before-totin_after),0.0)

TT=TPN
QQ=QPN
DO KR=1,NKR
   FF1(KR)=PSI1(KR)
   DO ICE=1,ICEMAX
      FF2(KR,ICE)=PSI2(KR,ICE)
   ENDDO
   FF3(KR)=PSI3(KR)
   FF4(KR)=PSI4(KR)
   FF5(KR)=PSI5(KR)
ENDDO

RETURN      
END SUBROUTINE ONECOND3

!##############################################################################
Subroutine COAL_BOTT_NEW (FF1R,FF2R,FF3R, &
   FF4R,FF5R,TT,QQ,PP,RHO,dthalf,TCRIT,TTCOAL, &
   kdir,idir,jdir)

use mem_micro
use mem_grid, only:ngrid

IMPLICIT NONE

INTEGER kdir,idir,jdir
INTEGER KR,ICE
INTEGER icol_drop,icol_snow,icol_graupel,icol_hail, &
 icol_column,icol_plate,icol_dendrite,icol_drop_brk
double precision  g1(nkr),g2(nkr,icemax),g3(nkr),g4(nkr),g5(nkr)
double precision gdumb(nkr),xl_dumb(0:nkr)
double precision g2_1(nkr),g2_2(nkr),g2_3(nkr)
real cont_fin_drop,dconc,conc_icempl,deldrop,t_new, &
 cont_fin_ice,conc_old,conc_new,cont_init_ice, &
 cont_init_drop,ALWC
REAL    FF1R(NKR),FF2R(NKR,ICEMAX),FF3R(NKR),FF4R(NKR),FF5R(NKR)
REAL DTHALF
REAL TCRIT,TTCOAL
INTEGER I,J,IT,NDIV
REAL RHO
REAL :: totalcloud,totalrain,totalprissnow
DOUBLE PRECISION break_drop_bef,break_drop_aft,dtbreakup
DOUBLE PRECISION break_drop_per
DOUBLE PRECISION TT,QQ,PP,prdkrn,prdkrn1
parameter (prdkrn1=1.d0)

icol_drop_brk=0
icol_drop=0
icol_snow=0
icol_graupel=0
icol_hail=0
icol_column=0
icol_plate=0
icol_dendrite=0
t_new=tt
CALL MISC1 (PP,cwll_1000mb,cwll_750mb,cwll_500mb,cwll,nkr)
! THIS IS FOR BREAKUP
DO I=1,NKR
   DO J=1,NKR
      CWLL(I,J)=ECOALMASSM(I,J)*CWLL(I,J)
   ENDDO
ENDDO

! THIS IS FOR TURBULENCE
IF (LIQTURB.EQ.1)THEN
   DO I=1,KRMAX_LL
      DO J=1,KRMAX_LL
         CWLL(I,J)=CTURBLL(I,J)*CWLL(I,J)
      END DO
   END DO
END IF
CALL modkrn (TT,QQ,PP,PRDKRN,TTCOAL)
DO KR=1,NKR
   G1(KR)=FF1R(KR)*3.*XL(KR)*XL(KR)*1.E3
   G2(KR,1)=FF2R(KR,1)*3*xi(KR,1)*XI(KR,1)*1.e3
   G2(KR,2)=FF2R(KR,2)*3.*xi(KR,2)*XI(KR,2)*1.e3
   G2(KR,3)=FF2R(KR,3)*3.*xi(KR,3)*XI(KR,3)*1.e3
   G3(KR)=FF3R(KR)*3.*xs(kr)*xs(kr)*1.e3
   G4(KR)=FF4R(KR)*3.*xg(kr)*xg(kr)*1.e3
   G5(KR)=FF5R(KR)*3.*xh(kr)*xh(kr)*1.e3
   g2_1(kr)=g2(KR,1)
   g2_2(KR)=g2(KR,2)
   g2_3(KR)=g2(KR,3)
   if(kr.gt.(nkr-jbreak).and.g1(kr).gt.1.e-17)icol_drop_brk=1
   if(g1(kr).gt.1.e-10)icol_drop=1
   if(g2_1(kr).gt.1.e-10)icol_column=1
   if(g2_2(kr).gt.1.e-10)icol_plate=1
   if(g2_3(kr).gt.1.e-10)icol_dendrite=1
   if(g3(kr).gt.1.e-10)icol_snow=1
   if(g4(kr).gt.1.e-10)icol_graupel=1
   if(g5(kr).gt.1.e-10)icol_hail=1
ENDDO 
! calculation of initial hydromteors content in g/cm**3 :
cont_init_drop=0.
cont_init_ice=0.
do kr=1,nkr
   cont_init_drop=cont_init_drop+g1(kr)
   cont_init_ice=cont_init_ice+g3(kr)+g4(kr)+g5(kr)
   do ice=1,icemax
      cont_init_ice=cont_init_ice+g2(kr,ice)
   enddo
enddo
cont_init_drop=col*cont_init_drop*1.e-3
cont_init_ice=col*cont_init_ice*1.e-3
! calculation of alwc in g/m**3
alwc=cont_init_drop*1.e6
! calculation of initial hydromteors content in g/cm**3 :

!----------LIQUID-LIQUID INTERACTIONS-----------------------
totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
if (icol_drop.eq.1)then 
  if (imbudget >= 1) then
     micro_g(ngrid)%cld2raint(kdir,idir,jdir) = micro_g(ngrid)%cld2raint(kdir,idir,jdir) &
         + totalcloud 
  endif
  CALL coll_xxx (G1,CWLL,XL_MG,CHUCM,IMA,NKR)
  if (imbudget >= 1) then
     totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
     micro_g(ngrid)%cld2raint(kdir,idir,jdir) = micro_g(ngrid)%cld2raint(kdir,idir,jdir) &
         - totalcloud 
  endif

!-----------RAINDROP BREAKUP------------------------------
  if(icol_drop_brk.eq.1)then
     ndiv=1
10 continue
     do it = 1,ndiv
       if (ndiv.gt.10000)stop
       dtbreakup = dthalf/ndiv
       if (it.eq.1)then
          do kr=1,nkr
             gdumb(kr)= g1(kr)*1.D-3
             xl_dumb(kr)=xl_mg(KR)*1.D-3
          end do
          break_drop_bef=0.d0
          do kr=1,nkr
             break_drop_bef=break_drop_bef+g1(kr)*1.D-3
          enddo
       end if
       CALL breakup (gdumb,xl_dumb,dtbreakup,brkweight, &
            pkij,qkj,nkr,jbreak)
     end do
     break_drop_aft=0.0d0
     do kr=1,nkr
        break_drop_aft=break_drop_aft+gdumb(kr)
     enddo
     break_drop_per=break_drop_aft/break_drop_bef
     if (break_drop_per.gt.1.001)then
         ndiv=ndiv*2
         GO TO 10
     else
        do kr=1,nkr
           g1(kr)=gdumb(kr)*1.D3
        end do
     end if
  end if
  totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
!-----------------SNOW-LIQUID COLLISIONS---------------------
 if (iceprocs .eq. 1) then
  if (imbudget >= 1) then
     micro_g(ngrid)%rimecldt(kdir,idir,jdir) = micro_g(ngrid)%rimecldt(kdir,idir,jdir) &
          + totalcloud
     micro_g(ngrid)%rain2icet(kdir,idir,jdir) = micro_g(ngrid)%rain2icet(kdir,idir,jdir) &
          + totalrain
  endif
  if (icol_snow.eq.1)then 
     if (imbudget >= 2) then
        micro_g(ngrid)%rimecldaggrt(kdir,idir,jdir) = micro_g(ngrid)%rimecldaggrt(kdir,idir,jdir) &
            + totalcloud 
        micro_g(ngrid)%rain2agt(kdir,idir,jdir) = micro_g(ngrid)%rain2agt(kdir,idir,jdir) &
            + totalrain
     endif
     if(tt.lt.tcrit) then
        CALL coll_xyz (g1,g3,g4,cwls,xl_mg,xs_mg, &
                       chucm,ima,prdkrn1,nkr,0)
     endif
     if(tt.ge.tcrit) then
        CALL coll_xyz (g1,g3,g5,cwls,xl_mg,xs_mg, &
                       chucm,ima,prdkrn1,nkr,0)
     endif
     if(alwc.lt.alcr) then
        CALL coll_xyx (g3,g1,cwsl,xs_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     endif
     if(alwc.ge.alcr) then
        CALL coll_xyz (g3,g1,g4,cwsl,xs_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     endif
     if (imbudget >= 2) then
        totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
        micro_g(ngrid)%rimecldaggrt(kdir,idir,jdir) = micro_g(ngrid)%rimecldaggrt(kdir,idir,jdir) &
            - totalcloud 
        totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
        micro_g(ngrid)%rain2agt(kdir,idir,jdir) = micro_g(ngrid)%rain2agt(kdir,idir,jdir) &
            - totalrain
     endif
  end if
!-----------------GRAUPEL-LIQUID COLLISIONS----------------
! water - graupel = graupel (t < tcrit ; xl_mg ge xg_mg)
! graupel - water = graupel (t < tcrit ; xg_mg > xl_mg)
! water - graupel = hail (t ge tcrit ; xl_mg ge xg_mg)
! graupel - water = hail (t ge tcrit ; xg_mg > xl_mg)
  if (icol_graupel.eq.1)then 
     if (imbudget >= 2) then
        micro_g(ngrid)%rimecldgraut(kdir,idir,jdir) = micro_g(ngrid)%rimecldgraut(kdir,idir,jdir) &
            + totalcloud 
        micro_g(ngrid)%rain2grt(kdir,idir,jdir) = micro_g(ngrid)%rain2grt(kdir,idir,jdir) &
            + totalrain
     endif
     if (alwc.lt.alcr_hail) then
        CALL coll_xyy (g1,g4,cwlg,xl_mg,xg_mg, &
                       chucm,ima,prdkrn1,nkr,0)
! if cwc .gt. alcr_hail, water+graupel - hail for nr .gt. kp_hail
     else
        CALL coll_xyyz (g1,g4,g5,cwlg,xl_mg,xg_mg,chucm, &
                       ima,prdkrn1,nkr,0,kp_hail)
     endif
 ! ice-multiplication preparation:
     conc_old=0.
     conc_new=0.
     do kr=kr_icempl,nkr
        conc_old=conc_old+col*g1(kr)/xl_mg(kr)
     enddo
     if (alwc.lt.alcr_hail) then
        CALL coll_xyx (g4,g1,cwgl,xg_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
     else
        CALL coll_xyxz (g4,g1,g5,cwlg,xg_mg,xl_mg,chucm, &
                       ima,prdkrn1,nkr,1,kp_hail)
     endif
     if (imbudget >= 2) then
         totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
         micro_g(ngrid)%rimecldgraut(kdir,idir,jdir) = micro_g(ngrid)%rimecldgraut(kdir,idir,jdir) &
             - totalcloud 
         totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
         micro_g(ngrid)%rain2grt(kdir,idir,jdir) = micro_g(ngrid)%rain2grt(kdir,idir,jdir) &
             - totalrain
     endif
!----ICE-MULTIPLICATION : Hallet-Mossop processes (1 per 250 collisions)
     if (icempl.eq.1) then
        if (tt.ge.265.15.and.tt.le.270.15) then
           do kr=kr_icempl,nkr
              conc_new=conc_new+col*g1(kr)/xl_mg(kr)
           enddo
           dconc=conc_old-conc_new
           if(tt.le.268.15) then
              conc_icempl=dconc*4.e-3*(265.15-tt)/(265.15-268.15)
           elseif(tt.gt.268.15) then
              conc_icempl=dconc*4.e-3*(tcrit-tt)/(tcrit-268.15)
           endif
           g2_2(1)=g2_2(1)+conc_icempl*xi2_mg(1)/col
         endif
     endif
  endif
!---------------HAIL-LIQUID COLLISIONS-----------------------------
! water - hail = hail (xl_mg ge xh_mg)    (kxyy=2)
! hail - water = hail (xh_mg > xl_mg)     (kxyx=3)
  if(icol_hail.eq.1) then
    if (imbudget >= 2) then
       micro_g(ngrid)%rimecldhailt(kdir,idir,jdir) = micro_g(ngrid)%rimecldhailt(kdir,idir,jdir) &
           + totalcloud 
       micro_g(ngrid)%rain2hat(kdir,idir,jdir) = micro_g(ngrid)%rain2hat(kdir,idir,jdir) &
           + totalrain
    endif
    CALL coll_xyy (g1,g5,cwlh,xl_mg,xh_mg, &
                       chucm,ima,prdkrn1,nkr,0)
    CALL coll_xyx (g5,g1,cwhl,xh_mg,xl_mg, &
                       chucm,ima,prdkrn1,nkr,1)
    if (imbudget >= 2) then
       totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
       micro_g(ngrid)%rimecldhailt(kdir,idir,jdir) = micro_g(ngrid)%rimecldhailt(kdir,idir,jdir) &
           - totalcloud 
       totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
       micro_g(ngrid)%rain2hat(kdir,idir,jdir) = micro_g(ngrid)%rain2hat(kdir,idir,jdir) &
           - totalrain
    endif
  endif
!--------------COLUMN-LIQUID COLLISIONS--------------------------
! interactions between water and crystals :
! interactions between water and columns :
! water - columns = graupel (t < tcrit ; xl_mg ge xi_mg)    (kxyz=6)
! water - columns = hail (t ge tcrit ; xl_mg ge xi_mg)      (kxyz=7)
! columns - water = columns/graupel (xi_mg > xl_mg)             (kxyx=4); kxyxz=2)
! now: columns - water = columns (xi_mg > xl_mg)             (kxyx=4); kxyxz=2)
  if (imbudget >= 2 .and. ipris > 0) then
    micro_g(ngrid)%rimecldsnowt(kdir,idir,jdir) = micro_g(ngrid)%rimecldsnowt(kdir,idir,jdir) &
        + totalcloud
    micro_g(ngrid)%rain2snt(kdir,idir,jdir) = micro_g(ngrid)%rain2snt(kdir,idir,jdir) &
        + totalrain
  endif
  if(icol_column.eq.1) then
     if(tt.lt.tcrit) then
        CALL coll_xyz (g1,g2_1,g4,cwli_1,xl_mg,xi1_mg, &
                         chucm,ima,prdkrn,nkr,0)
     endif
     if(tt.ge.tcrit) then
        CALL coll_xyz (g1,g2_1,g5,cwli_1,xl_mg,xi1_mg, &
                        chucm,ima,prdkrn,nkr,0)
     endif
     CALL coll_xyx (g2_1,g1,cwil_1,xi1_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
!------------PLATE-LIQUID COLLISIONS---------------------------------
! water - plates = graupel (t < tcrit ; xl_mg ge xi2_mg)    (kxyz=8)
! water - plates = hail (t ge tcrit ; xl_mg ge xi2_mg)      (kxyz=9)
! plates - water = plates/graupel (xi2_mg > xl_mg)              (kxyx=5; kxyxz=3)
!now: plates - water = plates (xi2_mg > xl_mg)              (kxyx=5; kxyxz=3)
  if(icol_plate.eq.1) then
    if(tt.lt.tcrit) then
      CALL coll_xyz (g1,g2_2,g4,cwli_2,xl_mg,xi2_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    if(tt.ge.tcrit) then
      CALL coll_xyz (g1,g2_2,g5,cwli_2,xl_mg,xi2_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    CALL coll_xyx (g2_2,g1,cwil_2,xi2_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
!--------------DENDRITE-LIQUID COLLISIONS------------------------------
! water - dendrites = graupel (t < tcrit ; xl_mg ge xi3_mg) (kxyz=10)
! water - dendrites = hail (t ge tcrit ; xl_mg ge xi3_mg)   (kxyz=11)
! dendrites - water = dendrites/graupel (xi3_mg > xl_mg)         (kxyx=6; kxyxz=4)
!now dendrites - water = dendrites (xi3_mg > xl_mg)         (kxyx=6; kxyxz=4)
  if(icol_dendrite.eq.1) then
    if(tt.lt.tcrit) then
      CALL coll_xyz (g1,g2_3,g4,cwli_3,xl_mg,xi3_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    if(tt.ge.tcrit) then
      CALL coll_xyz (g1,g2_3,g5,cwli_3,xl_mg,xi3_mg, &
                         chucm,ima,prdkrn,nkr,0)
    endif
    CALL coll_xyx (g2_3,g1,cwil_3,xi3_mg,xl_mg, &
                         chucm,ima,prdkrn,nkr,1)
  endif
  if (imbudget >= 1) then
     totalcloud = sum(g1(1:krdrop-1))*col*1e-3/rho
     totalrain = sum(g1(krdrop:nkr))*col*1e-3/rho
     micro_g(ngrid)%rimecldt(kdir,idir,jdir) = micro_g(ngrid)%rimecldt(kdir,idir,jdir) &
         - totalcloud
     micro_g(ngrid)%rain2icet(kdir,idir,jdir) = micro_g(ngrid)%rain2icet(kdir,idir,jdir) &
         - totalrain
  endif
  if (imbudget >= 2 .and. ipris > 0) then
     micro_g(ngrid)%rimecldsnowt(kdir,idir,jdir) = micro_g(ngrid)%rimecldsnowt(kdir,idir,jdir) &
         - totalcloud
     micro_g(ngrid)%rain2snt(kdir,idir,jdir) = micro_g(ngrid)%rain2snt(kdir,idir,jdir) &
         - totalrain
  endif
 endif !endif ice present
endif !endif liquid collisions

!------------CRYSTAL-CRYSTAL INTERACTIONS--------------------------
! if(t.le.TTCOAL) - no interactions between crystals
if(tt.gt.TTCOAL .and. iceprocs == 1) then
  if (imbudget >= 1) then
     totalprissnow = sum(g3)*col*1e-3/rho
     micro_g(ngrid)%aggregatet(kdir,idir,jdir) = micro_g(ngrid)%aggregatet(kdir,idir,jdir) &
         - totalprissnow
  endif
!------------COLUMN-ICE COLLISIONS--------------------------------
if(icol_column.eq.1) then
  ! columns - columns = snow
  CALL coll_xxy (g2_1,g3,cwii_1_1,xi1_mg, &
     chucm,ima,prdkrn,nkr)
  ! columns - plates = snow (xi1_mg ge xi2_mg)                (kxyz=12)
  ! plates - columns = snow (xi2_mg > xi1_mg)                 (kxyz=13)
  if(icol_plate.eq.1) then     
   CALL coll_xyz (g2_1,g2_2,g3,cwii_1_2,xi1_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,0)
   CALL coll_xyz (g2_2,g2_1,g3,cwii_2_1,xi2_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
  ! columns - dendrites = snow (xi1_mg ge xi3_mg)             (kxyz=14)
  ! dendrites - columns = snow (xi3_mg > xi1_mg)              (kxyz=15)
  if(icol_dendrite.eq.1) then
     CALL coll_xyz (g2_1,g2_3,g3,cwii_1_3,xi1_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,0)
     CALL coll_xyz (g2_3,g2_1,g3,cwii_3_1,xi3_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
  ! columns - snow = snow (xi1_mg ge xs_mg)                   (kxyy=3)
  ! snow - columns = snow (xs_mg > xi1_mg)                    (kxyx=7)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_1,g3,cwis_1,xi1_mg,xs_mg, &
!                      chucm,ima,prdkrn,nkr,0)
   CALL coll_xyx (g3,g2_1,cwsi_1,xs_mg,xi1_mg, &
     chucm,ima,prdkrn,nkr,1)
  endif          
endif
!-----------------PLATES-ICE COLLISIONS-------------------------
 ! plates - plates = snow
if(icol_plate.eq.1) then
   CALL coll_xxy (g2_2,g3,cwii_2_2,xi2_mg, &
     chucm,ima,prdkrn,nkr)
  ! plates - dendrites = snow (xi2_mg ge xi3_mg)              (kxyz=17)
  ! dendrites - plates = snow (xi3_mg > xi2_mg)               (kxyz=18)
  if(icol_dendrite.eq.1) then
   CALL coll_xyz (g2_2,g2_3,g3,cwii_2_3,xi2_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,0)
   CALL coll_xyz (g2_3,g2_2,g3,cwii_3_2,xi3_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
 ! plates - snow = snow (xi2_mg ge xs_mg)                    (kxyy=4)
 ! snow - plates = snow (xs_mg > xi2_mg)                     (kxyx=12)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_2,g3,cwis_2,xi2_mg,xs_mg, &
!                      chucm,ima,prdkrn,nkr,0)
   CALL coll_xyx (g3,g2_2,cwsi_2,xs_mg,xi2_mg, &
     chucm,ima,prdkrn,nkr,1)
   end if
endif
!------------DENDRITES-ICE COLLISIONS---------------------------
 ! dendrites - dendrites = snow
if(icol_dendrite.eq.1) then
  CALL coll_xxy (g2_3,g3,cwii_3_3,xi3_mg, &
                 chucm,ima,prdkrn,nkr)
  ! dendrites - snow = snow (xi3_mg ge xs_mg)                 (kxyy=5)
  ! snow - dendrites = snow (xs_mg > xi3_mg)                  (kxyx=17)
  if(icol_snow.eq.1) then
!       call coll_xyy (g2_3,g3,cwis_3,xi3_mg,xs_mg, &
!                    chucm,ima,prdkrn,nkr,0)
   CALL coll_xyx (g3,g2_3,cwsi_3,xs_mg,xi3_mg, &
     chucm,ima,prdkrn,nkr,1)
  end if
endif
!-------------SNOW-SNOW COLLISIONS-------------------------------
  if(icol_snow.ne.0) then
   ! snow - snow = snow
   CALL coll_xxx_prd (g3,cwss,xs_mg,chucm,ima,prdkrn,nkr)
   ! snow - graupel = snow (xs_mg > xg_mg)                     (kxyx=22)
   ! graupel - snow = graupel (xg_mg ge xs_mg)                 (kxyx=23)

   !JF4 no snow-graupel interactions
   !         if(icol_graupel.eq.1) then
   !          call coll_xyx (g3,g4,cwsg,xs_mg,xg_mg, &
   !     &                chucm,ima,prdkrn,nkr,1)
   !         endif
   !end JF4
  endif
  if (imbudget >= 1) then
     totalprissnow = sum(g3)*col*1e-3/rho
     micro_g(ngrid)%aggregatet(kdir,idir,jdir) = micro_g(ngrid)%aggregatet(kdir,idir,jdir) &
         + totalprissnow
  endif
endif

!-------------END OF COLLISION CALCULATIONS-------------------

! calculation of finish hydrometeors contents in g/cm**3 :
cont_fin_drop=0.
cont_fin_ice=0.
do kr=1,nkr
   g2(kr,1)=g2_1(kr)
   g2(kr,2)=g2_2(kr)
   g2(kr,3)=g2_3(kr)
   cont_fin_drop=cont_fin_drop+g1(kr)
   cont_fin_ice=cont_fin_ice+g3(kr)+g4(kr)+g5(kr)
   do ice=1,icemax
      cont_fin_ice=cont_fin_ice+g2(kr,ice)
   enddo
enddo
cont_fin_drop=col*cont_fin_drop*1.e-3
cont_fin_ice=col*cont_fin_ice*1.e-3
deldrop=cont_init_drop-cont_fin_drop
! deldrop in g/cm**3
! resulted value of temperature (rob in g/cm**3) :
if(t_new.le.273.15) then
  if(deldrop.ge.0.) then
    t_new=t_new+320.*deldrop/rho
  else
! if deldrop < 0
    if(abs(deldrop).gt.cont_init_drop*0.05) then
      print *, '*** deldrop < 0 in coal_bott ***'
      print *,'*** coal_bott ***'
      print*,'cont_fin_drop = ',cont_fin_drop,cont_fin_ice
      print*,'cont_init_drop = ',cont_init_drop,cont_init_ice
      stop
    endif
  endif
endif

! recalculation of density function f1,f2,f3,f4,f5 in 1/(g*cm**3) :  
DO KR=1,NKR
   FF1R(KR)=G1(KR)/(3.*XL(KR)*XL(KR)*1.E3)
   FF2R(KR,1)=G2(KR,1)/(3*xi(KR,1)*XI(KR,1)*1.e3)
   FF2R(KR,2)=G2(KR,2)/(3.*xi(KR,2)*XI(KR,2)*1.e3)
   FF2R(KR,3)=G2(KR,3)/(3.*xi(KR,3)*XI(KR,3)*1.e3)
   FF3R(KR)=G3(KR)/(3.*xs(kr)*xs(kr)*1.e3)
   FF4R(KR)=G4(KR)/(3.*xg(kr)*xg(kr)*1.e3)
   FF5R(KR)=G5(KR)/(3.*xh(kr)*xh(kr)*1.e3)
ENDDO

tt=t_new

RETURN
END SUBROUTINE COAL_BOTT_NEW

!##############################################################################
Subroutine MISC1 (PP,cwll_1000mb,cwll_750mb,cwll_500mb,cwll,nkr)

IMPLICIT NONE
INTEGER kr1,kr2,NKR
DOUBLE PRECISION PP
REAL P_Z
double precision cwll(nkr,nkr),cwll_1,cwll_2,cwll_3 &
     ,cwll_1000mb(nkr,nkr),cwll_750mb(nkr,nkr),cwll_500mb(nkr,nkr)

P_Z=PP
do kr1=1,nkr
   do kr2=1,nkr
      cwll_1=cwll_1000mb(kr1,kr2)
      cwll_2=cwll_750mb(kr1,kr2)
      cwll_3=cwll_500mb(kr1,kr2)
      if(p_z.ge.p1) cwll(kr1,kr2)=cwll_1
      if(p_z.eq.p2) cwll(kr1,kr2)=cwll_2
      if(p_z.eq.p3) cwll(kr1,kr2)=cwll_3
      if(p_z.lt.p1.and.p_z.gt.p2) &
         cwll(kr1,kr2)=cwll_2+ &
         (cwll_1-cwll_2)*(p_z-p2)/(p1-p2) 
      if(p_z.lt.p2.and.p_z.gt.p3) &
         cwll(kr1,kr2)=cwll_3+ &
         (cwll_2-cwll_3)*(p_z-p3)/(p2-p3)
      if(p_z.lt.p3) cwll(kr1,kr2)=cwll_3
   enddo
enddo

RETURN
END SUBROUTINE MISC1

!##############################################################################
Subroutine coll_xxx (g,ckxx,x,chucm,ima,nkr)

implicit double precision (a-h,o-z)
dimension g(nkr),ckxx(nkr,nkr),x(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(g(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(g(i).gt.gmin) exit
enddo

do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*g(i)*g(j)
      x0=min(x0,g(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,g(j)*x(i))
      endif
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      g(i)=g(i)-gsi
      if(g(i).lt.0.d0) g(i)=0.d0
      g(j)=g(j)-gsj
      gk=g(k)+gsk
      if(g(j).lt.0.d0.and.gk.lt.gmin) then
         g(j)=0.d0
         g(k)=g(k)+gsi
      endif
      flux=0.d0

      if(gk.gt.gmin) then
         x1=dlog(g(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if

         g(k)=gk-flux
         if(gk .lt. flux) flux=gk
         if(g(k).lt.0.d0) g(k)=0.d0
         g(kp)=g(kp)+flux
      endif
   end do
end do

return
END SUBROUTINE coll_xxx

!##############################################################################
Subroutine coll_xxx_prd (g,ckxx,x,chucm,ima,prdkrn,nkr)

implicit double precision (a-h,o-z)
dimension g(nkr),ckxx(nkr,nkr),x(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(g(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(g(i).gt.gmin) exit
enddo

do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*g(i)*g(j)*prdkrn
      x0=min(x0,g(i)*x(j))
      if(j.ne.k) then
        x0=min(x0,g(j)*x(i))
      endif
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      g(i)=g(i)-gsi
      if(g(i).lt.0.d0) g(i)=0.d0
      g(j)=g(j)-gsj
      gk=g(k)+gsk
      if(g(j).lt.0.d0.and.gk.lt.gmin) then
         g(j)=0.d0
         g(k)=g(k)+gsi
      endif
      flux=0.d0

      if(gk.gt.gmin) then
         x1=dlog(g(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         g(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(g(k).lt.0.d0) g(k)=0.d0
          g(kp)=g(kp)+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xxx_prd 

!##############################################################################
Subroutine modkrn (TT,QQ,PP,PRDKRN,TTCOAL)

implicit none
real epsf,tc,ttt1,ttt,factor,qs2,qq1,dele,f,factor_t
double precision TT,QQ,PP,satq2,t,p
double precision prdkrn
REAL at,bt,ct,dt,temp,a,b,c,d,tc_min,tc_max
real factor_max,factor_min
REAL TTCOAL
data at,bt,ct,dt/0.88333,0.0931878,0.0034793,4.5185186e-05/

satq2(t,p)=3.80e3*(10**(9.76421-2667.1/t))/p
temp(a,b,c,d,tc)=d*tc*tc*tc+c*tc*tc+b*tc+a

IF (QQ.LE.0)QQ=1.0E-12
epsf=.5
tc=tt-273.15
if(tc.le.0) then
  ttt1  =temp(at,bt,ct,dt,tc)
  ttt   =ttt1
  qs2   =satq2(tt,pp)
  qq1   =qq*(0.622+0.378*qs2)/(0.622+0.378*qq)/qs2
  dele  =ttt*qq1
  if(tc.ge.-6.) then
    factor = dele
    if(factor.lt.epsf) factor=epsf
    if(factor.gt.1.) factor=1.
  endif                        
  factor_t=factor
  if(tc.ge.-12.5.and.tc.lt.-6.) factor_t=0.5
  if(tc.ge.-17.0.and.tc.lt.-12.5) factor_t=1.
  if(tc.ge.-20.0.and.tc.lt.-17.) factor_t=0.4
  if(tc.lt.-20.) then
    tc_min=ttcoal-273.15
    tc_max=-20.
    factor_max=0.25
    factor_min=0.
    f=factor_min+(tc-tc_min)*(factor_max-factor_min)/  &
                                    (tc_max-tc_min)
    factor_t=f
  endif
  if (factor_t.lt.0)factor_t=0.01
  prdkrn=factor_t
else
  prdkrn=1.d0
end if

RETURN
END SUBROUTINE modkrn 

!##############################################################################
Subroutine coll_xxy (gx,gy,ckxx,x,chucm,ima,prdkrn,nkr)

implicit double precision (a-h,o-z)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
dimension  gx(nkr),gy(nkr),ckxx(nkr,nkr),x(0:nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! collisions
do i=ix0,ix1
   do j=i,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxx(i,j)*gx(i)*gx(j)*prdkrn
      x0=min(x0,gx(i)*x(j))
      x0=min(x0,gx(j)*x(i))
      gsi=x0/x(j)
      gsj=x0/x(i)
      gsk=gsi+gsj
      gx(i)=gx(i)-gsi
      if(gx(i).lt.0.d0) gx(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
         x1=dlog(gy(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gy(k)=gk-flux
         if(gk .lt. flux) flux=gk
         if(gy(k).lt.0.d0) gy(k)=0.d0
         gy(kp)=gy(kp)+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xxy

!##############################################################################
Subroutine coll_xyy (gx,gy,ckxy,x,y,chucm,ima, &
  prdkrn,nkr,indc)

implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1         
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      x0=min(x0,gx(j)*y(i))
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
        x1=dlog(gy(kp)/gk+1.d-15)
        if (x1.eq.0)then
           flux=0  
        else
           flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
           flux=min(flux,gsk)
        end if
        gy(k)=gk-flux
        if(gk .lt. flux) flux=gk
        if(gy(k).lt.0.d0) gy(k)=0.d0
        gy(kp)=gy(kp)+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xyy

!##############################################################################
Subroutine coll_xyx (gx,gy,ckxy,x,y,chucm,ima,prdkrn,nkr,indc)

implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,gx(j)*y(i))
      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      gk=gx(k)+gsk
      if(gx(j).lt.0.d0.and.gk.lt.gmin) then
         gx(j)=0.d0
         gx(k)=gx(k)+gsi
      endif
      flux=0.d0            
      if(gk.gt.gmin) then
         x1=dlog(gx(kp)/gk+1.d-15)
         if (x1.eq.0)then
             flux=0  
         else
             flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
             flux=min(flux,gsk)
         end if
         gx(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(gx(k).lt.0.d0) gx(k)=0.d0
         gx(kp)=gx(kp)+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xyx

!##############################################################################
Subroutine coll_xyxz (gx,gy,gz,ckxy,x,y,chucm,ima, &
 prdkrn,nkr,indc, kp_bound)

implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
integer kp_bound

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo

! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      if(j.ne.k) then
         x0=min(x0,gx(j)*y(i))
      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      gk=gx(k)+gsk
      if(gx(j).lt.0.d0.and.gk.lt.gmin) then
        gx(j)=0.d0
        gx(k)=gx(k)+gsi
      endif
      flux=0.d0
      if(kp.lt.kp_bound) gkp=gx(kp)
      if(kp.ge.kp_bound) gkp=gz(kp)
      if(gk.gt.gmin) then
         x1=dlog(gkp/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gx(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(gx(k).lt.0.d0) gx(k)=0.d0
         if(kp.lt.kp_bound) gx(kp)=gkp+flux
         if(kp.ge.kp_bound) gz(kp)=gkp+flux

      endif
   enddo
enddo

return
END SUBROUTINE coll_xyxz

!##############################################################################
Subroutine coll_xyz (gx,gy,gz,ckxy,x,y,chucm,ima, &
                    prdkrn,nkr,indc)

implicit double precision (a-h,o-z)
dimension gx(nkr),gy(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo
! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
    do j=jmin+indc,ix1         
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      x0=min(x0,gx(j)*y(i))
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gz(k)+gsk
      flux=0.d0
      if(gk.gt.gmin) then
         x1=dlog(gz(kp)/gk+1.d-15)
         if (x1.eq.0)then
            flux=0  
         else
            flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
            flux=min(flux,gsk)
         end if
         gz(k)=gk-flux
         if(gk .lt. flux) flux=gk

         if(gz(k).lt.0.d0) gz(k)=0.d0
         gz(kp)=gz(kp)+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xyz

!##############################################################################
Subroutine coll_xyyz (gx,gy,gz,ckxy,x,y,chucm,ima,prdkrn,nkr,indc,kp_bound)

implicit double precision (a-h,o-z)
dimension gy(nkr),gx(nkr),gz(nkr),ckxy(nkr,nkr),x(0:nkr),y(0:nkr)
dimension chucm(nkr,nkr)
double precision ima(nkr,nkr)
integer kp_bound

gmin=1.d-60
! lower and upper integration limit ix0,ix1
do i=1,nkr-1
   ix0=i
   if(gx(i).gt.gmin) exit
enddo

if(ix0.eq.nkr-1) return

do i=nkr-1,1,-1
   ix1=i
   if(gx(i).gt.gmin) exit
enddo
! lower and upper integration limit iy0,iy1
do i=1,nkr-1
   iy0=i
   if(gy(i).gt.gmin) exit
enddo

if(iy0.eq.nkr-1) return

do i=nkr-1,1,-1
   iy1=i
   if(gy(i).gt.gmin) exit
enddo

! collisions :
do i=iy0,iy1
   jmin=i
   if(jmin.eq.(nkr-1)) return
   if(i.lt.ix0) jmin=ix0-indc
   do j=jmin+indc,ix1
      k=ima(i,j)
      kp=k+1
      x0=ckxy(j,i)*gy(i)*gx(j)*prdkrn
      x0=min(x0,gy(i)*x(j))
      x0=min(x0,gx(j)*y(i))
!      if(j.ne.k) then
!        x0=min(x0,gx(j)*y(i))
!      endif
      gsi=x0/x(j)
      gsj=x0/y(i)
      gsk=gsi+gsj
      gy(i)=gy(i)-gsi
      if(gy(i).lt.0.d0) gy(i)=0.d0
      gx(j)=gx(j)-gsj
      if(gx(j).lt.0.d0) gx(j)=0.d0
      gk=gy(k)+gsk
      flux=0.d0
      if(kp.lt.kp_bound) gkp=gy(kp)
      if(kp.ge.kp_bound) gkp=gz(kp)
      if(gk.gt.gmin) then
        x1=dlog(gkp/gk+1.d-15)
        if (x1.eq.0)then
           flux=0  
        else
           flux=gsk/x1*(dexp(0.5*x1)-dexp(x1*(0.5-chucm(i,j))))
           flux=min(flux,gsk)
        end if
        gy(k)=gk-flux
        if(gk .lt. flux) flux=gk

        if(gy(k).lt.0.d0) gy(k)=0.d0
        if(kp.lt.kp_bound) gy(kp)=gkp+flux
        if(kp.ge.kp_bound) gz(kp)=gkp+flux
      endif
   enddo
enddo

return
END SUBROUTINE coll_xyyz

!##############################################################################
Subroutine breakup (GT_MG,XT_MG,DT,BRKWEIGHT,PKIJ,QKJ,JMAX,JBREAK)

!.....INPUT VARIABLES
!
! GT    : MASS DISTRIBUTION function
! XT_MG : MASS OF BIN IN MG
! JMAX  : NUMBER OF BINS
! DT    : TIMESTEP IN S

INTEGER JMAX

!.....LOCAL VARIABLES

INTEGER JBREAK,AP,IA,JA,KA,IE,JE,KE
DOUBLE PRECISION EPS,NEGSUM

PARAMETER (AP = 1)
PARAMETER (IA = 1)
PARAMETER (JA = 1)
PARAMETER (KA = 1)
PARAMETER (EPS = 1.D-20)

INTEGER I,J,K,JDIFF
DOUBLE PRECISION GT_MG(JMAX),XT_MG(0:JMAX),DT
DOUBLE PRECISION BRKWEIGHT(JBREAK),PKIJ(JBREAK,JBREAK,JBREAK), &
                 QKJ(JBREAK,JBREAK)
DOUBLE PRECISION HLP(JMAX)
DOUBLE PRECISION FT(JMAX),FA(JMAX)
DOUBLE PRECISION DG(JMAX),DF(JMAX),DBREAK(JBREAK),GAIN,LOSS
REAL PI
PARAMETER (PI = 3.1415927)
IE = JBREAK
JE = JBREAK
KE = JBREAK
!.....IN CGS

!.....SHIFT BETWEEN COAGULATION AND BREAKUP GRID

JDIFF = JMAX - JBREAK

!.....INITIALIZATION

!.....TRANSFORMATION FROM G(LN X) = X**2 F(X) TO F(X)
DO J=1,JMAX
   FT(J) = GT_MG(J) / XT_MG(J)**2
ENDDO

!.....SHIFT TO BREAKUP GRID

DO K=1,KE
   FA(K) = FT(K+JDIFF)
ENDDO

!.....BREAKUP: BLECK'S FIRST ORDER METHOD
!
!     PKIJ: GAIN COEFFICIENTS
!     QKJ : LOSS COEFFICIENTS
!

DO K=1,KE
   GAIN = 0.0
   DO I=1,IE
      DO J=1,I
         GAIN = GAIN + FA(I)*FA(J)*PKIJ(K,I,J)
      ENDDO
   ENDDO
   LOSS = 0.0
   DO J=1,JE
      LOSS = LOSS + FA(J)*QKJ(K,J)
   ENDDO
   DBREAK(K) = BRKWEIGHT(K) * (GAIN - FA(K)*LOSS)
ENDDO

!.....SHIFT RATE TO COAGULATION GRID

DO J=1,JDIFF
   DF(J) = 0.0
ENDDO
DO J=1,KE
   DF(J+JDIFF) = DBREAK(J)
ENDDO
!.....TRANSFORMATION TO MASS DISTRIBUTION function G(LN X)

DO J=1,JMAX
   DG(J) = DF(J) * XT_MG(J)**2
ENDDO

!.....TIME INTEGRATION

DO J=1,JMAX
HLP(J) = 0.0
NEGSUM = 0.0
GT_MG(J) = GT_MG(J) + DG(J) * DT
IF (GT_MG(J).LT.0) THEN
   HLP(J) = MIN(GT_MG(J),HLP(J))
   GT_MG(J) = EPS
ENDIF
ENDDO

RETURN
END SUBROUTINE breakup

!##############################################################################
Subroutine ccn_reg (fccn0,fccnin,fccnout,nkr,ccnreg, &
                    ff1r,xl,ff2r,xi,ff3r,xs,ff4r,xg)

!!! added for ccn regeneration (J. Fan Oct 2007)

implicit none

real delccn0,delccn,dn0,total,totalcn,totalhydro,coeff
real ccnreg
integer  nkr, kr
real fccn0(nkr), fccnin(nkr)
real ff1r(nkr), ff2r(nkr,3),ff3r(nkr),ff4r(nkr)
real xl(nkr) ,xi(nkr,3),xs(nkr),xg(nkr)
real ccn0(nkr), regcn0(nkr), regcn(nkr)

real fccnout(nkr)

!-----------------------------------
!       total number of ccn in the "ideal" spectrum
!---------------------------------------------------
dn0=0.
do kr=1,nkr
   dn0=dn0+fccn0(kr)*col
end do
!-----------------------------------
!       total ccn in the current spectrum if all drops evaporated
!----------------------------------
totalcn=0.
do kr=1,nkr
   totalcn=totalcn+fccnin(kr)*col
end do
totalhydro=0.
do kr=1,nkr
   totalhydro=totalhydro+ff1r(kr)*col3*xl(kr)           &
                        +ff2r(kr,1)*col3*xi(kr,1)       &
                        +ff2r(kr,2)*col3*xi(kr,2)       &
                        +ff2r(kr,3)*col3*xi(kr,3)       &
                        +ff3r(kr)*col3*xs(kr)           &
                        +ff4r(kr)*col3*xg(kr)
end do
totalhydro = totalhydro + ccnreg
total=totalhydro+totalcn
coeff=total/dn0
!------------------------------------------
!        spectrum of ccn if all drops evaporated
!---------------------------------------------
delccn0=0.
do kr=1,nkr
   delccn0=delccn0+amax1(0.,(fccn0(kr)*col*coeff-fccnin(kr)*col))
end do
if(abs(delccn0).gt.1.e-5) then
   do kr=1,nkr
      regcn0(kr)=amax1(0.,(fccn0(kr)*col*coeff-fccnin(kr)*col))/delccn0
      ccn0(kr)=fccnin(kr)+totalhydro*regcn0(kr)/col
   end do
else
   do kr=1,nkr
      ccn0(kr)=fccnin(kr)
   end do
end if  ! abs(delccn)...

if (ccnreg > 0.) then
   delccn=0.
   do kr=1,nkr
      delccn=delccn+amax1(0.,(ccn0(kr)*col-fccnin(kr)*col))
   end do
   if(abs(delccn).gt.1.e-5) then
      do kr=1,nkr
         regcn(kr)=amax1(0.,(ccn0(kr)*col-fccnin(kr)*col))/delccn
         fccnout(kr)=fccnin(kr)+ccnreg*regcn(kr)/col
      end do
   end if    ! abs(delccn)...
end if  ! cdreg.gt...

return
END SUBROUTINE ccn_reg

!##############################################################################
Subroutine evapfrz (ccnreg,frzfract,NKR,ICEMAX, TT, DT, xi,ff2r,rndrop)
 
! for dropet freezing

IMPLICIT NONE
DOUBLE PRECISION TT
integer nkr, icemax
REAL,DIMENSION (nkr,icemax) ::FF2R, xi
real rndrop
real ccnreg
real dt
real TPC, dx
integer ITYPE, ice
real frzfract, nfreez
 
! Type of ice
TPC=TT-273.15
ITYPE=0                                                                                  
IF((TPC.GT.-4.0).OR.(TPC.LE.-8.1.AND.TPC.GT.-12.7).OR.&
   (TPC.LE.-17.8.AND.TPC.GT.-22.4)) THEN
  ITYPE=2
ELSE
  IF((TPC.LE.-4.0.AND.TPC.GT.-8.1).OR.(TPC.LE.-22.4)) THEN
    ITYPE=1
  ELSE
    ITYPE=3
  ENDIF
ENDIF
 
ice=ITYPE
 
nfreez = frzfract*ccnreg        !cm-3
! put it in the smallest bin of ice crystal
DX=3.*XI(1,ICE)*COL
ff2r(1,ICE) = ff2r(1,ICE) + nfreez/DX
rndrop= nfreez/DT             !cm-3 s-1
 
return
END SUBROUTINE evapfrz

END MODULE module_hujisbm
