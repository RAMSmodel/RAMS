!##############################################################################
Module micro_prm

use grid_dims
implicit none

! variable to contain path to HUCM input files
character(len=strl1) :: hucmfile

! number of bins 
integer, parameter :: nkr=33

! for bulk nucleation
integer, parameter :: BULKNUC=1

! COL=Ln2/3 (for xl(k+1)=xl(k)*2.)
real :: COL = 0.231049060186648 !=log(2.)/3.
real :: COL3 = 0.693147180559945  !=log(2.)

! Aerosols
!  add for diagnostic CCN. If diagCCN=.false., CCN is prognostic
LOGICAL, PARAMETER :: diagCCN=.false.

! Parameters are used for calculation of 

! RO_SOLUTE - density of aerosol+h20
! ROCCN0, g/cm^3 - density of aerosol
! Molecular weight of dry aerosol
! Number of ions
real :: RO_SOLUTE,ROCCN0,MWAERO,IONS

! ACCN in 1/cm^3, BCCN - coefficients for calculation of FCCNR(KR)
! FCCNR(KR)=1.5*ACCN*BCCN*S_KR**BCCN,
! where S_KR - water supersaturation, %
real, parameter :: ACCN=1.0000E02, BCCN=0.4620E00

! ICCN - set maximum mass of aerosol
! XL(ICCN) - maximum mass of aerosol, is determined as minimum
! mass drop in bin with number ICCN
integer, parameter :: ICCN=1

! if ICEPROCS=1 it is ice microphysics
! ICEPROCS is set to one by setting any of the ice flags in RAMSIN to a value above 0
integer :: ICEPROCS,IFASTSBM

! ICEFLAG=0 USE MEYERS ICE NUCLEATION SCHEME; ICEFLAG=1 USE Classical Theory
integer :: ICEFLAG 

! ICEMAX - number of ice crystal types
integer, parameter :: ICEMAX=3

! flags for turbulance in collection kernals
integer, parameter :: LIQTURB=1
integer, parameter :: ICETURB=1
integer, parameter :: NHYDR=5,NHYDRO=7    &
 ,K0_LL=8,KRMIN_LL=1,KRMAX_LL=19,L0_LL=6  &
 ,IEPS_400=1,IEPS_800=0,IEPS_1600=0  & !For turbulent enhancement to collection kernal
 ,K0L_GL=16,K0G_GL=16                & !Bins that have turbulent enhancement
 ,KRMINL_GL=1,KRMAXL_GL=24           &
 ,KRMING_GL=1,KRMAXG_GL=33           &
 ,KRDROP=14                          & !First Bin number that is rain in FFCD
 ,JMAX=33,JBREAK = 18                  !For rain breakup routines 
integer,dimension(ICEMAX), parameter::    &
  KRPRIS=(/14,16,16/)   !First Bin number that is RAMS snow for each
                        !of the 3 ice crystal types

!Coefficient for deposition-condensation freezing, and contact nuclei 
!formulations from Meyers, respectively. The code came to me (Adele) 
!with C2_MEY=0, corresponding to no contact nuclei. This should stay zero 
!unless the nucleation code is fixed. It currently does not treat the contact 
!nuclei correctly.
real, parameter :: C1_MEY=0.00012,C2_MEY=0.

! For drop freezing (iceform=1) or drop evaporation nuclei (iceform = 2) 
integer, parameter :: iceform = 2 

! KRFREEZ set bin boundary between plates and hail with freezing
! AFREEZE, BFREEZE, BFREEZEMAX - coefficients in formula of 
!                                homogeneous freezing
integer, PARAMETER :: KRFREEZ=21
REAL, PARAMETER :: BFREEZMAX=0.66E0

! other parameters and thresholds 
real, parameter :: AFREEZMY=0.3333E-04,BFREEZMY=0.6600E00

! Parameters are used in algorithm of diffusional growth
! NCOND determine timestep (DT/NCOND) with diffusional growth
integer, parameter :: NCOND=1
real :: DTCOND

! Coeffients for diffusional growth
real, parameter :: A1_MYN=2.53,BB1_MYN=5.42, A2_MYN=3.41E1,BB2_MYN=6.13
real, parameter :: AA1_MY=2.53E12,BB1_MY=5.42E3, AA2_MY=3.41E13,BB2_MY=6.13E3 
real, parameter :: DEL_BB=BB2_MY-BB1_MY, DEL_BBN=BB2_MYN-BB1_MYN, DEL_BBR=BB1_MYN/DEL_BBN

! COAGULATION
!Do COAL_BOTT_NEW every NDTCOLL time steps
integer :: NDTCOLL
real :: DTCOLL

! Kernels (collisions), depend on heights
! P1=1000mb=1000000 dynes/cm^2
real, parameter :: p1=1000000.0,p2=750000.0,p3=500000.0

! Parameters used into subroutine coal_bott (collisions):
! if temperature less than TTCOAL then no calculations of
! temperature dependence of collision kernels
real, parameter :: TTCOAL=233.15E00

! TCRIT - temperature boundary between graupel and hail with coagulation
real, parameter :: TCRIT=270.15E00

! alcr=1.0 g/m**3 :threshold for graupel formation snow-drop collisions
real, parameter :: ALCR=1.5

!  threshold for hail formation graupel-drop collisions
integer, parameter :: alcr_hail=3

! threshold bin number for hail formation from graupel-drop collisions
integer, parameter :: kp_hail=25

! scal=1.0 is valid under condition of mass doubling : 
! xl(k+1)=xl(k)*2 
real, parameter :: SCAL=1.0E00

! Parameters used for ice multiplication :
! if icempl=0 no ice multiplication; 
! if icempl=1 - ice multiplication is included
integer, parameter :: ICEMPL=1
integer, parameter :: kr_icempl=9

! 3-point remapping
integer, parameter :: I3POINT=1

! For new ice nucleation from J. Comstock (ICEFLAG =1)

!       version
!               adapted from aerosol_prop_mks.h; from J.Comstock
!       purpose
!               -define properties for aerosol composition used in
!               heterogeneous nucleation (KC1999)
!               -define aerosol size distribution parameters
!       NOTE MKS units

!       =========================================================================
! aerosol solubility in terms of volume fraction (Qv)
real, parameter :: qvaero=0.4
! aerosol parameter describing composition
real, parameter :: betafr=0.5
! molecular weight of dry aerosol (g/mol) (NH4)2SO4 Ammonium Sulfate KC1999
real, parameter :: Ms2=0.132
! density of dry aerosol (kg/m3) (NH4)2SO4 Ammonium Sulate KC1999
real, parameter :: rhos2=1770.0
! wettability
real, parameter ::  wetcoef = 0.90
! relative area of active sites
real, parameter :: alf = 0.1e-5
! misfit strain parameter
real, parameter :: epsil = 2.5e-2
! IN fraction
real :: fracin
!       ==========================================================================
!!! For ccn regeneration
real :: ccnreg
real :: inreg
!       ================================================
! For idealized ccn distributions
real :: sig_g, rp_g !cm

! For radiation:
real, dimension(nzpmax,11) :: rxb2, cxb2

END MODULE micro_prm
