!##############################################################################
real Function MCKPP_CPSW (S,T1,P0)

  ! UNITS:      
  !       PRESSURE        P0       DECIBARS
  !       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
  !       SALINITY        S        (IPSS-78)
  !       SPECIFIC HEAT   CPSW     J/(KG DEG C)
  ! ***
  ! REF: MILLERO ET AL,1973,JGR,78,4499-4507
  !       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
  ! PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
  ! DEVELOPED BY FOFONOFF 1980.
  ! ***
  ! CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
  ! T = 40 DEG C, P0= 10000 DECIBARS
  !
  !   check that temperature is above -2
  T = T1
  if(T.lt.-2.) T = -2.
  !
  !   SCALE PRESSURE TO BARS
  P=P0/10.
  ! ***
  ! SQRT SALINITY FOR FRACTIONAL TERMS
  SR = SQRT(ABS(S))
  ! SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
  A = (-1.38385E-3*T+0.1072763)*T-7.643575
  B = (5.148E-5*T-4.07718E-3)*T+0.1770383
  C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T-3.720283)*T+4217.4
  CP0 = (B*SR + A)*S + C
  ! CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
  A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T-0.49592
  B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T+2.4931E-4
  C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
  CP1 = ((C*P+B)*P+A)*P
  ! CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
  A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T+4.9247E-3
  B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
  A = (A+B*SR)*S
  B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
  B = (B+9.971E-8*SR)*S
  C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
  C = (C-1.4300E-12*T*SR)*S
  CP2 = ((C*P+B)*P+A)*P
! SPECIFIC HEAT RETURN
  MCKPP_CPSW = CP0 + CP1 + CP2

return
END FUNCTION MCKPP_CPSW

!##############################################################################
Subroutine MCKPP_ABK80 (S,T1,P,Alpha,Beta,Kappa,Sig0,Sig)

implicit none

! Routine to Coordinate computation of the expansion coefficients of 
! seawater with respect to temperature (alpha), salinity (beta), and
! pressure (kappa) using the 1980 equation of state for seawater. 
! 
! Routines Bet80,Alf80,Kap80 are called IN THAT ORDER to assure a
! minimum of repetitive work between the routines, as well as to 
! minimize storage in common (alpha overwrites coefficients that beta 
! obtains through common from Sig80). 
! 
! On entry, the calling program should specify a nonzero value for
! any or all of: alpha,beta,kappa in the argument list.  This indicates 
! to THIS routine which values are desired.
! 
! If either alpha or beta (or both) are requested, then Sig80 is called 
! to compute the density at S,T,P.  Intermediate quantities in that 
! computation are passed in common for usage by Bet80 (FIRST!) for very 
! little extra work.  Alf80 doesn't gain too much from Sig80, and must
! redo many temperature power series.  To save space, the intermediate
! coefficients in common are overwritten as Alpha is computed since 
! beta will have already used them (if it was computed at all). 
! 
! Since Kappa does not require density values to be known, but does
! require the Bulk Modulus, an entry point to compute only K in Sig80 
! has been set up.  However, if alpha or beta are computed then K and 
! other necessary terms to compute kappa will be in common and the
! call to BLKMOD in Sig80 will not be necessary.  There is, however, a 
! special case when we must still call BlkMod in Sig80: if P=0.  This 
! results because the Bulk Modulus won't be needed to compute density 
! (sigma-t) but they WILL be needed to compute kappa.  In this case 
! we again go to the entry point BlkMod in Sig80. 
! 
! Coordination of whether or not density has been computed (Sig80 
! called) is handled by logical flag KapFlg.  This is used to control 
! entry BlkMod in Sig80, as well as to minimize extra work in Kappa if
! density is already known, so that K and the other terms needed for 
! Kappa are already in common.
! 
! Units used in this routine and external routines:
! Sig80,Alf80,Bet80,Kap80 are:
!
! S(1978 practical salinity),T(degC),P(dbar),
! Alpha(degC**-1),Beta(nondim. {divided by 10**-3 for 1978 prac. sal.}),
! Kappa(bar**-1),Sig & Sig0(kg/m**3).
!
! Check Values:  S=35,T=15(degC),P=0(dbar)-->Alpha=2.14136e-4,
!                                            Beta =7.51638e-4,
!                                            Kappa=4.32576e-5.
!
!                S=40,T=0(degC),P=10,000(dbar)-->Alpha=2.69822e-4,
!                                                Beta =6.88317e-4,
!                                                Kappa=3.55271e-5.
!
! References: UNESCO(1981): Background papers and supporting data on
!                           the International Equation of State of
!                           Seawater, 1980.  UNESCO Tech. Pap. in Mar.
!                           Sci., No. 38, 192 pp.
!
!             Fofonoff, N.P., and R.C. Millard, Jr. (1983): Algorithms
!                           for computation of fundamental properties
!                           of seawater.  UNESCO Tech. Pap. in Mar.
!                           Sci., No. 44, 53 pp.
!
!             Lillibridge, J.L. (1988): Computing the seawater expansion
!                           coefficients directly from the 1980 equation
!                           of state.  Jour. Atm. Ocean. Tech., in press.
!
! John L. Lillibridge @ URI/GSO: 25 March 1987 
! 

  !    *,Exp. Coeff. of SeaWater 1980 <870330.1557> 
 
  Real P,P0,T,S,SR,Sig,Sig0,R1,R2,R3,R4,T1 
  Real PK,A,B,Alpha,Beta,Kappa  
  Real A1,B1,C,D,E,K
  Real Rho,Rho0,ABFac 
  Logical KapFlg,ABFlg
  
  Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0,ABFac,ABFlg 

  ! Check that temperature is above freezing
  T = T1
  if(T.lt.-2.) T = -2.
  
  ! Set the KapFlg to .TRUE. to indicate Sig80 not called yet 
  ! and set the ABFlg to .TRUE. to be reset only by Beta to 
  ! save alpha a little work 
  KapFlg=.True. 
  ABFlg=.True.

  ! Begin with Beta 
  ! Note that this MUST be called before Alpha due to destruction 
  ! of Intermediate sums used between Sig80 and Beta
  If(Beta.ne.0)Then      
     !Main Program wants Beta 
     CALL MCKPP_Sig80 (S,T,P,KapFlg,Sig0,Sig) 
     CALL MCKPP_Bet80 (S,T,P,Beta)
  EndIf

  ! Compute Alpha next overwriting many common variables used by Sig80  
  If(Alpha.ne.0)Then     
     !Main Program wants Alpha
     If(KapFlg)Then        
        CALL MCKPP_Sig80 (S,T,P,KapFlg,Sig0,Sig) 
     EndIf
     CALL MCKPP_Alf80 (S,T,P,Alpha) 
  EndIf
 
! If we have not called Sig80 yet, zero out Sig,Sig0 to be safe.  
  If(KapFlg)Then
     Sig=0.
     Sig0=0. 
  EndIf

! Lastly compute Kappa using quantities in common from Sig80 if 
! possible, otherwise go to entry point BlkMod in Sig80 from
! within Kap80 
  If(Kappa.ne.0)Then
     CALL MCKPP_Kap80 (S,T,P,KapFlg,Kappa)
  EndIf

! All Done, Return with appropriate desired values calculated  

return
END SUBROUTINE MCKPP_ABK80

!##############################################################################
Subroutine MCKPP_Bet80 (S,T,P,Beta)

implicit none

! Coefficient of Haline Contraction using algebraic derivation of 
! the formulae for the 1980 Equation of State.
! 
! Units: P(db),T(DegC),S(1978 Practical Salinity),Sig & Sg80(kg/m**3)
!        Rho & Rho0(kg/m**3), Beta (nondim. {divided by 10**-3 for
!        1978 prac. sal.})
! 
! Common variables are polynomials of T evaluated by Sig80 in 
! computing density, which are then used again by Beta.  In addition
! the sigma-t value (sig) and secant bulk modulus (K) are stored for
! use in the overall expression for Beta. 
! 
! John L. Lillibridge at URI/GSO: 12/22/86  

  Real P,P0,T,S,SR,R1,R2,R3,R4 
  Real Beta,PK,A,B,SR5
  Real A1,B1,C,D,E,K
  Real Rho,Rho0 
  Real DRho,DK,DK0,DA,DB
  Real ABFac
  Logical ABFlg 
  Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0,ABFac,ABFlg

! Compute the Sigma-t derivative term 
! Note that SR is the Sq.Root of Salinity from Sig80 
  SR5=SR*1.5
  DRho=R2+SR5*R3+(S+S)*R4 
  If(P.eq.0)Then
     Beta=DRho/Rho 
     Return
  EndIf

! Next Compute the Derivative Terms for the Bulk Modulus 
  DK0=A1+SR5*B1 
  DA=C+SR5*D
  DB=E 

! Derivative DK/DS 
  DK=(DB*P0+DA)*P0+DK0 

! Assemble Beta from all the terms (Rho0,K,PK from Sig80)  
  ABFac=Rho0*P0/((K-P0)*(K-P0)) 
  ABFlg=.False. 
  Beta=DRho/(1.-PK) - ABFac*DK
  Beta=Beta/Rho 

return
END SUBROUTINE MCKPP_Bet80
 
!##############################################################################
Subroutine MCKPP_Alf80 (S,T,P,Alpha)

implicit none

  ! Thermal Expansion Coefficient for Sea Water  
  !    *,1980 Thermal Exp. Coeff. <870330.1549> 
  ! 
  ! Derivative of the 1980 equation of state with respect to
  ! Temperature.  Constants in the EOS have been replaced simply
  ! with c(i)-->i*c(i), and the powers in T are reduced by one
  ! in the Horner's Rule power series sums. 
  ! 
  ! UNITS: P(DBAR), T(DEG C), S(1978 Practical Salinity),
  !        P0,K(BAR), Alpha (DegC-1)
  ! 
  ! Updated for use with alpha,beta,kappa package by JLL March 24, 1987 
  ! Including use of extensive common to save intermediate summations.
  ! 
  ! Note that Rho,Rho0,P0(Bars),K,PK(=P0/K) & SR will be defined on 
  ! entry from a previous evaluation of Rho 
  ! 
  ! John L. Lillibridge @ URI/GSO: March 24, 1987
  !

  Real P,P0,PK                             
  !Pressure Related Terms
  Real T,S,SR                              
  !Temp & Salinity Terms 
  Real R1,R2,R3,R4                         
  !Atm. Prs. Density Terms  
  REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW       
  !Bulk Modulus Terms
  Real Rho,Rho0                            
  !Sigma + 1000
  Real Alph0,AlphaA,AlphB,AlphK,Alpha      
  !Alpha Terms 
  Real ABFac                               
  !Common to Alpha/Beta
  Logical ABFlg                            
 
! Common Block for the Equation Of State routines Alf80,Bet80,Kap80
! Which can utilize Intermediate sums to minimize overhead once density 
! computed by routine Sig80.  
  COMMON /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0,ABFac,ABFlg 
 
! COMPUTE Alpha PURE WATER AT ATM PRESS 
 
  R1=(((.3268166E-7*T-.4480332e-5)*T+.3005055e-3)*T-.1819058E-1)*T+6.793952E-2 

! SEAWATER Alpha AT ATM PRESS  
  R2=((.215500E-7*T-.247401E-5)*T+.152876E-3)*T-4.0899E-3 
  R3=-.33092E-5*T+1.0227E-4 
  Alph0=(R3*SR+R2)*S+R1 
 
! Alpha AT ATM PRESS 
  IF(P.EQ.0.0)Then
     Alpha=-Alph0/Rho
     Return
  EndIf
 
! COMPUTE COMPRESSION TERMS  
  B1=-.106018E-2*T+1.6483E-2
  A1=(-.18501E-3*T+.219974E-1)*T-0.603459 
  KW=((-.2062115E-3*T+.4081431E-1)*T-.4654210E+1)*T+148.4206 
  K0=(B1*SR+A1)*S+KW
 
! Pressure Terms in Bulk Modulus 
  E=.183394E-8*T+2.0816E-8
  BW=.105574E-6*T-6.12293E-6
  AlphB=BW+E*S 
  C=-.32156E-5*T-1.0981E-5
  AW=(-.1733715E-5*T+.232184E-3)*T+1.43713E-3 
  AlphaA=C*S+AW 
 
! EVALUATE PRESS POLYNOMIAL AND RETURN 
  AlphK=(AlphB*P0+AlphaA)*P0+K0 
  If(ABFlg) Then
     ABFac=Rho0*P0/((K-P0)*(K-P0)) 
  EndIf
  Alpha=Alph0/(1.-PK) - ABFac*AlphK 
  Alpha=-Alpha/Rho

return
END SUBROUTINE MCKPP_Alf80

!##############################################################################
Subroutine MCKPP_Kap80 (S,T,P,KapFlg,Kappa)

implicit none

! Coefficient of Compressibility using algebraic derivation of
! the formulae for the 1980 Equation of State.
! 
! Units: P(db),T(degC),S(1978 Practical Salinity),Sig & Sg80(kg/m**3)
!        Rho & Rho0(kg/m**3), Kappa (bar-1) 
! 
! If density has already been computed, for alpha or beta, then 
! the coefficients A,B and the secant bulk modulus K will have
! already been computed so simply evaluate the formula for kappa. 
! However, if P=0, then neither A,B, nor K will be computed in
! Sig80, therefore we must go to that routine at entry point 
! BlkMod to obtain those coefficients.  Similarly, if density 
! has not been computed, we go to entry BlkMod in routine
! Sig80 to obtain A,B,K for kappa.
! 
! John L. Lillibridge @ URI/GSO: 25 March 1987 
 
  Real P,P0,T,S,SR,R1,R2,R3,R4 
  Real PK,A,B,Kappa 
  Real A1,B1,C,D,E,K,DelK 
  Real Rho,Rho0,ABFac
  Logical KapFlg, ABFlg                    
  ! Signals we need entry to BlkMod 

  Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0,ABFac,ABFlg 

  ! If P=0 must go to Sig80 routine at entry BlkMod
  If(P.eq.0) Then 
     KapFlg=.True. 
     CALL MCKPP_BlkMod (S,T,P,KapFlg) 
     Kappa=1.0/K                   
     ! K in Common (Bars), Kappa in bar-1 
     Return
  EndIf
   
! If Density has not been computed must also go to BlkMod 
  If(KapFlg) Then 
     CALL MCKPP_BlkMod (S,T,P,KapFlg) 
  EndIf
   
! Full nonzero pressure formula for Kappa given A,B,K  
  DelK=A+(P0+P0)*B
  Kappa=(1.-PK*DelK)/(K-P0) 

return
END SUBROUTINE MCKPP_Kap80
 
!##############################################################################
Subroutine MCKPP_Sig80 (S,T,P,KapFlg,Sig0,Sig)

implicit none

! EQUATION OF STATE FOR SEA WATER  
!    *,1980 EQ. OF STATE <870330.1544>  
! 
! EQUATION OF STATE FOR SEAWATER PROPOSED BY JPOTS 1980 
! 
! REF: MILLERO, ET. AL. (1980) D.S.R.,27A,255-264.
! UNITS: P(DBAR), T(DEG C), S(P.S.U.), Sig,Sig0(KG/M**3)  
!        P0,K(BAR)
! 
! NPF,OCT 7 80; JLL,FEB 24 82.
! UPDATED TO USE ONLY SIGMA UNITS BY JLL MAR 10 82. 
! 
! Updated for use with alpha,beta,kappa package by JLL March 24, 1987 
! Including use of extensive common to save intermediate summations.
! 
! John L. Lillibridge @ URI/GSO: March 24 1987

  Real P,P0,PK                             
  !Pressure Related Terms
  Real T,S,SR                              
  !Temp & Salinity Terms 
  Real Sig0,R1,R2,R3,R4                    
  !Atm. Prs. Density Terms 
  REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW       
  !Bulk Modulus Terms
  Real BlkMod                              
  !entry Point for Kappa 
  Real Sig,Rho,Rho0,ABFac                    
  !Sigma + 1000
  Logical KapFlg, ABFlg                           
  !True when entry at BlkMod 

! Common Block for the Equation Of State routines Alf80,Bet80,Kap80
! Which can utilize Intermediate sums to minimize overhead once density 
! computed by this routine.  
  Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0,ABFac,ABFlg 
 
! CONVERT PRESSURE TO BARS AND SQ ROOT SALINITY 
! And Set the Logical Flag used to control entry Point Below for Kappa 
  P0=P/10.0 
  SR=SQRT(ABS(S)) 
  KapFlg=.False.
 
  ! COMPUTE SIGMA PURE WATER AT ATM PRESS  
  R1=((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T-9.095290E-3)&
       *T+6.793952E-2)*T-.157406 
  !Note Constant for Sigma vs. Rho 
  !SEAWATER SIGMA AT ATM PRESS  
  R2=(((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T+8.24493E-1 
  R3=(-1.6546E-6*T+1.0227E-4)*T-5.72466E-3
  R4=4.8314E-4
  Sig0=(R4*S+R3*SR+R2)*S+R1 
  Rho0=1000.0 + Sig0
 
  ! SIGMA AT ATM PRESS 
  IF(P.EQ.0.0)Then
     Sig=Sig0  
     Rho=Rho0
     Return
  EndIf
 
  ! COMPUTE COMPRESSION TERMS 
  
  ! This is the entry point when the Bulk Modulus K is desired, which 
  ! does not require knowledge of density. entry from Kappa will have
  ! KapFlg set .True. so that we exit before computing Sig.  
  Entry MCKPP_BlkMod (S,T,P,KapFlg)
  
  ! For Computation of Kappa  
  
  ! If entry without density computed need these two common terms 
  
  If(KapFlg)Then
     P0=P/10.0 
     SR=SQRT(ABS(S)) 
  EndIf
  
  B1=(-5.3009E-4*T+1.6483E-2)*T+7.944E-2
  A1=((-6.1670E-5*T+1.09987E-2)*T-0.603459)*T+54.6746 
  KW=(((-5.155288E-5*T+1.360477E-2)*T-2.327105)*T+148.4206)*T+19652.21 
  K0=(B1*SR+A1)*S+KW
  
  ! If computing Kappa and P=0, return here 
  If(P.eq.0.0)Then
     K=K0
     Return
  EndIf

! Pressure Terms in Bulk Modulus 
  E=(9.1697E-10*T+2.0816E-8)*T-9.9348E-7
  BW=(5.2787E-8*T-6.12293E-6)*T+8.50935E-5
  B=BW+E*S
  
  D=1.91075E-4
  C=(-1.6078E-6*T-1.0981E-5)*T+2.2838E-3
  AW=((-5.77905E-7*T+1.16092E-4)*T+1.43713E-3)*T+3.239908 
  A=(D*SR+C)*S+AW 

! EVALUATE PRESS POLYNOMIAL AND RETURN 
  K=(B*P0+A)*P0+K0
   
! If entry at BlkMod Exit Now 
  PK=P0/K 
  If (KapFlg) Return
  Sig=(1000.0*PK+Sig0)/(1.0-PK)   
  Rho=1000.0 + Sig  

return
END SUBROUTINE MCKPP_Sig80

