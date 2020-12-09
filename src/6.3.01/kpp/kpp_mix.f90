!##############################################################################
Subroutine mckpp_physics_ocnstep (kpp_1d_fields,kpp_const_fields,i,j)

use mem_kpp, only:kpp_1d_type,kpp_const_type
use kpp_parameters, only:kppNZ,kppNZP1,kppNVEL,kppNSCLR,kppitermax,kpprnt &
                        ,iprnt,jprnt
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid,time

implicit none

  !-----------------------------------------------------------------------  
  !     Main driver for ocean module.
  !     Integration is performed only on the permanent grid
  !     Written   3 Mar 1991 - WGL
  !     Modified  5 Jun 1992 - jan : implicit scheme
  !              16 Nov      - jan : latest version
  !              16 Nov 1994 - wgl : new KPP codes no temporary grid
  
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  real hmixe,hmixn,tol
  real Uo(kppNZP1,kppNVEL),Xo(kppNZP1,kppNSCLR)
  real Ux(kppNZP1,kppNVEL),XX(kppNZP1,kppNSCLR) 
  ! Additional variables to provide smoothing in the iteration.

  real Ui                   ! Ui used in damping (LH 8/08/2013)
  real dampU(kppNVEL)       ! dampU used to flag which Ui chosen in damping
  real lambda               ! Factor to control smoothing
  integer iter,iconv        ! number of iterations
  integer kmixe,kmixn
  
  ! More Local Variables (to make implicit none)
  real deltaz,a,b
  integer k,l,n,i,j

  !gravity constant
  real, parameter :: grav=9.816

  !For energetics computations (not done here; could uncomment 
  !if you want to compute these diagnostics (Saleeby 2018)
  !real rhonot,eflx,esnk,Tmke,Ptke,rmke(kppNZP1)  

  ! Number of iterations for computational instability
  integer comp_iter_max
  real rmsd(4),rmsd_threshold(4)
  data comp_iter_max /10/
  ! Critical depth-integrated RMS difference between old and new profiles
  ! for repeating integration, for (/U,V,T,S/). Typical (stable) 
  ! values are O(10^-2) for U and V, O(10^-3) for T and O(10^-4) for S.
  ! NPK 17/5/13
  data rmsd_threshold /1,1,1,1/
  
  data lambda /0.5/
  
  Uo=kpp_1d_fields%U(:,:)
  Xo=kpp_1d_fields%X(:,:)
  kpp_1d_fields%comp_flag=.TRUE.
  kpp_1d_fields%reset_flag=0
  
  DO WHILE (kpp_1d_fields%comp_flag .and. kpp_1d_fields%reset_flag &
       .le. comp_iter_max)
     !Print some things for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
       print*,'ITERATION-0',kpp_1d_fields%comp_flag,kpp_1d_fields%reset_flag &
         ,comp_iter_max
     endif
     ! Estimate new profiles by extrapolation
     DO k=1,kppNZP1
      DO l=1,kppNVEL
       IF (nint(kpp_1d_fields%old) .lt. 0 .or. nint(kpp_1d_fields%old) .gt. 1) THEN
        WRITE(*,*) 'Dodgy value of old at k=',k,'l=',l,'old=',nint(kpp_1d_fields%old)
        kpp_1d_fields%old=nint(kpp_1d_fields%new)
        stop
       ENDIF
       IF (nint(kpp_1d_fields%new) .lt. 0 .or. nint(kpp_1d_fields%new) .gt. 1) THEN
        WRITE(*,*) 'Dodgy value of new at k=',k,'l=',l,'new=',nint(kpp_1d_fields%new)
        kpp_1d_fields%new=nint(kpp_1d_fields%old)
        stop
       ENDIF        
       kpp_1d_fields%U(k,l)=2.*kpp_1d_fields%Us(k,l,nint(kpp_1d_fields%new))- &
            kpp_1d_fields%Us(k,l,nint(kpp_1d_fields%old))
       Ux(k,l)=kpp_1d_fields%U(k,l)
      ENDDO
      DO l=1,kppNSCLR
       kpp_1d_fields%X(k,l)=2.*kpp_1d_fields%Xs(k,l,nint(kpp_1d_fields%new))- &
            kpp_1d_fields%Xs(k,l,nint(kpp_1d_fields%old))
       Xx(k,l)=kpp_1d_fields%X(k,l)
      ENDDO
     ENDDO
     
     ! Iteration loop for semi-implicit integration
     ! Reset iteration counter
     iter=0
     iconv=0
     
     ! This loop controls the number of compulsory iterations
     ! The original value was 2, using (the upper value of the loop +2)
     ! added by SJW (17 Jan 03) to try an alleviate some non-convergences
     DO iter=0,2
       DO k=1,kppNZP1
         DO l=1,kppNVEL
          kpp_1d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*kpp_1d_fields%U(k,l)
          Ux(k,l)=kpp_1d_fields%U(k,l)
         ENDDO
         DO l=1,kppNSCLR
          kpp_1d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*kpp_1d_fields%X(k,l)
          Xx(k,l)=kpp_1d_fields%X(k,l)
         ENDDO
       ENDDO
       CALL MCKPP_PHYSICS_VERTICALMIXING (kpp_1d_fields,kpp_const_fields &
                                         ,hmixe,kmixe,i,j)
       CALL MCKPP_PHYSICS_OCNINT (kpp_1d_fields,kpp_const_fields,kmixe,Uo,Xo,i,j)
       !Print some things for debugging (Saleeby2018)
       if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
         print*,'ITERATION-1',iter,hmixe,kmixe
       endif
     ENDDO
     ! The original code can be restored by reseting iter=1 and removing the  
     ! above loop  
     ! iter=1
         
     IF (kpp_const_fields%LKPP) THEN
45      continue
        DO k=1,kppNZP1
           DO l=1,kppNVEL
            kpp_1d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*kpp_1d_fields%U(k,l)
            Ux(k,l)=kpp_1d_fields%U(k,l)
           ENDDO
           DO l=1,kppNSCLR
            kpp_1d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*kpp_1d_fields%X(k,l)
            Xx(k,l)=kpp_1d_fields%X(k,l)
           ENDDO
        ENDDO
        CALL MCKPP_PHYSICS_VERTICALMIXING (kpp_1d_fields,kpp_const_fields &
                                          ,hmixn,kmixn,i,j)       
        CALL MCKPP_PHYSICS_OCNINT (kpp_1d_fields,kpp_const_fields,kmixn,Uo,Xo,i,j)
        iter = iter + 1
         
        ! prepare tolerances to check for convergence
        tol = kpp_const_fields%hmixtolfrac*kpp_const_fields%hm(kmixn)
        if(kmixn.eq.kppNZP1) tol = kpp_const_fields%hmixtolfrac &
                              *kpp_const_fields%hm(kppNZ)

        !Print some things for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
          print*,'ITERATION-2',iter,iconv,hmixn,kmixn
          print*,'ITERATION-3',abs(hmixn-hmixe),tol,time,hmixn
        endif

        ! check iteration for convergence
        if(abs(hmixn-hmixe).gt.tol)  then
         ! Uncommeting the following the lines iconv=0 to IF (iconv ...)
         ! will make the model do two consecutive tests for convergence of the 
         ! hmix (added by SJW 17 Jan 03). This did not work well in testing for
         ! long timestep, high resolution (the model generally failed to satisfy 
         ! convergence test on two consecutive iterations.
         iconv=0
        ELSE
         iconv=iconv+1
        ENDIF
        IF (iconv .lt. 3) THEN
           if (iter.lt.kppitermax) then
              hmixe = hmixn
              kmixe = kmixn
              goto 45
           else
              ! use shallower hmix
              if(hmixn.gt.hmixe) then
                 hmixe = hmixn
                 kmixe = kmixn
                 goto 45
              endif
           endif
        endif

        if( iter.gt.(kppitermax+1) ) then 
         write(*,1009) i+mi0(ngrid),j+mj0(ngrid),hmixe,hmixn,hmixn-hmixe &
             ,kmixn,iter
1009     format('i=',i5,' j=',i5,' hmixest=',f7.2,' hmixnew=',f7.2 &
           ,' diff=',f6.1,' kmixn=',i3,' iteration=',i3)
        endif
     ENDIF

     ! Trap for profiles that are very different from original profile
     ! or clearly erroneous, to detect rare instances of instability
     ! in the semi-implicit integration.  Reset to original profile,
     ! add some noise via changing Coriolis term slightly, and try
     ! integration again.
     ! NPK 16/5/2013
     kpp_1d_fields%comp_flag=.FALSE.
     DO k=1,kppNZ
        IF (ABS(kpp_1d_fields%U(k,1)).ge. 10 &
            .or. ABS(kpp_1d_fields%U(k,2)).ge.10 .or. &
             ABS(kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1)) .ge. 10) THEN 
           kpp_1d_fields%comp_flag=.TRUE.
           kpp_1d_fields%f=kpp_1d_fields%f*1.01
        ENDIF
     END DO
     IF (.NOT. kpp_1d_fields%comp_flag) THEN
        rmsd(:)=0.
        DO k=1,kppNZP1
           rmsd(1)=rmsd(1)+(kpp_1d_fields%U(k,1)-Uo(k,1))*(kpp_1d_fields%U(k,1) &
               -Uo(k,1))*kpp_const_fields%hm(k)/kpp_const_fields%dm(kppNZ)
           rmsd(2)=rmsd(2)+(kpp_1d_fields%U(k,2)-Uo(k,2))*(kpp_1d_fields%U(k,2) &
               -Uo(k,2))*kpp_const_fields%hm(k)/kpp_const_fields%dm(kppNZ)
           rmsd(3)=rmsd(3)+(kpp_1d_fields%X(k,1)-Xo(k,1))*(kpp_1d_fields%X(k,1) &
               -Xo(k,1))*kpp_const_fields%hm(k)/kpp_const_fields%dm(kppNZ)
           rmsd(4)=rmsd(4)+(kpp_1d_fields%X(k,2)-Xo(k,2))*(kpp_1d_fields%X(k,2) &
               -Xo(k,2))*kpp_const_fields%hm(k)/kpp_const_fields%dm(kppNZ)
        ENDDO
        DO k=1,4
           rmsd(k)=SQRT(rmsd(k))
           IF (rmsd(k).ge.rmsd_threshold(k)) THEN
              kpp_1d_fields%comp_flag=.TRUE.
              kpp_1d_fields%f=kpp_1d_fields%f*1.01
           ENDIF
        ENDDO
     ENDIF
     kpp_1d_fields%reset_flag=kpp_1d_fields%reset_flag+1

     IF (kpp_1d_fields%reset_flag .gt. comp_iter_max) THEN
        print*,'Failed to find a reasonable solution in the semi-implicit &
             integration after ',comp_iter_max,' iterations.'
        print*,'At point =',i+mi0(ngrid),j+mj0(ngrid),kpp_1d_fields%X(1,1)
        WRITE(*,*) 'U = ',kpp_1d_fields%U(:,1)
        WRITE(*,*) 'V = ',kpp_1d_fields%U(:,2)
        WRITE(*,*) 'T = ',kpp_1d_fields%X(:,1)
        WRITE(*,*) 'S = ',kpp_1d_fields%X(:,2)
        stop            
     ENDIF
     ! End of trapping code.

  ENDDO !end DO WHILE
  
  ! Output Results from permanent grid iterations
  ! Compute diagnostic fluxes for writing
  do k=1,kppNZ
     deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
     do n=1,kppNSCLR
        kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*((kpp_1d_fields%X(k,n) &
            -kpp_1d_fields%X(k+1,n))/deltaz &
            -kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
        !Print some info for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1.and.n==1 &
          .and.k<=5) then
         print*,'Steve-WX',kpp_1d_fields%wX(0,n),kpp_1d_fields%wX(k,n) &
          ,kpp_1d_fields%difs(k),kpp_1d_fields%ghat(k)
        endif
     enddo
     if(kpp_const_fields%LDD) then
      kpp_1d_fields%wX(k,1)= -kpp_1d_fields%dift(k)*((kpp_1d_fields%X(k,1) &
         - kpp_1d_fields%X(k+1,1))/deltaz-kpp_1d_fields%ghat(k) &
         * kpp_1d_fields%wX(0,1))
     endif
     kpp_1d_fields%wB(k)= grav * (kpp_1d_fields%talpha(k) &
         *kpp_1d_fields%wX(k,1) - kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
     do n=1,kppNVEL
        kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*(kpp_1d_fields%U(k,n) &
           -kpp_1d_fields%U(k+1,n))/deltaz
     enddo
  enddo

  ! Compute energetics: Saleeby(2018):could uncomment this if we want these
  ! diagnostic for energy
!  rhonot = 1026.
!  Eflx = 0.5 * ( (Uo(1,1) + kpp_1d_fields%U(1,1)) * &
!       kpp_1d_fields%sflux(1) + (Uo(1,2) + kpp_1d_fields%U(1,2)) * &
!       kpp_1d_fields%sflux(2) )
!  Esnk = -0.5*rhonot* ( (Uo(kppNZ,1) + kpp_1d_fields%U(kppNZ,1)) &
!        * kpp_1d_fields%wU(kppNZ,1) + &
!         (Uo(kppNZ,2) + kpp_1d_fields%U(kppNZ,2)) * kpp_1d_fields%wU(kppNZ,2) )
!  Ptke = 0.0
!  ! use "amax1" to prevent "underflow" in single precision
!  do k=1,kppNZ-1
!   Ptke = Ptke - 0.5*( amax1(kpp_1d_fields%wU(k,1),1.E-10)* &
!    (rhonot   * (Uo(k,1) + kpp_1d_fields%U(k,1)) - rhonot * &
!    (Uo(k+1,1) + kpp_1d_fields%U(k+1,1)) ) +amax1(kpp_1d_fields%wU(k,2),1.E-10)*&
!    (rhonot   * (Uo(k,2) + kpp_1d_fields%U(k,2)) -rhonot * (Uo(k+1,2) &
!     + kpp_1d_fields%U(k+1,2)) ) )
!  ENDDO
!  Tmke = 0.0
!  do k=1,kppNZP1
!     rmke(k) = 0.5 * rhonot * (kpp_1d_fields%U(k,1)**2 &
!             + kpp_1d_fields%U(k,2)**2) * kpp_const_fields%hm(k)
!     Tmke = Tmke + rmke(k)
!  ENDDO
            
!  check heat and salt budgets
!  call budget(Xo,kpp_1d_fields,kpp_const_fields)
!            
!  Set new profiles
!          do k=1,kppNZP1     ! values at kppNZP1 only change for slab ocean
!             do n=1,kppNVEL
!                kpp_1d_fields%U(k,n) = kpp_1d_fields%U(k,n)
!             enddo
!             do n=1,kppNSCLR
!                kpp_1d_fields%X(k,n) = kpp_1d_fields%X(k,n)
!             enddo
!          enddo
!  Set correct surface values, and CP and rho profiles for new profiles
!  Get latest profiles
!  Removed to ensure that hmix,diff,cp,rho diagnostics are consistent with 
!  those used in the final iteration of the timestep 
!  (CP,RHO are not quite correct for updated values, but are the values
!  by the integration) (SJW 16/01/03)
!  Need to consider improving convergence test!!! (SJW 16/01/03)
!                 
! The final call to vmix is removed to ensure that the diffusion and 
! boundary layer profiles in the diagnostics are the ones used to calculate
! the fluxes, as it stands at the moment this means that the CP and rho are
! also the values used in the timestepping not the values appropriate to the
! S,T at the new time level.  
! call vmix(Un,Xn,hmixe,kmixe)  
! write(*,*) time,iter, hmixn,hmixe,kmixn,kmixe

  kpp_1d_fields%hmix = hmixn
  kpp_1d_fields%kmix = kmixn
  kpp_1d_fields%uref = kpp_1d_fields%U(1,1)
  kpp_1d_fields%vref = kpp_1d_fields%U(1,2)
  kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)
  
  ! Damping currents, Added LH (06/08/2013)
  IF (kpp_const_fields%L_DAMP_CURR) THEN 
   dampu(:)=0.
   do k=1,kppNZP1
    do l=1,kppNVEL
     a=0.99*ABS(kpp_1d_fields%U(k,l))
     b=kpp_1d_fields%U(k,l)**2 * kpp_const_fields%dt_uvdamp &
        * kpp_const_fields%dto
     Ui=MIN(a,b)
     ! LH (29/08/2013) Add Flags to check which Ui (a or b) is chosen, 
     ! dtuvdamp=360 (specified in namelist). 
     ! The flags for u and v can be requested as diagnostics dampu_flag, 
     ! dampv_flag. Note that the value of the flag is equal to 
     ! the *fraction* of levels at that point where (U**2)/r < alpha*ABS(U), 
     ! 1.0=all Ui are (U**2)/r.

     ! U**2/r < alpha*U =1, U**2/r > alpha*U =-1. Output is total of all 
     ! levels. If dampu_flag =kppNZ all levels U**2/r<alpha*U.
     ! (alpha = 0.99, r=tau*(86400/dto), tau=360).
     IF (b .lt. a) THEN  
        dampU(l)=dampU(l)+1.0/REAL(kppNZP1)
     ENDIF
     ! Apply damping
     kpp_1d_fields%U(k,l)= kpp_1d_fields%U(k,l) - SIGN(Ui,kpp_1d_fields%U(k,l))
    enddo
   enddo
   kpp_1d_fields%dampu_flag=dampU(1)
   kpp_1d_fields%dampv_flag=dampU(2)
  ENDIF

  ! Save variables for next timestep
  kpp_1d_fields%old = nint(kpp_1d_fields%new)
  kpp_1d_fields%new = 1 - nint(kpp_1d_fields%old)
  do k=1,kppNZP1
     do l=1,kppNVEL
        kpp_1d_fields%Us(k,l,nint(kpp_1d_fields%new))=kpp_1d_fields%U(k,l)
     enddo
     do l=1,kppNSCLR
        kpp_1d_fields%Xs(k,l,nint(kpp_1d_fields%new))=kpp_1d_fields%X(k,l)
     enddo
  enddo
  
return
END SUBROUTINE mckpp_physics_ocnstep

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING (kpp_1d_fields,kpp_const_fields &
                                        ,hmixn,kmixn,i,j)
use mem_kpp, only: kpp_1d_type,kpp_const_type
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid
use kpp_parameters, only:kppNZ,kppNZP1,kpprnt,iprnt,jprnt

implicit none

  !  Interface between 1-d model and vertical mixing

! inputs including those from common.inc and parameter.inc
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields
  real B0,B0sol,ustar

! outputs
  real hmixn                ! boundary layer depth (m)
  integer kmixn

! local
  real dVsq(kppNZP1)    ! (velocity shear re sfc)^2      (m/s)^2
  real Ritop(kppNZ)     ! numerator of bulk Richardson Number (m/s)^2
  real alphaDT(kppNZ)   ! alpha * DT across interfaces
  real betaDS(kppNZ)    ! beta  * DS across interfaces
  real epsilon,epsln
  real alpha,beta,exppr
  real sigma,sigma0
  real, external :: MCKPP_CPSW
  real tau
  real zref,wz,bref
  integer k,n,kl,i,j
  real del
  real rhob           !density of ocean corrected for seaice salinity
  real rhoh2o         !density of ocean corrected for temperature
  real dlimit,vlimit
  integer jerl(12)

  !gravity constant
  real, parameter :: grav=9.816

  ! month  1   2   3   4   5   6   7   8   9   10  11  12
  data jerl / 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /

  epsilon = 0.1
  epsln   = 1.e-20

  ! calculate density of fresh water and brine in surface layer
  alpha  = 1.
  beta   = 1.
  exppr  = 0.0
  sigma0 = 0.0
  sigma  = 0.0

  CALL MCKPP_ABK80 (0.0,kpp_1d_fields%X(1,1),-kpp_const_fields%zm(1),alpha &
                   ,beta,exppr,sigma0,sigma)
  rhoh2o = 1000. + sigma0

  CALL MCKPP_ABK80 (kpp_const_fields%SICE,kpp_1d_fields%X(1,1) &
                   ,-kpp_const_fields%zm(1),alpha,beta,exppr,sigma0,sigma)
  rhob = 1000. + sigma0

  ! calculate temperature and salt contributions of buoyancy gradients  
  ! calculate buoyancy profile (m/s**2) on gridlevels
  do k=1,kppNZP1
     CALL MCKPP_ABK80 (kpp_1d_fields%X(k,2)+kpp_1d_fields%Sref &
         ,kpp_1d_fields%X(k,1),-kpp_const_fields%zm(k) &
         ,alpha,beta,exppr,sigma0,sigma)
     kpp_1d_fields%rho(k)= 1000. + sigma0
     kpp_1d_fields%CP(k) = MCKPP_CPSW(kpp_1d_fields%X(k,2) &
          +kpp_1d_fields%Sref,kpp_1d_fields%X(k,1),-kpp_const_fields%zm(k))
     kpp_1d_fields%talpha(k) = alpha
     kpp_1d_fields%sbeta(k)  = beta
     kpp_1d_fields%buoy(k) = -grav * sigma0 / 1000.

     !Print some things for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1.and.k<=5) then
       print*,'Steve-Sigma',kpp_1d_fields%talpha(k),kpp_1d_fields%sbeta(k) &
         ,sigma0,sigma
     endif
  ENDDO        

  ! Call to ntflx, put here to allow removal of diagnostic call to vmix
  ! and to ensure the most recent cp,rho used. Assume surface rho and cp
  ! are equal to rho(1) and cp(1). No need to create rho(0) or cp(0).
  do k=0,kppNZ
        !Large et al.Eqn.A3c, right term, sort of
        kpp_1d_fields%wXNT(k)=-kpp_1d_fields%sflux(3) &
             *kpp_1d_fields%swdk_opt(k)&
             /(kpp_1d_fields%rho(1)*kpp_1d_fields%CP(1))
  enddo

  ! calculate kinematic surface momentum fluxes(Large et al.Eqn.A2a,A2b)
  kpp_1d_fields%wU(0,1) = -kpp_1d_fields%sflux(1) / kpp_1d_fields%rho(1)
  kpp_1d_fields%wU(0,2) = -kpp_1d_fields%sflux(2) / kpp_1d_fields%rho(1)
  tau = sqrt( kpp_1d_fields%sflux(1)**2 + kpp_1d_fields%sflux(2)**2 ) &
        + 1.e-16
  !  1.e-16 added to stop subsequent division by zero if tau=0.0
  ustar = sqrt( tau / kpp_1d_fields%rho(1) )

  ! total turbulent kinematic temperature flux (C m/s)(Large et al.Eqn.A2c)
  kpp_1d_fields%wX(0,1)  = -kpp_1d_fields%sflux(4) / kpp_1d_fields%rho(1) &
     / kpp_1d_fields%CP(1)
 
  ! total turbulent kinematic salinity flux (o/oo m/s)(Large et al.Eqn.A2d)
  kpp_1d_fields%wX(0,2) = kpp_1d_fields%Sref*kpp_1d_fields%sflux(6)/&
       rhoh2o+(kpp_1d_fields%Sref-kpp_const_fields%SICE)*&
       kpp_1d_fields%sflux(5)/rhob

  ! calculate total kinematic surface buoy flux (m**2/s**3)
  ! (Large et al.Eqn.A3d, left term)
  B0 = -grav*(kpp_1d_fields%talpha(1)*kpp_1d_fields%wX(0,1) &
       - kpp_1d_fields%sbeta(1)*kpp_1d_fields%wX(0,2) )

  ! calculate total kinematic surface buoy flux (m**2/s**3)
  ! (Large et al.Eqn.A3b)
  kpp_1d_fields%wB(0) =  - B0

  ! radiative contribution to the sfc buoy forcing (Large et al.Eqn.A3c left term)
  B0sol = grav * kpp_1d_fields%talpha(1) &
    * kpp_1d_fields%sflux(3) / (kpp_1d_fields%rho(1) * kpp_1d_fields%CP(1))

  ! calculate temperature and salt contributions of buoyancy gradients on 
  ! interfaces for double diffusion      
  do n = 1,kppNZ
     alphaDT(n) =0.5 *(kpp_1d_fields%talpha(n)+kpp_1d_fields%talpha(n+1)) * &
          (kpp_1d_fields%X(n,1) - kpp_1d_fields%X(n+1,1))
     betaDS(n)  =0.5 *(kpp_1d_fields%sbeta(n) + kpp_1d_fields%sbeta(n+1)) * &
          (kpp_1d_fields%X(n,2) - kpp_1d_fields%X(n+1,2))
  ENDDO

  ! compute buoyancy and shear profiles
  DO n = 1,kppNZ
     zref =  epsilon * kpp_const_fields%zm(n)
     ! compute reference buoyancy and velocity
     wz = AMAX1(kpp_const_fields%zm(1),zref) 
     kpp_1d_fields%uref  = kpp_1d_fields%U(1,1) * wz / zref
     kpp_1d_fields%vref  = kpp_1d_fields%U(1,2) * wz / zref
     bref  = kpp_1d_fields%buoy(1)* wz / zref
     do kl = 1,kppNZ
        IF(zref.ge.kpp_const_fields%zm(kl)) go to 126
        wz = AMIN1(kpp_const_fields%zm(kl)-kpp_const_fields%zm(kl+1) &
             ,kpp_const_fields%zm(kl)-zref) 
        del = 0.5 * wz / (kpp_const_fields%zm(kl) - kpp_const_fields%zm(kl+1))
        kpp_1d_fields%uref = kpp_1d_fields%uref - wz*( kpp_1d_fields%U(kl,1) &
          + del * (kpp_1d_fields%U(kl+1,1)- kpp_1d_fields%U(kl,1))) /zref
        kpp_1d_fields%vref = kpp_1d_fields%vref - wz*( kpp_1d_fields%U(kl,2) &
          + del * (kpp_1d_fields%U(kl+1,2)- kpp_1d_fields%U(kl,2))) /zref
        bref=bref -wz*(kpp_1d_fields%buoy(kl) + del * &
             (kpp_1d_fields%buoy(kl+1)-kpp_1d_fields%buoy(kl))) /zref
     ENDDO
126  continue
     !Ritop is top of Eqn.21 in Large et al.
     Ritop(n) = (zref - kpp_const_fields%zm(n)) * (bref - kpp_1d_fields%buoy(n))
     ! NPK Additions (25/9/2008). Prevent Ritop from going negative.
     ! IF (Ritop(ipt,n) .lt. 0) Ritop(ipt,n) = epsln
     kpp_1d_fields%dbloc(n) = kpp_1d_fields%buoy(n) - kpp_1d_fields%buoy(n+1)
     dVsq(n)  = (kpp_1d_fields%Uref - kpp_1d_fields%U(n,1))**2 &
         + (kpp_1d_fields%Vref - kpp_1d_fields%U(n,2))**2
     kpp_1d_fields%shsq(n)  = (kpp_1d_fields%U(n,1)-kpp_1d_fields%U(n+1,1))**2 &
         + (kpp_1d_fields%U(n,2)-kpp_1d_fields%U(n+1,2))**2         
  ENDDO

  !Print some things for debugging (Saleeby2018)
  if(KPPRNT==1)then
   DO n = 1,kppNZ
    if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.n<=5) then
     print*,'Steve-vertmix',kpp_1d_fields%dbloc(n),kpp_1d_fields%buoy(n) &
       ,Ritop(n),kpp_1d_fields%Shsq(n)
    endif
   ENDDO
  endif

  CALL MCKPP_PHYSICS_VERTICALMIXING_KPPMIX (kppNZ,kppNZP1,dVsq,ustar,B0,B0sol &
    ,alphaDT,betaDS,Ritop,hmixn,kmixn,kpp_1d_fields,kpp_const_fields,i,j)

  dlimit = 0.00001
  vlimit = 0.0001
 
  do k=kppNZ,kppNZP1
     kpp_1d_fields%difm(k) = vlimit
     kpp_1d_fields%difs(k) = dlimit
     kpp_1d_fields%dift(k) = dlimit
  enddo
  kpp_1d_fields%ghat(kppNZ) = 0.0

return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_KPPMIX (km,kmp1,dVsq,ustar,Bo,Bosol &
    ,alphaDT, betaDS, Ritop, hbl, kbl, kpp_1d_fields, kpp_const_fields,i,j)

use mem_kpp, only: kpp_1d_type,kpp_const_type
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid
use kpp_parameters, only:kpprnt,iprnt,jprnt

implicit none

!.......................................................................
!
!     Main driver routine for kpp vertical mixing scheme and 
!     interface to greater ocean model
!
!     written by : bill large,   june,  6, 1994
!     modified by: jan morzel,   june, 30, 1994
!                  bill large, august, 11, 1994
!                  bill large, november 1994,   for 1d code
!
!.......................................................................

  integer km,kmp1,mdiff,ki,i,j
  
  !parameter (km = kppNZ, kmp1 = kppNZP1)!, imt = 1) !NX*NY)
  parameter (mdiff = 3)  ! number of diffusivities for local arrays

  ! input
  real dVsq(kmp1)           ! (velocity shear re sfc)^2      (m/s)^2
  real ustar                ! surface friction velocity       (m/s)
  real Bo                   ! surface turbulent buoy. forcing (m^2/s^3)
  real Bosol                ! radiative buoyancy forcing      (m^2/s^3)
  real alphaDT(kmp1)        ! alpha * DT  across interfaces
  real betaDS(kmp1)         ! beta  * DS  across interfaces
  real Ritop(km)            ! numerator of bulk Richardson Number (m/s)^2
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields

  ! output
  ! visc replaced by kpp_1d_fields%difm (NPK 8/2/2013)
  real hbl                  ! boundary layer depth (m)
  integer kbl               ! index of first grid level below hbl     
  
  ! local
  real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
  real ws                   ! momentum velocity scale
  real wm                   ! scalar   velocity scale 
  real caseA                ! = 1 in case A; =0 in case B
  real stable               ! = 1 in stable forcing; =0 in unstable
  real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
  real gat1(mdiff)          ! shape function at sigma=1
  real dat1(mdiff)          ! derivative of shape function at sigma=1
  real blmc(km,mdiff)       ! boundary layer mixing coefficients
  real sigma                ! normalized depth (d / hbl)
  real Rib(2)               ! bulk Richardson number
      
  ! zero the mixing coefficients 
  DO ki=0,km
     kpp_1d_fields%difm(ki) = 0.0
     kpp_1d_fields%difs(ki) = 0.0
     kpp_1d_fields%dift(ki) = 0.0
  END DO

  ! compute RI and IW interior diffusivities everywhere
  IF(kpp_const_fields%LRI) THEN
    CALL MCKPP_PHYSICS_VERTICALMIXING_RIMIX (km,kmp1,kpp_1d_fields &
                                            ,kpp_const_fields,i,j)
  ENDIF

  ! add double diffusion if desired
  IF(kpp_const_fields%LDD) THEN
    CALL MCKPP_PHYSICS_VERTICALMIXING_DDMIX (km,kmp1,alphaDT,betaDS &
                                            ,kpp_1d_fields)
  ENDIF

  ! fill the bottom kmp1 coefficients for blmix      
  kpp_1d_fields%difm(kmp1) = kpp_1d_fields%difm(km)
  kpp_1d_fields%difs(kmp1) = kpp_1d_fields%difs(km)
  kpp_1d_fields%dift(kmp1) = kpp_1d_fields%dift(km)

  ! compute boundary layer mixing coefficients
  IF(kpp_const_fields%LKPP) THEN

     !Print some info for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
      DO ki=0,4
        print*,'Steve-KPPM1',kpp_1d_fields%difs(ki),kpp_1d_fields%dift(ki) &
        ,kpp_1d_fields%difm(ki)
      enddo
     endif

     ! diagnose the new boundary layer depth
     CALL MCKPP_PHYSICS_VERTICALMIXING_BLDEPTH (km, kmp1, dVsq, Ritop, ustar, &
          Bo, Bosol, hbl, bfsfc, stable, caseA, kbl, Rib, sigma, wm, &
          ws, kpp_1d_fields, kpp_const_fields,i,j)

     ! compute boundary layer diffusivities
     CALL MCKPP_PHYSICS_VERTICALMIXING_BLMIX (km, mdiff, ustar, bfsfc, hbl, &
          stable, caseA, kbl, gat1 , dat1 , dkm1, blmc, sigma, wm, ws, &
          kpp_1d_fields,kpp_const_fields,i,j)

     !Print some info for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
      DO ki=0,4
        print*,'Steve-KPPM2',kpp_1d_fields%difm(ki),blmc(ki,1) &
           ,kpp_1d_fields%difs(ki),blmc(ki,2)
      enddo
     endif

     ! enhance diffusivity at interface kbl - 1
     CALL MCKPP_PHYSICS_VERTICALMIXING_ENHANCE  (km, mdiff, dkm1, hbl, kbl, &
          caseA, blmc, kpp_1d_fields,kpp_const_fields)

     ! combine interior and boundary layer coefficients and nonlocal term
     do ki= 1,km
        if(ki.lt.kbl) then
           kpp_1d_fields%difm(ki)=blmc(ki,1)
           kpp_1d_fields%difs(ki)=blmc(ki,2)
           kpp_1d_fields%dift(ki)=blmc(ki,3)
        else
           kpp_1d_fields%ghat(ki)=0.
        endif
     enddo

     !Print some info for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
      DO ki=0,4
        print*,'Steve-KPPM3',kpp_1d_fields%difm(ki),kpp_1d_fields%difs(ki) &
           ,kpp_1d_fields%ghat(ki),kbl
      enddo
     endif

!     NPK 25/9/08.  Trap for negative values of diffusivities.
!     If negative, set to a background value of 1E-05.
!     
!     DO 205 ki= 1,km
!     DO 206 i = ipt,ipt
!     IF (dift(i,ki) .LT. 0) dift(i,ki)=1E-05
!     IF (difs(i,ki) .LT. 0) difs(i,ki)=1E-05
!     IF (visc(i,ki) .LT. 0) visc(i,ki)=1E-05
!     206     continue
!     205  continue
     
  ENDIF ! of LKPP

return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_KPPMIX

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_RIMIX (km,kmp1,kpp_1d_fields &
                                              ,kpp_const_fields,ii,jj)
use mem_kpp, only: kpp_1d_type,kpp_const_type
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid
use kpp_parameters, only:kpprnt,iprnt,jprnt

implicit none

  ! compute interior viscosity diffusivity coefficients due to
  ! shear instability (dependent on a local richardson number)
  ! and due to background internal wave activity.
  ! This routine related to Large et al.(1994) Eqns.28 & 29

  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  integer km,kmp1           ! number of vertical levels

  ! local variables
  real Rigg
  real fri,fcon             ! function of Rig
  real ratio
  real epsln,Riinfty,Ricon,difm0,difs0,difmiw,difsiw,difmcon,difscon,c1,c0
  INTEGER ki,mRi,j,ii,jj
      
  data  epsln   / 1.e-16 /  ! a small number          
  data  Riinfty /  0.8     / ! LMD default was = 0.7
  data  Ricon   / -0.2    / ! note: exp was repl by multiplication
  data  difm0   / 0.005  /  ! max visc due to shear instability
  data  difs0   / 0.005  /  ! max diff ..  .. ..    ..
  data  difmiw  / 0.0001  / ! background/internal waves visc(m^2/s)
  data  difsiw  / 0.00001 / ! ..         ..       ..    diff(m^2/s)
  data  difmcon / 0.0000   / ! max visc for convection  (m^2/s)
  data  difscon / 0.0000   / ! max diff for convection  (m^2/s)
  data  c1/ 1.0/
  data  c0/ 0.0/
  data  mRi/ 1 /            ! number of vertical smoothing passes
     
  !     compute interior gradient Ri at all interfaces, except surface     
  !-----------------------------------------------------------------------
  !     compute interior gradient Ri at all interfaces ki=1,km, (not surface)
  !     use visc(imt,ki=1,km) as temporary storage to be smoothed
  !     use dift(imt,ki=1,km) as temporary storage of unsmoothed Ri
  !     use difs(imt,ki=1,km) as dummy in smoothing call
      
  do ki = 1, km
     kpp_1d_fields%Rig(ki)  = kpp_1d_fields%dbloc(ki) &
        * (kpp_const_fields%zm(ki)-kpp_const_fields%zm(ki+1))/&
          (kpp_1d_fields%Shsq(ki) + epsln)
     kpp_1d_fields%dift(ki) = kpp_1d_fields%Rig(ki)
     kpp_1d_fields%difm(ki) = kpp_1d_fields%dift(ki)
  ENDDO

  !-----------------------------------------------------------------------
  !     vertically smooth Ri mRi times
  do j = 1,mRi
     CALL MCKPP_PHYSICS_VERTICALMIXING_Z121 (kmp1,c0,Riinfty &
          ,kpp_1d_fields%difm,kpp_1d_fields%difs)
  enddo

  !-----------------------------------------------------------------------
  !                           after smoothing loop
  DO ki = 1, km
     ! evaluate f of unsmooth Ri (fri) for convection        store in fcon
     ! evaluate f of   smooth Ri (fri) for shear instability store in fri
 
     Rigg  = AMAX1( kpp_1d_fields%dift(ki) , Ricon )
     ratio = AMIN1( (Ricon-Rigg)/Ricon , c1 )
     fcon  = (c1 - ratio*ratio)
     fcon  = fcon * fcon * fcon
     
     Rigg  = AMAX1( kpp_1d_fields%difm(ki) , c0 )
     ratio = AMIN1( Rigg/Riinfty , c1 )
     fri   = (c1 - ratio*ratio)
     fri   = fri * fri * fri
         	     
     ! ************************   Overwrite with Gent's PP **********
     !           fcon = 0.0
     !           Rigg  = AMAX1( dift(i,ki) , c0 )
     !           fri   = c1 / (c1 + 10. * Rigg )
     !           difm0 = 0.1 * fri
     !           difs0 = 0.1 * fri * fri

     !  ************************   Overwrite with original PP
     !           fcon = 0.0
     !           Rigg  = AMAX1( dift(i,ki) , c0 )
     !           fri   = c1 / (c1 +  5. * Rigg )
     !           difm0 = 0.01 * fri
     !           difs0 = (difmiw + fri * difm0)

     ! ----------------------------------------------------------------------
     !            evaluate diffusivities and viscosity
     ! mixing due to internal waves, and shear and static instability
 
     kpp_1d_fields%difm(ki) = (difmiw + fcon * difmcon + fri * difm0)
     kpp_1d_fields%difs(ki) = (difsiw + fcon * difscon + fri * difs0)
     kpp_1d_fields%dift(ki) = kpp_1d_fields%difs(ki)
  END DO

  ! ------------------------------------------------------------------------
  !         set surface values to 0.0

  kpp_1d_fields%difm(0) = c0
  kpp_1d_fields%dift(0) = c0
  kpp_1d_fields%difs(0) = c0

  !Print some info for debugging (Saleeby2018)
  if(KPPRNT==1)then
   DO ki = 0, km
    if(ii+mi0(ngrid)==iprnt.and.jj+mj0(ngrid)==jprnt.and.ki<=5) then
     print*,'Steve-rimix',kpp_1d_fields%Rig(ki+1),kpp_1d_fields%dbloc(ki+1) &
     ,kpp_1d_fields%difm(ki),kpp_1d_fields%difs(ki)
    endif
   ENDDO
  endif

return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_RIMIX

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_DDMIX (km,kmp1,alphaDT,betaDS &
                                              ,kpp_1d_fields)
use mem_kpp, only: kpp_1d_type

implicit none
  
  ! Rho dependent interior flux parameterization.
  ! Add double-diffusion diffusivities to Ri-mix values at blending
  ! interface and below.
 
  ! Necessary for IMPLICIT NONE (NPK 11/2/13)
  integer km,kmp1,ki
  real dsfmax,rrho0
  
  real alphaDT(kmp1)  ! alpha * DT  across interfaces
  real betaDS(kmp1)   ! beta  * DS  across interfaces

  TYPE(kpp_1d_type) :: kpp_1d_fields      

  !local
  real Rrho              ! dd parameter
  real diffdd            ! double diffusion diffusivity scale
  real prandtl           ! prandtl number

  data Rrho0  /  1.9   / ! Rp=(alpha*delT)/(beta*delS)
  data dsfmax / 1.0e-4 / ! .0001 m2/s

  DO ki=1,km           
     ! salt fingering case
     if((alphaDT(ki).gt.betaDS(ki)).and.(betaDS(ki).gt.0.)) then
        Rrho  = MIN(alphaDT(ki) / betaDS(ki) , Rrho0)
        diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
        diffdd     = dsfmax*diffdd*diffdd*diffdd
        kpp_1d_fields%dift(ki) = kpp_1d_fields%dift(ki) + diffdd * 0.8 / Rrho
        kpp_1d_fields%difs(ki) = kpp_1d_fields%difs(ki) + diffdd
        
        ! diffusive convection
     else if ((alphaDT(ki).lt.0.0).and.(betaDS(ki).lt.0.0).and. &
          (alphaDT(ki).lt.betaDS(ki)) ) then
        Rrho    = alphaDT(ki) / betaDS(ki) 
        diffdd  = 1.5e-6*9.0*0.101*exp(4.6*exp(-0.54*(1/Rrho-1)))
        prandtl = 0.15*Rrho
        if (Rrho.gt.0.5) prandtl = (1.85-0.85/Rrho)*Rrho
        kpp_1d_fields%dift(ki) = kpp_1d_fields%dift(ki) + diffdd
        kpp_1d_fields%difs(ki) = kpp_1d_fields%difs(ki) + prandtl*diffdd
        
     endif
  ENDDO
  
return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_DDMIX

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_BLDEPTH (km, kmp1, dVsq, Ritop &
  , ustar, Bo, Bosol, hbl, bfsfc, stable, caseA, kbl, Rib, sigma, wm &
  , ws, kpp_1d_fields, kpp_const_fields,i,j)

use mem_kpp, only:kpp_1d_Type,kpp_const_type
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid,time
use kpp_parameters, only:kpprnt,iprnt,jprnt

implicit none
 
  !     the oceanic planetary boundary layer depth, hbl, is determined as
  !     the shallowest depth where the bulk richardson number is
  !     equal to the critical value, Ricr.
  !     
  !     bulk richardson numbers are evaluated by computing velocity and
  !     buoyancy differences between values at zgrid(kl) < 0 and surface
  !     reference values.
  !     in this configuration, the reference values are equal to the
  !     values in the surface layer.  
  !     when using a very fine vertical grid, these values should be 
  !     computed as the vertical average of velocity and buoyancy from 
  !     the surface down to epsilon*zgrid(kl).
  !
  !     when the bulk richardson number at k exceeds Ricr, hbl is
  !     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
  !
  !     The water column and the surface forcing are diagnosed for 
  !     stable/ustable forcing conditions, and where hbl is relative 
  !     to grid points (caseA), so that conditional branches can be 
  !     avoided in later routines.

  ! Necessary for IMPLICIT NONE
  real bvsq,cekman,cmonob,cs,cv,epsilon,fekman,fmonob,&
       hbf,hekman,hmin,hmin2,hmonob,hri,Ricr,vtc,vtsq,epsln
  integer ka,ksave,ku,kl,i,j

  integer km,kmp1           ! number of vertical levels
  
  ! input
  real dVsq(kmp1)           ! (velocity shear re sfc)^2      (m/s)^2
  real Ritop(km)            ! numerator of bulk Richardson Number (m/s)^2
  !     Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
  real ustar                ! surface friction velocity         (m/s)
  real Bo                   ! surface turbulent buoyancy forcing(m^2/s^3)
  real Bosol                ! radiative buoyancy forcing        (m^2/s^3)
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields
     
  ! output
  real hbl                  ! boundary layer depth              (m)
  real bfsfc                ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
  real stable               ! =1 in stable forcing; =0 unstable
  real caseA                ! =1 in case A, =0 in case B 
  integer kbl               ! index of first grid level below hbl 
       
  ! local
  real Rib(2)               ! Bulk Richardson number
  real sigma                ! normalized depth (d/hbl)
  real wm,ws                ! turbulent velocity scales         (m/s)
  real dmo(2)               ! Monin-Obukhov Depth
  real hek                  ! Ekman depth

  !von karman constant
  real, parameter :: vonk=0.4

  data epsln           /  1.e-16 /
  data Ricr            /  0.30   /
  data epsilon         /  0.1    /
  data cekman          /  0.7    /
  data cmonob          /  1.0    /
  data cs              / 98.96   /
  data cv              /  1.6    /
  data hbf             /  1.0    /     
  
  ! find bulk Richardson number at every grid level find implied hri
  ! find Monin-Obukvov depth at every grid level find L
  ! Compare hri, L and hek to give hbl and kbl
  !   
  ! note: the reference depth is -epsilon/2.*zgrid(k), but the reference
  ! u,v,t,s values are simply the surface layer values,
  ! and not the averaged values from 0 to 2*ref.depth,
  ! which is necessary for very fine grids(top layer < 2m thickness)
  ! note: max values when Ricr never satisfied are
  ! kbl(i)=km and hbl(i) -zgrid(km)
  ! min values are                 kbl(i)=2      hbl(i) -zgrid(1)

  !Terms in Large et al. Eqn.23      
  Vtc =  cv * sqrt(0.2/cs/epsilon) / vonk**2 / Ricr
      
  ! indices for array Rib(i,k), the bulk Richardson number.
  ka = 1
  ku = 2

  ! initialize hbl and kbl to bottomed out values
  Rib(ka) = 0.0
  dmo(ka) = -kpp_const_fields%zm(kmp1)
  kbl = km
  hbl = -kpp_const_fields%zm(km)
  !Large et al. Eqn.24 (Saleeby 2018)
  hek =  cekman * ustar / (abs(kpp_1d_fields%f) + epsln)

  do kl = 2,km
                    
     IF(kbl.ge.km) THEN
        ! use caseA as temporary array for next call to wscale
        caseA = -kpp_const_fields%zm(kl)
            
        ! compute bfsfc= Bo + radiative contribution down to hbf * hbl
        ! bfsfc is Large et al.Eqn.A3d
        bfsfc  = Bo + Bosol * (1. - kpp_1d_fields%swfrac(kl))
        ! if bfsfc>0, stable=1 | if bfsfc<=, stable=0 
        stable = 0.5 + SIGN( 0.5, bfsfc+epsln )
        sigma  = stable * 1. + (1.-stable) * epsilon

        !Print some info here for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
            .and.kl<=20) then
          print*,'Steve-bldepth0',kbl,bfsfc,sigma,stable
        endif

        ! compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
        CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE (sigma,caseA,ustar,bfsfc &
                                             ,wm,ws,kpp_const_fields,i,j)
           
        ! compute the turbulent shear contribution to Rib
        bvsq =0.5*(kpp_1d_fields%dbloc(kl-1) / (kpp_const_fields%zm(kl-1) &
             -kpp_const_fields%zm(kl))+ kpp_1d_fields%dbloc(kl) &
             / (kpp_const_fields%zm(kl)-kpp_const_fields%zm(kl+1)))
        !Large et al. Eqn.23
        Vtsq = -kpp_const_fields%zm(kl) * ws * sqrt(abs(bvsq)) * Vtc
        !compute bulk Richardson number at new level, dunder.
        !Large et al. Eqn.21
        Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsln)
        Rib(ku) = MAX( Rib(ku), Rib(ka) + epsln)
        !     linear interpolate to find hbl where Rib = Ricr
        hri   = -kpp_const_fields%zm(kl-1) + (kpp_const_fields%zm(kl-1) &
                -kpp_const_fields%zm(kl)) * &
             (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))

        !hri can become very large/small (+/- Infinity) depending on layer
        !stability/shear differences, so limit this. Especially cannot have
        !(-Infinity) due to the "MIN" function below for determining mixed
        !layer depth.
        if(hri > 1.e10) hri=1.e10
        if(hri < 0.) hri=1.e10

        !Print some info here for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
            .and.kl<=20) then
          print*,'Steve-bldepth1',ku,Rib(ku),ka,Rib(ka)
        endif

        ! compute the Monin Obukov length scale
        ! Possibly use M-O depth with stable forcing, otherwise
        ! M-O depth is the KPP model bottom
        fmonob    = stable * 1.0
        dmo(ku) = cmonob * ustar * ustar * ustar &
             / vonk / (abs(bfsfc) + epsln)
        dmo(ku) = fmonob * dmo(ku) - (1.-fmonob) * kpp_const_fields%zm(kmp1) 
        if(dmo(ku).le.(-kpp_const_fields%zm(kl))) then
           hmonob =(dmo(ku)-dmo(ka))/(kpp_const_fields%zm(kl-1) &
                   -kpp_const_fields%zm(kl))
           hmonob =(dmo(ku)+hmonob*kpp_const_fields%zm(kl)) / (1.-hmonob)
        else
           hmonob = -kpp_const_fields%zm(kmp1)
        endif
        !Print some info here for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
            .and.kl<=20) then
          print*,'Steve-bldepth2',fmonob,dmo(ku),hmonob
        endif

        ! compute the Ekman depth
        ! Possibly use Ekman depth with stable forcing, otherwise
        ! ekman depth is the KPP model bottom
        fekman  = stable * 1.0
        hekman  = fekman * hek - (1.-fekman) * kpp_const_fields%zm(kmp1)
            
        ! compute boundary layer depth
        hmin  = MIN(hri, hmonob,  hekman, -kpp_1d_fields%ocdepth)

        !Print some info here for debugging (Saleeby2018)
        if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
            .and.kl<=20) then
          print*,'Steve-bldepth3',Rib(ku),hri,hmonob,hekman
        endif

        ! recheck boundary layer depth
        if(hmin .lt. -kpp_const_fields%zm(kl) ) then
     
           ! Code below added by SJW 09/07/04 to solve problems where hek 
           ! less than zgrid(kl-1) giving negative diffusions
           ! if this occurs too often then we need to rethink this fix
     
           ! Modified by NPK 25/09/08 to include Monin-Obukov depth as well.
           ! Richardson depth is sometimes calculated to be huge (> 1E10) when 
           ! Ritop is negative or very small and so is not always helpful in 
           ! this scenario.

           ! This "if" statement used to not be run on initialization step,
           ! but it seems to give better initial values of fluxes when allowed
           ! to run at initialization as well. Keep watch on this (Saleeby2018)
           if (hmin .lt. -kpp_const_fields%zm(kl-1)) then
             hmin2=MIN(hri,hmonob,-kpp_1d_fields%ocdepth)
             if (hmin2 .lt. -kpp_const_fields%zm(kl)) THEN
                hmin=hmin2
             endif
           endif
           
           hbl = hmin
           kbl = kl
        endif
        
     ENDIF !kbl

     !Print some info here for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
         .and.kl<=20) then
      print*,'Steve-bldepth4',kl,hbl,-kpp_const_fields%zm(kl-1) &
        ,-kpp_const_fields%zm(kl)
     endif

     ksave = ka
     ka    = ku
     ku    = ksave     
  ENDDO

  !Saleeby(2018): Do not let model continue if mixed layer depth goes
  !to the ocean model lowest level.
  if(hbl .eq. -kpp_const_fields%zm(km))then
    print*,'HMIX mixed to the bottom: ',hbl,i+mi0(ngrid),j+mj0(ngrid)
    stop
  endif

  !Compute shortwave fraction (bfsfc) at boundary layer depth (hbl)
  CALL MCKPP_PHYSICS_SWFRAC (-1.0,hbl,nint(kpp_1d_fields%jerlov),bfsfc,i,j)

  bfsfc  = Bo + Bosol * (1. - bfsfc)
  stable = 0.5 + SIGN( 0.5, bfsfc)
  bfsfc  = bfsfc + stable * epsln !ensures bfsfc never=0
  
  ! determine caseA and caseB
  caseA  = 0.5 + SIGN( 0.5,-kpp_const_fields%zm(kbl) -0.5 &
         * kpp_const_fields%hm(kbl) -hbl)

  !Print some info here for debugging (Saleeby2018)
  if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1) then
    print*,'Steve-bldepth5',bfsfc,stable,caseA,kbl
  endif

return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_BLDEPTH

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_BLMIX (km, mdiff, ustar, bfsfc, hbl &
    , stable, caseA, kbl,gat1, dat1, dkm1, blmc, sigma, wm, ws &
    , kpp_1d_fields,kpp_const_fields,i,j)

use mem_kpp, only:kpp_1d_type,kpp_const_type
use node_mod, only:mi0,mj0
use mem_grid, only:ngrid
use kpp_parameters, only:kpprnt,iprnt,jprnt

  ! mixing coefficients within boundary layer depend on surface
  ! forcing and the magnitude and gradient of interior mixing below
  ! the boundary layer ("matching").

implicit none

  ! CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
  ! fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)
  
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  integer km                !,kmp1        ! number of vertical levels
  integer mdiff             ! number of viscosities + diffusivities
  
  !    input
  real ustar                ! surface friction velocity         (m/s)
  real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
  real hbl                  ! boundary layer depth              (m)
  real stable               ! = 1 in stable forcing
  real caseA                ! = 1 in case A
  integer kbl               ! index of first grid level below hbl
     
  ! output
  real gat1(mdiff)
  real dat1(mdiff)
  real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
  real blmc(km,mdiff)       ! boundary layer mixing coefficients(m^2/s)
  
  !  local
  real sigma                ! normalized depth (d / hbl)
  real ws, wm               ! turbulent velocity scales         (m/s)
  !  None of these were previously declared ... (NPK 6/2/13)
  real a1,a2,a3,am,as,c1,c2,c3,cg,cm,cs,cstar,delhat,difsh,&
       difsp,difth,diftp,dvdzup,epsln,f1,gm,gs,&
       dvdzdn,epsilon,gt,r,visch,viscp,zetam,sig,&
       zetas
  integer ki,kn,i,j

  !von karman constant
  real, parameter :: vonk=0.4

  data epsln             /   1.e-20 /
  data epsilon           /   0.1    /
  data c1                /   5.0    /
  data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
  data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
  data cstar             /    5.    /     

  cg = cstar * vonk * (cs * vonk * epsilon)**(1./3.)
      
  ! compute velocity scales at hbl
  sigma = stable * 1.0 + (1.-stable) * epsilon
      
  CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE (sigma, hbl, ustar, bfsfc,wm,ws &
     ,kpp_const_fields,i,j)

  kn = ifix(caseA+epsln) *(kbl -1) + (1-ifix(caseA+epsln)) * kbl
      
  ! find the interior viscosities and derivatives at hbl(i) 
  delhat = 0.5*kpp_const_fields%hm(kn)-kpp_const_fields%zm(kn) - hbl
  R      = 1.0 - delhat / kpp_const_fields%hm(kn)
  dvdzup = (kpp_1d_fields%difm(kn-1) - kpp_1d_fields%difm(kn)) &
           / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%difm(kn)   - kpp_1d_fields%difm(kn+1)) &
           / kpp_const_fields%hm(kn+1)
  viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
      
  dvdzup = (kpp_1d_fields%difs(kn-1) - kpp_1d_fields%difs(kn)) &
           / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%difs(kn)   - kpp_1d_fields%difs(kn+1)) &
           / kpp_const_fields%hm(kn+1)
  difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
      
  dvdzup = (kpp_1d_fields%dift(kn-1) - kpp_1d_fields%dift(kn)) &
           / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%dift(kn)   - kpp_1d_fields%dift(kn+1)) &
           / kpp_const_fields%hm(kn+1)
  diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
     
  visch  = kpp_1d_fields%difm(kn) + viscp * delhat
  difsh  = kpp_1d_fields%difs(kn) + difsp * delhat
  difth  = kpp_1d_fields%dift(kn) + diftp * delhat
      
  f1 = stable * c1 * bfsfc / (ustar**4+epsln) 
  gat1(1) = visch / hbl / (wm+epsln)
  dat1(1) = -viscp / (wm+epsln) + f1 * visch
  dat1(1) = min(dat1(1),0.) 
  
  gat1(2) = difsh  / hbl / (ws+epsln)
  dat1(2) = -difsp / (ws+epsln) + f1 * difsh 
  dat1(2) = min(dat1(2),0.) 
  
  gat1(3) = difth /  hbl / (ws+epsln)
  dat1(3) = -diftp / (ws+epsln) + f1 * difth 
  dat1(3) = min(dat1(3),0.) 
      
!     Turn off interior matching here
!     gat1(i,1) = 0.0001
!     gat1(i,2) = 0.00001
!     gat1(i,3) = 0.00001
!     do m=1,3
!       dat1(i,m) = 0.0
!       enddo     

  DO ki = 1,km       
     ! compute turbulent velocity scales on the interfaces
     sig = (-kpp_const_fields%zm(ki) + 0.5 * kpp_const_fields%hm(ki)) / hbl
     sigma   = stable*sig + (1.-stable)*AMIN1(sig,epsilon)

     CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE (sigma, hbl, ustar, bfsfc, wm &
          ,ws,kpp_const_fields,i,j)

     ! compute the dimensionless shape functions at the interfaces     
     sig = (-kpp_const_fields%zm(ki) + 0.5 * kpp_const_fields%hm(ki)) / hbl
     a1 = sig - 2.
     a2 = 3.-2.*sig
     a3 = sig - 1.
   
     Gm = a1 + a2 * gat1(1) + a3 * dat1(1) 
     Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
     Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
   
     !  compute boundary layer diffusivities at the interfaces 
     blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
     blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
     blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)

     ! nonlocal transport term = ghats * <ws>o
     kpp_1d_fields%ghat(ki) = (1.-stable) * cg / (ws*hbl+epsln)

     !Print some info here for debugging (Saleeby2018)
     if(i+mi0(ngrid)==iprnt.and.j+mj0(ngrid)==jprnt.and.KPPRNT==1 &
         .and.ki<=5) then
       print*,'Steve-blmix1',blmc(ki,2),ws,sig,Gs
       print*,'Steve-blmix2',cg,ws,hbl,kpp_1d_fields%ghat(ki)
     endif
  ENDDO

  ! find diffusivities at kbl-1 grid level 
  sig   =  -kpp_const_fields%zm(kbl-1)  / hbl
  sigma =  stable * sig + (1.-stable) * AMIN1(sig,epsilon)
  
  CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE (sigma, hbl, ustar, bfsfc, wm, ws &
       , kpp_const_fields,i,j)

  sig = -kpp_const_fields%zm(kbl-1) / hbl
  a1= sig - 2.
  a2 = 3.-2.*sig
  a3 = sig - 1.
  Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
  Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
  Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
  dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
  dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
  dkm1(3) = hbl * ws * sig * (1. + sig * Gt)
  
return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_BLMIX

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_ENHANCE (km,mdiff,dkm1,hbl,kbl &
                           ,caseA,blmc,kpp_1d_fields,kpp_const_fields)

use mem_kpp, only:kpp_1d_type,kpp_const_type

implicit none

  ! enhance the diffusivity at the kbl-.5 interface

  !     Necessary for IMPLICIT NONE (NPK 11/2/13)
  real dkmp5,dstar
  integer ki
  integer km                !,kmp1           ! number of vertical levels  
  integer mdiff             ! number of viscosities + diffusivities
  integer kbl               ! grid above hbl
  real hbl                  ! boundary layer depth             (m)
  real dkm1(mdiff)          ! bl diffusivity at kbl-1 grid level
  real caseA                ! = 1 in caseA, = 0 in case B
 
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  ! output
  real blmc(km,mdiff)       ! enhanced bound. layer mixing coeff.

  ! local
  real delta                ! fraction hbl lies beteen zgrid neighbors

  do ki=1,km-1         
     if(ki .eq. (kbl - 1) ) then            
        delta = (hbl+kpp_const_fields%zm(ki)) / &
           (kpp_const_fields%zm(ki)-kpp_const_fields%zm(ki+1))
            
        dkmp5 = caseA * kpp_1d_fields%difm(ki) + (1.-caseA) * blmc(ki,1)
        dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5      
        blmc(ki,1) = (1.-delta) * kpp_1d_fields%difm(ki) + delta * dstar
        
        dkmp5 = caseA * kpp_1d_fields%difs(ki) + (1.-caseA) * blmc(ki,2)
        dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5    
        blmc(ki,2) = (1.-delta) * kpp_1d_fields%difs(ki) + delta * dstar
        
        dkmp5 = caseA * kpp_1d_fields%dift(ki) + (1.-caseA) * blmc(ki,3)
        dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5     
        blmc(ki,3) = (1.-delta) * kpp_1d_fields%dift(ki) + delta * dstar
        
        kpp_1d_fields%ghat(ki) = (1.-caseA) * kpp_1d_fields%ghat(ki)            
     endif
  enddo
  
return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_ENHANCE

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_Z121 (kmp1,vlo,vhi,V,w)

implicit none

  INTEGER kmp1
  REAL vlo,vhi,tmp,wait
  INTEGER k,km
  
  ! Apply 121 smoothing in k to 2-d array V(i,k=1,km)
  ! top (0) value is used as a dummy
  ! bottom (kmp1) value is set to input value from above.
 
  real V(0:kmp1)  ! 2-D array to be smoothed in kmp1 direction
  real w(0:kmp1)  ! 2-D array of internal weights to be computed
  
  km  = kmp1 - 1
  
  w(0)    =   0.0         
  w(kmp1) =   0.0             
  V(0)    =   0.0      
  V(kmp1) =   0.0
  
  do k=1,km  
     if((V(k).lt.vlo).or.(V(k).gt.vhi)) then
        w(k) = 0.0
        !     w(i,k) = 1.0
     else 
        w(k) = 1.0
     endif
     
  enddo
  
  do k=1,km
     tmp    = V(k)
     V(k) = w(k-1)*V(0)+2.*V(k)+w(k+1)*V(k+1)
     wait   = w(k-1) + 2.0 + w(k+1)
     V(k) = V(k) / wait             
     V(0) = tmp
  enddo
  
return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_Z121

!##############################################################################
Subroutine MCKPP_PHYSICS_SWFRAC (fact, z, jwtype, swdk, i,j)

implicit none

! compute fraction of solar short-wave flux penetrating to specified
! depth (times fact) due to exponential decay in  Jerlov water type
! reference : two band solar absorption model of simpson and paulson (1977)
      
  integer nwtype,i,j
  parameter(nwtype=5) ! max number of different water types 

  ! input
  real fact      ! scale factor to apply to depth array
  real z         ! vertical height ( <0.) for desired sw 
  !                           fraction (m)
  integer jwtype ! index for jerlov water type

  !  output
  real swdk      !  short wave (radiation) fractional decay

  !  local
  real  rfac(nwtype),a1(nwtype),a2(nwtype)
  real rmin,r1,r2 

  ! jerlov water type :  I       IA      IB      II      III
  !            jwtype    1       2       3       4       5
  data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
  data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
  data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
  data rmin         / -80. /
  
        
  r1   = MAX(z*fact/a1(jwtype), rmin)
  r2   = MAX(z*fact/a2(jwtype), rmin)
  swdk = rfac(jwtype)  * exp(r1) + (1.-rfac(jwtype)) * exp(r2)

return
END SUBROUTINE MCKPP_PHYSICS_SWFRAC

!##############################################################################
Subroutine MCKPP_PHYSICS_VERTICALMIXING_WSCALE (sigma, hbl, ustar, bfsfc, wm &
                                              , ws, kpp_const_fields,i,j)
use mem_kpp, only: kpp_const_type

implicit none

  ! compute turbulent velocity scales.
  ! use a 2D-lookup table for wm and ws as functions of ustar and
  ! zetahat (=vonk*sigma*hbl*bfsfc).

  ! Necessary for IMPLICIT NONE (NPK 11/2/13)
  INTEGER ni,nj,iz,izp1,ju,jup1,i,j
  REAL am,as,c1,c2,c3,cm,cs,epsln,fzfrac,ucube,udiff,ufrac,&
       usta,wam,was,wbm,wbs,zdiff,zetas,zfrac,zetam
  
  TYPE(kpp_const_type) :: kpp_const_fields

  ! lookup table
  parameter ( ni = 890,&   ! number of values for zehat
       nj = 48)            ! number of values for ustar

  real deltaz               ! delta zehat in table
  real deltau               ! delta ustar in table
  real zmin,zmax            ! zehat limits for table
  real umin,umax            ! ustar limits for table

  data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
  data umin,umax  /  0.   , .04   / ! m/s
      
  ! input
  real sigma                ! normalized depth (d/hbl)
  real hbl                  ! boundary layer depth (m)
  real ustar                ! surface friction velocity         (m/s)
  real bfsfc                ! total surface buoyancy flux       (m^2/s^3)
      
  ! output
  real wm,ws                ! turbulent velocity scales at sigma
      
  ! local
  real zehat                ! = zeta *  ustar**3
  real zeta                 ! = stability parameter d/L

  !von karman constant
  real, parameter :: vonk=0.4

  data epsln           /   1.0e-20/
  data c1              /   5.0   /
  data am,cm,c2,zetam  /   1.257 ,  8.380 , 16.0 , - 0.2  /
  data as,cs,c3,zetas  / -28.86  , 98.96  , 16.0 , - 1.0  /

  deltaz = (zmax-zmin)/(ni+1) 
  deltau = (umax-umin)/(nj+1)
  
  ! use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
  zehat = vonk * sigma * hbl * bfsfc
  
  IF (zehat .le. zmax) THEN
     zdiff  = zehat-zmin
     iz = int( zdiff/deltaz )
     iz = min( iz , ni )
     iz = max( iz , 0  )
     izp1=iz+1
     
     udiff  = ustar-umin
     ju = int( udiff/deltau)
     ju = min( ju , nj )
     ju = max( ju , 0  )
     jup1=ju+1
     
     zfrac = zdiff/deltaz - float(iz)
     ufrac = udiff/deltau - float(ju)
     
     fzfrac= 1.-zfrac
     wam   = (fzfrac)  * kpp_const_fields%wmt(iz,jup1) + &
          zfrac*kpp_const_fields%wmt(izp1,jup1)
     wbm   = (fzfrac)  * kpp_const_fields%wmt(iz,ju  ) + &
          zfrac*kpp_const_fields%wmt(izp1,ju  )
     wm    = (1.-ufrac)* wbm          + ufrac*wam
         
     was   = (fzfrac)  * kpp_const_fields%wst(iz,jup1) + &
          zfrac*kpp_const_fields%wst(izp1,jup1)
     wbs   = (fzfrac)  * kpp_const_fields%wst(iz,ju  ) + &
          zfrac*kpp_const_fields%wst(izp1,ju  )
     ws    = (1.-ufrac)* wbs          + ufrac*was       
  ELSE
     ucube = ustar**3
     wm = vonk * ustar * ucube / (ucube + c1 * zehat)
     ws = wm     
  ENDIF

return
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_WSCALE

