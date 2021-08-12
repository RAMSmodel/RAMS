!##############################################################################
Subroutine MCKPP_PHYSICS_OCNINT (kpp_1d_fields,kpp_const_fields &
                                ,kmixe,Uo,Xo,ii,jj)

use mem_kpp, only:kpp_1d_type,kpp_const_type
use kpp_parameters, only:kppNZ,kppNZP1,kppNVEL,kppNSCLR

implicit none

  ! Integrate the ocn model by backwards Euler(implicit)discretization
  ! On input : Un,Xn are estimated profiles which are used
  !            to estimate diffusivity profiles at new time.
  !          : Updated diffusivities from Un Xn are in common
  ! On output: Un,Xn are new profiles after integration.

  ! Written  19 March 1991 - jan

  ! Input
  REAL Uo(kppNZP1,kppNVEL),Xo(kppNZP1,kppNSCLR)
  
  ! Output
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  ! Local
  real cu (kppNZ), &! upper coeff for (k-1) on k line of trid.matrix
       cc (kppNZ), &! central ...     (k  ) ..
       cl (kppNZ), &! lower .....     (k-1) ..
       rhs(kppNZ)   ! right-hand-side terms
  real diff(0:kppNZ),gcap(kppNZ),ntflx(0:kppNZ,kppNSCLR)
  ! More local variables to make implicit none
  integer kmixe,i,imode,n,k,ii,jj
  real ftemp,ghatflux,sturflux
  integer adv_mode
  real adv_mag

  ! ********************************************************************
  ! U and V solution of tridiagonal matrix
  !                set f = 0 for equatorial application
  ftemp = kpp_1d_fields%f

  DO k=0,kppNZ
     diff(k)=kpp_1d_fields%difm(k)
  ENDDO
  CALL MCKPP_PHYSICS_SOLVERS_TRIDCOF (diff,kppNZ,cu,cc,cl,kpp_const_fields)

  ! U right hand side and solution
  rhs(1)= Uo(1,1) + kpp_const_fields%dto* (ftemp*.5*(Uo(1,2) &
      +kpp_1d_fields%U(1,2)) - &
       kpp_1d_fields%wU(0,1)/kpp_const_fields%hm(1)) 
  do i=2,kppNZ-1
     rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2) &
      +kpp_1d_fields%U(i,2))      
  enddo
  i=kppNZ                      ! bottom
  rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2) &
       + kpp_1d_fields%U(i,2)) &
       + kpp_const_fields%tri(i,1)*kpp_1d_fields%difm(i)*Uo(i+1,1)

  CALL MCKPP_PHYSICS_SOLVERS_TRIDMAT (cu,cc,cl,rhs,Uo(:,1),kppNZ &
                                     ,kpp_1d_fields%U(:,1))

  ! V rhs and solution
  rhs(1)= Uo(1,2) - kpp_const_fields%dto* (ftemp*.5*(Uo(1,1) &
       +kpp_1d_fields%U(1,1))+kpp_1d_fields%wU(0,2)/kpp_const_fields%hm(1))
  do i=2,kppNZ-1
     rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1) &
            +kpp_1d_fields%U(i,1)) 
  enddo
  i=kppNZ
  rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1) &
       +kpp_1d_fields%U(i,1)) &
       +kpp_const_fields%tri(i,1)*kpp_1d_fields%difm(i)*Uo(i+1,2)

  CALL MCKPP_PHYSICS_SOLVERS_TRIDMAT (cu,cc,cl,rhs,Uo(:,2),kppNZ &
                                     ,kpp_1d_fields%U(:,2))

! *******************************************************************
! Scalar solutions of tridiagonal matrix
!     Temperature (different from other scalars because of ghat-term
!                  and double diffusion)
!     ghatflux = wX(0,1) - (1-SWDK(-hmixe,real(time)))
!                         * sflux(3) / rho(ipt,0) / CP(ipt,0)
!     ghatflux = wX(0,1) - (1-SWDK(-d(1) ,real(time)))
!                         * sflux(3) / rho(ipt,0) / CP(ipt,0)

  ! **********************************************************************
  ! Ocean Temperature
  ! **********************************************************************
  ghatflux = kpp_1d_fields%wX(0,1)
  sturflux = kpp_1d_fields%wX(0,1)
  diff(0)=kpp_1d_fields%dift(0)
  ntflx(0,1)=kpp_1d_fields%wXNT(0)
  DO k=1,kppNZ
     diff(k)=kpp_1d_fields%dift(k)
     gcap(k)=kpp_1d_fields%ghat(k)
     ntflx(k,1)=kpp_1d_fields%wXNT(k)
  ENDDO
  CALL MCKPP_PHYSICS_SOLVERS_TRIDCOF (diff,kppNZ,cu,cc,cl,kpp_const_fields)
  
  CALL MCKPP_PHYSICS_SOLVERS_TRIDRHS (Xo(:,1),ntflx(:,1) &
  ,diff,gcap,sturflux,ghatflux,kppNZ,rhs,kpp_const_fields)  

  ! Sfc relaxation incompatible with flux corrections at depth (NPK 12/02/08)
  IF (kpp_const_fields%relax_sst .GT. 1.e-10) THEN     
     ! Relax the Mixed layer temperature back to SST0
     ! By using a flux correction at the surface
     ! Added by SJW (06/04)
     ! dm(kmixe) = hmix = mixed layer depth
     ! hm(1) = depth of top layer
      rhs(1)=rhs(1)+kpp_const_fields%dto*kpp_const_fields%relax_sst*&
            (kpp_1d_fields%SST0-Xo(1,1))*kpp_const_fields%dm(kmixe) &
            /kpp_const_fields%hm(1)
  ENDIF
      
  kpp_1d_fields%tinc_fcorr(:)=0.

  ! Relax the temperature at each layer in the model by computing a flux 
  ! correction at each layer.  Requires a three-dimensional (x,y,z) data
  ! of ocean temperatures
  IF (kpp_const_fields%relax_ocnT .GT. 1.e-10) THEN
     DO k=1,kppNZP1
        ! Store the relaxation term as tinc_fcorr so that, on output,
        ! that field contains the actual correction applied in K/timestep.
        kpp_1d_fields%tinc_fcorr(k)=kpp_1d_fields%tinc_fcorr(k)+&
             kpp_const_fields%dto*kpp_const_fields%relax_ocnT*&
             (kpp_1d_fields%ocnT_clim(k)-Xo(k,1))
     ENDDO
  ENDIF

  !The following do loop is to compute "ocntcorr" when performing
  !interactive relaxation (relax_ocnt). It is a flux of heat in W/m3.
  !This diagnostic is provided so the user can derive the heat
  !fluxes necessary to nudge the model when NOT doing interactive 
  !relaxation, but rather, a flux correction. This is done in multi-year
  !climate simulations. We have removed the flux correction method from
  !this version in RAMS, but could reinstall if needed. Commenting out
  !"ocntcorr" below for now.
  !In climate mode, as is done for KPP in CAM3, they use a 2-step process
  !whereby they first run the model with "interactive relaxation", then
  !compute a mean seasonal cycle of the nudging terms, then secondly, impose
  !thoses flux corrections in a second simulation.
  DO k=1,kppNZP1
     rhs(k) = rhs(k) + kpp_1d_fields%tinc_fcorr(k)
     ! Modify the correction field so that, when output, it is in
     ! the correct units to be input as a flux correction         
  !   kpp_1d_fields%ocnTcorr(k)=kpp_1d_fields%tinc_fcorr(k) &
  !       *kpp_1d_fields%rho(k) &
  !       *kpp_1d_fields%cp(k)/kpp_const_fields%dto
  ENDDO

  CALL MCKPP_PHYSICS_SOLVERS_TRIDMAT (cu,cc,cl,rhs,Xo(:,1) &
                                     ,kppNZ,kpp_1d_fields%X(:,1))

  ! **********************************************************************
  ! Salinity
  ! **********************************************************************
  ghatflux = kpp_1d_fields%wX(0,2) 
  sturflux = kpp_1d_fields%wX(0,2)
  diff(0)=kpp_1d_fields%difs(0)
  ntflx(0,2)=0.0
  DO k=1,kppNZ
     diff(k)=kpp_1d_fields%difs(k)
     gcap(k)=kpp_1d_fields%ghat(k)
     ntflx(k,2)=0.0
  ENDDO
  CALL MCKPP_PHYSICS_SOLVERS_TRIDCOF (diff,kppNZ,cu,cc,cl,kpp_const_fields)

  CALL MCKPP_PHYSICS_SOLVERS_TRIDRHS (Xo(:,2),ntflx(:,2) &
  ,diff,gcap,sturflux,ghatflux,kppNZ,rhs,kpp_const_fields)

  kpp_1d_fields%sinc_fcorr(:)=0.

  ! Relax the salinity at each layer in the model by computing a flux 
  ! correction at each layer.  Requires a three-dimensional(x,y,z) input 
  ! file of salinity
  IF (kpp_const_fields%relax_sal .GT. 1.e-10) THEN
    DO k=1,kppNZP1
       ! Store the relaxation term as sinc_fcorr so that, on output,
       ! that field contains the actual correction applied in psu/timestep
       kpp_1d_fields%sinc_fcorr(k) = kpp_1d_fields%sinc_fcorr(k)+&
            kpp_const_fields%dto*kpp_const_fields%relax_sal* &
            (kpp_1d_fields%sal_clim(k)-Xo(k,2))
    ENDDO
  ENDIF

  !The following do loop is to compute "scorr" when performing            
  !interactive relaxation (relax_sal). It is a flux of salt in W/m3.
  !This diagnostic is provided so the user can derive the salt         
  !fluxes necessary to nudge the model when NOT doing interactive 
  !relaxation, but rather, a flux correction. This is done in multi-year
  !climate simulations. We have removed the flux correction method from
  !this version in RAMS, but could reinstall if needed. Commenting out
  !"scorr" below for now.
  !In climate mode, as is done for KPP in CAM3, they use a 2-step process
  !whereby they first run the model with "interactive relaxation", then
  !compute a mean seasonal cycle of the nudging terms, then secondly, impose
  !thoses flux corrections in a second simulation.
  DO k=1,kppNZP1
     rhs(k) = rhs(k) + kpp_1d_fields%sinc_fcorr(k)
     ! Modify the correction field so that, when output, it is in
     ! the correct units to be input as a flux correction
  !   kpp_1d_fields%scorr(k)=kpp_1d_fields%sinc_fcorr(k) &
  !                         /kpp_const_fields%dto
  ENDDO

  CALL MCKPP_PHYSICS_SOLVERS_TRIDMAT (cu,cc,cl,rhs,Xo(:,2) &
                                     ,kppNZ,kpp_1d_fields%X(:,2))

return
END SUBROUTINE MCKPP_PHYSICS_OCNINT

!##############################################################################
Subroutine MCKPP_PHYSICS_SOLVERS_TRIDCOF (diff,nzi,cu,cc,cl,kpp_const_fields)

use mem_kpp, only:kpp_const_type

  ! Compute coefficients for tridiagonal matrix (dimension=nzi).
  !     Note: cu(1) = 0. and cl(nzi) = 0. are necessary conditions.
  
  IMPLICIT NONE

  ! Input
  TYPE(kpp_const_type) :: kpp_const_fields
  integer nzi             ! dimension of field
  real diff(0:nzi) ! diffusivity profile on interfaces
  
  ! Output
  real cu(nzi),&    ! upper coeff. for (k-1) on k line of trid.matrix
       cc(nzi),&    ! central ...      (k  ) ..
       cl(nzi)      ! lower .....      (k-1) ..
  integer i

  ! In the surface layer
  cu(1) = 0.
  cc(1) = 1. + kpp_const_fields%tri(1,1)*diff(1) ! 1.+ dto/h(1)/dzb(1)*diff(1)
  cl(1) =    - kpp_const_fields%tri(1,1)*diff(1) !   - dto/h(1)/dzb(1)*diff(1)

  ! Inside the domain
  do i=2,nzi
     cu(i) =    - kpp_const_fields%tri(i,0)*diff(i-1)
     cc(i) = 1. + kpp_const_fields%tri(i,1)*diff(i) &
                + kpp_const_fields%tri(i,0)*diff(i-1)
     cl(i) =    - kpp_const_fields%tri(i,1)*diff(i)
  ENDDO

  !In the bottom layer
  cl(nzi)= 0.

return
END SUBROUTINE MCKPP_PHYSICS_SOLVERS_TRIDCOF

!##############################################################################
Subroutine MCKPP_PHYSICS_SOLVERS_TRIDRHS (yo,ntflux,diff,ghat,sturflux &
    ,ghatflux,nzi,rhs,kpp_const_fields)

use mem_kpp, only:kpp_const_type

  ! Compute right hand side of tridiagonal matrix for scalar fields:
  !  =  yo (old field) 
  !     + flux-divergence of ghat
  !     + flux-divergence of non-turbulant fluxes
  ! Note: surface layer needs +dto/h(1) * surfaceflux
  ! bottom  ..... ..... +dto/h(nzi)*diff(nzi)/dzb(nzi)*yo(nzi+1)

  IMPLICIT NONE
  
  !  Input
  TYPE(kpp_const_type) :: kpp_const_fields

  integer nzi         ! dimension of field
  real yo(nzi+1),&    ! old profile
       ntflux(0:nzi),&! non-turbulent flux = wXNT(0:nzi)
       diff(0:nzi),&  ! diffusivity profile on interfaces
       ghat(nzi),&    ! ghat turbulent flux   
       sturflux,&     ! surface turbulent (kinematic) flux = wX(0,n)
       ghatflux       ! surface flux for ghat: includes solar flux      

  ! Output
  real rhs(nzi)      ! right hand side

  ! more local variables to make implicit none
  integer i
  
  ! In the surface layer (dto/h(1)=tri(0,1)
  rhs(1)= yo(1) + kpp_const_fields%dto/kpp_const_fields%hm(1) &
    * (ghatflux*diff(1)*ghat(1)- sturflux + &
       ntflux(1) - ntflux( 0 ) )

  ! Inside the rest of the domain
  do i=2,nzi-1
     rhs(i)= yo(i) + kpp_const_fields%dto/kpp_const_fields%hm(i) &
         * ( ghatflux*(diff(i)*ghat(i) &
           - diff(i-1)*ghat(i-1)) + ntflux(i) - ntflux(i-1) )
  ENDDO
      
  !     In the bottom layer     
  if(nzi.gt.1) then   ! not for slab ocean
     i=nzi
     rhs(i)= yo(i) + kpp_const_fields%dto/kpp_const_fields%hm(i) &
           * ( ghatflux*(diff(i)*ghat(i) &
             - diff(i-1)*ghat(i-1)) + ntflux(i) - ntflux(i-1) ) &
             + yo(i+1)*kpp_const_fields%tri(i,1)*diff(i)
  endif
  
return
END SUBROUTINE MCKPP_PHYSICS_SOLVERS_TRIDRHS

!##############################################################################
Subroutine MCKPP_PHYSICS_SOLVERS_TRIDMAT (cu,cc,cl,rhs,yo,nzi,yn)

  ! Solve tridiagonal matrix for new vector yn, given right hand side
  ! vector rhs. Note: yn(nzi+1) = yo(nzi+1).
  IMPLICIT NONE

  ! Input
  integer nzi               ! dimension of matrix
  real cu (nzi),&            ! upper coeff. for (k-1) on k line of tridmatrix
       cc (nzi),&            ! central ...      (k  ) ..
       cl (nzi),&            ! lower .....      (k-1) ..
       rhs(nzi),&           ! right hand side
       yo(nzi+1),yni            ! old field

  ! Output
  real yn(nzi+1)    ! new field
  
  ! Local 
  real gam(nzi),&    ! temporary array for tridiagonal solver
       bet           ! ...
  ! more local for implicit none
  integer i

  ! Solve tridiagonal matrix.
  bet   = cc(1)
  yn(1) =  rhs(1) / bet    ! surface
  DO i=2,nzi
     gam(i)= cl(i-1)/bet
     bet   = cc(i) - cu(i)*gam(i)
     if(bet.eq.0.) then
        write(*,*)'* algorithm for solving tridiag matrix fails'
        write(*,*)'* bet=',bet
        write(*,*)'*i-1=',i-1,' cc=',cc(i-1),'cl=',cl(i-1)
        write(*,*)'*i=',i,' cc=',cc(i),' cu=',cu(i),' gam=',gam(i)
        stop
        bet=1.0E-12
        !     Pause 3
     endif
     ! to avoid "Underflow" at single precision on the sun
     yn(i) =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet     
  ENDDO

  do i=nzi-1,1,-1
     yn(i)  = yn(i) - gam(i+1)*yn(i+1)
  ENDDO
  yn(nzi+1) = yo(nzi+1)
  
return
END SUBROUTINE MCKPP_PHYSICS_SOLVERS_TRIDMAT
