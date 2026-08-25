! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Provides cloud optical properties as a function of effective radius for the RRTMGP bands
!   Based on Mie calculations for liquid
!     and results from doi:10.1175/JAS-D-12-039.1 for ice with variable surface roughness
!   Can use either look-up tables or Pade approximates according to which data has been loaded
!   Mike Iacono (AER) is the original author
!
! The class can be used as-is but is also intended as an example of how to extend the RTE framework
! -------------------------------------------------------------------------------------------------

module mo_aerosol_optics
  use mo_rte_kind,      only: wp, wl
  use mo_rte_config,    only: check_values, check_extents
  use mo_rte_util_array,only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  implicit none
  private
  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_aerosol_optics
    private
    ! Lookup table information
    !
    ! Table size information
    integer :: naerbins=17
    integer :: nrh=20

    !Drop radius at bin center in centimeters
    real(wp),dimension(:), allocatable :: rma, radlim 
    real(wp),dimension(:,:), allocatable :: radparam


    ! The tables themselves.
    !
    real(wp), dimension(:,:,:,:), allocatable :: lut_qext, lut_qscat, lut_gasym
    real(wp), dimension(:,:,:), allocatable :: lut_qext_carb, lut_qscat_carb, lut_gasym_carb
    real(wp), dimension(:,:,:), allocatable :: lut_qext_dust, lut_qscat_dust, lut_gasym_dust
    real(wp), dimension(:,:,:), allocatable :: lut_growth

  contains
    generic,   public :: load  => load_lut
    ! Internal procedures
    procedure, private :: load_lut
    procedure, public :: aerosol_optics
  end type ty_aerosol_optics

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routines to load data needed for aerosol optics calculations.
  !
  ! ------------------------------------------------------------------------------
  function load_lut(this, band_lims_wvn, rma, radlim, radparam, &
                    lut_qext, lut_qext_carb, lut_qext_dust,    &
                    lut_qscat, lut_qscat_carb, lut_qscat_dust, &
                    lut_gasym, lut_gasym_carb, lut_gasym_dust, &
                    lut_growth) result(error_msg)

    class(ty_aerosol_optics),     intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn,radparam
    real(wp), dimension(:), intent(in   ) :: rma, radlim

    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter 
    real(wp), dimension(:,:,:,:), intent(in) :: lut_qext, lut_qscat, lut_gasym
    real(wp), dimension(:,:,:), intent(in) :: lut_qext_carb, lut_qscat_carb, lut_gasym_carb
    real(wp), dimension(:,:,:), intent(in) :: lut_qext_dust, lut_qscat_dust, lut_gasym_dust
    real(wp), dimension(:,:,:), intent(in) :: lut_growth
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer :: nbnd, ntype, nradbins, nrh, ncarb, ndust

    error_msg = this%init(band_lims_wvn, name="RRTMGP aerosol optics")
    !
    ! LUT coefficient dimensions
    !
    nbnd      = size(lut_qext,dim=2)
    ntype     = size(lut_qext,dim=1)
    nradbins  = size(lut_qext,dim=3)
    nrh       = size(lut_qext,dim=4)
    ncarb     = size(lut_qext_carb,dim=1)
    ndust     = size(lut_qext_dust,dim=1)

    allocate(this%rma(nradbins))
    this%rma = rma

    allocate(this%radlim(size(radlim)))
    this%radlim = radlim

    allocate(this%radparam(size(radlim),size(rma)))
    this%radparam = radparam

    ! Allocate LUT coefficients
    allocate(this%lut_qext(ntype, nbnd, nradbins,nrh), &
             this%lut_qext_carb(ncarb, nbnd, nradbins), &
             this%lut_qext_dust(ndust, nbnd, nradbins), &
             this%lut_qscat(ntype, nbnd, nradbins,nrh), &
             this%lut_qscat_carb(ncarb, nbnd, nradbins), &
             this%lut_qscat_dust(ndust, nbnd,nradbins), &
             this%lut_gasym(ntype, nbnd, nradbins, nrh), &
             this%lut_gasym_carb(ncarb, nbnd, nradbins), &
             this%lut_gasym_dust(ndust, nbnd,nradbins), &
             this%lut_growth(ntype, nradbins, nrh))


    ! Load LUT coefficients
    this%lut_qext = lut_qext
    this%lut_qext_carb = lut_qext_carb
    this%lut_qext_dust = lut_qext_dust
    this%lut_qscat = lut_qscat
    this%lut_qscat_carb = lut_qscat_carb
    this%lut_qscat_dust = lut_qscat_dust
    this%lut_gasym = lut_gasym
    this%lut_gasym_carb = lut_gasym_carb
    this%lut_gasym_dust = lut_gasym_dust
    this%lut_growth = lut_growth
  end function load_lut
  ! ------------------------------------------------------------------------------
  !
  ! Derive aerosol optical properties from provided aerosol physical properties
  !
  ! ------------------------------------------------------------------------------
  !
  ! Compute single-scattering properties
  !
  function aerosol_optics(this, &
                        aerocon, aeroradius, aerotype, rh_lay, &
                        optical_props, aodt) result(error_msg)
    class(ty_aerosol_optics), &
              intent(in   ) :: this
    real(wp), intent(in   ) :: aeroradius  (:,:,:), &    ! aerosol radius (m)
                               aerocon (:,:,:), &     ! aerosol concentration (#/m^2)
                               rh_lay(:,:)         ! RH
    integer,  intent(in   ) :: aerotype(:,:)
    real(wp), intent(inout), optional :: aodt(:)
    class(ty_optical_props_arry), &
              intent(inout) :: optical_props

    character(len=128)      :: error_msg
    ! ------- Local -------
                ! Optical properties: tau, tau*ssa, tau*ssa*g
    real(wp),    dimension(size(aerocon,1), size(aerocon,2), this%get_nband()) :: &
                tau, taussa, taussag
    real(wp), dimension(size(aerocon,2),this%naerbins) :: rnaer
    real(wp) :: tau_sum, om_sum, asym_sum, con_sum, om, gg
    real(wp) :: gfact, qext, qscat, gasym, dm, pi=3.141593
 
    integer  :: ncol, nlay, nbnd, aerocat, nradbins, nrh 
    integer  :: icol, ilay, ibnd, acat, atype, irh, ibns
    ! scalars for total tau, tau*ssa
    real(wp) :: tau_temp, taussa_temp
    ! ----------------------------------------
    !
    ! Error checking
    !
    ! ----------------------------------------

    error_msg = ''
    if(.not.(allocated(this%lut_qext))) then
      error_msg = 'aerosol optics: no data has been initialized'
      return
    end if

    ncol = size(aerocon,1)
    nlay = size(aerocon,2)
    aerocat = size(aerocon,3)
    nradbins  = size(this%lut_qext,dim=3)
    nrh       = size(this%lut_qext,dim=4)
    nbnd = this%get_nband()

    tau = 0.
    taussa = 0.
    taussag = 0.

    !
    ! Spectral consistency
    !
    if(check_values) then
      if(.not. this%bands_are_equal(optical_props)) &
        error_msg = "aerosol optics: optical properties don't have the same band structure"
      if(optical_props%get_nband() /= optical_props%get_ngpt() ) &
        error_msg = "aerosol optics: optical properties must be requested by band not g-points"
      if(error_msg /= "") return
    end if

    error_msg = ""
    if(error_msg == "") then
     do icol=1,ncol
      DO acat=1,aerocat
        IF (aerotype(icol,acat).gt.0) then
          atype = aerotype(icol,acat)
          ! For each aerosol, the code requires a total number concentration of
          ! the aerosol at each grid point, along with either a total mass at the
          ! grid point or a median radius.
          CALL aero_bin (nlay,nradbins,aerocon(icol,:,acat),aeroradius(icol,:,acat), &
                  this%radlim, this%radparam, rnaer)
!print*, aerocon(icol,10,acat),aeroradius(icol,10,acat),aerocon(icol,10,acat)/sum(rnaer(10,:))
          ! This is the routine that calculates the optical properties of the
          ! aerosols within the radiation routine:
          ! Doing over all vertical levels (not including radiation levels)
          DO ilay = 1,nlay
            ! Locate RH as a percentage for table lookup
            irh=INT(100.*(rh_lay(icol,ilay)-0.80)) + 1
            irh=MAX(1,MIN(irh,nrh) )
       
            ! Loop over radiation bands
            DO ibnd=1,nbnd
  
              ! Resetting temporary storage variabes to zero for next set of calcs:
              tau_sum = 0.0
              om_sum = 0.0
              asym_sum = 0.0
              con_sum = 0.0
              om  = 0.0
              gg  = 0.0
  
              ! Calculating optical properties at each aerosol bin:
              DO ibns = 1,this%naerbins
                ! Retrieving values for the growth factor, the extinction coefficient,
                ! the scattering coefficient and the asymmetry parameter, from their
                ! appropriate look up tables:
                if (atype<=2) then
                  gfact = this%lut_growth(atype,ibns,irh)
                  qext = this%lut_qext(atype,ibnd,ibns,irh)
                  qscat = this%lut_qscat(atype,ibnd,ibns,irh)
                  gasym = this%lut_gasym(atype,ibnd,ibns,irh)
                elseif (atype==4 .or. atype==5) then
                  gfact = 1.0
                  qext = this%lut_qext_carb(atype-3,ibnd,ibns)
                  qscat = this%lut_qscat_carb(atype-3,ibnd,ibns)
                  gasym = this%lut_gasym_carb(atype-3,ibnd,ibns)
                elseif (atype>=31) then
                  gfact = 1.0
                  qext = this%lut_qext_dust(atype-30,ibnd,ibns)
                  qscat = this%lut_qscat_dust(atype-30,ibnd,ibns)
                  gasym = this%lut_gasym_dust(atype-30,ibnd,ibns)
                endif
  
                ! Determining mean diameter of bin (meters), with deliquescence growth factor
                dm = 2. * this%rma(ibns) * gfact
                ! Updating temporary storage variables:
                if (qext>0.0 .and. rnaer(ilay,ibns)>0.0)then
                  tau_sum  = tau_sum  + (pi/4. * dm**2 * rnaer(ilay,ibns) * qext)  !Bext*dz
                  om_sum   = om_sum   + (pi/4. * dm**2 * rnaer(ilay,ibns) * qscat) !Bscat*dz
                  asym_sum = asym_sum + (rnaer(ilay,ibns) * gasym) !Asymmetry parameter
                  con_sum  = con_sum  + rnaer(ilay,ibns) !Total aerosol number
                endif
              ENDDO ! Aerosol bins (ibns)

              if(tau_sum>0.0) then
                om   = om_sum / tau_sum  !Omega = Bscat / Bext
                gg   = asym_sum/(con_sum+1.E-30) !normalize asym_sum by total number
              else
                tau_sum = 0.0
              endif
              tau(icol,ilay,ibnd)   = tau(icol,ilay,ibnd)    + tau_sum
              taussa(icol,ilay,ibnd) = taussa(icol,ilay,ibnd)  + om * tau_sum
              taussag(icol,ilay,ibnd)   = taussag(icol,ilay,ibnd)    + gg * om * tau_sum
            ENDDO ! Radiation band loop
          ENDDO ! Vertical level loop
        ENDIF !Aerosol type if
      ENDDO ! Aerotype loop
     enddo

      select type(optical_props)
        type is (ty_optical_props_1scl)
          do ibnd = 1, nbnd
            do ilay = 1, nlay
              do icol = 1,ncol
                ! Absorption optical depth  = (1-ssa) * tau = tau - taussa
                optical_props%tau(icol,ilay,ibnd) = tau(icol,ilay,ibnd) - taussa(icol,ilay,ibnd)
              end do
            end do
          end do
        type is (ty_optical_props_2str)
          do ibnd = 1, nbnd
            do ilay = 1, nlay
              do icol = 1,ncol
                tau_temp    = tau(icol,ilay,ibnd)
                taussa_temp = taussa(icol,ilay,ibnd)
                optical_props%g  (icol,ilay,ibnd) = taussag(icol,ilay,ibnd) / &
                                                      max(epsilon(tau_temp), taussa_temp)
                optical_props%ssa(icol,ilay,ibnd) = taussa_temp/max(epsilon(tau_temp), tau_temp)
                optical_props%tau(icol,ilay,ibnd) = tau_temp
              !AOD typically reported at ~550nm (MODIS). This is 19,000 cm^-1
                if (present(aodt)) then
                  if (this%band_lims_wvn(1,ibnd)<19000 .and. this%band_lims_wvn(2,ibnd)>19000) then
                    aodt(icol)=aodt(icol)+tau_temp
                  endif
                endif
              end do
            end do
          end do
        type is (ty_optical_props_nstr)
          error_msg = "aerosol optics: n-stream calculations not yet supported"
      end select
    else
      error_msg = "aerosol optics: no method to calculate cloud optical properties"
    endif

  end function aerosol_optics
!##############################################################################
Subroutine aero_bin (m1,naerbins,numconc,medianrad,radlim,radparam,rnaer)

! The purpose of this routine was to take a number concentration
! and median radius for any given aerosol and turn it into a binned
! distribution
! Input Variables: m1 - Integer, number of vertical levels
!    naerbins  - Integer number of bins used in aerosol radiation
!                calculations.  Set to 17 in radcalc3.
!    numconc   - Number concentration of aerosol (#/m^3)
!    medianrad - Median radius of the given aerosol type (m)
!    interp    - Determines beginning point for interpolation calculations
!    radlim    - Radii at which aerosol distribution was explicitly
!                  binned (in cm)
!    radparam  - Percentage of aerosol mass in each of the 17 bins,
!                as a func of median radius.
!    rg        - Median radius of aerosol (m)
!    L(4)      - Values used in the 3rd order Lagrangian polynomial
!                interpolation scheme to determine percentage of aerosol
!                within a given bin
!    m         - Used to define which four median radii should
!                 be used in interpolation
! Output variable: rnaer(m1,naerbins) - Number concentrations for
!                     all 17 bins at every vertical level (#/m3)
!**************************************************************************

implicit none

! Incoming variable declaration:
        INTEGER m1,naerbins
        REAL numconc(m1),medianrad(m1),dn0(m1)
        REAL radlim(:),radparam(:,:)

! Internal variable declaraion:
        INTEGER bin,i,k,n,m,interp
        REAL rg,L(4)

! Outgoing variable declaration:
        REAL rnaer(m1,naerbins)
! Vertical level loop:
  DO k = 1,m1

          rg = medianrad(k)

          !Determining which four points to use in the interpolation scheme
          interp = 0
          DO n=1,25
            IF(rg.GE.radlim(n).AND.rg.LE.radlim(n+1)) THEN
              interp = n
            ENDIF
          ENDDO
          if(interp == 0) then
            IF(rg.LT.radlim(1))   interp = 1
            IF(rg.GT.radlim(25)) interp = 25
          endif
          m = interp - 1
          IF(interp.EQ.1) m = 1
          IF(interp.EQ.25) m = 23

          ! Determining the coefficients to be used in the 3rd order
          !Lagrangian polynomial interpolation scheme:
          L(1) = (rg-radlim(m+3))*(rg-radlim(m+2))*  &
                 (rg-radlim(m+1))/(radlim(m)-radlim(m+3))/  &
                 (radlim(m)-radlim(m+2))/(radlim(m)-radlim(m+1))
          L(2) = (rg-radlim(m+3))*(rg-radlim(m+2))*  &
                 (rg-radlim(m))/(radlim(m+1)-radlim(m+3))/  &
                 (radlim(m+1)-radlim(m+2))/(radlim(m+1)-radlim(m))
          L(3) = (rg-radlim(m+3))*(rg-radlim(m+1))*  &
                 (rg-radlim(m))/(radlim(m+2)-radlim(m+3))/  &
                 (radlim(m+2)-radlim(m+1))/(radlim(m+2)-radlim(m))
          L(4) = (rg-radlim(m+2))*(rg-radlim(m+1))*  &
                 (rg-radlim(m))/(radlim(m+3)-radlim(m+2))/  &
                 (radlim(m+3)-radlim(m+1))/(radlim(m+3)-radlim(m))

          !At each bin, determine the number concentration of aerosol
          DO i = 1,naerbins
            rnaer(k,i) = L(1)*radparam(m,i) + L(2)*radparam(m+1,i) +  &
                         L(3)*radparam(m+2,i) + L(4)*radparam(m+3,i)
            IF(rnaer(k,i).LT.0.) rnaer(k,i) = 0.
            rnaer(k,i) = rnaer(k,i) * numconc(k) 
          ENDDO ! Bin loop (i)

  ENDDO   ! Vertical level loop (k)

return
END SUBROUTINE aero_bin
end module mo_aerosol_optics
