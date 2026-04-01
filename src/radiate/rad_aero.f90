!##############################################################################
Subroutine aerorad (i,j,mb,nb,nrad,m1,dzl,relh,tp,omgp,gp,dn0,aodt)

! ib .......... band number
! mb .......... maximum number of bands (8)
! nb .......... total number of bands
! m1 .......... number of vertical levels
! dzl .......... delta z in each level (m)
! nrad - number of vertical levels in radiation calculations
! naerbins - number of aerosol bins for radiative calculations (17)
! nrh - number of relative humidity levels for radiative calculations (20)
! relh(m1) - relative humidity at each vertical level
! iaerorad - flag to turn on aerosol radiation
! aerocon(m1,aerocat) - aerosol category number concentration (#/kg)
! aeromas(m1,aerocat) - aerosol category mass (kg/kg)
! irh     - the relative humidity level (1=80%,20=99%)
! aerocat  - number of aerosol species in the model
! aerotype - aerosol type to be used through the radiation code
! rnaer(m1,naerbins) - aerosol number concentration in a given bin (#/m3)
!***************************************************************************
use micphys

implicit none

! Variable declaration:
INTEGER mb,nb,ib,nrad,m1,ibns,naerbins,irh,nrh,k,i,j
REAL dzl(nrad),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb),relh(m1),aodt
REAL rma(17),ext_sum,om_sum,asym_sum,con_sum,ext,om,gg
REAL gfact,qext,qscat,gasym,dm,pi

PARAMETER(naerbins=17,nrh=20)

INTEGER aerotype(aerocat)
REAL rnaer(m1,naerbins),aeroradius(m1,aerocat),dn0(m1)

PARAMETER (pi=3.141593)

! Drop radius at bin center in centimeters
DATA rma /2.0358125E-08,3.2316509E-08,5.1299262E-08                &
         ,8.1432499E-08,1.2926605E-07,2.0519707E-07,3.2573003E-07  &
         ,5.1706418E-07,8.2078826E-07,1.3029201E-06,2.0682569E-06  &
         ,3.2831531E-06,5.2116811E-06,8.2730281E-06,1.3132613E-05  &
         ,2.0846725E-05,3.3092112E-05/

! Initializing aerotype:
! Although there are multiple aerosol variables, there are currently only
! three types that are coded. Aerotype indicies are as follows:
! 1 - Ammonium sulfate
! 2 - Sea Salt
! 3 - Mineral Dust
! 4 - Absorbing Carbon 1
! 5 - Absorbing Carbon 2
do acat=1,aerocat
  aerotype(acat) = 0
enddo

if(iaerosol>0) then
   aerotype(1) = 1   ! CCN-mode-1
   aerotype(2) = 1   ! CCN-mode-2
endif
if(idust>0) then
   aerotype(3) = 3   ! Small dust mode
   aerotype(4) = 3   ! Large dust mode
endif
if(isalt>0) then
   aerotype(5) = 2   ! Salt film mode
   aerotype(6) = 2   ! Salt jet mode
   aerotype(7) = 2   ! Salt spume mode
endif
if(iabcarb>0) then
   aerotype(8) = 4   ! Absorbing carbon mode-1 (1% BC, 99% OC)
   aerotype(9) = 5   ! Absorbing carbon mode-2 (2% BC, 98% OC)
endif
if(iccnlev>=2) then
   aerotype(aerocat-1) = 1   ! Small regenerated aerosol
   aerotype(aerocat)   = 1   ! Large regenerated aerosol
endif

!Loop thru all aerosol species and set number and size
!Set aeroradius() at all levels, 1 through m1, since aerobin() looks
!at all levels.
do acat=1,aerocat
 if(aerotype(acat)>0) then
  do k = 1,m1
   if(aerocon(k,acat)>mincon .and. aeromas(k,acat)>=minmas) then
    aeroradius(k,acat)=((0.23873/aero_rhosol(acat) &
        *aeromas(k,acat)/aerocon(k,acat))**(1./3.))/aero_rg2rm(acat)
   else
    aeroradius(k,acat)=0.005e-6
   endif
  enddo
 endif
enddo

! Looping over each aerosol species. If the aerosol is turned off, the
! code checks for the next active aerosol:
aodt=0.0
DO acat=1,aerocat
 IF(aerotype(acat).gt.0) then
  ! For each aerosol, the code requires a total number concentration of
  ! the aerosol at each grid point, along with either a total mass at the
  ! grid point or a median radius.
  CALL aerobin (m1,naerbins,aerocon(1,acat),aeroradius(1,acat),rnaer,dn0)
  ! This is the routine that calculates the optical properties of the
  ! aerosols within the radiation routine:
  ! Doing over all vertical levels (not including radiation levels)
  DO k = 2,m1-1
    ! Locate RH as a percentage for table lookup
    irh=INT(100.*(relh(k)-0.80)) + 1
    irh=MAX(1,MIN(irh,20) )

    DO ib=1,nb

      ! Resetting temporary storage variabes to zero for next set of calcs:
      ext_sum = 0.0
      om_sum = 0.0
      asym_sum = 0.0
      con_sum = 0.0
      ext = 0.0
      om  = 0.0
      gg  = 0.0

      ! Calculating optical properties at over each aerosol bin:
      DO ibns = 1,naerbins
        ! Retrieving values for the growth factor, the extinction coefficient,
        ! the scattering coefficient and the asymmetry parameter, from their
        ! appropriate look up tables:
        CALL aerogrowth (aerotype(acat),ibns,irh,gfact)
        CALL aeroqext   (aerotype(acat),ib,ibns,irh,qext)
        CALL aeroqscat  (aerotype(acat),ib,ibns,irh,qscat)
        CALL aerogasym  (aerotype(acat),ib,ibns,irh,gasym)
        ! Determining mean diameter of bin (meters), with deliquescence growth factor
        dm = 2. * rma(ibns) * gfact
        ! Updating temporary storage variables:
        if(qext>0.0 .and. rnaer(k,ibns)>0.0)then
          ext_sum  = ext_sum  + (pi/4. * dm**2 * rnaer(k,ibns) * qext)  !Bext
          om_sum   = om_sum   + (pi/4. * dm**2 * rnaer(k,ibns) * qscat) !Bscat
          asym_sum = asym_sum + (rnaer(k,ibns) * gasym) !Asymmetry parameter
          con_sum  = con_sum  + rnaer(k,ibns) !Total aerosol number
        endif
      ENDDO ! Aerosol bins (ibns)
      !Set up ext, om, gg same as Harrington's cloud_top for consistency
      ext = ext_sum * dzl(k) !Tau = Bext * dz
      if(ext>0.0) then
        om   = om_sum / ext_sum  !Omega = Bscat / Bext
        gg   = asym_sum/(con_sum+1.E-30) !normalize asym_sum by total number
      else
        ext = 0.0
      endif
      ! Obtain values of optical depth (tp), single scatter albedo (omgp),
      ! asymmetry parameter (gp) in a form to interact with cloud_opt.
      ! Some values of dzl may be <0; in this case set tp to 0.
      ! Factor in layer depth for optical depth extinction
      aodt       = aodt        + ext
      tp(k,ib)   = tp(k,ib)    + ext
      omgp(k,ib) = omgp(k,ib)  + om * ext
      gp(k,ib)   = gp(k,ib)    + gg * om * ext
    ENDDO ! Radiation band loop
  ENDDO ! Vertical level loop
 ENDIF !Aerosol type if
ENDDO ! Aerotype loop

return
END SUBROUTINE aerorad

!##############################################################################
Subroutine aerobin (m1,naerbins,numconc,medianrad,rnaer,dn0)

! The purpose of this routine was to take a number concentration
! and median radius for any given aerosol and turn it into a binned
! distribution
! Input Variables: m1 - Integer, number of vertical levels
!    naerbins  - Integer number of bins used in aerosol radiation
!                calculations.  Set to 17 in radcalc3.
!    numconc   - Number concentration of aerosol (#/kg)
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

! Internal variable declaraion:
        INTEGER bin,i,k,n,m,interp
        REAL radlim(26),radparam(26,17),rg,L(4)

! Outgoing variable declaration:
        REAL rnaer(m1,naerbins)

! Data statement containing the radii for which the binning routine was
! run off-line...i.e., explicit binning done at these radii:
        DATA radlim/0.00945E-06, 0.01191E-06, 0.01500E-06, 0.01890E-06  &
                   ,0.02382E-06, 0.03000E-06, 0.03780E-06, 0.04760E-06  &
                   ,0.06000E-06, 0.07560E-06, 0.09520E-06, 0.12000E-06  &
                   ,0.15120E-06, 0.19050E-06, 0.24000E-06, 0.30240E-06  &
                   ,0.38100E-06, 0.48000E-06, 0.60500E-06, 0.76200E-06  &
                   ,0.96000E-06, 1.20900E-06, 1.52400E-06, 1.92000E-06  &
                   ,2.41900E-06, 3.04800E-06/ ! in meters

! Percentage of mass within each bin at the above given radii:
        DATA (radparam(1,bin),bin=1,17)  & ! at 9.45nm
        /0.1580, 0.0489, 0.0084, 0.0008, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(2,bin),bin=1,17)  & ! at 11.91nm
        /0.2280, 0.0947, 0.0218, 0.0028, 0.0002, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(3,bin),bin=1,17)  & ! at 15nm
        /0.2839, 0.1580, 0.0488, 0.0083, 0.0008, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(4,bin),bin=1,17)  & ! at 18.9nm
        /0.3055, 0.2280, 0.0945, 0.0217, 0.0028, 0.0002  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(5,bin),bin=1,17)  & ! at 23.82nm
        /0.2839, 0.2841, 0.1581, 0.0487, 0.0084, 0.0008  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(6,bin),bin=1,17)  & ! at 30nm
        /0.2279, 0.3056, 0.2280, 0.0943, 0.0217, 0.0028  &
        ,0.0002, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(7,bin),bin=1,17)  & ! at 37.8nm
        /0.1579, 0.2840, 0.2841, 0.1577, 0.0488, 0.0084  &
        ,0.0008, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(8,bin),bin=1,17)  & ! at 47.6nm
        /0.0946, 0.2280, 0.3057, 0.2276, 0.0943, 0.0217  &
        ,0.0028, 0.0002, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(9,bin),bin=1,17)  & ! at 60nm
        /0.0488, 0.1579, 0.2481, 0.2838, 0.1579, 0.0488  &
        ,0.0084, 0.0008, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(10,bin),bin=1,17)  & ! at 75.6nm
        /0.0218, 0.0945, 0.2280, 0.3055, 0.2278, 0.0946  &
        ,0.0218, 0.0028, 0.0002, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(11,bin),bin=1,17)  & ! at 95.2nm
        /0.0084, 0.0489, 0.1581, 0.2841, 0.2838, 0.1578  &
        ,0.0488, 0.0083, 0.0008, 0.0000, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(12,bin),bin=1,17)  & ! at 120nm
        /0.0028, 0.0217, 0.0945, 0.2280, 0.3055, 0.2279  &
        ,0.0946, 0.0217, 0.0028, 0.0002, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(13,bin),bin=1,17)  & ! at 151.2nm
        /0.0008, 0.0084, 0.0488, 0.1580, 0.2839, 0.2840  &
        ,0.1581, 0.0487, 0.0083, 0.0008, 0.0000, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(14,bin),bin=1,17)  & ! at 190.5nm
        /0.0002, 0.0028, 0.0217, 0.0945, 0.2279, 0.3056  &
        ,0.2281, 0.0944, 0.0217, 0.0028, 0.0002, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(15,bin),bin=1,17)  & ! at 240nm
        /0.0000, 0.0008, 0.0083, 0.0488, 0.1580, 0.2840  &
        ,0.2841, 0.1578, 0.0487, 0.0084, 0.0008, 0.0000  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(16,bin),bin=1,17)  & ! at 302.4 nm
        /0.0000, 0.0002, 0.0028, 0.0218, 0.0945, 0.2279  &
        ,0.3057, 0.2278, 0.0944, 0.0218, 0.0028, 0.0002  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(17,bin),bin=1,17)  & ! at 381nm
        /0.0000, 0.0000, 0.0008, 0.0084, 0.0488, 0.1579  &
        ,0.2841, 0.2840, 0.1578, 0.0488, 0.0084, 0.0008  &
        ,0.0000, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(18,bin),bin=1,17)  & ! at 480nm
        /0.0000, 0.0000, 0.0002, 0.0028, 0.0215, 0.0945  &
        ,0.2279, 0.3056, 0.2277, 0.0945, 0.0218, 0.0028  &
        ,0.0002, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(19,bin),bin=1,17)  & ! at 605nm
        /0.0000, 0.0000, 0.0000, 0.0008, 0.0084, 0.0487  &
        ,0.1587, 0.2840, 0.2839, 0.1580, 0.0489, 0.0084  &
        ,0.0008, 0.0000, 0.0000, 0.0000, 0.0000/

        DATA (radparam(20,bin),bin=1,17)  & ! at 762nm
        /0.0000, 0.0000, 0.0000, 0.0002, 0.0028, 0.0217  &
        ,0.0944, 0.2280, 0.3055, 0.2279, 0.0946, 0.0218  &
        ,0.0028, 0.0002, 0.0000, 0.0000, 0.0000/

        DATA (radparam(21,bin),bin=1,17)  & ! at 960nm
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0008, 0.0084  &
        ,0.0488, 0.1580, 0.2840, 0.2839, 0.1580, 0.0488  &
        ,0.0083, 0.0008, 0.0000, 0.0000, 0.0000/

        DATA (radparam(22,bin),bin=1,17)  & ! at 1.209 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0002, 0.0028  &
        ,0.0218, 0.0946, 0.2281, 0.3055, 0.2278, 0.0945  &
        ,0.0216, 0.0028, 0.0002, 0.0000, 0.0000/

        DATA (radparam(23,bin),bin=1,17)  & ! at 1.524 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0008  &
        ,0.0083, 0.0488, 0.1580, 0.2839, 0.2840, 0.1581  &
        ,0.0486, 0.0084, 0.0008, 0.0000, 0.0000/

        DATA (radparam(24,bin),bin=1,17)  & ! at 1.92 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0002  &
        ,0.0028, 0.0217, 0.0945, 0.2279, 0.3056, 0.2281  &
        ,0.0943, 0.0217, 0.0028, 0.0002, 0.0000/

        DATA (radparam(25,bin),bin=1,17)  & ! at 2.419 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0008, 0.0084, 0.0488, 0.1580, 0.2840, 0.2842  &
        ,0.1577, 0.0487, 0.0084, 0.0008, 0.0000/

        DATA (radparam(26,bin),bin=1,17)  & ! at 3.048 microns
        /0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000  &
        ,0.0002, 0.0028, 0.0218, 0.0945, 0.2279, 0.3058  &
        ,0.2277, 0.0944, 0.0218, 0.0028, 0.0002/

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
            ! Conversion from #/kg to #/m3:
            rnaer(k,i) = rnaer(k,i) * numconc(k) * dn0(k)
          ENDDO ! Bin loop (i)

  ENDDO   ! Vertical level loop (k)

return
END SUBROUTINE aerobin

!##############################################################################
Subroutine aerogrowth (aerotype,bin,rh,value)

implicit none

        INTEGER aerotype,bin,rh,ib
        REAL growth(5,17,20),value

! Ammonium sulfate with 20% insoluble inclusion (dust low R(Im)))
        DATA(growth(1,ib,1),ib=1,17)  & ! Amm. Sulf. at  80% RH
        / 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000  &
        , 1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000  &
        , 1.0000,  1.0000,  1.0000,  1.0000,  1.0000/

        DATA(growth(1,ib,2),ib=1,17)  & ! Amm. Sulf. at  81% RH
        / 1.0100,  1.0094,  1.0099,  1.0102,  1.0104,  1.0105  &
        , 1.0106,  1.0106,  1.0107,  1.0107,  1.0107,  1.0107  &
        , 1.0107,  1.0107,  1.0107,  1.0107,  1.0107/

        DATA(growth(1,ib,3),ib=1,17)  & ! Amm. Sulf. at  82% RH
        / 1.0195,  1.0195,  1.0205,  1.0212,  1.0216,  1.0219  &
        , 1.0220,  1.0221,  1.0222,  1.0222,  1.0223,  1.0223  &
        , 1.0223,  1.0223,  1.0223,  1.0223,  1.0223/

        DATA(growth(1,ib,4),ib=1,17)  & ! Amm. Sulf. at  83% RH
        / 1.0295,  1.0305,  1.0320,  1.0331,  1.0337,  1.0342  &
        , 1.0344,  1.0346,  1.0347,  1.0348,  1.0348,  1.0349  &
        , 1.0349,  1.0349,  1.0349,  1.0349,  1.0349/

        DATA(growth(1,ib,5),ib=1,17)  & ! Amm. Sulf. at  84% RH
        / 1.0404,  1.0423,  1.0445,  1.0460,  1.0469,  1.0476  &
        , 1.0480,  1.0482,  1.0484,  1.0485,  1.0485,  1.0486  &
        , 1.0486,  1.0486,  1.0486,  1.0486,  1.0486/

        DATA(growth(1,ib,6),ib=1,17)  & ! Amm. Sulf. at  85% RH
        / 1.0522,  1.0551,  1.0581,  1.0601,  1.0614,  1.0622  &
        , 1.0628,  1.0631,  1.0633,  1.0635,  1.0635,  1.0636  &
        , 1.0636,  1.0637,  1.0637,  1.0637,  1.0637/

        DATA(growth(1,ib,7),ib=1,17)  & ! Amm. Sulf. at  86% RH
        / 1.0650,  1.0692,  1.0731,  1.0757,  1.0774,  1.0785  &
        , 1.0792,  1.0796,  1.0799,  1.0800,  1.0802,  1.0802  &
        , 1.0803,  1.0803,  1.0803,  1.0803,  1.0803/

        DATA(growth(1,ib,8),ib=1,17)  & ! Amm. Sulf. at  87% RH
        / 1.0790,  1.0847,  1.0897,  1.0930,  1.0951,  1.0965  &
        , 1.0974,  1.0980,  1.0983,  1.0985,  1.0987,  1.0988  &
        , 1.0988,  1.0989,  1.0989,  1.0989,  1.0989/

        DATA(growth(1,ib,9),ib=1,17)  & ! Amm. Sulf. at  88% RH
        / 1.0945,  1.1019,  1.1081,  1.1123,  1.1150,  1.1168  &
        , 1.1179,  1.1186,  1.1191,  1.1194,  1.1195,  1.1197  &
        , 1.1197,  1.1198,  1.1198,  1.1198,  1.1198/

        DATA(growth(1,ib,10),ib=1,17)  & ! Amm. Sulf. at  89% RH
        / 1.1117,  1.1213,  1.1290,  1.1342,  1.1376,  1.1398  &
        , 1.1412,  1.1422,  1.1427,  1.1431,  1.1433,  1.1435  &
        , 1.1435,  1.1436,  1.1436,  1.1437,  1.1437/

        DATA(growth(1,ib,11),ib=1,17)  & ! Amm. Sulf. at  90% RH
        / 1.1310,  1.1432,  1.1528,  1.1593,  1.1636,  1.1664  &
        , 1.1682,  1.1693,  1.1700,  1.1705,  1.1708,  1.1709  &
        , 1.1711,  1.1711,  1.1712,  1.1712,  1.1712/

        DATA(growth(1,ib,12),ib=1,17)  & ! Amm. Sulf. at  91% RH
        / 1.1529,  1.1684,  1.1804,  1.1886,  1.1940,  1.1975  &
        , 1.1997,  1.2011,  1.2020,  1.2026,  1.2030,  1.2032  &
        , 1.2034,  1.2034,  1.2035,  1.2035,  1.2036/

        DATA(growth(1,ib,13),ib=1,17)  & ! Amm. Sulf. at  92% RH
        / 1.1782,  1.1978,  1.2129,  1.2232,  1.2301,  1.2345  &
        , 1.2374,  1.2392,  1.2403,  1.2411,  1.2415,  1.2418  &
        , 1.2420,  1.2421,  1.2422,  1.2422,  1.2422/

        DATA(growth(1,ib,14),ib=1,17)  & ! Amm. Sulf. at  93% RH
        / 1.2079,  1.2329,  1.2521,  1.2653,  1.2740,  1.2797  &
        , 1.2833,  1.2856,  1.2871,  1.2880,  1.2886,  1.2890  &
        , 1.2892,  1.2894,  1.2895,  1.2895,  1.2896/

        DATA(growth(1,ib,15),ib=1,17)  & ! Amm. Sulf. at  94% RH
        / 1.2435,  1.2758,  1.3005,  1.3175,  1.3288,  1.3362  &
        , 1.3409,  1.3440,  1.3459,  1.3471,  1.3479,  1.3484  &
        , 1.3487,  1.3489,  1.3490,  1.3491,  1.3491/

        DATA(growth(1,ib,16),ib=1,17)  & ! Amm. Sulf. at  95% RH
        / 1.2874,  1.3298,  1.3623,  1.3847,  1.3996,  1.4094  &
        , 1.4157,  1.4197,  1.4222,  1.4238,  1.4249,  1.4255  &
        , 1.4259,  1.4262,  1.4263,  1.4264,  1.4265/

        DATA(growth(1,ib,17),ib=1,17)  & ! Amm. Sulf. at  96% RH
        / 1.3435,  1.4006,  1.4443,  1.4747,  1.4951,  1.5083  &
        , 1.5169,  1.5224,  1.5258,  1.5280,  1.5294,  1.5303  &
        , 1.5309,  1.5312,  1.5314,  1.5316,  1.5317/

        DATA(growth(1,ib,18),ib=1,17)  & ! Amm. Sulf. at  97% RH
        / 1.4184,  1.4982,  1.5597,  1.6028,  1.6324,  1.6517  &
        , 1.6643,  1.6723,  1.6774,  1.6806,  1.6826,  1.6839  &
        , 1.6847,  1.6852,  1.6856,  1.6858,  1.6859/

        DATA(growth(1,ib,19),ib=1,17)  & ! Amm. Sulf. at  98% RH
        / 1.5256,  1.6446,  1.7397,  1.8075,  1.8536,  1.8840  &
        , 1.9037,  1.9164,  1.9245,  1.9296,  1.9328,  1.9349  &
        , 1.9362,  1.9370,  1.9375,  1.9379,  1.9381/

        DATA(growth(1,ib,20),ib=1,17)  & ! Amm. Sulf. at  99% RH
        / 1.6978,  1.9026,  2.0745,  2.2054,  2.2983,  2.3613  &
        , 2.4028,  2.4298,  2.4471,  2.4581,  2.4651,  2.4695  &
        , 2.4723,  2.4740,  2.4751,  2.4758,  2.4763/

! Sea salt with no insoluble inclusion
        DATA(growth(2,ib,1),ib=1,17)  & ! Sea Salt at  80% RH
        /1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000  &
        ,1.0000, 1.0000, 1.0000, 1.0000, 1.0000/

        DATA(growth(2,ib,2),ib=1,17)  & ! Sea Salt at  81% RH
        /1.0151, 1.0142, 1.0134, 1.0136, 1.0137, 1.0137  &
        ,1.0138, 1.0138, 1.0138, 1.0138, 1.0138, 1.0138  &
        ,1.0138, 1.0138, 1.0138, 1.0138, 1.0138/

        DATA(growth(2,ib,3),ib=1,17)  & ! Sea Salt at  82% RH
        /1.0295, 1.0281, 1.0276, 1.0279, 1.0281, 1.0282  &
        ,1.0283, 1.0283, 1.0283, 1.0284, 1.0284, 1.0284  &
        ,1.0284, 1.0284, 1.0284, 1.0284, 1.0284/

        DATA(growth(2,ib,4),ib=1,17)  & ! Sea Salt at  83% RH
        /1.0437, 1.0428, 1.0425, 1.0430, 1.0433, 1.0435  &
        ,1.0437, 1.0438, 1.0438, 1.0439, 1.0439, 1.0439  &
        ,1.0439, 1.0439, 1.0439, 1.0439, 1.0439/

        DATA(growth(2,ib,5),ib=1,17)  & ! Sea Salt at  84% RH
        /1.0588, 1.0584, 1.0584, 1.0592, 1.0597, 1.0600  &
        ,1.0602, 1.0603, 1.0604, 1.0605, 1.0605, 1.0605  &
        ,1.0606, 1.0606, 1.0606, 1.0606, 1.0606/

        DATA(growth(2,ib,6),ib=1,17)  & ! Sea Salt at  85% RH
        /1.0748, 1.0750, 1.0756, 1.0767, 1.0774, 1.0778  &
        ,1.0781, 1.0783, 1.0784, 1.0785, 1.0785, 1.0786  &
        ,1.0786, 1.0786, 1.0786, 1.0786, 1.0786/

        DATA(growth(2,ib,7),ib=1,17)  & ! Sea Salt at  86% RH
        /1.0920, 1.0931, 1.0942, 1.0957, 1.0966, 1.0973  &
        ,1.0977, 1.0979, 1.0981, 1.0982, 1.0982, 1.0983  &
        ,1.0983, 1.0983, 1.0983, 1.0983, 1.0983/

        DATA(growth(2,ib,8),ib=1,17)  & ! Sea Salt at  87% RH
        /1.1107, 1.1127, 1.1146, 1.1166, 1.1179, 1.1187  &
        ,1.1192, 1.1196, 1.1198, 1.1199, 1.1200, 1.1201  &
        ,1.1201, 1.1201, 1.1201, 1.1201, 1.1201/

        DATA(growth(2,ib,9),ib=1,17)  & ! Sea Salt at  88% RH
        /1.1312, 1.1345, 1.1372, 1.1398, 1.1415, 1.1426  &
        ,1.1433, 1.1437, 1.1440, 1.1442, 1.1443, 1.1444  &
        ,1.1444, 1.1444, 1.1444, 1.1445, 1.1445/

        DATA(growth(2,ib,10),ib=1,17)  & ! Sea Salt at  89% RH
        /1.1539, 1.1588, 1.1626, 1.1659, 1.1681, 1.1695  &
        ,1.1704, 1.1710, 1.1713, 1.1715, 1.1717, 1.1718  &
        ,1.1718, 1.1719, 1.1719, 1.1719, 1.1719/

        DATA(growth(2,ib,11),ib=1,17)  & ! Sea Salt at  90% RH
        /1.1794, 1.1862, 1.1914, 1.1956, 1.1984, 1.2002  &
        ,1.2013, 1.2020, 1.2025, 1.2028, 1.2029, 1.2030  &
        ,1.2031, 1.2032, 1.2032, 1.2032, 1.2032/

        DATA(growth(2,ib,12),ib=1,17)  & ! Sea Salt at  91% RH
        /1.2083, 1.2177, 1.2245, 1.2298, 1.2333, 1.2356  &
        ,1.2370, 1.2379, 1.2385, 1.2388, 1.2390, 1.2392  &
        ,1.2393, 1.2393, 1.2394, 1.2394, 1.2394/

        DATA(growth(2,ib,13),ib=1,17)  & ! Sea Salt at  92% RH
        /1.2417, 1.2541, 1.2631, 1.2698, 1.2742, 1.2770  &
        ,1.2788, 1.2799, 1.2807, 1.2811, 1.2814, 1.2816  &
        ,1.2817, 1.2818, 1.2818, 1.2819, 1.2819/

        DATA(growth(2,ib,14),ib=1,17)  & ! Sea Salt at  93% RH
        /1.2808, 1.2971, 1.3088, 1.3173, 1.3229, 1.3265  &
        ,1.3288, 1.3302, 1.3311, 1.3317, 1.3321, 1.3323  &
        ,1.3324, 1.3325, 1.3326, 1.3326, 1.3327/

        DATA(growth(2,ib,15),ib=1,17)  & ! Sea Salt at  94% RH
        /1.3275, 1.3489, 1.3642, 1.3751, 1.3823, 1.3869  &
        ,1.3899, 1.3917, 1.3929, 1.3937, 1.3941, 1.3944  &
        ,1.3946, 1.3947, 1.3948, 1.3948, 1.3949/

        DATA(growth(2,ib,16),ib=1,17)  & ! Sea Salt at  95% RH
        /1.3846, 1.4131, 1.4334, 1.4478, 1.4572, 1.4632  &
        ,1.4671, 1.4696, 1.4711, 1.4721, 1.4727, 1.4731  &
        ,1.4734, 1.4735, 1.4736, 1.4737, 1.4737/

        DATA(growth(2,ib,17),ib=1,17)  & ! Sea Salt at  96% RH
        /1.4570, 1.4959, 1.5235, 1.5431, 1.5559, 1.5642  &
        ,1.5696, 1.5729, 1.5751, 1.5764, 1.5773, 1.5778  &
        ,1.5782, 1.5784, 1.5785, 1.5786, 1.5786/

        DATA(growth(2,ib,18),ib=1,17)  & ! Sea Salt at  97% RH
        /1.5534, 1.6089, 1.6489, 1.6770, 1.6956, 1.7077  &
        ,1.7154, 1.7204, 1.7235, 1.7255, 1.7267, 1.7275  &
        ,1.7280, 1.7283, 1.7285, 1.7287, 1.7287/

        DATA(growth(2,ib,19),ib=1,17)  & ! Sea Salt at  98% RH
        /1.6927, 1.7790, 1.8425, 1.8875, 1.9175, 1.9371  &
        ,1.9497, 1.9578, 1.9629, 1.9662, 1.9682, 1.9695  &
        ,1.9704, 1.9709, 1.9712, 1.9714, 1.9715/

        DATA(growth(2,ib,20),ib=1,17)  & ! Sea Salt at  99% RH
        /1.9269, 2.0871, 2.2119, 2.3031, 2.3656, 2.4071  &
        ,2.4341, 2.4515, 2.4626, 2.4696, 2.4741, 2.4769  &
        ,2.4787, 2.4798, 2.4805, 2.4809, 2.4812/

        value = growth(aerotype,bin,rh)

! No growth factor for dust or absorbing carbon at RH under 100%:
        IF(AEROTYPE.EQ.3) value = 1.0 ! Mineral dust
        IF(AEROTYPE.EQ.4) value = 1.0 ! Absorbing carbon (1% BC, 99% OC)
        IF(AEROTYPE.EQ.5) value = 1.0 ! Absorbing carbon (2% BC, 98% OC)
return
END SUBROUTINE aerogrowth

!##############################################################################
Subroutine aeroqext (aerotype,radband,bin,rh,value)

use rrad3, only:dust_ref_im

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL qext(2,8,17,20),value
        REAL qext_dust(3,8,17)
        REAL qext_carb(2,8,17)

! Qext Ammonium sulfate with 20% insoluble inclusion (dust low R(Im)))
        DATA(qext(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        / 0.00243,  0.00399,  0.00676,  0.01291,  0.03101,  0.11152  &
        , 0.39298,  1.20280,  2.59550,  3.42510,  2.50370,  2.36510  &
        , 2.25050,  2.18720,  2.14660,  2.10300,  2.07280/

        DATA(qext(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        / 0.00247,  0.00405,  0.00686,  0.01315,  0.03175,  0.11430  &
        , 0.40095,  1.22000,  2.61630,  3.42130,  2.49150,  2.34430  &
        , 2.24030,  2.17830,  2.14220,  2.09950,  2.07250/

        DATA(qext(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        / 0.00251,  0.00411,  0.00698,  0.01342,  0.03257,  0.11734  &
        , 0.40973,  1.23880,  2.63900,  3.41530,  2.47780,  2.34960  &
        , 2.23690,  2.18440,  2.13900,  2.09960,  2.07340/

        DATA(qext(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        / 0.00254,  0.00418,  0.00710,  0.01371,  0.03347,  0.12070  &
        , 0.41949,  1.25960,  2.66360,  3.40690,  2.47100,  2.32250  &
        , 2.24880,  2.17460,  2.14160,  2.10060,  2.07520/

        DATA(qext(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        / 0.00258,  0.00425,  0.00724,  0.01403,  0.03447,  0.12441  &
        , 0.43040,  1.28260,  2.69010,  3.39640,  2.46560,  2.33550  &
        , 2.23730,  2.18500,  2.14040,  2.10210,  2.07500/

        DATA(qext(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        / 0.00263,  0.00433,  0.00739,  0.01438,  0.03559,  0.12856  &
        , 0.44268,  1.30810,  2.71830,  3.38460,  2.45840,  2.30650  &
        , 2.24240,  2.17900,  2.13650,  2.10230,  2.07490/

        DATA(qext(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        / 0.00267,  0.00441,  0.00755,  0.01477,  0.03685,  0.13322  &
        , 0.45665,  1.33650,  2.74850,  3.37330,  2.44610,  2.30040  &
        , 2.24060,  2.17600,  2.13390,  2.10200,  2.07110/

        DATA(qext(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        / 0.00272,  0.00450,  0.00773,  0.01522,  0.03829,  0.13851  &
        , 0.47267,  1.36830,  2.78070,  3.36300,  2.43130,  2.30430  &
        , 2.24940,  2.18380,  2.13010,  2.09760,  2.07030/

        DATA(qext(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        / 0.00278,  0.00461,  0.00793,  0.01572,  0.03995,  0.14459  &
        , 0.49126,  1.40410,  2.81570,  3.34940,  2.42830,  2.28660  &
        , 2.24990,  2.18560,  2.12480,  2.09710,  2.07020/

        DATA(qext(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        / 0.00284,  0.00472,  0.00815,  0.01630,  0.04190,  0.15167  &
        , 0.51309,  1.44470,  2.85490,  3.32690,  2.42480,  2.27210  &
        , 2.24610,  2.17890,  2.12900,  2.09430,  2.06870/

        DATA(qext(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        / 0.00291,  0.00485,  0.00841,  0.01698,  0.04423,  0.16002  &
        , 0.53905,  1.49120,  2.90080,  3.29250,  2.41490,  2.28060  &
        , 2.24940,  2.16770,  2.12270,  2.09130,  2.06700/

        DATA(qext(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        / 0.00298,  0.00499,  0.00871,  0.01779,  0.04706,  0.17009  &
        , 0.57036,  1.54560,  2.95660,  3.25190,  2.41590,  2.26050  &
        , 2.24660,  2.16400,  2.13220,  2.09290,  2.06840/

        DATA(qext(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        / 0.00307,  0.00516,  0.00907,  0.01878,  0.05060,  0.18246  &
        , 0.60865,  1.61090,  3.02240,  3.21000,  2.41410,  2.25750  &
        , 2.24660,  2.15570,  2.12580,  2.08810,  2.06440/

        DATA(qext(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        / 0.00317,  0.00536,  0.00950,  0.02004,  0.05516,  0.19810  &
        , 0.65619,  1.69260,  3.09320,  3.14050,  2.42580,  2.26270  &
        , 2.22030,  2.16260,  2.11180,  2.08570,  2.06260/

        DATA(qext(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        / 0.00329,  0.00561,  0.01005,  0.02168,  0.06127,  0.21852  &
        , 0.71628,  1.79930,  3.16720,  3.05150,  2.44450,  2.29230  &
        , 2.20420,  2.16790,  2.12050,  2.08490,  2.06120/

        DATA(qext(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        / 0.00344,  0.00591,  0.01076,  0.02394,  0.06987,  0.24643  &
        , 0.79436,  1.93990,  3.25950,  2.92110,  2.46390,  2.31160  &
        , 2.20510,  2.15580,  2.10930,  2.08210,  2.05930/

        DATA(qext(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        / 0.00362,  0.00632,  0.01173,  0.02725,  0.08275,  0.28716  &
        , 0.90224,  2.12220,  3.34370,  2.75870,  2.47640,  2.30890  &
        , 2.20160,  2.14420,  2.10500,  2.07860,  2.05570/

        DATA(qext(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        / 0.00387,  0.00687,  0.01318,  0.03260,  0.10391,  0.35406  &
        , 1.07260,  2.38460,  3.40860,  2.54840,  2.42030,  2.26990  &
        , 2.18920,  2.13010,  2.10390,  2.07290,  2.05390/

        DATA(qext(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        / 0.00421,  0.00772,  0.01567,  0.04297,  0.14416,  0.48717  &
        , 1.37270,  2.76550,  3.33610,  2.38330,  2.29050,  2.22400  &
        , 2.17340,  2.11910,  2.09050,  2.06590,  2.04670/

        DATA(qext(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        / 0.00476,  0.00927,  0.02126,  0.07066,  0.24497,  0.81293  &
        , 1.99920,  3.29440,  2.81910,  2.45820,  2.32640,  2.20730  &
        , 2.14340,  2.10410,  2.07720,  2.05600,  2.04400/

        DATA(qext(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        / 0.00027,  0.00152,  0.00939,  0.05618,  0.24749,  0.91117  &
        , 2.28680,  3.57540,  2.87780,  2.38950,  2.34530,  2.22210  &
        , 2.16090,  2.12210,  2.08480,  2.05860,  2.04530/

        DATA(qext(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        / 0.00028,  0.00156,  0.00967,  0.05772,  0.25225,  0.92267  &
        , 2.30370,  3.57880,  2.86350,  2.38570,  2.34270,  2.22560  &
        , 2.16300,  2.12240,  2.08490,  2.06410,  2.04570/

        DATA(qext(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        / 0.00028,  0.00161,  0.00997,  0.05942,  0.25746,  0.93529  &
        , 2.32200,  3.58150,  2.84940,  2.37620,  2.33370,  2.22180  &
        , 2.16210,  2.11920,  2.08300,  2.05990,  2.04340/

        DATA(qext(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        / 0.00029,  0.00166,  0.01030,  0.06129,  0.26323,  0.94922  &
        , 2.34190,  3.58310,  2.83290,  2.36780,  2.33150,  2.22060  &
        , 2.16130,  2.11580,  2.07840,  2.05860,  2.04320/

        DATA(qext(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        / 0.00030,  0.00171,  0.01068,  0.06339,  0.26963,  0.96473  &
        , 2.36360,  3.58360,  2.81400,  2.35970,  2.32380,  2.22240  &
        , 2.15620,  2.10910,  2.08300,  2.06230,  2.04360/

        DATA(qext(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        / 0.00030,  0.00177,  0.01111,  0.06574,  0.27681,  0.98216  &
        , 2.38760,  3.58320,  2.79190,  2.35370,  2.33020,  2.21150  &
        , 2.15640,  2.10680,  2.08120,  2.06070,  2.04310/

        DATA(qext(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        / 0.00031,  0.00185,  0.01159,  0.06840,  0.28494,  1.00200  &
        , 2.41450,  3.58260,  2.76650,  2.34540,  2.31160,  2.20840  &
        , 2.15290,  2.11080,  2.07660,  2.05610,  2.04320/

        DATA(qext(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        / 0.00032,  0.00193,  0.01215,  0.07144,  0.29424,  1.02470  &
        , 2.44510,  3.58330,  2.74170,  2.33360,  2.31470,  2.21930  &
        , 2.14640,  2.10710,  2.07700,  2.05940,  2.04150/

        DATA(qext(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        / 0.00033,  0.00202,  0.01279,  0.07496,  0.30502,  1.05120  &
        , 2.48070,  3.58620,  2.71940,  2.31750,  2.31040,  2.21110  &
        , 2.14470,  2.10510,  2.07910,  2.05890,  2.04050/

        DATA(qext(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        / 0.00035,  0.00213,  0.01356,  0.07909,  0.31774,  1.08250  &
        , 2.52280,  3.58930,  2.69160,  2.30190,  2.30670,  2.20850  &
        , 2.14200,  2.10380,  2.07760,  2.05620,  2.04060/

        DATA(qext(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        / 0.00036,  0.00226,  0.01448,  0.08402,  0.33303,  1.12000  &
        , 2.57280,  3.58780,  2.65470,  2.29190,  2.29270,  2.21020  &
        , 2.13970,  2.10570,  2.07410,  2.05800,  2.04180/

        DATA(qext(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        / 0.00038,  0.00243,  0.01562,  0.09000,  0.35180,  1.16550  &
        , 2.63190,  3.57800,  2.62000,  2.26910,  2.28690,  2.21430  &
        , 2.13330,  2.10200,  2.07340,  2.05390,  2.03890/

        DATA(qext(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        / 0.00041,  0.00263,  0.01705,  0.09743,  0.37547,  1.22160  &
        , 2.70050,  3.56240,  2.58150,  2.25730,  2.28720,  2.20910  &
        , 2.13460,  2.10530,  2.07250,  2.05230,  2.03870/

        DATA(qext(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        / 0.00044,  0.00289,  0.01891,  0.10692,  0.40620,  1.29140  &
        , 2.78000,  3.54560,  2.53900,  2.24550,  2.28160,  2.19530  &
        , 2.13030,  2.10140,  2.07330,  2.05160,  2.03710/

        DATA(qext(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        / 0.00048,  0.00324,  0.02142,  0.11943,  0.44743,  1.37940  &
        , 2.87770,  3.50280,  2.49830,  2.24030,  2.26540,  2.18070  &
        , 2.13580,  2.09560,  2.07170,  2.04990,  2.03640/

        DATA(qext(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        / 0.00053,  0.00372,  0.02500,  0.13665,  0.50455,  1.49420  &
        , 3.00770,  3.43500,  2.46080,  2.25430,  2.24370,  2.16100  &
        , 2.13350,  2.09020,  2.06810,  2.04760,  2.03490/

        DATA(qext(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        / 0.00061,  0.00445,  0.03042,  0.16163,  0.58614,  1.65500  &
        , 3.16140,  3.31230,  2.42970,  2.29580,  2.19720,  2.16000  &
        , 2.12040,  2.08600,  2.06200,  2.04570,  2.03230/

        DATA(qext(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        / 0.00073,  0.00565,  0.03942,  0.20094,  0.70828,  1.89860  &
        , 3.34780,  3.10740,  2.40100,  2.34250,  2.18720,  2.15740  &
        , 2.11290,  2.08090,  2.05860,  2.04270,  2.03000/

        DATA(qext(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        / 0.00094,  0.00790,  0.05687,  0.27466,  0.91940,  2.27830  &
        , 3.51960,  2.77250,  2.32350,  2.31360,  2.20490,  2.13130  &
        , 2.09910,  2.07630,  2.05460,  2.03850,  2.02740/

        DATA(qext(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        / 0.00138,  0.01347,  0.10048,  0.46808,  1.40640,  2.96310  &
        , 3.40340,  2.42910,  2.25600,  2.21130,  2.14870,  2.11730  &
        , 2.08700,  2.06300,  2.04550,  2.03180,  2.01970/

        DATA(qext(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        / 0.00380,  0.02378,  0.12766,  0.52465,  1.52200,  3.05850  &
        , 3.54690,  2.43980,  2.39540,  2.20440,  2.19480,  2.12910  &
        , 2.09410,  2.07630,  2.05230,  2.03650,  2.02610/

        DATA(qext(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        / 0.00391,  0.02440,  0.13047,  0.53312,  1.53760,  3.07170  &
        , 3.53820,  2.42730,  2.39430,  2.19750,  2.19920,  2.13160  &
        , 2.09740,  2.07660,  2.05020,  2.03670,  2.02550/

        DATA(qext(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        / 0.00401,  0.02509,  0.13354,  0.54239,  1.55460,  3.08630  &
        , 3.52830,  2.41610,  2.39450,  2.20910,  2.19270,  2.13770  &
        , 2.09600,  2.06980,  2.05380,  2.04010,  2.02570/

        DATA(qext(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        / 0.00413,  0.02585,  0.13693,  0.55257,  1.57340,  3.10260  &
        , 3.51680,  2.40150,  2.39370,  2.21000,  2.19240,  2.13450  &
        , 2.09480,  2.07560,  2.05320,  2.03730,  2.02530/

        DATA(qext(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        / 0.00426,  0.02669,  0.14068,  0.56384,  1.59430,  3.12030  &
        , 3.50320,  2.38940,  2.39040,  2.21310,  2.20050,  2.13800  &
        , 2.10330,  2.07170,  2.05110,  2.03750,  2.02550/

        DATA(qext(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        / 0.00440,  0.02764,  0.14486,  0.57641,  1.61750,  3.13970  &
        , 3.48780,  2.37570,  2.39410,  2.20880,  2.20150,  2.13440  &
        , 2.09920,  2.06930,  2.05110,  2.03550,  2.02350/

        DATA(qext(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        / 0.00456,  0.02870,  0.14956,  0.59054,  1.64350,  3.16050  &
        , 3.47130,  2.36060,  2.39440,  2.21970,  2.19750,  2.13040  &
        , 2.09770,  2.06650,  2.05100,  2.03830,  2.02370/

        DATA(qext(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        / 0.00474,  0.02992,  0.15490,  0.60661,  1.67280,  3.18300  &
        , 3.45240,  2.34370,  2.39860,  2.21830,  2.20080,  2.13920  &
        , 2.09390,  2.06500,  2.05000,  2.03660,  2.02380/

        DATA(qext(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        / 0.00495,  0.03131,  0.16102,  0.62509,  1.70620,  3.20750  &
        , 3.43030,  2.32310,  2.39800,  2.22000,  2.19350,  2.13120  &
        , 2.09140,  2.07120,  2.04770,  2.03450,  2.02320/

        DATA(qext(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        / 0.00519,  0.03294,  0.16814,  0.64667,  1.74460,  3.23450  &
        , 3.40460,  2.29920,  2.39990,  2.23010,  2.18150,  2.13690  &
        , 2.08880,  2.06530,  2.04750,  2.03460,  2.02300/

        DATA(qext(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        / 0.00547,  0.03486,  0.17656,  0.67227,  1.78930,  3.26520  &
        , 3.37210,  2.27500,  2.40140,  2.23280,  2.18590,  2.13280  &
        , 2.09000,  2.06900,  2.04800,  2.03370,  2.02150/

        DATA(qext(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        / 0.00581,  0.03718,  0.18671,  0.70319,  1.84230,  3.30060  &
        , 3.32810,  2.25370,  2.39380,  2.24220,  2.17810,  2.11870  &
        , 2.08930,  2.06380,  2.04640,  2.03340,  2.02180/

        DATA(qext(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        / 0.00623,  0.04004,  0.19923,  0.74133,  1.90610,  3.34050  &
        , 3.27340,  2.23570,  2.39620,  2.24860,  2.16670,  2.11360  &
        , 2.08780,  2.06180,  2.04560,  2.03240,  2.02090/

        DATA(qext(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        / 0.00674,  0.04366,  0.21512,  0.78938,  1.98420,  3.38530  &
        , 3.21630,  2.20870,  2.39290,  2.24580,  2.15810,  2.10800  &
        , 2.08570,  2.06170,  2.04450,  2.03170,  2.01820/

        DATA(qext(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        / 0.00741,  0.04839,  0.23599,  0.85138,  2.08130,  3.43700  &
        , 3.13090,  2.19610,  2.38270,  2.24940,  2.14270,  2.11470  &
        , 2.08570,  2.05900,  2.04290,  2.02980,  2.01730/

        DATA(qext(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        / 0.00831,  0.05484,  0.26458,  0.93364,  2.20610,  3.49060  &
        , 3.01700,  2.19370,  2.35740,  2.23330,  2.13230,  2.11460  &
        , 2.08050,  2.05880,  2.04130,  2.02950,  2.01640/

        DATA(qext(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        / 0.00958,  0.06413,  0.30577,  1.04760,  2.37580,  3.54390  &
        , 2.86490,  2.22020,  2.30260,  2.19730,  2.14710,  2.10040  &
        , 2.07160,  2.05450,  2.03930,  2.02640,  2.01320/

        DATA(qext(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        / 0.01152,  0.07852,  0.36898,  1.21830,  2.61230,  3.55660  &
        , 2.65080,  2.28380,  2.21140,  2.15770,  2.13710,  2.10120  &
        , 2.06750,  2.05040,  2.03600,  2.02490,  2.01130/

        DATA(qext(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        / 0.01481,  0.10356,  0.47826,  1.50800,  2.95860,  3.46980  &
        , 2.37300,  2.36310,  2.19460,  2.17800,  2.12050,  2.09160  &
        , 2.06530,  2.04620,  2.03300,  2.02110,  2.00650/

        DATA(qext(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        / 0.02152,  0.15804,  0.72090,  2.06640,  3.41330,  3.02450  &
        , 2.18440,  2.33350,  2.21490,  2.13210,  2.10390,  2.07600  &
        , 2.05410,  2.03850,  2.02670,  2.01270,  1.99470/

        DATA(qext(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        / 0.00371,  0.00604,  0.00974,  0.01560,  0.02351,  0.03755  &
        , 0.06009,  0.09697,  0.16031,  0.28213,  0.56198,  1.20610  &
        , 2.15550,  2.80300,  2.70990,  2.47520,  2.35740/

        DATA(qext(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        / 0.00379,  0.00615,  0.00992,  0.01591,  0.02399,  0.03833  &
        , 0.06136,  0.09905,  0.16388,  0.28886,  0.57604,  1.22980  &
        , 2.17360,  2.80180,  2.70080,  2.47170,  2.35470/

        DATA(qext(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        / 0.00386,  0.00627,  0.01012,  0.01623,  0.02451,  0.03917  &
        , 0.06271,  0.10128,  0.16772,  0.29612,  0.59121,  1.25500  &
        , 2.19260,  2.80040,  2.69130,  2.46790,  2.35180/

        DATA(qext(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        / 0.00393,  0.00640,  0.01033,  0.01657,  0.02507,  0.04007  &
        , 0.06417,  0.10368,  0.17185,  0.30398,  0.60764,  1.28190  &
        , 2.21280,  2.79890,  2.68150,  2.46380,  2.34880/

        DATA(qext(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        / 0.00401,  0.00653,  0.01055,  0.01694,  0.02567,  0.04104  &
        , 0.06574,  0.10628,  0.17633,  0.31253,  0.62553,  1.31080  &
        , 2.23410,  2.79710,  2.67120,  2.45960,  2.34550/

        DATA(qext(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        / 0.00409,  0.00668,  0.01080,  0.01734,  0.02632,  0.04209  &
        , 0.06745,  0.10909,  0.18122,  0.32191,  0.64512,  1.34190  &
        , 2.25670,  2.79510,  2.66040,  2.45500,  2.34210/

        DATA(qext(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        / 0.00418,  0.00683,  0.01106,  0.01778,  0.02703,  0.04324  &
        , 0.06932,  0.11218,  0.18658,  0.33227,  0.66674,  1.37560  &
        , 2.28070,  2.79270,  2.64900,  2.45000,  2.33840/

        DATA(qext(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        / 0.00428,  0.00700,  0.01135,  0.01826,  0.02781,  0.04450  &
        , 0.07137,  0.11558,  0.19253,  0.34381,  0.69078,  1.41220  &
        , 2.30620,  2.79000,  2.63720,  2.44470,  2.33440/

        DATA(qext(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        / 0.00438,  0.00719,  0.01167,  0.01879,  0.02867,  0.04590  &
        , 0.07365,  0.11937,  0.19917,  0.35680,  0.71778,  1.45240  &
        , 2.33350,  2.78670,  2.62460,  2.43900,  2.33010/

        DATA(qext(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        / 0.00450,  0.00740,  0.01202,  0.01939,  0.02963,  0.04747  &
        , 0.07621,  0.12364,  0.20670,  0.37161,  0.74843,  1.49670  &
        , 2.36260,  2.78250,  2.61140,  2.43260,  2.32540/

        DATA(qext(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        / 0.00463,  0.00763,  0.01243,  0.02006,  0.03072,  0.04925  &
        , 0.07913,  0.12851,  0.21533,  0.38873,  0.78368,  1.54610  &
        , 2.39380,  2.77710,  2.59720,  2.42570,  2.32020/

        DATA(qext(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        / 0.00478,  0.00790,  0.01288,  0.02083,  0.03198,  0.05130  &
        , 0.08249,  0.13415,  0.22539,  0.40889,  0.82481,  1.60160  &
        , 2.42730,  2.76980,  2.58180,  2.41790,  2.31450/

        DATA(qext(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        / 0.00494,  0.00820,  0.01342,  0.02172,  0.03345,  0.05371  &
        , 0.08644,  0.14081,  0.23737,  0.43310,  0.87365,  1.66470  &
        , 2.46360,  2.76010,  2.56510,  2.40910,  2.30790/

        DATA(qext(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        / 0.00513,  0.00856,  0.01404,  0.02279,  0.03520,  0.05659  &
        , 0.09119,  0.14883,  0.25195,  0.46292,  0.93280,  1.73730  &
        , 2.50320,  2.74740,  2.54650,  2.39910,  2.30050/

        DATA(qext(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        / 0.00536,  0.00899,  0.01481,  0.02410,  0.03734,  0.06011  &
        , 0.09702,  0.15878,  0.27025,  0.50078,  1.00610,  1.82200  &
        , 2.54670,  2.73060,  2.52590,  2.38750,  2.29180/

        DATA(qext(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        / 0.00563,  0.00951,  0.01576,  0.02573,  0.04004,  0.06458  &
        , 0.10445,  0.17153,  0.29408,  0.55070,  1.09940,  1.92280  &
        , 2.59380,  2.70790,  2.50280,  2.37380,  2.28150/

        DATA(qext(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        / 0.00597,  0.01019,  0.01699,  0.02788,  0.04358,  0.07047  &
        , 0.11431,  0.18867,  0.32675,  0.61987,  1.22220,  2.04570  &
        , 2.64200,  2.67610,  2.47670,  2.35730,  2.26910/

        DATA(qext(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        / 0.00642,  0.01110,  0.01868,  0.03084,  0.04855,  0.07880  &
        , 0.12839,  0.21359,  0.37563,  0.72384,  1.39280,  2.20050  &
        , 2.68500,  2.63040,  2.44600,  2.33640,  2.25330/

        DATA(qext(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        / 0.00704,  0.01242,  0.02124,  0.03545,  0.05630,  0.09194  &
        , 0.15093,  0.25461,  0.45926,  0.89825,  1.64410,  2.39340  &
        , 2.70740,  2.56430,  2.40790,  2.30840,  2.23210/

        DATA(qext(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        / 0.00800,  0.01467,  0.02584,  0.04412,  0.07141,  0.11840  &
        , 0.19813,  0.34574,  0.65501,  1.26650,  2.06340,  2.62200  &
        , 2.65390,  2.46840,  2.35130,  2.26500,  2.19920/

        DATA(qext(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        / 0.00662,  0.01083,  0.01749,  0.02809,  0.04193,  0.06717  &
        , 0.10817,  0.17717,  0.30271,  0.55931,  1.07340,  1.84510  &
        , 2.52640,  2.67300,  2.44090,  2.33010,  2.24740/

        DATA(qext(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        / 0.00678,  0.01106,  0.01788,  0.02874,  0.04296,  0.06885  &
        , 0.11090,  0.18168,  0.31045,  0.57266,  1.09100,  1.85480  &
        , 2.51940,  2.65780,  2.43530,  2.32630,  2.24480/

        DATA(qext(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        / 0.00692,  0.01131,  0.01830,  0.02942,  0.04407,  0.07064  &
        , 0.11382,  0.18651,  0.31876,  0.58695,  1.10960,  1.86520  &
        , 2.51240,  2.64240,  2.42940,  2.32230,  2.24220/

        DATA(qext(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        / 0.00708,  0.01158,  0.01875,  0.03016,  0.04526,  0.07257  &
        , 0.11695,  0.19171,  0.32770,  0.60228,  1.12930,  1.87630  &
        , 2.50530,  2.62670,  2.42320,  2.31810,  2.23940/

        DATA(qext(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        / 0.00724,  0.01187,  0.01923,  0.03095,  0.04655,  0.07465  &
        , 0.12034,  0.19734,  0.33739,  0.61881,  1.15040,  1.88820  &
        , 2.49820,  2.61070,  2.41670,  2.31370,  2.23640/

        DATA(qext(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        / 0.00742,  0.01218,  0.01976,  0.03181,  0.04794,  0.07691  &
        , 0.12402,  0.20346,  0.34793,  0.63673,  1.17300,  1.90110  &
        , 2.49110,  2.59430,  2.40980,  2.30910,  2.23330/

        DATA(qext(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        / 0.00762,  0.01252,  0.02032,  0.03275,  0.04946,  0.07937  &
        , 0.12804,  0.21015,  0.35948,  0.65626,  1.19730,  1.91500  &
        , 2.48400,  2.57760,  2.40250,  2.30420,  2.23000/

        DATA(qext(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        / 0.00782,  0.01288,  0.02094,  0.03378,  0.05112,  0.08208  &
        , 0.13247,  0.21754,  0.37223,  0.67771,  1.22360,  1.93010  &
        , 2.47690,  2.56030,  2.39470,  2.29910,  2.22650/

        DATA(qext(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        / 0.00805,  0.01328,  0.02163,  0.03492,  0.05297,  0.08509  &
        , 0.13740,  0.22576,  0.38644,  0.70147,  1.25230,  1.94680  &
        , 2.46990,  2.54240,  2.38640,  2.29360,  2.22270/

        DATA(qext(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        / 0.00830,  0.01373,  0.02239,  0.03619,  0.05504,  0.08846  &
        , 0.14292,  0.23501,  0.40246,  0.72802,  1.28380,  1.96510  &
        , 2.46310,  2.52390,  2.37750,  2.28770,  2.21870/

        DATA(qext(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        / 0.00858,  0.01423,  0.02325,  0.03762,  0.05738,  0.09229  &
        , 0.14920,  0.24554,  0.42073,  0.75803,  1.31890,  1.98560  &
        , 2.45640,  2.50480,  2.36810,  2.28140,  2.21430/

        DATA(qext(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        / 0.00889,  0.01479,  0.02422,  0.03926,  0.06007,  0.09669  &
        , 0.15645,  0.25773,  0.44189,  0.79237,  1.35820,  2.00850  &
        , 2.45000,  2.48510,  2.35790,  2.27460,  2.20940/

        DATA(qext(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        / 0.00925,  0.01544,  0.02536,  0.04118,  0.06320,  0.10184  &
        , 0.16495,  0.27207,  0.46682,  0.83228,  1.40300,  2.03460  &
        , 2.44350,  2.46470,  2.34690,  2.26710,  2.20410/

        DATA(qext(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        / 0.00965,  0.01620,  0.02670,  0.04345,  0.06694,  0.10799  &
        , 0.17512,  0.28932,  0.49682,  0.87942,  1.45450,  2.06430  &
        , 2.43700,  2.44340,  2.33500,  2.25880,  2.19820/

        DATA(qext(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        / 0.01013,  0.01711,  0.02831,  0.04622,  0.07149,  0.11551  &
        , 0.18762,  0.31058,  0.53384,  0.93622,  1.51500,  2.09880  &
        , 2.42980,  2.42090,  2.32200,  2.24960,  2.19150/

        DATA(qext(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        / 0.01072,  0.01823,  0.03033,  0.04968,  0.07722,  0.12500  &
        , 0.20345,  0.33768,  0.58094,  1.00630,  1.58700,  2.13900  &
        , 2.42190,  2.39730,  2.30770,  2.23930,  2.18390/

        DATA(qext(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        / 0.01144,  0.01966,  0.03293,  0.05419,  0.08470,  0.13746  &
        , 0.22437,  0.37374,  0.64333,  1.09530,  1.67470,  2.18590  &
        , 2.41250,  2.37250,  2.29170,  2.22740,  2.17500/

        DATA(qext(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        / 0.01238,  0.02157,  0.03647,  0.06039,  0.09511,  0.15496  &
        , 0.25404,  0.42543,  0.73173,  1.21460,  1.78580,  2.24010  &
        , 2.39980,  2.34550,  2.27320,  2.21320,  2.16410/

        DATA(qext(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        / 0.01368,  0.02433,  0.04179,  0.06996,  0.11123,  0.18237  &
        , 0.30113,  0.50833,  0.86910,  1.38640,  1.93270,  2.30080  &
        , 2.37900,  2.31510,  2.25030,  2.19530,  2.15020/

        DATA(qext(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        / 0.01569,  0.02899,  0.05126,  0.08778,  0.14237,  0.23705  &
        , 0.39812,  0.68081,  1.13100,  1.67790,  2.14450,  2.35560  &
        , 2.33920,  2.27490,  2.21730,  2.16870,  2.12920/

        DATA(qext(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        / 0.01346,  0.02088,  0.03271,  0.05151,  0.08608,  0.13671  &
        , 0.21852,  0.35207,  0.56488,  0.86102,  1.23390,  1.64360  &
        , 2.01200,  2.22170,  2.21970,  2.17010,  2.13690/

        DATA(qext(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        / 0.01322,  0.02053,  0.03213,  0.05058,  0.08447,  0.13415  &
        , 0.21451,  0.34601,  0.55666,  0.85245,  1.22920,  1.64620  &
        , 2.02210,  2.23150,  2.22150,  2.17070,  2.13690/

        DATA(qext(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        / 0.01300,  0.02017,  0.03154,  0.04963,  0.08280,  0.13151  &
        , 0.21037,  0.33977,  0.54824,  0.84378,  1.22480,  1.64990  &
        , 2.03380,  2.24210,  2.22330,  2.17150,  2.13690/

        DATA(qext(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        / 0.01277,  0.01979,  0.03093,  0.04864,  0.08107,  0.12877  &
        , 0.20608,  0.33334,  0.53961,  0.83504,  1.22080,  1.65470  &
        , 2.04740,  2.25380,  2.22500,  2.17250,  2.13700/

        DATA(qext(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        / 0.01253,  0.01940,  0.03029,  0.04762,  0.07928,  0.12593  &
        , 0.20164,  0.32668,  0.53076,  0.82625,  1.21740,  1.66110  &
        , 2.06310,  2.26650,  2.22650,  2.17360,  2.13710/

        DATA(qext(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        / 0.01229,  0.01899,  0.02963,  0.04656,  0.07742,  0.12297  &
        , 0.19703,  0.31980,  0.52170,  0.81744,  1.21470,  1.66940  &
        , 2.08140,  2.28030,  2.22780,  2.17480,  2.13720/

        DATA(qext(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        / 0.01203,  0.01857,  0.02895,  0.04544,  0.07548,  0.11990  &
        , 0.19223,  0.31267,  0.51240,  0.80868,  1.21300,  1.67990  &
        , 2.10270,  2.29530,  2.22890,  2.17620,  2.13740/

        DATA(qext(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        / 0.01177,  0.01813,  0.02822,  0.04428,  0.07345,  0.11668  &
        , 0.18722,  0.30526,  0.50288,  0.80004,  1.21240,  1.69340  &
        , 2.12760,  2.31170,  2.22980,  2.17770,  2.13760/

        DATA(qext(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        / 0.01149,  0.01767,  0.02747,  0.04306,  0.07132,  0.11330  &
        , 0.18198,  0.29754,  0.49313,  0.79166,  1.21340,  1.71050  &
        , 2.15680,  2.32910,  2.23030,  2.17920,  2.13770/

        DATA(qext(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        / 0.01120,  0.01718,  0.02667,  0.04178,  0.06907,  0.10974  &
        , 0.17646,  0.28950,  0.48320,  0.78372,  1.21660,  1.73230  &
        , 2.19130,  2.34740,  2.23040,  2.18070,  2.13790/

        DATA(qext(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        / 0.01089,  0.01667,  0.02583,  0.04042,  0.06670,  0.10598  &
        , 0.17065,  0.28110,  0.47315,  0.77654,  1.22280,  1.76050  &
        , 2.23220,  2.36620,  2.23020,  2.18190,  2.13790/

        DATA(qext(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        / 0.01057,  0.01612,  0.02494,  0.03897,  0.06417,  0.10198  &
        , 0.16450,  0.27235,  0.46311,  0.77064,  1.23320,  1.79690  &
        , 2.28130,  2.38500,  2.22980,  2.18260,  2.13780/

        DATA(qext(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        / 0.01022,  0.01554,  0.02398,  0.03742,  0.06147,  0.09772  &
        , 0.15799,  0.26327,  0.45333,  0.76687,  1.24970,  1.84460  &
        , 2.34080,  2.40160,  2.22970,  2.18260,  2.13750/

        DATA(qext(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        / 0.00985,  0.01491,  0.02295,  0.03576,  0.05857,  0.09317  &
        , 0.15110,  0.25395,  0.44429,  0.76673,  1.27540,  1.90780  &
        , 2.41190,  2.41310,  2.23070,  2.18170,  2.13700/

        DATA(qext(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        / 0.00945,  0.01424,  0.02185,  0.03398,  0.05547,  0.08832  &
        , 0.14386,  0.24464,  0.43687,  0.77290,  1.31580,  1.99350  &
        , 2.49510,  2.41490,  2.23410,  2.18050,  2.13620/

        DATA(qext(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        / 0.00901,  0.01350,  0.02065,  0.03207,  0.05216,  0.08319  &
        , 0.13639,  0.23589,  0.43279,  0.79046,  1.38010,  2.11290  &
        , 2.59160,  2.39880,  2.24160,  2.17960,  2.13450/

        DATA(qext(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        / 0.00853,  0.01270,  0.01937,  0.03005,  0.04866,  0.07788  &
        , 0.12904,  0.22896,  0.43550,  0.82935,  1.48410,  2.27840  &
        , 2.68730,  2.35520,  2.25150,  2.17880,  2.13120/

        DATA(qext(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        / 0.00800,  0.01184,  0.01802,  0.02796,  0.04509,  0.07263  &
        , 0.12264,  0.22695,  0.45277,  0.91014,  1.65580,  2.51070  &
        , 2.75170,  2.28400,  2.24970,  2.17240,  2.12640/

        DATA(qext(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        / 0.00742,  0.01093,  0.01665,  0.02595,  0.04177,  0.06835  &
        , 0.12001,  0.23916,  0.50645,  1.07370,  1.96710,  2.82060  &
        , 2.68290,  2.24560,  2.21250,  2.16010,  2.11780/

        DATA(qext(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        / 0.00679,  0.01006,  0.01557,  0.02475,  0.04025,  0.06909  &
        , 0.13484,  0.30639,  0.70846,  1.51040,  2.59390,  3.05600  &
        , 2.28920,  2.31640,  2.19940,  2.13920,  2.10200/

        DATA(qext(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        / 0.01070,  0.01682,  0.02658,  0.04213,  0.06873,  0.11075  &
        , 0.18330,  0.32338,  0.64074,  1.34360,  2.38250,  3.07040  &
        , 2.81750,  2.42450,  2.36310,  2.26590,  2.19850/

        DATA(qext(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        / 0.01062,  0.01670,  0.02637,  0.04179,  0.06820,  0.10988  &
        , 0.18187,  0.32090,  0.63566,  1.32980,  2.35900,  3.05560  &
        , 2.81850,  2.42020,  2.36160,  2.26400,  2.19670/

        DATA(qext(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        / 0.01054,  0.01657,  0.02615,  0.04142,  0.06762,  0.10894  &
        , 0.18032,  0.31823,  0.63026,  1.31540,  2.33450,  3.03990  &
        , 2.81960,  2.41620,  2.35990,  2.26200,  2.19490/

        DATA(qext(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        / 0.01046,  0.01643,  0.02591,  0.04102,  0.06699,  0.10791  &
        , 0.17863,  0.31534,  0.62451,  1.30030,  2.30900,  3.02340  &
        , 2.82070,  2.41220,  2.35800,  2.26000,  2.19300/

        DATA(qext(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        / 0.01037,  0.01627,  0.02565,  0.04059,  0.06630,  0.10679  &
        , 0.17679,  0.31222,  0.61837,  1.28440,  2.28230,  3.00600  &
        , 2.82200,  2.40820,  2.35580,  2.25770,  2.19090/

        DATA(qext(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        / 0.01027,  0.01610,  0.02536,  0.04012,  0.06555,  0.10557  &
        , 0.17478,  0.30884,  0.61182,  1.26770,  2.25460,  2.98760  &
        , 2.82330,  2.40400,  2.35320,  2.25530,  2.18870/

        DATA(qext(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        / 0.01017,  0.01592,  0.02505,  0.03962,  0.06473,  0.10423  &
        , 0.17260,  0.30517,  0.60482,  1.25020,  2.22550,  2.96810  &
        , 2.82480,  2.39920,  2.35030,  2.25270,  2.18640/

        DATA(qext(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        / 0.01006,  0.01572,  0.02472,  0.03906,  0.06382,  0.10275  &
        , 0.17020,  0.30119,  0.59732,  1.23160,  2.19510,  2.94750  &
        , 2.82630,  2.39370,  2.34710,  2.24980,  2.18390/

        DATA(qext(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        / 0.00993,  0.01551,  0.02435,  0.03845,  0.06282,  0.10113  &
        , 0.16756,  0.29685,  0.58926,  1.21210,  2.16320,  2.92560  &
        , 2.82780,  2.38760,  2.34360,  2.24660,  2.18120/

        DATA(qext(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        / 0.00980,  0.01527,  0.02394,  0.03778,  0.06171,  0.09933  &
        , 0.16465,  0.29210,  0.58061,  1.19150,  2.12970,  2.90240  &
        , 2.82890,  2.38140,  2.33940,  2.24300,  2.17830/

        DATA(qext(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        / 0.00965,  0.01500,  0.02350,  0.03704,  0.06047,  0.09733  &
        , 0.16142,  0.28690,  0.57132,  1.16980,  2.09470,  2.87770  &
        , 2.82920,  2.37590,  2.33440,  2.23890,  2.17510/

        DATA(qext(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        / 0.00948,  0.01471,  0.02300,  0.03622,  0.05909,  0.09509  &
        , 0.15784,  0.28121,  0.56135,  1.14700,  2.05830,  2.85170  &
        , 2.82800,  2.37150,  2.32870,  2.23440,  2.17170/

        DATA(qext(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        / 0.00930,  0.01438,  0.02244,  0.03530,  0.05754,  0.09258  &
        , 0.15385,  0.27498,  0.55073,  1.12330,  2.02080,  2.82450  &
        , 2.82460,  2.36650,  2.32240,  2.22920,  2.16780/

        DATA(qext(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        / 0.00909,  0.01401,  0.02182,  0.03427,  0.05578,  0.08977  &
        , 0.14941,  0.26822,  0.53958,  1.09920,  1.98330,  2.79670  &
        , 2.81920,  2.35840,  2.31470,  2.22330,  2.16360/

        DATA(qext(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        / 0.00886,  0.01359,  0.02111,  0.03311,  0.05380,  0.08661  &
        , 0.14449,  0.26100,  0.52825,  1.07600,  1.94770,  2.76970  &
        , 2.81300,  2.35000,  2.30620,  2.21650,  2.15870/

        DATA(qext(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        / 0.00859,  0.01312,  0.02032,  0.03182,  0.05158,  0.08309  &
        , 0.13914,  0.25362,  0.51755,  1.05610,  1.91800,  2.74750  &
        , 2.80670,  2.33960,  2.29610,  2.20880,  2.15300/

        DATA(qext(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        / 0.00828,  0.01258,  0.01943,  0.03039,  0.04912,  0.07927  &
        , 0.13356,  0.24684,  0.50936,  1.04480,  1.90210,  2.73730  &
        , 2.79430,  2.32440,  2.28380,  2.19970,  2.14620/

        DATA(qext(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        / 0.00792,  0.01197,  0.01845,  0.02886,  0.04651,  0.07534  &
        , 0.12835,  0.24265,  0.50811,  1.05370,  1.91470,  2.75120  &
        , 2.77160,  2.30170,  2.26740,  2.18840,  2.13770/

        DATA(qext(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        / 0.00751,  0.01131,  0.01745,  0.02740,  0.04406,  0.07208  &
        , 0.12565,  0.24734,  0.52697,  1.10910,  1.99660,  2.81460  &
        , 2.71050,  2.27250,  2.24140,  2.17260,  2.12630/

        DATA(qext(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        / 0.00705,  0.01068,  0.01674,  0.02675,  0.04332,  0.07310  &
        , 0.13566,  0.28890,  0.63118,  1.31960,  2.29670,  2.94500  &
        , 2.48900,  2.27270,  2.20000,  2.14740,  2.10820/

        DATA(qext(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        / 0.00235,  0.00373,  0.00592,  0.00945,  0.01527,  0.02557  &
        , 0.04786,  0.11442,  0.32899,  0.93060,  2.06240,  3.20160  &
        , 2.76950,  2.40340,  2.31150,  2.19560,  2.14520/

        DATA(qext(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        / 0.00235,  0.00372,  0.00592,  0.00945,  0.01527,  0.02560  &
        , 0.04816,  0.11609,  0.33443,  0.94474,  2.08670,  3.21330  &
        , 2.75590,  2.40410,  2.30920,  2.19490,  2.14400/

        DATA(qext(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        / 0.00234,  0.00372,  0.00592,  0.00945,  0.01527,  0.02565  &
        , 0.04851,  0.11795,  0.34041,  0.96009,  2.11300,  3.22650  &
        , 2.73930,  2.40370,  2.30590,  2.19400,  2.14280/

        DATA(qext(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        / 0.00234,  0.00372,  0.00592,  0.00945,  0.01527,  0.02571  &
        , 0.04891,  0.12004,  0.34704,  0.97683,  2.14140,  3.24110  &
        , 2.71880,  2.40080,  2.30320,  2.19360,  2.14160/

        DATA(qext(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        / 0.00234,  0.00372,  0.00592,  0.00945,  0.01528,  0.02578  &
        , 0.04938,  0.12241,  0.35442,  0.99523,  2.17200,  3.25690  &
        , 2.69460,  2.39880,  2.29860,  2.19300,  2.14050/

        DATA(qext(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        / 0.00234,  0.00372,  0.00592,  0.00946,  0.01529,  0.02588  &
        , 0.04992,  0.12509,  0.36271,  1.01560,  2.20510,  3.27260  &
        , 2.66850,  2.39890,  2.29190,  2.19250,  2.13950/

        DATA(qext(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        / 0.00234,  0.00372,  0.00593,  0.00948,  0.01532,  0.02600  &
        , 0.05056,  0.12818,  0.37212,  1.03850,  2.24100,  3.28710  &
        , 2.64420,  2.39780,  2.28540,  2.19180,  2.13860/

        DATA(qext(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        / 0.00234,  0.00373,  0.00593,  0.00949,  0.01535,  0.02615  &
        , 0.05132,  0.13176,  0.38293,  1.06440,  2.28010,  3.29890  &
        , 2.62070,  2.39130,  2.27580,  2.19050,  2.13760/

        DATA(qext(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        / 0.00234,  0.00373,  0.00595,  0.00952,  0.01540,  0.02634  &
        , 0.05223,  0.13596,  0.39549,  1.09430,  2.32330,  3.30760  &
        , 2.59190,  2.38650,  2.26500,  2.18840,  2.13650/

        DATA(qext(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        / 0.00235,  0.00374,  0.00596,  0.00956,  0.01547,  0.02658  &
        , 0.05334,  0.14097,  0.41036,  1.12930,  2.37210,  3.31470  &
        , 2.55380,  2.38130,  2.25260,  2.18510,  2.13460/

        DATA(qext(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        / 0.00235,  0.00375,  0.00599,  0.00960,  0.01556,  0.02689  &
        , 0.05473,  0.14704,  0.42830,  1.17110,  2.42880,  3.32350  &
        , 2.51440,  2.37030,  2.23830,  2.18020,  2.13170/

        DATA(qext(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        / 0.00236,  0.00377,  0.00602,  0.00967,  0.01568,  0.02729  &
        , 0.05649,  0.15453,  0.45045,  1.22210,  2.49670,  3.33420  &
        , 2.47860,  2.35840,  2.22380,  2.17410,  2.12810/

        DATA(qext(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        / 0.00237,  0.00379,  0.00606,  0.00976,  0.01585,  0.02783  &
        , 0.05879,  0.16401,  0.47858,  1.28540,  2.57910,  3.33490  &
        , 2.43220,  2.34280,  2.21240,  2.16840,  2.12510/

        DATA(qext(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        / 0.00238,  0.00382,  0.00612,  0.00988,  0.01609,  0.02857  &
        , 0.06187,  0.17631,  0.51544,  1.36480,  2.67610,  3.31470  &
        , 2.39010,  2.32340,  2.20740,  2.16580,  2.12360/

        DATA(qext(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        / 0.00240,  0.00386,  0.00621,  0.01005,  0.01643,  0.02962  &
        , 0.06618,  0.19283,  0.56548,  1.46490,  2.78560,  3.28670  &
        , 2.34810,  2.30200,  2.21280,  2.16740,  2.12100/

        DATA(qext(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        / 0.00243,  0.00392,  0.00634,  0.01030,  0.01694,  0.03117  &
        , 0.07250,  0.21592,  0.63566,  1.59230,  2.91590,  3.22030  &
        , 2.31660,  2.28570,  2.22680,  2.16380,  2.11500/

        DATA(qext(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        / 0.00246,  0.00401,  0.00653,  0.01069,  0.01774,  0.03361  &
        , 0.08240,  0.24990,  0.73660,  1.76420,  3.08010,  3.10240  &
        , 2.31300,  2.28550,  2.22580,  2.14680,  2.11030/

        DATA(qext(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        / 0.00252,  0.00416,  0.00683,  0.01130,  0.01907,  0.03785  &
        , 0.09961,  0.30467,  0.88741,  2.02580,  3.25040,  2.88680  &
        , 2.34170,  2.29750,  2.18910,  2.14370,  2.10350/

        DATA(qext(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        / 0.00262,  0.00440,  0.00736,  0.01243,  0.02165,  0.04658  &
        , 0.13456,  0.40898,  1.14490,  2.41730,  3.38420,  2.55430  &
        , 2.35810,  2.24140,  2.17890,  2.12900,  2.09420/

        DATA(qext(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        / 0.00280,  0.00489,  0.00848,  0.01501,  0.02845,  0.07264  &
        , 0.23193,  0.70597,  1.73200,  3.07620,  3.14340,  2.29680  &
        , 2.28220,  2.22450,  2.14580,  2.10950,  2.07990/

! Qext for Sea Salt:
        DATA(qext(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00349, 0.00585, 0.01063, 0.02441, 0.07376, 0.26299  &
        ,0.84952, 2.05200, 3.35470, 2.85810, 2.47760, 2.33190  &
        ,2.20670, 2.15340, 2.10590, 2.08120, 2.05990/

        DATA(qext(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00357, 0.00599, 0.01089, 0.02518, 0.07646, 0.27096  &
        ,0.86962, 2.08480, 3.36910, 2.82990, 2.48150, 2.33470  &
        ,2.21820, 2.14500, 2.11100, 2.07930, 2.06000/

        DATA(qext(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00365, 0.00611, 0.01117, 0.02600, 0.07936, 0.27951  &
        ,0.89098, 2.11900, 3.38140, 2.80200, 2.47430, 2.32240  &
        ,2.21670, 2.14250, 2.10880, 2.08250, 2.05990/

        DATA(qext(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00372, 0.00625, 0.01146, 0.02689, 0.08251, 0.28878  &
        ,0.91393, 2.15500, 3.39170, 2.76950, 2.48360, 2.31850  &
        ,2.21360, 2.13720, 2.11040, 2.08130, 2.05790/

        DATA(qext(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00380, 0.00639, 0.01177, 0.02786, 0.08596, 0.29893  &
        ,0.93888, 2.19320, 3.40050, 2.73050, 2.48410, 2.31380  &
        ,2.23340, 2.13990, 2.11120, 2.07980, 2.05830/

        DATA(qext(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00389, 0.00655, 0.01211, 0.02894, 0.08979, 0.31019  &
        ,0.96636, 2.23450, 3.40910, 2.68780, 2.47120, 2.30470  &
        ,2.22760, 2.14050, 2.10520, 2.07710, 2.05760/

        DATA(qext(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00397, 0.00671, 0.01248, 0.03014, 0.09408, 0.32282  &
        ,0.99705, 2.27990, 3.41940, 2.65090, 2.46890, 2.30380  &
        ,2.22180, 2.14070, 2.10180, 2.07730, 2.05520/

        DATA(qext(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00407, 0.00689, 0.01290, 0.03150, 0.09895, 0.33721  &
        ,1.03180, 2.33070, 3.43100, 2.61480, 2.46280, 2.28200  &
        ,2.21590, 2.14180, 2.10560, 2.07630, 2.05630/

        DATA(qext(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00418, 0.00709, 0.01336, 0.03307, 0.10453, 0.35383  &
        ,1.07160, 2.38810, 3.43980, 2.56950, 2.44850, 2.27090  &
        ,2.20540, 2.13770, 2.10780, 2.07520, 2.05480/

        DATA(qext(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00429, 0.00731, 0.01389, 0.03490, 0.11102, 0.37334  &
        ,1.11760, 2.45290, 3.44060, 2.52130, 2.42610, 2.26750  &
        ,2.19970, 2.13440, 2.10170, 2.07460, 2.05400/

        DATA(qext(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00442, 0.00757, 0.01451, 0.03707, 0.11868, 0.39663  &
        ,1.17140, 2.52490, 3.43240, 2.48520, 2.39120, 2.23590  &
        ,2.16620, 2.12900, 2.09940, 2.07080, 2.05470/

        DATA(qext(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00457, 0.00786, 0.01524, 0.03968, 0.12784, 0.42494  &
        ,1.23440, 2.60350, 3.42150, 2.43800, 2.35710, 2.22370  &
        ,2.17250, 2.13740, 2.09600, 2.07170, 2.05300/

        DATA(qext(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00473, 0.00819, 0.01611, 0.04291, 0.13902, 0.45998  &
        ,1.30850, 2.68960, 3.40300, 2.40630, 2.31340, 2.21730  &
        ,2.16670, 2.13510, 2.09370, 2.06860, 2.05170/

        DATA(qext(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00493, 0.00859, 0.01718, 0.04701, 0.15294, 0.50417  &
        ,1.39610, 2.78890, 3.35590, 2.37480, 2.28100, 2.21960  &
        ,2.17270, 2.12740, 2.09420, 2.06590, 2.05000/

        DATA(qext(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00516, 0.00908, 0.01855, 0.05238, 0.17077, 0.56086  &
        ,1.50230, 2.90860, 3.29510, 2.36550, 2.23540, 2.23320  &
        ,2.16420, 2.11760, 2.08870, 2.06420, 2.04830/

        DATA(qext(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00544, 0.00970, 0.02037, 0.05975, 0.19450, 0.63487  &
        ,1.63810, 3.03920, 3.18750, 2.38200, 2.23630, 2.23270  &
        ,2.15860, 2.11720, 2.08660, 2.06270, 2.04640/

        DATA(qext(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00579, 0.01051, 0.02294, 0.07055, 0.22793, 0.73404  &
        ,1.82270, 3.18160, 3.02990, 2.41900, 2.28590, 2.20170  &
        ,2.15960, 2.11370, 2.08370, 2.05950, 2.04270/

        DATA(qext(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00626, 0.01167, 0.02696, 0.08791, 0.27982, 0.87708  &
        ,2.07510, 3.32470, 2.79050, 2.45920, 2.32460, 2.21040  &
        ,2.14170, 2.10430, 2.07820, 2.05680, 2.04480/

        DATA(qext(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00695, 0.01356, 0.03432, 0.12030, 0.37724, 1.12820  &
        ,2.46260, 3.40740, 2.49130, 2.38910, 2.24180, 2.17780  &
        ,2.12330, 2.09340, 2.07070, 2.05270, 2.04050/

        DATA(qext(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00814, 0.01752, 0.05311, 0.20048, 0.63541, 1.67050  &
        ,3.07500, 3.11800, 2.39820, 2.25790, 2.21080, 2.15350  &
        ,2.11390, 2.08310, 2.06090, 2.04480, 2.03340/

        DATA(qext(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00073, 0.00478, 0.02966, 0.15050, 0.53722, 1.54670  &
        ,3.05970, 3.42670, 2.44540, 2.26250, 2.22890, 2.15640  &
        ,2.12460, 2.08730, 2.06720, 2.04990, 2.03650/

        DATA(qext(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00077, 0.00503, 0.03101, 0.15593, 0.55377, 1.57760  &
        ,3.09020, 3.40500, 2.43890, 2.27260, 2.22590, 2.15460  &
        ,2.12160, 2.09360, 2.06470, 2.04870, 2.03610/

        DATA(qext(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00081, 0.00527, 0.03248, 0.16174, 0.57138, 1.61060  &
        ,3.12080, 3.37910, 2.43450, 2.27830, 2.21480, 2.16280  &
        ,2.12660, 2.09620, 2.06490, 2.04830, 2.03510/

        DATA(qext(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00085, 0.00555, 0.03409, 0.16803, 0.59022, 1.64620  &
        ,3.15150, 3.34960, 2.43090, 2.29100, 2.19980, 2.15830  &
        ,2.11930, 2.08830, 2.06620, 2.04760, 2.03490/

        DATA(qext(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00090, 0.00585, 0.03587, 0.17489, 0.61054, 1.68490  &
        ,3.18270, 3.31850, 2.42300, 2.30130, 2.19480, 2.16140  &
        ,2.12130, 2.08860, 2.06390, 2.04680, 2.03520/

        DATA(qext(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00095, 0.00618, 0.03785, 0.18244, 0.63264, 1.72720  &
        ,3.21510, 3.28760, 2.42060, 2.30780, 2.19060, 2.16150  &
        ,2.11790, 2.08930, 2.06310, 2.04700, 2.03460/

        DATA(qext(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00100, 0.00656, 0.04009, 0.19087, 0.65691, 1.77390  &
        ,3.25010, 3.25270, 2.41610, 2.31980, 2.18760, 2.16470  &
        ,2.11480, 2.08460, 2.06290, 2.04650, 2.03430/

        DATA(qext(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00107, 0.00700, 0.04265, 0.20038, 0.68382, 1.82560  &
        ,3.28870, 3.20910, 2.40750, 2.33030, 2.18670, 2.16050  &
        ,2.11590, 2.08410, 2.06090, 2.04530, 2.03320/

        DATA(qext(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00114, 0.00751, 0.04562, 0.21127, 0.71403, 1.88290  &
        ,3.33030, 3.15620, 2.40130, 2.33820, 2.18590, 2.15390  &
        ,2.11220, 2.08320, 2.06180, 2.04550, 2.03290/

        DATA(qext(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00123, 0.00810, 0.04911, 0.22392, 0.74842, 1.94670  &
        ,3.37220, 3.10120, 2.39480, 2.35070, 2.18800, 2.15340  &
        ,2.11070, 2.08280, 2.06080, 2.04420, 2.03280/

        DATA(qext(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00133, 0.00882, 0.05327, 0.23887, 0.78821, 2.01810  &
        ,3.41150, 3.03950, 2.38500, 2.34780, 2.19880, 2.15390  &
        ,2.11260, 2.08450, 2.06070, 2.04320, 2.03230/

        DATA(qext(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00145, 0.00971, 0.05833, 0.25688, 0.83522, 2.09940  &
        ,3.44850, 2.96190, 2.37200, 2.35100, 2.20660, 2.15880  &
        ,2.11170, 2.08000, 2.05730, 2.04270, 2.03170/

        DATA(qext(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00161, 0.01081, 0.06459, 0.27911, 0.89214, 2.19430  &
        ,3.48840, 2.88390, 2.35320, 2.33630, 2.21150, 2.15930  &
        ,2.10720, 2.07590, 2.05600, 2.04210, 2.03100/

        DATA(qext(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00180, 0.01223, 0.07253, 0.30737, 0.96298, 2.30810  &
        ,3.52580, 2.78820, 2.32890, 2.31120, 2.20410, 2.14200  &
        ,2.10330, 2.07500, 2.05530, 2.04220, 2.03110/

        DATA(qext(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00206, 0.01413, 0.08288, 0.34459, 1.05340, 2.44450  &
        ,3.54420, 2.68990, 2.29220, 2.28850, 2.19630, 2.13620  &
        ,2.10310, 2.07400, 2.05410, 2.03940, 2.02820/

        DATA(qext(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00242, 0.01676, 0.09691, 0.39573, 1.17090, 2.60460  &
        ,3.54950, 2.58520, 2.25350, 2.27360, 2.19760, 2.13130  &
        ,2.09650, 2.06900, 2.05110, 2.03880, 2.02800/

        DATA(qext(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00293, 0.02066, 0.11694, 0.46932, 1.32710, 2.80500  &
        ,3.50020, 2.49110, 2.22740, 2.26050, 2.17740, 2.12660  &
        ,2.09100, 2.06630, 2.04970, 2.03680, 2.02690/

        DATA(qext(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00374, 0.02698, 0.14786, 0.58054, 1.55110, 3.05690  &
        ,3.36840, 2.42170, 2.26450, 2.20950, 2.14900, 2.11740  &
        ,2.08720, 2.06380, 2.04620, 2.03400, 2.02590/

        DATA(qext(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00520, 0.03881, 0.20254, 0.76450, 1.91460, 3.35310  &
        ,3.06610, 2.38510, 2.34530, 2.19040, 2.15060, 2.10930  &
        ,2.08170, 2.05850, 2.04290, 2.03050, 2.02300/

        DATA(qext(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00854, 0.06789, 0.33643, 1.18010, 2.58590, 3.53320  &
        ,2.54050, 2.23250, 2.26170, 2.19240, 2.12710, 2.09350  &
        ,2.06850, 2.05000, 2.03640, 2.02810, 2.02020/

        DATA(qext(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.01390, 0.07864, 0.33227, 1.04460, 2.28450, 3.49860  &
        ,2.98590, 2.19520, 2.34920, 2.22500, 2.14170, 2.11100  &
        ,2.08200, 2.06070, 2.04250, 2.03250, 2.02220/

        DATA(qext(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.01464, 0.08201, 0.34275, 1.06840, 2.31620, 3.50960  &
        ,2.95390, 2.20050, 2.34220, 2.22620, 2.13860, 2.10850  &
        ,2.08410, 2.05850, 2.04410, 2.03120, 2.02240/

        DATA(qext(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.01537, 0.08541, 0.35400, 1.09360, 2.34930, 3.52010  &
        ,2.92280, 2.20660, 2.32180, 2.21840, 2.14490, 2.10860  &
        ,2.07680, 2.05670, 2.04130, 2.03020, 2.02270/

        DATA(qext(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.01612, 0.08908, 0.36614, 1.12070, 2.38480, 3.53090  &
        ,2.89170, 2.21720, 2.31400, 2.20760, 2.14590, 2.10120  &
        ,2.07690, 2.05470, 2.04240, 2.03040, 2.02310/

        DATA(qext(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.01694, 0.09310, 0.37934, 1.14990, 2.42310, 3.53880  &
        ,2.85880, 2.22360, 2.29950, 2.20100, 2.14770, 2.09830  &
        ,2.07260, 2.05660, 2.04210, 2.02990, 2.02290/

        DATA(qext(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.01784, 0.09754, 0.39382, 1.18200, 2.46470, 3.54480  &
        ,2.82240, 2.23140, 2.28320, 2.19710, 2.14640, 2.10400  &
        ,2.07280, 2.05520, 2.04110, 2.02970, 2.02150/

        DATA(qext(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.01885, 0.10249, 0.40984, 1.21740, 2.50980, 3.54840  &
        ,2.78170, 2.24250, 2.26770, 2.18710, 2.15200, 2.09650  &
        ,2.07650, 2.05450, 2.03980, 2.02960, 2.02250/

        DATA(qext(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.01999, 0.10807, 0.42771, 1.25700, 2.55930, 3.54920  &
        ,2.73530, 2.25490, 2.25220, 2.18260, 2.14880, 2.10060  &
        ,2.07430, 2.05380, 2.03970, 2.03000, 2.02120/

        DATA(qext(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.02130, 0.11445, 0.44784, 1.30140, 2.61360, 3.55050  &
        ,2.68520, 2.27430, 2.23320, 2.16970, 2.14060, 2.10390  &
        ,2.07340, 2.05340, 2.03930, 2.02810, 2.02020/

        DATA(qext(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.02282, 0.12184, 0.47072, 1.35140, 2.67260, 3.54930  &
        ,2.63090, 2.28640, 2.21800, 2.16130, 2.12830, 2.09900  &
        ,2.07080, 2.05160, 2.03850, 2.02900, 2.02030/

        DATA(qext(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.02461, 0.13049, 0.49706, 1.40850, 2.73990, 3.54170  &
        ,2.57320, 2.30470, 2.20080, 2.15640, 2.13060, 2.09570  &
        ,2.06830, 2.05220, 2.03740, 2.02830, 2.02120/

        DATA(qext(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.02676, 0.14077, 0.52789, 1.47490, 2.81590, 3.52420  &
        ,2.51200, 2.32930, 2.18840, 2.16050, 2.12120, 2.09040  &
        ,2.06970, 2.05110, 2.03810, 2.02760, 2.01930/

        DATA(qext(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.02940, 0.15318, 0.56477, 1.55310, 2.90210, 3.49850  &
        ,2.44450, 2.34220, 2.18330, 2.16790, 2.11720, 2.08720  &
        ,2.06750, 2.05040, 2.03660, 2.02660, 2.02010/

        DATA(qext(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.03270, 0.16840, 0.61014, 1.64650, 2.99490, 3.45990  &
        ,2.36890, 2.36290, 2.18670, 2.17930, 2.12060, 2.09160  &
        ,2.06580, 2.04780, 2.03540, 2.02620, 2.01870/

        DATA(qext(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.03695, 0.18755, 0.66800, 1.75840, 3.09880, 3.39210  &
        ,2.30260, 2.37550, 2.20480, 2.17800, 2.12780, 2.08800  &
        ,2.06190, 2.04560, 2.03420, 2.02420, 2.01850/

        DATA(qext(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.04262, 0.21252, 0.74478, 1.89610, 3.21790, 3.30210  &
        ,2.23270, 2.37720, 2.22460, 2.16150, 2.11120, 2.08180  &
        ,2.06340, 2.04430, 2.03370, 2.02530, 2.01850/

        DATA(qext(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.05057, 0.24710, 0.84932, 2.07510, 3.34410, 3.16020  &
        ,2.18690, 2.36780, 2.23490, 2.14050, 2.10870, 2.07700  &
        ,2.05890, 2.04380, 2.03110, 2.02410, 2.01650/

        DATA(qext(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.06256, 0.29951, 0.99765, 2.31660, 3.47200, 2.94570  &
        ,2.19220, 2.32250, 2.21080, 2.13720, 2.10380, 2.07590  &
        ,2.05420, 2.04030, 2.03010, 2.02140, 2.01610/

        DATA(qext(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.08290, 0.38942, 1.23930, 2.67350, 3.53380, 2.61790  &
        ,2.28460, 2.20780, 2.15350, 2.12660, 2.09820, 2.06710  &
        ,2.04930, 2.03730, 2.02680, 2.01970, 2.01480/

        DATA(qext(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.12530, 0.57315, 1.71420, 3.21310, 3.30290, 2.21790  &
        ,2.38380, 2.23300, 2.15730, 2.10010, 2.08190, 2.05830  &
        ,2.04370, 2.03130, 2.02260, 2.01610, 2.01320/

        DATA(qext(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00639, 0.01033, 0.01658, 0.02650, 0.04061, 0.06489  &
        ,0.10424, 0.17010, 0.28871, 0.52953, 1.02380, 1.76560  &
        ,2.43410, 2.64750, 2.49550, 2.36480, 2.27530/

        DATA(qext(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00652, 0.01052, 0.01687, 0.02697, 0.04136, 0.06610  &
        ,0.10622, 0.17347, 0.29498, 0.54260, 1.04820, 1.79500  &
        ,2.45320, 2.64710, 2.49060, 2.36200, 2.27320/

        DATA(qext(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00664, 0.01071, 0.01718, 0.02747, 0.04214, 0.06736  &
        ,0.10829, 0.17701, 0.30162, 0.55646, 1.07370, 1.82530  &
        ,2.47220, 2.64620, 2.48570, 2.35910, 2.27100/

        DATA(qext(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00676, 0.01091, 0.01750, 0.02798, 0.04297, 0.06869  &
        ,0.11047, 0.18076, 0.30869, 0.57126, 1.10070, 1.85660  &
        ,2.49110, 2.64460, 2.48050, 2.35600, 2.26870/

        DATA(qext(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00689, 0.01112, 0.01784, 0.02854, 0.04384, 0.07011  &
        ,0.11281, 0.18477, 0.31628, 0.58721, 1.12950, 1.88930  &
        ,2.51000, 2.64230, 2.47520, 2.35290, 2.26630/

        DATA(qext(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00702, 0.01134, 0.01821, 0.02913, 0.04478, 0.07163  &
        ,0.11532, 0.18911, 0.32454, 0.60458, 1.16030, 1.92360  &
        ,2.52880, 2.63920, 2.46960, 2.34950, 2.26370/

        DATA(qext(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00716, 0.01158, 0.01860, 0.02977, 0.04580, 0.07328  &
        ,0.11804, 0.19384, 0.33359, 0.62366, 1.19380, 1.95990  &
        ,2.54750, 2.63530, 2.46380, 2.34590, 2.26100/

        DATA(qext(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00731, 0.01183, 0.01902, 0.03046, 0.04691, 0.07509  &
        ,0.12104, 0.19905, 0.34363, 0.64486, 1.23020, 1.99850  &
        ,2.56620, 2.63030, 2.45760, 2.34200, 2.25800/

        DATA(qext(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00748, 0.01212, 0.01949, 0.03123, 0.04814, 0.07709  &
        ,0.12436, 0.20485, 0.35490, 0.66865, 1.27040, 2.03990  &
        ,2.58470, 2.62410, 2.45100, 2.33780, 2.25480/

        DATA(qext(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00766, 0.01243, 0.02002, 0.03209, 0.04951, 0.07933  &
        ,0.12808, 0.21139, 0.36770, 0.69565, 1.31500, 2.08430  &
        ,2.60300, 2.61660, 2.44400, 2.33320, 2.25130/

        DATA(qext(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00787, 0.01278, 0.02061, 0.03306, 0.05106, 0.08186  &
        ,0.13231, 0.21886, 0.38243, 0.72666, 1.36490, 2.13220  &
        ,2.62080, 2.60740, 2.43640, 2.32810, 2.24740/

        DATA(qext(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00810, 0.01318, 0.02128, 0.03417, 0.05283, 0.08476  &
        ,0.13716, 0.22749, 0.39963, 0.76273, 1.42120, 2.18400  &
        ,2.63770, 2.59620, 2.42810, 2.32250, 2.24320/

        DATA(qext(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00836, 0.01364, 0.02205, 0.03545, 0.05488, 0.08813  &
        ,0.14283, 0.23766, 0.42008, 0.80532, 1.48540, 2.23990  &
        ,2.65320, 2.58300, 2.41920, 2.31620, 2.23840/

        DATA(qext(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00867, 0.01418, 0.02296, 0.03696, 0.05731, 0.09213  &
        ,0.14958, 0.24986, 0.44492, 0.85650, 1.55930, 2.30050  &
        ,2.66600, 2.56720, 2.40920, 2.30910, 2.23290/

        DATA(qext(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00903, 0.01482, 0.02406, 0.03878, 0.06024, 0.09697  &
        ,0.15782, 0.26494, 0.47599, 0.91943, 1.64560, 2.36600  &
        ,2.67460, 2.54850, 2.39810, 2.30080, 2.22660/

        DATA(qext(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00947, 0.01561, 0.02541, 0.04106, 0.06390, 0.10306  &
        ,0.16825, 0.28427, 0.51639, 0.99916, 1.74810, 2.43730  &
        ,2.67710, 2.52640, 2.38520, 2.29120, 2.21930/

        DATA(qext(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.01001, 0.01661, 0.02716, 0.04401, 0.06870, 0.11107  &
        ,0.18212, 0.31048, 0.57190, 1.10440, 1.87310, 2.51460  &
        ,2.67010, 2.49990, 2.36980, 2.27940, 2.21030/

        DATA(qext(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.01074, 0.01797, 0.02957, 0.04812, 0.07542, 0.12242  &
        ,0.20207, 0.34912, 0.65467, 1.25180, 2.03180, 2.59460  &
        ,2.64740, 2.46830, 2.35040, 2.26450, 2.19900/

        DATA(qext(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.01178, 0.01999, 0.03324, 0.05452, 0.08603, 0.14061  &
        ,0.23485, 0.41494, 0.79538, 1.47760, 2.24250, 2.66530  &
        ,2.59810, 2.42900, 2.32400, 2.24430, 2.18350/

        DATA(qext(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.01349, 0.02361, 0.04019, 0.06707, 0.10742, 0.17855  &
        ,0.30676, 0.56758, 1.10240, 1.87870, 2.52380, 2.67720  &
        ,2.50420, 2.37300, 2.28190, 2.21220, 2.15910/

        DATA(qext(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.01209, 0.01960, 0.03151, 0.05042, 0.07698, 0.12332  &
        ,0.19914, 0.32792, 0.55770, 0.95268, 1.49610, 2.03990  &
        ,2.35980, 2.37690, 2.29660, 2.23220, 2.17940/

        DATA(qext(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.01237, 0.02001, 0.03214, 0.05143, 0.07859, 0.12591  &
        ,0.20340, 0.33514, 0.57015, 0.97128, 1.51680, 2.05430  &
        ,2.36200, 2.37290, 2.29380, 2.23010, 2.17790/

        DATA(qext(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.01263, 0.02042, 0.03279, 0.05248, 0.08026, 0.12863  &
        ,0.20787, 0.34272, 0.58320, 0.99062, 1.53800, 2.06890  &
        ,2.36410, 2.36880, 2.29090, 2.22800, 2.17620/

        DATA(qext(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.01288, 0.02084, 0.03348, 0.05359, 0.08203, 0.13149  &
        ,0.21258, 0.35074, 0.59698, 1.01090, 1.56000, 2.08370  &
        ,2.36610, 2.36470, 2.28800, 2.22580, 2.17460/

        DATA(qext(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.01315, 0.02128, 0.03420, 0.05477, 0.08390, 0.13453  &
        ,0.21760, 0.35929, 0.61167, 1.03220, 1.58290, 2.09880  &
        ,2.36780, 2.36050, 2.28510, 2.22360, 2.17280/

        DATA(qext(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.01344, 0.02175, 0.03498, 0.05602, 0.08591, 0.13779  &
        ,0.22299, 0.36850, 0.62746, 1.05490, 1.60700, 2.11440  &
        ,2.36940, 2.35620, 2.28200, 2.22120, 2.17100/

        DATA(qext(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.01374, 0.02226, 0.03581, 0.05738, 0.08807, 0.14132  &
        ,0.22884, 0.37851, 0.64458, 1.07920, 1.63230, 2.13040  &
        ,2.37070, 2.35170, 2.27880, 2.21880, 2.16910/

        DATA(qext(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.01407, 0.02281, 0.03671, 0.05886, 0.09044, 0.14518  &
        ,0.23525, 0.38950, 0.66334, 1.10560, 1.65930, 2.14710  &
        ,2.37180, 2.34700, 2.27550, 2.21610, 2.16710/

        DATA(qext(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.01442, 0.02341, 0.03771, 0.06049, 0.09305, 0.14944  &
        ,0.24234, 0.40169, 0.68408, 1.13430, 1.68820, 2.16430  &
        ,2.37260, 2.34210, 2.27190, 2.21330, 2.16490/

        DATA(qext(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.01481, 0.02407, 0.03881, 0.06230, 0.09596, 0.15420  &
        ,0.25028, 0.41536, 0.70723, 1.16580, 1.71940, 2.18230  &
        ,2.37300, 2.33690, 2.26810, 2.21030, 2.16250/

        DATA(qext(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.01524, 0.02481, 0.04005, 0.06434, 0.09923, 0.15957  &
        ,0.25926, 0.43087, 0.73334, 1.20080, 1.75320, 2.20110  &
        ,2.37280, 2.33140, 2.26400, 2.20710, 2.15990/

        DATA(qext(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.01572, 0.02565, 0.04146, 0.06667, 0.10297, 0.16570  &
        ,0.26955, 0.44870, 0.76313, 1.24000, 1.79020, 2.22070  &
        ,2.37200, 2.32560, 2.25960, 2.20350, 2.15710/

        DATA(qext(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.01628, 0.02661, 0.04309, 0.06935, 0.10729, 0.17282  &
        ,0.28153, 0.46949, 0.79755, 1.28440, 1.83080, 2.24120  &
        ,2.37030, 2.31920, 2.25470, 2.19950, 2.15400/

        DATA(qext(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.01692, 0.02774, 0.04499, 0.07251, 0.11238, 0.18123  &
        ,0.29574, 0.49422, 0.83796, 1.33520, 1.87590, 2.26250  &
        ,2.36750, 2.31230, 2.24920, 2.19510, 2.15050/

        DATA(qext(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.01767, 0.02907, 0.04727, 0.07631, 0.11851, 0.19140  &
        ,0.31302, 0.52435, 0.88637, 1.39450, 1.92640, 2.28440  &
        ,2.36330, 2.30460, 2.24300, 2.19000, 2.14640/

        DATA(qext(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.01858, 0.03071, 0.05008, 0.08102, 0.12614, 0.20413  &
        ,0.33475, 0.56232, 0.94596, 1.46520, 1.98390, 2.30640  &
        ,2.35720, 2.29590, 2.23580, 2.18420, 2.14180/

        DATA(qext(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.01972, 0.03279, 0.05369, 0.08712, 0.13609, 0.22081  &
        ,0.36348, 0.61250, 1.02220, 1.55190, 2.05010, 2.32730  &
        ,2.34830, 2.28570, 2.22720, 2.17710, 2.13610/

        DATA(qext(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.02122, 0.03558, 0.05864, 0.09558, 0.14998, 0.24434  &
        ,0.40440, 0.68366, 1.12520, 1.66280, 2.12760, 2.34560  &
        ,2.33560, 2.27300, 2.21630, 2.16820, 2.12900/

        DATA(qext(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.02334, 0.03972, 0.06616, 0.10869, 0.17179, 0.28185  &
        ,0.47054, 0.79670, 1.27790, 1.81370, 2.21820, 2.35610  &
        ,2.31650, 2.25560, 2.20140, 2.15600, 2.11930/

        DATA(qext(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.02685, 0.04709, 0.08028, 0.13427, 0.21559, 0.35945  &
        ,0.60968, 1.02060, 1.54900, 2.04310, 2.31850, 2.34380  &
        ,2.28470, 2.22730, 2.17750, 2.13660, 2.10400/

        DATA(qext(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00213, 0.00345, 0.00556, 0.00893, 0.01376, 0.02276  &
        ,0.04070, 0.08782, 0.23711, 0.63424, 1.54460, 2.85800  &
        ,3.34830, 2.16410, 2.41000, 2.23520, 2.15250/

        DATA(qext(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00218, 0.00353, 0.00567, 0.00911, 0.01405, 0.02327  &
        ,0.04172, 0.09039, 0.24392, 0.65093, 1.57260, 2.88340  &
        ,3.32840, 2.15330, 2.39860, 2.22840, 2.15190/

        DATA(qext(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00222, 0.00360, 0.00579, 0.00930, 0.01436, 0.02381  &
        ,0.04279, 0.09314, 0.25116, 0.66887, 1.60200, 2.90880  &
        ,3.30480, 2.14640, 2.38540, 2.22100, 2.15110/

        DATA(qext(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00227, 0.00367, 0.00591, 0.00950, 0.01468, 0.02438  &
        ,0.04394, 0.09611, 0.25891, 0.68830, 1.63310, 2.93380  &
        ,3.27840, 2.13970, 2.37090, 2.21320, 2.14980/

        DATA(qext(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00232, 0.00375, 0.00604, 0.00971, 0.01503, 0.02498  &
        ,0.04518, 0.09934, 0.26730, 0.70949, 1.66610, 2.95830  &
        ,3.25100, 2.13040, 2.35460, 2.20510, 2.14820/

        DATA(qext(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00237, 0.00384, 0.00618, 0.00993, 0.01540, 0.02563  &
        ,0.04652, 0.10289, 0.27645, 0.73278, 1.70100, 2.98230  &
        ,3.22320, 2.12250, 2.33660, 2.19700, 2.14620/

        DATA(qext(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00242, 0.00393, 0.00633, 0.01018, 0.01580, 0.02634  &
        ,0.04800, 0.10684, 0.28653, 0.75853, 1.73830, 3.00700  &
        ,3.19320, 2.12200, 2.31720, 2.18910, 2.14370/

        DATA(qext(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00248, 0.00402, 0.00649, 0.01044, 0.01624, 0.02712  &
        ,0.04965, 0.11129, 0.29776, 0.78720, 1.77840, 3.03400  &
        ,3.15690, 2.12510, 2.29540, 2.18190, 2.14100/

        DATA(qext(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00254, 0.00413, 0.00666, 0.01074, 0.01672, 0.02799  &
        ,0.05151, 0.11635, 0.31041, 0.81926, 1.82210, 3.06570  &
        ,3.11110, 2.12550, 2.27320, 2.17590, 2.13810/

        DATA(qext(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00261, 0.00425, 0.00686, 0.01106, 0.01726, 0.02897  &
        ,0.05363, 0.12220, 0.32482, 0.85527, 1.87090, 3.10360  &
        ,3.05750, 2.13550, 2.24940, 2.17170, 2.13520/

        DATA(qext(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00269, 0.00438, 0.00708, 0.01143, 0.01787, 0.03008  &
        ,0.05609, 0.12905, 0.34147, 0.89586, 1.92680, 3.14610  &
        ,2.99950, 2.15290, 2.22670, 2.16980, 2.13260/

        DATA(qext(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00277, 0.00453, 0.00733, 0.01185, 0.01857, 0.03137  &
        ,0.05899, 0.13720, 0.36101, 0.94179, 1.99290, 3.18670  &
        ,2.93040, 2.17320, 2.20670, 2.17040, 2.13020/

        DATA(qext(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00287, 0.00470, 0.00762, 0.01234, 0.01939, 0.03289  &
        ,0.06247, 0.14709, 0.38439, 0.99421, 2.07290, 3.21980  &
        ,2.84340, 2.20760, 2.19170, 2.17280, 2.12780/

        DATA(qext(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00299, 0.00490, 0.00796, 0.01291, 0.02035, 0.03472  &
        ,0.06675, 0.15938, 0.41314, 1.05520, 2.16900, 3.25110  &
        ,2.74560, 2.24990, 2.18570, 2.17450, 2.12470/

        DATA(qext(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00312, 0.00514, 0.00837, 0.01360, 0.02152, 0.03697  &
        ,0.07220, 0.17510, 0.44983, 1.12920, 2.28080, 3.29000  &
        ,2.62430, 2.29930, 2.19220, 2.17170, 2.12060/

        DATA(qext(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00328, 0.00543, 0.00887, 0.01446, 0.02299, 0.03987  &
        ,0.07946, 0.19601, 0.49931, 1.22580, 2.40990, 3.30810  &
        ,2.48670, 2.35200, 2.21070, 2.16180, 2.11610/

        DATA(qext(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00348, 0.00580, 0.00952, 0.01558, 0.02494, 0.04383  &
        ,0.08976, 0.22538, 0.57139, 1.36440, 2.57810, 3.29950  &
        ,2.33350, 2.38660, 2.22800, 2.14920, 2.11090/

        DATA(qext(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00375, 0.00629, 0.01041, 0.01714, 0.02771, 0.04973  &
        ,0.10589, 0.27002, 0.68707, 1.57060, 2.79420, 3.22210  &
        ,2.18910, 2.36130, 2.21230, 2.14280, 2.10420/

        DATA(qext(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00413, 0.00703, 0.01176, 0.01959, 0.03222, 0.06002  &
        ,0.13567, 0.34740, 0.88707, 1.87840, 3.06290, 2.99080  &
        ,2.15230, 2.23840, 2.16860, 2.12980, 2.09540/

        DATA(qext(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00475, 0.00834, 0.01432, 0.02449, 0.04203, 0.08539  &
        ,0.21228, 0.53671, 1.29000, 2.47370, 3.27290, 2.40570  &
        ,2.36280, 2.21670, 2.15200, 2.11140, 2.08170/

        DATA(qext(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00239, 0.00388, 0.00624, 0.01001, 0.01536, 0.02502  &
        ,0.04264, 0.08196, 0.19158, 0.47030, 1.13170, 2.22660  &
        ,3.17280, 2.62430, 2.28710, 2.21190, 2.16120/

        DATA(qext(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00245, 0.00396, 0.00637, 0.01021, 0.01568, 0.02557  &
        ,0.04364, 0.08411, 0.19678, 0.48106, 1.15000, 2.24760  &
        ,3.17130, 2.59950, 2.29250, 2.21170, 2.15940/

        DATA(qext(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00250, 0.00404, 0.00650, 0.01043, 0.01602, 0.02614  &
        ,0.04468, 0.08640, 0.20231, 0.49261, 1.16970, 2.26990  &
        ,3.17000, 2.57430, 2.29710, 2.21150, 2.15750/

        DATA(qext(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00255, 0.00413, 0.00664, 0.01065, 0.01638, 0.02674  &
        ,0.04580, 0.08885, 0.20822, 0.50511, 1.19100, 2.29390  &
        ,3.16900, 2.54940, 2.30220, 2.21120, 2.15560/

        DATA(qext(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00261, 0.00422, 0.00678, 0.01089, 0.01676, 0.02739  &
        ,0.04699, 0.09149, 0.21460, 0.51877, 1.21430, 2.32010  &
        ,3.16780, 2.52400, 2.30690, 2.21080, 2.15380/

        DATA(qext(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00266, 0.00431, 0.00694, 0.01114, 0.01717, 0.02808  &
        ,0.04828, 0.09438, 0.22156, 0.53386, 1.24010, 2.34870  &
        ,3.16550, 2.49650, 2.31050, 2.21010, 2.15190/

        DATA(qext(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00272, 0.00441, 0.00711, 0.01141, 0.01762, 0.02883  &
        ,0.04969, 0.09757, 0.22922, 0.55071, 1.26860, 2.38000  &
        ,3.16080, 2.46670, 2.31410, 2.20910, 2.15010/

        DATA(qext(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00279, 0.00452, 0.00729, 0.01171, 0.01810, 0.02966  &
        ,0.05125, 0.10113, 0.23774, 0.56974, 1.30050, 2.41400  &
        ,3.15320, 2.43670, 2.31650, 2.20770, 2.14820/

        DATA(qext(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00286, 0.00464, 0.00749, 0.01204, 0.01863, 0.03057  &
        ,0.05299, 0.10516, 0.24732, 0.59146, 1.33640, 2.45080  &
        ,3.14320, 2.40660, 2.31750, 2.20550, 2.14630/

        DATA(qext(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00294, 0.00478, 0.00771, 0.01240, 0.01923, 0.03160  &
        ,0.05496, 0.10977, 0.25821, 0.61649, 1.37670, 2.49050  &
        ,3.13200, 2.37340, 2.31690, 2.20250, 2.14420/

        DATA(qext(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00302, 0.00492, 0.00796, 0.01281, 0.01990, 0.03276  &
        ,0.05723, 0.11513, 0.27075, 0.64563, 1.42240, 2.53400  &
        ,3.11790, 2.34020, 2.31380, 2.19850, 2.14190/

        DATA(qext(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00312, 0.00509, 0.00824, 0.01328, 0.02067, 0.03410  &
        ,0.05986, 0.12145, 0.28538, 0.67986, 1.47450, 2.58300  &
        ,3.09590, 2.30790, 2.30720, 2.19340, 2.13910/

        DATA(qext(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00323, 0.00529, 0.00857, 0.01383, 0.02156, 0.03566  &
        ,0.06299, 0.12904, 0.30273, 0.72038, 1.53480, 2.63990  &
        ,3.06460, 2.27550, 2.29670, 2.18730, 2.13570/

        DATA(qext(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00336, 0.00551, 0.00895, 0.01446, 0.02261, 0.03752  &
        ,0.06678, 0.13839, 0.32377, 0.76877, 1.60660, 2.70490  &
        ,3.02570, 2.24770, 2.28130, 2.18060, 2.13180/

        DATA(qext(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00351, 0.00578, 0.00940, 0.01523, 0.02388, 0.03980  &
        ,0.07152, 0.15024, 0.35002, 0.82727, 1.69520, 2.77490  &
        ,2.96670, 2.22730, 2.26110, 2.17410, 2.12750/

        DATA(qext(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00369, 0.00610, 0.00997, 0.01619, 0.02547, 0.04269  &
        ,0.07770, 0.16587, 0.38424, 0.89983, 1.80760, 2.85230  &
        ,2.88670, 2.21860, 2.23760, 2.16810, 2.12280/

        DATA(qext(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00392, 0.00652, 0.01069, 0.01743, 0.02755, 0.04656  &
        ,0.08625, 0.18771, 0.43198, 0.99545, 1.95140, 2.93730  &
        ,2.76920, 2.22880, 2.21580, 2.16160, 2.11720/

        DATA(qext(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00422, 0.00708, 0.01168, 0.01915, 0.03049, 0.05219  &
        ,0.09923, 0.22079, 0.50632, 1.13760, 2.14720, 3.01460  &
        ,2.60090, 2.26170, 2.20290, 2.15090, 2.11010/

        DATA(qext(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00464, 0.00791, 0.01320, 0.02183, 0.03520, 0.06163  &
        ,0.12230, 0.27801, 0.64258, 1.37710, 2.43090, 3.03960  &
        ,2.37750, 2.29060, 2.19190, 2.13760, 2.10070/

        DATA(qext(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00534, 0.00938, 0.01605, 0.02713, 0.04504, 0.08336  &
        ,0.17931, 0.41017, 0.94219, 1.85680, 2.85380, 2.80770  &
        ,2.22470, 2.22210, 2.16160, 2.11770, 2.08630/

        DATA(qext(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00177, 0.00288, 0.00466, 0.00761, 0.01236, 0.02417  &
        ,0.06347, 0.21225, 0.66898, 1.71760, 3.12360, 3.25790  &
        ,2.27840, 2.29480, 2.24140, 2.15630, 2.11720/

        DATA(qext(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00182, 0.00294, 0.00476, 0.00777, 0.01264, 0.02484  &
        ,0.06566, 0.21914, 0.68848, 1.74790, 3.14920, 3.22830  &
        ,2.28360, 2.29320, 2.23880, 2.15360, 2.11640/

        DATA(qext(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00185, 0.00300, 0.00485, 0.00793, 0.01295, 0.02556  &
        ,0.06803, 0.22651, 0.70921, 1.77990, 3.17580, 3.19910  &
        ,2.28710, 2.29300, 2.23500, 2.15120, 2.11530/

        DATA(qext(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00189, 0.00307, 0.00496, 0.00811, 0.01327, 0.02634  &
        ,0.07060, 0.23445, 0.73136, 1.81430, 3.20270, 3.16830  &
        ,2.28770, 2.29420, 2.23000, 2.14930, 2.11380/

        DATA(qext(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00193, 0.00313, 0.00507, 0.00830, 0.01361, 0.02718  &
        ,0.07342, 0.24308, 0.75519, 1.85160, 3.22930, 3.13290  &
        ,2.29180, 2.29520, 2.22300, 2.14840, 2.11200/

        DATA(qext(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00197, 0.00320, 0.00518, 0.00850, 0.01399, 0.02810  &
        ,0.07655, 0.25256, 0.78103, 1.89270, 3.25500, 3.09140  &
        ,2.30210, 2.29600, 2.21520, 2.14840, 2.11020/

        DATA(qext(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00202, 0.00328, 0.00531, 0.00871, 0.01440, 0.02913  &
        ,0.08006, 0.26308, 0.80925, 1.93820, 3.27980, 3.04470  &
        ,2.31160, 2.29960, 2.20630, 2.14920, 2.10870/

        DATA(qext(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00207, 0.00336, 0.00544, 0.00895, 0.01485, 0.03029  &
        ,0.08405, 0.27490, 0.84036, 1.98920, 3.30450, 2.99600  &
        ,2.31780, 2.29990, 2.19720, 2.15010, 2.10740/

        DATA(qext(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00212, 0.00345, 0.00560, 0.00921, 0.01535, 0.03161  &
        ,0.08865, 0.28834, 0.87498, 2.04630, 3.33100, 2.94370  &
        ,2.32970, 2.30200, 2.18880, 2.15000, 2.10600/

        DATA(qext(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00218, 0.00354, 0.00576, 0.00951, 0.01592, 0.03313  &
        ,0.09402, 0.30384, 0.91397, 2.10990, 3.36090, 2.88030  &
        ,2.34080, 2.30070, 2.18270, 2.14770, 2.10440/

        DATA(qext(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00224, 0.00365, 0.00595, 0.00984, 0.01658, 0.03493  &
        ,0.10038, 0.32200, 0.95855, 2.18040, 3.39140, 2.80790  &
        ,2.35260, 2.29680, 2.18000, 2.14270, 2.10270/

        DATA(qext(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00231, 0.00378, 0.00617, 0.01023, 0.01735, 0.03707  &
        ,0.10805, 0.34366, 1.01060, 2.25850, 3.41560, 2.73660  &
        ,2.35920, 2.28670, 2.18150, 2.13670, 2.10110/

        DATA(qext(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00239, 0.00392, 0.00641, 0.01068, 0.01826, 0.03970  &
        ,0.11747, 0.37015, 1.07290, 2.34660, 3.42990, 2.64860  &
        ,2.36740, 2.27010, 2.18470, 2.13250, 2.09850/

        DATA(qext(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00249, 0.00409, 0.00671, 0.01121, 0.01937, 0.04300  &
        ,0.12934, 0.40351, 1.14980, 2.45100, 3.44240, 2.55890  &
        ,2.36480, 2.24720, 2.18500, 2.13050, 2.09550/

        DATA(qext(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00260, 0.00429, 0.00706, 0.01186, 0.02077, 0.04730  &
        ,0.14474, 0.44726, 1.24790, 2.57970, 3.44160, 2.46190  &
        ,2.35110, 2.22250, 2.17720, 2.12690, 2.09300/

        DATA(qext(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00273, 0.00453, 0.00749, 0.01269, 0.02260, 0.05319  &
        ,0.16556, 0.50753, 1.37520, 2.73180, 3.40190, 2.37440  &
        ,2.32520, 2.20960, 2.16710, 2.12280, 2.08970/

        DATA(qext(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00290, 0.00484, 0.00806, 0.01379, 0.02516, 0.06182  &
        ,0.19530, 0.59537, 1.54120, 2.90450, 3.32350, 2.30460  &
        ,2.29590, 2.22100, 2.16390, 2.11660, 2.08560/

        DATA(qext(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00312, 0.00526, 0.00884, 0.01538, 0.02911, 0.07586  &
        ,0.24142, 0.73009, 1.77310, 3.12150, 3.14720, 2.29110  &
        ,2.28460, 2.22470, 2.14560, 2.11030, 2.08040/

        DATA(qext(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00344, 0.00589, 0.01006, 0.01804, 0.03634, 0.10299  &
        ,0.32442, 0.95126, 2.14830, 3.34280, 2.80740, 2.34140  &
        ,2.29160, 2.17850, 2.13990, 2.10020, 2.07370/

        DATA(qext(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00396, 0.00701, 0.01245, 0.02405, 0.05555, 0.17662  &
        ,0.54577, 1.45180, 2.80730, 3.34810, 2.33220, 2.30030  &
        ,2.21140, 2.16340, 2.11850, 2.08590, 2.06320/

! Qext for Mineral Dust R(imaginary) = 0.0015 (LOW absorption)(No RH effect)
        DATA(qext_dust(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00053,  0.00086,  0.00144,  0.00276,  0.00728,  0.02876  &
        , 0.12569,  0.44984,  1.21290,  2.29250,  2.64450,  2.66060  &
        , 2.54140,  2.19160,  2.20710,  2.13410,  2.09580/ 
 
        DATA(qext_dust(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00050,  0.00128,  0.00509,  0.02746,  0.15370,  0.67428  &
        , 2.01810,  3.46660,  3.06260,  2.65140,  2.32890,  2.28180  &
        , 2.19870,  2.14360,  2.10470,  2.07500,  2.05570/ 
 
        DATA(qext_dust(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00246,  0.01192,  0.06727,  0.33158,  1.24500,  2.93860  &
        , 3.77620,  2.62540,  2.43390,  2.24850,  2.21360,  2.15880  &
        , 2.11950,  2.09170,  2.06560,  2.04800,  2.03390/ 
 
        DATA(qext_dust(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00232,  0.00368,  0.00584,  0.00927,  0.01472,  0.02341  &
        , 0.03732,  0.05997,  0.09838,  0.17106,  0.34448,  0.88580  &
        , 2.31920,  3.26560,  2.89620,  2.60710,  2.45250/ 
 
        DATA(qext_dust(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00154,  0.00244,  0.00387,  0.00615,  0.00978,  0.01558  &
        , 0.02498,  0.04082,  0.07056,  0.14147,  0.36788,  1.06280  &
        , 2.46800,  3.53330,  2.71720,  2.50290,  2.33060/ 
 
        DATA(qext_dust(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00319,  0.00506,  0.00804,  0.01277,  0.02033,  0.03253  &
        , 0.05279,  0.08919,  0.16752,  0.38114,  0.93498,  1.87920  &
        , 2.68110,  2.74410,  2.47570,  2.30710,  2.22650/ 
 
        DATA(qext_dust(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00347,  0.00550,  0.00874,  0.01389,  0.02215,  0.03557  &
        , 0.05837,  0.10191,  0.20923,  0.57176,  1.73950,  3.22860  &
        , 3.29090,  2.47190,  2.45570,  2.32820,  2.24460/ 
 
        DATA(qext_dust(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00147,  0.00234,  0.00371,  0.00591,  0.00942,  0.01519  &
        , 0.02535,  0.04659,  0.10365,  0.25753,  0.65195,  1.44010  &
        , 2.51080,  2.84480,  2.33600,  2.25070,  2.18000/ 

! Qext for Mineral Dust R(imaginary) = 0.003 (MID absorption)(No RH effect)
        DATA(qext_dust(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00070,  0.00112,  0.00185,  0.00343,  0.00839,  0.03067  &
        , 0.12895,  0.45450,  1.21520,  2.28620,  2.63730,  2.65780  &
        , 2.54130,  2.19440,  2.20660,  2.13270,  2.09540/ 
 
        DATA(qext_dust(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00090,  0.00193,  0.00616,  0.02932,  0.15704,  0.67954  &
        , 2.01910,  3.45710,  3.05920,  2.64770,  2.33040,  2.28090  &
        , 2.19820,  2.14350,  2.10420,  2.07570,  2.05580/ 
 
        DATA(qext_dust(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00327,  0.01328,  0.06968,  0.33582,  1.24930,  2.93440  &
        , 3.76520,  2.62660,  2.43290,  2.25180,  2.21400,  2.15940  &
        , 2.12010,  2.09120,  2.06570,  2.04800,  2.03400/ 
 
        DATA(qext_dust(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00232,  0.00368,  0.00584,  0.00927,  0.01472,  0.02341  &
        , 0.03732,  0.05997,  0.09838,  0.17106,  0.34448,  0.88580  &
        , 2.31920,  3.26560,  2.89620,  2.60710,  2.45250/ 
 
        DATA(qext_dust(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00154,  0.00244,  0.00387,  0.00615,  0.00978,  0.01558  &
        , 0.02498,  0.04082,  0.07056,  0.14147,  0.36788,  1.06280  &
        , 2.46800,  3.53330,  2.71720,  2.50290,  2.33060/ 
 
        DATA(qext_dust(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00319,  0.00506,  0.00804,  0.01277,  0.02033,  0.03253  &
        , 0.05279,  0.08919,  0.16752,  0.38114,  0.93498,  1.87920  &
        , 2.68110,  2.74410,  2.47570,  2.30710,  2.22650/ 
 
        DATA(qext_dust(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00347,  0.00550,  0.00874,  0.01389,  0.02215,  0.03557  &
        , 0.05837,  0.10191,  0.20923,  0.57176,  1.73950,  3.22860  &
        , 3.29090,  2.47190,  2.45570,  2.32820,  2.24460/ 
 
        DATA(qext_dust(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00147,  0.00234,  0.00371,  0.00591,  0.00942,  0.01519  &
        , 0.02535,  0.04659,  0.10365,  0.25753,  0.65195,  1.44010  &
        , 2.51080,  2.84480,  2.33600,  2.25070,  2.18000/ 

! Qext for Mineral Dust R(imaginary) = 0.008 (HIGH absorption)(No RH effect)
        DATA(qext_dust(3,1,ib),ib=1,17)  & ! Band 1, All RH's
        /0.00124, 0.00198, 0.00323, 0.00566, 0.01209, 0.03703  &
        ,0.13977, 0.46990, 1.22290, 2.26570, 2.61390, 2.64980  &
        ,2.54470, 2.19980, 2.20590, 2.13130, 2.09540/

        DATA(qext_dust(3,2,ib),ib=1,17)  & ! Band 2, All RH's
        /0.00224, 0.00409, 0.00974, 0.03551, 0.16816, 0.69694  &
        ,2.02260, 3.42610, 3.04780, 2.63640, 2.33260, 2.27750  &
        ,2.19670, 2.14290, 2.10390, 2.07640, 2.05630/

        DATA(qext_dust(3,3,ib),ib=1,17)  & ! Band 3, All RH's
        /0.00599, 0.01846, 0.08179, 0.36834, 1.28590, 2.91150  &
        ,3.72340, 2.62500, 2.42800, 2.26110, 2.21600, 2.16130  &
        ,2.12100, 2.09050, 2.06620, 2.04880, 2.03590/

        DATA(qext_dust(3,4,ib),ib=1,17)  & ! Band 4, All RH's
        /0.00232, 0.00368, 0.00584, 0.00927, 0.01472, 0.02341  &
        ,0.03733, 0.05997, 0.09838, 0.17106, 0.34448, 0.88581  &
        ,2.31920, 3.26570, 2.89620, 2.60710, 2.45250/

        DATA(qext_dust(3,5,ib),ib=1,17)  & ! Band 5, All RH's
        /0.00154, 0.00240, 0.00387, 0.00615, 0.00978, 0.01558  &
        ,0.02498, 0.04083, 0.07056, 0.14147, 0.36788, 1.06280  &
        ,2.46800, 3.53330, 2.71720, 2.50290, 2.33060/

        DATA(qext_dust(3,6,ib),ib=1,17)  & ! Band 6, All RH's
        /0.00319, 0.00506, 0.00804, 0.01277, 0.02033, 0.03253  &
        ,0.05279, 0.08919, 0.16752, 0.38114, 0.93498, 1.87920  &
        ,2.68110, 2.74410, 2.47570, 2.30710, 2.22650/

        DATA(qext_dust(3,7,ib),ib=1,17)  & ! Band 7, All RH's
        /0.00347, 0.00550, 0.00874, 0.01389, 0.02215, 0.03557  &
        ,0.05837, 0.10191, 0.20923, 0.57176, 1.73950, 3.22860  &
        ,3.29090, 2.47200, 2.45570, 2.32820, 2.24460/

        DATA(qext_dust(3,8,ib),ib=1,17)  & ! Band 8, All RH's
        /0.00147, 0.00234, 0.00371, 0.00591, 0.00942, 0.01519  &
        ,0.02535, 0.04659, 0.10365, 0.25753, 0.65195, 1.44010  &
        ,2.51080, 2.84480, 2.33600, 2.25070, 2.18000/

! Qext for Absorbing Carbon (1% BC, 99% OC)(No RH effect)
        DATA(qext_carb(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00245,  0.00392,  0.00643,  0.01137,  0.02481,  0.07891  &
        , 0.31208,  1.10800,  2.52240,  3.47310,  2.68080,  2.52040  &
        , 2.29430,  2.23790,  2.17310,  2.12700,  2.09340/

        DATA(qext_carb(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00547,  0.00931,  0.01846,  0.05122,  0.20027,  0.76872  &
        , 2.14280,  3.52800,  3.06150,  2.48770,  2.40170,  2.26170  &
        , 2.19620,  2.14160,  2.10360,  2.07600,  2.05560/

        DATA(qext_carb(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00903,  0.02314,  0.08790,  0.37122,  1.29850,  2.93160  &
        , 3.68880,  2.61620,  2.43400,  2.26980,  2.21970,  2.16330  &
        , 2.12150,  2.08990,  2.06580,  2.04790,  2.03400/

        DATA(qext_carb(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00207,  0.00329,  0.00522,  0.00828,  0.01316,  0.02092  &
        , 0.03338,  0.05370,  0.08851,  0.15621,  0.32870,  0.92364  &
        , 2.49150,  3.30090,  2.84960,  2.61440,  2.45410/

        DATA(qext_carb(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00213,  0.00338,  0.00536,  0.00852,  0.01354,  0.02158  &
        , 0.03461,  0.05658,  0.09779,  0.19599,  0.51048,  1.50240  &
        , 2.95000,  3.32070,  2.51700,  2.45020,  2.32250/

        DATA(qext_carb(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00454,  0.00721,  0.01145,  0.01824,  0.02920,  0.04748  &
        , 0.08070,  0.15532,  0.40043,  1.55040,  3.49400,  3.14430  &
        , 2.69580,  2.55150,  2.40690,  2.30370,  2.22700/

        DATA(qext_carb(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00286,  0.00454,  0.00721,  0.01147,  0.01833,  0.02967  &
        , 0.04984,  0.09332,  0.22860,  0.84677,  2.86630,  3.66030  &
        , 2.88190,  2.56760,  2.41750,  2.32800,  2.24340/

        DATA(qext_carb(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00161,  0.00256,  0.00407,  0.00648,  0.01042,  0.01720  &
        , 0.03094,  0.06965,  0.21644,  0.70260,  1.89450,  3.28730  &
        , 3.07130,  2.40340,  2.32590,  2.25300,  2.17850/

! Qext for Absorbing Carbon (2% BC, 98% OC)(No RH effect)
        DATA(qext_carb(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00342,  0.00547,  0.00892,  0.01543,  0.03173,  0.09177  &
        , 0.33775,  1.14950,  2.54460,  3.42600,  2.67820,  2.50580  &
        , 2.29910,  2.23690,  2.17270,  2.12680,  2.09320/

        DATA(qext_carb(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00750,  0.01260,  0.02395,  0.06106,  0.21973,  0.80684  &
        , 2.17380,  3.49700,  3.02970,  2.48110,  2.39390,  2.26130  &
        , 2.19440,  2.14100,  2.10330,  2.07590,  2.05550/

        DATA(qext_carb(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.01321,  0.03033,  0.10142,  0.39842,  1.33920,  2.93800  &
        , 3.63570,  2.60750,  2.43470,  2.27660,  2.22210,  2.16420  &
        , 2.12160,  2.08970,  2.06580,  2.04790,  2.03400/

        DATA(qext_carb(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00211,  0.00334,  0.00531,  0.00842,  0.01338,  0.02128  &
        , 0.03395,  0.05463,  0.09002,  0.15880,  0.33355,  0.93306  &
        , 2.49280,  3.29170,  2.84880,  2.61350,  2.45390/

        DATA(qext_carb(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00222,  0.00353,  0.00561,  0.00890,  0.01415,  0.02255  &
        , 0.03617,  0.05912,  0.10209,  0.20388,  0.52654,  1.53030  &
        , 2.95560,  3.29540,  2.52100,  2.44760,  2.32290/

        DATA(qext_carb(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00463,  0.00736,  0.01170,  0.01863,  0.02982,  0.04848  &
        , 0.08234,  0.15815,  0.40567,  1.55670,  3.48360,  3.14260  &
        , 2.69610,  2.55100,  2.40700,  2.30380,  2.22710/

        DATA(qext_carb(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00297,  0.00471,  0.00749,  0.01192,  0.01905,  0.03081  &
        , 0.05172,  0.09650,  0.23444,  0.85764,  2.86680,  3.64720  &
        , 2.87810,  2.56980,  2.42000,  2.32800,  2.24340/

        DATA(qext_carb(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00196,  0.00312,  0.00496,  0.00790,  0.01268,  0.02086  &
        , 0.03705,  0.08057,  0.23796,  0.74521,  1.94130,  3.27520  &
        , 3.01750,  2.41950,  2.32940,  2.24830,  2.17860/

! Returning value for qext:
        if (aerotype <= 2) value = qext(aerotype,radband,bin,rh)
        if (aerotype == 3) value = qext_dust(dust_ref_im,radband,bin)
        if (aerotype == 4) value = qext_carb(1,radband,bin)
        if (aerotype == 5) value = qext_carb(2,radband,bin)
return
END SUBROUTINE aeroqext

!##############################################################################
Subroutine aeroqscat (aerotype,radband,bin,rh,value)

use rrad3, only:dust_ref_im

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL qscat(2,8,17,20),value
        REAL qscat_dust(3,8,17)
        REAL qscat_carb(2,8,17)

! Qscat Ammonium sulfate with 20% insoluble inclusion (dust low R(Im)))
        DATA(qscat(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        / 0.00001,  0.00006,  0.00041,  0.00262,  0.01488,  0.08427  &
        , 0.34594,  1.12330,  2.47130,  3.25330,  2.28770,  2.11400  &
        , 1.97460,  1.89150,  1.81820,  1.72860,  1.63710/

        DATA(qscat(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        / 0.00001,  0.00006,  0.00042,  0.00270,  0.01536,  0.08657  &
        , 0.35310,  1.13920,  2.49080,  3.24820,  2.27460,  2.09340  &
        , 1.96480,  1.88360,  1.81340,  1.72510,  1.63570/

        DATA(qscat(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        / 0.00001,  0.00007,  0.00044,  0.00280,  0.01588,  0.08910  &
        , 0.36102,  1.15670,  2.51210,  3.24080,  2.25990,  2.09770  &
        , 1.96380,  1.88850,  1.81160,  1.72280,  1.63510/

        DATA(qscat(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        / 0.00001,  0.00007,  0.00045,  0.00291,  0.01647,  0.09190  &
        , 0.36985,  1.17600,  2.53510,  3.23100,  2.25210,  2.07070  &
        , 1.97100,  1.87960,  1.81280,  1.72460,  1.63610/

        DATA(qscat(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        / 0.00001,  0.00007,  0.00047,  0.00303,  0.01713,  0.09502  &
        , 0.37976,  1.19740,  2.55990,  3.21890,  2.24540,  2.08230  &
        , 1.96320,  1.89000,  1.80980,  1.72550,  1.63600/

        DATA(qscat(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        / 0.00001,  0.00007,  0.00049,  0.00316,  0.01788,  0.09850  &
        , 0.39096,  1.22120,  2.58640,  3.20540,  2.23690,  2.05350  &
        , 1.96770,  1.88250,  1.80680,  1.72330,  1.63480/

        DATA(qscat(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        / 0.00001,  0.00008,  0.00051,  0.00332,  0.01874,  0.10245  &
        , 0.40374,  1.24780,  2.61470,  3.19230,  2.22350,  2.04680  &
        , 1.96580,  1.88060,  1.80380,  1.72360,  1.62980/

        DATA(qscat(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        / 0.00001,  0.00008,  0.00054,  0.00350,  0.01973,  0.10694  &
        , 0.41845,  1.27750,  2.64480,  3.18020,  2.20750,  2.04960  &
        , 1.97200,  1.88710,  1.79840,  1.72010,  1.62920/

        DATA(qscat(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        / 0.00001,  0.00009,  0.00057,  0.00371,  0.02089,  0.11213  &
        , 0.43560,  1.31100,  2.67760,  3.16480,  2.20340,  2.03160  &
        , 1.97410,  1.88860,  1.79410,  1.71760,  1.62710/

        DATA(qscat(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        / 0.00001,  0.00009,  0.00061,  0.00397,  0.02228,  0.11819  &
        , 0.45580,  1.34910,  2.71430,  3.14030,  2.19870,  2.01660  &
        , 1.97110,  1.88150,  1.79690,  1.71240,  1.62450/

        DATA(qscat(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        / 0.00001,  0.00010,  0.00066,  0.00427,  0.02396,  0.12540  &
        , 0.47992,  1.39290,  2.75750,  3.10390,  2.18730,  2.02420  &
        , 1.97270,  1.87050,  1.78990,  1.70920,  1.62060/

        DATA(qscat(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        / 0.00002,  0.00011,  0.00071,  0.00466,  0.02604,  0.13412  &
        , 0.50910,  1.44410,  2.81030,  3.06100,  2.18660,  2.00370  &
        , 1.97050,  1.86580,  1.79790,  1.70880,  1.61930/

        DATA(qscat(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        / 0.00002,  0.00012,  0.00078,  0.00514,  0.02869,  0.14490  &
        , 0.54490,  1.50590,  2.87280,  3.01630,  2.18330,  2.00000  &
        , 1.96910,  1.85610,  1.79000,  1.70200,  1.61310/

        DATA(qscat(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        / 0.00002,  0.00013,  0.00088,  0.00579,  0.03218,  0.15862  &
        , 0.58947,  1.58360,  2.93980,  2.94330,  2.19310,  2.00420  &
        , 1.94380,  1.85870,  1.77470,  1.69710,  1.60760/

        DATA(qscat(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        / 0.00002,  0.00014,  0.00100,  0.00668,  0.03696,  0.17667  &
        , 0.64588,  1.68580,  3.00920,  2.85090,  2.20960,  2.03300  &
        , 1.92560,  1.86600,  1.78090,  1.69270,  1.60170/

        DATA(qscat(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        / 0.00002,  0.00017,  0.00119,  0.00798,  0.04384,  0.20153  &
        , 0.71929,  1.82120,  3.09630,  2.71670,  2.22650,  2.05100  &
        , 1.92650,  1.85290,  1.76580,  1.68450,  1.59420/

        DATA(qscat(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        / 0.00003,  0.00020,  0.00147,  0.01001,  0.05440,  0.23822  &
        , 0.82090,  1.99710,  3.17410,  2.54890,  2.23630,  2.04750  &
        , 1.92190,  1.83630,  1.75680,  1.67450,  1.58250/

        DATA(qscat(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        / 0.00003,  0.00026,  0.00196,  0.01354,  0.07223,  0.29942  &
        , 0.98224,  2.25040,  3.23060,  2.33270,  2.17640,  2.00390  &
        , 1.90610,  1.81800,  1.74810,  1.65870,  1.56840/

        DATA(qscat(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        / 0.00004,  0.00037,  0.00295,  0.02099,  0.10710,  0.42365  &
        , 1.26870,  2.61800,  3.14690,  2.15930,  2.04100,  1.95550  &
        , 1.88310,  1.79710,  1.72150,  1.63510,  1.54220/

        DATA(qscat(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        / 0.00006,  0.00065,  0.00577,  0.04291,  0.19708,  0.73191  &
        , 1.87400,  3.12550,  2.61160,  2.22190,  2.06750,  1.93050  &
        , 1.84030,  1.76180,  1.68050,  1.59160,  1.50230/

        DATA(qscat(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        / 0.00022,  0.00144,  0.00926,  0.05594,  0.24703,  0.91027  &
        , 2.28510,  3.57260,  2.87320,  2.38200,  2.33280,  2.20410  &
        , 2.13570,  2.08510,  2.02930,  1.97580,  1.91940/

        DATA(qscat(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        / 0.00023,  0.00148,  0.00953,  0.05748,  0.25179,  0.92179  &
        , 2.30210,  3.57610,  2.85900,  2.37840,  2.33090,  2.20810  &
        , 2.13850,  2.08560,  2.03000,  1.98170,  1.92190/

        DATA(qscat(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        / 0.00023,  0.00153,  0.00984,  0.05918,  0.25702,  0.93442  &
        , 2.32040,  3.57880,  2.84500,  2.36900,  2.32270,  2.20510  &
        , 2.13780,  2.08380,  2.02990,  1.97930,  1.92210/

        DATA(qscat(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        / 0.00024,  0.00158,  0.01017,  0.06106,  0.26279,  0.94837  &
        , 2.34030,  3.58050,  2.82860,  2.36080,  2.32070,  2.20460  &
        , 2.13710,  2.08060,  2.02640,  1.98000,  1.92520/

        DATA(qscat(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        / 0.00025,  0.00163,  0.01055,  0.06316,  0.26920,  0.96390  &
        , 2.36210,  3.58100,  2.80980,  2.35290,  2.31340,  2.20600  &
        , 2.13330,  2.07530,  2.03180,  1.98500,  1.92680/

        DATA(qscat(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        / 0.00026,  0.00170,  0.01098,  0.06552,  0.27639,  0.98135  &
        , 2.38610,  3.58070,  2.78780,  2.34710,  2.31900,  2.19650  &
        , 2.13340,  2.07350,  2.03110,  1.98570,  1.92930/

        DATA(qscat(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        / 0.00027,  0.00177,  0.01147,  0.06818,  0.28452,  1.00120  &
        , 2.41300,  3.58020,  2.76250,  2.33890,  2.30190,  2.19370  &
        , 2.13120,  2.07800,  2.02800,  1.98290,  1.93170/

        DATA(qscat(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        / 0.00028,  0.00185,  0.01202,  0.07123,  0.29383,  1.02390  &
        , 2.44370,  3.58100,  2.73780,  2.32730,  2.30480,  2.20370  &
        , 2.12550,  2.07520,  2.02980,  1.98740,  1.93390/

        DATA(qscat(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        / 0.00029,  0.00195,  0.01268,  0.07476,  0.30463,  1.05040  &
        , 2.47940,  3.58390,  2.71570,  2.31140,  2.30070,  2.19670  &
        , 2.12420,  2.07450,  2.03270,  1.98930,  1.93560/

        DATA(qscat(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        / 0.00030,  0.00206,  0.01345,  0.07889,  0.31735,  1.08170  &
        , 2.52150,  3.58710,  2.68800,  2.29600,  2.29730,  2.19480  &
        , 2.12230,  2.07430,  2.03320,  1.98900,  1.93880/

        DATA(qscat(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        / 0.00032,  0.00220,  0.01437,  0.08382,  0.33265,  1.11930  &
        , 2.57150,  3.58560,  2.65120,  2.28620,  2.28410,  2.19680  &
        , 2.12090,  2.07720,  2.03170,  1.99290,  1.94290/

        DATA(qscat(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        / 0.00034,  0.00236,  0.01551,  0.08981,  0.35143,  1.16480  &
        , 2.63070,  3.57600,  2.61670,  2.26380,  2.27860,  2.20140  &
        , 2.11490,  2.07450,  2.03220,  1.99180,  1.94460/

        DATA(qscat(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        / 0.00037,  0.00256,  0.01694,  0.09724,  0.37511,  1.22090  &
        , 2.69930,  3.56040,  2.57830,  2.25210,  2.27900,  2.19690  &
        , 2.11730,  2.07920,  2.03310,  1.99270,  1.94810/

        DATA(qscat(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        / 0.00040,  0.00283,  0.01880,  0.10674,  0.40586,  1.29070  &
        , 2.77890,  3.54370,  2.53590,  2.24060,  2.27390,  2.18430  &
        , 2.11360,  2.07670,  2.03600,  1.99530,  1.95130/

        DATA(qscat(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        / 0.00044,  0.00318,  0.02132,  0.11926,  0.44710,  1.37880  &
        , 2.87670,  3.50110,  2.49540,  2.23570,  2.25850,  2.17010  &
        , 2.12030,  2.07230,  2.03640,  1.99680,  1.95540/

        DATA(qscat(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        / 0.00050,  0.00367,  0.02490,  0.13648,  0.50423,  1.49370  &
        , 3.00670,  3.43340,  2.45820,  2.25000,  2.23680,  2.15100  &
        , 2.11880,  2.06880,  2.03560,  1.99800,  1.95940/

        DATA(qscat(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        / 0.00058,  0.00440,  0.03033,  0.16147,  0.58584,  1.65450  &
        , 3.16050,  3.31080,  2.42720,  2.29180,  2.19090,  2.15090  &
        , 2.10700,  2.06620,  2.03200,  2.00000,  1.96280/

        DATA(qscat(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        / 0.00070,  0.00560,  0.03934,  0.20080,  0.70801,  1.89810  &
        , 3.34690,  3.10610,  2.39880,  2.33890,  2.18180,  2.14890  &
        , 2.10080,  2.06290,  2.03130,  2.00120,  1.96690/

        DATA(qscat(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        / 0.00091,  0.00785,  0.05679,  0.27452,  0.91914,  2.27780  &
        , 3.51890,  2.77130,  2.32140,  2.31030,  2.19990,  2.12400  &
        , 2.08820,  2.05990,  2.03000,  2.00130,  1.97120/

        DATA(qscat(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        / 0.00135,  0.01342,  0.10041,  0.46794,  1.40620,  2.96270  &
        , 3.40270,  2.42790,  2.25410,  2.20830,  2.14420,  2.11080  &
        , 2.07720,  2.04800,  2.02300,  1.99800,  1.96980/

        DATA(qscat(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        / 0.00370,  0.02361,  0.12737,  0.52409,  1.52090,  3.05660  &
        , 3.54360,  2.43440,  2.38660,  2.19150,  2.17570,  2.10150  &
        , 2.05360,  2.01450,  1.95950,  1.89610,  1.81840/

        DATA(qscat(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        / 0.00381,  0.02424,  0.13018,  0.53257,  1.53640,  3.06980  &
        , 3.53500,  2.42200,  2.38590,  2.18550,  2.18050,  2.10420  &
        , 2.05760,  2.01580,  1.95890,  1.90010,  1.82170/

        DATA(qscat(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        / 0.00392,  0.02493,  0.13326,  0.54184,  1.55350,  3.08450  &
        , 3.52510,  2.41090,  2.38610,  2.19630,  2.17520,  2.11070  &
        , 2.05690,  2.01050,  1.96450,  1.90570,  1.82690/

        DATA(qscat(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        / 0.00404,  0.02569,  0.13665,  0.55204,  1.57240,  3.10070  &
        , 3.51370,  2.39650,  2.38550,  2.19750,  2.17540,  2.10910  &
        , 2.05660,  2.01770,  1.96570,  1.90660,  1.83050/

        DATA(qscat(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        / 0.00417,  0.02654,  0.14041,  0.56332,  1.59330,  3.11860  &
        , 3.50020,  2.38450,  2.38270,  2.20140,  2.18330,  2.11300  &
        , 2.06530,  2.01560,  1.96610,  1.90890,  1.83520/

        DATA(qscat(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        / 0.00431,  0.02749,  0.14460,  0.57590,  1.61650,  3.13790  &
        , 3.48490,  2.37100,  2.38660,  2.19740,  2.18370,  2.11060  &
        , 2.06340,  2.01500,  1.96880,  1.91180,  1.83810/

        DATA(qscat(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        / 0.00447,  0.02856,  0.14930,  0.59004,  1.64250,  3.15880  &
        , 3.46850,  2.35610,  2.38700,  2.20830,  2.18180,  2.10750  &
        , 2.06230,  2.01440,  1.97180,  1.91770,  1.84330/

        DATA(qscat(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        / 0.00466,  0.02977,  0.15464,  0.60612,  1.67190,  3.18140  &
        , 3.44970,  2.33930,  2.39140,  2.20770,  2.18500,  2.11610  &
        , 2.06020,  2.01440,  1.97320,  1.91970,  1.85000/

        DATA(qscat(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        / 0.00487,  0.03117,  0.16077,  0.62462,  1.70530,  3.20590  &
        , 3.42770,  2.31890,  2.39110,  2.20940,  2.17870,  2.10990  &
        , 2.05930,  2.02190,  1.97380,  1.92200,  1.85550/

        DATA(qscat(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        / 0.00511,  0.03280,  0.16790,  0.64622,  1.74370,  3.23300  &
        , 3.40210,  2.29520,  2.39340,  2.21990,  2.16730,  2.11600  &
        , 2.05750,  2.01800,  1.97630,  1.92700,  1.86160/

        DATA(qscat(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        / 0.00540,  0.03473,  0.17633,  0.67183,  1.78850,  3.26370  &
        , 3.36980,  2.27110,  2.39520,  2.22330,  2.17140,  2.11280  &
        , 2.06020,  2.02410,  1.98000,  1.93110,  1.86730/

        DATA(qscat(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        / 0.00574,  0.03706,  0.18648,  0.70277,  1.84150,  3.29920  &
        , 3.32590,  2.25000,  2.38820,  2.23320,  2.16480,  2.10000  &
        , 2.06100,  2.02150,  1.98210,  1.93590,  1.87530/

        DATA(qscat(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        / 0.00615,  0.03992,  0.19901,  0.74093,  1.90540,  3.33920  &
        , 3.27130,  2.23230,  2.39080,  2.24020,  2.15510,  2.09600  &
        , 2.06130,  2.02240,  1.98540,  1.94110,  1.88370/

        DATA(qscat(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        / 0.00667,  0.04354,  0.21491,  0.78901,  1.98350,  3.38410  &
        , 3.21430,  2.20550,  2.38790,  2.23870,  2.14650,  2.09180  &
        , 2.06140,  2.02490,  1.98850,  1.94670,  1.89080/

        DATA(qscat(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        / 0.00734,  0.04828,  0.23580,  0.85103,  2.08060,  3.43590  &
        , 3.12910,  2.19330,  2.37820,  2.24250,  2.13320,  2.09930  &
        , 2.06350,  2.02560,  1.99210,  1.95240,  1.90110/

        DATA(qscat(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        / 0.00825,  0.05474,  0.26441,  0.93331,  2.20550,  3.48960  &
        , 3.01550,  2.19110,  2.35350,  2.22710,  2.12380,  2.10170  &
        , 2.06060,  2.02880,  1.99580,  1.96050,  1.91330/

        DATA(qscat(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        / 0.00952,  0.06404,  0.30561,  1.04730,  2.37530,  3.54300  &
        , 2.86350,  2.21810,  2.29920,  2.19190,  2.13960,  2.08910  &
        , 2.05470,  2.02890,  2.00010,  1.96710,  1.92490/

        DATA(qscat(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        / 0.01147,  0.07844,  0.36883,  1.21810,  2.61190,  3.55590  &
        , 2.64970,  2.28200,  2.20870,  2.15350,  2.13070,  2.09220  &
        , 2.05370,  2.02950,  2.00410,  1.97700,  1.94070/

        DATA(qscat(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        / 0.01477,  0.10349,  0.47814,  1.50780,  2.95820,  3.46930  &
        , 2.37220,  2.36180,  2.19230,  2.17480,  2.11580,  2.08470  &
        , 2.05500,  2.03080,  2.00970,  1.98670,  1.95780/

        DATA(qscat(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        / 0.02148,  0.15799,  0.72082,  2.06630,  3.41310,  3.02420  &
        , 2.18390,  2.33270,  2.21370,  2.13040,  2.10130,  2.07200  &
        , 2.04820,  2.02980,  2.01440,  1.99660,  1.97740/

        DATA(qscat(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00006,  0.00040,  0.00258,  0.01637,  0.09695,  0.40877  &
        , 0.99543,  1.43050,  1.31110,  1.17900,  1.17950/

        DATA(qscat(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00007,  0.00042,  0.00270,  0.01708,  0.10063,  0.41801  &
        , 1.00260,  1.42490,  1.30340,  1.17870,  1.17940/

        DATA(qscat(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00007,  0.00044,  0.00282,  0.01786,  0.10471,  0.42801  &
        , 1.01030,  1.41930,  1.29560,  1.17840,  1.17920/

        DATA(qscat(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00007,  0.00046,  0.00297,  0.01875,  0.10925,  0.43886  &
        , 1.01860,  1.41350,  1.28770,  1.17820,  1.17910/

        DATA(qscat(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00008,  0.00049,  0.00313,  0.01975,  0.11432,  0.45070  &
        , 1.02780,  1.40750,  1.27980,  1.17790,  1.17890/

        DATA(qscat(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00008,  0.00052,  0.00331,  0.02090,  0.12004,  0.46370  &
        , 1.03770,  1.40140,  1.27170,  1.17760,  1.17880/

        DATA(qscat(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00009,  0.00055,  0.00353,  0.02222,  0.12654,  0.47806  &
        , 1.04860,  1.39500,  1.26360,  1.17740,  1.17870/

        DATA(qscat(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00009,  0.00059,  0.00378,  0.02376,  0.13401,  0.49404  &
        , 1.06060,  1.38840,  1.25530,  1.17710,  1.17860/

        DATA(qscat(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00010,  0.00064,  0.00408,  0.02560,  0.14268,  0.51197  &
        , 1.07390,  1.38150,  1.24700,  1.17690,  1.17850/

        DATA(qscat(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00011,  0.00070,  0.00444,  0.02781,  0.15289,  0.53231  &
        , 1.08850,  1.37410,  1.23850,  1.17670,  1.17840/

        DATA(qscat(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00012,  0.00077,  0.00489,  0.03052,  0.16509,  0.55564  &
        , 1.10470,  1.36610,  1.23000,  1.17660,  1.17830/

        DATA(qscat(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00013,  0.00086,  0.00546,  0.03395,  0.17991,  0.58274  &
        , 1.12280,  1.35710,  1.22140,  1.17650,  1.17830/

        DATA(qscat(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00015,  0.00097,  0.00620,  0.03838,  0.19830,  0.61470  &
        , 1.14300,  1.34680,  1.21280,  1.17650,  1.17820/

        DATA(qscat(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00003  &
        , 0.00018,  0.00113,  0.00721,  0.04432,  0.22164,  0.65301  &
        , 1.16570,  1.33460,  1.20430,  1.17650,  1.17820/

        DATA(qscat(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00003  &
        , 0.00021,  0.00135,  0.00864,  0.05262,  0.25205,  0.69977  &
        , 1.19120,  1.32000,  1.19600,  1.17670,  1.17820/

        DATA(qscat(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00004  &
        , 0.00026,  0.00169,  0.01080,  0.06483,  0.29286,  0.75805  &
        , 1.21990,  1.30230,  1.18830,  1.17700,  1.17820/

        DATA(qscat(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00005  &
        , 0.00035,  0.00225,  0.01431,  0.08406,  0.34964,  0.83237  &
        , 1.25070,  1.27990,  1.18180,  1.17750,  1.17820/

        DATA(qscat(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00008  &
        , 0.00051,  0.00329,  0.02086,  0.11762,  0.43356,  0.93051  &
        , 1.27950,  1.25120,  1.17730,  1.17810,  1.17810/

        DATA(qscat(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00013  &
        , 0.00087,  0.00571,  0.03568,  0.18469,  0.56813,  1.06060  &
        , 1.29550,  1.21580,  1.17580,  1.17870,  1.17760/

        DATA(qscat(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00005,  0.00033  &
        , 0.00223,  0.01470,  0.08736,  0.35672,  0.82620,  1.22560  &
        , 1.26460,  1.18190,  1.17780,  1.17880,  1.17600/

        DATA(qscat(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00005  &
        , 0.00032,  0.00205,  0.01278,  0.07358,  0.30598,  0.76319  &
        , 1.24540,  1.33480,  1.16390,  1.15450,  1.15530/

        DATA(qscat(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00005  &
        , 0.00033,  0.00212,  0.01319,  0.07550,  0.30950,  0.76195  &
        , 1.23310,  1.32040,  1.16170,  1.15330,  1.15440/

        DATA(qscat(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00005  &
        , 0.00035,  0.00219,  0.01364,  0.07762,  0.31339,  0.76107  &
        , 1.22100,  1.30590,  1.15940,  1.15210,  1.15360/

        DATA(qscat(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00006  &
        , 0.00036,  0.00228,  0.01415,  0.08000,  0.31770,  0.76063  &
        , 1.20900,  1.29150,  1.15700,  1.15100,  1.15270/

        DATA(qscat(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00006  &
        , 0.00038,  0.00238,  0.01473,  0.08267,  0.32251,  0.76070  &
        , 1.19720,  1.27720,  1.15450,  1.14980,  1.15190/

        DATA(qscat(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00006  &
        , 0.00039,  0.00249,  0.01539,  0.08570,  0.32790,  0.76136  &
        , 1.18560,  1.26290,  1.15200,  1.14860,  1.15110/

        DATA(qscat(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00006  &
        , 0.00041,  0.00262,  0.01616,  0.08916,  0.33398,  0.76273  &
        , 1.17420,  1.24860,  1.14950,  1.14740,  1.15020/

        DATA(qscat(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00007  &
        , 0.00044,  0.00277,  0.01706,  0.09315,  0.34090,  0.76494  &
        , 1.16300,  1.23440,  1.14690,  1.14630,  1.14940/

        DATA(qscat(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00007  &
        , 0.00047,  0.00295,  0.01813,  0.09782,  0.34885,  0.76817  &
        , 1.15220,  1.22020,  1.14420,  1.14520,  1.14860/

        DATA(qscat(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00008  &
        , 0.00050,  0.00318,  0.01942,  0.10334,  0.35807,  0.77265  &
        , 1.14180,  1.20610,  1.14160,  1.14410,  1.14780/

        DATA(qscat(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00009  &
        , 0.00055,  0.00345,  0.02100,  0.10998,  0.36892,  0.77869  &
        , 1.13200,  1.19210,  1.13890,  1.14300,  1.14700/

        DATA(qscat(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00009  &
        , 0.00060,  0.00380,  0.02300,  0.11810,  0.38185,  0.78669  &
        , 1.12280,  1.17840,  1.13640,  1.14200,  1.14630/

        DATA(qscat(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00011  &
        , 0.00067,  0.00425,  0.02559,  0.12824,  0.39756,  0.79719  &
        , 1.11460,  1.16510,  1.13400,  1.14110,  1.14570/

        DATA(qscat(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00012  &
        , 0.00077,  0.00487,  0.02906,  0.14122,  0.41704,  0.81093  &
        , 1.10740,  1.15240,  1.13180,  1.14040,  1.14510/

        DATA(qscat(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00014  &
        , 0.00091,  0.00575,  0.03389,  0.15829,  0.44180,  0.82891  &
        , 1.10160,  1.14070,  1.13010,  1.13990,  1.14470/

        DATA(qscat(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00003,  0.00017  &
        , 0.00112,  0.00706,  0.04098,  0.18143,  0.47424,  0.85251  &
        , 1.09740,  1.13020,  1.12900,  1.13970,  1.14440/

        DATA(qscat(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00023  &
        , 0.00147,  0.00920,  0.05210,  0.21396,  0.51833,  0.88359  &
        , 1.09510,  1.12190,  1.12890,  1.13990,  1.14420/

        DATA(qscat(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00005,  0.00033  &
        , 0.00211,  0.01313,  0.07138,  0.26245,  0.58175,  0.92502  &
        , 1.09500,  1.11670,  1.13000,  1.14070,  1.14430/

        DATA(qscat(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00008,  0.00055  &
        , 0.00359,  0.02190,  0.10983,  0.34045,  0.67795,  0.97980  &
        , 1.09730,  1.11630,  1.13290,  1.14220,  1.14450/

        DATA(qscat(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        / 0.00000,  0.00000,  0.00000,  0.00003,  0.00019,  0.00135  &
        , 0.00889,  0.05146,  0.20959,  0.49717,  0.83848,  1.04940  &
        , 1.10350,  1.12400,  1.13850,  1.14430,  1.14420/

        DATA(qscat(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00021  &
        , 0.00129,  0.00772,  0.03990,  0.14049,  0.31797,  0.58959  &
        , 0.88305,  1.07650,  1.09900,  1.09250,  1.10260/

        DATA(qscat(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00021  &
        , 0.00130,  0.00775,  0.04003,  0.14069,  0.32013,  0.59639  &
        , 0.89451,  1.08540,  1.09920,  1.09150,  1.10110/

        DATA(qscat(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00021  &
        , 0.00131,  0.00780,  0.04026,  0.14114,  0.32298,  0.60446  &
        , 0.90762,  1.09530,  1.09920,  1.09060,  1.09960/

        DATA(qscat(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00021  &
        , 0.00132,  0.00788,  0.04058,  0.14188,  0.32666,  0.61404  &
        , 0.92265,  1.10620,  1.09920,  1.08990,  1.09820/

        DATA(qscat(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00022  &
        , 0.00134,  0.00798,  0.04103,  0.14298,  0.33135,  0.62541  &
        , 0.93994,  1.11810,  1.09900,  1.08930,  1.09690/

        DATA(qscat(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00003,  0.00022  &
        , 0.00136,  0.00811,  0.04163,  0.14451,  0.33725,  0.63897  &
        , 0.95988,  1.13110,  1.09860,  1.08890,  1.09560/

        DATA(qscat(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00022  &
        , 0.00139,  0.00828,  0.04243,  0.14658,  0.34466,  0.65519  &
        , 0.98300,  1.14540,  1.09790,  1.08870,  1.09440/

        DATA(qscat(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00023  &
        , 0.00143,  0.00851,  0.04348,  0.14932,  0.35396,  0.67472  &
        , 1.00990,  1.16090,  1.09700,  1.08860,  1.09340/

        DATA(qscat(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00024  &
        , 0.00148,  0.00881,  0.04485,  0.15294,  0.36565,  0.69844  &
        , 1.04140,  1.17750,  1.09570,  1.08870,  1.09240/

        DATA(qscat(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00025  &
        , 0.00155,  0.00919,  0.04664,  0.15770,  0.38040,  0.72754  &
        , 1.07850,  1.19520,  1.09410,  1.08890,  1.09150/

        DATA(qscat(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00026  &
        , 0.00164,  0.00971,  0.04901,  0.16400,  0.39919,  0.76364  &
        , 1.12240,  1.21330,  1.09220,  1.08900,  1.09070/

        DATA(qscat(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00028  &
        , 0.00176,  0.01041,  0.05219,  0.17240,  0.42335,  0.80902  &
        , 1.17480,  1.23130,  1.09030,  1.08890,  1.09010/

        DATA(qscat(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00005,  0.00031  &
        , 0.00192,  0.01137,  0.05652,  0.18383,  0.45493,  0.86679  &
        , 1.23780,  1.24770,  1.08890,  1.08860,  1.08960/

        DATA(qscat(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00005,  0.00034  &
        , 0.00216,  0.01275,  0.06256,  0.19974,  0.49704,  0.94132  &
        , 1.31340,  1.25900,  1.08910,  1.08810,  1.08940/

        DATA(qscat(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00006,  0.00040  &
        , 0.00251,  0.01478,  0.07125,  0.22270,  0.55489,  1.03920  &
        , 1.40190,  1.26110,  1.09250,  1.08810,  1.08960/

        DATA(qscat(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00007,  0.00049  &
        , 0.00306,  0.01793,  0.08422,  0.25748,  0.63760,  1.17120  &
        , 1.50380,  1.24580,  1.10130,  1.09000,  1.09000/

        DATA(qscat(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00010,  0.00063  &
        , 0.00399,  0.02318,  0.10448,  0.31369,  0.76106,  1.35240  &
        , 1.60720,  1.20440,  1.11490,  1.09430,  1.09030/

        DATA(qscat(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00014,  0.00091  &
        , 0.00576,  0.03292,  0.13868,  0.41255,  0.95369,  1.60010  &
        , 1.67930,  1.13790,  1.12120,  1.09590,  1.09110/

        DATA(qscat(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        / 0.00000,  0.00000,  0.00000,  0.00004,  0.00023,  0.00156  &
        , 0.00986,  0.05426,  0.20376,  0.59500,  1.28030,  1.92210  &
        , 1.61610,  1.10790,  1.09990,  1.09610,  1.09140/

        DATA(qscat(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        / 0.00000,  0.00000,  0.00001,  0.00008,  0.00056,  0.00390  &
        , 0.02450,  0.11934,  0.39030,  1.02070,  1.89770,  2.14610  &
        , 1.22040,  1.19920,  1.11660,  1.09590,  1.09070/

        DATA(qscat(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00055  &
        , 0.00349,  0.02210,  0.13083,  0.54750,  1.31360,  1.84030  &
        , 1.53760,  1.19340,  1.22340,  1.20190,  1.19400/

        DATA(qscat(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00055  &
        , 0.00350,  0.02215,  0.13060,  0.54204,  1.30020,  1.83330  &
        , 1.54310,  1.18950,  1.22110,  1.19890,  1.19070/

        DATA(qscat(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00055  &
        , 0.00351,  0.02220,  0.13040,  0.53638,  1.28630,  1.82600  &
        , 1.54890,  1.18600,  1.21870,  1.19570,  1.18730/

        DATA(qscat(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00055  &
        , 0.00353,  0.02226,  0.13021,  0.53052,  1.27200,  1.81830  &
        , 1.55500,  1.18270,  1.21610,  1.19250,  1.18370/

        DATA(qscat(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00056  &
        , 0.00354,  0.02234,  0.13005,  0.52445,  1.25720,  1.81030  &
        , 1.56160,  1.17980,  1.21330,  1.18910,  1.17990/

        DATA(qscat(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00056  &
        , 0.00356,  0.02244,  0.12993,  0.51816,  1.24210,  1.80190  &
        , 1.56860,  1.17710,  1.21030,  1.18550,  1.17600/

        DATA(qscat(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00056  &
        , 0.00359,  0.02256,  0.12988,  0.51166,  1.22640,  1.79320  &
        , 1.57610,  1.17450,  1.20690,  1.18170,  1.17190/

        DATA(qscat(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00057  &
        , 0.00362,  0.02271,  0.12990,  0.50496,  1.21030,  1.78420  &
        , 1.58420,  1.17190,  1.20320,  1.17760,  1.16760/

        DATA(qscat(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00058  &
        , 0.00366,  0.02290,  0.13003,  0.49807,  1.19380,  1.77490  &
        , 1.59290,  1.16900,  1.19920,  1.17330,  1.16310/

        DATA(qscat(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00058  &
        , 0.00371,  0.02314,  0.13033,  0.49106,  1.17690,  1.76540  &
        , 1.60210,  1.16620,  1.19470,  1.16860,  1.15830/

        DATA(qscat(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00059  &
        , 0.00377,  0.02346,  0.13085,  0.48401,  1.15980,  1.75590  &
        , 1.61160,  1.16390,  1.18980,  1.16350,  1.15330/

        DATA(qscat(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00010,  0.00061  &
        , 0.00385,  0.02389,  0.13171,  0.47712,  1.14270,  1.74660  &
        , 1.62110,  1.16260,  1.18430,  1.15790,  1.14800/

        DATA(qscat(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00010,  0.00062  &
        , 0.00396,  0.02448,  0.13307,  0.47071,  1.12630,  1.73790  &
        , 1.63020,  1.16210,  1.17830,  1.15190,  1.14230/

        DATA(qscat(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00010,  0.00065  &
        , 0.00412,  0.02531,  0.13521,  0.46546,  1.11160,  1.73040  &
        , 1.63890,  1.16070,  1.17170,  1.14530,  1.13620/

        DATA(qscat(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00011,  0.00069  &
        , 0.00434,  0.02655,  0.13859,  0.46260,  1.10080,  1.72560  &
        , 1.64760,  1.15870,  1.16470,  1.13830,  1.12980/

        DATA(qscat(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00012,  0.00074  &
        , 0.00469,  0.02846,  0.14405,  0.46462,  1.09800,  1.72720  &
        , 1.65680,  1.15630,  1.15700,  1.13090,  1.12290/

        DATA(qscat(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00013,  0.00083  &
        , 0.00528,  0.03161,  0.15316,  0.47666,  1.11120,  1.74310  &
        , 1.66270,  1.15010,  1.14860,  1.12310,  1.11560/

        DATA(qscat(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        / 0.00000,  0.00000,  0.00000,  0.00002,  0.00015,  0.00101  &
        , 0.00636,  0.03743,  0.16944,  0.51027,  1.15510,  1.78530  &
        , 1.65830,  1.13830,  1.13840,  1.11500,  1.10780/

        DATA(qscat(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        / 0.00000,  0.00000,  0.00000,  0.00003,  0.00021,  0.00139  &
        , 0.00881,  0.05001,  0.20201,  0.59089,  1.26440,  1.87270  &
        , 1.61450,  1.12150,  1.12270,  1.10590,  1.09980/

        DATA(qscat(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        / 0.00000,  0.00000,  0.00001,  0.00006,  0.00040,  0.00274  &
        , 0.01728,  0.08856,  0.29694,  0.80396,  1.56750,  2.00610  &
        , 1.39930,  1.13830,  1.10260,  1.09580,  1.09100/

        DATA(qscat(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        / 0.00000,  0.00000,  0.00001,  0.00003,  0.00020,  0.00127  &
        , 0.00795,  0.04691,  0.21397,  0.73696,  1.76500,  2.77830  &
        , 2.19980,  1.68950,  1.49460,  1.30100,  1.19360/

        DATA(qscat(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00021,  0.00131  &
        , 0.00823,  0.04847,  0.21909,  0.75028,  1.78770,  2.78750  &
        , 2.18260,  1.68670,  1.48990,  1.29800,  1.19090/

        DATA(qscat(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00021,  0.00137  &
        , 0.00855,  0.05020,  0.22472,  0.76471,  1.81230,  2.79780  &
        , 2.16260,  1.68310,  1.48390,  1.29490,  1.18820/

        DATA(qscat(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00022,  0.00142  &
        , 0.00891,  0.05212,  0.23092,  0.78043,  1.83890,  2.80940  &
        , 2.13900,  1.67730,  1.47790,  1.29190,  1.18550/

        DATA(qscat(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00023,  0.00149  &
        , 0.00932,  0.05428,  0.23781,  0.79767,  1.86770,  2.82180  &
        , 2.11150,  1.67120,  1.47050,  1.28870,  1.18270/

        DATA(qscat(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00025,  0.00156  &
        , 0.00978,  0.05672,  0.24552,  0.81676,  1.89880,  2.83410  &
        , 2.08150,  1.66620,  1.46080,  1.28560,  1.18010/

        DATA(qscat(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        / 0.00000,  0.00000,  0.00001,  0.00004,  0.00026,  0.00165  &
        , 0.01031,  0.05951,  0.25423,  0.83811,  1.93240,  2.84510  &
        , 2.05180,  1.66040,  1.45040,  1.28190,  1.17740/

        DATA(qscat(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        / 0.00000,  0.00000,  0.00001,  0.00005,  0.00027,  0.00175  &
        , 0.01093,  0.06272,  0.26420,  0.86231,  1.96910,  2.85330  &
        , 2.02220,  1.65020,  1.43740,  1.27760,  1.17470/

        DATA(qscat(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        / 0.00000,  0.00000,  0.00001,  0.00005,  0.00029,  0.00187  &
        , 0.01166,  0.06647,  0.27576,  0.89016,  2.00950,  2.85810  &
        , 1.98780,  1.63960,  1.42260,  1.27220,  1.17160/

        DATA(qscat(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        / 0.00000,  0.00000,  0.00001,  0.00005,  0.00031,  0.00201  &
        , 0.01254,  0.07091,  0.28940,  0.92279,  2.05500,  2.86030  &
        , 1.94450,  1.62840,  1.40570,  1.26520,  1.16780/

        DATA(qscat(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        / 0.00000,  0.00000,  0.00001,  0.00006,  0.00034,  0.00219  &
        , 0.01361,  0.07627,  0.30579,  0.96181,  2.10760,  2.86290  &
        , 1.89810,  1.61160,  1.38670,  1.25650,  1.16290/

        DATA(qscat(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        / 0.00000,  0.00000,  0.00001,  0.00006,  0.00038,  0.00241  &
        , 0.01495,  0.08286,  0.32599,  1.00950,  2.17050,  2.86600  &
        , 1.85310,  1.59260,  1.36720,  1.24610,  1.15710/

        DATA(qscat(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        / 0.00000,  0.00000,  0.00001,  0.00007,  0.00042,  0.00269  &
        , 0.01669,  0.09116,  0.35160,  1.06870,  2.24680,  2.85910  &
        , 1.79820,  1.56910,  1.34960,  1.23570,  1.15180/

        DATA(qscat(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        / 0.00000,  0.00000,  0.00001,  0.00008,  0.00048,  0.00307  &
        , 0.01901,  0.10191,  0.38515,  1.14330,  2.33680,  2.83050  &
        , 1.74450,  1.54060,  1.33760,  1.22790,  1.14800/

        DATA(qscat(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        / 0.00000,  0.00000,  0.00001,  0.00009,  0.00056,  0.00361  &
        , 0.02225,  0.11631,  0.43070,  1.23760,  2.43810,  2.78930  &
        , 1.68920,  1.50920,  1.33490,  1.22350,  1.14300/

        DATA(qscat(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        / 0.00000,  0.00000,  0.00002,  0.00011,  0.00068,  0.00440  &
        , 0.02701,  0.13640,  0.49474,  1.35740,  2.55680,  2.70870  &
        , 1.64240,  1.48130,  1.33920,  1.21320,  1.13460/

        DATA(qscat(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        / 0.00000,  0.00000,  0.00002,  0.00015,  0.00088,  0.00569  &
        , 0.03457,  0.16596,  0.58704,  1.51790,  2.70440,  2.57020  &
        , 1.61780,  1.46720,  1.32650,  1.18850,  1.12780/

        DATA(qscat(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        / 0.00000,  0.00000,  0.00003,  0.00020,  0.00123,  0.00802  &
        , 0.04794,  0.21362,  0.72497,  1.76140,  2.85170,  2.32710  &
        , 1.61980,  1.46160,  1.27480,  1.17700,  1.11950/

        DATA(qscat(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        / 0.00000,  0.00001,  0.00004,  0.00032,  0.00201,  0.01322  &
        , 0.07593,  0.30470,  0.96070,  2.12400,  2.94730,  1.95050  &
        , 1.59920,  1.38150,  1.24510,  1.15310,  1.11010/

        DATA(qscat(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        / 0.00000,  0.00001,  0.00009,  0.00070,  0.00463,  0.03081  &
        , 0.15702,  0.56957,  1.50370,  2.72180,  2.63180,  1.61510  &
        , 1.46860,  1.32380,  1.18470,  1.12540,  1.09980/

! Qscat for Sea Salt:
        DATA(qscat(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00003, 0.00023, 0.00149, 0.00939, 0.04968, 0.22147  &
        ,0.77965, 1.94100, 3.20120, 2.66530, 2.25360, 2.08700  &
        ,1.94680, 1.87010, 1.79000, 1.71870, 1.63600/

        DATA(qscat(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00004, 0.00024, 0.00156, 0.00984, 0.05183, 0.22852  &
        ,0.79833, 1.97240, 3.21400, 2.63550, 2.25630, 2.08880  &
        ,1.95830, 1.86100, 1.79370, 1.71480, 1.63280/

        DATA(qscat(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00004, 0.00025, 0.00164, 0.01034, 0.05417, 0.23611  &
        ,0.81821, 2.00500, 3.22470, 2.60610, 2.24820, 2.07650  &
        ,1.95460, 1.85800, 1.78920, 1.71450, 1.63000/

        DATA(qscat(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00004, 0.00027, 0.00173, 0.01088, 0.05672, 0.24437  &
        ,0.83960, 2.03930, 3.23330, 2.57200, 2.25660, 2.07190  &
        ,1.95030, 1.85180, 1.78940, 1.71120, 1.62450/

        DATA(qscat(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00004, 0.00028, 0.00182, 0.01148, 0.05954, 0.25345  &
        ,0.86288, 2.07580, 3.24040, 2.53150, 2.25600, 2.06590  &
        ,1.96800, 1.85280, 1.78620, 1.70630, 1.62150/

        DATA(qscat(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00004, 0.00030, 0.00193, 0.01216, 0.06268, 0.26355  &
        ,0.88858, 2.11520, 3.24720, 2.48740, 2.24190, 2.05630  &
        ,1.96270, 1.85210, 1.77880, 1.70090, 1.61750/

        DATA(qscat(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00005, 0.00032, 0.00206, 0.01293, 0.06622, 0.27494  &
        ,0.91733, 2.15860, 3.25540, 2.44900, 2.23810, 2.05370  &
        ,1.95530, 1.85060, 1.77410, 1.69820, 1.61130/

        DATA(qscat(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00005, 0.00034, 0.00220, 0.01382, 0.07027, 0.28797  &
        ,0.94994, 2.20710, 3.26490, 2.41120, 2.23060, 2.03180  &
        ,1.94810, 1.84910, 1.77440, 1.69360, 1.60810/

        DATA(qscat(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00005, 0.00036, 0.00237, 0.01487, 0.07495, 0.30310  &
        ,0.98740, 2.26210, 3.27140, 2.36400, 2.21510, 2.01950  &
        ,1.93600, 1.84180, 1.77430, 1.68870, 1.60190/

        DATA(qscat(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00006, 0.00039, 0.00256, 0.01611, 0.08042, 0.32095  &
        ,1.03090, 2.32420, 3.26950, 2.31350, 2.19120, 2.01160  &
        ,1.92730, 1.83760, 1.76450, 1.68410, 1.59630/

        DATA(qscat(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00006, 0.00043, 0.00281, 0.01762, 0.08693, 0.34238  &
        ,1.08180, 2.39320, 3.25850, 2.27510, 2.15490, 1.98210  &
        ,1.89400, 1.82910, 1.75840, 1.67620, 1.59180/

        DATA(qscat(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00007, 0.00047, 0.00310, 0.01948, 0.09478, 0.36857  &
        ,1.14150, 2.46850, 3.24460, 2.22580, 2.11910, 1.96730  &
        ,1.89710, 1.83400, 1.75170, 1.67160, 1.58420/

        DATA(qscat(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00008, 0.00052, 0.00348, 0.02183, 0.10442, 0.40118  &
        ,1.21190, 2.55100, 3.22300, 2.19160, 2.07350, 1.95850  &
        ,1.88940, 1.82850, 1.74540, 1.66310, 1.57630/

        DATA(qscat(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00009, 0.00060, 0.00398, 0.02488, 0.11652, 0.44251  &
        ,1.29540, 2.64610, 3.17240, 2.15710, 2.03900, 1.96080  &
        ,1.89320, 1.81750, 1.74010, 1.65410, 1.56660/

        DATA(qscat(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00010, 0.00069, 0.00464, 0.02899, 0.13215, 0.49579  &
        ,1.39680, 2.76130, 3.10740, 2.14490, 1.99140, 1.97110  &
        ,1.88210, 1.80270, 1.72890, 1.64440, 1.55620/

        DATA(qscat(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00011, 0.00082, 0.00559, 0.03478, 0.15313, 0.56556  &
        ,1.52720, 2.88630, 2.99490, 2.15750, 1.98940, 1.96780  &
        ,1.87020, 1.79730, 1.71910, 1.63360, 1.54360/

        DATA(qscat(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00014, 0.00102, 0.00704, 0.04349, 0.18298, 0.65921  &
        ,1.70550, 3.02210, 2.83180, 2.19030, 2.03520, 1.93300  &
        ,1.86740, 1.78680, 1.70660, 1.61860, 1.52670/

        DATA(qscat(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00018, 0.00135, 0.00949, 0.05791, 0.22988, 0.79447  &
        ,1.94980, 3.15690, 2.58510, 2.22550, 2.06930, 1.93750  &
        ,1.84350, 1.76800, 1.68880, 1.60050, 1.51170/

        DATA(qscat(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.00025, 0.00200, 0.01443, 0.08564, 0.31960, 1.03340  &
        ,2.32480, 3.22770, 2.27650, 2.14850, 1.98140, 1.89820  &
        ,1.81540, 1.74170, 1.66240, 1.57390, 1.48310/

        DATA(qscat(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.00041, 0.00372, 0.02853, 0.15654, 0.56295, 1.55470  &
        ,2.91620, 2.92000, 2.16930, 2.00640, 1.94160, 1.85940  &
        ,1.78410, 1.70410, 1.61750, 1.52570, 1.43480/

        DATA(qscat(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.00072, 0.00477, 0.02963, 0.15046, 0.53715, 1.54660  &
        ,3.05940, 3.42630, 2.44470, 2.26150, 2.22730, 2.15380  &
        ,2.12090, 2.08190, 2.05890, 2.03700, 2.01700/

        DATA(qscat(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.00076, 0.00501, 0.03099, 0.15589, 0.55370, 1.57750  &
        ,3.09000, 3.40460, 2.43820, 2.27160, 2.22420, 2.15200  &
        ,2.11780, 2.08800, 2.05630, 2.03580, 2.01610/

        DATA(qscat(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.00080, 0.00526, 0.03246, 0.16170, 0.57131, 1.61050  &
        ,3.12060, 3.37870, 2.43390, 2.27720, 2.21310, 2.15960  &
        ,2.12290, 2.09040, 2.05620, 2.03510, 2.01490/

        DATA(qscat(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.00084, 0.00553, 0.03407, 0.16799, 0.59015, 1.64610  &
        ,3.15130, 3.34920, 2.43020, 2.28990, 2.19800, 2.15550  &
        ,2.11540, 2.08230, 2.05730, 2.03420, 2.01430/

        DATA(qscat(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.00089, 0.00583, 0.03584, 0.17484, 0.61047, 1.68470  &
        ,3.18240, 3.31810, 2.42230, 2.30020, 2.19300, 2.15860  &
        ,2.11740, 2.08240, 2.05500, 2.03310, 2.01400/

        DATA(qscat(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.00094, 0.00617, 0.03782, 0.18240, 0.63257, 1.72710  &
        ,3.21490, 3.28710, 2.41990, 2.30660, 2.18880, 2.15850  &
        ,2.11360, 2.08300, 2.05380, 2.03290, 2.01300/

        DATA(qscat(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.00099, 0.00655, 0.04006, 0.19083, 0.65683, 1.77380  &
        ,3.24990, 3.25220, 2.41540, 2.31860, 2.18580, 2.16170  &
        ,2.11070, 2.07830, 2.05340, 2.03210, 2.01230/

        DATA(qscat(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.00106, 0.00698, 0.04262, 0.20034, 0.68374, 1.82540  &
        ,3.28840, 3.20860, 2.40670, 2.32910, 2.18470, 2.15730  &
        ,2.11170, 2.07770, 2.05120, 2.03050, 2.01070/

        DATA(qscat(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.00113, 0.00749, 0.04559, 0.21122, 0.71395, 1.88270  &
        ,3.33000, 3.15580, 2.40060, 2.33690, 2.18390, 2.15110  &
        ,2.10790, 2.07670, 2.05180, 2.03040, 2.00970/

        DATA(qscat(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.00122, 0.00809, 0.04908, 0.22387, 0.74833, 1.94650  &
        ,3.37200, 3.10070, 2.39400, 2.34940, 2.18590, 2.15040  &
        ,2.10630, 2.07600, 2.05060, 2.02860, 2.00890/

        DATA(qscat(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.00132, 0.00881, 0.05324, 0.23881, 0.78812, 2.01800  &
        ,3.41120, 3.03900, 2.38420, 2.34650, 2.19670, 2.15080  &
        ,2.10790, 2.07760, 2.05020, 2.02720, 2.00780/

        DATA(qscat(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.00144, 0.00969, 0.05830, 0.25683, 0.83512, 2.09920  &
        ,3.44820, 2.96140, 2.37110, 2.34970, 2.20440, 2.15580  &
        ,2.10690, 2.07290, 2.04620, 2.02620, 2.00630/

        DATA(qscat(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.00160, 0.01079, 0.06456, 0.27906, 0.89204, 2.19410  &
        ,3.48810, 2.88340, 2.35240, 2.33490, 2.20920, 2.15610  &
        ,2.10200, 2.06850, 2.04480, 2.02490, 2.00460/

        DATA(qscat(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.00179, 0.01221, 0.07250, 0.30731, 0.96288, 2.30790  &
        ,3.52540, 2.78770, 2.32800, 2.30970, 2.20170, 2.13850  &
        ,2.09820, 2.06730, 2.04360, 2.02430, 2.00360/

        DATA(qscat(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.00205, 0.01411, 0.08285, 0.34453, 1.05330, 2.44430  &
        ,3.54380, 2.68930, 2.29120, 2.28690, 2.19390, 2.13260  &
        ,2.09760, 2.06590, 2.04180, 2.02050, 1.99920/

        DATA(qscat(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.00241, 0.01674, 0.09687, 0.39566, 1.17080, 2.60430  &
        ,3.54910, 2.58460, 2.25250, 2.27190, 2.19490, 2.12740  &
        ,2.09050, 2.06060, 2.03810, 2.01890, 1.99730/

        DATA(qscat(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.00292, 0.02064, 0.11690, 0.46925, 1.32690, 2.80480  &
        ,3.49980, 2.49040, 2.22630, 2.25870, 2.17460, 2.12250  &
        ,2.08490, 2.05700, 2.03570, 2.01540, 1.99400/

        DATA(qscat(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.00373, 0.02695, 0.14781, 0.58045, 1.55090, 3.05660  &
        ,3.36790, 2.42100, 2.26320, 2.20750, 2.14610, 2.11300  &
        ,2.08050, 2.05380, 2.03080, 2.01050, 1.98990/

        DATA(qscat(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.00518, 0.03878, 0.20249, 0.76440, 1.91440, 3.35280  &
        ,3.06560, 2.38420, 2.34390, 2.18810, 2.14700, 2.10430  &
        ,2.07410, 2.04700, 2.02540, 2.00370, 1.98210/

        DATA(qscat(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.00852, 0.06786, 0.33636, 1.18000, 2.58570, 3.53280  &
        ,2.53990, 2.23140, 2.25990, 2.18940, 2.12290, 2.08720  &
        ,2.05920, 2.03590, 2.01470, 1.99480, 1.96930/

        DATA(qscat(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.01390, 0.07864, 0.33227, 1.04460, 2.28450, 3.49860  &
        ,2.98590, 2.19520, 2.34920, 2.22500, 2.14170, 2.11090  &
        ,2.08200, 2.06070, 2.04250, 2.03250, 2.02210/

        DATA(qscat(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.01464, 0.08201, 0.34275, 1.06840, 2.31620, 3.50960  &
        ,2.95390, 2.20050, 2.34220, 2.22620, 2.13850, 2.10850  &
        ,2.08410, 2.05840, 2.04400, 2.03110, 2.02240/

        DATA(qscat(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.01537, 0.08541, 0.35400, 1.09360, 2.34930, 3.52010  &
        ,2.92280, 2.20660, 2.32180, 2.21840, 2.14490, 2.10860  &
        ,2.07680, 2.05670, 2.04120, 2.03020, 2.02260/

        DATA(qscat(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.01612, 0.08908, 0.36614, 1.12070, 2.38480, 3.53090  &
        ,2.89170, 2.21720, 2.31400, 2.20760, 2.14590, 2.10120  &
        ,2.07690, 2.05470, 2.04240, 2.03030, 2.02300/

        DATA(qscat(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.01694, 0.09310, 0.37934, 1.14990, 2.42310, 3.53880  &
        ,2.85880, 2.22360, 2.29950, 2.20100, 2.14770, 2.09830  &
        ,2.07260, 2.05660, 2.04210, 2.02990, 2.02280/

        DATA(qscat(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.01784, 0.09754, 0.39382, 1.18200, 2.46470, 3.54480  &
        ,2.82230, 2.23140, 2.28320, 2.19710, 2.14640, 2.10400  &
        ,2.07280, 2.05520, 2.04100, 2.02960, 2.02140/

        DATA(qscat(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.01885, 0.10249, 0.40984, 1.21740, 2.50980, 3.54840  &
        ,2.78170, 2.24250, 2.26770, 2.18710, 2.15200, 2.09650  &
        ,2.07650, 2.05450, 2.03970, 2.02960, 2.02240/

        DATA(qscat(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.01999, 0.10807, 0.42771, 1.25700, 2.55930, 3.54920  &
        ,2.73530, 2.25490, 2.25220, 2.18260, 2.14880, 2.10060  &
        ,2.07430, 2.05380, 2.03970, 2.02990, 2.02110/

        DATA(qscat(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.02130, 0.11445, 0.44784, 1.30140, 2.61360, 3.55050  &
        ,2.68520, 2.27430, 2.23320, 2.16970, 2.14060, 2.10390  &
        ,2.07340, 2.05340, 2.03930, 2.02810, 2.02010/

        DATA(qscat(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.02282, 0.12184, 0.47072, 1.35140, 2.67260, 3.54930  &
        ,2.63090, 2.28640, 2.21800, 2.16130, 2.12830, 2.09900  &
        ,2.07080, 2.05150, 2.03840, 2.02900, 2.02020/

        DATA(qscat(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.02461, 0.13049, 0.49706, 1.40850, 2.73990, 3.54170  &
        ,2.57320, 2.30470, 2.20080, 2.15640, 2.13060, 2.09570  &
        ,2.06830, 2.05220, 2.03740, 2.02830, 2.02120/

        DATA(qscat(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.02676, 0.14077, 0.52789, 1.47490, 2.81590, 3.52420  &
        ,2.51200, 2.32930, 2.18840, 2.16050, 2.12120, 2.09040  &
        ,2.06970, 2.05100, 2.03810, 2.02760, 2.01930/

        DATA(qscat(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.02940, 0.15318, 0.56477, 1.55310, 2.90210, 3.49850  &
        ,2.44450, 2.34220, 2.18330, 2.16790, 2.11720, 2.08720  &
        ,2.06750, 2.05040, 2.03660, 2.02660, 2.02000/

        DATA(qscat(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.03270, 0.16840, 0.61014, 1.64650, 2.99490, 3.45990  &
        ,2.36890, 2.36290, 2.18670, 2.17930, 2.12060, 2.09150  &
        ,2.06570, 2.04780, 2.03540, 2.02620, 2.01870/

        DATA(qscat(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.03695, 0.18755, 0.66800, 1.75840, 3.09880, 3.39210  &
        ,2.30260, 2.37550, 2.20480, 2.17790, 2.12780, 2.08800  &
        ,2.06190, 2.04560, 2.03420, 2.02410, 2.01840/

        DATA(qscat(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.04262, 0.21252, 0.74478, 1.89610, 3.21790, 3.30210  &
        ,2.23270, 2.37720, 2.22460, 2.16150, 2.11120, 2.08170  &
        ,2.06330, 2.04430, 2.03360, 2.02530, 2.01840/

        DATA(qscat(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.05057, 0.24710, 0.84932, 2.07510, 3.34410, 3.16020  &
        ,2.18690, 2.36780, 2.23490, 2.14050, 2.10870, 2.07700  &
        ,2.05890, 2.04380, 2.03100, 2.02410, 2.01650/

        DATA(qscat(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.06256, 0.29951, 0.99765, 2.31660, 3.47200, 2.94570  &
        ,2.19220, 2.32250, 2.21080, 2.13720, 2.10380, 2.07590  &
        ,2.05420, 2.04030, 2.03000, 2.02130, 2.01610/

        DATA(qscat(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.08290, 0.38942, 1.23930, 2.67350, 3.53380, 2.61790  &
        ,2.28460, 2.20780, 2.15350, 2.12650, 2.09820, 2.06710  &
        ,2.04930, 2.03730, 2.02680, 2.01960, 2.01470/

        DATA(qscat(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.12530, 0.57315, 1.71420, 3.21310, 3.30290, 2.21790  &
        ,2.38380, 2.23300, 2.15730, 2.10010, 2.08180, 2.05830  &
        ,2.04370, 2.03130, 2.02260, 2.01610, 2.01310/

        DATA(qscat(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00024, 0.00152, 0.00959, 0.05679, 0.25333, 0.66013  &
        ,1.11390, 1.27140, 1.18610, 1.16220, 1.16420/

        DATA(qscat(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00025, 0.00162, 0.01019, 0.06010, 0.26415, 0.67730  &
        ,1.12620, 1.27060, 1.18430, 1.16280, 1.16480/

        DATA(qscat(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00004  &
        ,0.00027, 0.00173, 0.01085, 0.06372, 0.27565, 0.69517  &
        ,1.13860, 1.26930, 1.18250, 1.16350, 1.16540/

        DATA(qscat(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00029, 0.00184, 0.01159, 0.06771, 0.28796, 0.71387  &
        ,1.15090, 1.26770, 1.18080, 1.16410, 1.16610/

        DATA(qscat(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00031, 0.00198, 0.01241, 0.07215, 0.30124, 0.73358  &
        ,1.16340, 1.26570, 1.17910, 1.16480, 1.16670/

        DATA(qscat(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00005  &
        ,0.00033, 0.00213, 0.01335, 0.07714, 0.31569, 0.75450  &
        ,1.17590, 1.26320, 1.17760, 1.16550, 1.16730/

        DATA(qscat(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00006  &
        ,0.00036, 0.00230, 0.01443, 0.08281, 0.33155, 0.77685  &
        ,1.18850, 1.26020, 1.17620, 1.16630, 1.16790/

        DATA(qscat(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00006  &
        ,0.00039, 0.00251, 0.01569, 0.08933, 0.34912, 0.80090  &
        ,1.20110, 1.25670, 1.17480, 1.16710, 1.16860/

        DATA(qscat(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00007  &
        ,0.00043, 0.00275, 0.01718, 0.09692, 0.36876, 0.82693  &
        ,1.21370, 1.25250, 1.17360, 1.16790, 1.16920/

        DATA(qscat(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00007  &
        ,0.00048, 0.00304, 0.01898, 0.10588, 0.39091, 0.85525  &
        ,1.22620, 1.24760, 1.17260, 1.16870, 1.16980/

        DATA(qscat(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00008  &
        ,0.00053, 0.00340, 0.02118, 0.11660, 0.41616, 0.88620  &
        ,1.23830, 1.24190, 1.17170, 1.16960, 1.17050/

        DATA(qscat(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00060, 0.00385, 0.02393, 0.12962, 0.44525, 0.92013  &
        ,1.24990, 1.23530, 1.17110, 1.17060, 1.17110/

        DATA(qscat(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00011  &
        ,0.00070, 0.00444, 0.02746, 0.14571, 0.47917, 0.95741  &
        ,1.26060, 1.22770, 1.17080, 1.17150, 1.17170/

        DATA(qscat(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00013  &
        ,0.00082, 0.00521, 0.03212, 0.16600, 0.51934, 0.99844  &
        ,1.26970, 1.21930, 1.17080, 1.17250, 1.17230/

        DATA(qscat(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00015  &
        ,0.00099, 0.00630, 0.03852, 0.19226, 0.56782, 1.04370  &
        ,1.27620, 1.21000, 1.17120, 1.17350, 1.17290/

        DATA(qscat(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00019  &
        ,0.00124, 0.00788, 0.04775, 0.22733, 0.62776, 1.09350  &
        ,1.27870, 1.20000, 1.17200, 1.17450, 1.17330/

        DATA(qscat(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00025  &
        ,0.00163, 0.01041, 0.06207, 0.27624, 0.70427, 1.14850  &
        ,1.27500, 1.19010, 1.17320, 1.17540, 1.17370/

        DATA(qscat(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00036  &
        ,0.00235, 0.01496, 0.08672, 0.34865, 0.80602, 1.20730  &
        ,1.26180, 1.18110, 1.17490, 1.17630, 1.17370/

        DATA(qscat(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00009, 0.00061  &
        ,0.00396, 0.02502, 0.13662, 0.46683, 0.94793, 1.26120  &
        ,1.23430, 1.17550, 1.17670, 1.17680, 1.17330/

        DATA(qscat(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00021, 0.00146  &
        ,0.00967, 0.05933, 0.27152, 0.70240, 1.14990, 1.27720  &
        ,1.19220, 1.17600, 1.17830, 1.17630, 1.17150/

        DATA(qscat(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00016  &
        ,0.00102, 0.00634, 0.03632, 0.15966, 0.42250, 0.78438  &
        ,1.05280, 1.11540, 1.11680, 1.12900, 1.13530/

        DATA(qscat(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00017  &
        ,0.00109, 0.00671, 0.03828, 0.16598, 0.43241, 0.79364  &
        ,1.05510, 1.11410, 1.11710, 1.12950, 1.13570/

        DATA(qscat(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00018  &
        ,0.00115, 0.00712, 0.04042, 0.17270, 0.44281, 0.80317  &
        ,1.05740, 1.11300, 1.11750, 1.13000, 1.13600/

        DATA(qscat(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00019  &
        ,0.00123, 0.00758, 0.04278, 0.17991, 0.45380, 0.81303  &
        ,1.05970, 1.11190, 1.11790, 1.13050, 1.13640/

        DATA(qscat(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00021  &
        ,0.00131, 0.00809, 0.04540, 0.18769, 0.46551, 0.82330  &
        ,1.06210, 1.11090, 1.11850, 1.13110, 1.13680/

        DATA(qscat(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00003, 0.00022  &
        ,0.00141, 0.00867, 0.04833, 0.19615, 0.47810, 0.83407  &
        ,1.06460, 1.11010, 1.11900, 1.13170, 1.13730/

        DATA(qscat(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00024  &
        ,0.00152, 0.00934, 0.05166, 0.20546, 0.49174, 0.84542  &
        ,1.06710, 1.10940, 1.11970, 1.13230, 1.13770/

        DATA(qscat(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00026  &
        ,0.00165, 0.01011, 0.05548, 0.21576, 0.50663, 0.85747  &
        ,1.06960, 1.10880, 1.12050, 1.13300, 1.13810/

        DATA(qscat(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00004, 0.00029  &
        ,0.00181, 0.01103, 0.05992, 0.22728, 0.52304, 0.87032  &
        ,1.07230, 1.10840, 1.12130, 1.13370, 1.13860/

        DATA(qscat(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00032  &
        ,0.00199, 0.01213, 0.06515, 0.24028, 0.54125, 0.88408  &
        ,1.07500, 1.10820, 1.12230, 1.13440, 1.13910/

        DATA(qscat(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00005, 0.00035  &
        ,0.00222, 0.01347, 0.07140, 0.25507, 0.56162, 0.89887  &
        ,1.07780, 1.10820, 1.12340, 1.13530, 1.13960/

        DATA(qscat(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00040  &
        ,0.00250, 0.01515, 0.07898, 0.27209, 0.58459, 0.91482  &
        ,1.08070, 1.10850, 1.12470, 1.13610, 1.14010/

        DATA(qscat(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00046  &
        ,0.00287, 0.01728, 0.08834, 0.29191, 0.61069, 0.93207  &
        ,1.08360, 1.10920, 1.12610, 1.13710, 1.14060/

        DATA(qscat(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00053  &
        ,0.00336, 0.02008, 0.10014, 0.31535, 0.64066, 0.95078  &
        ,1.08660, 1.11030, 1.12760, 1.13810, 1.14110/

        DATA(qscat(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00010, 0.00064  &
        ,0.00403, 0.02390, 0.11542, 0.34367, 0.67546, 0.97112  &
        ,1.08960, 1.11190, 1.12940, 1.13910, 1.14160/

        DATA(qscat(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00080  &
        ,0.00502, 0.02938, 0.13589, 0.37894, 0.71650, 0.99317  &
        ,1.09270, 1.11410, 1.13150, 1.14030, 1.14210/

        DATA(qscat(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00105  &
        ,0.00657, 0.03779, 0.16456, 0.42484, 0.76588, 1.01690  &
        ,1.09610, 1.11720, 1.13380, 1.14150, 1.14250/

        DATA(qscat(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00150  &
        ,0.00934, 0.05212, 0.20727, 0.48860, 0.82691, 1.04210  &
        ,1.10020, 1.12140, 1.13650, 1.14270, 1.14280/

        DATA(qscat(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00000, 0.00000, 0.00001, 0.00006, 0.00038, 0.00248  &
        ,0.01535, 0.08081, 0.27711, 0.58578, 0.90485, 1.06790  &
        ,1.10620, 1.12710, 1.13960, 1.14380, 1.14280/

        DATA(qscat(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00000, 0.00000, 0.00002, 0.00014, 0.00088, 0.00586  &
        ,0.03514, 0.15840, 0.41680, 0.75471, 1.00690, 1.09230  &
        ,1.11740, 1.13510, 1.14320, 1.14450, 1.14170/

        DATA(qscat(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00068  &
        ,0.00424, 0.02545, 0.12815, 0.43942, 1.22370, 2.36600  &
        ,2.65160, 1.26560, 1.38750, 1.18780, 1.11420/

        DATA(qscat(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00071  &
        ,0.00444, 0.02659, 0.13253, 0.45157, 1.24550, 2.38380  &
        ,2.62430, 1.24700, 1.37230, 1.17960, 1.11400/

        DATA(qscat(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00075  &
        ,0.00466, 0.02783, 0.13723, 0.46478, 1.26870, 2.40190  &
        ,2.59320, 1.23150, 1.35570, 1.17110, 1.11360/

        DATA(qscat(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00078  &
        ,0.00491, 0.02919, 0.14231, 0.47924, 1.29350, 2.41980  &
        ,2.55870, 1.21730, 1.33760, 1.16220, 1.11280/

        DATA(qscat(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00013, 0.00083  &
        ,0.00518, 0.03071, 0.14783, 0.49521, 1.31980, 2.43720  &
        ,2.52170, 1.20200, 1.31800, 1.15320, 1.11180/

        DATA(qscat(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00014, 0.00088  &
        ,0.00549, 0.03241, 0.15391, 0.51297, 1.34800, 2.45390  &
        ,2.48310, 1.18690, 1.29660, 1.14420, 1.11040/

        DATA(qscat(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00015, 0.00094  &
        ,0.00584, 0.03435, 0.16065, 0.53286, 1.37810, 2.47020  &
        ,2.44160, 1.17650, 1.27370, 1.13550, 1.10870/

        DATA(qscat(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00100  &
        ,0.00625, 0.03659, 0.16822, 0.55527, 1.41060, 2.48740  &
        ,2.39430, 1.17070, 1.24900, 1.12770, 1.10680/

        DATA(qscat(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00017, 0.00108  &
        ,0.00673, 0.03920, 0.17680, 0.58066, 1.44600, 2.50700  &
        ,2.33790, 1.16440, 1.22330, 1.12120, 1.10500/

        DATA(qscat(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00018, 0.00118  &
        ,0.00731, 0.04230, 0.18665, 0.60951, 1.48530, 2.53050  &
        ,2.27220, 1.16350, 1.19670, 1.11660, 1.10330/

        DATA(qscat(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00020, 0.00129  &
        ,0.00802, 0.04602, 0.19811, 0.64239, 1.53030, 2.55720  &
        ,2.19910, 1.17110, 1.17100, 1.11460, 1.10200/

        DATA(qscat(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00144  &
        ,0.00890, 0.05059, 0.21164, 0.67993, 1.58350, 2.58260  &
        ,2.11420, 1.18160, 1.14790, 1.11540, 1.10130/

        DATA(qscat(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00025, 0.00162  &
        ,0.01002, 0.05630, 0.22793, 0.72299, 1.64830, 2.60080  &
        ,2.01080, 1.20470, 1.13040, 1.11820, 1.10070/

        DATA(qscat(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00000, 0.00000, 0.00001, 0.00005, 0.00029, 0.00186  &
        ,0.01150, 0.06362, 0.24811, 0.77306, 1.72740, 2.61330  &
        ,1.89310, 1.23470, 1.12230, 1.12080, 1.09980/

        DATA(qscat(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00000, 0.00000, 0.00001, 0.00006, 0.00034, 0.00220  &
        ,0.01352, 0.07331, 0.27412, 0.83336, 1.82040, 2.62670  &
        ,1.75080, 1.27270, 1.12720, 1.11960, 1.09830/

        DATA(qscat(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00000, 0.00000, 0.00001, 0.00007, 0.00042, 0.00269  &
        ,0.01642, 0.08668, 0.30970, 0.91153, 1.92660, 2.61830  &
        ,1.58800, 1.31280, 1.14470, 1.11200, 1.09690/

        DATA(qscat(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00000, 0.00000, 0.00001, 0.00009, 0.00054, 0.00345  &
        ,0.02094, 0.10614, 0.36288, 1.02410, 2.06130, 2.57390  &
        ,1.40530, 1.33520, 1.16200, 1.10310, 1.09550/

        DATA(qscat(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00000, 0.00000, 0.00002, 0.00013, 0.00075, 0.00481  &
        ,0.02879, 0.13681, 0.45204, 1.19610, 2.23820, 2.45290  &
        ,1.22660, 1.29830, 1.14790, 1.10240, 1.09390/

        DATA(qscat(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00000, 0.00000, 0.00003, 0.00020, 0.00121, 0.00776  &
        ,0.04513, 0.19164, 0.61562, 1.45610, 2.44340, 2.15890  &
        ,1.14680, 1.16680, 1.10830, 1.09820, 1.09190/

        DATA(qscat(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00000, 0.00001, 0.00006, 0.00044, 0.00273, 0.01760  &
        ,0.09380, 0.32978, 0.95046, 1.95710, 2.54660, 1.47480  &
        ,1.30660, 1.14630, 1.10370, 1.09440, 1.08890/

        DATA(qscat(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00038  &
        ,0.00236, 0.01426, 0.07517, 0.27122, 0.80264, 1.72410  &
        ,2.46210, 1.70540, 1.23030, 1.13150, 1.10550/

        DATA(qscat(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00006, 0.00039  &
        ,0.00247, 0.01488, 0.07782, 0.27758, 0.81510, 1.73760  &
        ,2.45240, 1.67410, 1.23210, 1.13080, 1.10410/

        DATA(qscat(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00041  &
        ,0.00259, 0.01556, 0.08066, 0.28450, 0.82862, 1.75220  &
        ,2.44260, 1.64180, 1.23350, 1.13020, 1.10280/

        DATA(qscat(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00043  &
        ,0.00272, 0.01631, 0.08375, 0.29209, 0.84342, 1.76800  &
        ,2.43250, 1.60920, 1.23490, 1.12960, 1.10160/

        DATA(qscat(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00007, 0.00046  &
        ,0.00286, 0.01714, 0.08713, 0.30051, 0.85977, 1.78540  &
        ,2.42190, 1.57600, 1.23620, 1.12900, 1.10050/

        DATA(qscat(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00048  &
        ,0.00303, 0.01807, 0.09087, 0.30996, 0.87799, 1.80470  &
        ,2.40990, 1.54100, 1.23670, 1.12820, 1.09960/

        DATA(qscat(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00000, 0.00000, 0.00000, 0.00001, 0.00008, 0.00052  &
        ,0.00322, 0.01914, 0.09504, 0.32068, 0.89846, 1.82600  &
        ,2.39550, 1.50360, 1.23680, 1.12730, 1.09870/

        DATA(qscat(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00009, 0.00055  &
        ,0.00344, 0.02037, 0.09975, 0.33301, 0.92161, 1.84960  &
        ,2.37760, 1.46490, 1.23600, 1.12600, 1.09790/

        DATA(qscat(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00009, 0.00059  &
        ,0.00370, 0.02181, 0.10512, 0.34734, 0.94794, 1.87540  &
        ,2.35620, 1.42550, 1.23380, 1.12430, 1.09730/

        DATA(qscat(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00010, 0.00064  &
        ,0.00401, 0.02352, 0.11132, 0.36421, 0.97796, 1.90340  &
        ,2.33180, 1.38350, 1.23010, 1.12200, 1.09660/

        DATA(qscat(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00011, 0.00071  &
        ,0.00439, 0.02559, 0.11856, 0.38428, 1.01230, 1.93400  &
        ,2.30340, 1.34000, 1.22390, 1.11890, 1.09590/

        DATA(qscat(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00012, 0.00078  &
        ,0.00486, 0.02813, 0.12713, 0.40838, 1.05190, 1.96800  &
        ,2.26680, 1.29710, 1.21470, 1.11510, 1.09500/

        DATA(qscat(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00000, 0.00000, 0.00000, 0.00002, 0.00014, 0.00088  &
        ,0.00547, 0.03133, 0.13744, 0.43759, 1.09800, 2.00720  &
        ,2.21910, 1.25350, 1.20180, 1.11070, 1.09380/

        DATA(qscat(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00016, 0.00101  &
        ,0.00626, 0.03546, 0.15012, 0.47326, 1.15310, 2.05200  &
        ,2.16040, 1.21390, 1.18450, 1.10610, 1.09230/

        DATA(qscat(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00000, 0.00000, 0.00000, 0.00003, 0.00019, 0.00119  &
        ,0.00735, 0.04098, 0.16617, 0.51722, 1.22160, 2.09940  &
        ,2.08070, 1.18030, 1.16300, 1.10240, 1.09100/

        DATA(qscat(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00000, 0.00000, 0.00001, 0.00004, 0.00023, 0.00145  &
        ,0.00892, 0.04868, 0.18739, 0.57241, 1.30980, 2.14940  &
        ,1.97550, 1.15820, 1.13940, 1.10010, 1.08990/

        DATA(qscat(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00000, 0.00000, 0.00001, 0.00005, 0.00029, 0.00186  &
        ,0.01136, 0.06010, 0.21759, 0.64532, 1.42370, 2.20050  &
        ,1.82860, 1.15440, 1.11910, 1.09840, 1.08870/

        DATA(qscat(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00000, 0.00000, 0.00001, 0.00007, 0.00040, 0.00259  &
        ,0.01563, 0.07859, 0.26627, 0.75393, 1.57760, 2.23340  &
        ,1.62510, 1.17320, 1.11030, 1.09430, 1.08710/

        DATA(qscat(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00000, 0.00000, 0.00002, 0.00011, 0.00065, 0.00417  &
        ,0.02463, 0.11278, 0.36187, 0.94196, 1.80160, 2.19720  &
        ,1.35710, 1.19040, 1.10810, 1.09080, 1.08520/

        DATA(qscat(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00000, 0.00000, 0.00003, 0.00023, 0.00146, 0.00946  &
        ,0.05239, 0.19656, 0.59175, 1.32850, 2.11640, 1.86710  &
        ,1.14870, 1.12170, 1.09590, 1.08710, 1.08250/

        DATA(qscat(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00000, 0.00000, 0.00002, 0.00014, 0.00080, 0.00507  &
        ,0.03067, 0.15390, 0.56158, 1.53360, 2.83190, 2.82510  &
        ,1.68210, 1.55970, 1.40520, 1.23670, 1.15300/

        DATA(qscat(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00000, 0.00000, 0.00002, 0.00015, 0.00084, 0.00532  &
        ,0.03209, 0.15943, 0.57857, 1.56040, 2.85230, 2.78920  &
        ,1.67970, 1.55290, 1.39780, 1.23060, 1.15060/

        DATA(qscat(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00000, 0.00000, 0.00002, 0.00016, 0.00089, 0.00559  &
        ,0.03365, 0.16536, 0.59667, 1.58870, 2.87350, 2.75320  &
        ,1.67690, 1.54710, 1.38900, 1.22480, 1.14790/

        DATA(qscat(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00000, 0.00000, 0.00003, 0.00016, 0.00094, 0.00590  &
        ,0.03536, 0.17177, 0.61606, 1.61920, 2.89480, 2.71520  &
        ,1.67180, 1.54220, 1.37890, 1.21970, 1.14500/

        DATA(qscat(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00000, 0.00000, 0.00003, 0.00017, 0.00099, 0.00623  &
        ,0.03727, 0.17878, 0.63697, 1.65220, 2.91570, 2.67260  &
        ,1.66840, 1.53760, 1.36690, 1.21540, 1.14180/

        DATA(qscat(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00000, 0.00000, 0.00003, 0.00018, 0.00105, 0.00662  &
        ,0.03941, 0.18650, 0.65968, 1.68870, 2.93550, 2.62400  &
        ,1.66980, 1.53250, 1.35380, 1.21200, 1.13860/

        DATA(qscat(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00000, 0.00000, 0.00003, 0.00020, 0.00112, 0.00705  &
        ,0.04184, 0.19511, 0.68454, 1.72920, 2.95400, 2.56950  &
        ,1.67150, 1.52920, 1.33940, 1.20930, 1.13570/

        DATA(qscat(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00000, 0.00000, 0.00003, 0.00021, 0.00120, 0.00756  &
        ,0.04465, 0.20482, 0.71198, 1.77480, 2.97200, 2.51210  &
        ,1.67000, 1.52340, 1.32450, 1.20670, 1.13310/

        DATA(qscat(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00000, 0.00001, 0.00003, 0.00023, 0.00130, 0.00816  &
        ,0.04793, 0.21591, 0.74258, 1.82600, 2.99080, 2.45020  &
        ,1.67190, 1.51810, 1.31020, 1.20300, 1.13050/

        DATA(qscat(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00000, 0.00001, 0.00004, 0.00025, 0.00141, 0.00888  &
        ,0.05181, 0.22875, 0.77709, 1.88330, 3.01190, 2.37720  &
        ,1.67350, 1.50950, 1.29770, 1.19700, 1.12770/

        DATA(qscat(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00000, 0.00001, 0.00004, 0.00027, 0.00155, 0.00976  &
        ,0.05649, 0.24387, 0.81663, 1.94700, 3.03290, 2.29370  &
        ,1.67410, 1.49720, 1.28840, 1.18830, 1.12490/

        DATA(qscat(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00000, 0.00001, 0.00005, 0.00030, 0.00173, 0.01085  &
        ,0.06221, 0.26202, 0.86288, 2.01740, 3.04690, 2.20850  &
        ,1.66990, 1.47860, 1.28270, 1.17830, 1.12210/

        DATA(qscat(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00000, 0.00001, 0.00005, 0.00034, 0.00195, 0.01224  &
        ,0.06936, 0.28434, 0.91848, 2.09650, 3.04930, 2.10790  &
        ,1.66450, 1.45260, 1.27850, 1.17020, 1.11870/

        DATA(qscat(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00000, 0.00001, 0.00006, 0.00039, 0.00225, 0.01407  &
        ,0.07852, 0.31268, 0.98759, 2.18970, 3.04670, 2.00080  &
        ,1.64720, 1.41920, 1.27050, 1.16420, 1.11500/

        DATA(qscat(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00000, 0.00001, 0.00007, 0.00046, 0.00266, 0.01657  &
        ,0.09061, 0.35017, 1.07630, 2.30470, 3.02940, 1.88620  &
        ,1.61820, 1.38250, 1.25390, 1.15650, 1.11200/

        DATA(qscat(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00000, 0.00001, 0.00008, 0.00056, 0.00325, 0.02017  &
        ,0.10725, 0.40242, 1.19210, 2.44110, 2.97010, 1.77580  &
        ,1.57480, 1.35600, 1.23420, 1.14830, 1.10860/

        DATA(qscat(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00000, 0.00002, 0.00011, 0.00072, 0.00417, 0.02575  &
        ,0.13143, 0.47960, 1.34390, 2.59430, 2.86590, 1.67950  &
        ,1.52410, 1.35150, 1.22040, 1.13830, 1.10490/

        DATA(qscat(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00000, 0.00002, 0.00014, 0.00099, 0.00579, 0.03540  &
        ,0.16962, 0.59953, 1.55510, 2.78250, 2.65520, 1.63000  &
        ,1.48860, 1.33560, 1.19030, 1.12860, 1.10110/

        DATA(qscat(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00000, 0.00003, 0.00022, 0.00157, 0.00928, 0.05539  &
        ,0.23956, 0.79819, 1.89870, 2.96230, 2.26500, 1.63460  &
        ,1.46420, 1.26460, 1.17180, 1.11640, 1.09710/

        DATA(qscat(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00001, 0.00006, 0.00046, 0.00342, 0.02065, 0.11394  &
        ,0.43157, 1.25490, 2.49680, 2.88860, 1.70360, 1.52330  &
        ,1.33570, 1.21560, 1.13740, 1.10350, 1.09250/

! Qscat for Mineral Dust R(imaginary) = 0.0015 (LOW absorption)(No RH effect)
        DATA(qscat_dust(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00000,  0.00002,  0.00010,  0.00060,  0.00379,  0.02295  &
        , 0.11574,  0.43216,  1.18230,  2.24080,  2.56070,  2.52740  &
        , 2.33820,  1.90440,  1.81000,  1.61960,  1.45920/ 
 
        DATA(qscat_dust(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00010,  0.00063,  0.00402,  0.02560,  0.15020,  0.66682  &
        , 2.00350,  3.44090,  3.01910,  2.58140,  2.22970,  2.14250  &
        , 2.00430,  1.87080,  1.72780,  1.56920,  1.40860/ 
 
        DATA(qscat_dust(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00166,  0.01056,  0.06481,  0.32669,  1.23460,  2.91900  &
        , 3.74250,  2.56930,  2.34830,  2.13070,  2.05040,  1.93030  &
        , 1.80160,  1.65710,  1.49460,  1.33930,  1.21530/ 
 
        DATA(qscat_dust(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000  &
        , 0.00003,  0.00020,  0.00125,  0.00808,  0.05258,  0.31342  &
        , 1.17700,  1.80920,  1.38670,  1.24760,  1.24330/ 
 
        DATA(qscat_dust(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00010,  0.00066,  0.00420,  0.02663,  0.15649,  0.65345  &
        , 1.77580,  2.55670,  1.53740,  1.30830,  1.20600/ 
 
        DATA(qscat_dust(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00009  &
        , 0.00054,  0.00344,  0.02155,  0.12273,  0.47129,  1.15740  &
        , 1.71170,  1.59480,  1.27650,  1.14810,  1.13330/ 
 
        DATA(qscat_dust(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00015  &
        , 0.00094,  0.00603,  0.03906,  0.24079,  1.02410,  2.13950  &
        , 1.96750,  1.11300,  1.22940,  1.20200,  1.19660/ 
 
        DATA(qscat_dust(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00003,  0.00018  &
        , 0.00112,  0.00676,  0.03661,  0.14457,  0.46121,  1.13810  &
        , 2.05970,  2.21310,  1.51910,  1.29570,  1.16240/ 

! Qscat for Mineral Dust R(imaginary) = 0.003 (MID absorption)(No RH effect)
        DATA(qscat_dust(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00000,  0.00002,  0.00010,  0.00060,  0.00379,  0.02294  &
        , 0.11555,  0.43023,  1.17330,  2.21550,  2.52310,  2.47760  &
        , 2.27270,  1.81910,  1.69680,  1.48260,  1.31240/ 
 
        DATA(qscat_dust(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00010,  0.00063,  0.00402,  0.02559,  0.15004,  0.66465  &
        , 1.99000,  3.40630,  2.97380,  2.51370,  2.14500,  2.02750  &
        , 1.85490,  1.68360,  1.50960,  1.34620,  1.21980/ 
 
        DATA(qscat_dust(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00166,  0.01056,  0.06477,  0.32607,  1.22860,  2.89550  &
        , 3.69870,  2.51770,  2.27170,  2.03600,  1.92120,  1.76420  &
        , 1.59940,  1.43280,  1.28320,  1.18080,  1.13040/ 
 
        DATA(qscat_dust(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000  &
        , 0.00003,  0.00020,  0.00125,  0.00808,  0.05258,  0.31342  &
        , 1.17700,  1.80920,  1.38670,  1.24760,  1.24330/ 
 
        DATA(qscat_dust(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00010,  0.00066,  0.00420,  0.02663,  0.15649,  0.65345  &
        , 1.77580,  2.55670,  1.53740,  1.30830,  1.20600/ 
 
        DATA(qscat_dust(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00001,  0.00009  &
        , 0.00054,  0.00344,  0.02155,  0.12273,  0.47129,  1.15740  &
        , 1.71170,  1.59480,  1.27650,  1.14810,  1.13330/ 
 
        DATA(qscat_dust(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00002,  0.00015  &
        , 0.00094,  0.00603,  0.03906,  0.24079,  1.02410,  2.13950  &
        , 1.96750,  1.11300,  1.22940,  1.20200,  1.19660/ 
 
        DATA(qscat_dust(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00003,  0.00018  &
        , 0.00112,  0.00676,  0.03661,  0.14457,  0.46121,  1.13810  &
        , 2.05970,  2.21310,  1.51910,  1.29570,  1.16240/ 

! Qscat for Mineral Dust R(imaginary) = 0.008 (HIGH absorption)(No RH effect)
        DATA(qscat_dust(3,1,ib),ib=1,17)  & ! Band 1, All RH's
        /0.00000, 0.00002, 0.00010, 0.00060, 0.00379, 0.02291  &
        ,0.11493, 0.42401, 1.14410, 2.13480, 2.40680, 2.33500  &
        ,2.10380, 1.61510, 1.47100, 1.27120, 1.15620/

        DATA(qscat_dust(3,2,ib),ib=1,17)  & ! Band 2, All RH's
        /0.00001, 0.00063, 0.00402, 0.02557, 0.14955, 0.65760  &
        ,1.94660, 3.29530, 2.83420, 2.32070, 1.91820, 1.74640  &
        ,1.53580, 1.35690, 1.22550, 1.15300, 1.12440/

        DATA(qscat_dust(3,3,ib),ib=1,17)  & ! Band 3, All RH's
        /0.00177, 0.01133, 0.06890, 0.34280, 1.23190, 2.81180  &
        ,3.55540, 2.36150, 2.06270, 1.79890, 1.62930, 1.44460  &
        ,1.29170, 1.19010, 1.14040, 1.12370, 1.11720/

        DATA(qscat_dust(3,4,ib),ib=1,17)  & ! Band 4, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000  &
        ,0.00003, 0.00020, 0.00125, 0.00808, 0.05258, 0.31342  &
        ,1.17700, 1.80920, 1.38670, 1.24760, 1.24330/

        DATA(qscat_dust(3,5,ib),ib=1,17)  & ! Band 5, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00002  &
        ,0.00010, 0.00066, 0.00420, 0.02663, 0.15649, 0.65345  &
        ,1.77580, 2.55670, 1.53740, 1.30830, 1.20600/

        DATA(qscat_dust(3,6,ib),ib=1,17)  & ! Band 6, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00001, 0.00009  &
        ,0.00054, 0.00344, 0.02155, 0.12273, 0.47129, 1.15740  &
        ,1.71170, 1.59490, 1.27650, 1.14810, 1.13330/

        DATA(qscat_dust(3,7,ib),ib=1,17)  & ! Band 7, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00002, 0.00015  &
        ,0.00094, 0.00603, 0.03906, 0.24079, 1.02410, 2.13950  &
        ,1.96750, 1.11300, 1.22940, 1.20200, 1.19660/

        DATA(qscat_dust(3,8,ib),ib=1,17)  & ! Band 8, All RH's
        /0.00000, 0.00000, 0.00000, 0.00000, 0.00003, 0.00018  &
        ,0.00112, 0.00676, 0.03661, 0.14457, 0.46121, 1.13810  &
        ,2.05970, 2.21310, 1.51910, 1.29570, 1.16240/

! Qscat for Absorbing Carbon (1% BC, 99% OC)(No RH effect)
        DATA(qscat_carb(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00000,  0.00003,  0.00019,  0.00122,  0.00777,  0.04849  &
        , 0.25418,  0.98689,  2.31090,  3.12880,  2.17480,  1.87020  &
        , 1.52560,  1.36780,  1.23920,  1.16940,  1.13780/

        DATA(qscat_carb(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00010,  0.00065,  0.00414,  0.02635,  0.15382,  0.67180  &
        , 1.95580,  3.21550,  2.57880,  1.83720,  1.63100,  1.38500  &
        , 1.24630,  1.16880,  1.14060,  1.13040,  1.12380/

        DATA(qscat_carb(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00168,  0.01071,  0.06546,  0.32685,  1.20600,  2.76310  &
        , 3.41210,  2.19830,  1.88470,  1.60900,  1.43450,  1.28050  &
        , 1.18550,  1.14350,  1.12910,  1.12250,  1.11730/

        DATA(qscat_carb(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00004,  0.00022,  0.00143,  0.00932,  0.06155,  0.37469  &
        , 1.36790,  1.85340,  1.34360,  1.26400,  1.25110/

        DATA(qscat_carb(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00013,  0.00086,  0.00548,  0.03512,  0.20966,  0.89224  &
        , 2.01100,  2.09660,  1.22770,  1.23610,  1.19210/

        DATA(qscat_carb(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00007,  0.00043  &
        , 0.00279,  0.01836,  0.12544,  0.81406,  2.21360,  1.62300  &
        , 1.29270,  1.30940,  1.28940,  1.27790,  1.26630/

        DATA(qscat_carb(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00025  &
        , 0.00159,  0.01038,  0.07008,  0.46661,  1.93220,  2.41540  &
        , 1.55580,  1.30580,  1.26040,  1.25670,  1.24290/

        DATA(qscat_carb(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00058  &
        , 0.00370,  0.02321,  0.13279,  0.54355,  1.60440,  2.82840  &
        , 2.40960,  1.55340,  1.37700,  1.24730,  1.16110/

! Qscat for Absorbing Carbon (2% BC, 98% OC)(No RH effect)
        DATA(qscat_carb(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00000,  0.00003,  0.00019,  0.00124,  0.00790,  0.04929  &
        , 0.25710,  0.98291,  2.25740,  2.97020,  2.03230,  1.71230  &
        , 1.40210,  1.27390,  1.18670,  1.14780,  1.13170/

        DATA(qscat_carb(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00010,  0.00066,  0.00421,  0.02674,  0.15562,  0.67351  &
        , 1.92130,  3.08420,  2.41380,  1.69000,  1.49700,  1.29130  &
        , 1.19550,  1.15300,  1.13860,  1.13090,  1.12450/

        DATA(qscat_carb(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.00170,  0.01086,  0.06629,  0.32918,  1.19620,  2.68280  &
        , 3.22740,  2.01860,  1.70390,  1.44330,  1.29720,  1.19620  &
        , 1.15210,  1.13680,  1.12920,  1.12310,  1.11790/

        DATA(qscat_carb(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00001  &
        , 0.00004,  0.00022,  0.00144,  0.00935,  0.06173,  0.37502  &
        , 1.35900,  1.83760,  1.34160,  1.26330,  1.25140/

        DATA(qscat_carb(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00002  &
        , 0.00014,  0.00087,  0.00555,  0.03555,  0.21176,  0.89399  &
        , 1.98550,  2.04520,  1.21630,  1.22730,  1.19130/

        DATA(qscat_carb(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00007,  0.00043  &
        , 0.00278,  0.01835,  0.12531,  0.81059,  2.19370,  1.61690  &
        , 1.29100,  1.30830,  1.28930,  1.27780,  1.26630/

        DATA(qscat_carb(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00004,  0.00025  &
        , 0.00159,  0.01039,  0.07007,  0.46552,  1.91190,  2.38080  &
        , 1.53550,  1.29950,  1.25930,  1.25600,  1.24300/

        DATA(qscat_carb(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00000,  0.00000,  0.00000,  0.00001,  0.00009,  0.00060  &
        , 0.00378,  0.02372,  0.13536,  0.54986,  1.59030,  2.72850  &
        , 2.24840,  1.46690,  1.30380,  1.19790,  1.14560/

! Returning value for qscat:
        if (aerotype <= 2) value = qscat(aerotype,radband,bin,rh)
        if (aerotype == 3) value = qscat_dust(dust_ref_im,radband,bin)
        if (aerotype == 4) value = qscat_carb(1,radband,bin)
        if (aerotype == 5) value = qscat_carb(2,radband,bin)

return
END SUBROUTINE aeroqscat

!##############################################################################
Subroutine aerogasym (aerotype,radband,bin,rh,value)

use rrad3, only:dust_ref_im

implicit none

        INTEGER aerotype,radband,bin,rh,ib
        REAL gasym(2,8,17,20),value
        REAL gasym_dust(3,8,17)
        REAL gasym_carb(2,8,17)

! gasym Ammonium sulfate with 20% insoluble inclusion (dust low R(Im)))
        DATA(gasym(1,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        / 0.00187,  0.00483,  0.01231,  0.03118,  0.07419,  0.19005  &
        , 0.49242,  0.69369,  0.80106,  0.81543,  0.73993,  0.78459  &
        , 0.82767,  0.84881,  0.87068,  0.88537,  0.89958/

        DATA(gasym(1,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        / 0.00191,  0.00491,  0.01255,  0.03179,  0.07566,  0.19405  &
        , 0.50021,  0.69794,  0.80303,  0.81632,  0.74110,  0.79164  &
        , 0.83038,  0.85811,  0.87193,  0.88738,  0.90009/

        DATA(gasym(1,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        / 0.00194,  0.00501,  0.01280,  0.03245,  0.07727,  0.19842  &
        , 0.50840,  0.70244,  0.80505,  0.81717,  0.74209,  0.78602  &
        , 0.83210,  0.85212,  0.87425,  0.88525,  0.89913/

        DATA(gasym(1,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        / 0.00198,  0.00511,  0.01308,  0.03317,  0.07903,  0.20322  &
        , 0.51700,  0.70716,  0.80714,  0.81790,  0.74146,  0.79431  &
        , 0.83074,  0.85907,  0.87398,  0.88879,  0.89923/

        DATA(gasym(1,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        / 0.00202,  0.00522,  0.01338,  0.03397,  0.08098,  0.20854  &
        , 0.52604,  0.71212,  0.80928,  0.81843,  0.74213,  0.78781  &
        , 0.83395,  0.85726,  0.87352,  0.88743,  0.90090/

        DATA(gasym(1,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        / 0.00206,  0.00535,  0.01371,  0.03485,  0.08314,  0.21447  &
        , 0.53552,  0.71732,  0.81148,  0.81866,  0.74338,  0.79562  &
        , 0.83605,  0.85718,  0.87415,  0.88609,  0.90053/

        DATA(gasym(1,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        / 0.00211,  0.00549,  0.01409,  0.03584,  0.08556,  0.22113  &
        , 0.54543,  0.72275,  0.81376,  0.81861,  0.74621,  0.79596  &
        , 0.83681,  0.86229,  0.87770,  0.88822,  0.90184/

        DATA(gasym(1,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        / 0.00216,  0.00564,  0.01451,  0.03695,  0.08829,  0.22870  &
        , 0.55577,  0.72841,  0.81613,  0.81857,  0.74846,  0.79043  &
        , 0.83876,  0.86304,  0.87596,  0.89334,  0.90353/

        DATA(gasym(1,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        / 0.00222,  0.00581,  0.01499,  0.03822,  0.09142,  0.23739  &
        , 0.56651,  0.73430,  0.81865,  0.81881,  0.74804,  0.79537  &
        , 0.84040,  0.86375,  0.88013,  0.89083,  0.90383/

        DATA(gasym(1,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        / 0.00229,  0.00601,  0.01553,  0.03968,  0.09504,  0.24750  &
        , 0.57761,  0.74046,  0.82143,  0.81910,  0.75037,  0.79868  &
        , 0.84272,  0.86589,  0.88100,  0.89136,  0.90473/

        DATA(gasym(1,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        / 0.00237,  0.00624,  0.01617,  0.04139,  0.09930,  0.25946  &
        , 0.58898,  0.74692,  0.82457,  0.81874,  0.75533,  0.79687  &
        , 0.84428,  0.86628,  0.88316,  0.89364,  0.90554/

        DATA(gasym(1,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        / 0.00246,  0.00651,  0.01693,  0.04343,  0.10441,  0.27386  &
        , 0.60056,  0.75380,  0.82818,  0.81702,  0.75808,  0.80788  &
        , 0.84672,  0.86604,  0.88057,  0.89393,  0.90560/

        DATA(gasym(1,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        / 0.00256,  0.00683,  0.01784,  0.04592,  0.11068,  0.29160  &
        , 0.61231,  0.76133,  0.83213,  0.81535,  0.76358,  0.81005  &
        , 0.84903,  0.86821,  0.88206,  0.89720,  0.90796/

        DATA(gasym(1,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        / 0.00269,  0.00722,  0.01898,  0.04903,  0.11858,  0.31397  &
        , 0.62438,  0.76992,  0.83582,  0.81372,  0.76950,  0.81855  &
        , 0.85340,  0.86906,  0.88584,  0.89869,  0.90871/

        DATA(gasym(1,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        / 0.00285,  0.00772,  0.02043,  0.05304,  0.12889,  0.34302  &
        , 0.63746,  0.78001,  0.83887,  0.80796,  0.77631,  0.82251  &
        , 0.85564,  0.87383,  0.88810,  0.89999,  0.91028/

        DATA(gasym(1,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        / 0.00305,  0.00837,  0.02237,  0.05845,  0.14294,  0.38183  &
        , 0.65352,  0.79180,  0.84249,  0.80249,  0.78699,  0.83201  &
        , 0.85554,  0.87796,  0.88866,  0.90134,  0.91190/

        DATA(gasym(1,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        / 0.00331,  0.00926,  0.02508,  0.06613,  0.16322,  0.43486  &
        , 0.67682,  0.80557,  0.84607,  0.79253,  0.79745,  0.84326  &
        , 0.86806,  0.87771,  0.89282,  0.90406,  0.91402/

        DATA(gasym(1,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        / 0.00368,  0.01057,  0.02916,  0.07790,  0.19523,  0.50685  &
        , 0.71291,  0.82166,  0.84702,  0.77743,  0.80848,  0.84398  &
        , 0.87134,  0.88430,  0.89481,  0.90729,  0.91677/

        DATA(gasym(1,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        / 0.00425,  0.01270,  0.03615,  0.09886,  0.25403,  0.59123  &
        , 0.75669,  0.83788,  0.84058,  0.76828,  0.80618,  0.85339  &
        , 0.87174,  0.88806,  0.90126,  0.91020,  0.91994/

        DATA(gasym(1,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        / 0.00524,  0.01694,  0.05117,  0.14742,  0.39603,  0.66970  &
        , 0.80693,  0.85297,  0.80655,  0.80205,  0.83901,  0.86184  &
        , 0.88081,  0.89545,  0.90640,  0.91621,  0.92543/

        DATA(gasym(1,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        / 0.00849,  0.02182,  0.05540,  0.14088,  0.35082,  0.62082  &
        , 0.75525,  0.79558,  0.72584,  0.72424,  0.78566,  0.80938  &
        , 0.82638,  0.84025,  0.84942,  0.85955,  0.86636/

        DATA(gasym(1,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        / 0.00865,  0.02221,  0.05644,  0.14368,  0.35783,  0.62460  &
        , 0.75782,  0.79707,  0.72587,  0.72466,  0.78675,  0.80852  &
        , 0.82729,  0.84041,  0.84940,  0.85797,  0.86736/

        DATA(gasym(1,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        / 0.00880,  0.02263,  0.05757,  0.14673,  0.36543,  0.62870  &
        , 0.76050,  0.79859,  0.72596,  0.72674,  0.78941,  0.80984  &
        , 0.82870,  0.84190,  0.85239,  0.85962,  0.86798/

        DATA(gasym(1,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        / 0.00897,  0.02310,  0.05881,  0.15008,  0.37369,  0.63317  &
        , 0.76331,  0.80010,  0.72639,  0.72878,  0.79067,  0.81150  &
        , 0.82858,  0.84291,  0.85303,  0.86128,  0.87007/

        DATA(gasym(1,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        / 0.00915,  0.02360,  0.06017,  0.15377,  0.38270,  0.63806  &
        , 0.76626,  0.80154,  0.72703,  0.73108,  0.79324,  0.81200  &
        , 0.82935,  0.84393,  0.85270,  0.86187,  0.86860/

        DATA(gasym(1,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        / 0.00934,  0.02416,  0.06167,  0.15787,  0.39256,  0.64345  &
        , 0.76938,  0.80288,  0.72787,  0.73200,  0.79067,  0.81600  &
        , 0.83041,  0.84330,  0.85329,  0.86303,  0.87003/

        DATA(gasym(1,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        / 0.00956,  0.02477,  0.06334,  0.16246,  0.40341,  0.64941  &
        , 0.77270,  0.80406,  0.72862,  0.73340,  0.79717,  0.81867  &
        , 0.83272,  0.84541,  0.85604,  0.86424,  0.86982/

        DATA(gasym(1,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        / 0.00980,  0.02546,  0.06522,  0.16766,  0.41540,  0.65605  &
        , 0.77625,  0.80511,  0.72851,  0.73570,  0.79608,  0.81595  &
        , 0.83564,  0.84451,  0.85734,  0.86296,  0.87265/

        DATA(gasym(1,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        / 0.01007,  0.02624,  0.06736,  0.17359,  0.42872,  0.66346  &
        , 0.78006,  0.80620,  0.72829,  0.73829,  0.79820,  0.82004  &
        , 0.83407,  0.84771,  0.85632,  0.86401,  0.87278/

        DATA(gasym(1,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        / 0.01038,  0.02713,  0.06982,  0.18046,  0.44358,  0.67176  &
        , 0.78420,  0.80757,  0.72890,  0.74090,  0.79938,  0.82306  &
        , 0.83898,  0.84893,  0.85943,  0.86579,  0.87313/

        DATA(gasym(1,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        / 0.01073,  0.02816,  0.07268,  0.18853,  0.46023,  0.68102  &
        , 0.78870,  0.80917,  0.72965,  0.74294,  0.80450,  0.82427  &
        , 0.84152,  0.84969,  0.86084,  0.86607,  0.87257/

        DATA(gasym(1,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        / 0.01113,  0.02936,  0.07607,  0.19821,  0.47891,  0.69131  &
        , 0.79355,  0.81054,  0.72852,  0.74808,  0.80753,  0.82454  &
        , 0.84188,  0.85114,  0.86152,  0.86849,  0.87534/

        DATA(gasym(1,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        / 0.01161,  0.03080,  0.08019,  0.21006,  0.49982,  0.70260  &
        , 0.79878,  0.81106,  0.72915,  0.75060,  0.80884,  0.82811  &
        , 0.84229,  0.85312,  0.86219,  0.87016,  0.87514/

        DATA(gasym(1,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        / 0.01218,  0.03257,  0.08530,  0.22497,  0.52301,  0.71482  &
        , 0.80449,  0.81119,  0.72925,  0.75556,  0.81199,  0.83185  &
        , 0.84535,  0.85566,  0.86332,  0.87107,  0.87711/

        DATA(gasym(1,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        / 0.01288,  0.03480,  0.09185,  0.24436,  0.54827,  0.72788  &
        , 0.81102,  0.81094,  0.73168,  0.76304,  0.81673,  0.83221  &
        , 0.84814,  0.85631,  0.86369,  0.87211,  0.87797/

        DATA(gasym(1,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        / 0.01378,  0.03772,  0.10059,  0.27059,  0.57502,  0.74201  &
        , 0.81866,  0.80786,  0.73423,  0.77264,  0.81789,  0.83398  &
        , 0.84901,  0.85963,  0.86757,  0.87395,  0.87872/

        DATA(gasym(1,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        / 0.01497,  0.04172,  0.11288,  0.30789,  0.60297,  0.75827  &
        , 0.82573,  0.80249,  0.73979,  0.78614,  0.82149,  0.83856  &
        , 0.85197,  0.86249,  0.86956,  0.87525,  0.88034/

        DATA(gasym(1,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        / 0.01664,  0.04760,  0.13151,  0.36409,  0.63499,  0.77808  &
        , 0.83202,  0.79066,  0.74897,  0.80479,  0.82780,  0.84291  &
        , 0.85719,  0.86657,  0.87276,  0.87712,  0.88175/

        DATA(gasym(1,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        / 0.01919,  0.05715,  0.16381,  0.45478,  0.68339,  0.80251  &
        , 0.83523,  0.76912,  0.75969,  0.81392,  0.83549,  0.85114  &
        , 0.85945,  0.86818,  0.87435,  0.87947,  0.88366/

        DATA(gasym(1,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        / 0.02368,  0.07617,  0.23497,  0.58298,  0.75193,  0.83176  &
        , 0.82175,  0.74894,  0.78778,  0.82813,  0.84416,  0.85887  &
        , 0.86722,  0.87317,  0.87828,  0.88227,  0.88582/

        DATA(gasym(1,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        / 0.03881,  0.09905,  0.25222,  0.51963,  0.69486,  0.77757  &
        , 0.77036,  0.68247,  0.75387,  0.79048,  0.81340,  0.82732  &
        , 0.84203,  0.84902,  0.85872,  0.86540,  0.87545/

        DATA(gasym(1,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        / 0.03954,  0.10080,  0.25684,  0.52471,  0.69808,  0.77971  &
        , 0.77143,  0.68320,  0.75641,  0.79439,  0.81318,  0.82900  &
        , 0.84282,  0.85014,  0.85809,  0.86731,  0.87556/

        DATA(gasym(1,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        / 0.04023,  0.10270,  0.26184,  0.53015,  0.70147,  0.78192  &
        , 0.77250,  0.68353,  0.75898,  0.79231,  0.81625,  0.82987  &
        , 0.84326,  0.85108,  0.85928,  0.86731,  0.87688/

        DATA(gasym(1,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        / 0.04098,  0.10477,  0.26729,  0.53598,  0.70505,  0.78421  &
        , 0.77357,  0.68490,  0.76188,  0.79408,  0.81866,  0.83220  &
        , 0.84343,  0.85129,  0.85919,  0.86824,  0.87659/

        DATA(gasym(1,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        / 0.04179,  0.10703,  0.27326,  0.54227,  0.70883,  0.78658  &
        , 0.77461,  0.68599,  0.76466,  0.79497,  0.81743,  0.83264  &
        , 0.84141,  0.85297,  0.86043,  0.86691,  0.87680/

        DATA(gasym(1,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        / 0.04268,  0.10953,  0.27981,  0.54907,  0.71283,  0.78902  &
        , 0.77548,  0.68739,  0.76668,  0.79851,  0.81886,  0.83545  &
        , 0.84651,  0.85485,  0.86156,  0.86968,  0.87685/

        DATA(gasym(1,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        / 0.04366,  0.11230,  0.28706,  0.55647,  0.71707,  0.79152  &
        , 0.77607,  0.68889,  0.77045,  0.79829,  0.82115,  0.83719  &
        , 0.84487,  0.85698,  0.86342,  0.86953,  0.87703/

        DATA(gasym(1,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        / 0.04475,  0.11540,  0.29514,  0.56455,  0.72159,  0.79407  &
        , 0.77659,  0.69039,  0.77262,  0.80182,  0.82109,  0.83560  &
        , 0.84850,  0.85714,  0.86362,  0.86947,  0.87844/

        DATA(gasym(1,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        / 0.04597,  0.11889,  0.30422,  0.57345,  0.72641,  0.79666  &
        , 0.77698,  0.69232,  0.77644,  0.80448,  0.82417,  0.84092  &
        , 0.84993,  0.85512,  0.86450,  0.87054,  0.87885/

        DATA(gasym(1,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        / 0.04734,  0.12287,  0.31450,  0.58328,  0.73158,  0.79929  &
        , 0.77736,  0.69335,  0.77864,  0.80592,  0.82731,  0.83964  &
        , 0.85028,  0.85616,  0.86505,  0.87162,  0.87874/

        DATA(gasym(1,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        / 0.04891,  0.12748,  0.32628,  0.59421,  0.73716,  0.80199  &
        , 0.77763,  0.69483,  0.78237,  0.80978,  0.82774,  0.84082  &
        , 0.85012,  0.85883,  0.86554,  0.87233,  0.87870/

        DATA(gasym(1,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        / 0.05072,  0.13290,  0.33993,  0.60639,  0.74324,  0.80486  &
        , 0.77774,  0.69703,  0.78937,  0.81249,  0.82981,  0.84308  &
        , 0.85095,  0.86045,  0.86652,  0.87300,  0.87904/

        DATA(gasym(1,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        / 0.05286,  0.13939,  0.35593,  0.62001,  0.74994,  0.80799  &
        , 0.77677,  0.70050,  0.79198,  0.81541,  0.83306,  0.84457  &
        , 0.85406,  0.86236,  0.86803,  0.87402,  0.88075/

        DATA(gasym(1,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        / 0.05543,  0.14734,  0.37493,  0.63518,  0.75744,  0.81134  &
        , 0.77535,  0.70488,  0.79739,  0.82165,  0.83267,  0.84650  &
        , 0.85453,  0.86263,  0.86872,  0.87410,  0.88064/

        DATA(gasym(1,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        / 0.05860,  0.15738,  0.39778,  0.65201,  0.76613,  0.81475  &
        , 0.77321,  0.71066,  0.80256,  0.82315,  0.83665,  0.84717  &
        , 0.85729,  0.86483,  0.87046,  0.87587,  0.88104/

        DATA(gasym(1,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        / 0.06263,  0.17052,  0.42555,  0.67066,  0.77667,  0.81810  &
        , 0.76780,  0.72151,  0.80737,  0.82529,  0.84103,  0.85183  &
        , 0.85790,  0.86570,  0.87231,  0.87641,  0.88180/

        DATA(gasym(1,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        / 0.06797,  0.18855,  0.45967,  0.69167,  0.78979,  0.82079  &
        , 0.76019,  0.73768,  0.81182,  0.82851,  0.84377,  0.85444  &
        , 0.86211,  0.86951,  0.87295,  0.87752,  0.88217/

        DATA(gasym(1,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        / 0.07548,  0.21493,  0.50250,  0.71662,  0.80504,  0.82069  &
        , 0.74812,  0.76111,  0.81510,  0.83221,  0.84585,  0.85836  &
        , 0.86614,  0.87098,  0.87500,  0.87880,  0.88283/

        DATA(gasym(1,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        / 0.08693,  0.25742,  0.56095,  0.74759,  0.82028,  0.81432  &
        , 0.73036,  0.78855,  0.81905,  0.83971,  0.85382,  0.86023  &
        , 0.86912,  0.87329,  0.87666,  0.88010,  0.88327/

        DATA(gasym(1,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        / 0.10706,  0.33768,  0.64779,  0.78698,  0.83186,  0.78599  &
        , 0.73641,  0.81766,  0.83410,  0.84890,  0.85956,  0.86749  &
        , 0.87155,  0.87598,  0.87874,  0.88128,  0.88351/

        DATA(gasym(1,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        / 0.00001,  0.00003,  0.00008,  0.00020,  0.00048,  0.00121  &
        , 0.00305,  0.00769,  0.01933,  0.04873,  0.12603,  0.34952  &
        , 0.62871,  0.78210,  0.85086,  0.89593,  0.91627/

        DATA(gasym(1,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        / 0.00001,  0.00003,  0.00008,  0.00020,  0.00049,  0.00123  &
        , 0.00311,  0.00784,  0.01972,  0.04972,  0.12880,  0.35683  &
        , 0.63342,  0.78454,  0.85248,  0.89674,  0.91654/

        DATA(gasym(1,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        / 0.00001,  0.00003,  0.00008,  0.00021,  0.00050,  0.00126  &
        , 0.00318,  0.00801,  0.02015,  0.05081,  0.13183,  0.36470  &
        , 0.63844,  0.78710,  0.85415,  0.89756,  0.91682/

        DATA(gasym(1,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        / 0.00001,  0.00003,  0.00008,  0.00021,  0.00051,  0.00129  &
        , 0.00325,  0.00819,  0.02061,  0.05201,  0.13518,  0.37319  &
        , 0.64380,  0.78979,  0.85588,  0.89841,  0.91711/

        DATA(gasym(1,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        / 0.00001,  0.00003,  0.00009,  0.00022,  0.00052,  0.00132  &
        , 0.00333,  0.00840,  0.02113,  0.05333,  0.13891,  0.38239  &
        , 0.64953,  0.79262,  0.85768,  0.89927,  0.91741/

        DATA(gasym(1,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        / 0.00001,  0.00003,  0.00009,  0.00022,  0.00053,  0.00135  &
        , 0.00342,  0.00862,  0.02170,  0.05481,  0.14308,  0.39238  &
        , 0.65568,  0.79562,  0.85956,  0.90016,  0.91773/

        DATA(gasym(1,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        / 0.00001,  0.00003,  0.00009,  0.00023,  0.00055,  0.00139  &
        , 0.00352,  0.00888,  0.02235,  0.05647,  0.14779,  0.40328  &
        , 0.66229,  0.79883,  0.86153,  0.90108,  0.91807/

        DATA(gasym(1,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        / 0.00001,  0.00004,  0.00009,  0.00024,  0.00057,  0.00144  &
        , 0.00364,  0.00917,  0.02308,  0.05836,  0.15317,  0.41523  &
        , 0.66942,  0.80225,  0.86361,  0.90203,  0.91842/

        DATA(gasym(1,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        / 0.00001,  0.00004,  0.00010,  0.00024,  0.00059,  0.00149  &
        , 0.00377,  0.00950,  0.02392,  0.06053,  0.15939,  0.42837  &
        , 0.67712,  0.80594,  0.86583,  0.90303,  0.91880/

        DATA(gasym(1,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        / 0.00001,  0.00004,  0.00010,  0.00025,  0.00061,  0.00155  &
        , 0.00392,  0.00988,  0.02490,  0.06307,  0.16668,  0.44291  &
        , 0.68548,  0.80991,  0.86822,  0.90408,  0.91921/

        DATA(gasym(1,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        / 0.00002,  0.00004,  0.00010,  0.00026,  0.00064,  0.00162  &
        , 0.00410,  0.01034,  0.02605,  0.06607,  0.17538,  0.45908  &
        , 0.69459,  0.81420,  0.87082,  0.90520,  0.91965/

        DATA(gasym(1,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        / 0.00002,  0.00004,  0.00011,  0.00028,  0.00067,  0.00170  &
        , 0.00431,  0.01089,  0.02744,  0.06971,  0.18597,  0.47715  &
        , 0.70458,  0.81885,  0.87364,  0.90639,  0.92014/

        DATA(gasym(1,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        / 0.00002,  0.00004,  0.00011,  0.00029,  0.00071,  0.00181  &
        , 0.00458,  0.01156,  0.02916,  0.07421,  0.19915,  0.49744  &
        , 0.71560,  0.82395,  0.87673,  0.90768,  0.92067/

        DATA(gasym(1,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        / 0.00002,  0.00005,  0.00012,  0.00031,  0.00076,  0.00194  &
        , 0.00491,  0.01242,  0.03133,  0.07995,  0.21605,  0.52035  &
        , 0.72787,  0.82965,  0.88013,  0.90908,  0.92127/

        DATA(gasym(1,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        / 0.00002,  0.00005,  0.00013,  0.00034,  0.00082,  0.00211  &
        , 0.00535,  0.01353,  0.03419,  0.08754,  0.23840,  0.54636  &
        , 0.74170,  0.83614,  0.88395,  0.91064,  0.92194/

        DATA(gasym(1,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        / 0.00002,  0.00005,  0.00014,  0.00037,  0.00091,  0.00234  &
        , 0.00595,  0.01506,  0.03810,  0.09805,  0.26912,  0.57618  &
        , 0.75738,  0.84362,  0.88828,  0.91238,  0.92272/

        DATA(gasym(1,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        / 0.00002,  0.00006,  0.00016,  0.00042,  0.00104,  0.00267  &
        , 0.00681,  0.01726,  0.04377,  0.11356,  0.31303,  0.61107  &
        , 0.77520,  0.85223,  0.89326,  0.91435,  0.92363/

        DATA(gasym(1,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        / 0.00002,  0.00007,  0.00019,  0.00050,  0.00123,  0.00319  &
        , 0.00817,  0.02077,  0.05286,  0.13895,  0.37883,  0.65380  &
        , 0.79606,  0.86259,  0.89905,  0.91667,  0.92474/

        DATA(gasym(1,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        / 0.00003,  0.00008,  0.00023,  0.00063,  0.00158,  0.00414  &
        , 0.01065,  0.02719,  0.06976,  0.18772,  0.47735,  0.70782  &
        , 0.82144,  0.87565,  0.90585,  0.91948,  0.92613/

        DATA(gasym(1,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        / 0.00003,  0.00011,  0.00033,  0.00094,  0.00243,  0.00647  &
        , 0.01690,  0.04366,  0.11487,  0.31858,  0.61698,  0.77836  &
        , 0.85425,  0.89391,  0.91422,  0.92331,  0.92807/

        DATA(gasym(1,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        / 0.00003,  0.00007,  0.00018,  0.00045,  0.00109,  0.00275  &
        , 0.00695,  0.01750,  0.04418,  0.11414,  0.31696,  0.64379  &
        , 0.80159,  0.87297,  0.90762,  0.92878,  0.93681/

        DATA(gasym(1,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        / 0.00003,  0.00007,  0.00018,  0.00046,  0.00111,  0.00280  &
        , 0.00707,  0.01782,  0.04501,  0.11645,  0.32349,  0.64852  &
        , 0.80445,  0.87478,  0.90887,  0.92928,  0.93721/

        DATA(gasym(1,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        / 0.00003,  0.00007,  0.00018,  0.00047,  0.00113,  0.00286  &
        , 0.00721,  0.01818,  0.04593,  0.11898,  0.33062,  0.65353  &
        , 0.80741,  0.87662,  0.91012,  0.92980,  0.93762/

        DATA(gasym(1,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        / 0.00003,  0.00007,  0.00019,  0.00048,  0.00115,  0.00292  &
        , 0.00737,  0.01857,  0.04693,  0.12178,  0.33841,  0.65883  &
        , 0.81048,  0.87851,  0.91137,  0.93032,  0.93803/

        DATA(gasym(1,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        / 0.00003,  0.00008,  0.00019,  0.00049,  0.00118,  0.00299  &
        , 0.00754,  0.01900,  0.04805,  0.12490,  0.34699,  0.66448  &
        , 0.81368,  0.88045,  0.91263,  0.93087,  0.93844/

        DATA(gasym(1,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        / 0.00003,  0.00008,  0.00020,  0.00050,  0.00121,  0.00306  &
        , 0.00773,  0.01948,  0.04930,  0.12839,  0.35647,  0.67049  &
        , 0.81703,  0.88244,  0.91390,  0.93143,  0.93887/

        DATA(gasym(1,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        / 0.00003,  0.00008,  0.00020,  0.00052,  0.00124,  0.00314  &
        , 0.00794,  0.02002,  0.05071,  0.13234,  0.36701,  0.67693  &
        , 0.82056,  0.88449,  0.91519,  0.93200,  0.93929/

        DATA(gasym(1,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        / 0.00003,  0.00008,  0.00021,  0.00053,  0.00128,  0.00324  &
        , 0.00818,  0.02064,  0.05231,  0.13685,  0.37884,  0.68387  &
        , 0.82429,  0.88661,  0.91650,  0.93260,  0.93973/

        DATA(gasym(1,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        / 0.00003,  0.00008,  0.00021,  0.00055,  0.00132,  0.00335  &
        , 0.00846,  0.02135,  0.05415,  0.14207,  0.39220,  0.69138  &
        , 0.82827,  0.88879,  0.91782,  0.93322,  0.94018/

        DATA(gasym(1,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        / 0.00003,  0.00009,  0.00022,  0.00057,  0.00137,  0.00347  &
        , 0.00878,  0.02217,  0.05631,  0.14820,  0.40742,  0.69958  &
        , 0.83254,  0.89107,  0.91917,  0.93386,  0.94063/

        DATA(gasym(1,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        / 0.00003,  0.00009,  0.00023,  0.00059,  0.00143,  0.00362  &
        , 0.00916,  0.02315,  0.05887,  0.15552,  0.42493,  0.70859  &
        , 0.83717,  0.89348,  0.92057,  0.93454,  0.94109/

        DATA(gasym(1,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        / 0.00003,  0.00009,  0.00024,  0.00062,  0.00150,  0.00380  &
        , 0.00962,  0.02433,  0.06197,  0.16444,  0.44527,  0.71860  &
        , 0.84221,  0.89606,  0.92201,  0.93524,  0.94156/

        DATA(gasym(1,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        / 0.00004,  0.00010,  0.00025,  0.00065,  0.00158,  0.00402  &
        , 0.01019,  0.02579,  0.06583,  0.17557,  0.46908,  0.72983  &
        , 0.84773,  0.89885,  0.92353,  0.93599,  0.94205/

        DATA(gasym(1,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        / 0.00004,  0.00010,  0.00027,  0.00070,  0.00169,  0.00430  &
        , 0.01091,  0.02765,  0.07075,  0.18989,  0.49717,  0.74259  &
        , 0.85380,  0.90186,  0.92513,  0.93678,  0.94255/

        DATA(gasym(1,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        / 0.00004,  0.00011,  0.00029,  0.00075,  0.00183,  0.00467  &
        , 0.01186,  0.03008,  0.07726,  0.20894,  0.53037,  0.75726  &
        , 0.86051,  0.90514,  0.92683,  0.93762,  0.94306/

        DATA(gasym(1,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        / 0.00004,  0.00012,  0.00031,  0.00083,  0.00202,  0.00516  &
        , 0.01315,  0.03343,  0.08628,  0.23541,  0.56936,  0.77429  &
        , 0.86804,  0.90876,  0.92868,  0.93854,  0.94359/

        DATA(gasym(1,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        / 0.00005,  0.00013,  0.00035,  0.00093,  0.00229,  0.00588  &
        , 0.01501,  0.03830,  0.09958,  0.27418,  0.61420,  0.79409  &
        , 0.87665,  0.91286,  0.93070,  0.93954,  0.94414/

        DATA(gasym(1,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        / 0.00005,  0.00015,  0.00041,  0.00110,  0.00271,  0.00701  &
        , 0.01798,  0.04609,  0.12129,  0.33555,  0.66459,  0.81720  &
        , 0.88663,  0.91761,  0.93300,  0.94068,  0.94473/

        DATA(gasym(1,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        / 0.00006,  0.00018,  0.00051,  0.00139,  0.00347,  0.00907  &
        , 0.02342,  0.06059,  0.16280,  0.43979,  0.72114,  0.84414  &
        , 0.89835,  0.92328,  0.93573,  0.94204,  0.94541/

        DATA(gasym(1,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        / 0.00007,  0.00024,  0.00071,  0.00205,  0.00531,  0.01417  &
        , 0.03724,  0.09892,  0.27594,  0.61726,  0.79653,  0.87787  &
        , 0.91348,  0.93070,  0.93934,  0.94387,  0.94630/

        DATA(gasym(1,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        / 0.00006,  0.00017,  0.00043,  0.00110,  0.00260,  0.00658  &
        , 0.01670,  0.04278,  0.11378,  0.32215,  0.67624,  0.82827  &
        , 0.90221,  0.93561,  0.94966,  0.95759,  0.96206/

        DATA(gasym(1,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        / 0.00007,  0.00017,  0.00044,  0.00113,  0.00266,  0.00675  &
        , 0.01713,  0.04389,  0.11673,  0.32994,  0.68261,  0.83173  &
        , 0.90423,  0.93696,  0.95067,  0.95865,  0.96307/

        DATA(gasym(1,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        / 0.00007,  0.00018,  0.00045,  0.00116,  0.00274,  0.00694  &
        , 0.01761,  0.04511,  0.11993,  0.33842,  0.68908,  0.83529  &
        , 0.90628,  0.93830,  0.95167,  0.95971,  0.96408/

        DATA(gasym(1,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        / 0.00007,  0.00018,  0.00047,  0.00119,  0.00282,  0.00714  &
        , 0.01813,  0.04643,  0.12345,  0.34768,  0.69564,  0.83893  &
        , 0.90839,  0.93962,  0.95265,  0.96078,  0.96508/

        DATA(gasym(1,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        / 0.00007,  0.00019,  0.00048,  0.00123,  0.00290,  0.00737  &
        , 0.01870,  0.04789,  0.12732,  0.35785,  0.70227,  0.84266  &
        , 0.91054,  0.94090,  0.95360,  0.96186,  0.96609/

        DATA(gasym(1,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        / 0.00007,  0.00019,  0.00049,  0.00127,  0.00300,  0.00761  &
        , 0.01933,  0.04952,  0.13162,  0.36911,  0.70895,  0.84648  &
        , 0.91276,  0.94215,  0.95453,  0.96294,  0.96709/

        DATA(gasym(1,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        / 0.00008,  0.00020,  0.00051,  0.00131,  0.00311,  0.00789  &
        , 0.02003,  0.05133,  0.13644,  0.38165,  0.71565,  0.85039  &
        , 0.91504,  0.94335,  0.95543,  0.96402,  0.96808/

        DATA(gasym(1,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        / 0.00008,  0.00020,  0.00053,  0.00136,  0.00323,  0.00820  &
        , 0.02083,  0.05337,  0.14189,  0.39574,  0.72236,  0.85440  &
        , 0.91735,  0.94451,  0.95631,  0.96510,  0.96906/

        DATA(gasym(1,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        / 0.00008,  0.00021,  0.00055,  0.00141,  0.00337,  0.00856  &
        , 0.02174,  0.05571,  0.14813,  0.41171,  0.72910,  0.85852  &
        , 0.91966,  0.94561,  0.95716,  0.96618,  0.97002/

        DATA(gasym(1,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        / 0.00008,  0.00022,  0.00057,  0.00148,  0.00352,  0.00897  &
        , 0.02279,  0.05841,  0.15537,  0.43000,  0.73592,  0.86277  &
        , 0.92191,  0.94660,  0.95797,  0.96726,  0.97097/

        DATA(gasym(1,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        / 0.00009,  0.00023,  0.00060,  0.00155,  0.00371,  0.00944  &
        , 0.02402,  0.06158,  0.16391,  0.45117,  0.74300,  0.86719  &
        , 0.92407,  0.94745,  0.95878,  0.96832,  0.97189/

        DATA(gasym(1,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        / 0.00009,  0.00024,  0.00063,  0.00164,  0.00393,  0.01002  &
        , 0.02549,  0.06538,  0.17418,  0.47592,  0.75065,  0.87184  &
        , 0.92618,  0.94812,  0.95957,  0.96937,  0.97277/

        DATA(gasym(1,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        / 0.00010,  0.00026,  0.00067,  0.00175,  0.00420,  0.01072  &
        , 0.02728,  0.07002,  0.18685,  0.50512,  0.75948,  0.87684  &
        , 0.92834,  0.94850,  0.96041,  0.97039,  0.97361/

        DATA(gasym(1,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        / 0.00010,  0.00027,  0.00072,  0.00188,  0.00453,  0.01159  &
        , 0.02954,  0.07588,  0.20292,  0.53968,  0.77047,  0.88248  &
        , 0.93051,  0.94854,  0.96134,  0.97137,  0.97439/

        DATA(gasym(1,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        / 0.00011,  0.00029,  0.00078,  0.00205,  0.00497,  0.01272  &
        , 0.03246,  0.08350,  0.22405,  0.58020,  0.78491,  0.88915  &
        , 0.93245,  0.94800,  0.96247,  0.97231,  0.97510/

        DATA(gasym(1,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        / 0.00012,  0.00032,  0.00087,  0.00228,  0.00555,  0.01425  &
        , 0.03642,  0.09386,  0.25313,  0.62591,  0.80357,  0.89656  &
        , 0.93401,  0.94670,  0.96398,  0.97321,  0.97571/

        DATA(gasym(1,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        / 0.00013,  0.00036,  0.00098,  0.00260,  0.00637,  0.01642  &
        , 0.04208,  0.10882,  0.29563,  0.67227,  0.82468,  0.90250  &
        , 0.93453,  0.94432,  0.96594,  0.97402,  0.97621/

        DATA(gasym(1,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        / 0.00014,  0.00041,  0.00115,  0.00309,  0.00764,  0.01981  &
        , 0.05097,  0.13264,  0.36369,  0.70969,  0.84410,  0.90976  &
        , 0.93341,  0.94151,  0.96817,  0.97458,  0.97659/

        DATA(gasym(1,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        / 0.00017,  0.00050,  0.00144,  0.00395,  0.00991,  0.02592  &
        , 0.06713,  0.17699,  0.48390,  0.73971,  0.86477,  0.91643  &
        , 0.92767,  0.94373,  0.97011,  0.97499,  0.97682/

        DATA(gasym(1,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        / 0.00021,  0.00067,  0.00206,  0.00592,  0.01532,  0.04092  &
        , 0.10808,  0.29582,  0.67548,  0.82072,  0.89362,  0.91938  &
        , 0.91039,  0.95937,  0.97268,  0.97574,  0.97700/

        DATA(gasym(1,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        / 0.00007,  0.00018,  0.00045,  0.00114,  0.00274,  0.00693  &
        , 0.01746,  0.04414,  0.11510,  0.32844,  0.62906,  0.76296  &
        , 0.81231,  0.86433,  0.90421,  0.91371,  0.91848/

        DATA(gasym(1,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        / 0.00007,  0.00018,  0.00046,  0.00116,  0.00279,  0.00705  &
        , 0.01777,  0.04494,  0.11724,  0.33453,  0.63375,  0.76682  &
        , 0.81585,  0.86606,  0.90653,  0.91583,  0.92055/

        DATA(gasym(1,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        / 0.00007,  0.00018,  0.00046,  0.00118,  0.00284,  0.00718  &
        , 0.01811,  0.04581,  0.11958,  0.34116,  0.63864,  0.77088  &
        , 0.81957,  0.86793,  0.90892,  0.91803,  0.92268/

        DATA(gasym(1,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        / 0.00007,  0.00019,  0.00047,  0.00121,  0.00290,  0.00733  &
        , 0.01848,  0.04677,  0.12214,  0.34840,  0.64376,  0.77514  &
        , 0.82350,  0.86996,  0.91137,  0.92031,  0.92489/

        DATA(gasym(1,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        / 0.00007,  0.00019,  0.00048,  0.00123,  0.00296,  0.00749  &
        , 0.01890,  0.04782,  0.12497,  0.35637,  0.64914,  0.77964  &
        , 0.82765,  0.87218,  0.91389,  0.92267,  0.92718/

        DATA(gasym(1,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        / 0.00007,  0.00019,  0.00049,  0.00126,  0.00303,  0.00767  &
        , 0.01935,  0.04899,  0.12812,  0.36518,  0.65483,  0.78441  &
        , 0.83207,  0.87459,  0.91649,  0.92513,  0.92955/

        DATA(gasym(1,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        / 0.00008,  0.00020,  0.00051,  0.00129,  0.00311,  0.00787  &
        , 0.01986,  0.05029,  0.13166,  0.37501,  0.66087,  0.78949  &
        , 0.83678,  0.87717,  0.91916,  0.92770,  0.93203/

        DATA(gasym(1,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        / 0.00008,  0.00020,  0.00052,  0.00133,  0.00320,  0.00810  &
        , 0.02044,  0.05178,  0.13568,  0.38605,  0.66734,  0.79493  &
        , 0.84180,  0.87991,  0.92193,  0.93037,  0.93461/

        DATA(gasym(1,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        / 0.00008,  0.00021,  0.00053,  0.00137,  0.00330,  0.00835  &
        , 0.02110,  0.05347,  0.14029,  0.39859,  0.67434,  0.80077  &
        , 0.84714,  0.88277,  0.92482,  0.93318,  0.93732/

        DATA(gasym(1,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        / 0.00008,  0.00021,  0.00055,  0.00141,  0.00342,  0.00865  &
        , 0.02186,  0.05544,  0.14567,  0.41297,  0.68201,  0.80709  &
        , 0.85281,  0.88582,  0.92785,  0.93612,  0.94015/

        DATA(gasym(1,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        / 0.00008,  0.00022,  0.00057,  0.00147,  0.00355,  0.00900  &
        , 0.02276,  0.05777,  0.15205,  0.42967,  0.69056,  0.81398  &
        , 0.85875,  0.88922,  0.93100,  0.93922,  0.94313/

        DATA(gasym(1,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        / 0.00009,  0.00023,  0.00060,  0.00153,  0.00371,  0.00943  &
        , 0.02385,  0.06057,  0.15976,  0.44931,  0.70028,  0.82151  &
        , 0.86491,  0.89322,  0.93428,  0.94250,  0.94628/

        DATA(gasym(1,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        / 0.00009,  0.00024,  0.00063,  0.00161,  0.00391,  0.00994  &
        , 0.02518,  0.06402,  0.16934,  0.47270,  0.71163,  0.82980  &
        , 0.87129,  0.89780,  0.93778,  0.94596,  0.94960/

        DATA(gasym(1,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        / 0.00009,  0.00025,  0.00066,  0.00171,  0.00416,  0.01060  &
        , 0.02686,  0.06841,  0.18157,  0.50087,  0.72525,  0.83900  &
        , 0.87808,  0.90255,  0.94150,  0.94964,  0.95311/

        DATA(gasym(1,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        / 0.00010,  0.00027,  0.00071,  0.00184,  0.00449,  0.01145  &
        , 0.02907,  0.07417,  0.19777,  0.53500,  0.74206,  0.84931  &
        , 0.88568,  0.90770,  0.94550,  0.95353,  0.95680/

        DATA(gasym(1,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        / 0.00011,  0.00029,  0.00077,  0.00202,  0.00493,  0.01261  &
        , 0.03209,  0.08209,  0.22027,  0.57604,  0.76305,  0.86103  &
        , 0.89432,  0.91391,  0.94980,  0.95763,  0.96066/

        DATA(gasym(1,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        / 0.00011,  0.00032,  0.00086,  0.00227,  0.00557,  0.01429  &
        , 0.03647,  0.09366,  0.25344,  0.62343,  0.78850,  0.87411  &
        , 0.90314,  0.92061,  0.95449,  0.96191,  0.96465/

        DATA(gasym(1,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        / 0.00013,  0.00036,  0.00099,  0.00265,  0.00656,  0.01695  &
        , 0.04344,  0.11229,  0.30736,  0.67273,  0.81658,  0.88776  &
        , 0.91221,  0.92840,  0.95966,  0.96635,  0.96875/

        DATA(gasym(1,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        / 0.00014,  0.00043,  0.00122,  0.00333,  0.00836,  0.02179  &
        , 0.05627,  0.14726,  0.40678,  0.71493,  0.84530,  0.90448  &
        , 0.92079,  0.93815,  0.96539,  0.97084,  0.97284/

        DATA(gasym(1,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        / 0.00018,  0.00057,  0.00171,  0.00491,  0.01270,  0.03387  &
        , 0.08916,  0.24118,  0.61217,  0.78988,  0.88617,  0.92219  &
        , 0.92601,  0.95470,  0.97202,  0.97544,  0.97686/

        DATA(gasym(1,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        / 0.00021,  0.00054,  0.00137,  0.00350,  0.00834,  0.02104  &
        , 0.05293,  0.13426,  0.35861,  0.65635,  0.79545,  0.84985  &
        , 0.82088,  0.84615,  0.90728,  0.93448,  0.95451/

        DATA(gasym(1,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        / 0.00021,  0.00055,  0.00140,  0.00357,  0.00851,  0.02147  &
        , 0.05403,  0.13712,  0.36650,  0.65965,  0.79787,  0.85059  &
        , 0.82071,  0.84788,  0.90841,  0.93533,  0.95503/

        DATA(gasym(1,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        / 0.00022,  0.00056,  0.00143,  0.00364,  0.00869,  0.02194  &
        , 0.05522,  0.14026,  0.37507,  0.66321,  0.80041,  0.85142  &
        , 0.82073,  0.84970,  0.90943,  0.93619,  0.95557/

        DATA(gasym(1,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        / 0.00022,  0.00057,  0.00146,  0.00372,  0.00890,  0.02246  &
        , 0.05653,  0.14371,  0.38442,  0.66710,  0.80309,  0.85237  &
        , 0.82079,  0.85193,  0.91050,  0.93711,  0.95613/

        DATA(gasym(1,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        / 0.00022,  0.00058,  0.00150,  0.00382,  0.00912,  0.02303  &
        , 0.05798,  0.14752,  0.39465,  0.67139,  0.80590,  0.85347  &
        , 0.82061,  0.85412,  0.91169,  0.93808,  0.95672/

        DATA(gasym(1,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        / 0.00023,  0.00060,  0.00153,  0.00392,  0.00937,  0.02366  &
        , 0.05959,  0.15178,  0.40591,  0.67618,  0.80885,  0.85470  &
        , 0.81987,  0.85634,  0.91279,  0.93910,  0.95735/

        DATA(gasym(1,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        / 0.00024,  0.00061,  0.00158,  0.00403,  0.00964,  0.02437  &
        , 0.06139,  0.15656,  0.41836,  0.68156,  0.81199,  0.85599  &
        , 0.81859,  0.85869,  0.91388,  0.94016,  0.95801/

        DATA(gasym(1,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        / 0.00024,  0.00063,  0.00163,  0.00416,  0.00995,  0.02517  &
        , 0.06343,  0.16199,  0.43220,  0.68768,  0.81532,  0.85720  &
        , 0.81754,  0.86155,  0.91507,  0.94127,  0.95870/

        DATA(gasym(1,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        / 0.00025,  0.00065,  0.00168,  0.00430,  0.01031,  0.02608  &
        , 0.06576,  0.16823,  0.44766,  0.69470,  0.81891,  0.85820  &
        , 0.81712,  0.86420,  0.91619,  0.94243,  0.95943/

        DATA(gasym(1,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        / 0.00026,  0.00067,  0.00174,  0.00447,  0.01072,  0.02714  &
        , 0.06847,  0.17551,  0.46503,  0.70280,  0.82280,  0.85889  &
        , 0.81654,  0.86710,  0.91744,  0.94362,  0.96019/

        DATA(gasym(1,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        / 0.00026,  0.00070,  0.00181,  0.00466,  0.01120,  0.02838  &
        , 0.07167,  0.18415,  0.48459,  0.71213,  0.82704,  0.85937  &
        , 0.81471,  0.87045,  0.91867,  0.94487,  0.96097/

        DATA(gasym(1,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        / 0.00027,  0.00073,  0.00190,  0.00490,  0.01178,  0.02987  &
        , 0.07551,  0.19460,  0.50665,  0.72280,  0.83163,  0.86005  &
        , 0.81373,  0.87385,  0.92031,  0.94620,  0.96179/

        DATA(gasym(1,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        / 0.00029,  0.00077,  0.00200,  0.00518,  0.01249,  0.03171  &
        , 0.08024,  0.20754,  0.53138,  0.73477,  0.83645,  0.86079  &
        , 0.81372,  0.87777,  0.92227,  0.94774,  0.96272/

        DATA(gasym(1,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        / 0.00030,  0.00081,  0.00213,  0.00554,  0.01338,  0.03401  &
        , 0.08622,  0.22405,  0.55867,  0.74775,  0.84127,  0.86054  &
        , 0.81267,  0.88223,  0.92496,  0.94963,  0.96379/

        DATA(gasym(1,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        / 0.00032,  0.00087,  0.00230,  0.00599,  0.01453,  0.03702  &
        , 0.09404,  0.24586,  0.58780,  0.76125,  0.84619,  0.85952  &
        , 0.81445,  0.88748,  0.92866,  0.95195,  0.96494/

        DATA(gasym(1,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        / 0.00034,  0.00094,  0.00252,  0.00661,  0.01609,  0.04110  &
        , 0.10474,  0.27599,  0.61720,  0.77514,  0.85240,  0.85770  &
        , 0.81915,  0.89450,  0.93332,  0.95442,  0.96610/

        DATA(gasym(1,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        / 0.00037,  0.00104,  0.00283,  0.00748,  0.01832,  0.04697  &
        , 0.12025,  0.31991,  0.64533,  0.79081,  0.85971,  0.85294  &
        , 0.82908,  0.90402,  0.93793,  0.95681,  0.96745/

        DATA(gasym(1,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        / 0.00041,  0.00119,  0.00329,  0.00882,  0.02178,  0.05618  &
        , 0.14502,  0.38896,  0.67568,  0.81096,  0.86462,  0.84310  &
        , 0.84766,  0.91541,  0.94207,  0.96024,  0.96890/

        DATA(gasym(1,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        / 0.00048,  0.00143,  0.00408,  0.01119,  0.02799,  0.07291  &
        , 0.19119,  0.50212,  0.72495,  0.83574,  0.86863,  0.82513  &
        , 0.87245,  0.92426,  0.94898,  0.96382,  0.97045/

        DATA(gasym(1,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        / 0.00059,  0.00191,  0.00579,  0.01660,  0.04286,  0.11456  &
        , 0.31283,  0.64580,  0.79281,  0.86243,  0.85747,  0.82773  &
        , 0.90425,  0.93944,  0.95806,  0.96833,  0.97207/

! gasym for Sea Salt:
        DATA(gasym(2,1,ib,1),ib=1,17)  & ! Band 1, RH = 80%
        /0.00372, 0.00962, 0.02456, 0.06213, 0.14883, 0.39219  &
        ,0.65438, 0.78972, 0.83671, 0.78626, 0.78165, 0.83022  &
        ,0.85802, 0.86761, 0.88382, 0.89449, 0.90512/

        DATA(gasym(2,1,ib,2),ib=1,17)  & ! Band 1, RH = 81%
        /0.00383, 0.00989, 0.02521, 0.06379, 0.15294, 0.40271  &
        ,0.65885, 0.79259, 0.83781, 0.78397, 0.78391, 0.83075  &
        ,0.85349, 0.86989, 0.88407, 0.89752, 0.90520/

        DATA(gasym(2,1,ib,3),ib=1,17)  & ! Band 1, RH = 82%
        /0.00394, 0.01016, 0.02590, 0.06557, 0.15735, 0.41380  &
        ,0.66368, 0.79553, 0.83882, 0.78226, 0.78828, 0.83500  &
        ,0.85327, 0.87120, 0.88561, 0.89556, 0.90657/

        DATA(gasym(2,1,ib,4),ib=1,17)  & ! Band 1, RH = 83%
        /0.00405, 0.01044, 0.02665, 0.06747, 0.16210, 0.42554  &
        ,0.66893, 0.79854, 0.83966, 0.78107, 0.78868, 0.83865  &
        ,0.86204, 0.87445, 0.88581, 0.89632, 0.90781/

        DATA(gasym(2,1,ib,5),ib=1,17)  & ! Band 1, RH = 84%
        /0.00417, 0.01075, 0.02745, 0.06955, 0.16728, 0.43802  &
        ,0.67468, 0.80166, 0.84029, 0.77988, 0.79020, 0.83869  &
        ,0.85811, 0.87449, 0.88442, 0.89631, 0.90763/

        DATA(gasym(2,1,ib,6),ib=1,17)  & ! Band 1, RH = 85%
        /0.00429, 0.01109, 0.02833, 0.07182, 0.17299, 0.45137  &
        ,0.68102, 0.80492, 0.84070, 0.77778, 0.79472, 0.84174  &
        ,0.85738, 0.87524, 0.88720, 0.89890, 0.90901/

        DATA(gasym(2,1,ib,7),ib=1,17)  & ! Band 1, RH = 86%
        /0.00443, 0.01146, 0.02930, 0.07433, 0.17934, 0.46568  &
        ,0.68803, 0.80832, 0.84098, 0.77429, 0.79622, 0.83775  &
        ,0.85839, 0.87731, 0.88912, 0.90004, 0.90992/

        DATA(gasym(2,1,ib,8),ib=1,17)  & ! Band 1, RH = 87%
        /0.00458, 0.01187, 0.03038, 0.07715, 0.18650, 0.48106  &
        ,0.69577, 0.81187, 0.84136, 0.77204, 0.79778, 0.84188  &
        ,0.85927, 0.87715, 0.88897, 0.90059, 0.91080/

        DATA(gasym(2,1,ib,9),ib=1,17)  & ! Band 1, RH = 88%
        /0.00474, 0.01233, 0.03160, 0.08035, 0.19467, 0.49758  &
        ,0.70426, 0.81554, 0.84186, 0.77074, 0.79876, 0.84134  &
        ,0.86155, 0.87740, 0.88984, 0.90067, 0.91104/

        DATA(gasym(2,1,ib,10),ib=1,17)  & ! Band 1, RH = 89%
        /0.00493, 0.01285, 0.03300, 0.08403, 0.20411, 0.51526  &
        ,0.71348, 0.81924, 0.84218, 0.76760, 0.80246, 0.83760  &
        ,0.86210, 0.87828, 0.89180, 0.90254, 0.91217/

        DATA(gasym(2,1,ib,11),ib=1,17)  & ! Band 1, RH = 90%
        /0.00515, 0.01346, 0.03463, 0.08832, 0.21519, 0.53402  &
        ,0.72331, 0.82292, 0.84183, 0.76457, 0.80603, 0.84276  &
        ,0.86685, 0.88147, 0.89038, 0.90348, 0.91312/

        DATA(gasym(2,1,ib,12),ib=1,17)  & ! Band 1, RH = 91%
        /0.00540, 0.01417, 0.03655, 0.09340, 0.22839, 0.55361  &
        ,0.73358, 0.82659, 0.84074, 0.76428, 0.80602, 0.84499  &
        ,0.86246, 0.88025, 0.89539, 0.90508, 0.91456/

        DATA(gasym(2,1,ib,13),ib=1,17)  & ! Band 1, RH = 92%
        /0.00570, 0.01502, 0.03885, 0.09953, 0.24443, 0.57363  &
        ,0.74409, 0.83046, 0.83972, 0.76163, 0.81173, 0.84697  &
        ,0.86430, 0.88460, 0.89655, 0.90649, 0.91561/

        DATA(gasym(2,1,ib,14),ib=1,17)  & ! Band 1, RH = 93%
        /0.00606, 0.01605, 0.04167, 0.10709, 0.26435, 0.59349  &
        ,0.75476, 0.83499, 0.83755, 0.76377, 0.81183, 0.84940  &
        ,0.86876, 0.88303, 0.89586, 0.90712, 0.91671/

        DATA(gasym(2,1,ib,15),ib=1,17)  & ! Band 1, RH = 94%
        /0.00651, 0.01735, 0.04522, 0.11669, 0.28982, 0.61262  &
        ,0.76583, 0.84017, 0.83335, 0.76662, 0.81802, 0.85219  &
        ,0.87428, 0.88685, 0.89767, 0.90890, 0.91803/

        DATA(gasym(2,1,ib,16),ib=1,17)  & ! Band 1, RH = 95%
        /0.00707, 0.01902, 0.04987, 0.12941, 0.32355, 0.63113  &
        ,0.77813, 0.84456, 0.82795, 0.77394, 0.82352, 0.85770  &
        ,0.87246, 0.88918, 0.89986, 0.91026, 0.91965/

        DATA(gasym(2,1,ib,17),ib=1,17)  & ! Band 1, RH = 96%
        /0.00782, 0.02128, 0.05627, 0.14722, 0.37024, 0.65138  &
        ,0.79257, 0.84813, 0.81822, 0.78486, 0.82798, 0.86211  &
        ,0.87779, 0.89235, 0.90157, 0.91211, 0.92155/

        DATA(gasym(2,1,ib,18),ib=1,17)  & ! Band 1, RH = 97%
        /0.00888, 0.02459, 0.06582, 0.17448, 0.43773, 0.68085  &
        ,0.80980, 0.85182, 0.80225, 0.80110, 0.84028, 0.86074  &
        ,0.87985, 0.89430, 0.90479, 0.91487, 0.92395/

        DATA(gasym(2,1,ib,19),ib=1,17)  & ! Band 1, RH = 98%
        /0.01054, 0.03001, 0.08207, 0.22284, 0.53416, 0.72946  &
        ,0.82963, 0.85041, 0.77904, 0.81410, 0.84773, 0.86815  &
        ,0.88704, 0.89885, 0.90826, 0.91792, 0.92692/

        DATA(gasym(2,1,ib,20),ib=1,17)  & ! Band 1, RH = 99%
        /0.01363, 0.04120, 0.11833, 0.33813, 0.63653, 0.78644  &
        ,0.84970, 0.82810, 0.78299, 0.82584, 0.86127, 0.87790  &
        ,0.89090, 0.90419, 0.91355, 0.92296, 0.93167/

        DATA(gasym(2,2,ib,1),ib=1,17)  & ! Band 2, RH = 80%
        /0.01683, 0.04331, 0.11050, 0.28850, 0.58369, 0.74438  &
        ,0.81836, 0.80447, 0.73055, 0.77459, 0.81733, 0.83274  &
        ,0.84863, 0.85917, 0.86360, 0.86843, 0.87242/

        DATA(gasym(2,2,ib,2),ib=1,17)  & ! Band 2, RH = 81%
        /0.01733, 0.04452, 0.11345, 0.29658, 0.58950, 0.74770  &
        ,0.82000, 0.80393, 0.73114, 0.77681, 0.81572, 0.83395  &
        ,0.84796, 0.85772, 0.86497, 0.86979, 0.87284/

        DATA(gasym(2,2,ib,3),ib=1,17)  & ! Band 2, RH = 82%
        /0.01781, 0.04572, 0.11661, 0.30519, 0.59529, 0.75111  &
        ,0.82154, 0.80320, 0.73219, 0.78002, 0.81807, 0.83258  &
        ,0.84681, 0.85543, 0.86543, 0.86948, 0.87340/

        DATA(gasym(2,2,ib,4),ib=1,17)  & ! Band 2, RH = 83%
        /0.01830, 0.04701, 0.11999, 0.31444, 0.60112, 0.75462  &
        ,0.82297, 0.80204, 0.73364, 0.78163, 0.82102, 0.83552  &
        ,0.85075, 0.86037, 0.86432, 0.87073, 0.87386/

        DATA(gasym(2,2,ib,5),ib=1,17)  & ! Band 2, RH = 84%
        /0.01882, 0.04839, 0.12366, 0.32443, 0.60705, 0.75828  &
        ,0.82429, 0.80027, 0.73598, 0.78498, 0.81995, 0.83651  &
        ,0.84904, 0.86013, 0.86590, 0.87096, 0.87431/

        DATA(gasym(2,2,ib,6),ib=1,17)  & ! Band 2, RH = 85%
        /0.01938, 0.04989, 0.12767, 0.33533, 0.61318, 0.76211  &
        ,0.82555, 0.79818, 0.73709, 0.79020, 0.82072, 0.83782  &
        ,0.85120, 0.85892, 0.86730, 0.87137, 0.87490/

        DATA(gasym(2,2,ib,7),ib=1,17)  & ! Band 2, RH = 86%
        /0.01999, 0.05154, 0.13211, 0.34732, 0.61962, 0.76613  &
        ,0.82680, 0.79634, 0.73924, 0.79358, 0.82072, 0.83703  &
        ,0.85285, 0.86039, 0.86798, 0.87179, 0.87546/

        DATA(gasym(2,2,ib,8),ib=1,17)  & ! Band 2, RH = 87%
        /0.02067, 0.05338, 0.13709, 0.36063, 0.62650, 0.77037  &
        ,0.82814, 0.79455, 0.74193, 0.79697, 0.82168, 0.84015  &
        ,0.85140, 0.86180, 0.86834, 0.87270, 0.87592/

        DATA(gasym(2,2,ib,9),ib=1,17)  & ! Band 2, RH = 88%
        /0.02142, 0.05545, 0.14272, 0.37552, 0.63402, 0.77485  &
        ,0.82964, 0.79204, 0.74373, 0.80062, 0.82417, 0.84145  &
        ,0.85418, 0.86192, 0.86840, 0.87302, 0.87660/

        DATA(gasym(2,2,ib,10),ib=1,17)  & ! Band 2, RH = 89%
        /0.02228, 0.05780, 0.14920, 0.39230, 0.64240, 0.77962  &
        ,0.83120, 0.78817, 0.74632, 0.80321, 0.82706, 0.84247  &
        ,0.85506, 0.86377, 0.86879, 0.87357, 0.87728/

        DATA(gasym(2,2,ib,11),ib=1,17)  & ! Band 2, RH = 90%
        /0.02325, 0.06052, 0.15673, 0.41130, 0.65193, 0.78473  &
        ,0.83252, 0.78453, 0.74917, 0.80812, 0.82854, 0.84418  &
        ,0.85751, 0.86304, 0.86906, 0.87449, 0.87800/

        DATA(gasym(2,2,ib,12),ib=1,17)  & ! Band 2, RH = 91%
        /0.02439, 0.06372, 0.16566, 0.43288, 0.66296, 0.79028  &
        ,0.83332, 0.78048, 0.75208, 0.80862, 0.83105, 0.84615  &
        ,0.85732, 0.86441, 0.87112, 0.87513, 0.87846/

        DATA(gasym(2,2,ib,13),ib=1,17)  & ! Band 2, RH = 92%
        /0.02574, 0.06753, 0.17642, 0.45740, 0.67579, 0.79633  &
        ,0.83375, 0.77445, 0.75429, 0.81214, 0.83277, 0.84482  &
        ,0.85741, 0.86581, 0.87219, 0.87577, 0.87908/

        DATA(gasym(2,2,ib,14),ib=1,17)  & ! Band 2, RH = 93%
        /0.02736, 0.07218, 0.18969, 0.48504, 0.69060, 0.80288  &
        ,0.83432, 0.76893, 0.75687, 0.81547, 0.83550, 0.84979  &
        ,0.85766, 0.86814, 0.87204, 0.87633, 0.87988/

        DATA(gasym(2,2,ib,15),ib=1,17)  & ! Band 2, RH = 94%
        /0.02936, 0.07799, 0.20652, 0.51569, 0.70721, 0.80974  &
        ,0.83399, 0.76095, 0.76034, 0.81779, 0.83762, 0.85008  &
        ,0.85854, 0.86772, 0.87350, 0.87779, 0.88050/

        DATA(gasym(2,2,ib,16),ib=1,17)  & ! Band 2, RH = 95%
        /0.03190, 0.08551, 0.22871, 0.54865, 0.72499, 0.81699  &
        ,0.83222, 0.75349, 0.76463, 0.81971, 0.83948, 0.85158  &
        ,0.86297, 0.87029, 0.87469, 0.87837, 0.88158/

        DATA(gasym(2,2,ib,17),ib=1,17)  & ! Band 2, RH = 96%
        /0.03528, 0.09574, 0.25954, 0.58256, 0.74341, 0.82556  &
        ,0.82763, 0.74827, 0.77301, 0.82388, 0.84146, 0.85594  &
        ,0.86309, 0.87077, 0.87590, 0.87981, 0.88261/

        DATA(gasym(2,2,ib,18),ib=1,17)  & ! Band 2, RH = 97%
        /0.04004, 0.11070, 0.30572, 0.61709, 0.76374, 0.83377  &
        ,0.81822, 0.74813, 0.78878, 0.82730, 0.84369, 0.85811  &
        ,0.86754, 0.87232, 0.87738, 0.88085, 0.88399/

        DATA(gasym(2,2,ib,19),ib=1,17)  & ! Band 2, RH = 98%
        /0.04746, 0.13545, 0.38245, 0.66046, 0.78895, 0.84005  &
        ,0.79758, 0.75669, 0.81058, 0.83376, 0.84920, 0.86142  &
        ,0.86888, 0.87447, 0.87917, 0.88238, 0.88522/

        DATA(gasym(2,2,ib,20),ib=1,17)  & ! Band 2, RH = 99%
        /0.06134, 0.18747, 0.51997, 0.73291, 0.82151, 0.83610  &
        ,0.75735, 0.77150, 0.82576, 0.84278, 0.85340, 0.86711  &
        ,0.87323, 0.87765, 0.88149, 0.88446, 0.88713/

        DATA(gasym(2,3,ib,1),ib=1,17)  & ! Band 3, RH = 80%
        /0.09028, 0.23390, 0.46961, 0.68055, 0.77385, 0.81365  &
        ,0.76262, 0.72097, 0.80397, 0.82291, 0.83335, 0.84851  &
        ,0.85630, 0.86067, 0.86496, 0.86700, 0.86879/

        DATA(gasym(2,3,ib,2),ib=1,17)  & ! Band 3, RH = 81%
        /0.09305, 0.24036, 0.47702, 0.68471, 0.77665, 0.81460  &
        ,0.76139, 0.72276, 0.80423, 0.82000, 0.83609, 0.84798  &
        ,0.85289, 0.86159, 0.86407, 0.86804, 0.86939/

        DATA(gasym(2,3,ib,3),ib=1,17)  & ! Band 3, RH = 82%
        /0.09572, 0.24669, 0.48469, 0.68890, 0.77955, 0.81541  &
        ,0.75979, 0.72599, 0.80798, 0.82118, 0.83633, 0.84952  &
        ,0.85691, 0.86219, 0.86498, 0.86836, 0.86999/

        DATA(gasym(2,3,ib,4),ib=1,17)  & ! Band 3, RH = 83%
        /0.09841, 0.25336, 0.49265, 0.69317, 0.78244, 0.81583  &
        ,0.75812, 0.72909, 0.80724, 0.82354, 0.83865, 0.84906  &
        ,0.85590, 0.86321, 0.86591, 0.86857, 0.87033/

        DATA(gasym(2,3,ib,5),ib=1,17)  & ! Band 3, RH = 84%
        /0.10129, 0.26041, 0.50094, 0.69756, 0.78538, 0.81662  &
        ,0.75656, 0.73421, 0.80790, 0.82365, 0.84079, 0.85161  &
        ,0.85758, 0.86258, 0.86634, 0.86915, 0.87093/

        DATA(gasym(2,3,ib,6),ib=1,17)  & ! Band 3, RH = 85%
        /0.10440, 0.26793, 0.50964, 0.70212, 0.78843, 0.81709  &
        ,0.75477, 0.73853, 0.80949, 0.82281, 0.84233, 0.84842  &
        ,0.85731, 0.86385, 0.86704, 0.86945, 0.87143/

        DATA(gasym(2,3,ib,7),ib=1,17)  & ! Band 3, RH = 86%
        /0.10780, 0.27598, 0.51878, 0.70693, 0.79167, 0.81744  &
        ,0.75302, 0.74210, 0.80940, 0.82438, 0.84202, 0.84991  &
        ,0.85835, 0.86463, 0.86746, 0.86990, 0.87171/

        DATA(gasym(2,3,ib,8),ib=1,17)  & ! Band 3, RH = 87%
        /0.11154, 0.28468, 0.52844, 0.71205, 0.79497, 0.81768  &
        ,0.75058, 0.74726, 0.80966, 0.82402, 0.84153, 0.85015  &
        ,0.85866, 0.86446, 0.86818, 0.87067, 0.87227/

        DATA(gasym(2,3,ib,9),ib=1,17)  & ! Band 3, RH = 88%
        /0.11573, 0.29412, 0.53869, 0.71758, 0.79852, 0.81734  &
        ,0.74694, 0.75228, 0.81027, 0.82536, 0.84402, 0.85250  &
        ,0.86083, 0.86555, 0.86870, 0.87120, 0.87285/

        DATA(gasym(2,3,ib,10),ib=1,17)  & ! Band 3, RH = 89%
        /0.12046, 0.30444, 0.54964, 0.72354, 0.80224, 0.81720  &
        ,0.74357, 0.75964, 0.81048, 0.82748, 0.84634, 0.85461  &
        ,0.86036, 0.86708, 0.86947, 0.87178, 0.87323/

        DATA(gasym(2,3,ib,11),ib=1,17)  & ! Band 3, RH = 90%
        /0.12590, 0.31578, 0.56142, 0.72983, 0.80579, 0.81699  &
        ,0.73982, 0.76480, 0.81142, 0.82817, 0.84298, 0.85640  &
        ,0.86317, 0.86679, 0.87007, 0.87243, 0.87392/

        DATA(gasym(2,3,ib,12),ib=1,17)  & ! Band 3, RH = 91%
        /0.13222, 0.32835, 0.57429, 0.73649, 0.80936, 0.81594  &
        ,0.73664, 0.77042, 0.81320, 0.82963, 0.84684, 0.85505  &
        ,0.86363, 0.86739, 0.87048, 0.87279, 0.87444/

        DATA(gasym(2,3,ib,13),ib=1,17)  & ! Band 3, RH = 92%
        /0.13971, 0.34241, 0.58859, 0.74353, 0.81298, 0.81400  &
        ,0.73250, 0.77908, 0.81474, 0.83287, 0.84788, 0.85759  &
        ,0.86347, 0.86788, 0.87134, 0.87374, 0.87492/

        DATA(gasym(2,3,ib,14),ib=1,17)  & ! Band 3, RH = 93%
        /0.14874, 0.35842, 0.60473, 0.75111, 0.81661, 0.81155  &
        ,0.72798, 0.78519, 0.81856, 0.83556, 0.84886, 0.85778  &
        ,0.86367, 0.86874, 0.87203, 0.87436, 0.87543/

        DATA(gasym(2,3,ib,15),ib=1,17)  & ! Band 3, RH = 94%
        /0.15990, 0.37718, 0.62322, 0.75973, 0.82003, 0.80769  &
        ,0.72424, 0.79287, 0.82129, 0.83943, 0.85214, 0.86049  &
        ,0.86633, 0.87005, 0.87296, 0.87494, 0.87633/

        DATA(gasym(2,3,ib,16),ib=1,17)  & ! Band 3, RH = 95%
        /0.17410, 0.40022, 0.64446, 0.76933, 0.82368, 0.80243  &
        ,0.72191, 0.80250, 0.82666, 0.84324, 0.85335, 0.86023  &
        ,0.86645, 0.87120, 0.87388, 0.87565, 0.87709/

        DATA(gasym(2,3,ib,17),ib=1,17)  & ! Band 3, RH = 96%
        /0.19289, 0.43017, 0.66811, 0.78087, 0.82737, 0.79334  &
        ,0.72616, 0.81112, 0.83311, 0.84324, 0.85502, 0.86432  &
        ,0.86835, 0.87185, 0.87464, 0.87662, 0.87780/

        DATA(gasym(2,3,ib,18),ib=1,17)  & ! Band 3, RH = 97%
        /0.21911, 0.47083, 0.69366, 0.79592, 0.82979, 0.77907  &
        ,0.73832, 0.81622, 0.83326, 0.84596, 0.85760, 0.86539  &
        ,0.87008, 0.87354, 0.87586, 0.87724, 0.87857/

        DATA(gasym(2,3,ib,19),ib=1,17)  & ! Band 3, RH = 98%
        /0.25817, 0.52647, 0.72383, 0.81342, 0.82758, 0.75498  &
        ,0.76788, 0.81513, 0.83501, 0.84981, 0.85871, 0.86711  &
        ,0.87169, 0.87480, 0.87680, 0.87850, 0.87944/

        DATA(gasym(2,3,ib,20),ib=1,17)  & ! Band 3, RH = 99%
        /0.32054, 0.60325, 0.76351, 0.82916, 0.80845, 0.72768  &
        ,0.80522, 0.83202, 0.84412, 0.85774, 0.86561, 0.86980  &
        ,0.87372, 0.87640, 0.87834, 0.87969, 0.88053/

        DATA(gasym(2,4,ib,1),ib=1,17)  & ! Band 4, RH = 80%
        /0.00002, 0.00006, 0.00015, 0.00039, 0.00093, 0.00236  &
        ,0.00596, 0.01501, 0.03787, 0.09733, 0.26651, 0.58059  &
        ,0.76358, 0.85096, 0.89449, 0.91873, 0.92950/

        DATA(gasym(2,4,ib,2),ib=1,17)  & ! Band 4, RH = 81%
        /0.00002, 0.00006, 0.00016, 0.00040, 0.00096, 0.00243  &
        ,0.00612, 0.01543, 0.03893, 0.10018, 0.27471, 0.58770  &
        ,0.76704, 0.85241, 0.89521, 0.91890, 0.92941/

        DATA(gasym(2,4,ib,3),ib=1,17)  & ! Band 4, RH = 82%
        /0.00002, 0.00006, 0.00016, 0.00041, 0.00099, 0.00250  &
        ,0.00630, 0.01587, 0.04006, 0.10323, 0.28346, 0.59493  &
        ,0.77054, 0.85388, 0.89595, 0.91907, 0.92933/

        DATA(gasym(2,4,ib,4),ib=1,17)  & ! Band 4, RH = 83%
        /0.00003, 0.00007, 0.00017, 0.00043, 0.00102, 0.00257  &
        ,0.00649, 0.01635, 0.04128, 0.10654, 0.29287, 0.60231  &
        ,0.77411, 0.85540, 0.89672, 0.91925, 0.92925/

        DATA(gasym(2,4,ib,5),ib=1,17)  & ! Band 4, RH = 84%
        /0.00003, 0.00007, 0.00017, 0.00044, 0.00105, 0.00265  &
        ,0.00670, 0.01688, 0.04262, 0.11017, 0.30307, 0.60991  &
        ,0.77778, 0.85696, 0.89753, 0.91944, 0.92917/

        DATA(gasym(2,4,ib,6),ib=1,17)  & ! Band 4, RH = 85%
        /0.00003, 0.00007, 0.00018, 0.00045, 0.00108, 0.00274  &
        ,0.00693, 0.01745, 0.04409, 0.11418, 0.31420, 0.61781  &
        ,0.78157, 0.85860, 0.89837, 0.91964, 0.92910/

        DATA(gasym(2,4,ib,7),ib=1,17)  & ! Band 4, RH = 86%
        /0.00003, 0.00007, 0.00018, 0.00047, 0.00112, 0.00284  &
        ,0.00718, 0.01809, 0.04572, 0.11868, 0.32645, 0.62609  &
        ,0.78551, 0.86032, 0.89926, 0.91985, 0.92903/

        DATA(gasym(2,4,ib,8),ib=1,17)  & ! Band 4, RH = 87%
        /0.00003, 0.00007, 0.00019, 0.00049, 0.00117, 0.00296  &
        ,0.00746, 0.01881, 0.04756, 0.12376, 0.34003, 0.63484  &
        ,0.78964, 0.86216, 0.90021, 0.92007, 0.92896/

        DATA(gasym(2,4,ib,9),ib=1,17)  & ! Band 4, RH = 88%
        /0.00003, 0.00008, 0.00020, 0.00051, 0.00122, 0.00308  &
        ,0.00779, 0.01963, 0.04967, 0.12960, 0.35519, 0.64417  &
        ,0.79402, 0.86412, 0.90122, 0.92032, 0.92891/

        DATA(gasym(2,4,ib,10),ib=1,17)  & ! Band 4, RH = 89%
        /0.00003, 0.00008, 0.00021, 0.00053, 0.00128, 0.00323  &
        ,0.00816, 0.02057, 0.05210, 0.13638, 0.37221, 0.65418  &
        ,0.79868, 0.86623, 0.90230, 0.92058, 0.92886/

        DATA(gasym(2,4,ib,11),ib=1,17)  & ! Band 4, RH = 90%
        /0.00003, 0.00008, 0.00022, 0.00056, 0.00134, 0.00340  &
        ,0.00859, 0.02168, 0.05495, 0.14440, 0.39140, 0.66501  &
        ,0.80369, 0.86852, 0.90347, 0.92087, 0.92882/

        DATA(gasym(2,4,ib,12),ib=1,17)  & ! Band 4, RH = 91%
        /0.00003, 0.00009, 0.00023, 0.00059, 0.00142, 0.00360  &
        ,0.00911, 0.02299, 0.05834, 0.15402, 0.41310, 0.67679  &
        ,0.80911, 0.87101, 0.90473, 0.92119, 0.92880/

        DATA(gasym(2,4,ib,13),ib=1,17)  & ! Band 4, RH = 92%
        /0.00004, 0.00009, 0.00025, 0.00063, 0.00152, 0.00385  &
        ,0.00974, 0.02458, 0.06247, 0.16580, 0.43768, 0.68968  &
        ,0.81498, 0.87378, 0.90609, 0.92155, 0.92880/

        DATA(gasym(2,4,ib,14),ib=1,17)  & ! Band 4, RH = 93%
        /0.00004, 0.00010, 0.00026, 0.00068, 0.00164, 0.00415  &
        ,0.01051, 0.02655, 0.06761, 0.18061, 0.46549, 0.70386  &
        ,0.82136, 0.87688, 0.90758, 0.92196, 0.92881/

        DATA(gasym(2,4,ib,15),ib=1,17)  & ! Band 4, RH = 94%
        /0.00004, 0.00011, 0.00029, 0.00074, 0.00179, 0.00454  &
        ,0.01150, 0.02906, 0.07421, 0.19980, 0.49688, 0.71954  &
        ,0.82832, 0.88039, 0.90922, 0.92244, 0.92886/

        DATA(gasym(2,4,ib,16),ib=1,17)  & ! Band 4, RH = 95%
        /0.00004, 0.00012, 0.00032, 0.00082, 0.00198, 0.00505  &
        ,0.01281, 0.03241, 0.08309, 0.22574, 0.53221, 0.73703  &
        ,0.83608, 0.88442, 0.91102, 0.92300, 0.92895/

        DATA(gasym(2,4,ib,17),ib=1,17)  & ! Band 4, RH = 96%
        /0.00005, 0.00013, 0.00036, 0.00093, 0.00226, 0.00577  &
        ,0.01466, 0.03715, 0.09581, 0.26272, 0.57203, 0.75674  &
        ,0.84497, 0.88915, 0.91304, 0.92369, 0.92911/

        DATA(gasym(2,4,ib,18),ib=1,17)  & ! Band 4, RH = 97%
        /0.00006, 0.00016, 0.00042, 0.00110, 0.00269, 0.00688  &
        ,0.01750, 0.04450, 0.11590, 0.31919, 0.61788, 0.77927  &
        ,0.85536, 0.89489, 0.91538, 0.92457, 0.92935/

        DATA(gasym(2,4,ib,19),ib=1,17)  & ! Band 4, RH = 98%
        /0.00007, 0.00019, 0.00052, 0.00140, 0.00343, 0.00885  &
        ,0.02260, 0.05782, 0.15349, 0.41168, 0.67448, 0.80640  &
        ,0.86838, 0.90202, 0.91820, 0.92575, 0.92975/

        DATA(gasym(2,4,ib,20),ib=1,17)  & ! Band 4, RH = 99%
        /0.00009, 0.00026, 0.00075, 0.00208, 0.00522, 0.01365  &
        ,0.03525, 0.09188, 0.25351, 0.56302, 0.75191, 0.84210  &
        ,0.88706, 0.91128, 0.92206, 0.92754, 0.93046/

        DATA(gasym(2,5,ib,1),ib=1,17)  & ! Band 5, RH = 80%
        /0.00005, 0.00013, 0.00034, 0.00087, 0.00209, 0.00528  &
        ,0.01332, 0.03371, 0.08674, 0.23542, 0.56870, 0.77651  &
        ,0.87177, 0.91307, 0.93289, 0.94285, 0.94791/

        DATA(gasym(2,5,ib,2),ib=1,17)  & ! Band 5, RH = 81%
        /0.00005, 0.00014, 0.00035, 0.00090, 0.00214, 0.00542  &
        ,0.01368, 0.03462, 0.08919, 0.24252, 0.57811, 0.78044  &
        ,0.87338, 0.91371, 0.93313, 0.94287, 0.94784/

        DATA(gasym(2,5,ib,3),ib=1,17)  & ! Band 5, RH = 82%
        /0.00005, 0.00014, 0.00036, 0.00092, 0.00220, 0.00557  &
        ,0.01406, 0.03560, 0.09182, 0.25015, 0.58765, 0.78446  &
        ,0.87502, 0.91437, 0.93337, 0.94289, 0.94777/

        DATA(gasym(2,5,ib,4),ib=1,17)  & ! Band 5, RH = 83%
        /0.00006, 0.00015, 0.00037, 0.00095, 0.00227, 0.00573  &
        ,0.01447, 0.03665, 0.09467, 0.25839, 0.59738, 0.78858  &
        ,0.87669, 0.91504, 0.93362, 0.94292, 0.94770/

        DATA(gasym(2,5,ib,5),ib=1,17)  & ! Band 5, RH = 84%
        /0.00006, 0.00015, 0.00038, 0.00098, 0.00233, 0.00591  &
        ,0.01492, 0.03780, 0.09779, 0.26739, 0.60733, 0.79283  &
        ,0.87842, 0.91574, 0.93387, 0.94295, 0.94763/

        DATA(gasym(2,5,ib,6),ib=1,17)  & ! Band 5, RH = 85%
        /0.00006, 0.00015, 0.00040, 0.00101, 0.00241, 0.00610  &
        ,0.01541, 0.03907, 0.10124, 0.27729, 0.61755, 0.79725  &
        ,0.88021, 0.91646, 0.93414, 0.94298, 0.94756/

        DATA(gasym(2,5,ib,7),ib=1,17)  & ! Band 5, RH = 86%
        /0.00006, 0.00016, 0.00041, 0.00104, 0.00250, 0.00631  &
        ,0.01596, 0.04048, 0.10510, 0.28831, 0.62809, 0.80185  &
        ,0.88208, 0.91722, 0.93442, 0.94302, 0.94749/

        DATA(gasym(2,5,ib,8),ib=1,17)  & ! Band 5, RH = 87%
        /0.00006, 0.00017, 0.00042, 0.00108, 0.00259, 0.00656  &
        ,0.01658, 0.04208, 0.10947, 0.30069, 0.63900, 0.80669  &
        ,0.88403, 0.91803, 0.93472, 0.94307, 0.94742/

        DATA(gasym(2,5,ib,9),ib=1,17)  & ! Band 5, RH = 88%
        /0.00007, 0.00017, 0.00044, 0.00113, 0.00270, 0.00683  &
        ,0.01728, 0.04390, 0.11449, 0.31472, 0.65032, 0.81177  &
        ,0.88609, 0.91888, 0.93504, 0.94312, 0.94735/

        DATA(gasym(2,5,ib,10),ib=1,17)  & ! Band 5, RH = 89%
        /0.00007, 0.00018, 0.00046, 0.00118, 0.00282, 0.00715  &
        ,0.01809, 0.04600, 0.12031, 0.33079, 0.66211, 0.81714  &
        ,0.88827, 0.91979, 0.93538, 0.94319, 0.94729/

        DATA(gasym(2,5,ib,11),ib=1,17)  & ! Band 5, RH = 90%
        /0.00007, 0.00019, 0.00048, 0.00124, 0.00297, 0.00752  &
        ,0.01904, 0.04847, 0.12718, 0.34936, 0.67443, 0.82283  &
        ,0.89059, 0.92077, 0.93575, 0.94326, 0.94722/

        DATA(gasym(2,5,ib,12),ib=1,17)  & ! Band 5, RH = 91%
        /0.00007, 0.00020, 0.00051, 0.00131, 0.00314, 0.00796  &
        ,0.02017, 0.05141, 0.13542, 0.37103, 0.68737, 0.82888  &
        ,0.89306, 0.92184, 0.93616, 0.94336, 0.94716/

        DATA(gasym(2,5,ib,13),ib=1,17)  & ! Band 5, RH = 92%
        /0.00008, 0.00021, 0.00054, 0.00139, 0.00335, 0.00849  &
        ,0.02153, 0.05499, 0.14552, 0.39654, 0.70108, 0.83536  &
        ,0.89571, 0.92301, 0.93661, 0.94347, 0.94711/

        DATA(gasym(2,5,ib,14),ib=1,17)  & ! Band 5, RH = 93%
        /0.00008, 0.00022, 0.00058, 0.00150, 0.00360, 0.00915  &
        ,0.02323, 0.05944, 0.15818, 0.42683, 0.71583, 0.84233  &
        ,0.89861, 0.92429, 0.93711, 0.94360, 0.94706/

        DATA(gasym(2,5,ib,15),ib=1,17)  & ! Band 5, RH = 94%
        /0.00009, 0.00024, 0.00063, 0.00163, 0.00393, 0.01000  &
        ,0.02540, 0.06516, 0.17461, 0.46304, 0.73202, 0.84991  &
        ,0.90180, 0.92572, 0.93768, 0.94377, 0.94702/

        DATA(gasym(2,5,ib,16),ib=1,17)  & ! Band 5, RH = 95%
        /0.00010, 0.00026, 0.00070, 0.00180, 0.00436, 0.01111  &
        ,0.02828, 0.07285, 0.19685, 0.50646, 0.75031, 0.85825  &
        ,0.90536, 0.92733, 0.93833, 0.94398, 0.94701/

        DATA(gasym(2,5,ib,17),ib=1,17)  & ! Band 5, RH = 96%
        /0.00011, 0.00030, 0.00078, 0.00205, 0.00496, 0.01268  &
        ,0.03236, 0.08384, 0.22883, 0.55822, 0.77157, 0.86758  &
        ,0.90937, 0.92918, 0.93911, 0.94425, 0.94701/

        DATA(gasym(2,5,ib,18),ib=1,17)  & ! Band 5, RH = 97%
        /0.00012, 0.00034, 0.00092, 0.00241, 0.00588, 0.01510  &
        ,0.03868, 0.10114, 0.27893, 0.61835, 0.79685, 0.87836  &
        ,0.91411, 0.93139, 0.94007, 0.94462, 0.94706/

        DATA(gasym(2,5,ib,19),ib=1,17)  & ! Band 5, RH = 98%
        /0.00015, 0.00042, 0.00114, 0.00305, 0.00751, 0.01941  &
        ,0.05012, 0.13328, 0.36749, 0.68528, 0.82712, 0.89117  &
        ,0.91994, 0.93414, 0.94132, 0.94514, 0.94719/

        DATA(gasym(2,5,ib,20),ib=1,17)  & ! Band 5, RH = 99%
        /0.00019, 0.00057, 0.00165, 0.00454, 0.01142, 0.02999  &
        ,0.07911, 0.21829, 0.54490, 0.76654, 0.86496, 0.90770  &
        ,0.92785, 0.93794, 0.94316, 0.94596, 0.94746/

        DATA(gasym(2,6,ib,1),ib=1,17)  & ! Band 6, RH = 80%
        /0.00015, 0.00039, 0.00099, 0.00253, 0.00603, 0.01522  &
        ,0.03838, 0.09715, 0.25615, 0.63324, 0.78852, 0.87209  &
        ,0.88964, 0.84120, 0.93843, 0.96156, 0.97038/

        DATA(gasym(2,6,ib,2),ib=1,17)  & ! Band 6, RH = 81%
        /0.00015, 0.00040, 0.00102, 0.00259, 0.00619, 0.01564  &
        ,0.03942, 0.09985, 0.26384, 0.64208, 0.79358, 0.87364  &
        ,0.89017, 0.84265, 0.93975, 0.96205, 0.97071/

        DATA(gasym(2,6,ib,3),ib=1,17)  & ! Band 6, RH = 82%
        /0.00016, 0.00041, 0.00105, 0.00267, 0.00636, 0.01608  &
        ,0.04053, 0.10274, 0.27211, 0.65060, 0.79844, 0.87498  &
        ,0.89055, 0.84462, 0.94109, 0.96251, 0.97103/

        DATA(gasym(2,6,ib,4),ib=1,17)  & ! Band 6, RH = 83%
        /0.00016, 0.00042, 0.00108, 0.00274, 0.00655, 0.01655  &
        ,0.04174, 0.10586, 0.28109, 0.65875, 0.80305, 0.87617  &
        ,0.89083, 0.84727, 0.94241, 0.96298, 0.97134/

        DATA(gasym(2,6,ib,5),ib=1,17)  & ! Band 6, RH = 84%
        /0.00017, 0.00043, 0.00111, 0.00283, 0.00675, 0.01707  &
        ,0.04305, 0.10928, 0.29091, 0.66651, 0.80743, 0.87730  &
        ,0.89111, 0.85019, 0.94369, 0.96344, 0.97163/

        DATA(gasym(2,6,ib,6),ib=1,17)  & ! Band 6, RH = 85%
        /0.00017, 0.00045, 0.00115, 0.00292, 0.00698, 0.01764  &
        ,0.04449, 0.11304, 0.30177, 0.67382, 0.81156, 0.87854  &
        ,0.89151, 0.85308, 0.94491, 0.96392, 0.97192/

        DATA(gasym(2,6,ib,7),ib=1,17)  & ! Band 6, RH = 86%
        /0.00018, 0.00046, 0.00119, 0.00302, 0.00723, 0.01827  &
        ,0.04609, 0.11724, 0.31391, 0.68062, 0.81549, 0.88011  &
        ,0.89199, 0.85652, 0.94615, 0.96443, 0.97220/

        DATA(gasym(2,6,ib,8),ib=1,17)  & ! Band 6, RH = 87%
        /0.00018, 0.00048, 0.00123, 0.00314, 0.00750, 0.01898  &
        ,0.04790, 0.12198, 0.32762, 0.68687, 0.81928, 0.88226  &
        ,0.89225, 0.86119, 0.94743, 0.96499, 0.97248/

        DATA(gasym(2,6,ib,9),ib=1,17)  & ! Band 6, RH = 88%
        /0.00019, 0.00050, 0.00128, 0.00327, 0.00782, 0.01978  &
        ,0.04996, 0.12739, 0.34328, 0.69254, 0.82306, 0.88515  &
        ,0.89195, 0.86634, 0.94871, 0.96562, 0.97276/

        DATA(gasym(2,6,ib,10),ib=1,17)  & ! Band 6, RH = 89%
        /0.00020, 0.00052, 0.00134, 0.00342, 0.00818, 0.02071  &
        ,0.05233, 0.13366, 0.36137, 0.69770, 0.82703, 0.88868  &
        ,0.89111, 0.87203, 0.95005, 0.96634, 0.97304/

        DATA(gasym(2,6,ib,11),ib=1,17)  & ! Band 6, RH = 90%
        /0.00021, 0.00054, 0.00140, 0.00359, 0.00861, 0.02180  &
        ,0.05510, 0.14104, 0.38249, 0.70253, 0.83145, 0.89226  &
        ,0.89044, 0.87947, 0.95156, 0.96718, 0.97335/

        DATA(gasym(2,6,ib,12),ib=1,17)  & ! Band 6, RH = 91%
        /0.00022, 0.00057, 0.00148, 0.00380, 0.00911, 0.02308  &
        ,0.05840, 0.14985, 0.40742, 0.70753, 0.83661, 0.89495  &
        ,0.89006, 0.88727, 0.95321, 0.96813, 0.97368/

        DATA(gasym(2,6,ib,13),ib=1,17)  & ! Band 6, RH = 92%
        /0.00023, 0.00061, 0.00158, 0.00405, 0.00972, 0.02464  &
        ,0.06240, 0.16061, 0.43714, 0.71363, 0.84271, 0.89630  &
        ,0.88865, 0.89695, 0.95518, 0.96918, 0.97402/

        DATA(gasym(2,6,ib,14),ib=1,17)  & ! Band 6, RH = 93%
        /0.00024, 0.00065, 0.00169, 0.00435, 0.01046, 0.02657  &
        ,0.06735, 0.17405, 0.47285, 0.72241, 0.84971, 0.89759  &
        ,0.88628, 0.90700, 0.95756, 0.97026, 0.97436/

        DATA(gasym(2,6,ib,15),ib=1,17)  & ! Band 6, RH = 94%
        /0.00026, 0.00070, 0.00184, 0.00474, 0.01142, 0.02902  &
        ,0.07368, 0.19144, 0.51580, 0.73629, 0.85764, 0.90077  &
        ,0.88342, 0.91813, 0.96038, 0.97128, 0.97471/

        DATA(gasym(2,6,ib,16),ib=1,17)  & ! Band 6, RH = 95%
        /0.00029, 0.00077, 0.00202, 0.00525, 0.01268, 0.03227  &
        ,0.08213, 0.21495, 0.56668, 0.75791, 0.86722, 0.90382  &
        ,0.87974, 0.92950, 0.96351, 0.97218, 0.97507/

        DATA(gasym(2,6,ib,17),ib=1,17)  & ! Band 6, RH = 96%
        /0.00032, 0.00086, 0.00229, 0.00596, 0.01444, 0.03685  &
        ,0.09410, 0.24885, 0.62348, 0.78727, 0.87840, 0.90516  &
        ,0.87546, 0.94040, 0.96647, 0.97308, 0.97546/

        DATA(gasym(2,6,ib,18),ib=1,17)  & ! Band 6, RH = 97%
        /0.00036, 0.00099, 0.00268, 0.00703, 0.01713, 0.04389  &
        ,0.11272, 0.30262, 0.67653, 0.81631, 0.88584, 0.90556  &
        ,0.87435, 0.94957, 0.96848, 0.97415, 0.97586/

        DATA(gasym(2,6,ib,19),ib=1,17)  & ! Band 6, RH = 98%
        /0.00043, 0.00121, 0.00334, 0.00890, 0.02188, 0.05646  &
        ,0.14667, 0.40115, 0.71016, 0.83952, 0.89962, 0.90192  &
        ,0.89116, 0.95652, 0.97015, 0.97506, 0.97630/

        DATA(gasym(2,6,ib,20),ib=1,17)  & ! Band 6, RH = 99%
        /0.00055, 0.00167, 0.00480, 0.01323, 0.03325, 0.08738  &
        ,0.23464, 0.60652, 0.78084, 0.87799, 0.90836, 0.88433  &
        ,0.93891, 0.96699, 0.97369, 0.97604, 0.97682/

        DATA(gasym(2,7,ib,1),ib=1,17)  & ! Band 7, RH = 80%
        /0.00012, 0.00032, 0.00081, 0.00206, 0.00492, 0.01242  &
        ,0.03134, 0.07927, 0.20624, 0.54493, 0.75160, 0.86569  &
        ,0.90725, 0.89738, 0.92787, 0.96496, 0.97336/

        DATA(gasym(2,7,ib,2),ib=1,17)  & ! Band 7, RH = 81%
        /0.00013, 0.00033, 0.00083, 0.00212, 0.00505, 0.01276  &
        ,0.03219, 0.08146, 0.21230, 0.55717, 0.75707, 0.86831  &
        ,0.90831, 0.89791, 0.93035, 0.96570, 0.97367/

        DATA(gasym(2,7,ib,3),ib=1,17)  & ! Band 7, RH = 82%
        /0.00013, 0.00033, 0.00085, 0.00217, 0.00519, 0.01312  &
        ,0.03310, 0.08381, 0.21882, 0.56967, 0.76288, 0.87101  &
        ,0.90938, 0.89839, 0.93276, 0.96644, 0.97397/

        DATA(gasym(2,7,ib,4),ib=1,17)  & ! Band 7, RH = 83%
        /0.00013, 0.00034, 0.00088, 0.00224, 0.00534, 0.01351  &
        ,0.03409, 0.08634, 0.22589, 0.58246, 0.76904, 0.87377  &
        ,0.91049, 0.89892, 0.93513, 0.96716, 0.97427/

        DATA(gasym(2,7,ib,5),ib=1,17)  & ! Band 7, RH = 84%
        /0.00014, 0.00035, 0.00091, 0.00231, 0.00551, 0.01393  &
        ,0.03516, 0.08911, 0.23363, 0.59553, 0.77554, 0.87657  &
        ,0.91164, 0.89960, 0.93757, 0.96787, 0.97458/

        DATA(gasym(2,7,ib,6),ib=1,17)  & ! Band 7, RH = 85%
        /0.00014, 0.00036, 0.00094, 0.00238, 0.00569, 0.01439  &
        ,0.03633, 0.09216, 0.24219, 0.60888, 0.78236, 0.87937  &
        ,0.91278, 0.90033, 0.94001, 0.96859, 0.97489/

        DATA(gasym(2,7,ib,7),ib=1,17)  & ! Band 7, RH = 86%
        /0.00015, 0.00038, 0.00097, 0.00247, 0.00589, 0.01491  &
        ,0.03765, 0.09556, 0.25177, 0.62247, 0.78947, 0.88212  &
        ,0.91387, 0.90095, 0.94245, 0.96930, 0.97520/

        DATA(gasym(2,7,ib,8),ib=1,17)  & ! Band 7, RH = 87%
        /0.00015, 0.00039, 0.00100, 0.00256, 0.00612, 0.01549  &
        ,0.03912, 0.09939, 0.26261, 0.63619, 0.79678, 0.88476  &
        ,0.91489, 0.90152, 0.94495, 0.97000, 0.97552/

        DATA(gasym(2,7,ib,9),ib=1,17)  & ! Band 7, RH = 88%
        /0.00016, 0.00041, 0.00104, 0.00267, 0.00638, 0.01615  &
        ,0.04080, 0.10377, 0.27503, 0.64990, 0.80419, 0.88729  &
        ,0.91592, 0.90236, 0.94742, 0.97069, 0.97584/

        DATA(gasym(2,7,ib,10),ib=1,17)  & ! Band 7, RH = 89%
        /0.00016, 0.00042, 0.00109, 0.00279, 0.00667, 0.01691  &
        ,0.04274, 0.10884, 0.28944, 0.66337, 0.81160, 0.88978  &
        ,0.91707, 0.90341, 0.94999, 0.97138, 0.97618/

        DATA(gasym(2,7,ib,11),ib=1,17)  & ! Band 7, RH = 90%
        /0.00017, 0.00044, 0.00114, 0.00293, 0.00702, 0.01779  &
        ,0.04501, 0.11478, 0.30638, 0.67631, 0.81890, 0.89248  &
        ,0.91832, 0.90460, 0.95251, 0.97206, 0.97651/

        DATA(gasym(2,7,ib,12),ib=1,17)  & ! Band 7, RH = 91%
        /0.00018, 0.00047, 0.00121, 0.00310, 0.00743, 0.01884  &
        ,0.04770, 0.12188, 0.32660, 0.68839, 0.82604, 0.89574  &
        ,0.91935, 0.90627, 0.95509, 0.97273, 0.97685/

        DATA(gasym(2,7,ib,13),ib=1,17)  & ! Band 7, RH = 92%
        /0.00019, 0.00049, 0.00128, 0.00330, 0.00793, 0.02012  &
        ,0.05096, 0.13052, 0.35111, 0.69940, 0.83315, 0.89978  &
        ,0.92009, 0.90837, 0.95765, 0.97340, 0.97719/

        DATA(gasym(2,7,ib,14),ib=1,17)  & ! Band 7, RH = 93%
        /0.00020, 0.00053, 0.00138, 0.00355, 0.00854, 0.02169  &
        ,0.05500, 0.14129, 0.38140, 0.70944, 0.84056, 0.90424  &
        ,0.92108, 0.91157, 0.96021, 0.97409, 0.97754/

        DATA(gasym(2,7,ib,15),ib=1,17)  & ! Band 7, RH = 94%
        /0.00021, 0.00057, 0.00150, 0.00387, 0.00931, 0.02369  &
        ,0.06017, 0.15518, 0.41954, 0.71947, 0.84878, 0.90832  &
        ,0.92171, 0.91604, 0.96280, 0.97481, 0.97790/

        DATA(gasym(2,7,ib,16),ib=1,17)  & ! Band 7, RH = 95%
        /0.00023, 0.00063, 0.00165, 0.00428, 0.01034, 0.02635  &
        ,0.06706, 0.17388, 0.46850, 0.73214, 0.85832, 0.91235  &
        ,0.92184, 0.92243, 0.96547, 0.97559, 0.97827/

        DATA(gasym(2,7,ib,17),ib=1,17)  & ! Band 7, RH = 96%
        /0.00026, 0.00070, 0.00186, 0.00486, 0.01178, 0.03010  &
        ,0.07680, 0.20073, 0.53176, 0.75268, 0.87005, 0.91715  &
        ,0.92147, 0.93142, 0.96834, 0.97643, 0.97867/

        DATA(gasym(2,7,ib,18),ib=1,17)  & ! Band 7, RH = 97%
        /0.00029, 0.00081, 0.00218, 0.00574, 0.01398, 0.03585  &
        ,0.09192, 0.24320, 0.61010, 0.78691, 0.88521, 0.92196  &
        ,0.92012, 0.94321, 0.97149, 0.97729, 0.97910/

        DATA(gasym(2,7,ib,19),ib=1,17)  & ! Band 7, RH = 98%
        /0.00035, 0.00099, 0.00272, 0.00726, 0.01786, 0.04611  &
        ,0.11935, 0.32192, 0.68801, 0.82792, 0.89941, 0.92626  &
        ,0.91876, 0.95653, 0.97455, 0.97825, 0.97957/

        DATA(gasym(2,7,ib,20),ib=1,17)  & ! Band 7, RH = 99%
        /0.00045, 0.00136, 0.00392, 0.01080, 0.02716, 0.07134  &
        ,0.18956, 0.51157, 0.74833, 0.86898, 0.91861, 0.92646  &
        ,0.93225, 0.96876, 0.97709, 0.97928, 0.98014/

        DATA(gasym(2,8,ib,1),ib=1,17)  & ! Band 8, RH = 80%
        /0.00042, 0.00108, 0.00277, 0.00703, 0.01674, 0.04215  &
        ,0.10638, 0.27864, 0.61643, 0.76936, 0.84399, 0.83886  &
        ,0.79308, 0.88396, 0.92410, 0.94751, 0.96205/

        DATA(gasym(2,8,ib,2),ib=1,17)  & ! Band 8, RH = 81%
        /0.00043, 0.00111, 0.00284, 0.00722, 0.01719, 0.04329  &
        ,0.10931, 0.28684, 0.62246, 0.77267, 0.84581, 0.83798  &
        ,0.79638, 0.88633, 0.92539, 0.94838, 0.96257/

        DATA(gasym(2,8,ib,3),ib=1,17)  & ! Band 8, RH = 82%
        /0.00044, 0.00114, 0.00292, 0.00742, 0.01767, 0.04451  &
        ,0.11246, 0.29563, 0.62835, 0.77607, 0.84764, 0.83719  &
        ,0.80023, 0.88856, 0.92659, 0.94927, 0.96307/

        DATA(gasym(2,8,ib,4),ib=1,17)  & ! Band 8, RH = 83%
        /0.00045, 0.00117, 0.00300, 0.00764, 0.01819, 0.04583  &
        ,0.11586, 0.30513, 0.63414, 0.77959, 0.84941, 0.83661  &
        ,0.80440, 0.89104, 0.92778, 0.95020, 0.96356/

        DATA(gasym(2,8,ib,5),ib=1,17)  & ! Band 8, RH = 84%
        /0.00047, 0.00121, 0.00309, 0.00787, 0.01876, 0.04726  &
        ,0.11956, 0.31548, 0.63991, 0.78328, 0.85105, 0.83614  &
        ,0.80821, 0.89378, 0.92893, 0.95117, 0.96405/

        DATA(gasym(2,8,ib,6),ib=1,17)  & ! Band 8, RH = 85%
        /0.00048, 0.00125, 0.00319, 0.00813, 0.01938, 0.04884  &
        ,0.12365, 0.32687, 0.64572, 0.78721, 0.85251, 0.83539  &
        ,0.81239, 0.89631, 0.93011, 0.95221, 0.96454/

        DATA(gasym(2,8,ib,7),ib=1,17)  & ! Band 8, RH = 86%
        /0.00050, 0.00129, 0.00330, 0.00841, 0.02006, 0.05059  &
        ,0.12821, 0.33953, 0.65170, 0.79140, 0.85379, 0.83394  &
        ,0.81771, 0.89904, 0.93126, 0.95330, 0.96505/

        DATA(gasym(2,8,ib,8),ib=1,17)  & ! Band 8, RH = 87%
        /0.00051, 0.00133, 0.00343, 0.00873, 0.02084, 0.05255  &
        ,0.13334, 0.35372, 0.65797, 0.79588, 0.85494, 0.83181  &
        ,0.82321, 0.90209, 0.93253, 0.95442, 0.96558/

        DATA(gasym(2,8,ib,9),ib=1,17)  & ! Band 8, RH = 88%
        /0.00053, 0.00139, 0.00357, 0.00910, 0.02171, 0.05479  &
        ,0.13921, 0.36977, 0.66473, 0.80066, 0.85614, 0.82981  &
        ,0.82836, 0.90478, 0.93393, 0.95555, 0.96613/

        DATA(gasym(2,8,ib,10),ib=1,17)  & ! Band 8, RH = 89%
        /0.00055, 0.00145, 0.00372, 0.00951, 0.02272, 0.05737  &
        ,0.14600, 0.38809, 0.67226, 0.80574, 0.85762, 0.82801  &
        ,0.83483, 0.90789, 0.93556, 0.95663, 0.96668/

        DATA(gasym(2,8,ib,11),ib=1,17)  & ! Band 8, RH = 90%
        /0.00058, 0.00151, 0.00391, 0.01000, 0.02390, 0.06038  &
        ,0.15397, 0.40914, 0.68090, 0.81113, 0.85946, 0.82520  &
        ,0.84143, 0.91059, 0.93750, 0.95768, 0.96724/

        DATA(gasym(2,8,ib,12),ib=1,17)  & ! Band 8, RH = 91%
        /0.00061, 0.00159, 0.00413, 0.01057, 0.02530, 0.06395  &
        ,0.16349, 0.43345, 0.69111, 0.81690, 0.86129, 0.82161  &
        ,0.84824, 0.91334, 0.93976, 0.95874, 0.96783/

        DATA(gasym(2,8,ib,13),ib=1,17)  & ! Band 8, RH = 92%
        /0.00064, 0.00169, 0.00439, 0.01126, 0.02698, 0.06827  &
        ,0.17509, 0.46158, 0.70337, 0.82320, 0.86247, 0.81832  &
        ,0.85560, 0.91596, 0.94228, 0.95995, 0.96841/

        DATA(gasym(2,8,ib,14),ib=1,17)  & ! Band 8, RH = 93%
        /0.00068, 0.00181, 0.00471, 0.01211, 0.02905, 0.07363  &
        ,0.18956, 0.49400, 0.71804, 0.83009, 0.86306, 0.81357  &
        ,0.86282, 0.91842, 0.94481, 0.96132, 0.96900/

        DATA(gasym(2,8,ib,15),ib=1,17)  & ! Band 8, RH = 94%
        /0.00073, 0.00195, 0.00511, 0.01319, 0.03169, 0.08045  &
        ,0.20822, 0.53074, 0.73500, 0.83727, 0.86354, 0.81046  &
        ,0.87059, 0.92103, 0.94724, 0.96273, 0.96963/

        DATA(gasym(2,8,ib,16),ib=1,17)  & ! Band 8, RH = 95%
        /0.00079, 0.00214, 0.00564, 0.01461, 0.03518, 0.08953  &
        ,0.23336, 0.57082, 0.75331, 0.84405, 0.86235, 0.80778  &
        ,0.87889, 0.92489, 0.94989, 0.96424, 0.97027/

        DATA(gasym(2,8,ib,17),ib=1,17)  & ! Band 8, RH = 96%
        /0.00088, 0.00240, 0.00637, 0.01658, 0.04007, 0.10235  &
        ,0.26934, 0.61137, 0.77194, 0.85153, 0.86012, 0.81003  &
        ,0.88820, 0.93077, 0.95327, 0.96583, 0.97092/

        DATA(gasym(2,8,ib,18),ib=1,17)  & ! Band 8, RH = 97%
        /0.00100, 0.00277, 0.00745, 0.01956, 0.04752, 0.12218  &
        ,0.32554, 0.64888, 0.79263, 0.86063, 0.85335, 0.82235  &
        ,0.90136, 0.93724, 0.95670, 0.96759, 0.97157/

        DATA(gasym(2,8,ib,19),ib=1,17)  & ! Band 8, RH = 98%
        /0.00119, 0.00339, 0.00930, 0.02474, 0.06069, 0.15804  &
        ,0.42407, 0.69110, 0.81948, 0.86729, 0.83686, 0.85109  &
        ,0.91734, 0.94330, 0.96135, 0.96946, 0.97224/

        DATA(gasym(2,8,ib,20),ib=1,17)  & ! Band 8, RH = 99%
        /0.00154, 0.00466, 0.01338, 0.03675, 0.09234, 0.24911  &
        ,0.59662, 0.76793, 0.85173, 0.86554, 0.81500, 0.88824  &
        ,0.93140, 0.95415, 0.96660, 0.97150, 0.97288/

! gasym for Mineral Dust R(imaginary) = 0.0015 (LOW absorption)(No RH effect)
        DATA(gasym_dust(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00107,  0.00268,  0.00676,  0.01699,  0.04263,  0.10735  &
        , 0.28213,  0.61701,  0.76023,  0.82950,  0.82749,  0.83301  &
        , 0.85459,  0.86072,  0.89047,  0.91127,  0.92998/ 
 
        DATA(gasym_dust(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00490,  0.01231,  0.03085,  0.07714,  0.19799,  0.50394  &
        , 0.66481,  0.72768,  0.67004,  0.69465,  0.73482,  0.78342  &
        , 0.80902,  0.83487,  0.85836,  0.88236,  0.90635/ 
 
        DATA(gasym_dust(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.02210,  0.05517,  0.13807,  0.34799,  0.58337,  0.69190  &
        , 0.70762,  0.60868,  0.70584,  0.75077,  0.78408,  0.81217  &
        , 0.83760,  0.86277,  0.88741,  0.91151,  0.93118/ 
 
        DATA(gasym_dust(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00001,  0.00002,  0.00005,  0.00011,  0.00029,  0.00073  &
        , 0.00183,  0.00461,  0.01159,  0.02909,  0.07352,  0.19723  &
        , 0.49168,  0.66569,  0.75801,  0.84784,  0.87389/ 
 
        DATA(gasym_dust(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00002,  0.00004,  0.00010,  0.00026,  0.00065,  0.00165  &
        , 0.00415,  0.01043,  0.02616,  0.06562,  0.16892,  0.46511  &
        , 0.67484,  0.78564,  0.77736,  0.89280,  0.91936/ 
 
        DATA(gasym_dust(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00004,  0.00011,  0.00027,  0.00067,  0.00169,  0.00426  &
        , 0.01070,  0.02690,  0.06781,  0.17736,  0.50452,  0.71737  &
        , 0.82994,  0.88285,  0.92003,  0.94431,  0.95420/ 
 
        DATA(gasym_dust(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00004,  0.00010,  0.00025,  0.00064,  0.00160,  0.00403  &
        , 0.01014,  0.02545,  0.06394,  0.16706,  0.49136,  0.67547  &
        , 0.75696,  0.82304,  0.89691,  0.91021,  0.91682/ 
 
        DATA(gasym_dust(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00011,  0.00027,  0.00069,  0.00173,  0.00436,  0.01097  &
        , 0.02761,  0.06953,  0.17772,  0.45807,  0.71820,  0.84603  &
        , 0.90347,  0.91302,  0.90597,  0.94409,  0.96733/ 

! gasym for Mineral Dust R(imaginary) = 0.003 (MID absorption)(No RH effect)
        DATA(gasym_dust(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00107,  0.00268,  0.00676,  0.01699,  0.04263,  0.10738  &
        , 0.28234,  0.61728,  0.76091,  0.83063,  0.83010,  0.83856  &
        , 0.86246,  0.86982,  0.90315,  0.92773,  0.94773/ 
 
        DATA(gasym_dust(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00490,  0.01231,  0.03085,  0.07716,  0.19809,  0.50410  &
        , 0.66564,  0.72934,  0.67460,  0.70263,  0.74644,  0.79914  &
        , 0.82973,  0.86076,  0.88912,  0.91493,  0.93512/ 
 
        DATA(gasym_dust(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.02210,  0.05517,  0.13810,  0.34806,  0.58381,  0.69307  &
        , 0.71020,  0.61538,  0.71611,  0.76483,  0.80281,  0.83625  &
        , 0.86712,  0.89610,  0.92049,  0.93780,  0.94630/ 
 
        DATA(gasym_dust(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00001,  0.00002,  0.00005,  0.00011,  0.00029,  0.00073  &
        , 0.00183,  0.00461,  0.01159,  0.02909,  0.07352,  0.19723  &
        , 0.49168,  0.66569,  0.75801,  0.84784,  0.87389/ 
 
        DATA(gasym_dust(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00002,  0.00004,  0.00010,  0.00026,  0.00065,  0.00165  &
        , 0.00415,  0.01043,  0.02616,  0.06562,  0.16892,  0.46511  &
        , 0.67484,  0.78564,  0.77736,  0.89280,  0.91936/ 
 
        DATA(gasym_dust(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00004,  0.00011,  0.00027,  0.00067,  0.00169,  0.00426  &
        , 0.01070,  0.02690,  0.06781,  0.17736,  0.50452,  0.71737  &
        , 0.82994,  0.88285,  0.92003,  0.94431,  0.95420/ 
 
        DATA(gasym_dust(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00004,  0.00010,  0.00025,  0.00064,  0.00160,  0.00403  &
        , 0.01014,  0.02545,  0.06394,  0.16706,  0.49136,  0.67547  &
        , 0.75696,  0.82304,  0.89691,  0.91021,  0.91682/ 
 
        DATA(gasym_dust(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00011,  0.00027,  0.00069,  0.00173,  0.00436,  0.01097  &
        , 0.02761,  0.06953,  0.17772,  0.45807,  0.71820,  0.84603  &
        , 0.90347,  0.91302,  0.90597,  0.94409,  0.96733/ 

! gasym for Mineral Dust R(imaginary) = 0.008 (HIGH absorption)(No RH effect)
        DATA(gasym_dust(3,1,ib),ib=1,17)  & ! Band 1, All RH's
        /0.00107, 0.00268, 0.00676, 0.01699, 0.04264, 0.10747  &
        ,0.28301, 0.61817, 0.76309, 0.83417, 0.83819, 0.85421  &
        ,0.88159, 0.89181, 0.92795, 0.95285, 0.96718/

        DATA(gasym_dust(3,2,ib),ib=1,17)  & ! Band 2, All RH's
        /0.00490, 0.01231, 0.03085, 0.07719, 0.19842, 0.50460  &
        ,0.66831, 0.73464, 0.68896, 0.72559, 0.77863, 0.83901  &
        ,0.87691, 0.91052, 0.93413, 0.94697, 0.95155/

        DATA(gasym_dust(3,3,ib),ib=1,17)  & ! Band 3, All RH's
        /0.02390, 0.05985, 0.15098, 0.36031, 0.58493, 0.69319  &
        ,0.71724, 0.63499, 0.74390, 0.80146, 0.84751, 0.88670  &
        ,0.91730, 0.93703, 0.94596, 0.94845, 0.94879/

        DATA(gasym_dust(3,4,ib),ib=1,17)  & ! Band 4, All RH's
        /0.00001, 0.00002, 0.00005, 0.00011, 0.00029, 0.00073  &
        ,0.00183, 0.00461, 0.01159, 0.02909, 0.07352, 0.19723  &
        ,0.49168, 0.66569, 0.75801, 0.84784, 0.87389/

        DATA(gasym_dust(3,5,ib),ib=1,17)  & ! Band 5, All RH's
        /0.00002, 0.00004, 0.00010, 0.00026, 0.00065, 0.00165  &
        ,0.00415, 0.01043, 0.02617, 0.06562, 0.16892, 0.46511  &
        ,0.67484, 0.78564, 0.77736, 0.89280, 0.91936/

        DATA(gasym_dust(3,6,ib),ib=1,17)  & ! Band 6, All RH's
        /0.00004, 0.00011, 0.00026, 0.00067, 0.00169, 0.00426  &
        ,0.01071, 0.02690, 0.06781, 0.17736, 0.50452, 0.71737  &
        ,0.82994, 0.88285, 0.92003, 0.94431, 0.95420/

        DATA(gasym_dust(3,7,ib),ib=1,17)  & ! Band 7, All RH's
        /0.00004, 0.00010, 0.00025, 0.00064, 0.00160, 0.00403  &
        ,0.01014, 0.02545, 0.06394, 0.16706, 0.49135, 0.67547  &
        ,0.75696, 0.82304, 0.89691, 0.91021, 0.91682/

        DATA(gasym_dust(3,8,ib),ib=1,17)  & ! Band 8, All RH's
        /0.00011, 0.00027, 0.00069, 0.00173, 0.00436, 0.01087  &
        ,0.02761, 0.06953, 0.17772, 0.45807, 0.71820, 0.84603  &
        ,0.90347, 0.91302, 0.90597, 0.94409, 0.96733/

! gasym for Absorbing Carbon (1% BC, 99% OC)(No RH effect)
        DATA(gasym_carb(1,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00107,  0.00270,  0.00679,  0.01705,  0.04269,  0.10730  &
        , 0.28527,  0.58816,  0.71748,  0.75661,  0.70471,  0.79154  &
        , 0.84682,  0.90083,  0.92901,  0.94486,  0.95114/

        DATA(gasym_carb(1,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00487,  0.01223,  0.03064,  0.07669,  0.19739,  0.50061  &
        , 0.67041,  0.73939,  0.70207,  0.74514,  0.84526,  0.89083  &
        , 0.92559,  0.94233,  0.94802,  0.94930,  0.94959/

        DATA(gasym_carb(1,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.02215,  0.05533,  0.13872,  0.34940,  0.58663,  0.69946  &
        , 0.72426,  0.65395,  0.77353,  0.83330,  0.88199,  0.91698  &
        , 0.93769,  0.94623,  0.94834,  0.94864,  0.94867/

        DATA(gasym_carb(1,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00001,  0.00002,  0.00005,  0.00013,  0.00032,  0.00079  &
        , 0.00200,  0.00503,  0.01266,  0.03182,  0.08080,  0.22206  &
        , 0.48738,  0.64771,  0.74894,  0.84326,  0.87006/

        DATA(gasym_carb(1,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00002,  0.00004,  0.00011,  0.00027,  0.00069,  0.00173  &
        , 0.00436,  0.01095,  0.02747,  0.06899,  0.17976,  0.49530  &
        , 0.67468,  0.76643,  0.81133,  0.89678,  0.91569/

        DATA(gasym_carb(1,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00006,  0.00016,  0.00039,  0.00099,  0.00248,  0.00625  &
        , 0.01573,  0.03964,  0.10192,  0.29703,  0.46423,  0.63679  &
        , 0.77588,  0.82940,  0.84988,  0.86036,  0.86578/

        DATA(gasym_carb(1,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00005,  0.00012,  0.00031,  0.00079,  0.00200,  0.00503  &
        , 0.01266,  0.03187,  0.08120,  0.22619,  0.48089,  0.62299  &
        , 0.68746,  0.81476,  0.85347,  0.87638,  0.88352/

        DATA(gasym_carb(1,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00012,  0.00030,  0.00074,  0.00187,  0.00471,  0.01185  &
        , 0.02972,  0.07450,  0.19173,  0.50983,  0.70185,  0.79465  &
        , 0.78166,  0.80473,  0.89598,  0.93090,  0.94976/

! gasym for Absorbing Carbon (2% BC, 98% OC)(No RH effect)
        DATA(gasym_carb(2,1,ib),ib=1,17)  & ! Band 1, All RH's
        / 0.00107,  0.00270,  0.00680,  0.01708,  0.04279,  0.10763  &
        , 0.28695,  0.58914,  0.72051,  0.76198,  0.72395,  0.81290  &
        , 0.87119,  0.91879,  0.93973,  0.94903,  0.95197/

        DATA(gasym_carb(2,2,ib),ib=1,17)  & ! Band 2, All RH's
        / 0.00487,  0.01225,  0.03070,  0.07689,  0.19829,  0.50174  &
        , 0.67298,  0.74409,  0.71617,  0.76916,  0.86766,  0.90935  &
        , 0.93573,  0.94547,  0.94811,  0.94879,  0.94904/

        DATA(gasym_carb(2,3,ib),ib=1,17)  & ! Band 3, All RH's
        / 0.02220,  0.05548,  0.13924,  0.35063,  0.58838,  0.70320  &
        , 0.73248,  0.67786,  0.80437,  0.86426,  0.90901,  0.93414  &
        , 0.94460,  0.94740,  0.94797,  0.94812,  0.94814/

        DATA(gasym_carb(2,4,ib),ib=1,17)  & ! Band 4, All RH's
        / 0.00001,  0.00002,  0.00005,  0.00013,  0.00032,  0.00079  &
        , 0.00200,  0.00503,  0.01266,  0.03181,  0.08078,  0.22199  &
        , 0.48671,  0.64839,  0.75039,  0.84308,  0.86982/

        DATA(gasym_carb(2,5,ib),ib=1,17)  & ! Band 5, All RH's
        / 0.00002,  0.00004,  0.00011,  0.00027,  0.00069,  0.00173  &
        , 0.00436,  0.01097,  0.02751,  0.06912,  0.18033,  0.49576  &
        , 0.67629,  0.76929,  0.81859,  0.89819,  0.91594/

        DATA(gasym_carb(2,6,ib),ib=1,17)  & ! Band 6, All RH's
        / 0.00006,  0.00015,  0.00039,  0.00098,  0.00248,  0.00624  &
        , 0.01570,  0.03955,  0.10168,  0.29595,  0.46527,  0.63866  &
        , 0.77685,  0.82951,  0.84996,  0.86044,  0.86586/

        DATA(gasym_carb(2,7,ib),ib=1,17)  & ! Band 7, All RH's
        / 0.00005,  0.00012,  0.00031,  0.00079,  0.00199,  0.00502  &
        , 0.01264,  0.03182,  0.08107,  0.22568,  0.48122,  0.62597  &
        , 0.69327,  0.81826,  0.85513,  0.87658,  0.88349/

        DATA(gasym_carb(2,8,ib),ib=1,17)  & ! Band 8, All RH's
        / 0.00012,  0.00030,  0.00074,  0.00188,  0.00472,  0.01187  &
        , 0.02978,  0.07468,  0.19258,  0.51186,  0.70457,  0.79784  &
        , 0.78999,  0.82788,  0.90998,  0.93963,  0.95278/

! Returning value for gasym:
        if (aerotype <= 2) value = gasym(aerotype,radband,bin,rh)
        if (aerotype == 3) value = gasym_dust(dust_ref_im,radband,bin)
        if (aerotype == 4) value = gasym_carb(1,radband,bin)
        if (aerotype == 5) value = gasym_carb(2,radband,bin)

return
END SUBROUTINE aerogasym
