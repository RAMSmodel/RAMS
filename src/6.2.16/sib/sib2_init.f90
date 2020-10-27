!##############################################################################
Subroutine init_sib1 (n1,n2,n3,rco2p,ifm)

use mem_sib
use mem_grid, only:zt,print_msg

implicit none

integer :: n1,n2,n3,i,j,k,ip,ifm
real, dimension(n1,n2,n3) :: rco2p
real :: f

!This is the conversion factor to translate CO2 from (ppm) to
!a CO2 mixing ratio that RAMS can advect as a conserved quantity.
!CO2 molar mass = 44g/mol, AIR molar mass = 29g/mol
!Multiply CO2 ppm by 44/29/1.e6 to get the mass fraction which is 
!the same as mixing ratio units need by RAMS for conservative advection.
data f /1.51724e-6/  ! 1.51724e-6= 44/29/1.e6

!Initialization notes:
! 1. SiB needs initial stomatal resistance (rst). LEAF3 provides this.
! 2. SiB/RAMS needs initial CO2 profile
!     (RAMS-rco2p in mass fraction, SiB-atm_co2_ppm in ppm)

!Set CO2 initial profile
co2_init(1) = max(0., co2_init(1))
do k = 2,n1
    if(co2_init(k) <= 0.) co2_init(k) = co2_init(k-1)
enddo

if(print_msg)then
 do k = 1,n1
  print*,'CO2-init on grid:',ifm,k,zt(k),co2_init(k)
 enddo
endif

!Store CO2 in sib tracer in mixing ratio units rather than ppm.
do j = 1,n3
 do i = 1,n2
  do k = 2,n1
    rco2p(k,i,j) = f * co2_init(k-1)
  enddo
  rco2p(1,i,j) = rco2p(2,i,j) !Copy model level 2 data to level 1.
 enddo
enddo

return
END SUBROUTINE init_sib1

!##############################################################################
Subroutine init_sib2 (n1,n2,n3,np,rco2p,pco2ap,pi0,pp)

use mem_sib
use rconstants

implicit none

integer :: n1,n2,n3,np,i,j,k,ip
real, dimension(n1,n2,n3) :: pi0,pp,rco2p
real, dimension(n2,n3,np) :: pco2ap
real :: pis2,hcpi,prss,f

!This is the conversion factor to translate CO2 from (ppm) to
!a CO2 mixing ratio that RAMS can advect as a conserved quantity.
!CO2 molar mass = 44g/mol, AIR molar mass = 29g/mol
!Multiply CO2 ppm by 44/29/1.e6 to get the mass fraction which is 
!the same as mixing ratio units need by RAMS for conservative advection.
data f /1.51724e-6/  ! 1.51724e-6= 44/29/1.e6

!Initialize pco2ap (canopy air space pCO2, Pascals)

!Ppm in a gas is normally expressed on a mole fraction basis, so 385
!ppm of CO2 (in the atmosphere) is also 385 µmol/mol. By Dalton's law
!of partial pressure, the partial pressure of CO2 to total pressure
!will be in the same ratio. "Real" atmospheric pressure varies with
!altitude above sea level, and day-to-day with the weather; however,
!"standard" atmospheric pressure is 101325 Pa. 
!For example, the partial pressure of CO2 is:
!385 µmol/mol x 1.e-6(unit conversion to get mol/mol) x 101325 Pa = 39 Pa, 
!(at least on a day when real and standard pressure are the same value).

hcpi = .5 * cpi
do j = 1,n3
 do i = 1,n2
   pis2 = (pp(1,i,j) + pi0(1,i,j) + pp(2,i,j) + pi0(2,i,j)) * hcpi
   prss = pis2 ** cpor * p00 !air pressure (Pa)
   do ip = 2,np !Only done for land patches
    pco2ap(i,j,ip) = rco2p(1,i,j) / f * 1.e-6 * prss
   enddo
 enddo
enddo

return
END SUBROUTINE init_sib2
