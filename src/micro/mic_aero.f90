!##############################################################################
Subroutine aerosols ()

use mem_basic
use mem_micro
use mem_grid
use mem_leaf
use node_mod
use micphys

implicit none

integer :: i,j

!Run the SEASALT and DUST Source model before the call to Micro
if(idust==2) then
  CALL dust_sources (mzp,mxp,myp,ia,iz,ja,jz      &
                    ,grid_g(ngrid)%rtgt           &
                    ,grid_g(ngrid)%glat           &
                    ,grid_g(ngrid)%glon           &
                    ,basic_g(ngrid)%up            &
                    ,basic_g(ngrid)%vp            &
                    ,basic_g(ngrid)%dn0           &
                    ,micro_g(ngrid)%md1np         &
                    ,micro_g(ngrid)%md2np         &
                    ,micro_g(ngrid)%md1mp         &
                    ,micro_g(ngrid)%md2mp         &
                    ,micro_g(ngrid)%dustfrac      &
                    ,leaf_g(ngrid)%soil_water     &
                    ,leaf_g(ngrid)%patch_area     &
                    ,leaf_g(ngrid)%leaf_class     &
                    ,leaf_g(ngrid)%soil_text      &
                    ,leaf_g(ngrid)%veg_rough      &
                    )
endif

if(isalt==2) then
  CALL salt_sources (mzp,mxp,myp,ia,iz,ja,jz      &
                    ,grid_g(ngrid)%rtgt           &
                    ,basic_g(ngrid)%up            &
                    ,basic_g(ngrid)%vp            &
                    ,basic_g(ngrid)%dn0           &
                    ,micro_g(ngrid)%salt_film_np  &
                    ,micro_g(ngrid)%salt_jet_np   &
                    ,micro_g(ngrid)%salt_spum_np  &
                    ,micro_g(ngrid)%salt_film_mp  &
                    ,micro_g(ngrid)%salt_jet_mp   &
                    ,micro_g(ngrid)%salt_spum_mp  &
                    ,leaf_g(ngrid)%patch_area     &
                    ,leaf_g(ngrid)%leaf_class     &
                    )
endif

! Aerosol dry and wet deposition call when micro LEVEL < 3
if(iaerodep==1 .and. level<3) then
 do j = ja,jz
  do i = ia,iz

   CALL aero_copy (1,mzp &
    ,micro_g(ngrid)%cn1np(1,i,j),micro_g(ngrid)%cn1mp(1,i,j) &
    ,micro_g(ngrid)%cn2np(1,i,j),micro_g(ngrid)%cn2mp(1,i,j) &
    ,micro_g(ngrid)%md1np(1,i,j),micro_g(ngrid)%md1mp(1,i,j) &
    ,micro_g(ngrid)%md2np(1,i,j),micro_g(ngrid)%md2mp(1,i,j) &
    ,micro_g(ngrid)%salt_film_np(1,i,j),micro_g(ngrid)%salt_film_mp(1,i,j) &
    ,micro_g(ngrid)%salt_jet_np(1,i,j) ,micro_g(ngrid)%salt_jet_mp(1,i,j)  &
    ,micro_g(ngrid)%salt_spum_np(1,i,j),micro_g(ngrid)%salt_spum_mp(1,i,j) &
    ,micro_g(ngrid)%abc1np(1,i,j),micro_g(ngrid)%abc1mp(1,i,j) &
    ,micro_g(ngrid)%abc2np(1,i,j),micro_g(ngrid)%abc2mp(1,i,j))

   CALL deposition_driver (i,j,mzp,zm         &
    ,grid_g(ngrid)%rtgt(i,j)                  &
    ,basic_g(ngrid)%rv(1,i,j)                 &
    ,basic_g(ngrid)%theta(1,i,j)              &
    ,basic_g(ngrid)%up(1,i,j)                 &
    ,basic_g(ngrid)%vp(1,i,j)                 &
    ,basic_g(ngrid)%dn0(1,i,j)                &
    ,basic_g(ngrid)%pi0(1,i,j)                &
    ,basic_g(ngrid)%pp(1,i,j)                 &
    ,leaf_g(ngrid)%leaf_class(i,j,1:npatch)   &
    ,leaf_g(ngrid)%patch_area(i,j,1:npatch)   &
    ,leaf_g(ngrid)%ustar(i,j,1:npatch)        &
    ,leaf_g(ngrid)%patch_rought(i,j,1:npatch)  &
    ,imonth1                                  &
    )

   CALL aero_copy (2,mzp &
    ,micro_g(ngrid)%cn1np(1,i,j),micro_g(ngrid)%cn1mp(1,i,j) &
    ,micro_g(ngrid)%cn2np(1,i,j),micro_g(ngrid)%cn2mp(1,i,j) &
    ,micro_g(ngrid)%md1np(1,i,j),micro_g(ngrid)%md1mp(1,i,j) &
    ,micro_g(ngrid)%md2np(1,i,j),micro_g(ngrid)%md2mp(1,i,j) &
    ,micro_g(ngrid)%salt_film_np(1,i,j),micro_g(ngrid)%salt_film_mp(1,i,j) &
    ,micro_g(ngrid)%salt_jet_np(1,i,j) ,micro_g(ngrid)%salt_jet_mp(1,i,j)  &
    ,micro_g(ngrid)%salt_spum_np(1,i,j),micro_g(ngrid)%salt_spum_mp(1,i,j) &
    ,micro_g(ngrid)%abc1np(1,i,j),micro_g(ngrid)%abc1mp(1,i,j) &
    ,micro_g(ngrid)%abc2np(1,i,j),micro_g(ngrid)%abc2mp(1,i,j))

  enddo
 enddo
endif

return
END SUBROUTINE aerosols

!##############################################################################
Subroutine aerosol_init ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

real :: weightfac

! Set aerosol density depending on chemistry and soluble fraction
! Pure quantity densities (kg/m3) are:
! NH42S04 = 1769. (ammonium sulfate)
! Clay Dust (smaller) = 2500.
! Silt Dust (larger) = 2650.
! NaCl = 2165. (sodium chloride)

! Set Aerosol density (kg/m3) based on weighted mixture of soluble
! and insoluble material. Assume insoluble core to be like that of 
! silt dust with density = 2650 kg/m3, except for acat=3 which is
! already set to small sized clay dust
! Also set vanthoff factors for given chemistry

if(iprntstmt>=1 .and. print_msg) print*,''
if(iprntstmt>=1 .and. print_msg) print*,'Setting up default aerosol densities:'

do acat=1,aerocat
 
 aero_rhosol(acat)   = 0.0
 aero_vanthoff(acat) = 0.0

 if(acat==3) then !if small dust (clay)
   weightfac = 2500. * (1.0-aero_epsilon(acat)) !clay dust core
 else !all other
   weightfac = 2650. * (1.0-aero_epsilon(acat)) !silt dust core
 endif

 if(iaero_chem(acat)==1) then !NH42S04
   aero_rhosol(acat) = 1769. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 3
 elseif(iaero_chem(acat)==2) then !NaCl
   aero_rhosol(acat) = 2165. * aero_epsilon(acat) + weightfac
   aero_vanthoff(acat) = 2
 endif

 if(iprntstmt>=1 .and. print_msg) print*,'acat,rg,rho,i:',acat &
     ,aero_medrad(acat),aero_rhosol(acat),aero_vanthoff(acat)

enddo

if(iprntstmt>=1 .and. print_msg) print*,''

return
END SUBROUTINE aerosol_init

!##############################################################################
Subroutine aero_copy (aflag,m1,cn1np,cn1mp,cn2np,cn2mp,md1np,md1mp &
                    ,md2np,md2mp,salt_film_np,salt_film_mp,salt_jet_np &
                    ,salt_jet_mp,salt_spum_np,salt_spum_mp &
                    ,abc1np,abc1mp,abc2np,abc2mp)

!This routine is called in the event that MICRO LEVEL=1,2 so that
!aerosols can still be allowed to impact radiation.

use micphys

implicit none

integer :: m1,k,aflag
real, dimension(m1) :: cn1np,cn1mp,cn2np,cn2mp,md1np,md1mp &
                    ,md2np,md2mp,salt_film_np,salt_film_mp,salt_jet_np &
                    ,salt_jet_mp,salt_spum_np,salt_spum_mp &
                    ,abc1np,abc1mp,abc2np,abc2mp
if(aflag==1)then
 !Zero out aerosol scratch arrays
 do acat = 1,aerocat
     do k = 1,m1
       aerocon(k,acat) = 0.0
       aeromas(k,acat) = 0.0
     enddo
 enddo
 !Fill scratch arrays for aerosol modes for level=1,2
 do k = 1,m1-1
   if (iaerosol > 0) then
     aerocon(k,1) = cn1np(k)
     aeromas(k,1) = cn1mp(k)
     aerocon(k,2) = cn2np(k)
     aeromas(k,2) = cn2mp(k)
   endif
   if (idust > 0) then
     aerocon(k,3) = md1np(k)
     aeromas(k,3) = md1mp(k)
     aerocon(k,4) = md2np(k)
     aeromas(k,4) = md2mp(k)
   endif
   if (isalt > 0) then
     aerocon(k,5) = salt_film_np(k)
     aeromas(k,5) = salt_film_mp(k)
     aerocon(k,6) = salt_jet_np(k)
     aeromas(k,6) = salt_jet_mp(k)
     aerocon(k,7) = salt_spum_np(k)
     aeromas(k,7) = salt_spum_mp(k)
   endif
   if (iabcarb > 0) then
     aerocon(k,8) = abc1np(k)
     aeromas(k,8) = abc1mp(k)
     aerocon(k,9) = abc2np(k)
     aeromas(k,9) = abc2mp(k)
   endif
 enddo

elseif(aflag==2)then
 !Copy back scratch arrays to aerosol modes for level=1,2
 do k = 1,m1-1
   if (iaerosol > 0) then
    cn1np(k) = aerocon(k,1)
    cn1mp(k) = aeromas(k,1)
    cn2np(k) = aerocon(k,2)
    cn2mp(k) = aeromas(k,2)
   endif
   if (idust > 0) then
    md1np(k) = aerocon(k,3)
    md1mp(k) = aeromas(k,3)
    md2np(k) = aerocon(k,4)
    md2mp(k) = aeromas(k,4)
   endif
   if (isalt > 0) then
    salt_film_np(k) = aerocon(k,5)
    salt_film_mp(k) = aeromas(k,5)
    salt_jet_np(k)  = aerocon(k,6)
    salt_jet_mp(k)  = aeromas(k,6)
    salt_spum_np(k) = aerocon(k,7)
    salt_spum_mp(k) = aeromas(k,7)
   endif
   if (iabcarb > 0) then
    abc1np(k) = aerocon(k,8)
    abc1mp(k) = aeromas(k,8)
    abc2np(k) = aerocon(k,9)
    abc2mp(k) = aeromas(k,9)
   endif
 enddo
endif

return
END SUBROUTINE aero_copy
