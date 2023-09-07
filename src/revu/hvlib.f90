!##############################################################################
Subroutine rams_varlib (cvar,n1,n2,n3,ngrd,a,b,flnm,cdname,cdunits)

use micphys
use rcommons
use mem_grid

implicit none

integer :: n1,n2,n3,ngrd,lv,idim_type,irecind,irecsize,irecsizep,ierr,kp
character(len=*) :: cvar,flnm,cdname,cdunits
real :: a(*),b(*)
real, allocatable, save :: c(:),d(:),e(:),f(:)
integer, external :: rams_getvar
integer, external :: lastchar
real, allocatable, dimension(:,:,:) :: pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8,pv9
real, allocatable, dimension(:,:) :: dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9

!Arrays for computing reflectivity
!tmp arrays for mixing ratios
real, allocatable, dimension(:,:,:) :: trmix,tgmix,thmix,tpmix,tsmix,tamix
!tmp arrays for # conc's
real, allocatable, dimension(:,:,:) :: trnt,tgnt,thnt,tpnt,tsnt,tant
real, allocatable, dimension(:,:,:) :: reflc  !array for rflctvty values
real, allocatable, dimension(:,:,:) :: tden   !array for grd pt density

if(maxmem > 0) then
   if (allocated(c)) deallocate (c); allocate(c(maxmem))
   if (allocated(d)) deallocate (d); allocate(d(maxmem))
   if (allocated(e)) deallocate (e); allocate(e(maxmem))
   if (allocated(f)) deallocate (f); allocate(f(maxmem))
endif
                       
lv=lastchar(cvar)
!print*,'===> varlib- ',cvar,n1,n2,n3,ngrd

ivar_type=0
ierr_getvar=0
ierr=0
ifound=0
iany=0
iptrans_lim=0

!###########################################################################
! EMPTY PLACEHOLDER VARIABLE WITH ZEROS - 1
!###########################################################################
if(cvar(1:lv).eq.'empty3d') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   CALL rams_comp_empty (n1,n2,n3,a)
   cdname='empty-3D-variable;'
   cdunits='none;'

!######################################################################
! AOD DATA FROM OFFLINE CODE
!######################################################################
elseif(cvar(1:lv).eq.'ccn_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('ccn_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='ccn_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'ccn_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('ccn_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='ccn_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'dust1_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('dust1_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='dust1_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'dust1_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('dust1_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='dust1_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'dust2_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('dust2_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='dust2_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'dust2_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('dust2_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='dust2_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'regen_aero1_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('regen_aero1_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero1_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'regen_aero1_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('regen_aero1_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero1_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'regen_aero2_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('regen_aero2_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero2_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'regen_aero2_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('regen_aero2_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero2_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_film_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_film_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_film_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_film_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_film_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_film_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_jet_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_jet_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_jet_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_jet_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_jet_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_jet_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_spume_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_spume_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_spume_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'salt_spume_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('salt_spume_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='salt_spume_wet_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'Total_dry_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('Total_dry_AOD_550',idim_type,ngrd,a,flnm)
   cdname='Total_dry_AOD_550;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'Total_wet_AOD_550') then
   ivar_type=2
   iptrans_lim=1
   ierr=rams_getvar('Total_wet_AOD_550',idim_type,ngrd,a,flnm)
   cdname='Total_wet_AOD_550;'
   cdunits='AOD;'

!######################################################################
! EXTINCTION COEFFICIENT DATA FROM OFFLINE CODE
!######################################################################
elseif(cvar(1:lv).eq.'ccn_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('ccn_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='ccn_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'ccn_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('ccn_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='ccn_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'dust1_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('dust1_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='dust1_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'dust1_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('dust1_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='dust1_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'dust2_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('dust2_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='dust2_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'dust2_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('dust2_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='dust2_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'regen_aero1_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('regen_aero1_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero1_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'regen_aero1_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('regen_aero1_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero1_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'regen_aero2_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('regen_aero2_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero2_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'regen_aero2_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('regen_aero2_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='regen_aero2_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_film_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_film_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_film_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_film_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_film_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_film_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_jet_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_jet_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_jet_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_jet_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_jet_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_jet_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_spume_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_spume_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_spume_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'salt_spume_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('salt_spume_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='salt_spume_wet_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'Total_dry_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('Total_dry_ext_550',idim_type,ngrd,a,flnm)
   cdname='Total_dry_ext_550;'
   cdunits='1/Mm;'

elseif(cvar(1:lv).eq.'Total_wet_ext_550') then
   ivar_type=3
   iptrans_lim=1
   ierr=rams_getvar('Total_wet_ext_550',idim_type,ngrd,a,flnm)
   cdname='Total_wet_ext_550;'
   cdunits='1/Mm;'

!###########################################################################
! 3D VELOCITY AND VORTICITY VARIABLES
!###########################################################################
elseif(cvar(1:lv).eq.'u') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   cdname='u;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'v') then
   ivar_type=3
   ierr=rams_getvar('VP',idim_type,ngrd,a,flnm)
   cdname='v;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'u_avg') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgu (n1,n2,n3,a)
   cdname='u_avg;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'v_avg') then
   ivar_type=3
   ierr=rams_getvar('VP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgv (n1,n2,n3,a)
   cdname='v_avg;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'ue') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   CALL rams_comp_rotate (n1,n2,n3,a,c,ngrd)
   cdname='ue;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'ve') then
   ivar_type=3
   ierr=rams_getvar('VP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('UP',idim_type,ngrd,c,flnm)
   CALL rams_comp_rotate (n1,n2,n3,c,a,ngrd)
   cdname='ve;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'ue_avg') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   CALL rams_comp_rotate (n1,n2,n3,a,c,ngrd)
   CALL rams_comp_avgu (n1,n2,n3,a)
   cdname='ue_avg;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'ve_avg') then
   ivar_type=3
   ierr=rams_getvar('VP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('UP',idim_type,ngrd,c,flnm)
   CALL rams_comp_rotate (n1,n2,n3,c,a,ngrd)
   CALL rams_comp_avgv (n1,n2,n3,a)
   cdname='ve_avg;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'w') then
   ivar_type=3
   ierr=rams_getvar('WP',idim_type,ngrd,a,flnm)
   cdname='w;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'wcms') then
   ivar_type=3
   ierr=rams_getvar('WP',idim_type,ngrd,a,flnm)
   CALL rams_comp_wcms (n1,n2,n3,a)
   cdname='w;'
   cdunits='cm/s;'

elseif(cvar(1:lv).eq.'w_avg') then
   ivar_type=3
   ierr=rams_getvar('WP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgw (n1,n2,n3,a)
   cdname='w_avg;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'speed') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgu (n1,n2,n3,a)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   CALL rams_comp_avgv (n1,n2,n3,c)
   CALL rams_comp_speed (n1,n2,n3,a,c)
   cdname='speed;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'speed_mph') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgu (n1,n2,n3,a)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   CALL rams_comp_avgv (n1,n2,n3,c)
   CALL rams_comp_speed (n1,n2,n3,a,c)
   CALL rams_comp_mults (n1,n2,n3,a,2.237)
   cdname='speed;'
   cdunits='mph;'

elseif(cvar(1:lv).eq.'speed10m') then
   ivar_type=2
   ierr=rams_getvar('UP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,d,flnm)
   CALL rams_comp_speed (n1,n2,n3,c,d)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),npatch))
   ierr=rams_getvar('USTAR',idim_type,ngrd,pv1,flnm)
   ierr=rams_getvar('PATCH_ROUGH',idim_type,ngrd,pv2,flnm)
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,pv3,flnm)
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,pv4,flnm)
   CALL rams_reduced_wind (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),npatch  &
                          ,a,c,pv1,10.,ztn(2,ngrd),pv2,pv4,pv3,d,f,e  &
                          ,zmn(nnzp(1)-1,1))
   deallocate (pv1,pv2,pv3,pv4)
   cdname='speed-10m-AGL;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'direction') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   CALL rams_comp_avgu (n1,n2,n3,a)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   CALL rams_comp_avgv (n1,n2,n3,a)
   CALL rams_comp_dir (n1,n2,n3,a,c,ngrd)
   cdname='direction;'
   cdunits='deg;'

elseif(cvar(1:lv).eq.'relvortx') then
   ivar_type=3
   ierr=rams_getvar('VP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('WP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   CALL rams_comp_relvortx (n1,n2,n3,a,c,b,d,ngrd)
   cdname='x-vorticity;'
   cdunits='rad/s;'

elseif(cvar(1:lv).eq.'relvorty') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('WP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   CALL rams_comp_relvorty (n1,n2,n3,a,c,b,d,ngrd)
   cdname='y-vorticity;'
   cdunits='rad/s;'

elseif(cvar(1:lv).eq.'relvortz') then
   ivar_type=3 
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   CALL rams_comp_relvortz (n1,n2,n3,a,c,b,d,ngrd)
   cdname='relative-z-vorticity;'
   cdunits='rad/s;'

elseif(cvar(1:lv).eq.'absvortz') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   CALL rams_comp_totvortz (n1,n2,n3,a,c,b,d,ngrd)
   cdname='absolute-z-vorticity;'
   cdunits='rad/s;'

elseif(cvar(1:lv).eq.'potvortz') then
   ivar_type=3
   ierr=rams_getvar('UP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   CALL rams_comp_totvortz (n1,n2,n3,a,c,b,d,ngrd)
   CALL rams_comp_dn0 (n1,n2,n3,e,b,c,d,ngrd)
   ierr=rams_getvar('THETA',idim_type,ngrd,b,flnm)
   CALL rams_comp_potvortz (n1,n2,n3,a,b,c,e,d,ngrd)
   CALL rams_comp_mults (n1,n2,n3,a, 9.80 )  
   cdname='potential-z-vorticity;'
   cdunits='rad/s;'

elseif(cvar(1:lv).eq.'horiz_div') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('WP',idim_type,ngrd,a,flnm)
   CALL rams_comp_horizdiv (n1,n2,n3,a,d,ngrd)
   cdname='horizontal-divergence;'
   cdunits='/s;'

!#####################################################################
! 3D THERMODYNAMIC PROPERTIES OF AIR
!#####################################################################
elseif(cvar(1:lv).eq.'pi') then
   ivar_type=3
   ierr=rams_getvar('PI',idim_type,ngrd,a,flnm)
   cdname='Exner-func;'
   cdunits='J/(kg*K);'

elseif(cvar(1:lv).eq.'press') then
   ivar_type=3
   ierr=rams_getvar('PI',idim_type,ngrd,a,flnm)
   CALL rams_comp_press (n1,n2,n3,a)
   cdname='pressure;'
   cdunits='mb;'

elseif(cvar(1:lv).eq.'pprime') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,a,flnm)
   CALL rams_comp_z (n1,n2,n3,c,a,ngrd)

   ierr=rams_getvar('PI',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   CALL rams_comp_slpress (n1,n2,n3,a,d,c,a)

   ierr=rams_getvar('GLAT',idim_type,ngrd,c,flnm)
   CALL rams_comp_pprime (n1,n2,n3,a,c)
   cdname='mslp-perturbation;'
   cdunits='mb;'

elseif(cvar(1:lv).eq.'theta_il') then
   ivar_type=3
   ierr=rams_getvar('THP',idim_type,ngrd,a,flnm)
   cdname='ice-liquid-potential-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'theta') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   cdname='potential-temperature;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'dn0') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,a,e,ngrd)
   cdname='reference-density;'
   cdunits='kg/m3;'

elseif(cvar(1:lv).eq.'pi0') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,a,b,c,e,ngrd)
   cdname='reference-Exner-func;'
   cdunits='J/(kg K);'

elseif(cvar(1:lv).eq.'th0') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,a,c,e,ngrd)
   cdname='reference-virtual-potential-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'pert_pressure') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,a,b,e,ngrd)
   ierr=rams_getvar('PI',idim_type,ngrd,a,flnm)
   if (ierr.eq.0) CALL rams_comp_ppress (n1,n2,n3,a,c)
   cdname='perturbation-pressure;'
   cdunits='mb;'

elseif(cvar(1:lv).eq.'tempk') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,a,c)
   cdname='temperature;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'tempc') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,a,c)
   CALL rams_comp_tempC (n1,n2,n3,a)
   cdname='temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'tempf') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,a,c)
   CALL rams_comp_tempF (n1,n2,n3,a)
   cdname='temperature;'
   cdunits='F;'

elseif(cvar(1:lv).eq.'theta_e') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_thete (n1,n2,n3,a,c,d)
   cdname='equivalent-potential-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'theta_v') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   CALL rams_comp_thetv (n1,n2,n3,a,c)
   cdname='virtual-potential-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'theta_rho') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('RCP',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('RRP',idim_type,ngrd,e,flnm)
   CALL rams_comp_thetrho (n1,n2,n3,a,c,d,e)
   cdname='density-potential-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'buoyancy_liquid') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('RCP',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('RRP',idim_type,ngrd,e,flnm)
   CALL rams_comp_thetrho_buoy (n1,n2,n3,a,c,d,e)
   cdname='buoyancy-liquid;'
   cdunits='m/s2;'

elseif(cvar(1:lv).eq.'tempf2m') then
   ivar_type=2
   ierr=rams_getvar('UP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,d,flnm)
   CALL rams_comp_speed (n1,n2,n3,c,d)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv5(nnxp(ngrd),nnyp(ngrd),npatch))
   ierr=rams_getvar('USTAR',idim_type,ngrd,pv1,flnm)
   ierr=rams_getvar('SOIL_ROUGH',idim_type,ngrd,pv2,flnm)
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,pv3,flnm)
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,pv4,flnm)
   ierr=rams_getvar('TSTAR',idim_type,ngrd,pv5,flnm)
   CALL rams_reduced_temp (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),npatch  &
                          ,a,c,pv1,pv5,2.,ztn(2,ngrd),pv2,pv4,pv3,d,f,e  &
                          ,zmn(nnzp(1)-1,1))
   deallocate (pv1,pv2,pv3,pv4,pv5)
   CALL rams_comp_tempK (n1,n2,1,a,f)
   CALL rams_comp_tempF (n1,n2,1,a)
   cdname='temp-2m-AGL;'
   cdunits='F;'

elseif(cvar(1:lv).eq.'tempc2m') then
   ivar_type=2
   ierr=rams_getvar('UP',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('VP',idim_type,ngrd,d,flnm)
   CALL rams_comp_speed (n1,n2,n3,c,d)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),npatch))
   allocate (pv5(nnxp(ngrd),nnyp(ngrd),npatch))
   ierr=rams_getvar('USTAR',idim_type,ngrd,pv1,flnm)
   ierr=rams_getvar('SOIL_ROUGH',idim_type,ngrd,pv2,flnm)
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,pv3,flnm)
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,pv4,flnm)
   ierr=rams_getvar('TSTAR',idim_type,ngrd,pv5,flnm)
   CALL rams_reduced_temp (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),npatch  &
                          ,a,c,pv1,pv5,2.,ztn(2,ngrd),pv2,pv4,pv3,d,f,e  &
                          ,zmn(nnzp(1)-1,1))
   deallocate (pv1,pv2,pv3,pv4,pv5)
   CALL rams_comp_tempK (n1,n2,1,a,f)
   CALL rams_comp_tempC (n1,n2,1,a)
   cdname='temp-2m-AGL;'
   cdunits='C;'

!#####################################################################
!3D HYDROMETEOR GAMMA DISTRIBUTION INFO
!#####################################################################
elseif(cvar(1:lv).eq.'cloud_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CCP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(1),pwmas(1),gnu(1),1)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='cloud-gamma-dm;'
   cdunits='microns;'
elseif(cvar(1:lv).eq.'cloud_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CCP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(1),pwmas(1),gnu(1),2)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='cloud-gamma-d0;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'rain_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(2),pwmas(2),gnu(2),1)
   cdname='rain-gamma-dm;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'rain_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(2),pwmas(2),gnu(2),2)
   cdname='rain-gamma-d0;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'rain_gam_lognw') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(2),pwmas(2),gnu(2),3)
   cdname='rain-gamma-log10nw;'
   cdunits='1/mm x 1/m3;'
elseif(cvar(1:lv).eq.'rain_gam_sigma') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(2),pwmas(2),gnu(2),4)
   cdname='rain-gamma-sigma;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'pris_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CPP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(3),pwmas(3),gnu(3),1)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='pris-gamma-dm;'
   cdunits='microns;'
elseif(cvar(1:lv).eq.'pris_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CPP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(3),pwmas(3),gnu(3),2)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='pris-gamma-d0;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'snow_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(4),pwmas(4),gnu(4),1)
   cdname='snow-gamma-dm;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'snow_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(4),pwmas(4),gnu(4),2)
   cdname='snow-gamma-d0;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'snow_gam_lognw') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(4),pwmas(4),gnu(4),3)
   cdname='snow-gamma-log10nw;'
   cdunits='1/mm x 1/m3;'
elseif(cvar(1:lv).eq.'snow_gam_sigma') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(4),pwmas(4),gnu(4),4)
   cdname='snow-gamma-sigma;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'aggr_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CAP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(5),pwmas(5),gnu(5),1)
   cdname='aggr-gamma-dm;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'aggr_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CAP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(5),pwmas(5),gnu(5),2)
   cdname='aggr-gamma-d0;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'aggr_gam_lognw') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CAP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(5),pwmas(5),gnu(5),3)
   cdname='aggr-gamma-log10nw;'
   cdunits='1/mm x 1/m3;'
elseif(cvar(1:lv).eq.'aggr_gam_sigma') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CAP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(5),pwmas(5),gnu(5),4)
   cdname='aggr-gamma-sigma;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'grau_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CGP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(6),pwmas(6),gnu(6),1)
   cdname='graup-gamma-dm;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'grau_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CGP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(6),pwmas(6),gnu(6),2)
   cdname='graup-gamma-d0;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'grau_gam_lognw') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CGP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(6),pwmas(6),gnu(6),3)
   cdname='graup-gamma-log10nw;'
   cdunits='1/mm x 1/m3;'
elseif(cvar(1:lv).eq.'grau_gam_sigma') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CGP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(6),pwmas(6),gnu(6),4)
   cdname='graup-gamma-sigma;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'hail_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CHP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(7),pwmas(7),gnu(7),1)
   cdname='hail-gamma-dm;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'hail_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CHP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(7),pwmas(7),gnu(7),2)
   cdname='hail-gamma-d0;'
   cdunits='mm;'
elseif(cvar(1:lv).eq.'hail_gam_lognw') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CHP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(7),pwmas(7),gnu(7),3)
   cdname='hail-gamma-log10nw;'
   cdunits='1/mm x 1/m3;'
elseif(cvar(1:lv).eq.'hail_gam_sigma') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CHP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(7),pwmas(7),gnu(7),4)
   cdname='hail-gamma-sigma;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'driz_gam_dm') then
   ivar_type=3
   ierr=rams_getvar('RDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CDP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(8),pwmas(8),gnu(8),1)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='driz-gamma-dm;'
   cdunits='microns;'
elseif(cvar(1:lv).eq.'driz_gam_d0') then
   ivar_type=3
   ierr=rams_getvar('RDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CDP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mult (n1,n2,n3,f,d)
   CALL rams_comp_hydrogamma (n1,n2,n3,a,f,cfmas(8),pwmas(8),gnu(8),2)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='driz-gamma-d0;'
   cdunits='microns;'


!#####################################################################
! 3D MOISTURE MASS MIXING RATIOS AND HUMIDITY
!#####################################################################
elseif(cvar(1:lv).eq.'vapr_press') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,a,flnm)
   CALL rams_comp_vaporpress (n1,n2,n3,a,c)
   cdname='vapor-pressure;'
   cdunits='mb;'

elseif(cvar(1:lv).eq.'rslf') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('RV',idim_type,ngrd,d,flnm)
   CALL rams_comp_rslf (n1,n2,n3,a,c,d)
   cdname='liquid-supersaturation;'
   cdunits='percent;'

elseif(cvar(1:lv).eq.'rsif') then
   ivar_type=3
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('RV',idim_type,ngrd,d,flnm)
   CALL rams_comp_rsif (n1,n2,n3,a,c,d)
   cdname='ice-supersaturation;'
   cdunits='percent;'

elseif(cvar(1:lv).eq.'vapor') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='vapor-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vapor_m3') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='vapor-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'cloud') then
   ivar_type=3
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='cloud-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'cloud_m3') then
   ivar_type=3
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='cloud-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'rain') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='rain-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'rain_m3') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='rain-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'pristine') then
   ivar_type=3
   ierr=rams_getvar('RPP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='pristine-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'pristine_m3') then
   ivar_type=3
   ierr=rams_getvar('RPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'snow') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='snow-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'snow_m3') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snow-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'aggregates') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aggregate-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'aggregates_m3') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='aggregate-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'graupel') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='graupel-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'graupel_m3') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='graupel-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'hail') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='hail-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'hail_m3') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='hail-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'drizzle') then
   ivar_type=3
   ierr=rams_getvar('RDP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='drizzle-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'drizzle_m3') then
   ivar_type=3
   ierr=rams_getvar('RDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='drizzle-mixing-ratio;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'prissnowagg') then
   ivar_type=3
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snowprisagg-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'grauphail') then
   ivar_type=3
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='grauphail-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'liquid') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)

   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      ierr=rams_getvar('Q6',idim_type,ngrd,d,flnm)
      if(ierr.eq.0) then
         CALL rams_comp_fracliq (n1,n2,n3,d)
         CALL rams_comp_mult (n1,n2,n3,c,d)
      endif
      CALL rams_comp_accum (n1,n2,n3,a,c)
   endif

   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      ierr=rams_getvar('Q7',idim_type,ngrd,d,flnm)
      if(ierr.eq.0) then
         CALL rams_comp_fracliq (n1,n2,n3,d)
         CALL rams_comp_mult (n1,n2,n3,c,d)
      endif
      CALL rams_comp_accum (n1,n2,n3,a,c)
    endif

   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='liquid-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'ice') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)

   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      ierr=rams_getvar('Q6',idim_type,ngrd,d,flnm)
      if(ierr.eq.0) then
         CALL rams_comp_fracice (n1,n2,n3,d)
         CALL rams_comp_mult (n1,n2,n3,c,d)
      endif
      CALL rams_comp_accum (n1,n2,n3,a,c)
   endif

   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      ierr=rams_getvar('Q7',idim_type,ngrd,d,flnm)
      if(ierr.eq.0) then
         CALL rams_comp_fracice (n1,n2,n3,d)
         CALL rams_comp_mult (n1,n2,n3,c,d)
      endif
      CALL rams_comp_accum (n1,n2,n3,a,c)
   endif

   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ice-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'total_cond' .or. cvar(1:lv).eq.'total_cond_m3') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)

   if(cvar(1:lv).eq.'total_cond_m3')then
     ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
     CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
     CALL rams_comp_mult (n1,n2,n3,a,d)
   endif

   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='total-condensate-mixing-ratio;'

   if(cvar(1:lv).eq.'total_cond')    cdunits='g/kg;'
   if(cvar(1:lv).eq.'total_cond_m3') cdunits='g/m3;'

elseif(cvar(1:lv).eq.'total_mixr') then
   ivar_type=3
   ierr=rams_getvar('RTP',idim_type,ngrd,a,flnm)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='total-water-mixing-ratio-RTP;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'total_mixr_m3') then
   ivar_type=3
   ierr=rams_getvar('RTP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegm (n1,n2,n3,a)  !##### "nonegm" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='total-water-mixing-ratio-RTP;'
   cdunits='g/m3;'

elseif(cvar(1:lv).eq.'ctop_tempc_sstbase') then
   ivar_type=2
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)

   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,e,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,e)
   CALL rams_comp_tempC (n1,n2,n3,d)

   ierr=rams_getvar('SOIL_ENERGY',idim_type,ngrd,c,flnm)
   kp = nzg
   CALL rams_fill_sst (n1,n2,nzg*npatch,kp,f,c)

   CALL rams_comp_cloudtop_sstbase (n1,n2,n3,a,d,f)

   cdname='cloud-top-temperature-sstsurface;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'ctop_tempc_nobase') then
   ivar_type=2
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)

   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,e,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,e)
   CALL rams_comp_tempC (n1,n2,n3,d)

   CALL rams_comp_cloudtop_nobase (n1,n2,n3,a,d)

   cdname='cloud-top-temperature-nosurface;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'dewptk') then
   ivar_type=3
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_dewK (n1,n2,n3,a,c,d)
   cdname='dewpoint-temperature;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'dewptf') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_dewK (n1,n2,n3,a,c,d)
   CALL rams_comp_tempF (n1,n2,n3,a)
   cdname='dewpoint-temperature;'
   cdunits='F;'

elseif(cvar(1:lv).eq.'dewptc') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_dewK (n1,n2,n3,a,c,d)
   CALL rams_comp_tempC (n1,n2,n3,a)
   cdname='dewpoint-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'relhum') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_rh (n1,n2,n3,a,c,d)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='relative-humidity;'
   cdunits='pct;'

elseif(cvar(1:lv).eq.'relhum_frac') then
   ivar_type=3
   ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_rh (n1,n2,n3,a,c,d)
   CALL rams_comp_mults (n1,n2,n3,a,.01)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='relative-humidity;'
   cdunits='frac;'

elseif(cvar(1:lv).eq.'clear_frac') then
   ivar_type=2
   ierr=rams_getvar('RV',idim_type,ngrd,b,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_rh (n1,n2,n3,b,c,d)
   CALL rams_comp_noneg (n1,n2,n3,b)
   CALL cldfraction (n1,n2,n3,a,c,b)
   cdname='clear-sky;'
   cdunits='frac;'

elseif(cvar(1:lv).eq.'cloud_frac') then
   ivar_type=2
   ierr=rams_getvar('RV',idim_type,ngrd,b,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   CALL rams_comp_rh (n1,n2,n3,b,c,d)
   CALL rams_comp_noneg (n1,n2,n3,b)
   CALL cldfraction (n1,n2,n3,a,c,b)
   CALL rams_comp_1minus (n1,n2,n3,a)
   cdname='cloud-cover;'
   cdunits='frac;'

!#####################################################################
!3D HYDROMETEOR NUMBER CONCENTRATIONS
!#####################################################################
elseif(cvar(1:lv).eq.'cloud_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RCP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='cloud-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'cloud_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RCP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='cloud-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'rain_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RRP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='rain-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'pris_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RPP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'pris_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RPP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'snow_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RSP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snow-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'agg_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RAP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='aggregate-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'graup_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RGP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='graupel-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'hail_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RHP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='hail-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'drizzle_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RDP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='drizzle-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'drizzle_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RDP',idim_type,ngrd,f,flnm)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='drizzle-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'cloud_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RCP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='cloud-concen;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'rain_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='rain-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'rain_concen_dm3') then
   ivar_type=3
   ierr=rams_getvar('CRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RRP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='rain-concen;'
   cdunits='#/dm3;'

elseif(cvar(1:lv).eq.'pris_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RPP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'pris_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RPP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-concen;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'snow_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snow-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'snow_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RSP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snow-concen;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'agg_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RAP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='aggregates-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'graup_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RGP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='graupel-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'hail_concen_m3') then
   ivar_type=3
   ierr=rams_getvar('CHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RHP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='hail-concen;'
   cdunits='#/m3;'

elseif(cvar(1:lv).eq.'drizzle_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('RDP',idim_type,ngrd,f,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_nonegn (n1,n2,n3,a,f)  !##### "nonegn" minimum
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='drizzle-concen;'
   cdunits='#/cm3;'

!#####################################################################
!HUCM-SBM EXTRA HYDROMETEOR QUANTITIES
!#####################################################################
elseif(cvar(1:lv).eq.'ice_plates') then
   ivar_type=3
   ierr=rams_getvar('RIPP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='plates-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'ice_columns') then
   ivar_type=3
   ierr=rams_getvar('RICP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='columns-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'ice_dendrites') then
   ivar_type=3
   ierr=rams_getvar('RIDP',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dendrites-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'plates_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CIPP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='plates-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'plates_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CIPP',idim_type,ngrd,a,flnm)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='plates-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'columns_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CICP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='columns-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'columns_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CICP',idim_type,ngrd,a,flnm)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='columns-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'dendrites_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CIDP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='dendrites-concen;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'dendrites_concen_kg') then
   ivar_type=3
   ierr=rams_getvar('CIDP',idim_type,ngrd,a,flnm)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='dendrites-concen;'
   cdunits='#/kg;'

elseif(cvar(1:lv).eq.'pcpvip') then
   ivar_type=3
   ierr=rams_getvar('PCPVIP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-iceplates-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvic') then
   ivar_type=3
   ierr=rams_getvar('PCPVIC',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-icecolumns-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvid') then
   ivar_type=3
   ierr=rams_getvar('PCPVID',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-icedendrites-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprip') then
   ivar_type=2
   ierr=rams_getvar('PCPRIP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='iceplates-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpric') then
   ivar_type=2
   ierr=rams_getvar('PCPRIC',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='icecolumns-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprid') then
   ivar_type=2
   ierr=rams_getvar('PCPRID',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='icedendrites-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'accpip') then
   ivar_type=2
   ierr=rams_getvar('ACCPIP',idim_type,ngrd,a,flnm)
   cdname='accum-iceplates;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpic') then
   ivar_type=2
   ierr=rams_getvar('ACCPIC',idim_type,ngrd,a,flnm)
   cdname='accum-icecolumns;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpid') then
   ivar_type=2
   ierr=rams_getvar('ACCPID',idim_type,ngrd,a,flnm)
   cdname='accum-icedendrites;'
   cdunits='kg/m2;'

!#####################################################################
!3D AEROSOLS NUMBER, MASS, SIZE, SOLUBILITY
!#####################################################################
elseif(cvar(1:lv).eq.'ifn_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CIFNP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ice-nuclei-concentration;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'ifn_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CIFNP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ice-nuclei-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ccn1_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CN1NP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ccn-mode-1-concentration;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'ccn1_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CN1NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ccn-mode-1-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ccn2_concen_mg') then
   ivar_type=3
   ierr=rams_getvar('CN2NP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ccn-mode-2-concentration;'
   cdunits='#/mg;'

elseif(cvar(1:lv).eq.'ccn2_concen_cm3') then
   ivar_type=3
   ierr=rams_getvar('CN2NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ccn-mode-2-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'dust1_concen') then
   ivar_type=3
   ierr=rams_getvar('MD1NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust1-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'dust2_concen') then
   ivar_type=3
   ierr=rams_getvar('MD2NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust2-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'abs_carbon1_concen') then
   ivar_type=3
   ierr=rams_getvar('ABC1NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='absorbing-carbon1-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'abs_carbon2_concen') then
   ivar_type=3
   ierr=rams_getvar('ABC2NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='absorbing-carbon2-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'salt_film_concen') then
   ivar_type=3
   ierr=rams_getvar('SALT_FILM_NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-film-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'salt_jet_concen') then
   ivar_type=3
   ierr=rams_getvar('SALT_JET_NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-jet-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'salt_spume_concen') then
   ivar_type=3
   ierr=rams_getvar('SALT_SPUM_NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-spume-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'regen_aero1_concen') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO1_NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regenerated-aero1-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'regen_aero2_concen') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO2_NP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regenerated-aero2-concentration;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ccn1_mass') then
   ivar_type=3
   ierr=rams_getvar('CN1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ccn-mode-1-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'ccn2_mass') then
   ivar_type=3
   ierr=rams_getvar('CN2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ccn-mode-2-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust1_mass') then
   ivar_type=3
   ierr=rams_getvar('MD1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust1-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust1_massd10') then
   ivar_type=3
   ierr=rams_getvar('MD1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e8)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust1-mass;'
   cdunits='micro-grams/m3/10.0;'

elseif(cvar(1:lv).eq.'dust2_mass') then
   ivar_type=3
   ierr=rams_getvar('MD2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust2-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust2_massd10') then
   ivar_type=3
   ierr=rams_getvar('MD2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e8)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust2-mass;'
   cdunits='micro-grams/m3/10.0;'

elseif(cvar(1:lv).eq.'abs_carbon1_mass') then
   ivar_type=3
   ierr=rams_getvar('ABC1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='absorbing-carbon1-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'abs_carbon2_mass') then
   ivar_type=3
   ierr=rams_getvar('ABC2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='absorbing-carbon2-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'salt_film_mass') then
   ivar_type=3
   ierr=rams_getvar('SALT_FILM_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-film-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'salt_jet_mass') then
   ivar_type=3
   ierr=rams_getvar('SALT_JET_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-jet-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'salt_spume_mass') then
   ivar_type=3
   ierr=rams_getvar('SALT_SPUM_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='salt-spume-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'regen_aero1_mass') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO1_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regenerated-aero1-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'regen_aero2_mass') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO2_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regenerated-aero2-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'resol_aero1_mass') then
   ivar_type=3
   ierr=rams_getvar('RESOL_AERO1_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regen-soluble-aero1-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'resol_aero2_mass') then
   ivar_type=3
   ierr=rams_getvar('RESOL_AERO2_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='regen-soluble-aero2-mass;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'regen1_epsilon') then
   ivar_type=3
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RESOL_AERO1_MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_zero (n1,n2,n3,d)
   ierr=rams_getvar('REGEN_AERO1_MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   CALL rams_comp_aeroepsilon (n1,n2,n3,a,d)
   cdname='regen1-solubility-fraction;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'regen2_epsilon') then
   ivar_type=3
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RESOL_AERO2_MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_zero (n1,n2,n3,d)
   ierr=rams_getvar('REGEN_AERO2_MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   CALL rams_comp_aeroepsilon (n1,n2,n3,a,d)
   cdname='regen2-solubility-fraction;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'ccn1_medrad') then
   ivar_type=3
   ierr=rams_getvar('CN1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CN2NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,1769.) !1769 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ccn1-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'ccn2_medrad') then
   ivar_type=3
   ierr=rams_getvar('CN2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CN2NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,1769.) !1769 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='ccn2-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'dust1_medrad') then
   ivar_type=3
   ierr=rams_getvar('MD1MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('MD1NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2500.) !2500 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='dust1-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'dust2_medrad') then
   ivar_type=3
   ierr=rams_getvar('MD2MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('MD2NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2650.) !2650 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='dust2-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'salt_film_medrad') then
   ivar_type=3
   ierr=rams_getvar('SALT_FILM_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('SALT_FILM_NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2165.) !2165 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='salt-film-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'salt_jet_medrad') then
   ivar_type=3
   ierr=rams_getvar('SALT_JET_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('SALT_JET_NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2165.) !2165 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='salt-jet-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'salt_spume_medrad') then
   ivar_type=3
   ierr=rams_getvar('SALT_SPUM_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('SALT_SPUM_NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2165.) !2165 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='salt-spume-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'regen_aero1_medrad') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO1_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('REGEN_AERO1_NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2165.) !2165 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='regenerated-aero1-median-radius;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'regen_aero2_medrad') then
   ivar_type=3
   ierr=rams_getvar('REGEN_AERO2_MP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('REGEN_AERO2_NP',idim_type,ngrd,c,flnm)
   CALL rams_comp_aeromedrad (n1,n2,n3,a,c,2165.) !2165 kg/m3 is solute density
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='regenerated-aero2-median-radius;'
   cdunits='microns;'

!#####################################################################
!3D AEROSOLS TRACKING VARIABLES
!#####################################################################
elseif(cvar(1:lv).eq.'aerosol_cloud_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-cloud-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_rain_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-rain-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_pris_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-pristineice;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_snow_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-snow;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_aggr_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-aggregates;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_grau_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-graupel;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_hail_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-hail;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_driz_mass') then
   ivar_type=3
   ierr=rams_getvar('CNMDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-drizzle;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aerosol_hydro_mass') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('CNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('CNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='aerosol-mass-in-hydrometeors;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_cloud_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-cloud-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_rain_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-rain-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_pris_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-pristineice;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_snow_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-snow;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_aggr_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-aggregates;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_grau_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-graupel;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_hail_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-hail;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_driz_mass') then
   ivar_type=3
   ierr=rams_getvar('SNMDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='soluble-mass-in-drizzle;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'soluble_hydro_mass') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('SNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='ccn-soluble-mass-in-hydro;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'aero_epsilon') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('SNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('SNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_zero (n1,n2,n3,d)
   ierr=rams_getvar('CNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   ierr=rams_getvar('CNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,d,c)
   CALL rams_comp_aeroepsilon (n1,n2,n3,a,d)
   cdname='solubility-fraction;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'dust_cloud_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-cloud-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_rain_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-rain-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_pris_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-pristineice;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_snow_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-snow;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_aggr_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-aggregates;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_grau_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-graupel;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_hail_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-hail;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_driz_mass') then
   ivar_type=3
   ierr=rams_getvar('DNMDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-drizzle;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dust_hydro_mass') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dust-mass-in-hydrometeors;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_cloud_mass') then
   ivar_type=3
   ierr=rams_getvar('DINCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-cloud-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_rain_mass') then
   ivar_type=3
   ierr=rams_getvar('DINRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-rain-drops;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_pris_mass') then
   ivar_type=3
   ierr=rams_getvar('DINPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-pristineice;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_snow_mass') then
   ivar_type=3
   ierr=rams_getvar('DINSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-snow;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_aggr_mass') then
   ivar_type=3
   ierr=rams_getvar('DINAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-aggregates;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_grau_mass') then
   ivar_type=3
   ierr=rams_getvar('DINGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-graupel;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_hail_mass') then
   ivar_type=3
   ierr=rams_getvar('DINHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-hail;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_driz_mass') then
   ivar_type=3
   ierr=rams_getvar('DINDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-drizzle;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'dustifn_hydro_mass') then
   ivar_type=3
   iany=1
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DINCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DINDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='dustifn-mass-in-hydrometeors;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'ifn_nuc_numtrack') then
   ivar_type=3
   ierr=rams_getvar('IFNNUCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='IFN-already-nucleated-DeMott;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ifn_incloud') then
   ivar_type=3
   ierr=rams_getvar('IMMERCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='IFN-within-cloud-DeMott;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ifn_indriz') then
   ivar_type=3
   ierr=rams_getvar('IMMERDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='IFN-dust-within-drizzle-DeMott;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'ifn_inrain') then
   ivar_type=3
   ierr=rams_getvar('IMMERRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='IFN-dust-within-rain-DeMott;'
   cdunits='#/cm3;'

!###############################################################################
!3D VERTICAL VELOCITY AND MICROPHYSICAL BUDGETS (INSTANTANEOUS)
!###############################################################################
elseif(cvar(1:lv).eq.'wp_advdif') then       ! 1 - instant
   ivar_type=3
   ierr=rams_getvar('WP_ADVDIF',idim_type,ngrd,a,flnm)
   cdname='W-advection-diffusion;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'wp_buoy_theta') then       ! 1 - instant
   ivar_type=3
   ierr=rams_getvar('WP_BUOY_THETA',idim_type,ngrd,a,flnm)
   cdname='W-theta-buoyancy;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'wp_buoy_cond') then       ! 1 - instant
   ivar_type=3
   ierr=rams_getvar('WP_BUOY_COND',idim_type,ngrd,a,flnm)
   cdname='W-theta-cond;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'latheatvap') then
   ivar_type=3
   ierr=rams_getvar('LATHEATVAP',idim_type,ngrd,a,flnm)
   cdname='Lat-Heat-Vap-ThetaChange-inst;'
   cdunits='dTheta;'

elseif(cvar(1:lv).eq.'latheatfrz') then
   ivar_type=3
   ierr=rams_getvar('LATHEATFRZ',idim_type,ngrd,a,flnm)
   cdname='Lat-Heat-Frz-ThetaChange-inst;'
   cdunits='dTheta;'

elseif(cvar(1:lv).eq.'nuccldr') then
   ivar_type=3
   ierr=rams_getvar('NUCCLDR',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Cloud-Nucleated-Mixing-Ratio-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'cld2rain') then
   ivar_type=3
   ierr=rams_getvar('CLD2RAIN',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Cloud-to-rain-water-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'ice2rain') then
   ivar_type=3
   ierr=rams_getvar('ICE2RAIN',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Ice-to-rain-water-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'nucicer') then
   ivar_type=3
   ierr=rams_getvar('NUCICER',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Ice-Nucleated-Mixing-Ratio-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vapliq') then
   ivar_type=3
   ierr=rams_getvar('VAPLIQ',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Liquid-Vapor-diff-evap-Mixing-Ratio-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vapice') then
   ivar_type=3
   ierr=rams_getvar('VAPICE',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Ice-Vapor-diff-evap-Mixing-Ratio-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'meltice') then
   ivar_type=3
   ierr=rams_getvar('MELTICE',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melting-of-ice-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'rimecld') then
   ivar_type=3
   ierr=rams_getvar('RIMECLD',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Rimed-Amount-from-Cloud-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'rain2ice') then
   ivar_type=3
   ierr=rams_getvar('RAIN2ICE',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Rain-Water-Collected-by-Ice-Species-inst;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'aggregate') then
   ivar_type=3
   ierr=rams_getvar('AGGREGATE',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Aggregation-of-Pris-Snow-inst;'
   cdunits='g/kg;'

!###############################################################################
!3D VERTICAL VELOCITY AND MICROPHYSICAL TOTAL BUDGETS
!These are the accumulated variables between analysis file output times
!###############################################################################

! Microphysics budget variables for IMBUDGET == 1

elseif(cvar(1:lv).eq.'nuccldrt') then
   ivar_type=3
   ierr=rams_getvar('NUCCLDRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Cloud-Nucleated-Mixing-Ratio-Total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'cld2raint') then
   ivar_type=3
   ierr=rams_getvar('CLD2RAINT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Cloud-to-rain-water-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'ice2raint') then
   ivar_type=3
   ierr=rams_getvar('ICE2RAINT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Ice-to-rain-water-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'nucicert') then
   ivar_type=3
   ierr=rams_getvar('NUCICERT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Ice-Nucleated-Mixing-Ratio-Total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapliqt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPLIQT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Liquid-Vapor-diff-evap-Mixing-Ratio-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapicet') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPICET',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Ice-Vapor-diff-evap-Mixing-Ratio-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'melticet') then
   ivar_type=3
   ierr=rams_getvar('MELTICET',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melting-of-ice-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rimecldt') then
   ivar_type=3
   ierr=rams_getvar('RIMECLDT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Rimed-Amount-from-Cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2icet') then
   ivar_type=3
   ierr=rams_getvar('RAIN2ICET',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Rain-Water-Collected-by-Ice-Species-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'aggregatet') then
   ivar_type=3
   ierr=rams_getvar('AGGREGATET',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Aggregation-of-Pris-Snow-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'latheatvapt') then
   ivar_type=3
   ierr=rams_getvar('LATHEATVAPT',idim_type,ngrd,a,flnm)
   cdname='Lat-Heat-Vap-ThetaChange-total;'
   cdunits='dTheta/time;'

elseif(cvar(1:lv).eq.'latheatfrzt') then
   ivar_type=3
   ierr=rams_getvar('LATHEATFRZT',idim_type,ngrd,a,flnm)
   cdname='Lat-Heat-Frz-ThetaChange-total;'
   cdunits='dTheta/time;'

! Extra microphysics budget variables for IMBUDGET == 2

elseif(cvar(1:lv).eq.'inuchomrt') then
   ivar_type=3
   ierr=rams_getvar('INUCHOMRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e6)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Homogeous-ice-nucleation-total;'
   cdunits='mg/kg/time;'

elseif(cvar(1:lv).eq.'inuccontrt') then
   ivar_type=3
   ierr=rams_getvar('INUCCONTRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e6)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Contact-ice-nucleation-total;'
   cdunits='mg/kg/time;'

elseif(cvar(1:lv).eq.'inucifnrt') then
   ivar_type=3
   ierr=rams_getvar('INUCIFNRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e6)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='IFN-ice-nucleation-total;'
   cdunits='mg/kg/time;'

elseif(cvar(1:lv).eq.'inuchazrt') then
   ivar_type=3
   ierr=rams_getvar('INUCHAZRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e6)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Haze-ice-nucleation-total;'
   cdunits='mg/kg/time;'

elseif(cvar(1:lv).eq.'vapcldt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPCLDT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapraint') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPRAINT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapprist') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPPRIST',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Pristine-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapsnowt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPSNOWT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Snow-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapaggrt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPAGGRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Aggregate-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapgraut') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPGRAUT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Graupel-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vaphailt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPHAILT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Hail-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'vapdrizt') then ! can be -/+ for evap vs vapdep
   ivar_type=3
   ierr=rams_getvar('VAPDRIZT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='Vapor-DepEvap-Drizzle-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'meltprist') then
   ivar_type=3
   ierr=rams_getvar('MELTPRIST',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melt-pristine-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'meltsnowt') then
   ivar_type=3
   ierr=rams_getvar('MELTSNOWT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melt-snow-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'meltaggrt') then
   ivar_type=3
   ierr=rams_getvar('MELTAGGRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melt-aggregates-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'meltgraut') then
   ivar_type=3
   ierr=rams_getvar('MELTGRAUT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melt-graupel-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'melthailt') then
   ivar_type=3
   ierr=rams_getvar('MELTHAILT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Melt-hail-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rimecldsnowt') then
   ivar_type=3
   ierr=rams_getvar('RIMECLDSNOWT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Snow-rime-cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rimecldaggrt') then
   ivar_type=3
   ierr=rams_getvar('RIMECLDAGGRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Aggr-rime-cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rimecldgraut') then
   ivar_type=3
   ierr=rams_getvar('RIMECLDGRAUT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Graupel-rime-cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rimecldhailt') then
   ivar_type=3
   ierr=rams_getvar('RIMECLDHAILT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Hail-rime-cloud-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2prt') then
   ivar_type=3
   ierr=rams_getvar('RAIN2PRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Pristine-rime-rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2snt') then
   ivar_type=3
   ierr=rams_getvar('RAIN2SNT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Snow-rime-rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2agt') then
   ivar_type=3
   ierr=rams_getvar('RAIN2AGT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Aggr-rime-rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2grt') then
   ivar_type=3
   ierr=rams_getvar('RAIN2GRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Graupel-rime-rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'rain2hat') then
   ivar_type=3
   ierr=rams_getvar('RAIN2HAT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Hail-rime-rain-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'aggrselfprist') then
   ivar_type=3
   ierr=rams_getvar('AGGRSELFPRIST',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Pristine-Selfcollect-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'aggrselfsnowt') then
   ivar_type=3
   ierr=rams_getvar('AGGRSELFSNOWT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Snow-Selfcollect-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'aggrprissnowt') then
   ivar_type=3
   ierr=rams_getvar('AGGRPRISSNOWT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Pristine-Snow-collect-total;'
   cdunits='g/kg/time;'

! Extra microphysics budget variables for IMBUDGET == 3

elseif(cvar(1:lv).eq.'dust1cldrt') then
   ivar_type=3
   ierr=rams_getvar('DUST1CLDRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Dust1-cloud-nucleation-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'dust2cldrt') then
   ivar_type=3
   ierr=rams_getvar('DUST2CLDRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Dust2-cloud-nucleation-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'dust1drzrt') then
   ivar_type=3
   ierr=rams_getvar('DUST1DRZRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Dust1-drizzle-nucleation-total;'
   cdunits='g/kg/time;'

elseif(cvar(1:lv).eq.'dust2drzrt') then
   ivar_type=3
   ierr=rams_getvar('DUST2DRZRT',idim_type,ngrd,a,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_nonegp (n1,n2,n3,a)  !##### "nonegp" minimum
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
      CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='Dust2-drizzle-nucleation-total;'
   cdunits='g/kg/time;'

! Vertically integrated variables for IMBUDGET == 1

elseif(cvar(1:lv).eq.'vt_nuccldrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('NUCCLDRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-nuccldrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_cld2raint') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('CLD2RAINT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-cld2raint;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_ice2raint') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('ICE2RAINT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-ice2raint;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_nucicert') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('NUCICERT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-nucicert;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapliqt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPLIQT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapliqt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapicet') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPICET',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapicet;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_melticet') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTICET',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-melticet;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rimecldt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RIMECLDT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rimecldt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2icet') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2ICET',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2icet;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_aggregatet') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('AGGREGATET',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-aggregatet;'
   cdunits='kg/m2/time;'

! Extra vertically integrated variables for IMBUDGET == 2

elseif(cvar(1:lv).eq.'vt_inuchomrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('INUCHOMRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-inuchomrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_inuccontrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('INUCCONTRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-inuccontrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_inucifnrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('INUCIFNRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-inucifnrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_inuchazrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('INUCHAZRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-inuchazrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapcldt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPCLDT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapcldt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapraint') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPRAINT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapraint;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapprist') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPPRIST',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapprist;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapsnowt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPSNOWT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapsnowt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapaggrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPAGGRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapaggrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapgraut') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPGRAUT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapgraut;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vaphailt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPHAILT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vaphailt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_vapdrizt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('VAPDRIZT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapdrizt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_meltprist') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTPRIST',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-meltprist;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_meltsnowt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTSNOWT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-meltsnowt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_meltaggrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTAGGRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-meltaggrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_meltgraut') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTGRAUT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-meltgraut;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_melthailt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MELTHAILT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-melthailt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rimecldsnowt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RIMECLDSNOWT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rimecldsnowt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rimecldaggrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RIMECLDAGGRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rimecldaggrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rimecldgraut') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RIMECLDGRAUT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rimecldgraut;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rimecldhailt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RIMECLDHAILT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rimecldhailt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2prt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2PRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2prt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2snt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2SNT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2snt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2agt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2AGT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2agt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2grt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2GRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2grt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_rain2hat') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAIN2HAT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain2hat;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_aggrselfprist') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('AGGRSELFPRIST',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-aggrselfprist;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_aggrselfsnowt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('AGGRSELFSNOWT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-aggrselfsnowt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_aggrprissnowt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('AGGRPRISSNOWT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-aggrprissnowt;'
   cdunits='kg/m2/time;'

! Extra vertically integrated variables for IMBUDGET == 3

elseif(cvar(1:lv).eq.'vt_dust1cldrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DUST1CLDRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-dust1cldrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_dust2cldrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DUST2CLDRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-dust2cldrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_dust1drzrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DUST1DRZRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-dust1drzrt;'
   cdunits='kg/m2/time;'

elseif(cvar(1:lv).eq.'vt_dust2drzrt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DUST2DRZRT',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-dust2drzrt;'
   cdunits='kg/m2/time;'

!######################################################################
! 3D HYDROMETEOR DIAMETERS
!######################################################################
elseif(cvar(1:lv).eq.'cloudtop_diam') then
   ivar_type=2
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CCP',idim_type,ngrd,c,flnm)
   CALL rams_comp_ctopdiam (n1,n2,n3,a,c,cfmas(1),pwmas(1))
   cdname='cloud-top-diam;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'cloud_diam') then
   ivar_type=3
   ierr=rams_getvar('RCP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CCP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(1),pwmas(1))
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='cloud-diam;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'rain_diam') then
   ivar_type=3
   ierr=rams_getvar('RRP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CRP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(2),pwmas(2))
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='rain-diam;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'pris_diam') then
   ivar_type=3
   ierr=rams_getvar('RPP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CPP',idim_type,ngrd,c,flnm)
! more general case: write habit to anal file for cfmas & pwmas index
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(3),pwmas(3))
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='pristine-diam;'
   cdunits='microns;'

elseif(cvar(1:lv).eq.'snow_diam') then
   ivar_type=3
   ierr=rams_getvar('RSP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CSP',idim_type,ngrd,c,flnm)
! more general case: write habit to anal file for cfmas & pwmas index
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(4),pwmas(4))
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='snow-diam;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'agg_diam') then
   ivar_type=3
   ierr=rams_getvar('RAP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CAP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(5),pwmas(5))
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='aggregates-diam;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'graup_diam') then
   ivar_type=3
   ierr=rams_getvar('RGP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CGP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(6),pwmas(6))
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='graupel-diam;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'hail_diam') then
   ivar_type=3
   ierr=rams_getvar('RHP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CHP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(7),pwmas(7))
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='hail-diam;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'drizzle_diam') then
   ivar_type=3
   ierr=rams_getvar('RDP',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('CDP',idim_type,ngrd,c,flnm)
   CALL rams_comp_hydrodiam (n1,n2,n3,a,c,cfmas(16),pwmas(16))
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='drizzle-diam;'
   cdunits='microns;'

!######################################################################
! 3D HYDROMETEOR TEMP, THERMAL ENERGY, LIQUID FRACTION
!######################################################################
elseif(cvar(1:lv).eq.'q2') then
   ivar_type=3
   ierr=rams_getvar('Q2',idim_type,ngrd,a,flnm)
   cdname='q2;'
   cdunits='J/kg;'

elseif(cvar(1:lv).eq.'q6') then
   ivar_type=3
   ierr=rams_getvar('Q6',idim_type,ngrd,a,flnm)
   cdname='q6;'
   cdunits='J/kg;'

elseif(cvar(1:lv).eq.'q7') then
   ivar_type=3
   ierr=rams_getvar('Q7',idim_type,ngrd,a,flnm)
   cdname='q7;'
   cdunits='J/kg;'

elseif(cvar(1:lv).eq.'rain_temp') then
   ivar_type=3
   ierr=rams_getvar('Q2',idim_type,ngrd,a,flnm)
   CALL rams_comp_raintemp (n1,n2,n3,a)
   cdname='rain-temperature;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'graup_temp') then
   ivar_type=3
   ierr=rams_getvar('Q6',idim_type,ngrd,a,flnm)
   CALL rams_comp_qtcpcp (n1,n2,n3,a)
   cdname='graupel-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'hail_temp') then
   ivar_type=3
   ierr=rams_getvar('Q7',idim_type,ngrd,a,flnm)
   CALL rams_comp_qtcpcp (n1,n2,n3,a)
   cdname='hail-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'rain_air_tempdif') then
   ivar_type=3
   ierr=rams_getvar('Q2',idim_type,ngrd,a,flnm)
   CALL rams_comp_raintemp (n1,n2,n3,a)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,c)
   CALL rams_comp_tempC (n1,n2,n3,d)
   CALL rams_comp_subt (n1,n2,n3,a,d)
   cdname='rain-air-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'graup_air_tempdif') then
   ivar_type=3
   ierr=rams_getvar('Q6',idim_type,ngrd,a,flnm)
   CALL rams_comp_qtcpcp (n1,n2,n3,a)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,c)
   CALL rams_comp_tempC (n1,n2,n3,d)
   CALL rams_comp_subt (n1,n2,n3,a,d)
   cdname='graupel-air-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'hail_air_tempdif') then
   ivar_type=3
   ierr=rams_getvar('Q7',idim_type,ngrd,a,flnm)
   CALL rams_comp_qtcpcp (n1,n2,n3,a)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,c)
   CALL rams_comp_tempC (n1,n2,n3,d)
   CALL rams_comp_subt (n1,n2,n3,a,d)
   cdname='hail-air-temp;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'graup_fracliq') then
   ivar_type=3
   ierr=rams_getvar('Q6',idim_type,ngrd,a,flnm)
   CALL rams_comp_fracliq (n1,n2,n3,a)
   cdname='graupel-liq-frac;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'hail_fracliq') then
   ivar_type=3
   ierr=rams_getvar('Q7',idim_type,ngrd,a,flnm)
   CALL rams_comp_fracliq (n1,n2,n3,a)
   cdname='hail-liq-frac;'
   cdunits='fraction;'

!######################################################################
! 3D MISCELLANEOUS FIELDS
!######################################################################
elseif(cvar(1:lv).eq.'geo') then
   ivar_type=3
   ierr=rams_getvar('TOPT',idim_type,ngrd,c,flnm)
   CALL rams_comp_z (n1,n2,n3,a,c,ngrd)
   cdname='geopotential-height;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'tke') then
   ivar_type=3
   ierr=rams_getvar('TKEP',idim_type,ngrd,a,flnm)
   CALL rams_comp_noneg (n1,n2,n3,a)
   cdname='turb-kinetic-energy;'
   cdunits='m2/s2;'

elseif(cvar(1:lv).eq.'pbl_ht') then
   ivar_type=2
   ierr=rams_getvar('TKEP',idim_type,ngrd,b,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,c,flnm)
   CALL rams_comp_pbl (n1,n2,n3,a,b,c,ngrd)
   cdname='PBL-height;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'reflect_all') then
   ivar_type=3
   iany=1
   allocate(trmix(n3,n1,n2),tgmix(n3,n1,n2),thmix(n3,n1,n2),tpmix(n3,n1,n2) &
           ,tsmix(n3,n1,n2),tamix(n3,n1,n2), trnt(n3,n1,n2), tgnt(n3,n1,n2) &
            ,thnt(n3,n1,n2), tpnt(n3,n1,n2), tsnt(n3,n1,n2), tant(n3,n1,n2) &
           ,reflc(n3,n1,n2), tden(n3,n1,n2))

   CALL reflectivity_zero (n3,n1,n2,trmix,trnt,tgmix,tgnt,thmix,thnt &
          ,tpmix,tpnt,tsmix,tsnt,tamix,tant)

   !Calculate actual density over gridpoints
   ierr=rams_getvar('PI',idim_type,ngrd,f,flnm)
   CALL rams_comp_press (n1,n2,n3,f)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,c) !d is temperature [K]
   CALL density4reflc (n3,n1,n2,tden,f,d) !after call, tden is density [kg/m^3]

   !Get rain mixing ratio and # conc
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is rain mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,trmix,1)
   ierr=rams_getvar('CRP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is rain # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,trnt,1)

   !Get graupel mixing ratio and # conc
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is graupel mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,tgmix,1)
   ierr=rams_getvar('CGP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is graupel # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,tgnt,1)

   !Get hail mixing ratio and # conc
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is hail mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,thmix,1)
   ierr=rams_getvar('CHP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is hail # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,thnt,1)

   !Get pristine ice mixing ratio and # conc
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is pris mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,tpmix,1)
   ierr=rams_getvar('CPP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is pris # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,tpnt,1)

   !Get snow ice mixing ratio and # conc
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is snow mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,tsmix,1)
   ierr=rams_getvar('CSP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is snow # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,tsnt,1)

   !Get aggregate ice mixing ratio and # conc
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,c) !c is aggr mix ratio [kg/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,c,tamix,1)
   ierr=rams_getvar('CAP',idim_type,ngrd,d,flnm)
   if(ierr.eq.0) CALL rams_comp_noneg (n1,n2,n3,d) !d is aggr # conc [#/kg]
   if(ierr.eq.0) CALL arrayswap (n3,n1,n2,d,tant,1)

!Could use this to isolate reflectivity for particular species
   !Set temporary arrays for graupel and hail to zero
!   CALL set2zero (n3,n1,n2,tgmix,thmix,tgnt,thnt)

   !Compute reflectivity values [dBZ] for combined rain, graupel and hail
   CALL reflectivity_all (n3,n1,n2,trmix,trnt,tgmix,tgnt,thmix,thnt &
          ,tpmix,tpnt,tsmix,tsnt,tamix,tant,tden,reflc)
   CALL arrayswap (n3,n1,n2,a,reflc,2) !a is now total reflctvty [dBZ]
   cdname='radar-reflectivity;'
   cdunits='dBZ;'
   deallocate(trmix,tgmix,thmix,tpmix,tsmix,tamix,trnt,tgnt,thnt,tpnt,tsnt &
             ,tant,reflc,tden)

!######################################################################
! CUMULUS PARAMETERIZATION - RADIATION - TURBULENCE
!######################################################################
elseif(cvar(1:lv).eq.'cuparm_thetasrc') then
   ivar_type=3
   ierr=rams_getvar('THSRC',idim_type,ngrd,a,flnm)
   cdname='conv-heat-rate;'
   cdunits='K/s;'

elseif(cvar(1:lv).eq.'cuparm_rtsrc') then
   ivar_type=3
   ierr=rams_getvar('RTSRC',idim_type,ngrd,a,flnm)
   cdname='conv-moist-rate;'
   cdunits='kg/kg/s;'

elseif(cvar(1:lv).eq.'khh') then
   ivar_type=3
   ierr=rams_getvar('HKH',idim_type,ngrd,a,flnm)
   cdname='horiz-diffusion-coeff;'
   cdunits='m2/s;'

elseif(cvar(1:lv).eq.'khv') then
   ivar_type=3
   ierr=rams_getvar('VKH',idim_type,ngrd,a,flnm)
   cdname='vert-diffusion-coeff;'
   cdunits='m2/s;'

elseif(cvar(1:lv).eq.'visibility') then
   ivar_type=3
   ierr=rams_getvar('BEXT',idim_type,ngrd,a,flnm)
   cdname='visibility;'
   cdunits='km;'

elseif(cvar(1:lv).eq.'swup') then
   ivar_type=3
   ierr=rams_getvar('SWUP',idim_type,ngrd,a,flnm)
   cdname='shortwave-up;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'swdn') then
   ivar_type=3
   ierr=rams_getvar('SWDN',idim_type,ngrd,a,flnm)
   cdname='shortwave-down;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'lwup') then
   ivar_type=3
   ierr=rams_getvar('LWUP',idim_type,ngrd,a,flnm)
   cdname='longwave-up;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'lwdn') then
   ivar_type=3
   ierr=rams_getvar('LWDN',idim_type,ngrd,a,flnm)
   cdname='longwave-down;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'aodt') then
   ivar_type=2
   ierr=rams_getvar('AODT',idim_type,ngrd,a,flnm)
   cdname='Visible-Band-AOD;'
   cdunits='AOD;'

elseif(cvar(1:lv).eq.'rad_thetasrc') then
   ivar_type=3
   ierr=rams_getvar('FTHRD',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,86400.)
   cdname='rad-heat-rate;'
   cdunits='K/day;'

elseif(cvar(1:lv).eq.'column_net_rad_flx') then
   ivar_type=2
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (dv1(nnxp(ngrd),nnyp(ngrd)))
   allocate (dv2(nnxp(ngrd),nnyp(ngrd)))
   allocate (dv3(nnxp(ngrd),nnyp(ngrd)))
   allocate (dv4(nnxp(ngrd),nnyp(ngrd)))
   ierr=rams_getvar('SWUP',   idim_type,ngrd,pv1,flnm) !atmos shortwave up
   ierr=rams_getvar('LWUP',   idim_type,ngrd,pv2,flnm) !atmos longwave up
   ierr=rams_getvar('SWDN',   idim_type,ngrd,pv3,flnm) !atmos shortwave down
   ierr=rams_getvar('RSHORT', idim_type,ngrd,dv1,flnm) !surface shortwave down
   ierr=rams_getvar('RLONG',  idim_type,ngrd,dv2,flnm) !surface longwave down
   ierr=rams_getvar('RLONGUP',idim_type,ngrd,dv3,flnm) !surface longwave up
   ierr=rams_getvar('ALBEDT', idim_type,ngrd,dv4,flnm) !surface albdeo
   CALL rams_net_rad_flx (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),a &
     ,pv1,pv2,pv3,dv1,dv2,dv3,dv4)
   deallocate (pv1,pv2,pv3,dv1,dv2,dv3,dv4)
   cdname='column-net-radiative-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'sum_rad_flx') then
   ivar_type=3
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   ierr=rams_getvar('SWUP',idim_type,ngrd,pv1,flnm) !atmos shortwave up
   ierr=rams_getvar('LWUP',idim_type,ngrd,pv2,flnm) !atmos longwave up
   ierr=rams_getvar('SWDN',idim_type,ngrd,pv3,flnm) !atmos shortwave down
   ierr=rams_getvar('LWDN',idim_type,ngrd,pv4,flnm) !atmos longwave down
   CALL rams_sum_rad_flx (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),a &
     ,pv1,pv2,pv3,pv4)
   deallocate (pv1,pv2,pv3,pv4)
   cdname='sum-rad-flux-up-down;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'sw_heat_rate') then
   ivar_type=3
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   ierr=rams_getvar('SWUP',idim_type,ngrd,pv1,flnm) !atmos shortwave up
   ierr=rams_getvar('SWDN',idim_type,ngrd,pv2,flnm) !atmos shortwave down
   ierr=rams_getvar('DN0' ,idim_type,ngrd,pv3,flnm)
   ierr=rams_getvar('PI'  ,idim_type,ngrd,pv4,flnm)
   CALL rams_heatrate (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),a &
     ,pv1,pv2,pv3,pv4,ngrd)
   deallocate (pv1,pv2,pv3,pv4)
   CALL rams_comp_mults (n1,n2,n3,a,86400.)
   cdname='sw_heat_rate;'
   cdunits='K/day;'

elseif(cvar(1:lv).eq.'lw_heat_rate') then
   ivar_type=3
   allocate (pv1(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv2(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv3(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   allocate (pv4(nnxp(ngrd),nnyp(ngrd),nnzp(ngrd)))
   ierr=rams_getvar('LWUP',idim_type,ngrd,pv1,flnm) !atmos longwave up
   ierr=rams_getvar('LWDN',idim_type,ngrd,pv2,flnm) !atmos longwave down
   ierr=rams_getvar('DN0' ,idim_type,ngrd,pv3,flnm)
   ierr=rams_getvar('PI'  ,idim_type,ngrd,pv4,flnm)
   CALL rams_heatrate (nnxp(ngrd),nnyp(ngrd),nnzp(ngrd),a &
     ,pv1,pv2,pv3,pv4,ngrd)
   deallocate (pv1,pv2,pv3,pv4)
   CALL rams_comp_mults (n1,n2,n3,a,86400.)
   cdname='lw_heat_rate;'
   cdunits='K/day;'

!######################################################################
! 2D SURFACE PRECIP and VERTICALLY INTEGRATED FIELDS
!######################################################################
elseif(cvar(1:lv).eq.'accpr') then
   ivar_type=2
   ierr=rams_getvar('ACCPR',idim_type,ngrd,a,flnm)
   cdname='accum-rain;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpp') then
   ivar_type=2
   ierr=rams_getvar('ACCPP',idim_type,ngrd,a,flnm)
   cdname='accum-pristine;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accps') then
   ivar_type=2
   ierr=rams_getvar('ACCPS',idim_type,ngrd,a,flnm)
   cdname='accum-snow;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpa') then
   ivar_type=2
   ierr=rams_getvar('ACCPA',idim_type,ngrd,a,flnm)
   cdname='accum-aggregates;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpg') then
   ivar_type=2
   ierr=rams_getvar('ACCPG',idim_type,ngrd,a,flnm)
   cdname='accum-graupel;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accph') then
   ivar_type=2
   ierr=rams_getvar('ACCPH',idim_type,ngrd,a,flnm)
   cdname='accum-hail;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpd') then
   ivar_type=2
   ierr=rams_getvar('ACCPD',idim_type,ngrd,a,flnm)
   cdname='accum-drizzle;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'accpaero') then
   ivar_type=2
   ierr=rams_getvar('ACCPAERO',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   cdname='accum-total-aerosol-mass;'
   cdunits='milli-grams/m2;'

elseif(cvar(1:lv).eq.'accpdust') then
   ivar_type=2
   ierr=rams_getvar('ACCPDUST',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,1.e6)
   cdname='accum-dust-aerosol-mass;'
   cdunits='milli-grams/m2;'

elseif(cvar(1:lv).eq.'dustfrac') then
   ivar_type=2
   ierr=rams_getvar('DUSTFRAC',idim_type,ngrd,a,flnm)
   cdname='dust-erodible-fraction;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'totpcp'    .or. &
       cvar(1:lv).eq.'totpcp_in' .or. &
       cvar(1:lv).eq.'precip'    .or. &
       cvar(1:lv).eq.'precip_in') then
   ivar_type=2
   iany=1
   CALL rams_comp_zero (n1,n2,1,a)
   ierr=rams_getvar('ACCPD',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPR',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPS',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPA',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPG',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('ACCPH',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)

   if (cvar(1:lv)=='precip'.or.cvar(1:lv)=='precip_in') then
      ierr=rams_getvar('ACONPR',idim_type,ngrd,c,flnm)
      if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
      cdname='total-accum-precip;'
   else
      cdname='total-resolved-precip;'
   endif

   if(cvar(1:lv)=='totpcp'.or.cvar(1:lv)=='precip') then
      cdunits='mm-liq;'
   else
      CALL rams_comp_mults (n1,n2,n3,a,.03937)
      cdunits='in-liq;'
   endif
   CALL rams_comp_noneg (n1,n2,1,a)

elseif(cvar(1:lv).eq.'pcprr') then
   ivar_type=2
   ierr=rams_getvar('PCPRR',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='rain-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvr') then
   ivar_type=3
   ierr=rams_getvar('PCPVR',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-rain-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprp') then
   ivar_type=2
   ierr=rams_getvar('PCPRP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='pristine-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvp') then
   ivar_type=3
   ierr=rams_getvar('PCPVP',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-pristine-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprs') then
   ivar_type=2
   ierr=rams_getvar('PCPRS',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='snow-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvs') then
   ivar_type=3
   ierr=rams_getvar('PCPVS',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-snow-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpra') then
   ivar_type=2
   ierr=rams_getvar('PCPRA',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='aggregates-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpva') then
   ivar_type=3
   ierr=rams_getvar('PCPVA',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-aggregates-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprg') then
   ivar_type=2
   ierr=rams_getvar('PCPRG',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='graupel-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvg') then
   ivar_type=3
   ierr=rams_getvar('PCPVG',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-graupel-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprh') then
   ivar_type=2
   ierr=rams_getvar('PCPRH',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='hail precip rate;'
   cdunits='mm/hr liq equiv;'

elseif(cvar(1:lv).eq.'pcpvh') then
   ivar_type=3
   ierr=rams_getvar('PCPVH',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-hail-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcprd') then
   ivar_type=2
   ierr=rams_getvar('PCPRD',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   cdname='drizzle-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpvd') then
   ivar_type=3
   ierr=rams_getvar('PCPVD',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,n3,a,3600.)
   cdname='3D-drizzle-precip-rate;'
   cdunits='mm/hr-liq-equiv;'

elseif(cvar(1:lv).eq.'pcpg') then
   ivar_type=2
   ierr=rams_getvar('PCPG',idim_type,ngrd,a,flnm)
   cdname='pcpg;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'qpcpg') then
   ivar_type=2
   ierr=rams_getvar('QPCPG',idim_type,ngrd,a,flnm)
   cdname='qpcpg;'
   cdunits='J/m2;'

elseif(cvar(1:lv).eq.'dpcpg') then
   ivar_type=2
   ierr=rams_getvar('DPCPG',idim_type,ngrd,a,flnm)
   cdname='dpdpg;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'pcprate'    .or. &
       cvar(1:lv).eq.'pcprate_in' .or. &
       cvar(1:lv).eq.'precipr'    .or. &
       cvar(1:lv).eq.'precipr_in') then
   ivar_type=2
   iany=1
   CALL rams_comp_zero (n1,n2,1,a)
   ierr=rams_getvar('PCPRD',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRR',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRS',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRA',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRG',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   ierr=rams_getvar('PCPRH',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
   CALL rams_comp_noneg (n1,n2,1,a)

   if (cvar(1:lv)=='precipr'.or.cvar(1:lv)=='precipr_in') then
      ierr=rams_getvar('CONPRR',idim_type,ngrd,c,flnm)
      if(ierr.eq.0) CALL rams_comp_accum (n1,n2,1,a,c)
      cdname='total-precip-rate;'
   else
      cdname='resolved-precip-rate;'
   endif

   if(cvar(1:lv)=='pcprate'.or.cvar(1:lv)=='precipr') then
      CALL rams_comp_mults (n1,n2,1,a,3600.)
      cdunits='mm/hr;'
   elseif(cvar(1:lv)=='pcprate_in'.or.cvar(1:lv)=='precipr_in') then
      CALL rams_comp_mults (n1,n2,1,a,141.732)
      cdunits='in/hr;'
   endif

elseif(cvar(1:lv).eq.'conpcp') then
   ivar_type=2
   ierr=rams_getvar('CONPRR',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,3600.)
   CALL rams_comp_noneg (n1,n2,1,a)
   cdname='convective-pcp-rate;'
   cdunits='mm/hr;'

elseif(cvar(1:lv).eq.'acccon') then
   ivar_type=2
   ierr=rams_getvar('ACONPR',idim_type,ngrd,a,flnm)
   CALL rams_comp_noneg (n1,n2,1,a)
   cdname='accum-convective-pcp;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertavg_w') then
   ivar_type=2
   ierr=rams_getvar('WP',idim_type,ngrd,c,flnm)
   CALL rams_comp_vertavg (n1,n2,n3,a,c)
   cdname='average-vertical-motion;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'vertmax_w') then
   ivar_type=2
   ierr=rams_getvar('WP',idim_type,ngrd,c,flnm)
   CALL rams_comp_vertmax (n1,n2,n3,a,c)
   cdname='maximum-vertical-motion;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'vertmin_w') then
   ivar_type=2
   ierr=rams_getvar('WP',idim_type,ngrd,c,flnm)
   CALL rams_comp_vertmin (n1,n2,n3,a,c)
   cdname='minimum-vertical-motion;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'vertmax_cloud') then
   ivar_type=2
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-cloud-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_rain') then
   ivar_type=2
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-rain-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_pris') then
   ivar_type=2
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-pristine-ice-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_snow') then
   ivar_type=2
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-snow-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_aggr') then
   ivar_type=2
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-aggregates-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_grau') then
   ivar_type=2
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-graupel-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_hail') then
   ivar_type=2
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-hail-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'vertmax_driz') then
   ivar_type=2
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) then
      CALL rams_comp_vertmax (n1,n2,n3,a,c)
      CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   endif
   cdname='maximum-drizzle-mixing-ratio;'
   cdunits='g/kg;'

! Vertically-integrated atmospheric moisture
elseif(cvar(1:lv).eq.'vertint_rt' .or. cvar(1:lv).eq.'vertint_cond') then
   ivar_type=2
   iany=1

   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)

   if (cvar(1:lv)=='vertint_rt') then
      ierr=rams_getvar('RV',idim_type,ngrd,a,flnm)
      cdname='vertically-integrated-total-water;'
   else
      CALL rams_comp_zero (n1,n2,n3,a)
      cdname='vertically-integrated-condensate;'
   endif

   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_orig') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RTP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_subt (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-condensate;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_vapor') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RV',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-vapor;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_liq') then
   ivar_type=2
   iany=1
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-liquid;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_ice') then
   ivar_type=2
   iany=1
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-ice;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_cloud') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-cloud-water;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_driz') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-drizzle;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_rain') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-rain;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_pris') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-pristine;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_snow') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-snow;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_aggr') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-aggregates;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_graupel') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-graupel;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_hail') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('RHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   cdname='vertically-integrated-hail;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'vertint_dust') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('MD1MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('MD2MP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   CALL rams_comp_mults (n1,n2,n3,a,1.e3)
   cdname='vertically-integrated-dust;'
   cdunits='grams/m2;'

elseif(cvar(1:lv).eq.'vertint_dust_hydro') then
   ivar_type=2
   iany=1
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   CALL rams_comp_dn0 (n1,n2,n3,c,b,d,e,ngrd)
   CALL rams_comp_zero (n1,n2,n3,a)
   ierr=rams_getvar('DNMCP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMRP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMPP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMSP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMAP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMGP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMHP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   ierr=rams_getvar('DNMDP',idim_type,ngrd,c,flnm)
   if(ierr.eq.0) CALL rams_comp_accum (n1,n2,n3,a,c)
   CALL rams_comp_mult (n1,n2,n3,a,d)
   CALL rams_comp_vertint (n1,n2,n3,a,e,ngrd)
   CALL rams_comp_mults (n1,n2,n3,a,1.e9)
   cdname='vertint-dust-in-hydromets;'
   cdunits='micro-grams/m2;'

!######################################################################
! 2D SEA ICE COVERAGE, DEPTH, ROUGHNESS, TEMP, SNOW COVER
!######################################################################
elseif(cvar(1:lv).eq.'snowdepthonice') then
   ivar_type=2
   ierr=rams_getvar('DEPSNOW',idim_type,ngrd,a,flnm)
   cdname='snow-depth-on-ice;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'cicedepth') then
   ivar_type=2
   ierr=rams_getvar('DEPICE',idim_type,ngrd,a,flnm)
   cdname='cice-depth;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'cicefract') then
   ivar_type=2
   ierr=rams_getvar('FRACICE',idim_type,ngrd,a,flnm)
   cdname='cice-fraction;'
   cdunits='frac;'

elseif(cvar(1:lv).eq.'cicetemp') then
   ivar_type=2
   ierr=rams_getvar('CICETP',idim_type,ngrd,a,flnm)
   cdname='cice-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'cicerough') then
   ivar_type=2
   ierr=rams_getvar('CICERUF',idim_type,ngrd,a,flnm)
   cdname='cice-roughness;'
   cdunits='#;'

!######################################################################
! 2D SURFACE HEAT, MOISTURE, MOMENTUM AND RADIATIVE FLUX
!######################################################################
elseif(cvar(1:lv).eq.'sens_flux') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_T',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,1004.)
   cdname='sfc-sens-heat-flx;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'lat_flux') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_R',idim_type,ngrd,a,flnm)
   CALL rams_comp_mults (n1,n2,1,a,2.5e6)
   cdname='sfc-lat-heat-flx;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'etrans') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_R',idim_type,ngrd,a,flnm)
!                 Divide by water density to get depth and 
!                   convert units from m/s to mm/hour (3600./1000.)
   CALL rams_comp_mults (n1,n2,1,a,3.6)
   cdname='evapo-transpiration;'
   cdunits='mm/hour;'

elseif(cvar(1:lv).eq.'etrans_in') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_R',idim_type,ngrd,a,flnm)
!                 Divide by water density to get depth and 
!                   convert units from m/s to in/hour (39.37 * 3600./1000.)
   CALL rams_comp_mults (n1,n2,n3,a,141.732)
   cdname='evapo-transpiration;'
   cdunits='in/hour;'

elseif(cvar(1:lv).eq.'umom_flx') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_U',idim_type,ngrd,a,flnm)
   cdname='sfc-u-momentum-flx;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'vmom_flx') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_V',idim_type,ngrd,a,flnm)
   cdname='sfc-v-momentum-flx;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'wmom_flx') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_W',idim_type,ngrd,a,flnm)
   cdname='sfc-w-momentum-flx;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'bowen') then
   ivar_type=2
   ierr=rams_getvar('SFLUX_T',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('SFLUX_R',idim_type,ngrd,c,flnm)
   CALL rams_comp_bowen (n1,n2,1,a,c)
   cdname='bowen-ratio;'
   cdunits=' ;'

elseif(cvar(1:lv).eq.'rshort') then
   ivar_type=2
   ierr=rams_getvar('RSHORT',idim_type,ngrd,a,flnm)
   cdname='rshort;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'rlong') then
   ivar_type=2
   ierr=rams_getvar('RLONG',idim_type,ngrd,a,flnm)
   cdname='rlong;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'rlongup') then
   ivar_type=2
   ierr=rams_getvar('RLONGUP',idim_type,ngrd,a,flnm)
   cdname='rlongup;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'albedt') then
   ivar_type=2
   ierr=rams_getvar('ALBEDT',idim_type,ngrd,a,flnm)
   cdname='albedt;'
   cdunits='fraction;'

!######################################################################
! 2D TOPOGRAPHY AND GEOGRAPHIC VALUES
!######################################################################
elseif(cvar(1:lv).eq.'topt') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,a,flnm)
   cdname='topography;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'lat') then
   ivar_type=2
   ierr=rams_getvar('GLAT',idim_type,ngrd,a,flnm)
   cdname='latitude;'
   cdunits='deg;'

elseif(cvar(1:lv).eq.'lon') then
   ivar_type=2
   ierr=rams_getvar('GLON',idim_type,ngrd,a,flnm)
   cdname='longitude;'
   cdunits='deg;'

!######################################################################
! 2D MISCELLANEOUS FIELDS
!######################################################################
elseif(cvar(1:lv).eq.'sea_press') then
   ivar_type=2
   ierr=rams_getvar('TOPT',idim_type,ngrd,a,flnm)
   CALL rams_comp_z (n1,n2,n3,c,a,ngrd)
   ierr=rams_getvar('PI',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('THETA',idim_type,ngrd,a,flnm)
   CALL rams_comp_slpress (n1,n2,n3,a,d,c,a)
   cdname='sea-level-pressure;'
   cdunits='mb;'

elseif(cvar(1:lv).eq.'sfc_div') then
   ivar_type=2
   ierr=rams_getvar('WP',idim_type,ngrd,a,flnm)
   CALL rams_comp_sfcdiv (n1,n2,n3,a,ngrd)
   cdname='surface-divergence;'
   cdunits='1/s;'

! Special use of sst: acquired for patch #1 even where no water exists
elseif(cvar(1:lv).eq.'sst') then
   ivar_type=2
   ierr=rams_getvar('SOIL_ENERGY',idim_type,ngrd,c,flnm)
   kp = nzg
   CALL rams_fill_sst (n1,n2,nzg*npatch,kp,a,c)
   cdname='water-temperature;'
   cdunits='C;'

!######################################################################
! LEAF3/SIB variables section
!######################################################################
elseif(cvar(1:lv).eq.'patch_area') then
   ivar_type=6
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   cdname='patch-fractional-area;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'water') then
   ivar_type=2
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a,flnm)
   CALL rams_comp_1plus (nnxp(ngrd),nnyp(ngrd),1,a)
   cdname='water-frac-area;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'land') then
   ivar_type=2
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a,flnm)
   CALL rams_comp_1minus (nnxp(ngrd),nnyp(ngrd),1,a)
   cdname='land-frac-area;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'snow_levels') then
   ivar_type=6
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('SFCWATER_NLEV',idim_type,ngrd,a(irecind),flnm)
   cdname='number-of-snow-levels;'
   cdunits='#;'

elseif(cvar(1:lv).eq.'snow_depth_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SFCWATER_DEPTH',idim_type,ngrd,a(irecind),flnm)
   CALL rams_sum_snowlayers_ps (nnxp(ngrd),nnyp(ngrd),nzs,npatch,a(irecind),a(1))
   cdname='snow-depth;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'snow_mass_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SFCWATER_MASS',idim_type,ngrd,a(irecind),flnm)
   CALL rams_sum_snowlayers_ps (nnxp(ngrd),nnyp(ngrd),nzs,npatch,a(irecind),a(1))
   cdname='snow-water-equivalent;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'snow_temp_ps') then
   ivar_type=4
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SFCWATER_ENERGY',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_qtcpcp_ps (nnxp(ngrd),nnyp(ngrd),nzs,npatch,a(irecind))
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),nzs,npatch,a(irecind),a(1),b)
   cdname='snowlayer-water-temp;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'topo_z0_ps') then
   ivar_type=2
   ierr=rams_getvar('TOPZO',idim_type,ngrd,a,flnm)
   cdname='topo-roughness;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'net_z0_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PATCH_ROUGH',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='net-roughness;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'soil_z0_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_ROUGH',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='soil-roughness;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'veg_z0_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_ROUGH',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-roughness;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'veg_ndvi_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar ('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar ('VEG_NDVIC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='veg-ndvi;'
   cdunits='#;'

elseif(cvar(1:lv).eq.'veg_class_bp') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('LEAF_CLASS',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_vegclass (irecsize,1,1,a(irecind))
   CALL rams_comp_bigpatch (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='dominant-vegetation-class;'
   cdunits='#;'

elseif(cvar(1:lv).eq.'veg_albedo_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_ALBEDO',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-albedo;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'veg_fracarea_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_FRACAREA',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-frac-area;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'veg_lai_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_LAI',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-area-index;'
   cdunits='#;'

elseif(cvar(1:lv).eq.'veg_disp_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_HEIGHT',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-displacement-height;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'canopy_mixrat_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAN_RVAP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_mults (n1,n2,npatch,a(irecind),1.e3)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'grnd_mixrat_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('GROUND_RSAT',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_mults (n1,n2,npatch,a(irecind),1.e3)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'soil_mixrat_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('GROUND_RVAP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_mults (n1,n2,npatch,a(irecind),1.e3)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='soil-mixing-ratio;'
   cdunits='g/kg;'

elseif(cvar(1:lv).eq.'veg_moist_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_WATER',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-moisture;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'veg_temp_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VEG_TEMP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_tempC (n1,n2,npatch,a(irecind))
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'canopy_tempc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_tempC (n1,n2,npatch,a(irecind))
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'canopy_tempf_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_tempF (n1,n2,npatch,a(irecind))
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-temperature;'
   cdunits='F;'

elseif(cvar(1:lv).eq.'ustar_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('USTAR',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ustar;'
   cdunits='m/s;'

elseif(cvar(1:lv).eq.'tstar_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('TSTAR',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='tstar;'
   cdunits='K;'

elseif(cvar(1:lv).eq.'rstar_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RSTAR',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='rstar;'
   cdunits='kg/kg;'

elseif(cvar(1:lv).eq.'sltex_bp') then
   ivar_type=5
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_TEXT',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_bigpatch (nnxp(ngrd),nnyp(ngrd),nzg,npatch,a(irecind),a(1),b)
   cdname='dominant-soil-textural-class;'
   cdunits='#;'

elseif(cvar(1:lv).eq.'soilq_ps') then
   ivar_type=5
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_ENERGY',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),nzg,npatch,a(irecind),a(1),b)
   cdname='soil-q;'
   cdunits='J/m3;'

elseif(cvar(1:lv).eq.'soil_temp_ps') then
   ivar_type=5
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_ENERGY',idim_type,ngrd,a(irecind),flnm)
   ierr=rams_getvar('SOIL_WATER',idim_type,ngrd,c,flnm)
   ierr=rams_getvar('SOIL_TEXT',idim_type,ngrd,d,flnm)
   CALL rams_comp_copysst (n1,n2,nzg,a(irecind))
   irecsizep = nnxp(ngrd) * nnyp(ngrd) * nzg
   CALL rams_comp_qwtc (n1,n2,nzg*(npatch-1),a(irecind+irecsizep)  &
      ,c(1+irecsizep),d(1+irecsizep))
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),nzg,npatch,a(irecind),a(1),b)
   cdname='soil/sea temp;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'soil_moist_ps') then
   ivar_type=5
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_WATER',idim_type,ngrd ,a(irecind),flnm)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),nzg,npatch,a(irecind),a(1),b)
   cdname='soil-moisture;'
   cdunits='m3/m3;'

elseif(cvar(1:lv).eq.'soil_moistfrac_ps') then
   ivar_type=5
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),flnm)
   ierr=rams_getvar('SOIL_TEXT',idim_type,ngrd,c,flnm)
   CALL rams_comp_slmstf (irecsize,1,1,a(irecind),c)
   CALL rams_comp_patchsum_l (nnxp(ngrd),nnyp(ngrd),nzg,npatch,a(irecind),a(1),b)
   cdname='soil-moisture-fraction;'
   cdunits='m3/m3;'

elseif(cvar(1:lv).eq.'5050_tempc_ps' .or. &
       cvar(1:lv).eq.'5050_tempf_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAN_TEMP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   ierr=rams_getvar('THETA',idim_type,ngrd,d,flnm)
   ierr=rams_getvar('PI',idim_type,ngrd,c,flnm)
   CALL rams_comp_tempK (n1,n2,n3,d,c)
   CALL rams_comp_5050 (n1,n2,n3,a,d)
   if(cvar(1:lv)=='5050_tempc_ps') then
      CALL rams_comp_tempC (n1,n2,n3,a)
      cdname='avg-canopy-airlev2-tempC;'
      cdunits='C;'
   else
      CALL rams_comp_tempF (n1,n2,n3,a)
      cdname='avg-canopy-airlev2-tempF;'
      cdunits='F;'
   endif

!######################################################################
! SIB variables section
!######################################################################
elseif(cvar(1:lv).eq.'co2_concen') then
   ivar_type=3
   ierr=rams_getvar('RCO2P',idim_type,ngrd,a,flnm)
   !See sfc_driver for SiB for multiplier 1.51724e-6 used to get ppm.
   CALL rams_comp_mults (n1,n2,n3,a,1.0/1.51724e-6)
   cdname='co2-concentration;'
   cdunits='ppm;'

elseif(cvar(1:lv).eq.'snow1_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SNOW1',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-snow;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'snow2_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SNOW2',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-surface-snow;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'capac1_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAPAC1',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='vegetation-liquid-store;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'capac2_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CAPAC2',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-surface-liquid-store;'
   cdunits='kg/m2;'

elseif(cvar(1:lv).eq.'pco2ap_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PCO2AP',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='CAS-co2-concen;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'co2flx_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('CO2FLX',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   CALL rams_comp_mults (n1,n2,1,a,1.e6)
   cdname='surface-co2-flux;'
   cdunits='umol/m2/s;'

elseif(cvar(1:lv).eq.'sfcswa_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('SFCSWA',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='surface-albedo;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'uplwrf_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('UPLWRF',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='surface-longwave-upward-rad;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'assimn_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('ASSIMN',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-uptake-of-co2;'
   cdunits='umol/m2/s;'

elseif(cvar(1:lv).eq.'respg_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RESPG',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-respiration-flux;'
   cdunits='umol/m2/s;'

elseif(cvar(1:lv).eq.'rstfac1_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RSTFAC1',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-surface-humidity-resistance-stress;'
   cdunits='#(0-1);'

elseif(cvar(1:lv).eq.'rstfac2_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RSTFAC2',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='soil-moisture-resistance-stress;'
   cdunits='#(0-1);'

elseif(cvar(1:lv).eq.'rstfac3_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RSTFAC3',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='temperature-resistance-stress;'
   cdunits='#(0-1);'

elseif(cvar(1:lv).eq.'ect_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('ECT',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='transpiration-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'eci_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('ECI',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-interception-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'egi_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('EGI',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-interception-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'egs_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('EGS',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-surface-layer-evaporation;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'hc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('HC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-sensible-heat-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'hg_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('HG',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-surface-sensible-heat-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'ra_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RA',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='CAS-to-atmos-aerodynamic-resistance;'
   cdunits='s/m;'

elseif(cvar(1:lv).eq.'rb_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RB',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-surfce-to-CAS-aerodynamic-resistance;'
   cdunits='s/m;'

elseif(cvar(1:lv).eq.'rc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='total-canopy-resistance;'
   cdunits='s/m;'

elseif(cvar(1:lv).eq.'rd_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RD',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ground-to-CAS-aerodynamic-resistance;'
   cdunits='s/m;'

elseif(cvar(1:lv).eq.'roff_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('ROFF',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='water-runoff;'
   cdunits='mm;'

elseif(cvar(1:lv).eq.'green_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('GREEN',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='greenness-fraction;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'apar_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('APAR',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='absorbed-fraction-of-PAR;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'ventmf_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('VENTMF',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='ventilation-mass-flux;'
   cdunits='kg/m2/s;'

elseif(cvar(1:lv).eq.'pco2c_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PCO2C',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-chloroplast-co2-concen;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'pco2i_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PCO2I',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-internal-co2-concen;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'pco2s_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PCO2S',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='leaf-surface-co2-concen;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'pco2m_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PCO2M',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='lowest-atmos-level-co2-concen;'
   cdunits='Pa;'

elseif(cvar(1:lv).eq.'ea_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('EA',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='canopy-water-vapor-pressure;'
   cdunits='hPa;'

elseif(cvar(1:lv).eq.'em_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('EM',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='reference-level-vapor-pressure;'
   cdunits='hPa;'

elseif(cvar(1:lv).eq.'rha_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RHA',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='CAS-relative-humidity;'
   cdunits='fraction;'

elseif(cvar(1:lv).eq.'radvbc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RADVBC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='visible-direct-radiation;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'radvdc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RADVDC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='visible-diffuse-radiation;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'radnbc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RADNBC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='NIR-direct-radiation;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'radndc_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('RADNDC',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='NIR-diffuse-radiation;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'psy_ps') then
   ivar_type=2
   irecind = 1
   irecsize = nnxp(ngrd) * nnyp(ngrd) * npatch
   ierr=rams_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),flnm)
   irecind = irecind + irecsize
   ierr=rams_getvar('PSY',idim_type,ngrd,a(irecind),flnm)
   CALL rams_comp_patchsum (nnxp(ngrd),nnyp(ngrd),1,npatch,a(irecind),a(1),b)
   cdname='psychrometric-constant;'
   cdunits='hPa/deg;'

!######################################################################
! KPP OCEAN MIXED LAYER MODEL FIELDS
!######################################################################
elseif(cvar(1:lv).eq.'kpp_hmix') then
   ivar_type=2
   ierr=rams_getvar('KPP_HMIX',idim_type,ngrd,a,flnm)
   cdname='kpp-mixed-layer-depth;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'kpp_ocdepth') then
   ivar_type=2
   ierr=rams_getvar('KPP_OCDEPTH',idim_type,ngrd,a,flnm)
   cdname='kpp-ocean-depth;'
   cdunits='m;'

elseif(cvar(1:lv).eq.'kpp_flx_ust') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_UST',idim_type,ngrd,a,flnm)
   cdname='kpp-uwind-stress;'
   cdunits='N/m2;'

elseif(cvar(1:lv).eq.'kpp_flx_vst') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_VST',idim_type,ngrd,a,flnm)
   cdname='kpp-vwind-stress;'
   cdunits='N/m2;'

elseif(cvar(1:lv).eq.'kpp_flx_nsw') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_NSW',idim_type,ngrd,a,flnm)
   cdname='kpp-shortwave-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'kpp_flx_nlw') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_NLW',idim_type,ngrd,a,flnm)
   cdname='kpp-longwave-flux;'
   cdunits='W/m2;'

elseif(cvar(1:lv).eq.'kpp_flx_ice') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_ICE',idim_type,ngrd,a,flnm)
   cdname='kpp-ice-flux;'
   cdunits='not-used;'

elseif(cvar(1:lv).eq.'kpp_flx_pcp') then
   ivar_type=2
   ierr=rams_getvar('KPP_FLX_PCP',idim_type,ngrd,a,flnm)
   cdname='kpp-freshwater-flux;'
   cdunits='mm/sec;'

elseif(cvar(1:lv).eq.'kpp_depth_temp') then
   ivar_type=2
   ierr=rams_getvar('KPP_X_T',idim_type,ngrd,a,flnm)
   cdname='kpp-depth-temperature;'
   cdunits='C;'

elseif(cvar(1:lv).eq.'kpp_depth_salinity') then
   ivar_type=2
   ierr=rams_getvar('KPP_X_S',idim_type,ngrd,a,flnm)
   cdname='kpp-depth-salinity;'
   cdunits='o/oo;'

!######################################################################
!SECTION FOR ADDED TRACERS. ADD AS MANY AS YOU NEED HERE IN ORDER.
!BY DEFAULT HERE, WE ASSUME TRACERS WERE ADDED TO IDENTICALLY INITIALIZE
!AS THE CCN, DUST1, DUST2 NUMBER AND MASS, SO OUTPUTS ARE CREATED
!ACCORDINGLY FOR COMPARISON BETWEEN MICROPHYSICALLY ACTIVE AEROSOLS AND
!THOSE THAT ACT AS PASSIVE TRACERS.
!######################################################################
elseif(cvar(1:lv).eq.'tracer001') then
   ivar_type=3
   ierr=rams_getvar('TRACERP001',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-001;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'tracer002') then
   ivar_type=3
   ierr=rams_getvar('TRACERP002',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-002;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'tracer003') then
   ivar_type=3
   ierr=rams_getvar('TRACERP003',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e-6)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-003;'
   cdunits='#/cm3;'

elseif(cvar(1:lv).eq.'tracer004') then
   ivar_type=3
   ierr=rams_getvar('TRACERP004',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-004;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'tracer005') then
   ivar_type=3
   ierr=rams_getvar('TRACERP005',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-005;'
   cdunits='micro-grams/m3;'

elseif(cvar(1:lv).eq.'tracer006') then
   ivar_type=3
   ierr=rams_getvar('TRACERP006',idim_type,ngrd,a,flnm)
   ierr=rams_getvar('TOPT',idim_type,ngrd,e,flnm)
   if(ierr.eq.0) then
    CALL rams_comp_dn0 (n1,n2,n3,b,c,d,e,ngrd)
    CALL rams_comp_mult (n1,n2,n3,a,d)
    CALL rams_comp_mults (n1,n2,n3,a,1.e9)
    CALL rams_comp_noneg (n1,n2,n3,a)
   endif
   cdname='tracer-006;'
   cdunits='micro-grams/m3;'

else

   print*,'Variable name not found in hvlib.f - ',cvar(1:lv)
   ivar_type=0

endif

if(iany.eq.1 .and. ifound.eq.0) ivar_type=0
if(iany.eq.0 .and. (ierr_getvar.eq.1.or.ifound.eq.0)) ivar_type=0

return
END SUBROUTINE rams_varlib

!##############################################################################
Subroutine rams_varinfo (num,ival,iplim)

! Use this to indicate whether the computed variable exists here - if
! missing, the plot will be skipped.

use rcommons

implicit none

integer :: num,ival,iplim

if(num.eq.1) then
   ival=ivar_type
   iplim=iptrans_lim
endif

return
END SUBROUTINE rams_varinfo
