!##############################################################################
Subroutine simdata (fname,fln)

use mem_grid
use node_mod, only:mmzp,mzp
use io_params, only:frqstate

implicit none

integer :: i
character(len=*) :: fname,fln
character(len=64) :: fullname

fullname=trim(fln)//trim(fname)

open(50,file=trim(fullname),status='unknown')
write(50,*)iyear1
write(50,*)imonth1
write(50,*)idate1
write(50,*)itime1
write(50,*)
write(50,*)dtlt
write(50,*)time
write(50,*)ngrid
write(50,*)mzp
write(50,*)
write(50,*)ngrids
do i=1,ngrids
write(50,*)mmzp(i)
enddo
do i=1,ngrids
write(50,*)frqstate(i)
enddo
close(50)

return
END SUBROUTINE simdata

!##############################################################################
Subroutine iofil (fname,x,ix,wf,fln)

use mem_grid, only:time,dtlt

implicit none

character(len=*) :: fname,fln
character(len=64) :: fullname
integer :: ix                    !number vertical levels
integer :: wf                    !read or write flag
real, dimension(ix) :: x         !variable to read/write

integer k,l

fullname=trim(fln)//trim(fname)

!READ SCM DATA
if(wf == 1) then 
  stop 'No reading of SCM data for full RAMS model'
!WRITE SCM DATA
elseif (wf >= 2) then
  open(50,file=trim(fullname),status='unknown')
  if(wf==2)write(50,*)'Elapsed time',time
  if(wf==3)write(50,*)'Elapsed time',time+dtlt
  do k=1,ix
    write(50,'(i5,e26.18)') k,x(k)
  enddo
  close(50)
endif

return
END SUBROUTINE iofil

!##############################################################################
Subroutine readwrite_scm (wf,fln)

use mem_basic
use mem_micro
use mem_grid
use mem_radiate
use mem_leaf
use micphys
use rconstants
use node_mod, only:mi0,mj0,my_rams_num,my_scm_num,mzp

implicit none

character(len=*) :: fln
integer :: m1,k,i,j,wf,ng

if(iprntstmt>=1 .and. my_rams_num==my_scm_num) &
 print*,'Writing out SCM variables at time:',time,' on grid:',ngrid

!Set the grid we are on to local variable
ng = ngrid
m1 = mzp

!Get sub-grid point associated with absolute point
i=iscmx-mi0(ng)
j=iscmy-mj0(ng)

!Compute pressure from Exner function
do k=1,m1
 pitot(k) = basic_g(ng)%pi0(k,i,j) + basic_g(ng)%pp(k,i,j)
 press(k) = p00 * (pitot(k) * cpi) ** cpor
enddo
call iofil('pressure.txt'   ,press(:),m1,wf,fln)

!Output all the potential variables in the chosen column
call iofil('soil_water.txt',leaf_g(ng)%soil_water(nzg,i,j,:),npatch,wf,fln)
call iofil('patch_area.txt',leaf_g(ng)%patch_area(i,j,:),npatch,wf,fln)
call iofil('leaf_class.txt',leaf_g(ng)%leaf_class(i,j,:),npatch,wf,fln)
call iofil('soil_text.txt' ,leaf_g(ng)%soil_text(nzg,i,j,:),npatch,wf,fln)

call iofil('zm.txt',zm,m1,wf,fln)
call iofil('zt.txt',zt,m1,wf,fln)
call iofil('glat.txt',grid_g(ng)%glat(i,j),1,wf,fln)
call iofil('glon.txt',grid_g(ng)%glon(i,j),1,wf,fln)
call iofil('topt.txt',grid_g(ng)%topt(i,j),1,wf,fln)

call iofil('up.txt'   ,basic_g(ng)%up   (:,i,j),m1,wf,fln)
call iofil('vp.txt'   ,basic_g(ng)%vp   (:,i,j),m1,wf,fln)
call iofil('wp.txt'   ,basic_g(ng)%wp   (:,i,j),m1,wf,fln)
call iofil('thp.txt'  ,basic_g(ng)%thp  (:,i,j),m1,wf,fln)
call iofil('theta.txt',basic_g(ng)%theta(:,i,j),m1,wf,fln)
call iofil('pp.txt'   ,basic_g(ng)%pp   (:,i,j),m1,wf,fln)
call iofil('rtp.txt'  ,basic_g(ng)%rtp  (:,i,j),m1,wf,fln)
call iofil('rv.txt'   ,basic_g(ng)%rv   (:,i,j),m1,wf,fln)
call iofil('dn0.txt'  ,basic_g(ng)%dn0  (:,i,j),m1,wf,fln)
call iofil('pi0.txt'  ,basic_g(ng)%pi0  (:,i,j),m1,wf,fln)

if(iswrtyp == 3 .or. ilwrtyp == 3) then
 call iofil('fthrd.txt'  ,radiate_g(ng)%fthrd  (:,i,j),m1,wf,fln)
 call iofil('fthrdp.txt' ,radiate_g(ng)%fthrdp (:,i,j),m1,wf,fln)
 call iofil('bext.txt'   ,radiate_g(ng)%bext   (:,i,j),m1,wf,fln)
 call iofil('swdn.txt'   ,radiate_g(ng)%swdn   (:,i,j),m1,wf,fln)
 call iofil('swup.txt'   ,radiate_g(ng)%swup   (:,i,j),m1,wf,fln)
 call iofil('lwdn.txt'   ,radiate_g(ng)%lwdn   (:,i,j),m1,wf,fln)
 call iofil('lwup.txt'   ,radiate_g(ng)%lwup   (:,i,j),m1,wf,fln)
 call iofil('rshort.txt' ,radiate_g(ng)%rshort (i,j),1,wf,fln)
 call iofil('rlong.txt'  ,radiate_g(ng)%rlong  (i,j),1,wf,fln)
 call iofil('rlongup.txt',radiate_g(ng)%rlongup(i,j),1,wf,fln)
 call iofil('albedt.txt' ,radiate_g(ng)%albedt (i,j),1,wf,fln)
 call iofil('cosz.txt'   ,radiate_g(ng)%cosz   (i,j),1,wf,fln)
 call iofil('aodt.txt'   ,radiate_g(ng)%aodt   (i,j),1,wf,fln)
endif

if(iaerosol > 0) then
 call iofil('cccnp.txt',micro_g(ng)%cccnp(:,i,j),m1,wf,fln)
 call iofil('cccmp.txt',micro_g(ng)%cccmp(:,i,j),m1,wf,fln)
 call iofil('gccnp.txt',micro_g(ng)%gccnp(:,i,j),m1,wf,fln)
 call iofil('gccmp.txt',micro_g(ng)%gccmp(:,i,j),m1,wf,fln)
endif

if(idust > 0) then
 call iofil('md1np.txt',micro_g(ng)%md1np(:,i,j),m1,wf,fln)
 call iofil('md2np.txt',micro_g(ng)%md2np(:,i,j),m1,wf,fln)
 call iofil('md1mp.txt',micro_g(ng)%md1mp(:,i,j),m1,wf,fln)
 call iofil('md2mp.txt',micro_g(ng)%md2mp(:,i,j),m1,wf,fln)
 if(idust == 2) call iofil('dustfrac.txt',micro_g(ng)%dustfrac(i,j),1,wf,fln)
endif

if(isalt > 0) then
 call iofil('salt_film_np.txt',micro_g(ng)%salt_film_np(:,i,j),m1,wf,fln)
 call iofil('salt_jet_np.txt' ,micro_g(ng)%salt_jet_np (:,i,j),m1,wf,fln)
 call iofil('salt_spum_np.txt',micro_g(ng)%salt_spum_np(:,i,j),m1,wf,fln)
 call iofil('salt_film_mp.txt',micro_g(ng)%salt_film_mp(:,i,j),m1,wf,fln)
 call iofil('salt_jet_mp.txt' ,micro_g(ng)%salt_jet_mp (:,i,j),m1,wf,fln)
 call iofil('salt_spum_mp.txt',micro_g(ng)%salt_spum_mp(:,i,j),m1,wf,fln)
endif

if(level == 2 .or. level == 3) then
 call iofil('rcp.txt',micro_g(ng)%rcp(:,i,j),m1,wf,fln)
endif

if(level == 3) then
 call iofil('pcpg.txt' ,micro_g(ng)%pcpg (i,j),1,wf,fln)
 call iofil('qpcpg.txt',micro_g(ng)%qpcpg(i,j),1,wf,fln)
 call iofil('dpcpg.txt',micro_g(ng)%dpcpg(i,j),1,wf,fln)
 if(ipris>=5 .and. (iifn==1.or.iifn==2)) then
  call iofil('cifnp.txt',micro_g(ng)%cifnp(:,i,j),m1,wf,fln)
 endif
 if(idriz >= 1) then
  call iofil('rdp.txt'  ,micro_g(ng)%rdp  (:,i,j),m1,wf,fln)
  call iofil('pcpvd.txt',micro_g(ng)%pcpvd(:,i,j),m1,wf,fln)
  call iofil('pcprd.txt',micro_g(ng)%pcprd(i,j),1,wf,fln)
  call iofil('accpd.txt',micro_g(ng)%accpd(i,j),1,wf,fln)
 endif
 if(irain >= 1) then 
  call iofil('rrp.txt'  ,micro_g(ng)%rrp  (:,i,j),m1,wf,fln)
  call iofil('pcpvr.txt',micro_g(ng)%pcpvr(:,i,j),m1,wf,fln)
  call iofil('q2.txt'   ,micro_g(ng)%q2   (:,i,j),m1,wf,fln)
  call iofil('pcprr.txt',micro_g(ng)%pcprr(i,j),1,wf,fln)
  call iofil('accpr.txt',micro_g(ng)%accpr(i,j),1,wf,fln)
 endif
 if(ipris >= 1) then
  call iofil('rpp.txt'  ,micro_g(ng)%rpp  (:,i,j),m1,wf,fln)
  call iofil('pcpvp.txt',micro_g(ng)%pcpvp(:,i,j),m1,wf,fln)
  call iofil('pcprp.txt',micro_g(ng)%pcprp(i,j),1,wf,fln)
  call iofil('accpp.txt',micro_g(ng)%accpp(i,j),1,wf,fln)
 endif
 if(isnow >= 1) then
  call iofil('rsp.txt'  ,micro_g(ng)%rsp  (:,i,j),m1,wf,fln)
  call iofil('pcpvs.txt',micro_g(ng)%pcpvs(:,i,j),m1,wf,fln)
  call iofil('pcprs.txt',micro_g(ng)%pcprs(i,j),1,wf,fln)
  call iofil('accps.txt',micro_g(ng)%accps(i,j),1,wf,fln)
 endif
 if(iaggr >= 1) then
  call iofil('rap.txt'  ,micro_g(ng)%rap  (:,i,j),m1,wf,fln)
  call iofil('pcpva.txt',micro_g(ng)%pcpva(:,i,j),m1,wf,fln)
  call iofil('pcpra.txt',micro_g(ng)%pcpra(i,j),1,wf,fln)
  call iofil('accpa.txt',micro_g(ng)%accpa(i,j),1,wf,fln)
 endif
 if(igraup >= 1) then
  call iofil('rgp.txt'  ,micro_g(ng)%rgp  (:,i,j),m1,wf,fln)
  call iofil('pcpvg.txt',micro_g(ng)%pcpvg(:,i,j),m1,wf,fln)
  call iofil('q6.txt'   ,micro_g(ng)%q6   (:,i,j),m1,wf,fln)
  call iofil('pcprg.txt',micro_g(ng)%pcprg(i,j),1,wf,fln)
  call iofil('accpg.txt',micro_g(ng)%accpg(i,j),1,wf,fln)
 endif
 if(ihail >=1) then
  call iofil('rhp.txt'  ,micro_g(ng)%rhp  (:,i,j),m1,wf,fln)
  call iofil('pcpvh.txt',micro_g(ng)%pcpvh(:,i,j),m1,wf,fln)
  call iofil('q7.txt'   ,micro_g(ng)%q7   (:,i,j),m1,wf,fln)
  call iofil('pcprh.txt',micro_g(ng)%pcprh(i,j),1,wf,fln)
  call iofil('accph.txt',micro_g(ng)%accph(i,j),1,wf,fln)
 endif
 if(jnmb(1) >= 5) call iofil('ccp.txt',micro_g(ng)%ccp(:,i,j),m1,wf,fln)
 if(jnmb(8) >= 5) call iofil('cdp.txt',micro_g(ng)%cdp(:,i,j),m1,wf,fln)
 if(jnmb(2) >= 5) call iofil('crp.txt',micro_g(ng)%crp(:,i,j),m1,wf,fln)
 if(jnmb(3) >= 5) call iofil('cpp.txt',micro_g(ng)%cpp(:,i,j),m1,wf,fln)
 if(jnmb(4) >= 5) call iofil('csp.txt',micro_g(ng)%csp(:,i,j),m1,wf,fln)
 if(jnmb(5) >= 5) call iofil('cap.txt',micro_g(ng)%cap(:,i,j),m1,wf,fln)
 if(jnmb(6) >= 5) call iofil('cgp.txt',micro_g(ng)%cgp(:,i,j),m1,wf,fln)
 if(jnmb(7) >= 5) call iofil('chp.txt',micro_g(ng)%chp(:,i,j),m1,wf,fln)
 if(iccnlev >= 2) then
  call iofil('regen_aero1_np.txt',micro_g(ng)%regen_aero1_np(:,i,j),m1,wf,fln)
  call iofil('regen_aero1_mp.txt',micro_g(ng)%regen_aero1_mp(:,i,j),m1,wf,fln)
  call iofil('regen_aero2_np.txt',micro_g(ng)%regen_aero2_np(:,i,j),m1,wf,fln)
  call iofil('regen_aero2_mp.txt',micro_g(ng)%regen_aero2_mp(:,i,j),m1,wf,fln)
  if(jnmb(1)>=1) call iofil('cnmcp.txt',micro_g(ng)%cnmcp(:,i,j),m1,wf,fln)
  if(jnmb(8)>=1) call iofil('cnmdp.txt',micro_g(ng)%cnmdp(:,i,j),m1,wf,fln)
  if(jnmb(2)>=1) call iofil('cnmrp.txt',micro_g(ng)%cnmrp(:,i,j),m1,wf,fln)
  if(jnmb(3)>=1) call iofil('cnmpp.txt',micro_g(ng)%cnmpp(:,i,j),m1,wf,fln)
  if(jnmb(4)>=1) call iofil('cnmsp.txt',micro_g(ng)%cnmsp(:,i,j),m1,wf,fln)
  if(jnmb(5)>=1) call iofil('cnmap.txt',micro_g(ng)%cnmap(:,i,j),m1,wf,fln)
  if(jnmb(6)>=1) call iofil('cnmgp.txt',micro_g(ng)%cnmgp(:,i,j),m1,wf,fln)
  if(jnmb(7)>=1) call iofil('cnmhp.txt',micro_g(ng)%cnmhp(:,i,j),m1,wf,fln)
  call iofil('accpaero.txt',micro_g(ng)%accpaero(i,j),1,wf,fln)
  call iofil('pcpraero.txt',micro_g(ng)%pcpraero(i,j),1,wf,fln)
  if(itrkdust==1 .and. idust>0) then
   call iofil('accpdust.txt',micro_g(ng)%accpdust(i,j),1,wf,fln)
   call iofil('pcprdust.txt',micro_g(ng)%pcprdust(i,j),1,wf,fln)
   if(jnmb(1)>=1) call iofil('dnmcp.txt',micro_g(ng)%dnmcp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) call iofil('dnmdp.txt',micro_g(ng)%dnmdp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) call iofil('dnmrp.txt',micro_g(ng)%dnmrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) call iofil('dnmpp.txt',micro_g(ng)%dnmpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) call iofil('dnmsp.txt',micro_g(ng)%dnmsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) call iofil('dnmap.txt',micro_g(ng)%dnmap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) call iofil('dnmgp.txt',micro_g(ng)%dnmgp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) call iofil('dnmhp.txt',micro_g(ng)%dnmhp(:,i,j),m1,wf,fln)
  endif
  if(itrkdustifn==1 .and. idust>0) then
   if(jnmb(1)>=1) call iofil('dincp.txt',micro_g(ng)%dincp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) call iofil('dindp.txt',micro_g(ng)%dindp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) call iofil('dinrp.txt',micro_g(ng)%dinrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) call iofil('dinpp.txt',micro_g(ng)%dinpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) call iofil('dinsp.txt',micro_g(ng)%dinsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) call iofil('dinap.txt',micro_g(ng)%dinap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) call iofil('dingp.txt',micro_g(ng)%dingp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) call iofil('dinhp.txt',micro_g(ng)%dinhp(:,i,j),m1,wf,fln)
  endif
  if(itrkepsilon==1) then
   call iofil('resol_aero1_mp.txt',micro_g(ng)%resol_aero1_mp(:,i,j),m1,wf,fln)
   call iofil('resol_aero2_mp.txt',micro_g(ng)%resol_aero2_mp(:,i,j),m1,wf,fln)
   if(jnmb(1)>=1) call iofil('snmcp.txt',micro_g(ng)%snmcp(:,i,j),m1,wf,fln)
   if(jnmb(8)>=1) call iofil('snmdp.txt',micro_g(ng)%snmdp(:,i,j),m1,wf,fln)
   if(jnmb(2)>=1) call iofil('snmrp.txt',micro_g(ng)%snmrp(:,i,j),m1,wf,fln)
   if(jnmb(3)>=1) call iofil('snmpp.txt',micro_g(ng)%snmpp(:,i,j),m1,wf,fln)
   if(jnmb(4)>=1) call iofil('snmsp.txt',micro_g(ng)%snmsp(:,i,j),m1,wf,fln)
   if(jnmb(5)>=1) call iofil('snmap.txt',micro_g(ng)%snmap(:,i,j),m1,wf,fln)
   if(jnmb(6)>=1) call iofil('snmgp.txt',micro_g(ng)%snmgp(:,i,j),m1,wf,fln)
   if(jnmb(7)>=1) call iofil('snmhp.txt',micro_g(ng)%snmhp(:,i,j),m1,wf,fln)
  endif
 endif
 if(iifn==3 .and. iccnlev>=1) then
  if(jnmb(1)>=5)call iofil('ifnnucp.txt',micro_g(ng)%ifnnucp(:,i,j),m1,wf,fln)
  if(jnmb(1)>=5)call iofil('immercp.txt',micro_g(ng)%immercp(:,i,j),m1,wf,fln)
  if(jnmb(8)>=5)call iofil('immerdp.txt',micro_g(ng)%immerdp(:,i,j),m1,wf,fln)
  if(jnmb(2)>=5)call iofil('immerrp.txt',micro_g(ng)%immerrp(:,i,j),m1,wf,fln)
 endif
 if(imbudget>=1) then !14
  call iofil('latheatvap.txt' ,micro_g(ng)%latheatvap (:,i,j),m1,wf,fln)
  call iofil('latheatfrz.txt' ,micro_g(ng)%latheatfrz (:,i,j),m1,wf,fln)
  call iofil('nuccldrt.txt'   ,micro_g(ng)%nuccldrt   (:,i,j),m1,wf,fln)
  call iofil('cld2raint.txt'  ,micro_g(ng)%cld2raint  (:,i,j),m1,wf,fln)
  call iofil('ice2raint.txt'  ,micro_g(ng)%ice2raint  (:,i,j),m1,wf,fln)
  call iofil('nucicert.txt'   ,micro_g(ng)%nucicert   (:,i,j),m1,wf,fln)
  call iofil('vapliqt.txt'    ,micro_g(ng)%vapliqt    (:,i,j),m1,wf,fln)
  call iofil('vapicet.txt'    ,micro_g(ng)%vapicet    (:,i,j),m1,wf,fln)
  call iofil('melticet.txt'   ,micro_g(ng)%melticet   (:,i,j),m1,wf,fln)
  call iofil('rimecldt.txt'   ,micro_g(ng)%rimecldt   (:,i,j),m1,wf,fln)
  call iofil('rain2icet.txt'  ,micro_g(ng)%rain2icet  (:,i,j),m1,wf,fln)
  call iofil('aggregatet.txt' ,micro_g(ng)%aggregatet (:,i,j),m1,wf,fln)
  call iofil('latheatvapt.txt',micro_g(ng)%latheatvapt(:,i,j),m1,wf,fln)
  call iofil('latheatfrzt.txt',micro_g(ng)%latheatfrzt(:,i,j),m1,wf,fln)
 endif
 if(imbudget>=2) then !29
  call iofil('inuchomrt.txt'    ,micro_g(ng)%inuchomrt    (:,i,j),m1,wf,fln)
  call iofil('inuccontrt.txt'   ,micro_g(ng)%inuccontrt   (:,i,j),m1,wf,fln)
  call iofil('inucifnrt.txt'    ,micro_g(ng)%inucifnrt    (:,i,j),m1,wf,fln)
  call iofil('inuchazrt.txt'    ,micro_g(ng)%inuchazrt    (:,i,j),m1,wf,fln)
  call iofil('vapcldt.txt'      ,micro_g(ng)%vapcldt      (:,i,j),m1,wf,fln)
  call iofil('vapraint.txt'     ,micro_g(ng)%vapraint     (:,i,j),m1,wf,fln)
  call iofil('vapprist.txt'     ,micro_g(ng)%vapprist     (:,i,j),m1,wf,fln)
  call iofil('vapsnowt.txt'     ,micro_g(ng)%vapsnowt     (:,i,j),m1,wf,fln)
  call iofil('vapaggrt.txt'     ,micro_g(ng)%vapaggrt     (:,i,j),m1,wf,fln)
  call iofil('vapgraut.txt'     ,micro_g(ng)%vapgraut     (:,i,j),m1,wf,fln)
  call iofil('vaphailt.txt'     ,micro_g(ng)%vaphailt     (:,i,j),m1,wf,fln)
  call iofil('vapdrizt.txt'     ,micro_g(ng)%vapdrizt     (:,i,j),m1,wf,fln)
  call iofil('meltprist.txt'    ,micro_g(ng)%meltprist    (:,i,j),m1,wf,fln)
  call iofil('meltsnowt.txt'    ,micro_g(ng)%meltsnowt    (:,i,j),m1,wf,fln)
  call iofil('meltaggrt.txt'    ,micro_g(ng)%meltaggrt    (:,i,j),m1,wf,fln)
  call iofil('meltgraut.txt'    ,micro_g(ng)%meltgraut    (:,i,j),m1,wf,fln)
  call iofil('melthailt.txt'    ,micro_g(ng)%melthailt    (:,i,j),m1,wf,fln)
  call iofil('rimecldsnowt.txt' ,micro_g(ng)%rimecldsnowt (:,i,j),m1,wf,fln)
  call iofil('rimecldaggrt.txt' ,micro_g(ng)%rimecldaggrt (:,i,j),m1,wf,fln)
  call iofil('rimecldgraut.txt' ,micro_g(ng)%rimecldgraut (:,i,j),m1,wf,fln)
  call iofil('rimecldhailt.txt' ,micro_g(ng)%rimecldhailt (:,i,j),m1,wf,fln)
  call iofil('rain2prt.txt'     ,micro_g(ng)%rain2prt     (:,i,j),m1,wf,fln)
  call iofil('rain2snt.txt'     ,micro_g(ng)%rain2snt     (:,i,j),m1,wf,fln)
  call iofil('rain2agt.txt'     ,micro_g(ng)%rain2agt     (:,i,j),m1,wf,fln)
  call iofil('rain2grt.txt'     ,micro_g(ng)%rain2grt     (:,i,j),m1,wf,fln)
  call iofil('rain2hat.txt'     ,micro_g(ng)%rain2hat     (:,i,j),m1,wf,fln)
  call iofil('aggrselfprist.txt',micro_g(ng)%aggrselfprist(:,i,j),m1,wf,fln)
  call iofil('aggrselfsnowt.txt',micro_g(ng)%aggrselfsnowt(:,i,j),m1,wf,fln)
  call iofil('aggrprissnowt.txt',micro_g(ng)%aggrprissnowt(:,i,j),m1,wf,fln)
 endif
 if(imbudget==3 .and. idust>0) then !4
  call iofil('dust1cldrt.txt',micro_g(ng)%dust1cldrt(:,i,j),m1,wf,fln)
  call iofil('dust2cldrt.txt',micro_g(ng)%dust2cldrt(:,i,j),m1,wf,fln)
  call iofil('dust1drzrt.txt',micro_g(ng)%dust1drzrt(:,i,j),m1,wf,fln)
  call iofil('dust2drzrt.txt',micro_g(ng)%dust2drzrt(:,i,j),m1,wf,fln)
 endif

endif !if level == 3

return
END SUBROUTINE readwrite_scm
