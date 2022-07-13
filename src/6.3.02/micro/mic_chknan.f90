!##############################################################################
Subroutine checkmicro (string)

use mem_basic
use mem_micro
use mem_grid
use mem_turb
use mem_leaf
use mem_tend
use node_mod
use micphys
use rconstants

implicit none

integer :: k,i,j,prtflg,mprtflg
logical, external :: isnanr
character(len=*) :: string
real, external :: valugp
real :: tempvar

prtflg=0
mprtflg=0

do j = ja,jz
do i = ia,iz
do k = 1,mzp

   !CHECK FOR NEGATIVE VAPOR OR MIXING RATIO
   if(level>=1)then
     if(basic_g(ngrid)%rv(k,i,j)<0.0.or.basic_g(ngrid)%rtp(k,i,j)<0.0) mprtflg=1
   endif
   if(level>=2)then
     if(micro_g(ngrid)%rcp(k,i,j)<0.0) mprtflg=1
   endif
   if(level==3)then
     if((idriz  >= 1 .and. micro_g(ngrid)%rdp(k,i,j)<0.0) .or. &
        (irain  >= 1 .and. micro_g(ngrid)%rrp(k,i,j)<0.0) .or. &
        (ipris  >= 1 .and. micro_g(ngrid)%rpp(k,i,j)<0.0) .or. &
        (isnow  >= 1 .and. micro_g(ngrid)%rsp(k,i,j)<0.0) .or. &
        (iaggr  >= 1 .and. micro_g(ngrid)%rap(k,i,j)<0.0) .or. &
        (igraup >= 1 .and. micro_g(ngrid)%rgp(k,i,j)<0.0) .or. &
        (ihail  >= 1 .and. micro_g(ngrid)%rhp(k,i,j)<0.0)) mprtflg=1
   endif

   !CHECK ICE NUCLEI
   if (jnmb(3)>=5 .and. (iifn==1.or.iifn==2)) then
    if(isnanr(micro_g(ngrid)%cifnp(k,i,j))) prtflg=1
    tempvar = valugp(mzp,mxp,myp,k,i,j,tend%cifnt(1))
    if(isnanr(tempvar)) prtflg=1
    if(tempvar>1.e30) prtflg=1
   endif
   !CHECK DUST MODES
   if(idust > 0)then
     if(isnanr(micro_g(ngrid)%md1np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md1mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md2np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%md2mp(k,i,j)) ) prtflg=1
   endif
   !CHECK ABSORBING CARBON MODES
   if(iabcarb > 0)then
     if(isnanr(micro_g(ngrid)%abc1np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%abc1mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%abc2np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%abc2mp(k,i,j)) ) prtflg=1
   endif
   !CHECK SEA SALT MODES
   if(isalt > 0)then
     if(isnanr(micro_g(ngrid)%salt_film_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_jet_np(k,i,j))  .or. &
        isnanr(micro_g(ngrid)%salt_spum_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_film_mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%salt_jet_mp(k,i,j))  .or. &
        isnanr(micro_g(ngrid)%salt_spum_mp(k,i,j)) ) prtflg=1
   endif
   !CHECK CCN AND GCCN
   if(iaerosol > 0)then
     if(isnanr(micro_g(ngrid)%cccnp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%cccmp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%gccnp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%gccmp(k,i,j)) ) prtflg=1
     tempvar = valugp(mzp,mxp,myp,k,i,j,tend%cccnt(1))
     if(isnanr(tempvar)) prtflg=1
     if(tempvar>1.e30) prtflg=1
     tempvar = valugp(mzp,mxp,myp,k,i,j,tend%cccmt(1))
     if(isnanr(tempvar)) prtflg=1
     if(tempvar>1.e30) prtflg=1
     tempvar = valugp(mzp,mxp,myp,k,i,j,tend%gccnt(1))
     if(isnanr(tempvar)) prtflg=1
     if(tempvar>1.e30) prtflg=1
     tempvar = valugp(mzp,mxp,myp,k,i,j,tend%gccmt(1))
     if(isnanr(tempvar)) prtflg=1
     if(tempvar>1.e30) prtflg=1
   endif
   !CHECK REGENERATED AEROSOL MODES
   if(iccnlev>=2)then
     if(isnanr(micro_g(ngrid)%regen_aero1_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero2_np(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero1_mp(k,i,j)) .or. &
        isnanr(micro_g(ngrid)%regen_aero2_mp(k,i,j)) ) prtflg=1
   endif

   !CHECK TOTAL AEROSOL TRACKING VARIABLES
   if(iccnlev>=2)then
     if(icloud>=1)then 
       if(isnanr(micro_g(ngrid)%cnmcp(k,i,j))) prtflg=1 
     endif
     if(irain>=1)then
       if(isnanr(micro_g(ngrid)%cnmrp(k,i,j))) prtflg=1
     endif
     if(ipris>=1)then
       if(isnanr(micro_g(ngrid)%cnmpp(k,i,j))) prtflg=1
     endif
     if(isnow>=1)then
       if(isnanr(micro_g(ngrid)%cnmsp(k,i,j))) prtflg=1
     endif
     if(iaggr>=1)then
       if(isnanr(micro_g(ngrid)%cnmap(k,i,j))) prtflg=1
     endif
     if(igraup>=1)then
       if(isnanr(micro_g(ngrid)%cnmgp(k,i,j))) prtflg=1
     endif
     if(ihail>=1)then
       if(isnanr(micro_g(ngrid)%cnmhp(k,i,j))) prtflg=1
     endif
     if(idriz>=1)then
       if(isnanr(micro_g(ngrid)%cnmdp(k,i,j))) prtflg=1
     endif
     if(itrkepsilon==1) then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%snmcp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%snmrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%snmpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%snmsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%snmap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%snmgp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%snmhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%snmdp(k,i,j))) prtflg=1
       endif
       if(isnanr(micro_g(ngrid)%resol_aero1_mp(k,i,j)) .or. &
          isnanr(micro_g(ngrid)%resol_aero2_mp(k,i,j)) ) prtflg=1
     endif
     if(itrkdust==1)then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%dnmcp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%dnmrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%dnmpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%dnmsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%dnmap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%dnmgp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%dnmhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%dnmdp(k,i,j))) prtflg=1
       endif
     endif
     if(itrkdustifn==1)then
       if(icloud>=1)then
         if(isnanr(micro_g(ngrid)%dincp(k,i,j))) prtflg=1
       endif
       if(irain>=1)then
         if(isnanr(micro_g(ngrid)%dinrp(k,i,j))) prtflg=1
       endif
       if(ipris>=1)then
         if(isnanr(micro_g(ngrid)%dinpp(k,i,j))) prtflg=1
       endif
       if(isnow>=1)then
         if(isnanr(micro_g(ngrid)%dinsp(k,i,j))) prtflg=1
       endif
       if(iaggr>=1)then
         if(isnanr(micro_g(ngrid)%dinap(k,i,j))) prtflg=1
       endif
       if(igraup>=1)then
         if(isnanr(micro_g(ngrid)%dingp(k,i,j))) prtflg=1
       endif
       if(ihail>=1)then
         if(isnanr(micro_g(ngrid)%dinhp(k,i,j))) prtflg=1
       endif
       if(idriz>=1)then
         if(isnanr(micro_g(ngrid)%dindp(k,i,j))) prtflg=1
       endif
     endif
   endif

   !CHECK IMMERSION FREEZING NUCLEI TRACKING VARIABLES
   if(iifn==3 .and. iccnlev>=1)then
     if(jnmb(1) >= 5)then
       if(isnanr(micro_g(ngrid)%ifnnucp(k,i,j))) prtflg=1
     endif
     if(jnmb(1) >= 5)then
       if(isnanr(micro_g(ngrid)%immercp(k,i,j))) prtflg=1
     endif
     if(jnmb(8) >= 5)then
       if(isnanr(micro_g(ngrid)%immerdp(k,i,j))) prtflg=1
     endif
     if(jnmb(2) >= 5)then
       if(isnanr(micro_g(ngrid)%immerrp(k,i,j))) prtflg=1
     endif
   endif

   !CHECK HYDROMETEOR NUMBER CONCENTRATIONS
   if(jnmb(1) >= 5)then
     if(isnanr(micro_g(ngrid)%ccp(k,i,j))) prtflg=1
   endif
   if(jnmb(8) >= 5)then
     if(isnanr(micro_g(ngrid)%cdp(k,i,j))) prtflg=1
   endif
   if(jnmb(2) >= 5)then
     if(isnanr(micro_g(ngrid)%crp(k,i,j))) prtflg=1
   endif
   if(jnmb(3) >= 5)then
     if(isnanr(micro_g(ngrid)%cpp(k,i,j))) prtflg=1
   endif
   if(jnmb(4) >= 5)then
     if(isnanr(micro_g(ngrid)%csp(k,i,j))) prtflg=1
   endif
   if(jnmb(5) >= 5)then
     if(isnanr(micro_g(ngrid)%cap(k,i,j))) prtflg=1
   endif
   if(jnmb(6) >= 5)then
     if(isnanr(micro_g(ngrid)%cgp(k,i,j))) prtflg=1
   endif
   if(jnmb(7) >= 5)then
     if(isnanr(micro_g(ngrid)%chp(k,i,j))) prtflg=1
   endif
   !CHECK HYDROMETEOR MIXING RATIOS
   if(level>=2)then
     if(isnanr(micro_g(ngrid)%rcp(k,i,j))) prtflg=1
   endif
   if(idriz>=1)then
     if(isnanr(micro_g(ngrid)%rdp(k,i,j))) prtflg=1
   endif
   if(irain>=1)then
     if(isnanr(micro_g(ngrid)%rrp(k,i,j))) prtflg=1
   endif   
   if(ipris>=1)then
     if(isnanr(micro_g(ngrid)%rpp(k,i,j))) prtflg=1
   endif
   if(isnow>=1)then
     if(isnanr(micro_g(ngrid)%rsp(k,i,j))) prtflg=1
   endif
   if(iaggr>=1)then
     if(isnanr(micro_g(ngrid)%rap(k,i,j))) prtflg=1
   endif
   if(igraup>=1)then
     if(isnanr(micro_g(ngrid)%rgp(k,i,j))) prtflg=1
   endif
   if(ihail>=1)then
     if(isnanr(micro_g(ngrid)%rhp(k,i,j))) prtflg=1
   endif
   !CHECK 3D PRECIPITATION RATES
   if(idriz>=1)then
     if(isnanr(micro_g(ngrid)%pcpvd(k,i,j))) prtflg=1
   endif
   if(irain>=1)then
     if(isnanr(micro_g(ngrid)%pcpvr(k,i,j))) prtflg=1
   endif   
   if(ipris>=1)then
     if(isnanr(micro_g(ngrid)%pcpvp(k,i,j))) prtflg=1
   endif
   if(isnow>=1)then
     if(isnanr(micro_g(ngrid)%pcpvs(k,i,j))) prtflg=1
   endif
   if(iaggr>=1)then
     if(isnanr(micro_g(ngrid)%pcpva(k,i,j))) prtflg=1
   endif
   if(igraup>=1)then
     if(isnanr(micro_g(ngrid)%pcpvg(k,i,j))) prtflg=1
   endif
   if(ihail>=1)then
     if(isnanr(micro_g(ngrid)%pcpvh(k,i,j))) prtflg=1
   endif

   !CHECK DYNAMIC FIELDS
   if(isnanr(basic_g(ngrid)%pp(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%pc(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%up(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%uc(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%vp(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%vc(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%wp(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%wc(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%theta(k,i,j))) prtflg=1
   if(isnanr(basic_g(ngrid)%thp(k,i,j)))   prtflg=1
   if(isnanr(basic_g(ngrid)%rtp(k,i,j)))   prtflg=1
   if(isnanr(basic_g(ngrid)%rv(k,i,j)))    prtflg=1
   if(isnanr(turb_g(ngrid)%sflux_t(i,j)))  prtflg=1
   if(isnanr(turb_g(ngrid)%sflux_r(i,j)))  prtflg=1
   if(isnanr(turb_g(ngrid)%sflux_u(i,j)))  prtflg=1
   if(isnanr(turb_g(ngrid)%sflux_v(i,j)))  prtflg=1
   if(isnanr(turb_g(ngrid)%sflux_w(i,j)))  prtflg=1

   tempvar = valugp(mzp,mxp,myp,k,i,j,tend%wt(1))
   if(isnanr(tempvar)) prtflg=1

 !*******************************************************************
 !IF NEGATIVE MIXING RATIO or NAN EXISTS THEN PRINT MICRO INFORMATION
 !*******************************************************************
 if(mprtflg==1 .and. string.eq.'MICRO')then
    print*,'Checked After: ',string
    print*,'Negative Condensate MICRO (ngrid,k,i,j):' &
           ,ngrid,k,i+mi0(ngrid),j+mj0(ngrid)
    if(level>=1)then
      print*,'vapor:         ',basic_g(ngrid)%rv(k,i,j)
      print*,'rtp:           ',basic_g(ngrid)%rtp(k,i,j)
    endif
    if(level>=2)then 
      print*,'cloud:         ',micro_g(ngrid)%rcp(k,i,j)
    endif
    if(level==3)then
      if(idriz  >= 1) print*,'driz:          ',micro_g(ngrid)%rdp(k,i,j)
      if(irain  >= 1) print*,'rain:          ',micro_g(ngrid)%rrp(k,i,j)
      if(ipris  >= 1) print*,'ice:           ',micro_g(ngrid)%rpp(k,i,j)
      if(isnow  >= 1) print*,'snow:          ',micro_g(ngrid)%rsp(k,i,j)
      if(iaggr  >= 1) print*,'aggr:          ',micro_g(ngrid)%rap(k,i,j)
      if(igraup >= 1) print*,'graup:         ',micro_g(ngrid)%rgp(k,i,j)
      if(ihail  >= 1) print*,'hail:          ',micro_g(ngrid)%rhp(k,i,j)
    endif
    print*,'Try shorter timestep. These can become negative due to'
    print*,' diffusion, advection, cumulus parameterization, etc.'
 endif

 if(prtflg==1)then
    print*,'Checked After: ',string
    print*,'NAN Check (ngrid,k,j,i)',ngrid,k,j+mj0(ngrid),i+mi0(ngrid)

    print*,'wp:            ',basic_g(ngrid)%wp(k,i,j)
    print*,'wc:            ',basic_g(ngrid)%wc(k,i,j)
    print*,'w-tend:        ',valugp(mzp,mxp,myp,k,i,j,tend%wt(1))
    print*,'pp:            ',basic_g(ngrid)%pp(k,i,j)
    print*,'pc:            ',basic_g(ngrid)%pc(k,i,j)
    print*,'up:            ',basic_g(ngrid)%up(k,i,j)
    print*,'uc:            ',basic_g(ngrid)%uc(k,i,j)
    print*,'vp:            ',basic_g(ngrid)%vp(k,i,j)
    print*,'vc:            ',basic_g(ngrid)%vc(k,i,j)
    print*,'sflux_t:       ',turb_g(ngrid)%sflux_t(i,j)
    print*,'sflux_r:       ',turb_g(ngrid)%sflux_r(i,j)
    print*,'sflux_u:       ',turb_g(ngrid)%sflux_u(i,j)
    print*,'sflux_v:       ',turb_g(ngrid)%sflux_v(i,j)
    print*,'sflux_w:       ',turb_g(ngrid)%sflux_w(i,j)
    print*,'topt:          ',grid_g(ngrid)%topt(i,j)
    print*,'sst:           ',leaf_g(ngrid)%soil_energy(nzg,i,j,1)
    print*,'thp:           ',basic_g(ngrid)%thp(k,i,j)
    print*,'rtp:           ',basic_g(ngrid)%rtp(k,i,j)
    print*,'theta:         ',basic_g(ngrid)%theta(k,i,j)
    print*,'rv:            ',basic_g(ngrid)%rv(k,i,j)

    if(jnmb(3)>=5 .and. (iifn==1.or.iifn==2)) then
      print*,'cifnp:         ',micro_g(ngrid)%cifnp(k,i,j)
      print*,'cifnp-tend:    ',valugp(mzp,mxp,myp,k,i,j,tend%cifnt(1))
    endif
    if(iaerosol > 0)then
      print*,'cccnp:         ',micro_g(ngrid)%cccnp(k,i,j)
      print*,'cccnp-tend:    ',valugp(mzp,mxp,myp,k,i,j,tend%cccnt(1))
      print*,'cccmp:         ',micro_g(ngrid)%cccmp(k,i,j)
      print*,'cccmp-tend:    ',valugp(mzp,mxp,myp,k,i,j,tend%cccmt(1))
      print*,'gccnp:         ',micro_g(ngrid)%gccnp(k,i,j)
      print*,'gccnp-tend:    ',valugp(mzp,mxp,myp,k,i,j,tend%gccnt(1))
      print*,'gccmp:         ',micro_g(ngrid)%gccmp(k,i,j)
      print*,'gccmp-tend:    ',valugp(mzp,mxp,myp,k,i,j,tend%gccmt(1))
    endif
    if(idust > 0)then
      print*,'md1np:         ',micro_g(ngrid)%md1np(k,i,j)
      print*,'md1mp:         ',micro_g(ngrid)%md1mp(k,i,j)
      print*,'md2np:         ',micro_g(ngrid)%md2np(k,i,j)
      print*,'md2mp:         ',micro_g(ngrid)%md2mp(k,i,j)
    endif
    if(iabcarb > 0)then
      print*,'abc1np:         ',micro_g(ngrid)%abc1np(k,i,j)
      print*,'abc1mp:         ',micro_g(ngrid)%abc1mp(k,i,j)
      print*,'abc2np:         ',micro_g(ngrid)%abc2np(k,i,j)
      print*,'abc2mp:         ',micro_g(ngrid)%abc2mp(k,i,j)
    endif
    if(isalt > 0)then
      print*,'salt_film_np:  ',micro_g(ngrid)%salt_film_np(k,i,j)
      print*,'salt_jet_np:   ' ,micro_g(ngrid)%salt_jet_np(k,i,j)
      print*,'salt_spum_np:  ',micro_g(ngrid)%salt_spum_np(k,i,j)
      print*,'salt_film_mp:  ',micro_g(ngrid)%salt_film_mp(k,i,j)
      print*,'salt_jet_mp:   ' ,micro_g(ngrid)%salt_jet_mp(k,i,j)
      print*,'salt_spum_mp:  ',micro_g(ngrid)%salt_spum_mp(k,i,j)
    endif
    if(iccnlev>=2)then
      print*,'regen1:        ',micro_g(ngrid)%regen_aero1_np(k,i,j)
      print*,'regen2:        ',micro_g(ngrid)%regen_aero2_np(k,i,j)
      print*,'regem1:        ',micro_g(ngrid)%regen_aero1_mp(k,i,j)
      print*,'regem2:        ',micro_g(ngrid)%regen_aero2_mp(k,i,j)
      if(icloud>=1)print*,'cnmcp:         ',micro_g(ngrid)%cnmcp(k,i,j)
      if(irain>=1) print*,'cnmrp:         ',micro_g(ngrid)%cnmrp(k,i,j)
      if(ipris>=1) print*,'cnmpp:         ',micro_g(ngrid)%cnmpp(k,i,j)
      if(isnow>=1) print*,'cnmsp:         ',micro_g(ngrid)%cnmsp(k,i,j)
      if(iaggr>=1) print*,'cnmap:         ',micro_g(ngrid)%cnmap(k,i,j)
      if(igraup>=1)print*,'cnmgp:         ',micro_g(ngrid)%cnmgp(k,i,j)
      if(ihail>=1) print*,'cnmhp:         ',micro_g(ngrid)%cnmhp(k,i,j)
      if(idriz>=1) print*,'cnmdp:         ',micro_g(ngrid)%cnmdp(k,i,j)
      if(itrkdust==1)then
       if(icloud>=1)print*,'dnmcp:         ',micro_g(ngrid)%dnmcp(k,i,j)
       if(irain>=1) print*,'dnmrp:         ',micro_g(ngrid)%dnmrp(k,i,j)
       if(ipris>=1) print*,'dnmpp:         ',micro_g(ngrid)%dnmpp(k,i,j)
       if(isnow>=1) print*,'dnmsp:         ',micro_g(ngrid)%dnmsp(k,i,j)
       if(iaggr>=1) print*,'dnmap:         ',micro_g(ngrid)%dnmap(k,i,j)
       if(igraup>=1)print*,'dnmgp:         ',micro_g(ngrid)%dnmgp(k,i,j)
       if(ihail>=1) print*,'dnmhp:         ',micro_g(ngrid)%dnmhp(k,i,j)
       if(idriz>=1) print*,'dnmdp:         ',micro_g(ngrid)%dnmdp(k,i,j)
      endif
      if(itrkdustifn==1)then
       if(icloud>=1)print*,'dincp:         ',micro_g(ngrid)%dincp(k,i,j)
       if(irain>=1) print*,'dinrp:         ',micro_g(ngrid)%dinrp(k,i,j)
       if(ipris>=1) print*,'dinpp:         ',micro_g(ngrid)%dinpp(k,i,j)
       if(isnow>=1) print*,'dinsp:         ',micro_g(ngrid)%dinsp(k,i,j)
       if(iaggr>=1) print*,'dinap:         ',micro_g(ngrid)%dinap(k,i,j)
       if(igraup>=1)print*,'dingp:         ',micro_g(ngrid)%dingp(k,i,j)
       if(ihail>=1) print*,'dinhp:         ',micro_g(ngrid)%dinhp(k,i,j)
       if(idriz>=1) print*,'dindp:         ',micro_g(ngrid)%dindp(k,i,j)
      endif
      if(itrkepsilon==1)then
       print*,'regensol1:     ',micro_g(ngrid)%resol_aero1_mp(k,i,j)
       print*,'regensol2:     ',micro_g(ngrid)%resol_aero2_mp(k,i,j)
       if(icloud>=1)print*,'snmcp:         ',micro_g(ngrid)%snmcp(k,i,j)
       if(irain>=1) print*,'snmrp:         ',micro_g(ngrid)%snmrp(k,i,j)
       if(ipris>=1) print*,'snmpp:         ',micro_g(ngrid)%snmpp(k,i,j)
       if(isnow>=1) print*,'snmsp:         ',micro_g(ngrid)%snmsp(k,i,j)
       if(iaggr>=1) print*,'snmap:         ',micro_g(ngrid)%snmap(k,i,j)
       if(igraup>=1)print*,'snmgp:         ',micro_g(ngrid)%snmgp(k,i,j)
       if(ihail>=1) print*,'snmhp:         ',micro_g(ngrid)%snmhp(k,i,j)
       if(idriz>=1) print*,'snmdp:         ',micro_g(ngrid)%snmdp(k,i,j)
      endif
    endif
    if(iifn==3 .and. iccnlev>=1)then
      if(icloud>=5) print*,'ifnnucp:       ',micro_g(ngrid)%ifnnucp(k,i,j)
      if(icloud>=5) print*,'immercp:       ',micro_g(ngrid)%immercp(k,i,j)
      if(idriz>=5)  print*,'immerdp:       ',micro_g(ngrid)%immerdp(k,i,j)
      if(irain>=5)  print*,'immerrp:       ',micro_g(ngrid)%immerrp(k,i,j)
    endif
    if(icloud>=5)print*,'ccp:           ',micro_g(ngrid)%ccp(k,i,j)
    if(idriz>=5) print*,'cdp:           ',micro_g(ngrid)%cdp(k,i,j)
    if(irain>=5) print*,'crp:           ',micro_g(ngrid)%crp(k,i,j)
    if(ipris>=5) print*,'cpp:           ',micro_g(ngrid)%cpp(k,i,j)
    if(isnow>=5) print*,'csp:           ',micro_g(ngrid)%csp(k,i,j)
    if(iaggr>=5) print*,'cap:           ',micro_g(ngrid)%cap(k,i,j)
    if(igraup>=5)print*,'cgp:           ',micro_g(ngrid)%cgp(k,i,j)
    if(ihail>=5) print*,'chp:           ',micro_g(ngrid)%chp(k,i,j)
    if(icloud>=1)print*,'rcp:           ',micro_g(ngrid)%rcp(k,i,j)
    if(idriz>=1) print*,'rdp:           ',micro_g(ngrid)%rdp(k,i,j)
    if(irain>=1) print*,'rrp:           ',micro_g(ngrid)%rrp(k,i,j)
    if(ipris>=1) print*,'rpp:           ',micro_g(ngrid)%rpp(k,i,j)
    if(isnow>=1) print*,'rsp:           ',micro_g(ngrid)%rsp(k,i,j)
    if(iaggr>=1) print*,'rap:           ',micro_g(ngrid)%rap(k,i,j)
    if(igraup>=1)print*,'rgp:           ',micro_g(ngrid)%rgp(k,i,j)
    if(ihail>=1) print*,'rhp:           ',micro_g(ngrid)%rhp(k,i,j)
    if(idriz>=1) print*,'pcpvd:         ',micro_g(ngrid)%pcpvd(k,i,j)
    if(irain>=1) print*,'pcpvr:         ',micro_g(ngrid)%pcpvr(k,i,j)
    if(ipris>=1) print*,'pcpvp:         ',micro_g(ngrid)%pcpvp(k,i,j)
    if(isnow>=1) print*,'pcpvs:         ',micro_g(ngrid)%pcpvs(k,i,j)
    if(iaggr>=1) print*,'pcpva:         ',micro_g(ngrid)%pcpva(k,i,j)
    if(igraup>=1)print*,'pcpvg:         ',micro_g(ngrid)%pcpvg(k,i,j)
    if(ihail>=1) print*,'pcpvh:         ',micro_g(ngrid)%pcpvh(k,i,j)
 endif

 if(mprtflg==1 .or. prtflg==1) then
  print*,"Bad data values found! Stopping model!"
  stop
 endif

enddo
enddo
enddo

return
END SUBROUTINE checkmicro

!##############################################################################
logical Function isnanr (x)

implicit none

real :: x

isnanr = .false. 
if( (x .ne. x) .or. (x .ne. 0.0 .and. x/x .ne. 1) .or. (x*0.0 .ne. 0.0))then
   isnanr = .true.
endif
    
return
END FUNCTION isnanr

