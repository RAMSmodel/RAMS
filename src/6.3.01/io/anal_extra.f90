!##############################################################################
Module anal_extra

implicit none

integer, parameter :: max_anextra=130
type anextra
   character(len=16) :: name
   integer :: idim_type, npointer,ind_comp
end type

integer :: num_extra_anal,max_extra_anal

type(anextra), dimension(max_anextra) :: an_extra

Contains

!##############################################################################
Subroutine anal_extra_init ()

use micphys, only:level,ipris,igraup,ihail,iifn
use micro_prm, only: iceprocs,iceflag

implicit none

!Note that if more extra variables need to be added, you will need to adjust
!the constants in the "if" statement below and then adjust the numbering 
!order for the level=4 microphysics variables. You will also need to adjust
!the "case" statements in the routine that follows this one.

!Note that when naming a new variable, it cannot share the name of a currently
!allocated variable that already exists in the HDF5 output file as a dataset
!name. This will result in a seg fault or error in writing to HDF5 output.

if (level/=4) then
   num_extra_anal=7
elseif (level == 4) then
   num_extra_anal=13 !Warm phase vars and CCN
   if (iceprocs == 1 .and. iifn == 2) num_extra_anal=num_extra_anal + 2 !IN
   if (ipris >0 .and. ipris <= 3) num_extra_anal=num_extra_anal + 6 !Pristine ice and snow
   if (ipris >= 4) num_extra_anal=num_extra_anal+10 !All three ice crystal types
   if (igraup > 0) num_extra_anal=num_extra_anal+2  !Graupel
   if (ihail > 0) num_extra_anal=num_extra_anal+2 !Hail
   if (iceprocs == 1) num_extra_anal=num_extra_anal+2 !Aggregates always included if ice is on
endif

if (num_extra_anal <= 13) then
   max_extra_anal = num_extra_anal
else
   max_extra_anal = 31
endif
an_extra(1)%name='PI';  an_extra(1)%idim_type=3; an_extra(1)%ind_comp=1
an_extra(2)%name='HKH'; an_extra(2)%idim_type=3; an_extra(2)%ind_comp=2
an_extra(3)%name='VKH'; an_extra(3)%idim_type=3; an_extra(3)%ind_comp=3
an_extra(4)%name='DN01'; an_extra(4)%idim_type=3; an_extra(4)%ind_comp=4 !placeholder
an_extra(5)%name='SOIL_WATER1';  an_extra(5)%idim_type=4; an_extra(5)%ind_comp=5 !placeholder
an_extra(6)%name='SFCWATER1';    an_extra(6)%idim_type=5; an_extra(6)%ind_comp=6 !placeholder
an_extra(7)%name='STOM_RESIST1'; an_extra(7)%idim_type=6; an_extra(7)%ind_comp=7 !placeholder
if (level==4) then
 an_extra(8)%name='RCP';     an_extra(8)%idim_type=3;   an_extra(8)%ind_comp=8
 an_extra(9)%name='RRP';     an_extra(9)%idim_type=3;   an_extra(9)%ind_comp=9
 an_extra(10)%name='CCP';    an_extra(10)%idim_type=3;  an_extra(10)%ind_comp=10 
 an_extra(11)%name='CRP';    an_extra(11)%idim_type=3;  an_extra(11)%ind_comp=11 
 an_extra(12)%name='CCCMP';  an_extra(12)%idim_type=3;  an_extra(12)%ind_comp=12 
 an_extra(13)%name='CCCNP';  an_extra(13)%idim_type=3;  an_extra(13)%ind_comp=13
 if (iceprocs == 1) then
  an_extra(14)%name='RICP';  an_extra(14)%idim_type=3; an_extra(14)%ind_comp=14
  an_extra(15)%name='RIPP';  an_extra(15)%idim_type=3; an_extra(15)%ind_comp=15
  an_extra(16)%name='RIDP';  an_extra(16)%idim_type=3; an_extra(16)%ind_comp=16
  an_extra(17)%name='RAP';   an_extra(17)%idim_type=3; an_extra(17)%ind_comp=17
  an_extra(18)%name='RGP';   an_extra(18)%idim_type=3; an_extra(18)%ind_comp=18
  an_extra(19)%name='RHP';   an_extra(19)%idim_type=3; an_extra(19)%ind_comp=19
  an_extra(20)%name='CICP';  an_extra(20)%idim_type=3; an_extra(20)%ind_comp=20
  an_extra(21)%name='CIPP';  an_extra(21)%idim_type=3; an_extra(21)%ind_comp=21
  an_extra(22)%name='CIDP';  an_extra(22)%idim_type=3; an_extra(22)%ind_comp=22
  an_extra(23)%name='CAP';   an_extra(23)%idim_type=3; an_extra(23)%ind_comp=23
  an_extra(24)%name='CGP';   an_extra(24)%idim_type=3; an_extra(24)%ind_comp=24
  an_extra(25)%name='CHP';   an_extra(25)%idim_type=3; an_extra(25)%ind_comp=25
  an_extra(26)%name='RPP';   an_extra(26)%idim_type=3; an_extra(26)%ind_comp=26
  an_extra(27)%name='RSP';   an_extra(27)%idim_type=3; an_extra(27)%ind_comp=27
  an_extra(28)%name='CPP';   an_extra(28)%idim_type=3; an_extra(28)%ind_comp=28
  an_extra(29)%name='CSP';   an_extra(29)%idim_type=3; an_extra(29)%ind_comp=29
  an_extra(30)%name='RIFNP'; an_extra(30)%idim_type=3; an_extra(30)%ind_comp=30
  an_extra(31)%name='CIFNP'; an_extra(31)%idim_type=3; an_extra(31)%ind_comp=31
 endif
endif

return
END SUBROUTINE anal_extra_init

!##############################################################################
Subroutine anal_extra_comp (nv,a,ngr,skip)

! Compute the "extra" analysis file variables. Note 3-d arrays must be returned
!  as (i,j,k)

use mem_basic
use mem_grid
use mem_turb
use node_mod
use mem_micro
use mem_leaf
use micro_prm, only:nkr, krdrop,krpris
use module_hujisbm, only: xl,xi,xs,xg,xh,xccn,iceflag,iceprocs
use micphys, only:ipris,igraup,ihail,iifn
! Compute the "extra" analysis file variables. Note 3-d arrays must be returned
!  as (i,j,k)

implicit none
integer :: nv,ngr,izero
real :: a(*)
logical :: skip

skip = .false.
izero = 0

select case (an_extra(nv)%ind_comp)
   case(1)
      ! Total exner function
      CALL ancomp_pi (mmxyzp(ngr)  &
                 ,basic_g(ngr)%pp,basic_g(ngr)%pi0,a)
   case(2)
      ! Horizontal heat eddy diffusivity
      CALL ancomp_hkh (mmxyzp(ngr),turb_g(ngr)%hkm  &
         ,turb_g(ngr)%vkh,basic_g(ngr)%dn0  &
         ,idiffk(ngr),xkhkm(ngr),a)
   case(3)
      ! Vertical heat eddy diffusivity
      CALL ancomp_vkh (mmxyzp(ngr),turb_g(ngr)%vkh  &
             ,basic_g(ngr)%dn0,a)
   case(4) 
      ! Check; 3D variable here
      skip = .true.
   case(5)
      ! Check; add new here 4D leaf
      skip = .true.
   case(6)
      ! Check; add new here 5D leaf
      skip = .true.
   case(7)
      ! Check; add new here 6D leaf
      skip = .true.
   case(8)
      ! Total cloud water
      CALL sum_bins (micro_g(ngr)%ffcd,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krdrop-1,izero)
   case(9)
      ! Total rain water
      CALL sum_bins (micro_g(ngr)%ffcd,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krdrop,nkr,izero)
   case(10) 
      ! Total cloud number concentration
      CALL sum_bins_conc (micro_g(ngr)%ffcd,a,xl &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krdrop-1,izero)
   case(11) 
      ! Total rain number concentration
      CALL sum_bins_conc (micro_g(ngr)%ffcd,a,xl &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krdrop,nkr,izero)
   case(12) 
      ! Total CCN mass
      CALL sum_bins (micro_g(ngr)%fncn,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
   case(13)
      ! Total CCN number concentration
      CALL sum_bins_conc (micro_g(ngr)%fncn,a,xccn &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
   case(14)
      ! Total ice column mass
      if (ipris == 1 .or. ipris >= 4) then
         CALL sum_bins (micro_g(ngr)%ffic,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif
   case(15)
      ! Total ice plate mass
      if (ipris == 2 .or. ipris >=4) then
         CALL sum_bins (micro_g(ngr)%ffip,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(16)
      ! Total ice dendrite mass
      if (ipris >= 3) then
         CALL sum_bins (micro_g(ngr)%ffid,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(17)
      ! Total aggregate mass
      CALL sum_bins (micro_g(ngr)%ffsn,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
   case(18)
      ! Total graupel mass
      if (igraup > 0) then
         CALL sum_bins (micro_g(ngr)%ffgl,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(19)
      ! Total hail mass
      if (ihail > 0) then
         CALL sum_bins (micro_g(ngr)%ffhl,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(20)
      ! Total ice column number
      if (ipris == 1 .or. ipris >= 4) then
         CALL sum_bins_conc (micro_g(ngr)%ffic,a,xi(:,1) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(21)
      ! Total ice plate number
      if (ipris == 2 .or. ipris >=4) then
         CALL sum_bins_conc (micro_g(ngr)%ffip,a,xi(:,2) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(22)
      ! Total ice dendrite number
      if (ipris >=3) then
         CALL sum_bins_conc (micro_g(ngr)%ffid,a,xi(:,3) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(23)
      ! Total aggregate number
      CALL sum_bins_conc (micro_g(ngr)%ffsn,a,xs &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
   case(24)
      ! Total graupel number
      if (igraup > 0) then
         CALL sum_bins_conc (micro_g(ngr)%ffgl,a,xg &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(25)
      ! Total hail number
      if (ihail > 0) then
         CALL sum_bins_conc (micro_g(ngr)%ffhl,a,xh &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(26)
      ! Total pristine ice mass (Bulk scheme equivalent)
      if (ipris > 0) then
         izero = 0 ! If 0, then initialize a to 0 and set izero = 1
         if (ipris == 1 .or. ipris >= 4) &
            CALL sum_bins (micro_g(ngr)%ffic,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(1)-1,izero)
         if (ipris == 2 .or. ipris >= 4) &
            CALL sum_bins (micro_g(ngr)%ffip,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(2)-1,izero)
         if (ipris >= 3) &
            CALL sum_bins (micro_g(ngr)%ffid,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(3)-1,izero)
      else
         skip = .true.
      endif      
   case(27)
      ! Total snow mass (Bulk scheme equivalent)
      if (ipris > 0) then
         izero = 0 ! If 0, then initialize a to 0 and set izero = 1
         if (ipris == 1 .or. ipris >= 4) &
            CALL sum_bins (micro_g(ngr)%ffic,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(1),nkr,izero)
         if (ipris == 2 .or. ipris >= 4) &
            CALL sum_bins (micro_g(ngr)%ffip,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(2),nkr,izero)
         if (ipris >= 3) &
            CALL sum_bins (micro_g(ngr)%ffid,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(3),nkr,izero)
      else
         skip = .true.
      endif      
   case(28)
      ! Total pristine ice number (Bulk scheme equivalent)
      if (ipris > 0) then
         izero = 0 ! If 0, then initialize a to 0 and set izero = 1
         if (ipris == 1 .or. ipris >= 4) &
            CALL sum_bins_conc (micro_g(ngr)%ffic,a,xi(:,1) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(1)-1,izero)
         if (ipris == 2 .or. ipris >= 4) &
            CALL sum_bins_conc (micro_g(ngr)%ffip,a,xi(:,2) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(2)-1,izero)
         if (ipris >= 3) &
            CALL sum_bins_conc (micro_g(ngr)%ffid,a,xi(:,3) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,krpris(3)-1,izero)
      else
         skip = .true.
      endif      
   case(29)
      ! Total snow number (Bulk scheme equivalent)
      if (ipris > 0) then
         izero = 0 ! If 0, then initialize a to 0 and set izero = 1
         if (ipris == 1 .or. ipris >= 4) &
            CALL sum_bins_conc (micro_g(ngr)%ffic,a,xi(:,1) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(1),nkr,izero)
         if (ipris == 2 .or. ipris >= 4) &
            CALL sum_bins_conc (micro_g(ngr)%ffip,a,xi(:,2) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(2),nkr,izero)
         if (ipris >= 3) &
            CALL sum_bins_conc (micro_g(ngr)%ffid,a,xi(:,3) &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),krpris(3),nkr,izero)
      else
         skip = .true.
      endif      
   case(30)
      ! Total IN mass
      if (iifn == 2 .and. iceprocs == 1) then
         CALL sum_bins (micro_g(ngr)%ffin,a &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      
   case(31)
      ! Total IN number
      if (iifn == 2 .and. iceprocs == 1) then
         CALL sum_bins_conc (micro_g(ngr)%ffin,a,xccn &
            ,mmzp(ngr),mmxp(ngr),mmyp(ngr),1,nkr,izero)
      else
         skip = .true.
      endif      

end select
return
END SUBROUTINE anal_extra_comp

!##############################################################################
Subroutine ancomp_pi (n1,b,c,a)

implicit none

integer :: n1
real, dimension(n1) :: a,b,c
integer :: i
 
do i=1,n1
   a(i)=b(i)+c(i)
enddo

return
END SUBROUTINE ancomp_pi

!##############################################################################
Subroutine ancomp_hkh (n1,hkm,vkh,dn0,idiffk,xkhkm,a)

implicit none

integer :: n1,idiffk
real :: xkhkm
real, dimension(n1) :: hkm,vkh,dn0,a
integer :: ind

! Convert to HKM to HKH (note that VKH is HKH for Deardorff)

if (idiffk <= 3) then
   do ind = 1,n1
      a(ind) = hkm(ind) * xkhkm / dn0(ind)
   enddo
elseif (idiffk == 4) then
   do ind = 1,n1
      a(ind) = vkh(ind) / dn0(ind)
   enddo
endif

return
END SUBROUTINE ancomp_hkh

!##############################################################################
Subroutine ancomp_vkh (n1,vkh,dn0,a)

implicit none

integer :: n1
real :: vkh(n1),dn0(n1),a(n1)
integer :: ind

! Un-density weight VKH

do ind = 1,n1
   a(ind) = vkh(ind) / dn0(ind)
enddo

return
END SUBROUTINE ancomp_vkh

END MODULE anal_extra
