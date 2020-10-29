!##############################################################################
Subroutine tend0 ()

use mem_grid
use mem_tend
use var_tables
use node_mod
use micphys

implicit none

integer :: n

!     This routine simply sets all tendency arrays to zero.

!     First u,v tendencies

CALL azero (mxyzp,tend%ut(1))
CALL azero (mxyzp,tend%vt(1))
CALL azero (mxyzp,tend%wt(1))
CALL azero (mxyzp,tend%pt(1))

!     Now sclrr tendencies

do n = 1,num_scalar(ngrid)
   CALL azero (mxyzp,scalar_tab(n,ngrid)%var_t)
enddo

!******************************************************************************
!Add tendency forcings here for current timestep

!Convergence forcing
if(iprntstmt>=1 .and. print_msg) &
  print*,'Ngrid, Time, Domain Max W: ',ngrid,time,vertvel_max(ngrid)
if(ICONV > 0 .and. ICONGR == ngrid) &
  CALL conv_forcing (tend%ut(1),tend%vt(1))

return
END SUBROUTINE tend0

!##############################################################################
Subroutine hadvance (iac)

use mem_grid, only:ngrid,icorflg,jdim,ngbegun
use node_mod, only:mxyzp
use mem_basic
use mem_scratch

implicit none

integer :: iac

!     It is here that the Asselin filter is applied.  For the velocities
!     and pressure, this must be done in two stages, the first when
!     IAC=1 and the second when IAC=2.

!     For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.

CALL predict (mxyzp,basic_g(ngrid)%uc(1,1,1)   &
   ,basic_g(ngrid)%up(1,1,1),scratch%vt3da(1),iac &
   ,scratch%u_the_d)

if (icorflg .eq. 1 .or. jdim .eq. 1) then
   CALL predict (mxyzp,basic_g(ngrid)%vc(1,1,1)  &
      ,basic_g(ngrid)%vp(1,1,1),scratch%vt3da(1),iac &
      ,scratch%v_the_d)
endif

CALL predict (mxyzp,basic_g(ngrid)%wc(1,1,1),basic_g(ngrid)%wp(1,1,1)  &
   ,scratch%vt3da(1),iac,scratch%w_the_d)
CALL predict (mxyzp,basic_g(ngrid)%pc(1,1,1),basic_g(ngrid)%pp(1,1,1)  &
   ,scratch%vt3da(1),iac,scratch%p_the_d)

return
END SUBROUTINE hadvance

!##############################################################################
Subroutine predict (npts,ac,ap,af,iac,the_d)

use mem_grid, only:ngrid,ngbegun

implicit none

integer :: npts,iac,m,iraw
real :: epsu,the_a
real, dimension(*) :: ac,ap,af,the_d

!     This routine moves the arrays AC and AP forward by
!     1 time level by adding in the prescribed tendency. It also
!     applies the Asselin filter given by:
!
!              {AC} = AC + EPS * (AP - 2 * AC + AF)
!
!     where AP,AC,AF are the past, current and future time levels of A.
!     All IAC=1 does is to perform the {AC} calculation without the AF
!     term present.  IAC=2 completes the calculation of {AC} by adding
!     the AF term only, and advances AC by filling it with input AP
!     values which were already updated in ACOUSTC.

! Flag to turn on/off new Robert-Asselin-Williams (RAW) filter (1=on).
! Paul Williams(2011,MWR) "The RAW Filter: An Improvement to the Robert–Asselin
! Filter in Semi-Implicit Integrations".
! (Saleeby(10-26-2020):Consider this update for iraw=1 experimental at this time)
iraw=0

! Aryeh changed from .2 to .1 for a nu value of 0.2 following RAW notation
! Aryeh subsequently reverted back to 0.2 for a nu value of 0.4
epsu = 0.2
the_a = 0.53

if (ngbegun(ngrid) .eq. 0) epsu = 0.5
if (ngbegun(ngrid) .eq. 0) the_a = 1.

if (iac .eq. 1) then
   do m = 1,npts
    if(iraw==1)then
      the_d(m) = epsu * (ap(m) - 2. * ac(m))
      ac(m) = ac(m) + the_a * the_d(m)
    else
      ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))  ! this = this + nu(last - 2*this)
    endif
   enddo
   return
elseif (iac .eq. 2) then
   do m = 1,npts
      af(m) = ap(m)                                ! next = last
      if(iraw==1)then
       the_d(m) = the_d(m) + epsu * af(m)
       ap(m) = ac(m) + the_a * (epsu * af(m))
      else
       ap(m) = ac(m) + epsu * af(m)                ! last = this + nu*next
      endif
   enddo
endif

do m = 1,npts
  if(iraw==1)then
   ac(m) = af(m) + (the_a - 1.) * the_d(m)
  else
   ac(m) = af(m)                                   ! this = next
  endif
enddo

return
END SUBROUTINE predict

!##############################################################################
Subroutine predtr ()

use mem_grid
use var_tables
use node_mod

implicit none

integer :: n

!   -  Step thermodynamic variables from  t  to  t+1.
!   -  Set top, lateral and bottom boundary conditions on some variables
!        if needed.
!   -  call adjustment to assure all positive definite quantities
!        remain positive.
!   -  Rediagnose some thermodynamic quantities for use on the small
!        timestep.

!     Update the scalars and apply lateral, top, and bottom boundary
!     conditions.

do n = 1,num_scalar(ngrid)
   CALL update (mxyzp,scalar_tab(n,ngrid)%var_p  &
                    ,scalar_tab(n,ngrid)%var_t, dtlt)
enddo

return
END SUBROUTINE predtr

!##############################################################################
Subroutine predthp ()

use mem_grid, only:ngrid,dtlt
use node_mod, only:mxyzp
use mem_basic
use mem_radiate

implicit none


! Step thermodynamic variable from  t  to  t+1 with radiative tendendcy
if(ilwrtyp+iswrtyp > 0) &
 CALL update (mxyzp,basic_g(ngrid)%thp,radiate_g(ngrid)%fthrdp,dtlt)

return
END SUBROUTINE predthp

