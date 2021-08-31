!##############################################################################
Subroutine ipaths_cyc_alloc (nxp,nyp,ibnd,jbnd,maxmach)

use cyclic_mod

implicit none

integer :: nxp,nyp,ibnd,jbnd,maxmach

npts_cyc = 4 * maxmach
if (ibnd .eq. 2) npts_cyc = npts_cyc + nyp * 4
if (jbnd .eq. 2) npts_cyc = npts_cyc + nxp * 4

! 7 coded entries (isn,is,js,idn,id,jd,avgflg)
! by npts by 2 directions
allocate (ipathst_cyc(7,npts_cyc,2))
allocate (ipathsu_cyc(7,npts_cyc,2))
allocate (ipathsv_cyc(7,npts_cyc,2))

! 7 coded entries * 2 directions * npts
CALL azero_int (14*npts_cyc,ipathst_cyc(1,1,1))
CALL azero_int (14*npts_cyc,ipathsu_cyc(1,1,1))
CALL azero_int (14*npts_cyc,ipathsv_cyc(1,1,1))

return
END SUBROUTINE ipaths_cyc_alloc

!##############################################################################
Subroutine node_cycinit (nzp,npvar,nmachs,my_rams_num)

use cyclic_mod
use micro_prm, only:nkr

implicit none

integer :: nmachs,my_rams_num,nzp,icycpts,mdn,msn,ndn,nsn,npvar,nm
integer, save, allocatable :: ijcount(:)
integer :: maxijrecv_cyc, maxijsendt_cyc, maxijsendu_cyc, maxijsendv_cyc, &
           maxijrecvt_cyc, maxijrecvu_cyc, maxijrecvv_cyc, idir, &
           max_ndns_cyc, max_nsns_cyc

! This routine is called by each node process.  It uses information from
! the master cyclic parallel array [ipaths_cyc] to construct integer scalars 
! and arrays that the node needs for parallel data sends and receives for 
! cyclic boundary conditions

! ldgrant(2012): This routine is called within the time loop in rams_node
!   (rams_node calls init_fields if sending information back to the master,
!    and init_fields calls this routine)
! Thus added deallocate statements; this fixes the issue where memory usage
! increases through the run when cyclic boundary conditions are used,
! sometimes causing the node to crash.

! if ndn_cyc is allocated, the rest are as well
if ( allocated(ndn_cyc) ) then
   deallocate (ndn_cyc,nsn_cyc,msn_cyc,mdn_cyc  &
      ,nijsendt_cyc,nijsendu_cyc,nijsendv_cyc  &
      ,nijrecvt_cyc,nijrecvu_cyc,nijrecvv_cyc  &
      ,ijcount,isend_req_cyc,irecv_req_cyc)
endif

allocate (ndn_cyc(nmachs,2),nsn_cyc(nmachs,2),msn_cyc(nmachs,2),mdn_cyc(nmachs,2)  &
   ,nijsendt_cyc(nmachs,2),nijsendu_cyc(nmachs,2),nijsendv_cyc(nmachs,2)  &
   ,nijrecvt_cyc(nmachs,2),nijrecvu_cyc(nmachs,2),nijrecvv_cyc(nmachs,2)  &
   ,ijcount(nmachs)  &
   ,isend_req_cyc(6,nmachs,2),irecv_req_cyc(6,nmachs,2))

do nm = 1,nmachs
   nijsendt_cyc(nm,1:2) = 0
   nijsendu_cyc(nm,1:2) = 0
   nijsendv_cyc(nm,1:2) = 0
   nijrecvt_cyc(nm,1:2) = 0
   nijrecvu_cyc(nm,1:2) = 0
   nijrecvv_cyc(nm,1:2) = 0
   isend_req_cyc(1:6,nm,1:2)=-999
   irecv_req_cyc(1:6,nm,1:2)=-999
enddo

do idir = 1, 2
  do icycpts = 1,npts_cyc
    msn = ipathst_cyc(1,icycpts,idir)
    mdn = ipathst_cyc(4,icycpts,idir)
 
    if (msn .eq. my_rams_num .and. mdn .gt. 0) then
       nijsendt_cyc(mdn,idir) = nijsendt_cyc(mdn,idir) + 1
    endif
    
    if (mdn .eq. my_rams_num .and. msn .gt. 0) then
       nijrecvt_cyc(msn,idir) = nijrecvt_cyc(msn,idir) + 1
    endif
  enddo

  do icycpts = 1,npts_cyc
    msn = ipathsu_cyc(1,icycpts,idir)
    mdn = ipathsu_cyc(4,icycpts,idir)
 
    if (msn .eq. my_rams_num .and. mdn .gt. 0) then
       nijsendu_cyc(mdn,idir) = nijsendu_cyc(mdn,idir) + 1
    endif
    if (mdn .eq. my_rams_num .and. msn .gt. 0) then
       nijrecvu_cyc(msn,idir) = nijrecvu_cyc(msn,idir) + 1
    endif
  enddo

  do icycpts = 1,npts_cyc
    msn = ipathsv_cyc(1,icycpts,idir)
    mdn = ipathsv_cyc(4,icycpts,idir)

    if (msn .eq. my_rams_num .and. mdn .gt. 0) then
       nijsendv_cyc(mdn,idir) = nijsendv_cyc(mdn,idir) + 1
    endif
    if (mdn .eq. my_rams_num .and. msn .gt. 0) then
       nijrecvv_cyc(msn,idir) = nijrecvv_cyc(msn,idir) + 1
    endif
  enddo
enddo

! figure out sizes for send lists
maxijsendt_cyc = 0
maxijsendu_cyc = 0
maxijsendv_cyc = 0
maxijrecvt_cyc = 0
maxijrecvu_cyc = 0
maxijrecvv_cyc = 0

max_ndns_cyc = 0
max_nsns_cyc = 0

do idir = 1, 2
  nsn = 0
  ndn = 0

  do nm = 1,nmachs
    if (nijsendt_cyc(nm,idir) .gt. 0) then
      ndn = ndn + 1
      ndn_cyc(nm,idir) = ndn
      mdn_cyc(ndn,idir) = nm
      nijsendt_cyc(ndn,idir) = nijsendt_cyc(nm,idir)
      nijsendu_cyc(ndn,idir) = nijsendu_cyc(nm,idir)
      nijsendv_cyc(ndn,idir) = nijsendv_cyc(nm,idir)
      maxijsendt_cyc = max(maxijsendt_cyc,nijsendt_cyc(ndn,idir))
      maxijsendu_cyc = max(maxijsendu_cyc,nijsendu_cyc(ndn,idir))
      maxijsendv_cyc = max(maxijsendv_cyc,nijsendv_cyc(ndn,idir))
    endif

    if (nijrecvt_cyc(nm,idir) .gt. 0) then
      nsn = nsn + 1
      nsn_cyc(nm,idir) = nsn
      msn_cyc(nsn,idir) = nm
      nijrecvt_cyc(nsn,idir) = nijrecvt_cyc(nm,idir)
      nijrecvu_cyc(nsn,idir) = nijrecvu_cyc(nm,idir)
      nijrecvv_cyc(nsn,idir) = nijrecvv_cyc(nm,idir)
      maxijrecvt_cyc = max(maxijrecvt_cyc,nijrecvt_cyc(nsn,idir))
      maxijrecvu_cyc = max(maxijrecvu_cyc,nijrecvu_cyc(nsn,idir))
      maxijrecvv_cyc = max(maxijrecvv_cyc,nijrecvv_cyc(nsn,idir))
    endif
  enddo

  ndns_cyc(idir) = ndn
  nsns_cyc(idir) = nsn
  max_ndns_cyc = max(max_ndns_cyc, ndns_cyc(idir))
  max_nsns_cyc = max(max_nsns_cyc, nsns_cyc(idir))

enddo

maxijrecv_cyc = max(maxijrecvt_cyc,maxijrecvu_cyc,maxijrecvv_cyc)

! Allocate buffers for send lists
! ldgrant(2012): added deallocate statements; see comment at top of this routine
if ( allocated(ijsendt_cyc) ) deallocate(ijsendt_cyc)
if ( allocated(ijsendu_cyc) ) deallocate(ijsendu_cyc)
if ( allocated(ijsendv_cyc) ) deallocate(ijsendv_cyc)
if ( allocated(ijrecv_cyc) )  deallocate(ijrecv_cyc)

allocate (ijsendt_cyc(5,maxijsendt_cyc,max_ndns_cyc,2) &
         ,ijsendu_cyc(5,maxijsendu_cyc,max_ndns_cyc,2) &
         ,ijsendv_cyc(5,maxijsendv_cyc,max_ndns_cyc,2) &
         ,ijrecv_cyc(5,maxijrecv_cyc))

! Fill in the send lists
do idir = 1,2
  do ndn = 1,ndns_cyc(idir)
    ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
    msn = ipathst_cyc(1,icycpts,idir)
    mdn = ipathst_cyc(4,icycpts,idir)
    if (msn .eq. my_rams_num .and. mdn .gt. 0) then
      ndn = ndn_cyc(mdn,idir)
      ijcount(ndn) = ijcount(ndn) + 1
      ijsendt_cyc(1,ijcount(ndn),ndn,idir) = ipathst_cyc(2,icycpts,idir)
      ijsendt_cyc(2,ijcount(ndn),ndn,idir) = ipathst_cyc(3,icycpts,idir)
      ijsendt_cyc(3,ijcount(ndn),ndn,idir) = ipathst_cyc(5,icycpts,idir)
      ijsendt_cyc(4,ijcount(ndn),ndn,idir) = ipathst_cyc(6,icycpts,idir)
      ijsendt_cyc(5,ijcount(ndn),ndn,idir) = ipathst_cyc(7,icycpts,idir)
    endif
  enddo

  do ndn = 1,ndns_cyc(idir)
   ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
   msn = ipathsu_cyc(1,icycpts,idir)
   mdn = ipathsu_cyc(4,icycpts,idir)
   if (msn .eq. my_rams_num .and. mdn .gt. 0) then
      ndn = ndn_cyc(mdn,idir)
      ijcount(ndn) = ijcount(ndn) + 1
      ijsendu_cyc(1,ijcount(ndn),ndn,idir) = ipathsu_cyc(2,icycpts,idir)
      ijsendu_cyc(2,ijcount(ndn),ndn,idir) = ipathsu_cyc(3,icycpts,idir)
      ijsendu_cyc(3,ijcount(ndn),ndn,idir) = ipathsu_cyc(5,icycpts,idir)
      ijsendu_cyc(4,ijcount(ndn),ndn,idir) = ipathsu_cyc(6,icycpts,idir)
      ijsendu_cyc(5,ijcount(ndn),ndn,idir) = ipathsu_cyc(7,icycpts,idir)
   endif
  enddo

  do ndn = 1,ndns_cyc(idir)
   ijcount(ndn) = 0
  enddo

  do icycpts = 1,npts_cyc
   msn = ipathsv_cyc(1,icycpts,idir)
   mdn = ipathsv_cyc(4,icycpts,idir)
   if (msn .eq. my_rams_num .and. mdn .gt. 0) then
      ndn = ndn_cyc(mdn,idir)
      ijcount(ndn) = ijcount(ndn) + 1
      ijsendv_cyc(1,ijcount(ndn),ndn,idir) = ipathsv_cyc(2,icycpts,idir)
      ijsendv_cyc(2,ijcount(ndn),ndn,idir) = ipathsv_cyc(3,icycpts,idir)
      ijsendv_cyc(3,ijcount(ndn),ndn,idir) = ipathsv_cyc(5,icycpts,idir)
      ijsendv_cyc(4,ijcount(ndn),ndn,idir) = ipathsv_cyc(6,icycpts,idir)
      ijsendv_cyc(5,ijcount(ndn),ndn,idir) = ipathsv_cyc(7,icycpts,idir)
   endif
  enddo

enddo ! idir = 1, 2

! allocate cyclic buffers

nbuffsend_cyc = max(maxijsendt_cyc * (max(2,npvar) * nzp * nkr + 6) + 1  &
                   ,maxijsendu_cyc * (nzp + 6) + 1                 &
                  + maxijsendv_cyc * (nzp + 6) + 1                 ) 

! ldgrant(2012): added deallocate statement
if ( allocated(buffsend_cyc) ) deallocate(buffsend_cyc)
allocate (buffsend_cyc(nbuffsend_cyc,max_ndns_cyc))

nbuffrecv_cyc = max(maxijrecvt_cyc * (max(2,npvar) * nzp *nkr + 6) + 1  &
                   ,maxijrecvu_cyc * (nzp + 6) + 1                 &
                  + maxijrecvv_cyc * (nzp + 6) + 1                 ) 

! ldgrant(2012): added deallocate statement
if ( allocated(buffrecv_cyc) ) deallocate(buffrecv_cyc)
allocate (buffrecv_cyc(nbuffrecv_cyc,max_nsns_cyc))

return
END SUBROUTINE node_cycinit

!##############################################################################
! node_sendcyclic()
!
! This routine along with node_getcyclic() handles the application of cyclic
! boundary conditions upon the domain boundaries. The intention of the parallel
! cyclic boundary application is to mimic the sequential application as seen
! in the routine cyclic(). The order of execution in parallel mode can vary
! according to how the nodes are assigned so the algorithm that exists in cyclic()
! is split up into pieces, for the parallel scheme, that will work independent of
! the order of execution. This code appears in cyclic_para().
!
Subroutine node_sendcyclic (isflag,idir)

use mem_grid
use var_tables
use cyclic_mod
use mem_basic
use node_mod

implicit none

integer :: nsn,icycpts,msn,mdn,isn,jsn,isflag,idir,ndn,nv,mtp
real, pointer :: scalarp
real, dimension(:), allocatable :: cyc_buff

if (ibnd .ne. 2 .and. jbnd .ne. 2) return

!   First, before we send anything, let's post the receives

do nsn = 1,nsns_cyc(idir)
   msn = msn_cyc(nsn,idir)
!                                IN              IN
   CALL par_get_noblock (buffrecv_cyc(1,nsn),nbuffrecv_cyc  &
!                      IN                            IN  
      ,CYCLIC_TAG_BASE+100*machnum(msn)+10*machnum(my_rams_num)+isflag,machnum(msn)  &
!              OUT
      ,irecv_req_cyc(isflag,msn,idir))
enddo

!   Now we can actually go on to sending the stuff
mtp = mmzp(1)
allocate(cyc_buff(mtp))

do ndn = 1,ndns_cyc(idir)

   mdn = mdn_cyc(ndn,idir)
   CALL par_init_put (buffsend_cyc(1,ndn),nbuffsend_cyc)
   
   if (isflag == LBC_ALL_SCALARS .or. isflag == LBC_PP .or. isflag == LBC_WP) then   
      CALL par_put_int (nijsendt_cyc(ndn,idir),1)
      CALL par_put_int (ijsendt_cyc(1,1,ndn,idir),5*nijsendt_cyc(ndn,idir))

      do icycpts = 1,nijsendt_cyc(ndn,idir)

         isn = ijsendt_cyc(1,icycpts,ndn,idir) - mi0(1)
         jsn = ijsendt_cyc(2,icycpts,ndn,idir) - mj0(1)

         if (isflag == LBC_ALL_SCALARS) then
            ! scalar vars (all scalars are 3D)
            ! Even the bin micro variables are 3D
            ! in the scalar table. In the scalar table
            ! each bin of a species distribution is a
            ! separate variable. Adele
            do nv = 1, num_scalar(1)
              scalarp => scalar_tab(nv,1)%var_p
              CALL mkcycbuff (mmzp(1),mmxp(1),mmyp(1)  &
                  ,scalarp,cyc_buff,isn,jsn)
              CALL par_put_float (cyc_buff,mtp)
            enddo
         elseif (isflag == LBC_PP) then
            ! pp
            CALL mkcycbuff (mmzp(1),mmxp(1),mmyp(1)  &
               ,basic_g(1)%pp(1,1,1),cyc_buff  &
               ,isn,jsn)
            CALL par_put_float (cyc_buff,mtp)

         elseif (isflag == LBC_WP) then
            ! wp
            CALL mkcycbuff (mmzp(1),mmxp(1),mmyp(1)   &
               ,basic_g(1)%wp(1,1,1),cyc_buff  &
               ,isn,jsn)
            CALL par_put_float (cyc_buff,mtp)

         endif
      enddo
            
   elseif (isflag == LBC_UP) then   
      ! up
      CALL par_put_int (nijsendu_cyc(ndn,idir),1)
      CALL par_put_int (ijsendu_cyc(1,1,ndn,idir),5*nijsendu_cyc(ndn,idir))

      do icycpts = 1,nijsendu_cyc(ndn,idir)

         isn = ijsendu_cyc(1,icycpts,ndn,idir) - mi0(1)
         jsn = ijsendu_cyc(2,icycpts,ndn,idir) - mj0(1)

         CALL mkcycbuff (mmzp(1),mmxp(1),mmyp(1)  &
            ,basic_g(1)%up(1,1,1),cyc_buff  &
            ,isn,jsn)
         CALL par_put_float (cyc_buff,mtp)

      enddo

   elseif (isflag == LBC_VP) then
      ! vp
      CALL par_put_int (nijsendv_cyc(ndn,idir),1)
      CALL par_put_int (ijsendv_cyc(1,1,ndn,idir),5*nijsendv_cyc(ndn,idir))

      do icycpts = 1,nijsendv_cyc(ndn,idir)

         isn = ijsendv_cyc(1,icycpts,ndn,idir) - mi0(1)
         jsn = ijsendv_cyc(2,icycpts,ndn,idir) - mj0(1)

         CALL mkcycbuff (mmzp(1),mmxp(1),mmyp(1)  &
            ,basic_g(1)%vp(1,1,1),cyc_buff  &
            ,isn,jsn)
         CALL par_put_float (cyc_buff,mtp)

      enddo

   endif

!                         IN          IN
   CALL par_send_noblock (machnum(mdn),CYCLIC_TAG_BASE+100*machnum(my_rams_num)+10*machnum(mdn)+isflag  &
!             OUT
      ,isend_req_cyc(isflag,mdn,idir))
      
enddo

deallocate(cyc_buff)

return
END SUBROUTINE node_sendcyclic

!##############################################################################
Subroutine node_getcyclic (isflag,idir)

use mem_grid
use var_tables
use cyclic_mod
use mem_basic
use node_mod

implicit none

integer :: ndn,nsn,mdn,msn,mtp,ibytes,msgid,ihostnum,isflag,idir,mijrecv,ijr,nv
real, pointer :: scalarp
real, dimension(:), allocatable :: cyc_buff

if (ibnd .ne. 2 .and. jbnd .ne. 2) return

!  First, let's make sure our sends are all finished and de-allocated

do ndn = 1,ndns_cyc(idir)
   mdn = mdn_cyc(ndn,idir)
      
!                    input                   out   out     out
   CALL par_wait (isend_req_cyc(isflag,mdn,idir),ibytes,msgid,ihostnum)
      
enddo

!  Now, let's wait on our receives

mtp = mmzp(1)
allocate(cyc_buff(mtp))

do nsn = 1,nsns_cyc(idir)
   msn = msn_cyc(nsn,idir)
      
!                    input                  out   out     out
   CALL par_wait (irecv_req_cyc(isflag,msn,idir),ibytes,msgid,ihostnum)
      
!  We got all our stuff.  Now unpack it into appropriate space.

   CALL par_assoc_buff (buffrecv_cyc(1,nsn),nbuffrecv_cyc)

   if (isflag == LBC_ALL_SCALARS .or. isflag == LBC_PP .or. isflag == LBC_WP) then   
      CALL par_get_int (mijrecv,1)
      CALL par_get_int (ijrecv_cyc(1,1),5*mijrecv)

      do ijr = 1,mijrecv

         if (isflag == LBC_ALL_SCALARS) then
            ! scalar vars (all scalars are 3D)
            do nv = 1, num_scalar(1)
              scalarp => scalar_tab(nv,1)%var_p
              CALL par_get_float (cyc_buff,mtp)
              CALL cyclic_para (mmzp(1),mmxp(1),mmyp(1)  &
                  ,scalarp,cyc_buff  &
                  ,ijr,mi0(1),mj0(1))
            enddo
         elseif (isflag == LBC_PP) then
            ! pp
            CALL par_get_float (cyc_buff,mtp)
            CALL cyclic_para (mmzp(1),mmxp(1),mmyp(1)  &
               ,basic_g(1)%pp(1,1,1),cyc_buff,ijr,mi0(1),mj0(1))
         elseif (isflag == LBC_WP) then
            ! wp
            CALL par_get_float (cyc_buff,mtp)
            CALL cyclic_para (mmzp(1),mmxp(1),mmyp(1)  &
               ,basic_g(1)%wp(1,1,1),cyc_buff,ijr,mi0(1),mj0(1))
         endif
      enddo
 
   elseif (isflag == LBC_UP) then
      ! up
      CALL par_get_int (mijrecv,1)
      CALL par_get_int (ijrecv_cyc(1,1),5*mijrecv)

      do ijr = 1,mijrecv
         CALL par_get_float (cyc_buff,mtp)
         CALL cyclic_para (mmzp(1),mmxp(1),mmyp(1)  &
            ,basic_g(1)%up(1,1,1),cyc_buff,ijr,mi0(1),mj0(1))
      enddo

   elseif (isflag .eq. LBC_VP) then
      ! vp
      CALL par_get_int (mijrecv,1)
      CALL par_get_int (ijrecv_cyc(1,1),5*mijrecv)

      do ijr = 1,mijrecv
         CALL par_get_float (cyc_buff,mtp)
         CALL cyclic_para (mmzp(1),mmxp(1),mmyp(1)  &
            ,basic_g(1)%vp(1,1,1),cyc_buff,ijr,mi0(1),mj0(1))
      enddo

   endif

enddo

deallocate(cyc_buff)

return
END SUBROUTINE node_getcyclic

!##############################################################################
Subroutine cyclic_para (m1,m2,m3,af,bf,ijr,i0,j0)

use cyclic_mod

implicit none

integer :: m1,m2,m3,ijr,i0,j0
real, dimension(m1,m2,m3) :: af
real, dimension(*) :: bf
integer :: idn,jdn,k,isn,jsn

! The receive buffer contains:
!   1 - source node i
!   2 - source node j
!   3 - destination node i
!   4 - destination node j
!   5 - average flag (0 - do not average, 1 - average)

isn = ijrecv_cyc(1,ijr)
jsn = ijrecv_cyc(2,ijr)
idn = ijrecv_cyc(3,ijr) - i0
jdn = ijrecv_cyc(4,ijr) - j0

do k = 1,m1
   if (ijrecv_cyc(5,ijr) .eq. 1) then
     ! replace with average of existing and received value
     af(k,idn,jdn) = (af(k,idn,jdn) + bf(k)) * 0.5
   else
     ! replace with received value
     af(k,idn,jdn) =  bf(k)
   endif
enddo

return
END SUBROUTINE cyclic_para

!##############################################################################
Subroutine mkcycbuff (n1,n2,n3,a,b,i,j)

implicit none

integer :: k,i,j,n1,n2,n3
real, dimension(n1,n2,n3) :: a
real, dimension(*) :: b


do k = 1,n1
   b(k) = a(k,i,j)
enddo

return
END SUBROUTINE mkcycbuff

