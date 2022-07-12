!##############################################################################
Subroutine par_decomp_domain (nxp,nyp,nodes,ixb,ixe,iyb,iye)

implicit none

integer :: nxp,nyp,nodes
integer, dimension(nodes) :: ixb,ixe,iyb,iye

real :: anodes,aslabs
integer :: nslabs,min_blocks,nbigslabs
integer, dimension(:), allocatable :: nblocks
integer :: islab, inode, iblock
integer, dimension(:), allocatable :: temp_ixb,temp_ixe,temp_iyb,temp_iye
real, dimension(:), allocatable :: weights

  ! This routine decomposes grid domains of size (nnxp,nnyp) into a number,
  ! specified by nodes, of rectangular subdomains.  The convention is followed
  ! that any internal boundaries (between subdomains) that are parallel to
  ! the x-axis run continuously across the full domain, while boundaries
  ! parallel to the y-axis may or may not run the full distance across the
  ! domain.  For convenience, regions of the domain bounded by adjacent
  ! east-west internal boundaries are termed "slabs", while smaller divisions
  ! within each slab are termed "blocks".  Each block is required to have
  ! a minimum dimension of 6 by 6 grid cells.  If this cannot be satisfied
  ! with the given input parameters, the routine stops.

  ! Assume that all nodes run at the same speed, and that all columns in
  ! the domain take the same number of compute cycles during simulation. This
  ! boils the problem down to dividing the domain into regions that all have
  ! the same number of columns (ie, the same horizontal area).

  ! Make an intial guess by attempting to chop up the domain into same
  ! sized rectangle that tile into the domain. 
  anodes = float(nodes)
  aslabs = sqrt(anodes * (float(nyp) / float(nxp)))
  nslabs = min(nodes, max(1, nint(aslabs)))

  allocate(nblocks(nslabs), temp_iyb(nslabs), temp_iye(nslabs)) 

  ! The above code won't typically find an exact solution so clean up by
  ! placing leftover blocks, one each, into the bottom slabs. This will
  ! keep the tile sizes close to the same (note that this works better
  ! as the number of nodes goes up).
  min_blocks = nodes / nslabs
  nbigslabs = mod(nodes, nslabs)
  do islab = 1, nslabs
    if (islab .le. nbigslabs) then
      nblocks(islab) = min_blocks + 1
    else
      nblocks(islab) = min_blocks
    endif
  enddo

  ! In order to keep the balance between slabs where we want each block to be
  ! the same area, we want to split the slabs where the ones with more blocks
  ! are taller than the ones with fewer blocks. Do this by biasing the height
  ! of each slab by the count of blocks within each slab divided by the total
  ! number of nodes.
  allocate(weights(nslabs))
  do islab = 1, nslabs
    weights(islab) = float(nblocks(islab)) / float(nodes)
  enddo

  ! First figure out the spacing of the slabs (along the y dimension). Every 
  ! sub-domain in a given slab will have the same iyb and iye values so just 
  ! compute them once and record them to be copied into the actual nodes that 
  ! get assigned to each slab.
  CALL calc_ibeg_iend_values (nyp,nslabs,weights,temp_iyb,temp_iye)
  deallocate(weights)

  ! Now walk through each slab and assign the x beginning and end values to 
  ! each block within a given slab. Start with node 1 in the "southwest" 
  ! corner, increment the node number as you go across (west to east) each 
  ! slab. Walk through the slabs from south to north.
  !
  ! For each slab, we want to bias the width of the blocks by 1 divided by 
  ! the number of blocks in the slab.
  inode = 0
  do islab = 1, nslabs
    allocate(temp_ixb(nblocks(islab)),temp_ixe(nblocks(islab)) &
            ,weights(nblocks(islab)))

    do iblock = 1, nblocks(islab)
      weights(iblock) = 1.0 / float(nblocks(islab))
    enddo

    CALL calc_ibeg_iend_values (nxp,nblocks(islab),weights,temp_ixb,temp_ixe)

    do iblock = 1, nblocks(islab)
      inode = inode + 1
      ixb(inode) = temp_ixb(iblock)
      ixe(inode) = temp_ixe(iblock)
      iyb(inode) = temp_iyb(islab)
      iye(inode) = temp_iye(islab)
    enddo
    deallocate(temp_ixb,temp_ixe,weights)
  enddo
  
  deallocate(nblocks,temp_iyb,temp_iye)

  ! Check to make sure that each subdomain has at least 2 interior 
  ! rows and columns. Don't do this check if there is only one node
  ! in order to allow a single node run to do a 2D simulation.
  if(nodes .gt. 1) then
   do inode = 1,nodes
    if(iye(inode)-iyb(inode) .lt. 1 .or. ixe(inode)-ixb(inode) .lt. 1) then
     print*,'grid:',nxp,nyp,'  subdomain too small on node ',inode
     print*,'(ixb,ixe,iyb,iye) = ',ixb(inode),ixe(inode),iyb(inode),iye(inode)
     stop 'small_nodes'
    endif
   enddo
  endif

return
END SUBROUTINE par_decomp_domain

!##############################################################################
Subroutine calc_ibeg_iend_values (npts,nintervals,weights,ibeg,iend)

implicit none

integer :: npts, nintervals
integer, dimension(nintervals) :: ibeg, iend
real, dimension(nintervals) :: weights

integer :: i
real, dimension(:), allocatable :: dist, ends

  ! The decomposition does not include the boundaries of the domain meaning 
  ! that the interval across the entire axis goes from 2 to npts-1.
  !
  !    dist(1) = ((npts-1) - 2) * weights(1)  = (npts-3) * weights(1)
  !    dist(2) = (npts-3) * weights(2)
  !    dist(3) = (npts-3) * weights(3)
  !    ...
  !
  ! Then the endpoint of each interval:
  !
  !    end(1) = 2 + dist
  !    end(2) = end(1) + dist
  !    end(3) = end(2) + dist
  !    ...
  !
  ! and the beginning point of each interval:
  !
  !    beg(1) = 2
  !    beg(2) = end(1) + 1
  !    beg(3) = end(2) + 1
  !    ...
  !
  allocate(ends(nintervals), dist(nintervals))

  do i = 1, nintervals
    dist(i) = float(npts-3) * weights(i)
  enddo

  ibeg(1) = 2
  ends(1) = 2.0 + dist(1)
  iend(1) = nint(ends(1))
  do i = 2, nintervals
    ends(i) = ends(i-1) + dist(i)
    iend(i) = nint(ends(i))
    ibeg(i) = iend(i-1) + 1
  enddo

  ! Just to be safe, force the final end point to match npts-1
  iend(nintervals) = npts - 1

  deallocate(ends)

  return
END SUBROUTINE calc_ibeg_iend_values

!##############################################################################
Subroutine par_node_paths (maxgrds,ngrids,nxpmax,nypmax,nnxp,nnyp  &
   ,nxtnest,mynode,maxmach,nodes,ipm,jpm,ixb,ixe,iyb,iye,ipaths &
   ,igetpaths,ibnd,jbnd,not_a_node)

use cyclic_mod

implicit none

integer :: maxgrds,ngrids,nxpmax,nypmax,mynode,maxmach,nodes,not_a_node

integer :: nnxp(maxgrds),nnyp(maxgrds),nxtnest(maxgrds) &
     ,ipm(nxpmax,maxgrds),jpm(nypmax,maxgrds) &
     ,ixb(maxmach,maxgrds),ixe(maxmach,maxgrds),iyb(maxmach,maxgrds) &
     ,iye(maxmach,maxgrds),ipaths(5,7,maxgrds,maxmach) &
     ,igetpaths(6,maxgrds,maxmach)

integer :: ngr,isend_type,idn,isn,i,j,ibnd,jbnd,info,nnn &
        ,id,jd,is,js,is_t(2),js_t(2),is_u(2),js_u(2),is_v(2),js_v(2) &
        ,indt(2),indu(2),indv(2),nxp,nyp,avgflg_t(2),avgflg_u(2),avgflg_v(2) &
        ,dir,i1d,i2d,j1d,j2d,i1s,i2s,j1s,j2s

! if using cyclic boundary conditions, allocate ipaths_cyc array

if (ibnd .eq. 2 .or. jbnd .eq. 2) then
   CALL ipaths_cyc_alloc (nnxp(1),nnyp(1),ibnd,jbnd,maxmach)
endif

! Put in zeros
! Set entries representing node numbers with "not_a_node"
do ngr=1,ngrids
   do isend_type = 1,6
      do i = 1,nodes
         igetpaths(isend_type,ngr,i) = not_a_node
         do info = 1,4
            ipaths(info,isend_type,ngr,i) = 0
         enddo
         ipaths(5,isend_type,ngr,i) = not_a_node
      enddo
   enddo
enddo


do ngr=1,ngrids
   nnn = nxtnest(ngr)

   ! Fill in the ipaths array which represents this node (mynode) 
   ! acting as the source node
   do idn = 1,nodes
      if (idn .ne. mynode) then
         ! if there is overlap between mynode and idn then record the path
         if(ixb(idn,ngr)-1 .le. ixe(mynode,ngr) .and.  &
            iyb(idn,ngr)-1 .le. iye(mynode,ngr) .and.  &
            ixe(idn,ngr)+1 .ge. ixb(mynode,ngr) .and.  &
            iye(idn,ngr)+1 .ge. iyb(mynode,ngr)) then
   
            ipaths(1,1,ngr,idn)=max(ixb(mynode,ngr),ixb(idn,ngr)-1)
            ipaths(2,1,ngr,idn)=min(ixe(mynode,ngr),ixe(idn,ngr)+1)
            ipaths(3,1,ngr,idn)=max(iyb(mynode,ngr),iyb(idn,ngr)-1)
            ipaths(4,1,ngr,idn)=min(iye(mynode,ngr),iye(idn,ngr)+1)
            ipaths(5,1,ngr,idn)=idn
   
            ! Expand ipaths to include [coarse] grid extern boundary points.
   
            if (ipaths(1,1,ngr,idn) .eq. 2)  &
                ipaths(1,1,ngr,idn) = 1
            if (ipaths(2,1,ngr,idn) .eq. nnxp(ngr)-1)  &
                ipaths(2,1,ngr,idn) = nnxp(ngr)
            if (ipaths(3,1,ngr,idn) .eq. 2)  &
                ipaths(3,1,ngr,idn) = 1
            if (ipaths(4,1,ngr,idn) .eq. nnyp(ngr)-1)  &
                ipaths(4,1,ngr,idn) = nnyp(ngr)
   
            ! Small timestep overlap regions for u
   
            if (ixb(idn,ngr)-1 .eq. ixe(mynode,ngr) .and.  &
                iye(idn,ngr)   .ge. iyb(mynode,ngr) .and.  &
                iyb(idn,ngr)   .le. iye(mynode,ngr)) then
   
               ipaths(1,2,ngr,idn)=ixe(mynode,ngr)
               ipaths(2,2,ngr,idn)=ixe(mynode,ngr)
               ipaths(3,2,ngr,idn)=max(iyb(mynode,ngr),iyb(idn,ngr))
               ipaths(4,2,ngr,idn)=min(iye(mynode,ngr),iye(idn,ngr))
               ipaths(5,2,ngr,idn)=idn
            endif
   
            ! Small timestep overlap regions for v
   
            if (iyb(idn,ngr)-1 .eq. iye(mynode,ngr) .and.  &
                ixe(idn,ngr)   .ge. ixb(mynode,ngr) .and.  &
                ixb(idn,ngr)   .le. ixe(mynode,ngr)) then
   
               ipaths(1,3,ngr,idn)=max(ixb(mynode,ngr),ixb(idn,ngr))
               ipaths(2,3,ngr,idn)=min(ixe(mynode,ngr),ixe(idn,ngr))
               ipaths(3,3,ngr,idn)=iye(mynode,ngr)
               ipaths(4,3,ngr,idn)=iye(mynode,ngr)
               ipaths(5,3,ngr,idn)=idn
            endif
   
            ! Small timestep overlap regions for pi'
   
            if (ixe(idn,ngr)+1 .eq. ixb(mynode,ngr) .and.  &
                iye(idn,ngr)   .ge. iyb(mynode,ngr) .and.  &
                iyb(idn,ngr)   .le. iye(mynode,ngr)) then
   
               ipaths(1,4,ngr,idn)=ixb(mynode,ngr)
               ipaths(2,4,ngr,idn)=ixb(mynode,ngr)
               ipaths(3,4,ngr,idn)=max(iyb(mynode,ngr),iyb(idn,ngr))
               ipaths(4,4,ngr,idn)=min(iye(mynode,ngr),iye(idn,ngr))
               ipaths(5,4,ngr,idn)=idn
   
            elseif (iye(idn,ngr)+1 .eq. iyb(mynode,ngr) .and.  &
                    ixe(idn,ngr)   .ge. ixb(mynode,ngr) .and.  &
                    ixb(idn,ngr)   .le. ixe(mynode,ngr)) then
   
               ipaths(1,4,ngr,idn)=max(ixb(mynode,ngr),ixb(idn,ngr))
               ipaths(2,4,ngr,idn)=min(ixe(mynode,ngr),ixe(idn,ngr))
               ipaths(3,4,ngr,idn)=iyb(mynode,ngr)
               ipaths(4,4,ngr,idn)=iyb(mynode,ngr)
               ipaths(5,4,ngr,idn)=idn
            endif
         endif
      endif

      ! If this is a grid that belongs to a parent grid, then
      ! record paths for doing the nesting interpolation and feedback.
      if (nnn .gt. 0) then

         ! Coarse grid to fine grid interpolation communication

         if(ipm(ixb(idn,ngr)-1,ngr)-2 .le. ixe(mynode,nnn).and.  &
            jpm(iyb(idn,ngr)-1,ngr)-2 .le. iye(mynode,nnn).and.  &
            ipm(ixe(idn,ngr)+1,ngr)+1 .ge. ixb(mynode,nnn).and.  &
            jpm(iye(idn,ngr)+1,ngr)+1 .ge. iyb(mynode,nnn)) then

            ipaths(1,5,ngr,idn) = max(ixb(mynode,nnn), ipm(ixb(idn,ngr)-1,ngr)-2)
            ipaths(2,5,ngr,idn) = min(ixe(mynode,nnn), ipm(ixe(idn,ngr)+1,ngr)+1)
            ipaths(3,5,ngr,idn) = max(iyb(mynode,nnn), jpm(iyb(idn,ngr)-1,ngr)-2)
            ipaths(4,5,ngr,idn) = min(iye(mynode,nnn), jpm(iye(idn,ngr)+1,ngr)+1)
            ipaths(5,5,ngr,idn) = idn
         endif

         ! Fine grid to coarse grid averaging communication

         if(ixb(idn,nnn)-1 .le. ipm(ixe(mynode,ngr),ngr).and.  &
            iyb(idn,nnn)-1 .le. jpm(iye(mynode,ngr),ngr).and.  &
            ixe(idn,nnn)+1 .ge. ipm(ixb(mynode,ngr),ngr).and.  &
            iye(idn,nnn)+1 .ge. jpm(iyb(mynode,ngr),ngr)) then

            ipaths(1,6,ngr,idn) = max(ipm(ixb(mynode,ngr),ngr), ixb(idn,nnn)-1)
            ipaths(2,6,ngr,idn) = min(ipm(ixe(mynode,ngr),ngr), ixe(idn,nnn)+1)
            ipaths(3,6,ngr,idn) = max(jpm(iyb(mynode,ngr),ngr), iyb(idn,nnn)-1)
            ipaths(4,6,ngr,idn) = min(jpm(iye(mynode,ngr),ngr), iye(idn,nnn)+1)
            ipaths(5,6,ngr,idn) = idn

            ! A second index value of 7 of the ipaths array is used to determine
            ! the loop limits in fdbackp for averaging the fm over the overlap
            ! between the cm node and fm node, rather than always over the full
            ! fm node.  It is not used for actually sending stuff.  The
            ! ipaths(*,6,*,*,*) part of the array is still used for sending the
            ! block of averaged cm points from the fm node to the cm node.

            do i = ixb(mynode,ngr),ixe(mynode,ngr)
               if (ipm(i,ngr) .eq. ipaths(1,6,ngr,idn)) then
                  ipaths(1,7,ngr,idn) = i
                  exit
               endif
            enddo

            do i = ixe(mynode,ngr),ixb(mynode,ngr),-1
               if (ipm(i,ngr) .eq. ipaths(2,6,ngr,idn)) then
                  ipaths(2,7,ngr,idn) = i
                  exit
               endif
            enddo

            do j = iyb(mynode,ngr),iye(mynode,ngr)
               if (jpm(j,ngr) .eq. ipaths(3,6,ngr,idn)) then
                  ipaths(3,7,ngr,idn) = j
                  exit
               endif
            enddo

            do j = iye(mynode,ngr),iyb(mynode,ngr),-1
               if (jpm(j,ngr) .eq. ipaths(4,6,ngr,idn)) then
                  ipaths(4,7,ngr,idn) = j
                  exit
               endif
            enddo 
         endif
      endif
   enddo

   ! Fill in the igetpaths array - this node (mynode) acting 
   ! as the destination node
   do isn = 1,nodes
      if (isn .ne. mynode) then
         ! if there is overlap between mynode and isn then record the paths
         if(ixb(mynode,ngr)-1 .le. ixe(isn,ngr) .and.  &
            iyb(mynode,ngr)-1 .le. iye(isn,ngr) .and.  &
            ixe(mynode,ngr)+1 .ge. ixb(isn,ngr) .and.  &
            iye(mynode,ngr)+1 .ge. iyb(isn,ngr)) then
   
            igetpaths(1,ngr,isn) = isn
   
            ! Small timestep overlap regions for u
   
            if (ixb(mynode,ngr)-1 .eq. ixe(isn,ngr) .and.  &
                iye(mynode,ngr)   .ge. iyb(isn,ngr) .and.  &
                iyb(mynode,ngr)   .le. iye(isn,ngr)) then
   
               igetpaths(2,ngr,isn) = isn
            endif
   
            ! Small timestep overlap regions for v
   
            if (iyb(mynode,ngr)-1 .eq. iye(isn,ngr) .and.  &
                ixe(mynode,ngr)   .ge. ixb(isn,ngr) .and.  &
                ixb(mynode,ngr)   .le. ixe(isn,ngr)) then
   
               igetpaths(3,ngr,isn) = isn
            endif
   
            ! Small timestep overlap regions for pi'
   
            if (ixe(mynode,ngr)+1 .eq. ixb(isn,ngr) .and.  &
                iye(mynode,ngr)   .ge. iyb(isn,ngr) .and.  &
                iyb(mynode,ngr)   .le. iye(isn,ngr)) then
   
               igetpaths(4,ngr,isn) = isn
   
            elseif (iye(mynode,ngr)+1 .eq. iyb(isn,ngr) .and.  &
                    ixe(mynode,ngr)   .ge. ixb(isn,ngr) .and.  &
                    ixb(mynode,ngr)   .le. ixe(isn,ngr)) then
   
               igetpaths(4,ngr,isn) = isn
            endif
         endif
      endif

      ! If this is a grid that belongs to a parent grid, then
      ! record paths for doing the nesting interpolation and feedback.
      if (nnn .gt. 0) then

         ! Coarse grid to fine grid interpolation communication

         if(ipm(ixb(mynode,ngr)-1,ngr)-2 .le. ixe(isn,nnn).and.  &
            jpm(iyb(mynode,ngr)-1,ngr)-2 .le. iye(isn,nnn).and.  &
            ipm(ixe(mynode,ngr)+1,ngr)+1 .ge. ixb(isn,nnn).and.  &
            jpm(iye(mynode,ngr)+1,ngr)+1 .ge. iyb(isn,nnn)) then

            igetpaths(5,ngr,isn) = isn
         endif

         ! Fine grid to coarse grid averaging communication

         if(ixb(mynode,nnn)-1 .le. ipm(ixe(isn,ngr),ngr).and.  &
            iyb(mynode,nnn)-1 .le. jpm(iye(isn,ngr),ngr).and.  &
            ixe(mynode,nnn)+1.ge. ipm(ixb(isn,ngr),ngr).and.  &
            iye(mynode,nnn)+1.ge. jpm(iyb(isn,ngr),ngr)) then

            igetpaths(6,ngr,isn) = isn
         endif
      endif
   enddo
enddo

! Cyclic boundary conditions (grid 1 only)

indt = 0
indu = 0
indv = 0

nxp = nnxp(1)
nyp = nnyp(1)

! Loop over all destination ij columns

do jd = 1,nyp
   do id = 1,nxp
      if (id <= 2 .or. id >= nxp-1 .or. jd <= 2 .or. jd >= nyp-1) then

         is_t = 0
         js_t = 0
         is_u = 0
         js_u = 0
         is_v = 0
         js_v = 0
         avgflg_t = 0
         avgflg_u = 0
         avgflg_v = 0

         ! Record source i,j for t, u, and v in the i direction
         ! (element 1 of is,js arrays)
         if (ibnd == 2) then
! t
            if (id <= 2) then
               is_t(1) = id+nxp-3
               js_t(1) = jd
            elseif (id >= nxp-1) then
               is_t(1) = id-nxp+3
               js_t(1) = jd
            endif
! u
            if (id == 1) then
               is_u(1) = id+nxp-3
               js_u(1) = jd
            elseif (id == nxp-1) then
               is_u(1) = id-nxp+3
               js_u(1) = jd
            endif
! v
            if (id <= 2) then
               is_v(1) = id+nxp-3
               js_v(1) = jd
            elseif (id >= nxp-1) then
               is_v(1) = id-nxp+3
               js_v(1) = jd
            endif

         endif
          
         ! Record source i,j for t, u, and v in the j direction
         ! (element 2 of is,js arrays)
         if (jbnd == 2) then
!t
            if (jd <= 2) then
               is_t(2) = id
               js_t(2) = jd+nyp-3
            elseif (jd >= nyp-1) then
               is_t(2) = id
               js_t(2) = jd-nyp+3
            endif
!u
            if (jd <= 2) then
               is_u(2) = id
               js_u(2) = jd+nyp-3
            elseif (jd >= nyp-1) then
               is_u(2) = id
               js_u(2) = jd-nyp+3
            endif
!v
            if (jd == 1) then
               is_v(2) = id
               js_v(2) = jd+nyp-3
            elseif (jd == nyp-1) then
               is_v(2) = id
               js_v(2) = jd-nyp+3
            endif

         endif

         ! Do averaging if on the "2" or "n-1" points
         if ((id .eq. 2) .or. (id .eq. nxp-1)) then
            avgflg_t(1) = 1
            avgflg_v(1) = 1
         endif
         if ((jd .eq. 2) .or. (jd .eq. nyp-1)) then
            avgflg_t(2) = 1
            avgflg_u(2) = 1
         endif

! Loop over all destination nodes and check if current ij destination column 
! is in it

         do idn = 1,nodes
            ! grab end points (compute columns)
            i1d = ixb(idn,1)
            i2d = ixe(idn,1)
            j1d = iyb(idn,1)
            j2d = iye(idn,1)
            ! include domain boundaries
            if (i1d .eq. 2) i1d = 1
            if (i2d .eq. nxp-1) i2d = nxp
            if (j1d .eq. 2) j1d = 1
            if (j2d .eq. nyp-1) j2d = nyp

            if (id >= i1d .and. id <= i2d .and. jd >= j1d .and. jd <= j2d) then
                
! Loop over all source nodes and check if current ij source column is in it
! The t paths for each corresponding direction will always contain source 
! nodes, however the u and v paths may not. (Eg, dir == 1, u paths only exist
! for id == 1 and nxp-1 (and not for id == 2 and nxp).

               do dir = 1,2
                  ! dir == 1 --> along the x-axis ("i" direction)
                  ! dir == 2 --> along the y-axis ("j" direction)
                  is = is_t(dir)
                  js = js_t(dir)

                  do isn = 1,nodes
                    ! grab end points (compute columns)
                    i1s = ixb(isn,1)
                    i2s = ixe(isn,1)
                    j1s = iyb(isn,1)
                    j2s = iye(isn,1)
                    ! include domain boundaries
                    if (i1s .eq. 2) i1s = 1
                    if (i2s .eq. nxp-1) i2s = nxp
                    if (j1s .eq. 2) j1s = 1
                    if (j2s .eq. nyp-1) j2s = nyp

          if (is >= i1s .and. is <= i2s .and. js >= j1s .and. js <= j2s) then

           indt(dir) = indt(dir) + 1
           ipathst_cyc(1,indt(dir),dir) = isn            ! source node #
           ipathst_cyc(2,indt(dir),dir) = is_t(dir)      ! source node i
           ipathst_cyc(3,indt(dir),dir) = js_t(dir)      ! source node j
           ipathst_cyc(4,indt(dir),dir) = idn            ! destination node #
           ipathst_cyc(5,indt(dir),dir) = id             ! destination node i
           ipathst_cyc(6,indt(dir),dir) = jd             ! destination node j
           ipathst_cyc(7,indt(dir),dir) = avgflg_t(dir)  ! average flag

           if (is_u(dir) .gt. 0) then
             indu(dir) = indu(dir) + 1
             ipathsu_cyc(1,indu(dir),dir) = isn            ! source node #
             ipathsu_cyc(2,indu(dir),dir) = is_u(dir)      ! source node i
             ipathsu_cyc(3,indu(dir),dir) = js_u(dir)      ! source node j
             ipathsu_cyc(4,indu(dir),dir) = idn            ! destination node #
             ipathsu_cyc(5,indu(dir),dir) = id             ! destination node i
             ipathsu_cyc(6,indu(dir),dir) = jd             ! destination node j
             ipathsu_cyc(7,indu(dir),dir) = avgflg_u(dir)  ! average flag
           endif
   
           if (is_v(dir) .gt. 0) then
             indv(dir) = indv(dir) + 1
             ipathsv_cyc(1,indv(dir),dir) = isn            ! source node #
             ipathsv_cyc(2,indv(dir),dir) = is_v(dir)      ! source node i
             ipathsv_cyc(3,indv(dir),dir) = js_v(dir)      ! source node j
             ipathsv_cyc(4,indv(dir),dir) = idn            ! destination node #
             ipathsv_cyc(5,indv(dir),dir) = id             ! destination node i
             ipathsv_cyc(6,indv(dir),dir) = jd             ! destination node j
             ipathsv_cyc(7,indv(dir),dir) = avgflg_v(dir)  ! average flag
           endif

           exit ! since we found the source node
          endif
                  enddo ! source nodes

               enddo  ! dir = 1, 2

               exit ! since we found the destination node
            endif
         enddo ! destination nodes

      endif ! id, jd along boundaries
   enddo ! id
enddo ! jd

return
END SUBROUTINE par_node_paths
