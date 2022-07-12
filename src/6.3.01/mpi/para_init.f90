!##############################################################################
Subroutine assign_node_subdomain ()

use mem_grid
use node_mod
use micphys, only:iscmx,iscmy,iscm

implicit none

  integer :: ngr
  integer :: jnode, ncols

  ! Decompose all grids into subdomains

  do ngr = 1,ngrids
    ! Decompose the domain into rectangular sub-domains that have horizontal
    ! areas (ie, number of columns) match as closely as possible.
    CALL par_decomp_domain (nnxp(ngr),nnyp(ngr),nmachs,ixb(1,ngr),ixe(1,ngr) &
                           ,iyb(1,ngr),iye(1,ngr))

    if (nmachs .gt. 1) then
      print*,'!------------------------------------------------'
      print*,'!            Domain decomposition'
      print*,'!   grid# node# x-beg x-end y-beg y-end #cols'
      print*,'!------------------------------------------------'
      do jnode = 1,nmachs
         ncols= (1+ixe(jnode,ngr)-ixb(jnode,ngr))  &
               *(1+iye(jnode,ngr)-iyb(jnode,ngr))
         print '('' ! '',7i6,f12.4)', ngr,jnode,ixb(jnode,ngr) &
               ,ixe(jnode,ngr),iyb(jnode,ngr),iye(jnode,ngr),ncols
      enddo
      print*,'!------------------------------------------------'
      print*,''
    endif

    !Determine if this node has the SCM point to print in "subroutine timestep".
    !Only do this if RAMS running a single grid; no nested grids for SCM output.
    my_scm_num = -999
    if (iscm==ngr .and. nmachs==1) my_scm_num = 1
    if (iscm==ngr .and. nmachs>1) then
      do jnode = 1,nmachs
         if( ixb(jnode,ngr) <= iscmx .and. iscmx <= ixe(jnode,ngr) .and. &
             iyb(jnode,ngr) <= iscmy .and. iscmy <= iye(jnode,ngr) ) then
           my_scm_num = jnode
           print*,'SCM POINT ON NODE: ',my_scm_num
         endif
      enddo
    endif

    ! Figure out and record which sub domains are in the corners
    ! of the full domain
    do jnode = 1, nmachs
      if ((ixb(jnode,ngr) .eq.             2) .and. (iyb(jnode,ngr) &
         .eq.             2)) fd_sw_num(ngr) = jnode
      if ((ixb(jnode,ngr) .eq.             2) .and. (iye(jnode,ngr) &
         .eq. (nnyp(ngr)-1))) fd_nw_num(ngr) = jnode
      if ((ixe(jnode,ngr) .eq. (nnxp(ngr)-1)) .and. (iyb(jnode,ngr) &
         .eq.             2)) fd_se_num(ngr) = jnode
      if ((ixe(jnode,ngr) .eq. (nnxp(ngr)-1)) .and. (iye(jnode,ngr) &
         .eq. (nnyp(ngr)-1))) fd_ne_num(ngr) = jnode
    enddo
  enddo

return
END SUBROUTINE assign_node_subdomain

!##############################################################################
Subroutine node_init ()

  use mem_grid
  use node_mod

  implicit none

  integer :: ngr, n

!!DEBUG
!  integer :: inode, itype, info
!!DEBUG

  ! finish off the domain decomposition by filling in boundary information
  CALL par_decomp_bounds ()

  ! figure out which nodes need to talk to each other and record these paths
  CALL par_node_paths (maxgrds,ngrids,nxpmax,nypmax,nnxp,nnyp  &
   ,nxtnest,my_rams_num,maxmach,nmachs,ipm,jpm,ixb,ixe,iyb,iye,ipaths &
   ,iget_paths,ibnd,jbnd,not_a_node)

!!DEBUG
!do ngr = 1,ngrids
!   do inode=1,nmachs
!      print'(a,i0,a)', 'DEBUG: NODE', my_rams_num, ':'
!      print'(a,i0,a)', 'DEBUG: NODE', my_rams_num, ':  ngr idn itype info ipaths '
!      do itype=1,7
!        print'(a,i0,a)', 'DEBUG: NODE', my_rams_num, ':'
!         do info=1,5
!            print'(a,i0,a,5i6)', 'DEBUG: NODE', my_rams_num, ': ',ngr,inode,itype,info &
!               ,ipaths(info,itype,ngr,inode)
!         enddo
!      enddo
!   enddo
!enddo
!
!do ngr = 1,ngrids
!   do inode = 1,nmachs
!      print'(a,i0,a)', 'DEBUG: NODE', my_rams_num, ':'
!      print'(a,i0,a)', 'DEBUG: NODE', my_rams_num, ':  ngr isn itype iget_paths '
!      do itype = 1,6
!         print'(a,i0,a,4i6)', 'DEBUG: NODE', my_rams_num, ': ',ngr,inode,itype &
!             ,iget_paths(itype,ngr,inode)
!      enddo
!   enddo
!enddo
!!DEBUG

  ! In the parallel nesting scheme, each subgrid on the fine
  ! mesh will keep a copy of the coarse mesh "tile", for each
  ! variable, that covers the fine mesh subgrid enough for doing
  ! interpolation. 
  !
  do ngr = 2,ngrids
    ! Determine size of coarse tile
    !
    ! Need extra grids cells beyond those recorded in [ij]pm arrays for
    ! the interpolation.
    !      2 extra parent grid cells at the "left" end
    !      1 extra parent grid cell at the "right" end
    ctxbeg(ngr) = ipm(mxbeg(ngr),ngr) - 2
    ctxend(ngr) = ipm(mxend(ngr),ngr) + 1 
    if(jdim == 1) then
      ctybeg(ngr) = jpm(mybeg(ngr),ngr) - 2
      ctyend(ngr) = jpm(myend(ngr),ngr) + 1 
    else
      ctybeg(ngr) = 1
      ctyend(ngr) = 1
    endif

    ctnxp(ngr) = (ctxend(ngr) - ctxbeg(ngr)) + 1
    ctnyp(ngr) = (ctyend(ngr) - ctybeg(ngr)) + 1
    ctnzp(ngr) = mmzp(ngr)

    ! Determine offset of coarse tile
    ! This is the distance into the coarse grid. The array holding the coarse tile will go from 1..nxc and 
    ! 1..nyc, so the offset is the distance from (1,1) in the coarse grid to where the coarse tile
    ! is located.
    cti0(ngr) = ctxbeg(ngr) - 1
    ctj0(ngr) = ctybeg(ngr) - 1

    ! Shift ipm and jpm to hold indices relative to the origin being at (1,1)
    ! Copy kpm, for now.
    do n = 1,mmxp(ngr)
      ctipm(n,ngr) = ipm(n+mi0(ngr),ngr) - cti0(ngr)
    enddo
    do n = 1,mmyp(ngr)
      ctjpm(n,ngr) = jpm(n+mj0(ngr),ngr) - ctj0(ngr)
    enddo
    do n = 1,mmzp(ngr)
      ctkpm(n,ngr) = kpm(n,ngr)
    enddo
  enddo

return
END SUBROUTINE node_init

!##############################################################################
Subroutine par_set_comm_buff_sizes ()

use mem_grid
use node_mod
use var_tables
use micro_prm, only: nkr
use kpp_parameters, only: nkppz

implicit none

  integer :: inode,ng,mp_nzp,icm,ifm,nestvar,nf,nc,nv
  integer :: isn,idn,itype
  integer :: i1,i2,j1,j2
  integer :: ixy,ixyz,num_lbc_buff,num_nest_buff,num_feed_buff
  integer :: npvar2,npvar3,npvar7
  integer :: nivar2,nivar3,nivars4,nivarw4,nivars3,nivarb4,nivark4
  integer :: num_init_buff

  integer, dimension(nmachs) :: recv_buff_sizes

  ! Compute  send and receive buffer sizes. These will be maximum of
  ! long timestep, initialization, turbulence, nest boundaries, and nest feedback.
  !
  ! Small timestep will use same buffers as they are always smaller.
  !
  ! node_buffs(*)%nsend - sizes when I am the source node
  ! node_buffs(*)%nrecv - sizes when I am the destination node

  do inode=1,nmachs
     node_buffs(inode)%nsend=0
     node_buffs(inode)%nrecv=0
  enddo

  mp_nzp=0
  do ng=1,ngrids
     mp_nzp=max(mp_nzp,mmzp(ng))
  enddo

  
  do ng=1,ngrids
  
  !  Find number of nested variables to be communicated.
     icm=ng
     ifm=ng
     if(ng /= 1) icm=nxtnest(ifm)
     nestvar = NUM_NEST_VARS
     do nf=1,num_scalar(ifm)
        do nc=1,num_scalar(icm)
           if(scalar_tab(nf,ifm)%name==scalar_tab(nc,icm)%name)  &
                    nestvar=nestvar+1
        enddo
     enddo
  
  !  Find number of lbc variables to be communicated.
     npvar3=0 ; npvar2=0 ; npvar7=0
     do nv = 1,num_var(ng)
        if(vtab_r(nv,ng)%impt1 == 1 ) then
           if (vtab_r(nv,ng)%idim_type==2) npvar2=npvar2+1
           if (vtab_r(nv,ng)%idim_type==3) npvar3=npvar3+1
           if (vtab_r(nv,ng)%idim_type==7) npvar7=npvar7+1
        endif
     enddo
  
  !  Find number of init variables to be communicated.
     nivar2 = 0
     nivar3 = 0
     nivars4 = 0
     nivarw4 = 0
     nivars3 = 0
     nivarb4 = 0
     nivark4 = 0
     do nv = 1,num_var(ng)
        if(vtab_r(nv,ng)%impti == 1 ) then
           if (vtab_r(nv,ng)%idim_type==2) nivar2=nivar2+1
           if (vtab_r(nv,ng)%idim_type==3) nivar3=nivar3+1
           if (vtab_r(nv,ng)%idim_type==4) nivars4=nivars4+1
           if (vtab_r(nv,ng)%idim_type==5) nivarw4=nivarw4+1
           if (vtab_r(nv,ng)%idim_type==6) nivars3=nivars3+1
           if (vtab_r(nv,ng)%idim_type==7) nivarb4=nivarb4+1
           if (vtab_r(nv,ng)%idim_type==10)nivark4=nivark4+1
        endif
     enddo

     ! Figure out buffer sizes when I am the sending node 
     do idn=1,nmachs
        num_lbc_buff=0
        num_init_buff=0
        num_nest_buff=0
        num_feed_buff=0
  
        itype=1
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)
        if(i1.ne.0) then
           ixy=(i2-i1+1)*(j2-j1+1)
           ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)

           !Buffer for LBC things
           num_lbc_buff=ixyz*npvar3 &
                        +ixy*npvar2 &
                       +ixyz*nkr*npvar7  &
                +2*(npvar3+npvar2+npvar7+100)

           !Buffer for all things initialized (should be all vars)
           num_init_buff =                &
               (ixy*nivar2)               & !Standard 2D vars
             + (ixyz*nivar3)              & !Standard 3D vars
             + (nivars4*ixy*nzg*npatch)   & !Sfc model soil-levels/patches
             + (nivarw4*ixy*nzs*npatch)   & !Sfc model snow-levels/patches
             + (nivars3*ixy*npatch)       & !Sfc model patches
             + (nivarb4*ixyz*nkr)         & !HUCM-BIN micro bins
             + (nivark4*ixyz*nkppz)       & !KPP ocean levels
             + 2*(nivar2+nivar3+nivars4+nivarw4 &
                 +nivars3+nivarb4+nivark4+100)
        endif
  
        itype=5
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)
        if(i1.ne.0) then
           ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
           num_nest_buff=ixyz*nestvar+2*(nestvar+100)
        endif
  
        itype=6
        i1=ipaths(1,itype,ng,idn)
        i2=ipaths(2,itype,ng,idn)
        j1=ipaths(3,itype,ng,idn)
        j2=ipaths(4,itype,ng,idn)
        if(i1.ne.0) then
           ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
           num_feed_buff=ixyz*nestvar+2*(nestvar+100)
        endif

        node_buffs(idn)%nsend = max(node_buffs(idn)%nsend  &
             ,num_lbc_buff,num_init_buff,num_nest_buff,num_feed_buff)
     enddo
  enddo

  ! At this point we have placed the max buffer sizes for sending data to 
  ! other nodes considering all grids. Ie, one send size will account for 
  ! what is needed for all grids.

  ! Now figure out buffer sizes when this node is the destination node.
  !
  ! We need information from the other nodes in order to accomplish
  ! this since our copy of ipaths only accounts for when we are the
  ! sending node (and when another node is sending to us, the buffer
  ! size could be different). Can use an MPI gather to collect the
  ! send buffer sizes from the other nodes and load them (in correct
  ! order) into our receive buffer size vars.
  !
  ! The basic idea is that for node 1, the receive buffer sizes have been 
  ! computed in the first entry of all the other nodes' send buffer size 
  ! arrays (since they are sending to node 1). For node 2, we want to 
  ! collect the second entry of all the other nodes' send buffer size 
  ! arrays. For node 3, collect the third entries, and so on.
  !
  ! The algorithm then becomes:
  !   For each node number
  !     Issue a gather sending your individual send buffer size and 
  !     collecting into a buffer large enough to hold all the numbers 
  !     from the other nodes
  !   Once the gather loop is finished, dole out the numbers in your 
  !   buffer into your receive buffer size entries.
  !
  ! Every node must issue a gather to all nodes including itself. This 
  ! makes the gathers match up across the collection of processes and 
  ! each node receives the entries it needs.
  !
  ! Note that this algorithm effectively transposes the send matrix 
  ! (formed by the send arrays from all the nodes) into the receive matrix.
    
  do isn = 1, nmachs
    CALL par_gather_ints (node_buffs(isn)%nsend, recv_buff_sizes, 1, machnum(isn))
  enddo

  do isn = 1, nmachs
    node_buffs(isn)%nrecv = recv_buff_sizes(isn)
  enddo

return
END SUBROUTINE par_set_comm_buff_sizes

!##############################################################################
Subroutine par_decomp_bounds ()

use mem_grid
use node_mod
use micro_prm, only: nkr
use kpp_parameters, only: nkppz

implicit none

  integer :: igrid, inode
  integer :: nxpts, nypts

  !Compute various subdomain boundary numbers for this node
  ! mxbeg,mybeg,mxend,myend - portions of full domain that node will have
  !                         - includes overlap region
  ! mi0,mj0  - subdomain offsets relative to full domain
  ! mia,miz,mja,mjz - subdomain "compute" points,
  !              or normal thermodynamic tendency points (2-nx, 2-ny for
  !              non-parallel run
  ! mibcon - flag denoting if real boundary is on subdomain
  !     bit 1=west, bit 2=east, bit 3=south, bit 4=north
  !
  ! Subdomain boundary numbers for parallel HDF5 file IO
  ! mem_read(igrid)%xblock, mem_read(igrid)%yblock     
  !     - memory block size for read operation
  ! mem_write(igrid)%xblock, mem_write(igrid)%yblock   
  !     - memory block size for write operation
  ! mem_read(igrid)%xoff, mem_read(igrid)%yoff         
  !     - memory subdomain offset relative to subdomain for read operation
  ! mem_write(igrid)%xoff, mem_write(igrid)%yoff       
  !     - memory subdomain offset relative to subdomain for write operation
  ! file_read(igrid)%xblock, file_read(igrid)%yblock   
  !     - file block size for read operation
  ! file_write(igrid)%xblock, file_write(igrid)%yblock 
  !     - file block size for write operation
  ! file_read(igrid)%xoff, file_read(igrid)%yoff       
  !     - file subdomain offset relative to full domain for read operation
  ! file_write(igrid)%xoff, file_write(igrid)%yoff     
  !     - file subdomain offset relative to full domain for write operation
  ! file_xchunk, file_ychunk             
  !     - chunk size write operation

  ! This routine assumes that boundaries are all 1 point wide

  ! need chunk size to be max(all mmxp) by max(all mmyp)
  do igrid = 1, ngrids
    mibcon(igrid) = 0

    file_xchunk(igrid) = 0
    file_ychunk(igrid) = 0

    do inode = 1, nmachs
      ! Find the size of inode's sub-domain
      nxpts = (ixe(inode,igrid) - ixb(inode,igrid)) + 3
      nypts = (iye(inode,igrid) - iyb(inode,igrid)) + 3

      ! Save these if they are the max so far
      if (nxpts > file_xchunk(igrid)) then
        file_xchunk(igrid) = nxpts
      endif
      if (nypts > file_ychunk(igrid)) then
        file_ychunk(igrid) = nypts
      endif

      ! If we are on the inode matched with this process, then record
      ! and calculate all of the other boundary information
      !
      ! Calculate the numbers as if the sub-domain doesn't share a full
      ! domain boundary. Then check for the full domain boundary and adjust
      ! numbers accordingly.
      !
      ! ixb,ixe,iyb,iye denote the compute region of the subdomain
      if (inode .eq. my_rams_num) then
        mmxp(igrid) = nxpts
        mmyp(igrid) = nypts
        mmzp(igrid) = nnzp(igrid)

        !print*,'Z-pts,node,nnzp,mmzp: ',inode,nnzp(igrid),mmzp(igrid)
        !print*,'X-pts,node,nnxp,mmxp: ',inode,nnxp(igrid),mmxp(igrid)
        !print*,'Y-pts,node,nnyp,mmyp: ',inode,nnyp(igrid),mmyp(igrid)
        !print*,''

        mmxyzp(igrid) = mmxp(igrid) * mmyp(igrid) * mmzp(igrid)
        mmxysp(igrid) = mmxp(igrid) * mmyp(igrid) * (nzg+nzs+3) * npatch
        mmxyp(igrid)  = mmxp(igrid) * mmyp(igrid)
        mmxyzbp(igrid) = mmxyzp(igrid) * nkr
        mmxyzkp(igrid) = mmxp(igrid) * mmyp(igrid) * nkppz

        ! mxbeg, mxend, mybeg, myend extend 1 beyond ixb,ixe,iyb,iye
        mxbeg(igrid) = ixb(inode,igrid) - 1
        mxend(igrid) = ixe(inode,igrid) + 1
        mybeg(igrid) = iyb(inode,igrid) - 1
        myend(igrid) = iye(inode,igrid) + 1

        ! offsets are 2 less than ixb,iyb
        mi0(igrid) = ixb(inode,igrid) - 2
        mj0(igrid) = iyb(inode,igrid) - 2

        ! mia,mja are 2
        ! miz,mjz are 2 + (i[xy]e-i[xy]b)
        mia(igrid) = 2
        miz(igrid) = 2 + (ixe(inode,igrid) - ixb(inode,igrid))
        mja(igrid) = 2
        mjz(igrid) = 2 + (iye(inode,igrid) - iyb(inode,igrid))

        ! generate io sizes and offsets
        ! read operation uses entire sub-domain (compute area plus boundaries)
        mem_read(igrid)%xblock = mmxp(igrid)
        mem_read(igrid)%yblock = mmyp(igrid)
        mem_read(igrid)%xoff   = 0
        mem_read(igrid)%yoff   = 0

        file_read(igrid)%xblock = mem_read(igrid)%xblock
        file_read(igrid)%yblock = mem_read(igrid)%yblock
        file_read(igrid)%xoff   = mi0(igrid)
        file_read(igrid)%yoff   = mj0(igrid)
 
        ! write operations uses sub-domain area minus interior boundaries
        ! set values here as if sub-domain is completely surrounded 
        ! by interior boundaries
        mem_write(igrid)%xblock = (miz(igrid) - mia(igrid)) + 1
        mem_write(igrid)%yblock = (mjz(igrid) - mja(igrid)) + 1
        mem_write(igrid)%xoff   = 1
        mem_write(igrid)%yoff   = 1

        file_write(igrid)%xblock = mem_write(igrid)%xblock
        file_write(igrid)%yblock = mem_write(igrid)%yblock
        file_write(igrid)%xoff   = mi0(igrid) + 1
        file_write(igrid)%yoff   = mj0(igrid) + 1

        ! check west boundary
        if (ixb(inode,igrid) .eq. 2) then
          mibcon(igrid) = mibcon(igrid) + 1
          ! add west boundary back into the io write region
          mem_write(igrid)%xblock = mem_write(igrid)%xblock + 1
          mem_write(igrid)%xoff = mem_write(igrid)%xoff - 1

          file_write(igrid)%xblock = file_write(igrid)%xblock + 1
          file_write(igrid)%xoff = file_write(igrid)%xoff - 1
        endif
        
        ! check east boundary
        if (ixe(inode,igrid) .eq. (nnxp(igrid)-1)) then
          mibcon(igrid) = mibcon(igrid) + 2
          ! add east boundary back into the io write region
          mem_write(igrid)%xblock = mem_write(igrid)%xblock + 1

          file_write(igrid)%xblock = file_write(igrid)%xblock + 1
        endif
        
        ! check south boundary
        if (iyb(inode,igrid) .eq. 2) then
          mibcon(igrid) = mibcon(igrid) + 4
          ! add south boundary back into the io write region
          mem_write(igrid)%yblock = mem_write(igrid)%yblock + 1
          mem_write(igrid)%yoff = mem_write(igrid)%yoff - 1

          file_write(igrid)%yblock = file_write(igrid)%yblock + 1
          file_write(igrid)%yoff = file_write(igrid)%yoff - 1
        endif
        
        ! check north boundary
        if (iye(inode,igrid) .eq. (nnyp(igrid)-1)) then
          mibcon(igrid) = mibcon(igrid) + 8
          ! add north boundary back into the io write region
          mem_write(igrid)%yblock = mem_write(igrid)%yblock + 1

          file_write(igrid)%yblock = file_write(igrid)%yblock + 1
        endif
        
      endif
    enddo
  enddo

  ! In the special case where we only have one machine (nmachs ==1), 
  ! We need to make sure that the chunk size is not too big. 
  ! If the chunk size ends up too big, the entire domain will need to 
  ! be read into RAM. That won't work for sufficiently large domain sizes.
  ! So, in the special case where number of machs is 1, set chunk size
  ! to be equal to the minimum between number of points and an arbitrary value.
  ! In this case, an arbitrary value of 50 is chosen to match the value set 
  ! in core/rams_model.f90. See more info there.
  if (nmachs ==1) then
    do igrid = 1, ngrids
      file_xchunk(igrid) = min(50, nnxp(igrid))
      file_ychunk(igrid) = min(50, nnyp(igrid))
    enddo
  endif

return
END SUBROUTINE par_decomp_bounds
