!##############################################################################
Module mem_scratch

use grid_dims

implicit none

   Type scratch_vars

   real, allocatable, dimension(:) :: &
                         scr1,scr2 &
                        ,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dg  &
                        ,vt3dh,vt3di,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn  &
                        ,vt3do,vt3dp,u_the_d,v_the_d,w_the_d &
                        ,p_the_d
   real, allocatable, dimension(:) :: &
                        vt2da,vt2db,vt2dc,vt2dd,vt2de,vt2df
                         
   End Type
   
   type (scratch_vars) :: scratch

   real, dimension(maxdim) ::  vctr1 ,vctr2 ,vctr3 ,vctr4 ,vctr5 ,vctr6   &
                              ,vctr7 ,vctr8 ,vctr9 ,vctr10,vctr11,vctr12  &
                              ,vctr13,vctr14,vctr15,vctr16,vctr17,vctr18  &
                              ,vctr19,vctr20,vctr21,vctr22,vctr23,vctr24  &
                              ,vctr25,vctr26,vctr27,vctr28,vctr29,vctr30  &
                              ,vctr31,vctr32,vctr33,vctr34,vctr35,vctr36  &
                              ,vctr37,vctr38,vctr39,vctr40,vctr41

Contains

!##############################################################################
Subroutine alloc_scratch (numz_sd,numx_sd,numy_sd,numz_fd,numx_fd,numy_fd  &
              ,ngrs,nzg,nzs,npatch,maxx,maxy,maxz)

use micro_prm, only: nkr
use micphys, only: level

implicit none
   
   integer, dimension (*) :: numz_sd,numx_sd,numy_sd,numz_fd,numx_fd,numy_fd
   ! num[xyz]_sd - number of x,y,z points in sub domain (compute node)
   ! num[xyz]_fd - number of x,y,z points in full domain
   integer :: ngrs,nzg,nzs,npatch
   integer :: ng,ntpts3,ntpts2,ntptsx,maxx,maxy,maxz

!print*, 'enter alloc_scratch'

!         Find the maximum number of grid points needed for any grid.
!           The max points in each direction are passed back for use by 
!           various nesting things.    

   maxx = 0
   maxy = 0
   maxz = 0
   ntpts3=0
   ntpts2=0
   do ng=1,ngrs
      maxx = max(maxx,numx_fd(ng))
      maxy = max(maxy,numy_fd(ng))
      maxz = max(maxz,numz_fd(ng))
      ntpts3=max( numx_sd(ng)*numy_sd(ng)*numz_sd(ng),ntpts3 )
      ntpts2=max( numx_sd(ng)*numy_sd(ng),ntpts2 )
   enddo
   ! scr1 and scr2 needs to be the max of a passed field
   if (level == 4) then
      ntptsx=max(maxx*maxy*maxz,ntpts2*nzg*npatch, &     
                                 ntpts2*nzs*npatch,ntpts3*nkr)+1000
   else
      ntptsx=max(maxx*maxy*maxz,ntpts2*nzg*npatch,ntpts2*nzs*npatch)+1000
   endif

! Allocate arrays based on options (if necessary).
!   -scr1 and scr2 need to be allocated to full domain (even on compute nodes)
!        to max(nx)*max(ny)*max(nz)
!   -do not need all these arrays if it is a master process in a parallel run,
!      so just allocate some to 1 word.

      allocate (scratch%scr1 (ntptsx))
      allocate (scratch%scr2 (ntptsx))
      allocate (scratch%vt3da(ntpts3))
      allocate (scratch%vt3db(ntpts3))     

      allocate (scratch%vt3dc(ntpts3))
      allocate (scratch%vt3dd(ntpts3))
      allocate (scratch%vt3de(ntpts3))
      allocate (scratch%vt3df(ntpts3))
      allocate (scratch%vt3dg(ntpts3))
      allocate (scratch%vt3dh(ntpts3))
      allocate (scratch%vt3di(ntpts3))
      allocate (scratch%vt3dj(ntpts3))
      allocate (scratch%vt3dk(ntpts3))
      allocate (scratch%vt3dl(ntpts3))
      allocate (scratch%vt3dm(ntpts3))
      allocate (scratch%vt3dn(ntpts3))
      allocate (scratch%vt3do(ntpts3))
      allocate (scratch%u_the_d(ntpts3))
      allocate (scratch%v_the_d(ntpts3))
      allocate (scratch%w_the_d(ntpts3))
      allocate (scratch%p_the_d(ntpts3))
      if (level<4) then
         allocate (scratch%vt3dp(ntpts3))
      else
         allocate (scratch%vt3dp(ntpts3*nkr))
      endif

      allocate (scratch%vt2da(ntpts2))
      allocate (scratch%vt2db(ntpts2))     
      allocate (scratch%vt2dc(ntpts2))
      allocate (scratch%vt2dd(ntpts2))
      allocate (scratch%vt2de(ntpts2))
      allocate (scratch%vt2df(ntpts2))

return
END SUBROUTINE alloc_scratch

!##############################################################################
Subroutine dealloc_scratch ()
   
implicit none

! Deallocate all scratch arrays

   if (allocated(scratch%scr1 ))  deallocate (scratch%scr1 )
   if (allocated(scratch%scr2 ))  deallocate (scratch%scr2 )
   if (allocated(scratch%vt3da))  deallocate (scratch%vt3da)
   if (allocated(scratch%vt3db))  deallocate (scratch%vt3db)
   if (allocated(scratch%vt3dc))  deallocate (scratch%vt3dc)
   if (allocated(scratch%vt3dd))  deallocate (scratch%vt3dd)
   if (allocated(scratch%vt3de))  deallocate (scratch%vt3de)
   if (allocated(scratch%vt3df))  deallocate (scratch%vt3df)
   if (allocated(scratch%vt3dg))  deallocate (scratch%vt3dg)
   if (allocated(scratch%vt3dh))  deallocate (scratch%vt3dh)
   if (allocated(scratch%vt3di))  deallocate (scratch%vt3di)
   if (allocated(scratch%vt3dj))  deallocate (scratch%vt3dj)
   if (allocated(scratch%vt3dk))  deallocate (scratch%vt3dk)
   if (allocated(scratch%vt3dl))  deallocate (scratch%vt3dl)
   if (allocated(scratch%vt3dm))  deallocate (scratch%vt3dm)
   if (allocated(scratch%vt3dn))  deallocate (scratch%vt3dn)
   if (allocated(scratch%vt3do))  deallocate (scratch%vt3do)
   if (allocated(scratch%vt3dp))  deallocate (scratch%vt3dp)
   if (allocated(scratch%u_the_d))  deallocate (scratch%u_the_d)
   if (allocated(scratch%v_the_d))  deallocate (scratch%v_the_d)
   if (allocated(scratch%w_the_d))  deallocate (scratch%w_the_d)
   if (allocated(scratch%p_the_d))  deallocate (scratch%p_the_d)
   if (allocated(scratch%vt2da))  deallocate (scratch%vt2da)
   if (allocated(scratch%vt2db))  deallocate (scratch%vt2db)
   if (allocated(scratch%vt2dc))  deallocate (scratch%vt2dc)
   if (allocated(scratch%vt2dd))  deallocate (scratch%vt2dd)
   if (allocated(scratch%vt2de))  deallocate (scratch%vt2de)
   if (allocated(scratch%vt2df))  deallocate (scratch%vt2df)
        
return
END SUBROUTINE dealloc_scratch

END MODULE mem_scratch
