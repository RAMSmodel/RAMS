!##############################################################################
Subroutine modsched (ischtab,ntab,ngrids,nxtnest,nndtrat,nsub)

use mem_grid, only:print_msg

implicit none

integer :: ntab,ngrids,nsub,ischtab(ntab,*),nxtnest(*),nndtrat(*)  &
          ,ncnt(50),nnfm(50)
integer :: ng,nt,is,itopgrd,ntt,ngrid,isstp,n,k

do ng=1,ngrids
   nnfm(ng)=1
   ncnt(ng)=0
enddo
do nt=1,4
   do ntt=1,ntab
      ischtab(ntt,nt)=0
   enddo
enddo

is=0
itopgrd = 0
do ng=1,ngrids
   if(nxtnest(ng).eq.0)then
      ngrid=ng

30         ISSTP=MOD(NCNT(NGRID),NNDTRAT(NGRID))+1
      is=is+1

!        Table value 1 is the grid number to be stepped forward
      ischtab(is,1)=ngrid
!        Table value 3 is the sub-timestep counter
      ischtab(is,3)=ISSTP

      NCNT(NGRID)=NCNT(NGRID)+1

!        Table value 5 is the total sub-timestep counter
      ischtab(is,5)=NCNT(NGRID)

40         NNFM(NGRID)=NNFM(NGRID)+1
      IF (NNFM(NGRID).GT.NGRIDS) GO TO 50
      IF (NXTNEST(NNFM(NGRID)).NE.NGRID) GO TO 40
      NGRID=NNFM(NGRID)

!        Table value 2 is the grid number to interpolate. 0 = no interpolate.
      ischtab(is,2)=ngrid

      ISSTP=MOD(NCNT(NGRID),NNDTRAT(NGRID))+1

      GO TO 30

50         NNFM(NGRID)=1
      IF (NXTNEST(NGRID).EQ.0) GO TO 80
      IF(NCNT(NGRID).LT.NCNT(NXTNEST(NGRID))*NNDTRAT(NGRID))THEN
         GO TO 30
      ENDIF

!        Table value 4 is the number of grids to be fedback
      ischtab(is,4)=ischtab(is,4)+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ngrid=nxtnest(ngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      GO TO 40

80         CONTINUE

   endif
enddo

nsub=is

if(print_msg) then
  print*,'=== Timestep Schedule ===='
  print 333,(n,(ischtab(n,k),k=1,5),n=1,is)
  print*,''
333 format (6i5)
endif

return
END SUBROUTINE modsched

!##############################################################################
Subroutine cfl (n1,n2,n3,i0,j0,mynum)

use mem_basic
use mem_grid

implicit none

integer :: n1,n2,n3,i0,j0,mynum

CALL cfll (n1,n2,n3,i0,j0,mynum  &
   ,basic_g(ngrid)%up   (1,1,1)  ,basic_g(ngrid)%vp   (1,1,1)  &
   ,basic_g(ngrid)%wp   (1,1,1)  ,grid_g(ngrid)%rtgt    (1,1)  &
   ,grid_g(ngrid)%f13t    (1,1)  ,grid_g(ngrid)%f23t    (1,1)  &
   ,grid_g(ngrid)%dxt     (1,1)  ,grid_g(ngrid)%dyt     (1,1)  )

return
END SUBROUTINE cfl

!##############################################################################
Subroutine cfll (n1,n2,n3,i0,j0,mynum,up,vp,wp,rtgt,f13t,f23t,dxt,dyt)

use mem_grid
use mem_scratch

implicit none

integer :: n1,n2,n3,i0,j0,mynum
real, dimension(n1,n2,n3) :: up,vp,wp
real, dimension(n2,n3)    :: rtgt,f13t,f23t,dxt,dyt
         
integer :: i,j,k,ifm,icm,innest,nprints
real :: c1x,c1y,c1z,cflnumh,cflnumv

!     This routine returns flags the model to bring itself down when the CFL
!     linear stability criteria on advection is exceeded.
!     (Actually check on 90% of CFL)

5    format('NODE',i0,': cflx,ngrid,k,i,j = ',f5.1,4i5)
6    format('NODE',i0,': cfly,ngrid,k,i,j = ',f5.1,4i5)
7    format('NODE',i0,': cflz,ngrid,k,i,j = ',f5.1,4i5)

nprints = 0
cflnumh = .90
cflnumv = .90
cflxy(ngrid) = 0.
cflz(ngrid) = 0.
vertvel(ngrid) = 0.

! Let's try a new thing... if we have a grid point that is on a 
!   coarse grid, but it is under a nested grid, we will ignore it
!   under the assumption that the fine grid values will overwrite
!   it, hence not allowing it to go numerically unstable.

jloop: do j = 1+jdim,n3-jdim
   iloop: do i = 2,n2-1
   
      ! See if this is under a fine grid horizontally... ignore vertical for now
      innest=0
      if (ngrids > ngrid) then
         do ifm=ngrid+1,ngrids
            icm=nxtnest(ifm)
            if(icm == ngrid .and. &
               i+i0 >= ipm(1,ifm) .and. i+i0 <= ipm(nnxp(ifm),ifm) .and. &
               j+j0 >= jpm(1,ifm) .and. j+j0 <= jpm(nnyp(ifm),ifm) ) then
               innest=1
               exit
            endif
         enddo
      endif
      
      if(innest == 1) then
         !print '(a,5i4)', 'cfl under grid-' &
         !         ,ngrid,ifm,i+i0,j+j0,mynum
         cycle iloop
      endif
   
      kloop: do k = 2,n1-1
      
         vctr1(k) = .5*(up(k,i,j)+up(k,i-1,j))*dtlt*dxt(i,j)
         vctr2(k) = .5*(vp(k,i,j)+vp(k,i,j-jdim))*dtlt*dyt(i,j)
         vctr3(k) = ((wp(k,i,j)+wp(k-1,i,j))  &
           +(up(k,i,j)+up(k,i-1,j))*f13t(i,j)*ht(k)*rtgt(i,j)  &
           +(vp(k,i,j)+vp(k,i,j-jdim))*f23t(i,j)*ht(k)*rtgt(i,j)  &
           )*.5*dtlt*dzt(k)
      enddo kloop
      
      do k = 2,n1-1
         c1x = abs(vctr1(k))
         c1y = abs(vctr2(k))
         c1z = abs(vctr3(k))

         if (nprints .le. 10) then
            if (c1x .gt. cflnumh) then
               nprints = nprints + 1
               print 5, mynum,c1x,ngrid,k,i+i0,j+j0
               print'(a,i0,a,5e15.4)','NODE',mynum,': ' &
                 ,up(k,i,j),up(k,i-1,j),dtlt,dxt(i,j),vctr1(k)
            endif
            if (c1y .gt. cflnumh) then
               nprints = nprints + 1
               print 6, mynum,c1y,ngrid,k,i+i0,j+j0
            endif
            if (c1z .gt. cflnumv) then
               nprints = nprints + 1
               print 7, mynum,c1z,ngrid,k,i+i0,j+j0
            endif
         endif

         if (c1x .gt. cflxy(ngrid)) cflxy(ngrid) = c1x
         if (c1y .gt. cflxy(ngrid)) cflxy(ngrid) = c1y
         if (c1z .gt. cflz(ngrid)) cflz(ngrid) = c1z
         if (wp(k,i,j) .gt. vertvel(ngrid)) vertvel(ngrid) = wp(k,i,j)
      enddo
   enddo iloop
enddo jloop

return
END SUBROUTINE cfll

!##############################################################################
Subroutine dtset ()

use mem_grid
use rconstants
use ref_sounding
use io_params
use node_mod

implicit none

real, dimension(maxgrds) :: sscourn
real, dimension(nzpmax) :: vctr1

integer, save :: icflfirst
real, save :: cflnumh,cflnumv,ssodx(maxgrds)

integer :: ifm,id,n2,n3,k,icm
real :: ssmax,tmax,dxtmax,cflxyz,timeleft,dxt_sw,dxt_nw,dxt_se,dxt_ne

! On the first call to this routine, initialize ssodx, dtlongn,
! nnacoust, and if required, nndtrat.
data icflfirst/0/

if (icflfirst .eq. 0) then

   icflfirst=1

   cflnumh = .90
   cflnumv = .90

   do ifm = 1,ngrids
      n2 = mmxp(ifm)
      n3 = mmyp(ifm)
      do k = 1,mmzp(ifm)
         vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
      enddo
      tmax = maxval(vctr1(1:mmzp(ifm)))
      ssmax = sqrt(cp / cv * rgas * tmax)

      ! First record dxt at the sub domain corners. Note that if
      ! there is only one node, then mmxp and mmyp will describe
      ! the full domain corners. get_fd_corner_vals will read
      ! in the sub domain corner values and return the full
      ! domain corner values.
      dxt_sw = grid_g(ifm)%dxt(1,1)
      dxt_nw = grid_g(ifm)%dxt(1,mmyp(ifm))
      dxt_se = grid_g(ifm)%dxt(mmxp(ifm),1)
      dxt_ne = grid_g(ifm)%dxt(mmxp(ifm),mmyp(ifm))
      CALL get_fd_corner_vals (ifm, dxt_sw, dxt_nw, dxt_se, dxt_ne)

      dxtmax = max(dxt_sw, dxt_nw, dxt_se, dxt_ne)
      ssodx(ifm) = ssmax * dxtmax
   enddo

   sspct = 1.
   do ifm = 1,ngrids
         icm = nxtnest(ifm)
         if (icm .eq. 0) then
            dtlongn(ifm) = dtlong
         else
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif
         nnacoust(ifm) = nacoust
         sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
         sspct = min(sspct,  &
            .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))            
   enddo

! Print out initial values of dtlongn, nndtrat, nnacoust, sscourn, and sspct
   if(print_msg) then
     write(6,120)
120        format(/,'Initial timestep info: ngrid, nndtrat, nnacoust,'  &
            ,' dtlongn, sscourn, sspct')

     do ifm = 1,ngrids

        write(6,121) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)  &
           ,sscourn(ifm),sspct
121        format(23x,i3,i8,i9,f13.3,2f8.3)
     enddo
     print*,''
   endif

endif

! check Courant numbers

do ifm = 1,ngrids
      cflxyz = max(cflxy_max(ifm)/cflnumh,cflz_max(ifm)/cflnumv)
      if (cflxyz > 1.) then
         iflag = 1
         print'(a,i0,a)', 'NODE', my_rams_num, &
          ': Model will stop because CFL limit exceeded.'
      endif
enddo

! For the coarse grid(s), adjust dtlongn(1) to not jump past an analysis
! write time or the end of the model run time.

sspct = 1.
do ifm = 1,ngrids
         icm = nxtnest(ifm)
         if (ifm == 1) then
            dtlongn(ifm) = dtlong
            timeleft =  &
               min (timmax - time, frqstate(1) - mod(time,frqstate(1)))
            if (dtlongn(1) > timeleft) dtlongn(1) = timeleft
         else
            dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
         endif
         nnacoust(ifm) = nacoust
         sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
         sspct = min(sspct,  &
            .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
enddo

return
END SUBROUTINE dtset
