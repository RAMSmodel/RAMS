!##############################################################################
Subroutine grid_setup (num)

use mem_grid
use node_mod

implicit none

integer :: num

integer :: ifm

if (num == 1) then

   CALL gridinit ()
   CALL gridset ()

else  

   do ifm = 1,ngrids
      CALL newgrid (ifm)

      if (ihtran <= 1) then

         CALL polarst (mxp, myp, i0, j0                        &
            ,grid_g(ifm)%glat  (1,1) ,grid_g(ifm)%glon  (1,1)  &
            ,grid_g(ifm)%fmapu (1,1) ,grid_g(ifm)%fmapv (1,1)  &
            ,grid_g(ifm)%fmapt (1,1) ,grid_g(ifm)%fmapm (1,1)  &
            ,grid_g(ifm)%fmapui(1,1) ,grid_g(ifm)%fmapvi(1,1)  &
            ,grid_g(ifm)%fmapti(1,1) ,grid_g(ifm)%fmapmi(1,1)  )

         CALL grdspc (mxp, myp, i0, j0                        &
            ,grid_g(ifm)%dxu  (1,1)  ,grid_g(ifm)%dxv  (1,1)  &
            ,grid_g(ifm)%dxt  (1,1)  ,grid_g(ifm)%dxm  (1,1)  &
            ,grid_g(ifm)%dyu  (1,1)  ,grid_g(ifm)%dyv  (1,1)  &
            ,grid_g(ifm)%dyt  (1,1)  ,grid_g(ifm)%dym  (1,1)  &
            ,grid_g(ifm)%fmapu(1,1)  ,grid_g(ifm)%fmapv(1,1)  &
            ,grid_g(ifm)%fmapt(1,1)  ,grid_g(ifm)%fmapm(1,1)  )
      
      else
         print*,'Unknown IHTRAN value:',ihtran
         stop 'grid_setup: bad IHTRAN'
        
      endif


! Define transformation Jacobians for all grids

      CALL fill_toptuvm (mxp, myp, i0, j0                 &
         ,grid_g(ifm)%topt (1,1) ,grid_g(ifm)%topu (1,1)  &
         ,grid_g(ifm)%topv (1,1) ,grid_g(ifm)%topm (1,1))

      CALL transfm (mxp, myp                            &
         ,grid_g(ifm)%topt(1,1) ,grid_g(ifm)%topu(1,1)  &
         ,grid_g(ifm)%topv(1,1) ,grid_g(ifm)%topm(1,1)  &
         ,grid_g(ifm)%rtgt(1,1) ,grid_g(ifm)%rtgu(1,1)  &
         ,grid_g(ifm)%rtgv(1,1) ,grid_g(ifm)%rtgm(1,1)  &
         ,grid_g(ifm)%f13u(1,1) ,grid_g(ifm)%f13v(1,1)  &
         ,grid_g(ifm)%f13t(1,1) ,grid_g(ifm)%f13m(1,1)  &
         ,grid_g(ifm)%f23u(1,1) ,grid_g(ifm)%f23v(1,1)  &
         ,grid_g(ifm)%f23t(1,1) ,grid_g(ifm)%f23m(1,1)  &
         ,grid_g(ifm)%dxu (1,1) ,grid_g(ifm)%dxv (1,1)  &
         ,grid_g(ifm)%dxt (1,1) ,grid_g(ifm)%dxm (1,1)  &
         ,grid_g(ifm)%dyu (1,1) ,grid_g(ifm)%dyv (1,1)  &
         ,grid_g(ifm)%dyt (1,1) ,grid_g(ifm)%dym (1,1)  )

   enddo
endif

return
END SUBROUTINE grid_setup

!##############################################################################
Subroutine gridinit ()

use mem_grid

implicit none

!     +----------------------------------------------------------------
!     !    initialize the domain sizes for all nest levels.
!     +----------------------------------------------------------------
!
!          grid sizes:  the user has specified all nnxp, nnyp, nnzp,
!                       nzg, nzs, and npatch.  fill other arrays that are a
!                       func of these.

jdim = 1

if (nnyp(1) == 1) jdim = 0
do ngrid = 1,ngrids

!         x - direction

   nnx(ngrid) = nnxp(ngrid) - 1
   nnx1(ngrid) = nnxp(ngrid) - 2
   nnx2(ngrid) = nnxp(ngrid) - 3

!         y - direction

  if (jdim == 1) then
    nny(ngrid) = nnyp(ngrid) - 1
    nny1(ngrid) = nnyp(ngrid) - 2
    nny2(ngrid) = nnyp(ngrid) - 3
  else
    nny(ngrid) = 1
    nny1(ngrid) = 1
    nny2(ngrid) = 1
  endif

!         z - direction

  nnz(ngrid) = nnzp(ngrid) - 1
  nnz1(ngrid) = nnzp(ngrid) - 2

enddo

do ngrid = 1,ngrids
  nnxyzp(ngrid) = nnxp(ngrid) * nnyp(ngrid) * nnzp(ngrid)
  nnxysp(ngrid) = nnxp(ngrid) * nnyp(ngrid) * (nzg+nzs+3) * npatch
  nnxyp(ngrid) = nnxp(ngrid) * nnyp(ngrid)
enddo

return
END SUBROUTINE gridinit

!##############################################################################
Subroutine polarst (n2,n3,ioff,joff,glat,glon,fmapu,fmapv,fmapt,fmapm  &
                  ,fmapui,fmapvi,fmapti,fmapmi)

use mem_grid
use node_mod
use rconstants

implicit none

integer :: n2,n3,ioff,joff
real, dimension(n2,n3) :: glat,glon,fmapu,fmapv,fmapt,fmapm  &
         ,fmapui,fmapvi,fmapti,fmapmi

integer :: i,j
real :: c1,xm2,xt2,ym2,yt2
!  Calculates map factors and inverse map factors at u,v,t,m-points and
!  geographical lat/lon at t-points for a given polar stereographic grid

c1 = (2. * erad) ** 2
do j = 1,n3
   do i = 1,n2
      xm2 = xm(i+ioff) * xm(i+ioff)
      xt2 = xt(i+ioff) * xt(i+ioff)
      ym2 = ym(j+joff) * ym(j+joff)
      yt2 = yt(j+joff) * yt(j+joff)

      fmapt(i,j) = 1. + (xt2 + yt2) / c1
      fmapu(i,j) = 1. + (xm2 + yt2) / c1
      fmapv(i,j) = 1. + (xt2 + ym2) / c1
      fmapm(i,j) = 1. + (xm2 + ym2) / c1

      fmapui(i,j) = 1.0 / fmapu(i,j)
      fmapvi(i,j) = 1.0 / fmapv(i,j)
      fmapti(i,j) = 1.0 / fmapt(i,j)
      fmapmi(i,j) = 1.0 / fmapm(i,j)

      CALL xy_ll (glat(i,j),glon(i,j),polelat,polelon  &
         ,xt(i+ioff),yt(j+joff))

!       write(6,344)i,j,fmapt(i,j),fmapm(i,j),glat(i,j),glon(i,j)
! 344   format('polst:i,j,fmt,fmm,glt,gln',2i4,4e12.3)

   enddo
enddo

if (ihtran == 0) then
   CALL ae0 (n2*n3,fmapu,1.)
   CALL ae0 (n2*n3,fmapv,1.)
   CALL ae0 (n2*n3,fmapt,1.)
   CALL ae0 (n2*n3,fmapm,1.)
   CALL ae0 (n2*n3,fmapui,1.)
   CALL ae0 (n2*n3,fmapvi,1.)
   CALL ae0 (n2*n3,fmapti,1.)
   CALL ae0 (n2*n3,fmapmi,1.)
endif

return
END SUBROUTINE polarst

!##############################################################################
Subroutine grdspc (n2,n3,ioff,joff,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
                 ,fmapu,fmapv,fmapt,fmapm)

use mem_grid

implicit none

integer :: n2,n3,ioff,joff
real, dimension(n2,n3) :: dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
                         ,fmapu,fmapv,fmapt,fmapm

integer :: i,j

do j = 1,n3
   do i = 1,n2-1
      dxu(i,j) = fmapu(i,j) / (xtn(i+ioff+1,ngrid)-xtn(i+ioff,ngrid))
      dxm(i,j) = fmapm(i,j) / (xtn(i+ioff+1,ngrid)-xtn(i+ioff,ngrid))
   enddo
   dxu(n2,j)=dxu(n2-1,j)*fmapu(n2,j)/fmapu(n2-1,j)
   dxm(n2,j)=dxm(n2-1,j)*fmapm(n2,j)/fmapm(n2-1,j)
   do i = 2,n2
      dxv(i,j)=fmapv(i,j)/(xmn(i+ioff,ngrid)-xmn(i+ioff-1,ngrid))
      dxt(i,j)=fmapt(i,j)/(xmn(i+ioff,ngrid)-xmn(i+ioff-1,ngrid))
   enddo
   dxv(1,j)=dxv(2,j)*fmapv(1,j)/fmapv(2,j)
   dxt(1,j)=dxt(2,j)*fmapt(1,j)/fmapt(2,j)
enddo

if (jdim == 1) then
   do i = 1,n2
      do j = 1,n3-1
         dyv(i,j)=fmapv(i,j)/(ytn(j+joff+1,ngrid)-ytn(j+joff,ngrid))
         dym(i,j)=fmapm(i,j)/(ytn(j+joff+1,ngrid)-ytn(j+joff,ngrid))
      enddo
      dyv(i,n3)=dyv(i,n3-1)*fmapv(i,n3)/fmapv(i,n3-1)
      dym(i,n3)=dym(i,n3-1)*fmapm(i,n3)/fmapm(i,n3-1)
      do j = 2,n3
         dyu(i,j)=fmapu(i,j)/(ymn(j+joff,ngrid)-ymn(j+joff-1,ngrid))
         dyt(i,j)=fmapt(i,j)/(ymn(j+joff,ngrid)-ymn(j+joff-1,ngrid))
      enddo
      dyu(i,1)=dyu(i,2)*fmapu(i,1)/fmapu(i,2)
      dyt(i,1)=dyt(i,2)*fmapt(i,1)/fmapt(i,2)
   enddo
else
   do i=1,n2
      do j=1,n3
         dyu(i,j)=1./deltaxn(ngrid)
         dyv(i,j)=1./deltaxn(ngrid)
         dyt(i,j)=1./deltaxn(ngrid)
         dym(i,j)=1./deltaxn(ngrid)
      enddo
   enddo
endif

return
END SUBROUTINE grdspc

!##############################################################################
Subroutine fill_toptuvm (n2,n3,ioff,joff,topt,topu,topv,topm)

use mem_grid

implicit none

integer :: n2,n3,ioff,joff
real, dimension(n2,n3) :: topt,topu,topv,topm

integer :: i,j
real :: terdev

terdev = 0.
do j = 1,n3
   do i = 1,n2-1
      topu(i,j) = topt(i,j) + (topt(i+1,j) - topt(i,j))  &
                * (xm(i+ioff) - xt(i+ioff)) / (xt(i+ioff+1) - xt(i+ioff))
      terdev = max(terdev,abs(topt(i,j)))
   enddo
   topu(n2,j) = topt(n2,j) + (topt(n2,j) - topt(n2-1,j))  &
               * (xm(n2+ioff) - xt(n2+ioff)) / (xt(n2+ioff) - xt(n2+ioff-1))
enddo

if (terdev < 1.e-6) then
   itopo = 0
else
   itopo = 1
endif

if (jdim == 1) then
   do i = 1,n2
      do j = 1,n3-1
         topv(i,j) = topt(i,j) + (topt(i,j+1) - topt(i,j))  &
                   * (ym(j+joff) - yt(j+joff)) / (yt(j+joff+1) - yt(j+joff))
         topm(i,j) = topu(i,j) + (topu(i,j+1) - topu(i,j))  &
                   * (ym(j+joff) - yt(j+joff)) / (yt(j+joff+1) - yt(j+joff))
      enddo
      topv(i,n3) = topt(i,n3) + (topt(i,n3) - topt(i,n3-1))  &
                  * (ym(n3+joff) - yt(n3+joff)) / (yt(n3+joff) - yt(n3+joff-1))
      topm(i,n3) = topu(i,n3) + (topu(i,n3) - topu(i,n3-1))  &
                  * (ym(n3+joff) - yt(n3+joff)) / (yt(n3+joff) - yt(n3+joff-1))
   enddo
else
   do j = 1,n3
      do i = 1,n2
         topv(i,j) = topt(i,j)
         topm(i,j) = topu(i,j)
      enddo
   enddo
endif

return
END SUBROUTINE fill_toptuvm

!##############################################################################
Subroutine transfm (n2,n3,topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
                  ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
                  ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym)

use mem_grid

implicit none

integer :: n2,n3
real, dimension(n2,n3) :: topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
         ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
         ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym

integer :: iztflag=0,i,j,k

!     this routine computes the coordinate transformation constants
!     based on the topographical values of TOPT.

ztop = zmn(nnzp(1)-1,1)
do k = 1,nzp
  htn(k,ngrid) = zt(k) / ztop - 1.
  hwn(k,ngrid) = zm(k) / ztop - 1.
enddo
do k = 1,nzp
  ht2n(k,ngrid) = .5 * htn(k,ngrid)
  ht4n(k,ngrid) = .25 * htn(k,ngrid)
  hw2n(k,ngrid) = .5 * hwn(k,ngrid)
  hw4n(k,ngrid) = .25 * hwn(k,ngrid)
enddo
do k = 1,nzp
  ht(k)  = htn(k,ngrid)
  hw(k)  = hwn(k,ngrid)
  ht2(k) = ht2n(k,ngrid)
  ht4(k) = ht4n(k,ngrid)
  hw2(k) = hw2n(k,ngrid)
  hw4(k) = hw4n(k,ngrid)
enddo

do j = 1,n3
   do i = 1,n2
      rtgt(i,j) = 1. - topt(i,j) / ztop
      rtgu(i,j) = 1. - topu(i,j) / ztop
      rtgv(i,j) = 1. - topv(i,j) / ztop
      rtgm(i,j) = 1. - topm(i,j) / ztop
      if (topt(i,j) > .5 * ztop) then
         print*, 'Terrain height is over half the model domain'
         print*, 'height.  Model will stop here to prevent this.'
         print*, 'ngrid, i, j, topt, ztop = ',ngrid,i,j,topt(i,j),ztop
         iztflag = 1
      endif
   enddo
enddo

!Saleeby: Write out some grid navigation info to a file
!open(unit=34,file='topo.rtgt.txt')
!do j = 1,n3
!   do i = 1,n2
!      WRITE(34,100)'ngrid,i,j,ztop,topt,rtgt',ngrid,i,j,ztop,topt(i,j),rtgt(i,j)
!      100 FORMAT(A,3I5,2X,F9.2,2X,F11.3,2X,F9.5)
!   enddo
!enddo
!close(34)

if (iztflag == 1) stop 'topt/ztop'

do j = 1,n3
   do i = 2,n2
      f13t(i,j) = (topu(i,j) - topu(i-1,j)) * dxt(i,j) / rtgt(i,j)
      f13v(i,j) = (topm(i,j) - topm(i-1,j)) * dxv(i,j) / rtgv(i,j)
   enddo
   do i = 1,n2-1
      f13u(i,j) = (topt(i+1,j) - topt(i,j)) * dxu(i,j) / rtgu(i,j)
      f13m(i,j) = (topv(i+1,j) - topv(i,j)) * dxm(i,j) / rtgm(i,j)
   enddo
   f13t(1,j)  = f13u(1,j)
   f13v(1,j)  = f13m(1,j)
   f13u(n2,j) = f13t(n2,j)
   f13m(n2,j) = f13v(n2,j)
enddo

do i = 1,n2
   do j = 2,n3
      f23t(i,j) = (topv(i,j) - topv(i,j-jdim)) * dyt(i,j) / rtgt(i,j)
      f23u(i,j) = (topm(i,j) - topm(i,j-jdim)) * dyu(i,j) / rtgu(i,j)
   enddo
   do j = 1,n3-1
      f23v(i,j) = (topt(i,j+jdim) - topt(i,j)) * dyv(i,j) / rtgv(i,j)
      f23m(i,j) = (topu(i,j+jdim) - topu(i,j)) * dym(i,j) / rtgm(i,j)
   enddo
   if (jdim == 1) then
      f23t(i,1)  = f23v(i,1)
      f23u(i,1)  = f23m(i,1)
      f23v(i,n3) = f23t(i,n3)
      f23m(i,n3) = f23u(i,n3)
   endif
enddo

return
END SUBROUTINE transfm

!##############################################################################
Subroutine newgrid (ngr)

use mem_grid
use node_mod

implicit none

integer :: ngr,i,j,k

!     +----------------------------------------------------------------
!     !    Fill the single and 1D variables that the rest of the model
!     !      uses from the nest arrays and change grid level in the I/O.
!     +----------------------------------------------------------------

ngrid = ngr

!         grid point references

!         x - direction

nxp=nnxp(ngr)
nx=nnx(ngr)
nx1=nnx1(ngr)
nx2=nnx2(ngr)
!
!         y - direction
!
nyp=nnyp(ngr)
ny=nny(ngr)
ny1=nny1(ngr)
ny2=nny2(ngr)
!
!         z - direction
!
nzp=nnzp(ngr)
nzpp=nzp+1
nz=nnz(ngr)
nz1=nnz1(ngr)


nxyzp=nnxyzp(ngr)
nxysp=nnxysp(ngr)
nxyp=nnxyp(ngr)

!          grid spacings

deltax=deltaxn(ngr)

do i=1,nxp
  xt(i)=xtn(i,ngr)
  xm(i)=xmn(i,ngr)
enddo
do j=1,nyp
  yt(j)=ytn(j,ngr)
  ym(j)=ymn(j,ngr)
enddo

deltaz=zmn(2,ngr)-zmn(1,ngr)

do k=1,nzp
  zt(k)=ztn(k,ngr)
  zm(k)=zmn(k,ngr)
  dzm(k)=dzmn(k,ngr)
  dzt(k)=dztn(k,ngr)
  dzm2(k)=dzm2n(k,ngr)
  dzt2(k)=dzt2n(k,ngr)
  ht(k)=htn(k,ngr)
  ht2(k)=ht2n(k,ngr)
  ht4(k)=ht4n(k,ngr)
  hw(k)=hwn(k,ngr)
  hw2(k)=hw2n(k,ngr)
  hw4(k)=hw4n(k,ngr)
enddo

!         option flags

nsttop=nnsttop(ngr)
nstbot=nnstbot(ngr)

!         timesteps

dtlt=dtlongn(ngr)
dtlv=2.*dtlt

!        node gridpoint info

mxp=mmxp(ngr)
myp=mmyp(ngr)
mzp=mmzp(ngr)
ia=mia(ngr)
iz=miz(ngr)
ja=mja(ngr)
jz=mjz(ngr)
i0=mi0(ngr)
j0=mj0(ngr)
ibcon=mibcon(ngr)

mxyzp=mmxyzp(ngr)
mxysp=mmxysp(ngr)
mxyp=mmxyp(ngr)

return
END SUBROUTINE newgrid
