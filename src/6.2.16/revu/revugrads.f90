!##############################################################################
Subroutine gradr (n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt  &
                 ,xm,xt,ym,yt,zm,zt,deltax,dzm,dzt,vctr1,vctr2,ztop  &
                 ,jdim,ihtran,polelat,polelon) 

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,ihtran,jdim,jaa,jzz,i,j,k &
    ,isizexy,isizeall,iglat,iglon,ifmapu,ifmapv,ifmapt,ifmapm &
    ,ifmapui,ifmapvi,ifmapti,ifmapmi,idxu,idxv,idxt,idxm,idyu &
    ,idyv,idyt,idym,itopu,itopv,itopm,irtgt,irtgu,irtgv,irtgm &
    ,if13u,if13v,if13t,if13m,if23u,if23v,if23t,if23m
real :: polelat,polelon,deltax,ztop
real :: vc3da(n1,n2,n3),vc3db(n1,n2,n3),topt(n1,n2)  &
       ,xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),dzm(*),dzt(*)  &
       ,vctr1(*),vctr2(*)
character(len=*) :: dir,gpnt
character(len=6) :: optyp
real, allocatable :: e(:),hw(:),ht(:)

if(allocated(ht)) deallocate (ht); allocate (ht(n3))
if(allocated(hw)) deallocate (hw); allocate (hw(n3))

optyp = 'gradnt'

jaa = ja
jzz = jz
if(jdim.eq.0) then
   jaa = 1
   jzz = 1
endif

do k=1,n3
   do j=1,n2
      do i=1,n1
         vc3db(i,j,k)=0.
      enddo
   enddo
enddo

!print*,'in grad-',dir,gpnt,optyp,'   !!!!!!!!!!!!!!!!'

isizexy = n1 * n2
isizeall = isizexy * 33

if(allocated(e)) deallocate (e)
allocate (e(isizeall))

iglat   = 1
iglon   = iglat   + isizexy
ifmapu  = iglon   + isizexy
ifmapv  = ifmapu  + isizexy
ifmapt  = ifmapv  + isizexy
ifmapm  = ifmapt  + isizexy
ifmapui = ifmapm  + isizexy
ifmapvi = ifmapui + isizexy
ifmapti = ifmapvi + isizexy
ifmapmi = ifmapti + isizexy
idxu    = ifmapmi + isizexy
idxv    = idxu    + isizexy
idxt    = idxv    + isizexy
idxm    = idxt    + isizexy
idyu    = idxm    + isizexy
idyv    = idyu    + isizexy
idyt    = idyv    + isizexy
idym    = idyt    + isizexy
itopu   = idym    + isizexy
itopv   = itopu   + isizexy
itopm   = itopv   + isizexy
irtgt   = itopm   + isizexy
irtgu   = irtgt   + isizexy
irtgv   = irtgu   + isizexy
irtgm   = irtgv   + isizexy
if13u   = irtgm   + isizexy
if13v   = if13u   + isizexy
if13t   = if13v   + isizexy
if13m   = if13t   + isizexy
if23u   = if13m   + isizexy
if23v   = if23u   + isizexy
if23t   = if23v   + isizexy
if23m   = if23t   + isizexy

if(ihtran.eq.1)then
   CALL revu_polarst (N1,N2,E(IGLAT),E(IGLON),E(IFMAPU),E(IFMAPV)  &
                     ,E(IFMAPT),E(IFMAPM),E(IFMAPUI),E(IFMAPVI)  &
                     ,E(IFMAPTI),E(IFMAPMI),xm,xt,ym,yt,polelat,polelon)

else
   CALL ae0 (N1*N2,E(IFMAPU),1.)
   CALL ae0 (N1*N2,E(IFMAPV),1.)
   CALL ae0 (N1*N2,E(IFMAPT),1.)
   CALL ae0 (N1*N2,E(IFMAPM),1.)
   CALL ae0 (N1*N2,E(IFMAPUI),1.)
   CALL ae0 (N1*N2,E(IFMAPVI),1.)
   CALL ae0 (N1*N2,E(IFMAPTI),1.)
   CALL ae0 (N1*N2,E(IFMAPMI),1.)
ENDIF

CALL revu_grdspc (N1,N2,E(IDXU),E(IDXV),E(IDXT),E(IDXM)  &
                 ,E(IDYU),E(IDYV),E(IDYT),E(IDYM)  &
                 ,E(IFMAPU),E(IFMAPV),E(IFMAPT),E(IFMAPM)  &
                 ,xm,xt,ym,yt,jdim,deltax)

CALL revu_transfm (n1,n2,n3,TOPT,E(ITOPU),E(ITOPV),E(ITOPM)  &
                  ,E(IRTGT),E(IRTGU),E(IRTGV),E(IRTGM),E(IF13U),E(IF13V)  &
                  ,E(IF13T),E(IF13M),E(IF23U),E(IF23V),E(IF23T)  &
                  ,E(IF23M),E(IDXU),E(IDXV)  &
                  ,E(IDXT),E(IDXM),E(IDYU),E(IDYV),E(IDYT),E(IDYM)  &
                  ,xm,xt,ym,yt,zm,zt,jdim,ztop,hw,ht)

if(dir.eq.'xdir')then

   if(gpnt.eq.'upnt')then
      CALL gradxur (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgu),e(irtgt),e(idxt),dzt  &
                   ,e(if13t),hw,vctr2,'t')
   elseif(gpnt.eq.'vpnt')then
      CALL gradxtr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgv),e(irtgm),e(idxm),dzt  &
                   ,e(if13m),hw,vctr2,'t')
   elseif(gpnt.eq.'wpnt')then
      CALL gradxtr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgt),e(irtgu),e(idxu),dzm  &
                  ,e(if13u),ht,vctr2,'w')
   elseif(gpnt.eq.'tpnt')then
      CALL gradxtr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgt),e(irtgu),e(idxu),dzt  &
                   ,e(if13u),hw,vctr2,'t')
   elseif(gpnt.eq.'npnt')then
      CALL gradxtr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgv),e(irtgm),e(idxm),dzm  &
                   ,e(if13m),ht,vctr2,'w')
   elseif(gpnt.eq.'opnt')then
      CALL gradxur (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgu),e(irtgt),e(idxt),dzm  &
                   ,e(if13t),ht,vctr2,'w')
   elseif(gpnt.eq.'ppnt')then
      CALL gradxur (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgm),e(irtgv),e(idxv),dzt  &
                   ,e(if13v),hw,vctr2,'t')
   elseif(gpnt.eq.'mpnt')then
      CALL gradxur (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgm),e(irtgv),e(idxv),dzm  &
                   ,e(if13v),ht,vctr2,'w')
   endif
   
elseif(dir.eq.'ydir')then

   if(gpnt.eq.'upnt')then
      CALL gradytr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgu),e(irtgm),e(idym),dzt  &
                   ,e(if23m),hw,vctr2,'t',jdim)
   elseif(gpnt.eq.'vpnt')then
      CALL gradyvr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgv),e(irtgt),e(idyt),dzt  &
                   ,e(if23t),hw,vctr2,'t',jdim)
   elseif(gpnt.eq.'wpnt')then
      CALL gradytr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgt),e(irtgv),e(idyv),dzm  &
                   ,e(if23v),ht,vctr2,'w',jdim)
   elseif(gpnt.eq.'tpnt')then
      CALL gradytr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgt),e(irtgv),e(idyv),dzt  &
                   ,e(if23v),hw,vctr2,'t',jdim)
   elseif(gpnt.eq.'npnt')then
      CALL gradyvr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgv),e(irtgt),e(idyt),dzm  &
                   ,e(if23t),ht,vctr2,'w',jdim)
   elseif(gpnt.eq.'opnt')then
      CALL gradytr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgu),e(irtgm),e(idym),dzm  &
                   ,e(if23m),ht,vctr2,'w',jdim)
   elseif(gpnt.eq.'ppnt')then
      CALL gradyvr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgm),e(irtgu),e(idyu),dzt  &
                   ,e(if23u),hw,vctr2,'t',jdim)
   elseif(gpnt.eq.'mpnt')then
      CALL gradyvr (n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1  &
                   ,e(irtgm),e(irtgu),e(idyu),dzm  &
                   ,e(if23u),ht,vctr2,'w',jdim)
   endif
   
elseif(dir.eq.'zdir')then

   if(gpnt.eq.'upnt')then
     CALL gradztr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgu),dzm)
   elseif(gpnt.eq.'vpnt')then
     CALL gradztr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgv),dzm)
   elseif(gpnt.eq.'wpnt')then
     CALL gradzwr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgt),dzt)
   elseif(gpnt.eq.'tpnt')then
     CALL gradztr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgt),dzm)
   elseif(gpnt.eq.'npnt')then
     CALL gradzwr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgv),dzt)
   elseif(gpnt.eq.'opnt')then
     CALL gradzwr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgu),dzt)
   elseif(gpnt.eq.'ppnt')then
     CALL gradztr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgm),dzm)
   elseif(gpnt.eq.'mpnt')then
     CALL gradzwr (n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgm),dzt)
   endif
   
endif

deallocate (e,hw,ht)

return
END SUBROUTINE gradr

!##############################################################################
Subroutine gradxur (n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da  &
   ,rtge,rtgc,dx,dz,fq,hq,hq4,lev)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz,vc1da,hq,hq4
real, dimension(n1,n2) :: rtgc,rtge,dx,fq
character(len=1) :: lev
character(len=6) :: optyp

if(optyp.eq.'gradnt')then
   do j = ja,jz
      do i = ia,iz
         do k = 1,n3
            vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j)  &
               - vc3da(i-1,j,k) * rtge(i-1,j))  &
               * dx(i,j) / rtgc(i,j)
         enddo
      enddo
   enddo
endif

if(optyp.ne.'divstr')then
   if(lev.eq.'w')then
      do k = 1,n3
         hq4(k) = 0.25 * hq(k)
      enddo
   else
      do k = 2,n3
         hq4(k) = 0.25 * hq(k-1)
      enddo
   endif

   do j = ja,jz
      do i = ia,iz
         do k = 2,n3
            vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1)  &
               + vc3da(i-1,j,k) + vc3da(i-1,j,k-1))
         enddo
         do k = 2,n3-1
            vc3db(i,j,k) = vc3db(i,j,k)  &
               + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
         enddo
         vc3db(i,j,1) = vc3db(i,j,2)
         if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
         if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
      enddo
   enddo
endif

return
END SUBROUTINE gradxur

!##############################################################################
Subroutine gradxtr (n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da  &
              ,rtge,rtgc,dx,dz,fq,hq,hq4,lev)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz,vc1da,hq,hq4
real, dimension(n1,n2) :: rtgc,rtge,dx,fq
character(len=1) :: lev
character(len=6) :: optyp

if(optyp.eq.'gradnt')then
   do j = ja,jz
      do i = ia,iz
         do k = 1,n3
            vc3db(i,j,k) = (vc3da(i+1,j,k) * rtge(i+1,j)  &
               - vc3da(i,j,k) * rtge(i,j))  &
               * dx(i,j) / rtgc(i,j)
         enddo
      enddo
   enddo
endif

if(optyp.ne.'divstr')then
   if(lev.eq.'w')then
      do k = 1,n3
         hq4(k) = 0.25 * hq(k)
      enddo
   else
      do k = 2,n3
         hq4(k) = 0.25 * hq(k-1)
      enddo
   endif

   do j = ja,jz
      do i = ia,iz
         do k = 2,n3
            vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1)  &
               + vc3da(i+1,j,k) + vc3da(i+1,j,k-1))
         enddo
         do k = 2,n3-1
            vc3db(i,j,k) = vc3db(i,j,k)  &
               + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
         enddo
         vc3db(i,j,1) = vc3db(i,j,2)
         if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
         if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
      enddo
   enddo
endif

return
END SUBROUTINE gradxtr

!##############################################################################
Subroutine gradyvr (n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da  &
              ,rtge,rtgc,dy,dz,fq,hq,hq4,lev,jd)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k,jd
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz,vc1da,hq,hq4
real, dimension(n1,n2) :: rtgc,rtge,dy,fq
character(len=1) :: lev
character(len=6) :: optyp

if(optyp.eq.'gradnt')then
   do j = ja,jz
      do i = ia,iz
         do k = 1,n3
            vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j)  &
               - vc3da(i,j-jd,k) * rtge(i,j-jd))  &
               * dy(i,j) / rtgc(i,j)
         enddo
      enddo
   enddo
endif

if(optyp.ne.'divstr')then
   if(lev.eq.'w')then
      do k = 1,n3
         hq4(k) = 0.25 * hq(k)
      enddo
   else
      do k = 2,n3
         hq4(k) = 0.25 * hq(k-1)
      enddo
   endif

   do j = ja,jz
      do i = ia,iz
         do k = 2,n3
            vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1)  &
               + vc3da(i,j-jd,k) + vc3da(i,j-jd,k-1))
         enddo
         do k = 2,n3-1
            vc3db(i,j,k) = vc3db(i,j,k)  &
               + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
         enddo
         vc3db(i,j,1) = vc3db(i,j,2)
         if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
         if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
      enddo
   enddo
endif

return
END SUBROUTINE gradyvr

!##############################################################################
Subroutine gradytr (n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da  &
              ,rtge,rtgc,dy,dz,fq,hq,hq4,lev,jd)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k,jd
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz,vc1da,hq,hq4
real, dimension(n1,n2) :: rtgc,rtge,dy,fq
character(len=1) :: lev
character(len=6) :: optyp

if(optyp.eq.'gradnt')then
   do j = ja,jz
      do i = ia,iz
         do k = 1,n3
            vc3db(i,j,k) = (vc3da(i,j+jd,k) * rtge(i,j+jd)  &
               - vc3da(i,j,k) * rtge(i,j))  &
               * dy(i,j) / rtgc(i,j)
         enddo
      enddo
   enddo
endif

if(optyp.ne.'divstr')then
   if(lev.eq.'w')then
      do k = 1,n3
         hq4(k) = 0.25 * hq(k)
      enddo
   else
      do k = 2,n3
         hq4(k) = 0.25 * hq(k-1)
      enddo
   endif

   do j = ja,jz
      do i = ia,iz
         do k = 2,n3
            vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1)  &
               + vc3da(i,j+jd,k) + vc3da(i,j+jd,k-1))
         enddo
         do k = 2,n3-1
            vc3db(i,j,k) = vc3db(i,j,k)  &
               + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
         enddo
         vc3db(i,j,1) = vc3db(i,j,2)
         if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
         if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
      enddo
   enddo
endif

return
END SUBROUTINE gradytr

!##############################################################################
Subroutine gradzwr (n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz
real, dimension(n1,n2) :: rtgc

do j = ja,jz
   do i = ia,iz
      do k = 2,n3
         vc3db(i,j,k) = (vc3da(i,j,k) - vc3da(i,j,k-1)) * dz(k) / rtgc(i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE gradzwr

!##############################################################################
Subroutine gradztr (n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

! this routine which computes a component of the
! gradient of vc3da and stores it in vc3db.

implicit none

integer :: n1,n2,n3,ia,iz,ja,jz,i,j,k
real, dimension(n1,n2,n3) :: vc3da,vc3db
real, dimension(n3) :: dz
real, dimension(n1,n2) :: rtgc

do j = ja,jz
   do i = ia,iz
      do k = 1,n3-1
         vc3db(i,j,k) = (vc3da(i,j,k+1) - vc3da(i,j,k)) * dz(k) / rtgc(i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE gradztr

!##############################################################################
Subroutine revu_polarst (N1,N2,GLAT,GLON,FMAPU,FMAPV,FMAPT,FMAPM  &
                        ,FMAPUI,FMAPVI,FMAPTI,FMAPMI,xm,xt,ym,yt  &
                        ,polelat,polelon)

implicit none

integer :: N1,N2,I,J
real :: c1,erad,xm2,xt2,ym2,yt2,polelat,polelon
real :: GLAT(N1,N2),GLON(N1,N2),FMAPU(N1,N2),FMAPUI(N1,N2)  &
       ,FMAPV(N1,N2),FMAPT(N1,N2),FMAPM(N1,N2),FMAPVI(N1,N2)  &
       ,FMAPTI(N1,N2),FMAPMI(N1,N2),xm(*),xt(*),ym(*),yt(*)

! Calculates map factors and inverse map factors at u,v,t,m-points and
! geographical lat/lon at t-points for a given polar stereographic grid

ERAD=6367000.
c1 = (2. * erad) ** 2
DO J = 1,N2
   DO I = 1,N1
      xm2 = xm(i) * xm(i)
      xt2 = xt(i) * xt(i)
      ym2 = ym(j) * ym(j)
      yt2 = yt(j) * yt(j)

      FMAPT(I,J) = 1. + (xt2 + yt2) / c1
      FMAPU(I,J) = 1. + (xm2 + yt2) / c1
      FMAPV(I,J) = 1. + (xt2 + ym2) / c1
      FMAPM(I,J) = 1. + (xm2 + ym2) / c1

      FMAPUI(I,J) = 1.0 / FMAPU(I,J)
      FMAPVI(I,J) = 1.0 / FMAPV(I,J)
      FMAPTI(I,J) = 1.0 / FMAPT(I,J)
      FMAPMI(I,J) = 1.0 / FMAPM(I,J)

      CALL xy_ll (GLAT(I,J),GLON(I,J),polelat,polelon,XT(I),YT(J))

      !write(6,344)i,j,fmapt(i,j),fmapm(i,j),glat(i,j),glon(i,j),XT(I),YT(J)
      !344 format('polst:i,j,fmt,fmm,glt,gln,x,y',2i4,6e15.6)

   ENDDO
ENDDO

return
END SUBROUTINE revu_polarst

!##############################################################################
Subroutine revu_grdspc (N1,N2,DXU,DXV,DXT,DXM,DYU,DYV,DYT,DYM  &
                       ,FMAPU,FMAPV,FMAPT,FMAPM,xm,xt,ym,yt,jdim,deltax)

implicit none

integer :: N1,N2,I,J,jdim
real :: deltax
real :: DXU(N1,N2),DXV(N1,N2),DXT(N1,N2),DXM(N1,N2)  &
       ,DYU(N1,N2),DYV(N1,N2),DYT(N1,N2),DYM(N1,N2)  &
       ,FMAPU(N1,N2),FMAPV(N1,N2),FMAPT(N1,N2),FMAPM(N1,N2)  &
       ,xm(*),xt(*),ym(*),yt(*)

DO J=1,N2
   DO I=1,N1-1
      DXU(I,J)=FMAPU(I,J)/(XT(I+1)-XT(I))
      DXM(I,J)=FMAPM(I,J)/(XT(I+1)-XT(I))
   ENDDO
   DXU(N2,J)=DXU(N2-1,J)*FMAPU(N2,J)/FMAPU(N2-1,J)
   DXM(N2,J)=DXM(N2-1,J)*FMAPM(N2,J)/FMAPM(N2-1,J)
   DO I=2,N1
      DXV(I,J)=FMAPV(I,J)/(XM(I)-XM(I-1))
      DXT(I,J)=FMAPT(I,J)/(XM(I)-XM(I-1))
   ENDDO
   DXV(1,J)=DXV(2,J)*FMAPV(1,J)/FMAPV(2,J)
   DXT(1,J)=DXT(2,J)*FMAPT(1,J)/FMAPT(2,J)
ENDDO

IF(JDIM.EQ.1)THEN
   DO I=1,N1
      DO J=1,N2-1
         DYV(I,J)=FMAPV(I,J)/(YT(J+1)-YT(J))
         DYM(I,J)=FMAPM(I,J)/(YT(J+1)-YT(J))
      ENDDO
      DYV(I,N2)=DYV(I,N2-1)*FMAPV(I,N2)/FMAPV(I,N2-1)
      DYM(I,N2)=DYM(I,N2-1)*FMAPM(I,N2)/FMAPM(I,N2-1)
      DO J=2,N2
         DYU(I,J)=FMAPU(I,J)/(YM(J)-YM(J-1))
         DYT(I,J)=FMAPT(I,J)/(YM(J)-YM(J-1))
      ENDDO
      DYU(I,1)=DYU(I,2)*FMAPU(I,1)/FMAPU(I,2)
      DYT(I,1)=DYT(I,2)*FMAPT(I,1)/FMAPT(I,2)
   ENDDO
ELSE
   DO I=1,N1
      DO J=1,N2
         DYU(I,J)=1./DELTAX
         DYV(I,J)=1./DELTAX
         DYT(I,J)=1./DELTAX
         DYM(I,J)=1./DELTAX
      ENDDO
   ENDDO
ENDIF

return
END SUBROUTINE revu_grdspc

!##############################################################################
Subroutine revu_transfm (N1,N2,N3,TOPT,TOPU,TOPV,TOPM,RTGT,RTGU  &
                        ,RTGV,RTGM,F13U,F13V,F13T,F13M,F23U,F23V,F23T,F23M  &
                        ,DXU,DXV,DXT,DXM,DYU,DYV,DYT,DYM  &
                        ,xm,xt,ym,yt,zm,zt,jdim,ztop,ht,hw)

implicit none

integer :: N1,N2,N3,I,J,K,ITOPO,JDIM
real :: ztop,terdev
real ::   TOPT(N1,N2),TOPU(N1,N2),TOPV(N1,N2),TOPM(N1,N2)  &
         ,RTGT(N1,N2),RTGU(N1,N2),RTGV(N1,N2),RTGM(N1,N2)  &
         ,F13U(N1,N2),F13V(N1,N2),F13T(N1,N2),F13M(N1,N2)  &
         ,F23U(N1,N2),F23V(N1,N2),F23T(N1,N2),F23M(N1,N2)  &
         ,DXU(N1,N2),DXV(N1,N2),DXT(N1,N2),DXM(N1,N2)  &
         ,DYU(N1,N2),DYV(N1,N2),DYT(N1,N2),DYM(N1,N2)  &
         ,xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),ht(*),hw(*)
DATA TERDEV/0./
SAVE

! This routine computes the coordinate transformation constants
! based on the topographical values of TOPT.

DO J=1,N2
   DO I=1,N1-1
      TOPU(I,J)=TOPT(I,J)+(TOPT(I+1,J)-TOPT(I,J))  &
               *(XM(I)-XT(I))/(XT(I+1)-XT(I))
      TERDEV=MAX(TERDEV,ABS(TOPT(I,J)))
   ENDDO
   TOPU(N1,J)=TOPT(N1,J)+(TOPT(N1,J)-TOPT(N1-1,J))  &
             *(XM(N1)-XT(N1))/(XT(N1)-XT(N1-1))
ENDDO
IF(TERDEV.LT.1.E-6)THEN
   ITOPO=0
ELSE
   ITOPO=1
ENDIF

IF(JDIM.EQ.1)THEN
   DO I=1,N1
      DO J=1,N2-1
         TOPV(I,J)=TOPT(I,J)+(TOPT(I,J+1)-TOPT(I,J))  &
                  *(YM(J)-YT(J))/(YT(J+1)-YT(J))
         TOPM(I,J)=TOPU(I,J)+(TOPU(I,J+1)-TOPU(I,J))  &
                  *(YM(J)-YT(J))/(YT(J+1)-YT(J))
      ENDDO
      TOPV(I,N2)=TOPT(I,N2)+(TOPT(I,N2)-TOPT(I,N2-1))  &
                *(YM(N2)-YT(N2))/(YT(N2)-YT(N2-1))
      TOPM(I,N2)=TOPU(I,N2)+(TOPU(I,N2)-TOPU(I,N2-1))  &
                *(YM(N2)-YT(N2))/(YT(N2)-YT(N2-1))
   ENDDO
ELSE
   DO I=1,N1
      DO J=1,N2
         TOPV(I,J)=TOPT(I,J)
         TOPM(I,J)=TOPU(I,J)
      ENDDO
   ENDDO
ENDIF

DO K=1,N3
  HT(K)=ZT(K)/ZTOP-1.0
  HW(K)=ZM(K)/ZTOP-1.0
ENDDO

DO J=1,N2
  DO I=1,N1
    RTGT(I,J)=1.-TOPT(I,J)/ZTOP
    RTGU(I,J)=1.-TOPU(I,J)/ZTOP
    RTGV(I,J)=1.-TOPV(I,J)/ZTOP
    RTGM(I,J)=1.-TOPM(I,J)/ZTOP
    F13T(I,J)=(TOPU(I,J)-TOPU(I-1,J))*DXT(I,J)/RTGT(I,J)
    F13U(I,J)=(TOPT(I+1,J)-TOPT(I,J))*DXU(I,J)/RTGU(I,J)
    F13V(I,J)=(TOPM(I,J)-TOPM(I-1,J))*DXV(I,J)/RTGV(I,J)
    F13M(I,J)=(TOPV(I+1,J)-TOPV(I,J))*DXM(I,J)/RTGM(I,J)
    F23T(I,J)=(TOPV(I,J)-TOPV(I,J-JDIM))*DYT(I,J)/RTGT(I,J)
    F23U(I,J)=(TOPM(I,J)-TOPM(I,J-JDIM))*DYU(I,J)/RTGU(I,J)
    F23V(I,J)=(TOPT(I,J+JDIM)-TOPT(I,J))*DYV(I,J)/RTGV(I,J)
    F23M(I,J)=(TOPU(I,J+JDIM)-TOPU(I,J))*DYM(I,J)/RTGM(I,J)
  ENDDO
ENDDO
DO J=1,N2
  F13T(1,J)=F13U(1,J)
  F13V(1,J)=F13M(1,J)
  F13U(N1,J)=F13T(N1,J)
  F13M(N1,J)=F13V(N1,J)
ENDDO
IF(JDIM.EQ.1)THEN
  DO I=1,N1
    F23T(I,1)=F23V(I,1)
    F23U(I,1)=F23M(I,1)
    F23V(I,N2)=F23T(I,N2)
    F23M(I,N2)=F23U(I,N2)
  ENDDO
ENDIF

return
END SUBROUTINE revu_transfm
