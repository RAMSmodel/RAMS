!##############################################################################
Subroutine gridset ()

use mem_grid
use rconstants

implicit none

integer :: ifm,icm,k,nestza,iinc,icnt,if1,jinc  &
   ,jcnt,jf,kinc,kcnt,kf,nrat,i,j,kcy,kcw,kk
real :: centx1,centy1,centx,centy,dzr,dsum,dzrcm,dzrfm,tsum
real :: zmnvc(-1:nzpmax+1,maxgrds)

!Grid set up at model initialization

! Fill NRZ with the variable nest ratios.

   do ifm = 1,ngrids
      do k = 1,nzpmax
         nrz(k,ifm) = 1
      enddo
   enddo

   nestza = abs(nestz)

   if (nestza > 1 .and. nestza <= ngrids) then
      do k = 2,nnzp(1)
         nrz(k,nestza) = max(1,nstratz(k))
      enddo
      nrz(1,nestza) = nstratz(2)
   endif

!     FILL ALL DELTAXN VALUES, AND FIND NINEST AND NJNEST VALUES
!     IF SET TO ZERO IN NAMELIST

   deltaxn(1) = deltax
   CALL ll_xy (centlat(1),centlon(1),polelat,polelon,centx1,centy1)
   xmn(1,1) = centx1 - 0.5 * float(nnxp(1)-2) * deltaxn(1)
   ymn(1,1) = centy1 - 0.5 * float(nnyp(1)-2) * deltaxn(1)
   if (nnyp(1) == 1) ymn(1,1) = centy1

   do ifm = 1,ngrids
      icm = nxtnest(ifm)
      if (icm >= 1) then
         deltaxn(ifm) = deltaxn(icm) / float(nstratx(ifm))
         if (ninest(ifm) <= 0 .or. njnest(ifm) <= 0)  &
            CALL ll_xy (centlat(ifm),centlon(ifm),polelat  &
               ,polelon,centx,centy)
         if (ninest(ifm) <= 0) then
            xmn(1,ifm) = centx - 0.5 * float(nnxp(ifm)-2) * deltaxn(ifm)
            ninest(ifm) = int((xmn(1,ifm) - xmn(1,icm)) / deltaxn(icm) + 1.5)
         endif
         xmn(1,ifm) = xmn(1,icm) + deltaxn(icm) * float(ninest(ifm)-1)
         if (njnest(ifm) <= 0) then
            ymn(1,ifm) = centy - 0.5 * float(nnyp(ifm)-2) * deltaxn(ifm)
            njnest(ifm) = int((ymn(1,ifm) - ymn(1,icm)) / deltaxn(icm) + 1.5)
         endif
         ymn(1,ifm) = ymn(1,icm) + deltaxn(icm) * float(njnest(ifm)-1)
      endif
   enddo

!     Fill IPM, JPM and KPM arrays with parent grid index values for
!     all fine grids.

do ifm = 2,ngrids
   icm = nxtnest(ifm)
   if (icm >= 1) then
      ipm(1,ifm) = ninest(ifm)
      iinc = 1
      icnt = 0
      do if1 = 2,nnxp(ifm)
         ipm(if1,ifm) = ipm(if1-1,ifm) + iinc
         icnt = icnt + 1
         if (icnt >= nstratx(ifm)) then
            icnt = 0
            iinc = 1
         else
            iinc = 0
         endif
      enddo

      jpm(1,ifm) = njnest(ifm)
      jinc = 1
      jcnt = 0
      do jf = 2,nnyp(ifm)
         jpm(jf,ifm) = jpm(jf-1,ifm) + jinc
         jcnt = jcnt + 1
         if (jcnt >= nstraty(ifm)) then
            jcnt = 0
            jinc = 1
         else
            jinc = 0
         endif
      enddo

      if (nknest(ifm) == 1) then
         kpm(1,ifm) = 2
         kinc = 0
      else
         kpm(1,ifm) = nknest(ifm)
         kinc = 1
      endif
      kcnt = 0
      do kf = 2,nnzp(ifm)
         kpm(kf,ifm) = kpm(kf-1,ifm) + kinc
         nrat = nrz(kpm(kf,ifm),ifm)
         kcnt = kcnt + 1
         if (kcnt >= nrat .and. (kf < nnzp(ifm)-1 .or.  &
            kpm(kf,ifm) < nnzp(icm)-1)) then
            kcnt = 0
            kinc = 1
         else
            kinc = 0
         endif
      enddo
   endif
enddo

! Calculate xmn and ymn for the coarse grid(s) for an initial start.

   do i = 2,nnxp(1)
      xmn(i,1) = xmn(i-1,1) + deltaxn(1)
   enddo

   do j = 2,nnyp(1)
      ymn(j,1) = ymn(j-1,1) + deltaxn(1)
   enddo

! compute xmn and ymn for any required nested grids, and xtn and ytn for
! any required grids.

do ifm = 1,ngrids
   icm = nxtnest(ifm)
   if (icm >= 1) then
      xmn(1,ifm) = xmn(ninest(ifm),icm)
      ymn(1,ifm) = ymn(njnest(ifm),icm)
   endif
   do i = 2,nnxp(ifm)
      xmn(i,ifm) = xmn(i-1,ifm) + deltaxn(ifm)
      xtn(i,ifm) = .5 * (xmn(i,ifm) + xmn(i-1,ifm))
   enddo
   xtn(1,ifm) = 1.5 * xmn(1,ifm) - .5 * xmn(2,ifm)

   if (jdim == 1) then
      do j = 2,nnyp(ifm)
         ymn(j,ifm) = ymn(j-1,ifm) + deltaxn(ifm)
         ytn(j,ifm) = .5 * (ymn(j,ifm) + ymn(j-1,ifm))
      enddo
      ytn(1,ifm) = 1.5 * ymn(1,ifm) - .5 * ymn(2,ifm)
   else
      ytn(1,ifm) = ymn(1,ifm)
   endif
enddo

! calculate zmn for the coarse grid(s).

if (deltaz == 0.) then
   do k = 1,nnzp(1)
      zmn(k,1) = zz(k)
   enddo
else
   zmn(1,1) = 0.
   zmn(2,1) = deltaz
   do k = 3,nnzp(1)
   if (zmn(k-1,1)<zdelay) then
      dzr = 1.
   else
      dzr = dzrat
   endif
   if (nestz < -1 ) then
      if( k > kpm(2,nestza) .and.  &
          k <= kpm(nnzp(nestza)-1,nestza) .and.  &
          nrz(k,nestza) /= nrz(k-1,nestza)) then
         if (max(nrz(k,nestza),nrz(k-1,nestza)) == 2)  &
            dzr = 1.325 ** (nrz(k,nestza) - nrz(k-1,nestza))
         if (max(nrz(k,nestza),nrz(k-1,nestza)) == 3)  &
            dzr = 1.124 ** (nrz(k,nestza) - nrz(k-1,nestza))
         if (max(nrz(k,nestza),nrz(k-1,nestza)) == 4)  &
            dzr = 1.066 ** (nrz(k,nestza) - nrz(k-1,nestza))
         if (max(nrz(k,nestza),nrz(k-1,nestza)) == 5)  &
            dzr = 1.042 ** (nrz(k,nestza) - nrz(k-1,nestza))
      endif
   endif
   zmn(k,1) = zmn(k-1,1) + min(dzr * (zmn(k-1,1) - zmn(k-2,1)),dzmax)
   enddo
   deltazn(1) = deltaz
endif

! fill scratch array znmvc with zmn for coarse grid(s)

ztop = zmn(nnzp(1)-1,1)

do k = 1,nnzp(1)
   zmnvc(k,1) = zmn(k,1)
enddo
zmnvc(0,1) = -(zmnvc(2,1) - zmnvc(1,1)) ** 2 / (zmnvc(3,1) - zmnvc(2,1))
zmnvc(-1,1) = zmnvc(0,1) - (zmnvc(1,1) - zmnvc(0,1)) ** 2  &
   / (zmnvc(2,1) - zmnvc(1,1))
zmnvc(nnzp(1)+1,1) = zmnvc(nnzp(1),1)  &
   + (zmnvc(nnzp(1),1) - zmnvc(nnzp(1)-1,1)) ** 2  &
   / (zmnvc(nnzp(1)-1,1) - zmnvc(nnzp(1)-2,1))

! get zmn coordinates for all nested grids.

do ifm = 2,ngrids
   icm = nxtnest(ifm)
   if (icm >= 1) then
      kcy = max(-1,nrz(kpm(2,ifm),ifm)-3)
      kcw = 0 !Variable initialized
      do k = -1,nnzp(ifm) + 1
         if (k <= 0) then
            nrat = nrz(kpm(2,ifm),ifm)
            kcw = nknest(ifm)-1
            if (k == -1 .and. nrat == 1) kcw = nknest(ifm) - 2
         elseif (k == nnzp(ifm)+1) then
            nrat = nrz(kpm(nnzp(ifm)-1,ifm),ifm)
         else
            nrat = nrz(kpm(k,ifm),ifm)
         endif
         kcy = mod(kcy+1,nrat)

         if (kcy == 0) then
            if (k >= 1) kcw = kcw + 1
            zmnvc(k,ifm) = zmnvc(kcw,icm)
         else
            tsum = 0.  !Variable initialized
            dzrfm = 0. !Variable initialized
            dsum = 0.  !Variable initialized
            if (kcy == 1 .or. k == -1) then
               dsum = 0.
               dzrcm = sqrt((zmnvc(kcw+2,icm) - zmnvc(kcw+1,icm))  &
                  * nrz(max(1,kcw),ifm)  &
                  / ((zmnvc(kcw,icm) - zmnvc(kcw-1,icm)) * nrz(kcw+2,ifm)))
               dzrfm = dzrcm ** (1. / float(nrat))
               tsum = 0.
               do kk = 1,nrat
                  tsum = tsum + dzrfm ** (kk-1)
               enddo
            endif
            if (k == -1) then
               do kk = 1,kcy-1
                  dsum = dsum + ((zmnvc(kcw+1,icm)-zmnvc(kcw,icm))  &
                     / tsum) * dzrfm ** (kk-1)
               enddo
            endif
            dsum = dsum + (zmnvc(kcw+1,icm) - zmnvc(kcw,icm))  &
               / tsum * dzrfm ** (kcy-1)
            zmnvc(k,ifm) = zmnvc(kcw,icm) + dsum
         endif
         if (k >= 1 .and. k <= nnzp(ifm)) zmn(k,ifm) = zmnvc(k,ifm)
      enddo
   endif
enddo

!     compute ztn values for all grids by geometric interpolation.

do ifm = 1,ngrids
   do k = 1,nnzp(ifm)
      dzrfm = sqrt(sqrt((zmnvc(k+1,ifm) - zmnvc(k,ifm))  &
         / (zmnvc(k-1,ifm) - zmnvc(k-2,ifm))))
      ztn(k,ifm) = zmnvc(k-1,ifm)  &
         + (zmnvc(k,ifm) - zmnvc(k-1,ifm)) / (1. + dzrfm)
   enddo

!     compute other arrays based on the vertical grid.

   do k = 1,nnzp(ifm)-1
      dzmn(k,ifm) = 1. / (ztn(k+1,ifm) - ztn(k,ifm))
   enddo
   do k = 2,nnzp(ifm)
      dztn(k,ifm) = 1. / (zmn(k,ifm) - zmn(k-1,ifm))
   enddo
   do k = 2,nnzp(ifm)-1
      dzm2n(k,ifm) = 1. / (zmn(k+1,ifm) - zmn(k-1,ifm))
      dzt2n(k,ifm) = 1. / (ztn(k+1,ifm) - ztn(k-1,ifm))
   enddo

   dzmn(nnzp(ifm),ifm) = dzmn(nnzp(ifm)-1,ifm)  &
      * dzmn(nnzp(ifm)-1,ifm) / dzmn(nnzp(ifm)-2,ifm)
   dztn(1,ifm) = dztn(2,ifm) * dztn(2,ifm) / dztn(3,ifm)

   dzm2n(1,ifm) = dzm2n(2,ifm)
   dzm2n(nnzp(ifm),ifm) = dzm2n(nnzp(ifm)-1,ifm)
   dzt2n(1,ifm) = dzt2n(2,ifm)
   dzt2n(nnzp(ifm),ifm) = dzt2n(nnzp(ifm)-1,ifm)

   deltazn(ifm) = zmn(2,ifm) - zmn(1,ifm)

enddo

CALL cofnest ()

return
END SUBROUTINE gridset
