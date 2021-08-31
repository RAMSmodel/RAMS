!##############################################################################
Subroutine fdback (ac,af,dc,df,nzc,nxc,nyc,nzf,nxf,nyf,nf,vnam,sumflg)

use mem_grid

implicit none

integer :: nzc,nxc,nyc,nzf,nxf,nyf,nf,ibeg,jbeg,kbeg,iend,jend,kend  &
   ,nc,iinc,jinc,kv,ifbcf,ic,jc,kc,if,jf,kf
real, dimension(nzc,nxc,nyc) :: ac,dc,sumflg
real, dimension(nzf,nxf,nyf) :: af,df
character(len=*) :: vnam

nc = nxtnest(nf)
CALL azero (nzc*nxc*nyc,sumflg)
ibeg = 2
jbeg = 1 + jdim
kbeg = 2
iend = nxf - 1
jend = nyf - jdim
kend = nzf - 1
iinc = 1
jinc = 1
kv = 0

if (vnam .eq. 'u') then
   ibeg = 1 + nstratx(nf)
   iend = nxf - 1 - nstratx(nf)
   iinc = nstratx(nf)
   ifbcf = 1
elseif (vnam .eq. 'v') then
   jbeg = 1 + nstraty(nf) * jdim
   jend = nyf - (1 + nstraty(nf)) * jdim
   jinc = nstraty(nf)
   ifbcf = 2
elseif (vnam .eq. 'w' .or. vnam .eq. 'terr') then
   if (vnam .eq. 'w') then
      kbeg = 1 + nrz(kpm(2,nf),nf)
      kend = nzf - 1 - nrz(kpm(nzf-1,nf),nf)
      kv = 1
   else
      kbeg = 1
      kend = 1
   endif
   ifbcf = 3
else
   ifbcf = 4
endif

!print*,'fdback:',vnam,':',ibeg,iend,ipm(ibeg,nf),ipm(iend,nf)

kf = kbeg
1 continue
   kc = kpm(kf,nf)
   if (vnam .eq. 'terr') kc = 1

   do jf = jbeg,jend,jinc
      jc = jpm(jf,nf)
      if (vnam .eq. 'p' .or. vnam .eq. 'terr') then
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
                         + af(kf,if,jf) * fbcf(kf,nf,ifbcf)
            sumflg(kc,ic,jc) = 1.
         enddo
      elseif (vnam .eq. 'w') then
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
               + af(kf,if,jf) * fbcf(kf,nf,ifbcf)  &
               * (df(kf,if,jf) + df(kf+kv,if,jf))
            sumflg(kc,ic,jc) = 1.
         enddo
      else
         do if = ibeg,iend,iinc
            ic = ipm(if,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) * sumflg(kc,ic,jc)  &
               + af(kf,if,jf) * fbcf(kf,nf,ifbcf) * df(kf,if,jf)
            sumflg(kc,ic,jc) = 1.
         enddo
      endif
   enddo
   kf = kf + 1
   if (vnam .eq. 'w') kf = kf + nrz(kpm(kf,nf),nf) - 1
if (kf .le. kend) go to 1

do kc = kpm(kbeg,nf),kpm(kend,nf)
   do jc = jpm(jbeg,nf),jpm(jend,nf)

      if (vnam .eq. 'w') then

         do ic = ipm(ibeg,nf),ipm(iend,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) / (dc(kc,ic,jc) + dc(kc+kv,ic,jc))
         enddo

      elseif (vnam .ne. 'p' .and. vnam .ne. 'terr') then

         do ic = ipm(ibeg,nf),ipm(iend,nf)
            ac(kc,ic,jc) = ac(kc,ic,jc) / dc(kc,ic,jc)
         enddo

      endif
   enddo
enddo

return
END SUBROUTINE fdback
