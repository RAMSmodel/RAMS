!##############################################################################
Subroutine fillscr (m1,m2,m3,n1,n2,n3,k1,k2,scr,var)

implicit none

integer :: m1,m2,m3,n1,n2,n3,k1,k2
real :: scr(m1,m2,m3),var(n1,n2,n3)
integer :: i,j,k

do j = 1,n3
   do i = 1,n2
      do k = k1,k2
         scr(k,i,j) = var(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE fillscr

!##############################################################################
Subroutine fillvar (m1,m2,m3,n1,n2,n3,k1,k2,scr,var)

implicit none

integer :: m1,m2,m3,n1,n2,n3,k1,k2
real :: scr(m1,m2,m3),var(n1,n2,n3)
integer :: i,j,k

do j = 1,n3
   do i = 1,n2
      do k = k1,k2
         var(k,i,j) = scr(k,i,j)
      enddo
   enddo
enddo

return
END SUBROUTINE fillvar

!##############################################################################
Subroutine dnswt2 (nzmx,nxmx,nymx,n1,n2,n3,var1,dn0x,vnam,idir)

implicit none

integer :: nzmx,nxmx,nymx,n1,n2,n3,idir
real :: var1(nzmx,nxmx,nymx),dn0x(n1,n2,n3)
character(len=*) :: vnam
integer :: i,j,k

if (idir .eq. 1) then

   if (vnam .eq. 'w') then
      do j = 1,n3
         do i = 1,n2
            do k = 1,n1-1
               var1(k,i,j) = var1(k,i,j)  &
                  * .5 * (dn0x(k,i,j) + dn0x(k+1,i,j))
            enddo
         enddo
      enddo
   else
      do j = 1,n3
         do i = 1,n2
            do k = 1,n1
               var1(k,i,j) = var1(k,i,j) * dn0x(k,i,j)
            enddo
         enddo
      enddo
   endif

elseif (idir .eq. 2) then

   if (vnam .eq. 'w') then
      do j = 1,n3
         do i = 1,n2
            do k = 1,n1-1
               var1(k,i,j) = var1(k,i,j)  &
                  / (.5 * (dn0x(k,i,j) + dn0x(k+1,i,j)))
            enddo
         enddo
      enddo
   else
      do j = 1,n3
         do i = 1,n2
            do k = 1,n1
               var1(k,i,j) = var1(k,i,j) / dn0x(k,i,j)
            enddo
         enddo
      enddo
   endif

endif

return
END SUBROUTINE dnswt2
