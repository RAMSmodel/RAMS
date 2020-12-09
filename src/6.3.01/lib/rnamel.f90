!##############################################################################
!------------------------------------------------------------------
!
!      Namelist reading routines
!
!------------------------------------------------------------------

!##############################################################################
Subroutine varsetf (varn,var,is1,maxsub,fnum,fmin,fmax)

implicit none

integer :: is1,maxsub
character(len=*) :: varn
real :: var,fnum,fmin,fmax

! Assign real value read in (FNUM) to the corresponding real variable
! in the NAMELIST and do a bounds check

if(is1.le.maxsub) then
  var=fnum
else
  print 9,trim(varn),is1,maxsub,fnum
  9 format(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,F20.6)
endif

if(fnum.lt.fmin.or.fnum.gt.fmax) then
  print 10,trim(varn),fnum,fmin,fmax
  10 format(' -- ERROR --   Input variable - ',A,' - set to ',F18.5  &
      ,/,'                 allowable range ',F15.5,' to ',F15.5)
  stop
endif

return
END SUBROUTINE varsetf

!##############################################################################
Subroutine varseti (varn,ivar,is1,maxsub,inum,imin,imax)

implicit none

integer :: ivar,is1,maxsub,inum,imin,imax
character(len=*) :: varn

! Assign integer value read in (INUM) to the corresponding integer var
! (IVAR) in the NAMELIST and do a bounds check

if(is1.le.maxsub) then
  ivar=inum
else
  print 9,trim(varn),is1,maxsub,inum
  9 format(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,I20)
endif

if(inum.lt.imin.or.inum.gt.imax) then
  print 10,trim(varn),inum,imin,imax
  10 format(' -- ERROR --   Input variable - ',A,' - set to ',I10  &
      ,/,'                 allowable range ',I10,' to ',I10)
  stop
endif

return
END SUBROUTINE varseti

!##############################################################################
Subroutine varsetc (varn,var,is1,maxsub,ch,imin,imax)

implicit none

integer :: is1,maxsub,imin,imax
character(len=*) :: varn,ch,var
integer :: lch

! Assign character value read in (CH) to the corresponding character
! variable (VAR) in the NAMELIST and do a bounds check

if(is1.le.maxsub) then
  var=ch
else
  print 9,trim(varn),is1,maxsub,trim(ch)
  9 FORMAT(' -- ERROR --   Input variable - ',A,' - attempting'  &
        ,' to read extra values, values ignored.',/,  &
         ' Subscript, max dimension, value ',2I6,A)
endif

lch=len(ch)

if(lch.lt.imin.or.lch.gt.imax) then
  print 10,trim(varn),trim(ch),imin,imax
  10 FORMAT(' -- ERROR --   Input variable - ',A,' - set to ',A  &
      ,/,'                 allowable length ',I10,' to ',I10)
endif

return
END SUBROUTINE varsetc

!##############################################################################
Subroutine varchk (varn,group,namelst,iname,nname,inrflg)

implicit none

integer :: nname,inrflg,icompatible
integer :: iname(nname)
character(len=*) :: varn,group,namelst(nname)
integer :: nn,inoset

! This routine checks that all member variables of the NAMELIST
! specified by GROUP have been assigned values.
!
! GROUP   - name of the NAMELIST being checked
! INNAME  - storage vector for counting number of times each member
!           of the NAMELIST has been assigned a value
! NAMELST - character vector containing names of members of NAMELIST
! NNAME   - number of variables in NAMELIST
! VARN    - name of NAMELIST variable being checked

INRFLG=0
ICOMPATIBLE=0

IF(VARN.NE.'$END')THEN

  ! Update INAME value corresponding to the NAMELIST variable which has
  ! been assigned a value in the call to NVFILL.

  DO NN=1,NNAME
    IF(VARN.EQ.NAMELST(NN)) THEN
      INAME(NN)=1
      RETURN
    ENDIF
  enddo

  INRFLG=1
  PRINT 20,GROUP(1:len_trim(GROUP)),VARN(1:len_trim(VARN))
  20 FORMAT(' Extra variable in namelist - ',A,'- variable name: ',A)
  ICOMPATIBLE=1

ELSE

  ! End of NAMELIST input has been reached - check that all variables have
  ! been given values

  INRFLG=1
  INOSET=0
  DO NN=1,NNAME
    IF(INAME(NN).EQ.0) THEN
      IF(INOSET.EQ.0) THEN
        PRINT 30,GROUP(1:len_trim(GROUP))
      ENDIF
      INOSET=1
      PRINT 31,NAMELST(NN)(1:len_trim(NAMELST(NN)))
      ICOMPATIBLE=1
    ENDIF
  enddo

  30 FORMAT(' Variables not set in the -- ',A,' -- namelist:')
  31 FORMAT(' Missing variable name: ',A)

ENDIF

if(ICOMPATIBLE == 1)then
 print*,''
 print*,'Your RAMSIN namelist file has extra flags or missing flags and is not' 
 print*,' compatible with this version of RAMS. Please update your RAMSIN.'
 stop
endif

return
END SUBROUTINE varchk

!##############################################################################
Subroutine ch2int (str,int)

implicit none

integer :: int
character(len=*) :: str
character(len=8) :: form
integer :: nc

! Read integer value INT from character string STR
nc=len_trim(str)
write(form,90)nc
90 format('(i',i2,')')
read(str,form)int
  
return
END SUBROUTINE ch2int

!##############################################################################
Subroutine ch2real (str,fnum)

implicit none

real :: fnum
character(len=*) :: str
character(len=8) :: form
integer :: nc

! Read real value FNUM from character string STR

nc=len_trim(str)
write(form,90)nc
90 format('(f',i2,'.0)')
read(str,form)fnum
  
return
END SUBROUTINE ch2real

!##############################################################################
Subroutine ch2ch (str,chval,ncw)

implicit none

character(len=*) :: str,chval

integer :: nc,ncw,ncstr

! Remove trailing blanks from character string STR and store remaining
! NCW characters in character string CHVAL

ncstr=len(str)
do 10 nc=ncstr,1,-1
  if(str(nc:nc).eq.' ') goto 10
  chval=str(2:nc)
  ncw=nc-2
  return
10 continue

return
END SUBROUTINE ch2ch

!##############################################################################
Subroutine findgr (iunit,group)

use grid_dims

implicit none

integer :: iunit
character(len=*) :: group
character(len=strl1) :: line
integer :: nr,ind

! This routine checks to see if any of the first MAXREC lines on input
! unit IUNIT contains the character string GROUP

do nr=1,maxrec
  read(iunit,'(a128)',end=100) line
  ind=index(line,group)
  if(ind.ne.0) return
enddo
100 continue
print *,' Namelist read error -- group not found -- ',group

return
END SUBROUTINE findgr

!##############################################################################
Subroutine strip (lin1,lin2,nc2)

implicit none

integer :: nc2
character(len=*) :: lin1,lin2
integer :: iquote,nc,nl

! This routine strips blank characters from character string LIN1
! (as well as comments beginning with an '!') and stores the stripped-
! down remainder in character string LIN1.  NC2 is the number of
! characters in LIN2.

nl=len(lin1)

nc2=0
iquote=0
do nc=1,nl
  if(iquote.eq.0) then
    if(lin1(nc:nc).ne.' ') then
      if(lin1(nc:nc).eq.'!') return
      nc2=nc2+1
      lin2(nc2:nc2)=lin1(nc:nc)
    endif
  else
    nc2=nc2+1
    lin2(nc2:nc2)=lin1(nc:nc)
  endif
  if(lin1(nc:nc).eq.'''') then
    if(iquote.eq.0) then
      iquote=1
    else
      iquote=0
    endif
  endif
enddo

return
END SUBROUTINE strip

!##############################################################################
Subroutine toknze (str,nch,tokens,ntok)

implicit none

integer :: nch,ntok
character(len=*) :: str,tokens(*)
integer, parameter :: nsep=4
character(len=1) :: toksep(nsep)
data toksep/'=',',','(',')'/

integer :: nc,ns,npt

! This routine "parses" character string STR into different pieces
! or tokens by looking for one of four possible token separators (TOKS
! STR contains NCH characters.  The number of tokens identified is NTO
! the character string tokens are stored in TOKENS.


ntok=0
npt=1
do 10 nc=1,nch
  do 5 ns=1,nsep
    if(str(nc:nc).eq.toksep(ns))then
      if(nc-npt.ge.1)then
        ntok=ntok+1
        tokens(ntok)=str(npt:nc-1)
      endif
      ntok=ntok+1
      tokens(ntok)=str(nc:nc)
      npt=nc+1
      goto 10
    endif
  5 continue
10 continue

return
END SUBROUTINE toknze
