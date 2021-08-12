!##############################################################################
Subroutine deblank (str1,str2,nch)

implicit none

character(len=*) :: str1,str2
integer :: n,ln,nch

! strips blanks from a string and returns number of chars

str2=' '
ln=len(str1)
nch=0
do n=1,ln
   if(str1(n:n).ne.' ') then
      nch=nch+1
      str2(nch:nch)=str1(n:n)
   endif
enddo

return
END SUBROUTINE deblank

!##############################################################################
integer Function lastslash (str)

implicit none

character(len=*) :: str
integer :: n,ln

! returns last slash character position from a string

ln=len(str)
do n=ln,1,-1
   if(str(n:n).eq.'/') then
      lastslash = n
      return
   endif
enddo
lastslash = 0

return
END FUNCTION lastslash

!##############################################################################
Subroutine char_strip_var (line,var,line2)

implicit none

character(len=*) :: line,var,line2
integer :: nn,ncl,nb

! removes instances of a substring from a string

nb=0
ncl=len(line)
do nn=1,ncl
   if(line(nn:nn).ne.' ') then
      nb=index(line(nn:),' ')
      var=line(nn:nn+nb-1)
      goto 25
   endif
enddo
25 continue
line2=line(nn+nb-1:)

return
END SUBROUTINE char_strip_var

!##############################################################################
Subroutine parse (str,tokens,ntok)

implicit none

integer :: ntok
character(len=*) :: str,tokens(*)
character(len=1) :: sep
integer, parameter :: ntokmax=100
integer :: n,nc,npt,nch,ntbeg,ntend

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

sep=' '
ntok=0
npt=1
nch=len_trim(str)
nc=1
ntbeg=0
ntend=0
do ntok=1,ntokmax
   do n=nc,nch
      if(str(n:n) /= sep) then
         ntbeg=n
         exit
      endif
   enddo
   
   do n=ntbeg,nch
      if(str(n:n) == sep .or. str(n:n) == char(10)) then ! look for \n also
         ntend=n-1
         goto 22
      endif
      if(n == nch) then
         ntend=n
         exit
      endif
   enddo
   22 continue
   if(ntbeg==0.or.ntend==0) stop 'Error in PARSE routine'
   tokens(ntok)=str(ntbeg:ntend)
   nc=ntend+1
   if(nc >= nch) goto 25
enddo

25 continue

return
END SUBROUTINE parse

!##############################################################################
Subroutine tokenize1 (str1,tokens,ntok,toksep)

use grid_dims

implicit none

integer :: ntok
character(len=*) :: str1,tokens(*)
character(len=1) :: toksep
character(len=strl1) :: str
integer :: nch,ist,npt,nc

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

CALL deblank (str1,str,nch)

ist=1
if(str(1:1) == toksep) ist=2
npt=ist
ntok=0
do nc=ist,nch
   if(str(nc:nc) == toksep .or. nc == nch) then
      if(nc-npt >= 0) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc-1)
         if(nc == nch .and. str(nc:nc) /= toksep) then
            tokens(ntok)=str(npt:nc)
            !print*,'ttttttttttt2:',ntok,npt,nc,tokens(ntok)
            exit
         endif
         npt=nc+1
      elseif(nc == nch) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc)
         exit
      endif
   endif
enddo

return
END SUBROUTINE tokenize1

!##############################################################################
Subroutine tokenize2 (str,tokens,ntok,toksep)

implicit none

integer :: ntok
character(len=*) :: str,tokens(*)
character(len=1) :: toksep
integer :: nch,ist,npt,nc

! same as tokenize1, but doesn't deblank the string

str = adjustl(str)
nch = len_trim(str)

ist=1
if(str(1:1) == toksep) ist=2
npt=ist
ntok=0
do nc=ist,nch
   if(str(nc:nc) == toksep .or. nc == nch) then
      if(nc-npt >= 0) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc-1)
         if(nc == nch .and. str(nc:nc) /= toksep) then
            tokens(ntok)=str(npt:nc)
            exit
         endif
         npt=nc+1
      elseif(nc == nch) then
         ntok=ntok+1
         tokens(ntok)=str(npt:nc)
         exit
      endif
   endif
enddo

return
END SUBROUTINE tokenize2

!##############################################################################
integer Function letter (str)

implicit none

character(len=*) :: str

! First character alpha check - test to see if the first character of
! the string STR is alphabetic: LETTER equals 0 if 'no', = 1 if 'yes'.

letter = 0
if((str(1:1).ge.'A'.and.str(1:1).le.'Z').or.  &
   (str(1:1).ge.'a'.and.str(1:1).le.'z')) letter = 1
   
return
END FUNCTION letter

!##############################################################################
integer Function numberchk (str)

implicit none

character(len=*) :: str

! First character number check - test to see if the first character of
! the string STR is numeric:  NUMBER equals 0 if 'no', = 1 if 'yes' (includ
! a decimal point or minus sign).

numberchk = 0
if(str(1:1).ge.'0'.and.str(1:1).le.'9') numberchk = 1
if(str(1:1).eq.'.'.or.str(1:1).eq.'-') numberchk = 1

return
END FUNCTION numberchk

!##############################################################################
integer Function letint (str)

implicit none

character(len=*) :: str

! First character integer variable check - test to see if the first
! character of STR is an I, J, K, L, M, or N, or i, j, k, l ,m or n:
! LETINT equals 0 if 'no', = 1 if 'yes'

letint = 0
if((str(1:1).ge.'I'.and.str(1:1).le.'N').or.  &
   (str(1:1).ge.'i'.and.str(1:1).le.'n')) letint = 1
   
return
END FUNCTION letint

!##############################################################################
integer Function letquo (str)

implicit none

character(len=*) :: str

! First character quote check - test to see if the first character
! of STR is a quote:  LETQUO equals 0 if 'no', = 1 if 'yes'.

letquo = 0
if(str(1:1).eq.'''') letquo = 1

return
END FUNCTION letquo
