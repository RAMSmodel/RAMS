!##############################################################################
Subroutine nvfills (group,v,iv,fv,cv,nva)

! Namelist reading routines for ISAN

! Keep the list of variables here in sync with the list in
! broadcast_config() in mpass_init.f90

! This routine places values read in for a particular NAMELIST
! variable in the storage location(s) for that variable
!
! GROUP   - NAMELIST group name
! VR      - NAMELIST variable name
! II      - vector of integer values
! FF      - vector of real values
! CC      - vector of character strings
! NV      - number of values to be inserted

use isan_coms

implicit none

character(len=*) :: group,v,cv
integer :: iv,nva
real :: fv

integer, parameter :: nvcont=14,nvinst3=23
integer :: icontrol(nvcont),iinst3(nvinst3)
character(len=16) :: control(nvcont),inst3(nvinst3)
data icontrol/nvcont*0/,iinst3/nvinst3*0/
data control/'IASRFCE'  &
     ,'ISAN_INC','I1ST_FLG','IAPR','IARAWI','VAR_HFILE'  &
     ,'VARPFX','IOFLGISZ','IOFLGVAR'  &
     ,'IUPA_FLG','ISFC_FLG','ISZSTAGE','IVRSTAGE','IDATAIN'/
data inst3/'NOTSTA','NOTID','IOBSWIN','WVLNTH','RESPON'  &
     ,'SWVLNTH','STASEP','IGRIDFL','GRIDWT'  &
     ,'GOBSEP','GOBRAD','MAXSTA','MAXSFC'  &
     ,'NISN','LEVTH','NIGRIDS','SIGZWT','TOPSIGZ'  &
     ,'HYBBOT','HYBTOP','SFCINF','NFEEDVAR','USED_FILE'/

integer :: inrflg

IF(GROUP.EQ.'$ISAN_CONTROL') THEN
   CALL varchk (V,GROUP,CONTROL,ICONTROL,NVCONT,INRFLG)
   IF(INRFLG.EQ.1) RETURN
   IF(V.EQ.'ISZSTAGE')  CALL varseti (V,ISZSTAGE,NVA,1,IV,0,1)
   IF(V.EQ.'IVRSTAGE')  CALL varseti (V,IVRSTAGE,NVA,1,IV,0,1)
   IF(V.EQ.'ISAN_INC')  CALL varseti (V,ISAN_INC,NVA,1,IV,0,9999)
   IF(V.EQ.'I1ST_FLG')  CALL varseti (V,I1ST_FLG,NVA,1,IV,1,2)
   IF(V.EQ.'IUPA_FLG')  CALL varseti (V,IUPA_FLG,NVA,1,IV,1,10)
   IF(V.EQ.'ISFC_FLG')  CALL varseti (V,ISFC_FLG,NVA,1,IV,1,10)
   IF(V.EQ.'IAPR')      CALL varsetc (V,IAPR,NVA,1,CV,0,80)
   IF(V.EQ.'IARAWI')    CALL varsetc (V,IARAWI,NVA,1,CV,0,80)
   IF(V.EQ.'IASRFCE')   CALL varsetc (V,IASRFCE,NVA,1,CV,0,80)
   IF(V.EQ.'VAR_HFILE') CALL varsetc (V,VAR_HFILE,NVA,1,CV,0,strl1)
   IF(V.EQ.'VARPFX')    CALL varsetc (V,VARPFX,NVA,1,CV,0,80)
   IF(V.EQ.'IOFLGISZ')  CALL varseti (V,IOFLGISZ,NVA,1,IV,0,1)
   IF(V.EQ.'IOFLGVAR')  CALL varseti (V,IOFLGVAR,NVA,1,IV,0,1)
   IF(V.EQ.'IDATAIN')   CALL varseti (V,IDATAIN,NVA,1,IV,0,1)
ENDIF

IF(GROUP.EQ.'$ISAN_ISENTROPIC') THEN
   CALL varchk (V,GROUP,INST3,IINST3,NVINST3,INRFLG)
   IF(INRFLG.EQ.1) RETURN
   IF(V.EQ.'NIGRIDS')  CALL varseti (V,NIGRIDS,NVA,1,IV,1,9)
   IF(V.EQ.'NFEEDVAR') CALL varseti (V,NFEEDVAR,NVA,1,IV,0,1)
   IF(V.EQ.'SIGZWT')   CALL varsetf (V,SIGZWT,NVA,1,FV,0.,1.)
   IF(V.EQ.'NISN')     CALL varseti (V,NISN,NVA,1,IV,0,MAXISN)
   IF(V.EQ.'LEVTH')    CALL varseti (V,LEVTH(NVA),NVA,MAXISN,IV,200,800)
   IF(V.EQ.'TOPSIGZ')  CALL varsetf (V,TOPSIGZ,NVA,1,FV,0.,100000.)
   IF(V.EQ.'HYBBOT')   CALL varsetf (V,HYBBOT,NVA,1,FV,0.,100000.)
   IF(V.EQ.'HYBTOP')   CALL varsetf (V,HYBTOP,NVA,1,FV,0.,100000.)
   IF(V.EQ.'SFCINF')   CALL varsetf (V,SFCINF,NVA,1,FV,0.,100000.)
   IF(V.EQ.'MAXSTA')   CALL varseti (V,MAXSTA,NVA,1,IV,0,100000)
   IF(V.EQ.'MAXSFC')   CALL varseti (V,MAXSFC,NVA,1,IV,0,100000)
   IF(V.EQ.'GOBSEP')   CALL varsetf (V,GOBSEP,NVA,1,FV,0.,100.)
   IF(V.EQ.'GOBRAD')   CALL varsetf (V,GOBRAD,NVA,1,FV,0.,100.)
   IF(V.EQ.'WVLNTH')   CALL varsetf (V,WVLNTH(NVA),NVA,maxagrds,FV,0.,1e5)
   IF(V.EQ.'SWVLNTH')  CALL varsetf (V,SWVLNTH(NVA),NVA,maxagrds,FV,0.,1E5)
   IF(V.EQ.'RESPON')   CALL varsetf (V,RESPON(NVA),NVA,maxagrds,FV,0.,1.)
   IF(V.EQ.'STASEP')   CALL varsetf (V,STASEP,NVA,1,FV,0.,100.)
   IF(V.EQ.'IGRIDFL')  CALL varseti (V,IGRIDFL,NVA,1,IV,0,4)
   IF(V.EQ.'GRIDWT')   CALL varsetf (V,GRIDWT(NVA),NVA,maxagrds,FV,0.,1.)
   IF(V.EQ.'NOTSTA')   CALL varseti (V,NOTSTA,NVA,1,IV,0,50)
   IF(V.EQ.'NOTID')    CALL varsetc (V,NOTID(NVA),NVA,50,CV,1,8)
   IF(V.EQ.'USED_FILE') CALL varsetc (V,USED_FILE,NVA,1,CV,0,strl1)
   IF(V.EQ.'IOBSWIN')  CALL varseti (V,IOBSWIN,NVA,1,IV,-21600,21600)
ENDIF

return
END SUBROUTINE nvfills

!##############################################################################
Subroutine namein_isan (IUNIT,GROUP)

! This routine is called by routines ISAN_DRIVER
! to input the values for the NAMELIST specified by the character string
! GROUP from input unit IUNIT.

use grid_dims

implicit none

character(len=*) :: group
integer :: iunit
integer, parameter :: linvar=10
character(len=strl1) :: varn,line,linew,value(maxvalues),tokens(50)
integer :: nvalue,nvarn,nr,ncw,nt,ntok,int
real :: fnum
integer, external :: letter
integer, external :: numberchk
integer, external :: letquo

rewind iunit

NVALUE=0
NVARN=0
CALL findgr (IUNIT,GROUP)
DO 10 NR=1,MAXREC
  READ(IUNIT,'(A80)',END=100)LINE
  CALL strip (LINE,LINEW,NCW)
  NCW=MAX(NCW,1)
  CALL toknze (LINEW,NCW,TOKENS,NTOK)
  IF(LINEW(1:NCW).EQ.'$END') GOTO 100

  NT=1
  20 CONTINUE
    IF(NT.GT.NTOK) GOTO 10
    IF(letter(TOKENS(NT)).EQ.1) THEN
      IF(NVARN.GT.0) CALL nvtran_isan (GROUP,VARN,VALUE,NVALUE)
      NVALUE=0
      VARN=TOKENS(NT)
      NVARN=NVARN+1
      NT=NT+1
      IF(TOKENS(NT).EQ.'(' .or. TOKENS(NT).EQ.',') THEN
        stop 'Error in RAMSIN syntax. Routine namein_isan'
      ENDIF
    ELSEIF(numberchk(TOKENS(NT)).EQ.1 .OR. letquo(TOKENS(NT)).EQ.1) THEN
      NVALUE=NVALUE+1
      VALUE(NVALUE)=TOKENS(NT)
    ENDIF
    NT=NT+1
  GOTO 20
10 CONTINUE

100 CONTINUE
CALL nvtran_isan (GROUP,VARN,VALUE,NVALUE)
VARN='$END'
CALL nvfills (GROUP,VARN,INT,FNUM,LINEW(1:NCW),NR)

return
END SUBROUTINE namein_isan

!##############################################################################
Subroutine nvtran_isan (GROUP,VARN,VALUE,NVALUE)

! This routine converts a variable value(s) (VALUE) read in as a
! character string(s) to the proper type (integer, real, or character)
! and assigns it (them) to the storage location(s) corresponding to
! the NAMELIST variable VARN.  GROUP is the NAMELIST
! group name.  NVALUE is the number of values to be assigned.

use grid_dims

implicit none

integer :: NVALUE
character(len=*) :: GROUP,VARN,VALUE(*)
character(len=strl1) :: CHVAL

integer :: nv,int,ncw
real :: fnum
integer, external :: letint
integer, external :: letquo

ncw=1

IF(letint(VARN).EQ.1) THEN
  DO NV=1,NVALUE
    IF(letquo(VALUE(NV)).EQ.0) THEN
      CALL ch2int (VALUE(NV),INT)
      CALL nvfills (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
    ELSE
      CALL ch2ch (VALUE(NV),CHVAL,NCW)
      CALL nvfills (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
    ENDIF
  ENDDO
ELSE
  DO NV=1,NVALUE
    IF(letquo(VALUE(NV)).EQ.0) THEN
      CALL ch2real (VALUE(NV),FNUM)
      CALL nvfills (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
    ELSE
      CALL ch2ch (VALUE(NV),CHVAL,NCW)
      CALL nvfills (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
    ENDIF
  ENDDO
ENDIF

return
END SUBROUTINE nvtran_isan
