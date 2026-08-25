!##############################################################################
Subroutine opspec1 ()

! This routine checks the option specifications in the $MODEL_GRIDS
!   namelist for consistency, and overrides settings of ICLOUD,
!   IRAIN, IPRIS, ISNOW, IAGGR, IGRAUP, and IHAIL, setting them
!   all to zero if LEVEL is less than 3.

! Each error will be assigned a severity.  Fatal errors will cause
!   run to stop immediately and warning errors and informative
!   messages will be listed.

use mem_grid
use mem_radiate
use mem_leaf
use micphys
use micro_prm, only:iceflag,iceprocs

implicit none

integer :: ierr,ifm,icm,ng,ifaterr,iwarerr,infoerr,nparent


IFATERR=0
IWARERR=0
INFOERR=0

! Check number of grids and grid points in data versus maximums in
!   configuration file. (severity - F )

IF(NGRIDS.GT.MAXGRDS) THEN
  PRINT*,' FATAL - NGRIDS in data greater than MAXGRDS in config.'
  IFATERR=IFATERR+1
  goto 1000
ENDIF

IERR=0
DO NGRID=1,NGRIDS
  IF(NNXP(NGRID).GT.NXPMAX)IERR=IERR+1
  IF(NNYP(NGRID).GT.NYPMAX)IERR=IERR+1
  IF(NNZP(NGRID).GT.NZPMAX)IERR=IERR+1
ENDDO
IF(NZG.GT.NZGMAX)IERR=IERR+1
IF(IERR.GT.0) THEN
  PRINT*,' FATAL - grid points in data exceeded maximums in',  &
         ' config ',IERR,' times.'
  IFATERR=IFATERR+IERR
ENDIF

! Atmospheric levels must exceed soil levels
!   (severity - F)

DO NGRID=1,NGRIDS
  IF(NNZP(NGRID).LT.NZG+NZS+1) THEN
    PRINT*,' FATAL  - nnzp  must be at least nzg + nzs + 1.'
    IFATERR=IFATERR+1
  ENDIF
ENDDO

! Minimum values of NGRIDS, NNXP, NNYP, NNZP, NSTRATX, NSTRATY
!   (severity - F)

IF(NGRIDS.LT.1) THEN
  PRINT*,' FATAL - NGRIDS in data less than 1.'
  IFATERR=IFATERR+1
ENDIF

DO IFM=1,NGRIDS
  icm = nxtnest(ifm)
  IF(NNXP(IFM).LT.4) THEN
    PRINT*,' FATAL - NNXP must be at least 4.'
    IFATERR=IFATERR+1
  ENDIF

  IF(NNYP(IFM).LT.1) THEN
    PRINT*,' FATAL - NNYP must be at least 1.'
    IFATERR=IFATERR+1
  ENDIF

  IF(NNZP(IFM).LT.11) THEN
    PRINT*,' FATAL - NNZP must be at least 11.'
    IFATERR=IFATERR+1
  ENDIF

  IF(NSTRATX(IFM).LT.1) THEN
    PRINT*,' FATAL - NSTRATX must be at least 1.'
    IFATERR=IFATERR+1
  ENDIF

  !Set NSTRATY based on NSTRATX and 2D vs 3D
  if(nnyp(ifm).eq.1) then
    NSTRATY(IFM) = 1
  else
    NSTRATY(IFM) = NSTRATX(IFM)
  endif

  IF(NSTRATY(IFM).LT.1) THEN
    PRINT*,' FATAL - NSTRATY must be at least 1.'
    IFATERR=IFATERR+1
  ENDIF

  if(nstraty(ifm).ne.1.and.nnyp(ifm).eq.1)then
    print*,' FATAL - nstraty must be set to 1 for 2-d run.'
    ifaterr=ifaterr+1
  endif

  if (icm .ge. 1 .and. nnyp(ifm) .eq. 1 .and.  &
     (ninest(ifm) .lt. 3 .or. njnest(ifm) .lt. 1)) then
     print*, ' FATAL - nested 2d grid must have ninest > 2 '  &
     ,'and njnest = 1 in namelist.'
      IFATERR=IFATERR+1
   ENDIF
ENDDO

! Allowable values of CENTLAT, CENTLON, POLELAT, POLELON
!   (severity - F)

if(polelat.lt.-90..or.polelat.gt.90.) then
   PRINT*,' FATAL - POLELAT outside of legal bounds.'
   IFATERR=IFATERR+1
ENDIF

if(polelon.lt.-180..or.polelon.gt.180.) then
   PRINT*,' FATAL - POLELON outside of legal bounds.'
   IFATERR=IFATERR+1
ENDIF

do ng=1,ngrids
   if(centlat(ng).lt.-90..or.centlat(ng).gt.90.) then
      PRINT*,' FATAL - CENTLAT outside of legal bounds.'
      IFATERR=IFATERR+1
   ENDIF

   if(centlon(ng).lt.-180..or.centlon(ng).gt.180.) then
      PRINT*,' FATAL - CENTLON outside of legal bounds.'
      IFATERR=IFATERR+1
   ENDIF
enddo

! Check nxtnest values for validity
!   (severity - F)

if (nxtnest(1) .ne. 0) then
   print*, ' FATAL - Grid # 1 must have its parent mesh'  &
      ,' designated as Grid # 0 (nxtnest(1) = 0)'
   IFATERR=IFATERR+1
ENDIF

nparent = 0
do ifm = 1,ngrids
   icm = nxtnest(ifm)

   IF (icm .GE. ifm) THEN
      PRINT 1, ifm
      1 FORMAT (' FATAL - Nest #',I3,' has specified parent'  &
               ,' mesh of equal or higher number')
      IFATERR=IFATERR+1
   ENDIF

   if (icm .lt. 0) then
      PRINT 2, ifm
      2 FORMAT (' FATAL - Nest #',I3,' has specified parent'  &
               ,' mesh of number less than 0')
      IFATERR=IFATERR+1
   endif

   if (icm .eq. 0) then
      nparent = nparent + 1
   endif

enddo

if (nparent .gt. 1) then
   print*, ' FATAL - more than 1 grid has grid # 0 specified'  &
         ,' as their parent'
   IFATERR=IFATERR+1
endif

! Check to make sure that top grids have NSTTOP and NSTBOT set to 1
!   (severity - F)

IF (NNSTTOP(1) .NE. 1) THEN
  PRINT*,'NSTTOP not set to 1'
  IFATERR=IFATERR+1
ENDIF

IF (NNSTBOT(1) .NE. 1) THEN
  PRINT*,'NSTBOT not set to 1'
  IFATERR=IFATERR+1
ENDIF

! Check to make sure vegetation patches are appropriate relative to
! total number of patches.
if (nvegpat < 1 .or. nvegpat > npatch-1) then
  print*,'FATAL - NVEGPAT must be between 1 and NPATCH-1'
  IFATERR=IFATERR+1
endif

! Require use of 2 surface patches per grid cell for SCM so as to 
! keep this simple for interface to general parent model. May
! consider allowing more patches in the future if feasible.
if (iscm > 0 .and. npatch > 2) then
  print*,'FATAL - NPATCH must equal 2 for single column output'
  IFATERR=IFATERR+1
endif

!********************************************************************
! MICROPHYSICS AND AEROSOL FLAG CHECKING SECTION
!********************************************************************

!********************************************************************
! CHECK FOR AEROSOL SOURCE AND AEROSOL RADIATION FLAGS
!********************************************************************
 if (iaerosol .lt. 0 .or. iaerosol .gt. 1) THEN
    print*,'FATAL - IAEROSOL OUT OF RANGE: MUST BE 0-1'
    IFATERR = IFATERR + 1
 endif
 if (idust .lt. 0 .or. idust .gt. 2) THEN
    print*,'FATAL - IDUST OUT OF RANGE: MUST BE 0-2'
    IFATERR = IFATERR + 1
 endif
 if (isalt .lt. 0 .or. isalt .gt. 2) THEN
    print*,'FATAL - ISALT OUT OF RANGE: MUST BE 0-2'
    IFATERR = IFATERR + 1
 endif
 if (iabcarb .lt. 0 .or. iabcarb .gt. 1) THEN
    print*,'FATAL - IABCARB OUT OF RANGE: MUST BE 0-1'
    IFATERR = IFATERR + 1
 endif
 if (iaerorad .lt. 0 .or. iaerorad .gt. 1) THEN
    print*,'FATAL - IAERORAD OUT OF RANGE: MUST BE 0-1'
    IFATERR = IFATERR + 1
 endif
 if (iaerorad .eq. 1 .and. (iswrtyp.eq.1 .or. ilwrtyp.eq.1 .or. &
                            iswrtyp.eq.2 .or. ilwrtyp.eq.2 .or. &
                            iswrtyp.eq.4 .or. ilwrtyp.eq.4)) THEN
    print*,'FATAL - Aerosol radiation turned on but will'
    print*,'        only impact shortwave and/or longwave'  
    print*,'        radiation if type is set to 3(Harrington) or 5(RTE-RRTMGP)'
    IFATERR = IFATERR + 1
 endif

! If LEVEL is less than 3, set microphysics parameters to zero.
! If LEVEL equals 3, check for values of microphysics parameters
! that are out of bounds.

if (level .eq. 0) then
   iaerorad = 0
endif

if (level .le. 2) then

   icloud = 0
   irain = 0
   ipris = 0
   isnow = 0
   iaggr = 0
   igraup = 0
   ihail = 0
   idriz = 0
   imbudget = 0
   iccnlev = 0
   iifn = 0

elseif (level .eq. 3) then

!Guidance on hydrometeor shape parameter GNU
 if (icloud .gt. 0 .and. gnu(1) < 4) then
  print*,'FATAL - GNU for Cloud droplets is strongly recommended >= 4'
  IFATERR = IFATERR + 1
 endif
 if (idriz .gt. 0 .and. gnu(8) < 4) then
  print*,'FATAL - GNU for Drizzle droplets is strongly recommended >= 4'
  IFATERR = IFATERR + 1
 endif

!ALLOWABLE RANGE OF MICRO / RADIATION FLAGS
 if (icloud .lt. 0 .or. icloud .gt. 5) THEN
  print*,'FATAL - ICLOUD OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (irain .lt. 0 .or. irain .gt. 5) THEN
  print*,'FATAL - IRAIN OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (ipris .lt. 0 .or. ipris .gt. 5) THEN
  print*,'FATAL - IPRIS OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (ipris .ge. 1 .and. ipris .le. 4) THEN
  print*,'FATAL - IPRIS MUST BE 0 or 5'
  IFATERR = IFATERR + 1
 endif
 if (isnow .lt. 0 .or. isnow .gt. 5) THEN
  print*,'FATAL - ISNOW OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (iifn .eq. 4 .and. itrkdustifn .eq. 1) then
  print*, 'FATAL - ITRKDUSTIFN must be 0 when using simple icenuc (iifn=4)'
  IFATERR = IFATERR + 1
 endif 
 if (iaggr .lt. 0 .or. iaggr .gt. 5) THEN
  print*,'FATAL - IAGGR OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (igraup .lt. 0 .or. igraup .gt. 5) THEN
  print*,'FATAL - IGRAUP OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (ihail .lt. 0 .or. ihail .gt. 5) THEN
  print*,'FATAL - IHAIL OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (idriz .lt. 0 .or. idriz .gt. 5) THEN
  print*,'FATAL - IDRIZ OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (iccnlev .lt. 0 .or. iccnlev .gt. 2) THEN
  print*,'FATAL - ICCNLEV OUT OF RANGE: MUST BE 0-2'
  IFATERR = IFATERR + 1
 endif
 if (imbudget .lt. 0 .or. imbudget .gt. 3) THEN
  print*,'FATAL - IMBUDGET OUT OF RANGE'
  IFATERR = IFATERR + 1
 endif
 if (iifn .lt. 0 .or. iifn .gt. 4) THEN
  print*,'FATAL - IIFN OUT OF RANGE: MUST BE 0-4'
  IFATERR = IFATERR + 1
 endif
 if (isedim .lt. 0 .or. isedim .gt. 1) THEN
  print*,'FATAL - ISEDIM OUT OF RANGE: MUST BE 0-1'
  IFATERR = IFATERR + 1
 endif
 if (itrkepsilon .lt. 0 .or. itrkepsilon .gt. 1 .or. &
     itrkdust    .lt. 0 .or. itrkdust    .gt. 1 .or. &
     itrkdustifn .lt. 0 .or. itrkdustifn .gt. 1) THEN
  print*,'FATAL - AEROSOL TRACKING FLAGS MUST BE 0-1'
  IFATERR = IFATERR + 1
 endif

 if (icloud .lt. 5) then
  if (idriz .ge. 5) then
    print*,'FATAL - Microphysics - IDRIZ must be < 5 if ICLOUD < 5'
    IFATERR = IFATERR + 1
  endif
  if (iccnlev .gt. 0) THEN
    print*,'WARNING - iccnlev>0 and icloud<5.'
    print*,' No aerosol nucleation for these settings.'
    print*,' Aerosol removal/tracking now only active via deposition.'
    IWARERR = IWARERR + 1
  endif
  if (iifn .eq. 3 .and. icloud < 5) THEN
    print*,'FATAL - IIFN=3 requires ICLOUD>=5 for there to be any'
    print*,' primary heterogeneous ice nucleation.'
    IFATERR = IFATERR + 1
  endif
 endif
 if (icloud .eq. 5) then
  if (iaerosol.eq.0 .and. idust.eq.0 .and. isalt.eq.0 .and. iabcarb.eq.0) then
    print*,'FATAL - If icloud==5, then either IAEROSOL, IDUST, IABCARB, or ISALT'
    print*,' must be turned on so that some aerosols are present in order'
    print*,' for cloud nucleation to occur.'
    IFATERR = IFATERR + 1
  endif
 endif

 if (imbudget .ge. 3 .and. iccnlev .eq. 0) THEN
   print*,'FATAL - If IMBUDGET=3 then ICCNLEV must be > 0.'
   print*,' Cannot track specific aerosol type nucleation'
   print*,' unless ICCNLEV>0 for nucleation scavenging.'
   IFATERR = IFATERR + 1
 endif
 if (imbudget .ge. 3 .and. idust .eq. 0) THEN
   print*,'FATAL - If IMBUDGET=3 then IDUST must be > 0.'
   print*,' Cannot track dust if there is no dust turned on.'
   IFATERR = IFATERR + 1
 endif

 if ((itrkepsilon.gt.0 .or. itrkdust.gt.0 .or. itrkdustifn.gt.0) &
      .and. iccnlev.lt.2) THEN
   print*,'FATAL - Aerosol tracking = 0 since ICCNLEV < 2.'
   print*,'Either set ICCNLEV>=2 or set ITRKEPSILON=ITRKDUST=ITRKDUSTIFN=0'
   IFATERR = IFATERR + 1
 endif
 if ((itrkdust.gt.0 .or. itrkdustifn.gt.0) .and. idust.eq.0) THEN
   print*,'FATAL - Dust tracking = 0 since IDUST=0.'
   print*,'Either set IDUST>0 or set ITRKDUST=ITRKDUSTIFN=0'
   IFATERR = IFATERR + 1
 endif
 if ((iaerodep.eq.1 .and. isfcl.eq.3)) THEN
    print*,'FATAL - The surface scheme must be turned on to run aerosol deposition'
    IFATERR = IFATERR + 1
 endif
elseif (level .eq. 4) then
   idriz=0
   if (ipris>0.or.isnow>0.or.iaggr>0.or.igraup>0.or.ihail>0) then
      iceprocs = 1
      if (ipris == 1) print*, 'ONE ICE CRYSTAL TYPE: COLUMNS'
      if (ipris == 2) print*, 'ONE ICE CRYSTAL TYPE: PLATES'
      if (ipris == 3) print*, 'ONE ICE CRYSTAL TYPE: DENDRITES'
      if (ipris == 0) print*, 'NO ICE CRYSTAL TYPES TURNED ON'
      if (igraup == 0 .and. ihail == 0) then
        print*, 'For SBM microphysics, either graupel or hail must be turned on!'
        stop
      endif
      if (igraup == 0) print*, 'HAIL on, GRAUPEL off'
      if (ihail == 0) print*, 'GRAUPEL on, HAIL off'
      ICEFLAG=IIFN-1
      IF(ICEFLAG .ne. 0 .and. ICEFLAG .ne. 1) THEN
        print*, 'IIFN needs to be set to 1 or 2. IIFN=',IIFN
        stop
      ENDIF
   ENDIF
   if (icloud      > 0 .or. &
       idriz       > 0 .or. &
       irain       > 0 .or. &
       isnow       > 0 .or. &
       iaggr       > 0 .or. &
       iaerosol    > 0 .or. &
       idust       > 0 .or. &
       isalt       > 0 .or. &
       iabcarb     > 0 .or. &
       iccnlev     > 0 .or. &
       iaerorad    > 0 .or. &
       iaerodep    > 0 .or. &
       iaeroprnt   > 0 .or. &
       itrkepsilon > 0 .or. &
       itrkdust    > 0 .or. &
       itrkdustifn > 0 .or. &
       icheckmic   > 0) THEN
    print*,'FATAL - The following microphysics ON/OFF switches are not used'
    print*,'in HUCM Bin Microphysics, so set them to zero for OFF.'
    print*,''
    print*,'ICLOUD, IDRIZ, IRAIN, ISNOW, IAGGR, IAEROSOL' 
    print*,'IABCARB, ISALT, IDUST, ICCNLEV, IAERORAD, IAERODEP, IAEROPRNT' 
    print*,'ITRKEPSILON, ITRKDUST, ITRKDUSTIFN, ICHECKMIC.'
    print*,''
    print*,'HUCM only uses the following microphysics flags:'
    print*,''
    print*,'IMBUDGET, IPRIS, IGRAUP, IHAIL, HUCMFILE, NDTCOLL, IIFN'
    print*,'IAEROHIST, CIN_MAX, CCN_MAX, IAERO_CHEM, AERO_MEDRAD'
    print*,''
    print*,'Other microphysics flags that are NOT on/off switches are ignored.'
    print*,''
    print*,'The LEVEL=3 microphysics aerosol categories are not used in HUCM.'
    print*,'So aerosol radiation and deposition routines are not used.'
    print*,'CCN and INP are specified in CCN_MAX and CIN_MAX for HUCM.'
    IFATERR = IFATERR + 1
   endif
   if (imbudget==3)then
    print*,'FATAL - No dust nucleation tracking because there is'
    print*,'no dust in HUCM-SBM. Set IMBUDGET < 3.'
    IFATERR = IFATERR + 1
   endif


endif
!********************************************************************
! END MICROPHYSICS AND AEROSOL FLAG CHECKING
!********************************************************************

1000 continue

! Stop the run if there are any fatal errors.  List how many
!   warning and informative errors.

PRINT*,' -----------opspec1--------------------------'
PRINT*,' FATAL     errors - ',IFATERR
PRINT*,' WARNING   errors - ',IWARERR
PRINT*,' INFORM  messages - ',INFOERR
PRINT*,' -----------------------------------------------'

IF(IFATERR.GT.0) STOP 'OPSPEC1'

return
END SUBROUTINE opspec1

!##############################################################################
Subroutine opspec2 ()

! Check that fine mesh is a valid subset of its coarser mesh.
!   and that top and bottom boundary flags are set correctly.
!   (severity - F)

use mem_grid
use mem_varinit

implicit none

integer :: icm,ifm,ifaterr,iwarerr,infoerr,ncx,ncy,nfx,nfxp,nfy,nfyp
integer :: nesta,nfz,kc

IFATERR=0
IWARERR=0
INFOERR=0

DO ifm=1,NGRIDS

   icm = nxtnest(ifm)
   if (icm .ge. 1) then

      NCX=(NNXP(ifm)-2)/NSTRATX(ifm)
      NCY=(NNYP(ifm)-2)/NSTRATY(ifm)

      IF ((NNYP(ifm).EQ.1.AND.NNYP(icm).NE.1).OR.  &
          (NNYP(ifm).NE.1.AND.NNYP(icm).EQ.1)) THEN
         PRINT*,' FATAL - Grids must be either all 3-D or all 2-D'
         IFATERR=IFATERR+1
      ENDIF

      IF (NINEST(ifm).LT.3) THEN
         PRINT 11, ifm
         11 FORMAT (' FATAL - Nest #',I3,' too close to western'  &
                   ,' boundary of coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

      IF (NJNEST(ifm).LT.3.AND.NNYP(ifm).GT.1) THEN
         PRINT 12, ifm
         12 FORMAT (' FATAL - Nest #',I3,' too close to southern'  &
                   ,' boundary of coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

      IF (NKNEST(ifm).LT.3.AND.NNSTBOT(ifm).EQ.0) THEN
         PRINT 13, ifm
         13 FORMAT (' FATAL - Nest #',I3,' too close to lower'  &
                   ,' boundary of coarser mesh or NNSTBOT incorrect')
         IFATERR=IFATERR+1
      ENDIF

      IF (NKNEST(ifm).NE.1.AND.NNSTBOT(ifm).EQ.1) THEN
         PRINT 14, ifm
         14 FORMAT (' FATAL - Nest #',I3,' not to lower boundary of'  &
                   ,' coarser mesh or NNSTBOT flag set incorrectly')
         IFATERR=IFATERR+1
      ENDIF

      IF (NINEST(ifm)+NCX.GT.NNXP(icm)-3) THEN
         PRINT 15, ifm
         15 FORMAT (' FATAL - Nest #',I3,' too close to eastern'  &
                   ,' boundary of coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

      IF (NJNEST(ifm)+NCY.GT.NNYP(icm)-3.AND.NNYP(ifm).GT.1) THEN
         PRINT 16, ifm
         16 FORMAT (' FATAL - Nest #',I3,' too close to northern'  &
                   ,' boundary of coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

      IF (NCX.LT.2.OR.(NCY.LT.2.AND.NNYP(ifm).NE.1)) THEN
         PRINT 17, ifm
         17 FORMAT (' FATAL - Nest #',I3,' dimensioned too small'  &
                   ,' in at least one horizontal direction')
         IFATERR=IFATERR+1
      ENDIF

      NFX=NCX*NSTRATX(ifm)+2
      IF (NNXP(ifm).NE.NFX) THEN
         NFXP=NFX+NSTRATX(ifm)
         PRINT 18, ifm,NFXP,NFX
         18 FORMAT (' FATAL - Nest #',I3,' NNXP incompatible with'  &
                   ,' NSTRATX:  May increase NNXP to',I5  &
                   ,' or decrease it to',I5)
         IFATERR=IFATERR+1
      ENDIF

      NFY=NCY*NSTRATY(ifm)+2
      IF (NNYP(ifm).NE.NFY.AND.NNYP(ifm).GT.1) THEN
         NFYP=NFY+NSTRATY(ifm)
         PRINT 19, ifm,NFYP,NFY
         19 FORMAT (' FATAL - Nest #',I3,' NNYP incompatible with'  &
                   ,' NSTRATY:  May increase NNYP to',I5  &
                   ,' or decrease it to',I5)
         IFATERR=IFATERR+1
      ENDIF

      IF (NNSTBOT(ifm).EQ.1.AND.NNSTBOT(icm).NE.1) THEN
         PRINT 20, ifm
         20 FORMAT (' FATAL - Nest #',I3,' NNSTBOT flag incompatible'  &
                   ,' with coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

      IF (NNSTTOP(ifm).EQ.1.AND.NNSTTOP(icm).NE.1) THEN
         PRINT 21, ifm
         21 FORMAT (' FATAL - Nest #',I3,' NNSTTOP flag incompatible'  &
                   ,' with coarser mesh')
         IFATERR=IFATERR+1
      ENDIF

   endif
ENDDO

nesta=abs(nestz)
if(nestz.ne.0.and.nesta.le.ngrids)then
   nfz=nnzp(nesta)-2
   kc=nknest(nesta)
     1002    continue
      kc=kc+1
      nfz=nfz-nrz(kc,nesta)
      if(nfz.lt.0)then
         print 195,nesta
         195 format(' FATAL - vertically nested grid #',i3,  &
                    ' has illegal number of levels for given NSTRATZ values')
         IFATERR=IFATERR+1
      ENDIF
   if(nfz.gt.0)go to 1002
   if(nfz.eq.0)then
      if(kc.gt.nnzp(nxtnest(nesta))-3.and.nnsttop(nesta).ne.1)then
         PRINT 22, nesta
         22 FORMAT (' FATAL - Nest #',I3,' too high'  &
                   ,' or NNSTTOP flag set incorrectly')
         IFATERR=IFATERR+1
      ENDIF

      if(kc.ne.nnzp(nxtnest(nesta))-1.and.nnsttop(nesta).eq.1)then
         PRINT 23, nesta
         23 FORMAT (' FATAL - Nest #',I3,' not to upper boundary of'  &
                   ,' coarser mesh or NNSTTOP flag set incorrectly')
         IFATERR=IFATERR+1
      ENDIF
   endif
endif

DO ifm = 1,NGRIDS
   icm = nxtnest(ifm)
   if (ifm < nesta .and. icm .ge. 1)then
      if(nnzp(ifm).gt.nnzp(icm)-nknest(ifm)-1.and.nnsttop(ifm).eq.0)then
         PRINT 24, ifm
         24 FORMAT (' FATAL - Nest #',I3,' NNZP incompatible with'  &
                   ,' parent grid or NNSTTOP flag set incorrectly')
         IFATERR=IFATERR+1
      ENDIF

      if(nnzp(ifm).ne.nnzp(icm)-nknest(ifm)+1.and.nnsttop(ifm).eq.1)then
         PRINT 25, ifm
         25 FORMAT (' FATAL - Nest #',I3,' not to upper boundary of'  &
                   ,' coarser mesh or NNSTTOP flag set incorrectly')
         IFATERR=IFATERR+1
      ENDIF
   endif
enddo

! This need to be done here since VARFILES are filled before OPSPEC3
if (VWAITTOT.lt.VWAIT1) then
   print*,'Total wait time must be <= individual VARFILE wait'
   print*,'      resetting VWAITTOT to ',VWAIT1
   VWAITTOT=VWAIT1
   IWARERR=IWARERR+1
endif

! Stop the run if there are any fatal errors.  List how many
!   warning and informative errors.

PRINT*,' -----------opspec2--------------------------'
PRINT*,' FATAL     errors - ',IFATERR
PRINT*,' WARNING   errors - ',IWARERR
PRINT*,' INFORM  messages - ',INFOERR
PRINT*,' -----------------------------------------------'

IF(IFATERR.GT.0) STOP 'OPSPEC2'

return
END SUBROUTINE opspec2

!##############################################################################
Subroutine opspec3 ()

use mem_varinit
use mem_grid
use micphys
use io_params
use mem_radiate
use mem_cuparm
use mem_turb
use mem_leaf
use mem_sib, only:nzg_sib
use kpp_parameters, only:IKPP,FRQKPP

implicit none

integer :: k,ifaterr,iwarerr,infoerr,ng,ngr

IFATERR=0
IWARERR=0
INFOERR=0

! Check that moisture is turned on if radiation is used.
!   (severity - F)

  IF(ilwrtyp+iswrtyp.GT.0.AND.LEVEL.EQ.0)THEN
    PRINT*,' FATAL  - radiation scheme must be run with moisture.'
    IFATERR=IFATERR+1
  ENDIF
  
! Microphysics flags and parameter settings

IF( ((icloud .GE. 2 .AND. icloud .LE. 4 .AND. cparm .LE. 0.)  &
 .OR.(idriz  .GE. 2 .AND. idriz  .LE. 4 .AND. dparm .LE. 0.)  &
 .OR.(irain  .GE. 2 .AND. irain  .LE. 4 .AND. rparm .LE. 0.)  &
 .OR.(ipris  .GE. 2 .AND. ipris  .LE. 4 .AND. pparm .LE. 0.)  &
 .OR.(isnow  .GE. 2 .AND. isnow  .LE. 4 .AND. Sparm .LE. 0.)  &
 .OR.(igraup .GE. 2 .AND. igraup .LE. 4 .AND. gparm .LE. 0.)  &
 .OR.(iaggr  .GE. 2 .AND. iaggr  .LE. 4 .AND. Aparm .LE. 0.)  &
 .OR.(ihail  .GE. 2 .AND. ihail  .LE. 4 .AND. hparm .LE. 0.)) &
 .AND. level <=3) THEN
   print*,' '
   print*,'FATAL - Microphysics - xPARM must be positive'  &
             ,'if micro flags are set to 2, 3, or 4:'
   print*,'cparm:',cparm
   print*,'dparm:',dparm
   print*,'rparm:',rparm
   print*,'pparm:',pparm
   print*,'sparm:',sparm
   print*,'aparm:',aparm
   print*,'gparm:',gparm
   print*,'hparm:',hparm
   print*,' '
   IFATERR=IFATERR+1
ENDIF

IF(ibubble.lt.0 .or. ibubble.gt.3) then
   PRINT*,' FATAL - IBUBBLE must be 0, 1, 2, or 3'
   PRINT*,'         0 = off'
   PRINT*,'         1 = RAMSIN-set square bubble'
   PRINT*,'         2 = RAMSIN-set gaussian bubble'
   PRINT*,'         3 = Random bubble in ruser'
   IFATERR=IFATERR+1
ENDIF

IF(iconv.lt.0 .or. iconv.gt.5) then
   PRINT*,' FATAL - ICONV must be 0 - 5'
   PRINT*,'         0 = off'
   PRINT*,'         1-5 = On with gaussian shape'
   IFATERR=IFATERR+1
ENDIF

! Convective parameterization flags and parameter settings

DO NG=1,NGRIDS
  IF(NNQPARM(NG).GT.0.AND.LEVEL.EQ.0) THEN
    PRINT 27
    27 FORMAT (' FATAL - LEVEL must be at least'  &
              ,' 1 for the cumulus parameterization')
    IFATERR=IFATERR+1
  ENDIF
ENDDO

DO NG=1,NGRIDS
  IF(NNQPARM(NG) == 2 .AND.LEVEL /= 3) THEN
    PRINT 127
    127 FORMAT (' FATAL - LEVEL must be set to'  &
              ,' 3 for the Kain-Fritsch cumulus parameterization')
    IFATERR=IFATERR+1
  ENDIF
ENDDO

! Check horizontal and vertical grid spacings.

IF(DZMAX.LT.DELTAZ)THEN
  PRINT*,' WARNING - DELTAZ is being reduced by a low value',  &
    ' of DZMAX.'
  IWARERR=IWARERR+1
ENDIF

IF(DZRAT.GT.1.2)THEN
  PRINT*,' WARNING - Large vertical stretch ratios sacrifice',  &
    ' second order accuracy in the vertical differencing.'
  IWARERR=IWARERR+1
ENDIF

! Check numerical schemes.

IF((SSPCT.LT.0.2.OR.SSPCT.GT.1.0).AND.SSPCT.NE.0.0)THEN
  PRINT*,' WARNING - SSPCT should normally range from 0.2 to 1.0',sspct
  IWARERR=IWARERR+1
ENDIF

IF(NFPT.GT.NNZP(1))THEN
  PRINT*,' FATAL - NFPT must be less than nnzp(1).'
  IFATERR=IFATERR+1
ENDIF

! Check turbulence parameterization.

DO NGR=1,NGRIDS
  IF(IDIFFK(NGR).LT.1.OR.IDIFFK(NGR).GT.4)THEN
    PRINT*,' FATAL - IDIFFK must be 1, 2, 3, or 4.'
    IFATERR=IFATERR+1
  ENDIF
ENDDO

! Check that diffusion flags are compatible if using ihorgrad=1

if(ihorgrad.eq.2)then
  if(IDIFFK(NGR) >= 3)then
    print*,' FATAL - Cant use IHORGRAD=2 if IDIFFK >= 3'
    IFATERR=IFATERR+1
  endif
endif

! Check that diffusion of perturbations relative to varfile
! state only runs if varfiles are used for nudging.

if(idiffperts==3 .and. nud_type==0) then
   print*,' FATAL - Cant use IDIFFPERTS=3 without varfiles'
   IFATERR=IFATERR+1
endif

! Check whether the soil model will be run and make sure that the
!   number of soil levels are correct.(severity - F,I )

if(isoildat < 0 .or. isoildat > 1)then
  PRINT*,' SOIL INITIALIZATION OPTION IS NEEDED: ISOILDAT = 0 or 1'
  PRINT*,' ISOILDAT = 0 (for default horizontal homogeneous run)'
  PRINT*,' ISOILDAT = 1 (for soil varfile initialization)'
  ifaterr = ifaterr + 1
endif
if(isoildat == 1)then
  if (initial == 1 .or. (initial == 3 .and. initorig == 1)) then
   PRINT*,' VARIABLE SOIL INITIALIZATION FOR INITIAL=2 ONLY'
   ifaterr = ifaterr + 1
  endif
endif
if(isnowdat < 0 .or. isnowdat > 1)then
  PRINT*,' SNOW INITIALIZATION OPTION IS NEEDED: ISNOWDAT = 0 or 1'
  PRINT*,' ISNOWDAT = 0 (for default horizontal homogeneous run)'
  PRINT*,' ISNOWDAT = 1 (for snow varfile initialization)'
  ifaterr = ifaterr + 1
endif
if(isnowdat == 1)then
  if (initial == 1 .or. (initial == 3 .and. initorig == 1)) then
   PRINT*,' VARIABLE SNOW INITIALIZATION FOR INITIAL=2 ONLY'
   ifaterr = ifaterr + 1
  endif
endif

IF(ISFCL.EQ.0.AND.NZG.GT.1)THEN
  PRINT*,' INFO  - more soil levels specified than needed.'
  INFOERR=INFOERR+1
ENDIF

if (isfcl == 0 .and. npatch /= 2) then
  print*, ' fatal  - When isfcl = 0, npatch must be 2. '
  ifaterr = ifaterr + 1
endif

IF(ISFCL.GT.0.AND.ISFCL.LT.3.AND.NZG.LE.2)THEN
  PRINT*,  &
    ' FATAL  - at least 3 soil levels are needed for soil'  &
   ,' model.'
  IFATERR=IFATERR+1
ENDIF

!If running SiB surface model, can only have 1 surface water layer
IF(ISFCL.EQ.2.AND.NZS.NE.1)THEN
  PRINT*,  &
    ' FATAL  - For running SiB, set NZS=1 to match SiB layers'  &
   ,' model.'
  IFATERR=IFATERR+1
ENDIF

!If running SiB surface model, NZG=NZG_SIB
IF(ISFCL.EQ.2.AND.NZG.NE.nzg_sib)THEN
  PRINT*,  &
    ' FATAL  - For running SiB, set NZG=7 to match SiB layers'  &
   ,' model.'
  IFATERR=IFATERR+1
ENDIF

do k=1,nzg
   if (slz(k) .gt. -.001) then
      print*, 'FATAL - Soil level',k,' not (enough) below ground'  &
         ,' level'
     IFATERR=IFATERR+1
   endif
enddo

do k=1,nzg-1
   if (slz(k)-slz(k+1) .gt. .001) then
      print*, 'FATAL - Soil level',k,' not (enough) deeper than'  &
      ,' soil level',k+1
     IFATERR=IFATERR+1
   endif
enddo

! If the soil model will be run with no radiation, make a suggestion
!   that the radiation be turned on. (severity - F )

DO NGR=1,NGRIDS
  IF(ISFCL.GT.0.AND.ilwrtyp+iswrtyp.EQ.0)THEN
    PRINT*,' FATAL  - radiation scheme must be run with soil',  &
        ' model.'
    IFATERR=IFATERR+1
  ENDIF
ENDDO

! If running KPP ocean model, make sure FRQKPP >= coarse grid timestep.
if (IKPP > 0 .and. FRQKPP < dtlongn(1)) then
  print*, 'KPP timestep shorter than DTLONG for coard grid.'
  print*, 'KPP timestep must be >= DTLONG.'
  IFATERR=IFATERR+1
endif

! Make sure that if nudging, nudging time scales are greater than
! the model coarse grid timestep, and that Rayleigh friction nudging
! is not done with variable initialization.

if (initial == 1 .or. (initial == 3 .and. initorig == 1)) then
   if (nfpt .gt. 0 .and. distim .gt. 0. .and.  &
                         distim .lt. dtlongn(1)) then
      print*, 'Rayleigh friction nudging is being done'
      print*, 'and DISTIM is less than DTLONGN(1).'
      print*, 'This nudging is too strong.'
      IFATERR=IFATERR+1
   endif
endif

if (initial == 2 .or. (initial == 3 .and. initorig == 2)) then

   if (nfpt .gt. 0 .and. distim .gt. 0.) then
      print*, 'Rayleigh friction nudging may not be used when'
      print*, 'variable initialization is used.'
      IFATERR=IFATERR+1
   endif

   if (nudlat .ge. 1 .and. tnudlat .gt. 0. .and.  &
                           tnudlat .lt. dtlongn(1)) then
      print*, 'Lateral boundary nudging is being done'
      print*, 'and TNUDLAT is less than DTLONGN(1).'
      print*, 'This nudging is too strong.'
      IFATERR=IFATERR+1
   endif

   if (tnudcent .gt. 0. .and. tnudcent .lt. dtlongn(1)) then
      print*, 'Center nudging is being done'
      print*, 'and TNUDCENT is less than DTLONGN(1).'
      print*, 'This nudging is too strong.'
      IFATERR=IFATERR+1
   endif

   if (tnudtop .gt. 0. .and. tnudtop .lt. dtlongn(1)) then
      print*, 'Top boundary nudging is being done'
      print*, 'and TNUDTOP is less than DTLONGN(1).'
      print*, 'This nudging is too strong.'
      IFATERR=IFATERR+1
   endif

   if (snudcent .gt. 0. .and. snudcent .lt. dtlongn(1)) then
      print*, 'Center nudging is being done for soil moisture'
      print*, 'and SNUDCENT is less than DTLONGN(1).'
      print*, 'This nudging is too strong.'
      IFATERR=IFATERR+1
   endif

endif


!     Check the averaging and analysis frequencies for consistency.

if (abs(AVGTIM).gt.0.0.and.FRQMEAN.le.0.0.and.FRQBOTH.le.0.) then
   print*,'Have FRQMEAN=0 & FRQBOTH=0 even though AVGTIM=',AVGTIM
   print*,'Respecifying AVGTIM=0.'
   AVGTIM=0.
   IWARERR=IWARERR+1
endif
if (FRQMEAN.gt.0.0.and.abs(AVGTIM).gt.0.) then
   if ( abs(AVGTIM).gt.FRQMEAN ) then
      print*,'AVGTIM must be <= FRQMEAN'
      IFATERR=IFATERR+1
   endif
endif
if (FRQBOTH.gt.0.0.and.abs(AVGTIM).gt.0.) then
   if ( abs(AVGTIM).gt.FRQBOTH ) then
      print*,'AVGTIM must be <= FRQBOTH'
      IFATERR=IFATERR+1
   endif
endif
if (FRQMEAN.gt.0.0.and.FRQBOTH.gt.0.0.and.abs(AVGTIM).gt.0.) then
   if ( (FRQMEAN.gt.FRQBOTH.and.mod(FRQMEAN,FRQBOTH).ne.0.).or.  &
        (FRQMEAN.lt.FRQBOTH.and.mod(FRQBOTH,FRQMEAN).ne.0.) ) then
      print*,'FRQMEAN must be a multiple of FRQBOTH or vice versa'
      IFATERR=IFATERR+1
   endif
endif

! Stop the run if there are any fatal errors.  List how many
!   warning and informative errors.

PRINT*,' -----------opspec3--------------------------'
PRINT*,' FATAL     errors - ',IFATERR
PRINT*,' WARNING   errors - ',IWARERR
PRINT*,' INFORM  messages - ',INFOERR
PRINT*,' -----------------------------------------------'

IF(IFATERR.GT.0) STOP 'OPSPEC3'

return
END SUBROUTINE opspec3

