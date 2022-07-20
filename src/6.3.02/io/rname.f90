!##############################################################################
Subroutine nvfillm (group,vr,ii,ff,cc,nv)

!Namelist reading routines for MODEL

! Keep the list of variables here in sync with the list in broadcast_config() in mpass_init.f90

use mem_all
use node_mod
use leaf_coms, only:ifreeslip

implicit none

! This routine places values read in for a particular NAMELIST
! variable in the storage location(s) for that variable
!
! GROUP   - NAMELIST group name
! VR      - NAMELIST variable name
! II      - vector of integer values
! FF      - vector of real values
! CC      - vector of character strings
! NV      - number of values to be inserted

character(len=*) :: group,vr,cc
real :: ff
integer :: ii,nv
integer :: inrflg
integer, parameter ::nvgrid=37,nvstrt=77,nvindat=145,nvsound=10
integer ::  igrids(nvgrid),istart(nvstrt),iindat(nvindat),isound(nvsound)
character(len=16) :: grids(nvgrid),start(nvstrt),indat(nvindat),sound(nvsound)
data igrids/nvgrid*0/,istart/nvstrt*0/,iindat/nvindat*0/,isound/nvsound*0/

DATA GRIDS/  &
      'EXPNME','RUNTYPE','TIMEUNIT','TIMMAX','IMONTH1','IDATE1','IYEAR1' &
     ,'ITIME1','NGRIDS','NNXP','NNYP','NNZP','NZG','NZS','NXTNEST'       &
     ,'IPRNTSTMT','IHTRAN','DELTAX','DELTAZ','DZRAT'                     &
     ,'DZMAX','ZZ','DTLONG','NACOUST','NSTRATX'                          &
     ,'NNDTRAT','NESTZ','NSTRATZ','POLELAT','POLELON','NINEST','NJNEST'  &
     ,'NKNEST','CENTLAT','CENTLON','NNSTTOP','NNSTBOT'/
DATA START/  &
      'INITIAL','NUD_TYPE','VARFPFX','VWAIT1','VWAITTOT'                 &
     ,'NUDLAT','TNUDLAT','TNUDCENT','TNUDTOP','ZNUDTOP','WT_NUDGE_G'	 &
     ,'WT_NUDGE_UV','WT_NUDGE_TH','WT_NUDGE_PI','WT_NUDGE_RT'            &
     ,'NUD_COND','TCOND_BEG','TCOND_END','T_NUDGE_RC','WT_NUDGEC'        &
     ,'IF_ODA','ODA_UPAPREF','ODA_SFCPREF','FRQODA','TODABEG'            &
     ,'TODAEND','TNUDODA','WT_ODA_GRID','WT_ODA_UV','WT_ODA_TH'          &
     ,'WT_ODA_PI','WT_ODA_RT','RODA_SFCE','RODA_SFC0','RODA_UPAE'        &
     ,'RODA_UPA0','RODA_HGT','RODA_ZFAC','ODA_SFC_TIL','ODA_SFC_TEL'     &
     ,'ODA_UPA_TIL','ODA_UPA_TEL','HFILIN','IPAST_SFC','ICLOBBER'        &
     ,'IOUTPUT','AFILEPREF','FRQSTATE','FRQST_KEEP','FRQLITE'            &
     ,'NLITE_VARS','LITE_VARS','AVGTIM','FRQMEAN','FRQBOTH','TOPFILES'   &
     ,'SFCFILES','SSTFPFX','NDVIFPFX','ITOPTFLG','ISSTFLG','IVEGTFLG'    &
     ,'ISOILFLG','NDVIFLG','IUPDNDVI','IUPDSST','ITOPTFN'                &
     ,'ISSTFN','IVEGTFN','ISOILFN','NDVIFN','ITOPSFLG','TOPTENH'         &
     ,'TOPTWVL','IZ0FLG','Z0MAX','Z0FACT'/
DATA INDAT/  &
      'ICORFLG','IBND','JBND','CPHAS','LSFLG','NFPT','DISTIM','ISWRTYP'  &
     ,'ILWRTYP','RADFRQ','LONRAD','NNQPARM','CONFRQ','WCLDBS','IKPP'     &
     ,'NKPPZ','FRQKPP','RELAX_SST','RELAX_OCNT','RELAX_SAL','DMAXKPP'    &
     ,'DSCALEKPP','KPPITERMAX','KPPRNT','UBMN_KPP','NPATCH','NVEGPAT'    &
     ,'ISFCL','IFREESLIP','SIBFILE','CO2_INIT','ISOILDAT','SNUDCENT'     &
     ,'ISNOWDAT','NVGCON','PCTLCON','NSLCON','ZROUGH','ALBEDO','SEATMP'  &
     ,'DTHCON','DRTCON','SLZ','SLMSTR','STGOFF','IDIFFK','IDIFFPERTS'    &
     ,'IHORGRAD','CSX','CSZ','XKHKM','ZKHKM','AKMIN','IBUBBLE','IBUBGRD' &
     ,'IBDXIA','IBDXIZ','IBDYJA','IBDYJZ','IBDZK1','IBDZK2','BTHP'       &
     ,'BRTP','ICONV','ICONGR','ICICENT','ICJCENT','CXRAD','CYRAD'        &
     ,'ICVERT','ICKMAX','CZRAD','ICKCENT','CDIVMAX','CTAU','CTMAX'       &
     ,'IRCE','RCE_SZEN','RCE_SOLC','RCE_UBMN','RCE_BUBL','LEVEL','ISCM'  &
     ,'ICHECKMIC','ITRACER','ITRACHIST','IMBUDGET','IRIME','IPLAWS'      &
     ,'ISEDIM','ICLOUD','IDRIZ','IRAIN','IPRIS','ISNOW','IAGGR'          &
     ,'IGRAUP','IHAIL','CPARM','DPARM','RPARM','PPARM','SPARM','APARM'   &
     ,'GPARM','HPARM','GNU','HUCMFILE','NDTCOLL','IAEROSOL','ISALT'      &
     ,'IDUST','IDUSTLOFT','DUSTFILE','ICCNLEV','IIFN','IIFN_FORMULA'     &
     ,'IAERORAD','IAERODEP','IAEROPRNT','IAEROHIST','CIN_MAX','CCN_MAX'  &
     ,'GCCN_MAX','DUST1_MAX','DUST2_MAX','SALTF_MAX','SALTJ_MAX'         &
     ,'SALTS_MAX','IAEROLBC','ICO2LBC','BCTAU','IAERO_CHEM'              &
     ,'AERO_EPSILON','AERO_MEDRAD','ITRKEPSILON','ITRKDUST'              &
     ,'ITRKDUSTIFN','SCMTIME','ISCMX','ISCMY','FRACSAT','IABCARB'        &
     ,'ABC1_MAX','ABC2_MAX'/
DATA SOUND/  &
      'IPSFLG','ITSFLG','IRTSFLG','IUSFLG','HS','PS','TS','RTS','US','VS'/

!##############################################################################
interface
!##############################################################################
Subroutine varsetf (c1,f1,i1,i2,f2,f3,f4)
   
implicit none

   character(len=*) :: c1
   real :: f1,f2,f3,f4
   integer :: i1,i2

END SUBROUTINE varsetf
!##############################################################################
Subroutine varseti (c1,i1,i2,i3,i4,i5,i6)

implicit none

   character(len=*) :: c1
   integer :: i1,i2,i3,i4,i5,i6

END SUBROUTINE varseti
!##############################################################################
End Interface
!##############################################################################

!print*,'vars:',GROUP,'---',trim(VR),'---',ii,ff,cc,nv

IF(GROUP.EQ.'$MODEL_GRIDS') THEN
 CALL varchk (VR,GROUP,GRIDS,IGRIDS,NVGRID,INRFLG)
 IF(INRFLG.EQ.1) RETURN
 IF(VR.EQ.'EXPNME')    CALL varsetc (VR,EXPNME,NV,1,CC,1,strl1)
 IF(VR.EQ.'RUNTYPE')   CALL varsetc (VR,RUNTYPE,NV,1,CC,1,16)
 IF(VR.EQ.'TIMEUNIT')  CALL varsetc (VR,TIMEUNIT,NV,1,CC,1,10)
 IF(VR.EQ.'TIMMAX')    CALL varsetf (VR,TIMMAX,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'IMONTH1')   CALL varseti (VR,IMONTH1,NV,1,II,1,12)
 IF(VR.EQ.'IDATE1')    CALL varseti (VR,IDATE1,NV,1,II,1,31)
 IF(VR.EQ.'IYEAR1')    CALL varseti (VR,IYEAR1,NV,1,II,0,9999)
 IF(VR.EQ.'ITIME1')    CALL varseti (VR,ITIME1,NV,1,II,0,2359)
 IF(VR.EQ.'NGRIDS')    CALL varseti (VR,NGRIDS,NV,1,II,1,9)
 IF(VR.EQ.'NNXP')      CALL varseti (VR,NNXP(NV),NV,MAXGRDS,II,1,10000)
 IF(VR.EQ.'NNYP')      CALL varseti (VR,NNYP(NV),NV,MAXGRDS,II,1,10000)
 IF(VR.EQ.'NNZP')      CALL varseti (VR,NNZP(NV),NV,MAXGRDS,II,1,10000)
 IF(VR.EQ.'NZG')       CALL varseti (VR,NZG,NV,1,II,1,10000)
 IF(VR.EQ.'NZS')       CALL varseti (VR,NZS,NV,1,II,1,10000)
 IF(VR.EQ.'NXTNEST')   CALL varseti (VR,NXTNEST(NV),NV,MAXGRDS,II,0,20)
 IF(VR.EQ.'IPRNTSTMT') CALL varseti (VR,IPRNTSTMT,NV,1,II,0,2)
 IF(VR.EQ.'IHTRAN')    CALL varseti (VR,IHTRAN,NV,1,II,0,1)
 IF(VR.EQ.'DELTAX')    CALL varsetf (VR,DELTAX,NV,1,FF,.001,1.E6)
 IF(VR.EQ.'DELTAZ')    CALL varsetf (VR,DELTAZ,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'DZRAT')     CALL varsetf (VR,DZRAT,NV,1,FF,.1,10.)
 IF(VR.EQ.'DZMAX')     CALL varsetf (VR,DZMAX,NV,1,FF,.001,1.E5)
 IF(VR.EQ.'ZZ')        CALL varsetf (VR,ZZ(NV),NV,NZPMAX,FF,0.,1.E6)
 IF(VR.EQ.'DTLONG')    CALL varsetf (VR,DTLONG,NV,1,FF,0.,1.E8)
 IF(VR.EQ.'NACOUST')   CALL varseti (VR,NACOUST,NV,1,II,1,20)
 IF(VR.EQ.'NSTRATX')   CALL varseti (VR,NSTRATX(NV),NV,MAXGRDS,II,1,20)
 IF(VR.EQ.'NNDTRAT')   CALL varseti (VR,NNDTRAT(NV),NV,MAXGRDS,II,0,20)
 IF(VR.EQ.'NESTZ')     CALL varseti (VR,NESTZ,NV,1,II,-100,100)
 IF(VR.EQ.'NSTRATZ')   CALL varseti (VR,NSTRATZ(NV),NV,NZPMAX,II,1,10)
 IF(VR.EQ.'POLELAT')   CALL varsetf (VR,POLELAT,NV,1,FF,-90.,90.)
 IF(VR.EQ.'POLELON')   CALL varsetf (VR,POLELON,NV,1,FF,-180.,180.)
 IF(VR.EQ.'NINEST')    CALL varseti (VR,NINEST(NV),NV,MAXGRDS,II,0,10000)
 IF(VR.EQ.'NJNEST')    CALL varseti (VR,NJNEST(NV),NV,MAXGRDS,II,0,10000)
 IF(VR.EQ.'NKNEST')    CALL varseti (VR,NKNEST(NV),NV,MAXGRDS,II,1,10000)
 IF(VR.EQ.'CENTLAT')   CALL varsetf (VR,CENTLAT(NV),NV,MAXGRDS,FF,-90.,90.)
 IF(VR.EQ.'CENTLON')   CALL varsetf (VR,CENTLON(NV),NV,MAXGRDS,FF,-180.,180.)
 IF(VR.EQ.'NNSTTOP')   CALL varseti (VR,NNSTTOP(NV),NV,MAXGRDS,II,0,1)
 IF(VR.EQ.'NNSTBOT')   CALL varseti (VR,NNSTBOT(NV),NV,MAXGRDS,II,0,1)
ENDIF

IF(GROUP.EQ.'$MODEL_FILE_INFO') THEN
 CALL varchk (VR,GROUP,START,ISTART,NVSTRT,INRFLG)
 IF(INRFLG.EQ.1) RETURN
 IF(VR.EQ.'INITIAL')     CALL varseti (VR,INITIAL,NV,1,II,1,3)
 IF(VR.EQ.'NUD_TYPE')    CALL varseti (VR,NUD_TYPE,NV,1,II,0,1)
 IF(VR.EQ.'VARFPFX')     CALL varsetc (VR,VARFPFX,NV,1,CC,1,strl1)
 IF(VR.EQ.'VWAIT1')      CALL varsetf (VR,VWAIT1,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'VWAITTOT')    CALL varsetf (VR,VWAITTOT,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'NUDLAT')      CALL varseti (VR,NUDLAT,NV,1,II,0,10000)
 IF(VR.EQ.'TNUDLAT')     CALL varsetf (VR,TNUDLAT,NV,1,FF,0.,1.e10)
 IF(VR.EQ.'TNUDCENT')    CALL varsetf (VR,TNUDCENT,NV,1,FF,0.,1.e10)
 IF(VR.EQ.'TNUDTOP')     CALL varsetf (VR,TNUDTOP,NV,1,FF,0.,1.e10)
 IF(VR.EQ.'ZNUDTOP')     CALL varsetf (VR,ZNUDTOP,NV,1,FF,0.,1.e10)
 IF(VR.EQ.'WT_NUDGE_G')  CALL varsetf (VR,WT_NUDGE_G(NV),NV,MAXGRDS,FF,0.,10.0)
 IF(VR.EQ.'WT_NUDGE_UV') CALL varsetf (VR,WT_NUDGE_UV,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_NUDGE_TH') CALL varsetf (VR,WT_NUDGE_TH,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_NUDGE_PI') CALL varsetf (VR,WT_NUDGE_PI,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_NUDGE_RT') CALL varsetf (VR,WT_NUDGE_RT,NV,1,FF,0.,10.0)
 IF(VR.EQ.'NUD_COND')    CALL varseti (VR,NUD_COND,NV,1,II,0,1)
 IF(VR.EQ.'TCOND_BEG')   CALL varsetf (VR,TCOND_BEG,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'TCOND_END')   CALL varsetf (VR,TCOND_END,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'T_NUDGE_RC')  CALL varsetf (VR,T_NUDGE_RC,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'WT_NUDGEC')   CALL varsetf (VR,WT_NUDGEC(NV),NV,MAXGRDS,FF,0.,1.e20)
 IF(VR.EQ.'IF_ODA')      CALL varseti (VR,IF_ODA,NV,1,II,0,1)
 IF(VR.EQ.'ODA_UPAPREF') CALL varsetc (VR,ODA_UPAPREF,NV,1,CC,1,strl1)
 IF(VR.EQ.'ODA_SFCPREF') CALL varsetc (VR,ODA_SFCPREF,NV,1,CC,1,strl1)
 IF(VR.EQ.'FRQODA')      CALL varsetf (VR,FRQODA,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'TODABEG')     CALL varsetf (VR,TODABEG,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'TODAEND')     CALL varsetf (VR,TODAEND,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'TNUDODA')     CALL varsetf (VR,TNUDODA,NV,1,FF,0.,1.e20)
 IF(VR.EQ.'WT_ODA_GRID') CALL varsetf (VR,WT_ODA_GRID(NV),NV,MAXGRDS,FF,0.,10.)
 IF(VR.EQ.'WT_ODA_UV')   CALL varsetf (VR,WT_ODA_UV,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_ODA_TH')   CALL varsetf (VR,WT_ODA_TH,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_ODA_PI')   CALL varsetf (VR,WT_ODA_PI,NV,1,FF,0.,10.0)
 IF(VR.EQ.'WT_ODA_RT')   CALL varsetf (VR,WT_ODA_RT,NV,1,FF,0.,10.0)
 IF(VR.EQ.'RODA_SFCE')   CALL varsetf (VR,RODA_SFCE(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'RODA_SFC0')   CALL varsetf (VR,RODA_SFC0(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'RODA_UPAE')   CALL varsetf (VR,RODA_UPAE(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'RODA_UPA0')   CALL varsetf (VR,RODA_UPA0(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'RODA_HGT')    CALL varsetf (VR,RODA_HGT(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'RODA_ZFAC')   CALL varsetf (VR,RODA_ZFAC(NV),NV,MAXGRDS,FF,1.,1.e20)
 IF(VR.EQ.'ODA_SFC_TIL') CALL varsetf (VR,ODA_SFC_TIL,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'ODA_SFC_TEL') CALL varsetf (VR,ODA_SFC_TEL,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'ODA_UPA_TIL') CALL varsetf (VR,ODA_UPA_TIL,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'ODA_UPA_TEL') CALL varsetf (VR,ODA_UPA_TEL,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'HFILIN')      CALL varsetc (VR,HFILIN,NV,1,CC,1,strl1)
 IF(VR.EQ.'IPAST_SFC')   CALL varseti (VR,IPAST_SFC,NV,1,II,0,1)
 IF(VR.EQ.'ICLOBBER')    CALL varseti (VR,ICLOBBER,NV,1,II,0,1)
 IF(VR.EQ.'IOUTPUT')     CALL varseti (VR,IOUTPUT,NV,1,II,0,1)
 IF(VR.EQ.'AFILEPREF')   CALL varsetc (VR,AFILEPREF,NV,1,CC,1,strl1)
 IF(VR.EQ.'FRQSTATE')    CALL varsetf (VR,FRQSTATE(NV),NV,MAXGRDS,FF,0.,1.E20)
 IF(VR.EQ.'FRQST_KEEP')  CALL varsetf (VR,FRQST_KEEP,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'FRQLITE')     CALL varsetf (VR,FRQLITE,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'NLITE_VARS')  CALL varseti (VR,NLITE_VARS,NV,1,II,0,MAXLITE)
 IF(VR.EQ.'LITE_VARS')   CALL varsetc (VR,LITE_VARS(NV),NV,MAXLITE,CC,1,32)
 IF(VR.EQ.'AVGTIM')      CALL varsetf (VR,AVGTIM,NV,1,FF,-1.E20,1.E20)
 IF(VR.EQ.'FRQMEAN')     CALL varsetf (VR,FRQMEAN,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'FRQBOTH')     CALL varsetf (VR,FRQBOTH,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'TOPFILES')    CALL varsetc (VR,TOPFILES,NV,1,CC,1,strl1)
 IF(VR.EQ.'SFCFILES')    CALL varsetc (VR,SFCFILES,NV,1,CC,1,strl1)
 IF(VR.EQ.'SSTFPFX')     CALL varsetc (VR,SSTFPFX,NV,1,CC,1,strl1)
 IF(VR.EQ.'NDVIFPFX')    CALL varsetc (VR,NDVIFPFX,NV,1,CC,1,strl1)
 IF(VR.EQ.'ITOPTFLG')    CALL varseti (VR,ITOPTFLG(NV),NV,MAXGRDS,II,0,2)
 IF(VR.EQ.'ISSTFLG')     CALL varseti (VR,ISSTFLG(NV),NV,MAXGRDS,II,0,2)
 IF(VR.EQ.'IVEGTFLG')    CALL varseti (VR,IVEGTFLG(NV),NV,MAXGRDS,II,0,2)
 IF(VR.EQ.'ISOILFLG')    CALL varseti (VR,ISOILFLG(NV),NV,MAXGRDS,II,0,2)
 IF(VR.EQ.'NDVIFLG')     CALL varseti (VR,NDVIFLG(NV),NV,MAXGRDS,II,0,3)
 IF(VR.EQ.'IUPDNDVI')    CALL varseti (VR,IUPDNDVI,NV,1,II,0,2)
 IF(VR.EQ.'IUPDSST')     CALL varseti (VR,IUPDSST,NV,1,II,0,2)
 IF(VR.EQ.'ITOPTFN')     CALL varsetc (VR,ITOPTFN(NV),NV,MAXGRDS,CC,1,strl1)
 IF(VR.EQ.'ISSTFN')      CALL varsetc (VR,ISSTFN(NV),NV,MAXGRDS,CC,1,strl1)
 IF(VR.EQ.'IVEGTFN')     CALL varsetc (VR,IVEGTFN(NV),NV,MAXGRDS,CC,1,strl1)
 IF(VR.EQ.'ISOILFN')     CALL varsetc (VR,ISOILFN(NV),NV,MAXGRDS,CC,1,strl1)
 IF(VR.EQ.'NDVIFN')      CALL varsetc (VR,NDVIFN(NV),NV,MAXGRDS,CC,1,strl1)
 IF(VR.EQ.'ITOPSFLG')    CALL varseti (VR,ITOPSFLG(NV),NV,MAXGRDS,II,0,3)
 IF(VR.EQ.'TOPTENH')     CALL varsetf (VR,TOPTENH(NV),NV,MAXGRDS,FF,0.,5.)
 IF(VR.EQ.'TOPTWVL')     CALL varsetf (VR,TOPTWVL(NV),NV,MAXGRDS,FF,2.,9.)
 IF(VR.EQ.'IZ0FLG')      CALL varseti (VR,IZ0FLG(NV),NV,MAXGRDS,II,0,1)
 IF(VR.EQ.'Z0MAX')       CALL varsetf (VR,Z0MAX(NV),NV,MAXGRDS,FF,0.,20.)
 IF(VR.EQ.'Z0FACT')      CALL varsetf (VR,Z0FACT,NV,1,FF,.0001,.1)
ENDIF

IF(GROUP.EQ.'$MODEL_OPTIONS') THEN
 CALL varchk (VR,GROUP,INDAT,IINDAT,NVINDAT,INRFLG)
 IF(INRFLG.EQ.1) RETURN
 IF(VR.EQ.'ICORFLG')      CALL varseti (VR,ICORFLG,NV,1,II,0,1)
 IF(VR.EQ.'IBND')         CALL varseti (VR,IBND,NV,1,II,1,2)
 IF(VR.EQ.'JBND')         CALL varseti (VR,JBND,NV,1,II,1,2)
 IF(VR.EQ.'CPHAS')        CALL varsetf (VR,CPHAS,NV,1,FF,.001,1.E8)
 IF(VR.EQ.'LSFLG')        CALL varseti (VR,LSFLG,NV,1,II,0,3)
 IF(VR.EQ.'NFPT')         CALL varseti (VR,NFPT,NV,1,II,0,10000)
 IF(VR.EQ.'DISTIM')       CALL varsetf (VR,DISTIM,NV,1,FF,0.,10000.)
 IF(VR.EQ.'ISWRTYP')      CALL varseti (VR,ISWRTYP,NV,1,II,0,3)
 IF(VR.EQ.'ILWRTYP')      CALL varseti (VR,ILWRTYP,NV,1,II,0,3)
 IF(VR.EQ.'RADFRQ')       CALL varsetf (VR,RADFRQ,NV,1,FF,.001,100000.)
 IF(VR.EQ.'LONRAD')       CALL varseti (VR,LONRAD,NV,1,II,0,1)
 IF(VR.EQ.'NNQPARM')      CALL varseti (VR,NNQPARM(NV),NV,MAXGRDS,II,0,2)
 IF(VR.EQ.'CONFRQ')       CALL varsetf (VR,CONFRQ,NV,1,FF,.001,100000.)
 IF(VR.EQ.'WCLDBS')       CALL varsetf (VR,WCLDBS,NV,1,FF,-10.,10.)
 IF(VR.EQ.'IKPP')         CALL varseti (VR,IKPP,NV,1,II,0,2)
 IF(VR.EQ.'NKPPZ')        CALL varseti (VR,NKPPZ,NV,1,II,5,200)
 IF(VR.EQ.'FRQKPP')       CALL varsetf (VR,FRQKPP,NV,1,FF,.001,3600.)
 IF(VR.EQ.'RELAX_SST')    CALL varsetf (VR,RELAX_SST,NV,1,FF,0.,360.)
 IF(VR.EQ.'RELAX_OCNT')   CALL varsetf (VR,RELAX_OCNT,NV,1,FF,0.,360.)
 IF(VR.EQ.'RELAX_SAL')    CALL varsetf (VR,RELAX_SAL,NV,1,FF,0.,360.)
 IF(VR.EQ.'DMAXKPP')      CALL varsetf (VR,DMAXKPP,NV,1,FF,100.,3000.)
 IF(VR.EQ.'DSCALEKPP')    CALL varsetf (VR,DSCALEKPP,NV,1,FF,1.,20.)
 IF(VR.EQ.'KPPITERMAX')   CALL varseti (VR,KPPITERMAX,NV,1,II,1,1000)
 IF(VR.EQ.'KPPRNT')       CALL varseti (VR,KPPRNT,NV,1,II,0,5)
 IF(VR.EQ.'UBMN_KPP')     CALL varsetf (VR,UBMN_KPP,NV,1,FF,0.,7.)
 IF(VR.EQ.'NPATCH')       CALL varseti (VR,NPATCH,NV,1,II,2,10)
 IF(VR.EQ.'NVEGPAT')      CALL varseti (VR,NVEGPAT,NV,1,II,1,10)
 IF(VR.EQ.'ISFCL')        CALL varseti (VR,ISFCL,NV,1,II,0,2)
 IF(VR.EQ.'IFREESLIP')    CALL varseti (VR,IFREESLIP,NV,1,II,0,1)
 IF(VR.EQ.'SIBFILE')      CALL varsetc (VR,SIBFILE,NV,1,CC,1,strl1)
 IF(VR.EQ.'CO2_INIT')     CALL varsetf (VR,CO2_INIT(NV),NV,NZPMAX,FF,0.,1000.)
 IF(VR.EQ.'ISOILDAT')     CALL varseti (VR,ISOILDAT,NV,1,II,0,1)
 IF(VR.EQ.'SNUDCENT')     CALL varsetf (VR,SNUDCENT,NV,1,FF,0.,1.e10)
 IF(VR.EQ.'ISNOWDAT')     CALL varseti (VR,ISNOWDAT,NV,1,II,0,1)
 IF(VR.EQ.'NVGCON')       CALL varseti (VR,NVGCON,NV,1,II,0,20)
 IF(VR.EQ.'PCTLCON')      CALL varsetf (VR,PCTLCON,NV,1,FF,0.,1.)
 IF(VR.EQ.'NSLCON')       CALL varseti (VR,NSLCON,NV,1,II,1,12)
 IF(VR.EQ.'ZROUGH')       CALL varsetf (VR,ZROUGH,NV,1,FF,.0001,100.)
 IF(VR.EQ.'ALBEDO')       CALL varsetf (VR,ALBEDO,NV,1,FF,0.,1.)
 IF(VR.EQ.'SEATMP')       CALL varsetf (VR,SEATMP,NV,1,FF,100.,500.)
 IF(VR.EQ.'DTHCON')       CALL varsetf (VR,DTHCON,NV,1,FF,-100.,100.)
 IF(VR.EQ.'DRTCON')       CALL varsetf (VR,DRTCON,NV,1,FF,-1.,1.)
 IF(VR.EQ.'SLZ')          CALL varsetf (VR,SLZ(NV),NV,NZGMAX,FF,-1.E5,0.)
 IF(VR.EQ.'SLMSTR')       CALL varsetf (VR,SLMSTR(NV),NV,NZGMAX,FF,0.,1.)
 IF(VR.EQ.'STGOFF')       CALL varsetf (VR,STGOFF(NV),NV,NZGMAX,FF,-50.,50.)
 IF(VR.EQ.'IDIFFK')       CALL varseti (VR,IDIFFK(NV),NV,MAXGRDS,II,1,4)
 IF(VR.EQ.'IDIFFPERTS')   CALL varseti (VR,IDIFFPERTS,NV,1,II,0,3)
 IF(VR.EQ.'IHORGRAD')     CALL varseti (VR,IHORGRAD,NV,1,II,1,2)
 IF(VR.EQ.'CSX')          CALL varsetf (VR,CSX(NV),NV,MAXGRDS,FF,0.,10.)
 IF(VR.EQ.'CSZ')          CALL varsetf (VR,CSZ(NV),NV,MAXGRDS,FF,0.,10.)
 IF(VR.EQ.'XKHKM')        CALL varsetf (VR,XKHKM(NV),NV,MAXGRDS,FF,0.,100.)
 IF(VR.EQ.'ZKHKM')        CALL varsetf (VR,ZKHKM(NV),NV,MAXGRDS,FF,0.,100.)
 IF(VR.EQ.'AKMIN')        CALL varsetf (VR,AKMIN(NV),NV,MAXGRDS,FF,0.,5.)
 IF(VR.EQ.'FRACSAT')      CALL varsetf (VR,FRACSAT,NV,1,FF,0.000,999.000)
 IF(VR.EQ.'IBUBBLE')      CALL varseti (VR,IBUBBLE,NV,1,II,0,3)
 IF(VR.EQ.'IBUBGRD')      CALL varseti (VR,IBUBGRD,NV,1,II,1,10)
 IF(VR.EQ.'IBDXIA')       CALL varseti (VR,IBDXIA,NV,1,II,1,3000)
 IF(VR.EQ.'IBDXIZ')       CALL varseti (VR,IBDXIZ,NV,1,II,1,3000)
 IF(VR.EQ.'IBDYJA')       CALL varseti (VR,IBDYJA,NV,1,II,1,3000)
 IF(VR.EQ.'IBDYJZ')       CALL varseti (VR,IBDYJZ,NV,1,II,1,3000)
 IF(VR.EQ.'IBDZK1')       CALL varseti (VR,IBDZK1,NV,1,II,1,300)
 IF(VR.EQ.'IBDZK2')       CALL varseti (VR,IBDZK2,NV,1,II,1,300)
 IF(VR.EQ.'BTHP')         CALL varsetf (VR,BTHP,NV,1,FF,-20.,20.)
 IF(VR.EQ.'BRTP')         CALL varsetf (VR,BRTP,NV,1,FF,-1.,10.)
 IF(VR.EQ.'ICONV')        CALL varseti (VR,ICONV,NV,1,II,0,5)
 IF(VR.EQ.'ICONGR')       CALL varseti (VR,ICONGR,NV,1,II,0,10)
 IF(VR.EQ.'ICICENT')      CALL varseti (VR,ICICENT,NV,1,II,1,3000)
 IF(VR.EQ.'ICJCENT')      CALL varseti (VR,ICJCENT,NV,1,II,1,3000)
 IF(VR.EQ.'CXRAD')        CALL varsetf (VR,CXRAD,NV,1,FF,0.,100000.)
 IF(VR.EQ.'CYRAD')        CALL varsetf (VR,CYRAD,NV,1,FF,0.,100000.)
 IF(VR.EQ.'ICVERT')       CALL varseti (VR,ICVERT,NV,1,II,1,2)
 IF(VR.EQ.'ICKMAX')       CALL varseti (VR,ICKMAX,NV,1,II,1,300)
 IF(VR.EQ.'CZRAD')        CALL varsetf (VR,CZRAD,NV,1,FF,0.,10000.)
 IF(VR.EQ.'ICKCENT')      CALL varseti (VR,ICKCENT,NV,1,II,1,3000)
 IF(VR.EQ.'CDIVMAX')      CALL varsetf (VR,CDIVMAX,NV,1,FF,-1.,1.)
 IF(VR.EQ.'CTAU')         CALL varsetf (VR,CTAU,NV,1,FF,1.,86400.)
 IF(VR.EQ.'CTMAX')        CALL varsetf (VR,CTMAX,NV,1,FF,-200.,86400.)
 IF(VR.EQ.'IRCE')         CALL varseti (VR,IRCE,NV,1,II,0,1)
 IF(VR.EQ.'RCE_SZEN')     CALL varsetf (VR,RCE_SZEN,NV,1,FF,0.,90.)
 IF(VR.EQ.'RCE_SOLC')     CALL varsetf (VR,RCE_SOLC,NV,1,FF,0.,3000.)
 IF(VR.EQ.'RCE_UBMN')     CALL varsetf (VR,RCE_UBMN,NV,1,FF,0.,7.)
 IF(VR.EQ.'RCE_BUBL')     CALL varsetf (VR,RCE_BUBL,NV,1,FF,0.,5.)
 IF(VR.EQ.'ITRACER')      CALL varseti (VR,ITRACER,NV,1,II,0,100)
 IF(VR.EQ.'ITRACHIST')    CALL varseti (VR,ITRACHIST,NV,1,II,0,1)
 IF(VR.EQ.'LEVEL')        CALL varseti (VR,LEVEL,NV,1,II,0,4)
 IF(VR.EQ.'ISCM')         CALL varseti (VR,ISCM,NV,1,II,0,10)
 IF(VR.EQ.'SCMTIME')      CALL varsetf (VR,SCMTIME,NV,1,FF,0.,1.E20)
 IF(VR.EQ.'ISCMX')        CALL varseti (VR,ISCMX,NV,1,II,0,90000)
 IF(VR.EQ.'ISCMY')        CALL varseti (VR,ISCMY,NV,1,II,0,90000)
 IF(VR.EQ.'ICHECKMIC')    CALL varseti (VR,ICHECKMIC,NV,1,II,0,2)
 IF(VR.EQ.'IMBUDGET')     CALL varseti (VR,IMBUDGET,NV,1,II,0,3)
 IF(VR.EQ.'IRIME')        CALL varseti (VR,IRIME,NV,1,II,0,1)
 IF(VR.EQ.'IPLAWS')       CALL varseti (VR,IPLAWS,NV,1,II,0,2)
 IF(VR.EQ.'ISEDIM')       CALL varseti (VR,ISEDIM,NV,1,II,0,1)
 IF(VR.EQ.'ICLOUD')       CALL varseti (VR,ICLOUD,NV,1,II,0,5)
 IF(VR.EQ.'IDRIZ')        CALL varseti (VR,IDRIZ,NV,1,II,0,5)
 IF(VR.EQ.'IRAIN')        CALL varseti (VR,IRAIN,NV,1,II,0,5)
 IF(VR.EQ.'IPRIS')        CALL varseti (VR,IPRIS,NV,1,II,0,7)
 IF(VR.EQ.'ISNOW')        CALL varseti (VR,ISNOW,NV,1,II,0,5)
 IF(VR.EQ.'IAGGR')        CALL varseti (VR,IAGGR,NV,1,II,0,5)
 IF(VR.EQ.'IGRAUP')       CALL varseti (VR,IGRAUP,NV,1,II,0,5)
 IF(VR.EQ.'IHAIL')        CALL varseti (VR,IHAIL,NV,1,II,0,5)
 IF(VR.EQ.'CPARM')        CALL varsetf (VR,CPARM,NV,1,FF,0.,1.E10)
 IF(VR.EQ.'DPARM')        CALL varsetf (VR,DPARM,NV,1,FF,0.,1.E10)
 IF(VR.EQ.'RPARM')        CALL varsetf (VR,RPARM,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'PPARM')        CALL varsetf (VR,PPARM,NV,1,FF,0.,1.E8)
 IF(VR.EQ.'SPARM')        CALL varsetf (VR,SPARM,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'APARM')        CALL varsetf (VR,APARM,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'GPARM')        CALL varsetf (VR,GPARM,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'HPARM')        CALL varsetf (VR,HPARM,NV,1,FF,0.,1.E5)
 IF(VR.EQ.'GNU')          CALL varsetf (VR,GNU(NV),NV,8,FF,0.,20.)
 IF(VR.EQ.'HUCMFILE')     CALL varsetc (VR,HUCMFILE,NV,1,CC,1,strl1)
 IF(VR.EQ.'NDTCOLL')      CALL varseti (VR,NDTCOLL,NV,1,II,1,10)
 IF(VR.EQ.'IAEROSOL')     CALL varseti (VR,IAEROSOL,NV,1,II,0,1)
 IF(VR.EQ.'IABCARB')      CALL varseti (VR,IABCARB,NV,1,II,0,1)
 IF(VR.EQ.'ISALT')        CALL varseti (VR,ISALT,NV,1,II,0,2)
 IF(VR.EQ.'IDUST')        CALL varseti (VR,IDUST,NV,1,II,0,2)
 IF(VR.EQ.'IDUSTLOFT')    CALL varseti (VR,IDUSTLOFT,NV,1,II,0,99)
 IF(VR.EQ.'DUSTFILE')     CALL varsetc (VR,DUSTFILE,NV,1,CC,1,strl1)
 IF(VR.EQ.'ICCNLEV')      CALL varseti (VR,ICCNLEV,NV,1,II,0,3)
 IF(VR.EQ.'IIFN')         CALL varseti (VR,IIFN,NV,1,II,0,3)
 IF(VR.EQ.'IIFN_FORMULA') CALL varseti (VR,IIFN_FORMULA,NV,1,II,1,2)
 IF(VR.EQ.'IAERORAD')     CALL varseti (VR,IAERORAD,NV,1,II,0,1)
 IF(VR.EQ.'IAERODEP')     CALL varseti (VR,IAERODEP,NV,1,II,0,1)
 IF(VR.EQ.'IAEROPRNT')    CALL varseti (VR,IAEROPRNT,NV,1,II,0,1)
 IF(VR.EQ.'IAEROHIST')    CALL varseti (VR,IAEROHIST,NV,1,II,0,1)
 IF(VR.EQ.'CIN_MAX')      CALL varsetf (VR,CIN_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'CCN_MAX')      CALL varsetf (VR,CCN_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'GCCN_MAX')     CALL varsetf (VR,GCCN_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'DUST1_MAX')    CALL varsetf (VR,DUST1_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'DUST2_MAX')    CALL varsetf (VR,DUST2_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'ABC1_MAX')     CALL varsetf (VR,ABC1_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'ABC2_MAX')     CALL varsetf (VR,ABC2_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'SALTF_MAX')    CALL varsetf (VR,SALTF_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'SALTJ_MAX')    CALL varsetf (VR,SALTJ_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'SALTS_MAX')    CALL varsetf (VR,SALTS_MAX,NV,1,FF,0.,1.E4)
 IF(VR.EQ.'IAEROLBC')     CALL varseti (VR,IAEROLBC(NV),NV,MAXGRDS,II,0,1)
 IF(VR.EQ.'ICO2LBC')      CALL varseti (VR,ICO2LBC(NV),NV,MAXGRDS,II,0,1)
 IF(VR.EQ.'BCTAU')        CALL varsetf (VR,BCTAU(NV),NV,MAXGRDS,FF,0.,1.E4)
 IF(VR.EQ.'IAERO_CHEM')   CALL varseti (VR,IAERO_CHEM(NV),NV,aerocat,II,0,2)
 IF(VR.EQ.'AERO_EPSILON') CALL varsetf (VR,AERO_EPSILON(NV),NV,aerocat,FF,0.,1.)
 IF(VR.EQ.'AERO_MEDRAD')  CALL varsetf (VR,AERO_MEDRAD(NV),NV,aerocat,FF,0.,20.)
 IF(VR.EQ.'ITRKEPSILON')  CALL varseti (VR,ITRKEPSILON,NV,1,II,0,1)
 IF(VR.EQ.'ITRKDUST')     CALL varseti (VR,ITRKDUST,NV,1,II,0,1)
 IF(VR.EQ.'ITRKDUSTIFN')  CALL varseti (VR,ITRKDUSTIFN,NV,1,II,0,1)
ENDIF

IF(GROUP.EQ.'$MODEL_SOUND') THEN
 CALL varchk (VR,GROUP,SOUND,ISOUND,NVSOUND,INRFLG)
 IF(INRFLG.EQ.1) RETURN
 IF(VR.EQ.'IPSFLG')  CALL varseti (VR,IPSFLG,NV,1,II,0,10)
 IF(VR.EQ.'ITSFLG')  CALL varseti (VR,ITSFLG,NV,1,II,0,10)
 IF(VR.EQ.'IRTSFLG') CALL varseti (VR,IRTSFLG,NV,1,II,0,10)
 IF(VR.EQ.'IUSFLG')  CALL varseti (VR,IUSFLG,NV,1,II,0,10)
 IF(VR.EQ.'HS')      CALL varsetf (VR,HS(NV),NV,MAXSNDG,FF,-1.e3,1.E5)
 IF(VR.EQ.'PS')      CALL varsetf (VR,PS(NV),NV,MAXSNDG,FF,0.,1.E7)
 IF(VR.EQ.'TS')      CALL varsetf (VR,TS(NV),NV,MAXSNDG,FF,-200.,1000.)
 IF(VR.EQ.'RTS')     CALL varsetf (VR,RTS(NV),NV,MAXSNDG,FF,-200.,1000.)
 IF(VR.EQ.'US')      CALL varsetf (VR,US(NV),NV,MAXSNDG,FF,-500.,500.)
 IF(VR.EQ.'VS')      CALL varsetf (VR,VS(NV),NV,MAXSNDG,FF,-500.,500.)
ENDIF

return
END SUBROUTINE nvfillm

!##############################################################################
Subroutine nameout ()

use mem_all
use isan_coms
use leaf_coms, only:ifreeslip

implicit none

integer :: ng,k,m

! This routine prints out a listing of the values of all variables
! in the NAMELISTS

WRITE(6,100)
100  FORMAT(/,'----------------------------NAMELIST VARIABLES-------'  &
      ,'------------------------',/)

WRITE(6,*)'Grid-dependent Integers:'
WRITE(6,101)(' ',NNXP(NG),NNYP(NG),NNZP(NG)  &
            ,NXTNEST(NG),NSTRATX(NG),NG=1,NGRIDS)
WRITE(6,102)(' ',NNDTRAT(NG),NINEST(NG) &
            ,NJNEST(NG),NKNEST(NG),NG=1,NGRIDS)
WRITE(6,103)(' ',NNSTTOP(NG),NNSTBOT(NG),ITOPTFLG(NG) &
            ,ISSTFLG(NG),IVEGTFLG(NG),NG=1,NGRIDS)
WRITE(6,104)(' ',ISOILFLG(NG),NDVIFLG(NG)  &
            ,NNQPARM(NG),IDIFFK(NG),NG=1,NGRIDS)
WRITE(6,105)(' ',IAEROLBC(NG),ICO2LBC(NG)  &
            ,NG=1,NGRIDS)

101  FORMAT(A1,'    NNXP=',I4,'       NNYP=',I4,'       NNZP=',I4  &
     ,'    NXTNEST=',I4,'    NSTRATX=',I4,999(A1,/,I14,4I16))
102  FORMAT(A1,' NNDTRAT=',I4,'     NINEST=',I4  &
     ,'     NJNEST=',I4,'     NKNEST=',I4,999(A1,/,I14,3I16))
103  FORMAT(A1,' NNSTTOP=',I4,'    NNSTBOT=',I4,'   ITOPTFLG=',I4  &
     ,'    ISSTFLG=',I4,'   IVEGTFLG=',I4,999(A1,/,I14,4I16))
104  FORMAT(A1,'ISOILFLG=',I4,'    NDVIFLG=',I4  &
     ,'    NNQPARM=',I4,'     IDIFFK=',I4,999(A1,/,I14,3I16))
105  FORMAT(A1,'IAEROLBC=',I4,'    ICO2LBC=',I4  &
     ,999(A1,/,I14,1I16))

PRINT*, ' '

WRITE(6,*)'Non-grid-dependent Integers:'
WRITE(6,'(100(3(A19,I5)/))')         &
 ,'IMONTH1=',IMONTH1                 &
 ,'IDATE1=',IDATE1                   &
 ,'IYEAR1=',IYEAR1                   &
 ,'ITIME1=',ITIME1                   &
 ,'NGRIDS=',NGRIDS                   &
 ,'NZG=',NZG                         &
 ,'NZS=',NZS                         &
 ,'IPRNTSTMT=',IPRNTSTMT             &
 ,'IHTRAN=',IHTRAN                   &
 ,'NACOUST=',NACOUST                 &
 ,'NESTZ=',NESTZ                     &
 ,'INITIAL=',INITIAL                 &
 ,'NUD_TYPE=',NUD_TYPE               &
 ,'NUDLAT=',NUDLAT                   &
 ,'NUD_COND=',NUD_COND               &
 ,'IF_ODA=',IF_ODA                   &
 ,'IPAST_SFC=',IPAST_SFC             &
 ,'ICLOBBER=',ICLOBBER               &
 ,'IOUTPUT=',IOUTPUT                 &
 ,'NLITE_VARS=',NLITE_VARS           &
 ,'IUPDNDVI=',IUPDNDVI               &
 ,'IUPDSST=',IUPDSST                 &
 ,'ICORFLG=',ICORFLG                 &
 ,'IBND=',IBND                       &
 ,'JBND=',JBND                       &
 ,'LSFLG=',LSFLG                     &
 ,'NFPT=',NFPT                       &
 ,'ISWRTYP=',ISWRTYP                 &
 ,'ILWRTYP=',ILWRTYP                 &
 ,'LONRAD=',LONRAD                   &
 ,'NPATCH=',NPATCH                   &
 ,'NVEGPAT=',NVEGPAT                 &
 ,'ISFCL=',ISFCL                     &
 ,'IFREESLIP=',IFREESLIP             &
 ,'IKPP=',IKPP                       &
 ,'NKPPZ=',NKPPZ                     &
 ,'KPPITERMAX=',KPPITERMAX           &
 ,'KPPRNT=',KPPRNT                   &
 ,'ISOILDAT=',ISOILDAT               &
 ,'ISNOWDAT=',ISNOWDAT               &
 ,'NVGCON=',NVGCON                   &
 ,'NSLCON=',NSLCON                   &
 ,'IDIFFPERTS=',IDIFFPERTS           &
 ,'IHORGRAD=',IHORGRAD               &
 ,'IBUBBLE=',IBUBBLE                 &
 ,'IBUBGRD=',IBUBGRD                 &
 ,'IBDXIA=',IBDXIA                   &
 ,'IBDXIZ=',IBDXIZ                   &
 ,'IBDYJA=',IBDYJA                   &
 ,'IBDYJZ=',IBDYJZ                   &
 ,'IBDZK1=',IBDZK1                   &
 ,'IBDZK2=',IBDZK2                   &
 ,'ICONV=',ICONV                     &
 ,'ICONGR=',ICONGR                   &
 ,'ICICENT=',ICICENT                 &
 ,'ICJCENT=',ICJCENT                 &
 ,'ICVERT=',ICVERT                   &
 ,'ICKMAX=',ICKMAX                   &
 ,'ICKCENT=',ICKCENT                 &
 ,'IRCE=',IRCE                       &
 ,'ITRACER=',ITRACER                 &
 ,'ITRACHIST=',ITRACHIST             &
 ,'LEVEL=',LEVEL                     &
 ,'ISCM=',ISCM                       &
 ,'ISCMX=',ISCMX                     &
 ,'ISCMY=',ISCMY                     &
 ,'ICHECKMIC=',ICHECKMIC             &
 ,'IMBUDGET=',IMBUDGET               &
 ,'IRIME=',IRIME                     &
 ,'IPLAWS=',IPLAWS                   &
 ,'ISEDIM=',ISEDIM                   &
 ,'ICLOUD=',ICLOUD                   &
 ,'IDRIZ=',IDRIZ                     &
 ,'IRAIN=',IRAIN                     &
 ,'IPRIS=',IPRIS                     &
 ,'ISNOW=',ISNOW                     &
 ,'IAGGR=',IAGGR                     &
 ,'IGRAUP=',IGRAUP                   &
 ,'IHAIL=',IHAIL                     &
 ,'NDTCOLL=',NDTCOLL                 &
 ,'IAEROSOL=',IAEROSOL               &
 ,'IABCARB=',IABCARB                 &
 ,'ISALT=',ISALT                     &
 ,'IDUST=',IDUST                     &
 ,'IDUSTLOFT=',IDUSTLOFT             &
 ,'ICCNLEV=',ICCNLEV                 &
 ,'IIFN=',IIFN                       &
 ,'IIFN_FORMULA=',IIFN_FORMULA       &
 ,'IAERORAD=',IAERORAD               &
 ,'IAERODEP=',IAERODEP               &
 ,'IAEROPRNT=',IAEROPRNT             &
 ,'IAEROHIST=',IAEROHIST             &
 ,'ITRKEPSILON=',ITRKEPSILON         &
 ,'ITRKDUST=',ITRKDUST               &
 ,'ITRKDUSTIFN=',ITRKDUSTIFN         &
 ,'IPSFLG=',IPSFLG                   &
 ,'ITSFLG=',ITSFLG                   &
 ,'IRTSFLG=',IRTSFLG                 &
 ,'IUSFLG=',IUSFLG                   &
 ,'IMPL=',IMPL

PRINT*, ' '
WRITE(6,*)'Non-grid-dependent Floats:'
WRITE(6,'(100(3(A15,E11.4)/))')      &
 ,'TIMMAX=',TIMMAX                   &
 ,'DELTAX=',DELTAX                   &
 ,'DELTAZ=',DELTAZ                   &
 ,'DZRAT=',DZRAT                     &
 ,'DZMAX=',DZMAX                     &
 ,'DTLONG=',DTLONG                   &
 ,'POLELAT=',POLELAT                 &
 ,'POLELON=',POLELON                 &
 ,'TNUDLAT=',TNUDLAT                 &
 ,'TNUDCENT=',TNUDCENT               &
 ,'SNUDCENT=',SNUDCENT               &
 ,'TNUDTOP=',TNUDTOP                 &
 ,'ZNUDTOP=',ZNUDTOP                 &
 ,'FRQSTATE=',frqstate(1)            &
 ,'Z0FACT=',Z0FACT                   &
 ,'CPHAS=',CPHAS                     &
 ,'DISTIM=',DISTIM                   &
 ,'RADFRQ=',RADFRQ                   &
 ,'CONFRQ=',CONFRQ                   &
 ,'FRQKPP=',FRQKPP                   &
 ,'RELAX_SST=',RELAX_SST             &
 ,'RELAX_OCNT=',RELAX_OCNT           &
 ,'RELAX_SAL=',RELAX_SAL             &
 ,'DMAXKPP=',DMAXKPP                 &
 ,'DSCALEKPP=',DSCALEKPP             &
 ,'UBMN_KPP=',UBMN_KPP               &
 ,'WCLDBS=',WCLDBS                   &
 ,'PCTLCON=',PCTLCON                 &
 ,'ZROUGH=',ZROUGH                   &
 ,'ALBEDO=',ALBEDO                   &
 ,'SEATMP=',SEATMP                   &
 ,'DTHCON=',DTHCON                   &
 ,'DRTCON=',DRTCON                   &
 ,'BTHP=',BTHP                       &
 ,'BRTP=',BRTP                       &
 ,'CXRAD=',CXRAD                     &
 ,'CYRAD=',CYRAD                     &
 ,'CZRAD=',CZRAD                     &
 ,'CDIVMAX=',CDIVMAX                 &
 ,'CTAU=',CTAU                       &
 ,'CTMAX=',CTMAX                     &
 ,'RCE_SZEN=',RCE_SZEN               &
 ,'RCE_SOLC=',RCE_SOLC               &
 ,'RCE_UBMN=',RCE_UBMN               &
 ,'RCE_BUBL=',RCE_BUBL               &
 ,'FRACSAT=',FRACSAT                 &
 ,'SCMTIME=',SCMTIME                 &
 ,'CPARM=',CPARM                     &
 ,'DPARM=',DPARM                     &
 ,'RPARM=',RPARM                     &
 ,'PPARM=',PPARM                     &
 ,'SPARM=',SPARM                     &
 ,'APARM=',APARM                     &
 ,'GPARM=',GPARM                     &
 ,'HPARM=',HPARM                     &
 ,'CIN_MAX=',CIN_MAX                 &
 ,'CCN_MAX=',CCN_MAX                 &
 ,'GCCN_MAX=',GCCN_MAX               &
 ,'DUST1_MAX=',DUST1_MAX             &
 ,'DUST2_MAX=',DUST2_MAX             &
 ,'ABC1_MAX=',ABC1_MAX               &
 ,'ABC2_MAX=',ABC2_MAX               &
 ,'SALTF_MAX=',SALTF_MAX             &
 ,'SALTJ_MAX=',SALTJ_MAX             &
 ,'SALTS_MAX=',SALTS_MAX             &
 ,'SSPCT=',SSPCT

PRINT*, ' '
WRITE(6,*)'Grid-dependent Floats:'
WRITE(6,301)(' ',TOPTENH(NG),TOPTWVL(NG),CENTLAT(NG),NG=1,NGRIDS)
WRITE(6,302)(' ',CENTLON(NG),CSX(NG),CSZ(NG),NG=1,NGRIDS)
WRITE(6,303)(' ',XKHKM(NG),ZKHKM(NG),AKMIN(NG),NG=1,NGRIDS)
WRITE(6,304)(' ',BCTAU(NG),NG=1,NGRIDS)
301  FORMAT(A1,'TOPTENH=',E12.5,'        TOPTWVL=',E12.5  &
   ,'        CENTLAT=',E12.5,999(A1,/,E21.5,2E28.5))
302  FORMAT(A1,'CENTLON=',E12.5,'            CSX=',E12.5  &
   ,'            CSZ=',E12.5,999(A1,/,E21.5,2E28.5))
303  FORMAT(A1,'  XKHKM=',E12.5,'          ZKHKM=',E12.5  &
   ,'          AKMIN=',E12.5,999(A1,/,E21.5,2E28.5))
304  FORMAT(A1,'  BCTAU=',E12.5,999(A1,/,E21.5))

PRINT*, ' '

WRITE(6,601)(' ',ITOPTFN(M),M=1,NGRIDS)
WRITE(6,602)(' ',ISSTFN(M),M=1,NGRIDS)
WRITE(6,603)(' ',IVEGTFN(M),M=1,NGRIDS)
WRITE(6,604)(' ',ISOILFN(M),M=1,NGRIDS)
WRITE(6,605)(' ',NDVIFN(M),M=1,NGRIDS)
WRITE(6,607) ' ',trim(VARFPFX)
WRITE(6,608) ' ',trim(SSTFPFX)
WRITE(6,609) ' ',trim(NDVIFPFX)
WRITE(6,610) ' ',trim(DUSTFILE)
WRITE(6,611) ' ',trim(SIBFILE)
WRITE(6,612) ' ',trim(HUCMFILE)
601  FORMAT(A1,'  ITOPTFN=',A40,999(A1,/,11X,A40))
602  FORMAT(A1,'   ISSTFN=',A40,999(A1,/,11X,A40))
603  FORMAT(A1,'  IVEGTFN=',A40,999(A1,/,11X,A40))
604  FORMAT(A1,'  ISOILFN=',A40,999(A1,/,11X,A40))
605  FORMAT(A1,'   NDVIFN=',A40,999(A1,/,11X,A40))
607  FORMAT(A1,'  VARFPFX=',A)
608  FORMAT(A1,'  SSTFPFX=',A)
609  FORMAT(A1,' NDVIFPFX=',A)
610  FORMAT(A1,' DUSTFILE=',A)
611  FORMAT(A1,'  SIBFILE=',A)
612  FORMAT(A1,' HUCMFILE=',A)

PRINT*, ' '
WRITE(6,701)EXPNME
PRINT*, ' '
WRITE(6,702)HFILIN
WRITE(6,704)AFILEPREF
PRINT*, ' '
WRITE(6,705)RUNTYPE,TIMEUNIT

701  FORMAT('  EXPNME=',A40)
702  FORMAT('  HFILIN=',A40)
704  FORMAT(' AFILEPREF=',A40)
705  FORMAT(' RUNTYPE=',A10,'      TIMEUNIT=',A3)

PRINT*, ' '
WRITE(6,901)(' ',SLMSTR(K),SLZ(K),STGOFF(K),K=1,NZG)
901  FORMAT(A1,' SLMSTR=',F6.2,'      SLZ=',F7.3,'      STGOFF=',F8.2  &
     ,999(A1,/,F15.2,F17.3,F21.2))

PRINT*, ' '
WRITE(6,1001)(ZZ(K),K=1,NNZP(1))
1001  FORMAT('ZZ=',8F9.1,/,(F12.1,7F9.1))

PRINT*, ' '
WRITE(6,1002)(CO2_INIT(K),K=1,NNZP(1))
1002  FORMAT('CO2_INIT=',8F9.1,/,(F12.1,7F9.1))

WRITE(6,1101)(NSTRATZ(K),K=1,NNZP(1))
1101 FORMAT(/,'NSTRATZ=',(t9,23i3) )

write(6,1301)(gnu(k),k=1,8)
1301 format(/,'GNU=',(t9,8f5.2))

PRINT*, ' '
WRITE(6,*)'IAERO_CHEM    AERO_EPSILON    AERO_MEDRAD-(default)'
WRITE(6,*)'(category)     (fraction)      (microns)'
WRITE(6,1404)(iaero_chem(K),aero_epsilon(K),aero_medrad(K)*1.e6,K=1,aerocat)
1404  FORMAT(I6,F17.3,F16.2)

PRINT*, ' '

! write out ISAN config if this is a makevfile or makehfile run
if (trim(RUNTYPE) .eq. 'MAKEVFILE' .or. trim(RUNTYPE) .eq. 'MAKEHFILE') then
  WRITE(6,*)'ISAN Integers:'
  WRITE(6,'(100(3(A19,I5)/))')         &
   ,'ISZSTAGE=',ISZSTAGE               &
   ,'IVRSTAGE=',IVRSTAGE               &
   ,'ISAN_INC=',ISAN_INC               &
   ,'I1ST_FLG=',I1ST_FLG               &
   ,'IUPA_FLG=',IUPA_FLG               &
   ,'ISFC_FLG=',ISFC_FLG               &
   ,'IOFLGISZ=',IOFLGISZ               &
   ,'IOFLGVAR=',IOFLGVAR               &
   ,'IDATAIN=',IDATAIN                 &
   ,'NIGRIDS=',NIGRIDS                 &
   ,'NFEEDVAR=',NFEEDVAR               &
   ,'NISN=',NISN                       &
   ,'MAXSTA=',MAXSTA                   &
   ,'MAXSFC=',MAXSFC                   &
   ,'IGRIDFL=',IGRIDFL                 &
   ,'NOTSTA=',NOTSTA                   &
   ,'IOBSWIN=',IOBSWIN
  
  PRINT*, ' '
  WRITE(6,1501)(LEVTH(K),K=1,NISN)
1501  FORMAT('LEVTH=',8I9,/,(I15,7I9))


  PRINT*, ' '
  WRITE(6,*)'ISAN Floats:'
  WRITE(6,'(100(3(A15,E11.4)/))')      &
   ,'SIGZWT',SIGZWT                    &
   ,'TOPSIGZ',TOPSIGZ                  &
   ,'HYBBOT',HYBBOT                    &
   ,'HYBTOP',HYBTOP                    &
   ,'SFCINF',SFCINF                    &
   ,'GOBSEP',GOBSEP                    &
   ,'GOBRAD',GOBRAD                    &
   ,'STASEP',STASEP

  PRINT*, ' '
  WRITE(6,1503)(' ',WVLNTH(NG),SWVLNTH(NG),RESPON(NG),GRIDWT(NG),NG=1,NIGRIDS)
1503  FORMAT(A1,'  WVLNTH=',F8.2,' SWVLNTH=',F8.2,'  RESPON=',F8.2,'  GRIDWT=',F8.2  &
      ,999(A1,/,F18.2,F17.2,F17.2,F17.2))

  PRINT*, ' '
  WRITE(6,*) '       IAPR=', trim(IAPR)
  WRITE(6,*) '     IARAWI=', trim(IARAWI)
  WRITE(6,*) '    IASRFCE=', trim(IASRFCE)
  WRITE(6,*) '  VAR_HFILE=', trim(VAR_HFILE)
  WRITE(6,*) '     VARPFX=', trim(VARPFX)
!  WRITE(6,*) '      NOTID=', trim(NOTID)
  WRITE(6,*) '  USED_FILE=', trim(USED_FILE)

  PRINT*, ' '
endif

return
END SUBROUTINE nameout

!##############################################################################
Subroutine namein (IUNIT,GROUP)

! This routine is called by routines READ_NL
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
         IF(NVARN.GT.0) CALL nvtran (GROUP,VARN,VALUE,NVALUE)
         NVALUE=0
         VARN=TOKENS(NT)
         NVARN=NVARN+1
         NT=NT+1
         IF(TOKENS(NT).EQ.'(' .or. TOKENS(NT).EQ.',') THEN
           stop 'Error in RAMSIN syntax. Routine namein'
         ENDIF
      ELSEIF(numberchk(TOKENS(NT)).EQ.1 .OR. letquo(TOKENS(NT)).EQ.1) THEN
          NVALUE=NVALUE+1
          VALUE(NVALUE)=TOKENS(NT)
      ENDIF
      NT=NT+1
   GOTO 20
10 CONTINUE

100 CONTINUE
CALL nvtran (GROUP,VARN,VALUE,NVALUE)
VARN='$END'
CALL nvfillm (GROUP,VARN,INT,FNUM,LINEW(1:NCW),NR)

return
END SUBROUTINE namein

!##############################################################################
Subroutine nvtran (GROUP,VARN,VALUE,NVALUE)

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
         CALL nvfillm (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
      ELSE
         CALL ch2ch (VALUE(NV),CHVAL,NCW)
         CALL nvfillm (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
      ENDIF
   ENDDO
ELSE
   DO NV=1,NVALUE
      IF(letquo(VALUE(NV)).EQ.0) THEN
         CALL ch2real (VALUE(NV),FNUM)
         CALL nvfillm (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
      ELSE
         CALL ch2ch (VALUE(NV),CHVAL,NCW)
         CALL nvfillm (GROUP,VARN,INT,FNUM,CHVAL(1:NCW),NV)
      ENDIF
   ENDDO
ENDIF

return
END SUBROUTINE nvtran
