!##############################################################################
Subroutine broadcast_config ()

! This routine will broadcast the namelist variables to all processes. It is important
! to keep the order of the par_put_* calls matching the order of the par_get_* calls.
! par_put_* stuffs the variables into a buffer (MPI_Pack) in the order of the calls,
! while par_get_* extracts the variables from the buffer (MPI_Unpack) in the order
! of the calls. Mixing up the order between the puts and the gets will scramble
! the variable values, and there is nothing to prevent this from happinging, nor
! detect that it has happened.
!
! Keep the list of variables being broadcast here in sync with the list of variables
! that are being extracted from the namelist file via read_nl() and set by eng_params().
! The list of variables appears in:
!    nvfillm() - rname.f90
!    nvfills() - isan_name.f90
!    eng_params() - ruser.f90

use mem_all
use isan_coms
use node_mod
use leaf_coms, only:ifreeslip

implicit none

  real, allocatable :: buff(:)
  integer :: nwords, nm

  !Saleeby(2016)
  !Increment memory buffer size here if you add RAMSIN Namelist variables.
  !Add to the appropriate section below as (#-of-them * arraysize).
  nwords = 219 * 1                 & !single values
         +   1 * 8                 & !micro (8-hydromet types for gnu)
         +   3 * aerocat           & !micro (number aerosol species)
         +  44 * maxgrds           & !grid-dependent (max grids)
         +   3 * nzpmax            & !max vertical levels
         +   3 * nzgmax            & !max soil levels
         +   1 * maxisn            & !max isentropic levels
         +   4 * maxagrds          & !max varfile grids
         +   6 * maxsndg           & !max input sounding levels
         +   2 * 32                & !32 character length strings
         +   1 * maxlite * 32      & !lite variables 32 char length strings
         +  17 *       1 * strl1   & !individual input strings
         +   1 *      50 * strl1   & !array of input strings
         +   5 * maxgrds * strl1   & !grid-dependent array of input strings
         + 100                       !extras so we have enough buffer

  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! If you add variables below to "put" and "get", you need to increment "nwords"
  ! above. But only increment 1 for a single put/get combination.

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! $MODEL_GRID namelist group
    CALL par_put_char  (EXPNME,strl1)
    CALL par_put_char  (RUNTYPE,32)
    CALL par_put_char  (TIMEUNIT,32)
    CALL par_put_float (TIMMAX,1)
    CALL par_put_int   (IMONTH1,1)
    CALL par_put_int   (IDATE1,1)
    CALL par_put_int   (IYEAR1,1)
    CALL par_put_int   (ITIME1,1)
    CALL par_put_int   (NGRIDS,1)
    CALL par_put_int   (NNXP,MAXGRDS)
    CALL par_put_int   (NNYP,MAXGRDS)
    CALL par_put_int   (NNZP,MAXGRDS)
    CALL par_put_int   (NZG,1)
    CALL par_put_int   (NZS,1)
    CALL par_put_int   (NXTNEST,MAXGRDS)
    CALL par_put_int   (IPRNTSTMT,1)
    CALL par_put_int   (IHTRAN,1)
    CALL par_put_float (DELTAX,1)
    CALL par_put_float (DELTAZ,1)
    CALL par_put_float (DZRAT,1)
    CALL par_put_float (DZMAX,1)
    CALL par_put_float (ZZ,NZPMAX)
    CALL par_put_float (DTLONG,1)
    CALL par_put_int   (NACOUST,1)
    CALL par_put_int   (NSTRATX,MAXGRDS)
    CALL par_put_int   (NSTRATY,MAXGRDS)
    CALL par_put_int   (NNDTRAT,MAXGRDS)
    CALL par_put_int   (NESTZ,1)
    CALL par_put_int   (NSTRATZ,NZPMAX)
    CALL par_put_float (POLELAT,1)
    CALL par_put_float (POLELON,1)
    CALL par_put_int   (NINEST,MAXGRDS)
    CALL par_put_int   (NJNEST,MAXGRDS)
    CALL par_put_int   (NKNEST,MAXGRDS)
    CALL par_put_float (CENTLAT,MAXGRDS)
    CALL par_put_float (CENTLON,MAXGRDS)
    CALL par_put_int   (NNSTTOP,MAXGRDS)
    CALL par_put_int   (NNSTBOT,MAXGRDS)
    
    ! $MODEL_FILE_INFO namelist group
    CALL par_put_int   (INITIAL,1)
    CALL par_put_int   (NUD_TYPE,1)
    CALL par_put_char  (VARFPFX,strl1)
    CALL par_put_float (VWAIT1,1)
    CALL par_put_float (VWAITTOT,1)
    CALL par_put_int   (NUDLAT,1)
    CALL par_put_float (TNUDLAT,1)
    CALL par_put_float (TNUDCENT,1)
    CALL par_put_float (TNUDTOP,1)
    CALL par_put_float (ZNUDTOP,1)
    CALL par_put_float (WT_NUDGE_G,maxgrds)
    CALL par_put_float (WT_NUDGE_UV,1)
    CALL par_put_float (WT_NUDGE_TH,1)
    CALL par_put_float (WT_NUDGE_PI,1)
    CALL par_put_float (WT_NUDGE_RT,1)
    CALL par_put_int   (NUD_COND,1)
    CALL par_put_float (TCOND_BEG,1)
    CALL par_put_float (TCOND_END,1)
    CALL par_put_float (T_NUDGE_RC,1)
    CALL par_put_float (WT_NUDGEC,maxgrds)
    CALL par_put_float (SNUDCENT,1)
    CALL par_put_int   (IF_ODA,1)
    CALL par_put_char  (ODA_UPAPREF,strl1)
    CALL par_put_char  (ODA_SFCPREF,strl1)
    CALL par_put_float (FRQODA,1)
    CALL par_put_float (TODABEG,1)
    CALL par_put_float (TODAEND,1)
    CALL par_put_float (TNUDODA,1)
    CALL par_put_float (WT_ODA_GRID,MAXGRDS)
    CALL par_put_float (WT_ODA_UV,1)
    CALL par_put_float (WT_ODA_TH,1)
    CALL par_put_float (WT_ODA_PI,1)
    CALL par_put_float (WT_ODA_RT,1)
    CALL par_put_float (RODA_SFCE,MAXGRDS)
    CALL par_put_float (RODA_SFC0,MAXGRDS)
    CALL par_put_float (RODA_UPAE,MAXGRDS)
    CALL par_put_float (RODA_UPA0,MAXGRDS)
    CALL par_put_float (RODA_HGT,MAXGRDS)
    CALL par_put_float (RODA_ZFAC,MAXGRDS)
    CALL par_put_float (ODA_SFC_TIL,1)
    CALL par_put_float (ODA_SFC_TEL,1)
    CALL par_put_float (ODA_UPA_TIL,1)
    CALL par_put_float (ODA_UPA_TEL,1)
    CALL par_put_char  (HFILIN,strl1)
    CALL par_put_int   (IPAST_SFC,1)
    CALL par_put_int   (ICLOBBER,1)
    CALL par_put_int   (IOUTPUT,1)
    CALL par_put_char  (AFILEPREF,strl1)
    CALL par_put_float (FRQSTATE,MAXGRDS)
    CALL par_put_float (FRQST_KEEP,1)
    CALL par_put_float (FRQLITE,1)
    CALL par_put_float (NLITE_VARS,1)
    do nm = 1, nlite_vars
       print*,'lite pack:',nm,trim(LITE_VARS(nm))
       CALL par_put_char (LITE_VARS(nm),32)
    enddo
    CALL par_put_float (AVGTIM,1)
    CALL par_put_float (FRQMEAN,1)
    CALL par_put_float (FRQBOTH,1)
    CALL par_put_char  (TOPFILES,strl1)
    CALL par_put_char  (SFCFILES,strl1)
    CALL par_put_char  (SSTFPFX,strl1)
    CALL par_put_char  (NDVIFPFX,strl1)
    CALL par_put_int   (ITOPTFLG,MAXGRDS)
    CALL par_put_int   (ISSTFLG,MAXGRDS)
    CALL par_put_int   (IVEGTFLG,MAXGRDS)
    CALL par_put_int   (ISOILFLG,MAXGRDS)
    CALL par_put_int   (NDVIFLG,MAXGRDS)
    CALL par_put_int   (IUPDNDVI,1)
    CALL par_put_int   (IUPDSST,1)
    CALL par_put_char  (ITOPTFN,MAXGRDS*strl1)
    CALL par_put_char  (ISSTFN,MAXGRDS*strl1)
    CALL par_put_char  (IVEGTFN,MAXGRDS*strl1)
    CALL par_put_char  (ISOILFN,MAXGRDS*strl1)
    CALL par_put_char  (NDVIFN,MAXGRDS*strl1)
    CALL par_put_int   (ITOPSFLG,MAXGRDS)
    CALL par_put_float (TOPTENH,MAXGRDS)
    CALL par_put_float (TOPTWVL,MAXGRDS)
    CALL par_put_int   (IZ0FLG,MAXGRDS)
    CALL par_put_float (Z0MAX,MAXGRDS)
    CALL par_put_float (Z0FACT,1)

    ! $MODEL_OPTIONS namelist group
    CALL par_put_int   (ICORFLG,1)
    CALL par_put_int   (IBND,1)
    CALL par_put_int   (JBND,1)
    CALL par_put_int   (ISPONGE_PTS,MAXGRDS)
    CALL par_put_float (SPONGE_TAU,MAXGRDS)
    CALL par_put_float (CPHAS,1)
    CALL par_put_int   (LSFLG,1)
    CALL par_put_int   (NFPT,1)
    CALL par_put_float (DISTIM,1)
    CALL par_put_int   (ILWRTYP,1)
    CALL par_put_int   (ISWRTYP,1)
    CALL par_put_float (RADFRQ,1)
    CALL par_put_int   (LONRAD,1)
    CALL par_put_int   (NNQPARM,MAXGRDS)
    CALL par_put_float (CONFRQ,1)
    CALL par_put_float (WCLDBS,1)
    CALL par_put_int   (IKPP,1)
    CALL par_put_int   (NKPPZ,1)
    CALL par_put_float (FRQKPP,1)
    CALL par_put_float (RELAX_SST,1)
    CALL par_put_float (RELAX_OCNT,1)
    CALL par_put_float (RELAX_SAL,1)
    CALL par_put_float (DMAXKPP,1)
    CALL par_put_float (DSCALEKPP,1)
    CALL par_put_int   (KPPITERMAX,1)
    CALL par_put_int   (KPPRNT,1)
    CALL par_put_float (UBMN_KPP,1)
    CALL par_put_int   (NPATCH,1)
    CALL par_put_int   (NVEGPAT,1)
    CALL par_put_int   (ISFCL,1)
    CALL par_put_int   (IFREESLIP,1)
    CALL par_put_char  (SIBFILE,strl1)
    CALL par_put_float (CO2_INIT,NZPMAX)
    CALL par_put_int   (ISOILDAT,1)
    CALL par_put_int   (ISNOWDAT,1)
    CALL par_put_int   (NVGCON,1)
    CALL par_put_float (PCTLCON,1)
    CALL par_put_int   (NSLCON,1)
    CALL par_put_float (ZROUGH,1)
    CALL par_put_float (ALBEDO,1)
    CALL par_put_float (SEATMP,1)
    CALL par_put_float (DTHCON,1)
    CALL par_put_float (DRTCON,1)
    CALL par_put_float (SLZ,NZGMAX)
    CALL par_put_float (SLMSTR,NZGMAX)
    CALL par_put_float (STGOFF,NZGMAX)
    CALL par_put_int   (IDIFFK,MAXGRDS)
    CALL par_put_int   (IDIFFPERTS,1)
    CALL par_put_int   (IHORGRAD,1)
    CALL par_put_float (CSX,MAXGRDS)
    CALL par_put_float (CSZ,MAXGRDS)
    CALL par_put_float (XKHKM,MAXGRDS)
    CALL par_put_float (ZKHKM,MAXGRDS)
    CALL par_put_float (AKMIN,MAXGRDS)
    CALL par_put_float (FRACSAT,1)
    CALL par_put_int   (IBUBBLE,1)
    CALL par_put_int   (IBUBGRD,1)
    CALL par_put_int   (IBDXIA,1)
    CALL par_put_int   (IBDXIZ,1)
    CALL par_put_int   (IBDYJA,1)
    CALL par_put_int   (IBDYJZ,1)
    CALL par_put_int   (IBDZK1,1)
    CALL par_put_int   (IBDZK2,1)
    CALL par_put_float (BTHP,1)
    CALL par_put_float (BRTP,1)
    CALL par_put_int   (ICONV,1)
    CALL par_put_int   (ICONGR,1)
    CALL par_put_int   (ICICENT,1)
    CALL par_put_int   (ICJCENT,1)
    CALL par_put_float (CXRAD,1)
    CALL par_put_float (CYRAD,1)
    CALL par_put_int   (ICVERT,1)
    CALL par_put_int   (ICKMAX,1)
    CALL par_put_float (CZRAD,1)
    CALL par_put_int   (ICKCENT,1)
    CALL par_put_float (CDIVMAX,1)
    CALL par_put_float (CTAU,1)
    CALL par_put_float (CTMAX,1)
    CALL par_put_int   (IRCE,1)
    CALL par_put_float (RCE_SZEN,1)
    CALL par_put_float (RCE_SOLC,1)
    CALL par_put_float (RCE_UBMN,1)
    CALL par_put_float (RCE_BUBL,1)
    CALL par_put_int   (ITRACER,1)
    CALL par_put_int   (ITRACHIST,1)
    CALL par_put_int   (LEVEL,1)
    CALL par_put_int   (ISCM,1)
    CALL par_put_float (SCMTIME,1)
    CALL par_put_int   (ISCMX,1)
    CALL par_put_int   (ISCMY,1)
    CALL par_put_int   (ICHECKMIC,1)
    CALL par_put_int   (IMBUDGET,1)
    CALL par_put_int   (IRIME,1)
    CALL par_put_int   (IPLAWS,1)
    CALL par_put_int   (ISEDIM,1)
    CALL par_put_int   (ICLOUD,1)
    CALL par_put_int   (IDRIZ,1)
    CALL par_put_int   (IRAIN,1)
    CALL par_put_int   (IPRIS,1)
    CALL par_put_int   (ISNOW,1)
    CALL par_put_int   (IAGGR,1)
    CALL par_put_int   (IGRAUP,1)
    CALL par_put_int   (IHAIL,1)
    CALL par_put_float (CPARM,1)
    CALL par_put_float (DPARM,1)
    CALL par_put_float (RPARM,1)
    CALL par_put_float (PPARM,1)
    CALL par_put_float (SPARM,1)
    CALL par_put_float (APARM,1)
    CALL par_put_float (GPARM,1)
    CALL par_put_float (HPARM,1)
    CALL par_put_float (GNU,8)
    CALL par_put_char  (HUCMFILE,strl1)
    CALL par_put_int   (NDTCOLL,1)
    CALL par_put_int   (IAEROSOL,1)
    CALL par_put_int   (ISALT,1)
    CALL par_put_int   (IABCARB,1)
    CALL par_put_int   (IDUST,1)
    CALL par_put_int   (IDUSTLOFT,1)
    CALL par_put_char  (DUSTFILE,strl1)
    CALL par_put_int   (ICCNLEV,1)
    CALL par_put_int   (IIFN,1)
    CALL par_put_int   (IIFN_FORMULA,1)
    CALL par_put_int   (IAERORAD,1)
    CALL par_put_int   (IAERODEP,1)
    CALL par_put_int   (IAEROPRNT,1)
    CALL par_put_int   (IAEROHIST,1)
    CALL par_put_float (CIN_MAX,1)
    CALL par_put_float (CCN1_MAX,1)
    CALL par_put_float (CCN2_MAX,1)
    CALL par_put_float (DUST1_MAX,1)
    CALL par_put_float (DUST2_MAX,1)
    CALL par_put_float (SALTF_MAX,1)
    CALL par_put_float (SALTJ_MAX,1)
    CALL par_put_float (SALTS_MAX,1)
    CALL par_put_float (ABC1_MAX,1)
    CALL par_put_float (ABC2_MAX,1)
    CALL par_put_int   (IAEROLBC,MAXGRDS)
    CALL par_put_int   (ICO2LBC,MAXGRDS)
    CALL par_put_float (BCTAU,MAXGRDS)
    CALL par_put_int   (IAERO_CHEM,aerocat)
    CALL par_put_float (AERO_EPSILON,aerocat)
    CALL par_put_float (AERO_MEDRAD,aerocat)
    CALL par_put_int   (ITRKEPSILON,1)
    CALL par_put_int   (ITRKDUST,1)
    CALL par_put_int   (ITRKDUSTIFN,1)
    CALL par_put_int   (ICEPROCS,1)

    ! $MODEL_SOUND namelist group
    CALL par_put_int   (IPSFLG,1)
    CALL par_put_int   (ITSFLG,1)
    CALL par_put_int   (IRTSFLG,1)
    CALL par_put_int   (IUSFLG,1)
    CALL par_put_float (HS,MAXSNDG)
    CALL par_put_float (PS,MAXSNDG)
    CALL par_put_float (TS,MAXSNDG)
    CALL par_put_float (RTS,MAXSNDG)
    CALL par_put_float (US,MAXSNDG)
    CALL par_put_float (VS,MAXSNDG)

    ! $ISAN_CONTROL namelist group
    CALL par_put_int   (ISZSTAGE,1)
    CALL par_put_int   (IVRSTAGE,1)
    CALL par_put_int   (ISAN_INC,1)
    CALL par_put_int   (I1ST_FLG,1)
    CALL par_put_int   (IUPA_FLG,1)
    CALL par_put_int   (ISFC_FLG,1)
    CALL par_put_char  (IAPR,strl1)
    CALL par_put_char  (IARAWI,strl1)
    CALL par_put_char  (IASRFCE,strl1)
    CALL par_put_char  (VAR_HFILE,strl1)
    CALL par_put_char  (VARPFX,strl1)
    CALL par_put_int   (IOFLGISZ,1)
    CALL par_put_int   (IOFLGVAR,1)
    CALL par_put_int   (IDATAIN,1)

    ! $ISAN_ISENTROPIC namelist group
    CALL par_put_int   (NIGRIDS,1)
    CALL par_put_int   (NFEEDVAR,1)
    CALL par_put_float (SIGZWT,1)
    CALL par_put_int   (NISN,1)
    CALL par_put_int   (LEVTH,MAXISN)
    CALL par_put_float (TOPSIGZ,1)
    CALL par_put_float (HYBBOT,1)
    CALL par_put_float (HYBTOP,1)
    CALL par_put_float (SFCINF,1)
    CALL par_put_int   (MAXSTA,1)
    CALL par_put_int   (MAXSFC,1)
    CALL par_put_float (GOBSEP,1)
    CALL par_put_float (GOBRAD,1)
    CALL par_put_float (WVLNTH,MAXAGRDS)
    CALL par_put_float (SWVLNTH,MAXAGRDS)
    CALL par_put_float (RESPON,MAXAGRDS)
    CALL par_put_float (STASEP,1)
    CALL par_put_int   (IGRIDFL,1)
    CALL par_put_float (GRIDWT,MAXAGRDS)
    CALL par_put_int   (NOTSTA,1)
    CALL par_put_char  (NOTID,50*strl1)
    CALL par_put_char  (USED_FILE,strl1)
    CALL par_put_int   (IOBSWIN,1)

    ! eng_params() vars
    CALL par_put_float (SSPCT,1)
    CALL par_put_int   (IMPL,1)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! $MODEL_GRID namelist group
    CALL par_get_char  (EXPNME,strl1)
    CALL par_get_char  (RUNTYPE,32)
    CALL par_get_char  (TIMEUNIT,32)
    CALL par_get_float (TIMMAX,1)
    CALL par_get_int   (IMONTH1,1)
    CALL par_get_int   (IDATE1,1)
    CALL par_get_int   (IYEAR1,1)
    CALL par_get_int   (ITIME1,1)
    CALL par_get_int   (NGRIDS,1)
    CALL par_get_int   (NNXP,MAXGRDS)
    CALL par_get_int   (NNYP,MAXGRDS)
    CALL par_get_int   (NNZP,MAXGRDS)
    CALL par_get_int   (NZG,1)
    CALL par_get_int   (NZS,1)
    CALL par_get_int   (NXTNEST,MAXGRDS)
    CALL par_get_int   (IPRNTSTMT,1)
    CALL par_get_int   (IHTRAN,1)
    CALL par_get_float (DELTAX,1)
    CALL par_get_float (DELTAZ,1)
    CALL par_get_float (DZRAT,1)
    CALL par_get_float (DZMAX,1)
    CALL par_get_float (ZZ,NZPMAX)
    CALL par_get_float (DTLONG,1)
    CALL par_get_int   (NACOUST,1)
    CALL par_get_int   (NSTRATX,MAXGRDS)
    CALL par_get_int   (NSTRATY,MAXGRDS)
    CALL par_get_int   (NNDTRAT,MAXGRDS)
    CALL par_get_int   (NESTZ,1)
    CALL par_get_int   (NSTRATZ,NZPMAX)
    CALL par_get_float (POLELAT,1)
    CALL par_get_float (POLELON,1)
    CALL par_get_int   (NINEST,MAXGRDS)
    CALL par_get_int   (NJNEST,MAXGRDS)
    CALL par_get_int   (NKNEST,MAXGRDS)
    CALL par_get_float (CENTLAT,MAXGRDS)
    CALL par_get_float (CENTLON,MAXGRDS)
    CALL par_get_int   (NNSTTOP,MAXGRDS)
    CALL par_get_int   (NNSTBOT,MAXGRDS)
    
    ! $MODEL_FILE_INFO namelist group
    CALL par_get_int   (INITIAL,1)
    CALL par_get_int   (NUD_TYPE,1)
    CALL par_get_char  (VARFPFX,strl1)
    CALL par_get_float (VWAIT1,1)
    CALL par_get_float (VWAITTOT,1)
    CALL par_get_int   (NUDLAT,1)
    CALL par_get_float (TNUDLAT,1)
    CALL par_get_float (TNUDCENT,1)
    CALL par_get_float (TNUDTOP,1)
    CALL par_get_float (ZNUDTOP,1)
    CALL par_get_float (WT_NUDGE_G,maxgrds)
    CALL par_get_float (WT_NUDGE_UV,1)
    CALL par_get_float (WT_NUDGE_TH,1)
    CALL par_get_float (WT_NUDGE_PI,1)
    CALL par_get_float (WT_NUDGE_RT,1)
    CALL par_get_int   (NUD_COND,1)
    CALL par_get_float (TCOND_BEG,1)
    CALL par_get_float (TCOND_END,1)
    CALL par_get_float (T_NUDGE_RC,1)
    CALL par_get_float (WT_NUDGEC,maxgrds)
    CALL par_get_float (SNUDCENT,1)
    CALL par_get_int   (IF_ODA,1)
    CALL par_get_char  (ODA_UPAPREF,strl1)
    CALL par_get_char  (ODA_SFCPREF,strl1)
    CALL par_get_float (FRQODA,1)
    CALL par_get_float (TODABEG,1)
    CALL par_get_float (TODAEND,1)
    CALL par_get_float (TNUDODA,1)
    CALL par_get_float (WT_ODA_GRID,MAXGRDS)
    CALL par_get_float (WT_ODA_UV,1)
    CALL par_get_float (WT_ODA_TH,1)
    CALL par_get_float (WT_ODA_PI,1)
    CALL par_get_float (WT_ODA_RT,1)
    CALL par_get_float (RODA_SFCE,MAXGRDS)
    CALL par_get_float (RODA_SFC0,MAXGRDS)
    CALL par_get_float (RODA_UPAE,MAXGRDS)
    CALL par_get_float (RODA_UPA0,MAXGRDS)
    CALL par_get_float (RODA_HGT,MAXGRDS)
    CALL par_get_float (RODA_ZFAC,MAXGRDS)
    CALL par_get_float (ODA_SFC_TIL,1)
    CALL par_get_float (ODA_SFC_TEL,1)
    CALL par_get_float (ODA_UPA_TIL,1)
    CALL par_get_float (ODA_UPA_TEL,1)
    CALL par_get_char  (HFILIN,strl1)
    CALL par_get_int   (IPAST_SFC,1)
    CALL par_get_int   (ICLOBBER,1)
    CALL par_get_int   (IOUTPUT,1)
    CALL par_get_char  (AFILEPREF,strl1)
    CALL par_get_float (FRQSTATE,MAXGRDS)
    CALL par_get_float (FRQST_KEEP,1)
    CALL par_get_float (FRQLITE,1)
    CALL par_get_float (NLITE_VARS,1)
    do nm = 1, nlite_vars
       CALL par_get_char (LITE_VARS(nm),32)
    enddo
    CALL par_get_float (AVGTIM,1)
    CALL par_get_float (FRQMEAN,1)
    CALL par_get_float (FRQBOTH,1)
    CALL par_get_char  (TOPFILES,strl1)
    CALL par_get_char  (SFCFILES,strl1)
    CALL par_get_char  (SSTFPFX,strl1)
    CALL par_get_char  (NDVIFPFX,strl1)
    CALL par_get_int   (ITOPTFLG,MAXGRDS)
    CALL par_get_int   (ISSTFLG,MAXGRDS)
    CALL par_get_int   (IVEGTFLG,MAXGRDS)
    CALL par_get_int   (ISOILFLG,MAXGRDS)
    CALL par_get_int   (NDVIFLG,MAXGRDS)
    CALL par_get_int   (IUPDNDVI,1)
    CALL par_get_int   (IUPDSST,1)
    CALL par_get_char  (ITOPTFN,MAXGRDS*strl1)
    CALL par_get_char  (ISSTFN,MAXGRDS*strl1)
    CALL par_get_char  (IVEGTFN,MAXGRDS*strl1)
    CALL par_get_char  (ISOILFN,MAXGRDS*strl1)
    CALL par_get_char  (NDVIFN,MAXGRDS*strl1)
    CALL par_get_int   (ITOPSFLG,MAXGRDS)
    CALL par_get_float (TOPTENH,MAXGRDS)
    CALL par_get_float (TOPTWVL,MAXGRDS)
    CALL par_get_int   (IZ0FLG,MAXGRDS)
    CALL par_get_float (Z0MAX,MAXGRDS)
    CALL par_get_float (Z0FACT,1)

    ! $MODEL_OPTIONS namelist group
    CALL par_get_int   (ICORFLG,1)
    CALL par_get_int   (IBND,1)
    CALL par_get_int   (JBND,1)
    CALL par_get_int   (ISPONGE_PTS,MAXGRDS)
    CALL par_get_float (SPONGE_TAU,MAXGRDS)
    CALL par_get_float (CPHAS,1)
    CALL par_get_int   (LSFLG,1)
    CALL par_get_int   (NFPT,1)
    CALL par_get_float (DISTIM,1)
    CALL par_get_int   (ILWRTYP,1)
    CALL par_get_int   (ISWRTYP,1)
    CALL par_get_float (RADFRQ,1)
    CALL par_get_int   (LONRAD,1)
    CALL par_get_int   (NNQPARM,MAXGRDS)
    CALL par_get_float (CONFRQ,1)
    CALL par_get_float (WCLDBS,1)
    CALL par_get_int   (IKPP,1)
    CALL par_get_int   (NKPPZ,1)
    CALL par_get_float (FRQKPP,1)
    CALL par_get_float (RELAX_SST,1)
    CALL par_get_float (RELAX_OCNT,1)
    CALL par_get_float (RELAX_SAL,1)
    CALL par_get_float (DMAXKPP,1)
    CALL par_get_float (DSCALEKPP,1)
    CALL par_get_int   (KPPITERMAX,1)
    CALL par_get_int   (KPPRNT,1)
    CALL par_get_float (UBMN_KPP,1)
    CALL par_get_int   (NPATCH,1)
    CALL par_get_int   (NVEGPAT,1)
    CALL par_get_int   (ISFCL,1)
    CALL par_get_int   (IFREESLIP,1)
    CALL par_get_char  (SIBFILE,strl1)
    CALL par_get_float (CO2_INIT,NZPMAX)
    CALL par_get_int   (ISOILDAT,1)
    CALL par_get_int   (ISNOWDAT,1)
    CALL par_get_int   (NVGCON,1)
    CALL par_get_float (PCTLCON,1)
    CALL par_get_int   (NSLCON,1)
    CALL par_get_float (ZROUGH,1)
    CALL par_get_float (ALBEDO,1)
    CALL par_get_float (SEATMP,1)
    CALL par_get_float (DTHCON,1)
    CALL par_get_float (DRTCON,1)
    CALL par_get_float (SLZ,NZGMAX)
    CALL par_get_float (SLMSTR,NZGMAX)
    CALL par_get_float (STGOFF,NZGMAX)
    CALL par_get_int   (IDIFFK,MAXGRDS)
    CALL par_get_int   (IDIFFPERTS,1)
    CALL par_get_int   (IHORGRAD,1)
    CALL par_get_float (CSX,MAXGRDS)
    CALL par_get_float (CSZ,MAXGRDS)
    CALL par_get_float (XKHKM,MAXGRDS)
    CALL par_get_float (ZKHKM,MAXGRDS)
    CALL par_get_float (AKMIN,MAXGRDS)
    CALL par_get_float (FRACSAT,1)
    CALL par_get_int   (IBUBBLE,1)
    CALL par_get_int   (IBUBGRD,1)
    CALL par_get_int   (IBDXIA,1)
    CALL par_get_int   (IBDXIZ,1)
    CALL par_get_int   (IBDYJA,1)
    CALL par_get_int   (IBDYJZ,1)
    CALL par_get_int   (IBDZK1,1)
    CALL par_get_int   (IBDZK2,1)
    CALL par_get_float (BTHP,1)
    CALL par_get_float (BRTP,1)
    CALL par_get_int   (ICONV,1)
    CALL par_get_int   (ICONGR,1)
    CALL par_get_int   (ICICENT,1)
    CALL par_get_int   (ICJCENT,1)
    CALL par_get_float (CXRAD,1)
    CALL par_get_float (CYRAD,1)
    CALL par_get_int   (ICVERT,1)
    CALL par_get_int   (ICKMAX,1)
    CALL par_get_float (CZRAD,1)
    CALL par_get_int   (ICKCENT,1)
    CALL par_get_float (CDIVMAX,1)
    CALL par_get_float (CTAU,1)
    CALL par_get_float (CTMAX,1)
    CALL par_get_int   (IRCE,1)
    CALL par_get_float (RCE_SZEN,1)
    CALL par_get_float (RCE_SOLC,1)
    CALL par_get_float (RCE_UBMN,1)
    CALL par_get_float (RCE_BUBL,1)
    CALL par_get_int   (ITRACER,1)
    CALL par_get_int   (ITRACHIST,1)
    CALL par_get_int   (LEVEL,1)
    CALL par_get_int   (ISCM,1)
    CALL par_get_float (SCMTIME,1)
    CALL par_get_int   (ISCMX,1)
    CALL par_get_int   (ISCMY,1)
    CALL par_get_int   (ICHECKMIC,1)
    CALL par_get_int   (IMBUDGET,1)
    CALL par_get_int   (IRIME,1)
    CALL par_get_int   (IPLAWS,1)
    CALL par_get_int   (ISEDIM,1)
    CALL par_get_int   (ICLOUD,1)
    CALL par_get_int   (IDRIZ,1)
    CALL par_get_int   (IRAIN,1)
    CALL par_get_int   (IPRIS,1)
    CALL par_get_int   (ISNOW,1)
    CALL par_get_int   (IAGGR,1)
    CALL par_get_int   (IGRAUP,1)
    CALL par_get_int   (IHAIL,1)
    CALL par_get_float (CPARM,1)
    CALL par_get_float (DPARM,1)
    CALL par_get_float (RPARM,1)
    CALL par_get_float (PPARM,1)
    CALL par_get_float (SPARM,1)
    CALL par_get_float (APARM,1)
    CALL par_get_float (GPARM,1)
    CALL par_get_float (HPARM,1)
    CALL par_get_float (GNU,8)
    CALL par_get_char  (HUCMFILE,strl1)
    CALL par_get_int   (NDTCOLL,1)
    CALL par_get_int   (IAEROSOL,1)
    CALL par_get_int   (ISALT,1)
    CALL par_get_int   (IABCARB,1)
    CALL par_get_int   (IDUST,1)
    CALL par_get_int   (IDUSTLOFT,1)
    CALL par_get_char  (DUSTFILE,strl1)
    CALL par_get_int   (ICCNLEV,1)
    CALL par_get_int   (IIFN,1)
    CALL par_get_int   (IIFN_FORMULA,1)
    CALL par_get_int   (IAERORAD,1)
    CALL par_get_int   (IAERODEP,1)
    CALL par_get_int   (IAEROPRNT,1)
    CALL par_get_int   (IAEROHIST,1)
    CALL par_get_float (CIN_MAX,1)
    CALL par_get_float (CCN1_MAX,1)
    CALL par_get_float (CCN2_MAX,1)
    CALL par_get_float (DUST1_MAX,1)
    CALL par_get_float (DUST2_MAX,1)
    CALL par_get_float (SALTF_MAX,1)
    CALL par_get_float (SALTJ_MAX,1)
    CALL par_get_float (SALTS_MAX,1)
    CALL par_get_float (ABC1_MAX,1)
    CALL par_get_float (ABC2_MAX,1)
    CALL par_get_int   (IAEROLBC,MAXGRDS)
    CALL par_get_int   (ICO2LBC,MAXGRDS)
    CALL par_get_float (BCTAU,MAXGRDS)
    CALL par_get_int   (IAERO_CHEM,aerocat)
    CALL par_get_float (AERO_EPSILON,aerocat)
    CALL par_get_float (AERO_MEDRAD,aerocat)
    CALL par_get_int   (ITRKEPSILON,1)
    CALL par_get_int   (ITRKDUST,1)
    CALL par_get_int   (ITRKDUSTIFN,1)
    CALL par_get_int   (ICEPROCS,1)

    ! $MODEL_SOUND namelist group
    CALL par_get_int   (IPSFLG,1)
    CALL par_get_int   (ITSFLG,1)
    CALL par_get_int   (IRTSFLG,1)
    CALL par_get_int   (IUSFLG,1)
    CALL par_get_float (HS,MAXSNDG)
    CALL par_get_float (PS,MAXSNDG)
    CALL par_get_float (TS,MAXSNDG)
    CALL par_get_float (RTS,MAXSNDG)
    CALL par_get_float (US,MAXSNDG)
    CALL par_get_float (VS,MAXSNDG)

    ! $ISAN_CONTROL namelist group
    CALL par_get_int   (ISZSTAGE,1)
    CALL par_get_int   (IVRSTAGE,1)
    CALL par_get_int   (ISAN_INC,1)
    CALL par_get_int   (I1ST_FLG,1)
    CALL par_get_int   (IUPA_FLG,1)
    CALL par_get_int   (ISFC_FLG,1)
    CALL par_get_char  (IAPR,strl1)
    CALL par_get_char  (IARAWI,strl1)
    CALL par_get_char  (IASRFCE,strl1)
    CALL par_get_char  (VAR_HFILE,strl1)
    CALL par_get_char  (VARPFX,strl1)
    CALL par_get_int   (IOFLGISZ,1)
    CALL par_get_int   (IOFLGVAR,1)
    CALL par_get_int   (IDATAIN,1)

    ! $ISAN_ISENTROPIC namelist group
    CALL par_get_int   (NIGRIDS,1)
    CALL par_get_int   (NFEEDVAR,1)
    CALL par_get_float (SIGZWT,1)
    CALL par_get_int   (NISN,1)
    CALL par_get_int   (LEVTH,MAXISN)
    CALL par_get_float (TOPSIGZ,1)
    CALL par_get_float (HYBBOT,1)
    CALL par_get_float (HYBTOP,1)
    CALL par_get_float (SFCINF,1)
    CALL par_get_int   (MAXSTA,1)
    CALL par_get_int   (MAXSFC,1)
    CALL par_get_float (GOBSEP,1)
    CALL par_get_float (GOBRAD,1)
    CALL par_get_float (WVLNTH,MAXAGRDS)
    CALL par_get_float (SWVLNTH,MAXAGRDS)
    CALL par_get_float (RESPON,MAXAGRDS)
    CALL par_get_float (STASEP,1)
    CALL par_get_int   (IGRIDFL,1)
    CALL par_get_float (GRIDWT,MAXAGRDS)
    CALL par_get_int   (NOTSTA,1)
    CALL par_get_char  (NOTID,50*strl1)
    CALL par_get_char  (USED_FILE,strl1)
    CALL par_get_int   (IOBSWIN,1)

    ! eng_params() vars
    CALL par_get_float (SSPCT,1)
    CALL par_get_int   (IMPL,1)

  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_config

!##########################################################################
Subroutine broadcast_grid ()

use mem_grid
use node_mod

implicit none

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 50 + (20 + 10*nxpmax + 10*nypmax + 19*nzpmax) * maxgrds
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! grid vars
    CALL par_put_int (nnx,maxgrds)
    CALL par_put_int (nnx1,maxgrds)
    CALL par_put_int (nnx2,maxgrds)
    CALL par_put_int (nny,maxgrds)
    CALL par_put_int (nny1,maxgrds)
    CALL par_put_int (nny2,maxgrds)
    CALL par_put_int (nnz,maxgrds)
    CALL par_put_int (nnz1,maxgrds)
    CALL par_put_int (nnxyzp,maxgrds)
    CALL par_put_int (nnxysp,maxgrds)
    CALL par_put_int (nnxyp,maxgrds)
    CALL par_put_int (jdim,1)

    ! coordinate values
    CALL par_put_int   (nrz,nzpmax*maxgrds)
    CALL par_put_int   (ipm,nxpmax*maxgrds)
    CALL par_put_int   (jpm,nypmax*maxgrds)
    CALL par_put_int   (kpm,nzpmax*maxgrds)
    CALL par_put_float (xmn,nxpmax*maxgrds)
    CALL par_put_float (ymn,nypmax*maxgrds)
    CALL par_put_float (zmn,nzpmax*maxgrds)
    CALL par_put_float (xtn,nxpmax*maxgrds)
    CALL par_put_float (ytn,nypmax*maxgrds)
    CALL par_put_float (ztn,nzpmax*maxgrds)
    CALL par_put_float (dzmn,nzpmax*maxgrds)
    CALL par_put_float (dzm2n,nzpmax*maxgrds)
    CALL par_put_float (dztn,nzpmax*maxgrds)
    CALL par_put_float (dzt2n,nzpmax*maxgrds)
    CALL par_put_float (deltaxn,maxgrds)
    CALL par_put_float (deltazn,maxgrds)
    CALL par_put_float (ztop,1)

    ! nesting coefficients
    CALL par_put_float (ei1,nxpmax*maxgrds)
    CALL par_put_float (ei2,nxpmax*maxgrds)
    CALL par_put_float (ei3,nxpmax*maxgrds)
    CALL par_put_float (ei4,nxpmax*maxgrds)
    CALL par_put_float (ei5,nxpmax*maxgrds)
    CALL par_put_float (ei6,nxpmax*maxgrds)
    CALL par_put_float (ei7,nxpmax*maxgrds)
    
    CALL par_put_float (ej1,nypmax*maxgrds)
    CALL par_put_float (ej2,nypmax*maxgrds)
    CALL par_put_float (ej3,nypmax*maxgrds)
    CALL par_put_float (ej4,nypmax*maxgrds)
    CALL par_put_float (ej5,nypmax*maxgrds)
    CALL par_put_float (ej6,nypmax*maxgrds)
    CALL par_put_float (ej7,nypmax*maxgrds)
    
    CALL par_put_float (ek1,nzpmax*maxgrds)
    CALL par_put_float (ek2,nzpmax*maxgrds)
    CALL par_put_float (ek3,nzpmax*maxgrds)
    CALL par_put_float (ek4,nzpmax*maxgrds)
    CALL par_put_float (ek5,nzpmax*maxgrds)
    CALL par_put_float (ek6,nzpmax*maxgrds)
    CALL par_put_float (ek7,nzpmax*maxgrds)
    CALL par_put_float (fbcf,nzpmax*maxgrds*4)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! grid vars
    CALL par_get_int (nnx,maxgrds)
    CALL par_get_int (nnx1,maxgrds)
    CALL par_get_int (nnx2,maxgrds)
    CALL par_get_int (nny,maxgrds)
    CALL par_get_int (nny1,maxgrds)
    CALL par_get_int (nny2,maxgrds)
    CALL par_get_int (nnz,maxgrds)
    CALL par_get_int (nnz1,maxgrds)
    CALL par_get_int (nnxyzp,maxgrds)
    CALL par_get_int (nnxysp,maxgrds)
    CALL par_get_int (nnxyp,maxgrds)
    CALL par_get_int (jdim,1)

    ! coordinate values
    CALL par_get_int   (nrz,nzpmax*maxgrds)
    CALL par_get_int   (ipm,nxpmax*maxgrds)
    CALL par_get_int   (jpm,nypmax*maxgrds)
    CALL par_get_int   (kpm,nzpmax*maxgrds)
    CALL par_get_float (xmn,nxpmax*maxgrds)
    CALL par_get_float (ymn,nypmax*maxgrds)
    CALL par_get_float (zmn,nzpmax*maxgrds)
    CALL par_get_float (xtn,nxpmax*maxgrds)
    CALL par_get_float (ytn,nypmax*maxgrds)
    CALL par_get_float (ztn,nzpmax*maxgrds)
    CALL par_get_float (dzmn,nzpmax*maxgrds)
    CALL par_get_float (dzm2n,nzpmax*maxgrds)
    CALL par_get_float (dztn,nzpmax*maxgrds)
    CALL par_get_float (dzt2n,nzpmax*maxgrds)
    CALL par_get_float (deltaxn,maxgrds)
    CALL par_get_float (deltazn,maxgrds)
    CALL par_get_float (ztop,1)

    ! nesting coefficients
    CALL par_get_float (ei1,nxpmax*maxgrds)
    CALL par_get_float (ei2,nxpmax*maxgrds)
    CALL par_get_float (ei3,nxpmax*maxgrds)
    CALL par_get_float (ei4,nxpmax*maxgrds)
    CALL par_get_float (ei5,nxpmax*maxgrds)
    CALL par_get_float (ei6,nxpmax*maxgrds)
    CALL par_get_float (ei7,nxpmax*maxgrds)
    
    CALL par_get_float (ej1,nypmax*maxgrds)
    CALL par_get_float (ej2,nypmax*maxgrds)
    CALL par_get_float (ej3,nypmax*maxgrds)
    CALL par_get_float (ej4,nypmax*maxgrds)
    CALL par_get_float (ej5,nypmax*maxgrds)
    CALL par_get_float (ej6,nypmax*maxgrds)
    CALL par_get_float (ej7,nypmax*maxgrds)
    
    CALL par_get_float (ek1,nzpmax*maxgrds)
    CALL par_get_float (ek2,nzpmax*maxgrds)
    CALL par_get_float (ek3,nzpmax*maxgrds)
    CALL par_get_float (ek4,nzpmax*maxgrds)
    CALL par_get_float (ek5,nzpmax*maxgrds)
    CALL par_get_float (ek6,nzpmax*maxgrds)
    CALL par_get_float (ek7,nzpmax*maxgrds)
    CALL par_get_float (fbcf,nzpmax*maxgrds*4)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_grid

!##########################################################################
Subroutine broadcast_decomp ()

use node_mod

implicit none

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 50 + (4 * maxgrds) + (10 * maxmach * maxgrds)
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! node vars
    CALL par_put_int (ixb,maxmach*maxgrds) ! these are arrays with size maxmach*maxgrds
    CALL par_put_int (ixe,maxmach*maxgrds) !  that were set in assign_node_subdomain
    CALL par_put_int (iyb,maxmach*maxgrds)
    CALL par_put_int (iye,maxmach*maxgrds)
    CALL par_put_int (fd_sw_num,maxgrds)
    CALL par_put_int (fd_nw_num,maxgrds)
    CALL par_put_int (fd_se_num,maxgrds)
    CALL par_put_int (fd_ne_num,maxgrds)

    ! broadcast node that contains SCM output point    
    CALL par_put_int (my_scm_num,1)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! node vars
    CALL par_get_int (ixb,maxmach*maxgrds)
    CALL par_get_int (ixe,maxmach*maxgrds)
    CALL par_get_int (iyb,maxmach*maxgrds)
    CALL par_get_int (iye,maxmach*maxgrds)
    CALL par_get_int (fd_sw_num,maxgrds)
    CALL par_get_int (fd_nw_num,maxgrds)
    CALL par_get_int (fd_se_num,maxgrds)
    CALL par_get_int (fd_ne_num,maxgrds)

    ! broadcast node that contains SCM output point
    CALL par_get_int (my_scm_num,1)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_decomp

!##########################################################################
Subroutine broadcast_bub_rand_nums (npts,bub_rand_nums)

! Random numbers for T perturbations are only generated by one processor

use node_mod

implicit none

  integer :: npts
  real, dimension(npts) :: bub_rand_nums

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 50 + npts
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)

    ! random numbers
    CALL par_put_float (bub_rand_nums,npts)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! random numbers
    CALL par_get_float (bub_rand_nums,npts)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_bub_rand_nums

!##############################################################################
Subroutine broadcast_hist_header (ngrids1,nnxp1,nnyp1,nnzp1,nzg1 &
                                 ,nzs1,nkppz1,npatch1)

use mem_grid
use node_mod
use ref_sounding
use an_header
use kpp_parameters, only:nkppz

implicit none

  integer :: ngrids1, nzg1, nzs1, nkppz1, npatch1
  integer, dimension(maxgrds) :: nnxp1, nnyp1, nnzp1

  real, allocatable :: buff(:)
  integer :: nwords
  integer :: atable_elem_size, atable_elem_words

  type (head_table) :: dummy_atable

  ! this will produce number of bytes

  ! get the size of one record (an element of atable)
  atable_elem_size = sizeof(dummy_atable) 
  atable_elem_words = atable_elem_size / 4

  !Saleeby(2018): Increment "nwords" buffer if added "par_put" and
  !"par_get" calls below.
  nwords = 17 + ((3 + (6 * nzpmax)) * maxgrds) + (6 * maxsndg) &
           + (maxvars * atable_elem_words)

  ! note that what got allocated was nwords*sizeof(real) bytes
  allocate (buff(nwords))

  ! The mainnum process sends data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! former grid structure
    CALL par_put_int   (ngrids1,1)
    CALL par_put_int   (nnxp1,maxgrds)
    CALL par_put_int   (nnyp1,maxgrds)
    CALL par_put_int   (nnzp1,maxgrds)
    CALL par_put_int   (npatch1,1)
    CALL par_put_int   (nzg1,1)
    CALL par_put_int   (nzs1,1)
    CALL par_put_int   (nkppz,1)
    CALL par_put_float (time,1)

    ! Flag to determine if we are past the first timestep
    CALL par_put_int   (ngbegun,ngrids)
    ! Flag to retain initial runtype if we are doing history initialization
    CALL par_put_int   (initorig,1)

    ! Reference sounding
    CALL par_put_int   (iref,1)
    CALL par_put_int   (jref,1)
    CALL par_put_float (topref,1)
    CALL par_put_int   (nsndg,1)
    CALL par_put_float (us,maxsndg)
    CALL par_put_float (vs,maxsndg)
    CALL par_put_float (ts,maxsndg)
    CALL par_put_float (thds,maxsndg)
    CALL par_put_float (ps,maxsndg)
    CALL par_put_float (hs,maxsndg)

    ! 1-D reference state for all grids
    CALL par_put_float (u01dn,nzpmax*maxgrds)
    CALL par_put_float (v01dn,nzpmax*maxgrds)
    CALL par_put_float (pi01dn,nzpmax*maxgrds)
    CALL par_put_float (th01dn,nzpmax*maxgrds)
    CALL par_put_float (dn01dn,nzpmax*maxgrds)
    CALL par_put_float (rt01dn,nzpmax*maxgrds)

    ! anal_table - array of structures
    ! use par_put_char since sizeof returns number of bytes
    CALL par_put_int   (nvbtab,1)
    CALL par_put_char  (anal_table,nvbtab*atable_elem_size)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! former grid structure
    CALL par_get_int   (ngrids1,1)
    ngridsh = ngrids1
    CALL par_get_int   (nnxp1,maxgrds)
    CALL par_get_int   (nnyp1,maxgrds)
    CALL par_get_int   (nnzp1,maxgrds)
    CALL par_get_int   (npatch1,1)
    CALL par_get_int   (nzg1,1)
    CALL par_get_int   (nzs1,1)
    CALL par_get_int   (nkppz1,1)
    CALL par_get_float (time,1)

    ! Flag to determine if we are past the first timestep
    CALL par_get_int   (ngbegun,ngrids)
    ! Flag to retain initial runtype if we are doing history initialization
    CALL par_get_int   (initorig,1)

    ! Reference sounding
    CALL par_get_int   (iref,1)
    CALL par_get_int   (jref,1)
    CALL par_get_float (topref,1)
    CALL par_get_int   (nsndg,1)
    CALL par_get_float (us,maxsndg)
    CALL par_get_float (vs,maxsndg)
    CALL par_get_float (ts,maxsndg)
    CALL par_get_float (thds,maxsndg)
    CALL par_get_float (ps,maxsndg)
    CALL par_get_float (hs,maxsndg)

    ! 1-D reference state for all grids
    CALL par_get_float (u01dn,nzpmax*maxgrds)
    CALL par_get_float (v01dn,nzpmax*maxgrds)
    CALL par_get_float (pi01dn,nzpmax*maxgrds)
    CALL par_get_float (th01dn,nzpmax*maxgrds)
    CALL par_get_float (dn01dn,nzpmax*maxgrds)
    CALL par_get_float (rt01dn,nzpmax*maxgrds)

    ! anal_table - array of structures
    ! use par_get_char since sizeof returns number of bytes
    CALL par_get_int   (nvbtab,1)
    allocate (anal_table(nvbtab))
    CALL par_get_char  (anal_table,nvbtab*atable_elem_size)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_hist_header

!##########################################################################
Subroutine broadcast_vfile_selection ()

use mem_grid
use node_mod
use mem_varinit

implicit none

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 20 + (maxnudfiles * (strl1 + 15))
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! var file selection
    CALL par_put_int   (nvarffiles,1)
    CALL par_put_int   (nvarffl,1)
    CALL par_put_char  (fnames_varf,maxnudfiles*strl1)
    CALL par_put_char  (itotdate_varf,maxnudfiles*14)
    CALL par_put_float (varf_times,maxnudfiles)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! var file selection
    CALL par_get_int   (nvarffiles,1)
    CALL par_get_int   (nvarffl,1)
    CALL par_get_char  (fnames_varf,maxnudfiles*strl1)
    CALL par_get_char  (itotdate_varf,maxnudfiles*14)
    CALL par_get_float (varf_times,maxnudfiles)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_vfile_selection

!##########################################################################
Subroutine broadcast_vfile_refstate ()

use mem_grid
use node_mod
use ref_sounding

implicit none

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 20 + (6 * nzpmax * maxgrds)
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! var 1D ref state
    CALL par_put_int   (iref,1)
    CALL par_put_int   (jref,1)
    CALL par_put_float (topref,1)
    CALL par_put_float (u01dn,nzpmax*maxgrds)
    CALL par_put_float (v01dn,nzpmax*maxgrds)
    CALL par_put_float (rt01dn,nzpmax*maxgrds)
    CALL par_put_float (th01dn,nzpmax*maxgrds)
    CALL par_put_float (pi01dn,nzpmax*maxgrds)
    CALL par_put_float (dn01dn,nzpmax*maxgrds)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! var 1D ref state
    CALL par_get_int   (iref,1)
    CALL par_get_int   (jref,1)
    CALL par_get_float (topref,1)
    CALL par_get_float (u01dn,nzpmax*maxgrds)
    CALL par_get_float (v01dn,nzpmax*maxgrds)
    CALL par_get_float (rt01dn,nzpmax*maxgrds)
    CALL par_get_float (th01dn,nzpmax*maxgrds)
    CALL par_get_float (pi01dn,nzpmax*maxgrds)
    CALL par_get_float (dn01dn,nzpmax*maxgrds)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_vfile_refstate

!############################################################################
Subroutine get_fd_corner_vals (igrid,sw_val,nw_val,se_val,ne_val)

! This routine will find the max of the dxt values sampled from the
! four corners of the domain.

use mem_grid
use node_mod

implicit none

  integer :: igrid
  real :: sw_val, nw_val, se_val, ne_val

  integer :: inode
  real, dimension(:), allocatable :: all_sw_val, all_nw_val, all_se_val, all_ne_val

  ! If we are running in parallel, we need to exchange
  ! our local corner values with those at the full domain
  ! corners.
  !
  if (nmachs .gt. 1) then
    ! Allocate MPI buffers
    allocate(all_sw_val(nmachs))
    allocate(all_nw_val(nmachs))
    allocate(all_se_val(nmachs))
    allocate(all_ne_val(nmachs))

    ! Mainnum collects the corner values, adjusts corner dxt values and
    ! broadcasts to other nodes.
    if (my_rams_num .eq. mainnum) then
      CALL par_gather_floats (sw_val, all_sw_val, 1, machnum(mainnum))
      CALL par_gather_floats (nw_val, all_nw_val, 1, machnum(mainnum))
      CALL par_gather_floats (se_val, all_se_val, 1, machnum(mainnum))
      CALL par_gather_floats (ne_val, all_ne_val, 1, machnum(mainnum))

      ! Find the domain corners and reset local dxt values accordingly
      do inode = 1, nmachs
        sw_val = all_sw_val(fd_sw_num(igrid))  ! south west corner
        nw_val = all_nw_val(fd_nw_num(igrid))  ! north west corner
        se_val = all_se_val(fd_se_num(igrid))  ! south east corner
        ne_val = all_ne_val(fd_ne_num(igrid))  ! north east corner
      enddo

      ! Broadcast the new dxt_* values to the other nodes
      CALL broadcast_fd_corners (sw_val,nw_val,se_val,ne_val)
    else
      CALL par_gather_floats (sw_val, all_sw_val, 1, machnum(mainnum))
      CALL par_gather_floats (nw_val, all_nw_val, 1, machnum(mainnum))
      CALL par_gather_floats (se_val, all_se_val, 1, machnum(mainnum))
      CALL par_gather_floats (ne_val, all_ne_val, 1, machnum(mainnum))

      ! Recieve the new dxt_* values from mainnum
      CALL broadcast_fd_corners (sw_val,nw_val,se_val,ne_val)
    endif

    ! Deallocate MPI buffers
    deallocate(all_sw_val)
    deallocate(all_nw_val)
    deallocate(all_se_val)
    deallocate(all_ne_val)
  endif

return
END SUBROUTINE get_fd_corner_vals

!##########################################################################
Subroutine broadcast_fd_corners (sw_val,nw_val,se_val,ne_val)

use node_mod

implicit none

  real :: sw_val, nw_val, se_val, ne_val

  real, allocatable :: buff(:)
  integer :: nwords

  nwords = 20 
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)
    
    ! domain corners
    CALL par_put_float (sw_val,1)
    CALL par_put_float (nw_val,1)
    CALL par_put_float (se_val,1)
    CALL par_put_float (ne_val,1)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff, nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! domain corners
    CALL par_get_float (sw_val,1)
    CALL par_get_float (nw_val,1)
    CALL par_get_float (se_val,1)
    CALL par_get_float (ne_val,1)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_fd_corners

!##########################################################################
Subroutine broadcast_dustsource (source,nx_source,ny_source)

use node_mod

implicit none

  real, allocatable :: buff(:)
  integer :: nwords,nx_source,ny_source
  real, dimension(nx_source,ny_source) :: source

  nwords = 50 + nx_source*ny_source
  allocate (buff(nwords)) ! note that what got allocated was nwords*sizeof(real) bytes

  ! The mainnum process will send data, all others will receive in the broadcast
  !
  ! Therefore, the algorithm is:
  !   For mainnum
  !     pack buffer
  !     broadcast buffer to other processes
  !   For the other processes
  !     receive buffer from mainnum
  !     unpack buffer
  if (my_rams_num .eq. mainnum) then
    ! sending node: pack data into the buffer, then send
    CALL par_init_put (buff,nwords)

    ! dust source values
    CALL par_put_float (source,nx_source*ny_source)

    CALL par_broadcast (machnum(mainnum))
  else
    ! receiving nodes
    CALL par_init_recv_bcast (buff,nwords) ! don't call this if you are mainnum
    CALL par_broadcast (machnum(mainnum))

    ! random numbers
    CALL par_get_float (source,nx_source*ny_source)
  endif

  deallocate (buff)

return
END SUBROUTINE broadcast_dustsource
