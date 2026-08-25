

! CVS:  $Id: driver_read.F,v 1.10 2003/11/11 21:55:13 norm Exp $
! CVS:  $Name:  $


!-----------------------------------------------------------------------
      subroutine bugs_driver(nlm,amu0,alvdr,pl,tl,ql,qcwl,ncwl,qcil &
      ,     qrwl,qril,o3l,asl,atl,fuswb,fdswb,fulwb,fdlwb)

      use bugs_kinds, only:  int_kind, dbl_kind
      use bugsrad_physconst, only:  gravity, cp_dry_air, sol_const

      implicit none

!-----------------------------------------------------------------------
! driver_read is the main routine for running the CSU radiative transfer
! code offline (that is, apart from the CSU GCM).  It reads a profile
! from a file and also specifies variables that are not read in.  It
! then calls BUGSrad to do the radiative transfer.  This driver is not
! used when the code is compiled online with the CSU GCM.

! REFERENCES:
! Phil Partain /wombat (04-04-00).

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     bugs_rad   :The radiative transfer code

! FUNCTIONS CALLED:
!     none.
 
! INCLUDED COMMON BLOCKS:
!     none.
 
! LOCAL VARIABLES:
      integer (kind=int_kind):: &
       nlen,&      !Length of total domain.
       len,&       !Length of sub domain.
       nlm,&       !Number of layers.
       i,l

      real (kind=dbl_kind), dimension(1):: &
        ts   ,&      !Surface temperature                             (K).
        amu0 ,&      !Cosine of solar zenith angle                    (-).
        slr  ,&      !Fraction of daylight                            (-).
        alvdr,&      !Visible direct surface albedo                   (-).
        alndr,&      !Near-IR direct surface albedo                   (-).
        alvdf,&      !Visible diffuse surface albedo                  (-).
        alndf,&      !Near-IR diffuse surface albedo                  (-).
        umco2,&      !Col-avg concentration CO2                     (ppm).
        umch4,&      !Col-avg concentration CH4                     (ppm).
        umn2o        !Col-avg concentration N2O                     (ppm).

      real (kind=dbl_kind), dimension(1,nlm):: &
        pl,&         !Layer pressure                                (hPa).
        dpl,&        !Layer thickness                               (hPa).
        tl ,&        !Temperature                                     (K).
        ql ,&        !Specific humidity                           (kg/kg).
        qcwl ,&      !Cloud water mixing ratio                    (kg/kg).
        ncwl ,&      !Cloud droplet concentration                  (#/kg).
        qcil,&       !Cloud ice mixing ratio                      (kg/kg).
        qrwl,&       !Rain mixing ratio                           (kg/kg).
        qril,&       !Snow mixing ratio                           (kg/kg).
        o3l,&        !Ozone mixing ratio                          (kg/kg).
        acld,&         !Radiative cloud fraction                        (-).
        atl ,&       !All-sky LW radiative heating rate             (K/s).
        asl ,&       !All-sky SW radiative heating rate             (K/s).
        fulwb,&       !All-sky LW upwelling flux                   (W/m^2).
        fdlwb,&       !All-sky LW downwelling flux                 (W/m^2).
        fuswb,&       !All-sky SW upwelling flux                   (W/m^2).
        fdswb         !All-sky SW downwelling flux                 (W/m^2).

      !Note that rain mixing ratio is unused by BUGSrad, but I've put it
      !into the std_profile.dat file for completeness

      real (kind=dbl_kind), dimension(1,nlm+1):: &
        pl2          !Level pressure                                (hPa).

      real (kind=dbl_kind), dimension(1,nlm+1):: &
        fulw,&       !All-sky LW upwelling flux                   (W/m^2).
        fdlw,&       !All-sky LW downwelling flux                 (W/m^2).
        fusw,&       !All-sky SW upwelling flux                   (W/m^2).
        fdsw         !All-sky SW downwelling flux                 (W/m^2).

!-----------------------------------------------------------------------
      !print*,'inside bugs_driver'
      !print*,tl(1,1:5)
      !print*,ql(1,1:5)
      !print*,qcwl(1,1:5)
!---- 1. READ PROFILE DATA FROM FILE:
      
      nlen = 1    !to do timing tests
      len = nlen
!----
! AJD changed: original:
!      do l=2,nlm+1
!       pl2(1,l)=(pl(1,l-1)+pl(1,l))/2.
!      enddo
!      pl2(1,1)=pl2(1,2)+abs(pl2(1,2)-pl2(1,3))
! Bug in original: pl(1,nlm+1) is not defined!
!
! AJD changed: modified:
      do l=2,nlm
       pl2(1,l)=(pl(1,l-1)+pl(1,l))/2.
      enddo
      pl2(1,1)=pl2(1,2)+abs(pl2(1,2)-pl2(1,3))
! AJD - extrapolate to the top-layer pressure:
      pl2(1,nlm+1)=pl2(1,nlm)-abs(pl2(1,nlm-1)-pl2(1,nlm))

      pl = pl/100. !convert from Pascals to millibars
      pl2 = pl2/100. 

      do l=1,nlm
        dpl(1,l) = pl2(1,l)-pl2(1,l+1)
        if (qcwl(1,l)>1.e-14 .or. qrwl(1,l)>1.e-14 &
           .or. qcil(1,l)>1.e-14 .or. qril(1,l)>1.e-14) then
           acld(1,l)=1.
        else
           acld(1,l)=0.
        endif
      enddo

      ts = tl(1,1) !Adele - fix later

      !BUGSRAD works with columns that are upside down.
      CALL flip_profile(nlm,pl)
      CALL flip_profile(nlm+1,pl2)
      CALL flip_profile(nlm,dpl)
      CALL flip_profile(nlm,tl)
      CALL flip_profile(nlm,ql)
      CALL flip_profile(nlm,qcwl)
      CALL flip_profile(nlm,ncwl)
      CALL flip_profile(nlm,qcil)
      CALL flip_profile(nlm,qrwl)
      CALL flip_profile(nlm,qril)
      CALL flip_profile(nlm,o3l)
      CALL flip_profile(nlm,acld)

      alndr=alvdr
      alvdf = alvdr
      alndf = alvdr

      umco2=420.
      umch4=0. !Methane
      umn2o=0. !N2O
      fdsw=0.
      fusw=0.
      fdlw=0.
      fulw=0.

      slr = 1.0

!---- 4. CALL THE RADIATIVE TRANSFER CODE:
      call bugs_rad(nlen,len,nlm,pl2,pl,dpl,tl,ql,qcwl,ncwl,qcil,qril, &
                    o3l,ts,amu0,slr,alvdf,alndf,alvdr,alndr,sol_const, &
                    gravity,cp_dry_air,asl,atl,fdsw,fusw,fdlw,fulw, &
                    acld, umco2, umch4, umn2o)
      CALL flip_profile(nlm,pl)
      CALL flip_profile(nlm,tl)
      CALL flip_profile(nlm,ql)
      CALL flip_profile(nlm,qcwl)
      CALL flip_profile(nlm,ncwl)
      CALL flip_profile(nlm,qcil)
      CALL flip_profile(nlm,qrwl)
      CALL flip_profile(nlm,qril)
      CALL flip_profile(nlm,o3l)
      CALL flip_profile(nlm,asl)
      CALL flip_profile(nlm,atl)
      CALL flip_profile(nlm+1,fusw)
      CALL flip_profile(nlm+1,fdsw)
      CALL flip_profile(nlm+1,fulw)
      CALL flip_profile(nlm+1,fdlw)

      fdswb(1,:)=fdsw(1,1:nlm)
      fuswb(1,:)=fusw(1,1:nlm)
      fdlwb(1,:)=fdlw(1,1:nlm)
      fulwb(1,:)=fulw(1,1:nlm)

      end subroutine bugs_driver

      subroutine flip_profile(nrad,prof)
      use bugs_kinds, only: dbl_kind

      implicit none

      integer::nrad,i
      real(kind=dbl_kind),dimension(1,nrad)::prof,temp

      temp=prof
      do i=1,nrad
         prof(1,i)=temp(1,nrad+1-i)
      enddo

      end subroutine flip_profile
      subroutine flip_profile_int(nrad,prof)
      use bugs_kinds, only: dbl_kind

      implicit none

      integer::nrad,i
      integer,dimension(1,nrad)::prof,temp

      temp=prof
      do i=1,nrad
         prof(1,i)=temp(1,nrad+1-i)
      enddo

      end subroutine flip_profile_int
!-----------------------------------------------------------------------
