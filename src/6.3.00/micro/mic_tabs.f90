!##############################################################################
Subroutine haznuc ()

use micphys

implicit none

integer :: ithz,irhhz,k
real :: denccn,gnuccn,dnccn,ddccn,rhhz,c1hz,c2hz,c3hz,bhz,dm,sum  &
       ,dccn,y,dum,thz
real, external :: gammln

!  Haze nucleation table

!Density of ammonium sulfate (g/cm3)
denccn = 1.769
!Saleeby(02-21-2007) Shape parameter; Originally = 1.
gnuccn = 2.
!Saleeby(02-21-2007) Originally = .075e-4 (cm)
!Convert median radius (m) to (cm) and to diameter
dnccn = 2. * aero_medrad(1) * 1.e2
!Bin increment (cm)
ddccn = .005e-4

do ithz = 1,nthz
   thz = -60. + dthz * float(ithz - 1)
   do irhhz = 1,nrhhz
      rhhz = 0.82 + drhhz * float(irhhz - 1)
      c1hz = (3.14159 * denccn / 6.) ** (-.333333)
      c2hz = -14.65 - 1.045 * thz
      c3hz = -492.35 - 8.34 * thz - 0.0608 * thz ** 2
      bhz = min(38., max(-38., c2hz + c3hz * (1. - rhhz)))
      dm = c1hz * 10 ** (-bhz/6.)

      sum = 0.
      dccn = 0.
      do k=1,200
         dccn = dccn + ddccn
         y=dccn / dnccn
         dum=min(50., (dccn / dm) ** 6)
         sum = sum + y ** (gnuccn-1.) * exp(-y) * (1. - exp(-dum))
      enddo
      frachz(irhhz,ithz) = sum*ddccn/(exp(gammln(gnuccn))*dnccn)
   enddo
enddo

return
END SUBROUTINE haznuc

!##############################################################################
Subroutine homfrzcl (dtlt,ngr)

use micphys

implicit none

integer :: itc,ngr,k,idnc
real :: ddc,ajlso,dnc,sum,dc,v1,tc,y,dtlt
real, external :: gammln

!  Make table for homogeneous freezing of cloud droplets
!  Uses characteristic diameter instead of true diameter

ddc = 0.5e-6
!ndnc = 11     | both preset in micphys.h
!ddnc = 2.e-6  |

do itc = 1,ntc
   tc = -50. + dtc * float(itc-1)
   y = -(606.3952+tc*(52.6611+tc*(1.7439+tc*(.0265+tc*1.536e-4))))
   ajlso = 1.e6 * 10. ** y

   do idnc = 1,ndnc
      dnc = ddnc * float(idnc)
      sum = 0.
      dc = 0.
      do k = 1,2000
         dc = dc + ddc
         v1 = 0.523599 * dc ** 3
         sum = sum + (dc / dnc) ** (gnu(1) - 1.) * exp(-dc / dnc)  &
            * (1. - exp(-ajlso * v1 * dtlt))
      enddo
      fracc(idnc,itc,ngr) = sum * ddc / (exp(gammln(gnu(1))) * dnc)
   enddo
 
enddo

return
END SUBROUTINE homfrzcl

!##############################################################################
Subroutine mksedim_tab (m1,ngr,zm,dzt,pcpfillc,pcpfillr,sfcpcp,allpcp,dtsed)

! As before, sedimentation is not yet designed to transfer hydrometeor mass
! between grids in the case where a nested grid does not reach the top and/or
! bottom of the model domain.  Thus, keep vertical nested grid boundaries away
! from where sedimentation occurs.

use micphys

implicit none

integer, parameter :: nbin=50
integer :: m1,iembs,lcat,lhcat,k,kkf,ibin,kk,jbin,ngr,idensrtgt,iband
integer, dimension(16) :: lcat0

real, dimension(16) :: dmbsed

real :: dmbodn,diam0,diam1,fac1,fac3,sumc,sumr,diam,fac2,fac4  &
       ,disp,ztopnew,zbotnew,fallin,delzsfc,dispemb,dispmax,dispmx
real :: dtsed,dispmax0,dispmax1
real, dimension(m1) :: zm,dzt
real, dimension(m1,maxkfall,nembfall,nhcat,ndensrtgt,nband) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat,ndensrtgt,nband) :: sfcpcp
real, dimension(m1,nembfall,nhcat,ndensrtgt,nband) :: allpcp
real, dimension(nbin) :: cbin,rbin,reldisp
real, external :: gammln
real, external :: gammp

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/

data dmbsed/40.e-6,7.e-3,750.e-6,10.e-3,15.e-3,15.e-3,30.e-3,750.e-6 &
   ,750.e-6,750.e-6,750.e-6,10.e-3,10.e-3,10.e-3,10.e-3,150.e-6/

! Because timestep may now be variable in time, define sedtime0 and sedtime1
! here as 0.1 seconds and 3000 seconds.  The former is supposed to be
! less than 0.7 of the shortest timestep on any grid (sqrt(dn0i) never exceeds
! 0.7) and the latter is the longest timestep expected to ever be used (300
! seconds) times a factor of 2 for the maximum inverse of rtgt times a factor
! of 5 for the largest value of sqrt(dn0i).

sedtime0 = .1
sedtime1 = 3000.
dispmax0 = 500.

! Loop over hydrometeor categories

do lhcat = 1,nhcat
   lcat = lcat0(lhcat)

   !Displacement max based on largest size and timestep (Saleeby 8-17-05)
   dispmax1 = dtsed * cfvt(lhcat) * dmbsed(lhcat) ** pwvt(lhcat)
   dispmax = min(dispmax0,dispmax1)

   dispemb0(lhcat,ngr) = sedtime0 * cfvt(lhcat)  &
      * (emb0(lcat) / cfmas(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

   dispemb1(lhcat,ngr) = sedtime1 * cfvt(lhcat)  &
      * (emb1(lcat) / cfmas(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

   !Limit dispemb1 to a maximum of dispmax
   if (dispemb1(lhcat,ngr) .gt. dispmax) dispemb1(lhcat,ngr) = dispmax

! Loop over bins, filling them with fractional number, fractional mass,
! and displacement quotient relative to emb.

   dmbodn = (exp(gammln(gnu(lcat)+pwmas(lhcat)) &
           -gammln(gnu(lcat))))**pwmasi(lhcat)
   diam0 = 0.06 * dmbodn
   diam1 = 1.0 * dmbodn
   fac1 = gammp(gnu(lcat),diam0)
   fac3 = gammp(gnu(lcat) + pwmas(lhcat),diam0)
   sumc = 0.
   sumr = 0.

   do jbin = 1,nbin

      diam = diam0 * (diam1 / diam0) ** (float(jbin)/float(nbin))
      fac2 = gammp(gnu(lcat),diam)
      fac4 = gammp(gnu(lcat) + pwmas(lhcat),diam)
      cbin(jbin) = fac2 - fac1
      rbin(jbin) = fac4 - fac3
      fac1 = fac2
      fac3 = fac4
      sumc = sumc + cbin(jbin)
      sumr = sumr + rbin(jbin)
      reldisp(jbin) = diam ** pwvt(lhcat)

   enddo

   do jbin = 1,nbin
      cbin(jbin) = cbin(jbin) / sumc
      rbin(jbin) = rbin(jbin) / sumr
   enddo

! Loop over displacement distance for size emb.

   do iembs = 1,nembfall
      dispemb = dispemb0(lhcat,ngr)  &
         * (dispemb1(lhcat,ngr) / dispemb0(lhcat,ngr))  &
         ** (float(iembs-1) / float(nembfall-1))

! Zero out concentration and mass fill arrays and surface precip array
! before accumulation.
      do iband=1,nband
       do idensrtgt=1,ndensrtgt
        do k = 1,m1
          do kkf = 1,maxkfall
             pcpfillc(k,kkf,iembs,lhcat,idensrtgt,iband) = 0.
             pcpfillr(k,kkf,iembs,lhcat,idensrtgt,iband) = 0.
          enddo
          allpcp(k,iembs,lhcat,idensrtgt,iband) = 0.
          if (k .le. maxkfall) sfcpcp(k,iembs,lhcat,idensrtgt,iband) = 0.
        enddo
       enddo
      enddo

! Loop over vertical grid index.

      do k = 2,m1-1

!Bob (10/24/00):  Limit disp distance to (maxkfall-1) levels

         dispmx = dispmax
         if (k .gt. maxkfall) then
            dispmx = min(dispmx,zm(k-1) - zm(k-maxkfall))
         endif

! Loop over bins

         do ibin = 1,nbin
            disp = dispemb * reldisp(ibin)
            if (disp .gt. dispmx) disp = dispmx
            ztopnew = zm(k) - disp
            zbotnew = zm(k-1) - disp

! Loop over grid cells that a parcel falls into.

            do kkf = 1,min(k-1,maxkfall)

               kk = k + 1 - kkf
               if (zbotnew .gt. zm(kk)) go to 50

               if (ztopnew .le. zm(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzt(kk) *  &
                     (min(zm(kk),ztopnew) - max(zm(kk-1),zbotnew))
               endif

               pcpfillc(k,kkf,iembs,lhcat,1,1) = pcpfillc(k,kkf,iembs,lhcat,1,1)  &
                  + fallin * cbin(ibin)

               pcpfillr(k,kkf,iembs,lhcat,1,1) = pcpfillr(k,kkf,iembs,lhcat,1,1)  &
                  + fallin * rbin(ibin)

            enddo
50          continue

! Compute precip rate at levels above ground
            allpcp(k,iembs,lhcat,1,1) = allpcp(k,iembs,lhcat,1,1) + disp * rbin(ibin)

! Compute surface precipitation.
            if (zbotnew .lt. 0.) then
               delzsfc = min(0.,ztopnew) - zbotnew
               if (k .le. maxkfall) sfcpcp(k,iembs,lhcat,1,1)  &
                  = sfcpcp(k,iembs,lhcat,1,1) + delzsfc * rbin(ibin)
            endif

         enddo

      enddo
   enddo
enddo

return
END SUBROUTINE mksedim_tab

!##############################################################################
Subroutine mksedim_tab_trubin (m1,zm,dzt,pcpfillc,pcpfillr,sfcpcp,allpcp,dtsed)

! As before, sedimentation is not yet designed to transfer hydrometeor mass
! between grids in the case where a nested grid does not reach the top and/or
! bottom of the model domain.  Thus, keep vertical nested grid boundaries away
! from where sedimentation occurs.

use micphys

implicit none

integer, parameter :: nbin=41
integer :: m1,iembs,lcat,lhcat,k,kkf,ibin,kk,idensrtgt,iband
integer, dimension(16) :: lcat0

real :: dmbodn,fac1,sumc,sumr,fac2  &
       ,disp,ztopnew,zbotnew,fallin,delzsfc,dispmax,dispmx,dtsed
real, dimension(m1) :: zm,dzt
real, dimension(m1,maxkfall,nembfall,nhcat,ndensrtgt,nband) :: pcpfillc,pcpfillr
real, dimension(maxkfall,nembfall,nhcat,ndensrtgt,nband) :: sfcpcp
real, dimension(m1,nembfall,nhcat,ndensrtgt,nband) :: allpcp
real, dimension(nbin) :: cbin,rbin,fall
real, dimension(nbin+1) :: diam,massx
real :: tmpa,hemb,meand,dn1,pi,densrtgt
real, external :: gammln
real, external :: gammp

data lcat0 /1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/
pi=3.141592654

! NEED TO COMPLETELY REDO SEDIMENTATION TABLES IF YOU USE VARIABLE TIMESTEP

! Loop over three size bands for Vt power laws
do iband=1,nband

! Loop over hydrometeor categories
do lhcat = 1,nhcat
   lcat = lcat0(lhcat)

   pwmasi(lhcat) = 1. / pwmas(lhcat)
   dmbodn = ( exp(gammln(gnu(lcat) + pwmas(lhcat))  & !ratio of diam assoc w/ mean mass
        -gammln(gnu(lcat))) ) ** pwmasi(lhcat)        !to char. diam [unitless]

   !Calculating the mass categories (m.k.s) with the lowest diameter
   !of 3.125 microns and mass doubling every bin
   diam(1)=1.5625 * 2.0 * 1.0e-6
   massx(1)=cfmas(lhcat)*diam(1)**pwmas(lhcat)
   do ibin=2,nbin+1
     massx(ibin)=2.0*massx(ibin-1)
     diam(ibin)=(massx(ibin) / cfmas(lhcat)) ** pwmasi(lhcat)
   enddo

! Loop over displacement distance for size emb.
   do iembs = 1,nembfall
      sumc  = 0.
      sumr  = 0.
      tmpa  = emb1(lcat)/emb0(lcat)
      hemb  = emb0(lcat)*(tmpa**(float(iembs-1)/float(nembfall-1)))
      meand = (hemb/cfmas(lhcat))**pwmasi(lhcat) !diam assoc w/ mean mass (NOT MEAN DIAM!)
      dn1   = meand/dmbodn        !char diam for current mean mass-gnu pair

      !Determine gamma distribution values for size bin
      do ibin=1,nbin
       cbin(ibin) = 0.
       rbin(ibin) = 0.
       fac1 = gammp(gnu(lcat),diam(ibin)/dn1)
       fac2 = gammp(gnu(lcat),diam(ibin+1)/dn1)
       cbin(ibin) = fac2-fac1
       sumc = sumc + cbin(ibin)
       fac1 = gammp(gnu(lcat)+3.,(diam(ibin)/dn1))
       fac2 = gammp(gnu(lcat)+3.,(diam(ibin+1)/dn1))
       rbin(ibin) = fac2-fac1
       sumr = sumr + rbin(ibin)
      enddo

      !Truncate top and bottom 5% of gamma distribution to prevent
      !tails from dominating the sedimentation response
      sumc=0.0
      sumr=0.0
      do ibin=1,nbin
        if(rbin(ibin)<0.05 .and. cbin(ibin)<0.05) then
          cbin(ibin) = 0.0
          rbin(ibin) = 0.0
        endif
        sumc    = sumc + cbin(ibin)
        sumr    = sumr + rbin(ibin)
      enddo

      !Normalize bin moment values
      do ibin=1,nbin
        cbin(ibin) = cbin(ibin) / sumc
        rbin(ibin) = rbin(ibin) / sumr
      enddo

! Loop over density/rtgt range to get adjusted fall speed 0.8 to 3.65

      do idensrtgt=1,ndensrtgt
       densrtgt = idensrtgt / 10.0

! Zero out concentration and mass fill arrays and surface precip array
! before accumulation.

       do k = 1,m1
          do kkf = 1,maxkfall
             pcpfillc(k,kkf,iembs,lhcat,idensrtgt,iband) = 0.
             pcpfillr(k,kkf,iembs,lhcat,idensrtgt,iband) = 0.
          enddo
          allpcp(k,iembs,lhcat,idensrtgt,iband) = 0.
          if (k .le. maxkfall) sfcpcp(k,iembs,lhcat,idensrtgt,iband) = 0.
       enddo

!Displacement max based on largest size and timestep (Saleeby 8-17-05)

       dispmax=0.
       do ibin=1,nbin
         !distance fallen in 1 sec
         if(iplaws==0 .or. iplaws==1) &
          fall(ibin) = cfvt(lhcat) * (diam(ibin)**pwvt(lhcat))
         if(iplaws==2) &
          fall(ibin) = bcfvt(iband,lhcat) * (diam(ibin)**bpwvt(iband,lhcat)) 
         if(lhcat==2)  fall(ibin) = min(13.0,fall(ibin)) !limit rain to 13m/s max
         if(lhcat==4)  fall(ibin) = min( 6.0,fall(ibin)) !limit snow to 6m/s max
         if(lhcat==14) fall(ibin) = min( 6.0,fall(ibin)) !limit snow ndl to 6m/s max
         if(lhcat==15) fall(ibin) = min( 6.0,fall(ibin)) !limit snow ros to 6m/s max
         fall(ibin) = fall(ibin) * densrtgt !adjust for air density and rtgt
         if(rbin(ibin)>0.) dispmax = fall(ibin) * dtsed !max timestep displacement
       enddo

! Loop over vertical grid index.

       do k = 2,m1-1
         dispmx = dispmax
         if (k .gt. maxkfall) then
            dispmx = min(dispmx,zm(k-1) - zm(k-maxkfall))
         endif

! Loop over bins

         do ibin = 1,nbin
            disp = dtsed * fall(ibin) !distance fallen in 1 timestep
            if (disp .gt. dispmx) disp = dispmx
            ztopnew = zm(k) - disp
            zbotnew = zm(k-1) - disp !drops @ z(k) is over layer z(k-1) to z(k)

! Loop over grid cells that a parcel falls into.

            do kkf = 1,min(k-1,maxkfall)

               kk = k + 1 - kkf
               if (zbotnew .gt. zm(kk)) go to 50

               if (ztopnew .le. zm(kk-1)) then
                  fallin = 0.
               else
                  fallin = dzt(kk) *  &
                     (min(zm(kk),ztopnew) - max(zm(kk-1),zbotnew))
               endif

               pcpfillc(k,kkf,iembs,lhcat,idensrtgt,iband) = &
                    pcpfillc(k,kkf,iembs,lhcat,idensrtgt,iband)  &
                  + fallin * cbin(ibin)

               pcpfillr(k,kkf,iembs,lhcat,idensrtgt,iband) = &
                    pcpfillr(k,kkf,iembs,lhcat,idensrtgt,iband)  &
                  + fallin * rbin(ibin)
  
            enddo
50          continue

! Compute precip rate at levels above ground
            allpcp(k,iembs,lhcat,idensrtgt,iband) = & 
                  allpcp(k,iembs,lhcat,idensrtgt,iband) &
                + disp * rbin(ibin)

! Compute surface precipitation.
            if (zbotnew .lt. 0.) then
               delzsfc = min(0.,ztopnew) - zbotnew
               if (k .le. maxkfall) sfcpcp(k,iembs,lhcat,idensrtgt,iband)  &
                  = sfcpcp(k,iembs,lhcat,idensrtgt,iband) + delzsfc * rbin(ibin)
            endif

         enddo !ibin
       enddo !k-index
      enddo !idensrtgt
   enddo !iembs
enddo !lhcat
enddo !iband

return
END SUBROUTINE mksedim_tab_trubin

!##############################################################################
Subroutine tabmelt ()

use micphys

implicit none

integer, parameter :: nbins=500

integer :: lhcat,lcat,ndns1,ibin,inc,iter,idns
integer, dimension(16) :: lcat0
real :: dn,gammaa,totfmg,totmass,vtx,fre,totqm,qmgoal,qmnow,totmdqdt,deltat  &
       ,pliqmass,picemass,critmass,vk

real, dimension(nbins) :: db,fmg,pmass,binmass,dqdt,q
real, dimension(ncat) :: dmean
real, external :: gammln
data lcat0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/
data dmean/20.e-6,500.e-6,30.e-6,500.e-6,500.e-6,500.e-6,8000.e-6,60.e-6/
data vk/0.2123e-04/

do lhcat = 1,nhcat
   lcat = lcat0(lhcat)
   dn = dmean(lcat) / gnu(lcat)
   gammaa = exp(gammln(gnu(lcat)))

   rmlttab(1,lhcat) = 0.0
   rmlttab(ninc,lhcat) = 1.0
   enmlttab(1,lhcat) = 0.0
   enmlttab(ninc,lhcat) = 1.0

   ndns1 = 1
   if (lcat .eq. 7) ndns1 = ndns

   do idns = 1,ndns1
      shedtab(1,idns) = 0.0
      shedtab(ninc,idns) = 0.0

      if (ndns1 .gt. 1) dn = 1.0e-3 * float(idns) / gnu(lcat)

      totfmg = 0.
      totmass = 0.
      do ibin = 1,nbins
         db(ibin) = 0.02 * dn * (float(ibin) - 0.5)
         fmg(ibin) = (db(ibin) / dn) ** (gnu(lcat) - 1.)  &
            / (dn * gammaa) * exp(-db(ibin) / dn)
         totfmg = totfmg + fmg(ibin)
         q(ibin) = 0.
         pmass(ibin) = cfmas(lhcat) * db(ibin) ** pwmas(lhcat)
         binmass(ibin) = pmass(ibin) * fmg(ibin)
         totmass = totmass + binmass(ibin)
         vtx = cfvt(lhcat) * db(ibin) ** pwvt(lhcat)
         fre = (1.0 + 0.229 * sqrt(vtx * db(ibin) / vk))  &
            * shape(lhcat)
         dqdt(ibin) = db(ibin) ** (1. - pwmas(lhcat)) * fre
      enddo
      totqm = totmass * 80.

      do inc = 2,ninc-1
         qmgoal = totqm * float(inc-1) / float(ninc-1)

         do iter = 1,2
            qmnow = 0.
            totmdqdt = 0.
            do ibin = 1,nbins
               if(q(ibin) .lt. 79.9999)then
                  totmdqdt = totmdqdt + binmass(ibin) * dqdt(ibin)
               endif
               qmnow = qmnow + q(ibin) * binmass(ibin)
            enddo
            deltat = max(0.,(qmgoal - qmnow) / totmdqdt)
            do ibin = 1,nbins
               q(ibin) = min(80.,q(ibin) + dqdt(ibin) * deltat)
            enddo
         enddo

!  For the current inc value (representing total liquid fraction), compute
!  melted mixing ratio (rmlttab) and number (enmlttab) from totally-melted
!  bins and compute shedded mixing ratio (shedtab) from partially-melted bins.

         if(idns .eq. 7 .or. ndns1 .eq. 1)then
            rmlttab(inc,lhcat) = 0.0
            do ibin = 1,nbins
               if(q(ibin) .gt. 79.9)then
                  rmlttab(inc,lhcat) = rmlttab(inc,lhcat) + binmass(ibin)
               endif
            enddo
            rmlttab(inc,lhcat) = rmlttab(inc,lhcat) / totmass
         endif

         if(idns .eq. 7 .or. ndns1 .eq. 1)then
            enmlttab(inc,lhcat) = 0.0
            do ibin = 1,nbins
               if(q(ibin) .gt. 79.9)then
                  enmlttab(inc,lhcat) = enmlttab(inc,lhcat)  &
                     + fmg(ibin)
               endif
            enddo
            enmlttab(inc,lhcat) = enmlttab(inc,lhcat) / totfmg
         endif

         if(lcat .eq. 7)then
            shedtab(inc,idns) = 0.0
            do ibin = 1,nbins
               if(q(ibin) .le. 79.9)then
                  pliqmass = pmass(ibin) * q(ibin) / 80.
                  picemass = pmass(ibin) - pliqmass
                  critmass = .268e-3 + .1389 * picemass
                  shedtab(inc,idns) = shedtab(inc,idns)  &
                     + max(0.0, pliqmass - critmass) * fmg(ibin)
               endif
            enddo
            shedtab(inc,idns) = shedtab(inc,idns) / totmass
         endif

      enddo
   enddo
enddo

return
END SUBROUTINE tabmelt

!##############################################################################
Subroutine mkcoltb ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

integer, parameter :: ndx=20
integer :: ihx,ix,ihy,iy,iemby,iembx,idx
integer, dimension(16) :: ix0,iy0
real :: gxm,dnminx,dnmaxx,dxlo,dxhi,gyn,gyn1,gyn2,gynp,gynp1,gynp2,gym  &
       ,dnminy,dnmaxy,dny,vny,dnx,ans
real, external :: gammln
real, external :: xj_col

real, dimension(ndx) :: dx,fx,gx
data ix0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/
data iy0/1,2,3,4,5,6,7,3,3,3,3,4,4,4,4,8/

if(iprntstmt>=1 .and. print_msg) &
 print*,'Making the 2-Moment Microphysics Collection Tables'

do ihx = 1,nhcat

   ix = ix0(ihx)

   gxm = exp(gammln(gnu(ix))-gammln(gnu(ix) + pwmas(ihx)))
   dnminx = ((emb0(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
   dnmaxx = ((emb1(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
   dxlo = .01 * dnminx
   dxhi = 10. * dnmaxx

do ihy = 1,nhcat

   !print*, 'ihx,ihy,ipairc,ipairr',ihx,ihy,ipairc(ihx,ihy)  &
   !   ,ipairr(ihx,ihy)

   iy = iy0(ihy)

   if (ipairc(ihx,ihy) .gt. 0 .or. ipairr(ihx,ihy) .gt. 0) then
      gyn   = exp(gammln(gnu(iy)))
      gyn1  = exp(gammln(gnu(iy) + 1.)) / gyn
      gyn2  = exp(gammln(gnu(iy) + 2.)) / gyn
      gynp  = exp(gammln(gnu(iy) + pwvt(ihy))) / gyn
      gynp1 = exp(gammln(gnu(iy) + pwvt(ihy) + 1.)) / gyn
      gynp2 = exp(gammln(gnu(iy) + pwvt(ihy) + 2.)) / gyn

      gym = exp(gammln(gnu(iy))-gammln(gnu(iy) + pwmas(ihy)))
      dnminy = ((emb0(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))
      dnmaxy = ((emb1(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))

      do iemby = 1,nembc
         dny = dnminy * (dnmaxy / dnminy) ** (float(iemby-1)  &
            / float(nembc-1))
         vny = cfvt(ihy) * dny ** pwvt(ihy)
         do iembx = 1,nembc

            dnx = dnminx * (dnmaxx / dnminx) ** (float(iembx-1)  &
               / float(nembc-1))
            do idx = 1,ndx
               dx(idx) = dxlo * (dxhi / dxlo)  &
                  ** (float(idx-1) / float(ndx-1))
               fx(idx) = xj_col(dx(idx),cfvt(ihx),pwvt(ihx),cfvt(ihy)  &
                  ,pwvt(ihy),vny,dnx,dny,gnu(ix),gnu(iy)  &
                  ,gyn1,gyn2,gynp,gynp1,gynp2)
               gx(idx) = fx(idx) * cfmas(ihx)  &
                  * dx(idx) ** pwmas(ihx)

            enddo
            if (ipairc(ihx,ihy) .gt. 0) then
               CALL avint (dx,fx,ndx,dxlo,dxhi,ans)
               coltabc(iembx,iemby,ipairc(ihx,ihy))=  &
                  -log10(max(1.0e-30,ans))
            endif
            if (ipairr(ihx,ihy) .gt. 0) then
               CALL avint (dx,gx,ndx,dxlo,dxhi,ans)
               coltabr(iembx,iemby,ipairr(ihx,ihy))=  &
                  -log10(max(1.0e-30,ans))
            endif
         enddo
      enddo
   endif
enddo
enddo

return
END SUBROUTINE mkcoltb

!##############################################################################
real Function xj_col (dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu  &
                 ,gyn1,gyn2,gynp,gynp1,gynp2)

implicit none

real :: dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu,gyn1,gyn2,gynp,gynp1,gynp2  &
       ,dnxi,rdx,vx,dxy,ynup
real, external :: gammln
real, external :: gammp
real, external :: gammq

dnxi = 1. / dnx
rdx = dx * dnxi
vx = cvx * dx ** pvx
dxy = (vx / cvy) ** (1. / pvy) / dny
ynup = ynu + pvy

if (rdx .lt. 38.) then
   xj_col = exp(-rdx-gammln(xnu)-gammln(ynu))*rdx**(xnu-1.)*dnxi*(  &
       vx*(dx*dx*(gammp(ynu,dxy)-gammq(ynu,dxy))  &
         +2.*dx*dny*gyn1*(gammp(ynu+1.,dxy)-gammq(ynu+1.,dxy))  &
         +dny*dny*gyn2*(gammp(ynu+2.,dxy)-gammq(ynu+2.,dxy)))  &
     -vny*(dx*dx*gynp*(gammp(ynup,dxy)-gammq(ynup,dxy))  &
         +2.*dx*dny*gynp1*(gammp(ynup+1.,dxy)-gammq(ynup+1.,dxy))  &
         +dny*dny*gynp2*(gammp(ynup+2.,dxy)-gammq(ynup+2.,dxy))))
else
   xj_col = 0.
endif

return
END FUNCTION xj_col

!##############################################################################
Subroutine make_autotab ()

use micphys

implicit none

integer, parameter :: ithresh=14,ithresh1=17,ibins=36,icutoff=15
integer :: i,k,idcc,id1cd,id2cd,idccr,irrcr,idrcr,lcat,cld,maxcld
real :: r1,r2,r3,ri,en1,en2,en3,enice,en1i,en1i2,d1,d2,d3,di &
       ,sum1,sum10,sun1,sun10,sum2,sum20,sun2,sun20,sum3,sum30,sun3,sun30
real :: collectormass

real, dimension(ibins+1) :: x,diam,xcs,diamcs,xca,diamca,xcg,diamcg,xch,diamch
real, dimension(ibins+1,4) :: xi,diami
real, dimension(ibins) :: ank0,amk0,ank,amk,ank1,amk1,ank2,amk2,ank3,amk3
real, dimension(ibins,ibins,6) :: akbarx
real, dimension(ibins,ibins,4) :: akbarci,akbarxci
real, dimension(ibins,ibins) :: akbar,akbarcs,akbarca,akbarcg,akbarch

! This routine works in cgs units.
! read in mass grid x(k+1)=2*x(k), diameters (diam) and collection kernel kbar
! GNF: kernels for ice with cloud?

CALL data (x,diam,akbar,ibins)
CALL data_cs (xcs,diamcs,akbarcs,ibins,cfmas(4),pwmas(4))
CALL data_ca (xca,diamca,akbarca,ibins,cfmas(5),pwmas(5))
CALL data_cg (xcg,diamcg,akbarcg,ibins,cfmas(6),pwmas(6))
CALL data_ch (xch,diamch,akbarch,ibins,cfmas(7),pwmas(7))

CALL azero (ibins*ibins*6,akbarx)
CALL azero (ibins*ibins*4,akbarci)
CALL azero (ibins*ibins*4,akbarxci)

do i=1,ibins+1
 xi(i,1)=xcs(i)
 xi(i,2)=xca(i)
 xi(i,3)=xcg(i)
 xi(i,4)=xch(i)
 diami(i,1)=diamcs(i)
 diami(i,2)=diamca(i)
 diami(i,3)=diamcg(i)
 diami(i,4)=diamch(i)
enddo
do i=1,ibins
 do k=1,ibins
  akbarci(i,k,1)=akbarcs(i,k)
  akbarci(i,k,2)=akbarca(i,k)
  akbarci(i,k,3)=akbarcg(i,k)
  akbarci(i,k,4)=akbarch(i,k)
 enddo
enddo

!Break up kernel in 6 different segments for cloud,drizzle,rain
do i=1,ibins
   do k=1,ibins
      !For CLOUD-CLOUD collection for 1 and 2 mode options
      if(i <  ithresh  .and. k <  ithresh  .and. jnmb(8).gt.0) &
         akbarx(i,k,1) = akbar(i,k)
      if(i <=  icutoff  .and. k <=  icutoff  .and. jnmb(8).eq.0) &
         akbarx(i,k,1) = akbar(i,k)

      !For RAIN-CLOUD collection
      if(i >  ithresh1 .and. k < ithresh  .and. jnmb(8).gt.0) &
         akbarx(i,k,3) = akbar(i,k)
      if(i >  icutoff  .and. k <= icutoff  .and. jnmb(8).eq.0) &
         akbarx(i,k,3) = akbar(i,k)

      if(jnmb(8).gt.0) then  !Dual-cloud mode option only
       !For DRIZZLE-CLOUD collection
       if(i >= ithresh  .and. i <= ithresh1 .and. k < ithresh) &
          akbarx(i,k,2) = akbar(i,k)

       !For DRIZZLE-DRIZZLE collection
       if(i<=ithresh1 .and. i>=ithresh .and. k<=ithresh1 .and. k>=ithresh) &
          akbarx(i,k,5) = akbar(i,k)

       !For RAIN-DRIZZLEcollection
       if(i >  ithresh1 .and. k <= ithresh1 .and. k >= ithresh) &
          akbarx(i,k,4) = akbar(i,k)
      endif

      !For RAIN-RAIN collection (not currently used)
      !This is done in subroutine "cols" in order to include drop breakup
      !if(i >  ithresh1 .and. k >  ithresh1) &
      !   akbarx(i,k,6) = akbar(i,k)
   enddo
enddo

!Break up kernel in 2 different segments for c1,c2,ice species
do lcat=4,7
 do i=1,ibins
   do k=1,ibins
      if(i > icutoff .and. k <= icutoff) &
         akbarxci(i,k,lcat-3) = akbarci(i,k,lcat-3)
   enddo
 enddo
enddo

!(d1min and d1max are equivalent to dmb0 and dmb1, but may have diff values)
!Diameters in cm (4.e-4 cm = 4 microns)
!Mixing ratios in g/cm3 (.01e-6g/cm3 = .01e-3kg/kg = .01 g/kg for dn0=1)
!Mixing ratios in g/cm3 (.01e-8g/cm3 = .01e-5kg/kg = .0001 g/kg for dn0=1)
d1min = 4.e-4
if(jnmb(8).eq.0) d1max = 40.e-4
if(jnmb(8).gt.0) d1max = 35.e-4
d2min = 65.e-4
d2max = 100.e-4
d3min = 2.e-2
d3max = 1.
r3min = .01e-8
r3max = 20.e-6
!For ICE SPECIES: SNOW, AGGREGATES, GRAUPEL, AND HAIL
dimin = 2.e-2
dimax = 1.
rimin = .01e-8
rimax = 20.e-6

d1ecr = log10 (d1max / d1min) / float(ndccr-1)
d2ecr = log10 (d2max / d2min) / float(ndccr-1)
r3ecr = log10 (r3max / r3min) / float(nrrcr-1)

!For SNOW, AGGREGATES, GRAUPEL, AND HAIL
rieci = log10 (rimax / rimin) / float(nrrcr-1)

en1 = 100.
en2 = 1.
en3 = 1.0e-6
!For SNOW, AGGREGATES, GRAUPEL, AND HAIL
enice = 1.0e-6

en1i = 1. / en1
en1i2 = en1i ** 2.

!**************************************************************************
!**************************** CLOUD-CLOUD *********************************
!**************************************************************************
! Start 1 cc loop for r(1) and c(1,2) [r1 approx = (-r2)]
d2 = 65.e-4
r2 = en2 * .5236 * d2 ** 3
r3 = .01e-6

do idcc = 1,ndcc
   d1 = d1min + (d1max - d1min) * float(idcc-1) / float(ndcc-1)
   r1 = en1 * .5236 * d1 ** 3

   if(jnmb(8).gt.0) then
    CALL initg2mode (r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam &
     ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
    CALL sumn (ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
    CALL sumn (ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
    CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,1))
    CALL sumn (ank,amk,1,ithresh-1,ibins,sun1,sum1)
    CALL sumn (ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
   endif
   if(jnmb(8).eq.0) then
    CALL initg1mode (r1,r3,en1,en3,gnu(1),gnu(2),diam &
     ,amk0,ank0,ank1,amk1,ank3,amk3,icutoff,ibins,2)
    CALL sumn (ank0,amk0,1,icutoff,ibins,sun10,sum10)
    CALL sumn (ank0,amk0,icutoff+1,ibins,ibins,sun20,sum20)
    CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,1))
    CALL sumn (ank,amk,1,icutoff,ibins,sun1,sum1)
    CALL sumn (ank,amk,icutoff+1,ibins,ibins,sun2,sum2)
   endif

   r1tabcc(idcc) = max(0.,max((sum10-sum1)*en1i2,(sum2-sum20)*en1i2))
   c1tabcc(idcc) = max(0.,(sun10-sun1)*en1i2)
   c2tabcc(idcc) = max(0.,(sun2-sun20)*en1i2)
   if(r1tabcc(idcc)==0.0 .or. c1tabcc(idcc)==0.0 .or. c2tabcc(idcc)==0.0) then
     r1tabcc(idcc) = 0.
     c1tabcc(idcc) = 0.
     c2tabcc(idcc) = 0.
   endif

!   write(*,'(a2,2X,f7.1,2X,i2,2X,e10.3,2X,e10.3,2X,e10.3)') &
!     'CC',d1*10000,idcc,r1tabcc(idcc),c1tabcc(idcc),c2tabcc(idcc)
enddo

!**************************************************************************
!**************************** RAIN-CLOUD **********************************
!**************************************************************************
! Start 3 cr loops for r(1) c(1,3)  [r1 approx = (-r3)]
d2 = 65.e-4
r2 = en2 * .5236 * d2 ** 3

do idccr = 1,ndccr
   d1 = d1min * 10. ** (d1ecr * float(idccr-1))
   r1 = en1 * .5236 * d1 ** 3

   do irrcr = 1,nrrcr
      r3 = r3min * 10. ** (r3ecr * float(irrcr-1))
      d3minx = max(d3min,(r3 / (.1 * .5236)) ** .333333)
      d3ecr = alog10(d3max / d3minx) / float(ndrcr-1)

      do idrcr = 1,ndrcr
       d3 = d3minx * 10. ** (d3ecr * float(idrcr-1))
       en3 = r3 / (.5236 * d3 ** 3)

       if(jnmb(8).gt.0) then
        CALL initg2mode (r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam &
         ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
        CALL sumn (ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
        CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,3))
        CALL sumn (ank,amk,1,ithresh-1,ibins,sun1,sum1)
       endif
       if(jnmb(8).eq.0) then
        CALL initg1mode (r1,r3,en1,en3,gnu(1),gnu(2),diam &
         ,amk0,ank0,ank1,amk1,ank3,amk3,icutoff,ibins,2)
        CALL sumn (ank0,amk0,1,icutoff,ibins,sun10,sum10)
        CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,3))
        CALL sumn (ank,amk,1,icutoff,ibins,sun1,sum1)
       endif

       r1tabcr(idccr,irrcr,idrcr) = alog10(max(1.0e-20,(sum10-sum1)*en1i))
       c1tabcr(idccr,irrcr,idrcr) = alog10(max(1.0e-20,(sun10-sun1)*en1i))

!      write(*,'(a2,2X,i2,2X,i2,2X,i2,2X,f7.1,2X,f10.2,2X,f10.2)') &
!        'CR',idccr,irrcr,idrcr,d1*10000,r1tabcr(idccr,irrcr,idrcr) &
!        ,c1tabcr(idccr,irrcr,idrcr)
     enddo
   enddo
enddo

if(jnmb(8).gt.0) then
!**************************************************************************
!*************************** DRIZZLE-CLOUD ********************************
!**************************************************************************
r3 = .01e-6

do id1cd = 1,ndcd
   d1 = d1min + (d1max - d1min) * float(id1cd-1) / float(ndcd-1)
   r1 = en1 * .5236 * d1 ** 3

   do id2cd = 1,ndcd
      d2 = d2min + (d2max - d2min) * float(id2cd-1) / float(ndcd-1)
      r2 = en2 * .5236 * d2 ** 3

      CALL initg2mode (r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam &
       ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
      CALL sumn (ank0,amk0,1,ithresh-1,ibins,sun10,sum10)
      CALL sumn (ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
      CALL sumn (ank0,amk0,ithresh1+1,ibins,ibins,sun30,sum30)
      CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,2))
      CALL sumn (ank,amk,1,ithresh-1,ibins,sun1,sum1)
      CALL sumn (ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
      CALL sumn (ank,amk,ithresh1+1,ibins,ibins,sun3,sum3)

      r1tabcd(id1cd,id2cd) = (sum10-sum1)*en1i
      c1tabcd(id1cd,id2cd) = (sun10-sun1)*en1i
      r2tabcd(id1cd,id2cd) = (sum20-sum2)*en1i
      c2tabcd(id1cd,id2cd) = (sun20-sun2)*en1i

!      write(*,'(a3,2X,i2,2X,i2,2X,f6.1,2X,f6.1,6(2X,e10.3))') &
!       ,'CD',id1cd,id2cd,d1*10000,d2*10000 &
!       ,r1tabcd(id1cd,id2cd),r2tabcd(id1cd,id2cd),(sum30-sum3)*en1i &
!       ,c1tabcd(id1cd,id2cd),c2tabcd(id1cd,id2cd),(sun30-sun3)*en1i
   enddo
enddo

!**************************************************************************
!**************************** RAIN-DRIZZLE ********************************
!**************************************************************************
! Start 3 dr loops for r(2) c(2,3)   [r2 approx = (-r3)]
d1 = 4.e-4
r1 = en1 * .5236 * d1 ** 3

do idccr = 1,ndccr
   d2 = d2min * 10. ** (d2ecr * float(idccr-1))
   r2 = en2 * .5236 * d2 ** 3

   do irrcr = 1,nrrcr
      r3 = r3min * 10. ** (r3ecr * float(irrcr-1))
      d3minx = max(d3min,(r3 / (.1 * .5236)) ** .333333)
      d3ecr = alog10(d3max / d3minx) / float(ndrcr-1)

      do idrcr = 1,ndrcr
       d3 = d3minx * 10. ** (d3ecr * float(idrcr-1))
       en3 = r3 / (.5236 * d3 ** 3)

       CALL initg2mode (r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam &
        ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
       CALL sumn (ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
       CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,4))
       CALL sumn (ank,amk,ithresh,ithresh1,ibins,sun2,sum2)

       r2tabcr(idccr,irrcr,idrcr) = alog10(max(1.0e-20,sum20-sum2))
       c2tabcr(idccr,irrcr,idrcr) = alog10(max(1.0e-20,sun20-sun2))

!       write(*,'(a3,2X,i2,2X,i2,2X,i2,2X,f10.2,2X,f10.2)') &
!        ,'DR',idccr,irrcr,idrcr,r2tabcr(idccr,irrcr,idrcr) &
!        ,c2tabcr(idccr,irrcr,idrcr)
      enddo
   enddo
enddo

!**************************************************************************
!**************************** DRIZZLE-DRIZZLE *****************************
!**************************************************************************
! Start 1 dd loop for r(2) and c(2,3) [r2=(-r3)]
d1 = 4.e-4
r1 = en1 * .5236 * d1 ** 3
r3 = .01e-6

do idcc = 1,ndcc
   d2 = d2min + (d2max - d2min) * float(idcc-1) / float(ndcc-1)
   r2 = en2 * .5236 * d2 ** 3

   CALL initg2mode (r1,r2,r3,en1,en2,en3,gnu(1),gnu(8),gnu(2),diam &
    ,amk0,ank0,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)
   CALL sumn (ank0,amk0,ithresh,ithresh1,ibins,sun20,sum20)
   CALL sumn (ank0,amk0,ithresh1+1,ibins,ibins,sun30,sum30)
   CALL sxy (x,amk0,ank0,amk,ank,akbarx(1,1,5))
   CALL sumn (ank,amk,ithresh,ithresh1,ibins,sun2,sum2)
   CALL sumn (ank,amk,ithresh1+1,ibins,ibins,sun3,sum3)

   r2tabdd(idcc) = max(0.,sum20-sum2)
   c2tabdd(idcc) = max(0.,sun20-sun2)
   c3tabdd(idcc) = max(0.,sun3-sun30)

!   write(*,'(a4,2X,f7.1,2X,i2,4(2X,e10.3))') ,'DD',d2*10000,idcc &
!    ,r2tabdd(idcc),sum30-sum3,c2tabdd(idcc),c3tabdd(idcc)
enddo

endif !IF CLOUD2 turned on

!*************************************************************************
!*********************** CLOUD1/CLOUD2 ICE SPECIES ***********************
!*************************************************************************
maxcld=1
if(jnmb(8).gt.0) maxcld=2

do cld=1,maxcld
do lcat=4,7

do idccr = 1,ndccr
   if(cld==1)then
     d1 = d1min * 10. ** (d1ecr * float(idccr-1))
     r1 = en1 * .5236 * d1 ** 3
   elseif(cld==2)then
     d2 = d2min * 10. ** (d2ecr * float(idccr-1))
     r2 = en2 * .5236 * d2 ** 3
   endif

   do irrcr = 1,nrrcr
    ri = rimin * 10. ** (rieci * float(irrcr-1))
    !Convert ri (g to kg) & switch diminx (m to cm) after
    !Multiply by 200 instead of 100 to bump up the min diameter
    diminx=max(dimin,((ri/1000./cfmas(lcat))**(1./pwmas(lcat)))*200.)
    dieci = alog10(dimax / diminx) / float(ndrcr-1)

     do idrcr = 1,ndrcr
      di = diminx * 10. ** (dieci * float(idrcr-1))
      !Here ri has units of (g/cm3)
      enice = ri / (1000. * cfmas(lcat) * (di/100.) ** pwmas(lcat))

      if(cld==1) then
       CALL initg1mode (r1,ri,en1,enice,gnu(1),gnu(lcat) &
        ,diami(1,lcat-3),amk0,ank0,ank1,amk1 &
        ,ank3,amk3,icutoff,ibins,lcat)
       CALL sxyice (amk1,ank1,amk3,ank3,akbarxci(1,1,lcat-3),sum1 &
        ,collectormass,sun1)
      endif
      if(cld==2) then
       CALL initg1mode (r2,ri,en2,enice,gnu(8),gnu(lcat) &
        ,diami(1,lcat-3),amk0,ank0,ank2,amk2 &
        ,ank3,amk3,ithresh1,ibins,lcat)
       CALL sxyice (amk2,ank2,amk3,ank3,akbarxci(1,1,lcat-3) &
        ,sum1,collectormass,sun1)
      endif

   if(cld==1)then
    r1tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,sum1*en1i))
    c1tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,sun1*en1i))
    r1rimer(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,collectormass*en1i))
!    write(*,'(a2,2X,a5,i2,3(2X,a6,i2),3(2X,f10.2))') &
!     'CI','lcat=',lcat,'idccr=',idccr,'irrcr=',irrcr,'idrcr=',idrcr &
!     ,r1tabci(idccr,irrcr,idrcr,lcat-3) &
!     ,c1tabci(idccr,irrcr,idrcr,lcat-3),r1rimer(idccr,irrcr,idrcr,lcat-3)
   endif

   if(cld==2)then
    r2tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,sum1))
    c2tabci(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,sun1))
    r2rimer(idccr,irrcr,idrcr,lcat-3) = alog10(max(1.0e-20,collectormass))
!    write(*,'(a3,2X,a5,i2,3(2X,a6,i2),3(2X,f10.2))') &
!     'DI','lcat=',lcat,'idccr=',idccr,'irrcr=',irrcr,'idrcr=',idrcr &
!     ,r2tabci(idccr,irrcr,idrcr,lcat-3) &
!     ,c2tabci(idccr,irrcr,idrcr,lcat-3),r2rimer(idccr,irrcr,idrcr,lcat-3)
   endif

   enddo
 enddo
enddo

enddo !for looping lcat=4,7
enddo !for looping maxcld=1,2

return
END SUBROUTINE make_autotab

!##############################################################################
Subroutine sxy (x,amkd,ankd,amk,ank,akbar)

implicit none

integer, parameter :: ibins=36
integer :: i,ik,k,l
real, dimension(ibins+1) :: x
real, dimension(ibins,ibins) :: akbar
real, dimension(ibins) :: xave,ankd,ank,amkd,amk,am2,am3,am4,psi,f
real :: ap,pi,dm,dn,sm1,sm2,sm3,sm4,sm5,sn1,sn2,sn3,sn4,dm4,dm2

data ap/1.062500000/
data pi/3.141592654/

do l=1,ibins
   if(ankd(l).gt.0)then
      xave(l)=amkd(l)/ankd(l)
   else
      xave(l)=0.
      amkd(l)=0. !If bin number is zero, make sure mass is zero.
   endif
enddo

!LOOP THRU ALL 36 BINS OF TRIANGULAR MATRIX
do k=1,ibins

! calculation of the 2nd, 3rd, and 4th moments of the mass distribution
! based on equation 8 in reference.

   am2(k)=ap*xave(k)*amkd(k)
   am3(k)=ap*ap*xave(k)*am2(k)
   am4(k)=ap*ap*ap*xave(k)*am3(k)

! these come out of the linear approximations used to integrate
! over partial bins.  they are defined:
!      psi(k) = nk(k+1)
!        f(k) = nk(k)
! where nk is the distribution func.  see equation 13 in reference.

   psi(k)=2./x(k)*(amkd(k)/x(k)-ankd(k))
   f(k)=2./x(k)*(2.*ankd(k)-amkd(k)/x(k))
   !print*,'psi1',k,psi(k),f(k),amkd(k),ankd(k)

! zeroing the tendencies on the moments.

   sm1=0.
   sm2=0.
   sm3=0.
   sm4=0.
   sm5=0.
   sn1=0.
   sn2=0.
   sn3=0.
   sn4=0.

! calculation of tendencies on moments

   !Sum tendencies from current k bin with all larger bins.
   !This is a LOSS of mass and number from bin k.
   do i=k,ibins
      dm=akbar(i,k)*(am2(k)*ankd(i)+amkd(k)*amkd(i))
      dn=akbar(i,k)*(ankd(k)*amkd(i)+amkd(k)*ankd(i))
      sm5=sm5+dm
      sn4=sn4+dn
   enddo

   !Get self-collection tendencies for same bin interactions k-1.
   !This is a gain of mass and number for bin k 
   if(k.gt.1)then
      sm3=akbar(k-1,k-1)*(am2(k-1)*ankd(k-1)+amkd(k-1)**2)
      sn2=akbar(k-1,k-1)*ankd(k-1)*amkd(k-1)
      dn=sn2
      dm=sm3
   endif

   !Sum contributions of mass and number to current bin K due
   !to interactions between bin k and all previous bins up to k-1.
   do i=1,k-1
      !Source of mass to add to current bin k.
      dm4=akbar(k,i)*(ankd(k)*am2(i)+amkd(k)*amkd(i))
      sm4=sm4+dm4
      !If the average starting mass in bin k is greater than bin-center
      !mass then compute amount of mass and number to REMOVE from the bin
      !due to interactions between bin k and all previous bins up to k-1.
      if(xave(k).ge.x(k))then
         dm2=akbar(k,i)*(4.*x(k)**2*psi(k)*amkd(i)  &
            +0.5*x(k)*(4.*psi(k)+f(k))*am2(i)  &
            -(psi(k)-f(k))*am3(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am4(i))
         sm2=sm2+dm2
         dn=akbar(k,i)*(2.*x(k)*psi(k)*amkd(i)  &
            +0.5*f(k)*am2(i)  &
            -0.5/(x(k))*(psi(k)-f(k))*am3(i))
         sn3=sn3+dn
      endif
   enddo

   !Sum contributions of mass and number to current bin k due
   !to interactions between bin k-1 and all previous bins up to k-2.
   do i=1,k-2
      ik=k-1
      !Compute this mass and number to ADD to bin k if the average starting
      !mass in bin k-1 is greater than bin-center mass of k-1. Otherwise, 
      !mass can just be added to k-1 without crossing over to the larger 
      !mass in k      
      if(xave(ik).ge.x(ik))then
         dm=akbar(ik,i)*(4.*x(ik)**2*psi(ik)*amkd(i)  &
            +x(ik)/2.*(4.*psi(ik)+f(ik))*am2(i)  &
            -(psi(ik)-f(ik))*am3(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am4(i))
         sm1=sm1+dm
         dn=akbar(ik,i)*(2.*x(ik)*psi(ik)*amkd(i)  &
            +0.5*f(ik)*am2(i)  &
            -0.5/(x(ik))*(psi(ik)-f(ik))*am3(i))
         sn1=sn1+dn
      endif
   enddo

   amk(k)=amkd(k)+sm1-sm2+sm3+sm4-sm5
   ank(k)=ankd(k)+sn1+sn2-sn3-sn4

enddo

return
END SUBROUTINE sxy

!##############################################################################
Subroutine sxyice (amkdw,ankdw,amkdi,ankdi,akbar,smw,smi,sn)

implicit none

integer, parameter :: ibins=36
integer :: i,ik,k,l
real, dimension(ibins,ibins) :: akbar
real, dimension(ibins) :: xavew,xavei,ankdw,ankdi,amkdw,amkdi &
,am2w,am2i
real :: ap,pi,dmw,dmi,dn,smw,smi,sn

data ap/1.062500000/
data pi/3.141592654/

do l=1,ibins
   !Water spectrum
   if(ankdw(l).gt.0)then
      xavew(l)=amkdw(l)/ankdw(l)
   else
      xavew(l)=0.
      amkdw(l)=0. !If bin number is zero, make sure mass is zero.
   endif
   !Ice spectrum
   if(ankdi(l).gt.0)then
      xavei(l)=amkdi(l)/ankdi(l)
   else
      xavei(l)=0.
      amkdi(l)=0. !If bin number is zero, make sure mass is zero.
   endif
enddo

! zeroing the tendencies on the moments.
   smw=0.
   smi=0.
   sn=0.

!LOOP THRU ALL 36 BINS OF TRIANGULAR MATRIX
do k=1,ibins !ice bins

! calculation of tendencies on moments

   !Sum tendencies from current k bin with all bins.
   do i=1,ibins !water bins
      dn=akbar(k,i)*(ankdi(k)*amkdw(i)+amkdi(k)*ankdw(i))

      !Assumption is that dn is the same for water and ice
      dmw=dn*xavew(i)
      dmi=dn*xavei(k)
      
      sn=sn+dn    !Number of drops participating in collection
      smw=smw+dmw !Water participating in collection
      smi=smi+dmi !Ice participating in collection
   enddo
enddo

return
END SUBROUTINE sxyice

!##############################################################################
Subroutine data (x,diam,akbar,ibins)

implicit none

integer :: l,i,j,ibins,n,kernel
real :: pi
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbar
real, dimension(36,36) :: aabar,ahbar

kernel=1

if(kernel==1) then  !Long's collection kernel
data (aabar( 1,n),n=1, 1) /-.47757E-01  /
data (aabar( 2,n),n=1, 2) /-.26460E+00,-.47965E-01 /
data (aabar( 3,n),n=1, 3) /-.82258E+00,-.26760E+00,-.20453E-01 /
data (aabar( 4,n),n=1, 4) /-.19050E+01,-.82072E+00,-.11992E+00, .78909E-01 /
data (aabar( 5,n),n=1, 5) /-.39171E+01,-.18915E+01,-.33270E+00, .41936E+00  &
   ,.34801E+00 /
data (aabar( 6,n),n=1, 6) /-.76415E+01,-.38808E+01,-.73737E+00, .14121E+01  &
   ,.18851E+01,.99793E+00 /
data (aabar( 7,n),n=1, 7) /-.14595E+02,-.75638E+01,-.14861E+01, .33598E+01  &
   ,.61219E+01, .54314E+01, .24751E+01 /
data (aabar( 8,n),n=1, 8) /-.27720E+02,-.14442E+02,-.28741E+01, .69895E+01  &
   ,.14394E+02, .17479E+02, .13500E+02, .57110E+01 /
data (aabar( 9,n),n=1, 9) /-.52737E+02,-.27428E+02,-.54729E+01, .13703E+02  &
   ,.29792E+02, .40971E+02, .43267E+02, .31185E+02, .12630E+02 /
data (aabar(10,n),n=1,10) /-.10083E+03,-.52188E+02,-.10391E+02, .26218E+02  &
   ,.58283E+02, .84686E+02, .10128E+03, .99726E+02, .69014E+02, .27176E+02 /
data (aabar(11,n),n=1,11) /-.19396E+03,-.99799E+02,-.19790E+02, .49801E+02  &
   ,.11143E+03, .16558E+03, .20922E+03, .23326E+03, .22039E+03, .14858E+03  &
   ,.57396E+02 /
data (aabar(12,n),n=1,12) /-.37536E+03,-.19200E+03,-.37896E+02, .94692E+02  &
   ,.21165E+03, .31650E+03, .40896E+03, .48169E+03, .51524E+03, .47402E+03  &
   ,.31389E+03, .11962E+03 /
data (aabar(13,n),n=1,13) /-.73047E+03,-.37164E+03,-.73015E+02, .18089E+03  &
   ,.40253E+03, .60115E+03, .78166E+03, .94143E+03, .10638E+04, .11078E+04  &
   ,.10008E+04, .65436E+03, .24691E+03 /
data (aabar(14,n),n=1,14) /-.14285E+04,-.72333E+03,-.14152E+03, .34764E+03  &
   ,.76925E+03, .11434E+04, .14846E+04, .17993E+04, .20789E+04, .22870E+04  &
   ,.23385E+04, .20854E+04, .13509E+04, .50600E+03 /
data (aabar(15,n),n=1,15) /-.41365E+04,-.20869E+04,-.40697E+03, .99310E+03  &
   ,.21878E+04, .32394E+04, .41995E+04, .51084E+04, .59888E+04, .68297E+04  &
   ,.75528E+04, .79583E+04, .76785E+04, .62489E+04, .76776E+03 /
data (aabar(16,n),n=1,16) / .63760E+04, .64739E+04, .65970E+04, .67516E+04  &
   ,.69451E+04, .71861E+04, .74835E+04, .78448E+04, .82709E+04, .87453E+04  &
   ,.92111E+04, .95276E+04, .94079E+04, .83797E+04, .26045E+04, .89777E+03 /
data (aabar(17,n),n=1,17) / .62974E+04, .63746E+04, .64717E+04, .65934E+04  &
   ,.67457E+04, .69355E+04, .71702E+04, .74571E+04, .78005E+04, .81957E+04  &
   ,.86163E+04, .89879E+04, .91399E+04, .87394E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(18,n),n=1,18) / .62353E+04, .62963E+04, .63729E+04, .64689E+04  &
   ,.65889E+04, .67383E+04, .69233E+04, .71502E+04, .74238E+04, .77446E+04  &
   ,.81009E+04, .84538E+04, .87067E+04, .86514E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(19,n),n=1,19) / .61862E+04, .62344E+04, .62949E+04, .63707E+04  &
   ,.64653E+04, .65831E+04, .67290E+04, .69080E+04, .71250E+04, .73819E+04  &
   ,.76742E+04, .79815E+04, .82491E+04, .83524E+04, .66125E+04, .59471E+04  &
   ,.46530E+04, .26045E+04, .89777E+03 /
data (aabar(20,n),n=1,20) / .61474E+04, .61855E+04, .62334E+04, .62932E+04  &
   ,.63679E+04, .64608E+04, .65759E+04, .67172E+04, .68887E+04, .70932E+04  &
   ,.73291E+04, .75856E+04, .78311E+04, .79911E+04, .68735E+04, .66125E+04  &
   ,.59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(21,n),n=1,21) / .61166E+04, .61468E+04, .61847E+04, .62320E+04  &
   ,.62910E+04, .63644E+04, .64552E+04, .65668E+04, .67023E+04, .68644E+04  &
   ,.70531E+04, .72625E+04, .74738E+04, .76415E+04, .69140E+04, .68735E+04  &
   ,.66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(22,n),n=1,22) / .60923E+04, .61162E+04, .61462E+04, .61836E+04  &
   ,.62303E+04, .62883E+04, .63600E+04, .64481E+04, .65553E+04, .66836E+04  &
   ,.68338E+04, .70027E+04, .71786E+04, .73330E+04, .68498E+04, .69140E+04  &
   ,.68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
data (aabar(23,n),n=1,23) / .60730E+04, .60919E+04, .61157E+04, .61453E+04  &
   ,.61823E+04, .62281E+04, .62848E+04, .63545E+04, .64392E+04, .65408E+04  &
   ,.66601E+04, .67953E+04, .69391E+04, .70729E+04, .67447E+04, .68498E+04  &
   ,.69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04  &
   ,.89777E+03 /
data (aabar(24,n),n=1,24) / .60577E+04, .60727E+04, .60915E+04, .61150E+04  &
   ,.61443E+04, .61806E+04, .62254E+04, .62805E+04, .63475E+04, .64279E+04  &
   ,.65225E+04, .66304E+04, .67467E+04, .68590E+04, .66311E+04, .67447E+04  &
   ,.68498E+04, .69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04  &
   ,.26045E+04, .89777E+03 /
data (aabar(25,n),n=1,25) / .77967E+04, .78122E+04, .78316E+04, .78560E+04  &
   ,.78863E+04, .79242E+04, .79713E+04, .80294E+04, .81008E+04, .81878E+04  &
   ,.82924E+04, .84158E+04, .85571E+04, .87104E+04, .86265E+04, .88325E+04  &
   ,.90719E+04, .93363E+04, .95996E+04, .98007E+04, .98157E+04, .94274E+04  &
   ,.83361E+04, .63023E+04, .57988E+03 /
data (aabar(26,n),n=1,26) / .69349E+04, .69458E+04, .69595E+04, .69766E+04  &
   ,.69979E+04, .70244E+04, .70573E+04, .70978E+04, .71473E+04, .72072E+04  &
   ,.72788E+04, .73623E+04, .74565E+04, .75566E+04, .74715E+04, .76064E+04  &
   ,.77647E+04, .79435E+04, .81311E+04, .82983E+04, .83827E+04, .82640E+04  &
   ,.77406E+04, .65488E+04, .15807E+04, .51662E+03 /
data (aabar(27,n),n=1,27) / .61704E+04, .61781E+04, .61877E+04, .61997E+04  &
   ,.62147E+04, .62333E+04, .62562E+04, .62843E+04, .63186E+04, .63598E+04  &
   ,.64086E+04, .64648E+04, .65271E+04, .65912E+04, .65100E+04, .65961E+04  &
   ,.66976E+04, .68135E+04, .69382E+04, .70574E+04, .71390E+04, .71194E+04  &
   ,.68816E+04, .62379E+04, .26526E+04, .14083E+04, .46025E+03 /
data (aabar(28,n),n=1,28) / .54916E+04, .54971E+04, .55038E+04, .55123E+04  &
   ,.55228E+04, .55357E+04, .55517E+04, .55712E+04, .55949E+04, .56232E+04  &
   ,.56562E+04, .56938E+04, .57345E+04, .57747E+04, .57001E+04, .57533E+04  &
   ,.58161E+04, .58880E+04, .59661E+04, .60426E+04, .61008E+04, .61062E+04  &
   ,.59940E+04, .56500E+04, .31742E+04, .23632E+04, .12546E+04, .41004E+03 /
data (aabar(29,n),n=1,29) / .48886E+04, .48924E+04, .48971E+04, .49031E+04  &
   ,.49104E+04, .49195E+04, .49306E+04, .49441E+04, .49604E+04, .49797E+04  &
   ,.50020E+04, .50269E+04, .50530E+04, .50774E+04, .50108E+04, .50422E+04  &
   ,.50792E+04, .51212E+04, .51667E+04, .52111E+04, .52447E+04, .52486E+04  &
   ,.51861E+04, .49913E+04, .32935E+04, .28279E+04, .21054E+04, .11177E+04  &
   ,.36530E+03 /
data (aabar(30,n),n=1,30) / .43524E+04, .43551E+04, .43585E+04, .43626E+04  &
   ,.43678E+04, .43741E+04, .43818E+04, .43912E+04, .44024E+04, .44155E+04  &
   ,.44304E+04, .44467E+04, .44631E+04, .44771E+04, .44188E+04, .44361E+04  &
   ,.44561E+04, .44786E+04, .45022E+04, .45241E+04, .45384E+04, .45339E+04  &
   ,.44893E+04, .43663E+04, .31847E+04, .29342E+04, .25193E+04, .18757E+04  &
   ,.99579E+03, .32545E+03 /
data (aabar(31,n),n=1,31) / .38756E+04, .38775E+04, .38799E+04, .38828E+04  &
   ,.38864E+04, .38908E+04, .38961E+04, .39026E+04, .39102E+04, .39191E+04  &
   ,.39290E+04, .39395E+04, .39494E+04, .39568E+04, .39066E+04, .39149E+04  &
   ,.39241E+04, .39340E+04, .39435E+04, .39507E+04, .39516E+04, .39392E+04  &
   ,.39006E+04, .38129E+04, .29707E+04, .28372E+04, .26141E+04, .22445E+04  &
   ,.16710E+04, .88715E+03, .28994E+03 /
data (aabar(32,n),n=1,32) / .30106E+04, .30118E+04, .30132E+04, .30149E+04  &
   ,.30171E+04, .30197E+04, .30229E+04, .30266E+04, .30309E+04, .30357E+04  &
   ,.30408E+04, .30456E+04, .30491E+04, .30494E+04, .30032E+04, .30013E+04  &
   ,.29981E+04, .29929E+04, .29844E+04, .29706E+04, .29480E+04, .29111E+04  &
   ,.28504E+04, .27503E+04, .20956E+04, .19717E+04, .17926E+04, .15245E+04  &
   ,.11243E+04, .55892E+03, .22135E+03, .00000E+00 /
data (aabar(33,n),n=1,33) / .23888E+04, .23895E+04, .23903E+04, .23914E+04  &
   ,.23927E+04, .23943E+04, .23962E+04, .23983E+04, .24007E+04, .24033E+04  &
   ,.24057E+04, .24075E+04, .24077E+04, .24049E+04, .23645E+04, .23582E+04  &
   ,.23497E+04, .23382E+04, .23225E+04, .23007E+04, .22699E+04, .22258E+04  &
   ,.21613E+04, .20655E+04, .15572E+04, .14495E+04, .13057E+04, .11057E+04  &
   ,.82055E+03, .41732E+03, .17636E+03, .00000E+00, .00000E+00 /
data (aabar(34,n),n=1,34) / .18955E+04, .18959E+04, .18964E+04, .18971E+04  &
   ,.18979E+04, .18988E+04, .18999E+04, .19011E+04, .19024E+04, .19036E+04  &
   ,.19045E+04, .19047E+04, .19033E+04, .18990E+04, .18647E+04, .18567E+04  &
   ,.18462E+04, .18326E+04, .18145E+04, .17905E+04, .17581E+04, .17138E+04  &
   ,.16525E+04, .15662E+04, .11695E+04, .10771E+04, .95995E+03, .80547E+03  &
   ,.59516E+03, .30433E+03, .13287E+03, .00000E+00, .00000E+00, .00000E+00 /
data (aabar(35,n),n=1,35) / .15041E+04, .15044E+04, .15047E+04, .15051E+04  &
   ,.15056E+04, .15061E+04, .15067E+04, .15074E+04, .15080E+04, .15084E+04  &
   ,.15085E+04, .15079E+04, .15058E+04, .15011E+04, .14725E+04, .14642E+04  &
   ,.14536E+04, .14399E+04, .14221E+04, .13988E+04, .13682E+04, .13274E+04  &
   ,.12725E+04, .11976E+04, .88682E+03, .80897E+03, .71340E+03, .59221E+03  &
   ,.43360E+03, .22074E+03, .97241E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00 /
data (aabar(36,n),n=1,36) / .11936E+04, .11938E+04, .11940E+04, .11942E+04  &
   ,.11945E+04, .11948E+04, .11951E+04, .11954E+04, .11957E+04, .11957E+04  &
   ,.11954E+04, .11943E+04, .11921E+04, .11876E+04, .11640E+04, .11563E+04  &
   ,.11464E+04, .11337E+04, .11174E+04, .10963E+04, .10689E+04, .10330E+04  &
   ,.98554E+03, .92214E+03, .67808E+03, .61344E+03, .53582E+03, .44015E+03  &
   ,.31885E+03, .16089E+03, .70536E+02, .00000E+00, .00000E+00, .00000E+00  &
   ,.00000E+00, .00000E+00 /
endif

if(kernel==2) then  !Hall's collection kernel
data (ahbar( 1,n),n=1, 1)/ &
0.19826E+00/
data (ahbar( 2,n),n=1, 2)/ &
0.11115E+01, 0.31472E+00/
data (ahbar( 3,n),n=1, 3)/ &
0.29798E+01, 0.17644E+01, 0.49959E+00/
data (ahbar( 4,n),n=1, 4)/ &
0.48803E+01, 0.47301E+01, 0.28009E+01, 0.79305E+00/
data (ahbar( 5,n),n=1, 5)/ &
0.64832E+01, 0.77470E+01, 0.75086E+01, 0.44461E+01, 0.12589E+01/
data (ahbar( 6,n),n=1, 6)/ &
0.75069E+01, 0.10291E+02, 0.12298E+02, 0.11919E+02, 0.70577E+01, 0.19984E+01/
data (ahbar( 7,n),n=1, 7)/ &
0.97323E+01, 0.11917E+02, 0.16337E+02, 0.19521E+02, 0.18921E+02, 0.11203E+02,&
0.31722E+01/
data (ahbar( 8,n),n=1, 8)/ &
0.11705E+02, 0.15449E+02, 0.18916E+02, 0.25933E+02, 0.30988E+02, 0.30034E+02,&
0.17784E+02, 0.50356E+01/
data (ahbar( 9,n),n=1, 9)/ &
0.28432E+01, 0.18581E+02, 0.24524E+02, 0.30028E+02, 0.41166E+02, 0.49190E+02,&
0.47677E+02, 0.28231E+02, 0.79935E+01/
data (ahbar(10,n),n=1,10)/ &
0.11264E+01, 0.54577E+01, 0.24506E+02, 0.33433E+02, 0.44838E+02, 0.64398E+02,&
0.80924E+02, 0.77938E+02, 0.41342E+02, 0.10151E+02/
data (ahbar(11,n),n=1,11)/ &
0.19148E+00, 0.30048E+01, 0.99008E+01, 0.32365E+02, 0.45873E+02, 0.67472E+02,&
0.10098E+03, 0.13218E+03, 0.12667E+03, 0.61080E+02, 0.12790E+02/
data (ahbar(12,n),n=1,12)/ &
0.23333E+00, 0.24125E+00, 0.63028E+01, 0.17275E+02, 0.43141E+02, 0.63748E+02,&
0.10244E+03, 0.15873E+03, 0.21451E+03, 0.20480E+03, 0.91229E+02, 0.16114E+02/
data (ahbar(13,n),n=1,13)/ &
0.96499E+00, 0.18995E+01, 0.31461E+01, 0.21422E+02, 0.49047E+02, 0.95732E+02,&
0.19602E+03, 0.44815E+03, 0.86833E+03, 0.11774E+04, 0.10705E+04, 0.54659E+03,&
0.20565E+03/
data (ahbar(14,n),n=1,14)/ &
0.82031E+00, 0.11610E+02, 0.26319E+02, 0.45949E+02, 0.14169E+03, 0.27697E+03,&
0.47057E+03, 0.76235E+03, 0.13258E+04, 0.22194E+04, 0.28021E+04, 0.24852E+04,&
0.13356E+04, 0.61752E+03/
data (ahbar(15,n),n=1,15)/ &
0.49547E+01, 0.63623E+01, 0.11954E+03, 0.27411E+03, 0.48091E+03, 0.12745E+04,&
0.23941E+04, 0.39885E+04, 0.53365E+04, 0.64055E+04, 0.75273E+04, 0.78726E+04,&
0.71201E+04, 0.55775E+04, 0.12227E+04/
data (ahbar(16,n),n=1,16)/ &
0.19925E+02, 0.25489E+02, 0.32725E+02, 0.70048E+03, 0.16119E+04, 0.28297E+04,&
0.37418E+04, 0.48244E+04, 0.57636E+04, 0.68214E+04, 0.77280E+04, 0.83723E+04,&
0.84671E+04, 0.74915E+04, 0.23147E+04, 0.20493E+04/
data (ahbar(17,n),n=1,17)/ &
0.29658E+03, 0.37826E+03, 0.48383E+03, 0.62105E+03, 0.12813E+04, 0.21736E+04,&
0.33625E+04, 0.42800E+04, 0.53873E+04, 0.64521E+04, 0.73410E+04, 0.79541E+04,&
0.83280E+04, 0.80576E+04, 0.42900E+04, 0.25057E+04, 0.29399E+04/
data (ahbar(18,n),n=1,18)/ &
0.72103E+03, 0.91733E+03, 0.11698E+04, 0.14961E+04, 0.19199E+04, 0.25718E+04,&
0.34340E+04, 0.45765E+04, 0.53658E+04, 0.62564E+04, 0.68926E+04, 0.76242E+04,&
0.81438E+04, 0.82870E+04, 0.57426E+04, 0.44930E+04, 0.25977E+04, 0.27291E+04/
data (ahbar(19,n),n=1,19)/ &
0.96659E+03, 0.12273E+04, 0.15613E+04, 0.19908E+04, 0.25456E+04, 0.32656E+04,&
0.38523E+04, 0.46057E+04, 0.55962E+04, 0.62746E+04, 0.69491E+04, 0.70313E+04,&
0.78366E+04, 0.83524E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04,&
0.89099E+03/
data (ahbar(20,n),n=1,20)/ &
0.97638E+03, 0.12378E+04, 0.15716E+04, 0.19991E+04, 0.25486E+04, 0.32579E+04,&
0.41778E+04, 0.46749E+04, 0.52921E+04, 0.60950E+04, 0.66107E+04, 0.71113E+04,&
0.71824E+04, 0.77993E+04, 0.68735E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04,&
0.26045E+04, 0.89099E+03/
data (ahbar(21,n),n=1,21)/ &
0.94823E+03, 0.12006E+04, 0.15220E+04, 0.19322E+04, 0.24575E+04, 0.31324E+04,&
0.40029E+04, 0.51305E+04, 0.55255E+04, 0.59887E+04, 0.65800E+04, 0.69248E+04,&
0.72481E+04, 0.73294E+04, 0.69140E+04, 0.68735E+04, 0.66125E+04, 0.59471E+04,&
0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(22,n),n=1,22)/ &
0.82817E+03, 0.10475E+04, 0.13263E+04, 0.16812E+04, 0.21341E+04, 0.27139E+04,&
0.34583E+04, 0.44175E+04, 0.56582E+04, 0.59651E+04, 0.62990E+04, 0.67126E+04,&
0.69632E+04, 0.72083E+04, 0.68455E+04, 0.69140E+04, 0.68735E+04, 0.66125E+04,&
0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(23,n),n=1,23)/ &
0.69439E+03, 0.87760E+03, 0.11100E+04, 0.14053E+04, 0.17812E+04, 0.22609E+04,&
0.28745E+04, 0.36617E+04, 0.46749E+04, 0.59830E+04, 0.62384E+04, 0.64949E+04,&
0.67994E+04, 0.69710E+04, 0.66897E+04, 0.68477E+04, 0.69140E+04, 0.68735E+04,&
0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(24,n),n=1,24)/ &
0.57838E+03, 0.73052E+03, 0.92325E+03, 0.11677E+04, 0.14782E+04, 0.18735E+04,&
0.23776E+04, 0.30221E+04, 0.38482E+04, 0.49098E+04, 0.62770E+04, 0.64812E+04,&
0.66606E+04, 0.68557E+04, 0.66311E+04, 0.67447E+04, 0.68498E+04, 0.69140E+04,&
0.68735E+04, 0.66125E+04, 0.59471E+04, 0.46530E+04, 0.26045E+04, 0.89099E+03/
data (ahbar(25,n),n=1,25)/ &
0.59084E+03, 0.74589E+03, 0.94211E+03, 0.11907E+04, 0.15060E+04, 0.19065E+04,&
0.24163E+04, 0.30665E+04, 0.38980E+04, 0.49639E+04, 0.63340E+04, 0.80991E+04,&
0.83645E+04, 0.85993E+04, 0.86224E+04, 0.88325E+04, 0.90719E+04, 0.93363E+04,&
0.95996E+04, 0.98007E+04, 0.98157E+04, 0.94274E+04, 0.83361E+04, 0.63023E+04,&
0.57550E+03/
data (ahbar(26,n),n=1,26)/ &
0.41712E+03, 0.52636E+03, 0.66448E+03, 0.83925E+03, 0.10606E+04, 0.13414E+04,&
0.16979E+04, 0.21515E+04, 0.27296E+04, 0.34680E+04, 0.44128E+04, 0.56235E+04,&
0.71759E+04, 0.73865E+04, 0.73763E+04, 0.76028E+04, 0.77647E+04, 0.79435E+04,&
0.81311E+04, 0.82983E+04, 0.83827E+04, 0.82640E+04, 0.77406E+04, 0.65488E+04,&
0.15807E+04, 0.51271E+03/
data (ahbar(27,n),n=1,27)/ &
0.29457E+03, 0.37160E+03, 0.46891E+03, 0.59194E+03, 0.74760E+03, 0.94473E+03,&
0.11947E+04, 0.15119E+04, 0.19153E+04, 0.24289E+04, 0.30837E+04, 0.39193E+04,&
0.49856E+04, 0.63431E+04, 0.63635E+04, 0.65119E+04, 0.66944E+04, 0.68135E+04,&
0.69382E+04, 0.70574E+04, 0.71390E+04, 0.71194E+04, 0.68816E+04, 0.62379E+04,&
0.26526E+04, 0.14083E+04, 0.45678E+03/
data (ahbar(28,n),n=1,28)/ &
0.20808E+03, 0.26243E+03, 0.33104E+03, 0.41773E+03, 0.52730E+03, 0.66592E+03,&
0.84143E+03, 0.10639E+04, 0.13461E+04, 0.17045E+04, 0.21602E+04, 0.27398E+04,&
0.34765E+04, 0.44108E+04, 0.54855E+04, 0.56238E+04, 0.57419E+04, 0.58852E+04,&
0.59661E+04, 0.60426E+04, 0.61008E+04, 0.61062E+04, 0.59940E+04, 0.56500E+04,&
0.31742E+04, 0.23632E+04, 0.12546E+04, 0.40694E+03/
data (ahbar(29,n),n=1,29)/ &
0.14702E+03, 0.18538E+03, 0.23379E+03, 0.29491E+03, 0.37212E+03, 0.46970E+03,&
0.59313E+03, 0.74934E+03, 0.94723E+03, 0.11981E+04, 0.15162E+04, 0.19199E+04,&
0.24314E+04, 0.30782E+04, 0.38273E+04, 0.48524E+04, 0.49649E+04, 0.50559E+04,&
0.51642E+04, 0.52111E+04, 0.52447E+04, 0.52486E+04, 0.51861E+04, 0.49913E+04,&
0.32935E+04, 0.28279E+04, 0.21054E+04, 0.11177E+04, 0.36254E+03/
data (ahbar(30,n),n=1,30)/ &
0.10389E+03, 0.13098E+03, 0.16515E+03, 0.20827E+03, 0.26271E+03, 0.33147E+03,&
0.41837E+03, 0.52824E+03, 0.66723E+03, 0.84316E+03, 0.10659E+04, 0.13479E+04,&
0.17045E+04, 0.21543E+04, 0.26789E+04, 0.33884E+04, 0.42884E+04, 0.43778E+04,&
0.44448E+04, 0.45220E+04, 0.45384E+04, 0.45339E+04, 0.44893E+04, 0.43663E+04,&
0.31847E+04, 0.29342E+04, 0.25193E+04, 0.18757E+04, 0.99579E+03, 0.32299E+03/
data (ahbar(31,n),n=1,31)/ &
0.73425E+02, 0.92555E+02, 0.11668E+03, 0.14712E+03, 0.18553E+03, 0.23402E+03,&
0.29525E+03, 0.37261E+03, 0.47038E+03, 0.59398E+03, 0.75026E+03, 0.94780E+03,&
0.11972E+04, 0.15111E+04, 0.18798E+04, 0.23734E+04, 0.29973E+04, 0.37859E+04,&
0.38548E+04, 0.39003E+04, 0.39497E+04, 0.39392E+04, 0.39006E+04, 0.38129E+04,&
0.29707E+04, 0.28372E+04, 0.26141E+04, 0.22445E+04, 0.16710E+04, 0.88715E+03,&
0.28775E+03/
data (ahbar(32,n),n=1,32)/ &
0.45270E+02, 0.57059E+02, 0.71923E+02, 0.90671E+02, 0.11432E+03, 0.14416E+03,&
0.18182E+03, 0.22936E+03, 0.28938E+03, 0.36518E+03, 0.46087E+03, 0.58157E+03,&
0.73358E+03, 0.92435E+03, 0.11470E+04, 0.14441E+04, 0.18176E+04, 0.22860E+04,&
0.28721E+04, 0.29038E+04, 0.29104E+04, 0.29097E+04, 0.28504E+04, 0.27503E+04,&
0.20956E+04, 0.19717E+04, 0.17926E+04, 0.15245E+04, 0.11243E+04, 0.55892E+03,&
0.21886E+03, 0.00000E+00/
data (ahbar(33,n),n=1,33)/ &
0.28509E+02, 0.35930E+02, 0.45286E+02, 0.57082E+02, 0.71958E+02, 0.90722E+02,&
0.11439E+03, 0.14425E+03, 0.18193E+03, 0.22946E+03, 0.28939E+03, 0.36488E+03,&
0.45978E+03, 0.57859E+03, 0.71675E+03, 0.90062E+03, 0.11306E+04, 0.14176E+04,&
0.17740E+04, 0.22141E+04, 0.22189E+04, 0.21974E+04, 0.21603E+04, 0.20655E+04,&
0.15572E+04, 0.14495E+04, 0.13057E+04, 0.11057E+04, 0.82055E+03, 0.41732E+03,&
0.17442E+03, 0.00000E+00, 0.00000E+00/
data (ahbar(34,n),n=1,34)/ &
0.17955E+02, 0.22627E+02, 0.28517E+02, 0.35941E+02, 0.45302E+02, 0.57105E+02,&
0.71988E+02, 0.90757E+02, 0.11442E+03, 0.14426E+03, 0.18184E+03, 0.22912E+03,&
0.28847E+03, 0.36263E+03, 0.44864E+03, 0.56281E+03, 0.70510E+03, 0.88179E+03,&
0.11001E+04, 0.13676E+04, 0.16919E+04, 0.16753E+04, 0.16315E+04, 0.15655E+04,&
0.11695E+04, 0.10771E+04, 0.95995E+03, 0.80547E+03, 0.59516E+03, 0.30433E+03,&
0.13143E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (ahbar(35,n),n=1,35)/ &
0.11309E+02, 0.14251E+02, 0.17959E+02, 0.22632E+02, 0.28524E+02, 0.35950E+02,&
0.45313E+02, 0.57115E+02, 0.71989E+02, 0.90728E+02, 0.11432E+03, 0.14397E+03,&
0.18114E+03, 0.22752E+03, 0.28119E+03, 0.35228E+03, 0.44063E+03, 0.54992E+03,&
0.68429E+03, 0.84805E+03, 0.10451E+04, 0.12774E+04, 0.12438E+04, 0.11823E+04,&
0.88640E+03, 0.80897E+03, 0.71340E+03, 0.59221E+03, 0.43360E+03, 0.22074E+03,&
0.96191E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (ahbar(36,n),n=1,36)/ &
0.71229E+01, 0.89755E+01, 0.11310E+02, 0.14253E+02, 0.17962E+02, 0.22636E+02,&
0.28527E+02, 0.35951E+02, 0.45304E+02, 0.57082E+02, 0.71899E+02, 0.90509E+02,&
0.11382E+03, 0.14286E+03, 0.17642E+03, 0.22080E+03, 0.27581E+03, 0.34365E+03,&
0.42675E+03, 0.52753E+03, 0.64803E+03, 0.78904E+03, 0.94844E+03, 0.90139E+03,&
0.66943E+03, 0.61315E+03, 0.53582E+03, 0.44015E+03, 0.31885E+03, 0.16089E+03,&
0.69774E+02, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/ 
endif

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,ibins+1
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo

!Long's or Halls's collection kernel as calculated
!by Walko or Saleeby (1-36) with (1-36) weighted for (x+y)
do i=1,ibins
   do j=1,i
      if(kernel==1) akbar(i,j)=aabar(i,j)
      if(kernel==2) akbar(i,j)=ahbar(i,j)
      if(akbar(i,j).lt.0.) akbar(i,j)=0.
      akbar(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbar(i,j)=akbar(j,i)
   enddo
enddo

return
END SUBROUTINE data

!##############################################################################
Subroutine data_cs (x,diam,akbarcs,ibins,cfmass,pwmass)

implicit none

integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarcs
real, dimension(36,36) :: abarcs

data (abarcs( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcs( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcs( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.53661E+03, 0.10639E+04, 0.17461E+04,&
0.24939E+04, 0.33922E+04, 0.42133E+04, 0.49714E+04, 0.55472E+04, 0.60605E+04,&
0.63566E+04, 0.61133E+04, 0.44610E+04, 0.00000E+00/
data (abarcs(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.78871E+03, 0.13744E+04, 0.21286E+04,&
0.28583E+04, 0.37032E+04, 0.43698E+04, 0.50716E+04, 0.56088E+04, 0.61575E+04,&
0.67201E+04, 0.69407E+04, 0.62550E+04, 0.00000E+00, 0.00000E+00/
data (abarcs(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.98141E+03, 0.15955E+04, 0.23830E+04,&
0.30290E+04, 0.37710E+04, 0.43724E+04, 0.49581E+04, 0.54050E+04, 0.59129E+04,&
0.64845E+04, 0.68763E+04, 0.68208E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.94193E+03, 0.15293E+04, 0.22799E+04,&
0.28902E+04, 0.35847E+04, 0.41352E+04, 0.46578E+04, 0.50365E+04, 0.54635E+04,&
0.59603E+04, 0.63647E+04, 0.65881E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.90457E+03, 0.14671E+04, 0.21839E+04,&
0.27625E+04, 0.34159E+04, 0.39236E+04, 0.43940E+04, 0.47161E+04, 0.50718E+04,&
0.54881E+04, 0.58466E+04, 0.61587E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.86909E+03, 0.14084E+04, 0.20940E+04,&
0.26442E+04, 0.32615E+04, 0.37332E+04, 0.41607E+04, 0.44372E+04, 0.47337E+04,&
0.50775E+04, 0.53735E+04, 0.56824E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.83531E+03, 0.13528E+04, 0.20093E+04,&
0.25339E+04, 0.31192E+04, 0.35602E+04, 0.39523E+04, 0.41924E+04, 0.44415E+04,&
0.47247E+04, 0.49598E+04, 0.52295E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.80308E+03, 0.12999E+04, 0.19294E+04,&
0.24303E+04, 0.29870E+04, 0.34017E+04, 0.37643E+04, 0.39754E+04, 0.41870E+04,&
0.44214E+04, 0.46040E+04, 0.48253E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.77228E+03, 0.12496E+04, 0.18535E+04,&
0.23328E+04, 0.28635E+04, 0.32551E+04, 0.35929E+04, 0.37809E+04, 0.39630E+04,&
0.41587E+04, 0.42986E+04, 0.44740E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.74279E+03, 0.12015E+04, 0.17813E+04,&
0.22404E+04, 0.27474E+04, 0.31187E+04, 0.34353E+04, 0.36048E+04, 0.37635E+04,&
0.39289E+04, 0.40351E+04, 0.41711E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.71453E+03, 0.11555E+04, 0.17125E+04,&
0.21526E+04, 0.26377E+04, 0.29909E+04, 0.32892E+04, 0.34436E+04, 0.35840E+04,&
0.37256E+04, 0.38055E+04, 0.39095E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.68743E+03, 0.11114E+04, 0.16467E+04,&
0.20691E+04, 0.25338E+04, 0.28705E+04, 0.31528E+04, 0.32948E+04, 0.34205E+04,&
0.35435E+04, 0.36034E+04, 0.36821E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.66142E+03, 0.10692E+04, 0.15838E+04,&
0.19894E+04, 0.24350E+04, 0.27566E+04, 0.30247E+04, 0.31565E+04, 0.32704E+04,&
0.33787E+04, 0.34234E+04, 0.34824E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.63643E+03, 0.10287E+04, 0.15235E+04,&
0.19131E+04, 0.23408E+04, 0.26485E+04, 0.29038E+04, 0.30270E+04, 0.31313E+04,&
0.32280E+04, 0.32611E+04, 0.33051E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.61242E+03, 0.98978E+03, 0.14657E+04,&
0.18402E+04, 0.22508E+04, 0.25457E+04, 0.27893E+04, 0.29051E+04, 0.30015E+04,&
0.30889E+04, 0.31134E+04, 0.31459E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.58935E+03, 0.95241E+03, 0.14102E+04,&
0.17702E+04, 0.21648E+04, 0.24475E+04, 0.26805E+04, 0.27898E+04, 0.28796E+04,&
0.29595E+04, 0.29775E+04, 0.30014E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcs(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.56716E+03, 0.91651E+03, 0.13569E+04,&
0.17031E+04, 0.20823E+04, 0.23537E+04, 0.25768E+04, 0.26805E+04, 0.27646E+04,&
0.28383E+04, 0.28516E+04, 0.28690E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcs(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.54583E+03, 0.88199E+03, 0.13057E+04,&
0.16387E+04, 0.20033E+04, 0.22639E+04, 0.24778E+04, 0.25764E+04, 0.26557E+04,&
0.27243E+04, 0.27339E+04, 0.27466E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.52531E+03, 0.84880E+03, 0.12565E+04,&
0.15768E+04, 0.19274E+04, 0.21778E+04, 0.23830E+04, 0.24771E+04, 0.25522E+04,&
0.26165E+04, 0.26234E+04, 0.26326E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.50557E+03, 0.81688E+03, 0.12092E+04,&
0.15174E+04, 0.18546E+04, 0.20953E+04, 0.22923E+04, 0.23822E+04, 0.24536E+04,&
0.25141E+04, 0.25191E+04, 0.25256E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcs(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.48657E+03, 0.78617E+03, 0.11637E+04,&
0.14602E+04, 0.17846E+04, 0.20160E+04, 0.22053E+04, 0.22913E+04, 0.23594E+04,&
0.24167E+04, 0.24203E+04, 0.24249E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Wang & Ji 2000) collection kernel by Saleeby 6/23/04 (1-36) with (1-36)
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarcs(i,j)=abarcs(i,j)
      if(akbarcs(i,j).lt.0.) akbarcs(i,j)=0.
      akbarcs(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarcs(i,j)=akbarcs(j,i)
   enddo
enddo

return
END SUBROUTINE data_cs

!##############################################################################
Subroutine data_ca (x,diam,akbarca,ibins,cfmass,pwmass)

implicit none

integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarca
real, dimension(36,36) :: abarca

data (abarca( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarca( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarca( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.28025E+03, 0.19753E+04, 0.44248E+04, 0.51978E+04, 0.25782E+04,&
0.17774E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15165E+02, 0.21553E+02, 0.30125E+02,&
0.95567E+02, 0.54009E+03, 0.27604E+04, 0.42918E+04, 0.52676E+04, 0.59298E+04,&
0.46823E+04, 0.58190E+03, 0.47702E+03, 0.00000E+00, 0.00000E+00/
data (abarca(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.80274E+02, 0.11345E+03, 0.15744E+03,&
0.49499E+03, 0.11319E+04, 0.27491E+04, 0.39822E+04, 0.47690E+04, 0.54275E+04,&
0.50150E+04, 0.28583E+04, 0.24930E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.14965E+03, 0.21054E+03, 0.29048E+03,&
0.90652E+03, 0.17427E+04, 0.28230E+04, 0.38189E+04, 0.44627E+04, 0.51297E+04,&
0.53979E+04, 0.48803E+04, 0.43975E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22645E+03, 0.34181E+03, 0.49202E+03,&
0.11955E+04, 0.20566E+04, 0.28031E+04, 0.35326E+04, 0.41222E+04, 0.47278E+04,&
0.52342E+04, 0.55333E+04, 0.51884E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.30376E+03, 0.49729E+03, 0.74755E+03,&
0.13160E+04, 0.20052E+04, 0.26270E+04, 0.30831E+04, 0.36595E+04, 0.41208E+04,&
0.44659E+04, 0.47916E+04, 0.47526E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.34992E+03, 0.58085E+03, 0.87785E+03,&
0.13542E+04, 0.19195E+04, 0.24082E+04, 0.27781E+04, 0.31993E+04, 0.35298E+04,&
0.37789E+04, 0.40148E+04, 0.40734E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.36413E+03, 0.59007E+03, 0.87939E+03,&
0.12668E+04, 0.17268E+04, 0.21466E+04, 0.25022E+04, 0.28012E+04, 0.30599E+04,&
0.32423E+04, 0.34123E+04, 0.34479E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.36375E+03, 0.56448E+03, 0.82063E+03,&
0.11622E+04, 0.15651E+04, 0.19284E+04, 0.22304E+04, 0.24484E+04, 0.26537E+04,&
0.27880E+04, 0.29102E+04, 0.29263E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.32438E+03, 0.50280E+03, 0.72984E+03,&
0.10315E+04, 0.13856E+04, 0.17017E+04, 0.19600E+04, 0.21406E+04, 0.23060E+04,&
0.24061E+04, 0.24944E+04, 0.24966E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.28950E+03, 0.44831E+03, 0.64996E+03,&
0.91719E+03, 0.12295E+04, 0.15060E+04, 0.17290E+04, 0.18807E+04, 0.20161E+04,&
0.20921E+04, 0.21568E+04, 0.21500E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.25852E+03, 0.40004E+03, 0.57943E+03,&
0.81665E+03, 0.10930E+04, 0.13361E+04, 0.15299E+04, 0.16588E+04, 0.17715E+04,&
0.18303E+04, 0.18785E+04, 0.18666E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.23096E+03, 0.35720E+03, 0.51698E+03,&
0.72793E+03, 0.97301E+03, 0.11875E+04, 0.13571E+04, 0.14678E+04, 0.15628E+04,&
0.16093E+04, 0.16461E+04, 0.16315E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20642E+03, 0.31910E+03, 0.46158E+03,&
0.64942E+03, 0.86721E+03, 0.10571E+04, 0.12061E+04, 0.13020E+04, 0.13831E+04,&
0.14206E+04, 0.14494E+04, 0.14341E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18455E+03, 0.28518E+03, 0.41232E+03,&
0.57977E+03, 0.77362E+03, 0.94208E+03, 0.10736E+04, 0.11572E+04, 0.12272E+04,&
0.12581E+04, 0.12813E+04, 0.12663E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16503E+03, 0.25495E+03, 0.36848E+03,&
0.51788E+03, 0.69063E+03, 0.84038E+03, 0.95683E+03, 0.10302E+04, 0.10911E+04,&
0.11170E+04, 0.11362E+04, 0.11222E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarca(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.14760E+03, 0.22798E+03, 0.32940E+03,&
0.46280E+03, 0.61689E+03, 0.75022E+03, 0.85358E+03, 0.91828E+03, 0.97164E+03,&
0.99373E+03, 0.10099E+04, 0.99737E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarca(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13204E+03, 0.20390E+03, 0.29455E+03,&
0.41371E+03, 0.55127E+03, 0.67013E+03, 0.76204E+03, 0.81931E+03, 0.86633E+03,&
0.88544E+03, 0.89946E+03, 0.88837E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11812E+03, 0.18239E+03, 0.26344E+03,&
0.36994E+03, 0.49280E+03, 0.59886E+03, 0.68073E+03, 0.73156E+03, 0.77319E+03,&
0.78991E+03, 0.80225E+03, 0.79263E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10569E+03, 0.16318E+03, 0.23565E+03,&
0.33086E+03, 0.44066E+03, 0.53536E+03, 0.60838E+03, 0.65360E+03, 0.69057E+03,&
0.70534E+03, 0.71635E+03, 0.70812E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarca(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.94570E+02, 0.14600E+03, 0.21082E+03,&
0.29596E+03, 0.39412E+03, 0.47873E+03, 0.54392E+03, 0.58422E+03, 0.61715E+03,&
0.63028E+03, 0.64020E+03, 0.63323E+03, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Cober & List) collection kernel by Saleeby 6/23/04 for graupel/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarca(i,j)=abarca(i,j)
      if(akbarca(i,j).lt.0.) akbarca(i,j)=0.
      akbarca(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarca(i,j)=akbarca(j,i)
   enddo
enddo

return
END SUBROUTINE data_ca

!##############################################################################
Subroutine data_cg (x,diam,akbarcg,ibins,cfmass,pwmass)

implicit none

integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarcg
real, dimension(36,36) :: abarcg

data (abarcg( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarcg( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarcg( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.20097E+04, 0.27318E+04, 0.37157E+04,&
0.46118E+04, 0.57679E+04, 0.71030E+04, 0.83257E+04, 0.90047E+04, 0.95292E+04,&
0.95098E+04, 0.82568E+04, 0.48994E+04, 0.00000E+00/
data (abarcg(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16987E+04, 0.23372E+04, 0.31963E+04,&
0.39602E+04, 0.49292E+04, 0.60116E+04, 0.71421E+04, 0.77369E+04, 0.82349E+04,&
0.84640E+04, 0.79619E+04, 0.59465E+04, 0.00000E+00, 0.00000E+00/
data (abarcg(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.14234E+04, 0.19958E+04, 0.27574E+04,&
0.34313E+04, 0.42664E+04, 0.51380E+04, 0.61732E+04, 0.67127E+04, 0.71282E+04,&
0.74157E+04, 0.72934E+04, 0.62066E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.12169E+04, 0.17370E+04, 0.24222E+04,&
0.29970E+04, 0.37072E+04, 0.44738E+04, 0.53043E+04, 0.58948E+04, 0.62214E+04,&
0.64880E+04, 0.65268E+04, 0.59820E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10146E+04, 0.14914E+04, 0.21143E+04,&
0.26354E+04, 0.32668E+04, 0.39109E+04, 0.46347E+04, 0.52403E+04, 0.54895E+04,&
0.57105E+04, 0.58021E+04, 0.55433E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.84352E+03, 0.12842E+04, 0.18560E+04,&
0.23127E+04, 0.28662E+04, 0.34522E+04, 0.40906E+04, 0.46506E+04, 0.48989E+04,&
0.50724E+04, 0.51697E+04, 0.50527E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.67913E+03, 0.10889E+04, 0.16176E+04,&
0.20412E+04, 0.25492E+04, 0.30732E+04, 0.36139E+04, 0.41457E+04, 0.44188E+04,&
0.45511E+04, 0.46361E+04, 0.45862E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.55392E+03, 0.93719E+03, 0.14293E+04,&
0.18106E+04, 0.22657E+04, 0.27380E+04, 0.32206E+04, 0.37126E+04, 0.40228E+04,&
0.41228E+04, 0.41911E+04, 0.41713E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.44246E+03, 0.80234E+03, 0.12627E+04,&
0.16159E+04, 0.20349E+04, 0.24655E+04, 0.28794E+04, 0.33378E+04, 0.36877E+04,&
0.37675E+04, 0.38197E+04, 0.38122E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.35167E+03, 0.69069E+03, 0.11230E+04,&
0.14523E+04, 0.18408E+04, 0.22369E+04, 0.26134E+04, 0.30271E+04, 0.34046E+04,&
0.34689E+04, 0.35074E+04, 0.35041E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.33314E+03, 0.65352E+03, 0.10611E+04,&
0.13698E+04, 0.17326E+04, 0.21001E+04, 0.24465E+04, 0.28245E+04, 0.31656E+04,&
0.32144E+04, 0.32421E+04, 0.32396E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.31588E+03, 0.61910E+03, 0.10041E+04,&
0.12944E+04, 0.16345E+04, 0.19773E+04, 0.22980E+04, 0.26460E+04, 0.29572E+04,&
0.29945E+04, 0.30139E+04, 0.30108E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.29975E+03, 0.58706E+03, 0.95122E+03,&
0.12249E+04, 0.15448E+04, 0.18657E+04, 0.21644E+04, 0.24869E+04, 0.27731E+04,&
0.28018E+04, 0.28150E+04, 0.28112E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.28462E+03, 0.55710E+03, 0.90203E+03,&
0.11605E+04, 0.14620E+04, 0.17636E+04, 0.20429E+04, 0.23435E+04, 0.26084E+04,&
0.26308E+04, 0.26395E+04, 0.26350E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.27038E+03, 0.52899E+03, 0.85603E+03,&
0.11006E+04, 0.13854E+04, 0.16695E+04, 0.19316E+04, 0.22129E+04, 0.24596E+04,&
0.24773E+04, 0.24827E+04, 0.24779E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.25696E+03, 0.50254E+03, 0.81287E+03,&
0.10445E+04, 0.13139E+04, 0.15821E+04, 0.18289E+04, 0.20931E+04, 0.23239E+04,&
0.23381E+04, 0.23412E+04, 0.23363E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarcg(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.24427E+03, 0.47760E+03, 0.77225E+03,&
0.99191E+03, 0.12471E+04, 0.15007E+04, 0.17336E+04, 0.19824E+04, 0.21990E+04,&
0.22106E+04, 0.22123E+04, 0.22074E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarcg(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.23227E+03, 0.45403E+03, 0.73393E+03,&
0.94238E+03, 0.11843E+04, 0.14245E+04, 0.16446E+04, 0.18794E+04, 0.20835E+04,&
0.20932E+04, 0.20938E+04, 0.20892E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.22090E+03, 0.43173E+03, 0.69773E+03,&
0.89565E+03, 0.11252E+04, 0.13529E+04, 0.15613E+04, 0.17834E+04, 0.19760E+04,&
0.19842E+04, 0.19842E+04, 0.19800E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.21012E+03, 0.41060E+03, 0.66347E+03,&
0.85150E+03, 0.10695E+04, 0.12855E+04, 0.14830E+04, 0.16933E+04, 0.18755E+04,&
0.18826E+04, 0.18822E+04, 0.18783E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarcg(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.19990E+03, 0.39057E+03, 0.63102E+03,&
0.80970E+03, 0.10168E+04, 0.12219E+04, 0.14092E+04, 0.16086E+04, 0.17811E+04,&
0.17875E+04, 0.17868E+04, 0.17834E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Cober & List) collection kernel by Saleeby 6/23/04 for graupel/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarcg(i,j)=abarcg(i,j)
      if(akbarcg(i,j).lt.0.) akbarcg(i,j)=0.
      akbarcg(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarcg(i,j)=akbarcg(j,i)
   enddo
enddo

return
END SUBROUTINE data_cg

!##############################################################################
Subroutine data_ch (x,diam,akbarch,ibins,cfmass,pwmass)

implicit none

integer :: l,i,j,ibins,n
real :: pi,cfmass,pwmass
real, dimension(ibins+1) :: x,diam
real, dimension(ibins,ibins) :: akbarch
real, dimension(36,36) :: abarch

data (abarch( 1,n),n=1, 1)/ &
0.00000E+00/
data (abarch( 2,n),n=1, 2)/ &
0.00000E+00, 0.00000E+00/
data (abarch( 3,n),n=1, 3)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 4,n),n=1, 4)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 5,n),n=1, 5)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 6,n),n=1, 6)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch( 7,n),n=1, 7)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch( 8,n),n=1, 8)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch( 9,n),n=1, 9)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(10,n),n=1,10)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(11,n),n=1,11)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(12,n),n=1,12)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(13,n),n=1,13)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(14,n),n=1,14)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(15,n),n=1,15)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(16,n),n=1,16)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18260E+05, 0.20193E+05, 0.22753E+05,&
0.25470E+05, 0.28811E+05, 0.32097E+05, 0.34658E+05, 0.37576E+05, 0.40716E+05,&
0.43471E+05, 0.44623E+05, 0.42478E+05, 0.00000E+00/
data (abarch(17,n),n=1,17)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.15640E+05, 0.17168E+05, 0.19175E+05,&
0.21303E+05, 0.23918E+05, 0.26556E+05, 0.28532E+05, 0.30665E+05, 0.33125E+05,&
0.35677E+05, 0.37715E+05, 0.38066E+05, 0.00000E+00, 0.00000E+00/
data (abarch(18,n),n=1,18)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.13413E+05, 0.14624E+05, 0.16204E+05,&
0.17895E+05, 0.19982E+05, 0.22178E+05, 0.23750E+05, 0.25262E+05, 0.27063E+05,&
0.29078E+05, 0.31045E+05, 0.32337E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(19,n),n=1,19)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11513E+05, 0.12479E+05, 0.13731E+05,&
0.15093E+05, 0.16787E+05, 0.18664E+05, 0.19978E+05, 0.21041E+05, 0.22316E+05,&
0.23794E+05, 0.25369E+05, 0.26737E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(20,n),n=1,20)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.99570E+04, 0.10768E+05, 0.11814E+05,&
0.12914E+05, 0.14272E+05, 0.15772E+05, 0.16933E+05, 0.17716E+05, 0.18604E+05,&
0.19649E+05, 0.20808E+05, 0.21934E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(21,n),n=1,21)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.85951E+04, 0.92827E+04, 0.10166E+05,&
0.11059E+05, 0.12156E+05, 0.13373E+05, 0.14464E+05, 0.15063E+05, 0.15676E+05,&
0.16400E+05, 0.17218E+05, 0.18052E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(22,n),n=1,22)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.74082E+04, 0.79966E+04, 0.87497E+04,&
0.94753E+04, 0.10365E+05, 0.11367E+05, 0.12415E+05, 0.12916E+05, 0.13337E+05,&
0.13834E+05, 0.14397E+05, 0.14983E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(23,n),n=1,23)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.64420E+04, 0.69498E+04, 0.75977E+04,&
0.81863E+04, 0.89094E+04, 0.97458E+04, 0.10607E+05, 0.11155E+05, 0.11443E+05,&
0.11782E+05, 0.12164E+05, 0.12561E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(24,n),n=1,24)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.55999E+04, 0.60408E+04, 0.66020E+04,&
0.70880E+04, 0.76832E+04, 0.83750E+04, 0.90909E+04, 0.96786E+04, 0.98885E+04,&
0.10118E+05, 0.10374E+05, 0.10637E+05, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(25,n),n=1,25)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.48735E+04, 0.52583E+04, 0.57470E+04,&
0.61660E+04, 0.66700E+04, 0.72217E+04, 0.78235E+04, 0.84162E+04, 0.85953E+04,&
0.87498E+04, 0.89198E+04, 0.90907E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(26,n),n=1,26)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.42480E+04, 0.45742E+04, 0.49877E+04,&
0.53505E+04, 0.57835E+04, 0.62430E+04, 0.67549E+04, 0.72897E+04, 0.75073E+04,&
0.76107E+04, 0.77222E+04, 0.78304E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(27,n),n=1,27)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.37128E+04, 0.39889E+04, 0.43384E+04,&
0.46543E+04, 0.50280E+04, 0.54099E+04, 0.58491E+04, 0.63067E+04, 0.65815E+04,&
0.66513E+04, 0.67234E+04, 0.67898E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(28,n),n=1,28)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.32507E+04, 0.34939E+04, 0.38016E+04,&
0.40743E+04, 0.43920E+04, 0.46976E+04, 0.50766E+04, 0.54796E+04, 0.57855E+04,&
0.58353E+04, 0.58809E+04, 0.59199E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(29,n),n=1,29)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.28356E+04, 0.30476E+04, 0.33154E+04,&
0.35379E+04, 0.38059E+04, 0.41076E+04, 0.44301E+04, 0.47688E+04, 0.50874E+04,&
0.51352E+04, 0.51633E+04, 0.51845E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(30,n),n=1,30)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.24670E+04, 0.26464E+04, 0.28729E+04,&
0.30809E+04, 0.33290E+04, 0.35938E+04, 0.38390E+04, 0.41414E+04, 0.44456E+04,&
0.45303E+04, 0.45470E+04, 0.45569E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(31,n),n=1,31)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.21568E+04, 0.23211E+04, 0.25284E+04,&
0.26978E+04, 0.29021E+04, 0.31352E+04, 0.33501E+04, 0.36060E+04, 0.38781E+04,&
0.40046E+04, 0.40138E+04, 0.40168E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00/
data (abarch(32,n),n=1,32)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.18847E+04, 0.20216E+04, 0.21943E+04,&
0.23437E+04, 0.25242E+04, 0.27298E+04, 0.29188E+04, 0.31405E+04, 0.33788E+04,&
0.35454E+04, 0.35499E+04, 0.35489E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00/
data (abarch(33,n),n=1,33)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.16557E+04, 0.17717E+04, 0.19180E+04,&
0.20500E+04, 0.22098E+04, 0.23916E+04, 0.25583E+04, 0.27365E+04, 0.29461E+04,&
0.31428E+04, 0.31444E+04, 0.31411E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(34,n),n=1,34)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.14489E+04, 0.15486E+04, 0.16743E+04,&
0.17911E+04, 0.19318E+04, 0.20881E+04, 0.22352E+04, 0.23881E+04, 0.25675E+04,&
0.27711E+04, 0.27885E+04, 0.27841E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(35,n),n=1,35)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.12611E+04, 0.13496E+04, 0.14613E+04,&
0.15650E+04, 0.16877E+04, 0.18146E+04, 0.19444E+04, 0.20938E+04, 0.22391E+04,&
0.24180E+04, 0.24751E+04, 0.24705E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
data (abarch(36,n),n=1,36)/ &
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.11018E+04, 0.11805E+04, 0.12797E+04,&
0.13718E+04, 0.14800E+04, 0.15887E+04, 0.17036E+04, 0.18296E+04, 0.19537E+04,&
0.21119E+04, 0.21986E+04, 0.21941E+04, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,&
0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/

data pi/3.141592654/

!calculating the mass categories (c.g.s) with the lowest diameter
!of 3.125 microns and mass doubling every bin
diam(1)=1.5625*2.e-04
x(1)=pi/6.*diam(1)**3.
do l=2,15
   x(l)=2.*x(l-1)
   diam(l)=(6./pi*x(l))**(0.333333333333)
enddo
do l=16,ibins+1
   x(l)=2.*x(l-1)
   !uses power law values in m.k.s. and converts back to c.g.s.
   diam(l)=100 * ((x(l) / 1000.0 / cfmass) ** (1./pwmass))
enddo

!(Greenan & List) collection kernel by Saleeby 12/07/04 for hail/cloud droplet
! weighted for (x+y)
do i=1,ibins
   do j=1,i
      akbarch(i,j)=abarch(i,j)
      if(akbarch(i,j).lt.0.) akbarch(i,j)=0.
      akbarch(36,i)=0.
   enddo
enddo

!Invert here for use in routine 'sxy'
do j=1,ibins
   do i=1,j
      akbarch(i,j)=akbarch(j,i)
   enddo
enddo

return
END SUBROUTINE data_ch

!##############################################################################
Subroutine initg2mode (r1,r2,r3,n1,n2,n3,gnu1,gnu2,gnu3,diam,amk,ank  &
       ,ank1,amk1,ank2,amk2,ank3,amk3,ithresh,ithresh1,ibins)

use micphys

implicit none

integer :: i,ibins,ithresh,ithresh1
real :: r1,r2,r3,n1,n2,n3,gnu1,gnu2,gnu3,dn1,dn2,dn3,trunc,fac1,fac2,pi  &
       ,dmean,ex1
real, dimension(ibins+1) :: diam
real, dimension(ibins) :: amk,ank,amk1,ank1,amk2,ank2,amk3,ank3

real, external :: gammln
real, external :: gammp

! * *
! * Initial triple gamma distribution: n(D) = n1(D) + n2(D) + n3(D)
! * *

data pi/3.141592654/
data ex1/0.333333333/

! * *
! * gamma spectrum
! * *

do i=1,ibins
  ank1(i)=0.
  amk1(i)=0.
  ank2(i)=0.
  amk2(i)=0.
  ank3(i)=0.
  amk3(i)=0.
enddo

!*******************************************
! CLOUD
!*******************************************
dmean = (6.*r1/(pi*n1))**ex1   !mass mean diam
dn1 = dmean * (exp(gammln(gnu1)-gammln(gnu1+3.))) ** ex1

do i=1,ithresh
 fac1 = gammp(gnu1,diam(i)/dn1)
 fac2 = gammp(gnu1,diam(i+1)/dn1)
 trunc=fac2-fac1
 ank1(i)=n1*trunc

 fac1 = gammp(gnu1+3.,(diam(i)/dn1))
 fac2 = gammp(gnu1+3.,(diam(i+1)/dn1))
 trunc=fac2-fac1
 amk1(i)=r1*trunc
enddo

!This is a single bin patch for cases of truncation of cloud1
ank1(ithresh-1) = ank1(ithresh-1) + ank1(ithresh)
amk1(ithresh-1) = amk1(ithresh-1) + amk1(ithresh)
ank1(ithresh) = 0.0
amk1(ithresh) = 0.0

!********************************************
! CLOUD 2
!********************************************
dmean = (6. * r2 / (pi * n2)) ** ex1   !mass mean diam
dn2 = dmean * (exp(gammln(gnu2)-gammln(gnu2+3.))) ** ex1

do i=ithresh,ithresh1+1
 fac1 = gammp(gnu2,diam(i)/dn2)
 fac2 = gammp(gnu2,diam(i+1)/dn2)
 trunc=fac2-fac1
 ank2(i)=n2*trunc

 fac1 = gammp(gnu2+3.,(diam(i)/dn2))
 fac2 = gammp(gnu2+3.,(diam(i+1)/dn2))
 trunc=fac2-fac1
 amk2(i)=r2*trunc
enddo

!This is a single bin patch for cases of truncation of cloud2
ank2(ithresh1) = ank2(ithresh1) + ank2(ithresh1+1)
amk2(ithresh1) = amk2(ithresh1) + amk2(ithresh1+1)
ank2(ithresh1+1) = 0.0
amk2(ithresh1+1) = 0.0

!********************************************
! RAIN
!********************************************
dmean = (6. * r3 / (pi * n3)) ** ex1  !mass mean diam
dn3 = dmean * (exp(gammln(gnu3)-gammln(gnu3+3.))) ** ex1

do i=ithresh1+1,ibins
 fac1 = gammp(gnu3,diam(i)/dn3)
 fac2 = gammp(gnu3,diam(i+1)/dn3)
 trunc=fac2-fac1
 ank3(i)=n3*trunc

 fac1 = gammp(gnu3+3.,(diam(i)/dn3))
 fac2 = gammp(gnu3+3.,(diam(i+1)/dn3))
 trunc=fac2-fac1
 amk3(i)=r3*trunc
enddo

!**********************************************
! SUM THE NUMBER AND MASS OF EACH DISTRIBUTION
!**********************************************
do i=1,ibins
   ank(i)=ank1(i)+ank2(i)+ank3(i)
   amk(i)=amk1(i)+amk2(i)+amk3(i)
enddo

return
END SUBROUTINE initg2mode

!##############################################################################
Subroutine sumn (ank,amk,imin,imax,ibins,sun,sum)

implicit none

integer :: imin,imax,ibins
real :: sun,sum
real, dimension(ibins) :: ank,amk

integer :: i

  sum=0.
  sun=0.

  do i=imin,imax
   sun=sun+ank(i)
   sum=sum+amk(i)
  enddo

return
END SUBROUTINE sumn

!##############################################################################
Subroutine initg1mode (r1,r3,n1,n3,gnu1,gnu3,diam,amk,ank  &
       ,ank1,amk1,ank3,amk3,ithresh,ibins,hyd)

use micphys

implicit none

integer :: i,ibins,ithresh,hyd
real :: r1,r3,n1,n3,gnu1,gnu3,dn1,dn3,trunc,fac1,fac2,pi,dmean,ex1
real, dimension(ibins+1) :: diam
real, dimension(ibins) :: amk,ank,amk1,ank1,amk3,ank3
real, external :: gammln
real, external :: gammp

data pi/3.141592654/
data ex1/0.333333333/

do i=1,ibins
  ank1(i)=0.
  amk1(i)=0.
  ank3(i)=0.
  amk3(i)=0.
enddo

!*******************************************
! CLOUD
!*******************************************
dmean = (6.*r1/(pi*n1))**ex1   !mass mean diam
dn1 = dmean * (exp(gammln(gnu1)-gammln(gnu1+3.))) ** ex1

do i=1,ithresh
 fac1 = gammp(gnu1,diam(i)/dn1)
 fac2 = gammp(gnu1,diam(i+1)/dn1)
 trunc=fac2-fac1
 ank1(i)=n1*trunc

 fac1 = gammp(gnu1+3.,(diam(i)/dn1))
 fac2 = gammp(gnu1+3.,(diam(i+1)/dn1))
 trunc=fac2-fac1
 amk1(i)=r1*trunc
enddo

!********************************************
! RAIN, SNOW, AGGREGATES, GRAUPEL, or HAIL
!********************************************
dmean = 100 * (r3 / n3 / 1000. / cfmas(hyd)) ** (1./pwmas(hyd))
dn3 = dmean * (exp(gammln(gnu3)-gammln(gnu3+pwmas(hyd)))) ** (1./pwmas(hyd))

do i=ithresh+1,ibins
 fac1 = gammp(gnu3,diam(i)/dn3)
 fac2 = gammp(gnu3,diam(i+1)/dn3)
 trunc=fac2-fac1
 ank3(i)=n3*trunc

 fac1 = gammp(gnu3+pwmas(hyd),(diam(i)/dn3))
 fac2 = gammp(gnu3+pwmas(hyd),(diam(i+1)/dn3))
 trunc=fac2-fac1
 amk3(i)=r3*trunc
enddo

!**********************************************
! SUM THE NUMBER AND MASS OF EACH DISTRIBUTION
!**********************************************
do i=1,ibins
   ank(i)=ank1(i)+ank3(i)
   amk(i)=amk1(i)+amk3(i)
enddo

return
END SUBROUTINE initg1mode

!##############################################################################
Subroutine tabhab ()

use micphys
use mem_grid, only:iprntstmt,print_msg

implicit none

integer, parameter :: nhab=0
integer :: it,is

if(iprntstmt>=1 .and. print_msg)then
 if (nhab .eq.  0) print*,'VARIABLE HABIT PREDICTION'
 if (nhab .eq.  3) print*,'ASSUMED HABIT IS COLUMNS'
 if (nhab .eq.  8) print*,'ASSUMED HABIT IS HEX PLATES'
 if (nhab .eq.  9) print*,'ASSUMED HABIT IS DENDRITES'
 if (nhab .eq. 10) print*,'ASSUMED HABIT IS NEEDLES'
 if (nhab .eq. 11) print*,'ASSUMED HABIT IS ROSETTES'
 !if (nhab .eq.  x) print*,'ASSUMED HABIT IS SPHERES'
endif

! nt is temp, ns = satur (liq)

do it = 1,31
   do is = 1,100
      if (nhab .eq. 0) then
         if (it .ge. 0 .and. it .le. 2) then
            if (is .le. 95) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it .gt. 2 .and. it .le. 4) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            endif
         else if(it .gt. 4 .and. it .le. 6) then
            if (is .lt. 85) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it .gt. 6 .and. it .le. 9) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         else if(it .gt. 9 .and. it .le. 22) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 8
               jhabtab(it,is,2) = 12
            else
               jhabtab(it,is,1) = 9
               jhabtab(it,is,2) = 13
            endif
         elseif(it .gt. 22 .and. it .le. 30) then
            if (is .lt. 80) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 10
               jhabtab(it,is,2) = 14
            endif
         elseif(it .gt. 30) then
            if (is .lt. 90) then
               jhabtab(it,is,1) = 3
               jhabtab(it,is,2) = 4
            else
               jhabtab(it,is,1) = 11
               jhabtab(it,is,2) = 15
            endif
         endif
      else
         jhabtab(it,is,1) = nhab
         jhabtab(it,is,2) = nhab + 4
         if (nhab .eq. 3) jhabtab(it,is,2) = 4
      endif
   enddo
enddo

return
END SUBROUTINE tabhab
