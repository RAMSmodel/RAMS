!##############################################################################
Subroutine swrad (nz,nb,ns,npsb,              &
           u,pl,tl,dz,                        &
           xp,alpha,beta,wght,prf,trf,ral,       &
           solar,ngas,                           &
           alb,amu0,                             &
           tp,omgp,asym,fu,fd,flxsu,flxsd,ulim)

! * New swrad (Jerry, March 8, 1996)

implicit none

integer :: nz,nb,ns

!
!     two-stream radiative transfer code.
!     This version will attempt the use of FESFT
!     as outlined in R&G (MWR 1991)
!      JH 5-95
!
integer, parameter :: mb=8,mg=3,mk=7
real, parameter :: top=1800.,tm=1800./296.,gma=0.002,cp=1004.

!     input

real    :: u(nz,mg),pl(nz),tl(nz),dz(nz)
real    :: ulim(mg,mb)
integer :: npsb(mg,mb),na(mg),ngas(mg)

!     input parameters

real    :: ral(mb)
real    :: xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
real    :: wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
real    :: solar(mb)

!     arrays used in fluxes

real    :: tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb)
real    :: alb,amu0,asym(nz,nb)
real    :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real    :: re(nz),vd(nz),td(nz),vu(nz)
real    :: fu(nz,6),fd(nz,6),flxsu(nz),flxsd(nz)

integer ib,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,ig3,ik3
real :: uu,fact,ufac,dfac,tu1,tu2,td1,td2,fbu,fbd,tu3,td3

!     output

!     remember, for this code it is assumed that u, the gaseous
!     absorber amounts, are specified in Pa
!----------------------------------------------------------------------------
!     zero out flxsu and flxsd

!     do iz = 1,nz
!      flxsu(iz) = 0.0
!      flxsd(iz) = 0.0
!     enddo

!     loop through each band

do ib=1,ns

!        calculate here the properties that are considered grey,
!        i.e. averaged values are used across the band...
!        rayleigh scatter...cloud is now done outside of this routine

!           get rayleigh scattering AND
!           zero out the local flux arrays

     do iz=1,nz
         tcr(iz) = ral(ib)*pl(iz)*dz(iz)/  &
                              tl(iz)
         fu(iz,2) = 0.
         fu(iz,3) = 0.
         fu(iz,4) = 0.
         fd(iz,2) = 0.
         fd(iz,3) = 0.
         fd(iz,4) = 0.
     enddo

!        determine if, and how many overlaps...also check
!        to see if the gas is used

   mxa = 0
   do ig=1,mg
      if (npsb(ig,ib).gt.0.and.ngas(ig).eq.1)then
         mxa = mxa+1
         na(mxa) = ig
      endif
   enddo

   if (mxa.eq.0) go to 111
   if (mxa.eq.1) then

!-------------------------------------------------------------------------
!           no overlapping gasses, single aborber

      ig = na(1)

      do ik=1,npsb(ig,ib)

          do iz=nz,2,-1
            uu=min(ulim(ig,ib),u(iz,ig))
            tg(iz) = xp(ig,ik,ib)*uu*  &
               (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
               (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
          enddo

!               now do rest of stuff in routine

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

!               add pseudo-band fluxes to the total flux

          do iz=1,nz
            flxsu(iz)=flxsu(iz)+wght(ig,ik,ib)*fu(iz,6)
            flxsd(iz)=flxsd(iz)+wght(ig,ik,ib)*fd(iz,6)
          enddo

      enddo
   else if (mxa.eq.2) then

!--------------------------------------------------------------------------
!           overlap of two gasses using the FESFT

      ig1 = na(1)
      ig2 = na(2)

!           do the gray fluxes first

          tg(1:nz) = 0.

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),asym(1,ib))

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo

       enddo

!            do the 2nd gas

       do ik2=1,npsb(ig2,ib)

          do iz=1,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
      enddo

!           add together the fluxes: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,1)
!             fbd=fd(iz,2)/dfac*fd(iz,3)
!    +         /dfac*fd(iz,1)
!             flxsu(iz)=flxsu(iz)+fbu
!             flxsd(iz)=flxsd(iz)+fbd
!           enddo

! this is the stuff that is fixed (old stuff above)

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
        td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
        fbu = tu1*tu2*fu(iz,1)
        fbd = td1*td2*fd(iz,1)

        flxsu(iz)=flxsu(iz)+fbu
        flxsd(iz)=flxsd(iz)+fbd
      enddo


      else if (mxa.eq.3) then

!--------------------------------------------------------------------------
!           overlap of three gasses

      ig1 = na(1)
      ig2 = na(2)
      ig3 = na(3)

!           do the gray fluxes first

          tg(1:nz) = 0.

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),asym(1,ib))

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
        enddo

!           do the 2nd gas

      do ik2=1,npsb(ig2,ib)

          do iz=2,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
        enddo

!          do the 3rd gas

      do ik3=1,npsb(ig3,ib)

          do iz=2,nz
            uu=min(ulim(ig3,ib),u(iz,ig3))
            tg(iz) = xp(ig3,ik3,ib)*uu*  &
               (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
               (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
          enddo

!               now do rest of stuff in routine

          CALL flxsw (nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib))

!               sum the pseudo-band fluxes to get total flux

          fact = wght(ig3,ik3,ib)
          do iz=1,nz
            fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
            fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
          enddo
      enddo

!           sum the fluxes together: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,4)/ufac*fu(iz,1)
!             fbd=fd(iz,2)/dfac*fd(iz,3)
!    +         /dfac*fd(iz,4)/dfac*fd(iz,1)
!             flxsu(iz)=flxsu(iz)+fbu
!             flxsd(iz)=flxsd(iz)+fbd
!           enddo

! this is the stuff that is fixed (old stuff above)

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
        td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
        td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
        td3 = max(0.0,min(1.1,fd(iz,4)/dfac))

        fbu = tu1*tu2*tu3*fu(iz,1)
        fbd = td1*td2*td3*fd(iz,1)
        flxsu(iz)=flxsu(iz)+fbu
        flxsd(iz)=flxsd(iz)+fbd
      enddo


!-------------------------------------------------------------------------

   else
      write(*,5000)
5000       format('Two is the maximum amount of overlapping',/  &
            'gasses allowed: if you really want to have ',/  &
            'more gasses overlap, they better have few',/  &
            'pseudo-bands for all of them, and in any',/  &
            'case it will cost you. To do that come into',/  &
            'the code and modify it here, look at previous',/  &
            'structures to see how it is done. Its easy',/  &
            'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE')
      stop
   endif
111     continue

enddo

return
END SUBROUTINE swrad

!##############################################################################
Subroutine lwrad (i,j,nz,nb,ns,npsb,nuum,          &
           u,pl,tl,vp,                             &
           xp,alpha,beta,wght,prf,trf,ral,         &
           a0,a1,a2,a3,                            &
           exptabc,ngas,                           &
           tp,omgp,asym,fu,fd,flxlu,flxld,ulim)

! New lwrad (March 8, 1996)

implicit none

integer, parameter :: mb=8,mg=3,mk=7
real, parameter :: top=1800.,tm=1800./296.,gma=0.002,cp=1004.
integer :: nz,nb,ns,i,j
!
!     two-stream radiative transfer code.
!     This version will attempt the use of FESFT
!     as outlined in R&G (MWR 1991)
!      JH 5-95
!
!     input

real    u(nz,mg),pl(nz),tl(nz),vp(nz)
real    exptabc(150),ulim(mg,mb)
integer npsb(mg,mb),na(mg),nuum(mb),ngas(mg)

!     input parameters

real    ral(mb)
real    xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
real    wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
real    a0(mb),a1(mb),a2(mb),a3(mb)

!     arrays used in fluxes

real    tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb),src(nz)
real    t(nz),r(nz),sigu(nz),sigd(nz)
real    re(nz),vd(nz),td(nz),vu(nz),asym(nz,nb)
real    fu(nz,6),fd(nz,6),flxlu(nz),flxld(nz)

integer :: ib,iflag,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,nir,ii,ig3,ik3
real :: tf,ewght,expp,uu,chck,fact,ufac,tu1,tu2,fbu,dx,adx,fdx &
       ,dxc,tn1c,tn2c,tn1,tn2,difflxb,fbd,dfac,tu3,tn3c,tn3

!     remember, for this code it is assumed that u, the gaseous
!     absorber amounts, are specified in g/m^3
!----------------------------------------------------------------------------
!     zero out flxsu and flxsd

!     do iz = 1,nz
!      flxlu(iz) = 0.0
!      flxld(iz) = 0.0
!     enddo

!     loop through each band

nir=ns+1
do ib=nir,nb

!        calculate the properties that are grey across the band
!        Planck func at the interfaces, continuum absorption,
!        ...cloud is now done outside of this routine

     iflag=1 !Saleeby(2009): Set this to 1 to turn on water vapor continuum

!           do surface and top of atmosphere first

     src(1) = a0(ib)+tl(1)*(a1(ib)+tl(1)*  &
                     (a2(ib)+a3(ib)*tl(1)))
     src(2) = a0(ib)+tl(2)*(a1(ib)+tl(2)*  &
                     (a2(ib)+a3(ib)*tl(2)))
     src(nz) =  a0(ib)+tl(nz)*(a1(ib)+tl(nz)*  &
                         (a2(ib)+a3(ib)*tl(nz)))

     tcr(1)=0.
     tcr(2)=0.
     tcr(nz)=0.

     do iz=nz-1,3,-1

!                 get sources at the interface, temp first
            tf = 0.5*(tl(iz)+tl(iz-1))
            src(iz) = a0(ib)+tf*(a1(ib)+tf*  &
                         (a2(ib)+a3(ib)*tf))
            tcr(iz) = 0.

      enddo


      if (nuum(ib).eq.1) then

!              this band has continuum

         do iz=2,nz
             ii=int(tl(iz)-179.)
             ewght=tl(iz)-(float(ii)+179.)
             if(ii < 1) then !Deal with temps below 180K
               ii=1
               ewght=0.
             endif
             expp=(1.-ewght)*exptabc(ii)+ewght*exptabc(ii+1)
             tcr(iz) = ral(ib)*u(iz,1)*  &
                    expp*(vp(iz)+  &
                    gma*(pl(iz)-vp(iz)))
         enddo

      endif


   do iz=1,nz
         fu(iz,2) = 0.
         fu(iz,3) = 0.
         fu(iz,4) = 0.
         fd(iz,2) = 0.
         fd(iz,3) = 0.
         fd(iz,4) = 0.
         tg(iz)=0.
   enddo

!        determine if, and how many overlaps

   mxa = 0
   do ig=1,mg
      if (npsb(ig,ib).gt.0.and.ngas(ig).eq.1)then
         mxa = mxa+1
         na(mxa) = ig
      endif
   enddo

   if (mxa.eq.0) go to 111
   if (mxa.eq.1) then

!-------------------------------------------------------------------------
!           no overlapping gasses, single aborber

      ig = na(1)

      do ik=1,npsb(ig,ib)

          do iz=nz,2,-1
            uu=min(ulim(ig,ib),u(iz,ig))
            tg(iz) = xp(ig,ik,ib)*uu*  &
               (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
               (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
          enddo

!               now do rest of stuff in routine

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),1.,asym(1,ib))

!               add pseudo-band fluxes to the total flux

          do iz=1,nz
            flxlu(iz)=flxlu(iz)+wght(ig,ik,ib)*fu(iz,6)
            flxld(iz)=flxld(iz)+wght(ig,ik,ib)*fd(iz,6)
          enddo
      enddo
   else if (mxa.eq.2) then

!--------------------------------------------------------------------------
!           overlap of two gasses using the FESFT

      ig1 = na(1)
      ig2 = na(2)

!           do the gray fluxes first

      chck=float(iflag)

          tg(1:nz) = 0.

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),chck,asym(1,ib))

!           if there is continuum abs. do it now since tg=0
!
     if(nuum(ib).eq.1)then

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,5),fd(1,5),1.,asym(1,ib))
     endif

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
       enddo

!            do the 2nd gas

       do ik2=1,npsb(ig2,ib)

          do iz=1,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
      enddo

!           add together the fluxes: FESFT method...
!          NOTE: Here we use FESFT on the net flux NOT
!                the individual fluxes, as noted in R and G.
!                 FESFT seems not only to work for the net flux
!               but also for the upwelling IR. Thus, we may still use
!              FESFT to determine each flux individually.

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             dx=dfac-ufac
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,1)*(fu(iz,5)/ufac
!    +         *float(nuum(ib))+1.-float(nuum(ib)))
!             dx1=(fd(iz,2)-fu(iz,2))/dx
!             dx2=(fd(iz,3)-fu(iz,3))/dx
!             dxc=(fd(iz,5)-fu(iz,5))/dx*float(nuum(ib))
!    +         +1.-float(nuum(ib))
!             difflxb=dx1*dx2*dxc*dx
!             fbd=difflxb+fbu
!             flxlu(iz)=flxlu(iz)+fbu
!             flxld(iz)=flxld(iz)+fbd
!           enddo

!  above is old section, below is the new bounded part:

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        fbu = tu1*tu2*fu(iz,1)

        dx=fd(iz,1)-fu(iz,1)
        !Saleeby(08-10-2008) - added check to make sure upward flux does not
        !equal downward flux
        if (dx .eq. 0.) then
           dx=0.001
        endif

        adx = abs(dx)
        fdx = dx/adx
        dxc = max(1.e-14,adx)*fdx
        tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
        tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
        tn1 = min(max(tn1c,0.),2.) - 1.
        tn2 = min(max(tn2c,0.),2.) - 1.
        difflxb = tn1*tn2*dx
        fbd=difflxb+fbu

        flxlu(iz)=flxlu(iz)+fbu
        flxld(iz)=flxld(iz)+fbd
      enddo


      else if (mxa.eq.3) then

!--------------------------------------------------------------------------
!           overlap of three gasses

      ig1 = na(1)
      ig2 = na(2)
      ig3 = na(3)

!           do the gray fluxes first

      chck=float(iflag)

          tg(1:nz) = 0.

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),chck,asym(1,ib))

!           if there is continuum abs. do it now since tg=0
!
     if(nuum(ib).eq.1)then

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,5),fd(1,5),1.,asym(1,ib))
     endif

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
        enddo

!           do the 2nd gas

      do ik2=1,npsb(ig2,ib)

          do iz=2,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
        enddo

!          do the 3rd gas

      do ik3=1,npsb(ig3,ib)

          do iz=2,nz
            uu=min(ulim(ig3,ib),u(iz,ig3))
            tg(iz) = xp(ig3,ik3,ib)*uu*  &
               (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
               (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
          enddo

!               now do rest of stuff in routine

          CALL flxlw (i,j,ib,nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

!               sum the pseudo-band fluxes to get total flux

          fact = wght(ig3,ik3,ib)
          do iz=1,nz
            fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
            fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
          enddo
      enddo

!           sum the fluxes together: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             dx=dfac-ufac
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,4)/ufac*fu(iz,1)*(fu(iz,5)
!    +         /ufac*float(nuum(ib))+1.-float(nuum(ib)))
!             dx1=(fd(iz,2)-fu(iz,2))/dx
!             dx2=(fd(iz,3)-fu(iz,3))/dx
!             dx3=(fd(iz,4)-fu(iz,4))/dx
!             dxc=(fd(iz,5)-fu(iz,5))/dx*float(nuum(ib))
!    +         +1.-float(nuum(ib))
!             difflxb=dx1*dx2*dx3*dxc*dx
!             fbd=difflxb+fbu
!             flxlu(iz)=flxlu(iz)+fbu
!             flxld(iz)=flxld(iz)+fbd
!           enddo

! above is the old section, below is the new bounded section:


      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
        fbu = tu1*tu2*tu3*fu(iz,1)

        dx=fd(iz,1)-fu(iz,1)
        !Saleeby(08-10-2008) - added check to make sure upward flux does not
        !equal downward flux
        if (dx .eq. 0.) then
           dx=0.001
        endif

        adx = abs(dx)
        fdx = dx/adx
        dxc = max(1.e-14,adx)*fdx
        tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
        tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
        tn3c = (fd(iz,4)-fu(iz,4))/dxc + 1.
        tn1 = min(max(tn1c,0.),2.) - 1.
        tn2 = min(max(tn2c,0.),2.) - 1.
        tn3 = min(max(tn3c,0.),2.) - 1.
        difflxb = tn1*tn2*tn3*dx
        fbd=difflxb+fbu

        flxlu(iz)=flxlu(iz)+fbu
        flxld(iz)=flxld(iz)+fbd
      enddo



!-------------------------------------------------------------------------

   else
      write(*,5000)
5000       format('Two is the maximum amount of overlapping',/  &
            'gasses allowed: if you really want to have ',/  &
            'more gasses overlap, they better have few',/  &
            'pseudo-bands for all of them, and in any',/  &
            'case it will cost you. To do that come into',/  &
            'the code and modify it here, look at previous',/  &
            'structures to see how it is done. Its easy',/  &
            'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE')
      stop
   endif
111     continue

enddo

return
END SUBROUTINE lwrad

!##############################################################################
Subroutine flxsw (nz,tg,tp,tcr,omgp,  &
              alb,slr,amu0,t,r,tc,  &
              sigu,sigd,re,vd,td,vu,  &
              fu,fd,asym)

implicit none

integer :: nz
real :: tg(nz),tp(nz),tcr(nz),omgp(nz),  &
     alb,slr,amu0,asym(nz)

!     local variables

real :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real :: re(nz),vd(nz),td(nz),vu(nz),expt1

!     output variables

real :: fu(nz),fd(nz)

!     at this stage we have gotten rid of all the stuff
!     that are specific absorber gas dependent. The effects
!     of all the absorbing gasses in this specific pseudo-
!     band overlap (or no-overlap) is accumulated in tg.

!     parameter(asym=0.85,diffac=1.66,eps=1.e-6,etmp=0.367879)
real, parameter :: diffac=1.66,eps=1.e-6,etmp=0.367879
real, parameter :: e2=0.5,e3=0.3333,e4=0.25,e5=0.2,e6=0.166667
real, parameter :: e7=0.142857

integer :: iz
real :: pi,tau,omg0,af,fact,beta0,g1,g2,gg,rinf,ggtau,expp1,expp2,denom  &
       ,cc,g3,g4,aa,bb,tcm2,exp1,exp2

!     get total optical depth, single scattering albedo
!     and assymetry parameter
pi = 3.14159

   re(nz) = 0.
   vd(nz) = 0.
   expt1 = 1.
   fd(nz) = amu0*slr
   tc(nz) = 0.

do iz=nz,2,-1

      tau=tg(iz)+tp(iz)+tcr(iz)
      omg0 = min(.999999,(tcr(iz)+  &
           omgp(iz)*tp(iz))/max(1.e-20,tau))
      af   = asym(iz)*omgp(iz)*tp(iz)/  &
           (omg0*max(1.e-20,tau))

!           do delta-m scaling (wiscombe)

      fact        = af*af
      tau  = (1.-omg0*fact)*tau
      omg0 = ((1.-fact)*omg0)/(1.-omg0*fact)

!           determine the ODE matrix coefficients (Ritter and Geleyn)

      beta0     = (4.+af)/(8.*(1.+af))
      g1 = diffac*(1.-omg0*(1.-beta0))
      g2 = diffac*omg0*beta0
      gg = sqrt(g1**2-g2**2)

!           determine the local (true) reflection and transmission coefficients

      rinf  = g2/(gg+g1)
      ggtau=tau*gg

!Saleeby(08-10-2008): Old RAMS computation
      !expp1=1.-ggtau*(1.-ggtau*e2*(1.-ggtau*e3*(1.-ggtau*e4  &
      !       *(1.-ggtau*e5))))
!Saleeby(08-10-2008): From OLAM
      expp1=exp(-ggtau)

      expp1=max(expp1,0.0)
      expp2=expp1**2
      denom=(1.-rinf**2*expp2)
      t(iz) = ((1.-rinf**2)*expp1)/denom
      r(iz) = rinf*(1.-expp2)/denom

!        get the source func, go from top down to accomodate solar terms

!Saleeby(08-11-2008) - Original checked for dividing by 0.
!Variable 'fact' may still need bounds checking to prevent division by small number
         !if (gg-1./amu0.lt.eps) then
         !   fact = 1./(gg**2-1./(amu0+eps)**2)
         !else
         !  fact = 1./(gg**2-1./amu0**2)
         !endif
!Saleeby(08-11-2008) Version of 'fact' from OLAM - checked for dividing by 0
         if (abs(amu0*gg-1.0) < eps) then
            fact = amu0**2 / (gg**2*amu0**2 - 1.0 + sign(eps,amu0*gg-1.0))
         else
            fact = amu0**2 / (gg**2*amu0**2 - 1.0)
         endif

         cc          = omg0*slr*fact
         g3          = 0.5-0.75*af*amu0/(1.+af)
         g4          = 1.-g3
         aa          = g3*(g1-1./amu0)+g4*g2
         bb          = g4*(g1+1./amu0)+g3*g2
         tc(iz-1) = tc(iz)+tau
         tcm2=tc(iz-1)/amu0

!Saleeby(08-10-2008): Old RAMS computation
!         if(tcm2.le.1.0.or.tcm2.gt.10.0)then
!          exp2=1.-tcm2*(1.-tcm2*e2*(1.-tcm2*e3*(1.-tcm2*e4  &
!             *(1.-tcm2*e5*(1.-tcm2*e6*(1.-tcm2*e7))))))
!          exp2=max(exp2,0.0)
!         else
!           exp2=exp(-tcm2)
!          endif
!Saleeby(08-10-2008): From OLAM
         exp2=exp(-tcm2)

         exp1=expt1
         expt1=exp2
         sigu(iz) = cc*((aa-r(iz)*bb)*exp1-  &
           aa*t(iz)*exp2)
         sigd(iz) = cc*(-bb*t(iz)*  &
           exp1+(bb-r(iz)*aa)*exp2)
         fd(iz-1) = amu0*slr*exp2
      td(iz)   = 1. - re(iz)*r(iz)
      re(iz-1) = r(iz) + t(iz)**2*re(iz) /  &
                                      td(iz)
      vd(iz-1) = sigd(iz) + ( t(iz)*vd(iz) +  &
                    t(iz)*re(iz)*sigu(iz) ) /  &
                          td(iz)
      vu(iz-1) = ( r(iz)*vd(iz) + sigu(iz) ) /  &
                                      td(iz)

enddo

!     specify the boundary conditions

 fu(1)  = alb*(vd(1)+slr*amu0*  &
           exp(-tc(1)/amu0))/(1.-alb*re(1))

!     do adding, going from top down


!     calculate fluxes going up through the layers

do iz=2,nz
      fd(iz-1) = re(iz-1)*fu(iz-1) + vd(iz-1)  &
           + fd(iz-1)
      fu(iz)   = t(iz)*fu(iz-1)/td(iz)+vu(iz-1)
enddo

return
END SUBROUTINE flxsw

!##############################################################################
Subroutine flxlw (i,j,ib,nz,tg,tp,tcr,omgp,  &
                  src,t,r,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu,fd,chck,asym)

use node_mod, only:mi0,mj0
use mem_grid, only:ngrid

implicit none

integer :: nz,k,i,j,ib
real :: tg(nz),tp(nz),tcr(nz),omgp(nz),src(nz),chck,asym(nz)

!     local variables

real :: t(nz),r(nz),sigu(nz),sigd(nz)
real :: re(nz),vd(nz),td(nz),vu(nz)

!     output variables

real :: fu(nz),fd(nz)

!     parameter(asym=0.85,diffac=1.66,eps=1.e-6)
real, parameter :: diffac=1.66,eps=1.e-6,etmp=0.367879
real, parameter :: e2=0.5,e3=0.3333,e4=0.25,e5=0.2,e6=0.166667
real, parameter :: e7=0.142857

integer :: iz
real :: pi,tau,omg0,af,fact,beta0,g1,g2,gg,rinf,ggtau,expp1,expp2,aa,bb,cc

!     get total optical depth, single scattering albedo
!     and assymetry parameter
pi = 3.14159

do iz=nz,2,-1

!           calculate the single scattering albedo from both gasses
!           and particulates

      tau=tg(iz)+tp(iz)+tcr(iz)*chck
      if ( tau.eq.0.) then

         omg0 = 0.
         af   = 0.

      else

         omg0 = min(.999999,  &
          omgp(iz)*tp(iz)/max(1.e-20,tau))
         af   = asym(iz)

      endif

!           do delta-m scaling (wiscombe)

      fact        = af*af
      tau  = (1.-omg0*fact)*tau
      omg0 = ((1.-fact)*omg0)/(1.-omg0*fact)

!           determine the ODE matrix coefficients (Ritter and Geleyn)

      beta0     = (4.+af)/(8.*(1.+af))
      g1 = diffac*(1.-omg0*(1.-beta0))
      g2 = diffac*omg0*beta0
      gg = sqrt(g1**2-g2**2)

!           determine the local (true) reflection and transmission coefficients

      rinf  = g2/(gg+g1)
      ggtau = gg*tau
      expp1=1.-ggtau*(1.-ggtau*e2*(1.-ggtau*e3*(1.-ggtau*e4  &
             *(1.-ggtau*e5))))
      expp1=max(expp1,0.0)
      expp2=expp1**2
      t(iz) = (1.-rinf**2)*expp1 /  &
                  (1.-rinf**2*expp2)
      r(iz) = rinf*(1.-expp2) /  &
                  (1.-rinf**2*expp2)

!        get the source func, go from top down to accomodate solar terms

!            if (tau.lt..8e-2) then
      if (tau.lt.4.e-2) then      !changed June 12 after Jerry's recom.
            sigu(iz) = 0.5*pi*(src(iz)+src(iz-1))*tau*diffac
            sigd(iz) = sigu(iz)
      else
            aa = (g1+g2)*(1.-r(iz)) - (1.+r(iz)-t(iz))/tau
            bb = -(g1+g2)*t(iz) + (1.+r(iz)-t(iz))/tau
            cc = diffac*pi*(1.-omg0)/gg**2
            sigu(iz) = cc*(aa*src(iz)+bb*src(iz-1))
            sigd(iz) = cc*(bb*src(iz)+aa*src(iz-1))
            if (sigu(iz).lt.0. .or. sigd(iz).lt.0.)then
              print*,'NEGATIVE SOURCE IN LW RADIATION'
              print*,'Grid point (ib,k,i,j): ',ib,iz,i+mi0(ngrid),j+mj0(ngrid)
              print*,aa,bb,cc,g1,g1
              print*,tau,g1,g2,r(iz),t(iz)
              print*,src(iz),src(iz-1),sigu(iz),sigd(iz)
              if(sigu(iz) < 0.)sigu(iz)=0.
              if(sigd(iz) < 0.)sigd(iz)=0.
            endif
      endif

enddo

!     do adding
!     initialize

re(nz) = 0.
vd(nz) = 0.
fd(nz) = 0.
fu(1)  = pi*src(1)


do iz=nz,2,-1
      td(iz)   = 1. - re(iz)*r(iz)
      re(iz-1) = r(iz) + t(iz)**2*re(iz) /  &
                                      td(iz)
      vd(iz-1) = sigd(iz) + ( t(iz)*vd(iz) +  &
                    t(iz)*re(iz)*sigu(iz) ) /  &
                          td(iz)
      vu(iz-1) = ( r(iz)*vd(iz) + sigu(iz) ) /  &
                                      td(iz)
enddo


!     calculate fluxes going up through the layers

do iz=2,nz
      fd(iz-1) = re(iz-1)*fu(iz-1) + vd(iz-1)
      fu(iz)   = t(iz)*fu(iz-1)/td(iz)+vu(iz-1)
enddo

return
END SUBROUTINE flxlw

!##############################################################################
Subroutine radinit (ng,nb,ns,npsb,nuum,prf,alpha,trf,beta,xp,  &
    wght,wlenlo,wlenhi,solar,ralcs,a0,a1,a2,a3,exptabc,ulim,  &
    npartob,npartg,ncog,ncb,ocoef,bcoef,gcoef,gnu)

use mem_radiate, only:irce,rce_solc
    
implicit none

!     read input parameters for the radiative transfer scheme
!     Original: J Verlinde March 1993
!
!     Rewritten and modified: J.H. 1995/1996
!
integer, parameter :: mb=8,mg=3,mk=7
real, parameter :: top=1800.,tm=1800./296.
!
!     mb   --   maximum number of bands allowed
!     mg   --   maximum active gasses in scheme (different from
!               the max number of active gasses in any given band)
!     mk   --   maximum number of pseudobands allowed for any gas
!     (H2O - gas 1;  CO2 plus other gases - gas 2;  Ozone - gas 3)

! Saleeby(2008): For drizzle mode we would need to modify arrays
!   gnu,ocoef1,bcoef1,gcoef1,kkat and add to the tables below

integer :: ng,ns,nb,npsb(mg,mb),nuum(mb),npartob,npartg,ncog,ncb
real :: prf(mg,mb),alpha(mg,mk,mb),trf(mg,mb),beta(mg,mk,mb),  &
     xp(mg,mk,mb),wght(mg,mk,mb),  &
     wlenlo(mb),wlenhi(mb),solar(mb),ralcs(mb),  &
     a0(mb),a1(mb),a2(mb),a3(mb),solc,  &
     exptabc(150),ulim(mg,mb)
real :: ocoef(ncog,mb,npartob),bcoef(ncb,mb,npartob),  &
     gcoef(ncog,mb,npartg),gnu(7)

integer :: npsb1(3,8),nuum1(8)
real :: prf1(3,8),alpha1(3,7,8),trf1(3,8),beta1(3,7,8),  &
     xp1(3,7,8),wght1(3,7,8),  &
     wlenlo1(8),wlenhi1(8),solar1(8),ralcs1(8),  &
     a01(8),a11(8),a21(8),a31(8),  &
     ulim1(3,8)
real :: ocoef1(5,8,13,2),bcoef1(2,8,13,2)  &
    ,gcoef1(5,8,7,2)

integer :: ib,icog,icb,ig,ik,ngb,ip,icat,ignu,ic1,ic2,ic3,i,is
real :: temp
real, external :: sunavg

integer :: kkat(13)
data kkat /1,2,3,3,3,4,4,4,5,5,5,6,7/

data ((ocoef1(icog,ib, 1,1),icog=1,5),ib=1,8)/  &
   .53502E+00,  .25038E+00, -.28307E-02,  .17203E+00, -.22921E-01,  &
   .87670E+00, -.14653E-02, -.22618E+01,  .12308E+00, -.92743E-03,  &
   .99986E+00, -.49475E+00, -.19301E+01,  .60422E+00, -.21285E+01,  &
   .51950E+00, -.14330E+00, -.68641E-01, -.14036E+01, -.18877E+01,  &
   .51500E+00, -.75507E+00, -.25649E+00,  .16039E+01, -.19427E+01,  &
   .51014E+00, -.78594E+00, -.29105E+01,  .21041E+00, -.34407E-01,  &
   .50916E+00,  .68010E-01, -.12197E-01, -.19787E+01, -.25706E+01,  &
   .50477E+00,  .26895E+00, -.50744E-01,  .39204E-01, -.65051E-02/
data ((ocoef1(icog,ib, 2,1),icog=1,5),ib=1,8)/  &
   .50246E+00,  .27495E+00, -.86761E-02,  .15534E+00, -.11628E-02,  &
   .77184E+00,  .21225E+01, -.46616E+00,  .20555E+00, -.26433E-03,  &
   .99864E+00,  .29610E+01, -.76111E+00,  .24400E+01, -.21687E+01,  &
   .50373E+00, -.56477E+00, -.18834E+00,  .27339E-01, -.95556E-03,  &
   .50422E+00, -.10954E+01, -.32495E+00,  .26231E+01, -.12727E+01,  &
   .50285E+00,  .11969E+00, -.13137E-01,  .29679E+01, -.35296E+00,  &
   .50250E+00,  .96104E+00, -.36542E+00,  .70906E-01, -.63027E-02,  &
   .50217E+00,  .85769E+00, -.15454E+00,  .35989E-01, -.46899E-02/
data ((ocoef1(icog,ib, 3,1),icog=1,5),ib=1,8)/  &
   .81912E+00,  .89934E-01, -.95431E-02,  .77550E-01, -.64635E-01,  &
   .99919E+00,  .21984E+01, -.26508E+01, -.25080E+01, -.27903E+01,  &
   .99995E+00,  .30675E+00, -.24977E+01, -.40199E+00, -.27822E+01,  &
   .66707E+00,  .21413E+00, -.19915E-01, -.14337E+00, -.39737E+00,  &
   .56341E+00,  .21768E+00, -.47805E-01, -.39828E+00, -.15905E+01,  &
   .56926E+00, -.28098E+00, -.12872E+01,  .35146E+00, -.39604E-01,  &
   .55241E+00, -.21384E+00, -.80189E+00,  .18383E+00, -.36890E-01,  &
   .52773E+00,  .28583E+00, -.69522E-01,  .13389E+00, -.10619E-01/
data ((ocoef1(icog,ib, 4,1),icog=1,5),ib=1,8)/  &
   .88678E+00,  .95164E-01, -.18882E-01, -.79279E-01, -.26767E+01,  &
   .99955E+00, -.95345E+00, -.27709E+01,  .95466E+00, -.27609E+01,  &
   .99995E+00,  .30675E+00, -.24977E+01, -.40199E+00, -.27822E+01,  &
   .65337E+00, -.29702E+00, -.86659E+00,  .24877E+00, -.63270E-02,  &
   .57824E+00,  .23396E+00, -.29201E-01, -.97145E+00, -.20002E+01,  &
   .60869E+00, -.54567E+00, -.16368E+01,  .32509E+00, -.19304E-01,  &
   .57423E+00,  .20302E+00, -.20766E-01, -.52963E+00, -.12669E+01,  &
   .62105E+00,  .31371E+00, -.24231E-01, -.83159E+00, -.28613E+01/
data ((ocoef1(icog,ib, 5,1),icog=1,5),ib=1,8)/  &
   .85390E+00,  .24811E-01, -.59938E-01,  .97555E-01, -.13096E-01,  &
   .99945E+00,  .19216E+00, -.22428E+01, -.38142E+00, -.29731E+01,  &
   .99995E+00,  .30675E+00, -.24977E+01, -.40199E+00, -.27822E+01,  &
   .59890E+00, -.18424E+00, -.35464E+00,  .26931E+00, -.61489E-02,  &
   .56129E+00,  .20357E+00, -.27843E-01, -.52470E+00, -.13424E+01,  &
   .57703E+00, -.27255E+00, -.88195E+00,  .33391E+00, -.23252E-01,  &
   .54768E+00, -.29256E+00, -.61564E+00,  .17540E+00, -.19055E-01,  &
   .58753E+00,  .32487E+00, -.27846E-01, -.80520E+00, -.28833E+01/
data ((ocoef1(icog,ib, 6,1),icog=1,5),ib=1,8)/  &
   .60131E+00,  .14296E+01, -.26494E-01,  .19873E+00, -.25000E-03,  &
   .99253E+00,  .10755E+00, -.23677E-01, -.91790E+00, -.15857E+00,  &
   .99997E+00, -.27667E+01, -.85050E+00, -.15914E+01, -.28007E+01,  &
   .50331E+00,  .20589E+01, -.24104E+01,  .25282E+00, -.27455E-02,  &
   .50003E+00, -.18764E+01, -.39166E-01,  .12933E+00, -.45107E-02,  &
   .49990E+00, -.16069E+01, -.52634E-01,  .11701E+00, -.44343E-02,  &
   .49979E+00, -.23946E+01, -.22541E+01,  .91240E-01, -.41169E-02,  &
   .50012E+00,  .79855E+00, -.16299E+01,  .10765E+00, -.41424E-02/
data ((ocoef1(icog,ib, 7,1),icog=1,5),ib=1,8)/  &
   .61674E+00,  .21205E+01, -.29820E-01,  .23706E+00, -.22189E-03,  &
   .99432E+00,  .13598E+01, -.44456E-01,  .12953E+01, -.14623E+01,  &
   .99998E+00,  .21184E+01, -.92367E-01,  .53557E+00, -.35783E+00,  &
   .50773E+00,  .27353E+01, -.35146E-01,  .25352E+00, -.14164E-02,  &
   .49986E+00, -.23090E+01, -.14960E+01,  .13935E+00, -.41030E-02,  &
   .50042E+00,  .22300E+00, -.41742E-02, -.22590E+01, -.22655E+00,  &
   .50006E+00, -.24034E+01, -.19547E+01,  .14781E+00, -.41357E-02,  &
   .50123E+00, -.29943E+01, -.29163E+01,  .18482E+00, -.32265E-02/
data ((ocoef1(icog,ib, 8,1),icog=1,5),ib=1,8)/  &
   .60552E+00,  .19524E+00, -.91893E-02,  .22414E+00, -.22649E-03,  &
   .99345E+00, -.25853E+01, -.11440E+01,  .76242E+00, -.37469E-01,  &
   .99998E+00,  .96513E+00, -.84851E-01,  .29527E+01, -.11622E+01,  &
   .50469E+00,  .21314E+00, -.17656E-02,  .33216E+00, -.13907E-01,  &
   .50001E+00, -.10675E+01, -.52958E-01,  .11642E+00, -.42821E-02,  &
   .50003E+00, -.18789E+01, -.15510E+01,  .16530E+00, -.47971E-02,  &
   .50005E+00,  .11219E+00, -.44481E-02,  .15162E+00, -.11005E+01,  &
   .50059E+00,  .15061E+00, -.40776E-02, -.19571E+01, -.19943E+01/
data ((ocoef1(icog,ib, 9,1),icog=1,5),ib=1,8)/  &
   .59111E+00,  .19384E+00, -.20844E-03,  .15061E+00, -.50790E-02,  &
   .99242E+00,  .11224E+01, -.50262E+00, -.11549E+01, -.24816E+01,  &
   .99997E+00,  .10566E+01, -.12959E+01,  .29739E+01, -.11756E+01,  &
   .50337E+00,  .96222E+00, -.23622E+00,  .25860E+00, -.28177E-02,  &
   .49988E+00,  .10651E+00, -.42909E-02,  .17299E+01, -.28869E+00,  &
   .49997E+00,  .10887E+01, -.17114E+00,  .12077E+00, -.46694E-02,  &
   .49987E+00,  .94228E-01, -.43431E-02,  .54459E+00, -.17733E+00,  &
   .50023E+00,  .10820E+00, -.41047E-02,  .28487E+01, -.26901E+00/
data ((ocoef1(icog,ib,10,1),icog=1,5),ib=1,8)/  &
   .59885E+00,  .12438E+00, -.37180E-02,  .23437E+00, -.17364E-03,  &
   .99431E+00,  .98683E+00, -.50705E+00, -.15639E+01, -.17063E+01,  &
   .99998E+00,  .24721E+01, -.18065E+01,  .23966E+01, -.14163E+01,  &
   .50310E+00,  .79901E-01, -.51442E-03,  .30571E+00, -.34180E-02,  &
   .49967E+00,  .98947E+00, -.20676E+00,  .13659E+00, -.41175E-02,  &
   .50059E+00,  .25137E+00, -.39431E-01,  .22353E+00, -.42388E-02,  &
   .50008E+00,  .15606E+00, -.43760E-02,  .19032E+01, -.30191E+00,  &
   .50149E+00,  .19328E+00, -.33401E-02,  .16638E+01, -.22392E+00/
data ((ocoef1(icog,ib,11,1),icog=1,5),ib=1,8)/  &
   .59114E+00,  .22326E+00, -.18499E-03,  .13978E+00, -.47685E-02,  &
   .99306E+00,  .24021E+01, -.91627E+00,  .42397E-02, -.39950E-03,  &
   .99998E+00, -.27872E+01, -.17821E+01,  .27485E+01, -.11735E+01,  &
   .50579E+00,  .28701E+01, -.36958E+00,  .27189E+00, -.22118E-02,  &
   .50008E+00,  .27605E+01, -.32258E+00,  .11675E+00, -.44468E-02,  &
   .49984E+00,  .29173E+01, -.27598E+00,  .16851E+00, -.48549E-02,  &
   .50011E+00,  .13082E+00, -.50924E-02,  .10954E+01, -.28379E+00,  &
   .50050E+00,  .73313E+00, -.13738E+00,  .15749E+00, -.42053E-02/
data ((ocoef1(icog,ib,12,1),icog=1,5),ib=1,8)/  &
   .50246E+00,  .27495E+00, -.86761E-02,  .15534E+00, -.11628E-02,  &
   .77184E+00,  .21225E+01, -.46616E+00,  .20555E+00, -.26433E-03,  &
   .99864E+00,  .29610E+01, -.76111E+00,  .24400E+01, -.21687E+01,  &
   .50373E+00, -.56477E+00, -.18834E+00,  .27339E-01, -.95556E-03,  &
   .50422E+00, -.10954E+01, -.32495E+00,  .26231E+01, -.12727E+01,  &
   .50285E+00,  .11969E+00, -.13137E-01,  .29679E+01, -.35296E+00,  &
   .50250E+00,  .96104E+00, -.36542E+00,  .70906E-01, -.63027E-02,  &
   .50217E+00,  .85769E+00, -.15454E+00,  .35989E-01, -.46899E-02/
data ((ocoef1(icog,ib,13,1),icog=1,5),ib=1,8)/  &
   .50246E+00,  .27495E+00, -.86761E-02,  .15534E+00, -.11628E-02,  &
   .77184E+00,  .21225E+01, -.46616E+00,  .20555E+00, -.26433E-03,  &
   .99864E+00,  .29610E+01, -.76111E+00,  .24400E+01, -.21687E+01,  &
   .50373E+00, -.56477E+00, -.18834E+00,  .27339E-01, -.95556E-03,  &
   .50422E+00, -.10954E+01, -.32495E+00,  .26231E+01, -.12727E+01,  &
   .50285E+00,  .11969E+00, -.13137E-01,  .29679E+01, -.35296E+00,  &
   .50250E+00,  .96104E+00, -.36542E+00,  .70906E-01, -.63027E-02,  &
   .50217E+00,  .85769E+00, -.15454E+00,  .35989E-01, -.46899E-02/

data ((ocoef1(icog,ib, 1,2),icog=1,5),ib=1,8)/  &
   .52438E+00,  .24743E+00, -.34625E-02,  .18970E+00, -.29927E-01,  &
   .86366E+00,  .28636E-01, -.22157E+01,  .13490E+00, -.10851E-02,  &
   .99979E+00, -.29198E+01, -.28868E+01,  .18287E+01, -.24188E+01,  &
   .49952E+00,  .31959E-01, -.12977E-02, -.39239E+00, -.15836E+00,  &
   .50631E+00, -.43099E+00, -.21672E+00,  .29032E-01, -.42007E-02,  &
   .50446E+00,  .18046E+00, -.64406E-01,  .25188E-01, -.51384E-02,  &
   .50732E+00, -.20182E+01, -.27063E+01,  .61180E-01, -.13667E-01,  &
   .50425E+00,  .37006E-01, -.90129E-02,  .29626E+00, -.79197E-01/
data ((ocoef1(icog,ib, 2,2),icog=1,5),ib=1,8)/  &
   .50131E+00,  .27610E+00, -.10296E-01,  .13972E+00, -.14901E-02,  &
   .75380E+00,  .21442E+00, -.29393E-03,  .28470E+00, -.20896E+00,  &
   .99820E+00, -.29567E+01, -.74150E+00,  .13607E+01, -.59503E+00,  &
   .50496E+00, -.29518E+01, -.87080E+00, -.21751E+01, -.40022E+00,  &
   .50338E+00, -.29793E+01, -.93078E+00, -.24036E+01, -.51415E+00,  &
   .50201E+00,  .27620E+01, -.52448E+00,  .16007E+00, -.18536E-01,  &
   .50190E+00,  .23171E+01, -.36209E+00,  .19042E-01, -.18868E-02,  &
   .50138E+00,  .18239E+00, -.21161E-01,  .21794E+01, -.45835E+00/
data ((ocoef1(icog,ib, 3,2),icog=1,5),ib=1,8)/  &
   .81133E+00,  .69349E-01, -.10299E+00,  .10795E+00, -.14262E-01,  &
   .99895E+00, -.20536E+01, -.25147E+01,  .17417E+01, -.23481E+01,  &
   .99995E+00, -.40557E+00, -.26487E+01,  .36965E+00, -.25690E+01,  &
   .64454E+00,  .23063E+00, -.22302E-01, -.15667E+00, -.67213E+00,  &
   .55383E+00, -.70583E+00, -.26693E+01,  .22107E+00, -.59426E-01,  &
   .61802E+00, -.34675E-01,  .61173E-02,  .33720E+00, -.60327E-01,  &
   .54433E+00, -.19635E+00, -.91939E+00,  .19782E+00, -.48602E-01,  &
   .52284E+00,  .11900E+00, -.13524E-01,  .31394E+00, -.91911E-01/
data ((ocoef1(icog,ib, 4,2),icog=1,5),ib=1,8)/  &
   .87233E+00,  .99292E-01, -.19262E-01,  .17117E-01, -.19705E+00,  &
   .99953E+00, -.16322E+00, -.29089E+01,  .84448E-01, -.21449E+01,  &
   .99995E+00,  .30675E+00, -.24977E+01, -.40199E+00, -.27822E+01,  &
   .63934E+00,  .26529E+00, -.85349E-02, -.21627E+00, -.78589E+00,  &
   .56400E+00,  .24216E+00, -.36364E-01, -.14425E+01, -.29162E+01,  &
   .58293E+00,  .34643E+00, -.23927E-01, -.13853E+01, -.29851E+01,  &
   .55724E+00,  .21398E+00, -.24772E-01, -.80610E+00, -.20961E+01,  &
   .54695E+00,  .21108E+00, -.11528E-01,  .18701E+00, -.59256E-01/
data ((ocoef1(icog,ib, 5,2),icog=1,5),ib=1,8)/  &
   .84401E+00,  .23811E-01, -.14551E+00,  .10861E+00, -.15419E-01,  &
   .99932E+00, -.69472E+00, -.23245E+01,  .67102E+00, -.22771E+01,  &
   .99995E+00,  .30675E+00, -.24977E+01, -.40199E+00, -.27822E+01,  &
   .99963E+00, -.16618E+00,  .53578E-02, -.22109E+00, -.10289E+01,  &
   .55164E+00,  .21466E+00, -.36566E-01, -.60797E+00, -.18416E+01,  &
   .56015E+00, -.32017E+00, -.13440E+01,  .34671E+00, -.29540E-01,  &
   .54148E+00, -.25292E+00, -.72859E+00,  .18398E+00, -.25635E-01,  &
   .53004E+00,  .22726E+00, -.60634E-01,  .16822E+00, -.12088E-01/
data ((ocoef1(icog,ib, 6,2),icog=1,5),ib=1,8)/  &
   .56528E+00,  .20898E+00, -.23015E-03,  .18113E+00, -.66666E-02,  &
   .99010E+00,  .51630E-01, -.10120E-01, -.23299E+01, -.68731E-01,  &
   .99997E+00,  .16679E+01, -.88622E-01,  .18730E+01, -.63071E+00,  &
   .50148E+00,  .13349E+00, -.25006E+01,  .23413E+00, -.32827E-02,  &
   .49998E+00, -.28983E+01, -.40580E-01,  .12123E+00, -.52029E-02,  &
   .49984E+00,  .79700E-01, -.42625E-02, -.63772E+00, -.92501E+00,  &
   .49990E+00,  .73620E-01, -.41556E-02,  .22763E+01, -.24832E+01,  &
   .49997E+00,  .26563E+01, -.15949E+01,  .78663E-01, -.42985E-02/
data ((ocoef1(icog,ib, 7,2),icog=1,5),ib=1,8)/  &
   .56310E+00,  .25170E+00, -.17797E-03,  .88188E-01, -.24857E-02,  &
   .99255E+00,  .26672E+01, -.47333E-01, -.13995E+01, -.24663E+01,  &
   .99998E+00,  .18766E+01, -.18795E+01, -.85618E-01, -.19384E+01,  &
   .50248E+00,  .89183E-01, -.93480E-03,  .26905E+00, -.41338E-02,  &
   .50020E+00, -.11036E+01, -.30119E-01,  .18127E+00, -.59058E-02,  &
   .50014E+00, -.19260E+01, -.13113E+01,  .17751E+00, -.49112E-02,  &
   .50003E+00,  .11324E+00, -.45075E-02, -.17264E+00, -.12135E+01,  &
   .50038E+00,  .15476E+00, -.39471E-02,  .14488E+01, -.13307E+01/
data ((ocoef1(icog,ib, 8,2),icog=1,5),ib=1,8)/  &
   .58277E+00,  .24129E+00, -.25302E-03,  .29461E+01, -.32957E-01,  &
   .99189E+00,  .60466E+00, -.35096E-01, -.26453E+01, -.83121E+00,  &
   .99997E+00, -.27667E+01, -.85050E+00, -.15914E+01, -.28007E+01,  &
   .50333E+00,  .22535E+00, -.24322E-02,  .29507E+01, -.37696E-01,  &
   .49987E+00, -.19453E+01, -.29072E+01,  .89687E-01, -.42578E-02,  &
   .50009E+00, -.25865E+00, -.20818E+01,  .12573E+00, -.53493E-02,  &
   .49994E+00, -.99434E+00, -.13396E+01,  .85224E-01, -.46679E-02,  &
   .50013E+00,  .11770E+00, -.45465E-02, -.85177E+00, -.27648E+01/
data ((ocoef1(icog,ib, 9,2),icog=1,5),ib=1,8)/  &
   .56278E+00,  .20772E+00, -.21767E-03,  .15660E+00, -.52773E-02,  &
   .99044E+00, -.29442E+01, -.13689E+01,  .24217E+01, -.55470E+00,  &
   .99997E+00,  .96718E+00, -.25363E+01,  .24388E+01, -.11205E+01,  &
   .50139E+00,  .23373E+00, -.32921E-02,  .15407E+01, -.26574E+00,  &
   .49978E+00,  .84373E-01, -.38509E-02,  .62946E+00, -.19821E+00,  &
   .49997E+00,  .14601E+01, -.20115E+00,  .83662E-01, -.43680E-02,  &
   .49976E+00,  .25618E+01, -.32259E+00,  .71141E-01, -.41767E-02,  &
   .49958E+00,  .73392E-01, -.39063E-02,  .23186E+01, -.25438E+00/
data ((ocoef1(icog,ib,10,2),icog=1,5),ib=1,8)/  &
   .57280E+00,  .13298E+00, -.53572E-02,  .25483E+00, -.20239E-03,  &
   .99265E+00,  .11151E+00, -.22915E+01,  .29903E+01, -.60607E+00,  &
   .99998E+00, -.17678E+01, -.15604E+01,  .95765E+00, -.10485E+01,  &
   .50455E+00,  .19910E+01, -.28865E+00,  .27178E+00, -.20250E-02,  &
   .49994E+00,  .24198E+01, -.29927E+00,  .11348E+00, -.44506E-02,  &
   .50015E+00,  .27335E+00, -.33081E-01,  .17208E+00, -.49023E-02,  &
   .50019E+00,  .11259E+00, -.44642E-02,  .24168E+00, -.70501E-01,  &
   .50052E+00,  .15972E+00, -.40628E-02,  .47993E+00, -.93927E-01/
data ((ocoef1(icog,ib,11,2),icog=1,5),ib=1,8)/  &
   .56714E+00,  .14168E+00, -.60114E-02,  .23864E+00, -.20326E-03,  &
   .99203E+00,  .27487E+00, -.35208E+00, -.11985E+00, -.22339E+01,  &
   .99997E+00,  .28339E+01, -.12505E+01, -.28704E+01, -.29825E+01,  &
   .50298E+00,  .60429E+00, -.18340E+00,  .25075E+00, -.25888E-02,  &
   .50004E+00,  .86340E-01, -.44164E-02,  .22731E+01, -.29232E+00,  &
   .49986E+00,  .12787E+00, -.54055E-02,  .24263E+01, -.24682E+00,  &
   .49980E+00,  .85550E-01, -.43925E-02,  .10047E+01, -.22929E+00,  &
   .50024E+00,  .10162E+01, -.16810E+00,  .11731E+00, -.45331E-02/
data ((ocoef1(icog,ib,12,2),icog=1,5),ib=1,8)/  &
   .50131E+00,  .27610E+00, -.10296E-01,  .13972E+00, -.14901E-02,  &
   .75380E+00,  .21442E+00, -.29393E-03,  .28470E+00, -.20896E+00,  &
   .99820E+00, -.29567E+01, -.74150E+00,  .13607E+01, -.59503E+00,  &
   .50496E+00, -.29518E+01, -.87080E+00, -.21751E+01, -.40022E+00,  &
   .50338E+00, -.29793E+01, -.93078E+00, -.24036E+01, -.51415E+00,  &
   .50201E+00,  .27620E+01, -.52448E+00,  .16007E+00, -.18536E-01,  &
   .50190E+00,  .23171E+01, -.36209E+00,  .19042E-01, -.18868E-02,  &
   .50138E+00,  .18239E+00, -.21161E-01,  .21794E+01, -.45835E+00/
data ((ocoef1(icog,ib,13,2),icog=1,5),ib=1,8)/  &
   .50131E+00,  .27610E+00, -.10296E-01,  .13972E+00, -.14901E-02,  &
   .75380E+00,  .21442E+00, -.29393E-03,  .28470E+00, -.20896E+00,  &
   .99820E+00, -.29567E+01, -.74150E+00,  .13607E+01, -.59503E+00,  &
   .50496E+00, -.29518E+01, -.87080E+00, -.21751E+01, -.40022E+00,  &
   .50338E+00, -.29793E+01, -.93078E+00, -.24036E+01, -.51415E+00,  &
   .50201E+00,  .27620E+01, -.52448E+00,  .16007E+00, -.18536E-01,  &
   .50190E+00,  .23171E+01, -.36209E+00,  .19042E-01, -.18868E-02,  &
   .50138E+00,  .18239E+00, -.21161E-01,  .21794E+01, -.45835E+00/

data ((bcoef1(icb,ib, 1,1),icb=1,2),ib=1,8)/  &
   .34438E-11,  .19866E+01,  &
   .34388E-11,  .19863E+01,  &
   .33441E-11,  .19905E+01,  &
   .32585E-11,  .20009E+01,  &
   .31041E-11,  .20064E+01,  &
   .30466E-11,  .20084E+01,  &
   .29263E-11,  .20152E+01,  &
   .31930E-11,  .20001E+01/
data ((bcoef1(icb,ib, 2,1),icb=1,2),ib=1,8)/  &
   .32688E-11,  .19956E+01,  &
   .32122E-11,  .19975E+01,  &
   .31842E-11,  .19985E+01,  &
   .35364E-11,  .19873E+01,  &
   .33857E-11,  .19920E+01,  &
   .34428E-11,  .19898E+01,  &
   .33757E-11,  .19922E+01,  &
   .33885E-11,  .19916E+01/
data ((bcoef1(icb,ib, 3,1),icb=1,2),ib=1,8)/  &
   .58083E-11,  .17438E+01,  &
   .55947E-11,  .17484E+01,  &
   .52795E-11,  .17598E+01,  &
   .34851E-11,  .19093E+01,  &
   .54453E-11,  .17832E+01,  &
   .51652E-11,  .17924E+01,  &
   .41981E-11,  .18363E+01,  &
   .58427E-11,  .17551E+01/
data ((bcoef1(icb,ib, 4,1),icb=1,2),ib=1,8)/  &
   .18179E-10,  .14989E+01,  &
   .17374E-10,  .15048E+01,  &
   .16215E-10,  .15185E+01,  &
   .89107E-11,  .17122E+01,  &
   .16339E-10,  .15501E+01,  &
   .15479E-10,  .15601E+01,  &
   .11922E-10,  .16147E+01,  &
   .18206E-10,  .15130E+01/
data ((bcoef1(icb,ib, 5,1),icb=1,2),ib=1,8)/  &
   .48019E-11,  .17647E+01,  &
   .45784E-11,  .17712E+01,  &
   .42466E-11,  .17864E+01,  &
   .20760E-11,  .20058E+01,  &
   .41150E-11,  .18264E+01,  &
   .38476E-11,  .18393E+01,  &
   .29229E-11,  .18972E+01,  &
   .47067E-11,  .17835E+01/
data ((bcoef1(icb,ib, 6,1),icb=1,2),ib=1,8)/  &
   .17941E-10,  .14738E+01,  &
   .17785E-10,  .14748E+01,  &
   .17697E-10,  .14754E+01,  &
   .20203E-10,  .14601E+01,  &
   .18985E-10,  .14673E+01,  &
   .18692E-10,  .14691E+01,  &
   .18672E-10,  .14692E+01,  &
   .18378E-10,  .14710E+01/
data ((bcoef1(icb,ib, 7,1),icb=1,2),ib=1,8)/  &
   .10442E-10,  .16066E+01,  &
   .10358E-10,  .16075E+01,  &
   .10310E-10,  .16081E+01,  &
   .11679E-10,  .15936E+01,  &
   .11018E-10,  .16004E+01,  &
   .10855E-10,  .16021E+01,  &
   .10848E-10,  .16022E+01,  &
   .10687E-10,  .16039E+01/
data ((bcoef1(icb,ib, 8,1),icb=1,2),ib=1,8)/  &
   .16713E-11,  .19739E+01,  &
   .16595E-11,  .19747E+01,  &
   .16533E-11,  .19752E+01,  &
   .18340E-11,  .19632E+01,  &
   .17481E-11,  .19687E+01,  &
   .17252E-11,  .19703E+01,  &
   .17254E-11,  .19702E+01,  &
   .17030E-11,  .19717E+01/
data ((bcoef1(icb,ib, 9,1),icb=1,2),ib=1,8)/  &
   .15886E-10,  .14881E+01,  &
   .15600E-10,  .14903E+01,  &
   .15437E-10,  .14915E+01,  &
   .17739E-10,  .14753E+01,  &
   .17223E-10,  .14787E+01,  &
   .17031E-10,  .14800E+01,  &
   .16499E-10,  .14837E+01,  &
   .16659E-10,  .14826E+01/
data ((bcoef1(icb,ib,10,1),icb=1,2),ib=1,8)/  &
   .11227E-10,  .15981E+01,  &
   .10998E-10,  .16005E+01,  &
   .10864E-10,  .16019E+01,  &
   .12060E-10,  .15899E+01,  &
   .12088E-10,  .15895E+01,  &
   .11974E-10,  .15906E+01,  &
   .11453E-10,  .15958E+01,  &
   .11794E-10,  .15924E+01/
data ((bcoef1(icb,ib,11,1),icb=1,2),ib=1,8)/  &
   .18718E-11,  .19607E+01,  &
   .18382E-11,  .19628E+01,  &
   .18167E-11,  .19642E+01,  &
   .19564E-11,  .19556E+01,  &
   .19861E-11,  .19538E+01,  &
   .19683E-11,  .19549E+01,  &
   .18873E-11,  .19598E+01,  &
   .19525E-11,  .19558E+01/
data ((bcoef1(icb,ib,12,1),icb=1,2),ib=1,8)/  &
   .32688E-11,  .19956E+01,  &
   .32122E-11,  .19975E+01,  &
   .31842E-11,  .19985E+01,  &
   .35364E-11,  .19873E+01,  &
   .33857E-11,  .19920E+01,  &
   .34428E-11,  .19898E+01,  &
   .33757E-11,  .19922E+01,  &
   .33885E-11,  .19916E+01/
data ((bcoef1(icb,ib,13,1),icb=1,2),ib=1,8)/  &
   .32688E-11,  .19956E+01,  &
   .32122E-11,  .19975E+01,  &
   .31842E-11,  .19985E+01,  &
   .35364E-11,  .19873E+01,  &
   .33857E-11,  .19920E+01,  &
   .34428E-11,  .19898E+01,  &
   .33757E-11,  .19922E+01,  &
   .33885E-11,  .19916E+01/

data ((bcoef1(icb,ib, 1,2),icb=1,2),ib=1,8)/  &
   .10388E-10,  .19854E+01,  &
   .10167E-10,  .19885E+01,  &
   .98859E-11,  .19927E+01,  &
   .99564E-11,  .19966E+01,  &
   .95604E-11,  .20013E+01,  &
   .92513E-11,  .20056E+01,  &
   .89816E-11,  .20107E+01,  &
   .97704E-11,  .19962E+01/
data ((bcoef1(icb,ib, 2,2),icb=1,2),ib=1,8)/  &
   .96951E-11,  .19969E+01,  &
   .95759E-11,  .19982E+01,  &
   .95177E-11,  .19989E+01,  &
   .10482E-10,  .19885E+01,  &
   .10072E-10,  .19928E+01,  &
   .10166E-10,  .19916E+01,  &
   .10082E-10,  .19926E+01,  &
   .99843E-11,  .19936E+01/
data ((bcoef1(icb,ib, 3,2),icb=1,2),ib=1,8)/  &
   .20024E-10,  .16731E+01,  &
   .18971E-10,  .16826E+01,  &
   .18042E-10,  .16927E+01,  &
   .15730E-10,  .17736E+01,  &
   .20721E-10,  .16862E+01,  &
   .20019E-10,  .16905E+01,  &
   .16809E-10,  .17295E+01,  &
   .21236E-10,  .16694E+01/
data ((bcoef1(icb,ib, 4,2),icb=1,2),ib=1,8)/  &
   .44143E-10,  .15092E+01,  &
   .41297E-10,  .15213E+01,  &
   .38861E-10,  .15337E+01,  &
   .29552E-10,  .16484E+01,  &
   .44903E-10,  .15278E+01,  &
   .43528E-10,  .15316E+01,  &
   .34547E-10,  .15827E+01,  &
   .47297E-10,  .15043E+01/
data ((bcoef1(icb,ib, 5,2),icb=1,2),ib=1,8)/  &
   .12295E-10,  .17893E+01,  &
   .11473E-10,  .18021E+01,  &
   .10745E-10,  .18156E+01,  &
   .75080E-11,  .19480E+01,  &
   .12090E-10,  .18150E+01,  &
   .11598E-10,  .18212E+01,  &
   .90997E-11,  .18752E+01,  &
   .12992E-10,  .17872E+01/
data ((bcoef1(icb,ib, 6,2),icb=1,2),ib=1,8)/  &
   .42978E-10,  .14781E+01,  &
   .42685E-10,  .14789E+01,  &
   .42514E-10,  .14793E+01,  &
   .47449E-10,  .14666E+01,  &
   .45052E-10,  .14726E+01,  &
   .44365E-10,  .14744E+01,  &
   .44464E-10,  .14741E+01,  &
   .43801E-10,  .14759E+01/
data ((bcoef1(icb,ib, 7,2),icb=1,2),ib=1,8)/  &
   .26517E-10,  .16100E+01,  &
   .26350E-10,  .16107E+01,  &
   .26252E-10,  .16112E+01,  &
   .28994E-10,  .15997E+01,  &
   .27688E-10,  .16050E+01,  &
   .27286E-10,  .16067E+01,  &
   .27348E-10,  .16064E+01,  &
   .26969E-10,  .16080E+01/
data ((bcoef1(icb,ib, 8,2),icb=1,2),ib=1,8)/  &
   .47812E-11,  .19789E+01,  &
   .47558E-11,  .19795E+01,  &
   .47427E-11,  .19798E+01,  &
   .51322E-11,  .19707E+01,  &
   .49508E-11,  .19749E+01,  &
   .48910E-11,  .19763E+01,  &
   .49022E-11,  .19760E+01,  &
   .48465E-11,  .19773E+01/
data ((bcoef1(icb,ib, 9,2),icb=1,2),ib=1,8)/  &
   .39775E-10,  .14872E+01,  &
   .39221E-10,  .14888E+01,  &
   .38924E-10,  .14897E+01,  &
   .45126E-10,  .14725E+01,  &
   .42953E-10,  .14782E+01,  &
   .42418E-10,  .14797E+01,  &
   .41699E-10,  .14817E+01,  &
   .41391E-10,  .14825E+01/
data ((bcoef1(icb,ib,10,2),icb=1,2),ib=1,8)/  &
   .28126E-10,  .16031E+01,  &
   .27668E-10,  .16050E+01,  &
   .27422E-10,  .16061E+01,  &
   .30896E-10,  .15922E+01,  &
   .30342E-10,  .15943E+01,  &
   .30081E-10,  .15953E+01,  &
   .29174E-10,  .15989E+01,  &
   .29387E-10,  .15980E+01/
data ((bcoef1(icb,ib,11,2),icb=1,2),ib=1,8)/  &
   .52306E-11,  .19684E+01,  &
   .51549E-11,  .19701E+01,  &
   .51143E-11,  .19710E+01,  &
   .55947E-11,  .19606E+01,  &
   .55687E-11,  .19611E+01,  &
   .55284E-11,  .19620E+01,  &
   .53706E-11,  .19654E+01,  &
   .54315E-11,  .19640E+01/
data ((bcoef1(icb,ib,12,2),icb=1,2),ib=1,8)/  &
   .96951E-11,  .19969E+01,  &
   .95759E-11,  .19982E+01,  &
   .95177E-11,  .19989E+01,  &
   .10482E-10,  .19885E+01,  &
   .10072E-10,  .19928E+01,  &
   .10166E-10,  .19916E+01,  &
   .10082E-10,  .19926E+01,  &
   .99843E-11,  .19936E+01/
data ((bcoef1(icb,ib,13,2),icb=1,2),ib=1,8)/  &
   .96951E-11,  .19969E+01,  &
   .95759E-11,  .19982E+01,  &
   .95177E-11,  .19989E+01,  &
   .10482E-10,  .19885E+01,  &
   .10072E-10,  .19928E+01,  &
   .10166E-10,  .19916E+01,  &
   .10082E-10,  .19926E+01,  &
   .99843E-11,  .19936E+01/

data ((gcoef1(icog,ib, 1,1),icog=1,5),ib=1,8)/  &
   .85004E+00, -.23070E+01, -.29873E+01, -.72250E-01, -.13559E-01,  &
   .85583E+00, -.39371E-01, -.53225E-01, -.10528E+01, -.27121E+01,  &
   .89734E+00, -.36884E+00, -.28945E+00,  .29238E+00, -.46323E+00,  &
   .86046E+00, -.23999E+01, -.16198E+01, -.37815E+00, -.69269E-01,  &
   .91260E+00, -.19434E+00, -.77833E-01, -.27120E+01, -.14493E+01,  &
   .97134E+00, -.13436E+01, -.85240E+00, -.15603E+00, -.60289E-01,  &
   .96978E+00, -.27259E+01, -.14606E+01, -.13315E+00, -.55870E-01,  &
   .92327E+00, -.62450E+00, -.51258E+00, -.15262E+00, -.57748E-01/
data ((gcoef1(icog,ib, 2,1),icog=1,5),ib=1,8)/  &
   .85256E+00, -.28551E+01, -.11673E+01, -.87298E-01, -.12332E-01,  &
   .85579E+00, -.58954E+00, -.34927E+00, -.29288E+01, -.58065E+00,  &
   .89934E+00, -.25957E+01, -.51369E+00,  .16269E+01, -.95447E+00,  &
   .86130E+00, -.28143E+00, -.32093E-01,  .15132E+01, -.10311E+01,  &
   .91312E+00, -.11166E+01, -.24026E+00,  .24531E+01, -.64002E+00,  &
   .97157E+00, -.23646E+01, -.32910E+00,  .20974E+01, -.29852E+01,  &
   .97053E+00, -.12651E+00, -.44707E-01, -.27114E+01, -.18854E+01,  &
   .92363E+00, -.93127E+00, -.23277E+00,  .12683E+01, -.10108E+01/
data ((gcoef1(icog,ib, 3,1),icog=1,5),ib=1,8)/  &
   .84652E+00, -.14792E+00, -.49920E+00, -.82844E-01, -.32344E-01,  &
   .86490E+00, -.24280E-01, -.27182E-01, -.10035E+00, -.25304E+00,  &
   .90438E+00, -.57383E-01, -.15491E+00, -.10861E+00, -.94753E+00,  &
   .88805E+00, -.81787E+00, -.17892E+00, -.17723E+00, -.23314E-01,  &
   .89817E+00, -.18697E+00, -.74539E-01, -.94747E+00, -.42580E+00,  &
   .96632E+00, -.17709E+00, -.72026E-01, -.11097E+01, -.70518E+00,  &
   .96598E+00, -.15428E+00, -.70893E-01, -.11102E+01, -.59795E+00,  &
   .91847E+00, -.19007E+00, -.88016E-01, -.10278E+01, -.11070E+01/
data ((gcoef1(icog,ib, 4,1),icog=1,5),ib=1,8)/  &
   .86492E+00,  .14010E+01, -.18127E+01, -.12305E+01, -.34913E-01,  &
   .86413E+00, -.89078E-02, -.40056E+00, -.11315E+01, -.66225E+00,  &
   .90703E+00, -.23060E+01, -.14108E+00,  .16011E+01, -.10697E+01,  &
   .88834E+00, -.10417E+01, -.39147E+00, -.52226E+00, -.29662E-01,  &
   .90180E+00, -.22150E+01, -.26668E+01, -.20787E+01, -.52022E-01,  &
   .96995E+00, -.40797E+00, -.63701E-01, -.48035E+00, -.45402E-01,  &
   .97021E+00, -.22262E+00, -.38171E-01, -.27977E+01, -.75078E+00,  &
   .92271E+00, -.27473E+01, -.60188E-01, -.38875E+00, -.18987E+01/
data ((gcoef1(icog,ib, 5,1),icog=1,5),ib=1,8)/  &
   .86542E+00, -.10243E+00, -.14467E-01,  .18090E+01, -.88272E+00,  &
   .86609E+00, -.33357E+00, -.23570E+01, -.14838E+01, -.40147E+00,  &
   .90692E+00, -.21077E+01, -.67812E+00, -.20880E+01, -.51111E+00,  &
   .88844E+00,  .24039E+01, -.14473E+00, -.12302E+01, -.39473E-01,  &
   .90167E+00, -.18168E+01, -.29193E+00, -.79883E+00, -.43798E+00,  &
   .96988E+00, -.30375E+00, -.12232E+00,  .16400E+01, -.24678E+01,  &
   .97021E+00,  .17677E+01, -.66322E+00, -.16524E+01, -.29933E+00,  &
   .92259E+00, -.16496E+01, -.39920E+00, -.68006E+00, -.26008E+00/
data ((gcoef1(icog,ib, 6,1),icog=1,5),ib=1,8)/  &
   .86542E+00, -.10243E+00, -.14467E-01,  .18090E+01, -.88272E+00,  &
   .86609E+00, -.33357E+00, -.23570E+01, -.14838E+01, -.40147E+00,  &
   .90692E+00, -.21077E+01, -.67812E+00, -.20880E+01, -.51111E+00,  &
   .88844E+00,  .24039E+01, -.14473E+00, -.12302E+01, -.39473E-01,  &
   .90167E+00, -.18168E+01, -.29193E+00, -.79883E+00, -.43798E+00,  &
   .96988E+00, -.30375E+00, -.12232E+00,  .16400E+01, -.24678E+01,  &
   .97021E+00,  .17677E+01, -.66322E+00, -.16524E+01, -.29933E+00,  &
   .92259E+00, -.16496E+01, -.39920E+00, -.68006E+00, -.26008E+00/
data ((gcoef1(icog,ib, 7,1),icog=1,5),ib=1,8)/  &
   .86542E+00, -.10243E+00, -.14467E-01,  .18090E+01, -.88272E+00,  &
   .86609E+00, -.33357E+00, -.23570E+01, -.14838E+01, -.40147E+00,  &
   .90692E+00, -.21077E+01, -.67812E+00, -.20880E+01, -.51111E+00,  &
   .88844E+00,  .24039E+01, -.14473E+00, -.12302E+01, -.39473E-01,  &
   .90167E+00, -.18168E+01, -.29193E+00, -.79883E+00, -.43798E+00,  &
   .96988E+00, -.30375E+00, -.12232E+00,  .16400E+01, -.24678E+01,  &
   .97021E+00,  .17677E+01, -.66322E+00, -.16524E+01, -.29933E+00,  &
   .92259E+00, -.16496E+01, -.39920E+00, -.68006E+00, -.26008E+00/

data ((gcoef1(icog,ib, 1,2),icog=1,5),ib=1,8)/  &
   .85428E+00, -.22824E+00, -.32737E+00, -.57281E-01, -.11823E-01,  &
   .85677E+00, -.15857E+01, -.54104E+00,  .15231E+01, -.63857E+00,  &
   .89818E+00, -.25075E+00, -.24189E+01, -.36307E-01, -.12669E+00,  &
   .85962E+00, -.46391E+00, -.11261E+01, -.62972E+00, -.15502E+00,  &
   .91238E+00, -.15845E+01, -.18563E+01, -.26871E+00, -.14749E+00,  &
   .97117E+00, -.15075E+01, -.15209E+01, -.15929E+00, -.87209E-01,  &
   .96977E+00, -.15247E+00, -.91890E-01, -.29261E+01, -.27680E+01,  &
   .92298E+00, -.16593E+00, -.88094E-01, -.81450E+00, -.89494E+00/
data ((gcoef1(icog,ib, 2,2),icog=1,5),ib=1,8)/  &
   .85485E+00,  .16591E+01, -.13538E+01, -.10826E+00, -.20853E-01,  &
   .85877E+00, -.40050E+00, -.11467E+01, -.22559E+01, -.82045E+00,  &
   .89872E+00, -.48360E+00, -.88272E+00, -.27211E+01, -.10538E+01,  &
   .86142E+00, -.40265E+00, -.68366E-01, -.90545E+00, -.82210E+00,  &
   .91314E+00, -.28791E+00, -.96056E-01,  .16624E+01, -.78796E+00,  &
   .97158E+00, -.17553E+01, -.56189E+00, -.11158E+00, -.14263E+01,  &
   .97046E+00, -.15502E+00, -.10433E+00, -.14530E+01, -.95320E+00,  &
   .92357E+00, -.50976E-01, -.13639E+01, -.21887E+00, -.13027E+00/
data ((gcoef1(icog,ib, 3,2),icog=1,5),ib=1,8)/  &
   .83645E+00, -.10657E+00, -.88227E-01, -.57680E+00, -.29838E+01,  &
   .86231E+00, -.94684E-01, -.28902E+00, -.17472E-01, -.48864E-01,  &
   .90517E+00, -.12820E+00, -.97785E+00, -.30679E-01, -.13764E+00,  &
  -.10789E+00, -.88173E+00, -.17572E+00,  .92677E+00,  .76584E-03,  &
   .90332E+00, -.98943E+00, -.45205E+00, -.75741E-01, -.51938E-01,  &
   .96769E+00, -.14666E+00, -.86753E-01, -.11482E+01, -.95418E+00,  &
   .96824E+00, -.11279E+01, -.76324E+00, -.11148E+00, -.78466E-01,  &
   .91817E+00, -.18217E+00, -.12214E+00, -.11635E+01, -.16800E+01/
data ((gcoef1(icog,ib, 4,2),icog=1,5),ib=1,8)/  &
   .86761E+00, -.25078E+00, -.33773E-01, -.53825E+00, -.92074E+00,  &
   .86637E+00,  .22550E+01, -.19380E+00, -.13276E+01, -.93580E-01,  &
   .90660E+00,  .15513E+01, -.13785E+01, -.18898E+00, -.76385E+00,  &
   .88836E+00, -.24290E+01, -.76275E-01,  .14582E+01, -.16409E+01,  &
   .90184E+00, -.16807E+00, -.64505E-01, -.20434E+01, -.12036E+01,  &
   .96993E+00,  .15686E+01, -.15885E+01, -.87090E+00, -.91118E-01,  &
   .97029E+00, -.19457E+01, -.10327E+00,  .12140E+01, -.16934E+01,  &
   .92260E+00, -.20223E+01, -.12144E+01, -.13366E+01, -.10022E+00/
data ((gcoef1(icog,ib, 5,2),icog=1,5),ib=1,8)/  &
   .86750E+00, -.11043E+00, -.22250E-01, -.26821E+01, -.16937E+01,  &
   .86628E+00,  .40213E+00, -.13349E+01, -.29534E+01, -.86907E+00,  &
   .90633E+00,  .10585E+01, -.20870E+01, -.11491E+01, -.80802E+00,  &
   .88830E+00, -.12621E+00, -.29155E-01, -.46961E+00, -.13071E+00,  &
   .90172E+00, -.71294E+00, -.45464E+00, -.11305E+00, -.65795E-01,  &
   .96998E+00,  .18700E+01, -.22688E+01, -.15064E+00, -.63017E-01,  &
   .97035E+00, -.17586E+00, -.10184E+00,  .62649E+00, -.18503E+01,  &
   .92260E+00,  .10520E+01, -.14577E+01, -.19297E+00, -.11617E+00/
data ((gcoef1(icog,ib, 6,2),icog=1,5),ib=1,8)/  &
   .86750E+00, -.11043E+00, -.22250E-01, -.26821E+01, -.16937E+01,  &
   .86628E+00,  .40213E+00, -.13349E+01, -.29534E+01, -.86907E+00,  &
   .90633E+00,  .10585E+01, -.20870E+01, -.11491E+01, -.80802E+00,  &
   .88830E+00, -.12621E+00, -.29155E-01, -.46961E+00, -.13071E+00,  &
   .90172E+00, -.71294E+00, -.45464E+00, -.11305E+00, -.65795E-01,  &
   .96998E+00,  .18700E+01, -.22688E+01, -.15064E+00, -.63017E-01,  &
   .97035E+00, -.17586E+00, -.10184E+00,  .62649E+00, -.18503E+01,  &
   .92260E+00,  .10520E+01, -.14577E+01, -.19297E+00, -.11617E+00/
data ((gcoef1(icog,ib, 7,2),icog=1,5),ib=1,8)/  &
   .86750E+00, -.11043E+00, -.22250E-01, -.26821E+01, -.16937E+01,  &
   .86628E+00,  .40213E+00, -.13349E+01, -.29534E+01, -.86907E+00,  &
   .90633E+00,  .10585E+01, -.20870E+01, -.11491E+01, -.80802E+00,  &
   .88830E+00, -.12621E+00, -.29155E-01, -.46961E+00, -.13071E+00,  &
   .90172E+00, -.71294E+00, -.45464E+00, -.11305E+00, -.65795E-01,  &
   .96998E+00,  .18700E+01, -.22688E+01, -.15064E+00, -.63017E-01,  &
   .97035E+00, -.17586E+00, -.10184E+00,  .62649E+00, -.18503E+01,  &
   .92260E+00,  .10520E+01, -.14577E+01, -.19297E+00, -.11617E+00/

data ((ulim1(ig,ib),ig=1,3),ib=1,8)/  &
   1000.,   10000.,   1000.,  &
   10000.,  1000.,    1000.,  &
   1000.,   1000.,    2.,  &
   1000.,   1000.,    1000.,  &
   1000.,   100.,     1000.,  &
   1000.,   1000.,    1000.,  &
   1000.,   1000.,    0.3,  &
   1000.,   110.,     1000./
data (wlenlo1(ib),ib=1,8)/  &
   1.529999, .6999999, .2451000, 20.00000,  &
   12.50000, 8.333299, 9.008999, 4.642000/
data (wlenhi1(ib),ib=1,8)/  &
   4.642000, 1.529999, .6999999, 104.5149,  &
   20.00000, 9.008999, 10.30930, 8.333299/
data (solar1(ib),ib=1,3)/ 162.8590546, 556.7694702, 640.6521606/
data (ralcs1(ib),ib=1,8)/  &
  0.2306902E-09, 0.4273218E-08, 0.7142404E-07, 0.5187161E-04,  &
  0.4707304E-05, 0.8822867E-06, 0.5892014E-06, 0.4331474E-06/

!Saleeby(2009):Modify nuum1 statement to turn on water vapor continuum
!data (nuum1(ib),ib=1,8)/ 0, 0, 0, 0, 0, 1, 1, 0/
data (nuum1(ib),ib=1,8)/ 0, 0, 0, 1, 1, 1, 1, 1/

data (a01(ib),ib=1,8)/  &
    .0000000E+00,   .0000000E+00,   .0000000E+00,   .2380230E+01,  &
    .3292058E+02,   .1123225E+02,   .2146121E+01,  -.8898108E+02/
data (a11(ib),ib=1,8)/  &
    .0000000E+00,   .0000000E+00,   .0000000E+00,  -.1055970E+00,  &
   -.4862469E+00,  -.7751060E-01,   .2143994E-01,   .1392107E+01/
data (a21(ib),ib=1,8)/  &
    .0000000E+00,   .0000000E+00,   .0000000E+00,   .1012857E-02,  &
    .2017528E-02,  -.3353665E-03,  -.4981068E-03,  -.7344463E-02/
data (a31(ib),ib=1,8)/  &
    .0000000E+00,   .0000000E+00,   .0000000E+00,  -.8937388E-06,  &
   -.8957686E-06,   .2563034E-05,   .1815480E-05,   .1314578E-04/

data ((npsb1(ig,ib),ig=1,3),ib=1,8)/  &
    5,    3,    0,  &
    5,    1,    0,  &
    1,    1,    5,  &
    5,    0,    0,  &
    4,    4,    0,  &
    2,    1,    0,  &
    2,    1,    3,  &
    5,    1,    0/

data ((prf1(ig,ib),ig=1,3),ib=1,8)/  &
   .1013250E+06,   .8612625E+05,   .0000000E+00,  &
   .1013250E+06,   .8612625E+05,   .0000000E+00,  &
   .1013250E+06,   .8612625E+05,   .3039750E+04,  &
   .5066250E+05,   .0000000E+00,   .0000000E+00,  &
   .8612625E+05,   .1013250E+05,   .0000000E+00,  &
   .8612625E+05,   .5066250E+05,   .0000000E+00,  &
   .8612625E+05,   .5066250E+05,   .1013250E+05,  &
   .7599375E+05,   .5066250E+05,   .0000000E+00/

data ((trf1(ig,ib),ig=1,3),ib=1,8)/  &
   .2817000E+03,   .2558000E+03,   .0000000E+00,  &
   .2817000E+03,   .2558000E+03,   .0000000E+00,  &
   .2817000E+03,   .2558000E+03,   .2299000E+03,  &
   .2558000E+03,   .0000000E+00,   .0000000E+00,  &
   .2817000E+03,   .2299000E+03,   .0000000E+00,  &
   .2817000E+03,   .2558000E+03,   .0000000E+00,  &
   .2817000E+03,   .2558000E+03,   .2040000E+03,  &
   .2817000E+03,   .2558000E+03,   .0000000E+00/

data ((xp1(ig,ik,1),ig=1,3),ik=1,7)/  &
   .1247157E-01,   .5996504E-04,   .0000000E+00,  &
   .2623545E+01,   .8041044E-02,   .0000000E+00,  &
   .8128411E-05,   .6263048E+00,   .0000000E+00,  &
   .1484021E+00,   .0000000E+00,   .0000000E+00,  &
   .7244654E-03,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,2),ig=1,3),ik=1,7)/  &
   .4487467E-02,   .2326860E-04,   .0000000E+00,  &
   .5147366E+00,   .0000000E+00,   .0000000E+00,  &
   .5490468E-05,   .0000000E+00,   .0000000E+00,  &
   .5409482E-03,   .0000000E+00,   .0000000E+00,  &
   .3205971E-01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,3),ig=1,3),ik=1,7)/  &
   .1106090E-04,   .1005030E-04,   .1914260E+05,  &
   .0000000E+00,   .0000000E+00,   .5794290E+03,  &
   .0000000E+00,   .0000000E+00,   .7177940E+02,  &
   .0000000E+00,   .0000000E+00,   .1870680E+01,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,4),ig=1,3),ik=1,7)/  &
   .1790956E+00,   .0000000E+00,   .0000000E+00,  &
   .1755604E+01,   .0000000E+00,   .0000000E+00,  &
   .4080156E-02,   .0000000E+00,   .0000000E+00,  &
   .3145180E+02,   .0000000E+00,   .0000000E+00,  &
   .2784468E-01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,5),ig=1,3),ik=1,7)/  &
   .5217378E-03,   .6755403E+00,   .0000000E+00,  &
   .1272084E+01,   .3978131E-01,   .0000000E+00,  &
   .5748921E-02,   .3571308E-03,   .0000000E+00,  &
   .5702444E-01,   .1850978E+02,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,6),ig=1,3),ik=1,7)/  &
   .1298981E-03,   .5129330E-04,   .0000000E+00,  &
   .4006586E-02,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,7),ig=1,3),ik=1,7)/  &
   .2897531E-04,   .5129330E-04,   .4995610E+02,  &
   .1949597E-02,   .0000000E+00,   .4341897E-03,  &
   .0000000E+00,   .0000000E+00,   .1086045E+02,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((xp1(ig,ik,8),ig=1,3),ik=1,7)/  &
   .6476305E-03,   .1053610E-03,   .0000000E+00,  &
   .6949131E+00,   .0000000E+00,   .0000000E+00,  &
   .7677455E-02,   .0000000E+00,   .0000000E+00,  &
   .8363310E+01,   .0000000E+00,   .0000000E+00,  &
   .8868708E-01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/

data ((wght1(ig,ik,1),ig=1,3),ik=1,7)/  &
   .1415609E+00,   .7785574E+00,   .0000000E+00,  &
   .7056218E-01,   .1693635E+00,   .0000000E+00,  &
   .4375873E+00,   .5207906E-01,   .0000000E+00,  &
   .1332501E+00,   .0000000E+00,   .0000000E+00,  &
   .2170395E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,2),ig=1,3),ik=1,7)/  &
   .1096637E+00,   .1000000E+01,   .0000000E+00,  &
   .3161782E-01,   .0000000E+00,   .0000000E+00,  &
   .5323901E+00,   .0000000E+00,   .0000000E+00,  &
   .2418084E+00,   .0000000E+00,   .0000000E+00,  &
   .8451997E-01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,3),ig=1,3),ik=1,7)/  &
   .1000000E+01,   .1000000E+01,   .5640000E-03,  &
   .0000000E+00,   .0000000E+00,   .1086900E-01,  &
   .0000000E+00,   .0000000E+00,   .1243200E-01,  &
   .0000000E+00,   .0000000E+00,   .1844170E+00,  &
   .0000000E+00,   .0000000E+00,   .7917180E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,4),ig=1,3),ik=1,7)/  &
   .2591315E+00,   .0000000E+00,   .0000000E+00,  &
   .3134250E+00,   .0000000E+00,   .0000000E+00,  &
   .5277705E-01,   .0000000E+00,   .0000000E+00,  &
   .1492298E+00,   .0000000E+00,   .0000000E+00,  &
   .2254367E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,5),ig=1,3),ik=1,7)/  &
   .3616211E+00,   .1671412E+00,   .0000000E+00,  &
   .5159410E-01,   .2359013E+00,   .0000000E+00,  &
   .3750388E+00,   .5389559E+00,   .0000000E+00,  &
   .2117461E+00,   .5800171E-01,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,6),ig=1,3),ik=1,7)/  &
   .8236111E+00,   .1000000E+01,   .0000000E+00,  &
   .1763889E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,7),ig=1,3),ik=1,7)/  &
   .8024495E+00,   .1000000E+01,   .1427190E+00,  &
   .1975506E+00,   .0000000E+00,   .5857329E+00,  &
   .0000000E+00,   .0000000E+00,   .2726962E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((wght1(ig,ik,8),ig=1,3),ik=1,7)/  &
   .1278349E+00,   .1000000E+01,   .0000000E+00,  &
   .2494413E+00,   .0000000E+00,   .0000000E+00,  &
   .2109432E+00,   .0000000E+00,   .0000000E+00,  &
   .1522931E+00,   .0000000E+00,   .0000000E+00,  &
   .2594875E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/

data ((alpha1(ig,ik,1),ig=1,3),ik=1,7)/  &
   .8242220E+00,   .9999895E+00,   .0000000E+00,  &
   .6627719E+00,   .5184072E+00,   .0000000E+00,  &
   .9999794E+00,   .6601071E+00,   .0000000E+00,  &
   .7781165E+00,   .0000000E+00,   .0000000E+00,  &
   .7519682E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,2),ig=1,3),ik=1,7)/  &
   .7693452E+00,   .1385630E+00,   .0000000E+00,  &
   .5131682E+00,   .0000000E+00,   .0000000E+00,  &
   .3494183E+00,   .0000000E+00,   .0000000E+00,  &
   .7714978E+00,   .0000000E+00,   .0000000E+00,  &
   .6381509E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,3),ig=1,3),ik=1,7)/  &
  -.2816690E+00,   .1000000E+01,  -.2044000E-02,  &
   .0000000E+00,   .0000000E+00,   .7768400E-01,  &
   .0000000E+00,   .0000000E+00,  -.2296670E+00,  &
   .0000000E+00,   .0000000E+00,   .9945000E-01,  &
   .0000000E+00,   .0000000E+00,   .1000000E+01,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,4),ig=1,3),ik=1,7)/  &
   .9998096E+00,   .0000000E+00,   .0000000E+00,  &
   .8934920E+00,   .0000000E+00,   .0000000E+00,  &
   .4805090E+00,   .0000000E+00,   .0000000E+00,  &
   .7776389E+00,   .0000000E+00,   .0000000E+00,  &
   .7963042E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,5),ig=1,3),ik=1,7)/  &
   .2766540E+00,   .7073100E+00,   .0000000E+00,  &
   .5560525E+00,   .7701662E+00,   .0000000E+00,  &
   .6394594E+00,   .2286618E+00,   .0000000E+00,  &
   .9728559E+00,   .4597031E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,6),ig=1,3),ik=1,7)/  &
   .9998614E+00,   .7833650E+00,   .0000000E+00,  &
   .4105807E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,7),ig=1,3),ik=1,7)/  &
   .9986984E+00,   .2954650E+00,   .2158372E+00,  &
   .3831264E+00,   .0000000E+00,   .3950504E-06,  &
   .0000000E+00,   .0000000E+00,   .4674741E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((alpha1(ig,ik,8),ig=1,3),ik=1,7)/  &
   .9913880E+00,   .5113000E-02,   .0000000E+00,  &
   .8563318E+00,   .0000000E+00,   .0000000E+00,  &
   .9979001E+00,   .0000000E+00,   .0000000E+00,  &
   .7047079E+00,   .0000000E+00,   .0000000E+00,  &
   .8926126E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/

data ((beta1(ig,ik,1),ig=1,3),ik=1,7)/  &
  -.1712392E+01,  -.3158334E+00,   .0000000E+00,  &
  -.3852223E+00,  -.5670782E+00,   .0000000E+00,  &
  -.1378367E+01,  -.7933998E+00,   .0000000E+00,  &
  -.7341006E+00,   .0000000E+00,   .0000000E+00,  &
  -.1588977E+01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,2),ig=1,3),ik=1,7)/  &
  -.3969185E+00,   .4752930E+00,   .0000000E+00,  &
  -.3413553E+00,   .0000000E+00,   .0000000E+00,  &
  -.1065748E+01,   .0000000E+00,   .0000000E+00,  &
  -.5795957E+00,   .0000000E+00,   .0000000E+00,  &
  -.3413572E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,3),ig=1,3),ik=1,7)/  &
   .4186570E+00,  -.1006700E+01,   .4999120E+00,  &
   .0000000E+00,   .0000000E+00,   .4854630E+00,  &
   .0000000E+00,   .0000000E+00,   .4645810E+00,  &
   .0000000E+00,   .0000000E+00,  -.2546340E+00,  &
   .0000000E+00,   .0000000E+00,   .1000000E+01,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,4),ig=1,3),ik=1,7)/  &
  -.1718928E+01,   .0000000E+00,   .0000000E+00,  &
  -.1005695E+01,   .0000000E+00,   .0000000E+00,  &
   .9999878E+00,   .0000000E+00,   .0000000E+00,  &
  -.7752110E+00,   .0000000E+00,   .0000000E+00,  &
  -.3295007E+01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,5),ig=1,3),ik=1,7)/  &
  -.2055048E+01,  -.2231248E+01,   .0000000E+00,  &
  -.4119874E+01,  -.4965271E+01,   .0000000E+00,  &
  -.2822644E+01,  -.1611550E+01,   .0000000E+00,  &
  -.3811633E+01,  -.1467305E+01,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,6),ig=1,3),ik=1,7)/  &
  -.6162320E+01,  -.4473330E+01,   .0000000E+00,  &
  -.2685515E+01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,7),ig=1,3),ik=1,7)/  &
  -.9999836E+01,  -.5629570E+01,   .6873658E-01,  &
  -.2214210E+01,   .0000000E+00,  -.1255250E+01,  &
   .0000000E+00,   .0000000E+00,  -.1283144E+01,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/
data ((beta1(ig,ik,8),ig=1,3),ik=1,7)/  &
  -.4327339E-01,   .4952220E+00,   .0000000E+00,  &
  -.2665102E+01,   .0000000E+00,   .0000000E+00,  &
  -.6325287E+01,   .0000000E+00,   .0000000E+00,  &
  -.1731381E+01,   .0000000E+00,   .0000000E+00,  &
  -.4994114E+01,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00,  &
   .0000000E+00,   .0000000E+00,   .0000000E+00/

!Default solar constant
solc = 1370.

!Solar Constant for RCE simulations
if(irce == 1) solc = rce_solc

!Array/table dimensions
nb = 8
ng = 3
ns = 3
ngb = 2

!Loop over hydrometeor types (with habits) (Saleeby: modify for drizzle)
do ip = 1,13
   icat = kkat(ip)
   ignu = max(1,min(2,nint(gnu(icat))))
   do ib = 1,8
      do ic1 = 1,5
         ocoef(ic1,ib,ip) = ocoef1(ic1,ib,ip,ignu)
      enddo
      do ic2 = 1,2
         bcoef(ic2,ib,ip) = bcoef1(ic2,ib,ip,ignu)
      enddo
   enddo
enddo

!Loop over hydrometeor types (without habits) (Saleeby: modify for drizzle)
do ip = 1,7
   ignu = max(1,min(2,nint(gnu(ip))))
   do ib = 1,8
      do ic3 = 1,5
         gcoef(ic3,ib,ip) = gcoef1(ic3,ib,ip,ignu)
      enddo
   enddo
enddo

do ib = 1,8
   do ig = 1,3
     ulim(ig,ib) = ulim1(ig,ib)
     prf(ig,ib) = prf1(ig,ib)
     trf(ig,ib) = trf1(ig,ib)
     npsb(ig,ib) = npsb1(ig,ib)
   enddo
enddo

do i = 1,nb
   wlenlo(i) = wlenlo1(i)
   wlenhi(i) = wlenhi1(i)
   ralcs(i) = ralcs1(i)
   nuum(i) = nuum1(i)
   a0(i) = a01(i)
   a1(i) = a11(i)
   a2(i) = a21(i)
   a3(i) = a31(i)
enddo

do is = 1,ns
   solar(is) = solar1(is)
enddo

do ib = 1,8
   do ig = 1,3
      do ik = 1,7
         xp(ig,ik,ib) = xp1(ig,ik,ib)
         wght(ig,ik,ib) = wght1(ig,ik,ib)
         alpha(ig,ik,ib) = alpha1(ig,ik,ib)
         beta(ig,ik,ib) = beta1(ig,ik,ib)
      enddo
   enddo
enddo

!     do the solar, rayleigh scatter, and continuum abs.
do ib=1,nb

   if (ib.le.ns) then
      solar(ib) = sunavg(1.e4/wlenhi(ib),1.e4/wlenlo(ib),solc)
      CALL rayleigh (wlenlo(ib),wlenhi(ib),ralcs(ib))
   else
      CALL csband (wlenlo(ib),wlenhi(ib),ralcs(ib))
   endif

enddo
!
!     make the tables for the exponential func
!
do i=1,150
 temp=float(i)+179.
 exptabc(i) = exp(top/temp-tm)
enddo

return
END SUBROUTINE radinit

!##############################################################################
Subroutine rayleigh (wlnlo,wlnhi,rayavg)

implicit none

real :: wlnlo,wlnhi,rayavg
real :: an0,ssum,h1,sum,wl1,f1,fac,f2,h2,wl2,al
real, external :: sunirr
integer :: num,i

!     this stuff comes from Houghton 1985 (book), see also Slingo
!     and Schrecker 1982 (QJ)
!     calculate constant part (flux weighted) of rayleigh scattering
!     rayleigh scattering = rayavg*delta z*press/temp
!     written JV May 1993

real, parameter :: an01=1.000064328,an02=0.0294981,an03=.2554e-3
real, parameter :: an0d1=146.,an0d2=41.

an0(al) = an01+an02/(an0d1-al**2)+an03/(an0d2-al**2)

num = int(wlnhi-wlnlo)+1
num = min(max(25,num),5000)
sum = 0.
ssum= 0.
wl1 = wlnlo
f1  = (an0(wl1)**2-1.)**2/wl1**4*sunirr(1.e4/wl1)
h1  = sunirr(1.e4/wl1)
do i=2,num
   fac = (real(i)-1)/(real(num)-1.)
   wl2 = wlnlo+(wlnhi-wlnlo)*fac
   f2  = (an0(wl2)**2-1.)**2/wl2**4*sunirr(1.e4/wl2)
   h2  = sunirr(1.e4/wl2)
   sum = sum + 0.5*(f1+f2)*(wl2-wl1)
   ssum= ssum + 0.5*(h1+h2)*(wl2-wl1)
   wl1 = wl2
   f1  = f2
   h1  = h2
enddo
rayavg = 0.945319e-2*sum/ssum

return
END SUBROUTINE rayleigh

!##############################################################################
Subroutine csband (wlnlo,wlnhi,csavg)

implicit none

real :: wlnlo,wlnhi,csavg
real :: cs,wnhi,wnlo,sum,wn1,f1,fac,f2,wn,wn2
real, external :: plkavg
integer :: num,i
!
!     calculate the blackbody flux (t=296.) weighted self-
!     broadening coefficient for water vapor. See Kniezys
!     et al 1980 (LOWTRAN 5). Units are converted to have
!     the water vapor content input a Pascal.

cs(wn)  = 4.18+5578.*exp(-7.87e-3*wn)

wnhi = 1.e4/wlnlo
wnlo = 1.e4/wlnhi
num  = int(wnhi-wnlo)+1
num  = min(max(25,num),5000)
sum  = 0.
wn1  = wnlo
f1   = cs(wn1)
do i = 2,num
   fac = (real(i)-1)/(real(num)-1.)
   wn2 = wnlo+(wnhi-wnlo)*fac
   f2  = cs(wn2)
   sum = sum + 0.5*(f1+f2) * plkavg(wn1,wn2,296.)
   wn1 = wn2
   f1  = f2
enddo
csavg = sum/(plkavg(wnlo,wnhi,296.) *1013250.*9.81)

return
END SUBROUTINE csband

!##############################################################################
real Function plkavg (WNUMLO,WNUMHI,T)

implicit none

!      COMPUTES PLANCK FUNC INTEGRATED BETWEEN TWO WAVENUMBERS
!
!      I N P U T
!
!         WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL INTERVAL
!         WNUMHI : UPPER WAVENUMBER
!         T       : TEMPERATURE (K)
!
!      O U T P U T
!
!         PLKAVG : INTEGRATED PLANCK FUNC ( WATTS/SQ M )
!
!      REFERENCES-- (1) HOUGHTON,PHYSICS OF ATMOSPHERES,APPENDIX 7
!                   (2) SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
!                       OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.
!                       JAN. 1974
!
!      METHOD-- HOUGHTON'S EXPONENTIAL SERIES IS USED FOR V.GT.VCUT
!               ( 'V' IS HOUGHTON'S NOTATION ) AND HIS POWER SERIES
!               IN V FOR V.LE.VCUT.  MORE TERMS ARE TAKEN IN THE
!               EXPONENTIAL SERIES, THE LARGER V IS.  ( NOTE THAT
!               HOUGHTON'S ASSESSMENT THAT THE POWER SERIES IS USEFUL
!               FOR  V.LT.2*PI  IS INCORRECT--VCUT MUST BE LESS THAN
!               2 JUST IN ORDER TO GET 4-5 SIGNIFICANT DIGITS. )
!
!      ACCURACY-- 6 SIGNIFICANT DIGITS
!
!          *** ARGUMENTS
REAL     T, WNUMLO, WNUMHI
!          *** LOCAL VARIABLES
!                 C2       :  SECOND RADIATION CONSTANT
!                 CONA     :  15 / PI**4 / 13305600
!                 CONB     :  440 / 9
!                 D        :  INTEGRAL OF NORMALIZED PLANCK FUNC
!                             FROM 0 TO CURRENT WAVENUMBER
!                 EX       :  EXP( - V )
!                 F15PI4   :  15 / PI**4
!                 MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES
!                 MV       :  MULTIPLES OF *V*
!                 SIGDPI   :  STEFAN-BOLTZMANN CONSTANT DIVIDED BY PI
!                 V        :  H*C*NU / (K*T), WHERE H=PLANCKS CONSTANT,
!                             C=SPEED OF LIGHT, K=BOLTZMANN CONSTANT,
!                             NU=WAVENUMBER
!                 VCUT     :  POWER-SERIES CUTOFF POINT
!                 VCP      :  EXPONENTIAL SERIES CUTOFF POINTS
REAL     C2, CONA, CONB, D( 2 ), EX, F15PI4, MV, SIGDPI,  &
         V, VCUT, VCP( 7 ), VSQ
DATA  C2     / 1.438786    /,  CONA    / 0.11573303E-7 /,  &
      CONB   / 48.888889   /,  F15PI4 / 0.15398973    /,  &
      SIGDPI / 1.804919E-8 /,  VCUT    / 1.5 /,  &
      VCP    / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      
integer :: i,mmax,m

IF( T.LE.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LE.0. ) stop 'ahah'
V = 0.0 !Variable initialized
DO  50  I = 1, 2
   IF( I.EQ.1 )  V = ( C2 / T ) * WNUMLO
   IF( I.EQ.2 )  V = ( C2 / T ) * WNUMHI
   IF( V.LT.VCUT )  THEN

!            *** HOUGHTON'S POWER SERIES (FACTORED USING HORNER'S RULE)
      VSQ = V**2
      D( I ) = 1.0 - CONA *VSQ * V * ( 4435200. + V * ( -1663200.  &
               + V * ( 221760. + VSQ * ( -2640. + VSQ * ( CONB  &
               - VSQ ) ) ) ) )
   ELSE

!            *** HOUGHTON'S EXPONENTIAL SERIES
      MMAX = 0
!            *** SET UPPER LIMIT OF SERIES DEPENDING ON VALUE OF V
20       MMAX = MMAX + 1
         IF ( V.LT.VCP( MMAX ) )  GO TO 20

      EX = EXP( -V )
      D( I ) = EX * ( 6. + V * ( 6. + V * ( 3. + V ) ) )

      DO  30  M = 2, MMAX
         MV = M * V
         D( I ) = D( I ) + EX**M * ( 6. + MV * ( 6. + MV *  &
                           ( 3. + MV ) ) ) / M**4
30       CONTINUE

      D( I ) = F15PI4 * D( I )
   END IF

50 CONTINUE

plkavg = SIGDPI * T**4 * ( D( 1 ) - D( 2 ) )

return
END FUNCTION plkavg

!##############################################################################
real Function sunavg (wnumlo,wnumhi,solcon)

implicit none

real :: wnumlo, wnumhi, solcon
real :: scln=1372.844
integer :: num,i
real :: v1,s1,v2,s2,fac2,sum
real, external :: sunirr

! Spectral solar flux  (after LOWTRAN7)
!
! On input:
! wnumlo --- lower wavelength of a spectral interval (wavenumbers)
! wnumhi --- upper wavelength of a spectral interval (wavenumbers)
! solcon - value of solar "constant" in Watts/ sq m
!
! On output:
! Total extraterrestrial solar flux (W/sq m) between
! wavelengths wnumlo, wnumhi
! Trapezoidal integration is used with the resolution 1 cm^-1
! or such that there is at least 25 points between wnumhi and
! wnumlo but no more than 5000 points.
!
num = int(wnumhi-wnumlo)+1
num = min(max(25, num),5000)
sum = 0.
v1 = wnumlo
s1 = sunirr(v1)
do i=2, num
  fac2 = (real(i)-1.)/(real(num)-1.)
  v2 =  wnumlo+(wnumhi-wnumlo)*fac2
  s2 = sunirr(v2)
  sum = sum + (s1+s2)/2.*(10000./v1-10000./v2)
  v1 = v2
  s1 = s2
enddo
sunavg = (solcon/scln)*sum

return
END FUNCTION sunavg

!##############################################################################
real Function sunirr (v)

implicit none

real, intent(in) :: v

!        EVALUATES THE EXTRA-TERRESTRIAL SOLAR IRRADIANCE
!
!        INPUT:  V  =  FREQUENCY   (CM-1)
!                VALID RANGE   0 TO 57490 (CM-1)
!                    (EQUIVALENT TO WAVELENGTHS > 0.174 MICROMETERS)
!
!        OUTPUT:  SUNIRR  =  SOLAR IRRADIANCE  (WATTS M-2 MICROMETER-1)
!
!        WRITES A WARNING MESSAGE TO TAPE6  &  RETURNS  SUNIRR = 0
!            IF THE INPUT FREQUENCY IS OUT OF RANGE

real, parameter ::  a=3.50187e-13, b=3.93281
real :: wm,w0,w1,w2,p
integer :: i
real, save :: solara(1440), solarb(2910)

! SOLAR SPECTRUM FROM      0 TO    800 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=1,41) /                                    &
 0.0000E+00, 4.5756E-08, 7.0100E-07, 3.4580E-06, 1.0728E-05,  &
 2.5700E-05, 5.2496E-05, 9.6003E-05, 1.6193E-04, 2.5766E-04,  &
 3.9100E-04, 5.6923E-04, 8.0203E-04, 1.1006E-03, 1.4768E-03,  &
 1.9460E-03, 2.5213E-03, 3.2155E-03, 4.0438E-03, 5.0229E-03,  &
 6.1700E-03, 7.5145E-03, 9.0684E-03, 1.0853E-02, 1.2889E-02,  &
 1.5213E-02, 1.7762E-02, 2.0636E-02, 2.3888E-02, 2.7524E-02,  &
 3.1539E-02, 3.5963E-02, 4.0852E-02, 4.6236E-02, 5.2126E-02,  &
 5.8537E-02, 6.5490E-02, 7.3017E-02, 8.1169E-02, 9.0001E-02,  &
 9.9540E-02 /

! SOLAR SPECTRUM FROM    820 TO   3680 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=42,185) /                                      &
 .10980, .12080, .13260, .14520, .15860, .17310, .18850, .20490,  &
 .22240, .24110, .26090, .28200, .30430, .32790, .35270, .37890,  &
 .40650, .43550, .46600, .49800, .53160, .56690, .60390, .64260,  &
 .68320, .72560, .76990, .81620, .86440, .91470, .96710, 1.0220,  &
 1.0780, 1.1370, 1.1990, 1.2630, 1.3290, 1.3990, 1.4710, 1.5460,  &
 1.6250, 1.7060, 1.7910, 1.8800, 1.9710, 2.0670, 2.1660, 2.2680,  &
 2.3740, 2.4840, 2.5970, 2.7140, 2.8350, 2.9600, 3.0890, 3.2210,  &
 3.3570, 3.4980, 3.6420, 3.7900, 3.9440, 4.1040, 4.2730, 4.4450,  &
 4.6150, 4.7910, 4.9830, 5.1950, 5.4210, 5.6560, 5.8930, 6.1270,  &
 6.3560, 6.5820, 6.8080, 7.0360, 7.2700, 7.5170, 7.7890, 8.0910,  &
 8.4070, 8.7120, 8.9900, 9.2490, 9.5000, 9.7550, 10.010, 10.250,  &
 10.480, 10.700, 10.950, 11.230, 11.550, 11.900, 12.250, 12.600,  &
 12.930, 13.250, 13.530, 13.780, 14.040, 14.320, 14.660, 15.070,  &
 15.530, 16.011, 16.433, 16.771, 17.077, 17.473, 17.964, 18.428,  &
 18.726, 18.906, 19.141, 19.485, 19.837, 20.160, 20.509, 21.024,  &
 21.766, 22.568, 23.190, 23.577, 23.904, 24.335, 24.826, 25.236,  &
 25.650, 26.312, 27.208, 27.980, 28.418, 28.818, 29.565, 30.533,  &
 31.247, 31.667, 32.221, 33.089, 33.975, 34.597, 35.004, 35.395 /

! SOLAR SPECTRUM FROM   3700 TO   6560 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=186,329) /                                     &
 36.026, 36.985, 37.890, 38.401, 38.894, 39.857, 40.926, 41.570,  &
 42.135, 43.083, 44.352, 45.520, 45.982, 46.281, 48.335, 51.987,  &
 54.367, 54.076, 52.174, 50.708, 52.153, 55.707, 56.549, 54.406,  &
 53.267, 56.084, 61.974, 64.406, 60.648, 55.146, 53.067, 57.476,  &
 64.645, 68.348, 69.055, 69.869, 70.943, 71.662, 72.769, 74.326,  &
 75.257, 74.883, 73.610, 73.210, 74.886, 78.042, 80.204, 80.876,  &
 82.668, 84.978, 86.244, 88.361, 91.998, 95.383, 98.121, 100.29,  &
 100.64, 99.997, 101.82, 105.06, 107.50, 109.99, 112.45, 113.90,  &
 113.79, 119.23, 121.96, 124.58, 127.14, 125.19, 124.37, 125.00,  &
 127.88, 130.67, 131.98, 133.74, 136.69, 136.18, 135.02, 137.44,  &
 138.44, 137.25, 136.35, 142.60, 144.54, 148.37, 151.90, 151.55,  &
 155.35, 157.59, 159.70, 162.28, 168.44, 171.43, 169.82, 170.33,  &
 172.28, 176.68, 181.92, 186.06, 187.85, 186.00, 189.82, 189.35,  &
 192.86, 202.00, 209.63, 205.76, 212.88, 215.63, 216.51, 219.20,  &
 220.29, 221.12, 227.12, 229.97, 233.23, 233.95, 234.52, 234.45,  &
 235.77, 239.80, 243.11, 241.19, 242.34, 243.69, 242.84, 246.19,  &
 246.11, 246.76, 251.75, 255.38, 258.74, 260.26, 263.40, 268.68,  &
 271.81, 272.95, 273.93, 274.74, 274.43, 279.69, 287.76, 287.72 /

! SOLAR SPECTRUM FROM   6580 TO   9440 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=330,473) /                                     &
 287.96, 290.01, 291.92, 295.28, 296.78, 300.46, 302.19, 299.14,  &
 301.43, 305.68, 309.29, 310.63, 313.24, 314.61, 309.58, 318.81,  &
 320.54, 321.62, 328.58, 331.66, 337.20, 345.62, 345.54, 342.96,  &
 344.38, 346.23, 349.17, 351.79, 354.71, 356.97, 358.29, 362.29,  &
 364.15, 364.97, 367.81, 368.98, 369.07, 372.17, 377.79, 381.25,  &
 384.22, 388.66, 393.58, 396.98, 398.72, 400.61, 404.06, 408.23,  &
 412.47, 415.58, 416.17, 416.53, 419.55, 425.88, 433.30, 437.73,  &
 438.13, 439.79, 441.51, 438.71, 434.25, 437.54, 448.95, 448.86,  &
 439.46, 437.10, 439.34, 444.33, 455.00, 467.05, 473.04, 469.64,  &
 467.53, 473.78, 477.50, 477.50, 480.96, 483.94, 482.19, 479.08,  &
 482.09, 493.43, 498.40, 492.05, 489.53, 493.34, 495.51, 496.52,  &
 499.57, 504.65, 509.68, 512.00, 512.05, 512.31, 515.00, 520.70,  &
 527.30, 531.88, 532.16, 530.48, 532.33, 539.26, 548.57, 553.00,  &
 548.96, 546.05, 551.00, 556.41, 557.21, 557.85, 560.95, 564.02,  &
 565.57, 566.38, 567.88, 571.48, 576.68, 581.54, 586.51, 593.62,  &
 600.70, 602.79, 601.39, 603.00, 606.88, 605.95, 600.97, 600.79,  &
 607.21, 612.87, 614.13, 614.39, 616.61, 620.53, 625.19, 629.78,  &
 633.79, 637.31, 640.47, 642.53, 642.62, 641.93, 643.11, 646.68 /

! SOLAR SPECTRUM FROM   9460 TO  12320 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=474,617) /                                     &
 650.57, 654.30, 660.95, 672.10, 682.31, 684.89, 682.20, 682.53,  &
 687.79, 691.42, 689.62, 688.14, 693.71, 703.25, 708.07, 706.22,  &
 704.64, 708.97, 717.35, 725.43, 731.08, 734.17, 735.41, 736.60,  &
 739.34, 742.90, 745.04, 744.29, 742.44, 749.53, 755.70, 758.82,  &
 766.31, 761.53, 762.09, 769.68, 764.18, 763.75, 768.88, 762.69,  &
 753.93, 762.38, 765.79, 772.19, 760.67, 762.10, 766.76, 766.98,  &
 769.35, 773.50, 766.84, 763.60, 773.82, 777.18, 779.61, 792.48,  &
 797.54, 787.81, 793.75, 805.96, 804.77, 806.62, 821.72, 830.28,  &
 827.54, 831.06, 830.20, 826.22, 823.28, 822.18, 833.92, 854.58,  &
 859.80, 862.56, 871.16, 875.16, 867.67, 863.87, 883.30, 893.40,  &
 897.74, 905.24, 905.38, 911.07, 930.21, 939.24, 934.74, 935.15,  &
 942.38, 948.13, 947.00, 951.88, 960.12, 951.88, 954.22, 959.07,  &
 963.36, 980.16, 983.66, 978.76, 979.38, 985.24, 977.08, 919.94,  &
 899.68, 962.91, 997.17, 999.93, 995.65, 999.93, 1014.9, 951.57,  &
 893.52, 955.14, 1003.1, 990.13, 978.79, 1011.2, 1034.7, 1031.9,  &
 1029.9, 1039.7, 1045.5, 1044.1, 1049.6, 1056.1, 1049.8, 1038.0,  &
 1051.9, 1072.2, 1075.5, 1077.0, 1079.3, 1078.0, 1075.7, 1079.7,  &
 1081.0, 1069.8, 1078.4, 1104.3, 1111.4, 1111.7, 1117.6, 1119.6 /

! SOLAR SPECTRUM FROM  12340 TO  15200 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=618,761) /                                     &
 1109.3, 1100.6, 1112.9, 1122.7, 1119.5, 1123.9, 1136.1, 1143.7,  &
 1140.5, 1141.2, 1151.5, 1148.7, 1138.3, 1141.0, 1150.6, 1160.1,  &
 1170.6, 1177.7, 1179.8, 1181.7, 1182.4, 1179.8, 1181.8, 1188.3,  &
 1190.0, 1191.4, 1197.0, 1196.0, 1192.2, 1200.6, 1210.4, 1209.1,  &
 1207.5, 1205.3, 1193.3, 1192.9, 1220.0, 1243.3, 1245.4, 1241.5,  &
 1240.2, 1241.1, 1244.0, 1248.5, 1253.2, 1257.1, 1259.9, 1261.9,  &
 1263.6, 1265.7, 1269.6, 1277.0, 1284.2, 1284.4, 1282.7, 1287.2,  &
 1286.8, 1272.3, 1262.2, 1270.7, 1288.8, 1304.8, 1311.8, 1312.2,  &
 1314.4, 1320.2, 1326.2, 1328.4, 1325.3, 1322.5, 1325.4, 1334.6,  &
 1346.4, 1354.0, 1353.7, 1347.3, 1338.3, 1331.0, 1329.7, 1338.0,  &
 1351.9, 1363.0, 1368.8, 1372.0, 1375.9, 1382.1, 1387.8, 1388.8,  &
 1388.2, 1392.2, 1401.7, 1412.9, 1418.2, 1410.7, 1395.9, 1385.7,  &
 1388.1, 1405.0, 1424.0, 1428.1, 1422.2, 1423.6, 1434.5, 1445.2,  &
 1450.7, 1451.8, 1451.5, 1453.9, 1459.9, 1466.9, 1471.3, 1469.4,  &
 1462.5, 1460.4, 1468.9, 1481.8, 1490.8, 1495.3, 1497.9, 1500.7,  &
 1505.2, 1510.0, 1512.3, 1512.7, 1515.6, 1521.6, 1524.2, 1520.7,  &
 1520.3, 1531.6, 1545.7, 1548.2, 1541.7, 1542.2, 1553.6, 1563.6,  &
 1563.6, 1559.9, 1561.3, 1569.9, 1581.6, 1577.6, 1529.7, 1447.0 /

! SOLAR SPECTRUM FROM  15220 TO  18080 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=762,905) /                                     &
 1396.9, 1428.7, 1506.4, 1567.1, 1594.0, 1606.1, 1613.5, 1609.0,  &
 1588.6, 1567.8, 1567.3, 1587.2, 1610.2, 1624.4, 1630.2, 1630.9,  &
 1628.1, 1622.3, 1616.9, 1618.9, 1631.6, 1648.1, 1658.2, 1659.7,  &
 1658.1, 1658.0, 1659.4, 1660.4, 1659.2, 1653.7, 1645.3, 1642.1,  &
 1652.7, 1674.2, 1694.1, 1700.6, 1703.4, 1697.6, 1654.5, 1644.4,  &
 1661.6, 1676.3, 1707.7, 1703.1, 1710.8, 1732.3, 1716.5, 1719.6,  &
 1729.6, 1683.1, 1628.5, 1683.5, 1727.0, 1707.8, 1689.4, 1698.4,  &
 1733.1, 1737.8, 1714.1, 1734.6, 1750.1, 1750.1, 1760.3, 1764.3,  &
 1765.3, 1769.4, 1779.9, 1793.0, 1765.1, 1729.4, 1745.9, 1753.4,  &
 1758.1, 1775.0, 1768.4, 1767.9, 1789.5, 1806.6, 1799.3, 1782.6,  &
 1779.3, 1792.1, 1809.7, 1808.0, 1794.4, 1818.6, 1774.2, 1648.5,  &
 1674.3, 1789.3, 1847.2, 1848.3, 1812.9, 1796.4, 1840.3, 1868.3,  &
 1864.6, 1873.2, 1872.2, 1856.0, 1845.0, 1842.4, 1823.9, 1795.1,  &
 1819.6, 1861.5, 1857.7, 1838.6, 1840.5, 1863.5, 1876.8, 1884.4,  &
 1894.9, 1875.2, 1821.2, 1779.4, 1810.2, 1855.3, 1831.8, 1837.3,  &
 1882.3, 1866.4, 1819.6, 1804.8, 1831.4, 1861.6, 1867.1, 1862.9,  &
 1851.9, 1834.7, 1835.2, 1845.1, 1831.9, 1803.6, 1792.5, 1821.8,  &
 1845.8, 1832.3, 1847.6, 1894.2, 1909.2, 1901.0, 1891.2, 1869.9 /

! SOLAR SPECTRUM FROM  18100 TO  20960 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=906,1049) /                                    &
 1854.4, 1865.8, 1873.7, 1868.8, 1881.7, 1897.1, 1884.2, 1856.2,  &
 1840.6, 1855.1, 1885.3, 1903.6, 1900.1, 1887.4, 1887.7, 1879.0,  &
 1844.5, 1844.1, 1877.1, 1847.3, 1785.1, 1792.6, 1848.7, 1894.4,  &
 1908.8, 1892.8, 1867.4, 1885.6, 1959.9, 1971.9, 1895.8, 1883.5,  &
 1917.6, 1853.8, 1793.0, 1875.6, 1974.0, 1975.7, 1943.9, 1926.4,  &
 1914.4, 1902.7, 1882.5, 1813.3, 1710.8, 1717.9, 1859.7, 1965.1,  &
 1970.1, 1941.4, 1902.5, 1852.0, 1836.3, 1879.3, 1901.6, 1862.9,  &
 1839.1, 1840.9, 1780.0, 1684.9, 1677.3, 1718.7, 1697.3, 1684.3,  &
 1784.5, 1898.0, 1910.3, 1877.2, 1866.6, 1862.6, 1860.3, 1899.7,  &
 1971.0, 1999.9, 1970.9, 1936.5, 1922.8, 1922.8, 1924.0, 1917.2,  &
 1912.0, 1926.2, 1959.7, 1995.4, 1995.9, 1938.8, 1883.5, 1894.7,  &
 1933.3, 1935.1, 1899.3, 1852.7, 1820.2, 1821.5, 1865.2, 1935.5,  &
 1966.1, 1919.6, 1881.2, 1931.5, 2015.6, 2050.0, 2021.4, 1960.8,  &
 1938.2, 1997.0, 2051.0, 2003.4, 1912.1, 1880.2, 1895.2, 1898.0,  &
 1898.8, 1938.3, 1994.2, 2010.0, 1982.4, 1948.8, 1927.3, 1911.6,  &
 1877.7, 1791.6, 1679.8, 1645.0, 1727.3, 1845.2, 1926.2, 1973.4,  &
 2005.2, 2021.6, 2021.8, 2025.7, 2054.3, 2086.5, 2082.6, 2052.9,  &
 2047.1, 2070.2, 2072.4, 2038.1, 2020.2, 2049.9, 2074.0, 2038.1 /

! SOLAR SPECTRUM FROM  20980 TO  23840 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=1050,1193) /                                   &
 1978.6, 1963.5, 1996.8, 2037.5, 2057.5, 2048.2, 2018.4, 1999.2,  &
 2011.4, 2039.5, 2056.0, 2040.2, 1981.8, 1911.4, 1891.8, 1938.3,  &
 1991.7, 2005.5, 2000.8, 2011.3, 2022.7, 1997.5, 1947.7, 1936.3,  &
 1986.6, 2037.9, 2032.8, 1995.7, 1984.0, 2012.0, 2055.5, 2091.6,  &
 2106.5, 2094.9, 2070.4, 2052.8, 2046.7, 2043.8, 2035.5, 2016.6,  &
 1988.4, 1973.3, 1999.0, 2057.4, 2103.8, 2109.4, 2089.4, 2068.5,  &
 2051.8, 2031.2, 2005.9, 1986.7, 1981.5, 1979.4, 1964.1, 1943.6,  &
 1951.8, 2007.3, 2083.2, 2139.1, 2158.0, 2143.3, 2103.2, 2050.9,  &
 2001.9, 1974.5, 1988.0, 2037.8, 2075.1, 2050.6, 1971.5, 1884.5,  &
 1828.5, 1820.9, 1866.4, 1935.3, 1974.2, 1958.7, 1925.1, 1920.2,  &
 1949.7, 1984.6, 1996.4, 1966.4, 1884.8, 1781.9, 1726.8, 1759.4,  &
 1817.4, 1800.4, 1692.6, 1593.2, 1598.6, 1700.3, 1823.8, 1909.7,  &
 1937.7, 1902.5, 1822.4, 1737.8, 1683.2, 1666.8, 1682.7, 1715.3,  &
 1734.1, 1712.4, 1668.2, 1655.0, 1698.1, 1727.2, 1636.9, 1415.7,  &
 1204.2, 1155.8, 1278.4, 1450.0, 1560.5, 1595.1, 1587.8, 1570.6,  &
 1565.8, 1590.3, 1640.5, 1688.4, 1708.1, 1703.6, 1700.7, 1718.5,  &
 1749.0, 1772.2, 1772.5, 1745.2, 1690.2, 1624.9, 1589.0, 1618.5,  &
 1701.3, 1783.2, 1816.4, 1800.7, 1765.0, 1734.1, 1714.6, 1705.0 /

! SOLAR SPECTRUM FROM  23860 TO  26720 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=1194,1337) /                                   &
 1701.6, 1696.6, 1682.0, 1661.4, 1657.2, 1693.0, 1763.2, 1826.5,  &
 1841.6, 1806.1, 1755.6, 1725.8, 1724.2, 1736.8, 1749.0, 1756.1,  &
 1759.5, 1762.1, 1770.2, 1791.7, 1826.8, 1848.9, 1819.6, 1720.7,  &
 1595.5, 1513.9, 1522.5, 1602.0, 1706.2, 1793.4, 1837.9, 1820.3,  &
 1738.3, 1631.1, 1553.1, 1539.2, 1574.3, 1623.9, 1660.6, 1676.8,  &
 1673.1, 1652.9, 1626.4, 1606.7, 1604.2, 1620.9, 1654.5, 1701.2,  &
 1752.2, 1796.2, 1822.8, 1827.4, 1808.5, 1767.0, 1713.9, 1667.3,  &
 1643.7, 1643.5, 1652.5, 1655.3, 1638.7, 1592.2, 1506.4, 1377.3,  &
 1209.5, 1010.5, 807.59, 666.84, 664.53, 835.23, 1099.6, 1330.7,  &
 1423.2, 1363.7, 1194.1, 961.77, 725.04, 551.29, 504.01, 596.30,  &
 775.15, 975.62, 1150.2, 1287.2, 1386.1, 1447.5, 1473.7, 1468.5,  &
 1435.2, 1376.9, 1296.0, 1195.5, 1085.3, 985.40, 917.25, 894.59,  &
 910.86, 951.53, 1001.7, 1046.4, 1070.7, 1061.2, 1021.2, 977.16,  &
 959.15, 982.06, 1020.5, 1032.6, 983.44, 879.83, 762.66, 675.28,  &
 643.33, 662.65, 721.49, 808.35, 913.24, 1027.0, 1139.9, 1236.2,  &
 1293.2, 1287.1, 1210.4, 1102.1, 1021.6, 1022.8, 1109.3, 1232.6,  &
 1337.0, 1383.1, 1372.8, 1324.7, 1257.7, 1188.8, 1133.5, 1106.5,  &
 1113.7, 1136.8, 1147.9, 1121.4, 1054.1, 968.10, 889.19, 837.87 /

! SOLAR SPECTRUM FROM  26740 TO  28780 CM-1,  IN STEPS OF 20 CM-
data (solara(i), i=1338,1440) /                                   &
 817.64, 823.72, 851.04, 896.53, 959.85, 1041.2, 1137.6, 1231.2,  &
 1294.4, 1299.9, 1241.2, 1155.0, 1092.0, 1097.1, 1170.2, 1263.5,  &
 1322.4, 1307.4, 1233.6, 1146.1, 1090.8, 1092.5, 1134.6, 1188.9,  &
 1228.9, 1245.5, 1248.5, 1250.3, 1260.5, 1274.6, 1279.5, 1261.8,  &
 1214.3, 1145.4, 1069.6, 1001.4, 952.52, 930.48, 941.68, 990.34,  &
 1064.4, 1135.2, 1171.5, 1149.1, 1076.3, 984.35, 906.25, 868.17,  &
 873.75, 915.33, 984.41, 1067.2, 1137.1, 1163.1, 1115.5, 990.55,  &
 830.93, 692.29, 627.44, 654.10, 739.24, 838.88, 911.69, 941.90,  &
 944.42, 939.58, 946.10, 970.23, 1005.2, 1042.4, 1073.8, 1097.0,  &
 1114.3, 1128.8, 1142.9, 1153.4, 1152.4, 1131.5, 1084.2, 1016.7,  &
 945.95, 890.37, 866.15, 876.54, 913.13, 966.10, 1025.4, 1080.2,  &
 1119.0, 1102.7, 1243.5, 1209.9, 1079.2, 852.20, 956.80, 842.31,  &
 897.44, 1081.8, 914.23, 993.09, 1049.8, 844.95, 839.16/

! SOLAR SPECTRUM FROM  28400 TO  29830 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1,144) /                                       &
 876.54, 892.17, 913.13, 938.18, 966.10, 995.62, 1025.4, 1054.1,  &
 1080.2, 1102.1, 1119.0, 1132.2, 1102.7, 1159.3, 1243.5, 1238.3,  &
 1209.9, 1196.2, 1079.2, 895.60, 852.20, 935.59, 956.80, 897.09,  &
 842.31, 821.15, 897.44, 1042.7, 1081.8, 988.79, 914.23, 929.38,  &
 993.09, 1041.9, 1049.8, 984.33, 844.95, 770.76, 839.16, 939.65,  &
 1026.1, 1121.1, 1162.6, 1142.6, 1077.9, 1027.3, 1078.2, 1094.3,  &
 969.83, 853.72, 849.91, 909.12, 995.68, 1095.0, 1146.9, 1086.3,  &
 1010.4, 1065.4, 1128.9, 1080.6, 987.93, 898.18, 835.20, 771.63,  &
 687.12, 614.52, 606.14, 737.09, 908.13, 997.64, 1080.6, 1126.3,  &
 1056.7, 1028.4, 1141.7, 1252.6, 1225.3, 1103.2, 1038.6, 1043.4,  &
 1002.9, 965.51, 1035.0, 1150.7, 1200.9, 1152.0, 1068.5, 995.84,  &
 889.52, 818.48, 907.01, 1042.2, 1055.6, 1000.6, 972.00, 985.72,  &
 1027.2, 1054.8, 1078.0, 1126.6, 1205.3, 1245.7, 1201.0, 1144.7,  &
 1097.5, 1030.1, 926.85, 836.71, 864.11, 993.50, 1075.3, 1032.6,  &
 1008.9, 1066.1, 1067.4, 1004.8, 971.54, 923.18, 815.71, 799.70,  &
 946.19, 1100.1, 1126.4, 1032.2, 895.14, 784.30, 734.77, 726.53,  &
 726.88, 765.54, 863.90, 992.24, 1070.9, 1028.1, 858.78, 647.15,  &
 563.18, 679.98, 906.40, 1094.3, 1155.3, 1124.3, 1098.4, 1109.5 /

! SOLAR SPECTRUM FROM  29840 TO  31270 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=145,288) /                                     &
 1076.2, 944.17, 849.20, 928.54, 1062.0, 1118.9, 1119.2, 1074.6,  &
 1005.8, 980.02, 999.11, 1002.4, 939.78, 838.12, 816.13, 908.73,  &
 1014.9, 1058.3, 1043.7, 987.54, 946.35, 981.40, 1055.8, 1094.3,  &
 1028.3, 916.41, 908.99, 991.83, 1049.6, 1076.2, 1093.5, 1076.3,  &
 1014.5, 949.61, 947.26, 1001.2, 1051.5, 1072.8, 1068.0, 1012.5,  &
 907.81, 866.30, 950.89, 1037.5, 1079.5, 1183.9, 1291.3, 1268.6,  &
 1199.3, 1188.6, 1188.0, 1186.6, 1198.2, 1171.3, 1132.6, 1131.6,  &
 1096.0, 971.10, 847.07, 836.62, 922.78, 990.99, 987.51, 969.24,  &
 981.46, 981.36, 971.95, 985.34, 1003.0, 1037.2, 1071.2, 1065.7,  &
 1026.7, 984.84, 1002.7, 1070.3, 1117.5, 1116.0, 1048.9, 965.34,  &
 972.27, 1045.7, 1096.6, 1127.5, 1133.5, 1099.6, 1079.3, 1082.9,  &
 1026.8, 927.50, 879.08, 858.83, 831.01, 807.82, 789.56, 813.75,  &
 893.46, 937.62, 901.56, 864.46, 873.35, 891.03, 862.46, 810.30,  &
 787.36, 752.93, 715.34, 708.07, 728.93, 786.79, 807.73, 736.28,  &
 645.08, 616.90, 649.17, 691.77, 749.18, 820.21, 820.68, 791.26,  &
 854.27, 940.56, 956.38, 909.42, 824.18, 767.17, 722.06, 653.42,  &
 624.67, 633.73, 655.14, 707.93, 784.94, 880.79, 961.15, 985.60,  &
 986.18, 966.53, 921.47, 888.89, 855.85, 851.66, 886.78, 850.97 /

! SOLAR SPECTRUM FROM  31280 TO  32710 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=289,432) /                                     &
 766.97, 738.95, 724.53, 657.61, 587.77, 616.86, 760.61, 903.23,  &
 917.27, 838.49, 784.80, 759.41, 719.61, 671.48, 624.63, 588.57,  &
 574.70, 596.68, 698.02, 866.39, 974.82, 960.37, 930.10, 962.65,  &
 1007.1, 1001.9, 926.29, 816.64, 763.25, 772.93, 762.66, 729.39,  &
 725.01, 727.16, 672.73, 581.42, 520.97, 488.80, 478.60, 542.08,  &
 663.71, 749.48, 785.87, 811.05, 818.19, 813.80, 824.54, 836.62,  &
 799.66, 728.00, 660.36, 559.28, 473.28, 550.16, 752.04, 885.84,  &
 906.80, 912.21, 929.32, 899.72, 830.20, 774.56, 736.42, 724.09,  &
 740.12, 754.11, 764.96, 780.76, 788.94, 784.87, 758.80, 725.91,  &
 751.84, 804.24, 777.73, 703.36, 665.27, 663.99, 679.36, 706.09,  &
 757.57, 836.09, 880.02, 881.18, 907.91, 929.26, 894.32, 874.01,  &
 918.56, 953.50, 922.32, 866.61, 836.54, 825.28, 752.54, 586.02,  &
 427.46, 374.05, 437.23, 534.32, 556.74, 563.11, 629.31, 631.26,  &
 518.76, 438.31, 460.31, 530.45, 608.50, 657.99, 662.08, 686.17,  &
 775.18, 843.11, 797.46, 685.33, 611.33, 628.74, 711.36, 754.94,  &
 728.80, 722.79, 726.38, 679.68, 665.83, 710.48, 723.10, 724.09,  &
 760.18, 784.01, 742.78, 634.33, 546.55, 563.54, 611.03, 623.16,  &
 665.36, 743.55, 764.46, 671.14, 513.18, 401.86, 405.77, 515.72 /

! SOLAR SPECTRUM FROM  32720 TO  34150 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=433,576) /                                     &
 639.90, 677.85, 679.55, 759.33, 848.11, 819.89, 751.75, 710.50,  &
 615.33, 525.09, 583.35, 715.23, 767.53, 739.10, 664.05, 580.57,  &
 572.85, 634.13, 648.77, 561.27, 497.72, 591.71, 737.83, 794.19,  &
 802.51, 799.33, 735.79, 658.41, 659.47, 718.18, 761.67, 697.24,  &
 545.14, 474.47, 526.96, 597.65, 584.74, 447.28, 291.35, 261.28,  &
 330.26, 401.96, 466.32, 531.26, 572.34, 584.86, 585.17, 569.46,  &
 558.27, 559.41, 512.02, 426.37, 378.14, 398.26, 473.49, 542.18,  &
 531.76, 437.48, 341.85, 305.82, 299.88, 328.12, 440.04, 586.46,  &
 660.32, 625.22, 510.26, 418.85, 447.36, 534.89, 605.86, 667.07,  &
 687.31, 636.79, 549.63, 472.88, 419.53, 370.06, 327.98, 320.49,  &
 354.00, 399.17, 450.98, 528.34, 608.25, 696.07, 774.28, 760.75,  &
 690.58, 648.20, 580.63, 477.96, 453.91, 488.74, 464.02, 421.59,  &
 444.32, 446.59, 375.95, 342.13, 397.49, 510.97, 646.38, 725.14,  &
 703.06, 639.06, 619.10, 654.66, 665.99, 611.40, 580.22, 607.29,  &
 591.05, 542.30, 583.82, 673.02, 673.21, 582.44, 465.73, 377.25,  &
 377.04, 487.27, 607.93, 617.52, 583.46, 601.68, 615.94, 575.47,  &
 541.63, 542.06, 522.28, 472.49, 423.29, 438.09, 556.72, 664.34,  &
 669.88, 657.45, 684.71, 705.70, 683.11, 600.81, 509.90, 497.64 /

! SOLAR SPECTRUM FROM  34160 TO  35590 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=577,720) /                                     &
 511.07, 496.07, 500.32, 518.70, 529.91, 563.00, 609.20, 626.49,  &
 622.11, 615.72, 600.44, 591.26, 598.12, 593.07, 590.94, 631.58,  &
 696.48, 718.48, 676.11, 631.56, 619.64, 620.53, 624.10, 636.56,  &
 658.02, 688.78, 724.81, 742.60, 722.31, 675.86, 665.96, 704.73,  &
 703.70, 645.00, 598.26, 587.77, 590.94, 575.93, 528.03, 477.92,  &
 457.52, 456.80, 454.91, 448.65, 445.47, 445.38, 444.43, 446.04,  &
 455.91, 468.02, 454.34, 393.32, 301.22, 211.44, 167.11, 193.99,  &
 254.01, 305.35, 353.03, 385.08, 387.03, 391.60, 406.20, 415.34,  &
 435.34, 469.77, 492.15, 472.73, 409.86, 353.25, 340.68, 355.27,  &
 379.77, 401.81, 409.67, 406.89, 393.16, 378.89, 375.20, 373.52,  &
 360.19, 322.79, 273.55, 237.76, 212.33, 184.80, 156.20, 127.75,  &
 96.269, 68.806, 62.047, 77.143, 100.47, 127.56, 159.88, 194.05,  &
 225.20, 254.64, 285.75, 300.14, 294.40, 308.92, 340.83, 346.26,  &
 336.29, 347.54, 373.81, 388.78, 372.68, 325.29, 294.40, 317.56,  &
 360.30, 378.08, 374.22, 374.03, 383.34, 387.88, 377.55, 356.96,  &
 340.67, 328.71, 314.00, 316.91, 344.51, 355.54, 335.66, 318.68,  &
 318.65, 322.43, 318.61, 304.92, 284.84, 268.13, 265.80, 273.55,  &
 274.18, 252.38, 215.04, 188.60, 181.31, 181.31, 180.78, 175.24 /

! SOLAR SPECTRUM FROM  35600 TO  37030 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=721,864) /                                     &
 162.06, 145.08, 128.76, 113.76, 98.078, 83.072, 76.222, 78.359,  &
 78.434, 74.235, 75.843, 80.321, 77.859, 70.298, 64.651, 67.049,  &
 77.810, 83.167, 75.286, 71.202, 80.549, 92.008, 100.17, 108.63,  &
 119.44, 130.78, 142.31, 158.94, 177.12, 186.40, 186.60, 181.47,  &
 175.30, 175.54, 179.00, 177.04, 172.60, 172.67, 178.98, 193.77,  &
 215.13, 233.62, 252.05, 277.68, 298.91, 298.40, 280.81, 274.21,  &
 286.52, 285.46, 259.71, 241.39, 246.98, 259.87, 274.27, 298.47,  &
 316.85, 303.19, 263.69, 229.31, 227.90, 256.12, 281.58, 300.19,  &
 310.56, 279.54, 211.93, 152.18, 129.94, 147.47, 181.62, 215.37,  &
 239.50, 233.12, 191.55, 139.41, 110.51, 118.93, 134.79, 129.05,  &
 124.39, 143.53, 158.29, 141.84, 116.32, 111.59, 128.93, 149.17,  &
 153.44, 145.63, 148.52, 159.25, 155.84, 154.17, 177.28, 203.40,  &
 207.35, 205.27, 222.85, 253.18, 271.28, 279.27, 302.17, 321.47,  &
 288.83, 230.14, 206.40, 213.22, 216.49, 207.46, 196.20, 195.21,  &
 202.03, 194.33, 164.86, 136.65, 123.87, 128.14, 161.89, 216.99,  &
 253.68, 249.26, 222.89, 213.11, 243.64, 293.10, 309.42, 286.40,  &
 269.61, 272.23, 271.67, 265.84, 265.61, 264.77, 266.03, 289.51,  &
 325.67, 337.34, 321.17, 300.30, 282.60, 287.14, 322.06, 335.79 /

! SOLAR SPECTRUM FROM  37040 TO  38470 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=865,1008) /                                    &
 297.22, 254.10, 243.47, 239.49, 219.32, 211.94, 239.28, 271.43,  &
 279.37, 272.26, 264.77, 250.52, 229.93, 222.15, 235.30, 256.79,  &
 275.28, 286.92, 284.85, 269.52, 255.05, 253.46, 263.22, 274.78,  &
 279.19, 270.17, 249.41, 229.04, 221.64, 231.38, 252.70, 280.64,  &
 310.06, 328.33, 325.01, 290.26, 238.97, 223.38, 257.24, 282.60,  &
 264.32, 243.34, 253.18, 272.89, 271.32, 256.12, 260.24, 271.35,  &
 257.11, 236.61, 238.72, 248.92, 255.90, 272.04, 291.78, 297.40,  &
 288.09, 283.28, 292.92, 301.74, 309.07, 322.05, 320.42, 295.43,  &
 269.65, 254.41, 240.88, 228.18, 221.23, 213.72, 201.23, 197.17,  &
 212.29, 233.39, 247.65, 261.74, 286.17, 322.49, 349.47, 338.28,  &
 297.06, 261.55, 252.28, 264.65, 286.92, 298.94, 280.45, 244.37,  &
 213.47, 193.03, 182.07, 168.54, 143.12, 114.10, 89.615, 73.589,  &
 73.990, 87.912, 96.265, 94.813, 96.604, 102.30, 102.15, 103.07,  &
 117.81, 137.41, 146.09, 144.28, 137.89, 128.11, 122.82, 128.19,  &
 130.66, 117.31, 98.912, 93.397, 105.63, 122.73, 126.39, 113.05,  &
 92.317, 76.340, 69.032, 66.324, 71.280, 87.431, 105.94, 114.02,  &
 107.91, 91.872, 75.208, 69.123, 75.930, 90.928, 109.71, 125.70,  &
 135.79, 141.14, 138.14, 121.33, 91.806, 63.497, 52.106, 59.555 /

! SOLAR SPECTRUM FROM  38480 TO  39910 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1009,1152) /                                   &
 81.015, 106.67, 118.97, 116.36, 110.82, 100.88, 89.056, 90.431,  &
 104.41, 114.95, 124.85, 148.87, 171.72, 167.22, 142.25, 118.42,  &
 98.653, 78.908, 68.133, 77.286, 100.93, 120.08, 125.49, 131.79,  &
 155.69, 180.75, 181.81, 166.77, 150.06, 133.24, 116.14, 97.728,  &
 81.629, 76.695, 87.607, 110.23, 134.88, 149.13, 147.64, 139.88,  &
 135.19, 135.07, 138.00, 136.73, 128.84, 122.22, 120.48, 121.98,  &
 123.08, 116.30, 101.43, 86.303, 74.719, 68.800, 71.327, 80.626,  &
 90.485, 96.739, 100.69, 100.81, 93.677, 84.740, 81.532, 82.893,  &
 84.564, 87.584, 91.780, 91.272, 87.014, 87.386, 90.149, 84.917,  &
 71.266, 57.873, 51.863, 53.876, 57.909, 58.508, 57.020, 57.432,  &
 60.671, 64.667, 67.362, 67.511, 64.233, 59.035, 55.697, 56.636,  &
 59.400, 59.070, 56.522, 55.834, 55.860, 54.039, 51.976, 52.344,  &
 54.667, 56.450, 56.751, 56.769, 58.002, 60.029, 59.602, 53.134,  &
 42.926, 35.588, 33.447, 35.171, 39.379, 44.371, 47.745, 46.933,  &
 42.441, 37.879, 35.595, 36.458, 41.048, 47.300, 51.098, 50.024,  &
 45.331, 41.282, 40.082, 40.000, 39.104, 37.329, 36.632, 37.792,  &
 39.189, 41.058, 45.214, 50.737, 54.281, 55.015, 56.138, 60.931,  &
 67.383, 69.534, 65.159, 56.372, 47.326, 44.322, 49.944, 59.696 /

! SOLAR SPECTRUM FROM  39920 TO  41350 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1153,1296) /                                   &
 67.929, 71.334, 69.905, 65.620, 59.303, 54.016, 55.880, 65.155,  &
 74.065, 76.217, 73.506, 71.406, 70.849, 69.749, 69.268, 71.380,  &
 72.721, 68.929, 61.665, 54.896, 47.420, 38.325, 32.219, 31.243,  &
 33.310, 35.358, 35.623, 36.840, 41.551, 47.499, 51.176, 50.344,  &
 45.362, 38.341, 33.130, 33.801, 40.140, 49.121, 55.385, 55.174,  &
 50.450, 46.511, 47.495, 51.883, 56.354, 59.603, 61.584, 63.215,  &
 64.603, 64.101, 59.027, 50.956, 47.633, 52.543, 58.883, 59.829,  &
 57.617, 56.727, 57.371, 57.898, 57.177, 55.129, 52.952, 52.018,  &
 52.186, 52.044, 50.269, 46.592, 42.515, 40.755, 41.887, 44.119,  &
 46.536, 48.858, 50.490, 51.919, 54.085, 54.707, 51.927, 49.449,  &
 49.865, 50.933, 50.496, 48.616, 46.717, 46.070, 46.263, 46.733,  &
 48.009, 50.187, 52.420, 53.536, 52.507, 51.380, 53.214, 56.985,  &
 60.614, 63.139, 63.999, 63.869, 65.100, 69.385, 74.743, 78.184,  &
 78.103, 74.113, 67.371, 60.849, 58.924, 62.682, 68.032, 69.117,  &
 64.604, 59.110, 55.998, 56.838, 61.778, 65.874, 65.079, 63.038,  &
 64.809, 69.911, 74.841, 76.439, 73.587, 68.853, 67.497, 72.675,  &
 80.602, 83.422, 78.957, 72.228, 66.737, 62.842, 61.535, 63.574,  &
 69.248, 76.577, 79.922, 77.755, 73.938, 70.518, 68.003, 66.339 /

! SOLAR SPECTRUM FROM  41360 TO  42790 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1297,1440) /                                   &
 63.979, 61.098, 59.421, 58.103, 55.741, 52.549, 48.079, 42.578,  &
 38.373, 37.297, 37.455, 34.861, 30.483, 29.634, 34.734, 42.460,  &
 47.066, 45.848, 40.157, 34.290, 31.584, 30.650, 29.054, 27.788,  &
 30.427, 37.570, 44.196, 46.880, 47.848, 49.166, 49.180, 45.002,  &
 38.135, 35.055, 38.095, 41.750, 40.899, 35.722, 28.884, 24.835,  &
 28.670, 39.646, 50.310, 55.725, 57.401, 58.110, 59.406, 59.360,  &
 53.420, 43.004, 34.787, 33.697, 39.682, 47.554, 52.605, 53.632,  &
 51.001, 45.266, 37.844, 31.030, 25.936, 22.799, 21.882, 23.484,  &
 27.857, 33.447, 37.319, 39.195, 42.826, 50.398, 58.752, 63.301,  &
 61.094, 53.532, 46.046, 41.118, 37.646, 36.304, 40.426, 50.893,  &
 61.553, 65.395, 62.680, 58.087, 54.622, 51.330, 46.874, 42.870,  &
 40.547, 39.760, 40.217, 40.359, 39.559, 40.667, 46.260, 53.413,  &
 56.041, 52.566, 46.674, 41.073, 35.511, 31.231, 31.082, 35.955,  &
 45.199, 55.464, 61.802, 63.505, 61.850, 56.412, 49.388, 46.369,  &
 50.058, 56.694, 60.884, 61.030, 58.107, 54.303, 51.940, 50.508,  &
 46.749, 39.155, 31.535, 28.959, 30.973, 32.670, 31.567, 29.340,  &
 27.275, 25.184, 24.264, 27.068, 34.296, 42.475, 47.230, 47.425,  &
 44.435, 40.538, 36.868, 33.020, 29.405, 28.753, 34.079, 44.246 /

! SOLAR SPECTRUM FROM  42800 TO  44230 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1441,1584) /                                   &
 53.780, 57.974, 56.376, 51.200, 45.308, 40.273, 35.900, 33.344,  &
 34.011, 36.858, 41.283, 47.374, 53.088, 56.201, 55.633, 50.843,  &
 43.997, 38.767, 36.248, 36.380, 40.762, 50.700, 63.371, 73.432,  &
 76.418, 70.373, 58.741, 47.034, 38.598, 34.664, 35.794, 42.084,  &
 49.973, 54.338, 53.956, 52.287, 52.778, 55.571, 59.034, 60.268,  &
 56.247, 47.362, 38.056, 32.889, 31.739, 31.734, 32.476, 35.060,  &
 39.091, 43.398, 48.131, 53.574, 58.749, 63.599, 68.971, 73.421,  &
 73.861, 69.003, 60.557, 51.865, 44.879, 42.060, 44.802, 47.950,  &
 46.882, 42.973, 39.293, 37.711, 37.137, 35.222, 32.243, 30.488,  &
 32.605, 40.429, 51.099, 57.710, 57.150, 52.992, 50.275, 49.986,  &
 49.778, 48.371, 46.421, 44.604, 42.730, 41.244, 41.565, 43.805,  &
 47.013, 48.992, 46.428, 40.595, 37.840, 42.353, 52.248, 60.529,  &
 61.566, 56.800, 52.041, 52.260, 57.077, 61.019, 60.712, 57.048,  &
 51.481, 46.352, 44.366, 44.947, 45.478, 44.944, 43.825, 42.105,  &
 39.466, 36.826, 35.907, 36.357, 35.661, 33.947, 33.690, 34.429,  &
 34.000, 32.645, 31.410, 30.281, 29.409, 29.127, 29.326, 29.869,  &
 30.601, 31.311, 32.099, 32.779, 32.757, 32.098, 31.975, 33.484,  &
 36.048, 39.169, 43.365, 47.244, 48.214, 45.786, 41.586, 38.775 /

! SOLAR SPECTRUM FROM  44240 TO  45670 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1585,1728) /                                   &
 40.753, 46.752, 51.684, 52.597, 51.449, 50.684, 49.450, 46.747,  &
 45.369, 47.685, 50.240, 48.961, 46.693, 48.600, 53.694, 56.465,  &
 54.341, 50.722, 49.877, 51.246, 52.088, 52.765, 56.254, 63.326,  &
 69.744, 71.066, 68.349, 65.123, 62.551, 59.195, 53.705, 48.161,  &
 46.236, 47.710, 49.660, 50.799, 51.836, 54.537, 59.647, 64.707,  &
 65.844, 61.634, 55.570, 54.083, 58.781, 64.888, 69.777, 74.008,  &
 76.492, 76.226, 74.746, 74.941, 77.801, 79.619, 76.190, 67.190,  &
 55.231, 45.813, 43.141, 45.647, 49.466, 52.231, 52.221, 48.886,  &
 44.716, 42.613, 43.385, 45.968, 48.121, 48.998, 49.885, 50.707,  &
 49.893, 48.319, 48.198, 50.280, 53.830, 55.914, 54.822, 52.939,  &
 51.944, 49.438, 42.956, 34.614, 28.100, 24.503, 24.203, 27.839,  &
 34.604, 41.615, 45.324, 45.444, 45.527, 47.179, 45.756, 36.862,  &
 26.037, 20.569, 20.329, 24.263, 30.863, 35.939, 36.711, 35.693,  &
 37.256, 40.862, 44.416, 48.800, 54.182, 57.655, 58.427, 59.965,  &
 63.940, 66.820, 65.465, 59.482, 49.396, 39.422, 34.182, 35.388,  &
 42.875, 52.034, 57.595, 59.093, 57.272, 52.172, 45.493, 39.419,  &
 35.581, 35.902, 40.354, 46.732, 53.309, 58.781, 61.785, 59.255,  &
 50.030, 41.567, 40.523, 43.584, 44.875, 42.754, 40.077, 39.941 /

! SOLAR SPECTRUM FROM  45680 TO  47110 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1729,1872) /                                   &
 40.977, 39.567, 34.955, 30.424, 31.039, 38.687, 47.480, 49.830,  &
 46.790, 44.829, 46.546, 50.415, 54.602, 57.656, 58.463, 57.276,  &
 55.621, 54.514, 53.338, 50.026, 42.817, 33.636, 27.134, 25.516,  &
 27.897, 31.392, 32.125, 29.463, 26.581, 25.956, 27.737, 31.175,  &
 34.959, 37.671, 38.641, 37.958, 36.733, 35.681, 33.877, 30.849,  &
 28.059, 27.615, 29.319, 29.375, 25.390, 20.659, 19.484, 22.297,  &
 27.282, 32.467, 35.906, 37.137, 37.895, 39.130, 39.777, 39.872,  &
 40.778, 42.317, 42.934, 40.430, 34.227, 27.701, 23.880, 22.174,  &
 21.639, 22.589, 25.184, 29.017, 32.981, 36.110, 38.580, 41.239,  &
 44.426, 46.939, 47.010, 44.165, 39.659, 35.556, 32.838, 31.546,  &
 32.676, 36.963, 42.333, 44.931, 43.704, 40.943, 37.973, 35.199,  &
 33.574, 33.339, 34.185, 36.347, 39.963, 43.964, 47.162, 48.987,  &
 48.976, 47.948, 48.004, 49.892, 51.065, 47.834, 40.489, 32.665,  &
 26.795, 24.461, 26.655, 31.928, 37.634, 41.345, 40.956, 36.827,  &
 32.110, 28.612, 26.482, 26.602, 28.831, 30.877, 30.976, 30.063,  &
 29.887, 30.305, 29.974, 28.265, 26.517, 27.066, 30.403, 34.539,  &
 37.104, 37.598, 37.252, 37.060, 36.498, 34.167, 29.814, 24.192,  &
 18.515, 15.086, 15.040, 17.158, 20.807, 25.682, 30.352, 34.203 /

! SOLAR SPECTRUM FROM  47120 TO  48550 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=1873,2016) /                                   &
 37.902, 42.531, 47.832, 50.509, 48.019, 42.616, 38.321, 37.370,  &
 40.172, 44.395, 46.132, 43.911, 38.396, 31.379, 26.275, 25.075,  &
 26.652, 28.963, 31.168, 34.168, 38.050, 40.231, 38.347, 32.741,  &
 26.199, 21.863, 20.249, 20.185, 21.726, 25.562, 30.318, 33.431,  &
 34.453, 34.959, 36.374, 37.870, 36.655, 31.966, 25.920, 21.264,  &
 20.663, 24.658, 30.263, 34.021, 34.336, 31.356, 26.926, 23.109,  &
 20.867, 20.684, 22.416, 24.878, 26.779, 27.334, 26.537, 25.210,  &
 24.013, 22.944, 21.800, 20.449, 19.290, 19.528, 21.742, 24.125,  &
 23.994, 21.559, 19.555, 18.915, 18.342, 17.335, 16.549, 16.479,  &
 17.211, 18.445, 19.294, 18.980, 17.912, 17.156, 17.103, 17.256,  &
 16.925, 15.842, 14.485, 13.683, 13.647, 13.914, 14.009, 13.770,  &
 13.456, 13.399, 13.547, 13.760, 14.060, 14.427, 14.644, 14.438,  &
 13.986, 13.749, 13.927, 14.390, 14.759, 14.822, 14.679, 14.448,  &
 14.186, 13.937, 13.754, 13.657, 13.540, 13.308, 13.053, 12.841,  &
 12.704, 12.742, 12.811, 12.662, 12.355, 12.100, 12.003, 12.014,  &
 12.067, 12.223, 12.444, 12.472, 12.164, 11.732, 11.515, 11.619,  &
 11.873, 12.028, 11.947, 11.722, 11.399, 10.930, 10.473, 10.205,  &
 10.224, 10.694, 11.468, 12.007, 12.083, 11.905, 11.498, 10.891 /

! SOLAR SPECTRUM FROM  48560 TO  49990 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2017,2160) /                                   &
 10.575, 10.846, 11.353, 11.612, 11.411, 10.876, 10.383, 10.305,  &
 10.695, 11.245, 11.636, 11.828, 11.918, 11.865, 11.674, 11.510,  &
 11.407, 11.303, 11.216, 11.143, 11.039, 10.983, 11.004, 10.900,  &
 10.653, 10.562, 10.781, 11.186, 11.605, 11.806, 11.582, 11.056,  &
 10.567, 10.335, 10.408, 10.729, 11.165, 11.540, 11.646, 11.372,  &
 10.933, 10.524, 9.9973, 9.3783, 8.9883, 9.0163, 9.4125, 9.9179,  &
 10.278, 10.472, 10.553, 10.575, 10.519, 10.216, 9.6821, 9.1499,  &
 8.7057, 8.3894, 8.3442, 8.6241, 9.1371, 9.7184, 10.191, 10.443,  &
 10.458, 10.289, 9.9772, 9.5829, 9.3097, 9.3195, 9.4694, 9.5182,  &
 9.4326, 9.2478, 8.8197, 7.9809, 6.9996, 6.4856, 6.7462, 7.5406,  &
 8.2813, 8.7258, 9.0682, 9.1665, 8.8637, 8.4638, 8.2393, 8.1656,  &
 8.1880, 8.3578, 8.6488, 8.8980, 9.0117, 9.0659, 9.1955, 9.4207,  &
 9.5526, 9.4237, 9.1290, 8.8441, 8.6138, 8.4237, 8.2979, 8.2598,  &
 8.2859, 8.3475, 8.4533, 8.6285, 8.8310, 8.8866, 8.6750, 8.3312,  &
 8.0091, 7.7296, 7.6239, 7.8692, 8.2725, 8.4086, 8.2515, 8.0914,  &
 8.0003, 7.9367, 7.9266, 7.9580, 8.0492, 8.2376, 8.4263, 8.4811,  &
 8.3309, 8.0263, 7.7632, 7.6987, 7.8124, 7.9390, 8.0183, 8.0816,  &
 8.0428, 7.8923, 7.6963, 7.4969, 7.4013, 7.4289, 7.4489, 7.4059 /

! SOLAR SPECTRUM FROM  50000 TO  51430 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2161,2304) /                                   &
 7.4198, 7.5261, 7.5252, 7.3239, 7.1263, 7.1423, 7.3340, 7.5049,  &
 7.5484, 7.5319, 7.5163, 7.4995, 7.5728, 7.8104, 8.0588, 8.0948,  &
 7.9140, 7.6978, 7.5116, 7.2138, 6.8063, 6.5430, 6.5232, 6.5869,  &
 6.5610, 6.3984, 6.1889, 6.0587, 6.0676, 6.1988, 6.3140, 6.2527,  &
 6.0929, 6.0277, 6.0941, 6.3031, 6.6594, 6.9398, 6.9566, 6.8310,  &
 6.7374, 6.6812, 6.6558, 6.8336, 7.2020, 7.4012, 7.2950, 7.0488,  &
 6.7966, 6.6293, 6.5868, 6.5980, 6.6007, 6.6501, 6.7627, 6.7853,  &
 6.6321, 6.4856, 6.5198, 6.6486, 6.7271, 6.7227, 6.6696, 6.6189,  &
 6.5979, 6.6188, 6.7110, 6.8343, 6.8750, 6.8250, 6.7885, 6.8266,  &
 6.8556, 6.8068, 6.8377, 7.0467, 7.2779, 7.4139, 7.4712, 7.4621,  &
 7.4071, 7.3592, 7.3372, 7.3220, 7.2938, 7.2531, 7.2052, 7.1335,  &
 7.0298, 6.8533, 6.5535, 6.2227, 6.0139, 5.9384, 5.9038, 5.8568,  &
 5.7909, 5.7326, 5.7745, 5.9608, 6.1865, 6.3681, 6.4997, 6.5437,  &
 6.4637, 6.2708, 6.0451, 5.9557, 6.0855, 6.2542, 6.2454, 6.0795,  &
 5.9102, 5.8447, 5.9218, 6.1063, 6.2895, 6.3271, 6.1097, 5.7421,  &
 5.4452, 5.2981, 5.3256, 5.4935, 5.6819, 5.8245, 5.8933, 5.9630,  &
 6.1703, 6.4525, 6.6325, 6.6965, 6.7185, 6.6238, 6.3107, 5.9241,  &
 5.6987, 5.6651, 5.7428, 5.8790, 5.9715, 5.9618, 5.9674, 6.0754 /

! SOLAR SPECTRUM FROM  51440 TO  52870 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2305,2448) /                                   &
 6.2541, 6.4300, 6.4968, 6.4564, 6.4082, 6.3024, 6.0135, 5.6431,  &
 5.3963, 5.2989, 5.2635, 5.2227, 5.1279, 4.9315, 4.6348, 4.3168,  &
 4.0151, 3.6625, 3.2906, 3.1028, 3.1349, 3.1994, 3.2596, 3.4144,  &
 3.5949, 3.6534, 3.6296, 3.6281, 3.5876, 3.4292, 3.2659, 3.2284,  &
 3.2576, 3.3002, 3.4535, 3.7372, 4.0573, 4.3558, 4.5999, 4.7781,  &
 4.8855, 4.8999, 4.8392, 4.7624, 4.7059, 4.6981, 4.7666, 4.8453,  &
 4.8236, 4.7293, 4.6861, 4.7132, 4.7725, 4.8713, 4.9596, 4.9527,  &
 4.8957, 4.9252, 5.0736, 5.2229, 5.2505, 5.1537, 5.0156, 4.8880,  &
 4.7686, 4.6549, 4.5534, 4.4828, 4.4661, 4.5040, 4.5905, 4.7033,  &
 4.7858, 4.8334, 4.9283, 5.0377, 5.0065, 4.8471, 4.6828, 4.5586,  &
 4.4812, 4.4314, 4.3903, 4.3830, 4.4066, 4.3900, 4.2973, 4.1978,  &
 4.1462, 4.1084, 4.1495, 4.3897, 4.6859, 4.8206, 4.7938, 4.6781,  &
 4.5222, 4.3959, 4.3358, 4.2947, 4.2259, 4.1452, 4.1060, 4.1462,  &
 4.2149, 4.2549, 4.3061, 4.3742, 4.3738, 4.2718, 4.1389, 4.0405,  &
 3.9457, 3.8127, 3.7099, 3.7344, 3.8589, 3.9598, 3.9525, 3.8377,  &
 3.6708, 3.5357, 3.4929, 3.5375, 3.6381, 3.7890, 3.9671, 4.0995,  &
 4.1421, 4.1302, 4.1235, 4.1623, 4.2506, 4.2948, 4.2231, 4.0993,  &
 3.9680, 3.9475, 4.1958, 4.5131, 4.6101, 4.5130, 4.3474, 4.1749 /

! SOLAR SPECTRUM FROM  52880 TO  54310 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2449,2592) /                                   &
 4.0467, 3.9956, 4.0078, 4.0374, 4.0255, 3.9379, 3.8192, 3.7529,  &
 3.7675, 3.8260, 3.8654, 3.8518, 3.8148, 3.8028, 3.8098, 3.7934,  &
 3.7660, 3.7944, 3.8689, 3.8978, 3.8856, 3.8923, 3.8570, 3.6940,  &
 3.4693, 3.3222, 3.2824, 3.2887, 3.3039, 3.3222, 3.3313, 3.3326,  &
 3.3482, 3.3807, 3.4188, 3.4602, 3.4972, 3.5151, 3.5155, 3.5165,  &
 3.5258, 3.5406, 3.5478, 3.5345, 3.5339, 3.5820, 3.6396, 3.6448,  &
 3.5872, 3.5112, 3.4804, 3.5257, 3.6238, 3.7290, 3.8023, 3.8024,  &
 3.7268, 3.6578, 3.6439, 3.6422, 3.6373, 3.6397, 3.6410, 3.6494,  &
 3.6608, 3.6251, 3.5212, 3.4020, 3.2845, 3.1230, 2.9483, 2.8515,  &
 2.8432, 2.8638, 2.8967, 2.9505, 3.0025, 3.0552, 3.1106, 3.1178,  &
 3.0596, 2.9854, 2.9316, 2.8903, 2.8590, 2.8500, 2.8450, 2.8121,  &
 2.7626, 2.7424, 2.7667, 2.8024, 2.8165, 2.8111, 2.8128, 2.8569,  &
 2.9659, 3.1062, 3.1990, 3.2128, 3.2088, 3.2391, 3.2661, 3.2364,  &
 3.1173, 2.9094, 2.6952, 2.5324, 2.3959, 2.2953, 2.2510, 2.2245,  &
 2.1811, 2.1301, 2.1482, 2.3257, 2.5856, 2.7226, 2.6495, 2.4508,  &
 2.2444, 2.0850, 1.9891, 1.9843, 2.0816, 2.2233, 2.3248, 2.3551,  &
 2.3479, 2.3606, 2.4296, 2.5361, 2.6128, 2.6216, 2.6069, 2.6196,  &
 2.6464, 2.6427, 2.5823, 2.4682, 2.3320, 2.2405, 2.2637, 2.3973 /

! SOLAR SPECTRUM FROM  54320 TO  55750 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2593,2736) /                                   &
 2.5524, 2.6891, 2.8508, 3.0103, 3.0681, 3.0064, 2.9114, 2.8609,  &
 2.8517, 2.8374, 2.7894, 2.7288, 2.7138, 2.7729, 2.8707, 2.9536,  &
 2.9953, 2.9911, 2.9398, 2.8550, 2.7732, 2.7303, 2.7366, 2.7650,  &
 2.7705, 2.7374, 2.6830, 2.6218, 2.5663, 2.5341, 2.5351, 2.5681,  &
 2.6124, 2.6305, 2.6024, 2.5431, 2.4840, 2.4546, 2.4684, 2.5100,  &
 2.5445, 2.5532, 2.5564, 2.5889, 2.6616, 2.7553, 2.8466, 2.9290,  &
 2.9958, 3.0175, 2.9774, 2.8990, 2.8001, 2.6927, 2.6171, 2.5931,  &
 2.5809, 2.5276, 2.4284, 2.3365, 2.3162, 2.3855, 2.4872, 2.5455,  &
 2.5773, 2.6809, 2.9720, 3.5757, 4.4006, 5.0044, 5.0295, 4.5135,  &
 3.7071, 2.9059, 2.3600, 2.1418, 2.1119, 2.0871, 2.0301, 2.0043,  &
 2.0361, 2.0963, 2.1520, 2.1878, 2.1955, 2.1864, 2.1899, 2.2170,  &
 2.2574, 2.2895, 2.2783, 2.2148, 2.1641, 2.2343, 2.4726, 2.8119,  &
 3.1288, 3.2984, 3.2206, 2.8859, 2.4473, 2.1436, 2.0729, 2.1391,  &
 2.2171, 2.2580, 2.2654, 2.2481, 2.2103, 2.1657, 2.1356, 2.1321,  &
 2.1438, 2.1461, 2.1396, 2.1460, 2.1588, 2.1581, 2.1481, 2.1343,  &
 2.1101, 2.0754, 2.0400, 2.0121, 1.9930, 1.9799, 1.9699, 1.9613,  &
 1.9537, 1.9454, 1.9312, 1.9058, 1.8726, 1.8470, 1.8465, 1.8693,  &
 1.8844, 1.8635, 1.8143, 1.7618, 1.7188, 1.6853, 1.6656, 1.6708 /

! SOLAR SPECTRUM FROM  55760 TO  57190 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2737,2880) /                                   &
 1.7036, 1.7519, 1.8120, 1.9015, 2.0124, 2.0980, 2.1385, 2.1481,  &
 2.1347, 2.1086, 2.0953, 2.1062, 2.1095, 2.0685, 2.0001, 1.9461,  &
 1.9194, 1.9088, 1.9023, 1.8977, 1.9049, 1.9300, 1.9588, 1.9635,  &
 1.9357, 1.9019, 1.8887, 1.8939, 1.9018, 1.9038, 1.8975, 1.8747,  &
 1.8289, 1.7716, 1.7303, 1.7330, 1.7900, 1.8782, 1.9548, 1.9907,  &
 1.9807, 1.9430, 1.9173, 1.9218, 1.9203, 1.8717, 1.7832, 1.6965,  &
 1.6389, 1.6077, 1.5924, 1.5818, 1.5583, 1.5142, 1.4616, 1.4237,  &
 1.4252, 1.4834, 1.5970, 1.7410, 1.8771, 1.9784, 2.0451, 2.0872,  &
 2.0909, 2.0384, 1.9573, 1.9002, 1.8824, 1.8663, 1.8193, 1.7540,  &
 1.6874, 1.6222, 1.5726, 1.5450, 1.5290, 1.5312, 1.5699, 1.6411,  &
 1.7186, 1.7678, 1.7546, 1.6623, 1.5115, 1.3588, 1.2605, 1.2348,  &
 1.2611, 1.3091, 1.3588, 1.3884, 1.3800, 1.3482, 1.3224, 1.3159,  &
 1.3437, 1.4142, 1.4950, 1.5443, 1.5521, 1.5282, 1.4902, 1.4606,  &
 1.4465, 1.4398, 1.4399, 1.4544, 1.4760, 1.4781, 1.4506, 1.4229,  &
 1.4185, 1.4221, 1.4119, 1.3908, 1.3779, 1.3813, 1.3933, 1.4087,  &
 1.4268, 1.4417, 1.4408, 1.4188, 1.3861, 1.3548, 1.3261, 1.2980,  &
 1.2769, 1.2731, 1.2856, 1.3002, 1.3056, 1.2987, 1.2817, 1.2590,  &
 1.2291, 1.1868, 1.1428, 1.1183, 1.1141, 1.1120, 1.1009, 1.0797 /

! SOLAR SPECTRUM FROM  57200 TO  57490 CM-1,  IN STEPS OF 10 CM-
data (solarb(i), i=2881,2910) /                                   &
 1.0523, 1.0284, 1.0251, 1.0577, 1.1195, 1.1791, 1.2061, 1.2013,  &
 1.1936, 1.2000, 1.2040, 1.1824, 1.1489, 1.1400, 1.1539, 1.1629,  &
 1.1617, 1.1586, 1.1564, 1.1572, 1.1565, 1.1399, 1.1037, 1.0627,  &
 1.0341, 1.0223, 1.0199, 1.0188, 1.0174, 1.0163  /

!       WM, W0, W1, W2  ARE STATEMENT FUNC USED BY
!       THE 4 POINT LAGRANGE INTERPOLATION

wm(p) = p*(p - 1)*(p - 2)
w0(p) = 3*(p**2 - 1)*(p - 2)
w1(p) = 3*p*(p + 1)*(p - 2)
w2(p) = p*(p**2 - 1)

sunirr = 0.0 !Variable initialized

if(v .lt. 0.0) then
   !         IF  V  IS TOO SMALL,  WRITE WARNING  +  RETURN SUNIRR = 0
   sunirr = 0.0
   write(*, 900) v

elseif( v .ge. 0.0  .and.  v .lt. 100.0 ) then
   !         FOR LOW FREQUENCIES USE A POWER LAW APPROXIMATION
  sunirr = a*v**b

elseif( v .ge. 100.0  .and.  v .lt. 28420.0 ) then
   !         USE  4 POINT INTERPOLATION  ON  ARRAY  SOLARA
   !         WHICH IS AT  20 CM-1  SPACING  FROM 0 TO 28720 CM-1
   i = 1 + int(v/20.0)
   p = mod(v, 20.0)/20.0
   sunirr = ( w2(p)*solara(i+2) - w1(p)*solara(i+1) +  &
           w0(p)*solara(i) - wm(p)*solara(i-1) ) / 6

elseif( v .ge. 28420.0  .and.  v .le. 57470.0 ) then
   !         USE  4 POINT INTERPOLATION  ON  ARRAY  SOLARB
   !         WHICH IS AT  10 CM-1  SPACING  FROM 28400 TO 57490 CM-1
   i = int(v/10.0) - 2839
   p = mod(v, 10.0)/10.0
   sunirr = ( w2(p)*solarb(i+2) - w1(p)*solarb(i+1) +  &
           w0(p)*solarb(i) - wm(p)*solarb(i-1) ) / 6

elseif( v .gt. 57470.0 ) then
   !         IF  V  IS TOO LARGE,  WRITE WARNING  +  RETURN SUNIRR = 0
   sunirr = 0.0
   write(*, 900) v

endif

900 format('0 *****  WARNING - INPUT FREQUENCY = ', 1PG12.5, 'CM-1',  &
  /, '   OUTSIDE VALID RANGE OF 0 TO 57470 CM-1    *******', / )

return
END FUNCTION sunirr
