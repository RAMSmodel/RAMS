!##############################################################################
Subroutine geodat (n2,n3,datr,hfn,ofn,vt2da,vt2db,ngr,vnam)

use mem_grid
use io_params
use node_mod

implicit none

integer :: n2,n3,ngr
real ::  vt2da(*),vt2db(*),datr(n2,n3)
character(len=strl1) :: hfn,ofn,title,h5name
character(len=3) :: vnam

integer :: lb,iblksizo,no,isbego,iwbego,iodim,mof,niq,njq,np
real :: offlat,offlon,deltallo,deltaxq,deltaxp

real,allocatable:: dato(:)

print*,'====================================================='
if(VNAM(1:2).eq.'TO') then
   print*,'starting topography on grid:',NGR
elseif(VNAM(1:2).eq.'ZO') then
   print*,'starting surface roughness on grid:',NGR
else
   print*,'starting '//vnam//' data on grid:',NGR
endif

LB=len_trim(HFN)
if(LB.le.0) then
   print*,'==================================================='
   print*,'|  Problem in GEODAT, Input data prefix incorrect !'
   print*,'|  Grid :',ngrid
   print*,'|  File prefix:',HFN
   print*,'==================================================='
   stop 'GEODAT-file'
endif

if((vnam(1:2).eq.'TO'.or.vnam(1:2).eq.'ZO').and.  &
   (ITOPSFLG(NGR).eq.1.and.TOPTENH(NGR).gt.1.)) then
   print*,'==================================================='
   print*,'|  Problem in GEODAT, Silhouette weight too high !'
   print*,'|  Grid :',NGR
   print*,'|  Weight (range TOPTENH=0-1):',TOPTENH(NGR)
   print*,'==================================================='
   stop 'GEODAT'
endif

!     configure grid specs for raw data to rams grid (R) transfer

!     raw data grid (O)
TITLE=HFN(1:LB)//'HEADER'
LB=len_trim(TITLE)
CALL rams_f_open (29,title(1:lb),'FORMATTED','OLD','READ',0)
READ(29,*,end=1)IBLKSIZO,NO,ISBEGO,IWBEGO,offlat,offlon,h5name
1 continue
CLOSE(29)
DELTALLO=FLOAT(IBLKSIZO)/FLOAT(NO-1)

iodim=max(100000,4*no*no)
MOF=IODIM/(NO*NO)

   allocate(dato(iodim+mof+mof))

!     temp grid (Q) - smoothing only applied to topo
if(vnam(1:2).eq.'TO') then
      DELTAXQ=0.5*TOPTWVL(NGR)*DELTAXN(NGR)
else
      DELTAXQ=DELTAXN(NGR)
endif
NIQ=INT(FLOAT(N2-1)*DELTAXN(NGR)/DELTAXQ)+4
NJQ=INT(FLOAT(N3-1)*DELTAXN(NGR)/DELTAXQ)+4

!     interpollated raw data grid (P)
NP=MIN(10,MAX(1,INT(DELTAXQ/(DELTALLO*111000.))))
DELTAXP=DELTAXQ/FLOAT(NP)

! Need to line up XTN and YTN with this sub-domain. Do this by setting the
! addresses for the xt and yt arguments to 1 plus the x and y offsets associated
! with the sub-domain.
CALL sfcopqr (NO,MOF,NP,NIQ,NJQ,N2,N3,XTN(MI0(NGR)+1,NGR),YTN(MJ0(NGR)+1,NGR)  &
     ,polelat,polelon  &
     ,DELTALLO,DELTAXP,DELTAXQ,IBLKSIZO  &
     ,ISBEGO,IWBEGO,DATO(1),VT2DA,VT2DB,DATR  &
     ,OFN,offlat,offlon,VNAM,NGR,itopsflg(ngr),iz0flg(ngr),h5name)

deallocate(dato)

return
END SUBROUTINE geodat

!##############################################################################
Subroutine sfcopqr (no,mof,np,niq,njq,n2,n3,xt,yt,polelat,polelon  &
     ,deltallo,deltaxp,deltaxq,iblksizo  &
     ,isbego,iwbego,dato,datp,datq,datr  &
     ,ofn,offlat,offlon,vnam,ngr,itopsflg,iz0flg,h5name)

use grid_dims
use node_mod
use hdf5_utils

implicit none

integer :: no,mof,np,niq,njq,n2,n3,iblksizo,isbego,iwbego,ngr  &
          ,itopsflg,iz0flg
real :: dato(no,no,mof),datp(np,np),datq(niq,njq),datr(n2,n3)  &
         ,xt(n2),yt(n3)
real :: deltallo,deltaxp,deltaxq,offlat,offlon  
character(len=strl1) :: ofn,title3
character(len=3) :: title1,vnam
character(len=4) :: title2
character(len=*) :: h5name
logical :: l1,l2
character(len=strl1) :: fnmiss(maxmiss)
real, allocatable :: sdq(:,:),shaq(:,:),sdr(:,:),datre(:,:)
real, allocatable :: iso(:),iwo(:)

real :: polelat,polelon,xcentr,ycentr,glatp,glonp,rio_full,rjo_full  &
       ,xq,yq,xp,yp,wio1,wio2,wjo1,wjo2,sha,rha,rh2,sh,rh &
       ,xq1,yq1,xr,yr,rval,diff,difflcl
integer :: nmiss,nono,nofr,iof,iq,jq,ip,jp,iwoc,isoc,io1,io2,jo1,jo2 &
          ,lb,nn,isocpt,isocpo,iwocpo,iwocph,iwocpt,io_full,jo_full &
          ,iofr,jofr,ir,jr,is,js,i,j

integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

allocate (sdq(niq,njq),shaq(niq,njq),sdr(n2,n3),datre(n2,n3))
allocate (iso(mof),iwo(mof))

nmiss=0

nono=no*no
XCENTR=0.5*(XT(1)+XT(N2))
YCENTR=0.5*(YT(1)+YT(N3))
NOFR=0
DO IOF=1,MOF
   ISO(IOF)=0
   IWO(IOF)=0
ENDDO
DO JQ=1,NJQ
   DO IQ=1,NIQ
      XQ=(FLOAT(IQ)-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
      YQ=(FLOAT(JQ)-0.5*FLOAT(NJQ+1))*DELTAXQ+YCENTR
      DO JP=1,NP
         DO IP=1,NP
            XP=XQ+(FLOAT(IP)-0.5*FLOAT(NP+1))*DELTAXP
            YP=YQ+(FLOAT(JP)-0.5*FLOAT(NP+1))*DELTAXP

            CALL xy_ll (GLATP,GLONP,polelat,polelon,xp,yp)

            glatp = max(-89.9999,min(89.9999,glatp - offlat))
            glonp = max(-179.999,min(179.999,glonp - offlon))

            if (glonp >=  180.) glonp = glonp - 360.
            if (glonp <= -180.) glonp = glonp + 360.

            rio_full = (glonp - float(iwbego)) / deltallo
            rjo_full = (glatp - float(isbego)) / deltallo

            io_full = int(rio_full)
            jo_full = int(rjo_full)

            iwoc = (io_full / (no-1)) * iblksizo + iwbego
            isoc = (jo_full / (no-1)) * iblksizo + isbego

            wio2 = rio_full - float(io_full)
            wjo2 = rjo_full - float(jo_full)
           
            wio1 = 1. - wio2
            wjo1 = 1. - wjo2

            io1 = mod(io_full,no-1) + 1
            jo1 = mod(jo_full,no-1) + 1

            io2 = io1 + 1
            jo2 = jo1 + 1
            
            DO IOFR=1,NOFR
               JOFR=IOFR
               IF(ISO(IOFR).EQ.ISOC.AND.IWO(IOFR).EQ.IWOC)GO TO 10
            ENDDO

            ISOCPT=ABS(ISOC)/10
            ISOCPO=ABS(ISOC)-ISOCPT*10
            IWOCPH=ABS(IWOC)/100
            IWOCPT=(ABS(IWOC)-IWOCPH*100)/10
            IWOCPO=ABS(IWOC)-IWOCPH*100-IWOCPT*10
            IF(ISOC.GE.0) THEN
               WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'N'
            ELSE
               WRITE(TITLE1,'(2I1,A1)')ISOCPT,ISOCPO,'S'
            ENDIF
            IF(IWOC.GE.0) THEN
              !Saleeby(2009) Put in statement to solve SST dateline issues
              if(vnam.eq.'SST') then
               WRITE(TITLE2,'(4A1)')'0','0','0','E'
              else
               WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'E'
              endif
            ELSE
               WRITE(TITLE2,'(3I1,A1)')IWOCPH,IWOCPT,IWOCPO,'W'
            ENDIF

            !Determine data file name and see if it exists
            LB=len_trim(OFN)
            TITLE3=OFN(1:LB)//TITLE1//TITLE2
            LB=len_trim(TITLE3)
            title3=title3(1:lb)//'.h5'
            inquire(file=trim(title3),exist=l1,opened=l2)
     
            IF(.NOT.L1)THEN
               do nn=1,nmiss
                  if(trim(TITLE3).eq.fnmiss(nn)) goto 302
               enddo
               nmiss=nmiss+1
               fnmiss(nmiss)=trim(TITLE3)
302                 continue
               DATP(IP,JP)=0.
               GOTO 20
            ENDIF

            IF(NOFR.GE.MOF) THEN
               DO IOF=1,MOF
                  ISO(IOF)=0
                  IWO(IOF)=0
               ENDDO
               NOFR=0
            ENDIF
            NOFR=NOFR+1
            JOFR=NOFR

            !Open h5 data file and retrieve data
            ! Set the hyperslab specs to read the entire dataset
            mem_select%ndims = 2
            mem_select%dims(1:2)   = (/ no, no /)
            mem_select%block(1:2)  =  (/ no, no /)
            mem_select%count(1:2)  = (/ 1, 1 /)
            mem_select%offset(1:2) = (/ 0, 0 /)
            mem_select%stride(1:2) = (/ 1, 1 /)

            file_select%ndims = 2
            file_select%dims(1:2)   = (/ no, no /)
            file_select%block(1:2)  =  (/ no, no /)
            file_select%count(1:2)  = (/ 1, 1 /)
            file_select%offset(1:2) = (/ 0, 0 /)
            file_select%stride(1:2) = (/ 1, 1 /)

            CALL shdf5_open (title3,'R',iphdf5,h5_fid)
            CALL shdf5_irec (h5_fid,iphdf5,trim(h5name),mem_select &
                            ,file_select,rvara=dato(1,1,nofr))
            CALL shdf5_close (h5_fid)

            ISO(NOFR)=ISOC
            IWO(NOFR)=IWOC

10          CONTINUE

            datp(ip,jp)=wio1*(wjo1*dato(io1,jo1,jofr)   &
                             +wjo2*dato(io1,jo2,jofr))  &
                       +wio2*(wjo1*dato(io2,jo1,jofr)   &
                             +wjo2*dato(io2,jo2,jofr))
                                                          
20          CONTINUE
         ENDDO
      ENDDO

!           std dev for envelope orog and topo based zo schemes
      SHA=0.
      RHA=0.
      RH2=0.
      DO JP=1,NP
         SH=0.
         RH=0.
         DO IP=1,NP
            SH=MAX(SH,DATP(IP,JP))
            RH=RH+DATP(IP,JP)
            RH2=RH2+DATP(IP,JP)**2
         ENDDO
         SHA=SHA+SH/(2.*FLOAT(NP))
         RHA=RHA+RH
      ENDDO
      DATQ(IQ,JQ)=RHA/FLOAT(NP*NP)
      SDQ(IQ,JQ)=SQRT(max(0.,RH2/NP**2-DATQ(IQ,JQ)**2))
      DO IP=1,NP
         SH=0.
         DO JP=1,NP
            SH=MAX(SH,DATP(IP,JP))
         ENDDO
         SHA=SHA+SH/(2.*FLOAT(NP))
      ENDDO
      SHAQ(IQ,JQ)=SHA
      
   ENDDO
!         print*,'finished sfcopqr row jq = ',jq
ENDDO

!     envelope and zo schemes

if((vnam.eq.'TOP').and.  &
   (ITOPSFLG.eq.2).and.  &
   NP*NP.lt.8) print*,'Warning - '  &
   ,'trying to calc a std dev for: ',NP*NP,' points'
if((vnam.eq.'ZOT').and.  &
   IZ0FLG.eq.1.and.NP*NP.lt.8) print*,'Warning - '  &
   ,'trying to calc a std dev for: ',NP*NP,' points'


!           envelope orog and zo schemes
if(vnam.eq.'TOP')  &
   CALL topoq (NIQ,NJQ,DATQ,SDQ,SHAQ,DATRE,NGR,N2,N3)
if(vnam.eq.'ZOT')  &
   CALL zoq (NIQ,NJQ,DATQ,SDQ,NGR)

XQ1=(1.-0.5*FLOAT(NIQ+1))*DELTAXQ+XCENTR
YQ1=(1.-0.5*FLOAT(NJQ+1))*DELTAXQ+YCENTR
DO JR=1,N3
   DO IR=1,N2
      XR=(XT(IR)-XQ1)/DELTAXQ+1.
      YR=(YT(JR)-YQ1)/DELTAXQ+1.
      CALL gdtost (DATQ,NIQ,NJQ,XR,YR,RVAL)

      if(vnam.eq.'TOP') then
         DATR(IR,JR)=RVAL
      elseif(vnam.eq.'ZOT') then
         DATR(IR,JR)=MAX(0.,RVAL)
      else
         DATR(IR,JR)=MAX(0.,RVAL)
      endif
      
      CALL gdtost (SDQ,NIQ,NJQ,XR,YR,RVAL)
      SDR(IR,JR)=MAX(0.,RVAL)
   ENDDO
ENDDO

if(nmiss.gt.0) then
   print*,'-----------------------------------------------------'
   print*,'Input physiographical data file processing:'
   print*,'-----------------------------------------------------'
   print*,'  Input data blocks not found '  &
        ,' (data assumed to be zero):'
   do nn=1,nmiss
      print*,fnmiss(nn)
   enddo
   print*,'-----------------------------------------------------'
endif

!     check to find the largest change in topo height

if(vnam.eq.'TOP') then
   diff=0.
   difflcl=0.
   is=-999
   js=-999
   do j=2,n3-1
      do i=2,n2-1
         difflcl=max(difflcl,abs(datr(i,j)-datr(i-1,j)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i+1,j)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i,j-1)))
         difflcl=max(difflcl,abs(datr(i,j)-datr(i,j+1)))
         if(abs(diff-difflcl).gt.1.) then
            is=i
            js=j
         endif
         diff=max(diff,difflcl)
      enddo
   enddo
   write(6,100) ' Max d(topo) on grid @i,j=',ngr,is,js,diff
100     format(a,3i4,f8.1)
endif

deallocate(SDQ,SHAQ,SDR,DATRE)
deallocate(ISO,IWO)

return
END SUBROUTINE sfcopqr

!##############################################################################
Subroutine topoq (niq,njq,datq,sdq,shaq,datre,ngr,n2,n3)
                
use io_params
                
implicit none

integer :: niq,njq,ngr,n2,n3
real :: datq(niq,njq),sdq(niq,njq),shaq(niq,njq),datre(n2,n3)

integer :: iq,jq,jmin,imin,ire,jre,imax,jmax
real :: rad,count,total,remax,remin,average

!     orographic schemes

if(ITOPSFLG(ngr).lt.0) then                         ! No orography
   do jq=1,njq
      do iq=1,niq
         datq(iq,jq)=0.
!            print*,'None',iq,jq,datq(iq,jq)
      enddo
   enddo
   print *,'No orography'

elseif(ITOPSFLG(ngr).lt.0) then                      ! Average
   print *,'No orography enhancement applied'

   elseif(ITOPSFLG(ngr).eq.1) then                   ! Silhouette
      do jq=1,njq
         do iq=1,niq
            datq(iq,jq)=SHAQ(IQ,JQ)*toptenh(ngr)  &
                       +DATQ(IQ,JQ)*(1.-toptenh(ngr))
!                  print*,'Silhouette',iq,jq,datq(iq,jq)
         enddo
      enddo
      print *,'Silhouette Orography applied with'
      print *,'weighting = ',toptenh(ngr)

   elseif(ITOPSFLG(ngr).eq.2) then                   ! Envelope
      do jq=1,njq
         do iq=1,niq
            datq(iq,jq)=datq(iq,jq)+toptenh(ngr)*sdq(iq,jq)
!                  print*,'EO',iq,jq,datq(iq,jq)
         enddo
      enddo
      print *,'Envelope Orography applied with'
      print *,'enhancement = ',toptenh(ngr),' x std dev'

   else if(ITOPSFLG(ngr).ge.3) then                  ! Reflected Envelope

!        the radius we want to search for the current pts relative
!        height should correspond well to half the filtering wavelength
!        used on the topo (toptwvl)

         Rad=toptwvl(ngr)/2
      do jq=1,njq
         do iq=1,niq
            datre(iq,jq)=datq(iq,jq)
         enddo
      enddo
      do jq=1,njq
         do iq=1,niq
            count=0.
            total=0.
            remax=datre(iq,jq)
            remin=datre(iq,jq)
            jmin=jq-nint(Rad)
            imin=iq-nint(Rad)
            jmax=jq+nint(Rad)
            imax=iq+nint(Rad)
            do jre=max(1,jmin),min(njq,jmax)
               do ire=max(1,imin),min(niq,imax)
                  if((float((iq-ire)))**2  &
                    +(float((jq-jre)))**2.le.Rad**2) then
                   count=count+1.
                  total=total+datre(ire,jre)
                  remax=max(remax,datre(ire,jre))
                  remin=min(remin,datre(ire,jre))
               endif
            enddo
            enddo
         average=total/count
         if(remax.ne.remin)  &
                  datq(iq,jq)=datre(iq,jq)+(datre(iq,jq)-average)/  &
            ((remax-remin)/2)*toptenh(ngr)*sdq(iq,jq)
!               print*,'REO',iq,jq,datre(iq,jq),sdq(iq,jq),datq(iq,jq)
!               print*,'avg,n',average,count,remax,remin
      enddo
   enddo
   print *,'Reflected Envelope Orography applied with'
   print *,'enhancement = ',toptenh(ngr),' x std dev'
   print *,'and search radius (grid points) = ',Rad
endif

return
END SUBROUTINE topoq

!##############################################################################
Subroutine zoq (NIQ,NJQ,DATQ,SDQ,NGR)

use io_params
use mem_leaf

implicit none

integer :: NIQ,NJQ,NGR
real :: DATQ(NIQ,NJQ),SDQ(NIQ,NJQ)

integer :: iq,jq

!     topo base roughness length.


do jq=1,njq
   do iq=1,niq
      if(ITOPSFLG(ngr).lt.0) then  ! No orography
         datq(iq,jq)=zrough
      elseif(iz0flg(ngr).eq.1) then
         datq(iq,jq)=min(z0fact*sdq(iq,jq),z0max(NGR))
      else
         datq(iq,jq)=zrough
      endif
   enddo
enddo
if(ITOPSFLG(ngr).lt.0) then  ! No orography
   print *,'No orography'
else
   print *,'Subgrid terrain roughness applied with'
   print *,'factor  = ',z0fact,' x std dev'
   print *,'maximum = ',z0max(NGR)
endif

return
END SUBROUTINE zoq

