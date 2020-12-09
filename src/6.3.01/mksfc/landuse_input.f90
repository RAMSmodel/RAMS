!##############################################################################
Subroutine patch_array_size (npq,deltax  &
   ,ivegtflg,ivegtfn,isoilflg,isoilfn,ndviflg,ndvifn)

use grid_dims

implicit none

integer :: npq,ivegtflg,isoilflg,ndviflg
character(len=*) :: ivegtfn,isoilfn,ndvifn
character(len=strl1) :: h5name
real :: deltax
integer :: iblksiz,no,isbeg,iwbeg
real :: offlat,offlon,deltall,deltallo_min

! Find highest resolution dataset among all that are to be used

deltallo_min = 1.e20
if (ivegtflg == 1) then
   CALL read_header (ivegtfn,iblksiz,no,isbeg  &
      ,iwbeg,offlat,offlon,deltall,'veg',h5name)
   deltallo_min = min(deltallo_min,deltall)
endif

if (isoilflg == 1) then
   CALL read_header (isoilfn,iblksiz,no,isbeg  &
      ,iwbeg,offlat,offlon,deltall,'soil',h5name)
   deltallo_min = min(deltallo_min,deltall)
endif

if (ndviflg == 1) then
   CALL read_header (ndvifn,iblksiz,no,isbeg  &
      ,iwbeg,offlat,offlon,deltall,'ndvi',h5name)
   deltallo_min = min(deltallo_min,deltall)
endif

npq = min(10,max(1,nint(deltax / (deltallo_min * 111000.))))

return
END SUBROUTINE patch_array_size

!##############################################################################
Subroutine patch_latlon (n2,n3,xt,yt,deltax,polelat,polelon)

use mem_mksfc

implicit none

integer :: n2,n3
real :: xt(n2),yt(n3),polelat,polelon
integer :: jr,jp,ip,ir
real :: yp,xp,deltax,deltaxp

! Fill arrays with offset latitudes and longitudes of all p points

   deltaxp = deltax / float(npq)
   do jr = 1,n3
      do jp = 1,npq
         yp = yt(jr) + (float(jp) - .5 * float(npq+1)) * deltaxp
         do ir = 1,n2
            do ip = 1,npq
               xp = xt(ir) + (float(ip) - .5 * float(npq+1)) * deltaxp
               CALL xy_ll (glatp(ip,jp,ir,jr),glonp(ip,jp,ir,jr)  &
                  ,polelat,polelon,xp,yp)
            enddo
         enddo
      enddo
   enddo

return
END SUBROUTINE patch_latlon

!##############################################################################
Subroutine landuse_opqr (n2,n3,mzg,npat,nvegpat  &
   ,ivegtfn,isoilfn,ndvifn,cndvifil,iaction  &
   ,soil_text,patch_area,leaf_class,veg_ndvif)

! Routine landuse_opr reads in one or more landuse data types and defines
! and fills subgrid patches from them.  Currently, landuse (vegetation) class,
! soil textural class, and ndvi value are implemented.  In this version of
! landuse_opr, landuse class must be used if any datasets are used.  Patch
! areas are determined from landuse class alone.  If soil textural class data
! is used, the most dominant type occurring in each landuse-class-defined
! patch is assigned to that patch.  If ndvi data is used, the average value
! for each landuse-class-defined patch is assigned to that patch.

use mem_mksfc
use rconstants

implicit none

integer :: n2,n3,mzg,npat,nvegpat
character(len=*) :: iaction,ivegtfn,isoilfn,ndvifn,cndvifil
integer :: nmiss
real, dimension(mzg,n2,n3,npat) :: soil_text
real, dimension(n2,n3,npat) :: patch_area,leaf_class,veg_ndvif
character(len=strl1) :: fnmiss(maxmiss),h5name 
integer, parameter :: maxdatq=32,nsoil=12
integer :: datq,ngrdpix(0:maxdatq,2),datq_pat  &
          ,datsoil,soil_count      &
          ,jr,jp,ir,ip,iqv,ing,maxdq,jng,jng1,ngstor1    &
          ,ngstor2,npatpixs,nwat,ipat,idatq,isoil,jsoil,k      &
          ,no_veg ,iblksizo_veg ,isbego_veg ,iwbego_veg   &
          ,no_soil,iblksizo_soil,isbego_soil,iwbego_soil  &
          ,no_ndvi,iblksizo_ndvi,isbego_ndvi,iwbego_ndvi
integer, dimension(0:maxdatq) :: sumpix
integer, dimension(0:maxdatq,nsoil) :: datq_soil
real :: deltallo_veg ,deltallo_soil ,deltallo_ndvi  &
       ,offlat_veg   ,offlat_soil   ,offlat_ndvi    &
       ,offlon_veg   ,offlon_soil   ,offlon_ndvi    &
       ,fracwat,plpp
real, dimension(0:maxdatq) :: sumndvi

nmiss=0

if (iaction .eq. 'veg') then

   CALL read_header (ivegtfn,iblksizo_veg,no_veg,isbego_veg  &
      ,iwbego_veg,offlat_veg,offlon_veg,deltallo_veg,'veg',h5name)

   CALL fill_datp (n2,n3,no_veg,iblksizo_veg,isbego_veg,iwbego_veg  &
      ,offlat_veg,offlon_veg,deltallo_veg,ivegtfn,iaction  &
      ,nmiss,fnmiss,h5name)

! 9/30/97:  Carry out the first translation of the input DATP values into a
! condensed set called DATQ_patch.  The range of DATQ_patch values represents
! the total variety of land surface conditions (patches) to be allowed for 
! the present simulation, and may be a broader class than the LEAF-2 
! vegetation classes for which all the vegetation physical parameters are 
! defined.  For example, two different DATQ_patch classes may be mapped to 
! the same LEAF-2 vegetation class, but be initialized with different soil 
! moistures or soil types, and therefore require different patches.

! Fill datq_patch (patch class) values from input landuse dataset.
! Currently, this data serves as the primary criterion for defining patches.

   do jr = 1,n3
      do ir = 1,n2

         do iqv = 0,maxdatq
            ngrdpix(iqv,1) = 0     ! Initialize counter for datq pixels
            ngrdpix(iqv,2) = iqv   ! Initialize array of consecutive datq values 
         enddo

         do jp = 1,npq
            do ip = 1,npq
               CALL datp_datq (datp(ip,jp,ir,jr),datq_patch(ip,jp,ir,jr))
               datq = datq_patch(ip,jp,ir,jr)
               ngrdpix(datq,1) = ngrdpix(datq,1) + 1
            enddo
         enddo

! Sort values of ngrdpix by prevalence for non-water patches (datq .ge. 2)
         jng1=2 !Variable initialized
         do ing = 2,maxdatq
            maxdq = -1
            do jng = ing,maxdatq
               if (ngrdpix(jng,1) .gt. maxdq) then
                  jng1 = jng
                  maxdq = ngrdpix(jng,1)
               endif
            enddo
            ngstor1 = ngrdpix(ing,1)
            ngstor2 = ngrdpix(ing,2)
            ngrdpix(ing,1) = ngrdpix(jng1,1)
            ngrdpix(ing,2) = ngrdpix(jng1,2)
            ngrdpix(jng1,1) = ngstor1
            ngrdpix(jng1,2) = ngstor2
         enddo

! Fill patches numbered 2 through nvegpat+1 with the nvegpat most prevalent
!    nonwater landuse types.  Count pixels for these patches for normalization
!    of total grid cell area

         npatpixs = 1
         nwat = ngrdpix(0,1) + ngrdpix(1,1)

         if (nwat .lt. npq*npq) then
            npatpixs = 0
            do ipat = 2,nvegpat+1
               datq = ngrdpix(ipat,2)
               leaf_class(ir,jr,ipat) = float(datq)
               npatpixs = npatpixs + ngrdpix(ipat,1)
            enddo
         endif

         fracwat = float(nwat) / float(npq * npq)
         plpp = (1. - fracwat) / float(npatpixs)
         patch_area(ir,jr,1) = fracwat

         if (ngrdpix(0,1) .ge. ngrdpix(1,1)) then
            leaf_class(ir,jr,1) = 0.
         else
            leaf_class(ir,jr,1) = 1.
         endif

         do ipat = 2,nvegpat+1
            patch_area(ir,jr,ipat) = plpp * float(ngrdpix(ipat,1))
         enddo

      enddo

   enddo

elseif (iaction .eq. 'soil') then

   CALL read_header (isoilfn,iblksizo_soil,no_soil,isbego_soil  &
      ,iwbego_soil,offlat_soil,offlon_soil,deltallo_soil,'soil',h5name)

   CALL fill_datp (n2,n3,no_soil,iblksizo_soil,isbego_soil,iwbego_soil  &
      ,offlat_soil,offlon_soil,deltallo_soil,isoilfn,iaction  &
      ,nmiss,fnmiss,h5name)

   do jr = 1,n3
      do ir = 1,n2
         if (patch_area(ir,jr,1) .le. .9999) then

            do idatq = 0,maxdatq
               do isoil = 1,nsoil
                  datq_soil(idatq,isoil) = 0     ! Initialize counter for datq pixels
               enddo
            enddo

            do jp = 1,npq
               do ip = 1,npq

! Fill datq_soil values as secondary criterion.  This is for finding dominant
! soil class for each datq class.

                  CALL datp_datsoil (datp(ip,jp,ir,jr),datsoil)
                  datq_pat = datq_patch(ip,jp,ir,jr)
                  datq_soil(datq_pat,datsoil)  &
                     = datq_soil(datq_pat,datsoil) + 1

               enddo
            enddo

            do ipat = 2,nvegpat+1
               if (patch_area(ir,jr,ipat) .ge. .0001) then   

                  datq_pat = nint(leaf_class(ir,jr,ipat))

! Find isoil value for which soil_tab(datq_pat,isoil) is a maximum
                  jsoil=1 !Variable initialized
                  soil_count = 0
                  do isoil = 1,nsoil
                     if (datq_soil(datq_pat,isoil) .gt. soil_count) then
                        soil_count = datq_soil(datq_pat,isoil)
                        jsoil = isoil
                     endif
                  enddo
                  
! For now, assume single level of input soil data (e.g., FAO) and
! fill all soil levels with this value.

                  do k = 1,mzg
                     soil_text(k,ir,jr,ipat) = float(jsoil)
                  enddo

               endif
            enddo

         endif
      enddo
   enddo

elseif (iaction .eq. 'ndvi') then

   CALL read_header (ndvifn,iblksizo_ndvi,no_ndvi,isbego_ndvi  &
      ,iwbego_ndvi,offlat_ndvi,offlon_ndvi,deltallo_ndvi,'ndvi',h5name)
      
   CALL fill_datp (n2,n3,no_ndvi,iblksizo_ndvi,isbego_ndvi,iwbego_ndvi  &
      ,offlat_ndvi,offlon_ndvi,deltallo_ndvi,cndvifil,iaction  &
      ,nmiss,fnmiss,h5name)

   do jr = 1,n3
      do ir = 1,n2
         if (patch_area(ir,jr,1) .le. .9999) then

            do idatq = 0,maxdatq
               sumndvi(idatq) = 0.  ! initialize ndvi sum
               sumpix(idatq) = 0    ! initialize ndvi pixel count
            enddo

            do jp = 1,npq
               do ip = 1,npq

! Fill datq_ndvi values to compute average value for each datq class.

                  datq_pat = datq_patch(ip,jp,ir,jr)
                  sumndvi(datq_pat) = sumndvi(datq_pat) + datp(ip,jp,ir,jr)
                  sumpix(datq_pat) = sumpix(datq_pat) + 1

               enddo
            enddo

            do ipat = 2,nvegpat+1
                  datq_pat = nint(leaf_class(ir,jr,ipat))
                  if(sumpix(datq_pat) > 0) then
                    veg_ndvif(ir,jr,ipat) =  max( .05, &
                            sumndvi(datq_pat)/ sumpix(datq_pat)  )
                  else
                    veg_ndvif(ir,jr,ipat) = .05
                  endif
            enddo

         endif
      enddo
   enddo

endif

if(nmiss.gt.0) then
   print*,'-----------------------------------------------------'
   print*,'Input surface characteristic data file processing:',iaction
   print*,'-----------------------------------------------------'
   print*,'  Some input data blocks not found (data assumed to be ocean or default):'
!   do nn=1,nmiss
!      print*,trim(fnmiss(nn))
!   enddo
   print*,trim(fnmiss(1))
   print*,'  plus a total of ',nmiss,'files'
   print*,'-----------------------------------------------------'
endif

return
END SUBROUTINE landuse_opqr

!##############################################################################
Subroutine read_header (ofn,iblksizo,no,isbego,iwbego,offlat,offlon,deltallo  &
   ,ifield,h5name)

use grid_dims
use node_mod

implicit none

integer :: iblksizo,no,isbego,iwbego,lb
real :: offlat,offlon,deltallo
character(len=*) :: ofn,ifield,h5name
character(len=strl1) :: title

lb = len_trim(ofn)
if (lb .le. 0) then
   print*,'| ',ifield,' input data prefix incorrect !'
   print*,'|  file prefix:',ofn(1:lb)
   print*,'====================================================='
   stop 'landuse-file'
endif

title = ofn(1:lb)//'HEADER'
lb = len_trim(title)

CALL rams_f_open (29,title(1:lb),'FORMATTED','OLD','READ',0)
read(29,*,end=1) iblksizo,no,isbego,iwbego,offlat,offlon,h5name
1 continue
close(29)
deltallo = float(iblksizo) / float(no-1)

print*, 'read_header1 ',ifield
print*, 'read_header2 ',iblksizo,no,isbego,iwbego,offlat,offlon,deltallo

return
END SUBROUTINE read_header

!##############################################################################
Subroutine fill_datp (n2,n3,no,iblksizo,isbego,iwbego  &
   ,offlat,offlon,deltallo,ofn,iaction,nmiss,fnmiss,h5name)

use mem_mksfc
use hdf5_utils
use mem_grid, only:iprntstmt
use node_mod

implicit none

character(len=*) :: ofn,iaction,fnmiss(*),h5name
integer :: nmiss
integer :: n2,n3,no,iblksizo,isbego,iwbego,isoc,iwoc  &
   ,isocpt,isocpo,iwocph,iwocpt,iwocpo,lb,io,jo  &
   ,ir,jr,ip,jp,ind,nc3,nc2,j3d,j2d,j1d,ind1,ind2,io1,jo1  &
   ,ifile_max,jfile_max,ifile,jfile,missing,ptab,ptab0,idatp,nn
real :: rio,rjo,rno,offlat,offlon  &
       ,glatp1,glonp1,deltallo,wio2,wjo2,wio1,wjo1
character(len=3) :: title1
character(len=4) :: title2
character(len=strl1) :: title3
logical :: l1,l2
integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select

!if (nmachs .gt. 1) then
!  iphdf5 = 1
!else
!  iphdf5 = 0
!endif

! The set of files that are being read by this routine have separate files for
! each data sample point (lat,lon). This makes is so that the sub domains will
! require different sets of files since they exist in different lat,lon locations.
! Because of this, we can't open the files in parallel mode (all processes have
! open the same set of files) without requiring that all processes open all files
! over the entire domain. It turns out that multiple processes can access the same
! HDF5 file simultaneously in non-parallel mode as long as they are only doing reads. 
!
! Figure out the set of files for each process and read just those files in non-
! parallel mode. This means to give patch_latlon() the dimensions and locations
! of the sub domain, and always set iphdf5 to zero.
iphdf5 = 0

! Compute number of files in input dataset that span all latitudes and
! longitudes on earth.  Allocate nump and numpind arrays to this size and 
! initialize to 0.
! Allocate ptable array.

rno = float(no)
ifile_max = 360 / iblksizo
jfile_max = 180 / iblksizo

if (iaction == 'ndvi') then
   allocate (dato(no,no))
else
   allocate (cdato(no,no),idato(no,no))
endif

allocate (nump    (ifile_max,jfile_max)  &
         ,numpind (ifile_max,jfile_max)  &
         ,numpind1(ifile_max,jfile_max)  &
         ,numpind2(ifile_max,jfile_max)  &
         ,ptable  (npq*npq*n2*n3))

do jfile = 1,jfile_max
   do ifile = 1,ifile_max
      nump(ifile,jfile) = 0
      numpind(ifile,jfile) = 0
   enddo
enddo

! Get file index (ifile,jfile) within full dataset and count number of p 
! points (nump) that occur in each file

do jr = 1,n3
   do jp = 1,npq
      do ir = 1,n2
         do ip = 1,npq

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))

            !Saleeby(2010):For glonp1, use the min/max func
            !glonp1 = glonp(ip,jp,ir,jr) - offlon
            glonp1 = max(-179.9999,min(179.9999,glonp(ip,jp,ir,jr) - offlon))

            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
            jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1

            nump(ifile,jfile) = nump(ifile,jfile) + 1

            !write(6,202) ip,jp,ir,jr,ifile,jfile,nump(ifile,jfile)  &
            !,glatp1,glonp1
            !202 format('at2',7i5,2f10.4)
            
         enddo
      enddo
   enddo
enddo

! Set up array index values for ptable array

ind = 1
do jfile = 1,jfile_max
   do ifile = 1,ifile_max
      numpind1(ifile,jfile) = ind
      numpind2(ifile,jfile) = ind
      ind = ind + nump(ifile,jfile)
   enddo
enddo

! Fill ptable array

nc3 = n2 * npq * npq
nc2 = npq * npq

do jr = 1,n3
   j3d = (jr - 1) * nc3
   do ir = 1,n2
      j2d = (ir - 1) * nc2
      do jp = 1,npq
         j1d = (jp - 1) * npq
         do ip = 1,npq

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))

            !Saleeby(2010):For glonp1, use the min/max func
            !glonp1 = glonp(ip,jp,ir,jr) - offlon
            glonp1 = max(-179.9999,min(179.9999,glonp(ip,jp,ir,jr) - offlon))

            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            ifile = int((glonp1 - float(iwbego)) / float(iblksizo)) + 1
            jfile = int((glatp1 - float(isbego)) / float(iblksizo)) + 1

            ind = numpind2(ifile,jfile)
            ptable(ind) = j3d + j2d + j1d + ip
            numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1
            
         enddo
      enddo
   enddo
enddo

! Read files and extract data

do jfile = 1,jfile_max
   do ifile = 1,ifile_max
   
      ind1 = numpind1(ifile,jfile)
      ind2 = numpind2(ifile,jfile)
   
      if (ind2 .gt. ind1) then
         isoc = (jfile - 1) * iblksizo + isbego
         iwoc = (ifile - 1) * iblksizo + iwbego

! Construct filename

         isocpt = abs(isoc) / 10
         isocpo = abs(isoc) - isocpt*10
         iwocph = abs(iwoc) / 100
         iwocpt = (abs(iwoc) - iwocph * 100) / 10
         iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
         
         if (isoc .ge. 0) then
            write(title1,'(2i1,a1)') isocpt,isocpo,'N'
         else
            write(title1,'(2i1,a1)') isocpt,isocpo,'S'
         endif
         
         if (iwoc .ge. 0) then
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
         else
            write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
         endif

         lb = len_trim(ofn)
         title3 = ofn(1:lb)//title1//title2
         lb = len_trim(title3)
         title3=trim(title3)//'.h5'
         inquire(file=trim(title3),exist=l1,opened=l2)

! Read file or set missing flag to 1

         if(iprntstmt>=1)print*,'inquire file: ',title3(1:lb),l1,l2,ir,jr,ip,jp
         if (l1) then
            missing = 0
            if(iprntstmt>=1)print*,'getting file: ',title3

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

            if (iaction == 'ndvi') then
                  CALL shdf5_open (title3,'R',iphdf5,h5_fid)
                  CALL shdf5_irec (h5_fid,iphdf5,trim(h5name),mem_select &
                                  ,file_select,rvara=dato)
                  CALL shdf5_close (h5_fid)
            else
                  CALL shdf5_open (title3,'R',iphdf5,h5_fid)
                  CALL shdf5_irec (h5_fid,iphdf5,trim(h5name),mem_select &
                                  ,file_select,ivara=idato)
                  CALL shdf5_close (h5_fid)
            endif
         else
            do nn=1,nmiss
               if(trim(title3(1:lb)) == trim(fnmiss(nn)) ) goto 302
            enddo
            nmiss=nmiss+1
            fnmiss(nmiss)=title3(1:lb)
302         continue
            missing = 1
         endif

         do ind = ind1,ind2-1

            ptab = ptable(ind)         
            ptab0 = ptab - 1
            jr = ptab0 / nc3 + 1
            j3d = (jr - 1) * nc3
            ir = (ptab0 - j3d) / nc2 + 1
            j2d = (ir - 1) * nc2
            jp = (ptab0 - j3d - j2d) / npq + 1
            j1d = (jp - 1) * npq
            ip = ptab - j3d - j2d - j1d

            glatp1 = max(-89.9999,min(89.9999,glatp(ip,jp,ir,jr) - offlat))

            !Saleeby(2010):For glonp1, use the min/max func
            !glonp1 = glonp(ip,jp,ir,jr) - offlon
            glonp1 = max(-179.9999,min(179.9999,glonp(ip,jp,ir,jr) - offlon))

            if (glonp1 .ge.  180.) glonp1 = glonp1 - 360.
            if (glonp1 .le. -180.) glonp1 = glonp1 + 360.

            rio = (glonp1 - float(iwoc)) / deltallo + 1.
            !if( abs(glonp1 - float(iwoc)) < 1.e-5 ) rio = 1.
            rjo = (glatp1 - float(isoc)) / deltallo + 1.
            !if( abs(glatp1 - float(isoc)) < 1.e-5 ) rjo = 1.

            if (rio .lt. .9 .or. rio .gt. rno+.1 .or.  &
                rjo .lt. .9 .or. rjo .gt. rno+.1) then
                print*, 'rio,rjo out of range ',rio,rjo,ip,jp,ir,jr
                stop 45
            endif

            if (missing .eq. 0) then

               if (iaction .eq. 'veg' .or. iaction .eq. 'soil') then
                  io = nint(rio)
                  jo = nint(rjo)
                  idatp = idato(io,jo)
                  datp(ip,jp,ir,jr) = float(mod(idatp+256,256))

               elseif (iaction .eq. 'ndvi') then
                  io1 = max(1,min(int(rio),no-1))
                  jo1 = max(1,min(int(rjo),no-1))
                  wio2 = rio - float(io1)
                  wjo2 = rjo - float(jo1)
                  wio1 = 1. - wio2
                  wjo1 = 1. - wjo2
                  if(jo1==0) then
                     print*,'jp0:',ip,jp,ir,jr,io1,jo1,rio,rjo,rno
                     print*,'jp0:',rjo, glatp1 , float(isoc)  &
                        ,glatp1-float(isoc), deltallo
                     stop
                  endif
                  datp(ip,jp,ir,jr) =  &
                     wio1 * (wjo1 * dato(io1  ,jo1  )   &
                          +  wjo2 * dato(io1  ,jo1+1))  &
                   + wio2 * (wjo1 * dato(io1+1,jo1  )   &
                          +  wjo2 * dato(io1+1,jo1+1))
               endif

            else

               datp(ip,jp,ir,jr) = 0.

            endif

         enddo
      endif
   enddo
enddo

if (iaction == 'ndvi') then
   deallocate (dato)
else
   deallocate (cdato,idato)
endif

deallocate(nump,numpind,numpind1,numpind2,ptable)

return
END SUBROUTINE fill_datp
