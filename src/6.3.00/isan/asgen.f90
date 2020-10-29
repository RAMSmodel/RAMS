!##############################################################################
Subroutine isan_driver ()

use node_mod
use isan_coms
use mem_grid
use io_params
use hdf5_utils
use mem_leaf
use mem_varinit, only:snudcent

implicit none

character(len=3) :: csuff
character(len=strl1) :: locfn
integer :: ifm,icm,ng,k
integer*8 :: h5_fid
integer :: iphdf5

if (nmachs .gt. 1) then
  iphdf5 = 1
else
  iphdf5 = 0
endif

! Topo is read from surface files.

do ng=1,ngrids
   CALL newgrid (ng)
   CALL top_read (ng) 
enddo

!  Setup RAMS horizontal and vertical grid structure.

CALL grid_setup (2)

!  Allocate RAMS grid arrays where data analysis will be put 
!     for output and feedback.  Old way was to use "A"

maxix=maxval(mmxp(1:nigrids))   
maxiy=maxval(mmyp(1:nigrids))   
maxiz=maxval(mmzp(1:nigrids))   

do ngrid=1,nigrids
   !For upper air grids
   allocate(is_grids(ngrid)%rr_u (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_v (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_t (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_r (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_p (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_ug (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_vg (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_tg (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_rg (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_pg (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_pi0 (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_th0 (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0 (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0u(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0v(mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

   !For condensate nudging option
   allocate(is_grids(ngrid)%rr_cond (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

   !For soil/snow initialization
   allocate(is_grids(ngrid)%rr_soilmoist1 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soilmoist2 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soiltemp1 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soiltemp2 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_snowmass  (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_snowdepth (mmxp(ngrid),mmyp(ngrid)))
enddo
allocate(rr_scr1(maxix*maxiy*maxiz))
allocate(rr_scr2(maxix*maxiy*maxiz))
allocate(rr_scr3(maxix*maxiy*maxiz))
allocate(rr_vt2da(maxix*maxiy))

! Do inventory of input file names, determine which times to process

CALL isan_file_inv (iyear1,imonth1,idate1,itime1,timmax)

! Start main processing loop for each data time

do natime=1,npdates
   if (iproc_flag(natime,1) == 0) cycle
   CALL date_unmake_big (iyear,imonth,idate,ihour,iproc_dates(natime))
   print*
   print*,'========================================================'
   print*,'ISAN natime: ',natime,npdates,iproc_dates(natime)
   print*,'ISAN processing time: ',iyear,imonth,idate,ihour
   print*
   ihour=ihour/100

   innpr=iproc_names(natime,1)

   !Read in and process reanalysis pressure level data from dprep "dp" files
   if(iproc_flag(natime,2).eq.1) then
      CALL pressure_stage (mmxp(1),mmyp(1),grid_g(1)%glat(1,1),grid_g(1)%glon(1,1))
      print*,'';print*,'Finished pressure stage';print*,''
   endif

! Isentropic/sigma-z analysis to all requested RAMS grids
   
   do ngrid=1,nigrids

      ! Find number of sigma-z levels

      do k=1,mmzp(ngrid)
         if(ztn(k,ngrid) > topsigz) goto 100
         sigz(k)=ztn(k,ngrid)
      enddo
      k=mmzp(ngrid)+1
      100  continue
      nsigz=k-1


!         Allocate memory for isentropic analysis
!         --------------------------------------------------------
      print*,'Allocating RAMS polar/isentropic grid-'  &
             ,ngrid,mmxp(ngrid),mmyp(ngrid),nisn

      allocate(pi_u(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_v(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_p(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_s(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_r(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_scra(mmxp(ngrid),mmyp(ngrid),nisn))
      allocate(pi_scrb(mmxp(ngrid),mmyp(ngrid),nisn))

!         Allocate memory for sigma-z analysis 
!         --------------------------------------------------------

      print*,'Allocating RAMS polar/sigmaz grid-'  &
             ,ngrid,mmxp(ngrid),mmyp(ngrid),nsigz

      allocate(ps_u(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_v(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_p(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_t(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_r(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_scra(mmxp(ngrid),mmyp(ngrid),nsigz))
      allocate(ps_scrb(mmxp(ngrid),mmyp(ngrid),nsigz))

!         Allocate memory for surface analysis 
!         --------------------------------------------------------

      print*,'Allocating RAMS polar/surface grid-'  &
             ,ngrid,mmxp(ngrid),mmyp(ngrid)

      allocate(rs_u(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_v(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_p(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_t(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_r(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_s(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_top(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_qual(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_soilmoist1(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_soilmoist2(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_soiltemp1(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_soiltemp2(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_snowmass(mmxp(ngrid),mmyp(ngrid)))
      allocate(rs_snowdepth(mmxp(ngrid),mmyp(ngrid)))


      ! Do isentropic and sigma-z analysis
      
      if(iszstage == 1) then

         print*,'Doing isentropic-sigz analysis'
         CALL isnstage ()

         ! Output isentropic file if desired

         if(ioflgisz == 1)then
            write(csuff,'(a1,i1)') 'g',ngrid
            CALL makefnam (locfn,varpfx,0,iyear,imonth,idate,  &
                           ihour*100,'I',csuff,'h5')
            CALL shdf5_open (locfn,'W',iphdf5,h5_fid,iclobber)
            CALL isenio (h5_fid,iphdf5,'OUT',mmxp(ngrid),mmyp(ngrid),ngrid)
            CALL sigzio (h5_fid,iphdf5,'OUT',mmxp(ngrid),mmyp(ngrid),ngrid)
            CALL shdf5_close (h5_fid)
         endif

      elseif(ivrstage == 1) then
         print*,'Doing varfile stage and reading in data'
         write(csuff,'(a1,i1)') 'g',ngrid
         CALL makefnam (locfn,varpfx,0,iyear,imonth,idate,  &
              ihour*100,'I',csuff,'h5')
         CALL shdf5_open (locfn,'R',iphdf5,h5_fid,iclobber)
         CALL isenio (h5_fid,iphdf5,'IN',mmxp(ngrid),mmyp(ngrid),ngrid)
         CALL sigzio (h5_fid,iphdf5,'IN',mmxp(ngrid),mmyp(ngrid),ngrid)
         CALL shdf5_close (h5_fid)
      endif
           
      ! Prepare "varfiles"

      if(ivrstage == 1) then
         CALL makevarf (ngrid)
      endif

      deallocate(pi_u,pi_v,pi_p,pi_s,pi_r,pi_scra,pi_scrb)
      deallocate(ps_u,ps_v,ps_p,ps_t,ps_r,ps_scra,ps_scrb)
      deallocate(rs_u,rs_v,rs_p,rs_t,rs_r,rs_s,rs_top,rs_qual)
      deallocate(rs_soilmoist1,rs_soilmoist2,rs_soiltemp1,rs_soiltemp2 &
                ,rs_snowmass,rs_snowdepth)

   enddo

   ! Do the nesting feedback and write out the "varfiles"

   if(ivrstage == 1) then

      if(nfeedvar == 1 .and. nigrids > 1) then

         ! fill reference states for all grids
         CALL varfile_refstate (mmzp(1),mmxp(1),mmyp(1)  &
              ,is_grids(1)%rr_t(1,1,1),is_grids(1)%rr_p(1,1,1)  &
              ,is_grids(1)%rr_pi0(1,1,1),is_grids(1)%rr_th0(1,1,1)  &
              ,is_grids(1)%rr_r(1,1,1),is_grids(1)%rr_dn0(1,1,1)  &
              ,is_grids(1)%rr_dn0u(1,1,1),is_grids(1)%rr_dn0v(1,1,1)  &
              ,grid_g(1)%topt(1,1),grid_g(1)%rtgt(1,1)  &
              ,ztn(1,1),ztop,piref(1,1),thref(1,1),dnref(1,1),rtref(1,1))
         is_grids(1)%rr_p(1:mmzp(1),1:mmxp(1),1:mmyp(1)) =  &
                   is_grids(1)%rr_p  (1:mmzp(1),1:mmxp(1),1:mmyp(1)) &
                  -is_grids(1)%rr_pi0(1:mmzp(1),1:mmxp(1),1:mmyp(1))

         ! Initialize fine grids by interpolating from the parent coarse grid
         do ifm = 1,nigrids
           icm = nxtnest(ifm)
           if (icm .ge. 1) then
             CALL newgrid (ifm)
             CALL interp_fine_grid_sfcvar (ifm,icm,IFG_MKVF_INIT)
           endif
         enddo

         ! Feed back u,v,pi,t,and rt fields
         do ifm=nigrids,2,-1
            icm=nxtnest(ifm)
            if (icm == 0) cycle
            CALL varfile_nstfeed (ifm,icm,mmzp(ifm),mmxp(ifm),mmyp(ifm) &
                  ,mmzp(icm),mmxp(icm),mmyp(icm)  &
                  ,nnstbot(icm),nnsttop(icm))
         enddo

         do ifm=1,nigrids
         is_grids(ifm)%rr_p(1:mmzp(ifm),1:mmxp(ifm),1:mmyp(ifm)) =  &
            is_grids(ifm)%rr_p  (1:mmzp(ifm),1:mmxp(ifm),1:mmyp(ifm)) &
           +is_grids(ifm)%rr_pi0(1:mmzp(ifm),1:mmxp(ifm),1:mmyp(ifm))
         enddo
      endif

      !Set condensate nudging variable to null value for
      !standard varfile creation.
      do ifm=1,nigrids
        is_grids(ifm)%rr_cond = -9999.
      enddo

      !Quality control on soil moisture and temperature to make sure all
      !land grid cells have valid soil values and not missing data.
      if(isoildat == 1 .or. snudcent > 0.) then
       do ifm = 1,nigrids
        CALL newgrid (ifm)
        CALL sfc_read (ifm)
        CALL soil_quality_control (ifm,mmxp(ifm),mmyp(ifm),npatch &
          ,is_grids(ifm)%rr_soilmoist1(1,1),is_grids(ifm)%rr_soilmoist2(1,1) &
          ,is_grids(ifm)%rr_soiltemp1(1,1),is_grids(ifm)%rr_soiltemp2(1,1) &
          ,grid_g(ifm)%glat(1,1),grid_g(ifm)%glon(1,1) &
          ,leaf_g(ifm)%patch_area(1,1,1),leaf_g(ifm)%leaf_class(1,1,1))
       enddo
      endif

      !Write out the varfiles now that all the prep is finished.
      if(ioflgvar == 1) then
         CALL write_varf ()
      endif

   endif

enddo

return
END SUBROUTINE isan_driver

!##############################################################################
Subroutine soil_quality_control (ifm,n2,n3,npat,soil_moist_bot,soil_moist_top &
   ,soil_temp_bot,soil_temp_top,glat,glon,patch_area,leaf_class)

use mem_grid, only:iprntstmt,print_msg

implicit none

integer :: n2,n3,i,j,ii,jj,minii,minjj,ifm,npat
real, dimension(n2,n3) :: glat,glon,soil_moist_top,soil_moist_bot      &
  ,soil_temp_top,soil_temp_bot,soil_moist_top_temp,soil_moist_bot_temp &
  ,soil_temp_top_temp,soil_temp_bot_temp
real, dimension(n2,n3,npat) :: patch_area,leaf_class
real :: disty,distx,dist,mindist,maxsoil

 if(iprntstmt>=1 .and. print_msg)then
  print*,'------------------------------------------------------------------'
  print*,'RUNNING QC ON SOIL MOISTURE AND SOIL TEMPERATURE, GRID:',IFM
  print*,'------------------------------------------------------------------'
 endif

 !Do some initial quality control
 do j = 1,n3
  do i = 1,n2
   !Set NULL as in the write to varfiles
   if(soil_moist_top(i,j) >= 1E30) soil_moist_top(i,j) = -9999.
   if(soil_moist_bot(i,j) >= 1E30) soil_moist_bot(i,j) = -9999.
   if(soil_temp_top(i,j) >= 1E30) soil_temp_top(i,j) = -9999.
   if(soil_temp_bot(i,j) >= 1E30) soil_temp_bot(i,j) = -9999.
   !Check for high values
   if(soil_moist_top(i,j) > 0.49 .or. soil_moist_bot(i,j) > 0.49) then
    soil_moist_top(i,j) = .49
    soil_moist_bot(i,j) = .49
   endif
   !Check for non-zero low values
   if((soil_moist_top(i,j) > 0.0 .and. soil_moist_top(i,j) < 0.01) .or. &
      (soil_moist_bot(i,j) > 0.0 .and. soil_moist_bot(i,j) < 0.01)) then
    soil_moist_top(i,j) = 0.0101
    soil_moist_bot(i,j) = 0.0101
   endif
   !Set default temp array values
   soil_moist_top_temp(i,j) = -9999.
   soil_moist_bot_temp(i,j) = -9999.
   soil_temp_top_temp(i,j) = -9999.
   soil_temp_bot_temp(i,j) = -9999.
  enddo
 enddo

 !Interpolate to what amounts to coastal point where data is null, but
 !RAMS says there is land that needs valid soil moisture values.
 do j = 1,n3
  do i = 1,n2
   !if(i==7.and.j==1)print*,'st1',patch_area(i,j,2),nint(leaf_class(i,j,2)) &
   !,soil_moist_top(i,j),soil_moist_bot(i,j)
   if((soil_moist_top(i,j) < 0.01 .or. soil_moist_top(i,j) > 0.49 .or. &
       soil_moist_bot(i,j) < 0.01 .or. soil_moist_bot(i,j) > 0.49 ) .and. &
       patch_area(i,j,2) > 0.001 .and. &
       nint(leaf_class(i,j,2)) >= 2 .and. nint(leaf_class(i,j,2)) <= 21) then
    mindist=9.99e12
    minii=i
    minjj=j
    do jj = 1,n3
     do ii = 1,n2
      if(soil_moist_top(ii,jj) > 0.01 .and. &
         nint(leaf_class(ii,jj,2)) >= 2 .and. nint(leaf_class(ii,jj,2)) <= 21) then
       disty=abs(glat(i,j)-glat(ii,jj))*111.12
       distx=abs(glon(i,j)-glon(ii,jj))*111.12*cos((glat(i,j)+glat(ii,jj))/2.0*3.14/180.)
       dist=sqrt(distx**2.0+disty**2.0)
       if(dist<mindist)then
        mindist=dist
        minii=ii
        minjj=jj
       endif
      endif
     enddo
    enddo
    ii=minii
    jj=minjj
    soil_moist_top_temp(i,j) = soil_moist_top(ii,jj)
    soil_moist_bot_temp(i,j) = soil_moist_bot(ii,jj)
    soil_temp_top_temp(i,j)  = soil_temp_top(ii,jj)
    soil_temp_bot_temp(i,j)  = soil_temp_bot(ii,jj)
   endif
  enddo
 enddo

 !Copy max values back to passed variables to send back to varfiles
 do j = 1,n3
  do i = 1,n2
    soil_moist_top(i,j) = max(soil_moist_top_temp(i,j),soil_moist_top(i,j))
    soil_moist_bot(i,j) = max(soil_moist_bot_temp(i,j),soil_moist_bot(i,j))
    soil_temp_top(i,j)  = max(soil_temp_top_temp(i,j),soil_temp_top(i,j))
    soil_temp_bot(i,j)  = max(soil_temp_bot_temp(i,j),soil_temp_bot(i,j))
  enddo
 enddo

 !Quality check our interpolation and look for bad points
 !This could indicate that varfile was not created.
 do j = 1,n3
  do i = 1,n2
   maxsoil=max(soil_moist_bot(i,j),soil_moist_top(i,j))
   if(patch_area(i,j,2) > 0.001 .and. (maxsoil >.50 .or. maxsoil <.01)) then
    print*,'WARNING: Some land areas have invalid soil moisture values:'
    print*,'Need varfiles for all grids at initial time for valid soil on all grids.'
    print*,ifm,i,j,soil_moist_bot(i,j),soil_moist_top(i,j),patch_area(i,j,2)
    stop
   endif
  enddo
 enddo

return
END SUBROUTINE soil_quality_control

!##############################################################################
Subroutine opspec4 ()

use mem_grid
use isan_coms

implicit none

integer :: ifaterr,iwarerr,infoerr

! This routine checks the option specifications in the $MODEL_GRIDS
!    $ISAN_ISENTROPIC namelists for consistency.

IFATERR=0
IWARERR=0
INFOERR=0

! Don't allow NIGRIDS <= NGRIDS

if(nigrids.gt.ngrids)then
  print*,' FATAL - NIGRIDS must be <= NGRIDS'
  IFATERR=IFATERR+1
endif


! Stop the run if there are any fatal errors.  List how many
!   warning and informative errors.

PRINT*,' -----------opspec4--------------------------'
PRINT*,' FATAL     errors - ',IFATERR
PRINT*,' WARNING   errors - ',IWARERR
PRINT*,' INFORM  messages - ',INFOERR
PRINT*,' -----------------------------------------------'

IF(IFATERR.GT.0) STOP 'OPSPEC4'

return
END SUBROUTINE opspec4

!##############################################################################
Subroutine nudh_driver ()

use node_mod
use isan_coms
use mem_grid
use mem_varinit, only:nnudfl,fnames_nud,nnudfiles
use io_params
use hdf5_utils

implicit none

integer :: ng,lnf,ifm,icm,checkhist

!Allocate grid data type since almost all model memory is not used for ISAN.
print*,'start HISTORY-VAR grid allocation'
! If grid_g has already been allocated, deallocate it so
! we can start over. See comments in isan_driver() above.
if (allocated(grid_g)) then
  do ng=1,ngrids
    CALL dealloc_grid (grid_g(ng)) 
  enddo
  deallocate(grid_g)
endif

allocate(grid_g(ngrids))
do ng=1,ngrids
   CALL alloc_grid (grid_g(ng),mmxp(ng),mmyp(ng)) 
enddo

!Topo is read from surface files.
do ng=1,ngrids
   CALL newgrid (ng)
   CALL top_read (ng) 
enddo

!Setup RAMS horizontal and vertical grid structure.
CALL grid_setup (2)

!Allocate RAMS grid arrays
maxix=maxval(mmxp(1:nigrids))   
maxiy=maxval(mmyp(1:nigrids))   
maxiz=maxval(mmzp(1:nigrids))   

do ngrid=1,nigrids  
   allocate(is_grids(ngrid)%rr_u    (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_v    (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_t    (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_r    (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_p    (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_pi0  (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0  (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_th0  (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0u (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_dn0v (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

   allocate(is_grids(ngrid)%rr_cond (mmzp(ngrid),mmxp(ngrid),mmyp(ngrid)))

   allocate(is_grids(ngrid)%rr_soilmoist1 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soilmoist2 (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soiltemp1  (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_soiltemp2  (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_snowmass   (mmxp(ngrid),mmyp(ngrid)))
   allocate(is_grids(ngrid)%rr_snowdepth  (mmxp(ngrid),mmyp(ngrid)))
enddo
allocate(rr_scr1(maxix*maxiy*maxiz))
allocate(rr_scr2(maxix*maxiy*maxiz))
allocate(rr_scr3(maxix*maxiy*maxiz))
allocate(rr_vt2da(maxix*maxiy))

!Do inventory of input file names. All will be processed.
CALL nud_file_inv (var_hfile,iyear1,imonth1,idate1,itime1)

!Start main processing loop for each data time
histloop: do nnudfl=1,nnudfiles

 !Read the history files and extract needed data for varfiles
 CALL hvfiles (nnudfl,checkhist)

 !Check to see if all history grids are present at a given time
 if(checkhist==0) cycle histloop

 ! Update varfile pressure variable for grid 1
  is_grids(1)%rr_p(1:mmzp(1),1:mmxp(1),1:mmyp(1)) =  &
      is_grids(1)%rr_p  (1:mmzp(1),1:mmxp(1),1:mmyp(1)) &
     +is_grids(1)%rr_pi0(1:mmzp(1),1:mmxp(1),1:mmyp(1))

 ! Initialize the fine grids via interpolation
 do ifm = 1,nigrids
   icm = nxtnest(ifm)
   if (icm .ge. 1) then
     CALL newgrid (ifm)
     CALL interp_fine_grid_sfcvar (ifm,icm,IFG_MKHF_INIT)
   endif
 enddo

 !Varfiles expect the 6 soil and snow fields, so we must put in dummy values
 do ngrid=1,nigrids
   is_grids(ngrid)%rr_soilmoist1 = -9999.
   is_grids(ngrid)%rr_soilmoist2 = -9999.
   is_grids(ngrid)%rr_soiltemp1  = -9999.
   is_grids(ngrid)%rr_soiltemp2  = -9999.
   is_grids(ngrid)%rr_snowmass   = -9999.
   is_grids(ngrid)%rr_snowdepth  = -9999.
 enddo

 !Set some filetime naming conventions prior to writing
 lnf=len_trim(fnames_nud(nnudfl))
 read(fnames_nud(nnudfl)(lnf-25:lnf-9),20) iyear,imonth,idate,ihour
 20 format(i4,1x,i2,1x,i2,1x,i6)
 ihour=ihour/100

 !Write the varfile
 CALL write_varf ()

enddo histloop

return
END SUBROUTINE nudh_driver

!##############################################################################
Subroutine hvfiles (nnud,checkhist)

use isan_coms
use an_header
use mem_grid
use ref_sounding
use micphys, only:level
use grid_struct
use mem_leaf, only:nvegpat
use hdf5_utils
use mem_varinit
use node_mod

implicit none

integer :: nnud,ngrids1,nzg1,nzs1,npatch1,nvegpat1,ierr,ng,nc &
           ,ie,maxarr1,maxarr2,ngr,maxx1,maxy1,maxz1,npts,nv &
           ,ihtran1,nzpg1,checkhist,nvloopstart,dim_type
integer, external :: cio_i,cio_f
integer,save :: iunhd=11
integer, allocatable, dimension(:) :: nnxp1,nnyp1,nnzp1
real :: ztop1,polelat1,polelon1
real, allocatable, dimension(:,:) :: xmn1,xtn1,ymn1,ytn1,zmn1,ztn1,topt1
real, allocatable, dimension(:) :: scr,scr1,scr2,scr3
character(len=strl1) :: hnameinh
character(len=2) :: cng
character(len=1) :: cgrid
real, allocatable, dimension(:) :: u01dn1,v01dn1,rt01dn1,th01dn1  &
                         ,pi01dn1,dn01dn1
logical :: exists ! File existence

integer*8 :: h5_fid
integer :: iphdf5
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

type (head_table), allocatable :: hr_table(:)

type(grid_def), allocatable :: grdefh(:)
type(grid_def), allocatable :: grdefn(:)

print*,'Doing HISTORY FILE NUDGING option.'

33 format(a30,2i5,3x,a18,i8)
34 format(a,i3,3x,a,i3,a)

!Get length of history file name
nc=len_trim(fnames_nud(nnud))

!Open RAMS history header file
CALL rams_f_open (iunhd,fnames_nud(nnud),'FORMATTED','OLD','READ',0)

!Get grid specs for grid comparison
ie=cio_i(iunhd,1,'ngrids',ngrids1,1)

allocate (nnxp1(ngrids1),nnyp1(ngrids1),nnzp1(ngrids1))

ie=cio_i(iunhd,1,'nnxp',nnxp1,ngrids1)
ie=cio_i(iunhd,1,'nnyp',nnyp1,ngrids1)
ie=cio_i(iunhd,1,'nnzp',nnzp1,ngrids1)
ie=cio_i(iunhd,1,'npatch',npatch1,1)
ie=cio_i(iunhd,1,'nvegpat',nvegpat1,1)
ie=cio_i(iunhd,1,'nzg',nzg1,1)
ie=cio_i(iunhd,1,'nzs',nzs1,1)
ie=cio_f(iunhd,1,'ztop',ztop1,1)
ie=cio_f(iunhd,1,'polelat',polelat1,1)
ie=cio_f(iunhd,1,'polelon',polelon1,1)
ie=cio_i(iunhd,1,'ihtran',ihtran1,1)

!Make sure new grid top does not exceed old grid top
if(ztop > ztop1) then
  if(print_msg) then
   print*,'New interpolated grid top higher than past grid top!'
   print*,'Cannot interpolate. Extrapolation not recommended.'
   print*,'Past top:',ztop1, ' New top:',ztop
  endif
  stop
endif

if (ihtran /= ihtran1) then
  print*,'GRID TYPE (Cartesian or Polar Stereographic) must be the same'
  print*,'between history grid and current grid. May need to consider'
  print*,'wind rotation and such if we allow a change.'
  stop
endif  

!Find maximum size of any array on history file. Allocate scratch arrays.
maxarr1=0
maxarr2=0
maxx1=maxval(nnxp1(1:ngrids1))
maxy1=maxval(nnyp1(1:ngrids1))
maxz1=maxval(nnzp1(1:ngrids1))
do ngr=1,ngrids1
   maxarr1=max(maxarr1,nnxp1(ngr)*nnyp1(ngr)*nnzp1(ngr))
   maxarr2=max(maxarr2,nnxp1(ngr)*nnyp1(ngr))
enddo
allocate (scr(maxarr1),scr1(maxarr1))

!Allocate and read in grid specifications
allocate(xmn1(maxx1,ngrids1),xtn1(maxx1,ngrids1))
allocate(ymn1(maxy1,ngrids1),ytn1(maxy1,ngrids1))
allocate(zmn1(maxz1,ngrids1),ztn1(maxz1,ngrids1))
do ngr=1,ngrids1
   write(cng,'(i2.2)') ngr
   ie=cio_f(iunhd,1,'xmn'//cng,xmn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'xtn'//cng,xtn1(1,ngr),nnxp1(ngr))
   ie=cio_f(iunhd,1,'ymn'//cng,ymn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'ytn'//cng,ytn1(1,ngr),nnyp1(ngr))
   ie=cio_f(iunhd,1,'zmn'//cng,zmn1(1,ngr),nnzp1(ngr))
   ie=cio_f(iunhd,1,'ztn'//cng,ztn1(1,ngr),nnzp1(ngr))
enddo

!Read variable header info
rewind(iunhd)

!Top of header file. Total number of variables on all grids.
read(iunhd,*) nvbtab
!Allocate space for header file variable information
allocate (hr_table(nvbtab))
!Reader info on each variables in header file.
!(ex. Name  IfAPointer  Dimension   Grid   #GridPoints)
!(ex.  UP       0           3        1        49000)
do nv=1,nvbtab
  read(iunhd,*)              &
     hr_table(nv)%string     &
    ,hr_table(nv)%npointer   &
    ,hr_table(nv)%idim_type  &
    ,hr_table(nv)%ngrid      &
    ,hr_table(nv)%nvalues
enddo

!Check to see that all history grids are present at this time
checkhist=1
do ngr=1,ngrids1
  write(cgrid,'(i1)') ngr
  hnameinh=fnames_nud(nnudfl)(1:nc-9)//'-g'//cgrid//'.h5'
  inquire(file=hnameinh,exist=exists)
  if(.not. exists)then
   checkhist=0
   return
  endif
enddo

!Allocate history topography and get TOPT
allocate (topt1(maxarr2,ngrids1))
do ngr=1,ngrids1
   write(cgrid,'(i1)') ngr
   hnameinh=fnames_nud(nnud)(1:nc-9)//'-g'//cgrid//'.h5'
   CALL update_for_hist_grid (1,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
   if (nmachs .gt. 1) then
     iphdf5 = 1
   else
     iphdf5 = 0
   endif
   CALL shdf5_open (hnameinh,'R',iphdf5,h5_fid)
   CALL shdf5_set_hs_select (2,'R',ngr &
          ,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'TOPT' &
          ,mem_select,file_select,rvara=topt1(1,ngr))
   CALL shdf5_close (h5_fid)
   CALL update_for_hist_grid (2,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
enddo

!********* START CHECK GRID MATCHING *****************************
! Set a flag array (for each grid on the history file) to determine:
!  = 1 = This grid is identical to a current grid
!  = 0 = This grid is different.

igrid_match(1:maxgrds,1:maxgrds)=0

!Allocate grid structures
allocate (grdefn(nigrids))
allocate (grdefh(ngrids1))

print*,'###################################################################'
print*,'# History Grids:',ngrids1
print*,'# Current Grids:',nigrids

!Define the History Grid
do ngr=1,ngrids1
   CALL alloc_grid_def ( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr) )
   CALL fill_grid_def  ( grdefh(ngr),nnxp1(ngr),nnyp1(ngr),nnzp1(ngr)   &
                       ,nzg1,nzs1,npatch1,nvegpat1,polelat1,polelon1  &
                       ,xtn1(1,ngr),xmn1(1,ngr),ytn1(1,ngr),ymn1(1,ngr) &
                       ,ztn1(1,ngr),zmn1(1,ngr),topt1(1,ngr) )
enddo
!Define the Current Initial Grid
do ngr=1,nigrids
   CALL alloc_grid_def ( grdefn(ngr),mmxp(ngr),mmyp(ngr),mmzp(ngr) )
   CALL fill_grid_def  ( grdefn(ngr),mmxp(ngr),mmyp(ngr),mmzp(ngr)  &
                       ,nzg,nzs,npatch,nvegpat,polelat,polelon    &
                       ,xtn(1,ngr),xmn(1,ngr),ytn(1,ngr),ymn(1,ngr) &
                       ,ztn(1,ngr),zmn(1,ngr),grid_g(ngr)%topt )
enddo

! See if the history grids match any of the new grids...assuming 1:1 grid number
!   correspondence for now
do ngr=1,ngrids1
 do ng=1,nigrids
   CALL compare_grid_def (grdefh(ngr),grdefn(ng),'nud_update',ierr,ngr,ng)
   if (ierr /= 0) then
      ! No match...
      print*,'These grids do NOT match (Hist,New):',ngr,ng
   else
      ! We have a match...
      igrid_match(ngr,ng)=1
      print*,'These grids DO match (Hist,New):',ngr,ng
   endif
 enddo
enddo
print*,'###################################################################'

!********* END CHECK GRID MATCHING *****************************

!********* START READING IN THE DATA ***************************
! Finally, process the fields...

!*************************************************************************
!The code below copies or interpolates all history grids to all current 
!grids in the model if a current grid fits fully inside a history grid. 
!This is only done for the 5 main nudging fields.
!*************************************************************************

!*************************************************************************
!Loop over history grids "ngr" (ngrids1 = number of grids in history file)
!*************************************************************************
do ngr = 1,ngrids1

  !Get history file name and open the file
  write(cgrid,'(i1)') ngr
  hnameinh = fnames_nud(nnud)(1:nc-9)//'-g'//cgrid//'.h5'

  !*************************************************************************
  !Loop over current run grids "ng" (nigrids = number of grids in current
  !run that we will make varfiles for.
  !*************************************************************************
  do ng=1,nigrids

   if (xtn(1,ng) < xtn1(1,ngr) .or. xtn(mmxp(ng),ng) > xtn1(nnxp1(ngr),ngr) .or. &
       ytn(1,ng) < ytn1(1,ngr) .or. ytn(mmyp(ng),ng) > ytn1(nnyp1(ngr),ngr) ) then

    print*,'###################################################################'
    print 34,' Current Grid:',ng,'outside bounds of History Grid:',ngr,': NO Interp'
    print*,'Will not interpolate to a new grid bigger than hist grid.'
    if(ngr==1)stop 'Current grid larger than all history grids.'
    print*,'###################################################################'

   else 

    print*,'###################################################################'
    print 34,' Current Grid:',ng,' within bounds of History Grid:',ngr,': Interp'
    print*,'###################################################################'

    !Number of grid points for the current run variable
    npts=mmzp(ng)*mmxp(ng)*mmyp(ng)
    allocate (scr2(npts),scr3(npts))

    !Loop over history variables
    nvloopstart=0
    varloop: do nv=1,nvbtab

     !Cycle if the history variable is not on the grid we are looping over
     if(ngr /= hr_table(nv)%ngrid) cycle varloop

     !Assign or interpolate primary nudging variables
     if(hr_table(nv)%string=="UP"    .or. &
        hr_table(nv)%string=="VP"    .or. &
        hr_table(nv)%string=="THETA" .or. &
        hr_table(nv)%string=="RTP"   .or. &
        hr_table(nv)%string=="PP") then

      !Copy current grid structure to temporary arrays and tell
      !HDF5 to look for history file grid structure to read.
      CALL update_for_hist_grid (1,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
      if (nmachs .gt. 1) then
        iphdf5 = 1
      else
        iphdf5 = 0
      endif
      !Open file, Read data from the history files, Close file
      CALL shdf5_open (hnameinh,'R',iphdf5,h5_fid)
      dim_type = hr_table(nv)%idim_type
      CALL shdf5_set_hs_select (dim_type,'R',ngr &
                 ,mem_select,file_select,file_chunks)
      CALL shdf5_irec (h5_fid,iphdf5,trim(hr_table(nv)%string) &
                 ,mem_select,file_select,rvara=scr1)
      CALL shdf5_close (h5_fid)
      !Copy current grid structure back to main grid arrays.
      CALL update_for_hist_grid (2,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)

      !Arrange data for RAMS to work with to copy or interpolate
      CALL unarrange (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),scr1,scr)

      !Copy or interpolate history grids to current grids
      if(igrid_match(ngr,ng)==1) then
          !If grids match, just copy fields over to nudging arrays
          if(iprntstmt>=1)print 33,'nud_update: copy: ' &
                         ,ngr,ng,hr_table(nv)%string,npts
          CALL atob (npts,scr(1),scr2(1))
      else
          !Otherwise do interpolation to different grid
          if(iprntstmt>=1)print 33,'nud_update: interp: ' &
                         ,ngr,ng,hr_table(nv)%string,npts
          CALL hi_interp (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,1,scr(1)   &
               ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)         &
               ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1              &
               ,mmzp(ng),mmxp(ng),mmyp(ng),1,1,scr2(1)                  &
               ,ng,ngr,hr_table(nv)%string,3)
      endif
     endif

     !Assign or interpolate sum of RCP,RDP,RRP,RPP,RSP,RAP,RGP,RHP
     if(hr_table(nv)%string=="RCP" .or. &
        hr_table(nv)%string=="RDP" .or. &
        hr_table(nv)%string=="RRP" .or. &
        hr_table(nv)%string=="RPP" .or. &
        hr_table(nv)%string=="RSP" .or. &
        hr_table(nv)%string=="RAP" .or. &
        hr_table(nv)%string=="RGP" .or. &
        hr_table(nv)%string=="RHP") then
      if(nvloopstart==0) then
        CALL azero (npts,scr2(1))
        nvloopstart=1
      endif

      !Copy current grid structure to temporary arrays and tell
      !HDF5 to look for history file grid structure to read.
      CALL update_for_hist_grid (1,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)
      if (nmachs .gt. 1) then
        iphdf5 = 1
      else
        iphdf5 = 0
      endif
      !Open file, Read data from the history files, Close file
      CALL shdf5_open (hnameinh,'R',iphdf5,h5_fid)
      dim_type = hr_table(nv)%idim_type
      CALL shdf5_set_hs_select (dim_type,'R',ngr &
                 ,mem_select,file_select,file_chunks)
      CALL shdf5_irec (h5_fid,iphdf5,trim(hr_table(nv)%string) &
                 ,mem_select,file_select,rvara=scr1)
      CALL shdf5_close (h5_fid)
      !Copy current grid structure back to main grid arrays.
      CALL update_for_hist_grid (2,ngrids1,nnzp1,nnxp1,nnyp1,nzg1,nzs1,npatch1)

      !Arrange data for RAMS to work with to copy or interpolate
      CALL unarrange (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),scr1,scr)

      !Copy or interpolate history grids to current grids
      if(igrid_match(ngr,ng)==1) then
          !If grids match, just copy fields over to nudging arrays
         if(iprntstmt>=1)print 33,'nud_update: accum: ' &
                         ,ngr,ng,hr_table(nv)%string,npts
          CALL nud_cond_accum (npts,scr(1),scr2(1))
      else
          !Otherwise do interpolation to different grid
          if(iprntstmt>=1)print 33,'nud_update: interp: ' &
                         ,ngr,ng,hr_table(nv)%string,npts
          CALL hi_interp (nnzp1(ngr),nnxp1(ngr),nnyp1(ngr),1,1,scr(1) &
               ,xmn1(1,ngr),xtn1(1,ngr),ymn1(1,ngr),ytn1(1,ngr)       &
               ,zmn1(1,ngr),ztn1(1,ngr),topt1(1,ngr),ztop1            &
               ,mmzp(ng),mmxp(ng),mmyp(ng),1,1,scr3(1)                &
               ,ng,ngr,hr_table(nv)%string,3)
          CALL nud_cond_accum (npts,scr3(1),scr2(1))
      endif
     endif

     if(hr_table(nv)%string=="UP")    CALL atob (npts,scr2(1),is_grids(ng)%rr_u(1,1,1))
     if(hr_table(nv)%string=="VP")    CALL atob (npts,scr2(1),is_grids(ng)%rr_v(1,1,1))
     if(hr_table(nv)%string=="THETA") CALL atob (npts,scr2(1),is_grids(ng)%rr_t(1,1,1))
     if(hr_table(nv)%string=="RTP")   CALL atob (npts,scr2(1),is_grids(ng)%rr_r(1,1,1))
     if(hr_table(nv)%string=="PP")    CALL atob (npts,scr2(1),is_grids(ng)%rr_p(1,1,1))

     !Copy the accumulation each time a hydrometeor species is accesssed since
     !we do not know which will be the last active species to copy.
     if(hr_table(nv)%string=="RCP" .or. &
        hr_table(nv)%string=="RDP" .or. &
        hr_table(nv)%string=="RRP" .or. &
        hr_table(nv)%string=="RPP" .or. &
        hr_table(nv)%string=="RSP" .or. &
        hr_table(nv)%string=="RAP" .or. &
        hr_table(nv)%string=="RGP" .or. &
        hr_table(nv)%string=="RHP") then
          CALL atob (npts,scr2(1),is_grids(ng)%rr_cond(1,1,1))
     endif

    enddo varloop ! do looping over nudged variables

    deallocate (scr2,scr3)

   endif !if current grid fully contained in history grid

  enddo !do looping over current run grids

enddo !do looping over history grids

!*************************************************************************
!The following code for reference state variables is set up for
!grid-1. If there is more than one grid, the model will interpolate
!from grid-1 to additional grids. This is done in rdint.f90 routine "initlz".
!*************************************************************************

!Prepare 1D reference sounding for history grid-1 and apply to new grid-1. 
!Really only need to do this to get the isan Pressure variable.
 nzpg1=nnzp1(1)
 allocate(u01dn1(nzpg1), v01dn1(nzpg1),rt01dn1(nzpg1)  &
        ,th01dn1(nzpg1),pi01dn1(nzpg1),dn01dn1(nzpg1) )

 cng='01'
 ie=cio_f(iunhd,1,'u01dn'//cng,  u01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'v01dn'//cng,  v01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'pi01dn'//cng,pi01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'th01dn'//cng,th01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'dn01dn'//cng,dn01dn1(1),nzpg1)
 ie=cio_f(iunhd,1,'rt01dn'//cng,rt01dn1(1),nzpg1)

 CALL htint (nzpg1,th01dn1,ztn1(1,1),mmzp(1),th01dn(1,1),ztn(1,1))
 CALL htint (nzpg1,u01dn1 ,ztn1(1,1),mmzp(1),u01dn(1,1) ,ztn(1,1))
 CALL htint (nzpg1,v01dn1 ,ztn1(1,1),mmzp(1),v01dn(1,1) ,ztn(1,1))
 CALL htint (nzpg1,pi01dn1,ztn1(1,1),mmzp(1),pi01dn(1,1),ztn(1,1))
 CALL htint (nzpg1,dn01dn1,ztn1(1,1),mmzp(1),dn01dn(1,1),ztn(1,1))
 if (level .ge. 1) then
    CALL htint (nzpg1,rt01dn1,ztn1(1,1),mmzp(1),rt01dn(1,1),ztn(1,1))
 else
    rt01dn(1:mmzp(1),1) = 0.
 endif

 u01dn(1,1) = u01dn(2,1)
 v01dn(1,1) = v01dn(2,1)
 rt01dn(1,1) = rt01dn(2,1)
 th01dn(1,1) = th01dn(2,1)

 !Close the input history header file
 close(iunhd)

 !Compute 3d reference state for grid 1
 CALL newgrid (1)
 CALL refs3d (mmzp(1),mmxp(1),mmyp(1)  &
  ,is_grids(1)%rr_pi0  (1,1,1)  ,is_grids(1)%rr_dn0  (1,1,1)  &
  ,is_grids(1)%rr_dn0u (1,1,1)  ,is_grids(1)%rr_dn0v (1,1,1)  &
  ,is_grids(1)%rr_th0  (1,1,1)  ,grid_g(1)%topt      (1,1)    &
  ,grid_g(1)%rtgt      (1,1)                                    ) 

 !Deallocate temporary arrays
 deallocate(u01dn1,v01dn1,rt01dn1,th01dn1,pi01dn1,dn01dn1)

 !Deallocate temporary arrays
 deallocate(topt1,hr_table,xmn1,ymn1,zmn1,xtn1,ytn1,ztn1,scr,scr1 &
   ,nnxp1,nnyp1,nnzp1)

 do ngr=1,ngrids1
    CALL dealloc_grid_def (grdefh(ngr))
 enddo
 do ngr=1,nigrids
    CALL dealloc_grid_def (grdefn(ngr))
 enddo

return
END SUBROUTINE hvfiles

!##############################################################################
Subroutine nud_cond_accum (n,a,b)

implicit none

integer :: n,i
real :: a(n),b(n)

do i=1,n
 b(i)=a(i)+b(i)
enddo

return
END SUBROUTINE nud_cond_accum

!##############################################################################
Subroutine nud_file_inv (hfilin,iyear1,imonth1,idate1,itime1)

use mem_varinit

implicit none

character(len=*) :: hfilin
integer :: iyear1,imonth1,idate1,itime1,nf,lnf,nhftot &
          ,inyear,inmonth,indate,inhour
character(len=strl1) :: fnames(maxnudfiles),hpref
character(len=14)  :: itotdate
real(kind=8) :: secs_init,secs_nud

! Get abs seconds of run start

CALL date_abs_secs2 (iyear1,imonth1,idate1,itime1*100,secs_init)

! Go through history files and make inventory

nhftot=-1
hpref=hfilin

CALL rams_filelist (fnames,trim(hpref)  &
         //'????-??-??-??????-head.txt',nhftot)

if(nhftot > maxnudfiles) then
   print*,'too many nud history files'
   stop 'lots_of_nud_history'
endif

nnudfiles=0
do nf=1,nhftot
   lnf=len_trim(fnames(nf))
   read(fnames(nf)(lnf-25:lnf-9),20) inyear,inmonth,indate,inhour
   20 format(i4,1x,i2,1x,i2,1x,i6)

   CALL date_make_big (inyear,inmonth,indate,inhour,itotdate)

   nnudfiles=nnudfiles+1
   fnames_nud(nnudfiles)=fnames(nf)
   itotdate_nud(nnudfiles)=itotdate
   
   CALL date_abs_secs2 (inyear,inmonth,indate,inhour,secs_nud)
   nud_times(nnudfiles)=secs_nud - secs_init

enddo

CALL rams_dintsort (nnudfiles,itotdate_nud,fnames_nud)

!  start printing section

print*,' '
print*,'-------------------------------------------------------------'
print*,'-----------  History Nudging Input File Inventory -----------'
print*,'-------------------------------------------------------------'
do nf=1,nnudfiles
   print*,  itotdate_nud(nf),'   ',nud_times(nf)  &
           ,trim(fnames_nud(nf))
enddo
print*,'------------------------------------------------------'

return
END SUBROUTINE nud_file_inv
