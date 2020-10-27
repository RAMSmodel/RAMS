!##############################################################################
Subroutine isenio (h5_fid,iphdf5,inout,n1,n2,igrid)

use isan_coms
use hdf5_utils

implicit none

integer :: n1,n2,igrid,npts,nx3,ny3,ninn,l,levnn(maxisn)
integer*8 :: h5_fid
integer :: iphdf5
character(len=*) :: inout
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if(inout == 'IN') THEN

   ! scalar values
   CALL shdf5_set_hs_select (1,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_year',mem_select,file_select,ivars=iyy)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_month',mem_select,file_select,ivars=imm)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_date',mem_select,file_select,ivars=idd)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_hour',mem_select,file_select,ivars=ihh)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_nx',mem_select,file_select,ivars=nx3)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_ny',mem_select,file_select,ivars=ny3)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_nisn',mem_select,file_select,ivars=ninn)

   ! Vector, length is ninn which was just read in with the prior shdf5_irec call
   ! To indicate to shdf5_set_h5_select that the variable is a vector of
   ! length ninn, use negative ninn.
   CALL shdf5_set_hs_select (-ninn,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_levth',mem_select,file_select,ivara=levnn)

   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nisn) then
      print*,'Isentropic stage grid dimensions do not match'
      print*,'   configuration file on read !'
      print*,' File dimens - ',nx3,ny3,ninn
      print*,' Run  dimens - ',n1,n2,nisn
      stop 'IO3-2'
   endif

   npts=n1*n2*nisn
   ! ISAN 3D isentropic variables
   CALL shdf5_set_hs_select (8,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_u',mem_select,file_select,rvara=pi_u)
   CALL vmissr (pi_u,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_v',mem_select,file_select,rvara=pi_v)
   CALL vmissr (pi_v,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_s',mem_select,file_select,rvara=pi_s)
   CALL vmissr (pi_p,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_p',mem_select,file_select,rvara=pi_p)
   CALL vmissr (pi_s,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_r',mem_select,file_select,rvara=pi_r)
   CALL vmissr (pi_r,npts,1e30,-9998.)

   npts=n1*n2

   ! ISAN 2D variables
   CALL shdf5_set_hs_select (2,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_u',mem_select,file_select,rvara=rs_u)
   CALL vmissr (rs_u,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_v',mem_select,file_select,rvara=rs_v)
   CALL vmissr (rs_v,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_p',mem_select,file_select,rvara=rs_p)
   CALL vmissr (rs_p,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_t',mem_select,file_select,rvara=rs_t)
   CALL vmissr (rs_t,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_r',mem_select,file_select,rvara=rs_r)
   CALL vmissr (rs_r,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_s',mem_select,file_select,rvara=rs_s)
   CALL vmissr (rs_s,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_topo',mem_select,file_select,rvara=rs_top)
   CALL vmissr (rs_top,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_qual',mem_select,file_select,rvara=rs_qual)
   CALL vmissr (rs_qual,npts,1e30,-9998.)

   CALL shdf5_irec (h5_fid,iphdf5,'sfc_soilmoist1',mem_select,file_select &
                   ,rvara=rs_soilmoist1)
   CALL vmissr (rs_soilmoist1,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_soilmoist2',mem_select,file_select &
                   ,rvara=rs_soilmoist2)
   CALL vmissr (rs_soilmoist2,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_soiltemp1',mem_select,file_select  &
                   ,rvara=rs_soiltemp1)
   CALL vmissr (rs_soiltemp1,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_soiltemp2',mem_select,file_select  &
                   ,rvara=rs_soiltemp2)
   CALL vmissr (rs_soiltemp2,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_snowmass',mem_select,file_select   &
                   ,rvara=rs_snowmass)
   CALL vmissr (rs_snowmass,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sfc_snowdepth',mem_select,file_select  &
                   ,rvara=rs_snowdepth)
   CALL vmissr (rs_snowdepth,npts,1e30,-9998.)

   print 201,' *****  Isentropic file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of isentropic levels - ',I5,/  &
        ,' *',7X,' Isentropic levels (K)       - '/,(32X,8I5))
   print '(a)',' **********************************************'

endif

if(inout == 'OUT') then

   ! scalar values
   CALL shdf5_set_hs_select (1,'W',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_year',mem_select,file_select  &
                   ,file_chunks,ivars=iyear)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_month',mem_select,file_select &
                   ,file_chunks,ivars=imonth)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_date',mem_select,file_select  &
                   ,file_chunks,ivars=idate)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_hour',mem_select,file_select  &
                   ,file_chunks,ivars=ihour)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_nx',mem_select,file_select    &
                   ,file_chunks,ivars=n1)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_ny',mem_select,file_select    &
                   ,file_chunks,ivars=n2)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_nisn',mem_select,file_select  &
                   ,file_chunks,ivars=nisn)

   ! Vector, length is nisn 
   ! To indicate to shdf5_set_h5_select that the variable is a vector of
   ! length nisn, use negative nisn.
   CALL shdf5_set_hs_select (-nisn,'W',igrid,mem_select,file_select &
                            ,file_chunks)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_levth',mem_select,file_select &
                   ,file_chunks,ivara=levth)

   npts=n1*n2*nisn
   ! ISAN 3D isentropic variables
   CALL shdf5_set_hs_select (8,'W',igrid,mem_select,file_select,file_chunks)
   CALL vmissw (pi_u,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_u',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (pi_v,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_v',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (pi_p,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_p',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (pi_s,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_s',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (pi_r,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'isen_r',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)

   npts=n1*n2
   ! ISAN 2D variables
   CALL shdf5_set_hs_select (2,'W',igrid,mem_select,file_select,file_chunks)
   CALL vmissw (rs_u,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_u',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_v,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_v',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_p,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_p',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_t,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_t',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_r,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_r',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_s,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_s',mem_select,file_select    &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_top,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_topo',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_qual,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_qual',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   
   CALL vmissw (rs_soilmoist1,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_soilmoist1',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_soilmoist2,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_soilmoist2',mem_select,file_select &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_soiltemp1,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_soiltemp1',mem_select,file_select  &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_soiltemp2,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_soiltemp2',mem_select,file_select  &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_snowmass,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_snowmass',mem_select,file_select   &
                   ,file_chunks,rvara=pi_scra)
   CALL vmissw (rs_snowdepth,npts,pi_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sfc_snowdepth',mem_select,file_select  &
                   ,file_chunks,rvara=pi_scra)

   print 201,' *****  Isentropic file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nisn  &
        ,(levth(l),l=1,nisn)

   print 303,igridfl,gobsep,gobrad
   303 format(/,  &
         ' Grid flag (IGRIDFL)               -',I4,/  &
        ,' Grid-obs separation in degrees    -',F5.2,/  &
        ,' Grid-obs radius influence degrees -',F5.2)

endif

return
END SUBROUTINE isenio

!##############################################################################
Subroutine sigzio (h5_fid,iphdf5,inout,n1,n2,igrid)

use isan_coms
use hdf5_utils

implicit none

integer :: n1,n2,igrid,npts,l,ninn,nx3,ny3
integer*8 :: h5_fid
integer :: iphdf5
character(len=*) :: inout
type (hdf5_select_type) :: mem_select,file_select
integer, dimension(HDF5_MAX_DIMS) :: file_chunks

if(inout == 'IN') then

   ! scalar values
   CALL shdf5_set_hs_select (1,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_year',mem_select,file_select,ivars=iyy)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_month',mem_select,file_select,ivars=imm)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_date',mem_select,file_select,ivars=idd)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_hour',mem_select,file_select,ivars=ihh)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_nx',mem_select,file_select,ivars=nx3)
   CALL shdf5_irec (h5_fid,iphdf5,'isen_ny',mem_select,file_select,ivars=ny3)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_nsigz',mem_select,file_select,ivars=ninn)

   ! Vector, length is ninn which was just read in with the prior shdf5_irec call
   ! To indicate to shdf5_set_h5_select that the variable is a vector of
   ! length ninn, use negative ninn.
   CALL shdf5_set_hs_select (-ninn,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_sigz',mem_select,file_select,rvara=sigz)

   if(nx3.ne.n1.or.ny3.ne.n2.or.ninn.ne.nsigz)then
      print*,'Sigma-z grid dimensions do not match'
      print*,'   input data on read !'
      print*,' File  dimensions - ',nx3,ny3,ninn
      print*,' Input dimensions - ',n1,n2,nsigz
      stop 'iO3-2'
   endif

   npts=n1*n2*nsigz
   ! ISAN 3D sigma-z variables
   CALL shdf5_set_hs_select (9,'R',igrid,mem_select,file_select,file_chunks)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_u',mem_select,file_select,rvara=ps_u)
   CALL vmissr (ps_u,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_v',mem_select,file_select,rvara=ps_v)
   CALL vmissr (ps_v,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_p',mem_select,file_select,rvara=ps_p)
   CALL vmissr (ps_p,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_t',mem_select,file_select,rvara=ps_t)
   CALL vmissr (ps_t,npts,1e30,-9998.)
   CALL shdf5_irec (h5_fid,iphdf5,'sigz_r',mem_select,file_select,rvara=ps_r)
   CALL vmissr (ps_r,npts,1e30,-9998.)

   print 201,' *****  Sigma-z file input *****************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz  &
        ,(sigz(l),l=1,nsigz)
   201 format(//,a,//  &
        ,' *',7X,' Date (year,month,day,hour)  - ',4I5,/  &
        ,' *',7X,' Number of X,Y points        - ',2I5,/  &
        ,' *',7X,' Number of sigma-z levels    - ',I5,/  &
        ,' *',7X,' Sigma-z levels (m)          - '/,(32X,7F8.1))
   print '(a)',' **********************************************'

endif

if(inout == 'OUT') then
   
   ! scalar values
   CALL shdf5_set_hs_select (1,'W',igrid,mem_select,file_select &
                            ,file_chunks)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_nsigz',mem_select,file_select &
                   ,file_chunks,ivars=nsigz)

   ! Vector, length is nsigz.
   ! To indicate to shdf5_set_h5_select that the variable is a vector of
   ! length nsigz, use negative nsigz.
   CALL shdf5_set_hs_select (-nsigz,'R',igrid,mem_select,file_select &
                            ,file_chunks)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_sigz',mem_select,file_select &
                   ,file_chunks,rvara=sigz)

   npts=n1*n2*nsigz
   ! ISAN 3D sigma-z variables
   CALL shdf5_set_hs_select (9,'W',igrid,mem_select,file_select &
                            ,file_chunks)
   CALL vmissw (ps_u,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_u',mem_select,file_select &
                   ,file_chunks,rvara=ps_scra)
   CALL vmissw (ps_v,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_v',mem_select,file_select &
                   ,file_chunks,rvara=ps_scra)
   CALL vmissw (ps_p,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_p',mem_select,file_select &
                   ,file_chunks,rvara=ps_scra)
   CALL vmissw (ps_t,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_t',mem_select,file_select &
                   ,file_chunks,rvara=ps_scra)
   CALL vmissw (ps_r,npts,ps_scra,1E30,-9999.)
   CALL shdf5_orec (h5_fid,iphdf5,'sigz_r',mem_select,file_select &
                   ,file_chunks,rvara=ps_scra)

   print 201,' *****  Sigma-z file written *************'  &
        ,iyear,imonth,idate,ihour,n1,n2,nsigz   &
        ,(sigz(l),l=1,nsigz)

endif

return
END SUBROUTINE sigzio

!##############################################################################
Subroutine vmissw (af,n,as,fm,fx)

implicit none

integer :: n
real :: af(*),as(*),fm,fx
integer :: i

do i=1,n
   as(i)=af(i)
   if(af(i).ge.fm) as(i)=fx
enddo

return
END SUBROUTINE vmissw

!##############################################################################
Subroutine vmissr (af,n,fm,fx)

implicit none

integer :: n
real :: af(*),fm,fx
integer :: i

do i=1,n
   if(af(i).le.fx) af(i)=fm
enddo

return
END SUBROUTINE vmissr

