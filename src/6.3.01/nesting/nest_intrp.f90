!##############################################################################
Subroutine nstbtnd (m1,m2,m3,ia,iz,ja,jz,ibcon  &
     ,scp,sct,bx,by,bz,vnam,tymeinv,nstbot,nsttop,jdim)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,nstbot,nsttop,jdim,i,j,k  &
   ,nzfm,nxfm,nyfm,incia,inciz,incja,incjz
real :: tymeinv
real, dimension(m1,m2,m3) :: scp,sct
real, dimension(m1,m3,2) :: bx
real, dimension(m1,m2,2) :: by
real, dimension(m2,m3,2) :: bz
character(len=*) :: vnam

nzfm = m1
nxfm = iz + 1
nyfm = jz + 1
if (vnam .eq. 'w') nzfm =nzfm - 1
if (vnam .eq. 'u') nxfm =nxfm - 1
if (vnam .eq. 'v') nyfm =nyfm - 1
if (jdim .eq. 0)  nyfm = 1

incia = 0
inciz = 0
if (iand(ibcon,1) .ne. 0) incia = 1
if (iand(ibcon,2) .ne. 0 .and. vnam .ne. 'u') inciz = 1
incja = 0
incjz = 0
if (iand(ibcon,4) .ne. 0) incja = jdim
if (iand(ibcon,8) .ne. 0 .and. vnam .ne. 'v') incjz = jdim

if (nstbot .eq. 0) then
  do j = ja,jz
    do i = ia,iz
      sct(1,i,j) = (bz(i,j,1) - scp(1,i,j)) * tymeinv
    enddo
  enddo
endif
if (nsttop .eq. 0) then
  do j = ja,jz
    do i = ia,iz
      sct(nzfm,i,j) = (bz(i,j,2) - scp(nzfm,i,j)) * tymeinv
    enddo
  enddo
endif

if (iand(ibcon,1) .ne. 0) then
   do j = ja-incja,jz+incjz
      do k = 1,m1
         sct(k,1,j) = (bx(k,j,1) - scp(k,1,j)) * tymeinv
      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = ja-incja,jz+incjz
      do k = 1,m1
         sct(k,nxfm,j) = (bx(k,j,2) - scp(k,nxfm,j)) * tymeinv
      enddo
   enddo
endif

if (jdim .eq. 1) then
   if (iand(ibcon,4) .ne. 0) then
      do i = ia-incia,iz+inciz
         do k = 1,m1
            sct(k,i,1) = (by(k,i,1) - scp(k,i,1)) * tymeinv
         enddo
      enddo
   endif
   if (iand(ibcon,8) .ne. 0) then
      do i = ia-incia,iz+inciz
         do k = 1,m1
            sct(k,i,nyfm) = (by(k,i,2) - scp(k,i,nyfm)) * tymeinv
         enddo
      enddo
   endif
endif

return
END SUBROUTINE nstbtnd

!##############################################################################
Subroutine cofnest ()

use mem_grid

implicit none

integer :: nf,nc,nrat,if,jf,kf,kc
real :: alpha,et,ev
real, dimension(:), allocatable :: vctr1, vctr2, vctr3

do nf = 2,ngrids
   nc = nxtnest(nf)

   allocate(vctr1(nnzp(nc)))
   allocate(vctr2(nnzp(nc)))
   allocate(vctr3(nnzp(nc)))

   if (nc .eq. 0) go to 50
   nrat = nstratx(nf)
   alpha = ((1. / float(nrat)) ** 2 - 1.) / 24.
   do if = 1,nnxp(nf)
      et = -.5 + float(2 * mod(if+nrat-2,nrat) + 1) / (2.0 * float(nrat))
      ev = -.5 + float(mod(if+nrat-2,nrat) + 1) / float(nrat)

      ei1(if,nf) = et * (et - 1.) / 2. + alpha
      ei2(if,nf) = (1. - et * et) - 2. * alpha
      ei3(if,nf) = et * (et + 1.) / 2. + alpha
      ei4(if,nf) = (ev * ev - 0.25) * (1.5 - ev) / 6.0
      ei5(if,nf) = (0.5 - ev) * (2.25 - ev * ev) * 0.5
      ei6(if,nf) = (0.5 + ev) * (2.25 - ev * ev) * 0.5
      ei7(if,nf) = (ev * ev - 0.25) * (1.5 + ev) / 6.0
   enddo

   if (jdim .eq. 1) then
      nrat = nstraty(nf)
      alpha = ((1. / float(nrat)) ** 2 - 1.) / 24.
      do jf = 1,nnyp(nf)
         et = -.5 + float(2 * mod(jf+nrat-2,nrat) + 1) / (2.0 * float(nrat))
         ev = -.5 + float(mod(jf+nrat-2,nrat) + 1) / float(nrat)

         ej1(jf,nf) = et * (et - 1.) / 2. + alpha
         ej2(jf,nf) = (1. - et * et) - 2. * alpha
         ej3(jf,nf) = et * (et + 1.) / 2. + alpha
         ej4(jf,nf) = (ev * ev - 0.25) * (1.5 - ev) / 6.0
         ej5(jf,nf) = (0.5 - ev) * (2.25 - ev * ev) * 0.5
         ej6(jf,nf) = (0.5 + ev) * (2.25 - ev * ev) * 0.5
         ej7(jf,nf) = (ev * ev - 0.25) * (1.5 + ev) / 6.0
      enddo
   endif

   do kc = 1,nnzp(nc)
      vctr1(kc) = 0.
      vctr2(kc) = 0.
      vctr3(kc) = 0.
   enddo

   do kf = 1,nnzp(nf)
      kc = kpm(kf,nf)
      ek1(kf,nf) = (ztn(kf,nf) - ztn(kc,nc)) * (ztn(kf,nf) - ztn(kc+1,nc))  &
          / ((ztn(kc-1,nc) - ztn(kc,nc)) * (ztn(kc-1,nc) - ztn(kc+1,nc)))
      ek2(kf,nf) = (ztn(kf,nf) - ztn(kc-1,nc)) * (ztn(kf,nf) - ztn(kc+1,nc  &
          )) / ((ztn(kc,nc) - ztn(kc-1,nc)) * (ztn(kc,nc) - ztn(kc+1,nc)))
      ek3(kf,nf) = (ztn(kf,nf) - ztn(kc-1,nc)) * (ztn(kf,nf) - ztn(kc,nc))  &
          / ((ztn(kc+1,nc) - ztn(kc-1,nc)) * (ztn(kc+1,nc) - ztn(kc,nc)))
      if (kf .ge. 2 .and. kf .le. nnzp(nf) - 1) then
         vctr1(kc) = vctr1(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek1(kf,nf)
         vctr2(kc) = vctr2(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek2(kf,nf)
         vctr3(kc) = vctr3(kc) + (zmn(kf,nf) - zmn(kf-1,nf)) * ek3(kf,nf)
      endif
   enddo
!c         if(kpm(2,nf).gt.2)then
!c            vctr1(1)=vctr1(1+nrz(kpm(2,nf),nc)
!c            vctr2(1)=vctr2(1)
!c            vctr3(1)=vctr3(1)
!c         endif
!c         if(kpm(nnzp(nf)-1,nf).lt.nnzp(nc)-1)then
!c            vctr1(nnzp(nf))=vctr1(nnzp(nf)-1)
!c            vctr2(nnzp(nf))=vctr2(nnzp(nf)-1)
!c            vctr3(nnzp(nf))=vctr3(nnzp(nf)-1)
!c         endif
   do kf = 2,nnzp(nf) - 1
      kc = kpm(kf,nf)
      ek1(kf,nf) = ek1(kf,nf) - dztn(kc,nc) * vctr1(kc)
      ek2(kf,nf) = ek2(kf,nf) - dztn(kc,nc) * vctr2(kc) + 1.
      ek3(kf,nf) = ek3(kf,nf) - dztn(kc,nc) * vctr3(kc)
   enddo
!cccccccccccccccc
!cc         ek1(1,nf)=ek1(2,nf)
!cc         ek2(1,nf)=ek2(2,nf)
!cc         ek3(1,nf)=ek3(2,nf)
!cccccccccccccccccc

   do kf = 1,nnzp(nf)-1
      kc = kpm(kf+1,nf)
      if (kf .eq. 1 .or. kf .eq. nnzp(nf)-1 .or. kc .ne. kpm(kf,nf)) then
         ek4(kf,nf) = 0.
         ek5(kf,nf) = 1.
         ek6(kf,nf) = 0.
         ek7(kf,nf) = 0.
      else
         ek4(kf,nf) = ek4(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
              * (-ek1(kf,nf)) * dztn(kc-1,nc)
         ek5(kf,nf) = ek5(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
              * (ek1(kf,nf) * dztn(kc-1,nc) - ek2(kf,nf) * dztn(kc,nc))
         ek6(kf,nf) = ek6(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
              * (ek2(kf,nf) * dztn(kc,nc) - ek3(kf,nf) * dztn(kc+1,nc))
         ek7(kf,nf) = ek7(kf-1,nf) + (zmn(kf,nf) - zmn(kf-1,nf))  &
              * ek3(kf,nf) * dztn(kc+1,nc)
      endif
   enddo
   do kf = 1,nnzp(nf)
      fbcf(kf,nf,1) = dztn(kpm(kf,nf),nc) / (float(nstraty(nf)) * dztn(kf,nf))
      fbcf(kf,nf,2) = dztn(kpm(kf,nf),nc) / (float(nstratx(nf)) * dztn(kf,nf))
      fbcf(kf,nf,3) = 1. / float(nstratx(nf) * nstraty(nf))
      fbcf(kf,nf,4) = dztn(kpm(kf,nf),nc)  &
                    / (float(nstratx(nf) * nstraty(nf)) * dztn(kf,nf))
   enddo
50      continue

   deallocate(vctr1)
   deallocate(vctr2)
   deallocate(vctr3)
enddo

return
END SUBROUTINE cofnest

!##############################################################################
Subroutine fill_density_tiles (nzc,nxc,nyc,dn01dc,dn0c,dn0uc,dn0vc, &
                              nzct,nxct,nyct,dn01dct,dn0ct,dn0uct,dn0vct, &
                              ic0,jc0,ict0,jct0,ifm,icm)

! This routine will call fill_coarse_tile() for the 1D and 3D density
! variables.

  use node_mod

  implicit none

  integer :: nzc, nxc, nyc, nzct, nxct, nyct
  integer :: ic0, jc0, ict0, jct0, ifm, icm
  real, dimension(nzc) :: dn01dc
  real, dimension(nzct) :: dn01dct
  real, dimension(nzc,nxc,nyc) :: dn0c, dn0uc, dn0vc
  real, dimension(nzct,nxct,nyct) :: dn0ct, dn0uct, dn0vct

  ! 1D Density
  CALL fill_coarse_tile (nzc,1,1,dn01dc,nzct,1,1,dn01dct,0,0,0,0,ifm,icm,1)

  ! 3D Density
  CALL fill_coarse_tile (nzc,nxc,nyc,dn0c,nzct,nxct,nyct,dn0ct &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,dn0uc,nzct,nxct,nyct,dn0uct &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)
  CALL fill_coarse_tile (nzc,nxc,nyc,dn0vc,nzct,nxct,nyct,dn0vct &
                        ,ic0,jc0,ict0,jct0,ifm,icm,3)

  return
END SUBROUTINE fill_density_tiles

!##############################################################################
Subroutine interp_refs_1d (nzc,dn01dc,th01dc,u01dc,v01dc,rt01dc,pi01dc,nzct &
         ,dn01dct,nzf,nxf,nyf,dn01df,th01df,u01df,v01df,rt01df,pi01df,toptf &
         ,toptf_int,ifm,icm,ctipm,ctjpm,ctkpm,ibconf,ifg_mode)

! This routine will interpolate the 1D reference state from the parent grid.

  use node_mod, only:IFG_HH_INIT,IFG_HIST_INIT,IFG_VF_INIT
  use rconstants

  implicit none

  integer :: nzc, nzct, nzf, nxf, nyf, ifm, icm, ibconf, ifg_mode
  real, dimension(nzc) :: dn01dc, th01dc, u01dc, v01dc, rt01dc, pi01dc
  real, dimension(nzct) :: dn01dct
  real, dimension(nzf) :: dn01df, th01df, u01df, v01df, rt01df, pi01df, toptf, toptf_int
  integer, dimension(nxf) :: ctipm
  integer, dimension(nyf) :: ctjpm
  integer, dimension(nzf) :: ctkpm

  integer :: k
  real :: c1, c2

  ! for calculation of pi (exner)
  c1 = rgas / (cp - rgas)
  c2 = cp * (rgas / p00) ** c1

  ! Interpolate dn01d, th01d, u01d, v01d, rt01d and use dn01d, th01d to calculate pi
  CALL interp_var (nzc,1,1,dn01dc,nzct,1,1,dn01dct,nzf,1,1,dn01df,dn01df,toptf,toptf_int, &
                  ifm,icm,ctipm,ctjpm,ctkpm,0,0,0,0,0,0,ibconf,'t',0,0,1)

  CALL interp_var (nzc,1,1,th01dc,nzct,1,1,dn01dct,nzf,1,1,th01df,dn01df,toptf,toptf_int, &
                  ifm,icm,ctipm,ctjpm,ctkpm,0,0,0,0,0,0,ibconf,'t',1,0,1)

  if ((ifg_mode .eq. IFG_HH_INIT) .or. (ifg_mode .eq. IFG_HIST_INIT) & 
                                  .or. (ifg_mode .eq. IFG_VF_INIT)) then
    CALL interp_var (nzc,1,1,u01dc,nzct,1,1,dn01dct,nzf,1,1,u01df,dn01df,toptf,toptf_int, &
                    ifm,icm,ctipm,ctjpm,ctkpm,0,0,0,0,0,0,ibconf,'t',1,0,1)

    CALL interp_var (nzc,1,1,v01dc,nzct,1,1,dn01dct,nzf,1,1,v01df,dn01df,toptf,toptf_int, &
                    ifm,icm,ctipm,ctjpm,ctkpm,0,0,0,0,0,0,ibconf,'t',1,0,1)
  endif

  CALL interp_var (nzc,1,1,rt01dc,nzct,1,1,dn01dct,nzf,1,1,rt01df,dn01df,toptf,toptf_int, &
                  ifm,icm,ctipm,ctjpm,ctkpm,0,0,0,0,0,0,ibconf,'t',1,0,1)

  do k = 1,nzf
    pi01df(k) = c2 * (dn01df(k) * th01df(k)) ** c1
  enddo

  return
END SUBROUTINE interp_refs_1d

!##############################################################################
Subroutine interp_refs_3d (nzc,nxc,nyc,dn0c,th0c,nzct,nxct,nyct,dn0ct, &
            nzf,nxf,nyf,dn0f,dn0uf,dn0vf,th0f,pi0f,pf,toptf,toptf_int, &
            ifm,icm,ctipm,ctjpm,ctkpm,ic0,jc0,ict0,jct0,if0,jf0,ibconf,ifg_mode)

! This routine will interpolate the 3D reference state from the parent grid.

  use node_mod, only:nmachs,IFG_MKVF_INIT,IFG_MKHF_INIT
  use rconstants

  implicit none

  integer :: nzc,nxc,nyc,nzct,nxct,nyct,nzf,nxf,nyf
  integer :: ifm,icm,ic0,jc0,ict0,jct0,if0,jf0,ibconf,ifg_mode
  real, dimension(nzc,nxc,nyc) :: dn0c,th0c
  real, dimension(nzct,nxct,nyct) :: dn0ct
  real, dimension(nzf,nxf,nyf) :: dn0f,dn0uf,dn0vf,th0f,pi0f,pf,toptf,toptf_int
  integer, dimension(nzf) :: ctipm,ctjpm,ctkpm

  integer :: i,j,k
  real :: c1,c2

  ! for calculation of pi (exner)
  c1 = rgas / (cp - rgas)
  c2 = cp * (rgas / p00) ** c1

  ! Interpolate dn0, th0, and calculate pi0
  !   1. Do dn0 and th0 first without terrain
  !   2. Use dn0 and th0 (interpolated without terraine) to calculate pi0
  !   3. Apply terrain interpolation to th0 and pi0 (leave dn0 alone for now)

  CALL interp_var (nzc,nxc,nyc,dn0c,nzct,nxct,nyct,dn0ct, &
         nzf,nxf,nyf,dn0f,dn0f,toptf,toptf_int, &
         ifm,icm,ctipm,ctjpm,ctkpm,ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',0,0,3)

  CALL interp_var (nzc,nxc,nyc,th0c,nzct,nxct,nyct,dn0ct, &
         nzf,nxf,nyf,th0f,dn0f,toptf,toptf_int, &
         ifm,icm,ctipm,ctjpm,ctkpm,ic0,jc0,ict0,jct0,if0,jf0,ibconf,'t',1,0,3)

  do j = 1,nyf
    do i = 1,nxf
      do k = 1,nzf
        pi0f(k,i,j) = c2 * (dn0f(k,i,j) * th0f(k,i,j)) ** c1
      enddo
    enddo
  enddo
  if (nmachs .gt. 1) then
    CALL update_lbc_var (nzf,nxf,nyf,ifm,pi0f,3)
  endif

  CALL interp_topo (nzf,nxf,nyf,th0f,toptf_int,toptf,ifm,'t',ibconf,3)
  CALL interp_topo (nzf,nxf,nyf,pi0f,toptf_int,toptf,ifm,'t',ibconf,3)

  ! Average dn0 (grid center) to the horizontal grid cell edges.
  CALL calc_dn0uv (nzf,nxf,nyf,dn0f,dn0uf,dn0vf,ibconf,ifm)

  ! Form the perturbation pressure if doing MAKEVFILE
  if ((ifg_mode .eq. IFG_MKVF_INIT) .or. (ifg_mode .eq. IFG_MKHF_INIT)) then
    if (ifg_mode .eq. IFG_MKVF_INIT) then
      do j = 1,nyf
        do i = 1,nxf
          do k = 1,nzf
            pf(k,i,j) =  pf(k,i,j) - pi0f(k,i,j)
          enddo
        enddo
      enddo
    elseif (ifg_mode .eq. IFG_MKHF_INIT) then
      do j = 1,nyf
        do i = 1,nxf
          do k = 1,nzf
            pf(k,i,j) =  pf(k,i,j) + pi0f(k,i,j)
          enddo
        enddo
      enddo
    endif
    if (nmachs .gt. 1) then
      CALL update_lbc_var (nzf,nxf,nyf,ifm,pf,3)
    endif
  endif

  return
END SUBROUTINE interp_refs_3d

!##############################################################################
Subroutine interp_var (nzc,nxc,nyc,varc,nzct,nxct,nyct,densc,nzf,nxf,nyf,varf &
   ,densf,toptf,toptf_int,ifm,icm,ctipm,ctjpm,ctkpm,ic0,jc0,ict0,jct0,if0,jf0 &
   ,ibcon,vloc,idwflg,irtgflg,dimtype)

! This routine will interpolate one variable. This routine will handle the
! application of density weighting. It can also deal with 1D, 2D and 3D fields.
! However, the arrays are always declared as 3D (z,x,y). Just set the corresponding
! dimension size to 1 if calling with a field of lesser rank than 3. Eg, if have
! 2D horizontal field (x,y), set nz to 1; if have a 2D vertical field (z,y), set
! nx to 1; etc.
!
! The argument vloc represents where the quantity is located on the Arakawa-C grid
!   'u' - west and east edges
!   'v' - north and south edges
!   'w' - top and bottom edges
!   't' - center (this can actually be anything but 'u', 'v', or 'w')
!
! For the density weighting, make sure corresponding density variable is used.
! Eg, for scalars, use dn0; for u, use dn0u; for v, use dn0v; etc.

  use node_mod, only: nmachs

  implicit none

  integer :: nzc, nxc, nyc, nzct, nxct, nyct, nzf, nxf, nyf, ifm, icm
  integer :: ic0, jc0, ict0, jct0, if0, jf0, ibcon, idwflg, irtgflg, dimtype
  character :: vloc
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzct,nxct,nyct) :: densc
  real, dimension(nzf,nxf,nyf) :: varf, densf
  real, dimension(nxf,nyf) :: toptf, toptf_int
  integer, dimension(nxf) :: ctipm
  integer, dimension(nyf) :: ctjpm
  integer, dimension(nzf) :: ctkpm

  integer :: avgflg
  ! make temporary vars allocatable so that they don't take up space on the stack
  real, dimension(:,:,:), allocatable :: varc_tile, tmp_varc, tmp_varf

  ! Extract the coarse tile which is the rectangular piece 
  allocate(varc_tile(nzct,nxct,nyct))
  CALL fill_coarse_tile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile,ic0,jc0 &
                        ,ict0,jct0,ifm,icm,dimtype)

  if (idwflg .eq. 1) then
    allocate(tmp_varc(nzct,nxct,nyct))
    allocate(tmp_varf(nzf,nxf,nyf))

    ! For w, average the density, located in grid centers, to the grid edges.
    if (vloc .eq. 'w') then
      avgflg = 1
    else
      avgflg = 0
    endif

    ! apply density weighting
    ! this is a tile which needs all of its i,j locations worked on. Trick 
    ! denswt into applying weights everywhere by telling it that the varible 
    ! location is in the grid center (vloc = 't'), and
    ! all four sides are domain boundaries (ibcon = 15).
    CALL denswt (nzct,nxct,nyct,varc_tile,densc,tmp_varc,avgflg,1,'t',15)

    ! interpolate
    CALL interp_field (nzct,nxct,nyct,tmp_varc,nzf,nxf,nyf,tmp_varf,ctipm &
                      ,ctjpm,ctkpm,if0,jf0,ifm,vloc,ibcon)

    ! remove density weighting
    CALL denswt (nzf,nxf,nyf,tmp_varf,densf,varf,avgflg,0,vloc,ibcon)

    deallocate(tmp_varc)
    deallocate(tmp_varf)
  else
    ! just interpolate
    CALL interp_field (nzct,nxct,nyct,varc_tile,nzf,nxf,nyf,varf,ctipm,ctjpm &
                      ,ctkpm,if0,jf0,ifm,vloc,ibcon)
  endif

  ! Vertical interpolation according to topography
  ! If in parallel mode, interp_topo will update the subdomain boundaries. 
  ! If irtgflg is not 1, then we need to update the boundaries in this routine.
  if (irtgflg .eq. 1) then
    CALL interp_topo (nzf,nxf,nyf,varf,toptf_int,toptf,ifm,vloc,ibcon,dimtype)
  else
    ! Skip the update if this is a 1D (dimtype == 1) variable since there is no
    ! need (no lateral boundary overlap).
    if ((nmachs .gt. 1) .and. (dimtype .ne. 1)) then
      CALL update_lbc_var (nzf,nxf,nyf,ifm,varf,dimtype)
    endif
  endif

  deallocate(varc_tile)

  return
END SUBROUTINE interp_var

!##############################################################################
Subroutine interp_varsfc (nzc,nxc,nyc,varc,nzct,nxct,nyct,nzf,nxf,nyf,varf &
   ,ifm,icm,ctipm,ctjpm,ctkpm,ic0,jc0,ict0,jct0,if0,jf0 &
   ,ibcon,vloc,dimtype)

! This routine is the counterpart to interp_var but for use with runtypes:
!  IFG_MKVF_INIT     initialization for MAKEVFILE
!  IFG_MKHF_INIT     initialization for MAKEHFILE
!  IFG_TOPT1_INIT    initialization for MAKESFC, topography
!  IFG_TOPT2_INIT    initialization for MAKESFC, topography
!  IFG_SST_INIT      initialization for MAKESFC, sea surface temperature
!
! These runtypes associated with make surface and make varfile do not require
! full memory allocation, and so these are streamlined to save memory since
! these runtypes can only be done on a single processor (nmachs=1) by default.
!
! The argument vloc represents where the quantity is located on the Arakawa-C grid
!   'u' - west and east edges
!   'v' - north and south edges
!   'w' - top and bottom edges
!   't' - center (this can actually be anything but 'u', 'v', or 'w')

  use node_mod, only: nmachs

  implicit none

  integer :: nzc, nxc, nyc, nzct, nxct, nyct, nzf, nxf, nyf, ifm, icm
  integer :: ic0, jc0, ict0, jct0, if0, jf0, ibcon, dimtype
  character :: vloc
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzf,nxf,nyf) :: varf
  integer, dimension(nxf) :: ctipm
  integer, dimension(nyf) :: ctjpm
  integer, dimension(nzf) :: ctkpm

  integer :: avgflg
  ! make temporary vars allocatable so that they don't take up space on the stack
  real, dimension(:,:,:), allocatable :: varc_tile

  ! Extract the coarse tile which is the rectangular piece 
  allocate(varc_tile(nzct,nxct,nyct))
  CALL fill_coarse_tile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile,ic0,jc0 &
                        ,ict0,jct0,ifm,icm,dimtype)

  ! just interpolate
  CALL interp_field (nzct,nxct,nyct,varc_tile,nzf,nxf,nyf,varf,ctipm,ctjpm &
                      ,ctkpm,if0,jf0,ifm,vloc,ibcon)

  deallocate(varc_tile)

  return
END SUBROUTINE interp_varsfc

!##############################################################################
Subroutine fill_coarse_tile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile &
                            ,ic0,jc0,ict0,jct0,ifm,icm,dimtype)

! This routine will fill in coarse variable with the appropriate tile region.
!
! For sequential mode, simply pull the tile directly out of the coarse variables.
!
! For parallel mode, go through the send/receive MPI exchange which will piece
! together the tile.
!
! Variable numbers are defined in node_mod.f90. If no match for any variable
! number, then vnum is the index into the scalar table.

  use node_mod
  use var_tables
  use mem_basic
  use mem_grid
  use mem_varinit
  use ref_sounding

  implicit none

  integer :: nzc,nxc,nyc,nzct,nxct,nyct,ic0,jc0,ict0,jct0,ifm,icm,dimtype
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzct,nxct,nyct) :: varc_tile

  if (nmachs > 1) then
    ! Parallel mode
    CALL update_ctile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile,ic0,jc0 &
                      ,ifm,icm,dimtype)
  else
    ! Sequential mode
    CALL copy_coarse_tile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile,ict0,jct0)
  endif 

  return
END SUBROUTINE fill_coarse_tile

!##############################################################################
Subroutine copy_coarse_tile (nzc,nxc,nyc,varc,nzct,nxct,nyct,varc_tile,ict0,jct0)

! This routine will pull out a tile from the main coarse variable and
! place a copy of this into ctile.

  implicit none

  integer :: nzc, nxc, nyc, nzct, nxct, nyct, ict0, jct0
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzct,nxct,nyct) :: varc_tile

  integer :: i, j, k

  do j = 1, nyct
    do i = 1, nxct
      do k = 1, nzct
        varc_tile(k,i,j) = varc(k,i+ict0,j+jct0)
      enddo
    enddo
  enddo

  return
END SUBROUTINE copy_coarse_tile

!##############################################################################
Subroutine denswt (nz,nx,ny,var,dens,vout,avgflg,idir,vloc,ibcon)

! This routine will apply density weighting as part of the interpolation
! process.
!
! The avgflg argument controls whether the density is averaged in
! the vertical
!     avgflg    action
!       1         use average of density at k and k+1 levels
!       0         use (k,i,j) value of density
!
! The idir argument controls whether the weights are being
! applied or removed.
!      idir     action
!       1        multiply by density
!       0        divide by density

  implicit none

  integer :: nz, nx, ny, avgflg, idir, ibcon
  real, dimension(nz,nx,ny) :: var, dens, vout
  character :: vloc

  integer :: i, j, k
  integer :: i1, i2, j1, j2

  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,vloc,ibcon)
  if (idir .eq. 1) then
    ! multiply by density
    if (avgflg .eq. 1) then
      ! vertically average density
      do j = j1, j2
        do i = i1, i2
          do k = 1, nz-1
            vout(k,i,j) = var(k,i,j) * 0.5 * (dens(k,i,j)+dens(k+1,i,j))
          enddo
        enddo
      enddo
    else
      ! use density as is
      do j = j1, j2
        do i = i1, i2
          do k = 1, nz
            vout(k,i,j) = var(k,i,j) * dens(k,i,j)
          enddo
        enddo
      enddo
    endif
  else
    ! divide by density
    if (avgflg .eq. 1) then
      ! vertically average density
      do j = j1, j2
        do i = i1, i2
          do k = 1, nz-1
            vout(k,i,j) = var(k,i,j) / (0.5 * (dens(k,i,j)+dens(k+1,i,j)))
          enddo
        enddo
      enddo
    else
      ! use density as is
      do j = j1, j2
        do i = i1, i2
          do k = 1, nz
            vout(k,i,j) = var(k,i,j) / dens(k,i,j)
          enddo
        enddo
      enddo
    endif
  endif

  return
END SUBROUTINE denswt

!##############################################################################
Subroutine interp_field (nzc,nxc,nyc,varc,nzf,nxf,nyf,varf,ctipm,ctjpm,ctkpm &
                        ,if0,jf0,ifm,vloc,ibcon)

! This routine will interpolate the coarse tile (varc) onto the fine tile
! (varf). It is assumed that the caller arranged the coarse tile so that
! its origin has the correct offset from the origin of the fine tile.
!
! This routine can also deal with 1D, 2D and 3D fields. However, the arrays are
! declared as 3D (z,x,y). Just set the corresponding dimension size to 1 if
! calling with a field of lesser rank than 3. Eg, if have 2D horizontal
! field (x,y), set nz to 1; if have a 2D vertical field (z,y), set nx to 1; etc.
!
! The argument vloc gives where the quantity is located on the Arakawa-C grid
!   'u' - west and east edges
!   'v' - north and south edges
!   'w' - top and bottom edges
!   't' - center (this can actually be anything but 'u', 'v', or 'w')

  use mem_grid

  implicit none

  integer :: nzc,nxc,nyc,nzf,nxf,nyf,if0,jf0,ifm,ibcon
  character :: vloc
  real, dimension(nzc,nxc,nyc) :: varc
  real, dimension(nzf,nxf,nyf) :: varf
  integer, dimension(nzf) :: ctkpm
  integer, dimension(nxf) :: ctipm
  integer, dimension(nyf) :: ctjpm

  integer :: ic, jc, kc, if, jf, kf, i1f, i2f, j1f, j2f
  real, dimension(:,:,:), allocatable :: temp_var1, temp_var2

  ! The idea here is to process one dimension at a time. First the x-dimension,
  ! then the y-dimension, and last the z-dimension. As each dimension is
  ! interpolated, that dimension goes from the coarse positions to the fine
  ! positions. There are two scratch variables allocated in order to
  ! hold the results after processing the x and y dimensions. The flow
  ! looks like:
  !   1. interpolate x-dim:      varc(nzc,nxc,nyc) --> temp_var1(nzc,nxf,nyc)
  !   2. interpolate y-dim: temp_var1(nzc,nxf,nyc) --> temp_var2(nzc,nxf,nyf)
  !   3. interpolate z-dim: temp_var2(nzc,nxf,nyf) --> varf(nzf,nxf,nyf)
  !
  ! Make temp_var1, temp_var2 allocatable so that the memory gets placed
  ! in the heap instead of the stack. This will keep the stack from growing
  ! too large.

  allocate(temp_var1(nzc,nxf,nyc))
  allocate(temp_var2(nzc,nxf,nyf))

  CALL set_horiz_limits (nxf,nyf,i1f,i2f,j1f,j2f,vloc,ibcon)
  ! interpolate x-dimension
  if (nxf .eq. 1) then
    ! x-dimension is not used, assume that nxc is also 1
    ! copy varc to temp_var1
    do jc = 1,nyc
      do kc = 1,nzc
        temp_var1(kc,1,jc) = varc(kc,1,jc)
      enddo
    enddo
  else
    if (vloc .eq. 'u') then
      ! u - interpolate to the grid cell edge
      ! leave off last x position since that is an "extra" u entry
      do jc = 1,nyc
        do if = i1f,i2f
          ic = ctipm(if)
          do kc = 1,nzc
            temp_var1(kc,if,jc)=ei4(if+if0,ifm)*varc(kc,ic-2,jc)  &
                               +ei5(if+if0,ifm)*varc(kc,ic-1,jc)  &
                               +ei6(if+if0,ifm)*varc(kc,ic  ,jc)  &
                               +ei7(if+if0,ifm)*varc(kc,ic+1,jc)
          enddo
        enddo
      enddo
    else
      ! not u - interpolate to the grid cell center
      do jc = 1,nyc
        do if = i1f,i2f
          ic = ctipm(if)
          do kc = 1,nzc
            temp_var1(kc,if,jc)=ei1(if+if0,ifm)*varc(kc,ic-1,jc)  &
                               +ei2(if+if0,ifm)*varc(kc,ic  ,jc)  &
                               +ei3(if+if0,ifm)*varc(kc,ic+1,jc)
          enddo
        enddo
      enddo
    endif
  endif

  ! interpolate y-dimension
  if (nyf .eq. 1) then
    ! y-dimension is not used, assume that nyc is also 1
    ! copy temp_var1 to temp_var2
    do if = i1f,i2f
      do kc = 1,nzc
        temp_var2(kc,if,1) = temp_var1(kc,if,1)
      enddo
    enddo
  else
    if (vloc .eq. 'v') then
      ! v - interpolate to the grid cell edge
      ! leave off last y position since that is an "extra" v entry
      do jf = j1f,j2f
        jc = ctjpm(jf)
        do if = i1f,i2f
          do kc = 1,nzc
            temp_var2(kc,if,jf)=ej4(jf+jf0,ifm)*temp_var1(kc,if,jc-2)  &
                               +ej5(jf+jf0,ifm)*temp_var1(kc,if,jc-1)  &
                               +ej6(jf+jf0,ifm)*temp_var1(kc,if,jc  )  &
                               +ej7(jf+jf0,ifm)*temp_var1(kc,if,jc+1)
  
          enddo
        enddo
      enddo
    else
      ! not v - interpolate to the grid cell center
      do jf = j1f,j2f
        jc = ctjpm(jf)
        do if = i1f,i2f
          do kc = 1,nzc
            temp_var2(kc,if,jf)=ej1(jf+jf0,ifm)*temp_var1(kc,if,jc-1)  &
                               +ej2(jf+jf0,ifm)*temp_var1(kc,if,jc  )  &
                               +ej3(jf+jf0,ifm)*temp_var1(kc,if,jc+1)
          enddo
        enddo
      enddo
    endif
  endif

  ! interpolate z-dimension
  if (nzf .eq. 1) then
    ! z-dimension is not used, assume that nzc is also 1
    ! copy temp_var2 to varf
    do jf = j1f,j2f
      do if = i1f,i2f
        varf(1,if,jf) = temp_var2(1,if,jf)
      enddo
    enddo
  else
    if (vloc .eq. 'w') then
      ! w - interpolate to the grid cell edge
      ! leave off last z position since that is an "extra" w entry
      do jf = j1f,j2f
        do if = i1f,i2f
          do kf = 1,nzf-1
            kc = ctkpm(kf+1)
            varf(kf,if,jf) = ek4(kf,ifm) * temp_var2(max(1    ,kc-2),if,jf)  &
                           + ek5(kf,ifm) * temp_var2(          kc-1 ,if,jf)  &
                           + ek6(kf,ifm) * temp_var2(          kc   ,if,jf)  &
                           + ek7(kf,ifm) * temp_var2(min(nzc-1,kc+1),if,jf)
          enddo
        enddo
      enddo
    else
      ! not w - interpolate to the grid cell center
      do jf = j1f,j2f
        do if = i1f,i2f
          do kf = 1,nzf
            kc = ctkpm(kf)
            varf(kf,if,jf) = ek1(kf,ifm) * temp_var2(kc-1,if,jf)  &
                           + ek2(kf,ifm) * temp_var2(kc  ,if,jf)  &
                           + ek3(kf,ifm) * temp_var2(kc+1,if,jf)
          enddo
        enddo
      enddo
    endif
  endif

  deallocate(temp_var1)
  deallocate(temp_var2)

  return
END SUBROUTINE interp_field

!##############################################################################
Subroutine interp_topo (nz,nx,ny,var,topt_interp,topt_assign,ngr,vloc &
                       ,ibcon,dimtype)

! This routine will perform a special vertical interpolation to cover the
! case where the terrain assigned to the fine grid is different from what
! would be the terrain interpolated from the coarse grid.

  use mem_grid, only: zmn,ztop,ztn
  use node_mod, only: nmachs

  implicit none

  integer :: nz, nx, ny, ngr, ibcon, dimtype
  real, dimension(nz,nx,ny) :: var
  real, dimension(nx,ny) :: topt_interp, topt_assign
  character :: vloc

  integer :: i, j, k, i1, i2, j1, j2
  real, dimension(:), allocatable :: hghts_interp, hghts_assign, temp_var

  allocate(hghts_interp(nz))
  allocate(hghts_assign(nz))
  allocate(temp_var(nz))

  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,vloc,ibcon)
  do j = j1,j2
   do i = i1,i2
    ! If doing 'w', interpolate to the grid cell edges; otherwise interpolate
    ! to the grid cell centers.
    if (vloc .eq. 'w') then
     do k = 1,nz
      hghts_interp(k) = zmn(k,ngr)*(1.0-topt_interp(i,j)/ztop)+topt_interp(i,j)
      hghts_assign(k) = zmn(k,ngr)*(1.0-topt_assign(i,j)/ztop)+topt_assign(i,j)
      temp_var(k)     = var(k,i,j)
     enddo
    else
     do k = 1,nz
      hghts_interp(k) = ztn(k,ngr)*(1.0-topt_interp(i,j)/ztop)+topt_interp(i,j)
      hghts_assign(k) = ztn(k,ngr)*(1.0-topt_assign(i,j)/ztop)+topt_assign(i,j)
      temp_var(k)     = var(k,i,j)
     enddo
    endif
    CALL htint (nz,temp_var,hghts_interp,nz,var(1,i,j),hghts_assign)
   enddo
  enddo

  ! Update subdomain boundaries if in parallel mode. Skip the update if this
  ! is a column variable (no lateral overlap), ie dimtype == 1.
  if ((nmachs .gt. 1) .and. (dimtype .ne. 1)) then
    CALL update_lbc_var (nz,nx,ny,ngr,var,dimtype)
  endif

  deallocate(hghts_interp)
  deallocate(hghts_assign)
  deallocate(temp_var)
  
  return
END SUBROUTINE interp_topo

!##############################################################################
Subroutine copy_nest_bounds (nz,nx,ny,var,bx,by,bz,vloc,nstbot,nsttop,ibcon)

! This routine will mimic what nstb() does, but will in addition handle
! a subdomain when running in parallel mode.

  implicit none

  integer :: nz, nx, ny, nstbot, nsttop, ibcon
  character :: vloc
  real, dimension(nz,nx,ny) :: var
  real, dimension(nz,ny,2) :: bx
  real, dimension(nz,nx,2) :: by
  real, dimension(nx,ny,2) :: bz

  integer :: i, j, k
  integer :: i1, i2, j1, j2, k1, k2

  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,vloc,ibcon)

  ! West side
  if (iand(ibcon,1) .ne. 0) then
    do j = j1, j2
      do k = 1, nz
        bx(k,j,1) = var(k,i1,j)
      enddo
    enddo
  endif
      
  ! East side
  if (iand(ibcon,2) .ne. 0) then
    do j = j1, j2
      do k = 1, nz
        bx(k,j,2) = var(k,i2,j)
      enddo
    enddo
  endif
      
  ! South side
  if (iand(ibcon,4) .ne. 0) then
    do i = i1, i2
      do k = 1, nz
        by(k,i,1) = var(k,i,j1)
      enddo
    enddo
  endif
      
  ! North side
  if (iand(ibcon,8) .ne. 0) then
    do i = i1, i2
      do k = 1, nz
        by(k,i,2) = var(k,i,j2)
      enddo
    enddo
  endif

  k1 = 1
  if (vloc .eq. 'w') then
    k2 = nz - 1
  else
    k2 = nz
  endif

  ! Bottom side
  if (nstbot .eq. 0) then
    do j = j1,j2
      do i = i1,i2
        bz(i,j,1) = var(k1,i,j)
      enddo
    enddo
  endif

  ! Top side
  if (nsttop .eq. 0) then
    do j = j1,j2
      do i = i1,i2
        bz(i,j,2) = var(k2,i,j)
      enddo
    enddo
  endif

END SUBROUTINE copy_nest_bounds

!##############################################################################
Subroutine calc_dn0uv (nz,nx,ny,denst,densu,densv,ibcon,ngr)

! This routine will calculate the density values at the u and v grid locations
! by averaging the density values at the t grid locations.

  use node_mod, only: nmachs

  implicit none

  integer :: nz, nx, ny, ibcon, ngr
  real, dimension(nz,nx,ny) :: denst, densu, densv

  integer :: i, j, k, i1, i2, j1, j2

  ! Average to u locations. 
  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,'u',ibcon)
  do j = j1,j2
    do i = i1,i2
      do k = 1,nz
        if (i .lt. nx) then
          densu(k,i,j) = .5 * (denst(k,i,j) + denst(k,i+1,j))
        else
          densu(k,i,j) = denst(k,i,j)
        endif
      enddo
    enddo
  enddo

  ! Average to v locations. 
  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,'v',ibcon)
  do j = j1,j2
    do i = i1,i2
      do k = 1,nz
        if (j .lt. ny) then
          densv(k,i,j) = .5 * (denst(k,i,j) + denst(k,i,j+1))
        else
          densv(k,i,j) = denst(k,i,j)
        endif
      enddo
    enddo
  enddo

  ! Update subdomain boundaries if in parallel mode.
  if (nmachs .gt. 1) then
    CALL update_lbc_var (nz,nx,ny,ngr,densu,3)
    CALL update_lbc_var (nz,nx,ny,ngr,densv,3)
  endif

  return
END SUBROUTINE calc_dn0uv

!##############################################################################
Subroutine init_3d_refstate (nz,nx,ny,pi0,dn0,dn0u,dn0v,th0,topt,rtgt,ibcon)

  use mem_grid, only:zt,nzp,ngrid,ztop,dzm
  use ref_sounding
  use rconstants
  use node_mod, only:nmachs,mibcon

  implicit none

  integer :: nz,nx,ny
  real, dimension(nz,nx,ny) :: pi0, dn0, dn0u, dn0v, th0
  real, dimension(nx,ny) :: topt, rtgt

  integer :: i,j,k,i1,i2,j1,j2,ibcon
  real :: c1,c2,c3
  real, dimension(:), allocatable :: zt_topo

  CALL set_horiz_limits (nx,ny,i1,i2,j1,j2,'t',ibcon)

  allocate(zt_topo(nz))

  do j=j1,j2
    do i=i1,i2

      do k = 1,nz
         zt_topo(k) = zt(k) * rtgt(i,j) + topt(i,j)
      enddo
      CALL htint (nzp,pi01dn(1,ngrid),zt,nzp,pi0(1,i,j),zt_topo)
      CALL htint (nzp,th01dn(1,ngrid),zt,nzp,th0(1,i,j),zt_topo)
      c1 = g * 2. * (1. - topt(i,j) / ztop)

      c2 = 1. - cpor
      c3 = cp ** c2
      do k = nz-1,1,-1
        pi0(k,i,j) = pi0(k+1,i,j)  &
                   + c1 / ((th0(k,i,j) + th0(k+1,i,j)) * dzm(k))
      enddo

      do k = 1,nz
         dn0(k,i,j) = (c3 * p00) / (rgas * th0(k,i,j) * pi0(k,i,j) ** c2)
      enddo

    enddo
  enddo

  ! Need to update boundaries on dn0 for the subsequent calc_dn0uv call
  if (nmachs .gt. 1) CALL update_lbc_var (nz,nx,ny,ngrid,dn0,3)
  CALL calc_dn0uv (nz,nx,ny,dn0,dn0u,dn0v,mibcon(ngrid),ngrid)

  ! Update boundaries for th0 and pi0 (calc_dn0uv already updated dn0u and dn0v)
  if (nmachs .gt. 1) then
    CALL update_lbc_var (nz,nx,ny,ngrid,th0,3)
    CALL update_lbc_var (nz,nx,ny,ngrid,pi0,3)
  endif

  deallocate(zt_topo)

  return
END SUBROUTINE init_3d_refstate

!##############################################################################
Subroutine set_horiz_limits (nx,ny,i1,i2,j1,j2,vloc,ibcon)

! This routine will consider 1D, 2D, 3D fields, the variable being worked on,
! and the domain boundaries in order to set the x and y loop limits. This is
! intended to be used by the interpolation routines.

  implicit none

  integer :: nx, ny, i1, i2, j1, j2, ibcon
  character :: vloc

  if (nx .le. 1) then
    i1 = 1
    i2 = 1
  else
    i1 = 2
    ! If subdomain contains east domain boundary
    if (iand(ibcon,1) .ne. 0) i1 = 1

    i2 = nx-1
    ! If subdomain contains west domain boundary
    if ((iand(ibcon,2) .ne. 0) .and. (vloc .ne. 'u')) i2 = nx
  endif

  if (ny .le. 1) then
    j1 = 1
    j2 = 1
  else
    j1 = 2
    ! If subdomain contains south domain boundary
    if (iand(ibcon,4) .ne. 0) j1 = 1

    j2 = ny-1
    ! If subdomain contains north domain boundary
    if ((iand(ibcon,8) .ne. 0) .and. (vloc .ne. 'v')) j2 = ny
  endif

  return
END SUBROUTINE set_horiz_limits
