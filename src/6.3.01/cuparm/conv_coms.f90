!##############################################################################
Module conv_coms

implicit none

integer, parameter :: nkp=100

integer ::      kcon,klcl,klfc,ketl,kct,igo,kmt  &
               ,icprtfl,icpltfl
               
real    ::      zmid,cdzmin,dzlow,dzhigh,plcl,tlcl,dzlcl,zlcl,garea  &
               ,wconmin,contim,preff,envshr,supply,cptime,cprecip

real, dimension(nkp) :: ucon,vcon,wcon,thtcon ,rvcon,prcon,picon,tmpcon  &
               ,dncon,zcon,zzcon  &
               ,upe,vpe,wpe,ze,te,the,pe,rte,pke,rhoe,thve,zc  &
               ,rve,thee,qvct1,qvct2,qvct3,qvct4  &
               ,vheat,vmois,vmdry,frcon,ftcon,tcon,rcon &
               ,theu,rsu,thu,tu,thd,wtd,thcon,rtcon

END MODULE conv_coms
