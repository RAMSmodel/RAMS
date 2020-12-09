!##############################################################################
Subroutine prtopt ()

use mem_grid
use mem_scratch
use ref_sounding
use rconstants
use mem_turb
use ref_sounding

implicit none

integer :: k

!PRINT INITIAL SOUNDING INPUT ONLY IF DOING HORIZONTALLY HOMOGENEOUS
!SIMULATION AT TIME=0 SINCE THE SOUNDING DATA IS NOT CONVERTED TO
!UNITS IN THE LIST BELOW UNLESS ROUTINE "INITHH" IS USED. THAT ONLY
!OCCURS IF "INITIAL=1".
IF(initial == 1)THEN
  CALL mrsl (nsndg,ps(1),ts(1),vctr5(1))
  do k=1,nsndg
     vctr1(k) = 100. * rts(k) / vctr5(k)
  enddo
  WRITE(*,41)
41      FORMAT(/,'------------------------------SOUNDING INPUT-------'  &
   ,'---------------------------',//,7X,'PS',9X,'HS',7X,'TS',6X  &
   ,'THDS',6X,'US',7X,'VS',7X,'RTS',5X,'REL HUM',/,6X  &
   ,'(Pa)',7X,'(m)',6X,'(K)',6X,'(K)',6X,'(m/s)',4X,'(m/s)'  &
   ,3X,'(g/kg)',5X,'(%)',/)
  WRITE(*,42)(PS(K),HS(K),TS(K),THDS(K),US(K),VS(K),RTS(K)*1000.  &
             ,VCTR1(K),K=1,NSndg)
42      FORMAT(1X,F11.1,F10.1,2F9.2,2F9.2,F10.5,F9.2)
ENDIF

!PRINT REFERENCE STATE FOR ANY TYPE OF SIMULATION
DO K=1,NNZP(1)
  VCTR1(K)=P00*(PI01DN(K,1)/CP)**CPOR
  VCTR2(K)=TH01DN(K,1)/(1.+.61*RT01DN(K,1))
ENDDO
WRITE(*,310)IREF,JREF,TOPREF,(ZTN(K,1),U01DN(K,1),V01DN(K,1)  &
  ,DN01DN(K,1),PI01DN(K,1),VCTR1(K),TH01DN(K,1),VCTR2(K)  &
  ,RT01DN(K,1),K=1,NNZP(1))
310   FORMAT(/,'--------REFERENCE STATE at I,J=(',I4,',',I4  &
      ,')   SFC ELEV (M)= ',F6.1,'-------------'  &
 ,//,4X,'Z',6X,'U01D',4X,'V01D',4X,'DN01D',4X  &
 ,'PI01D',4X,'PRESS',4X,'TH01D',4X,'THD',6X,'RT01D'  &
 ,/,3X,'(m)',5X,'(m/s)',3X,'(m/s)',2X,'(kg/m3)',2X  &
 ,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X,'(K)',5X,'(kg/kg)'  &
 ,//,(1X,F7.1,2F8.2,F8.3,F10.2,F10.1,2F8.2,F10.5))

return
END SUBROUTINE prtopt
