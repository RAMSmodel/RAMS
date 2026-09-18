

! CVS:  $Id: bugsrad_planck.F90,v 1.1 2003/08/06 19:41:13 norm Exp $
! CVS:  $Name:  $

module bugsrad_planck

   use bugs_kinds
   implicit none

   integer, private, parameter :: MBIR = 12  !Number of IR bands
   real (kind=dbl_kind), private, dimension(6, MBIR) ::  &
      b = reshape ((/ -25.889132,     0.75038381,     -0.87074567E-02,  0.50701144E-04, -0.14856755E-06,  0.17579587E-09, &     !Band 1    Coefficients for fitted polynomial
                       25.397471,    -0.59596460,      0.53117737E-02, -0.21681758E-04,  0.36630792E-07, -0.11541419E-10, &     !Band 2    which computes the blackbody flux emission
                       57.891546,    -1.4745788,       0.14577775E-01, -0.68637478E-04,  0.14707480E-06, -0.98862337E-10, &     !Band 3    as a function of temperature and band.
                       21.837317,    -0.63194381,      0.71338812E-02, -0.38569394E-04,  0.95685257E-07, -0.76188561E-10, &     !Band 4
                      0.83155466,    -0.15281669,      0.31020500E-02, -0.23768837E-04,  0.74605666E-07, -0.67494167E-10, &     !Band 5
                      -19.432674,     0.37744942,     -0.22166529E-02,  0.11663914E-05,  0.22128830E-07, -0.28943829E-10, &     !Band 6 
                      -51.844021,      1.2280373,     -0.10600353E-01,  0.38135251E-04, -0.45111018E-07,  0.16679671E-10, &     !Band 7
                      -31.210771,     0.85737498,     -0.87947387E-02,  0.39416747E-04, -0.67469797E-07,  0.43711306E-10, &     !Band 8
                      -5.4417604,     0.28970317,     -0.44571665E-02,  0.26395273E-04, -0.52111967E-07,  0.37627129E-10, &     !Band 9
                       14.646543,    -0.25202253,      0.67234738E-03,  0.67552180E-05, -0.19815201E-07,  0.17221281E-10, &     !Band 10
                       12.218584,    -0.31591213,      0.26032011E-02, -0.58878366E-05,  0.73276694E-08, -0.38798834E-11, &     !Band 11
                       1.0183416,    -0.79710154E-01,  0.13753393E-02, -0.40247214E-05,  0.63186167E-08, -0.41250652E-11 /), &  !Band 12
                       (/ 6, MBIR /))

   contains
      subroutine planck  &
         ( NCOL,         &     ! Input:  Number of columns
           NLM,          &     ! Input:  Number of model layers
           nbir,         &     ! Input:  IR band number
           ts,           &     ! Input:  surface temperature, K
           tt,           &     ! Input:  atmospheric layer temperature, K
           bf   )              ! Output: blackbody emission, W/m^2

         integer (kind=int_kind), intent(in) ::   &
            NCOL,   &
            NLM,    &
            nbir
         real (kind=dbl_kind), intent(in), dimension(:) :: &
            ts
         real (kind=dbl_kind), intent(in), dimension(:,:) :: &
            tt
         real (kind=dbl_kind), intent(out), dimension(:,:) :: &
            bf

         
         !Local variables
         integer (kind=int_kind) :: &
            i_lay     !Layer index
         real(kind=dbl_kind), dimension(NCOL) :: &
             tmp      !Tmp var to hold interface temperature

    
         !Blackbody emission at top-of-model                        
         bf(:,1) = b(1,nbir)+tt(:,1)*(b(2,nbir)+tt(:,1)*(b(3,nbir)+tt(:,1)*(b(4,nbir)+tt(:,1)*(b(5,nbir)+tt(:,1)*b(6,nbir)))))
         !Emission at remaining interfaces
         do i_lay = 2,NLM
              tmp(:) = 0.5*(tt(:,i_lay-1)+tt(:,i_lay))
              bf(:,i_lay)= b(1,nbir)+tmp(:)*(b(2,nbir)+tmp(:)*(b(3,nbir)+tmp(:)*(b(4,nbir)+tmp(:)*(b(5,nbir)+tmp(:)*b(6,nbir)))))
           !enddo
         enddo
         !Surface emission
         bf(:,NLM+1) = b(1,nbir)+ts(:)*(b(2,nbir)+ts(:)*(b(3,nbir)+ts(:)*(b(4,nbir)+ts(:)*(b(5,nbir)+ts(:)*b(6,nbir)))))
         return
      end subroutine planck

end module bugsrad_planck
