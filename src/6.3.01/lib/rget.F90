!##############################################################################
!-------------------------------------------------------------------
!                 MACHINE/OS DEPENDENT ROUTINES
!  This contains short utility routines that are not
!  of the FORTRAN 77 standard and may differ from machine to machine
!  or OS to OS. These include bit manipulation, I/O, JCL calls, and 
!  vector funcs.
!-------------------------------------------------------------------

!##############################################################################
Subroutine ugetarg (i,arg)

!MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE

implicit none

      integer :: i
      character(len=*) :: arg
      
      !Routine to get command line argument

#if defined (PC_LINUX1)
      !Use lowercase "call" since this is a system call
      call getarg (i,arg)
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

return
END SUBROUTINE ugetarg

!##############################################################################
Subroutine usystem (arg)

!MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE

implicit none

      character(len=*) :: arg

      !Routine to get command line argument

#if defined (PC_LINUX1)
      !Use lowercase "call" since this is a system call
      call system (trim(arg))
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop
#endif

return
END SUBROUTINE usystem

