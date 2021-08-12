!##############################################################################
Module an_header

implicit none

type head_table
   character(len=32) :: string
   integer :: npointer,idim_type,ngrid,nvalues
end type

type (head_table), allocatable :: anal_table(:)
integer :: nvbtab

END MODULE an_header
