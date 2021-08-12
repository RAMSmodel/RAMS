!##############################################################################
Module var_tables

use grid_dims

implicit none

    !    Define data type for main variable table

type var_tables_r
   
   real, pointer :: var_p,var_m
   integer :: npts, idim_type
   integer :: ianal,imean,ilite,impti,impt1,irecycle_sfc
   character(len=32) :: name
   
end type

    !    Main variable table
type(var_tables_r) :: vtab_r(maxvars,maxgrds)

    !    "nvgrids" is "ngrids", for convenience
integer :: nvgrids

    !    number of variables for each grid
integer :: num_var(maxgrds)



    !    Define data type for scalar variable table

type scalar_table
   
   real, pointer :: var_p,var_t
   character(len=32) :: name
   
end type

    !    Scalar variable table
type(scalar_table) :: scalar_tab(maxsclr,maxgrds)


    !    number of scalars for each grid
integer :: num_scalar(maxgrds)


END MODULE var_tables
