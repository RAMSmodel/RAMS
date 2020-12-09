#if defined (PC_LINUX1)

/* LINUX APPENDS AN UNDERSCORE TO C ROUTINES */
/* Other machines/OS might require no-underscore or all capitals */

#define c_listfile          c_listfile_
#define irsleep             irsleep_
#define par_init_fortran    par_init_fortran_
#define par_init_put        par_init_put_
#define par_init_recv_bcast par_init_recv_bcast_
#define par_broadcast       par_broadcast_
#define par_allreduce_sum_doubles par_allreduce_sum_doubles_
#define par_gather_doubles  par_gather_doubles_
#define par_gather_floats   par_gather_floats_
#define par_gather_ints     par_gather_ints_
#define par_send            par_send_
#define par_put_int         par_put_int_
#define par_put_double      par_put_double_
#define par_put_float       par_put_float_
#define par_put_char        par_put_char_
#define par_send_noblock    par_send_noblock_
#define par_get_noblock     par_get_noblock_
#define par_assoc_buff      par_assoc_buff_
#define par_wait            par_wait_
#define par_get_new         par_get_new_
#define par_get_int         par_get_int_
#define par_get_float       par_get_float_
#define par_get_double      par_get_double_
#define par_get_char        par_get_char_
#define par_exit            par_exit_
#define par_pause           par_pause_
#define par_ready           par_ready_
#define fh5f_open           fh5f_open_
#define fh5f_create         fh5f_create_
#define fh5f_close          fh5f_close_
#define fh5d_info           fh5d_info_
#define fh5d_read           fh5d_read_
#define fh5d_write          fh5d_write_
#define fh5_prepare_write   fh5_prepare_write_
#define fh5_write           fh5_write_
#define fh5_close_write     fh5_close_write_

// Single versus Double Precision variables.
// Specify ifdef statement in include.mk file.
#ifdef RAMS_DOUBLE_PREC
#define HDF5_FLOAT_TYPE H5T_NATIVE_DOUBLE
#define MPI_FLOAT_TYPE MPI_DOUBLE
#define C_FLOAT_TYPE double
#else
#define HDF5_FLOAT_TYPE H5T_NATIVE_FLOAT
#define MPI_FLOAT_TYPE MPI_FLOAT
#define C_FLOAT_TYPE float
#endif

#else

   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop

#endif
