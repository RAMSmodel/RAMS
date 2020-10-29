#if defined(PC_LINUX1)

/* LINUX APPENDS AN UNDERSCORE TO C ROUTINES */ 
/* Other machines/OS might require no-underscore or all capitals */

#define grib_query grib_query_
#define grib_queryc grib_queryc_
#define grib_queryf grib_queryf_
#define grib_queryi grib_queryi_
#define grib_get grib_get_

#else

   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify filelist.F90 to add your machine/OS"
   stop

#endif
