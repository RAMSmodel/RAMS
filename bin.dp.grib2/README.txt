1. Type "make" to compile after you have setup your ../include.mk file.

2. Run sample executable like "dgrib-6.2.11" to see instructions on how to
   pre-process a single Grib-2 file.

3. Not all Grib reanalysis and forecast gridded datasets are accounted for.

   You can view the degribbing code in the files such as:
     src/dprep/dgrib2_main.f90   and
     src/lib/griber_grb2.c

   In these routines are the specfications for adding a new dataset. Grib
   names and numbers and labels can vary among datasets and each variation
   has to be customized as a new datatype.
