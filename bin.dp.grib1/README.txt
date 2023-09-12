1. Type "make.sc" to compile after you have setup your ../include.mk file.

2. Run sample executable like "dgrib-6.2.11" to see instructions on how to
   pre-process a single Grib-1 file.

3. Not all Grib reanalysis and forecast gridded datasets are accounted for. 

   You can view the degribbing code in the files such as:
     src/dprep/dgrib1_main.f90   and
     src/lib/griber_grb1.c

   In these routines are the specfications for adding a new dataset. Grib
   names and numbers and labels can vary among datasets and each variation
   has to be customized as a new datatype.

4. Note that rams_grib.c is just wgrib.c with some modifications for output.
   Currently using wgrib version 1.8.2 downloaded Sep 18, 2019.
