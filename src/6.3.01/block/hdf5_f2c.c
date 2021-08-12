#include "utils_sub_names.h"

#include "hdf5.h"

hid_t fileid, dsetid, dspcid, mspcid, propid;

/*************************** File routines ***********************************/

/*****************************************************************************/
 void fh5f_open (char*locfn,int*iaccess,int*hdferr)
{

unsigned flags;
hid_t access_id;
extern hid_t fileid;

access_id = H5P_DEFAULT;
if(*iaccess == 1) flags = H5F_ACC_RDONLY;
if(*iaccess == 2) flags = H5F_ACC_RDWR;

fileid = H5Fopen (locfn,flags,access_id);

/*printf("fh5f_open_ - fileid: %d\n",fileid);*/

*hdferr = fileid;

return;
}

/*****************************************************************************/
 void fh5f_create (char*locfn,int*iaccess,int*hdferr)
{

unsigned flags;
hid_t access_id,create_id;
extern hid_t fileid;

access_id = H5P_DEFAULT;
create_id = H5P_DEFAULT;
if(*iaccess == 1) flags = H5F_ACC_TRUNC;
if(*iaccess == 2) flags = H5F_ACC_EXCL ;

fileid = H5Fcreate (locfn,flags,create_id,access_id);

/*printf("fh5f_open_ - fileid: %d\n",fileid);*/

*hdferr = fileid;

return;
}

/*****************************************************************************/
 void fh5f_close (int*hdferr)
{

extern hid_t fileid;
herr_t herr;

herr = H5Fclose (fileid);
/*printf("fh5f_close: %d\n",herr);*/

/*
The following causes problems with REVU since this command forces everything 
to close and flushes all data to disk. This is bad if you have multiple 
processes with hdf5 open. I do not think we need to use this H5close() cmd.
*/
/*herr = H5close ();*/
/*printf("H5_close: %d\n",herr);*/

*hdferr = herr;

return;
}

/**************************** Writing routines *******************************/

/*****************************************************************************/
 void fh5_prepare_write (int*ndims,int*dims,int*hdferr)
{

extern hid_t fileid, dsetid, dspcid, mspcid, propid;
int i;
herr_t herr;

hsize_t dimsc[7]  = {1,1,1,1,1,1,1};
hsize_t maxdims[7]  = {1,1,1,1,1,1,1};
hsize_t chunk_size[7]  = {0,0,0,0,0,0,0};

for (i = 0; i < *ndims; i++) {  dimsc[i] = dims[i]; chunk_size[i] = dims[i];  
                                maxdims[i] = dims[i];}

/*  Create the data space for the dataset. */
mspcid = H5Screate_simple (*ndims,dimsc,maxdims);
/*printf("fh5_prepw - create 1: %d\n",mspcid);*/

/* Create properties for gzip compression.*/
propid = H5Pcreate (H5P_DATASET_CREATE);
/*printf("fh5_prepw - propid: %d\n",propid);*/

herr = H5Pset_chunk (propid,*ndims,chunk_size);
herr = H5Pset_shuffle (propid);
herr = H5Pset_deflate (propid,6);

*hdferr = herr;
/*printf("fh5_prepw - compress: %d\n",mspcid);*/

return;
}

/*****************************************************************************/
 void fh5_write (int*h5type,void*buf,char*dname,int*hdferr)
{

extern hid_t fileid, dsetid, dspcid, mspcid, propid;
herr_t herr;
hid_t memtype;

if (*h5type == 1) memtype=H5T_NATIVE_INT32;
if (*h5type == 2) memtype=HDF5_FLOAT_TYPE;
if (*h5type == 3) memtype=H5T_STRING;
if (*h5type == 4) memtype=H5T_NATIVE_DOUBLE;
if (*h5type == 5) memtype=H5T_NATIVE_HBOOL;
/*printf("fh5_write - start: %d \n",memtype);*/

/* HDF5 1.8 API */
dsetid = H5Dcreate (fileid,dname,memtype,mspcid,H5P_DEFAULT,propid,H5P_DEFAULT);
/*printf("fh5_write - create 1: %d\n",dsetid);*/

herr = H5Dwrite (dsetid,memtype,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf);
/*printf("fh5_write - write 1: %d\n",herr);*/

*hdferr = herr;

return;
}

/*****************************************************************************/
 void fh5_close_write (int*hdferr)
{

extern hid_t dsetid, dspcid, mspcid;
herr_t herr;

herr = H5Sclose (mspcid);
herr = H5Pclose (propid);
herr = H5Dclose (dsetid);
/*printf("fh5_close_write: %d\n",herr);*/

*hdferr = herr;

return;
}
