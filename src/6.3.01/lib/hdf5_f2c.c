#include <string.h>
#include <stdint.h>
#include <stdio.h>

#include "utils_sub_names.h"

#include "hdf5.h"

// Size of temporary hsize_t arrays
#define HDF5_MAX_DIMS 10

// Prototypes for internal routines
hid_t  fh5_set_hdf5_dtype (int type_code);
 void fh5_convert_array (int n,int *i_array,hsize_t *h_array,char *dir);

//************************** File routines ***********************************

//****************************************************************************
 void fh5f_open (char *locfn,int *iaccess,int *iphdf5,int64_t *fileid,int *hdferr)
  {
  unsigned flags;
  hid_t access_id;
  hid_t fileid_h5 = (hid_t)(*fileid);
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    // open a parallel file collectively
    access_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (access_id,MPI_COMM_WORLD,MPI_INFO_NULL);
#else
    access_id = H5P_DEFAULT;
#endif
    }
  else
    {
    access_id = H5P_DEFAULT;
    }
  
  if(*iaccess == 1) flags = H5F_ACC_RDONLY;
  if(*iaccess == 2) flags = H5F_ACC_RDWR;
  
  fileid_h5 = H5Fopen (locfn,flags,access_id);
  *fileid = (int64_t)fileid_h5;
  //printf("fh5f_open_ - fileid: %d\n",*fileid);
  if (*fileid < 0) { *hdferr = *fileid; return;}
  
  // Parallel mode created a property which needs to be closed
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    *hdferr = H5Pclose (access_id);
#else
    *hdferr = 0;
#endif
    }
  else
    {
    // no errors if we made it here
    *hdferr = 0;
    }
  
  return;
  }

//*****************************************************************************
 void fh5f_create (char *locfn,int *iaccess,int *iphdf5,int64_t *fileid,int *hdferr)
  {
  unsigned flags;
  hid_t access_id,create_id;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    // open a parallel file collectively
    access_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (access_id,MPI_COMM_WORLD,MPI_INFO_NULL);
#else
    access_id = H5P_DEFAULT;
#endif
    }
  else
    {
    access_id = H5P_DEFAULT;
    }
  create_id = H5P_DEFAULT;
  
  if(*iaccess == 1) flags = H5F_ACC_TRUNC;
  if(*iaccess == 2) flags = H5F_ACC_EXCL ;
  
  fileid_h5 = H5Fcreate (locfn,flags,create_id,access_id);
  *fileid = (int64_t)fileid_h5;
  //printf("fh5f_create: fileid: %d\n",*fileid);
  if (*fileid < 0) { *hdferr = *fileid; return;}
  
  // Parallel mode created a property which needs to be closed
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    *hdferr = H5Pclose (access_id);
#else
    *hdferr = 0;
#endif
    }
  else
    {
    // no errors if we made it here
    *hdferr = 0;
    }
  
  return;
  }

//****************************************************************************
 void fh5f_close (int64_t *fileid,int *hdferr)
  {
  herr_t herr;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);
  
  // Don't use H5close() to close everything since there might be more than one
  // file open (due to REVU for example).
  
  herr = H5Fclose (fileid_h5);
  
  *hdferr = herr;
  
  return;
  }

//******************************* Dataset routines ****************************

//*****************************************************************************
 void fh5d_info (int64_t *fileid,char *dname,int *ndims,int *dims,int *hdferr)
  {
  hid_t dsetid, dspcid;
  herr_t herr;
  hsize_t dimsc[HDF5_MAX_DIMS], maxdimsc[HDF5_MAX_DIMS];
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);

  int i;
  
  // Note that the code is the same whether or not accessing a
  // parallel HDF5 file.
  
  // Open the dataset, HDF5 1.8 API
  dsetid = H5Dopen (fileid_h5,dname,H5P_DEFAULT);
  //printf("fh5d_info: %d\n",dsetid);
  if (dsetid < 0) { *hdferr = dsetid; return;}
  
  // Get the description of the dataset
  dspcid = H5Dget_space (dsetid);
  //printf("fh5d_info: %d\n",dspcid);
  if (dspcid < 0) { *hdferr = dspcid; return;}
  
  // Get the number of dimensions and the sizes of each dimension
  // Place array values in temp space and then copy to output arguments
  *ndims = H5Sget_simple_extent_dims (dspcid, dimsc, maxdimsc);
  //printf("fh5d_info: %d %d %d %d\n",*ndims,dimsc[0],dimsc[1],dimsc[2]);
  if (*ndims < 0) { *hdferr = *ndims; return;}
  
  // Copy hsize_t arrays into FORTRAN arrays (passed in as type int)
  // since HDF5 routines expect hsize_t arrays and the lengths of
  // int and hsize_t may differ.
  fh5_convert_array (*ndims, dims, dimsc, "h2i");
  //printf("fh5d_info: %d %d %d %d\n",*ndims,dimsc[0],dimsc[1],dimsc[2]);
  
  // Close the dataset
  herr = H5Sclose (dspcid);
  herr = H5Dclose (dsetid);
  //printf("fh5d_info: %d\n",herr);
  
  *hdferr = herr;
  
  return;
  }

//****************************** Reading routines *****************************

//*****************************************************************************
// This routine will read a dataset from an HDF5 file. It accommodates hyperslab
// selection on both the input file space and output memory space. The selection
// of memory is contained in the input arguments m_ndims and m_select; ditto for
// the file selection and f_ndims, f_select.
//
// Since passing data in from FORTRAN, the hyperslab selection data has been
// joined together in the linear arrays m_select and f_select. This prevents the
// issue with column-major versus row-major. The data are stored in these
// arrays in groups that are the same size and the corresponding ndims value
// (m_ndims and f_ndims), and go in the order: dims, block, count, offset, stride.
//
// Another complication of the caller being FORTRAN is that the HDF5 routines
// want hsize_t arrays (instead of int arrays). FORTRAN has no idea what hsize_t
// is, so integer arrays are used instead. The lengths of int and hsize_t are
// likely to be different so it is necessary to convert between the argument int
// arrays and internal hsize_t arrays.
//
 void fh5d_read (int64_t *fileid,char *dname,int *h5type,int *m_ndims,int *m_select
               ,int *f_ndims,int *f_select,void *buf,int *hdferr)
  {
  hid_t dsetid, mspcid, fspcid;
  herr_t herr;
  hid_t memtype;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);

  int i;
  
  hsize_t m_dims[HDF5_MAX_DIMS];   // memory selection description
  hsize_t m_block[HDF5_MAX_DIMS];
  hsize_t m_count[HDF5_MAX_DIMS];
  hsize_t m_offset[HDF5_MAX_DIMS];
  hsize_t m_stride[HDF5_MAX_DIMS];

  hsize_t f_dims[HDF5_MAX_DIMS];   // file selection description
  hsize_t f_block[HDF5_MAX_DIMS];
  hsize_t f_count[HDF5_MAX_DIMS];
  hsize_t f_offset[HDF5_MAX_DIMS];
  hsize_t f_stride[HDF5_MAX_DIMS];
  
  // Note that the code is the same whether or not accessing a
  // parallel HDF5 file.
  
  //******************* CONVERT INPUT ARGUMENTS *************************

  // Extract and convert the input int arrays into hsize_t arrays
  // Note: &(array[i]) gives the address of the ith element of array[] making
  // it appear to fh5_convert_array() that this is the zeroth element of
  // its argument array.
  fh5_convert_array (*m_ndims, &(m_select[0]),            m_dims,   "i2h");
  fh5_convert_array (*m_ndims, &(m_select[*m_ndims]),     m_block,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*2]), m_count,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*3]), m_offset, "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*4]), m_stride, "i2h");

  fh5_convert_array (*f_ndims, &(f_select[0]),            f_dims,   "i2h");
  fh5_convert_array (*f_ndims, &(f_select[*m_ndims]),     f_block,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*2]), f_count,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*3]), f_offset, "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*m_ndims)*4]), f_stride, "i2h");

  // Set HDF5 data type
  memtype = fh5_set_hdf5_dtype(*h5type);

  //******************* SET UP SELECTION ******************************
  
  // Create descriptions of the file and memory data spaces

  dsetid = H5Dopen (fileid_h5,dname,H5P_DEFAULT); // file data space
  if (dsetid < 0) { *hdferr = dsetid; return;}
  fspcid = H5Dget_space (dsetid);
  if (fspcid < 0) { *hdferr = fspcid; return;}
  
  mspcid = H5Screate_simple (*m_ndims,m_dims,NULL); // memory data space
  if (mspcid < 0) { *hdferr = mspcid; return;}
  
  // Make the hyperslab selections

  // file
  herr = H5Sselect_hyperslab (fspcid,H5S_SELECT_SET,f_offset,f_stride
                             ,f_count,f_block);
  if (herr < 0) { *hdferr = herr; return;}

  // memory
  herr = H5Sselect_hyperslab (mspcid,H5S_SELECT_SET,m_offset,m_stride
                             ,m_count,m_block);
  if (herr < 0) { *hdferr = herr; return;}
  
  //******************* DO THE READ ******************************

  herr = H5Dread (dsetid,memtype,mspcid,fspcid,H5P_DEFAULT,buf);
  if (herr < 0) { *hdferr = herr; return;}
  
  // Close

  herr = H5Sclose (mspcid);
  herr = H5Sclose (fspcid);
  herr = H5Dclose (dsetid);
  
  *hdferr = herr;
  
  return;
  }

//**************************** Writing routines *******************************

//*****************************************************************************
// This routine will write a dataset into an HDF5 file. It accommodates hyperslab
// selection on both the input file space and output memory space. The selection
// of memory is contained in the input arguments m_ndims and m_select; ditto for
// the file selection and f_ndims, f_select.
//
// See the comments for fh5d_read() for info about the purpose of the hsize_t arrays.
//
 void fh5d_write (int64_t *fileid,char *dname,int *h5type,int *iphdf5,int *m_ndims
                ,int *m_select,int *f_ndims,int *f_select,int *f_csize,void *buf
                ,int *hdferr)
  {
  hid_t dsetid, fspcid, mspcid, propid;
  int i;
  herr_t herr;
  hid_t memtype;
  hid_t fileid_h5;
  fileid_h5 = (hid_t)(*fileid);


  hsize_t m_dims[HDF5_MAX_DIMS];   // memory selection description
  hsize_t m_block[HDF5_MAX_DIMS];
  hsize_t m_count[HDF5_MAX_DIMS];
  hsize_t m_offset[HDF5_MAX_DIMS];
  hsize_t m_stride[HDF5_MAX_DIMS];

  hsize_t f_dims[HDF5_MAX_DIMS];   // file selection description
  hsize_t f_block[HDF5_MAX_DIMS];
  hsize_t f_count[HDF5_MAX_DIMS];
  hsize_t f_offset[HDF5_MAX_DIMS];
  hsize_t f_stride[HDF5_MAX_DIMS];

  hsize_t chunk_size[HDF5_MAX_DIMS]; // for the file data space
  
  //******************* CONVERT INPUT ARGUMENTS *************************

  // Extract and convert the input int arrays into hsize_t arrays
  // Note: &(array[i]) gives the address of the ith element of array[] making
  // it appear to fh5_convert_array() that this is the zeroth element of
  // its argument array.
  fh5_convert_array (*m_ndims, &(m_select[0]),            m_dims,   "i2h");
  fh5_convert_array (*m_ndims, &(m_select[*m_ndims]),     m_block,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*2]), m_count,  "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*3]), m_offset, "i2h");
  fh5_convert_array (*m_ndims, &(m_select[(*m_ndims)*4]), m_stride, "i2h");

  fh5_convert_array (*f_ndims, &(f_select[0]),            f_dims,   "i2h");
  fh5_convert_array (*f_ndims, &(f_select[*f_ndims]),     f_block,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*2]), f_count,  "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*3]), f_offset, "i2h");
  fh5_convert_array (*f_ndims, &(f_select[(*f_ndims)*4]), f_stride, "i2h");

  fh5_convert_array (*f_ndims,f_csize,chunk_size,"i2h");

  // Set HDF5 data type
  memtype = fh5_set_hdf5_dtype (*h5type);

  //******************* CREATE THE FILE DATASET ***********************

  // This must be done before the selection section

  // Create description of data space for the file dataset
  fspcid = H5Screate_simple (*f_ndims,f_dims,f_dims);
  if (fspcid < 0) { 
    *hdferr = fspcid; 
    printf("Error in H5Screate_simple, fspcid: %d", fspcid);
    return;}
  
  // Create properties for chunking, gzip compression. Can't do
  // compression for PHDF5 file.
  propid = H5Pcreate (H5P_DATASET_CREATE);
  //printf("fh5d_write - propid: %d\n",propid);
  if (propid < 0) { 
    *hdferr = propid; 
    printf("Error in H5Pcreate, propid: %d", propid);
    return;}
  
  herr = H5Pset_chunk (propid,*f_ndims,chunk_size);
  #if H5_VERSION_GE(1,10,3)
  #ifdef ENABLE_PARALLEL_COMPRESSION 
  
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
  
  #else
  if (*iphdf5 != 1)
    {
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
    }
  #endif//enable parallel compression ifdef
  #else//H5_VERSION_GE
  if (*iphdf5 != 1)
    {
    herr = H5Pset_shuffle (propid);
    herr = H5Pset_deflate (propid,6);
    }
    #endif //H5_VERSION_GE if
  if (herr < 0) { 
    *hdferr = herr; 
    printf("Error in chunking, herr: %d", herr);
    return;}
  
  // Create the dataset
  dsetid = H5Dcreate (fileid_h5,dname,memtype,fspcid,H5P_DEFAULT,propid,H5P_DEFAULT);
  if (dsetid < 0) { 
    *hdferr = dsetid; 
    printf("Error in H5Dcreate, dsetid: %lld \n",dsetid);
    printf("Fileid_h5: %lld, dname: %s \n",fileid_h5,dname);
    printf("memtype: %d, fspcid: %d, propid: %d \n",memtype,fspcid, propid);
    for(int i=0; i<*f_ndims; i++){
      printf("Chunk size [%d], %lld \n",i,chunk_size[i]);
    }
    return;}

    

  // Close the items used for this section, but keep dsetid open for
  // subsequent use
  H5Pclose (propid);
  H5Sclose (fspcid);
  
  //******************* SET UP SELECTION ******************************
  
  // Create descriptions of the file and memory data spaces

  fspcid = H5Dget_space (dsetid);
  if (fspcid < 0) { *hdferr = fspcid; return;}
  
  mspcid = H5Screate_simple (*m_ndims,m_dims,NULL); // memory data space
  if (mspcid < 0) { *hdferr = mspcid; return;}
  
  // Make the hyperslab selections

  // file
  herr = H5Sselect_hyperslab (fspcid,H5S_SELECT_SET,f_offset,f_stride
                             ,f_count,f_block);
  if (herr < 0) { *hdferr = herr; return;}

  // memory
  herr = H5Sselect_hyperslab (mspcid,H5S_SELECT_SET,m_offset,m_stride
                             ,m_count,m_block); 
  if (herr < 0) { *hdferr = herr; return;}
  
  //******************* DO THE WRITE ******************************

  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    propid = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (propid,H5FD_MPIO_COLLECTIVE);
#else
    propid = H5P_DEFAULT;
#endif
    }
  else
    {
    propid = H5P_DEFAULT;
    }
  
  herr = H5Dwrite (dsetid,memtype,mspcid,fspcid,propid,buf);
  if (herr < 0) { *hdferr = herr; return;}
  
  // Close
  if (*iphdf5 == 1)
    {
#if defined (RAMS_MPI)
    H5Pclose (propid);
#endif
    }
  herr = H5Sclose (mspcid);
  herr = H5Sclose (fspcid);
  herr = H5Dclose (dsetid);
  
  *hdferr = herr;
  
  return;
  }

//*************************** Utilities **************************************

//****************************************************************************
// This routine will translate the type codes from FORTAN to the built in
// HDF5 type codes.

hid_t fh5_set_hdf5_dtype (int type_code)
  {
  hid_t dtype;

  if (type_code == 1) dtype = H5T_NATIVE_INT32;
  if (type_code == 2) dtype = HDF5_FLOAT_TYPE;
  if (type_code == 3) dtype = H5T_STRING;
  if (type_code == 4) dtype = H5T_NATIVE_DOUBLE;
  if (type_code == 5) dtype = H5T_NATIVE_HBOOL;

  return dtype;
  }

//****************************************************************************
// This routine will copy the contents between an integer array and an hsize_t
// array. It is likely to encounter a system where hsize_t and int are not the same
// length. As an example, this would necessitate the need to convert back and
// forth when the integer array is used as an argument from FORTRAN for example.
//
// dir == 'i2h' --> convert int to hsize_t
// dir == 'h2i' --> convert hsize_t to int

 void fh5_convert_array (int n,int *i_array,hsize_t *h_array,char *dir)
  {
  int i;

  if (strcmp(dir, "i2h") == 0)
    {
    // copy int array into the hsize_t array
    for (i = 0; i < n; i++)
      {
      h_array[i] = i_array[i];
      }
    }
  else if (strcmp(dir, "h2i") == 0)
    {
    // copy hsize_t array into the int array
    for (i = 0; i < n; i++)
      {
      i_array[i] = h_array[i];
      }
    }
  else
    {
    printf ("ERROR: fh5_convert_array:");
    printf ("Do not recognize dir argument value: %d\n", dir);
    }

  return;
  }
