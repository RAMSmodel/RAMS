/*
 * Using these C wrapper routines to be able to utilize the C API for HDF5.
 * Unfortunately, the fortran .mod files that came with our installation of 
 * HDF5 do not work with the pgf90 compiler.
 */

#include "revu_sub_names.h"
#include <string.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "utils_sub_names.h"
#include <stdio.h>

/*
 * Limits for arrays, strings, etc. Need to keep these in sync with like
 * named parameters in rhdf5_utils.f90
 */
#define RHDF5_MAX_DIMS    10

/*
 * Integer coding for HDF5 types. Need to keep these in sync with like named
 * parameters in rhdf5_utils.f90
 * These codes need to start with zero and be contiguous (0..n with no gaps)
 * for the rhdf5_type_names array. Make sure the entries in rhdf5_type_names
 * are consistent with the type numbers.
 */
#define RHDF5_NUM_TYPES     4

#define RHDF5_TYPE_STRING   0
#define RHDF5_TYPE_INTEGER  1
#define RHDF5_TYPE_FLOAT    2
#define RHDF5_TYPE_CHAR     3

char *rhdf5_type_names[RHDF5_NUM_TYPES] = {"STRING","INTEGER","FLOAT","CHAR"};

/* Prototypes for internal routines */

 void rf2c_get_mem_type (hid_t typid,int *type,int *size);
hid_t rf2c_set_mem_type (int dtype,int ssize);

/* Routines for HDF5 REVU IO. */

//############################################################################
//************************* DATASET ROUTINES *********************************
//############################################################################

//****************************************************************************
 void rh5d_open (int64_t *id,char *name,int64_t *dsetid,int *hdferr)
{
  hid_t fileid_h5 = (hid_t)(*id);
  hid_t dsetid_h5 = (hid_t)(*dsetid);
  dsetid_h5 = H5Dopen (fileid_h5,name,H5P_DEFAULT);
  *dsetid = (int64_t) dsetid_h5;
if (*dsetid < 0)
  {
  *hdferr = *dsetid;
  }
else
  {
  *hdferr = 0;
  }
return;
}

//****************************************************************************
 void rh5d_setup_and_write (int64_t *id,char *name,int *dtype,int *ssize
                         ,int *ndims,int *dims,int *ext_dims,int *chunk_dims
                         ,int *deflvl,int64_t *dsetid,void *data,int *hdferr)
{
/*
 * This routine will either create or extend a dataset and then
 * perform the write into that dataset.
 */
int i;
int cur_ndims;
hid_t fileid_h5 = (hid_t)(*id);
hid_t dsetid_h5 = (hid_t)(*dsetid);
hsize_t cur_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t max_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t new_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t mem_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t chunk_sizes[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
hsize_t hs_count[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t hs_start[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hid_t mtypid;
hid_t mspcid;
hid_t fspcid;
hid_t dcplid;

*hdferr = 0;

/*
 * Set the memory type (translate integer code in dtype to HDF5 type)
 */
mtypid = rf2c_set_mem_type (*dtype, *ssize);
if (mtypid < 0)
  {
  printf("ERROR: rh5d_setup_and_write: unrecognized dtype: %d\n", *dtype);
  *hdferr = -1;
  return;
  }

/*
 * Check to see if the dataset exists aleady. If so extend it and do
 * the write. If not, then create it and do the write.
 */
  
if ( H5Lexists (fileid_h5,name,H5P_DEFAULT) )
  {
  /*
   * dataset exists --> extend and write
   *
   * Open the existing dataset
   */

  dsetid_h5 = H5Dopen (fileid_h5,name,H5P_DEFAULT);
  *dsetid = (int64_t)(dsetid_h5);

  /*
   * Get the currently existing dimension information from the dataset
   * Close the dataspce after grabbing the current dimension information
   * since we need to reopen the dataspace after extending the dataset.
   */
  fspcid = H5Dget_space (dsetid_h5);
  if (fspcid < 0)
    {
    printf("ERROR: rh5d_setup_and_write: Cannot get file space id\
           from dataset\n");
    *hdferr = -1;
    return;
    }
  cur_ndims = H5Sget_simple_extent_dims (fspcid, cur_dims, max_dims);
  H5Sclose (fspcid);
  
  if (cur_ndims != *ndims)
    {
    printf("ERROR: rh5d_setup_and_write: number of dimensions in dataset\
           (%d) does not match specified number of dimensions (ndims): \
           %d\n",cur_ndims, *ndims);
    *hdferr = -1;
    return;
    }
  
  /*
   * Get the new dimension information from the input arguments
   */
  for (i = 0; i<*ndims; i++)
    {
    new_dims[i] = dims[i];
    /* hs_* vars are for the hyperslab selection later on */
    hs_count[i] = cur_dims[i];
    hs_start[i] = 0;
    if (ext_dims[i] == 1)
      {
      /* this is the dimension we are extended and assume for now */
      /* that we are extending by just one element                */
      mem_dims[i] = 1;

      if (new_dims[i] < cur_dims[i])
        {
        printf("ERROR: rh5d_setup_and_write: Shrinking an extendible \
                dimension is not supported\n");
        printf("ERROR:   Current dimension number: %d, Current dimension \
                size: %d\n", i+1, (int) cur_dims[i]);
        printf("ERROR:   New dimension number:     %d, New dimension size: \
                %d\n", i+1, (int) new_dims[i]);
        *hdferr = -1;
        return;
        }
      }
    else
      {
      mem_dims[i] = dims[i];
      }
    }
  
  /*
   * Create a memory dataspace to control the selection from the data 
   * buffer (input arg "data"). We are assuming that the entries in the
   * data buffer are organized so that the dimension that we are extending
   * is a size of one so if we set up the memory space with that in mind
   * the selection will happen automatically without the need to apply any 
   * selection mechanism on the memory space.
   */
  mspcid = H5Screate_simple (cur_ndims,mem_dims,NULL);
  
  /*
   * Extend the dataset and re-open the dataspace
   */
  H5Dset_extent (dsetid_h5,new_dims);
  fspcid = H5Dget_space (dsetid_h5);
  *dsetid = (int64_t)(dsetid_h5);
  
  if ((mspcid < 0) || (fspcid < 0))
    {
    *hdferr = -1;
    return;
    }
  
  /*
   * Create a hyperslab to select just the region where the new data will be
   * stored into the file. That is, we are appending the new data just after
   * the end of the existing data so select out the existing data. To do this
   * first select the entire dataspace (extended) and subtract out the 
   * dataspace prior to the extension. The select operator H5S_SELECT_NOTB
   * is used to do the subtraction.
   */
  H5Sselect_all(fspcid);
  *hdferr = H5Sselect_hyperslab (fspcid,H5S_SELECT_NOTB,hs_start
                                ,NULL,hs_count,NULL);
  
  if (*hdferr != 0)
    {
    return;
    }

  /*
   * write and free up resources
   */
  *hdferr = H5Dwrite (dsetid_h5,mtypid,mspcid,fspcid,H5P_DEFAULT,data);

  if (*dtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (mtypid);
    }
  H5Sclose (mspcid);
  H5Sclose (fspcid);
  }
else
  {
  /*
   * dataset does not exist --> create and write
   *
   * H5Screate_simple wants all of the elements of the dims and maxdim 
   * arrays (2nd, 3rd args) to be set to a non-zero value, even the elements
   * that are not being used (since they are beyond what ndims specifies).
   * Copy the used contents of dims (the first ndims entries) into an internal 
   * array which is initialized to all 1's so that the caller doesn't need
   * to bother with this detail.
   *
   * H5Pset_chunk wants elements not used in the chunk size array (2nd arg) 
   * to be set to zero and the used elements set to non-zero values.
   */
  for (i=0 ; i<*ndims; i++)
    {
    cur_dims[i] = dims[i];
    if (ext_dims[i] == 0)
      {
      max_dims[i] = dims[i];
      }
    else
      {
      max_dims[i] = H5S_UNLIMITED;
      }
    chunk_sizes[i] = chunk_dims[i];
    }
  
  /*
   * create the file dataspace
   */
  mspcid = H5Screate_simple (*ndims,cur_dims,max_dims);
  if (mspcid < 0)
    {
    printf("ERROR: rh5d_create: Cannot create memory space for dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * create the dataset
   * dataset creation properties
   */
  dcplid = H5Pcreate (H5P_DATASET_CREATE);
  if (dcplid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create properties for dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * Set up for data compression
   *   enable chunking of data
   *   enable shuffling of data
   *   enable deflation (compression) of data
   */
  H5Pset_chunk (dcplid,*ndims,chunk_sizes);
  H5Pset_shuffle (dcplid);
  H5Pset_deflate (dcplid,*deflvl);
  
  /*
   * dataset
   */

  /* 3rd arg is type of the data elements - always using float for now
   * 5th arg is the link creation property list (eg, automatically build 
   *   intermediate groups)
   * 7th arg is the data access property list
   */
  dsetid_h5 = H5Dcreate (fileid_h5,name,mtypid,mspcid,H5P_DEFAULT,dcplid,H5P_DEFAULT);
  *dsetid = (int64_t)(dsetid_h5);
  if (*dsetid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * write and free up resources
   */
  *hdferr = H5Dwrite (dsetid_h5,mtypid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

  if (*dtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (mtypid);
    }
  H5Pclose (dcplid);
  H5Sclose (mspcid);
  }

return;
}

//****************************************************************************
 void rh5d_read_get_dims (int64_t *dsetid,int *ndims,int *dims,int *hdferr)
{
/*
 * This routine will open the dataset given by id (file id) and name, 
 * read in and return the dimension information.
 *
 * The dataset is left open so that the caller can read the data and
 * retrieve attributes so it is up to the caller to close the dataset.
 */
  hid_t fspcid;
  hid_t dtypid;
  hid_t dsetid_h5 = (hid_t)(*dsetid);
  hsize_t dset_dims[RHDF5_MAX_DIMS];
  int i;

  *hdferr = 0;
  
  /*
   * Get the information about the dimensions
   * Need to read dimension sizes, from H5Sget_simple_extent_dims, into
   * an (hsize_t *) type since it's length is different than (int *) type.
   */
  fspcid = H5Dget_space (dsetid_h5);
  *ndims = H5Sget_simple_extent_dims (fspcid,dset_dims,NULL);
  for (i=0; i<*ndims; i++)
    {
    dims[i] = dset_dims[i];
    }

  H5Sclose (fspcid);
  if (*ndims < 0)
    {
    *hdferr = -1;
    }

  return;
}

//****************************************************************************
 void rh5d_close (int64_t *dsetid,int *hdferr)
{
*hdferr = H5Dclose ((hid_t)*dsetid);

return;
}

//############################################################################
//******************* ATTRIBUTE ROUTINES *************************************
//############################################################################

//****************************************************************************
 void rh5a_write_anyscalar (int64_t *id,char *name,void *value,int *vtype,int *hdferr)
{
/* This routine will write a scalar attribute into the HDF5 file, given the 
 * object id and attribute name.
 *
 * value is a VOID pointer so the caller can use any type (char, int, float 
 * for now). vtype is used to denote what type value is. The encoding for 
 * vtype is contained in the the defines: RHDF5_TYPE_* (see top of this file).
 */
  hid_t a_id, amem_id, atype;
  hid_t id_h5 = (hid_t)(*id);
  htri_t attr_exists;
  int ssize;
  
  *hdferr = 0;
  
  /*
   * Set up memory type
   */
  if (*vtype == RHDF5_TYPE_STRING)
    {
    ssize = strlen(value) + 1;
    }
  else
    {
    ssize = 0;
    }
  
  atype = rf2c_set_mem_type (*vtype, ssize);
  
  if (atype < 0)
    {
    printf("ERROR: rh5a_write_anyscalar: unrecognized vtype: %d\n", *vtype);
    *hdferr = -1;
    return;
    }
  
  if (*hdferr == 0)
    {
    /*
     * got a recognized type so keep going
     */
  
    attr_exists = H5Aexists (id_h5,name);
  
    if (attr_exists)
      {
      /*
       * attribute exists already --> update mode
       */
  
      a_id = H5Aopen (id_h5,name,H5P_DEFAULT);
      }
    else
      {
      /*
       * attribute does not exist --> create mode
       *
       *
       * scalar dataspace to hold the string for the file write
       */
      amem_id = H5Screate (H5S_SCALAR);
    
      /*
       * create the attribute and save the result for returning the status in hdferr
       * set hdferr here since we will be closing the a_id pointer at the end
       */

      /* HDF5 1.8 API
       * 5th arg is attribute creation property list (must be H5P_DEFAULT as of 4/19/12)
       * 6th arg is attribute access property list (must be H5P_DEFAULT as of 4/19/12)
       */
      a_id = H5Acreate (id_h5,name,atype,amem_id,H5P_DEFAULT,H5P_DEFAULT);
      }
  
    if (a_id < 0)
      {
      *hdferr = a_id;
      }
  
    /*
     * If we have an error don't attempt the write, but do continue on to
     * the close statments.
     */
  
    if (*hdferr == 0)
      {
      /*
       * The third arg to H5Awrite is a pointer to the data. value is that already,
       * but this does put the onus on the caller to make sure that value is really
       * pointing to the type that matches what's in vtype
       */
      H5Awrite (a_id,atype,value);
      }
  
    /*
     * clean up the attribute related pointers
     */
    if (*vtype == RHDF5_TYPE_STRING)
      {
      /*
       * Only call the close when using the string type (which is a complex structure)
       */
      H5Tclose (atype);
      }
    if (attr_exists == 0)
      {
      /*
       * if created the attribute, then close the memory space id
       */
      H5Sclose (amem_id);
      }
    H5Aclose (a_id);
    }
  
  return;
}

//****************************************************************************
 void rh5a_read_anyscalar (int64_t *dsetid,char *name,void *value,int *vtype,int *hdferr)
{
/* This routine will read in the attribute given by the object id and
 * name. The caller requests a certain type (arg: dtype) and that is checked
 * against the type found in the hdf5 file.
 */
  hid_t attrid;
  hid_t atypid;
  hid_t dsetid_h5 = (hid_t) *dsetid;
  int atype;
  int asize;
  
  *hdferr = 0;
  
  /*
   * First open the attribute
   */
  attrid = H5Aopen (dsetid_h5,name,H5P_DEFAULT);
  if (attrid < 0)
    {
    *hdferr = -1;
    return;
    }

  /*
   * Get the information about the data type
   */
  atypid = H5Aget_type (attrid);
  if (atypid < 0)
    {
    *hdferr = -1;
    return;
    }
  rf2c_get_mem_type (atypid,&atype,&asize);
  if (atype < 0)
    {
    *hdferr = -1;
    return;
    }
  if (atype != *vtype)
    {
    printf("ERROR: rh5a_read_anyscalar: Attribute type requested by caller");
    printf("       does not match type in file\n");
    printf("ERROR: Attribute naame: %s\n", name);
    printf("ERROR: Requested attribute type: %s\n",rhdf5_type_names[*vtype]);
    printf("ERROR: Attribute type found in HDF5 file: %s\n",rhdf5_type_names[atype]);
    *hdferr = -1;
    return;
    }

  /*
   * read the dataset, clean up and return
   */
  *hdferr = H5Aread (attrid,atypid,value);

  if (*vtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (atypid);
    }
  H5Aclose (attrid);

  return;
}

//############################################################################
//*************** DIMENSION SCALE ROUTINES ***********************************
//############################################################################

//****************************************************************************
 void rh5ds_set_scale (int64_t *dsetid,char *dimname,int *hdferr)
  {
  *hdferr = H5DSset_scale ((hid_t)*dsetid,dimname);

  return;
  }

//****************************************************************************
 void rh5ds_attach_scale (int64_t *dsetid,int64_t *dsclid,int *index,int *hdferr)
  {
  *hdferr = H5DSattach_scale ((hid_t)*dsetid,(hid_t)*dsclid,*index);

  return;
  }

//############################################################################
//*************** INTERNAL ROUTINES ******************************************
//############################################################################

//****************************************************************************
hid_t rf2c_set_mem_type (int dtype,int ssize)
{
/* This routine translates the integer coded memory type to the
 * corresponding HDF5 type.
 *
 * Note that ssize input argument is only used for string types.
 */
hid_t mtype;

switch (dtype)
  { 
  case RHDF5_TYPE_STRING:
    /*
     * define type for the string
     * C style string, null byte terminated, size taken from input argument ssize
     */
    mtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (mtype,ssize);
    H5Tset_strpad (mtype,H5T_STR_NULLTERM);
    break;

  case RHDF5_TYPE_INTEGER:
    /*
     * integer
     */
    mtype = H5T_NATIVE_INT;
    break;

  case RHDF5_TYPE_FLOAT:
    /*
     * float
     */
    mtype = HDF5_FLOAT_TYPE;
    break;

  case RHDF5_TYPE_CHAR:
    /*
     * char
     */
    mtype = H5T_NATIVE_UCHAR;
    break;

  default:
    /*
     * unrecognized type, send back -1 to inform caller
     */
    mtype = -1;
    break;
  }

return(mtype);
}

//****************************************************************************
 void rf2c_get_mem_type (hid_t typid,int *type,int *size)
{
/* This routine finds and translates the HDF5 type to
 * the integer coded memory type.
 */

  /* First check the class to see if this is a string type
   * since string type is not a "native" type. All the other
   * types (integer, float, char) are native types which can
   * checked using H5Tequal. H5Tequal returns positive value
   * for "true", zero for "false" and negative value for "error".
   */

  if ( H5Tget_class (typid) == H5T_STRING )
    {
    *type = RHDF5_TYPE_STRING;
    }
  else
    {
    /*
     * Not a string type so check the native types
     *
     * Check the character types first since they are a subset of
     * the integer types.
     */

    if (( H5Tequal (typid, H5T_NATIVE_UCHAR) > 0) || 
        ( H5Tequal (typid, H5T_NATIVE_SCHAR) > 0 ))
      {
      *type = RHDF5_TYPE_CHAR;
      }
    else if ( H5Tequal (typid,H5T_NATIVE_INT) > 0 )
      {
      *type = RHDF5_TYPE_INTEGER;
      }
    else if ( H5Tequal (typid,H5T_NATIVE_FLOAT) > 0 )
      {
      *type = RHDF5_TYPE_FLOAT;
      }
    else
      {
      /*
       * unrecognized type
       */
      *type = -1;
      }
    }

  /*
   * grab the size of this type
   */
  *size = H5Tget_size (typid);

  return;
}
