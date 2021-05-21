/**
 * @substitutes for the needed functionality previously obtained from minc
 */
/*
 * Overhaul Author: Bevin Brett
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

// HEAVILY BASED on minc-1.5.1/volume_io/Volumes/multidim_arrays.c

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "minc_multidim_arrays.h"

#include  <limits.h>
#include  <float.h>
#include  <stdio.h>
#include  <stdbool.h>
#include  <string.h>
#include "minc_internals.h"
	// various prototypes
	
#define for_less(VAR,INIT,LIM) for (VAR=INIT; VAR < LIM; VAR++)

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_empty_multidim_array
@INPUT      : n_dimensions
              data_type
@OUTPUT     : array
@RETURNS    : 
@DESCRIPTION: Creates a multidimensional array, without allocating its data.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void   create_empty_multidim_array(
    VIO_multidim_array  *array,
    int             n_dimensions,
    VIO_Data_types  data_type )
{
    if( n_dimensions < 1 || n_dimensions > VIO_MAX_DIMENSIONS )
    {
        fprintf(stderr,
     "create_empty_multidim_array(): n_dimensions (%d) not in range 1 to %d.\n",
               n_dimensions, VIO_MAX_DIMENSIONS );
    }

    array->n_dimensions = n_dimensions;
    array->data_type = data_type;
    array->data = (void *) 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_multidim_data_type
@INPUT      : array
@OUTPUT     : 
@RETURNS    : data type
@DESCRIPTION: Returns the data type of the multidimensional array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIO_Data_types  get_multidim_data_type(
    VIO_multidim_array       *array )
{
    return( array->data_type );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_multidim_data_type
@INPUT      : array
              data_type
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the data type of the array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  set_multidim_data_type(
    VIO_multidim_array       *array,
    VIO_Data_types        data_type )
{
    array->data_type = data_type;
}

/* ----------------------------------------------------------------------------
@NAME       : get_type_size
@INPUT      : type
@OUTPUT     : 
@RETURNS    : size of the type
@DESCRIPTION: Returns the size of the given type.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June, 1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

int  get_type_size(
    VIO_Data_types   type )
{
    int   size;

    switch( type )
    {
    case  UNSIGNED_BYTE:    size = sizeof( unsigned char );   break;
    case  SIGNED_BYTE:      size = sizeof( signed   char );   break;
    case  UNSIGNED_SHORT:   size = sizeof( unsigned short );  break;
    case  SIGNED_SHORT:     size = sizeof( signed   short );  break;
    case  UNSIGNED_INT:     size = sizeof( unsigned int );    break;
    case  SIGNED_INT:       size = sizeof( signed   int );    break;
    case  FLOAT:            size = sizeof( float );           break;
    case  DOUBLE:           size = sizeof( double );          break;
    default: size = 0; fprintf(stderr, "%s:%d no default\n", __FILE__,__LINE__); exit(1);
    }

    return( size );
}

void  get_type_range(
    VIO_Data_types   type,
    double         *min_value,
    double         *max_value )
{
    switch( type )
    {
    case UNSIGNED_BYTE:
        *min_value = 0.0;
        *max_value = (double) UCHAR_MAX;     break;
    case SIGNED_BYTE:
        *min_value = (double) SCHAR_MIN;
        *max_value = (double) SCHAR_MAX;     break;
    case UNSIGNED_SHORT:
        *min_value = 0.0;
        *max_value = (double) USHRT_MAX;     break;
    case SIGNED_SHORT:
        *min_value = (double) SHRT_MIN;
        *max_value = (double) SHRT_MAX;      break;
    case UNSIGNED_INT:
        *min_value = 0.0;
        *max_value = (double) UINT_MAX;     break;
    case SIGNED_INT:
        *min_value = (double) INT_MIN;
        *max_value = (double) INT_MAX;      break;
    case FLOAT:
        *min_value = (double) -FLT_MAX;
        *max_value = (double) FLT_MAX;       break;
    case DOUBLE:
        *min_value = (double) -DBL_MAX;
        *max_value = (double) DBL_MAX;       break;
    default: fprintf(stderr, "%s:%d no default\n", __FILE__,__LINE__); exit(1);
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_multidim_sizes
@INPUT      : array
              sizes
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the sizes of the array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  set_multidim_sizes(
    VIO_multidim_array   *array,
    int              sizes[] )
{
    int    dim;

    for_less( dim, 0, array->n_dimensions )
        array->sizes[dim] = sizes[dim];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_multidim_sizes
@INPUT      : array
@OUTPUT     : sizes
@RETURNS    : 
@DESCRIPTION: Passes back the sizes of the multidimensional array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  get_multidim_sizes(
    VIO_multidim_array   *array,
    int              sizes[] )
{
    int   i;

    for_less( i, 0, array->n_dimensions )
        sizes[i] = array->sizes[i];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : multidim_array_is_alloced
@INPUT      : array
@OUTPUT     : 
@RETURNS    : TRUE if array is allocated
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

bool  multidim_array_is_alloced(
    VIO_multidim_array   *array )
{
    return( array->data != NULL );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : alloc_multidim_array
@INPUT      : array
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Allocates the data for the multidimensional array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static void alloc_multidim_array_debug() {
}

void  alloc_multidim_array(
    VIO_multidim_array   *array )
{
    alloc_multidim_array_debug();

    if( multidim_array_is_alloced( array ) )
        delete_multidim_array( array );

    if( array->data_type == NO_DATA_TYPE )
    {
       fprintf(stderr,
           "Error: cannot allocate array data until size specified at %s:%d.\n",__FILE__,__LINE__ );
	exit(1);
        return;
    }

    // One malloc call allocates both the pointers and the data
    // The pointers are a tree, the leaves pointing to vectors of n-1 dimension number of elements
    //
    // For instance, a 4D array size[0] = 2, size[1] = 3, size[2] = 5, size[3] = 7
    // is set up like this...
    //
    // array->data  ->  0:  -> 2:  ->  8:      -> elements  0.. 6
    //				       9:      ->           7..
    //				    ..12:      ->
    //                         3:  -> 13:
    //                              ..17:
    //                         4:  ->
    //
    //                  1:  -> 5:  ->
    //			       6:  ->
    //			       7:  ->
    //				    ..37:      -> elements  ..209
    //
    int dim;

    // Calculate the number of elements needed
    //
    size_t numberOfElements = 1;
    for_less( dim, 0, array->n_dimensions )
        numberOfElements *= (size_t) array->sizes[dim];			// e.g.: 2*3*5*7 = 210

    // Calculate the number of pointers needed, 
    // not including the root pointer  
    //
    size_t numberOfPointers       = 0;
    size_t numberOfPointersPerDim = 1;
    for_less( dim, 0, array->n_dimensions - 1 )
        numberOfPointers += (numberOfPointersPerDim *= (size_t) array->sizes[dim]);
    									// e.g.: 2 + 6 + 30 = 38
    
    // Calculate the size of the pointers, rounding up to a multiple of
    // the cache line size to keep the elements cache line aligned
    //
    size_t dataTypeSize    = (size_t) get_type_size( array->data_type );
    size_t bytesOfPointers = (numberOfPointers*sizeof(void*) + 64) & ~63;
    size_t bytesOfData     = (numberOfElements*dataTypeSize  + 64) & ~63;
    
    // Malloc the space, and separate it into pointers and elements
    //
    void* data = malloc(bytesOfPointers + bytesOfData);
    char* elts = (char*)data + bytesOfPointers;
    
    // Fill in the dim 0, and get ready to fill in dim 1
    //
    array->data = data;    
    size_t prevNumberOfPtrsToFillIn = 1;
    void** ptrsToFillIn = (void**)data;

    // Fill in dim 1 .. n-2
    //
    for_less( dim, 1, array->n_dimensions - 1 ) {
    	// e.g: dim is 1 and 2
	
        size_t numberOfPtrsToFillIn = prevNumberOfPtrsToFillIn*(size_t) array->sizes[dim - 1];
		// e.g.: when dim = 1, this is 1 * 2 = 2
		//	      dim = 2, this is 2 * 3 = 6
		
	void** nextLayerOfPointers  = ptrsToFillIn + numberOfPtrsToFillIn;
		// e.g.: when dim = 1, this is data + 2
		//	      dim = 2, this is data + 2 + 6
	
	size_t i;
	for (i = 0; i < numberOfPtrsToFillIn; i++) {
	    ptrsToFillIn[i] = nextLayerOfPointers + i*(size_t) array->sizes[dim];
		// e.g.: when dim = 1, this is data + 2     + {0,3}
		//	      dim = 2, this is data + 2 + 6 + {0,5,...}
	}
	
	prevNumberOfPtrsToFillIn = numberOfPtrsToFillIn;
	ptrsToFillIn             = nextLayerOfPointers;
		// e.g.: when dim = 1, this is data + 2    
		//	      dim = 2, this is data + 2 + 6
    } 
    
    // Fill in the pointers to the elements
    //
    dim = array->n_dimensions - 1;
    {
    	size_t numberOfPtrsToFillIn = prevNumberOfPtrsToFillIn*(size_t) array->sizes[dim - 1];
	size_t i;
	for (i = 0; i < numberOfPtrsToFillIn; i++) {
	    ptrsToFillIn[i] = elts + i*(size_t) array->sizes[dim]*dataTypeSize;
	}
    }

    alloc_multidim_array_debug();
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_multidim_array
@INPUT      : array
              n_dimensions
              sizes
              data_type
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a multidimensional array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void   create_multidim_array(
    VIO_multidim_array  *array,
    int             n_dimensions,
    int             sizes[],
    VIO_Data_types      data_type )
{
    create_empty_multidim_array( array, n_dimensions, data_type );
    set_multidim_sizes( array, sizes );
    alloc_multidim_array( array );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_multidim_array
@INPUT      : array
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the multidimensional array.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  delete_multidim_array(
    VIO_multidim_array   *array )
{
    if( array->data == NULL )
    {
        fprintf(stderr, "Cannot free NULL multidim data.\n" );
	exit(1);
        return;
    }

    free( array->data );

    array->data = NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_multidim_n_dimensions
@INPUT      : array
@OUTPUT     : 
@RETURNS    : number of dimensions
@DESCRIPTION: Returns the number of dimensions of the array
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June, 1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

int  get_multidim_n_dimensions(
    VIO_multidim_array   *array )
{
    return( array->n_dimensions );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_multidim_data_reordered
@INPUT      : type_size
              void_dest_ptr
              n_dest_dims
              dest_sizes
              void_src_ptr
              n_src_dims
              src_sizes
              counts
              to_dest_index
              use_src_order   - whether to step through arrays in the
                                reverse order of src or dest
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Copies any type of multidimensional data from the src array
              to the destination array.  to_dest_index is a lookup that
              converts src indices to destination indices, to allow arbitrary
              reordering of array data.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : Feb. 27, 1996   D. MacDonald  - made more efficient
---------------------------------------------------------------------------- */

void  copy_multidim_data_reordered(
    int                 type_size,
    void                *void_dest_ptr,
    int                 n_dest_dims,
    int                 dest_sizes[],
    void                *void_src_ptr,
    int                 n_src_dims,
    int                 src_sizes[],
    int                 counts[],
    int                 to_dest_index[],
    BOOLEAN             use_src_order )
{
    char      *src_ptr, *dest_ptr;
    int       d;
    int       dest_offsets[VIO_MAX_DIMENSIONS], src_offsets[VIO_MAX_DIMENSIONS];
    int       dest_offset0, dest_offset1, dest_offset2, dest_offset3;
    int       dest_offset4;
    int       src_offset0, src_offset1, src_offset2, src_offset3;
    int       src_offset4;
    int       dest_steps[VIO_MAX_DIMENSIONS], src_steps[VIO_MAX_DIMENSIONS];
    int       dest_index;
    int       n_transfer_dims;
    int       src_axis[VIO_MAX_DIMENSIONS], dest_axis[VIO_MAX_DIMENSIONS];
    int       transfer_counts[VIO_MAX_DIMENSIONS];
    int       v0, v1, v2, v3, v4;
    int       size0, size1, size2, size3, size4;
    BOOLEAN   full_count_used;

    /*--- initialize dest */

    dest_ptr = (char *) void_dest_ptr;
    dest_steps[n_dest_dims-1] = type_size;
    for( d = n_dest_dims-2 ; d >= 0; d-- )
        dest_steps[d] = dest_steps[d+1] * dest_sizes[d+1];

    /*--- initialize src */

    src_ptr = (char *) void_src_ptr;
    src_steps[n_src_dims-1] = type_size;
    for( d = n_src_dims-2; d >= 0; d-- )
        src_steps[d] = src_steps[d+1] * src_sizes[d+1];

    n_transfer_dims = 0;

    if( getenv( "VOLUME_IO_SRC_ORDER" ) )
        use_src_order = true;
    else if( getenv( "VOLUME_IO_DEST_ORDER" ) )
        use_src_order = false;

    if( use_src_order )
    {
        for_less( d, 0, n_src_dims )
        {
            dest_index = to_dest_index[d];
            if( dest_index >= 0 )
            {
                src_axis[n_transfer_dims] = d;
                dest_axis[n_transfer_dims] = dest_index;
                src_offsets[n_transfer_dims] = src_steps[d];
                dest_offsets[n_transfer_dims] = dest_steps[dest_index];
                transfer_counts[n_transfer_dims] = counts[d];
                ++n_transfer_dims;
            }
        }
    }
    else
    {
        for_less( dest_index, 0, n_dest_dims )
        {
            for_less( d, 0, n_src_dims )
                if( to_dest_index[d] == dest_index )
                    break;

            if( d < n_src_dims )
            {
                src_axis[n_transfer_dims] = d;
                dest_axis[n_transfer_dims] = dest_index;
                src_offsets[n_transfer_dims] = src_steps[d];
                dest_offsets[n_transfer_dims] = dest_steps[dest_index];
                transfer_counts[n_transfer_dims] = counts[d];
                ++n_transfer_dims;
            }
        }
    }

    /*--- check if we can transfer more than one at once */

    full_count_used = true;

    while( n_transfer_dims > 0 &&
           src_axis[n_transfer_dims-1] == n_src_dims-1 &&
           dest_axis[n_transfer_dims-1] == n_dest_dims-1 && full_count_used )
    {
        if( transfer_counts[n_transfer_dims-1] != src_sizes[n_src_dims-1] ||
            transfer_counts[n_transfer_dims-1] != dest_sizes[n_dest_dims-1] )
        {
            full_count_used = false;
        }

        type_size *= transfer_counts[n_transfer_dims-1];
        --n_src_dims;
        --n_dest_dims;
        --n_transfer_dims;
    }

    for_less( d, 0, n_transfer_dims-1 )
    {
        src_offsets[d] -= src_offsets[d+1] * transfer_counts[d+1];
        dest_offsets[d] -= dest_offsets[d+1] * transfer_counts[d+1];
    }

    /*--- slide the transfer dims to the last of the 5 dimensions */

    for( d = n_transfer_dims-1; d >= 0; d-- )
    {
        src_offsets    [d+VIO_MAX_DIMENSIONS-n_transfer_dims] = src_offsets    [d];
        dest_offsets   [d+VIO_MAX_DIMENSIONS-n_transfer_dims] = dest_offsets   [d];
        transfer_counts[d+VIO_MAX_DIMENSIONS-n_transfer_dims] = transfer_counts[d];
    }

    for_less( d, 0, VIO_MAX_DIMENSIONS-n_transfer_dims )
    {
        transfer_counts[d] = 1;
        src_offsets[d] = 0;
        dest_offsets[d] = 0;
    }

    size0 = transfer_counts[0];
    size1 = transfer_counts[1];
    size2 = transfer_counts[2];
    size3 = transfer_counts[3];
    size4 = transfer_counts[4];

    src_offset0 = src_offsets[0];
    src_offset1 = src_offsets[1];
    src_offset2 = src_offsets[2];
    src_offset3 = src_offsets[3];
    src_offset4 = src_offsets[4];

    dest_offset0 = dest_offsets[0];
    dest_offset1 = dest_offsets[1];
    dest_offset2 = dest_offsets[2];
    dest_offset3 = dest_offsets[3];
    dest_offset4 = dest_offsets[4];

    for_less( v0, 0, size0 )
    {
        for_less( v1, 0, size1 )
        {
            for_less( v2, 0, size2 )
            {
                for_less( v3, 0, size3 )
                {
                    for_less( v4, 0, size4 )
                    {
                        (void) memcpy( dest_ptr, src_ptr, (size_t) type_size );
                        src_ptr += src_offset4;
                        dest_ptr += dest_offset4;
                    }
                    src_ptr += src_offset3;
                    dest_ptr += dest_offset3;
                }
                src_ptr += src_offset2;
                dest_ptr += dest_offset2;
            }
            src_ptr += src_offset1;
            dest_ptr += dest_offset1;
        }
        src_ptr += src_offset0;
        dest_ptr += dest_offset0;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_multidim_reordered
@INPUT      : dest
              dest_ind
              src
              src_ind
              counts
              to_dest_index
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Copies data from src array to dest array, with dimension
              translation given by to_dest_index[].
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 1, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

void  copy_multidim_reordered(
    VIO_multidim_array      *dest,
    int                 dest_ind[],
    VIO_multidim_array      *src,
    int                 src_ind[],
    int                 counts[],
    int                 to_dest_index[] )
{
    int       n_src_dims, n_dest_dims, type_size;
    int       src_sizes[VIO_MAX_DIMENSIONS], dest_sizes[VIO_MAX_DIMENSIONS];
    char      *dest_ptr, *src_ptr;
    void      *void_ptr;

    type_size = get_type_size( get_multidim_data_type(dest) );

    /*--- initialize dest */

    n_dest_dims = get_multidim_n_dimensions( dest );
    get_multidim_sizes( dest, dest_sizes );
    GET_MULTIDIM_PTR( void_ptr, *dest, dest_ind[0], dest_ind[1], dest_ind[2],
                      dest_ind[3], dest_ind[4] );
    dest_ptr = void_ptr;

    /*--- initialize src */

    n_src_dims = get_multidim_n_dimensions( src );
    get_multidim_sizes( src, src_sizes );
    GET_MULTIDIM_PTR( void_ptr, *src, src_ind[0], src_ind[1], src_ind[2],
                      src_ind[3], src_ind[4] );
    src_ptr = void_ptr;

    copy_multidim_data_reordered( type_size,
                                  dest_ptr, n_dest_dims, dest_sizes,
                                  src_ptr, n_src_dims, src_sizes,
                                  counts, to_dest_index, true );
}


#ifdef TEST_MULTIDIM_ARRAY_C
int main() {
    int n_dimensions = 5;

    int  sizes[5] = {2,3,5,7,11};
    long limit=1;
    {
        int i;
	for (i = 0; i < 5; i++) limit *= sizes[i];
    }
    
    VIO_multidim_array array;
    create_multidim_array(&array, n_dimensions, sizes, UNSIGNED_INT);

    int mode;
    for (mode = 0; mode<2; mode++) {
       long showLimit = 1;
       long i;
       size_t errors = 0;
       for (i = 0; (i < limit) && (errors < 10); i++) {
           int indices[5];
	   long tmp = i;
	   int j;
	   for (j = 5; j > 0;) {
	       j--;
	       indices[j] = tmp % sizes[j];
	       tmp             /= sizes[j];
	   }
	   void* voidPtr;  GET_MULTIDIM_PTR( voidPtr, array, indices[0], indices[1], indices[2], indices[3], indices[4] );
	   unsigned int* ptr = (unsigned int*)voidPtr;
	   if (!mode) *ptr = i; 
	   else {
	       bool show = (i < 32) || (i >= showLimit-1);
	       if (i != *ptr) {
	       	   errors++;
	           printf("Error expected %ld != gotten %ld pointer %p at ", (long)i, (long)*ptr, ptr);
	       } else if (i >= showLimit-1) {
	       	   showLimit *= 2;
	       }
	       if (show) {
	          printf("i:%ld ", i);
	          const char* sep="["; 
	          for (j = 0; j < 5; j++) { printf("%s%d", sep, indices[j]); sep = ","; }
		  printf("]\n");
	       }
	   }
       }
    }  
    return 0;
}
#endif
