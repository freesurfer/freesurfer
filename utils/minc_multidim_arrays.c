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

#include  <internal_volume_io.h>
#include  <limits.h>
#include  <float.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/Volumes/multidim_arrays.c,v 1.15.2.1 2004/10/04 20:20:14 bert Exp $";
#endif

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

VIOAPI   void   create_empty_multidim_array(
    multidim_array  *array,
    int             n_dimensions,
    Data_types      data_type )
{
    if( n_dimensions < 1 || n_dimensions > MAX_DIMENSIONS )
    {
        print_error(
     "create_empty_multidim_array(): n_dimensions (%d) not in range 1 to %d.\n",
               n_dimensions, MAX_DIMENSIONS );
    }

    array->n_dimensions = n_dimensions;
    array->data_type = data_type;
    array->data = (void *) NULL;
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

VIOAPI  Data_types  get_multidim_data_type(
    multidim_array       *array )
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

VIOAPI  void  set_multidim_data_type(
    multidim_array       *array,
    Data_types           data_type )
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

VIOAPI  int  get_type_size(
    Data_types   type )
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
    }

    return( size );
}

VIOAPI  void  get_type_range(
    Data_types   type,
    Real         *min_value,
    Real         *max_value )
{
    switch( type )
    {
    case UNSIGNED_BYTE:
        *min_value = 0.0;
        *max_value = (Real) UCHAR_MAX;     break;
    case SIGNED_BYTE:
        *min_value = (Real) SCHAR_MIN;
        *max_value = (Real) SCHAR_MAX;     break;
    case UNSIGNED_SHORT:
        *min_value = 0.0;
        *max_value = (Real) USHRT_MAX;     break;
    case SIGNED_SHORT:
        *min_value = (Real) SHRT_MIN;
        *max_value = (Real) SHRT_MAX;      break;
    case UNSIGNED_INT:
        *min_value = 0.0;
        *max_value = (Real) UINT_MAX;     break;
    case SIGNED_INT:
        *min_value = (Real) INT_MIN;
        *max_value = (Real) INT_MAX;      break;
    case FLOAT:
        *min_value = (Real) -FLT_MAX;
        *max_value = (Real) FLT_MAX;       break;
    case DOUBLE:
        *min_value = (Real) -DBL_MAX;
        *max_value = (Real) DBL_MAX;       break;
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

VIOAPI  void  set_multidim_sizes(
    multidim_array   *array,
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

VIOAPI  void  get_multidim_sizes(
    multidim_array   *array,
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

VIOAPI  BOOLEAN  multidim_array_is_alloced(
    multidim_array   *array )
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

VIOAPI  void  alloc_multidim_array(
    multidim_array   *array )
{
    int     dim;
    size_t  type_size, sizes[5];
    void    *p1, **p2, ***p3, ****p4, *****p5;

    if( multidim_array_is_alloced( array ) )
        delete_multidim_array( array );

    if( array->data_type == NO_DATA_TYPE )
    {
        print_error(
           "Error: cannot allocate array data until size specified.\n" );
        return;
    }

    for_less( dim, 0, array->n_dimensions )
        sizes[dim] = (size_t) array->sizes[dim];

    type_size = (size_t) get_type_size( array->data_type );

    switch( array->n_dimensions )
    {
    case  1:
        ASSIGN_PTR(p1) = alloc_memory_1d( sizes[0], type_size
                                          _ALLOC_SOURCE_LINE );
        array->data = (void *) p1;
        break;
    case  2:
        ASSIGN_PTR(p2) = alloc_memory_2d( sizes[0], sizes[1], type_size
                                          _ALLOC_SOURCE_LINE);
        array->data = (void *) p2;
        break;
    case  3:
        ASSIGN_PTR(p3) = alloc_memory_3d( sizes[0], sizes[1], sizes[2],
                                          type_size _ALLOC_SOURCE_LINE );
        array->data = (void *) p3;
        break;
    case  4:
        ASSIGN_PTR(p4) = alloc_memory_4d( sizes[0], sizes[1],
                             sizes[2], sizes[3], type_size _ALLOC_SOURCE_LINE );
        array->data = (void *) p4;
        break;
    case  5:
        ASSIGN_PTR(p5) = alloc_memory_5d( sizes[0], sizes[1],
                              sizes[2], sizes[3], sizes[4], type_size
                              _ALLOC_SOURCE_LINE );
        array->data = (void *) p5;
        break;
    }
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

VIOAPI   void   create_multidim_array(
    multidim_array  *array,
    int             n_dimensions,
    int             sizes[],
    Data_types      data_type )
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

VIOAPI  void  delete_multidim_array(
    multidim_array   *array )
{
    if( array->data == NULL )
    {
        print_error( "Warning: cannot free NULL multidim data.\n" );
        return;
    }

    switch( array->n_dimensions )
    {
    case  1:  free_memory_1d( (void **) &array->data _ALLOC_SOURCE_LINE );
              break;
    case  2:  free_memory_2d( (void ***) &array->data _ALLOC_SOURCE_LINE );
              break;
    case  3:  free_memory_3d( (void ****) &array->data _ALLOC_SOURCE_LINE );
              break;
    case  4:  free_memory_4d( (void *****) &array->data _ALLOC_SOURCE_LINE );
              break;
    case  5:  free_memory_5d( (void ******) &array->data _ALLOC_SOURCE_LINE );
              break;
    }

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

VIOAPI  int  get_multidim_n_dimensions(
    multidim_array   *array )
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

VIOAPI  void  copy_multidim_data_reordered(
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
    int       dest_offsets[MAX_DIMENSIONS], src_offsets[MAX_DIMENSIONS];
    int       dest_offset0, dest_offset1, dest_offset2, dest_offset3;
    int       dest_offset4;
    int       src_offset0, src_offset1, src_offset2, src_offset3;
    int       src_offset4;
    int       dest_steps[MAX_DIMENSIONS], src_steps[MAX_DIMENSIONS];
    int       dest_index;
    int       n_transfer_dims;
    int       src_axis[MAX_DIMENSIONS], dest_axis[MAX_DIMENSIONS];
    int       transfer_counts[MAX_DIMENSIONS];
    int       v0, v1, v2, v3, v4;
    int       size0, size1, size2, size3, size4;
    BOOLEAN   full_count_used;

    /*--- initialize dest */

    dest_ptr = (char *) void_dest_ptr;
    dest_steps[n_dest_dims-1] = type_size;
    for_down( d, n_dest_dims-2, 0 )
        dest_steps[d] = dest_steps[d+1] * dest_sizes[d+1];

    /*--- initialize src */

    src_ptr = (char *) void_src_ptr;
    src_steps[n_src_dims-1] = type_size;
    for_down( d, n_src_dims-2, 0 )
        src_steps[d] = src_steps[d+1] * src_sizes[d+1];

    n_transfer_dims = 0;

    if( getenv( "VOLUME_IO_SRC_ORDER" ) )
        use_src_order = TRUE;
    else if( getenv( "VOLUME_IO_DEST_ORDER" ) )
        use_src_order = FALSE;

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

    full_count_used = TRUE;

    while( n_transfer_dims > 0 &&
           src_axis[n_transfer_dims-1] == n_src_dims-1 &&
           dest_axis[n_transfer_dims-1] == n_dest_dims-1 && full_count_used )
    {
        if( transfer_counts[n_transfer_dims-1] != src_sizes[n_src_dims-1] ||
            transfer_counts[n_transfer_dims-1] != dest_sizes[n_dest_dims-1] )
        {
            full_count_used = FALSE;
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

    for_down( d, n_transfer_dims-1, 0 )
    {
        src_offsets[d+MAX_DIMENSIONS-n_transfer_dims] = src_offsets[d];
        dest_offsets[d+MAX_DIMENSIONS-n_transfer_dims] = dest_offsets[d];
        transfer_counts[d+MAX_DIMENSIONS-n_transfer_dims] = transfer_counts[d];
    }

    for_less( d, 0, MAX_DIMENSIONS-n_transfer_dims )
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

VIOAPI  void  copy_multidim_reordered(
    multidim_array      *dest,
    int                 dest_ind[],
    multidim_array      *src,
    int                 src_ind[],
    int                 counts[],
    int                 to_dest_index[] )
{
    int       n_src_dims, n_dest_dims, type_size;
    int       src_sizes[MAX_DIMENSIONS], dest_sizes[MAX_DIMENSIONS];
    char      *dest_ptr, *src_ptr;
    void      *void_ptr;

    type_size = get_type_size( get_multidim_data_type(dest) );

    /*--- initialize dest */

    n_dest_dims = get_multidim_n_dimensions( dest );
    get_multidim_sizes( dest, dest_sizes );
    GET_MULTIDIM_PTR( void_ptr, *dest, dest_ind[0], dest_ind[1], dest_ind[2],
                      dest_ind[3], dest_ind[4] );
    ASSIGN_PTR( dest_ptr ) = void_ptr;

    /*--- initialize src */

    n_src_dims = get_multidim_n_dimensions( src );
    get_multidim_sizes( src, src_sizes );
    GET_MULTIDIM_PTR( void_ptr, *src, src_ind[0], src_ind[1], src_ind[2],
                      src_ind[3], src_ind[4] );
    ASSIGN_PTR( src_ptr ) = void_ptr;

    copy_multidim_data_reordered( type_size,
                                  dest_ptr, n_dest_dims, dest_sizes,
                                  src_ptr, n_src_dims, src_sizes,
                                  counts, to_dest_index, TRUE );
}
