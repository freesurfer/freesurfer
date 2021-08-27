#include "minc_volume_io.h"
#include "minc_internals.h"

/**
 * @brief substitutes for the needed functionality previously obtained from minc
 */
/*
 * Original Author: Bevin Brett
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

/*------------------------------------------------------------------------
  HEADERS
  ------------------------------------------------------------------------*/

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "minc_internals.h"
#include "minc_files.h"


static VIO_Status  input_transform(
    FILE                *file,
    const char*          filename,
    General_transform   *transform );
    

char* statusToString(VIO_Status status) {
    static char* table[] = {
               "OK",
               "ERROR",
               "INTERNAL_ERROR",
               "END_OF_FILE",
               "QUIT"};
    return status >= 0 && status < VIO_Status__end 
        ? table[status] 
	: "Bad status";
};


// Based on minc-1.5.1/volume_io/Geometry/colour.c
// which requires the following...
/* ----------------------------------------------------------------------------
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
static VIO_Colour make_rgba_Colour(
    int    r,
    int    g,
    int    b,
    int    a )
{
    return 
       (((VIO_Colour)(unsigned char)a) <<  0)
     + (((VIO_Colour)(unsigned char)b) <<  8)
     + (((VIO_Colour)(unsigned char)g) << 16)
     + (((VIO_Colour)(unsigned char)r) << 24);
}


VIO_Colour  make_rgba_Colour_0_1(
    double   r,
    double   g,
    double   b,
    double   a )
{
    return( make_rgba_Colour( (int) (r * 255.0 + 0.5),
                              (int) (g * 255.0 + 0.5),
                              (int) (b * 255.0 + 0.5),
                              (int) (a * 255.0 + 0.5) ) );
}


// Based on minc-1.5.1/volume_io/Geometry/transforms.c
// which requires the following...
//
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
void  make_identity_transform( Transform   *transform )
{
    Transform_elem( *transform, 0, 0 ) = 1.0;
    Transform_elem( *transform, 0, 1 ) = 0.0;
    Transform_elem( *transform, 0, 2 ) = 0.0;
    Transform_elem( *transform, 0, 3 ) = 0.0;
    Transform_elem( *transform, 1, 0 ) = 0.0;
    Transform_elem( *transform, 1, 1 ) = 1.0;
    Transform_elem( *transform, 1, 2 ) = 0.0;
    Transform_elem( *transform, 1, 3 ) = 0.0;
    Transform_elem( *transform, 2, 0 ) = 0.0;
    Transform_elem( *transform, 2, 1 ) = 0.0;
    Transform_elem( *transform, 2, 2 ) = 1.0;
    Transform_elem( *transform, 2, 3 ) = 0.0;
    Transform_elem( *transform, 3, 0 ) = 0.0;
    Transform_elem( *transform, 3, 1 ) = 0.0;
    Transform_elem( *transform, 3, 2 ) = 0.0;
    Transform_elem( *transform, 3, 3 ) = 1.0;
}

bool close_to_identity(
    Transform   *transform )
{
    int i,j;
    for( i = 0; i < 4; i++ )
    {
        for( j = 0; j < 4; j++ )
        {
            double expected = ( i == j ) ? 1.0 : 0.0;
	    double actual   = Transform_elem(*transform,i,j);
            if( fabsf(expected - actual) > 0.001 ) return false;
        }
    }

    return( true );
}

void   concat_transforms(
    Transform   *result,
    Transform   *t1,
    Transform   *t2 )
{
    bool        result_is_also_an_arg;
    Transform   tmp, *t;

    /*--- check if the result transform is same as one of the arguments */

    if( result == t1 || result == t2 )
    {
        result_is_also_an_arg = true;
        t = &tmp;
    }
    else
    {
        result_is_also_an_arg = false;
        t = result;
    }

    /*--- perform multiplication */
    int i;
    for( i=0; i < 4; i++ )
    {
        int j;
        for( j=0; j < 4; j++ )
        {
            double sum = 0.0;
	    int k;
            for( k=0; k < 4; k++ )
            {
                sum += Transform_elem( *t2, i, k ) *
                       Transform_elem( *t1, k, j );
            }
            Transform_elem( *t, i, j ) = sum;
        }
    }

    if( result_is_also_an_arg )
        *result = tmp;
}

static void  homogenous_transform_point(
    Transform  *transform,
    double     x,
    double     y,
    double     z,
    double     w,
    double     *x_trans,
    double     *y_trans,
    double     *z_trans )
{
    double     w_trans;

    *x_trans = Transform_elem(*transform,0,0) * x +
               Transform_elem(*transform,0,1) * y +
               Transform_elem(*transform,0,2) * z +
               Transform_elem(*transform,0,3) * w;

    *y_trans = Transform_elem(*transform,1,0) * x +
               Transform_elem(*transform,1,1) * y +
               Transform_elem(*transform,1,2) * z +
               Transform_elem(*transform,1,3) * w;

    *z_trans = Transform_elem(*transform,2,0) * x +
               Transform_elem(*transform,2,1) * y +
               Transform_elem(*transform,2,2) * z +
               Transform_elem(*transform,2,3) * w;

    w_trans =  Transform_elem(*transform,3,0) * x +
               Transform_elem(*transform,3,1) * y +
               Transform_elem(*transform,3,2) * z +
               Transform_elem(*transform,3,3) * w;

    if( w_trans != 0.0 && w_trans != 1.0 )
    {
        *x_trans /= w_trans;
        *y_trans /= w_trans;
        *z_trans /= w_trans;
    }
}

void  transform_point(
    Transform  *transform,
    double     x,
    double     y,
    double     z,
    double     *x_trans,
    double     *y_trans,
    double     *z_trans )
{
    homogenous_transform_point( transform, x, y, z, 1.0,
                                x_trans, y_trans, z_trans );
}

void  get_transform_origin_real(
    Transform   *transform,
    double       origin[] )
{
    origin[X] = Transform_elem(*transform,0,3);
    origin[Y] = Transform_elem(*transform,1,3);
    origin[Z] = Transform_elem(*transform,2,3);
}

void  get_transform_x_axis_real(
    Transform   *transform,
    double       x_axis[] )
{
    x_axis[X] = Transform_elem(*transform,0,0);
    x_axis[Y] = Transform_elem(*transform,1,0);
    x_axis[Z] = Transform_elem(*transform,2,0);
}

void  get_transform_y_axis_real(
    Transform   *transform,
    double       y_axis[] )
{
    y_axis[X] = Transform_elem(*transform,0,1);
    y_axis[Y] = Transform_elem(*transform,1,1);
    y_axis[Z] = Transform_elem(*transform,2,1);
}

void  get_transform_z_axis_real(
    Transform   *transform,
    double       z_axis[] )
{
    z_axis[X] = Transform_elem(*transform,0,2);
    z_axis[Y] = Transform_elem(*transform,1,2);
    z_axis[Z] = Transform_elem(*transform,2,2);
}


// Based on minc-1.5.1/volume_io/Geometry/inverse.c
// which requires the following...
//
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
static bool compute_transform_inverse(
    Transform  *transform,
    Transform  *inverse )
{
    int        i, j;

    Double4x4  t, inv;

    for( i = 0; i < 4; i++ )
    {
        for( j = 0; j < 4; j++ )
            t[Index4x4(i,j)] = Transform_elem(*transform,i,j);
    }

    bool success = invert_4x4_matrix( (Double4x4 const*)&t, &inv );

    if( success )
    {
        /* --- copy the resulting numerical recipes matrix to the
               output argument */

        for( i = 0; i < 4; i++ )
        {
            for( j = 0; j < 4; j++ )
            {
                Transform_elem(*inverse,i,j) = inv[Index4x4(i,j)];
            }
        }

#ifdef  DEBUG
        /* --- check if this really is an inverse, by multiplying */

        {
            Transform  ident;

            concat_transforms( &ident, transform, inverse );

            if( !close_to_identity(&ident) )
            {
                print_error( "Error in compute_transform_inverse\n" );
            }
        }
#endif
    }
    else
        make_identity_transform( inverse );

    return( success );
}


// Based on minc-1.5.1/volume_io/MNI_formats/gen_xfs.c
// which requires the following...
//
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
int  get_n_concated_transforms(
    General_transform   *transform )
{
    if( transform->type == CONCATENATED_TRANSFORM )
        return( transform->n_transforms );
    else
        return( 1 );
}

General_transform  *get_nth_general_transform(
    General_transform   *transform,
    int                 n )
{
    if( n < 0 || n >= get_n_concated_transforms( transform ) )
    {
        fprintf(stderr, "get_nth_general_transform bad n\n" );
	exit(1);
        return( (General_transform *) NULL );
    }
    else if( transform->type == CONCATENATED_TRANSFORM )
        return( &transform->transforms[n] );
    else
        return( transform );
}

Transform  *get_linear_transform_ptr(
    General_transform   *transform )
{
    if( transform->type == LINEAR )
    {
        if( transform->inverse_flag )
            return( transform->inverse_linear_transform );
        else
            return( transform->linear_transform );
    }
    else
    {
        fprintf(stderr, "get_linear_transform_ptr\n" );
	exit(1);
        return( (Transform *) NULL );
    }
}

Transform  *get_inverse_linear_transform_ptr(
    General_transform   *transform )
{
    if( transform->type == LINEAR )
    {
        if( transform->inverse_flag )
            return( transform->linear_transform );
        else
            return( transform->inverse_linear_transform );
    }
    else
    {
        fprintf(stderr, "get_inverse_linear_transform_ptr\n" );
	exit(1);
        return( (Transform *) NULL );
    }
}

static void  general_transform_point(
    General_transform   *transform,
    double              x,
    double              y,
    double              z,
    double              *x_transformed,
    double              *y_transformed,
    double              *z_transformed );
    
static void  general_inverse_transform_point(
    General_transform   *transform,
    double              x,
    double              y,
    double              z,
    double              *x_transformed,
    double              *y_transformed,
    double              *z_transformed );
    
static  void  transform_or_invert_point(
    General_transform   *transform,
    bool                inverse_flag,
    double              x,
    double              y,
    double              z,
    double              *x_transformed,
    double              *y_transformed,
    double              *z_transformed )
{
    switch( transform->type )
    {
    case LINEAR:
        if( inverse_flag )
            transform_point( transform->inverse_linear_transform,
                             x, y, z,
                             x_transformed, y_transformed, z_transformed );
        else
            transform_point( transform->linear_transform,
                             x, y, z,
                             x_transformed, y_transformed, z_transformed );
        break;

#if defined(BEVIN_ALL_TYPES_SUPPORTED)
    case THIN_PLATE_SPLINE:
        if( inverse_flag )
        {
            thin_plate_spline_inverse_transform( transform->n_dimensions,
                                                 transform->n_points,
                                                 transform->points,
                                                 transform->displacements,
                                                 x, y, z,
                                                 x_transformed, y_transformed,
                                                 z_transformed );
        }
        else
        {
            thin_plate_spline_transform( transform->n_dimensions,
                                         transform->n_points,
                                         transform->points,
                                         transform->displacements,
                                         x, y, z,
                                         x_transformed, y_transformed,
                                         z_transformed );
        }
        break;

    case GRID_TRANSFORM:
        if( inverse_flag )
        {
            grid_inverse_transform_point( transform,
                                          x, y, z,
                                          x_transformed, y_transformed,
                                          z_transformed );
        }
        else
        {
            grid_transform_point( transform,
                                  x, y, z,
                                  x_transformed, y_transformed,
                                  z_transformed );
        }
        break;

    case USER_TRANSFORM:
        if( inverse_flag )
        {
            transform->user_inverse_transform_function(
                           transform->user_data, x, y, z,
                           x_transformed, y_transformed, z_transformed );
        }
        else
        {
            transform->user_transform_function(
                           transform->user_data, x, y, z,
                           x_transformed, y_transformed, z_transformed );
        }
        break;
#endif

    case CONCATENATED_TRANSFORM:
        *x_transformed = x;
        *y_transformed = y;
        *z_transformed = z;

        if( inverse_flag )
        {
    	    int trans;
            for( trans = transform->n_transforms-1;  trans >= 0;  --trans )
            {
                general_inverse_transform_point( &transform->transforms[trans],
                             *x_transformed, *y_transformed, *z_transformed,
                             x_transformed, y_transformed, z_transformed );
            }
        }
        else
        {
    	    int trans;
            for( trans = 0; trans < transform->n_transforms; trans++ )
            {
                general_transform_point( &transform->transforms[trans],
                             *x_transformed, *y_transformed, *z_transformed,
                             x_transformed, y_transformed, z_transformed );
            }
        }
        break;

    default:
        fprintf(stderr, "%s:%d transform_or_invert_point type NYI\n", __FILE__,__LINE__ );
	exit(1);
        break;
    }
}

static void  general_transform_point(
    General_transform   *transform,
    double              x,
    double              y,
    double              z,
    double              *x_transformed,
    double              *y_transformed,
    double              *z_transformed )
{

    transform_or_invert_point( transform, transform->inverse_flag, x, y, z,
                               x_transformed, y_transformed, z_transformed );
}

static void  general_inverse_transform_point(
    General_transform   *transform,
    double              x,
    double              y,
    double              z,
    double              *x_transformed,
    double              *y_transformed,
    double              *z_transformed )
{

    transform_or_invert_point( transform, !transform->inverse_flag, x, y, z,
                               x_transformed, y_transformed, z_transformed );
}

static  void  alloc_linear_transform(
    General_transform   *transform )
{
    transform->type = LINEAR;
    transform->inverse_flag = false;

    transform->linear_transform         = (Transform*)malloc(sizeof(Transform));
    transform->inverse_linear_transform = (Transform*)malloc(sizeof(Transform));
}

static void create_linear_transform(
    General_transform   *transform,
    Transform           *linear_transform )
{
    alloc_linear_transform( transform );

    if( linear_transform != (Transform *) NULL &&
        compute_transform_inverse( linear_transform,
                                   transform->inverse_linear_transform ) )
    {
        *(transform->linear_transform) = *linear_transform;
    }
    else
    {
        make_identity_transform( transform->linear_transform );
        make_identity_transform( transform->inverse_linear_transform );
    }
}


static  void  copy_and_invert_transform(
    General_transform  *transform,
    bool		invert_it,
    General_transform  *copy )
{
    Transform      *swap;
    int            trans;

    *copy = *transform;

    switch( transform->type )
    {
    case LINEAR:
        alloc_linear_transform( copy );
        *(copy->linear_transform) = *(transform->linear_transform);
        *(copy->inverse_linear_transform) =
                                       *(transform->inverse_linear_transform);

        if( transform->inverse_flag )
            invert_it = !invert_it;

        if( invert_it )
        {
            swap = copy->linear_transform;
            copy->linear_transform = copy->inverse_linear_transform;
            copy->inverse_linear_transform = swap;
        }
        copy->inverse_flag = false;
        break;

#if defined(BEVIN_ALL_TYPES_SUPPORTED)
    case THIN_PLATE_SPLINE:
        ALLOC2D( copy->points, copy->n_points, copy->n_dimensions);
        ALLOC2D( copy->displacements, copy->n_points + copy->n_dimensions + 1,
                 copy->n_dimensions);

        for_less( i, 0, copy->n_points )
            for_less( j, 0, copy->n_dimensions )
                copy->points[i][j] = transform->points[i][j];

        for_less( i, 0, copy->n_points + copy->n_dimensions + 1 )
            for_less( j, 0, copy->n_dimensions )
                copy->displacements[i][j] = transform->displacements[i][j];

        if( invert_it )
            copy->inverse_flag = !copy->inverse_flag;
        break;

    case GRID_TRANSFORM:
        copy->displacement_volume = (void *) copy_volume(
                                    (Volume) transform->displacement_volume );

        if( invert_it )
            copy->inverse_flag = !copy->inverse_flag;

        break;

    case USER_TRANSFORM:
        ALLOC( byte_ptr, copy->size_user_data );
        copy->user_data = byte_ptr;
        (void) memcpy( copy->user_data, transform->user_data,
                       copy->size_user_data );
        if( invert_it )
            copy->inverse_flag = !copy->inverse_flag;
        break;
#endif

    case CONCATENATED_TRANSFORM:
        copy->transforms = (General_transform*)malloc(
		sizeof(General_transform)*copy->n_transforms);
        for( trans=0; trans < copy->n_transforms; trans++ )
        {
            copy_general_transform( &transform->transforms[trans],
                                    &copy->transforms[trans] );
        }
        if( invert_it )
            copy->inverse_flag = !copy->inverse_flag;
        break;

    default:
        fprintf(stderr, "copy_and_invert_transform bad kind\n" );
	exit(1);
        break;
    }
}

static void create_inverse_general_transform(
    General_transform   *transform,
    General_transform   *inverse )
{
    copy_and_invert_transform( transform, true, inverse );
}

void  copy_general_transform(
    General_transform   *transform,
    General_transform   *copy )
{
    copy_and_invert_transform( transform, false, copy );
}

void  delete_general_transform(
    General_transform   *transform )
{

    switch( transform->type )
    {
    case LINEAR:
        free( transform->linear_transform );
        free( transform->inverse_linear_transform );
        break;
#if defined(BEVIN_ALL_TYPES_SUPPORTED)
    case THIN_PLATE_SPLINE:
        if( transform->n_points > 0 && transform->n_dimensions > 0 )
        {
            FREE2D( transform->points );
            FREE2D( transform->displacements );
        }
        break;

    case GRID_TRANSFORM:
        delete_volume( (Volume) transform->displacement_volume );
        break;

    case USER_TRANSFORM:
        if( transform->size_user_data )
            FREE( transform->user_data );
        break;
#endif
    case CONCATENATED_TRANSFORM: {
        int   trans;
        for( trans = 0; transform->n_transforms > 0; trans++)
            delete_general_transform( &transform->transforms[trans] );

        if( transform->n_transforms > 0 )
            free( transform->transforms );
    }   break;

    default:
        fprintf(stderr, "delete_general_transform\n" );
	exit(1);
        break;
    }
}


// Based on minc-1.5.1/volume_io/MNI_formats/mni_io.c
// which requires the following...
//
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
static const char COMMENT_CHAR1 = '%';
static const char COMMENT_CHAR2 = '#';

static VIO_Status  mni_get_nonwhite_character(
    FILE   *file,
    char   *ch )
{

    bool in_comment = false;

    VIO_Status status;
    do
    {
        status = input_character( file, ch );
        if( status == OK ) {
            if ( *ch == COMMENT_CHAR1 || *ch == COMMENT_CHAR2 ) 
	        in_comment = true;
	    else if( *ch == '\n' )
                in_comment = false;
	}
    }
    while( status == OK &&
           (in_comment || *ch == ' ' || *ch == '\t' || *ch == '\n' || 
            *ch == '\r') );     /* ignore carriage returns */

    if( status == ERROR )
        status = END_OF_FILE;

    return( status );
}

static VIO_Status  mni_skip_expected_character(
    FILE   *file,
    char   expected_ch )
{
    char       ch;
    VIO_Status status = mni_get_nonwhite_character( file, &ch );

    if( status == OK )
    {
        if( ch != expected_ch )
        {
            fprintf(stderr, "Expected '%c', found '%c'.\n", expected_ch, ch );
            status = ERROR;
        }
    }
    else
    {
        fprintf(stderr, "Expected '%c', found end of file.\n", expected_ch );
    }

    return( status );
}


static VIO_Status  mni_input_string(
    FILE  *file,
    char* *string,
    char   termination_char1,
    char   termination_char2 )
{
    *string = NULL;
    
    size_t buffer_capacity = 512;
    size_t buffer_size     = 1;		// trailing 0
    char*  buffer = (char*)malloc(buffer_capacity);
 
    char ch;
    VIO_Status status = mni_get_nonwhite_character( file, &ch );

    bool quoted = false;
    if( status == OK && ch == '"' )
    {
        quoted = true;
        status = mni_get_nonwhite_character( file, &ch );
        termination_char1 = '"';
        termination_char2 = '"';
    }

    while( status == OK &&
           ch != termination_char1 && ch != termination_char2 && ch != '\n' )
    {
        if (ch != '\r') {       /* Always ignore carriage returns */
	   buffer_size++;
           if ( buffer_capacity == buffer_size ) {
	       buffer_capacity += buffer_capacity/4;
	       buffer = (char*)realloc(buffer,buffer_capacity);
	   }
	   buffer[buffer_size - 2] = ch;
        }
        status = input_character( file, &ch );
    }

    if( !quoted )
        (void) unget_character( file, ch );

    while (buffer_size > 1 && buffer[buffer_size - 2] == ' ') buffer_size--;
    buffer[buffer_size - 1] = 0;
     
    if( status != OK )
    {
        free( buffer );
    } else { 
        *string = buffer;
    }

    return( status );
}

static VIO_Status  mni_input_double(
    FILE    *file,
    double  *d )
{
    char* str;
    VIO_Status status = mni_input_string( file, &str, (char) ' ', (char) ';' );

    if( status == OK && sscanf( str, "%lf", d ) != 1 )
    {
        fprintf(stderr, "%s is not a valid double\n", str);
        status = ERROR;
    }

    free( str );

    return( status );
}

static VIO_Status  mni_input_keyword_and_equal_sign(
    FILE         *file,
    const char*  keyword,
    bool         print_error_message )
{
    char* str;
    VIO_Status status = mni_input_string( file, &str, (char) '=', (char) 0 );

    if( status == END_OF_FILE )
        return( status );

    if( status != OK || strcmp( str, keyword ) ||
        mni_skip_expected_character( file, (char) '=' ) != OK )
    {
        if( print_error_message )
            fprintf(stderr, "Expected \"%s =\" got str:%s\n", keyword, str );
        status = ERROR;
    }

    free( str );

    return( status );
}


// Based on minc-1.5.1/volume_io/MNI_formats/gen_xf_io.c
// which requires the following...
//
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


void concat_general_transforms(
    General_transform   *first,
    General_transform   *second,
    General_transform   *result )
{
    int                  first_start, first_end, first_step;
    int                  second_start, second_end, second_step;
    int                  i, trans;
    BOOLEAN              crunching_linear;
    BOOLEAN              first_inverted_concat, second_inverted_concat;
    Transform            *first_transform, *first_inverse;
    Transform            *second_transform, *second_inverse;
    General_transform    result_tmp, *result_ptr;
    General_transform    *transform;

    if( result == first || result == second )
        result_ptr = &result_tmp;
    else
        result_ptr = result;
    

    first_inverted_concat = first->type == CONCATENATED_TRANSFORM &&
                            first->inverse_flag;
    second_inverted_concat = second->type == CONCATENATED_TRANSFORM &&
                             second->inverse_flag;

    if( first->inverse_flag )
    {
        first_start = get_n_concated_transforms( first ) - 1;
        first_end = 0;
        first_step = -1;
    }
    else
    {
        first_start = 0;
        first_end = get_n_concated_transforms( first ) - 1;
        first_step = 1;
    }

    if( second->inverse_flag )
    {
        second_start = get_n_concated_transforms( second ) - 1;
        second_end = 0;
        second_step = -1;
    }
    else
    {
        second_start = 0;
        second_end = get_n_concated_transforms( second ) - 1;
        second_step = 1;
    }

    result_ptr->n_transforms = fabs( first_end - first_start ) + 1 +
                               fabs( second_end - second_start ) + 1;

    crunching_linear = false;
    if( get_nth_general_transform( first, first_end )->type == LINEAR &&
        get_nth_general_transform( second, second_start )->type == LINEAR )
    {
        --result_ptr->n_transforms;
        crunching_linear = true;
        first_end -= first_step;
        second_start += second_step;
    }

    if( result_ptr->n_transforms == 1 )
        result_ptr->type = LINEAR;
    else
    {
        result_ptr->type = CONCATENATED_TRANSFORM;
        result_ptr->transforms = 
	    (General_transform*)malloc(
		sizeof(General_transform)*result_ptr->n_transforms );
    }

    result_ptr->inverse_flag = false;

    trans = 0;
    for( i = first_start;  i != first_end + first_step;  i += first_step )
    {
        copy_and_invert_transform( get_nth_general_transform( first, i ),
                                   first_inverted_concat,
                                   get_nth_general_transform(result_ptr,trans));
        ++trans;
    }

    if( crunching_linear )
    {
        transform = get_nth_general_transform( result_ptr, trans );
        alloc_linear_transform( transform );

        if( first_inverted_concat )
        {
            first_inverse = get_linear_transform_ptr(
                      get_nth_general_transform(first,first_end+first_step));
            first_transform = get_inverse_linear_transform_ptr(
                      get_nth_general_transform(first,first_end+first_step));
        }
        else
        {
            first_transform = get_linear_transform_ptr(
                      get_nth_general_transform(first,first_end+first_step));
            first_inverse = get_inverse_linear_transform_ptr(
                      get_nth_general_transform(first,first_end+first_step));
        }

        if( second_inverted_concat )
        {
            second_inverse = get_linear_transform_ptr(
                   get_nth_general_transform(second,second_start-second_step));
            second_transform = get_inverse_linear_transform_ptr(
                   get_nth_general_transform(second,second_start-second_step));
        }
        else
        {
            second_transform = get_linear_transform_ptr(
                   get_nth_general_transform(second,second_start-second_step));
            second_inverse = get_inverse_linear_transform_ptr(
                   get_nth_general_transform(second,second_start-second_step));
        }

        concat_transforms( get_linear_transform_ptr(transform),
                           first_transform, second_transform );

        concat_transforms( get_inverse_linear_transform_ptr(transform),
                           second_inverse, first_inverse );

        ++trans;
    }

    for( i = second_start;  i != second_end + second_step;  i += second_step )
    {
        copy_and_invert_transform( get_nth_general_transform( second, i ),
                                   second_inverted_concat,
                                   get_nth_general_transform(result_ptr,trans));
        ++trans;
    }

    if( result == first || result == second )
        *result = *result_ptr;
}


// input_transform_file
// Based on minc-1.5.1/volume_io/MNI_formats/gen_xf_io.c
//
static const char* const TRANSFORM_FILE_HEADER   = "MNI Transform File";
static const char* const TYPE_STRING             = "Transform_Type";
static const char* const LINEAR_TYPE             = "Linear";
static const char* const LINEAR_TRANSFORM_STRING = "Linear_Transform";
static const char* const INVERT_FLAG_STRING      = "Invert_Flag";
static const char* const TRUE_STRING 		 = "True";
static const char* const FALSE_STRING 		 = "False";

static VIO_Status input_one_transform(
    FILE                *file,
    const char*          filename,
    General_transform   *transform )
{
    VIO_Status        status;
    int               i, j;
    char* 	      type_name;
    char* 	      str;
    Transform         linear_transform;
    Transform_types   type;
    BOOLEAN           inverse_flag;
    General_transform inverse;
#ifdef BEVIN_ALL_TYPES_SUPPORTED
    int 	      n_points, n_dimensions;
    double            **points, **displacements;
    double	      *points_1d;
    minc_input_options options;
    char* 	      volume_filename;
    char* 	      directory;
    char* 	      tmp_filename;
    Volume            volume;
#endif

    inverse_flag = false;

    /* --- read the type of transform */

    status = mni_input_keyword_and_equal_sign( file, TYPE_STRING, true );

    if( status != OK )
        return( status );

    if( mni_input_string( file, &type_name, (char) ';', (char) 0 ) != OK )
    {
        fprintf(stderr, "input_transform(): missing transform type.\n");
        return( ERROR );
    }
    if( mni_skip_expected_character( file, (char) ';' ) != OK )
        return( ERROR );

    if( !strcmp( type_name, LINEAR_TYPE ) )
        type = LINEAR;
#ifdef BEVIN_ALL_TYPES_SUPPORTED
    else if( equal_strings( type_name, THIN_PLATE_SPLINE_STRING ) )
        type = THIN_PLATE_SPLINE;
    else if( equal_strings( type_name, GRID_TRANSFORM_STRING ) )
        type = GRID_TRANSFORM;
#endif
    else
    {
        fprintf(stderr, "input_transform(): invalid transform type %s.\n", type_name);
        free( type_name );
        return( ERROR );
    }

    free( type_name );

    /* --- read the next string */

    if( mni_input_string( file, &str, (char) '=', (char) 0 ) != OK )
        return( ERROR );

    if( !strcmp( str, INVERT_FLAG_STRING ) )
    {
        free( str );

        if( mni_skip_expected_character( file, (char) '=' ) != OK )
            return( ERROR );
        if( mni_input_string( file, &str, (char) ';', (char) 0 ) != OK )
            return( ERROR );
        if( mni_skip_expected_character( file, (char) ';' ) != OK )
        {
            free( str );
            return( ERROR );
        }

        if( !strcmp( str, TRUE_STRING ) )
            inverse_flag = true;
        else if( !strcmp( str, FALSE_STRING ) )
            inverse_flag = false;
        else
        {
            free( str );
            fprintf(stderr, "Expected %s or %s after %s =\n",
                         TRUE_STRING, FALSE_STRING, INVERT_FLAG_STRING );
            return( ERROR );
        }

        free( str );

        if( mni_input_string( file, &str, (char) '=', (char) 0 ) != OK )
            return( ERROR );
    }

    switch( type )
    {
    case LINEAR:
        if( strcmp( str, LINEAR_TRANSFORM_STRING ) )
        {
            fprintf(stderr, "Expected %s =\n", LINEAR_TRANSFORM_STRING );
            free( str );
            return( ERROR );
        }

        free( str );

        if( mni_skip_expected_character( file, (char) '=' ) != OK )
            return( ERROR );

        make_identity_transform( &linear_transform );

        /* now read the 3 lines of transforms */

        for( i=0; i < 3; i++ )
        {
            for( j=0; j < 4; j++ )
            {
    		double value;
                if( mni_input_double( file, &value ) != OK )
                {
                    fprintf(stderr,
                        "input_transform(): error reading transform elem [%d,%d]\n",
                        i+1, j+1 );
                    return( ERROR );
                }

                Transform_elem(linear_transform,i,j) = value;
            }
        }

        if( mni_skip_expected_character( file, (char) ';' ) != OK )
            return( ERROR );

        create_linear_transform( transform, &linear_transform );

        break;
#if defined(BEVIN_ALL_TYPES_SUPPORTED)
    case THIN_PLATE_SPLINE:

        /* --- read Number_Dimensions = 3; */

        if( !equal_strings( str, N_DIMENSIONS_STRING ) )
        {
            print_error( "Expected %s =\n", N_DIMENSIONS_STRING );
            delete_string( str );
            return( ERROR );
        }

        delete_string( str );

        if( mni_skip_expected_character( file, (char) '=' ) != OK )
            return( ERROR );
        if( mni_input_int( file, &n_dimensions ) != OK )
            return( ERROR );
        if( mni_skip_expected_character( file, (char) ';' ) != OK )
            return( ERROR );

        /* --- read Points = x y z x y z .... ; */

        if( mni_input_keyword_and_equal_sign( file, POINTS_STRING, TRUE ) != OK)
            return( ERROR );
        if( mni_input_reals( file, &n_points, &points_1d ) != OK )
            return( ERROR );

        if( n_points % n_dimensions != 0 )
        {
            print_error(
        "Number of points (%d) must be multiple of number of dimensions (%d)\n",
                  n_points, n_dimensions );
            return( ERROR );
        }

        n_points = n_points / n_dimensions;

        ALLOC2D( points, n_points, n_dimensions );
        for_less( i, 0, n_points )
        {
            for_less( j, 0, n_dimensions )
            {
                points[i][j] = points_1d[IJ(i,j,n_dimensions)];
            }
        }

        FREE( points_1d );

        /* --- allocate and input the displacements */

        ALLOC2D( displacements, n_points + n_dimensions + 1, n_dimensions );

        if( mni_input_keyword_and_equal_sign( file, DISPLACEMENTS_STRING, TRUE )
                                                                       != OK )
            return( ERROR );

        for_less( i, 0, n_points + n_dimensions + 1 )
        {
            for_less( j, 0, n_dimensions )
            {
                if( mni_input_real( file, &value ) != OK )
                {
                    print_error( "Expected more displacements.\n" );
                    return( ERROR );
                }
                displacements[i][j] = value;
            }
        }

        if( mni_skip_expected_character( file, (char) ';' ) != OK )
            return( ERROR );

        create_thin_plate_transform_real( transform, n_dimensions,
                                          n_points, points, displacements );


        FREE2D( points );
        FREE2D( displacements );

        break;

    case GRID_TRANSFORM:

        /*--- read the displacement volume filename */

        if( !equal_strings( str, DISPLACEMENT_VOLUME ) )
        {
            print_error( "Expected %s =\n", DISPLACEMENT_VOLUME );
            delete_string( str );
            return( ERROR );
        }

        delete_string( str );

        if( mni_skip_expected_character( file, (char) '=' ) != OK )
            return( ERROR );

        if( mni_input_string( file, &volume_filename,
                              (char) ';', (char) 0 ) != OK )
            return( ERROR );

        if( mni_skip_expected_character( file, (char) ';' ) != OK )
        {
            delete_string( volume_filename );
            return( ERROR );
        }

        /*--- if the volume filename is relative, add the required directory */

        if( volume_filename[0] != '/' && filename != NULL )
        {
            directory = extract_directory( filename );

            if( string_length(directory) > 0 )
            {
                tmp_filename = concat_strings( directory, "/" );
                concat_to_string( &tmp_filename, volume_filename );
                replace_string( &volume_filename, tmp_filename );
            }

            delete_string( directory );
        }

        /*--- input the displacement volume */

        set_default_minc_input_options( &options );
        set_minc_input_vector_to_scalar_flag( &options, FALSE );

        if( input_volume( volume_filename, 4, NULL, 
                          MI_ORIGINAL_TYPE, FALSE, 0.0, 0.0, 
                          TRUE, &volume, &options ) != OK )
        {
            delete_string( volume_filename );
            return( ERROR );
        }

        delete_string( volume_filename );

        /*--- create the transform */

        create_grid_transform_no_copy( transform, volume );

        break;
#endif
    default:
        fprintf(stderr, "Unhandled type at %s:%d\n", __FILE__, __LINE__);
	exit(1);
    }

    if( inverse_flag )
    {
        create_inverse_general_transform( transform, &inverse );
        delete_general_transform( transform );
        *transform = inverse;
    }

    return( OK );
}

static VIO_Status  input_transform(
    FILE                *file,
    const char*          filename,
    General_transform   *transform )
{
    VIO_Status          status;
    int                 n_transforms;
    char*         	line;
    General_transform   next, concated;

    /* parameter checking */

    if( !file )
    {
        fprintf(stderr, "input_transform(): passed NULL FILE ptr.\n");
        return( ERROR );
    }

    /* okay read the header */

    if( mni_input_string( file, &line, (char) 0, (char) 0 ) != OK )
    {
        free( line );
        fprintf(stderr, "input_transform(): could not read header in file.\n");
        return( ERROR );
    }

    if( strcmp( line, TRANSFORM_FILE_HEADER ) )
    {
        free( line );
        fprintf(stderr, "input_transform(): invalid header in file.\n");
        return( ERROR );
    }

    free( line );

    n_transforms = 0;
    while( (status = input_one_transform( file, filename, &next )) == OK )
    {
        if( n_transforms == 0 )
            *transform = next;
        else
        {
            concat_general_transforms( transform, &next, &concated );
            delete_general_transform( transform );
            delete_general_transform( &next );
            *transform = concated;
        }
        ++n_transforms;
    }

    if( status == ERROR )
    {
        fprintf(stderr, "input_transform: error reading transform.\n" );
        return( ERROR );
    }
    else if( n_transforms == 0 )
    {
        fprintf(stderr, "input_transform: no transform present.\n" );
        return( ERROR );
    }

    return( OK );
}

VIO_Status input_transform_file(
    const char* filename,
    General_transform   *transform ) {

	    
    if (false) printf("%s:%d input_transform_file %s\n", __FILE__, __LINE__,
      filename);
      
    VIO_Status status = ERROR;

    FILE* file;
    status = open_file_with_default_suffix( filename,
                      "xfm",
                      READ_FILE, ASCII_FORMAT, &file );
    if( status == OK )
        status = input_transform( file, filename, transform );
    if( status == OK )
        status = close_file( file );

    if (false) printf("%s:%d input_transform_file %s returning %s\n", __FILE__, __LINE__,
      filename, statusToString(status));
    return( status );
}

// Based on minc-1.5.1/volume_io/Volumes/evaluate.c
// which requires the following...
//
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
static double convert_voxel_to_value(
    Volume  volume,
    double  voxel )
{
    if( volume->real_range_set )
        return( volume->real_value_scale * voxel +
                volume->real_value_translation );
    else
        return( voxel );
}


void set_volume_voxel_value(
    Volume   volume,
    int      x,
    int      y,
    int      z,
    int      t,
    int      v,
    double   voxel )
{
//     SET_VOXEL( volume, v0, v1, v2, v3, v4, voxel );
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
    if( (volume)->is_cached_volume )
        set_cached_volume_voxel( volume, x, y, z, t, v, voxel );
    else
#endif
      	SET_MULTIDIM( (volume)->array, x, y, z, t, v, voxel );
}

double get_volume_voxel_value(
    Volume   volume,
    int      x,
    int      y,
    int      z,
    int      t,
    int      v ) 
{
    double voxel;

#if defined(BEVIN_ALL_VOLUME_MEMBERS)
    if( (volume)->is_cached_volume )
        voxel = (double) get_cached_volume_voxel( volume, x, y, z, t, v );
    else
#endif
        GET_MULTIDIM( voxel, (double), (volume)->array, x, y, z, t, v );

    return( voxel );
}


// Based on minc-1.5.1/volume_io/Volumes/input_volume.c
// which requires the following...
//
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
void  delete_volume_input(
    volume_input_struct   *input_info )
{
    switch( input_info->file_format )
    {
    case  MNC_FORMAT:
        (void) close_minc_input( input_info->minc_file );
        break;

    case  FREE_FORMAT:
        delete_free_format_input( input_info );
        break;
    }
}

// Based on minc-1.5.1/volume_io/Volumes/volumes.c
// which requires the following...
//
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
static double dot_vectors(
    int    n,
    double   v1[],
    double   v2[] )
{
    int   i;
    double  d;

    d = 0.0;
    for( i = 0; i <  n; i++ )
        d += v1[i] * v2[i];

    return( d );
}

static  void   cross_3D_vector(
    double   v1[],
    double   v2[],
    double   cross[] )
{
    cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
    cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
    cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static  void   normalize_vector(
    double   v1[],
    double   v1_normalized[] )
{
    int    d;
    double   mag;

    mag = dot_vectors( VIO_N_DIMENSIONS, v1, v1 );
    if( mag <= 0.0 )
        mag = 1.0;

    mag = sqrt( mag );

    for( d = 0; d < VIO_N_DIMENSIONS; d++ )
        v1_normalized[d] = v1[d] / mag;
}

static  void  assign_voxel_to_world_transform(
    Volume             volume,
    General_transform  *transform )
{
    delete_general_transform( &volume->voxel_to_world_transform );

    volume->voxel_to_world_transform = *transform;
}

static void  reorder_voxel_to_xyz(
    Volume   volume,
    double   voxel[],
    double   xyz[] )
{
    int   c;
    for( c = 0; c < VIO_N_DIMENSIONS; c++ )
    {
        int axis = volume->spatial_axes[c];
        if( axis >= 0 )
            xyz[c] = voxel[axis];
        else
            xyz[c] = 0.0;
    }
}

void compute_world_transform(
    int                 spatial_axes[VIO_N_DIMENSIONS],
    double              separations[],
    double              direction_cosines[][VIO_N_DIMENSIONS],
    double              starts[],
    General_transform   *world_transform )
{
    Transform           transform;
    double              separations_3D[VIO_N_DIMENSIONS];
    double              directions[VIO_N_DIMENSIONS][VIO_N_DIMENSIONS];
    double              starts_3D[VIO_N_DIMENSIONS];
    double              normal[VIO_N_DIMENSIONS];
    int                 dim, c, a1, a2, axis, n_axes;
    int                 axis_list[VIO_N_DIMENSIONS];

    /*--- find how many direction cosines are specified, and set the
          3d separations and starts */

    n_axes = 0;

    for( c = 0; c < VIO_N_DIMENSIONS; c++ )
    {
        axis = spatial_axes[c];
        if( axis >= 0 )
        {
            separations_3D[c] = separations[axis];
            starts_3D[c] = starts[axis];
            directions[c][X] = direction_cosines[axis][X];
            directions[c][Y] = direction_cosines[axis][Y];
            directions[c][Z] = direction_cosines[axis][Z];
            axis_list[n_axes] = c;
            ++n_axes;
        }
        else
        {
            separations_3D[c] = 1.0;
            starts_3D[c] = 0.0;
        }
    }

    if( n_axes == 0 )
    {
        fprintf(stderr, "error compute_world_transform:  no axes.\n" );
        return;
    }

    /*--- convert 1 or 2 axes to 3 axes */

    if( n_axes == 1 )
    {
        a1 = (axis_list[0] + 1) % VIO_N_DIMENSIONS;
        a2 = (axis_list[0] + 2) % VIO_N_DIMENSIONS;

        /*--- create an orthogonal vector */

        directions[a1][X] = directions[axis_list[0]][Y] +
                            directions[axis_list[0]][Z];
        directions[a1][Y] = -directions[axis_list[0]][X] -
                            directions[axis_list[0]][Z];
        directions[a1][Z] = directions[axis_list[0]][Y] -
                            directions[axis_list[0]][X];

        cross_3D_vector( directions[axis_list[0]], directions[a1],
                         directions[a2] );
        normalize_vector( directions[a1], directions[a1] );
        normalize_vector( directions[a2], directions[a2] );
    }
    else if( n_axes == 2 )
    {
        a2 = VIO_N_DIMENSIONS - axis_list[0] - axis_list[1];

        cross_3D_vector( directions[axis_list[0]], directions[axis_list[1]],
               directions[a2] );

        normalize_vector( directions[a2], directions[a2] );
    }

    /*--- check to make sure that 3 axes are not a singular system */

    for( dim = 0; dim < VIO_N_DIMENSIONS ; dim++)
    {
        cross_3D_vector( directions[dim], directions[(dim+1)%VIO_N_DIMENSIONS],
                         normal );
        if( normal[0] == 0.0 && normal[1] == 0.0 && normal[2] == 0.0 )
            break;
    }

    if( dim < VIO_N_DIMENSIONS )
    {
        directions[0][0] = 1.0;
        directions[0][1] = 0.0;
        directions[0][2] = 0.0;
        directions[1][0] = 0.0;
        directions[1][1] = 1.0;
        directions[1][2] = 0.0;
        directions[2][0] = 0.0;
        directions[2][1] = 0.0;
        directions[2][2] = 1.0;
    }

    /*--- make the linear transformation */

    make_identity_transform( &transform );

    for( c = 0; c < VIO_N_DIMENSIONS; c++ )
    {
        for( dim = 0; dim < VIO_N_DIMENSIONS; dim ++ )
        {
            Transform_elem(transform,dim,c) = directions[c][dim] *
                                              separations_3D[c];

            Transform_elem(transform,dim,3) += directions[c][dim] *
                                               starts_3D[c];
        }
    }

    create_linear_transform( world_transform, &transform );
}


void  check_recompute_world_transform(
    Volume  volume )
{
    General_transform        world_transform;

    if( !volume->voxel_to_world_transform_uptodate )
    {
        volume->voxel_to_world_transform_uptodate = true;

        compute_world_transform( volume->spatial_axes,
                                 volume->separations,
                                 volume->direction_cosines,
                                 volume->starts,
                                 &world_transform );

        assign_voxel_to_world_transform( volume, &world_transform );
    }
}

bool convert_dim_name_to_spatial_axis(
    const char* name,
    int     *axis )
{
    *axis = -1;

    if( !strcmp( name, MIxspace ) )
        *axis = X;
    else if( !strcmp( name, MIyspace ) )
        *axis = Y;
    else if( !strcmp( name, MIzspace ) )
        *axis = Z;

    return( *axis >= 0 );
}

static  void  convert_transform_origin_to_starts(
    double      origin[],
    int       	n_volume_dimensions,
    int         spatial_axes[],
    double      dir_cosines[][VIO_N_DIMENSIONS],
    double      starts[] )
{
    int         axis, dim, which[VIO_N_DIMENSIONS], n_axes;

    for( dim = 0; dim < n_volume_dimensions; dim++ )
        starts[dim] = 0.0;

    /*--- get the list of valid axes (which) */

    n_axes = 0;
    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        axis = spatial_axes[dim];
        if( axis >= 0 )
        {
            which[n_axes] = axis;
            ++n_axes;
        }
    }

    /*--- get the starts: computed differently for 1, 2, or 3 axes */

    if( n_axes == 1 )
    {
        double o_dot_c = dot_vectors( VIO_N_DIMENSIONS, origin,                dir_cosines[which[0]] );
        double c_dot_c = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[0]], dir_cosines[which[0]] );

        if( c_dot_c != 0.0 )
            starts[which[0]] = o_dot_c / c_dot_c;
    }
    
    else if( n_axes == 2 )
    {
        // The inverse of matrix [a b  is  1/(ad-bc) [d -b
	//			  c d]                -c a]
	//
	// Set up   [ xdx xdy    [ sx    = [ xdv
	//            ydx ydy ]    sy ]      ydv ]
	// and solve for sx sy by multiplying by the inverse
	//
	//                       [ sx    = 1/(xdx*ydy - xdy*ydx) [ ydy -xdy   [ xdv
	//                         sy ]                            -ydx xdx ]   ydv ]
	//
        double x_dot_v = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[0]], origin );
        double y_dot_v = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[1]], origin );
	
        double x_dot_x = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[0]], dir_cosines[which[0]] );
        double y_dot_y = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[1]], dir_cosines[which[1]] );
	
        double x_dot_y = dot_vectors( VIO_N_DIMENSIONS, dir_cosines[which[0]], dir_cosines[which[1]] );
	double y_dot_x = x_dot_y;
	
        double bottom = x_dot_x * y_dot_y - x_dot_y * y_dot_x;

        if( bottom != 0.0 )
        {
            starts[which[0]] = (x_dot_v * y_dot_y - x_dot_y * y_dot_v) / bottom;
            starts[which[1]] = (y_dot_v * x_dot_x - x_dot_y * x_dot_v) / bottom;
        }
    }
    else if( n_axes == 3 && VIO_N_DIMENSIONS == 3)
    {
        /*--- this is the usual case, solve the equations to find what
              starts give the desired origin */

    	Double4x4 matrix;	// only need 3x3 but have a 4x4

	int i,j;
        for( i = 0; i < VIO_N_DIMENSIONS; i++ )
        for( j = 0; j < VIO_N_DIMENSIONS; j++ )
        {
            matrix[Index4x4(i,j)] = dir_cosines[which[j]][i];		// I don't know why this transposes the matrix
        }
	
	for( i = 0; i < 4; i++ )					// Extend to 4x4
	    matrix[Index4x4(3,i)] = matrix[Index4x4(i,3)] = 0.0;
	matrix[Index4x4(3,3)] = 1.0;

	//                                     t         t
	// Find a solution to 	matrix solution  = origin
	//
    	Double4x4 inverse;	// only need 3x3 but have a 4x4
	if (!invert_4x4_matrix((const Double4x4*)&matrix, &inverse))
	{
	    fprintf(stderr,"%s:%d Could not invert matrix\n", __FILE__, __LINE__);
	    exit(1); 
	}
	
        for (i = 0; i < 3; i++)
	{
	    double solution_i = 0;
	    for (j = 0; j < 3; j++)
	    {
	        solution_i += inverse[Index4x4(i,j)] * origin[j];
	    }
            starts[which[i]] = solution_i;
        }
    }
    else
    {
        fprintf(stderr,
          "Invalid number of axes in convert_transform_origin_to_starts\n");
	exit(1);
    }
}

void  convert_transform_to_starts_and_steps(
    General_transform  *transform,
    int                n_volume_dimensions,
    double             step_signs[],
    int                spatial_axes[],
    double             starts[],
    double             steps[],
    double             dir_cosines[][VIO_N_DIMENSIONS] )
{
    double      sign, mag;
    int         axis, dim;
    double      axes[VIO_N_DIMENSIONS][VIO_N_DIMENSIONS];
    double      origin[VIO_N_DIMENSIONS];
    Transform   *linear_transform;

    if( transform->type != LINEAR )
    {
        fprintf(stderr, "convert_transform_to_starts_and_steps(): non-linear transform found.\n" );
        exit(1);
    }

    linear_transform = get_linear_transform_ptr( transform );

    get_transform_origin_real( linear_transform, origin );
    get_transform_x_axis_real( linear_transform, &axes[X][0] );
    get_transform_y_axis_real( linear_transform, &axes[Y][0] );
    get_transform_z_axis_real( linear_transform, &axes[Z][0] );

    /*--- assign default steps */

    for( dim = 0; dim < n_volume_dimensions; dim++ )
        steps[dim] = 1.0;

    /*--- assign the steps and dir_cosines for the spatial axes */

    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        axis = spatial_axes[dim];
        if( axis >= 0 )
        {
            mag = dot_vectors( VIO_N_DIMENSIONS, axes[dim], axes[dim] );

            if( mag <= 0.0 )
                mag = 1.0;
            mag = sqrt( mag );

            if( step_signs == NULL )
            {
                if( axes[dim][dim] < 0.0 )
                    sign = -1.0;
                else
                    sign = 1.0;
            }
            else  /*--- make the sign of steps match the step_signs passed in */
            {
                if( step_signs[axis] < 0.0 )
                    sign = -1.0;
                else
                    sign = 1.0;
            }

            steps[axis] = sign * mag;
            dir_cosines[axis][X] = axes[dim][X] / steps[axis];
            dir_cosines[axis][Y] = axes[dim][Y] / steps[axis];
            dir_cosines[axis][Z] = axes[dim][Z] / steps[axis];
        }
    }

    /*--- finally, get the starts */

    convert_transform_origin_to_starts( origin, n_volume_dimensions,
                                        spatial_axes, dir_cosines, starts );
}

static const char* default_dimension_names[VIO_MAX_DIMENSIONS][VIO_MAX_DIMENSIONS] =
{
    { MIxspace },
    { MIyspace, MIxspace },
    { MIzspace, MIyspace, MIxspace },
    { "", MIzspace, MIyspace, MIxspace },
    { "", "", MIzspace, MIyspace, MIxspace }
};

void  set_volume_space_type(
    Volume   volume,
    const char*  name )
{
    free( (void*)volume->coordinate_system_name );
    volume->coordinate_system_name = strdup( name );
}

int set_volume_irregular_starts(Volume volume, int idim, int count, double *starts)
{
#if !defined(BEVIN_ALL_VOLUME_MEMBERS)
    fprintf(stderr, "%s:%d set_volume_irregular_starts NYI\n", __FILE__, __LINE__);
    exit(1);
    return 0;
#else
    int i;

    if (idim >= volume->array.n_dimensions) {
        return (0);
    }

    if (volume->irregular_starts[idim] != NULL) {
        free(volume->irregular_starts[idim]);
    }

    if (starts == NULL) {
        return (0);
    }

    if (count > volume->array.sizes[idim]) {
        count = volume->array.sizes[idim];
    }

    volume->irregular_starts[idim] = malloc(count * sizeof (double));
    if (volume->irregular_starts[idim] == NULL) {
        return (0);
    }

    for (i = 0; i < count; i++) {
        volume->irregular_starts[idim][i] = starts[i];
    }

    return (count);
#endif
}

int set_volume_irregular_widths(Volume volume, int idim, int count, double *widths)
{
#if !defined(BEVIN_ALL_VOLUME_MEMBERS)
    fprintf(stderr, "%s:%d set_volume_irregular_widths NYI\n", __FILE__, __LINE__);
    exit(1);
    return 0;
#else
    int i;

    if (idim >= volume->array.n_dimensions) {
        return (0);
    }

    if (volume->irregular_widths[idim] != NULL) {
        free(volume->irregular_widths[idim]);
    }

    if (widths == NULL) {
        return (0);
    }

    if (count > volume->array.sizes[idim]) {
        count = volume->array.sizes[idim];
    }

    volume->irregular_widths[idim] = malloc(count * sizeof (double));
    if (volume->irregular_widths[idim] == NULL) {
        return (0);
    }

    for (i = 0; i < count; i++) {
        volume->irregular_widths[idim][i] = widths[i];
    }

    return (count);
#endif
}


static double  get_volume_real_min(
    Volume     volume )
{
    double   real_min = volume->voxel_min;

    if( volume->real_range_set )
        real_min = convert_voxel_to_value( volume, real_min );

    return( real_min );
}

static double  get_volume_real_max(
    Volume     volume )
{
    double   real_max = volume->voxel_max;

    if( volume->real_range_set )
        real_max = convert_voxel_to_value( volume, real_max );

    return( real_max );
}

void  get_volume_voxel_range(
    Volume volume,
    double *voxel_min,
    double *voxel_max )
{
    *voxel_min = volume->voxel_min;
    *voxel_max = volume->voxel_max;
}

void  set_volume_real_range(
    Volume volume,
    double real_min,
    double real_max )
{
    double    voxel_min, voxel_max;

    if( get_volume_data_type(volume) == FLOAT ||
        get_volume_data_type(volume) == DOUBLE )
    {
        set_volume_voxel_range( volume, real_min, real_max );
        volume->real_value_scale = 1.0;
        volume->real_value_translation = 0.0;
    }
    else
    {
        get_volume_voxel_range( volume, &voxel_min, &voxel_max );

        if( voxel_min < voxel_max )
        {
            volume->real_value_scale = (real_max - real_min) /
                                       (voxel_max - voxel_min);
            volume->real_value_translation = real_min -
                                       voxel_min * volume->real_value_scale;
        }
        else
        {
	    // FIXME: is scale = 0 correct??
            volume->real_value_scale = 0.0;
            volume->real_value_translation = real_min;
        }

        volume->real_range_set = true;
    }

#if 0	// cache_volume_range_has_changed is basically NYI in minc
    if( volume->is_cached_volume )
        cache_volume_range_has_changed( volume );
#endif
}

void  get_volume_real_range(
    Volume     volume,
    double       *min_value,
    double       *max_value )
{
    *min_value = get_volume_real_min( volume );
    *max_value = get_volume_real_max( volume );
}


void  set_volume_voxel_range(
    Volume   volume,
    double   voxel_min,
    double  voxel_max )
{
    double real_min = 0.0, real_max = 0.0;

    if( voxel_min >= voxel_max )
    {
        switch( get_volume_data_type( volume ) )
        {
        case UNSIGNED_BYTE:
            voxel_min = 0.0;
            voxel_max = (double) UCHAR_MAX;     break;
        case SIGNED_BYTE:
            voxel_min = (double) SCHAR_MIN;
            voxel_max = (double) SCHAR_MAX;     break;
        case UNSIGNED_SHORT:
            voxel_min = 0.0;
            voxel_max = (double) USHRT_MAX;     break;
        case SIGNED_SHORT:
            voxel_min = (double) SHRT_MIN;
            voxel_max = (double) SHRT_MAX;      break;
        case UNSIGNED_INT:
            voxel_min = 0.0;
            voxel_max = (double) UINT_MAX;     break;
        case SIGNED_INT:
            voxel_min = (double) INT_MIN;
            voxel_max = (double) INT_MAX;      break;
        case FLOAT:
            voxel_min = (double) -FLT_MAX;
            voxel_max = (double) FLT_MAX;       break;
        case DOUBLE:
            voxel_min = (double) -DBL_MAX;
            voxel_max = (double) DBL_MAX;       break;
	default:
	    fprintf(stderr, "%s:%d no default\n",__FILE__,__LINE__);
	    exit(1);
        }
    }

    if( volume->real_range_set )
        get_volume_real_range( volume, &real_min, &real_max );

    volume->voxel_min = voxel_min;
    volume->voxel_max = voxel_max;

    if( volume->real_range_set )
        set_volume_real_range( volume, real_min, real_max );
#if 0	// cache_volume_range_has_changed is basically NYI in minc
    else
        cache_volume_range_has_changed( volume );
#endif
}

void  set_volume_type(
    Volume       volume,
    nc_type      nc_data_type,
    bool         signed_flag,
    double       voxel_min,
    double       voxel_max )
{
    VIO_Data_types      data_type;

    if( nc_data_type != (nc_type)0 /* MI_ORIGINAL_TYPE */ )
    {
        switch( nc_data_type )
        {
        case  NC_BYTE:
            if( signed_flag )
                data_type = SIGNED_BYTE;
            else
                data_type = UNSIGNED_BYTE;
            break;

        case  NC_SHORT:
            if( signed_flag )
                data_type = SIGNED_SHORT;
            else
                data_type = UNSIGNED_SHORT;
            break;

        case  NC_LONG:
            if( signed_flag )
                data_type = SIGNED_INT;
            else
                data_type = UNSIGNED_INT;
            break;

        case  NC_FLOAT:
            data_type = FLOAT;
            break;

        case  NC_DOUBLE:
            data_type = DOUBLE;
            break;
	    
	default:
	    fprintf(stderr, "%s:%d bad nc_data_type\n", __FILE__, __LINE__);
            exit(1);
        }

        set_multidim_data_type( &volume->array, data_type );

        volume->signed_flag = signed_flag;

        set_volume_voxel_range( volume, voxel_min, voxel_max );
    }

    volume->nc_data_type = nc_data_type;
}

VIO_Data_types  get_volume_data_type(
    Volume       volume )
{
    return( get_multidim_data_type( &volume->array ) );
}

int get_volume_n_dimensions(
    Volume volume )
{
    return( get_multidim_n_dimensions( &volume->array ) );
}

void get_volume_sizes(
    Volume 	volume,
    int      	sizes[] )
{
    get_multidim_sizes( &volume->array, sizes );
}


nc_type get_volume_nc_data_type(
    Volume      volume,
    bool*	signed_flag ) 
{
    if( signed_flag != (BOOLEAN *) NULL )
        *signed_flag = volume->signed_flag;
    return( volume->nc_data_type );
}

void  set_rgb_volume_flag(
    Volume   volume,
    bool     flag )
{
    if( !flag || get_volume_data_type(volume) == UNSIGNED_INT )
        volume->is_rgba_data = flag;
}

Volume create_volume(
    int          n_dimensions,
    /*const*/ char*  dimension_names[],		// need compat with minc
    nc_type      nc_data_type,
    bool         signed_flag,
    double	 voxel_min,
    double       voxel_max ) 
{
    int             i;
    VIO_Status      status;
    volume_struct   *volume;
    Transform       identity;

    status = OK;

    if( n_dimensions < 1 || n_dimensions > VIO_MAX_DIMENSIONS )
    {
        fprintf(stderr,
            "create_volume(): n_dimensions (%d) not in range 1 to %d.\n",
               n_dimensions, VIO_MAX_DIMENSIONS );
        status = ERROR;
    }

    if( status == ERROR )
    {
        return( (Volume) NULL );
    }

    volume = (volume_struct*)malloc(sizeof(volume_struct));

    volume->is_rgba_data     = false;
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
    volume->is_cached_volume = false;
#endif

    volume->real_range_set         = false;
    volume->real_value_scale       = 1.0;
    volume->real_value_translation = 0.0;

    for( i = 0; i < VIO_N_DIMENSIONS; i++ )
        volume->spatial_axes[i] = -1;

    int sizes[VIO_MAX_DIMENSIONS];
    for( i = 0; i < n_dimensions; i++ )
    {
        volume->starts[i] = 0.0;
        volume->separations[i] = 1.0;
        volume->direction_cosines[i][X] = 0.0;
        volume->direction_cosines[i][Y] = 0.0;
        volume->direction_cosines[i][Z] = 0.0;
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
        volume->irregular_starts[i]     = NULL;
        volume->irregular_widths[i]     = NULL;
#endif

        sizes[i] = 0;

        const char* name;
        if( dimension_names != NULL )
            name = dimension_names[i];
        else
            name = default_dimension_names[n_dimensions-1][i];

        volume->dimension_names[i] = strdup( name );

        int axis;
        if( convert_dim_name_to_spatial_axis( name, &axis ) )
        {
            volume->spatial_axes        [axis] = i;
            volume->direction_cosines[i][axis] = 1.0;
        }

    }

    create_empty_multidim_array( &volume->array, n_dimensions, NO_DATA_TYPE );

    set_volume_type( volume, nc_data_type, signed_flag, voxel_min, voxel_max );
    set_volume_sizes( volume, sizes );

    make_identity_transform( &identity );
    create_linear_transform( &volume->voxel_to_world_transform, &identity );
    volume->voxel_to_world_transform_uptodate = true;

    volume->coordinate_system_name = strdup( MI_UNKNOWN_SPACE );

    return( volume );
}

bool volume_is_alloced(
    Volume   volume )
{
    return( 
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
            volume->is_cached_volume &&
            volume_cache_is_alloced( &volume->cache ) ||
            !volume->is_cached_volume &&
#endif
            multidim_array_is_alloced( &volume->array ) );
}

void  free_volume_data(
    Volume   volume )
{
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
    if( volume->is_cached_volume )
        delete_volume_cache( &volume->cache, volume );
    else 
#endif
    if( volume_is_alloced( volume ) )
        delete_multidim_array( &volume->array );
}

void delete_volume(
    Volume volume ) {

    if( volume == (Volume) NULL )
    {
        fprintf(stderr, "delete_volume():  cannot delete a null volume.\n" );
	exit(1);
        return;
    }

    free_volume_data( volume );

    delete_general_transform( &volume->voxel_to_world_transform );

    int   d;
    for( d = 0; d < get_volume_n_dimensions(volume); d++ )
        free( (void*)volume->dimension_names[d] );

    free( (void*)volume->coordinate_system_name );

    free( volume );
}

void  set_volume_sizes(
    Volume       volume,
    int          sizes[] )
{
    set_multidim_sizes( &volume->array, sizes );
}

void alloc_volume_data(
    Volume      volume ) 
{
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
    unsigned long   data_size;

    data_size = (unsigned long) get_volume_total_n_voxels( volume ) *
                (unsigned long) get_type_size( get_volume_data_type( volume ) );

    if( get_n_bytes_cache_threshold() >= 0 &&
        data_size > (unsigned long) get_n_bytes_cache_threshold() )
    {
        volume->is_cached_volume = true;
        initialize_volume_cache( &volume->cache, volume );
    }
    else
#endif
    {
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
        volume->is_cached_volume = false;
#endif
        alloc_multidim_array( &volume->array );
    }
}


void set_volume_separations(
    Volume      volume,
    double      separations[] ) 
{
    int i;
    for( i = 0; i < get_volume_n_dimensions( volume ); i++ )
        volume->separations[i] = separations[i];

    volume->voxel_to_world_transform_uptodate = false;
}

void  set_volume_direction_unit_cosine(
    Volume   volume,
    int      axis,
    double   dir[] )
{
    if( axis < 0 || axis >= get_volume_n_dimensions(volume) )
    {
        fprintf(stderr,
         "set_volume_direction_cosine:  cannot set dir cosine for axis %d\n",
          axis );
	exit(1);
        return;
    }

    /*--- check if this is a spatial axis */

    int    dim;

    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        if( volume->spatial_axes[dim] == axis )
            break;
    }

    if( dim == VIO_N_DIMENSIONS )   /* this is not a spatial axis, ignore the dir */
        return;

    volume->direction_cosines[axis][X] = dir[X];
    volume->direction_cosines[axis][Y] = dir[Y];
    volume->direction_cosines[axis][Z] = dir[Z];

    volume->voxel_to_world_transform_uptodate = false;
}

void set_volume_direction_cosine(
    Volume   	volume,
    int      	axis,
    double	dir[] )
{
    double len, unit_vector[VIO_N_DIMENSIONS];

    len = dir[X] * dir[X] + dir[Y] * dir[Y] + dir[Z] * dir[Z];

    if( len == 0.0 )
    {
        fprintf(stderr, "Warning: zero length direction cosine in set_volume_direction_cosine()\n" );
	exit(1);
        return;
    }

    if( len <= 0.0 )
        len = 1.0;

    len = sqrt( len );

    unit_vector[X] = dir[X] / len;
    unit_vector[Y] = dir[Y] / len;
    unit_vector[Z] = dir[Z] / len;

    set_volume_direction_unit_cosine( volume, axis, unit_vector );
}


void set_volume_starts(
    Volume  volume,
    double  starts[] )
{
    int  c;
    for( c =  0; c < get_volume_n_dimensions( volume ); c++ )
        volume->starts[c] = starts[c];
    volume->voxel_to_world_transform_uptodate = false;
}

void set_volume_translation(
    Volume  	volume,
    double    	voxel[],
    double    	world_space_voxel_maps_to[] )
{
    int         dim, dim2, axis, n_axes, a1, a2;
    double      world_space_origin[VIO_N_DIMENSIONS], len;
    double      starts[VIO_MAX_DIMENSIONS], starts_3d[VIO_N_DIMENSIONS];
    Transform   transform, inverse;

    /*--- find the world position where ( 0, 0, 0 ) maps to by taking
          the world position - voxel[x_axis] * Xaxis - voxel[y_axis] * Yaxis
          ..., and fill in the transform defined by Xaxis, Yaxis, Zaxis */

    make_identity_transform( &transform );

    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        world_space_origin[dim] = world_space_voxel_maps_to[dim];

        for( dim2 = 0; dim2 < VIO_N_DIMENSIONS; dim2++ )
        {
            axis = volume->spatial_axes[dim2];
            if( axis >= 0 )
            {
                world_space_origin[dim] -= volume->separations[axis] *
                           volume->direction_cosines[axis][dim] * voxel[axis];

                Transform_elem( transform, dim, dim2 ) =
                                           volume->direction_cosines[axis][dim];
            }
        }
    }

    n_axes = 0;

    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        axis = volume->spatial_axes[dim];
        if( axis >= 0 )
            ++n_axes;
    }

    /*--- if only one spatial axis, make a second orthogonal vector */

    if( n_axes == 1 )
    {
        /*--- set dim to the spatial axis */

        if( volume->spatial_axes[0] >= 0 )
            dim = 0;
        else if( volume->spatial_axes[1] >= 0 )
            dim = 1;
        else if( volume->spatial_axes[2] >= 0 )
            dim = 2;

        /*--- set a1 to the lowest occuring non-spatial axis, and create
              a unit vector normal to that of the spatial axis */

        if( dim == 0 )
            a1 = 1;
        else
            a1 = 0;

        Transform_elem( transform, 0, a1 ) = Transform_elem(transform,1,dim) +
                                             Transform_elem(transform,2,dim);
        Transform_elem( transform, 1, a1 ) = -Transform_elem(transform,0,dim) -
                                              Transform_elem(transform,2,dim);
        Transform_elem( transform, 2, a1 ) = Transform_elem(transform,1,dim) -
                                             Transform_elem(transform,0,dim);

        len = Transform_elem(transform,0,a1)*Transform_elem(transform,0,a1) +
              Transform_elem(transform,1,a1)*Transform_elem(transform,1,a1) +
              Transform_elem(transform,2,a1)*Transform_elem(transform,2,a1);

        if( len == 0.0 )
            len = 1.0;
        else
            len = sqrt( len );

        Transform_elem(transform,0,a1) /= len;
        Transform_elem(transform,1,a1) /= len;
        Transform_elem(transform,2,a1) /= len;
    }

    /*--- if only two spatial axis, make a third orthogonal vector */

    if( n_axes == 1 || n_axes == 2 )
    {
        /*--- set dim to the one axis that does not have a vector associated
              with it yet, and make one that is the unit cross product of 
              the other two */

        if( volume->spatial_axes[2] < 0 )
            dim = 2;
        else if( volume->spatial_axes[1] < 0 )
            dim = 1;
        else if( volume->spatial_axes[0] < 0 )
            dim = 0;

        a1 = (dim + 1) % VIO_N_DIMENSIONS;
        a2 = (dim + 2) % VIO_N_DIMENSIONS;

        /*--- take cross product */

        Transform_elem( transform, 0, dim ) = Transform_elem(transform,1,a1) *
                                              Transform_elem(transform,2,a2) -
                                              Transform_elem(transform,1,a2) *
                                              Transform_elem(transform,2,a1);
        Transform_elem( transform, 1, dim ) = Transform_elem(transform,2,a1) *
                                              Transform_elem(transform,0,a2) -
                                              Transform_elem(transform,2,a2) *
                                              Transform_elem(transform,0,a1);
        Transform_elem( transform, 2, dim ) = Transform_elem(transform,0,a1) *
                                              Transform_elem(transform,1,a2) -
                                              Transform_elem(transform,0,a2) *
                                              Transform_elem(transform,1,a1);

        /*--- normalize vector */

        len = Transform_elem(transform,0,dim)*Transform_elem(transform,0,dim) +
              Transform_elem(transform,1,dim)*Transform_elem(transform,1,dim) +
              Transform_elem(transform,2,dim)*Transform_elem(transform,2,dim);

        if( len == 0.0 )
            len = 1.0;
        else
            len = sqrt( len );

        Transform_elem(transform,0,dim) /= len;
        Transform_elem(transform,1,dim) /= len;
        Transform_elem(transform,2,dim) /= len;
    }

    /*--- find the voxel that maps to the world space origin, when there is
          no translation, and this is the starts */

    compute_transform_inverse( &transform, &inverse );

    transform_point( &inverse, world_space_origin[X],
                               world_space_origin[Y],
                               world_space_origin[Z],
                               &starts_3d[X], &starts_3d[Y], &starts_3d[Z] );

    /*--- map the X Y Z starts into the arbitrary axis ordering of the volume */

    for( dim = 0; dim < get_volume_n_dimensions(volume); dim++ )
        starts[dim] = 0.0;

    for( dim = 0; dim < VIO_N_DIMENSIONS; dim++ )
    {
        axis = volume->spatial_axes[dim];
        if( axis >= 0 )
            starts[axis] = starts_3d[dim];
    }

    set_volume_starts( volume, starts );
}


General_transform  *get_voxel_to_world_transform(
    Volume   volume )
{
    check_recompute_world_transform( volume );

    return( &volume->voxel_to_world_transform );
}


static const char* convert_spatial_axis_to_dim_name(
    int   axis )
{
    switch( axis )
    {
    case X:  return( MIxspace );
    case Y:  return( MIyspace );
    case Z:  return( MIzspace );
    default: 
    	fprintf(stderr, "convert_spatial_axis_to_dim_name, bad axis\n" );
        exit(1);
    }
    return( NULL );
}

char** get_volume_dimension_names(
    Volume   volume )
{
    char*   *names = (char**) calloc(sizeof(char*), get_volume_n_dimensions(volume) );

    int      i;
    for( i = 0; i < get_volume_n_dimensions(volume); i++ )
        names[i] = 
	    strdup( 
	    	(volume->spatial_axes[i] >= 0)
		? convert_spatial_axis_to_dim_name(i)
		: volume->dimension_names[i] );

    return( names );
}


void  delete_dimension_names(
    Volume   volume,
    char*    dimension_names[] )
{
    int   i;
    for( i = 0; i < get_volume_n_dimensions(volume); i++ )
        free( dimension_names[i] );

    free( dimension_names );
}



void get_volume_separations(
    Volume   	volume,
    double 	separations[] )
{
    int   i;

    for( i = 0; i < get_volume_n_dimensions( volume ); i++ )
        separations[i] = volume->separations[i];
}


void  convert_voxel_to_world(
    Volume   volume,
    double   voxel[],
    double  *x_world,
    double  *y_world,
    double  *z_world )
{
    double xyz[VIO_N_DIMENSIONS];

    check_recompute_world_transform( volume );

    reorder_voxel_to_xyz( volume, voxel, xyz );

    /* apply linear transform */

    general_transform_point( &volume->voxel_to_world_transform,
                             xyz[X], xyz[Y], xyz[Z],
                             x_world, y_world, z_world );
}

VIO_Status start_volume_input(
    const char*          filename,
    int                  n_dimensions,
    /*const*/ char*      dim_names[],			// need compat with minc
    nc_type              volume_nc_data_type,
    bool                 volume_signed_flag,
    double               volume_voxel_min,
    double               volume_voxel_max,
    bool                 create_volume_flag,
    Volume              *volume,
    minc_input_options  *options,
    volume_input_struct *input_info ) 
{
    VIO_Status status = OK;

    if( create_volume_flag || *volume == (Volume) NULL )
    {
        if( n_dimensions < 1 || n_dimensions > VIO_MAX_DIMENSIONS )
            n_dimensions = get_minc_file_n_dimensions( filename );

        if( n_dimensions < 1 )
            return( ERROR );

        if( dim_names == NULL )
           dim_names = (char**)default_dimension_names[ n_dimensions-1 ];

        *volume = create_volume( n_dimensions, dim_names, volume_nc_data_type,
                                 volume_signed_flag,
                                 volume_voxel_min, volume_voxel_max );
    }
    else if( n_dimensions != get_volume_n_dimensions( *volume ) &&
             volume_is_alloced( *volume ) )
        free_volume_data( *volume );

    const char* expanded_filename = expand_filename( filename );

    if( !filename_extension_matches( expanded_filename, "fre" ) )
        input_info->file_format = MNC_FORMAT;
    else
        input_info->file_format = FREE_FORMAT;

    switch( input_info->file_format )
    {
    case  MNC_FORMAT:
        if( !file_exists( expanded_filename ) )
        {
	    const char* old_filename = expanded_filename;
	    expanded_filename = file_exists_as_compressed( old_filename );
	    free( (void*)old_filename );
        }

        input_info->minc_file = 
	    initialize_minc_input( expanded_filename, *volume, options );
	    
        if( input_info->minc_file == (Minc_file) NULL )
            status = ERROR;
        else
        {
	    int d;
            for( d = 0; d < VIO_MAX_DIMENSIONS; d++ )
                input_info->axis_index_from_file[d] = d;
        }

        break;

    case  FREE_FORMAT:
        status = initialize_free_format_input( expanded_filename,
                                               *volume, input_info );
        break;
    }

    if( status != OK && create_volume_flag )
        delete_volume( *volume );

    free( (void*)expanded_filename );

    return( status );
}


bool input_more_of_volume(
    Volume                volume,
    volume_input_struct  *input_info,
    double               *fraction_done )
{
    bool more_to_do = false;

    switch( input_info->file_format )
    {
    case  MNC_FORMAT:
        more_to_do = input_more_minc_file( input_info->minc_file,
                                           fraction_done );
        break;

    case  FREE_FORMAT:
        more_to_do = input_more_free_format_file( volume, input_info,
                                                  fraction_done );
        break;
    }

    return( more_to_do );
}
