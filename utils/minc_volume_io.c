/**
 * @file  minc_volume_io.c
 * @brief substitutes for the needed functionality previously obtained from minc
 */
/*
 * Original Author: Bevin Brett
 * CVS Revision Info:
 *    $Author: Brett $
 *    $Date: 2017/11/01 12:46:00 $
 *    $Revision: 1.00 $
 *
 * Copyright © 2017 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "minc_volume_io.h"
#include "files.h"


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

#if 0
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
            in_comment = ( *ch == COMMENT_CHAR1 || *ch == COMMENT_CHAR2 );
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
    VIO_Status  status = mni_input_string( file, &str, (char) '=', (char) 0 );

    if( status == END_OF_FILE )
        return( status );

    if( status != OK || strcmp( str, keyword ) ||
        mni_skip_expected_character( file, (char) '=' ) != OK )
    {
        if( print_error_message )
            fprintf(stderr, "Expected \"%s =\"\n", keyword );
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


static void nyi(const char* file, int line, const char* name) {
    fprintf(stderr, "%s:%d NYI %s\n", file, line, name);
    exit(1);
}

#define NYI(NAME, VALUE) { nyi(__FILE__, __LINE__, NAME); return VALUE; }

Transform* get_linear_transform_ptr(
    General_transform   *transform ) NYI("get_linear_transform_ptr",NULL)

Transform* get_inverse_linear_transform_ptr(
    General_transform   *transform ) NYI("get_inverse_linear_transform_ptr",NULL)

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

void copy_general_transform(
    General_transform   *transform,
    General_transform   *copy ) NYI("copy_general_transform",)

void delete_general_transform(
    General_transform   *transform ) NYI("delete_general_transform",)


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

    status = mni_input_keyword_and_equal_sign( file, TYPE_STRING, false );

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

	    
    printf("%s:%d input_transform_file %s\n", __FILE__, __LINE__,
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

    printf("%s:%d input_transform_file %s returning %s\n", __FILE__, __LINE__,
      filename, statusToString(status));
    return( status );
}


void transform_point(
    Transform  *transform,
    double 	x,
    double	y,
    double 	z,
    double	*x_trans,
    double	*y_trans,
    double	*z_trans ) NYI("transform_point",)

Volume create_volume(
    int          n_dimensions,
    char*        dimension_names[],
    nc_type      nc_data_type,
    bool         signed_flag,
    double	 voxel_min,
    double       voxel_max ) NYI("create_volume",NULL)

void delete_volume(
    Volume volume ) NYI("",)

void delete_volume_input(
    volume_input_struct   *input_info ) NYI("delete_volume",)

void set_volume_sizes(
    Volume   	volume,
    int         sizes[] ) NYI("set_volume_sizes",)

void alloc_volume_data(
    Volume      volume ) NYI("alloc_volume_data",)

void set_volume_separations(
    Volume      volume,
    double      separations[] ) NYI("set_volume_separations",)

void set_volume_direction_cosine(
    Volume   	volume,
    int      	axis,
    double	dir[] ) NYI("set_volume_direction_cosine",)

void set_volume_translation(
    Volume  	volume,
    double    	voxel[],
    double    	world_space_voxel_maps_to[] ) NYI("set_volume_translation",)

void set_volume_voxel_value(
    Volume      volume,
    int         v0,
    int         v1,
    int         v2,
    int         v3,
    int         v4,
    double      voxel ) NYI("set_volume_voxel_value",)

VIO_Status  output_volume(
    const char*		  filename,
    nc_type		  file_nc_data_type,
    bool              	  file_signed_flag,
    double                file_voxel_min,
    double                file_voxel_max,
    Volume                volume,
    const char*	  	  history,
    minc_output_options  *options ) NYI("output_volume",OK)

VIO_Status start_volume_input(
    char*                filename,
    int                  n_dimensions,
    char*                dim_names[],
    nc_type              volume_nc_data_type,
    bool                 volume_signed_flag,
    double               volume_voxel_min,
    double               volume_voxel_max,
    bool                 create_volume_flag,
    Volume              *volume,
    minc_input_options  *options,
    volume_input_struct *input_info ) NYI("start_volume_input",OK)

bool input_more_of_volume(
    Volume                volume,
    volume_input_struct  *input_info,
    double               *fraction_done ) NYI("input_more_of_volume",0)

int get_volume_n_dimensions(
    Volume volume ) NYI("get_volume_n_dimensions",0)

void get_volume_sizes(
    Volume 	volume,
    int      	sizes[] ) NYI("get_volume_sizes",)

nc_type get_volume_nc_data_type(
    Volume      volume,
    bool*	signed_flag ) NYI("get_volume_nc_data_type",NC_UNSPECIFIED)

General_transform* get_voxel_to_world_transform(
    Volume   volume ) NYI("get_voxel_to_world_transform",NULL)

double get_volume_voxel_value(
    Volume   volume,
    int      v0,
    int      v1,
    int      v2,
    int      v3,
    int      v4 ) NYI("get_volume_voxel_value",0)
    
void get_volume_separations(
    Volume   	volume,
    double 	separations[] ) NYI("get_volume_separations",)

void convert_voxel_to_world(
    Volume   	volume,
    double     	voxel[],
    double     *x_world,
    double     *y_world,
    double     *z_world ) NYI("convert_voxel_to_world",)

