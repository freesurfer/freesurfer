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

    result_ptr->n_transforms = ABS( first_end - first_start ) + 1 +
                               ABS( second_end - second_start ) + 1;

    crunching_linear = FALSE;
    if( get_nth_general_transform( first, first_end )->type == LINEAR &&
        get_nth_general_transform( second, second_start )->type == LINEAR )
    {
        --result_ptr->n_transforms;
        crunching_linear = TRUE;
        first_end -= first_step;
        second_start += second_step;
    }

    if( result_ptr->n_transforms == 1 )
        result_ptr->type = LINEAR;
    else
    {
        result_ptr->type = CONCATENATED_TRANSFORM;
        ALLOC( result_ptr->transforms, result_ptr->n_transforms );
    }

    result_ptr->inverse_flag = FALSE;

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

static VIO_Status  input_transform(
    FILE                *file,
    const char*          filename,
    General_transform   *transform )
{
    VIO_Status          status;
    int                 n_transforms;
    const char*         line;
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
        delete_string( line );
        print_error( "input_transform(): could not read header in file.\n");
        return( ERROR );
    }

    if( !equal_strings( line, TRANSFORM_FILE_HEADER ) )
    {
        delete_string( line );
        print_error( "input_transform(): invalid header in file.\n");
        return( ERROR );
    }

    delete_string( line );

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
        print_error( "input_transform: error reading transform.\n" );
        return( ERROR );
    }
    else if( n_transforms == 0 )
    {
        print_error( "input_transform: no transform present.\n" );
        return( ERROR );
    }

    return( OK );
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

