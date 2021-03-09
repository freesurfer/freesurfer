/**
 * @brief Wrapper for MNI's volume_io.h, to decouple from MNI lib
 *
 */
/*
 * Original Author: Nick Schmansky
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

#ifndef MINC_VOLUME_IO_H
#define MINC_VOLUME_IO_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "netcdf.h"
#include "minc_multidim_arrays.h"

typedef bool BOOLEAN;


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/basic.h
//
// As such, it needs the following Copyright notice
/*
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
*/
typedef enum { 
               OK,
               ERROR,
               INTERNAL_ERROR,
               END_OF_FILE,
               QUIT,
	       VIO_Status__end
             } VIO_Status;

char* statusToString(VIO_Status status);


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/multidim.h
//
// As such, it needs the following Copyright notice
/*
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
*/

#define  VIO_MAX_DIMENSIONS     5


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/geom_structs.h
//
// As such, it needs the following Copyright notice
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
@VERSION    : $Header: /private-cvsroot/minc/volume_io/Include/volume_io/geom_structs.h,v 1.20.2.3 2006/11/30 09:15:13 rotor Exp $
---------------------------------------------------------------------------- */

typedef  unsigned  int     VIO_Colour;


#define VIO_N_DIMENSIONS 3

typedef struct
{
    double m[4][4];
} Transform;

#define  Transform_elem( t, i, j ) ((t).m[j][i])


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/transforms.h
//
// As such, it needs the following Copyright notice
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
@VERSION    : $Header: /private-cvsroot/minc/volume_io/Include/volume_io/transforms.h,v 1.12.2.2 2005/03/31 17:39:49 bert Exp $
---------------------------------------------------------------------------- */

typedef  void   (*User_transform_function)( void  *user_data,
                                            double  x,
                                            double  y,
                                            double  z,
                                            double  *x_trans,
                                            double  *y_trans,
                                            double  *z_trans );

typedef enum { 
	LINEAR,
#if defined(BEVIN_ALL_TYPES_SUPPORTED)
	THIN_PLATE_SPLINE, 
	USER_TRANSFORM,
#endif
        CONCATENATED_TRANSFORM
#if defined(BEVIN_ALL_TYPES_SUPPORTED)
	,
	GRID_TRANSFORM 
#endif
} Transform_types;


typedef struct General_transform
{
    Transform_types             type;
    bool                    	inverse_flag;

    /* --- linear transform */

    Transform               *linear_transform;
    Transform               *inverse_linear_transform;

    /* --- non-linear transform */

    int                         n_points;
    int                         n_dimensions;
    double                    **points;
    double                    **displacements;   /* n_points + n_dim + 1 by */
                                                   /* n_dim */

    /* --- grid transform */

    void                        *displacement_volume;

    /* --- user_defined */

    void                        *user_data;
    size_t                      size_user_data;
    User_transform_function     user_transform_function;
    User_transform_function     user_inverse_transform_function;

    /* --- concatenated transform */

    int                         n_transforms;
    struct General_transform    *transforms;

} General_transform;


// The following is a replacement for some portions of
// mni/1.5/include/minc.h
//
// As such, it needs the following Copyright notice
//
/*
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/
#define MI_MAX_IMGDIMS 		100

#define MIxspace           "xspace"
#define MIyspace           "yspace"
#define MIzspace           "zspace"
#define MItime             "time"

/* NC_UNSPECIFIED is defined here for backwards compatibility. With 
   NetCDF 2.x, NC_UNSPECIFIED may already be defined either through a macro
   or an enum. In the latter case, this macro will override the enum. */
#ifndef NC_UNSPECIFIED
#  define NC_UNSPECIFIED MI_ORIGINAL_TYPE
#endif

#define MI_ORIGINAL_TYPE ((nc_type) 0)


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/volume.h
// As such, it needs the following Copyright notice
/*
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
*/
typedef  enum  { MNC_FORMAT, FREE_FORMAT }       Volume_file_formats;


typedef struct volume_struct
{
        int                     spatial_axes[VIO_N_DIMENSIONS];

        nc_type                 nc_data_type;
        double                  voxel_min;
        double                  voxel_max;

        bool                    real_range_set;
        double                  real_value_scale;
        double                  real_value_translation;
	
        double                  separations[VIO_MAX_DIMENSIONS];
	double			starts[VIO_MAX_DIMENSIONS];
        double                  direction_cosines[VIO_MAX_DIMENSIONS][VIO_N_DIMENSIONS];

        bool                    voxel_to_world_transform_uptodate;
	General_transform       voxel_to_world_transform;

        VIO_multidim_array      array;
        bool                    signed_flag;

        const char*             dimension_names[VIO_MAX_DIMENSIONS];
        bool                    is_rgba_data;
  
        const char*             coordinate_system_name;
  

//#define   BEVIN_ALL_VOLUME_MEMBERS
#if defined(BEVIN_ALL_VOLUME_MEMBERS)
	//
	// When moving one of these out, make sure the current conditionalized uses
	// have the condition removed!
	//
        bool                    is_cached_volume;
        //VIO_volume_cache_struct cache;
  
        double                 *irregular_starts[VIO_MAX_DIMENSIONS];
        double                 *irregular_widths[VIO_MAX_DIMENSIONS];
#endif

} volume_struct;

typedef volume_struct* Volume;


#include "minc_structures.h"

#define MI_UNKNOWN_SPACE        "unknown___"
#define ANY_SPATIAL_DIMENSION   "any_spatial_dimension"


typedef  struct
{
    int         arent_any_yet;
} volume_creation_options;

typedef  struct
{
    bool        promote_invalid_to_zero_flag;
    bool        convert_vector_to_scalar_flag;
    bool        convert_vector_to_colour_flag;
    int         dimension_size_for_colour_data;
    int         max_dimension_size_for_colour_data;
    int         rgba_indices[4];
    double      user_real_range[2];
} minc_input_options;

typedef  struct
{
    bool               file_is_being_read;

    /* input and output */

    int                cdfid;
    int                img_var;
    int                n_file_dimensions;
    long               sizes_in_file[MAX_VAR_DIMS];
    long               indices[MAX_VAR_DIMS];
    const char*        dim_names[MAX_VAR_DIMS];
    Volume             volume;
    int                to_volume_index[MAX_VAR_DIMS];
    int                to_file_index[VIO_MAX_DIMENSIONS];
    int                minc_icv;
    const char*        filename;

    /* input only */

    bool               end_volume_flag;
    bool               converting_to_colour;
    int                rgba_indices[4];
    int                n_volumes_in_file;

    int                valid_file_axes[VIO_MAX_DIMENSIONS];

    int                n_slab_dims;

    int                spatial_axes[VIO_N_DIMENSIONS];
    General_transform  voxel_to_world_transform;
    minc_input_options original_input_options;

    /* output only */

    int                img_var_id;
    int                min_id;
    int                max_id;
    double             image_range[2];
    bool               end_def_done;
    bool               ignoring_because_cached;
    bool               variables_written;
    int                dim_ids[MAX_VAR_DIMS];
    bool               outputting_in_order;
    bool               entire_file_written;
    nc_type            nc_data_type;
    bool               signed_flag;
    double             valid_range[2];
    int                image_dims[MAX_VAR_DIMS];
    int                src_cdfid;
    int                src_img_var;
} minc_file_struct;

typedef  minc_file_struct  *Minc_file;


typedef struct
{
      Volume_file_formats  file_format;
  
      Minc_file            minc_file;
  
      /* for free format files only */
  
      FILE                 *volume_file;
      int                  slice_index;
      long                 sizes_in_file[VIO_MAX_DIMENSIONS];
      int                  axis_index_from_file[VIO_MAX_DIMENSIONS];
      VIO_Data_types       file_data_type;
      bool                 one_file_per_slice;
      const char*          directory;
      const char*          *slice_filenames;
      int                  *slice_byte_offsets;
      unsigned char        *byte_slice_buffer;
      unsigned short       *short_slice_buffer;
} volume_input_struct;



typedef struct
{
    double global_image_range[2];
    char*  dimension_names[VIO_MAX_DIMENSIONS];
    bool   use_starts_set;
    bool   use_volume_starts_and_steps;
} minc_output_options;


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/vol_io_prototypes.h
// which did not have its own Copyright notice
//

void  make_identity_transform( Transform   *transform );

bool close_to_identity(
    Transform   *transform );
    
void   concat_transforms(
    Transform   *result,
    Transform   *t1,
    Transform   *t2 );

Transform* get_linear_transform_ptr(
    General_transform   *transform );

Transform* get_inverse_linear_transform_ptr(
    General_transform   *transform );

void concat_general_transforms(
    General_transform   *first,
    General_transform   *second,
    General_transform   *result );
    
void copy_general_transform(
    General_transform   *transform,
    General_transform   *copy );

void delete_general_transform(
    General_transform   *transform );

VIO_Status input_transform_file(
    const char* filename,
    General_transform   *transform );

void transform_point(
    Transform  *transform,
    double 	x,
    double	y,
    double 	z,
    double	*x_trans,
    double	*y_trans,
    double	*z_trans );


int get_volume_n_dimensions(
    Volume volume );

void get_volume_sizes(
    Volume 	volume,
    int      	sizes[] );

nc_type get_volume_nc_data_type(
    Volume      volume,
    bool*	signed_flag );

void get_volume_separations(
    Volume   	volume,
    double 	separations[] );

void convert_voxel_to_world(
    Volume   	volume,
    double     	voxel[],
    double     *x_world,
    double     *y_world,
    double     *z_world );

void delete_volume(
    Volume volume );

void delete_volume_input(
    volume_input_struct   *input_info );

Volume create_volume(
    int          n_dimensions,
    /*const*/ char*  dimension_names[],		// need compat with minc
    nc_type      nc_data_type,
    bool         signed_flag,
    double	 voxel_min,
    double       voxel_max );

void  set_volume_space_type(
    Volume   volume,
    const char* name );
    
void set_volume_sizes(
    Volume   	volume,
    int         sizes[] );

void alloc_volume_data(
    Volume      volume );

void  free_volume_data(
    Volume   volume );

bool volume_is_alloced(
    Volume   volume );

void set_volume_separations(
    Volume      volume,
    double      separations[] );

void set_volume_direction_unit_cosine(
    Volume   volume,
    int      axis,
    double   dir[] );

void set_volume_starts(
    Volume  volume,
    double  starts[] );

int set_volume_irregular_starts(Volume volume, int idim, int count, double *starts);
int set_volume_irregular_widths(Volume volume, int idim, int count, double *widths);

void  set_volume_voxel_range(
    Volume volume,
    double voxel_min,
    double voxel_max );
    
void  get_volume_voxel_range(
    Volume volume,
    double *voxel_min,
    double *voxel_max );

void  set_volume_real_range(
    Volume volume,
    double real_min,
    double real_max );

void  get_volume_real_range(
    Volume     volume,
    double       *min_value,
    double       *max_value );
       
void  set_volume_type(
    Volume       volume,
    nc_type      nc_data_type,
    bool         signed_flag,
    double       voxel_min,
    double       voxel_max );

void set_volume_direction_cosine(
    Volume   	volume,
    int      	axis,
    double	dir[] );

void set_volume_translation(
    Volume  	volume,
    double    	voxel[],
    double    	world_space_voxel_maps_to[] );

void set_volume_voxel_value(
    Volume      volume,
    int         v0,
    int         v1,
    int         v2,
    int         v3,
    int         v4,
    double      voxel );


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
    volume_input_struct *input_info );

bool input_more_of_volume(
    Volume                volume,
    volume_input_struct  *input_info,
    double               *fraction_done );

double get_volume_voxel_value(
    Volume   volume,
    int      v0,
    int      v1,
    int      v2,
    int      v3,
    int      v4 );
    
General_transform* get_voxel_to_world_transform(
    Volume   volume );

VIO_Status  output_volume(
    const char*		  filename,
    nc_type		  file_nc_data_type,
    bool              	  file_signed_flag,
    double                file_voxel_min,
    double                file_voxel_max,
    Volume                volume,
    const char*	  	  history,
    minc_output_options  *options );

#endif // MINC_VOLUME_IO_H
