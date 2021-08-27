/**
 * @brief Wrapper for MNI's volume_io.h, to decouple from MNI lib
 *
 */
/*
 * Original Author: Bevin Brett, extracted and enhanced from various minc files
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

#ifndef MINC_INTERNALS_H
#define MINC_INTERNALS_H

#include "minc_volume_io.h"

enum {X=0, Y=1, Z=2};


// The following is a replacement for some portions of
// mni/1.5/include/volume_io/basic.h
//                 volume_io/progress.h
//	   volume_io/Prog_utils/time.c
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

double current_realtime_seconds( );

char* format_time(		// free() the result
    const char* format,
    double 	seconds );

    
typedef  struct
{
    bool       force_one_line;
    bool       first_msg_displayed;
    bool       one_line_flag;
    int        n_steps;
    int        n_dots_so_far;
    int        total_n_dots;
    double     start_time;
    double     previous_time;
    double     update_rate;
    double     sum_xy;
    double     sum_xx;
    const char* title;

    double     last_check_time;
    int        check_every;
    int        next_check_step;
    int        last_check_step;
} VIO_progress_struct;


void initialize_progress_report(
    VIO_progress_struct   *progress,
    bool              one_line_only,
    int               n_steps,
    const char*       title );

void  update_progress_report(
    VIO_progress_struct   *progress,
    int               current_step );

void  terminate_progress_report(
    VIO_progress_struct   *progress );
        
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
#define MI_EMPTY_STRING ""

#define MI_NOERROR 0
#define MI_ERROR (-1)

#define MI_ERR_BAD_STDVAR       1336  /* Not a standard variable */

#define MI_NUM_SPACE_DIMS 3


/* NetCDF standard attributes */
#define MIunits       "units"
#define MIlong_name   "long_name"
#define MIvalid_range "valid_range"
#define MIvalid_max   "valid_max"
#define MIvalid_min   "valid_min"
#define MI_FillValue  "_FillValue"
#define MItitle       "title"
#define MIhistory     "history"

/* General variable attributes */
#define MIvartype  "vartype"
#define MIvarid    "varid"
#define MIsigntype "signtype"
#define MIparent   "parent"
#define MIchildren "children"
#define MIcomments "comments"
#define MIversion  "version"

/* General attribute constants */
/*    Prefix for identifying a variable attribute pointer */
#define MI_VARATT_POINTER_PREFIX "--->"
/*    Separator for elements of MIchildren */
#define MI_CHILD_SEPARATOR "\n"
/*    MIvartype values */
#define MI_GROUP     "group________"
#define MI_DIMENSION "dimension____"
#define MI_DIM_WIDTH "dim-width____"
#define MI_VARATT    "var_attribute"
/*    MIvarid value */
#define MI_STDVAR "MINC standard variable"
/*    MIsigntype values */
#define MI_SIGNED   "signed__"
#define MI_UNSIGNED "unsigned"
/*    MIversion value */
#define MI_VERSION_1_0 "MINC Version    1.0"
#define MI_CURRENT_VERSION MI_VERSION_1_0
/* Generally useful values for boolean attributes */
#define MI_TRUE  "true_"
#define MI_FALSE "false"

/* Dimension names and names of associated variables */
#define MItfrequency       "tfrequency"
#define MIxfrequency       "xfrequency"
#define MIyfrequency       "yfrequency"
#define MIzfrequency       "zfrequency"
#define MIvector_dimension "vector_dimension"
#define MIxspace_width     "xspace-width"
#define MIyspace_width     "yspace-width"
#define MIzspace_width     "zspace-width"
#define MItime_width       "time-width"
#define MItfrequency_width "tfrequency-width"
#define MIxfrequency_width "xfrequency-width"
#define MIyfrequency_width "yfrequency-width"
#define MIzfrequency_width "zfrequency-width"

/* Dimension variable attribute names */
/* For dimension variables (MIspacing is also for dimension width vars) */
#define MIspacing           "spacing"
#define MIstep              "step"
#define MIstart             "start"
#define MIspacetype         "spacetype"
#define MIalignment         "alignment"
#define MIdirection_cosines "direction_cosines"
/* For dimension width variables */
#define MIwidth             "width"
#define MIfiltertype        "filtertype"

/* Dimension attribute constants */
/*    MIgridtype values */
#define MI_REGULAR   "regular__"
#define MI_IRREGULAR "irregular"
/*    MIspacetype values */
#define MI_NATIVE    "native____"
#define MI_TALAIRACH "talairach_"
#define MI_CALLOSAL  "callosal__"
/*    MIalignment values */
#define MI_START  "start_"
#define MI_CENTRE "centre"
#define MI_END    "end___"
#define MI_CENTER MI_CENTRE
/*    MIfiltertype values */
#define MI_SQUARE     "square____"
#define MI_GAUSSIAN   "gaussian__"
#define MI_TRIANGULAR "triangular"

/* The root variable */
#define MIrootvariable "rootvariable"

/* The image variable and its attributes */
#define MIimage    "image"
#define MIimagemax "image-max"
#define MIimagemin "image-min"
#define MIcomplete "complete"

/* The patient variable and its attributes */
#define MIpatient        "patient"
#define MIfull_name      "full_name"
#define MIother_names    "other_names"
#define MIidentification "identification"
#define MIother_ids      "other_ids"
#define MIbirthdate      "birthdate"
#define MIsex            "sex"
#define MIage            "age"
#define MIweight         "weight"
#define MIsize           "size"
#define MIaddress        "address"
#define MIinsurance_id   "insurance_id"

/* Patient attribute constants */
#define MI_MALE   "male__"
#define MI_FEMALE "female"
#define MI_OTHER  "other_"

/* The study variable and its attributes */
#define MIstudy               "study"
#define MIstart_time          "start_time"
#define MIstart_year          "start_year"
#define MIstart_month         "start_month"
#define MIstart_day           "start_day"
#define MIstart_hour          "start_hour"
#define MIstart_minute        "start_minute"
#define MIstart_seconds       "start_seconds"
#define MImodality            "modality"
#define MImanufacturer        "manufacturer"
#define MIdevice_model        "device_model"
#define MIinstitution         "institution"
#define MIdepartment          "department"
#define MIstation_id          "station_id"
#define MIreferring_physician "referring_physician"
#define MIattending_physician "attending_physician"
#define MIradiologist         "radiologist"
#define MIoperator            "operator"
#define MIadmitting_diagnosis "admitting_diagnosis"
#define MIprocedure           "procedure"
#define MIstudy_id            "study_id"

/* Study attribute constants */
#define MI_PET   "PET__"
#define MI_SPECT "SPECT"
#define MI_GAMMA "GAMMA"
#define MI_MRI   "MRI__"
#define MI_MRS   "MRS__"
#define MI_MRA   "MRA__"
#define MI_CT    "CT___"
#define MI_DSA   "DSA__"
#define MI_DR    "DR___"
#define MI_LABEL "label"

/* The acquisition variable and its attributes */
#define MIacquisition           "acquisition"
#define MIprotocol              "protocol"
#define MIscanning_sequence     "scanning_sequence"
#define MIrepetition_time       "repetition_time"
#define MIecho_time             "echo_time"
#define MIinversion_time        "inversion_time"
#define MInum_averages          "num_averages"
#define MIimaging_frequency     "imaging_frequency"
#define MIimaged_nucleus        "imaged_nucleus"
#define MIradionuclide          "radionuclide"
#define MIcontrast_agent        "contrast_agent"
#define MIradionuclide_halflife "radionuclide_halflife"
#define MItracer                "tracer"
#define MIinjection_time        "injection_time"
#define MIinjection_year        "injection_year"
#define MIinjection_month       "injection_month"
#define MIinjection_day         "injection_day"
#define MIinjection_hour        "injection_hour"
#define MIinjection_minute      "injection_minute"
#define MIinjection_seconds     "injection_seconds"
#define MIinjection_length      "injection_length"
#define MIinjection_dose        "injection_dose"
#define MIdose_units            "dose_units"
#define MIinjection_volume      "injection_volume"
#define MIinjection_route       "injection_route"






#define MItime    "time"
#define MIimage   "image"

#define MI_DEFAULT_MAX 1.0
#define MI_DEFAULT_MIN 0.0

#define MI_MAX_ATTSTR_LEN  64






/* For converting data type */
#define MI_ICV_TYPE             1
#define MI_ICV_SIGN             2
#define MI_ICV_DO_RANGE         3
#define MI_ICV_VALID_MAX        4
#define MI_ICV_VALID_MIN        5
/* For doing normalization */
#define MI_ICV_DO_NORM          6
#define MI_ICV_USER_NORM        7
#define MI_ICV_IMAGE_MAX        8
#define MI_ICV_IMAGE_MIN        9
/* Values actually used in normalization - read-only */
#define MI_ICV_NORM_MAX        10
#define MI_ICV_NORM_MIN        11
/* For doing dimension conversions */
#define MI_ICV_DO_DIM_CONV     12
/* For converting vector fields to scalar */
#define MI_ICV_DO_SCALAR       13
/* For flipping axis direction */
#define MI_ICV_XDIM_DIR        14
#define MI_ICV_YDIM_DIR        15
#define MI_ICV_ZDIM_DIR        16
/* For changing size of first two dimensions (excluding MIvector_dimension) */
#define MI_ICV_ADIM_SIZE       17
#define MI_ICV_BDIM_SIZE       18
#define MI_ICV_KEEP_ASPECT     19
/* The pixel size and location of first two dimensions (these are readonly) */
#define MI_ICV_ADIM_STEP       20
#define MI_ICV_BDIM_STEP       21
#define MI_ICV_ADIM_START      22
#define MI_ICV_BDIM_START      23
/* Number of image dimensions for dimension conversion */
#define MI_ICV_NUM_IMGDIMS     24
/* Number of dimensions of image variable taking into account vector/scalar
   data (read-only property) */
#define MI_ICV_NUM_DIMS        25
/* Id of file and image variable (read-only properties) */
#define MI_ICV_CDFID           26
#define MI_ICV_VARID           27
/* Names of MIimagemax and MIimagemin variables */
#define MI_ICV_MAXVAR          28
#define MI_ICV_MINVAR          29
/* For setting input values to a specified fillvalue */
#define MI_ICV_DO_FILLVALUE    30
#define MI_ICV_FILLVALUE       31
/* Image dimension properties. For each dimension, add the dimension 
   number (counting from fastest to slowest). */
#define MI_ICV_DIM_SIZE        1000
#define MI_ICV_DIM_STEP        1100
#define MI_ICV_DIM_START       1200


/* Constants that can be used as values for the above properties. */
/* Possible values for MI_ICV_?DIM_DIR */
#define MI_ICV_POSITIVE         1
#define MI_ICV_NEGATIVE       (-1)
#define MI_ICV_ANYDIR           0
/* Possible value for MI_ICV_?DIM_SIZE */
#define MI_ICV_ANYSIZE        (-1)


#define MI_MAX_NUM_ICV 		32
#define MI_ICV_NUM_IMGDIMS      24

int miicv_attach(int icvid, int cdfid, int varid);
int miicv_create(void);
int miicv_detach(int icvid);
int miicv_free(int icvid);
int miicv_inqdbl(int icvid, int icv_property, double *value);
int miicv_setdbl(int icvid, int icv_property, double value);
int miicv_setint(int icvid, int icv_property, int value);
int miicv_setstr(int icvid, int icv_property, char *value);
int miicv_get(int icvid, long start[], long count[], void *values);
int miicv_put(int icvid, long start[], long count[], void *values);

int miget_datatype(int cdfid, int imgid, 
                          nc_type *datatype, int *is_signed);
			  
int miget_default_range(nc_type datatype, int is_signed, 
                               double default_range[]);
int miget_image_range(int cdfid, double image_range[]);
int micreate_std_variable(int cdfid, const char *name, nc_type datatype, 
                                 int ndims, int dim[]);
int miset_valid_range(int cdfid, int imgid, double valid_range[]);
long *miset_coords(int nvals, long value, long coords[]);
long *mitranslate_coords(int cdfid, 
                                int invar,  long incoords[],
                                int outvar, long outcoords[]);

int mivarget1(int cdfid, int varid, long mindex[],
                     nc_type datatype, char *sign, void *value);				int miattget_with_sign(int cdfid, int varid, char *name, 
                              char *insign, nc_type datatype, char *outsign,
                              int max_length, void *value, int *att_length);
int miattget(int cdfid, int varid, char *name, nc_type datatype,
                    int max_length, void *value, int *att_length);
int miattget1(int cdfid, int varid, char *name, nc_type datatype,
                     void *value);
int mivarget(int cdfid, int varid, long start[], long count[],
                    nc_type datatype, char *sign, void *values);
char *miattgetstr(int cdfid, int varid, char *name,
                         int maxlen, char *value);
int miget_valid_range(int cdfid, int imgid, double valid_range[]);			     

// The following is a replacement for some portions of
// mni/1.5/include/volume_io/vol_io_prototypes.h
// which did not have its own Copyright notice
//
const char* get_date();

VIO_Colour  make_rgba_Colour_0_1(
    double r,
    double g,
    double b,
    double a );

void   create_empty_multidim_array(
    VIO_multidim_array* array,
    int                 n_dimensions,
    VIO_Data_types      data_type );

void   create_multidim_array(
    VIO_multidim_array  *array,
    int             n_dimensions,
    int             sizes[],
    VIO_Data_types      data_type );

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
    bool                use_src_order );
    
int  get_multidim_n_dimensions(
    VIO_multidim_array   *array );

void  get_multidim_sizes(
    VIO_multidim_array   *array,
    int              sizes[] );
       
void  delete_multidim_array(
    VIO_multidim_array   *array );

void  alloc_multidim_array(
    VIO_multidim_array   *array );
    
bool multidim_array_is_alloced(
    VIO_multidim_array   *array );

void  set_multidim_data_type(
    VIO_multidim_array       *array,
    VIO_Data_types        data_type );

VIO_Data_types  get_multidim_data_type(
    VIO_multidim_array       *array );

int  get_type_size(
    VIO_Data_types   type );

void  set_multidim_sizes(
    VIO_multidim_array   *array,
    int              sizes[] );
    
void compute_world_transform(
    int                 spatial_axes[VIO_N_DIMENSIONS],
    double              separations[],
    double              direction_cosines[][VIO_N_DIMENSIONS],
    double              starts[],
    General_transform   *world_transform );

bool convert_dim_name_to_spatial_axis(
    const char* name,
    int     *axis );

void  convert_transform_to_starts_and_steps(
    General_transform  *transform,
    int                n_volume_dimensions,
    double             step_signs[],
    int                spatial_axes[],
    double             starts[],
    double             steps[],
    double             dir_cosines[][VIO_N_DIMENSIONS] );

void  get_transform_origin_real(
    Transform   *transform,
    double       origin[] );
    
void  get_transform_x_axis_real(
    Transform   *transform,
    double       x_axis[] );
    
void  get_transform_y_axis_real(
    Transform   *transform,
    double       y_axis[] );
    
void  get_transform_z_axis_real(
    Transform   *transform,
    double       z_axis[] );

void  set_rgb_volume_flag(
    Volume   volume,
    bool     flag );

char** get_volume_dimension_names(
    Volume   volume );

void  delete_dimension_names(
    Volume   volume,
    char*    dimension_names[] );
    
VIO_Data_types  get_volume_data_type(
    Volume       volume );

void  set_default_minc_output_options(
    minc_output_options  *options           );

Minc_file  initialize_minc_output(
    const char*            filename,
    int                    n_dimensions,
    const char*            dim_names[],
    int                    sizes[],
    nc_type                file_nc_data_type,
    bool                   file_signed_flag,
    double                 file_voxel_min,
    double                 file_voxel_max,
    General_transform      *voxel_to_world_transform,
    Volume                 volume_to_attach,
    minc_output_options    *options );

VIO_Status  copy_auxiliary_data_from_minc_file(
    Minc_file   file,
    const char* filename,
    const char* history_string );

VIO_Status  add_minc_history(
    Minc_file   file,
    const char *history_string );

void  set_minc_output_real_range(
    minc_output_options  *options,
    double                real_min,
    double                real_max );

void  set_minc_output_use_volume_starts_and_steps_flag(
    minc_output_options  *options,
    bool                  flag );

VIO_Status  output_minc_volume(
    Minc_file   file );

VIO_Status  close_minc_output(
    Minc_file   file );  
         
// The following are in minc_files.c
//
bool file_exists(
    const char* expandedFilename );
    
char* file_exists_as_compressed(	// returns the pathname of the compressed file, use free() to delete
    const char* expandedFilename );

bool filename_extension_matches(
    const char*   expanded_filename_possibly_with_z,
    const char*   extension );
    
  
// The following are in minc_input_free.c
//
VIO_Status  initialize_free_format_input(
    const char*          filename,
    Volume               volume,
    volume_input_struct  *volume_input );
    
void  delete_free_format_input(
    volume_input_struct   *volume_input );

bool input_more_free_format_file(
    Volume                volume,
    volume_input_struct  *volume_input,
    double               *fraction_done );

// The following are in minc_input_mnc.c
//
Minc_file  initialize_minc_input(
    const char*          filename,
    Volume               volume,
    minc_input_options   *options );
    
int   get_minc_file_n_dimensions(
    const char*   filename );
    
VIO_Status  close_minc_input(
    Minc_file   file );

bool input_more_minc_file(
    Minc_file   file,
    double     *fraction_done );
    

// The following is a replacement for some portions of
// minc-1.5.1/libsrc/minc_basic.h
// As such, it needs the following Copyright notice
/*
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
#define MI_MAX_VAR_BUFFER_SIZE 10000

#define MI_WIDTH_SUFFIX "-width"

// The following is a replacement for some portions of
// minc-1.5.1/libsrc/minc_routines.h
// As such, it needs the following Copyright notice
/*
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

int MI_get_sign_from_string(nc_type datatype, char *sign);

int MI_convert_type(long number_of_values,
                                nc_type intype,  int insign,  void *invalues,
                                nc_type outtype, int outsign, void *outvalues,
                                mi_icv_type *icvp);
				
int MI_var_loop(int ndims, long start[], long count[],
                            int value_size, int *bufsize_step,
                            long max_buffer_size,
                            void *caller_data,
                            int (*action_func) (int, long [], long [], 
                                                long, void *, void *));


int MI_varaccess(int operation, int cdfid, int varid, 
                             long start[], long count[],
                             nc_type datatype, int sign, void *values,
                             int *bufsize_step, mi_icv_type *icvp);
			     						
typedef double Double4x4[4*4];
#define Index4x4(I,J) (4*(I)+(J))
bool invert_4x4_matrix( 
	const Double4x4 * mat, 	// doubles are read, not changed
	Double4x4 * 	  inv );	// doubles are written

#endif
