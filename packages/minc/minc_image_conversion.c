/* ----------------------------- MNI Header -----------------------------------
@NAME       : image_conversion.c
@DESCRIPTION: File of functions to manipulate image conversion variables
              (icv). These variables allow conversion of netcdf variables
              (the MINC image variable, in particular) to a form more
              convenient for a program.
@METHOD     : Routines included in this file :
              public :
                 miicv_create
                 miicv_free
                 miicv_setdbl
                 miicv_setint
                 miicv_setlong
                 miicv_setstr
                 miicv_inqdbl
                 miicv_inqint
                 miicv_inqlong
                 miicv_inqstr
                 miicv_ndattach
                 miicv_detach
                 miicv_get
                 miicv_put
              semiprivate :
                 MI_icv_chkid
              private :
                 MI_icv_get_type
                 MI_icv_get_vrange
                 MI_get_default_range
                 MI_icv_get_norm
                 MI_icv_access
                 MI_icv_zero_buffer
                 MI_icv_coords_tovar
                 MI_icv_calc_scale
@CREATED    : July 27, 1992. (Peter Neelin, Montreal Neurological Institute)
@MODIFIED   : 
 * $Log: image_conversion.c,v $
 * Revision 6.10.2.1  2004/09/28 20:23:40  bert
 * Minor portability fixes for Windows
 *
 * Revision 6.10  2003/09/18 16:17:00  bert
 * Correctly cast double to nc_type
 *
 * Revision 6.9  2001/11/28 15:38:07  neelin
 * Removed limit on number of icvs that can exist at one time.
 *
 * Revision 6.8  2001/11/13 21:00:24  neelin
 * Modified icv scaling calculations for no normalization. When the icv
 * type is double, normalization is always done, regardless of the
 * normalization setting. When the external type is floating point,
 * normalization to the slice real range is done (essentially a valid
 * range scaling, but where the valid range for a float is the slice real
 * range).
 *
 * Revision 6.7  2001/11/13 14:15:17  neelin
 * Added functions miget_image_range and mivar_exists
 *
 * Revision 6.6  2001/08/20 13:16:53  neelin
 * Removed extraneous variables from MI_icv_get_vrange.
 *
 * Revision 6.5  2001/08/16 19:24:11  neelin
 * Fixes to the code handling valid_range values.
 *
 * Revision 6.4  2001/08/16 16:41:31  neelin
 * Added library functions to handle reading of datatype, sign and valid range,
 * plus writing of valid range and setting of default ranges. These functions
 * properly handle differences between valid_range type and image type. Such
 * difference can cause valid data to appear as invalid when double to float
 * conversion causes rounding in the wrong direction (out of range).
 * Modified voxel_loop, volume_io and programs to use these functions.
 *
 * Revision 6.3  2001/08/16 13:32:18  neelin
 * Partial fix for valid_range of different type from image (problems
 * arising from double to float conversion/rounding). NOT COMPLETE.
 *
 * Revision 6.2  2001/04/17 18:40:12  neelin
 * Modifications to work with NetCDF 3.x
 * In particular, changed NC_LONG to NC_INT (and corresponding longs to ints).
 * Changed NC_UNSPECIFIED to NC_NAT.
 * A few fixes to the configure script.
 *
 * Revision 6.1  1999/10/19 14:45:07  neelin
 * Fixed Log subsitutions for CVS
 *
 * Revision 6.0  1997/09/12 13:24:54  neelin
 * Release of minc version 0.6
 *
 * Revision 5.0  1997/08/21  13:25:53  neelin
 * Release of minc version 0.5
 *
 * Revision 4.0  1997/05/07  20:07:52  neelin
 * Release of minc version 0.4
 *
 * Revision 3.3  1997/04/21  17:32:04  neelin
 * Fixed calculation of scale for icv so that values are not re-scaled
 * from real values to file floating-point values.
 *
 * Revision 3.2  1997/04/10  19:22:18  neelin
 * Removed redefinition of NULL and added pointer casts in appropriate places.
 *
 * Revision 3.1  1997/04/10  18:14:50  neelin
 * Fixed handling of invalid data when icv scale is zero.
 *
 * Revision 3.0  1995/05/15  19:33:12  neelin
 * Release of minc version 0.3
 *
 * Revision 2.3  1995/02/08  19:14:44  neelin
 * More changes for irix 5 lint.
 *
 * Revision 2.2  1995/02/08  19:01:06  neelin
 * Moved private function declarations from minc_routines.h to appropriate file.
 *
 * Revision 2.1  1994/12/09  09:12:30  neelin
 * Added test in miicv_detach to make sure that icv is attached before
 * detaching it.
 *
 * Revision 2.0  94/09/28  10:37:55  neelin
 * Release of minc version 0.2
 * 
 * Revision 1.18  94/09/28  10:37:06  neelin
 * Pre-release
 * 
 * Revision 1.17  93/08/11  12:59:31  neelin
 * We need only increment the chunk pointer (see previous fix) if we are
 * not doing dimension conversion (dimension conversion handles the 
 * offsets itself).
 * 
 * Revision 1.16  93/08/11  11:49:36  neelin
 * Added RCS logging in source.
 * Fixed bug in MI_icv_access so that pointer to values buffer is incremented
 * as we loop through the chunks. This affected calls to miicv_get/put that
 * had MIimagemax/min varying over the values read in one call (ie. reading
 * or writing a volume with MIimagemax/min varying over slices will give
 * incorrect results if the volume is read with one call).
 * 
              January 22, 1993 (P.N.)
                 - Modified handling of icv properties with miicv_set<type>.
                   Removed routine miicv_set. Use routines miicv_setdbl,
                   miicv_setint, miicv_setlong, miicv_setstr instead (this
                   gives type checking at compile time).
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
---------------------------------------------------------------------------- */


#ifndef lint
//static char rcsid[] = "$Header: /private-cvsroot/minc/libsrc/image_conversion.c,v 6.10.2.1 2004/09/28 20:23:40 bert Exp $ MINC (MNI)";
#endif

//BEVIN #include "minc_private.h"
//BEVIN #include "type_limits.h"
#include <limits.h>	//BEVIN 
#include <float.h>	//BEVIN 

//++
// BEVIN definitions to avoid big changes below
//
#include "minc_internals.h"
#include "minc_basic.h"
#include "minc_structures.h"

void* minc_image_conversion_dummy = &longMIN;	// just to avoid an error message

#include <string.h>
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <fcntl.h>
//#include <unistd.h>


#define MNCAPI 

#define SEMIPRIVATE static
#define PRIVATE static
#define FALSE 0
#define TRUE 1

static const char* saved_routine_name;
#define MI_SAVE_ROUTINE_NAME(X) (saved_routine_name = (X))
#define MALLOC(COUNT,TYPE)      ((TYPE*)malloc (      (COUNT)*sizeof(TYPE)))
#define REALLOC(PTR,COUNT,TYPE) ((TYPE*)realloc((PTR),(COUNT)*sizeof(TYPE)))
#define FREE(X) free((void*)(X))
#define MI_RETURN(X) return (X)
#define MI_RETURN_ERROR(X) return (X)
#define MI_LOG_PKG_ERROR2(CODE,MSG) { fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, (MSG)); exit(1); }
#define MI_LOG_SYS_ERROR1(MSG)      { fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, (MSG)); exit(1); }

#define MI_CHK_ERR(expr) { if ((expr)==MI_ERROR) MI_RETURN_ERROR(MI_ERROR); } 

SEMIPRIVATE mi_icv_type *MI_icv_chkid(int icvid);

//
//--

/* Private functions */

PRIVATE int MI_icv_get_type(mi_icv_type *icvp, int cdfid, int varid);
PRIVATE int MI_icv_get_vrange(mi_icv_type *icvp, int cdfid, int varid);
PRIVATE double MI_get_default_range(char *what, nc_type datatype, int sign);
PRIVATE int MI_icv_get_norm(mi_icv_type *icvp, int cdfid, int varid);
PRIVATE int MI_icv_access(int operation, mi_icv_type *icvp, long start[], 
                          long count[], void *values);
PRIVATE int MI_icv_zero_buffer(mi_icv_type *icvp, long count[], void *values);
PRIVATE int MI_icv_coords_tovar(mi_icv_type *icvp, 
                                long icv_start[], long icv_count[],
                                long var_start[], long var_count[]);
PRIVATE int MI_icv_calc_scale(int operation, mi_icv_type *icvp, long coords[]);

/* Array of pointers to image conversion structures */
static int minc_icv_list_nalloc = 0;
static mi_icv_type **minc_icv_list = NULL;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_create
@INPUT      : (none)
@OUTPUT     : (none)
@RETURNS    : icv id or MI_ERROR when an error occurs
@DESCRIPTION: Creates an image conversion variable (icv) and returns
              a handle to it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_create()
{
   int new_icv;       /* Id of newly created icv */
   mi_icv_type *icvp;  /* Pointer to new icv structure */
   int idim;
   int new_nalloc;

   MI_SAVE_ROUTINE_NAME("miicv_create");

   /* Look for free slot */
   for (new_icv=0; new_icv<minc_icv_list_nalloc; new_icv++)
      if (minc_icv_list[new_icv]==NULL) break;

   /* If none, then extend the list */
   if (new_icv>=minc_icv_list_nalloc) {

      /* How much space will be needed? */
      new_nalloc = minc_icv_list_nalloc + MI_MAX_NUM_ICV;

      /* Check for first allocation */
      if (minc_icv_list_nalloc == 0) {
         minc_icv_list = MALLOC(new_nalloc, mi_icv_type *);
      }
      else {
         minc_icv_list = REALLOC(minc_icv_list, new_nalloc, mi_icv_type *);
      }

      /* Check that the allocation was successful */
      if (minc_icv_list == NULL) {
         MI_LOG_SYS_ERROR1("miicv_create");
         MI_RETURN_ERROR(MI_ERROR);
      }
      /* Put in NULL pointers */
      for (new_icv=minc_icv_list_nalloc; new_icv<new_nalloc; new_icv++)
         minc_icv_list[new_icv] = NULL;

      /* Use the first free slot and update the list length */
      new_icv = minc_icv_list_nalloc;
      minc_icv_list_nalloc = new_nalloc;

   }

   /* Allocate a new structure */
   if ((minc_icv_list[new_icv]=MALLOC(1, mi_icv_type))==NULL) {
      MI_LOG_SYS_ERROR1("miicv_create");
      MI_RETURN_ERROR(MI_ERROR);
   }
   icvp=minc_icv_list[new_icv];

   /* Fill in defaults */

   /* Stuff for calling MI_varaccess */
   icvp->do_scale = FALSE;
   icvp->do_dimconvert = FALSE;
   icvp->do_fillvalue = FALSE;
   icvp->fill_valid_min = -DBL_MAX;
   icvp->fill_valid_max = DBL_MAX;

   /* User defaults */
   icvp->user_type = NC_SHORT;
   icvp->user_typelen = nctypelen(icvp->user_type);
   icvp->user_sign = MI_PRIV_SIGNED;
   icvp->user_do_range = TRUE;
   icvp->user_vmax = MI_get_default_range(MIvalid_max, icvp->user_type,
                                            icvp->user_sign);
   icvp->user_vmin = MI_get_default_range(MIvalid_min, icvp->user_type,
                                            icvp->user_sign);
   icvp->user_do_norm = FALSE;
   icvp->user_user_norm = FALSE;
   icvp->user_maxvar = strdup(MIimagemax);
   icvp->user_minvar = strdup(MIimagemin);
   icvp->user_imgmax = MI_DEFAULT_MAX;
   icvp->user_imgmin = MI_DEFAULT_MIN;
   icvp->user_do_dimconv = FALSE;
   icvp->user_do_scalar = TRUE;
   icvp->user_xdim_dir = MI_ICV_POSITIVE;
   icvp->user_ydim_dir = MI_ICV_POSITIVE;
   icvp->user_zdim_dir = MI_ICV_POSITIVE;
   icvp->user_num_imgdims = 2;
   icvp->user_keep_aspect = TRUE;
   icvp->user_do_fillvalue = FALSE;
   icvp->user_fillvalue = -DBL_MAX;
   for (idim=0; idim<MI_MAX_IMGDIMS; idim++) {
      icvp->user_dim_size[idim]=MI_ICV_ANYSIZE;
   }

   /* Variable values */
   icvp->cdfid = MI_ERROR;            /* Set so that we can recognise an */
   icvp->varid = MI_ERROR;            /* unattached icv */

   /* Values that can be read by user */
   icvp->derv_imgmax = MI_DEFAULT_MAX;
   icvp->derv_imgmin = MI_DEFAULT_MIN;
   for (idim=0; idim<MI_MAX_IMGDIMS; idim++) {
      icvp->derv_dim_step[idim] = 0.0;
      icvp->derv_dim_start[idim] = 0.0;
   }

   MI_RETURN(new_icv);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_free
@INPUT      : icvid
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Frees the image conversion variable (icv)
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_free(int icvid)
{
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_free");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Detach the icv if it is attached */
   if (icvp->cdfid != MI_ERROR)
      {MI_CHK_ERR(miicv_detach(icvid))}

   /* Free anything allocated at creation time */
   FREE(icvp->user_maxvar);
   FREE(icvp->user_minvar);

   /* Free the structure */
   FREE(icvp);
   minc_icv_list[icvid]=NULL;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_setdbl
@INPUT      : icvid        - icv id
              icv_property - property of icv to set
              value        - value to set it to
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets a property of an icv to a given double value
              Properties cannot be modified while the icv is attached to a 
              cdf file and variable (see miicv_attach and miicv_detach).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_setdbl(int icvid, int icv_property, double value)
{
   int ival, idim;
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_setdbl");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Check that the icv is not attached to a file */
   if (icvp->cdfid != MI_ERROR) {
      MI_LOG_PKG_ERROR2(MI_ERR_ICVATTACHED, 
         "Attempt to modify an attached image conversion variable");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Set the property */
   switch (icv_property) {
   case MI_ICV_TYPE:
      icvp->user_type   = (nc_type) (int) value;
      icvp->user_typelen= nctypelen(icvp->user_type);
      icvp->user_vmax   = MI_get_default_range(MIvalid_max, icvp->user_type,
                                               icvp->user_sign);
      icvp->user_vmin   = MI_get_default_range(MIvalid_min, icvp->user_type,
                                               icvp->user_sign);
      break;
   case MI_ICV_DO_RANGE:
      icvp->user_do_range = value; break;
   case MI_ICV_VALID_MAX:
      icvp->user_vmax   = value; break;
   case MI_ICV_VALID_MIN:
      icvp->user_vmin   = value; break;
   case MI_ICV_DO_NORM:
      icvp->user_do_norm = value; break;
   case MI_ICV_USER_NORM:
      icvp->user_user_norm = value; break;
   case MI_ICV_IMAGE_MAX:
      icvp->user_imgmax = value; break;
   case MI_ICV_IMAGE_MIN:
      icvp->user_imgmin = value; break;
   case MI_ICV_DO_FILLVALUE:
      icvp->user_do_fillvalue = value; break;
   case MI_ICV_FILLVALUE:
      icvp->user_fillvalue = value; break;
   case MI_ICV_DO_DIM_CONV:
      icvp->user_do_dimconv = value; break;
   case MI_ICV_DO_SCALAR:
      icvp->user_do_scalar = value; break;
   case MI_ICV_XDIM_DIR: 
      ival = value;
      icvp->user_xdim_dir = ((ival==MI_ICV_POSITIVE) || 
                             (ival==MI_ICV_NEGATIVE)) ? ival : MI_ICV_ANYDIR;
      break;
   case MI_ICV_YDIM_DIR:
      ival = value;
      icvp->user_ydim_dir = ((ival==MI_ICV_POSITIVE) || 
                             (ival==MI_ICV_NEGATIVE)) ? ival : MI_ICV_ANYDIR;
      break;
   case MI_ICV_ZDIM_DIR:
      ival = value;
      icvp->user_zdim_dir = ((ival==MI_ICV_POSITIVE) || 
                             (ival==MI_ICV_NEGATIVE)) ? ival : MI_ICV_ANYDIR;
      break;
   case MI_ICV_NUM_IMGDIMS:
      ival = value;
      if ((ival<0) || (ival>MI_MAX_IMGDIMS)) {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "MI_ICV_NUM_IMGDIMS out of range");
         MI_RETURN_ERROR(MI_ERROR);
      }
      icvp->user_num_imgdims = ival;
      break;
   case MI_ICV_ADIM_SIZE:
      icvp->user_dim_size[0] = value; break;
   case MI_ICV_BDIM_SIZE:
      icvp->user_dim_size[1] = value; break;
   case MI_ICV_KEEP_ASPECT:
      icvp->user_keep_aspect = value; break;
   case MI_ICV_SIGN:
   case MI_ICV_MAXVAR:
   case MI_ICV_MINVAR:
      MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                        "Tried to set icv string property to a number");
      MI_RETURN_ERROR(MI_ERROR);
      break;
   default:
      /* Check for image dimension properties */
      if ((icv_property>=MI_ICV_DIM_SIZE) && 
          (icv_property<MI_ICV_DIM_SIZE+MI_MAX_IMGDIMS)) {
         idim = icv_property - MI_ICV_DIM_SIZE;
         icvp->user_dim_size[idim] = value;
      }
      else {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "Tried to set unknown or illegal icv property");
         MI_RETURN_ERROR(MI_ERROR);
      }
      break;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_setint
@INPUT      : icvid        - icv id
              icv_property - property of icv to set
              value        - value to set it to
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets a property of an icv to a given integer value.
              Properties cannot be modified while the icv is attached to a 
              cdf file and variable (see miicv_attach and miicv_detach).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 7, 1992 (Peter Neelin)
@MODIFIED   : January 22, 1993 (P.N.)
                 - modified handling of icv properties
---------------------------------------------------------------------------- */
MNCAPI int miicv_setint(int icvid, int icv_property, int value)
{

   MI_SAVE_ROUTINE_NAME("miicv_setint");

   {MI_CHK_ERR(miicv_setdbl(icvid, icv_property, (double) value))}

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_setlong
@INPUT      : icvid        - icv id
              icv_property - property of icv to set
              value        - value to set it to
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets a property of an icv to a given long integer value.
              Properties cannot be modified while the icv is attached to a 
              cdf file and variable (see miicv_attach and miicv_detach).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_setlong(int icvid, int icv_property, long value)
{

   MI_SAVE_ROUTINE_NAME("miicv_setlong");

   {MI_CHK_ERR(miicv_setdbl(icvid, icv_property, (double) value))}

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_setstr
@INPUT      : icvid        - icv id
              icv_property - property of icv to set
              value        - value to set it to
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets a property of an icv to a given string value. 
              Properties cannot be modified while the icv is attached to a 
              cdf file and variable (see miicv_attach and miicv_detach).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_setstr(int icvid, int icv_property, char *value)
{
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_setstr");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Check that the icv is not attached to a file */
   if (icvp->cdfid != MI_ERROR) {
      MI_LOG_PKG_ERROR2(MI_ERR_ICVATTACHED, 
         "Attempt to modify an attached image conversion variable");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Set the property */
   switch (icv_property) {
   case MI_ICV_SIGN:
      icvp->user_sign   = MI_get_sign_from_string(icvp->user_type, value);
      icvp->user_vmax   = MI_get_default_range(MIvalid_max, icvp->user_type,
                                               icvp->user_sign);
      icvp->user_vmin   = MI_get_default_range(MIvalid_min, icvp->user_type,
                                               icvp->user_sign);
      break;
   case MI_ICV_MAXVAR:
      if (value!=NULL) {
         FREE(icvp->user_maxvar);
         icvp->user_maxvar = strdup(value);
      }
      break;
   case MI_ICV_MINVAR:
      if (value!=NULL) {
         FREE(icvp->user_minvar);
         icvp->user_minvar = strdup(value);
      }
      break;
   case MI_ICV_TYPE:
   case MI_ICV_DO_RANGE:
   case MI_ICV_VALID_MAX:
   case MI_ICV_VALID_MIN:
   case MI_ICV_DO_NORM:
   case MI_ICV_USER_NORM:
   case MI_ICV_IMAGE_MAX:
   case MI_ICV_IMAGE_MIN:
   case MI_ICV_DO_DIM_CONV:
   case MI_ICV_DO_SCALAR:
   case MI_ICV_XDIM_DIR: 
   case MI_ICV_YDIM_DIR:
   case MI_ICV_ZDIM_DIR:
   case MI_ICV_NUM_IMGDIMS:
   case MI_ICV_ADIM_SIZE:
   case MI_ICV_BDIM_SIZE:
   case MI_ICV_KEEP_ASPECT:
      MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                        "Tried to set icv numeric property to a string");
      MI_RETURN_ERROR(MI_ERROR);
      break;
   default:
      /* Check for image dimension properties */
      if ((icv_property>=MI_ICV_DIM_SIZE) && 
          (icv_property<MI_ICV_DIM_SIZE+MI_MAX_IMGDIMS)) {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "Tried to set icv numeric property to a string");
         MI_RETURN_ERROR(MI_ERROR);
      }
      else {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "Tried to set unknown or illegal icv property");
         MI_RETURN_ERROR(MI_ERROR);
      }
      break;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_inqdbl
@INPUT      : icvid        - icv id
              icv_property - icv property to get
@OUTPUT     : value        - value returned
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the value of an icv property
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_inqdbl(int icvid, int icv_property, double *value)
{
   int idim;
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_inqdbl");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Set the property */
   switch (icv_property) {
   case MI_ICV_TYPE:
      *value = icvp->user_type; break;
   case MI_ICV_DO_RANGE:
      *value = icvp->user_do_range; break;
   case MI_ICV_VALID_MAX:
      *value = icvp->user_vmax; break;
   case MI_ICV_VALID_MIN:
      *value = icvp->user_vmin; break;
   case MI_ICV_DO_NORM:
      *value = icvp->user_do_norm; break;
   case MI_ICV_USER_NORM:
      *value = icvp->user_user_norm; break;
   case MI_ICV_IMAGE_MAX:
      *value = icvp->user_imgmax; break;
   case MI_ICV_IMAGE_MIN:
      *value = icvp->user_imgmin; break;
   case MI_ICV_NORM_MAX:
      *value = icvp->derv_imgmax; break;
   case MI_ICV_NORM_MIN:
      *value = icvp->derv_imgmin; break;
   case MI_ICV_DO_FILLVALUE:
      *value = icvp->user_do_fillvalue; break;
   case MI_ICV_FILLVALUE:
      *value = icvp->user_fillvalue; break;
   case MI_ICV_DO_DIM_CONV:
      *value = icvp->user_do_dimconv; break;
   case MI_ICV_DO_SCALAR:
      *value = icvp->user_do_scalar; break;
   case MI_ICV_XDIM_DIR: 
      *value = icvp->user_xdim_dir; break;
   case MI_ICV_YDIM_DIR:
      *value = icvp->user_ydim_dir; break;
   case MI_ICV_ZDIM_DIR:
      *value = icvp->user_zdim_dir; break;
   case MI_ICV_NUM_IMGDIMS:
      *value = icvp->user_num_imgdims; break;
   case MI_ICV_NUM_DIMS:
      *value = icvp->var_ndims;
      if (icvp->var_is_vector && icvp->user_do_scalar) (*value)--;
      break;
   case MI_ICV_CDFID:
      *value = icvp->cdfid; break;
   case MI_ICV_VARID:
      *value = icvp->varid; break;
   case MI_ICV_ADIM_SIZE:
      *value = icvp->user_dim_size[0]; break;
   case MI_ICV_BDIM_SIZE:
      *value = icvp->user_dim_size[1]; break;
   case MI_ICV_ADIM_STEP:
      *value = icvp->derv_dim_step[0]; break;
   case MI_ICV_BDIM_STEP:
      *value = icvp->derv_dim_step[1]; break;
   case MI_ICV_ADIM_START:
      *value = icvp->derv_dim_start[0]; break;
   case MI_ICV_BDIM_START:
      *value = icvp->derv_dim_start[1]; break;
   case MI_ICV_KEEP_ASPECT:
      *value = icvp->user_keep_aspect; break;
   case MI_ICV_SIGN:
   case MI_ICV_MAXVAR:
   case MI_ICV_MINVAR:
      MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                   "Tried to inquire about icv string property as a number");
      MI_RETURN_ERROR(MI_ERROR);
      break;
   default:
      /* Check for image dimension properties */
      if ((icv_property>=MI_ICV_DIM_SIZE) && 
          (icv_property<MI_ICV_DIM_SIZE+MI_MAX_IMGDIMS)) {
         idim = icv_property - MI_ICV_DIM_SIZE;
         *value = icvp->user_dim_size[idim];
      }
      else if ((icv_property>=MI_ICV_DIM_STEP) && 
               (icv_property<MI_ICV_DIM_STEP+MI_MAX_IMGDIMS)) {
         idim = icv_property - MI_ICV_DIM_STEP;
         *value = icvp->derv_dim_step[idim];
      }
      else if ((icv_property>=MI_ICV_DIM_START) && 
               (icv_property<MI_ICV_DIM_START+MI_MAX_IMGDIMS)) {
         idim = icv_property - MI_ICV_DIM_START;
         *value = icvp->derv_dim_start[idim];
      }
      else {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "Tried to inquire about unknown icv property");
         MI_RETURN_ERROR(MI_ERROR);
      }
      break;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_inqint
@INPUT      : icvid        - icv id
              icv_property - icv property to get
@OUTPUT     : value        - value returned
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the value of an icv property
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_inqint(int icvid, int icv_property, int *value)
{
   double dvalue;

   MI_SAVE_ROUTINE_NAME("miicv_inqint");

   {MI_CHK_ERR(miicv_inqdbl(icvid, icv_property, &dvalue))}
   *value = dvalue;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_inqlong
@INPUT      : icvid        - icv id
              icv_property - icv property to get
@OUTPUT     : value        - value returned
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the value of an icv property
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : January 22, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_inqlong(int icvid, int icv_property, long *value)
{
   double dvalue;

   MI_SAVE_ROUTINE_NAME("miicv_inqlong");

   {MI_CHK_ERR(miicv_inqdbl(icvid, icv_property, &dvalue))}
   *value = dvalue;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_inqstr
@INPUT      : icvid        - icv id
              icv_property - icv property to get
@OUTPUT     : value        - value returned. Caller must allocate enough 
                 space for return string.
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the value of an icv property
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_inqstr(int icvid, int icv_property, char *value)
{
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_inqstr");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Set the property */
   switch (icv_property) {
   case MI_ICV_SIGN:
      if (icvp->user_sign==MI_PRIV_SIGNED)
         (void) strcpy(value, MI_SIGNED);
      else if (icvp->user_sign==MI_PRIV_UNSIGNED)
         (void) strcpy(value, MI_UNSIGNED);
      else
         (void) strcpy(value, MI_EMPTY_STRING);
      break;
   case MI_ICV_MAXVAR:
      (void) strcpy(value, icvp->user_maxvar);
      break;
   case MI_ICV_MINVAR:
      (void) strcpy(value, icvp->user_minvar);
      break;
   case MI_ICV_TYPE:
   case MI_ICV_DO_RANGE:
   case MI_ICV_VALID_MAX:
   case MI_ICV_VALID_MIN:
   case MI_ICV_DO_NORM:
   case MI_ICV_USER_NORM:
   case MI_ICV_IMAGE_MAX:
   case MI_ICV_IMAGE_MIN:
   case MI_ICV_NORM_MAX:
   case MI_ICV_NORM_MIN:
   case MI_ICV_DO_DIM_CONV:
   case MI_ICV_DO_SCALAR:
   case MI_ICV_XDIM_DIR: 
   case MI_ICV_YDIM_DIR:
   case MI_ICV_ZDIM_DIR:
   case MI_ICV_NUM_IMGDIMS:
   case MI_ICV_ADIM_SIZE:
   case MI_ICV_BDIM_SIZE:
   case MI_ICV_ADIM_STEP:
   case MI_ICV_BDIM_STEP:
   case MI_ICV_ADIM_START:
   case MI_ICV_BDIM_START:
   case MI_ICV_KEEP_ASPECT:
   case MI_ICV_NUM_DIMS:
   case MI_ICV_CDFID:
   case MI_ICV_VARID:
      MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                  "Tried to inquire about icv numeric property as a string");
      MI_RETURN_ERROR(MI_ERROR);
      break;
   default:
      /* Check for image dimension properties */
      if (((icv_property>=MI_ICV_DIM_SIZE) && 
           (icv_property<MI_ICV_DIM_SIZE+MI_MAX_IMGDIMS)) ||
          ((icv_property>=MI_ICV_DIM_STEP) && 
           (icv_property<MI_ICV_DIM_STEP+MI_MAX_IMGDIMS)) ||
          ((icv_property>=MI_ICV_DIM_START) && 
           (icv_property<MI_ICV_DIM_START+MI_MAX_IMGDIMS))) {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                  "Tried to inquire about icv numeric property as a string");
         MI_RETURN_ERROR(MI_ERROR);
      }
      else {
         MI_LOG_PKG_ERROR2(MI_ERR_BADPROP,
                           "Tried to inquire about unknown icv property");
         MI_RETURN_ERROR(MI_ERROR);
      }
      break;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_ndattach
@INPUT      : icvid - icv id
              cdfid - cdf file id
              varid - cdf variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Attaches an open cdf file and variable to an image conversion
              variable for subsequent access through miicvget and miicvput.
              File must be in data mode. This routine differs from 
              miicv_attach in that no dimension conversions will be made
              on the variable (avoids linking in a significant amount
              of code).
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 9, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_ndattach(int icvid, int cdfid, int varid)
{
   mi_icv_type *icvp;         /* Pointer to icv structure */
   int idim;

   MI_SAVE_ROUTINE_NAME("miicv_ndattach");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* If the icv is attached, then detach it */
   if (icvp->cdfid != MI_ERROR)
      {MI_CHK_ERR(miicv_detach(icvid))}

   /* Inquire about the variable's type, sign and number of dimensions */
   {MI_CHK_ERR(MI_icv_get_type(icvp, cdfid, varid))}

   /* If not doing range calculations, just set derv_firstdim for
      MI_icv_access, otherwise, call routines to calculate range and 
      normalization */
   if (!icvp->user_do_range) {
      icvp->derv_firstdim = -1;
   }
   else {
      /* Get valid range */
      {MI_CHK_ERR(MI_icv_get_vrange(icvp, cdfid, varid))}
          
      /* Get normalization info */
      {MI_CHK_ERR(MI_icv_get_norm(icvp, cdfid, varid))}
   }

   /* Set other fields to defaults */
   icvp->var_is_vector = FALSE;
   icvp->var_vector_size = 1;
   icvp->derv_do_zero = FALSE;
   icvp->derv_do_bufsize_step = FALSE;
   icvp->derv_var_pix_off = NULL;
   icvp->derv_usr_pix_off = NULL;
   for (idim=0; idim<icvp->user_num_imgdims; idim++) {
      icvp->derv_dim_flip[idim] = FALSE;
      icvp->derv_dim_grow[idim] = TRUE;
      icvp->derv_dim_scale[idim] = 1;
      icvp->derv_dim_off[idim] = 0;
      icvp->derv_dim_step[idim] = 0.0;
      icvp->derv_dim_start[idim] = 0.0;
   }

   /* Set the do_scale and do_dimconvert fields of icv structure
      We have to scale only if do_range is TRUE. If ranges don't
      match, or we have to do user normalization, or if we are normalizing
      and MIimagemax or MIimagemin vary over the variable. We don't have
      to scale if input and output are both floating point. */

   icvp->do_scale = 
      (icvp->user_do_range && 
       ((icvp->user_vmax!=icvp->var_vmax) ||
        (icvp->user_vmin!=icvp->var_vmin) ||
        (icvp->user_do_norm && icvp->user_user_norm) ||
        (icvp->user_do_norm && (icvp->derv_firstdim>=0))) );

   if ((icvp->derv_usr_float && icvp->derv_var_float))
      icvp->do_scale = FALSE;

   icvp->do_dimconvert = FALSE;

   /* Set the cdfid and varid fields */
   icvp->cdfid = cdfid;
   icvp->varid = varid;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_type
@INPUT      : icvp  - pointer to icv structure
              cdfid - cdf file id
              varid - variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the type and sign of a variable for miicv_attach.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_get_type(mi_icv_type *icvp, int cdfid, int varid)
{
   int oldncopts;            /* For saving value of ncopts */
   char stringa[MI_MAX_ATTSTR_LEN];
   char *string=stringa;     /* String for sign info */

   MI_SAVE_ROUTINE_NAME("MI_icv_get_type");

   /* Inquire about the variable */
   MI_CHK_ERR(ncvarinq(cdfid, varid, NULL, &(icvp->var_type), 
                       &(icvp->var_ndims), icvp->var_dim, NULL))

   /* Check that the variable type is numeric */
   if (icvp->var_type==NC_CHAR) {
      MI_LOG_PKG_ERROR2(MI_ERR_NONNUMERIC,"Non-numeric datatype");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Try to find out the sign of the variable using MIsigntype. */
   oldncopts = ncopts; ncopts = 0;
   string=miattgetstr(cdfid, varid, MIsigntype, MI_MAX_ATTSTR_LEN, string);
   ncopts = oldncopts;
   icvp->var_sign  = MI_get_sign_from_string(icvp->var_type, string);

   /* Get type lengths */
   icvp->var_typelen = nctypelen(icvp->var_type);
   icvp->user_typelen = nctypelen(icvp->user_type);

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_vrange
@INPUT      : icvp  - pointer to icv structure
              cdfid - cdf file id
              varid - variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the valid range of a variable
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_get_vrange(mi_icv_type *icvp, int cdfid, int varid)
{
   double vrange[2];         /* Valid range buffer */

   MI_SAVE_ROUTINE_NAME("MI_icv_get_vrange");

   if (miget_valid_range(cdfid, varid, vrange) == MI_ERROR) {
      MI_RETURN(MI_ERROR);
   }
   icvp->var_vmin = vrange[0];
   icvp->var_vmax = vrange[1];

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_get_default_range
@INPUT      : what     - MIvalid_min means get default min, MIvalid_min means 
                 get default min
              datatype - type of variable
              sign     - sign of variable
@OUTPUT     : (none)
@RETURNS    : default maximum or minimum for the datatype
@DESCRIPTION: Return the defaults maximum or minimum for a given datatype
              and sign.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE double MI_get_default_range(char *what, nc_type datatype, int sign)
{
   double range[2];

   MI_SAVE_ROUTINE_NAME("MI_get_default_range");

   (void) miget_default_range(datatype, (sign == MI_PRIV_SIGNED), range);

   if (!strcmp(what, MIvalid_max)) {
      MI_RETURN(range[1]);
   }
   else if (!strcmp(what, MIvalid_min)) {
      MI_RETURN(range[0]);
   }
   else {
      ncopts = NC_VERBOSE | NC_FATAL;
      MI_LOG_PKG_ERROR2(-1,"MINC bug - this line should never be printed");
   }

   MI_RETURN(MI_DEFAULT_MIN);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_norm
@INPUT      : icvp  - pointer to icv structure
              cdfid - cdf file id
              varid - variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets the normalization info for a variable
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_get_norm(mi_icv_type *icvp, int cdfid, int varid)
     /* ARGSUSED */
{
   int oldncopts;             /* For saving value of ncopts */
   int vid[2];                /* Variable ids for max and min */
   int ndims;                 /* Number of dimensions for image max and min */
   int dim[MAX_VAR_DIMS];     /* Dimensions */
   int imm;                   /* Counter for looping through max and min */
   double image_range[2];
   int idim, i;

   MI_SAVE_ROUTINE_NAME("MI_icv_get_norm");

   /* Check for floating point or double precision values for user or
      in variable - set flag to not do normalization if needed */
   icvp->derv_var_float = ((icvp->var_type == NC_DOUBLE) ||
                           (icvp->var_type == NC_FLOAT));
   icvp->derv_usr_float = ((icvp->user_type == NC_DOUBLE) ||
                           (icvp->user_type == NC_FLOAT));

   /* Initialize first dimension over which MIimagemax or MIimagemin
      vary - assume that they don't vary at all */
   icvp->derv_firstdim=(-1);

   /* Look for image max, image min variables */
   oldncopts=ncopts; ncopts=0;
   icvp->imgmaxid=ncvarid(cdfid, icvp->user_maxvar);
   icvp->imgminid=ncvarid(cdfid, icvp->user_minvar);
   ncopts = oldncopts;

   /* Check to see if normalization to variable max, min should be done */
   if (!icvp->user_do_norm) {
      icvp->derv_imgmax = MI_DEFAULT_MAX;
      icvp->derv_imgmin = MI_DEFAULT_MIN;
   }
   else {

      /* Get the image min and max, either from the user definition or 
         from the file. */
      if (icvp->user_user_norm) {
         icvp->derv_imgmax = icvp->user_imgmax;
         icvp->derv_imgmin = icvp->user_imgmin;
      }
      else {
         MI_CHK_ERR(miget_image_range(cdfid, image_range))
         icvp->derv_imgmin = image_range[0];
         icvp->derv_imgmax = image_range[1];
      }

      /* Check each of the dimensions of image-min/max variables to see
         which is the fastest varying dimension of the image variable. */
      vid[0]=icvp->imgminid;
      vid[1]=icvp->imgmaxid;
      if ((vid[0] != MI_ERROR) && (vid[1] != MI_ERROR)) {
         for (imm=0; imm < 2; imm++) {
            MI_CHK_ERR(ncvarinq(cdfid, vid[imm], NULL, NULL, 
                                &ndims, dim, NULL))
            for (idim=0; idim<ndims; idim++) {
               for (i=0; i<icvp->var_ndims; i++) {
                  if (icvp->var_dim[i]==dim[idim])
                     icvp->derv_firstdim = longMAX(icvp->derv_firstdim, i);
               }
            }
         }
      }

   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_detach
@INPUT      : icvid - icv id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Detaches the cdf file and variable from the image conversion
              variable, allowing modifications to the icv.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_detach(int icvid)
{
   mi_icv_type *icvp;
   int idim;

   MI_SAVE_ROUTINE_NAME("miicv_detach");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Check that the icv is in fact attached */
   if (icvp->cdfid == MI_ERROR)
      MI_RETURN(MI_NOERROR);

   /* Free the pixel offset arrays */
   if (icvp->derv_var_pix_off != NULL) FREE(icvp->derv_var_pix_off);
   if (icvp->derv_usr_pix_off != NULL) FREE(icvp->derv_usr_pix_off);

   /* Reset values that are read-only (and set when attached) */
   icvp->derv_imgmax = MI_DEFAULT_MAX;
   icvp->derv_imgmin = MI_DEFAULT_MIN;
   for (idim=0; idim<MI_MAX_IMGDIMS; idim++) {
      icvp->derv_dim_step[idim] = 0.0;
      icvp->derv_dim_start[idim] = 0.0;
   }

   /* Set cdfid field to MI_ERROR to indicate that icv is detached */
   icvp->cdfid = MI_ERROR;
   icvp->varid = MI_ERROR;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_get
@INPUT      : icvid  - icv id
              start  - coordinates of start of hyperslab (see ncvarget)
              count  - size of hyperslab (see ncvarget)
@OUTPUT     : values - array of values returned
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets a hyperslab of values from a netcdf variable through
              the image conversion variable (icvid) 
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_get(int icvid, long start[], long count[], void *values)
{
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_get");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);

   /* Get the data */
   MI_CHK_ERR(MI_icv_access(MI_PRIV_GET, icvp, start, count, values))

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_put
@INPUT      : icvid  - icv id
              start  - coordinates of start of hyperslab (see ncvarput)
              count  - size of hyperslab (see ncvarput)
              values - array of values to store
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Stores a hyperslab of values in a netcdf variable through
              the image conversion variable (icvid)
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miicv_put(int icvid, long start[], long count[], void *values)
{
   mi_icv_type *icvp;

   MI_SAVE_ROUTINE_NAME("miicv_put");

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) MI_RETURN_ERROR(MI_ERROR);


   MI_CHK_ERR(MI_icv_access(MI_PRIV_PUT, icvp, start, count, values))

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_access
@INPUT      : operation - MI_PRIV_GET or MI_PRIV_PUT
              icvid     - icv id
              start     - coordinates of start of hyperslab (see ncvarput)
              count     - size of hyperslab (see ncvarput)
              values    - array of values to put
@OUTPUT     : values    - array of values to get
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Does the work of getting or putting values from an icv.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 11, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_access(int operation, mi_icv_type *icvp, long start[], 
                          long count[], void *values)
{
   int *bufsize_step;                /* Pointer to array giving increments
                                        for allocating variable buffer
                                        (NULL if we don't care) */
   long chunk_count[MAX_VAR_DIMS];   /* Number of elements to get for chunk */
   long chunk_start[MAX_VAR_DIMS];   /* Starting index for getting a chunk */
   long chunk_size;                  /* Size of chunk in bytes */
   void *chunk_values;               /* Pointer to next chunk to get */
   long var_start[MAX_VAR_DIMS];     /* Coordinates of first var element */
   long var_count[MAX_VAR_DIMS];     /* Edge lengths in variable */
   long var_end[MAX_VAR_DIMS];       /* Coordinates of last var element */
   int firstdim;
   int idim, ndims;

   MI_SAVE_ROUTINE_NAME("MI_icv_access");

   /* Check that icv is attached to a variable */
   if (icvp->cdfid == MI_ERROR) {
      MI_LOG_PKG_ERROR2(MI_ERR_ICVNOTATTACHED, 
                        "ICV is not attached to an image variable");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Zero the user's buffer if needed */
   if ((operation == MI_PRIV_GET) && (icvp->derv_do_zero))
      {MI_CHK_ERR(MI_icv_zero_buffer(icvp, count, values))}

   /* Translate icv coordinates to variable coordinates */
   {MI_CHK_ERR(MI_icv_coords_tovar(icvp, start, count, var_start, var_count))}

   /* Save icv coordinates for future reference (for dimension conversion
      routines) */
   ndims = icvp->var_ndims;
   if (icvp->var_is_vector && icvp->user_do_scalar)
      ndims--;
   for (idim=0; idim < ndims; idim++) {
      icvp->derv_icv_start[idim] = start[idim];
      icvp->derv_icv_count[idim] = count[idim];
   }

   /* Do we care about getting variable in convenient increments ? 
      Only if we are getting data and the icv structure wants it */
   if ((operation==MI_PRIV_GET) && (icvp->derv_do_bufsize_step))
      bufsize_step = icvp->derv_bufsize_step;
   else
      bufsize_step = NULL;

   /* Set up variables for looping through variable. The biggest chunk that
      we can get in one call is determined by the subscripts of MIimagemax
      and MIimagemin. These must be constant over the chunk that we get if
      we are doing normalization. */
   for (idim=0; idim<icvp->var_ndims; idim++) {
      chunk_start[idim] = var_start[idim];
      var_end[idim]=var_start[idim]+var_count[idim];
   }
   (void) miset_coords(icvp->var_ndims, 1L, chunk_count);
   /* Get size of chunk in user's buffer. Dimension conversion routines
      don't need the buffer pointer incremented - they do it themselves */
   if (!icvp->do_dimconvert)
      chunk_size = nctypelen(icvp->user_type);
   else
      chunk_size = 0;
   for (idim=longMAX(icvp->derv_firstdim+1,0); idim < icvp->var_ndims; idim++) {
      chunk_count[idim]=var_count[idim];
      chunk_size *= chunk_count[idim];
   }
   firstdim = longMAX(icvp->derv_firstdim, 0);

   /* Loop through variable */
   chunk_values = values;
   while (chunk_start[0] < var_end[0]) {

      /* Set the do_fillvalue flag if the user wants it and we are doing
         a get. We must do it inside the loop since the scale factor 
         calculation can change it if the scale is zero. (Fillvalue checking
         is always done if the the scale is zero.) */
      icvp->do_fillvalue = 
         icvp->user_do_fillvalue && (operation == MI_PRIV_GET);
      icvp->fill_valid_min = icvp->var_vmin;
      icvp->fill_valid_max = icvp->var_vmax;

      /* Calculate scale factor */
      if (icvp->do_scale) {
         MI_CHK_ERR(MI_icv_calc_scale(operation, icvp, chunk_start))
      }

      /* Get the values */
      MI_CHK_ERR(MI_varaccess(operation, icvp->cdfid, icvp->varid,
                              chunk_start, chunk_count,
                              icvp->user_type, icvp->user_sign,
                              chunk_values, bufsize_step, icvp))

      /* Increment the start counter */
      chunk_start[firstdim] += chunk_count[firstdim];
      for (idim=firstdim; 
           (idim>0) && (chunk_start[idim]>=var_end[idim]); idim--) {
         chunk_start[idim]=var_start[idim];
         chunk_start[idim-1]++;
      }

      /* Increment the pointer to values */
      chunk_values = (void *) ((char *) chunk_values + (size_t) chunk_size);

   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_zero_buffer
@INPUT      : icvp      - icv structure pointer
              count     - count vector
              values    - pointer to user's buffer
@OUTPUT     : 
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Zeros the user's buffer, with a size given by the vector count.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 9, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_zero_buffer(mi_icv_type *icvp, long count[], void *values)
{
   double zeroval;
   union {
   	unsigned char  uc;
	signed char    sc;
	unsigned short us;
	signed short   ss;
	unsigned int   ui;
	signed int     si;
	float          f;
	double         d;
   } volatile zerobuf = { 0 };
   volatile char *zerostart;
   int zerolen, idim, ndims;
   char *bufptr, *bufend, *zeroptr, *zeroend;
   long buflen;

   MI_SAVE_ROUTINE_NAME("MI_icv_zero_buffer");

   /* Create a zero pixel and get its size */
   zerostart = (volatile char *) (&zerobuf);
   if (icvp->do_scale)
      zeroval = icvp->offset;
   else
      zeroval = 0.0;

   zerolen = icvp->user_typelen;
   
   /* Get the buffer size */
   ndims = icvp->var_ndims;
   if (icvp->var_is_vector && icvp->user_do_scalar)
      ndims--;
   buflen = zerolen;
   for (idim=0; idim<ndims; idim++)
      buflen *= count[idim];

   /* Loop through the buffer, copying the zero pixel */
   bufend = (char *) values + buflen;
   zeroend = (char *) zerostart + zerolen;
   for (bufptr = (char *) values, zeroptr = (char *) zerostart;
        bufptr < bufend; bufptr++, zeroptr++) {
      if (zeroptr >= zeroend)
         zeroptr = (char *) zerostart;
      *bufptr = *zeroptr;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_coords_tovar
@INPUT      : icvp      - icv structure pointer
              icv_start - start vector for icv
              icv_count - count vector for icv
@OUTPUT     : var_start - start vector for variable
              var_count - count vector for variable
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Converts a start and count vector for referencing an icv
              to the corresponding vectors for referencing a NetCDF variable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 1, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_coords_tovar(mi_icv_type *icvp, 
                                long icv_start[], long icv_count[],
                                long var_start[], long var_count[])
{
   int i, j;
   int num_non_img_dims;
   long coord, last_coord, icv_dim_size;

   MI_SAVE_ROUTINE_NAME("MI_icv_coords_tovar");

   /* Do we have to worry about dimension conversions? If not, then
      just copy the vectors and return. */
   if (!icvp->do_dimconvert) {
      for (i=0; i < icvp->var_ndims; i++) {
         var_count[i] = icv_count[i];
         var_start[i] = icv_start[i];
      }
      MI_RETURN(MI_NOERROR);
   }

   /* Get the number of non image dimensions */
   num_non_img_dims=icvp->var_ndims-icvp->user_num_imgdims;
   if (icvp->var_is_vector)
      num_non_img_dims--;

   /* Go through first, non-image dimensions */
   for (i=0; i < num_non_img_dims; i++) {
      var_count[i] = icv_count[i];
      var_start[i] = icv_start[i];
   }

   /* Go through image dimensions */
   for (i=num_non_img_dims, j=icvp->user_num_imgdims-1; 
        i < num_non_img_dims+icvp->user_num_imgdims; i++, j--) {
      /* Check coordinates. */
      icv_dim_size = (icvp->user_dim_size[j] > 0) ?
            icvp->user_dim_size[j] : icvp->var_dim_size[j];
      last_coord = icv_start[i] + icv_count[i] - 1;
      if ((icv_start[i]<0) || (icv_start[i]>=icv_dim_size) ||
          (last_coord<0) || (last_coord>=icv_dim_size) ||
          (icv_count[i]<0)) {
         MI_LOG_PKG_ERROR2(MI_ERR_ICV_INVCOORDS,
                           "Invalid icv coordinates");
         MI_RETURN_ERROR(MI_ERROR);
      }
      /* Remove offset */
      coord = icv_start[i]-icvp->derv_dim_off[j];
      /* Check for growing or shrinking */
      if (icvp->derv_dim_grow[j]) {
         var_count[i] = (icv_count[i]+icvp->derv_dim_scale[j]-1)
            /icvp->derv_dim_scale[j];
         coord /= icvp->derv_dim_scale[j];
      }
      else {
         var_count[i] = icv_count[i]*icvp->derv_dim_scale[j];
         coord *= icvp->derv_dim_scale[j];
      }
      /* Check for flipping */
      if (icvp->derv_dim_flip[j])
         coord = icvp->var_dim_size[j] - coord -
            ((icv_count!=NULL) ? var_count[i] : 0L);
      var_start[i] = coord;
      /* Check for indices out of variable bounds (but in icv bounds) */
      last_coord = var_start[i] + var_count[i];
      if ((var_start[i]<0) || (last_coord>=icvp->var_dim_size[j])) {
         if (var_start[i]<0) var_start[i] = 0;
         if (last_coord>=icvp->var_dim_size[j]) 
            last_coord = icvp->var_dim_size[j] - 1;
         var_count[i] = last_coord - var_start[i] + 1;
      }
   }

   /* Check for vector dimension */
   if (icvp->var_is_vector) {
      if (icvp->user_do_scalar) {
         var_count[icvp->var_ndims-1] = icvp->var_vector_size;
         var_start[icvp->var_ndims-1] = 0;
      }
      else {
         var_count[icvp->var_ndims-1] = icv_count[icvp->var_ndims-1];
         var_start[icvp->var_ndims-1] = icv_start[icvp->var_ndims-1];
      }
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_calc_scale
@INPUT      : operation - MI_PRIV_GET or MI_PRIV_PUT
              icvp      - icv structure pointer
              coords    - coordinates of first value to get or put
@OUTPUT     : icvp      - fields scale and offset set
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Calculates the scale and offset needed for getting or putting
              values, starting at index coords (assumes that scale is constant
              over that range).
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_calc_scale(int operation, mi_icv_type *icvp, long coords[])
{
   long mmcoords[MAX_VAR_DIMS];   /* Coordinates for max/min variable */
   double usr_imgmax, usr_imgmin;
   double var_imgmax, var_imgmin;
   double var_imgmax_true, var_imgmin_true;
   double usr_vmax, usr_vmin;
   double var_vmax, var_vmin;
   double slice_imgmax, slice_imgmin;
   double usr_scale;
   double denom;

   MI_SAVE_ROUTINE_NAME("MI_icv_calc_scale");

   /* Set variable valid range */
   var_vmax = icvp->var_vmax;
   var_vmin = icvp->var_vmin;

   /* Set image max/min for user and variable values depending on whether
      normalization should be done or not. Whenever floating-point values
      are involved, some type of normalization is done. When the icv type
      is floating point, normalization is always done. When the file type
      is floating point and the icv type is integer, slices are normalized
      to the real range of the slice (or chunk being read). */
   if (!icvp->derv_var_float && !icvp->derv_usr_float && !icvp->user_do_norm) {
      usr_imgmax = var_imgmax = MI_DEFAULT_MAX;
      usr_imgmin = var_imgmin = MI_DEFAULT_MIN;
   }
   else {

      /* Get the real range for the slice or chunk that is being examined */
      slice_imgmax = MI_DEFAULT_MAX;
      slice_imgmin = MI_DEFAULT_MIN;
      if ((!icvp->derv_var_float || !icvp->user_do_norm) &&
          (icvp->imgmaxid!=MI_ERROR) && (icvp->imgminid!=MI_ERROR)) {
         if (mitranslate_coords(icvp->cdfid, icvp->varid, coords, 
                                icvp->imgmaxid, mmcoords) == NULL)
            MI_RETURN_ERROR(MI_ERROR);
         {MI_CHK_ERR(mivarget1(icvp->cdfid, icvp->imgmaxid, mmcoords,
                               NC_DOUBLE, NULL, &slice_imgmax))}
         if (mitranslate_coords(icvp->cdfid, icvp->varid, coords, 
                                icvp->imgminid, mmcoords) == NULL)
            MI_RETURN_ERROR(MI_ERROR);
         {MI_CHK_ERR(mivarget1(icvp->cdfid, icvp->imgminid, mmcoords,
                               NC_DOUBLE, NULL, &slice_imgmin))}
      }

      /* Get the user real range */
      if (icvp->user_do_norm) {
         usr_imgmax = icvp->derv_imgmax;
         usr_imgmin = icvp->derv_imgmin;
      }
      else {
         usr_imgmax = slice_imgmax;
         usr_imgmin = slice_imgmin;
      }

      /* Get the file real range */
      if (icvp->derv_var_float) {
         var_imgmax = var_vmax;
         var_imgmin = var_vmin;
      }
      else {
         var_imgmax = slice_imgmax;
         var_imgmin = slice_imgmin;
      }
   }

   /* Prevent scaling between file floats and real value */
   if (icvp->derv_var_float) {
      var_imgmax = var_vmax;
      var_imgmin = var_vmin;
   }

   /* Get user valid range */
   if (icvp->derv_usr_float) {
      usr_vmax = usr_imgmax;
      usr_vmin = usr_imgmin;
   }
   else {
      usr_vmax = icvp->user_vmax;
      usr_vmin = icvp->user_vmin;
   }

   /* Save real var_imgmin/max for fillvalue checking later */
   var_imgmax_true = var_imgmax;
   var_imgmin_true = var_imgmin;

   /* Even though we have already carefully set the vmax/min and imgmax/min
      values to handle the floating point case, we can still have problems
      with the scale calculations (rounding errors) if full range max/min 
      are used (-FLT_MAX to FLT_MAX). To avoid this, we just force the
      values to 0 and 1 which will give the correct scale. That is why 
      we save the true values above. */

   if (icvp->derv_usr_float) {
      usr_imgmax = usr_vmax = MI_DEFAULT_MAX;
      usr_imgmin = usr_vmin = MI_DEFAULT_MIN;
   }
   if (icvp->derv_var_float) {
      var_imgmax = var_vmax = MI_DEFAULT_MAX;
      var_imgmin = var_vmin = MI_DEFAULT_MIN;
   }

   /* Calculate scale and offset for MI_PRIV_GET */

   /* Scale */
   denom = usr_imgmax - usr_imgmin;
   if (denom!=0.0)
      usr_scale=(usr_vmax - usr_vmin) / denom;
   else
      usr_scale=0.0;
   denom = var_vmax - var_vmin;
   if (denom!=0.0)
      icvp->scale = usr_scale * (var_imgmax - var_imgmin) / denom;
   else
      icvp->scale = 0.0;

   /* Offset */
   icvp->offset = usr_vmin - icvp->scale * var_vmin
                + usr_scale * (var_imgmin - usr_imgmin);

   /* If we want a MI_PRIV_PUT, invert */
   if (operation==MI_PRIV_PUT) {
      if (icvp->scale!=0.0) {
         icvp->offset = (-icvp->offset) / icvp->scale;
         icvp->scale  = 1.0/icvp->scale;
      }
      else {
         icvp->offset = var_vmin;
         icvp->scale  = 0.0;
      }
   }

   /* Do fill value checking if scale is zero */
   if (icvp->scale == 0.0) {

      /* Check for floating point on both sides of conversion. We should
         not be doing scaling in this case, but we will check to be safe. */
      if (icvp->derv_var_float && icvp->derv_usr_float) {
         icvp->do_scale = FALSE;
         icvp->do_fillvalue = FALSE;
      }

      else {      /* Not pure floating point */

         icvp->do_fillvalue = TRUE;

         /* For output, set the range properly depending on whether the user
            type is floating point or not */
         if (operation == MI_PRIV_PUT) {
            if (icvp->derv_usr_float) {
               icvp->fill_valid_min = var_imgmin_true;
               icvp->fill_valid_max = var_imgmax_true;
            }
            else if (usr_scale != 0.0) {
               icvp->fill_valid_min = 
                  usr_vmin + (var_imgmin_true - usr_imgmin) / usr_scale;
               icvp->fill_valid_max = 
                  usr_vmin + (var_imgmax_true - usr_imgmin) / usr_scale;
            }
            else {
               icvp->fill_valid_min = usr_vmin;
               icvp->fill_valid_max = usr_vmax;
            }
         }        /* If output operation */

      }        /* If not pure floating point */

   }       /* If scale == 0.0 */

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_chkid
@INPUT      : icvid  - icv id
@OUTPUT     : (none)
@RETURNS    : Pointer to icv structure if it exists, otherwise NULL.
@DESCRIPTION: Checks that an icv id is valid and returns a pointer to the
              structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 7, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
SEMIPRIVATE mi_icv_type *MI_icv_chkid(int icvid)
{
   MI_SAVE_ROUTINE_NAME("MI_icv_chkid");

   /* Check icv id */
   if ((icvid<0) || (icvid>=minc_icv_list_nalloc) || 
       (minc_icv_list[icvid]==NULL)) {
      MI_LOG_PKG_ERROR2(MI_ERR_BADICV,"Illegal icv identifier");
      MI_RETURN_ERROR((void *) NULL);
   }

   MI_RETURN(minc_icv_list[icvid]);
}

// Based on minc-1.5.1/libsrc/dim_conversion.c
// which requires the following...
//
/* ----------------------------------------------------------------------------
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
---------------------------------------------------------------------------- */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_dimconv_init
@INPUT      : operation  - MI_PRIV_GET or MI_PRIV_PUT
              icvp       - icv structure pointer
              dcp        - dimconvert structure pointer
              start      - start passed by user
              count      - count passed by user
              values     - pointer to user's data area (for put)
              bufstart   - start of variable buffer
              bufcount   - count of variable buffer
              buffer     - pointer to variable buffer (for get)
@OUTPUT     : 
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets up stuff for MI_icv_dimconvert.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 4, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_dimconv_init(int operation, mi_icv_type *icvp,
                              mi_icv_dimconv_type *dcp,
                              long start[], long count[], void *values,
                              long bufstart[], long bufcount[], void *buffer)
     /* ARGSUSED */
{
   long buffer_len, values_len; /* Buffer lengths, offsets and indices */
   long buffer_off, values_off;
   long buffer_index, values_index;
   int imgdim_high, imgdim_low; /* Range of subscripts of image dimensions */
   int scale, offset;            /* Dimension scale and offset */
   int idim, jdim;
   int fastdim;
   /* Variables for calculating pixel offsets for compress/expand */
   long var_dcount[MI_MAX_IMGDIMS+1], var_dend[MI_MAX_IMGDIMS+1];
   long usr_dcount[MI_MAX_IMGDIMS+1], usr_dend[MI_MAX_IMGDIMS+1];
   long pixcount;
   int var_fd, usr_fd, dshift;
   long ipix;

   MI_SAVE_ROUTINE_NAME("MI_icv_dimconv_init");

   /* Check to see if any compression or expansion needs to be done.
      Work it out for a GET and then swap if a PUT. */
   if (operation==MI_PRIV_GET) {
      dcp->do_compress = icvp->derv_var_compress;
      dcp->do_expand   = icvp->derv_usr_compress;
   }
   else {
      dcp->do_expand   = icvp->derv_var_compress;
      dcp->do_compress = icvp->derv_usr_compress;
   }

   fastdim = icvp->derv_dimconv_fastdim;

   /* Get the indices of high and low image dimensions */
   imgdim_high=icvp->var_ndims-1;
   if (icvp->var_is_vector) imgdim_high--;
   imgdim_low = imgdim_high - icvp->user_num_imgdims + 1;

   /* Get the buffer sizes */
   buffer_len = icvp->var_typelen;
   values_len = icvp->user_typelen;
   for (idim=0; idim < icvp->var_ndims; idim++) {
      buffer_len *= bufcount[idim];
      if (idim<=fastdim)
         values_len *= icvp->derv_icv_count[idim];
   }

   /* Calculate step size for variable and user buffers. This does not
      allow for growing or shrinking pixels. That correction is done below. */
   if (icvp->var_is_vector && icvp->user_do_scalar) {
      dcp->buf_step[fastdim+1] = icvp->var_typelen;
      dcp->buf_step[fastdim] = dcp->buf_step[fastdim+1] * bufcount[fastdim+1];
   }
   else {
      dcp->buf_step[fastdim] = icvp->var_typelen;
   }
   dcp->usr_step[fastdim] = icvp->user_typelen;
   for (idim=fastdim-1; idim>=0; idim--) {
      dcp->buf_step[idim] = dcp->buf_step[idim+1] * bufcount[idim+1];
      dcp->usr_step[idim] = dcp->usr_step[idim+1]
         * icvp->derv_icv_count[idim+1];
   }

   /* Set sign of user steps for flipping, if needed */
   for (idim=imgdim_low; idim <= imgdim_high; idim++) {
      if (icvp->derv_dim_flip[imgdim_high-idim])
         dcp->usr_step[idim] *= (-1);
   }

   /* Get the pointers to the start of buffers and the number of pixels
      in each dimension (count a pixel as one expansion/compression -
      one time through the loop below) */
   buffer_off = 0;
   values_off = 0;
   for (idim=0; idim <= fastdim; idim++) {
      if ((idim < imgdim_low) || (idim > imgdim_high)) {
         dcp->end[idim] = bufcount[idim];
         buffer_index = 0;
         values_index = bufstart[idim] - icvp->derv_icv_start[idim];
      }
      else {
         jdim = imgdim_high - idim;
         scale = icvp->derv_dim_scale[jdim];
         offset = icvp->derv_dim_off[jdim];
         if (icvp->derv_dim_grow[jdim]) {
            dcp->end[idim] = bufcount[idim];
            buffer_index = 0;
            if (!icvp->derv_dim_flip[jdim])
               values_index = bufstart[idim]*scale
                  - icvp->derv_icv_start[idim] + offset;
            else
               values_index = 
                  (icvp->var_dim_size[jdim] - bufstart[idim])*scale
                  - 1 - icvp->derv_icv_start[idim] + offset;
         }
         else {
            dcp->end[idim] = (bufcount[idim] - 1 + bufstart[idim]%scale) 
                                     / scale + 1; 
            buffer_index = -(bufstart[idim] % scale);
            if (!icvp->derv_dim_flip[jdim])
               values_index = bufstart[idim]/scale
                  - icvp->derv_icv_start[idim] + offset;
            else
               values_index = 
                  (icvp->var_dim_size[jdim] - bufstart[idim] - 1)/scale
                  - icvp->derv_icv_start[idim] + offset;
         }
      }
      buffer_off += buffer_index * labs(dcp->buf_step[idim]);
      values_off += values_index * labs(dcp->usr_step[idim]);
   }

   /* Calculate arrays of offsets for compress/expand. */
   if (dcp->do_compress || dcp->do_expand) {
      /* Initialize counters */
      var_fd = icvp->user_num_imgdims-1;
      usr_fd = icvp->user_num_imgdims-1;
      if (icvp->var_is_vector && icvp->user_do_scalar) {
         var_fd++;
         var_dcount[var_fd]=0;
         var_dend[var_fd]=icvp->var_vector_size;
      }
      for (jdim=0; jdim<icvp->user_num_imgdims; jdim++) {
         idim=icvp->user_num_imgdims - jdim - 1;
         var_dcount[idim] = 0;
         usr_dcount[idim] = 0;
         var_dend[idim] = (icvp->derv_dim_grow[jdim] ?
                           1 : longMIN(icvp->var_dim_size[jdim],
                                   icvp->derv_dim_scale[jdim]));
         usr_dend[idim] = (icvp->derv_dim_grow[jdim] ?
                           icvp->derv_dim_scale[jdim] : 1);
      }

      /* Loop through variable buffer pixels */
      pixcount=0;
      dshift = imgdim_low;
      for (ipix=0; ipix<icvp->derv_var_pix_num; ipix++) {
         icvp->derv_var_pix_off[ipix] = pixcount;
         pixcount += dcp->buf_step[var_fd+dshift];
         if ((++var_dcount[var_fd]) >= var_dend[var_fd]) {
            idim=var_fd;
            while ((idim>0) && (var_dcount[idim]>=var_dend[idim])) {
               var_dcount[idim]=0;
               idim--;
               var_dcount[idim]++;
            }
            for (idim=0, pixcount=0; idim<=var_fd; idim++) {
               pixcount += var_dcount[idim] * dcp->buf_step[idim+dshift];
            }
         }
      }

      /* Loop through user buffer pixels */
      pixcount=0;
      dshift = imgdim_low;
      for (ipix=0; ipix<icvp->derv_usr_pix_num; ipix++) {
         icvp->derv_usr_pix_off[ipix] = pixcount;
         pixcount += dcp->usr_step[usr_fd+dshift];
         if ((++usr_dcount[usr_fd]) >= usr_dend[usr_fd]) {
            idim=usr_fd;
            while ((idim>0) && (usr_dcount[idim]>=usr_dend[idim])) {
               usr_dcount[idim]=0;
               idim--;
               usr_dcount[idim]++;
            }
            for (idim=0, pixcount=0; idim<=var_fd; idim++) {
               pixcount += usr_dcount[idim] * dcp->usr_step[idim+dshift];
            }
         }
      }

      /* Correct buffer steps for compress/expand */
      for (idim=imgdim_low; idim <= imgdim_high; idim++) {
         jdim = imgdim_high-idim;
         if (icvp->derv_dim_grow[jdim])
            dcp->usr_step[idim] *= icvp->derv_dim_scale[jdim];
         else
            dcp->buf_step[idim] *= icvp->derv_dim_scale[jdim];
      }

   }           /* if compress/expand */

   /* Set input and output variables */
   if (operation==MI_PRIV_GET) {          /* For a GET */
      dcp->in_pix_num = icvp->derv_var_pix_num;
      dcp->in_pix_off = icvp->derv_var_pix_off;
      dcp->in_pix_first = buffer;
      dcp->in_pix_last = (void *) ((char *)buffer + buffer_len - 1);
      dcp->out_pix_num = icvp->derv_usr_pix_num;
      dcp->out_pix_off = icvp->derv_usr_pix_off;
      dcp->out_pix_first = values;
      dcp->out_pix_last = (void *) ((char *)values + values_len - 1);
      dcp->intype = icvp->var_type;
      dcp->insign = icvp->var_sign;
      dcp->outtype = icvp->user_type;
      dcp->outsign = icvp->user_sign;
      dcp->istep = dcp->buf_step;
      dcp->ostep = dcp->usr_step;
      dcp->istart = (void *) ((char *) buffer + buffer_off);
      dcp->ostart = (void *) ((char *) values + values_off);
   }                   /* if GET */
   else {                                 /* For a PUT */
      dcp->out_pix_num = icvp->derv_var_pix_num;
      dcp->out_pix_off = icvp->derv_var_pix_off;
      dcp->out_pix_first = buffer;
      dcp->out_pix_last = (void *) ((char *)buffer + buffer_len - 1);
      dcp->in_pix_num = icvp->derv_usr_pix_num;
      dcp->in_pix_off = icvp->derv_usr_pix_off;
      dcp->in_pix_first = values;
      dcp->in_pix_last = (void *) ((char *)values + values_len - 1);
      dcp->outtype = icvp->var_type;
      dcp->outsign = icvp->var_sign;
      dcp->intype = icvp->user_type;
      dcp->insign = icvp->user_sign;
      dcp->ostep = dcp->buf_step;
      dcp->istep = dcp->usr_step;
      dcp->ostart = (void *) ((char *) buffer + buffer_off);
      dcp->istart = (void *) ((char *) values + values_off);
   }                   /* if PUT */

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_dimconvert
@INPUT      : operation  - MI_PRIV_GET or MI_PRIV_PUT
              icvp       - icv structure pointer
              start      - start passed by user
              count      - count passed by user
              values     - pointer to user's data area (for put)
              bufstart   - start of variable buffer
              bufcount   - count of variable buffer
              buffer     - pointer to variable buffer (for get)
@OUTPUT     : values     - pointer to user's data area (for get)
              buffer     - pointer to variable buffer (for put)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Converts values and dimensions from an input buffer to the 
              user's buffer. Called by MI_var_action.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 27, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_dimconvert(int operation, mi_icv_type *icvp,
                              long start[], long count[], void *values,
                              long bufstart[], long bufcount[], void *buffer)
{
   mi_icv_dimconv_type dim_conv_struct;
   mi_icv_dimconv_type *dcp;
   double sum0, sum1;           /* Counters for averaging values */
   double dvalue = 0.0;         /* Pixel value */
   long counter[MAX_VAR_DIMS];  /* Dimension loop counter */
   void *ptr, *iptr, *optr;     /* Pointers for looping through fastest dim */
   void *ivecptr[MAX_VAR_DIMS]; /* Pointers to start of each dimension */
   void *ovecptr[MAX_VAR_DIMS];
   long *end;                   /* Pointer to array of dimension ends */
   int fastdim;                 /* Dimension that varies fastest */
   long ipix;                   /* Buffer subscript */
   int idim;                    /* Dimension subscript */
   int notmodified;             /* First dimension not reset */
   int out_of_range;            /* Flag indicating one pixel of sum out of 
                                   range */
   double dmin, dmax, epsilon;  /* Range limits */

   MI_SAVE_ROUTINE_NAME("MI_icv_dimconvert");

   /* Initialize variables */
   dcp = &dim_conv_struct;
   {MI_CHK_ERR(MI_icv_dimconv_init(operation, icvp, dcp, start, count, values,
                                   bufstart, bufcount, buffer))}

   /* Initialize local variables */
   iptr    = dcp->istart;
   optr    = dcp->ostart;
   end     = dcp->end;
   fastdim = icvp->derv_dimconv_fastdim;
   dmax = icvp->fill_valid_max;
   dmin = icvp->fill_valid_min;
   epsilon = (dmax - dmin) * FILLVALUE_EPSILON;
   epsilon = fabs(epsilon);
   dmax += epsilon;
   dmin -= epsilon;

   /* Initialize counters */
   for (idim=0; idim<=fastdim; idim++) {
      counter[idim] = 0;
      ivecptr[idim] = iptr;
      ovecptr[idim] = optr;
   }

   /* Loop through data */

   while (counter[0] < end[0]) {

      /* Compress data by averaging if needed */
      if (!dcp->do_compress) {
         {MI_TO_DOUBLE(dvalue, dcp->intype, dcp->insign, iptr)}
         out_of_range = (icvp->do_fillvalue && 
                         ((dvalue < dmin) || (dvalue > dmax)));
      }
      else {
         sum1 = 0.0;
         sum0 = 0.0;
         out_of_range=FALSE;
         for (ipix=0; ipix<dcp->in_pix_num; ipix++) {
            ptr=(void *) ((char *)iptr + dcp->in_pix_off[ipix]);

            /* Check if we are outside the buffer.
               If we are looking before the buffer, then we need to
               add in the previous result to do averaging properly. If
               we are looking after the buffer, then break. */
            if (ptr<dcp->in_pix_first) {
               /* Get the output value and re-scale it */
               {MI_TO_DOUBLE(dvalue, dcp->outtype, dcp->outsign, optr)}
               if (icvp->do_scale) {
                  dvalue = ((icvp->scale==0.0) ?
                            0.0 : (dvalue - icvp->offset) / icvp->scale);
               }
            }
            else if (ptr>dcp->in_pix_last) {
               continue;
            }
            else {
               {MI_TO_DOUBLE(dvalue, dcp->intype, dcp->insign, ptr)}
            }

            /* Add in the value, checking for range if needed */
            if (icvp->do_fillvalue && ((dvalue < dmin) || (dvalue > dmax))) {
               out_of_range = TRUE;
            }
            else {
               sum1 += dvalue;
               sum0++;
            }
         }         /* Foreach pixel to compress */

         /* Average values */
         if (sum0!=0.0)
            dvalue = sum1/sum0;
         else
            dvalue = 0.0;
      }           /* If compress */

      /* Check for out of range values and scale result */
      if (out_of_range) {
         dvalue = icvp->user_fillvalue;
      }
      else if (icvp->do_scale) {
         dvalue = icvp->scale * dvalue + icvp->offset;
      }

      /* Expand data if needed */
      if (!dcp->do_expand) {
         {MI_FROM_DOUBLE(dvalue, dcp->outtype, dcp->outsign, optr)}
      }
      else {
         for (ipix=0; ipix<dcp->out_pix_num; ipix++) {
            ptr=(void *) ((char *)optr + dcp->out_pix_off[ipix]);

            /* Check if we are outside the buffer. */
            if ((ptr>=dcp->out_pix_first) && (ptr<=dcp->out_pix_last)) {
               {MI_FROM_DOUBLE(dvalue, dcp->outtype, dcp->outsign, ptr)}
            }

         }         /* Foreach pixel to expand */
      }         /* if expand */

      /* Increment the counter and the pointers */
      if ((++counter[fastdim]) < end[fastdim]) {
         optr = (void *) ((char *) optr + dcp->ostep[fastdim]);
         iptr = (void *) ((char *) iptr + dcp->istep[fastdim]);
      }
      else {
         /* If we reach the end of fastdim, then reset the counter and
            increment the next dimension down - keep going as needed.
            The vectors ovecptr and ivecptr give the starting values of optr 
            and iptr for that dimension. */
         idim = fastdim;
         while ((idim>0) && (counter[idim] >= end[idim])) {
            counter[idim] = 0;
            idim--;
            counter[idim]++;
            ovecptr[idim] = (void *)((char *)ovecptr[idim]+dcp->ostep[idim]);
            ivecptr[idim] = (void *)((char *)ivecptr[idim]+dcp->istep[idim]);
         }
         notmodified = idim;

         /* Copy the starting index up the vector */
         for (idim=notmodified+1; idim<=fastdim; idim++) {
            ovecptr[idim]=ovecptr[notmodified];
            ivecptr[idim]=ivecptr[notmodified];
         }

         optr = ovecptr[fastdim];
         iptr = ivecptr[fastdim];
      }      /* if at end of row */

   }      /* while more pixels to process */


   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_dim_conversion
@INPUT      : icvp  - icv pointer
              subsc - array of dimension subscripts for fastest varying 
                 dimensions
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets the variables for dimensions converions
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 8, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_get_dim_conversion(mi_icv_type *icvp, int subsc[])
     /* ARGSUSED */
{
   int idim;

   MI_SAVE_ROUTINE_NAME("MI_icv_get_dim_conversion");

   /* Find out whether we need to compress variable or user buffer */
   icvp->derv_var_compress = (icvp->var_is_vector && icvp->user_do_scalar);
   icvp->derv_usr_compress = FALSE;
   for (idim=0; idim<icvp->user_num_imgdims; idim++) {
      if (icvp->derv_dim_scale[idim]!=1) {
         if (icvp->derv_dim_grow[idim])
            icvp->derv_usr_compress = TRUE;
         else
            icvp->derv_var_compress = TRUE;
      }
   }

   /* Get the fastest varying dimension in user's buffer */
   icvp->derv_dimconv_fastdim = icvp->var_ndims-1;
   if (icvp->var_is_vector && icvp->user_do_scalar)
      icvp->derv_dimconv_fastdim--;

   /* Find out how many pixels to compress/expand for variable and user
      buffers and allocate arrays */
   if (icvp->var_is_vector && icvp->user_do_scalar)
      icvp->derv_var_pix_num=icvp->var_vector_size;
   else
      icvp->derv_var_pix_num=1;
   icvp->derv_usr_pix_num=1;
   for (idim=0; idim<icvp->user_num_imgdims; idim++) {
      if (icvp->derv_dim_grow[idim])
         icvp->derv_usr_pix_num *= icvp->derv_dim_scale[idim];
      else
         icvp->derv_var_pix_num *= longMIN(icvp->var_dim_size[idim], 
                                       icvp->derv_dim_scale[idim]);
   }
   icvp->derv_var_pix_off = MALLOC(icvp->derv_var_pix_num, long);
   icvp->derv_usr_pix_off = MALLOC(icvp->derv_usr_pix_num, long);
   if ((icvp->derv_var_pix_off == NULL) || (icvp->derv_usr_pix_off == NULL)) {
      MI_LOG_SYS_ERROR1("MI_icv_get_dim_conversion");
      MI_RETURN_ERROR(MI_ERROR);
   }

   MI_RETURN(MI_NOERROR);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_dim_bufsize_step
@INPUT      : icvp  - icv pointer
              subsc - array of dimension subscripts for fastest varying 
                 dimensions
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Sets the variables giving variable buffer size
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 3, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_get_dim_bufsize_step(mi_icv_type *icvp, int subsc[])
{
   int idim;

   MI_SAVE_ROUTINE_NAME("MI_get_dim_bufsize_step");

   /* Set default buffer size step */
   for (idim=0; idim < MAX_VAR_DIMS; idim++)
      icvp->derv_bufsize_step[idim]=1;

   /* Check for converting vector to scalar */
   icvp->derv_do_bufsize_step = (icvp->var_is_vector && icvp->user_do_scalar);
   if (icvp->derv_do_bufsize_step)
      icvp->derv_bufsize_step[icvp->var_ndims-1] = icvp->var_vector_size;
      
   /* Check each dimension to see if we need to worry about the variable
      buffer size. This only occurs if we are shrinking the dimension from 
      the variable buffer to the user buffer. */
   for (idim=0; idim < icvp->user_num_imgdims; idim++) {
      if (!icvp->derv_dim_grow[idim])
         icvp->derv_bufsize_step[subsc[idim]]=icvp->derv_dim_scale[idim];
      if (icvp->derv_bufsize_step[subsc[idim]] != 1)
         icvp->derv_do_bufsize_step = TRUE;
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_dim_flip
@INPUT      : icvp  - icv pointer
              cdfid - cdf file id
              dimvid - variable id
              subsc - array of dimension subscripts for fastest varying 
                 dimensions
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Checks for flipping of icv.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 1, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_get_dim_flip(mi_icv_type *icvp, int cdfid, int dimvid[], 
                           int subsc[])
{
   int oldncopts;             /* For saving value of ncopts */
   char dimname[MAX_NC_NAME]; /* Dimensions name */
   int dim_dir;               /* Desired direction for current dimension */
   double dimstep;            /* Dimension step size (and direction) */
   int idim;

   MI_SAVE_ROUTINE_NAME("MI_get_dim_flip");

   /* Loop through fast dimensions */

   for (idim=0; idim < icvp->user_num_imgdims; idim++) {

      /* Get name of the dimension */
      {MI_CHK_ERR(ncdiminq(cdfid, icvp->var_dim[subsc[idim]], dimname, 
                           NULL))}

      /* Should we look for dimension flipping? */
      icvp->derv_dim_flip[idim]=FALSE;
      if (!strcmp(dimname, MIxspace) || 
          !strcmp(dimname, MIxfrequency))
         dim_dir = icvp->user_xdim_dir;
      else if (!strcmp(dimname, MIyspace) || 
               !strcmp(dimname, MIyfrequency))
         dim_dir = icvp->user_ydim_dir;
      else if (!strcmp(dimname, MIzspace) || 
               !strcmp(dimname, MIzfrequency))
         dim_dir = icvp->user_zdim_dir;
      else
         dim_dir = MI_ICV_ANYDIR;

      /* Look for variable corresponding to dimension */
      if (dim_dir != MI_ICV_ANYDIR) {   /* Look for flipping? */

         /* Get the value of the MIstep attribute to determine whether flipping
            is needed. Assume that direction is positive if no step is
            provided. */
         dimstep = 1.0;
         if (dimvid[idim] != MI_ERROR) {   /* if dimension exists */
            oldncopts = ncopts; ncopts = 0;
            (void) miattget1(cdfid, dimvid[idim], MIstep, NC_DOUBLE, &dimstep);
            ncopts = oldncopts;
         }                           /* if dimension exists */
         if (dim_dir == MI_ICV_POSITIVE)
            icvp->derv_dim_flip[idim] = (dimstep<0.0);
         else if (dim_dir == MI_ICV_NEGATIVE)
            icvp->derv_dim_flip[idim] = (dimstep>0.0);
      }                          /* if look for flipping */

   }                           /* for each dimension */

   MI_RETURN(MI_NOERROR);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_dim_scale
@INPUT      : icvp  - icv pointer
              cdfid - cdf file id
              dimvid - dimension variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Checks for scaling of images
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 1, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_get_dim_scale(mi_icv_type *icvp, int cdfid, int dimvid[])
{
   int oldncopts;             /* For saving value of ncopts */
   int min_grow, dim_grow;
   int min_scale, dim_scale;
   double dimstep, dimstart;
   int idim;
   long user_dim_size;

   MI_SAVE_ROUTINE_NAME("MI_get_dim_scale");

   /* Loop through dimensions, calculating scale and looking for smallest 
      one. For each dimension, check to see if we need to shrink or grow the 
      image. (This is get-oriented: grow is TRUE if the variable dimension 
      has to be expanded to fit the user's dimensions). */

   for (idim=0; idim < icvp->user_num_imgdims; idim++) {

      /* Check to see if we user wants resize */
      if (icvp->user_dim_size[idim] <= 0) {
         icvp->derv_dim_grow[idim]=TRUE;
         icvp->derv_dim_scale[idim]=1;
      }
      else {

         /* Check for growing or shrinking */
         icvp->derv_dim_grow[idim] =
            (icvp->var_dim_size[idim] <= icvp->user_dim_size[idim]);

         /* If growing */
         if (icvp->derv_dim_grow[idim]) {
            /* Get scale so that whole image fits in user array */
            icvp->derv_dim_scale[idim] = 
               icvp->user_dim_size[idim] / icvp->var_dim_size[idim];
         }

         /* Otherwise, shrinking. Things are complicated by the fact that
            the external variable must fit completely in the user's array */
         else {

            icvp->derv_dim_scale[idim] = 1 +
               (icvp->var_dim_size[idim] - 1) / icvp->user_dim_size[idim];
         }
      }           /* if user wants resizing */

      /* Check for smallest scale */
      if (idim==0) {
         min_grow = icvp->derv_dim_grow[idim];
         min_scale = icvp->derv_dim_scale[idim];
      }
      else {
         dim_grow  = icvp->derv_dim_grow[idim];
         dim_scale = icvp->derv_dim_scale[idim];
         /* Check for one of three conditions :
               (1) smallest so far is growing, but this dim is shrinking
               (2) both are growing and this dim has smaller scale
               (3) both are shrinking and this dim has larger scale */
         if ((min_grow && !dim_grow) ||
             ((min_grow && dim_grow) &&
              (min_scale > dim_scale)) ||
             ((!min_grow && !dim_grow) &&
              (min_scale < dim_scale))) {
            min_grow = dim_grow;
            min_scale = dim_scale;
         }
      }           /* if not first dim */

   }           /* for each dimension, get scale */

   /* Loop through dimensions, resetting scale if needed, setting offset 
      and pixel step and start */

   for (idim=0; idim < icvp->user_num_imgdims; idim++) {

      /* Check for aspect ratio */
      if (icvp->user_keep_aspect) {
         icvp->derv_dim_grow[idim]  = min_grow;
         icvp->derv_dim_scale[idim] = min_scale;
      }

      /* Get user's buffer size */
      user_dim_size = ((icvp->user_dim_size[idim]<=0) ? 
                       icvp->var_dim_size[idim] :
                       icvp->user_dim_size[idim]);

      /* Set offset of variable into user's image */

      /* If growing */
      if (icvp->derv_dim_grow[idim]) {
         /* Calculate remainder and split it in half */
         icvp->derv_dim_off[idim] = 
            ( user_dim_size -
             icvp->var_dim_size[idim] * icvp->derv_dim_scale[idim] )
                                      / 2;
      }
      /* Otherwise, shrinking. Things are complicated by the fact that
         the external variable must fit completely in the user's array */
      else {
         /* Calculate remainder and split it in half */
         icvp->derv_dim_off[idim] = 
            ( user_dim_size - 1 -
             (icvp->var_dim_size[idim] - 1) 
                              / icvp->derv_dim_scale[idim] ) / 2 ;
      }

      /* Get pixel step and start for variable and calculate for user.
         Look for them in the dimension variable (if MIstep is not
         there, then use defaults step = 1.0, start = 0.0 */
      oldncopts = ncopts; ncopts = 0;
      dimstep = 1.0;
      (void) miattget1(cdfid, dimvid[idim], MIstep, NC_DOUBLE, &dimstep);
      /* Flip dimstep if needed */
      if (icvp->derv_dim_flip[idim]) 
         dimstep *= (-1);
      /* Get step size for user's image */
      icvp->derv_dim_step[idim] = icvp->derv_dim_grow[idim] ?
         dimstep / icvp->derv_dim_scale[idim] :
            dimstep * icvp->derv_dim_scale[idim];
      /* Get start position for user's image - if no MIstart for
         dimension, then assume 0.0 */
      dimstart = 0.0;
      (void) miattget1(cdfid, dimvid[idim], MIstart, NC_DOUBLE, &dimstart);
      /* Flip dimstart if needed */
      if (icvp->derv_dim_flip[idim])
         dimstart -= dimstep * (icvp->var_dim_size[idim]-1);
      /* Calculate start position */
      icvp->derv_dim_start[idim] = dimstart + 
         (icvp->derv_dim_step[idim] - dimstep) / 2.0 -
            icvp->derv_dim_off[idim] * icvp->derv_dim_step[idim];
      ncopts = oldncopts;

   }                 /* for each dimension */

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_icv_get_dim
@INPUT      : icvp  - pointer to icv structure
              cdfid - cdf file id
              varid - variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Gets dimension info for the icv
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 10, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_icv_get_dim(mi_icv_type *icvp, int cdfid, int varid)
     /* ARGSUSED */
{
   int oldncopts;             /* For saving value of ncopts */
   char dimname[MAX_NC_NAME]; /* Dimensions name */
   int idim;                  /* Looping counter for fastest image dims */
   int subsc[MI_MAX_IMGDIMS]; /* Subscripts for fastest image dims */
   int dimvid[MI_MAX_IMGDIMS]; /* Variable ids for dimensions */

   MI_SAVE_ROUTINE_NAME("MI_icv_get_dim");

   /* Check that the variable has at least icvp->user_num_imgdims dimensions */
   if (icvp->var_ndims < icvp->user_num_imgdims) {
      MI_LOG_PKG_ERROR2(MI_ERR_TOOFEWDIMS, 
                       "Variable has too few dimensions to be an image");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Check the first dimensions of the variable */
   MI_CHK_ERR(ncdiminq(cdfid, icvp->var_dim[icvp->var_ndims-1], dimname, 
                       &(icvp->var_vector_size)))
   icvp->var_is_vector = !strcmp(dimname, MIvector_dimension);

   /* Check that the variable has at least icvp->user_num_imgdims+1 
      dimensions if it is a vector field */
   if (icvp->var_is_vector && (icvp->var_ndims < icvp->user_num_imgdims+1)) {
      MI_LOG_PKG_ERROR2(MI_ERR_TOOFEWDIMS, 
                        "Variable has too few dimensions to be an image");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Check for dimension flipping and get dimension sizes */

   /* Get subscripts for first icvp->user_num_imgdims dimensions */
   subsc[0] = (icvp->var_is_vector) ? icvp->var_ndims-2 : icvp->var_ndims-1;
   for (idim=1; idim < icvp->user_num_imgdims; idim++)
      subsc[idim]=subsc[idim-1]-1;

   /* Get dimension variable ids */
   for (idim=0; idim < icvp->user_num_imgdims; idim++) {
      {MI_CHK_ERR(ncdiminq(cdfid, icvp->var_dim[subsc[idim]], dimname, 
                       &(icvp->var_dim_size[idim])))};
      oldncopts = ncopts; ncopts = 0;
      dimvid[idim] = ncvarid(cdfid, dimname);
      ncopts = oldncopts;
   }

   /* Check for flipping */
   {MI_CHK_ERR(MI_get_dim_flip(icvp, cdfid, dimvid, subsc))}

   /* Check for scaling of dimension */
   {MI_CHK_ERR(MI_get_dim_scale(icvp, cdfid, dimvid))}

   /* Check for variable buffer size increments */
   {MI_CHK_ERR(MI_get_dim_bufsize_step(icvp, subsc))}

   /* Get information for dimension conversion */
   {MI_CHK_ERR(MI_icv_get_dim_conversion(icvp, subsc))}

   MI_RETURN(MI_NOERROR);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miicv_attach
@INPUT      : icvid - icv id
              cdfid - cdf file id
              varid - cdf variable id
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Attaches an open cdf file and variable to an image conversion
              variable for subsequent access through miicvget and miicvput.
              File must be in data mode.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : September 9, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int miicv_attach(int icvid, int cdfid, int varid)
{
   mi_icv_type *icvp;         /* Pointer to icv structure */
   long size_diff, user_dim_size;
   int idim;

   /* Check icv id */
   if ((icvp=MI_icv_chkid(icvid)) == NULL) return(MI_ERROR);

   /* Call routine to set variables for everything except dimension 
      conversion */
   {MI_CHK_ERR(miicv_ndattach(icvid, cdfid, varid))}

   /* Check to see if we need to worry about dimension conversion */
   if (!icvp->user_do_dimconv) {
      return (MI_NOERROR);
   }

   /* Reset cdfid and varid in icv structure in case something goes wrong 
      in dimension conversion calculations */
   icvp->cdfid = MI_ERROR;
   icvp->varid = MI_ERROR;

   /* Get dimensions info */
   {MI_CHK_ERR(MI_icv_get_dim(icvp, cdfid, varid))}

   /* Set the do_dimconvert field of icv structure
      We do dimension conversion if any dimension needs flipping, scaling
      or offset, or if we have to convert from vector to scalar. */

   icvp->do_dimconvert = (icvp->user_do_scalar && icvp->var_is_vector);
   for (idim=0; idim<icvp->user_num_imgdims; idim++) {
      if (icvp->derv_dim_flip[idim] || (icvp->derv_dim_scale[idim] != 1) ||
          (icvp->derv_dim_off[idim] != 0))
         icvp->do_dimconvert = TRUE;
   }
   icvp->dimconvert_func = MI_icv_dimconvert;

   /* Check if we have to zero user's buffer on GETs */
   icvp->derv_do_zero = FALSE;
   for (idim=0; idim<icvp->user_num_imgdims; idim++) {
      user_dim_size = ((icvp->user_dim_size[idim]<=0) ? 
                       icvp->var_dim_size[idim] :
                       icvp->user_dim_size[idim]);
      if (icvp->derv_dim_grow[idim])
         size_diff = user_dim_size -
                     icvp->var_dim_size[idim] * icvp->derv_dim_scale[idim];
      else
         size_diff = user_dim_size - 1 -
                     (icvp->var_dim_size[idim] - 1)
                        / icvp->derv_dim_scale[idim];
      if ((icvp->derv_dim_off[idim]!=0) || (size_diff!=0))
         icvp->derv_do_zero = TRUE;
   }

   /* Set the cdfid and varid fields */
   icvp->cdfid = cdfid;
   icvp->varid = varid;

   return (MI_NOERROR);
}
