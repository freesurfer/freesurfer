/*
 * Original Author: (assumed to be) Peter Neelin, modified to compile within freesurfer/utils by Bevin Brett
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


/* ----------------------------- MNI Header -----------------------------------
@NAME       : minc_convenience.c
@DESCRIPTION: File of convenience functions following the minc standard.
@METHOD     : Routines included in this file :
              public :
                 miget_datatype
                 miget_default_range
                 miget_valid_range
                 miset_valid_range
                 miget_image_range
                 mivar_exists
                 miattput_pointer
                 miattget_pointer
                 miadd_child
                 micreate_std_variable
                 micreate_group_variable
              private :
                 MI_create_dim_variable
                 MI_create_dimwidth_variable
                 MI_create_image_variable
                 MI_create_imaxmin_variable
                 MI_verify_maxmin_dims
                 MI_create_root_variable
                 MI_create_simple_variable
                 MI_add_stdgroup
                 MI_is_in_list
@CREATED    : July 27, 1992. (Peter Neelin, Montreal Neurological Institute)
@MODIFIED   : 
 * $Log: minc_convenience.c,v $
 * Revision 6.12.2.1  2004/09/28 20:23:40  bert
 * Minor portability fixes for Windows
 *
 * Revision 6.12  2004/03/24 20:53:48  bert
 * Increase att_length by one in miappend_history() in order to read the entire attribute
 *
 * Revision 6.11  2004/02/02 18:22:46  bert
 * Added miget_version() and miappend_history()
 *
 * Revision 6.10  2001/12/06 14:09:07  neelin
 * Corrected return from mivar_exists to use minc macro MI_RETURN so that
 * ncopts is properly restored.
 *
 * Revision 6.9  2001/11/13 14:15:18  neelin
 * Added functions miget_image_range and mivar_exists
 *
 * Revision 6.8  2001/10/17 14:32:20  neelin
 * Modified miset_valid_range to write out valid_range as double in all
 * cases except float. Unfortunately, writing out values in a type that
 * matched the type of the image variable caused problems with programs
 * linked against old minc libraries.
 *
 * Revision 6.7  2001/09/18 15:44:27  neelin
 * When output type is NC_BYTE, valid_range attribute should have type NC_SHORT.
 *
 * Revision 6.6  2001/08/20 13:19:14  neelin
 * Added function miattget_with_sign to allow the caller to specify the sign
 * of the input attribute since this information is ambiguous. This is
 * necessary for the valid_range attribute which should have the same sign
 * as the image data. Modified miget_valid_range to make use of this function.
 *
 * Revision 6.5  2001/08/16 19:24:11  neelin
 * Fixes to the code handling valid_range values.
 *
 * Revision 6.4  2001/08/16 16:41:32  neelin
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
 * Revision 6.2  2001/04/17 18:40:13  neelin
 * Modifications to work with NetCDF 3.x
 * In particular, changed NC_LONG to NC_INT (and corresponding longs to ints).
 * Changed NC_UNSPECIFIED to NC_NAT.
 * A few fixes to the configure script.
 *
 * Revision 6.1  1999/10/19 14:45:09  neelin
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
 * Revision 3.0  1995/05/15  19:33:12  neelin
 * Release of minc version 0.3
 *
 * Revision 2.1  1995/02/08  19:01:06  neelin
 * Moved private function declarations from minc_routines.h to appropriate file.
 *
 * Revision 2.0  1994/09/28  10:38:02  neelin
 * Release of minc version 0.2
 *
 * Revision 1.18  94/09/28  10:37:12  neelin
 * Pre-release
 * 
 * Revision 1.17  93/08/11  12:06:19  neelin
 * Added RCS logging in source.
 * 
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
//static char rcsid[] = "$Header: /private-cvsroot/minc/libsrc/minc_convenience.c,v 6.12.2.1 2004/09/28 20:23:40 bert Exp $ MINC (MNI)";
#endif

//#include "minc_private.h"
//#include "type_limits.h"
//#include "minc_varlists.h"
#include <limits.h>	//BEVIN 
#include <float.h>	//BEVIN 

//++
// BEVIN definitions to avoid big changes below
//
#include "minc_internals.h"
#include "minc_basic.h"
#include "minc_structures.h"
#include "minc_varlists.h"
#include "minc_netcdf_convenience.h"

void* minc_image_conversion_dummy1 = &longMAX;	// just to avoid an error message
void* minc_image_conversion_dummy2 = &longMIN;	// just to avoid an error message

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
#define MI_RETURN(X)       return (X)
#define MI_RETURN_ERROR(X) return (X)
#define MI_LOG_PKG_ERROR2(CODE,MSG)    { fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, (MSG)); exit(1); }
#define MI_LOG_PKG_ERROR3(CODE,FMT,P1) { fprintf(stderr, "%s:%d ",     __FILE__, __LINE__); \
                                         fprintf(stderr, (FMT), (P1));  \
					 fprintf(stderr, "\n"); \
					 exit(1); \
				       }
#define MI_LOG_SYS_ERROR1(MSG)         { fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, (MSG)); exit(1); }

#define MI_CHK_ERR(expr) { if ((expr)==MI_ERROR) MI_RETURN_ERROR(MI_ERROR); } 

static bool STRINGS_EQUAL(const char* rhs, const char* lhs) {
    return !strcmp(lhs,rhs);
}

//
//--

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_is_in_list
@INPUT      : string    - string for which to look
              list      - list in which to look (must be NULL terminated)
@OUTPUT     : (none)
@RETURNS    : TRUE if found, FALSE if not
@DESCRIPTION: Searches a list of character strings for string and returns
              TRUE if the string is in the list.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 5, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
static int MI_is_in_list(const char *string, const char *list[])
{
   int i;

   MI_SAVE_ROUTINE_NAME("MI_is_in_list");

   for (i=0; list[i] != NULL; i++) {
      if (!strcmp(string, list[i])) MI_RETURN(TRUE);
   }

   MI_RETURN(FALSE);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miget_image_range
@INPUT      : cdfid    - cdf file id
@OUTPUT     : image_range - array containing min and max of image-min/max
                 for the entire file
@RETURNS    : MI_ERROR when an error occurs.
@DESCRIPTION: Gets the image range for a file - that is, the maximum 
              image-max value and the minimum image-min value.
              For float images, ensures that values are cast to float to 
              ensure that they are correctly truncated. If the range cannot 
              be found, then use the default values for int images and 
              valid_range for floating-point images.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : October 19, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miget_image_range(int cdfid, double image_range[])
{
   int oldncopts;                /* For saving value of ncopt */
   int vid[2];                   /* Variable ids for min and max */
   int imgid;                    /* Image variable id */
   nc_type datatype;             /* Type of image variable */
   int is_signed;                /* Indicates if image variable is signed */
   int is_float, no_range_found; /* Flags */
   int imm;                      /* For looping over min and max */
   int ndims;                    /* Number of dimensions of variable */
   int idim;                     /* For looping over dimensions */
   int dim[MAX_VAR_DIMS];        /* Dimension ids of variable */
   long ientry;                  /* For stepping through values */
   long size;                    /* Size of min and max variables */
   long start[MAX_VAR_DIMS];     /* Start of variable */
   long count[MAX_VAR_DIMS];     /* Dimension sizes */
   double *buffer;               /* Pointer to buffer for min/max values */

   MI_SAVE_ROUTINE_NAME("miget_image_range");

   /* Set default values for image_range */
   image_range[0] = MI_DEFAULT_MIN;
   image_range[1] = MI_DEFAULT_MAX;

   /* Get the image-min/max variable ids */
   oldncopts=ncopts; ncopts=0;
   vid[0] = ncvarid(cdfid, MIimagemin);
   vid[1] = ncvarid(cdfid, MIimagemax);
   ncopts = oldncopts;

   /* Get the type information for the image variable */
   if ( ((imgid = ncvarid(cdfid, MIimage)) == MI_ERROR) ||
        (miget_datatype(cdfid, imgid, &datatype, &is_signed) == MI_ERROR) )
      MI_RETURN(MI_ERROR);

   /* No max/min variables, so use valid_range values for floats 
      if it is set and defaults otherwise */
   if ((vid[0] == MI_ERROR) || (vid[1] == MI_ERROR)) {

      /* Check for a floating-point type - if it is, try to get the
         valid_range. If the valid_range was set to full range
         for the type, then that means that the valid range was probably
         not set (and if it was, it was not particularly reasonable). */
      is_float = (datatype == NC_FLOAT || datatype == NC_DOUBLE);
      no_range_found = FALSE;
      if (is_float) {
         if (miget_valid_range(cdfid, imgid, image_range) == MI_ERROR)
            MI_RETURN(MI_ERROR);
         no_range_found = 
            (datatype == NC_FLOAT  && image_range[1] == FLT_MAX) ||
            (datatype == NC_DOUBLE && image_range[1] == DBL_MAX);
      }

      /* If it is not a float, or if the valid range was not set, then use 
         the default. */ 
      if (!is_float || no_range_found) {
         image_range[0] = MI_DEFAULT_MIN;
         image_range[1] = MI_DEFAULT_MAX;
      }

   }

   /* If the variables are there then get the max and min and fastest 
      varying dimension */
   else {

      /* Set initial values */
      image_range[0] = DBL_MAX;
      image_range[1] = -DBL_MAX;

      /* Loop over min and max */
      for (imm=0; imm<2; imm++) {

         /* Get dimension list */
         MI_CHK_ERR(ncvarinq(cdfid, vid[imm], NULL, NULL, 
                             &ndims, dim, NULL))

         /* Loop through dimensions, getting dimension sizes and 
            total min/max variable size */
         size=1;     /* Size of MIimagemin/max variable */
         for (idim=0; idim<ndims; idim++) {
            MI_CHK_ERR(ncdiminq(cdfid, dim[idim], NULL, &(count[idim])))
            size *= count[idim];
         }

         /* Get space */
         if ((buffer=MALLOC(size, double))==NULL) {
            MI_LOG_SYS_ERROR1("miget_image_range");
            MI_RETURN_ERROR(MI_ERROR);
         }

         /* Get values */
         if (mivarget(cdfid, vid[imm], 
                      miset_coords(ndims, 0L, start),
                      count, NC_DOUBLE, NULL, buffer)==MI_ERROR) {
            FREE(buffer);
            MI_RETURN_ERROR(MI_ERROR);
         }

         /* Loop through values, getting max/min */
         for (ientry=0; ientry<size; ientry++) {
            image_range[0] = doubleMIN(image_range[0], buffer[ientry]);
            image_range[1] = doubleMAX(image_range[1], buffer[ientry]);
         }
         FREE(buffer);

      }         /* End for (imm=0; imm<2; imm++) */

   }         /* End if {} else {} no min/max vars */


   /* Handle possible rounding errors in having double image-min/max for
      float image type */
   if (datatype == NC_FLOAT) {
      image_range[0] = (float) image_range[0];
      image_range[1] = (float) image_range[1];
   }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miget_datatype
@INPUT      : cdfid    - cdf file id
              imgid    - image variable id
@OUTPUT     : datatype
              is_signed - TRUE if type is signed
@RETURNS    : MI_ERROR when an error occurs.
@DESCRIPTION: Gets the datatype and sign of the image variable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : August 15, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miget_datatype(int cdfid, int imgid, 
                          nc_type *datatype, int *is_signed)
{
   int old_ncopts;
   int use_default_sign;
   char attstr[MI_MAX_ATTSTR_LEN];

   MI_SAVE_ROUTINE_NAME("miget_datatype");

   /* Get the type information for the variable */
   if (ncvarinq(cdfid, imgid, NULL, datatype, NULL, NULL, NULL) == MI_ERROR)
      MI_RETURN(MI_ERROR);

   /* Save the ncopts value */
   old_ncopts = ncopts;
   ncopts = 0;

   /* Get the sign information */
   if ((miattgetstr(cdfid, imgid, MIsigntype, 
                    MI_MAX_ATTSTR_LEN, attstr) != NULL)) {
      
      use_default_sign = FALSE;
      if (strcmp(attstr, MI_SIGNED) == 0)
         *is_signed = TRUE;
      else if (strcmp(attstr, MI_UNSIGNED) == 0)
         *is_signed = FALSE;
      else
         use_default_sign = TRUE;
   }
   else {
      use_default_sign = TRUE;
   }

   /* Set a default sign if needed */
   if (use_default_sign) {
      if (*datatype == NC_BYTE)
         *is_signed = FALSE;
      else
         *is_signed = TRUE;
   }

   /* Restore ncopts */
   ncopts = old_ncopts;

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miget_default_range
@INPUT      : datatype
              is_signed - TRUE if type is signed
@OUTPUT     : default_range - array containing default range for variable
@RETURNS    : MI_NOERROR
@DESCRIPTION: Gets the default range for a data type.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : August 15, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miget_default_range(nc_type datatype, int is_signed, 
                               double default_range[])
{
   MI_SAVE_ROUTINE_NAME("miget_default_range");

   switch (datatype) {
   case NC_INT:
      default_range[0] = (is_signed) ? INT_MIN : 0;
      default_range[1] = (is_signed) ? INT_MAX : UINT_MAX; 
      break;
   case NC_SHORT:
      default_range[0] = (is_signed) ? SHRT_MIN : 0;
      default_range[1] = (is_signed) ? SHRT_MAX : USHRT_MAX; 
      break;
   case NC_BYTE:
      default_range[0] = (is_signed) ? SCHAR_MIN : 0;
      default_range[1] = (is_signed) ? SCHAR_MAX : UCHAR_MAX; 
      break;
   case NC_FLOAT:
      default_range[0] = -FLT_MAX;
      default_range[1] = FLT_MAX;
      break;
   case NC_DOUBLE:
      default_range[0] = -DBL_MAX;
      default_range[1] = DBL_MAX;
      break;
   default: 
      default_range[0]= MI_DEFAULT_MIN;
      default_range[1]= MI_DEFAULT_MAX;
      break;
   }

   MI_RETURN(MI_NOERROR);
}
/* ----------------------------- MNI Header -----------------------------------
@NAME       : miget_valid_range
@INPUT      : cdfid    - cdf file id
              imgid    - image variable id
@OUTPUT     : valid_range - array containing valid min and max of image
@RETURNS    : MI_ERROR when an error occurs.
@DESCRIPTION: Gets the valid range for an image variable. Ensures that
              the values are cast to the appropriate type to ensure
              that they are correctly truncated. This is particularly
              important for float images. If the range cannot be found,
              then use the default values.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : August 15, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miget_valid_range(int cdfid, int imgid, double valid_range[])
{
   int old_ncopts;
   int status;
   int length;
   nc_type datatype;
   int is_signed;
   char *att_sign;
   double temp;

   MI_SAVE_ROUTINE_NAME("miget_valid_range");

   /* Get the type information for the variable */
   if (miget_datatype(cdfid, imgid, &datatype, &is_signed) == MI_ERROR)
      MI_RETURN(MI_ERROR);

   /* Save the ncopts value */
   old_ncopts = ncopts;
   ncopts = 0;

   /* Get the sign string for the attribute */
   if (is_signed)
      att_sign = MI_SIGNED;
   else
      att_sign = MI_UNSIGNED;

   /* Get valid range */
   status=miattget_with_sign(cdfid, imgid, MIvalid_range, 
                             att_sign, NC_DOUBLE, NULL, 
                             2, valid_range, &length);

   /* If not there, look for the max and min */
   if ((status==MI_ERROR) || (length!=2)) {

      /* Get the default range for the type */
      (void) miget_default_range(datatype, is_signed, valid_range);

      /* Try to read the valid max */
      (void) miattget_with_sign(cdfid, imgid, MIvalid_max, 
                                att_sign, NC_DOUBLE, NULL, 
                                1, &valid_range[1], NULL);

      /* Try to read the valid min */
      (void) miattget_with_sign(cdfid, imgid, MIvalid_min, 
                                att_sign, NC_DOUBLE, NULL,
                                1, &valid_range[0], NULL);

   }

   /* Restore the ncopts value */
   ncopts = old_ncopts;

   /* Make sure that the first element is the minimum */
   if (valid_range[1] < valid_range[0]) {
      temp = valid_range[0];
      valid_range[0] = valid_range[1];
      valid_range[1] = temp;
   }

   /* Cast to the appropriate type and back to make sure that things are
      rounded/truncated properly. This is only really needed for floats */
   switch (datatype) {
   case NC_INT:
   case NC_SHORT:
   case NC_BYTE:
      if (is_signed) {
         valid_range[0] = (int) valid_range[0];
         valid_range[1] = (int) valid_range[1];
      }
      else {
         valid_range[0] = (unsigned int) valid_range[0];
         valid_range[1] = (unsigned int) valid_range[1];
      }
      break;
   case NC_FLOAT:
      valid_range[0] = (float) valid_range[0];
      valid_range[1] = (float) valid_range[1];
      break;
   default:
       fprintf(stderr, "%s:%d bad data type\n", __FILE__, __LINE__);
       exit(1);
   }

   MI_RETURN(MI_NOERROR);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_dim_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
              datatype - type of data to store
              ndims    - number of dimensions - must be 0 or 1
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC dimension variable by calling ncvardef
              and then sets default attributes. The standard variables are 
              identified by name, so an unrecognised name produces an error.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_dim_variable(int cdfid, const char *name, 
                                   nc_type datatype, int ndims)
{
   int dimid;                /* Dimension id (for dimensions variables) */
   int varid;                /* Created variable id */

   MI_SAVE_ROUTINE_NAME("MI_create_dim_variable");

   /* Check for MIvector_dimension - no associated variable */
   if (STRINGS_EQUAL(name, MIvector_dimension)) {
      MI_LOG_PKG_ERROR3(MI_ERR_BAD_STDVAR,
                        "%s is not a standard MINC variable", name);
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Check for ndims being 0 or 1 */
   if (ndims>1) {
      MI_LOG_PKG_ERROR2(MI_ERR_WRONGNDIMS,
                        "Too many dimensions for a dimension variable");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Look for dimension and create the variable */
   MI_CHK_ERR(dimid=ncdimid(cdfid, name))
   MI_CHK_ERR(varid=ncvardef(cdfid, name, datatype, ndims, &dimid))

   /* Standard attributes */
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvarid, MI_STDVAR))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvartype, MI_DIMENSION))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIversion, MI_CURRENT_VERSION))

   /* Add comments for spatial dimensions */ 
   if (STRINGS_EQUAL(name, MIxspace))
      {MI_CHK_ERR(miattputstr(cdfid, varid, MIcomments,
                     "X increases from patient left to right"))}
   else if (STRINGS_EQUAL(name, MIyspace))
      {MI_CHK_ERR(miattputstr(cdfid, varid, MIcomments,
                     "Y increases from patient posterior to anterior"))}
   else if (STRINGS_EQUAL(name, MIzspace))
      {MI_CHK_ERR(miattputstr(cdfid, varid, MIcomments,
                     "Z increases from patient inferior to superior"))}

   /* Dimension attributes */
   if (ndims==0) {
      MI_CHK_ERR(miattputstr(cdfid, varid, MIspacing, MI_REGULAR))
   }
   else {
      MI_CHK_ERR(miattputstr(cdfid, varid, MIspacing, MI_IRREGULAR))
   }
   if (STRINGS_EQUAL(name, MItime))
      MI_CHK_ERR(miattputstr(cdfid, varid, MIalignment, MI_START))
   else
      MI_CHK_ERR(miattputstr(cdfid, varid, MIalignment, MI_CENTRE))

   MI_RETURN(varid);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_dimwidth_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
              datatype - type of data to store
              ndims    - number of dimensions - must be 0 or 1
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC dimension width variable by calling 
              ncvardef and then sets default attributes. The standard 
              variables are identified by name, so an unrecognised name 
              produces an error.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_dimwidth_variable(int cdfid, const char *name, 
                                        nc_type datatype, int ndims)
{
   int dimid;                /* Dimension id (for dimensions variables) */
   int varid;                /* Created variable id */
   char string[MAX_NC_NAME]; /* String for dimension name */
   char *str;

   MI_SAVE_ROUTINE_NAME("MI_create_dimwidth_variable");

   /* Look for dimension name in name (remove width suffix) */
   if ((str=strstr(strcpy(string, name),MI_WIDTH_SUFFIX)) == NULL) {
      MI_LOG_PKG_ERROR2(MI_ERR_BADSUFFIX,"Bad dimension width suffix");
      MI_RETURN_ERROR(MI_ERROR);
   }
   *str='\0';

   /* Check for ndims being 0 or 1 */
   if (ndims>1) {
      MI_LOG_PKG_ERROR2(MI_ERR_WRONGNDIMS,
                        "Too many dimensions for a dimension variable");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Look for the dimension */
   MI_CHK_ERR(dimid=ncdimid(cdfid, string))
   /* Create the variable and set defaults */
   MI_CHK_ERR(varid=ncvardef(cdfid, name, datatype, ndims, &dimid))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvarid, MI_STDVAR))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvartype, MI_DIM_WIDTH))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIversion, MI_CURRENT_VERSION))
   if (ndims==0) {
      MI_CHK_ERR(miattputstr(cdfid, varid, MIspacing, MI_REGULAR))
   }
   else {
      MI_CHK_ERR(miattputstr(cdfid, varid, MIspacing, MI_IRREGULAR))
   }
   MI_CHK_ERR(miattputstr(cdfid, varid, MIfiltertype, MI_SQUARE))

   MI_RETURN(varid);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_verify_maxmin_dims
@INPUT      : cdfid        - cdf file id
              image_ndims  - number of MIimage dimensions
              image_dim    - image dimensions
              maxmin_ndims - number of MIimagemax or MIimagemin dimensions
              maxmin_dim   - max/min dimensions
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if dimensions don't agree
@DESCRIPTION: Verifies that MIimage dimensions and MIimagemax/MIimagemin
              dimensions agree. MIimagemax/MIimagemin cannot vary over the
              two fastest varying (last) dimensions of MIimage - three
              fastest dimensions if MIvector_dimension is the fastest varying
              dimension of MIimage (this maintains the image nature of MIimage
              and its dimensional attributes MIimagemax and MIimagemin).
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 7, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_verify_maxmin_dims(int cdfid,
                                  int image_ndims,  int image_dim[],
                                  int maxmin_ndims, int maxmin_dim[])
{
   char dimname[MAX_NC_NAME];
   int i,j;
   int nbaddims = 2;         /* Number of dimension over which max/min
                                should not vary */

   MI_SAVE_ROUTINE_NAME("MI_verify_maxmin_dims");

   /* Check to see if last dimension is MIvectordimension */
   MI_CHK_ERR(ncdiminq(cdfid, image_dim[image_ndims-1], dimname, NULL))
   if (STRINGS_EQUAL(dimname, MIvector_dimension))
      nbaddims++;

   /* Loop through illegal image dimensions (last nbaddims) checking 
      dimensions against maxmin_dim */
   for (i=longMAX(0,image_ndims-nbaddims); i<image_ndims; i++)
      for (j=0; j<maxmin_ndims; j++)
         if (image_dim[i]==maxmin_dim[j]) {
            MI_LOG_PKG_ERROR2(MI_ERR_MAXMIN_DIMS,
                        "Imagemax/min dimensions vary over image dimensions");
            MI_RETURN_ERROR(MI_ERROR);
         }

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_root_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC root variable by calling ncvardef
              and then sets default attributes.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 6, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_root_variable(int cdfid, const char *name)
{
   int varid;                /* Created variable id */

   MI_SAVE_ROUTINE_NAME("MI_create_root_variable");

   /* Create the variable */
   MI_CHK_ERR(varid=ncvardef(cdfid, name, NC_INT, 0, NULL))

   /* Standard attributes */
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvarid, MI_STDVAR))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvartype, MI_GROUP))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIversion, MI_CURRENT_VERSION))

   /* Add empty parent pointer */
   MI_CHK_ERR(miattputstr(cdfid, varid, MIparent, MI_EMPTY_STRING))

   MI_RETURN(varid);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miadd_child
@INPUT      : cdfid        - cdf file id
              parent_varid - variable id of parent variable
              child_varid  - variable id of child variable
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs.
@DESCRIPTION: Adds the name of child_varid to the children attribute of
              parent_varid and sets the parent attribute of child_varid
              to the name of parent_varid.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miadd_child(int cdfid, int parent_varid, int child_varid)
{
   char *child_list;           /* Pointer to string list of children */
   int child_list_size;        /* Length of child list string */
   int oldncopts;              /* To set and reset ncopts */
   nc_type datatype;           /* Type of attribute */
   int status;                 /* Status of function call */
   char *new_child;            /* String containing name of new child */

   MI_SAVE_ROUTINE_NAME("miadd_child");

   /* Get the size of the child list in the parent. Handle the case where the
      child list does not exist. */
   oldncopts=ncopts; ncopts=0;
   status=ncattinq(cdfid, parent_varid, MIchildren, &datatype, 
                   &child_list_size);
   ncopts=oldncopts;
   if ((status == MI_ERROR) || (datatype != NC_CHAR)) 
      child_list_size=0;

   /* Allocate space for new child list */
   if ((child_list = MALLOC(child_list_size+MAX_NC_NAME+
                            strlen(MI_CHILD_SEPARATOR), char)) == NULL) {
      MI_LOG_SYS_ERROR1("miadd_child");
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Get the old list if needed and check for terminating null character 
      (child_list_size should point to the next spot in child_list, 
      overwriting any null character) */
   if (child_list_size>0) {
      if (ncattget(cdfid, parent_varid, MIchildren, child_list) == MI_ERROR) {
         FREE(child_list);
         MI_RETURN_ERROR(MI_ERROR);
      }
      if (child_list[child_list_size-1] == '\0')
         child_list_size--;

      /* Copy the child list element separator (only if there are other
         elements in the list */
      (void) strcpy(&child_list[child_list_size], MI_CHILD_SEPARATOR);
      child_list_size += strlen(MI_CHILD_SEPARATOR);
   }

   /* Get pointer to name of new child */
   new_child = &child_list[child_list_size];

   /* Add the new child name to the list */
   if (ncvarinq(cdfid, child_varid, new_child, NULL,
                NULL, NULL, NULL) == MI_ERROR) {
      FREE(child_list);
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Check for multiple copies of child */
   if (strstr(child_list, new_child) != new_child) {
      child_list_size -= strlen(MI_CHILD_SEPARATOR);
      child_list[child_list_size] = '\0';
   }

   /* Put the attribute MIchildren */
   if (miattputstr(cdfid, parent_varid, MIchildren, child_list) == MI_ERROR) {
      FREE(child_list);
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Get the parent variable name */
   if (ncvarinq(cdfid, parent_varid, child_list, NULL, NULL, NULL, NULL)
                       == MI_ERROR) {
      FREE(child_list);
      MI_RETURN_ERROR(MI_ERROR);
   }

   /* Put the attribute MIparent */
   if (miattputstr(cdfid, child_varid, MIparent, child_list) == MI_ERROR) {
      FREE(child_list);
      MI_RETURN_ERROR(MI_ERROR);
   }

   FREE(child_list);
   MI_RETURN(MI_NOERROR);

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_add_stdgroup
@INPUT      : cdfid    - cdf file id
              varid    - id of variable
@OUTPUT     : (none)
@RETURNS    : MI_ERROR if an error occurs
@DESCRIPTION: Adds an MI standard variable to the MIchildren list of 
              MIrootvariable and sets some standard attributes. If 
              MIrootvariable does not exist, it is created.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 6, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_add_stdgroup(int cdfid, int varid)
{
   int root_varid;          /* Id of root variable */
   int oldncopts;           /* Old value of ncopts */

   MI_SAVE_ROUTINE_NAME("MI_add_stdgroup");

   /* Check for root variable, and add it if it is not there */
   oldncopts=ncopts; ncopts=0;
   root_varid=ncvarid(cdfid, MIrootvariable);
   ncopts=oldncopts;
   if (root_varid==MI_ERROR) {
      MI_CHK_ERR(root_varid=MI_create_root_variable(cdfid, MIrootvariable))
   }

   /* Add group as child of root */
   MI_CHK_ERR(miadd_child(cdfid, root_varid, varid))

   /* Standard attributes */
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvarid, MI_STDVAR))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvartype, MI_GROUP))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIversion, MI_CURRENT_VERSION))

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : miattput_pointer
@INPUT      : cdfid    - cdf file id
              varid    - variable id
              name     - name of attribute to point to variable
              ptrvarid - variable id of existing variable to which name 
                 should point
@OUTPUT     : (none)
@RETURNS    : MI_ERROR when an error occurs
@DESCRIPTION: Creates an variable attribute which points to another variable
              (generally a multi-dimensional attribute that must be stored
              as a variable). The variable to which the attribute points must
              already exist in the file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miattput_pointer(int cdfid, int varid, const char *name, int ptrvarid)
{
   /* String to hold pointer to variable */
   char pointer_string[MAX_NC_NAME+sizeof(MI_VARATT_POINTER_PREFIX)];
   int index;           /* Index into string */

   MI_SAVE_ROUTINE_NAME("miattput_pointer");

   /* Set the first part of the string */
   index=strlen(strcpy(pointer_string,MI_VARATT_POINTER_PREFIX));

   /* Get the name of the variable to which we should point */
   MI_CHK_ERR(ncvarinq(cdfid, ptrvarid, &(pointer_string[index]), NULL,
                       NULL, NULL, NULL))

   /* Set the attribute of the parent */
   MI_CHK_ERR(miattputstr(cdfid, varid, name, pointer_string))

   /* Get the name of the parent variable */
   MI_CHK_ERR(ncvarinq(cdfid, varid, pointer_string, NULL,
                       NULL, NULL, NULL))

   /* Set the attribute of the variable to which we point */
   MI_CHK_ERR(miattputstr(cdfid, ptrvarid, MIparent, pointer_string))

   /* Set the MIvartype attribute for ptrvarid */
   MI_CHK_ERR(miattputstr(cdfid, ptrvarid, MIvartype, MI_VARATT))

   MI_RETURN(MI_NOERROR);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_image_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
              datatype - type of data to store (see ncvardef)
              ndims    - number of dimensions of variable (see ncvardef)
              dim      - vector of variable dimensions (see ncvardef)
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC image variable by calling ncvardef
              and then sets default attributes.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 6, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_image_variable(int cdfid, const char *name, nc_type datatype,
                                     int ndims, int dim[])
{
   int varid;                /* Created variable id */
   int max_varid;            /* Variable id for dimensional attribute */
   int min_varid;            /* Variable id for dimensional attribute */
   int maxmin_ndims;         /* Number of dimensions in max/min variable */
   int maxmin_dim[MAX_VAR_DIMS];  /* Dimensions of max/min variable */
   int oldncopts;            /* For saving and restoring ncopts */

   MI_SAVE_ROUTINE_NAME("MI_create_image_variable");

   /* Look to see if MIimagemax or MIimagemin exist for dimension checking 
      and pointers */
   oldncopts=ncopts; ncopts=0;
   max_varid=ncvarid(cdfid, MIimagemax);
   min_varid=ncvarid(cdfid, MIimagemin);
   ncopts=oldncopts;
   if (max_varid != MI_ERROR) {
      /* Get MIimagemax dimensions */
      MI_CHK_ERR(ncvarinq(cdfid, max_varid, NULL, NULL, &maxmin_ndims,
                          maxmin_dim, NULL))
      MI_CHK_ERR(MI_verify_maxmin_dims(cdfid, ndims, dim, 
                                       maxmin_ndims, maxmin_dim))
   }
   if (min_varid != MI_ERROR) {
      /* Get MIimagemin dimensions */
      MI_CHK_ERR(ncvarinq(cdfid, min_varid, NULL, NULL, &maxmin_ndims,
                          maxmin_dim, NULL))
      MI_CHK_ERR(MI_verify_maxmin_dims(cdfid, ndims, dim, 
                                       maxmin_ndims, maxmin_dim))
   }

   /* Create the variable */
   MI_CHK_ERR(varid=ncvardef(cdfid, name, datatype, ndims, dim))

   /* Standard attributes */
   MI_CHK_ERR(MI_add_stdgroup(cdfid, varid))

   /* Create pointers to MIimagemax and MIimagemin if they exist */
   if (max_varid!=MI_ERROR) 
      MI_CHK_ERR(miattput_pointer(cdfid, varid, MIimagemax, max_varid))
   if (min_varid!=MI_ERROR) 
      MI_CHK_ERR(miattput_pointer(cdfid, varid, MIimagemin, min_varid))

   MI_RETURN(varid);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_imaxmin_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
              datatype - type of data to store (see ncvardef)
              ndims    - number of dimensions of variable (see ncvardef)
              dim      - vector of variable dimensions (see ncvardef)
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC image maximum or minimum dimensional
              attribute variable by calling ncvardef and then sets default 
              attributes. If MIimage exists, then dimensions are checked
              (MIimagemax and MIimagemin cannot vary over the first two
              dimensions of MIimage (or first three if the first is
              MIvector_dimension)), and a pointer attribute is added to
              MIimage.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 6, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_imaxmin_variable(int cdfid, const char *name, nc_type datatype,
                                       int ndims, int dim[])
{
   int varid;                /* Created variable id */
   int image_varid;          /* Variable id for image */
   int image_ndims;          /* Number of image dimensions */
   int image_dim[MAX_VAR_DIMS]; /* Image dimensions */
   void *fillp;              /* Pointer to fill value */
   int oldncopts;            /* For saving and restoring ncopts */
   int index;
   static char fill_b[]={0,1};
   static short fill_s[]={0,1};
   static int fill_i[]={0,1};
   static float fill_f[]={0.0,1.0};
   static double fill_d[]={0.0,1.0};

   MI_SAVE_ROUTINE_NAME("MI_create_imaxmin_variable");

   /* Look to see if MIimage exists for dimension checking and pointers */
   oldncopts=ncopts; ncopts=0;
   image_varid=ncvarid(cdfid, MIimage);
   ncopts=oldncopts;
   if (image_varid != MI_ERROR) {
      /* Get image dimensions */
      MI_CHK_ERR(ncvarinq(cdfid, image_varid, NULL, NULL, &image_ndims,
                          image_dim, NULL))
      MI_CHK_ERR(MI_verify_maxmin_dims(cdfid, image_ndims, image_dim, 
                                       ndims, dim))
   }

   /* Create the variable */
   MI_CHK_ERR(varid=ncvardef(cdfid, name, datatype, ndims, dim))

   /* Standard attributes */
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvarid, MI_STDVAR))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIvartype, MI_VARATT))
   MI_CHK_ERR(miattputstr(cdfid, varid, MIversion, MI_CURRENT_VERSION))

   /* Attribute for setting default values to something reasonable */
   index = STRINGS_EQUAL(name, MIimagemax) ? 1 : 0;
   fillp = ((datatype==NC_BYTE) ?   (void *) &fill_b[index] :
            (datatype==NC_SHORT) ?  (void *) &fill_s[index] :
            (datatype==NC_INT) ?    (void *) &fill_i[index] :
            (datatype==NC_FLOAT) ?  (void *) &fill_f[index] :
            (datatype==NC_DOUBLE) ? (void *) &fill_d[index] :
                                    (void *) NULL);
   if (fillp != NULL) {
      MI_CHK_ERR(ncattput(cdfid, varid, MI_FillValue, datatype, 1, fillp))
   }

   /* Create pointer from MIimage to max or min if MIimage exists */
   if (image_varid != MI_ERROR) 
      MI_CHK_ERR(miattput_pointer(cdfid, image_varid, name, varid))

   MI_RETURN(varid);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : MI_create_simple_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC variable by calling ncvardef
              and then sets default attributes (only standard ones)
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 6, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
PRIVATE int MI_create_simple_variable(int cdfid, const char *name)
{
   int varid;                /* Created variable id */

   MI_SAVE_ROUTINE_NAME("MI_create_simple_variable");

   /* Create the variable */
   MI_CHK_ERR(varid=ncvardef(cdfid, name, NC_INT, 0, NULL))

   /* Standard attributes */
   MI_CHK_ERR(MI_add_stdgroup(cdfid, varid))

   MI_RETURN(varid);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : micreate_std_variable
@INPUT      : cdfid    - cdf file id
              name     - name of standard variable to create
              datatype - type of data to store (see ncvardef)
              ndims    - number of dimensions of variable (see ncvardef)
              dim      - vector of variable dimensions (see ncvardef)
@OUTPUT     : (none)
@RETURNS    : id of created variable, or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC variable by calling ncvardef
              and then sets default attributes. The standard variables are 
              identified by name, so an unrecognised name produces an error.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int micreate_std_variable(int cdfid, const char *name, nc_type datatype, 
                                 int ndims, int dim[])
{
   int varid;                /* Created variable id */

   MI_SAVE_ROUTINE_NAME("micreate_std_variable");

   /* Check to see if it is a standard dimension */
   if (MI_is_in_list(name, dimvarlist)) {
      MI_CHK_ERR(varid=MI_create_dim_variable(cdfid, name, datatype, ndims))
   }

   /* Check for a dimension width */
   else if (MI_is_in_list(name, dimwidthlist)) {
      MI_CHK_ERR(varid=MI_create_dimwidth_variable(cdfid, name, 
                                                   datatype, ndims))
   }

   /* Check for a standard variable or group */
   else if (MI_is_in_list(name, varlist)) {
      if (STRINGS_EQUAL(name, MIimage))
         MI_CHK_ERR(varid=MI_create_image_variable(cdfid, name, datatype,
                                                   ndims, dim))
      else if ((STRINGS_EQUAL(name, MIimagemax)) ||
               (STRINGS_EQUAL(name, MIimagemin)))
         MI_CHK_ERR(varid=MI_create_imaxmin_variable(cdfid, name, datatype,
                                                     ndims, dim))
      else if (STRINGS_EQUAL(name, MIrootvariable))
         MI_CHK_ERR(varid=MI_create_root_variable(cdfid, name))
      else if (STRINGS_EQUAL(name, MIpatient))
         MI_CHK_ERR(varid=MI_create_simple_variable(cdfid, name))
      else if (STRINGS_EQUAL(name, MIstudy))
         MI_CHK_ERR(varid=MI_create_simple_variable(cdfid, name))
      else if (STRINGS_EQUAL(name, MIacquisition))
         MI_CHK_ERR(varid=MI_create_simple_variable(cdfid, name))
      else {
         MI_LOG_PKG_ERROR3(MI_ERR_BAD_STDVAR, 
                           "%s is not recognised as a standard MINC variable",
                           name);
         MI_RETURN_ERROR(MI_ERROR);
      }
   }

   /* If not in any list, then return an error */
   else {
      MI_LOG_PKG_ERROR3(MI_ERR_BAD_STDVAR,
                        "%s is not a standard MINC variable", name);
      MI_RETURN_ERROR(MI_ERROR);
   }

   MI_RETURN(varid);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miset_valid_range
@INPUT      : cdfid    - cdf file id
              imgid    - image variable id
              valid_range - array containing valid min and max of image
@OUTPUT     : (none)
@RETURNS    : MI_ERROR when an error occurs.
@DESCRIPTION: Sets the valid range for an image variable. Ensures that
              the attribute types match the image variable type.
              This is particularly important for float images because of 
              potential rounding when going from double to float.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : August 15, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miset_valid_range(int cdfid, int imgid, double valid_range[])
{
   nc_type datatype;
   int is_signed;
   int status;
   char *attname;
   float fval[2];

   MI_SAVE_ROUTINE_NAME("miset_valid_range");

   /* Get the type information for the variable */
   if (miget_datatype(cdfid, imgid, &datatype, &is_signed) == MI_ERROR)
      MI_RETURN(MI_ERROR);

   /* Cast to the appropriate type and save. Originally, it was thought
      that casting to the type of the image variable would be a good idea 
      because NetCDF says that it should. Unfortunately, this breaks 
      compatibility in some cases with programs linked with old minc 
      libraries, so we only cast for floats to avoid rounding problems
      (cast from double to float can put values out of range) - for
      everything else double is used. */
   attname = MIvalid_range;
   switch (datatype) {
   case NC_FLOAT:
      fval[0] = valid_range[0];
      fval[1] = valid_range[1];
      status = ncattput(cdfid, imgid, attname, datatype, 2, fval);
      break;
   default:
      status = ncattput(cdfid, imgid, attname, NC_DOUBLE, 2, valid_range);
      break;
   }

   MI_RETURN(status);

}


#ifdef BEVIN_UNSUPPRESS

/* Private functions */
PRIVATE int MI_create_dim_variable(int cdfid, char *name, 
                                   nc_type datatype, int ndims);
PRIVATE int MI_create_dimwidth_variable(int cdfid, char *name, 
                                        nc_type datatype, int ndims);
PRIVATE int MI_create_image_variable(int cdfid, char *name, nc_type datatype,
                                     int ndims, int dim[]);
PRIVATE int MI_create_imaxmin_variable(int cdfid, char *name, nc_type datatype,
                                       int ndims, int dim[]);
PRIVATE int MI_verify_maxmin_dims(int cdfid,
                                  int image_ndims,  int image_dim[],
                                  int maxmin_ndims, int maxmin_dim[]);
PRIVATE int MI_create_root_variable(int cdfid, char *name);
PRIVATE int MI_create_simple_variable(int cdfid, char *name);
PRIVATE int MI_add_stdgroup(int cdfid, int varid);







/* ----------------------------- MNI Header -----------------------------------
@NAME       : mivar_exists
@INPUT      : cdfid    - cdf file id
              varname  - name of variable
@OUTPUT     : (none)
@RETURNS    : TRUE if variable exists, false otherwise
@DESCRIPTION: Checks for the existence of a variable.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines.
@CREATED    : October 22, 2001
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int mivar_exists(int cdfid, char *varname)
{
   int oldncopts;                /* For saving value of ncopt */
   int exists;                   /* Flag */

   MI_SAVE_ROUTINE_NAME("mivar_exists");

   oldncopts = ncopts;
   ncopts = 0;
   exists = (ncvarid(cdfid, varname) != MI_ERROR);
   ncopts = oldncopts;

   MI_RETURN(exists);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miattget_pointer
@INPUT      : cdfid - cdf file id
              varid - variable id
              name  - attribute name that should contain a pointer to
                 a variable
@OUTPUT     : (none)
@RETURNS    : variable id pointed to by name, MI_ERROR if an error occurs.
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF and MINC routines
@CREATED    : August 5, 1992
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int miattget_pointer(int cdfid, int varid, char *name)
{
   /* Character string to hold attribute */
   char pointer_string[MAX_NC_NAME+sizeof(MI_VARATT_POINTER_PREFIX)];
   int index;           /* Index into string */
   char *prefix_string=MI_VARATT_POINTER_PREFIX;  /* Prefix string */
   int ptrvarid;        /* Id of variable pointed to by name */

   MI_SAVE_ROUTINE_NAME("miattget_pointer");

   /* Get the attribute */
   if (miattgetstr(cdfid, varid, name, sizeof(pointer_string), 
                   pointer_string) == NULL)
      MI_RETURN_ERROR(MI_ERROR);

   /* Check for the prefix */
   for (index=0; prefix_string[index]!='\0'; index++) {
      if (pointer_string[index]!=prefix_string[index]) {
         MI_LOG_PKG_ERROR3(MI_ERR_NOTPOINTER,
                           "Attribute %s is not a pointer to a variable",
                           name);
         MI_RETURN_ERROR(MI_ERROR);
      }
   }

   /* Get the variable id */
   {MI_CHK_ERR((ptrvarid=ncvarid(cdfid, &pointer_string[index])))}

   MI_RETURN(ptrvarid);
}











/* ----------------------------- MNI Header -----------------------------------
@NAME       : micreate_group_variable
@INPUT      : cdfid - cdf file id
              name  - name of standard variable to create
@OUTPUT     : (none)
@RETURNS    : id of created variable or MI_ERROR if an error occurs
@DESCRIPTION: Creates a standard MINC variable whose values and dimensions
              are unimportant by calling ncvardef and then sets default 
              attributes. The standard variables are identified by name, so
              an unrecognised name produces an error.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : August 6, 1992 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int micreate_group_variable(int cdfid, char *name)
{
   int varid;

   MI_SAVE_ROUTINE_NAME("micreate_group_variable");

   MI_CHK_ERR(varid=micreate_std_variable(cdfid, name, NC_INT, 0, NULL))

   MI_RETURN(varid);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miappend_history
@INPUT      : id - cdf file id
              tm_stamp  - timestamp as returned by time_stamp() function.
@OUTPUT     : (none)
@RETURNS    : MI_NOERROR if successfuly
@DESCRIPTION: Appends the string (which should be in the format returned
              by the time_stamp() function) to the global "history" 
              attribute.
@METHOD     : 
@GLOBALS    : 
@CALLS      : NetCDF routines
@CREATED    : January 1, 2004 (Bert Vincent)
@MODIFIED   : 
---------------------------------------------------------------------------- */
MNCAPI int
miappend_history(int fd, const char *tm_stamp)
{
    nc_type att_type;
    int att_length;
    int r;
    char *att_value;

    r = ncattinq(fd, NC_GLOBAL, MIhistory, &att_type, &att_length);
    if (r < 0 || att_type != NC_CHAR) {
        att_length = 0;
    }

    /* For some reason, miattgetstr() needs to receive a value
     * one larger than the value returned by ncattinq() in order
     * to account for the terminating null character.
     */
    att_length++;

    /* Allocate enough bytes for the existing attribute, the string which 
     * will be appended, a terminating null character, and a possible
     * additional newline.
     */
    att_value = malloc(att_length + strlen(tm_stamp) + 1);
    if (att_value == NULL) {
        return (MI_ERROR);
    }
    if (miattgetstr(fd, NC_GLOBAL, MIhistory, att_length, att_value) == NULL) {
        return (MI_ERROR);
    }

    if (att_value[att_length-1] == '\0') {
        att_length--;
    }

    if (att_value[att_length-1] != '\n') {
        att_value[att_length] = '\n';
        att_length++;
    }

    /* Append the new history.
     */
    strcpy(att_value + att_length, tm_stamp);

    r = miattputstr(fd, NC_GLOBAL, MIhistory, att_value);

    free(att_value);

    return (r);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : miget_version
@INPUT      : (none)
@OUTPUT     : const char *
@RETURNS    : A string describing the MINC library version.
@DESCRIPTION: Just returns a fixed string.
              
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : December 8 2003
@MODIFIED   : 
---------------------------------------------------------------------------- */
const char * miget_version(void)
{
    return (VERSION);
}


#endif 	// BEVIN_UNSUPPRESS
