#ifndef  MINC_NETCDF_CONVENIENCE_H
#define  MINC_NETCDF_CONVENIENCE_H

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

// Based on code that had the following

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

#include "minc_internals.h"

int miopen(const char *path, int mode);
int miclose(int cdfid);

int miattputdbl(int cdfid, int varid, const char *name, double value);
int miattputstr(int cdfid, int varid, const char *name, const char *value);

int mivarput(int cdfid, int varid, long start[], long count[], nc_type datatype, char *sign, void *values);
int mivarput1(int cdfid, int varid, long mindex[], nc_type datatype, char *sign, void *value);

int micopy_all_atts(int incdfid, int invarid, int outcdfid, int outvarid);
int micopy_all_var_defs(int incdfid, int outcdfid, int nexclude, int excluded_vars[]);
int micopy_all_var_values(int incdfid, int outcdfid, int nexclude, int excluded_vars[]);

int micreate(const char *path, int cmode);
   
#endif
