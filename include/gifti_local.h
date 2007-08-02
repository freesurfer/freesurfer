/**
 * @file  gifti_local.h
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti.c and gifti_xml.c
 * 
 */
/*
 * Original Author: Kevin Teich 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/08/02 21:14:07 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef GIFTI_LOCAL_H
#define GIFTI_LOCAL_H

#include "gifti.h"


/* Create a new DataArray in an image and return a pointer to it. */
DataArray* gifti_alloc_and_add_darray (gifti_image* image);

/* Get or set a single value out of a data array. */
double gifti_get_DA_value_2D (DataArray* da, int row, int col);
void   gifti_set_DA_value_2D (DataArray* da, int row, int col, double value);

#endif
