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
 *    $Author: nicks $
 *    $Date: 2008/02/26 01:02:51 $
 *    $Revision: 1.2 $
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

#include "gifti_io.h"

/* Create a new DataArray in an image and return a pointer to it. */
giiDataArray* gifti_alloc_and_add_darray (gifti_image* image);

/* Get or set a single value out of a data array. */
double gifti_get_DA_value_2D (giiDataArray* da, int row, int col);
void gifti_set_DA_value_2D (giiDataArray* da, int row, int col, double value);

#endif
