/**
 * @file  gifti_local.h
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti_io.c and gifti_xml.c
 * 
 */
/*
 * Original Author: Kevin Teich 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/03/07 00:11:57 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2007-2008,
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
#include "mrisurf.h"

MRI_SURFACE * mrisReadGIFTIfile(char *fname);
int MRISwriteGIFTI(MRIS* mris, char *fname);
int MRISwriteScalarGIFTI(MRIS* mris, char *fname, char *scalar_fname);

#endif
