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
 *    $Date: 2010/03/02 19:38:02 $
 *    $Revision: 1.9 $
 *
 * Copyright (C) 2007-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef GIFTI_LOCAL_H
#define GIFTI_LOCAL_H

#include "gifti_io.h"
#include "mrisurf.h"

MRI_SURFACE * mrisReadGIFTIfile(const char *fname);
int mrisReadScalarGIFTIfile(MRI_SURFACE *mris, const char *fname);
int mrisReadLabelTableGIFTIfile(MRI_SURFACE *mris, const char *fname);
MRI * MRISreadGiftiAsMRI(const char *fname, int read_volume);
int MRISwriteGIFTI(MRIS* mris, const char *fname);
int MRISwriteScalarGIFTI(MRIS* mris,
                         const char *fname,
                         const char *scalar_fname);
int mriWriteGifti(MRI* mri, const char *fname);
int MRISwriteLabelTableGIFTI(MRI_SURFACE *mris, const char *fname);

#endif
