/**
 * @brief local utilities for GIFTI library
 *
 * This file has some some extra functions for use with the GIFTI
 * utilities. The official utilities reside in gifti_io.c and gifti_xml.c
 * 
 */
/*
 * Original Author: Kevin Teich 
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

#ifndef GIFTI_LOCAL_H
#define GIFTI_LOCAL_H

#include "gifti_io.h"
#include "mrisurf.h"

MRIS* mrisReadGIFTIfile(const char *fname, MRIS *mris);
MRIS* mrisReadGIFTIdanum(const char *fname, MRIS *mris, int daNum);
MRI* MRISreadGiftiAsMRI(const char *fname, int read_volume);
int MRISwriteGIFTI(MRIS* mris, int intent_code, const char *out_fname, const char *curv_fname);
int mriWriteGifti(MRI* mri, const char *out_fname);

#endif
