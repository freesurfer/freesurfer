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
#include "MRISurfOverlay.h"

MRIS* mrisReadGIFTIfile(const char *fname, MRIS *mris, MRI *outmri=NULL, int *frame=NULL);
MRIS* mrisReadGIFTIdanum(const char *fname, MRIS *mris, int daNum, MRI *outmri=NULL, int *frame=NULL);
MRI* MRISreadGiftiAsMRI(const char *fname, int read_volume);
int MRISwriteGIFTI(MRIS* mris, int intent_code, const char *out_fname, const char *curv_fname);
int mriWriteGifti(MRI* mri, const char *out_fname);

// GIFTI intent related functions extracted from MRISwriteGIFTI()
int MRISwriteGIFTIIntent(MRIS *mris, int intent_code, gifti_image *image, const char *out_fname, const char *curv_fname);
int MRISwriteGIFTIShape(MRIS *mris, gifti_image *image, int intent_code, const char *curv_fname);
int MRISwriteGIFTIStats(MRIS *mris, gifti_image *image, int intent_code);
int MRISwriteGIFTILabel(MRIS *mris, gifti_image *image, int intent_code);
int MRISwriteGIFTISurface(MRIS *mris, gifti_image *image, const char *out_fname);

// function to output multiple overlays
int MRISwriteGIFTICombined(MRIS *mris, MRISurfOverlay *poverlays, const char *out_fname);

// overloaded functions to handle combined GIFTI with multiple intents
int MRISwriteGIFTIIntent(MRIS *mris, const MRI *mri, int stframe, int endframe, gifti_image *image, int intent_code, const char *out_fname, const char *curv_fname, const char *datatype);
int MRISwriteGIFTIShape(MRIS *mris, const MRI *mri, int stframe, int endframe, gifti_image *image, int intent_code, const char *curv_fname, const char *datatype);
int MRISwriteGIFTIStats(MRIS *mris, const MRI *mri, int stframe, int endframe, gifti_image *image, int intent_code, const char *curv_fname, const char *datatype);

// function to return SHAPE and <STATS> intent counts
int getShapeStatIntentCount(const char *fgifti);



#endif
