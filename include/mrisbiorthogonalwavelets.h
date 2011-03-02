/**
 * @file  mrisbiorthogonalwavelets.c
 * @brief Routines for biothogonal wavelets transformation.
 *
 * References:
 * P. Schroder and W. Sweldens. Spherical wavelets: Texture processing. 
 * In Rendering Techniques '95. Springer Verlag, 1995.
 * P. Schroder and W. Sweldens. Spherical wavelets: Efficiently representing 
 * functions on the sphere. Computer Graphics Proceedings (SIGGRAPH 95), 
 * pages 161-172, 1995.
 */
/*
 * Original Author: Peng Yu
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef MRISBIORTHOGONALWAVELETS_H
#define MRISBIORTHOGONALWAVELETS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "icosahedron.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "matrix.h"

static MRI_SURFACE  *wavelet_analysis_curv(MRI_SURFACE *mris_out, int order) ;
static MRI_SURFACE  *wavelet_analysis_vec(MRI_SURFACE *mris_out, int order);
static MRI_SURFACE  *wavelet_synthesis_curv(MRI_SURFACE *mris_out, int order) ;
static MRI_SURFACE  *wavelet_synthesis_vec(MRI_SURFACE *mris_out, int order);

#endif
