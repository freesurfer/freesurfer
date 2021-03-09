/**
 * @brief Utilities for retinotopy analysis.
 *
 */
/*
 * Original Author: Doug Greve (and Marty and Anders, for now)
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



#ifndef RETINOTOPY_INC
#define RETINOTOPY_INC

#include "mri.h"
#include "mrisurf.h"

#ifdef X
#undef X
#endif

void RETcompute_fieldsign(MRIS *mris);
void RETcompute_angles(MRIS *mris, double EccenRotAngleRad, double PolarRotAngleRad);
int RETlogMap(MRIS *surf, double k, double a, double xc0, double yc0);
int RETinvLogMapFunc(double xc, double yc, double xc0, double yc0, 
		     double a, double k, double *r, double *theta);
float RETcircsubtract(float a,float b) ;
int RETreverseSign(MRI *mri);

#endif
