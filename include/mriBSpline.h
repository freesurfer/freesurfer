/**
 * @brief Compute B-Spline Coefficients and Interpolate MRI
 *
 * This Code is based on:
 *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *		IEEE Transactions on Medical Imaging,
 *		vol. 19, no. 7, pp. 739-758, July 2000.
 *
 * and code (for 2D images) obtained from Philippe Thevenaz at
 *   http://bigwww.epfl.ch/thevenaz/interpolation/
 *
 */
/*
 * Original Author: Martin Reuter
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

#ifndef MRIBSPLINE_H
#define MRIBSPLINE_H

#include "mri.h"
#include "transform.h"

/** BSpline structure has degree, coefficient image (float), source type and
    if the source had negative values */
typedef struct
{
   int degree;
   MRI * coeff;
   int srctype;
   int srcneg;
}
MRI_BSPLINE;

/**
  Allocate BSpline (mainly coefficients).
  Usually not called from outside.
  Has to be same size as the MRI for which this will be computed. */
MRI_BSPLINE* MRIallocBSpline(int width, int height, int depth, int nframes);

/** Free BSpline */
int MRIfreeBSpline(MRI_BSPLINE **pbspline) ;

/** Compute B-spline coefficients from image (step 1)
    bspline needs to be same dim as mri_src or NULL. */
MRI_BSPLINE* MRItoBSpline (const MRI *mri_src, MRI_BSPLINE *bspline, int degree);

/** Based on pre-computed B-spline coefficients interpolate image */
int MRIsampleBSpline(const MRI_BSPLINE * bspline, double x, double y, double z, int frame, double *pval);

/** Based on pre-computed B-spline coefficients interpolate image sequence */
int MRIsampleSeqBSpline(const MRI_BSPLINE *bspline, double x, double y, double z,
                        float *valvect, int firstframe, int lastframe);

/** Based on pre-computed B-spline coefficients interpolate image using voxel map MATRIX*/
MRI *MRIlinearTransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, MATRIX *mA);

/** Based on pre-computed B-spline coefficients interpolate image using RAS map MATRIX*/
MRI *MRIapplyRASlinearTransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, MATRIX *mA) ;


/** Based on pre-computed B-spline coefficients interpolate image using LTA*/
MRI *LTAtransformBSpline(const MRI_BSPLINE *bspline, MRI *mri_dst, LTA *lta) ;


/** Direct methods for downsample (based on simplified algorithm) */
MRI *MRIdownsample2BSpline(const MRI* mri_src, MRI *mri_dst) ;

/** Direct methods for upsample (based on simplified algorithm) */
MRI *MRIupsample2BSpline(const MRI* mri_src, MRI *mri_dst) ;

#endif
