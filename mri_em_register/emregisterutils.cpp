/**
 * @brief linear registration to a gca atlas
 *
 * Various utilities
 */
/*
 * Original Author: Bruce Fischl
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

#include <iostream>

#include "macros.h"
#include "diag.h"
#include "cma.h"
#include "error.h"

#include "emregisterutils.h"

extern int use_variance ;

// ===========================================

double local_GCAcomputeLogSampleProbability( GCA *gca,
    GCA_SAMPLE *gcas,
    MRI *mri,
    MATRIX *m_L,
    int nsamples,
    int exvivo,
    double clamp)
{
  static TRANSFORM *transform = NULL ;

  if (!transform)
  {
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  }
  ((LTA *)transform->xform)->xforms[0].m_L = m_L ;

  if (exvivo)
  {
    double gm, wm, fluid ;

    compute_tissue_modes(mri, gca, gcas, transform, nsamples,
                         &wm, &gm, &fluid ) ;
    return( (SQR(gm - wm) +
             SQR(gm-fluid) +
             SQR(fluid - wm) +
             SQR(gm) -
             SQR(fluid))) ;
  }


  double result;

  if (robust)
  {
    // Defined 0 at the top of the file
    result = GCAcomputeNumberOfGoodFittingSamples( gca, gcas, mri,
             transform, nsamples );
  }
  else
  {
    if (use_variance)
      result = GCAcomputeLabelIntensityVariance( gca, gcas, mri,
						 transform, nsamples );
    else
      result = GCAcomputeLogSampleProbability( gca, gcas, mri,
					       transform, nsamples, clamp );
  }


  return( result );
}


// ===========================================


int compute_tissue_modes( MRI *mri_inputs,
                          GCA *gca,
                          GCA_SAMPLE *gcas,
                          TRANSFORM *transform,
                          int nsamples,
                          double *pwm, double *pgm, double *pfluid )
{
  int        x, y, z, i, xp, yp, zp ;
  float      vals[MAX_GCA_INPUTS] ;
  int        countOutside = 0, ngm, nwm, nfluid;
  double     gm, wm, fluid ;

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */

  // store inverse transformation .. forward:input->gca template,
  // inv: gca template->input
  TransformInvert(transform, mri_inputs) ;

  // go through all sample points
  for (ngm = nwm = nfluid = 0, wm = gm = fluid = 0.0, i = 0 ;
       i < nsamples ;
       i++)
  {
    /////////////////// diag code /////////////////////////////
    if (i == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (Gdiag_no == gcas[i].label)
    {
      DiagBreak() ;
    }
    if (i == Gdiag_no ||
        (gcas[i].xp == Gxp && gcas[i].yp == Gyp && gcas[i].zp == Gzp))
    {
      DiagBreak() ;
    }
    ///////////////////////////////////////////////////////////

    // get prior coordinates
    xp = gcas[i].xp ;
    yp = gcas[i].yp ;
    zp = gcas[i].zp ;
    // if it is inside the source voxel
    if (!GCApriorToSourceVoxel(gca, mri_inputs, transform,
                               xp, yp, zp, &x, &y, &z))
    {
      if (x == Gx && y == Gy && z == Gz)
      {
        DiagBreak() ;
      }

      // (x,y,z) is the source voxel position
      gcas[i].x = x ;
      gcas[i].y = y ;
      gcas[i].z = z ;

      // get values from all inputs
      load_vals(mri_inputs, x, y, z, vals, gca->ninputs) ;
      if (FZERO(vals[0]) && gcas[i].label == Gdiag_no)
      {
        DiagBreak() ;
      }

      if (gcas[i].tissue_class == GM_CLASS)
      {
        ngm++ ;
        gm += vals[0] ;
      }
      else if (gcas[i].tissue_class == WM_CLASS)
      {
        nwm++ ;
        wm += vals[0] ;
      }
      else if (gcas[i].tissue_class == FLUID_CLASS)
      {
        nfluid++ ;
        fluid += vals[0] ;
      }

      if (!FZERO(vals[0]))
      {
        DiagBreak() ;
      }
      if (gcas[i].label != Unknown)
      {
        DiagBreak() ;
      }
      if (i == Gdiag_no)
      {
        DiagBreak() ;
      }
    }
    else  // outside the volume
    {
      countOutside++;
    }
  }

  if (nfluid == 0)
  {
    nfluid = 1 ;
  }
  if (ngm == 0)
  {
    ngm = 1 ;
  }
  if (nwm == 0)
  {
    nwm = 1 ;
  }
  wm /= nwm ;
  gm /= ngm ;
  fluid /= nfluid ;
  G_wm_mean = *pwm = wm ;
  G_gm_mean = *pgm = gm ;
  G_fluid_mean = *pfluid = fluid ;

  return(NO_ERROR) ;
}
