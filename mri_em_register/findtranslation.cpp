/**
 * @file  findtranslation.cpp
 * @brief linear registration to a gca atlas
 *
 * Find the best translation for mri_em_register
 */
/*
 * Original Author: Bruce Fischl
 * CUDA version : Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/03/22 15:47:05 $
 *    $Revision: 1.1 $
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

#include "diag.h"
#include "macros.h"
#include "proto.h"

#ifdef FS_CUDA
#define FAST_TRANSLATION 0
#include "devicemanagement.h"
#include "em_register_cuda.h"
#endif


#include "emregisterutils.h"

#include "findtranslation.h"

// ------------------------------------------------------------

double find_optimal_translation( GCA *gca,
                                 GCA_SAMPLE *gcas,
                                 MRI *mri,
                                 int nsamples,
                                 MATRIX *m_L,
                                 float min_trans,
                                 float max_trans,
                                 float trans_steps,
                                 int nreductions ) {
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta,
           log_p, max_log_p, mean_trans ;
  int      i ;

  log_p = 0;
  x_trans = 0;
  y_trans = 0;
  z_trans = 0;


#ifdef FS_CUDA
  CUDA_em_register_Prepare( gca, gcas, mri, nsamples );
#endif // FS_CUDA

  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p = local_GCAcomputeLogSampleProbability
              (gca, gcas, mri, m_L,nsamples, exvivo) ;

  for (i = 0 ; i <= nreductions ; i++)
  {
    delta = (max_trans-min_trans) / trans_steps ;
    if (FZERO(delta))
    {
      return(max_log_p) ;
    }
    if (Gdiag & DIAG_SHOW)
    {
      printf(
        "scanning translations %2.2f->%2.2f (step %2.1f) ",
        min_trans,max_trans, delta) ;
      fflush(stdout) ;
    }

#if defined(FS_CUDA) && FAST_TRANSLATION
    unsigned int nTrans = 1+((max_trans-min_trans)/delta);
    float myMaxLogP, mydx, mydy, mydz;
    CUDA_FindOptimalTranslation( m_L, min_trans, max_trans, nTrans,
                                 &myMaxLogP, &mydx, &mydy, &mydz );
    if( myMaxLogP > max_log_p )
    {
      max_log_p = myMaxLogP;
      x_max = mydx;
      y_max = mydy;
      z_max = mydz;
    }
#else
    for (x_trans = min_trans ; x_trans <= max_trans ; x_trans += delta)
    {
      *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
      for (y_trans = min_trans ; y_trans <= max_trans ; y_trans += delta)
      {
        *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
        for (z_trans= min_trans ;
             z_trans <= max_trans ;
             z_trans += delta)
        {
          *MATRIX_RELT(m_trans, 3, 4) = z_trans ;
          if (nint((x_trans)) == -9 && nint((y_trans)) == -5 &&
              nint((z_trans)) == -7)
          {
            DiagBreak() ;
          }
          // get the transform
          m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
          // calculate the LogSample probability
#ifdef FS_CUDA
          log_p = CUDA_ComputeLogSampleProbability( m_L_tmp );
#else
          log_p =
            local_GCAcomputeLogSampleProbability
            (gca, gcas, mri, m_L_tmp,nsamples, exvivo) ;
#endif

#if 0
          printf( "%s: %8.3f %8.3f %8.3f %12.6f\n",
                  __FUNCTION__,
                  x_trans, y_trans, z_trans,
                  log_p );
#endif

          if (log_p > max_log_p)
          {
            max_log_p = log_p ;
            x_max = x_trans ;
            y_max = y_trans ;
            z_max = z_trans ;
#if 0
            printf("new max p %2.1f found at "
                   "(%2.1f, %2.1f, %2.1f)\n",
                   max_log_p, x_trans, y_trans, z_trans) ;
#endif
          }
        }
      }
    }
#endif

    if( Gdiag & DIAG_SHOW )
    {
      printf(
        "max log p = %12.6f @ (%4.3f, %4.3f, %4.3f)\n",
        max_log_p, x_max, y_max, z_max) ;
    }


    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max ;
    // create a new transform by multiplying the previous one.
    m_L_tmp = MatrixMultiply(m_trans, m_L, m_L_tmp) ;
    MatrixCopy(m_L_tmp, m_L) ;
#ifdef FS_CUDA
    max_log_p = CUDA_ComputeLogSampleProbability( m_L_tmp );
#else
    max_log_p = local_GCAcomputeLogSampleProbability
                (gca, gcas, mri, m_L,nsamples, exvivo) ;
#endif

#if 0
    // Repeat for debugging
    printf(
      "max log p = %12.6f @ (%4.3f, %4.3f, %4.3f)\n",
      max_log_p, x_max, y_max, z_max) ;
#endif


    x_max = y_max = z_max = 0.0 ;
    /* we've translated transform by old maxs */

    mean_trans = (max_trans + min_trans) / 2 ;
    delta = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta ;
    max_trans = mean_trans + delta ;
  }

  MatrixFree(&m_trans) ;

#ifdef FS_CUDA
  CUDA_em_register_Release();
#endif

  return(max_log_p) ;
}
