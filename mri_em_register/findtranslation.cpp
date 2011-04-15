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
 *    $Date: 2011/04/15 17:33:17 $
 *    $Revision: 1.7 $
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

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "diag.h"
#include "macros.h"
#include "proto.h"

#ifdef FS_CUDA
#define FAST_TRANSLATION 0
#include "devicemanagement.h"
#include "em_register_cuda.h"
#endif


//#define OUTPUT_STAGES


#ifdef OUTPUT_STAGES
#ifdef FS_CUDA
const std::string stem( "TransGPU" );
#else
const std::string stem( "TransCPU" );
#endif

const std::string stern( ".output" );


//! See if an affine matrix really inverts a translation
void CheckInverseTranslation( const MATRIX* mat,
                              const float dx,
                              const float dy,
                              const float dz ) {

  bool exactInverse = true;

  // Check 'identity' portion
  for( int i=1; i<4; i++ ) {
    for( int j=1; j<3; j++ ) {
      const float mij = mat->rptr[i][j];
      if( (i==j) ) {
        exactInverse = (exactInverse && (mij==1));
      } else {
        exactInverse = (exactInverse && (mij==0));
      }
    }
  }

  // Check the translation itself
  exactInverse = (exactInverse && ( (-dx) == mat->rptr[1][4] ) );
  exactInverse = (exactInverse && ( (-dy) == mat->rptr[2][4] ) );
  exactInverse = (exactInverse && ( (-dz) == mat->rptr[3][4] ) );
  exactInverse = (exactInverse && ( (1) == mat->rptr[4][4] ) );

  if( !exactInverse ) {
    std::cout << "Inexact inverse" << std::endl;
    MatrixPrint( stdout, mat );
    std::cout << std::setprecision(12) << std::setw(20) << dx
              << std::setprecision(12) << std::setw(20) << dy
              << std::setprecision(12) << std::setw(20) << dz
              << std::endl;
  }
}

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
                                 int nreductions ,
                                 double clamp ) {
  MATRIX   *m_trans, *m_L_tmp ;
  double   x_trans, y_trans, z_trans, x_max, y_max, z_max, delta,
           log_p, max_log_p, mean_trans ;
  int      i ;

  log_p = 0;
  x_trans = 0;
  y_trans = 0;
  z_trans = 0;

#ifdef FS_CUDA
  if( exvivo ) {
    std::cerr << __FUNCTION__
              << ": exvivo not supported!"
              << std::endl;
    exit( EXIT_FAILURE );
  }
  CUDA_em_register_Prepare( gca, gcas, mri, nsamples );
#endif // FS_CUDA

  delta = (max_trans-min_trans) / trans_steps ;
  m_L_tmp = NULL ;
  m_trans = MatrixIdentity(4, NULL) ;
  x_max = y_max = z_max = 0.0 ;
  max_log_p = local_GCAcomputeLogSampleProbability
    (gca, gcas, mri, m_L,nsamples, exvivo, clamp) ;

  for (i = 0 ; i <= nreductions ; i++)
  {
#ifdef OUTPUT_STAGES
    std::stringstream numString;
    numString << std::setw(2)
              << std::setfill('0')
              << i;
    
    std::ofstream outFile( (stem + numString.str() + stern).c_str() );
#endif

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
    unsigned int nTrans = static_cast<unsigned int>( 1+((max_trans-min_trans)/delta) );
    float myMaxLogP, mydx, mydy, mydz;
    CUDA_FindOptimalTranslation( m_L, min_trans, max_trans, nTrans,
                                 clamp,
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
          log_p = CUDA_ComputeLogSampleProbability( m_L_tmp, clamp );
#else
          log_p =
            local_GCAcomputeLogSampleProbability
            (gca, gcas, mri, m_L_tmp,nsamples, exvivo, clamp) ;
#endif


#if 0
          printf( "%s: %8.3f %8.3f %8.3f %12.6f\n",
                  __FUNCTION__,
                  x_trans, y_trans, z_trans,
                  log_p );
#endif

#ifdef OUTPUT_STAGES
          outFile << std::setw(20) << std::setprecision(12) << x_trans << ",";
          outFile << std::setw(20) << std::setprecision(12) << y_trans << ",";
          outFile << std::setw(20) << std::setprecision(12) << z_trans << ",";
          outFile << std::setw(20) << std::setprecision(12) << log_p;
          outFile << "\n";

          MATRIX *inv_m_L = NULL;
          inv_m_L = MatrixInverse( (MATRIX*)m_L_tmp, inv_m_L );
#if 0
          // Check against original matrix
          CheckInverseTranslation( inv_m_L,
                                   m_L_tmp->rptr[1][4],
                                   m_L_tmp->rptr[2][4],
                                   m_L_tmp->rptr[3][4] );
#else
          // Check against base matrix + translation
          CheckInverseTranslation( inv_m_L,
                                   m_L->rptr[1][4] + x_trans,
                                   m_L->rptr[2][4] + y_trans,
                                   m_L->rptr[3][4] + z_trans );
#endif
          MatrixFree( &inv_m_L );
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
    max_log_p = CUDA_ComputeLogSampleProbability( m_L_tmp, clamp );
#else
    max_log_p = local_GCAcomputeLogSampleProbability
      (gca, gcas, mri, m_L,nsamples, exvivo, clamp) ;
#endif

#if 1
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

#ifdef OUTPUT_STAGES
  std::cerr << __FUNCTION__
            << ": Early exit for debugging"
            << std::endl;
  exit( 0 );
#endif

  return(max_log_p) ;
}
