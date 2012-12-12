/**
 * @file  gcamfots.cu
 * @brief Implementation of gcamFindOptimalTimeStep for the GPU
 *
 * Reference:
  * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
  * Structures in the Human Brain", Fischl et al.
  * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.12 $
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

#ifdef GCAMORPH_ON_GPU

#include "macros.h"

#include "chronometer.hpp"

#include "mriframegpu.hpp"
#include "gcamorphgpu.hpp"
#include "gcamorphenergy.hpp"

#include "gcamfots_cuda.hpp"
#include "gcamfots_cuda.h"

#define SHOW_TIMERS 1

// ========================================================================

// Stolen from gcamorph.c
const unsigned int MAX_SAMPLES = 100;


static int finitep( const float f )
{
  if( !finite(f) )
  {
    return 0;
  }
  else
  {
    return 1;
  }
}


// ========================================================================

template<typename T>
float FindOptimalTimestep( GPU::Classes::GCAmorphGPU& gcam,
                           GCA_MORPH_PARMS *parms,
                           const GPU::Classes::MRIframeGPU<T>& mri )
{

  GPU::Algorithms::GCAmorphEnergy gcamEnergy;
#if SHOW_TIMERS
  SciGPU::Utilities::Chronometer tFOTS;

  tFOTS.Start();
#endif

  float min_dt = 0;
  float min_rms, rms;
  float orig_dt, start_dt, max_dt;
  float dt_in[MAX_SAMPLES], rms_out[MAX_SAMPLES];


  gcam.ClearMomentum();

  min_rms = gcamEnergy.ComputeRMS( gcam, mri, parms );

  // Find right order of magnitude for time step
  orig_dt = parms->dt;
  if( DZERO( orig_dt ) )
  {
    orig_dt = 1e-6;
    std::cerr << "orig_dt = 0 in FindOptimalTimestep" << std::endl;
  }

  // Define a pretty broad search range initially
  start_dt = (sqrt(parms->navgs)+1)*orig_dt / (10*16.0*16.0);
  max_dt =   (sqrt(parms->navgs)+1)*orig_dt * (10*16.0*16.0);


  int i = 0;
  do
  {

    if( i++ > 0 )
    {
      start_dt *= 0.01;
    }

    if( DEQUAL( start_dt, 0 ) )
    {
      break;
    }

    for( parms->dt = start_dt; parms->dt <= max_dt; parms->dt *= 4 )
    {

      gcam.ApplyGradient( parms );
      rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
      gcam.UndoGradient();

      if( (gcam.neg > 0) && 0 )
      {
        std::cerr << __FUNCTION__
                  << ": How did you get here?" << std::endl;
        std::cerr << "Exiting on line " << __LINE__ << std::endl;
        exit( EXIT_FAILURE );
      }
      else
      {
        //prev_neg = suppressed = 0;
      }


      if( (gcam.neg>0) && DEQUAL(parms->dt, start_dt) )
      {
        start_dt *= 0.1;
        parms->dt /= 4;
        continue ;
      }

      if( (gcam.neg>0) && parms->noneg == True )
      {
        break;
      }

      if( rms < min_rms )
      {
        min_rms = rms;
        min_dt = parms->dt;
      }

      if( DEQUAL( min_dt, max_dt ) )
      {
        max_dt *= 10;
      }
    }

    if( i > 10 )
    {
      break;
    }

    max_dt = start_dt*10 ; /* will only iterate if min was at start_dt */


  }
  while( DEQUAL( min_dt, start_dt ) );


  // I believe that we now have a start point to do something more advanced

  // The next few blocks of code look very similar.....
  dt_in[0] = min_dt ;
  rms_out[0] = min_rms ;
  parms->dt = dt_in[1] = dt_in[0] - dt_in[0]*.2 ;
  gcam.ApplyGradient( parms ) ;
  rms = rms_out[1] = gcamEnergy.ComputeRMS( gcam, mri, parms );
  gcam.UndoGradient() ;

  if( (rms < min_rms) &&
      ( (gcam.neg==0) || (parms->noneg==False) ) )
  {
    min_rms = rms;
    min_dt = parms->dt;
  }


  parms->dt = dt_in[2] = dt_in[0] - dt_in[0]*.4 ;
  gcam.ApplyGradient( parms );
  rms = rms_out[2] = gcamEnergy.ComputeRMS( gcam, mri, parms );
  gcam.UndoGradient();

  if( (rms < min_rms) &&
      ( (gcam.neg==0) || (parms->noneg==False) ) )
  {
    min_rms = rms;
    min_dt = parms->dt;
  }


  parms->dt = dt_in[3] = dt_in[0] + dt_in[0]*.2 ;
  gcam.ApplyGradient( parms );
  rms = rms_out[3] = gcamEnergy.ComputeRMS( gcam, mri, parms );
  gcam.UndoGradient() ;

  if( (rms < min_rms) &&
      ( (gcam.neg==0) || (parms->noneg==False) ) )
  {
    min_rms = rms;
    min_dt = parms->dt;
  }


  parms->dt = dt_in[4] = dt_in[0] + dt_in[0]*.4 ;
  gcam.ApplyGradient( parms );
  rms = rms_out[4] = gcamEnergy.ComputeRMS( gcam, mri, parms );
  gcam.UndoGradient() ;

  if( (rms < min_rms) &&
      ( (gcam.neg==0) || (parms->noneg==False) ) )
  {
    min_rms = rms;
    min_dt = parms->dt;
  }


  /* now compute location of minimum of best quadratic fit */
  MATRIX   *mX, *m_xTx, *m_xTx_inv, *m_xTy, *mP, *m_xT;
  VECTOR   *vY;
  int N = 5 ;  /* min_dt +- .1*min_dt and .2*min_dt */
  mX = MatrixAlloc(N, 3, MATRIX_REAL) ;
  vY = VectorAlloc(N, MATRIX_REAL) ;

  for( i=1; i<=N; i++ )
  {
    *MATRIX_RELT(mX, i, 1) = dt_in[i-1] * dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 2) = 2*dt_in[i-1] ;
    *MATRIX_RELT(mX, i, 3) = 1.0f ;

    VECTOR_ELT(vY, i) = rms_out[i-1] ;
  }

  m_xT = MatrixTranspose(mX, NULL) ;
  m_xTx = MatrixMultiply(m_xT, mX, NULL) ;
  m_xTx_inv = MatrixInverse(m_xTx, NULL) ;

  if (m_xTx_inv)
  {
    float a, b;

    m_xTy = MatrixMultiply(m_xT, vY, NULL) ;
    mP = MatrixMultiply(m_xTx_inv, m_xTy, NULL) ;
    a = RVECTOR_ELT(mP, 1) ;
    b = RVECTOR_ELT(mP, 2) ;
    //c = RVECTOR_ELT(mP, 3);

    MatrixFree(&mP) ;
    MatrixFree(&m_xTx_inv) ;
    MatrixFree(&m_xTy) ;

    if (finitep(a) && !FZERO(a))
    {
      parms->dt = -b/a ;
      gcam.ApplyGradient( parms );
      rms = gcamEnergy.ComputeRMS( gcam, mri, parms );
      gcam.UndoGradient();

      if( (rms<min_rms) &&
          ( (gcam.neg==0) || (parms->noneg == False) ) )
      {
        min_rms = rms ;
        min_dt = parms->dt ;
      }
    }
  }

  MatrixFree(&m_xT) ;
  MatrixFree(&m_xTx) ;
  MatrixFree(&mX) ;
  VectorFree(&vY) ;

  int invalid;
  gcam.ComputeMetricProperties( invalid );

  parms->dt = orig_dt ;

#if SHOW_TIMERS
  tFOTS.Stop();
  std::cout << __FUNCTION__
            << ": Complete in "
            << tFOTS << std::endl;
#endif

  return( min_dt );
}





// ============================================

template<typename T>
float FindOptimalTimestepDispatch( GCA_MORPH *gcam,
                                   GCA_MORPH_PARMS *parms,
                                   MRI *mri )
{

#if SHOW_TIMERS
  SciGPU::Utilities::Chronometer tFOTSd;

  tFOTSd.Start();
#endif

  GPU::Classes::GCAmorphGPU myGCAM;
  myGCAM.CheckIntegrity(); // Shouldn't be necessary....
  myGCAM.SendAll( gcam );

  GPU::Classes::MRIframeGPU<T> myMRI;
  myMRI.Allocate( mri );
  myMRI.Send( mri, 0 );


  float dt = FindOptimalTimestep( myGCAM, parms, myMRI );

  myGCAM.RecvAll( gcam );

#if SHOW_TIMERS
  tFOTSd.Stop();
  std::cout << __FUNCTION__
            << ": Complete in "
            << tFOTSd << std::endl;
#endif

  return( dt );
}

// ============================================


float gcamFindOptimalTimestepGPU( GCA_MORPH *gcam,
                                  GCA_MORPH_PARMS *parms,
                                  MRI *mri )
{

  float dt;

  switch( mri->type )
  {

  case MRI_UCHAR:
    dt = FindOptimalTimestepDispatch<unsigned char>( gcam, parms, mri );
    break;

  case MRI_FLOAT:
    dt = FindOptimalTimestepDispatch<float>( gcam, parms, mri );
    break;

  default:
    std::cerr << __FUNCTION__
              << ": Unrecognised MRI type" << std::endl;
    exit( EXIT_FAILURE );
  }

  return( dt );

}


#endif
